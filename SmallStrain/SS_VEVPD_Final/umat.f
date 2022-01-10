C     UMAT - Isotropic Linear Viscoelastic - Perznya type J2 Viscoplastic 
C     (Linear Isotropic, Kinematic Hardening) coupled with Lamitre type Ductile Damage   

C     Mohib Mustafa - IMDEA 4 MAR 2021      

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,
     2           DPRED,CMNAME,
     3           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     3           DROT,PNEWDT,
     4           CELENT,DFGRD0,DFGRD1,NOEL,NPT,
     5           LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC' 

      DIMENSiON STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),DFGRD0(3,3),
     4     DFGRD1(3,3)


      call VE_VP_D(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &         NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,TIME,PNEWDT)

      return
      end     

C     !-------------------------------------------------------------------------
C     !   UMAT SUBROUTINE FOR ISOTROPIC LINEAR VISCOELASTIC - VISCOPLASTIC 
C     !    MATERIAL WITH LINEAR ISOTROPIC AND KINEMATIC HARDENING COUPLED 
C     !                   WITH LAMITRE TYPE DUCTILE DAMAGE
C     !-------------------------------------------------------------------------
      SUBROUTINE VE_VP_D(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,
     &        PROPS,NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,TIME,PNEWDT)
     
     
        IMPLICIT NONE

        INTERFACE
          FUNCTION dotProduct(voit_array)
            INTEGER, PARAMETER :: double=kind(1.d0)
            INTEGER, PARAMETER :: prec=double

            REAL(prec), INTENT(IN) :: voit_array(6)
            REAL(prec) :: dotProduct
          END FUNCTION dotProduct
        END INTERFACE
      
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
       
        INTEGER, INTENT(IN)      :: ntens,nprops,nstatv
        INTEGER, INTENT(IN)      :: kinc,kstep,noel
        REAL(prec), INTENT(IN)   :: stran(ntens), dstran(ntens)
        REAL(prec), INTENT(IN)   :: props(nprops), dtime, time(2)
        REAL(prec), INTENT(INOUT) :: statev(nstatv)
        REAL(prec), INTENT(OUT)   :: stress(ntens)
        REAL(prec), INTENT(OUT)   :: ddsdde(ntens,ntens), PNEWDT
      
        !List of internal variables
        INTEGER    :: rap_iter, K1, K2
        REAL(prec) :: xioi(6, 6), xii(6, 6), xpp(6, 6) 
        REAL(prec) :: retoret(6, 6), sigoret(6, 6)
        REAL(prec) :: EPS(6), alpha(6, 3), EPS_VP(6), xi, D
        REAL(prec) :: treps, dev(6)
        REAL(prec) :: BETA_VE_TR(6, 3), SIGMA_DEV_TR(6), beta_sc_tr
        REAL(prec) :: BETA_TR(6), zeta_tr(6), mod_zeta_tr, ret(6)
        REAL(prec) :: SIGMA_TR(6), SIGMA(6)
        REAL(prec) :: CHI_TR
        REAL(prec) :: GAM_i, alpha_i(6, 3), EPS_VP_i(6), xi_i
        REAL(prec) :: E_ve(6), E_e(6, 3)
        REAL(prec) :: Y_i, Y_ve, Y_e, Y_vphh
        REAL(prec) :: G_i, res_i
        REAL(prec) :: facve(3), fac1, fac2, fac3, fac4, fac5, fac6 
        REAL(prec) :: fac7, BB1, BB2, BB3, BB4, BB5, BB6   
        REAL(prec) :: ka, mu0, mu(3), eta(3) 
        REAL(prec) :: h, hh, y0, eta0, m, s, ss, dflag, dump_flag
       
        REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
        REAL(prec), PARAMETER :: THREE=3.D0, SIX=6.D0, ENUMAX=.4999D0
        REAL(prec), PARAMETER :: NEWTON=10, TOLER=1.0D-6

        !Define 4th order identity tensor
        data xioi(1,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
        data xioi(2,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
        data xioi(3,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
        data xioi(4,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
        data xioi(5,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
        data xioi(6,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/

        !Define 4th order symmetric identity tensor
        data xii(1, :) /ONE, ZERO, ZERO, ZERO, ZERO, ZERO/
        data xii(2, :) /ZERO, ONE, ZERO, ZERO, ZERO, ZERO/
        data xii(3, :) /ZERO, ZERO, ONE, ZERO, ZERO, ZERO/
        data xii(4, :) /ZERO, ZERO, ZERO, 0.5D0, ZERO, ZERO/
        data xii(5, :) /ZERO, ZERO, ZERO, ZERO, 0.5D0, ZERO/
        data xii(6, :) /ZERO, ZERO, ZERO, ZERO, ZERO, 0.5D0/
       
        !Compute deviatoric projection tensor
        xpp(:, :) = xii(:, :) - (ONE / THREE) * xioi(:, :)
       
        !Get material properties
        ka=PROPS(1)
        mu0=PROPS(2)
        mu(1)=PROPS(3)
        mu(2)=PROPS(4)
        mu(3)=PROPS(5)
        eta(1)=PROPS(6)
        eta(2)=PROPS(7)
        eta(3)=PROPS(8)
        h=PROPS(9)
        hh=PROPS(10)
        y0=PROPS(11)
        eta0=PROPS(12)
        m=PROPS(13)
        s=PROPS(14)
        ss=PROPS(15)
        dflag=PROPS(16)
                  
        !DEFINE STATE VARIABLES      
        EPS(:) = DSTRAN(:) + STRAN(:)
        alpha(:, 1) = STATEV(1 : 6)
        alpha(:, 2) = STATEV(7 : 12)
        alpha(:, 3) = STATEV(13 : 18)
        EPS_VP(:) = STATEV(19 : 24)
        xi = STATEV(25)
        D = STATEV(26)
        dump_flag = STATEV(27)
            
        !Trace of strain at time step n+1    
        treps = EPS(1) + EPS(2) + EPS(3)
      
        !Deviator of strain at time step n+1
        dev(1 : 3) = EPS(1 : 3) - treps / THREE
        dev(4 : 6) = EPS(4 : 6) / TWO
              
        !Define algorithmic constants
        facve(:) = (TWO * mu(:)) 
     1           / (ONE + (TWO * mu(:) * DTIME) / eta(:))
        fac1 = TWO * mu0 + facve(1) + facve(2) + facve(3)
        fac2 = fac1 + hh
        fac3 = fac2 + (TWO/THREE) * h
        fac4 = fac3 / y0
        fac5 = eta0 / (DTIME * y0)
        fac6 = fac1/y0
        fac7 = ((fac5 * s) / ss) * (TWO/THREE)**0.5D0

        !-----------------------------------------------------------------------
        ! Compute ViscoElastic Predictor (Trial) State for matrix material:     
        !-----------------------------------------------------------------------
      
        !Calculate BETA_VE at Trial for material matrix
        DO K1 = 1, 3
          BETA_VE_TR(:, K1) = facve(K1) 
     1                      * (dev(:) - EPS_VP(:) - alpha(:, K1))
        END DO

        !Calculate deviatoric Sigma at trial for matrix material
        SIGMA_DEV_TR(:) = TWO * mu0 * (dev(:) - EPS_VP(:)) 
     1                  + BETA_VE_TR(:, 1) + BETA_VE_TR(:, 2) 
     2                  + BETA_VE_TR(:, 3)

        !Calculate Beta scaler at trial for matrix material
        beta_sc_tr = h * xi

        !Calculate Beta at trial for matrix material
        BETA_TR(:) = hh * EPS_VP(:)
      
        !Calculate zeta(Relative stress) at trial for matrix material
        zeta_tr(:) = SIGMA_DEV_TR(:) - BETA_TR(:)

        !Calculate euler norm of zeta at trial
        CALL NORM(zeta_tr, mod_zeta_tr)

        !Calculate radial return tensor N
        ret(:) = zeta_tr(:) /  mod_zeta_tr
        
        !Calculate N o N
        DO K1 = 1, 6
          DO K2 = 1, 6
            retoret(K1, K2) = ret(K1) * ret(K2)
          END DO
        END DO

        !Calculate Stress at trial for matrix material
        SIGMA_TR(1 : 3) = ka * treps + SIGMA_DEV_TR(1 : 3)
        SIGMA_TR(4 : 6) = SIGMA_DEV_TR(4 : 6)
      
        !Calculate Yield surface at trial
        CHI_TR = mod_zeta_tr 
     1         - (y0 + beta_sc_tr) * (TWO/THREE)**0.5D0

        !-----------------------------------------------------------------------      
        !         IF Viscoelastic Domain at trial state is not breached 
        !                         Return trial state 
        !-----------------------------------------------------------------------
        IF ((CHI_TR .LE. ZERO) .OR. (dump_flag .GE. 5.D0)) THEN
      
          !Update Internal Variables - alphas need to be updated at this point
          !(Retundent for eps_vp, xi and D)
          STATEV(1 : 6) = STATEV(1 : 6) 
     1                  + (DTIME/eta(1)) * BETA_VE_TR(:, 1)
        
          STATEV(7 : 12) = STATEV(7 : 12) 
     2                   + (DTIME/eta(2)) * BETA_VE_TR(:, 2)

          STATEV(13 : 18) = STATEV(13 : 18) 
     3                    + (DTIME/eta(3)) * BETA_VE_TR(:, 3)

          STATEV(19 : 24) = EPS_VP(:)
          STATEV(25) = xi
          STATEV(26) = D
        
          !Update Stress and Consistent Tangent modulus to trial values 
          !after converting them to effective values for the RVE        
          STRESS(:) = SIGMA_TR(:) * (ONE - D)
          DDSDDE(:, :) = (ka * xioi(:, :) + fac1 * xpp(:, :)) 
     1                 * (ONE - D)
            
        !-----------------------------------------------------------------------      
        !             IF Viscoelastic Domain at trial state is breached 
        !             Find Plastic multipler Gamma and Correct trial state 
        !-----------------------------------------------------------------------     
        ELSE
          
          !----------------------------------------------------------------
          !           START OF NEWTON RAPHSON SOLVER FOR GAM_i
          !----------------------------------------------------------------

          !Initialize counter and initial values
          rap_iter = 0
          GAM_i = ZERO
        
          !Set values of internal variables for first newton raphson iteration
          DO K1 = 1,3
            alpha_i(:, K1) = alpha(:, K1) 
     1                     + (DTIME/eta(K1)) * BETA_VE_TR(:, K1)
          END DO

          EPS_VP_i(:) = EPS_VP(:)

          xi_i = xi

          !Compute algorithmic constants for GAM_i
          E_ve(:) = dev(:) - EPS_VP_i(:)
        
          DO K1 = 1, 3
            E_e(:, K1) = E_ve(:) - alpha_i(:, K1)
          END DO
        
          !Compute various portions of effective free energy for Gam_i.
          !(Helps in coding. Can be done in one line too)         
          Y_ve = mu0 * dotProduct(E_ve(:))
          Y_e = mu(1) * dotProduct(E_e(:,1)) 
     1        + mu(2) * dotProduct(E_e(:,2))
     2        + mu(3) * dotProduct(E_e(:,3))
          Y_vphh = 0.5d0 * hh * dotProduct(EPS_VP_i(:))

          !Compute effective free energy
          Y_i = 0.5d0 * treps ** TWO + Y_ve + Y_e 
     1        + 0.5d0 * h * xi_i ** TWO + Y_vphh          
        
          !Compute Algorithmic Constants for Gam_i
          BB4 = (CHI_TR / y0) - fac4 * GAM_i
          BB5 = dflag * GAM_i 
     1        * ((TWO/THREE) ** 0.5d0) * ((Y_i/ss) ** s)
          BB6 = ONE - D - BB5

          ! IF (BB6 .LE. 0.00001D0) THEN
          !   BB6 = 0.00001D0
          ! END IF


          !Compute resuidual and tangent for newton raphson
          G_i = m * fac4 * (BB4) ** (m - ONE) 
     1        - fac5 * BB5 + fac5 * BB6         
          res_i = BB4 ** m - fac5 * GAM_i * BB6
        
          DO WHILE((res_i .GT. TOLER) .AND. (rap_iter .LT. 25))
            rap_iter = rap_iter + 1
            GAM_i = GAM_i + (ONE/G_i)*res_i

            !Set values of internal variables for updated GAM_i
            DO K1 = 1,3
              alpha_i(:, K1) = alpha(:, K1) 
     1                     + (DTIME/eta(K1)) * BETA_VE_TR(:, K1)
     2                     - (DTIME/eta(K1)) * GAM_i * ret(:)
            END DO

            EPS_VP_i(:) = EPS_VP(:) + GAM_i * ret(:)
            xi_i = xi + GAM_i * (TWO/THREE)**0.5d0

            !Compute algorithmic constants for GAM_i
            E_ve(:) = dev(:) - EPS_VP_i(:)
        
            DO K1 = 1, 3
              E_e(:, K1) = E_ve(:) - alpha_i(:, K1)
            END DO

            !Compute various portions of effective free energy.
            !(Helps in coding. Can be done in one line too)         
            Y_ve = mu0 * dotProduct(E_ve(:))
            Y_e = mu(1) * dotProduct(E_e(:,1)) 
     1          + mu(2) * dotProduct(E_e(:,2))
     2          + mu(3) * dotProduct(E_e(:,3))
            Y_vphh = 0.5d0 * hh * dotProduct(EPS_VP_i(:))

            !Compute effective free energy
            Y_i = 0.5d0 * treps ** TWO + Y_ve + Y_e 
     1          + 0.5d0 * h * xi_i ** TWO + Y_vphh

            !Compute Algorithmic Constants
            BB4 = (CHI_TR / y0) - fac4 * GAM_i
            BB5 = dflag * GAM_i 
     1          * ((TWO/THREE) ** 0.5d0) * ((Y_i/ss) ** s)
            BB6 = ONE - D - BB5

            ! IF (BB6 .LE. 0.00001D0) THEN
            !   BB6 = 0.00001D0
            ! END IF

            !Compute resuidual and tangent for newton raphson
            G_i = m * fac4 * (BB4) ** (m - ONE) 
     1           - fac5 * BB5 + fac5 * BB6        
            res_i = BB4 ** m - fac5 * GAM_i * BB6

          END DO
        
          !----------------------------------------------------------------
          !           END OF NEWTON RAPHSON SOLVER FOR GAM_i
          !----------------------------------------------------------------

          IF (rap_iter .GE. 25) THEN
            PNEWDT = 0.5d0
            RETURN
          
          ELSE IF (D + BB5 .GE. 0.99999D0) THEN
            
            STATEV(27) = 6.D0
            ! WRITE(*,*) 'Dump : ', STATEV(27)

            STATEV(1 : 6) = alpha_i(:, 1)        
            STATEV(7 : 12) = alpha_i(:, 2)
            STATEV(13 : 18) = alpha_i(:, 3)
            STATEV(19 : 24) = EPS_VP_i(:)
            STATEV(25) = xi_i
            STATEV(26) = 0.99999D0

          !Update Stress and Consistent Tangent modulus to trial values 
          !after converting them to effective values for the RVE        
            STRESS(:) = SIGMA_TR(:) * (ONE - STATEV(26))
            DDSDDE(:, :) = (ka * xioi(:, :) + fac1 * xpp(:, :)) 
     1                   * (ONE - STATEV(26))


            RETURN

          ELSE
            !Pass converged internal variables from Newton Raphson solver
            STATEV(1 : 6) = alpha_i(:, 1)        
            STATEV(7 : 12) = alpha_i(:, 2)
            STATEV(13 : 18) = alpha_i(:, 3)
            STATEV(19 : 24) = EPS_VP_i(:)
            STATEV(25) = xi_i
            STATEV(26) = D + BB5

            ! IF (STATEV(26) .GE. 0.999999D0) THEN
            !   STATEV(27) = 6.D0
            ! WRITE(*,*) 'Dump : ', STATEV(27)
            ! END IF

            ! WRITE(*,*) 'alpha_1 : ', STATEV(1 : 6)

            ! WRITE(*,*) ' '
            ! WRITE(*,*) 'alpha_2 : ', STATEV(7 : 12)

            ! WRITE(*,*) ' '
            ! WRITE(*,*) 'alpha_3 : ', STATEV(13 : 18)
            
            ! WRITE(*,*) ' '
            ! WRITE(*,*) 'eps_vp : ', STATEV(19 : 24)
            
            ! WRITE(*,*) ' '
            ! WRITE(*,*) 'xi_i : ', STATEV(25)
            
            ! WRITE(*,*) ' '
            ! WRITE(*,*) 'D : ', STATEV(26)

            
            ! WRITE(*, *) 'R : ', res_i

            ! WRITE(*,*) ' '

            ! WRITE(*, *) 'G : ', G_i

            

            ! WRITE(*, *) '-------------------------------------------'

            !Correct Sigma for matrix material
            SIGMA(:) = SIGMA_TR(:) - fac1 * GAM_i * ret(:)

            !Update Stress for ABAQUS
            STRESS(:) = SIGMA(:) * (ONE - STATEV(26))

            !Compute algorithmic Constants
            BB1 = ONE - (fac1 * GAM_i) / mod_zeta_tr
            BB2 = ((fac6 * m) / G_i) * (BB4) ** (m - ONE) 
     1          - (fac1 * GAM_i) / mod_zeta_tr
            BB3 = dflag * (fac1 / (G_i)) * ((Y_i/ss) ** (s - ONE)) 
     1          * (GAM_i ** TWO)
            
            !Compute sigoret 4th order tensor
            DO K1 = 1, 6
              DO K2 = 1, 6
                sigoret(K1, K2) = SIGMA(K1) * ret(K2)
              END DO
            END DO

            DDSDDE(:, :) = (ka * xioi(:, :)
     1                   + fac1 * BB1 * xpp(:, :)
     2                   - fac1 * BB2 * retoret(:, :)
     3                   - fac1 * BB3 * sigoret(:, :)) 
     4                   * (ONE - STATEV(26))

          END IF
              
        END IF
        RETURN
      END SUBROUTINE VE_VP_D
      
      !-------------------------------------------------------------------------
      !           Helper function to compute Euler norm
      !-------------------------------------------------------------------------
      SUBROUTINE NORM(voit_array, output)
        IMPLICIT NONE
        
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: voit_array(6)
        REAL(prec), INTENT(OUT) :: output
        INTEGER :: K1
      
        output=0.d0
        DO K1 = 1, 3
            output = output + voit_array(K1)**2.0d0
            output = output + 2.d0 * voit_array(K1 + 3)**2.d0
        END DO
        
        output = output**0.5d0
        RETURN
      END SUBROUTINE NORM
      
      !-------------------------------------------------------------------------
      !       Helper function to compute dot product in voit notation
      !-------------------------------------------------------------------------
      FUNCTION dotProduct(voit_array)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: voit_array(6)
        REAL(prec) :: dotProduct
        INTEGER :: K1

        dotProduct = 0.d0
        DO K1 = 1, 3
          dotProduct = dotProduct + voit_array(K1)**2.0d0
          dotProduct = dotProduct + 2.d0 * voit_array(K1 + 3)**2.d0
        END DO
      END FUNCTION dotProduct
