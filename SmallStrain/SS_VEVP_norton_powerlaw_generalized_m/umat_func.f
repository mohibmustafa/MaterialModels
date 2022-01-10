C Mohib Mustafa - IMDEA 4 MAR 2021      
C UMAT - Isotropic Viscoplastic - Perznya with M Exponent
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


      call VE_VP(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,NPROPS,
     &         DTIME,DSTRAN,KINC,KSTEP,NOEL,TIME)

      return
      end      

      SUBROUTINE VE_VP(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &        NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL, TIME)
     
     
       IMPLICIT NONE
      
       INTEGER, PARAMETER :: double=kind(1.d0)
       INTEGER, PARAMETER :: prec=double
       
       INTEGER, INTENT(IN)      :: ntens,nprops,nstatv,kinc,kstep,noel
       REAL(prec), INTENT(IN)   :: stran(ntens), dstran(ntens)
       REAL(prec), INTENT(IN)   :: props(nprops), dtime, time(2)
       REAL(prec), INTENT(INOUT) :: statev(nstatv)
       REAL(prec), INTENT(OUT)   :: stress(ntens),ddsdde(ntens,ntens)
      
      
C      List of internal variables:
       INTEGER    :: rap_iter, K1, K2
       REAL(prec) :: xioi(6, 6), xii(6, 6), xpp(6, 6), retoret(6, 6)
       REAL(prec) :: ret(6), dev(6), BETA_TR(6), BETA_VE_TR(6, 3)
       REAL(prec) :: SIGMA_BAR_TR(6), SIGMA_BAR(6), EPS_VP(6)
       REAL(prec) :: EPS(6), SIGMA_TR(6), TANG_TR(6, 6), alpha(6, 3)
       REAL(prec) :: treps, xi, CHI_TR, GAM_i, BETA_VE(6, 3)
       REAL(prec) :: fac1, zeta_trial, zeta_tr(6), mod_zeta_tr
       REAL(prec) :: fac2, fac3, res_i, G_i, beta_sc_tr
       REAL(prec) :: ka, mu0, mu(3), h, hh, y0, eta0, eta(3), m
       REAL(prec) :: BB1, BB2, facve(3), fac4, fac5, fac6, BB3
       
       REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
       REAL(prec), PARAMETER :: THREE=3.D0, SIX=6.D0, ENUMAX=.4999D0
       REAL(prec), PARAMETER :: NEWTON=10, TOLER=1.0D-6

C      Define 4th order identity tensor
       data xioi(1,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
       data xioi(2,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
       data xioi(3,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
       data xioi(4,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
       data xioi(5,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
       data xioi(6,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/

C      Define 4th order symmetric identity tensor
       data xii(1, :) /ONE, ZERO, ZERO, ZERO, ZERO, ZERO/
       data xii(2, :) /ZERO, ONE, ZERO, ZERO, ZERO, ZERO/
       data xii(3, :) /ZERO, ZERO, ONE, ZERO, ZERO, ZERO/
       data xii(4, :) /ZERO, ZERO, ZERO, 0.5D0, ZERO, ZERO/
       data xii(5, :) /ZERO, ZERO, ZERO, ZERO, 0.5D0, ZERO/
       data xii(6, :) /ZERO, ZERO, ZERO, ZERO, ZERO, 0.5D0/
       
       
c       open(unit = 10, file = "debug.txt")
       
C       write(*, *) "-----------------------------------------"
C       write(*, *) "Time step : ", TIME(2)
C       write(*, *) "-----------------------------------------"
       
C     Compute deviatoric projection tensor
       xpp(:, :) = xii(:, :) - (ONE / THREE) * xioi(:, :)
       
c       write(10, *) "xioi"
c       DO K1 = 1, 6
c         write(10, *) xioi(K1, :)
c       END DO
       
c       write(10, *) "xii"
c       DO K1 = 1, 6
c         write(10, *) xii(K1, :)
c       END DO
       
c       write(10, *) "xpp"
c       DO K1 = 1, 6
c         write(10, *) xpp(K1, :)
c       END DO
       
       
       
C     Get material properties
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
      
C      write(*, *) "Material Properties"
C      write(*, *) " "
C      write(*, *) "ka : ", ka
C      write(*, *) "mu0 : ", mu0
C      write(*, *) "h : ", h
C      write(*, *) "hh : ", hh
C      write(*, *) "y0 : ", y0
C      write(*, *) "eta0 : ", eta0
C      write(*, *) "m : ", m
      
      
C     DEFINE STATE VARIABLES      
      EPS(:) = DSTRAN(:) + STRAN(:)
      alpha(:, 1) = STATEV(1 : 6)
      alpha(:, 2) = STATEV(7 : 12)
      alpha(:, 3) = STATEV(13 : 18)
      EPS_VP(:) = STATEV(19 : 24)
      xi = STATEV(25)
      
C      write(*, *) " "
C      write(*, *) "State Variables"
C      write(*, *) " "
C      write(*, *) "EPS"
C      write(*,*) EPS(:)
C      write(*, *) "EPS_VP"
C      write(*,*) EPS_VP(:)
C      write(*, *) "xi : ", xi
C      WRITE(*,*) "DTIME", DTIME
      
C     Trace of strain at time step n+1    
      treps = EPS(1) + EPS(2) + EPS(3)
      
C     Deviator of strain at time step n+1
      dev(1 : 3) = EPS(1 : 3) - treps / THREE
      dev(4 : 6) = EPS(4 : 6) / TWO
        
C      write(*,*) "treps : ", treps
C      write(*,*) "dev"
C      write(*, *) dev(:)
      
C     Define algorithmic constants
      facve(:) = (TWO * mu(:)) / (ONE + (TWO * mu(:) * DTIME) / eta(:))
      fac1 = TWO * mu0 + facve(1) + facve(2) + facve(3)
      fac2 = fac1 + hh
      fac3 = fac2 + (TWO/THREE) * h
      fac4 = fac3/y0
      fac5 = eta0 / (DTIME * y0)
      fac6 = fac1/y0

C     Compute Elastic Predictor (Trial) State:     

C     Calculate BETA_VE at trial
      DO K1 = 1, 3
        BETA_VE_TR(:, K1) = facve(K1) * (dev(:) - EPS_VP(:) - alpha(:, K1))
      END DO

C     Calculate Sigma_bar at trial
      SIGMA_BAR_TR(:) = TWO * mu0 * (dev(:) - EPS_VP(:)) 
     1                + BETA_VE_TR(:, 1) + BETA_VE_TR(:, 2) 
     2                + BETA_VE_TR(:, 3)

C      write(*,*) " "
C      write(*,*) "Sigma_bar_trial"
C      write(*,*) SIGMA_BAR_TR(:)

C     Calculate Beta scaler at trial
      beta_sc_tr = h * xi

C     Calculate Beta_VP at trial
      BETA_TR(:) = hh * EPS_VP(:)
C      write(*,*) " "
C      write(*,*) "beta_vp_trial"
C      write(*,*) BETA_TR(:)
      
C     Calculate xi at trial
      zeta_tr(:) = SIGMA_BAR_TR(:) - BETA_TR(:)

C     Calculate euler norm of xi at trial
      CALL NORM(zeta_tr, mod_zeta_tr)
C      write(*,*) "||zeta_tr|| : ", mod_zeta_tr

C     Calculate N_NP
      ret(:) = zeta_tr(:) /  mod_zeta_tr
        
C      write(*,*) ""
C      write(*,*) "N"
C      write(*,*) ret(:)


C     Calculate N_NP diadic N_NP
      DO K1 = 1, 6
        DO K2 = 1, 6
          retoret(K1, K2) = ret(K1) * ret(K2)
        END DO
      END DO
C      write(*,*) ""
C      write(*,*) "N x N"

C      DO K1 = 1, 6
C      write(*,*) retoret(K1, :)
C      End Do


C     Calculate Stress
      SIGMA_TR(1 : 3) = ka * treps + SIGMA_BAR_TR(1 : 3)
      SIGMA_TR(4 : 6) =  SIGMA_BAR_TR(4 : 6)
      
C      write(*,*) " "
C      write(*,*) "Sigma_tr"
C      write(*,*) SIGMA_TR(:)
      
C     Calculate consistent tangent modulus at trial
      TANG_TR(:, :) = ka * xioi(:, :) + fac1 * xpp(:, :)
c      write(*,*) " "
c      write(*,*) "tang_tr"
c      DO K1 = 1, 6
c       write(*,*) TANG_TR(K1, :)
c      END DO
      
      
C     Calculate Yield surface at trial
      CHI_TR = mod_zeta_tr 
     1        - (y0 + beta_sc_tr) * (TWO/THREE)**0.5D0
      
C      write(*,*) " "
C      write(*,*) "Yield Trial : ", CHI_TR
      
C     Check for breach of Yield Surface
      IF (CHI_TR .LE. ZERO) THEN
      

C       Update Internal Variables (alphas need to be updated at this point. retundent for eps_vp and xi)
        STATEV(1 : 6) = STATEV(1 : 6) 
     1                + (DTIME/eta(1)) * BETA_VE_TR(:, 1)
        
        STATEV(7 : 12) = STATEV(7 : 12) 
     2                 + (DTIME/eta(2)) * BETA_VE_TR(:, 2)

        STATEV(13 : 18) = STATEV(13 : 18) 
     3                  + (DTIME/eta(3)) * BETA_VE_TR(:, 3)

        STATEV(19 : 24) = EPS_VP(:)
        STATEV(25) = xi
        
C       If witin the elastic limit update Stress and Consistent Tangent modulus to trial values
        
        STRESS(:) = SIGMA_TR(:)
        DDSDDE(:, :) = TANG_TR(:, :)
 
C     Corrector steps for Radial return in case Breach of Yield surface     
      ELSE
        
C       ----------------------------------------------------------------
C                  NEWTON RAPHSON SOLVER FOR GAM_i
C       ----------------------------------------------------------------

C       Initialize counter and initial values
        rap_iter = 0
        GAM_i = ZERO
        
C       Set values for first newton raphson iteration
        
        BB3 = (ONE / y0) * CHI_TR - fac4 * GAM_i
        
        res_i = BB3 ** m - fac5 * GAM_i
        G_i = m * fac4 * (BB3) ** (m - ONE) + fac5

C        Write(*,*) "-------First Estimate Newton Raphson------------"
C        write(*,*) "For Gamma  = 0 "
C        write(*,*) "Residual : ", res_i
C        write(*,*) "G : ", G_i
C        Write(*,*) "------------------------------------------------"
       
        DO WHILE(res_i .GT. TOLER)
          rap_iter = rap_iter + 1
          GAM_i = GAM_i + (ONE/G_i)*res_i

          BB3 = (ONE / y0) * CHI_TR - fac4 * GAM_i
        
          G_i = m * fac4 * (BB3) ** (m - ONE) + fac5

          res_i = BB3 ** m - fac5 * GAM_i
        END DO
       
        Write(*,*) "-------End of Newton Raphson--------------------"
        write(*,*) "iterations : ", rap_iter
        write(*,*) "Residual : ", res_i
        write(*,*) "G : ", G_i
        write(*,*) "Gamma : ", GAM_i
        Write(*,*) "------------------------------------------------"
        
C       ----------------------------------------------------------------
C                  NEWTON RAPHSON SOLVER FOR GAM_i
C       ----------------------------------------------------------------

C       Correct BETA_VE at Yielding
        DO K1 = 1, 3
          BETA_VE(:, K1) = BETA_VE_TR(:, K1) - facve(K1) * GAM_i * ret(:)
        END DO


C       Correct Stress at Yielding
        STRESS(:) = SIGMA_TR(:) - fac1 * GAM_i * ret(:)
        
C        WRITE(*,*) " "
C        Write(*,*) "STRESS"
C        WRITE(*,*) STRESS(:)

C       Update Internal Variables
        STATEV(1 : 6) = STATEV(1 : 6) 
     1                + (DTIME/eta(1)) * BETA_VE(:, 1)
   
        STATEV(7 : 12) = STATEV(7 : 12) 
     2                 + (DTIME/eta(2)) * BETA_VE(:, 2)

        STATEV(13 : 18) = STATEV(13 : 18) 
     3                  + (DTIME/eta(3)) * BETA_VE(:, 3)

        STATEV(19 : 24) = STATEV(19 : 24) + GAM_i * ret(:)
        STATEV(25) = STATEV(25) + GAM_i * (TWO/THREE)**0.5d0
        
C        WRITE(*,*) " "
C        Write(*,*) "EPS_VP"
C        WRITE(*,*) STATEV(19 : 24)
       
C        WRITE(*,*) " "
C        Write(*,*) "xi : ", STATEV(25)

      
C       Compute algorithmic Constants
        BB1 = ONE - (fac1 * GAM_i) / mod_zeta_tr
        
        BB2 = ((fac6 * m) / G_i) * (BB3) ** (m - ONE) 
     1      - (fac1 * GAM_i) / mod_zeta_tr
       
        DDSDDE(:, :) = ka * xioi(:, :)
     1              + fac1 * BB1 * xpp(:, :)
     2              - fac1 * BB2 * retoret(:, :)
       
C        write(*,*) ""
C        write(*,*) "TANG"
        
C        DO K1 = 1, 6
C          write(*,*) DDSDDE(K1, :)
C        End Do
       
      
      END IF
      
      
      
C      close(10)
      
      END SUBROUTINE VE_VP
      
      
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
      
