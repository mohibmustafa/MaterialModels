C     UMAT - Finite strain J2 Plasticity
C     Formulation based on Computational Inelasticity by Simo

C     Mohib Mustafa - IMDEA 26 JULY 2021
            
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,
     3           DPRED,CMNAME,
     4           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     5           DROT,PNEWDT,
     6           CELENT,DFGRD0,DFGRD1,NOEL,NPT,
     7           LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),DFGRD0(3,3), DFGRD1(3,3)
      
      call FinJ2(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &         NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1)

      return
      end

C     !--------------------------------------------------------------
C     !   UMAT Finite Strain J2 Plasticity
C     !--------------------------------------------------------------

      SUBROUTINE FinJ2(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,
     1   PROPS,NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1)


        IMPLICIT NONE

        INTERFACE
          FUNCTION determinant(matrix)
            INTEGER, PARAMETER :: double=kind(1.d0)
            INTEGER, PARAMETER :: prec=double

            REAL(prec), INTENT(IN) :: matrix(3, 3)
            REAL(prec) :: determinant
          END FUNCTION determinant
        END INTERFACE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
       
        INTEGER, INTENT(IN)      :: ntens,nprops,nstatv
        INTEGER, INTENT(IN)      :: kinc,kstep,noel
        REAL(prec), INTENT(IN)   :: stran(ntens), dstran(ntens)
        REAL(prec), INTENT(IN)   :: props(nprops), dtime
        REAL(prec), INTENT(IN)   :: DFGRD0(3,3), DFGRD1(3,3)
        REAL(prec), INTENT(INOUT) :: statev(nstatv)
        REAL(prec), INTENT(OUT)   :: stress(ntens)
        REAL(prec), INTENT(OUT)   :: ddsdde(ntens,ntens)

        !List of internal variables
        INTEGER    :: index, ii, jj, mm, nn, K1, K2, ll
        REAL(prec) :: I(6), m2v(3, 3), Be_bar_tr(3, 3), tau_tr(3, 3)
        REAL(prec) :: F_n_inv(3, 3), ff(3, 3), BB1(6, 6), C_n(3, 3)
        REAL(prec) :: c_strike(6, 6), xioi(6, 6), xii(6, 6)
        REAL(prec) :: c_strike_tr(6, 6), F_bar_inv(3, 3), alpha_n
        REAL(prec) :: dev_Be_bar_tr(3, 3), trBe_bar_tr, dev_tau(6)
        REAL(prec) :: ff_bar(3, 3), F_bar_inv_n(3, 3), C_bar(3, 3)
        REAL(prec) :: J, kappa, mu, Y0, K, J_n, mu_bar, chi_tr
        REAL(prec) :: s_tr(3, 3), norm_s_tr, ret(6), retsq(6)
        REAL(prec) :: H(3, 3), p, I_mat(3, 3), dev_retsq(6)
        REAL(prec) :: Be_bar_n(3, 3), C_bar_n(3, 3), BB4(6, 6)
        REAL(prec) :: C(3, 3), tau_tr_v(6), beta3, beta4, BB3(6, 6)
        REAL(prec) :: s_tr_v(6), BB2(6, 6), beta0, beta1, beta2
        REAL(prec) :: c_strike_bar_tr(6, 6), xpp(6, 6), retoret(6, 6)
        REAL(prec) :: c_strike_ep_tr(6, 6), gamma, s_v(6), tau_v(6)
        REAL(prec) :: c_strike_ep(6, 6), ret_mat(3, 3),retsq_mat(3,3)
        REAL(prec) :: dev_retsq_mat(3, 3)
        !Decleration of constants
        REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
        REAL(prec), PARAMETER :: THREE=3.D0, FOUR= 4.D0, SIX=6.D0
        REAL(prec), PARAMETER :: NINE=9.D0


        !Declare matrix to voit mapping
        data m2v/ 1, 4, 5,
     1            4, 2, 6,
     2            5, 6, 3/

        !2nd Order Identity
        data I/ ONE, ONE, ONE, ZERO, ZERO, ZERO/
        
        !Define 4th order identity tensor
        data I_mat(1,:) /ONE, ZERO, ZERO/
        data I_mat(2,:) /ZERO, ONE, ZERO/
        data I_mat(3,:) /ZERO, ZERO, ONE/
      
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
        kappa=PROPS(1)
        mu=PROPS(2)
        Y0=PROPS(3)
        K=PROPS(4)
        
        write(*,*) 'dt : ', DTIME
        ! write(*,*) 'kappa : ', kappa
        ! write(*,*) 'mu : ', mu
        ! write(*,*) 'Y0 : ', Y0
        ! write(*,*) 'K : ', K
        
        
        !Calculate determinant of deformation gradient
        J = determinant(DFGRD1(:,:))
        J_n = determinant(DFGRD0(:,:))

        !Calculate inverse of deformation gradient
        CALL matInv(DFGRD0(:, :), F_n_inv(:, :))
        
        ff(:, :) = ZERO
        DO K1 = 1, 3
          DO ii = 1, 3
            DO jj = 1, 3
              ff(ii, jj) = ff(ii, jj) 
     1                   + DFGRD1(ii, K1) * F_n_inv(K1, jj)
            END DO
          END DO
        END DO

        ff_bar(:, :) = ff(:, :) 
     1               * (determinant(ff(:, :))) ** (-ONE / THREE)


        WRITE(*,*) 'F : '
        DO K1 = 1, 3
          WRITE(*,*) DFGRD1(K1, :)
        END DO

        WRITE(*,*) 'F_n : '
        DO K1 = 1, 3
          WRITE(*,*) DFGRD0(K1, :)
        END DO

        WRITE(*,*) 'J : ', J
        WRITE(*,*) ''

        WRITE(*,*) 'J_n : ', J_n
        WRITE(*,*) ''

        WRITE(*,*) 'F_n_inv : '
        DO K1 = 1, 3
          WRITE(*,*) F_n_inv(K1, :)
        END DO

        WRITE(*,*) 'ff : '
        DO K1 = 1, 3
          WRITE(*,*) ff(K1, :)
        END DO

        WRITE(*,*) 'ff_bar : '
        DO K1 = 1, 3
          WRITE(*,*) ff_bar(K1, :)
        END DO
        

        
        Be_bar_n(1, 1) = STATEV(1)
        Be_bar_n(2, 2) = STATEV(2)
        Be_bar_n(3, 3) = STATEV(3)
        Be_bar_n(1, 2) = STATEV(4)
        Be_bar_n(1, 3) = STATEV(5)
        Be_bar_n(2, 1) = STATEV(4)
        Be_bar_n(2, 3) = STATEV(6)
        Be_bar_n(3, 1) = STATEV(5)
        Be_bar_n(3, 2) = STATEV(6)

        alpha_n = STATEV(7)

        ! WRITE(*,*) 'Be_bar_n : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) Be_bar_n(K1, :)
        ! END DO
        ! WRITE(*,*) 'alpha_n : ', alpha_n
        
    
        
        !Calculate Be_bar_tr
        Be_bar_tr(:, :) = ZERO
        DO ii = 1, 3
          DO jj = 1, 3
            DO mm = 1, 3
              DO nn = 1, 3
                Be_bar_tr(ii, jj) = Be_bar_tr(ii, jj) 
     1          + ff_bar(ii, mm) * Be_bar_n(mm, nn) * ff_bar(jj, nn)
              END DO
            END DO
         END DO
        END DO

        ! WRITE(*,*) 'Be_bar_tr : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) Be_bar_tr(K1, :)
        ! END DO

        

        CALL dev(Be_bar_tr(:, :), dev_Be_bar_tr(:, :))

        trBe_bar_tr = Be_bar_tr(1, 1) + Be_bar_tr(2, 2) 
     1              + Be_bar_tr(3, 3)
      
        ! WRITE(*,*) 'dev_Be_bar_tr : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) dev_Be_bar_tr(K1, :)
        ! END DO

        ! Write(*,*) 'trBe_bar_tr : ', trBe_bar_tr
        
        
        
        !Calculate initial Kirchoff stress dev
        s_tr(:, :) = mu * dev_Be_bar_tr(:, :)

        DO K1 = 1, 3
          s_tr_v(K1) = s_tr(K1, K1)
        END DO
        s_tr_v(4) = s_tr(1, 2)
        s_tr_v(5) = s_tr(1, 3)
        s_tr_v(6) = s_tr(2, 3)
       
        ! WRITE(*,*) 's_tr : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) s_tr(K1, :)
        ! END DO

        ! WRITE(*,*) 's_tr_v : ', s_tr_v(:)
        

        norm_s_tr = ZERO
        DO ii = 1, 3
          DO jj = 1, 3
            norm_s_tr = norm_s_tr + s_tr(ii, jj) * s_tr(ii, jj)
          END DO
        END DO
        norm_s_tr = norm_s_tr ** (0.5D0)

        IF (norm_s_tr .EQ. ZERO) THEN
          ret(:) = ZERO
        ELSE
          ret(:) = s_tr_v(:) / norm_s_tr
        END IF

        ! WRITE(*,*) '||s_tr|| : ', norm_s_tr

        ! WRITE(*,*) 'n_n+1 : ', ret(:)
        
        p = 0.5D0 * kappa * (J - (ONE / J))

        tau_tr(:, :) = J * p * I_mat(:, :) + s_tr(:, :)

        DO K1 = 1, 3
          tau_tr_v(K1) = tau_tr(K1, K1)
        END DO
        tau_tr_v(4) = tau_tr(1, 2)
        tau_tr_v(5) = tau_tr(1, 3)
        tau_tr_v(6) = tau_tr(2, 3)

        mu_bar = mu * (ONE / THREE) * trBe_bar_tr
        beta0 = ONE + K / (THREE * mu_bar)

        ! write(*,*) 'p : ', p

        ! write(*,*) 'tau_tr: '
        ! DO K1 = 1, 3
        !   WRITE(*,*) tau_tr(K1, :)
        ! END DO

        ! write(*,*) 'tau_tr_v: ', tau_tr_v(:)

        !Algorithmic constants for the Tangent
        
        ret_mat(1, 1) = ret(1)
        ret_mat(2, 2) = ret(2)
        ret_mat(3, 3) = ret(1)
        ret_mat(1, 2) = ret(4)
        ret_mat(1, 3) = ret(5)
        ret_mat(2, 3) = ret(6)
        ret_mat(2, 1) = ret(4)
        ret_mat(3, 1) = ret(5)
        ret_mat(3, 2) = ret(6)

        retsq_mat(:,:) = ZERO
        DO K1 = 1, 3
          DO ii = 1, 3
            DO jj = 1, 3
                retsq_mat(ii, jj) = retsq_mat(ii, jj) 
     1                  + ret_mat(ii, K1) * ret_mat(K1, jj)
            END DO
          END DO
        END DO

        CALL dev(retsq_mat(:, :), dev_retsq_mat(:, :))
        
        dev_retsq(1) = dev_retsq_mat(1, 1)
        dev_retsq(2) = dev_retsq_mat(2, 2)
        dev_retsq(3) = dev_retsq_mat(3, 3)
        dev_retsq(4) = dev_retsq_mat(1, 2)
        dev_retsq(5) = dev_retsq_mat(1, 3)
        dev_retsq(6) = dev_retsq_mat(2, 3)

        ! write(*,*) 'retsq : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) retsq_mat(K1, :)
        ! END DO
        ! write(*,*) 'devretsq : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) dev_retsq_mat(K1, :)
        ! END DO
        ! write(*,*) 'devretsq_v : ', dev_retsq(:)

        
        
        BB1(:, :) = ZERO
        BB2(:, :) = ZERO
        retoret(:, :) = ZERO
        BB3(:, :) = ZERO
        DO ii = 1, 6
          DO jj = 1, 6
            BB1(ii, jj) = ret(ii) * I(jj)
            BB2(ii, jj) = I(ii) * ret(jj)
            retoret(ii, jj) = ret(ii) * ret(jj)
            BB3(ii, jj) = ret(ii) * dev_retsq(jj)
          END DO
        END DO

        BB4 = ZERO
        DO ii = 1, 6
          DO jj = 1, 6
            BB4(ii, jj) = BB4(ii, jj) 
     1                  + 0.5D0 * (BB3(ii, jj) + BB3(jj, ii))
          END DO
        END DO
        
        

        c_strike_bar_tr(:, :) = TWO * mu_bar * xpp(:, :) 
     1      - (TWO / THREE) * norm_s_tr * (BB1(:, :) + BB2(:, :))

        c_strike_tr(:, :) = kappa * (J ** TWO) * xioi(:, :) 
     2                    - kappa * ((J ** TWO) - ONE) * xii(:, :) 
     3                    + c_strike_bar_tr(:, :)

        chi_tr = norm_s_tr 
     1         - ((TWO / THREE) ** (0.5D0)) * (Y0 + K * alpha_n)


        ! write(*,*) 'chi_trial : ', chi_tr
        write(*,*) '----------------------------------'
        ! IF Yield surface isnt breached (Hyper Elastic Step)
        IF (chi_tr .LE. ZERO) THEN

          ! Return Cauchy stress 
          STRESS(:) = tau_tr_v(:) / J

          ! Calculate algorithmic constants for trial state

          
          beta3 = ONE / beta0
          beta4 = beta3 * norm_s_tr / mu_bar

          STATEV(1 : 6) = (s_tr_v(:) / mu) 
     1                  + (trBe_bar_tr/THREE) * I(:)

          c_strike_ep_tr(:, :) = c_strike_tr(:,:) 
     1                 - TWO * mu_bar * beta3 * retoret(:, :) 
     2                 - TWO * mu_bar * beta4 * BB4(:, :)

          !Tangent for ABAQUS
          DO ii = 1, 3
            DO jj = 1, 3
              DO ll = 1, 3
                DO mm = 1, 3
                  DDSDDE(m2v(ii, jj), m2v(ll, mm))   
     1        = DDSDDE(m2v(ii, jj), m2v(ll, mm)) 
     2        + (ONE / J) * (c_strike_ep_tr(m2v(ii,jj), m2v(ll, mm)) 
     3        + I(m2v(ii, ll)) * tau_tr_v(m2v(jj, mm)) 
     4        + I(m2v(jj, mm)) * tau_tr_v(m2v(ii, ll)))
                END DO
              END DO
            END DO
          END DO
            

        ELSE

          gamma = (chi_tr / (TWO * mu_bar)) / beta0

          s_v = s_tr_v - TWO * mu_bar * gamma * ret(:)

          STATEV(7) = STATEV(7) + ((TWO/THREE) ** (0.5D0)) * gamma

          tau_v(:) = J * p * I(:) + s_v(:)

          STRESS(:) =  tau_v(:) / J

          STATEV(1 : 6) = (s_v(:) / mu) 
     1                  + (trBe_bar_tr / THREE) * I(:)


          beta1 = (TWO * mu_bar * gamma) / norm_s_tr
          beta2 = (ONE - (ONE / beta0)) 
     1          * (TWO / THREE) * (norm_s_tr / mu_bar) * gamma
          beta3 = (ONE / beta0) - beta1 + beta2
          beta4 = ((ONE / beta0) - beta1) * norm_s_tr / mu_bar

          c_strike_ep(:, :) = c_strike_tr(:, :) 
     1                      - beta1 * c_strike_bar_tr(:, :) 
     2                      - TWO * mu_bar * beta3 * retoret(:, :) 
     3                      - TWO * mu_bar * beta4 * BB4(:, :)

                !Tangent for ABAQUS
          DO ii = 1, 3
            DO jj = 1, 3
              DO ll = 1, 3
                DO mm = 1, 3
                  DDSDDE(m2v(ii, jj), m2v(ll, mm))   
     1        = DDSDDE(m2v(ii, jj), m2v(ll, mm)) 
     2        + (ONE / J) * (c_strike_ep(m2v(ii,jj), m2v(ll, mm)) 
     3        + I(m2v(ii, ll)) * tau_v(m2v(jj, mm)) 
     4        + I(m2v(jj, mm)) * tau_v(m2v(ii, ll)))
                END DO
              END DO
            END DO
          END DO
          

        END IF
     
      RETURN
      END SUBROUTINE FinJ2

      !--------------------------------------------------------------
      !     Helper function to compute Determinant of 3x3 matrix 
      !--------------------------------------------------------------
      FUNCTION determinant(matrix)
        IMPLICIT NONE
                
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec):: determinant

        determinant = matrix(1,1) * (matrix(2,2) * matrix(3,3) 
     1   - matrix(3,2) * matrix(2, 3)) - matrix(1,2) * (matrix(2,1)
     2   * matrix(3,3) - matrix(2,3) * matrix(3,1)) + matrix(1,3) 
     3   * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))

      END FUNCTION determinant


      !--------------------------------------------------------------
      !      Helper function to compute inverse of 3x3 matrix
      !--------------------------------------------------------------
      SUBROUTINE matInv(matrix, inv)
        IMPLICIT NONE
        
        INTERFACE
          FUNCTION determinant(mat)
            INTEGER, PARAMETER :: double=kind(1.d0)
            INTEGER, PARAMETER :: prec=double

            REAL(prec), INTENT(IN) :: mat(3, 3)
            REAL(prec) :: determinant
          END FUNCTION determinant
        END INTERFACE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec), INTENT(OUT) :: inv(3, 3)

        inv(1,1) = matrix(2,2)*matrix(3,3) - matrix(3,2)*matrix(2,3)
        inv(1,2) = matrix(3,2)*matrix(1,3) - matrix(1,2)*matrix(3,3)
        inv(1,3) = matrix(1,2)*matrix(2,3) - matrix(2,2)*matrix(1,3)
        inv(2,1) = matrix(3,1)*matrix(2,3) - matrix(2,1)*matrix(3,3)
        inv(2,2) = matrix(1,1)*matrix(3,3) - matrix(3,1)*matrix(1,3)
        inv(2,3) = matrix(2,1)*matrix(1,3) - matrix(1,1)*matrix(2,3)
        inv(3,1) = matrix(2,1)*matrix(3,2) - matrix(3,1)*matrix(2,2)
        inv(3,2) = matrix(3,1)*matrix(1,2) - matrix(1,1)*matrix(3,2)
        inv(3,3) = matrix(1,1)*matrix(2,2) - matrix(2,1)*matrix(1,2)
      
        inv(:, :) = inv(:, :)/determinant(matrix(:, :))

        RETURN
      END SUBROUTINE matInv

      !--------------------------------------------------------------
      !      Helper function to compute dev of voit array
      !--------------------------------------------------------------
      SUBROUTINE devVoit(array, dev)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: array(6)
        REAL(prec), INTENT(OUT) :: dev(6)

        dev(1 : 3) = array(1 : 3) 
     1             - (1.D0/3.D0) * (array(1) + array(2) + array(3))

        dev(4 : 6) = array(4 : 6)
        RETURN
      END SUBROUTINE devVoit

      !--------------------------------------------------------------
      !      Helper function to compute dev of voit array
      !--------------------------------------------------------------
      SUBROUTINE dev(matrix, out)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec), INTENT(OUT) :: out(3, 3)
        REAL(prec) :: I_mat(3, 3), tr
        INTEGER :: K1

        !Define 2nd order identity tensor
        data I_mat(1,:) /1.D0, 0.D0, 0.D0/
        data I_mat(2,:) /0.D0, 1.D0, 0.D0/
        data I_mat(3,:) /0.D0, 0.D0, 1.D0/

        ! write(*,*) '-------- start of dev ----------'

        ! WRITE(*,*) 'matrix : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) matrix(K1, :)
        ! END DO

        ! WRITE(*,*) 'I : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) I_mat(K1, :)
        ! END DO

        tr = matrix(1, 1) + matrix(2, 2) + matrix(3, 3)

        ! WRITE(*,*) 'trace : ', tr

        out(:, :) = matrix(:, :) 
     1            - (1.D0/3.D0) * tr * I_mat(:, :)

        ! WRITE(*,*) 'dev : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) out(K1, :)
        ! END DO

        ! write(*,*) '-------- start of dev ----------'

        RETURN
      END SUBROUTINE dev



      !--------------------------------------------------------------
      !      Helper function to compute DEV of voit array
      !--------------------------------------------------------------
      SUBROUTINE ddevVoit(array, C, dev)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: array(6), C(6)
        REAL(prec), INTENT(OUT) :: dev(6)
        REAL(prec) :: C_mat(3, 3), C_inv_mat(3, 3), m2v(3, 3)
        REAL(prec) :: C_inv(6), Ctr
        INTEGER :: K1, K2

        !Declare matrix to voit mapping
        data m2v/ 1, 4, 5,
     1            4, 2, 6,
     2            5, 6, 3/


        C_mat(1, 1) = C(1)
        C_mat(2, 2) = C(2)
        C_mat(3, 3) = C(3)
        C_mat(1, 2) = C(4)
        C_mat(1, 3) = C(5)
        C_mat(2, 1) = C(4)
        C_mat(2, 3) = C(6)
        C_mat(3, 1) = C(5)
        C_mat(3, 2) = C(6)

        CALL matInv(C_mat, C_inv_mat)

        DO K1 = 1, 3
          DO K2 = 1, 3
            C_inv(m2v(K1, K2)) = C_inv_mat(K1, K2)
          END DO
        END DO

        Ctr = 0.D0
        DO K1 = 1, 3
          DO K2 = 1, 3
            Ctr = Ctr + array(m2v(K1, K2)) * C(m2v(K1, K2))
          END DO
        END DO

        dev(:) = array(:) - (1.D0 / 3.D0) * Ctr * C_inv(:)
        RETURN
      END SUBROUTINE ddevVoit

      !--------------------------------------------------------------
      !      Helper function to compute DEV of matrix
      !--------------------------------------------------------------
      SUBROUTINE ddev(matrix, C, dev)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3), C(3, 3)
        REAL(prec), INTENT(OUT) :: dev(3, 3)
        REAL(prec) :: C_mat(3, 3), C_inv(3, 3), Ctr
        INTEGER :: ii, jj, kk, ll

        ! WRITE(*,*) '------------ ddev starts here-------------'

        CALL matInv(C, C_inv)
        
        ! WRITE(*,*) 'matrix : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) matrix(K1, :)
        ! END DO

        ! WRITE(*,*) 'C : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) C(K1, :)
        ! END DO

        ! WRITE(*,*) 'C_inv : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) C_inv(K1, :)
        ! END DO
        
        Ctr = 0.D0
        DO ii = 1, 3
          DO jj = 1, 3
            Ctr = Ctr + matrix(ii, jj) * C(ii, jj)
          END DO
        END DO 

        ! write(*,*) 'Ctr : ', Ctr

        dev(:, :) = matrix(:, :) - (1.D0 / 3.D0) * Ctr * C_inv(:, :)

        ! WRITE(*,*) 'dev : '
        ! DO K1 = 1, 3
        !   WRITE(*,*) dev(K1, :)
        ! END DO

        ! WRITE(*,*) '------------ ddev ends here-------------'

        RETURN
      END SUBROUTINE ddev