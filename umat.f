C     UMAT - Isotropic Hyperelastic from Computational inealsticity
C     Goal is to implement just the initail stored energy formulation
C     Later it will be modified to incorporate Kinematically non linear VE

C     Mohib Mustafa - IMDEA 14 JULY 2021
            
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
      
      call VE(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &         NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1)

      return
      end

C     !--------------------------------------------------------------
C     !   UMAT SUBROUTINE FOR ISOTROPIC NEO HOOKIAN MATERIAL MODEL
C     !--------------------------------------------------------------

      SUBROUTINE VE(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,
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
        INTEGER    :: index, ii, jj, mm, ll, K1, K2
        REAL(prec) :: I(6), m2v(3, 3), B_bar(3, 3), tau(3, 3)
        REAL(prec) :: F_inv(3, 3), F_bar(3, 3), BB1(6, 6), C_n(3, 3)
        REAL(prec) :: c_strike(6, 6), xioi(6, 6), xii(6, 6)
        REAL(prec) :: c_strike_bar(6, 6), F_bar_inv(3, 3), g_star
        REAL(prec) :: dev_B_bar(3, 3), trB_bar, dev_tau(6), J_n
        REAL(prec) :: F_bar_n(3, 3), F_bar_inv_n(3, 3), C_bar(3, 3)
        REAL(prec) :: J, kappa, mu0, mu1, eta1, muinf, ginf, g1, tau1
        REAL(prec) :: tau_o_bar(3, 3), S_til_o_n(3, 3), H_til_n(3, 3)
        REAL(prec) :: S_til_o(3, 3), H(3, 3), p, I_mat(3, 3)
        REAL(prec) :: H_n(3, 3), C_bar_n(3, 3), D(3, 3), xpp(6, 6)
        REAL(prec) :: hh_bar_n(3, 3), h_bar_n, C(3, 3), tau_v(6)
        REAL(prec) :: tau_o_bar_v(6), BB2(6, 6), BB3(6, 6), BB4(6, 6)
        REAL(prec) :: hh_bar_n_v(6), c_strike_o_bar(6, 6)
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
        mu0=PROPS(2)
        mu1=PROPS(3)
        eta1=PROPS(4)
        
        write(*,*) 'dt : ', DTIME
        ! write(*,*) 'kappa : ', kappa
        ! write(*,*) 'mu0 : ', mu0
        ! write(*,*) 'mu1 : ', mu1
        ! write(*,*) 'eta1 : ', eta1

        !Calculate other derived material parameters
        muinf=mu0 - mu1
        ginf=muinf/mu0
        g1=mu1/mu0
        tau1=eta1/mu1
        g_star = ginf + g1 * EXP((-DTIME/(TWO * tau1))) 

        ! write(*,*) 'muinf : ', muinf
        ! write(*,*) 'ginf : ', ginf
        ! write(*,*) 'g1 : ', g1
        ! write(*,*) 'tau1 : ', tau1
        ! write(*,*) 'g_star : ', g_star

        !Calculate determinant of deformation gradient
        J = determinant(DFGRD1(:,:))
        J_n = determinant(DFGRD0(:,:))

        !Calculate inverse of deformation gradient
        !CALL matInv(DFGRD1(:,:), F_inv(:, :))
        
        F_bar(:, :) = (J ** (-ONE / THREE)) * DFGRD1(:, :)
        CALL matInv(F_bar(:, :), F_bar_inv(:, :))

        F_bar_n(:, :) = (J_n ** (-ONE / THREE)) * DFGRD0(:, :)
        !CALL matInv(F_bar_n(:, :), F_bar_inv_n(:, :))

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

        WRITE(*,*) 'F_bar : '
        DO K1 = 1, 3
          WRITE(*,*) F_bar(K1, :)
        END DO

        WRITE(*,*) 'F_bar_n : '
        DO K1 = 1, 3
          WRITE(*,*) F_bar_n(K1, :)
        END DO


        H_n(1, 1) = STATEV(1)
        H_n(2, 2) = STATEV(2)
        H_n(3, 3) = STATEV(3)
        H_n(1, 2) = STATEV(4)
        H_n(1, 3) = STATEV(5)
        H_n(2, 1) = STATEV(4)
        H_n(2, 3) = STATEV(6)
        H_n(3, 1) = STATEV(5)
        H_n(3, 2) = STATEV(6)

        WRITE(*,*) 'H_n : '
        DO K1 = 1, 3
          WRITE(*,*) H_n(K1, :)
        END DO
        

    !     !Calculate C
    !     C(:, :) = ZERO
    !     DO K1 = 1, 3
    !       DO ii = 1, 3
    !         DO jj = 1, 3
    !             C(ii, jj) = C(ii, jj) 
    !  1                    + DFGRD1(K1, ii) * DFGRD1(K1, jj)
    !         END DO
    !      END DO
    !     END DO

        !Calculate C_n
        C_n(:, :) = ZERO
        DO K1 = 1, 3
          DO ii = 1, 3
            DO jj = 1, 3
                C_n(ii, jj) = C_n(ii, jj) 
     1                           + DFGRD0(K1, ii) * DFGRD0(K1, jj)
            END DO
         END DO
        END DO

        WRITE(*,*) 'C_n : '
        DO K1 = 1, 3
          WRITE(*,*) C_n(K1, :)
        END DO



    !     !Calculate C_bar
    !     C_bar(:, :) = ZERO
    !     DO K1 = 1, 3
    !       DO ii = 1, 3
    !         DO jj = 1, 3
    !             C_bar(ii, jj) = C_bar(ii, jj) 
    !  1                             + F_bar(K1, ii) * F_bar(K1, jj)
    !         END DO
    !      END DO
    !     END DO

    !     !Calculate C_bar_n
    !     C_bar_n(:, :) = ZERO
    !     DO K1 = 1, 3
    !       DO ii = 1, 3
    !         DO jj = 1, 3
    !             C_bar_n(ii, jj) = C_bar_n(ii, jj) 
    !  1                            + F_bar_n(K1, ii) * F_bar_n(K1, jj)
    !         END DO
    !      END DO
    !     END DO
        
        !Calculate B_bar
        B_bar(:, :) = ZERO
        DO K1 = 1, 3
          DO ii = 1, 3
            DO jj = 1, 3
                B_bar(ii, jj) = B_bar(ii, jj) 
     1                             + F_bar(ii, K1) * F_bar(jj, K1)
            END DO
         END DO
        END DO

        WRITE(*,*) 'B_bar : '
        DO K1 = 1, 3
          WRITE(*,*) B_bar(K1, :)
        END DO


        CALL dev(B_bar, dev_B_bar)

        trB_bar = B_bar(1, 1) + B_bar(2, 2) + B_bar(3, 3)
      
        WRITE(*,*) 'dev_B_bar : '
        DO K1 = 1, 3
          WRITE(*,*) dev_B_bar(K1, :)
        END DO

        Write(*,*) 'tr[B_bar] : ', trB_bar
        
        !Calculate initial Kirchoff stress dev
        tau_o_bar(:, :) = mu0 * dev_B_bar(:, :)

        DO K1 = 1, 3
          tau_o_bar_v(K1) = tau_o_bar(K1, K1)
        END DO
        tau_o_bar_v(4) = tau_o_bar(1, 2)
        tau_o_bar_v(5) = tau_o_bar(1, 3)
        tau_o_bar_v(6) = tau_o_bar(2, 3)
       
        WRITE(*,*) 'tau_o_bar : '
        DO K1 = 1, 3
          WRITE(*,*) tau_o_bar(K1, :)
        END DO

        WRITE(*,*) 'tau_o_bar_v : ', tau_o_bar_v(:)

        CALL ddev(mu0 * I_mat(:, :), C_n(:, :), S_til_o_n(:, :))
        
        WRITE(*,*) 'S_til_o_n : '
        DO K1 = 1, 3
          WRITE(*,*) S_til_o_n(K1, :)
        END DO
        
        H_til_n(:, :) = EXP(-DTIME / tau1) * H_n(:, :) 
     1             - EXP(-DTIME / (TWO * tau1)) * S_til_o_n(:, :)

        WRITE(*,*) 'H_til_n : '
        DO K1 = 1, 3
          WRITE(*,*) H_til_n(K1, :)
        END DO

         
        !Calculate S_til_o
        S_til_o(:, :) = ZERO
        DO K1 = 1, 3
          DO K2 = 1, 3
            DO ii = 1, 3
              DO jj = 1, 3

                  S_til_o(ii, jj) = S_til_o(ii, jj) 
     1             + F_bar_inv(ii, K1) 
     2             * tau_o_bar(K1, K2) * F_bar_inv(jj, K2)

              END DO
            END DO
          END DO
        END DO

        WRITE(*,*) 'S_til_o : '
        DO K1 = 1, 3
          WRITE(*,*) S_til_o(K1, :)
        END DO

        
        
        H(:, :) = H_til_n(:, :) 
     1          + EXP(-DTIME/(TWO * tau1)) * S_til_o(:, :)

        DO K1 = 1, 3
          STATEV(K1) = H(K1, K1)
        END DO
        STATEV(4) = H(1, 2)
        STATEV(5) = H(1, 3)
        STATEV(6) = H(2, 3)

        WRITE(*,*) 'H : '
        DO K1 = 1, 3
          WRITE(*,*) H(K1, :)
        END DO

        WRITE(*,*) 'STAT : ', STATEV(1 : 6)
        write(*,*) '------------------------------------------------'
        p = 0.5D0 * kappa * (J - (ONE / J))

        

        !CALL devVoit(tau, dev_tau)
        
        D(:, :) = ZERO
        DO K1 = 1, 3
          DO K2 = 1, 3
            DO ii = 1, 3
              DO jj = 1, 3
                D(ii, jj) = D(ii, jj) 
     1                    + F_bar(ii, K1) * H_til_n(K1, K2) 
     2                    * F_bar(jj, K2)
              END DO
            END DO
          END DO
        END DO
        
        
        
        CALL dev(g1 * D(:, :), hh_bar_n(:, :))

        DO K1 = 1, 3
          hh_bar_n_v(K1) = hh_bar_n(K1, K1)
        END DO
        hh_bar_n_v(4) = hh_bar_n(1, 2)
        hh_bar_n_v(5) = hh_bar_n(1, 3)
        hh_bar_n_v(6) = hh_bar_n(2, 3)

        h_bar_n = g1 * (D(1, 1) + D(2, 2) + D(3, 3))

        tau(:, :) = J * p * I_mat(:, :) 
     1            + g_star * tau_o_bar(:, :) + hh_bar_n(:, :)

        
        
        DO K1 = 1, 3
          tau_v(K1) = tau(K1, K1)
        END DO
        tau_v(4) = tau(1, 2)
        tau_v(5) = tau(1, 3)
        tau_v(6) = tau(2, 3)
        
        !Return Cauchy stress to Abaqus
        STRESS(:) = (ONE / J) * tau_v(:)
        
        !Algorithmic constants for the Tangent
        BB1(:, :) = ZERO
        DO ii = 1, 6
          DO jj = 1, 6
            BB1(ii, jj) = tau_o_bar_v(ii) * I(jj)
            BB2(ii, jj) = I(ii) * tau_o_bar_v(jj)
            BB3(ii, jj) = hh_bar_n_v(ii) * I(jj)
            BB4(ii, jj) = I(ii) * hh_bar_n_v(jj)
          END DO
        END DO

        c_strike_o_bar(:, :) = (-TWO / THREE) 
     1     * (BB1(:, :) + BB2(:, :)) 
     2     + (TWO / THREE) * mu0 * (B_bar(1, 1) + B_bar(2, 2) 
     3     + B_bar(3, 3)) * xpp(:, :)

        c_strike_bar(:, :) = g_star * c_strike_o_bar(:, :) 
     1                     - (TWO / THREE) * (BB3(:, :) + BB4(:, :)) 
     2                     + (TWO / THREE) * h_bar_n * xpp(:, :)

        ! Tangennt in current config
        c_strike(:, :) = c_strike_bar(:, :) 
     1                 + J * p * (xioi(:, :) - TWO * xii(:, :)) 
     2               + 0.5D0 * kappa * (J**TWO + ONE) * xioi(:, :)  
       
        !Tangent for ABAQUS
        DO ii = 1, 3
          DO jj = 1, 3
            DO ll = 1, 3
              DO mm = 1, 3
                DDSDDE(m2v(ii, jj), m2v(ll, mm))   
     1        = DDSDDE(m2v(ii, jj), m2v(ll, mm)) 
     2        + (ONE / J) * (c_strike(m2v(ii,jj), m2v(ll, mm)) 
     3        + I(m2v(ii, ll)) * tau_v(m2v(mm, jj)) 
     4        + I(m2v(jj, mm)) * tau_v(m2v(ii, ll)))
              END DO
            END DO
          END DO
        END DO

        

      RETURN
      END SUBROUTINE VE

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