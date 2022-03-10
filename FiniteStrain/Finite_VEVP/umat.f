C     UMAT - Finite strain VEVP
C     Based on Paper by Nguyen and Noels

C     Note : Computations are done in standard matrix notation and
C            at the end changed to voit where needed for return to main 
C            routine. I understand this isnt the most elegant way of 
C            doing things but when dealing with 2nd order tensors 
C            with 4 or more dummy indices, I feel I can code much 
C            faster this way. If you know the algebra and have the 
C            time to work in voit notation change it and send me an 
C            email at mohibmustafa@gmail.com
C
C

C     Note : All consistent tangents are calculated numerically. 

C     Mohib Mustafa - IMDEA 1 Dec 2021
            
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
      
      call VEVP(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &         NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1)

      return
      end

C     !--------------------------------------------------------------
C     !                     UMAT SUBROUTINE
C     !--------------------------------------------------------------

      SUBROUTINE VEVP(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,
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
        INTEGER    :: ii, jj, mm, ll, O1, O5, order
        REAL(prec) :: tau_tr(3,3), J_n, J, KK_inf, KK_1
        REAL(prec) :: k_1, GG_inf, GG_1, g_1, KK_e, GG_e
        REAL(prec) :: E_ve_tr(3,3), E_ve_n(3,3), I_mat(3,3), AA(3,3)
        REAL(prec) :: A_1_n(3,3), F_hat(6,3,3), J_hat(6), beta
        REAL(prec) :: F_vp_n(3,3), B_1_n, C_ve_tr(3,3),kappa_tr(3,3)
        REAL(prec) :: A_1(3,3), F_ve_tr(3, 3)
        REAL(prec) :: F_ve_tr_hat(6, 3, 3), C_ve_tr_hat(6, 3, 3)
        REAL(prec) :: E_ve_tr_hat(6, 3, 3)
        REAL(prec) :: kappa_tr_hat(6, 3, 3)
        REAL(prec) :: AA_hat(6, 3, 3), tr_phi_tr_hat(6)
        REAL(prec) :: tau_tr_hat(6, 3, 3), tau_tr_hat_v(6, 6)
        REAL(prec) :: S_tr(3, 3), S_tr_hat(6, 3, 3)
        REAL(prec) :: logAA(3, 3)
        REAL(prec) :: dlogAA(3, 3, 3, 3), ddlogAA(3, 3, 3, 3, 3, 3)
        REAL(prec) :: II_mat(3,3,3,3), temp1(3,3,3,3), temp2(3,3,3,3)
        REAL(prec) :: B_1
        REAL(prec) :: logAA_hat(6, 3, 3), dHHbdgma
        REAL(prec) :: dlogAA_hat(6,3,3,3,3)
        REAL(prec) :: ddlogAA_hat(6, 3, 3, 3, 3, 3, 3)
        REAL(prec) :: A_1_hat(6, 3, 3)
        REAL(prec) :: B_1_hat(6), phi_tr(3, 3), b_tr(3, 3), phi_p_tr
        REAL(prec) :: dev_phi_tr(3, 3), phi_e_tr, sigma_t, sigma_c
        REAL(prec) :: m, a0, a1, a2, F_tr, alpha, nu_p
        REAL(prec) :: eta, p_exp, k, sigma_c0, h_c1, h_c2, h_b1
        REAL(prec) :: h_b2, h_b0, HHc, h_cexp
        REAL(prec) :: dev_phi(3, 3), Q(3, 3), tr_Q, dev_Q(3, 3)
        REAL(prec) :: kappa(3, 3), exp_GQ(3, 3)
        REAL(prec) :: gma_n, u, v, tr_phi_tr, HHb
        REAL(prec) :: F_vp(3, 3), C_ve(3, 3)
        REAL(prec) :: AA_upd(3, 3), logAA_upd(3, 3), S(3, 3)
        REAL(prec) :: dlogAA_upd(3, 3, 3, 3), E_ve(3, 3), F_ve(3, 3)
        REAL(prec) :: ddlogAA_upd(3, 3, 3, 3, 3, 3), tau(3, 3)
        REAL(prec) :: b(3, 3), phi_tr_hat(6, 3, 3)
        REAL(prec) :: phi_p_tr_hat(6), dev_phi_tr_hat(6, 3, 3)
        REAL(prec) :: phi_e_tr_hat(6), F_tr_hat(6), u_hat(6),v_hat(6)
        REAL(prec) :: tr_Q_hat(6)
        REAL(prec) :: dev_phi_hat(6, 3, 3), Q_hat(6, 3, 3)
        REAL(prec) :: dev_Q_hat(6, 3, 3), exp_GQ_hat(6,3,3)
        REAL(prec) :: kappa_hat(6, 3, 3)
        REAL(prec) :: F_vp_hat(6, 3, 3)
        REAL(prec) :: C_ve_hat(6, 3, 3), AA_upd_hat(6,3,3)
        REAL(prec) :: logAA_upd_hat(6, 3, 3), E_ve_hat(6,3,3)
        REAL(prec) :: dlogAA_upd_hat(6,3,3,3,3), S_hat(6,3,3)
        REAL(prec) :: ddlogAA_upd_hat(6,3,3,3,3,3,3), b_e
        REAL(prec) :: F_ve_hat(6,3,3)
        REAL(prec) :: tau_hat(6,3,3), tau_hat_v(6, 6)
        REAL(prec) :: temp6(3, 3), temp6_hat(6, 3, 3)
        REAL(prec) :: GAMMA, GG_til, KK_til, A, HHt
        REAL(prec) :: dDgammaDGamma, DfDGamma
        REAL(prec) :: ptilde, PhiEq
        REAL(prec) :: F_vp_inv(3, 3), ptilde_hat(6), PhiEq_hat(6)
        REAL(prec) :: GAMMA_hat(6), A_hat(6), F_vp_inv_hat(6, 3, 3)

        !Decleration of constants
        REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
        REAL(prec), PARAMETER :: THREE=3.D0, FOUR= 4.D0, SIX=6.D0
        REAL(prec), PARAMETER :: NINE=9.D0, TOLL=0.001D0    ! Be carefull this toll is used by perturb_F Subroutine
        REAL(prec), PARAMETER :: TOLL_N=0.0000001D0         ! This toll is used by newton raphson routine
        REAL(prec), PARAMETER :: TOLL_G=0.999D-7
        !Define 2nd order identity tensor
        data I_mat(1,:) /ONE, ZERO, ZERO/
        data I_mat(2,:) /ZERO, ONE, ZERO/
        data I_mat(3,:) /ZERO, ZERO, ONE/

        ! Define 4th order symetric identity tensor 
        II_mat(:,:,:,:) = 0.D0
        temp1(:,:,:,:) = 0.D0
        temp2(:,:,:,:) = 0.D0
        DO ii = 1, 3
          DO jj = 1, 3
            DO mm = 1, 3
              DO ll = 1, 3
                temp1(ii,jj,mm,ll) = temp1(ii,jj,mm,ll) 
     1                             + I_mat(ii,mm) * I_mat(jj,ll)
                temp2(ii,jj,mm,ll) = temp2(ii,jj,mm,ll) 
     1                             + I_mat(ii,ll) * I_mat(jj,mm)
              END DO
            END DO
          END DO
        END DO
        II_mat(:,:,:,:) = 0.5D0 * (temp1(:,:,:,:) + temp2(:,:,:,:))

        ! State Variables at previous time step

        ! This is a workaround to set the diagnol values of C_vp = 1
        ! in FFTMAD. Can be removed for ABAQUS, IN abaqus use the inp
        ! file to set initail values for C_vp
        IF ((KINC .EQ. 1) .AND. (KSTEP .EQ. 1)) THEN
          STATEV(1) = ONE
          STATEV(2) = ONE
          STATEV(3) = ONE
        END IF

        F_vp_n(1, 1) = STATEV(1)
        F_vp_n(2, 2) = STATEV(2)
        F_vp_n(3, 3) = STATEV(3)
        F_vp_n(1, 2) = STATEV(4)
        F_vp_n(1, 3) = STATEV(5)
        F_vp_n(2, 3) = STATEV(6)
        F_vp_n(2, 1) = STATEV(7)
        F_vp_n(3, 1) = STATEV(8)
        F_vp_n(3, 2) = STATEV(9)

        B_1_n = STATEV(10)

        A_1_n(1, 1) = STATEV(11)
        A_1_n(2, 2) = STATEV(12)
        A_1_n(3, 3) = STATEV(13)
        A_1_n(1, 2) = STATEV(14)
        A_1_n(1, 3) = STATEV(15)
        A_1_n(2, 3) = STATEV(16)
        A_1_n(2, 1) = STATEV(17)
        A_1_n(3, 1) = STATEV(18)
        A_1_n(3, 2) = STATEV(19)

        E_ve_n(1, 1) = STATEV(20)
        E_ve_n(2, 2) = STATEV(21)
        E_ve_n(3, 3) = STATEV(22)
        E_ve_n(1, 2) = STATEV(23)
        E_ve_n(1, 3) = STATEV(24)
        E_ve_n(2, 3) = STATEV(25)
        E_ve_n(2, 1) = STATEV(26)
        E_ve_n(3, 1) = STATEV(27)
        E_ve_n(3, 2) = STATEV(28)

        b_tr(1, 1) = STATEV(29)
        b_tr(2, 2) = STATEV(30)
        b_tr(3, 3) = STATEV(31)
        b_tr(1, 2) = STATEV(32)
        b_tr(1, 3) = STATEV(33)
        b_tr(2, 3) = STATEV(34)
        b_tr(2, 1) = STATEV(35)
        b_tr(3, 1) = STATEV(36)
        b_tr(3, 2) = STATEV(37)

        gma_n = STATEV(38)

        !Get material properties
        KK_inf   =PROPS(1)
        KK_1     =PROPS(2)
        k_1      =PROPS(3)
        GG_inf   =PROPS(4)
        GG_1     =PROPS(5)
        g_1      =PROPS(6)
        order    =PROPS(7)
        m        =PROPS(8)
        alpha    =PROPS(9)
        nu_p     =PROPS(10)
        eta      =PROPS(11)
        p_exp    =PROPS(12)
        sigma_c0 =PROPS(13)
        h_c1     =PROPS(14)
        h_c2     =PROPS(15)
        h_cexp   =PROPS(16)
        h_b0     =PROPS(17)
        h_b1     =PROPS(18)
        h_b2     =PROPS(19)

        !Calculate other derived material parameters
        KK_e = KK_inf + KK_1 * EXP(-DTIME / (TWO * k_1))
        GG_e = GG_inf + GG_1 * EXP(-DTIME / (TWO * g_1))

        beta = (9.D0 / 2.D0) 
     1         * ((1.D0 - 2.D0 * nu_p) / (nu_p + 1.D0))

        k = 1.D0 / (1.D0 + 2.D0 * nu_p ** 2.D0)**(0.5D0)
        
        CALL perturb_F(DFGRD1(:, :), F_hat(:, :, :))

        !Calculate determinant of deformation gradient
        J = determinant(DFGRD1(:,:))
        J_n = determinant(DFGRD0(:,:))

        J_hat(:) = ZERO
        DO O1 = 1, 6
          J_hat(O1) = determinant(F_hat(O1, :, :))
        END DO
     
         CALL vevpSplit(DFGRD1(:, :), F_vp_n(:, :),
     1    F_ve_tr(:, :), C_ve_tr(:, :))

        
        DO O5 = 1, 6
          CALL vevpSplit(F_hat(O5, :, :), F_vp_n(:, :),
     1    F_ve_tr_hat(O5, :, :), C_ve_tr_hat(O5, :, :))
        END DO

       AA(:, :) = C_ve_tr(:, :) - I_mat(:, :)
       CALL approx_log(AA, order, logAA, dlogAA, ddlogAA)
       E_ve_tr(:, :) = 0.5D0 * logAA(:, :)
       
       DO O5 = 1, 6
        AA_hat(O5, :, :) = C_ve_tr_hat(O5, :, :) - I_mat(:, :)
        CALL approx_log(AA_hat(O5,:,:), order, logAA_hat(O5,:,:), 
     1   dlogAA_hat(O5,:,:,:,:), ddlogAA_hat(O5, :, :, :, :, :, :))
        E_ve_tr_hat(O5, :, :) = 0.5D0 * logAA_hat(O5, :, :)
       END DO

       CALL vePredictor_log(E_ve_tr(:,:), E_ve_n(:, :), A_1_n(:,:),
     1   B_1_n, DTIME, GG_inf, GG_1, g_1, KK_inf, KK_1, k_1,
     2   kappa_tr(:, :), A_1(:, :), B_1)

       DO O5  = 1, 6
        CALL vePredictor_log(E_ve_tr_hat(O5,:,:), E_ve_n(:, :)
     1   , A_1_n(:,:), B_1_n, DTIME, GG_inf, GG_1, g_1, KK_inf, KK_1
     2   , k_1, kappa_tr_hat(O5, :, :), A_1_hat(O5, :, :)
     3   , B_1_hat(O5))
       END DO

       

       ! Calculate S as per paper
       CALL mat24ddot(kappa_tr(:, :), dlogAA(:, :, :, :), S_tr(:, :))

       DO O5 = 1, 6
        CALL mat24ddot(kappa_tr_hat(O5, :, :),
     1            dlogAA_hat(O5, :, :, :, :), S_tr_hat(O5, :, :))
       END DO

    !    ! Calculate kirchoff stress
       CALL mat2dot(F_ve_tr(:, :), S_tr(:, :), temp6(:, :))
       CALL mat2dotT(temp6(:, :), F_ve_tr(:, :), tau_tr(:, :))

      DO O5 = 1, 6
       CALL mat2dot(F_ve_tr_hat(O5,:, :), S_tr_hat(O5,:, :),
     1   temp6_hat(O5, :, :))
       CALL mat2dotT(temp6_hat(O5, :, :), F_ve_tr_hat(O5, :, :),
     1   tau_tr_hat(O5, :, :))
      END DO

       ! Calculate yield function at VE trial state

        ! Get phi_e_tr, phi_p_tr
        phi_tr(:, :) = kappa_tr(:,:) - b_tr(:,:)
        
        CALL tr_dev_split(phi_tr(:,:), dev_phi_tr(:,:), tr_phi_tr)
        phi_p_tr = tr_phi_tr / 3.D0

        CALL mat2ddot(dev_phi_tr(:, :), dev_phi_tr(:, :), phi_e_tr)
        phi_e_tr = (3.D0 / 2.D0) * phi_e_tr
        phi_e_tr = SQRT(phi_e_tr)


        ptilde = phi_p_tr
        PhiEq = phi_e_tr

        GAMMA = 0.0D0


        DO O5 = 1, 6
          phi_tr_hat(O5, :, :) = kappa_tr_hat(O5, :,:) - b_tr(:,:)
        
          CALL tr_dev_split(phi_tr_hat(O5, :,:), 
     1     dev_phi_tr_hat(O5,:,:), tr_phi_tr_hat(O5))
          phi_p_tr_hat(O5) = tr_phi_tr_hat(O5) / 3.D0
  
          CALL mat2ddot(dev_phi_tr_hat(O5, :, :),
     1      dev_phi_tr_hat(O5, :, :), phi_e_tr_hat(O5))
          phi_e_tr_hat(O5) = (3.D0 / 2.D0) * phi_e_tr_hat(O5)
          phi_e_tr_hat(O5) = SQRT(phi_e_tr_hat(O5))
  
  
          ptilde_hat(O5) = phi_p_tr_hat(O5)
          PhiEq_hat(O5) = phi_e_tr_hat(O5)
  
          GAMMA_hat(O5) = 0.0D0
        END DO

        CALL getC(sigma_c0, h_c1, h_c2, h_cexp, gma_n, sigma_c,HHc)
        sigma_t = m * sigma_c
        HHt = HHc
        CALL geta(alpha, m, sigma_c, a0, a1, a2)
        
        CALL getB(h_b0, h_b1, h_b2, gma_n, b_e, HHb, dHHbdgma)

        GG_til = GG_e + (k/2.D0)*HHb
        KK_til = KK_e + (k/3.D0)*HHb

        F_tr = a2 * PhiEq**alpha  - a1 * ptilde - a0

        ! DfDGamma = 0.D0
        

        u = 1.D0
        v = 1.D0

         A = (6.D0 * (PhiEq**2.D0)
     1        + (4.D0 / 3.D0) * (beta**2.D0) 
     2        * (ptilde ** 2.D0)) ** 0.5D0

        ! dDgammaDGamma = 0.D0
        ! Dgamma = 0.D0

        DO O5 = 1, 6
          F_tr_hat(O5) = a2 * PhiEq_hat(O5)**alpha  
     1                 - a1 * ptilde_hat(O5) - a0

        ! DfDGamma = 0.D0

        u_hat(O5) = 1.D0
        v_hat(O5) = 1.D0

         A_hat(O5) = (6.D0 * (PhiEq_hat(O5)**2.D0)
     1        + (4.D0 / 3.D0) * (beta**2.D0) 
     2        * (ptilde_hat(O5) ** 2.D0)) ** 0.5D0

        ! dDgammaDGamma = 0.D0
        ! Dgamma = 0.D0
        END DO

        
        ! Check for yielding
        WRITE(*,*) "F_tr : ", F_tr
        IF(F_tr .LE. 9.99D-6) THEN
          ! Update VE internal variables
          STATEV(10) = B_1
          ! Store A_1 in the statev vector
          DO O1 = 1, 3
            STATEV(10 + O1) = A_1(O1, O1)
          END DO
          STATEV(14) = A_1(1, 2)
          STATEV(15) = A_1(1, 3)
          STATEV(16) = A_1(2, 3)
          STATEV(17) = A_1(2, 1)
          STATEV(18) = A_1(3, 1)
          STATEV(19) = A_1(3, 2)

          ! Store E_ve in the statev vector
          STATEV(20) = E_ve_tr(1, 1)
          STATEV(21) = E_ve_tr(2, 2)
          STATEV(22) = E_ve_tr(3, 3)
          STATEV(23) = E_ve_tr(1, 2)
          STATEV(24) = E_ve_tr(1, 3)
          STATEV(25) = E_ve_tr(2, 3)
          STATEV(26) = E_ve_tr(2, 1)
          STATEV(27) = E_ve_tr(3, 1)
          STATEV(28) = E_ve_tr(3, 2)

          !Return cauchy stress for abaqus
          DO ii = 1, 3
            STRESS(ii) = (ONE / J) * tau_tr(ii, ii)
          END DO
          STRESS(4) = (ONE / J) * tau_tr(1, 2)
          STRESS(5) = (ONE / J) * tau_tr(1, 3)
          STRESS(6) = (ONE / J) * tau_tr(2, 3)

          !Turn perturbated stress into voit
          tau_tr_hat_v(:, :) = ZERO
          DO O5 = 1, 6
            DO ii = 1, 3
              tau_tr_hat_v(O5, ii) = tau_tr_hat(O5, ii, ii)
            END DO
            tau_tr_hat_v(O5, 4) = tau_tr_hat(O5, 1, 2)
            tau_tr_hat_v(O5, 5) = tau_tr_hat(O5, 1, 3)
            tau_tr_hat_v(O5, 6) = tau_tr_hat(O5, 2, 3)
          END DO

          !Tangent for Abaqus
          DO ii = 1, 6
            DO jj = 1, 6
              DDSDDE(ii, jj) = (ONE / (J * TOLL))
     1                     *  (tau_tr_hat_v(jj, ii) - J*STRESS(ii)) 
            END DO
          END DO

        ELSE

          ! VEVP  Here

          CALL nlinSolver(F_tr, eta, DTIME, GG_til, PhiEq, u,
     1      KK_til, beta, ptilde, v, A, k, GAMMA, HHt, sigma_c, HHc,
     2    sigma_t, alpha, m, a0, a1, a2, p_exp, gma_n, sigma_c0,
     3    h_c1, h_c2, h_cexp, h_b0, h_b1, h_b2, GG_e, KK_e, phi_e_tr,
     4    phi_p_tr)

          dev_phi(:,:) = dev_phi_tr(:,:) / u
          ptilde = phi_p_tr / v

          dev_Q(:, :) = dev_phi(:,:)
          dev_Q(:, :) = dev_Q(:, :) * 3.D0
          tr_Q = 2.D0 * beta * ptilde
          Q(:, :) = dev_Q(:,:)
          Q(1, 1) = Q(1, 1) + tr_Q / 3.D0 
          Q(2, 2) = Q(2, 2) + tr_Q / 3.D0
          Q(3, 3) = Q(3, 3) + tr_Q / 3.D0

          CALL approx_exp(GAMMA * Q(:,:), order, exp_GQ(:,:))

          CALL mat2dot(exp_GQ(:,:), F_vp_n(:,:), F_vp(:,:))

          STATEV(38) = gma_n

          CALL matInv(F_vp(:,:), F_vp_inv(:,:))
          CALL mat2dot(DFGRD1(:,:), F_vp_inv(:,:), F_ve(:, :))
          CALL mat2Tdot(F_ve(:,:), F_ve(:,:), C_ve(:,:))

          AA_upd(:,:) = C_ve(:,:) - I_mat(:,:)
          CALL approx_log(AA_upd(:,:), order, logAA_upd(:,:),
     1        dlogAA_upd(:,:,:,:), ddlogAA_upd(:,:,:,:,:,:))

          E_ve(:,:) = 0.5D0 * logAA_upd(:,:)

          CALL vePredictor_log(E_ve(:,:), E_ve_n(:, :), A_1_n(:,:),
     1   B_1_n, DTIME, GG_inf, GG_1, g_1, KK_inf, KK_1, k_1,
     2   kappa(:, :), A_1(:, :), B_1)

          STATEV(1) = F_vp(1, 1)
          STATEV(2) = F_vp(2, 2)
          STATEV(3) = F_vp(3, 3)
          STATEV(4) = F_vp(1, 2)
          STATEV(5) = F_vp(1, 3)
          STATEV(6) = F_vp(2, 3)
          STATEV(7) = F_vp(2, 1)
          STATEV(8) = F_vp(3, 1)
          STATEV(9) = F_vp(3, 2)

          STATEV(10) = B_1
          DO O1 = 1, 3
            STATEV(10 + O1) = A_1(O1, O1)
          END DO
          STATEV(14) = A_1(1, 2)
          STATEV(15) = A_1(1, 3)
          STATEV(16) = A_1(2, 3)
          STATEV(17) = A_1(2, 1)
          STATEV(18) = A_1(3, 1)
          STATEV(19) = A_1(3, 2)

          STATEV(20) = E_ve(1, 1)
          STATEV(21) = E_ve(2, 2)
          STATEV(22) = E_ve(3, 3)
          STATEV(23) = E_ve(1, 2)
          STATEV(24) = E_ve(1, 3)
          STATEV(25) = E_ve(2, 3)
          STATEV(26) = E_ve(2, 1)
          STATEV(27) = E_ve(3, 1)
          STATEV(28) = E_ve(3, 2)

         b(:,:) = b_tr(:,:) + k*HHb * GAMMA *Q(:,:)

         STATEV(29) = b(1, 1)
         STATEV(30) = b(2, 2)
         STATEV(31) = b(3, 3)
         STATEV(32) = b(1, 2)
         STATEV(33) = b(1, 3)
         STATEV(34) = b(2, 3)
         STATEV(35) = b(2, 1)
         STATEV(36) = b(3, 1)
         STATEV(37) = b(3, 2)

         ! Calculate S as per Ludovic's paper
         CALL mat24ddot(kappa(:, :), dlogAA_upd(:, :, :, :), S(:, :))
         ! Calculate kirchoff stress
         CALL mat2dot(F_ve(:, :), S(:, :), temp6(:, :))
         CALL mat2dotT(temp6(:, :), F_ve(:, :), tau(:, :))

          DO O5 = 1, 6

            CALL nlinSolver(F_tr_hat(O5), eta, DTIME, GG_til,
     1        PhiEq_hat(O5), u_hat(O5), KK_til, beta, ptilde_hat(O5),
     2        v_hat(O5), A_hat(O5), k, GAMMA_hat(O5), HHt, sigma_c,
     3        HHc, sigma_t, alpha, m, a0, a1, a2, p_exp, gma_n,
     4        sigma_c0, h_c1, h_c2, h_cexp, h_b0, h_b1, h_b2, GG_e,
     5        KK_e, phi_e_tr_hat(O5), phi_p_tr_hat(O5))

            
            dev_phi_hat(O5,:,:) = dev_phi_tr_hat(O5,:,:) / u_hat(O5)
            ptilde_hat(O5) = phi_p_tr_hat(O5) / v_hat(O5)

            dev_Q_hat(O5,:, :) = dev_phi_hat(O5, :,:)
            dev_Q_hat(O5,:, :) = dev_Q_hat(O5, :, :) * 3.D0
            tr_Q_hat(O5) = 2.D0 * beta * ptilde_hat(O5)
            Q_hat(O5, :, :) = dev_Q_hat(O5, :,:)
            Q_hat(O5,1, 1) = Q_hat(O5,1, 1) + tr_Q_hat(O5) / 3.D0 
            Q_hat(O5,2, 2) = Q_hat(O5,2, 2) + tr_Q_hat(O5) / 3.D0
            Q_hat(O5,3, 3) = Q_hat(O5,3, 3) + tr_Q_hat(O5) / 3.D0
          
            CALL approx_exp(GAMMA_hat(O5) * Q_hat(O5,:,:), order,
     1        exp_GQ_hat(O5,:,:))

            CALL mat2dot(exp_GQ_hat(O5,:,:), F_vp_n(:,:),
     1        F_vp_hat(O5,:,:))

            CALL matInv(F_vp_hat(O5,:,:), F_vp_inv_hat(O5,:,:))
            CALL mat2dot(F_hat(O5,:,:), F_vp_inv_hat(O5,:,:),
     1        F_ve_hat(O5,:, :))
            CALL mat2Tdot(F_ve_hat(O5, :,:), F_ve_hat(O5,:,:),
     1        C_ve_hat(O5,:,:))

            AA_upd_hat(O5,:,:) = C_ve_hat(O5,:,:) - I_mat(:,:)

            CALL approx_log(AA_upd_hat(O5, :,:), order, 
     1       logAA_upd_hat(O5,:,:), dlogAA_upd_hat(O5,:,:,:,:),
     2       ddlogAA_upd_hat(O5,:,:,:,:,:,:))

            E_ve_hat(O5,:,:) = 0.5D0 * logAA_upd_hat(O5,:,:)

            CALL vePredictor_log(E_ve_hat(O5,:,:), E_ve_n(:, :),
     1        A_1_n(:,:), B_1_n, DTIME, GG_inf, GG_1, g_1, KK_inf,
     2        KK_1, k_1, kappa_hat(O5, :, :), A_1(:, :), B_1)

            ! Calculate S as per Ludovic's paper
            CALL mat24ddot(kappa_hat(O5, :, :),
     1           dlogAA_upd_hat(O5, :, :, :, :), S_hat(O5, :, :))
         ! Calculate kirchoff stress
            CALL mat2dot(F_ve_hat(O5,:, :), S_hat(O5, :, :),
     1        temp6_hat(O5,:, :))
            CALL mat2dotT(temp6_hat(O5,:, :), F_ve_hat(O5,:, :),
     1        tau_hat(O5,:, :))

          END DO

          ! Start Returning stuff here
          !Return cauchy stress for abaqus
          DO ii = 1, 3
            STRESS(ii) = (ONE / J) * tau(ii, ii)
          END DO
          STRESS(4) = (ONE / J) * tau(1, 2)
          STRESS(5) = (ONE / J) * tau(1, 3)
          STRESS(6) = (ONE / J) * tau(2, 3)


          !Turn perturbated stress into voit
          tau_hat_v(:, :) = ZERO
          DO O5 = 1, 6
            DO ii = 1, 3
              tau_hat_v(O5, ii) = tau_hat(O5, ii, ii)
            END DO
            tau_hat_v(O5, 4) = tau_hat(O5, 1, 2)
            tau_hat_v(O5, 5) = tau_hat(O5, 1, 3)
            tau_hat_v(O5, 6) = tau_hat(O5, 2, 3)
          END DO

          !Tangent for Abaqus
          DO ii = 1, 6
            DO jj = 1, 6
              DDSDDE(ii, jj) = (ONE / (J * TOLL))
     1                     *  (tau_hat_v(jj, ii) - J*STRESS(ii)) 
            END DO
          END DO

        END IF
        

      RETURN
      END SUBROUTINE VEVP



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
      !      Helper function to compute tr-dev split of matrix
      !--------------------------------------------------------------

      SUBROUTINE tr_dev_split(matrix, dev, tr)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec), INTENT(OUT) :: dev(3, 3), tr
        REAL(prec) :: I_mat(3, 3)
        INTEGER :: ii, jj, kk, ll

        !Define 2nd order identity tensor
        data I_mat(1,:) /1.D0, 0.D0, 0.D0/
        data I_mat(2,:) /0.D0, 1.D0, 0.D0/
        data I_mat(3,:) /0.D0, 0.D0, 1.D0/
        
        tr = matrix(1, 1) + matrix(2, 2) + matrix(3, 3)

        dev(:, :) = matrix(:, :) - (1.D0 / 3.D0) * tr * I_mat(:, :)
        RETURN
      END SUBROUTINE tr_dev_split

      !--------------------------------------------------------------
      !      Function to compute perturb of matrix
      !--------------------------------------------------------------

      SUBROUTINE perturb_F(F, F_hat)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: F(3, 3)
        REAL(prec), INTENT(OUT) :: F_hat(6, 3, 3)
        REAL(prec) :: DD1(3, 3, 3, 3), DD2(3, 3, 3, 3), I_mat(3, 3)
        REAL(prec) :: delF(3, 3, 3, 3), F_hat_mat(3, 3, 3, 3)
        INTEGER :: ii, jj, ll, mm, nn

        !Decleration of constants
        REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
        REAL(prec), PARAMETER :: THREE=3.D0, FOUR= 4.D0, SIX=6.D0
        REAL(prec), PARAMETER :: NINE=9.D0, TOLL=0.001D0

        !Define 4th order identity tensor
        data I_mat(1,:) /ONE, ZERO, ZERO/
        data I_mat(2,:) /ZERO, ONE, ZERO/
        data I_mat(3,:) /ZERO, ZERO, ONE/

        ! Calculate Perturbation of F
        DD1(:, :, :, :) = ZERO
        DD2(:, :, :, :) = ZERO
        DO ii = 1, 3
          DO jj = 1, 3
            DO ll = 1, 3
              DO mm = 1, 3
                DO nn = 1, 3
          DD1(ii, jj, ll, mm) = DD1(ii, jj, ll, mm) 
     1              + (TOLL/TWO) * (I_mat(ll, ii) * I_mat(nn, jj) 
     2              * F(nn, mm))
          DD2(ii, jj, ll, mm) = DD2(ii, jj, ll, mm) 
     1              + (TOLL/TWO) * (I_mat(ll, jj) * I_mat(nn, ii) 
     2              * F(nn, mm))
                END  DO
              END DO
            END DO
          END DO
        END DO

        delF(:, :, :, :) = ZERO
        DO ii = 1, 3
          DO jj = 1, 3
            delF(ii, jj, :, :) = DD1(ii, jj, :, :) 
     1                             + DD2(ii, jj, :, :)
          END DO
        END DO

        F_hat_mat(:, :, :, :) = ZERO
        DO ii = 1, 3
          DO jj = 1, 3
          F_hat_mat(ii, jj, :, :) = F(:, :)+ delF(ii, jj, :, :)
          END DO
        END DO

        F_hat(1, : , :) = F_hat_mat(1, 1, :, :)
        F_hat(2, : , :) = F_hat_mat(2, 2, :, :)
        F_hat(3, : , :) = F_hat_mat(3, 3, :, :)
        F_hat(4, : , :) = F_hat_mat(1, 2, :, :)
        F_hat(5, : , :) = F_hat_mat(1, 3, :, :)
        F_hat(6, : , :) = F_hat_mat(2, 3, :, :)

        RETURN

      END SUBROUTINE perturb_F

      !--------------------------------------------------------------
      !      Function to compute log of matrix
      !--------------------------------------------------------------

      SUBROUTINE approx_log(AA, order, logAA, dlogAA, ddlogAA)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        INTEGER, INTENT(IN)      :: order
        REAL(prec), INTENT(IN)   :: AA(3, 3)
        REAL(prec), INTENT(OUT)  :: logAA(3, 3), dlogAA(3, 3, 3, 3)
        REAL(prec), INTENT(OUT)  :: ddlogAA(3, 3, 3, 3, 3, 3)

        ! List of internal variables
        INTEGER    :: ii, jj, nn, mm, ll, rr, ss, pp, qq
        REAL(prec) :: coeffs(order + 1), I_mat(3, 3)
        REAL(prec) :: FF(3, 3), dFF(3, 3, 3, 3)
        REAL(prec) :: ddFF(3, 3, 3, 3, 3, 3), II_mat(3, 3, 3, 3)
        REAL(prec) :: temp1(3, 3, 3, 3), temp2(3, 3, 3, 3)

        ! Define 2nd order identity tensor
        data I_mat(1,:) /1.D0, 0.D0, 0.D0/
        data I_mat(2,:) /0.D0, 1.D0, 0.D0/
        data I_mat(3,:) /0.D0, 0.D0, 1.D0/

        ! Define 4th order symetric identity tensor 
        II_mat(:,:,:,:) = 0.D0
        temp1(:,:,:,:) = 0.D0
        temp2(:,:,:,:) = 0.D0
        DO ii = 1, 3
          DO jj = 1, 3
            DO mm = 1, 3
              DO ll = 1, 3
                temp1(ii,jj,mm,ll) = temp1(ii,jj,mm,ll) 
     1                             + I_mat(ii,mm) * I_mat(jj,ll)
                temp2(ii,jj,mm,ll) = temp2(ii,jj,mm,ll) 
     1                             + I_mat(ii,ll) * I_mat(jj,mm)
              END DO
            END DO
          END DO
        END DO
        II_mat(:,:,:,:) = 0.5D0 * (temp1(:,:,:,:) + temp2(:,:,:,:))


        ! Calculate coeffs for the power law approxiamation of log
        DO ii = 1, order+1
          IF (ii .EQ. 1) THEN
           coeffs(ii) = 0.D0
          ELSE IF (modulo(ii, 2) .EQ. 0) THEN
           coeffs(ii) = 1.D0 / (real(ii, 8) - 1.D0)
          ELSE
           coeffs(ii) = -1.D0 / (real(ii, 8) - 1.D0)
          END IF
        END DO

        ! Implementation of polynomial func DG3D
        logAA(:, :) = 0.D0
        DO ii = 1, 3
          DO jj = 1, 3
            logAA(ii, jj) = I_mat(ii, jj) * coeffs(order + 1)
          END DO
        END DO

        dlogAA(:, :, :, :) = 0.D0
        ddlogAA(:, :, :, :, :, :) = 0.D0

        DO ii = 1, order
        FF(:, :) = 0.D0
        dFF(:, :, :, :) = 0.D0
        ddFF(:, :, :, :, :, :) = 0.D0

        nn = order + 1 - ii

        DO jj = 1, 3
          DO mm = 1, 3
            FF(jj, mm) = coeffs(nn)*I_mat(jj, mm)

            DO ll = 1, 3
              FF(jj, mm) = FF(jj, mm) + logAA(jj, ll) * AA(ll, mm)

              DO rr = 1, 3
                DO ss = 1, 3
                  dFF(jj, mm, rr, ss) = dFF(jj, mm, rr, ss) 
     1             + dlogAA(jj, ll, rr, ss) * AA(ll, mm) 
     2             + logAA(jj, ll) * II_mat(ll, mm, rr, ss)
                  
                  DO pp = 1, 3
                    DO qq = 1, 3
            ddFF(jj, mm, rr, ss, pp, qq) = ddFF(jj, mm, rr, ss, pp, qq)
     1               + ddlogAA(jj, ll, rr, ss, pp, qq) * AA(ll, mm)
     2               + dlogAA(jj, ll, rr, ss) * II_mat(ll, mm, pp, qq)
     3               + dlogAA(jj, ll, pp, qq) * II_mat(ll, mm, rr, ss)
                  END DO
                END DO
              
              END DO
            END DO
            
          END DO

          END DO
        END DO
        logAA(:, :) = FF(:, :)
        dlogAA(:,:,:,:) = dFF(:,:,:,:)
        ddlogAA(:,:,:,:,:,:) = ddFF(:, :, :, :, :, :)
      END DO

      END SUBROUTINE approx_log

      !--------------------------------------------------------------
      !      Function to compute exp of matrix
      !--------------------------------------------------------------

      SUBROUTINE approx_exp(Q, order, expQ)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        INTEGER, INTENT(IN)      :: order
        REAL(prec), INTENT(IN)   :: Q(3, 3)
        REAL(prec), INTENT(OUT)  :: expQ(3, 3)


        INTEGER    :: ii, jj, nn, mm, O5
        REAL(prec) :: Q_temp(order,  3, 3), I_mat(3, 3)
        REAL(prec) :: coeffs(order - 1)

        !Define 2nd order identity tensor
        data I_mat(1,:) /1.D0, 0.D0, 0.D0/
        data I_mat(2,:) /0.D0, 1.D0, 0.D0/
        data I_mat(3,:) /0.D0, 0.D0, 1.D0/


        Q_temp(:, :, :) = 0.D0
        Q_temp(1, :, :) = I_mat(:, :)
        Q_temp(2, :, :) = Q(:, :)

        DO O5 = 3, order
          DO ii = 1, 3
            DO jj = 1, 3
              DO mm = 1, 3
                Q_temp(O5, ii, jj) = Q_temp(O5, ii, jj) 
     1           + Q_temp(O5 - 1, ii, mm) * Q(mm, jj) 
              END DO
            END DO
          END DO          
        END DO

        coeffs(1) = 1.d0
        DO O5 = 2, order
          coeffs(O5) = coeffs(O5 - 1) * REAL(O5, 8) 
        END  DO

        expQ(:, :) = 0.D0
        expQ(:, :) = Q_temp(1, :, :)

        DO O5 = 2, order
          expQ(:, :) = expQ(:, :) 
     1     + Q_temp(O5, :, :) / coeffs(O5 - 1)
        END DO
      END SUBROUTINE approx_exp



      SUBROUTINE getC(sigma_c0, h_c1, h_c2, h_cexp, gma, sigma_c,HHc)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: sigma_c0, h_c1, h_c2, h_cexp, gma
        REAL(prec), INTENT(OUT) :: sigma_c, HHc

        REAL(prec) :: temp

        temp = EXP(-h_cexp * gma)

        sigma_c = sigma_c0 + h_c1 * gma + h_c2*(1.D0 - temp)
        HHc = h_c1  + h_c2 * h_cexp * temp
      END SUBROUTINE getC

      SUBROUTINE getB(h_b0, h_b1, h_b2, gma, b_e, HHb, dHHbdgma)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.D0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: h_b0, h_b1, h_b2, gma
        REAL(prec), INTENT(OUT) :: b_e, HHb, dHHbdgma

        b_e = h_b0 + h_b1 * gma
        HHb = h_b1 
        dHHbdgma = 0.D0
      END SUBROUTINE getB

      SUBROUTINE geta(alpha, m, sigma_c, a0, a1, a2)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.D0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: alpha, m, sigma_c
        REAL(prec), INTENT(OUT) :: a0, a1, a2

        REAL(prec) :: mral, mp1, sigcral

        mral = m**alpha
        mp1 = m + 1.D0
        sigcral = sigma_c**alpha

        a2 = 1.D0 / sigcral

        a1 = 3.D0 * ( (mral - 1.D0)/(mp1) ) * (1.D0 / sigma_c)

        a0 = (mral + m) / mp1

      END SUBROUTINE geta


      SUBROUTINE mat2Tdot(a, b, c)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: a(3, 3), b(3, 3)
        REAL(prec), INTENT(OUT) :: c(3, 3)

        INTEGER  :: O1, ii, jj

        c(:, :) = 0.D0

        DO O1 = 1, 3
          DO ii = 1, 3
            DO jj = 1, 3
              c(ii, jj) = c(ii, jj) 
     1                    + a(O1, ii) * b(O1, jj)
            END DO
          END DO
        END DO

      END SUBROUTINE mat2Tdot

      SUBROUTINE mat2dotT(a, b, c)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: a(3, 3), b(3, 3)
        REAL(prec), INTENT(OUT) :: c(3, 3)

        INTEGER  :: O1, ii, jj

        c(:, :) = 0.D0

        DO O1 = 1, 3
          DO ii = 1, 3
            DO jj = 1, 3
              c(ii, jj) = c(ii, jj) 
     1                    + a(ii, O1) * b(jj, O1)
            END DO
          END DO
        END DO

      END SUBROUTINE mat2dotT

      SUBROUTINE mat2dot(a, b, c)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: a(3, 3), b(3, 3)
        REAL(prec), INTENT(OUT) :: c(3, 3)

        INTEGER  :: O1, ii, jj

        c(:, :) = 0.D0

        DO ii = 1, 3
          DO jj = 1, 3
            DO O1 = 1, 3
              c(ii, jj) = c(ii, jj) 
     1                  + a(ii, O1) * b(O1, jj)
            END DO
          END DO
        END DO

      END SUBROUTINE mat2dot

      SUBROUTINE mat2ddot(a, b, c)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: a(3, 3), b(3, 3)
        REAL(prec), INTENT(OUT) :: c

        INTEGER  :: O1, ii, jj

        c = 0.D0

        DO ii = 1, 3
          DO jj = 1, 3
              c = c + a(ii, jj) * b(ii, jj)
          END DO
        END DO

      END SUBROUTINE mat2ddot


      SUBROUTINE mat24ddot(a, b, c)
        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: a(3, 3), b(3, 3, 3, 3)
        REAL(prec), INTENT(OUT) :: c(3, 3)

        INTEGER  :: ii, jj, mm, nn

        DO ii = 1, 3
        DO jj= 1, 3
          c(ii, jj) = 0.D0
          DO mm = 1, 3
            DO nn = 1, 3
              c(ii, jj) = c(ii, jj) 
     1          + a(mm, nn) * b(mm, nn, ii, jj)
            END DO
          END DO
        END DO
       END DO

      END SUBROUTINE mat24ddot

      SUBROUTINE vevpSplit(F, F_vp, F_ve, C_ve)

        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: F(3, 3), F_vp(3, 3)
        REAL(prec), INTENT(OUT) :: F_ve(3, 3), C_ve(3, 3)
        REAL(prec) :: C(3, 3), F_vp_inv(3, 3), temp1(3, 3)


        ! Compute C
        CALL mat2Tdot(F(:, :), F(:, :), C(:, :))

        ! Compute inverse of int variable F_vp
        CALL matInv(F_vp(:, :), F_vp_inv(:, :))

        ! Compute F_ve
        CALL mat2dot(F(:, :), F_vp_inv(:, :), F_ve(:, :))

        ! Compute C_ve
        CALL mat2Tdot(F_vp_inv(:, :), C(:, :), temp1(:, :))
        CALL mat2dot(temp1(:, :), F_vp_inv(:, :), C_ve(:, :))

      END SUBROUTINE vevpSplit

      SUBROUTINE vePredictor_log(E_ve, E_ve_n, Ai_n, Bi_n,
     1    dt, GGinf, GG1, g1, KKinf, KK1, k1, kappa, Ai, Bi)

         IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN)  :: E_ve(3, 3), E_ve_n(3, 3)
        REAL(prec), INTENT(IN)  :: Ai_n(3, 3), Bi_n, dt, GG1
        REAL(prec), INTENT(IN)  :: g1, KK1, k1, GGinf, KKinf
        REAL(prec), INTENT(OUT) :: kappa(3, 3), Ai(3, 3), Bi

        INTEGER  :: O1, ii, jj
        REAL(prec) :: dev_E_ve(3, 3), tr_E_ve, DE_ve(3, 3)
        REAL(prec) :: dev_DE_ve(3, 3), tr_DE_ve, dtg, expmdtg, ztag
        REAL(prec) :: dtk, expmdtk, ztak, dev_kappa(3, 3), GGe, KKe
        REAL(prec) :: p, I_mat(3, 3)

        !Define 2nd order identity tensor
        data I_mat(1,:) /1.D0, 0.D0, 0.D0/
        data I_mat(2,:) /0.D0, 1.D0, 0.D0/
        data I_mat(3,:) /0.D0, 0.D0, 1.D0/

        CALL tr_dev_split(E_ve(:,:), dev_E_ve(:,:), tr_E_ve)
        
        DE_ve(:, :) = E_ve(:, :) - E_ve_n(:, :)
        CALL tr_dev_split(DE_ve(:,:), dev_DE_ve(:,:), tr_DE_ve)

         dtg = dt / g1
         expmdtg = EXP(-dtg)
         ztag = EXP(-dtg/2.D0)
         GGe = GGinf + GG1 * EXP(-dt / (2.D0 * g1))

    !      Ai(:,:) = expmdtg*Ai_n(:, :) 
    !  1           + 2.D0 * GG1 * ztag * dev_DE_ve(:, :)
         Ai(:,:) = expmdtg*Ai_n(:, :) 
     1           + ztag * dev_DE_ve(:, :)
         
         dtk = dt / k1
         expmdtk = EXP(-dtk)
         ztak = exp(-dtk/2.D0)
         KKe = KKinf + KK1 * EXP(-dt / (2.D0 * k1))

        !  Bi = Bi_n * expmdtk + KK1 * ztak * tr_DE_ve
         Bi = Bi_n * expmdtk + ztak * tr_DE_ve

         dev_kappa(:, :) = dev_E_ve(:, :) * 2.D0 * GGinf 
     1                   + 2 * GG1 * Ai(:, :)
         p = tr_E_ve * KKinf + Bi * KK1
         kappa(:, :) = dev_kappa(:, :) + p * I_mat(:, :)

        END SUBROUTINE vePredictor_log


        SUBROUTINE nlinSolver(F_tr, eta, DTIME, GG_til, PhiEq, u,
     1      KK_til, beta, ptilde, v, A, k, GAMMA, HHt, sigma_c, HHc,
     2    sigma_t, alpha, m, a0, a1, a2, p_exp, gma_n, sigma_c0,
     3    h_c1, h_c2, h_cexp, h_b0, h_b1, h_b2, GG_e, KK_e, phi_e_tr,
     4    phi_p_tr)
          
          IMPLICIT NONE
          INTEGER, PARAMETER :: double=kind(1.d0)
          INTEGER, PARAMETER :: prec=double

          REAL(prec), INTENT(IN)    :: eta, DTIME, beta, k, alpha, m
          REAL(prec), INTENT(IN)    :: p_exp, sigma_c0, h_c1, h_c2
          REAL(prec), INTENT(IN)    :: h_cexp, h_b0, h_b1, h_b2
          REAL(prec), INTENT(IN)    :: GG_e, KK_e, phi_e_tr, phi_p_tr
          REAL(prec), INTENT(INOUT) :: F_tr, GG_til, PhiEq, u, KK_til
          REAL(prec), INTENT(INOUT) :: ptilde, v, A, GAMMA, HHt
          REAL(prec), INTENT(INOUT) :: sigma_c, HHc, sigma_t
          REAL(prec), INTENT(INOUT) :: a0, a1, a2, gma_n

          REAL(prec), PARAMETER :: TOLL_G=0.999D-6

          REAL(prec) :: iter_G, etaOverDt, dAdGamma, dDgammaDGamma,Dm
          REAL(prec) :: Da1Dm, H2, H1, H0, dfdDgamma, DfDGamma,Dgamma
          REAL(prec) :: b_e, HHb, dHHbdgma, viscoterm

          iter_G = 0

          DO WHILE (((ABS(F_tr) .GE. TOLL_G) .OR. (iter_G .LT. 1))
     1        .AND. (iter_G .LE. 250))
              etaOverDt = eta / DTIME
              dAdGamma = -(72.D0 * GG_til * PhiEq * PhiEq /u 
     1         + 16.D0 * KK_til * beta*beta*beta*ptilde*ptilde
     2         /(3.D0*v)) / (2.D0 * A)

              dDgammaDGamma = k*(A + GAMMA * dAdGamma)

              Dm = (HHt*sigma_c - HHc*sigma_t) / (sigma_c * sigma_c)

              Da1Dm  = (3.D0 / sigma_c) 
     1             * (alpha * (m**(alpha - 1.D0))/(m + 1.D0) 
     2             - (((m**alpha) - 1.D0)/(m + 1.D0))/(m + 1.D0)) 
            
              H2 = -alpha * (sigma_c**(-alpha - 1.D0)) * HHc

              H1 = Da1Dm * Dm 
     1           - 3.D0 * ((((m**alpha) - 1.D0) / (m + 1.D0))
     2           /(sigma_c * sigma_c))*HHc

              H0 =((alpha * (m**(alpha - 1.D0)) + 1.D0) / (m + 1.D0) 
     1            - (((m**alpha) + m)/(m + 1.D0))/(m + 1.D0))*Dm

              dfdDgamma = H2 * (PhiEq**alpha) -H1*ptilde - H0

              DfDGamma = (dfdDgamma * dDgammaDGamma) 
     1       - (alpha * a2 * 6.D0 * GG_til) * (PhiEq**alpha) / u
     2       + a1 * ptilde * 2.D0 * beta * KK_til / v

              IF ((GAMMA .GT. 0.D0) .AND. (etaOverDt .GT. 0.D0)) THEN
               DfDGamma = DfDGamma - (etaOverDt**p_exp) * p_exp 
     1                  * (GAMMA**(p_exp - 1.D0))
              END IF

              dGamma = -F_tr / DfDGamma

              IF ((GAMMA + dGamma) .LE. 0.D0) THEN
                GAMMA = GAMMA / 2.D0

              ELSE
                GAMMA = GAMMA + dGAMMA
              END IF
            
              u = 1.D0 + 6.D0 * GG_til * GAMMA
              v = 1.D0 + 2.D0 * beta * KK_til * GAMMA

             PhiEq = phi_e_tr / u
             ptilde = phi_p_tr / v

             A = (6.D0 * (PhiEq**2.D0)
     1         + (4.D0 / 3.D0) * (beta**2.D0) 
     2         * (ptilde ** 2.D0)) ** 0.5D0

             Dgamma = k * GAMMA * A

             gma_n = gma_n + Dgamma

             CALL getC(sigma_c0, h_c1, h_c2, h_cexp, gma_n, sigma_c,HHc)
             sigma_t = m * sigma_c
             HHt = HHc
             CALL geta(alpha, m, sigma_c, a0, a1, a2)
            
             CALL getB(h_b0, h_b1, h_b2, gma_n, b_e, HHb, dHHbdgma)

             GG_til = GG_e + (k/2.D0)*HHb
             KK_til = KK_e + (k/3.D0)*HHb

             F_tr = a2 * PhiEq**alpha  - a1 * ptilde - a0

             viscoterm = etaOverDt * GAMMA

             IF ((GAMMA .GT. 0.D0) .AND. (etaOverDt .GT. 0.D0)) THEN
               F_tr = F_tr - viscoterm ** p_exp
             END IF

             iter_G = iter_G + 1

            

          END DO
          WRITE(*,*) "Iter : ", iter_G
          WRITE(*,*) "GAMMA : ", GAMMA
          WRITE(*,*) "RESI : ", ABS(F_tr)
          WRITE(*,*) "-------------------------------------------"
        END SUBROUTINE nlinSolver