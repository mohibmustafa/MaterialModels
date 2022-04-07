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

C     Mohib Mustafa - IMDEA 15 March 2022


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
     &         NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1,
     &         PNEWDT)

      return
      end



C     !--------------------------------------------------------------
C     !                     UMAT SUBROUTINE
C     !--------------------------------------------------------------

      SUBROUTINE VEVP(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,
     1   PROPS,NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1,
     2   PNEWDT)


         IMPLICIT NONE

         INTEGER, PARAMETER :: double=kind(1.d0)
         INTEGER, PARAMETER :: prec=double
       
         INTEGER, INTENT(IN)      :: ntens,nprops,nstatv
         INTEGER, INTENT(IN)      :: kinc,kstep,noel
         REAL(prec), INTENT(IN)   :: stran(ntens), dstran(ntens)
         REAL(prec), INTENT(IN)   :: props(nprops), dtime
         REAL(prec), INTENT(IN)   :: DFGRD0(3,3), DFGRD1(3,3)
         REAL(prec), INTENT(INOUT) :: statev(nstatv)
         REAL(prec), INTENT(OUT)   :: stress(ntens)
         REAL(prec), INTENT(OUT)   :: ddsdde(ntens,ntens), PNEWDT    

         !List of local variables
         INTEGER    :: ii, jj, O6, order
         REAL(prec) :: I_mat(3,3)
         ! Internal Variables 
         REAL(prec) :: F_vp_n(3, 3), E_ve_n(3, 3), gma_n, b_n(3, 3)
         REAL(prec) :: AA_n(3, 3, 3), BB_n(3)
         ! Material Params
         REAL(prec) :: KK_inf, GG_inf, alpha, nu_p, eta, p_exp
         REAL(prec) :: sigmac0, hc1, hc2, hcexp
         REAL(prec) :: sigmat0, ht1, ht2, htexp
         REAL(prec) :: hb0, hb1, hb2
         REAL(prec) :: KK(3), k(3), GG(3), g(3), beta, k_plast
         ! VE  trial variables
         REAL(prec) :: F_ve_tr(3, 3), C_ve_tr(3, 3), D_ve_tr(3, 3)
         REAL(prec) :: E_ve_tr(3, 3), dlogC_ve_tr(3,3,3,3)
         REAL(prec) :: ddlogC_ve_tr(3,3,3,3,3,3), AA_tr(3,3,3)
         REAL(prec) :: BB_tr(3),GGe, KKe, kappa_tr(3, 3), S_tr(3, 3)
         REAL(prec) :: temp1(3,3), tau_tr(3, 3)
         ! variables for trial yield
         REAL(prec) :: phi_tr(3, 3), tr_phi_tr, dev_phi_tr(3, 3)
         REAL(prec) :: phi_p_tr, phi_e_tr, F_tr
         ! trial VE variables perturbation
         REAL(prec) :: F_hat(6, 3, 3), J, F_ve_tr_hat(6, 3, 3)
         REAL(prec) :: C_ve_tr_hat(6, 3, 3), D_ve_tr_hat(6,3,3)
         REAL(prec) :: E_ve_tr_hat(6,3,3), dlogC_ve_tr_hat(6,3,3,3,3)
         REAL(prec) :: ddlogC_ve_tr_hat(6,3,3,3,3,3,3)
         REAL(prec) :: AA_tr_hat(6,3,3,3)
         REAL(prec) :: BB_tr_hat(6, 3), GGe_hat(6), KKe_hat(6)
         REAL(prec) :: kappa_tr_hat(6, 3, 3), S_tr_hat(6, 3, 3)
         REAL(prec) :: temp1_hat(6, 3, 3), tau_tr_hat(6, 3, 3)
         REAL(prec) :: tau_tr_hat_v(6, 6)
         ! VP variables
         REAL(prec) :: ptilde, PhiEq, sigma_c, sigma_t, HHc, HHt, HHb
         REAL(prec) :: m, a0, a1, a2, GAMMA, u, v, dev_phi(3, 3)
         REAL(prec) :: dev_Q(3, 3), tr_Q, Q(3, 3), GQ(3,3)
         REAL(prec) :: exp_GQ(3, 3), F_vp(3,3), F_vp_inv(3,3)
         ! VE corrected  variables
         REAL(prec) :: F_ve(3, 3), C_ve(3, 3), D_ve(3, 3), E_ve(3, 3)
         REAL(prec) :: dlogC_ve(3,3,3,3), ddlogC_ve(3,3,3,3,3,3)
         REAL(prec) :: AA(3,3,3), BB(3), kappa(3, 3), b(3, 3)
         REAL(prec) :: S(3, 3), temp2(3, 3), tau(3, 3)
         ! VP variables perturbated
         REAL(prec) :: phi_tr_hat(6,3,3), dev_phi_tr_hat(6, 3, 3)
         REAL(prec) :: tr_phi_tr_hat(6), phi_p_tr_hat(6)
         REAL(prec) :: phi_e_tr_hat(6), ptilde_hat(6), PhiEq_hat(6)
         REAL(prec) :: gma_n_hat(6), sigma_c_hat(6), sigma_t_hat(6)
         REAL(prec) :: HHc_hat(6), HHt_hat(6), HHb_hat(6), m_hat(6)
         REAL(prec) :: a0_hat(6), a1_hat(6), a2_hat(6), F_tr_hat(6)
         REAL(prec) :: GAMMA_hat(6), u_hat(6), v_hat(6)
         REAL(prec) :: dev_phi_hat(6, 3, 3), dev_Q_hat(6, 3, 3)
         REAL(prec) :: tr_Q_hat(6), Q_hat(6,3,3)
         REAL(prec) :: GQ_hat(6,3,3), exp_GQ_hat(6,3,3)
         REAL(prec) :: F_vp_hat(6,3,3), F_vp_inv_hat(6,3,3)
         ! perturbated VE corrected variables
         REAL(prec) :: F_ve_hat(6,3,3), C_ve_hat(6,3,3)
         REAL(prec) :: D_ve_hat(6,3,3), E_ve_hat(6,3,3)
         REAL(prec) :: dlogC_ve_hat(6,3,3,3,3)
         REAL(prec) :: ddlogC_ve_hat(6,3,3,3,3,3,3), AA_hat(6,3,3,3)
         REAL(prec) :: BB_hat(6,3), kappa_hat(6,3,3), S_hat(6,3,3)
         REAL(prec) :: temp2_hat(6,3,3), tau_hat(6,3,3)
         REAL(prec) :: tau_hat_v(6,6)

        ! Tollerances
         REAL(prec), PARAMETER :: TOLL=0.001D0, TOLL_G=0.999D-11
         INTEGER, PARAMETER :: MAX_i=500

        !Define 2nd order identity tensor
         data I_mat(1,:) /1.D0, 0.D0, 0.D0/
         data I_mat(2,:) /0.D0, 1.D0, 0.D0/
         data I_mat(3,:) /0.D0, 0.D0, 1.D0/

        ! This is a workaround to set the diagnol values of F_vp = 1
        ! in FFTMAD. Can be removed for ABAQUS, IN abaqus use the inp
        ! file to set initail values for C_vp
         IF ((KINC .EQ. 1) .AND. (KSTEP .EQ. 1)) THEN
          STATEV(1) = 1.D0
          STATEV(2) = 1.D0
          STATEV(3) = 1.D0
         END IF

         ! State variables at previous time step
         CALL voit2mat(STATEV(1:9), F_vp_n(:,:))
         CALL voit2mat(STATEV(10:18), E_ve_n(:,:))
         gma_n = STATEV(19)
         CALL voit2mat(STATEV(20:28), b_n(:,:))
         CALL voit2mat(STATEV(29:37), AA_n(1,:,:))
         CALL voit2mat(STATEV(38:46), AA_n(2,:,:))
         CALL voit2mat(STATEV(47:55), AA_n(3,:,:))
         BB_n(1) = STATEV(56)
         BB_n(2) = STATEV(57)
         BB_n(3) = STATEV(58)

         ! Get Material Properties
         order         =PROPS(1)    
         KK_inf        =PROPS(2)     
         GG_inf        =PROPS(3)     
         alpha         =PROPS(4)    
         nu_p          =PROPS(5)   
         eta           =PROPS(6)  
         p_exp         =PROPS(7)    
         sigmac0       =PROPS(8) 

         hc1           =PROPS(9)  
         hc2           =PROPS(10)  
         hcexp         =PROPS(11)    
         sigmat0       =PROPS(12)      
         ht1           =PROPS(13)  
         ht2           =PROPS(14)  
         htexp         =PROPS(15)    
         hb0           =PROPS(16)

         hb1           =PROPS(17)  
         hb2           =PROPS(18)
         KK(1)         =PROPS(19)
         KK(2)         =PROPS(20)
         KK(3)         =PROPS(21)
         k(1)          =PROPS(22)
         k(2)          =PROPS(23)
         k(3)          =PROPS(24)

         GG(1)         =PROPS(25)
         GG(2)         =PROPS(26)
         GG(3)         =PROPS(27)
         g(1)          =PROPS(28)
         g(2)          =PROPS(29)
         g(3)          =PROPS(30)

         ! compute det(F)
         CALL determinant(DFGRD1(:,:), J)
        ! Calculate VE Right Cauchy Green Tensor (C_ve) at trial State
         CALL vevpSplit(DFGRD1(:, :), F_vp_n(:, :),
     1                  F_ve_tr(:, :), C_ve_tr(:, :))

         ! Approximate VE log Strain at trial state 
         D_ve_tr(:,:) = C_ve_tr(:, :) - I_mat(:,:)
         CALL approx_log(D_ve_tr(:,:), order, E_ve_tr(:, :),
     1               dlogC_ve_tr(:,:,:,:), ddlogC_ve_tr(:,:,:,:,:,:))
         E_ve_tr(:, :) = 0.5D0 * E_ve_tr(:, :)

         ! Calculate the corotated kirchoff stress at trial state along 
         ! with VE internal variables
         CALL vePredictor_log(E_ve_tr(:,:), E_ve_n(:,:), AA_n(:,:,:),
     1             BB_n(:), DTIME, PROPS, AA_tr(:,:,:), BB_tr(:), GGe,
     2             KKe, kappa_tr(:,:))

         ! Check yielding at trial state => gma = gma_n,  u=1, v=1, GAMMA = 0
           
         ! Get phi_e_tr, phi_p_tr
         phi_tr(:, :) = kappa_tr(:,:) - b_n(:,:)
             
         CALL tr_dev_split(phi_tr(:,:), dev_phi_tr(:,:), tr_phi_tr)
         phi_p_tr = tr_phi_tr / 3.D0
    
         CALL mat2ddot(dev_phi_tr(:, :), dev_phi_tr(:, :), phi_e_tr)
         phi_e_tr = (3.D0 / 2.D0) * phi_e_tr
         phi_e_tr = SQRT(phi_e_tr)
    
         ptilde = phi_p_tr
         PhiEq = phi_e_tr
    
         ! Update hardening variables on trial state
         CALL hardn(PROPS, gma_n, sigma_c, sigma_t, HHc, HHt, HHb)
         ! Update Drucker Pragger coefficients along with ecc factor m
         CALL DPcoeff(alpha, sigma_c, sigma_t, m, a0, a1, a2)
         
         F_tr = a2 * PhiEq**alpha  - a1 * ptilde - a0

          ! Numeric tangent
          CALL perturb_F(DFGRD1(:, :), TOLL, F_hat(:, :, :))
          DO O6 = 1, 6

            CALL vevpSplit(F_hat(O6, :, :), F_vp_n(:, :),
     1                  F_ve_tr_hat(O6, :, :), C_ve_tr_hat(O6, :, :))

            ! Approximate VE log Strain at trial state 
            D_ve_tr_hat(O6, :,:) = C_ve_tr_hat(O6, :, :) - I_mat(:,:)
            CALL approx_log(D_ve_tr_hat(O6, :,:), order,
     1        E_ve_tr_hat(O6, :, :), dlogC_ve_tr_hat(O6,:,:,:,:),
     2        ddlogC_ve_tr_hat(O6,:,:,:,:,:,:))
            E_ve_tr_hat(O6,:, :) = 0.5D0 * E_ve_tr_hat(O6, :, :)

            ! Calculate the corotated kirchoff stress at trial state
            ! along with VE internal variables
            CALL vePredictor_log(E_ve_tr_hat(O6, :,:), E_ve_n(:,:),
     1        AA_n(:,:,:), BB_n(:), DTIME, PROPS, AA_tr_hat(O6, :,:,:),
     2        BB_tr_hat(O6, :), GGe_hat(O6), KKe_hat(O6),
     3        kappa_tr_hat(O6, :,:)) 

            ! Calculate 2Pk as per paper
            CALL mat24ddot(kappa_tr_hat(O6, :, :),
     1        dlogC_ve_tr_hat(O6, :, :, :, :), S_tr_hat(O6, :, :))

            ! Calculate kirchoff stress
            CALL mat2dot(F_ve_tr_hat(O6, :, :), S_tr_hat(O6, :, :),
     1                    temp1_hat(O6, :, :))
            CALL mat2dotT(temp1_hat(O6, :, :), F_ve_tr_hat(O6, :, :),
     1                    tau_tr_hat(O6, :, :))
          
          END DO


         IF(F_tr .LE. TOLL_G) THEN
         
          ! Calculate 2Pk as per paper
          CALL mat24ddot(kappa_tr(:, :), dlogC_ve_tr(:, :, :, :),
     1                  S_tr(:, :))

          ! Calculate kirchoff stress
          CALL mat2dot(F_ve_tr(:, :), S_tr(:, :), temp1(:, :))
          CALL mat2dotT(temp1(:, :), F_ve_tr(:, :), tau_tr(:, :))

          ! Return cauchy stress for abaqus
          STRESS(1) = (1.D0 / J) * tau_tr(1, 1)
          STRESS(2) = (1.D0 / J) * tau_tr(2, 2)
          STRESS(3) = (1.D0 / J) * tau_tr(3, 3)
          STRESS(4) = (1.D0 / J) * tau_tr(1, 2)
          STRESS(5) = (1.D0 / J) * tau_tr(1, 3)
          STRESS(6) = (1.D0 / J) * tau_tr(2, 3)

          ! Update internal variables for trial state
          CALL mat2voit(E_ve_tr(:,:),STATEV(10:18))
          CALL mat2voit(AA_tr(1,:,:),STATEV(29:37))
          CALL mat2voit(AA_tr(2,:,:),STATEV(38:46))
          CALL mat2voit(AA_tr(3,:,:),STATEV(47:55))
          STATEV(56)=BB_tr(1)
          STATEV(57)=BB_tr(2)
          STATEV(58)=BB_tr(3)
          
          

          !Turn perturbated stress into voit
          tau_tr_hat_v(:, :) = 0.D0
          DO O6 = 1, 6
            DO ii = 1, 3
              tau_tr_hat_v(O6, ii) = tau_tr_hat(O6, ii, ii)
            END DO
            tau_tr_hat_v(O6, 4) = tau_tr_hat(O6, 1, 2)
            tau_tr_hat_v(O6, 5) = tau_tr_hat(O6, 1, 3)
            tau_tr_hat_v(O6, 6) = tau_tr_hat(O6, 2, 3)
          END DO

          !Tangent for Abaqus
          DO ii = 1, 6
            DO jj = 1, 6
              DDSDDE(ii, jj) = (1.D0 / (J * TOLL))
     1                     *  (tau_tr_hat_v(jj, ii) - J*STRESS(ii)) 
            END DO
          END DO

         

        ELSE
          ! Solve Newton Raphson to get GAMMA, gma and related variables
          CALL nlinSolver(PROPS, DTIME, TOLL_G, MAX_i,
     1      GGe, KKe,
     1      F_tr, gma_n, ptilde, PhiEq, GAMMA, u, v, beta, k_plast)

          ! Update dev_phi and ptilde (already returned by the subrout)
          dev_phi(:,:) = dev_phi_tr(:,:) / u

          ! Update flow normal Q
          dev_Q(:, :) = 3.D0 * dev_phi(:,:)
          tr_Q = 2.D0 * beta * ptilde / 3.D0
          Q(:,:) = dev_Q(:,:) + tr_Q * I_mat(:,:)
          
          ! Make corrections to F_vp, F_ve and C_ve
          GQ(:,:) = GAMMA * Q(:,:)
          CALL approx_exp(GQ(:,:), order, exp_GQ(:,:))
          CALL mat2dot(exp_GQ(:,:), F_vp_n(:,:), F_vp(:,:))
          CALL matInv(F_vp(:,:), F_vp_inv(:,:))
          CALL mat2dot(DFGRD1(:,:), F_vp_inv(:,:), F_ve(:, :))
          CALL mat2Tdot(F_ve(:,:), F_ve(:,:), C_ve(:,:))

          ! Approximate VE log Strain at the corrected state
          D_ve(:,:) = C_ve(:, :) - I_mat(:,:)
          CALL approx_log(D_ve(:,:), order, E_ve(:, :),
     1                     dlogC_ve(:,:,:,:), ddlogC_ve(:,:,:,:,:,:))
          E_ve(:, :) = 0.5D0 * E_ve(:, :)

          ! Calculate the corotated kirchoff stress at corrected state along 
          ! with VE internal variables
          CALL vePredictor_log(E_ve(:,:), E_ve_n(:,:), AA_n(:,:,:),
     1             BB_n(:), DTIME, PROPS, AA(:,:,:), BB(:), GGe, KKe,
     2             kappa(:,:))

          ! Calculate 2PK stress at  corrected state
          CALL mat24ddot(kappa(:, :), dlogC_ve(:, :, :, :), S(:, :))
          ! Calculate kirchoff stress at corrected state
          CALL mat2dot(F_ve(:, :), S(:, :), temp2(:, :))
          CALL mat2dotT(temp2(:, :), F_ve(:, :), tau(:, :))

          ! Return cauchy stress for abaqus
          STRESS(1) = (1.D0 / J) * tau(1, 1)
          STRESS(2) = (1.D0 / J) * tau(2, 2)
          STRESS(3) = (1.D0 / J) * tau(3, 3)
          STRESS(4) = (1.D0 / J) * tau(1, 2)
          STRESS(5) = (1.D0 / J) * tau(1, 3)
          STRESS(6) = (1.D0 / J) * tau(2, 3)

          ! Update kin hardening variable for the corrected state
          CALL hardn(PROPS, gma_n, sigma_c, sigma_t, HHc, HHt, HHb)
          b(:,:) = b_n(:,:) + k_plast*HHb * GQ(:,:)

          ! Return updated internal variables
          CALL mat2voit(F_vp(:,:),STATEV(1:9))
          CALL mat2voit(E_ve(:,:),STATEV(10:18))
          STATEV(19) = gma_n
          CALL mat2voit(b(:,:),STATEV(20:28))
          CALL mat2voit(AA(1,:,:),STATEV(29:37))
          CALL mat2voit(AA(2,:,:),STATEV(38:46))
          CALL mat2voit(AA(3,:,:),STATEV(47:55))
          STATEV(56)=BB(1)
          STATEV(57)=BB(2)
          STATEV(58)=BB(3)

          ! Numeric Tangent
          DO O6  =  1, 6
           
          gma_n_hat(O6) = gma_n
          
          ! Check yielding at Perturbated state 
           
          ! Get phi_e_tr, phi_p_tr
          phi_tr_hat(O6, :, :) = kappa_tr_hat(O6, :,:) - b_n(:,:)
             
          CALL tr_dev_split(phi_tr_hat(O6, :,:),
     1         dev_phi_tr_hat(O6, :,:), tr_phi_tr_hat(O6))
          phi_p_tr_hat(O6) = tr_phi_tr_hat(O6) / 3.D0
    
          CALL mat2ddot(dev_phi_tr_hat(O6, :, :),
     1             dev_phi_tr_hat(O6, :, :), phi_e_tr_hat(O6))
          phi_e_tr_hat(O6) = (3.D0 / 2.D0) * phi_e_tr_hat(O6)
          phi_e_tr_hat(O6) = SQRT(phi_e_tr_hat(O6))
    
          ptilde_hat(O6) = phi_p_tr_hat(O6)
          PhiEq_hat(O6) = phi_e_tr_hat(O6)
    
          ! Update hardening variables on puturb state
          CALL hardn(PROPS, gma_n_hat(O6), sigma_c_hat(O6),
     1      sigma_t_hat(O6), HHc_hat(O6), HHt_hat(O6), HHb_hat(O6))
          ! Update Drucker Pragger coefficients along with ecc factor m
          CALL DPcoeff(alpha, sigma_c_hat(O6), sigma_t_hat(O6),
     1        m_hat(O6), a0_hat(O6), a1_hat(O6), a2_hat(O6))
          ! Calculate yield func at purturbed state
          F_tr_hat(O6) = a2_hat(O6) * PhiEq_hat(O6)**alpha 
     1        - a1_hat(O6) * ptilde_hat(O6) - a0_hat(O6)

          ! Solve Newton raphson to get gma and GAMMA at purturbed state
          CALL nlinSolver(PROPS, DTIME, TOLL_G, MAX_i,
     1      GGe, KKe,
     1      F_tr_hat(O6), gma_n_hat(O6), ptilde_hat(O6),
     2      PhiEq_hat(O6), GAMMA_hat(O6), u_hat(O6), v_hat(O6),
     3      beta, k_plast)

          ! Update dev_phi_hat and ptilde_hat (already returned by the subrout)
          dev_phi_hat(O6, :,:) = dev_phi_tr_hat(O6,:,:) / u_hat(O6)

          ! Update flow normal Q_hat
          dev_Q_hat(O6, :, :) = 3.D0 * dev_phi_hat(O6, :,:)
          tr_Q_hat(O6) = 2.D0 * beta * ptilde_hat(O6) / 3.D0
          Q_hat(O6,:,:) = dev_Q_hat(O6,:,:) 
     1                  + tr_Q_hat(O6) * I_mat(:,:)
          
          ! Make corrections to F_vp_hat, F_ve_hat and C_ve_hat
          GQ_hat(O6,:,:) = GAMMA_hat(O6) * Q_hat(O6,:,:)
          CALL approx_exp(GQ_hat(O6,:,:), order, exp_GQ_hat(O6,:,:))
          CALL mat2dot(exp_GQ_hat(O6, :,:), F_vp_n(:,:),
     1                  F_vp_hat(O6,:,:))
          CALL matInv(F_vp_hat(O6,:,:), F_vp_inv_hat(O6,:,:))
          CALL mat2dot(F_hat(O6, :,:), F_vp_inv_hat(O6,:,:),
     1      F_ve_hat(O6, :, :))
          CALL mat2Tdot(F_ve_hat(O6,:,:), F_ve_hat(O6,:,:),
     1                  C_ve_hat(O6,:,:))


          ! Approximate VE log Strain at the purturbated state
          D_ve_hat(O6,:,:) = C_ve_hat(O6,:, :) - I_mat(:,:)
          CALL approx_log(D_ve_hat(O6,:,:), order, E_ve_hat(O6,:, :),
     1                     dlogC_ve_hat(O6,:,:,:,:),
     2                ddlogC_ve_hat(O6,:,:,:,:,:,:))
          E_ve_hat(O6,:, :) = 0.5D0 * E_ve_hat(O6,:, :)

          ! Calculate the corotated kirchoff stress at purturbed state along 
          ! with VE internal variables
          CALL vePredictor_log(E_ve_hat(O6,:,:), E_ve_n(:,:),
     1      AA_n(:,:,:), BB_n(:), DTIME, PROPS, AA_hat(O6,:,:,:),
     2      BB_hat(O6,:), GGe, KKe, kappa_hat(O6,:,:))

          ! Calculate 2PK stress at  purturbed state
          CALL mat24ddot(kappa_hat(O6,:, :),
     1      dlogC_ve_hat(O6,:, :, :, :), S_hat(O6,:, :))

          ! Calculate kirchoff stress at purturbed state
          CALL mat2dot(F_ve_hat(O6,:, :), S_hat(O6,:, :),
     1      temp2_hat(O6,:, :))
          CALL mat2dotT(temp2_hat(O6,:, :), F_ve_hat(O6,:, :),
     1      tau_hat(O6,:, :))

          END DO

           !Turn perturbated stress into voit
          tau_hat_v(:, :) = 0.D0
          DO O6 = 1, 6
            DO ii = 1, 3
              tau_hat_v(O6, ii) = tau_hat(O6, ii, ii)
            END DO
            tau_hat_v(O6, 4) = tau_hat(O6, 1, 2)
            tau_hat_v(O6, 5) = tau_hat(O6, 1, 3)
            tau_hat_v(O6, 6) = tau_hat(O6, 2, 3)
          END DO

          !Tangent for Abaqus
          DO ii = 1, 6
            DO jj = 1, 6
              DDSDDE(ii, jj) = (1.D0 / (J * TOLL))
     1                     *  (tau_hat_v(jj, ii) - J*STRESS(ii)) 
            END DO
          END DO
          
        END IF

          RETURN
       END SUBROUTINE VEVP

      !--------------------------------------------------------------
      !     Helper function to convert a voit array to matrix 
      !--------------------------------------------------------------

      SUBROUTINE voit2mat(voit, mat)
        IMPLICIT NONE
                
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: voit(9)
        REAL(prec), INTENT(OUT):: mat(3, 3)

        mat(1, 1) = voit(1)
        mat(2, 2) = voit(2)
        mat(3, 3) = voit(3)
        mat(1, 2) = voit(4)
        mat(1, 3) = voit(5)
        mat(2, 3) = voit(6)
        mat(2, 1) = voit(7)
        mat(3, 1) = voit(8)
        mat(3, 2) = voit(9)
        RETURN
      END SUBROUTINE

      !--------------------------------------------------------------
      !     Helper function to convert a  matrix to voit array
      !--------------------------------------------------------------

      SUBROUTINE mat2voit(mat, voit)
        IMPLICIT NONE
                
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) ::  mat(3, 3)
        REAL(prec), INTENT(OUT):: voit(9)

        voit(1) = mat(1, 1)  
        voit(2) = mat(2, 2)  
        voit(3) = mat(3, 3)  
        voit(4) = mat(1, 2)  
        voit(5) = mat(1, 3)  
        voit(6) = mat(2, 3)  
        voit(7) = mat(2, 1)  
        voit(8) = mat(3, 1)  
        voit(9) = mat(3, 2)  
        RETURN
      END SUBROUTINE

      !--------------------------------------------------------------
      !     Helper function to compute Determinant of 3x3 matrix 
      !--------------------------------------------------------------
      SUBROUTINE determinant(matrix, det)
        IMPLICIT NONE
                
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec), INTENT(OUT):: det

        det = matrix(1,1) * (matrix(2,2) * matrix(3,3) 
     1   - matrix(3,2) * matrix(2, 3)) - matrix(1,2) * (matrix(2,1)
     2   * matrix(3,3) - matrix(2,3) * matrix(3,1)) + matrix(1,3) 
     3   * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))
        RETURN
      END SUBROUTINE determinant


      !--------------------------------------------------------------
      !      Helper function to compute inverse of 3x3 matrix
      !--------------------------------------------------------------
      SUBROUTINE matInv(matrix, inv)
        IMPLICIT NONE
        
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec), INTENT(OUT) :: inv(3, 3)
        REAL(prec) :: J

        inv(1,1) = matrix(2,2)*matrix(3,3) - matrix(3,2)*matrix(2,3)
        inv(1,2) = matrix(3,2)*matrix(1,3) - matrix(1,2)*matrix(3,3)
        inv(1,3) = matrix(1,2)*matrix(2,3) - matrix(2,2)*matrix(1,3)
        inv(2,1) = matrix(3,1)*matrix(2,3) - matrix(2,1)*matrix(3,3)
        inv(2,2) = matrix(1,1)*matrix(3,3) - matrix(3,1)*matrix(1,3)
        inv(2,3) = matrix(2,1)*matrix(1,3) - matrix(1,1)*matrix(2,3)
        inv(3,1) = matrix(2,1)*matrix(3,2) - matrix(3,1)*matrix(2,2)
        inv(3,2) = matrix(3,1)*matrix(1,2) - matrix(1,1)*matrix(3,2)
        inv(3,3) = matrix(1,1)*matrix(2,2) - matrix(2,1)*matrix(1,2)
        
        CALL determinant(matrix, J)

        inv(:, :) = inv(:, :)/ J

        RETURN
      END SUBROUTINE matInv

      
      !------------------------------------------------------------
      !     Helper function to for c = a.b  
      !------------------------------------------------------------

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

      !------------------------------------------------------------
      !     Helper function to for c = a^T.b  
      !------------------------------------------------------------
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

      !------------------------------------------------------------
      !     Helper function
      !------------------------------------------------------------
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

      !------------------------------------------------------------
      !     Helper function
      !------------------------------------------------------------
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

       !--------------------------------------------------------------
      !     Helper function to for 2nd order 4th order contraction
      !--------------------------------------------------------------
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

      !--------------------------------------------------------------
      !     Helper function to get additive dev trace split  
      !--------------------------------------------------------------
      SUBROUTINE tr_dev_split(matrix, dev, tr)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec), INTENT(OUT) :: dev(3, 3), tr
        REAL(prec) :: I_mat(3, 3)

        !Define 2nd order identity tensor
        data I_mat(1,:) /1.D0, 0.D0, 0.D0/
        data I_mat(2,:) /0.D0, 1.D0, 0.D0/
        data I_mat(3,:) /0.D0, 0.D0, 1.D0/
        
        tr = matrix(1, 1) + matrix(2, 2) + matrix(3, 3)

        dev(:, :) = matrix(:, :) - (1.D0 / 3.D0) * tr * I_mat(:, :)
        RETURN
      END SUBROUTINE tr_dev_split

      !--------------------------------------------------------------
      !     Helper function to get VEVP Multiplicative Split
      !--------------------------------------------------------------
      SUBROUTINE vevpSplit(F, F_vp, F_ve, C_ve)

        IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN) :: F(3, 3), F_vp(3, 3)
        REAL(prec), INTENT(OUT) :: F_ve(3, 3), C_ve(3, 3)
        REAL(prec) :: F_vp_inv(3, 3)


        ! Compute inverse of int variable F_vp
        CALL matInv(F_vp(:, :), F_vp_inv(:, :))

        ! Compute F_ve
        CALL mat2dot(F(:, :), F_vp_inv(:, :), F_ve(:, :))

        ! Compute C_ve
        CALL mat2Tdot(F_ve(:, :), F_ve(:, :), C_ve(:, :))
      END SUBROUTINE vevpSplit

       !--------------------------------------------------------------
      !      Function to compute perturb of matrix
      !--------------------------------------------------------------

      SUBROUTINE perturb_F(F, TOLL, F_hat)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: F(3, 3), TOLL
        REAL(prec), INTENT(OUT) :: F_hat(6, 3, 3)
        REAL(prec) :: DD1(3, 3, 3, 3), DD2(3, 3, 3, 3), I_mat(3, 3)
        REAL(prec) :: delF(3, 3, 3, 3), F_hat_mat(3, 3, 3, 3)
        INTEGER :: ii, jj, ll, mm, nn

        !Decleration of constants
        REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
        REAL(prec), PARAMETER :: THREE=3.D0, FOUR= 4.D0, SIX=6.D0
        REAL(prec), PARAMETER :: NINE=9.D0

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
      !      Function to compute apprrox log strains and derivatives
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


      !--------------------------------------------------------------
      !      VE Predictor routine . 3 Branches
      !--------------------------------------------------------------

      SUBROUTINE vePredictor_log(E_ve, E_ve_n, AA_n, BB_n,
     1    dt, PROPS, AA, BB, GGe, KKe, kappa)

         IMPLICIT NONE
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double

        REAL(prec), INTENT(IN)  :: E_ve(3, 3), E_ve_n(3, 3)
        REAL(prec), INTENT(IN)  :: AA_n(3, 3, 3), BB_n(3), dt
        REAL(prec), INTENT(IN)  :: PROPS(30)
        REAL(prec), INTENT(OUT) :: AA(3, 3, 3), BB(3), GGe, KKe
        REAL(prec), INTENT(OUT) :: kappa(3, 3)

        INTEGER  :: ii
        REAL(prec) :: dev_E_ve(3, 3), tr_E_ve, DE_ve(3, 3)
        REAL(prec) :: dev_DE_ve(3, 3), tr_DE_ve, dtg, expmdtg, ztag
        REAL(prec) :: dtk, expmdtk, ztak, dev_kappa(3, 3)
        REAL(prec) :: p, I_mat(3, 3)

        ! Material Params
        REAL(prec) :: KK_inf, GG_inf
        REAL(prec) :: KK(3), k(3), GG(3), g(3)

        !Define 2nd order identity tensor
        data I_mat(1,:) /1.D0, 0.D0, 0.D0/
        data I_mat(2,:) /0.D0, 1.D0, 0.D0/
        data I_mat(3,:) /0.D0, 0.D0, 1.D0/

        ! Get Material Properties   
        KK_inf        =PROPS(2)     
        GG_inf        =PROPS(3)     
        KK(1)         =PROPS(19)
        KK(2)         =PROPS(20)
        KK(3)         =PROPS(21)
        k(1)          =PROPS(22)
        k(2)          =PROPS(23)
        k(3)          =PROPS(24)
        GG(1)         =PROPS(25)
        GG(2)         =PROPS(26)
        GG(3)         =PROPS(27)
        g(1)          =PROPS(28)
        g(2)          =PROPS(29)
        g(3)          =PROPS(30)

        CALL tr_dev_split(E_ve(:,:), dev_E_ve(:,:), tr_E_ve)
        
        DE_ve(:, :) = E_ve(:, :) - E_ve_n(:, :)
        CALL tr_dev_split(DE_ve(:,:), dev_DE_ve(:,:), tr_DE_ve)

         GGe = GG_inf
         AA(:,:,:)=0.D0

         KKe = KK_inf
         BB(:)=0.D0

         DO ii = 1,3
            dtg = dt / g(ii)
            expmdtg = EXP(-dtg)
            ztag = EXP(-dtg/2.D0)

            GGe = GGe + GG(ii) * ztag
            AA(ii,:,:) = expmdtg*AA_n(ii,:, :) 
     1                 + ztag * dev_DE_ve(:, :)
         
         
            dtk = dt / k(ii)
            expmdtk = EXP(-dtk)
            ztak = EXP(-dtk/2.D0)

            KKe = KKe + KK(ii) * ztak
            BB(ii) = BB_n(ii) * expmdtk + ztak  * tr_DE_ve
         END DO


         dev_kappa(:, :) = dev_E_ve(:, :) * 2.D0 * GG_inf 
         p = tr_E_ve * KK_inf

         DO ii = 1, 3
            dev_kappa(:, :) = dev_kappa(:, :) 
     1       + 2.D0 * GG(ii) * AA(ii, :, :)

            p = p + BB(ii) * KK(ii)
         END DO

         kappa(:, :) = dev_kappa(:, :) + p * I_mat(:, :)

            RETURN
        END SUBROUTINE vePredictor_log

         !--------------------------------------------------------------
         !      Iso Hardening Function (Linear exponential)
         !--------------------------------------------------------------

        SUBROUTINE hardn(PROPS, gma, sigma_c, sigma_t, HHc, HHt, HHb)
          IMPLICIT NONE
          INTEGER, PARAMETER :: double=kind(1.d0)
          INTEGER, PARAMETER :: prec=double
  
          REAL(prec), INTENT(IN) :: PROPS(30), gma
          REAL(prec), INTENT(OUT) :: sigma_c, sigma_t, HHc, HHt,HHb
  
          REAL(prec) :: sigmac0, hc1, hc2, hcexp, tempc
          REAL(prec) :: sigmat0, ht1, ht2, htexp, tempt
          REAL(prec) :: hb0, hb1, hb2

          sigmac0       =PROPS(8) 
          hc1           =PROPS(9)  
          hc2           =PROPS(10)  
          hcexp         =PROPS(11)    
          sigmat0       =PROPS(12)      
          ht1           =PROPS(13)  
          ht2           =PROPS(14)  
          htexp         =PROPS(15)
          hb0           =PROPS(16)
          hb1           =PROPS(17)  
          hb2           =PROPS(18)

          tempc = EXP(-hcexp * gma)
          sigma_c = sigmac0 + hc1 * gma + hc2*(1.D0 - tempc)
          HHc = hc1  + hc2 * hcexp * tempc

          tempt = EXP(-htexp * gma)
          sigma_t = sigmat0 + ht1 * gma + ht2*(1.D0 - tempt)
          HHt = ht1  + ht2 * htexp * tempt

          HHb = hb1
        END SUBROUTINE hardn


        !--------------------------------------------------------------
         !      Drucker - Prager coefficient update subroutine
         !--------------------------------------------------------------
        SUBROUTINE DPcoeff(alpha, sigma_c, sigma_t, m, a0, a1, a2)
          IMPLICIT NONE
          INTEGER, PARAMETER :: double=kind(1.D0)
          INTEGER, PARAMETER :: prec=double
  
          REAL(prec), INTENT(IN) :: alpha, sigma_c, sigma_t
          REAL(prec), INTENT(OUT) ::  m, a0, a1, a2
  
          REAL(prec) :: mral, mp1, sigcral
  
          m = sigma_t / sigma_c
  
          mral = m**alpha
          mp1 = m + 1.D0
          sigcral = sigma_c**alpha
  
          a2 = 1.D0 / sigcral
  
          a1 = 3.D0 * ( (mral - 1.D0)/(mp1) ) * (1.D0 / sigma_c)
  
          a0 = (mral + m) / mp1
  
        END SUBROUTINE DPcoeff

        SUBROUTINE nlinSolver(PROPS, DTIME, TOLL, MAX_i,
     1      GGe, KKe, F, gma,
     1      ptilde, PhiEq,  GAMMA, u, v, beta, k_plast)

          IMPLICIT NONE
          INTEGER, PARAMETER :: double=kind(1.d0)
          INTEGER, PARAMETER :: prec=double

          REAL(prec), INTENT(IN)    :: PROPS(30), DTIME, TOLL
          INTEGER, INTENT(IN)       ::  MAX_i
          REAL(prec), INTENT(IN)    :: GGe, KKe
          REAL(prec), INTENT(OUT)   :: GAMMA, beta, k_plast, u, v
          REAL(prec), INTENT(INOUT) :: F, gma, ptilde, PhiEq

          INTEGER    :: iter_G
          REAL(prec) :: gma_i, A, nu_p, GG_til, KK_til, sigma_c
          REAL(prec) :: phi_e_tr, phi_p_tr, alpha, eta, p_exp
          REAL(prec) :: sigma_t, HHc, HHt, HHb
          REAL(prec) :: dAdGamma, dDgmadGamma, Dm, Da1Dm
          REAL(prec) :: m, a0, a1, a2, H2, H1, H0, etaOverDt
          REAL(prec) :: dfdDgma, DfDGamma, dGamma, Dgma, viscoterm

          ! Material Props
          alpha  = PROPS(4)
          nu_p   = PROPS(5)
          eta    = PROPS(6)
          p_exp  = PROPS(7) 

          ! Calc derived material props
          beta = (9.D0 / 2.D0) 
     1          * ((1.D0 - 2.D0 * nu_p) / (nu_p + 1.D0))

          k_plast = 1.D0 / ((1.D0 + 2.D0 * nu_p ** 2.D0)**(0.5D0))
          
          

          phi_e_tr = PhiEq
          phi_p_tr = ptilde
          ! Initialize for first step
          iter_G = 0
          
          GAMMA = 0.0
          u = 1.D0
          v = 1.D0

          A = (6.D0 * (PhiEq**2.D0)
     1        + (4.D0 / 3.D0) * (beta**2.D0) 
     2        * (ptilde ** 2.D0)) ** 0.5D0

          gma_i = gma
          
          ! Update hardening variables on trial state
          CALL hardn(PROPS, gma_i, sigma_c, sigma_t, HHc, HHt, HHb)
          ! Update Drucker Pragger coefficients along with ecc factor m
          CALL DPcoeff(alpha, sigma_c, sigma_t, m, a0, a1, a2)

          GG_til = GGe + (k_plast/2.D0)*HHb
          KK_til = KKe + (k_plast/3.D0)*HHb

          DO WHILE (((ABS(F) .GE. TOLL) .OR. (iter_G .LT. 1))
     1        .AND. (iter_G .LE. MAX_i))

            etaOverDt = eta / DTIME

            dAdGamma = -(72.D0 * GG_til * PhiEq * PhiEq /u 
     1         + 16.D0 * KK_til * beta*beta*beta*ptilde*ptilde
     2         /(3.D0*v)) / (2.D0 * A)

            dDgmadGamma = k_plast*(A + GAMMA * dAdGamma)

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

            
            dfdDgma = H2 * (PhiEq**alpha) -H1*ptilde - H0

            DfDGamma = (dfdDgma * dDgmadGamma) 
     1       - (alpha * a2 * 6.D0 * GG_til) * (PhiEq**alpha) / u
     2       + a1 * ptilde * 2.D0 * beta * KK_til / v

            IF ((GAMMA .GT. 0.D0) .AND. (etaOverDt .GT. 0.D0)) THEN
               DfDGamma = DfDGamma - (etaOverDt**p_exp) * p_exp 
     1                  * (GAMMA**(p_exp - 1.D0))
            END IF

            dGamma = -F / DfDGamma

            IF ((GAMMA + dGamma) .LE. 0.D0) THEN
                GAMMA = GAMMA / 2.D0

            ELSE
                GAMMA = GAMMA + dGamma
            END IF

            u = 1.D0 + 6.D0 * GG_til * GAMMA
            v = 1.D0 + 2.D0 * beta * KK_til * GAMMA

            PhiEq = phi_e_tr / u
            ptilde = phi_p_tr / v

            A = (6.D0 * (PhiEq**2.D0)
     1         + (4.D0 / 3.D0) * (beta**2.D0) 
     2         * (ptilde ** 2.D0)) ** 0.5D0

            Dgma = k_plast * GAMMA * A
            gma_i = gma + Dgma

            ! Update hardening variables on gma_i
            CALL hardn(PROPS, gma_i, sigma_c, sigma_t, HHc, HHt, HHb)
            ! Update Drucker Pragger coefficients along with ecc factor m
            CALL DPcoeff(alpha, sigma_c, sigma_t, m, a0, a1, a2)

            GG_til = GGe + (k_plast/2.D0)*HHb
            KK_til = KKe + (k_plast/3.D0)*HHb

            F = a2 * PhiEq**alpha  - a1 * ptilde - a0
            viscoterm = etaOverDt * GAMMA

            IF ((GAMMA .GT. 0.D0) .AND. (etaOverDt .GT. 0.D0)) THEN
              F = F - viscoterm ** p_exp
            END IF

            iter_G = iter_G + 1

          END DO

          

          gma = gma_i

          ! WRITE(*,*) "Iter : ", iter_G
          ! WRITE(*,*) "GAMMA : ", GAMMA
          ! WRITE(*,*) "RESI : ", ABS(F)
          ! WRITE(*,*) "MAX_i : ", MAX_i
          ! WRITE(*,*) "-------------------------------------------"
            
        END SUBROUTINE nlinSolver