C Isotropic non linear J2 Plasticity UMAT
C           
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

C     List of internal variables:
      DIMENSION xioi(6, 6), xii(6, 6), xpp(6, 6), nponp(6, 6)
      DIMENSION ret_vec(6), dev(6), BETA_NP_TR(6), BETA_NP(6)
      DIMENSION SIGMA_dev_NP_TR(6), SIGMA_dev_NP(6), xi(6)
      DIMENSION STRAIN_NP(6), SIGMA_NP_TR(6), MOD_NP_TR(6, 6)
      INTEGER K1, K2, K3, newton_iter
      DOUBLE PRECISION treps, mod_beta_np_tr, CHI_NP_TR, GAM_i
      DOUBLE PRECISION BB1, BB2, mod_n, resi, local_tan, over_stress
      DOUBLE PRECISION beta_np_sc_tr, mod_xi_tr, rap_iter, res_i
      DOUBLE PRECISION fac1, fac2, fac3 
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1          FOUR=4.D0, ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-9)

c     define 4th order identity tensor
      data xioi/ ONE, ONE, ONE, ZERO, ZERO, ZERO,
     1           ONE, ONE, ONE, ZERO, ZERO, ZERO,
     2           ONE, ONE, ONE, ZERO, ZERO, ZERO,
     3           ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,
     4           ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,
     5           ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/

c     define 4th order symmetric identity tensor
      data xii/   ONE, ZERO, ZERO, ZERO, ZERO, ZERO,
     1           ZERO, ONE, ZERO, ZERO, ZERO, ZERO,
     2           ZERO, ZERO, ONE, ZERO, ZERO, ZERO,
     3           ZERO, ZERO, ZERO, 0.5D0, ZERO, ZERO,
     4           ZERO, ZERO, ZERO, ZERO, 0.5D0, ZERO,
     5           ZERO, ZERO, ZERO, ZERO, ZERO, 0.5D0/

c     compute deviatoric projection tensor
      DO K1 = 1, 6
       DO K2 = 1, 6
        xpp(K1, K2) = xii(K1, K2) - (ONE / THREE) * xioi(K1, K2)
       END DO
      END DO

c     Get material properties
      P_KA=PROPS(1)
      P_MU=PROPS(2)
      P_h=PROPS(3)
      P_HH=PROPS(4)
      P_Y0=PROPS(5)
      P_eta=PROPS(6)
      P_m=PROPS(7)

      fac2 = P_eta / (DTIME * P_Y0)
      fac3 = (TWO * P_MU + P_HH + (TWO/THREE) * P_h) / P_Y0
    
C      write(*,*) "Total time mine", TIME(2) + DTIME
C      write(*,*) "dt : ", DTIME

C     Get Strain at time step n+1
      DO K1 = 1, 6
        STRAIN_NP(K1) = DSTRAN(K1) + STRAN(K1)

        write(*,*) "EPSILON(", K1, ") : ", STRAIN_NP(K1) 
      END DO

C      Write(*,*) "------------------------------------------"

c     Trace of strain at time step n+1    
      treps = STRAIN_NP(1) + STRAIN_NP(2) + STRAIN_NP(3)

c     Deviator of strain at time step n+1
      DO K1 = 1, 3
        dev(K1)= STRAIN_NP(K1) - treps / THREE
        dev(K1 + 3) = STRAIN_NP(K1 + 3) / TWO
      END DO
      
c     Calculate sigma_bar trial and beta_trial
      DO K1 = 1, 6
        SIGMA_dev_NP_TR(K1) = TWO * P_MU * (dev(K1) - STATEV(K1))
        BETA_NP_TR(K1) = P_HH * STATEV(K1)     
      END DO

C     Calculate sigma_bar trial and beta_trial
      DO K1 = 1, 6
        xi(K1) = SIGMA_dev_NP_TR(K1) - BETA_NP_TR(K1)
      END DO

c     beta trial scaler
      beta_np_sc_tr = P_h * STATEV(7)

c     Calculate mod of trial forces
      mod_xi_tr = ZERO
      DO K1 = 1, 3
        mod_xi_tr = mod_xi_tr + xi(K1)**TWO   
        mod_xi_tr = mod_xi_tr + 2 * xi(K1 + 3)**TWO
      END DO
      mod_xi_tr = (mod_xi_tr)**0.5D0

c     Calculate N_NP
      DO K1 = 1, 6
        ret_vec(K1) = xi(K1) / mod_xi_tr 
      END DO

c     Calculate N_NP diadic N_NP
      DO K1 = 1, 6
        DO K2 = 1, 6
          nponp(K1, K2) = ret_vec(K1) *  ret_vec(K2)
        END DO
      END DO

C     Elastic Predictor step
c     Trial stress
      DO K1 = 1, 3
        SIGMA_NP_TR(K1) = P_KA * treps + SIGMA_dev_NP_TR(K1)
        SIGMA_NP_TR(K1 + 3) = SIGMA_dev_NP_TR(K1 + 3)
      END DO

c     trial modulus
      DO K1 = 1, 6
        DO K2 = 1, 6
          MOD_NP_TR(K1, K2) = P_KA * xioi(K1, K2) 
     1                      + TWO * P_MU * xpp(K1, K2)
        END DO 
      END DO

c     Trial yield function
C     write(*,*) "Radius of yield surface : ", over_stress
      CHI_NP_TR = mod_xi_tr - ((TWO / THREE)**0.5D0)
     1          * (P_Y0 + beta_np_sc_tr)

C     Check for plastic loading
      IF(CHI_NP_TR .LT. ZERO) THEN
        
C        write(*, *) "Elastic Step Chi: ", CHI_NP_TR
C       Update Internal Variables
        DO K1 = 1, 7
          STATEV(K1) = STATEV(K1)
        END DO
C       Update Stress 
        DO K1 = 1, 6
            STRESS(K1) = SIGMA_NP_TR(K1)
        END DO
C       Update Modulus
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = MOD_NP_TR(K1, K2)
          END DO
        END DO

C     Plastic Corrector Step
      ELSE
      
C      write(*, *) "Timestep : ", TIME(2)
      
C    ----------------------------------------------------------------
C                  NEWTON RAPHSON SOLVER FOR GAM_i
C    ----------------------------------------------------------------

C      Initialize counter and initial values
        rap_iter = 0
        GAM_i = ZERO
 
C      Set residual to trial value of elastic domain boundary

        fac1 = (ONE / P_Y0) * CHI_NP_TR - fac3 * GAM_i

        res_i = fac1 ** P_m - fac2 * GAM_i

        G_i = P_m * fac3 * (fac1) ** (P_m - ONE) + fac2    
        
C        Write(*,*) "Gam_i before first iter", GAM_i
C        write(*,*) "Res before 1st iter : ", res_i
C        write(*,*) "G before 1st iter : ", G_i

        DO WHILE(res_i .GT. TOLER)
          rap_iter = rap_iter + 1
          GAM_i = GAM_i + (ONE/G_i)*res_i

          fac1 = (ONE / P_Y0) * CHI_NP_TR - fac3 * GAM_i
            
          G_i = P_m * fac3 * (fac1) ** (P_m - ONE) + fac2
    
          res_i = fac1 ** P_m - fac2 * GAM_i
          
C          write(*,*) "Gam_i at ", rap_iter, " : ", GAM_i
C          write(*,*) "Res at ", rap_iter, " : ", res_i
C          write(*,*) "G at ", rap_iter, " : ", G_i

        END DO
        
C        write(*, *) "Chi Trial : ", CHI_NP_TR
C        write(*, *) "Gam_i : ", GAM_i

C       Calculate intermediate values
        BB1 = ONE - TWO * (P_MU * GAM_i) / mod_xi_tr
        BB2 = ((TWO * P_MU * P_m) / (G_i * P_Y0)) * fac1 ** (P_m - ONE) 
     3      - (TWO * P_MU * GAM_i) / mod_xi_tr

C       Update Plastic Strain
        DO K1 = 1, 6
          STATEV(K1) = STATEV(K1) + GAM_i * ret_vec(K1) 
        END DO

C       Update Plastic Arc Length
        STATEV(7) = STATEV(7) + ((TWO / THREE)**0.5D0) * GAM_i

C       Correct Beta
        DO K1 = 1, 6
          SIGMA_dev_NP(K1) = SIGMA_dev_NP_TR(K1)
     1                - TWO * P_MU * GAM_i * ret_vec(K1)      
        END DO

C       Correct Stress
        DO K1 = 1, 6
          STRESS(K1) = SIGMA_NP_TR(K1) 
     1               - TWO * P_MU * GAM_i * ret_vec(K1)
        END DO

C       Correct Modulus
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = P_KA * xioi(K1, K2)
     1                     + TWO * P_MU * BB1 * xpp(K1, K2)
     2                     - TWO * P_MU * BB2 * nponp(K1, K2)
          END DO
        END DO
      END IF

      RETURN
      END
