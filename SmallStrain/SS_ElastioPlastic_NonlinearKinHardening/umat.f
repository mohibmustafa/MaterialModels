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
      DIMENSION STRAIN_NP(6), SIGMA_NP_TR(6), MOD_NP_TR(6, 6)
      INTEGER K1, K2, K3, newton_iter
      DOUBLE PRECISION treps, mod_beta_np_tr, CHI_NP_TR, GAMMA_NP
      DOUBLE PRECISION CC1, CC2, mod_n, resi, local_tan, over_stress
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1          FOUR=4.D0, ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)

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
      P_Y0=PROPS(3)
      P_Yinf=PROPS(4)
      P_OMEGA=PROPS(5) 

C     Get Strain at time step n+1
      DO K1 = 1, 6
        STRAIN_NP(K1) = DSTRAN(K1) + STRAN(K1)
      END DO


c     Trace of strain at time step n+1    
      treps = STRAIN_NP(1) + STRAIN_NP(2) + STRAIN_NP(3)

c     Deviator of strain at time step n+1
      DO K1 = 1, 3
        dev(K1)= STRAIN_NP(K1) - treps / THREE
        dev(K1 + 3) = STRAIN_NP(K1 + 3) / TWO
      END DO
      
c     Calculate Beta trial
      DO K1 = 1, 6
        BETA_NP_TR(K1) = TWO * P_MU * (dev(K1) - STATEV(K1))    
      END DO

c     Calculate mod of beta trial
      mod_beta_np_tr = ZERO
      DO K1 = 1, 3
        mod_beta_np_tr = mod_beta_np_tr + BETA_NP_TR(K1)**TWO
        mod_beta_np_tr = mod_beta_np_tr + 2 * BETA_NP_TR(K1 +3)**TWO
      END DO
      mod_beta_np_tr = (mod_beta_np_tr)**0.5D0

C     Elastic Predictor step

c     Trial stress
      DO K1 = 1, 3
        SIGMA_NP_TR(K1) = P_KA * treps + BETA_NP_TR(K1)
        SIGMA_NP_TR(K1 + 3) = BETA_NP_TR(K1 + 3)
      END DO

c     trial modulus
      DO K1 = 1, 6
        DO K2 = 1, 6
          MOD_NP_TR(K1, K2) = P_KA * xioi(K1, K2) 
     1                      + TWO * P_MU * xpp(K1, K2)
        END DO 
      END DO

c     Trial yield function
      over_stress = ((TWO/THREE)**0.5D0) 
     1          * (P_Y0 + (P_Yinf - P_Y0) 
     2          * (ONE - EXP(-P_OMEGA * STATEV(7))))

C     write(*,*) "Radius of yield surface : ", over_stress
      CHI_NP_TR = mod_beta_np_tr - over_stress
      
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
C       write(*, *) "Plastic Step Chi: ", CHI_NP_TR
c       Calculate N_NP
        DO K1 = 1, 6
          ret_vec(K1) = BETA_NP_TR(K1) /  mod_beta_np_tr 
      
        END DO

c       Calculate N_NP diadic N_NP
        DO K1 = 1, 6
          DO K2 = 1, 6
           nponp(K1, K2) = ret_vec(K1) *  ret_vec(K2)
          END DO
        END DO

C       Newton Raphson Loop to get Gamma 
        GAMMA_NP = ZERO

        resi = mod_beta_np_tr - TWO * P_MU * GAMMA_NP
     1       - ((TWO / THREE) ** 0.5D0)
     2       * (P_Y0 + (P_Yinf - P_Y0) * (ONE 
     3       - EXP(-P_OMEGA * (STATEV(7) + GAMMA_NP 
     4       * (TWO / THREE)**0.5D0))))
        newton_iter = 0
        DO WHILE(resi .GT. TOLER)

          local_tan = TWO * P_MU + (TWO / THREE) * (P_Yinf - P_Y0)
     1              * P_OMEGA * EXP( -P_OMEGA * STATEV(7) - P_OMEGA 
     2              * ((TWO / THREE)**0.5D0) * GAMMA_NP)
          
          GAMMA_NP = GAMMA_NP + (1 / local_tan) * resi

          resi = mod_beta_np_tr - TWO * P_MU * GAMMA_NP
     1         - ((TWO / THREE) ** 0.5D0)
     2         * (P_Y0 + (P_Yinf - P_Y0) * (ONE 
     3         - EXP(-P_OMEGA * (STATEV(7) + GAMMA_NP 
     4         * (TWO / THREE)**0.5D0))))

          newton_iter = newton_iter + 1
C         write(*, *) "Newton iteration : ", newton_iter
C         write(*, *) "local residual : ", resi
        END DO
      
C       Calculate intermediate values
        CC1 = 1 / local_tan
        CC2 = GAMMA_NP / mod_beta_np_tr


C       Update Plastic Strain
        DO K1 = 1, 6
          STATEV(K1) = STATEV(K1) + GAMMA_NP * ret_vec(K1) 
        END DO

C       Update Plastic Arc Length
        STATEV(7) = STATEV(7) + ((TWO / THREE)**0.5D0) * GAMMA_NP

C       Correct Beta
        DO K1 = 1, 6
          BETA_NP(K1) = BETA_NP_TR(K1)
     1                - TWO * P_MU * GAMMA_NP * ret_vec(K1)      
        END DO


C       Correct Stress
        DO K1 = 1, 6
          STRESS(K1) = SIGMA_NP_TR(K1) 
     1               - TWO * P_MU * GAMMA_NP * ret_vec(K1)
        END DO

C       Correct Modulus
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = P_KA * xioi(K1, K2)
     1                     + TWO * P_MU * (ONE - TWO * P_MU  * CC2)
     2                     * xpp(K1, K2) + FOUR * (P_MU ** 2)
     3                     *(CC2 - CC1) * nponp(K1, K2)
          END DO
        END DO
      END IF

      RETURN
      END
