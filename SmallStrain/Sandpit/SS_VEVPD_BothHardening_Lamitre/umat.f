C UMAT - Linear Isotropic Viscoelastic Viscoplastic including Lamitre Damage
C Mohib Mustafa - IMDEA 8 Feb 2021      
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
      DIMENSION xioi(6, 6), xii(6, 6), xpp(6, 6), retoret(6, 6)
      DIMENSION ret(6), dev(6), BETA_TR(6), BETA(6), BETA_VP_TR(6)
      DIMENSION BETA_VP_PAR_TR(6), eps_vp_i(6), alpha_i(6)
      DIMENSION STRAIN(6), SIGMA_TR(6), TANG_TR(6, 6), E1(6), E2(6)
      DOUBLE PRECISION treps, mod_beta_vp_tr, CHI_TR, GAM_i, d_arc
      DOUBLE PRECISION fac2, fac4, mod_n, zeta_trial, y_i, sq_eps_e
      DOUBLE PRECISION fac1, fac5, res_i, xi_i, ee1, ee2, evp, G_i
      DOUBLE PRECISION fac3
      INTEGER rap_iter
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1          ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)

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
      P_MU0=PROPS(2)
      P_MU1=PROPS(3)
      P_ETA1=PROPS(4)
      P_Y0=PROPS(5)
      P_H=PROPS(6)
      P_HH=PROPS(7)
      P_ETA0=PROPS(8)
      P_ESS=PROPS(9)
      P_ES=PROPS(10)

C      write(*,*) "S0 : ", P_ESS
C      write(*,*) "s0 : ", P_ES


C     Get Strain at time step n+1
      DO K1 = 1, 6
        STRAIN(K1) = DSTRAN(K1) + STRAN(K1)
      END DO

C     Trace of strain at time step n+1    
      treps = STRAIN(1) + STRAIN(2) + STRAIN(3)

C     Deviator of strain at time step n+1
      DO K1 = 1, 3
        dev(K1)= STRAIN(K1) - treps / THREE
        dev(K1 + 3) = STRAIN(K1 + 3) / TWO
      END DO

C     Algorithmic factors
      fac1 = (TWO * P_MU1) / (ONE + (TWO * P_MU1 * DTIME) / P_ETA1)
      fac5 = TWO * P_MU0 + fac1
      
C     Calculate Beta trial
      DO K1 = 1, 6
        BETA_TR(K1) = fac1 * (dev(K1) - STATEV(K1) - STATEV(K1 + 6))
      END DO

C     Calculate Beta VP trial, Partial
      DO K1 = 1, 6
        BETA_VP_PAR_TR(K1) = TWO * P_MU0 
     1           * (dev(K1) - STATEV(K1 + 6)) + BETA_TR(K1)
      END DO      

C     Calculate Beta VP trial
      DO K1 = 1, 6
        BETA_VP_TR(K1) = BETA_VP_PAR_TR(K1) - P_HH * STATEV(K1 + 6)
      END DO      

C     Calculate mod of Beta VP trial
      mod_beta_vp_tr = ZERO
      DO K1 = 1, 3
        mod_beta_vp_tr = mod_beta_vp_tr + BETA_VP_TR(K1)**TWO
        mod_beta_vp_tr = mod_beta_vp_tr 
     1                 + 2 * BETA_VP_TR(K1 + 3)**TWO
      END DO
      mod_beta_vp_tr = (mod_beta_vp_tr)**0.5D0

C     Calculate zeta trial
      zeta_trial = P_H * STATEV(13)

C     Elastic Predictor step

C     Trial stress
      DO K1 = 1, 3
        SIGMA_TR(K1) = P_KA * treps + BETA_VP_PAR_TR(K1)
        SIGMA_TR(K1 + 3) = BETA_VP_PAR_TR(K1 + 3)
      END DO

C     Trial tangent modulus
      DO K1 = 1, 6
        DO K2 = 1, 6
          TANG_TR(K1, K2) = P_KA * xioi(K1, K2) 
     1                      + fac5 * xpp(K1, K2)
        END DO 
      END DO

C     Trial yield function
      CHI_TR = mod_beta_vp_tr 
     1       - (P_Y0 + zeta_trial) * (TWO/THREE)**0.5D0
      
C      write(*,*) "CHI_TR :", CHI_TR
C     Check for plastic loading
      IF(CHI_TR .LT. ZERO) THEN

C       Update Internal variable for viscoelasticity
        DO K1 = 1, 6
          STATEV(K1) = STATEV(K1) + (DTIME / P_ETA1) * BETA_TR(K1)
        END DO

C       Update Stress 
        DO K1 = 1, 6
            STRESS(K1) = SIGMA_TR(K1) * (ONE - STATEV(14))
        END DO

C       Update tangent modulus
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = TANG_TR(K1, K2) * (ONE - STATEV(14))
          END DO
        END DO

C     Plastic Corrector Step
      ELSE

C       Calculate N_NP
        DO K1 = 1, 6
          ret(K1) = BETA_VP_TR(K1) /  mod_beta_vp_tr       
        END DO

C       Calculate mod of n (Only for debugging remove later !!!!)
C        mod_n = ZERO
C        DO K1 = 1, 3
C          mod_n = mod_n + ret(K1)**TWO
C          mod_n = mod_n + 2 * ret(K1 +3)**TWO
C        END DO
C        mod_n = (mod_n)**0.5D0

c      Calculate N_NP diadic N_NP
       DO K1 = 1, 6
       DO K2 = 1, 6
         retoret(K1, K2) = ret(K1) * ret(K2)
       END DO
       END DO

C       Newton Raphson loop for Plastic Parameter Gamma
        rap_iter = 0
        GAM_i = ZERO
C       Set residual to trial value of elastic domain boundary
        res_i = CHI_TR

        DO WHILE(res_i .GT. TOLER)
        
          DO K1 = 1, 6
            alpha_i(K1) = STATEV(K1) - (DTIME / P_ETA1) 
     1                  * fac1 * GAM_i * ret(K1)

            eps_vp_i(K1) = STATEV(K1 + 6) + GAM_i * ret(K1)
          END DO
          
          DO K1 = 1, 6
            alpha_i(K1) = STATEV(K1) - (DTIME / P_ETA1) 
     1                  * fac1 * GAM_i * ret(K1)

            eps_vp_i(K1) = STATEV(K1 + 6) + GAM_i * ret(K1)
          END DO


          xi_i = STATEV(13) + GAM_i * (TWO / THREE)**0.5D0

          DO K1 = 1, 6
            E1(K1) = dev(K1) - eps_vp_i(K1)
          END DO

          DO K1 = 1, 6
            E2(K1) = E1(K1) - alpha_i(K1)
          END DO

          ee1 = ZERO
          ee2 = ZERO
          evp = ZERO
          DO K1 = 1, 3
            ee1 = ee1 + E1(K1)**TWO
            ee2 = ee2 + E2(K1)**TWO
            evp = evp + eps_vp_i(K1)**TWO

            ee1 = ee1 + TWO * E1(K1 + 3)**TWO
            ee2 = ee2 + TWO * E2(K1 + 3)**TWO
            evp = evp + TWO * eps_vp_i(K1 + 3)**TWO
          END DO

          y_i = 0.5D0 * P_KA * treps**TWO + P_MU0 * ee1 
     1        + P_MU1 * ee2 + 0.5D0 * P_HH * evp 
     2        + 0.5D0 * P_H * xi_i**TWO

          G_i = fac5 + P_HH + (TWO/THREE)*P_H 
     1        + (P_ETA0/DTIME) * ((TWO/THREE)**0.5D0) 
     2        * ((y_i / P_ESS)**P_ES) * GAM_i
     3        + (P_ETA0/DTIME) * (ONE - STATEV(14) 
     4        - GAM_i*((TWO/THREE)**0.5D0) * ((y_i / P_ESS)**P_ES))

          res_i = CHI_TR - (fac5 + P_HH + (TWO/THREE)*P_H)*GAM_i 
     1          - (P_ETA0/DTIME) * (ONE - STATEV(14) 
     2          - GAM_i*((TWO/THREE)**0.5D0) 
     3          * ((y_i / P_ESS)**P_ES))*GAM_i
          
          GAM_i = GAM_i + (ONE/G_i)*res_i
          rap_iter = rap_iter + 1
        END DO

C       Update Internal Variables
        DO K1 = 1, 6
            STATEV(K1) = alpha_i(K1)
            STATEV(K1 + 6) = eps_vp_i(K1) 
        END DO
        STATEV(13) = xi_i
        STATEV(14) = STATEV(14) + GAM_i * ((TWO/THREE)**0.5D0) * (y_i/P_ESS)**P_ES

        IF (STATEV(14) .GT. ONE) THEN
          STATEV(14) = ONE
        END IF

C        write(*,*) "Guass Iter : ", rap_iter
C        write(*,*) "Residual : ", res_i
C        write(*,*) "D : ", STATEV(14)
C        write(*,*) "Y : ", y_i
C        write(*,*) "S0 : ", P_ESS
C        write(*,*) "s0 : ", P_ES


C       Correct Stress
        DO K1 = 1, 3
            STRESS(K1) = (SIGMA_TR(K1) 
     2                 - fac5 * GAM_i * ret(K1)) 
     3                 * (ONE - STATEV(14))
        END DO

C       Correct tangent modulus
        fac2 = fac5 / mod_beta_vp_tr
        fac4 = GAM_i * fac2
        fac3 = (fac5 / G_i) - fac4
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = (TANG_TR(K1, K2) 
     1                     - fac5 * (fac3 * retoret(K1, K2) 
     2                     + fac4 *xpp(K1, K2)))
     3                     * (ONE - STATEV(14))
          END DO
        END DO
      END IF

      RETURN
      END
