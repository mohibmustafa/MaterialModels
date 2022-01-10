C Mohib Mustafa - IMDEA 23 Feb 2021      
C UMAT - Linear Isotropic Viscoelastic Viscoplastic - 2 Maxwell Branches
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
      DIMENSION ret(6), dev(6), BETA_TR(6, 2), BETA(6, 2)
      DIMENSION BETA_VP_TR(6), eps_ve_trial(6), BETA_TR_SUM(6)
      DIMENSION BETA_VP_PAR_TR(6)
      DIMENSION STRAIN(6), SIGMA_TR(6), TANG_TR(6, 6)
      DOUBLE PRECISION treps, mod_beta_vp_tr, CHI_TR, gam
      DOUBLE PRECISION fac4, zeta_trial
      DIMENSION fac1(2), ETA(2), xMU(2)
      DOUBLE PRECISION fac3, fac5, fac6
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
      xMU(1)=PROPS(3)
      xMU(2)=PROPS(4)
      ETA(1)=PROPS(5)
      ETA(2)=PROPS(6)
      P_Y0=PROPS(7)
      P_H=PROPS(8)
      P_HH=PROPS(9)
      P_ETA0=PROPS(10)

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

C     Difference of dev strain tensor at t n+1 and vp strain at n  
      DO K1 = 1, 6
        eps_ve_trial(K1)= dev(K1) - STATEV(12 + K1)
      END DO

C     Algorithmic factors
      DO K1 = 1, 2
        fac1(K1) = (TWO * xMU(K1)) / (ONE + (TWO * xMU(K1) * DTIME) / ETA(K1))
      END DO
      
      fac3 = TWO * P_MU0 + fac1(1) + fac1(2)
      
C     Calculate Beta trial
      x_iter = 0
      DO K2 = 1, 2
        DO K1 = 1, 6
          BETA_TR(K1, K2) = fac1(K2) * (eps_ve_trial(K1) - STATEV(K1 + x_iter))
        END DO
        x_iter = x_iter + 6
      END DO

C     Calculate Beta trial Sum
      DO K1 = 1, 6
        BETA_TR_SUM(K1) = ZERO
      END DO

      DO K2 = 1, 2
        DO K1 = 1, 6
          BETA_TR_SUM(K1) = BETA_TR_SUM(K1) + BETA_TR(K1, K2)
        END DO 
      END DO

C     Calculate Beta VP trial, Partial
      DO K1 = 1, 6
        BETA_VP_PAR_TR(K1) = TWO * P_MU0 * eps_ve_trial(K1) + BETA_TR_SUM(K1)
      END DO      

C     Calculate Beta VP trial
      DO K1 = 1, 6
        BETA_VP_TR(K1) = BETA_VP_PAR_TR(K1) - P_HH * STATEV(K1 + 12)
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
      zeta_trial = P_H * STATEV(19)

C     Elastic Predictor step

C     Trial stress
      DO K1 = 1, 3
        SIGMA_TR(K1) = P_KA * treps + BETA_VP_PAR_TR(K1)
        SIGMA_TR(K1 + 3) = BETA_VP_PAR_TR(K1 + 3)
      END DO

C     Trial tangent modulus
      DO K1 = 1, 6
        DO K2 = 1, 6
          TANG_TR(K1, K2) = P_KA * xioi(K1, K2) + fac3 * xpp(K1, K2)
        END DO 
      END DO

C     Trial yield function
      CHI_TR = mod_beta_vp_tr 
     1       - (P_Y0 + zeta_trial) * (TWO/THREE)**0.5D0
      
C      write(*,*) "CHI_TR :", CHI_TR
C     Check for plastic loading
      IF(CHI_TR .LT. ZERO) THEN

C       Update Internal variable for viscoelasticity
        x_iter = 0
        DO K2 = 1, 2
          DO K1 = 1, 6
            STATEV(K1 + x_iter) = STATEV(K1 + x_iter) 
     1                          + (DTIME / ETA(K2)) * BETA_TR(K1, K2)
          END DO
          x_iter = x_iter + 6
        END DO

C       Update Stress 
        DO K1 = 1, 6
            STRESS(K1) = SIGMA_TR(K1)
        END DO

C       Update tangent modulus
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = TANG_TR(K1, K2)
          END DO
        END DO

C     Plastic Corrector Step
      ELSE
C       Calculate plastic strain parameter
        gam = CHI_TR / 
     1        (fac3 + P_HH + (TWO / THREE)*P_H + (P_ETA0 / DTIME))

C       Calculate N_NP
        DO K1 = 1, 6
          ret(K1) = BETA_VP_TR(K1) /  mod_beta_vp_tr       
        END DO

c      Calculate N_NP diadic N_NP
       DO K1 = 1, 6
         DO K2 = 1, 6
           retoret(K1, K2) = ret(K1) * ret(K2)
         END DO
       END DO

C      Correct Beta
       DO K2 = 1, 2
         DO K1 = 1, 6
           BETA(K1, K2) = BETA_TR(K1, K2) - fac1(K2) * gam * ret(K1)      
         END DO
       END DO

C       Update Internal Variables
        x_iter = 0
        DO K2 = 1, 2
          DO K1 = 1, 6
            STATEV(K1 + x_iter) = STATEV(K1 + x_iter) 
     1                          + (DTIME/ETA(K2)) * BETA(K1, K2)             
          END DO
          x_iter = x_iter + 6
        END DO

        DO K1 = 1, 6
          STATEV(K1 + 12) = STATEV(K1 + 12) + gam * ret(K1)
        END DO

        STATEV(19) = STATEV(19) + gam * (TWO / THREE)**0.5D0

C       Correct Stress
        DO K1 = 1, 6
            STRESS(K1) = SIGMA_TR(K1) - fac3 * gam * ret(K1)     
        END DO

C       Correct tangent modulus
        fac4 = fac3 / mod_beta_vp_tr
        fac6 = gam * fac4
        fac5 = ((fac3) / (fac3 + P_HH + (TWO/THREE) * P_H 
     1       + P_ETA0/DTIME)) - fac6
       
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = TANG_TR(K1, K2) - fac3 
     1                      * (fac5 * retoret(K1, K2) + fac6 *xpp(K1, K2))
          END DO
        END DO
      END IF

      RETURN
      END