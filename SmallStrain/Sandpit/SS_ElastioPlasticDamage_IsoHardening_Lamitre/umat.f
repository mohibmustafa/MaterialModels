C UMAT - Linear Isotropic Elastic Plastic Iso Hardening Lamitre Damage
C Mohib Mustafa - IMDEA 28 Jan 2021      
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
      DIMENSION ret(6), dev(6), BETA_TR(6), BETA(6)
      DIMENSION STRAIN(6), SIGMA_TR(6), TANG_TR(6, 6)
      DOUBLE PRECISION treps, mod_beta_tr, CHI_TR, GAM, d_arc
      DOUBLE PRECISION BB1, BB2, mod_n, xi_trial, my, sq_eps_e
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
      P_MU=PROPS(2)
      P_Y0=PROPS(3)
      P_h=PROPS(4)
      P_SS=PROPS(5)
      P_S=PROPS(6)

C     Get Strain at time step n+1
      DO K1 = 1, 6
        STRAIN(K1) = DSTRAN(K1) + STRAN(K1)
      END DO

c     Trace of strain at time step n+1    
      treps = STRAIN(1) + STRAIN(2) + STRAIN(3)

c     Deviator of strain at time step n+1
      DO K1 = 1, 3
        dev(K1)= STRAIN(K1) - treps / THREE
        dev(K1 + 3) = STRAIN(K1 + 3) / TWO
      END DO
      
c     Calculate Beta trial
      DO K1 = 1, 6
        BETA_TR(K1) = TWO * P_MU * (dev(K1) - STATEV(K1))    
      END DO

c     Calculate mod of beta trial
      mod_beta_tr = ZERO
      DO K1 = 1, 3
        mod_beta_tr = mod_beta_tr + BETA_TR(K1)**TWO
        mod_beta_tr = mod_beta_tr + 2 * BETA_TR(K1 + 3)**TWO
      END DO
      mod_beta_tr = (mod_beta_tr)**0.5D0

C     Calculate xi trial
      xi_trial = P_h * STATEV(7)


C     Elastic Predictor step

c     Trial stress
      DO K1 = 1, 3
        SIGMA_TR(K1) = P_KA * treps + BETA_TR(K1)
        SIGMA_TR(K1 + 3) = BETA_TR(K1 + 3)
      END DO

c     Trial tangent modulus
      DO K1 = 1, 6
        DO K2 = 1, 6
          TANG_TR(K1, K2) = P_KA * xioi(K1, K2) 
     1                      + TWO * P_MU * xpp(K1, K2)
        END DO 
      END DO

c     Trial yield function
      CHI_TR = mod_beta_tr - (P_Y0 + xi_trial) * (TWO/THREE)**0.5D0
      
C     Check for plastic loading
      IF(CHI_TR .LE. ZERO) THEN

C       Update Stress 
        DO K1 = 1, 6
            STRESS(K1) = SIGMA_TR(K1) * (ONE - STATEV(8))
        END DO
C       Update tangent modulus
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = TANG_TR(K1, K2) * (ONE - STATEV(8))
          END DO
        END DO

C     Plastic Corrector Step
      ELSE
C       Calculate plastic strain parameter
        GAM = CHI_TR / (TWO * P_MU + (TWO / THREE) * P_h)
        
c       Calculate N_NP
        DO K1 = 1, 6
          ret(K1) = BETA_TR(K1) /  mod_beta_tr       
        END DO
      
c       Calculate mod of n
        mod_n = ZERO
        DO K1 = 1, 3
          mod_n = mod_n + ret(K1)**TWO
          mod_n = mod_n + 2 * ret(K1 +3)**TWO
        END DO
        mod_n = (mod_n)**0.5D0

c      Calculate N_NP diadic N_NP
       DO K1 = 1, 6
         DO K2 = 1, 6
           retoret(K1, K2) = ret(K1) * ret(K2)
         END DO
      END DO

C       Update Plastic Strain
        DO K1 = 1, 6
            STATEV(K1) = STATEV(K1) + GAM * ret(K1) 
        END DO

C      Update Plastic arc length
        d_arc = GAM * (TWO / THREE)**0.5D0
        STATEV(7) = STATEV(7) + d_arc

C     Calculate -Y (Driving force for damage parameter)   
        sq_eps_e = 0.0D0
        DO K1 = 1, 6
          sq_eps_e = sq_eps_e +  (dev(K1) - STATEV(K1))**TWO
        END DO

        my = 0.5D0 * P_KA * treps**TWO + P_MU * sq_eps_e 
     1       + 0.5D0 * P_h * STATEV(7)**TWO 

C        write(*,*) "-Y : ", my
C        write(*,*) "S0 : ", P_SS
C        write(*,*) "s0 : ", P_S
C        write(*,*) "del arc : ", d_arc
C     Update D
        STATEV(8) = STATEV(8) + ((my/P_SS)**P_S) * d_arc

        IF (STATEV(8) .GT. ONE) THEN
          STATEV(8) = ONE
        END IF

        write(*,*) "D : ", STATEV(8)

C     Correct Beta
        DO K1 = 1, 6
          BETA(K1) = BETA_TR(K1) - TWO * P_MU * GAM * ret(K1)
        END DO

C       Correct Stress
        DO K1 = 1, 3
            STRESS(K1) = (P_KA * treps + BETA(K1)) * (ONE - STATEV(8))
            STRESS(K1 + 3) = BETA_TR(K1 + 3) * (ONE - STATEV(8))
        END DO

C       Correct tangent modulus
        BB2 = (TWO * P_MU * GAM) / mod_beta_tr
        BB1 = (TWO * P_MU) / ((TWO / THREE) * P_h + TWO * P_MU) - BB2
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = (TANG_TR(K1, K2) 
     1                     - TWO * P_MU * (BB1 * retoret(K1, K2)
     2                     + BB2 * xpp(K1, K2))) * (ONE - STATEV(8))
          END DO
        END DO
      END IF

      RETURN
      END
