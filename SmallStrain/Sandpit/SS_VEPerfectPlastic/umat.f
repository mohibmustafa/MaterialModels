C Isotropic linear visco elastic with perfect J2 Plasticity UMAT
C Mohib Mustafa - IMDEA 7 Jan 2021
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
      DIMENSION ret(6), dev(6), BETA_TR(6), BETA(6), BETA_P_TR(6)
      DIMENSION STRAIN(6), SIGMA_TR(6), TAN_TR(6, 6)
      DOUBLE PRECISION treps, mod_beta_p_tr, CHI_TR, gam
      DOUBLE PRECISION fac1, fac2, fac3, fac4, mod_ret, eff
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

      eff = (TWO * P_MU1) / (ONE + (TWO * P_MU1 * DTIME) / P_ETA1)
      fac1 = TWO * P_MU0 + eff

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
      
c     Calculate Beta trial - VE
      DO K1 = 1, 6
        BETA_TR(K1) = eff * (dev(K1) - STATEV(K1) - STATEV(K1 + 6))    
      END DO

c     Calculate Beta trial - Plastic
      DO K1 = 1, 6
        BETA_P_TR(K1) = TWO * P_MU0 * (dev(K1) - STATEV(K1 + 6))
     1                + BETA_TR(K1)    
      END DO

c     Trial stress
      DO K1 = 1, 3
        SIGMA_TR(K1) = P_KA * treps + BETA_P_TR(K1)
        SIGMA_TR(K1 + 3) = BETA_P_TR(K1 + 3)
      END DO

c     Trial Consistent Tangent
      DO K1 = 1, 6
        DO K2 = 1, 6
          TAN_TR(K1, K2) = P_KA * xioi(K1, K2) + fac1 * xpp(K1, K2)
        END DO 
      END DO

c     Calculate norm of Beta trial - Plastic
      mod_beta_p_tr = ZERO
      DO K1 = 1, 3
        mod_beta_p_tr = mod_beta_p_tr + BETA_P_TR(K1)**TWO
        mod_beta_p_tr = mod_beta_p_tr + 2 * BETA_P_TR(K1 + 3)**TWO
      END DO
      mod_beta_p_tr = (mod_beta_p_tr)**0.5D0

c     Trial yield function
      CHI_TR = mod_beta_p_tr - P_Y0 * (TWO/THREE)**0.5D0

C     Viscoelastic Predictor step
      IF(CHI_TR .LE. ZERO) THEN

C       Update VE Internal Variable
        DO K1 = 1, 6
          STATEV(K1) = STATEV(K1) + (DTIME / P_ETA1) * BETA_TR(K1)
        END DO

C       Update Stress 
        DO K1 = 1, 6
            STRESS(K1) = SIGMA_TR(K1)
        END DO
C       Update Modulus
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = TAN_TR(K1, K2)
          END DO
        END DO

C     Plastic Corrector Step
      ELSE
        write(*, *) "Plastic Step: Chi_TR = ", CHI_TR
C       Calculate plastic strain parameter
        gam = CHI_TR / fac1

c       Calculate N
        DO K1 = 1, 6
          ret(K1) = BETA_P_TR(K1) /  mod_beta_p_tr 
        END DO

c       Calculate N diadic N
        DO K1 = 1, 6
          DO K2 = 1, 6
            retoret(K1, K2) = ret(K1) *  ret(K2)
          END DO
        END DO

C       Correct Beta
        DO K1 = 1, 6
          BETA(K1) = BETA_TR(K1) - eff * gam * ret(K1)      
        END DO

C       Update internal variables
        DO K1 = 1, 6
            STATEV(K1) = STATEV(K1) + (DTIME / P_ETA1) * BETA(K1)
            STATEV(K1 + 6) = STATEV(K1 + 6) + gam * ret(K1) 
        END DO

C       Correct Stress
        DO K1 = 1, 6
          STRESS(K1) = SIGMA_TR(K1) - fac1 * gam * ret(K1)
        END DO

C       Correct Consistent Tangent
        fac2 = fac1 / mod_beta_p_tr
        fac4 = gam * fac2
        fac3 = ONE - fac4
        
        DO K1 = 1, 6
          DO K2 = 1, 6
            DDSDDE(K1, K2) = TAN_TR(K1, K2) 
     1                     - fac1 * (fac3 * retoret(K1, K2)
     2                     + fac4 * xpp(K1, K2))
          END DO
        END DO
      END IF

      RETURN
      END
