C UMAT - Linear Isotropic Viscoelastic
C Mohib Mustafa - IMDEA 2 Dec 2020
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,
     2           DPRED,CMNAME,
     3           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     3           DROT,PNEWDT,
     4           CELENT,DFGRD0,DFGRD1,NOEL,NPT,
     5           LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC' 

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),DFGRD0(3,3),
     4     DFGRD1(3,3)

C     List of internal variables:
      DIMENSION xioi(6, 6), xii(6, 6), xpp(6, 6), dev(6)
      DOUBLE PRECISION treps
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1          ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)

c     Define 4th order identity tensor
      data xioi/ ONE, ONE, ONE, ZERO, ZERO, ZERO,
     1           ONE, ONE, ONE, ZERO, ZERO, ZERO,
     2           ONE, ONE, ONE, ZERO, ZERO, ZERO,
     3           ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,
     4           ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,
     5           ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/

c     Define 4th order symmetric identity tensor
      data xii/   ONE, ZERO, ZERO, ZERO, ZERO, ZERO,
     1           ZERO, ONE, ZERO, ZERO, ZERO, ZERO,
     2           ZERO, ZERO, ONE, ZERO, ZERO, ZERO,
     3           ZERO, ZERO, ZERO, 0.5D0, ZERO, ZERO,
     4           ZERO, ZERO, ZERO, ZERO, 0.5D0, ZERO,
     5           ZERO, ZERO, ZERO, ZERO, ZERO, 0.5D0/

c     Compute deviatoric projection tensor
      DO K1 = 1, 6
        DO K2 = 1, 6
          xpp(K1, K2) = xii(K1, K2) - ONE / THREE * xioi(K1, K2)
        END DO
      END DO

c     Get material properties
      P_KA=PROPS(1)
      P_MU0=PROPS(2)
      P_MU1=PROPS(3)
      P_ETA1=PROPS(4)
      
c     Compute trace of strain tensor     
      treps = (DSTRAN(1) + DSTRAN(2) + DSTRAN(3)) + (STRAN(1) 
     1       + STRAN(2) + STRAN(3))

c     Compute deviatoric part of strain tensor
      DO K1 = 1, 3
            dev(K1)= DSTRAN(K1) + STRAN(K1) - treps / THREE
            dev(K1 + 3) = (DSTRAN(K1 + 3) + STRAN(K1 + 3)) / TWO
      END DO
      
c     Compute algorithmic parameters
      fac1 = TWO * P_MU1 * DTIME / P_ETA1
      fac2 = ONE/(ONE + fac1)

c     Update history variables
      DO K1 = 1,6
        STATEV(K1) = (STATEV(K1) + fac1 * dev(K1)) / (ONE + fac1)
      END DO

c     Compute Stress Tensor
      DO K1=1, 3
        STRESS(K1) = P_KA * treps + TWO * P_MU0 * dev(K1)
     1             + TWO * P_MU1 * (dev(K1) - STATEV(K1))
            
        STRESS(K1 + 3) = TWO * P_MU0 * dev(K1 + 3)
     1                 + TWO * P_MU1 * (dev(K1 + 3) - STATEV(K1 + 3)) 
      END DO

c     Compute Tangent Modulus
      DO K1=1, 6
        DO K2=1, 6
          DDSDDE(K1 ,K2) = P_KA * xioi(K1, K2) + TWO 
     1                   * (P_MU0 + P_MU1 * fac2) * xpp(K1, K2)
        END DO
      END DO

      RETURN
      END