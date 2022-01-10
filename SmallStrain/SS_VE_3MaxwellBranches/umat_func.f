C Mohib Mustafa - IMDEA 4 MAR 2021      
C UMAT - Isotropic Viscoplastic - Perznya with M Exponent
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


      call VE3(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,NPROPS,
     &         DTIME,DSTRAN,KINC,KSTEP,NOEL,TIME)

      return
      end      

      SUBROUTINE VE3(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &        NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL, TIME)
     
     
       IMPLICIT NONE
      
       INTEGER, PARAMETER :: double=kind(1.d0)
       INTEGER, PARAMETER :: prec=double
       
       INTEGER, INTENT(IN)      :: ntens,nprops,nstatv,kinc,kstep,noel
       REAL(prec), INTENT(IN)   :: stran(ntens), dstran(ntens)
       REAL(prec), INTENT(IN)   :: props(nprops), dtime, time(2)
       REAL(prec), INTENT(INOUT) :: statev(nstatv)
       REAL(prec), INTENT(OUT)   :: stress(ntens),ddsdde(ntens,ntens)
      
      
C      List of internal variables:
       INTEGER    :: K1, K2
       REAL(prec) :: xioi(6, 6), xii(6, 6), xpp(6, 6)
       REAL(prec) :: ret(6), dev(6), BETA_VE(6, 3)
       REAL(prec) :: SIGMA_BAR(6), alpha(6, 3)
       REAL(prec) :: EPS(6), fac_ve(3)
       REAL(prec) :: treps
       REAL(prec) :: ka, mu0, mu(3), eta(3)
       
       REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
       REAL(prec), PARAMETER :: THREE=3.D0, SIX=6.D0, ENUMAX=.4999D0
       REAL(prec), PARAMETER :: NEWTON=10, TOLER=1.0D-6

C      Define 4th order identity tensor
       data xioi(1,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
       data xioi(2,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
       data xioi(3,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
       data xioi(4,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
       data xioi(5,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
       data xioi(6,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/

C      Define 4th order symmetric identity tensor
       data xii(1, :) /ONE, ZERO, ZERO, ZERO, ZERO, ZERO/
       data xii(2, :) /ZERO, ONE, ZERO, ZERO, ZERO, ZERO/
       data xii(3, :) /ZERO, ZERO, ONE, ZERO, ZERO, ZERO/
       data xii(4, :) /ZERO, ZERO, ZERO, 0.5D0, ZERO, ZERO/
       data xii(5, :) /ZERO, ZERO, ZERO, ZERO, 0.5D0, ZERO/
       data xii(6, :) /ZERO, ZERO, ZERO, ZERO, ZERO, 0.5D0/
       
C     Compute deviatoric projection tensor
       xpp(:, :) = xii(:, :) - (ONE / THREE) * xioi(:, :)
         
C     Get material properties
      ka=PROPS(1)
      mu0=PROPS(2)
      mu(1)=PROPS(3)
      mu(2)=PROPS(4)
      mu(3)=PROPS(5)
      eta(1)=PROPS(6)
      eta(2)=PROPS(7)
      eta(3)=PROPS(8)
      
C     DEFINE STATE VARIABLES      
      EPS(:) = DSTRAN(:) + STRAN(:)
      alpha(:, 1) = STATEV(1 : 6)
      alpha(:, 2) = STATEV(7 : 12)
      alpha(:, 3) = STATEV(13 : 18)
            
C     Trace of strain at time step n+1    
      treps = EPS(1) + EPS(2) + EPS(3)
      
C     Deviator of strain at time step n+1
      dev(1 : 3) = EPS(1 : 3) - treps / THREE
      dev(4 : 6) = EPS(4 : 6) / TWO

C     Define algorithmic constants
      fac_ve(:) = (TWO * mu(:)) / (ONE + (TWO * mu(:) * DTIME) / eta(:))
      
C      write(*,*)"fac1", fac_ve(1)
C      write(*,*)"fac2", fac_ve(2)
C      write(*,*)"fac3", fac_ve(3)

C     Calculate Beta_VE
      DO K1 = 1, 3 
        BETA_VE(:, K1) = fac_ve(K1) * (dev(:) - alpha(:, K1)) 
      END DO
      
        
C     Update Internal Variables
      STATEV(1 : 6) = STATEV(1 : 6) + (DTIME/eta(1)) * BETA_VE(:, 1)
      STATEV(7 : 12) = STATEV(7 : 12) + (DTIME/eta(2)) * BETA_VE(:, 2)
      STATEV(13 : 18) = STATEV(13 : 18) + (DTIME/eta(3)) * BETA_VE(:, 3)

C     Calculate Stress
      STRESS(:) = ka * treps + TWO * mu0 * dev(:) + BETA_VE(:,1) 
     1          + BETA_VE(:,2) + BETA_VE(:,3)
              
C     Calculate Consistent Tangent
      DDSDDE(:, :) = ka * xioi(:, :)
     1              + (TWO * mu0 + fac_ve(1) + fac_ve(2) 
     2              + fac_ve(3)) * xpp(:, :)
      
      END SUBROUTINE VE3
      
