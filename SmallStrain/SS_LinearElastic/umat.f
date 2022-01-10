C     UMAT - Isotropic Linear Elastic
C     
C     Mohib Mustafa - IMDEA 5 Dec 2020      

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


      call LIN_ELAST(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,NPROPS,
     &         DTIME,DSTRAN,KINC,KSTEP,NOEL,TIME)

      return
      end    


C     !-------------------------------------------------------------------------
C     !          UMAT SUBROUTINE FOR ISOTROPIC LINEAR ELASTIC MATERIAL 
C     !-------------------------------------------------------------------------
      SUBROUTINE LIN_ELAST(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     1        NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL, TIME)

        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
       
        INTEGER, INTENT(IN)      :: ntens,nprops,nstatv,kinc,kstep,noel
        REAL(prec), INTENT(IN)   :: stran(ntens), dstran(ntens)
        REAL(prec), INTENT(IN)   :: props(nprops), dtime, time(2)
        REAL(prec), INTENT(INOUT) :: statev(nstatv)
        REAL(prec), INTENT(OUT)   :: stress(ntens),ddsdde(ntens,ntens)

        !List of internal variables:
        REAL(prec) :: xioi(6, 6), xii(6, 6), xpp(6, 6)
        REAL(prec) :: EPS(6), dev(6), treps
        REAL(prec) :: ka, mu

        REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
        REAL(prec), PARAMETER :: THREE=3.D0
        
        !Define 4th order identity tensor
        data xioi(1,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
        data xioi(2,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
        data xioi(3,:) /ONE, ONE, ONE, ZERO, ZERO, ZERO/
        data xioi(4,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
        data xioi(5,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
        data xioi(6,:) /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/

        !Define 4th order symmetric identity tensor
        data xii(1, :) /ONE, ZERO, ZERO, ZERO, ZERO, ZERO/
        data xii(2, :) /ZERO, ONE, ZERO, ZERO, ZERO, ZERO/
        data xii(3, :) /ZERO, ZERO, ONE, ZERO, ZERO, ZERO/
        data xii(4, :) /ZERO, ZERO, ZERO, 0.5D0, ZERO, ZERO/
        data xii(5, :) /ZERO, ZERO, ZERO, ZERO, 0.5D0, ZERO/
        data xii(6, :) /ZERO, ZERO, ZERO, ZERO, ZERO, 0.5D0/

        !Compute deviatoric projection tensor
        xpp(:, :) = xii(:, :) - (ONE / THREE) * xioi(:, :)

        !Get Material Properties
        ka = PROPS(1)
        mu = PROPS(2)
        
        !DEFINE STATE VARIABLES      
        EPS(:) = DSTRAN(:) + STRAN(:)

        !Trace of strain at time step n+1    
        treps = EPS(1) + EPS(2) + EPS(3)
      
        !Deviator of strain at time step n+1
        dev(1 : 3) = EPS(1 : 3) - treps / THREE
        dev(4 : 6) = EPS(4 : 6) / TWO
        
        !Compute Stresses
        STRESS(1 : 3) = ka * treps + TWO * mu * dev(1 : 3)
        STRESS(4 : 6) = TWO * dev(4 : 6)

        !Compute Tangent Modulus
        DDSDDE(:, :) = ka * xioi(:, :) + TWO * mu * xpp(:, :) 

        RETURN
      END SUBROUTINE LIN_ELAST

