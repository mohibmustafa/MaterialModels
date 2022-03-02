      PROGRAM MAIN
      IMPLICIT NONE

      INTEGER, PARAMETER :: double=kind(1.d0)
      INTEGER, PARAMETER :: prec=double
      
      INTEGER    :: KINC, KSTEP
      REAL(prec) :: STRESS(6),STATEV(38)
      REAL(prec) :: DDSDDE(6,6),DDSDDT(6),DRPLDE(6)
      REAL(prec) :: STRAN(6),DSTRAN(6),TIME(2),PREDEF(1),DPRED(1)
      REAL(prec) :: PROPS(19),DFGRD0(3,3),DFGRD1(3,3),DTIME(1),NTENS(1)
      REAL(prec) :: NSTATV(1), NPROPS(1), NOEL(1)

      STRESS(:) = 0.D0
      STATEV(:) = 0.D0
      DDSDDE(:, :) = 0.D0
      DDSDDT(:) = 0.D0
      DRPLDE(:) = 0.D0
      STRAN(:) = 0.D0
      DSTRAN(:) = 0.D0
      TIME(:) = 0.D0
      PREDEF(:) = 0.D0
      DPRED(:) = 0.D0
      NTENS = 0.D0
      NSTATV = 0.D0
      NPROPS = 0.D0
      KINC = 1
      KSTEP = 1
      NOEL = 0.D0


      DTIME = 0.01


      data PROPS/2250.0, 1150.0, 7202.0, 1038.4615384615383,
     1  530.76923076923072, 7202.0, 11, 0.78, 3.5,
     2  0.42105263157894735, 30000.0, 0.21, 100.0, 300.0, 0.0, 10.0,
     3  0.0, 370.0, 0.0/

      data DFGRD0(1,:) /1.D0, 0.D0, 0.D0/
      data DFGRD0(2,:) /0.D0, 1.D0, -0.0052253236793590712D0/
      data DFGRD0(3,:) /0.D0, 0.D0, 1.0183333333333333D0/

      data DFGRD1(1,:) /1.D0, 0.D0 , 0.D0/
      data DFGRD1(2,:) /0.D0, 1.D0, -0.0052253236793590712D0/
      data DFGRD1(3,:) /0.D0, 0.D0, 1.02D0/

      KINC = 2
      KSTEP = 2

      DTIME = 0.01

      STATEV(1) = 0.99999999989585497D0
      STATEV(2) = 0.99999999989567556D0
      STATEV(3) = 1.0000000003850633D0
      STATEV(4) = 0.D0
      STATEV(5) = 0.D0
      STATEV(6) = -6.9057616628254854D-11
      STATEV(7) = 0.D0
      STATEV(8) = 0.D0
      STATEV(9) = -6.9057616628254828D-11

      STATEV(10) = 0.018167164623473293D0

      STATEV(11) = -0.0060557214370127004D0
      STATEV(12) = -0.0060623835261151161D0
      STATEV(13) = 0.012118104963127813D0
      STATEV(14) = 0.D0
      STATEV(15) = 0.D0
      STATEV(16) = -0.0025654529455137139D0
      STATEV(17) = 0.D0
      STATEV(18) = 0.D0
      STATEV(19) = -0.0025654529455137144D0

      STATEV(20) = 1.0414513694793137D-10
      STATEV(21) = -6.6620193775470514D-6
      STATEV(22) = 0.018173965694087944D0
      STATEV(23) = 0.D0
      STATEV(24) = 0.D0
      STATEV(25) = -0.0025654727161421968D0
      STATEV(26) = 0.D0
      STATEV(27) = 0.D0
      STATEV(28) = -0.002565472716142197D0

      STATEV(29) = -3.3108509812844912D-8
      STATEV(30) = -3.3165520744013637D-8
      STATEV(31) = 1.2241454892020397D-7
      STATEV(32) = 0.D0
      STATEV(33) = 0.D0
      STATEV(34) = -2.195394040455492D-8
      STATEV(35) = 0.D0
      STATEV(36) = 0.D0
      STATEV(37) = -2.195394040455492D-8

      STATEV(38) = 3.64067175811993D-10


      
      call VEVP(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &         NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1)

      END PROGRAM MAIN