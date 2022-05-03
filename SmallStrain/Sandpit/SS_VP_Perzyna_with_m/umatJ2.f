C==========================================================================
C                             SUBROUTINE UMAT
C==========================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4           CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c

       INCLUDE 'ABA_PARAM.INC' 
 

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),
     4     DFGRD1(3,3)

     
      call J2(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,NPROPS,
     &         DTIME,DSTRAN,KINC,KSTEPNOEL)

      return
      end                       

C==========================================================================
C                           SUBROUTINE J2
C==========================================================================

      SUBROUTINE J2(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,NPROPS,
     &         DTIME,DSTRAN,KINC,KSTEP,NOEL)

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: double=kind(1.d0)
      INTEGER, PARAMETER :: prec=double

      INTEGER, INTENT(IN)                  :: ntens,nprops,nstatv,kinc,kstep,noel
      REAL(prec), INTENT(IN)               :: stran(ntens), dstran(ntens)
      REAL(prec), INTENT(IN)               :: props(nprops), dtime
      REAL(prec), INTENT(INOUT)            :: statev(nstatv)
      REAL(prec), INTENT(OUT)              :: stress(ntens),ddsdde(ntens,ntens)

      !local variables
      INTEGER    :: i,j
      REAL(prec) :: shear,bulk,sY,Kbar,Dt
      REAL(prec) :: ii(6),iii(6),Id(6,6),Is(6,6),IdT(6,6),Pdev(6,6),Cel(6,6)
      REAL(prec) :: eps(6),epsn(6),epsPn(6),sigman(6),alphan
      REAL(prec) :: epsP(6),alpha
      REAL(prec) :: sigma_trial(6),s_trial(6),p_trial,q_trial,n_trial(6),phi_trial
      REAL(prec) :: dgamma,theta,thetabar

      
      !define material constants
      shear=props(1)/(2.d0*(1.d0+props(2)))
      bulk=props(1)/(3.d0*(1.d0-2.d0*props(2)))
      sY=props(3)
      Kbar=props(4)
      Dt=dtime

      !define useful matrices
C----------------------------------------------------------------------------      

      !vector to get trace of a tensor
      ii(:)=0.d0
      ii(1)=1.d0
      ii(2)=1.d0
      ii(3)=1.d0
      
      !identity vector
      iii(:)=1.d0
      
      !identity matrix (stress)
      Id(:,:) = 0.d0
      do i=1,6
        Id(i,i) = 1.d0
      end do
       
      !identity matrix (strain)
      Is(:,:)=0.d0  
      Is(1,1)=1.d0
      Is(2,2)=1.d0
      Is(3,3)=1.d0
      Is(4,4)=0.5d0
      Is(5,5)=0.5d0
      Is(6,6)=0.5d0
      
      !to get volumetric part of a tensor
      IdT(:,:)=0.d0
      do i=1,3
        do j=1,3
          IdT(i,j)=1.d0
        end do
      end do
      
      !to get deviatoric part of a tensor
      Pdev =Id - 1.d0/3.d0*IdT
            
      !define matrix of elastic constants
      do i=1,6
        do j=1,6
         Cel(i,j)=(bulk-2.d0/3.d0*shear)*IdT(i,j)+2.d0*shear*Is(i,j)
        end do
      end do
      
C----------------------------------------------------------------------------        
      
      !introduce state variables 
      eps(1:6) = stran(1:6)+dstran(1:6)
      epsn(1:6) = stran(1:6)
      epsPn(1:6) = statev(1:6)
      sigman(1:6) = statev(7:12)
      alphan=statev(13)
      
      !initialize variables
      epsP = epsPn
      alpha = alphan
      
C----------------------------------------------------------------------------  

      !trial state
      sigma_trial=sigman
      do i=1,6
        do j=1,6
        sigma_trial(i)=sigma_trial(i)+Cel(i,j)*(eps(j)-epsn(j))
        end do
      end do
      
      s_trial(:)=0.d0
      do i=1,6
        do j=1,6
        s_trial(i)=s_trial(i)+Pdev(i,j)*sigma_trial(j)
        end do
      end do
      
      p_trial=0.d0
      do i=1,6
        p_trial=p_trial-1.d0/3.d0*sigma_trial(i)*ii(i)
      end do
      
      
      q_trial=sqrt(s_trial(1)*s_trial(1)+
     &                        s_trial(2)*s_trial(2)+
     &                        s_trial(3)*s_trial(3)+
     &                   2.d0*s_trial(4)*s_trial(4)+
     &                   2.d0*s_trial(5)*s_trial(5)+
     &                   2.d0*s_trial(6)*s_trial(6))
      
      n_trial=1.d0/q_trial*s_trial

      n_trialm=0.d0
      do i=1,6
        do j=1,6
          n_trialm(i)=n_trialm(i)+Is(i,j)*n_trial(j)
        end do
      end do

      phi_trial = q_trial - sqrt(2.d0/3.d0)*(sY+Kbar*alpha)
      
C---------------------------------------------------------------------------- 
      
      !elastic trial
      
      if (phi_trial.lt.0.d0) then
      
        stress=sigma_trial
        ddsdde=Cel
      
        statev(1:6)=epsP(1:6)
        statev(7:12)=stress(1:6)
        statev(13)=alpha
      
      else
      
C----------------------------------------------------------------------------  

      !plastic corrector
      Dgamma=(q_trial-sqrt(2.d0/3.d0)*(sY+Kbar*alphan))/(2.d0*shear+2.d0/3.d0*Kbar)
      
      !consistent tangent operator
      theta=1.d0-2.d0*shear*Dgamma/q_trial
      thetabar=1.d0/(1.d0+Kbar/3.d0/shear)-(1.d0-theta)
      ddsdde=bulk*IdT+2.d0*shear*theta*(Is-1.d0/3.d0*IdT)
      do i=1,6
        do j=1,6
          ddsdde(i,j)=ddsdde(i,j)-2.d0*shear*thetabar*n_trial(i)*n_trial(j)
        end do
      end do
   

      !update stress and internal variables
      stress=sigma_trial-2.d0*shear*Dgamma*n_trial
      epsP=epsPn+Dgamma*n_trial
      alpha=alphan+sqrt(2.d0/3.d0)*Dgamma
      
      statev(1:6)=epsP(1:6)
      statev(7:12)=stress(1:6)
      statev(13)=alpha

      end if
      
      END SUBROUTINE J2
      
C==========================================================================
C                           SUBROUTINE NORM
C==========================================================================

      SUBROUTINE norm(dnorm,MAT,m,n)
      REAL*8 MAT(m,n),dnorm
      integer ii,jj,m,n

      dnorm=0d0
      DO,ii=1,m
         do,jj=1,n
            dnorm=dnorm+MAT(ii,jj)*MAT(ii,jj)
         enddo
      enddo
      dnorm=sqrt(dnorm)
      RETURN 
      END SUBROUTINE norm
 
C==========================================================================
C                           SUBROUTINE LUDCMP
C==========================================================================
      
      Subroutine LUDCMP(A,N,INDX,D,CODE)
      PARAMETER(NMAX=100,TINY=1d-16)
      REAL*8  AMAX,DUM, SUM, A(N,N),VV(NMAX)
      INTEGER CODE, D, INDX(N)
      
      D=1; CODE=0

      DO I=1,N
         AMAX=0.d0
         DO J=1,N
            IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
         END DO                 ! j loop
         IF(AMAX.LT.TINY) THEN
            CODE = 1
            RETURN
         END IF
         VV(I) = 1.d0 / AMAX
      END DO                    ! i loop
      
      DO J=1,N
         DO I=1,J-1
            SUM = A(I,J)
            DO K=1,I-1
               SUM = SUM - A(I,K)*A(K,J) 
            END DO              ! k loop
            A(I,J) = SUM
         END DO                 ! i loop
         AMAX = 0.d0
         DO I=J,N
            SUM = A(I,J)
            DO K=1,J-1
               SUM = SUM - A(I,K)*A(K,J) 
            END DO              ! k loop
            A(I,J) = SUM
            DUM = VV(I)*DABS(SUM)
            IF(DUM.GE.AMAX) THEN
               IMAX = I
               AMAX = DUM
            END IF
         END DO                 ! i loop  
         
         IF(J.NE.IMAX) THEN
            DO K=1,N
               DUM = A(IMAX,K)
               A(IMAX,K) = A(J,K)
               A(J,K) = DUM
            END DO              ! k loop
            D  = -D
            VV(IMAX) = VV(J)
         END IF
         
         INDX(J) = IMAX
         IF(DABS(A(J,J)) < TINY) A(J,J) = TINY
         
         IF(J.NE.N) THEN
            DUM = 1.d0 / A(J,J)
            DO I=J+1,N
               A(I,J) = A(I,J)*DUM
            END DO              ! i loop
         END IF 
      END DO                    ! j loop
      
      RETURN
      END subroutine LUDCMP

C==========================================================================
C                           SUBROUTINE GAUSS_3
C==========================================================================

      SUBROUTINE gauss_3(AA,BB,xx,nmax,n,iflag)
      integer iflag
      real*8 AA(nmax,nmax),bb(nmax),xx(nmax)
      real*8 A(n,n),b(n),d
      integer n,nmax,i,j,indx(n)
      det=0.
     
      do,i=1,n
         b(i)=bb(i)
         do,j=1,n
            A(i,j)=AA(i,j)
         enddo
      enddo
      call ludcmp(a,n,indx,d,iflag)
      call LUBKSB(A,N,INDX,B)
      do,i=1,n
         xx(i)=b(i)
      enddo
      return
      end
      
C==========================================================================
C                           SUBROUTINE LUBKSB
C==========================================================================

      Subroutine LUBKSB(A,N,INDX,B)
      REAL*8  SUM, A(N,N),B(N)
      INTEGER INDX(N)
      
      II = 0     

      DO I=1,N
         LL = INDX(I)
         SUM = B(LL)
         B(LL) = B(I)
         IF(II.NE.0) THEN
            DO J=II,I-1
               SUM = SUM - A(I,J)*B(J)
            END DO              ! j loop
         ELSE IF(SUM.NE.0.d0) THEN
            II = I
         END IF
         B(I) = SUM
      END DO                    ! i loop
      
      DO I=N,1,-1
         SUM = B(I)
         IF(I < N) THEN
            DO J=I+1,N
               SUM = SUM - A(I,J)*B(J)
            END DO              ! j loop
         END IF
         B(I) = SUM / A(I,I)
      END DO                    ! i loop
          

      RETURN
      END subroutine LUBKSB






	      

      

