C     UMAT - Isotropic Hyperelastic from Computational inealsticity
C     Goal is to implement just the initail stored energy formulation
C     Later it will be modified to incorporate Kinematically non linear VE

C     Mohib Mustafa - IMDEA 14 JULY 2021
            
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,
     3           DPRED,CMNAME,
     4           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     5           DROT,PNEWDT,
     6           CELENT,DFGRD0,DFGRD1,NOEL,NPT,
     7           LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),DFGRD0(3,3), DFGRD1(3,3)
      
      call HyperElast(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &         NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1)

      return
      end

C     !--------------------------------------------------------------
C     !   UMAT SUBROUTINE FOR ISOTROPIC NEO HOOKIAN MATERIAL MODEL
C     !--------------------------------------------------------------

      SUBROUTINE HyperElast(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,
     1   PROPS,NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1)


        IMPLICIT NONE

        INTERFACE
          FUNCTION determinant(matrix)
            INTEGER, PARAMETER :: double=kind(1.d0)
            INTEGER, PARAMETER :: prec=double

            REAL(prec), INTENT(IN) :: matrix(3, 3)
            REAL(prec) :: determinant
          END FUNCTION determinant
        END INTERFACE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
       
        INTEGER, INTENT(IN)      :: ntens,nprops,nstatv
        INTEGER, INTENT(IN)      :: kinc,kstep,noel
        REAL(prec), INTENT(IN)   :: stran(ntens), dstran(ntens)
        REAL(prec), INTENT(IN)   :: props(nprops), dtime
        REAL(prec), INTENT(IN)   :: DFGRD0(3,3), DFGRD1(3,3)
        REAL(prec), INTENT(INOUT) :: statev(nstatv)
        REAL(prec), INTENT(OUT)   :: stress(ntens)
        REAL(prec), INTENT(OUT)   :: ddsdde(ntens,ntens)

        !List of internal variables
        INTEGER    :: index, ii, jj, mm, ll, K1
        REAL(prec) :: I(6), m2v(3, 3), B_bar(6), tau(6)
        REAL(prec) :: F_inv(3, 3), F_bar(3, 3), BB(6, 6)
        REAL(prec) :: c_strike(6, 6), xioi(6, 6), xii(6, 6), xpp(6,6)
        REAL(prec) :: c_strike_bar(6, 6)
        REAL(prec) :: J, kappa, mu, dev_B_bar(6), trB_bar, dev_tau(6)
        
        !Decleration of constants
        REAL(prec), PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0
        REAL(prec), PARAMETER :: THREE=3.D0, FOUR= 4.D0, SIX=6.D0
        REAL(prec), PARAMETER :: NINE=9.D0


        !Declare matrix to voit mapping
        data m2v/ 1, 4, 5,
     1            4, 2, 6,
     2            5, 6, 3/

        !2nd Order Identity
        data I/ ONE, ONE, ONE, ZERO, ZERO, ZERO/
        
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

        !Get material properties
        kappa=PROPS(1)
        mu=PROPS(2)
        
      !   WRITE(*,*) 'cE, cnu : ', cE, ' : ', cnu

        !Calculate determinant of deformation gradient
        J = determinant(DFGRD1(:,:))

        !Calculate inverse of deformation gradient
        CALL matInv(DFGRD1(:,:), F_inv(:, :))
        
        F_bar(:, :) = (J ** (-ONE / THREE)) * DFGRD1(:, :)

        WRITE(*,*) 'F : '
        DO K1 = 1, 3
          WRITE(*,*) DFGRD1(K1, :)
        END DO

        WRITE(*,*) 'J : ', J
        WRITE(*,*) ''

        WRITE(*,*) 'F_bar : '
        DO K1 = 1, 3
          WRITE(*,*) F_bar(K1, :)
        END DO
        
       
        !Calculate Deviatoric finger tensor B_bar
        B_bar(:) = ZERO
        DO K1 = 1, 3
          DO ii = 1, 3
            DO jj = 1, 3
              IF (ii .EQ. jj) THEN
                B_bar(m2v(ii, jj)) = B_bar(m2v(ii, jj)) 
     1                             + F_bar(ii, K1) * F_bar(jj, K1)
              ELSE
                B_bar(m2v(ii, jj)) = B_bar(m2v(ii, jj)) + 0.5d0 
     1                             * F_bar(ii, K1) * F_bar(jj, K1)
              END IF
            END DO
         END DO
        END DO

        WRITE(*,*) 'B_bar : ', B_bar(:)

        CALL devVoit(B_bar, dev_B_bar)

        trB_bar = B_bar(1) + B_bar(2) + B_bar(3)

        WRITE(*,*) 'dev_B_bar : ', dev_B_bar(:)
        Write(*,*) 'tr B_bar : ', trB_bar
        write(*,*) '------------------------------------------------'

        !Calculate Kirchoff stress
        tau(:) = 0.5D0 * kappa * (J**TWO - ONE) * I(:) 
     1         + mu * dev_B_bar(:)

        CALL devVoit(tau, dev_tau)
        
        !Return Cauchy stress to Abaqus
        STRESS(:) = (ONE / J) * tau(:)

        !Algorithmic constants for the Tangent
        BB(:, :) = ZERO
        DO ii = 1, 3
          DO jj = 1, 3
            DO ll = 1, 3
              DO mm = 1, 3
        BB(m2v(ii, jj), m2v(ll,mm)) = BB(m2v(ii, jj), m2v(ll,mm))
     1                      + dev_tau(m2v(ii, jj)) * I(m2v(ll, mm))
     2                      + I(m2v(ii, jj)) * dev_tau(m2v(ll, mm))
              END DO
            END DO
          END DO
        END DO
        BB(:, :) = (-TWO/THREE) * BB(:, :)

        c_strike_bar(:, :) = BB(:, :) 
     1                     + (TWO/THREE) * mu * trB_bar * xpp(:,:)



        ! Tangennt in current config
        c_strike(:, :) = 0.5D0 * kappa * (J**TWO + ONE) * xioi(:, :)
     1                 + 0.5D0 * kappa * (J**TWO - ONE) 
     2                 * (xioi(:, :) - TWO * xii(:, :)) 
     3                 + c_strike_bar(:, :)  

        !Tangent for ABAQUS
        DO ii = 1, 3
          DO jj = 1, 3
            DO ll = 1, 3
              DO mm = 1, 3
                DDSDDE(m2v(ii, jj), m2v(ll, mm))   
     1        = DDSDDE(m2v(ii, jj), m2v(ll, mm)) 
     2        + (ONE / J) * (c_strike(m2v(ii,jj), m2v(ll, mm)) 
     3        + I(m2v(ii, ll)) * tau(m2v(mm, jj)) 
     4        + I(m2v(jj, mm)) * tau(m2v(ii, ll)))
              END DO
            END DO
          END DO
        END DO

        index = 0
        DO ii = 1, 3
          DO jj = 1, 3
            index = index + 1
            STATEV(index) = DFGRD1(ii,jj)
          END DO
        END DO

      RETURN
      END SUBROUTINE HyperElast

      !--------------------------------------------------------------
      !     Helper function to compute Determinant of 3x3 matrix 
      !--------------------------------------------------------------
      FUNCTION determinant(matrix)
        IMPLICIT NONE
                
        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec):: determinant

        determinant = matrix(1,1) * (matrix(2,2) * matrix(3,3) 
     1   - matrix(3,2) * matrix(2, 3)) - matrix(1,2) * (matrix(2,1)
     2   * matrix(3,3) - matrix(2,3) * matrix(3,1)) + matrix(1,3) 
     3   * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))

      END FUNCTION determinant


      !--------------------------------------------------------------
      !      Helper function to compute inverse of 3x3 matrix
      !--------------------------------------------------------------
      SUBROUTINE matInv(matrix, inv)
        IMPLICIT NONE
        
        INTERFACE
          FUNCTION determinant(mat)
            INTEGER, PARAMETER :: double=kind(1.d0)
            INTEGER, PARAMETER :: prec=double

            REAL(prec), INTENT(IN) :: mat(3, 3)
            REAL(prec) :: determinant
          END FUNCTION determinant
        END INTERFACE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: matrix(3, 3)
        REAL(prec), INTENT(OUT) :: inv(3, 3)

        inv(1,1) = matrix(2,2)*matrix(3,3) - matrix(3,2)*matrix(2,3)
        inv(1,2) = matrix(3,2)*matrix(1,3) - matrix(1,2)*matrix(3,3)
        inv(1,3) = matrix(1,2)*matrix(2,3) - matrix(2,2)*matrix(1,3)
        inv(2,1) = matrix(3,1)*matrix(2,3) - matrix(2,1)*matrix(3,3)
        inv(2,2) = matrix(1,1)*matrix(3,3) - matrix(3,1)*matrix(1,3)
        inv(2,3) = matrix(2,1)*matrix(1,3) - matrix(1,1)*matrix(2,3)
        inv(3,1) = matrix(2,1)*matrix(3,2) - matrix(3,1)*matrix(2,2)
        inv(3,2) = matrix(3,1)*matrix(1,2) - matrix(1,1)*matrix(3,2)
        inv(3,3) = matrix(1,1)*matrix(2,2) - matrix(2,1)*matrix(1,2)
      
        inv(:, :) = inv(:, :)/determinant(matrix(:, :))

        RETURN
      END SUBROUTINE matInv

      !--------------------------------------------------------------
      !      Helper function to compute dev of voit array
      !--------------------------------------------------------------
      SUBROUTINE devVoit(array, dev)
        IMPLICIT NONE

        INTEGER, PARAMETER :: double=kind(1.d0)
        INTEGER, PARAMETER :: prec=double
        
        REAL(prec), INTENT(IN) :: array(6)
        REAL(prec), INTENT(OUT) :: dev(6)

        dev(1 : 3) = array(1 : 3) 
     1             - (1.D0/3.D0) * (array(1) + array(2) + array(3))

        dev(4 : 6) = array(4 : 6)
        RETURN
      END SUBROUTINE devVoit