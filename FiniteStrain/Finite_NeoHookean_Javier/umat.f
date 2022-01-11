C     UMAT - Isotropic Neohookian
C     Mohib Mustafa - IMDEA 18 MAY 2021
            
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
      
      call NEOHOOK(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,PROPS,
     &         NPROPS,DTIME,DSTRAN,KINC,KSTEP,NOEL,DFGRD0,DFGRD1)

      return
      end

C     !--------------------------------------------------------------
C     !   UMAT SUBROUTINE FOR ISOTROPIC NEO HOOKIAN MATERIAL MODEL
C     !--------------------------------------------------------------

      SUBROUTINE NEOHOOK(STRESS,STATEV,DDSDDE,STRAN,NTENS,NSTATV,
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
        INTEGER    :: index, ii, jj, mm, K1, K2, ll
        REAL(prec) :: I(6), m2v(3, 3), B(6), B_bar(6), dev_B_bar(6)
        REAL(prec) :: F_inv(3, 3)
        REAL(prec) :: J, I_bar, fac1, fac2, fac3, fac4, cE, cnu
        
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
                   
        !Get material properties
        cE=PROPS(1)
        cnu=PROPS(2)
        
      !   WRITE(*,*) 'cE, cnu : ', cE, ' : ', cnu

        !Calculate determinant of deformation gradient
        J = determinant(DFGRD1(:,:))

        !Calculate inverse of deformation gradient
        CALL matInv(DFGRD1(:,:), F_inv(:, :))

        WRITE(*,*) 'F : '
        DO K1 = 1, 3
          WRITE(*,*) DFGRD1(K1, :)
        END DO

        WRITE(*,*) 'F_inv : '
        DO K1 = 1, 3
          WRITE(*,*) F_inv(K1, :)
        END DO

        WRITE(*,*) 'J : ', J
        WRITE(*,*) ''
        !Calculate finger tensor B
        B(:) = ZERO
        DO K1 = 1, 3
          DO ii = 1, 3
            DO jj = 1, 3
              IF (ii .EQ. jj) THEN
                B(m2v(ii, jj)) = B(m2v(ii, jj)) +  DFGRD1(ii, K1) 
     1                         * DFGRD1(jj, K1)
              ELSE
                B(m2v(ii, jj)) = B(m2v(ii, jj)) + 0.5d0 
     1                         * DFGRD1(ii, K1) * DFGRD1(jj, K1)
              END IF
            END DO
         END DO
        END DO

        B_bar(:) = B(:) *  J ** (-TWO / THREE)
        I_bar = B_bar(1) + B_bar(2) + B_bar(3)
        dev_B_bar(:) = B_bar(:) - (ONE / THREE) * I_bar * I(:)

        !Return Cauchy stress to Abaqus
        STRESS(:) = 2 * cE * (J - ONE) * J * I(:)
     1          + 2 * cnu * dev_B_bar(:)

        !Algorithmic constants for the Tangent
        fac1 = -(FOUR / THREE) * cnu * J ** (-ONE)
        fac2 = -(TWO * cE * (J - ONE) * J 
     1       - (TWO/THREE) * cnu * I_bar) * J ** (-ONE)
        fac3 =  (TWO * cE * (J - ONE) * J 
     1       + (FOUR / NINE) * cnu * I_bar + TWO * cE * J ** TWO) 
     2       * J ** (-ONE)
        fac4 = TWO * cnu * J ** (-ONE)

        !Tangent for ABAQUS
        DDSDDE(:,:) = ZERO
        DO ii = 1, 3
          DO jj = 1, 3
            DO ll = 1, 3
              DO mm = 1, 3
                DDSDDE(m2v(ii, jj), m2v(ll, mm)) = 
     1            DDSDDE(m2v(ii, jj), m2v(ll, mm))
     2            + fac1 * (B_bar(m2v(ii, jj)) * I(m2v(ll, mm)) 
     3            + I(m2v(ii, jj)) * B_bar(m2v(ll, mm)))
     4            + (fac4 - TWO * cnu * J**(-ONE)) 
     5            * B_bar(m2v(ii, ll)) * I(m2v(jj, mm))
     6            + TWO * cnu * (dev_B_bar(m2v(ii, ll)) 
     7            * I(m2v(jj, mm)) 
     8            + I(m2v(ii, ll)) * dev_B_bar(m2v(jj, mm)))
     9            + fac2 * I(m2v(ii, mm)) * I(m2v(ll, jj))
     1            + fac3 * I(m2v(ii, jj)) * I(m2v(ll, mm))
     2            + (fac2 + FOUR * cE * (J - ONE)) 
     3            * I(m2v(ii, ll)) * I(m2v(jj, mm))
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
      END SUBROUTINE NEOHOOK

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
