! PROGRAM TO TEST THE OPTIONAL ARGUMENTS TO RESHAPE
  INTERFACE
     SUBROUTINE WRITE_MATRIX(A)
         REAL, DIMENSION(:,:) :: A
     END SUBROUTINE  WRITE_MATRIX
  END INTERFACE

  REAL, DIMENSION (1:9) :: B = (/ 11, 12, 13, 14, 15, 16, 17, 18, 19 /)
  REAL, DIMENSION (1:9) :: A = (/ 1, 1, 1, 2, 2, 2, 3, 3, 3 /)
  REAL, DIMENSION (1:9) :: Bpack,Apack,Cfill

  REAL, DIMENSION (1:3, 1:3) :: C, D, E
  REAL, DIMENSION (1:4, 1:4) :: F, G, H

  INTEGER, DIMENSION (1:2) :: ORDER1 = (/ 1, 2 /)
  INTEGER, DIMENSION (1:2) :: ORDER2 = (/ 2, 1 /)
  REAL, DIMENSION (1:16)   :: PAD1 = (/ -1, -2, -3, -4, -5, -6, -7, -8, &
                                 &   -9, -10, -11, -12, -13, -14, -15, -16 /)

  Cfill=0
  Bpack = pack(B,A==2,Cfill)

  print*, Bpack
  ! C = RESHAPE( B, (/ 3, 3 /) )
  ! CALL WRITE_MATRIX(C)

  ! D = RESHAPE( B, (/ 3, 3 /), ORDER = ORDER1)
  ! CALL WRITE_MATRIX(D)

  ! E = RESHAPE( B, (/ 3, 3 /), ORDER = ORDER2)
  ! CALL WRITE_MATRIX(E)

  ! F = RESHAPE( B, (/ 4, 4 /), PAD = PAD1)
  ! CALL WRITE_MATRIX(F)

  ! G = RESHAPE( B, (/ 4, 4 /), PAD = PAD1, ORDER = ORDER1)
  ! CALL WRITE_MATRIX(G)

  ! H = RESHAPE( B, (/ 4, 4 /), PAD = PAD1, ORDER = ORDER2)
  ! CALL WRITE_MATRIX(H)

  END

  SUBROUTINE WRITE_MATRIX(A)
  REAL, DIMENSION(:,:) :: A
  WRITE(*,*)
  DO I = LBOUND(A,1), UBOUND(A,1)
     WRITE(*,*) (A(I,J), J = LBOUND(A,2), UBOUND(A,2))
  END DO
  END SUBROUTINE WRITE_MATRIX
