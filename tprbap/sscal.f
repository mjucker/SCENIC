      SUBROUTINE  SSCAL(N,SA,SX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      REAL SA,SX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF(N.LE.0)RETURN
!      IF(INCX.EQ.1)GO TO 20
      IF(INCX.ne.1)then
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        SX(I) = SA*SX(I)
   10 end do
      RETURN
      end if
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
!      IF( M .EQ. 0 ) GO TO 40
      IF( M .ne. 0 ) then
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 end do
      IF( N .LT. 5 ) RETURN
      end if
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 end do
      RETURN
      END
