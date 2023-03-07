!      MODULE FUNCTIONS
!      CONTAINS
      function isamax (n,sx,incx)
c  ******************************************************************
c
c 1. function
c  isamax returns the index of element absolute mumimum value 
c 2.  description
c  n     ; in   : Number of elements to process in the vector to be searched
c                  (n=vector length if incx=1;n=vector length/2 if incx=2;
c                    and so on).
c                 if n<=0,isamax returns 0 .
c  sx    ; in   : Real vector to be searched
c  incx  ; in   : Increment between elements of sx
c  ******************************************************************
c
      real :: sx(*), vmax
      integer :: n, incx, ind, j,  i, isamax
      intrinsic abs
c
      ind=0
!      if(n.le.0) goto 20
      if(n.gt.0) then
        if(incx.ge.0) then
          j=1
        else
          j=1-incx*(n-1)
        endif
        vmax=abs(sx(j))
        ind=1
        do 10 i=2,n
          j=j+incx
!        if(vmax.ge.abs(sx(j))) goto 10
        if(vmax.lt.abs(sx(j))) then
          vmax=abs(sx(j))
          ind=i
        end if
   10   end do
        end if
   20   isamax=ind
      return
      end function isamax
!
!
      real function ssum (dim,ary,dummy)                                    1
C...Translated by  fopp     4.02E36 17:59:18  12/14/93
C...Switches: -eadejlpuvx18 -dchimr07 -e1 -gi
      integer   k,dummy,dim                                                 2
      real, dimension(dim) ::  ary                                          3
c     ssum = sum(ary,dim)                                                   4
      ssum=0.                                                               5
*VDIR NODEP
      do k=1,dim,dummy                                                      6
         ssum=ssum+ary(k)                                                   7
      enddo                                                                 8
      return                                                               10
      end function ssum                                                    10
                                                                           11
      real function sdot (dim,v1,d1,v2,d2)
C...Translated by  fopp     4.02E36 17:59:18  12/14/93
C...Switches: -eadejlpuvx18 -dchimr07 -e1 -gi
      integer j1
      integer   dim,d1,d2 
      real :: v1(dim),v2(dim)
      real*8 d3
      d3 = 0
*VDIR NODEP
      do j1 = 1, dim 
         d3 = d3 + v1(1+(j1-1)*d1)*v2(1+(j1-1)*d2)
      end do
      sdot = d3
      return
      end function sdot
 
      subroutine scopy(nn,x,incx,y,incy)
      real :: x(1),y(1)
      integer nn, incx, incy, i
      do i = 1, nn
      y(1+(i-1)*incy) = x(1+(i-1)*incx)
      end do
      return
      end subroutine scopy

      SUBROUTINE SAXPY(n,sa,sx,incx,sy,incy)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      REAL, DIMENSION (1) ::  sx, sy
      REAL :: sa
      INTEGER :: I,incx,incy,IX,IY,M,MP1,n
C
      IF(n.LE.0)RETURN
      IF (sa .EQ. 0.0) RETURN
      IF(incx.NE.1.OR.incy.NE.1)THEN
C   
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(incx.LT.0)IX = (-N+1)*incx + 1
      IF(incy.LT.0)IY = (-N+1)*incy + 1
      DO 10 I = 1,n
        sy(IY) = sy(IY) + sa*sx(IX)
        IX = IX + incx
        IY = IY + incy
   10 END DO
      RETURN
      END if
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(n,4)
      IF( M .NE. 0 )THEN
      DO 30 I = 1,M
        sy(I) = sy(I) + sa*sx(I)
   30 END DO
      IF( n .LT. 4 ) RETURN
      END if
   40 MP1 = M + 1
      DO 50 I = MP1,n,4
        sy(I) = sy(I) + sa*sx(I)
        sy(I + 1) = sy(I + 1) + sa*sx(I + 1)
        sy(I + 2) = sy(I + 2) + sa*sx(I + 2)
        sy(I + 3) = sy(I + 3) + sa*sx(I + 3)
   50 END DO
      RETURN
      END SUBROUTINE SAXPY

      SUBROUTINE SSCAL(n,sa,sx,incx)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      real :: sa
      real, dimension (1) :: sx
      integer :: I,incx,M,MP1,Nincx,n
C
      IF(N.LE.0)return
!      IF(incx.EQ.1)GO TO 20
      IF(incx.ne.1)then
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      Nincx = n*incx
      DO 10 I = 1,Nincx,incx
        sx(I) = sa*sx(I)
   10 END DO
      RETURN
      end if
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(n,5)
!      IF( M .EQ. 0 ) GO TO 40
      IF( M .ne. 0 ) then
      DO 30 I = 1,M
        sx(I) = sa*sx(I)
   30 end do
      IF( n .LT. 5 ) return
      end if
   40 MP1 = M + 1
      DO 50 I = MP1,n,5
        sx(I) = sa*sx(I)
        sx(I + 1) = sa*sx(I + 1)
        sx(I + 2) = sa*sx(I + 2)
        sx(I + 3) = sa*sx(I + 3)
        sx(I + 4) = sa*sx(I + 4)
   50 END DO
      RETURN
      END SUBROUTINE SSCAL
      SUBROUTINE SSWAP (n,sx,incx,sy,incy)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      real :: sx(1),sy(1),STEMP
      integer :: I,incx,incy,IX,IY,M,MP1,n
C
      IF(n.LE.0)return
!      IF(incx.EQ.1.AND.incy.EQ.1)GO TO 20
      IF(incx.ne.1.or.incy.ne.1)then
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(incx.LT.0)IX = (-N+1)*incx + 1
      IF(incy.LT.0)IY = (-N+1)*incy + 1
      DO 10 I = 1,n
        STEMP = sx(IX)
        sx(IX) = sy(IY)
        sy(IY) = STEMP
        IX = IX + incx
        IY = IY + incy
   10 end do
      RETURN
      end if
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(n,3)
!      IF( M .EQ. 0 ) GO TO 40
      IF( M .ne. 0 ) then
      DO 30 I = 1,M
        STEMP = sx(I)
        sx(I) = sy(I)
        sy(I) = STEMP
   30 end do
      IF( n .LT. 3 ) return
      end if
   40 MP1 = M + 1
      DO 50 I = MP1,n,3
        STEMP = sx(I)
        sx(I) = sy(I)
        sy(I) = STEMP
        STEMP = sx(I + 1)
        sx(I + 1) = sy(I + 1)
        sy(I + 1) = STEMP
        STEMP = sx(I + 2)
        sx(I + 2) = sy(I + 2)
        sy(I + 2) = STEMP
   50 end do
      RETURN
      END SUBROUTINE SSWAP
      real function cvmgp (TR,FR,MR)
      real :: TR, FR, MR
      cvmgp = merge(TR,FR,MR>=0.0)
      return
      end function cvmgp
      real function cvmgz (TR,FR,MR)
      real :: TR, FR, MR
      cvmgz = merge(TR,FR,MR==0.0)
      return
      end function cvmgz
C****************************************************************************A
      subroutine second(tt)
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
C.. In/Out Status: Not Read, Overwritten ..
      real tt
C
C.. Local Scalars ..
!      integer i
!      real tt1
C
C.. Local Arrays ..
      real ra(2)
C
C.. External Functions ..
!      real etime
!      external etime
C
C ... Executable Statements ...
C
!      do i=1,2
!         ra(i) = 0.0
!      end do
C      call etime(ra)
      call cpu_time(tt)
!      tt1 = etime(ra)
!      tt = ra(1)
      end subroutine second
C****************************************************************************A
      subroutine fdate(t)
C
C.. Implicits ..
      implicit none
C
C.. External Functions..
      external date
      character t*24
      call date(t)
      end subroutine fdate
!
!      END MODULE FUNCTIONS
