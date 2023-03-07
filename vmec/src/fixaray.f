      subroutine fixaray(ierflag)
!
C.. Implicits ..
      implicit none
!
      include 'name2.inc'
!
      integer :: ierflag
!
C.. Local Scalars ..
      integer :: i, m, j, n
      real :: argi, arg, argj, power, t1, rtest, ztest
!
C.. External Functions ..
      real ssum
      external  ssum
!
!
C.. Intrinsic Functions ..
      intrinsic cos, sin, atan
!

!************
!                 INDEX OF LOCAL VARIABLES
!************
!         mscale  array for norming theta-trig functions (internal use only)
!         nscale  array for norming zeta -trig functions (internal use only)
!************
!                 COMPUTE TRIGONOMETRIC FUNCTION ARRAYS
!************
        twopi = c8p0*atan(c1p0)
        dnorm = c1p0/real(nzeta*(ntheta2-1))

      do 15 i = 1,ntheta2
       argi = twopi*real(i-1)/real(ntheta1)
        do 10 m = 0,mpol1
          arg = argi*real(m)
          cosmui(i,m) = cos(arg)*mscale(m)
          if(i.eq.1.or.i.eq.ntheta2)cosmui(i,m) = cp5*cosmui(i,m)
          sinmu (i,m) = sin(arg)*mscale(m)
          cosmu (i,m) = cos(arg)*mscale(m)
          sinmum(i,m) =-sinmu(i,m)*real(m)
          cosmum(i,m) = cosmu(i,m)*real(m)
          cosmumi(i,m)= cosmui(i,m)*real(m)
 10     end do
 15   end do

      do 25 j = 1,nzeta
       argj = twopi*real(j-1)/real(nzeta)
        do 20 n = 0,nmax
          arg = argj*real(n)
          cosnv (j,n) = cos(arg)*nscale(n)
          sinnv (j,n) = sin(arg)*nscale(n)
          cosnvn(j,n) = cosnv(j,n)*real(n*nfp)
          sinnvn(j,n) =-sinnv(j,n)*real(n*nfp)
 20     end do
 25   end do

        do 35 m = 0,mpol1
         power = -cp5*real(m)
         xmpq(m,1) = real(max0(0,m**2-1))
         xmpq(m,2) = real(m**4)
         xmpq(m,3) = real(m**5)
         t1 = -cp075*real(isigng)*dnorm/real((1+m)**2*max0(m,1)**2)
          do 30 n = 0,nmax
            faccon(n,m) = t1/(mscale(m)*nscale(n))**2
            if(m.eq.0.or.n.eq.0)faccon(n,m) = cp5*faccon(n,m)
            xrz3(n,m) = c2p0**(power + c1p0)
            xrz4(n,m) =-c3p0**power
 30     end do
 35   end do
!************
!                 CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS ISIGNG)
!************
        rtest = ssum(nmax1,rb(0,1,1),1)
        ztest = ssum(nmax1,zb(0,1,2),1)
        if( (rtest*ztest*real(isigng)) .ge. czero)ierflag = 5
        return
        end subroutine fixaray
