        subroutine totzsp(rmncc,rmnss,zmncs,zmnsc,lmncs,lmnsc,
     &  r1,ru,rv,z1,zu,zv,lu,lv,rcon,zcon,
     &  work1,work2,work3,realsp)
!--------0---------0---------0---------0---------0---------0---------0-c
!
C.. Implicits ..
      implicit none
!
      include 'name2.inc'
       real, dimension(ns,0:nmax,0:mpol1)    :: rmncc,rmnss,zmncs,zmnsc
       real, dimension(ns,0:nmax,0:mpol1)    :: lmncs,lmnsc
       real, dimension(ns*nzeta,ntheta2,0:1) :: r1,ru,rv,z1,zu,zv,lu,lv
       real, dimension(ns*nzeta,ntheta2,0:1) :: rcon,zcon
       real, dimension(ns*nzeta*12)          :: work1
       real, dimension(ns*nzeta,12)          :: work2
       real, dimension(ns,nzeta,12)          :: work3
       real, dimension(16*nrztd)             :: realsp
!
C.. Local Scalars ..
      integer :: l, n, m, k, js, i, jk, mparity
      real :: cosmux, sinmux
!************
!                 THIS ROUTINE ASSUMES THE FOLLOWING STACKING OF R, Z, 
!                 LAMBDA ARRAYS:
!                 rmncc(ns,0:nmax,0:mpol1),rmnss,zmncs,zmncc,lmncs,lmnsc
!************
!                 INITIALIZATION BLOCK
!************
      do 10 l = 1,16*nrztd
         realsp(l) = czero
   10 end do
!************
!                 EXTRAPOLATION AT JS=1 FOR M=1 MODES
!************
      do 30 n = 0,nmax
        rmncc(1,n,1) = c2p0*rmncc(2,n,1) - rmncc(3,n,1)
        rmnss(1,n,1) = c2p0*rmnss(2,n,1) - rmnss(3,n,1)
        zmncs(1,n,1) = c2p0*zmncs(2,n,1) - zmncs(3,n,1)
        zmnsc(1,n,1) = c2p0*zmnsc(2,n,1) - zmnsc(3,n,1)
        lmncs(1,n,1) = c2p0*lmncs(2,n,1) - lmncs(3,n,1)
        lmnsc(1,n,1) = c2p0*lmnsc(2,n,1) - lmnsc(3,n,1)
   30 end do
!************
!                 COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!                 BEGIN INVERSE TRANSFORM IN N-ZETA
!************
      do 70 m = 0,mpol1
       mparity = mod(m,2)
        do 35 l = 1,12*ns*nzeta
           work1(l) = czero
   35   end do
       do 50 n = 0,nmax
        do 45 k = 1,nzeta
!*VDIR NODEP
CDIR$ IVDEP
         do 40 js= jmin1(m),ns
          work3(js,k,1) = work3(js,k,1) + rmncc(js,n,m)*cosnv (k,n)
          work3(js,k,2) = work3(js,k,2) + rmnss(js,n,m)*sinnv (k,n)
          work3(js,k,3) = work3(js,k,3) + rmncc(js,n,m)*sinnvn(k,n)
          work3(js,k,4) = work3(js,k,4) + rmnss(js,n,m)*cosnvn(k,n)
          work3(js,k,5) = work3(js,k,5) + zmncs(js,n,m)*sinnv (k,n)
          work3(js,k,6) = work3(js,k,6) + zmnsc(js,n,m)*cosnv (k,n)
          work3(js,k,7) = work3(js,k,7) + zmncs(js,n,m)*cosnvn(k,n)
          work3(js,k,8) = work3(js,k,8) + zmnsc(js,n,m)*sinnvn(k,n)
          work3(js,k,9) = work3(js,k,9) + lmncs(js,n,m)*sinnv (k,n)
          work3(js,k,10)= work3(js,k,10)+ lmnsc(js,n,m)*cosnv (k,n)
          work3(js,k,11)= work3(js,k,11)+ lmncs(js,n,m)*cosnvn(k,n)
          work3(js,k,12)= work3(js,k,12)+ lmnsc(js,n,m)*sinnvn(k,n)
   40    end do
   45   end do
   50  end do
!************
!                 INVERSE TRANSFORM IN M-THETA
!************
        do 65 i = 1,ntheta2
          cosmux = xmpq(m,1)*cosmu(i,m)
          sinmux = xmpq(m,1)*sinmu(i,m)
        do 60 jk = 1, nzeta*ns
          r1(jk,i,mparity) = r1(jk,i,mparity) +
     &        work2(jk,01)*cosmu (i,m) + work2(jk,02)*sinmu (i,m)
          ru(jk,i,mparity) = ru(jk,i,mparity) +
     &        work2(jk,02)*cosmum(i,m) + work2(jk,01)*sinmum(i,m)
          rv(jk,i,mparity) = rv(jk,i,mparity) +
     &        work2(jk,03)*cosmu (i,m) + work2(jk,04)*sinmu (i,m)
          z1(jk,i,mparity) = z1(jk,i,mparity) +
     &        work2(jk,05)*cosmu (i,m) + work2(jk,06)*sinmu (i,m)
          zu(jk,i,mparity) = zu(jk,i,mparity) +
     &        work2(jk,06)*cosmum(i,m) + work2(jk,05)*sinmum(i,m)
          zv(jk,i,mparity) = zv(jk,i,mparity) +
     &        work2(jk,07)*cosmu (i,m) + work2(jk,08)*sinmu (i,m)
          lu(jk,i,mparity) = lu(jk,i,mparity) +
     &        work2(jk,10)*cosmum(i,m) + work2(jk,09)*sinmum(i,m)
          lv(jk,i,mparity) = lv(jk,i,mparity) -
     &       (work2(jk,11)*cosmu (i,m) + work2(jk,12)*sinmu (i,m))
          rcon(jk,i,mparity) = rcon(jk,i,mparity) +
     &        work2(jk,01)*cosmux      + work2(jk,02)*sinmux
          zcon(jk,i,mparity) = zcon(jk,i,mparity) +
     &        work2(jk,05)*cosmux      + work2(jk,06)*sinmux
   60   end do
   65   end do
   70 end do
        return
        end subroutine totzsp
