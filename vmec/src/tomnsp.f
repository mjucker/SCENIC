        subroutine tomnsp(frcc,frss,fzcs,fzsc,flcs,flsc,armn,brmn,
     &  crmn,azmn,bzmn,czmn,blmn,clmn,arcon,azcon,work1,work2,work3)
!--------0---------0---------0---------0---------0---------0---------0-c
!
C.. Implicits ..
      implicit none
!
      include 'name2.inc'
      real, dimension(ns,0:nmax,0:mpol1) :: frcc,frss,fzcs,fzsc,flcs
      real, dimension(ns,0:nmax,0:mpol1) :: flsc
      real, dimension(ns*nzeta,ntheta2,0:1) :: armn,brmn,crmn,azmn,bzmn
      real, dimension(ns*nzeta,ntheta2,0:1) :: czmn,blmn,clmn,arcon
      real, dimension(ns*nzeta,ntheta2,0:1) :: azcon
      real, dimension(ns,nzeta,12)          :: work3
      real, dimension(ns*nzeta,12)          :: work2
      real, dimension(ns*nzeta*12)          :: work1
!
C.. Local Scalars ..
      integer :: jmax, m, mparity, l, i, jk, n, k, js
      real    :: temp1, temp3
!
        jmax = ns
        if( ivac.lt.1 )jmax = ns-1
!************
!                 BEGIN INVERSE FOURIER TRANSFORM
!                 DO THETA (U) INTEGRATION FIRST
!************
      do 60 m = 0,mpol1
        mparity = mod(m,2)
        do 10 l = 1,12*ns*nzeta
          work1(l) = czero
   10   end do
        do 30 i = 1,ntheta2
         do 20 jk= 1,ns*nzeta
          temp1 = armn(jk,i,mparity) + xmpq(m,1)*arcon(jk,i,mparity)
          temp3 = azmn(jk,i,mparity) + xmpq(m,1)*azcon(jk,i,mparity)
          work2(jk,01) = work2(jk,01) + temp1               *cosmui(i,m)
     &                                + brmn(jk,i,mparity)*sinmum(i,m)
          work2(jk,02) = work2(jk,02) - crmn(jk,i,mparity)*cosmui(i,m)
          work2(jk,03) = work2(jk,03) + temp1               *sinmu (i,m)
     &                                + brmn(jk,i,mparity)*cosmumi(i,m)
          work2(jk,04) = work2(jk,04) - crmn(jk,i,mparity)*sinmu (i,m)
          work2(jk,05) = work2(jk,05) + temp3               *cosmui(i,m)
     &                                + bzmn(jk,i,mparity)*sinmum(i,m)
          work2(jk,06) = work2(jk,06) - czmn(jk,i,mparity)*cosmui(i,m)
          work2(jk,07) = work2(jk,07) + temp3               *sinmu (i,m)
     &                                + bzmn(jk,i,mparity)*cosmumi(i,m)
          work2(jk,08) = work2(jk,08) - czmn(jk,i,mparity)*sinmu (i,m)
          work2(jk,09) = work2(jk,09) + blmn(jk,i,mparity)*sinmum(i,m)
          work2(jk,10) = work2(jk,10) - clmn(jk,i,mparity)*cosmui(i,m)
          work2(jk,11) = work2(jk,11) + blmn(jk,i,mparity)*cosmumi(i,m)
          work2(jk,12) = work2(jk,12) - clmn(jk,i,mparity)*sinmu (i,m)
   20    end do
   30   end do
!************
!                 NEXT, DO ZETA (V) INTEGRATION
!************
        do 57 n = 0,nmax
         do 56 k = 1,nzeta
          do 50 js= jmin2(m),jmax
           frcc(js,n,m) = frcc(js,n,m) + work3(js,k,01)*cosnv (k,n)
     &                                 + work3(js,k,02)*sinnvn(k,n)
           frss(js,n,m) = frss(js,n,m) + work3(js,k,03)*sinnv (k,n)
     &                                 + work3(js,k,04)*cosnvn(k,n)
           fzcs(js,n,m) = fzcs(js,n,m) + work3(js,k,05)*sinnv (k,n)
     &                                 + work3(js,k,06)*cosnvn(k,n)
           fzsc(js,n,m) = fzsc(js,n,m) + work3(js,k,07)*cosnv (k,n)
     &                                 + work3(js,k,08)*sinnvn(k,n)
   50     end do
          do 55 js= jlam(m),ns
           flcs(js,n,m) = flcs(js,n,m) + work3(js,k,09)*sinnv (k,n)
     &                                 + work3(js,k,10)*cosnvn(k,n)
           flsc(js,n,m) = flsc(js,n,m) + work3(js,k,11)*cosnv (k,n)
     &                                 + work3(js,k,12)*sinnvn(k,n)
   55     end do
   56    end do
   57   end do
   60 end do
        return
        end subroutine tomnsp
