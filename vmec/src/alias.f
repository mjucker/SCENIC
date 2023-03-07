!
        subroutine alias(gcon,zcon,work1,work2,work3,gcs,gsc)
        include 'name2.inc'
        real, dimension(ns*nzeta,ntheta2)  :: gcon,zcon
        real, dimension(ns,0:nmax,0:mpol1) :: gcs,gsc
        real, dimension(ns,nzeta,4)        :: work3
        real, dimension(ns*nzeta,4)        :: work2
        real, dimension(4*ns*nzeta)        :: work1
C.. Local Scalars ..
      integer :: m, l, i, jk, n, k, js
      real :: fm
*************
*                 BEGIN DE-ALIASING (TRUNCATION OF GCON IN FOURIER-SPACE)
*************
        do 60 m = 1,mpol1-1
          do 20 l = 1,4*ns*nzeta
            work1(l) = czero
 20       end do
        do 26 i = 1,ntheta2
          do 25 jk = 1,ns*nzeta
           work2(jk,01) = work2(jk,01) + zcon(jk,i)*cosmui(i,m)
           work2(jk,02) = work2(jk,02) + zcon(jk,i)*sinmu (i,m)
 25       end do
 26     end do
        do 32 n = 0,nmax
         fm = faccon(n,m)
          do 31 k = 1,nzeta
           do 30 js= 2,ns
        gcs(js,n,m) =gcs(js,n,m) +fm*tcon(js)*work3(js,k,01)*sinnv(k,n)
        gsc(js,n,m) =gsc(js,n,m) +fm*tcon(js)*work3(js,k,02)*cosnv(k,n)
 30        end do
 31       end do
 32     end do
*************
*                 RECONSTRUCT DE-ALIASED GCON
*************
        do 42 n = 0,nmax
          do 41 k = 1,nzeta
           do 40 js= 2,ns
              work3(js,k,03) = work3(js,k,03) + gcs(js,n,m)*sinnv(k,n)
              work3(js,k,04) = work3(js,k,04) + gsc(js,n,m)*cosnv(k,n)
 40        end do
 41       end do
 42     end do
        do 51 i = 1,ntheta2
          do 50 jk= 1,ns*nzeta
            gcon(jk,i) = gcon(jk,i) + work2(jk,03)*cosmu(i,m)
     >                          + work2(jk,04)*sinmu(i,m)
 50       end do
 51     end do
!
 60    end do
        return
        end subroutine alias
