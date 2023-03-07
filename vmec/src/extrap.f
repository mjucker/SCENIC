      subroutine extrap (rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc, x3,
     & x4, ns)
!***********************************************************************
C.. Implicits ..
      implicit none

      include 'name1.inc'  

      integer :: ns
!
      real, dimension(ns,0:mnd1) :: rmncc,rmnss,zmncs,zmnsc,lmncs,lmnsc
      real, dimension(0:mnd1)    :: x3, x4
C.. Local Scalars ..
      integer :: mn, n
!
      do 10 mn = 2 * nmax1, mnd1  
       rmncc (2, mn) = x3 (mn) * rmncc (3, mn) + x4 (mn) * rmncc (4, mn)
       rmnss (2, mn) = x3 (mn) * rmnss (3, mn) + x4 (mn) * rmnss (4, mn)
       zmncs (2, mn) = x3 (mn) * zmncs (3, mn) + x4 (mn) * zmncs (4, mn)
       zmnsc (2, mn) = x3 (mn) * zmnsc (3, mn) + x4 (mn) * zmnsc (4, mn)
   10 end do
      do 20 mn = nmax1, mnd1  
       lmncs (2, mn) = x3 (mn) * lmncs (3, mn) + x4 (mn) * lmncs (4, mn)
       lmnsc (2, mn) = x3 (mn) * lmnsc (3, mn) + x4 (mn) * lmnsc (4, mn)  
   20 end do
      do 30 n = 1, nmax  
        lmncs (1, n) = x3 (0) * lmncs (2, n) - lmncs (3, n)  
   30 end do
      return  
      end subroutine extrap
