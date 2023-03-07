
        subroutine getfsq(gcr,gcz,gnormr,gnormz,gnorm,mprecon)
!
C.. Implicits ..
      implicit none
      include 'name2.inc'
        real, dimension(ns,0:mnd1,2) ::  gcr, gcz
        real                         :: gnormr,gnormz,gnorm
        integer :: mprecon
C.. Local Scalars ..
        integer :: js,mn,jsmax
        gnormr = czero
        gnormz = czero
        jsmax = (ns-1) + mprecon
        do 15 mn = 0,mnd1
        do 10 js = 1,jsmax
        gnormr = gnormr + gnorm*(gcr(js,mn,1)**2 + gcr(js,mn,2)**2)
        gnormz = gnormz + gnorm*(gcz(js,mn,1)**2 + gcz(js,mn,2)**2)
 10     end do
 15     end do
        return
        end subroutine getfsq
