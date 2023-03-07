        subroutine vsetup
!--------0---------0---------0---------0---------0---------0---------0-c
!
C.. Implicits ..
      implicit none
!
C.. Local Scalars ..
      integer :: l,i
!
      include 'name2.inc'
      include 'name3.inc'

      do 10 l = 1,nrztd
        rcon0(l) = czero
        zcon0(l) = czero
   10 end do
      do 20 i = 0,10
        timer(i) = czero
   20 end do
        iequi = 0
        itfsq = 0
        bscale = c1p0                                                   vac
        delbsq = c1p0                                                   vac
        ivac   = 0                                                      vac
        return
        end subroutine vsetup
