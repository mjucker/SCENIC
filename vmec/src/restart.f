        subroutine restart
!--------0---------0---------0---------0---------0---------0---------0-
!
C.. Implicits ..
      implicit none
      include 'name2.inc'
C.. Local Scalars ..
      integer :: l
!
C.. External Calls ..
      external scopy
        goto (10,20,20)irst
 10     call scopy(neqs,xc,1,xstore,1)
        return
 20   do 30 l = 1,neqs
        xcdot(l) = czero
        xc(l) = xstore(l)
   30 end do
        delt = delt*(cp96*(irst-2) + cp9*(3-irst))
        irst = 1
        return
        end subroutine restart
