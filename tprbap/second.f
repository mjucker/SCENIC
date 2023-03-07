C****************************************************************************A
      subroutine second(tt)
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
C.. In/Out Status: Not Read, Overwritten ..
!      real :: tt, tt1, ra(2)
       real :: tt
C
C.. External Functions ..
!      real etime
!      external etime
C
C ... Executable Statements ...
C
!      call etime(tt)
      call cpu_time(tt)
!      tt1=etime(ra)
!      tt=ra(1)
      end subroutine second
