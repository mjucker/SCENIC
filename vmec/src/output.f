        subroutine output(ierflag)
!--------0---------0---------0---------0---------0---------0---------0-c
!
C.. Implicits ..
      implicit none
!
      include 'name2.inc'
!
      integer :: ierflag
!
        iequi = 1
        call funct3d(ierflag)
        print 10, ijacob
        write(3,10)ijacob
 10     format(/,'  NUMBER OF JACOBIAN RESETS = ',i4)
        print 20, timer(0),timer(1)
        write(3,20)timer(0),timer(1)
 20     format (/,
     &  '  TOTAL COMPUTATIONAL TIME :       ',1pe10.2,' SECONDS',/
     &  '  TIME IN VACUUM LOOP :            ',1pe10.2,' SECONDS',/)
        return
        end subroutine output
