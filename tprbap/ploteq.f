C-----------------------------------------------------------------------
C  22.06.88        LAST MODIFICATION 23.11.05     UHS: PLOTEQ DATA D
C-----------------------------------------------------------------------
C
      subroutine PLOTEQ(LMAP)
C 
      use tprb_param
      use triggs
      use vmcfou
      use boofou
C.. Implicits .. 
      implicit none
C
C     WRITES EQUILIBRIUM COEFFICIENTS TO UNIT 17
C     IF LMAP =1: ORIGINAL DATA ARE USED,
C     IF LMAP =2: BC AFTER THE MAPPING ARE USED
C
!      include 'parbal.inc'
C 
C.. Formal Arguments .. 
C.. In/Out Status: Read, Not Written ..
      integer LMAP
C 
C.. Local Scalars .. 
      integer I,L
!      include 'tprcom.bal'
C 
C ... Executable Statements ...
C 
C
      if (LMAP .eq. 1) then
C
 1000   format (
     &  /' TERPSICHORE EQUILIBRIUM DATA'//,' LMAP = ',i2,' NI = ',i3,
     &  ' LMN =',i3,//)
C
        write (17,1000) LMAP, NI, LMNV !
                                       !
C
C
        do I = 0,NI
C
          do L = 1,LMNV
 1001       format (1x,4i5,1p2e25.14)
C
            write (17,1001) I, L, MX(L), NX(L), FRV(L,I), FZV(L,I)
C
          end do
C
        end do
      else
C
        write (17,1000) LMAP, NI, LMNB
C
        do I = 0,NI
C
          do L = 1,LMNB
C
            write (17,1001) I, L, MB(L), NB(L), FRI(L,I), FZI(L,I)
C
          end do
C
        end do
      end if
      end
