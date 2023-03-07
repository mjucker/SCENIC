C-----------------------------------------------------------------------
C  18.12.89        LAST MODIFICATION 01.02.94     UHS: trgfun DATA   D
C-----------------------------------------------------------------------
C
      subroutine trgfun(njk,lmnl,ml,nl,pol,tor,tcos,tsin,tsss)
C
C 
C.. Implicits .. 
      implicit none
C 
C.. Formal Arguments .. 
      integer njk,lmnl,ml(*),nl(*)
      real pol(*),tor(*),tcos(njk,*),tsin(njk,*),tsss(njk,*)
C 
C.. Local Scalars .. 
      integer JK,l
      real Tc
C 
C.. Intrinsic Functions .. 
      intrinsic COS, sin
C 
C ... Executable Statements ...
C 
C
C     COMPUTE TRIGFCT FOR LAMBDA
C
      do l = 1,lmnl
        do JK = 1,NJK
          Tc = COS(ML(L)*POL(JK)-NL(L)*TOR(JK))
          TCOS(JK,L) = +ML(L)*TC
          TSIN(JK,L) = -NL(L)*Tc
          TSss(JK,L) = sin(ml(l)*pol(jk)-NL(L)*Tor(jk))
        end do
      end do
      end
