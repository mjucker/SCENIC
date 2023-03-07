C-----------------------------------------------------------------------
C  15.02.88          LAST MODIFICATION 28.02.06  UHS:  LAMNEW DATA D
C-----------------------------------------------------------------------
C
      subroutine LAMNEW(I,ni,njk,lmnl,lss,nper,llampr,ftp,fpp,ml,nl,fvl,
     &                  fvli,tcos,tsin,tsc,lsx,pgttv,pgtpv,pgppv,lmnb,
     &                  al,gla,dna,dnb)
Cpara      SUBROUTINE LAMNEW ( AL, I )
C
C-----------------------------------------------------------------------
C
C     SOLVES LINEAR 2D-PDE FOR COEFFICIENTS OF LAMBDA(FL)
C
C     VARIABLES: =THETA*FTP(S) -PHI* FPP(S) + LAM(U,V,S)
C
C-----------------------------------------------------------------------
C
C 
C.. Implicits .. 
      implicit none
C 
C.. Formal Arguments .. 
      integer lmnb,njk,lss,lmnl,I,ni,nper,llampr,ml(*),nl(*),
     &        lsx(-36:72,-72:72)
      real ftp(*),fpp(*),fvl(lmnl,*),fvli(lmnl,0:*),tcos(njk,*),
     &     tsin(njk,*),tsc(lss,*),pgttv(*),pgtpv(*),pgppv(*),al(lmnl,*),
     &     gla(3,*),dna(lmnl,2,*),dnb(lmnl,2,*)
C 
C.. Local Scalars .. 
      integer jk,L,lc,ld,lp,lr,ls,MODE
      real DET,EPSMNV,FPP2,FTP2
C 
C.. External Calls .. 
      external XMINV
C 
C ... Executable Statements ...
C 
C
C     COMPUTE GEOMETRICAL FACTORS
C
      FTP2 = 2. * FTP(I)
      FPP2 = 2. * FPP(I)
C
C
C     COMPUTE MATRIX ELEMENTS AND RHS
C
      do lr = 1,lss
        do lc = 1,3
          gla(lc,lr) = 0.0
        end do
      end do
C--sh-000831      lmnl2 = 2*lmnl
      do lc = 1,lmnl
        do lr = 1,lmnl
          dna(lc,1,lr) = 0.
          dna(lc,2,lr) = 0.
          dnb(lc,1,lr) = 0.
          dnb(lc,2,lr) = 0.
        end do
      end do
C
C      do 190 lc = 1,lmnl
C      do 190 jk = 1,njk
C      gla(lc,1,jk) = pgppv(jk) * tcos(jk,lc)
C      gla(lc,2,jk) = pgtpv(jk) * tsin(jk,lc)
C 190  continue
C
C      call mxm(gla,lmnl2,tcos,njk,dna,lmnl)
C
C      do 191 lc = 1,lmnl
C      do 191 jk = 1,njk
C      gla(lc,1,jk) = pgtpv(jk) * tcos(jk,lc)
C      gla(lc,2,jk) = pgttv(jk) * tsin(jk,lc)
C 191  continue
C
C      call mxm(gla,lmnl2,tsin,njk,dnb,lmnl)
C
      do ls = 1,lss
!CDIR UNROLL=8
        do jk = 1,njk
          gla(1,ls) = gla(1,ls) + pgppv(jk)*tsc(ls,jk)
          gla(2,ls) = gla(2,ls) + pgtpv(jk)*tsc(ls,jk)
          gla(3,ls) = gla(3,ls) + pgttv(jk)*tsc(ls,jk)
        end do
      end do
C
      do lc = 1,lmnl
        do lr = 1,lmnl
          ld = lsx(ml(lc)-ml(lr),(nl(lc)-nl(lr))/nper)
          lp = lsx(ml(lc)+ml(lr),(nl(lc)+nl(lr))/nper)
          dna(lc,1,lr) = ml(lc) * ml(lr) * (gla(1,ld)+gla(1,lp))
          dna(lc,2,lr) = -ml(lc)*nl(lr)*(gla(2,ld)+gla(2,lp))
          dnb(lr,1,lc) = dna(lc,2,lr)
          dnb(lc,2,lr) = nl(lc) * nl(lr) * (gla(3,ld)+gla(3,lp))
        end do
      end do
C
      do lc = 1,lmnl
        do Lr = 1,LMNL+3
C
          AL(LC,LR) = 0.0
C
        end do
      end do
C
      do lc = 1,lmnl
        do lr = 1,lmnl
          al(lc,lr) =
     &      al(lc,lr) + dna(lc,1,lr) + dna(lc,2,lr) + dnb(lc,1,lr) +
     &      dnb(lc,2,lr)
        end do
      end do
C
      do lr = 1,lmnl
!CDIR UNROLL=8
        do JK = 1,NJK
C
          AL(LR,LMNL+1) =
     &      AL(LR,LMNL+1) - TCOS(JK,LR)*(+FTP2*PGPPV(JK)-FPP2*PGTPV(JK))
     &      - TSIN(JK,LR)*(-FPP2*PGTTV(JK)+FTP2*PGTPV(JK))
C
        end do
      end do
C
C--sh-000831c           WRITE( 6,1002) I,LR,ML(LR),NL(LR)
C--sh-000831c           WRITE( 6,1003)(AL(LR,LW),LW=1,LMNL)
C--sh-000831C           WRITE(16,1002) I,LR,ML(LR),NL(LR)
C--sh-000831C           WRITE(16,1003)(AL(LR,LW),LW=1,LMNL)
C--sh-000831C
C--sh-000831  300       CONTINUE
C--sh-000831 1002 FORMAT (1X,'AL(LR,LW):',4I3)
C--sh-000831 1003 FORMAT (1X,1P7E10.3)
C
C     SOLVE SYSTEM
C
      MODE = 0
C--sh-000831      EPSMNV = 1.E-14
C
C-CRAY/YMP      CALL MINV ( AL(1,1), LMNL, LMNL, AL(1,LMNL+2)
      call XMINV(AL(1,1),LMNL,LMNL,AL(1,LMNL+2),DET,EPSMNV,1,MODE)
C
C     STORE NEW LAMBDA
C
      do L = 1,LMNL
C
        FVL(L,I) = AL(L,LMNL+1)
        if (I.le.NI .and. LLAMPR.eq.1) then
C
 1000     format (1x'FVL',4i3,1p2e12.4)
          write (16,1000) I, L, ML(L), NL(L), FVL(L,I), FVLI(L,I)
        end if
      end do
      end
