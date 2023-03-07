C-----------------------------------------------------------------------
C  17.11.89        LAST MODIFICATION 28.02.06     UHS: lamcal DATA   D
C-----------------------------------------------------------------------
C
      subroutine lamcal(ni,njk,lmnl,lss,nper,nprocs,llampr,dtdp,s,ftp,
     &                  fpp,civ,cjv,vvp,wmagv,ml,nl,fvl,fvli,tcos,tsin,
     &                  tsss,zv,zvt,rvq,rvqt,vjac,vl,gttv,gtpv,gppv,vbp,
     &                  vbt,vbsq,tsc,lsx,lmnb,vlt,vlp,psivt,psivp,
     &                  sigmav,parpv,parpav,pgppv,pgttv,pgtpv,al,gla,
     &                  dna,dnb)
C
C      common /lamcom/    PGPPV(njk ),PGTTV(njk ),PGTPV(njk )
Cl1t      parameter (njkloc=4320)
C      parameter (njkloc=2880)
C
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
      integer lmnb,njk,lss,lmnl,ni,nper,nprocs,llampr,ml(*),nl(*),
     &        lsx(-36:72,-72:72)
      real, dimension(njk,0:ni) :: sigmav  ,parpv
      real, dimension(ni)       :: parpav
      real dtdp,s(0:*),ftp(*),fpp(*),civ(*),cjv(*),vvp(*),wmagv(*),
     &     fvl(lmnl,*),fvli(lmnl,0:*),tcos(njk,*),tsin(njk,*),
     &     tsss(njk,*),zv(njk,0:*),zvt(njk,0:*),rvq(njk,0:*),
     &     rvqt(njk,0:*),vjac(njk,0:*),vl(njk,*),gttv(njk,0:*),
     &     gtpv(njk,0:*),gppv(njk,0:*),vbp(njk,0:*),vbt(njk,0:*),
     &     vbsq(njk,0:*),tsc(lss,*),vlt(*),vlp(*),psivt(*),psivp(*),
     &     pgppv(*),pgttv(*),pgtpv(*),al(lmnl,*),gla(3,*),dna(lmnl,2,*),
     &     dnb(lmnl,2,*)
C
C.. Local Scalars ..
      integer i,JK,L
      real :: dvol,FPPI,FTPI,RVQS,tsi,ZVS, sigovjac
C
C.. External Calls ..
      external LAMNEW
C
C ... Executable Statements ...
C
C
C     COMPUTE VJAC ETC and geometrical factors ON HALFGRID
C
Cpara      DO 200 I = NI,1,-1
Cpara      DO 300 I = npr,NI,nprocs
      do i = 1,ni
C
        tsi = 1. / (s(i)-s(i-1))
        do JK = 1,NJK
C
          sigovjac = sigmav(jk,i) / vjac(jk,i)
          ZVS = TSI * (ZV(JK,I)-ZV(JK,I-1))
          RVQS = TSI * (RVQ(JK,I)-RVQ(JK,I-1))
C           VJAC(JK,I) = 0.25* ( (RVQT(JK,I)+RVQT(JK,I-1))*ZVS
C    $                           -RVQS      *(ZVT(JK,I)+ZVT(JK,I-1)))
C    $                       * TPONP
          pGTTV(JK) = 0.5 * (GTTV(JK,I)+GTTV(JK,I-1)) * sigovjac
          pGTPV(JK) =-0.5 * (GTPV(JK,I)+GTPV(JK,I-1)) * sigovjac
          pGPPV(JK) = 0.5 * (GPPV(JK,I)+GPPV(JK,I-1)) * sigovjac
C
        end do
C--sh-000831  200    CONTINUE
C
C     COMPUTE NEW COEFFICIENTS FOR LAMBDA, EXTRAPOLATE FOR I=1,2
C
        call LAMNEW(I,ni,njk,lmnl,lss,nper,llampr,ftp,fpp,ml,nl,fvl,
     &              fvli,tcos,tsin,tsc,lsx,pgttv,pgtpv,pgppv,lmnb,al,
     &              gla,dna,dnb)
C
C     COMPUTE DERIVATIVES AND COMBIN. WITH LAMBDA
C     INITIALIZE R AND Z ON HALF INTEGER GRID.
C
        do JK = 1,NJK
C
          VL(JK,I) = 0.
          VLT(JK) = 0.
          VLP(JK) = 0.
        end do
C
        do L = 1,LMNL
          do JK = 1,NJK
            VL(JK,I) = VL(JK,I) + FVL(L,I)*tsss(jk,l)
            VLT(JK) = VLT(JK) + FVL(L,I)*TCOS(JK,L)
            VLP(JK) = VLP(JK) + FVL(L,I)*TSIN(JK,L)
          end do
        end do
C
        FTPI = FTP(I)
        FPPI = FPP(I)
C
        do JK = 1,NJK
C
          PSIVT(JK) = FTPI + VLT(JK)
          PSIVP(JK) = -FPPI + VLP(JK)
        end do
C
C     COMPUTE sigma*B(SUB)T=VBT, sigma*B(SUB)P=VBP, J, I, ETC.
C
        CIV(I) = 0.
        CJV(I) = 0.
        VVP(I) = 0.
        WMAGV(I) = 0.
        dvol = (s(i)-s(i-1)) * dtdp
        parpav(i) = 0.
C
        do JK = 1,NJK
C
          VBP(JK,I) = (pGPPV(JK)*PSIVT(JK)+pgtpv(jk)*psivp(jk))
C
          VBT(JK,I) = -(pGTPV(JK)*PSIVT(JK)+pgttv(jk)*psivp(jk))
C
          VBSQ(JK,I) =
     &      (vbp(jk,i)*psivt(jk)-vbt(jk,i)*psivp(jk)) / vjac(jk,i)
     &      / sigmav(jk,i)
C
          CIV(I) = CIV(I) - DTDP*VBP(JK,I)
          CJV(I) = CJV(I) + DTDP*VBT(JK,I)
C
          VVP(I) = VVP(I) - DTDP*VJAC(JK,I)
          WMAGV(I) = WMAGV(I) + 0.5*VBSQ(JK,I)*VJAC(JK,I)*DVOL
!
          PARPAV(I) = PARPAV(I) + DTDP * VJAC(JK,I) * PARPV(JK,I)
C
        end do
C
      end do
      end subroutine lamcal
