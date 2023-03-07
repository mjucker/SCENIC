C-----------------------------------------------------------------------
C  29.03.88        LAST MODIFICATION 17.03.06     UHS: VMTOBO DATA D
C-----------------------------------------------------------------------
C
      subroutine VMTOBO(ni,njk,nj,lmnb,lmnl,nper,lmnb0,lvmtpr,tpi,dtdp,
     &                  dph,dth,fpp,ftp,cjv,civ,mb,nb,ml,nl,tsin,tcos,
     &                  vbp,vbt,vbsq,vjac,r,z,vl,phv,fvl,fr,fz,fphv,
     &                  fphvt,fphvp,fbsq,fbjac,fvalf,fvgam,fvb,fvq,
     &                  valft,vgamt,valf,vgam,vjtobj,thv,vq,sqgbsq,
     &                  sigmav,parpv,perpv,tauv,fsigmb,ftaub,fgparp,
     &                  fgperp)
C
C     COMPUTES VALF, FVALF, VGAM, FVGAM, VQ, FVQ AND VJTOBJ FOR MAPPING
C        GIVES FR, FZ, FPHV, R, Z, PHV IN BC
C              BY DIRECT INTEGRATION
!     FGPARP, FGPERP ARE FOURIER AMPLITUDES OF JACOBIAN TIMES HOT PRESSURE
!     RADIAL GRADIENTS.
C
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
      integer njk,lmnb,lmnb0,ni,nj,lmnl,nper,lvmtpr,mb(*),nb(*),ml(*),
     &        nl(*)
      real tpi,dtdp,dph,dth,fpp(*),ftp(*),cjv(*),civ(*),tsin(njk,*),
     &     tcos(njk,*),vbp(njk,0:*),vbt(njk,0:*),vbsq(njk,0:*),
     &     vjac(njk,0:*),r(njk,0:*),z(njk,0:*),vl(njk,*),phv(njk,*),
     &     fvl(lmnl,*),fr(lmnb,*),fz(lmnb,*),fphv(lmnb,*),fphvt(lmnb,*),
     &     fphvp(lmnb,*),fbsq(lmnb,*),fbjac(lmnb,0:*),fvalf(*),fvgam(*),
     &     fvb(*),fvq(*),valft(*),vgamt(*),valf(*),vgam(*),vjtobj(*),
     &     thv(*),vq(*),sqgbsq(*)
      real, dimension(njk,0:ni) :: sigmav ,parpv  ,perpv  ,tauv
      real, dimension(lmnb,ni)  :: fsigmb ,fgparp ,fgperp ,ftaub
C
C.. Local Scalars ..
      integer I,J,JK,K,L,LB,LL,NBPFP
      real AT,cosat,FPPI,FTPI,TBZT,TWO
C
C.. Intrinsic Functions ..
      intrinsic MOD, SIN, cos, MERGE
C
C ... Executable Statements ...
C
C
C     if (i.eq.0) go to 300
C     if (i.eq.1)WRITE( 6,1001)
C     if (i.eq.1)WRITE(16,1001)
C 1001 FORMAT(/1X,'VMTOBO: I,J,VL,VQ,VALF,VGAM'/)
C
c..DETERMINE THE INDEX FOR WHICH mb=0 AND nb=0
      do l=1,lmnb
        if (mb(l).eq.0 .and. nb(l).eq.0) then
        lmnb0 = l
        endif
      end do
C
C     COMPUTE FVQ AND VQ FROM VBT AND VBP
C
      do I = 1,NI
C
        do L = 1,LMNB
          FVB(L) = 0.
          FVQ(L) = 0.
C
          if (NB(L) .ne. 0) then
            do JK = 1,NJK
              FVB(L) = FVB(L) + 2.*VBP(JK,I)*TCOS(JK,L)
            end do
CPERF   FVQ(L) = - DTDP / NB(L) * FVB(L) / TPI
            FVQ(L) = -DTDP/NB(L)*FVB(L)
C
          elseif (NB(L).eq.0 .and. MB(L).ne.0) then
            do JK = 1,NJK
              FVB(L) = FVB(L) + 2.*VBT(JK,I)*TCOS(JK,L)
            end do
CPERF   FVQ(L) =   DTDP / MB(L) * FVB(L) / TPI
            FVQ(L) = DTDP / MB(L) * FVB(L)
          end if
C
        end do
C--sh-000831  140   CONTINUE
C
Cpara      DO 150 I=1,NI
        do JK = 1,NJK
C
          VQ(JK) = 0.
C
          do L = 1,LMNB
C
            VQ(JK) = VQ(JK) + FVQ(L)*TSIN(JK,L)
          end do
C
        end do
C
C     ALL MAPPING FUNCTIONS COMPUTED ON HALF GRID
C
Cpara      DO 200 I=1,NI
C
        FPPI = FPP(I)
        FTPI = FTP(I)
        SQGBSQ(I) = FPPI*CJV(I) - FTPI*CIV(I)
C
        do LL = 1,LMNL
!          do LB = 1,LMNB
!            if (ML(LL).eq.MB(LB) .and. NL(LL).eq.NB(LB)) then
              LB = MERGE(LL,LL+1,LL<lmnb0)
              FVALF(LB) = (FPPI*FVQ(LB)-CIV(I)*FVL(LL,I)) / sqgbsq(i)
              FVGAM(LB) = (FTPI*FVQ(LB)-CJV(I)*FVL(LL,I)) / sqgbsq(i)
!            end if
!          end do
        end do
            FVALF(lmnb0) = 0.
            FVGAM(lmnb0) = 0.
C
        do JK = 1,NJK
C
          VALF(JK) = (FPPI*VQ(JK)-CIV(I)*VL(JK,I)) / SQGBSQ(I)
          VGAM(JK) = (FTPI*VQ(JK)-CJV(I)*VL(JK,I)) / SQGBSQ(I)
          VJTOBJ(JK) = SIGMAV(JK,I) * VJAC(JK,I) * VBSQ(JK,I)/SQGBSQ(I)
C        WRITE (16,1002)I,JK,VL(JK,I),VQ(JK),VALF(JK),VGAM(JK)
        end do
C 1002 FORMAT (1X,2I3,1P4E15.6)
C
C        IF (I.LE.8) THEN
C        DO 208 L = 1,LMNB
C 208    WRITE (16, *  ) I,ML(L),MB(L),FVL(L),FVQ(L)
C
C
C        END IF
        do JK = 1,NJK
C
C     PRESET PHV ETC
C
          K = ((JK-1)/NJ)
          J = MOD((JK-1),NJ)
CPERF  PHV(JK,I) = DPH * K
CPERF  THV(JK  ) = DTH * J
          PHV(JK,I) = DPH * K * TPI / NPER
          THV(JK) = DTH * J * TPI
          VGAMT(JK) = 0.
          VALFT(JK) = 0.
C
        end do
C
C     COMPUTE VALFT, VGAMT
C
        do L = 1,LMNB
          do JK = 1,NJK
CPERF    TBZT = TPI * MB(L) * TCOS(JK,L)
            TBZT = MB(L) * TCOS(JK,L)
            VGAMT(JK) = VGAMT(JK) + FVGAM(L)*TBZT
            VALFT(JK) = VALFT(JK) + FVALF(L)*TBZT
          end do
        end do
C
C
C     DIRECT FOURIER INTEGRALS FOR BC COEFFICIENTS
C
        do L = 1,LMNB
C
          FR(L,I) = 0.
          FZ(L,I) = 0.
          FBSQ(L,I) = 0.
          fbjac(l,i) = 0.
          FPHV(L,I) = 0.
          FPHVT(L,I) = 0.
          FPHVP(L,I) = 0.
          FTAUB(L,I) = 0.
          FSIGMB(L,I)= 0.
          FGPARP(L,I)= 0.
          FGPERP(L,I)= 0.
C
          TWO = 2.
          if (MB(L).eq.0 .and. NB(L).eq.0) then
            TWO = 1.
          end if
          do JK = 1,NJK
C
CPERF       AT = TPI * (MB(L)*( THV(JK  ) + VALF(JK  ) )
            AT = (MB(L)*(THV(JK)+VALF(JK))-NB(L)*(PHV(JK,I)+VGAM(JK)))
            cosat = cos(at)
            FR(L,I)     = FR    (L,I) + R(JK,I)*VJTOBJ(JK)*COSat
            FZ(L,I)     = FZ    (L,I) + Z(JK,I)*VJTOBJ(JK)*SIN(AT)
            FBSQ(L,I)   = FBSQ  (L,I) + VBSQ(JK,I)*VJTOBJ(JK)*COSat
            fbjac(l,i)  = fbjac (l,i) + vjac(jk,i)*cosat
            FPHVT(L,I)  = FPHVT (L,I) - VGAMT(JK)*COSat
            FPHVP(L,I)  = FPHVP (L,I) + (1.+VALFT(JK))*COSat
            FSIGMB(L,I) = FSIGMB(L,I) + SIGMAV(JK,I)*VJTOBJ(JK)*COSat
            FTAUB(L,I)  = FTAUB (L,I) + TAUV  (JK,I)*VJTOBJ(JK)*COSat
            FGPARP(L,I) = FGPARP(L,I) + vjac(jk,i)*PARPV (JK,I)*COSat
            FGPERP(L,I) = FGPERP(L,I) + vjac(jk,i)*PERPV (JK,I)*COSat
C
          end do
C
          FR(L,I) = DTDP * FR(L,I) * TWO
          FZ(L,I) = DTDP * FZ(L,I) * TWO
          FBSQ(L,I)   = DTDP * FBSQ(L,I)   * TWO
          Fbjac(l,i)  = dtdp * fbjac(l,i)  * two
          FPHVT(L,I)  = DTDP * FPHVT(L,I)  * TWO
          FPHVP(L,I)  = DTDP * FPHVP(L,I)  * TWO
          FSIGMB(L,I) = DTDP * FSIGMB(L,I) * TWO
          FTAUB(L,I)  = DTDP * FTAUB(L,I)  * TWO
          FGPARP(L,I) = DTDP * FGPARP(L,I) * TWO
          FGPERP(L,I) = DTDP * FGPERP(L,I) * TWO
C
        end do
C
        do L = 1,LMNB
C
          if (MB(L).ne.0 .or. NB(L).ne.0) then
            if (MB(L) .ne. 0) then
CPERF       FPHV(L,I) = + 1./(TPI*MB(L)) * FPHVT(L,I)
              FPHV(L,I) = +1./(MB(L))*FPHVT(L,I)
            else
CPERF       FPHV(L,I) = - 1./(TPI*NB(L)) * FPHVP(L,I)
              FPHV(L,I) = -1./(NB(L))*FPHVP(L,I)
            end if
          end if
C
        end do
C
        if ((I.le.4.or.I.gt.NI-5) .and. LVMTPR.eq.1) then
C         IF(I.GT.NI-5          ) THEN
          do L = 1,LMNB
            NBPFP = NB(L) / NPER
 1000       format (1x,'BOFOU',4i3,1p4e15.6)
CPERF     WRITE (16,1000) I,L,MB(L),NB(L),FR(L,I),FZ(L,I),FPHV(L,I)
            write (16,1000)
     &            I, L, MB(L), NBPFP, FR(L,I), FZ(L,I), FPHV(L,I),
     &            FBSQ(L,I)
CPERF     WRITE ( 6,1000) I,L,MB(L),NB(L),FR(L,I),FZ(L,I),FPHV(L,I)
            write (6,1000)
     &            I, L, MB(L), NBPFP, FR(L,I), FZ(L,I), FPHV(L,I),
     &            FBSQ(L,I)
CPERF     WRITE ( 6,1000) I,L,MB(L),NB(L),FR(L,I),FZ(L,I),FPHV(L,I)
          end do
        end if
      end do
      end subroutine vmtobo
