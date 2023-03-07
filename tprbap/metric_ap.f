C-----------------------------------------------------------------------
C  17.03.06        LAST MODIFICATION 17.03.06     UHS: METRICT DATA D
C-----------------------------------------------------------------------
C
      subroutine METRIC(ni,njk,lmnb,s,fpp,ftp,cjv,civ,mb,nb,tsin,tcos,r,
     &                  z,phv,rt,zt,rp,zp,rs,zs,rsq,ri,zi,phvi,bsq,bjac,
     &                  vjac,gttl,gtpl,gppl,gssl,gstl,fr,fz,fbsq,fphv,
     &                  phvp,phvt,fgparp,gparp,fgperp,gperp,
     &                  fsigmb,sigmab)
C
C     COMPUTES METRIC ELEMENTS IN BOOZERS COORDINATES, TOGETHER
C              WITH BJAC, ETC. (PRELIMINARY VERSION)
C
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
      integer njk,lmnb,ni,mb(*),nb(*)
      real s(0:*),fpp(*),ftp(*),cjv(*),civ(*),tsin(njk,*),tcos(njk,*),
     &     r(njk,0:*),z(njk,0:*),phv(njk,*),rt(njk,*),zt(njk,*),
     &     rp(njk,*),zp(njk,*),rs(njk,*),zs(njk,*),rsq(njk,*),
     &     ri(njk,0:*),zi(njk,0:*),phvi(njk,0:*),bsq(njk,0:*),
     &     bjac(njk,0:*),vjac(njk,0:*),gttl(njk,0:*),gtpl(njk,0:*),
     &     gppl(njk,0:*),gssl(njk,0:*),gstl(njk,0:*),fr(lmnb,*),
     &     fz(lmnb,*),fbsq(lmnb,*),fphv(lmnb,*),phvp(*),phvt(*)
      real, dimension(lmnb,ni)  :: fgparp ,fgperp ,fsigmb
      real, dimension(njk,ni)   ::  gparp , gperp ,sigmab
C
C.. Local Scalars ..
      integer i,JK,L
      real PHVSJK,RSJK,TBRP,TBRT,TBZP,TBZT,tsi,ZSJK
C
C ... Executable Statements ...
C
C
C     JACOBIAN AND METRIC ELEMENTS ON HALF GRID
C
C     EQUILIBRIUM QUANTITIES ON HALF INTEGER BOOZER COORDINATE GRID.
C     THETA/PHI - DERIVATIVES
C     INITIALIZE QUANTITIES. NOTE PHV INITIALIZED IN VMTOBO.
C
      do i = 1,ni
        tsi = 1. / (s(i)-s(i-1))
C
        do JK = 1,NJK
          R(JK,I) = 0.
          Z(JK,I) = 0.
          bsq(JK,I) = 0.
          RT(JK,I) = 0.
          ZT(JK,I) = 0.
          RP(JK,I) = 0.
          ZP(JK,I) = 0.
          phvp(jk) = 1.
          PHVt(JK) = 0.
          GPARP(JK,I)  = 0.
          GPERP(JK,I)  = 0.
          SIGMAB(JK,I) = 0.
        end do
C
        do L = 1,LMNB
C
          do JK = 1,NJK
C
            TBRT = -MB(L)*TSIN(JK,L)
            TBZT = MB(L) * TCOS(JK,L)
            TBRP = NB(L) * TSIN(JK,L)
            TBZP = -NB(L)*TCOS(JK,L)
            R(JK,I) = R(JK,I) + FR(L,I)*TCOS(JK,L)
            Z(JK,I) = Z(JK,I) + FZ(L,I)*TSIN(JK,L)
            BSQ(JK,I)    = BSQ(JK,I)    + FBSQ  (L,I)*TCOS(JK,L)
            GPARP(JK,I)  = GPARP(JK,I)  + FGPARP(L,I)*TCOS(JK,L)
            GPERP(JK,I)  = GPERP(JK,I)  + FGPERP(L,I)*TCOS(JK,L)
            SIGMAB(JK,I) = SIGMAB(JK,I) + FSIGMB(L,I)*TCOS(JK,L)
            PHV(JK,I) = PHV(JK,I) + FPHV(L,I)*TSIN(JK,L)
            RT(JK,I) = RT(JK,I) + FR(L,I)*TBRT
            ZT(JK,I) = ZT(JK,I) + FZ(L,I)*TBZT
            RP(JK,I) = RP(JK,I) + FR(L,I)*TBRP
            ZP(JK,I) = ZP(JK,I) + FZ(L,I)*TBZP
            PHVT(JK) = PHVT(JK) + FPHV(L,I)*TBZT
            PHVP(JK) = PHVP(JK) + FPHV(L,I)*TBZP
          end do
C
        end do
C
        do JK = 1,NJK
C
C     S - DIFFERENCING
C
          RSJK = (RI(JK,I)-RI(JK,I-1)) * TSI
          ZSJK = (ZI(JK,I)-ZI(JK,I-1)) * TSI
          PHVSJK = (PHVI(JK,I)-PHVI(JK,I-1)) * TSI
          RS(JK,I) = RSJK
          ZS(JK,I) = ZSJK
          RSQ(JK,I) = R(JK,I)**2
C--sh-000831         rsqsjk     = ( rsq(jk,i) - rsq(jk,i-1)) * tsi
C--sh-000831         rsqth      = r(jk,i) * rt(jk,i)
C--sh-000831         rsqph      = r(jk,i) * rp(jk,i)
C
C     BJAC COMPUTED FROM BSQ GOT FROM DIRECT INTEGRATION:
C
          BJAC(JK,I) = (FPP(I)*CJV(I)-FTP(I)*CIV(I)) / BSQ(JK,I)
     &                 / SIGMAB(JK,I)
C
C     L(OWER) M.E.
C
          GTTL(JK,I) =
     &      (rt(jk,i)*rt(jk,i)+rsq(jk,i)*phvt(jk)*phvt(jk)+
     &       zt(jk,i)*zt(jk,i)) / bjac(jk,i)
C
C     GPPL AND GTPL COMPUTED USING EQUILIBRIUM STUFF AND GTTL
C
          GTPL(JK,I) = CJV(I)/(FTP(I) * SIGMAB(JK,I))
     &               - FPP(I)/FTP(I)*GTTL(JK,I)
C
          GPPL(JK,I) =
     &      -(FPP(I)*CJV(I)+FTP(I)*CIV(I))/(SIGMAB(JK,I) * FTP(I)**2)
     &       +     (FPP(I)/FTP(I))**2*GTTL(JK,I)
C
C     FROM HERE NO STRICT SQUARE BEFORE AVERAGING!
          GSTL(JK,I) =
     &      (rsjk*rt(jk,i)+rsq(jk,i)*phvsjk*phvt(jk)+zsjk*zt(jk,i)) /
     &      bjac(jk,i)
          GSSL(JK,I) =
     &      (rsjk*rsjk+rsq(jk,i)*phvsjk*phvsjk+zsjk*zsjk) / bjac(jk,i)
C
C--THE BOOZER JACOBIAN COMPUTED FROM THE GEOMETRY STORED IN vjac.
          vjac(jk,i) =
     &      r(jk,i) *
     &      (rsjk*(phvt(jk)*zp(jk,i)-phvp(jk)*zt(jk,i))+
     &       phvsjk*(zt(jk,i)*rp(jk,i)-zp(jk,i)*rt(jk,i))+
     &       zsjk*(rt(jk,i)*phvp(jk)-rp(jk,i)*phvt(jk)))
C          vjac(jk,i) = 0.5 *
C     $ (rsqsjk * (phvt(jk) * zp(jk,i) - phvp(jk) * zt(jk,i))
C     $ +phvsjk * (  zt(jk,i) * rsqph   - zp(jk,i) * rsqth)
C     $ +zsjk * (   rsqth * phvp(jk) - rsqph * phvt(jk)))
Cbound      if (i.ge.ni-2 .and. jk.le.nj)
Cbound     $ write (16,1009)i,vjac(jk,i),r(jk,i),rsjk,zsjk,phvsjk,rt(jk,i)
Cbound     $,zt(jk,i)
C
C     IF(I.LE.6.AND.JK.LE.13)WRITE(16,1002)
C    $   I,JK,RSJK**2, RSQ(JK,I)*PHVSJK**2, ZSJK**2
C1002 FORMAT(1X,'GSSL',2I3,1P3E15.6)
C
        end do
C
      end do
      end subroutine metric
