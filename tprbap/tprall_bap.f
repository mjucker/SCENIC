CDIRD NOLIST
C-----------------------------------------------------------------------
C  07.07.04        LAST MODIFICATION 24.03.06   UHS: tprall_bap DATA   D
C-----------------------------------------------------------------------
C
      subroutine TPRALL_BAP
C
      use tprb_param
C
      use physeq
      use radfct
      use triggs
      use paprpl
      use vmcfou
      use vmcrsp
      use boofou
      use boorsp
      use colamn
      use lambda
      use mappin
      use boomet
      use bderiv
      use comerc
      use merfou
      use balfou
C
      save
C
C.. ALLOCATE BLOCKS ..
C
C
C... FORMERLY Variables in Common Block /BDERIV/ ...
      allocate (BT(NJK,0:NI),BP(NJK,0:NI),BSQ(NJK,0:NI),SIGBS(NJK,0:NI),
     &          stat=ierr2)
!      common /BDERIV/ BT,BP,BSQ,SIGBS
C... In/Out Status of Common:
C     For FR to FZ: Maybe Read, Maybe Written
C     For FPHV: Not Read, Not Written
C     For FRI to FZI: Maybe Read, Not Written
C     For FPHVI to FPARK: Not Read, Not Written
C
C... FORMERLY Variables in Common Block /BOOFOU/ ...
      allocate(FR(MLMNB,NI),FZ(MLMNB,NI),FPHV(MLMNB,NI),FRI(MLMNB,0:NI),
     &         FZI(MLMNB,0:NI) ,FPHVI(MLMNB,0:NI)      ,FPHVT(MLMNB,NI),
     &         FPHVP(MLMNB,NI) ,fbsq (mlmnb,ni)      ,FBJAC(MLMNB,0:NI),
     &         FSIGMB(MLMNB,NI),FTAUB(MLMNB,NI)       ,FGPARP(MLMNB,NI),
     &         FGPERP(MLMNB,NI),FSIGBS(MLMNB,0:NI)    ,FPARJP(MLMNB,NI),
     &         FPARK(MLMNB,NI) ,stat=ierr3)
!      common /BOOFOU/ FR,FZ,FPHV,FRI,FZI,FPHVI,FPHVT,FPHVP,fbsq,FBJAC,
!     &                FBS
C
C... FORMERLY Variables in Common Block /BOOMET/ ...
      allocate    (BJAC(NJK,0:NI),BJACS(NJK,0:NI),GTTL(NJK,0:NI),
     &     GTPL(NJK,0:NI),GPPL(NJK,0:NI),GSSL(NJK,0:NI),
     &     GSTL(NJK,0:NI),GSsu(NJK,0:NI),stat=ierr4)
!      common /BOOMET/ BJAC,BJACS,GTTL,GTPL,GPPL,GSSL,GSTL,GSsu
C... In/Out Status of Common:
C     For R to Z: Maybe Read, Maybe Written
C     For phv: Not Read, Not Written
C     For RI to phvi: Not Read, Maybe Written
C     For RT to PHVT: Not Read, Not Written
C
C... FORMERLY Variables in Common Block /BOORSP/ ...
      allocate ( R(NJK,0:NI)  ,Z(NJK,0:NI)   ,phv(njk,ni) ,RI(NJK,0:NI),
     &          ZI(NJK,0:NI)  ,phvi(njk,0:ni),RT(NJK,NI)  ,ZT(NJK,NI),
     &          rs(njk,ni)    ,RP(NJK,NI)    ,ZP(NJK,NI)  ,zs(njk,ni),
     &          SIGMAB(NJK,NI),GPARP(NJK,NI) ,GPERP(NJK,NI),
     &          PARJP(NJK,NI) ,PARKUR(NJK,NI),TAUB(NJK,NI) ,
     &          RSq(NJK,NI)   ,PHVP(njk)     ,PHVT(njk)    , stat=ierr5)
!      common /BOORSP/ R,Z,phv,RI,ZI,phvi,RT,ZT,rs,RP,ZP,zs,RSq,
!     &                PHVP,PHVT
C... In/Out Status of Common:
C     For AL to dnb: Maybe Read, Maybe Written
C
C... FORMERLY Variables in Common Block /COLAMN/ ...
      allocate  (AL(mlmnl,mlmnl+3),gla(3,lssl),dna(mlmnl,2,mlmnl),
     &     dnb(mlmnl,2,mlmnl),stat=ierr6)
!     common /COLAMN/ AL,gla,dna,dnb
C
C... FORMERLY Variables in Common Block /COMERC/ ...
      allocate(BGRADG(njk),SIGBXG(njk),BALDP(NJK),BALDS(njk),BALCP(njk),
     &         SHEARL(njk),BALCS(njk) ,BALCQ(njk),BALNP(njk),BALNS(njk),
     &         CURVDS(njk),CURVDN(njk),BALDN(njk),BALNQ(njk),stat=ierr7)
!     common /COMERC/ BGRADG
C
C...FORMERLY Variables in Common Block /MERFOU/
      allocate
     &    (fcurvs(mlmnb), flshea(mlmnb), fcurvn(mlmnb), fgssu(mlmnb)
     &    ,stat=ierr8)
!      COMMON /MERFOU/ fcurvs,flshea,fcurvn,fgssu
C
C...FORMERLY Variables in Common Block /BALFOU/
       allocate
     &    ( fbaldp(mlmnb), fbalds(mlmnb), fbaldn(mlmnb) ,fbaldt(mlmnb) ,
     &      fbalcp(mlmnb), fbalcs(mlmnb), fbalcq(mlmnb) ,
     &      fbalnp(mlmnb), fbalns(mlmnb), fbalnq(mlmnb), stat=ierr9)
!     COMMON /BALFOU/ fbaldp, fbalds, fbalcp, fbalcs, fbalcq
C
C... In/Out Status of Common:
C     For VL to PGTPV: Maybe Read, Maybe Written
C
C... FORMERLY Variables in Common Block /LAMBDA/ ...
      allocate (VL(NJK,NI),vlt(NJK),vlp(NJK),PSIVT(NJK),PSIVP(NJK),
     &     PGPPV(NJK),PGTTV(NJK),PGTPV(NJK),stat=ierr10)
!      common /LAMBDA/ VL,vlt,vlp,PSIVT,PSIVP,PGPPV,PGTTV,PGTPV
C
C... FORMERLY Variables in Common Block /MAPPIN/ ...
      allocate (VQ(NJK),sqgbsq(0:ni),VALF(NJK),VGAM(NJK),THV(NJK),
     &     VALFT(NJK),VGAMT(NJK),VJTOBJ(NJK),stat=ierr11)
!      common /MAPPIN/ VQ,sqgbsq,VALF,VGAM,THV,VALFT,VGAMT,VJTOBJ
C... In/Out Status of Common:
C     For AM: Not Read, Not Written
C     For PTH: Maybe Read, Not Written
C     For PVP to PVPI: Maybe Read, Maybe Written
C     For FPP to FTP: Maybe Read, Not Written
C     For PPI to AIOTA: Not Read, Not Written
C     For VVP: Maybe Read, Maybe Written
C     For f0 to ftppp: Not Read, Not Written
C     For WMAGV: Maybe Read, Maybe Written
C     For WMAG: Not Read, Not Written
C     For CIV to CJVP: Maybe Read, Maybe Written
C     For CI to VPP: Not Read, Not Written
C     For EQUIV: Maybe Read, Maybe Written
C     For EQUI: Not Read, Not Written
C     For PARPVI to PARPAV: Not Read, Maybe Written
C     For s: Maybe Read, Maybe Written
C
C... FORMERLY Variables in Common Block /RADFCT/ ...
      allocate(          AM(  NI),   PTH(  NI),  PVP(  NI),PVPI(  NI)
     &              ,   FPP(NI+1),  FTP(NI+1),  PPI(  NI),  PP(  NI)
     &              ,  FPPP(NI+1),  FTPP(NI+1),AIOTA(  NI), VVP(  NI)
     &              ,    f0(  ni),    f1(  ni),   f2(  ni),  f3(  ni)
     &              , fpppp(  ni), ftppp(  ni),WMAGV(  NI),WMAG(  NI)
     &              ,   CIV(  NI),   CJV(  NI), CIVP(  NI),CJVP(  NI)
     &              ,   CI (  NI),   CJ (  NI), CIP (  NI),CJP (  NI)
     &              ,  CIPI(  NI),  CJPI(  NI), ppp (  ni),glm (ni,4)
     &              ,    VP(  NI),   VPP(  NI),EQUIV(  NI),EQUI(  NI)
     &              ,PARPVI(  NI),PARPAV( NI) , s (0:ni),stat=ierr12)
!      common /RADFCT/ AM,PTH,PVP,PVPI,FPP,FTP,PPI,PP,FPPP,FTPP,AIOTA,
!     &                VVP,f0,f1,f2,f3,fpppp,ftppp,WMAGV,WMAG,CIV,CJV,
!     &                CIVP,CJVP,CI,CJ,CIP,CJP,CIPI,CJPI,ppp,glm,VP,VPP,
!     &                EQUIV,EQUI,s
C
C... In/Out Status of Common:
C     For TCOS to TOR: Maybe Read, Maybe Written
C     For TSS: Not Read, Not Written
C     For TSC: Maybe Read, Maybe Written
C     For MX to LMNV: Maybe Read, Not Written
C     For LMNL: Read, Maybe Written
C     For LMNB: Read, Not Written
C     For MS to LMNS: Not Read, Not Written
C     For lfrz: Not Read, Maybe Written
C     For lfx: Not Read, Not Written
C     For lsx to nsx: Maybe Read, Not Written
C     For lss: Read, Not Written
C     For cospar to parity: Not Read, Not Written
C
C... FORMERLY Variables in Common Block /TRIGGS/ ...
      allocate(MX(MLMNV)  ,NX(MLMNV)  ,ML(MLMNB) ,NL(MLMNB),MB(MLMNB)
     &        ,NB(MLMNB)  ,lfrz(0:36,-36:36)     ,lfx(-36:72,-72:72)
     &        ,lsx(-36:72,-72:72)     ,msx(lssl) ,nsx(lssl),stat=ierr13)
      allocate( TCOS(NJK,MLMNB),TSIN(NJK,MLMNB)  ,tsss(njk,mlmnb-1)
     &        , POL(NJK)       ,TOR(NJK)         ,TSC(lssl,njk)
     &        , stat=ierr14)
!      common /TRIGGS/ TCOS,TSIN,tsss,POL,TOR,TSS,TSC,MX,NX,ML,NL,MB,NB,
!     &                LMNV,LMNL,LMNB,MS,NS,QL,LMNS,lfrz,lfx,lsx,msx,nsx,
!     &                lss,cospar,cospr2,parity
C
C... In/Out Status of Common:
C     For FRV to FZV: Maybe Read, Not Written
C     For FVL: Maybe Read, Maybe Written
C     For FVLI to FVJAC: Maybe Read, Not Written
C     For FSIGMV to FPERPV: Maybe Read, Not Written
C     For FVGAM to fvq: Not Read, Not Written
C
C... FORMERLY Variables in Common Block /VMCFOU/ ...
      allocate( FRV(MLMNV,0:NI)   ,FZV(MLMNV,0:NI)   ,FVL(MLMNL,NI)
     &        , FVLI(MLMNL,0:NI)  ,FVJAC(MLMNB,0:NI) ,FSIGMV(MLMNV,0:NI)
     &        , FTAUV(MLMNV,0:NI),FPARPV(MLMNV,0:NI) ,FPERPV(MLMNV,0:NI)
     &        , FVGAM(MLMNB)     ,FVALF(MLMNB)       ,FVB(MLMNB)
     &        , fvq(mlmnb)       ,stat=ierr15)
!      common /VMCFOU/ FRV,FZV,FVL,FVLI,FVJAC,FVGAM,FVALF,FVB,fvq
C
C... In/Out Status of Common:
C     For RV to RVQ: Maybe Read, Maybe Written
C     For RVQT: Not Read, Maybe Written
C     For GTTV to PERPV: Maybe Read, Maybe Written
C
C... FORMERLY Variables in Common Block /VMCRSP/ ...
      allocate(    RV(NJK,0:NI),   ZV(NJK,0:NI),  RVT(NJK,0:NI)
     &        ,   ZVT(NJK,0:NI),  RVP(NJK,0:NI),  ZVP(NJK,0:NI)
     &        ,   RVQ(NJK,0:NI), RVQT(NJK,0:NI), GTTV(NJK,0:NI)
     &        ,  GTPV(NJK,0:NI), GPPV(NJK,0:NI), VJAC(NJK,0:NI)
     &        ,SIGMAV(NJK,0:NI),  VBT(NJK,0:NI),  VBP(NJK,0:NI)
     &        ,  VBSQ(NJK,0:NI), TAUV(NJK,0:NI),PARPV(NJK,0:NI)
     &        , PERPV(NJK,0:NI),   stat=ierr16)
!      common /VMCRSP/ RV,ZV,RVT,ZVT,RVP,ZVP,RVQ,RVQT,GTTV,GTPV,GPPV,
!     &                VJAC,VBT,VBP,VBSQ
C
      end subroutine TPRALL_BAP
