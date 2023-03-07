c----------------------------------------------------------------------c
C  17.03.06        LAST MODIFICATION 04.07.07     UHS: mtaskb DATA   D
c-----------------------------------------------------------------------
C
      subroutine mtaskb(tim)
C
      use tprb_param
      use physeq
      use radfct
      use triggs
      use paprpl
      use vmcfou
      use vmcrsp
      use boofou
      use boorsp
      use lambda
      use mappin
      use boomet
      use bderiv
      use comerc
      use merfou
      use balfou
!     use coboot
C
C.. Implicits ..
      implicit none
C
!      include 'parbal.inc'
C
C.. External Calls ..
      external BOPHYS, extint, metric, vmtobo, mercap, second
!      include 'tprcom.bal'
C
C.. Internal scalars
      integer :: i, l
      real :: tim(*), shalf
C ... Executable Statements ...
C
C
C     PERFORM MAPPING TO BOOZER COORDINATES.
c
c      do 10 i = 1,ni
      call vmtobo
     &                 (ni,njk,nj,lmnb,lmnl,nper,lmnb0,lvmtpr,tpi,dtdp,
     &                  dph,dth,fpp,ftp,cjv,civ,mb,nb,ml,nl,tsin,tcos,
     &                  vbp,vbt,vbsq,vjac,r,z,vl,phv,fvl,fr,fz,fphv,
     &                  fphvt,fphvp,fbsq,fbjac,fvalf,fvgam,fvb,fvq,
     &                  valft,vgamt,valf,vgam,vjtobj,thv,vq,sqgbsq,
     &                  sigmav,parpv,perpv,tauv,fsigmb,ftaub,fgparp,
     &                  fgperp)
C 10   continue
      call second (tim(5))
c
C     COMPUTE R/Z/PHI ON THE HALF INTEGER GRID.
c
c      do 20 l = 1,lmnb
      call extint
     &          (ni,njk,lmnb,mm,nmin,nmax,nper,lmetpr,rplmin,mb,nb,
     &             lfrz,tsin,tcos,ri,zi,phvi,fr,fz,fphv,fri,fzi,fphvi,
     &             nbalmn)
C 20   continue
c
c     calculate metric elements, jacobian and geometry (r/z/phi) on
c     the half integer grid in real boozer coordinate space.
c
c      do 30 i = 1,ni
      call metric
     &                 (ni,njk,lmnb,s,fpp,ftp,cjv,civ,mb,nb,tsin,tcos,r,
     &                  z,phv,rt,zt,rp,zp,rs,zs,rsq,ri,zi,phvi,bsq,bjac,
     &                  vjac,gttl,gtpl,gppl,gssl,gstl,fr,fz,fbsq,fphv,
     &                  phvp,phvt,fgparp,gparp,fgperp,gperp,
     &                  fsigmb,sigmab)
C 30   continue
      call second (tim(6))
c
      call bophys
     &                 (ni,njk,lmnb,lmnb0,mm,nmin,nmax,lcurrf,lvmtpr,
     &                  lpress,rplmin,djp,dtdp,s,ftp,fpp,ftpp,fppp,ci,
     &                  cj,cip,cjp,pp,wmag,vp,cipi,cjpi,pvpi,pth,ppi,
     &                  equi,civ,cjv,vvp,mb,nb,tsin,tcos,lfrz,gttl,gtpl,
     &                  gppl,gssu,bjac,bjacs,vjac,sigbs,bp,bt,bsq,
     &                  fbjac,fgparp,fsigbs,fphv,glm,bgradg,gparp,
     &                  parpav,parpvi,sigmab,taub,gperp,ftaub,fgperp,
     &                  fparjp,fpark,parjp,parkur,curfac)
      call second(tim(7))
!           print *,' lmnb0=',lmnb0,'   fr(lmnb0,ni)=',fr(lmnb0,ni)
C
      write (6, 889)
  889 format (/2x,"i",8x,"q(si)",6x,"mercier criterion: fluid and non-in
     &teracting hot particle models"/)
!
      if (lpress.ne.-9)
     &call mercap
     &                 (ni,njk,lmnb,lmnb0,nper,dtdp,s,ftp,fpp,ftpp,fppp,
     &                  ci,cj,pp,cjp,cip,mb,nb,tcos,tsin,gssu,gssl,gstl,
     &                  fbjac,bjac,bjacs,bsq,fr,fz,fbsq,fpark,fparjp,
     &                  fcurvs,fcurvn,flshea,fgssu,fbalcp,fbalcs,fbalcq,
     &                  fbaldp,fbalds,fbalnp,fbalns,fbalnq,fbaldn,
     &                  fbaldt,parkur,parjp,sigmab,taub,sigbs,gparp,
     &                  gperp,glm,bgradg,sigbxg,baldp,balds,balcp,balcs,
     &                  balcq,balnp,balns,balnq,baldn,curvds,curvdn,
     &                  shearl,vp)
      call second(tim(8))
c
C..OUTPUT IN UNIT 37 FOR GUIDING CENTRE CALCULATIONS
      write(37,820) ni, lmnb
      do i=1,ni
       do l=1,lmnb
        write(37,821) mb(l),nb(l),fr(l,i),fz(l,i),fphv(l,i),fbsq(l,i),
     &                fsigmb(l,i),ftaub(l,i),fgperp(l,i)-fgparp(l,i),
     &                fsigbs(l,i)
       enddo
      enddo
      do i=1,ni
       shalf=0.5*(s(i)+s(i-1))
       write(37,822)shalf,-ftp(i),-fpp(i),-ci(i),-cj(i),-cip(i),-cjp(i)
      enddo
 820  format(1x,2i6)
 821  format(1x,2i8,1p8e14.6)
 822  format(1x,1p7e14.6)
!.. BOOTSTRAP CURRENT CALCULATION
!     if (lbootj.eq.1) then
!     call bootsj(ni,njk,lmnb,npitch,dtdp,pitch,bnorm,fbjac,bjac,
!    &bsq,mb,nb,tcos,tsin,fc,ft,gb,g1mnin,g4mn,g1av,g4av,g2,g4,g2av,
!    &s,vp,ftp,fpp,ci,cj,pth,pp,ppi,tempp,tpi,boot,bjav,sqbbsq,
!    &dboot,rl31,rl32e,rl32i,jkmax,bsqav,bsqmax,djp,bdgradg,bdgradb,
!    &bgradbmn,bgradgmn,bdglnb,gbplat,rlampl,rnue)
!     end if
!
      return
      end subroutine mtaskb
