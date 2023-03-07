C-----------------------------------------------------------------------
C  17.11.89        LAST MODIFICATION 28.02.06     UHS: mtaskl DATA   D
C-----------------------------------------------------------------------
C
      subroutine mtaskl(ni,njk,lmnv,lmnl,lss,nbalmn,nper,nprocs,llampr,
     &                  dtdp,s,ftp,fpp,civ,cjv,vvp,wmagv,mx,ml,nl,pol,
     &                  tor,fr,fz,frv,fzv,fvl,fvli,tcos,tsin,tsss,r,z,
     &                  zv,zvt,rvq,rvqt,vjac,vl,gttv,gtpv,gppv,vbp,vbt,
     &                  vbsq,tsc,lsx,lmnb,vlt,vlp,psivt,psivp,sigmav,
     &                  parpv,parpav,pgppv,pgttv,pgtpv,al,gla,dna,dnb)
C
C
C --- DIMENSIONING FOR INTEXT
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
      integer lmnv,lmnb,njk,lss,lmnl,ni,nper,nprocs,llampr,mx(*),ml(*),
     &        nl(*),lsx(-36:72,-72:72),nbalmn
      real dtdp,s(0:*),ftp(*),fpp(*),civ(*),cjv(*),vvp(*),wmagv(*),
     &     pol(*),tor(*),fr(lmnv,*),fz(lmnv,*),frv(lmnv,0:*),
     &     fzv(lmnv,0:*),fvl(lmnl,*),fvli(lmnl,0:*),tcos(njk,*),
     &     tsin(njk,*),tsss(njk,*),r(njk,0:*),z(njk,0:*),zv(njk,0:*),
     &     zvt(njk,0:*),rvq(njk,0:*),rvqt(njk,0:*),vjac(njk,0:*),
     &     vl(njk,*),gttv(njk,0:*),gtpv(njk,0:*),gppv(njk,0:*),
     &     vbp(njk,0:*),vbt(njk,0:*),vbsq(njk,0:*),tsc(lss,*),vlt(*),
     &     vlp(*),psivt(*),psivp(*),pgppv(*),pgttv(*)
      real, dimension(njk,0:ni) :: sigmav ,parpv
      real, dimension(ni)       :: parpav
      real pgtpv(*),al(lmnl,*),gla(3,*),dna(lmnl,2,*),dnb(lmnl,2,*)
C
C.. External Calls ..
      external intext, lamcal, trgfun
C
C ... Executable Statements ...
C
C
C     COMPUTE R/Z ON THE HALF INTEGER GRID.
C
C      do 10 l = 1,lmnv
      call intext(ni,njk,lmnv,nbalmn,mx,fr,fz,frv,fzv,tcos,tsin,r,z)
C 10   continue
C
C     PREPARE COEFFICIENTS FOR LAMBDA CALCULATION.
C     CALL LAMNEW TO DETERMINE FOURIER AMPLITUDES OF LAMBDA.
C     CALCULATE LAMBDA, COVARIANT B-FIELDS, B SQUARED AND
C     CURRENT FLUXES.
C      do 20 l = 1,lmnl
      call trgfun(njk,lmnl,ml,nl,pol,tor,tcos,tsin,tsss)
C 20   continue
C
C      do 30 i = 1,ni
C 30   continue
C
      call lamcal(ni,njk,lmnl,lss,nper,nprocs,llampr,dtdp,s,ftp,fpp,civ,
     &            cjv,vvp,wmagv,ml,nl,fvl,fvli,tcos,tsin,tsss,zv,zvt,
     &            rvq,rvqt,vjac,vl,gttv,gtpv,gppv,vbp,vbt,vbsq,tsc,lsx,
     &            lmnb,vlt,vlp,psivt,psivp,sigmav,parpv,parpav,
     &            pgppv,pgttv,pgtpv,al,gla,dna,dnb)
      end subroutine mtaskl
