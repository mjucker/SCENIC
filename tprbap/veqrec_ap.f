C-----------------------------------------------------------------------
C  15.12.89        LAST MODIFICATION 20.03.06     UHS: veqrec DATA D
C-----------------------------------------------------------------------
C
      subroutine veqrec(tim)
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
      use colamn
      use lambda
C
C.. Implicits ..
      implicit none
C
C.. Local Scalars ..
      integer lmap
      real :: tim(*)
C
C.. External Calls ..
      external LGIKVM, cospol, mtaskl, vforce, second
!      include 'tprcom.inc'
C
C ... Executable Statements ...
C
C
C-----------------------------------------------------------------------
C
C     COSPOL: COMPUTES ANGLES AND TRIGONOMETRIC FUNCTIONS.
C
C-----------------------------------------------------------------------
      call cospol(nper,nj,njk,lmnv,lss,lmap,tpi,dth,dph,mx,nx,msx,nsx,
     &            tor,pol,tsin,tcos,tsc)
C-----------------------------------------------------------------------
C
C     LGIKVM: COMPUTES METRIC ELEMENTS IN EQUILIBRIUM SYSTEM,
C             COMPUTES LAMBDA AND, USING A STRAIGHTFORWARD DIS-
C                      CRETIZATION, EQUILIBRIUM QUANTITIES (S.A.).
C     SURPLO: PLOTS INPUT EQUILIBRIUM GRID
C
C-----------------------------------------------------------------------
C      do 10 i = 0,ni
      call LGIKVM(ni,nj,njk,mlmnv,nper,tpi,dph,mx,nx,tcos,tsin,frv,fzv,
     &            fvjac,fsigmv,ftauv,fparpv,fperpv,r,z,ri,zi,phvi,vjac,
     &            rv,zv,rvt,zvt,rvp,zvp,rvq,rvqt,sigmav,tauv,parpv,
     &            perpv,gttv,gtpv,gppv,mlmnb)
C 10   continue
      call second (tim(2))
C-----------------------------------------------------------------------
C
C     MTASKL: COMPUTES LAMBDA AND, USING A STRAIGHTFORWARD DIS-
C                      CRETIZATION, EQUILIBRIUM QUANTITIES (S.A.).
C             INTERPOLATES R/Z TO HALF INTEGER GRID.
C
C-----------------------------------------------------------------------
      call mtaskl(ni,njk,mlmnv,lmnl,lss,nbalmn,nper,nprocs,llampr,dtdp,
     &            s,ftp,fpp,civ,cjv,vvp,wmagv,mx,ml,nl,pol,tor,fr,fz,
     &            frv,fzv,fvl,fvli,tcos,tsin,tsss,r,z,zv,zvt,rvq,rvqt,
     &            vjac,vl,gttv,gtpv,gppv,vbp,vbt,vbsq,tsc,lsx,mlmnb,vlt,
     &            vlp,psivt,psivp,sigmav,parpv,parpav,pgppv,pgttv,pgtpv,
     &            al,gla,dna,dnb)
      call second (tim(3))
C
C-----------------------------------------------------------------------
C
C     VFORCE: COMPUTES MHD EQUILIBRIUM FORCE BALANCE IN VMEC COORDINATES
C
C-----------------------------------------------------------------------
C
C--------0---------0---------0---------0---------0---------0---------0-C
C
      call vforce(ni,njk,lmnb,nper,lnbaln,lnbalx,mm,lmxbal,nmin,nmax,
     &            lskewn,s,vvp,ftp,fpp,civ,cjv,civp,cjvp,pth,pvpi,equiv,
     &            pvp,mb,nb,pol,tor,tcos,tsin,lfrz,parpvi,parpav)
      end subroutine veqrec
