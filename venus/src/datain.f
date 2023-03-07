      subroutine datain

      use particle
      use plasma
      use load
      use para
      use mpi_var
      use scatter
      use kicks
      use dstfun
      use leman
      use pival
      use terp
      implicit none

      integer j

      namelist /background/
     $     density,n0,n1,n2,n3,n4,n5,current,j1,j2,j3,j4,j5
     $     ,TiTe,mplasma,qplasma,Zeff,vcurr

      namelist /particule/
     $     Qel, mp, me, vc, epso,signelambda,
     $     Apart, Zpart, mo, Qe, sinit,thetainit,
     $     xiinit,tfin,muin, pitchangle, qp, Es, startpos,
     $     lengthpas,Tmax,navg,splitfact
      namelist /loading/
     $     Tbc,Tpc

      namelist /collisions/
     $     coulomb, icrh, Powph, DepPow, SELFO, Ecplus, Ecminus, kperpc

      namelist /diagnostics/
     $     ndiagnostics,nparts,nbinsvpar,nbinsvperp,
     $     nbinslam,nbinsen,
     $     nbinspol,nbinsrad,nbinstor,
     $     profexp,polexp,tens,constbins,precond,stat
     $     ,writerad,writepol,writetor
     $     ,autostop,dragcurr

      pi=2.*acos(0.)
      twopi=2.*pi
      twopiPer=twopi/Lper

c initialize values here if you want
      Qel=1.6022e-19      ! charge elementaire
      mp=1.6726e-27       ! masse proton
      me=9.1095e-31       ! masse electron
      vc=2.99792e+8       ! vitesse de la lumiere
      epso=8.8542e-12
      signelambda=1.      ! signe de v_par/v initial
      Apart=1             ! A particule (masse de la particule)
      Zpart=1             ! Z numero atomique
      density=3.6e+17
      n0=1.               ! background density expansion coefficients
      n1=0.
      n2=0.
      n3=0.
      n4=0.
      n5=0.
      current=1.
      j1=0.               ! if vcurr=1, toroidal ohmic current expansion coefficients
      j2=0.               ! note that j0=0. always!
      j3=0.
      j4=0.
      j5=0.
c      Tp=1.0e+03          ! Thermal particle temperature (Maxwellian)
      mplasma=1.           ! A thermal ions
      Zeff=-1.             ! plasma Zeff for drag current, locally computed without impurities if <0
      qp=0.0
      startpos=0          ! read initial positions from file: 0 no, 1 yes, re-write positions, 2 yes, do not re-write initial.out
      sinit=0.
      thetainit=0.
      xiinit=0.
      coulomb=1           ! 1 -> Boozer & Kuo-Petravic; 2 -> Strumberger
      icrh=1              ! ICRH kicks on/off
      Powph=2.            ! Power repartition. Max. 6 components for 6 modes
      SELFO=0             ! k||=n/R (>0), E=cst (2)
      Ecplus=1.           ! E+ if (SELFO>2
      Ecminus=1.          ! E- if (SELFO)>2
      kperpc=30.          ! kperp if SELFO>1
      lengthpas=0.        ! use variable time step if lengthpas=0.
      Tmax=1.e7           ! defines maximum energy for loading [eV]
      Es=0.               ! defines minimum energy for loading [eV]
      splitfact=5.        ! defines where to split into hot and thermal
      Tbc=1.              ! coefficient for loading Tpar
      Tpc=1.              ! coefficient for loading Tperp
      ndiagnostics=1      ! # of subintervals for diagnostics
      nparts=1            ! # of particles per cell
      nbinsvpar=14        ! # bins vpar
      nbinsvperp=7        ! # bins vperp
      nbinslam=50
      nbinsen=50
      nbinspol=35         ! # bins theta
      nbinsrad=48         ! # bins radial
      nbinstor=55         ! # bins toroidal
      profexp=1.          ! exponent for s bin profile
      polexp=1.           ! exponent for theta bin profile
      tens=1.e-5          ! tension for cubic splines
      constbins=0         ! wheter to use same v-bins for all ndiag
      vcurr=1             ! iota profile if =0, current profile if =1
      dragcurr=1          ! whether to correct current density with drag current
      precond=0           ! whether comment-out integration
      writerad=0          ! radial index to write fort.39
      writepol=0          ! poloidal index to write fort.39
      writetor=1          ! toroidal index to write fort.39
      autostop=0.         ! stop simulation at x slowing down times
      stat=0              ! Power deposition diagnostics
      

c read values from inputfile
      if(me2.eq.0)then
         read(5,background)
         read(5,particule)
         read(5,loading)
         read(5,collisions)
         read(5,diagnostics)
      endif
      
      call MPI_BCAST(lengthpas,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(Qel,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(mp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(me,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(vc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(epso,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(signelambda,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(Apart,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(Zpart,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(density,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(n1,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(n2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(n3,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(n4,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(n5,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(current,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(j1,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(j2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(j3,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(j4,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(j5,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(TiTe,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(mplasma,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(qplasma,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(Zeff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(qp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(sinit,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(muin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(pitchangle,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(tfin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(startpos,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(Tmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(Es,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(splitfact,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(Tbc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(Tpc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(navg,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(thetainit,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(xiinit,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(coulomb,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(icrh,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(Powph,6,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      call MPI_BCAST(DepPow,1,MPI_DOUBLE_PRECISION,0,comm,ierr)  
      call MPI_BCAST(SELFO,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(Ecplus,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(Ecminus,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      call MPI_BCAST(kperpc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)   
      call MPI_BCAST(ndiagnostics,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nparts,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(nbinsvpar,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nbinsvperp,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nbinslam,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nbinsen,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nbinspol,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nbinsrad,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nbinstor,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(profexp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(polexp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(tens,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(constbins,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(vcurr,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(dragcurr,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(precond,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(writerad,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(writepol,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(writetor,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(stat,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(autostop,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
        
      npitch=nprocs     ! one random number group / processor
      ne=nint(1.*nparts/nprocs)
      if(sinit.ne.0.)ne=1
      ntot=npitch*ne      ! always divideable by nprocs 

      if(me2.eq.0)then
         if(writerad.gt.nbinsrad)then
            print*,'WRITERAD .GT. NBINSRAD'
            stop
         elseif(writepol.gt.nbinspol)then
            print*,'WRITEPOL .GT. NBINSPOL'
            stop
         elseif(writetor.gt.nbinstor)then
            print*,'WRITETOR .GT. NBINSTOR'
            stop
         endif
      endif
      mplasma=mplasma*mp
      qplasma=qplasma*Qel
      Qe=Zpart*Qel     ! charge de la particule
      mo=Apart*mp            ! masse de la particule
      Es=Es*Qel
c      if(precond.eq.1)stat=2
      if(stat.gt.0)then
         allocate(Encollip(0:nbinsrad,2,ntot/nprocs))
         allocate(Encoll(0:nbinsrad,2),Encolltot(0:nbinsrad,2))
         allocate(Encollall(0:nbinsrad,nbinspol,nbinstor)
     $        ,Encollalltot(0:nbinsrad,nbinspol,nbinstor))
         allocate(lostE(nbinspol),lostEtot(nbinspol))
         Encollip=0.
         Encollall=0.
      else
         allocate(lostE(1))
      endif
      lostE=0.

      return
      end subroutine datain
