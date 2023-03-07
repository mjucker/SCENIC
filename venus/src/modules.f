
      module particle
      real Qel, mp, me, vc, muin,dels, 
     >     signelambda,Apart,Zpart,Zeff,mo,Qe,
     $     sinit,tfin,pitchangle,qp,Es,gamma,Tmax
      real smin,smax,splitfact
      integer startpos,runctrl,precond
      real thetainit,xiinit,lostweight!totalweight
      real plasmaV
      real,allocatable :: CIv(:,:),savg(:),Ecrit(:)
      end module particle

      module load
      real Tbc,Tpc
      end module load

      module pival
      real pi,twopi,twopiPer
      end module pival

      module para
      integer ::
     $     npitch,!=1   ,! # of random number generator groups for compatibility
     $     ne,!=1      ,  ! # to complete total number of particles
     $     ntot !=npitch*ne, ! total # of particles
      real lengthpas
      end module para

      module dimarray
      integer,parameter :: nxia= 72,ntha=1000
      integer :: nsa
      end module dimarray

      module terp
      integer,parameter :: Lper=2
      integer :: ni,lmnb
      end module terp

      module terpall
      use terp
      integer,allocatable :: mb(:),nm(:)
      real,allocatable,dimension(:) :: si,psi,chi,psip,chip,
     >     ci,cj,cip,cjp
      real,allocatable,dimension(:,:) :: bsmn,sigmn,taumn,dpdsmn,
     $     sigBsmn
      end module terpall

      module plasma
      real epso,qplasma
      real mplasma
      real density,current,Te,Ti,TiTe
      real n0,n1,n2,n3,n4,n5
      real j1,j2,j3,j4,j5
      integer vcurr
      end module plasma

      module machine
      real r00,a00,b00,Bc
      end module machine

      module mpi_var
      use mpi
      integer me2,ierr,status(MPI_STATUS_SIZE)
      integer comm,nprocs
      character*40 filename      
c      integer ipstart,ipstop
      end module mpi_var

      module arrays
      use dimarray
      REAL,allocatable :: tabarray(:,:,:,:)
      REAL,allocatable :: tabarrays(:,:)
c      REAL,allocatable :: pertarray(:,:,:,:)
      end module arrays
      
      module check
      use para
      real thin,phin
      real thend,phend
      integer tin,tend
      real yn2b4
      end module check

      module interp
      real val(10),vals(8),pert(8)
      end module interp

      module sinitial
      integer sinind
      real starts,ds,ds2
      end module sinitial
      
      module scatter
      integer coulomb,icrh,nph,stat,SELFO
      integer,allocatable :: ntor(:)
      real,dimension(:),allocatable ::  pp,weight
      double complex Eplus,Eminus
      real kperpend,kparallel,omega,Ecplus,Ecminus,kperpc
      real kparavg,kperpavg,kparavgtot,kperpavgtot,Ravg,Ravgtot
      real,allocatable :: Enkicks(:,:,:),Enkickstot(:,:,:),RZ(:,:,:,:)
     $     ,Encollip(:,:,:),Enkicksip(:,:,:,:),Enkicksradip(:,:,:)
     $     ,Encoll(:,:),Encolltot(:,:)
     $     ,lostE(:),lostEtot(:)
     $     ,Encollall(:,:,:),Encollalltot(:,:,:)
     $     ,Enkicksrad(:,:),Enkicksradtot(:,:)
      real RZbox(4),DepPow,Powph(6)
      real,allocatable :: Efact(:)
      end module scatter

      module kicks
      use para
      use mpi_var
      integer harmonics
      integer(8) nkicks,nkicklocal,nkickstot
      parameter(harmonics=1)
      real,allocatable ::  ifarg(:,:,:),nuDDbef(:,:,:),dtbef(:)
      real airyf,kickfactor,kickfactot
      real kickfactor2,kickfactot2
      end module kicks

      module rand_var
      use para
      real random
      integer,allocatable :: ran_index(:)
      double precision,allocatable :: ran_array(:,:)
      external random
      character*(34) random_seed_str
      end module rand_var

      module dstfun
      integer :: constbins,dragcurr
      integer :: nbinsvpar,nbinsvperp
      integer :: nbinspol,nbinsrad,nbinstor,ndiagnostics
      integer :: nbinslam,nbinsen
      integer :: writerad,writepol,writetor
      integer,allocatable :: nbinspolr(:),writepolr(:)
     $     ,nmarks(:,:,:),nmarksloc(:,:,:),nmarksall(:,:,:)
      real autostop
      real nparts,distfn,navg,totnavg,nall,tens
      real vperpmax,vparmax,dvperp,dvpar,Eavg,profexp,polexp
      real corr,lossall
      real,allocatable :: sl(:),s(:),intervs(:),Arad(:),trap(:)
      real,allocatable,dimension(:,:,:) :: V,BpB
     $     ,densang,densthang,pparang,pperpang
      real,allocatable :: distmatrixstat(:,:),distmatrixEstat(:)
     $     ,distmatrixEthstat(:)
      real,allocatable :: currt(:),currp(:),currtot(:),dsb(:),tensall(:)
      real,allocatable :: densth(:),pth(:)
      real :: totweight,totalweight
      real,allocatable :: dens_s(:),anis_s(:),Tpar_s(:),densi(:,:,:)

      end module dstfun

      module leman
      integer nsL,njL,nkL
      real,allocatable :: cR_hp(:,:),R_hp(:,:,:,:)
     $     ,cZ_hp(:,:),Z_hp(:,:,:,:)
     $     ,kpar(:,:,:,:),kperp(:,:,:,:)
     $     ,sL(:)
      double complex,allocatable :: E(:,:,:,:,:)
      end module leman

      module lkup
      real,allocatable :: loadarray(:,:,:,:)
      real Emax,Emin
      integer nrad,npol,nla,nen
      end module lkup
      
