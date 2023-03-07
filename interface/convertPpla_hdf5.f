      program convertPpla
      use futils
      implicit none

      integer ns,nj,nk,i,j,AP,N,ints(4),fid,n_ant
      double precision,allocatable :: s_hp(:),cR_hp(:,:),cZ_hp(:,:)
      double precision,allocatable :: ReE(:,:,:),ImE(:,:,:)
     $     ,ReP(:,:),ImP(:,:),PowHot(:,:)
     $     ,kpar(:,:),kperp(:,:)
     $     ,dummy(:,:),dummy3(:,:,:)
      double complex,allocatable :: Ep(:,:),Em(:,:)
      double precision ImPtot,RePtot,ImPhot,ImPth,Pwrite(2)
      double complex,allocatable ::  E1(:,:),E2(:,:)
      integer OutputSelect(8),NbDistr,AP_nb
      integer,allocatable :: k_lkp(:),Pres_Sp(:)


      open(unit=20,file='Ppla3D.dat',form='unformatted',action='read')!,access='direct')
      read(20) ns, nj, nk  ! # elements radial, poloidal,toroidal
      allocate(s_hp(ns),cR_hp(ns,nj*nk),cZ_hp(ns,nj*nk),dummy(ns,nj*nk))
      allocate(PowHot(ns,nj*nk),kpar(ns,nj*nk),kperp(ns,nj*nk))
      allocate(ReE(3,ns,nj*nk),ImE(3,ns,nj*nk),dummy3(3,ns,nj*nk))
c     $     ,ReP(ns,nj*nk),ImP(ns,nj*nk))
c      allocate(dummy(3,ns,nj*nk))
      read(20) s_hp             ! s half-mesh point
      read(20) cR_hp            ! R half-mesh point
      read(20) cZ_hp            ! Z half-mesh point
      read(20) OutputSelect     ! what does Ppla3D.dat contain?
      read(20) NbDistr          ! number of species considered for k_par
      allocate(k_lkp(NbDistr-1))
      k_lkp(1)=1
      do i=2,NbDistr-1
         read(20)k_lkp(i)
      enddo
      read(20) AP_nb            ! number of species for power absorption
      allocate(Pres_Sp(AP_nb))
      do i=1,AP_nb
         read(20) Pres_Sp(i)
      enddo

      read(20) RePtot           ! Puissance deposee: All species
      read(20) ImPtot
      if(Pres_Sp(AP_nb-1).eq.3)read(20)ImPth
      if(Pres_Sp(AP_nb).eq.4)then
         read(20) ImPhot        ! Dep. power hot species
      else
         print*,'NO HOT PARTICLES IN Ppla3D.dat'
         stop
      endif
     
      if(OutputSelect(1).eq.1)then
         read(20) ReE           ! Champ Electrique
         read(20) ImE
      else
         print*,'NO ELECTRIC FIELD IN Ppla3D.dat!'
         stop
      endif
 
      if(OutputSelect(2).eq.1)then
         read(20)dummy3         ! Mag. field
         read(20)dummy3
      endif
      if(OutputSelect(3).eq.1)then
         read(20)dummy3         ! Potential
         read(20)dummy3
      endif
      if(OutputSelect(4).eq.1)then
         read(20)dummy         ! Total Power dep
         read(20)dummy
         do i=2,AP_nb+1
            read(20)dummy
            read(20)PowHot
         enddo
         PowHot=abs(PowHot)
         print*,'PowHot/ImPhot',sum(PowHot)/ImPhot
      endif
      do i=1,sum(OutputSelect(5:7))
         read(20)dummy         ! dielectric tensor
         read(20)dummy
      enddo

      if(OutputSelect(8).eq.1)then
         read(20) n_ant
         do i=1,NbDistr
            read(20)kpar
         enddo
         read(20)kperp
      endif
         
      
      close(20)

      write(*,'(1x,a,f4.1,a,f4.1,a,f4.1,a)')
     $     'Absorption on minority: Thermal ',ImPth/ImPtot*100.
     $     ,'%  Hot: ',ImPhot/ImPtot*100.
     $     ,'%  Total: ',(ImPhot+ImPth)/ImPtot*100.,'%'


      if(ImPtot.ne.ImPtot.or.ImPtot*2.eq.ImPtot)then
         print*,'LEMAN OUTPUT NOT PHYSICAL'
         open(unit=99,file='LEMAN_ERROR')
         close(99)
      else

      allocate(Ep(ns,nj*nk),Em(ns,nj*nk))

c$$$      Ep(:,:)=(ReE(1,:,:)-ImE(2,:,:))*(ReE(1,:,:)-ImE(2,:,:)) + 
c$$$     $      (ImE(1,:,:)+ReE(2,:,:))*(ImE(1,:,:)+ReE(2,:,:))
c$$$      Ep(:,:)=sqrt(0.5*Ep(:,:))
c$$$      Em(:,:)=(ReE(1,:,:)+ImE(2,:,:))*(ReE(1,:,:)+ImE(2,:,:)) +
c$$$     $      (ImE(1,:,:)-ReE(2,:,:))*(ImE(1,:,:)-ReE(2,:,:))
c$$$      Em(:,:)=sqrt(0.5*Em(:,:))

      allocate(E1(ns,nj*nk),E2(ns,nj*nk))

      E1=cmplx(ReE(1,:,:),ImE(1,:,:))
      E2=cmplx(ReE(2,:,:),ImE(2,:,:))

      Ep=0.5*(E1 + cmplx(0,1)*E2)
      Em=0.5*(E1 - cmplx(0,1)*E2)


      ints(1)=ns
      ints(2)=nj
      ints(3)=nk
      ints(4)=n_ant
      Pwrite(1)=ImPtot
      Pwrite(2)=ImPhot+ImPth
      
      call creatf('fort.27',fid)
      call putarr(fid,'/ints',ints)
      call putarr(fid,'/coordinateR',cR_hp)
      call putarr(fid,'/coordinateZ',cZ_hp)
      call creatg(fid,'/2D_Quants')
      call putarr(fid,'/2D_Quants/Eplus',Ep)
      call putarr(fid,'/2D_Quants/Eminus',Em)
      call putarr(fid,'/2D_Quants/PowHot',PowHot)
      call putarr(fid,'/2D_Quants/kpar',kpar)
      call putarr(fid,'/2D_Quants/kperp',kperp)
      call putarr(fid,'/Dep_Power',Pwrite)
      call closef(fid)

      endif
 

      end program convertPpla
