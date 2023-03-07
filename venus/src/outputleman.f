      subroutine outputleman
      use leman
      use scatter
      use machine
      use mpi_var
      use futils
      implicit none
      
      real,allocatable :: SmP(:),ShP(:),PowHotL(:,:,:)
     $     ,kparL(:,:,:),kperpL(:,:,:)
     $     ,PowHot(:,:,:,:),factor(:)
      complex,allocatable :: EpL(:,:,:),EmL(:,:,:)
      integer,dimension(:),allocatable :: n_ant
      integer i,k,j,itp,d,nsLb
      real Eps,Ems,ShT,SmT,factot,Pwrite(2)
      character(10) file_name
      character(20) form
      logical filecheck
      integer fid,ints(4)
      

      r00=0.
      
      if(me2.eq.0)then
         d=0
         nph=0
         do while (d.eq.0)
            if(nph+1.lt.10)then
               write (file_name,"('fort.27_',i1.1)") nph+1
            else
               write (file_name,"('fort.27_',i2.1)") nph+1
            endif

            INQUIRE(FILE=file_name,EXIST=filecheck)
            if(filecheck.and.Powph(nph+1).gt.0.)then    
               nph=nph+1   
            else                ! for file fort.27_d
               if(filecheck)
     $              print*,'Omitting ',file_name,'and further files'
               d=1
            endif
         enddo
      endif

      call MPI_BCAST(nph,1,MPI_INTEGER,0,comm,ierr)
      allocate(SmP(nph),ShP(nph),n_ant(nph))
      if(sum(Powph).gt.size(Powph))Powph=1./nph
      

      if(me2.eq.0)then
         do d=1,nph  !loop over toroidal modes

            if(d.lt.10)then
               write (file_name,"('fort.27_',i1.1)") d
            else
               write (file_name,"('fort.27_',i2.1)") d
            endif


            call openf(file_name,fid,'r')
            call getarr(fid,'/ints',ints)
            nsL=ints(1)
            njL=ints(2)
            nkL=ints(3)
            n_ant(d)=ints(4)
           
            if(d.eq.1)allocate(EpL(nph,nsL,njL*nkL)
     $           ,EmL(nph,nsL,njL*nkL),PowHotL(nph,nsL,njL*nkL)
     $           ,kparL(nph,nsL,njL*nkL),kperpL(nph,nsL,njL*nkL))
            if(d.eq.1)allocate(cR_hp(nsL,njL*nkL),cZ_hp(nsL,njL*nkL))
            call getarr(fid,'/coordinateR',cR_hp)
            call getarr(fid,'/coordinateZ',cZ_hp)
            call getarr(fid,'/2D_Quants/Eplus',EpL(d,:,:))
            call getarr(fid,'/2D_Quants/Eminus',EmL(d,:,:))
            call getarr(fid,'/2D_Quants/kpar',kparL(d,:,:))
            call getarr(fid,'/2D_Quants/kperp',kperpL(d,:,:))
            call getarr(fid,'/Dep_Power',Pwrite)
            SmP(d)=Pwrite(1)
            ShP(d)=Pwrite(2)
            call closef(fid)
      enddo
      endif
       
      call MPI_BCAST(nsL,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(njL,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nkL,1,MPI_INTEGER,0,comm,ierr)
      if(me2.ne.0)then
         allocate(EpL(nph,nsL,njL*nkL),EmL(nph,nsL,njL*nkL)
     $        ,PowHotL(nph,nsL,njL*nkL)
     $        ,kparL(nph,nsL,njL*nkL),kperpL(nph,nsL,njL*nkL))
         allocate(cR_hp(nsL,njL*nkL),cZ_hp(nsL,njL*nkL))
      endif

      call MPI_BCAST(cR_hp,nsL*njL*nkL,MPI_DOUBLE_PRECISION
     $     ,0,comm,ierr)
      call MPI_BCAST(cZ_hp,nsL*njL*nkL,MPI_DOUBLE_PRECISION
     $     ,0,comm,ierr)
      call MPI_BCAST(EpL,nph*nsL*njL*nkL,MPI_DOUBLE_COMPLEX
     $     ,0,comm,ierr)
      call MPI_BCAST(EmL,nph*nsL*njL*nkL,MPI_DOUBLE_COMPLEX
     $     ,0,comm,ierr)
      call MPI_BCAST(PowHotL,nph*nsL*njL*nkL,MPI_DOUBLE_PRECISION
     $     ,0,comm,ierr)
      call MPI_BCAST(kparL,nph*nsL*njL*nkL,MPI_DOUBLE_PRECISION
     $     ,0,comm,ierr)
      call MPI_BCAST(kperpL,nph*nsL*njL*nkL,MPI_DOUBLE_PRECISION
     $     ,0,comm,ierr)
      call MPI_BCAST(SmP,nph,MPI_DOUBLE_PRECISION
     $     ,0,comm,ierr) 
      call MPI_BCAST(ShP,nph,MPI_DOUBLE_PRECISION
     $     ,0,comm,ierr)
      call MPI_BCAST(n_ant,nph,MPI_INTEGER
     $     ,0,comm,ierr)
      allocate(E(0:nsL,njL,nkL,2,nph)
     $     ,PowHot(0:nsL,njL,nkL,nph)
     $     ,R_hp(0:nsL,njL,nkL,nph)
     $     ,Z_hp(0:nsL,njL,nkL,nph)
     $     ,kpar(0:nsL,njL,nkL,nph)
     $     ,kperp(0:nsL,njL,nkL,nph)
     $     ,ntor(nph))
         

      ShT=sum(abs(ShP(:)))            ! total hot absorbed power
      SmT=sum(abs(SmP(:)))            ! total total deposited power
c         SmT=SmT+SmP(d)*SmP(d)

      do d=1,nph
         do k=1,nkL
            do j=1,njL
               itp=(k-1)*njL+j
               do i=1,nsL
                  E(i,j,k,1,d)=EpL(d,i,itp) ! Eplus
                  E(i,j,k,2,d)=EmL(d,i,itp) ! Eminus
                  PowHot(i,j,k,d)=PowHotL(d,i,itp) ! Power on hot ptcles
                  kpar(i,j,k,d)=kparL(d,i,itp) ! kpar 
                  kperp(i,j,k,d)=kperpL(d,i,itp) !kperp
                  R_hp(i,j,k,d)=cR_hp(i,itp) ! R coordinate
                  Z_hp(i,j,k,d)=cZ_hp(i,itp) ! Z coordinate
               enddo
               E(0,j,k,1,d)=1.5*E(1,j,k,1,d)-0.5*E(2,j,k,1,d)
               E(0,j,k,2,d)=1.5*E(1,j,k,2,d)-0.5*E(2,j,k,2,d)
               PowHot(0,j,k,d)=1.5*PowHot(1,j,k,d)-0.5*PowHot(2,j,k,d)
               kpar(0,j,k,d)=1.5*kpar(1,j,k,d)-0.5*kpar(2,j,k,d)
               kperp(0,j,k,d)=1.5*kperp(1,j,k,d)-0.5*kperp(2,j,k,d)
               R_hp(0,j,k,d)=1.5*R_hp(1,j,k,d)-0.5*R_hp(2,j,k,d)!r00
               Z_hp(0,j,k,d)=0.
               r00=r00+R_hp(0,j,k,d)
            enddo
         enddo
      enddo

      allocate(sL(0:nsL))
      if(me2.eq.0)then
         open(unit=25,file='grid.dat')
         read(25,'(i3)')nsLb
         if(nsL.ne.nsLb)then
            print*,'nsL,nsL grid',nsL,nsLb
            stop
         endif
         do i=0,nsL
            read(25,*)sL(i)
         enddo
      endif
      call MPI_BCAST(sL,nsL+1,MPI_DOUBLE_PRECISION
     $     ,0,comm,ierr)

      r00=r00/(njL*nkL*nph)
c      DepPow=DepPow*ShT/SmT ! input DepPow = total, including electrons and background
      
      allocate(factor(nph))
      factot=0.
      do d=1,nph
c         factor(d)=sqrt(abs(SmP(d))*DepPow/SmT) ! normalise each field to its deposited power, s.t. Power(n1)+Power(n2)=DepPow
         factor(d)=sqrt(DepPow/(nph*abs(SmP(d)))) ! equal power absorbed by each field
c         factor(d)=sqrt(DepPow/SmT) ! normalise to total dep. Power input as in Choi PoP 12,072505
c         factor(d)=sqrt(DepPow/ShT) ! normalise to hot dep. Power input as in Choi PoP 12,072505
c$$$c      factor=sqrt(DepPow*ShT)/SmT ! normalize to input Deposited Power, but only hot particles
! given power repartition in Powph
         factor(d)=sqrt(Powph(d)*DepPow/abs(SmP(d)))
         factot=factot+factor(d)
      enddo

      if(SELFO.gt.0)then
         ntor=n_ant
         write(form,"('(1x,a,',i2,'i5,a)')")nph
         if(me2.eq.0)write(*,form)'Forcing k||=',ntor,'/R'
         if(SELFO.gt.2)then
            E(:,:,:,1,:)=cmplx(Ecplus,0.)
            E(:,:,:,2,:)=cmplx(Ecminus,0.)
            factor=1.
            if(SELFO.gt.1)then
               if(me2.eq.0)write(*,'(a,f6.2)')' Forcing k_|_=',kperpc
            endif
         endif
      endif

      Eps=0.
      Ems=0.
      do d=1,nph
         do k=1,nkL
            do j=1,njL
               do i=1,nsL
                  E(i,j,k,1,d)=factor(d)*E(i,j,k,1,d)
                  E(i,j,k,2,d)=factor(d)*E(i,j,k,2,d)
                  Eps=Eps+abs(E(i,j,k,1,d))
                  Ems=Ems+abs(E(i,j,k,2,d))
               enddo
            enddo
         enddo
      enddo
      

      if(me2.eq.0.and.icrh.ne.0)then
         if(SELFO.eq.0)then
            write(form,"('(1x,a,',i2,'i5)')")nph
            write(*,form)'Toroidal modes',n_ant
            if(nph.gt.1)
c$$$     $           write(*,form)'E scaling [%] ',nint(factor/factot*100)
     $           write(*,form)'Power prop [%]',nint(Powph(1:nph)*100)
         endif
c         write(*,form)'P contrib [%] ',nint(ShP/SmT*100)
         write(*,'(a9,2es8.1E1,f6.2,a)')'<E+->,P '
     $        ,Eps/(nsL*njL*nkL*nph)!sqrt(real(nph)))
     $        ,Ems/(nsL*njL*nkL*nph)!sqrt(real(nph)))
     $        ,DepPow/1e6,'MW'
      endif

      deallocate(SmP,ShP,PowHotL,kparL,kperpL,PowHot,factor)
      deallocate(EpL,EmL,n_ant)

      end subroutine outputleman
