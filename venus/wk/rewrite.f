      program rewrite
      use mpi
      implicit none

      integer me2,ierr,status(MPI_STATUS_SIZE)
      integer comm,nprocs,npnt,tmp
      logical file_exists
      character*40 filename,endian 
      integer ntot,i,nparts,start,stp,k,ne,j
      real nall,wtmp,lostweight
      real,allocatable :: CIv(:,:),lambda(:),Ec(:),weight(:)
      real,allocatable :: Cv(:,:),lbd(:),E(:),wet(:)

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,me2,ierr)
      comm=MPI_COMM_WORLD
      call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
      
      if(me2.eq.0)then
         do i=1,10000
            write(filename,'(a8,i4.4)')'state.vs',i
            INQUIRE(FILE=filename,EXIST=file_exists)
            if(.not.file_exists)exit
         enddo
         ntot=i  
         write(*,'(a)')'Number of particles:'
         read*,nparts 
         write(*,'(a)')'Endian (b/l):'
         read*,endian
         write(*,'(a,i4,a,i4)')'Reducing ',ntot
     $        ,' files into ',nprocs
         if(endian.eq.'b')then
            write(*,'(a)')'Reading big endian files'
            endian='big_endian'
         else
            endian='little_endian'
         endif
      endif
      call MPI_BCAST(ntot,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nparts,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(endian,40,MPI_CHARACTER,0,comm,ierr)

      ne=nint(1.*nparts/ntot)
      if(ne*ntot.ne.nparts)then
         if(me2.eq.0)print*,'Real number of particles: ',ne*ntot
      endif
      nparts=ne*ntot

      allocate(Cv(4,nparts/ntot),lbd(nparts/ntot),E(nparts/ntot)
     $     ,wet(nparts/ntot))
      allocate(CIv(4,nparts/nprocs),lambda(nparts/nprocs)
     $     ,Ec(nparts/nprocs),weight(nparts/nprocs))

      npnt=nparts/ntot
      lostweight=0.
      k=0
      do i=me2,ntot-1,nprocs
         start=k*npnt+1
         stp=(k+1)*npnt
         write(filename,'(a8,i4.4)')'state.vs',i
         open(unit=11,file=filename,form='unformatted',convert=endian)
         read(11)wtmp
         lostweight=lostweight+wtmp
         read(11)nall
         read(11)Cv
         read(11)lbd
         read(11)E
         read(11)wet
         close(11,status='delete')
         do j=1,4
            CIv(j,start:stp)=Cv(j,:)
         enddo
         lambda(start:stp)=lbd
         Ec(start:stp)=E
         weight(start:stp)=wet
         k=k+1
      enddo   

      write(filename,'(a8,i4.4)')'state.vs',me2
      open(unit=12,file=filename,form='unformatted')
      write(12)lostweight
      write(12)nall
      write(12)CIv
      write(12)lambda
      write(12)Ec
      write(12)weight
      close(12)
         

      call mpi_finalize(ierr)

      end
