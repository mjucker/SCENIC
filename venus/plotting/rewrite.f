      program rewrite
      use mpi
      implicit none

      integer me2,ierr,status(MPI_STATUS_SIZE)
      integer comm,nprocs
      logical file_exists
      character*40 filename   
      integer ntot,i,nparts,start,stp,k
      real nall
      real,allocatable :: CIv(:,:),lambda(:),Ec(:),weight(:)

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
         write(*,'(a,i4,a,i4)')'Reducing ',ntot
     $        ,' files into ',nprocs
      endif
      call MPI_BCAST(ntot,1,MPI_INTEGER,0,comm,ierr)
      call MPI_BCAST(nparts,1,MPI_INTEGER,0,comm,ierr)

      allocate(CIv(4,nparts/nprocs),lambda(nparts/nprocs)
     $     ,Ec(nparts/nprocs),weight(nparts/nprocs))

      k=0
      do i=me2,ntot-1,nprocs
         start=k*nparts/ntot+1
         stp=(k+1)*nparts/ntot
         write(filename,'(a8,i4.4)')'state.vs',i
         open(unit=11,file=filename,form='unformatted')
         read(11)nall
         read(11)CIv(:,start:stp)
         read(11)lambda(start:stp)
         read(11)Ec(start:stp)
         read(11)weight(start:stp)
         close(11,status='delete')
         if(i.eq.me2)then
            write(filename,'(a8,i4.4)')'state.vs',me2
            open(unit=12,file=filename,form='unformatted')
            write(12)nall
         endif         
         k=k+1
      enddo   

      write(12)CIv
      write(12)lambda
      write(12)Ec
      write(12)weight
      close(12)
         

      call mpi_finalize(ierr)

      end
