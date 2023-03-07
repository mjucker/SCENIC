      program main
      use para
      use mpi_var
      implicit none
c      external datain,empero
      
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,me2,ierr)
      comm=MPI_COMM_WORLD
      call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)

      call walltime

      call datain

      if(ntot.lt.nprocs)then
         print*,'synchronize ntot,nprocs!'
         print*,'ntot,nprocs',ntot,nprocs
         stop
      endif

      call emperor  

      call MPI_BARRIER(comm,ierr)
      if(me2.eq.0)write(6,*) ' TERMINATED NORMALLY'

      call walltime
      call mpi_finalize(ierr)
      end program main
