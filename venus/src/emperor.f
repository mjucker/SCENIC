      subroutine emperor
      use particle
      use para
      use scatter
      use mpi_var
      use dstfun
      implicit none

      real,allocatable :: mu(:),Ec(:),lambda(:)
      integer ndiagtmp
      logical file_exists

C--------------------------------------------------------------
C      Initialisation
C--------------------------------------------------------------

      call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)


      call init_random

      open(unit=1,file='fort.37')
      
      call outputterp                 ! read output de terpsichore
      
      call outputleman                ! read output LEMan
      
      call array                      ! construction de tabarray(s),pertarray
     
      allocate(mu(ntot/nprocs),Ec(ntot/nprocs),lambda(ntot/nprocs))

      allocate(CIv(4,ntot/nprocs))
      
      if(me2.eq.0)write(*,'(a13,i5,es8.1E1)')
     $     ' nprocs,ntot ',nprocs,real(ntot)
c      call walltime
C-------------------------------------------------------------
C      Solver
C-------------------------------------------------------------
      call MPI_BARRIER(comm,ierr)
      runctrl=-1
      filename='state.vs'                        ! check if new start or 
      INQUIRE(FILE=filename, EXIST=file_exists)  ! restart from previous run
      if(file_exists)then
         runctrl=1              !runctrl+1
      else
         INQUIRE(FILE='loadarray.in',EXIST=file_exists)
         if(file_exists)runctrl=runctrl+1
      endif

c$$$      write(filename,'(a8,i4.4)')'state.vs',me2  ! only for backward compatibility reasons
               
         ndiagtmp=ndiagnostics
      do while(runctrl.lt.2)
         
         if(precond.eq.1)ndiagnostics=1
         call makeinputrandom(Ec,lambda,mu)
         
         if(runctrl.ge.0.and.precond.eq.0)then
            call solverRK(Ec,lambda,mu)
         endif
         
         call MPI_BARRIER(comm,ierr)
         
         icrh=0
         coulomb=0
         if(stat.ne.2)startpos=-2
         constbins=1
         ndiagnostics=ndiagtmp
         if(runctrl.ne.-1)then
            runctrl=2
c     ndiagnostics=0
            if(me2.eq.0)print*,' POSTRUN'
         else
c     ndiagnostics=0
            runctrl=-2
            ndiagnostics=0
            if(me2.eq.0)print*,' MAKE LOADARRAY'
         endif
         
         call solverRK(Ec,lambda,mu)

         if(me2.eq.0.and.runctrl.eq.-2)then
            print*,'loadarray.in created, stopping'
            stop
         endif
         call MPI_BARRIER(comm,ierr)
         runctrl=runctrl+1
      enddo

      end
