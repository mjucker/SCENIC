      subroutine walltime
      implicit none
      include 'mpif.h'
      real, save :: stime
      real :: eltime, eltime_max
      integer, save :: icall
      integer :: me2, ierr
      integer, parameter :: fn = 2
      integer :: hours,mins,secs
!
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      eltime =  mpi_wtime()
!
      if (icall == 0) then      ! Initial call
         stime = eltime
         icall = 1
      end if
      eltime=eltime-stime
!     
      call mpi_comm_rank(MPI_COMM_WORLD, me2, ierr)
      call mpi_allreduce(eltime, eltime_max, 1, MPI_DOUBLE_PRECISION,
     $     MPI_MAX, MPI_COMM_WORLD, ierr)
      hours=floor(eltime_max/3600)
      mins=floor((eltime_max-hours*3600)/60.)
      secs=nint((eltime_max-hours*3600-mins*60))
      if(me2.eq.0.and.eltime_max.ne.0.) then
         write(*,'(/a,i4,a,i3.2,a,i3.2,a)')
     $        'CPU TIME ',hours,' h'
     $        ,mins,' m',secs,' s'
         open(fn,file='time.txt')
         write (fn, '(1x, 1e19.10)') eltime_max
         close(fn)
      end if
      
      end subroutine walltime
