subroutine walltime
  implicit none
  include 'mpif.h'
  real*8, save :: stime
  real*8 :: eltime, eltime_max
  integer, save :: icall
  integer :: me, ierr
  integer, parameter :: fn = 2
!
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  eltime =  mpi_wtime()
!
  if (icall == 0) then  ! Initial call
     stime = eltime
     icall = 1
  end if
  eltime=eltime-stime
!
  call mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
  call mpi_reduce(eltime, eltime_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0,  &
       &  MPI_COMM_WORLD, ierr)
  if(me.eq.0) then
     write(*,'(/a,1pe14.6,a)') 'CPU TIME USED SO FAR =',eltime_max,' SECS'
     open(fn,file='time.txt')
     write (fn, '(1x, 1e19.10)') eltime_max
     close(fn)
  end if

end subroutine walltime
