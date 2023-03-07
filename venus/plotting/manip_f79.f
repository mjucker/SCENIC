      program manip_f79
      implicit none

      integer nstep,nskip
      real,allocatable,dimension(:) :: lambda,Ep,mu
      integer t,it
      real twopi,rnstep
      

      open(unit=78,file='data.79')
      twopi=8.*atan(1.)
c   Read initial stuff
      write(*,'("Number of time steps (..x.x.. and in thousands): ")',
     $     advance='no')
      read(*,710) rnstep
      write(*,'("Reading interval : ")',
     $     advance='no')
      read(*,711) nskip
 711  format(i8)
      nstep=nint(rnstep*1000)
      
      allocate(lambda(nstep),Ep(nstep),mu(nstep))

      write(78,'(a)')'lambda Energy mu'

      do t=1,nstep
         read(79,10)lambda(t),Ep(t),mu(t)
      enddo
      do t=1,nstep,nskip
         write(78,10)lambda(t),Ep(t),mu(t)
      enddo
      
 710  format(f9.5)
 10   format(1p3e14.6)
      end
c----------------------------------------------------------------------c
