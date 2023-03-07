      subroutine mastersolv(initall)
      use particle
      use para
      use mpi_var
      implicit none

      real mu(ntot/nprocs),Ec(ntot/nprocs)
      real lambda(ntot/nprocs),initall(3,ntot/nprocs)
     

      if(ntot.lt.nprocs)then
         print*,'synchronize ntot,nprocs!'
         print*,'ntot,nprocs',ntot,nprocs
         stop
      endif
      
      Ec(:)=initall(1,:)
      lambda(:)=initall(2,:)
      mu(:)=initall(3,:)

      call solverRK(Ec,lambda,mu,CIv)
c      Integration de l'EDO

      


      end                       ! fin master
