      subroutine findsigma(yn,vpar,vperp
     $     ,sigmadens,sigmappar,sigmapperp
     $     ,sigmadensang,sigmapparang,sigmapperpang)
      
      use particle
      use machine
      use dimarray
      use para
      use sinitial
      use dstfun
      use scatter
      use mpi_var
      use pival
      implicit none

      integer ip,i,j,k
      real yn(4,ntot/nprocs),vpar(ntot/nprocs),vperp(ntot/nprocs)
      real vperpmaxloc,vparmaxloc,intpol1,intpol2,norm,
     $     vparatbin(nbinsvpar),vperpatbin(nbinsvperp)
      integer vperpbin,vparbin
      real distmatrix(nbinsvpar,nbinsvperp),
     $     Nvperp(nbinsvperp),Nvpar(nbinsvpar),
     $     distmatrixtot(nbinsvpar,nbinsvperp)
      real pperp(nbinsrad),ppar(nbinsrad)
     $     ,dens(nbinsrad)
      real TpoTb(nsa),thetaper,phiper
      integer,allocatable,dimension(:) :: indxr(:),indxp(:),indxt(:)
      real sigmadens(nbinsrad),sigmadensang(nbinstor,nbinspol,nbinsrad)
     $     ,sigmappar(nbinsrad),sigmapparang(nbinstor,nbinspol,nbinsrad)
     $     ,sigmapperp(nbinsrad)
     $     ,sigmapperpang(nbinstor,nbinspol,nbinsrad)
      real temptot(nbinstor,nbinspol,nbinsrad)


      allocate(indxr(ntot/nprocs),indxp(ntot/nprocs),indxt(ntot/nprocs))
      

      do ip=1,ntot/nprocs
         do i=0,nbinsrad-1
            if(yn(1,ip).ge.intervs(i).and.yn(1,ip).lt.intervs(i+1))
     $           indxr(ip)=i+1
         enddo
         if(yn(1,ip).ge.intervs(nbinsrad))indxr(ip)=nbinsrad
         if(yn(1,ip).lt.intervs(0))indxr(ip)=0
         thetaper=
     $        mod(mod(yn(2,ip)-pi/nbinspolr(indxr(ip))
     $        ,twopi)+twopi,twopi)/twopi
         
         indxp(ip)=floor(thetaper*nbinspolr(indxr(ip)))
         indxp(ip)=mod(indxp(ip)+(nbinspolr(indxr(ip))+1)/2
     $        ,nbinspolr(indxr(ip)))+1
         phiper=mod(mod(yn(3,ip),twopiPer)+twopiPer,twopiPer)/twopiPer
         indxt(ip)=floor(phiper*nbinstor)+1
      enddo
  
      sigmadensang=0.
      sigmapparang=0.
      sigmapperpang=0.
      

      do ip=1,ntot/nprocs

         i=indxr(ip)
         j=indxp(ip)
         k=indxt(ip)
         

         if(i*weight(ip).gt.0)then
            sigmadensang(k,j,i)=sigmadensang(k,j,i)
     $           +(weight(ip)/V(k,j,i)-densang(k,j,i))**2.   
            sigmapparang(k,j,i)=sigmapparang(k,j,i)
     $           +(mo*vpar(ip)**2.*weight(ip)/V(k,j,i)
     $           -pparang(k,j,i))**2.
            sigmapperpang(k,j,i)=sigmapperpang(k,j,i)
     $           +(0.5*mo*vperp(ip)**2.*weight(ip)/V(k,j,i)
     $           -pperpang(k,j,i))**2.
         endif

      enddo
      temptot=0.
      call MPI_ALLREDUCE(sigmadensang,temptot,
     $     nbinstor*nbinspol*nbinsrad,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     comm,ierr)
      do i=1,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               sigmadensang(k,j,i)=temptot(k,j,i)
     $              /max(1,nmarksall(k,j,i))
               temptot(k,j,i)=0.
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(sigmapparang,temptot,
     $     nbinstor*nbinspol*nbinsrad,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     comm,ierr)
      do i=1,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               sigmapparang(k,j,i)=temptot(k,j,i)
     $              /max(1,nmarksall(k,j,i))
               temptot(k,j,i)=0.
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(sigmapperpang,temptot,
     $     nbinstor*nbinspol*nbinsrad,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     comm,ierr)
      do i=1,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               sigmapperpang(k,j,i)=temptot(k,j,i)
     $              /max(1,nmarksall(k,j,i))
            enddo
         enddo
      enddo
      sigmadens=0.
      sigmappar=0.
      sigmapperp=0.
      do i=1,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               sigmadens(i)=sigmadens(i)+sigmadensang(k,j,i)
               sigmappar(i)=sigmappar(i)+sigmapparang(k,j,i)
               sigmapperp(i)=sigmapperp(i)+sigmapperpang(k,j,i)
            enddo
         enddo
      enddo
      sigmadensang=sqrt(sigmadensang)
      sigmapparang=sqrt(sigmapparang)
      sigmapperpang=sqrt(sigmapperpang)
      sigmadens=sqrt(sigmadens/nbinspolr/nbinstor)
      sigmappar=sqrt(sigmappar/nbinspolr/nbinstor)
      sigmapperp=sqrt(sigmapperp/nbinspolr/nbinstor)
      
      deallocate(indxr,indxp,indxt)

      end subroutine findsigma
