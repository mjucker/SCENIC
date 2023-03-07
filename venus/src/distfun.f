      subroutine distfun(d,ttot,vpar,vperp)
      
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

      integer ip,i,n,m,d,j,k
      real vpar(ntot/nprocs),vperp(ntot/nprocs),ttot
      real vperpmaxloc,vparmaxloc,intpol1,intpol2,norm,
     $     vparatbin(nbinsvpar),vperpatbin(nbinsvperp)
      integer vperpbin,vparbin
      integer nbinse,wndiag,linbin
      real distmatrix(nbinsvpar,nbinsvperp),
     $     distmatrixtot(nbinsvpar,nbinsvperp)
      real,allocatable,dimension(:) :: distmatrixE
     $     ,distmatrixEtot,distmatrixEth,distmatrixEthtot
     $     ,Eatbin,logEbin,intere
      real val(10),vals(8),pert(8),TpoTb(nsa)
      real theta,thetaper,phi,phiper,VV,VE
      real temptot(nbinsrad),radloc(ntot/nprocs)
      integer,allocatable,dimension(:) :: indxr(:),indxt(:),indxp(:)
     $     ,indxe(:)
      real Edmax,Eloc,dlogE
 
      linbin=0 ! linear or logarithmic energy bin

      allocate(indxr(ntot/nprocs),indxt(ntot/nprocs),indxp(ntot/nprocs)
     $     ,indxe(ntot/nprocs))
      nbinse=nbinsvpar*nbinsvperp
      allocate(distmatrixE(nbinse),distmatrixEtot(nbinse),Eatbin(nbinse)
     $     ,logEbin(nbinse),intere(0:nbinse))
      allocate(distmatrixEth(nbinse),distmatrixEthtot(nbinse))
      if(runctrl.eq.2.and.d.eq.0)then
         allocate(distmatrixstat(nbinsvpar,nbinsvperp))
         allocate(distmatrixEstat(nbinse))
         allocate(distmatrixEthstat(nbinse))
         distmatrixstat=0.
         distmatrixEstat=0.
         distmatrixEthstat=0.
      endif

      if(d.eq.-1)then    
         if(me2.eq.0)then
            open(unit=39,file='fort.39',form='unformatted')
            open(unit=59,file='fort.59',form='unformatted')
            open(unit=69,file='fort.69',form='unformatted')
            wndiag=ndiagnostics+2
            if(precond.eq.1)wndiag=ndiagnostics+1
            write(39)nbinsvpar,nbinsvperp
     $           ,nbinstor,nbinspolr(writerad),nbinsrad,
     $           wndiag,tfin
            write(69)nbinse,wndiag
            write(69)tfin,nall,totalweight
         endif
      endif
      

      pi=2.*acos(0.)
      twopi=2.*pi

! #################################################################
!            velocity distribution function: velocity bins
! #################################################################
      if(constbins*d.le.0)then   
         vperpmaxloc=0.
         vparmaxloc=0.
         do i=1,ntot/nprocs
            if(weight(i).ne.0.)then
               if(vperp(i).gt.vperpmaxloc)vperpmaxloc=vperp(i)
               if(abs(vpar(i)).gt.vparmaxloc)vparmaxloc=abs(vpar(i))
            endif
         enddo
         call MPI_ALLREDUCE(vperpmaxloc,vperpmax,1,
     $        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
         call MPI_ALLREDUCE(vparmaxloc,vparmax,1,
     $        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
         vperpmax=min(vperpmax,sqrt(2.*Tmax*Qel/mo))
         vparmax=min(vparmax,sqrt(2.*Tmax*Qel/mo))
         dvperp=vperpmax/nbinsvperp
         dvpar=2.*vparmax/nbinsvpar
      endif
      do i=1,nbinsvpar
         vparatbin(i)=-vparmax+real(i-0.5)*dvpar
      enddo
      do i=1,nbinsvperp
         vperpatbin(i)=real(i-0.5)*dvperp
      enddo
      do ip=1,ntot/nprocs
         vpar(ip)=vpar(ip)+vparmax
      enddo
! #################################################################
!           energy distribution function: bins
! #################################################################
!      Edmax=0.5*mo*(vparmax**2 + vperpmax**2)/Qel
      Edmax=Tmax
      if(linbin.eq.0)then
c$$$         logEbin(nbinse)=ceiling(log10(Edmax))
         logEbin(nbinse)=log10(Edmax)
         logEbin(1)=log10(Es/Qel) !max(log10(Es),logEbin(nbinse)-4) ! max 4 decades
         dlogE=(logEbin(nbinse)-logEbin(1))
     $        /(nbinse-1)
         Eatbin(1)=10.**logEbin(1)
         do i=2,nbinse
            logEbin(i)=logEbin(1) + dlogE*(i-1)
            Eatbin(i)=10.**logEbin(i)
         enddo
         intere(0)=0.
         do i=1,nbinse-1
            intere(i)=10.**(0.5*(logEbin(i+1)+logEbin(i)))
         enddo
         intere(nbinse)=10.**logEbin(nbinse)
      else
         intere(0)=Es/Qel
         do i=1,nbinse
            Eatbin(i)=Es/Qel+(i-0.5)*(Edmax-Es/Qel)/nbinse
            intere(i)=intere(0)+i*(Edmax-Es/Qel)/nbinse
         enddo
      endif
      if(me2.eq.0)then!.and.writerad*writepolr(writerad).ne.0)then
         if(d.eq.ndiagnostics.or.runctrl.ne.2)then
            write(39)-vparatbin ! due to negative Jacobian
            write(39)vperpatbin
            write(69)Eatbin
         endif
      endif
! #################################################################
!            distribution function: spatial and energy bins
! #################################################################
      do ip=1,ntot/nprocs
         do i=0,nbinsrad-1
c            if(runctrl.eq.2)then
c               radloc(ip)=savg(ip)
c            else
               radloc(ip)=CIv(1,ip)
c            endif
            if(radloc(ip).ge.intervs(i).and.radloc(ip).lt.intervs(i+1))
     $           indxr(ip)=i+1
         enddo
         if(radloc(ip).ge.intervs(nbinsrad))then
            indxr(ip)=nbinsrad!0
            weight(ip)=0.
         endif
         if(radloc(ip).lt.intervs(0))indxr(ip)=0
         thetaper=
     $        mod(mod(CIv(2,ip)-pi/nbinspolr(indxr(ip)),twopi)+twopi
     $        ,twopi)/twopi
c     $        mod(mod(CIv(2,ip),twopi)+twopi,twopi)/twopi
         indxp(ip)=floor(thetaper*nbinspolr(indxr(ip)))
         indxp(ip)=mod(indxp(ip)+(nbinspolr(indxr(ip))+1)/2
     $        ,nbinspolr(indxr(ip)))+1
         phiper=mod(mod(CIv(3,ip),twopiPer)+twopiPer,twopiPer)/twopiPer
         indxt(ip)=floor(phiper*nbinstor)+1
         
         Eloc=0.5*mo*((vpar(ip)-vparmax)**2. + vperp(ip)**2.)/Qel
         do i=0,nbinse-1
            if((Eloc.gt.intere(i)).and.(Eloc.le.intere(i+1)))then
               indxe(ip)=i+1
            endif
         enddo
         if(Eloc.ge.intere(nbinse))indxe(ip)=nbinse   
         if(Eloc.le.intere(0))indxe(ip)=1
      enddo
       
! #################################################################
!            create distribution functions
! #################################################################      

      vparbin=0
      vperpbin=0
      distmatrix=0.
      distmatrixE=0.
      distmatrixEth=0.
      do ip=1,ntot/nprocs
         
         i=indxr(ip)
         j=indxp(ip)
         k=indxt(ip)

         if(i.gt.-1)then

         if(i.eq.writerad.and.j.eq.writepolr(i).and.k.eq.writetor)then
            vperpbin=floor(vperp(ip)/dvperp)+1
            if(vperp(ip).ge.vperpmax)vperpbin=nbinsvperp
            vparbin=floor(vpar(ip)/dvpar)+1
            if(vpar(ip).ge.2.*vparmax)vparbin=nbinsvpar
            if(vpar(ip).le.0.)vparbin=1
            
            VV=vperpatbin(vperpbin)*V(k,j,i) ! total jacobian
            VV=VV*dvperp*dvpar  ! total phase space volume.
            
            distmatrix(vparbin,vperpbin)=
     $           distmatrix(vparbin,vperpbin)+weight(ip)/VV     
         endif
c         if(i.eq.writerad)then ! energy distribution
         VE=(intere(indxe(ip)) - intere(indxe(ip)-1))
         VE=VE*totalweight*sqrt(Eatbin(indxe(ip)))
         if(0.5*mo*((vpar(ip)-vparmax)**2. + vperp(ip)**2.)
     $        .lt.Ecrit(i))then
            distmatrixEth(indxe(ip))=distmatrixEth(indxe(ip))
     $           +weight(ip)/VE
         else
            distmatrixE(indxe(ip))=distmatrixE(indxe(ip))
     $           +weight(ip)/VE
         endif
c         endif ! energy distribution
         endif
      enddo

      distmatrixtot=0.
      call MPI_ALLREDUCE(distmatrix,distmatrixtot,
     $     nbinsvperp*nbinsvpar,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     comm,ierr)
      distmatrixEtot=0.
      call MPI_ALLREDUCE(distmatrixE,distmatrixEtot,
     $     nbinse,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     comm,ierr)
      distmatrixEthtot=0.
      call MPI_ALLREDUCE(distmatrixEth,distmatrixEthtot,
     $     nbinse,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     comm,ierr)
      
      if(runctrl.eq.2)then
         distmatrixstat=distmatrixstat+distmatrixtot
         distmatrixEstat=distmatrixEstat+distmatrixEtot
         distmatrixEthstat=distmatrixEthstat+distmatrixEthtot
      endif
      
      if(runctrl.ne.2.or.d.eq.ndiagnostics)then
         if(runctrl.eq.2)then
            distmatrixtot=distmatrixstat/(d+1)
            distmatrixEtot=distmatrixEstat/(d+1)
            distmatrixEthtot=distmatrixEthstat/(d+1)
         endif
         if(me2.eq.0)then
            write(39)distmatrixtot !f(v_par,v_perp,s,theta)
            write(69)distmatrixEtot !f(E)
            write(69)distmatrixEthtot !f(E)
         endif
      endif



      if(d.eq.ndiagnostics.and.runctrl.eq.2)
     $     deallocate(distmatrixstat,distmatrixEstat
     $     ,distmatrixEthstat)

      deallocate(indxr,indxt,indxe)
      deallocate(distmatrixE,distmatrixEtot,Eatbin,logEbin,intere)
      deallocate(distmatrixEth,distmatrixEthtot)


      end subroutine distfun


  
