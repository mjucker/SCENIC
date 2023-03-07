      subroutine moments(d,ttot,vpar,vperp)
      
      use particle
      use machine
      use plasma
      use dimarray
      use para
      use sinitial
      use dstfun
      use scatter
      use mpi_var
      use pival
      implicit none

      integer ip,i,d,j,k,ndummy,angsum
      real ttot,vip2,pos(4)
      real intpol1,intpol2,norm,Eavgtot
     $     ,denssum
      real,allocatable ::  pperp(:),ppar(:),dens(:),currpp(:)
      real vpar(ntot/nprocs),vperp(ntot/nprocs)
c      real,intent(inout)::vpar(*)
      real val(10),vals(8),pert(8),TpoTb(nbinsrad),TpoTb_d(nsa)
      real theta,thetaper,phi,phiper
      integer,allocatable :: indxr(:),indxp(:),indxt(:)
      real sigmadens(nbinsrad),sigmappar(nbinsrad),sigmapperp(nbinsrad)
     $     ,sigmaanis(nbinsrad)
     $     ,sigmadensang(nbinstor,nbinspol,nbinsrad)
     $     ,sigmapparang(nbinstor,nbinspol,nbinsrad)
     $     ,sigmapperpang(nbinstor,nbinspol,nbinsrad)
      real,allocatable ::  temptot(:,:,:),ssum(:)
      real dummy,radloc(ntot/nprocs)
      real nud,tsi,tse,dnuedee,dnuedei,val1(4)
      double complex ef(2)
      logical file

      allocate(indxr(ntot/nprocs),indxt(ntot/nprocs),indxp(ntot/nprocs)) 
      allocate(pperp(0:nbinsrad),ppar(0:nbinsrad)
     $     ,dens(0:nbinsrad),currpp(0:nbinsrad)
     $     ,temptot(nbinstor,nbinspol,0:nbinsrad)
     $     ,ssum(0:nbinsrad))


      if(runctrl.eq.2.and.d.eq.ndiagnostics.and.me2.eq.0)then
         print*,' '
         print*,' *** statistics ***'
      endif

      if(d.eq.-1)call initialise

  
      if(d.le.0.or.runctrl.ne.2)then
         do i=0,nbinsrad
            sl(i)=0.
            do j=1,nbinspol
               do k=1,nbinstor
                  nmarksloc(k,j,i)=0
               enddo
            enddo
         enddo
      endif

      totweight=0.
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
     $        mod(mod(CIv(2,ip)-pi,twopi)+twopi
     $        ,twopi)/twopi
c     $        mod(mod(CIv(2,ip),twopi)+twopi,twopi)/twopi
         indxp(ip)=floor(thetaper*nbinspolr(indxr(ip)))+1
c         indxp(ip)=mod(indxp(ip)+(nbinspolr(indxr(ip))+1)/2
c     $        ,nbinspolr(indxr(ip)))+1
         phiper=mod(mod(CIv(3,ip),twopiPer)+twopiPer,twopiPer)/twopiPer
         indxt(ip)=floor(phiper*nbinstor)+1
         if(indxt(ip).eq.nbinstor+1)indxt(ip)=1
c   find <s> of every bin 
         if(weight(ip).ne.0.)then
            sl(indxr(ip))=sl(indxr(ip))+radloc(ip)
            nmarksloc(indxt(ip),indxp(ip),indxr(ip))
     $           =nmarksloc(indxt(ip),indxp(ip),indxr(ip))+1
            if(indxr(ip).gt.0)totweight=totweight+weight(ip)
         endif
      enddo

      call MPI_REDUCE(lostweight,lossall,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,0,comm,ierr)
c$$$      call MPI_ALLREDUCE(totweight,totalweight,1
c$$$     $     ,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
c$$$      dummy=totalweight
      call MPI_ALLREDUCE(sum(weight),totalweight,1
     $     ,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
c$$$      corr=totalweight/(totalweight+lossall)
      corr=1.
c      print*,dummy/totalweight


      if(d.eq.ndiagnostics.or.runctrl.ne.2)then      
         do i=0,nbinsrad
            do j=1,nbinspol
               do k=1,nbinstor
                  nmarksall(k,j,i)=0
               enddo
            enddo
         enddo
         call MPI_ALLREDUCE(nmarksloc,nmarksall,
     $        nbinstor*nbinspol*(nbinsrad+1),MPI_INTEGER,MPI_SUM,
     $        comm,ierr)
                  
         call MPI_ALLREDUCE(sl,ssum,nbinsrad+1,MPI_DOUBLE_PRECISION
     $        ,MPI_SUM,comm,ierr)
         do i=0,nbinsrad
            angsum=0
            do j=1,nbinspol
               do k=1,nbinstor
                  angsum=angsum+nmarksall(k,j,i)
               enddo
            enddo
            if(angsum.gt.0)then
               s(i)=ssum(i)/angsum
            elseif(i.gt.0)then
c            if(i.gt.0)then
               s(i)=0.5*(intervs(i-1)+intervs(i))
            else
               s(i)=0.5*intervs(0)
            endif
         enddo
      endif

      call volume

      if(d.le.0.or.runctrl.lt.2)then
         Eavg=0.

         do i=0,nbinsrad
            do j=1,nbinspol
               do k=1,nbinstor
                  nmarks(k,j,i)=0
                  densang(k,j,i)=0.
                  pparang(k,j,i)=0.
                  pperpang(k,j,i)=0.
                  densthang(k,j,i)=0.
               enddo
            enddo
            currt(i)=0.
            currp(i)=0.
            if(vcurr.eq.1)currtot(i)=0.
            densth(i)=0.
            pth(i)=0.
         enddo
      endif
!***************************************************************
! INTEGRATION OF MOMENTS      
!***************************************************************
      do ip=1,ntot/nprocs

         i=indxr(ip)
         j=indxp(ip)
         k=indxt(ip)
        
         if(i.gt.-1)then
            
            vip2=vpar(ip)**2.+vperp(ip)**2.

            Eavg=Eavg+vip2*weight(ip)

            if(0.5*mo*vip2.lt.Ecrit(i))then
               densthang(k,j,i)=densthang(k,j,i)+weight(ip)/V(k,j,i)
               densth(i)=densth(i)+weight(ip)/V(k,j,i)
               pth(i)=pth(i)+vip2*weight(ip)/V(k,j,i)
            else
               
               if(weight(ip).ne.0.)nmarks(k,j,i)=nmarks(k,j,i)+1

               densang(k,j,i)=densang(k,j,i) + weight(ip)/V(k,j,i)

               pparang(k,j,i)=pparang(k,j,i)
     $              + vpar(ip)*vpar(ip)*weight(ip)/V(k,j,i)
               pperpang(k,j,i)=pperpang(k,j,i)
     $              + vperp(ip)*vperp(ip)*weight(ip)/V(k,j,i) 
            endif
           
            if(runctrl.eq.2)then
               if(trap(ip).gt.0)then ! trapped particles
                  currt(i)=currt(i)
     $                 + vpar(ip)*BpB(k,j,i)*weight(ip)/V(k,j,i)
c     $              + vpar(ip)*weight(ip)/V(j,i)
               else
                  currp(i)=currp(i) ! passing particles
     $                 + vpar(ip)*BpB(k,j,i)*weight(ip)/V(k,j,i)
c     $              + vpar(ip)*weight(ip)/V(j,i)
               endif
            else
               currt(i)=currt(i)
     $              + vpar(ip)*BpB(k,j,i)*weight(ip)/V(k,j,i)
            endif

         endif
      enddo
! END INTEGRATION PER PROCESSOR  
!--------------------------------------------------------
! Complete moment normalisations
!--------------------------------------------------------
    
      if(d.eq.ndiagnostics.or.runctrl.ne.2)then

! average particle energy
      call MPI_ALLREDUCE(Eavg,Eavgtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     $           comm,ierr)
      Eavgtot=Eavgtot/Qel  ![eV]
      Eavg=0.5*mo*Eavgtot/totalweight

      if(d.eq.-1.and.me2.eq.0)then
         pos=0.
         pos(1)=0.5/nsa
         call updatenu(Eavg*Qel,pos,nud,tsi,tse
     $        ,dnuedee,dnuedei)
         write(*,'(a12,es8.1E1,a,es8.1E1,a,es8.1E1,a)')
     $        ' tfin,ts,<E>',tfin,'s',tse,'s',Eavg,'eV'
      endif
 
! number of hot particle markers
      do i=0,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               temptot(k,j,i)=0.
               nmarksall(k,j,i)=0
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(nmarks,nmarksall,
     $     nbinstor*nbinspol*(nbinsrad+1),MPI_INTEGER,MPI_SUM,
     $     comm,ierr)

! hot particle density 
      call MPI_ALLREDUCE(densang,temptot,
     $     nbinstor*nbinspol*(nbinsrad+1)
     $     ,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
c      denspol=navg*twopi*temptot/totalweight
      do i=0,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               densang(k,j,i)=nall*temptot(k,j,i)/totalweight
               temptot(k,j,i)=0.
            enddo
         enddo
      enddo

! hot particle pressures 
c      temptot=0.
      call MPI_ALLREDUCE(pparang,temptot,
     $     nbinstor*nbinspol*(nbinsrad+1)
     $     ,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
c      pparpol=navg*twopi*mo*temptot/totalweight
      do i=0,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               pparang(k,j,i)=nall*mo*temptot(k,j,i)/totalweight
               temptot(k,j,i)=0.
            enddo
         enddo
      enddo

c      temptot=0.  
      call MPI_ALLREDUCE(pperpang,temptot,
     $     nbinstor*nbinspol*(nbinsrad+1),MPI_DOUBLE_PRECISION,MPI_SUM,
     $     comm,ierr)
c      pperppol=navg*twopi*0.5*mo*temptot/totalweight
      do i=0,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               pperpang(k,j,i)=nall*0.5*mo*temptot(k,j,i)/totalweight
               temptot(k,j,i)=0.
            enddo
         enddo
      enddo
            
! hot particle currents
      do i=0,nbinsrad
         currpp(i)=0.
      enddo
      call MPI_ALLREDUCE(currt,currpp,
     $           nbinsrad+1,MPI_DOUBLE_PRECISION,MPI_SUM,
     $           comm,ierr)
c      currpol=Qe*navg*twopi*mo*temptot/totalweight
      do i=0,nbinsrad
         currt(i)=nall*Qe*currpp(i)/totalweight
     $        /nbinspolr(i)/nbinstor
         currpp(i)=0.
      enddo
c      currt=currt/twopi ! flux surface average
c      currt=currt/nbinspolr ! flux surface average

c      currpp=0.
      call MPI_ALLREDUCE(currp,currpp,
     $           nbinsrad+1,MPI_DOUBLE_PRECISION,MPI_SUM,
     $           comm,ierr)
c      currpol=Qe*navg*twopi*mo*temptot/totalweight
      do i=0,nbinsrad
         currp(i)=nall*Qe*currpp(i)/totalweight
     $        /nbinspolr(i)/nbinstor
         currpp(i)=0.
      enddo
c      currp=currp/twopi ! flux surface average
c      currp=currp/nbinspolr ! flux surface average

! minority thermal density
c      temptot=0.
      call MPI_ALLREDUCE(densthang,temptot,
     $     nbinstor*nbinspol*(nbinsrad+1),MPI_DOUBLE_PRECISION,MPI_SUM,
     $     comm,ierr)
c      pperppol=navg*twopi*0.5*mo*temptot/totalweight
      do i=0,nbinsrad
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               densthang(k,j,i)=nall*temptot(k,j,i)/totalweight
               temptot(k,j,i)=0.
            enddo
         enddo
      enddo

c      currpp=0.
      call MPI_ALLREDUCE(densth,currpp,
     $           nbinsrad+1,MPI_DOUBLE_PRECISION,MPI_SUM,
     $           comm,ierr)
      do i=0,nbinsrad
         densth(i)=nall*currpp(i)/totalweight/nbinspolr(i)/nbinstor ! flux surface average
         currpp(i)=0.
      enddo

! minority thremal pressure
c      currpp=0.
      call MPI_ALLREDUCE(pth,currpp,
     $           nbinsrad+1,MPI_DOUBLE_PRECISION,MPI_SUM,
     $           comm,ierr)
      do i=0,nbinsrad
         pth(i)=nall*mo*currpp(i)/totalweight/nbinspolr(i)/nbinstor/3.
      enddo
c      pth=pth/nbinspolr/3. ! flux surface average, p=(ppar+2pperp)/3
!*****************************************************************
! END INTEGRATION
!*****************************************************************
      if(runctrl.eq.2.and.d.gt.0)then
         Eavg=Eavg/(d+1)
         do i=0,nbinsrad
            do j=1,nbinspolr(i)
               do k=1,nbinstor
                  densang(k,j,i)=densang(k,j,i)/(d+1)
                  pparang(k,j,i)=pparang(k,j,i)/(d+1)
                  pperpang(k,j,i)=pperpang(k,j,i)/(d+1)
                  densthang(k,j,i)=densthang(k,j,i)/(d+1)
                  nmarksall(k,j,i)=nmarksall(k,j,i)/(d+1)
               enddo
            enddo
            currt(i)=currt(i)/(d+1)
            currp(i)=currp(i)/(d+1)
            if(vcurr.eq.1)currtot(i)=currt(i)+currp(i)
            densth(i)=densth(i)/(d+1)
            pth(i)=pth(i)/(d+1)
         enddo
      endif
      denssum=0.
      do i=0,nbinsrad
         pperp(i)=0.
         ppar(i)=0.
         dens(i)=0.
         do j=1,nbinspolr(i)
            theta=-pi+twopi*(j-0.5)/nbinspolr(i)
            do k=1,nbinstor
               ppar(i)=ppar(i)+pparang(k,j,i)
               pperp(i)=pperp(i)+pperpang(k,j,i)
               dens(i)=dens(i)+densang(k,j,i)
               if(d.eq.ndiagnostics.and.runctrl.eq.2)then
                  phi=twopiPer*(k-0.5)/nbinstor
                  ndummy=1
                  call interpE(ndummy,s(i),theta,phi,val1,ef)
                  write(13,'(7es11.3,i9)')val1(1),val1(2),phi
     $                 ,pparang(k,j,i),pperpang(k,j,i)
     $                 ,densang(k,j,i),densthang(k,j,i)
     $                 ,nmarksall(k,j,i)
               endif
            enddo
         enddo
         dens(i)=dens(i)/nbinspolr(i)/nbinstor
         ppar(i)=ppar(i)/nbinspolr(i)/nbinstor
         pperp(i)=pperp(i)/nbinspolr(i)/nbinstor
         denssum=denssum+dens(i)+densth(i)
      enddo
c$$$      corr=navg/denssum*nbinsrad
      
 
      call findsigma(CIv,vpar,vperp
     $     ,sigmadens,sigmappar,sigmapperp
     $     ,sigmadensang,sigmapparang,sigmapperpang)
      

      if(me2.eq.0)then
         if(nbinsvpar*nbinsvperp.ne.0)then
            write(59)ttot
            write(59)s!(1:nbinsrad)
            write(59)dens
            write(59)sigmadens
            write(59)densth
            write(59)ppar
            write(59)sigmappar
            write(59)pperp      !n(s),p_par(s),p_perp(s)
            write(59)sigmapperp
            write(59)currt
            write(59)currp
         endif
         
         if(d.eq.ndiagnostics.and.runctrl.eq.2)then
            call inputvmec(sigmapparang,sigmapperpang,sigmaanis
     $           ,TpoTb,TpoTb_d)
            call inputleman(TpoTb,TpoTb_d,sigmadens,sigmaanis)
         endif
      endif    


      if(me2.eq.0.and.d.eq.-1.and.runctrl.ge.0)
     $     write(*,'(a,i3,a,i3)')' nbinspolr'
     $     ,nbinspolr(1),' ->',nbinspolr(nbinsrad)
      pos=0.
      pos(1)=0.5/nsa
      call updatenu(Eavg*Qel,pos,nud,tsi,tse
     $     ,dnuedee,dnuedei)
      if(me2.eq.0.and.runctrl.ne.2.and.runctrl*d.ne.-1)
     $     write(79,'(4es14.6)')ttot
     $     ,tsi*tse/(tsi+tse),1./nud,Eavg 
      if(autostop.gt.0..and.ttot.gt.autostop*tsi*tse/(tsi+tse))then
         if(me2.eq.0)write(*,100)
     $        ' STOPPING AT ',autostop,' SLOWING DOWN TIMES'
         d=max(d,ndiagnostics)
      endif
      INQUIRE(FILE='checkpoint',EXIST=file)
      if(file.and.d.ge.0.and.runctrl.lt.2)then
         if(me2.eq.0)write(*,'(a)')
     $        ' STOPPING FOR CHECKPOINT'
         d=max(d,ndiagnostics)
      endif
         
 100  format(a,f4.1,a)
      
      if(me2.eq.0.and.d.eq.ndiagnostics.and.runctrl.eq.2)
     $     write(*,'(a,f6.2,a)')' Hot minority',(1.-sum(densth)
     $     /(sum(dens) + sum(densth)))*100,'%'

      endif

      


      deallocate(indxr,indxp,indxt)
      deallocate(pperp,ppar,dens,currpp,temptot,ssum)

      end subroutine moments
         

c------------------------------------------------------------------

      subroutine initialise
      use dimarray
      use dstfun
      use particle
      use plasma
      use interp
      use pival
      implicit none

      integer i,j,k
      real theta,phi
      integer nes,nth,nph
      real ess,Vp,sm
      real profile,pden,ptherm,np
c     real dummV(0:nbinsrad)


      

      allocate(intervs(0:nbinsrad),sl(0:nbinsrad),s(0:nbinsrad)
     $     ,V(nbinstor,nbinspol,0:nbinsrad),Arad(0:nbinsrad))
      allocate(nbinspolr(0:nbinsrad),writepolr(0:nbinsrad))
      allocate(nmarks(nbinstor,nbinspol,0:nbinsrad)
     $     ,nmarksloc(nbinstor,nbinspol,0:nbinsrad)
     $     ,nmarksall(nbinstor,nbinspol,0:nbinsrad))
      allocate(densang(nbinstor,nbinspol,0:nbinsrad)
     $     ,densthang(nbinstor,nbinspol,0:nbinsrad)
     $     ,pparang(nbinstor,nbinspol,0:nbinsrad)
     $     ,pperpang(nbinstor,nbinspol,0:nbinsrad))
      allocate(BpB(nbinstor,nbinspol,0:nbinsrad))
      allocate(currt(0:nbinsrad),currp(0:nbinsrad))
      if(vcurr.eq.1)allocate(currtot(0:nbinsrad))
      allocate(dsb(0:nbinsrad),densth(0:nbinsrad),pth(0:nbinsrad))
      allocate(tensall(nbinsrad+2))
      allocate(Ecrit(0:nbinsrad))
      

      intervs(0)=0.!5/nsa!smin  ! no smaller than 0.5/nsa because of modB
      if(nbinsrad.eq.nsa)intervs(0)=0.
      dsb(0)=intervs(0)
      sm=intervs(0)
      do i=1,nbinsrad
         intervs(i)=sm+(smax-sm)*(real(i)/nbinsrad)**profexp
         dsb(i)=intervs(i)-intervs(i-1)

c         if(polexp.ne.1.)then
c            nbinspolr(i)=nbinspol + nint((1.-polexp)*(nbinsrad-i))
c         else
c            nbinspolr(i)=nbinspol - nint((1.-polexp)*(i-1))
c         endif
         if(polexp.ge.1.)then
            nbinspolr(i)=nbinspol/polexp
     $           + nbinspol*(1.-1./polexp)*(i-1)/(nbinsrad-1)
         else
            nbinspolr(i)=nbinspol
     $           +nbinspol*(polexp-1.)*(i-1)/(nbinsrad-1)
         endif            
         nbinspolr(i)=max(nbinspolr(i),3)
         nbinspolr(i)=min(nbinspolr(i),nbinspol)

         if(writepol.eq.0)then
            writepolr(i)=(nbinspolr(i)+1)/2
         else
            writepolr(i)=min(writepol,nbinspolr(i))
         endif
         
         ess=0.5*(intervs(i)+intervs(i-1))
         profile=n0+n1*ess+n2*ess**2.+n3*ess**3.+n4*ess**4.+n5*ess**5. 
         pden=max(0.,density*profile)
         if(pden.eq.0.)print*,ess,profile
         np=sum(densi(:,:,ceiling(ess*nsa)))/nbinspol/nbinstor
         np=(pden-Zpart*np)*Qel/qplasma
         call interpPP(ess,ptherm)
         Te=ptherm/(np*TiTe+pden) ! [J]
c         Ecrit(i)=14.8*Apart*Te/((mplasma/mp)**(2./3.))
         Ecrit(i)=splitfact*Te
      enddo
      nbinspolr(0)=nbinspolr(1)
      writepolr(0)=0
      Ecrit(0)=Ecrit(1)
c$$$      intervs(0)=0.!smin!0.5/nsa
c$$$      do i=1,nbinsrad
c$$$c$$$         intervs(i)=smin+i*(smax-smin)/nbinsrad
c$$$         intervs(i)=real(i)/nbinsrad
c$$$         dsb(i)=intervs(i)-intervs(i-1)
c$$$      enddo

      plasmaV=0.
      nes=100
      nth=100
      nph=100
      do i=1,nes
         ess=(i-0.5)/nes
         do j=1,nth
            theta=twopi*(j-0.5)/nth
            do k=1,nph
               phi=twopiPer*(k-0.5)/nph
               call interp3DRGC(ess,theta,phi,val,vals,pert)
               Vp=(-vals(6)*vals(2)+vals(5)*vals(1))/val(1)/val(5) ! jacobian
               plasmaV=plasmaV+Vp
            enddo
         enddo
      enddo
      plasmaV=plasmaV/nes*twopi/nth*twopi/nph ! Lper times (2pi/Lper)/nph
      nall=navg*plasmaV
 
      end

c---------------------------------------------------------------
      subroutine volume
      use dimarray
      use dstfun
      use particle
      use interpos
      use pival
      use mpi_var !debug
      implicit none

      integer i,j,k,l,ndum
      real val(10),vals(8),pert(8),val1(4)
      real theta,phi,Vtmp,stmp,dstmp
      integer niter
      complex ef

c      smin=0.5/nsa
c      smax=1.-smin
c      smin=0.
      niter=100
      ndum=1
      

      do i=0,nbinsrad
         Arad(i)=0.
         do j=1,nbinspolr(i)
            theta=-pi+twopi*(j-0.5)/nbinspolr(i)
            do k=1,nbinstor
               phi=twopiPer*(k-0.5)/nbinstor
               V(k,j,i)=0.
c     theta=twopi*(j-0.5)/nbinspol
               Vtmp=0.
               if(i.eq.0)then
                  stmp=0.5*intervs(i)
                  dstmp=dsb(i)
                  call interp3DRGC(stmp,theta,phi,val,vals,pert)
                  Vtmp=(-vals(6)*vals(2)+vals(5)*vals(1))/val(1)/val(5) ! jacobian
                  Vtmp=Vtmp*dstmp*twopi/nbinspolr(i)*twopiPer/nbinstor ! real 3D space volume
                  V(k,j,i)=V(k,j,i)+Vtmp
               else
                  dstmp=intervs(i-1)
                  do l=1,niter
                     stmp=intervs(i-1)+real(l)/niter*dsb(i)
                     dstmp=stmp-dstmp
                     call interp3DRGC(stmp,theta,phi,val,vals,pert)
                     Vtmp=(-vals(6)*vals(2)+vals(5)*vals(1))
     $                    /val(1)/val(5) ! jacobian
                     Vtmp=Vtmp*dstmp*twopi/nbinspolr(i) ! real 3D space volume
     $                    *twopi/nbinstor ! Lper times twopiPer/nbinstor
                     dstmp=stmp
                     V(k,j,i)=V(k,j,i)+Vtmp
                  enddo
               endif
c     V(j,i)=1. ! weights already in boozer coords
c     V(j,i)=V(j,i)*dsb(i)*twopi/nbinspol ! real 2D space volume
               call interpE(ndum,s(i),theta,phi,val1,ef)
               Arad(i)=Arad(i)+V(k,j,i)/(twopi*val1(1)) ! A=V/2*pi*R
c     BpB(j,i)=-vals(1)/val(5)/sqrt(val(1)) ! sigmaB_phi=-mu0*I(s)
c     BpB(j,i)=vals(5)/(-vals(6)*vals(2)+vals(5)*vals(1)) ! part of J_tor which depends on s only
c     BpB(j,i)=BpB(j,i)*sqrt(val(1)) 
c     BpB(j,i)=BpB(j,i)*twopi/nbinspol !remove pol bin for <j> 
               BpB(k,j,i)=vals(5)/(val(5)*sqrt(val(1)))
            enddo
         enddo
      enddo      
      end subroutine volume

