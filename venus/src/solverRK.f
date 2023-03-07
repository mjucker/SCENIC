      subroutine solverRK(Ec,lambdaip,Cmu)
      use futils
      use particle
      use dimarray
      use machine
      use mpi_var
      use sinitial
      use scatter
      use kicks
      use interp
      use rand_var
      use dstfun
      use pival
      implicit none

      integer N,i,j,k,ip,step,dummi,fid
      integer vstep,madecol,ncoulcol,d,ipall
      parameter (N=4)
      real f(N),past,Tpall,tlastcol,eps(0:nbinsrad),dummy,dummy2
      real,allocatable :: wn(:,:),wnall(:,:),zyn(:,:)
      real,allocatable,dimension(:) :: zmu,
     $     Ecall,vperp,vparall
      real Ep,jac,collpast,Cmu(ntot/nprocs),mu,ynt(N),yn(N),
     $     norm,B,slim,normp,rholim,lambda,
     $     lambdaip(ntot/nprocs),ppar,Ec(ntot/nprocs),
     $     interp3DB,smax3,nud,nudpast
     $     ,tsi,tse,nuepast,dnuedee,dnuedei,dnuede
      real vpar,deltas,Epall,thtmp,thtmp1,fs,fs1,Z1
      real tfinloc,ttot,realt,Bp,savgloc,redstep
      integer jt,js,nlost,nlostall,dstart,rctrlloc
      integer check,chk,ind,indp,icol,ic,pickip,Zck,cntcol
      real val1(4),lamchk,rand,Ptot,tavg,hdfarr(2)
      complex ef(2)
      character(20) :: form,dummc
      real,allocatable :: CIvall(:,:),lambdaall(:),Eall(:)
     $     ,weightall(:)
      logical file_exists


      allocate(zmu(ntot/nprocs),
     $     zyn(N,ntot/nprocs),Ecall(ntot/nprocs),
     $     vperp(ntot/nprocs),vparall(ntot/nprocs))

      if(runctrl.eq.-2)allocate(savg(ntot/nprocs))
      if(runctrl.eq.2)allocate(trap(ntot/nprocs))


      slim=min(sinit+nsa/(2.*nsa),smax)
c      deltas=slim-starts   !   ds  annulus size
      deltas=1.


      norm=abs(Qe)*b00/mo
      eps=a00*sqrt(s)/r00

      

      rholim=1.e-5   ! init de rholim pour rkmixte (ion~2-4e-5.electron~2-4e-6)
                     ! a modifier aussi dans rungekuttamixte.f
      icol=1

c$$$      if(me2.eq.0)then
c$$$      ip=0
c$$$      do ipall=me2+1,ntot,nprocs
c$$$         ip=ip+1
c$$$         Ec(ip)=2e6*Qel
c$$$         CIv(1,ip)=0.22
c$$$         i=mod(ipall-1+npitch,npitch)+1
c$$$         rand=random(ran_index(i),ran_array(1,i))
c$$$         CIv(2,ip)=0.!rand*twopi
c$$$         CIv(3,ip)=0.
c$$$         call interp3DRGC(CIv(1,ip),CIv(2,ip),CIv(3,ip),
c$$$     $     val,vals,pert)
c$$$         B=sqrt(val(1))
c$$$         rand=random(ran_index(i),ran_array(1,i))
c$$$         lambdaip(ip)=-0.5*sqrt(1.-B/Bc)!2.*(rand-0.5)!2.948e-1
c$$$         Cmu(ip)=Ec(ip)*(1.-lambdaip(ip)**2.)/B
c$$$         normp=sqrt(2.*Ec(ip)*mo)
c$$$         CIv(4,ip)=lambdaip(ip)*normp/(Qe*B*val(5))
c$$$      enddo
c$$$         print*,B/Bc
c$$$      endif
         
         
c------------------------------------------------------
c DEBUT INTEGRATION
c------------------------------------------------------
c initial  condition
      nkicks=0
      nkicklocal=0
      nlost=0
      realt=0.
      if(runctrl.eq.0)lostweight=0.

C######################################################################C
C     DIAGNOSTICS LOOP                                                 C
C######################################################################C
c$$$      if(me2.eq.0.and.runctrl.eq.2)write(*,'(a25,i3)')
c$$$     $     ' Time averaging, ndiag = ',ndiagnostics
      if(abs(runctrl).eq.2)then
         dstart=0   !find out if particle trapped or passing, or find savg
      else
         dstart=1
      endif
      d=dstart-1
      do while(d.lt.ndiagnostics)
      d=d+1
C######################################################################C
C     PARTICLE LOOP                                                    C
C######################################################################C
      ip=0
      do ipall=me2+1,ntot,nprocs   ! main particle loop
      ip=ip+1
      
      madecol=0
      ncoulcol=0
      ttot=0.
      tlastcol=0.
      redstep=1. 
      cntcol=icol
      if(d.eq.0)then
         if(runctrl.eq.2)then
            tavg=0.
            trap(ip)=0
            Zck=0
         else
            Zck=0
            tavg=0.
c            savg(ip)=CIv(1,ip)
            savg(ip)=0.
         endif            
      endif

      Ep=Ec(ip)
      normp=sqrt(2.*Ep*mo)
      call interp3DRGC(CIv(1,ip),CIv(2,ip),CIv(3,ip),
     $     val,vals,pert)
      B=sqrt(val(1))
      lambda=lambdaip(ip)    
      if(Cmu(ip)*B/Ep.gt.1.and.d.eq.1)then
         Cmu(ip)=Ep/B
      endif
     
      zmu(ip)=Cmu(ip)
      do i=1,N                  ! initialisation
        zyn(i,ip)=CIv(i,ip)
      enddo 
      
      
c   determine time step
      if(lengthpas.eq.0.)then
         if(ipall.eq.1.and.d.eq.1.and.runctrl.ne.2)
     $        print*,'variable time step'
c          ind=max(1,floor(zyn(1,ip)*nbinsrad))
c         past=eps(ind)*mo/normp*(1.+1.*zmu(ip)*(B/Ep)) !kovanen, J.Comp.Phys. 105,14
         past=0.2*mo/normp*(1.+1.*zmu(ip)*(B/Ep)) !kovanen, J.Comp.Phys. 105,14
         if(coulomb.eq.1)then
            call updatenu(Ep,zyn(:,ip),nud,tsi,tse
     $           ,dnuedee,dnuedei)
             past=min(0.1/max(nud,1./tsi+1./tse),past)
         endif
         past=max(1.,past*norm)
         vstep=1
      elseif(ip.eq.1.and.d.eq.dstart)then
         past=tfin*norm/lengthpas 
         vstep=0
      endif
      


c end initial condition
C##########################################################################C
C           TIME LOOP BETWEEN TWO DIAGNOSTICS                              C
C##########################################################################C
      step=1
      if(runctrl.eq.2.and.d.gt.0)tfin=ndiagnostics*5.e-7
      tfinloc=tfin/ndiagnostics
      if(d.eq.0)tfinloc=10.!10000*past/norm ! max num steps for savg
      if(icrh.eq.1.and.(runctrl.eq.0.or.runctrl.eq.1))then
         if(tfin/ndiagnostics.gt.2.e-04)then
            tfinloc=(tfin-2.e-04)/(ndiagnostics-2)
            if(d.le.2)tfinloc=1.e-04
         else
            tfinloc=tfin/ndiagnostics
         endif
      endif
!      past=min(past,abs(tfinloc-ttot)*norm) !in order to finish at t=tfin
      if(weight(ip).eq.0.)then
         ttot=tfinloc
         do i=1,N
            yn(i)=zyn(i,ip)
         enddo
         yn(2)=mod(mod(yn(2),twopi)+twopi,twopi)
         yn(3)=mod(mod(yn(3),twopi)+twopi,twopi)
      endif

      do while (ttot.lt.tfinloc)!(step.le.lengthpas(ip))       !   main time  loop

         
         do i=1,N
            yn(i)=zyn(i,ip)
         enddo
         yn(2)=mod(mod(yn(2),twopi)+twopi,twopi)  !0->2pi
         yn(3)=mod(mod(yn(3),twopi)+twopi,twopi)  !0->2pi
         
         mu=zmu(ip)
c#####################################################################
c$$$        if(ipall.eq.1.and.runctrl.lt.2)
c$$$c     $        .and.(d.eq.3.or.d.eq.3))
c$$$     $       write(99,'(1p5e16.8)')
c$$$     $       (d-1)*tfinloc+ttot,yn(1),yn(2),yn(3),yn(4)
c$$$c     $       mu,yn(1),yn(2)/pi,yn(3)/pi,yn(4)
c SOLVER   make step
        call rungekuttamixte(past,yn,ynt,Ep,mu,
     >       b00,rholim,smax3,check)
         ttot=ttot+past/norm    !real simulated time

c#####################################################################
         
c########## LOST PARTICLES ###########################################
cj  particle lost on wall or low energy - duplicate existing ptcle
         if(weight(ip).eq.0.)then
            ttot=tfinloc
         elseif(yn(1).ge.smax)then
            lostweight=lostweight+weight(ip)
            indp=1
            if(stat.gt.0)indp=floor(mod(mod(yn(2),twopi)+twopi,twopi)
     $           /twopi*nbinspol)+1
            lostE(indp)=lostE(indp)+Ep*weight(ip)
            nlost=nlost+1
cc           weight(ip)=0.
cc           ttot=tfinloc
cc           yn(1)=smax
cc        elseif(runctrl.ge.0.and.Ep.lt.0.*Es)then
            chk=0
            do while(chk.eq.0)
               i=mod(ipall-1+npitch,npitch)+1
               rand=random(ran_index(i),ran_array(:,i))
c           pickip=ip+ceiling(rand*(ntot/nprocs-ip)) ! pick ptcle at ttot=0
               pickip=ceiling(rand*ntot/nprocs)
               Ep=Ec(pickip)
               if(Ep.gt.Es.and.weight(pickip).gt.0.)chk=1
            enddo
            yn(:)=CIv(:,pickip)
            lambda=lambdaip(pickip)
            mu=Cmu(pickip)
            weight(ip)=0.5*weight(pickip)
            weight(pickip)=weight(ip)
            if(pickip.gt.ip)then
               ttot=past/norm   ! ptcle not integrated to tfinloc yet
               tlastcol=0.
               Enkicksradip(:,ip,:)=0.
               if(stat.gt.0)then
                  Enkicksip(:,:,:,ip)=0.
                  Encollip(:,:,ip)=0.
               endif
            else
               ttot=tfinloc     ! ptcle already integrated to tfinloc
               Enkicksradip(:,pickip,:)=0.5*Enkicksradip(:,pickip,:)
               Enkicksradip(:,ip,:)=Enkicksradip(:,pickip,:)
               if(stat.gt.0)then
                  Enkicksip(:,:,:,pickip)=0.5*Enkicksip(:,:,:,pickip)
                  Encollip(:,:,pickip)=0.5*Encollip(:,:,pickip)
                  Enkicksip(:,:,:,ip)=Enkicksip(:,:,:,pickip)
                  Encollip(:,:,ip)=Encollip(:,:,pickip)
               endif
            endif
         endif
         call interp3DRGC(yn(1),yn(2),yn(3),val,vals,pert)
         B=sqrt(val(1))

c################# FIND <S> ###########################################
        if(d.eq.0)then
           if(step.eq.1)then
              call interpE(d+1,yn(1),yn(2),yn(3),val1,ef) !val1(1)=R,val1(2)=Z
              Z1=val1(2)
           else
              call interpE(d+1,yn(1),yn(2),yn(3),val1,ef) !val1(1)=R,val1(2)=Z
              if(Z1*val1(2).lt.0.)Zck=Zck+1
              Z1=val1(2)
           endif

           if(Zck.gt.0)then
c$$$              tavg=tavg+past
              if(Zck.eq.3)then
                 ttot=tfinloc
c$$$                 write(17,'(3es14.6)')val1(1),mu*Bc/Ep,tavg/norm
              endif
              if(runctrl.eq.-2)then
                 savg(ip)=savg(ip)+yn(1)*past
                 tavg=tavg+past
              endif
           endif
           
           if(yn(4)*CIv(4,ip).lt.0..and.runctrl.eq.2)then
              trap(ip)=1
              ttot=tfinloc
           endif

           if(stat.eq.2.and.runctrl.eq.2)ttot=tfinloc

           if(step.eq.1e5)then
              print*,ipall,'<s> not terminated'
              ttot=tfinloc
           endif
        endif
c#####################################################################
c MAKE DIFFUSION. MONTE CARLO COULOMB (PITCH,E) & ICRH (VPERP,VPAR)  #
c#####################################################################
       normp=sqrt(2.*Ep*mo)
       lambda=yn(4)*Qe*B*val(5)/normp ! v_par/v     
       if(abs(lambda).gt.1..and.weight(ip).ne.0.)then
!          print*,'ip,d,step,lambda,icol',ipall,d,step,lambda,icol
c          lambda=sign(.999,lambda)    ! avoid numerical imprecision 
          lambda=sign(2.,lambda)-lambda   ! avoid numerical imprecision 
          yn(4)=lambda*normp/(Qe*B*val(5))
c          weight(ip)=0.
       endif
c pitch angle and energy scattering with Lorentz operator
       if(coulomb*weight(ip).gt.0.and.ttot.ne.tfinloc)then
          if(cntcol.ge.icol)then
            call updatenu(Ep,yn,nud,tsi,tse
     $           ,dnuedee,dnuedei)
             collpast=ttot-tlastcol
             i=mod(ipall-1+npitch,npitch)+1
             nudpast=nud*collpast
             call changepitch(
     $            yn(1),yn(2),yn(3),lambda,nudpast,Ep,collpast
     $            ,tsi,tse,dnuedee,dnuedei,i,ip)
             madecol=1
             ncoulcol=ncoulcol+1
             if(Ep.lt.0.)then
                weight(ip)=0.
                Ep=abs(Ep)
             endif
             tlastcol=ttot      !-past/norm + past/norm/icol*ic
             cntcol=1
          else
             cntcol=cntcol+1  ! collision subsampling
          endif
       endif
       
c vperp scattering due to ICRH
       if(icrh*weight(ip).gt.0..and.ttot.ne.tfinloc)then
          i=mod(ipall-1+npitch,npitch)+1
          call makekicks(yn,Ep,lambda,B,val(5),past/norm,
     $       ip,i,madecol,redstep)
       endif
       

c################# END DIFFUSION ######################################c
       if(madecol.eq.1)then
          if(abs(lambda).gt.1.)lambda=sign(2.,lambda)-lambda
          normp=sqrt(2.*Ep*mo)
          ppar=normp*lambda  
          yn(4)=ppar/Qe/B/val(5)
       endif

       do i=1,N
          zyn(i,ip)=yn(i)
       enddo
       mu=abs(Ep*(1.-lambda*lambda)/B)
       zmu(ip)=mu

c   determine new time step
       if(vstep.eq.1.and.weight(ip).ne.0.)then
          past=0.2*mo/normp*(1.+1.*mu*(B/Ep))
          if(coulomb.gt.0)then
            call updatenu(Ep,yn,nud,tsi,tse
     $           ,dnuedee,dnuedei)
             past=min(0.1*icol/max(nud,1./tsi+1./tse),past)
          endif
          past=past*redstep
          if(runctrl.eq.2.and.d.gt.0)past=min(past,abs(tfinloc-ttot)) !in order to finish at t=tfin
          past=past*norm
          madecol=0
          redstep=1.
       endif

       step=step+1
      enddo
C#######################################################################C
C        END TIME LOOP PER DIAGNOSTICS TIME                             C
C#######################################################################C
      if(d.eq.0)then
         do i=1,N
            yn(i)=CIv(i,ip)
         enddo
         zmu(ip)=Cmu(ip)
         lambda=lambdaip(ip)
         Ep=Ec(ip)
         call interp3DRGC(yn(1),yn(2),yn(3),val,vals,pert)
         B=sqrt(val(1))
      endif
         
      Ec(ip)=Ep
      if(weight(ip).eq.0..and.abs(lambda).gt.1.)lambda=sign(.99,lambda)
      lambdaip(ip)=lambda
      Cmu(ip)=zmu(ip) 
      do i=1,N
         CIv(i,ip)=yn(i)
      enddo 
!      if(weight(ip).eq.0.)CIv(1,ip)=smax
      if(runctrl.eq.-2.and.weight(ip).eq.0.)tavg=ttot
      if(d.eq.0.and.runctrl.eq.-2)then
         savg(ip)=savg(ip)/tavg
         if(savg(ip).gt.1..and.weight(ip).ne.0.)then
            print*,'ipall,savg,tavg,step',ipall,savg(ip),tavg,step
            weight(ip)=0.
         endif
      endif

      vparall(ip)=lambda*sqrt(2.*Ep/mo) 
      vperp(ip)=sqrt(2.*Ep/mo-vparall(ip)*vparall(ip))


      if(zmu(ip)*B/Ep.gt.1.)then
         if(weight(ip).gt.0.)print*,'pitch>1',ipall,d,zmu(ip)*B/Ep !(vperp/v)^2 
         mu=0.99*Ep/B
         weight(ip)=0.
      endif
      

      enddo
C#######################################################################C
C      END PARTICLE LOOP                                                C
C#######################################################################C
      realt=realt+ttot
      if(runctrl.ge.0)then
c      if(runctrl.eq.2)then
         call moments(d,realt,vparall,vperp)
         if(nbinsvpar*nbinsvperp.ne.0)
     $        call distfun(d,realt,vparall,vperp)


         INQUIRE(FILE=filename, EXIST=file_exists) ! create state.vs if not existing
         if(runctrl.ne.2.or..not.file_exists)then      
            hdfarr(1)=lostweight
            hdfarr(2)=nall
            allocate(CIvall(4,ntot),lambdaall(ntot)
     $           ,Eall(ntot),weightall(ntot))
            do i=1,4
               call MPI_GATHER(CIv(i,:),ntot/nprocs
     $              ,MPI_DOUBLE_PRECISION
     $              ,CIvall(i,:),ntot/nprocs
     $              ,MPI_DOUBLE_PRECISION,0,comm,ierr)
            enddo
            call MPI_GATHER(lambdaip,ntot/nprocs,MPI_DOUBLE_PRECISION
     $           ,lambdaall,ntot/nprocs
     $           ,MPI_DOUBLE_PRECISION,0,comm,ierr)
            call MPI_GATHER(Ec,ntot/nprocs,MPI_DOUBLE_PRECISION
     $           ,Eall,ntot/nprocs
     $           ,MPI_DOUBLE_PRECISION,0,comm,ierr)
            call MPI_GATHER(weight,ntot/nprocs,MPI_DOUBLE_PRECISION
     $           ,weightall,ntot/nprocs
     $           ,MPI_DOUBLE_PRECISION,0,comm,ierr)
            if(me2.eq.0)then
               call creatf(filename,fid)
               call putarr(fid,'/CIv',CIvall)
               call putarr(fid,'/hdfarr',hdfarr)
               call putarr(fid,'/lambda',lambdaall)
               call putarr(fid,'/Ec',Eall)
               call putarr(fid,'/weight',weightall)
               call closef(fid)
            endif
            deallocate(CIvall,lambdaall,Eall,weightall)
            
            if(stat.ne.0.or.icrh.eq.1)call writestat(d,tfinloc,ttot)
         else
            if(d.eq.0.and.me2.eq.0)then
               print*,' trap/pass determined'
            elseif(mod(d,max(1,ndiagnostics/10)).eq.9.and.me2.eq.0)then
               write(*,'(a,$)'),' #'
            endif
         endif
      endif
      enddo
C#######################################################################C
C      END DIAGNOSTICS LOOP                                             C
C#######################################################################C
      if(runctrl.eq.0.or.runctrl.eq.1)then
         if(icrh.gt.0)then
            call MPI_REDUCE(kparavg,kparavgtot,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,comm,ierr)
            call MPI_REDUCE(kperpavg,kperpavgtot,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,comm,ierr)
            call MPI_REDUCE(Ravg,Ravgtot,1,MPI_DOUBLE_PRECISION,
     $           MPI_SUM,0,comm,ierr)
            call MPI_REDUCE(kickfactor,kickfactot,1
     $           ,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,ierr)
            call MPI_REDUCE(kickfactor2,kickfactot2,1
     $           ,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,ierr)
            call MPI_REDUCE(nkicks,nkickstot,1,MPI_INTEGER8
     $           ,MPI_SUM,0,comm,ierr)
         endif

         
         nlost=0
         do ip=1,ntot/nprocs
            if(weight(ip).eq.0.)then
               nlost=nlost+1
            endif
         enddo
         call MPI_REDUCE(nlost,nlostall,1,MPI_INTEGER,
     $        MPI_SUM,0,comm,ierr)
         
         
         if(me2.eq.0)then
            if(icrh.gt.0)then
               write(*,'(a21,3f9.2)')
     $              ' <R>,<k_par>,<k_perp> ',
     $              Ravgtot/nkickstot
     $              ,kparavgtot/nkickstot,kperpavgtot/nkickstot
               write(*,'(a24,3es8.1E1)')
     $              ' <nkicks>,<kick _|_ ||> ',real(nkickstot)/ntot
     $              ,kickfactot/nkickstot,kickfactot2/nkickstot
            endif
               
c$$$  open(unit=5,file='result.out')
            write(*,'(1x,a,f6.2,a,es6.1E1,a)')'Lost particles '
     $           ,lossall/(totalweight+lossall)*100.
     $           ,'% (',real(nlostall),' ptcles)'
         endif
         
         
 101     format(6es24.15e3)     !,i12)
   


      elseif(runctrl.eq.-2)then
         
         deallocate(CIv,zmu,zyn,Ecall,vperp,vparall)
         call writearray
         deallocate(weight)

      endif


      if(abs(startpos).le.1)then
         allocate(wn(6,ntot/nprocs))
         allocate(wnall(6,ntot))
         do ip=1,ntot/nprocs
            wn(1:3,ip)=CIv(1:3,ip)
            wn(4,ip)=lambdaip(ip)
            wn(5,ip)=Ec(ip)/Qel
            if(weight(ip).eq.0.)then
               wn(6,ip)=0.
            else
               wn(6,ip)=weight(ip)!Cmu(ip)*b00/Ec(ip)
            endif
         enddo
         
         do i=1,6
            call MPI_GATHER(wn(i,:),ntot/nprocs,MPI_DOUBLE_PRECISION,
     $           wnall(i,:),ntot/nprocs
     $           ,MPI_DOUBLE_PRECISION,0,comm,ierr)
         enddo
         
         
         deallocate(wn)
         
         if(me2.eq.0)then
            open(unit=3,file='positions.out')
            do ip=1,ntot
               write(3,101)wnall(1:6,ip)
            enddo
            close(3)
            deallocate(wnall)
         endif
      endif
         

      if(runctrl.ge.0)then
         if(runctrl.eq.2)deallocate(trap)
         if(runctrl.eq.-2)deallocate(savg)
         deallocate(zmu)
         deallocate(zyn,Ecall)
         deallocate(vperp,vparall)
      endif

      end


c--------------------------------------------------------------------
c$$$      subroutine reload(yn,mu,Ep,B,w,i)
c$$$      use dstfun
c$$$      use particle
c$$$      use dimarray
c$$$      use para
c$$$      use rand_var
c$$$      use interp
c$$$      use machine
c$$$      use sinitial
c$$$      implicit none
c$$$
c$$$      real yn(4),mu,Ep,B,w
c$$$      real pi,twopi,jac,radcr,vpar,vperp
c$$$      real vmax,dstfn,x,x2,x3,dsi,si2
c$$$      real Tloc,Tpar,Tperp,densround,fjacmax,fjac
c$$$      integer is,check,i
c$$$      
c$$$
c$$$      pi=2.*acos(0.)
c$$$      twopi=2.*pi
c$$$      
c$$$      vmax=sqrt(2.*Tmax*Qel/mo)
c$$$      distfn=4.*pi*vmax*vmax*(smax-smin)  ! 2D distribution (s,theta)
c$$$      distfn=ntot/distfn
c$$$
c$$$      x=random(ran_index(i),ran_array(1,i)) 
c$$$      yn(2)=twopi*x         ! theta
c$$$      x=random(ran_index(i),ran_array(1,i)) 
c$$$      yn(3)   =twopi*x ! phi
c$$$      x=random(ran_index(i),ran_array(1,i))
c$$$      yn(1)=smin+x*(smax-smin)      ! radially inside grid  ! s
c$$$      
c$$$      is=floor(x*nsa)
c$$$      si2=is*ds+ds2
c$$$      dsi=(si2-x)/ds
c$$$      if(dsi.lt.0.)then
c$$$         is=is+1
c$$$         si2=is*ds+ds2
c$$$         dsi=(si2-x)/ds
c$$$      endif
c$$$      dsi=1.-dsi
c$$$      
c$$$      Tloc=dsi*Tpar_s(is+1)+(1.-dsi)*Tpar_s(is)
c$$$      Tpar=Tloc*Qel             ! [J]
c$$$      Tperp=(dsi*anis_s(is+1)+(1.-dsi)*anis_s(is))*Tpar
c$$$      densround=dsi*dens_s(is+1)+(1.-dsi)*dens_s(is)
c$$$      
c$$$      call interp3DRGC(yn(1),yn(2),yn(3),
c$$$     $     val,vals,pert)
c$$$      B=sqrt(val(1))
c$$$      jac=-(vals(6)*vals(2) -  vals(5)*vals(1))/val(1)/val(5) ! jacobian   
c$$$      
c$$$      fjacmax=densround*((mo/twopi)**1.5)/(sqrt(Tpar)*Tperp)
c$$$      radcr=1.
c$$$c      radcr=1.-CIv(1,ip)
c$$$      check=0
c$$$      do while(check.eq.0)
c$$$         x2=random(ran_index(i),ran_array(1,i)) ! for vperp
c$$$         x3=random(ran_index(i),ran_array(1,i)) ! for vpar
c$$$         vperp=x2*vmax*radcr
c$$$         x3=2.*x3 - 1.          ! x3=vpar/v in (-1,1)
c$$$         vpar=x3*vmax*radcr
c$$$               
c$$$         fjac=fjacmax
c$$$     $        *exp(-mo/2.*(Bc/B*vperp*vperp/Tperp
c$$$     $        +abs(vpar*vpar+(1.-Bc/B)*vperp*vperp)/Tpar))
c$$$         
c$$$         w=fjac*jac*vperp/distfn*radcr*radcr
c$$$         if(w.gt.0.)check=1
c$$$      enddo
c$$$               
c$$$
c$$$c other quantities needed
c$$$      Ep=0.5*mo*(vpar*vpar+vperp*vperp)
c$$$      mu=0.5*mo*vperp*vperp/B      
c$$$      yn(4)=mo*vpar/(abs(Qe)*B*val(5))
c$$$ 
c$$$      end

c--------------------------------------------------------------
      subroutine writearray
      use particle
      use machine
      use dstfun
      use mpi_var
      use para
      use scatter
      use dimarray
      use interp
      use futils
      implicit none


      real,allocatable,dimension(:) :: so,lambda,Ec,th
      integer ila,ien,isa,ith,fid,dims(4)
      character*40 fname
      real Emax,twopi,Ewrite(2)
      real,allocatable :: ss(:,:,:,:),st(:,:,:,:)
      integer,allocatable :: markers(:,:,:,:),mtot(:,:,:,:)
      integer ip

      allocate(so(ntot/nprocs),lambda(ntot/nprocs),Ec(ntot/nprocs)
     $     ,th(ntot/nprocs))
      allocate(ss(nbinslam,nbinsen,nbinsrad,nbinspol)
     $     ,markers(nbinslam,nbinsen,nbinsrad,nbinspol))

      write(fname,'(a4,i4.4)')'load',me2
      open(unit=12,file=fname,form='unformatted')
      read(12)so
      read(12)th
      read(12)lambda
      read(12)Ec
      close(12,status='delete')

      twopi=4.*acos(0.)
      Emax=Tmax*Qel
      lambda=lambda+1.

      markers=0
      ss=0.
      do ip=1,ntot/nprocs
         if(weight(ip).gt.0.)then
            ila=floor(lambda(ip)/2.*nbinslam)+1
            ien=floor((Ec(ip)-Es)/(Emax-Es)*nbinsen)+1
            isa=floor(savg(ip)*nbinsrad)+1
            ith=floor(th(ip)/twopi*nbinspol)+1
            
            ss(ila,ien,isa,ith)=ss(ila,ien,isa,ith) + so(ip)*weight(ip)
            markers(ila,ien,isa,ith)
     $           =markers(ila,ien,isa,ith) + nint(weight(ip))
         endif
      enddo
      deallocate(so,lambda,Ec,th)
      allocate(st(nbinslam,nbinsen,nbinsrad,nbinspol)
     $     ,mtot(nbinslam,nbinsen,nbinsrad,nbinspol))
      mtot=0
      call MPI_REDUCE(markers,mtot,nbinslam*nbinsen*nbinsrad*nbinspol
     $     ,MPI_INTEGER,MPI_SUM,
     $     0,comm,ierr)
      st=0.
      call MPI_REDUCE(ss,st,nbinslam*nbinsen*nbinsrad*nbinspol
     $     ,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     0,comm,ierr)
      
      if(me2.eq.0)then
         do ith=1,nbinspol
            do isa=1,nbinsrad
               do ien=1,nbinsen
                  do ila=1,nbinslam
                     if(mtot(ila,ien,isa,ith).eq.0)then
                        mtot(ila,ien,isa,ith)=-1
                        st(ila,ien,isa,ith)=2.
                     endif
                     st(ila,ien,isa,ith)
     $                 =st(ila,ien,isa,ith)/mtot(ila,ien,isa,ith)
                  enddo
               enddo
            enddo
         enddo

         dims(1)=nbinslam
         dims(2)=nbinsen
         dims(3)=nbinsrad
         dims(4)=nbinspol
         Ewrite(1)=Es
         Ewrite(2)=Emax

         call creatf('loadarray.in',fid)
         call putarr(fid,'/dims',dims)
         call putarr(fid,'/Erange',Ewrite)
         call putarr(fid,'/st',st)
         call closef(fid)
      endif
      deallocate(ss,st,markers,mtot)

      end
      

c---------------------------------------------------------------------
      subroutine writestat(d,tfinloc,ttot)
      use particle
      use mpi_var
      use scatter
      use dstfun
      use pival
      use terp
      implicit none

      integer d,i,j,ip,k
      real tfinloc,ttot,cosph,sinph
      real dummy,dummth,dummph,val1(4)
      real RZC(0:nbinsrad,nbinspol,nbinstor,3)
      double complex ef(2)
      character(80) form

      write(form,"('(a,i3,a,f6.2,a,f6.2,a,'
     $,i2,'i3,a,f6.2,a,f6.2,a)')")nph

      if(d.eq.1)then
         do k=1,nbinstor
            dummph=(k-0.5)/nbinstor*twopi/Lper
            cosph=cos(dummph)
            sinph=sin(dummph)
            do j=1,nbinspol
               dummth=(j-0.5)/nbinspol*twopi
               do i=0,nbinsrad
                  call interpE(1,s(i),dummth,dummph,val1,ef)
                  RZC(i,j,k,1)=val1(1)*cosph
                  RZC(i,j,k,2)=val1(1)*sinph
                  RZC(i,j,k,3)=val1(2)
               enddo
            enddo
         enddo
      endif

      if(icrh.eq.1)then
         if(stat.gt.0)then
            do j=1,nph
               do i=0,nbinsrad
                  Enkicksrad(i,j)=sum(Enkicksradip(i,:,j))
               enddo
            enddo
         else
            do j=1,nph
               Enkicksrad(1,j)=sum(Enkicksradip(1,:,j))
            enddo
         endif
         Enkicksrad=Enkicksrad/totalweight
      endif
      if(stat*coulomb*icrh.gt.0)then
         Encoll=sum(Encollip,3)/totalweight
         Encollall=Encollall/totalweight
         Enkicks=sum(Enkicksip,4)/totalweight
         lostE=lostE/totalweight
         call MPI_REDUCE(Encoll,Encolltot,(nbinsrad+1)*2
     $        ,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,ierr)
         call MPI_REDUCE(Encollall,Encollalltot
     $        ,(nbinsrad+1)*nbinspol*nbinstor
     $        ,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,ierr)
         call MPI_REDUCE(Enkicks,Enkickstot
     $        ,(nbinsrad)*nbinspol*nbinstor
     $        ,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,ierr)
         call MPI_ALLREDUCE(Enkicksrad,Enkicksradtot
     $        ,(nbinsrad+1)*nph
     $        ,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
         call MPI_REDUCE(lostE,lostEtot
     $        ,nbinspol,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,ierr)
         if(me2.eq.0)then
            if(d.eq.1)then
               open(unit=89,file='fort.89'
     $              ,form='unformatted')
               write(89)ndiagnostics,nbinsrad,nbinspol,nbinstor,Lper,nph
               write(89)s
               write(89)RZC !(x,y,z) for Power plots
            endif
            write(89)(d-1)*tfinloc+ttot
            write(89)Encolltot*nall/tfinloc
            write(89)Encollalltot*nall/tfinloc
            write(89)Enkicksradtot*nall/tfinloc
            write(89)Enkickstot*nall/tfinloc
            write(89)lostEtot*nall/tfinloc
            write(*,form)
     $           ' Power lost//ICRH//Coul t(e), d =',d,': '
     $           ,(-sum(lostEtot))*nall/tfinloc/1e6
     $           ,' /// ',sum(Enkicksradtot)*nall/tfinloc/1e6
     $           ,' ('
     $           ,nint(
     $           sum(Enkicksradtot,1)/sum(Enkicksradtot)*100)
     $           ,') // '
     $           ,sum(Encolltot)*nall/tfinloc/1e6
     $           ,' ('
     $           ,sum(Encolltot(:,2))*nall/tfinloc/1e6,') MW'
            if(d.eq.ndiagnostics)close(89)
         endif
         Encollip=0.
         Enkicksip=0.
         lostE=0.
      elseif(icrh.eq.1)then
         if(stat.eq.1)then
            call MPI_ALLREDUCE(Enkicksrad,Enkicksradtot
     $           ,(nbinsrad+1)*nph
     $           ,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
         else
            call MPI_ALLREDUCE(Enkicksrad,Enkicksradtot
     $           ,1*nph
     $           ,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
         endif
         if(me2.eq.0)then
            write(*,'(a,i3,a,f5.2,a)')
     $           ' Power ICRH, d =',d,': '
     $           ,sum(Enkicksradtot)*nall/tfinloc/1e6,' MW'
         endif 
      endif
      if(icrh.eq.1)then
         do j=1,nph
            if(sum(Enkicksradtot(:,j)).gt.0..and.d.ge.1)then
               dummy=sqrt(Powph(j)*DepPow/(sum(Enkicksradtot(:,j))
     $              *nall/tfinloc))
               Efact(j)=dummy*Efact(j)
            endif 
         enddo  
         Enkicksradip=0.
      endif
      
      end subroutine writestat
