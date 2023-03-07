      subroutine makeinputrandom(Ec,lambda,mu)
      use futils
      use particle
      use terp
      use dimarray
      use para
      use machine
      use mpi_var
      use sinitial
      use rand_var
      use scatter
      use interp
      use dstfun
      use pival
      implicit none


      real Tloc,dsi,dth,si2
      integer  i,j,ip,ipall,js,jt,jk,check,is,reload,reldtot
      REAL p,val1(3),B  ,yn(4)
      real x,x2,x3,x1,x4,fjacmax,gjacmax
      real,allocatable,dimension(:) ::  sall,thall,xiall,
     $     lam,Eall,weightall!,totweight
      real,allocatable,dimension(:) ::  Etemp,vparall,vperpall
      real Ec(ntot/nprocs),lambda(ntot/nprocs),mu(ntot/nprocs)
      real Etemptot,Epall,jac,vmax,lamax,lamin
      real radcr,dummy,dummy2,hdfarr(2)
      character*40 fname,dummc
      integer fid,ndum,isnan
      real,allocatable :: wn(:,:),wnall(:,:)
      real ss,thta,A,BoBc,C,dC,theta,phi
      


      if(startpos.gt.0)allocate(sall(ntot),thall(ntot),xiall(ntot),
     $     lam(ntot),Eall(ntot))
      if(runctrl.ge.0)allocate(Etemp(ntot/nprocs),vparall(ntot/nprocs)
     $     ,vperpall(ntot/nprocs))!,totweight(ntot/nprocs))
      allocate(dens_s(0:nsa),anis_s(0:nsa),Tpar_s(0:nsa)
     $     ,densi(nbinstor,nbinspol,nsa))

      smin=0.5/nsa 
      smax=1.-smin
c      smin=1./30
      smin=0.

      lamin=1./nbinslam
      lamax=2.-lamin

      gamma=1.
      
      if(me2.eq.0)then
c  read in density dens(s,theta) and pressures ppar(s,theta), pperp(s,theta) from VMEC/MATLAB
            
         open(unit=10,file='moments.in')
         read(10,*)Bc        
         do js=1,nsa
            read(10,*)dens_s(js),Tpar_s(js),anis_s(js)
            ss=(js-0.5)/nsa
            A=1./sqrt(anis_s(js))
            do jt=1,nbinspol
               theta=(jt-0.5)*twopi/nbinspol
               do jk=1,nbinstor
                  phi=twopiPer*(jk-0.5)/nbinstor
                  call interp3DRGC(ss,theta,phi,val,vals,pert)
                  BoBc=sqrt(val(1))/Bc
                  call calcC(BoBc,anis_s(js),C,dC)
                  densi(jk,jt,js)=dens_s(js)*A*C
               enddo
            enddo
         enddo
         dens_s(0)=1.5*dens_s(1)-0.5*dens_s(2)
         Tpar_s(0)=1.5*Tpar_s(1)-0.5*Tpar_s(2)
         anis_s(0)=1.5*anis_s(1)-0.5*anis_s(2)
         close(10)
      endif
      if(me2.eq.0)write(*,'(1x,a,es9.1)')'<n>',navg
c     distribute read quantities to all processors
      call MPI_BCAST(
     $     navg,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(
     $     Bc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(
     $     Tpar_s,(nsa+1),
     $     MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(
     $     dens_s,nsa+1,
     $     MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(
     $     anis_s,nsa+1,
     $     MPI_DOUBLE_PRECISION,0,comm,ierr)
      call MPI_BCAST(
     $     densi,nbinstor*nbinspol*nsa,
     $     MPI_DOUBLE_PRECISION,0,comm,ierr)
c-------------------------------------------------------------
c  PRODUCTION RUN IF RUNCTRL.GE.0, LOOK-UP TABLE IF -1
c-------------------------------------------------------------
      if(runctrl.ge.0)then
c-------------------------------------------------------------
c  READ INITIAL POSITIONS FROM FILE IF STARTPOS.GT.0
c-------------------------------------------------------------
      if(startpos.gt.0)then
         if(me2.eq.0)then
            open(unit=10,file='moments.in')
            read(10,*)Bc
            close(10)
         endif

         allocate(weight(ntot/nprocs),weightall(ntot))
         if(me2.eq.0)then
            open(unit=3,file='positions.out')
            totalweight=0.
            do i=1,ntot
               read(3,'(6es24.15e2)')sall(i),thall(i),xiall(i),
     $              lam(i),Eall(i),weightall(i)
               totalweight=totalweight+weightall(i)
            enddo
            close(3)
         endif
         
         call MPI_BCAST(sall,ntot,MPI_DOUBLE_PRECISION,0,comm,ierr)
         call MPI_BCAST(thall,ntot,MPI_DOUBLE_PRECISION,0,comm,ierr)
         call MPI_BCAST(xiall,ntot,MPI_DOUBLE_PRECISION,0,comm,ierr)
         call MPI_BCAST(lam,ntot,MPI_DOUBLE_PRECISION,0,comm,ierr)
         call MPI_BCAST(Eall,ntot,MPI_DOUBLE_PRECISION,0,comm,ierr)
         call MPI_BCAST(weightall,ntot,MPI_DOUBLE_PRECISION,0,comm,ierr)
         
         ip=0            
c     distribute quantities among processors
         do ipall=me2+1,ntot,nprocs
            ip=ip+1
            CIv(1,ip)=sall(ipall)
            CIv(2,ip)=thall(ipall)
            CIv(3,ip)=xiall(ipall)
            Ec(ip)=Eall(ipall)*Qel
            weight(ip)=weightall(ipall)
            
            
            call interp3DcoefRGC2(CIv(1,ip),CIv(2,ip),CIv(3,ip),B)
            B=sqrt(B)
            lambda(ip)=lam(ipall)
            mu(ip)=(1.-lambda(ip)*lambda(ip))*Ec(ip)/B ! init de mu
            if(lambda(ip).gt.1.)then
               lambda(ip)=sign(.99,lambda(ip))
            endif
            p=sqrt(2.*Ec(ip)*mo)         
            CIv(4,ip)=p*lambda(ip)/abs(Qe)/B
         enddo
         
         call MPI_BCAST(
     $        Bc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
         
         call MPI_REDUCE(Ec,Etemp,ntot/nprocs,MPI_DOUBLE_PRECISION
     $        ,MPI_SUM,0,comm,ierr)
         
         if(me2.eq.0)then
            Etemptot=0.
            do i=1,ntot/nprocs
               Etemptot=Etemptot+Etemp(i)
            enddo
            write(*,'(a,2es9.1)')
     $           ' Initialisation from file: <E>, tloc  ',
     $           Etemptot/Qel/ntot,tfin
         endif
         
         deallocate(sall,thall,xiall,
     $        lam,Eall)
c--------------------------------------------------------------
c  CREATE NEW INPUT IF STARTPOS.EQ.0
c--------------------------------------------------------------
      else    

         allocate(weight(ntot/nprocs))

         vmax=sqrt(2.*Tmax*Qel/mo)


         if(runctrl.eq.0)then
            if(me2.eq.0)print*,' NEW INPUT'

            call lookup
            call findmaxf(fjacmax,gjacmax)

            ip=0
            do ipall=me2+1,ntot,nprocs
               ip=ip+1
               
               i=mod(ipall-1+npitch,npitch)+1
               call loadptcle(i,ip,fjacmax,gjacmax
     $              ,yn,Ec(ip),mu(ip),lambda(ip),weight(ip))

               CIv(1,ip)=yn(1)
               CIv(2,ip)=yn(2)
               CIv(3,ip)=yn(3)
               CIv(4,ip)=yn(4)
            enddo
         else
            if(me2.eq.0)print*,' LOADING PTCLES'

            call openf(filename,fid,'r',mpicomm=comm)
            call getarr(fid,'/CIv',CIv,pardim=2)
            call getarr(fid,'/hdfarr',hdfarr)
            lostweight=hdfarr(1)
            nall=hdfarr(2)
            call getarr(fid,'/lambda',lambda,pardim=1)
            call getarr(fid,'/Ec',Ec,pardim=1)
            call getarr(fid,'/weight',weight,pardim=1)
            call closef(fid)

            reload=0
            ip=0
            do ipall=me2+1,ntot,nprocs
               ip=ip+1
            
               if(weight(ip).eq.0.)then
                  reload=reload+1
               endif
               
               call interp3DRGC(CIv(1,ip),CIv(2,ip),CIv(3,ip)
     $              ,val,vals,pert)
               mu(ip)=Ec(ip)*(1.-lambda(ip)**2.)/sqrt(val(1))

            enddo
            call MPI_REDUCE(reload,reldtot,1,MPI_INTEGER,
     $           MPI_SUM,0,comm,ierr)
            
            if(me2.eq.0)write(*,'(a,es7.1,a)')
     $           ' Missing ',real(reldtot),' ptcles'
         endif
      endif
      deallocate(dens_s,anis_s,Tpar_s)
c---------Compute initial f(s,theta) and moments-------------------
c  first diagnostics for initial conditions
      ip=0
      do ipall=me2+1,ntot,nprocs
         ip=ip+1
         vparall(ip)=lambda(ip)*sqrt(2.*Ec(ip)/mo)
         vperpall(ip)=sqrt(abs(2.*Ec(ip)/mo-vparall(ip)*vparall(ip)))
         if(isnan(vparall(ip)*vperpall(ip)))then
            print*,ipall,CIv(:,ip),lambda(ip),Ec(ip),weight(ip)
     $           ,vparall(ip),vperpall(ip)
            stop
         endif
      enddo
c      call MPI_REDUCE(weight,totweight,ntot/nprocs
c     $     ,MPI_DOUBLE_PRECISION,
c     $     MPI_SUM,0,comm,ierr)
c      if(me2.eq.0.and.startpos.lt.1)then
c         do ip=1,ntot/nprocs
c            totalweight=totalweight+totweight(ip)
c         enddo
c      endif
      
c      call MPI_BCAST(totalweight,1,MPI_DOUBLE_PRECISION,0,comm,ierr)  
      ndum=-1
      dummy=0.
      call moments(ndum,dummy,vparall,vperpall)
      if(nbinsvpar*nbinsvperp.ne.0)
     $     call distfun(ndum,dummy,vparall,vperpall)

      deallocate(Etemp,vparall,vperpall)!,totweight)

c------write positions.out of loading----------------
      if(abs(startpos).le.1)then
         allocate(wn(6,ntot/nprocs))
         allocate(wnall(6,ntot))
         do ip=1,ntot/nprocs
            wn(1:3,ip)=CIv(1:3,ip)
            wn(4,ip)=lambda(ip)
            wn(5,ip)=Ec(ip)/Qel
            if(weight(ip).eq.0.)then
               wn(6,ip)=0.
            else
               wn(6,ip)=weight(ip)!mu(ip)*b00/Ec(ip)
            endif
         enddo
         
         do i=1,6
            call MPI_GATHER(wn(i,:),ntot/nprocs,MPI_DOUBLE_PRECISION,
     $           wnall(i,:),ntot/nprocs
     $           ,MPI_DOUBLE_PRECISION,0,comm,ierr)
         enddo
         
         
         deallocate(wn)
         
         if(me2.eq.0)then
            open(unit=3,file='initial.out')
            do ip=1,ntot
               write(3,101)wnall(1:6,ip)
            enddo
            close(3)
            deallocate(wnall)
         endif
      endif
 101  format(6es24.15e3)        !,i12)
c--------------------------------------------------------------------
      else ! create look-up table for <s>
c--------------------------------------------------------------------  
      allocate(s(nbinsrad))
      allocate(sall(ntot/nprocs),thall(ntot/nprocs),xiall(ntot/nprocs))
      do i=1,nbinsrad
         s(i)=(i-0.5)/nsa
      enddo
      allocate(weight(ntot/nprocs))
      ip=0
      do ipall=me2+1,ntot,nprocs
         ip=ip+1
         weight(ip)=1.

         i=mod(ipall-1+npitch,npitch)+1
         x=random(ran_index(i),ran_array(:,i))
         x=smin+x*(smax-smin)   ! radially inside grid  ! s
         sall(ip)=x
         x=random(ran_index(i),ran_array(:,i))
         thall(ip)=x*twopi
         xiall(ip)=0.
         call interp3DRGC(sall(ip),thall(ip),xiall(ip),
     $        val,vals,pert)
         B=sqrt(val(1))

c         vmax=sqrt(2.*Tmax*Qel/mo)
         x1=random(ran_index(i),ran_array(:,i)) ! for Ep
         x2=random(ran_index(i),ran_array(:,i)) ! for lambda
         Ec(ip)=Es+x1*(Tmax*Qel-Es)
         lambda(ip)=2.*(x2-0.5)

         mu(ip)=Ec(ip)*(1.-lambda(ip)**2.)/B 
         dummy=lambda(ip)*sqrt(2.*Ec(ip)/mo)
         CIv(4,ip)=mo*dummy/(abs(Qe)*B*val(5))
      enddo
      
      write(fname,'(a4,i4.4)')'load',me2
      open(unit=12,file=fname,form='unformatted')
      write(12)sall
      write(12)thall
      write(12)lambda
      write(12)Ec
      close(12)
      
      CIv(1,:)=sall
      CIv(2,:)=thall
      CIv(3,:)=xiall

      deallocate(anis_s,Tpar_s)
      deallocate(sall,thall,xiall)

      endif

      
      end


c-------------------------------------------------------------------------
c#########################################################################
c-------------------------------------------------------------------------
      subroutine loadptcle(i,ip,fjacmax,gjacmax,yn,Ep,mu,lambda,wgt)
      use sinitial
      use terp
      use particle
      use dimarray
      use dstfun
      use interp
      use rand_var
      use machine
      use lkup
      use mpi_var
      use pival
      implicit none

      real yn(4),Ep,mu,lambda
      integer i,ip,check,is,check2
      real x,x1,x2,x3,x4,dsi,si2,Tloc,Tpar,Tperp
      real fjacmax,densround,jac,B,radcr,fjac,gjacmax,max,gjac
      real lamin,lamax,emx,emn,smx,smn
      real vmax,dstfn,wgt,vpar,vperp,densg


c      smn=1./nrad
c      smx=1.-1./nrad
      smn=smin
      smx=smax
      emx=Emax!*(1.-1./nen)
      emn=Es             ! must correspond to findmaxf
c      emx=Emax
c      lamin=1./nla
      lamin=0.
      lamax=2.-lamin

            
      check=0
      do while(check.eq.0)
         x=random(ran_index(i),ran_array(:,i))
         x=smn+x*(smx-smn)      ! radially inside grid  ! s
!         x=0.49  !!!!!CCCCCCCAAAAAAAAAAREFULLLLLLLLLLL
!         x=0.09  !!!!!CCCCCCCAAAAAAAAAAREFULLLLLLLLLLL
         is=floor(x*nsa)
         si2=is*ds+ds2
         dsi=(si2-x)/ds
         if(dsi.lt.0.)then
            is=is+1
            si2=is*ds+ds2
            dsi=(si2-x)/ds
         endif
         dsi=1.-dsi
      
         densround=dsi*dens_s(is+1)+(1.-dsi)*dens_s(is)
         Tloc=dsi*Tpar_s(is+1)+(1.-dsi)*Tpar_s(is)
         Tpar=Tloc*Qel          ! [J]
         Tperp=(dsi*anis_s(is+1)+(1.-dsi)*anis_s(is))*Tpar 
         
c     initial velocities
         x1=random(ran_index(i),ran_array(:,i)) ! for theta
!         x1=0.  !!!!!CCCCCCCAAAAAAAAAAREFULLLLLLLLLLL
         x2=random(ran_index(i),ran_array(:,i)) ! for lambda
         x3=random(ran_index(i),ran_array(:,i)) ! for Ep
         yn(2)=twopi*x1
         x1=random(ran_index(i),ran_array(:,i))  ! for phi
         yn(3)=twopi*x1
         lambda=lamin+x2*(lamax-lamin)-1.       ! lambda in (0,2)-1
c         Ep=Emin+x3*(emx-Emin)
         Ep=emn+(emx-emn)*x3
!         Ep=1.e4*Qel  !!!!!CCCCCCCAAAAAAAAAAREFULLLLLLLLLLL

c         if(Ep.gt.Emin)then
c            call interpS(lambda,Ep,yn(2),x,yn(1))
c         else
            yn(1)=x
c         endif
         if(yn(1).gt.0..and.yn(1).lt.smax)then
            call interp3DRGC(yn(1),yn(2),yn(3),
     $           val,vals,pert)
            jac=-(vals(6)*vals(2)-vals(5)*vals(1))/val(1)/val(5)
            B=sqrt(val(1))
            mu=Ep*(1.-lambda**2.)/B

            call computef(fjac,gjac,jac,Ep,densround,Tpar,Tperp,mu,Bc)
           
            if(gjac.gt.gjacmax)then
               print*,'gjac>gjacmax!', gjac,gjacmax
               print*,'Ep,lambda,s,th',Ep,lambda,yn(1),yn(2)/pi
               stop
            endif
            
            x4=random(ran_index(i),ran_array(:,i)) ! for f
            if(x4*gjacmax.lt.gjac)then
               wgt=fjac/gjac
               check=1
            endif
         endif
         
         
      enddo
      
      vpar=lambda*sqrt(2.*Ep/mo)
      yn(4)=mo*vpar/(abs(Qe)*B*val(5))

      if(ip.eq.ntot/nprocs)deallocate(loadarray)
                
      end subroutine loadptcle
c-------------------------------------------------------------------------
c#########################################################################
c-------------------------------------------------------------------------
      subroutine lookup
      use lkup
      use particle
      use interp
      use dstfun
      use dimarray
      use mpi_var
      use futils
      implicit none

      real Erange(2)
      integer fid,dims(4)
            
      call openf('loadarray.in',fid,'r')
      call getarr(fid,'/dims',dims)
      nla=dims(1)
      nen=dims(2)
      nrad=dims(3)
      npol=dims(4)
      if(me2.eq.0)!.and.runctrl.eq.-1)
     $     write(*,'(a,4i4)')' nla,nen,nrad,npol '
     $     ,nla,nen,nrad,npol

      allocate(loadarray(nla,nen,nrad,npol))
      call getarr(fid,'/Erange',Erange)
      call getarr(fid,'/st',loadarray)
      call closef(fid)

      if(Erange(1).gt.Es)then
         if(me2.eq.0)write(*,'(a,es8.1E1,a)')' Esl>Es'
     $        ,Erange(1)/Qel,'eV'
      endif
      if(Erange(2).lt.Emax)then
         if(me2.eq.0)write(*,'(a,es8.1E1,a)')' Eml<Emax'
     $        ,Erange(2)/Qel,'eV'
      endif
      Emin=Erange(1) !0.5*Emax/nen
      Emax=Tmax*Qel
c      Emin=0.
      

      end subroutine lookup
c-------------------------------------------------------------------------
c#########################################################################
c-------------------------------------------------------------------------
      subroutine interpS(lambda,Ep,theta,sa,so)
      use particle
      use dimarray
      use machine
      use dstfun
      use lkup
      implicit none

      real lambda,Ep,theta,sa,so
      real coef(2),twopi,thtmp
      real f22,f21,f12,f11
      real g111,g112,g121,g122,g211,g212,g221,g222
      real ds,ds2,si2,dsi,dla,dla2,lai2,dlai,den,den2,eni2,deni
     $     ,dth,dth2,thi2,dthi
      integer ila,ien,is,ith,j
      integer ila1,ien1,is1,ith1
      
      twopi=4.*acos(0.)
      lambda=lambda+1.
      thtmp=theta
      
      ila=floor(lambda/2.*nla)
      ien=floor((Ep-Es)/(Emax-Es)*nen)
      is=floor(sa*nrad)
      ith=floor(theta/twopi*npol)

      if(ila.eq.nla)ila=ila-1
      if(ien.eq.nen)ien=ien-1
      if(is.eq.nrad)is=is-1
      if(ila.eq.0)ila=1
      if(ien.eq.0)ien=1
      if(is.eq.0)is=1
      if(ith.eq.0)then ! theta periodic
         ith=npol
         theta=theta+twopi
      endif

      ds=1./nrad
      ds2=0.5*ds
      dth=twopi/npol
      dth2=0.5*dth
      dla=2./nla
      dla2=0.5*dla
      den=Emax/nen
      den2=0.5*den

      si2=ds*is+ds2
      thi2=dth*ith+dth2
      lai2=dla*ila+dla2
      eni2=den*ien+den2

      dsi=(si2-sa)/ds
      dthi=(thi2-theta)/dth
      dlai=(lai2-lambda)/dla
      deni=(eni2-Ep)/den
      if(dsi.lt.0.)then
         is=is+1
         si2=ds*is+ds2
         dsi=(si2-sa)/ds
      endif
      if(dthi.lt.0.)then
         ith=ith+1
         thi2=dth*ith+dth2
         dthi=(thi2-theta)/dth
      endif
      if(dlai.lt.0..and.ila.lt.nla-1)then
         ila=ila+1
         lai2=dla*ila+dla2
         dlai=(lai2-lambda)/dla
      endif
      if(deni.lt.0.)then
         ien=ien+1
         eni2=den*ien+den2
         deni=(eni2-Ep)/den
      endif
      dsi=1.-dsi
      dthi=1.-dthi
      dlai=1.-dlai
      deni=1.-deni

      is=min(is,nrad-1)
      is1=is+1
      if(ith.eq.npol+1)ith=1
      ith1=ith+1
      if(ith1.eq.npol+1)ith1=1
      ila1=ila+1
      ien=min(ien,nen-1)
      ien1=ien+1

c      if(is.ge.nrad.or.is.eq.0)print*,'sa',sa,is,is1
c      if(ila.ge.nla.or.ila.eq.0)print*,'lambda',lambda/2.,ila,ila1
c      if(ien.ge.nen.or.ien.eq.0)print*,'en',Ep/Emax,ien,ien1
      
      g111=loadarray(ila,ien,is,ith)*(1.-dthi)
     $     + loadarray(ila,ien,is,ith1)*dthi
      g112=loadarray(ila,ien,is1,ith)*(1.-dthi)
     $     + loadarray(ila,ien,is1,ith1)*dthi
      g121=loadarray(ila,ien1,is,ith)*(1.-dthi)
     $     + loadarray(ila,ien1,is,ith1)*dthi
      g211=loadarray(ila1,ien,is,ith)*(1.-dthi)
     $     + loadarray(ila1,ien,is,ith1)*dthi
      g122=loadarray(ila,ien1,is1,ith)*(1.-dthi)
     $     + loadarray(ila,ien1,is1,ith1)*dthi
      g212=loadarray(ila1,ien,is1,ith)*(1.-dthi)
     $     + loadarray(ila1,ien,is1,ith1)*dthi
      g221=loadarray(ila1,ien1,is,ith)*(1.-dthi)
     $     + loadarray(ila1,ien1,is,ith1)*dthi
      g222=loadarray(ila1,ien1,is1,ith)*(1.-dthi)
     $     + loadarray(ila1,ien1,is1,ith1)*dthi

      f11=g111*(1.-dsi) + g112*dsi
      f12=g121*(1.-dsi) + g122*dsi
      f21=g211*(1.-dsi) + g212*dsi
      f22=g221*(1.-dsi) + g222*dsi

      

c$$$      f11=loadarray(ila,ien,is)*(1.-dsi)
c$$$     $     + loadarray(ila,ien,is1)*dsi
c$$$      f12=loadarray(ila,ien1,is)*(1.-dsi)
c$$$     $     + loadarray(ila,ien1,is1)*dsi
c$$$      f21=loadarray(ila1,ien,is)*(1.-dsi)
c$$$     $     + loadarray(ila1,ien,is1)*dsi
c$$$      f22=loadarray(ila1,ien1,is)*(1.-dsi)
c$$$     $     +loadarray(ila1,ien1,is1)*dsi

      coef(1)=f11*(1.-dlai) + f12*dlai
      coef(2)=f21*(1.-dlai) + f22*dlai

      
      so=coef(1)*(1.-deni) + coef(2)*deni
     
         
      if(min(loadarray(ila,ien,is,ith),loadarray(ila,ien,is,ith1)
     $     ,loadarray(ila,ien,is1,ith),loadarray(ila,ien,is1,ith1)
     $     ,loadarray(ila,ien1,is,ith),loadarray(ila,ien1,is,ith1)
     $     ,loadarray(ila1,ien,is,ith),loadarray(ila1,ien,is,ith1)
     $     ,loadarray(ila,ien1,is1,ith),loadarray(ila,ien1,is1,ith1)
     $     ,loadarray(ila1,ien1,is,ith),loadarray(ila1,ien1,is,ith1)
     $     ,loadarray(ila1,ien1,is1,ith),loadarray(ila1,ien1,is1,ith1))
     $     .lt.0.)so=-2.      
         
      lambda=lambda-1.
      theta=thtmp
      


      end subroutine interpS
c-------------------------------------------------------------------------
c#########################################################################
c-------------------------------------------------------------------------
      subroutine interp3DcoefRGC2(s,theta,xi,val)!,ntheta)
      use particle
      use terp
      use arrays
      use sinitial
      implicit none

      integer  kxi,kxi1,is,jth,is1,jth1,l,kxim1,ntheta
      REAL  twopi,dxi,dth,xiperiode,thperiode,s,theta,xi,
     >    xik,xik1,si2,
     >     twopiPer,thetaj,dsi,dthj,f11,f12,f21,f22,xikm1,pxi,
     >     coef(4),br(3),val
c     $     ,var1(0:nsa,0:ntheta),var2(0:nsa,0:ntheta),


c      smax=1.00-0.5/nsa  

      twopi=4.*acos(0.)
      twopiPer=twopi/Lper

      if (s.lt.0.) then
         s=abs(s)
c         theta=theta+twopi/2.
      endif
      

      dxi=twopiPer/nxia
      dth=twopi/ntha

c      smax=starts+1.*nsat/nsa

      xiperiode=mod(mod(xi,twopiPer)+twopiPer,twopiPer) ! xi --> [0,2pi/L]
      thperiode=mod(mod(theta,twopi)+twopi,twopi)    ! theta --> [0,2pi]

      
c determination de xii et xi1 tq xi in [xii,xi1] 
      kxi=floor(xiperiode/dxi)+1
      if (kxi.eq.(nxia+1)) then
         kxi1=kxi
         kxim1=kxi-2
         kxi=kxi-1
      elseif (kxi.eq.1) then
         kxi1=kxi+1
         kxim1=nxia
c         kxi=kxi+1
      else
         kxi1=kxi+1
         kxim1=kxi-1
      endif

      xik=(kxi-1)*dxi
      xik1=(kxi1-1)*dxi
      xikm1=(kxim1-1)*dxi
      pxi=(xiperiode-xik)/dxi

c determination de la boite pour interpolation 2-D dans plan (s,u)
c c'eat-a-dire les (si2,thetaj) et DSi et DTHj
c$$$      if (s.gt.smax.or.(s-starts.lt.0.and.starts.ne.0)) then
c$$$c         print*,'PTCLE OUT OF BOUNDARY (INT3DRGC)'
c$$$c         print*,'s',s
c$$$         s=smax
c$$$      endif

      is=floor(s*nsa)

      si2=is*ds+ds2
      dsi=(si2-s)/ds
      if(dsi.lt.0.)then
         is=is+1
         si2=is*ds+ds2
         dsi=(si2-s)/ds
      endif      
      dsi=1-dsi


c controle de la taille du maillage (s-theta)
      is1=is+1
      if (is1-sinind.gt.nsa.or.(is-sinind.lt.1.and.starts.ne.0)) then
         if(is.lt.nsa)then
            print*,'s,is,is1',s,is-sinind,is1-sinind
            print*,'ADJUST NSAT TO HIGHER VALUE'
            stop
         else
            is=nsa-1
         endif
      endif

      jth=floor(thperiode/dth)+1
      jth1=jth+1
      if (jth.eq.ntha+1) then
         jth1=ntha+1
         jth=ntha
      endif

      thetaj=(jth-1)*dth
      dthj=(thperiode-thetaj)/dth
      f22=dsi*dthj
      f21=dsi-f22
      f12=dthj-f22
      f11=1+f22-(dsi+dthj)

c     if(abs(s-0.1355).lt.0.0005)print*,s,dsi,si2



      coef(1) = pxi*(pxi-1)/2.
      coef(2) = (1.0-pxi*pxi)
      coef(3) = pxi*(pxi+1)/2.

c interpolation 2-D lineaire sur les plan xik et xik1 et xikm1
      
      is=is-sinind
      is1=is+1


      l=1
      br(1)=
     >     f11 * tabarray(jth,kxi,is,l)+f22*tabarray(jth1,kxi,is1,l)
     >     +f21 *tabarray(jth,kxi,is1,l)+f12*tabarray(jth1,kxi,is,l)
      br(2)=
     >     f11*tabarray(jth,kxi1,is,l)+f22*tabarray(jth1,kxi1,is1,l)
     >     +f21*tabarray(jth,kxi1,is1,l)+f12*tabarray(jth1,kxi1,is,l)
      br(3)=
     >     f11*tabarray(jth,kxim1,is,l)+f22*tabarray(jth1,kxim1,is1,l)
     >     +f21*tabarray(jth,kxim1,is1,l)+f12*tabarray(jth1,kxim1,is,l)
c     interpolation 2-D lineaire entre xik et xik1 et xikm1

      val=coef(1)*br(3) + coef(2)*br(1) +
     >     coef(3)*br(2)


      end

c-------------------------------------------------------------------------
c#########################################################################
c-------------------------------------------------------------------------
      subroutine findmaxf(fmax,gmax)
      use machine
      use particle
      use terp
      use interp
      use dimarray
      use dstfun
      use lkup
      use mpi_var
      implicit none

      real fmax,gmax,fjac,gjac,fmaxall,gmaxall
      integer i,j,k,m,e,npolo,ntoro,nm,nengy
      real ss,thta,dens,Ep,mu,twopi,Tpar,Tperp
      real phi
      real jac,emn

      twopi=4.*acos(0.)
      emn=Es
      
      npolo=ceiling(100./nprocs)*nprocs
      ntoro=2*nxia
      nm=100
      nengy=1000

      fmax=0.
      gmax=0.

      do i=1,nsa
         ss=(i-0.5)/nsa
         dens=dens_s(i)
         Tpar=Tpar_s(i)*Qel
         Tperp=anis_s(i)*Tpar
         do j=1,npolo,nprocs
            thta=(j-0.5)/npolo*twopi
            do k=1,ntoro
               phi=(k-0.5)/ntoro*twopi/Lper
               call interp3DRGC(ss,thta,phi,val,vals,pert)
               jac=-(vals(6)*vals(2)-vals(5)*vals(1))/val(1)/val(5)
               do e=1,nengy
c     Ep=Emin+(e-1)*(Emax-Emin)/(nengy-1)
                  Ep=emn+(e-1)*(Emax-emn)/(nengy-1)
                  do m=1,nm
                     mu=(m-1)*Ep/sqrt(val(1))/(nm-1)
                     
                     call computef(fjac,gjac,jac,Ep,dens
     $                    ,Tpar,Tperp,mu,Bc)
                  
                     if(fjac.gt.fmax)fmax=fjac
                     if(gjac.gt.gmax)gmax=gjac
                  enddo
               enddo
            enddo
         enddo
      enddo

      call MPI_ALLREDUCE(fmax,fmaxall,1,MPI_DOUBLE_PRECISION,MPI_MAX
     $     ,comm,ierr)
      call MPI_ALLREDUCE(gmax,gmaxall,1,MPI_DOUBLE_PRECISION,MPI_MAX
     $     ,comm,ierr)


      fmax=2.*fmaxall
      gmax=2.*gmaxall

      end

c-------------------------------------------------------------------------
c#########################################################################
c-------------------------------------------------------------------------
      subroutine computef(fout,gout,jac,Ep,dens,Tpar,Tperp,mu,Bc)
      use load
      implicit none

      real jac,Ep,dens,Tpar,Tperp,mu,Bc
      real fout,gout
      real Tbg,Tpg

      
      fout=jac*sqrt(Ep)
     $     *dens
!     $     /(sqrt(Tpar)*Tperp) ! dens=n_c
     $     /(Tperp**1.5)      ! dens=N(s)
     $     *exp(-(mu*Bc/Tperp
     $     +abs(Ep-mu*Bc)/Tpar)) 

      Tbg=Tbc*Tpar
      Tpg=Tpc*Tperp

      gout=jac*sqrt(Ep)
     $     *dens
!     $     /(sqrt(Tbg)*Tpg) ! dens=n_c
     $     /(Tpg**1.5)      ! dens=N(s)
     $     *exp(-(mu*Bc/Tpg
     $     +abs(Ep-mu*Bc)/Tbg)) 

      end
