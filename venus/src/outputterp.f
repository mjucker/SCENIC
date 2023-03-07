      subroutine outputterp

      use terpall
      use machine
      use mpi_var
      use scatter
      use particle
      implicit none
      

      integer i,l,jth,jxi,ir(2)
      REAL pas,theta,xi,twopi,Bsqare,volume
      real,allocatable :: Fct(:),rmn(:,:),zmn(:,:),phimn(:,:)
      real rdump,rdumm,zdum
      logical f_exist
      real,allocatable,dimension(:) :: pp_b,nth_o,TT_o

      twopi=4.*acos(0.)

      read(1,*) ni,lmnb
      call allocateEq
      allocate(rmn(lmnb,ni),zmn(lmnb,ni),phimn(lmnb,ni))
      allocate(pp(ni))

      
      do i=1,ni
         do l=1,lmnb
            read(1,*) mb(l),nm(l),
     >           rmn(l,i),zmn(l,i),phimn(l,i),bsmn(l,i), ! m,n,r_mn,z_mn,phi_mn,B^2_mn
     >           sigmn(l,i),taumn(l,i),dpdsmn(l,i),sigBsmn(l,i) ! sigma_mn,tau_,m,Dpperp/Ds-Dppar/Ds,sigmaBs
         enddo
      enddo
      rdump=0.
      rdumm=0.
      zdum=0.
      do l=1,lmnb
         rdump=rdump+rmn(l,ni)
         rdumm=rdumm+rmn(l,ni)*cos(0.5*twopi*mb(l))
         zdum=zdum+zmn(l,ni)*sin(0.25*twopi*mb(l))
      enddo
      
c      call extrap(bsmn,b00)
c      b00=sqrt(bsmn(1,1)) ! not exact: bsmn is at first half grid point!
      b00=1.5*bsmn(1,1)-0.5*bsmn(1,2)
      b00=sqrt(b00)
      r00=rmn(1,1) !R0 will be computed in outputleman.f
      a00=0.5*(0.5*(rdump-rdumm)+zdum)
c       print*,'nm(2)',nm(2)
      if(nm(2).ne.0.and.Lper.ne.nm(2))then
         print*,'CHECK Lper'
         stop
      endif

      
      INQUIRE(FILE='thermal.in', EXIST=f_exist)  ! thermal.in existing?
      if(f_exist)then
         allocate(pp_b(ni),nth_o(ni),TT_o(ni))
         open(unit=21,file='thermal.in')
      endif

      pas=1./ni
c VERSION MOI a verifier signe + pas
      do i=1,ni
         read(1,*) si(i),psip(i),chip(i),pp(i),
     >        ci(i),cj(i),cip(i),cjp(i) ! si,psi',chi',thermal pressure,muoI,muoJ,muoI',muoJ'
                                ! flux toro',flux polo',current polo,current toro
c         write(77,'(3es14.6)')si(i),ci(i),cj(i)
c     3D pas de Iplasma => J=0
c     cj(i)=0.
c     cjp(i)=0.
c     3D pas de Iplasma
c     si(i)=si(i)-pas/2.  ! on remet sur demi-point
c     on change de signe
c     psip(i)=-psip(i)
c            chip(i)=-chip(i)
c     ci(i)=-ci(i)
c     cj(i)=-cj(i)
         if(f_exist)then
            read(21,*)pp_b(i),nth_o(i),TT_o(i)
            pp(i)=pp_b(i)/(2.*twopi*1.e-07)
         else
            pp(i)=pp(i)/(2.*twopi*1.e-07)
         endif
      enddo
      cip(1)=cip(2)
      cjp(1)=cjp(2)
      
c     FIN LECTURE
      
      cip(ni)=cip(ni-1)
      cjp(ni)=cjp(ni-1)

      psi(1)=0.                 ! construction de psi,chi
      chi(1)=0.
      do i=2,ni
         psi(i)=psip(i)*pas+psi(i-1)
         chi(i)=chip(i)*pas+chi(i-1)
      enddo


      if(me2.eq.0)then
         write(*,'(1a5,2a4,3f7.3)')' b00,','r00,','a00 ',b00,r00,a00
c         write(*,'(a12,1es9.1)')' Tp on axis ',pp(1)/density/Qel
      endif

      if(f_exist)then
         deallocate(pp_b,nth_o,TT_o)
         close(21)
      endif
      deallocate(rmn,zmn,phimn)
      end


c---------------------------
      subroutine allocateEq
      use terpall
      use arrays
      use mpi_var
      implicit none

      nsa=ni
      allocate(mb(lmnb),nm(lmnb))
      allocate(si(ni),psi(ni),chi(ni),psip(ni),chip(ni),
     >     ci(ni),cj(ni),cip(ni),cjp(ni),bsmn(lmnb,ni),
     >     sigmn(lmnb,ni),taumn(lmnb,ni),dpdsmn(lmnb,ni),
     $     sigBsmn(lmnb,ni))
      allocate(tabarray(ntha+1,nxia+1,0:nsa,10),
     $     tabarrays(0:nsa,8))
c     $     pertarray(ntha+1,nxia+1,0:nsa,8))

      if(me2.eq.0)write(*,'(a,2i3)')' nsa,lmnb',nsa,lmnb

      end
