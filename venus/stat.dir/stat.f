      program stat

      implicit none

      integer ip,i,n,m,d,j,ntot,nbinsrad,nbinspol,nbinsvpar,nbinsvperp
      integer nsa,constbins
      real ttot,pi,twopi,Qel,mo,vparmax,vperpmax,dvperp,dvpar
      real,allocatable ::  yn(:,:),vpar(:),vperp(:),
     $     vparatbin(:),vperpatbin(:)
      real vperpmaxloc,vparmaxloc,intpol1,intpol2,norm
      integer vperpbin,vparbin
      integer,allocatable :: indxr(:),indxt(:)
      real,allocatable ::  distmatrix(:,:),
     $     Nvperp(:),Nvpar(:),
     $     distmatrixtot(:,:)
      real,allocatable ::  pperp(:),ppar(:)
     $     ,s(:),th(:),dens(:)
      real,allocatable ::  pperppol(:,:),pparpol(:,:)
     $     ,denspol(:,:)
     $     ,H(:,:),ph(:),polfac(:)
      real,allocatable :: TpoTb(:)
      real val(10),vals(8),pert(8),BoBc,ptherm,onembc
      real,allocatable :: thetaper(:),intervs(:),dsb(:)
      real,allocatable :: pchangle(:),lambda(:),Ec(:),weight(:)
      real Esum

      pi=2.*acos(0.)
      twopi=2.*pi
      Qel=1.6022e-19
      mo=1.6726e-27
      d=1
      constbins=0

      write(*,'("Number of particles: ")',
     $     advance='no')
      read(*,710) ntot
      write(*,'("Number of total radial grid points: ")',
     $     advance='no')
      read(*,710) nsa
       write(*,'("nbinsrad: ")',
     $     advance='no')
      read(*,710) nbinsrad
       write(*,'("nbinspol: ")',
     $     advance='no')
      read(*,710) nbinspol
       write(*,'("nbinsvpar: ")',
     $     advance='no')
      read(*,710) nbinsvpar
       write(*,'("nbinsvperp: ")',
     $     advance='no')
      read(*,710) nbinsvperp
 710  format(i9)

      allocate(yn(4,ntot),vpar(ntot),vperp(ntot))
      allocate(vparatbin(nbinsvpar),vperpatbin(nbinsvperp))
      allocate(indxr(ntot),indxt(ntot))
      allocate(distmatrix(nbinsvpar,nbinsvperp),
     $     Nvperp(nbinsvperp),Nvpar(nbinsvpar),
     $     distmatrixtot(nbinsvpar,nbinsvperp))
      allocate(pperp(nbinsrad),ppar(nbinsrad)
     $     ,s(nbinsrad),th(nbinspol),dens(nbinsrad))
      allocate(pperppol(nbinspol,nbinsrad),pparpol(nbinspol,nbinsrad)
     $     ,denspol(nbinspol,nbinsrad)
     $     ,H(nbinspol,nbinsrad),ph(nbinsrad),polfac(nbinsrad))
      allocate(TpoTb(nbinsrad)
     $     ,thetaper(ntot),intervs(0:nbinsrad),dsb(nbinsrad))
      allocate(pchangle(ntot),lambda(ntot),Ec(ntot),weight(ntot))


      open(unit=11,file='positions.out')
      write(39)nbinsvpar,nbinsvperp,nbinspol,nbinsrad,
     $     1,0.

      Esum=0.
      do ip=1,ntot
         read(11,101)yn(1,ip),yn(2,ip),yn(3,ip)
     $        ,lambda(ip),Ec(ip),weight(ip)
         Esum=Esum+Ec(ip)
         Ec(ip)=Ec(ip)*Qel
         vpar(ip)=sqrt(2.*Ec(ip)/mo)*lambda(ip)
         vperp(ip)=sqrt(abs(2.*Ec(ip)/mo-vpar(ip)*vpar(ip)))
c         if(weight(ip).ne.0.)write(88,'(2es14.6)')vpar(ip),vperp(ip)
      enddo
101   format(6es23.15e3)

      write(*,'(a,es9.1)')' <E> ',Esum/ntot

      open(unit=39,file='fort.39',form='unformatted')
      open(unit=59,file='fort.59',form='unformatted')

      intervs(0)=0.
      s(1)=real(1-0.5)/nsa
      intervs(1)=2.*s(1)
      dsb(1)=intervs(1)
c      V(1)=2.*pi*pi*r00*s(1)*a00*a00
      s(nbinsrad-1)=2.*s(1)
     $     +(1.-2.*s(1))*real(nbinsrad-2.5)/(nbinsrad-1)
c      s(nbinsrad)=real(sinind+nsa-0.5)/nsa
      s(nbinsrad)=2.*s(1)
     $     +(1.-2.*s(1))*real(nbinsrad-1.5)/(nbinsrad-1)
c      intervs(nbinsrad-1)=1.-1./nsa
      intervs(nbinsrad-1)=1.-1./(nbinsrad-1)
      intervs(nbinsrad)=1.
      dsb(nbinsrad)=1.-intervs(nbinsrad-1)
!$OMP PARALLEL
!$OMP DO
      do i=2,nbinsrad-2
         s(i)=2.*s(1)+(1.-2.*s(1))*real(i-1.5)/(nbinsrad-1)
c         V(i)=2.*pi*pi*r00*s(i)*a00*a00-V(i-1)
      enddo
!$OMP END DO
!$OMP DO
      do i=2,nbinsrad-2
         intervs(i)=0.5*(s(i+1)+s(i))
         dsb(i)=intervs(i)-intervs(i-1)
      enddo
!$OMP END DO
      dsb(nbinsrad-1)=intervs(nbinsrad-1)-intervs(nbinsrad-2)
!$OMP DO
      do i=1,nbinspol
         th(i)=mod(twopi*real(i-0.5)/nbinspol+3.*pi,twopi)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c       endif
      if(d.eq.1.or.constbins.eq.0)then   
         vperpmaxloc=vperp(1)
         vparmaxloc=abs(vpar(1))
         do i=2,ntot
            if(weight(i).ne.0.)then
            if(vperp(i).gt.vperpmaxloc) vperpmaxloc=vperp(i)
            if(abs(vpar(i)).gt.vparmaxloc)vparmaxloc=abs(vpar(i))
            endif
         enddo
         vperpmax=vperpmaxloc/4.
         vparmax=vparmaxloc/2.
         dvperp=vperpmax/nbinsvperp
         dvpar=2.*vparmax/nbinsvpar
      endif
      do i=1,nbinsvpar
         vparatbin(i)=-vparmax+real(i-0.5)*dvpar
      enddo
      do i=1,nbinsvperp
         vperpatbin(i)=real(i-0.5)*dvperp
      enddo
!$OMP PARALLEL
!$OMP DO
      do ip=1,ntot
         vpar(ip)=vpar(ip)+vparmax
c         vpar(ip)=max(vpar(ip),0.)   !new faster particles included in max bin
c         vpar(ip)=min(vpar(ip),2.*vparmax) !idem
c         vperp(ip)=min(vperp(ip),vperpmax) !idem
      enddo
!$OMP END DO
!$OMP END PARALLEL
                  
      write(39)vparatbin
      write(39)vperpatbin
!$OMP PARALLEL
!$OMP DO
      do ip=1,ntot
         do i=0,nbinsrad-1
            if(yn(1,ip).ge.intervs(i).and.yn(1,ip).lt.intervs(i+1))
     $           indxr(ip)=i+1
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      thetaper(:)=mod(mod(yn(2,:)+pi,twopi)+twopi,twopi)/twopi
      indxt(:)=floor(real(thetaper(:))*nbinspol)+1
      do i=1,nbinsrad
         pperp(i)=0.
         ppar(i)=0.
         dens(i)=0.
         do j=1,nbinspol
            pperppol(j,i)=0.
            pparpol(j,i)=0.
            denspol(j,i)=0.
!$OMP PARALLEL
!$OMP DO
            do n=1,nbinsvpar
               Nvpar(n)=0
            enddo
!$OMP ENDDO
!$OMP DO
            do m=1,nbinsvperp
               Nvperp(m)=0
            enddo
!$OMP ENDDO
!$OMP DO
            do m=1,nbinsvperp
               do n=1,nbinsvpar
                  distmatrix(n,m)=0.
               enddo
            enddo
!$OMP ENDDO
!$OMP DO
            do ip=1,ntot
               vparbin=0
               vperpbin=0
               if(indxr(ip).eq.i.and.indxt(ip).eq.j)then
                  vperpbin=floor(vperp(ip)/dvperp)+1
                  if(vperp(ip).ge.vperpmax)vperpbin=nbinsvperp
                  vparbin=floor(vpar(ip)/dvpar)+1
                  if(vpar(ip).ge.2.*vparmax)vparbin=nbinsvpar
                  if(vpar(ip).le.0.)vparbin=1
                  distmatrix(vparbin,vperpbin)=
     $                 distmatrix(vparbin,vperpbin)+weight(ip)
               endif
            enddo
!$OMP ENDDO
!$OMP END PARALLEL
            distmatrixtot=distmatrix
            norm=dvperp*dvpar*dsb(i)/nbinspol !d^2v*ds*dtheta
c     $              *0.5*a00**3.*sqrt(s(i))*cos(th(j)) !spatial jacobian
            write(39)distmatrixtot/norm !f(v_par,v_perp,s,theta)*vperpatbin(n)
            do n=1,nbinsvperp
               do m=1,nbinsvpar
                  Nvperp(n)=Nvperp(n)+distmatrixtot(m,n) !integrate over vpar
               enddo
               pperppol(j,i)=pperppol(j,i)+ !integrate over vperp
     $              vperpatbin(n)*vperpatbin(n)* Nvperp(n)
            enddo
            pperp(i)=pperp(i)+0.5*pperppol(j,i)
            do n=1,nbinsvpar
               do m=1,nbinsvperp
                  Nvpar(n)=Nvpar(n)+distmatrixtot(n,m) !integrate over vperp
                  denspol(j,i)=denspol(j,i)+distmatrixtot(n,m)
               enddo           
               pparpol(j,i)=pparpol(j,i)+ !integrate over vpar
     $              vparatbin(n)*vparatbin(n)* Nvpar(n)
            enddo
            write(12,'(5es14.6)')
     $           s(i),th(j)
     $           ,denspol(j,i)/(4.*pi*1.e-7*mo)
     $           ,pparpol(j,i)/(4.*pi*1.e-7)
     $           ,0.5*pperppol(j,i)/(4.*pi*1.e-7)

            dens(i)=dens(i)+denspol(j,i) !average over theta
            ppar(i)=ppar(i)+pparpol(j,i) !average over theta
         enddo
      enddo

      dens=dens/(4.*pi*1.e-7*mo) !correct for mu0 and mass/number density
      ppar=ppar/(4.*pi*1.e-7)
      pperp=pperp/(4.*pi*1.e-7)
      write(59)ttot
      write(59)s
      write(59)dens
      write(59)ppar
      write(59)pperp            !n(s),p_par(s),p_perp(s)
      
c      call inputvmec(s,pparpol,pperppol,TpoTb)
c      call inputleman(s,pparpol,TpoTb,denspol)



      end program stat
         
c-------------------------------------------


  
