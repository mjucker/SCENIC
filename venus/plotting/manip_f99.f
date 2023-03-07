      program manip_f98
      implicit none

      integer nnstep,nstep,snstep,nskip
      real,allocatable,dimension(:) :: r,z,phi,rho  
      real,allocatable,dimension(:) :: n,s,th,ph
      integer t,it
      real twopi,rnstep,rsnstep,phioff
      

      open(unit=97,file='orbits.99')
      open(unit=96,file='RZ.99')
      open(unit=95,file='energy.99')
      twopi=8.*atan(1.)
c   Read initial stuff
      write(*,'("Number of time steps (..x.x.. and in thousands): ")',
     $     advance='no')
      read(*,710) rnstep
      write(*,'("Starting time step (..x.x.. and in thousands): ")',
     $     advance='no')
      read(*,710) rsnstep
      write(*,'("Reading interval : ")',
     $     advance='no')
      read(*,711) nskip
      write(*,'("Toroidal angle offset (in units of pi): ")',
     $     advance='no')
      read(*,710) phioff
 710  format(f9.5)
 711  format(i8)
      nnstep=nint(rnstep*1000)
      snstep=nint(rsnstep*1000)
      nstep=nnstep+snstep

      phioff=phioff*twopi/2.
      
      allocate(r(nstep),z(nstep),phi(nstep),
     $     n(nstep),s(nstep),th(nstep),ph(nstep),rho(nstep))

      do t=1,nstep
         read(99,12)n(t),s(t),th(t),ph(t),rho(t)
      enddo
      call boozer2RZ(nstep,n,s,th,ph,r,z,phi)
      
      write(95,'(a)')'time   energy'
      it=0
      do t=snstep,nstep,nskip
         it=it+1
         phi(t)=phi(t)+phioff
         write(97,11)r(t)*cos(phi(t)),r(t)*sin(phi(t)),z(t),rho(t),it
         write(96,13)r(t),z(t),rho(t),it
         write(95,13)n(t),rho(t)
      enddo
      
 10   format(1p3e16.8)
 11   format(1p4e16.8,i6)
 12   format(1p5e16.8)!i9,1p4e16.8)
 13   format(1p3e16.8,i6)
      end
c----------------------------------------------------------------------c
      subroutine boozer2RZ(nstep,n,s,thetab,zeta,r,z,phi)
c
!
C.. Implicits ..
!      implicit none
c
      real :: r(nstep)  , z(nstep)  , phi(nstep)
c
      real, allocatable,dimension(:,:) :: rmn , zmn , phmn,
     $     sigmb,taub,fperp,fparp
      real, allocatable,dimension(:,:) :: rmnv, zmnv, phmnv
      real :: s(nstep),thetab(nstep),zeta(nstep),time(nstep)
      real, allocatable,dimension(:,:) :: tcos, tsin
      integer, allocatable,dimension(:)    :: xm  , xn
      integer :: x(nstep)
      
!
!
C.. Internal scalars
      integer ::  mnt, mn, iarg
      integer :: t,i
      real :: twopi, pi, arg, indbef,indaf,rnstep

c   Read initial stuff
      read (37,711) ns,mnt
 711  format(1x,2i6)
c   Allocate
C
      allocate(rmn(mnt,0:ns),zmn(mnt,0:ns),phmn(mnt,0:ns),
     $     sigmb(mnt,ns),taub(mnt,ns),
     $     fperp(mnt,nstep),fparp(mnt,nstep))
      allocate(rmnv(mnt,nstep),zmnv(mnt,nstep),phmnv(mnt,nstep))
c      allocate(time(nstep),s(nstep),thetab(nstep),zeta(nstep),
c     $     rhop(nstep),r(nstep),z(nstep),phi(nstep),x(nstep))
      allocate(xm(mnt),xn(mnt))
      allocate(tcos(nstep,mnt),tsin(nstep,mnt))
c
c   Calculations
      twopi = 8. * atan(1.) 
      pi=0.5*twopi
c
      do 10 i = 1,ns
        do 8 mn = 1,mnt
          read (37,810) xm(mn),xn(mn),rmn(mn,i),zmn(mn,i),phmn(mn,i),
     $          sigmb(mn,i),taub(mn,i),fperp(mn,i),fparp(mn,i)
 8      end do
 10   end do
 810  format (1x,2i8,1p7e14.6)
c
c   Extrapolation to magnetic axis from half integer meshpoints
      rmn(1,0)=1.5*rmn(1,1)-0.5*rmn(1,2)
      zmn(1,0)=1.5*zmn(1,1)-0.5*zmn(1,2)
      phmn(1,0)=1.5*phmn(1,1)-0.5*phmn(1,2)
      do 30 mn=2,mnt
         rmn(mn,0)=0.
         zmn(mn,0)=0.
         phmn(mn,0)=0.
 30   end do

      do 11 t=1,nstep
         x(t)=floor(s(t)*ns)
 11   end do
      
      
C   Interpolation of fort.99 orbit values onto the fort.37 equilibrium grid
      do 12 t=1,nstep
         do 13 mn=1,mnt
            indbef=s(t)*ns-x(t)
            indaf=x(t)+1.-s(t)*ns
            rmnv(mn,t)=indbef*rmn(mn,x(t)+1)
     $           +indaf*rmn(mn,x(t))
            zmnv(mn,t)=indbef*zmn(mn,x(t)+1)
     $           +indaf*zmn(mn,x(t))
            phmnv(mn,t)=indbef*phmn(mn,x(t)+1)
     $           +indaf*phmn(mn,x(t))
 13      end do
 12   end do
      
C   Preparing conversion into real space    
      do 25 mn = 1,mnt
         do 22 t = 1,nstep
             arg = xm(mn) * thetab(t)  -  xn(mn) * zeta(t)
c             iarg = arg
c             arg = twopi * (arg - iarg)
             tcos(t,mn) = cos(arg)
             tsin(t,mn) = sin(arg)
 22      end do
 25   end do
c
C Conversion into R,Z,phi
      do 50 t = 1,nstep
            r(t) = 0.
            z(t) = 0.
            phi(t) = 0.
c
        do 45 mn = 1,mnt
              r(t) = r(t)+rmnv(mn,t)*tcos(t,mn) !!!!!!!!!!!!!!
              z(t) = z(t)+zmnv(mn,t)*tsin(t,mn) !!!!!!!!!!!!!!
              phi(t)= phi(t)+phmnv(mn,t)*tsin(t,mn)
 45     end do
        phi(t)=phi(t)+zeta(t)
 50   end do
c
      end subroutine boozer2RZ
C----------------------------------------------------------------------C
