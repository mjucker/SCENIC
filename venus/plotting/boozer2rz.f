c--------0---------0---------0---------0---------0---------0---------0-c
c---TRANSFORMS VMECAP OUTPUT FROM CYLINDRICAL COORDINATE            ---c
c---INFORMATION TO CARTESIAN COORDINATES. THE OUTPUT CORRESPONDS TO ---c
c---THE CARTESIAN COORDINATES OF EACH GRID POINT OF A FLUX SURFACE. ---c
c---'b' CONTAINS mod-B, p CONTAINS p_par and u CONTAINS p_perp      ---c
c---'d' CONTAINS THE HOT PARTICLE DENSITY                           ---c
c---            Last Modification:                         ---c
c----------------------------------------------------------------------c
      program boozer2RZ
c
!
C.. Implicits ..
!      implicit none
c
      real, allocatable, dimension(:) :: r  , z  , ph
c
      real, allocatable,dimension(:,:) :: rmn , zmn , phmn,
     $     bs,sigmb,taub,dsigds,sigBs
      real, allocatable,dimension(:,:) :: rmnv, zmnv, phmnv
      real, allocatable,dimension(:)     :: s,thetab,zeta,rhop,time
      real, allocatable,dimension(:,:) :: tcos, tsin
      integer, allocatable,dimension(:)    :: xm  , xn
      integer, allocatable,dimension(:)  :: x
      
!
!
C.. Internal scalars
      integer ::  mnt, mn, iarg
      integer :: t,tt,i,nstep,nstart
      real :: twopi, pi, arg, indbef,indaf,rnstep,nrstart,sshift

c   Read initial stuff
      write(*,'("Number of time steps (..x.x.. and in thousands): ")',
     $     advance='no')
      read(*,710) rnstep
      write(*,'("Starting time step (..x.x.. and in thousands): ")',
     $     advance='no')
      read(*,710) nrstart
      write(*,'("Shift radial position (in terms of s): ")',
     $     advance='no')
      read(*,710) sshift
      nrstart=max(0.001,nrstart)
      read (37,711) ns,mnt
 710  format(f9.5)
 711  format(1x,2i6)
      nstep=nint(rnstep*1000)
      nstart=nint(nrstart*1000)-1
c   Allocate
C
      allocate(rmn(mnt,0:ns),zmn(mnt,0:ns),phmn(mnt,0:ns),
     $     bs(mnt,ns),sigmb(mnt,ns),taub(mnt,ns),
     $     dsigds(mnt,nstep),sigBs(mnt,nstep))
      allocate(rmnv(mnt,nstep+nstart),zmnv(mnt,nstep+nstart),
     $     phmnv(mnt,nstep+nstart))
      allocate(time(nstep+nstart),s(nstep+nstart),thetab(nstep+nstart)
     $     ,zeta(nstep+nstart),rhop(nstep+nstart),x(nstep+nstart)
     $     ,r(nstep),z(nstep),ph(nstep))
      allocate(xm(mnt),xn(mnt))
      allocate(tcos(nstep+nstart,mnt),tsin(nstep+nstart,mnt))
c
c   Calculations
      twopi = 8. * atan(1.) 
      pi=0.5*twopi
c
      do 10 i = 1,ns
        do 8 mn = 1,mnt
          read (37,810) xm(mn),xn(mn),rmn(mn,i),zmn(mn,i),phmn(mn,i),
     $          bs(mn,i),sigmb(mn,i),taub(mn,i),dsigds(mn,i),
     $          sigBs(mn,i)
 8      end do
 10   end do
 810  format (1x,2i8,1p8e14.6)
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
      
      do 11 tt=1,nstart+nstep
         read(99,812)time(tt),s(tt),thetab(tt),zeta(tt),rhop(tt)
         s(tt)=s(tt)+sshift
         if(s(tt).gt.1.)then
            print*,'SSHIFT TOO LARGE, s(tt)',s(tt)
            stop
         endif
         x(tt)=floor(s(tt)*ns)
 11   end do
 812  format (1p5e16.8)
      
      
C   Interpolation of fort.99 orbit values onto the fort.37 equilibrium grid
      do 12 tt=1,nstep+nstart
         do 13 mn=1,mnt
            indbef=s(tt)*ns-x(tt)
            indaf=x(tt)+1.-s(tt)*ns
            rmnv(mn,tt)=indbef*rmn(mn,x(tt)+1)
     $           +indaf*rmn(mn,x(tt))
            zmnv(mn,tt)=indbef*zmn(mn,x(tt)+1)
     $           +indaf*zmn(mn,x(tt))
            phmnv(mn,tt)=indbef*phmn(mn,x(tt)+1)
     $           +indaf*phmn(mn,x(tt))
 13      end do
 12   end do
      
C   Preparing conversion into real space    
      do 25 mn = 1,mnt
         do 22 tt = 1,nstep+nstart
             arg = xm(mn) * thetab(tt)  -  xn(mn) * zeta(tt)
c             iarg = arg
c             arg = twopi * (arg - iarg)
             tcos(tt,mn) = cos(arg)
             tsin(tt,mn) = sin(arg)
 22      end do
 25   end do
c
C Conversion into R,Z,phi
      do 50 t = 1,nstep!-nstart
            r(t) = 0.
            z(t) = 0.
            ph(t) = 0.
c
            tt=nstart+t
        do 45 mn = 1,mnt
              r(t) = r(t)+rmnv(mn,tt)*tcos(tt,mn) !!!!!!!!!!!!!!
              z(t) = z(t)+zmnv(mn,tt)*tsin(tt,mn) !!!!!!!!!!!!!!
              ph(t)= ph(t)+phmnv(mn,tt)*tsin(tt,mn)
 45     end do
        ph(t)=ph(t)+zeta(tt)
        write (98,1008) r(t),z(t),ph(t)
 50   end do
 1008 format(1p3e16.8)
c
C----------------------------------------------------------------------C

      end program boozer2RZ
C-----------------------------------------------------------------------
C  
C----------------------------------------------------------------------C
