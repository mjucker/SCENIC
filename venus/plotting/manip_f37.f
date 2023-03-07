      program manip_f37
      implicit none

      integer ni,nm
      integer,allocatable :: m(:),n(:)
      integer i,j,k,t,np,nt,nstep
      real,allocatable :: f(:,:,:),g(:,:,:),h(:,:)
      real phi,twopi,theta
      
      open(unit=38,file='geometry.37')
      open(unit=39,file='coils.37')
      open(unit=40,file='modB.37')
      twopi=8.*atan(1.)
      write(*,'("Number of poloidal grid points: ")',
     $     advance='no')
      read(*,710) np
      write(*,'("Number of toroidal grid points points: ")',
     $     advance='no')
      read(*,710) nt      
      write(*,'("Step size for radial points: ")',
     $     advance='no')
      read(*,710) nstep      
 710  format(i4)

      read(37,*)ni,nm
      allocate(f(nm,ni,8),g(np,ni,3),h(np,2),m(nm),n(nm))

      do i=1,ni
         do j=1,nm
            read(37,10)m(j),n(j),f(j,i,1:8)
         enddo
      enddo

    ! Write Coils 
      write(39,11)np
      do i=1,2
         do j=1,np
            h(j,i)=0.
         enddo
      enddo
      do j=1,np
         theta=(twopi+0.001)/(np-1)*(j-1)
         do i=1,nm
            if(i.eq.2)then
               h(j,1)=h(j,1)+
     $              (f(i,ni,1)+f(i,ni,1)/10.)*
     $              cos(m(i)*theta)
               h(j,2)=h(j,2)+
     $              (f(i,ni,2)+f(i,ni,2)/10.)*
     $              sin(m(i)*theta)
            else
               h(j,1)=h(j,1)+
     $              (f(i,ni,1))*cos(m(i)*theta)
               h(j,2)=h(j,2)+
     $              (f(i,ni,2))*sin(m(i)*theta)
            endif
         enddo
         write(39,12)
     $        h(j,1),
     $        h(j,2),
     $        phi,h(j,2)
      enddo

      call convert(m,n,f,g,nm,np,ni)

    ! Write LCFS
      write(38,14)nt,np
      do j=1,np
         do k=1,nt
            phi=twopi/(nt-1)*(k-1)
            write(38,13)
     $           g(j,ni,1)*cos(phi),g(j,ni,1)*sin(phi),g(j,ni,2),
     $           k
         enddo
      enddo

    ! Write ModB in (R,Z)
      write(40,15)np,ni/nstep
      do i=1,ni,nstep
         do j=1,np
            write(40,16)g(j,i,1),g(j,i,2),g(j,i,3)
         enddo
      enddo


 10   format(1x,2i8,1p8e14.6)
 11   format(i5)
 12   format(1p4e14.6)!,i8)
 13   format(1p3e14.6,i8)
 14   format(2i5)
 15   format(3i5)
 16   format(1p3e14.6)
      end

c-------------------------------------------
      subroutine convert(m,n,f,g,nm,np,ni)
      implicit none

      real f(nm,ni,8),g(np,ni,3)
      real twopi
      real theta,dth
      integer i,j,k,t,np,ni
      integer nm,m(nm),n(nm)

      twopi=8.*atan(1.)

    ! Conversion for wall, but including all flux surfaces
    ! g(:,:,1)=R(theta,r),g(:,:,2)=Z(theta,r),g(:,:,3)=B(theta,r)
      dth=twopi/(np-1)
      do t=1,3
         do i=1,ni
            do j=1,np
               g(j,i,t)=0.
            enddo
         enddo
      enddo
      do i=1,ni
         do j=1,np
            theta=dth*(j-1)
            do k=1,nm
               g(j,i,1)=g(j,i,1)+f(k,i,1)*cos(m(k)*theta)
               g(j,i,2)=g(j,i,2)+f(k,i,2)*sin(m(k)*theta)
               g(j,i,3)=g(j,i,3)+f(k,i,4)*cos(m(k)*theta)
            enddo
            g(j,i,3)=sqrt(g(j,i,3))
         enddo
      enddo


      end
      


