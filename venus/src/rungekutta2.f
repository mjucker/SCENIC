      subroutine rungekutta2(past,yn,ynt,frk2,mu,
     >     Btilt)
      use particle
      use dimarray
      implicit none

c      include 'include.dir/dimarray.inc'
c      include 'include.dir/particule.inc'

      integer i,N
      parameter (N=4)
      REAL past,past2,
     >     mu,Btilt,yn(N),ynt(N),k1(N),
     >     k2(N),f(N),frk2(N)
      real chk

c INPUT   yn a t=t
c OUTPUT  yn a t=t+past
c         ynt c'est yn a t=t
      past2=past/2.

       
      call systEDO(mu,Btilt,yn,f) ! f(to,yo)
     

      do i=1,N
         ynt(i)=yn(i)           ! save yn(tn)
         k1(i)=past*f(i)
         yn(i)=ynt(i)+k1(i)/2.
         frk2(i)=f(i)                      ! save f en to pour rk4mixte.f
      enddo 
c$$$      yn(1)=abs(yn(1))
c$$$      yn(1)=min(yn(1),1.0)

      call systEDO(mu,Btilt,yn,f) ! f(to+h/2,yo+k1/2)

      do i=1,N
         k2(i)=past*f(i)
         yn(i)=ynt(i)+k2(i)     ! y(n) -->y(n+1) 
      enddo
c$$$      yn(1)=abs(yn(1))
c$$$      yn(1)=min(yn(1),1.0)


      end
