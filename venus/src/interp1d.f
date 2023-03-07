      subroutine interp1d(xi,fxi,nx,x,fx,eps)
      implicit none

      integer ns,i,m,nx,nxmax
      parameter (nxmax=15)
      REAL dif,dift,ho,hp,w,den,eps,x,fx,xi(nx),fxi(nx),
     >     c(nxmax),d(nxmax)

      ns=1
      dif=abs(x-xi(1))
      do i=1,nx
         dift=abs(x-xi(i))
         if (dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=fxi(i)
         d(i)=fxi(i)
      enddo
      fx=fxi(ns)
      ns=ns-1
      do m=1,nx-1
         do i=1,nx-m
            ho=xi(i)-x
            hp=xi(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         enddo
         if (2*ns.lt.nx-m) then
            eps=c(ns+1)
         else
            eps=d(ns)
            ns=ns-1
         endif
         fx=fx+eps
      enddo

      end
