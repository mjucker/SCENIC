      subroutine rungekutta4mixte(past,yn,ynt,frk2,
     >     Ep,mu,Btilt,signerho,smax3,check)
      use particle
      use interp
      use mpi_var !debugging only
      implicit none

      
      integer N,i,NN,check
      parameter (N=4,NN=3)
      REAL past,past2,Ep,Btilt,
     >     twomo,signerho,modB,interp3DB,
     >     mu,ppar,
     >     yn(N),ynt(N),f(N),frk2(N),
     >     k1(NN),k2(NN),k3(NN),k4(NN)
      real smax3,dummy

      past2=past/2.
      twomo=2.*mo
c INPUT   yn a t=t
c OUTPUT  yn a t=t+past
c         ynt c'est yn a t=t

c------------------------------------------------------
c FIRST STEP
c------------------------------------------------------
      do i=1,NN
         ynt(i)=yn(i)
         k1(i)=past*frk2(i)
         yn(i)=ynt(i)+k1(i)/2.
      enddo
      
      if(yn(1).lt.0.)then
         check=2
         return
      endif
      call interp3DRGC(yn(1),yn(2),yn(3),val,vals,pert)
      modB=sqrt(val(1))
      yn(4)=twomo*(Ep-mu*modB)/(Qe*Qe*val(1)*val(5)*val(5)) !rho^2
      if(yn(4).lt.0.)then
         check=1
         return
      endif
      yn(4)=sign(sqrt(yn(4)),signerho)

c------------------------------------------------------
c SECOND STEP
c------------------------------------------------------
      call systEDOmixte(mu,Btilt,yn,f) ! f(to+h/2,yo+k1/2)
      do i=1,NN
         k2(i)=past*f(i)
         yn(i)=ynt(i)+k2(i)/2. 
      enddo

      if(yn(1).lt.0.)then
         check=2
         return
      endif      
      call interp3DRGC(yn(1),yn(2),yn(3),val,vals,pert)
      modB=sqrt(val(1))
      yn(4)=twomo*(Ep-mu*modB)/(Qe*Qe*val(1)*val(5)*val(5)) !rho^2
      if(yn(4).lt.0.)then
         check=1
         return
      endif 
      yn(4)=sign(sqrt(yn(4)),signerho)

c------------------------------------------------------
c THIRD STEP
c------------------------------------------------------
      call systEDOmixte(mu,Btilt,yn,f) ! f(to+h/2,yo+k2/2)
      do i=1,NN    
         k3(i)=past*f(i)
         yn(i)=ynt(i)+k3(i)    
      enddo

      if(yn(1).lt.0.)then
         check=2
         return
      endif
      call interp3DRGC(yn(1),yn(2),yn(3),val,vals,pert)
      modB=sqrt(val(1))
      yn(4)=twomo*(Ep-mu*modB)/(Qe*Qe*val(1)*val(5)*val(5)) !rho^2
      if(yn(4).lt.0.)then
         check=1
         return
      endif
      yn(4)=sign(sqrt(yn(4)),signerho)

c------------------------------------------------------
c FOURTH STEP
c------------------------------------------------------
      call systEDOmixte(mu,Btilt,yn,f) ! f(to+h,yo+k3)
      do i=1,NN
         k4(i)=past*f(i)     
      enddo

c------------------------------------------------------
c FINAL STEP
c------------------------------------------------------
      do i=1,NN
         yn(i) = ynt(i) + 1./6.*(k1(i)+2.*k2(i)+2.*k3(i)+k4(i)) ! y(n)-->y(n+1)
      enddo

      if(yn(1).lt.0.)then
         check=2
         return
      endif
      call interp3DRGC(yn(1),yn(2),yn(3),val,vals,pert)
      modB=sqrt(val(1))
      yn(4)=twomo*(Ep-mu*modB)/(Qe*Qe*val(1)*val(5)*val(5)) !rho^2
      if(yn(4).lt.0.)then
         check=1
         return
      endif
      yn(4)=sign(sqrt(yn(4)),signerho)


      end
