      subroutine rungekuttamixte(past,yn,ynt,Ep,mu,
     >     Btilt,rholim,smax3,check)
      use particle
      use mpi_var !debugging only
      implicit none


      integer N,i,j,reg,inreg,check
      parameter (N=4)
      REAL past,past0,past1,Ep,mu,Btilt,B,rholim,
     >     signerho,frk2(N),yn(N),ynt(N),yntt(N),yn0(N),
     >     interp3DB
      real smax3,const

c INPUT   yn a t=t
c OUTPUT  yn a t=t+past
c         ynt c'est yn a t=t
      check=0 

      call rungekutta2(past,yn,ynt,frk2,mu,Btilt) ! t_n-->t_n+1 avec RK2 with N=4

c particule ne change pas de signe v// ==> RK4 with N=3
      if((yn(4)*ynt(4).ge.0.).and.abs(yn(4)).gt.rholim)then
         do i=1,N
            yntt(i)=ynt(i)
         enddo

         signerho=sign(1.,yn(4))
         call rungekutta4mixte(past,ynt,yntt,frk2,Ep,
     >        mu,Btilt,signerho,smax3,check) ! t_n-->t_n+1 avec RK4 with N=3

c         rholim=4.e-5 ! same as in solverRK
         if(check.eq.1)then
            if(yn(1).lt.0.)then
               yn(1)=abs(yn(1))!0.
               yn(2)=-yn(2)
            endif
            ynt=yn                  !problem in RK4, take RK2 result
         elseif(check.eq.2)then     !s<0 in RK4
c$$$            ynt(1)=abs(yn(1))!0.
c$$$            ynt(2)=-ynt(2)
            ynt=yntt
            ynt(2)=-yntt(2)
         endif
 
         do i=1,N
            yn(i)=ynt(i)
            ynt(i)=yntt(i)
         enddo
      else
         check=1
         if(yn(1).lt.0.)then
            yn(1)=-yn(1)!0.
            yn(2)=-yn(2)
         endif
      endif

      
      end
