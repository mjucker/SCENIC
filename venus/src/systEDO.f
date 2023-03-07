      subroutine systEDO(mu,Btilt,x,f)
      use particle
      use interp
      implicit none

c      include 'include.dir/dimarray.inc'
c      include 'include.dir/particule.inc'


      integer N
      parameter (N=4)
      REAL rho,D,mu,Btilt,
     >     x(N),f(N),coef(8)
      integer i

c update 04.07.07 with sigmaBs

c********************
c x(1) --->s        *
c x(2) --->theta    *
c x(3) --->xi       *
c x(4) --->rho      *
c t=tau             *
c with t=m/eB *tau  *
c                   *
c Ds/Dtau     = f(1)*
c Dtheta/Dtau = f(2)*
c Dxi/Dtau    = f(3)*
c Drho/Dtau   = f(4)*
c********************

      call interp3DRGC(x(1),x(2),x(3),val,vals,pert)

      rho=x(4)

      D=vals(6)*vals(2) - vals(5)*vals(1) + (rho + pert(1))*
     >     (vals(2)*vals(3) - vals(1)*vals(4))
c   add sigmaBs terms
     $     + (rho + pert(1))*(vals(1)*val(9) + vals(2)*val(10))

      coef(5)=1.0/D/Btilt      ! 1/(D Btilt)
      coef(6)=val(5)*val(6)   ! sigma*tau
      coef(1)=val(5)*val(5)*val(1)*rho/gamma ! sigma^2 B^2 rho/(gamma)
      coef(2)=1./gamma*(mo*mu/Qe/Qe + 
     >     coef(6)*sqrt(val(1))*rho*rho)
      coef(3)=vals(1)*coef(5)    ! muoI/D/Btilt
      coef(4)=vals(2)*coef(5)    ! muoJ/D/Btilt
      coef(7)=val(1)*val(5)*rho*rho/gamma !sigma*B^2 rho^2/(gamma)
      coef(8)=mo/Qe



      f(1) = coef(3)*(coef(8)*pert(6) + coef(2)*val(3)
     $     - coef(1)*pert(2))
     >     + coef(4)*(coef(8)*pert(7) + coef(2)*val(4)
     $     - coef(1)*pert(3))

 
      f(2) = coef(1)*coef(5)*(vals(6)
     $     + (rho+pert(1))*vals(3)+vals(1)*pert(4))
     >     - coef(3)*(coef(8)*pert(5) + coef(2)*val(2)
     >     + coef(7)*val(7)) ! anisotropy
c   add sigmaBs terms
     $     + coef(5)*(coef(1)*(rho + pert(1))*val(10) !sigmaBs terms
     $     - val(8)*(coef(8)*pert(7)+coef(2)*val(4)-coef(1)*pert(3)))


      f(3) = coef(1)*coef(5)*(vals(5)
     $     + (rho+pert(1))*vals(4)+vals(2)*pert(4))
     >     - coef(4)*(coef(8)*pert(5) + coef(2)*val(2)
     >     + coef(7)*val(7)) ! anisotropy
c   add sigmaBs terms
     $     + coef(5)*(-coef(1)*(rho + pert(1))*val(9) !sigmaBs terms
     $     + val(8)*(coef(8)*pert(6) + coef(2)*val(3)-coef(1)*pert(2)))


      f(4) = -pert(8)
     $     -coef(5)*(vals(6)
     $     +(rho+pert(1))*vals(3)+vals(1)*pert(4))
     >     * (coef(8)*pert(6)+coef(2)*val(3))
     >     - coef(5)*(vals(5)
     $     + (rho+pert(1))*vals(4)+vals(2)*pert(4))
     >     * (coef(8)*pert(7)+coef(2)*val(4))
     >     + coef(5)*(vals(1)*pert(2) + vals(2)*pert(3))
     $     * (coef(8)*pert(5)+coef(2)*val(2)
     >     + coef(7)*val(7)) !anisotropy
cc   add sigmaBs terms
     $     - coef(5)*((rho + pert(1))*val(10) + val(8)*pert(3))*   !sigmaBs terms
     $     (coef(8)*pert(6) + coef(2)*val(3))
     $     + coef(5)*((rho + pert(1))*val(9) + val(8)*pert(2))* 
     $     (coef(8)*pert(7) + coef(2)*val(4))



      end
