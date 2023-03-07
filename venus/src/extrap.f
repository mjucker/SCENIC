      subroutine extrap(bmn,b00)
      
      use terpall
      implicit none


      REAL bmn(lmnb,ni),rmn(lmnb,ni)
      REAL b00,r00
      REAL dth,dph,twopi,arg,iarg
      integer nj,nk,nst,ks
      parameter (nj=98,nk=1)
      REAL tcos(lmnb,nj*nk),cosp
      REAL bt(nj*nk,ni),bax(ni),b(nj*nk*Lper,ni,nj)
      REAL rt(nj*nk,ni),rax(ni),r(nj*nk*Lper,ni,nj)
      integer mn,j,k,lk,njk,kf,dt,jk,i

      intrinsic atan,cos,sin
c conversion into real space of first component
      njk=nj*nk
      twopi=8.*atan(1.)/Lper
      dth=1./nj   
      dph=1./(nk*Lper)
      

      do 20 mn = 1,lmnb
         do 21 k=1,nk
            do 22 j = 1,nj
               jk = j + (k-1) * nj
               arg = mb(mn) * (j-1) * dth  -  nm(mn) * (k-1) * dph
               iarg = 0.!arg
               arg = twopi * (arg - iarg)
               tcos(mn,jk) = cos(arg)
c               tsin(jk) = sin(arg)
 22         end do
 21      end do
 20   end do
      do 40 i=1,ni
         do lk=1,njk
            bt(lk,i)=0.
c           rt(lk,i)=0.
         end do
         do 30 lk=1,njk
            do 31 mn=1,lmnb
               bt(lk,i)=bt(lk,i)+bmn(mn,i)*tcos(mn,lk)
c              rt(lk,i)=rt(lk,i)+rmn(mn,i)*tcos(mn,lk)
 31         end do
 30      end do
 40   end do
      do i=1,ni
         do nst=1,Lper
            do k=1,nk
               ks=(nst-1)*nk+k
               do j=1,nj
                  jk=(k-1)*nj+j
                  b(ks,i,j)=bt(jk,i)
c                  r(ks,i,j)=rt(jk,i)
               end do
            end do
         end do
      end do

c extrapolation to axis
      
      kf=Lper*nk+1
      do i=1,2
         bax(i)=0.
c         rax(i)=0. 
      end do
      do i=1,2
         do lk=1,kf-1
            do j=1,nj
               bax(i)=bax(i)+b(lk,i,j)
c               rax(i)=rax(i)+r(lk,i,j)
            end do
         end do
      end do
      
      b00=sqrt((1.5*bax(1)-0.5*bax(2))/nj)
c      r00=(1.5*rax(1)-0.5*rax(2))/nj
      
      end subroutine extrap
