      subroutine array

      use particle
      use terpall
      use arrays
      use mpi_var
      use sinitial
      implicit none

      integer ixi,jth,ks,l,ii
      REAL r00,twopi,dxi,dth,before,tt,after,xi,theta,fact
      REAL jac,cosarg(lmnb),sinarg(lmnb)


c tabarray: equilibrium functions of (s,theta,xi)
c
c****************TABARRAY****************************
c (1) ---> B^2              (5)---> sigma           *
c (2) ---> DB/Ds            (6)---> tau             *
c (3) ---> DB/Dtheta        (7)---> Dsigma/Ds       *
c (4) ---> DB/Dxi           (8)---> sigmaBs         *
c                           (9)---> DsigmaBs/Dtheta *
c                          (10)---> DsigmaBs/Dxi    *
c****************************************************
c
c tabarrays: equilitbrium functions of s only
c
c****************TABARRAYS***************************
c (1) ---> muo I            (5) ---> phi' (toro)    *
c (2) ---> muo J            (6) ---> psi' (polo)    *
c (3) ---> muo I'           (7) ---> phi            *
c (4) ---> muo J'           (8) ---> psi            *
c****************************************************
c
c pertarray: perturbations
c
c****************PERTARRAY***************************
c (1) ---> gamma            (5) ---> Dchi/Ds        *
c (2) ---> Dgamma/Dtheta    (6) ---> Dchi/Dtheta    *
c (3) ---> Dgamma/Dxi       (7) ---> Dchi/Dxi       *
c (4) ---> Dgamma/Ds        (8) ---> Dgamma/Dtau    *
c****************************************************
    
      twopi=4.*acos(0.)

      dxi=twopi/Lper/nxia
      dth=twopi/ntha
      ds=1.0/nsa
      ds2=ds/2.
      

      do l=1,8
         do ks=0,nsa
            tabarrays(ks,l)=0.
            do ixi=1,nxia+1
               do jth=1,ntha+1
                  tabarray(jth,ixi,ks,l)=0.
c                  pertarray(jth,ixi,ks,l)=0.
               enddo
            enddo
         enddo
      enddo
      do l=9,10
        do ks=0,nsa
            do ixi=1,nxia+1
               do jth=1,ntha+1
                  tabarray(jth,ixi,ks,l)=0.
               enddo      
            enddo
         enddo
      enddo
      sinind=floor(sinit*nsa)
      sinind=max(sinind,0)
      sinind=min(sinind,nsa)

c      print *,'
c      print *,'********COMPUTE TABRGC********'       
      do ks=1,nsa         
         do ixi=1,nxia+1
            xi=(ixi-1)*dxi ! xi in [0,twopi/Lper]
            do jth=1,ntha+1
               theta=(jth-1)*dth   ! theta in [0,twopi]

               do l=1,lmnb
                  cosarg(l)=cos(mb(l)*theta - nm(l)*xi) !calculate cosine for all transformations
                  tabarray(jth,ixi,ks,1) = tabarray(jth,ixi,ks,1) + 
     >                 bsmn(l,sinind+ks)*cosarg(l) ! B^2
               enddo

               fact=2.*sqrt(tabarray(jth,ixi,ks,1))

               do l=1,lmnb
                  sinarg(l)=sin(mb(l)*theta - nm(l)*xi) !calculate sine for all transformations
                  tabarray(jth,ixi,ks,3) = tabarray(jth,ixi,ks,3) -
     >                 mb(l)*bsmn(l,sinind+ks)*sinarg(l) ! dB^2/dtheta
                  tabarray(jth,ixi,ks,4) = tabarray(jth,ixi,ks,4) +
     >                 nm(l)*bsmn(l,sinind+ks)*sinarg(l) ! dB^2/dxi
               enddo
               tabarray(jth,ixi,ks,3) = tabarray(jth,ixi,ks,3)/fact
               tabarray(jth,ixi,ks,4) = tabarray(jth,ixi,ks,4)/fact

               if (ks.ne.1)
     >              tabarray(jth,ixi,ks-1,2) = 
     >              (sqrt(tabarray(jth,ixi,ks,1))-
     >              sqrt(tabarray(jth,ixi,ks-1,1)))/ds ! dB/ds

               do l=1,lmnb
                  tabarray(jth,ixi,ks,5)=tabarray(jth,ixi,ks,5) +
     >                 sigmn(l,sinind+ks)*cosarg(l) ! sigma 
                  tabarray(jth,ixi,ks,6)=tabarray(jth,ixi,ks,6) +
     >                 taumn(l,sinind+ks)*cosarg(l) ! tau
                  tabarray(jth,ixi,ks,7)=tabarray(jth,ixi,ks,7) +
     >                 dpdsmn(l,sinind+ks)*cosarg(l) ! Dpperp/Ds-Dppar/Ds
                  tabarray(jth,ixi,ks,8)=tabarray(jth,ixi,ks,8) +
     >                 sigBsmn(l,sinind+ks)*cosarg(l) ! sigmaBs
               enddo

               do l=1,lmnb
                  tabarray(jth,ixi,ks,9) = tabarray(jth,ixi,ks,9) -
     >                 mb(l)*sigBsmn(l,sinind+ks)*sinarg(l) ! DsigBs/Dtheta
                  tabarray(jth,ixi,ks,10) = tabarray(jth,ixi,ks,10) +
     >                 nm(l)*sigBsmn(l,sinind+ks)*sinarg(l) ! DsigBs/Dxi
               enddo

            enddo
         enddo
         tabarrays(ks,1) =ci(sinind+ks) ! muoI
         tabarrays(ks,2) =cj(sinind+ks) ! muoJ

c isaev
c              print *,'cip, cjp',cip(ks),cjp(ks)

         tabarrays(ks,3) =cip(sinind+ks) ! muoI'
         tabarrays(ks,4) =cjp(sinind+ks) ! muoJ'
         tabarrays(ks,5) =psip(sinind+ks) ! phi'
         tabarrays(ks,6)=chip(sinind+ks) ! psi'
         tabarrays(ks,7)=psi(sinind+ks) ! phi
         tabarrays(ks,8)=chi(sinind+ks) ! psi


         jac=tabarrays(ks,6)*tabarrays(ks,2) -
     >        tabarrays(ks,5)*tabarrays(ks,1)

         do ixi=1,nxia+1
            do jth=1,ntha+1
               tabarray(jth,ixi,ks,7)=tabarray(jth,ixi,ks,5)*
     >              tabarray(jth,ixi,ks,7)/jac !Dsigma/Ds @B=cst
            enddo
         enddo 
c        print *,'  Plan xi',ixi,'    TIME=',after-before
c         print*,'sinind,sinind+ks',sinind,sinind+ks
      enddo

      do ixi=1,nxia+1
         do jth=1,ntha+1
            tabarray(jth,ixi,nsa,2)=tabarray(jth,ixi,nsa-1,2) ! dB/ds at s=stops
         enddo
      enddo

      if(sinind.eq.0)call getaxis

      call deallocateEq

      end

c--------------------------------------------------------------------------------------
c
c     interpolation/extrapolation of quantities on magnetic axis
c
c-------------------------------------------------------------------------------------

      subroutine getaxis

      use arrays
      implicit none

      integer n,ixi,jth
      
    ! interpolation of quantities depending on theta
c      print*,'GETAXIS'
      do n=1,8
         tabarrays(0,n)=1.5*tabarrays(1,n)-.5*tabarrays(2,n)
         if(n.eq.7.or.n.eq.8)tabarrays(0,n)=0. !fluxes are zero on axis
         do ixi=1,nxia+1
            do jth=1,ntha+1
              tabarray(jth,ixi,0,n)=1.5*tabarray(jth,ixi,1,n)-
     $              .5*tabarray(jth,ixi,2,n)
!              pertarray(jth,ixi,0,n)=1.5*pertarray(jth,ixi,1,n)-
!     $             .5*pertarray(jth,ixi,2,n)
c               tabarray(jth,ixi,0,n)=.5*
c     $              (tabarray(jth,ixi,1,n)+tabarray(jth+ntha/2,ixi,1,n))
c               pertarray(jth,ixi,0,n)=.5*(pertarray(jth,ixi,1,n)+
c     $              pertarray(jth+ntha/2,ixi,1,n))
c               if(n.eq.2)then ! dBds zero on axis (anti-symmetric)
c                  tabarray(jth,ixi,0,n)=0.
c               endif
            enddo
c            do jth=ntha/2+2,ntha+1
c               tabarray(jth,ixi,0,n)=.5*
c     $              (tabarray(jth,ixi,1,n)+tabarray(jth-ntha/2,ixi,1,n))
c               pertarray(jth,ixi,0,n)=.5*(pertarray(jth,ixi,1,n)+
c     $              pertarray(jth-ntha/2,ixi,1,n))
c              tabarray(jth,ixi,0,n)=1.5*tabarray(jth,ixi,1,n)+
c    $              .5*tabarray(jth,ixi,2,n)
c              pertarray(jth,ixi,0,n)=1.5*pertarray(jth,ixi,1,n)+
c    $              .5*pertarray(jth,ixi,2,n)
c               if(n.eq.2)then ! dBds zero on axis (anti-symmetric)
c                  tabarray(jth,ixi,0,n)=0.
c               endif

c            enddo
               
         enddo
      enddo
      do n=9,10
         do ixi=1,nxia+1
            do jth=1,ntha+1
               tabarray(jth,ixi,0,n)=1.5*tabarray(jth,ixi,1,n)-
     $              .5*tabarray(jth,ixi,2,n)               
c               tabarray(jth,ixi,0,n)=.5*
c     $              (tabarray(jth,ixi,1,n)+tabarray(jth+ntha/2,ixi,1,n))
            enddo
c            do jth=ntha/2+2,ntha+1
c              tabarray(jth,ixi,0,n)=1.5*tabarray(jth,ixi,1,n)+
c    $              .5*tabarray(jth,ixi,2,n)               
c               tabarray(jth,ixi,0,n)=.5*
c     $              (tabarray(jth,ixi,1,n)+tabarray(jth-ntha/2,ixi,1,n))
c            enddo            
         enddo
      enddo
    ! extrapolation of quantities depending on s only
c      do n=1,8
c         tabarrays(0,n)=1.5*tabarrays(1,n)-.5*tabarrays(2,n)
c      enddo

      end
               

c---------------------------
      subroutine deallocateEq
      use terpall
      use arrays
      implicit none

      deallocate(mb,nm)
      deallocate(si,psi,chi,psip,chip,
     >     ci,cj,cip,cjp,bsmn,sigmn,taumn,dpdsmn,sigBsmn)

      end
