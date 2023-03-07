      subroutine makekicks(yn,Ep,lambda,B,val5,dt,ip,i,madecol,redstep)
c     $     ,tstep,ipall,ttot,step)
      use particle
      use machine
      use scatter
      use kicks
      use rand_var
      use para
      use interp
      use dstfun
      use pival
      use terp
      implicit none


      real ran,cosa,dt
      real yn(4),ynt(4),ynr(4),Ep,lambda,B,val5,val1(4)
      complex ef(2)
      real tau,tau1,denom,x
      real ds,dBds,R1,Z1,Bt
      real normp,ppar,mu,Omegac,vasq
      real vperp,vpar,arg,Besselp,Besselm
      real deltavperp(nph),deltavpar(nph),deltavperpmean
      integer n,ip,i,j,k,d,chk,madecol,inds,indp,indt
      real check,modenum(2),vtot,dummy,thetaper
      real argnow,redstep,nudot,nuDD,nuDDD
      real const,airconst,aairy
      real Pphi,f(4),tht,ph

      airconst=1./69.654248904899973
c      call random_seed
      normp=sqrt(2.*Ep*mo)
      ppar=normp*lambda
      mu=abs(normp*normp-ppar*ppar)/(2.*mo*B)
      vpar=ppar/mo
      vperp=sqrt(abs(normp*normp-ppar*ppar))/mo
      Omegac=Qe*B/mo
c$$$      call interpE(1,yn(1),yn(2),yn(3),val1,ef)
c$$$      Omegac=34.5*3./val1(1)*Qe/mo !B=B0*R0/R

      deltavpar=0.
      deltavperp=0.
      chk=0

      if(nkicklocal.eq.0)call InitKicks

c loop over toroidal mode numbers
      do d=1,nph    
         call interpE(d,yn(1),yn(2),yn(3),val1,ef)
         kparallel=val1(3)
         if(SELFO.gt.0)kparallel=ntor(d)/val1(1)
         
c vperp ICRH scattering after Kovanen et al. Nucl. Fus. 32,787 (1992)
         do n=1,harmonics
            nudot=omega-(n*Omegac+kparallel*vpar)
            check=nudot*ifarg(n,d,ip)
            if(check.le.0.0.and.abs(nkicklocal*ifarg(n,d,ip)).gt.0.)then

               madecol=1
               chk=1
               nkicks=nkicks+1

               Eplus=Efact(d)*ef(1)
               Eminus=Efact(d)*ef(2)

               kperpend=val1(4) ! from LEMan
               if(SELFO.gt.1)kperpend=kperpc ! constant kperp
               
               kparavg=kparavg+kparallel
               kperpavg=kperpavg+kperpend
               Ravg=Ravg+val1(1)
               
               nuDD=(nudot-ifarg(n,d,ip))/dt !nu..
               tau=sqrt(twopi/abs(nuDD))
               nuDDD=(nuDD-nuDDbef(n,d,ip))/dt!+dtbef(ip)) !nu...
               if(nuDDD**2..gt.airconst*abs(nuDD)**3.)then
c                  x=-nuDD*nuDD/(2.**(2./3.)*nuDDD**(4./3.))
! large number manipulation
                  x=abs(nuDDD)**(1./3)
                  x=x**4.
                  x=x*(2.**(2./3))
                  x=nuDD*nuDD/x
                  airyf=aairy(x)
                  tau=twopi*((2./abs(nuDDD))**(1./3.))*airyf
               endif
               ran=random(ran_index(i),ran_array(:,i))
c               cosa=sign(1.,-1.+2.*ran)   ! random gyrophase
               cosa=sqrt(3.)*2.*(ran-0.5)
               arg=kperpend*vperp/Omegac
      
               call Bessel(n-1,arg,Besselm)
               call Bessel(n+1,arg,Besselp)

               deltavperp(d)=sqrt(2.)*Qe/mo*tau*  !Murakami
     $              abs(Eplus*Besselm
     $              + Eminus*Besselp)

               deltavperp(d)=0.25*deltavperp(d)**2./vperp
     $              + cosa*sqrt(0.5)*deltavperp(d)


               deltavpar(d)=kparallel*vperp*deltavperp(d)/(n*Omegac)
               
               const=vperp**2.+(vpar-omega/kparallel)**2.
               dummy=(vperp+deltavperp(d))**2.
     $              +(vpar+deltavpar(d)-omega/kparallel)**2.
               if(abs(dummy/const-1.).gt.0.1)then
c$$$                  write(*,'(a,f5.2,i3,f5.2)')
c$$$     $                 'KICK NOT CONSERVATIVE: C/C0,d,kpar', dummy/const
c$$$     $                 ,d,kparallel
                  deltavperp(d)=0.
                  deltavpar(d)=0.
               endif
               

            endif
            nuDDbef(n,d,ip)=(nudot -ifarg(n,d,ip))/dt
            ifarg(n,d,ip)=nudot
         enddo
      enddo
      
      
c$$$      vperp=vperp+sum(deltavperp(:))
c$$$      vpar=vpar+sum(deltavpar(:))
      

      kickfactor=kickfactor+abs(sum(deltavperp(:))/vperp)
      kickfactor2=kickfactor2+abs(sum(deltavpar(:))/vpar)
      
      dtbef(ip)=dt
      
      vtot=(vperp+sum(deltavperp))**2.+(vpar+sum(deltavpar))**2. ! new v^2
      

      if(chk.eq.1)then
         if(stat.gt.0)then

            inds=0
            do j=0,nbinsrad-1
               if(yn(1).ge.intervs(j).and.yn(1).le.intervs(j+1))inds=j+1
            enddo
            tht=mod(mod(yn(2),twopi)+twopi,twopi)
            indp=floor(tht/twopi*nbinspol)+1
            indp=min(indp,nbinspol)
            ph=mod(mod(yn(3),twopiPer)+twopiPer,twopiPer)
            indt=floor(ph/twopiper*nbinstor)+1
            indt=min(indt,nbinstor)
         
            Enkicksip(inds,indp,indt,ip)=Enkicksip(inds,indp,indt,ip)
     $           + (0.5*mo*vtot - Ep)*weight(ip)

         else
            inds=1         
         endif
         do d=1,nph
            Enkicksradip(inds,ip,d)=Enkicksradip(inds,ip,d)
     $           + (0.5*mo*((vperp+deltavperp(d))**2.
     $           +(vpar+deltavpar(d))**2.) - Ep)*weight(ip)
         enddo
      endif

      vperp=vperp+sum(deltavperp)
      vpar=vpar+sum(deltavpar)

      Ep=0.5*mo*vtot
      lambda=vpar/sqrt(vtot)

      nkicklocal=1

      end

c--------------------------------------------------------------------
      subroutine InitKicks
      use particle
      use machine
      use terp
      use scatter
      use kicks
      use mpi_var
      use rand_var
      use para
      use interp
      use dstfun
      use pival
      implicit none

      integer j,k,l
      real rztmp,cosphi,sinphi

      
      allocate(ifarg(harmonics,nph,ntot/nprocs))
      allocate(nuDDbef(harmonics,nph,ntot/nprocs))
      allocate(dtbef(ntot/nprocs))
      allocate(Efact(nph))
      if(stat.eq.1)then
         allocate(Enkicks(nbinsrad,nbinspol,nbinstor))
         allocate(Enkicksip(nbinsrad,nbinspol,nbinstor,ntot/nprocs))
         allocate(Enkickstot(nbinsrad,nbinspol,nbinstor))
         allocate(RZ(3,nbinsrad,nbinspol,nbinstor))
         allocate(Enkicksrad(0:nbinsrad,nph))
         allocate(Enkicksradip(0:nbinsrad,ntot/nprocs,nph))
         allocate(Enkicksradtot(0:nbinsrad,nph))
      else
         allocate(Enkicksrad(1,nph))
         allocate(Enkicksradip(1,ntot/nprocs,nph))
         allocate(Enkicksradtot(1,nph))
      endif
      ifarg=0.
      kparavg=0.
      kperpavg=0.
      Ravg=0.
      kickfactor=0.
      kickfactor2=0.
      if(stat.eq.1)Enkicksip=0.
      Enkicksradip=0.
      Efact=1.
      omega=Qe*Bc/mo            ! ICRH heating @ B=Bc
      if(me2.eq.0)
     $     write(*,'(1x,a,1es9.2,a,f4.2,a)')'ICRH frequency: '
     $     ,omega,
     $     ' (',omega*mo/(Qe*b00),' x Omegac_o)'
      if(stat.eq.1)then
         RZbox(1)=1.4!0.8*b00*r00/Bc
         RZbox(2)=3.8!1.25*b00*r00/Bc
         RZbox(3)=-1.5!-1.3*a00
         RZbox(4)=1.5!1.3*a00
         do l=1,nbinstor
            cosphi=cos((l-0.5)*twopi/nbinstor/Lper)
            sinphi=sin((l-0.5)*twopi/nbinstor/Lper)
            do j=1,nbinspol
               do k=1,nbinsrad
                  rztmp=RZbox(1)   ! x
     $                 +(k-0.5)*(RZbox(2)-RZbox(1))/nbinsrad
                  RZ(1,k,j,l)=rztmp*cosphi
                  RZ(2,k,j,l)=rztmp*sinphi  ! y
                  RZ(3,k,j,l)=RZbox(3)     ! z
     $                 +(j-0.5)*(RZbox(4)-RZbox(3))/nbinspol
               enddo
            enddo
         enddo
      endif
      end
      
c--------------------------------------------------------------------
      subroutine Bessel(n,x,y)
      real  A(0:100)
      real  x,y
      integer m, n

      m=5
      call Bessel_coeff(m,n,A)

      y = A(0) + A(1) * x;
      if (m > 1)  y = y + A(2) * x * x;
      if (m > 2)  y = y + A(3) * x * x * x;
      if (m > 3)  y = y + A(4) * x**4
      do i = 4, m
         if (m > i)  y = y + A(i + 1) * x**(i+1)
      end do

      end

!*************************************************************
!* Bessel function series coefficient evaluation subroutine  *
!* m+1 is the number of coefficients desired, n is the order *
!* of the Bessel function. The coefficients are returned in  *
!* A(i).                                                     *
!*************************************************************
      subroutine Bessel_coeff(m,n,A)
      real A(0:m+n), B(0:m+n)  
      real a1,b1;
      integer i,m,n
      a1 = 1.d0; b1 = 1.d0;
      do i = 1, n
         B(i - 1) = 0.d0; b1 = b1 * i; a1 = a1 / 2.d0
      end do
      b1 = a1 / b1; a1 = 1.d0
      i = 0
      do while (i <= m)  
         A(i) = a1 * b1; A(i + 1) = 0
         a1 = -a1 / ((i + 2) * (n + n + i + 2))
         i = i + 2
      end do
      a1 = a1 / 2.d0;
      do i = 0, m
         B(i + n) = A(i)
      end do
      A = B
      return
      end
