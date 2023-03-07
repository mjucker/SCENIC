      subroutine changepitch(ss,theta,xi,lambda,nudpast,Ep,past
     $     ,tsi,tse,dnuedee,dnuedei,i,ip)
      use particle
      use plasma
      use scatter
      use mpi_var
      use rand_var
      use para
      use dstfun
      use pival
      use terp
      implicit none


      real nudpast,lambda,dlambda
      real Ep,past,tsi,tse
      real nuepast,dnuedee,dnuedei,nuipast
      real plusminus,flam
      integer ip,i,j,choice,ind,indp,indt
      real dummy,ss,theta,dE,dEe,dEi
      real delT,phi,tht,ph,xi



      choice=1

c debut pitch angle scattering. Boozer Kuo-Petravic. Phys. Fluids 24,851(1981) 
c make random first random call

      flam=1.0-lambda*lambda
      flam=max(flam,0.)
      
       
      if(choice.eq.1)then
c pitch angle scattering Boozer Kuo-Petravic
         plusminus=random(ran_index(i),ran_array(:,i))
c      call random_gauss(plusminus,1,ran_index(i),ran_array(:,i))
         plusminus=sign(1.,plusminus-0.5)
c         plusminus=sqrt(3.)*2.*(plusminus-0.5)
         lambda=lambda*(1.0-nudpast) + 
     $        plusminus*sqrt(flam*nudpast) ! cos(pitch). Lorentz scattering

      elseif(choice.eq.2)then
c pitch angle scattering as proposed by Vernay
         call random_gauss(plusminus,1.,ran_index(i),ran_array(:,i))
         delT=2.*plusminus*sqrt(nudpast)
         phi=random(ran_index(i),ran_array(:,i))
         phi=sin(phi*twopi)
         lambda=sin(delT)*phi*sqrt(flam) + lambda*cos(delT)
      endif


c     scattering Energy. Only NRel.  Boozer and Kuo-Petravic. Phys. Fluids 24, 851 (1981)
      ind=0
      do j=0,nbinsrad-1
         if(ss.ge.intervs(j).and.ss.le.intervs(j+1))ind=j+1
      enddo
      tht=mod(mod(theta,twopi)+twopi,twopi)
      indp=floor(tht/twopi*nbinspol)+1
      indp=min(indp,nbinspol)
      ph=mod(mod(xi,twopiPer)+twopiPer,twopiPer)
      indt=floor(ph/twopiPer*nbinstor)+1
      indt=min(indt,nbinstor)

c complete scattering
      nuepast=past/tsi
      nuipast=nuepast

c      plusminus=random(ran_index(i),ran_array(:,i))
c      plusminus=sign(1.,plusminus-0.5)

      dEi= - 2.*nuepast*(Ep - dnuedei)
c     $     +plusminus*2.*sqrt(Ep*Ti*Qel*nuepast)

      nuepast=past/tse

      plusminus=random(ran_index(i),ran_array(:,i))
      plusminus=sign(1.,plusminus-0.5)

      dEe= - 2.*nuepast*(Ep - dnuedee)
c     $     +plusminus*2.*sqrt(Ep*Te*Qel*nuepast)


c total scattering
      dE=dEi+dEe
     $     +plusminus*2.*sqrt(Ep*Qel*(Ti*nuipast+Te*nuepast))
      Ep=Ep+dE

c statistics
      if(stat.eq.1)then
         Encollip(ind,1,ip)=Encollip(ind,1,ip)
     $        +dE*tse/(tsi+tse)*weight(ip)
         Encollip(ind,2,ip)=Encollip(ind,2,ip)
     $        +dE*tsi/(tsi+tse)*weight(ip)
         Encollall(ind,indp,indt)=Encollall(ind,indp,indt)
     $        +dE*weight(ip)
      endif

      end


