      subroutine updatenu(Ep,yn,nud,tsi,tse,dnuedee,dnuedei)
      use particle
      use plasma
      use terp
      use scatter
      use dimarray
      use dstfun
      use pival
      implicit none

      REAL Ep,nud,ss,yn(4),erfp,profile,
     >     logCi,logCe,phi,psi,xboo,
     >     QeQp2,QeQel2,nude,nudi,nub,
     >     pden,logC,nuee,nuei,F,dnuedee,dnuedei,psip,
     >     vpart,twopiesp2m2,ptherm
     $     ,tsi,tse,Ad,tdi,tde
     $     ,ne,np,nh,thta,xi
      real nue,dnue,dnuede
      real dE,nudde,dtsi,dtse
      real dummy
      integer i,j,k

      
      QeQp2=Qe*Qe*qplasma*qplasma !q(part)^2 * q(background)^2
      QeQel2=Qe*Qe*Qel*Qel        !q(part)^2 * q(electrons)^2
      
      twopiesp2m2=twopi*epso*epso*mo*mo
    
      ss=yn(1)
      thta=mod(mod(yn(2),twopi)+twopi,twopi)
      xi=mod(mod(yn(3),twopiPer)+twopiPer,twopiPer)
c      profile=1.
c      profile=(1.-0.8*ss)**2.
c      profile=1. -0.795*ss +0.115*ss**2 -0.305*ss**3
c      profile=1. -0.06*ss -0.32*ss**2 -3.11*ss**3 +7.8*ss**4 -5.07*ss**5
      profile=n0 +n1*ss +n2*ss**2. +n3*ss**3. +n4*ss**4. +n5*ss**5. 
     
      pden=max(0.,density*profile)
      if(pden.le.0.)print*,ss,profile
      ne=pden
      i=floor(ss*nsa)+1
      j=floor(thta/twopi*nbinspol)+1
      k=floor(xi/twopiPer*nbinstor)+1
      nh=densi(k,j,i)
      np=(ne-Zpart*nh)*Qel/qplasma
      if(np.lt.0.)then
         print*,'np<0: ne,nh',ne,nh
         stop
      endif
      call interpPP(ss,ptherm)
      Te=ptherm/(np*TiTe+ne)/Qel ! [eV]
      Ti=TiTe*Te

c      logCi=30.-log(Zpart*Zpart*Zpart*sqrt(pden)*Tp**(-1.5))
c      logCe=31.3 - log(sqrt(pden)/Tp)
cj Coulomb logarithm as in Wesson p. 726
      logCi=30.3 + 1.5*log(Ti) - log(sqrt(np)*Zpart*qplasma/Qel)
c      logCe=30.3 + 1.5*log(Te) - log(sqrt(ne)*Zpart) 
cj SELFO for Z=1, Wesson p. 727
c      logCi=17.3 - 0.5*log(ne/(10.**20.)) + 1.5*log(Ti/1000.)
c      logCe=15.2 - 0.5*log(ne/(10.**20.)) + log(Te/1000.)
cj NRL p. 34: ne[cm-3]=ne[m-3]/1e6, 31.3-24=ln(1e6), same as SELFO
      logCe=31.3 - 0.5*log(ne) + log(Te)
c      logC=23-log(sqrt(pden*(1.e-06))/T0)

cRel      normp=sqrt(gamma*gamma-1.0)*mo*vc ! p=sqrt(gamma^2-1)*mo*c


c pitch angle scattering: nu estimated by Strumberger. Nucl. Fusion 40, 1697 (2000) and as in Wesson p.65
      vpart=sqrt(2.*Ep*mo)/gamma/mo

c  on ions (plasma ion mass=mplasma) [before= ion plasma mass=2.*mo]

      xboo=vpart/sqrt(2.*Ti*Qel/mplasma)
      phi=erf(xboo)
      erfp=exp(-xboo*xboo)*2./sqrt(pi)
      psi=(phi-xboo*erfp)/(2.*xboo*xboo)
c ion on ion slowing-down time
      Ad=np*QeQp2*logCi/twopiesp2m2
c      tsi=2.*Tp*Qel/mplasma*vpart/((1+mo/mplasma)*Ad*psi)
c      tsi=vpart**3./((1+mo/mplasma)*Ad*psi*xboo*xboo)
      tsi=vpart**3./((mo/mplasma)*Ad*psi*xboo*xboo)
c      tsi=0.5*vpart**3./(Ad*psi*xboo*xboo)
c ion on ion deflection time
      tdi=vpart**3./(0.5*Ad*(phi-psi))
      nudi=1./tdi

      nuei=1./tsi
      
c On electrons
      xboo=vpart/sqrt(2.*Te*Qel/me)
      phi=erf(xboo)
      erfp=exp(-xboo*xboo)*2./sqrt(pi)
      psi=(phi-xboo*erfp)/(2.*xboo*xboo)
c ion on electron slowing-down time
      Ad=ne*QeQel2*logCe/twopiesp2m2
c      tse=2.*Tp*Qel/me*vpart/((1+mo/me)*Ad*psi)
c      tse=vpart**3./((1+mo/me)*Ad*psi*xboo*xboo)
      tse=vpart**3./((mo/me)*Ad*psi*xboo*xboo)
c      tse=0.5*vpart**3./(Ad*psi*xboo*xboo)
c ion on electron deflection time
      tde=vpart**3./(0.5*Ad*(phi-psi))
      nude=1./tde

      nuee=1./tse

      nud=nudi+nude

      
      nue=nuei+nuee

c      nue=1./tse
c      print*,'tse,tsi',tse,tsi
c      tse=6.27e8*Apart*Tp**1.5/(Zpart**2*pden*1e-6*logCe) !Spitzer time
c      print*,'ts',tse
c      nue=1./tse      !Spitzer time

c compute dnue/dE
      dE=0.001*Ep
      vpart=sqrt(2.*(Ep+dE)*mo)/gamma/mo

      xboo=vpart/sqrt(2.*Ti*Qel/mplasma)
      phi=erf(xboo)
      erfp=exp(-xboo*xboo)*2./sqrt(pi)
      psi=(phi-xboo*erfp)/(2.*xboo*xboo)
      Ad=np*QeQp2*logCi/twopiesp2m2
c      tsi=2.*Tp*Qel/mplasma*vpart/((1+mo/mplasma)*Ad*psi)
c      dtsi=vpart**3./((1+mo/mplasma)*Ad*psi*xboo*xboo)
      dtsi=vpart**3./((mo/mplasma)*Ad*psi*xboo*xboo)


      xboo=vpart/sqrt(2.*Te*Qel/me)
      phi=erf(xboo)
      erfp=exp(-xboo*xboo)*2./sqrt(pi)
      psi=(phi-xboo*erfp)/(2.*xboo*xboo)
      Ad=ne*QeQel2*logCe/twopiesp2m2
c      tse=2.*Tp*Qel/me*vpart/((1+mo/me)*Ad*psi)
c      dtse=vpart**3./((1+mo/me)*Ad*psi*xboo*xboo)
      dtse=vpart**3./((mo/me)*Ad*psi*xboo*xboo)
      
c      tse=6.27e8*Apart*Tp**1.5/(Zpart**2*pden*1e-6*logCe) !Spitzer time


      dnue=1./dtsi
      dnuedei=(dnue-nuei)/dE

      dnue=1./dtse
      dnuedee=(dnue-nuee)/dE
           

      dnuedei=(1.5+Ep*dnuedei*tsi)*Ti*Qel
      dnuedee=(1.5+Ep*dnuedee*tse)*Te*Qel


 
      end
