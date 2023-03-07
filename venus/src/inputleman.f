      subroutine inputleman(anis,anis_p,sigmadens,sigmaanis)
      use particle
      use dimarray
      use dstfun
      use machine
      use interpos
      use para
      implicit none

      real,dimension(nbinsrad) ::  Tpa,N,anis
     $     ,N_int,p,pp,Tpa_int
      real,dimension(nsa) :: dens_p,Tpa_p,N_p,anis_p
     $     ,an_v,N_v,Tpa_v
      integer i
      real c(0:nbinsrad),dd,ca(0:nbinsrad)
      integer l,j,la,isnan
      real sv(nsa),tns(nsa),Noutp(nsa),Noutpp(nsa),ybc(6)
     $     ,Nin(nbinsrad+2),tensl(nbinsrad+2),sigmadens(nbinsrad)
     $     ,sigmaanis(nbinsrad),sigmaN(nbinsrad)
     $     ,Tin(nbinsrad+2),Toutp(nsa),Toutpp(nsa)
     $     ,Nip(nbinsrad),Nipp(nbinsrad)
      integer nbc(2),iflag,markmin
      real esin(nbinsrad+2),sgdensmin,dummy


c find density n(s) same way as anis(s)
      call findNT(Bc,anis,N,Tpa,sigmaanis,sigmadens,sigmaN)

c      Tpa=corr*Tpa
c      N=corr*N
      
      open(unit=31,file='moments.in')
      read(31,*)dummy ! Bcrit
      do i=1,nsa
         read(31,*)N_v(i),Tpa_v(i),an_v(i)
      enddo
      close(31)

      
      do i=1,nsa
         sv(i)=real(i-0.5)/nsa
      enddo
      tns=1.e-10
      nbc=0
      ybc=0.
      call cbsplgen0(sv,N_v,nsa,s(1:nbinsrad),N_int,Nip,Nipp,
     $     nbinsrad,0,0,tns,nbc,ybc,iflag)
      call cbsplgen0(sv,Tpa_v,nsa,s(1:nbinsrad),Tpa_int,Nip,Nipp,
     $     nbinsrad,0,0,tns,nbc,ybc,iflag)
      do i=1,nbinsrad
c         markmin=ntot
c         do j=1,nbinspolr(i)
c            if(nmarksall(j,i).lt.markmin)markmin=nmarksall(j,i)
c         enddo 
         N_int(i)=max(N_int(i),1.e-10)
         Tpa_int(i)=max(Tpa_int(i),1.e-10)
         write(11,'(3es14.6)')s(i),N(i),Tpa(i)        
c         if(isnan(N(i)).or.N(i).le.0.)N(i)=N_int(i)
c         if(isnan(Tpa(i)).or.Tpa(i).le.0.)Tpa(i)=Tpa_int(i)
c         if(markmin.lt.30)then!.and.i.gt.1)then
cc            N(i)=N_int(i)!N(i-1)
c            sigmaN(i)=0.
c            Tpa(i)=Tpa_int(i)!Tpa(i-1)
c         endif
      enddo        

      sgdensmin=huge(sgdensmin)!N_int(1)!sigmaN(1)
      do i=1,nbinsrad
         if(isnan(sigmaN(i)).or.sigmaN(i).le.0.)
     $        sigmaN(i)=N(i)!N_int(i) ! NaN
         if(sigmaN(i).lt.sgdensmin.and.sigmaN(i).gt.0.)
     $        sgdensmin=sigmaN(i)
      enddo

      esin(1)=0.
      Nin(1)=N(1)
      Tin(1)=Tpa(1)
      esin(nbinsrad+2)=1.
      Nin(nbinsrad+2)=0.!N_v(nsa)
      Tin(nbinsrad+2)=Tpa(nbinsrad)!Tpa_v(nsa)
      tensl=tens
c      tensl(1)=tens
c      tensl(nbinsrad+2)=tens
      do i=1,nbinsrad
         esin(i+1)=s(i)
         Nin(i+1)=N(i)
         Tin(i+1)=Tpa(i)
         tensl(i+1)=tens*sigmaN(i)/sgdensmin
c         if(N(i).eq.N_int(i))tensl(i+1)=10.*tensl(1)
         if(sigmaN(i).eq.0.)tensl(i+1)=10.*tensl(1)
      enddo      
c      tensl(1)=100*tensl(1)
c      tensl(nbinsrad+2)=tensl(1)!0.
      nbc(1)=0!2
      nbc(2)=0
      ybc=0.
      ybc(1)=0.!N(1)
c      tensl=0.
      call cbsplgen0(esin(2:nbinsrad+2),Nin(2:nbinsrad+2),nbinsrad+1
     $     ,sv,N_p,Noutp,Noutpp,
     $     nsa,0,1,tensl,nbc,ybc,iflag)
c      tensl=tens*100.
      call cbsplgen0(esin(2:nbinsrad+2),Tin(2:nbinsrad+2),nbinsrad+1
     $     ,sv,Tpa_p,Toutp,Toutpp,
     $     nsa,0,1,tensl,nbc,ybc,iflag)
      

      open(unit=31,file='moments.in')      
      write(31,'(es14.6)')Bc
      
      do i=1,nsa
         N_p(i)=max(N_p(i),1.e-10)
         Tpa_p(i)=max(Tpa_p(i),1.e-10)
         anis_p(i)=max(anis_p(i),1.e-10)
         write(31,'(3es14.6)')N_p(i),Tpa_p(i),anis_p(i)
      enddo

      end


c-----------------------------------------

  
      subroutine findNT(Bc,anis,N,Tpa,sigmaanis,sigmadens,sigmaN)
      use dstfun
      use particle
      use pival
      implicit none

      real Bc,N_d(nbinsrad),IsqA,sigloc
      real val(10),vals(8),pert(8),BoBc,onembc,
     $     intpol1,intpol2,C,H,dC
      integer i,j,k,check
      integer l,npol
      real sigmaanis(nbinsrad),sigmadens(nbinsrad),sigmaN(nbinsrad)
      real N(nbinsrad),anis(nbinsrad),Tpa(nbinsrad)


      do i=1,nbinsrad
         N_d(i)=0.
         Tpa(i)=0.
         sigmaN(i)=0.
         npol=0
!         IsqA=1.!/sqrt(anis(i))
         IsqA=sqrt(anis(i))
         do j=1,nbinspolr(i)
            intpol1=mod(twopi*real(j-0.5)/nbinspol-pi,twopi)
            do k=1,nbinstor
               if(nmarksall(k,j,i).gt.0)then
                  intpol2=twopiPer*real(k-0.5)/nbinstor
                  call interp3DRGC(s(i),intpol1,intpol2,val,vals,pert)
                  BoBc=sqrt(val(1))/Bc
                  call calcC(BoBc,anis(i),C,dC)
                  N_d(i)=N_d(i)+densang(k,j,i)/C*IsqA !N(s)=n_c*sqrt(A)
                  sigmaN(i)=sigmaN(i) + (sigmadens(i)/C)**2.
     $                 + (densang(k,j,i)/C*dC*sigmaanis(i))**2.
                  call calcH(BoBc,anis(i),H)
                  Tpa(i)=Tpa(i)+pparang(k,j,i)/densang(k,j,i)*C/H !*IsqA

                  npol=npol+1
               endif
            enddo
         enddo
         if(npol.gt.0)then
            Tpa(i)=Tpa(i)/Qel/npol  ! Tpa [eV]
            N(i)=N_d(i)/npol
            sigmaN(i)=sqrt(sigmaN(i))/npol
         else
            N(i)=0.
         endif
            
      enddo

      
      end


c--------------------------------------------------------
      subroutine calcC(BoBc,A,C,dC)
      implicit none

      real BoBc,A,C,dC

      if(BoBc.ge.1.)then
         C=BoBc/(1-A*(1.-BoBc))
         dC=(1.-BoBc)/(1-A*(1.-BoBc))
      else
         C=BoBc*(1.+A*(1.-BoBc)-2.*(A*(1.-BoBc))**1.5)
         C=C/((1.-A*(1.-BoBc))*(1.+A*(1.-BoBc)))
         dC=(2.-BoBc-3.*sqrt(A)*(1.-BoBc)**1.5)
     $        /(1.+A*(1.-BoBc)-2.*(A*(1.-BoBc))**1.5)
         dC=dC + 2.*A*(1.-BoBc)**2./(1.-(A*(1-BoBc))**2.)         
      endif

      end
c--------------------------------------------------------
      subroutine calcH(BoBc,A,H)
      implicit none

      real BoBc,A,H

      if(BoBc.ge.1.)then
         H=BoBc/(1.-A*(1.-BoBc))
      else
         H=BoBc*(1.+A*(1.-BoBc)-2.*(A*(1.-BoBc))**2.5)
         H=H/((1.-A*(1.-BoBc))*(1.+A*(1.-BoBc)))
      endif

      end

c--------------------------------------------------------
