      subroutine inputvmec(sigmapparang,sigmapperpang,sigmaanis
     $     ,TpoTb,TT_o)
      use dimarray
      use machine
      use plasma
      use scatter
      use dstfun
      use interpos
      use para
      use pival
      implicit none

      real ph(nbinsrad),polfac(nbinsrad)
      real pp_v(nsa),ph_v(nsa),A_v(nsa),iota_v(nsa)
     $     ,pp_b(nsa)
     $     ,TT_o(nsa),ph_o(nsa),pp_o(nsa),nth_o(nsa),Tth_o(nsa)
      real val(10),vals(8),pert(8),BoBc,TpoTb(nbinsrad),ptherm,onembc
     $     ,intpol1,intpol2
      real TpoTb_d(nbinsrad),polfacloc,H
      integer i,j,k,l,nbins,isnan
      real c(0:nbinsrad),dd
      real sv(nsa),tns(nsa),php(nsa),phpp(nsa),ybc(6)
     $     ,phin(nbinsrad+2),tensa(nbinsrad+2),tensph(nbinsrad+2)
     $     ,TTin(nbinsrad+2),TToutp(nbinsrad),TToutpp(nbinsrad)
     $     ,TTop(nsa),TTopp(nsa)
     $     ,ppin(nbinsrad+2),nthin(nbinsrad+2)
      real,allocatable :: io_int(:),currin(:),curr_o(:)
      integer nbc(2),iflag,markmin,npol
      real esin(nbinsrad+2),sigmina,sigminph
     $     ,sigmapparang(nbinstor,nbinspol,nbinsrad)
     $     ,sigmapperpang(nbinstor,nbinspol,nbinsrad)
     $     ,sigmaanis(nbinsrad),sigmaph(nbinsrad)
      real A_int(nbinsrad),P_int(nbinsrad),Pp_int(nbinsrad)
      real mu0,Qel
      logical f_exist

      if(vcurr.eq.1)allocate(io_int(nbinsrad)
     $     ,currin(nbinsrad+2),curr_o(nsa))
      
      mu0=2.*twopi*1.e-7
      Qel=1.6022e-19

      open(unit=21,file='profiles.in')

      do i=1,nsa
         read(21,*)pp_v(i),ph_v(i),A_v(i),iota_v(i)
      enddo
      close(21)

      INQUIRE(FILE='thermal.in', EXIST=f_exist)  ! thermal.in existing?
      if(f_exist)then
         open(unit=21,file='thermal.in')
         do i=1,nsa
            read(21,*)pp_b(i),nth_o(i),TT_o(i)
         enddo
         close(21)
      else
         pp_b=pp_v
      endif
      
         
      
c prepare interpolations
      do i=1,nsa
         sv(i)=real(i-0.5)/nsa
      enddo
      tns=1.e-10
      nbc=0
      ybc=0.
      call cbsplgen0(sv,A_v,nsa,s(1:nbinsrad),A_int,TToutp,TToutpp,
     $     nbinsrad,0,1,tns,nbc,ybc,iflag)
      call cbsplgen0(sv,ph_v,nsa,s(1:nbinsrad),P_int,TToutp,TToutpp,
     $     nbinsrad,0,1,tns,nbc,ybc,iflag)
      call cbsplgen0(sv,pp_b,nsa,s(1:nbinsrad),Pp_int,TToutp,TToutpp,
     $     nbinsrad,0,1,tns,nbc,ybc,iflag)
      if(vcurr.eq.1)
     $     call cbsplgen0(sv,iota_v,nsa,s(1:nbinsrad),io_int
     $     ,TToutp,TToutpp,
     $     nbinsrad,0,1,tns,nbc,ybc,iflag)

c find Anisotropy
      call FindA(sigmapparang,sigmapperpang,Bc,TpoTb,sigmaanis)

      
c use found anisotropy for ph
      polfac=0.
      sigmaph=0.
      do i=1,nbinsrad
         P_int(i)=max(P_int(i),1.e-10)
         A_int(i)=max(A_int(i),1.e-10)
         Pp_int(i)=max(Pp_int(i),1.e-10)
         if(isnan(TpoTb(i)).or.TpoTb(i).le.0.)TpoTb(i)=A_int(i) ! NaN
         npol=0
         do j=1,nbinspolr(i)
            intpol1=mod(twopi*real(j-0.5)/nbinspolr(i)-pi,twopi)
            do k=1,nbinstor
               if(nmarksall(k,j,i).gt.0)then
                  intpol2=twopiPer*real(k-0.5)/nbinstor
                  call interp3DRGC(s(i),intpol1,intpol2,val,vals,pert)
                  BoBc=sqrt(val(1))/Bc
                  call calcH(BoBc,TpoTb(i),H)
                  polfacloc=pparang(k,j,i)/H
                  polfac(i)=polfac(i)+polfacloc
               
                  sigmaph(i)=sigmaph(i)+
     $                 (1./H)**2.*sigmapparang(k,j,i)**2.
                  onembc=1.-BoBc
                  sigmaph(i)=sigmaph(i)
     $                 +(polfacloc*onembc/(1.-TpoTb(i)*onembc))**2.
     $                 *sigmaanis(i)**2.
                  npol=npol+1
               endif
            enddo
         enddo

         if(npol.gt.0)then
            polfac(i)=polfac(i)/npol ! average over poloidal bins
            sigmaph(i)=sigmaph(i)/npol

            call interpPP(s(i),ptherm)
            ph(i)=polfac(i)/ptherm
            if(isnan(ph(i)).or.ph(i).le.0.)ph(i)=P_int(i)
            sigmaph(i)=sqrt(sigmaph(i))/ptherm
         else
            ph(i)=0.!P_int(i)
            sigmaph(i)=0.!P_int(i)!ph(i)
         endif

      enddo
      
c      ph=corr*ph

      if(vcurr.eq.1)call findcurr
      

      do i=1,nbinsrad
         markmin=ntot
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               if(nmarksall(k,j,i).lt.markmin)markmin=nmarksall(k,j,i)
            enddo
         enddo
         write(14,'(6es14.6)')s(i),ph(i),pth(i)
     $        ,TpoTb(i),densth(i),mu0*pth(i)+Pp_int(i)
c         if(markmin.lt.30.and.i.gt.1)then
c            TpoTb(i)=A_int(i)!TpoTb(i-1)
c$$$            ph(i)=ph(i-1)!P_int(i)
c$$$c            write(*,'(i3,a12,i3)')markmin,' ptcls in i=',i
c         endif
      enddo         

      esin(1)=0.
      TTin(1)=TpoTb(1)          !A_int(1)
      phin(1)=ph(1)
      if(densth(1).gt.0.)then
         ppin(1)=pth(1)/densth(1)
      else
         ppin(1)=0.
      endif
      nthin(1)=densth(1)
      if(vcurr.eq.1)currin(1)=0.!currtot(1)
      esin(nbinsrad+2)=1.
      TTin(nbinsrad+2)=TpoTb(nbinsrad)!A_v(nsa)  !TpoTb(nbinsrad)
      phin(nbinsrad+2)=ph(nbinsrad)!ph_v(nsa)
      if(densth(nbinsrad).gt.0.)then
         ppin(nbinsrad+2)=pth(nbinsrad)/densth(nbinsrad)
      else
         ppin(nbinsrad+2)=0.
      endif
      nthin(nbinsrad+2)=densth(nbinsrad)
      if(vcurr.eq.1)currin(nbinsrad+2)=currtot(nbinsrad)
!      if(isnan(sigmaanis(1)).or.sigmaanis(1).eq.0.)
!     $     sigmaanis(1)=TpoTb(1)
!      if(isnan(sigmaph(1)).or.sigmaph(1).eq.0.)
!     $     sigmaph(1)=ph(1)aph
      do i=1,nbinsrad
         esin(i+1)=s(i)
         TTin(i+1)=TpoTb(i)
         phin(i+1)=ph(i)
         if(densth(i).gt.0.)then
            ppin(i+1)=pth(i)/densth(i)
         else
            ppin(i+1)=0.
         endif
         nthin(i+1)=densth(i)
         if(vcurr.eq.1)currin(i+1)=currtot(i)
         if(isnan(sigmaanis(i)).or.sigmaanis(i).le.0.)
     $        sigmaanis(i)=TpoTb(i)!10.*sigmina ! rid of NaN
         if(isnan(sigmaph(i)).or.sigmaph(i).le.0.)
     $        sigmaph(i)=ph(i)!10.*sigminph ! get rid of NaN
      enddo
c      sigmina=A_int(1)!sigmaanis(1)
c      sigminph=P_int(1)!sigmaph(1)
      sigmina=huge(sigmina)
      sigminph=sigmina
      do i=1,nbinsrad
         if(sigmaanis(i).lt.sigmina.and.sigmaanis(i).gt.0.)
     $        sigmina=sigmaanis(i)
         if(sigmaph(i).lt.sigminph.and.sigmaph(i).gt.0.)
     $        sigminph=sigmaph(i)
      enddo     
      
      nbc=0
      ybc=0.
      tensa=tens
      tensph=tens
c      tensa(1)=tens
c      tensph(1)=tens
      do i=1,nbinsrad
         tensa(i+1)=tens*sigmaanis(i)/sigmina
         tensph(i+1)=tens*sigmaph(i)/sigminph
c$$$         if(ph(i).eq.P_int(i))then
c$$$            tensa(i+1)=10.*tensa(1)
c$$$            tensph(i+1)=10.*tensph(1)
c$$$         endif
      enddo
      tensall=tensa
c      tensa(1)=100*tensa(1)
c      tensph(1)=100*tensph(1)
c      tensa(nbinsrad+2)=0.!tensa(1)
c      tensph(nbinsrad+2)=0.!tensph(1)

      nbc(1)=0
c      nbc(1)=1
      nbc(2)=0
      ybc=0.
c      ybc(2)=ph(nbinsrad+2)
c      tens=0.e-4
      call cbsplgen0(esin,phin,nbinsrad+2,sv,ph_o,php,phpp,
     $     nsa,0,1,tensph,nbc,ybc,iflag)
c      nbc=0
c      ybc(2)=TTin(nbinsrad+2)
c      nbc(2)=1
      call cbsplgen0(esin(2:nbinsrad+2),TTin(2:nbinsrad+2),nbinsrad+1
     $     ,sv,TT_o,TTop,TTopp,
     $     nsa,0,1,tensa(2:nbinsrad+2),nbc,ybc,iflag)
c      nbc(2)=0
      call cbsplgen0(esin(2:nbinsrad+2),ppin(2:nbinsrad+2),nbinsrad+1
     $     ,sv,pp_o,TTop,TTopp,
     $     nsa,0,1,tensa(2:nbinsrad+2),nbc,ybc,iflag)

      call cbsplgen0(esin(2:nbinsrad+2),nthin(2:nbinsrad+2),nbinsrad+1
     $     ,sv,nth_o,TTop,TTopp,
     $     nsa,0,1,tensa(2:nbinsrad+2),nbc,ybc,iflag)
      
      if(vcurr.eq.1)then
         nbc(1)=0
         tensa(1)=tens*0.01!0.
         call cbsplgen0(esin,currin
     $        ,nbinsrad+1,sv
     $        ,curr_o,TTop,TTopp,
     $        nsa,0,1,tensa,nbc,ybc,iflag)
      endif


c       if(s(1).gt.1./2./nsa.or.s(nbinsrad).lt.1.-1./2./nsa)print*,
c     $     'CAUTION: radial profiles being extrapolated!'

c      call extraProf(s,ph,TpoTb,ph_o,TT_o)
      
      open(unit=21,file='profiles.in')
      do i=1,nsa
         pp_v(i)=pp_b(i)+mu0*pp_o(i)
         ph_o(i)=max(ph_o(i),1.e-10)
         TT_o(i)=max(TT_o(i),1.e-10)
         if(vcurr.eq.0)then
            write(21,'(4es14.6)')pp_v(i),ph_o(i),TT_o(i),iota_v(i)
         else
            curr_o(i)=mu0*curr_o(i)
            write(21,'(4es14.6)')pp_v(i),ph_o(i),TT_o(i),curr_o(i) ! sign(jacobian of VMEC)
         endif            
      enddo
      close(21)
      
      open(unit=21,file='thermal.in')
      do i=1,nsa
         nth_o(i)=max(nth_o(i),1.e-10)
         ph_o(i)=max(pp_o(i)/Qel,1.e-10) ! Tth
         write(21,'(3es14.6)')pp_b(i),nth_o(i),ph_o(i)
      enddo
      close(21)
      
      


      end subroutine inputvmec

c-----------------------------------------

  
      subroutine FindA(sigmapparang,sigmapperpang,Bc,TpoTb_d,sigmaanis)
      use dstfun
      use pival
      implicit none

      real Bc,TpoTb_d(nbinsrad)
      real fact,val(10),vals(8),pert(8),BoBc,onembc
     $     ,intpol1,intpol2
      integer i,j,k,check
      integer l,nbins
      real c(0:nbinsrad),dd,deltaB
     $     ,sigmapparang(nbinstor,nbinspol,nbinsrad)
     $     ,sigmapperpang(nbinstor,nbinspol,nbinsrad)
     $     ,sigmaanis(nbinsrad)
      real dB,B,Tpp,Tpb,TpoTbloc,sigmaloc
      real val1(1) !debug
      double complex ef(2) !debug


      TpoTb_d=0.
      sigmaanis=0.
      do i=1,nbinsrad
         nbins=0
         do j=1,nbinspolr(i)
            intpol1=mod(twopi*real(j-0.5)/nbinspolr(i)-pi,twopi)
            do k=1,nbinstor
               if(nmarksall(k,j,i).gt.0)then
                  intpol2=twopiPer*(k-0.5)/nbinstor
                  call interp3DRGC(s(i),intpol1,intpol2,val,vals,pert)
                  BoBc=sqrt(val(1))/Bc
                  onembc=1.-BoBc
                  fact=pperpang(k,j,i)/pparang(k,j,i)
                  if(BoBc.ge.1.)then !0.and.BoBc+fact*onembc.gt.0.)then
c                  TpoTbloc=fact/(BoBc+fact*onembc)  
                     call findzeroA(i,j,k,fact,BoBc,TpoTbloc)      
                  
                     sigmaloc=((1./pperpang(k,j,i)
     $                    -onembc/(pparang(k,j,i)
     $                    *(1.+1./BoBc*fact*onembc))))**2.
     $                    *sigmapperpang(k,j,i)**2.
                     sigmaloc=sigmaloc
     $                    +((TpoTbloc*onembc-1.)/pparang(k,j,i))**2.
     $                    *sigmapparang(k,j,i)**2.
                     sigmaloc=sqrt(sigmaloc)*abs(TpoTbloc)
                     
                  else
                     call findzeroA(i,j,k,fact,BoBc,TpoTbloc)         
                     
                  
                     sigmaloc=((1./pperpang(k,j,i)
     $                    -onembc/(pparang(k,j,i)
     $                    *(1.+1./BoBc*fact*onembc))))**2.
     $                    *sigmapperpang(k,j,i)**2.
                     sigmaloc=sigmaloc
     $                    +((TpoTbloc*onembc-1.)/pparang(k,j,i))**2.
     $                    *sigmapparang(k,j,i)**2.
                     sigmaloc=sqrt(sigmaloc)*abs(TpoTbloc)
c$$$               
                  endif

                  if(isnan(TpoTbloc))then
                     sigmaloc=0.
                     TpoTbloc=0.
c     print*,'NO CONTRIBUTION FROM (I,J)=',i,j
                  elseif(TpoTbloc.lt.0.)then
                     sigmaloc=0.
                     TpoTbloc=0.
c     print*,'ANIS NEGATIVE (I,J)=',i,j 
c     print*,'TpoTbloc,fact,BoBc',TpoTbloc,fact,BoBc
                  else
                     nbins=nbins+1
                     TpoTb_d(i)=TpoTb_d(i)+TpoTbloc
                     sigmaanis(i)=sigmaanis(i)+sigmaloc
                  endif
               endif
            enddo
         enddo
         if(nbins.gt.0)then
            TpoTb_d(i)=TpoTb_d(i)/nbins
            sigmaanis(i)=sigmaanis(i)/nbins
         else
            TpoTb_d(i)=0.
            sigmaanis(i)=0.
         endif
      enddo

      end


c--------------------------------------------------------
      subroutine calcM(BoBc,fact,A,M)
      implicit none

      real BoBc,A,M,onembc,fact

      onembc=1.-BoBc

      if(BoBc.ge.1.)then
c         print*,'PROBLEM: B<>Bc',BoBc
c         M=0.
         M=A*(BoBc+fact*onembc)
      else
         M=A*BoBc/((1.-A*onembc)*(1.+A*onembc))
         M=M*((1.+A*onembc)**2.-5.*(A*onembc)**1.5+(A*onembc)**3.5)
         M=M/(1.+A*onembc-2.*(A*onembc)**2.5)
      endif

      end

c--------------------------------------------------------
      subroutine findzeroA(i,j,k,fact,BoBc,A)
      implicit none

      real BoBc,onembc,A,fact
      real f,f1,f2,A0,Ainit,A1,A2,M,fp,eps
      integer check,count,cmax
      integer i,j,k

      A0=0.
      Ainit=100.
      eps=1.e-2
      cmax=100

      onembc=1.-BoBc

c bisection method
      check=0
      count=0
      do while (check.eq.0)
         count=count+1
         A=A0+(count-1)*(Ainit-A0)/cmax
         call calcM(BoBc,fact,A,M)
         f=fact-M
         A1=A+(Ainit-A0)/cmax
         call calcM(BoBc,fact,A1,M)
         f1=fact-M

         if(f*f1.lt.0..or.count.gt.cmax)then
            check=1
         endif
      enddo
      
      if(count.gt.cmax)then
c            write(*,'(a,f5.1)')'A larger than ',A1
c            write(*,'(f5.2,2i3)')BoBc,i,j
            A=-1.
      else
         if(f1.gt.0.)then
            A=Ainit
            A1=0.
         endif

         check=0
         count=0
         do while (check.eq.0)
            count=count+1
            A2=A+0.5*(A1-A)
            call calcM(BoBc,fact,A2,M)
            f2=fact-M
            if(abs(A-A1).lt.eps)then
               A=A2
               check=1
            elseif(f2.gt.0.)then
               A=A2
            else
               A1=A2
            endif
            
            if(count.ge.cmax)then
c               write(*,'(a,2f6.2)')'max iterations: A,f',A2,f2
c               write(*,'(f5.2,2i3)')BoBc,i,j
               A=-1.
               check=1
            endif
         enddo
      endif
      

      end

c-----------------------------------------------------------
      subroutine findcurr
      use particle
      use machine
      use dstfun
      use plasma
      implicit none

      real bkcurr,hotcurr(0:nbinsrad),Zf
      integer i,j,k,Zi
      real dragcoeff,ni,nh(0:nbinsrad)

      
      do i=0,nbinsrad
         nh(i)=0
         do j=1,nbinspolr(i)
            do k=1,nbinstor
               nh(i)=nh(i)+densang(k,j,i)
            enddo
         enddo
      enddo
      nh=nh/nbinspolr/nbinstor
      hotcurr=0.
      dragcoeff=1.
      do i=1,nbinsrad
         bkcurr=j1*s(i)+j2*s(i)**2.+j3*s(i)**3.+j4*s(i)**4.+j5*s(i)**5. ! ohmic total current
         bkcurr=bkcurr*current
         if(dragcurr.eq.1)then!drag current, Graves NF 50,052002 
            Zi=nint(qplasma/Qel)
            ni=density*(n0 +n1*s(i) +n2*s(i)**2. +n3*s(i)**3.
     $           +n4*s(i)**4. +n5*s(i)**5.) ! electron density
            ni=(ni-Zpart*nh(i))/Zi
            if(Zeff.le.0.)then
               Zf=(ni*Zi*Zi+nh(i)*Zpart*Zpart)/(ni*Zi+nh(i)*Zpart) !local Zeff
            else
               Zf=Zeff ! constant input
            endif
            dragcoeff=1-(Zpart/Zf
     $           +(mo*Zi*ni*(1.-Zi/Zf))/(Zpart*ni*mplasma)
     $           -1.46*sqrt(sqrt(s(i))*a00/r00)*2.0   !also: Connor NF 14,185: approx. A(Zeff)=2.0
     $           *(Zpart/Zf-mo*ni*Zi*Zi/(Zpart*Zf*ni*mplasma)))
c$$$            (1.-Zpart/Zf      !drag current, Carlsson PoP 5,2885
c$$$     $           +1.46*sqrt(sqrt(s(i))*a00/r00) !also: Connor NF 14,185
c$$$     $           *(Zpart/Zf-mo*qplasma/(mplasma*Zpart*Qel))*1.5) ! approximate A(Zeff)=1.6f
c$$$            dragcoeff=1.
c$$$            write(15,'(3es14.6)')
c$$$     $           sqrt(s(i)),dragcoeff,currtot(i)*dragcoeff
         endif
         hotcurr(i)=hotcurr(i-1) + currtot(i)*Arad(i)*dragcoeff ! jacobian sign in VMEC
         currtot(i)=hotcurr(i)+bkcurr
      enddo
      end
c-----------------------------------------------------------
      subroutine NaN2Real(ind,V_v,out)
      use dstfun
      use dimarray
      use sinitial
      implicit none
         
      integer ind
      real V_v(nsa),out
      integer is
      real si2,dsi,x

      x=s(ind)
      is=floor(x*nsa)
      si2=is*ds+ds2
      dsi=(si2-x)/ds
      if(dsi.lt.0.)then
         is=is+1
         si2=is*ds+ds2
         dsi=(si2-x)/ds
      endif
      dsi=1.-dsi
      out=dsi*V_v(is+1)+(1.-dsi)*V_v(is)

      end
      
c-----------------------------------------------------------
      function isnan(x)
      implicit none
      
      real x
      integer isnan
      
      isnan=0
      if(2.*x.eq.x.and.x.ne.0.)then
         isnan=1
      elseif(x.ne.x)then
         isnan=1
      endif
      
      end

