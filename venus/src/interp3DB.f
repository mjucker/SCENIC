      function interp3DB(s,theta,xi,smax3)
      use particle
      use arrays
      use terp
      use sinitial
      implicit none

c      include 'include.dir/terp.inc'
c      include 'include.dir/dimarray.inc'
c      include 'include.dir/particule.inc'

      integer  kxi,kxi1,is,jth,is1,jth1,kxim1
      REAL  twopi,dxi,dth,xiperiode,thperiode,s,theta,
     >     xi,xik,xik1,si2,
     >     interp3DB,twopiPer,thetaj,dsi,dthj,f11,f12,f21,f22,
     >     smax,xikm1,pxi,coef(3),br(3)
      real smax3

c      allocate(tabarray(ntha+1,nxia+1,nsat,4))

      smax=1.00-1.0/nsa/2.0 
c      smax=1.00 

      twopi=4.*acos(0.)
      twopiPer=twopi/Lper

      if (s.le.0.) then
         s=abs(s)

       print*,'s,theta,xi',s,theta,xi
        print*,'interp3DB, s.le.0!!!'
         stop

         theta=theta+twopi/2.
      endif

     
      if(abs(s-sinit).gt.smax3)smax3=abs(s-sinit)
      dxi=twopiPer/nxia
      dth=twopi/ntha

      smax=(sinind+nsat)*ds

      xiperiode=mod(mod(xi,twopiPer)+twopiPer,twopiPer) ! xi --> [0,2pi/L]
      thperiode=mod(mod(theta,twopi)+twopi,twopi)    ! theta --> [0,2pi]

c determination de xii et xi1 tq xi in [xii,xi1] 
      kxi=nint(xiperiode/dxi)+1
      if (kxi.eq.(nxia+1)) then
         kxi1=kxi
         kxim1=kxi-2
         kxi=kxi-1
      elseif (kxi.eq.1) then
         kxi1=kxi+2
         kxim1=kxi
         kxi=kxi+1
      else
         kxi1=kxi+1
         kxim1=kxi-1
      endif

      xik=(kxi-1)*dxi
      xik1=(kxi1-1)*dxi
      xikm1=(kxim1-1)*dxi
      pxi=(xiperiode-xik)/dxi

c determination de la boite pour interpolation 2-D dans plan (s,u)
c c'eat-a-dire les (si2,thetaj) et DSi et DTHj

      is=floor(s*nsa)
      jth=nint(thperiode/dth)+1
      

      si2=is*ds+ds2
      dsi=(si2-s)/ds
      if(dsi.lt.0.)then
         is=is+1
         si2=is*ds+ds2
         dsi=(si2-s)/ds
      endif
      dsi=1-dsi

c controle de la taille du maillage
      is1=is+1
      if (is1-sinind.gt.nsat.or.(is-sinind.lt.1.and.starts.ne.0)) then
         print*,'s,is,is1',s,is-sinind,is1-sinind
         print*,'ADJUST NSAT TO HIGHER VALUE'
         stop
      endif

c      if (s.gt.(smax-ds2)) then
c         is=nsat-1
c         is1=nsat
c      endif

      jth1=jth+1
      if (jth.eq.ntha+1) then
         jth1=ntha+1
         jth=ntha
      endif

      thetaj=(jth-1)*dth
      dthj=(thperiode-thetaj)/dth
      f22=dsi*dthj
      f21=dsi-f22
      f12=dthj-f22
      f11=1+f22-(dsi+dthj)




      coef(1) = pxi*(pxi-1)/2.
      coef(2) = (1.0-pxi*pxi)
      coef(3) = pxi*(pxi+1)/2.      
c interpolation 2-D lineaire sur les plan xik et xik1 et xikm1

      if(jth.gt.ntha+1)print*,'jth.gt.ntha+1'


      if(kxi1.gt.nxia+1.or.kxim1.gt.nxia+1)then
         print *,'kxi1.gt.nxia+1'
         stop
      endif

      is=is-sinind
      is1=is+1

      br(1)=
     >        f11 * tabarray(jth,kxi,is,1)+f22*tabarray(jth1,kxi,is1,1)
     >        +f21 *tabarray(jth,kxi,is1,1)+f12*tabarray(jth1,kxi,is,1)
      br(2)=
     >        f11*tabarray(jth,kxi1,is,1)+f22*tabarray(jth1,kxi1,is1,1)
     >        +f21*tabarray(jth,kxi1,is1,1)+f12*tabarray(jth1,kxi1,is,1)

      br(3)=
     >       f11*tabarray(jth,kxim1,is,1)+f22*tabarray(jth1,kxim1,is1,1)
     >      +f21*tabarray(jth,kxim1,is1,1)+f12*tabarray(jth1,kxim1,is,1)




c interpolation 2-D lineaire entre xik et xik1 et xikm1

      interp3DB=coef(1)*br(3) + coef(2)*br(1) +
     >        coef(3)*br(2)
      
      if (interp3DB.le.0.0) then
          write(6,*) 'interp3DB, STOP. B<0',interp3DB,s
          stop
c isaev
c        interp3DB=sqrt(abs(interp3DB))
      else
         interp3DB=sqrt(interp3DB)
      endif
         
      end
