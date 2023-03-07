      subroutine interp3DRGC(s,theta,xi,val,vals,pert)
      use particle
      use terp
      use arrays
      use sinitial
      implicit none

      integer  kxi,kxi1,is,jth,is1,jth1,l,kxim1
      REAL  twopi,dxi,dth,xiperiode,thperiode,s,theta,xi,
     >    xik,xik1,si2,
     >     twopiPer,thetaj,dsi,dthj,f11,f12,f21,f22,xikm1,pxi,
     >     smax,
     >     coef(4),br(3),br2(3),val(10),vals(8),pert(8)


      smax=1.00-1./nsa/2.  

      twopi=4.*acos(0.)
      twopiPer=twopi/Lper

      if (s.lt.0.) then
         s=abs(s)
         theta=theta+twopi/2.
      endif
      


      dxi=twopiPer/nxia
      dth=twopi/ntha

c      smax=starts+1.*nsat/nsa

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


c controle de la taille du maillage (s-theta)
      is1=is+1
      if (is1-sinind.gt.nsat.or.(is-sinind.lt.1.and.starts.ne.0)) then
         if(is.lt.nsa)then
            print*,'s,is,is1',s,is-sinind,is1-sinind
            print*,'ADJUST NSAT TO HIGHER VALUE'
            stop
         else
            is=nsa-1
         endif
      endif

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

c     if(abs(s-0.1355).lt.0.0005)print*,s,dsi,si2



      coef(1) = pxi*(pxi-1)/2.
      coef(2) = (1.0-pxi*pxi)
      coef(3) = pxi*(pxi+1)/2.

c interpolation 2-D lineaire sur les plan xik et xik1 et xikm1
      
      is=is-sinind
      is1=is+1


      do l=1,10
         br(1)=
     >        f11 * tabarray(jth,kxi,is,l)+f22*tabarray(jth1,kxi,is1,l)
     >        +f21 *tabarray(jth,kxi,is1,l)+f12*tabarray(jth1,kxi,is,l)
         br(2)=
     >        f11*tabarray(jth,kxi1,is,l)+f22*tabarray(jth1,kxi1,is1,l)
     >        +f21*tabarray(jth,kxi1,is1,l)+f12*tabarray(jth1,kxi1,is,l)
         br(3)=
     >       f11*tabarray(jth,kxim1,is,l)+f22*tabarray(jth1,kxim1,is1,l)
     >      +f21*tabarray(jth,kxim1,is1,l)+f12*tabarray(jth1,kxim1,is,l)
c interpolation 2-D lineaire entre xik et xik1 et xikm1

         val(l)=coef(1)*br(3) + coef(2)*br(1) +
     >        coef(3)*br(2)
      enddo

c jucker perturbations
      do l=1,8
         br(1)=
     >        f11 * pertarray(jth,kxi,is,l)
     $        +f22*pertarray(jth1,kxi,is1,l)
     >        +f21 *pertarray(jth,kxi,is1,l)
     $        +f12*pertarray(jth1,kxi,is,l)
         br(2)=
     >        f11*pertarray(jth,kxi1,is,l)
     $        +f22*pertarray(jth1,kxi1,is1,l)
     >        +f21*pertarray(jth,kxi1,is1,l)
     $        +f12*pertarray(jth1,kxi1,is,l)
         br(3)=
     >       f11*pertarray(jth,kxim1,is,l)
     $        +f22*pertarray(jth1,kxim1,is1,l)
     >      +f21*pertarray(jth,kxim1,is1,l)
     $        +f12*pertarray(jth1,kxim1,is,l)
c interpolation 2-D lineaire entre xik et xik1 et xikm1

         pert(l)=coef(1)*br(3) + coef(2)*br(1) +
     >        coef(3)*br(2)
         
      enddo


      if (val(1).le.0.0) then
         write(6,*) 'interp3DRGC, STOP. B<0',val(1),s
c isaev
         stop

c        val(1)=tabarray(jth,kxi,is,1)
      endif

      coef(4)=1.0-dsi
      do l=1,8
         vals(l) =coef(4)*tabarrays(is,l) + dsi*tabarrays(is1,l)
      end do
      

      end
