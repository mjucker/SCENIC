c--------0---------0---------0---------0---------0---------0---------0-c
c---TRANSFORMS VMECAP OUTPUT FROM CYLINDRICAL COORDINATE            ---c
c---INFORMATION TO CARTESIAN COORDINATES. THE OUTPUT CORRESPONDS TO ---c
c---THE CARTESIAN COORDINATES OF EACH GRID POINT OF A FLUX SURFACE. ---c
c---'b' CONTAINS mod-B, p CONTAINS p_par and u CONTAINS p_perp      ---c
c---'d' CONTAINS THE HOT PARTICLE DENSITY                           ---c
c---            Last Modification:   15.12.05                       ---c
c----------------------------------------------------------------------c
c        parameter(nit=48,njt=73,nkt=32,npert=4,mnt=50,nkst=nkt*npert+1
c        parameter(nit=24,njt=65,nkt=16,npert=1,mnt=50,nkst=nkt*npert+1
c     >            ,njkt=(njt-1)*nkt,njkst=(njt-1)*nkt*npert)
        program vma3dp
c
        include 'plovma.inc'
!
C.. Implicits ..
!      implicit none
c
        real, dimension(njkt,nit)       :: rt  , zt  , bt  , pt  , ut
     &                                    ,dn
        real, dimension(nkst,0:nit,njt) :: x   , y   , z   , b   , p , u
     &                                    ,d
        real, dimension(2)              :: xax , yax , zax , pax , bax
     &                                    ,uax ,dax
        real, dimension(mnt,0:nit) :: rmn , zmn , bmn , pparmn , ppermn
     &                                                         , hotdmn
        real, dimension(mnt)       :: xm  , xn
!
!
C.. Internal scalars
      integer :: nimax, nskipr, nskipp, nskipt, kstart, kstop, nistop  ,
     &           ni   , nj    , nk    , i     , j     , k    , nst,
     &           ks   , jk    , jks   , kf    , nper  , nrp  , npp, ntp,
     &           mns  , njk
      real    :: argc , arg   , argp  , cosp  , sinp  , dt   , twopif
!
C.. External Calls ..
      external realsp
!
C.. Intrinsic Functions ..
      intrinsic MOD, SIN, cos, atan, min
!
        nimax  = 192
        nskipr = 1 
        nskipp = 1
        nskipt = 1
        kstart = 1
        kstop  = nkst
        mns    = mnt
        ni     = nit
        nj     = njt - 1
        nk     = nkt
        njk    = nj*nk
c
        read (5,1000) nimax,nskipr,nskipp,nskipt,kstart,kstop
c        write (6,1000) nimax,nskipr,nskipp,nskipt,kstart,kstop
 1000   format (/6i7)
c
 1001   format (5i3)
c
        call realsp(ns,ni,nj,nk,njk,mns,nper,rt,zt,bt,pt,ut,dn)
        if (ns.gt.ni .or. nper.gt.npert) stop
c
        if (mod(nj,nskipp).ne.0) stop
 1002   format (1p3e16.8,/1p5e16.8)
c
c--    COMPUTE QUANTITIES ON CARTESIAN GRID.                     ------c
        twopif = 8. * atan(1.) / nper
c
        do 305 i=1,ns
          do 303 nst=1,nper
           argc = twopif * (nst-1.)
            do 302 k=1,nk
              ks = (nst-1) * nk + k
              argp = twopif * (k-1)/nk
               do 300 j=1,nj
                 jk = (k-1)*nj + j
                 jks = (nst-1)*nj*nk + jk
                 arg  = argp + argc
                 cosp = cos(arg)
                 sinp = sin(arg)
                 x(ks,i,j)  = rt(jk,i) * cosp
                 y(ks,i,j)  = rt(jk,i) * sinp
                 z(ks,i,j)  = zt(jk,i)
                 p(ks,i,j)  = pt(jk,i)
                 u(ks,i,j)  = ut(jk,i)
                 d(ks,i,j)  = dn(jk,i)
                 b(ks,i,j)  = bt(jk,i)
 300           end do
 302        end do
 303      end do
 305   end do
c
c--- EXTRAPOLATE TO THE MAGNETIC AXIS                               ---c
c
        kf = nper * nk + 1
        dt = 1./nj
c
        do 350 ks=1,kf-1
          do 320 i=1,2
            xax(i) = 0.
            yax(i) = 0.
            zax(i) = 0.
            pax(i) = 0.
            uax(i) = 0.
            dax(i) = 0.
            bax(i) = 0.
 320      end do
          do 335 i=1,2
            do 330 j=1,nj
              xax(i)  = xax(i)  + x(ks,i,j)
              yax(i)  = yax(i)  + y(ks,i,j)
              zax(i)  = zax(i)  + z(ks,i,j)
              pax(i)  = pax(i)  + p(ks,i,j)
              uax(i)  = uax(i)  + u(ks,i,j)
              dax(i)  = dax(i)  + d(ks,i,j)
              bax(i)  = bax(i)  + b(ks,i,j)
 330        end do
 335      end do
             x(ks,0,1) = (1.5*xax(1)  - 0.5*xax(2) ) * dt
             y(ks,0,1) = (1.5*yax(1)  - 0.5*yax(2) ) * dt
             z(ks,0,1) = (1.5*zax(1)  - 0.5*zax(2) ) * dt
             p(ks,0,1) = (1.5*pax(1)  - 0.5*pax(2) ) * dt
             u(ks,0,1) = (1.5*uax(1)  - 0.5*uax(2) ) * dt
             d(ks,0,1) = (1.5*dax(1)  - 0.5*dax(2) ) * dt
             b(ks,0,1) = (1.5*bax(1)  - 0.5*bax(2) ) * dt
          do 345 j=2,nj
            x(ks,0,j) = x(ks,0,1)
            y(ks,0,j) = y(ks,0,1) 
            z(ks,0,j) = z(ks,0,1) 
            p(ks,0,j) = p(ks,0,1)   
            u(ks,0,j) = u(ks,0,1)   
            d(ks,0,j) = d(ks,0,1)   
            b(ks,0,j) = b(ks,0,1) 
 345      end do
 350    end do
c
c---    CLOSE THE TORUS POLOIDALLY AND TOROIDALLY                ------c
c
        do 405 i=0,ns
           do 402 ks=1,kf-1
             x(ks,i,nj+1)  = x(ks,i,1)
             y(ks,i,nj+1)  = y(ks,i,1)
             z(ks,i,nj+1)  = z(ks,i,1)
             p(ks,i,nj+1)  = p(ks,i,1)
             u(ks,i,nj+1)  = u(ks,i,1)
             d(ks,i,nj+1)  = d(ks,i,1)
             b(ks,i,nj+1)  = b(ks,i,1)
 402       end do
c      
           do 404 j=1,nj+1
             x(kf,i,j)  = x(1,i,j)
             y(kf,i,j)  = y(1,i,j)
             z(kf,i,j)  = z(1,i,j)
             p(kf,i,j)  = p(1,i,j)
             u(kf,i,j)  = u(1,i,j)
             d(kf,i,j)  = d(1,i,j)
             b(kf,i,j)  = b(1,i,j)
 404       end do
 405    end do
c
c--     WRITE OUTPUT.                                           -------c
c
        nistop = min(ns,nimax)
        nrp = 1 + nistop / nskipr
        npp = 1 + nj / nskipp  
        ntp = 1 + (kstop - kstart)/nskipt
        write (45,1006) ntp, nrp, npp
        write (46,1006) ntp, nrp, npp
        write (47,1006) ntp, nrp, npp
        write (48,1006) ntp, nrp, npp
 1006   format (3i4)
        xx = 0.
        xy = 0.
        xz = 0.
        do 505 ks = kstart,kstop,nskipt
           do 503 i=0,nistop,nskipr
              do 500 j=1,nj+1,nskipp
cbaspl           write (6,1007) 1+(ks-kstart)/nskipt,i/nskipr,1+(j-1)/nskipp
cbaspl     >       ,x(ks,i,j),y(ks,i,j),z(ks,i,j),p(ks,i,j),b(ks,i,j)
cbaspl     >       ,xx,xy,xz
                 write (45,1008) x(ks,i,j),y(ks,i,j),z(ks,i,j),b(ks,i,j)
                 write (46,1008) x(ks,i,j),y(ks,i,j),z(ks,i,j),p(ks,i,j)
                 write (47,1008) x(ks,i,j),y(ks,i,j),z(ks,i,j),u(ks,i,j)
                 write (48,1008) x(ks,i,j),y(ks,i,j),z(ks,i,j),d(ks,i,j)
 500          end do
 503       end do
 505    end do
c
 1007   format(3i4,1p4e16.8,/1p4e16.8)
 1008   format(1p4e16.8)
        stop
        end program vma3dp
C-----------------------------------------------------------------------
C  04.12.91       LAST MODIFICATION 15.12.05        UHS:  realsp1 DATA D
C----------------------------------------------------------------------C
C
        subroutine realsp(ns,ni,nj,nk,njk,mns,nper,r,z,b,p,u,d)
C
C----------------------------------------------------------------------C
c       RECONSTRUCTION IN REAL SPACE OF THE GEOMETRY, MAGNITUDE OF THE
c       MAGNETIC FIELD AND THE TOTAL PRESSURE.
C----------------------------------------------------------------------C
C
c
        include 'plovma.inc'
c
        real, dimension(mnt,0:nit) :: rmn , zmn , bmn , pparmn , ppermn
     &                                                         , hotdmn
        real, dimension(njkt,mnt ) :: tcos, tsin
        real, dimension(njkt,*)    :: r   , z   , b   , p      , u
     &                                    , d
        real, dimension(mnt)       :: xm  , xn
C
C.. Local Scalars ..
      integer ::  nper, mnmax, nsp, ns, mns, i, mn, k, j, jk, iarg, lk
      real :: twopi, dth, dph, arg
C
C.. Intrinsic Functions ..
      intrinsic atan, cos, sin
C
C
      twopi = 8. * atan(1.)
      read (23,711) nper,mnmax,nsp
 711  format(3i6)
      if (mnmax.gt.mns) stop
      ns = nsp - 1
      dth = 1. / nj
      dph = 1. / (nk * nper)
c
      do 10 i = 0,ns
        do 8 mn = 1,mnmax
!          read (23,710) xm(mn),xn(mn),rmn(mn,i),zmn(mn,i),bmn(mn,i)
!     &                 ,pparmn(mn,i),ppermn(mn,i)
          read (23,710) xm(mn),xn(mn),rmn(mn,i),zmn(mn,i),bmn(mn,i)
          read (23,720) pparmn(mn,i),ppermn(mn,i),hotdmn(mn,i)
c jucker Fourier amplitudes for shaping coefficients
          write(49,721) xm(mn),xn(mn),rmn(mn,i),zmn(mn,i)
 8      end do
 10   end do
 710  format (1x,1p5e22.14)
 720  format (1x,1p3e22.14)
 721  format (1x,1p4e22.14)
c
      do 25 mn = 1,mnmax
         do 22 k = 1,nk
           do 20 j = 1,nj
             jk = j + (k-1) * nj
             arg = xm(mn) * (j-1) * dth  -  xn(mn) * (k-1) * dph
             iarg = arg
             arg = twopi * (arg - iarg)
             tcos(jk,mn) = cos(arg)
             tsin(jk,mn) = sin(arg)
 20        end do
 22      end do
 25   end do
c
      do 50 i = 1,ns
        do 30 lk = 1,njk
            r(lk,i) = 0.
            z(lk,i) = 0.
            p(lk,i) = 0.
            u(lk,i) = 0.
            d(lk,i) = 0.
            b(lk,i) = 0.
 30     end do
c
        do 45 mn = 1,mnmax
           do 40 lk = 1,njk
      r(lk,i) = r(lk,i) + 0.5*(rmn(mn,i)+rmn(mn,i-1))   * tcos(lk,mn)
      z(lk,i) = z(lk,i) + 0.5*(zmn(mn,i)+zmn(mn,i-1))   * tsin(lk,mn)
      p(lk,i) = p(lk,i) + pparmn(mn,i)                  * tcos(lk,mn)
      u(lk,i) = u(lk,i) + ppermn(mn,i)                  * tcos(lk,mn)
      d(lk,i) = d(lk,i) + hotdmn(mn,i)                  * tcos(lk,mn)
      b(lk,i) = b(lk,i) + bmn(mn,i)                     * tcos(lk,mn)
 40        end do
 45     end do
 50   end do
c
      return
      end subroutine realsp
c
C----------------------------------------------------------------------C
