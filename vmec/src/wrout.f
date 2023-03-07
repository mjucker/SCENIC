
        subroutine wrout(bsq,gsqrt,bsubu,bsubv,bsubs,br,bphi,bz,lu,lv,
     &                   ppar,pperp,sigma,tau,onembc)
!--------0---------0---------0---------0---------0---------0---------0-c
!
C.. Implicits ..
      implicit none
      include 'name2.inc'
      include 'name3.inc'
        real, dimension(nznt) :: bmod, ppprim, pbprim, densit
        real, dimension(ns,*) :: bsq,gsqrt,bsubu,bsubv,bsubs,ppar,pperp
        real, dimension(ns,*) :: lu,lv,onembc,sigma,tau
        real, dimension(*)    :: br,bphi,bz
!
C.. Intrinsic Functions ..
      intrinsic sqrt, abs
C.. External Calls ..
      real ssum
      external convert, ssum
!
C.. Local Scalars ..
      integer :: js, nmin0, lk,j,k, m, n, mn,iprint,ntskip,nzskip,nl,l,i
     &          ,kpres
      real :: helsym,dmult,gmn,bmn,bsubumn,bsubvmn,bsubsmn,bsupumn,
     &  bsupvmn,tcosi,tsini,fac,zeta,pparmn,pperpmn,sigmamn,taumn
     & ,pme ,pde ,pmh ,pdh ,omtbc   ,optbc ,ppprmn ,pbprmn
     & ,hotdam, hotdmn
      real, dimension(nsd)    :: phi
      real, dimension(mnmax)  :: rmnc,zmns,lmns,xm,xn,bmodmn
!*************
!                 THIS SUBROUTINE CREATES THE BINARY FILE WOUT WHICH 
!                 CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
!                 COEFFICIENTS RMNC,ZMNS,LMNS (FULL MESH)
!*************
        write(8,701)voli,gam,c1p0/nfp,mnmax,nmax1,1,
     &  itfsq,niter/100+1
        kpres = 0
        helsym = 0.
        WRITE(18,40)VOLI,GAM,c1p0/NFP,helsym,MNMAX,NS,MPOL,nmax,NMAX1,1,
     &              ITFSQ,NITER/100+1,kpres
        write(23,711)nfp,mnmax,ns
  701   format(f20.13,f10.3,f20.13,5i6)
  711   format(3i6)
!... CALCULATE RADIAL DERIVATIVES OF HOT PARTICLE PRESSURE TERMS
!... STORE IN ARRAYS pm AND pd PREVIOUSLY USED IN PRESSURE AND EQFOR
        do 5 js=2,ns-1
            pd(js) = ohs * (pres(js+1) * ph(js+1) - pres(js) * ph(js))
            pm(js) = ohs * (tpotb(js+1) - tpotb(js))
    5   enddo
!... INTERPOLATE (EXTRAPOLATE) TO HALF INTEGER MESH
         pdh = c1p5 * pd(2) - cp5 * pd(3)
         pmh = c1p5 * pm(2) - cp5 * pm(3)
         pde = c1p5 * pd(ns-1) - cp5 * pd(ns-2)
         pme = c1p5 * pm(ns-1) - cp5 * pm(ns-2)
       do 6 js=ns-2,2,-1
         pd(js+1) = cp5 * (pd(js+1) + pd(js)) / (pres(js+1)*ph(js+1))
         pm(js+1) = cp5 * (pm(js+1) + pm(js))
    6  enddo
         pd(2)  = pdh / (pres(2)*ph(2))
         pd(ns) = pde / (pres(ns)*ph(ns)) 
         pm(2)  = pmh
         pm(ns) = pme
!ALTERNATE EXTRAPOLATION
         pd(2) = c2p0*pd(3) - pd(4)
         pd(ns) = c2p0*pd(ns-1) - pd(ns-2) 
        WRITE(19,51) ohs
 51     FORMAT(/,' PRESSURE ARRAYS;  OHS =',1p1e14.6)
        WRITE(19,50)(PD(js),PM(js),PRES(js),PH(js),TPOTB(js),js=1,NS)
!*VDIR noassume
        do 29 js=1,ns
        call convert(rmnc,zmns,lmns,xm,xn,js,xc,xc(1+mns),
     &  xc(1+2*mns),xc(1+3*mns),xc(1+4*mns),xc(1+5*mns))
        mn = 0
        hotdam = pres(js) * ph(js) / sqrt(tpotb(js))
        do 10 lk=1,nznt
           bmod(lk)=sqrt(c2p0*abs(bsq(js,lk)/bscale**2))
           omtbc = c1p0 - tpotb(js) * onembc(js,lk)
           optbc = c1p0 + tpotb(js) * onembc(js,lk)
        if (onembc(js,lk) .le. czero) then
           densit(lk)= (ppar(js,lk) - pres(js))*hotdam/(pres(js)*ph(js))
           pbprim(lk) =  (ppar(js,lk) -pres(js)) *
     &             (pd(js) + onembc(js,lk) * pm(js) / omtbc)
           ppprim(lk) =  (pperp(js,lk)-pres(js)) *
     &             (pd(js) + optbc * pm(js) / (omtbc * tpotb(js)))
        else
          densit(lk) = hotdam * (c1p0 - onembc(js,lk)) *
     &    (optbc - c2p0*(tpotb(js)*onembc(js,lk))**1.5) / (omtbc*optbc)
          pbprim(lk) =  ((ppar(js,lk) -pres(js)) * pd(js) +
     &    ( c2p0 * tpotb(js) * onembc(js,lk)**2 * (ppar(js,lk)-pres(js))
     &   + pres(js)*ph(js)*(c1p0-onembc(js,lk))*onembc(js,lk)*(c1p0 -5.0
     &        *(tpotb(js)*onembc(js,lk))**1.5))* pm(js) / (omtbc*optbc))
          ppprim(lk) =  ((pperp(js,lk)-pres(js)) * pd(js) +
     & ((pperp(js,lk)-pres(js))*(c1p0+c3p0*(tpotb(js)*onembc(js,lk))**2)
     &   / tpotb(js) + pres(js)*ph(js)*tpotb(js)*(c1p0-onembc(js,lk))**2
     &   * onembc(js,lk)*(c2p0*optbc-sqrt(tpotb(js)*onembc(js,lk))*(7.5
     &   - 3.5*(tpotb(js)*onembc(js,lk))**2))/(omtbc*optbc))
     &   * pm(js)/ (omtbc * optbc))
        endif
   10   end do
        do 25 m = 0,mpol1
        nmin0 = -nmax
        if(m.eq.0)nmin0 = 0
        do 20 n = nmin0,nmax
        mn = mn+1
!        dmult = c2p0*dnorm/(mscale(m)*nscale(iabs(n)))
        dmult = c2p0*dnorm/(mscale(m)*nscale(abs(n)))
        if(m.eq.0.and.n.eq.0)dmult=cp5*dmult
        gmn = czero
        bmn = czero
        bsubumn = czero
        bsubvmn = czero
        bsubsmn = czero
        bsupumn = czero
        bsupvmn = czero
        pperpmn = czero
        pparmn  = czero
        taumn   = czero
        sigmamn = czero
        ppprmn  = czero
        pbprmn  = czero
        hotdmn  = czero
        if(js.ne.1)then
        do 201 j = 1,ntheta2
        do 200 k = 1,nzeta
        lk = k + nzeta*(j-1)
        l  = js + ns*(lk-1)
        if(n.ge.0)then
        tcosi = dmult*(cosmui(j,m)*cosnv(k,n)+sinmu(j,m) *sinnv(k,n))
        tsini = dmult*(sinmu (j,m)*cosnv(k,n)-cosmui(j,m)*sinnv(k,n))
        else
        tcosi = dmult*(cosmui(j,m)*cosnv(k,-n)-sinmu(j,m) *sinnv(k,-n))
        tsini = dmult*(sinmu (j,m)*cosnv(k,-n)+cosmui(j,m)*sinnv(k,-n))
        endif
        bmn = bmn + tcosi*bmod(lk)
        gmn = gmn + tcosi*gsqrt(js,lk)
        bsubumn = bsubumn + tcosi*bsubu(js,lk)
        bsubvmn = bsubvmn + tcosi*bsubv(js,lk)
        bsubsmn = bsubsmn + tsini*bsubs(js,lk)
        bsupumn = bsupumn + tcosi*lv(js,lk)
        bsupvmn = bsupvmn + tcosi*lu(js,lk)
        pperpmn = pperpmn + tcosi*pperp(js,lk)  !CAREFUL: total pressure
        pparmn  = pparmn  + tcosi*ppar(js,lk)   !CAREFUL: total pressure
!        pperpmn = pperpmn + tcosi*(pperp(js,lk)-pres(js))  !CAREFUL: hot particle pressure only
!        pparmn  = pparmn  + tcosi*(ppar(js,lk)-pres(js))   !CAREFUL: hot particle pressure only
        sigmamn = sigmamn + tcosi*sigma(js,lk)
        taumn   = taumn   + tcosi*tau(js,lk)
        pbprmn  = pbprmn  + tcosi*pbprim(lk)
        ppprmn  = ppprmn  + tcosi*ppprim(lk)
        hotdmn  = hotdmn  + tcosi*densit(lk)
 200    end do
 201  end do
        if(js.eq.ns/2)bmodmn(mn) = bmn
        end if
        WRITE(18,50)XM(MN),XN(MN),RMNC(MN),ZMNS(MN),GMN
        write(18,50)sigmamn, taumn, pbprmn, ppprmn
        write(23,50)xm(mn),xn(mn),rmnc(mn),zmns(mn),bmn,pparmn,pperpmn
     &             ,hotdmn
!        write(23,50)xm(mn),xn(mn),rmnc(mn),zmns(mn),bmn,-pbprmn,-ppprmn
!     &             ,taumn
 20     write(8,703)xm(mn),xn(mn),rmnc(mn),zmns(mn),lmns(mn),
     &  bmn,gmn,bsubumn/bscale,bsubvmn/bscale,bsubsmn/bscale,
     &  bsupumn/bscale,bsupvmn/bscale
 25     end do
 29   end do
  703    format(4e20.13)
        phi(1) = czero
        do 30 js = 2,ns
           phi(js) = twopi*hs*ssum(js-1,phips(2),1)
   30   end do
        fac = abs(bscale)**(gam-c2p0)
        write(8,703)(iotas(js),mass(js)*fac,pres(js),phips(js),buco(js),
     &  bvco(js),phi(js),vp(js),ju(js),jv(js),specw(js),js=2,ns)
        write(8,703)(fsqt(i),wdot(i),i=1,100)
        WRITE(18,50)(IOTAS(js),MASS(js),PRES(js),-PHIPS(js),VP(js),
     &               js=1,NS)
        write(44,501)PHIPS(2)
        write(44,502)(iotas(js),bvco(js),js=2,ns)
 40     FORMAT(1x,1pe22.12,1x,1p3e12.5,8i4,i1)
! 50     FORMAT(1x,1p5e14.6)
 501    format(1p1e22.14)
 502    format(1p2e22.14)
 50     FORMAT(1x,1p5e22.14)
        if(nvac.eq.0)return                                             vac
        write(3,60)rbtor/bscale,ctor/bscale,bscale                      vac
 60     format(/,' RBTOR = ',1pe16.8,' NET TOROIDAL CURRENT = ',        vac
     &  1pe16.8,' BSCALE = ',1pe16.8)                                   vac
        do 69 iprint=1,2                                                vac
        if(iprint.eq.1)write(3,70)                                      vac
        if(iprint.eq.2)write(3,75)                                      vac
        ntskip=1+ntheta1/12                                             vac
        nzskip=1+nzeta/6                                                vac
        do 67 l=1,nzeta,nzskip                                          vac
        zeta = real(360*(l-1))/real(nzeta)                              vac
        do 65 k=1,ntheta2,ntskip                                        vac
        lk=l+nzeta*(k-1)                                                vac
        nl = ns*lk                                                      vac
        if(iprint.eq.1)write(3,90)zeta,r1(nl)+r1(nl+nrzt),              vac
     &  z1(nl)+z1(nl+nrzt),(bsqsav(lk,n)/bscale**2,n=1,3),              vac
     &  bsqvac(lk)/bscale**2                                            vac
        if(iprint.eq.2)write(3,95)zeta,r1(nl)+r1(nl+nrzt),              vac
     &  z1(nl)+z1(nl+nrzt),(1.5*br(nl) - 0.5*br(nl-1))/bscale,          vac
     &  (1.5*bphi(nl) - 0.5*bphi(nl-1))/bscale,                         vac
     &  (1.5*bz(nl) - 0.5*bz(nl-1))/bscale,brv(lk)/bscale,              vac
     &  bphiv(lk)/bscale,bzv(lk)/bscale                                 vac
 65     end do                                                          vac
 67    end do                                                           vac
 69   end do                                                            vac
 70     format(/,4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,                      vac
     &  'BSQMHDI',5x,'BSQVACI',5x,'BSQMHDF',5x,'BSQVACF',/)             vac
 75     format(/,4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,                      vac
     &  'BR',8x,'BPHI',6x,'BZ',8x,'BRv',7x,'BPHIv',5x,'BZv',/)          vac
 90     format(1pe10.2,1p6e12.4)                                        vac
 95     format(1pe10.2,1p2e12.4,1p6e10.2)                               vac
        write(3,100)                                                    vac
 100    format(//,3x,'mb',2x,'nb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)',     vac
     &  6x,'mb',2x,'nb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)'/)              vac
        do 110 mn=1,mnmax,2                                             vac
        write(3,115)nint(xm(mn)),nint(xn(mn)/nfp),rmnc(mn),zmns(mn),    vac
     &  bmodmn(mn),nint(xm(mn+1)),nint(xn(mn+1)/nfp),rmnc(mn+1),        vac
     &  zmns(mn+1),bmodmn(mn+1)                                         vac
  110   end do                                                          vac
 115    format(i5,i4,1p3e12.4,3x,i5,i4,1p3e12.4)                        vac
        write(3,120)                                                    vac
 120    format(/,3x,'mf',2x,'nf',5x,'potvacs',6x,'mf',2x,'nf',5x,       vac
     &  'potvacs',6x,'mf',2x,'nf',5x,'potvacs'/)                        vac
        do 130 mn=1,mpmax,3                                             vac
        write(3,135)nint(xmpot(mn)),nint(xnpot(mn)/nfp),                vac
     &  potvac(mn)/bscale,nint(xmpot(mn+1)),nint(xnpot(mn+1)/nfp),      vac
     &  potvac(mn+1)/bscale,nint(xmpot(mn+2)),nint(xnpot(mn+2)/nfp),    vac
     &  potvac(mn+2)/bscale                                             vac
  130   end do                                                          vac
 135    format(i5,i4,1pe12.4,3x,i5,i4,1pe12.4,3x,i5,i4,1pe12.4)         vac
        return
        end subroutine wrout
