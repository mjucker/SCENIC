C--------0---------0---------0---------0---------0---------0---------0-c
C  13.01.06          LAST MODIFICATION 31.03.06       UHS: MERCIE DATA D
C--------0---------0---------0---------0---------0---------0---------0-c
      subroutine mercap(ni,njk,lmnb,lmnb0,nper,dtdp,s,ftp,fpp,ftpp,fppp,
     &                  ci,cj,pp,cjp,cip,mb,nb,tcos,tsin,gssu,gssl,gstl,
     &                  fbjac,bjac,bjacs,bsq,fr,fz,fbsq,fpark,fparjp,
     &                  fcurvs,fcurvn,flshea,fgssu,fbalcp,fbalcs,fbalcq,
     &                  fbaldp,fbalds,fbalnp,fbalns,fbalnq,fbaldn,
     &                  fbaldt,parkur,parjp,sigmab,taub,sigbs,gparp,
     &                  gperp,glm,bgradg,sigbxg,baldp,balds,balcp,balcs,
     &                  balcq,balnp,balns,balnq,baldn,curvds,curvdn,
     &                  shearl,vp)
C
C
C.. Implicits ..
      implicit none
C
C.. Parameters ..
      real tiny
      parameter (tiny = 1e-300)
C
C.. Formal Arguments ..
      integer                   :: njk,ni,lmnb,lmnb0,nper
      integer, dimension(lmnb)  :: mb    ,nb
      real                      :: dtdp
!
      real, dimension(lmnb)     :: fbalcp ,fbalcs ,fbalcq ,fbaldp,fbalds
     &                            ,fbalnp ,fbalns ,fbalnq ,fbaldn,fbaldt
     &                            ,flshea ,fcurvs ,fcurvn ,fgssu
      real, dimension(lmnb,0:ni):: fbjac
      real, dimension(lmnb,ni)  :: fpark  ,fparjp ,fr     ,fz    ,fbsq
      real, dimension(njk,lmnb) :: tcos   ,tsin
      real, dimension(njk,0:ni) :: sigbs  ,bjacs  ,bjac   ,bsq
     &                            ,gssu   ,gssl   ,gstl
      real, dimension(njk,ni)   :: parkur,parjp ,sigmab,taub  ,gparp
     &                            ,gperp
      real, dimension(njk)      :: bgradg,sigbxg,baldp ,balds ,balcp
     &                            ,balcs ,balcq ,shearl,curvds,curvdn
     &                            ,baldn ,balnp ,balns ,balnq
      real, dimension(ni+1)     :: ftp   ,fpp   ,fppp  ,ftpp
      real, dimension(ni)       :: ci    ,cj    ,cjp   ,cip   ,pp    ,vp
      real, dimension(0:ni)     :: s
      real, dimension(ni,4)     :: glm
C
C.. Local Scalars ..
      integer ::   i,    jk,     l
      real :: dper  ,dsca  ,dsec  ,hdsc,  hdsq  ,hq    ,fintls,tbalin,
     &        dperni,dscani,dsecni,hdscni,hdsqni,hqni  ,vp0   ,b0r0  ,
     &        t0    ,t1    ,t2    ,t3    ,qsurf ,dpamp ,dsamp ,cpamp ,
     &        csamp ,cqamp ,cuamp ,cwamp ,dnamp ,bdtemp,dsampk,dsampn,
     &        sigtau   ,pgradt    ,sqgbsq,pmnjac,cmnjac,gspl  ,sqrbsq
C
C.. Intrinsic Functions ..
      intrinsic abs, sqrt, merge
C
C ... Executable Statements ...
C
C
! --- OUTPUT FOR BALLOONING CALCULATION
      b0r0 = fbsq (lmnb0, ni) / (fr (lmnb0, ni) * fr (lmnb0, ni) )
      vp0 = 1.5 * vp (1) - 0.5 * vp (2)
      write (25, 850) ni, lmnb, b0r0, vp0
      write (25, 860) (mb (l), l = 1, lmnb), (nb (l), l = 1, lmnb)
      write (26, 850) ni, lmnb, b0r0, vp0
      write (26, 860) (mb (l), l = 1, lmnb), (nb (l), l = 1, lmnb)
      write(66,858)nper,lmnb,ni-2
      write(64,858)nper,lmnb,ni-2
!
      do 700 i = 2,ni-1
C
C..SURFACE QUANTITIES
        sqgbsq = fpp(i) * cj(i) - ftp(i) * ci(i)
        qsurf  = ftp(i) / fpp(i)
        t0 = ftpp(i) - ftp(i)*fppp(i)/fpp(i)    ! q-prime x psi-prime
        t1 = sqgbsq / (t0*t0)
        t2 = cj(i) * fppp(i) - ci(i) * ftpp(i)
        t3 = t2 + cjp(i) * fpp(i) - cip(i) * ftp(i)
!        dsamp = - t0 / (sqgbsq * fpp(i) * fpp(i))
        dpamp = 2.0 * dtdp / (fpp(i) * fpp(i) * sqgbsq)
        dnamp = dpamp * pp(i)
        dsamp = t0 / (fpp(i) * fpp(i))          ! q-prime / psi-prime
        dsampn = pp(i) * dsamp / sqgbsq
        dsampk = 2.0 * dtdp * dsamp / sqgbsq
        cpamp = 2.*dtdp
        csamp = 4.*dtdp*t0/sqgbsq
        cqamp = 2.*dtdp*t0*t0/sqgbsq
!00        cqamp = 2.*dtdp*t0
        cuamp = 2.*dtdp/sqgbsq
        cwamp = cuamp * pp(i)
C
C..INITIALISATION
        hq = 0.
        hqni = 0.
        dper = 0.
        dsec = 0.
        hdsc = 0.
        hdsq = 0.
        dperni = 0.
        dsecni = 0.
        hdscni = 0.
        hdsqni = 0.
        do jk = 1,njk
          bgradg(jk) = 0.
          sigbxg(jk) = 0.
          baldp(jk)  = 0.
          balds(jk)  = 0.
          balcp(jk)  = 0.
          balcs(jk)  = 0.
          balcq(jk)  = 0.
          curvds(jk) = 0.
        end do
C
C..THE DERIVATIVES OF THE JACOBIAN ALONG AND ACROSS THE FIELD LINES.
C..SPECIFICALLY Jacob B.grad(Jacob) AND Jacob sigmaBxgrad(s).grad(Jacob)
        do l = 1,lmnb
          pmnjac = (mb(l) * fpp(i) - nb(l) * ftp(i)) * fbjac(l,i)
          cmnjac = (mb(l) * ci(i)  - nb(l) * cj(i))  * fbjac(l,i)
          do jk = 1,njk
            bgradg(jk) = bgradg(jk) - pmnjac * tsin(jk,l)
            sigbxg(jk) = sigbxg(jk) + cmnjac * tsin(jk,l)
!            bst(jk) = bst(jk) + mb(l)*fbs(l,i)*tcos(jk,l)
          end do
        end do
C
C..THE MERCIER COEFFICIENTS FOR FULL KO AND NONINTERACTING HOT PARTICLE MODEL
! --- AND CORRESPONDING DETERMINATION OF THE BALLOONING COEFFICIENTS
        do 50 jk = 1,njk
          hq = hq + 1./gssu(jk,i)
          hqni = hqni + sigmab(jk,i)/gssu(jk,i)
          sigtau = sigmab(jk,i) / taub(jk,i)
          sqrbsq = sqrt(bsq(jk,i))
          gspl   = sigbs(jk,i)/(ftp(i)*sigmab(jk,i)) - gstl(jk,i) /qsurf
          balns(jk) = (ci(i)*gstl(jk,i) + cj(i)*gspl)
          balcs(jk) = balns(jk) * sigmab(jk,i)
          shearl(jk) = - balcs(jk) / gssu(jk,i)            ! this is h_s
          pgradt = bjac(jk,i) * pp(i)* (1. + sigtau) + gparp(jk,i)
     &             + sigtau * gperp(jk,i)
          bdtemp = - sqgbsq * bjacs(jk,i)
     &             + sigbs(jk,i) * bgradg(jk) / bjac(jk,i) + t3
          baldp(jk) = 1./(1.+sigtau) * pgradt * (pgradt + bdtemp)
          balds(jk) = 1./(1.+sigtau) * pgradt * sigbxg(jk) / bjac(jk,i)
          baldn(jk) = bjac(jk,i)* (bjac(jk,i)*(1.+sigmab(jk,i))*pp(i)+
     &                             gparp(jk,i) + bdtemp)
          curvds(jk)= baldp(jk) + 
     &             shearl(jk)*sigbxg(jk)*pgradt/((1.+sigtau)*bjac(jk,i))
          curvdn(jk)= baldn(jk) + shearl(jk) * sigbxg(jk)
          balcp(jk) = sigmab(jk,i) * gssl(jk,i) -
     &                                   sigbs(jk,i)*sigbs(jk,i)/sqgbsq
          balnq(jk) = gssu(jk,i) / sigmab(jk,i)
!00          balcp(jk) =  sigmab(jk,i)**3 / (bjac(jk,i)*gssu(jk,i))
!00          balcs(jk) =(ci(i)*gstl(jk,i)+cj(i)*gspl)/(sigmab(jk,i)*sqrbsq)
!00          balcq(jk) = gssu(jk,i) / (sigmab(jk,i)*sigmab(jk,i)*sqrbsq)
          balnp(jk) = balcp(jk) / sigmab(jk,i)
C--------0---------0---------0---------0---------0---------0---------0-c
          dper   = dper + baldp(jk)
          dperni = dperni + bjac(jk,i) * (sigmab(jk,i)*bjac(jk,i) *pp(i)
     &            - sqgbsq * bjacs(jk,i) + t2)
!          dper = dper - bjacs(jk,i) + (pp(i)*bjac(jk,i)+t3)/bsq(jk,i)
!          dsca = cjp(i) - bst(jk) + t5/bsq(jk,i)
          dsec   = dsec   - parkur(jk,i)
          dsecni = dsecni - parjp(jk,i)
          hdsc = hdsc - parkur(jk,i) / gssu(jk,i)
          hdsq = hdsq + parkur(jk,i) * parkur(jk,i) / gssu(jk,i)
          hdscni = hdscni - parjp(jk,i) * sigmab(jk,i) / gssu(jk,i)
          hdsqni = hdsqni + parjp(jk,i) * sigmab(jk,i) / gssu(jk,i) *
     &                      parjp(jk,i)
 50     end do
        do 60 l=1,lmnb
         fcurvs(l) = 0.
         fcurvn(l) = 0.
         flshea(l) = 0.
         fbaldp(l) = 0.
         fbalds(l) = 0.
         fbaldn(l) = 0.
         fbalnp(l) = 0.
         fbalns(l) = 0.
         fbalnq(l) = 0.
         fbalcp(l) = 0.
         fbalcs(l) = 0.
 60      fbalcq(l) = 0.
        do 67 l=1,lmnb
         do 65 jk=1,njk
          fcurvs(l) = fcurvs(l) + curvds(jk)  * tcos(jk,l)
          fcurvn(l) = fcurvn(l) + curvdn(jk)  * tcos(jk,l)
          fbaldp(l) = fbaldp(l) + baldp(jk)   * tcos(jk,l)
          fbaldn(l) = fbaldn(l) + baldn(jk)   * tcos(jk,l)
          fbalcp(l) = fbalcp(l) + balcp(jk)   * tcos(jk,l)
          fbalcq(l) = fbalcq(l) +  gssu(jk,i) * tcos(jk,l)
          fbalnq(l) = fbalnq(l) + balnq(jk)   * tcos(jk,l)
!00        fbalcq(l) = fbalcq(l) + gssl(jk,3)  * tcos(jk,l)
!00        fbalcp(l) = fbalcp(l) + gssl(jk,1)  * tcos(jk,l)
          fbalnp(l) = fbalnp(l) + balnp(jk)   * tcos(jk,l)
          flshea(l) = flshea(l) + shearl(jk)  * tsin(jk,l)
!00 65     fbalcs(l) = fbalcs(l) + gssl(jk,2)  * tsin(jk,l)
          fbalcs(l) = fbalcs(l) + balcs(jk)   * tsin(jk,l)
          fbalns(l) = fbalns(l) + balns(jk)   * tsin(jk,l)
          fbalds(l) = fbalds(l) + balds(jk)   * tsin(jk,l)
 65      end do
 67     end do
        do 70 l=1,lmnb
          pmnjac = mb(l) * fpp(i) - nb(l) * ftp(i)
        flshea(l) = -cpamp * (mb(l)*fpp(i) - nb(l)*ftp(i))*flshea(l)
        fbaldp(l) = - dpamp * fbaldp(l)            ! d_p with opposite sign
        fbaldn(l) = - dnamp * fbaldn(l)            ! d_p with opposite sign
!        fbalds(l) = - dsamp * pmnjac * fpark(l,i)  ! d_s with opposite sign
!        fbaldt(l) = - dsamp * pmnjac * fparjp(l,i) ! d_s with opposite sign
        fbalds(l) = -dsampk * fbalds(l)            ! d_s with opposite sign
        fbaldt(l) = -dsampn * (mb(l) * ci(i) - nb(l) * cj(i))*fbjac(l,i)
!                                                  ! d_s with opposite sign
        fgssu (l) = cpamp * fbalcq(l)
        fbalcp(l) = cpamp * fbalcp(l)
        fbalnp(l) = cpamp * fbalnp(l)
        fbalcs(l) = csamp * fbalcs(l)
        fbalns(l) = csamp * fbalns(l)
        fbalcq(l) = cqamp * fbalcq(l)
        fbalnq(l) = cqamp * fbalnq(l)
        fcurvs(l) = cuamp * fcurvs(l)        
        fcurvn(l) = cwamp * fcurvn(l)        
 70     end do
      fbalcp(lmnb0) = 0.5 * fbalcp(lmnb0)
      fbalcs(lmnb0) = 0.
      fbalcq(lmnb0) = 0.5 * fbalcq(lmnb0) 
      fbaldp(lmnb0) = 0.5 * fbaldp(lmnb0)
      fbaldn(lmnb0) = 0.5 * fbaldn(lmnb0)
      fbalnp(lmnb0) = 0.5 * fbalnp(lmnb0)
      fbalns(lmnb0) = 0.
      fbalds(lmnb0) = 0.
      fbaldt(lmnb0) = 0.
      fbalnq(lmnb0) = 0.5 * fbalnq(lmnb0)
      fcurvs(lmnb0) = 0.5 * fcurvs(lmnb0)
      fcurvn(lmnb0) = 0.5 * fcurvn(lmnb0)
      flshea(lmnb0) = fpp(i)*ftpp(i) - ftp(i)*fppp(i)
      fgssu (lmnb0) = 0.5 * fgssu (lmnb0)
!
!...MERCIER CRITERION COEFFICIENTS
!
        hq     = dtdp * t1 * hq
        hqni   = dtdp * t1 * hqni
        dper   = dtdp * dper           / (fpp(i) * fpp(i) * sqgbsq)
        dperni = dtdp * dperni * pp(i) / (fpp(i) * fpp(i) * sqgbsq)
        dsec   = dtdp * t0 * dsec   / fpp(i)
        dsecni = dtdp * t0 * dsecni / fpp(i)
!        write(25,2111) hq, hqni, dper, dperni, dsec, dsecni
        hdsc   = dtdp * t1 * t0 / fpp(i) * hdsc
        hdscni = dtdp * t1 * t0 / fpp(i) * hdscni
        hdsq   = dtdp * sqgbsq / (fpp(i)*fpp(i)) * hdsq
        hdsqni = dtdp * sqgbsq / (fpp(i)*fpp(i)) * hdsqni
        dsca   = hq   * (dper  -dsec  +hdsq)
        dscani = hqni * (dperni-dsecni+hdsqni)
        dper   = hq   * (dper  -dsec)   + hdsc
        dperni = hqni * (dperni-dsecni) + hdscni
        hdsq   = (hq  *hdsq  -hdsc  *hdsc)
        hdsqni = (hqni*hdsqni-hdscni*hdscni)
C--sh-000719      if (abs(hdsq).lt.1.e-800)hdsq=1.e-800
C--sh+000719>>>>>>>>>>
!        if (abs(hdsq) .lt. tiny) then
!          hdsq = tiny
!        end if
         hdsq  =merge(hdsq  ,tiny,abs(hdsq)>=tiny)
         hdsqni=merge(hdsqni,tiny,abs(hdsq)>=tiny)
C--sh+000719<<<<<<<<<<
        hdsq   = 1. / hdsq
        hdsqni = 1. / hdsqni
        glm(i,2) = sqrt((dper*dper*hdsq+1.)*hdsq)
        glm(i,2) = 1. + 0.5*(dper*hdsq-glm(i,2))
        glm(i,1) = (0.5-hdsc)  *(0.5-hdsc)   - dsca
        glm(i,4) = sqrt((dperni*dperni*hdsqni+1.)*hdsqni)
        glm(i,4) = 1. + 0.5*(dperni*hdsqni-glm(i,4))
        glm(i,3) = (0.5-hdscni)*(0.5-hdscni) - dscani
C
      write (6,890) i, qsurf, glm(i,1), -glm(i,2), glm(i,3), -glm(i,4)
      write(16,890) i, qsurf, glm(i,1), -glm(i,2), glm(i,3), -glm(i,4)
!
!.. PRINT OUT BALLOONING FOURIER AMPLITUDES
!
      t1 = 0.5 * (s(i) + s(i-1)) 
!      write (25,870)t1,qsurf, pp(i),vp(i),tbalin,sqgbsq, ftp(i),pth(i)
      write (25,870)t1,qsurf, pp(i),vp(i),tbalin,sqgbsq, ftp(i),pp(i)
      write (25,99) (fbalcp(l),l=1,lmnb)
      write (25,99) (fbalcs(l),l=1,lmnb)
      write (25,99) (fbalcq(l),l=1,lmnb)
      t1 = 1. / (fpp(i) * fpp(i))
      write (25,99) (t1*fbjac(l,i),l=1,lmnb)
      write (25,99) (fbaldp(l),l=1,lmnb)
      write (25,99) (fbalds(l),l=1,lmnb)
      t1 = 0.5 * (s(i) + s(i-1)) 
      write (26,870)t1,qsurf, pp(i),vp(i),tbalin,sqgbsq, ftp(i),pp(i)
      write (26,99) (fbalnp(l),l=1,lmnb)
      write (26,99) (fbalns(l),l=1,lmnb)
      write (26,99) (fbalnq(l),l=1,lmnb)
      t1 = 1. / (fpp(i) * fpp(i))
      write (26,99) (t1*fbjac(l,i),l=1,lmnb)
      write (26,99) (fbaldn(l),l=1,lmnb)
      write (26,99) (fbaldt(l),l=1,lmnb)
c..print the normal curvature, the local shear, the parallel currrent
c..density, |sigma*grad s|^2, B^2, the integrated local shear and the Jacobian
c..fourier components for plotting purposes.
      do 80 l=1,lmnb
      pmnjac    = mb(l)*fpp(i) - nb(l)*ftp(i)
      fintls = flshea(l)
      if(l.ne.lmnb0)fintls = fintls/pmnjac
      write(66,859) mb(l),nb(l),fr(l,i),fz(l,i),fcurvs(l),flshea(l)
     >             ,fpark(l,i),fgssu(l),fbsq(l,i),fintls,fbjac(l,i)
      write(64,859) mb(l),nb(l),fr(l,i),fz(l,i),fcurvn(l),flshea(l)
     >             ,fparjp(l,i),fgssu(l),fbsq(l,i),fintls,fbjac(l,i)
 80   end do
 700  end do
 99   format (1p6e22.14)
 850  format (2i5,1p2e15.6)
 858  format(3i6)
 859  format(1x,2i6,1p9e13.6)
 860  format (12i6)
 870   format (1p8e16.8)
 890  format(i4,1p1e15.7,1p4e12.4)
 2111   format(1p6e14.6)
      end subroutine mercap
