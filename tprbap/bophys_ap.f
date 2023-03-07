C--------0---------0---------0---------0---------0---------0---------0-c
C  23.09.88          LAST MODIFICATION 24.03.06  UHS: BOPHYST DATA D
C-----------------------------------------------------------------------
C
      subroutine BOPHYS(ni,njk,lmnb,lmnb0,mm,nmin,nmax,lcurrf,lvmtpr,
     &                  lpress,rplmin,djp,dtdp,s,ftp,fpp,ftpp,fppp,ci,
     &                  cj,cip,cjp,pp,wmag,vp,cipi,cjpi,pvpi,pth,ppi,
     &                  equi,civ,cjv,vvp,mb,nb,tsin,tcos,lfrz,gttl,gtpl,
     &                  gppl,gssu,bjac,bjacs,vjac,sigbs,bp,bt,bsq,
     &                  fbjac,fgparp,fsigbs,fphv,glm,bgradg,gparp,
     &                  parpav,parpvi,sigmab,taub,gperp,ftaub,fgperp,
     &                  fparjp,fpark,parjp,parkur,curfac)
C
C-----------------------------------------------------------------------
C     RECONSTRUCTION OF EQUILIBRIUM IN BC
C-----------------------------------------------------------------------
C
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
      integer :: njk,nmin,nmax,lmnb,lmnb0,ni,mm,lcurrf,lvmtpr,lpress
      integer, dimension(0:36,nmin:nmax) :: lfrz
      real, dimension(ni)  :: ci   ,cj    ,cip    ,cjp    ,pp    ,wmag
     &                       ,vp   ,cipi  ,cjpi   ,pvpi   ,pth   ,ppi
     &                       ,equi ,civ   ,cjv    ,vvp
      real, dimension(njk,lmnb)  :: tsin  ,tcos
      real, dimension(lmnb,ni)   :: ftaub ,fpark  ,fparjp ,fgparp,fgperp
     &                             ,fphv
      real, dimension(lmnb,0:ni) :: fsigbs,fbjac
      real, dimension(njk,ni)    :: gparp ,sigmab ,taub   ,gperp
     &                             ,parjp ,parkur
      real, dimension(njk,0:ni)  :: sigbs ,bsq    ,gttl   ,gtpl  ,gppl
     &                             ,gssu  ,bjac   ,bjacs  ,vjac  ,bp ,bt
      real, dimension(njk)       :: bgradg
      real, dimension(ni)        :: parpav,parpvi
      real, dimension(0:ni)      :: s
      real, dimension(ni+1)      :: ftp   ,fpp    ,fppp   ,ftpp
      real, dimension(ni,4)      :: glm
      real                       :: curfac,b0r0   ,rplmin ,djp   ,dtdp
      integer, dimension(lmnb)   :: mb    ,nb
C
C.. Local Scalars ..
      integer I,JK,L,m,mboo,n
      real DIFJAC,dvol,EQUIN,FPPI,FTPI,PMN,T7,T8,tds,tsi,VPI,djpboo
C
C.. Intrinsic Functions ..
      intrinsic ABS, MAX
C
      data djpboo /1.0e-16/
!
C ... Executable Statements ...
C
      djpboo = max(djp,djpboo)
C
C
      if (LCURRf .ne. 2) then
 2000   format (/'************** CI/J FROM BP/T *******************')
        write (16,2000)
      end if
      if (LCURRf .eq. 2) then
 2001   format (/'************** CI/J FROM VMEC *******************')
        write (16,2001)
      end if
      do 115 I = 1,NI
C
        CI(I) = 0.
        CJ(I) = 0.
        VP(I) = 0.
        pvpi(i) = 0.
        WMAG(I) = 0.
        PARPAV(i) = 0.
        dvol = (s(i)-s(i-1)) * dtdp
C
!..BP AND BT CORRESPOND TO THE ACTUAL MAGNETIC FIELDS IN THE COVARIANT
!..REPRESENTATION. THE CURRENT FLUX FUNCTIONS CI AND CJ CORRESPOND TO SIGMA
!..TIMES THE MAGNETIC FIELD COMPONENTS.
        do JK = 1,NJK
C
          BP(JK,I) = (GPPL(JK,I)*FTP(I)+GTPL(JK,I)*FPP(I))
C
          BT(JK,I) = (GTPL(JK,I)*FTP(I)+GTTL(JK,I)*FPP(I))
C
C     FBSQ(L,I) AND BSQ(JK,I) IS COMPUTED DIRECTLY IN VMTOBO (LIKE R/Z)
C
C          BSQ(JK,I) = (  GPPL(JK,I) * FTP(I)**2
C    $                +2.*GTPL(JK,I) * FPP(I) * FTP(I)
C    $                  + GTTL(JK,I) * FPP(I)**2  ) / BJAC(JK,I)
C
          if (LCURRf .ne. 2) then
            CI(I) = CI(I) - DTDP*BP(JK,I) * SIGMAB(JK,I)
            CJ(I) = CJ(I) + DTDP*BT(JK,I) * SIGMAB(JK,I)
          end if
          if (LCURRf .eq. 2) then
            CI(I) = CIV(I)
            CJ(I) = CJV(I)
          end if
          VP(I) = VP(I) - DTDP*BJAC(JK,I)
          pvpi(i) = pvpi(i) - dtdp*vjac(jk,i)
C
        end do
C
!...MAGNETIC FIELD ENERGY AND FLUX SURFACE AVERAGE OF PARALLEL
!...ENERGETIC SPECIES PRESSURE
        do JK = 1,NJK
          WMAG(I) = WMAG(I) + 0.5*BSQ(JK,I)*BJAC(JK,I)*DVOL
          PARPAV(I) = PARPAV(I) + DTDP  * GPARP(JK,I)
C
        end do
C
 115  end do
 1033 format (
     &'DIFFERENTIAL VOLUMES FROM GEOMETRIC AND MAGNETIC BOOZER J',
     &'ACOBIANS'
     &)
C
C      write (6,1033)
      write (16,1033)
      do i = 1,ni
 1034   format (23x,1p1e13.6,5x,1p1e13.6)
C      write (6,1034) pvpi(i),vp(i)
        write (16,1034) pvpi(i), vp(i)
      end do
C
      mboo = 0
      do m = 0,mm
        do n = nmin,nmax
          mboo = mboo + lfrz(m,n)
        end do
      end do
 1036 format (
     &/
     &'Modes in the Boozer Table in which R,Z,Phi Fourier Amplit',
     &'ude exceeds:'
     &,/25x,1p1e13.5,3x,i4,/)
      write (16,1036) rplmin, mboo
      do n = nmin,nmax
 1037   format (5x,37i2,i3)
        write (16,1037) (lfrz(m,n), m = 0,36), n
      end do
C
      do I = 1,NI-1
C
        tsi = 2.0 / (s(i+1)-s(i-1))
        VPI = 0.5 * (VP(I+1)+VP(I))
        FTPI = 0.5 * (FTP(I+1)+FTP(I))
        FPPI = 0.5 * (FPP(I+1)+FPP(I))
        CIPI(I) = TSI * (CI(I+1)-CI(I))
        CJPI(I) = TSI * (CJ(I+1)-CJ(I))
        PVPI(I) = TSI * (PTH(I+1)-PTH(I))
        PARPVI(I) = 0.5 * (PARPAV(I+1) + PARPAV(I))
        EQUIN = ABS(CJPI(I)*FPPI) + ABS(CIPI(I)*FTPI) + ABS(VPI*PVPI(I))
     &        + ABS(PARPVI(I))
        EQUI(I) = ((CJPI(I)*FPPI-CIPI(I)*FTPI)+PARPVI(I)-VPI*PVPI(I))
     &          / EQUIN
        PPI(I) = (CJPI(I)*FPPI-CIPI(I)*FTPI+PARPVI(I)) / VPI
      end do
      pvpi(ni) = 0.
 1002 format (
     &//1x,'***** BOPHYS, RECONSTRUCTED EQUILIBRIUM *****',//3x,
     &'I     VP(I)    PPI(I)     PVPI(I)   PARPVI(I)   CJPI(I)    CIPI',
     &'(I)     EQUI'
     &,/)
C
      write (6,1002)
C     WRITE(16,1002)
      do I = 1,NI
C-----------------------------------------------------------------------
 1003   format (1x,i3,1p7e11.4)
        write (6,1003)
     &        I, VP(I), PPI(I), PVPI(I),-PARPVI(I)/VPI, CJPI(I),CIPI(I),
     &         EQUI(I)
C        WRITE(16,1003)I,VP(I),PPI(I),PVPI(I),CJPI(I),CIPI(I),EQUI(I)
      end do
C
C     INTER/EXTRAPOLATION TO HALF GRID
C
      if (lpress .eq. 2) then
        do i = 1,ni-1
          ppi(i) = pvpi(i)
        end do
      end if
C
      do 145 I = 2,NI-1
C
        tds = 2.0 / (s(i+1)+s(i)-s(i-1)-s(i-2))
        FPPP(I) = TDS * (FPP(I+1)-FPP(I-1))
        FTPP(I) = TDS * (FTP(I+1)-FTP(I-1))
        PP(I) = 0.5 * (PPI(I)+PPI(I-1))
        CIP(I) = 0.5 * (CIPI(I)+CIPI(I-1))
        CJP(I) = 0.5 * (CJPI(I)+CJPI(I-1))
C
!...NOTE ARRAY BJACS CONTAINS delta(ln(JACOBIAN))/delta(s)
        do JK = 1,NJK
          BJACS(JK,I) = TDS * (BJAC(JK,I+1)-BJAC(JK,I-1)) / BJAC(JK,I)
        end do
C
 145  end do
      do 200 I = 2,NI-1
C
        do JK = 1,NJK
          VJAC(JK,I) = (FPP(I)*CJ(I)-FTP(I)*CI(I)) / BSQ(JK,I)
        end do
C
!...COMPUTE FOURIER AMPLITUDES OF JACOBIAN (for comparison purposes),
!...sigma B_s AND PART OF THE A7 COEFFICIENT.
        do 185 L = 1,LMNB
C
Cbjac         FBJAC(L,I) = 0.
          FPHV(L,I)   = 0.
          FSIGBS(L,I)   = 0.
C
C
          do JK = 1,NJK
Cbjac         FBJAC(L,I) = FBJAC(L,I)
Cbjac     $              + BJAC(JK,I) *TCOS(JK,L)
            FPHV(L,I) = FPHV(L,I) + VJAC(JK,I)*TCOS(JK,L)
          end do
C
Cbjac      FBJAC(L,I) = 2. * DTDP * FBJAC(L,I)
          FPHV(L,I) = 2. * DTDP * FPHV(L,I)
            FPHV(lmnb0,i) = 0.5 * FPHV(lmnb0,i)

!          if (mb(l).eq.0 .and. nb(l).eq.0) then
!            fphv(l,i) = 0.5 * fphv(l,i)
!          end if
C-----------------------------------------------------------------------
C
C     COMPARE (DISCRETIZED) BJAC COEFFICIENTS WITH PHV(=BJAC) COMPUTED
C     BY DIRECT INTEGRATION OVER BSQ
C
          DIFJAC = 2. * (FBJAC(L,I)-FPHV(L,I)) / (VP(I)+VVP(I))
C
          if (ABS(DIFJAC).ge.1.e-06 .and. LVMTPR.eq.1) then
 1006       format (1x,3i4,1p2e20.8,1p3e10.2)
            write (16,1006)
     &            I, MB(L), NB(L), VP(I), VVP(I), FBJAC(L,I), FPHV(L,I),
     &            DIFJAC
          end if
C
C..DETERMINATION OF THE FOURIER AMPLITUDES OF (sigma B_s) STORED IN
C..fsigbs USING THE RADIAL FORCE BALANCE RELATION.
          if (MB(L).ne.0 .or. NB(L).ne.0) then
            PMN         = MB(L) * FPP(I) - NB(L) * FTP(I)
            PMN = PMN / (PMN*PMN + DJPBOO * (PMN+FPP(I)) * (PMN+FPP(I)))
*******
            FSIGBS(L,I) = (PP(I) * FBJAC(L,I) + FGPARP(L,I) ) * PMN
!     &                   * PMN / (PMN*PMN+DJP)
*******
            FPARK(L,I)  = (MB(L)*CI(I)-NB(L)*CJ(I))*FSIGBS(L,I) /
     &                     (CJ(I) * FPP(I) - CI(I) * FTP(I))
            FPARJP(L,I) = (MB(L)*CI(I)-NB(L)*CJ(I))* PP(I) *FBJAC(L,I) /
     &                   (CJ(I)*FPP(I)-CI(I)*FTP(I)) * PMN
!    &                   /(PMN*PMN+DJP)
          end if
 185    end do
           FPARK(LMNB0,I) = (CJ(I) * CIP(I) - CI(I) * CJP(I)) /
     &                    (CJ(I) * FPP(I) - CI(I) * FTP(I))
           FPARJP(LMNB0,I) = curfac * FPARK(LMNB0,I)
C
C
!..GSSU CORRESPONDS TO (sigma)^2 TIMES GRAD(s).GRAD(s)
        T7 = (CJ(I)*FPP(I)-CI(I)*FTP(I)) / (FTP(I)*FTP(I))
        T8 = CJ(I) * CJ(I) / (FTP(I)*FTP(I))
        do JK = 1,NJK
          GSSU(JK,I) = T7*SIGMAB(JK,I)*GTTL(JK,I) - T8
          SIGBS(JK,I) = 0.
          TAUB(JK,I) = 0.
          PARKUR(JK,I) = 0.
          PARJP (JK,I) = 0.
        end do
C
!..sigbs CALCULATED HERE CORRESPONDS TO (sigma B_s).
        do L = 1,LMNB
          do JK = 1,NJK
            SIGBS(JK,I) = SIGBS(JK,I) + FSIGBS(L,I)     * TSIN(JK,L)
            TAUB(JK,I)  = TAUB(JK,I)  + FTAUB(L,I)      * TCOS(JK,L)
            PARKUR(JK,I)= PARKUR(JK,I)+ FPARK(L,I)      * TCOS(JK,L)
            PARJP(JK,I) = PARJP(JK,I) + FPARJP(L,I)     * TCOS(JK,L)
          end do
        end do
!
C      WRITE (16,1004)
C1004 FORMAT(/1X,'FBS(1...4),CJP,CIP'/)
C      WRITE (16,1005) I,(FBS(L,I),L=1,4),CJP(I),CIP(I)
C1005 FORMAT(I4,1P4E11.4,1X,1P2E11.4)
C
C     CI,FPP AND VP MULTIPLIED BY NSTA:
C
CPERF      CI(I) = CI(I) * NSTA
CPERF     CIP(I) =CIP(I) * NSTA
CPERF      VP(I) = VP(I) * NSTA
CPERF     FPP(I) = FPP(I) * NSTA
CPERF    FPPP(I) =FPPP(I) * NSTA
C
 200  end do

      T7 = (CJ(1)*FPP(1)-CI(1)*FTP(1)) / (FTP(1)*FTP(1))
      T8 = CJ(1) * CJ(1) / (FTP(1)*FTP(1))
      do JK = 1,NJK
         GSSU(JK,1) = T7*SIGMAB(JK,1)*GTTL(JK,1) - T8
      enddo
      T7 = (CJ(ni)*FPP(ni)-CI(ni)*FTP(ni)) / (FTP(ni)*FTP(ni))
      T8 = CJ(ni) * CJ(ni) / (FTP(ni)*FTP(ni))
      do JK = 1,NJK
         GSSU(JK,ni) = T7*SIGMAB(JK,ni)*GTTL(JK,ni) - T8
      enddo
C
!      if (lpress .ne. -9) then
!        call mercie(ni,njk,lmnb,dtdp,ftp,fpp,ftpp,fppp,ci,cj,pp,
!     &              cjp,cip,mb,nb,tsin,gssu,fbjac,bjac,bjacs,parkur,
!     &              parjp,sigmab,taub,sigbs,gparp,gperp,glm,bgradg)
!      end if
C
C     CALL RADPLO(NI,   CI,'   CI$',1,0)
C     CALL RADPLO(NI,  CIP,'  CIP$',1,0)
C     CALL RADPLO(NI,   CJ,'   CJ$',1,0)
C     CALL RADPLO(NI,  CJP,'  CJP$',1,0)
C     CALL RADPLO(NI,  FPP,'  FPP$',1,0)
C     CALL RADPLO(NI, FPPP,' FPPP$',1,0)
C     CALL RADPLO(NI,  FBS,'   VP$',1,0)
      end subroutine bophys
