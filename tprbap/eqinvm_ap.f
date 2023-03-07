CDIRD NOLIST
C-----------------------------------------------------------------------
C /20.09/89/ UHS        LAST MODIFIED /22.03.06/ UHS: EQINVMF DATA D
C-----------------------------------------------------------------------
C
      subroutine EQINVM
C
      use tprb_param
      use physeq
      use radfct
      use triggs
      use paprpl
      use vmcfou
C
C.. Implicits ..
      implicit none
C
C.. Local Scalars ..
      integer i,IFSQ,insol,J,K,KPRES,L,lc,LCL,LCX,lr,M0,MN,MN0,MPNT,
     &        MPOL,MRZ,MSOLOV,mtbl,mxdif,N0,NIT,NPERIODS,NRHO,
     &        NRZ,NSI,NSIN,ntbl,NTHSOL,NTOR,NTOR0,nxdif,nzero,m,n
     &       ,NBALMX ,nowall   ,npwall
      character zdat1*24,zdat3*24
      real AIOTAX,AR,ASP,BTOR0,DVDSJ,ELONG,
     &     ENM1,EP2,FLAMP,GMN,PHIPAM,PRAX,PSIB,qonax,
     &     RHO,rhohaf,RSOL,signz,THETA,TWO,VOLI,xlc,XM,XN,XRC,XZC
     &    ,SMNC,TMNC,PBMNC,PPMNC
     &    ,awall ,ewall ,dwall ,gwall ,drwal ,dzwal ,wct
C
C.. External Calls ..
      external fdate
C
C.. Intrinsic Functions ..
      intrinsic ATAN, COS, NINT, REAL, abs, sqrt, max
!      include 'tprcom.inc'
C
C.. Local Arrays ..
      integer, parameter ::  IP_193=193
      integer, parameter ::  IP_36 =36
      integer, parameter ::  IP_28 =28
!      real, allocatable  :: AIOTAH(:)  ,AMASS(:)  ,PRES(:)  ,PHIP1(:)  ,
!     &                      VOL(:)     ,VPMOM(:)  ,PPMNF(:,:,:)        ,
!     &                      FLI(:,:,:) ,RF(:,:,:) ,ZF(:,:,:),VJB(:,:,:),
!     &                      SMNF(:,:,:),TMNF(:,:,:),PBMNF(:,:,:)
      real, dimension(IP_193) :: AIOTAH ,AMASS ,PRES ,PHIP1,VOL  ,VPMOM
!
      real, dimension(0:IP_36,-IP_28:IP_28,IP_193) :: FLI  ,RF   ,ZF
     &                                               ,VJB  ,SMNF ,TMNF
     &                                               ,PBMNF,PPMNF
C
C ... Executable Statements ...
!      allocate
!     & (AIOTAH(IP_193) ,RF(0:IP_36,-IP_28:IP_28,IP_193)   ,
!     &  AMASS(IP_193)  ,ZF(0:IP_36,-IP_28:IP_28,IP_193)   ,
!     &  PRES(IP_193)   ,FLI(0:IP_36,-IP_28:IP_28,IP_193)  , 
!     &  PHIP1(IP_193)  ,SMNF(0:IP_36,-IP_28:IP_28,IP_193) , 
!     &  VOL(IP_193)    ,TMNF(0:IP_36,-IP_28:IP_28,IP_193) , 
!     &  VPMOM(IP_193)  ,PBMNF(0:IP_36,-IP_28:IP_28,IP_193), 
!     &                  PPMNF(0:IP_36,-IP_28:IP_28,IP_193))
      MPNT = 0
      MSOLOV = 0
      NRHO = 0
      NTHSOL = 0
      zdat1 = ' '
      AR = 0.0
      ASP = 0.0
      BTOR0 = 0.0
      ELONG = 0.0
      ENM1 = 0.0
      EP2 = 0.0
      signz = 1.0
      do i=1,IP_193
         AIOTAH(i) = 0.0
         AMASS(i) = 0.0
         PHIP1(i) = 0.0
         PRES(i) = 0.0
         AMASS(i) = 0.0
         VOL(i) = 0.0
         VPMOM(i) = 0.0
      end do
      do k=1,IP_193
         do j=-IP_28,IP_28
            do i=0,IP_36
               FLI(i,j,k)   = 0.0
               RF(i,j,k)    = 0.0
               VJB(i,j,k)   = 0.0
               ZF(i,j,k)    = 0.0
               SMNF(i,j,k)  = 0.0
               TMNF(i,j,k)  = 0.0
               PBMNF(i,j,k) = 0.0
               PPMNF(i,j,k) = 0.0
            end do
         end do
      end do
C
C
C
C-----------------------------------------------------------------------
C     ---- INPUT FORMATTED -----
C     READ TABLE OF BC - COEFF AND PRESET MB AND NB
C     READ VMEC OUTPUT FROM 18 AND PRESET MX,NX,ML AND NL
C     ---- VMEC89 ----- /03.11.89/ -------
C-----------------------------------------------------------------------
C
C     initialise parameters
      nzero = 0
      dth = 1. / nj
      dph = 1. / nk
      dtdp = dth * dph
C     PI   = 3.1415926535898  D0
      PI = 4.0 * ATAN(1.0)
      TPI = 2 * PI
      alpha = 0.
      djp = 0.
      lnbaln = 0
C
      insol = 1
      mm = 3
      nmin = 0
      nmax = 0
      nprocs = 1
 1000 format (15x,a8,//)
C      READ (5,1000) zdat(3)
      read (5,1000) zdat3
 1008 format (5x,8i6,/////)
      read (5,1008) mm, nmin, nmax, mms, nsmin, nsmax, nprocs, insol
C     INPUT FOR SOLOV'EV EQUILIBRIA
C--sh-000831      QONAX  = 0.35
      if (INSOL .gt. 0) then
        if (insol .eq. 1) then
 1010     format (
     &    5x,'Solov ev Equilibria with radial variable as radius')
          write (6,1010)
        end if
        if (insol .eq. 2) then
 1011     format (
     &    5x,'Solov ev Equilibria with radial variable as volume')
          write (6,1011)
        end if
        NRHO = NI + 1
        ENM1 = 1.0
        MSOLOV = MLMNV - 1
        NTHSOL = 200
        AR = 1.0
        BTOR0 = 1.0
        ASP = 3.
        ELONG = 1.0
        GAM = 0.
        wct = 0.0001
      end if
C
      rewind (18)
      if (INSOL .eq. 0) then
C
C--sh-000831      IF (LIOTPL.EQ.1) THEN
C--sh-000831      DO 37 J=1,NSIN
C--sh-000831 37   AIOTAH(J)=AIOTA(J)/NPER
Cplot       CALL RADPLO (NI,AIOTAH,        'IOTA$', 1,0)
C--sh-000831      ENDIF
C
 40     format (1x,1pe22.12,1x,1p3e12.5,8i4,i1)
        read (18,40)
     &       VOLI, GAM, ENM1, alpha, MPNT, NRHO, MPOL, NTOR, NTOR0, MN0,
     &       NIT, IFSQ, KPRES
      end if
!       WRITE (16,40)VOLI,GAM,ENM1,alpha  ,MPNT,NRHO,MPOL,
!     &  NTOR,NTOR0,MN0,NIT,IFSQ,KPRES
C ....   ALPHA .NE. 0  CORRESPONDS ONLY TO HELICAL SYMMETRY.
      NPERIODS = NINT(1/ENM1)
      NPER = NPERIODS
      NSIN = NRHO - 1
      NSI = NSIN
C       call fdate(zdat(1))
CNEC-SX4      call fdate(zdat1)
      call fdate(zdat1)
 1001 format (
     &'1'//1x,
     &'###### TERPSICHORE #########  EQINVM  ######### TERPSICHORE ##',
     &'## '
     &,/8x,a24,3x,a24,/,//4x,'NIM',3x,'IVAC',3x,'NJ',4x,'NK',4x,'MM',2x,
     &'NMIN',2x,'NMAX',//1x,i6,7i6,//4x,'TABLE OF R AND Z COEFFICIENTS',
     &/)
C-CRAY      call date(zdat1)
C      call clock(zdat(2))
C      CALL DATJOB (ZDAT)
C      WRITE (16,1001) (ZDAT(K),K=1,3),NSI,ivac,NJ,NK,MM,NMIN,NMAX
C      WRITE ( 6,1001) (ZDAT(K),K=1,3),NSI,ivac,NJ,NK,MM,NMIN,NMAX
      write (16,1001) ZDAT1, ZDAT3, NSI, ivac, NJ, NK, MM, NMIN, NMAX
      write (6,1001) ZDAT1, ZDAT3, NSI, ivac, NJ, NK, MM, NMIN, NMAX
 1021 format (//4x,'TABLE OF R/Z (PR/Z) COEFFICIENTS',/)
      write (16,1021)
      write (6,1021)
 1002 format (3x,'M=  0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 ...N',/)
C
C     READ MATRIX, CONTAINING THE WANTED MODE NUMBERS FOR R AND Z
C     IF (LFRZ(M,N) = 1, LMNB   = LMNB  +1, M/N(LMNB) = M/N
C
      write (16,1002)
      write (6,1002)
C
      do N = NMIN,NMAX
C
 1003   format (5x,37i2,i3)
C
        read (5,1003) (LFRZ(M,N), M = 0,36)
 1004   format (6x,37i2,i3)
        write (16,1004) (LFRZ(M,N), M = 0,36), N
        write (6,1004) (LFRZ(M,N), M = 0,36), N
      end do
 1006 format (//,4i12,//,4i12,//,4i12,//,4i12,//,4i12)
C
C     READ PRINT AND PLOT PARAMETERS
C
C
C     PAPRPL: PRINT & PLOT SWITCHES
C
      read (5,1006)
     &     LLAMPR, LVMTPR, LMETPR, LRIPPL, LRADMN, LRADMX, LNBALN,
     &     LNBALX, LXYZPR, LIOTPL, LSKEWN, LMXBAL, LCURRF, LMESHP,
     &     LMESHV, LITERS, LXYZPL, LEFPLS, LBOOTJ, LPRESS
      write (16,1006)
     &      LLAMPR, LVMTPR, LMETPR, LRIPPL, LRADMN, LRADMX, LNBALN,
     &      LNBALX, LXYZPR, LIOTPL, LSKEWN, LMXBAL, LCURRF, LMESHP,
     &      LMESHV, LITERS, LXYZPL, LEFPLS, LBOOTJ, LPRESS
 1007 format (//,1p6e12.4,5x,i2)
C
      read (5,1007) pvac, parfac, qonax, qn, dsvac, qvac, nowall
      read (5,1007) awall, ewall, dwall, gwall, drwal, dzwal, npwall
 1015 format (//,1p5e12.4)
      read (5,1015) rplmin, xplo, djp, wct, curfac
      AIOTAX = 1.0 / QONAX
      if (ivac.ne.0.and.pvac.le.1.0 .or. qvac.le.1.0) stop
 1016 format (
     &/
     &'   QN   QONAX  AIOTAX  PARFAC SVMAX - 1    PVAC      QVAC',
     &'      NOWALL'
     &,/1x,1f4.1,3f8.4,1p3e10.3,i8/)
C
      write (6,1016)
     &      QN, QONAX, aiotax, parfac, dsvac, pvac, qvac, nowall
      write (16,1016)
     &      QN, QONAX, aiotax, parfac, dsvac, pvac, qvac, nowall
 1018 format (
     &/
     &' npwall    awall     ewall     dwall     gwall     drwal ',
     &'    dzwal'
     &,/1x,i4,3x,1p6e10.3,/)
      write (6,1018) npwall, awall, ewall, dwall, gwall, drwal, dzwal
      write (16,1018) npwall, awall, ewall, dwall, gwall, drwal, dzwal
 1019 format ('rplmin =',1p1e12.4,' xplo =',1p1e12.4,' djp =',1p1e12.4,
     &' wct =',1p1e12.4,' curfac=',1p1e12.4/)
      write (6,1019) rplmin, xplo, djp, wct ,curfac
      write (16,1019) rplmin, xplo, djp, wct ,curfac
C
      LMNB = 0
      LMNL = 0
C
      do M = 0,MM
        do N = NMIN,NMAX
C
C        WRITE(6,*) M,N,LFRZ(M,N)
C       WRITE(16,*) M,N,LFRZ(M,N)
          if (LFRZ(M,N) .ne. 0) then
            LMNB = LMNB + 1
C        WRITE(6,*)M,N,LMNB,LFRZ(M,N)
C       WRITE(16,*)M,N,LMNB,LFRZ(M,N)
            MB(LMNB) = M
            NB(LMNB) = N * NPER
C
C     COMPUTE LAMBDA WITH BC SET
            if (M.ne.0 .or. N.ne.0) then
              LMNL = LMNL + 1
              ML(LMNL) = M
              NL(LMNL) = N * NPER
            end if
          end if
C     WRITE (16,*)LMNB,LMNL
        end do
      end do
C
C     WRITE (16,*)LMNB,LMNL
C
      M0 = 0
      N0 = 0
      LCX = 0
      LCL = 0
C
      if (INSOL .eq. 0) then
        xlc = 0. !
        do J = 1,NRHO
C--sh-000831        LJ=(J-1)*MPNT
          do MN = 1,MPNT
C--sh-000831 50     FORMAT(1X,1P5E14.6,I3)
C--sh-000831 51     FORMAT(1X,1P2E10.2,1P2E22.14)
 52         format (1x,1p5e22.14)
 53         format (1x,1p4e22.14)
! 52         format (1x,1p5e14.6)
            read (18,52) XM, XN, XRC, XZC, GMN
            read (18,53) SMNC, TMNC, PBMNC, PPMNC
C       IF ( J.LE.2.OR.J.GE.NRHO-2)
C    $  WRITE(16,50)  XM, XN, XRC, XZC, GMN*VFVM, J
C    $  WRITE( 6,50)  XRC, XZC, XLC, XM, XN, J
C--sh-000831        NXM = XN/ NPER
            NRZ = nINT(XN/NPER)
            MRZ = nINT(XM)
C
C     PRESET INDIRECT INDEX ARRAYS FOR EQ. COOR.
C
            if (J .eq. 2) then
C
              LCX = LCX + 1
              MX(LCX) = MRZ
              NX(LCX) = NRZ * NPER
C
              if (MRZ.ne.0 .or. NRZ.ne.0) then
C
                LCL = LCL + 1
              end if
            end if
C
            RF(MRZ,NRZ,J) = XRC
            ZF(MRZ,NRZ,J) = XZC
            FLI(MRZ,NRZ,J) = XLC
            VJB(MRZ,NRZ,J) = GMN
            SMNF(MRZ,NRZ,J) = SMNC
            TMNF(MRZ,NRZ,J) = TMNC
            PBMNF(MRZ,NRZ,J) = PBMNC
            PPMNF(MRZ,NRZ,J) = PPMNC
            if (M0 .lt. MRZ) then
              M0 = MRZ
            end if
            if (N0 .lt. ABS(NRZ)) then
              N0 = ABS(NRZ)
            end if
          end do
C--sh-000831        LJ=(J-1)*MPNT
        end do
        signz = zf(1,0,nrho) / abs(zf(1,0,nrho))
C
        LMNV = MPNT
CBC     LMNL = MPNT - 1
C
        NBALMN = MAX(1,LRADMN)
        NBALMX = MAX(NBALMN+2,LRADMX)
        if (nbalmx .ge. nrho .or. nbalmx-nbalmn .ge. ni) then
           print*,'ni,nrho,lradmx',ni,nrho,lradmx
           stop
        endif
!
        EPS = 1. / RF(0,0,NRHO)
        EL = TPI / EPS * ENM1
        do J = nbalmn,nbalmx+1
          do M = 0,MM
            do N = NMIN,NMAX
              do L = 1,LMNV
                if (MX(L).eq.M .and. NX(L).eq.N*NPER) then
                  FRV(L,J-nbalmn) = RF(M,N,J)
                  FZV(L,J-nbalmn) = signz * ZF(M,N,J)
                  FVJAC(L,J-nbalmn) = signz * VJB(M,N,J)
                  FSIGMV(L,J-nbalmn) = SMNF(M,N,J)
                  FTAUV(L,J-nbalmn)  = TMNF(M,N,J)
                  FPARPV(L,J-nbalmn) = PBMNF(M,N,J)
                  FPERPV(L,J-nbalmn) = PPMNF(M,N,J)
                end if
                if (MX(L).gt.0 .and. J.eq.1) then
                  FRV(L,0) = 0.
                  FZV(L,0) = 0.
                end if
              end do
            end do
          end do
        end do
C
        do J = nbalmn,nbalmx+1
          do M = 0,MM
            do N = NMIN,NMAX
              do L = 1,LMNV
                if (ML(L).eq.M .and. NL(L).eq.N*NPER) then
                  FVLI(L,J-nbalmn) = FLI(M,N,J)
                end if
              end do
            end do
          end do
        end do
      else
C
        LCX = 0
        do M = 0,MSOLOV
          LCX = LCX + 1
          MX(LCX) = M
          NX(LCX) = 0
        end do
        LMNV = LCX
C
C      FOURIER ANALYSIS OF R AND Z ON INTEGER GRID POINTS (SOLOV'EV).
C
        EPS = 1. / ASP
        EP2 = 2. * EPS
        EL = TPI / EPS * ENM1
        do M = 0,MSOLOV
          FRV(M+1,0) = 0.
          FZV(M+1,0) = 0.
        end do
        FRV(1,0) = AR * ASP
C
        do J = 1,NSI
          RHO = EP2 * J / NSI
          if (insol .eq. 2) then
            rho = ep2 * sqrt(REAL(j)/REAL(nsi))
          end if
          do M = 0,MSOLOV
            FRV(M+1,J) = 0.
            TWO = 2.0 * FRV(1,0) / NTHSOL
            if (M .eq. 0) then
              TWO = 0.5 * TWO
            end if
            do K = 1,NTHSOL
              THETA = TPI * (K-1) / NTHSOL
              FRV(M+1,J) =
     &          FRV(M+1,J) + SQRT(1.0+RHO*COS(THETA))*COS(MX(M+1)*THETA)
            end do
            FRV(M+1,J) = TWO * FRV(M+1,J)
            FZV(M+1,J) = MX(M+1) * ELONG * FRV(M+1,J)
          end do
        end do
C
C      FOURIER AMPLITUDES OF LAMBDA ON HALF-INTEGER GRID (SOLOV'EV).
C
        FLAMP = -ELONG*AR*AR*BTOR0/nthsol
        if (insol .eq. 1) then
          flamp = 2. * flamp
        end if
        do J = 1,NSI
          RHO = EP2 * (J-0.5) / NSI
          if (insol .eq. 2) then
            rho = ep2 * sqrt((j-0.5)/REAL(nsi))
          end if
          FVJAC(1,J) = 0.
          do M = 1,LMNL
            FVJAC(M+1,J) = 0.
            FVLI(M,J) = 0.
            two = flamp / ml(m)
            if (insol .eq. 1) then
              TWO = two * (J-0.5) / REAL(nsi)
            end if
            do K = 1,NTHSOL
              THETA = TPI * (K-1) / NTHSOL
              RSOL = 1.0 / SQRT(1.0+RHO*COS(THETA))
              FVJAC(M+1,J) = FVJAC(M+1,J) + RSOL*COS(ML(M)*THETA)
              FVLI(M,J) = FVLI(M,J) + RSOL*RSOL*RSOL*COS(ML(M)*THETA)
            end do
            FVJAC(M+1,J) = ML(M) * FRV(1,0) * TWO * FVJAC(M+1,J) / BTOR0
            FVLI(M,J) = TWO * FVLI(M,J)
          end do
        end do
      end if
C
C--sh-000831       DO 29 I = 0,NI
C--sh-000831       DO 29 L = 1,LMNV
C
C         WRITE (16,1005) I,L,MX(L),NX(L), FRV (L,I), FZV(L,I)
C         WRITE ( 6,1005) I,L,MX(L),NX(L), FRV (L,I), FZV(L,I)
C--sh-000831   29     CONTINUE
C--sh-000831 1005 FORMAT (1X,'EQIN', 4I3,1P2E20.10)
C
C--sh-000831        IF ( N0 .LE. 1 )  N0 = 2
C--sh-000831 1111   FORMAT(4I3,1P5E12.4)
      if (INSOL .eq. 0) then
        read (18,52)
     &       (AIOTAH(J), AMASS(J), PRES(J), PHIP1(J), VOL(J), J = 1,NRHO
     &       )
      end if
C--sh-000831 60     FORMAT(  ' I, FTP, FPP',(I3,1P3E15.6))
 70   format (
     &/,'      IOTA           M           VP          P ',/,
     &(1p2e10.2,1pe20.8,1pe10.2))
C
      write (16,70) (AIOTAH(J), AMASS(J), VOL(J), PRES(J), J = 1,NRHO)
      write (6,70) (AIOTAH(J), AMASS(J), VOL(J), PRES(J), J = 1,NRHO)
C
C--sh-000831        GAM1 = GAM - 1
C
      if (INSOL .eq. 0) then
        s(0) = (nbalmn-1)/real(NSIN)
!
!      SURFACE FUNCTIONS FOR SOLOV'EV EQUILIBRIA.
        do J = nbalmn+1,nbalmx+1
          s(j-NBALMN) = (j-1) / REAL(nsin)
          FTP(J-NBALMN) = PHIP1(J)
          FPP(J-NBALMN) = AIOTAH(J) * FTP(J-NBALMN)
          AIOTA(J-NBALMN) = FPP(J-NBALMN) / FTP(J-NBALMN)
          PTH(J-NBALMN) = PRES(J)
          VPMOM(J-NBALMN) = VOL(J)
          DVDSJ = VPMOM(J-NBALMN)
          AM(J-NBALMN) = PTH(J-NBALMN) * DVDSJ
          if (GAM .ne. 0) then
            AM(J-NBALMN) = (PTH(J-NBALMN))**(1/GAM) * DVDSJ
          end if
C          WRITE(16,60) J-NBALMN, FTP(J-NBALMN), FPP( J-NBALMN)
        end do
      else
        PHIPAM = -ELONG*AR*AR*BTOR0/NTHSOL
        if (insol .eq. 2) then
          phipam = 0.5 * phipam
        end if
        s(0) = 0.
        do J = 1,NSI
          VPMOM(J) = 0.
          FTP(J) = 0.
          s(j) = j / REAL(nsi)
          rhohaf = (j-0.5) / REAL(nsi)
          RHO = EP2 * rhohaf
          if (insol .eq. 2) then
            rho = ep2 * sqrt(rhohaf)
          end if
          do K = 1,NTHSOL
            THETA = TPI * (K-1) / NTHSOL
            RSOL = 1. / SQRT(1.0+RHO*COS(THETA))
            VPMOM(J) = VPMOM(J) + RSOL
            FTP(J) = FTP(J) + RSOL*RSOL*RSOL
          end do
          vpmom(j) = -phipam*frv(1,0)*vpmom(j)/btor0
          if (insol .eq. 1) then
            vpmom(j) = vpmom(j) * rhohaf
          end if
          ftp(j) = phipam * ftp(j)
          if (insol .eq. 1) then
            ftp(j) = ftp(j) * rhohaf
          end if
          FVJAC(1,J) = -VPMOM(J)
        end do
C
        PSIB = -0.5*ELONG*AR*AR*BTOR0*AIOTAX
        PRAX = 2. * (1.+ELONG*ELONG) * PSIB * PSIB /
     &         (AR**2*FRV(1,0)*FRV(1,0)*ELONG*ELONG)
        if (insol .eq. 1) then
          do J = 1,NSI
            RHO = (J-0.5) / NSI
            FPP(J) = 2. * PSIB * RHO
            AIOTA(J) = FPP(J) / FTP(J)
            PTH(J) = PRAX * (1.0-RHO*RHO)
C     WRITE (16,60) J, FTP(J), FPP(J), VPMOM(J)
          end do
        else
          do J = 1,NSI
            RHO = (J-0.5) / NSI
            FPP(J) = PSIB
            AIOTA(J) = FPP(J) / FTP(J)
            PTH(J) = PRAX * (1.0-RHO)
C     WRITE (16,60) J, FTP(J), FPP(J), VPMOM(J)
          end do
        end if
      end if
 1022 format (
     &/1x,' EPS',1pe10.3,'  EL',1pe10.3,' NPER',i5,' GAM',1pe10.3,
     &'  ALPHA',1pe10.2,//1x,'LMNV',i3,' LMNL',i3,' LMNB',i3,//)
C--sh-000831 90     FORMAT((6E12.4))
C
      write (16,1022) EPS, EL, NPER, GAM, ALPHA, LMNV, LMNL, LMNB
      write (6,1022) EPS, EL, NPER, GAM, ALPHA, LMNV, LMNL, LMNB
C
C     AUXILIARY TABLE OF MODES FOR LAMBDA CALCULATION
      do mtbl = -36,72
        do ntbl = -72,72
          lfx(mtbl,ntbl) = 0
        end do
      end do
      lss = 0
      do lc = 1,lmnl
        do lr = 1,lmnl
          mxdif = ml(lc) - ml(lr)
          nxdif = (nl(lc)-nl(lr)) / nper
          if (lfx(mxdif,nxdif) .le. 0) then
            lfx(mxdif,nxdif) = 1
            lss = lss + 1
            msx(lss) = mxdif
            nsx(lss) = nxdif * nper
            lsx(mxdif,nxdif) = lss
          end if
        end do
      end do
      do lc = 1,lmnl
        do lr = 1,lmnl
          mxdif = ml(lc) + ml(lr)
          nxdif = (nl(lc)+nl(lr)) / nper
          if (lfx(mxdif,nxdif) .le. 0) then
            lfx(mxdif,nxdif) = 1
            lss = lss + 1
            msx(lss) = mxdif
            nsx(lss) = nxdif * nper
            lsx(mxdif,nxdif) = lss
          end if
        end do
      end do
 1027 format (/,'Number of Modes in the Lambda Reconstruction Table:',i4
     &)
C
      write (6,1027) lss
      write (16,1027) lss
      if (lss > lssl) then
       write(6,1028)lssl,lss
      stop
      end if
 1028 format(' PARAMETER LSSL = ',i6,' IS TOO SMALL; MODIFY IT TO: ',i6)
 1029 format(' PARAMETER NI = ',i6,' MUST BE EQUAL TO: ',i6)
 1030 format(' PARAMETER MLMNB = ',i6,' MUST BE GREATER OR EQUAL TO: ',i
     &6)
 1031 format(' PARAMETER MLMNV = ',i6,' IS TOO SMALL; MODIFY IT TO: ',i6
     &)
      if (NI > NSI) then
       write(6,1029)ni,nsi
      stop
      end if
      if (lmnb > mlmnb) then
       write(6,1030)mlmnb,lmnb
      stop
      end if
      if (lmnv > mlmnv) then
       write(6,1031)mlmnv,lmnv
      stop
      end if
      end subroutine eqinvm
