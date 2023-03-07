C-----------------------------------------------------------------------
C  15.12.89        LAST MODIFICATION 20.03.06     UHS: vforce DATA D
C-----------------------------------------------------------------------
C
      subroutine vforce(ni,njk,lmnb,nper,lnbaln,lnbalx,mm,lmxbal,nmin,
     &            nmax,lskewn,s,vvp,ftp,fpp,civ,cjv,civp,cjvp,pth,pvpi,
     &            equiv,pvp,mb,nb,pol,tor,tcos,tsin,lfrz,parpvi,parpav)
C
C     CALCULATE THE EQUILIBRIUM FORCE BALANCE IN VMEC COORDINATES.
C     PRINT THESE AND OTHER FLUX SURFACE QUANTITIES.
C
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
      integer njk,nmin,nmax,ni,lmnb,ivac,lskewn,nper,lnbaln,lnbalx,
     &        mm,lmxbal,mb(*),nb(*),lfrz(0:36,nmin:nmax)
      real dsvac,s(0:*),vvp(*),ftp(*),fpp(*),civ(*),cjv(*),civp(*),
     &     cjvp(*),pth(*),pvpi(*),equiv(*),pvp(*),pol(*),tor(*),
     &     tcos(njk,*),tsin(njk,*)
      real, dimension(ni) :: parpav ,parpvi
C
C.. Local Scalars ..
      integer I,jk,l,m,n, mcent, nminb, nmaxb
      real EQUIN,FPPI,FTPI,tsi,VVPI,rmcent
C
C.. Intrinsic Functions ..
      intrinsic ABS, COS, SIN,  MIN
C
C ... Executable Statements ...
C
C
C
      do 320 I = 1,NI-1
C
        tsi = 2. / (s(i+1)-s(i-1))
        VVPI = 0.5 * (VVP(I+1)+VVP(I))
        FTPI = 0.5 * (FTP(I+1)+FTP(I))
        FPPI = 0.5 * (FPP(I+1)+FPP(I))
        CIVP(I)  = TSI * (CIV(I+1)-CIV(I))
        CJVP(I)  = TSI * (CJV(I+1)-CJV(I))
        PVPI(I)  = TSI * (PTH(I+1)-PTH(I))
        PARPVI(I)= 0.5 * (PARPAV(I+1)+PARPAV(I))
        EQUIN    = ABS(CJVP(I)*FPPI) + ABS(CIVP(I)*FTPI) +
     &             ABS(VVPI*PVPI(I)) + ABS(PARPVI(I))
        EQUIV(I) = ((CJVP(I)*FPPI-CIVP(I)*FTPI)+PARPVI(I)-VVPI*PVPI(I))
     &           / EQUIN
        PVP(I) = (CJVP(I)*FPPI-CIVP(I)*FTPI+PARPVI(I)) / VVPI
 320  end do
 1002 format (
     &//1x,'***** LGIKVM, RECONSTRUCTED EQUILIBRIUM *****',//3x,'I',9x,
     &'VVP(I)',5x,'PVP(I)',4x,'PVPI(I)',4x,'PARPVI(I)',2x,'CJVP(I)',4x,
     &'CIVP(I)',5x,'EQUI',/)
C
      write (6,1002)
      write (16,1002)
      do I = 1,NI
 1003   format (1x,i3,1pe16.7,1p6e11.4)
        write (6,1003)
     & I, VVP(I),PVP(I),PVPI(I),-PARPVI(I)/VVPI,CJVP(I),CIVP(I),EQUIV(I)
        write (16,1003)
     & I, VVP(I),PVP(I),PVPI(I),-PARPVI(I)/VVPI,CJVP(I),CIVP(I),EQUIV(I)
      end do
C
C     COMPUTE COS/SIN FOR BC SET
C
C     EXTRA MODES FOR BALLOONING
      nminb = lnbaln
      nmaxb = lnbalx
      do 113 m=mm+1,lmxbal
        if (lskewn > 0) then
        lskewn = min(lskewn,ni)
        rmcent = m * fpp(lskewn) / (ftp(lskewn)*nper)
        mcent  = nint(rmcent)
        nminb = mcent + lnbaln
        nmaxb = mcent + lnbalx
        endif
         do 112 n=nminb,nmaxb
           lmnb = lmnb + 1
           mb(lmnb) = m
           nb(lmnb) = n * nper
 112     end do
 113  end do    
C
      do l = 1,lmnb
        do jk = 1,njk
          TCOS(JK,L) = COS(MB(L)*POL(JK)-NB(L)*TOR(JK))
          TSIN(JK,L) = SIN(MB(L)*POL(JK)-NB(L)*TOR(JK))
        end do
      end do
C
C     reinitialise lfrz array
C
      do m = 0,36
        do n = nmin,nmax
          lfrz(m,n) = 0
        end do
      end do
!
      return
      end subroutine vforce
