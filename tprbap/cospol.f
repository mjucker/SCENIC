C-----------------------------------------------------------------------
C  15.12.89        LAST MODIFICATION 01.02.94     UHS: cospol DATA D
C-----------------------------------------------------------------------
C
      subroutine cospol(nper,nj,njk,lmnv,lss,lmap,tpi,dth,dph,mx,nx,msx,
     &                  nsx,tor,pol,tsin,tcos,tsc)
C
C
C 
C.. Implicits .. 
      implicit none
C 
C.. Formal Arguments .. 
      integer njk,lss,nper,nj,lmnv,lmap,mx(*),nx(*),msx(*),nsx(*)
      real tpi,dth,dph,tor(*),pol(*),tsin(njk,*),tcos(njk,*),tsc(lss,*)
C 
C.. Local Scalars .. 
      integer JK,L,LS
C 
C.. External Calls .. 
      external ploteq
C 
C.. Intrinsic Functions .. 
      intrinsic COS, MOD, SIN
C 
C ... Executable Statements ...
C 
C     COMPUTES POLOIDAL AND TOROIDAL ANGLES AND TRIGONOMETRIC FUNCTIONS.    
C
C     WRITE(6,1005) NI,DS,TSI,TDS,TPONP
C1005 FORMAT(/1X'LGIKVM: NI, DS, TSI, TDS, TPONP =',I3,1P4E12.4,/)
C
C
C     PRECOMPUTE COS/SIN FOR EQUILIBRIUM M/N'S
C
      do JK = 1,NJK
C
C        POL(JK) = TPI * DTH * (MOD((JK+NJ),NJ)-1)
        POL(JK) = TPI * DTH * MOD((JK-1),NJ)
        TOR(JK) = TPI * DPH * ((JK-1)/NJ) / NPER
        do L = 1,LMNV
          TCOS(JK,L) = COS(MX(L)*POL(JK)-NX(L)*TOR(JK))
          TSIN(JK,L) = SIN(MX(L)*POL(JK)-NX(L)*TOR(JK))
        end do
      end do
C
C
C     WRITE EQUILIBRIUM COEFFICIENTS TO UNIT 17 FOR PLOTTING
C
      lmap = 1
      call ploteq(lmap)
C
      do LS = 1,LsS
C
        do JK = 1,NJK
          TSC(ls,jk) = COS(MSx(LS)*POL(jk)-NSx(LS)*TOR(jk))
        end do
C
      end do
      end
