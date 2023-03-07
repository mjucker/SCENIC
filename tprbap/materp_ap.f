!--------0---------0---------0---------0---------0---------0---------0-c
! /19.05.89/ UHS        LAST MODIFIED /17.03.06/ UHS:  TPRM    DATA    D
!--------0---------0---------0---------0---------0---------0---------0-c
!
!                 TERPSICHORE ANISOTROPIC BALLOONING VERSION
!
!                       MASTER VERSION FOR EPFL ITANIUM MACHINE
!-----------------------------------------------------------------------
! /11.01.89/ UHS        LAST MODIFIED /17.03.06/ UHS:  MATERP1 DATA    D
!-----------------------------------------------------------------------
!
      program materp  
!
      use tprb_param
C.. Implicits ..
      implicit none
!
!     MAIN PROGRAM FOR TERPSICHORE
!
!     BY
!
!     D. V. ANDERSON, W. A. COOPER, R. GRUBER AND U. SCHWENN
!
!
!     LINEAR STABILITY ANALYSIS OF 3D EQUILIBRIA
!
!     EQUILIBRIA ARE GIVEN BY A FOURIER REPRESENTATION OF THE GEOMETRY
!     AND THE INVARIANT FUNCTIONS FT(TOROIDAL FLUX), IOTA=FP'/FT',
!     AND THE PRESSURE (USED TO CHECK THE ACCURACY OF THE RECONSTRUCTED
!     EQUILIBRIUM QUANTITIES LIKE P', I', J' AND EQUI (LIKE IN FIT/VMEC)
!
!-----------------------------------------------------------------------
!
!      PARAMETERS FOR GRIDPOINTS, NUMBER OF POINTS ETC.
!
!      THE $$ PARTPR DATA D STATEMENT OR THE EQUIVALENT UPDATE
!      MUST REFER TO A LOWERLEVEL $$ ... OR UPDATE, CONTAINING
!      THE PARAMETERS FOR A SPECIFIC CONFIGURATION
!
!      I.E. $$ PARBAL DATA D ETC. FOR THE SOLOV'EV TESTCASE
!
!      NI     : NUMBER OF RADIAL INTERVALS, NI = NRHO+1 IN VMEC
!      NJ     :   "    "  POLOIDAL INTEGRATION POINTS, NJ >= 3*MAX(M)
!      NK     :   "    "  TOROIDAL      "         "     K          N
!
!      MLMNV  : NUMBER OF R/Z MODES IN THE EQUIL. INPUT, >= LMNV(INPUT)
!      MLMNB  :   "    "        "   FOR THE MAPPING, >= LMNB(COMPUTED
!                                                       FROM THE TABLE)
!
!      BALLOONING PARAMETERS AND DATA:
!      LRADMN : INDEX OF INNER RADIAL SURFACE FOR BALLOONING CALCULATION
!      LRADMX : INDEX OF OUTER RADIAL SURFACE FOR BALLOONING CALCULATION
!      LMXBAL : MAX. POLOIDAL MODE NUMBER OF EXTRA MODES FOR BALLOONING.
!      LSKEWN : SWITCH TO SKEW MODE SELECTION PATTERN IN (M,N) SPACE.
!      LNBALN : MIN. TOROIDAL MODE NUMBER OF EXTRA MODES--DIFF. FOR SKEW
!      LNBALX : MAX. TOROIDAL MODE NUMBER OF EXTRA MODES--DIFF. FOR SKEW
!
!-----------------------------------------------------------------------
!     include'parbal.inc'  
!--------0---------0---------0---------0---------0---------0---------072
      integer :: i
      real :: tim (20)  
!
C
C.. External Calls ..
      external EQINVM, mtaskb, second, veqrec, TPRALL_BAP
!     include 'common_paprpl.inc'
      do i=1,20
         tim(i) = 0.0
      end do
!      close (6)  
C-----------------------------------------------------------------------
C
C    TPRALL_BAL:   INITIALISE BY CALLING THE ALLOCATE ROUTINE
C
C
      call TPRALL_BAP
C
!
!-----------------------------------------------------------------------
!
!     EQINVM: READS R/Z FOURIERCOEFFICIENTS ON INTEGER GRID,
!             READS INVARIANTS AND ALL INFORMATION OF THE EQUILIBRIUM.
!             READS INDEX TABLE FOR FOURIERCOEFFICIENTS USED TO
!                   RECOMPUTE LAMBDA WITH THE SET USED IN BOOZER'S
!                   COORDINATES (BC).
!
!-----------------------------------------------------------------------
!
      call eqinvm  
      call second (tim (1) )  
!-----------------------------------------------------------------------
!
!     VEQREC: Vmec EQuilibrium REConstruction.
!
!-----------------------------------------------------------------------
      call veqrec (tim)  
      call second (tim (4) )  
!-----------------------------------------------------------------------
!
!     MTASKB: PERFORMS THE MAPPING FRM VMEC TO BOOZER COORDINATES.
!             CALL ROUTINES VMTOBO, EXTINT AND METRIC.
!
!     VMTOBO: COMPUTES AC'S Q AND THE MAPPING FUNCTIONS ALFA/GAMMA,
!             PERFORMES THE DIRECT FOURIER INTEGRALS TO OBTAIN THE
!             BOOZER COEFFICIENTS FOR R, Z AND PHI.
!             COMPUTES R, Z AND PHI OM THE HALF GRID
!
!     EXTINT: INTERPOLATES BOOZER FOURIER AMPLITUDES OF R/Z/PHI ONTO
!             INTEGER GRID. CALCULATES THESE VARIABLES IN REAL SPACE.
!
!     PLAVAC: CALCULATES R, Z, PHI AND THEIR DERIVATIVES AT THE
!             PLASMA-VACUUM INTERFACE (PVI).
!
!     CONWAL: CALCULATES R, Z, PHI AND THEIR DERIVATIVES WITH RESPECT
!             TO THE BOOZER ANGLES AT A PRECRIBED CONDUCTING WALL.
!
!     METRIC: COMPUTES UPPER AND LOWER METRIC ELEMENTS IN BC,
!             USING HALF GRID VALUES AND INTEGER GRID VALUES
!             OBTAINED BY INTERPOLATION.
!             COMPUTES S-DERIVATIVES OF BJAC, AT THE BOUNDARY AGAIN
!             BY INTERPOLATION.
!
!     VACMET: CONSTRUCTS THE JACOBIAN AND ALL THE LOWER METRIC ELEMENTS
!             ON A HALF-GRID RADIAL MESH IN THE VACUUM REGION.
!
!     BOPHYS: COMPUTES EQUILIBRIUM QUANTITIES IN BC AND ALL RADIAL
!             FUNCTIONS NEEDED FOR STABILITY.
!
!     MERCIE: COMPUTES MERCIER CRITERION AND COEFFICIENTS OF THE
!             BALLOONING MODE EQUATION.
!
!-----------------------------------------------------------------------
      call mtaskb (tim)  
      call second (tim (9) )  
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  100 continue  
!
      write (6, 1020) tim (1)  
      write (6, 1021) tim (2) - tim (1)  
      write (6, 1022) tim (3) - tim (2)  
      write (6, 1023) tim (4) - tim (1)  
      write (6, 1024) tim (5) - tim (4)  
      write (6, 1025) tim (6) - tim (5)  
      write (6, 1026) tim (7) - tim (6)  
      write (6, 1027) tim (8) - tim (7)  
      write (6, 1028) tim (9) - tim (8)  
      write (6, 1030) tim (9)  
 1020 format   (" timing in subroutine eqinvm =",1p1e15.6, " seconds")  
 1021 format   (" timing in subroutine lgikvm =",1p1e15.6, " seconds")  
 1022 format   (" timing in subroutine mtaskl =",1p1e15.6, " seconds")  
 1023 format   (" timing in subroutine veqrec =",1p1e15.6, " seconds")  
 1024 format   (" timing in subroutine vmtobo =",1p1e15.6, " seconds")  
 1025 format   (" timing in subroutine metric =",1p1e15.6, " seconds")  
 1026 format   (" timing in subroutine bophys =",1p1e15.6, " seconds")  
 1027 format   (" timing in subroutine mercie =",1p1e15.6, " seconds")  
 1028 format   (" timing in subroutine bootsj =",1p1e15.6, " seconds")  
 1030 format   (" total running time          =",1p1e15.6, " seconds")  
!
!end      call endfr
!-----------------------------------------------------------------------
!
      end program materp
