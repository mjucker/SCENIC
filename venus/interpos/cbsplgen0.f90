      SUBROUTINE cbsplgen0(XIN,YIN,nin,xout,yout,youtp,youtpp, &
     &  nout,ioptder,iextrapo,psig,nbc,ybc,IFLAG)
!
      USE prec_rkind
      implicit none
      integer LENREAL
      parameter(LENREAL=8)
      integer nin,nout, nbc(2)
      REAL(RKIND) ::  psig(nin)
      REAL(RKIND) ::  xin(nin), yin(nin), ybc(6)
      REAL(RKIND) ::  xout(nout),yout(nout), youtp(nout), youtpp(nout)
      INTEGER ioptder, NBCLFT,NBCRGT, iextrapo
      REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
      REAL(RKIND) :: PXEXP0,PXEXPDL
!      pointer(iptr_pynew,pynew)
      REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pynew
!      REAL(RKIND) :: pynew(1)
!      pointer(iptr_pyinpp,pyinpp)
      REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pyinpp
!      REAL(RKIND) :: pyinpp(1)
!      pointer(iptr_pamat,pamat)
      REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pamat
!      REAL(RKIND) :: pamat(1)
!      pointer(iptr_pwork,pwork)
      REAL(RKIND), DIMENSION(:), ALLOCATABLE :: pwork
!      REAL(RKIND) :: pwork(1)
!      pointer(iptr_kwork,kwork)
      INTEGER, DIMENSION(:), ALLOCATABLE :: kwork
!      integer kwork(1)
!
      integer iflag, idamat, mdamat, nbytes
!OS      integer*4 malloc_f
!
!
!%OS      NBCLFT = 1
!%OS      NBCRGT = 1
!%OS      YBCLFT = 1.E32_RKIND
!%OS      YBCRGT = 1.E32_RKIND
      NBCLFT = nbc(1)
      NBCRGT = nbc(2)
      YBCLFT = ybc(1)
      YBCRGT = ybc(2)
      if (NBCLFT .ge. 10) XBCLFT = ybc(3)
      if (NBCRGT .ge. 10) XBCRGT = ybc(4)
!
      IDAMAT = 3
      IF (PSIG(1) .EQ. 0._RKIND) IDAMAT = 2
      IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
     &  .OR. NBCRGT.EQ.2) IDAMAT = 3*(IDAMAT-1)+1
      IF (NBCLFT .GE. 10) IDAMAT = IDAMAT + 1
      IF (NBCRGT .GE. 10) IDAMAT = IDAMAT + 2
      mdamat = IDAMAT*nin
      IF( .NOT. ALLOCATED(pynew) ) ALLOCATE(pynew(nin))
      IF( .NOT. ALLOCATED(pyinpp) ) ALLOCATE(pyinpp(nin))
      IF( .NOT. ALLOCATED(pamat) ) ALLOCATE(pamat(mdamat))
      IF( .NOT. ALLOCATED(pwork) ) ALLOCATE(pwork(2*nin))
      IF( .NOT. ALLOCATED(kwork) ) ALLOCATE(kwork(5*nin))
      PXEXP0 = ybc(5)
      PXEXPDL= ybc(6)
      CALL CBSPLGEN(XIN,YIN,PYNEW,PYINPP,Nin,XOUT,YOUT,YOUTP,YOUTPP, &
           &    nout,ioptder,PSIG,KWORK,PWORK,PAMAT,mdamat,NBCLFT, &
           &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL &
           &    ,IFLAG)
!
      END SUBROUTINE cbsplgen0

