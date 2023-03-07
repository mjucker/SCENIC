module erf

contains

COMPLEX FUNCTION CW(CZ)
!------------------------------------------------------------------------------
!     EVALUATES w(z)=exp(-z**2)erfc(-i*z) FOR COMPLEX z.
!     THE ALGORITHM IS PRESENTED IN :
!     "efficient computation of the complex error function", W.Gautschi, 
!     SIAM J. Numer. Anal. Vol 7, No.1, March 1970
!     AUTHOR OF PROGRAM : S. Brunner
!------------------------------------------------------------------------------
      PARAMETER (X0=5.33,Y0=4.29,H0=1.6,N0=6,N1=23,NU0=9, &
                NU1=21,ZERO=0.,ONE=1.,HALF=0.5)
      DOUBLE COMPLEX CI,CZ,CZT,R
      DIMENSION HPJ(0:NU0+NU1)
      DATA CI/(0.,1.)/,SQRTPI/1.7724538509055160D0/
!
!-----defining corresponding value in first quadrant z ---> zt ----------------
!
      CZR=REAL(CZ)
      CZI=AIMAG(CZ)
      X=ABS(CZR)
      Y=ABS(CZI)
      CZT=CMPLX(X,Y)
!
!-----evaluating the parameters H, N, NU as a function of zt-------------------
!
      IF (X .LE. X0 .AND. Y .LE. Y0) THEN
         S=(ONE-Y/Y0)*SQRT(ONE-(X/X0)**2)
         H=H0*S
         N=NINT(N0+N1*S)
         NU=NINT(NU0+NU1*S)
      ELSE
         H=ZERO
         N=0
         NU=8
      ENDIF
!
!-----aplying Gautschi's algorithm---------------------------------------------
!
      R=ZERO
      CW=ZERO
      HPJ(0)=1.
      DO 15 J=1,NU
        HPJ(J)=2.*H*HPJ(J-1)
 15   CONTINUE
      DO 20 J=NU,0,-1
        R=HALF/(H-CI*CZT+(J+1)*R)
        IF (J .LE. N) THEN
          CW=R*(HPJ(J)+CW)
        ENDIF
 20   CONTINUE
      CW=2.*CW/SQRTPI
!
!-----evaluating w(z) depending on the quadrant of z
!
      IF (CZR .GE. ZERO) THEN
         IF (CZI .LT. ZERO) THEN
            CW=CMPLX(REAL(CW),-AIMAG(CW))
            CW=2.*EXP(-CZ*CZ)-CW
         ENDIF
      ELSEIF (CZI .GE. ZERO) THEN
         CW=CMPLX(REAL(CW),-AIMAG(CW))
      ELSE 
         CW=2.*EXP(-CZ*CZ)-CW
      ENDIF
!
      END FUNCTION

end module
