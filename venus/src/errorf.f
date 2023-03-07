      SUBROUTINE ERRORF(X,ERR)
C
C       =========================================
C       Purpose: Compute error function erf(x)
C       Input:   x   --- Argument of erf(x)
C       Output:  ERR --- erf(x)
C       =========================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EPS=1.0D-15
      PI=3.141592653589793D0
      X2=X*X
      IF (DABS(X).LT.3.5D0) THEN
         ER=1.0D0
         R=1.0D0
         DO 10 K=1,50
            R=R*X2/(K+0.5D0)
            ER=ER+R
            IF (DABS(R).LE.DABS(ER)*EPS) GO TO 15
 10      CONTINUE
 15      C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
         ERR=C0*ER
      ELSE
         ER=1.0D0
         R=1.0D0
         DO 20 K=1,12
            R=-R*(K-0.5D0)/X2
 20         ER=ER+R
            C0=DEXP(-X2)/(DABS(X)*DSQRT(PI))
            ERR=1.0D0-C0*ER
            IF (X.LT.0.0) ERR=-ERR
         ENDIF
         END
      
