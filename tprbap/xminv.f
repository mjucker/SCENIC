C****************************************************************************A
      subroutine XMINV(AB,N,ND,SCR,DET,EPS,M,MODE)
C
C     LINPACK SIMULATION OF CRAY-MINV
C
C     AB      AUGMENTED MATRIX  FIRST DIMENSION IS ND
C     N       ORDER OF MATRIX
C     ND      FIRST DIMENSION OF MATRIX AB
C     SCR     SCRATCH VECTOR OF DIMENSION AT LEAST 2*N
C     DET     DETERMINANT
C     EPS     TOLERANCE FOR THE PRODUCT OF PIVOT ELEMENTS
C                         ----> NOT NEEDED HERE
C     M       > 0 NUMBER OF SYSTEMS OF LINEAR EQUATIONS TO SOLVE
C             = 0 DETERMINANT IS COMPUTED
C     MODE    = 1 AB IS OVERWRITTEN WITH INV (AB)
C             = 0 DETERMINANT AS WELL AS INVERSE ARE NOT COMPUTED
C                 ONLY EQUATIONS ARE SOLVED  (IF ANY)
C                 MATRIX IS NOT SAVED IN THIS CASE
C
C--sh-000829      REAL AB(ND,1), SCR(1), DT(2)
C--sh+000829>>>>>>>>>>
C 
C.. Implicits .. 
      implicit none
C 
C.. Formal Arguments .. 
      integer ND,N,SCR(*),M,MODE
      real AB(ND,*),DET,EPS
C 
C.. Local Scalars .. 
      integer I,INFO
C 
C.. Local Arrays .. 
      real DT(2),WORK(ND)
C 
C.. External Calls .. 
      external SGEDI, SGEFA, SGESL
C 
C ... Executable Statements ...
C 
      do i=1,2
         DT(i) = 0.0
      end do
C--sh+000829<<<<<<<<<<      
C
C     PERFORM L-U DECOMPOSITION
C
C     CALL SGECO (AB(1,1), ND, N, SCR, RCOND, SCR(N+1))
      call SGEFA(AB(1,1),ND,N,SCR,INFO)
C
C     SOLVE LINEAR EQUATIONS IF ANY
C
      if (M .ne. 0) then
        do I = 1,M
          call SGESL(AB(1,1),ND,N,SCR,AB(1,N+I),0)
        end do
      end if
      if (MODE .ne. 0) then
C
C     COMPUTE DETERMINANT AND INVERSE
C
C--sh-000831      CALL SGEDI(AB, ND, N, SCR(1), DT, SCR(N+1), 11)
C--sh+000831>>>>>>>>>>
        call SGEDI(AB,ND,N,SCR(1),DT,WORK(1),11)
C--sh+000831<<<<<<<<<<                                        
        DET = DT(1) * (10.0**DT(2))
      end if
      end
