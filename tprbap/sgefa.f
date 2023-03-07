C****************************************************************************A
      subroutine SGEFA(A,LDA,N,IPVT,INFO)
C--sh-000829      INTEGER LDA,N,IPVT(1),INFO
C--sh-000829      REAL A(LDA,1)
C--sh+000829>>>>>>>>>>
C 
C.. Implicits .. 
      implicit none
C 
C.. Formal Arguments .. 
C.. In/Out Status: Read, Not Written ..
      integer LDA
C.. In/Out Status: Maybe Read, Maybe Written ..
      real A(LDA,*)
C.. In/Out Status: Read, Not Written ..
      integer N
C.. In/Out Status: Not Read, Maybe Written ..
      integer IPVT(*)
C.. In/Out Status: Not Read, Overwritten ..
      integer INFO
C 
C.. Local Scalars .. 
      integer i,J,K,KP1,L,NM1
      real T
C 
C.. External Functions .. 
      integer ISAMAX
      external ISAMAX
C 
C.. External Calls .. 
      external SSCAL
C 
C ... Executable Statements ...
C 
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      if (NM1 .ge. 1) then
        do K = 1,NM1
          KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
          L = ISAMAX(N-K+1,A(K,K),1) + K - 1
          IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
          if (A(L,K) .eq. 0.0e0) then
            INFO = K
          else
C
C           INTERCHANGE IF NECESSARY
C
            if (L .ne. K) then
              T = A(L,K)
              A(L,K) = A(K,K)
              A(K,K) = T
            end if
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0e0/A(K,K)
            call SSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            do J = KP1,N
              T = A(L,J)
              if (L .ne. K) then
                A(L,J) = A(K,J)
                A(K,J) = T
              end if
C--sh-000725               CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
C--sh+000725>>>>>>>>>>
              if (t .ne. 0.0) then
                do i = k+1,n
                  a(i,j) = a(i,j) + t*a(i,k)
                end do
              end if
C--sh+000725<<<<<<<<<<
C
            end do
          end if
        end do
      end if
      IPVT(N) = N
      if (A(N,N) .eq. 0.0e0) then
        INFO = N
      end if
      end
