C****************************************************************************A
      subroutine SGESL(A,LDA,N,IPVT,B,JOB)
C--sh-000829      INTEGER LDA,N,IPVT(1),JOB
C--sh-000829      REAL A(LDA,1),B(1)
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
C.. In/Out Status: Maybe Read, Not Written ..
      integer IPVT(*)
C.. In/Out Status: Maybe Read, Maybe Written ..
      real B(*)
C.. In/Out Status: Read, Not Written ..
      integer JOB
C 
C.. Local Scalars .. 
      integer K,KB,L,NM1
      real T
C 
C.. External Functions .. 
      real SDOT
      external SDOT
C 
C.. External Calls .. 
      external SAXPY
C 
C ... Executable Statements ...
C 
C
      NM1 = N - 1
      if (JOB .ne. 0) then
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
        do K = 1,N
          T = SDOT(K-1,A(1,K),1,B(1),1)
          B(K) = (B(K)-T) / A(K,K)
        end do
C
C        NOW SOLVE TRANS(L)*X = Y
C
        if (NM1 .lt. 1) return
        do KB = 1,NM1
          K = N - KB
          B(K) = B(K) + SDOT(N-K,A(K+1,K),1,B(K+1),1)
          L = IPVT(K)
          if (L .ne. K) then
            T = B(L)
            B(L) = B(K)
            B(K) = T
          end if
        end do
      else
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
        if (NM1 .ge. 1) then
          do K = 1,NM1
            L = IPVT(K)
            T = B(L)
            if (L .ne. K) then
              B(L) = B(K)
              B(K) = T
            end if
            call SAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
          end do
        end if
C
C        NOW SOLVE  U*X = Y
C
        do KB = 1,N
          K = N + 1 - KB
          B(K) = B(K) / A(K,K)
          T = -B(K)
          call SAXPY(K-1,T,A(1,K),1,B(1),1)
        end do
      end if
      end
