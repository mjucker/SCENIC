C****************************************************************************A
      subroutine SGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
C--sh-000829      INTEGER LDA,N,IPVT(1),JOB
C--sh-000829      REAL A(LDA,1),DET(2),WORK(1)
C--sh+000829>>>>>>>>>>
C 
C.. Implicits .. 
      implicit none
C 
C.. Formal Arguments .. 
      integer LDA,N,IPVT(*),JOB
      real A(LDA,*),DET(2),WORK(*)
C 
C.. Local Scalars .. 
      integer I,J,K,KB,KP1,L,NM1
      real T,TEN
C 
C.. External Calls .. 
      external SAXPY, SSCAL, SSWAP
C 
C.. Intrinsic Functions .. 
      intrinsic ABS, MOD
C 
C ... Executable Statements ...
C 
C
C
C     COMPUTE DETERMINANT
C
      if (JOB/10 .ne. 0) then
        DET(1) = 1.0e0
        DET(2) = 0.0e0
        TEN = 10.0e0
        do I = 1,N
          if (IPVT(I) .ne. I) then
            DET(1) = -DET(1)
          end if
          DET(1) = A(I,I) * DET(1)
C        ...EXIT
          if (DET(1) .eq. 0.0e0) goto 1000
          do while (ABS(DET(1)) .lt. 1.0e0)
            DET(1) = TEN * DET(1)
            DET(2) = DET(2) - 1.0e0
          end do
          do while (ABS(DET(1)) .ge. TEN)
            DET(1) = DET(1) / TEN
            DET(2) = DET(2) + 1.0e0
          end do
        end do
      end if
C
C     COMPUTE INVERSE(U)
C
 1000 if (MOD(JOB,10) .eq. 0) return
      do K = 1,N
        A(K,K) = 1.0e0 / A(K,K)
        T = -A(K,K)
        call SSCAL(K-1,T,A(1,K),1)
        KP1 = K + 1
        if (N .ge. KP1) then
          do J = KP1,N
            T = A(K,J)
            A(K,J) = 0.0e0
            call SAXPY(K,T,A(1,K),1,A(1,J),1)
          end do
        end if
      end do
C
C        FORM INVERSE(U)*INVERSE(L)
C
      NM1 = N - 1
      if (NM1 .lt. 1) return
      do KB = 1,NM1
        K = N - KB
        KP1 = K + 1
        do I = KP1,N
          WORK(I) = A(I,K)
          A(I,K) = 0.0e0
        end do
        do J = KP1,N
          T = WORK(J)
          call SAXPY(N,T,A(1,J),1,A(1,K),1)
        end do
        L = IPVT(K)
        if (L .ne. K) then
          call SSWAP(N,A(1,K),1,A(1,L),1)
        end if
      end do
      end
