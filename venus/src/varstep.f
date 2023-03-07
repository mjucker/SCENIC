
      subroutine varstep(ip,past,v,iota,omega_c,rhocrit,madecol)

      use machine
      use particle
      use kicks
      use scatter

      implicit none

      real omegab,omega_c,rhocrit
      real v,iota,past,twopi,argloc(harmonics)
      integer N,argcheck,i,ip,madecol,check


      argcheck=0
      check=0
      rhocrit=abs(rhocrit)

      twopi=4.*acos(0.)

      N=50  ! not really number of integration steps per orbit

      omegab=v*iota/r00/twopi !estimate of bounce frequency
      
c      if(icrh.gt.0)then
         do i=1,harmonics
            if(ifarg(i,2,ip).eq.0.)then
               argloc(i)=ifarg(i,1,ip)
            else
               argloc(i)=min(ifarg(i,1,ip),ifarg(i,2,ip))
            endif
         enddo


         do i=1,harmonics
            if(abs(argloc(i)).lt.0.01*omega)then
               argcheck=1
               check=1
            endif
         enddo
c      endif

      if(madecol.eq.1.or.rhocrit.lt.1.)check=1

      if(check.eq.1)then
         past=1./(omegab*N)     !time step [s] s.t. one bounce=N steps
      
         past=past*omega_c      !normalisation to non-dimensional time step
      endif


      if (argcheck.eq.1)then
         past=past/5.
      elseif(rhocrit.lt.1.0)then
         past=past/5.           !smaller time step for banana tips
      endif
      
      past=max(past,1.e-7*omega_c)   ! very good results are obtained with 1.e-7*omega_c

      end







C###### variable step including pitch angle & aspect ratio ############C
C###### not continuous over passin-trapped boundary        ############C

c$$$
c$$$      lengthpas=int(tfin*10*normp/mo*iota/
c$$$     $     r00)
c$$$      eps=sqrt(yn(1))*a00/r00
c$$$      pch=(1.-lambda*lambda)*b00/B
c$$$      ksq=2.*eps*pch/(1.-(1.-eps)*pch)
c$$$      corr=sqrt(2.*eps/(2.*eps+(1.-eps)*ksq))
c$$$      ksq=sqrt(ksq)
c$$$      if(ksq.le.1.)then
c$$$         call comelp(ksq,ellipk,ellipe)
c$$$         corr=corr/4./ellipk
c$$$      else
c$$$         call comelp(1./ksq,ellipk,ellipe)
c$$$         corr=corr*ksq/4./ellipk
c$$$      endif
c$$$      corr=1.;
c$$$      lengthpas=max(int(lengthpas*corr),2)
c$$$      
c$$$      if(coulomb.eq.1)then
c$$$         past=min(tfin*norm/lengthpas,0.6*norm/nud) !nudpast<0.9
c$$$      else
c$$$         past=tfin*norm/lengthpas
c$$$      endif
c$$$C########################################################################
c$$$        SUBROUTINE COMELP(HK,CK,CE)
c$$$C
c$$$C       ==================================================
c$$$C       Purpose: Compute complete elliptic integrals K(k)
c$$$C                and E(k)
c$$$C       Input  : K  --- Modulus k ( 0 %GÃ¯Â¿Â½%@k %GÃ¯Â¿Â½%@1 )
c$$$C       Output : CK --- K(k)
c$$$C                CE --- E(k)
c$$$C       ==================================================
c$$$C
c$$$        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$        PK=1.0D0-HK*HK
c$$$        IF (HK.EQ.1.0) THEN
c$$$           CK=1.0D+300
c$$$           CE=1.0D0
c$$$        ELSE
c$$$           AK=(((.01451196212D0*PK+.03742563713D0)*PK
c$$$     &        +.03590092383D0)*PK+.09666344259D0)*PK+
c$$$     &        1.38629436112D0
c$$$           BK=(((.00441787012D0*PK+.03328355346D0)*PK+
c$$$     &        .06880248576D0)*PK+.12498593597D0)*PK+.5D0
c$$$           CK=AK-BK*DLOG(PK)
c$$$           AE=(((.01736506451D0*PK+.04757383546D0)*PK+
c$$$     &        .0626060122D0)*PK+.44325141463D0)*PK+1.0D0
c$$$           BE=(((.00526449639D0*PK+.04069697526D0)*PK+
c$$$     &        .09200180037D0)*PK+.2499836831D0)*PK
c$$$           CE=AE-BE*DLOG(PK)
c$$$        ENDIF
c$$$        RETURN
c$$$        END
