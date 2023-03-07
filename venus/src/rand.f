!----------------------------------------------------------------------
      
      SUBROUTINE init_random
      use mpi_var
      use rand_var
      
      IMPLICIT NONE
      
      INTEGER :: i,seed(8)
      CHARACTER*(34) :: string
      character cdate*8,ctime*10,czone*5
      integer time(8)
      
!     Determine and print the seed

      allocate(ran_index(npitch),ran_array(100,npitch))

      IF (me2.eq.0) THEN
cj  for reproducable results, same seed at every run (random_seed_str)
c       random_seed_str= "3.7182813"
c       print*,'Fixed seed ',random_seed_str
c       CALL decimal_to_seed(random_seed_str, seed)
c       CALL seed_to_decimal(seed, string)
cj  for different seed at every run
        call date_and_time(cdate,ctime,czone,time)
        call set_random_seed(time,seed)
        call seed_to_decimal(seed,string)
        write(*,'(1x,a,8i4)')'Arbitrary seed ',time
      END IF
	      
      CALL mpi_bcast(seed, 8, mpi_integer, 0, comm, ierr)  
      
!     Each processor initializes all random number series relative to 
!     different bunches
cj  # bunches = npitch
      CALL next_seed(me2*npitch, seed)
      
      DO i = 1, npitch
         CALL random_init(seed, ran_index(i), ran_array(:,i))
         CALL next_seed(1, seed)
      END DO
      
      RETURN
      
      END SUBROUTINE init_random
      

!----------------------------------------------------------------------
* Version 1.0 of random number routines
* Author: Charles Karney <karney@princeton.edu>
* Date: 1999-08-05 10:44:51 -0400
*
      function random(p,r)
      implicit none
      DOUBLE PRECISION ulp2
      parameter(ulp2=2.0d0**(-47-1))
      REAL ulps
      DOUBLE PRECISION mult
      parameter(ulps=2.0**(-23),mult=2.0d0**23)
      DOUBLE PRECISION random
      REAL srandom
      integer p
      DOUBLE PRECISION r(0:100-1)
      external rand_batch
      if(p.GE.100)then
      call rand_batch(p,r(0))
      end if
      random=r(p)+ulp2
      p=p+1
      return
      entry srandom(p,r)
      if(p.GE.100)then
      call rand_batch(p,r(0))
      end if
      srandom=(int(mult*r(p))+0.5)*ulps
      p=p+1
      return
      end
*
      subroutine random_array(y,n,p,r)
      implicit none
      DOUBLE PRECISION ulp2
      parameter(ulp2=2.0d0**(-47-1))
      REAL ulps
      DOUBLE PRECISION mult
      parameter(ulps=2.0**(-23),mult=2.0d0**23)
      integer n
      DOUBLE PRECISION y(0:n-1)
      REAL ys(0:n-1)
      integer p
      DOUBLE PRECISION r(0:100-1)
      integer i,k,j
      external rand_batch
      if(n.LE.0)return
      k=min(n,100-p)
      do i=0,k-1
      y(i)=r(i+p)+ulp2
      end do
      p=p+(k)
      do j=k,n-1,100
      call rand_batch(p,r(0))
      do i=j,min(j+100,n)-1
      y(i)=r(i-j+p)+ulp2
      end do
      p=p+(min(100,n-j))
      end do
      return
      entry srandom_array(ys,n,p,r)
      if(n.LE.0)return
      k=min(n,100-p)
      do i=0,k-1
      ys(i)=(int(mult*r(i+p))+0.5)*ulps
      end do
      p=p+(k)
      do j=k,n-1,100
      call rand_batch(p,r(0))
      do i=j,min(j+100,n)-1
      ys(i)=(int(mult*r(i-j+p))+0.5)*ulps
      end do
      p=p+(min(100,n-j))
      end do
      return
      end
*
      subroutine rand_batch(p,r)
      implicit none
      integer p
      DOUBLE PRECISION r(0:100-1)
      integer i
      DOUBLE PRECISION w(0:1009-100-1)
      DOUBLE PRECISION tmp
      do i=0,63-1
      tmp=r(i)+r(i+100-63)
      w(i)=tmp-int(tmp)
      end do
      do i=63,100-1
      tmp=r(i)+w(i-63)
      w(i)=tmp-int(tmp)
      end do
      do i=100,1009-100-1
      tmp=w(i-100)+w(i-63)
      w(i)=tmp-int(tmp)
      end do
      do i=1009-100,1009-100+63-1
      tmp=w(i-100)+w(i-63)
      r(i-1009+100)=tmp-int(tmp)
      end do
      do i=1009-100+63,1009-1
      tmp=w(i-100)+r(i-1009+100-63)
      r(i-1009+100)=tmp-int(tmp)
      end do
      p=0
      return
      end
*
      subroutine random_init(seed,p,r)
      implicit none
      integer b
      DOUBLE PRECISION del,ulp
      parameter(del=2.0d0**(-14),ulp=2.0d0**(-47),b=2**14)
      integer seed(0:8-1)
      integer p
      DOUBLE PRECISION r(0:100-1)
      integer i,j,t,s(0:8-1)
      logical even
      integer z(0:8-1)
      integer a0,a1,a2,a3,a4,a5,a6,c0
      data a0,a1,a2,a3,a4,a5,a6,c0/15661,678,724,5245,13656,11852,29,1/
      do i=0,8-1
      s(i)=seed(i)
      end do
      even=.TRUE.
      even=even.AND.(mod(s(7),2).EQ.0)
      r(0)=(((s(7)*del+s(6))*del+s(5))*del+int(s(4)/512))*512*del
      do j=1,100-1
      t=c0+a0*s(0)
      z(0)=mod(t,b)
      t=int(t/b)+a0*s(1)+a1*s(0)
      z(1)=mod(t,b)
      t=int(t/b)+a0*s(2)+a1*s(1)+a2*s(0)
      z(2)=mod(t,b)
      t=int(t/b)+a0*s(3)+a1*s(2)+a2*s(1)+a3*s(0)
      z(3)=mod(t,b)
      t=int(t/b)+a0*s(4)+a1*s(3)+a2*s(2)+a3*s(1)+a4*s(0)
      z(4)=mod(t,b)
      t=int(t/b)+a0*s(5)+a1*s(4)+a2*s(3)+a3*s(2)+a4*s(1)+a5*s(0)
      z(5)=mod(t,b)
      t=int(t/b)+a0*s(6)+a1*s(5)+a2*s(4)+a3*s(3)+a4*s(2)+a5*s(1)+a6*s(0)
      z(6)=mod(t,b)
      t=int(t/b)+a0*s(7)+a1*s(6)+a2*s(5)+a3*s(4)+a4*s(3)+a5*s(2)+a6*s(1)
      z(7)=mod(t,b)
      do i=0,8-1
      s(i)=z(i)
      end do
      even=even.AND.(mod(s(7),2).EQ.0)
      r(j)=(((s(7)*del+s(6))*del+s(5))*del+int(s(4)/512))*512*del
      end do
      if(even)then
      t=c0+a0*s(0)
      z(0)=mod(t,b)
      t=int(t/b)+a0*s(1)+a1*s(0)
      z(1)=mod(t,b)
      t=int(t/b)+a0*s(2)+a1*s(1)+a2*s(0)
      z(2)=mod(t,b)
      t=int(t/b)+a0*s(3)+a1*s(2)+a2*s(1)+a3*s(0)
      z(3)=mod(t,b)
      t=int(t/b)+a0*s(4)+a1*s(3)+a2*s(2)+a3*s(1)+a4*s(0)
      z(4)=mod(t,b)
      t=int(t/b)+a0*s(5)+a1*s(4)+a2*s(3)+a3*s(2)+a4*s(1)+a5*s(0)
      z(5)=mod(t,b)
      t=int(t/b)+a0*s(6)+a1*s(5)+a2*s(4)+a3*s(3)+a4*s(2)+a5*s(1)+a6*s(0)
      z(6)=mod(t,b)
      t=int(t/b)+a0*s(7)+a1*s(6)+a2*s(5)+a3*s(4)+a4*s(3)+a5*s(2)+a6*s(1)
      z(7)=mod(t,b)
      do i=0,8-1
      s(i)=z(i)
      end do
      j=int((s(8-1)*100)/b)
      r(j)=r(j)+(ulp)
      end if
      p=100
      return
      end
*
      subroutine decimal_to_seed(decimal,seed)
      implicit none
      integer b
      parameter(b=2**14)
      character*(*)decimal
      integer seed(0:8-1)
      integer i,ten(0:8-1),c(0:8-1),ch
      data ten/10,7*0/
      external rand_axc
      do i=0,8-1
      seed(i)=0
      c(i)=0
      end do
      do i=1,len(decimal)
      ch=ichar(decimal(i:i))
      if(ch.GE.ichar('0').AND.ch.LE.ichar('9'))then
      c(0)=ch-ichar('0')
      call rand_axc(ten,seed,c)
      end if
      end do
      return
      end
*
      subroutine string_to_seed(string,seed)
      implicit none
      integer b
      parameter(b=2**14)
      character*(*)string
      integer seed(0:8-1)
      integer t,i,k,unity(0:8-1),c(0:8-1),ch
      data unity/1,7*0/
      external rand_axc
      do i=0,8-1
      seed(i)=0
      c(i)=0
      end do
      do i=1,len(string)
      ch=ichar(string(i:i))
      if(ch.GT.ichar(' ').AND.ch.LT.127)then
      t=mod(seed(0),2)*(b/2)
      do k=0,8-1
      seed(k)=int(seed(k)/2)
      if(k.LT.8-1)then
      seed(k)=seed(k)+(mod(seed(k+1),2)*(b/2))
      else
      seed(k)=seed(k)+(t)
      end if
      end do
      c(0)=ch
      call rand_axc(unity,seed,c)
      end if
      end do
      return
      end
*
      subroutine set_random_seed(time,seed)
      implicit none
      integer time(8)
      integer seed(0:8-1)
      character*21 c
      external decimal_to_seed
      write(c(1:8),'(i4.4,2i2.2)')time(1),time(2),time(3)
      write(c(9:12),'(i1.1,i3.3)') (1-sign(1,time(4)))/2,abs(time(4))
      write(c(13:21),'(3i2.2,i3.3)')time(5),time(6),time(7),time(8)
      call decimal_to_seed(c,seed)
      return
      end
*
      subroutine seed_to_decimal(seed,decimal)
      implicit none
      integer pow,decbase,b
      parameter(pow=4,decbase=10**pow,b=2**14)
      character*(*)decimal
      integer seed(0:8-1)
      integer z(0:8-1),i,t,j,k
      character*36 str
      k=-1
      do i=0,8-1
      z(i)=seed(i)
      if(z(i).GT.0)k=i
      end do
      str=' '
      i=9
90000 continue
      i=i-1
      t=0
      do j=k,0,-1
      z(j)=z(j)+t*b
      t=mod(z(j),decbase)
      z(j)=int(z(j)/decbase)
      end do
      if(z(max(0,k)).EQ.0)k=k-1
      if(k.GE.0)then
      write(str(pow*i+1:pow*(i+1)),'(i4.4)')t
      else
      write(str(pow*i+1:pow*(i+1)),'(i4)')t
      end if
      if(k.GE.0)goto 90000
      if(len(decimal).GT.len(str))then
      decimal(:len(decimal)-len(str))=' '
      decimal(len(decimal)-len(str)+1:)=str
      else
      decimal=str(len(str)-len(decimal)+1:)
      end if
      return
      end
*
      subroutine rand_next_seed(n,ax,cx,y)
      implicit none
      integer n,ax(0:8-1),cx(0:8-1)
      integer y(0:8-1)
      integer a(0:8-1),c(0:8-1),z(0:8-1),t(0:8-1),m,i
      data z/8*0/
      external rand_axc
      if(n.EQ.0)return
      m=n
      do i=0,8-1
      a(i)=ax(i)
      c(i)=cx(i)
      end do
90000 continue
      if(mod(m,2).GT.0)then
      call rand_axc(a,y,c)
      end if
      m=int(m/2)
      if(m.EQ.0)return
      do i=0,8-1
      t(i)=c(i)
      end do
      call rand_axc(a,c,t)
      do i=0,8-1
      t(i)=a(i)
      end do
      call rand_axc(t,a,z)
      goto 90000
      end
*
      subroutine next_seed3(n0,n1,n2,seed)
      implicit none
      integer n0,n1,n2
      integer seed(0:8-1)
      integer af0(0:8-1),cf0(0:8-1)
      integer ab0(0:8-1),cb0(0:8-1)
      integer af1(0:8-1),cf1(0:8-1)
      integer ab1(0:8-1),cb1(0:8-1)
      integer af2(0:8-1),cf2(0:8-1)
      integer ab2(0:8-1),cb2(0:8-1)
      data af0/15741,8689,9280,4732,12011,7130,6824,12302/
      data cf0/16317,10266,1198,331,10769,8310,2779,13880/
      data ab0/9173,9894,15203,15379,7981,2280,8071,429/
      data cb0/8383,3616,597,12724,15663,9639,187,4866/
      data af1/8405,4808,3603,6718,13766,9243,10375,12108/
      data cf1/13951,7170,9039,11206,8706,14101,1864,15191/
      data ab1/6269,3240,9759,7130,15320,14399,3675,1380/
      data cb1/15357,5843,6205,16275,8838,12132,2198,10330/
      data af2/445,10754,1869,6593,385,12498,14501,7383/
      data cf2/2285,8057,3864,10235,1805,10614,9615,15522/
      data ab2/405,4903,2746,1477,3263,13564,8139,2362/
      data cb2/8463,575,5876,2220,4924,1701,9060,5639/
      external rand_next_seed
      if(n2.GT.0)then
      call rand_next_seed(n2,af2,cf2,seed)
      else if(n2.LT.0)then
      call rand_next_seed(-n2,ab2,cb2,seed)
      end if
      if(n1.GT.0)then
      call rand_next_seed(n1,af1,cf1,seed)
      else if(n1.LT.0)then
      call rand_next_seed(-n1,ab1,cb1,seed)
      end if
      entry next_seed(n0,seed)
      if(n0.GT.0)then
      call rand_next_seed(n0,af0,cf0,seed)
      else if(n0.LT.0)then
      call rand_next_seed(-n0,ab0,cb0,seed)
      end if
      return
      end
*
      subroutine rand_axc(a,x,c)
      implicit none
      integer b
      parameter(b=2**14)
      integer a(0:8-1),c(0:8-1)
      integer x(0:8-1)
      integer z(0:8-1),i,j,t
      t=0
      do i=0,8-1
      t=c(i)+int(t/b)
      do j=0,i
      t=t+(a(j)*x(i-j))
      end do
      z(i)=mod(t,b)
      end do
      do i=0,8-1
      x(i)=z(i)
      end do
      return
      end
*
!-----------------------------------------------------------------------
* Version 1.1 of random number routines
* Author: Charles Karney <karney@princeton.edu>
* Date: 1999-08-06 14:18:03 -0400
*
      subroutine random_gauss(y,n,p,r)
      implicit none
      integer p
      DOUBLE PRECISION r(0:100-1)
      integer n
      DOUBLE PRECISION y(0:n-1)
      integer i
      DOUBLE PRECISION pi,theta,z
      DOUBLE PRECISION random
      external random_array,random
      data pi/3.14159265358979323846264338328d0/
      REAL ys(0:n-1)
      REAL spi,stheta,sz
      REAL srandom
      external srandom_array,srandom
      data spi/3.14159265358979323846264338328/
      if(n.LE.0)return
      call random_array(y,n,p,r(0))
      do i=0,int(n/2)*2-1,2
      theta=pi*(2.0d0*y(i)-1.0d0)
      z=sqrt(-2.0d0*log(y(i+1)))
      y(i)=z*cos(theta)
      y(i+1)=z*sin(theta)
      end do
      if(mod(n,2).EQ.0)return
      theta=pi*(2.0d0*y(n-1)-1.0d0)
      z=sqrt(-2.0d0*log(random(p,r(0))))
      y(n-1)=z*cos(theta)
      return
      entry srandom_gauss(ys,n,p,r)
      if(n.LE.0)return
      call srandom_array(ys,n,p,r(0))
      do i=0,int(n/2)*2-1,2
      stheta=spi*(2.0*ys(i)-1.0)
      sz=sqrt(-2.0*log(ys(i+1)))
      ys(i)=sz*cos(stheta)
      ys(i+1)=sz*sin(stheta)
      end do
      if(mod(n,2).EQ.0)return
      stheta=spi*(2.0*ys(n-1)-1.0)
      sz=sqrt(-2.0*srandom(p,r(0)))
      ys(n-1)=sz*cos(stheta)
      return
      end
*
      subroutine random_isodist(v,n,p,r)
      implicit none
      integer p
      DOUBLE PRECISION r(0:100-1)
      integer n
      DOUBLE PRECISION v(0:3*n-1)
      integer i
      DOUBLE PRECISION pi,costheta,phi
      external random_array
      data pi/3.14159265358979323846264338328d0/
      REAL vs(0:3*n-1)
      REAL spi,scostheta,sphi
      external srandom_array
      data spi/3.14159265358979323846264338328/
      if(n.LE.0)return
      call random_array(v(n),2*n,p,r(0))
      do i=0,n-1
      costheta=2.0d0*v(n+2*i)-1.0d0
      phi=pi*(2.0d0*v(n+2*i+1)-1.0d0)
      v(3*i)=cos(phi)*sqrt(1.0d0-costheta**2)
      v(3*i+1)=sin(phi)*sqrt(1.0d0-costheta**2)
      v(3*i+2)=costheta
      end do
      return
      entry srandom_isodist(vs,n,p,r)
      if(n.LE.0)return
      call srandom_array(vs(n),2*n,p,r(0))
      do i=0,n-1
      scostheta=2.0*vs(n+2*i)-1.0
      sphi=spi*(2.0*vs(n+2*i+1)-1.0)
      vs(3*i)=cos(sphi)*sqrt(1.0-scostheta**2)
      vs(3*i+1)=sin(sphi)*sqrt(1.0-scostheta**2)
      vs(3*i+2)=scostheta
      end do
      return
      end
*
      subroutine random_cosdist(v,n,p,r)
      implicit none
      integer p
      DOUBLE PRECISION r(0:100-1)
      integer n
      DOUBLE PRECISION v(0:3*n-1)
      integer i
      DOUBLE PRECISION pi,costheta2,phi
      external random_array
      data pi/3.14159265358979323846264338328d0/
      REAL vs(0:2*n-1)
      REAL spi,scostheta2,sphi
      external srandom_array
      data spi/3.14159265358979323846264338328/
      if(n.LE.0)return
      call random_array(v(n),2*n,p,r(0))
      do i=0,n-1
      costheta2=v(n+2*i)
      phi=pi*(2.0d0*v(n+2*i+1)-1.0d0)
      v(3*i)=cos(phi)*sqrt(1.0d0-costheta2)
      v(3*i+1)=sin(phi)*sqrt(1.0d0-costheta2)
      v(3*i+2)=sqrt(costheta2)
      end do
      return
      entry srandom_cosdist(vs,n,p,r)
      if(n.LE.0)return
      call srandom_array(vs(n),2*n,p,r(0))
      do i=0,n-1
      scostheta2=vs(n+2*i)
      sphi=spi*(2.0*vs(n+2*i+1)-1.0)
      vs(3*i)=cos(sphi)*sqrt(1.0-scostheta2)
      vs(3*i+1)=sin(sphi)*sqrt(1.0-scostheta2)
      vs(3*i+2)=sqrt(scostheta2)
      end do
      return
      end
*

