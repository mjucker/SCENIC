      SUBROUTINE airy(x,ai,bi,aip,bip)
      REAL ai,aip,bi,bip,x
C USES bessik,bessjy
c     Returns Airy functions Ai(x), Bi(x), and their derivatives Ai (x), Bi (x).
      REAL absx,ri,rip,rj,rjp,rk,rkp,rootx,ry,ryp,z,
     $     PI,THIRD,TWOTHR,ONOVRT
      PARAMETER (PI=3.1415927,THIRD=1./3.,TWOTHR=2.*THIRD,
     $     ONOVRT=.57735027)
      absx=abs(x)
      rootx=sqrt(absx)
      z=TWOTHR*absx*rootx
      if(x.gt.0.)then
         call bessik(z,THIRD,ri,rk,rip,rkp)
         ai=rootx*ONOVRT*rk/PI
         bi=rootx*(rk/PI+2.*ONOVRT*ri)
         call bessik(z,TWOTHR,ri,rk,rip,rkp)
         aip=-x*ONOVRT*rk/PI
         bip=x*(rk/PI+2.*ONOVRT*ri)
      else if(x.lt.0.)then
         call bessjy(z,THIRD,rj,ry,rjp,ryp)
         ai=.5*rootx*(rj-ONOVRT*ry)
         bi=-.5*rootx*(ry+ONOVRT*rj)
         call bessjy(z,TWOTHR,rj,ry,rjp,ryp)
         aip=.5*absx*(ONOVRT*ry+rj)
         bip=.5*absx*(ONOVRT*rj-ry)
      else                      !Case x = 0.
         ai=.35502805
         bi=ai/ONOVRT
         aip=-.25881940
         bip=-aip/ONOVRT
      endif
      return
      END


c-------------------------------------------------------------------------------
      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi) 
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8)
C USES chebev
c      Evaluates Γ1 and Γ2 by Chebyshev expansion for |x| ≤ 1/2. Also returns 1/Γ(1 + x) and
c      1/Γ(1 − x). If converting to double precision, set NUSE1 = 7, NUSE2 = 8.
      REAL xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371168d0,6.5165112670737d-3,
     $     3.087090173086d-4,-3.4706269649d-6,6.9437664d-9,
     $     3.67795d-11,-1.356d-13/
      DATA c2/1.843740587300905d0,-7.68528408447867d-2,
     $     1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8,
     $     2.423096d-10,-1.702d-13,-1.49d-15/
      xx=8.d0*x*x-1.d0                         !Multiply x by 2 to make range be −1 to 1, and then
      gam1=chebev(-1.,1.,c1,NUSE1,xx)           !   apply transformation for evaluating even Cheby-
      gam2=chebev(-1.,1.,c2,NUSE2,xx)            !  shev series.
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END
c-------------------------------------------------------------------------------
      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      INTEGER MAXIT
      REAL rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-16,FPMIN=1.e-30,MAXIT=10000,XMIN=2.,
     $     PI=3.141592653589793d0)
C USES beschb
c     Returns the Bessel functions rj = Jν , ry = Yν and their derivatives rjp = Jν , ryp = Yν ,
c     for positive x and for xnu = ν ≥ 0. The relative accuracy is within one or two signiﬁcant
c     digits of EPS, except near a zero of one of the functions, where EPS controls its absolute
c     accuracy. FPMIN is a number close to the machine’s smallest ﬂoating-point number. All
c     internal arithmetic is in double precision. To convert the entire routine to double precision,
c     change the REAL declaration above and decrease EPS to 10−16. Also convert the subroutine
c     beschb.
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,
     $     dr,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,
     $     p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,
     $     rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'x<0 in bessjy'
      if(x.lt.XMIN)then         !nl is the number of downward recurrences of the J’s and
         nl=int(xnu+.5d0)       !upward recurrences of Y ’s. xmu lies between −1/2 and
      else                      !1/2 for x < XMIN, while it is chosen so that x is greater
         nl=max(0,int(xnu-x+1.5d0)) !than the turning point for x ≥ XMIN.
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI                  !The Wronskian.
      isign=1                   ! Evaluate CF1 by modiﬁed Lentz’s method (§5.2). isign keeps
      h=xnu*xi                  !     track of sign changes in the denominator.
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
         b=b+xi2
         d=b-d
         if(abs(d).lt.FPMIN)d=FPMIN
         c=b-1.d0/c
         if(abs(c).lt.FPMIN)c=FPMIN
         d=1.d0/d
         del=c*d
         h=del*h
         if(d.lt.0.d0)isign=-isign
         if(abs(del-1.d0).lt.EPS)goto 1
 11   enddo
      pause 'x too large in bessjy; try asymptotic expansion'
 1    continue
      rjl=isign*FPMIN           !Initialize Jν and Jν for downward recurrence.
      rjpl=h*rjl
      rjl1=rjl                  !Store values for later rescaling.
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
         rjtemp=fact*rjl+rjpl
         fact=fact-xi
         rjpl=fact*rjtemp-rjl
         rjl=rjtemp
 12   enddo
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl                !Now have unnormalized Jµ and Jµ .
      if(x.lt.XMIN) then        !Use series.
         x2=.5d0*x
         pimu=PI*xmu
         if(abs(pimu).lt.EPS)then
            fact=1.d0
         else
            fact=pimu/sin(pimu)
         endif
         d=-log(x2)
         e=xmu*d
         if(abs(e).lt.EPS)then
            fact2=1.d0
         else
            fact2=sinh(e)/e
         endif
         call beschb(xmu,gam1,gam2,gampl,gammi) !Chebyshev evaluation of Γ1 and Γ2 .
         ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d) !f0 .
         e=exp(e)
         p=e/(gampl*PI)         !p0 .
         q=1.d0/(e*PI*gammi)    !q0 .
         pimu2=0.5d0*pimu
         if(abs(pimu2).lt.EPS)then
            fact3=1.d0
         else
            fact3=sin(pimu2)/pimu2
         endif
         r=PI*pimu2*fact3*fact3
         c=1.d0
         d=-x2*x2
         sum=ff+r*q
         sum1=p
         do 13 i=1,MAXIT
            ff=(i*ff+p+q)/(i*i-xmu2)
            c=c*d/i
            p=p/(i-xmu)
            q=q/(i+xmu)
            del=c*(ff+r*q)
            sum=sum+del
            del1=c*p-i*del
            sum1=sum1+del1
            if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
 13      enddo
         pause 'bessy series failed to converge'
 2       continue
         rymu=-sum
         ry1=-sum1*xi2
         rymup=xmu*xi*rymu-ry1
         rjmu=w/(rymup-f*rymu)  !Equation (6.7.13).
      else                      !Evaluate CF2 by modiﬁed Lentz’s method
         a=.25d0-xmu2           !    (§5.2).
         p=-.5d0*xi
         q=1.d0
         br=2.d0*x
         bi=2.d0
         fact=a*xi/(p*p+q*q)
         cr=br+q*fact
         ci=bi+p*fact
         den=br*br+bi*bi
         dr=br/den
         di=-bi/den
         dlr=cr*dr-ci*di
         dli=cr*di+ci*dr
         temp=p*dlr-q*dli
         q=p*dli+q*dlr
         p=temp
         do 14 i=2,MAXIT
            a=a+2*(i-1)
            bi=bi+2.d0
            dr=a*dr+br
            di=a*di+bi
            if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
            fact=a/(cr*cr+ci*ci)
            cr=br+cr*fact
            ci=bi-ci*fact
            if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
            den=dr*dr+di*di
            dr=dr/den
            di=-di/den
            dlr=cr*dr-ci*di
            dli=cr*di+ci*dr
            temp=p*dlr-q*dli
            q=p*dli+q*dlr
            p=temp
            if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
 14      enddo
         pause 'cf2 failed in bessjy'
 3       continue
         gam=(p-f)/q            !Equations (6.7.6) – (6.7.10).
         rjmu=sqrt(w/((p-f)*gam+q))
         rjmu=sign(rjmu,rjl)
         rymu=rjmu*gam
         rymup=rymu*(p+q/gam)
         ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact              !Scale original Jν and Jν .
      rjp=rjp1*fact
      do 15 i=1,nl              !Upward recurrence of Yν .
         rytemp=(xmu+i)*xi2*ry1-rymu
         rymu=ry1
         ry1=rytemp
 15   enddo
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END
c--------------------------------------------------------------------------------------------
      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL chebev,a,b,x,c(m)
c    Chebyshev evaluation: All arguments are input. c(1:m) is an array of Chebyshev coeﬃ-
c    cients, the ﬁrst m elements of c output from chebft (which must have been called with
c    the same a and b). The Chebyshev polynomial m ck Tk−1(y) − c1 /2 is evaluated at a
c                                                     k=1
c    point y = [x − (b + a)/2]/[(b − a)/2], and the result is returned as the function value.
      INTEGER j
      REAL d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)        !Change of variable.
      y2=2.*y
      do 11 j=m,2,-1            !Clenshaw’s recurrence.
         sv=d
         d=y2*d-dd+c(j)
         dd=sv
 11   enddo
      chebev=y*d-dd+0.5*c(1) !Last step is diﬀerent.
      return
      END
c----------------------------------------------------------------
      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
      INTEGER MAXIT
      REAL ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-16,FPMIN=1.e-30,MAXIT=10000,XMIN=2.,
     $     PI=3.141592653589793d0)
C USES beschb
c     Returns the modiﬁed Bessel functions ri = Iν , rk = Kν and their derivatives rip = Iν ,
c     rkp = Kν , for positive x and for xnu = ν ≥ 0. The relative accuracy is within one or
c     two signiﬁcant digits of EPS. FPMIN is a number close to the machine’s smallest ﬂoating-
c     point number. All internal arithmetic is in double precision. To convert the entire routine
c     to double precision, change the REAL declaration above and decrease EPS to 10−16 . Also
c     convert the subroutine beschb.
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,
     $     fact2,ff,gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,
     $     qnew,ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,
     $     rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessik'
      nl=int(xnu+.5d0)          !nl is the number of downward recurrences
      xmu=xnu-nl                !    of the I’s and upward recurrences
      xmu2=xmu*xmu              !   of K’s. xmu lies between −1/2 and
      xi=1.d0/x                 !  1/2.
      xi2=2.d0*xi
      h=xnu*xi                  !Evaluate CF1 by modiﬁed Lentz’s method
      if(h.lt.FPMIN)h=FPMIN     !    (§5.2).
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
         b=b+xi2
         d=1.d0/(b+d)           !Denominators cannot be zero here, so no
         c=b+1.d0/c             !    need for special precautions.
         del=c*d
         h=del*h
         if(abs(del-1.d0).lt.EPS)goto 1
 11   enddo
      pause 'x too large in bessik; try asymptotic expansion'
 1    continue
      ril=FPMIN                                    !Initialize Iν and Iν for downward recur-
      ripl=h*ril                ! rence.
      ril1=ril                  !Store values for later rescaling.
      rip1=ripl
      fact=xnu*xi
      do 12 l=nl,1,-1
         ritemp=fact*ril+ripl
         fact=fact-xi
         ripl=fact*ritemp+ril
         ril=ritemp
 12   enddo 
      f=ripl/ril                                  ! Now have unnormalized Iµ and Iµ .
      if(x.lt.XMIN) then        !Use series.
         x2=.5d0*x
         pimu=PI*xmu
         if(abs(pimu).lt.EPS)then
            fact=1.d0
         else
            fact=pimu/sin(pimu)
         endif
         d=-log(x2)
         e=xmu*d
         if(abs(e).lt.EPS)then
            fact2=1.d0
         else
            fact2=sinh(e)/e
         endif
         call beschb(xmu,gam1,gam2,gampl,gammi)  !Chebyshev evaluation of Γ1 and Γ2 .
         ff=fact*(gam1*cosh(e)+gam2*fact2*d) !f0 .
         sum=ff
         e=exp(e)
         p=0.5d0*e/gampl        !p0 .
         q=0.5d0/(e*gammi)      !q0 .
         c=1.d0
         d=x2*x2
         sum1=p
         do 13 i=1,MAXIT
            ff=(i*ff+p+q)/(i*i-xmu2)
            c=c*d/i
            p=p/(i-xmu)
            q=q/(i+xmu)
            del=c*ff
            sum=sum+del
            del1=c*(p-i*ff)
            sum1=sum1+del1
            if(abs(del).lt.abs(sum)*EPS)goto 2
 13      enddo
         pause 'bessk series failed to converge'
 2       continue
         rkmu=sum
         rk1=sum1*xi2
      else                                          !Evaluate CF2 by Steed’s algorithm (§5.2),
         b=2.d0*(1.d0+x)        !which is OK because there can be no
         d=1.d0/b               !zero denominators.
         delh=d
         h=delh
         q1=0.d0                                  !Initializations for recurrence (6.7.35).
         q2=1.d0
         a1=.25d0-xmu2
         c=a1
         q=c                    !First term in equation (6.7.34).
         a=-a1
         s=1.d0+q*delh
         do 14 i=2,MAXIT
            a=a-2*(i-1)
            c=-a*c/i
            qnew=(q1-b*q2)/a
            q1=q2
            q2=qnew
            q=q+c*qnew
            b=b+2.d0
            d=1.d0/(b+a*d)
            delh=(b*d-1.d0)*delh
            h=h+delh
            dels=q*delh
            s=s+dels
            if(abs(dels/s).lt.EPS)goto 3        !Need only test convergence of sum since
 14      enddo                  !    CF2 itself converges more quickly.
         pause 'bessik: failure to converge in cf2'
 3       continue
         h=a1*h
         rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s         !Omit the factor exp(−x) to scale all the
       rk1=rkmu*(xmu+x+.5d0-h)*xi                  !  returned functions by exp(x) for x ≥
      endif                     ! XMIN.
      rkmup=xmu*xi*rkmu-rk1
      rimu=xi/(f*rkmu-rkmup)    !Get Iµ from Wronskian.
      ri=(rimu*ril1)/ril        !Scale original Iν and Iν .
      rip=(rimu*rip1)/ril
      do 15 i=1,nl              !Upward recurrence of Kν .
         rktemp=(xmu+i)*xi2*rk1+rkmu
         rkmu=rk1
         rk1=rktemp
 15   enddo
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1
      return
      END
