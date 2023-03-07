!  Waves  3D
!
! Variational form + matrix composition
!
! Hermite cubics
!
 module matrix

 use erf
 use math
 
 implicit none


 complex(dp) CA(Njk,73),     &  ! theta&phi integration element and gauss point & radial (one element)
             CF(Njk,73),     &  ! theta&phi integration element and gauss point & radial (one element)
             CAF(NoVarTerms) ! coefficients at the variational form terms

 complex(dp) Acoeff_k(Njk,NGaussR,(DIM+1)*(DIM+1))


 real(dp) FTransformError



  real(dp) loopcount

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 contains
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!-------------------------------------------
 ! Returns the global index (in KK matrix) for every unknown
 ! s: 1..Nmesh,  m: mAF(1..mAFTot),  n: nAF(1..nAFTot),  k: 1..NDIM
 ! e - radial element, a - local index (1 or 2 for "hat" basic functions)

! this function is commented because a static variable (lookup table) is
! used insted -- for the sake of vectorization

! function GlobalIndex(e,a,n,m,k)
! integer e,a,m,n,k, r,  GlobalIndex
!
!  r = LMat(a,e);
!  GlobalIndex = (r-1)*nAFTot*mAFTot*NDIM + (n-1)*mAFTot*NDIM + (m-1)*NDIM + k
!
! return
! end function
!-------------------------------------------
! ! Opposite to the GlobalIndex.  r - integer radial grid point (A) number.
! subroutine GetIndex(rnmk, r,n,m,k)
! integer rnmk, r,m,n,k, nmk,mk
!
!  r    = (rnmk-1) / (nAFTot*mAFTot*NDIM) + 1
!  nmk  = mod( rnmk-1, nAFTot*mAFTot*NDIM ) + 1
!  n    = (nmk-1)/ (mAFTot*NDIM) + 1
!  mk   = mod( nmk-1, mAFTot*NDIM ) + 1
!  m    = (mk-1)/ (NDIM) + 1
!  k    = mod( mk-1, NDIM ) + 1
! 
! return
! end subroutine
!-------------------------------------------
subroutine Init_dPsiR(e)
  integer e,gR, ij,a,b,da,df, AModified, FModified, iAF
  
  do a = 1, Nbf
  do b = 1, Nbf
  do gR = 1, NGaussR
  do AModified = 0, 1
  do FModified = 0, 1
  do ij = 1, (DIM+1)*(DIM+1)
    da = da_ij_lkp(ij)
    df = df_ij_lkp(ij)

!    dPsiRAF(gR,ij,a,b,0) = Psi_FE(da,a,e, e,gR, 0)*Psi_FE(df,b,e, e,gR, 0)
!    dPsiRAF(gR,ij,a,b,1) = Psi_FE(da,a,e, e,gR, 0)*Psi_FE(df,b,e, e,gR, 1)
!    dPsiRAF(gR,ij,a,b,2) = Psi_FE(da,a,e, e,gR, 1)*Psi_FE(df,b,e, e,gR, 0)
!    dPsiRAF(gR,ij,a,b,3) = Psi_FE(da,a,e, e,gR, 1)*Psi_FE(df,b,e, e,gR, 1)

    iAF = AModified*2 + FModified
    dPsiRAF(gR,iAF,ij,b,a) = Psi_FE(da,a,e, e,gR, AModified)*Psi_FE(df,b,e, e,gR, FModified)
  enddo

  enddo
  enddo
  enddo
  enddo
  enddo
          ! S & M -- standart & modified (the first FE is modified for some poloidal
          ! modes to be consistent with the asymptotics on the axis)

end subroutine
!-------------------------------------------
subroutine InitialWaveVector
  integer itpr,itp,eR,gR

  k_par_it = 0

  do eR = 1, NelR

    if (CylMetric) then
      call Init_Cyl_Equilibruim(eR)
    else
      call Init_Equilibruim(eR)
    end if

    do itpr = 1, Njkr
      itp = itp_itpr_lkp(itpr)
      gR  = gR_itpr_lkp(itpr)

      select case (kparApprox)

       case(1) ! 1 -  k|| = n/R
        k_par_it(:,eR,itpr) = abs(mftp(gR)*max(abs(real(n_ant)),0.01) &
          /(bjac(gR,itp)*B02*bmod(gR,itp)))
       case(2) ! 2 -  k|| = 1/(2qR)
        k_par_it(:,eR,itpr) = 1/(2*mftp(gR)/mfpp(gR)*cR(gR,itp))
       case(3) ! 3 -  k|| = omega/c_A
        k_par_it(:,eR,itpr) = Om_SI*sqrt(fnp(eR,gR))/(CA_SI*bmod(gR,itp))
      end select

    enddo

  enddo

end subroutine
!-------------------------------------------
subroutine Init_Cyl_Equilibruim(eR) ! Cylindrical equilibrium
 integer eR,gR, itp
 real(dp) rw, r, phi, theta, ampl
 real(dp) M3(3,3), iM3(3,3), detM
 real(dp) TL3(3,3), iTU3(3,3)
 integer  me, ne



  rw = 1.0_dp
  
  ! Attention! Lower metric coefficients are divided by bjac
  ! P.ex.  tgspl = gspl_real / bjac

  mTU  = 0;   mTL  = 0;   bjac  = 0;   bmod = 0;
  mTUs = 0;   mTLs = 0;   bjacs = 0;
  mTUt = 0;   mTLt = 0;   bjact = 0;
  mTUp = 0;   mTLp = 0;   bjacp = 0;

  do gR = 1,NGaussR
    do itp = 1, Njk
      r = sqrt(s_gauss(eR,gR))*rw
!      if (r<0.1) r = 0.1
      
!----  s = r^2   ------------
!      mTU(1,1,gR,itp) = 2.*r/(rw*rw)
!      mTU(2,2,gR,itp) = 1./r
!      mTU(3,3,gR,itp) = -1./RRcyl
!
!      mTL(1,1,gR,itp) = (rw*rw)/(2.*r)
!      mTL(2,2,gR,itp) = r
!      mTL(3,3,gR,itp) = -RRcyl
!
!
!      mTUs(1,1,gR,itp) = 1./r
!      mTUs(2,2,gR,itp) = -rw*rw/(2.*r*r*r)
!
!      mTLs(1,1,gR,itp) = -(rw*rw*rw*rw)/(4.*r*r*r)
!      mTLs(2,2,gR,itp) = rw*rw/(2.*r)
!
!
!      bjac(gR,itp)    = -RRcyl/2.*(rw*rw)
!      bmod(gR,itp)    = 1.0


!----  s = r^2,  Normalised metric  ------------
!      mTU(1,1,gR,itp) = 2._dp/(rw*rw)
!      mTU(2,2,gR,itp) = 1._dp
!      mTU(3,3,gR,itp) = -1._dp/RRcyl
!
!      mTL(1,1,gR,itp) = (rw*rw)/(2._dp)
!      mTL(2,2,gR,itp) = 1._dp
!      mTL(3,3,gR,itp) = -RRcyl
!
!
!      bjac(gR,itp)    = -RRcyl/2._dp*(rw*rw)
!      bmod(gR,itp)    = 1.0_dp



!----  s = r   ------------
!      r = s_gauss(eR,gR)
!
!      mTU(1,1,gR,itp) = 1
!      mTU(2,2,gR,itp) = 1._dp/r
!      mTU(3,3,gR,itp) = -1._dp/RRcyl
!
!      mTL(1,1,gR,itp) = 1
!      mTL(2,2,gR,itp) = r
!      mTL(3,3,gR,itp) = -RRcyl
!
!
!      mTUs(1,1,gR,itp) = 0
!      mTUs(2,2,gR,itp) = -1._dp/(r*r)
!
!      mTLs(1,1,gR,itp) = 0
!      mTLs(2,2,gR,itp) = 1
!
!
!      bjac(gR,itp)    = -RRcyl*r
!      bjacs(gR,itp)   = -RRcyl
!      bmod(gR,itp)    = 1.0_dp


      phi   = ph_itp(itp)
      theta = th_itp(itp)


!      mTU(1,1,gR,itp) = 0.9_dp! + 0.3_dp*sin(theta)
!
!      mTU(2,1,gR,itp) = 0.2_dp*(1._dp + 0.8_dp*sin(0*phi))
!      mTU(2,2,gR,itp) = 1.1_dp*(1._dp + 0.8_dp*sin(0*phi))
!      mTU(2,3,gR,itp) = 0.1_dp*(1._dp + 0.8_dp*sin(0*phi))
!
!      mTU(3,1,gR,itp) = 0.32_dp*(1._dp + 0.4_dp*sin(0*theta))
!      mTU(3,2,gR,itp) = 0.24_dp*(1._dp + 0.4_dp*sin(0*theta))
!      mTU(3,3,gR,itp) = 1.2_dp*(1._dp + 0.4_dp*sin(0*theta))




!      mTU(1,1,gR,itp) = 1.1_dp! + 0.3_dp*sin(theta)
!
!!      mTU(2,1,gR,itp) = 1.3_dp*(1._dp + 0.8_dp*sin(0*phi))
!      mTU(2,2,gR,itp) = 1.3_dp
!!      mTU(2,3,gR,itp) = 1.78_dp*(1._dp + 0.8_dp*sin(0*phi))
!
!!      mTU(3,1,gR,itp) = 0.32_dp*(1._dp + 0.4_dp*cos(0*theta))
!!      mTU(3,2,gR,itp) = 0.29_dp*(1._dp + 0.4_dp*cos(0*theta))
!      mTU(3,3,gR,itp) = 1.7_dp

!      mTU(1,1,gR,itp) = 1.1_dp
!      mTU(2,2,gR,itp) = -0.8_dp
!      mTU(3,3,gR,itp) = 1.13_dp

!      me = 1
!      ne = 0
!      ampl = 0.00_dp
!
!      mTU(1,1,gR,itp) = 1.1_dp  * (1._dp + ampl*sin(me*theta + ne*phi)*s_gauss(eR,gR))
!                                                       
!      mTU(2,1,gR,itp) = 0.2_dp  * (1._dp + ampl*sin(me*theta + ne*phi)*s_gauss(eR,gR))*s_gauss(eR,gR)
!      mTU(2,2,gR,itp) = -0.8_dp * (1._dp + ampl*sin(me*theta + ne*phi)*s_gauss(eR,gR))
!      mTU(2,3,gR,itp) = 0.4_dp  * (1._dp + ampl*sin(me*theta + ne*phi)*s_gauss(eR,gR))*s_gauss(eR,gR)
!                                                 
!      mTU(3,1,gR,itp) = 0.1_dp  * (1._dp + ampl*sin(me*theta + ne*phi)*s_gauss(eR,gR))*s_gauss(eR,gR)
!      mTU(3,2,gR,itp) = 0.4_dp  * (1._dp + ampl*sin(me*theta + ne*phi)*s_gauss(eR,gR))*s_gauss(eR,gR)
!      mTU(3,3,gR,itp) = 1.13_dp * (1._dp + ampl*sin(me*theta + ne*phi)*s_gauss(eR,gR))
!
!
!
!
!      M3(:,:) = mTU(:,:,gR,itp)
!      call InverseMatrix3(M3,iM3)
!      iTU3 = iM3
!      call TransposeMatrix3(iM3,M3)
!
!
!      mTL(:,:,gR,itp) = M3(:,:)
!
!      TL3 = mTL(:,:,gR,itp)
!
!      call DetMultMatrix3(TL3,iTU3,detM)
!
!      bjac(gR,itp)    =  - sqrt(detM)
!      bmod(gR,itp)    = 1.0




     
!      theta = 2*pi*(it-1)/nj
!      trx(ip,it,is)   = RRcyl + r*cos(theta)
!      trz(ip,it,is)   = r*sin(theta)
      
    enddo
  enddo


  call Init_Equilibruim(eR)

end subroutine
!-------------------------------------------
subroutine Init_Equilibruim(eR) ! Interpolate the equilibrium (TU,TL,B,etc) on the gauss radial
                               ! grid for the radial element "e"
 integer is,ism,isp, eR,gR, it,ip,itp, ii,jj
 real(dp) rw, RR, s, r
 integer  it_m,it_p, ip_m,ip_p
 
 real(dp) M3(3,3), iM3(3,3), detM
 real(dp) TL3(3,3), iTU3(3,3)

  do gR = 1,NGaussR

    ! ism,isp - radial surfaces on the TERPSICHORE grid, neighbour to the s_gauss(eR,gR) surface
    ism = min(niT-1, max(1, int(s_gauss(eR,gR)*niT + 0.5_dp)))
    isp = ism + 1

    mftp(gR) = (ftp(ism)*(s_T(isp)-s_gauss(eR,gR)) &
              + ftp(isp)*(s_gauss(eR,gR)-s_T(ism)))/(s_T(isp)-s_T(ism))
    mfpp(gR) = (fpp(ism)*(s_T(isp)-s_gauss(eR,gR)) &
              + fpp(isp)*(s_gauss(eR,gR)-s_T(ism)))/(s_T(isp)-s_T(ism))

    mci(gR) = (ci(ism)*(s_T(isp)-s_gauss(eR,gR)) &
             + ci(isp)*(s_gauss(eR,gR)-s_T(ism)))/(s_T(isp)-s_T(ism))
    mcj(gR) = (cj(ism)*(s_T(isp)-s_gauss(eR,gR)) &
             + cj(isp)*(s_gauss(eR,gR)-s_T(ism)))/(s_T(isp)-s_T(ism))


    call InterpEqSurf(ism)

    do is = 1, 2
      tTU(1,1,:,is) = tTU(1,1,:,is) / setTU(1,is+ism-1)
      tTU(2,2,:,is) = tTU(2,2,:,is) / setTU(2,is+ism-1)
      tTL(1,1,:,is) = tTL(1,1,:,is) / setTL(1,is+ism-1)
      tTL(2,2,:,is) = tTL(2,2,:,is) / setTL(2,is+ism-1)
    enddo

  do itp = 1, Njk  ! 1..N <-- used for poloidal derivative calculation
    ! Linear interpolation radially
    do ii=1,DIM;  do jj=1,DIM
    mTU(ii,jj,gR,itp) = ( tTU(ii,jj,itp,1)*(s_T(isp)-s_gauss(eR,gR)) + &
                          tTU(ii,jj,itp,2)*(s_gauss(eR,gR)-s_T(ism))      )/(s_T(isp)-s_T(ism))
!    mTL(ii,jj,gR,itp) = ( tTL(ii,jj,itp,1)*(s_T(isp)-s_gauss(eR,gR)) + &
!                          tTL(ii,jj,itp,2)*(s_gauss(eR,gR)-s_T(ism))      )/(s_T(isp)-s_T(ism))
    enddo; enddo


    M3(:,:) = mTU(:,:,gR,itp)
    M3(1,1) = M3(1,1) * seTU(1,eR,gR)
    M3(2,2) = M3(2,2) * seTU(2,eR,gR)
    call InverseMatrix3(M3,iM3)
    iTU3 = iM3
    call TransposeMatrix3(iM3,M3)

    mTL(:,:,gR,itp) = M3(:,:)
    TL3 = mTL(:,:,gR,itp)

    call DetMultMatrix3(TL3,iTU3,detM)

    bjac(gR,itp)    = - sqrt(detM)
    mTL(1,1,gR,itp) = mTL(1,1,gR,itp) / seTL(1,eR,gR)
    mTL(2,2,gR,itp) = mTL(2,2,gR,itp) / seTL(2,eR,gR)

!    bjac(gR,itp) = ( tbjac(itp,1)*(s_T(isp)-s_gauss(eR,gR)) + &
!                     tbjac(itp,2)*(s_gauss(eR,gR)-s_T(ism))      )/(s_T(isp)-s_T(ism))
    bmod(gR,itp) = ( tbmod(itp,1)*(s_T(isp)-s_gauss(eR,gR)) + &
                     tbmod(itp,2)*(s_gauss(eR,gR)-s_T(ism))      )/(s_T(isp)-s_T(ism))

     mgssu(gR,itp) = (tgssu(itp,1)*(s_T(isp)-s_gauss(eR,gR)) + &
                     tgssu(itp,2)*(s_gauss(eR,gR)-s_T(ism))      )/(s_T(isp)-s_T(ism))

    ip = (itp-1)/njT + 1
    cR(gR,itp) = (( (trx(itp,1)-MAxisR(ip))/sqrt(s_T(ism))*(s_T(isp)-s_gauss(eR,gR)) + &
  (trx(itp,2)-MAxisR(ip))/sqrt(s_T(isp))*(s_gauss(eR,gR)-s_T(ism)))/(s_T(isp)-s_T(ism)))* &
          sqrt(s_gauss(eR,gR)) + MAxisR(ip)
    cZ(gR,itp) = (( (trz(itp,1)-MAxisZ(ip))/sqrt(s_T(ism))*(s_T(isp)-s_gauss(eR,gR)) + &
  (trz(itp,2)-MAxisZ(ip))/sqrt(s_T(isp))*(s_gauss(eR,gR)-s_T(ism)))/(s_T(isp)-s_T(ism)))* &
          sqrt(s_gauss(eR,gR)) + MAxisZ(ip)
    phiv(gR,itp) = ( tphiv(itp,1)*(s_T(isp)-s_gauss(eR,gR)) + &
                     tphiv(itp,2)*(s_gauss(eR,gR)-s_T(ism))      )/(s_T(isp)-s_T(ism))


    ! Calculating radial derivatives
    ! linear interpolation - derivatives supposed to be step-wise constant
    do ii=1,DIM;  do jj=1,DIM
    mTLs(ii,jj,gR,itp) = (tTL(ii,jj,itp,2)-tTL(ii,jj,itp,1))/(s_T(isp)-s_T(ism))
    mTUs(ii,jj,gR,itp) = (tTU(ii,jj,itp,2)-tTU(ii,jj,itp,1))/(s_T(isp)-s_T(ism))
    enddo; enddo
    bjacs(gR,itp) = (tbjac(itp,2)-tbjac(itp,1))/(s_T(isp)-s_T(ism))
  enddo ! itp

  enddo ! gR


  do gR = 1,NGaussR
    do itp = 1, Njk

    it = mod(itp-1,njT) + 1
    ip = (itp-1)/njT + 1


    ! Calculating poloidal derivatives from the newly interpolated values

    it_m = it - 1
    if (it_m==0)     it_m = njT
    it_p = it + 1
    if (it_p==njT+1) it_p = 1

    do ii=1,DIM;  do jj=1,DIM
    mTLt(ii,jj,gR,itp) = (mTL(ii,jj,gR,itp_it_ip_lkp(it_p,ip)) - mTL(ii,jj,gR,itp_it_ip_lkp(it_m,ip))) &
                        /(theta_T(it+1)-theta_T(it-1))
    mTUt(ii,jj,gR,itp) = (mTU(ii,jj,gR,itp_it_ip_lkp(it_p,ip)) - mTU(ii,jj,gR,itp_it_ip_lkp(it_m,ip))) &
                        /(theta_T(it+1)-theta_T(it-1))
    enddo; enddo
    bjact(gR,itp) = (bjac(gR,itp_it_ip_lkp(it_p,ip)) - bjac(gR,itp_it_ip_lkp(it_m,ip)))/(theta_T(it+1)-theta_T(it-1))


    ! Calculating toroidal derivatives from the newly interpolated values
    ip_m = ip - 1
    if (ip_m==0)     ip_m = nkT
    ip_p = ip + 1
    if (ip_p==nkT+1) ip_p = 1

    do ii=1,DIM;  do jj=1,DIM
    mTLp(ii,jj,gR,itp) = (mTL(ii,jj,gR,itp_it_ip_lkp(it,ip_p)) - mTL(ii,jj,gR,itp_it_ip_lkp(it,ip_m))) &
                        /(phi_T(ip+1)-phi_T(ip-1))
    mTUp(ii,jj,gR,itp) = (mTU(ii,jj,gR,itp_it_ip_lkp(it,ip_p)) - mTU(ii,jj,gR,itp_it_ip_lkp(it,ip_m))) &
                        /(phi_T(ip+1)-phi_T(ip-1))
    enddo; enddo
    bjacp(gR,itp) = (bjac(gR,itp_it_ip_lkp(it,ip_p)) - bjac(gR,itp_it_ip_lkp(it,ip_m)))/(phi_T(ip+1)-phi_T(ip-1))

    enddo ! itp
  enddo ! gR


! contains

!  ! 1D linear interpolation: f(x) calculated from a,b,f(a),f(b),x
!  function LinI(a,b,fa,fb,x)
!    real(dp) a,b,fa,fb,x, LinI
!    LinI = (fa*(b-x) + fb*(x-a))/(b-a)
!    return
!  end function


!  ! Poloidal derivative + radial interpolation
!  ! a,b - radial surfaces (x - point of interpolation - lies between a and b)
!  ! tm, tp - theta surfaces. The derivative is calculated at (tm+tp)/2
!  ! fatm,fatp,fbtm,fbtp - function values at (a,tm), (a,tp), (b,tm), (b,tp)
!  function DiffTheta(a,b,tm,tp,fatm,fatp,fbtm,fbtp,x)
!    real(dp) a,b,tm,tp,fatm,fatp,fbtm,fbtp,x, DiffTheta,   dfa, dfb
!    dfa = (fatp - fatm)/(tp-tm)
!    dfb = (fatp - fatm)/(tp-tm)
!    DiffTheta = LinI(a,b,dfa,dfb,x)
!    return
!  end function



end subroutine
!-------------------------------------------
!subroutine Init_Equilibruim(eR) ! Interpolate the equilibrium (TU,TL,B,etc) on the gauss radial
!                                ! grid for the radial element "e"
! integer  eR,gR, it,ip,itp, ii,jj
! integer  it_m,it_p, ip_m,ip_p
! 
! real(dp) vint_tmp, vint_tmpy1
! real(dp) M3(3,3), iM3(3,3), detM
! real(dp) TL3(3,3), iTU3(3,3)
!
!
!
!  do gR = 1,NGaussR
!  do itp = 1, Njk  ! 1..N <-- used for poloidal derivative calculation
!
!    ! Cubic splines radially
!    do ii=1,DIM;  do jj=1,DIM
!      vect_tmp(:)   = tTU(ii,jj,itp,:)
!      vect_tmpy2(:) = tTUy2(ii,jj,itp,:)
!      call SplInt(s_T,vect_tmp,vect_tmpy2,niT,s_gauss(eR,gR),vint_tmp,vint_tmpy1)
!      mTU(ii,jj,gR,itp)  = vint_tmp
!      mTUs(ii,jj,gR,itp) = vint_tmpy1
!   
!      vect_tmp(:)   = tTL(ii,jj,itp,:)
!      vect_tmpy2(:) = tTLy2(ii,jj,itp,:)
!      call SplInt(s_T,vect_tmp,vect_tmpy2,niT,s_gauss(eR,gR),vint_tmp,vint_tmpy1)
!!      mTL(ii,jj,gR,itp)  = vint_tmp
!      mTLs(ii,jj,gR,itp) = vint_tmpy1
!    enddo;  enddo
!
!
!!    print *, s_gauss(eR,gR), mTU(1,1,gR,itp), 2._dp*sqrt(s_gauss(eR,gR))
!!    print *, s_gauss(eR,gR), mTU(2,2,gR,itp), 1._dp/sqrt(s_gauss(eR,gR))
!!    print *, s_gauss(eR,gR), mTU(3,3,gR,itp), -1._dp
!!    stop
!
!!    mTU(1,1,gR,itp) = 2._dp*sqrt(s_gauss(eR,gR))
!!    mTU(2,2,gR,itp) = 1._dp/sqrt(s_gauss(eR,gR))
!!    mTU(3,3,gR,itp) = -1._dp
!
!!    mTUs(1,1,gR,itp) = 1._dp/sqrt(s_gauss(eR,gR))
!!    mTUs(2,2,gR,itp) = -0.5_dp/(sqrt(s_gauss(eR,gR))*s_gauss(eR,gR))
!!    mTUs(3,3,gR,itp) = 0
!
!!    print *, 'Diff: '
!!    print *, mTLs(1,1,gR,itp), (-0.25_dp/(sqrt(s_gauss(eR,gR))*s_gauss(eR,gR)))
!!    mTLs(1,1,gR,itp) = -0.25_dp/(sqrt(s_gauss(eR,gR))*s_gauss(eR,gR))
!!    mTLs(2,2,gR,itp) = 0.5_dp/sqrt(s_gauss(eR,gR))
!!    mTLs(3,3,gR,itp) = 0
!
!    M3(:,:) = mTU(:,:,gR,itp)
!    call InverseMatrix3(M3,iM3)
!    iTU3 = iM3
!    call TransposeMatrix3(iM3,M3)
!    mTL(:,:,gR,itp) = M3(:,:)
!    TL3 = mTL(:,:,gR,itp)
!    call DetMultMatrix3(TL3,iTU3,detM)
!    bjac(gR,itp)    = - sqrt(detM)
!
!    vect_tmp(:)   = tbjac(itp,:)
!    vect_tmpy2(:) = tbjacy2(itp,:)
!    call SplInt(s_T,vect_tmp,vect_tmpy2,niT,s_gauss(eR,gR),vint_tmp,vint_tmpy1)
!    bjacs(gR,itp) = vint_tmpy1
!
!    vect_tmp(:)   = tbmod(itp,:)
!    vect_tmpy2(:) = tbmody2(itp,:)
!    call SplInt(s_T,vect_tmp,vect_tmpy2,niT,s_gauss(eR,gR),vint_tmp,vint_tmpy1)
!    bmod(gR,itp)  = vint_tmp
!
!
!!    bjacs(gR,itp) = 0
!  enddo ! itp
!  enddo ! gR
!
!
!  do gR = 1,NGaussR
!    do itp = 1, Njk
!
!    it = mod(itp-1,njT) + 1
!    ip = (itp-1)/njT + 1
!
!
!    ! Calculating poloidal derivatives from the newly interpolated values
!    do ii=1,DIM;  do jj=1,DIM
!
!    it_m = it - 1
!    if (it_m==0)     it_m = njT
!    it_p = it + 1
!    if (it_p==njT+1) it_p = 1
!
!
!    mTLt(ii,jj,gR,itp) = (mTL(ii,jj,gR,itp_it_ip_lkp(it_p,ip)) - mTL(ii,jj,gR,itp_it_ip_lkp(it_m,ip))) &
!                        /(theta_T(it+1)-theta_T(it-1))
!    mTUt(ii,jj,gR,itp) = (mTU(ii,jj,gR,itp_it_ip_lkp(it_p,ip)) - mTU(ii,jj,gR,itp_it_ip_lkp(it_m,ip))) &
!                        /(theta_T(it+1)-theta_T(it-1))
!    enddo; enddo
!    bjact(gR,itp) = (bjac(gR,itp_it_ip_lkp(it_p,ip)) - bjac(gR,itp_it_ip_lkp(it_m,ip)))/(theta_T(it+1)-theta_T(it-1))
!
!
!    ! Calculating toroidal derivatives from the newly interpolated values
!    do ii=1,DIM;  do jj=1,DIM
!
!    ip_m = ip - 1
!    if (ip_m==0)     ip_m = nkT
!    ip_p = ip + 1
!    if (ip_p==nkT+1) ip_p = 1
!
!
!    mTLp(ii,jj,gR,itp) = (mTL(ii,jj,gR,itp_it_ip_lkp(it,ip_p)) - mTL(ii,jj,gR,itp_it_ip_lkp(it,ip_m))) &
!                        /(phi_T(ip+1)-phi_T(ip-1))
!    mTUp(ii,jj,gR,itp) = (mTU(ii,jj,gR,itp_it_ip_lkp(it,ip_p)) - mTU(ii,jj,gR,itp_it_ip_lkp(it,ip_m))) &
!                        /(phi_T(ip+1)-phi_T(ip-1))
!    enddo; enddo
!    bjacp(gR,itp) = (bjac(gR,itp_it_ip_lkp(it,ip_p)) - bjac(gR,itp_it_ip_lkp(it,ip_m)))/(phi_T(ip+1)-phi_T(ip-1))
!
!    enddo ! itp
!  enddo ! gR
!
!
!! contains
!
!!  ! 1D linear interpolation: f(x) calculated from a,b,f(a),f(b),x
!!  function LinI(a,b,fa,fb,x)
!!    real(dp) a,b,fa,fb,x, LinI
!!    LinI = (fa*(b-x) + fb*(x-a))/(b-a)
!!    return
!!  end function
!
!
!!  ! Poloidal derivative + radial interpolation
!!  ! a,b - radial surfaces (x - point of interpolation - lies between a and b)
!!  ! tm, tp - theta surfaces. The derivative is calculated at (tm+tp)/2
!!  ! fatm,fatp,fbtm,fbtp - function values at (a,tm), (a,tp), (b,tm), (b,tp)
!!  function DiffTheta(a,b,tm,tp,fatm,fatp,fbtm,fbtp,x)
!!    real(dp) a,b,tm,tp,fatm,fatp,fbtm,fbtp,x, DiffTheta,   dfa, dfb
!!    dfa = (fatp - fatm)/(tp-tm)
!!    dfb = (fatp - fatm)/(tp-tm)
!!    DiffTheta = LinI(a,b,dfa,dfb,x)
!!    return
!!  end function
!
!
!
!end subroutine
!-------------------------------------------

! Calculate variational form coefficients on the "e" radial element, for given frequency
subroutine Init_a_eq(eR,gR,iR)

  integer eR,gR,iR, ii,jj,ij,  &
          it,ip,         & ! poloidal & toroidal index
          itp,           & ! combined toroidal & poloidal index
          itpr,          & ! combined toroidal & poloidal & radial gauss point (not element!) index
          k,kp, a, ni, nj, Ncf, Nsb, i, j, ij_i_j_lkp(0:DIM,0:DIM)

  complex(dp) K0,om,A1,A2,t11,t12,t13,t21,t22,t23,t31,t32,t33, &
              OmCIe,OmCIi1,OmCIi2,OmCIih

  complex(dp) eS(4),eD(4),eP(4)
 
  real(dp) s, TL(DIM,DIM),TLs(DIM,DIM),TLt(DIM,DIM),TLp(DIM,DIM), &
              TU(DIM,DIM),TUs(DIM,DIM),TUt(DIM,DIM),TUp(DIM,DIM), &
              jac, jacs, jact, jacp,                              &
              sf1, sf2, sf3, sf4, dsf1, dsf2, dsf3, dsf4,         &
              seTU1, seTU2, seTL1, seTL2, dseTU1, dseTU2, dseTL1, dseTL2, &
              OmPlaesqr, OmPlai1sqr, OmPlai2sqr, OmPlaihsqr
  complex(dp) Z_1i1,Z_1i2,Z_1ih,Z_1ip,Z_1e,Z0i1,Z0i2,Z0ih,Z0ip,Z0e,Z1i1,Z1i2,Z1ih,Z1ip,Z1e,cw_arg
  real(dp) k_par(8), v_the, v_thi1, v_thi2, v_thih, CF_Anis_A, CF_Anis_B, CF_Anis_C

 do i=0,DIM
   do j=0,DIM
     ij_i_j_lkp(i,j)=i*(DIM+1)+j+1
   enddo
 enddo


 CA = 0;  CF = 0
 Acoeff = 0

 K0 = Om_SI/c_SI
 om = Om_SI/OmCi1_SI


!!! Testing
 CAF(1) = -1
 CAF(2) = -1
 CAF(3) = K0*K0
 CAF(4) = iim*K0
 CAF(5) = 1
 CAF(6) = -iim*K0
 CAF(7) = 0
 CAF(8) = 0

! Thermal velocities
 v_the =sqrt(2*Te0*fte(eR,gR)*1.6e-19/me_SI)
 v_thi1=sqrt(2*Ti0*fti(eR,gR)*1.6e-19/mi1_SI)
 v_thi2=sqrt(2*Ti0*ftm(eR,gR)*1.6e-19/mi2_SI)
 v_thih=sqrt(2*Th0*fth(eR,gR)*1.6e-19/mih_SI)


 do itp = 1, Njk
 
   s  = s_gauss(eR,gR)
 
   itpr = itpr_itp_gR_lkp(itp,gR)
   it  = it_itpr_lkp(itpr)
   ip  = ip_itpr_lkp(itpr)


! Dielectric tensor
  OmPlaesqr  = 4*pi*(ne_SI*fnp(eR,gR)+nh0_SI*iZh*fnhp(eR,gR))/1.e6_dp*qe_Gauss*qe_Gauss/(me_SI*1.e3_dp)
  OmPlai1sqr = 4*pi*n0_SI*ifr1/1.e6_dp*fnp(eR,gR)*iZ1*iZ1*qe_Gauss*qe_Gauss/(mi1_SI*1.e3_dp)
  OmPlai2sqr = 4*pi*n0_SI*ifr2/1.e6_dp*fnm(eR,gR)*iZ2*iZ2*qe_Gauss*qe_Gauss/(mi2_SI*1.e3_dp)
  OmPlaihsqr = 4*pi*nh0_SI/1.e6_dp*fnhp(eR,gR)*iZh*iZh*qe_Gauss*qe_Gauss/(mih_SI*1.e3_dp)
  OmCIe  = qe_Gauss*B0*1.e4_dp*bmod(gR,itp)/(c_SI*100*me_SI*1.e3_dp)
  OmCIi1 = iZ1*qe_Gauss*B0*1.e4_dp*bmod(gR,itp)/(c_SI*100*mi1_SI*1.e3_dp)
  OmCIi2 = iZ2*qe_Gauss*B0*1.e4_dp*bmod(gR,itp)/(c_SI*100*mi2_SI*1.e3_dp)
  OmCIih = iZh*qe_Gauss*B0*1.e4_dp*bmod(gR,itp)/(c_SI*100*mih_SI*1.e3_dp)


! Approximated k||

  k_par(1) = k_par_it(1,eR,itpr)
  j = 2
  do i = 2, 8
   if (k_par_lkp(i-1)==1) then
    k_par(i) = k_par_it(j,eR,itpr)
    j = j+1
   else
    k_par(i) = k_par_it(1,eR,itpr)
   end if
  enddo


 select case (DielTensor)

 case (1)  ! 1 - full cold plasma dielectric tensor, "LION"

      eS(1) = - OmPlaesqr/((Om_SI*Om_SI - OmCIe*OmCIe))
      eS(2) = - OmPlai1sqr/((Om_SI*Om_SI - OmCIi1*OmCIi1))
      eS(3) = - OmPlai2sqr/((Om_SI*Om_SI - OmCIi2*OmCIi2))
      eS(4) = 0.0_dp

      eD(1) = -iim * OmCIe/Om_SI*OmPlaesqr/ &
                  (Om_SI*Om_SI - OmCIe*OmCIe)
      eD(2) =  iim * OmCIi1/Om_SI*OmPlai1sqr/ &
                  (Om_SI*Om_SI - OmCIi1*OmCIi1)
      eD(3) =  iim *  OmCIi2/Om_SI*OmPlai2sqr/ &
                  (Om_SI*Om_SI - OmCIi2*OmCIi2)
      eD(4) = 0.0_dp

      eP(1) = - OmPlaesqr/(Om_SI*Om_SI)
      eP(2) = - OmPlai1sqr/(Om_SI*Om_SI)
      eP(3) = - OmPlai2sqr/(Om_SI*Om_SI)
      eP(4) = 0.0_dp

      eP(:) = eP(:)*(1.0_dp + iim*om_im*nu_f(eR,gR))
                  

 case (2)  ! 2 - 'Andre' model (imaginary eps_parallel)

      eS(1) = - OmPlaesqr/((Om_SI*Om_SI - OmCIe*OmCIe))
      eS(2) = - OmPlai1sqr/((Om_SI*Om_SI - OmCIi1*OmCIi1))
      eS(3) = - OmPlai2sqr/((Om_SI*Om_SI - OmCIi2*OmCIi2))
      eS(4) = 0.0_dp

      eD(1) = -iim * OmCIe/Om_SI*OmPlaesqr/ &
                  (Om_SI*Om_SI - OmCIe*OmCIe)
      eD(2) =  iim * OmCIi1/Om_SI*OmPlai1sqr/ &
                  (Om_SI*Om_SI - OmCIi1*OmCIi1)
      eD(3) =  iim *  OmCIi2/Om_SI*OmPlai2sqr/ &
                  (Om_SI*Om_SI - OmCIi2*OmCIi2)
      eD(4) = 0.0_dp

      eP(1) = - OmPlaesqr/(Om_SI*iim*om_im)
      eP(2) = - OmPlai1sqr/(Om_SI*Om_SI)
      eP(3) = - OmPlai2sqr/(Om_SI*Om_SI)
      eP(4) = 0.0_dp

 case (3)  ! 3 - vacuum  E = I

      eS(1) = iim*om_im
      eS(2) = 0.0_dp
      eS(3) = 0.0_dp
      eS(4) = 0.0_dp

      eS(1) = 0.0_dp
      eS(2) = 0.0_dp
      eS(3) = 0.0_dp
      eS(4) = 0.0_dp

      eP(1) = iim*om_im
      eP(2) = 0.0_dp
      eP(3) = 0.0_dp
      eP(4) = 0.0_dp

 case (5)  ! 5 - full cold plasma, om_e = om_e*(1 + i*delta)

      eS(1) = - OmPlaesqr/(1.0_dp/(1.0_dp + iim*om_im )*  &
           (Om_SI*Om_SI*(1.0_dp + iim*om_im )*(1.0_dp + iim*om_im ) - OmCIe*OmCIe))
      eS(2) = - OmPlai1sqr/(1.0_dp/(1.0_dp + iim*om_imi)*  &
           (Om_SI*Om_SI*(1.0_dp + iim*om_imi)*(1.0_dp + iim*om_imi) - OmCIi1*OmCIi1))
      eS(3) = - OmPlai2sqr/(1.0_dp/(1.0_dp + iim*om_imi)*  &
           (Om_SI*Om_SI*(1.0_dp + iim*om_imi)*(1.0_dp + iim*om_imi) - OmCIi2*OmCIi2))
      eS(4) = 0.0_dp

      eD(1) = -iim * OmCIe/Om_SI*OmPlaesqr/ &
                  (Om_SI*Om_SI*(1.0_dp + iim*om_im )*(1.0_dp + iim*om_im ) - OmCIe*OmCIe)
      eD(2) =  iim * OmCIi1/Om_SI*OmPlai1sqr/ &
                  (Om_SI*Om_SI*(1.0_dp + iim*om_imi)*(1.0_dp + iim*om_imi) - OmCIi1*OmCIi1)
      eD(3) =  iim * OmCIi2/Om_SI*OmPlai2sqr/ &
                  (Om_SI*Om_SI*(1.0_dp + iim*om_imi)*(1.0_dp + iim*om_imi) - OmCIi2*OmCIi2)
      eD(4) = 0.0_dp

      eP(1) = - OmPlaesqr/(Om_SI*Om_SI*(1.0_dp + iim*om_im ))
      eP(2) = - OmPlai1sqr/(Om_SI*Om_SI*(1.0_dp + iim*om_imi))
      eP(3) = - OmPlai2sqr/(Om_SI*Om_SI*(1.0_dp + iim*om_imi))
      eP(4) = 0.0_dp

 case (6)  ! 6 - warm model

   cw_arg = (Om_SI*cmplx(1.00,om_imi)-real(OmCIi1))/(k_par(6)*v_thi1)
   Z1i1=-iim*(sqrt(pi)*OmPlai1sqr)/(k_par(6)*v_thi1)*cw(cw_arg)
   cw_arg = (Om_SI*cmplx(1.00,om_imi)-real(OmCIi2))/(k_par(7)*v_thi2)
   Z1i2=-iim*(sqrt(pi)*OmPlai2sqr)/(k_par(7)*v_thi2)*cw(cw_arg)
   cw_arg = (Om_SI*cmplx(1.00,om_imi)-real(OmCIih))/(k_par(8)*v_thih)
   Z1ih=-iim*(sqrt(pi)*OmPlaihsqr)/(k_par(8)*v_thih)*cw(cw_arg)
   cw_arg = (Om_SI*cmplx(1.00,om_im)+real(OmCIe))/(k_par(1)*v_the)
   Z1e=-iim*(sqrt(pi)*OmPlaesqr)/(k_par(1)*v_the)*cw(cw_arg)
   cw_arg = (Om_SI*cmplx(1.00,om_imi)+real(OmCIi1))/(k_par(1)*v_thi1)
   Z_1i1=-iim*(sqrt(pi)*OmPlai1sqr)/(k_par(1)*v_thi1)*cw(cw_arg)
   cw_arg = (Om_SI*cmplx(1.00,om_imi)+real(OmCIi2))/(k_par(1)*v_thi2)
   Z_1i2=-iim*(sqrt(pi)*OmPlai2sqr)/(k_par(1)*v_thi2)*cw(cw_arg)
   cw_arg = (Om_SI*cmplx(1.00,om_imi)+real(OmCIih))/(k_par(1)*v_thih)
   Z_1ih=-iim*(sqrt(pi)*OmPlaihsqr)/(k_par(1)*v_thih)*cw(cw_arg)
   cw_arg = (Om_SI*cmplx(1.00,om_im)-real(OmCIe))/(k_par(1)*v_the)
   Z_1e=-iim*(sqrt(pi)*OmPlaesqr)/(k_par(1)*v_the)*cw(cw_arg)
   cw_arg = Om_SI*cmplx(1.00,om_imi)/(k_par(3)*v_thi1)
   Z0i1=-iim*(sqrt(pi)*OmPlai1sqr)/(k_par(3)*v_thi1)*cw(cw_arg)
   cw_arg = Om_SI*cmplx(1.00,om_imi)/(k_par(4)*v_thi2)
   Z0i2=-iim*(sqrt(pi)*OmPlai2sqr)/(k_par(4)*v_thi2)*cw(cw_arg)
   cw_arg = Om_SI*cmplx(1.00,om_imi)/(k_par(5)*v_thih)
   Z0ih=-iim*(sqrt(pi)*OmPlaihsqr)/(k_par(5)*v_thih)*cw(cw_arg)
   cw_arg = Om_SI*cmplx(1.00,om_im)/(k_par(2)*v_the)
   Z0e=-iim*(sqrt(pi)*OmPlaesqr)/(k_par(2)*v_the)*cw(cw_arg)

   eS(1) = - 1.0_dp/(2.0_dp*Om_SI)*(Z1e+Z_1e)
   eS(2) = - 1.0_dp/(2.0_dp*Om_SI)*(Z1i1+Z_1i1)
   eS(3) = - 1.0_dp/(2.0_dp*Om_SI)*(Z1i2+Z_1i2)

   eD(1) = iim/(2.0_dp*Om_SI)*(Z1e-Z_1e)
   eD(2) = iim/(2.0_dp*Om_SI)*(Z1i1-Z_1i1)
   eD(3) = iim/(2.0_dp*Om_SI)*(Z1i2-Z_1i2)


   if(((abs(k_par(2))*v_the)/Om_SI)**2<8.0e-5) then
    eP(1) = - OmPlaesqr/(Om_SI*Om_SI*(1.0_dp + iim*om_im ))
   else
    eP(1) = 2.0_dp/(k_par(2)*v_the)**2*(OmPlaesqr-Om_SI*cmplx(1.00,om_im)*Z0e)*cmplx(1.00,om_im)
   end if

   if(((abs(k_par(3))*v_thi1)/Om_SI)**2<8.0e-5) then
    eP(2) = - OmPlai1sqr/(Om_SI*Om_SI*(1.0_dp + iim*om_imi))
   else
    eP(2) = 2.0_dp/(k_par(3)*v_thi1)**2*(OmPlai1sqr-Om_SI*cmplx(1.00,om_imi)*Z0i1)*cmplx(1.00,om_imi)
   end if

   if(((abs(k_par(4))*v_thi2)/Om_SI)**2<8.0e-5) then
    eP(3) = - OmPlai2sqr/(Om_SI*Om_SI*(1.0_dp + iim*om_imi))
   else
    eP(3) = 2.0_dp/(k_par(4)*v_thi2)**2*(OmPlai2sqr-Om_SI*cmplx(1.00,om_imi)*Z0i2)*cmplx(1.00,om_imi)
   end if

   eS(4) = 0.0_dp
   eD(4) = 0.0_dp
   eP(4) = 0.0_dp


   if (nh0_SI>1.0) then

    Cf_Anis_A = Bc/(B0*bmod(gR,itp)) + anis(eR,gR)*(1-Bc/(B0*bmod(gR,itp)))
    Cf_Anis_B = Bc/(B0*bmod(gR,itp)) + anis(eR,gR)*(Bc/(B0*bmod(gR,itp))-1)

    if (B0*bmod(gR,itp)>Bc) then
      eS(4) = - 1.0_dp/(2.0_dp*Om_SI*Cf_Anis_A)*(Z1ih+Z_1ih)
      eD(4) = iim/(2.0_dp*Om_SI*Cf_Anis_A)*(Z1ih-Z_1ih)
      eP(4) = 2.0_dp/(k_par(5)*v_thih)**2/Cf_Anis_A &
              * (OmPlaihsqr-Om_SI*cmplx(1.00,om_imi)*Z0ih)*cmplx(1.00,om_imi)

    else

      Cf_Anis_C = sqrt((Bc-B0*bmod(gR,itp))/Bc)

      cw_arg = (Om_SI*cmplx(1.00,om_imi)-real(OmCIih))/(Cf_Anis_C*k_par(8)*v_thih*sqrt(anis(eR,gR)))
      Z1ip=-iim*(sqrt(pi)*OmPlaihsqr)/(Cf_Anis_C*k_par(8)*v_thih*sqrt(anis(eR,gR)))*cw(cw_arg)
      cw_arg = (Om_SI*cmplx(1.00,om_imi)+real(OmCIih))/(Cf_Anis_C*k_par(1)*v_thih*sqrt(anis(eR,gR)))
      Z_1ip=-iim*(sqrt(pi)*OmPlaihsqr)/(Cf_Anis_C*k_par(1)*v_thih*sqrt(anis(eR,gR)))*cw(cw_arg)
      cw_arg = Om_SI*cmplx(1.00,om_imi)/(Cf_Anis_C*k_par(5)*v_thih*sqrt(anis(eR,gR)))
      Z0ip=-iim*(sqrt(pi)*OmPlaihsqr)/(Cf_Anis_C*k_par(5)*v_thih*sqrt(anis(eR,gR)))*cw(cw_arg)

      eS(4) = - 1.0_dp/(2.0_dp*Om_SI) * 1/(Cf_Anis_A*CF_Anis_B)*( Cf_Anis_B *(Z1ih+Z_1ih)  &
              + 2.0*anis(eR,gR)**(3./2)*(1-Bc/(B0*bmod(gR,itp)))*Cf_Anis_C*(Z1ip+Z_1ip))

      eD(4) = iim/(2.0_dp*Om_SI) * 1/(Cf_Anis_A*CF_Anis_B)*( Cf_Anis_B*(Z1ih-Z_1ih)    &
              + 2.0*anis(eR,gR)**(3./2)*(1-Bc/(B0*bmod(gR,itp)))*Cf_Anis_C*(Z1ip-Z_1ip))

      eP(4) = 2.0_dp/(k_par(5)*v_thih)**2 * (1/Cf_Anis_A*                   &
               (OmPlaihsqr-Om_SI*cmplx(1.00,om_imi)*Z0ih)*cmplx(1.00,om_imi)    &
              - 2.0 * Bc/(B0*bmod(gR,itp)) * Cf_Anis_C/(Cf_Anis_A*Cf_Anis_B*sqrt(anis(eR,gR))) *  &
               (OmPlaihsqr-Om_SI*cmplx(1.00,om_imi)*Z0ip)*cmplx(1.00,om_imi))

    end if

   end if

   endselect

   t11 = 1.0_dp + eS(1) + eS(2) + eS(3) + eS(4)
   t12 = - eD(1) - eD(2) - eD(3) - eD(4)
   t13 = 0.0_dp
   t21 = eD(1) + eD(2) + eD(3) + eD(4)
   t22 = 1.0_dp + eS(1) + eS(2) + eS(3) + eS(4)
   t23 = 0.0_dp
   t31 = 0.0_dp
   t32 = 0.0_dp
   t33 = 1.0_dp + eP(1) + eP(2) + eP(3) + eP(4)

! Saving the dielectric tensor components for the diagnostics 
! (used for absorbed power calculations in fields.f)
 t11m(iR,itpr,1) = t11;    t12m(iR,itpr,1) = t12
 t21m(iR,itpr,1) = t21;    t22m(iR,itpr,1) = t22
                                                t33m(iR,itpr,1) = t33

 do i = 1, AP_Nb-1
  t11m(iR,itpr,i+1) = eS(AP_lkp(i));   t12m(iR,itpr,i+1) = - eD(AP_lkp(i))
  t21m(iR,itpr,i+1) = eD(AP_lkp(i));   t22m(iR,itpr,i+1) = eS(AP_lkp(i))
                                       t33m(iR,itpr,i+1) = eP(AP_lkp(i))
 enddo

! Optimization. Making a temporary static copy to avoid large dynamic array access

 TL(:,:)  = mTL (:,:,gR,itp)
 TLs(:,:) = mTLs(:,:,gR,itp)
 TLt(:,:) = mTLt(:,:,gR,itp)
 TLp(:,:) = mTLp(:,:,gR,itp)
 TU(:,:)  = mTU (:,:,gR,itp)
 TUs(:,:) = mTUs(:,:,gR,itp)
 TUt(:,:) = mTUt(:,:,gR,itp)
 TUp(:,:) = mTUp(:,:,gR,itp)
 jac = bjac(gR,itp); jacs = bjacs(gR,itp); jact = bjact(gR,itp); jacp = bjacp(gR,itp); 


! The following lines are generated in Mathematica
!   itpr Nvf i  k


! ----------------------------------------------------------------------------------------


! Ar     = ksi1 * sf1
! Atheta = ksi2 * sf2
! Aphi   = ksi3 * sf3
! phi4   = ksi4 * sf4

 sf1  = sf_v(1,eR,gR)
 sf2  = sf_v(2,eR,gR)
 sf3  = sf_v(3,eR,gR)
 sf4  = sf_v(4,eR,gR)

 dsf1 = dsf_v(1,eR,gR)
 dsf2 = dsf_v(2,eR,gR)
 dsf3 = dsf_v(3,eR,gR)
 dsf4 = dsf_v(4,eR,gR)

 ! Equlibrium asymptotics 
 seTU1  = seTU(1,eR,gR)
 seTU2  = seTU(2,eR,gR)
 seTL1  = seTL(1,eR,gR)
 seTL2  = seTL(2,eR,gR)

 dseTU1 = dseTU(1,eR,gR)
 dseTU2 = dseTU(2,eR,gR)
 dseTL1 = dseTL(1,eR,gR)
 dseTL2 = dseTL(2,eR,gR)


! Full eq. EQ + A asymptotic.

CF(itp,1) = -((seTL1*seTL2*sf2*TL(1,1)*TLp(2,2))/jac**2) + (seTL1*sf2*TL(1,1)*TLt(3,2))/jac**2
CF(itp,2) = -((seTL1*sf3*TL(1,1)*TLp(2,3))/jac**2) + (seTL1*sf3*TL(1,1)*TLt(3,3))/jac**2
CF(itp,3) = (seTL1*sf2*TL(1,1)*TL(3,2))/jac**2
CF(itp,4) = -((seTL1*seTL2*sf2*TL(1,1)*TL(2,2))/jac**2)
CF(itp,5) = (seTL1*sf3*TL(1,1)*TL(3,3))/jac**2
CF(itp,6) = -((seTL1*sf3*TL(1,1)*TL(2,3))/jac**2)
CF(itp,7) = (seTL1*seTL2*sf1*TL(2,2)*TLp(1,1))/jac**2 - (seTL1*sf1*TL(3,2)*TLt(1,1))/jac**2
CF(itp,8) = (dseTL2*sf2*TL(2,2)*TL(3,2))/jac**2 + (seTL2*sf2*TL(2,2)*TLp(1,2))/jac**2 &
- (seTL2*sf2*TL(1,2)*TLp(2,2))/jac**2 + (seTL2*sf2*TL(3,2)*TLs(2,2))/jac**2 &
- (seTL2*sf2*TL(2,2)*TLs(3,2))/jac**2 - (sf2*TL(3,2)*TLt(1,2))/jac**2 &
+ (sf2*TL(1,2)*TLt(3,2))/jac**2
CF(itp,9) = (dsf3*TL(2,3)*TL(3,2))/jac**2 - (dsf3*seTL2*TL(2,2)*TL(3,3))/jac**2 &
+ (seTL2*sf3*TL(2,2)*TLp(1,3))/jac**2 - (sf3*TL(1,2)*TLp(2,3))/jac**2 &
+ (sf3*TL(3,2)*TLs(2,3))/jac**2 - (seTL2*sf3*TL(2,2)*TLs(3,3))/jac**2 &
- (sf3*TL(3,2)*TLt(1,3))/jac**2 + (sf3*TL(1,2)*TLt(3,3))/jac**2
CF(itp,10) = -((seTL1*sf1*TL(1,1)*TL(3,2))/jac**2)
CF(itp,11) = (seTL1*seTL2*sf1*TL(1,1)*TL(2,2))/jac**2
CF(itp,12) = (sf3*TL(2,3)*TL(3,2))/jac**2 - (seTL2*sf3*TL(2,2)*TL(3,3))/jac**2
CF(itp,13) = -((sf3*TL(1,3)*TL(3,2))/jac**2) + (sf3*TL(1,2)*TL(3,3))/jac**2
CF(itp,14) = (seTL2*sf3*TL(1,3)*TL(2,2))/jac**2 - (sf3*TL(1,2)*TL(2,3))/jac**2
CF(itp,15) = (seTL1*sf1*TL(2,3)*TLp(1,1))/jac**2 - (seTL1*sf1*TL(3,3)*TLt(1,1))/jac**2
CF(itp,16) = -((dsf2*TL(2,3)*TL(3,2))/jac**2) + (dsf2*seTL2*TL(2,2)*TL(3,3))/jac**2 &
+ (dseTL2*sf2*TL(2,2)*TL(3,3))/jac**2 + (sf2*TL(2,3)*TLp(1,2))/jac**2 &
- (seTL2*sf2*TL(1,3)*TLp(2,2))/jac**2 + (seTL2*sf2*TL(3,3)*TLs(2,2))/jac**2 &
- (sf2*TL(2,3)*TLs(3,2))/jac**2 - (sf2*TL(3,3)*TLt(1,2))/jac**2 + (sf2*TL(1,3)*TLt(3,2))/jac**2
CF(itp,17) = (sf3*TL(2,3)*TLp(1,3))/jac**2 - (sf3*TL(1,3)*TLp(2,3))/jac**2 &
+ (sf3*TL(3,3)*TLs(2,3))/jac**2 - (sf3*TL(2,3)*TLs(3,3))/jac**2 &
- (sf3*TL(3,3)*TLt(1,3))/jac**2 + (sf3*TL(1,3)*TLt(3,3))/jac**2
CF(itp,18) = -((seTL1*sf1*TL(1,1)*TL(3,3))/jac**2)
CF(itp,19) = (seTL1*sf1*TL(1,1)*TL(2,3))/jac**2
CF(itp,20) = -((sf2*TL(2,3)*TL(3,2))/jac**2) + (seTL2*sf2*TL(2,2)*TL(3,3))/jac**2
CF(itp,21) = (sf2*TL(1,3)*TL(3,2))/jac**2 - (sf2*TL(1,2)*TL(3,3))/jac**2
CF(itp,22) = -((seTL2*sf2*TL(1,3)*TL(2,2))/jac**2) + (sf2*TL(1,2)*TL(2,3))/jac**2
CF(itp,23) = (dsf1*seTU1*TU(1,1))/jac + (dseTU1*sf1*TU(1,1))/jac &
+ (jacs*seTU1*sf1*TU(1,1))/jac**2 + (jact*sf1*TU(2,1))/jac**2 &
+ (jacp*sf1*TU(3,1))/jac**2 + (sf1*TUp(3,1))/jac + (seTU1*sf1*TUs(1,1))/jac &
+ (sf1*TUt(2,1))/jac
CF(itp,24) = (jact*seTU2*sf2*TU(2,2))/jac**2 + (jacp*sf2*TU(3,2))/jac**2 &
+ (sf2*TUp(3,2))/jac + (seTU2*sf2*TUt(2,2))/jac
CF(itp,25) = (jact*sf3*TU(2,3))/jac**2 + (jacp*sf3*TU(3,3))/jac**2 + (sf3*TUp(3,3))/jac &
+ (sf3*TUt(2,3))/jac
CF(itp,26) = (seTU1*sf1*TU(1,1))/jac
CF(itp,27) = (sf1*TU(2,1))/jac
CF(itp,28) = (sf1*TU(3,1))/jac
CF(itp,29) = (seTU2*sf2*TU(2,2))/jac
CF(itp,30) = (sf2*TU(3,2))/jac
CF(itp,31) = (sf3*TU(2,3))/jac
CF(itp,32) = (sf3*TU(3,3))/jac
CF(itp,33) = sf1
CF(itp,34) = 0
CF(itp,35) = 0
CF(itp,36) = sf2
CF(itp,37) = sf3
CF(itp,38) = sf1
CF(itp,39) = 0
CF(itp,40) = 0
CF(itp,41) = 0
CF(itp,42) = 0
CF(itp,43) = sf2
CF(itp,44) = 0
CF(itp,45) = 0
CF(itp,46) = 0
CF(itp,47) = 0
CF(itp,48) = sf3
CF(itp,49) = 0
CF(itp,50) = 0
CF(itp,51) = dsf4*seTU1*TU(1,1)
CF(itp,52) = seTU1*sf4*TU(1,1)
CF(itp,53) = sf4*TU(2,1)
CF(itp,54) = sf4*TU(3,1)
CF(itp,55) = 0
CF(itp,56) = 0
CF(itp,57) = seTU2*sf4*TU(2,2)
CF(itp,58) = sf4*TU(3,2)
CF(itp,59) = sf4*TU(2,3)
CF(itp,60) = sf4*TU(3,3)
CF(itp,61) = 0
CF(itp,62) = 0
CF(itp,63) = dsf4*seTU1*TU(1,1)
CF(itp,64) = seTU1*sf4*TU(1,1)
CF(itp,65) = sf4*TU(2,1)
CF(itp,66) = sf4*TU(3,1)
CF(itp,67) = 0
CF(itp,68) = 0
CF(itp,69) = seTU2*sf4*TU(2,2)
CF(itp,70) = sf4*TU(3,2)
CF(itp,71) = 0
CF(itp,72) = sf4*TU(2,3)
CF(itp,73) = sf4*TU(3,3)
CA(itp,1) = CAF(1)*(-(seTL1*seTL2*sf2*TL(1,1)*TLp(2,2)) + seTL1*sf2*TL(1,1)*TLt(3,2))
CA(itp,2) = CAF(1)*(-(seTL1*sf3*TL(1,1)*TLp(2,3)) + seTL1*sf3*TL(1,1)*TLt(3,3))
CA(itp,3) = seTL1*sf2*CAF(1)*TL(1,1)*TL(3,2)
CA(itp,4) = -(seTL1*seTL2*sf2*CAF(1)*TL(1,1)*TL(2,2))
CA(itp,5) = seTL1*sf3*CAF(1)*TL(1,1)*TL(3,3)
CA(itp,6) = -(seTL1*sf3*CAF(1)*TL(1,1)*TL(2,3))
CA(itp,7) = CAF(1)*(seTL1*seTL2*sf1*TL(2,2)*TLp(1,1) - seTL1*sf1*TL(3,2)*TLt(1,1))
CA(itp,8) = CAF(1)*(dseTL2*sf2*TL(2,2)*TL(3,2) + seTL2*sf2*TL(2,2)*TLp(1,2) &
- seTL2*sf2*TL(1,2)*TLp(2,2) + seTL2*sf2*TL(3,2)*TLs(2,2) - seTL2*sf2*TL(2,2)*TLs(3,2) &
- sf2*TL(3,2)*TLt(1,2) + sf2*TL(1,2)*TLt(3,2))
CA(itp,9) = CAF(1)*(dsf3*TL(2,3)*TL(3,2) - dsf3*seTL2*TL(2,2)*TL(3,3) &
+ seTL2*sf3*TL(2,2)*TLp(1,3) - sf3*TL(1,2)*TLp(2,3) + sf3*TL(3,2)*TLs(2,3) &
- seTL2*sf3*TL(2,2)*TLs(3,3) - sf3*TL(3,2)*TLt(1,3) + sf3*TL(1,2)*TLt(3,3))
CA(itp,10) = -(seTL1*sf1*CAF(1)*TL(1,1)*TL(3,2))
CA(itp,11) = seTL1*seTL2*sf1*CAF(1)*TL(1,1)*TL(2,2)
CA(itp,12) = CAF(1)*(sf3*TL(2,3)*TL(3,2) - seTL2*sf3*TL(2,2)*TL(3,3))
CA(itp,13) = CAF(1)*(-(sf3*TL(1,3)*TL(3,2)) + sf3*TL(1,2)*TL(3,3))
CA(itp,14) = CAF(1)*(seTL2*sf3*TL(1,3)*TL(2,2) - sf3*TL(1,2)*TL(2,3))
CA(itp,15) = CAF(1)*(seTL1*sf1*TL(2,3)*TLp(1,1) - seTL1*sf1*TL(3,3)*TLt(1,1))
CA(itp,16) = CAF(1)*(-(dsf2*TL(2,3)*TL(3,2)) + dsf2*seTL2*TL(2,2)*TL(3,3) &
+ dseTL2*sf2*TL(2,2)*TL(3,3) + sf2*TL(2,3)*TLp(1,2) - seTL2*sf2*TL(1,3)*TLp(2,2) &
+ seTL2*sf2*TL(3,3)*TLs(2,2) - sf2*TL(2,3)*TLs(3,2) - sf2*TL(3,3)*TLt(1,2) &
+ sf2*TL(1,3)*TLt(3,2))
CA(itp,17) = CAF(1)*(sf3*TL(2,3)*TLp(1,3) - sf3*TL(1,3)*TLp(2,3) + sf3*TL(3,3)*TLs(2,3) &
- sf3*TL(2,3)*TLs(3,3) - sf3*TL(3,3)*TLt(1,3) + sf3*TL(1,3)*TLt(3,3))
CA(itp,18) = -(seTL1*sf1*CAF(1)*TL(1,1)*TL(3,3))
CA(itp,19) = seTL1*sf1*CAF(1)*TL(1,1)*TL(2,3)
CA(itp,20) = CAF(1)*(-(sf2*TL(2,3)*TL(3,2)) + seTL2*sf2*TL(2,2)*TL(3,3))
CA(itp,21) = CAF(1)*(sf2*TL(1,3)*TL(3,2) - sf2*TL(1,2)*TL(3,3))
CA(itp,22) = CAF(1)*(-(seTL2*sf2*TL(1,3)*TL(2,2)) + sf2*TL(1,2)*TL(2,3))
CA(itp,23) = CAF(2)*(dsf1*jac*seTU1*TU(1,1) + dseTU1*jac*sf1*TU(1,1) &
+ jacs*seTU1*sf1*TU(1,1) + jact*sf1*TU(2,1) + jacp*sf1*TU(3,1) + jac*sf1*TUp(3,1) &
+ jac*seTU1*sf1*TUs(1,1) + jac*sf1*TUt(2,1))
CA(itp,24) = CAF(2)*(jact*seTU2*sf2*TU(2,2) + jacp*sf2*TU(3,2) + jac*sf2*TUp(3,2) &
+ jac*seTU2*sf2*TUt(2,2))
CA(itp,25) = CAF(2)*(jact*sf3*TU(2,3) + jacp*sf3*TU(3,3) + jac*sf3*TUp(3,3) + jac*sf3*TUt(2,3))
CA(itp,26) = jac*seTU1*sf1*CAF(2)*TU(1,1)
CA(itp,27) = jac*sf1*CAF(2)*TU(2,1)
CA(itp,28) = jac*sf1*CAF(2)*TU(3,1)
CA(itp,29) = jac*seTU2*sf2*CAF(2)*TU(2,2)
CA(itp,30) = jac*sf2*CAF(2)*TU(3,2)
CA(itp,31) = jac*sf3*CAF(2)*TU(2,3)
CA(itp,32) = jac*sf3*CAF(2)*TU(3,3)
CA(itp,33) = sf1*t11*CAF(3)
CA(itp,34) = sf2*t12*CAF(3)
CA(itp,35) = sf1*t21*CAF(3)
CA(itp,36) = sf2*t22*CAF(3)
CA(itp,37) = sf3*t33*CAF(3)
CA(itp,38) = 0
CA(itp,39) = dsf4*seTU1*t11*CAF(4)*TU(1,1)
CA(itp,40) = seTU1*sf4*t11*CAF(4)*TU(1,1)
CA(itp,41) = CAF(4)*(sf4*t11*TU(2,1) + seTU2*sf4*t12*TU(2,2))
CA(itp,42) = CAF(4)*(sf4*t11*TU(3,1) + sf4*t12*TU(3,2))
CA(itp,43) = 0
CA(itp,44) = dsf4*seTU1*t21*CAF(4)*TU(1,1)
CA(itp,45) = seTU1*sf4*t21*CAF(4)*TU(1,1)
CA(itp,46) = CAF(4)*(sf4*t21*TU(2,1) + seTU2*sf4*t22*TU(2,2))
CA(itp,47) = CAF(4)*(sf4*t21*TU(3,1) + sf4*t22*TU(3,2))
CA(itp,48) = 0
CA(itp,49) = sf4*t33*CAF(4)*TU(2,3)
CA(itp,50) = sf4*t33*CAF(4)*TU(3,3)
CA(itp,51) = dsf4*seTU1*t11*CAF(5)*TU(1,1)
CA(itp,52) = seTU1*sf4*t11*CAF(5)*TU(1,1)
CA(itp,53) = CAF(5)*(sf4*t11*TU(2,1) + seTU2*sf4*t12*TU(2,2))
CA(itp,54) = CAF(5)*(sf4*t11*TU(3,1) + sf4*t12*TU(3,2))
CA(itp,55) = dsf4*seTU1*t21*CAF(5)*TU(1,1)
CA(itp,56) = seTU1*sf4*t21*CAF(5)*TU(1,1)
CA(itp,57) = CAF(5)*(sf4*t21*TU(2,1) + seTU2*sf4*t22*TU(2,2))
CA(itp,58) = CAF(5)*(sf4*t21*TU(3,1) + sf4*t22*TU(3,2))
CA(itp,59) = sf4*t33*CAF(5)*TU(2,3)
CA(itp,60) = sf4*t33*CAF(5)*TU(3,3)
CA(itp,61) = sf1*t11*CAF(6)
CA(itp,62) = sf2*t12*CAF(6)
CA(itp,63) = 0
CA(itp,64) = 0
CA(itp,65) = 0
CA(itp,66) = 0
CA(itp,67) = sf1*t21*CAF(6)
CA(itp,68) = sf2*t22*CAF(6)
CA(itp,69) = 0
CA(itp,70) = 0
CA(itp,71) = sf3*t33*CAF(6)
CA(itp,72) = 0
CA(itp,73) = 0

! End of Mathematica output

 enddo ! itp


 Ncf=1
 do while(Ncf<73)
     
   Nsb=CFCAINFO(4,Ncf)
   do ni=1,Nsb
     do nj=1,Nsb


       k=CFCAINFO(3,Ncf+ni-1)
       kp=CFCAINFO(3,Ncf+nj-1)
       ii=CFCAINFO(2,Ncf+ni-1)
       jj=CFCAINFO(2,Ncf+nj-1)
       ij=ij_i_j_lkp(ii,jj)

       do itp = 1, Njk

         Acoeff(itp,k,kp,ij) = Acoeff(itp,k,kp,ij) + CA(itp,Ncf+ni-1) * CF(itp,Ncf+nj-1) * bjac(gR,itp) * wgauss(NGaussR,gR)

       enddo
     enddo
   enddo

   Ncf=Ncf+Nsb

 enddo

end subroutine
!-------------------------------------------
subroutine TransformEq2Fourier(eR,gR,k,kp)
 integer ij, eR,gR, ieq, k,kp
 real(dp) norm

 complex(dp) X(Njk), Y(Njk)
 complex(dp) Y0(0:Njk) ! used to avoid 'if' in the loop -- coefficients
                               ! outside mEqMax,nEqMax are automatically zero


  do ij = 1, (DIM+1)*(DIM+1)

      X = Acoeff(:,k,kp,ij)/sqrt(real(Njk))

      call fftwnd_f77_one(trans,X,Y)

      Y0(0) = 0
      Y0(1:Njk) = Y
      ! ieq: Get the position of (m,n) in the X,Y table
      do ieq = 1, mnEqTot
        Acoeff_k_f(ieq, ij) = Y0(mnFFT_lkp(ieq))
      enddo !ieq
  enddo !ij

  do ieq = 1, mnEqTot
  do ij = 1, (DIM+1)*(DIM+1)

      Acoeff_k_f(ieq, ij) = Acoeff_k_f(ieq, ij) / sqrt(real(Njk)) ! (2*pi*2*pi)/Njk/nperiod

      ! Fill the table of maximums for eq. Fourier amplitudes
      norm =  Acoeff_k_f(ieq, ij) * conjg(Acoeff_k_f(ieq, ij)) 
      if (norm > maxFEqAmp(ieq)) maxFEqAmp(ieq) = norm

  enddo
  enddo

!  do ieq = 1, mnEqTot
!    print *, mEq(ieq), nEq(ieq),  maxFEqAmp(ieq)
!  enddo
!  stop

if (FourierCheck)  call CheckFourier(k,kp)


end subroutine
!-------------------------------------------
subroutine CheckFourier(k,kp)
 integer ij, gR, itp,itpr, ieq, k,kp
 complex(dp) f1, f2
 real(dp) dist, norm

  do ij = 1, (DIM+1)*(DIM+1)
    do itpr = 1, Njkr
      gR  = gR_itpr_lkp(itpr)
      itp = itp_itpr_lkp(itpr)

      f1 = 0
      do ieq = 1, mnEqTot
        f1 = f1 + Acoeff_k_f(ieq, ij) / EqExp(ieq,itp)
      enddo
      f2 = Acoeff_k(itp,gR,ij)

      norm = f2*conjg(f2)
      if (norm>maxEqValue) maxEqValue = norm

      dist = (f1-f2)*conjg(f1-f2)

      if (dist>FTransformError)  FTransformError = dist
      if (FTransformError>1.0e-10_dp) then
        print *, 'Discrepancy in eq. Fourier transform:'
        print *, k,kp,ij,gR,itp,f1,f2  
        stop
      endif
    enddo
  enddo

end subroutine
!-------------------------------------------

!-------------------------------------------
! Main subroutine for the matrix construction!
 subroutine Init_Matrices(e,i)
 integer k,kp,e,a,b,i,m,mp,n,np,mn,mnp,iKK,jKK, iKK_mn,jKK_mn,  percent,pper
 complex(dp) ke(Nbf,Nbf), tmp
 integer mnmnp,ab

 ! Objects used in CalcLocMatrix_ke

 ! a,b - index of basic functions on one radial element (FE). For "hat" functions - 0 or 1.
 ! i,j - derivative index (0 - no derivative, 1,2,3 -- d/ds,d/dtheta,d/dphi)
 ! NoVarTerms -  number of additive terms in the variational form (6)
 ! NoSPTerms - number of additive terms in the scalar products (3 - dimension of space)
 ! Njk - combined poloidal & toroidal index, nj*nk
 ! gR - radial gauss point index (over the current (one!) radial element)

  
 integer AModified, FModified, iAF

 complex(dp) I_co, I_gR
 integer ij, gR, ieq, ip, status(MPI_STATUS_SIZE)


    ! Interpolate the equilibrium (TU,TL,B,etc) on the gauss radial
    ! grid for the radial element "e" (from TERPSICHORE grid)
    if (CylMetric) then
      call Init_Cyl_Equilibruim(e)
    else
      call Init_Equilibruim(e)
    endif
    cR_hp(e,:) = cR(NGaussR/2+1,:)
    cZ_hp(e,:) = cZ(NGaussR/2+1,:)
    phiv_hp(e,:) = phiv(NGaussR/2+1,:)

    call Init_dPsiR(e)


  do gR = 1, NGaussR

    call Init_a_eq(e,gR,i) ! Calculate the coefficients of the variational form on one
                      ! radial interval (gauss points&poloidal&toroidal grid)


      Acoeff = Acoeff * (s_nod(e+1)-s_nod(e)) ! for radial integration


      do kp = 1, NDIM
      do k  = 1, NDIM

      if ( ((PrecalcMatrix/=1).or.(OmSteps/=0)) .and. (EmptyRun==0) ) &
        call TransformEq2Fourier(e,gR,k,kp)


      do ij = 1, (DIM+1)*(DIM+1)
      do mnmnp = 1, mnAFTot*mnAFTot

        Acoeff_k_f_Dmult_imn(mnmnp,ij) = Acoeff_k_f(mneq0_lkp(mnmnp), ij) &
                                    &  * DMult_imn(mnmnp,ij,(nindex-1)/np_npar+1)

      enddo


!if (ij==5) then
!print *, Acoeff_k_f_Dmult_imn(1,5), Acoeff_k_f_Dmult_imn(2,5), Acoeff_k_f_Dmult_imn(3,5)
!stop
!end if

      enddo


!$omp parallel do private(mnmnp,mn,mnp,iAF,jKK,iKK,jKK_mn,iKK_mn,tmp,a,b)
      do ab = 1, Nbf*Nbf

        a = a_ab_lkp(ab)
        b = b_ab_lkp(ab)

        ! Part of the index in the KK matrix (without mn,mnp)
        jKK_mn = mod(a+1,2)*mnAFTot*NDIM + (k-1)*mnAFTot 
        iKK_mn = mod(b+1,2)*mnAFTot*NDIM + (kp-1)*mnAFTot

!CDIR LISTVEC, NODEP

      do mnmnp = 1, mnAFTot*mnAFTot

        mn  = mn_mnmnp_lkp(mnmnp)
        mnp = mnp_mnmnp_lkp(mnmnp)


        iAF = iAF_lkp(mnmnp,k,kp)

        tmp = 0
        tmp = tmp + Acoeff_k_f_Dmult_imn(mnmnp,1)  * dPsiRAF(gR,iAF,1,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,2)  * dPsiRAF(gR,iAF,2,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,3)  * dPsiRAF(gR,iAF,3,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,4)  * dPsiRAF(gR,iAF,4,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,5)  * dPsiRAF(gR,iAF,5,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,6)  * dPsiRAF(gR,iAF,6,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,7)  * dPsiRAF(gR,iAF,7,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,8)  * dPsiRAF(gR,iAF,8,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,9)  * dPsiRAF(gR,iAF,9,b,a)  &
                  + Acoeff_k_f_Dmult_imn(mnmnp,10) * dPsiRAF(gR,iAF,10,b,a) &
                  + Acoeff_k_f_Dmult_imn(mnmnp,11) * dPsiRAF(gR,iAF,11,b,a) &
                  + Acoeff_k_f_Dmult_imn(mnmnp,12) * dPsiRAF(gR,iAF,12,b,a) &
                  + Acoeff_k_f_Dmult_imn(mnmnp,13) * dPsiRAF(gR,iAF,13,b,a) &
                  + Acoeff_k_f_Dmult_imn(mnmnp,14) * dPsiRAF(gR,iAF,14,b,a) &
                  + Acoeff_k_f_Dmult_imn(mnmnp,15) * dPsiRAF(gR,iAF,15,b,a) &
                  + Acoeff_k_f_Dmult_imn(mnmnp,16) * dPsiRAF(gR,iAF,16,b,a)

        jKK = jKK_mn + mn   ! column number
        iKK = iKK_mn + mnp  ! row number


        if (ab_lkp(b,a)==1) then
          KKE(iKK,jKK)=KKE(iKK,jKK) + tmp * 2*pi*2*pi/nperiod
        elseif (ab_lkp(b,a)==2) then
          KK_tmp(iKK,jKK)=KK_tmp(iKK,jKK) + tmp * 2*pi*2*pi/nperiod
        elseif (ab_lkp(b,a)==3) then
          KKA(iKK,jKK)=KKA(iKK,jKK) + tmp * 2*pi*2*pi/nperiod
        else
          KKD(iKK,jKK)=KKD(iKK,jKK) + tmp * 2*pi*2*pi/nperiod
        endif

      enddo

      enddo
!$omp end parallel do

      enddo
      enddo

  enddo !gR

!if (sgme==0.and.m_npar==1) print *, "step : ", e


 end subroutine
!-------------------------------------------
 subroutine Init_RHS
 integer b,iFF,i,  eR, kp,mnp,  itpr,itp,gR
 complex(dp) tmp
 real(dp) tmpr, tmpi, radial_dep
 integer FModified, status(MPI_STATUS_SIZE)
 real(dp), dimension(:), allocatable   :: jext_r, jext_i
 real(dp), dimension(:,:,:), allocatable :: Exp_mnp_r, Exp_mnp_i
 complex(dp), dimension(:), allocatable :: FF_local

 allocate(jext_r(Njkr), jext_i(Njkr))
 allocate(Exp_mnp_r(Njk,mnAFTot,(nindex-1)/np_npar+1))
 allocate(Exp_mnp_i(Njk,mnAFTot,(nindex-1)/np_npar+1))
 allocate(FF_local(MSize2))

 Exp_mnp_r=real(Exp_mnp)
 Exp_mnp_i=imag(Exp_mnp)

  call walltime
  if (me.eq.0) then
    print *, 'RHS:'
  endif

  FF_local = 0._dp

  do eR = emin(msep), emax(msep)
!  do eR = 1, NelR

    ! Interpolate the equilibrium (TU,TL,B,etc) on the gauss radial
    ! grid for the radial element "e" (from TERPSICHORE grid)
    if (CylMetric) then
      call Init_Cyl_Equilibruim(eR)
    else
      call Init_Equilibruim(eR)
    endif

    call Init_jext3D(eR)

    do kp = 1, NDIM

    jext_r(:)=real(jext_3D(:,kp))
    jext_i(:)=imag(jext_3D(:,kp))

    do b = 1, Nbf

!CDIR LISTVEC, NODEP

!$omp parallel do private(FModified,tmpr,tmpi,gR,itp,itpr,radial_dep,tmp,i,iFF)
      do mnp = 1, mnAFTot

        FModified = iF_lkp(mnp,kp)

        tmpr = 0._dp
        tmpi = 0._dp
        do gR = 1, NGaussR

          do itp = 1, Njk

            itpr=(gR-1)*Njk+itp

            radial_dep = Psi_FE(0,b,eR, eR,gR,FModified) * sf_v(kp,eR,gR)* &
                       wgauss(NGaussR,gR)

            tmpr = tmpr + radial_dep * bjac(gR,itp) * &
                          (jext_r(itpr)*Exp_mnp_r(itp,mnp,(nindex-1)/np_npar+1)- &
                          jext_i(itpr)*Exp_mnp_i(itp,mnp,(nindex-1)/np_npar+1))

            tmpi = tmpi + radial_dep * bjac(gR,itp) * &
                          (jext_r(itpr)*Exp_mnp_i(itp,mnp,(nindex-1)/np_npar+1)+ &
                          jext_i(itpr)*Exp_mnp_r(itp,mnp,(nindex-1)/np_npar+1))

          enddo
        enddo
!        tmp = - tmp * (s_nod(eR+1)-s_nod(eR)) * (2*pi*2*pi)/Njk/nperiod

        tmp = tmpr+iim*tmpi

        tmp = - tmp * (s_nod(eR+1)-s_nod(eR)) * (2*pi*2*pi)/Njk/nperiod * &
                mu0_SI

        i = LMat(b,eR);
        iFF = (i-1)*mnAFTot*NDIM + (kp-1)*mnAFTot + mnp  ! row number
        FF_local(iFF) = FF_local(iFF) + tmp
      enddo !mnp
!$omp end parallel do

    enddo !b      
    enddo !kp

  enddo !eR

  call MPI_ALLREDUCE(FF_local,FF,MSize2,MPI_DOUBLE_COMPLEX,MPI_SUM, &
    sep_comm(m_npar),ierr)

  call walltime
  if (me.eq.0) print *, 'Matrix construction and factorisation:'

 deallocate(FF_local)
 deallocate(Exp_mnp_r,Exp_mnp_i)
 deallocate(jext_r,jext_i)

 end subroutine
!-------------------------------------------
 subroutine RHS_BConditions
 integer e,b,m,n,mp,np,k,kp,i,j,jj, mn,mnp,ii

! Axis

  e = 1;  b = 1
  do mn = 1, mnAFTot   ! boundary conditions: Ak = 0 at r=br0, r=br1
    mnp = mn
    select case (mAF(mn))
      case (0)
        do k = 1, 2
          kp = k
          i = GlobalIndex_lkp(e,b,mnp,kp)      ! row number (mp)
          FF(i) = 0._dp
        enddo

      case default
        do k = 1, NDIM
          kp = k
          i = GlobalIndex_lkp(e,b,mnp,kp)      ! row number (mp)
          FF(i) = 0._dp
        enddo
    end select
  enddo

! Outer boundary
  e = NelR;  b = Nbf-1
  do mn = 1, mnAFTot   ! boundary conditions: Ak = 0 at r=br0, r=br1
    mnp = mn
    do k = 2, NDIM  ! NDIM=4
      kp = k
      i = GlobalIndex_lkp(e,b,mnp,kp)      ! row number (mp)
      FF(NelR*LBKK+i) = 0._dp
    enddo
  enddo

 end subroutine
!-------------------------------------------
 subroutine Axis_BConditions
 integer e,b,m,n,mp,np,k,kp,i,j,jj, mn,mnp,ii

  e = 1;  b = 1
  do mn = 1, mnAFTot   ! boundary conditions: Ak = 0 at r=br0, r=br1
    mnp = mn
    select case (mAF(mn))
      case (0)
        do k = 1, 2
          kp = k
          j = GlobalIndex_lkp(e,b,mn,k)        ! column number (m)
          i = GlobalIndex_lkp(e,b,mnp,kp)      ! row number (mp)
          do jj = 1, LBKK
            KKE(i,jj) = 0._dp
            KK_tmp(i,jj) = 0._dp
          enddo
          KKE(i,j) = 1._dp
        enddo

      case default
        do k = 1, NDIM
          kp = k
          j = GlobalIndex_lkp(e,b,mn,k)        ! column number (m)
          i = GlobalIndex_lkp(e,b,mnp,kp)      ! row number (mp)
          do jj = 1, LBKK
            KKE(i,jj) = 0._dp
            KK_tmp(i,jj) = 0._dp
          enddo
          KKE(i,j) = 1._dp
        enddo
    end select
  enddo

end subroutine
!-------------------------------------------
 subroutine Border_BConditions
 integer e,b,m,n,mp,np,k,kp,i,j,jj, mn,mnp,ii

  e = NelR;  b = Nbf-1

  if (role<(min(npsep,4)+1)/2) then

  do mn = 1, mnAFTot   ! boundary conditions: Ak = 0 at r=br0, r=br1
    mnp = mn
    do k = 2, NDIM  ! NDIM=4
      kp = k
      j = GlobalIndex_lkp(e,b,mn,k)        ! column number (m)
      i = GlobalIndex_lkp(e,b,mnp,kp)      ! row number (mp)
      do jj = 1, LBKK
        KKA(i,jj) = 0._dp
        KKD(i,jj) = 0._dp
      enddo
      KKD(i,j) = 1._dp
    enddo
  enddo

  else

  do mn = 1, mnAFTot   ! boundary conditions: Ak = 0 at r=br0, r=br1
    mnp = mn
    do k = 2, NDIM  ! NDIM=4
      kp = k
      j = GlobalIndex_lkp(e,b,mn,k)        ! column number (m)
      i = GlobalIndex_lkp(e,b,mnp,kp)      ! row number (mp)
      do jj = 1, LBKK
        KK_tmp(i,jj) = 0._dp
        KKE(i,jj) = 0._dp
      enddo
      KKE(i,j) = 1._dp
    enddo
  enddo

  end if

 end subroutine
!-------------------------------------------
 subroutine GaussElimination(e,i)
  integer e,i,ip,ifile
  integer  INFO, status(MPI_STATUS_SIZE)
  integer, dimension(:), allocatable :: IPIV
  integer, parameter :: fn = 2
  complex(dp) alpha, beta
  real(dp) detime, artime
  character(5) stmp

  allocate(IPIV(LBKK))

  alpha = (1.0, 0.0)
  beta = (0.0, 0.0)

    KKC(:,:)=KK_tmp(:,sgme*LBKK/same+1:(sgme+1)*LBKK/same)


    if (npsep>2.and.e/=1.and.e/=NelR) then
      if (role<(min(npsep,4)+1)/2) then
        call MPI_RECV(KK_tmp,LBKK*LBKK,MPI_DOUBLE_COMPLEX, &
          &  me_lkp(mod(e/2,same),mod(e,2),m_npar),0,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(FF((e-1)*LBKK+1:e*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
          &  me_lkp(sgme,mod(role+1,2),m_npar),0,MPI_COMM_WORLD,status,ierr)
      else
        call MPI_RECV(KK_tmp,LBKK*LBKK,MPI_DOUBLE_COMPLEX, &
          &  me_lkp(mod((NelR+1-e)/2,same),mod(NelR-e,2)+2,m_npar),0,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(FF(e*LBKK+1:(e+1)*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
          &  me_lkp(sgme,mod(role+1,2)+2,m_npar),0,MPI_COMM_WORLD,status,ierr)
      end if
      KK_tmp(:,:)=KK_tmp(:,:)+KKE(:,:)
    else
      KK_tmp=KKE
    end if


    call ZGETRF(LBKK,LBKK,KK_tmp,LBKK,IPIV,INFO)

    call ZGETRS('N',LBKK,LBKK/same,KK_tmp,LBKK,IPIV, &
      &  KKC(1,1),LBKK,INFO)


    if (role<(min(npsep,4)+1)/2) then

      call ZGETRS('N',LBKK,1,KK_tmp,LBKK,IPIV, &
        &  FF((e-1)*LBKK+1),LBKK,INFO)

      call ZGEMV('N',LBKK,LBKK,alpha,KKA(1,1),LBKK, &
        &  FF((e-1)*LBKK+1),1,beta,FF_tmp,1)

      FF(e*LBKK+1:(e+1)*LBKK)=FF(e*LBKK+1:(e+1)*LBKK)-FF_tmp

    else

      call ZGETRS('N',LBKK,1,KK_tmp,LBKK,IPIV, &
        &  FF(e*LBKK+1),LBKK,INFO)

      call ZGEMV('N',LBKK,LBKK,alpha,KKA(1,1),LBKK, &
        &  FF(e*LBKK+1),1,beta,FF_tmp,1)

      if(e==NelR/2+1) then
        FF((e-1)*LBKK+1:e*LBKK)=-FF_tmp
      else
        FF((e-1)*LBKK+1:e*LBKK)=FF((e-1)*LBKK+1:e*LBKK)-FF_tmp
      end if

    end if

    call ZGEMM('N','N',LBKK,LBKK/same,LBKK,alpha,KKA(1,1), &
      &  LBKK,KKC,LBKK,beta,KK_mult,LBKK)


  if (sgme==mod(min((NelR-e)/2,(e-1)/2),same)) then
    do ip = 0, same-1
      if (sgme==ip) then
        KKE(:,sgme*LBKK/same+1:(sgme+1)*LBKK/same)=KKC(:,:)
      else
        call MPI_RECV(KKE(1,ip*LBKK/same+1),LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
          &  me_lkp(ip,role,m_npar),0,MPI_COMM_WORLD,status,ierr)
      end if
    enddo
    if (i<=indexfile) then
      KKB(:,:,i)=KKE(:,:)
    else
      ifile=0
      do while(ifile*LBKK_red < LBKK)
        write(stmp,'(i5)') nprocs*(NelR/2+1)*ifile+me*(NelR/2+1)+i
        stmp = adjustl(stmp)
 open(fn, file=STORAGEPATH//'/BL'//trim(stmp)//'.dat',form='unformatted')
        write(fn) KKE(:,ifile*LBKK_red+1:min(LBKK,(ifile+1)*LBKK_red))
        close(fn)
        ifile=ifile+1
      end do
    endif
    esup = e
  else
    call MPI_SEND(KKC,LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
      &  me_lkp(mod(min((NelR-e)/2,(e-1)/2),same),role,m_npar),0,MPI_COMM_WORLD,ierr)
  end if


  KKD(:,sgme*LBKK/same+1:(sgme+1)*LBKK/same) = &
    KKD(:,sgme*LBKK/same+1:(sgme+1)*LBKK/same) - KK_mult(:,:)


  if (sgme==mod(min((NelR-e+2)/2,(e+1)/2),same)) then
    do ip = 0, same-1
      if (sgme/=ip) then
        call MPI_RECV(KKD(1,ip*LBKK/same+1),LBKK*LBKK/same, &
          &  MPI_DOUBLE_COMPLEX,me_lkp(ip,role,m_npar),0,MPI_COMM_WORLD,status,ierr)
      end if
    enddo
  else
    call MPI_SEND(KKD(1,sgme*LBKK/same+1),LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
      &  me_lkp(mod(min((NelR-e+2)/2,(e+1)/2),same),role,m_npar),0,MPI_COMM_WORLD,ierr)
  end if


  if (npsep>2.and.e/=NelR/2.and.e/=NelR/2+1) then
    if (sgme==mod(min((NelR-e+2)/2,(e+1)/2),same)) then
      do ip = 0, same-1
        call MPI_SEND(KKD(1,1),LBKK*LBKK,MPI_DOUBLE_COMPLEX, &
          &  me_lkp(ip,role+1-2*mod(role,2),m_npar),0,MPI_COMM_WORLD,ierr)
      enddo
    end if


    if (role<(min(npsep,4)+1)/2) then
      call MPI_SEND(FF(e*LBKK+1:(e+1)*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
        me_lkp(sgme,mod(role+1,2),m_npar),0,MPI_COMM_WORLD,ierr)
    else
      call MPI_SEND(FF((e-1)*LBKK+1:e*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
        me_lkp(sgme,mod(role+1,2)+2,m_npar),0,MPI_COMM_WORLD,ierr)
    end if
  end if



    if(e==NelR/2+1.and.npsep>1) then
      if (sgme==0) then
        do ip = 0, same-1
          if (npsep==2) then
            call MPI_RECV(KK_tmp(1,ip*LBKK/same+1),LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
            &  me_lkp(ip,0,m_npar), &
            &  0,MPI_COMM_WORLD,status,ierr)
          else
            call MPI_RECV(KK_tmp(1,ip*LBKK/same+1),LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
            &  me_lkp(ip,mod(NelR/2+1,2),m_npar), &
            &  0,MPI_COMM_WORLD,status,ierr)
          end if
          if (sgme/=ip) then
            call MPI_RECV(KKD(1,ip*LBKK/same+1),LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
            &  me_lkp(ip,role,m_npar),0,MPI_COMM_WORLD,status,ierr)
          end if
        enddo

        if (npsep==2) then
          call MPI_RECV(FF_tmp,LBKK,MPI_DOUBLE_COMPLEX, &
            me_lkp(0,0,m_npar),0,MPI_COMM_WORLD,status,ierr)
        else
          call MPI_RECV(FF_tmp,LBKK,MPI_DOUBLE_COMPLEX, &
            me_lkp(0,mod(NelR/2+1,2),m_npar),0,MPI_COMM_WORLD,status,ierr)
        end if

        KK_tmp(:,:)=KK_tmp(:,:)+KKD(:,:)
        FF((e-1)*LBKK+1:e*LBKK)=FF((e-1)*LBKK+1:e*LBKK)+FF_tmp(:)

        call ZGETRF(LBKK,LBKK,KK_tmp,LBKK,IPIV,INFO)
        call ZGETRS('N',LBKK,1,KK_tmp,LBKK,IPIV, &
          &  FF((e-1)*LBKK+1),LBKK,INFO)

      else
        call MPI_SEND(KKD(1,sgme*LBKK/same+1),LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
          &  me_lkp(0,role,m_npar),0,MPI_COMM_WORLD,ierr)
      end if
    end if





    if (e==NelR/2.and.npsep>1) then
      if (npsep==2) then
        call MPI_SEND(KKD(1,sgme*LBKK/same+1),LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
          me_lkp(0,1,m_npar),0,MPI_COMM_WORLD,ierr)
      else
        call MPI_SEND(KKD(1,sgme*LBKK/same+1),LBKK*LBKK/same,MPI_DOUBLE_COMPLEX, &
          me_lkp(0,3-mod(NelR-NelR/2-1,2),m_npar),0,MPI_COMM_WORLD,ierr)
      end if
      if (npsep==2) then
        call MPI_SEND(FF(e*LBKK+1:(e+1)*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
          me_lkp(0,1,m_npar),0,MPI_COMM_WORLD,ierr)
      elseif (sgme==0) then
        call MPI_SEND(FF(e*LBKK+1:(e+1)*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
          me_lkp(0,3-mod(NelR-NelR/2-1,2),m_npar),0,MPI_COMM_WORLD,ierr)
      end if
    end if

  deallocate(IPIV)

 end subroutine
!-------------------------------------------
 subroutine SolveEq(i)

 integer e,a,NN,k,r,j,mn,i,l,ip,ifile

 integer  INFO, status(MPI_STATUS_SIZE)
 integer, dimension(:), allocatable :: IPIV
 integer, parameter :: fn = 2
 character(1) TRANS
 complex(dp) alpha, beta
 character(5) stmp
 
  allocate(IPIV(LBKK))

  alpha = (1.0, 0.0)
  beta = (0.0, 0.0)

  if (npsep==1) then
    bmin = gmin
    bmax = gmax
  else
    if (role<(min(npsep,4))/2) then
      bmin = gmin+2*sgme
      bmax = NelR/2
    else
      bmin = esup
      bmax = NelR
    end if
  end if


  if (EmptyRun==0) then

  call walltime
  if (me.eq.0) print *, 'Backsolve:'

  if (npsep==1) then

    KK_tmp=KKD(:,:)
    call ZGETRF(LBKK,LBKK,KK_tmp,LBKK,IPIV,INFO)
    call ZGETRS('N',LBKK,1,KK_tmp,LBKK,IPIV, &
      &  FF(NelR*LBKK+1),LBKK,INFO)

! Backsolve

    do i = NelR, 1, -1

      if (i<=indexfile) then
        KKE(:,:)=KKB(:,:,i)
      else
        ifile = 0
        do while(ifile*LBKK_red < LBKK)
          write(stmp,'(i5)') nprocs*(NelR/2+1)*ifile+me*(NelR/2+1)+i
          stmp = adjustl(stmp)
 open(fn, file=STORAGEPATH//'BL'//trim(stmp)//'.dat',form='unformatted',action='read')
          read(fn) KKE(:,ifile*LBKK_red+1:min(LBKK,(ifile+1)*LBKK_red))
          close(fn)
          ifile=ifile+1
        end do
      end if

      call ZGEMV('N',LBKK,LBKK,alpha,KKE(1,1), &
        &  LBKK,FF(i*LBKK+1),1,beta,FF_tmp,1)

      FF((i-1)*LBKK+1:i*LBKK)=FF((i-1)*LBKK+1:i*LBKK)-FF_tmp

    enddo

  else

    i=i-1

    if (npsep==2) then
      call MPI_BCAST(FF(NelR/2*LBKK+1:(NelR/2+1)*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
        &  me_lkp(0,1,1),sep_comm(m_npar),ierr) ! m_npar set to 1 because it is not the full
                                                ! communicator
    else
      call MPI_BCAST(FF(NelR/2*LBKK+1:(NelR/2+1)*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
        &  me_lkp(0,3-mod(NelR-NelR/2-1,2),1),sep_comm(m_npar),ierr)
    end if

! Burn At Both Ends (BABE)

    if (role<(min(npsep,4)+1)/2) then

! First processor

    do e = esup, 1, -incr*same

      i=i-1

      if (npsep>2.and.e/=NelR/2) then
        call MPI_RECV(FF(e*LBKK+1:(e+1)*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
          me_lkp(mod(sgme+role,same),mod(role+1,2),m_npar),0,MPI_COMM_WORLD,status,ierr)
      end if

      if (i<indexfile) then
        KKE(:,:)=KKB(:,:,i+1)
      else
        ifile=0
        do while(ifile*LBKK_red < LBKK)
          write(stmp,'(i5)') nprocs*(NelR/2+1)*ifile+me*(NelR/2+1)+i+1
          stmp = adjustl(stmp)
 open(fn, file=STORAGEPATH//'BL'//trim(stmp)//'.dat',form='unformatted',action='read')
          read(fn) KKE(:,ifile*LBKK_red+1:min(LBKK,(ifile+1)*LBKK_red))
          close(fn)
          ifile=ifile+1
        end do
      end if

      call ZGEMV('N',LBKK,LBKK,alpha,KKE, &
      &  LBKK,FF(e*LBKK+1),1,beta,FF_tmp,1)

      FF((e-1)*LBKK+1:e*LBKK)=FF((e-1)*LBKK+1:e*LBKK)-FF_tmp

      if (npsep>2.and.e/=1) then
        call MPI_SEND(FF((e-1)*LBKK+1:e*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
          me_lkp(mod(sgme+role-1+same,same),mod(role+1,2),m_npar),0,MPI_COMM_WORLD,ierr)
      end if
 
    enddo

    else

    do e = esup, NelR, -incr*same

      i=i-1

      if (npsep>2.and.e/=NelR/2+1) then
        call MPI_RECV(FF((e-1)*LBKK+1:e*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
          me_lkp(mod(sgme+mod(role+1,2),same),mod(role+1,2)+2,m_npar),0,MPI_COMM_WORLD,status,ierr)
      end if

      if (i<indexfile) then
        KKE(:,:)=KKB(:,:,i+1)
      else
        ifile = 0
        do while(ifile*LBKK_red < LBKK)
          write(stmp,'(i5)') nprocs*(NelR/2+1)*ifile+me*(NelR/2+1)+i+1
          stmp = adjustl(stmp)
 open(fn, file=STORAGEPATH//'BL'//trim(stmp)//'.dat',form='unformatted',action='read')
          read(fn) KKE(:,ifile*LBKK_red+1:min(LBKK,(ifile+1)*LBKK_red))
          close(fn)
          ifile=ifile+1
        end do
      end if

      call ZGEMV('N',LBKK,LBKK,alpha,KKE, &
      &  LBKK,FF((e-1)*LBKK+1),1,beta,FF_tmp,1)


      FF(e*LBKK+1:(e+1)*LBKK)=FF(e*LBKK+1:(e+1)*LBKK)-FF_tmp

      if (npsep>2.and.e/=NelR) then
        call MPI_SEND(FF(e*LBKK+1:(e+1)*LBKK),LBKK,MPI_DOUBLE_COMPLEX, &
          me_lkp(mod(sgme-mod(role,2)+same,same),mod(role+1,2)+2,m_npar),0,MPI_COMM_WORLD,ierr)
      end if

    enddo

    endif

    ! Transmission of the results to the first processor


    if (msep.eq.0) then
      do e = 1, NelR
        if (e>NelR/2) then
          if (npsep>2) then
            ip=me_lkp(mod((NelR-e)/2,same),mod(NelR+1-e,2)+2,m_npar)
          else
            ip=me_lkp(0,1,m_npar)
          end if
          call MPI_RECV(FF((e-1)*LBKK+1),LBKK,MPI_DOUBLE_COMPLEX,ip,0, &
            MPI_COMM_WORLD,status,ierr)
        elseif (mod(e,incr*same)/=1.and.npsep>2) then
          ip=me_lkp(mod((e-1)/2,same),mod(e-1,2),m_npar)
          call MPI_RECV(FF((e-1)*LBKK+1),LBKK,MPI_DOUBLE_COMPLEX,ip,0, &
            MPI_COMM_WORLD,status,ierr)
        end if
      enddo
      call MPI_RECV(FF(NelR*LBKK+1),LBKK,MPI_DOUBLE_COMPLEX, &
        me_lkp(0,min(npsep-1,3),m_npar),0,MPI_COMM_WORLD,status,ierr)
    else
      do e = 1, NelR
        if(mod(e,abs(same*incr))==mod(bmin,abs(same*incr)) &
          .and.e>=bmin.and.e<=bmax) then
          call MPI_SEND(FF((e-1)*LBKK+1),LBKK,MPI_DOUBLE_COMPLEX, &
            me_lkp(0,0,m_npar),0,MPI_COMM_WORLD,ierr)
        end if
      enddo
      if (me.eq.me_lkp(0,min(npsep-1,3),m_npar)) then
        call MPI_SEND(FF(NelR*LBKK+1),LBKK,MPI_DOUBLE_COMPLEX, &
          me_lkp(0,0,m_npar),0,MPI_COMM_WORLD,ierr)
      end if
    end if

  endif

endif

    do mn = 1, mnAFTot
      do e = 1, NelR
      do a = 1, Nbf
      do k = 1, NDIM
        r = LMat(a,e)
        NN = (r-1)*mnAFTot*NDIM + (k-1)*mnAFTot + mn
        XX(r,mn,k) = FF(NN)
      enddo
      enddo
      enddo
    enddo

  do k = 1, NDIM
    if (m_npar==1) then
      do ip = 2, np_npar
        call MPI_RECV(XXc(1,(nindex+ip-2)*mnAFTot+1,k),Nunknowns*mnAFTot,MPI_DOUBLE_COMPLEX, &
          me_lkp(sgme,role,ip),0,MPI_COMM_WORLD,status,ierr)
      enddo
    else
      call MPI_SEND(XX(1,1,k),Nunknowns*mnAFTot,MPI_DOUBLE_COMPLEX, &
        me_lkp(sgme,role,1),0,MPI_COMM_WORLD,ierr)
    end if
  enddo

  XXc(:,(nindex-1)*mnAFTot+1:nindex*mnAFTot,:) = XX(:,1:mnAFTot,:)

  deallocate(IPIV)


 end subroutine
!-------------------------------------------


!-------------------------------------------

 end module matrix
