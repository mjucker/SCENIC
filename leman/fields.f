!
!
!  E/M waves propagation in 3D plasmas  -- Appendix.
!
!  Calculation of the E/M fields from the potentials (in Fourier space)
!
!  Version 3.0: 3D, general geometry, fixed boundary
!
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!==========================================================================
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


 module Fields

 use math
 use matrix
 
 implicit none


 ! grid for E&B
 integer, parameter ::  NtotR = NelR * NGaussR

 ! Fields
 complex(dp) PdivCum,                            &
             Pant(Nmesh), PantSum(Nmesh),        &
             Pdiv(Nmesh), PdivSum(Nmesh),        &
             Ppoynt(NelR), PpoyntSum(NelR),      &
             PAdivA(NelR), PAdivASum(NelR)

 complex(dp), dimension(:,:), allocatable :: Ppla, PplaSum
 complex(dp), dimension(:), allocatable   :: PplaCum, EeE

 complex(dp), dimension(:), allocatable   :: scalar_temp
 complex(dp), dimension(:,:), allocatable :: vector_temp, quadri_temp
 real(dp), dimension(:), allocatable      :: real_temp
 complex(dp), dimension(:,:,:), allocatable :: k_par_smooth

 real(dp) divAint(NelR,NGaussR), divAintR,  TotalVolume, &
          avrgA, avrgEdgeA,  avrgFluxDiff, maxFluxDiff

 real(dp) TotalVolumeSum, divAintRSum, avrgASum, &
          avrgEdgeASum

 complex(dp), dimension(:,:,:,:,:,:), allocatable :: vdAph_dui_g   ! (3,..,3)=(j,i)  ui:  1 - d/ds,  2 - d/dtheta,  3 - d/dphi

 real(dp) maxDiv, maxA


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 contains
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


!-------------------------------------------
subroutine EvaluateFields
  integer k

  call walltime
  if (me.eq.0) print *, 'Treatment of the solution:'

  vA_g = 0; vAp_g = 0
  do k = 1, NDIM
    call A_IntervalValues(k)
  enddo
  
  call CalculatePowerFlux
  if (me.eq.0) call CoulombGauge
!  call FourierTransformEB

end subroutine
!-------------------------------------------
subroutine A_IntervalValues(kk)
  integer kk,mn
  integer  i, eR, gR, AModified, a, j
  real(dp)  s, sf, dsf

  j = 0
  do nindex = m_npar, nwidth, np_npar
  j = j+1

  do mn = 1, mnAFTot

  AModified = iA_lkp(mn,kk)

  ! Gauss points
  do eR = bmin, bmax, abs(same*incr)
  do gR = 1, NGaussR

    sf  = sf_v(kk,eR,gR)
    dsf = dsf_v(kk,eR,gR)

    vA_g (kk,mn,eR,gR,j)  = (XXc(LMat(1,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(0,1,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(2,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(0,2,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(3,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(0,3,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(4,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(0,4,eR,s_gauss(eR,gR),AModified)) * sf
                                             
    vAp_g(kk,mn,eR,gR,j)  = (XXc(LMat(1,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(1,1,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(2,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(1,2,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(3,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(1,3,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(4,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(1,4,eR,s_gauss(eR,gR),AModified)) * sf + &
                                             
                            (XXc(LMat(1,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(0,1,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(2,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(0,2,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(3,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(0,3,eR,s_gauss(eR,gR),AModified) + &
                             XXc(LMat(4,eR),mn+(nindex-1)*mnAFTot,kk)*Psi_al(0,4,eR,s_gauss(eR,gR),AModified)) * dsf

  enddo
  enddo

  ! Half-mesh points
  gR = NGaussR/2+1  ! NGaussR should be odd!!!
  do eR = bmin, bmax, abs(same*incr)
    vA_hp (kk, mn, eR) = vA_g (kk, mn, eR, gR, j) ! This quantity is not used afterwards
    vAp_hp(kk, mn, eR) = vAp_g(kk, mn, eR, gR, j)
  enddo

  enddo
  enddo


  do mn = 1, nwidth*mnAFTot

  ! FE grid
  do eR = bmin, bmax, abs(same*incr)

    sf  = sf_v_A(kk,eR)
    dsf = dsf_v_A(kk,eR)

    vA_fe (kk, mn, eR)  = (XXc(LMat(1,eR),mn,kk)*Psi_al(0,1,eR,s_nod(eR),AModified) + &
                           XXc(LMat(2,eR),mn,kk)*Psi_al(0,2,eR,s_nod(eR),AModified)) * sf
                                           
    vAp_fe(kk, mn, eR)  = (XXc(LMat(1,eR),mn,kk)*Psi_al(1,1,eR,s_nod(eR),AModified) + &
                           XXc(LMat(2,eR),mn,kk)*Psi_al(1,2,eR,s_nod(eR),AModified)) * sf + &
                                           
                          (XXc(LMat(1,eR),mn,kk)*Psi_al(0,1,eR,s_nod(eR),AModified) + &
                           XXc(LMat(2,eR),mn,kk)*Psi_al(0,2,eR,s_nod(eR),AModified)) * dsf

  enddo

  eR = NelR+1

  sf  = sf_v_A(kk,NelR+1)
  dsf = dsf_v_A(kk,NelR+1)
  vA_fe (kk, mn, NelR+1)  = (XXc(LMat(3,NelR),mn,kk)*Psi_al(0,3,NelR,s_nod(eR),AModified) + &
                             XXc(LMat(4,NelR),mn,kk)*Psi_al(0,4,NelR,s_nod(eR),AModified)) * sf
                                               
  vAp_fe(kk, mn, NelR+1)  = (XXc(LMat(3,NelR),mn,kk)*Psi_al(1,3,NelR,s_nod(eR),AModified) + &
                             XXc(LMat(4,NelR),mn,kk)*Psi_al(1,4,NelR,s_nod(eR),AModified)) * sf + &
                                               
                            (XXc(LMat(3,NelR),mn,kk)*Psi_al(0,3,NelR,s_nod(eR),AModified) + &
                             XXc(LMat(4,NelR),mn,kk)*Psi_al(0,4,NelR,s_nod(eR),AModified)) * dsf

  enddo

end subroutine
!-------------------------------------------
subroutine CalculatePowerFlux
 integer, parameter :: fn = 2

 integer eR,gR, mn, ip, status(MPI_STATUS_SIZE)
 integer itp, itpr, iR, aR, i, j, k, v_div, ieq
 complex(dp) B2, Sn, Ej, &
             PantCum, PpoyntCum, PAdivACum


 complex(dp) A3D(NDIM,Njkr),            & ! <- n,b,||,phi components, one radial element
             A3D_edgemodes(NDIM,Njkr),  & ! <- n,b,||,phi components, one radial element, edge Fourier modes
             dAdu(Njkr,DIM,DIM), Rk(Njkr,DIM), &
             dAph_dui_3D(3,Njkr,3)  ! (3,..,3)=(j,i)  ui:  1 - d/ds,  2 - d/dtheta,  3 - d/dphi
 complex(dp) Spoynt_ph(3,Njk)  ! (3,..,3)=(j,i)  ui:  1 - d/ds,  2 - d/dtheta,  3 - d/dphi

 complex(dp) AA, Fj

 real(dp) EpsTens(3,3,3), s,  TUtmp(3,3),TLtmp(3,3), TUstmp(3,3),TUttmp(3,3),TUptmp(3,3)

 real(dp) dTLdui(3,3,Njkr,3) ! (3,3,...,3)=(iii,jjj,ui)

 integer iii,jjj,kkk
 real(dp) maxFluxTmp, maxw_norm, sum_node, sum_total

 complex(dp) X(Njk), Y(Njk)
 complex(dp) Y0(0:Njk) ! used to avoid 'if' in the loop -- coefficients
                               ! outside mEqMax,nEqMax are automatically zero

 real(dp) v_th, v_the, v_thi1, v_thi2, v_thih, OmCIe, OmCIi1, OmCIi2, OmCIih
 real(dp) ApproxVls(Njk,2), v_par
 real(dp) FFT_Vls(0:Njk,2)
 integer mnp, mnmnp, INFO
 real(dp), dimension(:,:), allocatable  :: Vls_Mat
 real(dp), dimension(:), allocatable    :: Vls_RHS
 integer, dimension(:), allocatable     :: IPIV

 real(dp), dimension(:,:), allocatable  :: k_par_new, k_perp
 real(dp), dimension(:), allocatable    :: k_par_tmp
 complex(dp), dimension(:), allocatable :: k_cmp_tmp
 real(dp), dimension(:), allocatable    :: k_perp_rad, k_perp_bin


 complex(dp) k_tmp, inv_temp(Njkr)

  allocate(Vls_Mat(mnDTTot,mnDTTot))
  allocate(Vls_RHS(mnDTTot))
  allocate(IPIV(mnDTTot))
  allocate(k_par_new(NbDistr-1,Njkr))
  allocate(k_par_tmp(Njkr))
  allocate(k_cmp_tmp(Njkr))
  allocate(k_perp_rad(Njkr))
  allocate(k_perp_bin(Njkr))
  allocate(k_perp(NelR,Njk))


! Antisymmetric tensor
  EpsTens = 0
  EpsTens(1,2,3) = 1
  EpsTens(2,3,1) = 1
  EpsTens(3,1,2) = 1
  EpsTens(2,1,3) = -1
  EpsTens(1,3,2) = -1
  EpsTens(3,2,1) = -1


! All three power fluxes in one radial loop

  Ppla   = 0;
  Pant   = 0;
  Ppoynt = 0;  PpoyntCum = 0
  PAdivA = 0;  PAdivACum = 0
  Pdiv   = 0;
  Ppla3D = 0

  divA3D = 0; divAint = 0; divAintR = 0; TotalVolume = 0;  avrgA = 0;  avrgEdgeA = 0
  vdAph_dui_g = 0

  do aR = 1, NelR

  eR = (aR+1)/2 + (NelR+1)/2*mod(aR+1,2) ! Different order for the radial surfaces treatment in order
                                         ! to avoid waiting time between processors


  if(mod(eR,abs(same*incr))==mod(bmin,abs(same*incr)) &
    .and.eR>=bmin.and.eR<=bmax) then

    if(role<(min(npsep,4)+1)/2) then
      iR=(eR-1)/((npsep+1)/2)+1
    else
      iR=(NelR-eR)/((npsep+1)/2)+1
    end if

    if (CylMetric) then
      call Init_Cyl_Equilibruim(eR)
    else
      call Init_Equilibruim(eR)
    endif
    cR_hp(eR,:) = cR(NGaussR/2+1,:)
    cZ_hp(eR,:) = cZ(NGaussR/2+1,:)
    phiv_hp(eR,:) = phiv(NGaussR/2+1,:)

    call Init_jext3D(eR)

  E3D = 0; B3D = 0; GradPhi3D = 0; A3D = 0; k_par_new = 0

  j = 0
  do nindex = m_npar, nwidth, np_npar
  j = j+1
!$omp parallel do private(itp,gR,mn)
  do itpr = 1, Njkr
    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    do mn = 1, mnAFTot

      GradPhi3D(1,itpr) = GradPhi3D(1,itpr) + &
                                   (mTU(1,1,gR,itp)*seTU(1,eR,gR)*vAp_g(4,mn,eR,gR,j) +                   &
                                    mTU(2,1,gR,itp)*iim*mAF(mn)*vA_g(4,mn,eR,gR,j) +    &
                                    mTU(3,1,gR,itp)*iim*nAF(mn,j)*vA_g(4,mn,eR,gR,j)) * Exp_mn(itp,mn,j) 
      GradPhi3D(2,itpr) = GradPhi3D(2,itpr) + (mTU(1,2,gR,itp)*vAp_g(4,mn,eR,gR,j) +                 &
                                    mTU(2,2,gR,itp)*seTU(2,eR,gR)*iim*mAF(mn)*vA_g(4,mn,eR,gR,j) +      &
                                    mTU(3,2,gR,itp)*iim*nAF(mn,j)*vA_g(4,mn,eR,gR,j)) * Exp_mn(itp,mn,j)
      GradPhi3D(3,itpr) = GradPhi3D(3,itpr) + (mTU(1,3,gR,itp)*vAp_g(4,mn,eR,gR,j) +                 &
                                    mTU(2,3,gR,itp)*iim*mAF(mn)*vA_g(4,mn,eR,gR,j) +    &
                                    mTU(3,3,gR,itp)*iim*nAF(mn,j)*vA_g(4,mn,eR,gR,j)) * Exp_mn(itp,mn,j)


      A3D(1,itpr) = A3D(1,itpr) + vA_g(1,mn,eR,gR,j) * Exp_mn(itp,mn,j)
      A3D(2,itpr) = A3D(2,itpr) + vA_g(2,mn,eR,gR,j) * Exp_mn(itp,mn,j)
      A3D(3,itpr) = A3D(3,itpr) + vA_g(3,mn,eR,gR,j) * Exp_mn(itp,mn,j)
      A3D(4,itpr) = A3D(4,itpr) + vA_g(4,mn,eR,gR,j) * Exp_mn(itp,mn,j)


      E3D(1,itpr) = E3D(1,itpr) + iim*Om_SI*vA_g(1,mn,eR,gR,j) * Exp_mn(itp,mn,j)
      E3D(2,itpr) = E3D(2,itpr) + iim*Om_SI*vA_g(2,mn,eR,gR,j) * Exp_mn(itp,mn,j)
      E3D(3,itpr) = E3D(3,itpr) + iim*Om_SI*vA_g(3,mn,eR,gR,j) * Exp_mn(itp,mn,j)

    enddo

  enddo
!$omp end parallel do
  enddo

  E3D(:,:) = - GradPhi3D(:,:)*c_SI + E3D(:,:)


! Calculation of k_par and k_perp based on scalar potential

   sum_total = sum(abs(vA_g(4,:,eR,NGaussR/2+1,j)))
   k_par_tmp = 0.0

   j = 1
!$omp parallel do private(itp,gR,mn)
   do itpr = 1, Njkr

    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    do mn = 1, mnDTTot
      k_par_tmp(itpr) = k_par_tmp(itpr)+abs(mfpp(gR)*mAF(mn_mnDT_lkp(mn))+mftp(gR)*nAF(mn_mnDT_lkp(mn),j)) &
         /abs(bjac(gR,itp)*B02*bmod(gR,itp))*abs(vA_g(4,mn,eR,NGaussR/2+1,j))/sum_total
!      k_cmp_tmp(itpr) = k_cmp_tmp(itpr)+abs(mfpp(gR)*mAF(mn_mnDT_lkp(mn))+mftp(gR)*nAF(mn_mnDT_lkp(mn),j)) &
!       /abs(bjac(gR,itp)*B02*bmod(gR,itp))*vA_g(4,mn,eR,NGaussR/2+1,j) * Exp_mn(itp,mn,j)
    enddo
!    k_par_tmp(itpr) = abs(k_cmp_tmp(itpr)/A3D(4,itpr))
   enddo
!$omp end parallel do
 
!  k_par_it(NbDistr,eR,:) = (1.0-RelaxParam)*k_par_it(NbDistr,eR,:) + RelaxParam*abs(k_par_tmp(:))
  k_par_it(NbDistr,eR,:) = abs(k_par_tmp(:))


   k_perp_bin = 0.0

   j = 1
!$omp parallel do private(itp,gR,mn)
   do itpr = 1, Njkr

    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    do mn = 1, mnDTTot
      k_perp_bin(itpr) = k_perp_bin(itpr)+abs(mci(gR)*mAF(mn_mnDT_lkp(mn))+mcj(gR)*nAF(mn_mnDT_lkp(mn),j)) &
         /abs(bjac(gR,itp)*B02*bmod(gR,itp))*abs(vA_g(4,mn,eR,NGaussR/2+1,j))/sum_total
    enddo
   enddo
!$omp end parallel do


   k_perp_rad = 0.0

   j = 1
!$omp parallel do private(itp,gR,mn)
   do itpr = 1, Njkr

    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    do mn = 1, mnDTTot
      k_perp_rad(itpr) = k_perp_rad(itpr)+abs(vAp_g(4,mn,eR,NGaussR/2+1,j))/sum_total
    enddo
   enddo
!$omp end parallel do


  j = 1
!$omp parallel do private(itp,gR,mn)
   do itpr = 1, Njkr

    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    if (gR==NGaussR/2+1) then
      k_perp(eR,itp) = sqrt(abs(mgssu(gR,itp)) * k_perp_rad(itpr)**2 + 1./abs(mgssu(gR,itp)) * k_perp_bin(itpr)**2)
    end if
   enddo
!$omp end parallel do



  do gR = 1, NGaussR

! Thermal velocities
  v_the =sqrt(2*Te0*fte(eR,gR)*1.6e-19/me_SI)
  v_thi1=sqrt(2*Ti0*fti(eR,gR)*1.6e-19/mi1_SI)
  v_thi2=sqrt(2*Ti0*ftm(eR,gR)*1.6e-19/mi2_SI)
  v_thih=sqrt(2*Th0*fth(eR,gR)*1.6e-19/mih_SI)

! All species

  k = 0

  do i = 1, 7

  select case (i)
   case(1)
    v_th = v_the
   case(2)
    v_th = v_thi1
   case(3)
    v_th = v_thi2
   case(4)
    v_th = v_thih
   case(5)
    v_th = v_thi1
   case(6)
    v_th = v_thi2
   case(7)
    v_th = v_thih
  end select

  if (k_par_lkp(i)==1) then
  k = k+1

  do itp = 1, Njk

   OmCIe  = qe_Gauss*B0*1.e4_dp*bmod(gR,itp)/(c_SI*100*me_SI*1.e3_dp)
   OmCIi1 = iZ1*qe_Gauss*B0*1.e4_dp*bmod(gR,itp)/(c_SI*100*mi1_SI*1.e3_dp)
   OmCIi2 = iZ2*qe_Gauss*B0*1.e4_dp*bmod(gR,itp)/(c_SI*100*mi2_SI*1.e3_dp)
   OmCIih = iZh*qe_Gauss*B0*1.e4_dp*bmod(gR,itp)/(c_SI*100*mih_SI*1.e3_dp)

   if (i==5) then
    ApproxVls(itp,1) = Om_SI-OmCIi1
   elseif (i==6) then
    ApproxVls(itp,1) = Om_SI-OmCIi2
   elseif (i==7) then
    ApproxVls(itp,1) = Om_SI-OmCIih
   else
    ApproxVls(itp,1) = Om_SI
   end if

   ApproxVls(itp,2) = 1/(bjac(gR,itp)*B02*bmod(gR,itp))

  enddo

  do j = 1, 2

   X = ApproxVls(:,j)
   call fftwnd_f77_one(trans,X,Y)
   Y0(0) = 0
   Y0(1:Njk) = Y

   do ieq = 1, mnEqDTTot
   FFT_Vls(ieq,j) = Y0(mnFFT_DT_lkp(ieq))
   enddo

  enddo

 maxw_norm = 0.0


 do v_div = 1, Nb_intv

  v_par = -2.0*v_range*v_th*(v_div-Nb_intv/2.0-0.5)/real(Nb_intv)

  if (i<=4) then
  maxw_norm = maxw_norm + (v_par/v_th)**2 * exp(-(v_par/v_th)**2)
  else
  maxw_norm = maxw_norm + exp(-(v_par/v_th)**2)
  end if

 enddo

 do v_div = 1, Nb_intv

  sum_node  = 0.0
  k_par_tmp = 0.0

  j = 0
  do nindex = m_npar, nwidth, np_npar
  j = j+1

   v_par = -2.0*v_range*v_th*(v_div-Nb_intv/2.0-0.5)/real(Nb_intv)


   do mnmnp = 1, mnDTTot*mnDTTot

    mnp = mn_mnmnpDT_lkp(mnmnp)
    mn = mnp_mnmnpDT_lkp(mnmnp)


    Vls_Mat(mnp,mn) = FFT_Vls(mneq0DT_lkp(mnmnp),1) + &
      v_par*(mftp(gR)*nDT(mnp,j)+mfpp(gR)*mDT(mnp)) * &
      FFT_Vls(mneq0DT_lkp(mnmnp),2)

   enddo

   if (i<=4) then
    X = E3D(3,(gR-1)*Njk+1:gR*Njk)
   else
    X = E3D(1,(gR-1)*Njk+1:gR*Njk)+iim*E3D(2,(gR-1)*Njk+1:gR*Njk)
   end if

   do mn = 1, mnAFTot
    if (mAF(mn)==m_ant.and.nAF(mn,j)==n_ant) then
      do itp = 1, Njk
       X(itp) = X(itp)*Exp_mnp(itp,mn,j)
      enddo
    end if
   enddo

   call fftwnd_f77_one(trans,X,Y)
   Y0(0) = 0
   Y0(1:Njk) = Y

   do mn = 1, mnDTTot
    Vls_RHS(mn) = abs(Y0(mnFFT_RHS_lkp(mn)))
   enddo

   call DGETRF(mnDTTot,mnDTTot,Vls_Mat,mnDTTot,IPIV,INFO)
   call DGETRS('N',mnDTTot,1,Vls_Mat,mnDTTot,IPIV,Vls_RHS,mnDTTot,INFO)

!  do mn = 1, mnDTTot
!   if (abs(Vls_RHS(mn))<0.04*abs(maxval(Vls_RHS))) then
!    Vls_RHS(mn)=0.0
!   end if
!  enddo

!  if (eR>NelR-15) then
!  do mn = 1, mnDTTot
!   if (abs(Vls_RHS(mn))<0.2*abs(maxval(Vls_RHS))) then
!    Vls_RHS(mn)=0.0
!   end if
!  enddo
!  end if

   sum_node = sum_node + sum(abs(Vls_RHS))

!$omp parallel do private(itpr,mn)
   do itp = 1, Njk

   itpr=itpr_itp_gR_lkp(itp,gR)

    do mn = 1, mnDTTot
      k_par_tmp(itpr) = k_par_tmp(itpr) + abs(mfpp(gR)*mDT(mn)+mftp(gR)*nDT(mn,j)) &
         /abs(bjac(gR,itp)*B02*bmod(gR,itp))*abs(Vls_RHS(mn))
    enddo
   enddo
!$omp end parallel do
  enddo ! j


  call MPI_REDUCE(sum_node,sum_total,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,nind_comm(msep),ierr)
  call MPI_REDUCE(k_par_tmp,real_temp,Njkr,MPI_DOUBLE_PRECISION,MPI_SUM,0,nind_comm(msep),ierr)
  if (m_npar==1) k_par_tmp(:) = real_temp

   if (i<=4) then
    k_par_new(k,:) = k_par_new(k,:) + &
         k_par_tmp(:) * (v_par/v_th)**2 * exp(-(v_par/v_th)**2)/(maxw_norm*sum_total)
   else
    k_par_new(k,:) = k_par_new(k,:) + &
         k_par_tmp(:) * exp(-(v_par/v_th)**2)/(maxw_norm*sum_total)
   end if

 enddo ! v_div

 end if

 enddo ! i

enddo

! Relaxation method and treatment of the border

 do i = 2, NbDistr-1
  if (eR<=NelR-BorderSurf) then
   k_par_it(i,eR,:) = (1.0-RelaxParam)*k_par_it(i,eR,:) + RelaxParam*abs(k_par_new(i-1,:))
  end if
 enddo

  call MPI_REDUCE(E3D,vector_temp,3*Njkr,MPI_DOUBLE_COMPLEX,MPI_SUM,0,nind_comm(msep),ierr)
  if (m_npar==1) E3D = vector_temp


  ! A3D edge modes (on the edge of the Fourier box)

  j = 0
  do nindex = m_npar, nwidth, np_npar
  j = j+1
  A3D_edgemodes = 0
  do mn = 1, EdgeModeNo
    do itpr = 1, Njkr
      itp = itp_itpr_lkp(itpr)
      gR  = gR_itpr_lkp(itpr)

      A3D_edgemodes(1,itpr) = A3D_edgemodes(1,itpr) + &
   vA_g(1,EdgeModes_lkp(mn),eR,gR,j) * Exp_mn(itp,EdgeModes_lkp(mn),j)
      A3D_edgemodes(2,itpr) = A3D_edgemodes(2,itpr) + &
   vA_g(2,EdgeModes_lkp(mn),eR,gR,j) * Exp_mn(itp,EdgeModes_lkp(mn),j)
      A3D_edgemodes(3,itpr) = A3D_edgemodes(3,itpr) + &
   vA_g(3,EdgeModes_lkp(mn),eR,gR,j) * Exp_mn(itp,EdgeModes_lkp(mn),j)
      A3D_edgemodes(4,itpr) = A3D_edgemodes(4,itpr) + &
   vA_g(4,EdgeModes_lkp(mn),eR,gR,j) * Exp_mn(itp,EdgeModes_lkp(mn),j)
    enddo
  enddo
  enddo



  ! B

  ! a. dAph_j/du_i in Fourier space

  j = 0
  do nindex = m_npar, nwidth, np_npar
  j = j+1
  do mn = 1, mnAFTot
    do gR = 1, NGaussR

      vdAph_dui_g(1,mn,eR,gR,1,j) = vAp_g(1,mn,eR,gR,j)
      vdAph_dui_g(1,mn,eR,gR,2,j) = iim*mAF(mn)*vA_g(1,mn,eR,gR,j)
      vdAph_dui_g(1,mn,eR,gR,3,j) = iim*nAF(mn,j)*vA_g(1,mn,eR,gR,j)
      
      vdAph_dui_g(2,mn,eR,gR,1,j) = vAp_g(2,mn,eR,gR,j)
      vdAph_dui_g(2,mn,eR,gR,2,j) = iim*mAF(mn)*vA_g(2,mn,eR,gR,j)
      vdAph_dui_g(2,mn,eR,gR,3,j) = iim*nAF(mn,j)*vA_g(2,mn,eR,gR,j)
      
      vdAph_dui_g(3,mn,eR,gR,1,j) = vAp_g(3,mn,eR,gR,j)
      vdAph_dui_g(3,mn,eR,gR,2,j) = iim*mAF(mn)*vA_g(3,mn,eR,gR,j)
      vdAph_dui_g(3,mn,eR,gR,3,j) = iim*nAF(mn,j)*vA_g(3,mn,eR,gR,j)

    enddo
  enddo
  enddo


  ! dAph_j/du_i in real space (3D)
  dAph_dui_3D = 0
  j = 0
  do nindex = m_npar, nwidth, np_npar
  j = j+1
  do mn = 1, mnAFTot
    do itpr = 1, Njkr
      itp = itp_itpr_lkp(itpr)
      gR  = gR_itpr_lkp(itpr)
      do iii = 1, 3
      do jjj = 1, 3
        dAph_dui_3D(jjj,itpr,iii) = dAph_dui_3D(jjj,itpr,iii) + vdAph_dui_g(jjj,mn,eR,gR,iii,j)*Exp_mn(itp,mn,j)
      enddo
      enddo

    enddo
  enddo
  enddo



  ! b. dA_j/du_i in real space
  dAdu = 0
  do itpr = 1, Njkr
    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    do iii = 1, 3
    do jjj = 1, 3
      dTLdui(iii,jjj,itpr,1) = mTLs(iii,jjj,gR,itp)
      dTLdui(iii,jjj,itpr,2) = mTLt(iii,jjj,gR,itp)
      dTLdui(iii,jjj,itpr,3) = mTLp(iii,jjj,gR,itp)
    enddo
    enddo
    dTLdui(1,1,itpr,1) = mTLs(1,1,gR,itp)*seTL(1,eR,gR) + mTL(1,1,gR,itp)*dseTL(1,eR,gR)
    dTLdui(2,2,itpr,1) = mTLs(2,2,gR,itp)*seTL(2,eR,gR) + mTL(2,2,gR,itp)*dseTL(2,eR,gR)
    dTLdui(1,1,itpr,2) = mTLt(1,1,gR,itp)*seTL(1,eR,gR)
    dTLdui(2,2,itpr,2) = mTLt(2,2,gR,itp)*seTL(2,eR,gR)
    dTLdui(1,1,itpr,3) = mTLp(1,1,gR,itp)*seTL(1,eR,gR)
    dTLdui(2,2,itpr,3) = mTLp(2,2,gR,itp)*seTL(2,eR,gR)

    TLtmp = mTL(:,:,gR,itp)
    TLtmp(1,1) = TLtmp(1,1)*seTL(1,eR,gR)
    TLtmp(2,2) = TLtmp(2,2)*seTL(2,eR,gR)

    do iii = 1, 3
    do jjj = 1, 3
    do kkk = 1, 3
      dAdu(itpr,jjj,iii) = dAdu(itpr,jjj,iii) + &
                    TLtmp(jjj,kkk)*dAph_dui_3D(kkk,itpr,iii) + A3D(kkk,itpr)*dTLdui(jjj,kkk,itpr,iii)
    enddo
    enddo
    enddo
  enddo


  ! c.  B = rot A in real space, 3D
  do itpr = 1, Njkr
    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    TLtmp = mTL(:,:,gR,itp)
    TLtmp(1,1) = TLtmp(1,1)*seTL(1,eR,gR)
    TLtmp(2,2) = TLtmp(2,2)*seTL(2,eR,gR)

    do iii = 1, 3
    do jjj = 1, 3
    do kkk = 1, 3
      B3D(1,itpr) = B3D(1,itpr) + EpsTens(iii,jjj,kkk)*dAdu(itpr,jjj,iii)*TLtmp(kkk,1)/bjac(gR,itp)
      B3D(2,itpr) = B3D(2,itpr) + EpsTens(iii,jjj,kkk)*dAdu(itpr,jjj,iii)*TLtmp(kkk,2)/bjac(gR,itp)
      B3D(3,itpr) = B3D(3,itpr) + EpsTens(iii,jjj,kkk)*dAdu(itpr,jjj,iii)*TLtmp(kkk,3)/bjac(gR,itp)
    enddo
    enddo
    enddo
  enddo


  call MPI_REDUCE(B3D,vector_temp,3*Njkr,MPI_DOUBLE_COMPLEX,MPI_SUM,0,nind_comm(msep),ierr)
  if (m_npar==1) B3D = vector_temp


  ! divA
  do kkk = 1, DIM
  do itpr = 1, Njkr
    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    TUtmp = mTU(:,:,gR,itp)
    TUtmp(1,1) = TUtmp(1,1)*seTU(1,eR,gR)
    TUtmp(2,2) = TUtmp(2,2)*seTU(2,eR,gR)

    TUstmp = mTUs(:,:,gR,itp)
    TUstmp(1,1) = mTUs(1,1,gR,itp)*seTU(1,eR,gR) + mTU(1,1,gR,itp)*dseTU(1,eR,gR)
    TUstmp(2,2) = mTUs(2,2,gR,itp)*seTU(2,eR,gR) + mTU(2,2,gR,itp)*dseTU(2,eR,gR)

    TUttmp = mTUt(:,:,gR,itp)
    TUttmp(1,1) = mTUt(1,1,gR,itp)*seTU(1,eR,gR)
    TUttmp(2,2) = mTUt(2,2,gR,itp)*seTU(2,eR,gR)

    TUptmp = mTUp(:,:,gR,itp)
    TUptmp(1,1) = mTUp(1,1,gR,itp)*seTU(1,eR,gR)
    TUptmp(2,2) = mTUp(2,2,gR,itp)*seTU(2,eR,gR)

    divA3D(itpr,eR) = divA3D(itpr,eR) + 1._dp/bjac(gR,itp)*( &

       ! sqrt(g) * TU_ik * dA_k/du_i
       bjac(gR,itp)*TUtmp(1,kkk)*dAph_dui_3D(kkk,itpr,1) + &
       bjac(gR,itp)*TUtmp(2,kkk)*dAph_dui_3D(kkk,itpr,2) + &
       bjac(gR,itp)*TUtmp(3,kkk)*dAph_dui_3D(kkk,itpr,3) + &
  
       ! sqrt(g) * dTU_ik/du_i * A_k
       bjac(gR,itp)*TUstmp(1,kkk)*A3D(kkk,itpr) + &
       bjac(gR,itp)*TUttmp(2,kkk)*A3D(kkk,itpr) + &
       bjac(gR,itp)*TUptmp(3,kkk)*A3D(kkk,itpr) + &
  
       ! dsqrt(g)/du_i * TU_ik * A_k
       bjacs(gR,itp)*TUtmp(1,kkk)*A3D(kkk,itpr) + &
       bjact(gR,itp)*TUtmp(2,kkk)*A3D(kkk,itpr) + &
       bjacp(gR,itp)*TUtmp(3,kkk)*A3D(kkk,itpr))
  enddo
  enddo

  call MPI_REDUCE(A3D,quadri_temp,4*Njkr,MPI_DOUBLE_COMPLEX,MPI_SUM,0,nind_comm(msep),ierr)
  if (m_npar==1) A3D = quadri_temp
  call MPI_REDUCE(divA3D(:,eR),scalar_temp,Njkr,MPI_DOUBLE_COMPLEX,MPI_SUM,0,nind_comm(msep),ierr)
  if (m_npar==1) divA3D(:,eR) = scalar_temp

!  ! divA
!  do kkk = 1, DIM
!  do itpr = 1, Njkr
!    itp = itp_itpr_lkp(itpr)
!    gR  = gR_itpr_lkp(itpr)
!    s   = s_gauss(eR,gR)
!
!    divA3D(eR,itpr) = divA3D(eR,itpr) + &
!   (-A3D(1,itpr)/2/Sqrt(s)) + dAph_dui_3D(3,itpr,3)/2 - &
!   dAph_dui_3D(2,itpr,2)/2/Sqrt(s) - sqrt(s)*dAph_dui_3D(1,itpr,1)
!
!  enddo
!  enddo


  do itpr = 1, Njkr
    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)
    divAint(eR,gR) = divAint(eR,gR) + abs(divA3D(itpr,eR))
!    if (s_gauss(eR,gR) < 0.18) &
    divAintR = divAintR + abs(divA3D(itpr,eR))*bjac(gR,itp)*wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))
    TotalVolume = TotalVolume + bjac(gR,itp)*wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))
    avrgA = avrgA + sqrt(abs( A3D(1,itpr)*conjg(A3D(1,itpr)) + A3D(2,itpr)*conjg(A3D(2,itpr)) + &
                              A3D(3,itpr)*conjg(A3D(3,itpr)) )) &
                         *bjac(gR,itp)*wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))
    avrgEdgeA = avrgEdgeA + sqrt(abs( A3D_edgemodes(1,itpr)*conjg(A3D_edgemodes(1,itpr)) + &
                                      A3D_edgemodes(2,itpr)*conjg(A3D_edgemodes(2,itpr)) + &
                                      A3D_edgemodes(3,itpr)*conjg(A3D_edgemodes(3,itpr)) )) &
                         *bjac(gR,itp)*wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))
  enddo


  ! E x B   --  Poynting flux in real space, center of the radial cell
  do itp = 1, Njk
    gR = NGaussR/2+1
    itpr = (gR-1)*Njk + itp

    ! Cross Product c = a x b
    ! c1 = a2b3 - a3b2
    ! c2 = a3b1 - a1b3
    ! c3 = a1b2 - a2b1
    Spoynt_ph(1,itp) = conjg(E3D(2,itpr))*B3D(3,itpr) - conjg(E3D(3,itpr))*B3D(2,itpr)
    Spoynt_ph(2,itp) = conjg(E3D(3,itpr))*B3D(1,itpr) - conjg(E3D(1,itpr))*B3D(3,itpr)
    Spoynt_ph(3,itp) = conjg(E3D(1,itpr))*B3D(2,itpr) - conjg(E3D(2,itpr))*B3D(1,itpr)


!    Spoynt_ph(1,itp) = conjg(A3D(2,itpr))*B3D(3,eR,itpr) - conjg(A3D(3,itpr))*B3D(2,eR,itpr)
!    Spoynt_ph(2,itp) = conjg(A3D(3,itpr))*B3D(1,eR,itpr) - conjg(A3D(1,itpr))*B3D(3,eR,itpr)
!    Spoynt_ph(3,itp) = conjg(A3D(1,itpr))*B3D(2,eR,itpr) - conjg(A3D(2,itpr))*B3D(1,eR,itpr)
!
!
!    Spoynt_ph(1,itp) = Spoynt_ph(1,itp) + conjg(A3D(1,itpr))*divA3D(eR,itpr)
!    Spoynt_ph(2,itp) = Spoynt_ph(2,itp) + conjg(A3D(2,itpr))*divA3D(eR,itpr)
!    Spoynt_ph(3,itp) = Spoynt_ph(3,itp) + conjg(A3D(3,itpr))*divA3D(eR,itpr)
  enddo

  ! Poynting flux
  PpoyntCum = 0
  do itp = 1, Njk
    gR = NGaussR/2+1
    itpr = (gR-1)*Njk + itp

    PpoyntCum = PpoyntCum + iim*0.5_dp/mu0_SI*Spoynt_ph(1,itp)*mTU(1,1,gR,itp)*seTU(1,eR,gR)*bjac(gR,itp)

!    do iii = 1, 3
!      ! E X B
!      PpoyntCum = PpoyntCum + iim*Spoynt_ph(iii,itp)*mTU(1,iii,gR,itp)*seTU(1,eR,gR)*bjac(gR,itp)
!!      PpoyntCum = PpoyntCum + Spoynt_ph(iii,itp)*mTU(1,iii,gR,itp)*bjac(gR,itp)
!
!!      ! F . V A
!!      PpoyntCum = PpoyntCum + conjg(A3D(iii,itpr))*divA3D(eR,itpr)*mTU(1,iii,gR,itp)*bjac(gR,itp)
!    enddo
  enddo


  
  ! F*divA - error (flux)
  PAdivACum = 0
  do itp = 1, Njk
    gR = NGaussR/2+1
    itpr = (gR-1)*Njk + itp

    PAdivACum = PAdivACum + conjg(A3D(1,itpr))*divA3D(itpr,eR)*mTU(1,1,gR,itp)*seTU(1,eR,gR)*bjac(gR,itp)
    PAdivACum = PAdivACum + conjg(A3D(2,itpr))*divA3D(itpr,eR)*mTU(1,2,gR,itp)*bjac(gR,itp)
    PAdivACum = PAdivACum + conjg(A3D(3,itpr))*divA3D(itpr,eR)*mTU(1,3,gR,itp)*bjac(gR,itp)
  enddo
  PAdivACum = PAdivACum * Om_SI/2/mu0_SI


  PplaCum = 0; PantCum = 0; PdivCum = 0;

  do itpr = 1, Njkr
    itp = itp_itpr_lkp(itpr)
    gR  = gR_itpr_lkp(itpr)

    Ej = conjg(E3D(1,itpr))*jext_3D(itpr,1) + conjg(E3D(2,itpr))*jext_3D(itpr,2) + &
         conjg(E3D(3,itpr))*jext_3D(itpr,3)

    Fj = conjg(A3D(1,itpr))*jext_3D(itpr,1) + conjg(A3D(2,itpr))*jext_3D(itpr,2) + &
         conjg(A3D(3,itpr))*jext_3D(itpr,3)

    B2 = B3D(1,itpr)*conjg(B3D(1,itpr)) + B3D(2,itpr)*conjg(B3D(2,itpr)) + B3D(3,itpr)*conjg(B3D(3,itpr))

    AA = A3D(1,itpr)*conjg(A3D(1,itpr)) + A3D(2,itpr)*conjg(A3D(2,itpr)) + &
         A3D(3,itpr)*conjg(A3D(3,itpr)) + A3D(4,itpr)*conjg(A3D(4,itpr))

    EeE(:) = conjg(E3D(1,itpr)) * (t11m(iR,itpr,:)*E3D(1,itpr) + t12m(iR,itpr,:)*E3D(2,itpr)) + &
             conjg(E3D(2,itpr)) * (t21m(iR,itpr,:)*E3D(1,itpr) + t22m(iR,itpr,:)*E3D(2,itpr)) + &
             conjg(E3D(3,itpr)) *  t33m(iR,itpr,:)*E3D(3,itpr)

!    PplaCum = PplaCum +  (B2 + AA + conjg(divA3D(eR,itpr))*divA3D(eR,itpr))*bjac(gR,itp)* &
!                         wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))
!
!
!    PantCum = PantCum +  4.0_dp*pi/c_SI * Fj*bjac(gR,itp)* &
!                         wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))


!!!Testing
!    PplaCum = PplaCum +  Om_SI/c_SI*(-B2 + EeE - conjg(divA3D(eR,itpr))*divA3D(eR,itpr) )*bjac(gR,itp)* &
!                         wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))

    PdivCum = PdivCum +  Om_SI/2/mu0_SI*(conjg(divA3D(itpr,eR))*divA3D(itpr,eR) )*bjac(gR,itp)* &
                         wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))


    PplaCum(:) = PplaCum(:) +  Om_SI/2/mu0_SI*(- 1.0_dp/c_SI/c_SI*EeE(:))*bjac(gR,itp)* &
                         wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))

    PplaCum(1) = PplaCum(1) +  Om_SI/2/mu0_SI*B2*bjac(gR,itp)* &
                         wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))


    PantCum = PantCum +  iim/2.0_dp * Ej*bjac(gR,itp)* &
                         wgauss(NGaussR,gR)*(s_nod(eR+1)-s_nod(eR))

    if (gR == NGaussR/2+1) then
      Ppla3D(:,eR,itp) = Om_SI/2/mu0_SI*(- 1.0_dp/c_SI/c_SI*EeE(:))*bjac(gR,itp)* &
                         (s_nod(eR+1)-s_nod(eR))*(2*pi*2*pi/Njk)/nperiod
      Ppla3D(1,eR,itp) = Ppla3D(1,eR,itp) + Om_SI/2/mu0_SI*B2*bjac(gR,itp)* &
                         (s_nod(eR+1)-s_nod(eR))*(2*pi*2*pi/Njk)/nperiod
      E3D_hp(:,eR,itp) = E3D(:,itpr)
      B3D_hp(:,eR,itp) = B3D(:,itpr)
      A3D_hp(:,eR,itp) = A3D(:,itpr)
      i = 1
      if (OutputSelect(5)==1) then
!        t_hp(i,eR,itp) = sum(t11m(iR,itpr,:))
        t_hp(i,eR,itp) = t11m(iR,itpr,3) ! first species
        i = i+1
      end if
      if (OutputSelect(6)==1) then
!        t_hp(i,eR,itp) = sum(t21m(iR,itpr,:))
        t_hp(i,eR,itp) = t11m(iR,itpr,4) !  second species
        i = i+1
      end if
      if (OutputSelect(7)==1) then
!        t_hp(i,eR,itp) = sum(t33m(iR,itpr,:))
        t_hp(i,eR,itp) = t11m(iR,itpr,5) ! fast ions
      end if
    end if

  enddo

    if(me.ne.0.and.m_npar==1) then
      call MPI_SEND(k_par_it(:,eR,:),NbDistr*Njkr, &
        &  MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(k_perp(eR,:),Njk, &
        &  MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(t_hp(:,eR,:),sum(OutputSelect(5:7))*Njk, &
        &  MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(Ppla3D(:,eR,:),AP_nb*Njk, &
        &  MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(E3D_hp(:,eR,:),DIM*Njk, &
        &  MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(B3D_hp(:,eR,:),DIM*Njk, &
        &  MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(A3D_hp(:,eR,:),DIM*Njk, &
        &  MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(divA3D(:,eR),Njkr, &
        &  MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(cR_hp(eR,:),Njk, &
        &  MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(cZ_hp(eR,:),Njk, &
        &  MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(phiv_hp(eR,:),Njk, &
        &  MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
      call MPI_SEND(vA_fe(:,:,eR),NDIM*nwidth*mnAFTot, &
        &  MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
    end if

  Ppla(:,eR+1) = PplaCum(:) * (2*pi*2*pi/Njk)/nperiod !   Ppla(1) = 0
  Pant(eR+1)   = PantCum * (2*pi*2*pi/Njk)/nperiod !   Pant(1) = 0
  Ppoynt(eR)   = PpoyntCum * (2*pi*2*pi/Njk)/nperiod
  PAdivA(eR)   = PAdivACum * (2*pi*2*pi/Njk)/nperiod
  Pdiv(eR+1)   = PdivCum * (2*pi*2*pi/Njk)/nperiod !   Ppla(1) = 0

  else

    if(me.eq.0) then

      if(eR<=NelR/2) then
        ip=mod(2*(eR-1)-mod(eR-1,2),npsep)
      else
        if(npsep>2) then
          ip=mod(2*(NelR-eR)+3*mod(NelR-eR+1,2),npsep)
        else
          ip=1
        end if
      end if

      call MPI_RECV(k_par_it(:,eR,:),NbDistr*Njkr, &
        &  MPI_DOUBLE_PRECISION,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(k_perp(eR,:),Njk, &
        &  MPI_DOUBLE_PRECISION,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(t_hp(:,eR,:),sum(OutputSelect(5:7))*Njk, &
        &  MPI_DOUBLE_PRECISION,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(Ppla3D(:,eR,:),AP_nb*Njk, &
        &  MPI_DOUBLE_COMPLEX,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(E3D_hp(:,eR,:),DIM*Njk, &
        &  MPI_DOUBLE_COMPLEX,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(B3D_hp(:,eR,:),DIM*Njk, &
        &  MPI_DOUBLE_COMPLEX,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(A3D_hp(:,eR,:),DIM*Njk, &
        &  MPI_DOUBLE_COMPLEX,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(divA3D(:,eR),Njkr, &
        &  MPI_DOUBLE_COMPLEX,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(cR_hp(eR,:),Njk, &
        &  MPI_DOUBLE_PRECISION,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(cZ_hp(eR,:),Njk, &
        &  MPI_DOUBLE_PRECISION,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(phiv_hp(eR,:),Njk, &
        &  MPI_DOUBLE_PRECISION,ip,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(vA_fe(:,:,eR),NDIM*nwidth*mnAFTot, &
        &  MPI_DOUBLE_COMPLEX,ip,0,MPI_COMM_WORLD,status,ierr)

    end if

  end if

  enddo ! eR


  if(me.eq.0) then

   do eR = 1, NelR
      do i = 2, NbDistr!-1
       if (eR>NelR-BorderSurf) then
        if (BorderMethod==2) then
         k_par_it(i,eR,:) = k_par_it(i,eR-1,:)
         k_perp(eR,:) = k_perp(eR-1,:)
        elseif (BorderMethod==3) then
         k_par_it(i,eR,:) = (1e-10+(NelR-eR)/real(BorderSurf))*k_par_it(i,NelR-BorderSurf,:)
         k_perp(eR,:) = (1e-10+(NelR-eR)/real(BorderSurf))*k_perp(NelR-BorderSurf,:)
        end if
       elseif (eR<BorderSurf) then
          k_perp(eR,:) = k_perp(BorderSurf,:)
       end if
      enddo

      if (OutputSelect(8)==1) then
       do itpr = 1, Njkr
        itp = itp_itpr_lkp(itpr)
        gR  = gR_itpr_lkp(itpr)
        if (gR == NGaussR/2+1) then
         k_par_hp(:,eR,itp) = k_par_it(:,eR,itpr)
        end if
       enddo
      end if
   enddo

  end if

  deallocate(Vls_Mat)
  deallocate(Vls_RHS)
  deallocate(IPIV)
  deallocate(k_par_new)
  deallocate(k_par_tmp)
  deallocate(k_cmp_tmp)
  deallocate(k_perp_rad,k_perp_bin)


  call MPI_BCAST(k_par_it, NbDistr*NelR*Njkr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  do eR = 1, NelR
    Ppla(:,eR+1) = Ppla(:,eR+1) + Ppla(:,eR)
    Pant(eR+1) = Pant(eR+1) + Pant(eR)
    Pdiv(eR+1) = Pdiv(eR+1) + Pdiv(eR)
!    PAdivA(eR) = PAdivA(eR) + PAdivA(eR-1)
  enddo

  if (m_npar==1) then
    call MPI_REDUCE(Ppla,PplaSum,AP_Nb*Nmesh,MPI_DOUBLE_COMPLEX,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(Pant,PantSum,Nmesh,MPI_DOUBLE_COMPLEX,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(Ppoynt,PpoyntSum,NelR,MPI_DOUBLE_COMPLEX,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(PAdivA,PAdivASum,NelR,MPI_DOUBLE_COMPLEX,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(Pdiv,PdivSum,Nmesh,MPI_DOUBLE_COMPLEX,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(TotalVolume,TotalVolumeSum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(divAintR,divAintRSum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(avrgA,avrgASum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(avrgEdgeA,avrgEdgeASum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,sep_comm(1),ierr)
    call MPI_REDUCE(maxFEqAmp,maxFEqAmpMax,mnEqTot,MPI_DOUBLE_PRECISION,MPI_MAX,0,sep_comm(1),ierr)
  end if


  if (me.eq.0) then

  divAint(:,:) = divAint(:,:)/Njk
  divAintRSum     = abs(divAintRSum) * (2*pi*2*pi/Njk)/nperiod
  TotalVolumeSum  = abs(TotalVolumeSum) * (2*pi*2*pi/Njk)/nperiod
  avrgASum        = abs(avrgASum) * (2*pi*2*pi/Njk)/nperiod
  avrgEdgeASum    = abs(avrgEdgeASum) * (2*pi*2*pi/Njk)/nperiod
  PdivCum = PdivSum(NMesh)
!  PdivCum = PdivCum * (2*pi*2*pi/Njk)/nperiod

  avrgFluxDiff = 0
  do eR = 1, NelR
    avrgFluxDiff = avrgFluxDiff + &
        abs( (PplaSum(1,eR)+PplaSum(1,eR+1))/2 - (PantSum(eR)+PantSum(eR+1))/2 - PpoyntSum(eR) )*(s_nod(eR+1)-s_nod(eR))
  enddo
  avrgFluxDiff = avrgFluxDiff / NelR

  maxFluxDiff = 0
  do eR = 1, NelR
    maxFluxTmp = abs( (PplaSum(1,eR)+PplaSum(1,eR+1))/2 - (PantSum(eR)+PantSum(eR+1))/2 - PpoyntSum(eR) )
    if (maxFluxTmp > maxFluxDiff)  maxFluxDiff = maxFluxTmp
  enddo

  print *, 'Global power balance'
  print *, 'Ppla(1)     = ', PplaSum(1,Nmesh)
  print *, 'Pant(1)     = ', PantSum(Nmesh)
  print *, 'diff        = ', PplaSum(1,Nmesh) - PantSum(Nmesh)
  print *, 'divCum      = ', PdivCum
  print *, 'AbsDivCum   = ', divAintRSum
  print *, 'diff+div    = ', PplaSum(1,Nmesh) - PantSum(Nmesh) + PdivCum
  print *, 'avrg.Div    = ', divAintRSum/TotalVolumeSum
  print *, 'avrg.A      = ', avrgASum/TotalVolumeSum
  print *, 'avrg.EdgeA  = ', avrgEdgeASum/TotalVolumeSum
  print *, 'av.Flux.Dif = ', avrgFluxDiff
  print *, 'maxFlux.Dif = ', maxFluxDiff
  print *, 'Tot. volume = ', TotalVolumeSum
  print *, ' '

  if (AP_nb>1) then
    print *, 'Power deposition by species:'
  end if

  i = 2
  if (AbsPower(1)==1) then
    print *, '     Electrons:', aimag(PplaSum(i,Nmesh))
    i = i+1
  end if
  if (AbsPower(2)==1) then
    print *, ' First species:', aimag(PplaSum(i,Nmesh))
    i = i+1
  end if
  if (AbsPower(3)==1) then
    print *, 'Second species:', aimag(PplaSum(i,Nmesh))
    i = i+1
  end if
  if (AbsPower(4)==1) then
    print *, '     Fast ions:', aimag(PplaSum(i,Nmesh))
    i = i+1
  end if

  if(AP_nb>1) print *, ' '


  open(fn, file='PowersDiff.txt')
  write (fn, *) Nmesh

  do eR = 1, Nmesh
    write(fn,'(1x, 9e25.16)') s_nod(eR), real(PplaSum(1,eR)), imag(PplaSum(1,eR)), &
     real(PantSum(eR)), imag(PantSum(eR)), real(PplaSum(1,eR)+PantSum(eR)), imag(PplaSum(1,eR)+PantSum(eR)), &
     real(PdivSum(eR)), imag(PdivSum(eR))
  enddo

  write (fn, *) NelR
  do eR = 1, NelR
    write(fn,'(1x, 3e25.16)') s_hp(eR), real(PpoyntSum(eR)), imag(PpoyntSum(eR))
  enddo

  write (fn, *) NelR
  do eR = 1, NelR
    write(fn,'(1x, 3e25.16)') s_hp(eR), real(PAdivASum(eR)), imag(PAdivASum(eR))
  enddo


  write (fn, *) NelR
  gR = NGaussR / 2 + 1
  do eR = 1, NelR
    write(fn,'(1x, 2e25.16)') s_gauss(eR,gR), divAint(eR,gR)
  enddo

  write (fn, *) AP_nb-1
  do i = 1, 4
  if (AbsPower(i)==1) write(fn, *) i
  enddo

  do i = 2, AP_nb
    do eR = 1, Nmesh
      write(fn,'(1x, 2e25.16)') real(PplaSum(i,eR)), aimag(PplaSum(i,eR))
    enddo
  enddo

  close(fn)


  open(fn, file='convergence.txt')
  write (fn, *) 'mAFTot, nAFTot, NelR'
  write (fn, *) mnAFTot, NelR
  write (fn, *) 'PplaSum(1), Ppla(1) + Pant(1), divCum'
  write (fn, '(1x, 6e19.10)') real(PplaSum(1,Nmesh)),               imag(PplaSum(1,Nmesh))              , & 
                              real(PplaSum(1,Nmesh) + PantSum(Nmesh)), imag(PplaSum(1,Nmesh) + PantSum(Nmesh)), &
                              real(PdivCum),                   imag(PdivCum)
  write (fn, *) 'avrg. div,  avrgASum,  avrgEdgeA, tot. volume'
  write (fn, '(1x, 4e19.10)') divAintRSum/TotalVolumeSum, avrgASum/TotalVolumeSum, &
                              avrgEdgeASum/TotalVolumeSum, TotalVolumeSum
  write (fn, *) 'avrg. and max flux diff.'
  write (fn, '(1x, 2e19.10)') avrgFluxDiff, maxFluxDiff
  write (fn, *) 'PplaSum+Pant+PdivCum'
  write (fn, '(1x, 2e19.10)') real(PplaSum(1,Nmesh) + PantSum(Nmesh) + PdivCum),  &
                              imag(PplaSum(1,Nmesh) + PantSum(Nmesh) + PdivCum)
  close(fn)


  open(fn, file=Trim(OUTPUTPATH)//'Ppla3D.dat',form='unformatted')
  write(fn) NelR, njT, nkT
  write(fn) s_hp
  write(fn) cR_hp
  write(fn) cZ_hp
  write(fn) OutputSelect
  write(fn) NbDistr
  do i = 1, 7
    if (k_par_lkp(i)==1) write(fn) i+1
  enddo
  write(fn) AP_nb-1
  do i = 1, 4
    if (AbsPower(i)==1) write(fn) i
  enddo
  write(fn) real(sum(PplaSum(2:AP_nb,Nmesh)))
  write(fn) imag(sum(PplaSum(2:AP_nb,Nmesh)))
  if (AbsPower(3)==1) write(fn) imag(PplaSum(AP_nb-1,Nmesh))
  if (AbsPower(4)==1) write(fn) imag(PplaSum(AP_nb,Nmesh))

  if (OutputSelect(1)==1) then
    write(fn) real(E3D_hp)
    write(fn) imag(E3D_hp)
  end if

  if (OutputSelect(2)==1) then
    write(fn) real(B3D_hp)
    write(fn) imag(B3D_hp)
  end if

  if (OutputSelect(3)==1) then
    write(fn) real(A3D_hp)
    write(fn) imag(A3D_hp)
  end if

  if (OutputSelect(4)==1) then
    write(fn) real(Ppla3D(1,:,:))
    write(fn) imag(Ppla3D(1,:,:))
    do i = 2, AP_nb
      write(fn) real(Ppla3D(i,:,:))
      write(fn) imag(Ppla3D(i,:,:))
    enddo
  end if

  do i = 1, sum(OutputSelect(5:7))
    write(fn) real(t_hp(i,:,:))
    write(fn) imag(t_hp(i,:,:))   
  enddo

  if (OutputSelect(8)==1) then
    write(fn) n_ant
    do i = 1, NbDistr
      write(fn) sign(1,n_ant)*k_par_hp(i,:,:)
    enddo
    write(fn) k_perp
  end if
 
  close(fn)

  end if

  deallocate(k_perp)

end subroutine
!-------------------------------------------
 subroutine CoulombGauge
  integer eR,gR,itp,itpr, j
  integer, parameter :: fn = 2

  
  gR = NGaussR/2+1
  open(fn, file='divA.dat')
  write (fn, *) NelR

  itp = 1
  itpr = (gR-1)*Njk + itp
  do eR = 1, NelR
    write(fn,'(1x, 3e17.10)') s_hp(eR), real(divA3D(itpr,eR)),  imag(divA3D(itpr,eR))
  enddo

  close(fn)


 end subroutine
!-------------------------------------------
 subroutine SaveAntenna
  integer, parameter :: fn = 2
  integer  eR, gR, itpr, itp

  do eR = 1, NelR

    if (CylMetric) then
      call Init_Cyl_Equilibruim(eR)
    else
      call Init_Equilibruim(eR)
    endif

    call Init_jext3D(eR)

    do itpr = 1, Njkr
      gR = gR_itpr_lkp(itpr)
      itp = itp_itpr_lkp(itpr)
    
      if (gR == NGaussR/2+1) then
        A3D_hp(1,eR,itp) = jext_3D(itpr,1)
        A3D_hp(2,eR,itp) = jext_3D(itpr,2)
        A3D_hp(3,eR,itp) = jext_3D(itpr,3)
        A3D_hp(4,eR,itp) = 0
      endif
    enddo

  enddo

  open(fn, file=Trim(OUTPUTPATH)//'Antenna3D.dat',form='unformatted')
  write(fn) NelR, njT, nkT
  write(fn) cR_hp
  write(fn) cZ_hp
  write(fn) real(A3D_hp)
  write(fn) imag(A3D_hp)
  close(fn)

 end subroutine
!-------------------------------------------

 end module fields
