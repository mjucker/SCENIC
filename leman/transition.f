!
!  Management of the memory
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!==========================================================================
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


 module transition

 use math
 use matrix
 use fields

 implicit none


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 contains
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


!-------------------------------------------
! Memory allocation before the matrix construction
 subroutine SolveStart
  integer i,itp

  IF( .NOT. ALLOCATED(Acoeff) ) allocate(Acoeff(Njk,NDIM,NDIM,(DIM+1)*(DIM+1)))
  IF( .NOT. ALLOCATED(Acoeff_k_f_Dmult_imn) ) allocate(Acoeff_k_f_Dmult_imn(mnAFtot*mnAFTot,(DIM+1)*(DIM+1)))
  IF( .NOT. ALLOCATED(Acoeff_k_f) ) allocate(Acoeff_k_f(mnEqTot,(DIM+1)*(DIM+1)))

  IF( .NOT. ALLOCATED(KKA) ) allocate(KKA(LBKK,LBKK))
  IF( .NOT. ALLOCATED(KKB) ) allocate(KKB(LBKK,LBKK,indexfile))
  IF( .NOT. ALLOCATED(KKC) ) allocate(KKC(LBKK,LBKK/same))
  IF( .NOT. ALLOCATED(KKD) ) allocate(KKD(LBKK,LBKK))
  IF( .NOT. ALLOCATED(KKE) ) allocate(KKE(LBKK,LBKK))

  IF( .NOT. ALLOCATED(KK_tmp) ) allocate(KK_tmp(LBKK,LBKK))
  IF( .NOT. ALLOCATED(KK_mult) ) allocate(KK_mult(LBKK,LBKK/same))
  IF( .NOT. ALLOCATED(FF_tmp) ) allocate(FF_tmp(LBKK))

  KKB = 0._dp
  KKD = 0._dp
  FTransformError = 0
  maxFEqAmp = 0
  maxEqValue = 0
  XXc = 0._dp

  if (FourierCheck) then
    IF( .NOT. ALLOCATED(EqExp) ) allocate(EqExp(mnEqTot, Njk))
    do i = 1, mnEqTot
      do itp = 1, Njk
        EqExp(i, itp) = Exp(-iim*(mEq(i)*th_itp(itp)+nEq(i)*ph_itp(itp)))
      enddo
    enddo
  end if

  if (npsep==1) then
    gmin = 1
    gmax = NelR
    incr = 1
  elseif (npsep==2) then
    if (role==0) then
      gmin = 1
      gmax = NelR/2
      incr = 1
    else
      gmin = NelR
      gmax = NelR/2+1
      incr = -1
    end if
  else
    if (role==0) then
      gmin = 1
      gmax = NelR/2
      incr = 2
    elseif (role==1) then
      gmin = 2
      gmax = NelR/2
      incr = 2
    elseif (role==2) then
      gmin = NelR-1
      gmax = NelR/2+1
      incr = -2
    else
      gmin = NelR
      gmax = NelR/2+1
      incr = -2
    end if
  end if

 end subroutine

!-------------------------------------------
! Memory allocation between the computation and the treatment of the results
 subroutine SolveToTreatment

  deallocate(Acoeff_k_f_Dmult_imn)
  deallocate(Acoeff_k_f)

  deallocate(KK_tmp)
  deallocate(KK_mult)
  deallocate(FF_tmp)

  if (FourierCheck) deallocate(EqExp)

  allocate(E3D(DIM,Njkr))
  allocate(B3D(DIM,Njkr))
  allocate(GradPhi3D(DIM,Njkr))

  allocate(vdAph_dui_g(3,mnAFTot,NelR,NGaussR,3,nwidth/np_npar+1))
  allocate(Ppla3D(AP_nb,NelR,Njk))
  allocate(E3D_hp(DIM,NelR,Njk))
  allocate(B3D_hp(DIM,NelR,Njk))
  allocate(A3D_hp(NDIM,NelR,Njk))
  allocate(divA3D(Njkr,NelR))
  allocate(t_hp(sum(OutputSelect(5:7)),NelR,Njk))
  allocate(k_par_hp(NbDistr,NelR,Njk))

  allocate(scalar_temp(Njkr))
  allocate(vector_temp(3,Njkr))
  allocate(quadri_temp(4,Njkr))
  allocate(real_temp(Njkr))

  allocate(Ppla(AP_nb,Nmesh))
  allocate(PplaSum(AP_nb,Nmesh))
  allocate(PplaCum(AP_nb))
  allocate(EeE(AP_nb))

  allocate(vA_hp(NDIM, mnAFTot, NelR))
  allocate(vAp_hp(NDIM, mnAFTot, NelR))
  allocate(vA_g(NDIM, mnAFTot, NelR, NGaussR, nwidth/np_npar+1))
  allocate(vAp_g(NDIM, mnAFTot, NelR, NGaussR, nwidth/np_npar+1))


 end subroutine

!-------------------------------------------
! Memory deallocation at the end of the treatment of the solution
 subroutine FinishTreatment

    deallocate(E3D)
    deallocate(B3D)
    deallocate(GradPhi3D)

    deallocate(Ppla3D)
    deallocate(E3D_hp)
    deallocate(B3D_hp)
    deallocate(A3D_hp)
    deallocate(divA3D)
    deallocate(t_hp)
    deallocate(k_par_hp)

    deallocate(vdAph_dui_g)

    deallocate(vA_hp);
    deallocate(vAp_hp);
    deallocate(vA_g);
    deallocate(vAp_g);

    deallocate(scalar_temp)
    deallocate(vector_temp)
    deallocate(quadri_temp)
    deallocate(real_temp)

    deallocate(Ppla)
    deallocate(PplaSum)
    deallocate(PplaCum)
    deallocate(EeE)

 end subroutine
!-------------------------------------------
! Final memory deallocation
 subroutine DeallocateAll


  deallocate(Acoeff)

  deallocate(KKA)
  deallocate(KKB)
  deallocate(KKC)
  deallocate(KKD,KKE)

  deallocate(mEq,nEq)
  deallocate(mEqDT,nEqDT)
  deallocate(maxFEqAmp,maxFEqAmpMax)
  deallocate(mnFFT_lkp,mnFFT_DT_lkp)

  deallocate(om_a)
  deallocate(GlobalIndex_lkp)
  deallocate(Psi_FE)

  deallocate(mn_mnmnp_lkp,mnp_mnmnp_lkp)
  deallocate(mn_mnmnpDT_lkp,mnp_mnmnpDT_lkp)

  deallocate(mneq0_lkp,mneq0DT_lkp)
  deallocate(mnFFT_RHS_lkp)

  deallocate(iF_lkp,iA_lkp)

  deallocate(k_par_it)

  if(UseDensCoeffs==3.or.UseTempCoeffs==3)deallocate(nT_fits)
  if(UseDensCoeffs==4.or.UseTempCoeffs==4)deallocate(thermal)
  if(UseHotDensCfs==3.or.UseHotTempCfs==3)deallocate(moments)

  deallocate(vA_fe)
  deallocate(vAp_fe)

  deallocate(tTU)
  deallocate(tTL)
  deallocate(tbjac)
  deallocate(tbmod)
  deallocate(trx)
  deallocate(trz)
  deallocate(tphiv)
  deallocate(tgssu)
 
 
  deallocate(mTL)
  deallocate(mTU)
  deallocate(mTLs)
  deallocate(mTUs)
  deallocate(mTLt)
  deallocate(mTUt)
  deallocate(mTLp)
  deallocate(mTUp)

  deallocate(mgssu)
  deallocate(bjac)
  deallocate(bjacs)
  deallocate(bjact)
  deallocate(bjacp)
  deallocate(bmod)
  deallocate(Response)

  deallocate(cR)
  deallocate(cZ)
  deallocate(phiv)
  deallocate(cR_hp)
  deallocate(cZ_hp)
  deallocate(phiv_hp)

  deallocate(mAF)
  deallocate(nAF)
  deallocate(mpAF)
  deallocate(npAF)

  deallocate(Exp_mnp)
  deallocate(Exp_mn)

  deallocate(t11m);
  deallocate(t12m);
  deallocate(t21m);
  deallocate(t22m);
  deallocate(t33m);

  deallocate(iAF_lkp);

  deallocate(emin);
  deallocate(emax);
  deallocate(FF);
  deallocate(XXc);


!  deallocate(tTUy2);
!  deallocate(tTLy2);
!  deallocate(tbjacy2);
!  deallocate(tbmody2);

  deallocate(tbjacTO);
  deallocate(tbmodTO);
  deallocate(tgssuTO)
  deallocate(tgpplTO)
  deallocate(tgtplTO)
  deallocate(tgttlTO)
  deallocate(tgstlTO)
  deallocate(tgsslTO)
  deallocate(tgsplTO)
  deallocate(tb2TO)
  deallocate(trxTO)
  deallocate(trzTO)
  deallocate(tphivTO)

 end subroutine

 end module
