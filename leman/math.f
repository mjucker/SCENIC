!  Waves  3D  
!
! Main constants definition, functions initialization, gauss coefficients etc.
!
!
 module math

 implicit none
 include 'mpif.h'
   
 integer, parameter :: dp = SELECTED_REAL_KIND(10)


! ---------------------------------------------------------------------------

 integer, parameter ::  NelR    =150


 integer, parameter ::  NGaussR =  5     ! NGaussR -- gauss points/interval (radial direction)
 ! Attention: NGaussR should be odd so that the middle gauss points
 ! coincide with half-mesh points (used in field calculations and diagnostics)

! Data directly from TERPSICHORE - number of radial, poloidal & toroidal cells
 integer, parameter :: niTO =  96,     & ! number of surfaces in TERPSICHORE equilibrium (T.eq.)
                       njTO = 100,     & ! number of poloidal nodes (cells) ! Only even!
                       nkTO =  2,      & ! number of toroidal nodes (cells) ! Only even!
                       NjkTO = njTO*nkTO ! total number of poloidal & toroidal cells (intervals) = nj*nk (TERPSICHORE)



! Interpolated on the new angular mesh equilibrium:
! needed more precision than in T. equilibrium to correctly describe the F. modes.
 integer, parameter :: niT  = niTO,    &
                       njT  = 450,     & ! number of poloidal nodes (cells)
                       nkT  =   2,     & ! number of toroidal nodes (cells)  ! even!
                       Njk  = njT*nkT, & ! total number of poloidal & toroidal cells (intervals) = nj*nk
                       Njkr = Njk * NGaussR
                       ! j <-> theta,  k <-> phi



 ! Perturbations
 integer mnAFTot, mnDTTot
 integer, dimension(:), allocatable   :: mAF    ! m mode lookup table, size: 1..mnAFTot
 integer, dimension(:,:), allocatable :: nAF    ! n mode lookup table, size: 1..mnAFTot
 integer, dimension(:), allocatable   :: mpAF   ! m mode lookup table for test functions
 integer, dimension(:,:), allocatable :: npAF   ! n mode lookup table for test functions
 integer, dimension(:), allocatable   :: mDT
 integer, dimension(:,:), allocatable   :: nDT
 integer, dimension(:), allocatable   :: mpDT
 integer, dimension(:,:), allocatable   :: npDT

! ---------------------------------------------------------------



! ---------------------------------------------------------------------------

 integer, parameter :: fnInput = 5  ! Standard input (for the Namelist)
 
 character (len=*), parameter :: EQFILE = 'fort.37.DAT.N' ! Equilibrium (produced with eqtransf.f)
 character (len=*), parameter :: AFTableFILE = 'AFtable.txt' ! "A" Fourier mode table
 character (len=*), parameter :: DTTableFILE = 'DTtable.txt' ! Dielectric Tensor Fourier mode table
 character (len=*), parameter :: STORAGEPATH = '/scratch/jucker/' ! Storage of the matrix elements

 character (30)   eqname        ! Equilibrium name

 character (200)  EQPATH,     & ! Path to the equilibrium file (EQFILE)
                  GRIDPATH,   & ! Path to the grid file ('grid.dat')
                  OUTPUTPATH    ! Path for the output

!  Mesh
 integer, parameter ::  &  !  NelR = above,   & ! number or radial FE intervals
                        Nnod = 2,   & ! number of nodes per element
                        Nbf  = 4,   & ! number of basis functions per element (2 for LFE, 4 for cubics)
                        Nmesh = NelR + 1,   & ! total number of nodes
                        Nunknowns = 2*(NelR-1) + 4 ! total number of unknowns (coefficients) (cubics)
 

 integer, parameter ::  &  !  NGaussR = above, & ! NGaussR -- gauss points/interval (radial direction), 
                       MaxNGauss = 6

 
 integer, parameter :: NDIM = 4   ! Variables: A1 <-> Ar,  A2 <-> Atheta,  A3 <-> Aphi,  A4 <-> phi
 integer, parameter :: DIM  = 3   ! Coordinates:  1 <-> r,  2 <-> theta,  3 <-> phi
 integer, parameter :: xxx = 1,  yyy = 2,  zzz = 3


! ---------------------------------------------------------------------------
! EQUILIBRIUM

 integer  mnEqTot,mnEqDTTot ! number of all equilibrium modes

 integer, dimension(:), allocatable   :: mEq ! list of all equilibrium modes
 integer, dimension(:), allocatable   :: nEq
 integer, dimension(:), allocatable   :: mEqDT
 integer, dimension(:), allocatable   :: nEqDT

 real(dp) maxEqValue
 real(dp), dimension(:), allocatable  :: maxFEqAmp, maxFEqAmpMax
 integer  mEqMax, nEqMax

 integer, dimension(:), allocatable :: mneq0_lkp, mneq0DT_lkp
 ! Lookup table for equilibrium modes that contribute to the variational form
 ! for given m,n,mp,np. Returns the index in the mEq, nEq table

 integer, dimension(:), allocatable :: mnFFT_lkp, mnFFT_DT_lkp
 integer, dimension(:), allocatable :: mnFFT_RHS_lkp
 ! Lookup table for equilibrium modes. Returns the index in the FFT mode
 ! table for a given index in the mEq,nEq table

 integer, dimension(:), allocatable :: mn_mnDT_lkp

 complex(dp), dimension(:,:), allocatable ::  EqExp

 integer, dimension(:,:,:), allocatable :: iAF_lkp ! lookup table: selects the type of basis
                                            ! functions according to the mode number
 integer, dimension(:,:), allocatable  :: iF_lkp ! same for RHS calculation
 integer, dimension(:,:), allocatable  :: iA_lkp ! ... and for the potentials in real space

 integer, parameter :: MaxEdgeModeNo = 100
 integer                  EdgeModeNo
 integer EdgeModes_lkp(MaxEdgeModeNo) ! lookup table: indices of the F. modes
                                      ! located on the edge of the Fourier box.
                                      ! Used for the convergence study
! ---------------------------------------------------------------------------
! FFTW

 integer*8 :: forw, backw, trans

! ---------------------------------------------------------------------------

 real(dp) pi
 complex(dp) iim

 
!  Geometry
 real(dp)  bx0,  bx1,  by0,  by1, bz0,  bz1  ! boundary (r, theta, phi)
 real(dp)  MajRadius, B02
 real(dp) ftp(niT),fpp(niT),ci(niT),cj(niT),cip(niT),cjp(niT)


! -----------------------------------------
 ! Equilibrium metric (TU,TL,bjac,bmod, derivatives) on the gauss grid for only
 ! one radial element (FE). The values are calculated "on the fly" by interpolation
 ! from the TERPSICHORE grid (tTU,tTL), on the each radial loop step.
 ! dimensions: (:,:,:,:) <-> ii,jj,NGaussR,njk (ii,jj -- transform matrix indices)
 real(dp), dimension(:,:,:,:), allocatable :: mTL,mTU,mTLs,mTUs,mTLt,mTUt,mTLp,mTUp
 real(dp) mftp(NGaussR),mfpp(NGaussR),mci(NGaussR),mcj(NGaussR)
 ! dimensions: (:,:) <-> NGaussR,njk
 real(dp), dimension(:,:), allocatable :: bjac,bjacs,bjact,bjacp,bmod,mgssu, &
                                     cR,cZ, cR_hp, cZ_hp, phiv, phiv_hp

   !       bjac            Jacobian(s,theta) (Boozer)
   !       bjacs           d/ds(Jacobian(s,theta))
   !       bjact           d/dt(Jacobian(s,theta))
   !       bmod            B(s,theta) <- normalized to 1.0 on the axis!
   !       cR              cR(s,theta) r in cylindrical coord.
   !       cZ              cZ(s,theta) z in cylindrical coord.
   !       cR/Z -- on one radial element, cR/Z_hp -- all elements, half-points




 ! Equilibrium directly from TERPSICHORE, before integpolation on the new abgular grid
 real(dp), dimension(:,:), allocatable :: tgssuTO, tgpplTO, tgtplTO, tgttlTO, tgstlTO, &
                                          tgsslTO, tgsplTO, trxTO, trzTO, tb2TO,       &
                                          tbjacTO, tbmodTO, tphivTO



 ! Same matrices (equilibrium), new angular grid grid, all radial surfaces stored simultaneously
 ! dimensions: (:,:,:,:) <-> ii,jj,njk,ni (ii,jj -- transform matrix indices)
 real(dp), dimension(:,:,:,:), allocatable :: tTU, tTL  !, tTUy2, tTLy2


 real(dp), dimension(:,:), allocatable :: tbjac, tbmod !, tbjacy2, tbmody2
! real(dp)  vect_tmp(niT), vect_tmpy2(niT)

 ! temporary arrays; are only used to calculate tTL, tTU matrices
 real(dp), dimension(:,:), allocatable :: tgssu, tgppl, tgtpl, tgttl, tgstl, &
                                          tgssl, tgspl, trx, trz, tb2, tphiv



 real(dp) MAxisR(nkT),MAxisZ(nkT)


 real(dp) s_T(niT), iota(niT), iaspct(niT), &
          theta_T(0:njT+1), phi_T(0:nkT+1),   &
          th_itp(Njk), & ! "lookup table" for poloidal angle
          ph_itp(Njk)    ! "lookup table" for toroidal angle

 integer itp_itpr_lkp(Njkr), it_itpr_lkp(Njkr), ip_itpr_lkp(Njkr), &
         gR_itpr_lkp(Njkr), itpr_itp_gR_lkp(Njk,NGaussR)

  ! Dielectric tensor components
  ! The tensor is stored globally, in all the grid points because it is not only
  ! used in the variational form construction, but also for power flux calculations (fields.f)
 complex(dp), dimension(:,:,:), allocatable ::  t11m, t12m, t21m, t22m, t33m



! -----------------------------------------
 integer, parameter :: NoVarTerms  = 8,  & ! number of additive terms in the variational form
                       NoSPTerms   = 3,  & ! number of additive terms in the scalar products
                       NoAllTerms  = NoVarTerms*NoSPTerms,  & ! total number of additive terms
                       
   ! number of variables
   !(3 vector potential components + scalar potential) * 4 - no derivation + d/ds,d/dtheta,d/dphi
                       NoVariables = 16

 integer  itp_it_ip_lkp(njT,nkT)

! -----------------------------------------


 ! NGaussR -- gauss points/interval (radial), NGaussT -- gauss points/interval (theta) 
 real(dp) xgauss(2:MaxNGauss, 1:MaxNGauss), wgauss(2:MaxNGauss, 1:MaxNGauss)
 real(dp) s_nod(1:Nmesh),              & ! FE radial nodes, integer points
          s_gauss(1:NelR, 1:NGaussR),  & ! radial gauss points
          s_hp(NelR)                     ! half-points (FE) radial mesh


 ! Multiplying factor ( ksi_i = A_i * sf_v_i )
 real(dp) sf_v(4,NelR,NGaussR), dsf_v(4,NelR,NGaussR)  ! values in gauss points
 real(dp) sf_v_A(4,Nmesh),      dsf_v_A(4,Nmesh)       ! values in FE grid points


 ! Multiplying factor, equilibrium ( TU/L ->  TU/L *  se_TU/L 
 !                                   1 <-> (1,1),  2 <-> (2,2)  )
 real(dp) seTU(2,NelR,NGaussR), dseTU(2,NelR,NGaussR)  ! values in gauss points
 real(dp) seTL(2,NelR,NGaussR), dseTL(2,NelR,NGaussR)
 real(dp) setTU(2,niT), dsetTU(2,niT)  ! values in TERPSICHORE radial points
 real(dp) setTL(2,niT), dsetTL(2,niT)


 integer i_ij_lkp((DIM+1)*(DIM+1)) /0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3/
 integer j_ij_lkp((DIM+1)*(DIM+1)) /0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3/

! integer da_i_lkp(0:DIM) / 0,1,0,0 /
 integer da_ij_lkp((DIM+1)*(DIM+1)) /0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0/
 integer df_ij_lkp((DIM+1)*(DIM+1)) /0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0/

! Lookup tables for matrix contruction
 integer, dimension(:), allocatable :: mn_mnmnp_lkp, mnp_mnmnp_lkp
 integer, dimension(:), allocatable :: mn_mnmnpDT_lkp, mnp_mnmnpDT_lkp

 integer ab_lkp(NDIM,NDIM), a_ab_lkp(Nbf*Nbf), b_ab_lkp(Nbf*Nbf)
 integer, dimension(:,:,:,:), allocatable :: GlobalIndex_lkp
 integer CFCAINFO(4,73)


 ! External antenna currents
 complex(dp) jext_3D(Njkr,NDIM)
 complex(dp), dimension(:,:,:), allocatable :: Exp_mnp ! table of Exp() for right-hand side integration
 complex(dp), dimension(:,:,:), allocatable :: Exp_mn ! table of Exp() for transformation to real space (potentials)

 complex(dp), dimension(:,:,:), allocatable :: Dmult_imn   ! mnmnp,i,j
          ! Multiplication factor (comes from the derivative of exponential (p.e. exp(i m theta)))
          ! i,j -- derivative index (0-no derivative, 1-d/ds, 2-d/dtheta, 3-d/dphi)

 real(dp) dPsiRAF(NGaussR,0:3,(DIM+1)*(DIM+1),1:Nbf,1:Nbf) 
          ! A & F -- potentials & test function
          ! derivative index (0-no derivative, 1-d/ds, 2-d/dtheta, 3-d/dphi)
          ! 1:4 -- AModified, FModified: 
          !     0: (A=0,F=0),   1: (A=0,F=1),   2: (A=1,F=0),   3: (A=1,F=1)


 !               da    a    e   s_gauss(eR,gR) Modified/n
 real(dp), dimension(:,:,:,:,:,:), allocatable :: Psi_FE


 integer VDATo16(0:DIM, NDIM) ! converts derivative index (0..DIM) and A component numb. (1..NDIM)
                              ! to variable number (1..16) 
 integer LMat(1:Nbf, 1:NelR+1)

!  ----------- Plasma parameters ----------- !
 real(dp) c_SI,mp_SI,me_SI,mi1_SI,mi2_SI,mih_SI,iA1,iA2,iAh,ifr1,ifr2, & ! constants
          iZ1,iZ2,iZh,n0_SI,ne_SI,Ti0,Te0,OmCi1_SI,OmCi2_SI,cA_SI, &
          mu0_SI,qe_Gauss,CCa,Ra,Rw,B0,om_im,om_imi,om_SI,a_plasma, &
          s_ant_min,s_ant_max,theta_ant1,theta_ant2,theta_ant3,theta_ant4
 integer  m_ant, n_ant, nperiod, ant_type, m_shape, n_shape
 real(dp) theta_antB
 real(dp) Th0, nh0_SI, An0

 ! Plasma density and temperature profile coefficients
 real(dp) np0, np1, np2, np3, np4, np5
 real(dp) nh0, nh1, nh2, nh3, nh4, nh5
 real(dp) Tip0, Tip1, Tip2, Tip3, Tip4, Tip5
 real(dp) Tep0, Tep1, Tep2, Tep3, Tep4, Tep5
 real(dp) ArtBorder
 real(dp) Thp0, Thp1, Thp2, Thp3, Thp4, Thp5
 real(dp) Anp0, Anp1, Anp2, Anp3, Anp4, Anp5
 integer  UseDensCoeffs
 integer  UseTempCoeffs
 integer  UseHotDensCfs
 integer  UseHotTempCfs

 real(dp) nu_f(1:NelR,1:NgaussR) ! Absorption profile in cold plasma  (values in the gauss points only)
 real(dp) fnp(1:NelR,1:NgaussR)  ! Plasma density  n(r)/n(0)  (values in the gauss points only)
 real(dp) fnm(1:NelR,1:NgaussR)  ! Minority density  n(r)/n(0)  (values in the gauss points only)
 real(dp) fnhp(1:NelR,1:NgaussR) ! Fast ions density  n(r)/n(0)  (values in the gauss points only)
 real(dp) fti(1:NelR,1:NgaussR)  ! Ion temperature Ti(r)/Ti(0)  (values in the gauss points only)
 real(dp) fte(1:NelR,1:NgaussR)  ! Electron temperature Te(r)/Te(0)  (values in the gauss points only)
 real(dp) ftm(1:NelR,1:NgaussR)  ! Minority temperature Tm(r)/Ti(0)  (values in the gauss points only)
 real(dp) fth(1:NelR,1:NgaussR)  ! Fast Ions parallel temperature Th(r)/Th(0)  (values in the gauss points only)
 real(dp) anis(1:NelR,1:NgaussR) ! Fast Ions anisotropy Thpar/Thperp  (values in the gauss points only)

 real(dp) dens_T(niT),Bc ! plasma density on the TERPSICHORE grid
 real(dp), dimension(:,:), allocatable ::  moments,thermal
!  ----------------------------------------- !
 real(dp), dimension(:,:), allocatable :: om_a  ! Frequency of the Alfven resonance in a cylindrical
                                ! approximation

 integer nif
 real(dp), dimension(:,:), allocatable :: nT_fits

 integer  OmSteps ! Number of steps for frequency scan (read from input)
 integer  NbIter  ! Number of iteration for a single frequency
 real(dp) om_SI1, om_SI2, om_SIv(1:100) ! min, max frequency, array of used frequencies values
 real(dp), allocatable :: Response(:)  ! plasma response (energy coupled in the antenna)
 

! Output switches
 logical OutputFlags(2)
   ! 1 - Output divA 3D
   ! 2 - Output EB   3D
 integer OutputSelect(8)

! Separated power deposition switches
 integer AbsPower(4)
   ! 1 - Electrons
   ! 2 - First species
   ! 3 - Second species
   ! 4 - Fast population
 integer AP_nb, AP_lkp(4)

! Density display parameters

 integer, parameter :: Gwidth  = 48
 integer, parameter :: Gheight = 16


 integer PrecalcMatrix, MaxResponseFreq, EmptyRun
 logical CylMetric, FourierCheck
 integer DielTensor, kparApprox, k_par_lkp(7), Nb_intv
 integer BorderMethod, BorderSurf, NbDistr
 real(dp) RRcyl, v_range, RelaxParam

 integer   meq_cyl,neq_cyl
 real(dp)  ampl_cyl,s1_cyl,s2_cyl

 integer qq, qiter


! variables for matrix.f

 integer  MSize,   &
          MSize2,  & ! number of rows of the matrix KK (total number of unknowns)
          MDiag      ! number of non-zero diagonals in KK
                     ! overlapping linear FE ('hat'): 3
                     ! overlapping Hermite cubics FE: 6
 ! parameters used in the solver
 integer  MBlock, KL, KU, LDAB, LBKK, LBKK_red

 complex(dp), dimension(:,:,:,:), allocatable :: Acoeff
 complex(dp), dimension(:,:), allocatable     :: Acoeff_k_f_Dmult_imn
 complex(dp), dimension(:,:), allocatable     :: Acoeff_k_f

 complex(dp), dimension(:),   allocatable :: FF     ! right-hand side vector (antenna)
 complex(dp), dimension(:,:,:), allocatable :: XX, XXc  ! solution (fourier harmonics stored separately)


! variables for fields.f

 ! Attention!  vA here is the potential (after multiplying by the factor "sfr")
 complex(dp), dimension(:,:,:), allocatable   :: vA_hp, vAp_hp ! fe mesh half points
 complex(dp), dimension(:,:,:,:,:), allocatable :: vA_g, vAp_g ! gauss points
 complex(dp), dimension(:,:,:), allocatable   :: vA_fe, vAp_fe ! nodes (FE)

 complex(dp), dimension(:,:), allocatable   :: E3D, B3D, GradPhi3D
 real(dp), dimension(:,:,:), allocatable   :: k_par_it

 complex(dp), dimension(:,:), allocatable   :: divA3D
 complex(dp), dimension(:,:,:), allocatable :: E3D_hp, B3D_hp, A3D_hp, Ppla3D
 complex(dp), dimension(:,:,:), allocatable :: t_hp
 real(dp), dimension(:,:,:), allocatable    :: k_par_hp


! Variables for MPI and OpenMp

 integer :: ierr, me, nprocs, role, same, posbuff, sgme, msep, npsep
 integer, parameter :: lbuff=2000
 character :: buff(lbuff)
 integer, dimension(:), allocatable :: emin, emax, sep_comm, nind_comm
 integer, dimension(:,:,:), allocatable :: me_lkp
 integer  :: gmin, gmax, bmin, bmax, incr
 integer  :: esup, indexfile
 real(dp) :: KKMem
 integer  :: ParaMethod, nindex, nwidth, m_npar, np_npar

 complex(dp), dimension(:,:), allocatable   :: KKA, KKC, KKD, KKE
 complex(dp), dimension(:,:,:), allocatable :: KKB

 complex(dp), dimension(:,:), allocatable   :: KK_tmp, KK_mult
 complex(dp), dimension(:), allocatable     :: FF_tmp


 namelist /InputNameList/ EQPATH,OUTPUTPATH, GRIDPATH, &
    Ra,Rw,mp_SI,me_SI,iA1,ifr1,iZ1,iA2,ifr2,iZ2,iAh,iZh,n0_SI,np0,np1,np2,np3, &
    np4,np5,nh0_SI,nh0,nh1,nh2,nh3,nh4,nh5,Ti0,Te0,Th0,Tip0,Tip1,Tip2,Tip3,Tip4, &
    Tip5,Tep0,Tep1,Tep2,Tep3,Tep4,Tep5,Thp0,Thp1,Thp2,Thp3,Thp4,Thp5,An0,Anp0,Anp1, &
    Anp2,Anp3,Anp4,Anp5,ArtBorder,UseTempCoeffs,UseDensCoeffs,UseHotDensCfs, &
    UseHotTempCfs,B0,Bc,Om_SI1,Om_SI2,OmSteps,NbIter,ant_type,n_shape,m_shape,n_ant, &
    m_ant,s_ant_min,s_ant_max,theta_ant1,theta_ant2,theta_ant3,theta_ant4,om_im, &
    om_imi,PrecalcMatrix,MaxResponseFreq,EmptyRun,OutputFlags,OutputSelect,CylMetric, &
    nperiod,meq_cyl,neq_cyl,ampl_cyl,s1_cyl,s2_cyl,RRcyl,DielTensor,kparApprox, &
    k_par_lkp,Nb_intv,v_range,BorderMethod,BorderSurf,RelaxParam,mEqMax,nEqMax, &
    FourierCheck,AbsPower,KKMem,ParaMethod,nwidth,np_npar

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 contains
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 

!-------------------------------------------
! Parallel initialization
subroutine Init_Para
  integer ip, ipa, ipb, ipc

  if (ParaMethod == 1) then
    nwidth  = 1
    np_npar = 1
  end if

  if (nprocs/=np_npar.and.nprocs/=2*np_npar.and.mod(nprocs,4*np_npar)/=0) then
    if (me.eq.0) then
      print *, "The number of processors must be 1,2,4,8,12,16,...&
                         & times np_npar."
      print *, "Exiting..."
    end if
    stop
  end if

  npsep = nprocs/np_npar
  msep = mod(me,npsep)
  m_npar = me/npsep+1

  role = mod(msep,4)
  same = (npsep-role+3)/4
  sgme = msep/4

  allocate(me_lkp(0:same-1,0:3,np_npar))

  do ipc = 1, np_npar
    do ipb = 0, 3
      do ipa = 0, same-1
        me_lkp(ipa,ipb,ipc) = npsep*(ipc-1) + 4*ipa + ipb
      enddo
    enddo
  enddo

   allocate(emin(0:npsep-1),emax(0:npsep-1))
   allocate(nind_comm(0:npsep-1))
   allocate(sep_comm(1:np_npar))


   if (me.eq.0) then
      print *, "*******************************************"
      print *, "* LEMan: Parallel version (MPI + Open MP) *"
      write (*,'(a2,i4,a12,a26)') "*", nprocs, " MPI task(s)", "*"
      print *, "* Solver using BABE (Burn at Both Ends)   *"
      print *, "*******************************************"
      print *, " "
   end if
   do ip = 0, npsep-1
     emin(ip) = 1+int(ip*(NelR/real(npsep)))
     emax(ip) = int((ip+1)*(NelR/real(npsep)))
     call MPI_COMM_SPLIT(MPI_COMM_WORLD, msep, m_npar-1, nind_comm(ip), ierr)
   enddo

  do ip = 1, np_npar
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, m_npar, msep, sep_comm(ip), ierr)
  enddo

end subroutine

!-------------------------------------------
! Location matrix LMat: 
! element number e + local basis function number a -->  global bf number LMat(a,e)
subroutine Init_LM
 integer a,e
  do a = 1, Nbf
    do e = 1, NelR+1
      LMat(a,e) = 2*(e-1) + a
    enddo
 enddo
end subroutine
!-------------------------------------------
subroutine Initialize_Gauss
  
   xgauss(2,:) = (/ (-0.577350269189626_dp + 1._dp)/2._dp, &
                  (+0.577350269189626_dp + 1._dp)/2._dp,  0._dp,0._dp,0._dp,0._dp  /)
   wgauss(2,:) = (/ 1.0_dp, 1.0_dp,  0._dp,0._dp,0._dp,0._dp /)

   xgauss(3,:) = (/ (-0.774596669241483_dp + 1._dp)/2._dp, &
               ( 0._dp                + 1._dp)/2._dp, &
               (+0.774596669241483_dp + 1._dp)/2._dp,  0._dp,0._dp,0._dp  /)
   wgauss(3,:) = (/ 0.555555555555556_dp, &
               0.888888888888889_dp, &
               0.555555555555556_dp,  0._dp,0._dp,0._dp /)

   xgauss(4,:) = (/ (-0.86113631_dp + 1._dp)/2._dp, &
               (-0.33998104_dp + 1._dp)/2._dp, &
               (+0.33998104_dp + 1._dp)/2._dp, &
               (+0.86113631_dp + 1._dp)/2._dp,  0._dp,0._dp  /)
   wgauss(4,:) = (/ 0.34785485_dp, &
               0.65214515_dp, &
               0.65214515_dp, &
               0.34785485_dp,  0._dp,0._dp /)

   xgauss(5,:) = (/ (-0.906179845938664_dp + 1._dp)/2._dp, &
               (-0.538469310105683_dp + 1._dp)/2._dp, &
               ( 0._dp + 1._dp)/2._dp, &
               (+0.538469310105683_dp + 1._dp)/2._dp, &
               (+0.906179845938664_dp + 1._dp)/2._dp,  0._dp  /)
   wgauss(5,:) = (/ 0.236926885056189_dp, &
               0.478628670499366_dp, &
               0.568888888888889_dp, &
               0.478628670499366_dp, &
               0.236926885056189_dp,  0._dp /)

   xgauss(6,:) = (/ (-0.932469514203152_dp + 1._dp)/2._dp, &
               (-0.661209386466265_dp + 1._dp)/2._dp, &
               (-0.238619186083197_dp + 1._dp)/2._dp, &
               (+0.238619186083197_dp + 1._dp)/2._dp, &
               (+0.661209386466265_dp + 1._dp)/2._dp, &
               (+0.932469514203152_dp + 1._dp)/2._dp  /)
   wgauss(6,:) = (/ 0.171324492379170_dp, &
               0.360761573048139_dp, &
               0.467913934572691_dp, &
               0.467913934572691_dp, &
               0.360761573048139_dp, &
               0.171324492379170_dp /)

  wgauss(:,:) = wgauss(:,:)/2._dp

  iim = cmplx(0._dp,1._dp,dp)
  pi = asin(1._dp)*2._dp

  bx0 = 0._dp
  bx1 = 1._dp
  by0 = 0._dp
  by1 = 2._dp*pi
  bz0 = 0._dp
  bz1 = 2._dp*pi

end subroutine
!-------------------------------------------
subroutine InverseMatrix3(A,B)
  real(dp) A(3,3), B(3,3)!, C(3,3)
  real(dp) detA

  detA = A(1,1) * (A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
         A(1,2) * (A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
         A(1,3) * (A(2,1)*A(3,2)-A(2,2)*A(3,1))

  B(1,1)= (A(2,2)*A(3,3)-A(2,3)*A(3,2)); B(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2)); B(1,3)= (A(1,2)*A(2,3)-A(1,3)*A(2,2))
  B(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1)); B(2,2)= (A(1,1)*A(3,3)-A(1,3)*A(3,1)); B(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
  B(3,1)= (A(2,1)*A(3,2)-A(2,2)*A(3,1)); B(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1)); B(3,3)= (A(1,1)*A(2,2)-A(1,2)*A(2,1))

  B = B/detA
  
end subroutine
!-------------------------------------------
subroutine TransposeMatrix3(A,B)
  real(dp) A(3,3), B(3,3)

  B = A;
  B(1,2) = A(2,1);  B(1,3) = A(3,1);
  B(2,1) = A(1,2);  B(2,3) = A(3,2);
  B(3,1) = A(1,3);  B(3,2) = A(2,3);

end subroutine
!-------------------------------------------
subroutine DetMultMatrix3(A,B,detAB)
  real(dp) A(3,3), B(3,3), AB(3,3), detAB
  integer i,j,k

  AB = 0
  do i = 1, 3
  do j = 1, 3
    do k = 1, 3
      AB(i,j) = AB(i,j) + A(i,k)*B(k,j)
    enddo
  enddo
  enddo

  detAB = AB(1,1) * (AB(2,2)*AB(3,3)-AB(2,3)*AB(3,2)) - &
          AB(1,2) * (AB(2,1)*AB(3,3)-AB(2,3)*AB(3,1)) + &
          AB(1,3) * (AB(2,1)*AB(3,2)-AB(2,2)*AB(3,1))
  
end subroutine
!-------------------------------------------
subroutine MultMatrixVect3(A,B,C)
  real(dp) A(3,3), B(3), C(3)
  integer i,j

  C = 0
  do i = 1, 3
  do j = 1, 3
    C(i) = C(i) + A(i,j)*B(j)
  enddo
  enddo

end subroutine
!-------------------------------------------
! Compute the product of two square matrices
! On entry: A,B two complex matrices,
!           ms the size of the matrices
! On exit: C the complex result
subroutine MXM(A,B,C,ms)

  integer j,k,l,ms
  complex(dp), dimension(ms,ms) :: A,B,C

  C=0
  do j = 1, ms
    do k = 1, ms
      do l = 1, ms
      C(l,j)=C(l,j)+A(l,k)*B(k,j)
      enddo
    enddo
  enddo

end subroutine
!-------------------------------------------
! Compute the product of a matrix and a vector
! On entry: A a square matrix, vec a vector,
!           ms the size of the matrix
! On exit: vecp the complex result
subroutine MXv(A,vec,vecp,ms)

  integer j,k,ms
  complex(dp), dimension(ms) :: vec, vecp
  complex(dp), dimension(ms,ms) :: A

  vecp=0
    do j = 1, ms
      do k = 1, ms
        vecp(j)=vecp(j)+A(j,k)*vec(k)
      enddo
    enddo

end subroutine
!-------------------------------------------
subroutine ReadNamelist

  ! 'param.in'

  if (me.eq.0) then
    read(fnInput, nml=InputNameList)
    posbuff=0
    call MPI_PACK(Ra, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Rw, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(mp_SI, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(me_SI, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(iA1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(ifr1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(iZ1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(iA2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(ifr2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(iZ2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(iAh, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(iZh, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(n0_SI, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(np0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(np1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(np2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(np3, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(np4, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(np5, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nh0_SI, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nh0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nh1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nh2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nh3, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nh4, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nh5, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Ti0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tip0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tip1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tip2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tip3, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tip4, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tip5, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Te0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tep0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tep1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tep2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tep3, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tep4, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Tep5, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Th0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Thp0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Thp1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Thp2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Thp3, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Thp4, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Thp5, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(An0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Anp0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Anp1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Anp2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Anp3, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Anp4, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Anp5, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(UseDensCoeffs, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(UseTempCoeffs, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(UseHotDensCfs, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(UseHotTempCfs, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(ArtBorder, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(B0, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Bc, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Om_SI1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Om_SI2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(OmSteps, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(NbIter, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(ant_type, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(n_shape, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(m_shape, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(n_ant, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(m_ant, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(s_ant_min, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(s_ant_max, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(theta_ant1, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(theta_ant2, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(theta_ant3, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(theta_ant4, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(om_im, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(om_imi, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(PrecalcMatrix, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(MaxResponseFreq, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(FourierCheck, 1, MPI_LOGICAL, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(EmptyRun, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(OutputFlags, 2, MPI_LOGICAL, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(OutputSelect, 8, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(CylMetric, 1, MPI_LOGICAL, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nperiod, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(meq_cyl, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(neq_cyl, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(ampl_cyl, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(s1_cyl, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(s2_cyl, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(RRcyl, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(DielTensor, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(kparApprox, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(k_par_lkp, 7, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(Nb_intv, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(v_range, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(BorderMethod, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(BorderSurf, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(RelaxParam, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(mEqMax, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nEqMax, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(AbsPower, 4, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(KKMem, 1, MPI_DOUBLE_PRECISION, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(ParaMethod, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(nwidth, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_PACK(np_npar, 1, MPI_INTEGER, buff, lbuff, posbuff, &
      MPI_COMM_WORLD, ierr)
    call MPI_BCAST(buff, lbuff, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
  else
    call MPI_BCAST(buff, lbuff, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
    posbuff=0
    call MPI_UNPACK(buff, lbuff, posbuff, Ra, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Rw, 1, MPI_DOUBLE_PRECISION,  &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, mp_SI, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, me_SI, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, iA1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, ifr1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, iZ1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, iA2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, ifr2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, iZ2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, iAh, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, iZh, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, n0_SI, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, np0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, np1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, np2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, np3, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, np4, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, np5, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nh0_SI, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nh0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nh1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nh2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nh3, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nh4, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nh5, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Ti0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tip0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tip1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tip2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tip3, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tip4, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tip5, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Te0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tep0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tep1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tep2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tep3, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tep4, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Tep5, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Th0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Thp0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Thp1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Thp2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Thp3, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Thp4, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Thp5, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, An0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Anp0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Anp1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Anp2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Anp3, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Anp4, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Anp5, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, UseDensCoeffs, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, UseTempCoeffs, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, UseHotDensCfs, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, UseHotTempCfs, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, ArtBorder, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, B0, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Bc, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Om_SI1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Om_SI2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, OmSteps, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, NbIter, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, ant_type, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, n_shape, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, m_shape, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, n_ant, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, m_ant, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, s_ant_min, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, s_ant_max, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, theta_ant1, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, theta_ant2, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, theta_ant3, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, theta_ant4, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, om_im, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, om_imi, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, PrecalcMatrix, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, MaxResponseFreq, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, FourierCheck, 1, MPI_LOGICAL, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, EmptyRun, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, OutputFlags, 2, MPI_LOGICAL, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, OutputSelect, 8, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, CylMetric, 1, MPI_LOGICAL, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nperiod, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, meq_cyl, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, neq_cyl, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, ampl_cyl, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, s1_cyl, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, s2_cyl, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, RRcyl, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, DielTensor, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, kparApprox, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, k_par_lkp, 7, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, Nb_intv, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, v_range, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, BorderMethod, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, BorderSurf, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, RelaxParam, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, mEqMax, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nEqMax, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, AbsPower, 4, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, KKMem, 1, MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, ParaMethod, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, nwidth, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
    call MPI_UNPACK(buff, lbuff, posbuff, np_npar, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, ierr)
  end if

end subroutine
!-------------------------------------------
! Reads the input file containing the radial grid
subroutine ReadGrid
 integer, parameter :: fn = 2
 integer :: eR, gR, i

if (me.eq.0) then
  open(fn, file=Trim(GRIDPATH)//'grid.dat')
  read(fn, *) i
  if (i/=NelR) then
    print *, 'Error: the radial grid doesn`t correspond to actual parameters.'
    print *, 'NelR(LEMan) = ',NelR, ',    NelR(grid.dat) = ',i,'        Exiting...'
    stop
  endif
  do i = 1, NelR+1
    read(fn, *) s_nod(i)
  enddo
  close(fn)
  call MPI_BCAST(s_nod, NelR+1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
else
  call MPI_BCAST(s_nod, NelR+1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end if

  do i = 1, NelR+1
    s_nod(i) = bx0 + (bx1-bx0)*s_nod(i)
  enddo

  do eR = 1, NelR  ! Initializing s in gauss points
  do gR = 1, NGaussR 
    s_gauss(eR,gR) = s_nod(eR) + (s_nod(eR+1)-s_nod(eR))*xgauss(NGaussR,gR)
  enddo
  enddo


  ! Half-mesh points
  do eR = 1, NelR
    s_hp(eR) = (s_nod(eR)+s_nod(eR+1))/2
  enddo

end subroutine
!-------------------------------------------
! Reads the perturbed mode table used for the Fourier decomposition (AFtable.txt)
! and calculates TU,TL,bjac,iota etc
subroutine ReadAFTable

 integer, parameter :: fn = 2
 integer mmin,mm,mmax,nmin,nn,nmax,mn
 integer, dimension(:,:), allocatable :: mnTable

 integer m,n,i,j
 character(253) s, stmp

  if (me.eq.0) then
    open(fn, file=AFTableFILE)
    read(fn, *)
    read(fn, *)
    read(fn, *)
    read(fn, '(1x,4i10)') mmin,mm,nmin,nn
    read(fn, *)
    read(fn, *)

    call MPI_BCAST(mmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    mmax = mmin + mm - 1
    nmax = nmin + (nn-1)*nperiod
    allocate(mnTable(1:mm,1:nn))

    write(stmp,'(i4)') mmax-mmin+1
    s = '(5x,' // trim(stmp) // 'i2,i5)'
    s = trim(s)

    do n = 1,nn
      read(fn, s) (mnTable(m,n),  m=1,mm)
    enddo

    call MPI_BCAST(mnTable, mm*nn, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    close(fn)

  else

    call MPI_BCAST(mmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    mmax = mmin + mm - 1
    nmax = nmin + (nn-1)*nperiod
    allocate(mnTable(1:mm,1:nn))

    call MPI_BCAST(mnTable, mm*nn, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  end if

  if (me.eq.0) then

    print *, ' '
    print *, 'AF mode table'
    s = ' '
    do m = 3, mm, 4
      write(stmp,'(i4)') mmin+m-1
      s = trim(s) // stmp
    enddo
    write(*,'(a)') trim(s)

    do n = 1, nn
    do i = 1, nwidth
      s=' '
      do m = 1, mm
       write(stmp,'(i1)') mnTable(m,n)
       s = trim(s) // stmp
      enddo
       write(stmp,'(i5)') nmin + (n+i-2)*nperiod
       s = ' ' // trim(s) // stmp
      if (sum(mnTable(:,n))/=0) write(*,'(a)') trim(s)
    enddo
    enddo
  end if

  mnAFTot = sum(mnTable)
  if (me.eq.0) then
    print *, 'Total: ',nwidth*mnAFTot, '  modes'
  end if

  allocate(mAF(mnAFTot))
  allocate(nAF(mnAFTot,nwidth/np_npar+1))
  allocate(mpAF(mnAFTot))
  allocate(npAF(mnAFTot,nwidth/np_npar+1))

  nAF = 0

  i = 0
  do m = 1, mm
  do n = 1, nn
    if (mnTable(m,n) == 1) then
      i = i + 1
      mAF(i) = mmin + m - 1
      mpAF(i) = -mAF(i)

      j = 0
      do nindex = m_npar, nwidth, np_npar
        j = j+1
        nAF(i,j) = nmin + (n+nindex-2)*nperiod
        npAF(i,j) = -nAF(i,j)
      enddo

! CAREFUL!!!
!      if ((abs(mAF(i))>njT/2).or.(abs(nAF(i)/real(nperiod))>nkT/2)) then
!        print *, 'Mode number too high for the selected angular grid'
!        print *, 'm =', mAF(i), ',  nj = ',njT,',   n =', nAF(i), ',  nk = ',nkT
!        stop
!      endif

!      if ((abs(nAF(i)/real(nperiod))>nkT/2).and.me.eq.0) then
!        print *, 'Mode number too high for the selected angular grid'
!        print *, 'm =', mAF(i), ',  nj = ',njT,',   n =', nAF(i), &
!           ',  nk =',nkT
!        stop
!      endif
    endif
  enddo
  enddo

!  do i = 1, mnAFTot
!    print '(1x,a,i4,a,i4,a,i4,a,i4)', 'm = ',mAF(i), ',  n = ',nAF(i), ',  mp = ',mpAF(i), ',  np = ',npAF(i)
!  enddo


  deallocate(mnTable)


  if (me.eq.0) then

    open(fn, file=DTTableFILE)
    read(fn, *)
    read(fn, *)
    read(fn, *)
    read(fn, '(1x,4i10)') mmin,mm,nmin,nn
    read(fn, *)
    read(fn, *)

    call MPI_BCAST(mmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    mmax = mmin + mm - 1
    nmax = nmin + (nn-1)*nperiod
    allocate(mnTable(1:mm,1:nn))

    write(stmp,'(i4)') mmax-mmin+1
    s = '(5x,' // trim(stmp) // 'i2,i5)'
    s = trim(s)

    do n = 1,nn
      read(fn, s) (mnTable(m,n),  m=1,mm)
    enddo

    call MPI_BCAST(mnTable, mm*nn, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    close(fn)

  else

    call MPI_BCAST(mmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    mmax = mmin + mm - 1
    nmax = nmin + (nn-1)*nperiod
    allocate(mnTable(1:mm,1:nn))

    call MPI_BCAST(mnTable, mm*nn, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  end if

  if (me.eq.0) then
    print *, ' '
    print *, 'DT mode table'
    s = ' '
    do m = 3, mm, 4
      write(stmp,'(i4)') mmin+m-1
      s = trim(s) // stmp
    enddo
    write(*,'(a)') trim(s)

    do n = 1, nn
    do i = 1, nwidth
      s=' '
      do m = 1, mm
       write(stmp,'(i1)') mnTable(m,n)
       s = trim(s) // stmp
      enddo
       write(stmp,'(i5)') nmin + (n+i-2)*nperiod
       s = ' ' // trim(s) // stmp
      if (sum(mnTable(:,n))/=0) write(*,'(a)') trim(s)
    enddo
    enddo
  end if


  mnDTTot = sum(mnTable)
  if (me.eq.0) then
    print *, 'Total: ',nwidth*mnDTTot, '  modes'
  end if

  allocate(mDT(mnDTTot))
  allocate(nDT(mnDTTot,nwidth/np_npar+1))
  allocate(mpDT(mnDTTot))
  allocate(npDT(mnDTTot,nwidth/np_npar+1))
  allocate(mn_mnDT_lkp(mnDTTot))

  mn_mnDT_lkp = 0

  i = 0
  do m = 1, mm
  do n = 1, nn
    if (mnTable(m,n) == 1) then
      i = i + 1
      mDT(i) = mmin + m - 1
      mpDT(i) = -mDT(i)

      j = 0
      do nindex = m_npar, nwidth, np_npar
        j = j+1
        nDT(i,j) = nmin + (n+nindex-2)*nperiod
        npDT(i,j) = -nDT(i,j)
      enddo

      do mn = 1, mnAFTot
        if(mDT(i)==mAF(mn).and.nDT(i,1)==nAF(mn,1)) then
          mn_mnDT_lkp(i) = mn
        end if
      enddo
      if (mn_mnDT_lkp(i)==0) then
        if (me.eq.0) &
          print *, "AF modes tables are not chosen properly. Exting...", i
        stop
      end if

    endif
  enddo
  enddo

  deallocate(mnTable)

end subroutine
!-------------------------------------------
! Reads the equilibrium on the TERPSICHORE grid (from fort.37), interpolates it to the new angular grid
! and calculates TU,TL,bjac,iota etc
subroutine ReadEquilibriumData
 integer, parameter :: fn = 2
 integer is,itTO,ipTO,itpTO, it,ip,itp, itpr, gR

 integer jkTO
 ! toroidal magn. flux, poloidal m. flux,  &
 ! toroidal current flux, poloidal cur. flux, corresponding radial derivatives
 real(dp) tmp(11)

 integer itp_it_ip_lkpTO(njTO,nkTO)

 real(dp) s, ksi, s_mult, theta,phi, angle_mult
 real(dp) dx,dy,thG
 real(dp) th_ant_rad(2), thG_closest(2)


  ! Initialize lookup table for poloidal & toroidal indices -- TERPSICHORE mesh
  do it = 1, njT
  do ip = 1, nkT
    itp_it_ip_lkp(it,ip) = (ip-1)*njT + it
  enddo
  enddo
  ! Initialize lookup table for poloidal & toroidal indices -- the new angular mesh
  do itTO = 1, njTO
  do ipTO = 1, nkTO
    itp_it_ip_lkpTO(itTO,ipTO) = (ipTO-1)*njTO + itTO
  enddo
  enddo

! TERPSICHORE mesh
! itpr = (gR-1)*Njk + itp
  do itpr = 1, Njkr
    itp_itpr_lkp(itpr) = mod(itpr-1,Njk) + 1
    ip_itpr_lkp(itpr)  = (itp_itpr_lkp(itpr)-1)/njT + 1
    it_itpr_lkp(itpr) = mod(itp_itpr_lkp(itpr)-1,njT) + 1
    gR_itpr_lkp(itpr) = (itpr - itp_itpr_lkp(itpr))/Njk + 1
  enddo
  do gR = 1, NGaussR
    do itp = 1, Njk
      itpr_itp_gR_lkp(itp,gR) = (gR-1)*Njk + itp
    enddo
  enddo


! New angular grid
  allocate(tTU(DIM,DIM, Njk,2));
  allocate(tTL(DIM,DIM, Njk,2));
  allocate(tbjac(Njk,2));
  allocate(tbmod(Njk,2));
  allocate(trx(Njk,2));
  allocate(trz(Njk,2));
  allocate(tphiv(Njk,2));
  allocate(tgssu(Njk,2));


 MajRadius = RRCyl   ! <- cylindrical metric
! nperiod   = 1  ! cylinder
 if (.not.(CylMetric)) then


! TERPSICHORE grid
  allocate(tbjacTO(NjkTO,niT));
  allocate(tbmodTO(NjkTO,niT));
  allocate(tgssuTO(NjkTO,niT));
  allocate(tgpplTO(NjkTO,niT));
  allocate(tgtplTO(NjkTO,niT));
  allocate(tgttlTO(NjkTO,niT));
  allocate(tgstlTO(NjkTO,niT));
  allocate(tgsslTO(NjkTO,niT));
  allocate(tgsplTO(NjkTO,niT));
  allocate(trxTO(NjkTO,niT));
  allocate(trzTO(NjkTO,niT));
  allocate(tb2TO(NjkTO,niT));
  allocate(tphivTO(NjkTO,niT));

  if (me.eq.0) then

  open(fn, file=Trim(EQPATH)//EQFILE, form='unformatted')
  read(fn) is,itTO,ipTO, nperiod
  if ((abs(is)>1000).or.(abs(itTO)>1000).or.(abs(ipTO)>1000)) then
    rewind(fn)
    read(fn) eqname
    read(fn) is,itTO,ipTO, nperiod
  else
    eqname = 'Not specified'
  endif
  if ((is/=niT).or.(itTO/=njTO).or.(ipTO/=nkTO)) then
    print *, 'TERPSICHORE grid: ( ni=',is ,',nj=',itTO,',nk=',ipTO ,')'
    print *, 'LEMan eq. grid:   ( ni=',niT,',nj=',njTO,',nk=',nkTO,')'
    print *, 'Error: the metric doesn`t correspond to actual parameters (ni,nj,nk). Exiting...'
    stop
  endif
  write(*,'(a24,i4,a6,i4,a6,i4,a1)') ' TERPSICHORE grid: (ni =',niT, &
                                     ', nj =',njTO,', nk =',nkTO,')'
  write(*,'(a39,i4,a6,i4,a1)') ' Interpolated LEMan angular grid: (nj =',njT,', nk =',nkT,')'
  if ((njT<njTO).or.(nkT<nkTO)) then 
    print *, 'The new interpolated angular grid is coarser than the original. Exiting...'
    stop
  endif
  print *, ' '

!  1D values f(s)
  read(fn) s_T,ftp,fpp,ci,cj,cip,cjp

! Scale the fluxes (in T. they are one toroidal period fluxes)
!  fpp = fpp*nperiod
!  ci  = ci*nperiod
!  cip = cip*nperiod

! Reading values in the real space (Boozer coordinates).  One toroidal period!!!

  read(fn) tbjacTO, tgssuTO, tgpplTO, tgtplTO, tgttlTO, &
           tgstlTO, tgsslTO, tgsplTO, trxTO, trzTO, tb2TO, tphivTO

  close(fn)

! Read density and ion temperature from file

  if (UseHotDensCfs==3.or.UseHotTempCfs==3) then
  open(fn,file='moments.in')
  read(fn,*) Bc
  allocate(moments(0:niT,3))
  do it=1,niT
      read(fn,*)moments(it,1:3)
  enddo
  moments(0,1:3)=1.5*moments(1,1:3)-0.5*moments(2,1:3)
  do it=1,3
      if(moments(0,it).lt.0.)moments(0,it)=moments(1,it)
  enddo
  end if


  if (UseDensCoeffs==3.or.UseTempCoeffs==3) then
  open(fn,file='ne_Te_fits.in')
  read(fn,*) nif
  allocate (nT_fits(nif,3))
  do it=1,nif
    read(fn,*) nT_fits(it,1:3)
  enddo
  end if

!j Read 2nd population density and ion temperature from file

  if (UseDensCoeffs==4.or.UseTempCoeffs==4) then
  open(fn,file='thermal.in')
  allocate(thermal(0:niT,3))
  do it=1,niT
      read(fn,*)thermal(it,1:3)
  enddo
  thermal(0,1:3)=1.5*thermal(1,1:3)-0.5*thermal(2,1:3)
  do it=1,3
      if(thermal(0,it).lt.0.)thermal(0,it)=thermal(1,it)
  enddo
  end if



  call MPI_BCAST(nperiod, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(s_T, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ftp, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(fpp, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ci, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(cj, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(cip, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(cjp, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tbjacTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgssuTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgpplTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgtplTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgttlTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgstlTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgsslTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgsplTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(trxTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(trzTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tb2TO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tphivTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (UseHotDensCfs==3.or.UseHotTempCfs==3) then
   call MPI_BCAST(moments, 3*(niT+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Bc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  endif
  if (UseDensCoeffs==3.or.UseTempCoeffs==3) then
   call MPI_BCAST(nif, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(nT_fits, 3*nif, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end if
  if (UseDensCoeffs==4.or.UseTempCoeffs==4) then
   call MPI_BCAST(thermal, 3*(niT+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end if

  else

  call MPI_BCAST(nperiod, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(s_T, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ftp, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(fpp, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ci, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(cj, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(cip, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(cjp, niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tbjacTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgssuTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgpplTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgtplTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgttlTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgstlTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgsslTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tgsplTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(trxTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(trzTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tb2TO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tphivTO, NjkTO*niT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (UseHotDensCfs==3.or.UseHotTempCfs==3) then
   allocate(moments(0:niT,3))
   call MPI_BCAST(moments, 3*(niT+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Bc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end if
  if (UseDensCoeffs==3.or.UseTempCoeffs==3) then
   call MPI_BCAST(nif, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   allocate (nT_fits(nif,3))
   call MPI_BCAST(nT_fits, 3*nif, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end if
  if (UseDensCoeffs==4.or.UseTempCoeffs==4) then
   allocate(thermal(0:niT,3))
   call MPI_BCAST(thermal, 3*(niT+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end if

  end if

!  if(UseHotTempCfs.eq.3)Th0=moments(0,2)
!  if(UseHotDensCfs.eq.3)nh0_SI=moments(0,1)

  if (me.eq.0) then
    write(*,'(a9,es11.4,a10,es11.4)') ' s_T(1) =', s_T(1), ', s_T(2) =', s_T(2)
    write(*,'(a17,es11.4)') ' gauss point(1) =', s_gauss(1,1)
    print *, ' '
  end if

  is = 1;  MajRadius = 0
  do ipTO = 1, nkTO
    itpTO = itp_it_ip_lkpTO(1,ipTO); tmp(1) = trxTO(itpTO,is); tmp(2) = tmp(1)
    do itTO = 1, njTO
      itpTO = itp_it_ip_lkpTO(itTO,ipTO)
      if (trxTO(itpTO,is)<tmp(1)) tmp(1) = trxTO(itpTO,is)
      if (trxTO(itpTO,is)>tmp(2)) tmp(2) = trxTO(itpTO,is)  
    enddo
    MajRadius = MajRadius + (tmp(1)+tmp(2))/2
  enddo
  MajRadius = MajRadius / nkTO
  if (me.eq.0) then
    write(*,'(a19,i7)')   ' nperiod          =', nperiod
    write(*,'(a19,f7.2,a2)') ' Major radius     =', MajRadius, ' m'
  end if

  ! Aspect ratio in the ipTO=1 toroidal plane  
  ipTO = 1
  do is = 1, niT
    itpTO = itp_it_ip_lkpTO(1,ipTO); tmp(1) = trxTO(itpTO,is); tmp(2) = tmp(1)
    do itTO = 1, njTO
      itpTO = itp_it_ip_lkpTO(itTO,ipTO)
      if (trxTO(itpTO,is)<tmp(1)) tmp(1) = trxTO(itpTO,is)
      if (trxTO(itpTO,is)>tmp(2)) tmp(2) = trxTO(itpTO,is)  
    enddo
    iaspct(is) = (tmp(2)-tmp(1))/(tmp(2)+tmp(1))
  enddo

  do is = 1, niT
    iota(is) = fpp(is)/ftp(is)  ! Total iota
  enddo
  if (me.eq.0) then
    write (*,'(a19,f7.2,a8,f7.2)') ' Total iota: i(1) =',iota(1), &
                                   ', i(ni)=',iota(niT)
    write (*,'(a19,f7.2,a8,f7.2)') ' Per period: i(1) =',iota(1)/nperiod, &
                                   ', i(ni)=',iota(niT)/nperiod
    write(*,'(a19,f7.2,a2)') ' R*q(total)_axis  =', &
                               MajRadius/iota(1), ' m'
  end if


  call InterpEqSurf(1)
  if (me.eq.0) write(*,'(a19,f7.2,a2)') ' B0               =', B0, ' T'

  do ip = 1, nkT
    itp = itp_it_ip_lkp(1,ip); 
    tmp(1) = trx(itp,1); tmp(2) = tmp(1)
    tmp(3) = trz(itp,1); tmp(4) = tmp(3)
    do it = 1, njT
      itp = itp_it_ip_lkp(it,ip)
      if (trx(itp,1)<tmp(1)) tmp(1) = trx(itp,1)
      if (trx(itp,1)>tmp(2)) tmp(2) = trx(itp,1)  
      if (trz(itp,1)<tmp(3)) tmp(3) = trz(itp,1)
      if (trz(itp,1)>tmp(4)) tmp(4) = trz(itp,1)  
    enddo
    MAxisR(ip) = (tmp(1)+tmp(2))/2
    MAxisZ(ip) = (tmp(3)+tmp(4))/2
  enddo

  else

   print *, '"Cylindrical" metric + 3D terms: (ni=',niT,',nj=',njT,& 
    &  ',nk=',nkT,')'
   print *, 'me=',meq_cyl,',ne=',neq_cyl,',   ampl=',ampl_cyl
   print *, ' '


!   MajRadius = 5._dp
   MAxisR(:) = MajRadius
   MAxisZ(:) = 0

   do is = 1, niT
     s_T(is) = (is-0.5_dp)/niT
   enddo

  end if ! CylMetric

!  call InterpolateTERPSMetric

  do it = 0, njT+1  ! Initializing theta vector
    theta_T(it) = 2*pi*(it-1)/njT
  enddo
  do ip = 0, nkT+1  ! Initializing phi vector
    phi_T(ip)   = 2*pi*(ip-1)/nkT/nperiod
  enddo  
  do it = 1,njT
  do ip = 1,nkT
    itp = itp_it_ip_lkp(it,ip)
    th_itp(itp) = theta_T(it)
    ph_itp(itp) = phi_T(ip)
  enddo
  enddo

  call InterpEqSurf(niT-1)


if (ant_type == 1 .and. m_shape == 0) then

  theta_ant1 = theta_ant1*pi/(180*nperiod)
  theta_ant2 = theta_ant2*pi/(180*nperiod)

elseif (ant_type /= 0) then

  th_ant_rad(1) = theta_ant1*pi/180
  th_ant_rad(2) = theta_ant2*pi/180
  thG_closest = 0
  ip = 1
  do it = 1, njT
    itp = itp_it_ip_lkp(it,ip)

    ! Get the geometrical angle 0<=th<=2*pi
    dx = (trx(itp,2)-MAxisR(ip))
    dy = (trz(itp,2)-MAxisZ(ip))
    thG = atan2(dy,dx)
    if (thG<0) thG = thG + 2*pi

    if (abs(thG - th_ant_rad(1))<abs(thG_closest(1) - th_ant_rad(1))) then
      thG_closest(1) = thG
      theta_ant1  = theta_T(it)
    endif

    if (abs(thG - th_ant_rad(2))<abs(thG_closest(2) - th_ant_rad(2))) then
      thG_closest(2) = thG
      theta_ant2  = theta_T(it)
    endif

  enddo

  theta_ant3 = theta_ant3*pi/(180*nperiod)
  theta_ant4 = theta_ant4*pi/(180*nperiod)

end if


if (ant_type == 1) then

  if (m_shape*m_ant+n_shape*n_ant/nperiod /= 0) then
    if (me==0) then
      print *, 'Antenna mode and shape are not compatible. Exiting...'
    endif
    stop
  endif

  if (me==0) then
    if (m_shape==0) then

      print *, 'Toroidally localised antenna.'
    write(*,'(a21,f5.1,a11,f5.1,a8)') ' Extension:  phi1  = ', &
      theta_ant1*(180*nperiod)/pi, ',  phi2  = ', theta_ant2*(180*nperiod)/pi, ' degrees'

    elseif (n_shape==0) then

      print *, 'Poloidally localised antenna.'
      print *, 'Extension: theta1(Boozer)  = ', theta_ant1*180/pi, &
        &  ',   theta2(Boozer)  = ', theta_ant2*180/pi

    else

      print *, 'Helically localised antenna with shape  m  = ', &
        &  m_shape, ',   n  = ', n_shape
      print *, 'Extension at phi = 0 :  theta1(Boozer)  = ', &
        &  theta_ant1*180/pi, ',   theta2(Boozer)  = ', theta_ant2*180/pi

    end if
  end if

endif

if (ant_type == 2) then

  if (me==0) then
    print *, 'Localised antenna.'
    print *, 'Poloidal extension: theta1(Boozer)  = ', theta_ant1*180/pi, &
        &  ',   theta2(Boozer)  = ', theta_ant2*180/pi, 'degrees'
  end if

end if


end subroutine

!--------------------------------------------------------------------------------
! Interpolate data in poloidal and toroidal direction for one magnetic surface
subroutine InterpEqSurf(ism)

  implicit none

  integer is, ism, it, ip, itp
  real(dp) M3(3,3), iM3(3,3), detM
  real(dp) TL3(3,3), iTU3(3,3)
  real(dp) s, ksi, s_mult, theta,phi, angle_mult


  if (.not.(CylMetric)) then

  allocate(tgppl(Njk,2));
  allocate(tgtpl(Njk,2));
  allocate(tgttl(Njk,2));
  allocate(tgstl(Njk,2));
  allocate(tgssl(Njk,2));
  allocate(tgspl(Njk,2));
  allocate(tb2(Njk,2));


  ! Interpolate the equilibrium to the new angular grid
  call InterpAnglesEq(ism,njTO,nkTO,tbjacTO,niT,njT,nkT,tbjac,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tgssuTO,niT,njT,nkT,tgssu,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tgpplTO,niT,njT,nkT,tgppl,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tgtplTO,niT,njT,nkT,tgtpl,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tgttlTO,niT,njT,nkT,tgttl,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tgstlTO,niT,njT,nkT,tgstl,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tgsslTO,niT,njT,nkT,tgssl,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tgsplTO,niT,njT,nkT,tgspl,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,trxTO,niT,njT,nkT,trx,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,trzTO,niT,njT,nkT,trz,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tb2TO,niT,njT,nkT,tb2,forw,backw)
  call InterpAnglesEq(ism,njTO,nkTO,tphivTO,niT,njT,nkT,tphiv,forw,backw)


!  print *
!  print *
!  print *
!  print *, 'COMPARE:'
!  print *, tbjacTO(1,50), tbjac(1,50)
!  print *, tbjacTO(1,30), tbjac(1,30)
!  print *
!  print *
!  print *
!  print *

  ! Calculation of the transformation matrices, new angular grid, real space, 3D
  ! Attention! Lower metric coefficients are divided by bjac
  ! P.ex.  tgspl = gspl_real / bjac

  tTU = 0
  do is = 1, 2
   do it = 1, njT
    do ip = 1, nkT
      itp = itp_it_ip_lkp(it,ip)
      tTU(1,1,itp,is) = sqrt(tgssu(itp,is))
      tTU(2,1,itp,is) = (tgspl(itp,is)*tgtpl(itp,is) - tgstl(itp,is)*tgppl(itp,is)) &
                          /sqrt(tgssu(itp,is))
      tTU(2,2,itp,is) = -ci(ism+is-1)/(sqrt(tb2(itp,is)*tgssu(itp,is))*tbjac(itp,is))
      tTU(2,3,itp,is) = fpp(ism+is-1)/(sqrt(tb2(itp,is))*tbjac(itp,is))
      tTU(3,1,itp,is) = (tgstl(itp,is)*tgtpl(itp,is) - tgspl(itp,is)*tgttl(itp,is)) &
                          /sqrt(tgssu(itp,is))
      tTU(3,2,itp,is) = -cj(ism+is-1)/(sqrt(tb2(itp,is)*tgssu(itp,is))*tbjac(itp,is))
      tTU(3,3,itp,is) = ftp(ism+is-1)/(sqrt(tb2(itp,is))*tbjac(itp,is))


!      if ((is == 50).and.(it == 3).and.(ip == 30)) then
!        print *, tTU(1,1,itp,is)
!        print *, tTU(2,1,itp,is), tTU(2,2,itp,is), tTU(2,3,itp,is)
!        print *, tTU(3,1,itp,is), tTU(3,2,itp,is), tTU(3,3,itp,is)
!      endif

      M3(:,:) = tTU(:,:,itp,is)
      call InverseMatrix3(M3,iM3)
      iTU3 = iM3
      call TransposeMatrix3(iM3,M3)
  
      tTL(:,:,itp,is) = M3(:,:)
      TL3 = tTL(:,:,itp,is)
  
      call DetMultMatrix3(TL3,iTU3,detM)

      tbjac(itp,is)    = - sqrt(detM)

    enddo
   enddo
  enddo


!  tTL = 0
!  do is = 1, niT
!    do it = 1, njT
!    do ip = 1, nkT
!      itp = itp_it_ip_lkp(it,ip)
!      tTL(1,1,itp,is) = 1/sqrt(tgssu(itp,is))
!      tTL(1,2,itp,is) = -(tgstl(itp,is)*ci(is)+tgspl(itp,is)*cj(is))  &
!                       /sqrt(tb2(itp,is)*tgssu(itp,is))
!      tTL(1,3,itp,is) =  (tgstl(itp,is)*fpp(is)+tgspl(itp,is)*ftp(is))  &
!                       /sqrt(tb2(itp,is))
!      tTL(2,2,itp,is) = ftp(is)*sqrt(tgssu(itp,is)/tb2(itp,is))
!      tTL(2,3,itp,is) = cj(is)/sqrt(tb2(itp,is))
!
!      tTL(3,2,itp,is) = -fpp(is)*sqrt(tgssu(itp,is)/tb2(itp,is))
!      tTL(3,3,itp,is) = -ci(is)/sqrt(tb2(itp,is))
!    enddo
!    enddo
!  enddo

  if (ism==1) then
    B02 = 0._dp;
    do ip = 1, nkT
    do it = 1, njT
      itp = itp_it_ip_lkp(it,ip)
      B02 = B02 + sqrt(tb2(itp,1))
    enddo
    enddo
    B02 = B02 / njT / nkT
  end if

  if(B0 < 0.)B0=B02     ! take TERP B0 if input B0<0

  do is = 1,2
    do it = 1,njT
      do ip = 1,nkT
        itp = itp_it_ip_lkp(it,ip)
        tbmod(itp,is) = sqrt(tb2(itp,is))/B02
      enddo
    enddo
  enddo

  deallocate(tgppl)
  deallocate(tgtpl)
  deallocate(tgttl)
  deallocate(tgstl)
  deallocate(tgssl)
  deallocate(tgspl)
  deallocate(tb2)

 else 
 ! CylMetric = True

   tTU = 0

   do is = 1, 2

     s = s_T(ism+is-1)
   
     ksi = 2._dp*(s-s1_cyl)/(s2_cyl-s1_cyl)-1._dp
     s_mult = 0._dp
     if ( (s>s1_cyl).and.(s<s2_cyl) ) &
       s_mult = (1._dp - ksi*ksi)*(1._dp - ksi*ksi)

     do it = 1, njT
       do ip = 1, nkT
         itp = itp_it_ip_lkp(it,ip)
   
         theta   = 2*pi*(it-1)/njT
         phi     = 2*pi*(ip-1)/nkT !/nperiod
         angle_mult = sin(meq_cyl*theta + neq_cyl*phi)


!       tTU(1,1,itp,is) =  0.6_dp              * (1._dp + ampl_cyl*s_mult*angle_mult)
!                                                         
!       tTU(2,1,itp,is) =  0.12e-01_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(2,2,itp,is) = -0.6_dp              * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(2,3,itp,is) =  0.87e-01_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!                                                  
!       tTU(3,1,itp,is) = -0.11e-04_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(3,2,itp,is) =  0.58e-03_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(3,3,itp,is) =  0.47e-02_dp             * (1._dp + ampl_cyl*s_mult*angle_mult)


!       tTU(1,1,itp,is) =  0.6_dp  !            * (1._dp + ampl_cyl*s_mult*angle_mult)
!                                                         
!!       tTU(2,1,itp,is) =  0.12e-01_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(2,2,itp,is) = -0.6_dp  !            * (1._dp + ampl_cyl*s_mult*angle_mult)
!!       tTU(2,3,itp,is) =  0.87e-01_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!                                                  
!!       tTU(3,1,itp,is) = -0.11e-04_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!!       tTU(3,2,itp,is) =  0.58e-03_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(3,3,itp,is) =  0.47e-02_dp   !          * (1._dp + ampl_cyl*s_mult*angle_mult)



!       tTU(1,1,itp,is) =  1.1_dp
!       tTU(2,2,itp,is) =  -1.0_dp
!       tTU(3,3,itp,is) =  0.7_dp


!       tTU(1,1,itp,is) =  2.0_dp*sqrt(s)     * (1._dp + ampl_cyl*s_mult*angle_mult)
!
!       tTU(2,1,itp,is) =  0.1_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(2,2,itp,is) =  1.0_dp/sqrt(s)     * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(2,3,itp,is) =  0*0.2_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!
!       tTU(3,1,itp,is) = (-0.01_dp) * s_mult   * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(3,2,itp,is) = 0.02_dp * s_mult   * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(3,3,itp,is) = -1.0_dp             * (1._dp + ampl_cyl*s_mult*angle_mult)



!         tTU(1,1,itp,is) =  2.0_dp*sqrt(s)
!         tTU(2,2,itp,is) =  1.0_dp/sqrt(s)
!         tTU(3,3,itp,is) = -1.0_dp

         tTU(1,1,itp,is) =  2.0_dp*sqrt(s)/Rw
         tTU(2,2,itp,is) =  1.0_dp/(Rw*sqrt(s))
         tTU(3,3,itp,is) = -1.0_dp/RRcyl

!       tTU(1,1,itp,is) =  1.1_dp              * (1._dp + ampl_cyl*s_mult*angle_mult)
!                                                         
!       tTU(2,1,itp,is) =  0.10_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(2,2,itp,is) =  0.8_dp              * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(2,3,itp,is) =  0.13_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!                                                  
!       tTU(3,1,itp,is) =  0.07_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(3,2,itp,is) =  0.11_dp * s_mult    * (1._dp + ampl_cyl*s_mult*angle_mult)
!       tTU(3,3,itp,is) = -1.13_dp             * (1._dp + ampl_cyl*s_mult*angle_mult)

       M3(:,:) = tTU(:,:,itp,is)
       call InverseMatrix3(M3,iM3)
       iTU3 = iM3
       call TransposeMatrix3(iM3,M3)
   
       tTL(:,:,itp,is) = M3(:,:)
       TL3 = tTL(:,:,itp,is)
   
       call DetMultMatrix3(TL3,iTU3,detM)
   
       tbjac(itp,is)    = - sqrt(detM)

!       if ((is == 100).and.(itp == 1)) then
!
!         print *, 's=',s
!         print *, 'TL:'
!         print *, tTL(1,1,itp,is),  tTL(1,2,itp,is),  tTL(1,3,itp,is)
!         print *, tTL(2,1,itp,is),  tTL(2,2,itp,is),  tTL(2,3,itp,is)
!         print *, tTL(3,1,itp,is),  tTL(3,2,itp,is),  tTL(3,3,itp,is)
!
!         print *, 'tbjac:'
!         print *, tbjac(itp,is)
!
!       endif

       tbmod(itp,is)    = 1._dp
       trx(itp,is)      = MajRadius + sqrt(s) * cos(theta)
       trz(itp,is)      = sqrt(s) * sin(theta)
       tphiv(itp,is)    = 2*pi*(ip-1)/nkT/nperiod

    enddo
   enddo
   enddo


 endif ! CylMetric

end subroutine

!-------------------------------------------
!subroutine InterpolateTERPSMetric
! integer ii,jj, itp
!
! allocate(tTUy2(DIM,DIM, Njk,niT));
! allocate(tTLy2(DIM,DIM, Njk,niT));
! allocate(tbjacy2(Njk,niT));
! allocate(tbmody2(Njk,niT));
!
! do itp = 1, Njk
!
!   do ii = 1, 3
!   do jj = 1, 3
!     vect_tmp(:) = tTU(ii,jj,itp,:)
!     call Spline(s_T,vect_tmp,niT,vect_tmpy2)
!     tTUy2(ii,jj,itp,:) = vect_tmpy2
!
!     vect_tmp(:) = tTL(ii,jj,itp,:)
!     call Spline(s_T,vect_tmp,niT,vect_tmpy2)
!     tTLy2(ii,jj,itp,:) = vect_tmpy2
!   enddo
!   enddo
!
!   vect_tmp(:) = tbjac(itp,:)
!   call Spline(s_T,vect_tmp,niT,vect_tmpy2)
!   tbjacy2(itp,:) = vect_tmpy2
!   vect_tmp(:) = tbmod(itp,:)
!   call Spline(s_T,vect_tmp,niT,vect_tmpy2)
!   tbmody2(itp,:) = vect_tmpy2
!
! enddo
!
!end subroutine
!-------------------------------------------
subroutine RenormalizeEquilibrium
 integer eR, gR
 real(dp) s
 integer is

  do eR = 1, NelR
  do gR = 1, NGaussR 
    s = s_gauss(eR,gR)

    seTU (1,eR,gR) =  Sqrt(s)
    seTU (2,eR,gR) =  1._dp/Sqrt(s)
    dseTU(1,eR,gR) =  0.5_dp/Sqrt(s)
    dseTU(2,eR,gR) = -0.5_dp/(Sqrt(s)*s)

    seTL (1,eR,gR) =  1._dp/Sqrt(s)
    seTL (2,eR,gR) =  Sqrt(s)
    dseTL(1,eR,gR) = -0.5_dp/(Sqrt(s)*s)
    dseTL(2,eR,gR) =  0.5_dp/Sqrt(s)

!    seTU (1,eR,gR) = 1
!    seTU (2,eR,gR) = 1
!    dseTU(1,eR,gR) = 0
!    dseTU(2,eR,gR) = 0
!                     
!    seTL (1,eR,gR) = 1
!    seTL (2,eR,gR) = 1
!    dseTL(1,eR,gR) = 0
!    dseTL(2,eR,gR) = 0
  enddo
  enddo

  do is = 1, niT
    s = s_T(is)

    setTU (1,is) =  Sqrt(s)
    setTU (2,is) =  1._dp/Sqrt(s)
    dsetTU(1,is) =  0.5_dp/Sqrt(s)
    dsetTU(2,is) = -0.5_dp/(Sqrt(s)*s)

    setTL (1,is) =  1._dp/Sqrt(s)
    setTL (2,is) =  Sqrt(s)
    dsetTL(1,is) = -0.5_dp/(Sqrt(s)*s)
    dsetTL(2,is) =  0.5_dp/Sqrt(s)

!    setTU (1,is) = 1
!    setTU (2,is) = 1
!    dsetTU(1,is) = 0
!    dsetTU(2,is) = 0
!
!    setTL (1,is) = 1
!    setTL (2,is) = 1
!    dsetTL(1,is) = 0
!    dsetTL(2,is) = 0
  enddo

end subroutine
!-------------------------------------------
subroutine InitializeVariables
 integer eR, gR, eR_FE, bf, a, b, e, ab
 real(dp) s, tsize
 integer m,n,mp,np,mn,mnp,mnmnp, it,ip, ieq, k,kp
 integer itp, itpr

 integer Dm,Dn, i,j,ij,l
 logical NewEntry
 integer AModified, FModified

 real(dp) pf

 double precision mem

 real(dp) DensDpl(Gwidth+1), TempDpl(Gwidth+1)
 integer Gx
 character(50) TitleChar

 integer, dimension(:), allocatable   :: mEq_tmp, nEq_tmp


  ! matrix.f

  MSize  = Nunknowns*NDIM
  MSize2 = MSize*mnAFTot
  MDiag  = NDIM*mnAFTot*6
  MBlock = NDIM*mnAFTot
  KL     = MBlock*4 - 1
  KU     = KL
  LDAB   = 2*KL + KU + 1
  LBKK   = 2*NDIM*mnAFTot

  ! Treatment of the memory

  tsize = (KKMem*1024**2)/(16*real(LBKK)**2)-4.0-2.0/real(same)

  if (tsize<0) then
    if (me.eq.0) print *, 'Space allocated for main matrix is too small.&
                        & Exiting...'
    stop
  end if
  indexfile = min(NelR/npsep+1,int(tsize))

  if (me.eq.0) print *, ' '
  if(NelR/npsep+1>int(tsize).and.me==0) then
    print *, '* Code using spare memory *'
    write(*,'(f8.3,a22,i4,a35)') KKMem/1.0e3, ' GB allocated to store', &
      indexfile, ' matrix blocks and temporary tables'
    write(*,'(f8.3,a57)') &
      (NelR/npsep+1-indexfile)*real(16*LBKK**2)/real(1024**3), &
      ' GB for the rest of the matrix to be stored on local disk'
  print *, ' '
  end if

  LBKK_red = int((1900.0*1024**2)/(16*real(LBKK)))

  NbDistr = sum(k_par_lkp)+2


  ! Power deposition output

  j = 1
  do i = 1, 4
    if (AbsPower(i)<0.or.AbsPower(i)>1) AbsPower(i)=1
    if (AbsPower(i)==1) then
      AP_lkp(j) = i
      j=j+1
    end if
  enddo

  AP_nb = sum(AbsPower)+1


  ! Table allocation

  allocate(FF(MSize2))
  allocate(XX(Nunknowns,mnAFTot, NDIM))
  allocate(XXc(Nunknowns, nwidth*mnAFTot, NDIM))

  allocate(Exp_mnp(Njk,mnAFTot,nwidth/np_npar+1))
  allocate(Exp_mn(Njk,mnAFTot,nwidth/np_npar+1))
  allocate(DMult_imn(mnAFTot*mnAFTot,(DIM+1)*(DIM+1),nwidth/np_npar+1))
  allocate(GlobalIndex_lkp(NelR,Nbf,mnAFTot,NDIM))
  allocate(Psi_FE(0:1,Nbf, NelR, NelR,NGaussR, 0:1))

  allocate(t11m(NelR/npsep+2,Njkr,AP_nb))
  allocate(t12m(NelR/npsep+2,Njkr,AP_nb))
  allocate(t21m(NelR/npsep+2,Njkr,AP_nb))
  allocate(t22m(NelR/npsep+2,Njkr,AP_nb))
  allocate(t33m(NelR/npsep+2,Njkr,AP_nb))

  allocate(iAF_lkp(mnAFTot*mnAFTot,NDIM,NDIM))

  allocate(vA_fe(NDIM, nwidth*mnAFTot, Nmesh))
  allocate(vAp_fe(NDIM, nwidth*mnAFTot, Nmesh))

  allocate(k_par_it(NbDistr,NelR,Njkr))

  allocate(mneq0_lkp(mnAFTot*mnAFTot))
  allocate(mneq0DT_lkp(mnDTTot*mnDTTot))

  allocate(mnFFT_RHS_lkp(mnDTTot))

  allocate(iF_lkp(mnAFTot,NDIM))
  allocate(iA_lkp(mnAFTot,NDIM))

  allocate(mn_mnmnp_lkp(mnAFTot*mnAFTot))
  allocate(mnp_mnmnp_lkp(mnAFTot*mnAFTot))
  allocate(mn_mnmnpDT_lkp(mnDTTot*mnDTTot))
  allocate(mnp_mnmnpDT_lkp(mnDTTot*mnDTTot))

  allocate(om_a(niT,mnAFTot))

  Rw   = Ra ! fixed-boundary version
  a_plasma = 1.0_dp

!  ----------- Constants ----------- !
!!!Testing
  c_SI     = 2.9979e8_dp       ! m/s
  mu0_SI   = 4*pi*1.0e-7_dp    ! H/m
  qe_Gauss = 4.8032e-10_dp     ! statcoulomb (Gauss)


  mi1_SI    = mp_SI*iA1
  mi2_SI    = mp_SI*iA2
  mih_SI    = mp_SI*iAh

!  ---------- Parameters ----------- !
  OmCi1_SI  = iZ1*qe_Gauss*B0*1.e4_dp/(mi1_SI*1.e3_dp*c_SI*100)     ! rad/s
  OmCi2_SI  = iZ2*qe_Gauss*B0*1.e4_dp/(mi2_SI*1.e3_dp*c_SI*100)     ! rad/s

  ne_SI    = n0_SI*(iZ1*ifr1+iZ2*ifr2)
  cA_SI    = B0*1.e4_dp/sqrt(4*pi*n0_SI* &
             (mi1_SI*ifr1+mi2_SI*ifr2+me_SI)/1.e6_dp*1.e3_dp)/100  ! m/s
  CCa      = c_SI/cA_SI


  ! Initialize density n(s) and temperature T(s) profiles
  do eR = 1,NelR
  do gR = 1,NGaussR
    s = s_gauss(eR,gR)
!    fnp(eR,gR) = 1.0000 + s*(-1.4574) + s**2*0.5430 + s**3*(-0.0943) + s**4*0.0086 + s**5*0.0001
    ! ^^^^--- pressure profile from LION (p.ex. TEST_GAP equilibrium (n = 1-s_LION^2 ) )

    fnp(eR,gR)  = PlasmaDensity(s)
    fnm(eR,gR)  = MinorityDensity(s)
    fnhp(eR,gR) = FastPartDensity(s)
    nu_f(eR,gR) = AbsorpProf(s)
    fti(eR,gR)  = IonTemperature(s)
    fte(eR,gR)  = ElectronTemperature(s)
    ftm(eR,gR)  = MinorityTemperature(s)
    fth(eR,gR)  = FastTemperature(s)
    anis(eR,gR) = FastAnisotropy(s)
  enddo
  enddo
  do eR = 1, niT
    s = s_T(eR)
    dens_T(eR) = max(PlasmaDensity(s),1.0e-15_dp)
  enddo

  
  do i = 0, 1  ! da
  do j = 0, 1  ! Modified/or not
  do eR_FE = 1, NelR
    do eR = 1, NelR
      do gR = 1, NGaussR
        do bf = 1, Nbf
!               da  a   e   s_gauss(eR,gR) Modified/n
          Psi_FE(i,bf,eR_FE, eR,gR,j) = Psi_al(i,bf,eR_FE, s_gauss(eR,gR),j)
        enddo
      enddo
    enddo
  enddo
  enddo
  enddo

  ! Lookup tables for matrix construction

  do e = 1, NelR
  do a = 1, Nbf
    do mn = 1, mnAFTot
      do k = 1, NDIM
        i = LMat(a,e);
        GlobalIndex_lkp(e,a,mn,k) = mod(a+1,2)*mnAFTot*NDIM + (k-1)*mnAFTot + mn
      enddo
    enddo
  enddo
  enddo

if (role<(min(npsep,4)+1)/2) then
  do a = 1, Nbf
    do b = 1, Nbf
      if(b<=2.and.a<=2) then
        ab_lkp(b,a)=1
      elseif(b<=2.and.a>2) then
        ab_lkp(b,a)=2
      elseif(b>2.and.a<=2) then
        ab_lkp(b,a)=3
      else
        ab_lkp(b,a)=4
      end if
    enddo
  enddo
else
  do a = 1, Nbf
    do b = 1, Nbf
      if(b<=2.and.a<=2) then
        ab_lkp(b,a)=4
      elseif(b<=2.and.a>2) then
        ab_lkp(b,a)=3
      elseif(b>2.and.a<=2) then
        ab_lkp(b,a)=2
      else
        ab_lkp(b,a)=1
      end if
    enddo
  enddo
endif

  do ab = 1, Nbf*Nbf
    a_ab_lkp(ab) = (ab-1)/Nbf+1
    b_ab_lkp(ab) = mod(ab-1,Nbf)+1
  enddo

  do mnmnp = 1, mnAFTot*mnAFTot
    mn_mnmnp_lkp(mnmnp)  = (mnmnp-1)/mnAFTot + 1
    mnp_mnmnp_lkp(mnmnp) = mod(mnmnp-1,mnAFTot) + 1
  enddo

  do mnmnp = 1, mnDTTot*mnDTTot
    mn_mnmnpDT_lkp(mnmnp)  = (mnmnp-1)/mnDTTot + 1
    mnp_mnmnpDT_lkp(mnmnp) = mod(mnmnp-1,mnDTTot) + 1
  enddo

  CFCAINFO(:,1) = (/ 1, 0, 2, 6 /)
  CFCAINFO(:,2) = (/ 1, 0, 3, 0 /)
  CFCAINFO(:,3) = (/ 1, 2, 2, 0 /)
  CFCAINFO(:,4) = (/ 1, 3, 2, 0 /)
  CFCAINFO(:,5) = (/ 1, 2, 3, 0 /)
  CFCAINFO(:,6) = (/ 1, 3, 3, 0 /)
  CFCAINFO(:,7) = (/ 2, 0, 1, 8 /)
  CFCAINFO(:,8) = (/ 2, 0, 2, 0 /)
  CFCAINFO(:,9) = (/ 2, 0, 3, 0 /)
  CFCAINFO(:,10) = (/ 2, 2, 1, 0 /)
  CFCAINFO(:,11) = (/ 2, 3, 1, 0 /)
  CFCAINFO(:,12) = (/ 2, 1, 3, 0 /)
  CFCAINFO(:,13) = (/ 2, 2, 3, 0 /)
  CFCAINFO(:,14) = (/ 2, 3, 3, 0 /)
  CFCAINFO(:,15) = (/ 3, 0, 1, 8 /)
  CFCAINFO(:,16) = (/ 3, 0, 2, 0 /)
  CFCAINFO(:,17) = (/ 3, 0, 3, 0 /)
  CFCAINFO(:,18) = (/ 3, 2, 1, 0 /)
  CFCAINFO(:,19) = (/ 3, 3, 1, 0 /)
  CFCAINFO(:,20) = (/ 3, 1, 2, 0 /)
  CFCAINFO(:,21) = (/ 3, 2, 2, 0 /)
  CFCAINFO(:,22) = (/ 3, 3, 2, 0 /)
  CFCAINFO(:,23) = (/ 4, 0, 1, 10 /)
  CFCAINFO(:,24) = (/ 4, 0, 2, 0 /)
  CFCAINFO(:,25) = (/ 4, 0, 3, 0 /)
  CFCAINFO(:,26) = (/ 4, 1, 1, 0 /)
  CFCAINFO(:,27) = (/ 4, 2, 1, 0 /)
  CFCAINFO(:,28) = (/ 4, 3, 1, 0 /)
  CFCAINFO(:,29) = (/ 4, 2, 2, 0 /)
  CFCAINFO(:,30) = (/ 4, 3, 2, 0 /)
  CFCAINFO(:,31) = (/ 4, 2, 3, 0 /)
  CFCAINFO(:,32) = (/ 4, 3, 3, 0 /)
  CFCAINFO(:,33) = (/ 7, 0, 1, 2 /)
  CFCAINFO(:,34) = (/ 7, 0, 2, 0 /)
  CFCAINFO(:,35) = (/ 8, 0, 1, 2 /)
  CFCAINFO(:,36) = (/ 8, 0, 2, 0 /)
  CFCAINFO(:,37) = (/ 9, 0, 3, 1 /)
  CFCAINFO(:,38) = (/ 10, 0, 1, 5 /)
  CFCAINFO(:,39) = (/ 10, 0, 4, 0 /)
  CFCAINFO(:,40) = (/ 10, 1, 4, 0 /)
  CFCAINFO(:,41) = (/ 10, 2, 4, 0 /)
  CFCAINFO(:,42) = (/ 10, 3, 4, 0 /)
  CFCAINFO(:,43) = (/ 11, 0, 2, 5 /)
  CFCAINFO(:,44) = (/ 11, 0, 4, 0 /)
  CFCAINFO(:,45) = (/ 11, 1, 4, 0 /)
  CFCAINFO(:,46) = (/ 11, 2, 4, 0 /)
  CFCAINFO(:,47) = (/ 11, 3, 4, 0 /)
  CFCAINFO(:,48) = (/ 12, 0, 3, 3 /)
  CFCAINFO(:,49) = (/ 12, 2, 4, 0 /)
  CFCAINFO(:,50) = (/ 12, 3, 4, 0 /)
  CFCAINFO(:,51) = (/ 13, 0, 4, 4 /)
  CFCAINFO(:,52) = (/ 13, 1, 4, 0 /)
  CFCAINFO(:,53) = (/ 13, 2, 4, 0 /)
  CFCAINFO(:,54) = (/ 13, 3, 4, 0 /)
  CFCAINFO(:,55) = (/ 14, 0, 4, 4 /)
  CFCAINFO(:,56) = (/ 14, 1, 4, 0 /)
  CFCAINFO(:,57) = (/ 14, 2, 4, 0 /)
  CFCAINFO(:,58) = (/ 14, 3, 4, 0 /)
  CFCAINFO(:,59) = (/ 15, 2, 4, 2 /)
  CFCAINFO(:,60) = (/ 15, 3, 4, 0 /)
  CFCAINFO(:,61) = (/ 16, 0, 1, 6 /)
  CFCAINFO(:,62) = (/ 16, 0, 2, 0 /)
  CFCAINFO(:,63) = (/ 16, 0, 4, 0 /)
  CFCAINFO(:,64) = (/ 16, 1, 4, 0 /)
  CFCAINFO(:,65) = (/ 16, 2, 4, 0 /)
  CFCAINFO(:,66) = (/ 16, 3, 4, 0 /)
  CFCAINFO(:,67) = (/ 17, 0, 1, 4 /)
  CFCAINFO(:,68) = (/ 17, 0, 2, 0 /)
  CFCAINFO(:,69) = (/ 17, 2, 4, 0 /)
  CFCAINFO(:,70) = (/ 17, 3, 4, 0 /)
  CFCAINFO(:,71) = (/ 18, 0, 3, 3 /)
  CFCAINFO(:,72) = (/ 18, 2, 4, 0 /)
  CFCAINFO(:,73) = (/ 18, 3, 4, 0 /)


  ! Allocating memory for storing equlibrium (gauss grid)
  ! Values are assigned in the matrix construction loop only
  allocate(mTL(DIM,DIM,  NGaussR,Njk))
  allocate(mTU(DIM,DIM,  NGaussR,Njk))
  allocate(mTLs(DIM,DIM, NGaussR,Njk))
  allocate(mTUs(DIM,DIM, NGaussR,Njk))
  allocate(mTLt(DIM,DIM, NGaussR,Njk))
  allocate(mTUt(DIM,DIM, NGaussR,Njk))
  allocate(mTLp(DIM,DIM, NGaussR,Njk))
  allocate(mTUp(DIM,DIM, NGaussR,Njk))
  allocate(mgssu(NGaussR,Njk))
  allocate(bjac (NGaussR,Njk))
  allocate(bjacs(NGaussR,Njk))
  allocate(bjact(NGaussR,Njk))
  allocate(bjacp(NGaussR,Njk))
  allocate(bmod (NGaussR,Njk))
  allocate(cR   (NGaussR,Njk))
  allocate(cZ   (NGaussR,Njk))
  allocate(phiv (NGaussR,Njk))
  allocate(cR_hp(NelR,Njk))
  allocate(cZ_hp(NelR,Njk))
  allocate(phiv_hp(NelR,Njk))
  if (OmSteps > 0) then
    allocate(Response(0:OmSteps))
  else
    allocate(Response(0:1))
  endif


 VDATo16(0,1) = 1;  VDATo16(0,2) = 2;  VDATo16(0,3) = 3;  VDATo16(0,4) = 4
 VDATo16(1,1) = 5;  VDATo16(2,1) = 6;  VDATo16(3,1) = 7  
 VDATo16(1,2) = 8;  VDATo16(2,2) = 9;  VDATo16(3,2) = 10  
 VDATo16(1,3) = 11; VDATo16(2,3) = 12; VDATo16(3,3) = 13  
 VDATo16(1,4) = 14; VDATo16(2,4) = 15; VDATo16(3,4) = 16  

 l = 0
 do nindex = m_npar, nwidth, np_npar
   l = l+1

   do ij = 1, (DIM+1)*(DIM+1)
     i = i_ij_lkp(ij)
     j = j_ij_lkp(ij)

     do mn  = 1, mnAFTot
     do mnp = 1, mnAFTot

       mnmnp = (mn-1)*mnAFTot + mnp
       DMult_imn(mnmnp,ij,l) = 1._dp

       if (i==yyy) DMult_imn(mnmnp,ij,l) = DMult_imn(mnmnp,ij,l)*iim*mAF(mn)
       if (i==zzz) DMult_imn(mnmnp,ij,l) = DMult_imn(mnmnp,ij,l)*iim*nAF(mn,l)

       if (j==yyy) DMult_imn(mnmnp,ij,l) = DMult_imn(mnmnp,ij,l)*iim*mpAF(mnp)
       if (j==zzz) DMult_imn(mnmnp,ij,l) = DMult_imn(mnmnp,ij,l)*iim*npAF(mnp,l)

     enddo
     enddo
   enddo
 enddo


 allocate(mEq_tmp(mnAFTot*mnAFTot),nEq_tmp(mnAFTot*mnAFTot))

 mnEqTot = 0;  mEq_tmp = 0;  nEq_tmp = 0;  mneq0_lkp = 0
 do mn  = 1, mnAFTot
 do mnp = 1, mnAFTot
   
   Dm = mAF(mn)+mpAF(mnp) ! mEq = -Dm
   Dn = nAF(mn,1)+npAF(mnp,1) ! nEq = -Dn

   NewEntry = .TRUE.
   ! check if already in table
   do i = 1, mnEqTot
     if ( (Dm == -mEq_tmp(i)).and.(Dn == -nEq_tmp(i)) ) then
       NewEntry = .False.
       mneq0_lkp((mn-1)*mnAFTot + mnp) = i
     endif
   enddo

   if (NewEntry) then
    mnEqTot = mnEqTot + 1
    mEq_tmp(mnEqTot) = -Dm
    nEq_tmp(mnEqTot) = -Dn
    mneq0_lkp((mn-1)*mnAFTot + mnp) = mnEqTot
   endif

 enddo
 enddo

 allocate(mEq(mnEqTot),nEq(mnEqTot))

 mEq = 0; nEq = 0
 do i = 1, mnEqTot
   mEq(i) = mEq_tmp(i)
   nEq(i) = nEq_tmp(i)
 enddo


 if (me.eq.0) &
   print *, 'Number of equilibrium Fourier modes: ', mnEqTot


 mnEqDTTot = 0;  mEq_tmp = 0;  nEq_tmp = 0;  mneq0DT_lkp = 0
 do mn  = 1, mnDTTot
 do mnp = 1, mnDTTot
   
   Dm = mDT(mn)+mpDT(mnp) ! mEq = -Dm
   Dn = nDT(mn,1)+npDT(mnp,1) ! nEq = -Dn

   NewEntry = .TRUE.
   ! check if already in table
   do i = 1, mnEqDTTot
     if ( (Dm == -mEq_tmp(i)).and.(Dn == -nEq_tmp(i)) ) then
       NewEntry = .False.
       mneq0DT_lkp((mn-1)*mnDTTot + mnp) = i
     endif
   enddo

   if (NewEntry) then
     mnEqDTTot = mnEqDTTot + 1
     mEq_tmp(mnEqDTTot) = -Dm
     nEq_tmp(mnEqDTTot) = -Dn
     mneq0DT_lkp((mn-1)*mnDTTot + mnp) = mnEqDTTot
   endif

 enddo
 enddo


 if (me.eq.0) &
  print *, 'Number of equilibrium Fourier modes for  k||: ', mnEqDTTot

 allocate(mEqDT(mnEqDTTot),nEqDT(mnEqDTTot))

 mEqDT = 0; nEqDT = 0
 do i = 1, mnEqDTTot
   mEqDT(i) = mEq_tmp(i)
   nEqDT(i) = nEq_tmp(i)
 enddo

 deallocate(mEq_tmp,nEq_tmp)

 allocate(maxFEqAmp(mnEqTot))
 allocate(maxFEqAmpMax(mnEqTot))
 allocate(mnFFT_lkp(mnEqTot))
 allocate(mnFFT_DT_lkp(mnEqDTTot))


 iAF_lkp = 0;  iF_lkp = 0;  iA_lkp = 0
 do k  = 1, NDIM
 do kp = 1, NDIM
 do mn  = 1, mnAFTot
 do mnp = 1, mnAFTot

   AModified = 0
   select case (abs(mAF(mn)))
     case (1)
       if ((k==3).or.(k==4)) AModified = 1
     case (0,2)
       if ((k==1).or.(k==2)) AModified = 1
   endselect
   FModified = 0
   select case (abs(mpAF(mnp)))
     case (1)
       if ((kp==3).or.(kp==4)) FModified = 1
     case (0,2)
       if ((kp==1).or.(kp==2)) FModified = 1
   endselect
   iAF_lkp((mn-1)*mnAFTot + mnp,k,kp) = AModified*2 + FModified
   iA_lkp(mn,k)   = AModified
   iF_lkp(mnp,kp) = FModified

 enddo
 enddo
 enddo
 enddo



 mnFFT_lkp = 0
 do ieq  = 1, mnEqTot

   m = mEq(ieq)
   n = nEq(ieq)/nperiod
!   if (   ((abs(m)>njT/2).and.(abs(m)<=mEqMax))  &
!     .or. ((abs(n)>nkT/2).and.(abs(n)<=nEqMax)) ) then
!     print *, 'Angular resolution insufficient for the required number of Fourier modes'
!     print *, 'mEq = ',m, '    nEq = ',n
!     print *, 'njT/2=',njT/2, '    nkT/2=',nkT/2
!     stop
!   endif
   if (m>=0) then
     it = m + 1
   else
     it = njT - abs(m) + 1
   endif
   if (n>=0) then
     ip = n + 1
   else
     ip = nkT - abs(n) + 1
   endif

!   mnFFT_lkp(ieq) = itp_it_ip_lkp(it,ip)

   if ((abs(mEq(ieq)) > njT/2).or.(abs(nEq(ieq)/nperiod) > nkT/2)) then
     mnFFT_lkp(ieq) = 0
   else
     mnFFT_lkp(ieq) = itp_it_ip_lkp(it,ip)
   endif

 enddo

 j = 0
 do nindex = m_npar, nwidth, np_npar
  j = j+1
  mnFFT_DT_lkp = 0
  do mn = 1, mnEqDTTot

   m = mEqDT(mn)
   n = nEqDT(mn)/nperiod

   if (m>=0) then
     it = m + 1
   else
     it = njT - abs(m) + 1
   endif
   if (n>=0) then
     ip = n + 1
   else
     ip = nkT - abs(n) + 1
   endif

   if ((abs(mEqDT(mn)) > njT/2).or.(abs(nEqDT(mn)/nperiod) > nkT/2)) then
     mnFFT_DT_lkp(mn) = 0
   else
     mnFFT_DT_lkp(mn) = itp_it_ip_lkp(it,ip)
   endif

  enddo

 enddo

j = 0
 do nindex = m_npar, nwidth, np_npar
  j = j+1
  mnFFT_RHS_lkp = 0
  do mn = 1, mnDTTot

   m = mDT(mn)
   n = nDT(mn,j)/nperiod

   if (m-m_ant>=0) then
     it = m - m_ant + 1
   else
     it = njT - abs(m-m_ant) + 1
   endif
   if (n-n_ant/nperiod>=0) then
     ip = n - n_ant/nperiod + 1
   else
     ip = nkT - abs(n-n_ant/nperiod) + 1
   endif
   mnFFT_RHS_lkp(mn) = itp_it_ip_lkp(it,ip)

  enddo

 enddo


! print *, 'Lookup table'
! do ieq  = 1, mnEqTot
!   print *, ieq, mEq(ieq), nEq(ieq), mnFFT_lkp(ieq)
! enddo

  ! Locate edge modes (on the edge of the Fourier box)
  EdgeModes_lkp(:) = 0
  EdgeModeNo = 0


!
!  ! Remember, doesn't work for m/nAFTot = 2!!! 
!  if ((mAFTot>1).and.(nAFTot>1)) then
!    do n = 1, nAFTot
!    do m = 1, mAFTot
!      if ((m==1).or.(m==mAFTot).or.(n==1).or.(n==nAFTot)) then
!        EdgeModeNo = EdgeModeNo + 1
!        EdgeModes_lkp(EdgeModeNo) = (n-1)*mAFTot + m
!      endif
!    enddo
!    enddo
!  endif
!
!  if ((mAFTot>1).and.(nAFTot==1)) then
!    n = 1;   m = mAFTot;  EdgeModeNo = EdgeModeNo + 1
!    EdgeModes_lkp(EdgeModeNo) = (n-1)*mAFTot + m
!    n = 1;   m = 1;       EdgeModeNo = EdgeModeNo + 1
!    EdgeModes_lkp(EdgeModeNo) = (n-1)*mAFTot + m
!  endif
!
!  if ((mAFTot==1).and.(nAFTot>1)) then
!    n = 1;        m = 1;   EdgeModeNo = EdgeModeNo + 1
!    EdgeModes_lkp(EdgeModeNo) = (n-1)*mAFTot + m
!    n = nAFTot;   m = 1;   EdgeModeNo = EdgeModeNo + 1
!    EdgeModes_lkp(EdgeModeNo) = (n-1)*mAFTot + m
!  endif
!
!  print *, 'Edge modes:', EdgeModeNo
!  do mn = 1, EdgeModeNo
!    print *, 'm = ', mAFmn(EdgeModes_lkp(mn)), ',      n = ', nAFmn(EdgeModes_lkp(mn))
!  enddo

! print *, 'N:',nperiod, n_ant,  nEq(1:mnEqTot)
 do i = 1, mnEqTot
   if ((mod(nEq(i),nperiod) /= 0).or.(mod(n_ant+npAF(1,1),nperiod) /= 0)) then
     print *, 'Error: Periodicity condition is not satisfied.'
     stop
   endif
 enddo
! if (me.eq.0) print *, 'Number of eq. Fourier modes: ', mnEqTot


 l = 0
 do nindex = m_npar, nwidth, np_npar
   l = l+1

   do itp = 1, Njk
    do mn = 1, mnAFTot
      Exp_mn(itp,mn,l)   = exp(iim*(mAF(mn)*th_itp(itp) + nAF(mn,l)*ph_itp(itp)))
    enddo
    do mnp = 1, mnAFTot
      Exp_mnp(itp,mnp,l) = exp(iim*(mpAF(mnp)*th_itp(itp) + npAF(mnp,l)*ph_itp(itp)))
    enddo
   enddo

 enddo

 call RenormalizeEquilibrium

!!!Testing  pf = -1.0_dp
  pf = -1.0_dp
!  pf = 0.5_dp

  do eR = 1, NelR
  do gR = 1, NGaussR 
    s = s_gauss(eR,gR)

    sf_v(1,eR,gR)  = sfr(1,s)
    sf_v(2,eR,gR)  = sfr(2,s)
    sf_v(3,eR,gR)  = sfr(3,s)
    sf_v(4,eR,gR)  = sfr(4,s)

    dsf_v(1,eR,gR) = dsfr(1,s)
    dsf_v(2,eR,gR) = dsfr(2,s)
    dsf_v(3,eR,gR) = dsfr(3,s)
    dsf_v(4,eR,gR) = dsfr(4,s)
  enddo
  enddo

  s = s_nod(1)
  sf_v_A(1,1) = 0
  sf_v_A(2,1) = 0
  sf_v_A(3,1) = sfr(3,s)
  sf_v_A(4,1) = sfr(4,s)
  do eR = 2, Nmesh
    s = s_nod(eR)

    sf_v_A(1,eR)  = sfr(1,s)
    sf_v_A(2,eR)  = sfr(2,s)
    sf_v_A(3,eR)  = sfr(3,s)
    sf_v_A(4,eR)  = sfr(4,s)
  enddo


  do eR = 2, Nmesh
    s = s_nod(eR)

    dsf_v_A(1,eR) = dsfr(1,s)
    dsf_v_A(2,eR) = dsfr(2,s)
    dsf_v_A(3,eR) = dsfr(3,s)
    dsf_v_A(4,eR) = dsfr(4,s)
  enddo
  dsf_v_A(1,1) = 0
  dsf_v_A(2,1) = 0
  dsf_v_A(3,1) = 0
  dsf_v_A(4,1) = 0


  ! Display of all the profiles
  if (me.eq.0) then

  ! Density
    do Gx = 0, Gwidth
      s = Gx/real(Gwidth)
      DensDpl(Gx+1) = PlasmaDensity(s)
      DensDpl(Gx+1) = FastPartDensity(s)
    enddo
    TitleChar = 'Density profile:'
    TitleChar = 'Hot density profile:'
    call ProfileDisplay(TitleChar,'n ',UseHotDensCfs,DensDpl,np0,np1,np2,np3,np4,np5)

  ! Ion Temperature
    do Gx = 0, Gwidth
      s = Gx/real(Gwidth)
      TempDpl(Gx+1) = IonTemperature(s)
      TempDpl(Gx+1) = FastTemperature(s)
    enddo
    TitleChar = 'Ion Temperature:'
    TitleChar = 'Hot ion Temperature:'
    call ProfileDisplay(TitleChar,'Ti',UseHotTempCfs,TempDpl,Tip0,Tip1,Tip2,Tip3,Tip4,Tip5)

  ! Electron Temperature
    do Gx = 0, Gwidth
      s = Gx/real(Gwidth)
      TempDpl(Gx+1) = ElectronTemperature(s)
    enddo
    TitleChar = 'Electron Temperature:'
    call ProfileDisplay(TitleChar,'Te',UseTempCoeffs,TempDpl,Tep0,Tep1,Tep2,Tep3,Tep4,Tep5)

    print *, " "

  end if

! print *, 'Fourier modes'
! do mn  = 1, mnAFTot
! do mnp = 1, mnAFTot
!   i = mneq0_lkp(mn,mnp)
!   print *, mAFmn(mn), nAFmn(mn), mpAFmn(mnp), npAFmn(mnp), i, mEq(i), nEq(i)
! enddo
! enddo
!
! print *, 'Eq. modes'
! do i = 1, mnEqTot
!   print *, mEq(i), nEq(i)
! enddo

 contains

 function sfr(i,s)  ! Multiplying factor
   real(dp) sfr, s
   integer  i
   select case(i)
     case (1);   sfr = s**pf !1./sqrt(s) !1./s
     case (2);   sfr = s**pf !1./sqrt(s) !1./s
     case (3);   sfr = 1._dp
     case (4);   sfr = 1._dp
   end select
 end function

 function dsfr(i,s)  ! Multiplying factor
   real(dp) dsfr, s
   integer  i
   select case(i)
     case (1);   dsfr = pf*s**(pf-1)  !-0.5/(s**1.5) !-1./s/s
     case (2);   dsfr = pf*s**(pf-1)  !-0.5/(s**1.5) !-1./s/s
     case (3);   dsfr = 0._dp
     case (4);   dsfr = 0._dp
   end select
 end function



 function PlasmaDensity(s)
 real(dp) s, s_mod, PlasmaDensity
 integer int_a, int_b, int_c
 real(dp) dsi

 PlasmaDensity = 1e-10

 if(UseDensCoeffs.eq.1.or.UseDensCoeffs.eq.4) then

   s_mod = s/ArtBorder

   if (s < ArtBorder) then
     PlasmaDensity = np0 + np1*s_mod + np2*s_mod**2 + np3*s_mod**3 + np4*s_mod**4 + np5*s_mod**5
   end if

 elseif(UseDensCoeffs.eq.2) then

!   PlasmaDensity = 1.0_dp - s !*0.9_dp

   if (sqrt(s)<0.95_dp) &
!      PlasmaDensity = 1.0_dp
     PlasmaDensity = (1.0_dp - (sqrt(s)/0.95_dp)**2.0_dp)**1.5_dp

!   if (sqrt(s)<0.95_dp) &
!     PlasmaDensity = 1.0_dp

!     PlasmaDensity = (1.0_dp - (sqrt(s)/0.95_dp)**2.0_dp)**1.5_dp

!     PlasmaDensity = (1.0_dp - (s/0.98_dp)**0.8_dp)**1.3_dp
!     PlasmaDensity = (1.0_dp - (s/0.9_dp)**1.5_dp)**3.5_dp
!%n = (1.0 - (r/0.95).^2.0).^1.5;


!   PlasmaDensity = 1.0_dp - s !*0.9_dp

!   if (s<0.9_dp) &
!     PlasmaDensity = 1.0_dp - s*s/0.81_dp


!    PlasmaDensity = 1.0_dp - s*s


!   PlasmaDensity = 1.0000 + s*(-1.0332) + s**2*0.0313 + s**3*0.0016 + s**4*0.0003
 ! TEST case: n = 1-s_LION^2,  QAS 2D, Soloviev eq.


!    fnp(eR,gR) = 1.0000 + s*(-1.4574) + s**2*0.5430 + s**3*(-0.0943) + s**4*0.0086 + s**5*0.0001
    ! ^^^^--- pressure profile from LION (p.ex. TEST_GAP equilibrium (n = 1-s_LION^2 ) )

 elseif(UseDensCoeffs.eq.3) then

 int_a = 1
 int_b = nif

 do while (int_b-int_a>1)
   int_c = (int_a+int_b)/2
   if (s==nT_fits(int_a,1)) then
    int_b = int_a+1
   elseif (s==nT_fits(int_b,1)) then
    int_a = int_b-1
   elseif ((s-nT_fits(int_a,1))*(s-nT_fits(int_c,1))<0.0) then
    int_b = int_c
   else
    int_a = int_c
   end if
 enddo

   dsi=(s-nT_fits(int_a,1))/(nT_fits(int_b,1)-nT_fits(int_a,1))
   PlasmaDensity = dsi*nT_fits(int_b,2) + (1.-dsi)*nT_fits(int_a,2)
   PlasmaDensity = PlasmaDensity/nT_fits(1,2)

 end if

 end function


 function MinorityDensity(s)
 real(dp) s, s_mod, MinorityDensity
 integer is
 real(dp) dsi

 MinorityDensity=1.e-10

 if(UseDensCoeffs.le.3) then
      MinorityDensity=PlasmaDensity(s)
 else
   is=floor(s*niT)
   is=min(is,niT-1)
   dsi=s*niT-is
   MinorityDensity = dsi*thermal(is+1,2) + (1.-dsi)*thermal(is,2)
   MinorityDensity = MinorityDensity/n0_SI/ifr2
 endif

 end function
      

 function FastPartDensity(s)
 real(dp) s, s_mod, FastPartDensity
 integer is
 real(dp) dsi

 FastPartDensity = 0.

 if(UseHotDensCfs.eq.1) then

   s_mod = s/ArtBorder

   if (s < ArtBorder) then
     FastPartDensity = nh0 + nh1*s_mod + nh2*s_mod**2 + nh3*s_mod**3 + nh4*s_mod**4 + nh5*s_mod**5
   end if

 elseif(UseHotDensCfs.eq.2)then

   if (sqrt(s)<0.95_dp) &
     FastPartDensity = (1.0_dp - (sqrt(s)/0.95_dp)**2.0_dp)**1.5_dp

 else
   is=floor(s*niT)
   is=min(is,niT-1)
   dsi=s*niT-is
   FastPartDensity = dsi*moments(is+1,1) + (1.-dsi)*moments(is,1)
   FastPartDensity = FastPartDensity/nh0_SI
 end if

 end function


! Modelling absorption profile in cold plasma
 function AbsorpProf(s)
  real(dp) AbsorpProf, s, r
  real(dp), parameter :: r1 = 0.3_dp,  r2 = 0.4_dp
  
!  r = sqrt(s)
  AbsorpProf = 1.0_dp ! - s
  
 end function

! Ion and electron temperature profiles
 function IonTemperature(s)
  real(dp) s, s_mod, IonTemperature
  integer int_a, int_b, int_c
  real(dp) dsi

  IonTemperature = 1e-10

  if(UseTempCoeffs.eq.1.or.UseTempCoeffs.eq.4) then

   s_mod = s/ArtBorder

   if (s < ArtBorder) then
    IonTemperature = Tip0 + Tip1*s_mod + Tip2*s_mod**2 + Tip3*s_mod**3 + Tip4*s_mod**4 + Tip5*s_mod**5
   end if

  elseif(UseTempCoeffs.eq.2) then

   IonTemperature = 1.0_dp

  elseif(UseTempCoeffs.eq.3) then

 int_a = 1
 int_b = nif

 do while (int_b-int_a>1)
   int_c = (int_a+int_b)/2
   if (s==nT_fits(int_a,1)) then
    int_b = int_a+1
   elseif (s==nT_fits(int_b,1)) then
    int_a = int_b-1
   elseif ((s-nT_fits(int_a,1))*(s-nT_fits(int_c,1))<0.0) then
    int_b = int_c
   else
    int_a = int_c
   end if
 enddo

   dsi=(s-nT_fits(int_a,1))/(nT_fits(int_b,1)-nT_fits(int_a,1))
   IonTemperature = dsi*nT_fits(int_b,3) + (1.-dsi)*nT_fits(int_a,3)
   IonTemperature = IonTemperature/nT_fits(1,3)

  end if

 end function

 function ElectronTemperature(s)
  real(dp) s, s_mod, ElectronTemperature
  integer int_a, int_b, int_c
  real(dp) dsi

  ElectronTemperature = 1e-10

  if(UseTempCoeffs.eq.1.or.UseTempCoeffs.eq.4) then

   s_mod = s/ArtBorder

   if (s < ArtBorder) then
    ElectronTemperature = Tep0 + Tep1*s_mod + Tep2*s_mod**2 + Tep3*s_mod**3 + Tep4*s_mod**4 + Tep5*s_mod**5
   end if

  elseif (UseTempCoeffs.eq.2) then

   ElectronTemperature = 1.0_dp

  elseif (UseTempCoeffs.eq.3) then

 int_a = 1
 int_b = nif

 do while (int_b-int_a>1)
   int_c = (int_a+int_b)/2
   if (s==nT_fits(int_a,1)) then
    int_b = int_a+1
   elseif (s==nT_fits(int_b,1)) then
    int_a = int_b-1
   elseif ((s-nT_fits(int_a,1))*(s-nT_fits(int_c,1))<0.0) then
    int_b = int_c
   else
    int_a = int_c
   end if
 enddo

   dsi=(s-nT_fits(int_a,1))/(nT_fits(int_b,1)-nT_fits(int_a,1))
   ElectronTemperature = dsi*nT_fits(int_b,3) + (1.-dsi)*nT_fits(int_a,3)
   ElectronTemperature = ElectronTemperature/nT_fits(1,3)

  end if

 end function


 function MinorityTemperature(s)
 real(dp) s, s_mod, MinorityTemperature
 integer is
 real(dp) dsi

 MinorityTemperature=1.e-10

 if(UseTempCoeffs.le.3) then
      MinorityTemperature=IonTemperature(s)
 else
   is=floor(s*niT)
   is=min(is,niT-1)
   dsi=s*niT-is
   MinorityTemperature = dsi*thermal(is+1,3) + (1.-dsi)*thermal(is,3)
   MinorityTemperature = MinorityTemperature/Ti0
 endif

 end function


 function FastTemperature(s)
  real(dp) s, s_mod, FastTemperature
  integer is
  real(dp) dsi

  FastTemperature = 0

  if(UseHotTempCfs.eq.1) then

   s_mod = s/ArtBorder

   if (s < ArtBorder) then
    FastTemperature = Thp0 + Thp1*s_mod + Thp2*s_mod**2 + Thp3*s_mod**3 + Thp4*s_mod**4 + Thp5*s_mod**5
   end if

  elseif(UseHotTempCfs.eq.2)then

   FastTemperature = 1.0_dp

  else
   is=floor(s*niT)
   is=min(is,niT-1)
   dsi=s*niT-is
   FastTemperature = dsi*moments(is+1,2) + (1.-dsi)*moments(is,2)
   FastTemperature = FastTemperature/Th0
  end if

 end function

 function FastAnisotropy(s)
  real(dp) s, s_mod, FastAnisotropy
  integer is
  real(dp) dsi

  FastAnisotropy = 0

  if(UseHotTempCfs.eq.1)then

   FastAnisotropy = An0*(Anp0 + Anp1*s_mod + Anp2*s_mod**2 + Anp3*s_mod**3 + Anp4*s_mod**4 + Anp5*s_mod**5)

  elseif(UseHotTempCfs.eq.2)then

   FastAnisotropy = 1.0_dp

  else
   is=floor(s*niT)
   dsi=s*niT-is
   FastAnisotropy = dsi*moments(is+1,3) + (1.-dsi)*moments(is,3)
  end if

 end function


end subroutine
!-------------------------------------------
subroutine ProfileDisplay(profile_title,axis_name,UseCoeffs,DplData,prf0,prf1,prf2,prf3,prf4,prf5)
 real(dp) prf0, prf1, prf2, prf3, prf4, prf5, DplData(Gwidth+1)
 real(dp) prfmax
 character(50) profile_title
 character(2) axis_name
 integer UseCoeffs
 integer  Gx, Gy


    print *, " "
    print *, profile_title
    if (UseCoeffs.eq.1.or.axis_name.eq.'Te') then

      write(*,'(a1,a2,a5,f4.1)',advance='no') ' ', axis_name, '(s) =', prf0

      if (prf1>=0) then
        write(*,'(a2,f4.1)',advance='no') ' +', prf1
      else
        write(*,'(a2,f4.1)',advance='no') ' -', -prf1
      endif

      if (prf2>=0) then
        write(*,'(a4,f4.1)',advance='no') ' s +', prf2
      else
        write(*,'(a4,f4.1)',advance='no') ' s -', -prf2
      endif

      if (prf3>=0) then
        write(*,'(a6,f4.1)',advance='no') ' s^2 +', prf3
      else
        write(*,'(a6,f4.1)',advance='no') ' s^2 -', -prf3
      endif

      if (prf4>=0) then
        write(*,'(a6,f4.1)',advance='no') ' s^3 +', prf4
      else
        write(*,'(a6,f4.1)',advance='no') ' s^3 -', -prf4
      endif

      if (prf5>=0) then
        write(*,'(a6,f4.1)',advance='no') ' s^4 +', prf5
      else
        write(*,'(a6,f4.1)',advance='no') ' s^4 -', -prf5
      endif

      write(*,'(a4)',advance='yes') ' s^5'

      write(*,'(a8,f4.1)') ' Border:', ArtBorder

    elseif(UseCoeffs.eq.2)then
      print *, 'Arbitrary profile'
    else
      prfmax=0.
      do Gx=1,Gwidth
         if(DplData(Gx).gt.prfmax)prfmax=DplData(Gx)
      enddo
      write(*,'(a,f5.2,a,f5.2)')' Profile read from moments.in, axis = ',DplData(1),',max = ',prfmax
      DplData=DplData/prfmax   
    end if
    print *, " "


    write(*,'(a2,a1)',advance='yes') axis_name, "^"   
 
    do Gy = Gheight, 1, -1
      if (DplData(1)<(Gy+0.5)/real(Gheight).and.DplData(1)>(Gy-0.5)/real(Gheight)) then
        write(*,'(a3)',advance='no') "  *"
      else
        write(*,'(a3)',advance='no') "  |"
      end if
      do Gx = 1, Gwidth
        if (DplData(Gx+1)<(Gy+0.5)/real(Gheight).and.DplData(Gx+1)>(Gy-0.5)/real(Gheight)) then
          write(*,'(a1)',advance='no') "*"
        else
          write(*,'(a1)',advance='no') " "
        end if
      enddo
      write(*,'(a1)',advance='yes') " "
    enddo

    write(*,'(a3)',advance='no') "   "    
    do Gx = 1, Gwidth
      if (DplData(Gx+1)<0.5/real(Gheight).and.DplData(Gx+1)>-0.5/real(Gheight)) then
        write(*,'(a1)',advance='no') "*"
      else
        write(*,'(a1)',advance='no') "-"
      end if
    enddo
    write(*,'(a1)',advance='yes') ">"

    write(*,'(a3)',advance='no') "   "    
    do Gx = 1, Gwidth
      write(*,'(a1)',advance='no') " "
    enddo
    write(*,'(a1)',advance='yes') "s"

end subroutine
!-------------------------------------------
subroutine Max_Resp(MR) ! Searches for the maximum response
real(dp) mmax
integer i, im, MR

  mmax = Response(0);  im = 0
  do i = 0, OmSteps
    if (Response(i)>mmax) then
      mmax = Response(i);  im = i
    endif
  enddo
  MR = im

end subroutine
!-------------------------------------------
subroutine CylinderAlvfenResonance
 integer, parameter :: fn = 2
 integer is,it,ip,itp
 real(dp) Bavg, k_par, v_a(niT)
 integer mn


  do is = 1,niT-1

    call InterpEqSurf(is)

    Bavg = 0._dp
    do it = 1,njT
    do ip = 1,nkT
      itp = itp_it_ip_lkp(it,ip)
      Bavg = Bavg + tbmod(itp,1)
    enddo
    enddo
    Bavg = Bavg / njT / nkT
    
    if (CylMetric) Bavg = 1.0_dp

    do mn = 1, mnAFTot
      k_par = (nAF(mn,1) + mAF(mn)*iota(is))/MajRadius
      v_a(is)   = cA_SI * Bavg/sqrt(dens_T(is))
      om_a(is,mn)  = v_a(is) * abs(k_par)
    enddo

  enddo


  Bavg = 0._dp
  do it = 1,njT
  do ip = 1,nkT
    itp = itp_it_ip_lkp(it,ip)
    Bavg = Bavg + tbmod(itp,2)
  enddo
  enddo
  Bavg = Bavg / njT / nkT
    
  if (CylMetric) Bavg = 1.0_dp

  do mn = 1, mnAFTot
    k_par = (nAF(mn,1) + mAF(mn)*iota(niT))/MajRadius
    v_a(niT)   = cA_SI * Bavg/sqrt(dens_T(niT))
    om_a(niT,mn)  = v_a(niT) * abs(k_par)
  enddo


  if (me.eq.0) then
    open(fn, file='AlfvRes.dat',form='unformatted')
    write(fn) niT
    write(fn) s_T
    write(fn) v_a
    write(fn) iota
    close(fn)
  end if

end subroutine


! ::::::::: FEM :::::::::
!-------------------------------------------
function Psi_al(p,a,e,x, Modified)   ! Functions used for FEM, a,e -- local indices, x -- global coordinate.
                                     ! d^p/dx psi(a,x)  on the element e
 real(dp) x, ksi, Psi_al,  s,a1,b1,a2,b2,h
 integer p,a, e, Modified

  Psi_al = 0._dp
  ksi = 0._dp
  if ((x>=s_nod(e)).and.(x<=s_nod(e+1))) then
    ksi = (2._dp*x-s_nod(e)-s_nod(e+1))/(s_nod(e+1)-s_nod(e))
  else
    return
  endif

!  if ((Modified == 1).and.((e==1).and.((a==1).or.(a==2)))) then
!
!    h = s_nod(e+1)-s_nod(e)
!    
!    a1 = -1.5/h**0.5
!    b1 =  0.5/h**1.5
!    a2 = -2.5/h**1.5
!    b2 =  1.5/h**2.5
!
!    select case (p)
!      case (0) ! basis function
!        if (a==1) Psi_al = 0.25_dp *( ksi*ksi*ksi - 3*ksi + 2)
!        if (a==2) Psi_al = 0.125_dp*( ksi*ksi*ksi - ksi*ksi - ksi + 1)*(s_nod(e+1)-s_nod(e))
!        if (a==3) Psi_al = 0.25_dp *(-ksi*ksi*ksi + 3*ksi + 2)
!        if (a==4) Psi_al = 0.125_dp*( ksi*ksi*ksi + ksi*ksi - ksi - 1)*(s_nod(e+1)-s_nod(e))
!      case (1) ! derivative of the basis function
!        if (a==1) Psi_al = 0.25_dp* ( 3*ksi*ksi - 3)*2/(s_nod(e+1)-s_nod(e))
!        if (a==2) Psi_al = 0.125_dp*( 3*ksi*ksi - 2*ksi - 1)*2
!        if (a==3) Psi_al = 0.25_dp* (-3*ksi*ksi + 3)*2/(s_nod(e+1)-s_nod(e))
!        if (a==4) Psi_al = 0.125_dp*( 3*ksi*ksi + 2*ksi - 1)*2
!    endselect
!
!
!  else

    select case (p)
      case (0) ! basis function
        if (a==1) Psi_al = 0.25_dp *( ksi*ksi*ksi - 3*ksi + 2)
        if (a==2) Psi_al = 0.125_dp*( ksi*ksi*ksi - ksi*ksi - ksi + 1)*(s_nod(e+1)-s_nod(e))
        if (a==3) Psi_al = 0.25_dp *(-ksi*ksi*ksi + 3*ksi + 2)
        if (a==4) Psi_al = 0.125_dp*( ksi*ksi*ksi + ksi*ksi - ksi - 1)*(s_nod(e+1)-s_nod(e))
      case (1) ! derivative of the basis function
        if (a==1) Psi_al = 0.25_dp* ( 3*ksi*ksi - 3)*2/(s_nod(e+1)-s_nod(e))
        if (a==2) Psi_al = 0.125_dp*( 3*ksi*ksi - 2*ksi - 1)*2
        if (a==3) Psi_al = 0.25_dp* (-3*ksi*ksi + 3)*2/(s_nod(e+1)-s_nod(e))
        if (a==4) Psi_al = 0.125_dp*( 3*ksi*ksi + 2*ksi - 1)*2
    endselect

!   endif

  return


!    select case (p)
!      case (0) ! basis function
!        if (a==1) Psi_al = 1.!(s**1.5 + a1*s*s + b1*s*s*s)
!        if (a==2) Psi_al = 0.125_dp*( ksi*ksi*ksi - ksi*ksi - ksi + 1)*(s_nod(e+1)-s_nod(e))
!     !Psi_al = s**0.5 + a2*s*s + b2*s*s*s
!      case (1) ! derivative of the basis function
!        if (a==1) Psi_al = 0.!(1.5*s**0.5 + 2*a1*s + 3*b1*s*s)
!        if (a==2) Psi_al = 0.125_dp*( 3*ksi*ksi - 2*ksi - 1)*2 !Psi_al = 0.5/s**0.5 + 2*a1*s + 3*b1*s*s
!    endselect


!    sh = 0.5_dp*(1._dp+ksi)
!    ah = -1.5
!    bh = 0.5
!    select case (p)
!      case (0) ! basis function
!        Psi_al =  sqrt(sh) + ah*sh + bh*sh*sh
!      case (1) ! derivative of the basis function
!        Psi_al = (0.5/sqrt(sh) + ah + 2*bh*sh)/sqrt(s_nod(e+1)-s_nod(e))
!    endselect



 end function
!-------------------------------------------


!-------------------------------------------
 subroutine Init_jext3D(eR)
 integer  eR,gR, itp,itpr,ip,it
 real(dp) sigma_s, r, s, ksi, ksi2, thG,thB, phB, sigma_theta, &
   dsigma_theta, sigma_phi, dsigma_phi
 complex(dp) sigma, dthsigma, dphsigma


  jext_3D = 0
  do itpr = 1, Njkr

    gR = gR_itpr_lkp(itpr)
    itp = itp_itpr_lkp(itpr)
    ip = ip_itpr_lkp(itpr)
    it = it_itpr_lkp(itpr)

    s = s_gauss(eR,gR)
!    r = sqrt(s)
    r = s

    sigma_s = 0
    ksi = 2._dp*(s-s_ant_min)/(s_ant_max-s_ant_min)-1
    if ( (s>s_ant_min).and.(s<s_ant_max) ) &
      sigma_s = (1 - ksi*ksi)*(1 - ksi*ksi)
!    sigma_s = 1


   select case (ant_type)

 case (0)  ! 0 -- helical antenna, one mode
    sigma = sigma_s * exp(iim*(m_ant*th_itp(itp) + n_ant*ph_itp(itp)))
    dthsigma = iim*m_ant * sigma
    dphsigma = iim*n_ant * sigma


 case (1)  ! 1 -- Directionaly localised antenna

    thB = th_itp(itp)
    phB = ph_itp(itp)
    sigma_theta=0
    dsigma_theta=0
    dsigma_phi=0

    if (theta_ant1<theta_ant2) then

      if (m_shape==0) then ! Mirror case
        ksi = 2*(phB-theta_ant1)/(theta_ant2-theta_ant1)-1
      elseif (abs(m_shape)>abs(n_shape).and.n_shape/=0) then
        ksi = 2*mod(thB+n_shape/real(m_shape)*phB-theta_ant1+2*pi &
         * (abs(n_shape/(real(m_shape)))+1),2*n_shape/real(m_shape)*pi) &
           /(theta_ant2-theta_ant1)-1
      else
        ksi = 2*mod(thB+n_shape/real(m_shape)*phB-theta_ant1+2*pi &
         * abs(n_shape/(real(m_shape))),2*pi) &
           /(theta_ant2-theta_ant1)-1
      end if

      sigma_theta = (1 - ksi*ksi)*(1 - ksi*ksi)

      if (ksi>-1.0.and.ksi<1.0) then
        if (m_shape==0) then
          dsigma_theta = 0
          dsigma_phi = -8*ksi*(1 - ksi*ksi)/(theta_ant2-theta_ant1)
        else
          dsigma_theta = -8*ksi*(1 - ksi*ksi)/(theta_ant2-theta_ant1)
          dsigma_phi = -8*n_shape/real(m_shape)*ksi*(1 - ksi*ksi) &
            &  /(theta_ant2-theta_ant1)
        end if
      else
        sigma_theta=0
        dsigma_theta=0
        dsigma_phi=0
      end if

    else

      if (m_shape==0) then !Mirror case

        if (phB>theta_ant1) then
          ksi = 2*(phB-theta_ant1)/(theta_ant2-theta_ant1+2*pi)-1
        else
          ksi = 2*(phB-theta_ant1+2*pi)/(theta_ant2-theta_ant1+2*pi)-1
        end if

        if (ksi>-1.0.and.ksi<1.0) then
          sigma_theta = (1 - ksi*ksi)*(1 - ksi*ksi)
          dsigma_theta = 0
          dsigma_phi = -8*ksi*(1 - ksi*ksi)/(theta_ant2-theta_ant1+2*pi)
        end if

      else

        if (abs(m_shape)>abs(n_shape)/nperiod.and.n_shape/=0) then
          ksi = 2*mod(thB+n_shape/real(m_shape)*phB-theta_ant1+2*pi &
           * (2+abs(n_shape/(real(m_shape)))),2*n_shape/real(m_shape)*pi) &
             /(theta_ant2-theta_ant1+2*pi)-1
        else
          ksi = 2*mod(thB+n_shape/real(m_shape)*phB-theta_ant1+2*pi &
           * (1+abs(n_shape/(real(m_shape)))),2*pi) &
             /(theta_ant2-theta_ant1+2*pi)-1
        end if

        if (ksi>-1.0.and.ksi<1.0) then
          sigma_theta = (1 - ksi*ksi)*(1 - ksi*ksi)
          dsigma_theta = -8*ksi*(1 - ksi*ksi)/(theta_ant2-theta_ant1+2*pi)
          dsigma_phi = -8*n_shape/real(m_shape)*ksi*(1 - ksi*ksi) &
            &  /(theta_ant2-theta_ant1+2*pi)
        end if


      end if

    end if

    sigma    = sigma_s * exp(iim*(m_ant*th_itp(itp)+n_ant*ph_itp(itp))) * sigma_theta
    dthsigma = sigma_s * exp(iim*(m_ant*th_itp(itp)+n_ant*ph_itp(itp))) &
               & * dsigma_theta + iim*m_ant * sigma
    dphsigma = sigma_s * exp(iim*(m_ant*th_itp(itp)+n_ant*ph_itp(itp))) &
               & * dsigma_phi + iim*n_ant * sigma


 case(2)  ! 2 -- Localised antenna

    thB = th_itp(itp)
    phB = ph_itp(itp)
    sigma_theta=0
    sigma_phi=0
    dsigma_theta=0
    dsigma_phi=0

    ksi  = 2*mod(thB-theta_ant1+2*pi,2*pi)/mod(theta_ant2-theta_ant1+2*pi,2*pi)-1
    ksi2 = 2*mod(phB-theta_ant3+2*pi,2*pi)/mod(theta_ant4-theta_ant3+2*pi,2*pi)-1

    if (ksi>-1.0.and.ksi<1.0.and.ksi2>-1.0.and.ksi2<1.0) then
      sigma_theta  = (1 - ksi*ksi)*(1 - ksi*ksi)
      sigma_phi    = (1 - ksi2*ksi2)*(1 - ksi2*ksi2)
      dsigma_theta = -8*ksi*(1 - ksi*ksi)/mod(theta_ant2-theta_ant1+2*pi,2*pi)
      dsigma_phi   = -8*ksi2*(1 - ksi2*ksi2)/mod(theta_ant4-theta_ant3+2*pi,2*pi)
    end if

    sigma    = sigma_s * sigma_theta * sigma_phi
    dthsigma = sigma_s * dsigma_theta * sigma_phi
    dphsigma = sigma_s * sigma_theta * dsigma_phi

   end select



  if (abs(bjac(gR,itp))<1e-10) &
    print *, 'BJAC: ', eR, gR, itp, bjac(gR,itp)


      jext_3D(itpr,1) = 0
      jext_3D(itpr,2) = jext_3D(itpr,2) + &
   ( mTL(3,2,gR,itp)*dthsigma - mTL(2,2,gR,itp)*seTL(2,eR,gR)*dphsigma ) / bjac(gR,itp)
      jext_3D(itpr,3) = jext_3D(itpr,3) + &
   ( mTL(3,3,gR,itp)*dthsigma - mTL(2,3,gR,itp)*dphsigma ) / bjac(gR,itp)


!      jext_3D(itpr,1) = 0
!      jext_3D(itpr,2) = 1
!      jext_3D(itpr,3) = 0


! Old expression for div-free j_ext
!      jext_3D(itpr,1) = 0
!      jext_3D(itpr,2) = jext_3D(itpr,2) + (- mTU(1,1,gR,itp)* ( mTU(2,3,gR,itp)*iim*ma*sigma + &
!                                             mTU(3,3,gR,itp)*iim*na*sigma ))
!      jext_3D(itpr,3) = jext_3D(itpr,3) + ( mTL(3,3,gR,itp)*iim*ma*sigma -                     &
!                                            mTL(2,3,gR,itp)*iim*na*sigma ) / bjac(gR,itp)



!   enddo

!    Just one mode in the antenna (m_ant,n_ant)
!
!    sigma = sigma_s * exp(iim*(m_ant*th_itp(itp) + n_ant*ph_itp(itp)))
!
!    jext_3D(itpr,1) = 0
!    jext_3D(itpr,2) = - mTU(1,1,gR,itp)* ( mTU(2,3,gR,itp)*iim*m_ant*sigma + &
!                        mTU(3,3,gR,itp)*iim*n_ant*sigma )
!    jext_3D(itpr,3) = ( mTL(3,3,gR,itp)*iim*m_ant*sigma - &
!                        mTL(2,3,gR,itp)*iim*n_ant*sigma ) / bjac(gR,itp)

  enddo ! itpr


 end subroutine
!-------------------------------------------


 end module math


! ==================================================================
subroutine Spline(x,y,N,y2)
 implicit none
 integer, parameter :: dp = SELECTED_REAL_KIND(10)

 integer N
 real(dp), dimension(:), allocatable :: U
! real(dp) U(1000)
 real(dp) x(N), y(N), y2(N)
 real(dp) sig, p, qn,un

 integer i,k



 allocate(U(N))

 y2(1) = 0._dp
 U(1)  = 0._dp

 do i = 2,N-1
   sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
   p   = sig*y2(i-1) + 2._dp
   y2(i) = (sig-1._dp)/p
   U(i)  = (6._dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
           /(x(i+1)-x(i-1))-sig*U(i-1))/p
 enddo

 qn = 0._dp
 un = 0._dp

 y2(N) = un
 do k = N-1,1,-1
   y2(k) = y2(k)*y2(k+1)+U(k)
 enddo

 deallocate(U)

end subroutine
! ==================================================================
subroutine SplInt(xa,ya,y2a,N,x,y,y1)
 implicit none
 integer, parameter :: dp = SELECTED_REAL_KIND(10)

 real(dp) x,y,y1
 real(dp) xa(N), ya(N), y2a(N)
 integer  klo,khi,k, N
 real(dp) a,b,h


 klo = 1
 khi = N

 do while ((khi-klo) > 1)
   k = (khi + klo)/2
   if (xa(k)>x) then
     khi = k
   else
     klo = k
   endif
 enddo

 h = xa(khi) - xa(klo)
 a = (xa(khi)-x)/h
 b = (x-xa(klo))/h

 y = a*ya(klo) + b*ya(khi) + ((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*(h*h)/6._dp

 y1 = (ya(khi)-ya(klo))/h - (3._dp*a*a-1)*h*ya(klo)/6._dp + (3._dp*b*b-1)*h*ya(khi)/6._dp

end subroutine
! ==================================================================
! Interpolates the TERPSICHORE equilibrium from the T. angular grid
! to a finer grid. Needed to be able to run LEMan for higher angular
! harmonics (e.g. in the ICRF domain)
! Uses FFT for the interpolation
subroutine InterpAnglesEq(ism,njTO,nkTO,v3DTO,ni,nj,nk,v3D,forw,backw)
 implicit none

 integer*8 :: forw, backw

 integer, parameter :: dp = SELECTED_REAL_KIND(10)

 integer ism,njTO,nkTO,ni,nj,nk
 real(dp) v3DTO(njTO*nkTO,ni),v3D(nj*nk,2)

 complex(dp), dimension(:), allocatable :: XTO, YTO
 complex(dp), dimension(:), allocatable :: X, Y

 integer NjkTO, Njk
 integer itTO,ipTO,itpTO, m,n
 integer it,ip,itp
 integer is

 NjkTO = njTO*nkTO
 Njk = nj*nk
 allocate(XTO(NjkTO))
 allocate(YTO(NjkTO))
 allocate(X(Njk))
 allocate(Y(Njk))


 do is = 1, 2

! 1. Fourier transform

  XTO(:) = cmplx(v3DTO(:,ism+is-1)/NjkTO,0._dp)

  call fftwnd_f77_one(forw,XTO,YTO)

  X(:) = 0
  do m = -njTO/2, +njTO/2
  do n = -nkTO/2, +nkTO/2

    if (m>=0) then
      itTO = m + 1
    else
      itTO = njTO - abs(m) + 1
    endif
    if (n>=0) then
      ipTO = n + 1
    else
      ipTO = nkTO - abs(n) + 1
    endif
    itpTO = (ipTO-1)*njTO + itTO  ! itp_it_ip_lkpTO(itTO,ipTO)

    if (m>=0) then
      it = m + 1
    else
      it = nj - abs(m) + 1
    endif
    if (n>=0) then
      ip = n + 1
    else
      ip = nk - abs(n) + 1
    endif
    itp = (ip-1)*nj + it  ! itp_it_ip_lkp(it,ip)

    X(itp) = YTO(itpTO)

  enddo
  enddo

! 2. Back Fourier transform on fine angular mesh

  ! Complex conjugate for the back Fourier transform
  call fftwnd_f77_one(backw,X,Y)

  v3D(:,is) = real(Y(:))
! Y should be 0 for a real value

 enddo  ! is


 deallocate(XTO)
 deallocate(YTO)
 deallocate(X)
 deallocate(Y)

end subroutine
