!  E/M waves propagation in 3D plasmas.  Main module.
!
!  Cold plasma
!  Potential formulation of the wave equation
!  Normalization:  All values are in SI
!
! For the comparison with LION the following normalization is used:
!                  freq    ->  1/(om.Larmor.ion*a/Ca)
!                  length  ->  1/a
!                  time    ->  1/(a/Ca),  Ca - Alfven velosity
!                  B(r)  = B0*Bp(r),  B0 = B(r=0)
!                  n(r)  = n0*np(r),  n0 = n(r=0)
!
!
!  Version 3.0: 3D, general geometry, cold plasma, Boozer coordinates, values in SI
!  Cubic FE (r)
!
!  theta,phi decomposition: exp(i*m*theta), exp(i*n*phi)
!
!
!
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!==========================================================================
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 program WAVES3D

 use math
 use matrix
 use inout
 use fields
 use unicity
 use transition

 implicit none

 integer, parameter :: FFTW_FORWARD  = -1
 integer, parameter :: FFTW_BACKWARD = +1
 integer, parameter :: FFTW_ESTIMATE = 64
 
 integer totmem, e, i
 double precision mem
 integer, parameter :: fid = 2

 call MPI_INIT(ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

 call fftw2d_f77_create_plan(forw,njTO,nkTO,FFTW_FORWARD,FFTW_ESTIMATE)
 call fftw2d_f77_create_plan(backw,njT,nkT,FFTW_BACKWARD,FFTW_ESTIMATE)
 call fftw2d_f77_create_plan(trans,njT,nkT,FFTW_FORWARD,FFTW_ESTIMATE)
 
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


   call ReadNameList
   call Init_Para
   call Init_LM
   call Initialize_Gauss
   ! Reading equilibrium quantities from the file generated with TERPSICHORE (transformed to binary)
   call ReadGrid
   call ReadEquilibriumData
   call ReadAFTable
   call InitializeVariables
 
 !  call ShowEqTable


   call CylinderAlvfenResonance

  if (me.eq.0) then 

  print *, '* Ver. 3.0:  General geometry (metric from TERPSICHORE), 3D'
  print *, '* Input files:  parameters  - from standard input'
  print *, '*               metric      - fort.37.DAT.N'
  write (*,'(a17,a30)') ' * Equilibrium:  ', eqname
  print *, ' '
  print *, 'Equilibrium parameters:'
  write (*,'(a18, 1e11.3, a)')    'n0     = ',n0_SI, ' m^(-3)'
  write (*,'(a18, 1e11.3, a)')    'B0     = ',B0   , ' T'
  write (*,'(a18, 1e11.3, a)')    'R0     = ',MajRadius,' m'
  write (*,'(a18, 1e11.3     )')  'ifr1   = ',ifr1
  write (*,'(a18, 1e11.3     )')  'ifr2   = ',ifr2
  write (*,'(a18, 1e11.3     )')  'C/Ca   = ',CCa
  write (*,'(a18, 1e11.3, a)')    'Om(ci1)= ',OmCI1_SI, ' [rad/s]'
  write (*,'(a18, 1e11.3, a)')    'Om(ci2)= ',OmCI2_SI, ' [rad/s]'
  write (*,'(a18, 1e11.3, a)')    'Om(ci)*a/Ca = ',OmCI1_SI*a_plasma/cA_SI, ' [0]'
  print *, 'Antenna:'
  write (*,'(a18, i5)')    'Toroidal mode = ',n_ant
  write (*,'(a18, i5)')    'Poloidal mode = ',m_ant
  write (*,'(a18, 1e11.3, a,1e9.3,a)')  'Omega  = ',Om_SI1,'..',Om_SI2, ' rad/s'
  write (*,'(a18, 1e11.3)')  'Imag. part  = ',om_im
 
!     print *, 'Fourier decomposition: m(theta)  =',mAF(:),  '   n(z)  =',nAF(:)
!     print *, '                       mp(theta) =',mpAF(:), '   np(z) =',npAF(:)
 
  print *, 'Integration: '
  print *, '  ', NelR, ' radial elements,     ', &
           & NGaussR, 'radial gauss points'
  print *, '  ', njT, ' poloidal elements   '
  print *, '  ', nkT, ' toroidal elements   '
  print *, '  MSize2=',MSize2
  print *, ' '
  print *, ' '

  end if


! Loop for each frequency

  do qq = 0, OmSteps
    if (OmSteps > 0) then
      Om_SI = Om_SI1 + (Om_SI2-Om_SI1)/OmSteps * qq
    elseif (Om_SI1 > 0.) then
      Om_SI = Om_SI1
    else
      Om_SI = iZh * 1.6e-19 * Bc / (iAh * mp_SI)
    endif

   call InitialWaveVector
   
   do qiter = 1, NbIter

    if (me.eq.0) then     
      write (*,'(a10, 1e9.3, a, 1e9.3, a, i3,a,i3,a)') 'Omega = ', Om_SI, &
      ' rad/s;  Omega/OmCI = ', Om_SI/OmCI1_SI,';   Step ',qq,'/',OmSteps,';'
      write (*,'(a12, i3,a,i3)') 'Iteration ', qiter,'/',NbIter
    end if

    ! Matrix contruction + Solver

    call SolveStart

  do nindex = m_npar, nwidth, np_npar

    call Init_RHS
    call RHS_BConditions


    KKD = 0._dp
    i = 1
    do e = gmin, gmax, incr
      KKA = 0._dp; KKE = 0._dp; KK_tmp = 0._dp
      if (npsep<=2) KKE(:,:) = KKD(:,:)
      KKD = 0._dp;
      call Init_Matrices(e,i)
      if (e==1) call Axis_BConditions
      if (e==NelR) call Border_BConditions
      call GaussElimination(e,i)
      if (nprocs==1.or.sgme==mod(min((NelR-e)/2,(e-1)/2),same)) i = i+1
    enddo

    call SolveEq(i)

  enddo

!    if (me.eq.0) then
!    open(fid,file='memory.txt')
!    write (fid, *) mnAFTot
!    write (fid, '(1x, 1e19.10)') mem()
!    close(fid)
!    end if

    call SolveToTreatment

!    call SaveAntenna
    call EvaluateFields

    if (me.eq.0) call SaveRawData

    if (me.eq.0.and.qiter==NbIter) then
      Response(qq) = abs(imag(PantSum(Nmesh)))

      if (qq==0) then
        call SaveResults
        call SaveRawData
      endif
    end if

   call FinishTreatment

    enddo
 
   enddo ! qq

   if ((OmSteps > 0).and.(me.eq.0)) call SaveResponse


! Additional loop for Maximum response frequency calculation
   
   if ((MaxResponseFreq == 1).and.(OmSteps > 0)) then

    call Max_Resp(qq)

    Om_SI = Om_SI1 + (Om_SI2-Om_SI1)/OmSteps * qq
!    Om_SI = Om_SI*(1._dp + iim*om_im)
    print *, ' '
    write (*,'(a, 1e9.3, a, 1e9.3)') 'Maximum response at: Omega = ', &
      real(Om_SI), ' rad/s;  Omega/OmCI = ', real(Om_SI)/OmCI1_SI

    ! Matrix contruction + Solver

    call Init_RHS
    call RHS_BConditions

    call SolveStart

    i = 1
    do e = gmin, gmax, incr
      KKA = 0._dp; KKE = 0._dp; KK_tmp = 0._dp 
      if (nprocs<=2) KKE(:,:) = KKD(:,:)
      KKD = 0._dp;
      call Init_Matrices(e,i)
      if (e==1) call Axis_BConditions
      if (e==NelR) call Border_BConditions
      call GaussElimination(e,i)
      if (nprocs==1.or.sgme==mod(min((NelR-e)/2,(e-1)/2),same)) i = i+1
    enddo

    call SolveEq(i)

    call SolveToTreatment

    call EvaluateFields
    call walltime
    if (me.eq.0) then
      call SaveResults
    end if

   call FinishTreatment

   endif


 
  if (me.eq.0) then

    call ShowEqTable

    totmem = 16*(size(KKA)+size(KKB)+size(KKC)+size(CA)+size(CF) &
           & +size(Acoeff))+8*(size(mTU))

    print *, 'Memory usage:'
    write (*,'(a18, i12)')    'KK     = ',size(KKA)+size(KKB)+size(KKC)
    write (*,'(a18, i12)')    'CA+CF  = ',size(CA)+size(CF)
    write (*,'(a18, i12)')    'Acoeff = ',size(Acoeff)
  if (FourierCheck) &
    write (*,'(a18, i12)')    'EqExp  = ',size(EqExp)
    write (*,'(a18, i12)')    'mEQUIL = ',size(mTU)*9
    print *, ' '
    write (*,'(a18, f12.4, a3)')    'Total  = ', totmem/1.0e6, " MB"

    call ShowAFTable

  end if
  call walltime

  call DeallocateAll
 
!-------------------------------------------

  call fftwnd_f77_destroy_plan(forw)
  call fftwnd_f77_destroy_plan(backw)
  call fftwnd_f77_destroy_plan(trans)


  call MPI_FINALIZE(ierr)

 end program
