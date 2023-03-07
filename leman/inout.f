!  Waves  3D  
!
! Response calculation, results output
!
!
 module inout
 
 use math
 use matrix
 use fields

 
 implicit none


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 contains
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!-------------------------------------------
 subroutine SaveResponse
 integer, parameter :: fn = 2
 integer i,k, mn, j

  open(fn, file=Trim(OUTPUTPATH)//'A.resp')  ! response

! --------- Eq. parameters ----------
  write (fn, *) 'Ra. Plasma radius'
  write (fn, '(1x, e17.10)') Ra
  write (fn, *) 'Rw.  Wall radius'
  write (fn, '(1x, e17.10)') Rw
  write (fn, *) 'Mie. Ion/electron mass ratio'
  write (fn, '(1x, e17.10)') ifr1
  write (fn, *) 'n0_SI.  Plasma density, 1/m^3'
  write (fn, '(1x, e17.10)') n0_SI
  write (fn, *) 'B0.  z-magnetic field on the axis, T'
  write (fn, '(1x, e17.10)') B0
  write (fn, *) 'OmCi_SI. ICF [1/rad]'
  write (fn, '(1x, e17.10)') OmCI1_SI
  write (fn, *) 'Imaginary part in the frequency'
  write (fn, '(1x, e17.10)') om_im
  write (fn, *) 'cA_SI. [m/s]'
  write (fn, '(1x, e17.10)') cA_SI
  write (fn, *) 'm_ant. Antenna: poloidal mode'
  write (fn, *) m_ant
  write (fn, *) 'n_ant. Antenna: toroidal mode'
  write (fn, *) n_ant
  write (fn, *) 'nperiod. Number of toroidal periods'
  write (fn, *) nperiod
  write (fn, *) 'Major radius'
  write (fn, '(1x, e17.10)') MajRadius
  write (fn, *) 'Fourier harmonics (mAFTot, nAFTot), NelR'
  write (fn, *) nwidth*mnAFTot, NelR
  do j = 1, nwidth
    do mn = 1, mnAFTot
      write (fn, *) mAF(mn), nAF(mn,1)+(j-1)*nperiod
    enddo
  enddo

  write(fn, '(1x, 1i5)') OmSteps
  write(fn, '(1x, 2e17.10)') Om_SI1, Om_SI2

  do i = 0, OmSteps
    write (fn, '(1x, e17.10)') Response(i)
  enddo

  close(fn)

 end subroutine
!-------------------------------------------
 subroutine SaveResults
 integer, parameter :: fn = 2
 integer mn,k, eR,gR, is, a,r, itp, itpr, it,ip,j
 complex(dp) ksi(4), dksi(4), dA(4)


  open(fn, file=Trim(OUTPUTPATH)//'A.3D')  ! solution vector XX

! --------- Eq. parameters ----------
  write (fn, *) 'Ra. Plasma radius'
  write (fn, '(1x, e17.10)') Ra
  write (fn, *) 'Rw.  Wall radius'
  write (fn, '(1x, e17.10)') Rw
  write (fn, *) 'Mie. Ion/electron mass ratio'
  write (fn, '(1x, e17.10)') ifr1
  write (fn, *) 'n0_SI.  Plasma density, 1/m^3'
  write (fn, '(1x, e17.10)') n0_SI
  write (fn, *) 'B0.  z-magnetic field on the axis, T'
  write (fn, '(1x, e17.10)') B0
  write (fn, *) 'Om_SI.  Excitation frequency, 1/rad'
  write (fn, '(1x, e17.10)') real(Om_SI)
  write (fn, *) 'Imaginary part in the frequency'
  write (fn, '(1x, e17.10)') om_im
  write (fn, *) 'OmCi_SI. ICF [1/rad]'
  write (fn, '(1x, e17.10)') OmCI1_SI
  write (fn, *) 'cA_SI. [m/s]'
  write (fn, '(1x, e17.10)') cA_SI
  write (fn, *) 'm_ant. Antenna: poloidal mode'
  write (fn, *) m_ant
  write (fn, *) 'n_ant. Antenna: toroidal mode'
  write (fn, *) n_ant
  write (fn, *) 'nperiod. Number of toroidal periods'
  write (fn, *) nperiod
  write (fn, *) 'Major radius'
  write (fn, '(1x, e17.10)') MajRadius

  
! --------- Radial grids ----------
  write (fn, *) 's. Gauss grid.'
  write (fn, *) NelR, NGaussR
  do eR = 1, NelR
  do gR = 1, NGaussR
    write (fn, '(1x, e17.10)') s_gauss(eR,gR)
  enddo
  enddo

  write(fn, *) 's. TERPSICHORE grid.'
  write (fn, *) niT
  do is = 1, niT
    write (fn, '(1x, e17.10)') s_T(is)
  enddo

  write(fn, *) 's. Finite Elements grid.'
  write(fn, *) Nmesh
  do eR = 1, Nmesh
    write (fn, '(1x, e17.10)') s_nod(eR)
  enddo


! --------- Eq. quantities ----------
  write(fn, *) 'Density profile. Gauss grid.'
  do eR = 1, NelR
  do gR = 1, NGaussR
    write (fn, '(1x, e17.10)') fnp(eR,gR)
  enddo
  enddo

  write(fn, *) 'Iota profile + inverse aspect ratio (on the TERPSICHORE grid).'
  do is = 1, niT
    write (fn, '(1x, 2e17.10)') iota(is), iaspct(is)
  enddo


! --------- Perturbations -----------
  write (fn, *) 'P&T fourier modes m,n.'
  write(fn, '(1x, 1i5)') nwidth*mnAFTot
  do j = 1, nwidth
    do mn = 1, mnAFTot
      write (fn, *) mAF(mn), nAF(mn,1)+(j-1)*nperiod
    enddo
  enddo


! --------- Potentials ----------
  write (fn, *) 'Potentials. An, Ab, A||, phi (re + im)'
  do mn = 1, mnAFTot
    do eR = 1, Nmesh
      a = 1
      r = LMat(a,eR)
      write (fn, '(1x, 8e17.10)') real(XXc(r,mn,1)*sf_v_A(1,eR)), aimag(XXc(r,mn,1)*sf_v_A(1,eR)), &
                                  real(XXc(r,mn,2)*sf_v_A(2,eR)), aimag(XXc(r,mn,2)*sf_v_A(2,eR)), &
                                  real(XXc(r,mn,3)*sf_v_A(3,eR)), aimag(XXc(r,mn,3)*sf_v_A(3,eR)), &
                                  real(XXc(r,mn,4)*sf_v_A(4,eR)), aimag(XXc(r,mn,4)*sf_v_A(4,eR))
    enddo ! eR
  enddo ! mn


  write (fn, *) 'Radial derivatives of potentials. An, Ab, A||, phi (re + im)'
  do mn = 1, mnAFTot
    do eR = 1, Nmesh
      a = 1;   r = LMat(a,eR)
      ksi(:)  = XXc(r,mn,:)
      a = 2;   r = LMat(a,eR)
      dksi(:) = XXc(r,mn,:)
    
      dA(:) = dksi(:)*sf_v_A(:,eR) + ksi(:)*dsf_v_A(:,eR)
    
      write (fn, '(1x, 8e17.10)') real(dA(1)), aimag(dA(1)), &
                                  real(dA(2)), aimag(dA(2)), &
                                  real(dA(3)), aimag(dA(3)), &
                                  real(dA(4)), aimag(dA(4))
    enddo ! eR
  enddo ! mn


  write(fn, *) 'Alfven resonance freq., cylinder, TERPSICHORE grid.'
  do is = 1, niT
    do mn = 1,mnAFTot
      write (fn, '(1x, e17.10)') om_a(is,mn)
    enddo
  enddo

! -------- Power fluxes and fields ----------
  write (fn, *) 'Half-points radial mesh'
  do eR = 1, NelR
    write (fn, '(1x, e15.8)') s_hp(eR)
  enddo

  write (fn, *) 'Poynting flux, half-mesh'
  do eR = 1, NelR
    write (fn, '(1x, 2e15.8)') real(PpoyntSum(eR)), aimag(PpoyntSum(eR))
  enddo
  write (fn, *) 'Ppla, FE mesh'
  do eR = 1, Nmesh
    write (fn, '(1x, 2e15.8)') real(PplaSum(1,eR)), aimag(PplaSum(1,eR))
  enddo
  write (fn, *) 'Pant, FE mesh'
  do eR = 1, Nmesh
    write (fn, '(1x, 2e15.8)') real(PantSum(eR)), aimag(PantSum(eR))
  enddo


  write (fn, *) 'TERPSICHORE angular grid size'
  write (fn, *) njT, nkT

  ! Divergence of A in real space
  if (OutputFlags(1)) then
    write (fn, *) 'Div(A), real space, 3D'
    gR = NGaussR/2+1
    do eR = 1, NelR
      do itp = 1, Njk
        itpr = (gR-1)*Njk + itp
        write (fn, '(2e11.3)') real(divA3D(itpr,eR)), imag(divA3D(itpr,eR))
      enddo
    enddo
  endif

  ! Electric and magnetic fields in real space
  if (OutputFlags(2)) then
    write (fn, *) 'E&B in 3D'
    do eR = 1, NelR
      do itp = 1, Njk
        write (fn, '(12e12.4)') real(E3D_hp(1,eR,itp)), imag(E3D_hp(1,eR,itp)), &
                                real(E3D_hp(2,eR,itp)), imag(E3D_hp(2,eR,itp)), &
                                real(E3D_hp(3,eR,itp)), imag(E3D_hp(3,eR,itp)), &
                                real(B3D_hp(1,eR,itp)), imag(B3D_hp(1,eR,itp)), &
                                real(B3D_hp(2,eR,itp)), imag(B3D_hp(2,eR,itp)), &
                                real(B3D_hp(3,eR,itp)), imag(B3D_hp(3,eR,itp))
      enddo
    enddo
  endif


  close(fn)


 end subroutine
!-------------------------------------------
 subroutine SaveRawData
 integer, parameter :: fn = 2
 integer, dimension(:), allocatable :: tmpAF
 real(dp) Om_SI_re
 integer j, mn

  allocate(tmpAF(nwidth*mnAFTot))
 
  open(fn, file=Trim(OUTPUTPATH)//'RawA.dat',form='unformatted')  ! solution vector XX

  ! Parameters
  Om_SI_re = real(Om_SI)
  write(fn) ifr1, n0_SI, B0, Om_SI_re, om_im, OmCI1_SI, cA_SI, MajRadius
  write(fn) m_ant, n_ant, nperiod

  ! Fourier mode numbers
  write(fn) Nmesh, nwidth*mnAFTot

  do j = 1, nwidth
    tmpAF((j-1)*mnAFTot+1:j*mnAFTot) = mAF(1:mnAFTot)
  enddo
  write(fn) tmpAF

  do j = 1, nwidth
    do mn = 1, mnAFTot
      tmpAF((j-1)*mnAFTot+mn) = nAF(mn,1) + (j-1)*nperiod
    enddo
  enddo
  write(fn) tmpAF


  ! Potentials
  write(fn) s_nod
  write(fn) real(XXc)
  write(fn) imag(XXc)
  write(fn) real(vA_fe(:, 1:mnAFTot, :))
  write(fn) imag(vA_fe(:, 1:mnAFTot, :))

  close(fn)

  deallocate(tmpAF)

 end subroutine
!-------------------------------------------
 subroutine ShowEqTable
  integer ieq, m,n
  integer maxm, maxn

  real(dp), dimension(:,:), allocatable :: vt
  logical, dimension(:,:), allocatable :: vtexist

  real(dp) maxAmp
  real(dp), parameter :: minEqAmp = 1.0e-3_dp

  character(250) s, stmp
  
  maxm = abs(mEq(1));  maxn = abs(nEq(1));
  do ieq = 1, mnEqTot
    if (abs(mEq(ieq))>maxm) maxm = abs(mEq(ieq))
    if (abs(nEq(ieq))>maxn) maxn = abs(nEq(ieq))
  enddo
  maxm = maxm + 1
  maxn = maxn/nperiod + 1
  
  print *, 'Maxm:', maxm
  print *, 'Maxn:', maxn
!  print *, nEq(:)

  allocate(vt(-maxm:maxm, -maxn:maxn))
  allocate(vtexist(-maxm:maxm, -maxn:maxn))

  vt = 0; vtexist = .FALSE.
  do ieq = 1, mnEqTot
    if ((mEq(ieq)>=-maxm).and.(mEq(ieq)<=maxm).and. &
        (nEq(ieq)/nperiod>=-maxn).and.(nEq(ieq)/nperiod<=maxn)) then
      vt(mEq(ieq),nEq(ieq)/nperiod)      = sqrt(maxFEqAmpMax(ieq))
      vtexist(mEq(ieq),nEq(ieq)/nperiod) = .TRUE.
!    if ((mEq(ieq)>=-maxm).and.(mEq(ieq)<=maxm).and. &
!        (nEq(ieq)/nperiod>=-maxn).and.(nEq(ieq)/nperiod<=maxn)) then
!      vt(mEq(ieq),nEq(ieq)/nperiod)      = sqrt(maxFEqAmp(ieq))
!      vtexist(mEq(ieq),nEq(ieq)/nperiod) = .TRUE.
    else
      print *, 'Insufficient storage for the equilibrium Fourier mode table'
      stop
    endif
  enddo

  print *, ' Equilibrium Fourier mode amplitude table'
  s = 'M'
  do m = -maxm, maxm
    write(stmp,'(i4)') m
    if (s(len_trim(s):len_trim(s))=='M') then
      s = trim(s) // stmp
    else
      s = trim(s) // '    ' // stmp
    endif
  enddo
  s = trim(s) // '       N'
  print *, s

  do n = -maxn, maxn
    s = ''
    do m = -maxm, maxm

      if (vtexist(m,n)) then
        if ((abs(m) > mEqMax).or.(abs(n) > nEqMax)) then
!          write(stmp,'(a5)') '+'
          write(stmp,'(e8.1)') vt(m,n)
        else
          write(stmp,'(e8.1)') vt(m,n)
        endif
      else
        write(stmp,'(a5)') '-'
      endif


      if ((s(len_trim(s):len_trim(s))=='-').or.(s(len_trim(s):len_trim(s))=='+')) then
        s = trim(s) // '   ' // stmp
      else
        s = trim(s) // stmp
      endif
    enddo
    write(stmp,'(i8)') n*nperiod
    s = trim(s) // stmp
    print *, trim(s)
    print *, " "
  enddo

  maxAmp = sqrt(maxFEqAmpMax(1))
  do ieq = 2, mnEqTot
    if (sqrt(maxFEqAmp(ieq))>maxAmp) maxAmp = sqrt(maxFEqAmp(ieq))
  enddo
  maxAmp = max(maxAmp,1.0e-30_dp)
  print *, maxAmp

  print *, ' Equilibrium Fourier mode amplitude table (trace)'
  s = 'M'
  do m = -maxm+3, maxm, 4
    write(stmp,'(i4)') m
    s = trim(s) // stmp
  enddo
  s = trim(s) // '   N'
  print *, trim(s)

  do n = -maxn, maxn
    s = ''
    do m = -maxm, maxm

      if (vtexist(m,n)) then

!        if ((abs(m) > mEqMax).or.(abs(n) > nEqMax)) then
!          write(stmp,'(a1)') '+'
!        else

          if (vt(m,n)/maxAmp > minEqAmp) then
            write(stmp,'(a1)') '*'
          else
            write(stmp,'(a1)') '0'
          endif

!        endif

      else
        write(stmp,'(a1)') '-'
      endif

      s = trim(s) // stmp

    enddo
    write(stmp,'(i6)') n*nperiod
    s = trim(s) // stmp
    print *, trim(s)
  enddo

  deallocate(vt)
  deallocate(vtexist)

  print *, 'Legend:  ''-'' - eq. mode is not included in the sum'
  print *, '         ''+'' - eq. mode is included, but forced to zero (m/nEqMax parameter)'
  print '(a,f8.1,a)', '         ''*'' - eq. mode amplitude is greater than ',minEqAmp,' maxAmp'
  print '(a,f8.1,a)', '         ''0'' - eq. mode amplitude is less than ',minEqAmp,' maxAmp'

 end subroutine
!-------------------------------------------
 subroutine ShowAFTable
  integer m,n,mn,mm,nn, kk,r
  integer minm,maxm, minn,maxn

  real(dp), dimension(:,:), allocatable :: vt
  logical, dimension(:,:), allocatable  :: vtexist
  real(dp), dimension(:), allocatable   :: mnAFAmp

  real(dp) maxAmp

  character(1000) s, stmp
  integer, parameter :: fn = 2
  
  minm = minval(mAF(:));  maxm = maxval(mAF(:));
  minn = minval(nAF(:,1));  maxn = maxval(nAF(:,1))+(nwidth-1)*nperiod;

  mm = maxm - minm + 1
  nn = (maxn-minn)/nperiod + 1

  allocate(vt(mm,nn))
  allocate(vtexist(mm,nn))
  allocate(mnAFAmp(mnAFTot))


  maxAmp = 1.0e-30_dp
  mnAFAmp(:) = 0
  do kk = 1,3
  do r = 1,Nmesh
  do mn = 1, mnAFTot
    mnAFAmp(mn) = max(mnAFAmp(mn),abs(vA_fe(kk, mn, r)))
    maxAmp = max(maxAmp,abs(vA_fe(kk, mn, r)))
  enddo
  enddo
  enddo


  vt = 0; vtexist = .FALSE.
  do mn = 1, mnAFTot
   m = mAF(mn) - minm + 1
   n = (nAF(mn,1) - minn)/nperiod + 1

   vt(m,n)      = mnAFAmp(mn)
   vtexist(m,n) = .TRUE.
  enddo


  open(fn, file=Trim(OUTPUTPATH)//'AFdisplay.dat',form='unformatted')
  write(fn) mm,nn
  write(fn) vt
  close(fn)


  write (16, *) ' Perturbed "A" Fourier mode amplitude table'
  s = 'M'
  do m = 1, mm
    write(stmp,'(i4)') minm + m - 1
    if (s(len_trim(s):len_trim(s))=='M') then
      s = trim(s) // stmp
    else
      s = trim(s) // '    ' // stmp
    endif
  enddo
  s = trim(s) // '       N'
  write (16, *) trim(s)

  do n = 1, nn
    s = ''
    do m = 1, mm

      if (vtexist(m,n)) then
        write(stmp,'(e8.1)') vt(m,n)
      else
        write(stmp,'(a5)') '-'
      endif


      if ((s(len_trim(s):len_trim(s))=='-').or.(s(len_trim(s):len_trim(s))=='+')) then
        s = trim(s) // '   ' // stmp
      else
        s = trim(s) // stmp
      endif
    enddo
    write(stmp,'(i8)') minn + (n-1)*nperiod
    s = trim(s) // stmp
    write (16, *) trim(s)
  enddo

  write (16, *) maxAmp

  write (16, *) ' Perturbed "A" Fourier mode amplitude table (trace)'
  s = 'M'
  m = 1
  write(stmp,'(i3)') minm + m - 1
  s = trim(s) // stmp
  do m = 1, mm
    write(stmp,'(i4)') minm + m - 1
    s = trim(s) // stmp
  enddo
  s = trim(s) // '   N'
  write (16, *) trim(s)

  do n = 1, nn
    s = ''
    do m = 1, mm

      if (vtexist(m,n)) then

          if (vt(m,n)/maxAmp > 4.e-3) then
            write(stmp,'(a4)') '*'
          else
            write(stmp,'(a4)') '0'
          endif

      else
        write(stmp,'(a4)') '-'
      endif

      s = trim(s) // stmp

    enddo
    write(stmp,'(i4)') minn + (n-1)*nperiod
    s = trim(s) // stmp
    write (16, *) trim(s)
  enddo

  deallocate(vt)
  deallocate(vtexist)
  deallocate(mnAFAmp)

  write (16, *) 'Legend:  ''-'' - eq. mode is not included in the sum'
  write (16, *) '         ''*'' - eq. mode amplitude is greater than 1.e-4*maxAmp'
  write (16, *) '         ''0'' - eq. mode amplitude is less than 1.e-4*maxAmp'

 end subroutine
!-------------------------------------------

!-------------------------------------------

 end module inout
