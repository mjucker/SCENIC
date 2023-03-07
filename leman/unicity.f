!
!
!  E/M waves propagation in 3D plasmas  -- Appendix.
!
!  Check the unicity of the solution on the axis
!
!  Version 3.0: 3D, general geometry, fixed boundary
!
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!==========================================================================
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


 module Unicity

 use math
 use matrix
 use fields
 
 implicit none

 ! dimensions: (:,:) <-> DIM,DIM,NelR,njk

 real(dp),  dimension(:,:,:,:), allocatable :: TU_hp, dRR_hp
 complex(dp), dimension(:,:,:), allocatable :: ARZ3D_hp, ERZ3D_hp, &
                                               BRZ3D_hp


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 contains
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!-------------------------------------------
subroutine Init_Equilibruim_Unicity(eR) 
 integer ism,isp, eR,gR, it,ip,itp
 integer  it_m,it_p, ip_m,ip_p
 real(dp) s, theta,phi,  R1,R1m,R1p, Z1,Z1m,Z1p
 

  gR = NGaussR/2+1

!  ! Cylinder
!  s = s_gauss(eR,gR)
!  do itp = 1, Njk
!    theta = th_itp(itp)
!    cR_hp(eR,itp) = sqrt(s)*cos(theta)
!    cZ_hp(eR,itp) = sqrt(s)*sin(theta)
!    phiv_hp(eR,itp) = ph_itp(itp)
!  enddo



    ! ism,isp - radial surfaces on the TERPSICHORE grid, neighbour to the s_gauss(eR,gR) surface
    ism = min(niT-1, max(1, int(s_gauss(eR,gR)*niT + 0.5_dp)))
    isp = ism + 1

  do itp = 1, Njk  ! 1..N <-- used for poloidal derivative calculation

    s     = s_gauss(eR,gR)
    ip = (itp-1)/njT + 1

    R1m = (trx(itp,ism)-MAxisR(ip))/sqrt(s_T(ism))
    R1p = (trx(itp,isp)-MAxisR(ip))/sqrt(s_T(isp))
    R1  = R1m + (s-s_T(ism))*(R1p-R1m)/(s_T(isp)-s_T(ism))

    Z1m = (trz(itp,ism)-MAxisZ(ip))/sqrt(s_T(ism))
    Z1p = (trz(itp,isp)-MAxisZ(ip))/sqrt(s_T(isp))
    Z1  = Z1m + (s-s_T(ism))*(Z1p-Z1m)/(s_T(isp)-s_T(ism))

    ! dR_ds
    dRR_hp(1,1,eR,itp) = 0.5_dp/sqrt(s)*R1 + sqrt(s)*(R1p-R1m)/(s_T(isp)-s_T(ism))
    ! dphi_ds
    dRR_hp(2,1,eR,itp) = (tphiv(itp,isp)-tphiv(itp,ism))/(s_T(isp)-s_T(ism))
    ! dZ_ds
    dRR_hp(3,1,eR,itp) = 0.5_dp/sqrt(s)*Z1 + sqrt(s)*(Z1p-Z1m)/(s_T(isp)-s_T(ism))

  enddo ! itp


  do itp = 1, Njk

    it = mod(itp-1,njT) + 1
    ip = (itp-1)/njT + 1

    ! Calculating poloidal derivatives from the radially interpolated values
    it_m = it - 1
    if (it_m==0)     it_m = njT
    it_p = it + 1
    if (it_p==njT+1) it_p = 1


!    R1m = cR_hp(eR,itp_it_ip_lkp(it_m,ip))/cos(theta_T(it-1))
!    R1p = cR_hp(eR,itp_it_ip_lkp(it_p,ip))/cos(theta_T(it+1))
!    R1  = R1m + (theta_T(it)-theta_T(it-1))*(R1p-R1m)/(theta_T(it+1)-theta_T(it-1))
!
!    Z1m = cZ_hp(eR,itp_it_ip_lkp(it_m,ip))/sin(theta_T(it-1))
!    Z1p = cZ_hp(eR,itp_it_ip_lkp(it_p,ip))/sin(theta_T(it+1))
!    Z1  = Z1m + (theta_T(it)-theta_T(it-1))*(Z1p-Z1m)/(theta_T(it+1)-theta_T(it-1))
!
!
!    ! dR_dtheta
!    dRR_hp(1,2,eR,itp) = -sin(theta_T(it))*R1 + &
!         cos(theta_T(it))*(R1p-R1m)/(theta_T(it+1)-theta_T(it-1))
!    ! dphi_dtheta
!    dRR_hp(2,2,eR,itp) = 0
!    ! dZ_dtheta
!    dRR_hp(3,2,eR,itp) = cos(theta_T(it))*Z1 + &
!         sin(theta_T(it))*(Z1p-Z1m)/(theta_T(it+1)-theta_T(it-1))


    ! dR_dtheta
    dRR_hp(1,2,eR,itp) = (cR_hp(eR,itp_it_ip_lkp(it_p,ip)) - &
                          cR_hp(eR,itp_it_ip_lkp(it_m,ip)))/(theta_T(it+1)-theta_T(it-1))
    ! dphi_dtheta
    dRR_hp(2,2,eR,itp) = (phiv_hp(eR,itp_it_ip_lkp(it_p,ip)) - &
                          phiv_hp(eR,itp_it_ip_lkp(it_m,ip)))/(theta_T(it+1)-theta_T(it-1))
    ! dZ_dtheta
    dRR_hp(3,2,eR,itp) = (cZ_hp(eR,itp_it_ip_lkp(it_p,ip)) - &
                          cZ_hp(eR,itp_it_ip_lkp(it_m,ip)))/(theta_T(it+1)-theta_T(it-1))



    ! Calculating toroidal derivatives from the newly interpolated values
    ip_m = ip - 1
    if (ip_m==0)     ip_m = nkT
    ip_p = ip + 1
    if (ip_p==nkT+1) ip_p = 1


    ! dR_dphi
    dRR_hp(1,3,eR,itp) = (cR_hp(eR,itp_it_ip_lkp(it,ip_p)) - &
                          cR_hp(eR,itp_it_ip_lkp(it,ip_m)))/(phi_T(ip+1)-phi_T(ip-1))
    ! dphi_dphi
    dRR_hp(2,3,eR,itp) = (phiv_hp(eR,itp_it_ip_lkp(it,ip_p)) - &
                          phiv_hp(eR,itp_it_ip_lkp(it,ip_m)))/(phi_T(ip+1)-phi_T(ip-1))
    ! dZ_dphi
    dRR_hp(3,3,eR,itp) = (cZ_hp(eR,itp_it_ip_lkp(it,ip_p)) - &
                          cZ_hp(eR,itp_it_ip_lkp(it,ip_m)))/(phi_T(ip+1)-phi_T(ip-1))

  enddo ! itp



!  ! Cylinder
!  s = s_gauss(eR,gR)
!  do itp = 1, Njk
!    theta = th_itp(itp)
!
!    ! dR_ds
!    dRR_hp(1,1,eR,itp) = 0.5/sqrt(s)*cos(theta)
!    ! dphi_ds
!    dRR_hp(2,1,eR,itp) = 0
!    ! dZ_ds
!    dRR_hp(3,1,eR,itp) = 0.5/sqrt(s)*sin(theta)
!
!    ! dR_dtheta
!    dRR_hp(1,2,eR,itp) = -sqrt(s)*sin(theta)
!    ! dphi_dtheta
!    dRR_hp(2,2,eR,itp) = 0
!    ! dZ_dtheta
!    dRR_hp(3,2,eR,itp) = sqrt(s)*cos(theta)
!
!    ! dR_dtheta
!    dRR_hp(1,3,eR,itp) = 0
!    ! dphi_dtheta
!    dRR_hp(2,3,eR,itp) = 1
!    ! dZ_dtheta
!    dRR_hp(3,3,eR,itp) = 0
!
!  enddo ! itp



end subroutine
!-------------------------------------------
subroutine InitTransformMatrices
 integer eR,gR, itp
 real(dp)  M3(3,3), iM3(3,3)
 real(dp)  x1,x2,x3, y(3), a(3)

 
  gR = NGaussR/2+1
  TU_hp = 0
  do eR = 1, NelR

    ! Interpolate the equilibrium (TU,TL,B,etc) on the gauss radial
    ! grid for the radial element "e" (from TERPSICHORE grid)
    if (CylMetric) then
      call Init_Cyl_Equilibruim(eR)
    else
      call Init_Equilibruim(eR)
    endif
    call Init_Equilibruim_Unicity(eR)

    TU_hp(1,1,eR,:) = mTU(1,1,gR,:) * seTU(1,eR,gR)
    TU_hp(2,2,eR,:) = mTU(2,2,gR,:) * seTU(2,eR,gR)
    TU_hp(3,3,eR,:) = mTU(3,3,gR,:)

!    TU_hp(1,1,eR,:) = 2.0_dp*sqrt(s_gauss(eR,NGaussR/2+1))
!    TU_hp(2,2,eR,:) = 1.0_dp/sqrt(s_gauss(eR,NGaussR/2+1))
!    TU_hp(3,3,eR,:) =-1.0_dp

  enddo

end subroutine
!-------------------------------------------
subroutine InitARZ
 integer eR, it,ip,itp, ii,kk,jj
 real(dp) Rot(3,3)
 real(dp) s, theta

  ARZ3D_hp = 0
  ERZ3D_hp = 0
  BRZ3D_hp = 0
  do eR = 1, NelR
  do itp = 1, Njk
    Rot = 0
    do ii = 1, 3
    do jj = 1, 3
    do kk = 1, 3

      Rot(ii,jj) = Rot(ii,jj) + dRR_hp(ii,kk,eR,itp)*TU_hp(kk,jj,eR,itp)

    enddo
    enddo
    enddo

    do ii = 1, 3
    do kk = 1, 3
      ARZ3D_hp(ii,eR,itp) = ARZ3D_hp(ii,eR,itp) + Rot(ii,kk)*A3D_hp(kk,eR,itp) 
      ERZ3D_hp(ii,eR,itp) = ERZ3D_hp(ii,eR,itp) + Rot(ii,kk)*E3D_hp(kk,eR,itp) 
      BRZ3D_hp(ii,eR,itp) = BRZ3D_hp(ii,eR,itp) + Rot(ii,kk)*B3D_hp(kk,eR,itp) 
    enddo
    enddo

  enddo
  enddo

  ! Extrapolation to the axis
  do itp = 1, Njk
  do ii = 1, DIM
    ARZ3D_hp(ii,0,itp) = ARZ3D_hp(ii,1,itp) - s_hp(1)/(s_hp(2)-s_hp(1))*(ARZ3D_hp(ii,2,itp)-ARZ3D_hp(ii,1,itp))
    ERZ3D_hp(ii,0,itp) = ERZ3D_hp(ii,1,itp) - s_hp(1)/(s_hp(2)-s_hp(1))*(ERZ3D_hp(ii,2,itp)-ERZ3D_hp(ii,1,itp))
    BRZ3D_hp(ii,0,itp) = BRZ3D_hp(ii,1,itp) - s_hp(1)/(s_hp(2)-s_hp(1))*(BRZ3D_hp(ii,2,itp)-BRZ3D_hp(ii,1,itp))
  enddo
  enddo


!  ARZ3D_hp = 0
!  do eR = 1, NelR
!  s = s_gauss(eR,NGaussR/2+1)
!  do itp = 1, Njk
!    theta = th_itp(itp)
!
!    ARZ3D_hp(1,eR,itp) = cos(theta)*A3D_hp(1,eR,itp) - sin(theta)*A3D_hp(2,eR,itp)
!    ARZ3D_hp(2,eR,itp) = sin(theta)*A3D_hp(1,eR,itp) + cos(theta)*A3D_hp(2,eR,itp)
!    ARZ3D_hp(3,eR,itp) = A3D_hp(3,eR,itp)
!
!  enddo
!  enddo


end subroutine
!-------------------------------------------
 subroutine SaveARZ
 integer, parameter :: fn = 2
 
  open(fn, file=Trim(OUTPUTPATH)//'RawARZ.dat',form='unformatted')

  write(fn) NelR, njT, nkT
  write(fn) cR_hp
  write(fn) cZ_hp
  write(fn) real(ARZ3D_hp)
  write(fn) imag(ARZ3D_hp)
  write(fn) real(ERZ3D_hp)
  write(fn) imag(ERZ3D_hp)
  write(fn) real(BRZ3D_hp)
  write(fn) imag(BRZ3D_hp)

  close(fn)

 end subroutine
!-------------------------------------------
subroutine CheckUnicity
 integer ii, it, ip, itp, eR
 integer, parameter :: fn = 2
 real(dp) UnA(DIM,0:NelR),  avrgA(DIM)
 real(dp) UnE(DIM,0:NelR),  avrgE(DIM)
 real(dp) UnB(DIM,0:NelR),  avrgB(DIM)
 complex(dp) avrg(3)

 allocate(TU_hp(DIM,DIM,NelR,Njk))
 allocate(dRR_hp(DIM,DIM,NelR,Njk))
 allocate(ARZ3D_hp(DIM,0:NelR,Njk))
 allocate(ERZ3D_hp(DIM,0:NelR,Njk))
 allocate(BRZ3D_hp(DIM,0:NelR,Njk))


  call InitTransformMatrices
  call InitARZ
  call SaveARZ


  ip = 1
  UnA = 0
  UnE = 0
  UnB = 0
  avrgA = 0
  avrgE = 0
  avrgB = 0
  do ii = 1, DIM
  do eR = 0, NelR
    do it = 1, njT
      itp = itp_it_ip_lkp(it,ip)
      avrgA(ii) = avrgA(ii) + ARZ3D_hp(ii,eR,itp)*conjg(ARZ3D_hp(ii,eR,itp))
      avrgE(ii) = avrgE(ii) + ERZ3D_hp(ii,eR,itp)*conjg(ERZ3D_hp(ii,eR,itp))
      avrgB(ii) = avrgB(ii) + BRZ3D_hp(ii,eR,itp)*conjg(BRZ3D_hp(ii,eR,itp))
    enddo
  enddo
  enddo
  avrgA = avrgA/njT
  avrgE = avrgE/njT
  avrgB = avrgB/njT

  do ii = 1, DIM
  do eR = 0, NelR
  
    avrg = 0
    do it = 1, njT
      itp = itp_it_ip_lkp(it,ip)
      avrg(1) = avrg(1) + ARZ3D_hp(ii,eR,itp)
      avrg(2) = avrg(2) + ERZ3D_hp(ii,eR,itp)
      avrg(3) = avrg(3) + BRZ3D_hp(ii,eR,itp)
    enddo
    avrg(:) = avrg(:)/njT
  
    do it = 1, njT
      itp = itp_it_ip_lkp(it,ip)
      UnA(ii,eR) = UnA(ii,eR) + &
        real((ARZ3D_hp(ii,eR,itp)-avrg(1))*conjg(ARZ3D_hp(ii,eR,itp)-avrg(1)))
      UnE(ii,eR) = UnE(ii,eR) + &
        real((ERZ3D_hp(ii,eR,itp)-avrg(2))*conjg(ERZ3D_hp(ii,eR,itp)-avrg(2)))
      UnB(ii,eR) = UnB(ii,eR) + &
        real((BRZ3D_hp(ii,eR,itp)-avrg(3))*conjg(BRZ3D_hp(ii,eR,itp)-avrg(3)))
    enddo
!    if (abs(avrg)<1.e-30) avrg = 1.e-30
    UnA(ii,eR) = UnA(ii,eR)/njT/avrgA(ii)
    UnE(ii,eR) = UnE(ii,eR)/njT/avrgE(ii)
    UnB(ii,eR) = UnB(ii,eR)/njT/avrgB(ii)
  
  enddo
  enddo

  open(fn, file=Trim(OUTPUTPATH)//'UnicityARZ.dat',form='unformatted')

  write(fn) NelR
  write(fn) s_hp
  write(fn) UnA
  write(fn) UnE
  write(fn) UnB

  close(fn)


  open(fn, file='convergence.txt',position='append')

  write(fn,*) 'Unicity:'

  write(fn,'(1x, 3e19.10)') UnA(1,0),UnA(2,0),UnA(3,0)
  write(fn,'(1x, 3e19.10)') UnE(1,0),UnE(2,0),UnE(3,0)
  write(fn,'(1x, 3e19.10)') UnB(1,0),UnB(2,0),UnB(3,0)

  close(fn)

  deallocate(TU_hp)
  deallocate(dRR_hp)
  deallocate(ARZ3D_hp)
  deallocate(ERZ3D_hp)
  deallocate(BRZ3D_hp)

end subroutine
!-------------------------------------------

 end module unicity



























!  ! Radial derivatives -------------------------
!  do eR = 2, NelR-1
!
!   x1 = s_hp(eR-1)
!   x2 = s_hp(eR)
!   x3 = s_hp(eR+1)
!   M3(1,1) = x1*x1;   M3(1,2) = x1;   M3(1,3) = 1._dp
!   M3(2,1) = x2*x2;   M3(2,2) = x2;   M3(2,3) = 1._dp
!   M3(3,1) = x3*x3;   M3(3,2) = x3;   M3(3,3) = 1._dp
!   call InverseMatrix3(M3,iM3)
!
!  do it = 2, NjT-1
!  do ip = 2, NkT-1
!
!   itp = itp_it_ip_lkp(it,ip)
!
!
!   ! dR/ds
!   y(1) = cR_hp(eR-1,itp);  y(2) = cR_hp(eR,itp);  y(3) = cR_hp(eR+1,itp)
!   call MultMatrixVect3(iM3,y,a)
!   dRR_hp(1,1) = 2._dp*a(1)*x2 + a(2)
!
!   ! dphi/ds
!   dRR_hp(2,1) = 0
!
!   ! dZ/ds
!   y(1) = cZ_hp(eR-1,itp);  y(2) = cZ_hp(eR,itp);  y(3) = cZ_hp(eR+1,itp)
!   call MultMatrixVect3(iM3,y,a)
!   dRR_hp(1,1) = 2._dp*a(1)*x2 + a(2)
!
!  enddo
!  enddo
!  enddo
!
!  ! Theta derivatives -------------------------
!  do it = 2, NjT-1
!
!   x1 = theta_T(it-1)
!   x2 = theta_T(it)
!   x3 = theta_T(it+1)
!   M3(1,1) = x1*x1;   M3(1,2) = x1;   M3(1,3) = 1._dp
!   M3(2,1) = x2*x2;   M3(2,2) = x2;   M3(2,3) = 1._dp
!   M3(3,1) = x3*x3;   M3(3,2) = x3;   M3(3,3) = 1._dp
!   call InverseMatrix3(M3,iM3)
!
!  do eR = 2, NelR-1
!  do ip = 2, NkT-1
!
!   itp = itp_it_ip_lkp(it,ip)
!
!
!   ! dR/dtheta
!   y(1) = cR_hp(eR,itp);  y(2) = cR_hp(eR,itp);  y(3) = cR_hp(eR+1,itp)
!   call MultMatrixVect3(iM3,y,a)
!   dRR_hp(1,1) = 2._dp*a(1)*x2 + a(2)
!
!   ! dphi/ds
!   dRR_hp(2,1) = 0
!
!   ! dZ/ds
!   y(1) = cZ_hp(eR-1,itp);  y(2) = cZ_hp(eR,itp);  y(3) = cZ_hp(eR+1,itp)
!   call MultMatrixVect3(iM3,y,a)
!   dRR_hp(1,1) = 2._dp*a(1)*x2 + a(2)
!
!  enddo
!  enddo
!  enddo
