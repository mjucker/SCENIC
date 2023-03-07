!-----------------------------------------------------------------------
!  31.05.88        LAST MODIFICATION 20.02.96     UHS: INTEXT DATA   D
!-----------------------------------------------------------------------
!
      subroutine intext (ni, njk, lmnv, nbalmn, mx, fr, fz, frv, fzv,
     & tcos, tsin, r, z)
!
!
!.. Implicits ..
      implicit none  
!
!.. Formal Arguments ..
      integer :: lmnv, njk, ni, nbalmn, mx ( * )  
      real :: fr (lmnv, * ), fz (lmnv, * ), frv (lmnv, 0: * ), fzv (
     & lmnv, 0: * ), tcos (njk, * ), tsin (njk, * ), r (njk, 0: * ),
     & z (njk, 0: * )
!
!.. Local Scalars ..
      integer :: I, JK, l, nbegn, inbmn  
      real :: fim, fip  
!
! ... Executable Statements ...
!
!     COMPUTE R/Z ON HALFGRID USING INTER/EXTRAPOLATED VALUES OF FR/Z
!
!     INTERPOLATE FOURIER AMPL.
!
      nbegn = 1  
      do 200 l = 1, lmnv  
        if (nbalmn.eq.1) then  
          nbegn = 2  
          if (mx (l) .eq.0) then  
          fr (l, 1) = 0.5 * (frv (l, 1) + frv (l, 0) )  
          fz (l, 1) = 0.5 * (fzv (l, 1) + fzv (l, 0) )  
         else  
          fim = 0.5** (mx (l) / 2.)  
          fr (l, 1) = fim * frv (l, 1)  
          fz (l, 1) = fim * fzv (l, 1)  
         endif  
        endif  
!
      do 110 i = nbegn, ni  
!
       inbmn = i + nbalmn - 1  
       fip = 0.5 * ( (inbmn - .5) / (inbmn) ) ** (mx (l) / 2.)  
       fim = 0.5 * ( (inbmn - .5) / (inbmn - 1) ) ** (mx (l) / 2.)  
       fr (l, i) = fip * frv (l, i) + fim * frv (l, i - 1)  
       fz (l, i) = fip * fzv (l, i) + fim * fzv (l, i - 1)  
  110   end do  
!
!     REAL SPACE
!
      do 120 i = 1, ni  
        do 120 jk = 1, njk  
          r (jk, i) = r (jk, i) + tcos (jk, l) * fr (l, i)  
  120   z (jk, i) = z (jk, i) + tsin (jk, l) * fz (l, i)  
!
  200 end do  
!
      return  
      end subroutine intext
