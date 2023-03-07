
        subroutine difvol(gsqrt,wint,vp,dnorm,nznt,ns)
        integer :: nznt, ns, js
        real :: dnorm
        real, dimension(ns,*) :: gsqrt,wint
        real, dimension (*)   :: vp
C.. External Functions ..
      real  sdot
      external  sdot
C.. Intrinsic Functions ..
      intrinsic  abs
c
        do 10 js = 2,ns
         vp(js) = sdot(nznt,gsqrt(js,1),ns,wint(js,1),ns)
 10     end do
        do 20 js = 2,ns
         vp(js) = dnorm*abs(vp(js))
 20     end do
        return
        end subroutine difvol
