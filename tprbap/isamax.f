      function isamax (n,sx,incx)
c  ******************************************************************
c
c 1. function
c  isamax returns the index of element absolute mumimum value 
c 2.  description
c  n     ; in   : Number of elements to process in the vector to be searched
c                  (n=vector length if incx=1;n=vector length/2 if incx=2;
c                    and so on).
c                 if n<=0,isamax returns 0 .
c  sx    ; in   : Real vector to be searched
c  incx  ; in   : Increment between elements of sx
c  ******************************************************************
c
      real :: sx(*), vmax
      integer :: n, incx, ind, j,  i, isamax
      intrinsic abs
c
      ind=0
!      if(n.le.0) goto 20
      if(n.gt.0) then
        if(incx.ge.0) then
          j=1
        else
          j=1-incx*(n-1)
        endif
        vmax=abs(sx(j))
        ind=1
        do 10 i=2,n
          j=j+incx
!        if(vmax.ge.abs(sx(j))) goto 10
        if(vmax.lt.abs(sx(j))) then
          vmax=abs(sx(j))
          ind=i
        end if
   10   end do
        end if
   20   isamax=ind
      return
      end
