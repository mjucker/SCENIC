        subroutine trid(a,d,b,c,gam,alf,jmin,nn,ns)
!--------0---------0---------0---------0---------0---------0---------0-c
!
C.. Implicits ..
      implicit none
!
      include 'name1.inc'
!************
!                 SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,NN
!                 AND RETURNS ANSWER IN C(I)
!************
        integer                    :: nn,ns
        real, dimension(ns,0:mnd1) :: a,b,c,d
        real, dimension(0:mnd1,ns) :: alf,gam
        integer, dimension(0:3)    :: jmin
        integer, dimension(3)      :: mupper, mlower
        common /bounds/ mupper,mlower
C.. Local Scalars ..
      integer :: imodes, in, in1, mn, i, n2, i1
!************
!                 SEPARATE M=0 (IMODES=1), M=1 (IMODES=2), M>1 (IMODES=3)
!************
      do 100 imodes = 1,3
        in = jmin(imodes-1)
        in1 = in + 1
        do 10 mn = mlower(imodes),mupper(imodes)
          gam(mn,in) = d(in,mn)
   10   end do
        do 25 i=in1,nn
          do 20 mn = mlower(imodes),mupper(imodes)
            alf(mn,i-1) = a(i-1,mn)/gam(mn,i-1)
            gam(mn,i)   = d(i,mn) - b(i,mn)*alf(mn,i-1)
   20     end do
   25   end do
        do 30 mn = mlower(imodes),mupper(imodes)
          c(in,mn) = c(in,mn)/gam(mn,in)
   30   end do
        do 45 i=in1,nn
          do 40 mn = mlower(imodes),mupper(imodes)
            c(i,mn) = (c(i,mn) - b(i,mn)*c(i-1,mn))/gam(mn,i)
   40     end do
   45   end do
        n2 = nn + in
        do 55 i=in1,nn
          i1 = n2 -i
            do 50 mn = mlower(imodes),mupper(imodes)
               c(i1,mn) = c(i1,mn) - alf(mn,i1)*c(i1+1,mn)
   50       end do
   55   end do
  100 end do
        return
        end subroutine trid
