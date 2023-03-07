      subroutine getiota(phipog,guu,guv,lu,lv,wint,iotas,jv,
     &  czero,ns,ncurr)
!
C.. Implicits ..
      implicit none
      include 'name1.inc'
      integer               :: ns   , ncurr
      real                  :: czero
      real, dimension(*)    :: iotas, jv
      real, dimension(ns,*) :: guu  , guv  , lu, lv, wint, phipog
C.. Local Scalars ..
      integer :: js, lk
      real :: top, bot
        if(ncurr .ne. 0) then
         do 10 js = 2,ns
          top = jv(js)
          bot = czero
            do 20 lk = 1,nznt
              top = top - wint(js,lk)*(guu(js,lk)*lv(js,lk)
     &            +                    guv(js,lk)*lu(js,lk))
              bot = bot + wint(js,lk)*phipog(js,lk)*guu(js,lk)
   20       end do
         iotas(js) = top/bot
   10   end do
!************
!                 ADD IOTA TO LAMBDA (ON HALF MESH NOW)
!************
        endif
        do 45 js = 2,ns
          do 40 lk = 1,nznt
            lv(js,lk) = lv(js,lk) + phipog(js,lk)*iotas(js)
   40     end do
   45   end do
        return
        end subroutine getiota
