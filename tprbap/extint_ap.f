!--------0---------0---------0---------0---------0---------0---------0--
!  15.12.89        LAST MODIFICATION 20.03.06     UHS: extint DATA D
!--------0---------0---------0---------0---------0---------0---------0--
!
      subroutine extint (ni, njk, lmnb, mm, nmin, nmax, nper, lmetpr,
     &rplmin, mb, nb, lfrz, tsin, tcos, ri, zi, phvi, fr, fz, fphv, fri,
     &fzi, fphvi, nbalmn)
!
!      INTERPOLATION TO INTEGER GRID FOR FRI, FZI AND FPHVI
!
C.. Implicits ..
      implicit none
C
C.. Formal Arguments .. 
      integer nmin,nmax,njk,lmnb,ni,mm,nper,lmetpr,mb(*),nb(*),
     &        lfrz(0:36,nmin:nmax),nbalmn
      real rplmin,tsin(njk,*),tcos(njk,*),ri(njk,0:*),zi(njk,0:*),
     &     phvi(njk,0:*),fr(lmnb,*),fz(lmnb,*),fphv(lmnb,*),
     &     fri(lmnb,0:*),fzi(lmnb,0:*),fphvi(lmnb,0:*)
C
C.. Local Scalars ..
      integer i,JK,l,m,n,inbmn
      real FIM,FIP,frzpmx
C
C.. Intrinsic Functions ..
      intrinsic MAX, abs
C
C ... Executable Statements ...
C
!
      do 80 l = 1, lmnb  
         do  m = 0, mm  
           do  n = nmin, nmax  
              if (mb (l) .eq.m.and.nb (l) .eq.n * nper) then  
              frzpmx = 1.e-99  
                   do i = 1, ni  
                   frzpmx = MAX(frzpmx,abs(fr(l,i)),abs(fz(l,i)),
     &                           abs(fphv(l,i)))
                   end do
                   if (frzpmx.gt.rplmin)then
                   lfrz (m, n) = 1  
                   end if  
              end if  
           end do
         end do
!
         if (mb (l) .eq.0) then  
           do  i = 1, ni - 1  
             fri (l, i) = 0.5 * (fr (l, i) + fr (l, i + 1) )  
             fzi (l, i) = 0.5 * (fz (l, i) + fz (l, i + 1) )  
           end do
         else  
           do  i = 1, ni - 1  
!
             inbmn = i + nbalmn - 1  
             fim = 0.5 * ( (inbmn) / (inbmn - 0.5) ) ** (mb(l)   / 2.)
             fip = 0.5 * ( (inbmn) / (inbmn + 0.5) ) ** (mb (l)  / 2.)
             fri (l, i) = fim * fr (l, i) + fip * fr (l, i + 1)  
             fzi (l, i) = fim * fz (l, i) + fip * fz (l, i + 1)  
           end do
         endif  
      do 20 i = 1, ni - 1  
         fphvi (l, i) = 0.5 * (fphv (l, i) + fphv (l, i + 1) )  
 20   end do
      if (mb (l) .eq.0.and.nb (l) .eq.0) then  
      do 25 i = 1, ni - 1  
         fphvi (l, i) = 0.  
 25   end do
      endif  
!--------0---------0---------0---------0---------0---------0---------0--
! extrapolation to origin and border
      IF(MB(L).EQ.0) THEN
         FRI(L,0) = 1.5 *  FR(L,1) - 0.5*  FR(L,2)
         FZI(L,0) = 1.5 *  FZ(L,1) - 0.5*  FZ(L,2)
         FPHVI(L,0) = 0.
      ELSE
         FRI(L,0) = 0.
         FZI(L,0) = 0.
         FPHVI(L,0) = 0.
      END IF

      FRI(L,NI) =       1.5 * FR(L,NI ) - 0.5 * FR(l,ni-1)
      FZI(L,NI) =       1.5 * FZ(L,NI ) - 0.5 * FZ(l,ni-1)
      FPHVI(L,NI) =     1.5 *FPHV(L,NI) - 0.5*FPHV(L,NI-1)
      

!
      if (lmetpr.eq.1) then  
      do 52 i = 0, 4  
          write (16, 1001) fri (l, i), fzi (l, i), fphvi (l, i),
     &       l, mb (l), nb (l), i
 52   end do
      do 53 i = ni - 4, ni  
          write (16, 1001) fri (l, i), fzi (l, i), fphvi (l, i),
     &       l, mb (l), nb (l), i
 53   end do
 1001 format ('fr/zi:' ,1p3e15.6,4i4)  
      endif  
!
!     RI, ZI, RSQI AND PHVI ON INTEGER GRID
!
       do 75 i = 0, ni   
!
        do jk = 1, njk  
           ri (jk, i) = ri (jk, i) + fri (l, i) * tcos (jk, l)  
           zi (jk, i) = zi (jk, i) + fzi (l, i) * tsin (jk, l)  
           phvi (jk, i) = phvi (jk, i) + fphvi (l, i) * tsin (jk, l)  
        end do
   75  end do  
   80 end do  
!
      return  
      end subroutine extint
