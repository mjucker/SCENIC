      program eqtransf
      implicit none
   
      integer, parameter :: dp = SELECTED_REAL_KIND(10)

      integer, parameter :: fn = 2


 ! dimensions: (:,:) <-> njk,ni
      real(dp), dimension(:,:), allocatable ::
     $     tbjac, tbmod, tgssu, tgppl,tgtpl, tgttl, tgstl, tgssl,
     $     tgspl, trx, trz, tb2, phv
     $     ,sigbs
      
      real(dp), dimension(:,:), allocatable ::
     $     rmn,zmn,phimn,bsmn,fsigbs

      real(dp), dimension(:), allocatable ::
     $     s_T, ftp, fpp, pp, ci, cj, cip, cjp

      integer,allocatable :: mb(:),nm(:)

      integer ni, nj,nk,Njk, nper,lmnb
      integer is, jk
      integer i,l,j,k
      real(dp) theta,xi,twopi,cosarg,sinarg
      character(30) eqname

c      intrinsic cos


c  write(*,'(a40)',advance='no') 'Enter configuration name (max 30 char): '
c  read *, eqname
      eqname='JET_ICRH'

      open(fn, file='fort.38')
      read(fn, *) nper
!  nper = 1
      read(fn, *) ni
  
      allocate(s_T(ni));
      allocate(ftp(ni));
      allocate(fpp(ni));
      allocate(ci(ni));
      allocate(cj(ni));
      allocate(cip(ni));
      allocate(cjp(ni));

      do is = 1, ni
         read(fn, *)
     $        s_T(is),ftp(is),fpp(is),ci(is),cj(is),cip(is),cjp(is)
      enddo
      cip(1)=cip(2)
      cjp(1)=cjp(2)
      cip(ni)=cip(ni-1)
      cjp(ni)=cjp(ni-1)


! Reading values in the real space (Boozer coordinates).  One toroidal period!!!
      read(fn, *) ni, nk, nj
      nk = nk/nj
      Njk = nj*nk
c      print *, 'TERPSICHORE grid: (ni=',ni,',nj=',nj,',nk=',nk,')'

      allocate(tbjac(Njk,ni));
      allocate(tbmod(Njk,ni));
      allocate(tgssu(Njk,ni));
      allocate(tgppl(Njk,ni));
      allocate(tgtpl(Njk,ni));
      allocate(tgttl(Njk,ni));
      allocate(tgstl(Njk,ni));
      allocate(tgssl(Njk,ni));
      allocate(tgspl(Njk,ni));
      allocate(trx  (Njk,ni));
      allocate(trz  (Njk,ni));
      allocate(tb2  (Njk,ni));
      allocate(phv  (Njk,ni));
      allocate(sigBs(Njk,ni));

      do is = 1, ni
         do jk = 1, Njk
            read(fn, *)
     $           tbjac(jk,is),tgssu(jk,is),tgppl(jk,is),tgtpl(jk,is)
     $           ,tgttl(jk,is),tgstl(jk,is),tgssl(jk,is)
c     $           ,trx (jk,is),trz (jk,is),tb2(jk,is),phv(jk,is)
         enddo
c         tgstl(1,is)=(4.*tgstl(2,is)-tgstl(3,is))/3.
c         tgstl(Njk/2+1,is)=(4.*tgstl(Njk/2,is)-tgstl(Njk/2-1,is))/3.
      enddo
c      do jk=1,Njk
c         tgssu(jk,1)=(4.*tgssu(jk,2)-tgssu(jk,3))/3.
c         tgssu(jk,ni)=(4.*tgssu(jk,ni-1)-tgssu(jk,ni-2))/3.
c      enddo


      read(fn, *) ni,lmnb

      allocate(mb(lmnb),nm(lmnb)
     $     ,rmn(lmnb,ni),zmn(lmnb,ni),phimn(lmnb,ni),bsmn(lmnb,ni)
     $     ,fsigbs(lmnb,ni))


      do i=1,ni
         do l=1,lmnb
            read(fn, *) mb(l),nm(l),
     >           rmn(l,i),zmn(l,i),phimn(l,i),bsmn(l,i),fsigbs(l,i)
         enddo
      enddo


      close(fn)

c      print*, 'Converting into real space...'

      twopi=4.*acos(0.)

      do i=1,ni
         jk=0
         do k=1,nk
            xi=(k-1)*twopi/nper/nk
            do j=1,nj
               theta=(j-1)*twopi/nj
               jk=jk+1
               tb2(jk,i)=0.
               trx(jk,i)=0.
               trz(jk,i)=0.
               phv(jk,i)=0.
               sigBs(jk,i)=0.
               do l=1,lmnb
                  cosarg=cos(mb(l)*theta - nm(l)*xi)
                  sinarg=sin(mb(l)*theta - nm(l)*xi)
                  tb2(jk,i)=tb2(jk,i)+bsmn(l,i)*cosarg
                  trx(jk,i)=trx(jk,i)+rmn(l,i)*cosarg
                  trz(jk,i)=trz(jk,i)+zmn(l,i)*sinarg
                  phv(jk,i)=phv(jk,i)+phimn(l,i)*cosarg
                  sigBs(jk,i)=sigBs(jk,i)+fsigbs(l,i)*cosarg
               enddo
               tgspl(jk,i)=sigBs(jk,i)/ftp(i)
     $              - fpp(i)/ftp(i)*tgstl(jk,i)
            enddo
         enddo
      enddo
               

c      print *, 'Writing to unformatted file...'

      open(fn, file='fort.37.DAT.N',form='unformatted')

      write(fn) eqname  
      write(fn) ni,nj,nk,nper
      write(fn) s_T, ftp, fpp, ci, cj, cip, cjp
      write(fn) tbjac, tgssu, tgppl, tgtpl, tgttl, tgstl, tgssl, tgspl
     $     , trx, trz, tb2, phv

      close(fn)

c      print *, 'Done'


      deallocate(s_T);
      deallocate(ftp);
      deallocate(fpp);
      deallocate(ci);
      deallocate(cj);
      deallocate(cip);
      deallocate(cjp);
      
      deallocate(tbjac);
      deallocate(tbmod);
      deallocate(tgssu);
      deallocate(tgppl);
      deallocate(tgtpl);
      deallocate(tgttl);
      deallocate(tgstl);
      deallocate(tgssl);
      deallocate(tgspl);
      deallocate(tb2  );
      deallocate(phv  );
      
      end
