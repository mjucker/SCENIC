c        cliche name2
	include 'name0.inc'
        include 'name1.inc'
c        use name0
c        use name1
        real, dimension(nsd) :: buco,bvco,ju,jv
        common/current/buco,bvco,ju,jv
        real, dimension(nznt) :: potvac,xmpot,xnpot,brv,bphiv,bzv       vac
        real :: bscale                                                  vac
        integer :: mpmax, ivac, ivac2                                   vac
        common/extfld/potvac,xmpot,xnpot,brv,bphiv,bzv,bscale,mpmax,    vac
     >                ivac,ivac2                                        vac
        real :: fnorm,fsqr,fsqz,fsql,fnorm1,fsqr1,fsqz1,fsql1,fsq
        real :: fedge,wb,wpar,wper
        real, dimension(100) :: fsqt,wdot
        real, dimension(nsd) :: equif
        common/fsqu/fnorm,fsqr,fsqz,fsql,fnorm1,fsqr1,fsqz1,fsql1,fsq,
     >         fedge,wb,wpar,wper,fsqt,wdot,equif
        real, dimension(0:10) ::am,ah,at,ai
        real, dimension(0:nmax) :: raxis,zaxis
        real, dimension(0:nmax,0:mpol1,2) :: rb,zb
        real :: ftol,gam,bcrit
        integer :: ncurr,nfp,nstep,niter,nvacskip
        common/inputdat/am,ah,at,ai,raxis,zaxis,rb,zb,kpress,
     >         ftol,gam,bcrit,ncurr,nfp,nstep,niter,nvacskip
        real, dimension(0:nmax,0:mpol1) :: xrz3,xrz4
        real, dimension(0:mpol1,3) :: xmpq
        real, dimension(0:mpol1) :: mscale
        real, dimension(0:nmax1) :: nscale
        integer, dimension(0:mpol) :: jmin1,jmin2,jlam
	integer :: kpress
        common/mnarray/xrz3,xrz4,xmpq,mscale,nscale,jmin1,jmin2,jlam
        real, dimension(nsd) :: iotas,mass,ph,pres,phot,pm,bmin,vp,pppr
     >                         ,tpotb,pd
        real, dimension(nsd1) :: phips
        real, dimension(nrztd) :: phip
        common/profs/iotas,mass,ph,phips,pres,phot,pm,pppr,bmin,vp,tpotb
     >        ,pd,phip
        real :: dnorm,hs,ohs,twopi,voli
        integer :: ijacob,itfsq,iequi,irst,iter1,iter2,isigng
        integer :: meven,modd,ndamp,ns,ns4,neqs,nrzt,mns
        common/scalars/dnorm,hs,ohs,twopi,voli,ijacob,itfsq,iequi,
     >         irst,iter1,iter2,isigng,meven,modd,ndamp,ns,ns4,neqs,
     >         nrzt,mns
        real, dimension(mnd2*nsd) :: faclam
        real, dimension(nrztd) :: shalf,sqrts,wint
        real, dimension(neq) :: scalxc
        common/scalefac/faclam,shalf,sqrts,scalxc,wint
        real, dimension(0:nmax,0:mpol1) :: faccon
        real, dimension(nsd) :: specw, tcon
        common/spectra/faccon,specw,tcon
        real, dimension(2*nsd) :: otau
        real, dimension(0:10) :: timer
        real :: delt, otav
        common/time/delt,otav,otau,timer
        real, dimension(ntheta2,0:mpol1) :: cosmu,sinmu,cosmum,sinmum,
     >                                      cosmui,cosmumi
        real, dimension(nzeta,0:nmax) :: cosnv,sinnv,cosnvn,sinnvn
        common/trignew/cosmu,sinmu,cosmum,sinmum,cosmui,cosmumi,
     >                 cosnv,sinnv,cosnvn,sinnvn
        real, dimension(nznt,3) :: bsqsav                               vac
        real, dimension(nznt) :: dbsq,bsqvac,rbsq                       vac
        real :: curpol,curtor,rbtor,ctor,phiedge,delbsq                 vac
        common/magfield/bsqsav,dbsq,bsqvac,rbsq,curpol,curtor,rbtor,    vac
     >                                          ctor,phiedge,delbsq     vac
        real, dimension(neq) :: xcdot,xstore,xc,gc
        common/xstuff/xcdot,xstore,xc,gc
c        endcliche
