c        cliche name1
        integer nsd,mpol,nmax,ntheta,nzeta,nvac,nsd1,ntheta1,ntheta2,
     &          nznt,nrztd,mpol1,mpol2,nmax1,nmax2,nmax3,mnd,mnd1,
     &          mnd2,mnmax,neq
        parameter(nsd=97, mpol=10, nmax=0, ntheta=39, nzeta=1, nvac=0,
     &  nsd1=nsd+1, ntheta1=2*(ntheta/2), ntheta2=1+ntheta1/2, nznt=
     &  nzeta*ntheta2, nrztd=nznt*nsd+1, mpol1 = mpol-1,
     &  nmax1=1+nmax, mnd=mpol*nmax1, mnd1=mnd-1, mnd2=2*mnd,
     &  mpol2 = 2*mpol1,nmax2=2*nmax,nmax3=nmax2+1,
     &  mnmax=nmax1+mpol1*(1+2*nmax), neq=6*nsd*mnd)
c        endcliche
