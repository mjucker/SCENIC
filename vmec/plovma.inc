c----6789012345678901234567890123456789012345678901234567890123456789012
      integer :: nit, njt, nkt, npert, mnt, nkst, njkt, njkst
        parameter(nit=128,njt=99,nkt=1,npert=1,mnt=10,nkst=nkt*npert+1
c        parameter(nit=64,njt=99,nkt=16,npert=4,mnt=72,nkst=nkt*npert+1
c        parameter(nit=96,njt=99,nkt=48,npert=03,mnt=196,nkst=nkt*npert+1
     &            ,njkt=(njt-1)*nkt,njkst=(njt-1)*nkt*npert)
c
