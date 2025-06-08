
      double precision ff,evec,ww,vv,zz
      common /se_work/ ff(msorb2),evec(msrorb2),ww(1+6*msorb+2*msorb2),
     &                 vv(msorb), zz(msorb2)
      integer iww, izz
      common /se_worki/ iww(10*msorb), izz(2*msorb)
