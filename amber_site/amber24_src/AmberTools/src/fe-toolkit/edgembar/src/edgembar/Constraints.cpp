#include <cstdlib>
#include <iostream>

#include "Constraints.hpp"

extern "C"
{

  void dgemv_( char const * TRANS,
               int const * M,
               int const * N,
               double const * ALP,
               double const * A,
               int const * LDA,
               double const * X,
               int const * INCX,
               double const * BETA,
               double * Y,
               int const * INCY );

  void dgesvd_( char const * JOBU,
                char const * JOBVT,
                int const * M,
                int const * N,
                double * A,
                int const * LDA,
                double * S,
                double * U,
                int const * LDU,
                double * VT,
                int const * LDVT,
                double * WORK,
                int const * LWORK,
                int * INFO );
}


edgembar::Constraints::Constraints()
  : mNumParam(0),
    mNumFreeParam(0),
    mNumCon(0)
{

}

edgembar::Constraints::Constraints
( int nparam )
  : mNumParam(nparam),
    mNumFreeParam(nparam),
    mNumCon(0)
{

}


void edgembar::Constraints::reset( int nparam )
{
  mNumParam = nparam;
  mNumFreeParam = nparam;
  mNumCon = 0;
  mConMat.resize(0);
  mConVals.resize(0);
  VT.resize(0);
  P0.resize(0);
}


void edgembar::Constraints::push_back
( std::vector<int> const & pidxs,
  std::vector<double> const & coefs,
  double const cval )
{
  if ( pidxs.size() > 0  )
    {
      
      if ( mNumFreeParam == 0 )
	{
	  std::cerr << "Error: Could not push restraint because there are "
		    << "already no free parameters\n";
	  std::exit(EXIT_FAILURE);
	}
      
      if ( pidxs.size() != coefs.size() )
	{
	  std::cerr << "Error in edgembar::Constraints::push_back : "
		    << " pidxs.size() != coefs.size() : "
		    << pidxs.size() << " " << coefs.size()
		    << "\n";
	  std::exit(EXIT_FAILURE);
	}
      
      mNumFreeParam -= 1;
      mNumCon += 1;

      for ( int i=0; i<mNumParam; ++i )
	{
	  mConMat.push_back( 0. );
	}
      mConVals.push_back( cval );

      int const icon = mNumCon - 1;
      
      for ( int k=0, nk=pidxs.size(); k<nk; ++k )
	{
	  
	  if ( pidxs[k] >= mNumParam or pidxs[k] < 0 )
	    {
	      std::cerr << "Error in edgembar::Constraints::push_back : "
			<< " pidxs[" << k << "]=" << pidxs[k]
			<< " is out of bounds [0," << mNumParam-1
			<< "\n";
	      std::exit(EXIT_FAILURE);
	    }
	  int i = pidxs[k];
	  mConMat[i+icon*mNumParam] = coefs[k];
	};
    };
}


void edgembar::Constraints::finalize()
{
  if ( mNumCon > 0 )
    {
      
      std::vector<double> A(mNumCon*mNumParam,0.);
      for ( int i=0; i<mNumCon; ++i )
        for ( int j=0; j<mNumParam; ++j )
          A[i+j*mNumCon] = mConMat[j+i*mNumParam];

      std::vector<double> S(mNumParam,0.);
      std::vector<double> U(mNumCon*mNumCon,0.);
      VT.assign(mNumParam*mNumParam,0.);

      std::vector<double> WORK(1,0.);
      int LWORK = -1;
      int INFO = 0;
      
      dgesvd_("A","A",&mNumCon,&mNumParam,
              A.data(),&mNumCon,
              S.data(),
              U.data(),&mNumCon,
              VT.data(),&mNumParam,
              WORK.data(),
              &LWORK,&INFO);
      
      LWORK = 0.5+WORK[0];
      WORK.resize(LWORK);
      
      dgesvd_("A","A",&mNumCon,&mNumParam,
              A.data(),&mNumCon,
              S.data(),
              U.data(),&mNumCon,
              VT.data(),&mNumParam,
              WORK.data(),
              &LWORK,&INFO);

      for ( int i=0; i<mNumCon; ++i )
        {
          if ( S[i] > 1.e-12 )
            {
              S[i] = 1./S[i];
            }
          else
            {
              S[i] = 0.;
            };
        };

      std::vector<double> t(mNumParam,0.);
      double alpha = 1.;
      double beta  = 0.;
      int inc = 1;
      dgemv_("T",&mNumCon,&mNumCon,&alpha,
             U.data(),&mNumCon,
             mConVals.data(),&inc,
             &beta,
             t.data(),&inc);
      for ( int i=0; i<mNumCon; ++i )
        t[i] *= S[i];
      
      
      P0.assign(mNumParam,0.);
      dgemv_("T",&mNumParam,&mNumParam,&alpha,
             VT.data(),&mNumParam,
             t.data(),&inc,
             &beta,
             P0.data(),&inc);

      /*
      for ( int i=0; i<mNumParam; ++i )
        {
          std::printf("P0 %20.10e\n",P0[i]);
        }
      for ( int i=0; i<mNumParam; ++i )
        {
          for ( int j=mNumCon; j<mNumParam; ++j )
            std::printf("%12.3e",VT[j+i*mNumParam]);
          std::printf("\n");
        }
      */
      
      /*
      for ( int ii=0; ii<mNumParam-mNumCon; ++ii )
        {
          std::fill(t.begin(),t.end(),0.);
          t[mNumCon+ii] = 1.;
          dgemv_("T",&mNumParam,&mNumParam,&alpha,
                 VT.data(),&mNumParam,
                 t.data(),&inc,
                 &alpha,
                 p0.data(),&inc);
      
          std::vector<double> cv(mNumCon,0.);
          for ( int i=0; i<mNumCon; ++i )
            for ( int j=0; j<mNumParam; ++j )
              cv[i] += mConMat[j+i*mNumParam] * p0[j];
        
          for ( int i=0; i<mNumCon; ++i )
            {
              std::printf("%20.10e %20.10e\n",cv[i],ConVals[i]);
            }
          std::printf("\n");
        };
      
      std::exit(1);
      */

      
    }

}


std::vector<double> edgembar::Constraints::GetParamsFromFree
( std::vector<double> const & qin ) const
{
  std::vector<double> p(mNumParam,0.);
  if ( mNumCon == 0 )
    {
      p = qin;
    }
  else
    {
      // In the way that we are doing thigs, we
      // didn't eliminate the first mNumCon rows of
      // VT, but we effectively do so by passing it
      // a q vector whose first mNumCon elements are
      // zero
      std::vector<double> q( mNumParam, 0. );
      for ( int i=0; i < mNumFreeParam; ++i )
        q[mNumCon+i] = qin[i];
      
      double alpha = 1.;
      double beta  = 0.;
      int inc = 1;
      dgemv_("T",&mNumParam,&mNumParam,&alpha,
             VT.data(),&mNumParam,
             q.data(),&inc,
             &beta,
             p.data(),&inc);
      for ( int i=0; i<mNumParam; ++i )
        p[i] += P0[i];
    }
  return p;
}


std::vector<double> edgembar::Constraints::GetFreeFromParams
( std::vector<double> const & p ) const
{
  std::vector<double> q(mNumFreeParam,0.);
  if ( mNumCon == 0 )
    {
      q = p;
    }
  else
    {
      double alpha = 1.;
      double beta  = 0.;
      int inc = 1;
      std::vector<double> tmp( mNumParam, 0. );
      dgemv_("N",&mNumParam,&mNumParam,&alpha,
             VT.data(),&mNumParam,
             p.data(),&inc,
             &beta,
             tmp.data(),&inc);
      for ( int i=0; i<mNumFreeParam; ++i )
        q[i] = tmp[mNumCon+i];
    }
  return q;
}

