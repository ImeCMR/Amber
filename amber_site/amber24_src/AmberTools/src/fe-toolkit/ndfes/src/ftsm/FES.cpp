#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <cstdio>
#include <random>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>

#include "FES.hpp"
#include "Bspline.hpp"
#include "MeshUtils.hpp"
#include "PeriodicUtils.hpp"
#include "Utils.hpp"
#include "TubeUtils.hpp"
#include "../xmlio/xmlio.hpp"



namespace strop
{
  inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    s.erase(0, s.find_first_not_of(t));
    return s;
  }
  
  // trim from right
  inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    s.erase(s.find_last_not_of(t) + 1);
    return s;
  }
  
  // trim from left & right
  inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    return ltrim(rtrim(s, t), t);
  } 
}




extern "C"
{
  void dsysv_(const char* uplo, const int* n,
              const int* nrhs,
              double* a, const int* lda, int* ipiv,
              double* b, const int* ldb,
              double* work, const int* lwork, int* info);

  void dgemv_(const char *trans, const int *m, const int *n,
	      double const *alpha, double const *a, const int *lda,
	      double const *x, const int *incx, double const *beta,
	      double *y, const int *incy);
  
  

/* DSYTRF - compute the factorization of a real symmetric matrix */
/* A using the Bunch-Kaufman diagonal pivoting method */


/** \private */
  void
  dsytrf_(const char* uplo, const int* n,
	  double* a, const int* lda, int* ipiv,
	  double* work, const int* lwork, int* info);

/* DSYTRI - compute the inverse of a real symmetric indefinite */
/* matrix A using the factorization A = U*D*U**T or A = L*D*L**T */
/* computed by DSYTRF */


/** \private */
  void
  dsytri_(const char* uplo, const int* n,
	  double* a, const int* lda, const int* ipiv,
	  double* work, int* info);

/* DSYTRS - solve a system of linear equations A*X = B with a */
/* real symmetric matrix A using the factorization A = U*D*U**T or */
/* A = L*D*L**T computed by DSYTRF */


  void
  dsymm_(const char *side, const char *uplo, const int *m,
  	 const int *n, double const *alpha,
  	 double const *a, const int *lda,
  	 double const *b, const int *ldb,
  	 double const *beta, double *c, const int *ldc);

  void
  dgemv_(const char *trans, const int *m, const int *n,
  	 double const *alpha, double const *a, const int *lda,
  	 double const *x, const int *incx, double const *beta,
  	 double *y, const int *incy);
  
}



int sym_inverse( int const nf,  double * A );

void ge_dot_sy
( int nfA, int nsA, double const * A, 
  int nfB, int nsB, double const * B, 
  double * AB,
  double alpha = 1.0, double beta = 0. );

void ge_dot_v
( int nfA, int nsA, double const * A, double const * x,double * Ax,
  double alpha = 1., double beta  = 0. );





int sym_inverse( int const nf, double * A )
{

  std::vector<int> ipiv(nf,0);
  std::vector<double> work(1,0);
  int lwork=-1;
  int info=0;
  dsytrf_( "U", &nf, A, &nf, ipiv.data(), work.data(), &lwork, &info );

  if ( info == 0 )
    {
      lwork = work[0];
      work.resize(lwork);
      dsytrf_( "U", &nf, A, &nf, ipiv.data(), work.data(), &lwork, &info );
      if ( info != 0 )
	{
	  std::printf("dsytrf failed with info=%i\n",info);
	}
      else
	{
	  work.resize(nf);
	  dsytri_( "U", &nf, A, &nf, ipiv.data(), work.data(), &info );
	  if ( info != 0 )
	    {
	      std::printf("dsytri failed with info=%i\n",info);
	    }
	}
    }
  else
    {
      std::printf("dsytrf failed to get lwork with info=%i\n",info);
    }

  if ( info == 0 )
    {
      for ( int j=1; j < nf; ++j )
	{
	  for ( int i=0; i < j; ++i )
	    {
	      A[j+i*nf] = A[i+j*nf];
	    }
	}
    }
  return info;
}


void ge_dot_sy
( int nfA, int , double const * A, 
  int nfB, int nsB, double const * B, 
  double * AB,
  double alpha, double beta )
{
  int nf = nfA;
  int ns = nsB;
  dsymm_("R","U",&nf,&ns,&alpha, B,&nfB, A,&nfA,&beta, AB,&nf);
}



void ge_dot_v
( int nfA, int nsA, double const * A, 
  double const * x,
  double * Ax,
  double alpha, double beta )
{
  int const inc = 1;
  dgemv_("N",&nfA,&nsA,&alpha,A,&nfA,x,&inc,&beta,Ax,&inc);
}







void RbfSolve( int const n, double * A, double * w )
{
  std::vector<int> ipiv( n, 0 );
  int LWORK = -1;
  int INFO = 0;
  //double alpha = 1.;
  //double beta  = 0.;
  int inc = 1;
  std::vector<double> WORK(1,0);
  
  dsysv_("U",&n,&inc,A,&n,ipiv.data(),w,
	 &n,WORK.data(),&LWORK,&INFO);
  
  LWORK = WORK[0]+1;
  WORK.resize( LWORK );
  
  dsysv_("U",&n,&inc,A,&n,ipiv.data(),w,
	 &n,WORK.data(),&LWORK,&INFO);
  
  if ( INFO > 0 )
    {
      std::cerr << "Failed to solve the RBF system equations using dsyev_."
		<< " Returned INFO=" << INFO << std::endl;
      std::exit(EXIT_FAILURE);
  }

}

void RbfSolveWithInverse( int const n, double * A, double * w )
{
  // on input, A is the RBF matrix
  // on output, A is the inverse of A
  // on input, w is the control values
  // on output, w is the weights Ainv.b
  sym_inverse(n,A);
  std::vector<double> b(w,w+n);
  ge_dot_v(n,n,A,b.data(),w);
}




static void GaussLegendreRule
( double const x1,
  double const x2,
  int const n,
  double * rVecX, 
  double * rVecW )
{
  int const MAXIT=10;
  double const EPS = 3.0e-14;
  int const m = (n+1)/2;
  double const xm = 0.5 * ( x2 + x1 );
  double const xl = 0.5 * ( x2 - x1 );  
  double const C1 = M_PI / (n+0.5);
  
  std::vector<double> z(m);
  std::vector<double> z1(m);
  std::vector<double> p1(m);
  std::vector<double> p2(m);
  std::vector<double> p3(m);
  std::vector<double> pp(m);
  std::vector<bool> unfinished(m);

  for ( int i=0; i<m; ++i )
    z[i] = std::cos( C1 * ( i + 0.75 ) );
  for ( int i=0; i<m; ++i )
    unfinished[i] = true;
  for ( int its=0; its < MAXIT; ++its )
    {
      for ( int i=0; i < m; ++i )
        {
          if ( unfinished[i] )
            {
              p1[i] = 1.0;
              p2[i] = 0.0;
            };
        };
      for ( int j=0; j < n; ++j )
        {
          for ( int i=0; i < m; ++i )
            {
              if ( unfinished[i] )
                {
                  p3[i]=p2[i];
                  p2[i]=p1[i];
                  p1[i]=((2.0*j+1.0)*z[i]*p2[i]-j*p3[i])/(j+1.0);
                };
            };
        };

      bool HasUnfinished = false;
      for ( int i=0; i < m; ++i )
        {
          if ( unfinished[i] )
            {
              pp[i]=n*(z[i]*p1[i]-p2[i])/(z[i]*z[i]-1.0);
              z1[i]=z[i];
              z[i]=z1[i]-p1[i]/pp[i];
              unfinished[i] = ( std::abs(z[i]-z1[i]) > EPS );
              if ( unfinished[i] )
                HasUnfinished = true;
            };
        };
      if ( ! HasUnfinished ){ break; };
      
      if ( its == MAXIT )
        std::cerr << "Too many iterations in GaussLegendreRule" << std::endl;
    };

  for ( int i=0; i < m; ++i )
    rVecX[i] = xm-xl*z[i];
  for ( int i=0; i < n-m; ++i )
    rVecX[n-i-1] = xm+xl*z[i];
  for ( int i=0; i < m; ++i )
    rVecW[i] = 2.0*xl/((1.-z[i]*z[i])*pp[i]*pp[i]);
  for ( int i=0; i < n-m; ++i )
    rVecW[n-i-1] = rVecW[i];
}







void RbfEval( int const ndim,
	      int const ncpts,
	      double const * cpts,
	      double const * wts,
	      double const shape,
	      int const npts,
	      double const * pts,
	      double * vals )
{
  if ( npts > 0 )
    {

      std::vector<double> A(npts*ncpts,0);

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif      
      for ( int i=0; i<ncpts; ++i )
	{
	  for ( int j=0; j<npts; ++j )
	    {
	      double r2 = 0;
	      for ( int k=0; k<ndim; ++k )
		{
		  double dx = pts[k+j*ndim] - cpts[k+i*ndim];
		  r2 += dx*dx;
		};
	      A[j+i*npts] = std::sqrt( 1 + shape * r2 );
	    }
	}
      
      double alpha=1;
      double beta=0;
      int inc=1;
      dgemv_( "N", &npts, &ncpts, &alpha, A.data(), &npts,
	      wts,  &inc, &beta, vals, &inc );
    };
  
}




void RbfEval( int const ndim,
	      int const ncpts,
	      double const * cpts,
	      double const * wts,
	      double const shape,
	      int const npts,
	      double const * pts,
	      double * vals,
	      double * grds )
{
  std::fill( vals, vals + npts, 0. );
  std::fill( grds, grds + npts*ndim, 0. );

  
  if ( npts > 0 )
    {

      //std::vector<double> A(npts*ncpts,0);
      //std::vector<double> D(ndim*npts*ncpts,0);

#ifdef WITH_OPENMP
#pragma omp parallel
#endif
      {
	std::vector<double> myvals( npts, 0 );
	std::vector<double> mygrds( ndim*npts, 0 );
#ifdef WITH_OPENMP
#pragma omp for schedule(dynamic)
#endif      
	for ( int i=0; i<ncpts; ++i )
	  {
	    for ( int j=0; j<npts; ++j )
	      {
		double r2 = 0;
		for ( int k=0; k<ndim; ++k )
		  {
		    double dx = pts[k+j*ndim] - cpts[k+i*ndim];
		    r2 += dx*dx;
		  };
		double f = std::sqrt( 1 + shape * r2 );
		double v = wts[i] * f;
		myvals[j] += v;
		
		double dvdr2 = 0.5 * wts[i] * shape / f;
		for ( int k=0; k<ndim; ++k )
		  {
		    double dx = pts[k+j*ndim] - cpts[k+i*ndim];
		    mygrds[k+j*ndim] += dvdr2 * 2 * dx;
		  }
	      }
	  }
#ifdef WITH_OPENMP
#pragma omp critical
#endif
	{
	  for ( int j=0; j<npts; ++j )
	    {
	      vals[j] += myvals[j];
	      for ( int k=0; k<ndim; ++k )
		{
		  grds[k+j*ndim] += mygrds[k+j*ndim];
		}
	    }
	}

      }
    }
  
}






ndfes::FES::FES()
  : mIsMBAR(false),
    mInterpMode(0),
    mNumBins(0),
    mOrder(0),
    mNBspl(0),
    mNumCorners(0),
    mRBFShape(100),
    mNumQuad(0),
    mSubBinSize(0),
    mARBFDelta(-1),
    mWithInverse(false)
{

}


void ndfes::FES::SetupHist()
{
  mInterpMode = 1;
}

void ndfes::FES::SetupRBF( double const shape, std::size_t const minsize, double const maxerr, bool const with_inverse )
{
  mInterpMode = 2;
  mRBFShape = shape;
  mWithInverse = with_inverse;
  std::vector<std::size_t> occbins;
  std::vector<std::size_t> unoccbins;
  for ( std::size_t ibin=0; ibin < mNumBins; ++ibin )
    {
      if ( mBinSizes[ibin] >= minsize )
	{
	  occbins.push_back( ibin );
	}
      else
	{
	  unoccbins.push_back( ibin );
	}
    }
  std::size_t nocc = occbins.size();
  std::size_t nunocc = unoccbins.size();
  std::size_t ndim = mDimInfo.GetNumDims();


  std::cout << "There are " << mNumBins << " occupied bins, "
	    << nunocc << " of which are excluded because they "
	    << "have fewer than " << minsize << " samples\n";
  

  mRBFCpts.resize( ndim*nocc );
  mRBFCwts.resize( nocc );
  mRBFCerr.resize( nocc );
  for ( std::size_t i=0; i<nocc; ++i )
    {
      double const * ci = mBins[ occbins[i] ].center.data();
      for ( std::size_t k=0; k<ndim; ++k )
	{
	  mRBFCpts[k+i*ndim] = ci[k];
	}
      mRBFCwts[i] = mBinVals[ occbins[i] ];
      mRBFCerr[i] = mBinErrs[ occbins[i] ];
    }


  

  for ( std::size_t d=0; d<ndim; ++d )
    {
      if ( mDimInfo.IsPeriodic(d) )
	{
	  std::vector<double> cp;
	  std::vector<double> cm;
	  std::vector<double> vp;
	  std::vector<double> vm;
	  std::vector<double> ep;
	  std::vector<double> em;
	  std::size_t ncpts = mRBFCwts.size();
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      mRBFCpts[d+i*ndim] = wraptorange(mRBFCpts[d+i*ndim],0.,360.);
	      
	      if ( mRBFCpts[d+i*ndim] < 0.07*360 )
		{
		  vp.push_back(mRBFCwts[i]);
		  ep.push_back(mRBFCerr[i]);
		  for ( std::size_t m=0; m<ndim; ++m )
		    {
		      if ( m == d )
			{
			  cp.push_back( mRBFCpts[m+i*ndim] + 360 );
			}
		      else
			{
			  cp.push_back( mRBFCpts[m+i*ndim] );
			}
		    }
		}
	      
	      if ( mRBFCpts[d+i*ndim] > 0.93*360 )
		{
		  vm.push_back(mRBFCwts[i]);
		  em.push_back(mRBFCerr[i]);
		  for ( std::size_t m=0; m<ndim; ++m )
		    {
		      if ( m == d )
			{
			  cm.push_back( mRBFCpts[m+i*ndim] - 360 );
			}
		      else
			{
			  cm.push_back( mRBFCpts[m+i*ndim] );
			}
		    }
		}
	    }
	  
	  std::copy( cp.begin(), cp.end(), std::back_inserter(mRBFCpts) );
	  std::copy( cm.begin(), cm.end(), std::back_inserter(mRBFCpts) );
	  std::copy( vp.begin(), vp.end(), std::back_inserter(mRBFCwts) );
	  std::copy( vm.begin(), vm.end(), std::back_inserter(mRBFCwts) );
	  std::copy( ep.begin(), ep.end(), std::back_inserter(mRBFCerr) );
	  std::copy( em.begin(), em.end(), std::back_inserter(mRBFCerr) );

	}
    }



  
  std::size_t ncpts = mRBFCwts.size();
  //std::printf("ncpts=%lu\n",ncpts);

  // for ( std::size_t i=0; i<ncpts; ++i )
  //   {
  //     std::printf("mRBFCwts %6lu %12.3e\n",i,mRBFCwts[i]);
  //   }
  
  // for ( std::size_t i=0; i<ncpts; ++i )
  //   {
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf("mRBFCpts %3lu %6lu %12.3e\n",dim,i,mRBFCpts[dim+i*ndim]);
  // 	};
  //   }
  // std::cout << "shape = " << shape
  // 	    << " with inverse = " << with_inverse
  // 	    << std::endl;

  
  std::vector<double> A( ncpts*ncpts, 0. );
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for ( std::size_t i=0; i<ncpts; ++i )
    {
      A[i+i*ncpts] = 1;
      
      double const * ci = mRBFCpts.data() + i*ndim;
      for ( std::size_t j=0; j<i; ++j )
	{
	  double const * cj = mRBFCpts.data() + j*ndim;

	  double r2 = 0;
	  for ( std::size_t k=0; k<ndim; ++k )
	    {
	      double dx = ci[k]-cj[k];
	      r2 += dx*dx;
	    }
	  // if ( r2 < 0.15*0.15-1.e-14 )
	  //   {
	  //     std::printf("%4lu %4lu r2 %14.5e\n",i,j,r2);
	  //   }
	  double f = std::sqrt( 1 + shape * r2 );
	  A[j+i*ncpts] = f;
	  A[i+j*ncpts] = f;
	}
    }

  if ( with_inverse )
    {
      mRBFAinv = A;
      RbfSolveWithInverse( ncpts, mRBFAinv.data(), mRBFCwts.data() );
    }
  else
    {
      RbfSolve( ncpts, A.data(), mRBFCwts.data() );
    }


  // for ( std::size_t i=0; i<ncpts; ++i )
  //   {
  //     std::printf("mRBFCwts %6lu %12.3e\n",i,mRBFCwts[i]);
  //   }
  

  
  if ( nunocc > 0 )
    {
      std::vector<double> unoccpts( ndim*nunocc, 0 );
      std::vector<double> unoccvals( nunocc, 0 );
      for ( std::size_t i=0; i<nunocc; ++i )
	{
	  unoccvals[i] = mBinVals[ unoccbins[i] ];
	  double const * ci = mBins[ unoccbins[i] ].center.data();
	  for ( std::size_t k=0; k<ndim; ++k )
	    {
	      unoccpts[k+i*ndim] = ci[k];
	    }
	}
      
      std::vector<double> pvals( nunocc, 0 );
      EvalRBF( nunocc, unoccpts.data(), pvals.data() );


      std::vector<std::size_t> glbskips;
      
      for ( std::size_t i=0; i<nunocc; ++i )
	{
	  double d = pvals[i]-unoccvals[i];
	  if ( std::abs(d) > maxerr )
	    {
	      glbskips.push_back( mBins[ unoccbins[i] ].glbidx );
	      //std::printf("%8lu %14.5e %14.5e %14.5e\n",i,pvals[i],unoccvals[i],d);
	    }
	}

      typedef std::vector<std::size_t>::iterator it;
      for ( it p=glbskips.begin(), pend=glbskips.end(); p != pend; ++p )
	{
	  std::size_t idx = mBins.size();
	  for ( std::size_t k=0; k<mBins.size(); ++k )
	    {
	      if ( mBins[k].glbidx == *p )
		{
		  idx = k;
		  break;
		}
	    }
	  
	  if ( idx != mBins.size() )
	    {
	      mNumBins -= 1;
	      mBins.erase( mBins.begin() + idx );
	      mBinVals.erase( mBinVals.begin() + idx );
	      mBinErrs.erase( mBinErrs.begin() + idx );
	      mBinREs.erase( mBinREs.begin() + idx );
	      mBinSizes.erase( mBinSizes.begin() + idx );
	    }
	  else
	    {
	      std::cerr << "Failed to remove bin " << *p << std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	}

      
      std::cout << "There are " << mNumBins
		<< " remaining bins that agree with the RBF within "
		<< std::fixed << std::setprecision(1) << maxerr << " kcal/mol"
		<< "\n";

    };
  
}



void ndfes::FES::EvalRBF( std::size_t const npts, double const * pts, double * vals ) const
{
  if ( npts > 0 )
    {
      int const ndim = mDimInfo.GetNumDims();
      int const ncpts = mRBFCwts.size();
      RbfEval( ndim, ncpts, mRBFCpts.data(), mRBFCwts.data(), mRBFShape, npts, pts, vals );
    };
}

void ndfes::FES::EvalRBF( std::size_t const npts, double const * pts, double * vals, double * grds ) const
{
  if ( npts > 0 )
    {
      int const ndim = mDimInfo.GetNumDims();
      int const ncpts = mRBFCwts.size();
      RbfEval( ndim, ncpts, mRBFCpts.data(), mRBFCwts.data(), mRBFShape, npts, pts, vals, grds );
    };
}


void ndfes::FES::EvalRBF( std::vector<double> const & cpts, std::vector<double> const & cwts, std::size_t const npts, double const * pts, double * vals ) const
{
  if ( npts > 0 )
    {
      int const ndim = mDimInfo.GetNumDims();
      int const ncpts = cwts.size();
      RbfEval( ndim, ncpts, cpts.data(), cwts.data(), mRBFShape, npts, pts, vals );
    };
}


void ndfes::FES::EvalRBF( std::vector<double> const & cpts, std::vector<double> const & cwts, std::size_t const npts, double const * pts, double * vals, double * grds ) const
{
  if ( npts > 0 )
    {
      int const ndim = mDimInfo.GetNumDims();
      int const ncpts = cwts.size();
      RbfEval( ndim, ncpts, cpts.data(), cwts.data(), mRBFShape, npts, pts, vals, grds );
    };
}




void ndfes::FES::SetupBspl()
{
  mInterpMode = 3;

  if ( mIsMBAR )
    {
      std::cerr << "Cannot use B-spline interpolation on MBAR surface"
		<< std::endl;
      std::exit(EXIT_FAILURE);
    }
  
}


void ndfes::FES::SetupARBF( int const delta, bool const with_inverse )
{
  mInterpMode = 4;
  mARBFDelta = delta;
  mWithInverse = with_inverse;
}


void ndfes::FES::SetupWAvg( int const order, int const niter )
{
  if ( mGlbIdxMap.size() == 0 )
    {
      for ( std::size_t ibin=0; ibin < mNumBins; ++ibin )
	{
	  mGlbIdxMap.insert( {mBins[ibin].glbidx,ibin} );
	}
    }
  
  mInterpMode = 5;
  mOrder = order;
  mNBspl = mOrder + (mOrder%2);

  if ( niter > 0 )
    {
      std::vector<double> tvals;
      std::vector<std::size_t> occbins;
      for ( std::size_t ibin=0; ibin<mNumBins; ++ibin )
	{
	  if ( mBinSizes[ibin] > 0 )
	    {
	      occbins.push_back(ibin);
	      tvals.push_back(mBinVals[ibin]);
	    }
	}

      std::size_t const nocc = occbins.size();
      for ( int it=0; it<niter+1; ++it )
	{
	  std::vector<double> ovals(nocc,0);
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	  for ( std::size_t i=0; i<nocc; ++i )
	    {
	      std::size_t ibin = occbins[i];
	      std::vector<double> c( mBins[ibin].center );
	      double v = 0;
	      EvalWAvgValue( c.data(), v );
	      ovals[i] = v;
	    }

	  double maxerr = 0;
	  double avgerr = 0;
	  for ( std::size_t i=0; i<nocc; ++i )
	    {
	      double d = std::abs(tvals[i]-ovals[i]);
	      maxerr = std::max(maxerr,d);
	      avgerr += d;
	    }
	  avgerr /= nocc;

	  std::printf("wavg iter %4i  MAE %10.2e  MaxE %10.2e\n",
		      it, avgerr, maxerr);

	  if ( it < niter )
	    {
	      for ( std::size_t i=0; i<nocc; ++i )
		{
		  std::size_t ibin = occbins[i];
		  mBinVals[ibin] += tvals[i] - ovals[i];
		}
	    }
	}
      
    }
  
}



void ndfes::FES::EvalWAvgValue( double const * pt, double & ene ) const
{

  ene = 0;
  // now calculate the bin center values
  // we can calculate a single set of bspline weights

  std::size_t ndim = mDimInfo.GetNumDims();
  double const * widths = mDimInfo.GetTargetWidths();
  int const * sizes = mDimInfo.GetDimSizes();
  int const * isper = mDimInfo.GetIsPeriodic();
  double const * xmins = mDimInfo.GetXmin();
  //double const * xmaxs = mDimInfo.GetXmax();
  std::vector<std::size_t> bidxs( mDimInfo.CptBinIdxs( pt ) );
  //std::vector<double> bcent( mDimInfo.CptBinCenter( bidxs.data() ) );

  std::vector< std::vector<double> > linsubwts(ndim);
  std::vector<double> gwts(mOrder,0.);
  std::vector<int> gidx(mOrder,0);
  
  std::vector<double> cpt(ndim,0);
  std::vector<int> cidxs(ndim,0);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      double h = widths[dim]/2;
      double xlo = xmins[dim];
      double xp = xlo+h;
      int idelta = 0;
      if ( not isper[dim] )
	{
	  double frac = std::max(0.,(pt[dim]-xp)/widths[dim]);
	  idelta = frac;
	}
      else
	{
	  idelta = std::floor( (pt[dim]-xp)/widths[dim] );
	};

      cidxs[dim] = idelta;
      cpt[dim] = xp + (idelta+0.5)*widths[dim];
      // if ( pt[dim] - cpt[dim] > h )
      // 	{
      // 	  cpt[dim] += widths[dim];
      // 	}
      // else if ( pt[dim] - cpt[dim] < -h )
      // 	{
      // 	  cpt[dim] -= widths[dim];
      // 	}
      //std::printf("%12.5f %12.5f %12.5f\n",xlo,xp,pt[dim]);
      //std::printf("%3i %12.5f %12.5f %12.5f\n",idelta,cpt[dim],pt[dim]-cpt[dim],frac);
    }
  //std::printf("\n");
  
  std::vector< std::vector<std::size_t> > lidxs(ndim);

  
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      linsubwts[dim].assign(mOrder,0.);
      lidxs[dim].assign(mOrder,0);
      std::fill( gwts.data(), gwts.data() + mOrder, 0. );
      double w = (mNBspl-1)*widths[dim];
      double c = w/2.;
      double dx = c + (pt[dim]-cpt[dim]);
      //std::printf("mOrder %i/%i %i %i %12.6f %12.6f\n",isub,subbinsize,mOrder,mNBspl,dx,w);

      ndfes::bspline_aperiodic( dx, w, mNBspl, mOrder,
				gwts.data(), gidx.data() );


      // double dx = pt[dim]-(xmins[dim]+widths[dim]/2);
      // if ( isper[dim] )
      // 	{
      // 	  dx = wraptorange(dx,xmins[dim]+widths[dim]/2,xmaxs[dim]+widths[dim]/2);
      // 	}
      // dx = std::max(0.,dx);
      
      // ndfes::bspline_aperiodic( dx, xmaxs[dim]-xmins[dim], sizes[dim], mOrder,
      // 				gwts.data(), gidx.data() );

      for ( std::size_t ib=0; ib<mOrder; ++ib )
	{
	  int idx = (int)(cidxs[dim] + gidx[ib]) - (mNBspl/2-1);
	  //int idx = gidx[ib];
	  //std::printf("%12.5f %5i\n",wt,idx);
	  if ( isper[dim] )
	    {
	      idx = INTWRAP(idx,sizes[dim]);
	    }
	  else if ( idx < 0 or idx >= sizes[dim] )
	    {
	      std::cout << "INVALID IDX " << idx << " " << bidxs[dim] << " " << gidx[ib] << " IN FES.cpp ndfes::FES::EvalWAvgValue" << std::endl;
	      idx = bidxs[dim];
	      //wt = 0.;
	    }
	  //lidxs[dim][ gidx[ib] ] = idx;
	  //linsubwts[dim][ gidx[ib] ] = gwts[ib];
	  lidxs[dim][ib] = idx;
	  linsubwts[dim][ib] = gwts[ib];
	}
    }
  
  std::vector<double> gridsubwts;
  ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
  std::vector<std::size_t> midxs;
  ndfes::LinearSpacingsToMeshgrid( lidxs, midxs );

  std::size_t nm = midxs.size() / ndim;

  std::vector<bool> isocc(nm,false);
  std::vector<std::size_t> gidxs(nm,0);
  std::vector<std::size_t> locidxs(nm,0);
  double wsum = 0;

  for ( std::size_t i=0; i<nm; ++i )
    {
      gidxs[i] = mDimInfo.CptGlbBinIdx( midxs.data() + i*ndim );



      std::unordered_map<std::size_t,std::size_t>::const_iterator
	p = mGlbIdxMap.find(gidxs[i]);
      if ( p != mGlbIdxMap.end() )
	{
	  isocc[i] = true;
	  locidxs[i] = p->second;
	  wsum += gridsubwts[i];
	}
      
      
      // for ( std::size_t ibin=0; ibin<mNumBins; ++ibin )
      // 	{
      // 	  if ( mBins[ibin].glbidx == gidxs[i] )
      // 	    {
      // 	      isocc[i] = true;
      // 	      locidxs[i] = ibin;
      // 	      wsum += gridsubwts[i];
      // 	      break;
      // 	    }
      // 	}
    }

  
  for ( std::size_t i=0; i<nm; ++i )
    {
      if ( isocc[i] )
	{
	  double w = (gridsubwts[i]/wsum);
	  //double w = (gridsubwts[i]);
	  ene +=  w * mBinVals[locidxs[i]];
	}
    }
  
  // {
  //   std::printf("P  ");
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %6lu",0);
  //     }
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %12.5f",pt[dim]);
  //     }
  //   std::printf("\n");
    
  //   std::vector<double> c( mDimInfo.CptBinCenter( bidxs.data() ) );
  //   std::printf("C %lu",0);
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %6lu",bidxs[dim]);
  //     }
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %12.5f",c[dim]);
  //     }
  //   std::printf("\n");
    
  // }
  // for ( std::size_t i=0; i<nm; ++i )
  //   {
  //     std::vector<std::size_t> midx( midxs.data() + i*ndim, midxs.data() + (i+1)*ndim );
  //     std::vector<double> c( mDimInfo.CptBinCenter( midx.data() ) );
  //     std::printf("C %lu",i+1);
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf(" %6lu",midx[dim]);
  // 	}
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf(" %12.5f",c[dim]);
  // 	}
  //     std::printf(" %12.7f",gridsubwts[i]);
  //     std::printf("\n");
  //   }
  // std::printf("ene: %20.10e\n",ene);
  // std::exit(0);
}




void ndfes::FES::EvalWAvgValueAndError( double const * pt, double & ene, double & err ) const
{

  ene = 0;
  err = 0;
  // now calculate the bin center values
  // we can calculate a single set of bspline weights

  std::size_t ndim = mDimInfo.GetNumDims();
  double const * widths = mDimInfo.GetTargetWidths();
  int const * sizes = mDimInfo.GetDimSizes();
  int const * isper = mDimInfo.GetIsPeriodic();
  double const * xmins = mDimInfo.GetXmin();
  //double const * xmaxs = mDimInfo.GetXmax();
  std::vector<std::size_t> bidxs( mDimInfo.CptBinIdxs( pt ) );
  std::vector<double> bcent( mDimInfo.CptBinCenter( bidxs.data() ) );

  std::vector< std::vector<double> > linsubwts(ndim);
  std::vector<double> gwts(mOrder,0.);
  std::vector<int> gidx(mOrder,0);
  
  std::vector<double> cpt(ndim,0);
  //int shift = 0;
  std::vector<int> cidxs(ndim,0);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      double h = widths[dim]/2;
      double xlo = xmins[dim];
      double xp = xlo+h;
      int idelta = 0;
      if ( not isper[dim] )
	{
	  double frac = std::max(0.,(pt[dim]-xp)/widths[dim]);
	  idelta = frac;
	}
      else
	{
	  idelta = std::floor( (pt[dim]-xp)/widths[dim] );
	};
      // if ( isper[dim] )
      // 	{
      // 	  if ( pt[dim] < xp )
      // 	    {
      // 	      idelta = sizes[dim]-1;
      // 	    }
      // 	}
      cidxs[dim] = idelta;
      cpt[dim] = xp + (idelta+0.5)*widths[dim];


      // if ( std::abs(pt[dim]-cpt[dim]) > h )
      // 	{
      // 	  std::printf("Whoops\n");
      // 	}
      
      // if ( pt[dim] - cpt[dim] > h )
      // 	{
      // 	  cpt[dim] += widths[dim];
      // 	  shift += 1;
      // 	}
      // else if ( pt[dim] - cpt[dim] < -h )
      // 	{
      // 	  cpt[dim] -= widths[dim];
      // 	  shift -= 1;
      // 	}
      //std::printf("%3i %12.5f %12.5f %12.5f %2i\n",idelta,cpt[dim],pt[dim]-cpt[dim],frac,shift);
    }
  //std::printf("\n");
  
  std::vector< std::vector<std::size_t> > lidxs(ndim);

  
  // {
  //   std::printf("P %6lu",0);
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %6lu",0);
  //     }
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %12.5f",pt[dim]);
  //     }
  //   std::printf("\n");
    
  //   std::vector<double> c( mDimInfo.CptBinCenter( bidxs.data() ) );
  //   std::printf("C %6lu",0);
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %6lu",bidxs[dim]);
  //     }
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %12.5f",c[dim]);
  //     }
  //   std::printf("\n");
  // }

  for ( std::size_t dim=0; dim<ndim; ++dim )
    {

      
      linsubwts[dim].assign(mOrder,0.);
      lidxs[dim].assign(mOrder,0);
      std::fill( gwts.data(), gwts.data() + mOrder, 0. );
      double w = (mNBspl-1)*widths[dim];
      double c = w/2.;
      double dx = c + (pt[dim]-cpt[dim]);
      //double dx = c + (pt[dim]-bcent[dim]);

      // std::printf("dim: %3lu pt: %20.8f\n",dim,pt[dim]);
      // std::printf("mOrder %i %i %12.6f %12.6f\n",mOrder,mNBspl,dx,w);


      
      ndfes::bspline_aperiodic( dx, w, mNBspl, mOrder,
				gwts.data(), gidx.data() );

      // double dx = pt[dim]-(xmins[dim]+widths[dim]/2);
      // if ( isper[dim] )
      // 	{
      // 	  dx = wraptorange(dx,xmins[dim]+widths[dim]/2,xmaxs[dim]+widths[dim]/2);
      // 	}
      // dx = std::max(0.,dx);
      
      // ndfes::bspline_aperiodic( dx, xmaxs[dim]-xmins[dim], sizes[dim], mOrder,
      // 				gwts.data(), gidx.data() );

      //std::printf("%5i %5i\n",cidxs[dim],bidxs[dim]);
      
      for ( std::size_t ib=0; ib<mOrder; ++ib )
	{
	  int idx = (int)(cidxs[dim] + gidx[ib]) - (mNBspl/2-1) ;
	  //int idx = gidx[ib];
	  //std::printf("%2lu %5lu %5i %5i %5i\n",dim,ib,gidx[ib],idx,sizes[dim]);
	  if ( isper[dim] )
	    {
	      idx = INTWRAP(idx,sizes[dim]);
	    }
	  else if ( idx < 0 or idx >= sizes[dim] )
	    {
	      idx = bidxs[dim];
	      //wt = 0.;
	    }
	  //lidxs[dim][ gidx[ib] ] = idx;
	  //linsubwts[dim][ gidx[ib] ] = gwts[ib];
	  lidxs[dim][ib] = idx;
	  linsubwts[dim][ib] = gwts[ib];
	}
    }
  //std::cout << "OUT" << std::endl;

  // std::printf("pt");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf("%12.5f",pt[dim]);
  //   }
  // for ( std::size_t dim=0; dim<ndim; ++dim)
  //   {
  //     std::printf(" %5i",cidxs[dim]);
  //   }
  // std::printf("\n");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf("%1lu : ",dim);
  //     for ( std::size_t ii=0; ii<lidxs[dim].size(); ++ii )
  // 	{
  // 	  std::printf("%6lu",lidxs[dim][ii]);
  // 	}
  //     std::printf("\n");
  //   }
  // std::printf("\n");
  // std::printf("xmin:  ");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf("%12.5f",mDimInfo.GetXmin()[dim]);
  //   }
  // std::printf("\n");
  // std::printf("xmax:  ");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf("%12.5f",mDimInfo.GetXmax()[dim]);
  //   }
  // std::printf("\n");
  // std::printf("size:  ");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf("%6i",mDimInfo.GetDimSizes()[dim]);
  //   }
  // std::printf("\n");
  // std::printf("width: ");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf("%12.5f",mDimInfo.GetTargetWidths()[dim]);
  //   }
  // std::printf("\n");
  //std::exit(0);
  
  std::vector<double> gridsubwts;
  ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
  std::vector<std::size_t> midxs;
  ndfes::LinearSpacingsToMeshgrid( lidxs, midxs );

  std::size_t nm = midxs.size() / ndim;

  std::vector<bool> isocc(nm,false);
  std::vector<std::size_t> gidxs(nm,0);
  std::vector<std::size_t> locidxs(nm,0);
  double wsum = 0;

  for ( std::size_t i=0; i<nm; ++i )
    {
      gidxs[i] = mDimInfo.CptGlbBinIdx( midxs.data() + i*ndim );
      //std::cout << "midx ";
      //for ( std::size_t dim=0; dim<ndim; ++dim )
      //{
      //  std::cout << std::setw(6) << midxs[dim+i*ndim];
      //}
      //std::cout << " : " << gidxs[i] << std::endl;



      std::unordered_map<std::size_t,std::size_t>::const_iterator
	p = mGlbIdxMap.find(gidxs[i]);
      if ( p != mGlbIdxMap.end() )
	{
	  isocc[i] = true;
	  locidxs[i] = p->second;
	  wsum += gridsubwts[i];
	}
      
      
      // for ( std::size_t ibin=0; ibin<mNumBins; ++ibin )
      // 	{
      // 	  //std::printf("gidx search %6lu %6lu ? %6lu  : %i\n",
      // 	  //	      ibin,mBins[ibin].glbidx,gidxs[i],mBins[ibin].glbidx == gidxs[i]);
      // 	  if ( mBins[ibin].glbidx == gidxs[i] )
      // 	    {
      // 	      isocc[i] = true;
      // 	      locidxs[i] = ibin;
      // 	      wsum += gridsubwts[i];
      // 	      break;
      // 	    }
      // 	}
      if ( not isocc[i] )
	{
	  std::cerr << "FAILED TO FIND " << gidxs[i]
		    << " in FES.cpp ndfes::FES::EvalWAvgValueAndError"
		    << std::endl;
	  //std::exit(0);
	}
    }


  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     for ( std::size_t ib=0; ib<mOrder; ++ib )
  // 	{
  // 	  std::printf("%6lu",lidxs[dim][ib]);
  // 	}
  //     std::printf("\n");
  //   }
  
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     for ( std::size_t ib=0; ib<mOrder; ++ib )
  // 	{
  // 	  std::printf("%6.2f",linsubwts[dim][ib]);
  // 	}
  //     std::printf("\n");
  //   }
  
  
  
  for ( std::size_t i=0; i<nm; ++i )
    {
      if ( isocc[i] )
	{
	  double w = (gridsubwts[i]/wsum);
	  //double w = (gridsubwts[i]);
	  ene +=  w * mBinVals[locidxs[i]];
	  double t = w * mBinErrs[locidxs[i]];
	  err += t*t;
	}
    }
  err = std::sqrt(err);



  
}



void ndfes::FES::EvalWAvgValueAndGrad( double const * pt, double & ene, double * grd ) const
{

  
  // now calculate the bin center values
  // we can calculate a single set of bspline weights

  std::size_t const ndim = mDimInfo.GetNumDims();

  std::fill( grd, grd+ndim, 0. );
  ene = 0;
  
  double const * widths = mDimInfo.GetTargetWidths();
  int const * sizes = mDimInfo.GetDimSizes();
  int const * isper = mDimInfo.GetIsPeriodic();
  double const * xmins = mDimInfo.GetXmin();
  //double const * xmaxs = mDimInfo.GetXmax();
  std::vector<std::size_t> bidxs( mDimInfo.CptBinIdxs( pt ) );
  //std::vector<double> bcent( mDimInfo.CptBinCenter( bidxs.data() ) );

  std::vector< std::vector<double> > linsubwts(ndim);
  std::vector< std::vector<double> > linsubdwts(ndim);
  std::vector<double> gwts(mOrder,0.);
  std::vector<double> gdwts(mOrder,0.);
  std::vector<int> gidx(mOrder,0);

  std::vector<double> cpt(ndim,0);
  std::vector<int> cidxs(ndim,0);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      double h = widths[dim]/2;
      double xlo = xmins[dim];
      double xp = xlo+h;
      int idelta = 0;
      if ( not isper[dim] )
	{
	  double frac = std::max(0.,(pt[dim]-xp)/widths[dim]);
	  idelta = frac;
	}
      else
	{
	  idelta = std::floor( (pt[dim]-xp)/widths[dim] );
	};
      cidxs[dim] = idelta;
      cpt[dim] = xp + (idelta+0.5)*widths[dim];
      // if ( pt[dim] - cpt[dim] > h )
      // 	{
      // 	  cpt[dim] += widths[dim];
      // 	}
      // else if ( pt[dim] - cpt[dim] < -h )
      // 	{
      // 	  cpt[dim] -= widths[dim];
      // 	}
      //std::printf("%12.5f %12.5f %12.5f\n",xlo,xp,pt[dim]);
      //std::printf("%3i %12.5f %12.5f %12.5f\n",idelta,cpt[dim],pt[dim]-cpt[dim],frac);
    }
  //std::printf("\n");
  
  std::vector< std::vector<std::size_t> > lidxs(ndim);
  
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      linsubwts[dim].assign(mOrder,0.);
      linsubdwts[dim].assign(mOrder,0.);
      lidxs[dim].assign(mOrder,0);
      std::fill( gwts.data(), gwts.data() + mOrder, 0. );
      double w = (mNBspl-1)*widths[dim];
      double c = w/2.;
      double dx = c + (pt[dim]-cpt[dim]);
      //std::printf("mOrder %i/%i %i %i %12.6f %12.6f\n",isub,subbinsize,mOrder,mNBspl,dx,w);
      ndfes::bspline_aperiodic( dx, w, mNBspl, mOrder,
				gwts.data(), gdwts.data(), gidx.data() );

      // double dx = pt[dim]-(xmins[dim]+widths[dim]/2);
      // if ( isper[dim] )
      // 	{
      // 	  dx = wraptorange(dx,xmins[dim]+widths[dim]/2,xmaxs[dim]+widths[dim]/2);
      // 	}
      // dx = std::max(0.,dx);
      
      // ndfes::bspline_aperiodic( dx, xmaxs[dim]-xmins[dim], sizes[dim], mOrder,
      // 				gwts.data(), gdwts.data(), gidx.data() );

      for ( std::size_t ib=0; ib<mOrder; ++ib )
	{
	  int idx = (int)(cidxs[dim] + gidx[ib]) - (mNBspl/2-1);
	  //int idx = gidx[ib];
	  //std::printf("%12.5f %5i\n",wt,idx);
	  if ( isper[dim] )
	    {
	      idx = INTWRAP(idx,sizes[dim]);
	    }
	  else if ( idx < 0 or idx >= sizes[dim] )
	    {
	      std::cout << "INVALID IDX " << idx << " " << bidxs[dim] << " " << gidx[ib] << " in FES.cpp ndfes::FES::EvalWAvgValueAndGrad" << std::endl;
	      idx = bidxs[dim];
	      //wt = 0.;
	    }
	  //lidxs[dim][ gidx[ib] ] = idx;
	  //linsubwts[dim][ gidx[ib] ] = gwts[ib];
	  lidxs[dim][ib] = idx;
	  linsubwts[dim][ib] = gwts[ib];
	  linsubdwts[dim][ib] = gdwts[ib];
	}
    }

  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     for ( std::size_t ib=0; ib<mOrder; ++ib )
  // 	{
  // 	  std::printf("%6lu",lidxs[dim][ib]);
  // 	}
  //     std::printf("\n");
  //   }
  
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     for ( std::size_t ib=0; ib<mOrder; ++ib )
  // 	{
  // 	  std::printf("%6.2f",linsubwts[dim][ib]);
  // 	}
  //     std::printf("\n");
  //   }
  
      
  std::vector<double> gridsubwts;
  ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
  std::vector<std::size_t> midxs;
  ndfes::LinearSpacingsToMeshgrid( lidxs, midxs );

  std::size_t nm = midxs.size() / ndim;

  std::vector<bool> isocc(nm,false);
  std::vector<std::size_t> gidxs(nm,0);
  std::vector<std::size_t> locidxs(nm,0);
  double wsum = 0;

  
  for ( std::size_t i=0; i<nm; ++i )
    {
      gidxs[i] = mDimInfo.CptGlbBinIdx( midxs.data() + i*ndim );

      std::unordered_map<std::size_t,std::size_t>::const_iterator
	p = mGlbIdxMap.find(gidxs[i]);
      if ( p != mGlbIdxMap.end() )
	{
	  isocc[i] = true;
	  locidxs[i] = p->second;
	  wsum += gridsubwts[i];
	}
      
      // for ( std::size_t ibin=0; ibin<mNumBins; ++ibin )
      // 	{
      // 	  if ( mBins[ibin].glbidx == gidxs[i] )
      // 	    {
      // 	      isocc[i] = true;
      // 	      locidxs[i] = ibin;
      // 	      wsum += gridsubwts[i];
      // 	      break;
      // 	    }
      // 	}
    }

  
  for ( std::size_t i=0; i<nm; ++i )
    {
      if ( isocc[i] )
	{
	  double w = (gridsubwts[i]/wsum);
	  //double w = (gridsubwts[i]);
	  ene += w * mBinVals[locidxs[i]];
	}
    }


  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      std::vector<double> dwts;
      std::vector< std::vector<double> > lwts( linsubwts );
      lwts[dim] = linsubdwts[dim];
      ndfes::LinearWtsToMeshWts(lwts,dwts);

      double dsum = 0;
      
      for ( std::size_t i=0; i<nm; ++i )
	{
	  if ( isocc[i] )
	    {
	      dsum += dwts[i];
	    }
	}
      
      for ( std::size_t i=0; i<nm; ++i )
	{
	  if ( isocc[i] )
	    {
	      double t1 = dwts[i]/wsum;
	      double t2 = (-gridsubwts[i]/(wsum*wsum)) * dsum;
	      grd[dim] += (t1+t2) * mBinVals[locidxs[i]];
	      //grd[dim] += dwts[i] * mBinVals[locidxs[i]];
	    }
	}
    }
  
  // {
  //   std::printf("P %6lu",0);
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %6lu",0);
  //     }
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %12.5f",pt[dim]);
  //     }
  //   std::printf("\n");
    
  //   std::vector<double> c( mDimInfo.CptBinCenter( bidxs.data() ) );
  //   std::printf("C %6lu",0);
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %6lu",bidxs[dim]);
  //     }
  //   for ( std::size_t dim=0; dim<ndim; ++dim )
  //     {
  // 	std::printf(" %12.5f",c[dim]);
  //     }
  //   std::printf("\n");
    
  // }
  // for ( std::size_t i=0; i<nm; ++i )
  //   {
  //     std::vector<std::size_t> midx( midxs.data() + i*ndim, midxs.data() + (i+1)*ndim );
  //     std::vector<double> c( mDimInfo.CptBinCenter( midx.data() ) );
  //     std::printf("C %6lu",i+1);
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf(" %6lu",midx[dim]);
  // 	}
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	  {
  // 	    std::printf(" %12.5f",c[dim]);
  // 	  }
  //     std::printf(" %12.7f",gridsubwts[i]);
  //     int o = 0;
  //     if ( isocc[i] )
  // 	{
  // 	  o = 1;
  // 	}
  //     std::printf(" %i",o);
  //     std::printf("\n");
  //   }
  
}










double ndfes::FES::CptAvgOccAtPt( double const * pt, std::size_t const order ) const
{

  double ene = 0;
  std::size_t nbspl = order + (order%2);
  // now calculate the bin center values
  // we can calculate a single set of bspline weights

  std::size_t ndim = mDimInfo.GetNumDims();
  double const * widths = mDimInfo.GetTargetWidths();
  int const * sizes = mDimInfo.GetDimSizes();
  int const * isper = mDimInfo.GetIsPeriodic();
  double const * xmins = mDimInfo.GetXmin();
  std::vector<std::size_t> bidxs( mDimInfo.CptBinIdxs( pt ) );

  std::vector< std::vector<double> > linsubwts(ndim);
  std::vector<double> gwts(order,0.);
  std::vector<int> gidx(order,0);
  
  std::vector<double> cpt(ndim,0);
  std::vector<int> cidxs(ndim,0);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      double h = widths[dim]/2;
      double xlo = xmins[dim];
      double xp = xlo+h;
      int idelta = 0;
      if ( not isper[dim] )
	{
	  double frac = std::max(0.,(pt[dim]-xp)/widths[dim]);
	  idelta = frac;
	}
      else
	{
	  idelta = std::floor( (pt[dim]-xp)/widths[dim] );
	};

      cidxs[dim] = idelta;
      cpt[dim] = xp + (idelta+0.5)*widths[dim];
    }
  
  std::vector< std::vector<std::size_t> > lidxs(ndim);

  
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      linsubwts[dim].assign(order,0.);
      lidxs[dim].assign(order,0);
      std::fill( gwts.data(), gwts.data() + order, 0. );
      double w = (nbspl-1)*widths[dim];
      double c = w/2.;
      double dx = c + (pt[dim]-cpt[dim]);

      ndfes::bspline_aperiodic( dx, w, nbspl, order,
				gwts.data(), gidx.data() );


      for ( std::size_t ib=0; ib<order; ++ib )
	{
	  int idx = (int)(cidxs[dim] + gidx[ib]) - (nbspl/2-1);
	  if ( isper[dim] )
	    {
	      idx = INTWRAP(idx,sizes[dim]);
	    }
	  else if ( idx < 0 or idx >= sizes[dim] )
	    {
	      //std::cout << "INVALID IDX " << idx << " " << bidxs[dim] << " " << gidx[ib] << " IN FES.cpp ndfes::FES::EvalWAvgValue" << std::endl;
	      idx = bidxs[dim];
	      gwts[ib] = 0;
	      //wt = 0.;
	    }
	  //lidxs[dim][ gidx[ib] ] = idx;
	  //linsubwts[dim][ gidx[ib] ] = gwts[ib];
	  lidxs[dim][ib] = idx;
	  linsubwts[dim][ib] = gwts[ib];
	}
    }
  
  std::vector<double> gridsubwts;
  ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
  std::vector<std::size_t> midxs;
  ndfes::LinearSpacingsToMeshgrid( lidxs, midxs );

  std::size_t nm = midxs.size() / ndim;

  for ( std::size_t i=0; i<nm; ++i )
    {
      std::size_t gidx = mDimInfo.CptGlbBinIdx( midxs.data() + i*ndim );

      std::unordered_map<std::size_t,std::size_t>::const_iterator
	p = mGlbIdxMap.find(gidx);
      if ( p != mGlbIdxMap.end() )
	{
	  ene += gridsubwts[i] * mBinSizes[p->second];
	}
    }
  
  return ene;
}
















void ndfes::FES::PreparevFEP()
{
  
  typedef std::vector<std::size_t>::iterator it;

  mNBspl = mOrder + (mOrder%2);
  
  mNumCorners = mDimInfo.GetNumCorners();
  std::size_t nc = mNumCorners;

  mLocCidxs.assign( mNumCorners * mNumBins, 0 );
  for ( std::size_t ibin=0; ibin<mNumBins; ++ibin )
    {
      for ( std::size_t ic=0; ic<nc; ++ic )
	{
	  it p = std::find( mGlbCidxs.begin(), mGlbCidxs.end(),
			    mBins[ibin].cidxs[ic] );
	  if ( p != mGlbCidxs.end() )
	    {
	      mLocCidxs[ic+ibin*nc] = std::distance(mGlbCidxs.begin(),p);
	    }
	  else
	    {
	      std::cerr << "Failed to locate corner index "
                        << mBins[ibin].cidxs[ic]
                        << " from the global array" << std::endl;
              std::exit(EXIT_FAILURE);
	    }
	}
    }


  // now calculate the bin center values
  // we can calculate a single set of bspline weights

  int ndim = mDimInfo.GetNumDims();
  double const * widths = mDimInfo.GetTargetWidths();
  
  std::vector< std::vector<double> > linsubwts(ndim);
  std::vector<double> gwts(mOrder,0.);
  std::vector<int> gidx(mOrder,0);
  
  for ( int dim=0; dim<ndim; ++dim )
    {
      linsubwts[dim].assign(mNBspl,0.);
      std::fill( gwts.data(), gwts.data() + mOrder, 0. );
      double w = (mNBspl-1)*widths[dim];
      double c = w/2.;
      double dx = c;
      //std::printf("mOrder %i/%i %i %i %12.6f %12.6f\n",isub,subbinsize,mOrder,mNBspl,dx,w);
      ndfes::bspline_aperiodic( dx, w, mNBspl, mOrder,
				gwts.data(), gidx.data() );
      
      for ( std::size_t ib=0; ib<mOrder; ++ib )
	{
	  linsubwts[dim][ gidx[ib] ] = gwts[ib];
	}
    }
  
  std::vector<double> gridsubwts;
  ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);

  for ( std::size_t ibin=0; ibin<mNumBins; ++ibin )
    {
      double v = 0;
      double e2 = 0;
      for ( std::size_t ic=0; ic<nc; ++ic )
	{
	  std::size_t loc = mLocCidxs[ic+ibin*nc];
	  v += gridsubwts[ic] * mGlbCvals[loc];
	  double e = gridsubwts[ic] * mGlbCerrs[loc];
	  e2 += e*e;
	}
      mBinVals[ibin] = v;
      if ( e2 > 0 )
	{
	  mBinErrs[ibin] = std::sqrt(e2);
	}
      else
	{
	  mBinErrs[ibin] = 0;
	}
      mBinREs[ibin] = 1;
      //std::printf("%8lu %20.10e\n",mBins[ibin].glbidx,mBinVals[ibin]);
    }

  
}




void ndfes::FES::SetupQuad( int nquad )
{
  std::size_t ndim = mDimInfo.GetNumDims();
  mNumQuad = std::max(1,nquad);
  std::vector< std::vector<double> > linpts(ndim);
  std::vector< std::vector<double> > linwts(ndim);
  double const * pwidths = mDimInfo.GetTargetWidths();
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      linpts[dim].assign( mNumQuad, 0 );
      linwts[dim].assign( mNumQuad, 0 );
      GaussLegendreRule( -pwidths[dim]/2, pwidths[dim]/2,
			 mNumQuad,
			 linpts[dim].data(),
			 linwts[dim].data() );
    };

  mSubBinSize = std::pow(mNumQuad,ndim);
  ndfes::LinearSpacingsToMeshgrid(linpts,mQuadPts);
  ndfes::LinearWtsToMeshWts(linwts,mQuadWts);

  if ( mInterpMode == 3 )
    {
      mC2SC.assign( mNumCorners*mSubBinSize, 0. );
      std::size_t bsplorder = mOrder;
      std::size_t nbspl = mNBspl;
      
      double const * widths = mDimInfo.GetTargetWidths();
      for ( std::size_t isub=0; isub<mSubBinSize; ++isub )
	{
	  std::vector< std::vector<double> > linsubwts(ndim);
	  std::vector<double> gwts(bsplorder,0.);
	  std::vector<int> gidx(bsplorder,0);
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      linsubwts[dim].assign(nbspl,0.);
	      std::fill( gwts.data(), gwts.data() + bsplorder, 0. );
	      double w = (nbspl-1)*widths[dim];
	      double c = w/2.;
	      double dx = mQuadPts[dim+isub*ndim] + c;
	      //std::printf("bsplorder %i/%i %i %i %12.6f %12.6f\n",isub,subbinsize,bsplorder,nbspl,dx,w);
	      ndfes::bspline_aperiodic( dx, w, nbspl, bsplorder,
					gwts.data(), gidx.data() );
	      for ( std::size_t ib=0; ib<bsplorder; ++ib )
		{
		  linsubwts[dim][ gidx[ib] ] = gwts[ib];
		}
	    }
	  std::vector<double> gridsubwts;
	  ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
	  for ( std::size_t icorner=0, nc=mNumCorners; icorner<nc; ++icorner )
	    {
	      mC2SC[icorner+isub*nc] = gridsubwts[icorner];
	    };
	}
      
    }

  

  std::size_t nbins = mBins.size();
  mMeshFE.resize(nbins);
  if ( mInterpMode == 4 )
    {
      mARBFCpts.resize( nbins );
      mARBFCwts.resize( nbins );
      if ( mWithInverse )
	{
	  mARBFAinv.resize( nbins );
	  mARBFCidx.resize( nbins );
	}
    }
  
#ifdef WITH_OPENMP
#pragma omp parallel
#endif
  {
    std::size_t nqpts = mSubBinSize;
    std::vector<double> work( ndim * nqpts, 0 );
#ifdef WITH_OPENMP
#pragma omp for schedule(dynamic)
#endif
    for ( std::size_t ibin=0; ibin<nbins; ++ibin )
      {
	mMeshFE[ibin].assign( nqpts, 0 );
	EvalQuadInBin( ibin, work, mMeshFE[ibin].data() );  
      }
  }

  
}



void ndfes::FES::SetupNoQuad()
{
  std::size_t nbins = mBins.size();
  std::size_t ndim = mDimInfo.GetNumDims();
  // if ( mInterpMode == 3 )
  //   {
  //     PreparevFEP();
  //   }
  // else
  if ( mInterpMode == 4 )
    {
      mARBFCpts.resize( nbins );
      mARBFCwts.resize( nbins );
      if ( mWithInverse )
	{
	  mARBFAinv.resize( nbins );
	  mARBFCidx.resize( nbins );
	}

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for ( std::size_t ibin=0; ibin<nbins; ++ibin )
	{
	  //int didx = mARBFDelta;

      
	  //int nd = 2*didx+1;
	  //std::size_t maxnbor = std::pow( nd, (int)ndim );
      
      
      
	  ////////////////////////////////////
	  std::vector<std::size_t> nbors;
	  std::vector<double> cwts;

	  { // get neighbor indexes and values
      
	    // std::vector< std::vector<std::size_t> > lidxs(ndim);
	    // std::size_t nmidxs = 1;
	    // for ( std::size_t dim=0; dim<ndim; ++dim )
	    //   {
	    // 	int dimsize =  mDimInfo.GetDimSizes()[dim];
	    // 	int bidx = mBins[ibin].bidxs[dim];
	    // 	if ( mDimInfo.IsPeriodic(dim) )
	    // 	  {
	    // 	    int idel = -didx;
	    // 	    for ( int i=0; i<nd; ++i, ++idel )
	    // 	      {
	    // 		int k = bidx + idel;
	    // 		lidxs[dim].push_back( INTWRAP( k, dimsize ) );
	    // 	      }
	    // 	  }
	    // 	else
	    // 	  {
	    // 	    int idel = -didx;
	    // 	    for ( int i=0; i<nd; ++i, ++idel )
	    // 	      {
	    // 		int k = bidx + idel;
	    // 		if ( k >= 0 and k < dimsize )
	    // 		  {
	    // 		    lidxs[dim].push_back( k );
	    // 		  };
	    // 	      }
	    // 	  }
	    // 	nmidxs *= lidxs[dim].size();
	    //   }
	    // std::vector<std::size_t> midxs;
	    // ndfes::LinearSpacingsToMeshgrid( lidxs, midxs );

	    std::vector<std::size_t> midxs
	      ( ndfes::GetMeshIdxs( mDimInfo,mARBFDelta,
				    mBins[ibin].bidxs.data(),false) );
	    std::size_t nmidxs = midxs.size()/ndim;
	    
	    // if ( nmidxs != midxs.size()/ndim )
	    //   {
	    // 	std::cerr << "nmidx size mismatch " << nmidxs
	    // 		  << " " << midxs.size()/ndim << std::endl;
	    // 	std::exit(EXIT_FAILURE);
	    //   }
	    if ( midxs.size() % ndim != 0 )
	      {
		std::cerr << "midxs not a multiple of " << ndim
			  << " " << midxs.size() << std::endl;
		std::exit(EXIT_FAILURE);
	      }
	
	    for ( std::size_t i=0; i<nmidxs; ++i )
	      {
		std::size_t gidx = mDimInfo.CptGlbBinIdx( midxs.data() + i*ndim );
	    
		std::size_t idx = mBins.size();
		for ( std::size_t k=0; k<mBins.size(); ++k )
		  {
		    if ( mBins[k].glbidx == gidx )
		      {
			idx = k;
			break;
		      }
		  }
	    
		if ( idx != mBins.size() )
		  {
		    nbors.push_back( idx );
		    cwts.push_back( mBinVals[ idx ] );
		  }
	    
	      }

	    //std::printf("%5lu %5lu\n",nbors.size(),maxnbor);
	  } // get neighbor indexes and values


      
	  if ( true )
	    { // orig
      
	      ////////////////////////////////////
	      std::size_t ncpts = nbors.size();
	      std::vector<double> cpts( ndim*ncpts, 0 );
	      { // calc cpts
		std::vector<double> c0( mBins[ibin].center );
		std::vector<double> halfwidths( ndim, 0 );
		for ( std::size_t dim=0; dim<ndim; ++dim )
		  {
		    halfwidths[dim] = 0.5 * ( mDimInfo.GetXmax()[dim] - mDimInfo.GetXmin()[dim] );
		  }
		for ( std::size_t i=0; i<ncpts; ++i )
		  {
		    for ( std::size_t dim=0; dim<ndim; ++dim )
		      {
			if ( mDimInfo.IsPeriodic(dim) )
			  {
			    double dx = mBins[ nbors[i] ].center[dim] - c0[dim];
			    dx = wraptorange( dx, -halfwidths[dim], halfwidths[dim] );
			    cpts[dim+i*ndim] = c0[dim] + dx;
			  }
			else
			  {
			    cpts[dim+i*ndim] = mBins[ nbors[i] ].center[dim];
			  }
		      }
		  }
	      } // calc cpts
	

	      ////////////////////////////////////
	      { // calc A
		std::vector<double> A( ncpts*ncpts, 0. );
		for ( std::size_t i=0; i<ncpts; ++i )
		  {
		    A[i+i*ncpts] = 1;
		    //cwts[i] = mBinVals[ nbors[i] ];
	      
		    double const * ci = cpts.data() + i*ndim;
		    for ( std::size_t j=0; j<i; ++j )
		      {
			double const * cj = cpts.data() + j*ndim;
		  
			double r2 = 0;
			for ( std::size_t k=0; k<ndim; ++k )
			  {
			    double dx = ci[k]-cj[k];
			    r2 += dx*dx;
			  }
			double f = std::sqrt( 1 + mRBFShape * r2 );
			A[j+i*ncpts] = f;
			A[i+j*ncpts] = f;
		      }
		  }

		if ( mWithInverse )
		  {
		    RbfSolveWithInverse( ncpts, A.data(), cwts.data() );
		    mARBFAinv[ibin] = A;
		    mARBFCidx[ibin] = nbors;
		  }
		else
		  {
		    RbfSolve( ncpts, A.data(), cwts.data() );
		  }
	      } // calc A
	
	      mARBFCpts[ibin] = cpts;
	      mARBFCwts[ibin] = cwts;
	    } // orig

	  
	}

    }

  if ( mGlbIdxMap.size() == 0 )
    {
      for ( std::size_t ibin=0; ibin < mNumBins; ++ibin )
	{
	  mGlbIdxMap.insert( {mBins[ibin].glbidx,ibin} );
	}
    }
  
}




void ndfes::FES::EvalQuadInBin
( std::size_t ibin,
  std::vector<double> & work,
  double * vals )
{
  std::size_t ndim = mDimInfo.GetNumDims();
  
  if ( mInterpMode == 1 )
    {
      std::fill( vals, vals+mSubBinSize, mBinVals[ibin] );
    }
  else if ( mInterpMode == 2 )
    {
      if ( work.size() < ndim * mSubBinSize )
	{
	  work.resize( ndim * mSubBinSize );
	}
      double const * c = mBins[ibin].center.data();
      for ( std::size_t i=0; i<mSubBinSize; ++i )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      work[dim+i*ndim] = mQuadPts[dim+i*ndim] + c[dim];
	    }
	}
      EvalRBF( mSubBinSize, work.data(), vals );
    }
  else if ( mInterpMode == 3 )
    {
      // CornerValsToSubCellVals
      if ( work.size() < mNumCorners )
	{
	  work.resize(mNumCorners);
	}
      for ( std::size_t ic=0; ic<mNumCorners; ++ic )
	{
	  work[ic] = mGlbCvals[ mLocCidxs[ic+ibin*mNumCorners] ];
	}
      double alp = 1.;
      double bet = 0.;
      int inc = 1;
      int ncorner = mNumCorners;
      int subbinsize = mSubBinSize;
      dgemv_( "T", &ncorner, &subbinsize, &alp, mC2SC.data(),
	      &ncorner, work.data(), &inc, &bet, vals, &inc );
    }
  else if ( mInterpMode == 4 )
    {
      //int didx = mARBFDelta;

      
      //int nd = 2*didx+1;
      //std::size_t maxnbor = std::pow( nd, (int)ndim );
      
      
      
      ////////////////////////////////////
      std::vector<std::size_t> nbors;
      std::vector<double> cwts;

      { // get neighbor indexes and values
      
	// std::vector< std::vector<std::size_t> > lidxs(ndim);
	// std::size_t nmidxs = 1;
	// for ( std::size_t dim=0; dim<ndim; ++dim )
	//   {
	//     int dimsize =  mDimInfo.GetDimSizes()[dim];
	//     int bidx = mBins[ibin].bidxs[dim];
	//     if ( mDimInfo.IsPeriodic(dim) )
	//       {
	// 	int idel = -didx;
	// 	for ( int i=0; i<nd; ++i, ++idel )
	// 	  {
	// 	    int k = bidx + idel;
	// 	    lidxs[dim].push_back( INTWRAP( k, dimsize ) );
	// 	  }
	//       }
	//     else
	//       {
	// 	int idel = -didx;
	// 	for ( int i=0; i<nd; ++i, ++idel )
	// 	  {
	// 	    int k = bidx + idel;
	// 	    if ( k >= 0 and k < dimsize )
	// 	      {
	// 		lidxs[dim].push_back( k );
	// 	      };
	// 	  }
	//       }
	//     nmidxs *= lidxs[dim].size();
	//   }
	// std::vector<std::size_t> midxs;
	// ndfes::LinearSpacingsToMeshgrid( lidxs, midxs );

	// if ( nmidxs != midxs.size()/ndim )
	//   {
	//     std::cerr << "nmidx size mismatch " << nmidxs
	// 	      << " " << midxs.size()/ndim << std::endl;
	//     std::exit(EXIT_FAILURE);
	//   }

	std::vector<std::size_t> midxs
	  ( ndfes::GetMeshIdxs( mDimInfo, mARBFDelta,
				mBins[ibin].bidxs.data(), false ) );
	std::size_t nmidxs = midxs.size()/ndim;
	
	if ( midxs.size() % ndim != 0 )
	  {
	    std::cerr << "midxs not a multiple of " << ndim
		      << " " << midxs.size() << std::endl;
	    std::exit(EXIT_FAILURE);
	  }
	
	for ( std::size_t i=0; i<nmidxs; ++i )
	  {
	    std::size_t gidx = mDimInfo.CptGlbBinIdx( midxs.data() + i*ndim );
	    
	    std::size_t idx = mBins.size();
	    for ( std::size_t k=0; k<mBins.size(); ++k )
	      {
		if ( mBins[k].glbidx == gidx )
		  {
		    idx = k;
		    break;
		  }
	      }
	    
	    if ( idx != mBins.size() )
	      {
		nbors.push_back( idx );
		cwts.push_back( mBinVals[ idx ] );
	      }
	    
	  }

	//std::printf("%5lu %5lu\n",nbors.size(),maxnbor);
      } // get neighbor indexes and values


      
      if ( true )
      { // orig
      
	////////////////////////////////////
	std::size_t ncpts = nbors.size();
	std::vector<double> cpts( ndim*ncpts, 0 );
	{ // calc cpts
	  std::vector<double> c0( mBins[ibin].center );
	  std::vector<double> halfwidths( ndim, 0 );
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      halfwidths[dim] = 0.5 * ( mDimInfo.GetXmax()[dim] - mDimInfo.GetXmin()[dim] );
	    }
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  if ( mDimInfo.IsPeriodic(dim) )
		    {
		      double dx = mBins[ nbors[i] ].center[dim] - c0[dim];
		      dx = wraptorange( dx, -halfwidths[dim], halfwidths[dim] );
		      cpts[dim+i*ndim] = c0[dim] + dx;
		    }
		  else
		    {
		      cpts[dim+i*ndim] = mBins[ nbors[i] ].center[dim];
		    }
		}
	    }
	} // calc cpts
	

	////////////////////////////////////
	{ // calc A
	  std::vector<double> A( ncpts*ncpts, 0. );
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      A[i+i*ncpts] = 1;
	      //cwts[i] = mBinVals[ nbors[i] ];
	      
	      double const * ci = cpts.data() + i*ndim;
	      for ( std::size_t j=0; j<i; ++j )
		{
		  double const * cj = cpts.data() + j*ndim;
		  
		  double r2 = 0;
		  for ( std::size_t k=0; k<ndim; ++k )
		    {
		      double dx = ci[k]-cj[k];
		      r2 += dx*dx;
		    }
		  double f = std::sqrt( 1 + mRBFShape * r2 );
		  A[j+i*ncpts] = f;
		  A[i+j*ncpts] = f;
		}
	    }

	  if ( mWithInverse )
	    {
	      RbfSolveWithInverse( ncpts, A.data(), cwts.data() );
	      mARBFAinv[ibin] = A;
	      mARBFCidx[ibin] = nbors;
	    }
	  else
	    {
	      RbfSolve( ncpts, A.data(), cwts.data() );
	    }
	} // calc A
	

	if ( work.size() < ndim * mSubBinSize )
	  {
	    work.resize( ndim * mSubBinSize );
	  }
	double const * c = mBins[ibin].center.data();
	for ( std::size_t i=0; i<mSubBinSize; ++i )
	  {
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		work[dim+i*ndim] = mQuadPts[dim+i*ndim] + c[dim];
	      }
	  }
	RbfEval( ndim, ncpts, cpts.data(), cwts.data(), mRBFShape, mSubBinSize, work.data(), vals );
	mARBFCpts[ibin] = cpts;
	mARBFCwts[ibin] = cwts;
      } // orig


      
    }
}



void ndfes::FES::BootstrapFES( double const booterror )
{
  std::random_device rd;
  std::mt19937 gen( rd() );
  
  for ( std::size_t ibin=0; ibin<mNumBins; ++ibin )
    {
      if ( booterror > 0 )
	{
	  std::normal_distribution<double> dist(mBinVals[ibin],booterror);
	  mBinVals[ibin] = dist(gen);
	}
      else
	{
	  std::normal_distribution<double> dist(mBinVals[ibin],mBinErrs[ibin]);
	  mBinVals[ibin] = dist(gen);
	}
    }
}



bool ndfes::FES::GetValue( double const * pt, double & val ) const
{
  val = 5000;
  std::vector<std::size_t> bidx( mDimInfo.CptBinIdxs( pt ) );
  std::size_t gidx = mDimInfo.CptGlbBinIdx( bidx.data() );
  std::size_t ndim = mDimInfo.GetNumDims();

  bool found = false;
  std::size_t ibin = 0;
  // for ( std::size_t i=0; i < mNumBins; ++i )
  //   {
  //     if ( mBins[i].glbidx == gidx )
  // 	{
  // 	  ibin = i;
  // 	  found = true;
  // 	  break;
  // 	}
  //   }
  {
    std::unordered_map<std::size_t,std::size_t>::const_iterator
      p = mGlbIdxMap.find(gidx);
    if ( p != mGlbIdxMap.end() )
      {
	ibin = p->second;
	found = true;
      }
  }


  

  if ( found )
    {
      if ( mInterpMode == 1 )
	{
	  val = mBinVals[ibin];
	}
      else if ( mInterpMode == 2 )
	{
	  EvalRBF( 1, pt, &val );
	}
      else if ( mInterpMode == 3 )
	{
	  std::vector<double> work( mNumCorners, 0 );
	  for ( std::size_t ic=0; ic<mNumCorners; ++ic )
	    {
	      work[ic] = mGlbCvals[ mLocCidxs[ic+ibin*mNumCorners] ];
	    }
	  
	  std::size_t bsplorder = mOrder;
	  std::size_t nbspl = mNBspl;
	  
	  double const * widths = mDimInfo.GetTargetWidths();
	  //for ( std::size_t isub=0; isub<mSubBinSize; ++isub )
	    {
	      std::vector< std::vector<double> > linsubwts(ndim);
	      std::vector<double> gwts(bsplorder,0.);
	      std::vector<int> gidx(bsplorder,0);
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  linsubwts[dim].assign(nbspl,0.);
		  std::fill( gwts.data(), gwts.data() + bsplorder, 0. );
		  double w = (nbspl-1)*widths[dim];
		  double c = w/2.;
		  double dx = (pt[dim] - mBins[ibin].center[dim]) + c;
		  ndfes::bspline_aperiodic( dx, w, nbspl, bsplorder,
					    gwts.data(), gidx.data() );
		  for ( std::size_t ib=0; ib<bsplorder; ++ib )
		    {
		      linsubwts[dim][ gidx[ib] ] = gwts[ib];
		    }
		}
	      std::vector<double> gridsubwts;
	      ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
	      
	      val = 0;
	      for ( std::size_t i=0; i<mNumCorners; ++i )
		{
		  val += gridsubwts[i] * work[i];
		}
	      
	    }
	}
      else if ( mInterpMode == 4 )
	{
	  std::size_t ncpts = mARBFCwts[ibin].size();
	  RbfEval( ndim, ncpts, mARBFCpts[ibin].data(), mARBFCwts[ibin].data(),
		   mRBFShape, 1, pt, &val );
	}
      else if ( mInterpMode == 5 )
	{
	  EvalWAvgValue(pt,val);
	}

    }
  return found;
}

bool ndfes::FES::GetValueErrorAndEntropy( double const * pt, double & val, double & err, double & re, std::size_t & ns ) const
{
  val = 5000;
  err = 0;
  re = 0;
  ns = 0;
  
  std::size_t ndim = mDimInfo.GetNumDims();
  
  std::vector<std::size_t> bidx( mDimInfo.CptBinIdxs( pt ) );

  // std::printf("pt");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf(" %20.10e",pt[dim]);
  //   }
  
  // std::printf(" bidx");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf(" %6lu",bidx[dim]);
  //   }
  // std::printf("\n");
  
  std::size_t gidx = mDimInfo.CptGlbBinIdx( bidx.data() );

  bool found = false;
  std::size_t ibin = 0;
  // for ( std::size_t i=0; i < mNumBins; ++i )
  //   {
  //     if ( mBins[i].glbidx == gidx )
  // 	{
  // 	  ibin = i;
  // 	  found = true;
  // 	  break;
  // 	}
  //   }

  {
    std::unordered_map<std::size_t,std::size_t>::const_iterator
      p = mGlbIdxMap.find(gidx);
    if ( p != mGlbIdxMap.end() )
      {
	ibin = p->second;
	found = true;
      }
  }
  
  if ( found )
    {
      //std::cout << "FOUND\n";
      // std::printf("pt ");
      // for ( std::size_t dim=0; dim<ndim; ++dim )
      // 	{
      // 	  std::printf("%9.5f",pt[dim]);
      // 	}
      // std::printf(" bidx ");
      // for ( std::size_t dim=0; dim<ndim; ++dim )
      // 	{
      // 	  std::printf("%6lu",bidx[dim]);
      // 	}
      // std::printf(" gidx ");
      // std::printf("%6lu",gidx);
      // std::printf("%6lu",mBinSizes[ibin]);
      // std::printf("\n");
      
      re = mBinREs[ibin];
      ns = mBinSizes[ibin];
      if ( mInterpMode == 1 )
	{
	  val = mBinVals[ibin];
	  err = mBinErrs[ibin];
	}
      else if ( mInterpMode == 2 )
	{
	  std::size_t ncpts = mRBFCwts.size();
	  std::vector<double> A(ncpts,0);
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      double r2 = 0;
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  double dx = pt[dim] - mRBFCpts[dim+i*ndim];
		  r2 += dx*dx;
		};
	      A[i] = std::sqrt( 1 + mRBFShape * r2 );
	    }
	  val = 0;
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      val += A[i] * mRBFCwts[i];
	    }
	  std::vector<double> T(ncpts,0);
	  double alpha=1;
	  double beta=0;
	  int inc=1;
	  int n = ncpts;
	  dgemv_( "T", &n, &n, &alpha, mRBFAinv.data(), &n,
		  A.data(),  &inc, &beta, T.data(), &inc );
	  err = 0;
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      double dx = T[i] * mRBFCerr[i];
	      err += dx*dx;
	    }
	  err = std::sqrt(err);
	}
      else if ( mInterpMode == 3 )
	{
	  std::vector<double> work( mNumCorners, 0 );
	  std::vector<double> dwork( mNumCorners, 0 );
	  for ( std::size_t ic=0; ic<mNumCorners; ++ic )
	    {
	      // std::size_t iloc = ic+ibin*mNumCorners;
	      // std::cout << "ic " << std::setw(6) << ic
	      // 		<< " ibin " << std::setw(6) << ibin
	      // 		<< " ncor " << std::setw(6) << mNumCorners
	      // 		<< " iloc " << std::setw(6) << iloc
	      // 		<< " lsize " << std::setw(6) << mLocCidxs.size()
	      // 		<< " gsize " << std::setw(6) << mGlbCvals.size()
	      // 		<< std::endl;
	      // std::cout << "gcidx " << std::setw(6) << mLocCidxs[iloc]
	      // 		<< std::endl;
	      work[ic]  = mGlbCvals[ mLocCidxs[ic+ibin*mNumCorners] ];
	      dwork[ic] = mGlbCerrs[ mLocCidxs[ic+ibin*mNumCorners] ];
	    }
	  
	  std::size_t bsplorder = mOrder;
	  std::size_t nbspl = mNBspl;
	  
	  double const * widths = mDimInfo.GetTargetWidths();
	  //for ( std::size_t isub=0; isub<mSubBinSize; ++isub )
	    {
	      std::vector< std::vector<double> > linsubwts(ndim);
	      std::vector<double> gwts(bsplorder,0.);
	      std::vector<int> gidx(bsplorder,0);
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  linsubwts[dim].assign(nbspl,0.);
		  std::fill( gwts.data(), gwts.data() + bsplorder, 0. );
		  double w = (nbspl-1)*widths[dim];
		  double c = w/2.;
		  double dx = (pt[dim] - mBins[ibin].center[dim]) + c;
		  ndfes::bspline_aperiodic( dx, w, nbspl, bsplorder,
					    gwts.data(), gidx.data() );
		  for ( std::size_t ib=0; ib<bsplorder; ++ib )
		    {
		      linsubwts[dim][ gidx[ib] ] = gwts[ib];
		    }
		}
	      std::vector<double> gridsubwts;
	      ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
	      
	      val = 0;
	      for ( std::size_t i=0; i<mNumCorners; ++i )
		{
		  val += gridsubwts[i] * work[i];
		}

	      err = 0;
	      for ( std::size_t i=0; i<mNumCorners; ++i )
		{
		  double dx = gridsubwts[i] * dwork[i];
		  err += dx*dx;
		}
	      err = std::sqrt(err);
	    }
	}
      else if ( mInterpMode == 4 )
	{
	  std::size_t ncpts = mARBFCwts[ibin].size();
	  std::vector<double> A(ncpts,0);
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      double r2 = 0;
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  double dx = pt[dim] - mARBFCpts[ibin][dim+i*ndim];
		  r2 += dx*dx;
		};
	      A[i] = std::sqrt( 1 + mRBFShape * r2 );
	    }
	  val = 0;
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      val += A[i] * mARBFCwts[ibin][i];
	    }
	  std::vector<double> T(ncpts,0);
	  double alpha=1;
	  double beta=0;
	  int inc=1;
	  int n = ncpts;
	  dgemv_( "T", &n, &n, &alpha, mARBFAinv[ibin].data(), &n,
		  A.data(),  &inc, &beta, T.data(), &inc );
	  err = 0;
	  for ( std::size_t i=0; i<ncpts; ++i )
	    {
	      double dx = T[i] * mBinErrs[ mARBFCidx[ibin][i] ];
	      err += dx*dx;
	    }
	  err = std::sqrt(err);
	}
      else if ( mInterpMode == 5 )
	{
	  if ( ns > 0 )
	    {
	      EvalWAvgValueAndError(pt,val,err);
	    };
	}
    }
  return found;
}



std::vector<double> ndfes::FES::GetModifiedBinValues() const
{
  return GetModifiedBinValues( mBinVals );
}


std::vector<double> ndfes::FES::GetModifiedBinValues( std::vector<double> const & binvals ) const
{
  std::vector<double> vs( binvals );

  double const minre = 1.e-20;
  double const maxpen = -std::log(minre);
  for ( std::size_t ibin=0; ibin<mNumBins; ++ibin )
    {
      double re = 1;
      if ( mIsMBAR )
	{
	  re = SwitchOn(mBinREs[ibin],0.3,0.5);
	}
      double penfact = 0;
      if ( re <= minre )
	{
	  penfact = 1;
	}
      else
	{
	  penfact = -std::log(re) / maxpen;
	}

      if ( mBinSizes[ibin] < 25 )
	{
	  penfact += 1. - (mBinSizes[ibin]/25.);
	}

      penfact = std::min(1.,std::max(0.,penfact));
      vs[ibin] += penfact * maxpen;
      
      // p = re * exp(-beta*U)
      // ln(p) = ln(re) - beta*U
      // ln(p) = - beta*(U-ln(re))
      

    }

  return vs;
}











bool ndfes::FES::GetValue( double const * inppt, double & val, double * grd ) const
{
  val = 5000;
  std::size_t const ndim = mDimInfo.GetNumDims();
  std::vector<double> pt( inppt, inppt + ndim );
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      if ( mDimInfo.IsPeriodic(dim) )
	{
	  double xmin = mDimInfo.GetXmin()[dim];
	  double xmax = mDimInfo.GetXmax()[dim];
	  pt[dim] = wraptorange(pt[dim],xmin,xmax);
	}
    }
  std::vector<std::size_t> bidx( mDimInfo.CptBinIdxs( pt.data() ) );
  std::size_t gidx = mDimInfo.CptGlbBinIdx( bidx.data() );


  std::fill( grd, grd+ndim, 0. );
  
  bool found = false;
  std::size_t ibin = 0;
  // for ( std::size_t i=0; i < mNumBins; ++i )
  //   {
  //     if ( mBins[i].glbidx == gidx )
  // 	{
  // 	  ibin = i;
  // 	  found = true;
  // 	  break;
  // 	}
  //   }
  {
    std::unordered_map<std::size_t,std::size_t>::const_iterator
      p = mGlbIdxMap.find(gidx);
    if ( p != mGlbIdxMap.end() )
      {
	ibin = p->second;
	found = true;
      }
  }

  if ( found )
    {
      val = 0;
      if ( mInterpMode == 1 )
	{
	  val = mBinVals[ibin];
	}
      else if ( mInterpMode == 2 )
	{
	  EvalRBF( 1, pt.data(), &val, grd );
	}
      else if ( mInterpMode == 3 )
	{
	  std::vector<double> work( mNumCorners, 0 );
	  for ( std::size_t ic=0; ic<mNumCorners; ++ic )
	    {
	      work[ic] = mGlbCvals[ mLocCidxs[ic+ibin*mNumCorners] ];
	    }
	  
	  std::size_t bsplorder = mOrder;
	  std::size_t nbspl = mNBspl;
	  
	  double const * widths = mDimInfo.GetTargetWidths();
	  //for ( std::size_t isub=0; isub<mSubBinSize; ++isub )
	    {
	      std::vector< std::vector<double> > linsubwts(ndim);
	      std::vector< std::vector<double> > linsubdwts(ndim);
	      std::vector<double> gwts(bsplorder,0.);
	      std::vector<double> gdwts(bsplorder,0.);
	      std::vector<int> gidx(bsplorder,0);
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  linsubwts[dim].assign(nbspl,0.);
		  linsubdwts[dim].assign(nbspl,0.);
		  std::fill( gwts.data(), gwts.data() + bsplorder, 0. );
		  double w = (nbspl-1)*widths[dim];
		  double c = w/2.;
		  double dx = (pt[dim] - mBins[ibin].center[dim]) + c;
		  
		  ndfes::bspline_aperiodic( dx, w, nbspl, bsplorder,
					    gwts.data(), gdwts.data(), gidx.data() );

		  // {
		  //   double const DEL = 5.e-5;
		  //   std::vector<double> whi(bsplorder,0.);
		  //   std::vector<double> wlo(bsplorder,0.);
		  //   std::vector<double> dw(bsplorder,0.);
		  //   ndfes::bspline_aperiodic( dx+DEL, w, nbspl, bsplorder,
		  // 			      whi.data(), gidx.data() );
		  //   ndfes::bspline_aperiodic( dx-DEL, w, nbspl, bsplorder,
		  // 			      wlo.data(), gidx.data() );
		    
		  // }
		  
		  for ( std::size_t ib=0; ib<bsplorder; ++ib )
		    {
		      linsubwts[dim][ gidx[ib] ] = gwts[ib];
		      linsubdwts[dim][ gidx[ib] ] = gdwts[ib];
		    }
		}
	      std::vector<double> gridsubwts;
	      ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
	      
	      val = 0;
	      for ( std::size_t i=0; i<mNumCorners; ++i )
		{
		  val += gridsubwts[i] * work[i];
		}

	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  std::vector< std::vector<double> > tmplin( linsubwts );
		  tmplin[dim] = linsubdwts[dim];
		  ndfes::LinearWtsToMeshWts(tmplin,gridsubwts);
		  for ( std::size_t i=0; i<mNumCorners; ++i )
		    {
		      grd[dim] += gridsubwts[i] * work[i];
		    }
		}
	      
	    }
	}
      else if ( mInterpMode == 4 )
	{
	  std::size_t ncpts = mARBFCwts[ibin].size();
	  RbfEval( ndim, ncpts, mARBFCpts[ibin].data(), mARBFCwts[ibin].data(),
		   mRBFShape, 1, pt.data(), &val, grd );
	}
      else if ( mInterpMode == 5 )
	{
	  EvalWAvgValueAndGrad(pt.data(),val,grd);
	}
    }
  return found;
}





bool ndfes::FES::GetPenalizedValue
( double const oobk,
  double const * inppt,
  double & val,
  double * grd ) const
{
  std::size_t const ndim = mDimInfo.GetNumDims();

  val = 5000;
  std::fill( grd, grd+ndim, 0. );
  double pen = 0;
  std::vector<double> dpen(ndim,0.);
  
  std::vector<double> pt( inppt, inppt + ndim );
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      if ( mDimInfo.IsPeriodic(dim) )
	{
	  double xmin = mDimInfo.GetXmin()[dim];
	  double xmax = mDimInfo.GetXmax()[dim];
	  pt[dim] = wraptorange(pt[dim],xmin,xmax);
	}
    }
  std::vector<std::size_t> bidx( mDimInfo.CptBinIdxs( pt.data() ) );
  std::size_t gidx = mDimInfo.CptGlbBinIdx( bidx.data() );


  bool found = false;
  std::size_t ibin = 0;

  //std::cout << "GetPenValue" << std::endl;
  // for ( std::size_t i=0; i < mNumBins; ++i )
  //   {
  //     //std::printf("%6lu %6lu %6lu %i\n",i,mBins[i].glbidx,gidx,mBins[i].glbidx == gidx);
  //     if ( mBins[i].glbidx == gidx )
  // 	{
  // 	  ibin = i;
  // 	  found = true;
  // 	  break; 
  // 	}
  //   }


  std::unordered_map<std::size_t,std::size_t>::const_iterator
    p = mGlbIdxMap.find( gidx );
  if ( p != mGlbIdxMap.end() )
    {
      ibin = p->second;
      found = true;
    }
  
  //std::cout << "ibin,found " << ibin << " " << found << std::endl;

  if ( mInterpMode == 5 )
    {

      std::size_t dbin = mNBspl/2;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  std::size_t dsize = mDimInfo.GetDimSizes()[dim];
	  if ( bidx[dim] < dbin or bidx[dim]+dbin >= dsize )
	    {
	      found = false;
	      break;
	    }
	};
    }

  // std::printf("pt ");
  // for ( std::size_t dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf("%10.5f",inppt[dim]);
  //   }
  // if ( found )
  //   {
  //     std::printf(" gidx: %6lu is  found, nbin: %6lu\n",gidx,mNumBins);
  //   }
  // else
  //   {
  //     std::printf(" gidx: %6lu NOT found, nbin: %6lu\n",gidx,mNumBins);
  //   }

  
  bool calc_pen = not found;
  if ( found )
    {
      if ( mBinSizes[ibin] == 0 )
	{
	  calc_pen = true;
	}
    }

  if ( calc_pen )
    {
      std::size_t jbin = 0;
      
      double mindist = 1.e+30;
      for ( std::size_t i=0; i<mNumBins; ++i )
	{
	  if ( mBinSizes[i] == 0 ) // here: only look at occupied bins
	    {
	      continue;
	    }
	  double dist = 0;
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      double dx = pt[dim] - mBins[i].center[dim];
	      if ( mDimInfo.IsPeriodic(dim) )
		{
		  dx = wrap(dx,360.);
		}
	      dist += dx*dx;
	    }
	  if ( dist < mindist )
	    {
	      mindist = dist;
	      jbin = i;
	    }
	}

      if ( not found )
	{
	  ibin = jbin;
	}

      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  double dx = pt[dim] - mBins[jbin].center[dim];
	  if ( mDimInfo.IsPeriodic(dim) )
	    {
	      dx = wrap(dx,360.);
	    }
	  pen += oobk * dx * dx;
	  dpen[dim] = 2 * oobk * dx;
	  // ok so this makes sure that pt lives in an occupied bin
	  pt[dim] = mBins[jbin].center[dim];
	}

      
    }

  
  //if ( found )
    {
      val = 0;
      if ( mInterpMode == 1 )
	{
	  val = mBinVals[ibin];
	}
      else if ( mInterpMode == 2 )
	{
	  EvalRBF( 1, pt.data(), &val, grd );
	}
      else if ( mInterpMode == 3 )
	{
	  std::vector<double> work( mNumCorners, 0 );
	  for ( std::size_t ic=0; ic<mNumCorners; ++ic )
	    {
	      work[ic] = mGlbCvals[ mLocCidxs[ic+ibin*mNumCorners] ];
	    }
	  
	  std::size_t bsplorder = mOrder;
	  std::size_t nbspl = mNBspl;
	  
	  double const * widths = mDimInfo.GetTargetWidths();
	  //for ( std::size_t isub=0; isub<mSubBinSize; ++isub )
	    {
	      std::vector< std::vector<double> > linsubwts(ndim);
	      std::vector< std::vector<double> > linsubdwts(ndim);
	      std::vector<double> gwts(bsplorder,0.);
	      std::vector<double> gdwts(bsplorder,0.);
	      std::vector<int> gidx(bsplorder,0);
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  linsubwts[dim].assign(nbspl,0.);
		  linsubdwts[dim].assign(nbspl,0.);
		  std::fill( gwts.data(), gwts.data() + bsplorder, 0. );
		  double w = (nbspl-1)*widths[dim];
		  double c = w/2.;
		  double dx = (pt[dim] - mBins[ibin].center[dim]) + c;
		  
		  ndfes::bspline_aperiodic( dx, w, nbspl, bsplorder,
					    gwts.data(), gdwts.data(), gidx.data() );

		  // {
		  //   double const DEL = 5.e-5;
		  //   std::vector<double> whi(bsplorder,0.);
		  //   std::vector<double> wlo(bsplorder,0.);
		  //   std::vector<double> dw(bsplorder,0.);
		  //   ndfes::bspline_aperiodic( dx+DEL, w, nbspl, bsplorder,
		  // 			      whi.data(), gidx.data() );
		  //   ndfes::bspline_aperiodic( dx-DEL, w, nbspl, bsplorder,
		  // 			      wlo.data(), gidx.data() );
		    
		  // }
		  
		  for ( std::size_t ib=0; ib<bsplorder; ++ib )
		    {
		      linsubwts[dim][ gidx[ib] ] = gwts[ib];
		      linsubdwts[dim][ gidx[ib] ] = gdwts[ib];
		    }
		}
	      std::vector<double> gridsubwts;
	      ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
	      
	      val = 0;
	      for ( std::size_t i=0; i<mNumCorners; ++i )
		{
		  val += gridsubwts[i] * work[i];
		}

	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  std::vector< std::vector<double> > tmplin( linsubwts );
		  tmplin[dim] = linsubdwts[dim];
		  ndfes::LinearWtsToMeshWts(tmplin,gridsubwts);
		  for ( std::size_t i=0; i<mNumCorners; ++i )
		    {
		      grd[dim] += gridsubwts[i] * work[i];
		    }
		}
	      
	    }
	}
      else if ( mInterpMode == 4 )
	{

	  // std::cout << "ibin,sizes" << ibin
	  // 	    << " " << mARBFCwts.size()
	  // 	    << " " << mARBFCpts.size()
	  // 	    << " " << mARBFCwts.size()
	  // 	    << " " << pt.size()
	  // 	    << std::endl;
	  
	  std::size_t ncpts = mARBFCwts[ibin].size();
	  RbfEval( ndim, ncpts, mARBFCpts[ibin].data(), mARBFCwts[ibin].data(),
		   mRBFShape, 1, pt.data(), &val, grd );
	}
      else if ( mInterpMode == 5 )
	{
	  EvalWAvgValueAndGrad(pt.data(),val,grd);
	}
    }


  val += pen;
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      grd[dim] += dpen[dim];
    }
  //std::cout << "Return" << std::endl;
  return found;
}





bool ndfes::FES::GetBiasedValue
( double const * rc,
  double const * fc,
  double const oobk,
  double const * inppt,
  double & val,
  double * grd ) const
{
  bool ok = GetPenalizedValue( oobk, inppt, val, grd );

  std::size_t const ndim = mDimInfo.GetNumDims();
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      double dx = inppt[dim] - rc[dim];
      if ( mDimInfo.IsPeriodic(dim) )
	{
	  dx = wrap(dx,360.);
	}
      val += fc[dim] * dx * dx;
      grd[dim] += 2. * fc[dim] * dx;
    }
  
  return ok;
}



bool ndfes::FES::PointIsOcc( double const * pt ) const
{
  std::vector<std::size_t> bidx( mDimInfo.CptBinIdxs( pt ) );
  std::size_t gidx = mDimInfo.CptGlbBinIdx( bidx.data() );
  //std::size_t ndim = mDimInfo.GetNumDims();

  bool found = false;
  //std::size_t ibin = 0;

  std::unordered_map<std::size_t,std::size_t>::const_iterator
    p = mGlbIdxMap.find(gidx);
  if ( p != mGlbIdxMap.end() )
    {
      found = true;
    }
  
  // for ( std::size_t i=0; i < mNumBins; ++i )
  //   {
  //     if ( mBins[i].glbidx == gidx )
  // 	{
  // 	  //ibin = i;
  // 	  found = true;
  // 	  break;
  // 	}
  //   }
  return found;
}

std::size_t ndfes::FES::BinOccAtPt( double const * pt ) const
{
  std::vector<std::size_t> bidx( mDimInfo.CptBinIdxs( pt ) );
  std::size_t gidx = mDimInfo.CptGlbBinIdx( bidx.data() );
  //std::size_t ndim = mDimInfo.GetNumDims();

  std::size_t n=0;
  
  std::unordered_map<std::size_t,std::size_t>::const_iterator
    p = mGlbIdxMap.find(gidx);
  if ( p != mGlbIdxMap.end() )
    {
      n = mBinSizes[ p->second ];
    }
  
  // for ( std::size_t i=0; i < mNumBins; ++i )
  //   {
  //     if ( mBins[i].glbidx == gidx )
  // 	{
  // 	  //ibin = i;
  // 	  found = true;
  // 	  break;
  // 	}
  //   }
  return n;
}


bool ndfes::FES::PointIsInRange( double const * pt ) const
{
  bool good = true;
  std::size_t ndim = mDimInfo.GetNumDims();
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      if ( pt[dim] < mDimInfo.GetXmin()[dim] or
	   pt[dim] > mDimInfo.GetXmax()[dim] )
	{
	  good = false;
	}
    }
  return good;
}

bool ndfes::FES::PointIsInRangeAndNonzeroOcc( double const * pt ) const
{
  bool good = PointIsInRange(pt);
  if ( good )
    {
      std::vector<std::size_t> bidx( mDimInfo.CptBinIdxs( pt ) );
      std::size_t gidx = mDimInfo.CptGlbBinIdx( bidx.data() );
      //std::size_t ndim = mDimInfo.GetNumDims();
      
      bool found = false;
      std::size_t ibin = 0;
      std::unordered_map<std::size_t,std::size_t>::const_iterator
	p = mGlbIdxMap.find(gidx);
      if ( p != mGlbIdxMap.end() )
	{
	  ibin = p->second;
	  found = true;
	}
      // for ( std::size_t i=0; i < mNumBins; ++i )
      // 	{
      // 	  if ( mBins[i].glbidx == gidx )
      // 	    {
      // 	      ibin = i;
      // 	      found = true;
      // 	      break;
      // 	    }
      // 	}
      if ( found )
	{
	  good = mBinSizes[ibin] > 0;
	}
      else
	{
	  good = false;
	}
    }
  return good;
}



std::shared_ptr<ndfes::FES>
ndfes::FES::AddPenaltyBuffer( std::size_t const minsize, std::size_t const nlayers ) const
{
  std::shared_ptr<ndfes::FES> xFES( AddPenaltyBuffer(minsize) );
  for ( std::size_t i=1; i<nlayers; ++i )
    {
      xFES = xFES->AddPenaltyBuffer(0);
    }
  return xFES;
}



///////

struct TripleT
{
  TripleT() : first(0), value(-1.e+30) {}
  
  TripleT( std::size_t const gidx, std::vector<std::size_t> const & bidx, double const fe )
    : first(gidx), second(bidx), value(fe) {}
  
  std::size_t first;
  std::vector<std::size_t> second;
  double value;
};

typedef std::vector< TripleT > vtriple;
typedef vtriple::iterator itriple;

typedef std::unordered_map< std::size_t, TripleT > mtriple;
typedef mtriple::iterator imtriple;

struct FirstEquals
{
  FirstEquals( std::size_t v ) : first(v) {}
  std::size_t first;
};

bool operator== ( TripleT const & a, FirstEquals const & b )
{
  return a.first == b.first;
}

bool operator== ( TripleT const & a, TripleT const & b )
{
  return a.first == b.first;
}

bool operator< ( TripleT const & a, TripleT const & b )
{
  return a.first < b.first;
}
///////




std::shared_ptr<ndfes::FES>
ndfes::FES::AddPenaltyBuffer( std::size_t const minsize ) const
{

  double const DELTAV = 0.5;
  
  std::size_t const ndim = mDimInfo.GetNumDims();
  
  std::shared_ptr<ndfes::FES> xFES( new ndfes::FES() );
  xFES->mDimInfo = ndfes::MakeBufferedGrid( mDimInfo, 1 );
  xFES->mIsMBAR = true;
  xFES->mInterpMode = 0;
  xFES->mOrder = mOrder;
  xFES->mNBspl = 0;
  xFES->mNumCorners = 0;
  xFES->mRBFShape = mRBFShape;
  xFES->mNumQuad = mNumQuad;
  xFES->mSubBinSize = mSubBinSize;
  xFES->mARBFDelta = mARBFDelta;
  xFES->mWithInverse = mWithInverse;


  //std::vector<std::size_t> oldgidxs;
  typedef std::unordered_set<std::size_t> useti_t;
  useti_t oldgidxs;

  //std::printf("\n");
  std::size_t nbins = mBins.size();
  for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      if ( mBinSizes[ibin] >= minsize )
	{
	  double const * c = mBins[ibin].center.data();
	  std::vector<std::size_t> bidxs( xFES->mDimInfo.CptBinIdxs( c ) );
	  std::size_t gidx = xFES->mDimInfo.CptGlbBinIdx( bidxs.data() );
	  xFES->mBins.push_back(ndfes::SpatialBin(xFES->mDimInfo,gidx,bidxs.data()));
	  xFES->mBinVals.push_back( mBinVals[ibin] );
	  xFES->mBinErrs.push_back( mBinErrs[ibin] );
	  xFES->mBinREs.push_back( mBinREs[ibin] );
	  xFES->mBinSizes.push_back( mBinSizes[ibin] );
	  //oldgidxs.push_back( gidx );
	  oldgidxs.insert( gidx );
	  
	  // std::printf("%6lu %6lu %5lu %6lu ",ibin,mBins[ibin].glbidx,gidx,mBinSizes[ibin]);
	  // for ( std::size_t dim=0; dim<ndim; ++dim )
	  //   {
	  //     std::printf("%12.5f",c[dim]);
	  //   }
	  // std::printf("\n");
	}
    }

  mtriple newidxs;
  
  //std::vector<std::size_t> newgidxs;
  //std::vector< std::vector<std::size_t> > newbidxs;

  std::size_t const nxbins = xFES->mBins.size();

#ifdef WITH_OPENMP
#pragma omp parallel
  {
    mtriple mynewidxs;
#pragma omp for schedule(dynamic)
#endif
    for ( std::size_t ibin=0; ibin < nxbins; ++ibin )
      {
	double const fevalue = xFES->mBinVals[ibin];
	std::vector<std::size_t> midxs
	  ( ndfes::GetMeshIdxs(xFES->mDimInfo,1,
			       xFES->mBins[ibin].bidxs.data(),true));
	std::size_t nmidxs = midxs.size()/ndim;
	for ( std::size_t i=0; i<nmidxs; ++i )
	  {
	    std::vector<std::size_t> bidx( midxs.data() + i*ndim,
					   midxs.data() + (i+1)*ndim );
	    std::size_t gidx = xFES->mDimInfo.CptGlbBinIdx( bidx.data() );
	    
	    //std::vector<std::size_t>::iterator q
	    // = std::find( oldgidxs.begin(), oldgidxs.end(), gidx );
	    useti_t::iterator q = oldgidxs.find(gidx);
	    
	    if ( q == oldgidxs.end() )
	      {
		TripleT val( gidx, bidx, fevalue );
#ifdef WITH_OPENMP
		//itriple p = std::find( mynewidxs.begin(), mynewidxs.end(), val );
		imtriple p = mynewidxs.find(gidx);
		if ( p == mynewidxs.end() )
		  {
		    //mynewidxs.push_back( val );
		    mynewidxs.insert( {gidx,val} );
		  }
		else
		  {
		    p->second.value = std::max( p->second.value, fevalue );
		  }
#else
		//itriple p = std::find( newidxs.begin(), newidxs.end(), val );
		imtriple p = newidxs.find(gidx);
		if ( p == newidxs.end() )
		  {
		    //newidxs.push_back( val );
		    newidxs.insert( {gidx,val} );
		  }
		else
		  {
		    p->second.value = std::max( p->second.value, fevalue );
		  }
#endif
	      }
	  }
      }
#ifdef WITH_OPENMP
#pragma omp critical
    {
      newidxs.reserve( newidxs.size() + std::distance( mynewidxs.begin(), mynewidxs.end() ) );
      // for ( itriple it=mynewidxs.begin(); it != mynewidxs.end(); ++it )
      // 	{
      // 	  itriple p = std::find( newidxs.begin(), newidxs.end(), *it );
      // 	  if ( p == newidxs.end() )
      // 	    {
      // 	      newidxs.push_back( *it );
      // 	    }
      // 	  else
      // 	    {
      // 	      p->value = std::max( p->value, it->value );
      // 	    }
      // 	}
      for ( imtriple it=mynewidxs.begin(); it != mynewidxs.end(); ++it )
	{
	  imtriple p = newidxs.find(it->first);
	  if ( p == newidxs.end() )
	    {
	      newidxs.insert( *it );
	    }
	  else
	    {
	      p->second.value = std::max( p->second.value, it->second.value );
	    }
	}
    }
  } // omp parallel
#endif

  
  // std::size_t const nnew = newidxs.size();
  // for ( std::size_t inew=0; inew<nnew; ++inew )
  //   {
  //     std::vector<std::size_t> bidxs( newidxs[inew].second );
  //     std::size_t gidx = newidxs[inew].first;
  //     double v = newidxs[inew].value;
  //     {
  // 	xFES->mBins.push_back(ndfes::SpatialBin(xFES->mDimInfo,gidx,bidxs.data()));
  // 	xFES->mBinVals.push_back( v + DELTAV );
  // 	xFES->mBinErrs.push_back( 0. );
  // 	xFES->mBinREs.push_back( 0. );
  // 	xFES->mBinSizes.push_back( 0 );
  //     }
  //   }

  std::cout << "Had " << xFES->mBins.size() << " bins; adding "
	    << newidxs.size() << " for a total of "
	    << xFES->mBins.size() + newidxs.size()
	    << std::endl;
  
  for ( imtriple p = newidxs.begin(); p != newidxs.end(); ++p )
    {
      std::vector<std::size_t> bidxs( p->second.second );
      std::size_t gidx = p->first;
      double v = p->second.value;
      {
	xFES->mBins.push_back(ndfes::SpatialBin(xFES->mDimInfo,gidx,bidxs.data()));
	xFES->mBinVals.push_back( v + DELTAV );
	xFES->mBinErrs.push_back( 0. );
	xFES->mBinREs.push_back( 0. );
	xFES->mBinSizes.push_back( 0 );
      }
    }

  
  
  xFES->mNumBins = xFES->mBins.size();
  
  return xFES;
}

std::shared_ptr<ndfes::FES>
ndfes::FES::DeleteBinsBySize( std::size_t const minsize ) const
{
  
  std::shared_ptr<ndfes::FES> xFES( new ndfes::FES() );
  xFES->mDimInfo = mDimInfo;
  xFES->mIsMBAR = true;
  xFES->mInterpMode = 0;
  xFES->mOrder = mOrder;
  xFES->mNBspl = 0;
  xFES->mNumCorners = 0;
  xFES->mRBFShape = mRBFShape;
  xFES->mNumQuad = mNumQuad;
  xFES->mSubBinSize = mSubBinSize;
  xFES->mARBFDelta = mARBFDelta;
  xFES->mWithInverse = mWithInverse;


  std::size_t nbins = mBins.size();
  for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      if ( mBinSizes[ibin] >= minsize )
	{
	  double const * c = mBins[ibin].center.data();
	  std::vector<std::size_t> bidxs( xFES->mDimInfo.CptBinIdxs( c ) );
	  std::size_t gidx = xFES->mDimInfo.CptGlbBinIdx( bidxs.data() );
	  xFES->mBins.push_back(ndfes::SpatialBin(xFES->mDimInfo,gidx,bidxs.data()));
	  xFES->mBinVals.push_back( mBinVals[ibin] );
	  xFES->mBinErrs.push_back( mBinErrs[ibin] );
	  xFES->mBinREs.push_back( mBinREs[ibin] );
	  xFES->mBinSizes.push_back( mBinSizes[ibin] );
	}
    }

  xFES->mNumBins = xFES->mBins.size();

  std::cout << "Pruning bins with fewer than "
	    << minsize << " samples. Before: "
	    << mNumBins << " After: "
	    << xFES->mNumBins << std::endl;
  
  return xFES;
}


void ndfes::FES::AddOccPen( double const occpen )
{
   std::size_t nbins = mBins.size();
   for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      mBinVals[ibin] += occpen * mBinSizes[ibin];
    }
}

std::shared_ptr<ndfes::FES>
ndfes::FES::AddPenaltyBuffer_ORIG( std::size_t const minsize ) const
{

  double const DELTAV = 0.5;
  
  std::size_t const ndim = mDimInfo.GetNumDims();
  
  std::shared_ptr<ndfes::FES> xFES( new ndfes::FES() );
  xFES->mDimInfo = ndfes::MakeBufferedGrid( mDimInfo, 1 );
  xFES->mIsMBAR = true;
  xFES->mInterpMode = 0;
  xFES->mOrder = mOrder;
  xFES->mNBspl = 0;
  xFES->mNumCorners = 0;
  xFES->mRBFShape = mRBFShape;
  xFES->mNumQuad = mNumQuad;
  xFES->mSubBinSize = mSubBinSize;
  xFES->mARBFDelta = mARBFDelta;
  xFES->mWithInverse = mWithInverse;


  std::vector<std::size_t> oldgidxs;
  
  std::size_t nbins = mBins.size();
  for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      if ( mBinSizes[ibin] >= minsize )
	{
	  double const * c = mBins[ibin].center.data();
	  std::vector<std::size_t> bidxs( xFES->mDimInfo.CptBinIdxs( c ) );
	  std::size_t gidx = xFES->mDimInfo.CptGlbBinIdx( bidxs.data() );
	  xFES->mBins.push_back(ndfes::SpatialBin(xFES->mDimInfo,gidx,bidxs.data()));
	  xFES->mBinVals.push_back( mBinVals[ibin] );
	  xFES->mBinErrs.push_back( mBinErrs[ibin] );
	  xFES->mBinREs.push_back( mBinREs[ibin] );
	  xFES->mBinSizes.push_back( mBinSizes[ibin] );
	  oldgidxs.push_back( gidx );
	}
    }

  std::vector<std::size_t> newgidxs;
  std::vector< std::vector<std::size_t> > newbidxs;

  for ( std::size_t ibin=0; ibin < xFES->mBins.size(); ++ibin )
    {
      std::vector<std::size_t> midxs
	( ndfes::GetMeshIdxs(xFES->mDimInfo,1,
			     xFES->mBins[ibin].bidxs.data(),true));
      std::size_t nmidxs = midxs.size()/ndim;
      for ( std::size_t i=0; i<nmidxs; ++i )
	{
	  std::vector<std::size_t> bidx( midxs.data() + i*ndim,
					 midxs.data() + (i+1)*ndim );
	  std::size_t gidx = xFES->mDimInfo.CptGlbBinIdx( bidx.data() );
	  std::vector<std::size_t>::iterator p
	    = std::find( newgidxs.begin(), newgidxs.end(), gidx );
	  std::vector<std::size_t>::iterator q
	    = std::find( oldgidxs.begin(), oldgidxs.end(), gidx );
	  
	  if ( p == newgidxs.end() and q == oldgidxs.end() )
	    {
	      newgidxs.push_back( gidx );
	      newbidxs.push_back( bidx );
	    }
	}
    }


  for ( std::size_t inew=0; inew<newgidxs.size(); ++inew )
    {
      std::vector<std::size_t> midxs
	( ndfes::GetMeshIdxs(xFES->mDimInfo,1,
			     newbidxs[inew].data(),false));
      std::size_t nmidxs = midxs.size()/ndim;

      //std::cout << "inew, nmidxs " << inew << " " << nmidxs << std::endl;
      std::size_t nfound = 0;
      double v = -1.e+30;
      for ( std::size_t i=0; i<nmidxs; ++i )
	{
	  std::vector<std::size_t> bidx( midxs.data() + i*ndim,
				    midxs.data() + (i+1)*ndim );
	  std::size_t gidx = xFES->mDimInfo.CptGlbBinIdx( bidx.data() );

	  std::vector<std::size_t>::iterator p
	    = std::find( oldgidxs.begin(), oldgidxs.end(), gidx );
	  
	  if ( p != oldgidxs.end() )
	    {
	      std::size_t ibin = std::distance( oldgidxs.begin(), p );
	      //v += xFES->mBinVals[ibin];
	      v = std::max(v,xFES->mBinVals[ibin]);
	      ++nfound;
	    }
	}
      if ( nfound == 0 )
	{
	  std::cerr << "There was an internal error creating the grid buffer" << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      //v /= nfound;
      std::vector<std::size_t> bidxs( newbidxs[inew] );
      std::size_t gidx = newgidxs[inew];
      xFES->mBins.push_back(ndfes::SpatialBin(xFES->mDimInfo,gidx,bidxs.data()));
      xFES->mBinVals.push_back( v + DELTAV );
      xFES->mBinErrs.push_back( 0. );
      xFES->mBinREs.push_back( 0. );
      xFES->mBinSizes.push_back( 0 );
    }
  
  xFES->mNumBins = xFES->mBins.size();


  // for ( std::size_t i=0; i<xFES->mBins.size(); ++i )
  //   {
  //     if ( i < oldgidxs.size() )
  // 	{
  // 	  std::printf("old %6lu",xFES->mBins[i].glbidx);
  // 	}
  //     else
  // 	{
  // 	  std::printf("new %6lu",xFES->mBins[i].glbidx);
  // 	}
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf("%9.4f",xFES->mBins[i].center[dim]);
  // 	}
  //     std::printf("%12.3e\n",xFES->mBinVals[i]);
  //   }
  
  
  
  return xFES;
}




int ndfes::FES::NumLayersOfOcc( double const * pt ) const
{
  //typedef std::vector<std::size_t>::iterator vit;
  
  int maxlayers = -1;
  std::size_t const ndim = mDimInfo.GetNumDims();
  
  std::vector<std::size_t> bidx( mDimInfo.CptBinIdxs( pt ) );
  std::size_t gidx = mDimInfo.CptGlbBinIdx( bidx.data() );

  // std::vector<std::size_t> gidxs(mNumBins);
  // for ( std::size_t i=0; i<mNumBins; ++i )
  //   {
  //     gidxs[i] = mBins[i].glbidx;
  //   }
  
  // bool found = false;
  // {
  //   vit p = std::find(gidxs.begin(),gidxs.end(),gidx);
  //   if ( p != gidxs.end() )
  //     {
  // 	found = true;
  // 	maxlayers = 0;
  //     }
  // }

  bool found = false;
  {
    std::unordered_map<std::size_t,std::size_t>::const_iterator
      p = mGlbIdxMap.find(gidx);
    if ( p != mGlbIdxMap.end() )
      {
	found = true;
      }
  }

  if ( found )
    {
      for ( std::size_t ilayer=1; ilayer < 4; ++ilayer )
	{
	  std::vector<std::size_t> midxs
	    ( ndfes::GetMeshIdxs(mDimInfo,ilayer,
				 bidx.data(),false));
	  std::size_t nm = midxs.size() / ndim;
	  if ( nm != std::pow(2*ilayer+1,ndim) )
	    {
	      break;
	    }
	  else
	    {
	      bool foundall = true;
	      for ( std::size_t i=0; i<nm; ++i )
		{
		  std::size_t gidx = mDimInfo.CptGlbBinIdx( midxs.data() + i*ndim );

		  std::unordered_map<std::size_t,std::size_t>::const_iterator
		    p = mGlbIdxMap.find(gidx);
		  if ( p == mGlbIdxMap.end() )
		    {
		      foundall = false;
		      break;
		    }
		  
		  // vit p = std::find(gidxs.begin(),gidxs.end(),gidx);
		  // if ( p == gidxs.end() )
		  //   {
		  //     foundall = false;
		  //     break;
		  //   }
		}
	      if ( foundall )
		{
		  ++maxlayers;
		}
	    }
	}
    }
  return maxlayers;
}







ndfes::FES::FES( std::string xmlname, int const imodel )
  : mIsMBAR(false),
    mInterpMode(0),
    mNumBins(0),
    mOrder(0),
    mNBspl(0),
    mNumCorners(0),
    mRBFShape(100),
    mNumQuad(0),
    mSubBinSize(0),
    mARBFDelta(-1),
    mWithInverse(false)
{

  std::string mstr;
  {
    std::stringstream ss;
    ss << imodel;
    mstr = ss.str();
  }
  
  XMLIO::xml_document xml;
  XMLIO::LoadXml( xmlname, xml );

  XMLIO::xml_node ndfesnode( XMLIO::ExtractNode(xml,"ndfes") );

  //bool foundmodel = false;

  XMLIO::xml_node modelnode = ndfesnode.find_child_by_attribute("model","idx",mstr.c_str());

  if ( ! modelnode )
    {
      std::cerr << "Failed to find <model idx=\"" << mstr << "\"> in "
   		<< xmlname << std::endl;
      std::exit(EXIT_FAILURE);
    }
  

  
  {
    std::string modeltype = "MBAR";
    XMLIO::ExtractValue(modelnode,"type",modeltype);
    strop::trim(modeltype);
    std::transform( modeltype.begin(), modeltype.end(),
		    modeltype.begin(), ::toupper);
    if ( modeltype.compare("MBAR") == 0 )
      {
	mIsMBAR = true;
      }
    else if ( modeltype.compare("VFEP") == 0 )
      {
	mIsMBAR = false;
      }
    else
      {
	std::cerr << "Invalid <type> value: " << modeltype
		  << " found in model " << imodel
		  << " within " << xmlname << std::endl;
	std::exit(EXIT_FAILURE);
      };
  }

  mOrder = 1;
  if ( not mIsMBAR )
    {
      int order = 0;
      XMLIO::ExtractValue(modelnode,"order",order);
      mOrder = order;
    }

  std::size_t ndim = 0;
  {
    XMLIO::xml_node gnode( XMLIO::ExtractNode(modelnode,"grid") );
    for ( XMLIO::xml_node dnode = gnode.first_child();
	  dnode; dnode = dnode.next_sibling() )
      {
	++ndim;
      }

    std::vector<double> xmins(ndim,0);
    std::vector<double> xmaxs(ndim,0);
    std::vector<int> sizes(ndim,0);
    std::vector<int> ispers(ndim,0);

    for ( XMLIO::xml_node dnode = gnode.first_child();
	  dnode; dnode = dnode.next_sibling() )
      {
	int idx = 0;
	XMLIO::ExtractAttribute(dnode,"idx",idx);
	if ( idx < 0 or idx >= (int)ndim )
	  {
	    std::cerr << "<grid><dim idx='" << idx
		      << "'> is out of bounds"
		      << std::endl;
	    std::exit(EXIT_FAILURE);
	  }
	XMLIO::ExtractValue(dnode,"xmin",xmins[idx]);
 	XMLIO::ExtractValue(dnode,"xmax",xmaxs[idx]);
 	XMLIO::ExtractValue(dnode,"size",sizes[idx]);
 	XMLIO::ExtractValue(dnode,"isper",ispers[idx]);
	//std::cout << idx << " " << ispers[idx] << "\n";
     }


    mDimInfo = ndfes::DimInfo
      ( xmins.size(), mOrder,
	ispers.data(), xmins.data(),
	xmaxs.data(), sizes.data() );
    
  }


  for ( XMLIO::xml_node bnode = modelnode.child("bin");
	bnode; bnode = bnode.next_sibling("bin") )
    {
      std::size_t gidx=0;
      XMLIO::ExtractAttribute(bnode,"idx",gidx);
      std::vector<std::size_t> bidxs(ndim,0);
      for ( XMLIO::xml_node tnode = bnode.child("bidx");
	    tnode; tnode = tnode.next_sibling("bidx") )
	{
	  int idim=0;
	  XMLIO::ExtractAttribute(tnode,"idx",idim);
	  XMLIO::ExtractValue(tnode,bidxs[idim]);
	}
      double value = 0;
      double stderr = 0;
      double entropy = 0;
      int size = 0;
      XMLIO::ExtractOptionalValue(bnode,"val",value);
      XMLIO::ExtractOptionalValue(bnode,"err",stderr);
      XMLIO::ExtractOptionalValue(bnode,"re",entropy);
      XMLIO::ExtractOptionalValue(bnode,"size",size);

      mBinVals.push_back( value );
      mBinErrs.push_back( stderr );
      mBinREs.push_back( entropy );
      mBinSizes.push_back( size );
      mBins.push_back( ndfes::SpatialBin( mDimInfo, gidx, bidxs.data() ) );
      mNumBins = mBins.size();
    }

  
  for ( XMLIO::xml_node cnode = modelnode.child("corner");
	cnode; cnode = cnode.next_sibling("corner") )
    {
      std::size_t cidx = 0;
      double value = 0;
      double stderr = 0;
      XMLIO::ExtractAttribute(cnode,"idx",cidx);
      XMLIO::ExtractValue(cnode,"val",value);
      XMLIO::ExtractValue(cnode,"err",stderr);
      mGlbCidxs.push_back( cidx );
      mGlbCvals.push_back( value );
      mGlbCerrs.push_back( stderr );
    }
  mNumCorners = mDimInfo.GetNumCorners();

  if ( not mIsMBAR )
    {
      PreparevFEP();
    }
  
}
