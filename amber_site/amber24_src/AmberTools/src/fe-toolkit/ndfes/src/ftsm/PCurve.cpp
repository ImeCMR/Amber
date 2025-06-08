#include <cmath>
#include <cstdio>

#include "PCurve.hpp"
#include "LinearInterp.hpp"

void ApproxTFromPts
( int const ndim,
  int const npts,
  double const * pts,
  double * ts )
{
  ts[0] = 0.;
  for ( int i=0; i<npts-1; ++i )
    {
      double dist = 0;
      for ( int k=0; k<ndim; ++k )
	{
	  double d = pts[k+(i+1)*ndim] - pts[k+i*ndim];
	  dist += d*d;
	}
      dist = std::sqrt(dist);
      ts[i+1] = ts[i] + dist;
    }
  double maxt = ts[npts-1];
  for ( int i=0; i<npts; ++i )
    {
      ts[i] /= maxt;
    };
}

void ApproxTFromSpl
( int const ndim,
  int const nseg,
  std::vector< ccdl::AkimaSpline > const & pspls,
  double * ts )
{
  std::vector<double> pts(ndim*nseg,0.);
  double nm1 = nseg-1;
  for ( int i=0; i<nseg; ++i )
    {
      double t = i/nm1;
      for ( int k=0; k<ndim; ++k )
	{
	  pts[k+i*ndim] = pspls[k].GetValue(t);
	}
    }
  ApproxTFromPts( ndim, nseg, pts.data(), ts );
}


void CptNewTs
( int const ndim,
  int const nseg,
  std::vector< ccdl::AkimaSpline > const & pspls,
  double * tnew )
{
  std::vector<double> tin(nseg,0);
  std::vector<double> tout(nseg,0);
  double nm1 = nseg-1;
  for ( int i=0; i<nseg; ++i )
    {
      tin[i] = i/nm1;
    }
  ApproxTFromSpl(ndim,nseg,pspls,tout.data() );
  
  ndfes::LinearInterp tspl(nseg,tin.data(),tout.data());
  
  int ncntrl = pspls[0].GetNumPts();
  double const * tcntrl = pspls[0].GetXValues();
  for ( int i=0; i<ncntrl; ++i )
    {
      tnew[i] = tspl.GetValue( tcntrl[i] );
    }
  tnew[0] = 0;
  tnew[ncntrl-1] = 1;
}



ndfes::PCurve::PCurve
( int myndim,
  int mnpts,
  double const * mpts,
  double const * mts,
  bool const akima )
  : mNumDim(myndim),
    mNumPts(mnpts),
    mT(mts,mts+mnpts),
    mX(mpts,mpts+myndim*mnpts),
    mAkima(akima)
{
  std::vector<double> data(mNumPts,0);
  if ( mAkima )
    {
      mASpl.resize(mNumDim);
      for ( int k=0; k<mNumDim; ++k )
	{
	  for ( int i=0; i<mNumPts; ++i )
	    {
	      data[i] = mX[k+i*mNumDim];
	    }
	  mASpl[k] = ccdl::AkimaSpline( mNumPts, mT.data(), data.data() );
	}
    }
  else
    {
      mLSpl.resize(mNumDim);
      for ( int k=0; k<mNumDim; ++k )
	{
	  for ( int i=0; i<mNumPts; ++i )
	    {
	      data[i] = mX[k+i*mNumDim];
	    }
	  mLSpl[k] = ndfes::LinearInterp( mNumPts, mT.data(), data.data() );
	}
    }
}



ndfes::PCurve::PCurve
( int myndim,
  int mnpts,
  double const * mpts,
  bool const akima,
  int const maxit, int const nseg )
  : mNumDim(myndim),
    mNumPts(mnpts),
    mT(mnpts,0),
    mX(mpts,mpts+myndim*mnpts),
    mAkima(akima)
{
  double const TOL = 1.e-8;

  int mynseg = nseg;
  int maxnseg = 16*nseg;

  std::vector<double> data(mNumPts*mNumDim,0);
      
  for ( int k=0; k<mNumDim; ++k )
    {
      for ( int i=0; i<mNumPts; ++i )
	{
	  data[i+k*mNumPts] = mX[k+i*mNumDim];
	}
    }
  
  ApproxTFromPts( mNumDim, mNumPts, mX.data(), mT.data() );

  
  if ( not mAkima )
    {
      mLSpl.resize( mNumDim );

      for ( int k=0; k<mNumDim; ++k )
	{
	  mLSpl[k] = ndfes::LinearInterp( mNumPts, mT.data(),
					  data.data()+k*mNumPts );
	}
      
    }
  else
    {

      mASpl.resize( mNumDim );
  
      
      for ( int k=0; k<mNumDim; ++k )
	{
	  mASpl[k] = ccdl::AkimaSpline( mNumPts, mT.data(),
					data.data()+k*mNumPts );
	}
      
      std::vector<double> tnew(mT);
      
      for ( int it=0; it<maxit; ++it )
	{
	  CptNewTs(mNumDim,mynseg,mASpl,tnew.data());
	  
	  bool bad = false;
	  for ( int i=1; i<mNumPts; ++i )
	    {
	      if ( std::abs(tnew[i]-tnew[i-1]) < 1.e-6 )
		{
		  bad = true;
		}
	    }
	  if ( bad )
	    {
	      if ( mynseg < maxnseg )
		{
		  mynseg *= 2;
		  CptNewTs(mNumDim,mynseg,mASpl,tnew.data());
		  bad = false;
		  for ( int i=1; i<mNumPts; ++i )
		    {
		      if ( std::abs(tnew[i]-tnew[i-1]) < 1.e-8 )
			{
			  bad = true;
			}
		    }
		}
	      if ( bad )
		{
		  break;
		};
	    }
       
	  
	  // std::printf("it=%3i nseg=%5i\n",it,mynseg);
	  // for ( int i=1; i<mNumPts; ++i )
	  //   {
	  //     //if ( std::abs(tnew[i]-tnew[i-1]) < 1.e-8 )
	  // 	{
	  // 	  std::printf("%4i %13.4e %13.4e %13.4e\n",
	  // 		      i,mT[i],tnew[i],tnew[i]-tnew[i-1]);
	  // 	};
	  //   }
	  
	  double err = 0;
	  for ( int i=0; i<mNumPts; ++i )
	    {
	      double d = mT[i]-tnew[i];
	      err += d*d;
	      mT[i] = tnew[i];
	    }
	  err = std::sqrt(err);

	  
	  for ( int k=0; k<mNumDim; ++k )
	    {
	      mASpl[k] = ccdl::AkimaSpline( mNumPts, mT.data(),
					    data.data()+k*mNumPts );
	    }
	  
	  if ( err < TOL )
	    {
	      break;
	    }
	  
	}
    }
  
  // std::printf("\n");
  // for ( std::size_t i=0; i<mNumPts; ++i )
  //   {
  //     std::printf("%6lu %17.8e",i,mT[i]);
  //     for ( int k=0; k<mNumDim; ++k )
  // 	{
  // 	  std::printf(" %17.8e",mX[k+i*mNumDim]);
  // 	}
  //     std::printf("\n");
  //   }
  // std::printf("\n");

}






void ndfes::PCurve::GetValue( double const t, double * vs ) const
{
  if ( mAkima )
    {
      for ( int k=0; k<mNumDim; ++k )
	{
	  vs[k] = mASpl[k].GetValue(t);
	}
    }
  else
    {
      for ( int k=0; k<mNumDim; ++k )
	{
	  vs[k] = mLSpl[k].GetValue(t);
	}
    }
}




void ndfes::PCurve::GetDeriv( double const t, double * vs ) const
{
  if ( mAkima )
    {
      for ( int k=0; k<mNumDim; ++k )
	{
	  vs[k] = mASpl[k].GetDeriv(t);
	}
    }
  else
    {
      for ( int k=0; k<mNumDim; ++k )
	{
	  vs[k] = mLSpl[k].GetDeriv(t);
	}
    }
}


namespace
{
   double distmetric( int ndim, double const * a, double const * b )
   {
       double d=0;
       for ( int i=0; i<ndim; ++i )
       {
          double e = a[i]-b[i];
	  d += e*e;
       }
       return std::sqrt(d);
   }
}

std::vector<double> ndfes::PCurve::GetPointClosestTo( double const * pt ) const
{
   int imin = 0;
   double dmin = 1.e+30;
   for ( int i=0; i<mNumPts; ++i )
   {
      double d = distmetric(mNumDim,pt,mX.data()+i*mNumDim);
      if ( d < dmin )
      {
         dmin=d;
	 imin=i;
      }
   }
   double tlo = mT[imin];
   double thi = mT[imin];
   if ( imin > 0 )
   {
      tlo = mT[imin-1];
   }
   if ( imin+1 < mNumPts )
   {
      thi = mT[imin+1];
   }
   double tmin = tlo;
   dmin=1.e+30;
   int ntry=201;
   double dt = (thi-tlo)/(ntry-1);
   for ( int itry=0; itry<ntry; ++itry )
   {
      double t = tlo + itry * dt;
      std::vector<double> crd( mNumDim, 0 );
      GetValue(t,crd.data());
      double d = distmetric(mNumDim,crd.data(),pt);
      //std::printf("t %12.5f %20.10e %20.10e\n",t,d,d-dmin);
      if ( d < dmin )
      {
         dmin=d;
	 tmin = t;
      }
      else
      {
	      break;
      }
   }
   std::vector<double> crd( mNumDim, 0 );
   GetValue(tmin,crd.data());
   return crd;
}


