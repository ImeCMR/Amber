#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "StatsMod.hpp"

#include "Stats.hpp"


/*
namespace ccdl
{
  
  void arth
  ( double const first, 
    double const increment,
    int const n,
    double * rVec );

  double gammln( double const xx );
  
  double beta_cfrac( double a, double b, double x );
  
  double incbeta( double a, double b, double x );

  // reject if small
  double Ttest( double t, double dof );
  
  // reject if small
  double Ttest( int N1, double mu1, double var1,
		int N2, double mu2, double var2 );

  // reject if small
  double Ftest( double F, double dof1, double dof2 );
  
  // reject if small
  double Ftest( int N1, double var1,int N2, double var2 );


  // utility function for NormalDistInvCDF
  double NormalDistInvCDF_RationalApprox( double t );

  // Inverse cumulative distribution function (aka the probit function)
  // https://www.johndcook.com/blog/cpp_phi_inverse/
  // https://www.johndcook.com/blog/normal_cdf_inverse/
  // To get the error bars for a 0.95% confidence interval, call
  // NormalDistInvCDF( 0.5*(0.95+1.0) )
  // The result is 1.96, which you'd scale the standard error by
  // The following call produces "1.0000"
  // NormalDistInvCDF( 0.5*(0.68270473138+1.0) )
  double NormalDistInvCDF( double p );
  
}


void ccdl::arth
( double const first, 
  double const increment,
  int const n,
  double * rVec )
{
  int const NPAR_ARTH = 16;
  int const NPAR2_ARTH = 8;
  
  if ( n > 0 )
    {
      rVec[0] = first;
    }
  
  if ( n <= NPAR_ARTH )
    {
      for ( int k=1; k < n; ++k )
	{
	  rVec[k] = rVec[k-1] + increment;
	}
    }
  else
    {
      for ( int k=1; k < NPAR2_ARTH; ++k )
	{
	  rVec[k] = rVec[k-1] + increment;
	}
      double temp = increment * NPAR2_ARTH;
      int k = NPAR2_ARTH;
      
      while (1)
	{
	  if ( k >= n-1 )
	    {
	      return;
	    }
	  else
	    {
	      int k2 = k + k;
	      int nk = n-k;
	      int m  = std::min(k,nk);
	      for ( int i=0; i<m; ++i )
		{
		  rVec[k+i] = temp + rVec[i];
		}
	      temp = temp + temp;
	      k = k2;
	    };
	};
    };
  
}



double ccdl::gammln( double const xx )
{
  double const stp = 2.5066282746310005;
  double const coef[] = { 76.18009172947146, -86.50532032941677,
			  24.01409824083091, -1.231739572450155,
			  0.1208650973866179e-2, -0.5395239384953e-5 };
  
  if ( xx <= 0.0 )
    {
      std::cerr << "ccdl::gammln xx <= 0.0" << std::endl;
    }
  
  double const x = xx;
  double tmp = x + 5.5;
  tmp = (x+0.5)*std::log(tmp)-tmp;
  
  double arth[6];
  
  ccdl::arth(x+1.0,1.0,6,arth);
  double sum = 0.0;
  for ( int i=0; i<6; ++i )
    {
      sum += coef[i] / arth[i];
    }
  
  double g = tmp + std::log( stp * ( 1.000000000190015 + sum ) / x );
  return g; 
}




double ccdl::beta_cfrac( double a, double b, double x )
{
  int const MAXIT = 100;
  double const EPS = 3.e-7;
  double const FPMIN = 1.e-30;
  
  double qab = a+b;
  double qap = a+1.;
  double qam = a-1.;
  double c = 1.;
  double d = 1.-qab*x/qap;
  if ( std::abs(d) < FPMIN )
    {
      d = FPMIN;
    }
  d = 1./d;
  double h=d;
  for ( int m=1; m<MAXIT+1; ++m )
    {
      double m2=2*m;
      double aa=m*(b-m)*x/((qam+m2)*(a+m2));
      d = 1.+aa*d;
      if ( std::abs(d) < FPMIN )
	{
	  d = FPMIN;
	}
      c = 1.+aa/c;
      if ( std::abs(c) < FPMIN )
	{
	  c = FPMIN;
	}
      d = 1./d;
      h *= d*c;
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      d = 1.+aa*d;
      if ( std::abs(d) < FPMIN )
	{
	  d = FPMIN;
	}
      c = 1.+aa/c;
      if ( std::abs(c) < FPMIN )
	{
	  c = FPMIN;
	}
      d = 1./d;
      double del=d*c;
      h *= del;
      if ( std::abs(del-1.) < EPS )
	{
	  break;
	};
    }
  return h;
}



double ccdl::incbeta( double a, double b, double x )
{
  double result = 0.;
  double bt = 0.;
  x = std::max(0.,std::min(1.,x));
  if ( x > 0. and x < 1. )
    {
      bt = std::exp( ccdl::gammln(a+b)
		     - ccdl::gammln(a)
		     - ccdl::gammln(b)
		     + a*std::log(x)
		     + b*std::log(1.-x) );
    }
  if ( x < (a+1.)/(a+b+2.) )
    {
      result = bt*ccdl::beta_cfrac(a,b,x)/a;
    }
  else
    {
      result = 1.-bt*ccdl::beta_cfrac(b,a,1.-x)/b;
    }
  return result;
}


double ccdl::Ttest( double t, double dof )
{
  return ccdl::incbeta( 0.5*dof, 0.5, dof/(dof+t*t) );
}

double ccdl::Ttest
( int N1, double mu1, double var1,
  int N2, double mu2, double var2 )
{
  double vn1 = var1/N1;
  double vn2 = var2/N2;
  double t = (mu1-mu2)/std::sqrt(vn1+vn2);
  double n1m1 = N1 > 1 ? N1-1 : 1;
  double n2m1 = N2 > 1 ? N2-1 : 1;
  double dof = (std::pow(vn1+vn2,2)) / ( vn1*vn1/n1m1 + vn2*vn2/n2m1 );
  return ccdl::Ttest( t, dof );
}


double ccdl::Ftest( double F, double dof1, double dof2 )
{
  double p = 2. * ccdl::incbeta( 0.5*dof2, 0.5*dof1, dof2/(dof2+dof1*F) );
  if ( p > 1. )
    {
      p = 2.-p;
    }
  return p;
}



double ccdl::Ftest
( int N1, double var1,
  int N2, double var2 )
{
  double f = 0;
  double df1 = 0;
  double df2 = 0;
  if ( var1 > var2 )
    {
      f = var1/var2;
      df1 = N1-1;
      df2 = N2-1;
    }
  else
    {
      f = var2/var1;
      df1 = N2-1;
      df2 = N1-1;
    }
  return ccdl::Ftest( f, df1, df2 );
}



bool edgembar::HaveTheSameMean
( int N1, double mu1, double serr1,
  int N2, double mu2, double serr2,
  double tol )
{
  bool res = true;
  double err = std::sqrt( serr1*serr1 + serr2*serr2 );
  if ( err > 1.e-8 )
    {
      double var1 = serr1*serr1*N1;
      double var2 = serr2*serr2*N2;
      double p = ccdl::Ttest(N1,mu1,var1,N2,mu2,var2);
      res = p > tol;
    }
  return res;
}



bool edgembar::HaveTheSameVar
( int N1, double serr1,
  int N2, double serr2,
  double tol )
{
  bool res = true;
  double err = std::sqrt( serr1*serr1 + serr2*serr2 );
  if ( err > 1.e-8 )
    {
      double var1 = serr1*serr1*N1;
      double var2 = serr2*serr2*N2;
      double p = ccdl::Ftest(N1,var1,N2,var2);
      res = p > tol;
    }
  return res;
}




double ccdl::NormalDistInvCDF_RationalApprox(double t)
{
  // Abramowitz and Stegun formula 26.2.23.
  // The absolute value of the error should be less than 4.5 e-4.
  double c[] = {2.515517, 0.802853, 0.010328};
  double d[] = {1.432788, 0.189269, 0.001308};
  return t - ((c[2]*t + c[1])*t + c[0]) / 
    (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}


double ccdl::NormalDistInvCDF( double p )
{
  if ( p <= 0 )
    {
      p = 1.e-10;
    }
  else if ( p >= 1 )
    {
      p = 1 - 1.e-10;
    }
  
  // See article above for explanation of this section.
  if (p < 0.5)
    {
      // F^-1(p) = - G^-1(p)
      return -ccdl::NormalDistInvCDF_RationalApprox
	( std::sqrt(-2.0*std::log(p)) );
    }
  else
    {
      // F^-1(p) = G^-1(1-p)
      return ccdl::NormalDistInvCDF_RationalApprox
	( std::sqrt(-2.0*std::log(1-p)) );
    }
}




void edgembar::CptMeanAndStderr
( int const N,
  double const * A,
  double & mu,
  double & err )
{
  mu = 0.;
  for ( int i=0; i<N; ++i )
    {
      mu += A[i];
    };
  mu /= N;
  err = 0.;
  for ( int i=0; i<N; ++i )
    {
      double x = A[i]-mu;
      err += x*x;
    }
  if ( N > 1 )
    {
      err = std::sqrt( err/( N*(N-1) ));
    };
}



// void edgembar::CptMeanAndStderr
// ( int const istart,
//   int const stride,
//   int const M,
//   double const * A,
//   double & mu,
//   double & err )
// {
//   mu  = 0.;
//   err = 0.;
//   if ( stride < 1 )
//     {
//       std::cerr << "edgembar::CptMeanAndStderr illegal stride value "
// 		<< stride << std::endl;
//       std::exit(EXIT_FAILURE);
//     }
//   int N = 0;
//   for ( int i=istart; i<M; i += stride )
//     {
//       N += 1;
//       mu += A[i];
//     }
//   if ( N > 0 )
//     {
//       mu /= N;
//       for ( int i=istart; i<M; i += stride )
// 	{
// 	  double x = A[i]-mu;
// 	  err += x*x;
// 	}
//       if ( N > 1 )
// 	{
// 	  err = std::sqrt( err/( N*(N-1) ));
// 	};
//     };
// }



int edgembar::CptStatIneff( int const N, double const * A )
{
  int const mintime = 3;
  
  
  //std::printf("dats[0:1] %.8f %.8f\n",A[0],A[1]);
  
  double g = 1.0;
  
  double mu = 0.;
  for ( int i=0; i<N; ++i )
    mu += A[i];
  mu /= N;
  
  std::vector<double> dA(A,A+N);
  for ( int i=0; i<N; ++i )
    dA[i] -= mu;
  
  double sigma2_AB = 0.;
  for ( int i=0; i<N; ++i )
    sigma2_AB += dA[i]*dA[i];
  sigma2_AB /= N;
  
  //std::printf("N,sig2 %i %.8f\n",N,sigma2_AB);
  
  if ( std::abs(sigma2_AB) > 1.e-8 )
    {
      for ( int t=1; t < N-1; t += 1 )
        {
          double C = 0.;
          for ( int j=0; j < N-t; ++j )
            C += dA[j]*dA[j+t];
          C /= ( sigma2_AB * ( N-t ) );

          if ( C <= 0.0 and t > mintime )
            break;

          g += 2.0 * C * (1.0 - ((double)t)/((double)N));
        }
    }
  //std::printf("g=%.6f sigma2_AB=%.6f\n",g,sigma2_AB);
  if ( g < 1.0 )
    g = 1.0;
  int blocksize = g+0.5;
  return blocksize;
}


int edgembar::CptSampleStride( int const N, double const * A )
{
  std::vector<double> x(A,A+N);
  int g = CptStatIneff(x.size(),x.data());
  int s = g;
  while ( g > 1 )
    {
      if ( s >= N/2 )
	{
	  break;
	}
      else
	{	  
	  x.resize(0);
	  for ( int i=0; i<N; i += s )
	    {
	      x.push_back( A[i] );
	    };
	  
	  g = CptStatIneff(x.size(),x.data());
	  if ( g > 1 )
	    {
	      s += std::max(g,1);
	    }
	};
    }
  return s;
}


*/



edgembar::TimeseriesSegment::TimeseriesSegment()
  : istart(0),
    istop(0),
    s(1),
    cn(0),
    cavg(0),
    cvar(0),
    cstd(0),
    cerr(0),
    un(0),
    uavg(0),
    uvar(0),
    ustd(0),
    uerr(0)
{}


edgembar::TimeseriesSegment::TimeseriesSegment
( int mistart,
  int mistop,
  int mn,
  double const * mxs )
{
  reset( mistart, mistop, mn, mxs );
}

void edgembar::TimeseriesSegment::reset
( int mistart,
  int mistop,
  int nn,
  double const * xs )
{
  istart = mistart;
  istop = std::min(mistop,nn);
  cdata.assign( xs + istart, xs + istop );
  s = ccdl::CptSampleStride(cdata.size(),cdata.data());
  udata.resize(0);
  for ( int i=istart; i<istop; i += s )
    {
      udata.push_back( xs[i] );
    }
  
  un = udata.size();
  ccdl::CptMeanAndVar( udata.size(), udata.data(), uavg, uvar );
  ustd=0.;
  uerr=0.;
  if ( uvar > 0 )
    {
      ustd = std::sqrt(uvar);
      if ( un > 0 )
	{
	  uerr = std::sqrt(uvar/un);
	}
    }
 
  cn = cdata.size();
  ccdl::CptMeanAndVar( cdata.size(), cdata.data(), cavg, cvar );
  cstd=0.;
  cerr=0.;
  if ( cvar > 0 )
    {
      cstd = std::sqrt(cvar);
      if ( un > 0 )
	{
	  cerr = std::sqrt(cvar/un);
	}
    };
}


bool edgembar::TimeseriesSegment::IsSimilarTo
( edgembar::TimeseriesSegment const & rhs,
  double const ptol ) const
{
  bool res = false;
  if ( ptol > 0 )
    {
      res = ccdl::HaveTheSameMean(un,cavg,cerr,rhs.un,rhs.cavg,rhs.cerr,ptol);
    }
  else
    {
      res = std::abs(cavg-rhs.cavg) < std::sqrt( cerr*cerr + rhs.cerr*rhs.cerr + 1.e-8 );
    }

  return res;
}

/*
edgembar::TimeseriesHalves::TimeseriesHalves()
  : istart(0),
    istop(0),
    s(1),
    n(0),
    n1(0),
    n2(0),
    avg1(0),
    err1(0),
    avg2(0),
    err2(0)
{}


edgembar::TimeseriesHalves::TimeseriesHalves
( int mistart,
  int mistop,
  int mn,
  double const * mxs )
{
  reset( mistart, mistop, mn, mxs );
}


void edgembar::TimeseriesHalves::reset
( int mistart,
  int mistop,
  int nn,
  double const * xs )
{
  istart = mistart;
  istop = std::min(mistop,nn);

  std::vector<double> x( xs + istart, xs + istop );
  s = ccdl::CptSampleStride(x.size(),x.data());

  int imid = x.size()/2;
  std::vector<double> a,b;
  for ( int i=0, nx=x.size(); i<nx; i += s )
    {
      if ( i < imid )
	{
	  a.push_back( x[i] );
	}
      else
	{
	  b.push_back( x[i] );
	};
    }
  data1 = a;
  data2 = b;
  n = a.size() + b.size();
  n1 = a.size();
  avg1 = 0.;
  err1 = 0.;
  ccdl::CptMeanAndStderr( a.size(), a.data(), avg1, err1 );
  n2 = b.size();
  avg2 = 0.;
  err2 = 0.;
  ccdl::CptMeanAndStderr( b.size(), b.data(), avg2, err2 );
}


bool edgembar::TimeseriesHalves::HasSameMeans( double const ptol ) const
{
  //bool close = false;
  //bool ttest = false;
  bool res = false;
  if ( ptol > 0 )
    {
      res = ccdl::HaveTheSameMean(n1,avg1,err1,n2,avg2,err2,std::abs(ptol));
    }
  else
    {
      res = std::abs(avg1-avg2) < std::sqrt( err1*err1 + err2*err2 + 1.e-8 );
    };
  return res;
}


bool edgembar::TimeseriesHalves::HasSameMeansAndDist( double const ptol ) const
{
  bool res = false;
  if ( ptol <= 0 )
    {
      res = HasSameMeans(ptol);
    }
  else
    {
      bool same_means = ccdl::HaveTheSameMean(n1,avg1,err1,
					      n2,avg2,err2,
					      ptol );
      bool same_dist = ccdl::HaveTheSameDist(n1,data1.data(),
					     n2,data2.data(),
					     ptol );
      res = same_means and same_dist;
    }
  return res;
}


bool edgembar::TimeseriesHalves::HasSameDist( double const ptol ) const
{
  bool res = false;
  if ( ptol <= 0 )
    {
      res = HasSameMeans(ptol);
    }
  else
    {
      //bool same_means = ccdl::HaveTheSameMean(n1,avg1,err1,
      //n2,avg2,err2,
      //				      ptol );
      bool same_dist = ccdl::HaveTheSameDist(n1,data1.data(),
					     n2,data2.data(),
					     ptol );
      res = same_dist;
    }
  return res;
}
*/







/*

bool edgembar::CptSampleStartAndStride
( int const N,
  double const * x,
  double const maxeq,
  int & istart,
  int & stride )
{
  int maxn = 0;
  istart = 0;
  stride = 1;

  bool isconv = false;
  
  int iskip = 1;

  int ncyc=0;

  int i0max = std::max(1,std::min(N,(int)(maxeq*N+0.5)));
  
  for ( int i0=0; i0 < i0max; i0 += std::max(1,iskip) )
    {
      
      ncyc += 1;


      std::vector<double> xprod(x+i0,x+N);
      int s = edgembar::CptSampleStride(xprod.size(),xprod.data());
      
      std::vector<double> xlo;
      std::vector<double> xhi;
      int kmid = (i0+N)/2;
      for ( int k=i0; k < N; k += s )
	{
	  if ( k < kmid )
	    {
	      xlo.push_back(x[k]);
	    }
	  else
	    {
	      xhi.push_back(x[k]);
	    };
	}
      
      int n = xlo.size() + xhi.size();

      double mulo = 0.;
      double muhi = 0.;
      double errlo = 0.;
      double errhi = 0.;
      edgembar::CptMeanAndStderr(xlo.size(),xlo.data(),mulo,errlo);
      edgembar::CptMeanAndStderr(xhi.size(),xhi.data(),muhi,errhi);
      
      //std::printf("%5i %4i %5i %5i   (%4lu) %12.4f +- %12.4f   (%4lu) %12.4f +- %12.4f\n",i0,s,n,maxn,xlo.size(),mulo,errlo,xhi.size(),muhi,errhi);

      
      if ( std::abs( mulo - muhi ) < std::sqrt( errlo*errlo + errhi*errhi ) + 1.e-7 )
	{
	  isconv=true;
	  
	  if ( n > maxn )
	    {
	      maxn = n;
	      istart = i0;
	      stride = s;
	      iskip = std::max(1,stride/2);
	      if ( stride < 2 )
		{
		  break;
		}
	    }
	  else if ( i0 > istart + 2*stride )
	    {
	      //std::printf("Early termination at %i cycle %i\n",i0,ncyc);
	      break;
	    }

	}
      else if ( (! isconv) and (i0 + std::max(1,iskip) >= i0max) )
	{
	  istart = std::max(0,i0max);
	  std::vector<double> xprod(x+istart,x+N);
	  stride = edgembar::CptSampleStride(xprod.size(),xprod.data());
	  break;
	}
      else if ( N > 100 and i0 + iskip + 1 < i0max-1 )
	{
	  iskip++;
	};

    }
  return isconv;
}




bool edgembar::CptSampleStartAndStrideByBlock
( int const N,
  double const * x,
  double const maxeq,
  int const Nblocks,
  int & istart,
  int & stride )
{
  int maxn = 0;
  istart = 0;
  stride = 1;

  bool isconv = false;
  
  int iskip = std::max(1,N/Nblocks);

  int ncyc=0;

  int i0max = std::max(1,std::min(N,(int)(maxeq*N+0.5)));
  //std::printf("i0max N %i %i\n",i0max,N);
  for ( int i0=0; i0 < i0max; i0 += std::max(1,iskip) )
    {
      //std::printf("i0: %i\n",i0);
      // if ( i0 >= i0max-1 or i0+std::max(1,iskip) >= i0max )
      // 	{
      // 	  istart = i0max-1;
      // 	  isconv = false;
      // 	  //std::printf("Setting not conv\n");
      // 	  break;
      // 	}
      ncyc += 1;


      std::vector<double> xprod(x+i0,x+N);
      int s = edgembar::CptSampleStride(xprod.size(),xprod.data());
      
      std::vector<double> xlo;
      std::vector<double> xhi;
      int kmid = (i0+N)/2;
      for ( int k=i0; k < N; k += s )
	{
	  if ( k < kmid )
	    {
	      xlo.push_back(x[k]);
	    }
	  else
	    {
	      xhi.push_back(x[k]);
	    };
	}

      int n = xlo.size() + xhi.size();
	
      
      double mulo = 0.;
      double muhi = 0.;
      double errlo = 0.;
      double errhi = 0.;
      edgembar::CptMeanAndStderr(xlo.size(),xlo.data(),mulo,errlo);
      edgembar::CptMeanAndStderr(xhi.size(),xhi.data(),muhi,errhi);
      
      //std::printf("%5i %4i %5i %5i   (%4lu) %12.4f +- %12.4f   (%4lu) %12.4f +- %12.4f\n",i0,s,n,maxn,xlo.size(),mulo,errlo,xhi.size(),muhi,errhi);

      
      if ( std::abs( mulo - muhi ) < std::sqrt( errlo*errlo + errhi*errhi ) + 1.e-7 )
	{
	  isconv=true;
	  
	  if ( n > maxn )
	    {
	      maxn = n;
	      istart = i0;
	      stride = s;
	      if ( stride < 2 )
		{
		  break;
		}
	    }
	  else if ( i0 > istart + 2*iskip )
	    {
	      //std::printf("Early termination at %i cycle %i %i\n",i0,ncyc,isconv);
	      break;
	    }

	}
      else if ( (! isconv) and (i0 + std::max(1,iskip) >= i0max) )
	{
	  istart = std::max(0,i0max);
	  std::vector<double> xprod(x+istart,x+N);
	  stride = edgembar::CptSampleStride(xprod.size(),xprod.data());
	  break;
	}
	

    }

  // if ( ! isconv )
  //   {
  //     std::printf("not conv\n");
  //   }
  
  return isconv;
}
*/







bool edgembar::CptSampleStartAndStrideByBlock_v3
( int const N,
  double const * x,
  double const maxeq,
  int const Nblocks,
  int & istart,
  int & stride,
  double const ptol )
{
  
  bool conv = false;
  istart = 0;
  stride = 1;
  //int imax = std::max(1,std::min(N-2,(int)(maxeq*N+0.5)));
  int imax = std::max(1,std::min(N,(int)(maxeq*N)));

  
  std::vector<int> blkSizes(Nblocks,0);
  for ( int i=0; i<N; ++i )
    {
      ++blkSizes[ i % Nblocks ];
    };
  std::vector<int> blkStart(Nblocks,0);
  std::vector<int> blkStop(Nblocks,0);
  for ( int i=0; i<Nblocks; ++i )
    {
      if ( i > 0 )
	{
	  blkStart[i] = blkStop[i-1];
	}
      blkStop[i] = blkStart[i] + blkSizes[i];
    }

  // qblk is the first block past the first quarter
  int qblk = (Nblocks*0.25 + 0.5);
  // mblk is the first block past the midway point
  int mblk = (Nblocks*0.5 + 0.5);
  // fblk is the first block past the last allowable block
  int fblk = (Nblocks*maxeq + 0.5);
  
  edgembar::TimeseriesSegment sec( blkStart[mblk], N, N, x );

  int ntarget = sec.un;
  
  // find the first proposed production block by walking
  // backwards from the midpoint and performing a KS test
  // to see if the additional samples appear to come from
  // the same distribution as the second half.
  int pblk = mblk;
  for ( int iblk=mblk; iblk --> 0; )
    {
      edgembar::TimeseriesSegment seg( blkStart[iblk], blkStop[mblk-1], N, x );
      
      bool same_dist = ccdl::HaveTheSameDist( sec.cn, sec.cdata.data(),
					      seg.cn, seg.cdata.data(),
					      ptol );
      
      bool same_mean = ccdl::HaveTheSameMean( sec.un, sec.cavg, sec.cerr,
					      seg.un, seg.cavg, seg.cerr,
					      ptol );

      bool match = false;
      if ( iblk >= qblk )
	{
	  match = same_dist or same_mean;
	}
      else
	{
	  match = same_dist and same_mean;
	}

      bool morepts = false;
      int segsec_n = 0;
      if ( match )
	{
	  edgembar::TimeseriesSegment segsec( blkStart[iblk], N, N, x );
	  segsec_n = segsec.un;
	  morepts = segsec_n >= ntarget;
	}

      // double dist_p = 1.;
      // double ks = 0.;
      // ccdl::TwoSampleKSTest( sec.cdata.size(), sec.cdata.data(),
      // 			     seg.cdata.size(), seg.cdata.data(),
      // 			     ks, dist_p );

      // double var1 = sec.cerr*sec.cerr*sec.un;
      // double var2 = seg.cerr*seg.cerr*seg.un;
      
      // double mean_p = ccdl::Ttest(sec.un,sec.cavg,var1,
      // 				  seg.un,seg.cavg,var2);

	
      // std::printf("%3i ",iblk);
      // std::printf("%5i %6i %4i ",blkStart[iblk],blkStop[mblk-1],seg.s);
      // std::printf("%6i ",seg.un);
      // std::printf("%9.3f %7.3f ",seg.cavg,seg.cerr);
      // std::printf("%5i %6i %4i ",blkStart[mblk],N,sec.s);
      // std::printf("%6i ",sec.un);
      // std::printf("%9.3f %7.3f ",sec.cavg,sec.cerr);
      // std::printf("%7.3f %7.3f ",mean_p,dist_p);
      // std::printf("%2i %2i\n",match,morepts);

      
      if ( match and morepts )
	{
	  pblk = iblk;
	  ntarget = segsec_n;
	}
      else if ( iblk < qblk )
	{
	  break;
	}
    };


  //std::printf("Final pblk: %i\n",pblk);

  for ( int iblk=pblk; iblk < Nblocks; ++iblk )
    {

      //edgembar::TimeseriesHalves h( blkStart[iblk], N, N, x );
      int imid = (blkStart[iblk] + N)/2 + 1;
      edgembar::TimeseriesSegment fhalf( blkStart[iblk], imid, N, x );
      edgembar::TimeseriesSegment lhalf( imid, N, N, x );
	
      bool same_mean = ccdl::HaveTheSameMean( fhalf.un, fhalf.cavg, fhalf.cerr,
					      lhalf.un, lhalf.cavg, lhalf.cerr,
					      ptol );


      //std::printf("iblk %3i %5i %3i (cur best size %4i)\n",iblk,h.n,h.s,size);
      //std::printf("    %12.3f +- %5.3f  %12.3f +- %5.3f %i\n",
      //	  h.avg1,h.err1,h.avg2,h.err2,(int)(h.IsSelfSimilar()));

      //std::printf("halves iblk: %2i %i\n",iblk,same_mean);
      
      if ( same_mean ) // or h.HasSameMeansAndDist(ptol/2) )
	{
	  conv = true;
	  //size = h.n;
	  istart = blkStart[iblk];
	  edgembar::TimeseriesSegment prod( blkStart[iblk], N, N, x );
	  stride = prod.s;
	  break;
	}

      if ( iblk == fblk )
	{
	  break;
	};
      
    }
  //std::printf("conv? %i\n",conv);
  
  if ( ! conv )
    {
      edgembar::TimeseriesSegment full( blkStart[pblk], N, N, x );
      edgembar::TimeseriesSegment last( imax, N, N, x );
      
      bool same_mean = ccdl::HaveTheSameMean( full.un, full.cavg, full.cerr,
					      last.un, last.cavg, last.cerr,
					      ptol );


      // double var1 = full.cerr*full.cerr*full.un;
      // double var2 = last.cerr*last.cerr*last.un;
     
      // double mean_p = ccdl::Ttest(full.un,full.cavg,var1,
      // 				  last.un,last.cavg,var2);
      
      // std::printf("not conv same mean? %i %7.3f %i %i %i\n",same_mean,
      // 		  mean_p,full.un,last.un,imax);
      
      // std::printf("%13.4e +- %13.4e  %13.4e %13.4e\n",
      // 		  full.cavg,full.cerr, last.cavg, last.cerr );

      
      if ( same_mean )
	{
	  istart = full.istart;
	  stride = full.s;
	  conv = true;
	}
      else
	{
	  istart = last.istart;
	  stride = last.s;
	}
    }
  //std::printf("final result: %i %i %i\n",istart,stride,conv);

  return conv;
}
