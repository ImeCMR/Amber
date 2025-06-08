#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include "StatsMod.hpp"

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
  //double Ttest( int N1, double mu1, double var1,
 //		int N2, double mu2, double var2 );

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
  // double NormalDistInvCDF( double p );


  double KSTestProb( double alam );

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



bool ccdl::HaveTheSameMean
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



bool ccdl::HaveTheSameVar
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




void ccdl::CptMeanAndStderr
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
  //std::printf("pre  %8i %13.4e %13.4e",N,mu,err);
  if ( N > 1 )
    {
      double f = N;
      err = std::sqrt( err/( f*(f-1) ));
    };
  //std::printf("post %8i %13.4e %13.4e",N,mu,err);
}



void ccdl::CptMeanAndVar
( int const N,
  double const * A,
  double & mu,
  double & var )
{
  mu = 0.;
  for ( int i=0; i<N; ++i )
    {
      mu += A[i];
    };
  mu /= N;
  var = 0.;
  for ( int i=0; i<N; ++i )
    {
      double x = A[i]-mu;
      var += x*x;
    }
  if ( N > 1 )
    {
      var /= (N-1.);
    };
}


double ccdl::KSTestProb( double alam )
{
  double const EPS1 = 0.001;
  double const EPS2 = 1.e-8;
  double a2 = -2.*alam*alam;
  double fac = 2.;
  double probks = 0.;
  double termbf = 0.;
  for ( int j=0; j<100; ++j )
    {
      double term = fac*std::exp(a2*(j+1)*(j+1));
      probks += term;
      if ( std::abs(term) <= EPS1*termbf or std::abs(term) <= EPS2*probks )
	{
	  return probks;
	}
      fac = -fac;
      termbf = std::abs(term);
    }
  return 1.;
}


void ccdl::TwoSampleKSTest
( int n1, double const * x1,
  int n2, double const * x2,
  double & ksstatistic,
  double & pvalue )
{
  ksstatistic = 0.;
  std::vector<double> data1(x1,x1+n1);
  std::vector<double> data2(x2,x2+n2);
  std::sort( data1.begin(), data1.end() );
  std::sort( data2.begin(), data2.end() );
  double en1=n1;
  double en2=n2;
  int j1=0;
  int j2=0;
  double fn1=0.;
  double fn2=0.;
  while ( j1 < n1 and j2 < n2 )
    {
      double d1 = data1[j1];
      double d2 = data2[j2];
      if ( d1 <= d2 )
	{
	  fn1 = (j1+1.)/en1;
	  ++j1;
	}
      if ( d2 <= d1 )
	{
	  fn2 = (j2+1.)/en2;
	  ++j2;
	}
      double dt = std::abs(fn2-fn1);
      ksstatistic = std::max(ksstatistic,dt);
    }
  double en = std::sqrt( en1*en2/(en1+en2));
  pvalue = KSTestProb( (en+0.12+0.11/en)*ksstatistic );
}



bool ccdl::HaveTheSameDist
( int n1, double const * x1,
  int n2, double const * x2,
  double tol )
{
  double ks=0;
  double pvalue=1;
  TwoSampleKSTest( n1, x1, n2, x2, ks, pvalue );
  return pvalue > tol;
}




int ccdl::CptStatIneff( int const N, double const * A )
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


int ccdl::CptSampleStride( int const N, double const * A )
{
  std::vector<double> x(A,A+N);
  int g = ccdl::CptStatIneff(x.size(),x.data());
  int s = std::max(g/2,1);
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
	  
	  g = ccdl::CptStatIneff(x.size(),x.data());
	  if ( g > 1 )
	    {
	      s += std::max(g/2,1);
	    }
	};
    }
  return s;
}


