#ifndef _StatsMod_hpp_
#define _StatsMod_hpp_

#include <vector>

namespace ccdl
{

  // Inverse cumulative distribution function (aka the probit function)
  // https://www.johndcook.com/blog/cpp_phi_inverse/
  // https://www.johndcook.com/blog/normal_cdf_inverse/
  // To get the error bars for a 0.95% confidence interval, call
  // NormalDistInvCDF( 0.5*(0.95+1.0) )
  // The result is 1.96, which you'd scale the standard error.
  //
  // The following call produces "1.0000". This is a 2-tailed 68% confidence
  // NormalDistInvCDF( 0.5*(0.68270473138+1.0) )
  
  double NormalDistInvCDF( double p );

  
  void CptMeanAndStderr( int const N, double const * A,
			 double & mu, double & err );

  void CptMeanAndVar( int const N, double const * A,
		      double & mu, double & var );

  
  int CptStatIneff( int const n, double const * x );

  
  int CptSampleStride( int const n, double const * x );


  bool HaveTheSameMean( int N1, double mu1, double serr1,
			int N2, double mu2, double serr2,
			double tol );

  bool HaveTheSameVar( int N1, double serr1,
		       int N2, double serr2,
		       double tol );
  
  double Ttest( int N1, double mu1, double var1,
                int N2, double mu2, double var2 );
  
  
  void TwoSampleKSTest( int n1, double const * x1,
			int n2, double const * x2,
			double & ksstatistic,
			double & pvalue );

  bool HaveTheSameDist( int n1, double const * x1,
			int n2, double const * x2,
			double tol );
}

#endif

