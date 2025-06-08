#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Convergence.hpp"

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

  double Ttest( double t, double dof );
  
  double Ttest
  ( double N1, double mu1, double var1,
    double N2, double mu2, double var2 );

  double Ttest_ind_from_stats( double N1, double mu1, double stddev1,
			       double N2, double mu2, double stddev2 );
  
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
( double N1, double mu1, double var1,
  double N2, double mu2, double var2 )
{
  double vn1 = var1/N1;
  double vn2 = var2/N2;
  double t = (mu1-mu2)/std::sqrt(vn1+vn2);
  double n1m1 = N1 > 1 ? N1-1 : 1;
  double n2m1 = N2 > 1 ? N2-1 : 1;
  double dof = (std::pow(vn1+vn2,2)) / ( vn1*vn1/n1m1 + vn2*vn2/n2m1 );
  return ccdl::Ttest( t, dof );
}


double ccdl::Ttest_ind_from_stats( double N1, double mu1, double stddev1,
				   double N2, double mu2, double stddev2 )
{
  return ccdl::Ttest(N1,mu1,stddev1*stddev1,N2,mu2,stddev2*stddev2);
}
  



void ndfes::FTSMCheckSameMeans
( ndfes::PathIter const & opath,
  ndfes::PathIter const & cpath,
  std::vector<double> const & deffc,
  double const beta,
  double const ptol,
  bool & conv,
  double & pmin,
  bool const verbose )
{
  conv = true;
  pmin = 2;

  std::size_t ndim = opath.mPath.mNumDim;
  std::size_t nsim = opath.mPath.mNumPts;

  if ( cpath.mPath.mNumPts != nsim )
    {
      conv = false;
      return;
    }
  
  //std::vector<double> cps( cpath.GetPathProjection( cpath ) );
  std::vector<double> pvals(ndim,0);
  std::vector<double> diff(ndim,0);

  // for ( std::size_t i=0; i<nsim; ++i )
  //   {
  //     std::printf("%2lu",i);
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf("%10.5f",opath.rcs[dim+i*ndim]);
  // 	}
  //      for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf("%10.5f",cpath.rcs[dim+i*ndim]);
  // 	}
  //      std::printf("\n");
  //   }

  int const nobs = 100;
  std::vector<double> stds(ndim,0);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      stds[dim] = 1. / ( 2*beta*deffc[dim] );
    }
  
  for ( std::size_t i=0; i<nsim; ++i )
    {
      std::string allconv = "T";
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  pvals[dim] = ccdl::Ttest_ind_from_stats
	    ( nobs,
	      opath.mPath.mPts[dim+i*ndim],
	      stds[dim],
	      nobs,
	      cpath.mPath.mPts[dim+i*ndim],
	      stds[dim] );
	  
	  diff[dim] = std::abs( opath.mPath.mPts[dim+i*ndim] -
				cpath.mPath.mPts[dim+i*ndim] );
	  
	  if ( pvals[dim] < pmin )
	    {
	      pmin = pvals[dim];
	    }
	  if ( pvals[dim] < ptol )
	    {
	      conv = false;
	      allconv="F";
	    };
	}

#define FMTF(w,p) std::fixed << std::setprecision((p)) << std::setw((w))
#define FMTE(w,p) std::scientific << std::setprecision((p)) << std::setw((w))

      if ( verbose )
	{
	  std::cout << "img " << std::setw(3) << i+1
		    << " " << allconv;
      
	  //<< " cos: " << FMTF(5,2) << cps[i];
      
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      std::string tf = "T";
	      if ( pvals[dim] < ptol ) { tf="F"; };
	      std::cout << " [d:" << FMTE(8,1) << diff[dim]
			<< " p: " << FMTF(5,3) << pvals[dim]
			<< ">" << FMTF(5,3) << ptol
			<< " ? " << tf << "]";
	    }

	  std::cout << std::endl;
	}
      
#undef FMTF
#undef FMTE
      
    }
  
}






double ndfes::GetSlope( int mlen, int nin, double const * vs )
{
  double M = 1.e+30;
  
  int n = std::min(mlen,nin);
  if ( n > 1 )
    {
      int i0 = nin - n;
      double xavg = 0;
      double yavg = 0;
      for ( int i=i0; i<nin; ++i )
	{
	  xavg += i;
	  yavg += vs[i];
	}
      xavg /= n;
      yavg /= n;
      double den = 0;
      double num = 0;
      for ( int i=i0; i<nin; ++i )
	{
	  double dx = i-xavg;
	  double dy = vs[i]-yavg;
	  den += dx*dx;
	  num += dx*dy;
	}
      M = num/den;
    }
  return M;
}







void ndfes::OptDisplacements
( std::shared_ptr<ndfes::PCurve> pspl,
  std::shared_ptr<ndfes::PCurve> ospl,
  std::size_t ndim,
  std::size_t nsim,
  std::vector<double> & maxdisps,
  std::vector<double> & avgdisps )
{
  //std::size_t nsim = paths.back()->rcs.size() / ndim;
   maxdisps.assign(ndim,0);
   avgdisps.assign(ndim,0);

   //std::shared_ptr<ndfes::PCurve> pspl = paths.back()->pspl;
   //std::shared_ptr<ndfes::PCurve> ospl = paths.back()->GetOutputSpline( popts );

   double dt = 1. / ( nsim - 1. );
   for ( std::size_t i=0; i<nsim; ++i )
   {
	std::vector<double> a(ndim,0);
	std::vector<double> b(ndim,0);
	pspl->GetValue( i*dt, a.data() );
	ospl->GetValue( i*dt, b.data() );
        for ( std::size_t dim=0; dim<ndim; ++dim )
	{
           double dx = std::abs(a[dim] - b[dim]);
           avgdisps[dim] += dx / nsim;
	   maxdisps[dim] = std::max(maxdisps[dim],dx);
	}
   }

}








double ndfes::GetMinSlope( int mlen, int nin, double const * vs )
{
   double Mmin = 1.e+30;
   for ( int len = nin; len --> mlen; )
     {
       double M = ndfes::GetSlope(len,nin,vs);
	Mmin = std::min(M,Mmin);
     }
   return Mmin;
}




void ndfes::FTSMSlopeTest
( ndfes::PathOpt const & simpaths,
  //ndfes::PathOptions const & popts,
  int const mlen,
  double const distol,
  double const angtol,
  double const rmsdtol,
  std::vector<bool> const tidxisangle,
  double & distolgood,
  double & angtolgood,
  double & rmsdgood,
  bool & maxdispisbig,
  bool & avgdispisbig,
  std::vector< std::ostream * > couts )
{
  
  distolgood = 0;
  angtolgood = 0;
  rmsdgood = 1.e+30;
  maxdispisbig = false;
  avgdispisbig = false;

  std::size_t npath = simpaths.mIters.size();
  std::size_t path0 = std::max(0,(int)(npath)-mlen);
  std::size_t ndata = npath - path0;
  
  std::size_t nsim = simpaths.mIters.back().mPath.mNumPts;
  std::size_t ndim = simpaths.mIters.back().mPath.mNumDim;
  double dt = 1. / (nsim - 1.);


  if ( npath == 1 )
    {
      distolgood = 1.e+30;
      angtolgood = 1.e+30;
      rmsdgood = 1.e+30;
      maxdispisbig = true;
      avgdispisbig = true;
    }
  

  std::size_t ncout = couts.size();
  
  for ( std::size_t iout=0; iout<ncout; ++iout )
    {
      std::ostream & fh = *(couts[iout]);
      for ( std::size_t cnt=0; cnt<80; ++cnt )
	{
	  fh << "=";
	};
      fh << "\n"
	 << "SLOPE TEST\n\n";
    }
      

  std::vector< std::shared_ptr<ndfes::PCurve> > ospls;
  for ( std::size_t j=0; j<npath; ++j )
    {
      ospls.push_back( simpaths.mIters[j].mPath.mSpl );
    }
  
  
  std::vector<double> data( ndata * ndim * nsim, 0 );

  for ( std::size_t i=0; i<nsim; ++i )
    {
      double t = i * dt;
      for ( std::size_t j=0; j<ndata; ++j )
	{
	  std::size_t ipath = j + path0;
	  std::vector<double> c(ndim,0);
	  //simpaths[ipath]->pspl->GetValue( t, c.data() );
	  
	  ospls[ipath]->GetValue( t, c.data() );
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      data[j+(dim+i*ndim)*ndata] = c[dim];
	    }
	}
    }

  std::vector<double> ms( ndim * nsim, 0 );
  for ( std::size_t i=0; i<nsim; ++i )
    {
      
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  ms[dim+i*ndim]
	    = ndfes::GetSlope( ndata, ndata, data.data() + (dim+i*ndim)*ndata );
	}

      std::string conv = "T";
      std::vector<std::string> convs;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  double m = std::abs(ms[dim+i*ndim]);
	  std::string myconv = "T";
	  if ( tidxisangle[dim] )
	    {
	      if ( m > angtol )
		{
		  conv   = "F";
		  myconv = "F";
		}
	      angtolgood = std::max(m,angtolgood);
	    }
	  else
	    {
	      if ( m > distol )
		{
		  conv   = "F";
		  myconv = "F";
		}
	      distolgood = std::max(m,distolgood);
	    }
	  convs.push_back( myconv );
	}
      

      if ( npath > 1 )
	{
      
	  for ( std::size_t iout=0; iout<ncout; ++iout )
	    {
	      std::ostream & fh = *(couts[iout]);
	      fh << std::setw(3) << i+1 << " " << conv << " ";
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  double tol = distol;
		  if ( tidxisangle[dim] )
		    {
		      tol = angtol;
		    }
		  fh << "|m=" << std::fixed << std::setw(7)
		     << std::setprecision(4) << ms[dim+i*ndim]
		     << "|<" << std::fixed << std::setw(6)
		     << std::setprecision(4) << tol
		     << " ? " << convs[dim] << "]";
		  if ( dim+1 < ndim )
		    {
		      fh << " ";
		    }
		}
	      fh << "\n";
	    }

	}
    }


  if ( (int)npath == 1 )
    {
      for ( std::size_t iout=0; iout<ncout; ++iout )
	{
	  std::ostream & fh = *(couts[iout]);
	  fh << "No slope to test because it is the initial iteration\n";
	}
    }
  
  for ( std::size_t iout=0; iout<ncout; ++iout )
    {
      std::ostream & fh = *(couts[iout]);
      fh << "Convergence would be met if distol >= "
	 << std::scientific << std::setw(11) << std::setprecision(2)
	 << distolgood
	 << " and angtol >= "
	 << std::scientific << std::setw(11) << std::setprecision(2)
	 << angtolgood
	 << "\n";
      if ( distolgood < distol and angtolgood < angtol )
	{
	  fh << "Slopes appear to be converged\n\n";
	}
      else
	{
	  fh << "Slopes are NOT converged\n\n";
	}
    }
  

  //////////////////////////////////////////////////////////////
  
  for ( std::size_t iout=0; iout<ncout; ++iout )
    {
      std::ostream & fh = *(couts[iout]);
      for ( std::size_t cnt=0; cnt<80; ++cnt )
	{
	  fh << "=";
	};
      fh << "\n"
	 << "RMSD TEST\n\n";
    }
  
  std::vector<double> Ds( npath, 0 );

  data.assign( ndim*nsim, 0 );
  for ( std::size_t i=0; i<nsim; ++i )
    {
      simpaths.mIters[0].mPath.GetValue(i*dt,data.data()+i*ndim);
    }
  
  for ( std::size_t ipath=0; ipath < npath; ++ipath )
    {
      double D = 0.;
      for ( std::size_t i=0; i<nsim; ++i )
	{
	  std::vector<double> c(ndim,0);
	  //simpaths[ipath]->pspl->GetValue(i*dt,c.data());
	  ospls[ipath]->GetValue(i*dt,c.data());
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      double dx = data[dim+i*ndim] - c[dim];
	      D += dx*dx;
	    }
	}
      Ds[ipath] = std::sqrt( D / (nsim*ndim) );
    }

  
  for ( std::size_t iout=0; iout<ncout; ++iout )
    {
      std::ostream & fh = *(couts[iout]);
      fh << "\n"
	 << "Check Path RMSD relative to initial guess\n"
	 << std::setw(5) << "Iter" << " " << std::setw(12) << "RMSD" << "\n";
      
      for ( std::size_t ipath=0; ipath<npath; ++ipath )
	{
	  // std::string dname = "init";
	  // if ( ipath > 0 )
	  //   {
	  //     std::stringstream ss;
	  //     ss << "it" << std::setfill('0') << std::setw(3) << ipath;
	  //     dname = ss.str();
	  //   }
	  fh << std::setw(5) << ipath << " "
	     << std::fixed << std::setw(12) << std::setprecision(4)
	     << Ds[ipath] << "\n";
	};
    };

  
  if ( (int)npath >= mlen )
    {
      rmsdgood = ndfes::GetMinSlope( mlen, npath, Ds.data() );

      for ( std::size_t iout=0; iout<ncout; ++iout )
	{
	  std::ostream & fh = *(couts[iout]);
	  fh << "\n"
	     << "\nMin slope in RMSD "
	     << std::scientific << std::setw(12) << std::setprecision(3)
	     << rmsdgood;
	  if ( rmsdgood < rmsdtol )
	    {
	      fh << " Path is effectively CONVERGED\n\n";
	    }
	  else
	    {
	      fh << " Path is NOT converged\n\n";
	    }
	};
    }
  else
    {
      for ( std::size_t iout=0; iout<ncout; ++iout )
	{
	  std::ostream & fh = *(couts[iout]);
	  fh << "\n"
	     << "RMSD slope not checked because there are fewer than "
	     << mlen << " paths.\n\n";
	}
    }
  
  
  //////////////////////////////////////////////////////////////
  
  for ( std::size_t iout=0; iout<ncout; ++iout )
    {
      std::ostream & fh = *(couts[iout]);
      for ( std::size_t cnt=0; cnt<80; ++cnt )
	{
	  fh << "=";
	};
      fh << "\n"
	 << "OPTIMIZATION DISPLACEMENT TEST\n\n";
    }


  
  std::vector<double> maxdisps(ndim,0);
  std::vector<double> avgdisps(ndim,0);
  //ndfes::OptDisplacements( simpaths, popts, ndim, maxdisps, avgdisps );

  if ( simpaths.mIters.size() == 1 )
    {
      ndfes::OptDisplacements( simpaths.mIters[0].mPath.mSpl, ospls[0], ndim, nsim, maxdisps, avgdisps );
    }
  else
    {
      int ilast = (int)ospls.size() - 1;
      ndfes::OptDisplacements( ospls[ilast-1], ospls[ilast], ndim, nsim, maxdisps, avgdisps );
    }
  
  for ( std::size_t icout=0; icout<ncout; ++icout )
    {
      std::ostream & fh = *(couts[icout]);
      
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  double tol = distol;
	  if ( tidxisangle[dim] )
	    {
	      tol = angtol;
	    }
	  fh << "Max disp in dim " << std::setw(2) << dim + 1 << " "
	     << std::scientific << std::setw(12) << std::setprecision(3) 
	     << maxdisps[dim] << " < " << tol 
	     << " ? " << (maxdisps[dim] < tol) << "\n";
	}
  
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  double tol = distol;
	  if ( tidxisangle[dim] )
	    {
	      tol = angtol;
	    }
	  fh << "Avg disp in dim " << std::setw(2) << dim + 1 << " "
	     << std::scientific << std::setw(12) << std::setprecision(3)
	     << avgdisps[dim] << " < " << tol
	     << " ? " << (avgdisps[dim] < tol) << "\n";
	}
    }

  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      double tol = distol;
      if ( tidxisangle[dim] )
	{
	  tol = angtol;
	}

      if ( avgdisps[dim] > tol )
	{
	  avgdispisbig = true;
	}
      
      if ( maxdisps[dim] > tol )
	{
	  maxdispisbig = true;
	}
      
    }
  
}
    
