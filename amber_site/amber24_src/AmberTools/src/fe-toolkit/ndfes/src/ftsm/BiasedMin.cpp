#include <iostream>
#include <iomanip>
#include <cstdlib>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include <nlopt.hpp>

#include "BiasedMin.hpp"
#include "PeriodicUtils.hpp"


namespace ndfes
{
  class OutOfBounds : public std::exception
  {
  public:
    virtual ~OutOfBounds() {};
    virtual char const * what () const noexcept
    {
      return "ndfes::OutOfBounds";
    }
  };
}


namespace ndfes
{
  class OptData_t
  {
  public:

    OptData_t
    ( std::size_t const my_ndim,
      double const * my_rc,
      double const * my_fc,
      double const my_oobk,
      std::shared_ptr<ndfes::FES> my_fes,
      bool my_verbose )
      :
      ndim(my_ndim),
      rc(my_rc, my_rc+my_ndim),
      fc(my_fc, my_fc+my_ndim),
      oobk(my_oobk),
      fes(my_fes),
      verbose(my_verbose),
      cnt(0),
      minval(1.e+30),
      relerr(0),
      isoob(false)
    {}

    bool evaluate( double const * x, double & v, double * g ) const
    {
      return fes->GetBiasedValue( rc.data(), fc.data(), oobk, x,
				  v, g );
    }

    void FindOcc( double * x ) const;
    
    std::size_t ndim;
    std::vector<double> rc;
    std::vector<double> fc;
    double oobk;
    std::shared_ptr<ndfes::FES> fes;
    bool verbose;
    std::size_t cnt;
    double minval;
    double relerr;
    bool isoob;
    std::vector<double> lastgoodpt;
  };


  double min_objective( std::vector<double> const & f,
			std::vector<double> & g,
			void * p );
  
}





void ndfes::OptData_t::FindOcc( double * pt ) const
{

  if ( lastgoodpt.size() > 0 )
    {
      std::copy( lastgoodpt.begin(), lastgoodpt.end(), pt );
      return;
    }
  
  std::size_t const ndim = fes->mDimInfo.GetNumDims();
  std::size_t const nbins = fes->mBins.size();
  std::vector<double> dvec(ndim,0);
  std::size_t ibin = 0;
  double mindist = 1.e+30;
  for ( std::size_t i=0; i<nbins; ++i )
    {
      if ( fes->mBinSizes[i] == 0 )
      	{
      	  continue;
      	};
      double dist = 0;
      std::vector<double> mydvec(ndim,0);
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  double dx = fes->mBins[i].center[dim] - pt[dim];
	  if ( fes->mDimInfo.IsPeriodic(dim) )
	    {
	      dx = wrap(dx,360.);
	    }
	  mydvec[dim] = dx;
	  dist += dx*dx;
	}
      if ( dist < mindist )
	{
	  mindist = dist;
	  ibin = i;
	  dvec = mydvec;
	}
    }

  double const * widths = fes->mDimInfo.GetTargetWidths();
  std::vector<double> mins( fes->mBins[ibin].center );
  std::vector<double> maxs( mins );
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      mins[dim] -= 0.5*widths[dim];
      maxs[dim] += 0.5*widths[dim];
    }

  
  std::size_t maxit = 101;
  double dt = 1./(maxit-1.);

  std::vector<double> tmp(pt,pt+ndim);

  for ( std::size_t it=0; it<maxit; ++it )
    {
      double t = it*dt;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  tmp[dim] = pt[dim] + t * dvec[dim];
	  //std::printf("%9.4f %9.4f %i %i\n",tmp[dim],fes->mBins[ibin].center[dim],tmp[dim]>mins[dim],tmp[dim]<maxs[dim]);
	}
      bool ok = true;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  if ( tmp[dim] < mins[dim] or tmp[dim] > maxs[dim] )
	    {
	      ok = false;
	      break;
	    }
	}
      if ( ok )
	{
	  //std::printf("ok %12.5f\n",t);
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      pt[dim] = tmp[dim];
	    }
	  break;
	}


	 
	
    };
}


double ndfes::min_objective
( std::vector<double> const & f,
  std::vector<double> & g,
  void * p )
{
  ndfes::OptData_t * data = (ndfes::OptData_t *)(p);
  std::size_t n = f.size();

  std::fill(g.begin(),g.end(),0);
  
  if ( n != data->ndim )
    {
      std::cerr << "ndfes::min_objective size mismatch "
		<< n << " " << data->ndim << std::endl;
    }
  
  double chisq = 0;
  std::vector<double> myg(n,0);
  
  bool ok =  data->evaluate( f.data(), chisq, myg.data() );

  // std::printf("%5lu ",data->cnt);
  // for ( std::size_t i=0; i<n; ++i )
  //   {
  //     std::printf("%12.6f",f[i]);
  //   }
  // std::printf("%13.4e",chisq);
  // for ( std::size_t i=0; i<n; ++i )
  //   {
  //     std::printf("%12.6f",myg[i]);
  //   }
  // std::printf(" %i\n",ok);
  
  
  if ( not ok )
    {
      data->isoob = true;
      throw ndfes::OutOfBounds();
    }
  else
    {
      data->isoob = false;
      data->lastgoodpt = f;
    }
  
  if ( g.size() > 0 )
    {
      std::copy( myg.begin(), myg.end(), g.begin() );
    }

  ++(data->cnt);
  double best = data->minval;
  double rel = (chisq-best)/std::abs(best);
  if ( data->cnt == 1 )
    {
      rel = 0.;
    }
  data->relerr = rel;
      
  if ( data->verbose )
    {
      std::cout << std::setw(5)
		<< data->cnt
		<< std::scientific << std::setw(20) << std::setprecision(10)
		<< chisq
		<< std::scientific << std::setw(14) << std::setprecision(2)
                << rel;

      for ( std::size_t i=0; i<n; ++i )
	{
	  std::cout << std::scientific << std::setw(13) << std::setprecision(4) << f[i];
	}
      
      if ( chisq < data->minval )
	{
	  std::cout << " *";
	  data->minval = chisq;
	}
      std::cout << "\n";
    }
  else if ( chisq < data->minval )
    {
      data->minval = chisq;
    }
  
  return chisq;
}




void ndfes::BiasedMins
( std::size_t const ndim,
  std::size_t const nsim,
  double const * rcs,
  double const * fcs,
  double * obsmeans,
  double const oobk,
  std::shared_ptr<ndfes::FES> fes,
  std::size_t const maxit,
  double const TOL,
  std::vector<double> const & minbounds,
  std::vector<double> const & maxbounds )
{
  bool verbose = false;
  
  std::copy( rcs, rcs+ndim*nsim, obsmeans );
  std::vector<double> ircs(rcs,rcs+ndim*nsim);
  
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      ndfes::OptData_t data( ndim, rcs+isim*ndim, fcs+isim*ndim,
			     oobk, fes, verbose );

      std::shared_ptr<nlopt::opt> Optimizer
	( new nlopt::opt( nlopt::LD_LBFGS, ndim ) );
      //( new nlopt::opt( nlopt::LN_COBYLA, ndim ) );

      Optimizer->set_min_objective( &ndfes::min_objective, &data );
      Optimizer->set_ftol_abs(1.e-13);
      Optimizer->set_ftol_rel(TOL);
      Optimizer->set_maxeval(maxit);

      std::vector<double> x(rcs+isim*ndim, rcs+(isim+1)*ndim);
      
      std::vector<double> xlo(ndim,0);
      std::vector<double> xhi(ndim,0);

      double const * ws = fes->mDimInfo.GetTargetWidths();
      int nlayers = std::max( 1, fes->NumLayersOfOcc( x.data() ) );
      //int nlayers = 1;
      //std::printf("isim nlayers %6lu %3i\n",isim,nlayers);
      
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  xlo[dim] = fes->mDimInfo.GetXmin()[dim];
	  xhi[dim] = fes->mDimInfo.GetXmax()[dim];
	  if ( nlayers > 0 )
	    {
	      xlo[dim] = std::max(xlo[dim],x[dim]-ws[dim]*nlayers);
	      xhi[dim] = std::min(xhi[dim],x[dim]+ws[dim]*nlayers);
	    }
	}
      if ( minbounds.size() == ndim*nsim )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      xlo[dim] = std::max(xlo[dim],minbounds[dim+isim*ndim]);
	    }
	}
      if ( maxbounds.size() == ndim*nsim )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      xhi[dim] = std::min(xhi[dim],maxbounds[dim+isim*ndim]);
	    }
	}
  

      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  x[dim] = std::max(std::min(x[dim],xhi[dim]),xlo[dim]);
	}

      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  ircs[dim+isim*ndim] = x[dim];
	}
      
      Optimizer->set_lower_bounds(xlo);
      Optimizer->set_upper_bounds(xhi);


      double chisq = 1.e+100;
      nlopt::result res = nlopt::SUCCESS;

      std::vector<double> x0(x);
      
      try
	{
	  res = Optimizer->optimize( x, chisq );
	}
      catch( nlopt::roundoff_limited & e )
	{
	  std::cout << e.what() << "\n";
	  std::cout << "Optimization halted because roundoff error limited progress\n";
	}
      catch( nlopt::forced_stop & e )
	{
	  std::cout << e.what() << "\n";
	  std::cout << "Optimization halted due to calling nlopt::opt::force_stop()\n";
	}
      
      catch( std::runtime_error & e )
	{
	  if ( ! data.isoob )
	    {
	      if ( data.relerr > 1.e-7 )
		{
		  std::cout << e.what() << "\n";
		  std::cout << "Encountered generic runtime exception\n";
		  std::cout << "relerr = "
			    << std::scientific << std::setw(12)
			    << std::setprecision(3)
			    << data.relerr << std::endl;
		};
	    }
	  else
	    {
	      data.FindOcc( x0.data() );
	      x = x0;
	    }
	}
      catch( std::bad_alloc & e )
	{
	  std::cout << e.what() << "\n";
	  std::cout << "Memory allocation error\n";
	}
      catch(...)
	{
	  std::cout << "NLopt threw an unknown exception. Aborting optimization.\n";
	}

      switch(res)
	{
	case nlopt::SUCCESS:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization completed successfully\n";
	      }
	    break;
	  }
	case nlopt::STOPVAL_REACHED:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization stopped because 'stopval' was reached\n";
	      }
	    break;
	  }
	case nlopt::FTOL_REACHED:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization completed successfully because 'ftol_rel' or 'ftol_abs' was reached\n";
	      }
	    break;
	  }
	case nlopt::XTOL_REACHED:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization completed successfully because 'xtol_rel' or 'xtol_abs' was reached\n";
	      }
	    break;
	  }
	case nlopt::MAXEVAL_REACHED:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization stopped because 'maxeval' was reached\n";
	      }
	    break;
	  }
	case nlopt::MAXTIME_REACHED:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization stopped because 'maxtime' was reached\n";
	      }
	    break;
	  }
	case nlopt::FAILURE:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization failed; generic failure\n";
	      }
	    break;
	  }
	case nlopt::INVALID_ARGS:
	  {
	    std::cout << "Optimization failed; invalid arguments (check bounds or unknown algorithm)\n";
	    break;
	  }
	case nlopt::OUT_OF_MEMORY:
	  {
	    std::cout << "Optimization failed; out of memory\n";
	    break;
	  }
	case nlopt::ROUNDOFF_LIMITED:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization failed; roundoff errors limited progress\n";
	      }
	    break;
	  }
	case nlopt::FORCED_STOP:
	  {
	    if ( verbose )
	      {
		std::cout << "Optimization failed; forced termination was issued\n";
	      }
	    break;
	  }
	default:
	  {
	    std::cout << "An unknown code was returned " << res << std::endl;
	    break;
	  }
	}

      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  obsmeans[dim+isim*ndim] = x[dim];
	}
      
    } // isim

  // for ( std::size_t i=0; i<nsim; ++i )
  //   {
  //     std::printf("%4lu",i);
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf("%13.4e",rcs[dim+i*ndim]);
  // 	}
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf("%13.4e",ircs[dim+i*ndim]);
  // 	}
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf("%13.4e",obsmeans[dim+i*ndim]);
  // 	}
  //     std::printf("\n");
  //   }
  
}
