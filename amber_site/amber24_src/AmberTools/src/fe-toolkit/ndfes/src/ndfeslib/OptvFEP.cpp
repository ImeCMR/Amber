#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <nlopt.hpp>

#include "OptvFEP.hpp"


namespace ndfes
{
  class OptvFEP_t
  {
  public:

    OptvFEP_t( ndfes::SystemInfo const & a,
	       std::vector< std::vector<std::size_t> > const & b,
	       bool const verb );
    
  public:

    ndfes::SystemInfo const * sysinfo;
    std::vector<double> cornerh;
    bool verbose;
    double best_chisq;
    int iteration;
    
  };


  double OptvFEP_callback( std::vector<double> const & f,
			   std::vector<double> & g,
			   void * p );
  
}

ndfes::OptvFEP_t::OptvFEP_t
( ndfes::SystemInfo const & a,
  std::vector< std::vector<std::size_t> > const & sbinsamples,
  bool const verb )
  : sysinfo( &a ),
    cornerh( a.CptvFEPLinearTerm(sbinsamples) ),
    verbose(verb),
    best_chisq( 1.e+100 ),
    iteration( 0 )
{
  
}


double ndfes::OptvFEP_callback
( std::vector<double> const & f,
  std::vector<double> & g,
  void * p )
{
  
  ndfes::OptvFEP_t * data = (ndfes::OptvFEP_t *)(p);

  double chisq = data->sysinfo->CptvFEPObjective(f,g,data->cornerh);
 

  if ( data->verbose )
    {
      data->iteration++;
      double best = data->best_chisq;
      double rel = (chisq-best)/std::abs(best);
      if ( data->iteration == 1 )
	{
	  rel = 0.;
	}
      std::cout << std::setw(5)
		<< data->iteration
		<< std::scientific << std::setw(20) << std::setprecision(10)
		<< chisq
		<< std::scientific << std::setw(14) << std::setprecision(2)
		<< rel;
      if ( chisq < best )
	{
	  std::cout << " *";
	  data->best_chisq = chisq;
	}
      std::cout << "\n";
    }

  return chisq;
}


void ndfes::OptvFEP
( std::size_t const nparam,
  ndfes::SystemInfo const & sysinfo,
  std::vector< std::vector<std::size_t> > const & sbinsamples,
  double const TOL,
  int const maxit,
  bool const verbose,
  std::vector<double> & state_fs )
{
  
  double const DEL = 100.;


  
  ndfes::OptvFEP_t data(sysinfo,sbinsamples,verbose);

  
  std::shared_ptr<nlopt::opt> Optimizer
    ( new nlopt::opt( nlopt::LD_LBFGS, nparam ) );
    //( new nlopt::opt( nlopt::LN_COBYLA, nparam ) );

  Optimizer->set_min_objective( &ndfes::OptvFEP_callback, &data );
  Optimizer->set_ftol_abs(1.e-13);
  Optimizer->set_ftol_rel(TOL);
  Optimizer->set_maxeval(maxit);

  
  std::vector<double> x(nparam,0.);
  if ( (std::size_t)state_fs.size() == nparam )
    {
      x=state_fs;
    };

  
  {
    std::vector<double> xlo(x);
    for ( std::size_t i=0; i<nparam; ++i )
      {
	xlo[i] -= DEL;
      }
    Optimizer->set_lower_bounds(xlo);
  }

  
  {
    std::vector<double> xhi(x);
    for ( std::size_t i=0; i<nparam; ++i )
      {
	xhi[i] += DEL;
      }
    Optimizer->set_upper_bounds(xhi);
  }


  if ( verbose )
    {
      std::cout << "\n\n";
      std::cout << "Nonlinear optimization of the vFEP equations\n"
		<< "performed with the L-BFGS algorithm implemented in NLopt"
		<< "\n\n"
		<< "     Giese, T. J.; Ekesan, S.; York, D. M. Extension of the Variational\n"
		<< "  Free Energy Profile and Multistate Bennett Acceptance Ratio Methods\n"
		<< "  for High-Dimensional Potential of Mean Force Profile Analysis.\n"
		<< "  J. Phys. Chem. A 2021, 125, 4216-4232.\n"
		<< "     Lee, T.-S.; Radak, B. K.; Pabis, A.; York, D. M. A new maximum\n"
		<< "  likelihood approach for free energy profile construction from\n"
		<< "  molecular simulations. J. Chem. Theory Comput. 2013, 9, 153-164.\n"
		<< "     Lee, T.-S.; Radak, B. K.; Huang, M.; Wong, K.-Y.; York, D. M.\n"
		<< "  Roadmaps through free energy landscapes calculated using the\n"
		<< "  multidimensional vFEP approach. J. Chem. Theory Comput. 2014, 10,\n"
		<< "  24-34.\n"
		<< "     Johnson, S. G. The NLopt nonlinear-optimization package.\n"
		<< "  http://github.com/stevengj/nlopt.\n"
		<< "     Nocedal, J. Updating quasi-Newton matrices with limited storage.\n"
		<< "  Math. Comput. 1980, 35, 773-782.\n"
		<< "     Liu, D. C.; Nocedal, J. On the limited memory BFGS method for large\n"
		<< "  scale optimization. Math. Programming, 1989, 45, 503-528.\n"
		<< "\n";
      std::cout << std::setw(4) << "Iter"
		<< std::setw(20) << "Objective"
		<< std::setw(14) << "Rel Change"
		<< "\n";
    }

  double chisq = 1.e+100;
  nlopt::result res = nlopt::SUCCESS;

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
      std::cout << e.what() << "\n";
      std::cout << "Encountered generic runtime exception\n";
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
        std::cout << "Optimization stopped because 'stopval' was reached\n";
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
        std::cout << "Optimization stopped because 'maxeval' was reached\n";
        break;
      }
    case nlopt::MAXTIME_REACHED:
      {
        std::cout << "Optimization stopped because 'maxtime' was reached\n";
        break;
      }
    case nlopt::FAILURE:
      {
        std::cout << "Optimization failed; generic failure\n";
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
        std::cout << "Optimization failed; roundoff errors limited progress\n";
        break;
      }
    case nlopt::FORCED_STOP:
      {
        std::cout << "Optimization failed; forced termination was issued\n";
        break;
      }
    default:
      {
        std::cout << "An unknown code was returned " << res << std::endl;
        break;
      }
    }

  if ( verbose )
    {
      std::cout << "\n";
    }
  
  state_fs = x;
  for ( std::size_t i=0; i<nparam; ++i )
    {
      state_fs[i] -= x[0];
    };
}



void ndfes::OptvFEP_CheckGradients
( std::size_t const nparam,
  ndfes::SystemInfo const & sysinfo,
  std::vector< std::vector<std::size_t> > const & sbinsamples,
  std::vector<double> & state_fs )
{
  double const DEL = 5.e-5;
  
  std::vector<double> b(nparam,0.);
  
  if ( (std::size_t)state_fs.size() == nparam )
    {
      b = state_fs;
    }
  else
    {
      for ( std::size_t i=0; i<nparam; ++i )
	{
	  b[i] = std::cos( ((double)i)/nparam );
	};
    }

  std::vector<double> tmp(nparam,0.);
  std::vector<double> g(nparam,0.);
  std::vector<double> cornerh( sysinfo.CptvFEPLinearTerm(sbinsamples) );
  sysinfo.CptvFEPObjective(b,g,cornerh);

  std::cout << "\n\n"
	    << "Checking if CptvFEPObjective gradients are correct\n\n";
  std::cout << std::setw(4) << "idx"
	    << std::setw(15) << "analytic grd"
	    << std::setw(15) << "numerical grd"
	    << std::setw(15) << "difference"
	    << "\n";
  
  for ( std::size_t i=0; i<nparam; ++i )
    {
      b[i] += DEL;
      double hi = sysinfo.CptvFEPObjective(b,tmp,cornerh);
      //double hi = sysinfo.CptFastvFEPObjective(b,samples);
      b[i] -= 2*DEL;
      double lo = sysinfo.CptvFEPObjective(b,tmp,cornerh);
      //double lo = sysinfo.CptFastvFEPObjective(b,samples);
      b[i] += DEL;

      double num = (hi-lo)/(2*DEL);
      double diff = g[i]-num;
      
#define FMTE15 std::scientific << std::setw(15) << std::setprecision(3)
      std::cout << std::setw(4) << i
		<< FMTE15 << g[i]
		<< FMTE15 << num
		<< FMTE15 << diff
		<< "\n";
#undef FMTE15
    }
  std::cout << "\n";

}
