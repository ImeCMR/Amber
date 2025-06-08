#include <nlopt.hpp>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <memory>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "Opt.hpp"
#include "Edge.hpp"


namespace edgembar
{
  
  double ObjFcn( std::vector<double> const & xfree,
		 std::vector<double> & gfree,
		 void * objdata );

  
  class Objective
  {
  public:
    
    Objective( edgembar::Edge const * e, int verb );

    int curit;
    double bestchisq;
    int verbosity;
    edgembar::Edge const * edge;
    std::vector< edgembar::Trial const * > trials;
  };


}




double edgembar::ObjFcn
( std::vector<double> const & xfree,
  std::vector<double> & gfree,
  void * objdata )
{
  double chisq = 0.;

  edgembar::Objective * P = (edgembar::Objective *)(objdata);
  edgembar::Edge const * edge = P->edge;
  std::vector<double> x( edge->GetParamsFromFree(xfree) );
  int const n = x.size();

  /*
   std::printf("Opt.cpp/ObjFcn\n");
   std::printf("Xfree\n");
   for ( int i=0, nfree=xfree.size(); i<nfree; ++i )
     {
       std::printf("%5i %12.3e\n",i,xfree[i]);
     }
  
   std::printf("X\n");
   for ( int i=0; i<n; ++i )
     {
       std::printf("%5i %12.3e\n",i,x[i]);
     }
  */

  
  std::vector<double> g( n, 0. );
  int const nt = P->trials.size();
  
#ifdef WITH_OPENMP
  //if ( omp_in_parallel() )
  if ( false )
    {
#pragma omp parallel
      {
	std::vector<double> myg( x.size(), 0. );
	double mychisq = 0.;
#pragma omp for
	for ( int i=0; i<nt; ++i )
	{
	  edgembar::Trial const * trial = P->trials[i];
	  int o = trial->GetParamOffset();
	  mychisq += trial->CptObjective( x.data() + o, myg.data() + o );
	} // omp for
#pragma omp critical
	{
	  chisq += mychisq;
	  for ( int i=0; i<n; ++i )
	    {
	      g[i] += myg[i];
	    }
	} // omp critical
      } // omp parallel
    }
  else // omp_in_parallel
#endif
    {
      for ( int i=0; i<nt; ++i )
	{
	  edgembar::Trial const * trial = P->trials[i];
	  int o = trial->GetParamOffset();
	  chisq += trial->CptObjective( x.data() + o, g.data() + o );
	}
    }
  gfree = edge->GetFreeFromParams(g);

  // for ( int i=0, ng=g.size(); i<ng; ++i )
  //   {
  //     std::printf("%4i %10.4f %12.3e",i,x[i],g[i]);
  //     if ( i < (int)gfree.size() )
  // 	{
  // 	  std::printf("  %10.4f %12.3e",xfree[i],gfree[i]);
  // 	}
  //     std::printf("\n");
  //   }
  
  // -------------------------------------------------------------------
  // Optional printing
  // -------------------------------------------------------------------
  if ( P->verbosity > 0 )
  //if ( true )
    {
      if ( P->curit > 0 )
	{
	  double rel = (chisq-P->bestchisq) / std::abs(P->bestchisq);
	  std::printf("mbar iter: %6i chisq: %19.10e rel: %13.4e",
		      P->curit,chisq,rel);
	  if ( chisq < P->bestchisq )
	    {
	      std::printf(" *");
	    }
	}
      else
	{
	  std::printf("mbar iter: %6i chisq: %19.10e",
		      P->curit,chisq);
	}
      std::printf("\n");
    }
  // -------------------------------------------------------------------

  if ( chisq < P->bestchisq )
    {
      P->bestchisq = chisq;
    };
  P->curit += 1;
  
  return chisq;
}



edgembar::Objective::Objective
( edgembar::Edge const * e,
  int verb )
  : curit(0),
    bestchisq(1.e+100),
    verbosity(verb),
    edge(e),
    trials( e->GetConstTrials() )
{
  
}



double edgembar::Optimize
( edgembar::Edge const * e,
  std::vector<double> const & pguess,
  std::vector<double> & popt,
  double const tol,
  int const verbosity )
{
  double chisq = 1.e+100;
  int nparam = e->GetNumParam();
  int nfree = e->GetNumFreeParam();
 
  if ( (int)pguess.size() != nparam )
    {
      std::cerr << "Error in edgembar::Optimize number of guess "
		<< "parameters does not match system size : "
		<< pguess.size() << " " << nparam << std::endl;
      std::exit(EXIT_FAILURE);
    }

  popt.assign(nparam, 0.);

  edgembar::Objective objdata( e, verbosity );
  
  std::vector<double> pfree( e->GetFreeFromParams(pguess) );

  if ( nfree == 0 )
  {
     std::vector<double> gfree;
     chisq = edgembar::ObjFcn( pfree, gfree, &objdata );
     popt = e->GetParamsFromFree( pfree );
     return chisq;
  }


  std::shared_ptr<nlopt::opt> Optimizer
    ( new nlopt::opt( nlopt::LD_LBFGS, nfree ) );



  // std::printf("Opt.cpp/Optimize\n");
  // std::printf("pguess\n");
  // for ( int i=0, ni=pguess.size(); i<ni; ++i )
  //   {
  //     std::printf("%5i %13.4e\n",i,pguess[i]);
  //   }
  // std::printf("pfree\n");
  // for ( int i=0, ni=pfree.size(); i<ni; ++i )
  //   {
  //     std::printf("%5i %13.4e\n",i,pfree[i]);
  //   }
  

  
  Optimizer->set_min_objective( &edgembar::ObjFcn, &objdata );
  Optimizer->set_ftol_abs(1.e-13);
  Optimizer->set_ftol_rel(tol);
  Optimizer->set_maxeval(10000);

  nlopt::result res = nlopt::SUCCESS;

  bool verbose = (verbosity > 0);
  
  try
    {
      res = Optimizer->optimize( pfree, chisq );
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

      // check the gradients. If they are all small, then it's probably fine
      std::vector<double> gfree(pfree.size(),0);
      edgembar::ObjFcn( pfree, gfree, &objdata );
      double maxg = -1.;
      for ( int i=0, nf=gfree.size(); i<nf; ++i )
	{
	  maxg = std::max(maxg,std::abs(gfree[i]));
	}
      if ( maxg > 1.e-5 )
	{
	  std::cout << e.what() << "\n";
	  std::cout << "LD_LBFGS Encountered generic runtime exception. Restarting optimization with LN_COBYLA\n";
	  
	  Optimizer.reset( new nlopt::opt( nlopt::LN_COBYLA, nfree ) );
	  Optimizer->set_min_objective( &edgembar::ObjFcn, &objdata );
	  Optimizer->set_ftol_abs(1.e-13);
	  Optimizer->set_ftol_rel(tol);
	  Optimizer->set_maxeval(10000);
	  bool hasexcept=false;
	  try
	    {
	      res = Optimizer->optimize( pfree, chisq );
	    }
	  catch(...)
	    {
	      hasexcept=true;
	      std::cout << "Failed a second time with LN_COBYLA\n";
	    }
	  if ( ! hasexcept )
	    {
	      std::cout << "LN_COBYLA completed without throwing an exception\n";
	    }
	};
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
	  };
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
	  };
        break;
      }
    case nlopt::MAXTIME_REACHED:
      {
	if ( verbose )
	  {
	    std::cout << "Optimization stopped because 'maxtime' was reached\n";
	  };
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
	if ( verbose )
	  {
	    std::cout << "Optimization failed; invalid arguments (check bounds or unknown algorithm)\n";
	  };
        break;
      }
    case nlopt::OUT_OF_MEMORY:
      {
	if ( verbose )
	  {
	    std::cout << "Optimization failed; out of memory\n";
	  };
        break;
      }
    case nlopt::ROUNDOFF_LIMITED:
      {
	if ( verbose )
	  {
	    std::cout << "Optimization failed; roundoff errors limited progress\n";
	  };
        break;
      }
    case nlopt::FORCED_STOP:
      {
	if ( verbose )
	  {
	    std::cout << "Optimization failed; forced termination was issued\n";
	  };
        break;
      }
    default:
      {
	if ( verbose )
	  {
	    std::cout << "An unknown code was returned " << res << std::endl;
	  };
        break;
      }
    }

  popt = e->GetParamsFromFree( pfree );
  return chisq;
}
