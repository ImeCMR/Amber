#include <iostream>
#include <cstdlib>
#include <cmath>

#include "StateValues.hpp"


void edgembar::StateValues::Assert
( std::string msg, std::vector<double> const & x ) const
{
  int n = x.size();
  if ( mResult.size() > 0 )
    {
      if ( (int)mResult.size() != n )
	{
	  std::cerr << "Error in " << msg
		    << "input array size does not match result size "
		    << x.size() << " " << mResult.size() << std::endl;
	  std::exit(EXIT_FAILURE);
	};
    }
  if ( mBootstrap.size() > 0 )
    {
      if ( (int)mBootstrap.size() != n )
	{
	  std::cerr << "Error in " << msg
		    << "input array size does not match bootstrap size "
		    << x.size() << " " << mBootstrap.size() << std::endl;
	  std::exit(EXIT_FAILURE);
	};
    }
}

void edgembar::StateValues::SetResult
( std::vector<double> const & x )
{
  Assert("edgembar::StateValues::SetResult",x);
  mResult = x;
  if ( (int)mBootstrap.size() == 0 )
    {
      mBootstrap.resize( x.size() );
      for ( int i=0, n=x.size(); i<n; ++i )
	{
	  mBootstrap.resize(0);
	}
    }
}

void edgembar::StateValues::PushBootstrap
( std::vector<double> const & x )
{
  Assert("edgembar::StateValues::PushBootstrap",x);
  if ( (int)mBootstrap.size() == 0 )
    {
      mBootstrap.resize( x.size() );
    }
  for ( int i=0, n=x.size(); i<n; ++i )
    {
      mBootstrap[i].push_back( x[i] );
    }
}

std::vector<double> edgembar::StateValues::GetStderr() const
{
  int n = mResult.size();
  
  std::vector<double> stderr(n,0.);

  if ( (int)mBootstrap.size() > 0 )
    {
      for ( int i=0; i<n; ++i )
	{
	  double mu = 0;
	  int nboot = mBootstrap[i].size();
	  for ( int k=0; k<nboot; ++k )
	    {
	      mu += mBootstrap[i][k];
	    }
	  mu /= nboot;
	  
	  double sq = 0.;
	  for ( int k=0; k<nboot; ++k )
	    {
	      double d = mBootstrap[i][k] - mu;
	      sq += d*d;
	    };
	  if ( nboot > 1 )
	    {
	      stderr[i] = std::sqrt( sq/(nboot-1) );
	    }
	  else
	    {
	      stderr[i] = 0.;
	    }
	};
    }
  
  return stderr;
}


void edgembar::StateValues::Join
( edgembar::StateValues const & o )
{
  if ( o.mResult.size() > 0 and
       mResult.size() > 0 )
    {
      std::cerr << "Error in edgembar::StateValues::Join both objects "
		<< "unepxectedly contain results (only bootstraps "
		<< "can be joined) " << o.mResult.size()
		<< " " << mResult.size() << std::endl;
      std::exit(EXIT_FAILURE);
    }
  else if ( (int)mResult.size() == 0 and
	    (int)o.mResult.size() > 0 )
    {
      mResult = o.mResult;
    }

  if ( (int)mBootstrap.size() == 0 and
       (int)o.mBootstrap.size() > 0 )
    {
      mBootstrap = o.mBootstrap;
      if ( (int)mResult.size() > 0 and
	   mBootstrap.size() != mResult.size() )
	{
	  std::cerr << "Error in edgembar::StateValues::Join bootstrap "
		    << "size doesn't match result "
		    << mBootstrap.size() << " " << mResult.size()
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}
    }
  else if ( (int)mBootstrap.size() > 0 and
       (int)o.mBootstrap.size() > 0 )
    {
      if ( (int)mBootstrap.size() != (int)o.mBootstrap.size() )
	{
	  std::cerr << "Error in edgembar::StateValues::Join objects "
		    << "have different bootstrap sizes "
		    << mBootstrap.size() << " " << o.mBootstrap.size()
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      else
	{
	  int n = mBootstrap.size();
	  for ( int i=0; i<n; ++i )
	    {
	      mBootstrap[i].reserve(mBootstrap[i].size()+o.mBootstrap.size());
	      mBootstrap[i].insert( mBootstrap[i].end(),
				    o.mBootstrap[i].begin(),
				    o.mBootstrap[i].end() );
	    }
	};
    }
  
}
