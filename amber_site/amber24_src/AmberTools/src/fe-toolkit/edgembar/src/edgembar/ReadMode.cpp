#include "ReadMode.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

edgembar::ReadMode edgembar::GetMode( std::string mode )
{
  std::transform( mode.begin(), mode.end(), mode.begin(), ::toupper );
  edgembar::ReadMode i = edgembar::AUTO;
  if ( mode.compare("AUTO") == 0 )
    {
      i = edgembar::AUTO; 
    }
  else if ( mode.compare("MBAR") == 0 )
    {
      i = edgembar::MBAR;
    }
  else if ( mode.compare("MBAREXP0") == 0 )
    {
      i = edgembar::MBAREXP0;
    }
  else if ( mode.compare("MBAREXP1") == 0 )
    {
      i = edgembar::MBAREXP1;
    }
  else if ( mode.compare("MBAREXP") == 0 )
    {
      i = edgembar::MBAREXP;
    }
  else if ( mode.compare("BAR") == 0 )
    {
      i = edgembar::BAR;
    }
  else if ( mode.compare("BAREXP0") == 0 )
    {
      i = edgembar::BAREXP0;
    }
  else if ( mode.compare("BAREXP1") == 0 )
    {
      i = edgembar::BAREXP1;
    }
  else if ( mode.compare("BAREXP") == 0 )
    {
      i = edgembar::BAREXP;
    }
  else
    {
      std::cerr << "Invalid mode '" << mode << "'\n";
      std::exit(EXIT_FAILURE);
    }
  return i;
}



std::string edgembar::GetModeString( edgembar::ReadMode mode )
{
  std::string name;
  switch ( mode )
    {
    case edgembar::INVALID : { name="INVALID"; break; }
    case edgembar::MBAR : { name="MBAR"; break; }
    case edgembar::MBAREXP0 : { name="MBAREXP0"; break; }
    case edgembar::MBAREXP1 : { name="MBAREXP1"; break; }
    case edgembar::MBAREXP : { name="MBAREXP"; break; }
    case edgembar::BAR : { name="BAR"; break; }
    case edgembar::BAREXP0 : { name="BAREXP0"; break; }
    case edgembar::BAREXP1 : { name="BAREXP1"; break; }
    case edgembar::BAREXP : { name="BAREXP"; break; }
    default:
      {
	std::cerr << "Programming error: Invalid mode in "
		  << "edgembar::GetModeString " << mode << "\n";
	std::exit(EXIT_FAILURE);
	break;
      }

    };
  return name;
}




bool edgembar::CanEvalMode
( int const ns,
  int const * avail,
  edgembar::ReadMode const mode )
{
  bool ok = true;

  switch ( mode )
    {
      
    case edgembar::AUTO:
      {
	break;
      }

    case edgembar::MBAR:
      {
	if ( ns < 2 )
	  {
	    ok = false;
	  }
	for ( int t=0; t<ns; ++t )
	  {
	    for ( int e=0; e<ns; ++e )
	      {
		if ( avail[e+t*ns] != 1 )
		  {
		    ok=false;
		  }
	      }
	  }
	break;
      }

      
    case edgembar::MBAREXP0:
      {
	if ( ns < 2 )
	  {
	    ok=false;
	  };
	for ( int t=1; t<ns; ++t )
	  {
	    for ( int e=1; e<ns; ++e )
	      {
		if ( avail[e+t*ns] != 1 )
		  {
		    ok=false;
		  }
	      }
	  }
	if ( avail[0+1*ns] != 1 )
	  {
	    ok = false;
	  }
	break;
      }

    case edgembar::MBAREXP1:
      {
	if ( ns < 2 )
	  {
	    ok=false;
	  };
	for ( int t=0; t<ns-1; ++t )
	  {
	    for ( int e=0; e<ns-1; ++e )
	      {
		if ( avail[e+t*ns] != 1 )
		  {
		    ok=false;
		  }
	      }
	  }
	if ( avail[(ns-1)+(ns-2)*ns] != 1 )
	  {
	    ok = false;
	  }
	break;
      }

    case edgembar::MBAREXP:
      {
	if ( ns < 3 )
	  {
	    ok=false;
	  };
	for ( int t=1; t<ns-1; ++t )
	  {
	    for ( int e=1; e<ns-1; ++e )
	      {
		if ( avail[e+t*ns] != 1 )
		  {
		    ok=false;
		  }
	      }
	  }
	if ( avail[0+1*ns] != 1 )
	  {
	    ok = false;
	  }
	if ( avail[(ns-1)+(ns-2)*ns] != 1 )
	  {
	    ok = false;
	  }
	break;
      }

   
    case edgembar::BAR:
      {
	if ( ns < 2 )
	  {
	    ok=false;
	  };
	for ( int t=0; t<ns; ++t )
	  {
	    for ( int e=std::max(t-1,0); e<std::min(t+2,ns); ++e )
	      {
		if ( avail[e+t*ns] != 1 )
		  {
		    ok=false;
		  }
	      }
	  }
	break;
      }

    case edgembar::BAREXP0:
      {
	if ( ns < 2 )
	  {
	    ok=false;
	  };
	for ( int t=1; t<ns; ++t )
	  {
	    for ( int e=std::max(t-1,1); e<std::min(t+2,ns); ++e )
	      {
		if ( avail[e+t*ns] != 1 )
		  {
		    ok=false;
		  }
	      }
	  }
	if ( avail[0+1*ns] != 1 )
	  {
	    ok = false;
	  }
	break;
      }

    case edgembar::BAREXP1:
      {
	if ( ns < 2 )
	  {
	    ok=false;
	  };
	for ( int t=0; t<ns-1; ++t )
	  {
	    for ( int e=std::max(t-1,0); e<std::min(t+2,ns-1); ++e )
	      {
		if ( avail[e+t*ns] != 1 )
		  {
		    ok=false;
		  }
	      }
	  }
	if ( avail[(ns-1)+(ns-2)*ns] != 1 )
	  {
	    ok = false;
	  }
	break;
      }

    case edgembar::BAREXP:
      {
	if ( ns < 3 )
	  {
	    ok=false;
	  };
	for ( int t=1; t<ns-1; ++t )
	  {
	    for ( int e=std::max(t-1,1); e<std::min(t+2,ns-1); ++e )
	      {
		if ( avail[e+t*ns] != 1 )
		  {
		    ok=false;
		  }
	      }
	  }
	if ( avail[0+1*ns] != 1 )
	  {
	    ok = false;
	  }
	if ( avail[(ns-1)+(ns-2)*ns] != 1 )
	  {
	    ok = false;
	  }
	break;
      }

    default:
      {
	std::cerr << "Programming error: Unimplemented mode in "
		  << "edgembar::CanEvalMode " << mode << "\n";
	std::exit(EXIT_FAILURE);
	break;
      }

    };

  return ok;
}


edgembar::ReadMode edgembar::AutoDetectMode
( int const ns,
  int const * avail )
{
  edgembar::ReadMode mode = edgembar::INVALID;
  if ( edgembar::CanEvalMode(ns,avail,edgembar::MBAR) )
    {
      mode = edgembar::MBAR;
    }
  else if ( edgembar::CanEvalMode(ns,avail,edgembar::MBAREXP0) )
    {
      mode = edgembar::MBAREXP0;
    }
  else if ( edgembar::CanEvalMode(ns,avail,edgembar::MBAREXP1) )
    {
      mode = edgembar::MBAREXP1;
    }
  else if ( edgembar::CanEvalMode(ns,avail,edgembar::MBAREXP) )
    {
      mode = edgembar::MBAREXP;
    }
  else if ( edgembar::CanEvalMode(ns,avail,edgembar::BAR) )
    {
      mode = edgembar::BAR;
    }
  else if ( edgembar::CanEvalMode(ns,avail,edgembar::BAREXP0) )
    {
      mode = edgembar::BAREXP0;
    }
  else if ( edgembar::CanEvalMode(ns,avail,edgembar::BAREXP1) )
    {
      mode = edgembar::BAREXP1;
    }
  else if ( edgembar::CanEvalMode(ns,avail,edgembar::BAREXP) )
    {
      mode = edgembar::BAREXP;
    }

  return mode;
}



std::vector< std::vector<int> > edgembar::GetSimIndexes
( int const ns,
  edgembar::ReadMode const mode )
{
  std::vector< std::vector<int> > idxs(ns);

  switch ( mode )
    {
      
    case edgembar::MBAR :
      {
	for ( int t=0; t<ns; ++t )
	  {
	    for ( int e=0; e<ns; ++e )
	      {
		idxs[t].push_back(e);
	      };
	  };
	break;
      }
      
    case edgembar::MBAREXP0 :
      {
	for ( int t=1; t<ns; ++t )
	  {
	    for ( int e=1; e<ns; ++e )
	      {
		idxs[t].push_back(e);
	      };
	  };
	idxs[1].push_back(0);
	break;
      }

    case edgembar::MBAREXP1 :
      {
	for ( int t=0; t<ns-1; ++t )
	  {
	    for ( int e=0; e<ns-1; ++e )
	      {
		idxs[t].push_back(e);
	      };
	  };
	idxs[ns-2].push_back(ns-1);
	break;
      }

    case edgembar::MBAREXP :
      {
	for ( int t=1; t<ns-1; ++t )
	  {
	    for ( int e=1; e<ns-1; ++e )
	      {
		idxs[t].push_back(e);
	      };
	  };
	idxs[1].push_back(0);
	idxs[ns-2].push_back(ns-1);
	break;
      }

    case edgembar::BAR :
      {
	for ( int t=0; t<ns; ++t )
	  {
	    for ( int e=std::max(0,t-1); e<std::min(ns,t+2); ++e )
	      {
		idxs[t].push_back(e);
	      };
	  };
	break;
      }

    case edgembar::BAREXP0 :
      {
	for ( int t=1; t<ns; ++t )
	  {
	    for ( int e=std::max(1,t-1); e<std::min(ns,t+2); ++e )
	      {
		idxs[t].push_back(e);
	      };
	  };
	idxs[1].push_back(0);
	break;
      }

      
    case edgembar::BAREXP1 :
      {
	for ( int t=0; t<ns-1; ++t )
	  {
	    for ( int e=std::max(0,t-1); e<std::min(ns-1,t+2); ++e )
	      {
		idxs[t].push_back(e);
	      };
	  };
	idxs[ns-2].push_back(ns-1);
	break;
      }

    case edgembar::BAREXP :
      {
	for ( int t=1; t<ns-1; ++t )
	  {
	    for ( int e=std::max(1,t-1); e<std::min(ns-1,t+2); ++e )
	      {
		idxs[t].push_back(e);
	      };
	  };
	idxs[1].push_back(0);
	idxs[ns-2].push_back(ns-1);
	break;
      }

    default:
      {
	std::cerr << "Programming error: Unimplemented mode in "
		  << "edgembar::GetSimIndexes " << mode << "\n";
	std::exit(EXIT_FAILURE);
	break;
      }
      
    };

  for ( int i=0; i<ns; ++i )
    {
      std::sort(idxs[i].begin(),idxs[i].end());
    }

  return idxs;
}
