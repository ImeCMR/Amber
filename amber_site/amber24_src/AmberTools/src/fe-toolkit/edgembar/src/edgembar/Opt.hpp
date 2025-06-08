#ifndef _edgembar_Opt_hpp_
#define _edgembar_Opt_hpp_

#include <vector>

namespace edgembar
{
  class Edge;
  
  double Optimize( edgembar::Edge const * e,
		   std::vector<double> const & pguess,
		   std::vector<double> & popt,
		   double const tol,
		   int const verbosity );
}

#endif
