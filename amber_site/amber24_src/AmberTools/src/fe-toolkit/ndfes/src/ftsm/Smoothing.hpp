#ifndef _Smoothing_hpp_
#define _Smoothing_hpp_

#include <vector>

namespace ndfes
{
  std::vector<double> GetWinAvg
  ( int const n,
    double const * xs,
    int wlen );
  
  std::vector<double> GetIterWinAvg
  ( int const n,
    double const * xs,
    int wlen,
    int wmax );

  std::vector<double> GetIterWinAvg
  ( int const ndim,
    int const n,
    double const * xs,
    int const wlen,
    int const wmax );
  
}

#endif

