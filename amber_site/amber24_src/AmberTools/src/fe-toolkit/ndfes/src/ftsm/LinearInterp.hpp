#ifndef _LinearInterp_hpp_
#define _LinearInterp_hpp_

#include <vector>

namespace ndfes
{
  class LinearInterp
  {
  public:

    LinearInterp() {}

    LinearInterp( int const n, double const * x, double const * y );

    void reset( int const n, double const * x, double const * y );

    bool push_back( double const x, double const y );
    
    double GetValue( double const x ) const;
    
    double GetDeriv( double const x ) const;

  private:

    int FindIdx( double const x ) const;

  public:
    std::vector<double> mX;
    std::vector<double> mY;
    std::vector<double> mB;

    
  };
  
}

#endif
