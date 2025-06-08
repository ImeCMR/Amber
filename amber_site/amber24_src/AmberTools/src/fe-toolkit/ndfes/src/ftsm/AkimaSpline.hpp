#ifndef _ccdl_AkimaSpline_hpp_
#define _ccdl_AkimaSpline_hpp_

#include <vector>

namespace ccdl
{
  class AkimaSpline
  {
  public:

    AkimaSpline() {};
    
    AkimaSpline( int const n, double const * xs, double const * ys );
    
    double GetValue( double const x ) const;
    
    double GetDeriv( double const x ) const;

    int GetNumPts() const { return mX.size(); }

    double const * GetXValues() const { return mX.data(); };
    
  private:
    
    int FindIdx( double const x ) const;
    
    std::vector<double> mX;
    std::vector<double> mY;
    std::vector<double> mB;
    std::vector<double> mC;
    std::vector<double> mD;
  };
  
}

#endif
