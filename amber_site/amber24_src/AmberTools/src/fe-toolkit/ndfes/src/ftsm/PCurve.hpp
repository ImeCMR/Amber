#ifndef _ndfes_PCurve_hpp_
#define _ndfes_PCurve_hpp_

#include <vector>
#include "AkimaSpline.hpp"
#include "LinearInterp.hpp"

namespace ndfes
{

  class PCurve
  {
  public:

    PCurve( int myndim, int mnpts, double const * pts, double const * ts, bool const akima );
    PCurve( int myndim, int mnpts, double const * pts, bool const akima, int const maxit, int const nseg );

    void GetValue( double const t, double * vs ) const;
    
    void GetDeriv( double const t, double * vs ) const;

    std::vector<double> GetPointClosestTo( double const * pt ) const;

    int mNumDim;
    int mNumPts;
    std::vector<double> mT;
    std::vector<double> mX;
    bool mAkima;
    std::vector<ccdl::AkimaSpline> mASpl;
    std::vector<ndfes::LinearInterp> mLSpl;

  };
  
}

#endif

