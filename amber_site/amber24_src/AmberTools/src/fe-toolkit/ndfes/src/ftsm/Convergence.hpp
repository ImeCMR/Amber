#ifndef _Convergence_hpp_
#define _Convergence_hpp_

#include <vector>
#include <memory>
#include <ostream>

#include "PathData.hpp"

namespace ndfes
{


  double GetSlope( int mlen, int nin, double const * vs );
  double GetMinSlope( int mlen, int nin, double const * vs );


  
  void FTSMCheckSameMeans
  ( ndfes::PathIter const & opath,
    ndfes::PathIter const & cpath,
    std::vector<double> const & deffc,
    double const beta,
    double const ptol,
    bool & samemeans,
    double & pmin,
    bool const verbose );

  
  void FTSMSlopeTest
  ( ndfes::PathOpt const & simpaths,
    //ndfes::PathOptions const & popts,
    int const mlen,
    double const distol,
    double const angtol,
    double const rmsdtol,
    std::vector<bool> const tidxisangle,
    double & distolgood,
    double & angtolgood,
    double & rmsdgood,
    bool & maxdispisbig,
    bool & avgdispisbig,
    std::vector<std::ostream *> couts );


  
  void OptDisplacements
  ( std::shared_ptr<ndfes::PCurve> pspl,
    std::shared_ptr<ndfes::PCurve> ospl,
    std::size_t ndim,
    std::size_t nsim,
    std::vector<double> & maxdisps,
    std::vector<double> & avgdisps );
    
}

#endif

