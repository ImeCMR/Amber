#ifndef _BiasedMin_hpp_
#define _BiasedMin_hpp_

#include <memory>
#include "FES.hpp"

namespace ndfes
{
  void BiasedMins
  ( std::size_t const ndim,
    std::size_t const nsim,
    double const * rcs,
    double const * fcs,
    double * obsmeans,
    double const oobk,
    std::shared_ptr<ndfes::FES> fes,
    std::size_t const maxit,
    double const TOL,
    std::vector<double> const & minbounds,
    std::vector<double> const & maxbounds );
}

#endif

