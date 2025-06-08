#ifndef _ndfes_PredictCentroids_hpp_
#define _ndfes_PredictCentroids_hpp_

#include <vector>
#include <memory>
#include "FES.hpp"
#include "PathOptions.hpp"

namespace ndfes
{

  std::vector<double> PredictCentroids
  ( std::size_t const ndim,
    std::size_t const npts,
    std::vector<double> const & rcs,
    std::vector<bool> const & done,
    std::shared_ptr<ndfes::FES> fes,
    ndfes::PathOptions const & popts );
  
  std::vector<double> PredictCentroids
  ( std::size_t const ndim,
    std::size_t const npts,
    std::vector<double> const & rcs,
    std::vector<double> const & fcs,
    std::shared_ptr<ndfes::FES> fes,
    ndfes::PathOptions const & popts );
  
  std::vector<double> PredictCentroids
  ( std::size_t const ndim,
    std::size_t const npts,
    std::vector<double> const & rcs,
    std::vector<double> const & fcs,
    std::vector<bool> const & done,
    std::shared_ptr<ndfes::FES> fes,
    ndfes::PathOptions const & popts );
  
}

#endif
