#ifndef _ndfes_Tube_hpp_
#define _ndfes_Tube_hpp_

#include <vector>
#include <memory>
#include <ostream>

#include "DimInfo.hpp"
#include "PCurve.hpp"
#include "FES.hpp"
#include "PathOptions.hpp"
#include "TubeUtils.hpp"

namespace ndfes
{

  void SamplingConvPrint
  ( std::vector<ndfes::AreaSummary> const & simsizes,
    std::size_t const conv_layers,
    std::size_t const conv_samples,
    std::ostream & fh );

  bool SamplingConvResult
  ( std::vector<ndfes::AreaSummary> const & simsizes,
    std::size_t const conv_samples );

  
  std::vector<ndfes::AreaSummary> SamplingConvTest
  ( std::size_t const ndim,
    std::size_t const nsim,
    std::shared_ptr<ndfes::PCurve> pspl,
    std::shared_ptr<ndfes::FES> fes,
    std::size_t const conv_layers );


  std::vector<ndfes::AreaSummary> SamplingConvOpt
  ( std::size_t const ndim,
    std::size_t const nsim,
    double * rcs,
    double * fcs,
    std::shared_ptr<ndfes::PCurve> pspl,
    std::shared_ptr<ndfes::FES> fes,
    ndfes::PathOptions const & popts,
    std::size_t const conv_layers,
    std::vector<double> const & reservedpts );


  std::vector<ndfes::AreaSummary> SamplingConvOptThrowable
  ( std::size_t const ndim,
    std::size_t const nsim,
    double * rcs,
    double * fcs,
    std::shared_ptr<ndfes::PCurve> pspl,
    std::shared_ptr<ndfes::FES> fes,
    ndfes::PathOptions const & popts,
    std::size_t const conv_layers,
    std::vector<double> const & reservedpts );

  
}


#endif

