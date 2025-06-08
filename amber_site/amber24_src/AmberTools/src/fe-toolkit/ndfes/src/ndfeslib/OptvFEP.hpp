#ifndef _ndfes_OptvFEP_hpp_
#define _ndfes_OptvFEP_hpp_

#include "SystemInfo.hpp"

namespace ndfes
{
  
  void OptvFEP
  ( std::size_t const nparam,
    ndfes::SystemInfo const & sysinfo,
    std::vector< std::vector<std::size_t> > const & sbinsamples,
    double const TOL,
    int const maxit,
    bool const verbose,
    std::vector<double> & state_fs );

  void OptvFEP_CheckGradients
  ( std::size_t const nparam,
    ndfes::SystemInfo const & sysinfo,
    std::vector< std::vector<std::size_t> > const & sbinsamples,
    std::vector<double> & state_fs );
  
}

#endif
