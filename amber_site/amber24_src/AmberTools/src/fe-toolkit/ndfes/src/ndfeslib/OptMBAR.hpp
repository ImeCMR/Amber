#ifndef _ndfes_OptMBAR_hpp_
#define _ndfes_OptMBAR_hpp_

#include "SystemInfo.hpp"

namespace ndfes
{
  
  void OptMBAR( ndfes::SystemInfo const & sysinfo,
		std::vector< std::vector<std::size_t> > const & state_samples,
		double const TOL,
		int const maxit,
		bool const verbose,
		std::vector<double> & state_fs );

  void OptMBAR_CheckGradients
  ( ndfes::SystemInfo const & sysinfo,
    std::vector< std::vector<std::size_t> > const & state_samples,
    std::vector<double> & state_fs );
  
}

#endif
