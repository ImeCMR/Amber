#ifndef _read_metafile_hpp_
#define _read_metafile_hpp_

#include <string>
#include <vector>
#include "SystemInfo.hpp"

namespace ndfes
{
  
  ndfes::pSystemInfo ReadMetafile
  ( std::string fname,
    std::size_t const NumHam,
    std::size_t const BsplOrder,
    std::vector<double> const & TargetWidths,
    std::vector<int> const & PeriodicDims,
    double const TargetT,
    double const StartFrac,
    double const StopFrac,
    double const deltaDoS,
    bool const sdosSingleRef,
    std::string const sdosHistPrefix,
    bool keeponlyuniq,
    int const maxuniq );
  
}

#endif
