#ifndef _ndfes_SpatialBin_hpp_
#define _ndfes_SpatialBin_hpp_

#include <string>
#include <cstddef>

#include "DimInfo.hpp"

namespace ndfes
{

  class SpatialBin
  {
  public:
    
    SpatialBin();

    SpatialBin( ndfes::DimInfo const & dinfo );
    
    SpatialBin( ndfes::DimInfo const & dinfo,
		std::size_t const gidx,
		std::size_t const * binidxs );

    bool operator==( ndfes::SpatialBin const & rhs ) const;
    
    bool operator<( ndfes::SpatialBin const & rhs ) const;

    
  public:

    std::size_t glbidx;
    std::vector<std::size_t> bidxs;
    std::vector<std::size_t> cidxs;
    std::vector<double> center;
    
  };
  
}


#endif
