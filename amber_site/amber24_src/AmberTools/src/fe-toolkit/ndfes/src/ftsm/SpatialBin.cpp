#include "SpatialBin.hpp"


ndfes::SpatialBin::SpatialBin()
  : glbidx(-1)
{
}


ndfes::SpatialBin::SpatialBin( ndfes::DimInfo const & info )
  : glbidx(-1),
    bidxs( info.GetNumDims(), 0 ),
    cidxs( info.GetNumCorners(), 0 ),
    center( info.GetNumDims(), 0 )
{
}

ndfes::SpatialBin::SpatialBin( ndfes::DimInfo const & info, std::size_t const gidx, std::size_t const * binidxs )
  : glbidx(gidx),
    bidxs( binidxs, binidxs + info.GetNumDims() ),
    cidxs( info.CptCornerGlbIdxs( binidxs ) ),
    center( info.CptBinCenter( binidxs ) )
{
}

