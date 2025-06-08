#ifndef _TubeUtils_hpp_
#define _TubeUtils_hpp_

#include <vector>
#include <memory>
#include <ostream>
#include "DimInfo.hpp"
#include "PCurve.hpp"

namespace ndfes
{

  class NearbyBin
  {
  public:
    
    NearbyBin()
      : gidx(0), size(0), simidx(0), reserved(false)
    {};
    
    NearbyBin( std::size_t const a )
      : gidx(a), size(0), simidx(0), reserved(false)
    {};

    NearbyBin( std::size_t const gidx,
	       std::vector<std::size_t> const & bidx,
	       std::vector<double> const & center )
      : gidx(gidx),
	bidx(bidx),
	center(center),
	size(0),
	simidx(0),
	reserved(false)
    {}

    bool operator==( ndfes::NearbyBin const & rhs ) const
    {
      return gidx == rhs.gidx;
    }

    double CalcMinDist( ndfes::DimInfo const & diminfo,
			std::vector<double> const & pts ) const;
    
    std::size_t gidx;
    std::vector<std::size_t> bidx;
    std::vector<double> center;
    std::size_t size;
    std::size_t simidx;
    bool reserved;
  };



  bool sortbysize( ndfes::NearbyBin const & a, ndfes::NearbyBin const & b );
  //bool operator==( ndfes::NearbyBin const & a, ndfes::NearbyBin const & b );


  
  class AreaSummary
  {
  public:
    AreaSummary() {}

    void push_back( std::vector<double> c, std::size_t s );
    
    void push_back( ndfes::NearbyBin const & bin );

    std::size_t size() const;

    std::size_t NumInclusive( std::size_t nlo, std::size_t nhi ) const;    

    std::size_t min() const;
    
    std::vector< std::vector<double> > mCenter;
    std::vector< std::size_t > mSize;
  };


  
  ndfes::DimInfo MakeBufferedGrid
  ( ndfes::DimInfo const & grid,
    std::size_t const conv_layers );
 
  ndfes::DimInfo MakeBufferedGrid
  ( ndfes::DimInfo const & grid,
    std::size_t const conv_layers,
    std::vector<double> rcs );
  
  ndfes::DimInfo MakeBufferedGrid
  ( ndfes::DimInfo const & grid,
    std::size_t const conv_layers,
    std::shared_ptr<ndfes::PCurve> pspl,
    std::size_t const npts );
  
  ndfes::DimInfo MakeBufferedGrid
  ( ndfes::DimInfo const & grid,
    std::size_t const conv_layers,
    std::vector<double> rcs,
    std::shared_ptr<ndfes::PCurve> pspl,
    std::size_t const npts );

  std::vector<std::size_t> GetMeshIdxs
  ( ndfes::DimInfo const & xDimInfo,
    std::size_t const conv_layers,
    std::size_t const * bidxs,
    bool const strict = true );

  
}


std::ostream & operator<<( std::ostream & cout, ndfes::AreaSummary const & s );
std::ostream & operator<<( std::ostream & cout, std::vector<ndfes::AreaSummary> const & s );


#endif
