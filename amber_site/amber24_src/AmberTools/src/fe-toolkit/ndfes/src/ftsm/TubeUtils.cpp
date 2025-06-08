#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "TubeUtils.hpp"
#include "MeshUtils.hpp"
#include "PeriodicUtils.hpp"




bool ndfes::sortbysize( ndfes::NearbyBin const & a, ndfes::NearbyBin const & b )
{
  return a.size < b.size;
}

// bool ndfes::operator==( ndfes::NearbyBin const & a, ndfes::NearbyBin const & b )
// {
//   return a.gidx == b.gidx;
// }



ndfes::DimInfo ndfes::MakeBufferedGrid
( ndfes::DimInfo const & grid,
  std::size_t const conv_layers )
{
  std::size_t ndim = grid.GetNumDims();
    
  std::vector<double> minrange(grid.GetXmin(), grid.GetXmin() + ndim);
  std::vector<double> maxrange(grid.GetXmax(), grid.GetXmax() + ndim);
  std::vector<double> widths( grid.GetTargetWidths(), grid.GetTargetWidths() + ndim );
  std::vector<int> sizes( grid.GetDimSizes(), grid.GetDimSizes() + ndim );
  std::vector<int> isper( grid.GetIsPeriodic(), grid.GetIsPeriodic() + ndim );
  {
    for ( std::size_t dim=0; dim<ndim; ++dim )
      {
	if ( not isper[dim] )
	  {
	    double w = widths[dim];
	    double xlo = minrange[dim] - conv_layers*w;
	    // double xhi = maxrange[dim] + conv_layers*w;
	    // int clo = ndfes::ifloor(xlo/w);
	    // int chi = ndfes::iceil(xhi/w);
	    // minrange[dim] = clo * w;
	    // maxrange[dim] = chi * w;
	    // sizes[dim] = chi-clo;
	    sizes[dim] += 2*conv_layers;
	    minrange[dim] = xlo;
	    maxrange[dim] = xlo + w * sizes[dim];
	  }
      }
  }
  ndfes::DimInfo xDimInfo( ndim, 2, isper.data(), minrange.data(), maxrange.data(), sizes.data() );
  return xDimInfo;
}



ndfes::DimInfo ndfes::MakeBufferedGrid
( ndfes::DimInfo const & grid,
  std::size_t const conv_layers,
  std::vector<double> rcs )
{
  std::size_t ndim = grid.GetNumDims();
  std::size_t nsim = rcs.size() / ndim;
  
  ndfes::WrapPath( grid, rcs );
  
  std::vector<double> minrange(grid.GetXmin(), grid.GetXmin() + ndim);
  std::vector<double> maxrange(grid.GetXmax(), grid.GetXmax() + ndim);
  std::vector<double> widths( grid.GetTargetWidths(), grid.GetTargetWidths() + ndim );
  std::vector<int> sizes( grid.GetDimSizes(), grid.GetDimSizes() + ndim );
  std::vector<int> isper( grid.GetIsPeriodic(), grid.GetIsPeriodic() + ndim );
  {
    for ( std::size_t i=0; i<nsim; ++i )
      {
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    if ( not isper[dim] )
	      {
		minrange[dim] = std::min(minrange[dim],rcs[dim+i*ndim]);
		maxrange[dim] = std::max(maxrange[dim],rcs[dim+i*ndim]);
	      }
	  }
      };
    for ( std::size_t dim=0; dim<ndim; ++dim )
      {
	if ( not isper[dim] )
	  {
	    double w = widths[dim];
	    double xlo = minrange[dim] - conv_layers*w;
	    double xhi = maxrange[dim] + conv_layers*w;
	    int clo = ndfes::ifloor(xlo/w);
	    int chi = ndfes::iceil(xhi/w);
	    minrange[dim] = clo * w;
	    maxrange[dim] = chi * w;
	    sizes[dim] = chi-clo;
	  }
      }
  }
  ndfes::DimInfo xDimInfo( ndim, 2, isper.data(), minrange.data(), maxrange.data(), sizes.data() );
  return xDimInfo;
}


ndfes::DimInfo ndfes::MakeBufferedGrid
( ndfes::DimInfo const & grid,
  std::size_t const conv_layers,
  std::shared_ptr<ndfes::PCurve> pspl,
  std::size_t const npts )
{
  std::size_t ndim = grid.GetNumDims();
  
  //ndfes::WrapPath( grid, rcs );
  
  std::vector<double> minrange(grid.GetXmin(), grid.GetXmin() + ndim);
  std::vector<double> maxrange(grid.GetXmax(), grid.GetXmax() + ndim);
  std::vector<double> widths( grid.GetTargetWidths(), grid.GetTargetWidths() + ndim );
  std::vector<int> sizes( grid.GetDimSizes(), grid.GetDimSizes() + ndim );
  std::vector<int> isper( grid.GetIsPeriodic(), grid.GetIsPeriodic() + ndim );
  {
    double dt = 1./(npts-1.);
    for ( std::size_t i=0; i<npts; ++i )
      {
	std::vector<double> c(ndim,0);
	pspl->GetValue( i*dt, c.data() );
	ndfes::WrapPath(grid,c);
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    if ( not isper[dim] )
	      {
		minrange[dim] = std::min(minrange[dim],c[dim]);
		maxrange[dim] = std::max(maxrange[dim],c[dim]);
	      }
	  }
      };
    for ( std::size_t dim=0; dim<ndim; ++dim )
      {
	if ( not isper[dim] )
	  {
	    double w = widths[dim];
	    double xlo = minrange[dim] - conv_layers*w;
	    double xhi = maxrange[dim] + conv_layers*w;
	    int clo = ndfes::ifloor(xlo/w);
	    int chi = ndfes::iceil(xhi/w);
	    minrange[dim] = clo * w;
	    maxrange[dim] = chi * w;
	    sizes[dim] = chi-clo;
	  }
      }
  }
  ndfes::DimInfo xDimInfo( ndim, 2, isper.data(), minrange.data(), maxrange.data(), sizes.data() );
  return xDimInfo;
}





ndfes::DimInfo ndfes::MakeBufferedGrid
( ndfes::DimInfo const & grid,
  std::size_t const conv_layers,
  std::vector<double> rcs,
  std::shared_ptr<ndfes::PCurve> pspl,
  std::size_t const npts )
{
  std::size_t ndim = grid.GetNumDims();
  std::size_t nsim = rcs.size()/ndim;
  
  ndfes::WrapPath( grid, rcs );
  
  std::vector<double> minrange(grid.GetXmin(), grid.GetXmin() + ndim);
  std::vector<double> maxrange(grid.GetXmax(), grid.GetXmax() + ndim);
  std::vector<double> widths( grid.GetTargetWidths(), grid.GetTargetWidths() + ndim );
  std::vector<int> sizes( grid.GetDimSizes(), grid.GetDimSizes() + ndim );
  std::vector<int> isper( grid.GetIsPeriodic(), grid.GetIsPeriodic() + ndim );
  {
    double dt = 1./(npts-1.);
    for ( std::size_t i=0; i<npts; ++i )
      {
	std::vector<double> c(ndim,0);
	pspl->GetValue( i*dt, c.data() );
	ndfes::WrapPath(grid,c);
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    if ( not isper[dim] )
	      {
		minrange[dim] = std::min(minrange[dim],c[dim]);
		maxrange[dim] = std::max(maxrange[dim],c[dim]);
	      }
	  }
      };

    for ( std::size_t i=0; i<nsim; ++i )
      {
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    if ( not isper[dim] )
	      {
		minrange[dim] = std::min(minrange[dim],rcs[dim+i*ndim]);
		maxrange[dim] = std::max(maxrange[dim],rcs[dim+i*ndim]);
	      }
	  }
      }
    
    for ( std::size_t dim=0; dim<ndim; ++dim )
      {
	if ( not isper[dim] )
	  {
	    double w = widths[dim];
	    double xlo = minrange[dim] - conv_layers*w;
	    double xhi = maxrange[dim] + conv_layers*w;
	    int clo = ndfes::ifloor(xlo/w);
	    int chi = ndfes::iceil(xhi/w);
	    minrange[dim] = clo * w;
	    maxrange[dim] = chi * w;
	    sizes[dim] = chi-clo;
	  }
      }
  }
  ndfes::DimInfo xDimInfo( ndim, 2, isper.data(), minrange.data(), maxrange.data(), sizes.data() );
  return xDimInfo;
}




std::vector<std::size_t> ndfes::GetMeshIdxs
( ndfes::DimInfo const & xDimInfo,
  std::size_t const conv_layers,
  std::size_t const * bidxs,
  bool const strict )
{
  std::size_t const ndim = xDimInfo.GetNumDims();
  int const * sizes = xDimInfo.GetDimSizes();
  //double const * widths = xDimInfo.GetTargetWidths();
  int const * isper = xDimInfo.GetIsPeriodic();
  //double const * xmins = xDimInfo.GetXmin();
  //double const * xmaxs = xDimInfo.GetXmax();

  int didx = conv_layers;
  int nd = 2*didx+1;

  std::vector< std::vector<std::size_t> > lidxs(ndim);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      int dimsize = sizes[dim];
      int bidx = bidxs[dim];
      if ( isper[dim] )
	{
	  int idel = -didx;
	  for ( int i=0; i<nd; ++i, ++idel )
	    {
	      int k = bidx + idel;
	      lidxs[dim].push_back( INTWRAP( k, dimsize ) );
	    }
	}
      else
	{
	  int idel = -didx;
	  for ( int i=0; i<nd; ++i, ++idel )
	    {
	      int k = bidx + idel;
	      if ( (k < 0 or k >= dimsize) and strict )
		{
		  std::cerr << "Programming error getting extended region "
			    << " dim=" << dim << " i=" << i
			    << " k=" << k << " dimsize=" << dimsize << std::endl;
		  std::exit(EXIT_FAILURE);
		}
	      else
		{
		  lidxs[dim].push_back( k );
		}
	    }
	}
      //nmidxs *= lidxs[dim].size();
    }
  std::vector<std::size_t> midxs;
  ndfes::LinearSpacingsToMeshgrid( lidxs, midxs );
  return midxs;
}


double ndfes::NearbyBin::CalcMinDist
( ndfes::DimInfo const & diminfo,
  std::vector<double> const & pts ) const
{
  std::size_t const ndim = center.size();
  std::size_t const npts = pts.size() / ndim;
  double mindist = 1.e+30;
  for ( std::size_t i=0; i<npts; ++i )
    {
      double dist = 0;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  double dx = pts[dim+i*ndim] - center[dim];
	  if ( diminfo.IsPeriodic(dim) )
	    {
	      dx = wrap(dx,360.);
	    }
	  dist += dx*dx;
	}
      dist = std::sqrt(dist);
      mindist = std::min(mindist,dist);
    }
  return mindist;
}




void ndfes::AreaSummary::push_back( std::vector<double> c, std::size_t s )
{
  mCenter.push_back( c );
  mSize.push_back( s );
}

std::size_t ndfes::AreaSummary::min() const
{
  std::size_t m = 0;
  std::vector<std::size_t>::const_iterator p = std::min_element( mSize.begin(), mSize.end() );
  if ( p != mSize.end() )
    {
      m = *p;
    };
  return m;
}
    
void ndfes::AreaSummary::push_back( ndfes::NearbyBin const & bin )
{
  push_back( bin.center, bin.size );
}

std::size_t ndfes::AreaSummary::size() const { return mSize.size(); }

std::size_t ndfes::AreaSummary::NumInclusive( std::size_t nlo, std::size_t nhi ) const
{
  std::size_t cnt = 0;
  for ( std::size_t a=0; a<mSize.size(); ++a )
    {
      if ( mSize[a] >= nlo and mSize[a] <= nhi )
	{
	  ++cnt;
	}
    }
  return cnt;
}

std::ostream & operator<<( std::ostream & cout, ndfes::AreaSummary const & s )
{
  cout << std::setw(10) << s.size()
       << std::setw(10) << s.NumInclusive(0,0)
       << std::setw(10) << s.NumInclusive(1,50)
       << std::setw(10) << s.NumInclusive(51,100)
       << std::setw(10) << s.NumInclusive(101,1000000);
  return cout;
}

std::ostream & operator<<( std::ostream & cout, std::vector<ndfes::AreaSummary> const & s )
{
  std::size_t n = 0;
  std::size_t n0 = 0;
  std::size_t n50 = 0;
  std::size_t n100 = 0;
  std::size_t nmany = 0;

  cout << std::setw(5) << "Sim"
       << std::setw(10) << "Nbin"
       << std::setw(10) << "N=0"
       << std::setw(10) << "0<N<=50"
       << std::setw(10) << "50<N<=100"
       << std::setw(10) << "N>100"
       << "\n";

  for ( std::size_t i=0; i<s.size(); ++i )
    {
      n     += s[i].size();
      n0    += s[i].NumInclusive(0,0);
      n50   += s[i].NumInclusive(1,50);
      n100  += s[i].NumInclusive(51,100);
      nmany += s[i].NumInclusive(101,1000000);
      cout << std::setw(5) << i+1 << s[i] << "\n";
    }
  
  cout << std::setw(5) << "all"
       << std::setw(10) << n
       << std::setw(10) << n0
       << std::setw(10) << n50
       << std::setw(10) << n100
       << std::setw(10) << nmany
       << "\n";
  return cout;
}
