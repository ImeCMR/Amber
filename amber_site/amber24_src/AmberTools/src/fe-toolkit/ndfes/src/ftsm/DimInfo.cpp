#include <cmath>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "DimInfo.hpp"
#include "MeshUtils.hpp"
#include "PeriodicUtils.hpp"

int ndfes::ifloor( double const x )
{
  int c = 0;
  if ( x < 0 )
    {
      c = -( (int)(-x) + 1 );
    }
  else
    {
      c = (int)(x);
    }
  return c;
}

int ndfes::iceil( double const x )
{
  return ndfes::ifloor(x)+1;
}





void ndfes::WrapPath
( ndfes::DimInfo const & info, std::vector<double> & pts )
{
  std::size_t const ndim = info.GetNumDims();
  std::size_t const npts = pts.size() / ndim;
  
  int const * isper = info.GetIsPeriodic();
  double const * xmins = info.GetXmin();
  double const * xmaxs = info.GetXmax();
  
  for ( std::size_t ipt=0; ipt<npts; ++ipt )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  if ( isper[dim] )
	    {
	      pts[dim+ipt*ndim] = wraptorange( pts[dim+ipt*ndim], xmins[dim],
					       xmaxs[dim] );
	    }
	}
    }
}


void ndfes::WrapPathToNearestOccBin
( ndfes::DimInfo const & info,
  std::vector< std::vector<std::size_t> > const & mbinidxs,
  std::vector<double> & pts )
{
  std::size_t const ndim = info.GetNumDims();
  std::size_t const npts = pts.size() / ndim;
  
  int const * isper = info.GetIsPeriodic();
  //double const * xmins = info.GetXmin();
  //double const * xmaxs = info.GetXmax();

  std::size_t nbins = mbinidxs.size();
  std::vector< std::vector<double> > centers( nbins );
  for ( std::size_t i=0; i<nbins; ++i )
    {
      centers[i] = info.CptBinCenter( mbinidxs[i].data() );
    }

  //std::printf("nbins: %8lu\n",nbins);
  //ndfes::WrapPath( info, pts );


  std::vector<std::size_t> obins;
  std::vector<double> opts;
  
  for ( std::size_t ipt=0; ipt<npts; ++ipt )
    {
      double mindist = 1.e+30;
      std::size_t ibin = 0;
      for ( std::size_t i=0; i<nbins; ++i )
	{
	  double dist = 0;
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      double dx = pts[dim+ipt*ndim] - centers[i][dim];
	      if ( isper[dim] )
		{
		  dx = wrap(dx,360.);
		}
	      dist += dx*dx;
	    }
	  dist = std::sqrt(dist);
	  if ( dist < mindist )
	    {
	      mindist = dist;
	      ibin = i;
	      if ( mindist < 1.e-8 )
		{
		  break;
		}
	    }
	}

      std::size_t gidx = info.CptGlbBinIdx( mbinidxs[ibin].data() );
      
      std::vector<std::size_t>::iterator p
	= std::find( obins.begin(), obins.end(), gidx );
      
      if ( p == obins.end() )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      opts.push_back( centers[ibin][dim] );
	    }
	  obins.push_back( gidx );
	}
    }
  
  pts = opts;
}



void ndfes::UnwrapPath
( ndfes::DimInfo const & info, std::vector<double> & pts )
{
  std::size_t const ndim = info.GetNumDims();
  std::size_t const npts = pts.size() / ndim;
  
  int const * isper = info.GetIsPeriodic();
  double const * xmins = info.GetXmin();
  double const * xmaxs = info.GetXmax();
  
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      if ( isper[dim] )
	{
	  pts[dim+0*ndim] = wraptorange( pts[dim+0*ndim], xmins[dim],
					 xmaxs[dim] );
	}
    }
  
  for ( std::size_t ipt=1; ipt<npts; ++ipt )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  double dx = pts[dim+ipt*ndim] - pts[dim+(ipt-1)*ndim];
	  if ( isper[dim] )
	    {
	      dx = wrap(dx,360.);
	    }
	  pts[dim+ipt*ndim] = dx + pts[dim+(ipt-1)*ndim];
	}
    }
}



void ndfes::UnwrapCentroids
( ndfes::DimInfo const & info,
  std::vector<double> & pts,
  std::vector<double> const & rcs )
{
  std::size_t const ndim = info.GetNumDims();
  std::size_t const npts = pts.size() / ndim;
  
  int const * isper = info.GetIsPeriodic();
  //double const * xmins = info.GetXmin();
  //double const * xmaxs = info.GetXmax();
  
  for ( std::size_t ipt=0; ipt<npts; ++ipt )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  double dx = pts[dim+ipt*ndim] - rcs[dim+ipt*ndim];
	  if ( isper[dim] )
	    {
	      dx = wrap(dx,360.);
	    }
	  pts[dim+ipt*ndim] = dx + rcs[dim+ipt*ndim];
	}
    }
}







ndfes::DimInfo::DimInfo()
  : ndim(0),
    bsplorder(1),
    nbspl(2),
    ncorner(2)
{
}


ndfes::DimInfo::DimInfo
( int const indim,
  int const ibspl,
  int const * iisperiodic,
  double const * itarget_widths )
  : ndim(indim),
    bsplorder( std::max(ibspl,1) ),
    isperiodic(iisperiodic,iisperiodic+indim),
    target_widths(itarget_widths,itarget_widths+indim),
    xmin(indim,0.),
    xmax(indim,0.),
    dimsizes(indim,0)
{
  nbspl = bsplorder + (bsplorder%2);
  ncorner = std::pow(nbspl,ndim);
  ResetRange();
}

ndfes::DimInfo::DimInfo
( int const indim,
  int const ibspl,
  int const * iisperiodic,
  double const * imins,
  double const * imaxs,
  int const * isizes )
  : ndim(indim),
    bsplorder( std::max(ibspl,1) ),
    isperiodic(iisperiodic,iisperiodic+indim),
    target_widths(indim,0.),
    xmin(imins,imins+indim),
    xmax(imaxs,imaxs+indim),
    dimsizes(isizes,isizes+indim)
{
  nbspl = bsplorder + (bsplorder%2);
  ncorner = std::pow(nbspl,ndim);
  for ( int idim=0; idim<ndim; ++idim )
    {
      target_widths[idim] = (xmax[idim]-xmin[idim])/dimsizes[idim];
    }
}



void ndfes::DimInfo::ResetRange()
{
  for ( int idim=0; idim<ndim; ++idim )
    {
      xmin[idim] =  1.e+60;
      xmax[idim] = -1.e+60;
      dimsizes[idim] = 0;
    }
}

void ndfes::DimInfo::ModifyRange( double const * pt )
{
  for ( int idim=0; idim<ndim; ++idim )
    {
      double x = pt[idim];
      xmin[idim] = std::min(xmin[idim],x);
      xmax[idim] = std::max(xmax[idim],x);
    }
}

void ndfes::DimInfo::FinalizeRange()
{
  int nbuf = (nbspl-2)/2;

  for ( int dim=0; dim<ndim; ++dim )
    {

      
      if ( IsPeriodic(dim) )
	{
	  xmin[dim] = 0.;
	  xmax[dim] = 360.;
	  double width = xmax[dim]-xmin[dim];
	  dimsizes[dim] = (int)(0.5 + width/target_widths[dim]);
	  target_widths[dim] = width / dimsizes[dim];
	}
      else
	{
	  double w = target_widths[dim];
	  double xlo = xmin[dim] - nbuf*w;
	  double xhi = xmax[dim] + nbuf*w;


	  
	  int clo = ndfes::ifloor(xlo/w);
	  int chi = ndfes::iceil(xhi/w);
	  
	  //std::printf("xlo,xhi %2i %12.4f %12.4f -> %12.4f %12.4f\n",
	  //dim,xmin[dim],xmax[dim],clo * w,chi * w);

	  xmin[dim] = clo * w;
	  xmax[dim] = chi * w;
	  dimsizes[dim] = chi-clo;


	}

     
    }
  
  /*
  
  for ( int dim=0; dim<ndim; ++dim )
    {
      if ( IsPeriodic(dim) )
	{
	  xmin[dim] = 0.;
	  xmax[dim] = 360.;
	}
      else
	{
	  double w = target_widths[dim];
	  
	  xmin[dim] -= w*nbuf;
	  xmax[dim] += w*nbuf;

	  // Added this to make the ranges "round" wrt the width
	  xmin[dim] = std::floor(xmin[dim]/w)*w;
	  xmax[dim] = std::ceil(xmax[dim]/w)*w;
	}
    };

  for ( int dim=0; dim<ndim; ++dim )
    {
      double width = xmax[dim]-xmin[dim];
      if ( IsPeriodic(dim) )
	{
	  dimsizes[dim] = (int)(0.5 + width/target_widths[dim]);
	  target_widths[dim] = width / dimsizes[dim];
	}
      else
	{
	  dimsizes[dim] = 2 + width / target_widths[dim];
	  xmin[dim] -= target_widths[dim];
	  xmax[dim] = xmin[dim] + dimsizes[dim] * target_widths[dim];
	}
    };

  */
}


std::vector<std::size_t> ndfes::DimInfo::CptBinIdxs( double const * crd ) const
{
  std::vector<std::size_t> idxs( ndim, 0 );
  for ( int dim=0; dim<ndim; ++dim )
    {
      double x = crd[dim];
      if ( IsPeriodic(dim) )
	{
	  x = wraptorange(x,xmin[dim],xmax[dim]);
	}

      double frac = (x-xmin[dim])/target_widths[dim];
      idxs[dim] = (std::size_t)frac;
    }

  return idxs;
}


std::size_t ndfes::DimInfo::CptGlbBinIdx( std::size_t const * idxs ) const
{
  std::size_t gidx = 0;
  //std::cout << "dimsizes.size() " << dimsizes.size() << std::endl;
  for ( int dim=ndim; dim --> 0; )
    {
      //std::cout << "idxs[" << dim << "] = " << idxs[dim] << std::endl;
      //std::cout << "dimsizes[" << dim << "] = " << dimsizes[dim] << std::endl;

      gidx = idxs[dim] + gidx*((std::size_t)dimsizes[dim]);
    }
  return gidx;
}


std::size_t ndfes::DimInfo::CptGlbBinIdx( double const * crd ) const
{
  std::size_t gidx = 0;
  for ( int dim=ndim; dim --> 0 ; )
    {
      double x = crd[dim];
      if ( IsPeriodic(dim) )
	{
	  x = wraptorange(x,xmin[dim],xmax[dim]);
	}

      double frac = (x-xmin[dim])/target_widths[dim];
      std::size_t idx = (std::size_t)frac;
      gidx = idx + gidx*((std::size_t)dimsizes[dim]);
    }

  return gidx;
}


std::vector< std::vector<std::size_t> > ndfes::DimInfo::CptCornerIdxs( std::size_t const * binidxs ) const
{
  std::vector< std::vector<std::size_t> > linidx(ndim);
  for ( int dim=0; dim<ndim; ++dim )
    {
      for ( int ib=0; ib<nbspl; ++ib )
	{
	  int lidx = (int)binidxs[dim] + ib - (nbspl/2-1);
	  if ( IsPeriodic(dim) )
	    {
	      lidx = INTWRAP(lidx,dimsizes[dim]);
	    }
	  linidx[dim].push_back( (std::size_t)lidx );
	}
    }
  return linidx;
}


std::vector<std::size_t> ndfes::DimInfo::CptCornerGlbIdxs( std::vector< std::vector<std::size_t> > const & linidx ) const
{
  std::vector<std::size_t> glbidx(ncorner,0);
  std::vector<std::size_t> grididx(ndim*ncorner);
  ndfes::LinearSpacingsToMeshgrid( linidx, grididx );
  for ( int ic=0,nc=ncorner; ic<nc; ++ic )
    {
      std::size_t cidx=0;
      for ( int dim=ndim; dim-->0; )
	{
	  if ( IsPeriodic(dim) )
	    {
	      cidx = grididx[dim+ic*ndim] + (std::size_t)(cidx*dimsizes[dim]);
	    }
	  else
	    {
	      cidx = grididx[dim+ic*ndim] + (std::size_t)(cidx*(dimsizes[dim]+1));
	    }
	}
      glbidx[ic] = cidx;
    };
  return glbidx;
}


std::vector<std::size_t> ndfes::DimInfo::CptCornerGlbIdxs( std::size_t const * binidx ) const
{
  std::vector< std::vector<std::size_t> > linidxs( CptCornerIdxs(binidx) );
  return CptCornerGlbIdxs( linidxs );
}


std::vector<double> ndfes::DimInfo::CptBinCenter( std::size_t const * binidxs ) const
{
  std::vector<double> c( ndim, 0. );
  for ( int dim=0; dim<ndim; ++dim )
    {
      c[dim] = xmin[dim] + (binidxs[dim]+0.5)*target_widths[dim];
    }
  return c;
}


void ndfes::DimInfo::CptBinCenterFromPt( double const * pt, double * c ) const
{
  for ( int dim=0; dim<ndim; ++dim )
    {
      double const w = target_widths[dim];
      double x = pt[dim];
      if ( IsPeriodic(dim) )
	{
	  x = wraptorange(x,xmin[dim],xmax[dim]);
	}

      double sign = 1.;
      if ( x < 0 )
	{
	  sign = -1.;
	}
      int clo = ndfes::ifloor(x/w);
      //int chi = ndfes::iceil(x/w);
      c[dim] = (clo + sign*0.5) * w;
    }
}


void ndfes::DimInfo::WriteChkpt
( std::ostream & cout,
  int const indent ) const
{
#define FMTE std::scientific << std::setw(22) << std::setprecision(13)

  cout << "\n";
  
  for ( int dim=0; dim<ndim; ++dim )
    {
      cout << std::setw(indent) << ""
	   << "dim" << dim+1
	   << " = ndfes.SpatialDim( "
	   << FMTE << xmin[dim] << ", "
	   << FMTE << xmax[dim] << ", "
	   << std::setw(6) << dimsizes[dim] << ", ";
      if ( IsPeriodic(dim) )
	{
	  cout << std::setw(5) << "True";
	}
      else
	{
	  cout << std::setw(5) << "False";
	};
      cout << " )\n";
    }

  cout << "\n";
  
  cout << std::setw(indent) << ""
       << "grid = ndfes.VirtualGrid([";
  for ( int dim=0; dim<ndim; ++dim )
    {
      cout << "dim" << dim+1;
      if ( dim < ndim-1 )
	{
	  cout << ",";
	}
    }
  cout << "])\n";

#undef FMTE
}



