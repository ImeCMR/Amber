#include "PredictCentroids.hpp"
#include "BiasedMin.hpp"


std::vector<double> ndfes::PredictCentroids
( std::size_t const ndim,
  std::size_t const npts,
  std::vector<double> const & rcs,
  std::vector<bool> const & done,
  std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts )
{
  std::vector<double> centroids(ndim*npts,0);

  if ( popts.acc )
    {
      std::size_t my_npts = 0;
      std::vector<double> my_rcs;
      for ( std::size_t i=0; i<npts; ++i )
	{
	  if ( not done[i] )
	    {
	      ++my_npts;
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  my_rcs.push_back( rcs[dim+i*ndim] );
		}
	    }
	}

      std::vector<double> my_centroids(my_npts*ndim,0.);
      std::vector<double> fcs( my_rcs.size(), 0 );
      for ( std::size_t i=0; i<my_npts; ++i )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      fcs[dim+i*ndim] = popts.deffc[dim];
	    }
	};
      
      std::vector<double> my_wrcs(my_rcs);
      ndfes::WrapPath(fes->mDimInfo, my_wrcs);

      ndfes::BiasedMins( ndim, my_npts, my_wrcs.data(), fcs.data(),
			 my_centroids.data(), popts.acc_oobk,
			 fes, popts.acc_maxit, 1.e-13,
			 popts.minbounds, popts.maxbounds );

      ndfes::UnwrapCentroids( fes->mDimInfo, my_centroids, my_rcs );
      
      for ( std::size_t i=0, ii=0; i<npts; ++i )
	{
	  if ( not done[i] )
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  centroids[dim+i*ndim] = my_centroids[dim+ii*ndim];
		}
	      ++ii;
	    }
	  else
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  centroids[dim+i*ndim] = rcs[dim+i*ndim];
		}
	    }
	};
      
    }

  return centroids;
}



std::vector<double> ndfes::PredictCentroids
( std::size_t const ndim,
  std::size_t const npts,
  std::vector<double> const & rcs,
  std::vector<double> const & fcs,
  std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts )
{
  std::vector<double> centroids(ndim*npts,0);

  if ( popts.acc )
    {

      std::vector<double> wrcs( rcs );
      ndfes::WrapPath( fes->mDimInfo, wrcs );
      
      ndfes::BiasedMins( ndim, npts, wrcs.data(), fcs.data(),
			 centroids.data(), popts.acc_oobk,
			 fes, popts.acc_maxit, 1.e-13,
			 popts.minbounds, popts.maxbounds );

      ndfes::UnwrapCentroids( fes->mDimInfo, centroids, rcs );
    }

  return centroids;
}

std::vector<double> ndfes::PredictCentroids
( std::size_t const ndim,
  std::size_t const npts,
  std::vector<double> const & rcs,
  std::vector<double> const & fcs,
  std::vector<bool> const & done,
  std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts )
{
  std::vector<double> centroids(ndim*npts,0);

  if ( popts.acc )
    {
      std::size_t my_npts = 0;
      for ( std::size_t i=0; i<npts; ++i )
	{
	  if ( not done[i] )
	    {
	      ++my_npts;
	    }
	}
      std::vector<double> my_rcs(my_npts*ndim,0);
      std::vector<double> my_fcs(my_npts*ndim,0);
      std::vector<double> my_centroids(my_npts*ndim,0);
      for ( std::size_t i=0, ii=0; i<npts; ++i )
	{
	  if ( not done[i] )
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  my_rcs[dim+ii*ndim] = rcs[dim+i*ndim];
		  my_fcs[dim+ii*ndim] = fcs[dim+i*ndim];
		}
	      ++ii;
	    }
	}

      std::vector<double> my_wrcs(my_rcs);
      ndfes::WrapPath( fes->mDimInfo, my_wrcs );

      ndfes::BiasedMins( ndim, my_npts, my_wrcs.data(), my_fcs.data(),
			 my_centroids.data(), popts.acc_oobk,
			 fes, popts.acc_maxit, 1.e-13,
			 popts.minbounds, popts.maxbounds );

      ndfes::UnwrapCentroids( fes->mDimInfo, my_centroids, my_rcs );
      
      for ( std::size_t i=0, ii=0; i<npts; ++i )
	{
	  if ( not done[i] )
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  centroids[dim+i*ndim] = my_centroids[dim+ii*ndim];
		}
	      ++ii;
	    }
	  else
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  centroids[dim+i*ndim] = rcs[dim+i*ndim];
		};
	    }
	}
    }

  return centroids;
}
