#include <iostream>
#include <iomanip>
#include <cstdio>
#include <algorithm>
#include <cmath>

#include "Tube.hpp"
#include "PeriodicUtils.hpp"
#include "MeshUtils.hpp"
#include "TubeUtils.hpp"
#include "PredictCentroids.hpp"







void ndfes::SamplingConvPrint
( std::vector<ndfes::AreaSummary> const & simsizes,
  std::size_t const conv_layers,
  std::size_t const conv_samples,
  std::ostream & fh )
{
  //bool isconv = true;
  //std::size_t const nsim = simsizes.size();

  fh << "\nFES sampling check around the new path.\n"
     << "What is the fewest number of bin samples within "
     << conv_layers << " layers of each control point?\n"
     << "If each nearby bin contains at least " << conv_samples
     << " samples, then we assume the FES is converged.\n\n";

  fh << simsizes;
  
  // std::size_t num_not_enough = 0;
  // for ( std::size_t isim=0; isim<nsim; ++isim )
  //   {
  //     fh << std::setw(4) << isim+1
  // 	      << std::setw(11) << simsizes[isim];
  //     if ( simsizes[isim] < conv_samples )
  // 	{
  // 	  isconv = false;
  // 	  fh << " F";
  // 	  num_not_enough += 1;
  // 	}
  //     else
  // 	{
  // 	  fh << " T";
  // 	};
  //     fh << std::endl;
  //   }
  
  // if ( isconv )
  //   {
  //     fh << "\nSampling appears to be converged" << std::endl;
  //   }
  // else
  //   {
  //     fh << "\nSampling is NOT converged. Areas near the path need more samples."
  // 	 << std::endl;
  //     fh << "There are " << num_not_enough
  // 	 << " simulations with insufficient sampling."
  // 	 << std::endl;
  //   }
}



bool ndfes::SamplingConvResult
( std::vector<ndfes::AreaSummary> const & simsizes,
  std::size_t const conv_samples )
{
  bool isconv = true;
  std::size_t const nsim = simsizes.size();
  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      if ( simsizes[isim].min() < conv_samples )
	{
	  isconv = false;
	}
    }
  return isconv;
}



std::vector<ndfes::AreaSummary> ndfes::SamplingConvTest
( std::size_t const ndim,
  std::size_t const nsim,
  std::shared_ptr<ndfes::PCurve> pspl,
  std::shared_ptr<ndfes::FES> fes,
  std::size_t const conv_layers )
{
  std::vector<ndfes::AreaSummary> simsizes(nsim);
  //typedef std::vector<std::size_t>::iterator vit;
  std::size_t const NINTERP = 2000;

  std::vector<double> rcs( ndim*nsim, 0 );
  {
    double dt = 1./(nsim-1.);
    for ( std::size_t i=0; i<nsim; ++i )
      {
	pspl->GetValue( i*dt, rcs.data() + i*ndim );
      }
  }
  ndfes::WrapPath( fes->mDimInfo, rcs );

  ndfes::DimInfo xDimInfo( ndfes::MakeBufferedGrid( fes->mDimInfo, conv_layers, rcs, pspl, NINTERP ) );


  
  
  

  // Figure out the occupied bin indexes on the extended grid
  
  std::size_t nbins = fes->mBins.size();
  //std::vector< std::vector<std::size_t> > xbidxs( nbins );
  //std::vector<std::size_t> xgidxs( nbins, 0 );
  typedef std::unordered_map<std::size_t,std::size_t> idxmap;
  typedef idxmap::iterator idxiter;
  idxmap xGlbIdxMap;
  xGlbIdxMap.reserve(nbins);
  for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      //xbidxs[ibin] = xDimInfo.CptBinIdxs( fes->mBins[ibin].center.data() );
      //xgidxs[ibin] = xDimInfo.CptGlbBinIdx( xbidxs[ibin].data() );
      std::vector<std::size_t> bidxs( xDimInfo.CptBinIdxs( fes->mBins[ibin].center.data() ) );
      std::size_t gidx = xDimInfo.CptGlbBinIdx( bidxs.data() );
      xGlbIdxMap.insert( { gidx, ibin } );
    }


  // Find ALL nearby bins
  //std::vector<ndfes::NearbyBin> nearbybins;
  typedef std::unordered_map<std::size_t,ndfes::NearbyBin> binmap;
  typedef binmap::iterator biniter;
  binmap nearbybins;

  std::vector<double> pathcrds(ndim*NINTERP,0);
  
  {
    std::size_t prev_gidx = 0;
    double const dt = 1./(NINTERP-1.);
    for ( std::size_t ipt=0; ipt<NINTERP; ++ipt )
      {
	std::vector<double> pt(ndim,0);
	pspl->GetValue( ipt*dt, pt.data() );
	ndfes::WrapPath( fes->mDimInfo, pt );
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    pathcrds[dim+ipt*ndim] = pt[dim];
	  }
	std::vector<std::size_t> bidxs( xDimInfo.CptBinIdxs( pt.data() ) );
	std::size_t gidx = xDimInfo.CptGlbBinIdx( bidxs.data() );
	if ( ipt > 0 and gidx == prev_gidx )
	  {
	    continue;
	  }
	else
	  {
	    prev_gidx = gidx;
	  };
	

	std::vector<std::size_t> midxs( ndfes::GetMeshIdxs( xDimInfo, conv_layers, bidxs.data() ) );
	std::size_t nmidxs = midxs.size() / ndim;
	
	for ( std::size_t i=0; i<nmidxs; ++i )
	  {
	    std::vector< std::size_t > mbidxs( midxs.data() + i*ndim,
					       midxs.data() + (i+1)*ndim );
	    
	    std::size_t mgidx = xDimInfo.CptGlbBinIdx( mbidxs.data() );
	    std::vector<double> mc = xDimInfo.CptBinCenter( mbidxs.data() );
	    ndfes::NearbyBin bin( mgidx, mbidxs, mc );
	    
	    //std::vector<ndfes::NearbyBin>::iterator p =
	    //std::find( nearbybins.begin(), nearbybins.end(), bin );

	    biniter pbin = nearbybins.find( mgidx );
	    
	    if ( pbin == nearbybins.end() )
	      {

		double r2min = 1.e+30;
		std::size_t imin = 0;
		for ( std::size_t isim=0; isim<nsim; ++isim )
		  {
		    double r2 = 0;
		    for ( std::size_t dim=0; dim<ndim; ++dim )
		      {
			double dx = rcs[dim+isim*ndim] - mc[dim];
			if ( xDimInfo.IsPeriodic(dim) )
			  {
			    dx = wrap(dx,360.);
			  };
			r2 += dx*dx;
		      }
		    if ( r2 < r2min )
		      {
			r2min = r2;
			imin = isim;
		      }
		  }
		bin.simidx = imin;

		std::size_t samples = 0;
		//vit p = std::find( xgidxs.begin(), xgidxs.end(), mgidx );
		idxiter p = xGlbIdxMap.find(mgidx);
		//if ( p != xgidxs.end() )
		if ( p != xGlbIdxMap.end() )
		  {
		    //std::size_t oidx = std::distance(xgidxs.begin(),p);
		    //samples = fes->mBinSizes[oidx];
		    samples = fes->mBinSizes[ p->second ];
		  }
		bin.size = samples;
		
		//nearbybins.push_back( bin );
		nearbybins.insert( {mgidx,bin} );
	      }
	  }
	
	
      }
  }



  //
  // Where are the target bins?
  //
  
  std::vector<double> simrcs( rcs.data(), rcs.data()+ndim*nsim );
  std::vector<std::size_t> simgidx(nsim,0);

  
  
  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      std::vector<std::size_t> bidxs( xDimInfo.CptBinIdxs( rcs.data() + isim*ndim ) );
      //std::size_t gidx = xDimInfo.CptGlbBinIdx( bidxs.data() );
      //std::vector<ndfes::NearbyBin>::iterator ptarget = nearbybins.end();
      biniter ptarget = nearbybins.end();

      std::size_t minsize = 1000000000;
      double mindist = 1.e+30;

      
      //for ( std::vector<ndfes::NearbyBin>::iterator
      for ( biniter
	      p=nearbybins.begin(); p != nearbybins.end(); ++p )
	{
	  if ( p->second.simidx == isim )
	    {
	      simsizes[isim].push_back( p->second );
	    
	      if ( not p->second.reserved )
		{
		  double dist = 0;
		  /////
		  dist = p->second.CalcMinDist(fes->mDimInfo,pathcrds);
		  /////
	      
		  if ( p->second.size < minsize or
		       (p->second.size == minsize and dist < mindist) )
		    {
		      ptarget = p;
		      minsize = p->second.size;
		      mindist = dist;
		    }
		}
	    }
	}
      

      
      if ( ptarget == nearbybins.end() )
	{
	  std::vector<std::size_t> midxs( ndfes::GetMeshIdxs( xDimInfo, conv_layers, bidxs.data() ) );
	  std::size_t nmidxs = midxs.size() / ndim;
	  
	  for ( std::size_t i=0; i<nmidxs; ++i )
	    {
	      std::vector< std::size_t > mbidxs( midxs.data() + i*ndim,
						 midxs.data() + (i+1)*ndim );
	      std::size_t mgidx = xDimInfo.CptGlbBinIdx( mbidxs.data() );
	      ndfes::NearbyBin bin(mgidx);
	      //std::vector<ndfes::NearbyBin>::iterator p = std::find(nearbybins.begin(), nearbybins.end(), bin );
	      biniter p = nearbybins.find(mgidx);
	      if ( p == nearbybins.end() )
		{
		  std::cerr << "Failed to locate gidx " << mgidx << " in nearbybins" << std::endl;
		  std::exit(EXIT_FAILURE);
		}
	      else if ( ! p->second.reserved )
		{
		  double dist = 0;
		  // for ( std::size_t dim=0; dim<ndim; ++dim )
		  //   {
		  //     double dx = rcs[dim+isim*ndim] - p->center[dim];
		  //     if ( fes->mDimInfo.IsPeriodic(dim) )
		  // 	{
		  // 	  dx = wrap(dx,360.);
		  // 	}
		  //     dist += dx*dx;
		  //   }
		  // dist = std::sqrt(dist);

		  ////////
		  dist = p->second.CalcMinDist(fes->mDimInfo,pathcrds);
		  ////////

		  
		  
		  
		  if ( p->second.size < minsize or
		       (p->second.size == minsize and dist < mindist) )
		    {
		      minsize = p->second.size;
		      mindist = dist;
		      ptarget = p;
		      //std::printf("isim %3lu size %5lu dist %12.4f gidx %8lu\n",
		      //isim,p->size,dist,p->gidx);
		    }
		}
	    }
	}
      
      if ( ptarget == nearbybins.end() )
	{
	  //std::cerr << "Failed to locate a target bin for isim " << isim << std::endl;
	  //std::exit(EXIT_FAILURE);
	  //simsizes[isim] = 99999;
	}
      else
        { 
          ptarget->second.reserved = true;
          simgidx[isim]  = ptarget->second.gidx;
          //simsizes[isim] = ptarget->size;
	};
      //double w = std::min(1.,(double)(ptarget->size) / DEF_FC_SIZE);
      //for ( std::size_t dim=0; dim<ndim; ++dim )
      //{
	  //simfcs[dim+isim*ndim] = (1-w) * minfcs[dim] + w*fcs[dim+isim*ndim];
      //  simrcs[dim+isim*ndim] = ptarget->center[dim];
      //}

    };

  return simsizes;
}











std::vector<ndfes::AreaSummary> ndfes::SamplingConvOpt
( std::size_t const ndim,
  std::size_t const nsim,
  double * rcs,
  double * fcs,
  std::shared_ptr<ndfes::PCurve> pspl,
  std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts,
  std::size_t const conv_layers,
  std::vector<double> const & reservedpts )
{
  std::vector<ndfes::AreaSummary> simsizes;
   try
   {
     simsizes = ndfes::SamplingConvOptThrowable
       (ndim,nsim,
	rcs,fcs,pspl,fes,popts,
	conv_layers,reservedpts);   
   }
   catch ( std::size_t simidx )
   {
      if ( conv_layers == 0 )
      {
         try
	 {
            simsizes = ndfes::SamplingConvOptThrowable
	      (ndim,nsim,
	       rcs,fcs,pspl,fes,popts,
	       conv_layers+1,reservedpts);
	 }
	 catch ( std::size_t newsimidx )
	 {
            std::cerr << "ndfes::SamplingConvOpt failed to locate a bin for simidx "
		      << simidx << " within " << conv_layers
		      << " layers of the path\n";
	    std::cerr << "Furthermore, it failed to locate a bin for simidx "
		      << newsimidx << " when the layers were increased to "
		      << conv_layers+1 << std::endl;
	    std::exit(EXIT_FAILURE);
	 }
      }
   }
   return simsizes;
}



std::vector<ndfes::AreaSummary> ndfes::SamplingConvOptThrowable
( std::size_t const ndim,
  std::size_t const nsim,
  double * rcs,
  double * fcs,
  std::shared_ptr<ndfes::PCurve> pspl,
  std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts,
  std::size_t const conv_layers,
  std::vector<double> const & reservedpts )
{
  std::vector<ndfes::AreaSummary> simsizes(nsim);
  //typedef std::vector<std::size_t>::iterator vit;
  std::size_t const NINTERP = 2000;
  //double const DEF_FC_SIZE = 25;

  //double const * minfcs = popts.minfc.data();
  double const * maxfcs = popts.maxfc.data();

  //std::size_t const nres = reservedpts.size();
  double minwidth = 1.e+30;
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      minwidth = std::min(minwidth,fes->mDimInfo.GetTargetWidths()[dim]);
    }
  
  // Create a new, extended DimInfo object
  
  std::vector<double> orig_xmins( fes->mDimInfo.GetXmin(), fes->mDimInfo.GetXmin() + ndim );
  std::vector<double> orig_xmaxs( fes->mDimInfo.GetXmax(), fes->mDimInfo.GetXmax() + ndim );

  std::vector<double> trcs(rcs,rcs+ndim*nsim);
  ndfes::WrapPath( fes->mDimInfo, trcs );
  std::copy(trcs.begin(),trcs.end(),rcs);
  
  ndfes::DimInfo xDimInfo( ndfes::MakeBufferedGrid( fes->mDimInfo, conv_layers, trcs, pspl, NINTERP ) );
  trcs.resize(0);
  

  // Figure out the occupied bin indexes on the extended grid
  
  std::size_t nbins = fes->mBins.size();
  //std::vector< std::vector<std::size_t> > xbidxs( nbins );
  //std::vector<std::size_t> xgidxs( nbins, 0 );
  //std::vector< std::vector<std::size_t> > occidxs;

  typedef std::unordered_map<std::size_t,std::size_t> idxmap;
  typedef idxmap::iterator idxiter;
  idxmap xGlbIdxMap;
  xGlbIdxMap.reserve(nbins);
  for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      //xbidxs[ibin] = xDimInfo.CptBinIdxs( fes->mBins[ibin].center.data() );
      //xgidxs[ibin] = xDimInfo.CptGlbBinIdx( xbidxs[ibin].data() );
      std::vector<std::size_t> bidxs( xDimInfo.CptBinIdxs( fes->mBins[ibin].center.data() ) );
      std::size_t gidx = xDimInfo.CptGlbBinIdx( bidxs.data() );
      //occidxs.push_back( fes->mBins[ibin].bidxs );
      xGlbIdxMap.insert( { gidx, ibin } );
    }


  // Find ALL nearby bins
  //std::vector<ndfes::NearbyBin> nearbybins;

  typedef std::unordered_map<std::size_t,ndfes::NearbyBin> binmap;
  typedef binmap::iterator biniter;
  binmap nearbybins;
  
  std::vector<double> pathcrds(ndim*NINTERP,0);

  
  {
    std::size_t prev_gidx = 0;
    double const dt = 1./(NINTERP-1.);
    for ( std::size_t ipt=0; ipt<NINTERP; ++ipt )
      {
	std::vector<double> pt(ndim,0);
	pspl->GetValue( ipt*dt, pt.data() );
	ndfes::WrapPath( fes->mDimInfo, pt );
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    pathcrds[dim+ipt*ndim] = pt[dim];
	  }
	std::vector<std::size_t> bidxs( xDimInfo.CptBinIdxs( pt.data() ) );
	std::size_t gidx = xDimInfo.CptGlbBinIdx( bidxs.data() );
	if ( ipt > 0 and gidx == prev_gidx )
	  {
	    continue;
	  }
	else
	  {
	    prev_gidx = gidx;
	  };
	

	std::vector<std::size_t> midxs( ndfes::GetMeshIdxs( xDimInfo, conv_layers, bidxs.data() ) );
	std::size_t nmidxs = midxs.size() / ndim;
	
	for ( std::size_t i=0; i<nmidxs; ++i )
	  {
	    std::vector< std::size_t > mbidxs( midxs.data() + i*ndim,
					       midxs.data() + (i+1)*ndim );
	    
	    std::size_t mgidx = xDimInfo.CptGlbBinIdx( mbidxs.data() );
	    std::vector<double> mc = xDimInfo.CptBinCenter( mbidxs.data() );
	    ndfes::NearbyBin bin( mgidx, mbidxs, mc );

	    
	    //std::vector<ndfes::NearbyBin>::iterator p =
	    //  std::find( nearbybins.begin(), nearbybins.end(), bin );

	    biniter pbin = nearbybins.find(mgidx);
	    
	    if ( pbin == nearbybins.end() )
	      {

		double r2min = 1.e+30;
		std::size_t imin = 0;
		for ( std::size_t isim=0; isim<nsim; ++isim )
		  {
		    double r2 = 0;
		    for ( std::size_t dim=0; dim<ndim; ++dim )
		      {
			double dx = rcs[dim+isim*ndim] - mc[dim];
			if ( xDimInfo.IsPeriodic(dim) )
			  {
			    dx = wrap(dx,360.);
			  };
			r2 += dx*dx;
		      }
		    if ( r2 < r2min )
		      {
			r2min = r2;
			imin = isim;
		      }
		  }
		bin.simidx = imin;

		std::size_t samples = 0;
		//vit p = std::find( xgidxs.begin(), xgidxs.end(), mgidx );
		idxiter p = xGlbIdxMap.find(mgidx);
		if ( p != xGlbIdxMap.end() )
		  {
		    //std::size_t oidx = std::distance(xgidxs.begin(),p);
		    //samples = fes->mBinSizes[oidx];
		    samples = fes->mBinSizes[ p->second ];
		  }
		bin.size = samples;

		// wait - why is this reserved?
		// oh, it's checking against a special set of external
		// reserved points that we are avoiding
		double minresdist = bin.CalcMinDist( xDimInfo, reservedpts );
		//std::cout << "reserved pts " << reservedpts.size() << " " << minresdist << std::endl;
		if ( minresdist < minwidth/2. )
		  {
		    bin.reserved = true;
		  }

		
		//nearbybins.push_back( bin );
		nearbybins.insert( {mgidx,bin} );
	      }
	  }
	
	
      }
  }

  //std::sort( nearbybins.begin(), nearbybins.end(), sortbysize );


  //
  // Where are the target bins?
  //
  
  std::vector<double> simrcs( rcs, rcs+ndim*nsim );
  std::vector<double> simfcs( fcs, fcs+ndim*nsim );
  std::vector<std::size_t> simgidx(nsim,0);

  
  
  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      std::vector<std::size_t> bidxs( xDimInfo.CptBinIdxs( rcs + isim*ndim ) );
      //std::vector<ndfes::NearbyBin>::iterator ptarget = nearbybins.end();
      biniter ptarget = nearbybins.end();

      std::size_t minsize = 1000000000;
      double mindist = 1.e+30;
      std::size_t mingidx = -1;

      
      //for ( std::vector<ndfes::NearbyBin>::iterator
      for ( biniter
	      p=nearbybins.begin(); p != nearbybins.end(); ++p )
	{

	  if ( p->second.simidx == isim )
	    {
	      simsizes[isim].push_back( p->second );
	    
	      if ( not p->second.reserved )
		{
		  double dist = 0;
		  
		  /////
		  dist = p->second.CalcMinDist(fes->mDimInfo,pathcrds);
		  /////
		  
		  if ( p->second.size < minsize or
		       (p->second.size == minsize and dist < mindist) or
		       (p->second.size == minsize and
			std::abs(dist - mindist) < 1.e-6 and
			p->second.gidx < mingidx ) )
		    {
		      ptarget = p;
		      minsize = p->second.size;
		      mindist = dist;
		      mingidx = p->second.gidx;
		    }
		}
	    }
	}
      

      
      if ( ptarget == nearbybins.end() )
	{
	  std::vector<std::size_t> midxs( ndfes::GetMeshIdxs( xDimInfo, conv_layers, bidxs.data() ) );
	  std::size_t nmidxs = midxs.size() / ndim;
	  
	  for ( std::size_t i=0; i<nmidxs; ++i )
	    {
	      std::vector< std::size_t > mbidxs( midxs.data() + i*ndim,
						 midxs.data() + (i+1)*ndim );
	      std::size_t mgidx = xDimInfo.CptGlbBinIdx( mbidxs.data() );
	      ndfes::NearbyBin bin(mgidx);
	      
	      //std::vector<ndfes::NearbyBin>::iterator p = std::find(nearbybins.begin(), nearbybins.end(), bin );
	      biniter p = nearbybins.find(mgidx);
	      
	      if ( p == nearbybins.end() )
		{
		  std::cerr << "Failed to locate gidx " << mgidx << " in nearbybins" << std::endl;
		  std::exit(EXIT_FAILURE);
		}
	      else if ( ! p->second.reserved )
		{
		  double dist = 0;
		  /////
		  dist = p->second.CalcMinDist(fes->mDimInfo,pathcrds);
		  /////

	      
		  if ( p->second.size < minsize or
		       (p->second.size == minsize and dist < mindist) or
		       (p->second.size == minsize and
			std::abs(dist - mindist) < 1.e-6 and
			p->second.gidx < mingidx) )
		    {
		      minsize = p->second.size;
		      mindist = dist;
		      ptarget = p;
		      mingidx = p->second.gidx;
		    }
		}
	    }
	}
      
      if ( ptarget == nearbybins.end() )
	{
	  //std::cerr << "Failed to locate a target bin for isim " << isim << std::endl;
	  //std::exit(EXIT_FAILURE);
	  throw isim;
	}
      
      ptarget->second.reserved = true;
      simgidx[isim] = ptarget->second.gidx;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  simfcs[dim+isim*ndim] = fcs[dim+isim*ndim];
	  simrcs[dim+isim*ndim] = ptarget->second.center[dim];
	}

    };

  

  {
    std::size_t mynres = 0;
    //for ( std::vector<ndfes::NearbyBin>::iterator
    for ( biniter
	    p=nearbybins.begin(); p!=nearbybins.end(); ++p )
      {
	if ( p->second.reserved )
	  {
	    ++mynres;
	  }
      }
    std::printf("Total reservations: %lu\n",mynres);
  }


  std::printf("\n\nSTAGE 1: Check if centroids of target bins remains\n");
  std::printf("         in the bin. If not, make a linear guess for\n");
  std::printf("         a proposed umbrella center.\n\n");

  std::vector<bool> isempty(nsim,false);
  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      if ( simsizes[isim].min() == 0 )
	{
	  isempty[isim] = true;
	}
      else
	{
	  isempty[isim] = false;
	}
    }



  
  std::vector<double> effmeans0(ndim*nsim,0);

  {
    effmeans0 = ndfes::PredictCentroids(ndim,nsim,simrcs,simfcs,isempty,fes,popts);
  }

  std::vector<std::size_t> propgidx( simgidx );
  std::vector<double> proprcs( simrcs );
  std::vector<bool> donefc(nsim,false);
  std::vector<bool> donerc(nsim,false);
  
  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      // std::vector<ndfes::NearbyBin>::iterator p =
      // 	std::find(nearbybins.begin(),nearbybins.end(),
      // 		  ndfes::NearbyBin(simgidx[isim]));

      biniter p = nearbybins.find( simgidx[isim] );
      
      if ( isempty[isim] )
	{
	  donefc[isim] = true;
	  donerc[isim] = true;
	  
	  std::printf("sim %3lu ",isim+1);
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      std::printf("%10.4f",p->second.center[dim]);
	    }
	  std::printf(" not moving: empty\n");
	}
      else
	{
	  std::vector<std::size_t> cen_bidxs( xDimInfo.CptBinIdxs( effmeans0.data() + isim*ndim ) );
	  std::size_t cen_gidx = xDimInfo.CptGlbBinIdx( cen_bidxs.data() );

	  std::vector<std::size_t> sim_bidxs( xDimInfo.CptBinIdxs( simrcs.data() + isim*ndim ) );
	  std::size_t sim_gidx = xDimInfo.CptGlbBinIdx( sim_bidxs.data() );


	  if ( sim_gidx != cen_gidx )
	    {
	      std::vector<double> est(ndim,0);
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  double disp = effmeans0[dim+isim*ndim] - simrcs[dim+isim*ndim];
		  est[dim] = simrcs[dim+isim*ndim] - disp;
		  est[dim] = std::min(orig_xmaxs[dim]-1.e-10,std::max(orig_xmins[dim]+1.e-10,est[dim]));
		}
	      cen_bidxs = xDimInfo.CptBinIdxs( est.data() );
	      cen_gidx = xDimInfo.CptGlbBinIdx( cen_bidxs.data() );
	    }


	  if ( sim_gidx == cen_gidx )
	    {
	      donefc[isim] = true;
	      donerc[isim] = true;
	      
	      std::printf("sim %3lu ",isim+1);
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  std::printf("%10.4f",p->second.center[dim]);
		}
	      std::printf(" not moving: centroid remains\n");
	      
	    }
	  else
	    {
	      	      
	      //std::vector<ndfes::NearbyBin>::iterator pcentroid =
	      //std::find(nearbybins.begin(),nearbybins.end(),ndfes::NearbyBin(cen_gidx));

	      biniter pcentroid = nearbybins.find( cen_gidx );
	      
	      if ( pcentroid == nearbybins.end() )
		{
		  std::size_t newsize = 0;
		  //vit pxg = std::find(xgidxs.begin(),xgidxs.end(),cen_gidx);
		  idxiter pxg = xGlbIdxMap.find(cen_gidx);
		  if ( pxg != xGlbIdxMap.end() )
		    {
		      //std::size_t oidx = std::distance(xgidxs.begin(),pxg);
		      //newsize = fes->mBinSizes[oidx];
		      newsize = fes->mBinSizes[ pxg->second ];
		    }
		  std::vector<double> c( xDimInfo.CptBinCenter( cen_bidxs.data() ) );

		  ndfes::NearbyBin bin(cen_gidx,cen_bidxs,c);
		  bin.size = newsize;
		  bin.simidx = isim;
		  bin.reserved = true;
		  //nearbybins.push_back( bin );
		  //pcentroid = nearbybins.begin() + nearbybins.size()-1;
		  pcentroid = nearbybins.insert( { cen_gidx, bin } ).first;

		  
		  donefc[isim] = false;
		  donerc[isim] = false;

		  if ( newsize == 0 )
		    {
		      donefc[isim] = true;
		    }
		      
		  propgidx[isim] = cen_gidx;
		  p->second.reserved = false;
		  
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      proprcs[dim+isim*ndim] = c[dim];
		    }
		  
		  std::printf("sim %3lu ",isim+1);
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      std::printf("%10.4f",p->second.center[dim]);
		    }
		  std::printf(" moving to: ");
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      std::printf("%10.4f",pcentroid->second.center[dim]);
		    }
		  std::printf(" (note: outside tube)\n");
		}
	      else
		{
		  if ( pcentroid->second.reserved )
		    {
		      donefc[isim] = false;
		      donerc[isim] = true;
		  
	
		      std::printf("sim %3lu ",isim+1);
		      for ( std::size_t dim=0; dim<ndim; ++dim )
			{
			  std::printf("%10.4f",p->second.center[dim]);
			}
		      std::printf(" not moving: reserved\n");
		      
		    }
		  else
		    {
		      donefc[isim] = false;
		      donerc[isim] = false;
		  
	
		      std::printf("sim %3lu ",isim+1);
		      for ( std::size_t dim=0; dim<ndim; ++dim )
			{
			  std::printf("%10.4f",p->second.center[dim]);
			}
		      std::printf(" moving to: ");
		      for ( std::size_t dim=0; dim<ndim; ++dim )
			{
			  std::printf("%10.4f",pcentroid->second.center[dim]);
			}
		      std::printf("\n");

		      propgidx[isim] = pcentroid->second.gidx;
		      pcentroid->second.reserved = true;
		      p->second.reserved = false;

		      for ( std::size_t dim=0; dim<ndim; ++dim )
			{
			  proprcs[dim+isim*ndim] = pcentroid->second.center[dim];
			}
		      
		    }
		}
	    };
	}
    }
	  


  {
    std::size_t mynres = 0;
    //for ( std::vector<ndfes::NearbyBin>::iterator
    for ( biniter
	    p=nearbybins.begin(); p!=nearbybins.end(); ++p )
      {
	if ( p->second.reserved )
	  {
	    ++mynres;
	  }
      }
    std::printf("Total reservations: %lu\n",mynres);
  }


  ////////////////////////////////////////////////////////////////////////////

  std::printf("\n\nSTAGE 2: If the simulation center has changed, make\n");
  std::printf("         sure its centroid is closer to the target\n");
  std::printf("         bin than the simply placing the center at\n");
  std::printf("         the bin. If not, revert the linear guess.\n\n");
  
  std::vector<double> effmeans1(ndim*nsim,0);

  {
    effmeans1 = ndfes::PredictCentroids(ndim,nsim,proprcs,simfcs,donerc,fes,popts);
  }
  

  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      if ( donerc[isim] )
	{
	  std::printf("sim %3lu ",isim+1);
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      std::printf("%10.4f",proprcs[dim+isim*ndim]);
	    }
	  std::printf(" untested\n");
	}
      else
	{

	  double r0 = 0;
	  double r1 = 0;
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      double d0 = effmeans0[dim+isim*ndim] - simrcs[dim+isim*ndim];
	      r0 += d0*d0;
	      double d1 = effmeans1[dim+isim*ndim] - simrcs[dim+isim*ndim];
	      r1 += d1*d1;
	    }
	  r0 = std::sqrt(r0);
	  r1 = std::sqrt(r1);

	  if ( r1 > r0 )
	    {
	      // std::vector<ndfes::NearbyBin>::iterator pold =
	      // 	std::find(nearbybins.begin(),nearbybins.end(),ndfes::NearbyBin(simgidx[isim]));
	      
	      // std::vector<ndfes::NearbyBin>::iterator pnew =
	      // 	std::find(nearbybins.begin(),nearbybins.end(),ndfes::NearbyBin(propgidx[isim]));

	      biniter pold = nearbybins.find( simgidx[isim] );
	      biniter pnew = nearbybins.find( propgidx[isim] );
	      
	      
	      if ( pold->second.reserved )
		{
		  std::printf("sim %3lu ",isim+1);
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      std::printf("%10.4f",proprcs[dim+isim*ndim]);
		    }
		  std::printf(" not reverting to ");
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      std::printf("%10.4f",simrcs[dim+isim*ndim]);
		    }
		  std::printf(" because it is reserved\n");
		}
	      else
		{

		  std::printf("sim %3lu ",isim+1);
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      std::printf("%10.4f",proprcs[dim+isim*ndim]);
		    }
		  std::printf(" reverting to ");
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      std::printf("%10.4f",simrcs[dim+isim*ndim]);
		    }
		  std::printf(" because %10.2e >= %10.2e\n",r1,r0);
		  
		  pold->second.reserved = true;
		  pnew->second.reserved = false;
		  propgidx[isim] = pold->second.gidx;
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      proprcs[dim+isim*ndim] = pold->second.center[dim];
		    }
		  
		}
	      
	    }
	  else
	    {
	      std::printf("sim %3lu ",isim+1);
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  std::printf("%10.4f",proprcs[dim+isim*ndim]);
		}
	      std::printf(" validated %10.2e <= %10.2e\n",r1,r0);
	    }
	  
	}
    }

  
     

  {
    std::size_t mynres = 0;
    //for ( std::vector<ndfes::NearbyBin>::iterator
    for ( biniter
	    p=nearbybins.begin(); p!=nearbybins.end(); ++p )
      {
	if ( p->second.reserved )
	  {
	    ++mynres;
	  }
      }
    std::printf("Total reservations: %lu\n",mynres);
  }
 

  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      //std::vector<ndfes::NearbyBin>::iterator p =
      //std::find(nearbybins.begin(),nearbybins.end(),ndfes::NearbyBin(propgidx[isim]));

      biniter p = nearbybins.find( propgidx[isim] );
      
      if ( p->second.size == 0 )
	{
	  donefc[isim] = true;
	}
    }



  {
  
    std::printf("\n\nSTAGE 3: Adjust force constants, if necessary\n\n");
    
    
    std::size_t ntry = 11;
    std::vector<double> ws(ntry,0);
    for ( std::size_t itry=0; itry < ntry; ++itry )
      {
	ws[itry] = (double)(itry) / ( (double)(ntry) - 1. );
      }
    
    for ( std::size_t itry=0; itry<ntry; ++itry )
      {
	std::size_t ncantry = 0;
	for ( std::size_t isim=0; isim<nsim; ++isim )
	  {
	    if ( not donefc[isim] )
	      {
		++ncantry;
	      }
	  }
	  
	std::printf("\nIteration %2lu ncantry %3lu scale %7.3f\n",itry+1,ncantry,ws[itry]);
	
	if ( ncantry == 0 )
	  {
	    break;
	  }
	
	std::vector<double> tryfcs(simfcs);
	for ( std::size_t isim=0; isim<nsim; ++isim )
	  {
	    if ( not donefc[isim] )
	      {
		for ( std::size_t dim=0; dim<ndim; ++dim )
		  {
		    tryfcs[dim+isim*ndim] = (1.-ws[itry]) * simfcs[dim+isim*ndim] + ws[itry] * maxfcs[dim];
		  }
	      }
	  }


	std::vector<double> ptry
	  ( ndfes::PredictCentroids(ndim,nsim,proprcs,tryfcs,donefc,fes,popts) );
	
	for ( std::size_t isim=0; isim<nsim; ++isim )
	  {
	    
	    if ( donefc[isim] )
	      {
		std::printf("isim %2lu fc not adjusted\n",isim+1);
	      }
	    else
	      {
		double maxf = -1.e+10;
		
		std::printf("isim %2lu fc ",isim+1);
		for ( std::size_t dim=0; dim<ndim; ++dim )
		  {
		    std::printf("%10.4f",tryfcs[dim+isim*ndim]);
		  }
		std::printf(" u: ");
		double const * widths = fes->mDimInfo.GetTargetWidths();
		for ( std::size_t dim=0; dim<ndim; ++dim )
		  {
		    double dx = std::abs( ptry[dim+isim*ndim] - simrcs[dim+isim*ndim] );
		    double u = dx / (0.5*widths[dim]) - 1.;
		    maxf = std::max(maxf,u);
		    std::printf("%8.4f",u);
		  }
		
		double r0 = 0;
		double r1 = 0;
		for ( std::size_t dim=0; dim<ndim; ++dim )
		  {
		    double dx = 0;
		    
		    dx = ptry[dim+isim*ndim] - simrcs[dim+isim*ndim];
		    r0 += dx*dx;
		    
		    dx = ptry[dim+isim*ndim] - proprcs[dim+isim*ndim];
		    r1 += dx*dx;
		  }
		r0 = std::sqrt(r0);
		r1 = std::sqrt(r1);
		std::printf(" r1: %10.2e r0: %10.2e",r1,r0);
		
		if ( maxf < 0 or itry == ntry - 1 )
		  {
		    std::printf(" FINAL\n");
		    donefc[isim] = true;
		    for ( std::size_t dim=0; dim<ndim; ++dim )
		      {
			simfcs[dim+isim*ndim] = tryfcs[dim+isim*ndim];
		      }
		  }
		else if ( r1 < r0 )
		  {
		    donefc[isim] = true;
		    std::printf(" reverting this iteration; moving away from target\n");
		    if ( itry > 0 )
		      {
			double w = ws[itry-1];
			for ( std::size_t dim=0; dim<ndim; ++dim )
			  {
			    simfcs[dim+isim*ndim] = (1.-w) * simfcs[dim+isim*ndim] + w * maxfcs[dim];
			  }
		      }
		  }
		else
		  {
		    std::printf("\n");
		  }
	      }
	  }
      }
  }
  

  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  rcs[dim+isim*ndim] = proprcs[dim+isim*ndim];
	  fcs[dim+isim*ndim] = simfcs[dim+isim*ndim];
	}
    };
  
  
  std::printf("\nFINAL RCS/FCS\n");
  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      std::printf("%2lu ",isim+1);
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  std::printf("%12.6f",rcs[dim+isim*ndim]);
	}
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  std::printf("%12.6f",fcs[dim+isim*ndim]);
	}
      std::printf("\n");
    }


  
  return simsizes;
}



