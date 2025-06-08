#define _USE_MATH_DEFINES
#include <cmath>

#include <algorithm>
#include <random>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "SystemInfo.hpp"
#include "GetBeta.hpp"
#include "PeriodicUtils.hpp"
#include "RunAvg.hpp"


namespace ndfes
{
  double logsumexp( std::size_t const n, double const * w, double const * z )
  {

    double maxz = -1.e+100;
    for ( std::size_t i=0; i<n; ++i )
      {
	if ( std::abs(w[i]) > 1.e-100 )
	  {
	    maxz = std::max(maxz,z[i]);
	  };
      };
    double s = 1.e-300;
    for ( std::size_t i=0; i<n; ++i )
      {
	if ( std::abs(w[i]) > 1.e-100 )
	  {
	    s += w[i]*std::exp(z[i]-maxz);
	  };
      }
    return std::log(s) + maxz;
  }

  double logsumexp( std::size_t const n, double const * z )
  {
    double maxz = -1.e+100;
    for ( std::size_t i=0; i<n; ++i )
      {
	maxz = std::max(maxz,z[i]);
      };
    double s = 1.e-300;
    for ( std::size_t i=0; i<n; ++i )
      {
	s += std::exp(z[i]-maxz);
      }
    return std::log(s) + maxz;
  }

  
}


namespace ndfes
{
  static std::random_device static_random_device;
  static std::mt19937 static_random_gen( ndfes::static_random_device() );

  void ResampleIdxs( std::size_t const g, std::vector<std::size_t> const & refa, std::vector<std::size_t> & a );

  std::size_t GetEneIdx( double const x, double const dx );

  std::vector<std::size_t> GetHashKeys( std::vector<std::size_t> const & seen );

  std::size_t GetHashIdx( std::size_t const key, std::vector<std::size_t> const & unique_keys );
  
  std::vector<std::size_t> CptHashTransform( std::vector<std::size_t> & seen );

  template<typename T>
  T HashLookup( std::size_t const key, std::vector<T> const & values, std::vector<std::size_t> const & unique_keys );
}


void ndfes::ResampleIdxs
( std::size_t const g,
  std::vector<std::size_t> const & refa,
  std::vector<std::size_t> & a )
{  
  a.resize(0);
 
  std::size_t n = refa.size() / 2;
  std::size_t nblks = n/g;
  if ( n%g > 0 )
    {
      nblks += 1;
    }
  
  if ( n > 0 )
    {
      std::uniform_int_distribution<std::size_t> dist(0,n-1);
      std::vector<std::size_t> vals;
      for ( std::size_t k=0; k<nblks; ++k )
	{
	  std::size_t refk = dist(ndfes::static_random_device);
	  for ( std::size_t ig=0; ig<g; ++ig, ++refk )
	    {
	      refk = refk % n;
	      if ( (std::size_t)vals.size() == n )
		{
		  break;
		}
	      else
		{
		  vals.push_back( refa[0+refk*2] );
		}
	    };
	}
      if ( (std::size_t)vals.size() != n )
	{
	  std::cerr << "Block resampling failed. Generated "
		    << vals.size() << " points, but expected "
		    << n << std::endl;
	  std::exit(1);
	}
      std::sort(vals.begin(),vals.end());

      a.push_back( vals[0] );
      a.push_back( 1 );
      for ( std::size_t k=1; k<n; ++k )
	{
	  if ( vals[k] != vals[k-1] )
	    {
	      a.push_back( vals[k] );
	      a.push_back( 1 );
	    }
	  else
	    {
	      a.back() += 1;
	    }
	}
    } 
}



//
// Return a signed index of the bin containing "x".
// The central bin is centered on 0 and has width "dx";
// that is, the central bin captures all x values in the
// range [-dx/2,+dx/2)
//
inline std::size_t ndfes::GetEneIdx
( double const x,
  double const dx )
{
  return x >= 0 ? x/dx+0.5 : -(std::size_t)(-x/dx+0.5);
}


//
// Return the list of sorted, unique elements in "seen"
//
inline std::vector<std::size_t> ndfes::GetHashKeys
( std::vector<std::size_t> const & seen )
{
  std::vector<std::size_t> unique_keys(seen);
  std::sort(unique_keys.begin(),unique_keys.end());
  std::vector<std::size_t>::iterator end=std::unique(unique_keys.begin(),unique_keys.end());
  unique_keys.resize(std::distance(unique_keys.begin(),end));
  return unique_keys;
}

//
// Return the index within "unique_keys" that matches "key"
//
inline std::size_t ndfes::GetHashIdx
( std::size_t const key,
  std::vector<std::size_t> const & unique_keys )
{
  return std::distance(unique_keys.begin(),std::find(unique_keys.begin(),unique_keys.end(),key));
}

//
// Return the value in "values" using the same index that matches the "key" within "unique_keys"
//
template<typename T>
inline T ndfes::HashLookup
( std::size_t const key,
  std::vector<T> const & values,
  std::vector<std::size_t> const & unique_keys )
{
  // assert values.size() == unique_keys.size()
  return values[GetHashIdx(key,unique_keys)];
}

//
// Returns a list of unique values within "seen"
// The values of the "seen" array are modified to point
// to the element within the list of unique values
//
inline std::vector<std::size_t> ndfes::CptHashTransform
( std::vector<std::size_t> & seen )
{
  std::vector<std::size_t> unique_keys(ndfes::GetHashKeys(seen));
  for ( std::size_t i=0, n=seen.size(); i<n; ++i )
    {
      seen[i] = ndfes::GetHashIdx(seen[i],unique_keys);
    }
  return unique_keys;
}


ndfes::SystemInfo::SystemInfo()
  : temperature(0.),
    beta(0.),
    deltaDoS(-1.),
    sdosFullSampling(false),
    sdosHistPrefix(""),
    nham(0)
{
}


ndfes::SystemInfo::SystemInfo
( ndfes::DimInfo const & dinfo,
  double const T,
  double const DeltaDoS,
  bool FullSampling,
  std::string histPrefix,
  std::size_t const nhams )
  : diminfo(dinfo),
    temperature(T),
    beta( ndfes::GetBeta(T) ),
    deltaDoS(DeltaDoS),
    sdosFullSampling(FullSampling),
    sdosHistPrefix(histPrefix),
    nham(nhams)
{
}


void ndfes::SystemInfo::SetTargetTemp( double const T )
{
  temperature = T;
  beta = ndfes::GetBeta(T);
}


void ndfes::SystemInfo::SetStatIneff
( std::size_t const istate,
  std::size_t const g )
{
  if ( istate < (std::size_t)states.size() )
    {
      states[istate].SetStatIneff(g);
    }
  else
    {
      std::cerr << "ndfes::SystemInfo::SetStatIneff invalid state index " << istate << " : " << states.size() << std::endl;
      std::exit(1);
    }
}


std::size_t ndfes::SystemInfo::InsertState
( std::size_t const iiham,
  double const * icenter,
  double const * ifconst,
  double const itemperature )
{
  ndfes::State s( diminfo.GetNumDims(),iiham,icenter,ifconst,itemperature);
  return InsertState( s );
}


std::size_t ndfes::SystemInfo::InsertState
( ndfes::State const & s )
{
  // Check to see if this is a new state or if we've already read this state
  std::size_t istate = 0;
  std::size_t ns = states.size();
  bool found = false;
  for ( std::size_t i=0; i<ns; ++i )
    {
      if ( s == states[i] )
	{
	  istate = i;
	  found = true;
	  break;
	}
    }

  if ( not found )
    {
      // This is a new state, so append it to the list of states
      istate = states.size();
      states.push_back( s );
    }

  // Return the state idx
  return istate;
}


void ndfes::SystemInfo::CptBiasEnergies( double * W, double const * pt ) const
{
  int const * dimisper = diminfo.GetIsPeriodic();
  std::size_t ns = states.size();
  for ( std::size_t i=0; i<ns; ++i )
    {
      W[i] = states[i].CptBiasEnergy( dimisper, pt );
    }
}


void ndfes::SystemInfo::CptSpatialHistogram()
{
  // Spatial bin objects
  sbins.resize(0);

  // array of sample idxs for each spatial bin
  sbinsamples.resize(0);

  // array of sample idxs for each state
  statesamples.resize(0);
  statesamples.resize( states.size() );
  
  // list of unique corner idxs
  ucidxs.resize(0);

  std::size_t npts = samples.size();
  
  // calculate min, max, and spacing of the grid
  diminfo.ResetRange();
  for ( std::size_t ipt=0; ipt<npts; ++ipt )
    {
      diminfo.ModifyRange( samples[ipt].GetPt() );
    }
  diminfo.FinalizeRange();
  
  
  typedef std::vector<ndfes::SpatialBin>::iterator iter;
  typedef std::vector<std::size_t>::iterator viter;
  for ( std::size_t ipt=0; ipt<npts; ++ipt )
    {

      {
	// Append the sample to the list of samples provided by state istate
	std::size_t istate = samples[ipt].GetStateIdx();
	statesamples[istate].push_back( ipt );
	statesamples[istate].push_back( 1 );
      }
      
      // The binidxs array has length ndim; the bin idx in each dimension
      std::vector<std::size_t> binidxs( diminfo.CptBinIdxs( samples[ipt].GetPt() ) );
      
      // The unique index of this bin
      std::size_t glbidx = diminfo.CptGlbBinIdx( binidxs.data() );
      
      // Determine if this is a new bin or a bin we already added
      iter p = FindBin( glbidx );
      
      if ( p == sbins.end() )
	{
	  // This is a new bin

	  // Modify the sample so it points at the new bin
	  samples[ipt].SetSpatialBinArrayIdx( sbins.size() );

	  // Create the new bin and push it to the list of bins
	  sbins.push_back( ndfes::SpatialBin( diminfo, glbidx, binidxs.data() ) );

	  // Append the sample to a the list of samples in this bin
	  // The list of samples has 2 values: the sample idx and the degeneracy
	  // When bootstrapping, we can include only the unique sample idxs and
	  // increase the degeneracy
	  std::vector<std::size_t> bsamples(2,0);
	  bsamples[0] = ipt;
	  bsamples[1] = 1;
	  
	  sbinsamples.push_back(bsamples);

	  // Append the bins corners to the list of corner idxs 
	  for ( viter q=sbins.back().cidxs.begin(),
		  qend=sbins.back().cidxs.end(); q!=qend; ++q )
	    {
	      ucidxs.push_back( *q );
	    }
	}
      else
	{

	  if ( binidxs[0] != p->bidxs[0] )
	    {
	      std::cout << "ERROR in ndfes::SystemInfo::CptSpatialHistogram() "
			<< binidxs[0] << " " << p->bidxs[0] << "\n";
	    }
	  
	  // This bin has already been created

	  // The bin idx in the array of bins
	  std::size_t ibin = std::distance( sbins.begin(), p );

	  // Modify the sample so it points at the bin
	  samples[ipt].SetSpatialBinArrayIdx( ibin );

	  // Append the sample to a the list of samples in this bin
	  sbinsamples[ibin].push_back( ipt );
	  sbinsamples[ibin].push_back( 1 );
	};
    }

  {
    // Make ucidxs a list of unique corner idxs
    std::sort( ucidxs.begin(), ucidxs.end() );
    viter p = std::unique( ucidxs.begin(), ucidxs.end() );
    ucidxs.resize( std::distance( ucidxs.begin(), p ) );
  }
}


std::vector<ndfes::SpatialBin>::iterator ndfes::SystemInfo::FindBin( std::size_t const gidx )
{
  typedef std::vector<ndfes::SpatialBin>::iterator iter;

  iter p = sbins.begin();
  for ( iter pend=sbins.end(); p!=pend; ++p )
    {
      if ( p->glbidx == gidx )
	break;
    }
  return p;
}


std::vector<ndfes::SpatialBin>::const_iterator ndfes::SystemInfo::FindBin( std::size_t const gidx ) const
{
  typedef std::vector<ndfes::SpatialBin>::const_iterator iter;

  iter p = sbins.begin();
  for ( iter pend=sbins.end(); p!=pend; ++p )
    {
      if ( p->glbidx == gidx )
	break;
    }
  return p;
}



void ndfes::SystemInfo::GetSampleIdxs
( std::vector< std::vector<std::size_t> > & ssamples,
  std::vector< std::vector<std::size_t> > & bsamples ) const
{
  ssamples = statesamples;
  bsamples = sbinsamples;
}


void ndfes::SystemInfo::ResampleStates
( std::vector< std::vector<std::size_t> > & ssamples,
  std::vector< std::vector<std::size_t> > & bsamples ) const
{
  std::size_t ns = statesamples.size();
  std::size_t nb = sbinsamples.size();
  ssamples.resize( ns );
  bsamples.resize(0);
  bsamples.resize( nb );
  for ( std::size_t i=0; i<ns; ++i )
    {
      std::size_t g = states[i].GetStatIneff();
      ndfes::ResampleIdxs( g, statesamples[i], ssamples[i] );
      std::size_t n = ssamples[i].size()/2;
      for ( std::size_t k=0; k<n; ++k )
	{
	  std::size_t isample = ssamples[i][0+k*2];
	  std::size_t ibin = samples[isample].GetSpatialBinArrayIdx();
	  bsamples[ibin].push_back(isample);
	  bsamples[ibin].push_back(ssamples[i][1+k*2]);
	}	  
    }
}


// void ndfes::SystemInfo::ResampleBins
// ( std::vector< std::vector<std::size_t> > & bsamples ) const
// {
//   std::size_t nb = sbinsamples.size();
//   bsamples.resize( nb );
//   for ( std::size_t i=0; i<nb; ++i )
//     {
//       ndfes::ResampleIdxs( sbinsamples[i], bsamples[i] );
//     }
// }


void ndfes::SystemInfo::ResampleBin
( std::size_t ibin,
  std::vector<std::size_t> & bsamples ) const
{
  if ( ibin < sbinsamples.size() )
    {
      ndfes::ResampleIdxs( 1, sbinsamples[ibin], bsamples );
    }
  else
    {
      std::cerr << "Error: ibin ("
		<< ibin
		<< ") out of range [0,"
		<< sbinsamples.size()
		<< "] in ndfes::SystemInfo::ResampleBin"
		<< std::endl;
      std::exit(1);
    }
}


void ndfes::SystemInfo::ShiftPotEnesToAvg()
{
  std::size_t nsample = samples.size();
  std::vector<double> avgs(nham,0.);

  for ( std::size_t ipt=0; ipt<nsample; ++ipt )
    {
      double const * pU = samples[ipt].GetPotEnes();
      for ( std::size_t iham=0; iham < nham; ++iham )
	{
	  avgs[iham] += pU[iham];
	};
    };
  for ( std::size_t iham=0; iham < nham; ++iham )
    {
      avgs[iham] = - avgs[iham] / nsample;
    }
  for ( std::size_t ipt=0; ipt<nsample; ++ipt )
    {
      samples[ipt].AddToPotEnes( avgs.data() );
    }
}





void ndfes::SystemInfo::PrintSummary( std::ostream & cout ) const
{
  std::size_t ns = GetNumStates();
  std::size_t ndim = diminfo.GetNumDims();
  std::size_t nbin = sbinsamples.size();

#define FMTW5 std::setw(5)
#define FMTW6 std::setw(6)
#define FMTW8 std::setw(8)
#define FMTW12 std::setw(12)
#define FMTF8 std::fixed << std::setw(8) << std::setprecision(2)
#define FMTF12 std::fixed << std::setw(12) << std::setprecision(6)
  
  cout << "\n";
  
  std::size_t totn = 0;

  for ( std::size_t i=0; i<ns; ++i )
    {
      std::size_t npts = statesamples[i].size()/2;
      totn += npts;
    };
  
  cout << "Num. dimensions:          " << ndim << "\n";
  cout << "Num. unique states:       " << ns << "\n";
  cout << "Num. occ. spatial bins:   " << nbin << "\n";
  cout << "Num. target Hamiltonians: " << nham << "\n";
  cout << "Num. Samples:             " << totn << "\n";

  cout << "Target temperature:   " << FMTF8 << temperature << "\n";

  
  cout << "\n";
  cout << FMTW5 << "State"
       << FMTW5 << "Ham"
       << FMTW8 << "T"
       << FMTW6 << "Ineff"
       << FMTW8 << "N";
    

  
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      std::stringstream ss;
      ss << "X" << dim+1;
      cout << FMTW8 << ss.str();
      
      ss.str("");
      ss.clear();
      ss << "K" << dim+1;
      cout << FMTW8 << ss.str();
      
      ss.str("");
      ss.clear();
      ss << "<X" << dim+1 << ">";
      cout << FMTW8 << ss.str();
      
      ss.str("");
      ss.clear();
      ss << "SD X" << dim+1;
      cout << FMTW8 << ss.str();
    }
  cout << "\n";
  cout << FMTW5 << "idx"
       << FMTW5 << "idx"
       << "\n";

  for ( std::size_t i=0; i<ns; ++i )
    {
      std::size_t npts = statesamples[i].size()/2;
      std::size_t g = states[i].GetStatIneff();
      
      std::vector<ndfes::RunAvg> CVs(ndim);
      for ( std::size_t ipt=0; ipt<npts; ++ipt )
	{
	  std::size_t idx = statesamples[i][0+ipt*2];
	  double const * xs = samples[idx].GetPt();
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      CVs[dim].push_back( xs[dim] );
	    };
	}
      
      cout << FMTW5 << i
	   << FMTW5 << states[i].GetHamIdx()
	   << FMTF8 << states[i].GetTemperature()
	   << FMTW6 << g
	   << FMTW8 << npts;
      
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  cout << FMTF8 << states[i].GetCenter()[dim]
	    //<< FMTF8 << states[i].ConvertKTToKcal( states[i].GetFConst()[dim] )
	       << FMTF8 << states[i].GetFConst()[dim]
	       << FMTF8 << CVs[dim].mean()
	       << FMTF8 << CVs[dim].stddev();
	}
      cout << "\n";
    }
  cout << "\n\n";

  cout << "Virtual Grid Information\n\n";
  cout << FMTW5 << "Dim"
       << FMTW12 << "Min Val"
       << FMTW12 << "Max Val"
       << FMTW12 << "Bin Width"
       << FMTW8 << "Size"
       << " " << FMTW8 << "Periodic"
       << "\n";
  for ( std::size_t dim=0; dim<ndim; ++dim)
    {
      double xmin = diminfo.GetXmin()[dim];
      double xmax = diminfo.GetXmax()[dim];
      double w = diminfo.GetTargetWidths()[dim];
      std::size_t nx = (std::size_t)( (xmax-xmin)/w + 0.5);
      cout << FMTW5 << dim+1
	   << FMTF12 << xmin
	   << FMTF12 << xmax
	   << FMTF12 << w
	   << FMTW8 << nx
	   << " ";
      if ( diminfo.IsPeriodic(dim) )
	{
	  cout << "T";
	}
      else
	{
	  cout << "F";
	}
      cout << "\n";
    }
  cout << "\n\n\n";

  cout << "Occupied Spatial Bin Information\n\n";

  cout << FMTW8 << "Bin"
       << FMTW12 << "Glb";
  for ( std::size_t dim=0; dim<ndim; ++dim)
    {
      std::stringstream ss;
      ss << "X" << dim+1;
      cout << FMTW5 << ss.str();
    }
  cout << FMTW8 << "N";
  for ( std::size_t dim=0; dim<ndim; ++dim)
    {
      std::stringstream ss;
      ss << "<X" << dim+1 << ">";
      cout << FMTW12 << ss.str();
    }
  
  cout << "\n";
  cout << FMTW8 << "idx"
       << FMTW12 << "idx";
  for ( std::size_t dim=0; dim<ndim; ++dim)
    {
      cout << FMTW5 << "idx";
    }
  cout << "\n";
  
  
  for ( std::size_t i=0; i<nbin; ++i )
    {
      cout << FMTW8 << i
	   << FMTW12 << sbins[i].glbidx;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  cout << FMTW5 << sbins[i].bidxs[dim];
	}
      cout << FMTW8 << sbinsamples[i].size()/2;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  cout << FMTF12 << sbins[i].center[dim];
	}
      cout << "\n";
    }
  cout << "\n";

  
#undef FMTW5
#undef FMTW6
#undef FMTW8
#undef FMTW12
#undef FMTF8
#undef FMTF12
}






double ndfes::SystemInfo::CptFastMBARObjective
( std::vector<double> const & bs,
  std::vector< std::vector<std::size_t> > const & idxs ) const
{
  double chisq = 0.;

  //
  // Evaluates Eq. (7) from
  // https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.8b01010
  // Ding, Vilseck, Brooks III. J. Chem. Theory. Comput. (2019) 15, 799-802.
  //

  //
  // \sum_j exp(z_j) = exp(zmax) \sum_j exp(z_j-zmax)
  //
  // Let, v = \sum_j exp(z_j-zmax), then
  //
  // \sum_j exp(z_j) = v exp(zmax)
  //
  // ln( \sum_j exp(z_j) ) = ln( v exp(zmax) )
  //                       = ln(v) + ln(exp(zmax))
  //                       = ln(v) + zmax
  //
  // In Eq. (7), z_j = -b_j - e_j
  // where "j" indexes a "state" and e_j is the potential+bias energy
  // of the state at the simulated sample.
  
  std::size_t const ndim = diminfo.GetNumDims();
  int const * dimisper = diminfo.GetIsPeriodic();
  std::size_t const totnpts = samples.size();
  std::size_t const nstate = states.size();
  
#ifdef WITH_OPENMP
#pragma omp parallel
  {
    double mychisq = 0.;
    std::vector<double> vs(nstate,0.);
#pragma omp for schedule(dynamic)
#else
    std::vector<double> vs(nstate,0.);
#endif
    for ( std::size_t i=0; i<nstate; ++i )
      {
	std::size_t const nupts = idxs[i].size()/2;
	
	for ( std::size_t iu=0; iu<nupts; ++iu )
	  {
	    std::size_t const ipt = idxs[i][0+iu*2];
	    std::size_t const degen = idxs[i][1+iu*2];
	    double const * pt = samples[ipt].GetPt();
	    double const * pU = samples[ipt].GetPotEnes();
	    double maxv = -1.e+100;
	    
	    for ( std::size_t j=0; j<nstate; ++j )
	      {
		double const * ks = states[j].GetFConst();
		double const * cs = states[j].GetCenter();
		std::size_t const enehamidx = states[j].GetHamIdx();
		double ene = pU[enehamidx]; // potential energy of state j
		for ( std::size_t dim=0; dim<ndim; ++dim )
		  {
		    double dx = pt[dim]-cs[dim];
		    if ( dimisper[dim] > 0 )
		      {
			dx = wrap(dx,360.);
		      }
		    ene += ks[dim] * dx * dx;
		  }
		// ene is the potential+bias energy of state j
		ene *= states[j].GetBeta();
		vs[j] = -bs[j]-ene;
		maxv = std::max(maxv,vs[j]);
	      }; // loop over energy states
	    double v = 0.;
	    for ( std::size_t j=0; j<nstate; ++j )
	      {
		v += std::exp( vs[j] - maxv );
	      }
#ifdef WITH_OPENMP
	    mychisq += degen * ((std::log(v+1.e-300)+maxv)/totnpts);
#else
	    chisq += degen * ((std::log(v+1.e-300)+maxv)/totnpts);
#endif
	  }; // loop over unique pts
	std::size_t const npts = statesamples[i].size()/2;
#ifdef WITH_OPENMP
	mychisq += (bs[i]*npts)/totnpts;
#else
	chisq += (bs[i]*npts)/totnpts;
#endif	
      }; // loop over simulated states
    
#ifdef WITH_OPENMP
#pragma omp critical
    {
      chisq += mychisq;
    }
  }
#endif

  return chisq;
}



double ndfes::SystemInfo::CptFastMBARObjective
( std::vector<double> const & bs,
  std::vector<double> & gs,
  std::vector< std::vector<std::size_t> > const & idxs ) const
{
  double chisq = 0.;

  //
  // Evaluates Eqs. (6) and (7) from
  // https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.8b01010
  // Ding, Vilseck, Brooks III. J. Chem. Theory. Comput. (2019) 15, 799-802.
  //

  
  std::size_t const ndim = diminfo.GetNumDims();
  int const * dimisper = diminfo.GetIsPeriodic();
  std::size_t const totnpts = samples.size();
  std::size_t const nstate = states.size();

  gs.assign( nstate, 0. );
  
#ifdef WITH_OPENMP
#pragma omp parallel
  {
    double mychisq = 0.;
    std::vector<double> mygs(nstate,0.);
    std::vector<double> vs(nstate,0.);
#pragma omp for schedule(dynamic)
#else
    std::vector<double> vs(nstate,0.);
#endif
    for ( std::size_t i=0; i<nstate; ++i )
      {
	std::size_t const nupts = idxs[i].size()/2;
	
	for ( std::size_t iu=0; iu<nupts; ++iu )
	  {
	    std::size_t const ipt = idxs[i][0+iu*2];
	    std::size_t const degen = idxs[i][1+iu*2];
	    double const * pt = samples[ipt].GetPt();
	    double const * pU = samples[ipt].GetPotEnes();
	    double maxv = -1.e+100;
	    
	    for ( std::size_t j=0; j<nstate; ++j )
	      {
		double const * ks = states[j].GetFConst();
		double const * cs = states[j].GetCenter();
		std::size_t const enehamidx = states[j].GetHamIdx();
		double ene = pU[enehamidx]; // potential energy of state j
		for ( std::size_t dim=0; dim<ndim; ++dim )
		  {
		    double dx = pt[dim]-cs[dim];
		    if ( dimisper[dim] > 0 )
		      {
			dx = wrap(dx,360.);
		      }
		    ene += ks[dim] * dx * dx;
		  }
		// ene is the potential+bias energy of state j
		ene *= states[j].GetBeta();
		vs[j] = -bs[j]-ene;
		maxv = std::max(maxv,vs[j]);
	      }; // loop over energy states
	    double v = 0.;
	    for ( std::size_t j=0; j<nstate; ++j )
	      {
		vs[j] = std::exp( vs[j] - maxv );
		v += vs[j];
	      }
#ifdef WITH_OPENMP
	    mychisq += degen * ((std::log(v+1.e-300)+maxv)/totnpts);
	    v *= totnpts;
	    for ( std::size_t j=0; j<nstate; ++j )
	      {
		mygs[j] -= degen*vs[j]/v;
	      }
#else
	    chisq += degen * ((std::log(v+1.e-300)+maxv)/totnpts);
	    v *= totnpts;
	    for ( std::size_t j=0; j<nstate; ++j )
	      {
		gs[j] -= degen*vs[j]/v;
	      }
#endif
	  }; // loop over unique pts
	std::size_t const npts = statesamples[i].size()/2;
	double const factor = ((double)npts)/totnpts;
#ifdef WITH_OPENMP
	mychisq += bs[i]*factor;
	mygs[i] += factor;
#else
	chisq += bs[i]*factor;
	gs[i] += factor;
#endif
      }; // loop over simulated states

#ifdef WITH_OPENMP
#pragma omp critical
    {
      chisq += mychisq;
      for ( std::size_t i=0; i<nstate; ++i )
	gs[i] += mygs[i];
    }
  }
#endif

    
  return chisq;
}



void ndfes::SystemInfo::CptWTP
( std::vector<double> const & Fin,
  std::vector< std::vector<std::size_t> > const & spts,
  std::vector< std::vector<std::size_t> > const & bpts,
  std::vector<double> & Fout,
  std::vector<double> & Fham,
  std::vector< std::vector<double> > & Fbins,
  std::vector< std::vector<double> > & Sbins ) const
{
  //
  //  f_t = -ln \sum_i { exp(-u_t) / \sum_j N_j exp(f_j-u_j) }
  //      = -ln exp(-C) \sum_i { exp(-u_t) / \sum_j N_j exp(f_j-u_j) exp(-C) }
  //      = -ln exp(-C) \sum_i { exp(-u_t) / \sum_j N_j exp(f_j-u_j-C) }
  //      = -ln exp(-C) \sum_i { exp(-u_t+D) / \sum_j N_j exp(f_j-u_j-C+D) }
  //      = C - ln(\sum_i { exp(-u_t+D) / \sum_j N_j exp(f_j-u_j-C+D) })
  //
  // Let C = max { f_i }
  // Let D = U_i (the potential energy of the hamiltonian
  //              producing the trajectory)


  //std::size_t simham_singleref = 1;
  
  std::size_t const ndim = diminfo.GetNumDims();
  int const * dimisper = diminfo.GetIsPeriodic();
  std::size_t const totnpts = samples.size();
  std::size_t const nstate = states.size();
  std::size_t const nbins = bpts.size();

  //
  // state2uham is a map that takes a biased-state index (corresponding to
  // an umbrella simulation) and returns a key-value to an petite array of unique
  // unbiased hamiltonians.  The "unique" index is not the hamiltonian index,
  // but an array index because the simulated hamiltonian indexes are not
  // necessarily 0-based on contiguous, whereas the unique indexes are.
  //
  // uham2simham takes a unique hamiltonian index and returns
  // the unbiased hamiltonian index
  //
  std::vector<std::size_t> state2uham;
  for ( std::size_t i=0; i<nstate; ++i )
    {
      state2uham.push_back( states[i].GetHamIdx() );
    }
  std::vector<std::size_t> const uham2simham( ndfes::CptHashTransform( state2uham ) );
  std::size_t const nsimham = uham2simham.size();
  // if ( sdosFullSampling )
  //   {
  //     nsimham = 1;
  //   }

  // for ( std::size_t uham=0; uham < nsimham; ++uham )
  //   {
  //     std::size_t simham=uham2simham[uham];
  //     std::prstd::size_tf("uham =%3i => simham=%3i\n",uham,simham);
  //   }
  
  // for ( std::size_t simham=0; simham < nstate; ++simham )
  //   {
  //     std::size_t uham=state2uham[simham];
  //     std::printf("state=%3i => uham  =%3i\n",simham,uham);
  //   }
  
  
  Fout.assign( nstate, 0. );
  Fham.assign( nham, 0. );
  Fbins.resize( nbins );
  for ( std::size_t i=0; i<nbins; ++i )
    {
      Fbins[i].assign( nham, 0. );
    }
  Sbins.resize( nbins );
  for ( std::size_t i=0; i<nbins; ++i )
    {
      Sbins[i].assign( nham, 0. );
    }
  
  // double C = -1.e+50;
  // for ( std::size_t i=0; i<nstate; ++i )
  //   {
  //     C = std::max(C,Fin[i]);
  //   }
  //std::vector<double> maxzs(totnpts,0.);
  //std::vector<double> Ps(totnpts,0.);
  std::vector<double> dPs(totnpts*nsimham,0.);

  //std::printf("BETA %20.10e %20.10e %20.10e\n",GetBeta(),states[0].GetBeta(),GetBeta()/states[0].GetBeta());

  // for ( std::size_t i=0; i<nstate; ++i )
  //   {
  //     std::printf("State Beta %6i %20.10e\n",i,states[i].GetBeta());
  //   }
  // std::printf("State Beta target %20.10e\n",GetBeta());

  // std::exit(0);

  std::vector<double> degens( totnpts, 0. );
  //std::vector<double> normexps( totnpts, 0. );
  std::vector<double> hamexps( totnpts*nham, 0. );
  std::vector<double> stateexps;
  std::size_t totnpts_nstate = totnpts * nstate;
  if ( totnpts_nstate < (std::size_t)(100000000) )
    {
      stateexps.assign( totnpts_nstate, 0. );
    };



#ifdef WITH_OPENMP
#pragma omp parallel
  {
#endif
    std::vector<double> BiasedUs(nstate,0.);
    std::vector<double> zs(nstate,0.);
    std::vector<double> ns(nstate,0.);
#ifdef WITH_OPENMP
#pragma omp for
#endif
    for ( std::size_t i=0; i<nstate; ++i )
      {
	std::size_t const nupts = spts[i].size()/2;
	//std::size_t const simhamidx = states[i].GetHamIdx();
	
	for ( std::size_t iu=0; iu<nupts; ++iu )
	  {
	    std::size_t const ipt = spts[i][0+iu*2];
	    std::size_t const degen = spts[i][1+iu*2];
	    double const * pt = samples[ipt].GetPt();
	    double const * pU = samples[ipt].GetPotEnes();
	    
	    degens[ipt] = degen;
	    double maxz = -1.e+30;
	    for ( std::size_t j=0; j<nstate; ++j )
	      {
		double const * ks = states[j].GetFConst();
		double const * cs = states[j].GetCenter();
		std::size_t const enehamidx = states[j].GetHamIdx();
		double bias = 0.;
		for ( std::size_t dim=0; dim<ndim; ++dim )
		  {
		    double dx = pt[dim]-cs[dim];
		    if ( dimisper[dim] > 0 )
		      {
			dx = wrap(dx,360.);
		      }
		    bias += ks[dim] * dx * dx;
		  }

		// {
		//   if ( enehamidx == 0 )
		//     {
		//       double t0 = (pU[enehamidx] + bias) * GetBeta();
		//       double t1 = (pU[enehamidx] + bias) * states[j].GetBeta();
		//       double tbias = (t1-t0)/GetBeta();
		//       std::printf("%8i %10.2f\n",iu,tbias);
		//     }
		// }
		
		BiasedUs[j] = (pU[enehamidx] + bias) * states[j].GetBeta();
		zs[j]       = Fin[j] - BiasedUs[j];
		maxz        = std::max(maxz,zs[j]);
		double const Nj = statesamples[j].size()/2;
		ns[j] = Nj;
		
		// std::printf("%5i %3i %3i %13.4e %13.4e %13.4e\n",
		// 	    ipt,j,enehamidx,
		// 	    pU[enehamidx]*states[j].GetBeta(),
		// 	    pU[enehamidx]* GetBeta(),
		// 	    pU[enehamidx]*(states[j].GetBeta()-GetBeta()));
		
		if ( sdosFullSampling )
		  {
		    for ( std::size_t uham=0; uham < nsimham; ++uham )
		      {
			std::size_t simham = uham2simham[uham];
			double ene = BiasedUs[j] - pU[simham] * GetBeta();
			dPs[uham+ipt*nsimham] += Nj * std::exp( Fin[j] - ene );
			//dPs[uham+ipt*nsimham] += Nj * std::exp( Fin[j] - bias - (pU[enehamidx]-pU[simham]) );
		      }
		  }
		else
		  {
		    std::size_t uham = state2uham[samples[ipt].GetStateIdx()];
		    std::size_t simham = uham2simham[uham];
		    if ( simham == enehamidx )
		      {
			double ene = BiasedUs[j] - pU[simham] * GetBeta();
			dPs[uham+ipt*nsimham] += Nj * std::exp( Fin[j] - ene );
			//dPs[uham+ipt*nsimham] += Nj * std::exp( Fin[j] - bias );
		      }
		  }

	      }; // energy states

	    

	    double normexp = ndfes::logsumexp( nstate, ns.data(), zs.data() );
	    if ( stateexps.size() > 0 )
	      {
		for ( std::size_t j=0; j<nstate; ++j )
		  {
		    stateexps[ipt+j*totnpts] = -BiasedUs[j] - normexp;
		  }
	      };
	    
	    for ( std::size_t j=0; j<nham; ++j )
	      {
		hamexps[ipt+j*totnpts] = -GetBeta()*pU[j] - normexp;
	      }
	    
	    //double den = 0.;
	    //for ( int j=0; j<nstate; ++j )
	    //  {
	    //	double const Nj = statesamples[j].size()/2;
	    //	den += Nj * std::exp( zs[j]-maxz );
	    //  }
	    //Ps[ipt] = 1./den;

	    for ( std::size_t uham=0; uham<nsimham; ++uham )
	      {
		if ( dPs[uham+ipt*nsimham] > 0. )
		  {
		    dPs[uham+ipt*nsimham] = 1./dPs[uham+ipt*nsimham];
		  }
		else
		  {
		    dPs[uham+ipt*nsimham] = 0.;
		  }
	      }

// 	    maxzs[ipt] = maxz;

// #ifdef WITH_OPENMP
// #pragma omp critical
// #endif
// 	    {
// 	      for ( std::size_t j=0; j<nstate; ++j )
// 		{
// 		  Fout[j] += degen * std::exp( -BiasedUs[j] - maxz ) / den;
// 		}
	      
// 	      for ( std::size_t j=0; j<nham; ++j )
// 		{
// 		  std::printf("U=%20.10e  %20.10e\n",-GetBeta()*pU[j] - MaxMinusUs[j],degen/den);
// 		  Fham[j] += degen * std::exp( -(GetBeta()*pU[j]) - MaxMinusUs[j] - maxz )/den;
// 		}
// 	    }
	    

	    
	  }; // unique pts
      }; // simulation states
#ifdef WITH_OPENMP
  }
#endif

  // X = \sum_i ( (x[i]+A) - y[i] )^2 = \sum_i (yp[i]-y[i])^2 = \sum_i yp[i]^2+y[i]^2-2*y[i]*yp[i]
  // dX/dA = \sum_i 2*yp[i]-2*y[i] = 0 = \sum_i yp[i]-y[i] = (\sum_i (x[i]-y[i])) + N*A
  // N*A = \sum_i y[i]-x[i]
  // A = <y-x>

  // dF = Fout-Fin = -ln(h)-Fin
  // exp(-dF) = exp(ln(h)+Fin) = h*exp(Fin)

  if ( stateexps.size() > 0 )
    {
      double mue = 0.;
      for ( std::size_t i=0; i<nstate; ++i )
	{
	  Fout[i] = - ndfes::logsumexp( totnpts, degens.data(), stateexps.data() + i*totnpts );
	  //Fout[i] = - std::log(Fout[i]+1.e-300);
	  mue += Fin[i]-Fout[i];
	};
      mue /= nstate;
      for ( std::size_t i=0; i<nstate; ++i )
	{
	  Fout[i] += mue;
	}
    }
  else
    {
      for ( std::size_t i=0; i<nstate; ++i )
	{
	  Fout[i] = Fin[i];
	}
    }
  

  for ( std::size_t i=0; i<nham; ++i )
    {
      //std::printf("Fham[i] = %20.10e %20.10e %20.10e\n",Fham[i],std::log(Fham[i]),MaxMinusUs[i]);
      //Fham[i] = - (std::log(Fham[i]) + MaxMinusUs[i]) ; // * (GetBeta()/states[0].GetBeta());
      Fham[i] = - ndfes::logsumexp( totnpts, degens.data(), hamexps.data() + i*totnpts );
      for ( std::size_t ipt=0; ipt<totnpts; ++ipt )
	{
	  hamexps[ipt+i*totnpts] += Fham[i];
	}
    };


  bool const do_dos = deltaDoS > 1.e-5;



#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
  for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      // Loop over spatial bins
      std::vector<std::size_t> const & upts = bpts[ibin];
      std::size_t const nupts = upts.size()/2;

      //std::printf("ibin = %i %i\n",ibin,nupts);

      std::vector< std::vector<std::size_t> > ebins(nsimham);
      for ( std::size_t uham=0; uham<nsimham; ++uham )
	{
	  ebins[uham].assign( nupts, 0 );
	}
      std::vector<double> Psums(nsimham,0.);
      std::vector<double> Uavgs(nsimham,0.);
      std::vector<double> Uvars(nsimham,0.);
      std::vector< std::vector<double> > rho_s(nsimham);
      std::vector< std::vector<double> > rho_g(nsimham);
      
      for ( std::size_t iham=0; iham<nham; ++iham )
	{
	  double S = 0.;
	  //double F = 0.;

	  if ( do_dos )
	    {
	      for ( std::size_t uham=0; uham<nsimham; ++uham )
		{
		  std::fill( ebins[uham].data(), ebins[uham].data() + nupts, 0 );
		};
	      std::fill( Psums.data(), Psums.data() + nsimham, 0. );
	      std::fill( Uavgs.data(), Uavgs.data() + nsimham, 0. );
	      std::fill( Uvars.data(), Uvars.data() + nsimham, 0. );

	      for ( std::size_t iu=0; iu<nupts; ++iu )
		{
		  std::size_t const ipt     = upts[0+iu*2];
		  std::size_t const degen   = upts[1+iu*2];
		  double const * pU = samples[ipt].GetPotEnes();

		  if ( sdosFullSampling )
		    {
		      for ( std::size_t uham=0; uham<nsimham; ++uham )
			{
			  std::size_t simham = uham2simham[uham];
			  double const dU = GetBeta()*(pU[iham]-pU[simham]);
			  Psums[uham] += degen * dPs[uham+ipt*nsimham];
			  Uavgs[uham] += degen * dPs[uham+ipt*nsimham] * dU;
			  ebins[uham][iu] = ndfes::GetEneIdx(dU,deltaDoS);
			};
		    }
		  else
		    {
		      std::size_t const stateidx= samples[ipt].GetStateIdx();
		      std::size_t const uham = state2uham[stateidx];
		      std::size_t const simham= uham2simham[uham];
		      double const dU = GetBeta()*(pU[iham]-pU[simham]);
		      Psums[uham] += degen * dPs[uham+ipt*nsimham];
		      Uavgs[uham] += degen * dPs[uham+ipt*nsimham] * dU;
		      ebins[uham][iu] = ndfes::GetEneIdx(dU,deltaDoS);
		    }

		  // if ( simham == 1 and iham == 3 )
		  // {
		  //   std::printf("dU %8i iham=%3i simham=%3i uham=%3i dPs=%14.5e dU=%14.5e Psums=%14.5e Uavgs=%14.5e\n",ipt,iham,simham,uham,dPs[ipt],dU,Psums[uham],Uavgs[uham]);
		  // }
		}

	      std::vector< std::vector<std::size_t> > ekeys(nsimham);
	      
	      for ( std::size_t uham=0; uham < nsimham; ++uham )
		{
		  ekeys[uham] = ndfes::CptHashTransform(ebins[uham]);
		  rho_s[uham].assign( ekeys[uham].size(), 0. );
		  rho_g[uham].assign( ekeys[uham].size(), 0. );

		  // if ( sbins[ibin].glbidx == 1 )
		  //   {
		  //     std::printf("Psum %2i %20.10e %20.10e %20.10e\n",uham,Psums[uham],Uavgs[uham],Uavgs[uham]/Psums[uham]);
		  //   }
		  if ( Psums[uham] > 1.e-300 )
		    {
		      Uavgs[uham] /= Psums[uham];
		    }
		  else
		    {
		      Uavgs[uham] = 0.;
		    }
		}
	      
	      for ( std::size_t iu=0; iu<nupts; ++iu )
		{
		  std::size_t const ipt     = upts[0+iu*2];
		  std::size_t const degen   = upts[1+iu*2];
		  double const * pU = samples[ipt].GetPotEnes();
		  if ( sdosFullSampling )
		    {
		      for ( std::size_t uham=0; uham < nsimham; ++uham )
			{
			  std::size_t const simham  = uham2simham[uham];
			  double const dU   = GetBeta()*(pU[iham]-pU[simham]);
			  double const dx   = dU - Uavgs[uham];
			  double const fact = degen * dPs[simham+ipt*nsimham]/Psums[uham];
			  Uvars[uham] += fact * dx*dx;
			  rho_s[uham][ebins[uham][iu]] += fact;
			}
		    }
		  else
		    {
		      std::size_t const stateidx= samples[ipt].GetStateIdx();
		      std::size_t const uham = state2uham[stateidx];
		      std::size_t const simham= uham2simham[uham];
		      double const dU = GetBeta()*(pU[iham]-pU[simham]);
		      double const dx   = dU - Uavgs[uham];
		      double const fact = degen * dPs[simham+ipt*nsimham]/Psums[uham];
		      Uvars[uham] += fact * dx*dx;
		      rho_s[uham][ebins[uham][iu]] += fact;
		    }
		  
		}

	      
	      
	      for ( std::size_t uham=0; uham<nsimham; ++uham )
		{
		  double const avg = Uavgs[uham];
		  std::size_t const nekeys = rho_g[uham].size();
		  if ( Uvars[uham] < 1.e-8 )
		    {
		      for ( std::size_t ukey=0; ukey<nekeys; ++ukey )
			{
			  rho_g[uham][ukey] = rho_s[uham][ukey];
			}
		    }
		  else
		    {
		      // if ( sbins[ibin].glbidx == 1 )
		      // 	{
		      // 	  std::printf("iham,uham,Uvars %3i %3i %14.4e %14.4e\n",iham,uham,Uavgs[uham],Uvars[uham]);
		      // 	};

		      double const stddev = std::sqrt( Uvars[uham] );
		      for ( std::size_t ukey=0; ukey<nekeys; ++ukey )
			{
			  double const a = (ekeys[uham][ukey] - 0.5) * deltaDoS;
			  double const b = a + deltaDoS;
			  double const t1 = std::erf( M_SQRT1_2 * (b-avg)/stddev );
			  double const t2 = std::erf( M_SQRT1_2 * (a-avg)/stddev );
			  rho_g[uham][ukey] = std::max( 0.5*(t1-t2), 1.e-10 );
			  // if ( sbins[ibin].glbidx == 1 )
			  //   {
			  //     std::printf("iham,uham,ukey %3i %3i %3i %14.4e %14.4e %14.4e %14.4e %14.4e\n",iham,uham,ukey,rho_g[uham][ukey],M_SQRT1_2 * (b-avg)/stddev,M_SQRT1_2 * (a-avg)/stddev,a,b);
			  //   };
			};
		    }
		}


	      if ( (std::size_t)(sdosHistPrefix.size()) > 0 )
		{
		  for ( std::size_t uham=0; uham<nsimham; ++uham )
		    {
		      std::size_t simham = uham2simham[uham];
		      
		      std::size_t const nekeys = rho_s[uham].size();
		      //if ( Uvars[uham] > 1.e-8 )
			
		      std::stringstream fname;
		      fname << sdosHistPrefix
			    << "hist"
			//<< std::fixed
			//<< std::setprecision(3)
			//<< sbins[ibin].center[0]
			    << "_bin_" << ibin
			    << "_target_" << iham
			    << "_ref_" << simham
			    << ".dat";
		      std::ofstream fout;
		      fout.open( fname.str().c_str() );
		      if ( ! fout )
			{
			  std::cerr << "Unable to open " << fname.str()
				    << " for writing" << std::endl;
			}
#define FMT std::scientific << std::setw(20) << std::setprecision(10)
		      fout << "# Uavg: " << FMT << Uavgs[uham]
			   << " Uvar: " << FMT << Uvars[uham]
			   << " RC: ";
		      for ( std::size_t idim=0; idim<sbins[ibin].center.size(); ++idim )
			{
			  fout << std::fixed << std::setprecision(3)
			       << sbins[ibin].center[idim];
			}
		      fout << std::endl;
		      
		      for ( std::size_t ukey=0; ukey<nekeys; ++ukey )
			{
			  double const a = ekeys[uham][ukey] * deltaDoS;
			  fout << FMT << a
			       << FMT << rho_s[uham][ukey]
			       << FMT << rho_g[uham][ukey]
			       << std::endl;
#undef FMT
			}
		      fout.close();
		    }
		  
		}

	      

	    };
	  
	  //F = 0.;
	  S = 0.;
	  std::size_t nbinpts = 0;
	  std::vector<double> zs(nupts,0.);
	  std::vector<double> ps(nupts,0.);
	  std::vector<double> as(nupts,0.);
	  //double maxz = -1.e+300;
	  for ( std::size_t iu=0; iu<nupts; ++iu )
	    {
	      std::size_t const ipt     = upts[0+iu*2];
	      std::size_t const degen   = upts[1+iu*2];
	      nbinpts += degen;
	      //std::size_t const simham  = samples[ipt].GetHamIdx();
	      //double const * pU = samples[ipt].GetPotEnes();
	      //zs[iu] = Fham[iham] - GetBeta()*pU[iham] - maxzs[ipt];
	      //maxz = std::max(maxz,zs[iu]);
	      zs[iu] = hamexps[ipt+iham*totnpts];

	      double alpha = 1.0;
	      if ( do_dos )
		{
		  std::size_t const stateidx  = samples[ipt].GetStateIdx();
		  std::size_t uham    = state2uham[stateidx];
		  std::size_t const ukey    = ebins[uham][iu];
		  if ( rho_s[uham][ukey] > 0. )
		    {
		      alpha = rho_g[uham][ukey] / rho_s[uham][ukey];
		    };
		  alpha = std::min(10.0,alpha);
		};
	      ps[iu] = alpha * degen;
	      as[iu] = alpha;
	    };
	  double norm = ndfes::logsumexp( nupts, ps.data(), zs.data() );
	  if ( nbinpts > 1 )
	    {
	      //double wsum = 0.;
	      for ( std::size_t iu=0; iu<nupts; ++iu )
		{
		  std::size_t const degen   = upts[1+iu*2];
		  double w = as[iu] * std::exp( zs[iu] - norm );
		  //wsum += degen * w;
		  if ( nbinpts > 1 and w > 1.e-100 )
		    {
		      S -= degen * w * std::log(w);
		    }
		};
	      S /= std::log((double)nbinpts);
	      //std::printf("Wsum = %20.10e\n",wsum);
	    }
	  
	  
	  // for ( std::size_t iu=0; iu<nupts; ++iu )
	  //   {
	  //     std::size_t const ipt     = upts[0+iu*2];
	  //     std::size_t const degen   = upts[1+iu*2];

	  //     double h = as[iu] * std::exp( zs[iu] - maxz ) * Ps[ipt];
          //     F += degen * h;
          //     S += degen * h * ( maxz + std::log(h+1.e-300) );
	      
	  //   };
	  // if ( nbinpts > 1 and std::abs(F) > 1.e-100 )
	  //   {
	  //     S = - ( S/F - ( maxz + std::log(F+1.e-300) ) ) / std::log((double)nbinpts);
	  //   }
	  // else
	  //   {
	  //     S = 0.;
	  //   }

	  Fbins[ibin][iham] = -norm;
	  Sbins[ibin][iham] = S;
	  
	} // iham
      
    } // ibin
  
}





void ndfes::SystemInfo::WriteChkpt_MBAR
( std::string const fname,
  std::vector< std::vector<double> > const & Fbins,
  std::vector< std::vector<double> > const & Ebins,
  std::vector< std::vector<double> > const & Sbins ) const
{

  std::size_t const nbin = Fbins.size();
  std::size_t const ndim = diminfo.GetNumDims();
  
  if ( Fbins.size() != sbinsamples.size() )
    {
      std::cerr << "Error: SystemInfo spatial bin size mismatch (Fbins) "
		<< Fbins.size() << " != " << sbinsamples.size()
		<< std::endl;
      std::exit(1);						
    };

  if ( Ebins.size() != sbinsamples.size() )
    {
      std::cerr << "Error: SystemInfo spatial bin size mismatch (Ebins) "
		<< Ebins.size() << " != " << sbinsamples.size()
		<< std::endl;
      std::exit(1);						
    };

  if ( Sbins.size() != sbinsamples.size() )
    {
      std::cerr << "Error: SystemInfo spatial bin size mismatch (Sbins) "
		<< Sbins.size() << " != " << sbinsamples.size()
		<< std::endl;
      std::exit(1);						
    };

  
  
  std::ofstream cout;
  cout.open( fname.c_str() );
  
  std::size_t indent = 4;

#define INDENT std::setw(indent) << ""
#define I6 std::setw(6)
#define FMTE std::scientific << std::setw(23) << std::setprecision(14)
  
  cout << "#!/usr/bin/env python3\n\n"
       << "import ndfes\n\n";

  diminfo.WriteChkpt(cout,0);

  
  for ( std::size_t iham=0; iham<nham; ++iham )
    {
      cout << "\nmodel" << iham
	   << " = ndfes.MBAR(grid,\n"
	   << INDENT << "{ ";
      std::size_t ilast = 0;
      for ( std::size_t ibin=0; ibin<nbin; ++ibin )
	{
	  if ( Ebins[ibin][iham] >= 0 )
	    {
	      ilast = ibin;
	    };
	};
      
      bool first = true;
      for ( std::size_t ibin=0; ibin<nbin; ++ibin )
	{
	  double avg = ConvertKTToKcal(Fbins[ibin][iham]);
	  double err = ConvertKTToKcal(Ebins[ibin][iham]);
	  double ent = Sbins[ibin][iham];
	  {
	    if ( ! first )
	      {
		cout << INDENT << "  ";
	      }
	    cout << I6 << sbins[ibin].glbidx << " :";	
	    first=false;
	    cout << " ndfes.SpatialBin([";
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		cout << I6 << sbins[ibin].bidxs[dim];
		if ( dim < ndim-1 )
		  {
		    cout << ",";
		  }
	      }
	    cout << " ],"
		 << FMTE << avg << ","
		 << FMTE << err << ","
		 << FMTE << ent << ","
		 << I6 << sbinsamples[ibin].size()/2 << " )";
	    if ( ibin < ilast )
	      {
		cout << ",\n";
	      }
	    else
	      {
		cout << "})\n";
	      }
	  }
	};
    }
  
  
  cout << "\nmodels = [";
  for ( std::size_t iham=0; iham<nham; ++iham )
    {
      cout << "model" << iham;
      if ( iham < nham-1 )
	{
	  cout << ",";
	}
    };
  cout << "]\n\n";
  
  
#undef INDENT
#undef I6
#undef FMTE
  
}






void ndfes::SystemInfo::WriteChkptXML_MBAR
( std::string const fname,
  std::vector< std::vector<double> > const & Fbins,
  std::vector< std::vector<double> > const & Ebins,
  std::vector< std::vector<double> > const & Sbins ) const
{

  std::size_t const nbin = Fbins.size();
  std::size_t const ndim = diminfo.GetNumDims();
  
  if ( Fbins.size() != sbinsamples.size() )
    {
      std::cerr << "Error: SystemInfo spatial bin size mismatch (Fbins) "
		<< Fbins.size() << " != " << sbinsamples.size()
		<< std::endl;
      std::exit(1);						
    };

  if ( Ebins.size() != sbinsamples.size() )
    {
      std::cerr << "Error: SystemInfo spatial bin size mismatch (Ebins) "
		<< Ebins.size() << " != " << sbinsamples.size()
		<< std::endl;
      std::exit(1);						
    };

  if ( Sbins.size() != sbinsamples.size() )
    {
      std::cerr << "Error: SystemInfo spatial bin size mismatch (Sbins) "
		<< Sbins.size() << " != " << sbinsamples.size()
		<< std::endl;
      std::exit(1);						
    };

  
  
  std::ofstream cout;
  cout.open( fname.c_str() );
  
  std::size_t indent = 2;

#define INDENT(LEV) std::setw(indent*LEV) << ""
#define I6 std::setw(6)
#define FMTE  std::scientific << std::setw(23) << std::setprecision(14)
#define FMTF  std::fixed << std::setw(13) << std::setprecision(8)

  //diminfo.WriteChkpt(cout,0);

  cout << "<ndfes>\n";
  
  
  for ( std::size_t iham=0; iham<nham; ++iham )
    {
      cout << INDENT(1)
	   << "<model idx=\"" << iham << "\">\n";
      cout << INDENT(2)
	   << "<type>MBAR</type>\n";
      {
	double const * xmins = diminfo.GetXmin();
	double const * xmaxs = diminfo.GetXmax();
	int const * sizes = diminfo.GetDimSizes();
	int const * isper = diminfo.GetIsPeriodic();
	
	cout << INDENT(2)
	     << "<grid>\n";
	
	for ( std::size_t idim=0; idim<ndim; ++idim )
	  {
	    cout << INDENT(3)
		 << "<dim idx=\"" << idim << "\">\n"
		 << INDENT(4)
		 << "<xmin>" << FMTF << xmins[idim] << "</xmin>\n"
		 << INDENT(4)
		 << "<xmax>" << FMTF << xmaxs[idim] << "</xmax>\n"
		 << INDENT(4)
		 << "<size>" << sizes[idim] << "</size>\n"
		 << INDENT(4)
		 << "<isper>" << isper[idim] << "</isper>\n"
		 << INDENT(3)
		 << "</dim>\n";
	  }
	cout << INDENT(2)
	     << "</grid>\n";
      }
      
      for ( std::size_t ibin=0; ibin<nbin; ++ibin )
	{
	  double avg = ConvertKTToKcal(Fbins[ibin][iham]);
	  double err = ConvertKTToKcal(Ebins[ibin][iham]);
	  double ent = Sbins[ibin][iham];
	  {
	    cout << INDENT(2)
		 << "<bin idx=\"" << sbins[ibin].glbidx << "\">\n";
	    for ( std::size_t idim=0; idim<ndim; ++idim )
	      {
		cout << INDENT(3)
		     << "<bidx idx=\"" << idim << "\">"
		     << sbins[ibin].bidxs[idim]
		     << "</bidx>\n";
	      }
	    cout << INDENT(3)
		 << "<val>" << FMTE << avg << "</val>\n"
		 << INDENT(3)
		 << "<err>" << FMTE << err << "</err>\n"
		 << INDENT(3)
		 << "<re> " << FMTE << ent << "</re>\n"
		 << INDENT(3)
		 << "<size>" << sbinsamples[ibin].size()/2 << "</size>\n";
	  }
	  cout << INDENT(2) << "</bin>\n";
	}
      cout << INDENT(1) << "</model>\n";
    }
  cout << "</ndfes>\n";

  
  
#undef INDENT
#undef I6
#undef FMTE
#undef FMTF

}




std::size_t ndfes::SystemInfo::InitvFEP( std::size_t const nsubbins )
{
  vfep.reset( new ndfes::vFEPData( diminfo, nsubbins, states, sbins ) );
  return vfep->glbcidxs.size();
}



std::vector<double> ndfes::SystemInfo::CptvFEPLinearTerm
( std::vector< std::vector<std::size_t> > const & bpts ) const
{
  if ( ! vfep )
    {
      std::cerr << "Called ndfes::SystemInfo::CptvFEPLinearTerm "
		<< "without first calling ndfes::SystemInfo::InitvFEP"
		<< std::endl;
      std::exit(EXIT_FAILURE);
    }
  return vfep->CptLinearTerm( samples, bpts, sbins );
}
    

double ndfes::SystemInfo::CptvFEPObjective
( std::vector<double> const & bs,
  std::vector<double> & gs,
  std::vector<double> const & cornerh ) const
{
  return vfep->CptChisq( bs, gs, cornerh );
}


std::size_t ndfes::SystemInfo::GetParamCornerGlbIdx( std::size_t ic ) const
{
  if ( ! vfep )
    {
      std::cerr << "Called ndfes::SystemInfo::GetParamCornerGlbIdx "
		<< "without first calling ndfes::SystemInfo::InitvFEP"
		<< std::endl;
      std::exit(EXIT_FAILURE);
    }
  else if ( ic >= (std::size_t)vfep->glbcidxs.size() )
    {
      std::cerr << "Error: ndfes::SystemInfo::GetParamCornerGlbIdx "
		<< "parameter index is out of bounds"
		<< std::endl;
      std::exit(EXIT_FAILURE);
    }
  return vfep->glbcidxs[ic];
}




void ndfes::SystemInfo::InterpvFEP
( double const * c,
  std::vector<double> const & p,
  std::vector<double> const & dp,
  double & val,
  double & err ) const
{
  if ( ! vfep )
    {
      std::cerr << "Called ndfes::SystemInfo::InterpvFEP "
		<< "without first calling ndfes::SystemInfo::InitvFEP"
		<< std::endl;
      std::exit(EXIT_FAILURE);
    }
  vfep->InterpvFEP(c,p,dp,val,err);
}




void ndfes::SystemInfo::WriteChkpt_vFEP
( std::string const fname,
  std::vector<double> const & Fs,
  std::vector<double> const & Es ) const
{

  std::size_t const nparam = Fs.size();
  std::size_t const ndim = diminfo.GetNumDims();
  std::size_t const nbins = sbins.size();

  if ( ! vfep )
    {
      std::cerr << "Error: ndfes::SystemInfo::WriteChkpt_vFEP "
		<< "called without initiating vFEP"
		<< std::endl;
      std::exit(1);
    }
  
  if ( Fs.size() != vfep->glbcidxs.size() )
    {
      std::cerr << "Error: SystemInfo vfep parameter size mismatch (Fs) "
		<< Fs.size() << " != " << vfep->glbcidxs.size()
		<< std::endl;
      std::exit(1);						
    };

  if ( Es.size() != Fs.size() )
    {
      std::cerr << "Error: SystemInfo vfep parameter size mismatch (Es) "
		<< Es.size() << " != " << Fs.size()
		<< std::endl;
      std::exit(1);						
    };


  
  
  std::ofstream cout;
  cout.open( fname.c_str() );
  
  std::size_t indent = 4;

#define INDENT std::setw(indent) << ""
#define I6 std::setw(6)
#define FMTE std::scientific << std::setw(23) << std::setprecision(14)
  
  cout << "#!/usr/bin/env python3\n\n"
       << "import ndfes\n\n";

  diminfo.WriteChkpt(cout,0);

  
  for ( std::size_t iham=0; iham<nham; ++iham )
    {
      cout << "\nmodel" << iham
	   << " = ndfes.vFEP(grid,\n"
	   << INDENT << "{ ";

      
      bool first = true;
      std::size_t ilast = nbins-1;
      for ( std::size_t ibin=0; ibin<nbins; ++ibin )
	{
	  //double avg = ConvertKTToKcal(Fbins[ibin][iham]);
	  //double err = ConvertKTToKcal(Ebins[ibin][iham]);
	  //double ent = Sbins[ibin][iham];
	  {
	    if ( ! first )
	      {
		cout << INDENT << "  ";
	      }
	    cout << I6 << sbins[ibin].glbidx << " :";	
	    first=false;
	    cout << " ndfes.SpatialBin([";
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		cout << I6 << sbins[ibin].bidxs[dim];
		if ( dim < ndim-1 )
		  {
		    cout << ",";
		  }
	      }
	    cout << " ]" << ", size="
		 << I6 << sbinsamples[ibin].size()/2 << ")";
	    if ( ibin < ilast )
	      {
		cout << ",\n";
	      }
	    else
	      {
		cout << "}, "
		     << diminfo.GetBsplOrder()
		     << ",\n";
	      }
	  }
	};
      cout << INDENT << "{ ";
      ilast = nparam-1;
      first=true;
      for ( std::size_t iparam=0; iparam<nparam; ++iparam )
	{
	  if ( ! first )
	    {
	      cout << INDENT << "  ";
	    }
	  first=false;
	  cout << I6 << vfep->glbcidxs[iparam]
	       << " : ( "
	       << FMTE << ConvertKTToKcal(Fs[iparam]) << ","
	       << FMTE << ConvertKTToKcal(Es[iparam]) << " )";
	  if ( iparam < ilast )
	    {
	      cout << ",\n";
	    }
	  else
	    {
	      cout << "})\n";
	    }
	};
    }
  
  
  cout << "\nmodels = [";
  for ( std::size_t iham=0; iham<nham; ++iham )
    {
      cout << "model" << iham;
      if ( iham < nham-1 )
	{
	  cout << ",";
	}
    };
  cout << "]\n\n";
  
  
#undef INDENT
#undef I6
#undef FMTE
  
}






void ndfes::SystemInfo::WriteChkptXML_vFEP
( std::string const fname,
  std::vector<double> const & Fs,
  std::vector<double> const & Es ) const
{

  std::size_t const nparam = Fs.size();
  std::size_t const ndim = diminfo.GetNumDims();
  std::size_t const nbins = sbins.size();

  if ( ! vfep )
    {
      std::cerr << "Error: ndfes::SystemInfo::WriteChkptXML_vFEP "
		<< "called without initiating vFEP"
		<< std::endl;
      std::exit(1);
    }
  
  if ( Fs.size() != vfep->glbcidxs.size() )
    {
      std::cerr << "Error: SystemInfo vfep parameter size mismatch (Fs) "
		<< Fs.size() << " != " << vfep->glbcidxs.size()
		<< std::endl;
      std::exit(1);						
    };

  if ( Es.size() != Fs.size() )
    {
      std::cerr << "Error: SystemInfo vfep parameter size mismatch (Es) "
		<< Es.size() << " != " << Fs.size()
		<< std::endl;
      std::exit(1);						
    };


  
  
  std::ofstream cout;
  cout.open( fname.c_str() );
  
  std::size_t indent = 2;

#define INDENT(LEV) std::setw(indent*LEV) << ""
#define I6 std::setw(6)
#define FMTE std::scientific << std::setw(23) << std::setprecision(14)
#define FMTF std::scientific << std::setw(13) << std::setprecision(8)

  cout << "<ndfes>\n";

  
  for ( std::size_t iham=0; iham<nham; ++iham )
    {
      cout << INDENT(1)
	   << "<model idx=\"" << iham << "\">\n";
      cout << INDENT(2)
	   << "<type>VFEP</type>\n"
	   << INDENT(2)
	   << "<order>" << diminfo.GetBsplOrder() << "</order>\n";

      {
	double const * xmins = diminfo.GetXmin();
	double const * xmaxs = diminfo.GetXmax();
	int const * sizes = diminfo.GetDimSizes();
	int const * isper = diminfo.GetIsPeriodic();
	
	cout << INDENT(2)
	     << "<grid>\n";
	
	for ( std::size_t idim=0; idim<ndim; ++idim )
	  {
	    cout << INDENT(3)
		 << "<dim idx=\"" << idim << "\">\n"
		 << INDENT(4)
		 << "<xmin>" << FMTF << xmins[idim] << "</xmin>\n"
		 << INDENT(4)
		 << "<xmax>" << FMTF << xmaxs[idim] << "</xmax>\n"
		 << INDENT(4)
		 << "<size>" << sizes[idim] << "</size>\n"
		 << INDENT(4)
		 << "<isper>" << isper[idim] << "</isper>\n"
		 << INDENT(3)
		 << "</dim>\n";
	  }
	cout << INDENT(2)
	     << "</grid>\n";
      }
      
      for ( std::size_t ibin=0; ibin<nbins; ++ibin )
	{
	  //double avg = ConvertKTToKcal(Fbins[ibin][iham]);
	  //double err = ConvertKTToKcal(Ebins[ibin][iham]);
	  //double ent = Sbins[ibin][iham];
	  {

	    cout << INDENT(2)
		 << "<bin idx=\"" << sbins[ibin].glbidx << "\">\n";
	    for ( std::size_t idim=0; idim<ndim; ++idim )
	      {
		cout << INDENT(3)
		     << "<bidx idx=\"" << idim << "\">"
		     << sbins[ibin].bidxs[idim]
		     << "</bidx>\n";
	      }
	    
	    cout << INDENT(3)
		 << "<size>" << sbinsamples[ibin].size()/2 << "</size>\n";

	  }
	  cout << INDENT(2) << "</bin>\n";
	}
      
      for ( std::size_t iparam=0; iparam<nparam; ++iparam )
	{
	  cout << INDENT(2)
	       << "<corner idx=\"" << vfep->glbcidxs[iparam] << "\">\n"
	       << INDENT(3)
	       << "<val>" 
	       << FMTE << ConvertKTToKcal(Fs[iparam])
	       << "</val>\n"
	       << INDENT(3)
	       << "<err>" 
	       << FMTE << ConvertKTToKcal(Es[iparam])
	       << "</err>\n"
	       << INDENT(2)
	       << "</corner>\n";
	};
      cout << INDENT(1) << "</model>\n";
    }
    cout << "</ndfes>\n";

  
#undef INDENT
#undef I6
#undef FMTE
  
}














void ndfes::SystemInfo::CheckForConsistentTemperatures( bool expert ) const
{                       
  std::size_t ns = states.size();
  std::vector<double> temps;
  for ( std::size_t i=0; i<ns; ++i )
    {
      double T = states[i].GetTemperature();
      if ( std::abs( T - GetTargetTemp() ) > 0.001 )
	{	  
	  std::vector<double>::iterator p =
	    std::find(temps.begin(),temps.end(),T);
	  
	  if ( p == temps.end() )
	    {
	      temps.push_back(T);
	    }
	}
    }

  if ( (std::size_t)temps.size() > 0 )
    {

      std::stringstream msg1;
      msg1 << "There are " << temps.size()
	   << " simulation temperatures in the metafile (";
      for ( std::size_t k=0, n=temps.size(); k<n; ++k )
	{
	  msg1 << std::setw(7) << std::setprecision(3) << std::fixed
	       << temps[k];
	  if ( k < n-1 )
	    {
	      msg1 << ",";
	    }
	}
      msg1 << ") that differ from the target temperature ("
	   << std::setw(7) << std::setprecision(3) << std::fixed
	   << GetTargetTemp()
	   << ").";

      std::stringstream msg2;
      msg2 << "This will likely give inaccurate "
	   << "results. This is an experimental feature,"
	   << " and its use is not recommended!\n\n";
      


      
      if ( vfep )
	{
	  std::cerr << "\nERROR: " << msg1.str()
		    << "\n--vfep does not currently support "
		    << "multiple temperatures."
		    << std::endl;
	  std::exit(1);
	}
      else
	{
	  if ( expert )
	    {
	      std::cout << "\nWARNING: " << msg1.str() << " " << msg2.str();
	    }
	  else
	    {
	      std::cerr << "\nERROR: " << msg1.str() << " " << msg2.str();
	      std::cerr << "You likely either:\n"
			<< "  1. set an incorrect temperature on the command line\n"
			<< "  2. you forgot to set the correct temperature on the command line\n"
			<< "  3. you set an incorrect temperature in the metafile\n"
			<< "But, if you're an expert and you REALLY want to analyze multiple\n"
			<< "temperatures, then you can do so by using the --expert flag\n";
	      std::exit(1);
	    }
	}
    }
  
}
