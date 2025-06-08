#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <cstdlib>

#include "ReadMetafile.hpp"
#include "GetBeta.hpp"

namespace ndfes
{
  ndfes::pDimInfo CreateDimInfo
  ( std::size_t const ndim,
    std::size_t const BsplOrder,
    std::vector<double> const & TargetWidths,
    std::vector<int> const & PeriodicDims );

  int ReadSamples
  ( std::string dfile,
    double StartFrac,
    double StopFrac,
    double beta,
    std::size_t const istate,
    std::size_t const ndim,
    std::size_t const nham,
    std::size_t const iham,
    std::vector<ndfes::Sample> & samples );
    
}


double CalcAvg( std::vector<double> const & a )
{
  double s=0.;
  std::size_t const n = a.size();
  for ( std::size_t i=0; i<n; ++i )
    {
      s += a[i];
    }
  return s/n;
}


std::vector<std::size_t> CalcIntAvgs( std::vector< std::vector<double> > const & a )
{
  std::size_t const n = a.size();
  std::vector<std::size_t> b(n,0);
  for ( std::size_t i=0; i<n; ++i )
    {
      double c = CalcAvg(a[i]);
      b[i] = c+0.5;
    };
      
  return b;
}



std::size_t CalcStatisticalInefficiency( std::size_t N, double const * A )
{
  std::size_t const mintime = 3;


  //std::printf("dats[0:1] %.8f %.8f\n",A[0],A[1]);
  
  double g = 1.0;

  double mu = 0.;
  for ( std::size_t i=0; i<N; ++i )
    mu += A[i];
  mu /= N;

  std::vector<double> dA(A,A+N);
  for ( std::size_t i=0; i<N; ++i )
    dA[i] -= mu;

  double sigma2_AB = 0.;
  for ( std::size_t i=0; i<N; ++i )
    sigma2_AB += dA[i]*dA[i];
  sigma2_AB /= N;

  //std::printf("N,sig2 %i %.8f\n",N,sigma2_AB);
  
  if ( std::abs(sigma2_AB) > 1.e-8 )
    {
      for ( std::size_t t=1; t < N-1; t += 1 )
        {
          double C = 0.;
          for ( std::size_t j=0; j < N-t; ++j )
            C += dA[j]*dA[j+t];
          C /= ( sigma2_AB * ( N-t ) );

          if ( C <= 0.0 and t > mintime )
            break;

          g += 2.0 * C * (1.0 - ((double)t)/((double)N));
        }
    }
  //std::printf("g=%.6f sigma2_AB=%.6f\n",g,sigma2_AB);
  if ( g < 1.0 )
    g = 1.0;
  std::size_t blocksize = g+0.5;
  return blocksize;
}



ndfes::pDimInfo ndfes::CreateDimInfo
( std::size_t const ndim,
  std::size_t const BsplOrder,
  std::vector<double> const & TargetWidths,
  std::vector<int> const & PeriodicDims )
{
  std::vector<int> dimisper(ndim,0);
  for ( std::size_t i=0, n=PeriodicDims.size(); i<n; ++i )
    {
      int k = PeriodicDims[i]-1;
      if ( k >= 0 and k < (int)ndim )
	{
	  dimisper[k] = 1;
	}
      else
	{
	  std::cerr << "Invalid command line argument\n";
	  std::cerr << "Periodic dim index "
		    << PeriodicDims[i]
		    << " is out of range [1," << ndim
		    << "]" << std::endl;
	  std::exit(1);
	};
    }

  std::vector<double> twidths;
  if ( (std::size_t)TargetWidths.size() > 0 )
    {
      twidths.assign(ndim,TargetWidths[0]);
      if ( (std::size_t)TargetWidths.size() == ndim )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      twidths[dim] = TargetWidths[dim];
	    }
	}
      else if ( (std::size_t)TargetWidths.size() == 1 )
	{
	  std::cout << "Using bin width "
		    << std::fixed
		    << std::setw(9)
		    << std::setprecision(4)
		    << twidths[0]
		    << "\n";
	}
      else
	{
	  std::cerr << "Invalid command line argument\n";
	  std::cerr << "The number of specified bin widths -w is "
		    << TargetWidths.size()
		    << ", but there are "
		    << ndim
		    << " dimensions"
		    << std::endl;
	  std::exit(1);
	}
    }
  else
    {
      std::cerr << "Error in ndfes ReadMetafile: "
		<< "TargetWidths has no elements"
		<< std::endl;
      std::exit(1);
    }
	      
  return ndfes::pDimInfo( new ndfes::DimInfo( ndim, BsplOrder, dimisper.data(), twidths.data() ) );
}




int ndfes::ReadSamples
( std::string dfile,
  double StartFrac,
  double StopFrac,
  double beta,
  std::size_t const istate,
  std::size_t const ndim,
  std::size_t const nham_read,
  std::size_t const iham,
  std::vector<ndfes::Sample> & samples )
{
  //std::cout << "ReadSamples nham = " << nham << " " << iham << std::endl;
  
  int retcode = 0;
  samples.resize(0);

  StartFrac = std::max(0.,StartFrac);
  StopFrac = std::min(1.,StopFrac);
  if ( StartFrac >= StopFrac )
    {
      std::cerr << "Error: StartFrac >= StopFrac will discard all samples"
		<< std::endl;
      std::exit(1);
    };
  
  std::size_t nham = std::max((int)nham_read,1);

  std::ifstream cin;
  cin.open( dfile.c_str() );
  if ( ! cin.good() )
    {
      retcode = -1;
    }
  else
    {
      
      //int nham_read = nham;
      
      /*
      if ( nham_read < 2 )
	{
	  nham_read = 0;
	};
      */

      std::size_t mincols = 1 + ndim + nham_read;
      std::size_t maxcols = 1 + ndim + nham_read;
      if ( nham_read == 0 )
	{
	  maxcols += 1;
	}

      //std::printf("ReadSamples nham_read %i\n",nham_read);
      
      //int nskip=0;
      double prevt = -1.;
      std::size_t tmp_npts = 0;
      std::vector<double> tmp_obspts;
      std::vector<double> tmp_potenes;
      std::string line;

      while ( std::getline( cin, line ) )
	{
	  std::vector<std::string> tokens;
	  std::istringstream iss(line);
	  std::copy( std::istream_iterator<std::string>(iss),
		     std::istream_iterator<std::string>(),
		     std::back_inserter(tokens) );
	  if ( tokens.size() == 0 )
	    {
	      continue;
	    }
	  else if ( (std::size_t)tokens.size() < mincols )
	    {
	      std::cerr << "Error: Not enough columns for ndim="
			<< ndim << " and nham=" << nham_read
			<< std::endl;
	      std::cerr << "In file " << dfile << std::endl;
	      std::cerr << "On line " << line << std::endl;
	      std::exit(1);
	    }
	  else if ( (std::size_t)tokens.size() > maxcols )
	    {
	      std::cerr << "Error: Too many columns for ndim="
			<< ndim << " and nham=" << nham_read
			<< std::endl; 
	      std::cerr << "In file " << dfile << std::endl;
	      std::cerr << "On line " << line << std::endl;
	      std::exit(1);
	    }
	  else
	    {
	      {
		double x;
		std::stringstream sstr(tokens[0]);
		sstr >> x;
		if ( std::abs(x-prevt) < 0.0001 )
		  {
		    //++nskip;
		    continue; // skip repeated time value
		  };
		prevt = x;
	      }
	      for ( std::size_t i=1; i < ndim+1; ++i )
		{
		  double x;
		  std::stringstream sstr(tokens[i]);
		  sstr >> x;
		  tmp_obspts.push_back(x);
		};
	      for ( std::size_t i=0; i<nham_read; ++ i )
		{
		  double x;
		  std::stringstream sstr(tokens[1+ndim+i]);
		  sstr >> x;
		  tmp_potenes.push_back(beta*x);
		}

	      tmp_npts++;
	    };
	};

      if ( tmp_npts == 0 )
	{
	  retcode = -2;
	}
      
      int istart = (int)(tmp_npts * StartFrac + 0.5);
      int istop = (int)(tmp_npts * StopFrac + 0.5);
      istart = std::max(0,istart);
      istop = std::min((int)tmp_npts,istop);

      std::vector<double> zeros(nham,0.);
      double const * pU = zeros.data();
      for ( int i=istart; i<istop; ++i )
	{
	  if ( nham_read > 0 )
	    {
	      pU = tmp_potenes.data() + i*nham;
	    }
	  samples.push_back( ndfes::Sample(istate,ndim,
					   tmp_obspts.data() + i*ndim,
					   nham, iham,
					   pU ) );
	}

      if ( (std::size_t)samples.size() == 0 )
	{
	  retcode = -3;
	}
      
    }
  return retcode;
}





ndfes::pSystemInfo ndfes::ReadMetafile
( std::string fname,
  std::size_t const NumHam,
  std::size_t const BsplOrder,
  std::vector<double> const & TargetWidths,
  std::vector<int> const & PeriodicDims,
  double const TargetT,
  double const StartFrac,
  double const StopFrac,
  double const DeltaDoS,
  bool const sdosFullSampling,
  std::string sdosHistPrefix,
  bool const keeponlyuniq,
  int const maxuniq )
{
  ndfes::pDimInfo diminfo;
  ndfes::pSystemInfo sysinfo;


  std::vector< std::vector<double> > statineff;

  std::size_t nham_read = NumHam;
  /*
  if ( nham_read < 2 )
    {
      nham_read = 0;
    }
  */

  //std::printf("nham_read=%i\n",nham_read);
  
  // find basename
  std::string basename;
  {
    std::size_t p = fname.find_last_of('/');
    if ( p != std::string::npos )
      {
	basename = fname.substr(0,p+1);
      }
  }

  
  std::ifstream cin;
  cin.open( fname.c_str() );

  if ( ! cin.good() )
    {
      std::cerr << "Error: Could not open "
		<< fname << std::endl;
      std::exit(1);
    };
   

  bool found_a_line = false;
  std::string line;
  while( std::getline(cin,line) )
    {
      std::vector<std::string> tokens;
      std::istringstream iss(line);
      std::copy( std::istream_iterator<std::string>(iss),
                 std::istream_iterator<std::string>(),
                 std::back_inserter(tokens) );
      if ( tokens.size() == 0 )
        {
          continue;
        }
      else
	{
	  found_a_line=true;
	  int ndim2 = (tokens.size()-3);
	  if ( ndim2 % 2 )
	    {
	      std::cerr << "Incorrect number of columns in file "
			<< fname << " on line:\n" << line << "\n";
	      std::exit(1);
	    };
	  
	  std::size_t ndim = ndim2/2;
	  std::vector<double> biasloc(ndim,0.);
	  std::vector<double> biasfcs(ndim,0.);

	  
	  std::size_t hamidx = std::atoi(tokens[0].c_str());

	  if ( hamidx >= NumHam and nham_read > 0 )
	    {
	      std::cerr << "Error: Hamiltonian index "
			<< hamidx << " is out of range [0,"
			<< NumHam-1 << "]\n";
	      std::cerr << "In file " << fname << std::endl;
	      std::cerr << "On line " << line << std::endl;
	      std::cerr << "Perhaps you meant to specify --nham="
			<< hamidx+1
			<< " on the command line?"
			<< std::endl;
	      std::cerr << "Or perhaps this is a typo in the metafile."
			<< std::endl;
	      std::exit(1);
	    }
	  
	  double simtemp = std::atof(tokens[1].c_str());
	  if ( std::abs(simtemp - TargetT) > 0.001 and nham_read < 1 )
	    {
	      std::cerr << "Error:  The simulation temperatures do not match the target\n"
			<< "temperature set on the command line.   You likely forgot to\n"
		        << "set the correct temperature with the --temp flag.  Or maybe\n"
		        << "you set the temperature incorrectly in the metafile. If you\n"
		        << "are  an  expert  and  you're  actually  trying  to  analyze\n"
			<< "simulations performed at different temperatures,  be  aware\n"
			<< "that  this  is an  experimental feature.  It  often  yields\n"
			<< "inaccurate results - but if you insist, then you can try it\n"
			<< "by including the unbiased energies of each potential  (even\n"
			<< "if you have only 1 potential)  as additional columns in the\n"
			<< "dumpave  files.  You'll  then need to set  --nham > 0.  You\n"
			<< "didn't  set  a  --nham  to  a  positive  integer,  however.\n\n";
	      std::exit(1);
	    }
	  //double simbeta = ndfes::GetBeta(simtemp);
	  std::string dumpave = tokens[2];
	  for ( std::size_t i=0; i<ndim; ++i )
	    {
	      biasloc[i] = std::atof(tokens[3+2*i].c_str());
	      biasfcs[i] = std::atof(tokens[4+2*i].c_str());
	    };

	  
	  if ( ! diminfo )
	    {
	      diminfo = ndfes::CreateDimInfo( ndim, BsplOrder, TargetWidths, PeriodicDims );
	    }
	  else
	    {
	      if ( diminfo->GetNumDims() != ndim )
		{
		  std::cerr << "Error: Expected "
			    << diminfo->GetNumDims()
			    << " dimensions, but found "
			    << ndim
			    << " on line:\n"
			    << line
			    << std::endl;
		  std::exit(1);
		};
	    }

	  if ( ! sysinfo )
	    {
	      sysinfo.reset( new ndfes::SystemInfo( *diminfo, TargetT, DeltaDoS, sdosFullSampling, sdosHistPrefix, std::max((int)NumHam,1) ) );
	    }

	  std::size_t istate = sysinfo->InsertState(hamidx,biasloc.data(),biasfcs.data(),simtemp);


	  
	  for ( std::size_t itry=0; itry<2; ++itry )
	    {	      
	      std::string dfile = dumpave;
	      if ( itry == 1 )
		{
		  dfile = basename + dumpave;
		}

	      std::vector<ndfes::Sample> samples;
	      
	      //int retcode = ndfes::ReadSamples(dfile,StartFrac,StopFrac,simbeta,istate,ndim,std::max(NumHam,1),hamidx, samples );
	      //std::printf("NumHam,max %i %i\n",NumHam,std::max(NumHam,1));
	      int retcode = ndfes::ReadSamples(dfile,StartFrac,StopFrac,1.,istate,ndim,NumHam,hamidx, samples );
	      
	      
	      if ( retcode == 0 )
		{
		  //
		  // extract statistically independent samples
		  //
		  ndfes::State state( sysinfo->GetState(istate) );
		  std::vector<double> biase(samples.size(),0.);
		  for ( std::size_t i=0, n=samples.size(); i<n; ++i )
		    {
		      biase[i] = state.CptBiasEnergy( diminfo->GetIsPeriodic(), samples[i].GetPt() );
		    }
		  std::size_t g = CalcStatisticalInefficiency( biase.size(), biase.data() );
		  //std::printf("%3i\n",g);

		  if ( keeponlyuniq )
		    {
		      std::vector<ndfes::Sample> usamples;
		      for ( std::size_t i=0, n=samples.size(); i<n; i += g )
			{
			  usamples.push_back( samples[i] );
			}
		      if ( maxuniq > 0 and (int)usamples.size() > maxuniq )
			{
			  usamples.resize(maxuniq);
			}
		      samples = usamples;
		      g=1;
		    }
		  //else
		  //{
		  //  int gold = state.GetStatIneff();
		  //  g = std::max(gold,g);
		  //};
		  //sysinfo->SetStatIneff(istate,g);

		  if ( istate+1 > statineff.size() )
		    {
		      statineff.resize(istate+1);
		    }
		  statineff[istate].push_back( (double)g );
		  
		  
		  for ( std::size_t i=0, n=samples.size(); i<n; ++i )
		    {
		      sysinfo->InsertSample(samples[i]);
		    }
		  break;
		}
	      else if ( retcode == -2 )
		{
		  std::cerr << "Error: No samples were read from file "
			    << dfile << std::endl;
		  std::exit(1);
		}
	      else if ( retcode == -3 )
		{
		  std::cerr << "Error: Samples were read from file "
			    << dfile << "\n"
			    << "but the user-specified start/stop range "
			    << "discarded all samples"
			    << std::endl;
		  std::exit(1);
		}
	      else if ( itry == 0 ) // retcode == -1
		{
		  if ( (int)basename.size() == 0 )
		    {
		      std::cerr << "Error: Could not read file "
				<< dumpave
				<< std::endl;
		      std::exit(1);
		    }
		}
	      else if ( itry == 1 )
		{
		  std::cerr << "Error: Could not read file "
			    << dumpave
			    << " nor "
			    << dfile
			    << std::endl;
		  std::exit(1);
		}
	    }
	  
	};

    };


  std::vector<std::size_t> gs( CalcIntAvgs( statineff ) );
  for ( std::size_t istate=0, ns=sysinfo->GetNumStates(); istate<ns; ++istate )
    {
      sysinfo->SetStatIneff(istate,gs[istate]);
    }

  if ( ! found_a_line )
    {
      std::cerr << "Error: No text within " << fname << "\n";
      std::exit(1);
    };

  sysinfo->ShiftPotEnesToAvg();
  sysinfo->CptSpatialHistogram();

  return sysinfo;
}
