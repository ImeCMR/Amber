#include <iostream>
#include <cstdlib>
#include <iomanip>

#include "../ndfeslib/ReadMetafile.hpp"
#include "../ndfeslib/OptMBAR.hpp"
#include "../ndfeslib/OptvFEP.hpp"
#include "../ndfeslib/RunAvg.hpp"

#include "options.hpp"



void RunMBAR( ndfes::options const & cli,
	      ndfes::pSystemInfo sysinfo );

void RunVFEP( ndfes::options const & cli,
	      ndfes::pSystemInfo sysinfo,
	      std::size_t const nparam );






int main( int argc, char * argv[] )
{  
  ndfes::options cli( ndfes::ReadOptions( argc, argv ) );

  std::cout << "Reading metafile: " << cli.metafile << std::endl;

  ndfes::pSystemInfo sysinfo =
    ndfes::ReadMetafile( cli.metafile, cli.nham, cli.order,
			 cli.target_widths, cli.perdims,
			 cli.temp, cli.startfrac, cli.stopfrac,
			 cli.sdos, cli.sdosFullSampling, cli.sdosHistPrefix,
			 cli.uniq, cli.maxuniq );


  sysinfo->CheckForConsistentTemperatures( cli.expert );
  
  
  if ( ! sysinfo )
    {
      std::cerr << "Error: Failed to read metafile" << std::endl;
      std::exit(1);
    }
  
  std::size_t const ns = sysinfo->GetNumStates();

  if ( ns < 1 )
    {
      std::cerr << "Error: No states have been read" << std::endl;
      std::exit(1);
    };

  sysinfo->PrintSummary( std::cout );

  if ( cli.mbar )
    {
      RunMBAR( cli, sysinfo );
    }
  else if ( cli.vfep )
    {
      std::size_t nparam = sysinfo->InitvFEP( cli.nquad );
      RunVFEP( cli, sysinfo, nparam );
    }
  
  
  return 0;
}










void RunMBAR
( ndfes::options const & cli,
  ndfes::pSystemInfo sysinfo )
{

  typedef std::vector< std::vector<std::size_t> > vvint;
  typedef std::vector< std::vector<double> > vvdouble;

  std::size_t const ns = sysinfo->GetNumStates();

  vvint StateSamples;
  vvint BinSamples;
  sysinfo->GetSampleIdxs( StateSamples, BinSamples );

      
  std::vector<double> StateFs(ns,0.);
  ndfes::OptMBAR( *sysinfo, StateSamples, cli.reltol,
		  cli.maxiter, true, StateFs );
      
      
  std::size_t const nham = std::max(cli.nham,1);
  std::size_t const nbin = sysinfo->GetNumBins();
  std::vector<double> refHamFs(nham,0.);
  vvdouble refBinFs(nbin);
  vvdouble refBinSs(nbin);

  std::vector<ndfes::RunAvg> avgStateFs(ns);
  std::vector<ndfes::RunAvg> avgHamFs(nham);
  std::vector< std::vector<ndfes::RunAvg> > avgBinFs(nbin);
  for ( std::size_t i=0; i<nbin; ++i )
    {
      avgBinFs[i].resize(nham);
    }

  for ( std::size_t iter=0; iter<(std::size_t)(1+cli.nboot); ++iter )
    {
      std::vector<double> tmpStateFs( StateFs );
      std::vector<double> myStateFs(ns,0.);
      std::vector<double> myHamFs(nham,0.);
      vvdouble myBinFs(nbin);
      vvdouble myBinSs(nbin);
      vvint myStateSamples(StateSamples);
      vvint myBinSamples(BinSamples);
	  
      if ( iter > 0 )
	{
	  sysinfo->SetHistPrefix("");
	  std::cout << "Bootstrap "
		    << std::setw(4) << iter << " / "
		    << std::setw(4) << cli.nboot << "\n";
	  
	  sysinfo->ResampleStates( myStateSamples,myBinSamples );
	  if ( cli.btol > 0 )
	    {
	      ndfes::OptMBAR( *sysinfo, myStateSamples, cli.btol,
			      cli.maxiter, false, tmpStateFs );
	    };
	}
	  
      sysinfo->CptWTP( tmpStateFs, myStateSamples, myBinSamples,
		       myStateFs, myHamFs,
		       myBinFs, myBinSs );

      if ( iter == 0 )
	{
	  refHamFs = myHamFs;
	  refBinFs = myBinFs;
	  refBinSs = myBinSs;
	}

      //#ifdef WITH_OPENMP
      //#pragma omp critical
      //#endif
      {
	for ( std::size_t i=0; i<ns; ++i )
	  {
	    avgStateFs[i].push_back( myStateFs[i] );
	  }
	for ( std::size_t i=0; i<nham; ++i )
	  {
	    avgHamFs[i].push_back( myHamFs[i] );
	  }
	for ( std::size_t ibin=0; ibin<nbin; ++ibin )
	  {
	    if ( (std::size_t)myBinSamples[ibin].size()/2 > 0 )
	      {
		for ( std::size_t iham=0; iham<nham; ++iham )
		  {
		    avgBinFs[ibin][iham].push_back( myBinFs[ibin][iham] );
		  }
	      };
	  }
      } // critical
	  
    }
      
	  
#define FMT std::fixed << std::setw(11) << std::setprecision(4)

  std::cout << "\nBiased free energy of each simulated state\n\n";
  std::cout << std::setw(6)  << "State"
	    << std::setw(11) << "F"
	    << std::setw(11) << "<F>"
	    << std::setw(11) << "dF"
	    << "\n";
  for ( std::size_t i=0; i<ns; ++i )
    {
      std::cout << std::setw(6) << i
		<< FMT << sysinfo->GetState(i).ConvertKTToKcal(StateFs[i])
		<< FMT << sysinfo->GetState(i).ConvertKTToKcal(avgStateFs[i].mean())
		<< FMT << sysinfo->GetState(i).ConvertKTToKcal(avgStateFs[i].stderr(true))
		<< "\n";
    };
      
  std::cout << "\n\n";

  std::cout << "Unbiased free energy of each Hamiltonian\n\n";
  std::cout << std::setw(6)  << "Ham"
	    << std::setw(11) << "F"
	    << std::setw(11) << "<F>"
	    << std::setw(11) << "dF"
	    << "\n";
  for ( std::size_t i=0; i<nham; ++i )
    {
      std::cout << std::setw(6) << i
		<< FMT << sysinfo->ConvertKTToKcal(refHamFs[i])
		<< FMT << sysinfo->ConvertKTToKcal(avgHamFs[i].mean())
		<< FMT << sysinfo->ConvertKTToKcal(avgHamFs[i].stderr(true))
		<< "\n";
    }
      
  // std::cout << "\n\n";
      
  // for ( std::size_t i=0; i<nbin; ++i )
  //   {
	    
  //     double const * c = sysinfo->GetBin(i).center.data();
  //     std::cout << std::fixed << std::setw(15) << std::setprecision(8)
  // 		<< c[0];
	  
  //     for ( std::size_t iham=0; iham < nham; ++iham )
  // 	{
  // 	  double stderr = avgBinFs[i][iham].stderr(true);
  // 	  std::cout << FMT << sysinfo->ConvertKTToKcal(refBinFs[i][iham])
  // 		    << FMT << sysinfo->ConvertKTToKcal(avgBinFs[i][iham].mean())
  // 		    << FMT << sysinfo->ConvertKTToKcal(stderr);
  // 	  std::cout << FMT << refBinSs[i][iham];
  // 	}
  //     std::cout << "\n";
  //   }
      
      
#undef FMT

  { // Write Checkpoint
    std::vector< std::vector<double> > BinEs(nbin);
    for ( std::size_t i=0; i<nbin; ++i )
      {
	BinEs[i].assign(nham,0.);
	for ( std::size_t iham=0; iham<nham; ++iham )
	  {
	    double err = avgBinFs[i][iham].stderr(true);
	    if ( cli.nboot > 1 and avgBinFs[i][iham].size() < 2 )
	      {
		err = sysinfo->ConvertKcalToKT( 1. );
	      };
	     
	    BinEs[i][iham] = err;
	  }
      }

    std::cout << "\n\nWriting bin free energies to " << cli.chkpt << "\n\n";
    
    //sysinfo->WriteChkpt_MBAR( cli.chkpt, refBinFs, BinEs, refBinSs );
    sysinfo->WriteChkptXML_MBAR( cli.chkpt, refBinFs, BinEs, refBinSs );
  }

      

}








void RunVFEP
( ndfes::options const & cli,
  ndfes::pSystemInfo sysinfo,
  std::size_t const nparam )
{
  std::size_t const nham = std::max(cli.nham,1);
  //std::size_t const ns = sysinfo->GetNumStates();
  //std::size_t const nbin = sysinfo->GetNumBins();

  if ( nham > 1 )
    {
      std::cerr << "Error: vFEP cannot be run with multiple Hamiltonians"
		<< std::endl;
      std::exit(1);
    }
  
  typedef std::vector< std::vector<std::size_t> > vvint;
  //typedef std::vector< std::vector<double> > vvdouble;

  vvint StateSamples;
  vvint BinSamples;
  sysinfo->GetSampleIdxs( StateSamples, BinSamples );

  std::vector<double> refParams(nparam,0.);

  //ndfes::OptvFEP_CheckGradients( nparam, *sysinfo, BinSamples, refParams );
  
  ndfes::OptvFEP( nparam, *sysinfo, BinSamples,
		  cli.reltol, cli.maxiter, true, refParams );
  
  std::vector<ndfes::RunAvg> avgParams(nparam);


  for ( std::size_t iter=0; iter<(std::size_t)(1+cli.nboot); ++iter )
    {
      std::vector<double> params( refParams );
      vvint myStateSamples(StateSamples);
      vvint myBinSamples(BinSamples);
	  
      if ( iter > 0 )
	{
	  std::cout << "Bootstrap "
		    << std::setw(4) << iter << " / "
		    << std::setw(4) << cli.nboot << "\n";
	  
	  sysinfo->ResampleStates( myStateSamples,myBinSamples );
	  if ( cli.btol > 0 )
	    {
	      ndfes::OptvFEP( nparam, *sysinfo, myBinSamples,
			      cli.btol, cli.maxiter,
			      false, params );
	    };
	}
	
      //#ifdef WITH_OPENMP
      //#pragma omp critical
      //#endif
      {
	for ( std::size_t i=0; i<nparam; ++i )
	  {
	    avgParams[i].push_back( params[i] );
	  }
      } // critical
	  
    }
      
	  

  std::vector<double> refErrors(nparam,0.);
  for ( std::size_t i=0; i<nparam; ++i )
    {
      refErrors[i] = avgParams[i].stderr(true);
    }

  
#define FMT std::fixed << std::setw(11) << std::setprecision(4)
  
  std::cout << "\nvFEP parameters, bootstrap means and standard errors\n\n";
  std::cout << std::setw(6)  << "Corner"
	    << std::setw(11) << "F"
	    << std::setw(11) << "<F>"
	    << std::setw(11) << "dF"
    //<< std::setw(6)  << "Nobs"
	    << "\n";
  std::cout << std::setw(6)  << "Idx"
	    << "\n";
  for ( std::size_t i=0; i<nparam; ++i )
    {
      std::cout << std::setw(6) << sysinfo->GetParamCornerGlbIdx(i)
		<< FMT << sysinfo->ConvertKTToKcal(refParams[i])
		<< FMT << sysinfo->ConvertKTToKcal(avgParams[i].mean())
		<< FMT << sysinfo->ConvertKTToKcal(refErrors[i])
	//<< std::setw(6) << avgParams[i].size()
		<< "\n";
    };


  /*
  std::size_t const nbin = sysinfo->GetNumBins();
  std::size_t const ndim = sysinfo->GetDimInfo().GetNumDims();
  std::cout << "\n\nvFEP potential at the bin centers\n\n"
	    << std::setw(6) << "Bin"
	    << std::setw(11) << "F"
	    << std::setw(11) << "dF"
	    << "\n"
	    << std::setw(6) << "Idx"
	    << "\n";
  for ( std::size_t ibin=0; ibin<nbin; ++ibin )
    {
      double const * c = sysinfo->GetBin(ibin).center.data();
      double val = 0.;
      double err = 0.;
      sysinfo->InterpvFEP( c, refParams, refErrors, val, err );
      std::cout << std::setw(6) << ibin;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  std::cout << FMT << c[dim];
	};
      std::cout << FMT << sysinfo->ConvertKTToKcal(val)
		<< FMT << sysinfo->ConvertKTToKcal(err)
		<< "\n";
    };
  */
      
      
#undef FMT

  
  std::cout << "\n\nWriting bin free energies to " << cli.chkpt << "\n\n";
  //sysinfo->WriteChkpt_vFEP( cli.chkpt, refParams, refErrors );
  sysinfo->WriteChkptXML_vFEP( cli.chkpt, refParams, refErrors );

  
}



