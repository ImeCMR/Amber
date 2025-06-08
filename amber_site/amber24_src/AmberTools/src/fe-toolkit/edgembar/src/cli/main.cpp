#include <iomanip>
#include <fstream>
#include <cmath>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "cli_options.hpp"
#include "../edgembar/ReadInput.hpp"
#include "../edgembar/Edge.hpp"
#include "../edgembar/StateValues.hpp"
#include "../edgembar/CalcIdx.hpp"

static double get_beta( double T )
{
  double const k_au = 3.16681534222374e-06;
  double const au_per_kcal = 1.59360145069066e-03;
  double const k_kcal = ( k_au / au_per_kcal );
  //std::printf("%20.10e\n",1./( k_kcal * T ) );
  return 1. / ( k_kcal * T );
}


std::string strip_suffix( std::string str, std::string const suffix )
{
  if ( str.size() >= suffix.size() )
    {
      if (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)
	{
	  //std::cout << "Match " << str << " " << suffix <<"\n";
	  str = str.substr( 0, str.size() - suffix.size() );
	  //std::cout << "Subst " << str << " " << suffix <<"\n";
	}
      else
	{
	  //std::cout << "No match " << str << " " << suffix <<"\n";
	}
    }
  return str;
}


void WriteOutput( std::ostream & cout,
		  std::string const & body )
{
  //std::string pystr = pdata->GetPython( beta, calcs, GlbValues );
  cout << "#!/usr/bin/env python3\n"
       << "import edgembar\n"
       << "import numpy as np\n\n"
       << "edge = " << body << ";\n\n"
       << "if __name__ == '__main__':\n"
       << "    import edgembar\n"
       << "    edgembar.ProcessEdgeCli(edge)\n";

}


void PrintParams( edgembar::Edge const * e,
		  std::vector<double> const & x,
		  double const beta )
{
  std::vector< edgembar::Trial const * > trials( e->GetConstTrials() );
  for ( int i=0, n=trials.size(); i<n; ++i )
    {
      int nstates = trials[i]->GetNumStates();
      int o = trials[i]->GetParamOffset();
      for ( int j=0; j<nstates; ++j )
	{
	  std::cout << std::setw(12) << std::setprecision(4)
		    << std::fixed << x[o+j] / beta
		    << " " << std::setw(3) << j << " "
		    << trials[i]->GetEdge()->GetName()
		    << " "
		    << trials[i]->GetEnv()->GetName()
		    << " "
		    << trials[i]->GetStage()->GetName()
		    << " "
		    << trials[i]->GetName()
		    << "\n";
	};
    };
}







int main( int argc, char * argv[] )
{
  typedef std::shared_ptr<edgembar::Edge> edge_t;
  typedef edgembar::CalcIdx calc_t;
  
  edgembar::cli_options cli( edgembar::read_options( argc, argv ) );
  
  double beta = get_beta( cli.temp );
  
  edge_t cdata( ReadInput( cli.inpfile, beta, cli.readmode,
			   cli.fstart, cli.fstop, cli.stride,
			   cli.ptol, cli.autoeqmode ) );

  if ( cli.uwts )
    {
      cdata->UseUniformWts();
    }
  
  // if ( cli.fstart > 0 or cli.fstop < 1 )
  //   {
  //     cdata = cdata->ExtractRange( cli.fstart, cli.fstop );
  //   }

  cdata->ApplyMinimalConstraints();

  edge_t pdata;
  if ( cli.autoequil )
    {
      pdata = cdata->ExtractProd();
    }
  else
    {
      pdata = cdata->ExtractCopy();
    }

  int nparam = cdata->GetNumParam();
  std::vector<double> pguess(nparam,0);
  if ( not cli.nullguess )
    {
      pguess = cdata->MakeExpAvgGuess();
    };
  
  std::vector<double> popt(nparam,0);
  std::vector<double> fopt(nparam,0);
  pdata->Optimize( pguess, popt, fopt, cli.tol, cli.verbosity );
  pguess = popt;
  double fe0 = pdata->GetFreeEnergy( fopt );


  std::vector<edge_t> bases;
  std::vector<calc_t> calcs;

  //
  // Setup the unconstrained production calculation
  //
  {
    int baseidx = bases.size();
    bases.push_back(pdata);
    calc_t idx;
    idx.BaseIdx = baseidx;
    idx.IsCon = false;
    idx.IsBootstrap = false;
    calcs.push_back(idx);
  }

  
  //
  // Setup spline calculation
  //
  {
    double df = cli.dcon * beta;
    double lo = fe0 - df;
    double hi = fe0 + df;
    int nbase = 2*cli.ncon+1;
    if ( nbase > 1 )
      {
	for ( int i=0; i < nbase; ++i )
	  {
	    double conval = lo + (hi-lo)*i/( (float)(nbase-1) );
	    edge_t e( pdata->ExtractCopy() );
	    e->ApplyConstraint( conval );

	    int baseidx = bases.size();
	    bases.push_back(e);
	    calc_t idx;
	    idx.BaseIdx = baseidx;
	    idx.IsCon = true;
	    idx.ConVal = conval;
	    idx.IsBootstrap = false;
	    calcs.push_back(idx);
	  }
      }
    /*
    else
      {
	calc_t idx;
	idx.BaseIdx = bases.size();
	idx.IsCon = true;
	idx.ConVal = fe0;
	idx.IsBootstrap = false;
	calcs.push_back(idx);
	bases.push_back(pdata);
      };
    */
  }

  
  //
  // Setup forward/reverse time series calculations
  //
  {
    int nana = cli.ntimes;
    if ( cli.fwdrev )
      {
	for ( int i=0; i<nana; ++i )
	  {
	    double tstart = 0.;
	    double tstop = (i+1)/((double)(nana));
	    //std::printf("fwd %14.8f %14.8f\n",tstart,tstop);
	    edge_t e( cdata->ExtractRange(tstart,tstop) );

	    int baseidx = bases.size();
	    bases.push_back(e);
	    calc_t idx;
	    idx.BaseIdx = baseidx;
	    idx.IsFwd = true;
	    idx.RangeLo = tstart;
	    idx.RangeHi = tstop;
	    idx.IsBootstrap = false;
	    calcs.push_back(idx);
	  }
	for ( int i=0; i<nana-1; ++i )
	  {
	    double tstart = 1. - (i+1)/((double)(nana));
	    double tstop = 1.;
	    //std::printf("rev %14.8f %14.8f\n",tstart,tstop);
	    edge_t e( cdata->ExtractRange(tstart,tstop) );

	    int baseidx = bases.size();
	    bases.push_back(e);
	    calc_t idx;
	    idx.BaseIdx = baseidx;
	    idx.IsRev = true;
	    idx.RangeLo = tstart;
	    idx.RangeHi = tstop;
	    idx.IsBootstrap = false;
	    calcs.push_back(idx);
	  }
      }
    if ( cli.halves )
      {
	for ( int i=0; i<nana; ++i )
	  {
	    double tstart = 1. - (i+1)/((double)(nana));
	    double tmid =  0.5*(1. + tstart);
	    double tstop = 1.;
	    //std::printf("rev %14.8f %14.8f\n",tstart,tstop);
	    edge_t e( cdata->ExtractRange(tstart,tmid) );

	    int baseidx = bases.size();
	    bases.push_back(e);
	    calc_t idx;
	    idx.BaseIdx = baseidx;
	    idx.IsHalf1 = true;
	    idx.RangeLo = tstart;
	    idx.RangeHi = tstop;
	    idx.IsBootstrap = false;
	    calcs.push_back(idx);
	  }
	for ( int i=0; i<nana; ++i )
	  {
	    double tstart = 1. - (i+1)/((double)(nana));
	    double tmid =  0.5*(1. + tstart);
	    double tstop = 1.;
	    //std::printf("rev %14.8f %14.8f\n",tstart,tstop);
	    edge_t e( cdata->ExtractRange(tmid,tstop) );

	    int baseidx = bases.size();
	    bases.push_back(e);
	    calc_t idx;
	    idx.BaseIdx = baseidx;
	    idx.IsHalf2 = true;
	    idx.RangeLo = tstart;
	    idx.RangeHi = tstop;
	    idx.IsBootstrap = false;
	    calcs.push_back(idx);
	  }
      }
  }


  //
  // Setup bootstrap calculations
  //
  int const nbase = bases.size();
  for ( int ibase=0; ibase<nbase; ++ibase )
    {
      if ( ! calcs[ibase].IsCon )
	{
	  for ( int iboot=0; iboot < cli.nbootstrap; ++iboot )
	    {
	      calc_t idx( calcs[ibase] );
	      idx.IsBootstrap = true;
	      calcs.push_back( idx );
	    }
	};
    }


  //
  // Run all calculations in parallel
  //
  int const ncalc = calcs.size();
  std::vector<edgembar::StateValues> GlbValues( nbase );
#ifdef WITH_OPENMP
#pragma omp parallel
#endif
  {
    std::vector<edgembar::StateValues> myValues( nbase );
    std::vector<double> myp( nparam, 0. );
    std::vector<double> myguess( pguess );
    std::vector<double> myf( nparam, 0. );
#ifdef WITH_OPENMP
    int myverbosity = cli.verbosity;
#pragma omp for schedule(dynamic)
#endif
    for ( int icalc=0; icalc<ncalc; ++icalc )
      {
	calc_t idx( calcs[icalc] );

	//std::printf("icalc=%i\n",icalc);
	//std::printf("basidx=%i\n",idx.BaseIdx);
	//std::printf("nbases=%lu\n",bases.size());

	if ( cli.verbosity > 0 )
	  {
	    idx.Write( std::cout );
	    std::cout << "\n";
	  };
	
	edge_t e = bases[idx.BaseIdx];
	double tol = cli.tol;
	if ( idx.IsBootstrap )
	  {
	    e = e->ExtractBootstrap();
	    tol = cli.btol;
	  }
	// -----------------------------------------------------------
#ifdef WITH_OPENMP
	double chisq = e->Optimize( myguess, myp, myf, tol, myverbosity );
	if ( idx.IsBootstrap )
	  {
	    myValues[idx.BaseIdx].PushBootstrap(myf);
	  }
	else
	  {
	    myValues[idx.BaseIdx].SetResult(myf);
	    if ( idx.IsCon )
	      {
		calcs[icalc].Chisq = chisq;
	      }
	  };
#else
	// -----------------------------------------------------------
	double chisq = e->Optimize( pguess, popt, fopt, tol, cli.verbosity );
	if ( idx.IsBootstrap )
	  {
	    GlbValues[idx.BaseIdx].PushBootstrap(fopt);
	  }
	else
	  {
	    GlbValues[idx.BaseIdx].SetResult(fopt);
	    if ( idx.IsCon )
	      {
		calcs[icalc].Chisq = chisq;
	      }
	  };
#endif
	// -----------------------------------------------------------
      }
#ifdef WITH_OPENMP
#pragma omp critical
    {
      for ( int ibase=0; ibase<nbase; ++ibase )
	{
	  //std::cout << "Join " << ibase << " / " << nbase << std::endl;
	  GlbValues[ibase].Join( myValues[ibase] );
	}
    }
#endif
  }

#define SPC(n) std::setw((n)) << ""

  calcs.resize( nbase );

  //std::cout << "Main GetPython" << std::endl;
  
  std::string pystr = pdata->GetPython( argc, argv, beta, calcs, GlbValues, cli.ptol );
  
  //std::cout << "Got Python " << pystr << std::endl;
  
  std::string fname = strip_suffix(cli.inpfile,".xml") + ".py";
  if ( cli.outfile.size() > 0 )
    {
      fname = cli.outfile;
    };
  std::ofstream cout;
  //std::cout << "Open " << fname << "\n";
  cout.open(fname.c_str());

  //std::cout << "WriteOutput" << "\n";
  WriteOutput(cout,pystr);
  //std::cout << "Main done" << "\n";

  
  /*
  std::cout << "#!/usr/bin/env python3\n"
	    << "import edgembar\n"
	    << "import numpy as np\n\n"
	    << "edge = " << pystr << ";\n\n"
	    << "if __name__ == '__main__':\n"
	    << SPC(4) << "import edgembar\n"
	    << SPC(4) << "edgembar.ProcessEdgeCli(edge)\n";
  */

  /*
	    << SPC(4) << "import argparse\n"
	    << SPC(4) << "import sys\n"
	    << SPC(4) << "from pathlib import Path\n"
	    << SPC(4) << "parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=\"\"\"Output analysis of the edge in the desired format\"\"\")\n"
	    << SPC(4) << "parser.add_argument(\"--html\",help=\"Write HTML path/{edge}.html where path is the same path as this script and edge is replaced with the edge's name\",action='store_true',required=False)\n"
	    << SPC(4) << "parser.add_argument(\"--xml\",help=\"Write XML to stdout\",action='store_true',required=False)\n"
    	    << SPC(4) << "parser.add_argument(\"--imgfmt\",help=\"Image format. If left unset, then HTML is generated with the google-charts API. The resulting HTML file sizes are small, but they may take a few seconds to load. One can set this option to \\\"png\\\" or \\\"jpeg\\\" to embed the images into the HTML as a base64-encoded string. The resulting HTML files load very quickly, but the file sizes can be large.\",required=False)\n\n"
	    << SPC(4) << "args = parser.parse_args()\n\n"
	    << SPC(4) << "if args.xml:\n"
	    << SPC(8) << "print(edge)\n"
	    << SPC(4) << "else:\n"
	    << SPC(8) << "fname = Path(sys.argv[0]).parent / (edge.name + \".html\")\n"
	    << SPC(8) << "edge.WriteHtmlFile(fname,imgfmt=args.imgfmt)\n";
  */
  
  return 0;
}

