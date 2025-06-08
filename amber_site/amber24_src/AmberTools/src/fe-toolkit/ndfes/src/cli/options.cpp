
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "options.hpp"

ndfes::options::options()
  : temp(298.0),
    nham(0),
    maxiter(99999),
    nquad(5),
    order(5),
    nboot(0),
    reltol(1.e-9),
    btol(1.e-5),
    vfep(false),
    mbar(false),
    sdosFullSampling(true),
    sdosHistPrefix(""),
    sdos(-1.),
    startfrac(0.),
    stopfrac(1.),
    uniq(false),
    maxuniq(0),
    expert(false)
{}

void ndfes::options::print_help()
{
  //                    1         2         3         4         5         6         7         8
  //           12345678901234567890123456789012345678901234567890123456789012345678901234567890
  std::printf("Calculate a N-dimensional potential of mean force profile using either MBAR or\n");
  std::printf("variational free energy profile (vFEP) within a cardinal B-spline basis.\n");
  std::printf("\n\n");
  std::printf("MBAR Usage: ndfes --mbar -c <chk> -w <dx> [-w <dy> [...]] \n");
  std::printf("                  [--sdos <nrand> [--sepdos] [--hist <prefix>]]\n");
  std::printf("                  [--nham <nham>] [-T <tol>] [-n <maxit>]\n");
  std::printf("                  [--nboot <nboot>] [--btol <btol>] [-u]\n");
  std::printf("                  metafile \n");
  std::printf("\n\n");
  std::printf("vFEP Usage: ndfes --vfep -c <chk> [-q <nq>] [-o <order>] -w <dx> [-w <dy> [...]]\n");
  std::printf("                  [-t <temp>] [-g <nrand>] [-T <tol>] [-n <maxit>] [-u]\n");
  std::printf("                  metafile\n");
  std::printf("\n\n");
 
  std::printf("\n");
  std::printf("Required arguments:\n");
  std::printf(" metafile               File containing umbrella window information\n");
  std::printf("                        More than one metafile can be specified\n");
  std::printf(" --chkpt,-c <chk>       This is the output of the program. It is read by\n");
  std::printf("                        the ndfes python library to generate images and\n");
  std::printf("                        search for minimum free energy paths\n");
  std::printf(" --width,-w <dx>        Bin width size. Use this option multiple times to\n");
  std::printf("                        define the bin width for each dimension.\n");
  std::printf("                        If only one instance is used, then all dimensions\n");
  std::printf("                        share the same bin width.\n");
  std::printf("\n");
  std::printf("vFEP-Specific Options:\n");
  std::printf(" --nquad,-q <nq>        Number of Gauss-Legendre quadrature points used to\n");
  std::printf("                        integrate the partition function within each bin\n");
  std::printf("                        nq is the number of points per dimension. Default: 5\n");
  std::printf(" --order,-o <order>     The B-spline order. A value of 3 is required to yield\n");
  std::printf("                        a smooth FES. Increasing the order may produce\n");
  std::printf("                        polynomial artifacts from incomplete distributions.\n");
  std::printf("                        Default 5.\n");
  std::printf("\n");
  std::printf("MBAR-Specific Options:\n");
  std::printf(" --nham, -G <nham>      Number of potential energies to expect in each tracefile\n");
  std::printf(" --sdos, -D <dU>        If dU > 0, smooth the MBAR density-of-states using the\n");
  std::printf("                        method in: https://doi.org/10.1021/acs.jctc.0c00794\n");
  std::printf("                        The value of dU is in reduced energy units. The value\n");
  std::printf("                        used in the paper was 0.2. Default is -1. (no smoothing)\n");
  std::printf(" --sepsdos, -P          If present, then perform density-of-states smoothing from\n");
  std::printf("                        from the samples of the individual reference potentials\n");
  std::printf("                        rather than the combined sampling from all potentials.\n");
  std::printf("                        It is arguably more correct to not use this option, but\n");
  std::printf("                        one may find it to yield more stable results in some\n");
  std::printf("                        circumstances.\n");
  std::printf(" --hist, -O <prefix>    Write the density-of-states energy histograms to file. The\n");
  std::printf("                        filenames will be:\n");
  std::printf("                            <prefix>hist_bin_<R>_target_<T>_ref_<H>.dat\n");
  std::printf("                        where <R> is the spatial bin index, <T> is the target\n");
  std::printf("                        potential index, and <H> is the reference potential index.\n");
  std::printf("\n");
  std::printf("Options Applicable to both methods:\n");
  std::printf("--temp,-t <temp>        Temperature. Default 298\n");
  std::printf("--tol,-T <tol>          Optimization stopping criteria. Default 1.e-9\n");
  std::printf("--maxit,-n <maxit>      Perform a maximum of maxit iterations. Default 99999\n");
  std::printf("--periodic,-p <dim>     Flags dimension <dim> as a periodic coordinate\n");
  std::printf("                        ranging from 0 to 360. The option can be used\n");
  std::printf("                        multiple times for multiple periodic coordinates.\n");
  std::printf("--start,-s <frac>       Discard the first <frac> fraction of data as\n");
  std::printf("                        equilibration. Range [0.0,1.0). Default: 0.0\n");
  std::printf("--stop,-S <frac>        Discard samples after <frac> fraction of data.\n");
  std::printf("                        Range (0.0,1.0]. Default: 1.0\n");
  std::printf("--uniq,-u               If present, then prune the inputs to keep only\n");
  std::printf("                        the statistically unique samples.  The default\n");
  std::printf("                        is to keep all samples and estimate errors from\n");
  std::printf("                        block bootstrap analysis, where the block size\n");
  std::printf("                        is chosen from the statistical inefficiency of the\n");
  std::printf("                        biasing energy.\n");
  std::printf("--maxuniq,-U <int>      Maximum number of statistically indepedent samples\n");
  std::printf("                        to use from each simulation. To use this option,\n");
  std::printf("                        you must also use --uniq.  The default value is 0,\n");
  std::printf("                        which analyzes all statistically unique samples.\n");
  std::printf("                        This option is only really useful when trying to\n");
  std::printf("                        determine the smallest subset of data that reproduces\n");
  std::printf("                        the results from the full set of sampling.\n");
  std::printf("--expert,-X             Don't terminate abnormally if the target and simulation\n");
  std::printf("                        temperatures disagree. It is recommended that all of\n");
  std::printf("                        your simulations are run at the same temperature and\n");
  std::printf("                        you analyze the FES with that same temperature. The\n");
  std::printf("                        support for mixing different temperatures is still\n");
  std::printf("                        experimental and likely won't give you good results.\n");
  std::printf(" --nboot, -M <nboot>    Number of bootstrap resamples to estimate errors. Default: 0\n");
  std::printf(" --btol, -B <btol>      Optimization tolerance for bootstrap solutions. When MBAR\n");
  std::printf("                        is used, a value <= 0 will estimate the solution of the\n");
  std::printf("                        MBAR/UWHAM equations by performing 1 iteration of the\n");
  std::printf("                        MBAR self-consistent equations, rather than a full nonlinear\n");
  std::printf("                        optimization. Default: 1.e-5\n");
  std::printf("\n");
  std::printf("Suggested vFEP options:\n");
  std::printf("    ndfes --vfep -w 0.15 -o 5 -q 5 -c vfep.15o5q5.chk\n");
  std::printf("Alternatively, one can try several sets of similar control parameters to estimate\n");
  std::printf("the uncertainty of the analysis. Here is an example set of calculations I've\n");
  std::printf("used for 2d FESs\n");
  std::printf("    ndfes --vfep -w 0.12 -o 4 -q 5 -c vfep.12o4q5.chk\n");
  std::printf("    ndfes --vfep -w 0.13 -o 4 -q 6 -c vfep.13o4q6.chk\n");
  std::printf("    ndfes --vfep -w 0.14 -o 4 -q 6 -c vfep.14o4q6.chk\n");
  std::printf("    ndfes --vfep -w 0.15 -o 4 -q 7 -c vfep.15o4q7.chk\n\n");
  std::printf("Suggested MBAR options:\n");
  std::printf("    ndfes --mbar -w 0.15 -c mbar.15.chk\n");
  std::printf("The recommended width for periodic calculations is 10 or 12 degrees.\n");
  std::printf("The recommended width for distances is 0.15 Angstroms.\n");

}


ndfes::options ndfes::ReadOptions( int argc, char * argv[] )
{
  // Tell getopt to not print
  //extern int opterr;
  //opterr = 0;

  ndfes::options cli;
  
  static struct option long_options[] =
    {
     { "help",    no_argument,       NULL, 'h'   },
     { "vfep",    no_argument,       NULL, 'z'   },
     { "mbar",    no_argument,       NULL, 'H'   },
     { "nham",    required_argument, NULL, 'G'   },
     { "chkpt",   required_argument, NULL, 'c'   },
     { "order",   required_argument, NULL, 'o'   },
     { "temp",    required_argument, NULL, 't'   },
     { "maxit",   required_argument, NULL, 'n'   },
     { "tol",     required_argument, NULL, 'T'   },
     { "nquad",   required_argument, NULL, 'q'   },
     { "width",   required_argument, NULL, 'w'   },
     { "periodic",required_argument, NULL, 'p'   },
     { "sdos",    required_argument, NULL, 'D'   },
     { "hist",    required_argument, NULL, 'O'   },
     { "sepsdos", no_argument,       NULL, 'P'   },
     { "start",   required_argument, NULL, 's'   },
     { "stop",    required_argument, NULL, 'S'   },
     { "nboot",   required_argument, NULL, 'M'   },
     { "btol",    required_argument, NULL, 'B'   },
     { "uniq",    no_argument,       NULL, 'u'   },
     { "maxuniq", required_argument, NULL, 'U'   },
     { "expert",  no_argument,       NULL, 'X'   },
   {NULL,0,NULL,0}
    };


  cli.target_widths.resize(0);
  int opt = 0;
  int long_index = 0;
  //char * subopts, * value;
  while ( (opt = getopt_long
          ( argc, argv, "hzPXHuG:c:o:t:n:T:q:w:p:D:s:S:M:B:O:U:", 
            long_options, &long_index )) != -1 )
    {

      switch (opt)
        {
	case 'X':  { cli.expert = true; break; }
	case 'z':  { cli.vfep = true; cli.mbar = false; break; }
        case 'H':  { cli.mbar = true; cli.vfep = false; break; }
	case 'P':  { cli.sdosFullSampling = false; break; }
	case 'u':  { cli.uniq = true; break; }
	case 'U':  { cli.maxuniq   = std::atoi(optarg); break; }
	case 'O':  { cli.sdosHistPrefix  = optarg; break; }
	case 'G':  { cli.nham = std::atoi(optarg); break; }
	case 'c':  { cli.chkpt  = optarg; break; }
	case 'o':  { cli.order = std::atoi(optarg); break; }
	case 't':  { cli.temp      = std::atof(optarg); break; }
	case 'n':  { cli.maxiter   = std::atoi(optarg); break; }
	case 'T':  { cli.reltol    = std::atof(optarg); break; }
	case 'q':  { cli.nquad = std::atoi(optarg); break; }
	case 'w':  { cli.target_widths.push_back( std::atof(optarg) ); break; }
	case 'p':  { cli.perdims.push_back( std::atoi(optarg) ); break; }
	case 'D':  { cli.sdos = std::atof(optarg); break; }
	case 's':  { cli.startfrac = std::atof(optarg); break; }
	case 'S':  { cli.stopfrac = std::atof(optarg); break; }
	case 'M':  { cli.nboot = std::atoi(optarg); break; }
	case 'B':  { cli.btol = std::atof(optarg); break; }
        case '?': // getopt_long already printed an error
	  std::printf("%s: use -h for usage\n", argv[0]);
          std::exit(EXIT_FAILURE);
          break;
        case 'h':
          cli.print_help();
          std::exit(EXIT_SUCCESS);
          break;
        default:
          std::printf("An error occured while reading command line options\n");
          std::exit(EXIT_FAILURE);
        }
    };

  if ( optind == argc )
    {
      std::printf("Error: Missing nonoptional argument(s): metafile\n");
      std::exit(EXIT_FAILURE);
    }
  else
    {
      for ( int iarg=0; optind<argc; ++optind, ++iarg )
	{
	  if ( iarg > 0 )
	    {
	      std::printf("Error: Unexpected command line argument: %s\n",argv[optind]);
	      std::exit(EXIT_FAILURE);
	    }
	  cli.metafile = argv[optind];
	}
    }


  if ( cli.maxuniq != 0 and not cli.uniq )
    {
      std::printf("Error: One must use --uniq if --maxuniq is specified\n");
      std::exit(EXIT_FAILURE);
    }
  
  if ( cli.nboot < 0 )
    {
      cli.nboot = 0;
    }

  if ( cli.btol < 1.e-15 )
    {
      cli.btol = -1.;
    }


  if ( cli.startfrac < 0. )
    {
      cli.startfrac = 0.;
    }
  else if ( cli.startfrac >= 1. )
    {
      std::printf("Error: --start value must be less than 1.0\n");
      std::exit(EXIT_FAILURE);
    };

  if ( cli.stopfrac > 1. )
    {
      cli.stopfrac = 1.;
    }
  else if ( cli.stopfrac <= 0. )
    {
      std::printf("Error: --stop value must be greater than 0.0\n");
      std::exit(EXIT_FAILURE);
    };

  if ( cli.stopfrac <= cli.startfrac )
    {
      std::printf("Error: --stop value must be greater than --start value\n");
      std::exit(EXIT_FAILURE);
    }

  
  
  if ( cli.chkpt.size() == 0 )
    {
      std::printf("Error: --chkpt is mandatory\n");
      std::exit(EXIT_FAILURE);
    }
  
  if ( cli.vfep and cli.order < 1 )
    {
      std::printf("Error: --order must be 1 or larger if --vfep is used\n");
      std::exit(EXIT_FAILURE);
    }

  if ( cli.mbar )
    {
      cli.order = 1;
    }
  
  if ( cli.nham < 0 )
    {
      cli.nham = 0;
    }
  /*
  else if ( cli.nham == 1 )
    {
      std::printf("Error: --nham must be =0 or >1\n");
      std::exit(EXIT_FAILURE);
    };
  */
  
  if ( cli.target_widths.size() < 1 )
    {
      cli.target_widths.push_back(0.15);
      std::printf("WARNING: bin widths were not specified with -w. Defaulting to 0.15 for each dimension. This is likely fine for distances in Angstroms, but too small for angles measured in degrees.\n");
    }

  if ( cli.sdos > 0. and cli.mbar )
    {
      std::printf("MBAR density-of-states smoothing will be applied with energy histogram bin widths %0.2f (reduced energy units)\n",cli.sdos);
    }

  if ( not cli.vfep and not cli.mbar )
    {
      std::printf("Error: You must specify --vfep or --mbar. See ndfes --help\n");
      std::exit(EXIT_FAILURE);
    }
  
  return cli;
}
