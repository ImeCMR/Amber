#include "cli_options.hpp"
#include "../edgembar/ReadMode.hpp"
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <sstream>
#include <string>
#include <algorithm>



edgembar::cli_options::cli_options()
  : temp(298.0),
    tol(1.e-13),
    btol(1.e-7),
    ptol(0.31729526862),
    nbootstrap(20),
    verbosity(0),
    nullguess(false),
    readmode(edgembar::AUTO),
    inpfile("inpfile"),
    outfile(""),
    ncon(2),
    dcon(2.),
    ntimes(4),
    fwdrev(false),
    halves(false),
    fstart(-1.),
    fstop(1.),
    stride(1),
    autoequil(true),
    uwts(false),
    autoeqmode(1)
{}


void edgembar::cli_options::print_help()
{
  //                    1         2         3         4         5         6         7         8
  //           12345678901234567890123456789012345678901234567890123456789012345678901234567890
  std::printf("Performs MBAR optimization of a \"graph edge\", generally defined as a\n");
  std::printf("difference in free energies in two environments.\n");
  std::printf("\n");
  std::printf("Usage: edgembar [-v <n>] [-t <f>] [--tol <f>] [--btol <f>] [--nboot <n>]\n");
  std::printf("                [--ptol <f>] [--ncon <n>] [--dcon <f>]\n");
  std::printf("                [--fstart <f>] [--fstop <f>] [--stride <i>]\n");
  std::printf("                [--times <n>] [--halves] [--fwdrev]\n");
  std::printf("                [--null] [--no-auto] [--auto-algo]\n");
  std::printf("                [--mode <s>] [--out <s>] inpfile\n");
  std::printf("\n");
  std::printf("Required arguments:\n");
  std::printf(" inpfile                input file\n");
  std::printf("\n");
  std::printf("Options:\n");
  std::printf(" --out,-o <s>           Output python file. If not specified, then '.py' is\n");
  std::printf("                        appended to the input filename\n");
  std::printf(" --temp,-t <f>          Temperature (K) [Default: 298]\n");
  std::printf(" --tol,-T <f>           Optimization relative tolerance [Default: 1.e-13]\n");
  std::printf(" --btol,-B <f>          Optimization relative tolerance during bootstrap\n");
  std::printf("                        calculations [Default: 1.e-7]\n");
  std::printf(" --ptol,-p <f>          The p-value (significance level) used in the Welch \n");
  std::printf("                        T-test for determining the production region of each\n");
  std::printf("                        simulation. Large values (p =~ 1.0) will cause more of\n");
  std::printf("                        the simulation to be discarded; small values (p =~ 0.0)\n");
  std::printf("                        will include more samples within the production region.\n");
  std::printf("                        [Default: 0.31729526862]. p-values less than ptol\n");
  std::printf("                        implies that the difference between two means is\n");
  std::printf("                        significant. If ptol <= 0, then T-tests are not\n");
  std::printf("                        performed; instead, it checks if the difference in\n");
  std::printf("                        means is less than 1 standard error of the propagated\n");
  std::printf("                        error. The default value is 1-0.68270473138, where\n");
  std::printf("                        0.6827... is a 68%% confidence interval (1 standard\n");
  std::printf("                        error), such that the T-test is very similar to\n");
  std::printf("                        accepting |mu1-mu2| < 1.0*sqrt(err1**2+err2**2).\n");
  std::printf(" --nboot,-b <n>         Number of bootstrap samples to estimate errors\n");
  std::printf("                        [Default: 20] (Recommended: 20 or higher)\n");
  std::printf(" --verbosity,-v <n>     Verbosity level. No printing is <= 0.\n");
  std::printf("                        [Default: 0] (0 is recommended if openmp parallel)\n");
  std::printf(" --null,-G              If present, then initial guess free energies are 0.\n");
  std::printf("                        The default is to use exponential averaging. Use\n");
  std::printf("                        the --null option if exp.avg. is numerically unstable\n");
  std::printf(" --ncon,-c <n>          Create a spline of the objective function using 2*n+1\n");
  std::printf("                        points centered about the unconstrained solution.\n");
  std::printf("                        [Default: 2, resulting in a 5-point spline]\n");
  std::printf("                        Setting this value to 0 will skip the spline generation\n");
  std::printf(" --dcon,-C <f>          The maximum displacement of the generated spline in\n");
  std::printf("                        kcal/mol. [Default: 2.0]\n");
  std::printf(" --ntimes,-A <n>        The number of time series calculations [Default: 4,\n");
  std::printf("                        which schedules: 0.25, 0.50, 0.75, 1.00\n");
  std::printf(" --halves,-H            First and last half analysis of the correlated samples\n");
  std::printf("                        after excluding someportion of equilibration\n");
  std::printf(" --fwdrev,-F            Forward and reverse analysis of the correlated samples\n");
  std::printf(" --fstart <f>           Skip the first <f> fraction of samples when reading\n");
  std::printf("                        the data files. [Default: 0.0]\n");
  std::printf(" --fstop <f>            Ignore the last <f> fraction of samples when reading\n");
  std::printf("                        the data files. [Default: 1.0]\n");
  std::printf(" --stride <i>           Stride through the data when reading the data files.\n");
  std::printf("                        [Default: 1]\n");
  std::printf(" --no-auto              The default behavior is to perform an automatic\n");
  std::printf("                        equilibration procedure that identifies the production\n");
  std::printf("                        region of each simulation and extract the statistically\n");
  std::printf("                        independent samples. If --no-auto is used, then all\n");
  std::printf("                        samples from --fstart to --fstop are analyzed and the\n");
  std::printf("                        errors are estimated from block bootstrap analysis,\n");
  std::printf("                        where the block size is chosen from the autocorrelation\n");
  std::printf("                        in the data.\n");
  std::printf(" --auto-algo <i>        Select the algorithm for detecting the equilibrated and\n");
  std::printf("                        production regions of each simulation\n");
  std::printf("                        At the present time, there is only 1 algorithm\n");
  std::printf("                        corresponding to --auto-algo=1. This option exists for\n");
  std::printf("                        future development.\n");
  std::printf(" --mode <s>             Defines the objective function and the expected\n");
  std::printf("                        filenames.\n");
  std::printf("                        --mode=AUTO (default) :\n");
  std::printf("                                  Automatically detect based on the existing\n");
  std::printf("                                  filenames.\n");
  std::printf("                        --mode=MBAR :\n");
  std::printf("                                  A full matrix of trajectories and energies\n");
  std::printf("                        --mode=MBAREXP0 :\n");
  std::printf("                                  MBAR between windows [2,N] and exponential\n");
  std::printf("                                  averaging from 2->1\n");
  std::printf("                        --mode=MBAREXP1 :\n");
  std::printf("                                  MBAR between windows [1,N-1] and exponential\n");
  std::printf("                                  averaging from N-1->N\n");
  std::printf("                        --mode=MBAREXP\n");
  std::printf("                                   MBAR between windows [2,N-1] and\n");
  std::printf("                                   exponential averaging from N-1->N and 2->1\n");
  std::printf("                        --mode=BAR\n");
  std::printf("                                   N-1 pairs of MBAR objectives connecting\n");
  std::printf("                                   adjacent states [i,i+1]\n");
  std::printf("                        --mode=BAREXP0\n");
  std::printf("                                   BAR between all adjacent pairs of states\n");
  std::printf("                                   in the range [2,N], and exponential\n");
  std::printf("                                   averaging from 2->1\n");
  std::printf("                        --mode=BAREXP1\n");
  std::printf("                                   BAR between all adjacent pairs of states\n");
  std::printf("                                   in the range [1,N-1], and exponential\n");
  std::printf("                                   averaging from N-1->N\n");
  std::printf("                        --mode=BAREXP\n");
  std::printf("                                   BAR between windows [2,N-1]\n");
  std::printf("                                   and exponential averaging from N-1->N\n");
  std::printf("                                   and exponential averaging from 2->1\n");
  std::printf(" --uwts                 Use uniform weights in the MBAR and BAR objective\n");
  std::printf("                        functions rather than scaling them to account for\n");
  std::printf("                        the number of states. This is a debug option to make\n");
  std::printf("                        comparison with the graphmbar program; you should not\n");
  std::printf("                        use this option.\n");

  std::printf("\n");

}


edgembar::cli_options edgembar::read_options( int argc, char * argv[] )
{
  // Tell getopt to not print
  //extern int opterr;
  //opterr = 0;

  edgembar::cli_options cli;


  
  static struct option long_options[] =
    {
     { "help",      no_argument,       NULL, 'h'   },
     { "out",       required_argument, NULL, 'o'   },
     { "temp",      required_argument, NULL, 't'   },
     { "tol",       required_argument, NULL, 'T'   },
     { "btol",      required_argument, NULL, 'B'   },
     { "ptol",      required_argument, NULL, 'p'   },
     { "nboot",     required_argument, NULL, 'b'   },
     { "mode",      required_argument, NULL, 'm'   },
     { "verbosity", required_argument, NULL, 'v'   },
     { "null",      no_argument,       NULL, 'G'   },
     { "ncon",      required_argument, NULL, 'c'   },
     { "dcon",      required_argument, NULL, 'C'   },
     { "ntimes",    required_argument, NULL, 'A'   },
     { "halves",    no_argument,       NULL, 'H'   },
     { "fwdrev",    no_argument,       NULL, 'F'   },
     { "no-auto",   no_argument,       NULL, 'Q'   },
     { "fstart",    required_argument, NULL, 's'   },
     { "fstop",     required_argument, NULL, 'S'   },
     { "stride",    required_argument, NULL, 'g'   },
     { "uwts",      no_argument,       NULL, 'U'   },
     { "auto-algo", required_argument, NULL, 'E'   },
    {NULL,0,NULL,0}
    };


  cli.readmode = edgembar::AUTO;
  
  int opt = 0;
  int long_index = 0;
  //char * subopts, * value;
  while ( (opt = getopt_long
          ( argc, argv, "ho:t:T:B:b:m:v:Gc:C:A:HFQs:S:g:Up:E:", 
            long_options, &long_index )) != -1 )
    {
      switch ( opt )
	{
	  
	case 'o':     { cli.outfile  = optarg; break; }
	case 't':     { cli.temp  = std::atof(optarg); break; }
 	case 'T':     { cli.tol   = std::atof(optarg); break; }
 	case 'B':     { cli.btol  = std::atof(optarg); break; }
 	case 'p':     { cli.ptol  = std::atof(optarg); break; }
 	case 'b':     { cli.nbootstrap = std::atoi(optarg); break; }
	case 'm':     { cli.readmode = edgembar::GetMode(optarg); break; }
 	case 'v':     { cli.verbosity = std::max(0,std::atoi(optarg)); break; }
	case 'G':     { cli.nullguess = true; break; }
 	case 'c':     { cli.ncon = std::max(0,std::atoi(optarg)); break; }
 	case 'C':     { cli.dcon = std::atof(optarg); break; }
 	case 'A':     { cli.ntimes = std::max(0,std::atoi(optarg)); break; }
	case 'H':     { cli.halves = true; break; }
	case 'F':     { cli.fwdrev = true; break; }
	case 'Q':     { cli.autoequil = false; break; }
 	case 's':     { cli.fstart = std::max(0.,std::min(1.,std::atof(optarg))); break; }
 	case 'S':     { cli.fstop = std::max(0.,std::min(1.,std::atof(optarg))); break; }
 	case 'g':     { cli.stride = std::max(1,std::atoi(optarg)); break; }
	case 'U':     { cli.uwts = true; break; }
	case 'E':     { cli.autoeqmode = std::atoi(optarg); break; }
     
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
      std::printf("Error: Missing nonoptional argument: inpfile\n");
      std::exit(EXIT_FAILURE);
    }
  else if ( optind == argc-1 )
    {
      cli.inpfile = argv[optind];
    }
  else
    {
      std::printf("Error: Expected only one non-optional argument\n");
      std::exit(EXIT_FAILURE);
    }


  if ( cli.autoeqmode < 1 or cli.autoeqmode > 1 )
    {
      std::printf("Error: --auto-algo must be 1\n");
      std::exit(EXIT_FAILURE);
    }
  
  if ( cli.ntimes < 1 )
    {
      cli.halves = false;
      cli.fwdrev = false;
    }

  if ( cli.fstart >= cli.fstop )
    {
      std::printf("Error: --fstart must be less than --fstop\n");
      std::exit(EXIT_FAILURE);
    }

  //cli.ptol = std::max(0.,std::min(1.,cli.ptol));

  if ( cli.readmode < edgembar::AUTO or cli.readmode > edgembar::BAREXP )
    {
      std::printf("Error: --mode value out of bounds %i\n",cli.readmode);
      std::exit(EXIT_FAILURE);
    }

  return cli;
}
