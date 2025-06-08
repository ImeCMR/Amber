#include "Options.hpp"

#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>


Options::Options()
  : disang(""),
    mdin(""),
    ipath(""),
    opath(""),
    chk(""),
    curit(-1),
    odir(""),
    //dacc(false),
    model(0),
    rbf(false),
    bspl(false),
    hist(false),
    arbf(-1),
    shape(100.),
    minsize(0),
    rbfmaxerr(1.0),
    maxit(-1),
    ptol(0.95),
    distol(0.05),
    angtol(0.5),
    mlen(5),
    //closest(false),
    temp(298.),
    fix0(false),
    fix1(false),
    //dampfactor(1.),
    //smooth_frac(1.),
    //smooth_nmin(3),
    //smooth_nmax(11),
    smoothfc(false),
    //scalefc(false),
    predfc(false),
    maxfc(300.),
    minfc(50.),
    maxfc_ang(0.5),
    minfc_ang(0.1),
    tdispfc(0.75),
    //pkl(false),
    //nquad(5),
    cp_nearest(false),
    boot(false),
    booterror(-1),
    npathpts(-1),
    nsplpts(-1),
    //nebk(-1),
    //conv_layers(-1),
    //conv_samples(50),
    //tetris(false),
    //tetris_smooth(false),
    //acc(false),
    acc_oobk(300.),
    acc_maxit(100),
    buffer(false),
    //sim_layers(-1),
    nlayers(-1),
    dry_run(false),
    //explore_once(false),
    neqit(-1),
    //disp_limit(false),
    wavg(0),
    wavg_niter(0),
    prefix(""),
    smooth(true),
    msm(false),
    gbsm(false),
    stype(-1),
    printfreq(50),
    pad(3),
    sasm(false) //,
    //gxsm(false)
{}

void Options::print_help() const
{
  //                    1         2         3         4         5         6         7         8
  //           12345678901234567890123456789012345678901234567890123456789012345678901234567890
  std::printf("usage: ndfes-path [--help]\n");
  std::printf("\n(options for creating a new directory of simulations)\n");
  std::printf("       [-d DISANG -i MDIN --curit CURIT] [--gbsm | --msm]\n");
  std::printf("       [--odir ODIR] [--dry-run] [--prefix PREFIX] [--extra-prefix XPREF]\n");
  std::printf("       [--cp-nearest] [--neqit N]\n");
  std::printf("\n(options for path projections and optimization)\n");
  std::printf("       [--ipath IPATH --chk CHK --opath OPATH [--nsplpts NSP]]\n");
  std::printf("\n(options that control the FES representation)\n");
  std::printf("       [-m MODEL] [--minsize MINSIZE] [--bspl] [--hist] \n");
  std::printf("       [--wavg N] [--wavg-niter N]\n");
  std::printf("       [--rbf] [--arbf ARBF] [-R SHAPE] [--rbfmaxerr RBFMAXERR]\n");
  std::printf("       [--buffer] [--boot] [--boot-error ERROR]\n");
  std::printf("\n(options that control the movement of the path)\n");
  std::printf("       [--maxit MAXIT] [--ptol PTOL] \n");
  std::printf("       [--distol DISTOL] [--angtol ANGTOL] [--mlen MLEN]\n");
  std::printf("       [--acc-oobk K] [--acc-maxit]\n");
  std::printf("\n(options that control the path representation)\n");
  std::printf("       [--npathpts NPP] [--linear] [--akima] [--no-smooth]\n");
  std::printf("       [--nlayers NS] [--fix0] [--fix1]\n");
  std::printf("       [--predfc [--smoothfc] [--tdispfc TDISP]]\n");
  std::printf("       [--maxfc MAXFC] [--minfc MINFC]\n");
  std::printf("       [--maxfc-ang MAXFC_ANG] [--minfc-ang MINFC_ANG]\n");
  std::printf("\n(misc options)\n");
  std::printf("       [-T TEMP] [--pad PAD] [--print-freq FREQ]\n");
  std::printf("\n");
  std::printf("\n");
  std::printf("--------------------------------------------------------------------------------\n");
  std::printf("**Required options for generating the next directory of simulations using\n");
  std::printf("the modified string method described in Rosta et. al., J. Am. Chem.\n");
  std::printf("Soc. (2011) 133, 8934:\n\n");
  std::printf("ndfes-path -d DISANG -i MDIN --curit CURIT --maxit=0 --msm --no-smooth\n\n");
  std::printf("You can also choose a spline representation, either: --akima or --linear\n");
  std::printf("By setting --maxit=0, the new spline will be constructed from the\n");
  std::printf("observed reaction coordinate averages obtained from simulation (by\n");
  std::printf("averaging the appropriate columns in the amber dumpave output file).\n");
  std::printf("The discretization of the path must be the same as the number of simulations\n");
  std::printf("\n");
  std::printf("**Required options for generating the next directory of simulations using\n");
  std::printf("the surface-accelerated string method:\n\n");
  std::printf("ndfes-path -d DISANG -i MDIN --curit CURIT --sasm --minsize=10\n\n");
  std::printf("By default --maxit=300 --npathpts=100 when setting --sasm; the new spline will\n");
  std::printf("be constructed by performing minimizations on the biased free energy surface\n");
  std::printf("in the reduced space of reaction coordinates. In this manner, the output path\n");
  std::printf("will be the minimum free energy pathway along the current estimate of the \n");
  std::printf("unbiased free energy surface. The value of npathpts sets the number of\n");
  std::printf("synthetic images.\n");
  std::printf("\n");
  std::printf("\n");

  std::printf("--------------------------------------------------------------------------------\n");
  std::printf("Options for printing a 1D projection through a FES\n\n");
  std::printf("ndfes-path --ipath IPATH --chk CHK --opath OPATH --maxit=0\n");
  std::printf("                       [--nsplpts NSP]\n");
  std::printf("\n");
  std::printf("Options for optimizing a path on a fixed FES without generating a new\n");
  std::printf("directory of simulations\n\n");
  std::printf("ndfes-path --ipath IPATH --chk CHK --opath OPATH --maxit MAXIT\n");
  std::printf("                       [--npathpts NPTS] [--nsplpts NSP]\n");
  std::printf("\n\n");
  std::printf("--------------------------------------------------------------------------------\n");
  std::printf("Recommended options for representing the free energy surface:\n");
  std::printf("To smooth out a noisey FES: --wavg=4 --wavg-niter=1 --minsize=10\n");
  std::printf("To better reproduce the noise: --wavg=4 --wavg-niter=5 --minsize=10\n");
  std::printf("To exactly reproduce the noise: --rbf --minsize=10\n");
  std::printf("--------------------------------------------------------------------------------\n");


  std::printf("\n");
  std::printf("optional arguments:\n");
  std::printf("  -h, --help            show this help message and exit\n");
  std::printf("\n");
  std::printf("\nDirectory Creation\n");
  std::printf("  -d DISANG, --disang DISANG\n");
  std::printf("                        template disang file\n");
  std::printf("  -i MDIN, --mdin MDIN  template mdin file\n");
  std::printf("  --curit CURIT         Index of the current iteration. The current directory\n");
  std::printf("                        is assumed to be itXXX, where XXX is a zero-padded\n");
  std::printf("                        integer, specified by --curit. If --curit=0, then the\n");
  std::printf("                        directory name is assumed to be init\n");
  std::printf("  --odir ODIR           Output directory. If unspecified, it is itXXX, where\n");
  std::printf("                        XXX is XXX + 1\n");
  std::printf("  --pad PAD             The with of the file and directory name padding.\n");
  std::printf("                        Default: 3\n");
  std::printf("  --dry-run             Do not write any output files to disk\n");
  
  std::printf("  --neqit N             Remove the first N iterations as equilibration.\n");
  std::printf("                        Default: -1. The value of N is the last iteration\n");
  std::printf("                        to be excluded, such that N=0 would only exclude\n");
  std::printf("                        the simulation data in the init directory. This\n");
  std::printf("                        option only effects the simulations restart files\n");
  std::printf("                        that can be copied to the next directory. You will\n");
  std::printf("                        need to set the --neqit option in \n");
  std::printf("                        ndfes-FTSM-PrepareAnalysis.py to effect the FES.\n");

  std::printf("  --cp-nearest          If present, then the restart of the next iteration\n");
  std::printf("                        image is copied from the image whose centroid was\n");
  std::printf("                        closest to the new window center. The default behavior\n");
  std::printf("                        is to copy the restarts from the corresponding image\n");

  std::printf("  --prefix SUBDIR       The subdirectory to read and write mdin and rst7 files.\n");
  std::printf("                        By default, these files are written to itXX/, but\n");
  std::printf("                        setting a prefix will write them to itXX/SUBDIR/. This\n");
  std::printf("                        is useful when optimizing several paths on a common FES\n");
  std::printf(" --extra-prefix SUBDIR  Extra subdirectories to search for rst7 files.\n");
  std::printf("                        This is only useful when tryin to optimize several paths\n");
  std::printf("                        at the same time. By using this option, one can copy the\n");
  std::printf("                        restart file from any of the paths.  This option can be\n");
  std::printf("                        used more than once.\n");

  std::printf("\nPath Analysis and Optimization Without Directory Creation\n");

  std::printf("  --chk CHK             Read free energy surface from the specified ndfes\n");
  std::printf("                        checkpoint file. This option must be used with --ipath\n");
  std::printf("                        and --opath, and it cannot be used with \n");
  std::printf("                        --curit --mdin and --disang\n");
  std::printf("  --ipath IPATH         Read the umbrella centers and force constants from the\n");
  std::printf("                        specified metafile. META can also be a table whose\n");
  std::printf("                        first 2 columns are ignored. Columns 3-to-3+Ndim\n");
  std::printf("                        are assumed to be umbrella centers. All remaining\n");
  std::printf("                        columns are ignored. In this case, all force constants\n");
  std::printf("                        are assumed to be 200 kcal/mol/Ang**2. META can also\n");
  std::printf("                        be a ndfes path pkl file, from which it will read the\n");
  std::printf("                        last path.\n");
  std::printf("                        This option must be used with --chk and --opath, and\n");
  std::printf("                        it cannot be used with --curit --mdin and --disang\n");
  std::printf("  --opath OPATH         Write the output path to OPATH. This file will contain\n");
  std::printf("                        2+2*Ndim columns: i, t[i], rc[i][dim], and fc[i][dim]\n");
  std::printf("                        Two additional files will be written:\n");
  std::printf("                          OPATH.xml : the path iterations\n");
  std::printf("                          OPATH.interp.model.dat : A file containing the path\n");
  std::printf("                        properties: i t[i] rc[i][dim] F[i] dF[i] RE[i] N[i]\n");
  std::printf("                        where t is the progress variable, rc is the reaction\n");
  std::printf("                        coordinate, F is the free energy, Err is the standard\n");
  std::printf("                        error in F, RE is the reweighting entropy, and N is\n");
  std::printf("                        the number of samples in the bin that this point lies\n");
  std::printf("                        within.\n");
  std::printf("                        This option must be used with --ipath and --opath, and\n");
  std::printf("                        it cannot be used with --curit --mdin and --disang\n");
  std::printf("\nThe Free Energy Surface\n");
  std::printf("  -m MODEL, --model MODEL\n");
  std::printf("                        model to read from the checkpoint file (default: 0)\n");
  std::printf("  --minsize MINSIZE\n");
  std::printf("                        If a FES bin has fewer than MINSIZE samples, then\n");
  std::printf("                        discard it. Default: 0 (although I recommend\n");
  std::printf("                        --minsize=10)\n");
  std::printf("\n*FES Interpolator: The B-spline Weighted Average (wavg)\n");
  std::printf("  --wavg N              Approximate the FES by a weighted average of nearby bin\n");
  std::printf("                        values using B-splines. N is the B-spline order.\n");
  std::printf("                        Default: 0\n");
  std::printf("  --wavg-niter N        Iterative apply corrections to the wavg histogram values\n");
  std::printf("                        such that the interpolated value match the values at the\n");
  std::printf("                        bin centers. This will potentially correct artifacts\n");
  std::printf("                        introduced from the averaging procedure, but it will\n");
  std::printf("                        also reproduce noise in the data. Default: 0.\n");
  std::printf("\n*FES Interpolator: Radial Basis Function (rbf)\n");
  std::printf("  -r, --rbf             Radial Basis Function interpolation. Valid for both\n");
  std::printf("                        vFEP and MBAR\n");
  std::printf("  -R SHAPE, --shape SHAPE\n");
  std::printf("                        RBF shape parameter used when --rbf is specified.\n");
  std::printf("                        small values tend to yield larger oscillations.\n");
  std::printf("                        Default=100\n");
  std::printf("  --rbfmaxerr RBFMAXERR\n");
  std::printf("                        If a bin was excluded from the RBF solution, and the\n");
  std::printf("                        interpolated free energy differs from the histogram\n");
  std::printf("                        value by more than rbfmaxerr, then eliminate the bin\n");
  std::printf("                        entirely -- effectively treating the free energy as\n");
  std::printf("                        positive infinity. Default: 1.0 kcal/mol\n");
  std::printf("\n*FES Interpolator: Approximate Radial Basis Function (arbf)\n");
  std::printf("  --arbf ARBF           Approximate radial basis function parametrized from\n");
  std::printf("                        the central bin and its nearest ARBF neighbors.\n");
  std::printf("                        Default: -1, which disables the arbf. A value of 1\n");
  std::printf("                        constructs a RBF from the central bin and first layer\n");
  std::printf("                        of nearby neighbors. This would be fast, but noticable\n");
  std::printf("                        artifacts would be present in the FES. A value of 2\n");
  std::printf("                        reduces the artifacts. Setting ARBF >= 3 is not\n");
  std::printf("                        reccommended because it would likely be more expensive\n");
  std::printf("                        than simply using --rbf.\n");
  std::printf("\n*FES Interpolator: vFEP B-splines (bspl)\n");
  std::printf("  -b, --bspl            Cardinal B-spline interpolation. Valid only for vFEP\n");
  std::printf("\n*FES Interpolator: Histogram (hist)\n");
  std::printf("  -H, --hist            No interpolation; only use histogram values\n");
  std::printf("\n*Other FES options:\n");
  std::printf("  --buffer              Adds an extra layers of fake bins to the FES that\n");
  std::printf("                        increase the energy by 2 kcal/mol. This is recommended\n");
  std::printf("                        when using --acc to perform biased minimizations.\n");
  std::printf("                        The penalty buffer will prevent the string optimization\n");
  std::printf("                        from 'sliding off' the region of occupied bins\n");
  std::printf("  --boot                If present, then modify the FES by adding random\n");
  std::printf("                        numbers to the free energy values by drawing from\n");
  std::printf("                        normal distributions centered about the input free\n");
  std::printf("                        energy and whose standard deviation is the standard\n");
  std::printf("                        error of the free energy estimate\n");
  std::printf("  --boot-error ERR      When modifying the free energy values with --boot,\n");
  std::printf("                        ignore the standard errors and instead draw from\n");
  std::printf("                        normal distributions with the specified standard\n");
  std::printf("                        deviation (kcal/mol)\n");


  
  std::printf("\nPath Representation\n");
  std::printf("  --no-smooth           Do not smooth the spline control points when generating\n");
  std::printf("                        a new path. The default is to apply smoothing procedure.\n");
  std::printf("  --linear              Use piecewise linear splines. This is the default when\n");
  std::printf("                        using --msm --maxit=0\n");
  std::printf("  --akima               Use Akima splines. This is usually the default.\n");
  
  std::printf("  --fix0                If present, then don't update the position of t=0\n");
  std::printf("  --fix1                If present, then don't update the position of t=1\n");

  std::printf("  --npathpts NPP        Reconstruct a new initial guess from IPATH using NPP\n");
  std::printf("                        points. This is useful when making a simple guess but\n");
  std::printf("                        you want to optimize a path with more points. The\n");
  std::printf("                        default is to optimize with the same number of points\n");
  std::printf("                        as found in the IPATH file. If this option is unused\n");
  std::printf("                        and one is creating a new directory of simulations\n");
  std::printf("                        with --sasm, then it will default to 100\n");
  std::printf("                        If --npathpts is used, you should not use --pred-fc\n");
  std::printf("  --nsplpts NSP         Given a spline of the final path, interpolate the spline\n");
  std::printf("                        with NSP points to write OPATH.interp.model.dat.\n");
  std::printf("                        To print the free energy values at the input points,\n");
  std::printf("                        set: --maxit=0 and do not use --nsplpts nor --npathpts\n");

  std::printf("\nPath Optimization\n");
  std::printf("  -I MAXIT, --maxit MAXIT\n");
  std::printf("                        maximum number of synthetic iterations used to find the\n");
  std::printf("                        minimum free energy path on the fixed surface. If one is\n");
  std::printf("                        generating a new directory of simulations, then the default\n");
  std::printf("                        is 300. If one is analyzing a 1D projection, then the default\n");
  std::printf("                        is 0.\n");
  std::printf("  --nlayers N           The number of layers around the path used to select\n");
  std::printf("                        the next set of simulations.\n");

  std::printf("  --ptol PTOL           T-test p-value tolerance used during synthetic iterations.\n");
  std::printf("                        Values close to (but less than) 1.0 represent a strict\n");
  std::printf("                        tolerance requiring many iterations to achieve converence.\n");
  std::printf("                        A value >= 1 will never converge. A value <= 0 will\n");
  std::printf("  --acc-oobk            The force constant used to push the path toward the\n");
  std::printf("                        nearest occupied bin, if it is out-of-bounds. Default: 300.\n");
  std::printf("  --acc-maxit           The number of nonlinear sythetic optimizations steps.\n");
  std::printf("                        Default: 100\n");
  
  std::printf("                        pass the convergence test. Default: 0.95\n");
  std::printf("  --distol DISTOL       FES-FTSM slope termination tolerance for distance\n");
  std::printf("                        coordinates (default: 0.05 Ang).\n");
  std::printf("  --angtol ANGTOL       FES-FTSM slope termination tolerance for angle\n");
  std::printf("                        coordinates (default: 0.5 Deg).\n");
  std::printf("  --mlen MLEN           The number of iterations to consider when calculating\n");
  std::printf("                        the slope of the path with respect to iteration\n");
  std::printf("                        (default: 5).\n");


  std::printf("\nForce Constants\n");
  std::printf("  --predfc              Adjust the simulation force constants and progress\n");
  std::printf("                        variable values to place simulations along the path in a\n");
  std::printf("                        manner such that the expected centroid means are evenly\n");
  std::printf("                        spaced (rather than spacing the umbrella window centers).\n");
  std::printf("                        By default, the force constants are not adjusted and \n");
  std::printf("                        the simulations are chosen to uniformly space the\n");
  std::printf("                        biasing potential centers.\n");
  
  std::printf("  --smoothfc            If present, apply smoothing to the predicted force\n");
  std::printf("                        constants.\n");

  std::printf("  --tdispfc TDISP       Adjust the force constants to achieve the target\n");
  std::printf("                        displacement. Only used if --pred-fc is present.\n");
  std::printf("                        A value of 1 will scale the force constants\n");
  std::printf("                        such that the predicted centroids will displace along\n");
  std::printf("                        the path by the same amount as the nearest neighbors\n");
  std::printf("                        of a uniformly discretized set of windows. A value of\n");
  std::printf("                        0.5 will adjust the force constants so that the\n");
  std::printf("                        centroid does not displace by more than 1/2 the\n");
  std::printf("                        distance between nearest neighbors. Setting tdisp to\n");
  std::printf("                        small values will generally cause the force constants\n");
  std::printf("                        to increase. Default: 0.75.\n");
  
  std::printf("  --maxfc MAXFC         Maximum allowable force constants for distance\n");
  std::printf("                        restraints. Only used if --pred-fc is present.\n");
  std::printf("                        This value is the largest force constant that\n");
  std::printf("                        you feel comfortable simulating. Default: 300\n");
  std::printf("                        kcal/mol/A^2\n");
  std::printf("  --minfc MINFC         Minimum allowable force constants for distance\n");
  std::printf("                        restraints. Only used if --pred-fc is present.\n");
  std::printf("                        This value is the smallest force constant\n");
  std::printf("                        that you feel comfortable simulating. Default: 50\n");
  std::printf("                        kcal/mol/A^2\n");
  std::printf("  --maxfc-ang MAXFC_ANG\n");
  std::printf("                        Maximum allowable force constants for angle\n");
  std::printf("                        restraints. Only used if --pred-fc is present\n");
  std::printf("                        This value is the largest force constant that\n");
  std::printf("                        you feel comfortable simulating. Default: 0.5\n");
  std::printf("                        kcal/mol/deg^2\n");
  std::printf("  --minfc-ang MINFC_ANG\n");
  std::printf("                        Minimum allowable force constants for angle\n");
  std::printf("                        restraints. Only used if --pred-fc is present.\n");
  std::printf("                        This value is the smallest force constant\n");
  std::printf("                        that you feel comfortable simulating. Default: 0.1\n");
  std::printf("                        kcal/mol/deg^2\n");



  std::printf("\nMisc\n");
  
  std::printf("  -T TEMP, --temp TEMP  Temperature, K (default: 298.)\n");
  std::printf("  --print-freq FREQ     Output frequency during optimization (default: 50)\n");

} 



Options ReadOptions( int argc, char * argv[] )
{
  // Tell getopt to not print
  //extern int opterr;
  //opterr = 0;

  Options cli;

  
  static struct option long_options[] =
    {
     { "help",       no_argument,       NULL, 'h'  },
     { "disang",     required_argument, NULL, 'd'  },
     { "mdin",       required_argument, NULL, 'i'  },
     { "curit",      required_argument, NULL, 1000 },
     { "odir",       required_argument, NULL, 1001 },
     //{ "dacc",       no_argument,       NULL, 1002 },
     { "model",      required_argument, NULL, 'm'  },
     { "rbf",        no_argument,       NULL, 'r'  },
     { "bspl",       no_argument,       NULL, 'b'  },
     { "hist",       no_argument,       NULL, 'H'  },
     { "shape",      required_argument, NULL, 'R'  },
     { "minsize", required_argument, NULL, 1003 },
     { "rbfmaxerr",  required_argument, NULL, 1004 },
     { "maxit",      required_argument, NULL, 'I'  },
     { "distol",     required_argument, NULL, 1005 },
     { "angtol",     required_argument, NULL, 1006 },
     { "mlen",       required_argument, NULL, 1007 },
     //{ "closest",    no_argument,       NULL, 1008 },
     { "temp",       required_argument, NULL, 'T'  },
     { "fix0",       no_argument,       NULL, 1009 },
     { "fix1",       no_argument,       NULL, 1010 },
     //{ "dampfactor", required_argument, NULL, 1011 },
     //{ "smooth-frac",required_argument, NULL, 1012 },
     //{ "smooth-nmin",required_argument, NULL, 1013 },
     //{ "smooth-nmax",required_argument, NULL, 1014 },
     { "smoothfc",  no_argument,       NULL, 1015 },
     { "scalefc",   no_argument,       NULL, 1016 },
     { "predfc",    no_argument,       NULL, 1017 },
     { "maxfc",      required_argument, NULL, 1018 },
     { "minfc",      required_argument, NULL, 1019 },
     { "maxfc-ang",  required_argument, NULL, 1020 },
     { "minfc-ang",  required_argument, NULL, 1021 },
     { "tdispfc",      required_argument, NULL, 1022 },
     //{ "pkl",        no_argument,       NULL, 1023 },
     //{ "nquad",      required_argument, NULL, 1024 },
     { "cp-nearest", no_argument,       NULL, 1025 },
     { "arbf",       required_argument, NULL, 1026 },
     { "ipath",       required_argument, NULL, 1027 },
     { "opath",      required_argument, NULL, 1028 },
     { "chk",        required_argument, NULL, 1029 },
     { "boot",       no_argument,       NULL, 1030 },
     { "boot-error", required_argument, NULL, 1031 },
     { "npathpts",   required_argument, NULL, 1032 },
     { "nsplpts",    required_argument, NULL, 1033 },
     //{ "nebk",       required_argument, NULL, 1035 },
     { "ptol",       required_argument, NULL, 1036 },
     //{ "conv-layers",required_argument, NULL, 1037 },
     //{ "conv-samples",required_argument,NULL, 1038 },
     //{ "tetris",     no_argument,       NULL, 1039 },
     //{ "tetris-smooth",no_argument,     NULL, 1040 },
     { "acc",        no_argument,       NULL, 1041 },
     { "acc-oobk",   required_argument, NULL, 1042 },
     { "acc-maxit",  required_argument, NULL, 1043 },
     { "buffer",     no_argument,       NULL, 1044 },
     { "nlayers",    required_argument, NULL, 1045 },
     { "dry-run",    no_argument,       NULL, 1046 },
     //{ "explore-once",no_argument,      NULL, 1047 },
     { "neqit",      required_argument, NULL, 1048 },
     //{ "disp-limit", no_argument,       NULL, 1049 },
     { "wavg",       required_argument, NULL, 1050 },
     { "wavg-niter", required_argument, NULL, 1051 },
     { "prefix",     required_argument, NULL, 1052 },
     { "extra-prefix", required_argument, NULL, 1053 },
     { "no-smooth",  no_argument, NULL, 1054 },
     { "msm",        no_argument, NULL, 1055 },
     { "gbsm",       no_argument, NULL, 1056 },
     { "linear",     no_argument, NULL, 1057 },
     { "akima",      no_argument, NULL, 1058 },
     { "pad",        required_argument, NULL, 1059 },
     { "print-freq", required_argument, NULL, 1060 },
     { "sasm",     no_argument, NULL, 1061 },

     {NULL,0,NULL,0}
    };


  int opt = 0;
  int long_index = 0;
  //char * subopts, * value;
  while ( (opt = getopt_long
          ( argc, argv, "hd:i:m:rbHR:I:T:", 
            long_options, &long_index )) != -1 )
    {

      switch (opt)
        {
	  
        case 'h':
          cli.print_help();
          std::exit(EXIT_SUCCESS);
          break;
	  
        case 'd':  { cli.disang = optarg; break; }
        case 'i':  { cli.mdin = optarg; break; }
        case 1000: { cli.curit = std::atoi(optarg); break; }
        case 1001: { cli.odir = optarg; break; }
	  //case 1002: { cli.dacc = true; break; }
        case 'm':  { cli.model = std::atoi(optarg); break; }
        case 'r':  { cli.rbf = true; break; }
        case 'b':  { cli.bspl = true; break; }
        case 'H':  { cli.hist = true; break; }
        case 'R':  { cli.shape = std::atof(optarg); break; }
        case 1003: { cli.minsize = std::atoi(optarg); break; }
        case 1004: { cli.rbfmaxerr = std::atof(optarg); break; }
        case 'I':  { cli.maxit = std::atoi(optarg); break; }
        case 1005: { cli.distol = std::atof(optarg); break; }
        case 1006: { cli.angtol = std::atof(optarg); break; }
        case 1007: { cli.mlen = std::atoi(optarg); break; }
	  //case 1008: { cli.closest = true; break; }
        case 'T':  { cli.temp = std::atof(optarg); break; }
        case 1009: { cli.fix0 = true; break; }
        case 1010: { cli.fix1 = true; break; }
	  //case 1011: { cli.dampfactor = std::atof(optarg); break; }
	  //case 1012: { cli.smooth_frac = std::atof(optarg); break; }
	  //case 1013: { cli.smooth_nmin = std::atoi(optarg); break; }
	  //case 1014: { cli.smooth_nmax = std::atoi(optarg); break; }
        case 1015: { cli.smoothfc = true; break; }
	  //case 1016: { cli.scalefc = true; break; }
        case 1017: { cli.predfc = true; break; }
        case 1018: { cli.maxfc = std::atof(optarg); break; }
        case 1019: { cli.minfc = std::atof(optarg); break; }
        case 1020: { cli.maxfc_ang = std::atof(optarg); break; }
        case 1021: { cli.minfc_ang = std::atof(optarg); break; }
        case 1022: { cli.tdispfc = std::atof(optarg); break; }
	  //case 1023: { cli.pkl = true; break; }
	  //case 1024: { cli.nquad = std::atoi(optarg); break; }
        case 1025: { cli.cp_nearest = true; break; }
        case 1026: { cli.arbf = std::atoi(optarg); break; }
        case 1027: { cli.ipath = optarg; break; }
        case 1028: { cli.opath = optarg; break; }
        case 1029: { cli.chk = optarg; break; }
	case 1030: { cli.boot = true; break; }
	case 1031: { cli.booterror = std::atof(optarg); break; }
	case 1032: { cli.npathpts = std::atoi(optarg); break; }
	case 1033: { cli.nsplpts = std::atoi(optarg); break; }
	  //case 1035: { cli.nebk = std::atof(optarg); break; }
	case 1036: { cli.ptol = std::atof(optarg); break; }
	  //case 1037: { cli.conv_layers = std::atoi(optarg); break; }
	  //case 1038: { cli.conv_samples = std::atoi(optarg); break; }
	  //case 1039: { cli.tetris = true; break; }
	  //case 1040: { cli.tetris_smooth = true; break; }
	case 1041: { cli.acc = true; break; }
	case 1042: { cli.acc_oobk = std::atof(optarg); break; }
	case 1043: { cli.acc_maxit = std::atoi(optarg); break; }
	case 1044: { cli.buffer = true; break; }
        case 1045: { cli.nlayers = std::atoi(optarg); break; }
	case 1046: { cli.dry_run = true; break; }
	  //case 1047: { cli.explore_once = true; break; }
        case 1048: { cli.neqit = std::atoi(optarg); break; }
	  //case 1049: { cli.disp_limit = true; break; }
        case 1050: { cli.wavg = std::atoi(optarg); break; }
        case 1051: { cli.wavg_niter = std::atoi(optarg); break; }
        case 1052: { cli.prefix = optarg; break; }
        case 1053:
	  {
	    std::string sdir = optarg;
	    bool found=false;
	    for ( std::size_t ii=0; ii<cli.extra_prefixes.size(); ++ii )
	      {
		if ( cli.extra_prefixes[ii].compare(sdir) == 0 )
		  {
		    found=true;
		    break;
		  }
	      };
	    if ( not found )
	      {
		cli.extra_prefixes.push_back(sdir);
	      };
	    break;
	  }
	case 1054: { cli.smooth = false; break; }
	case 1055: { cli.msm = true; break; }
	case 1056: { cli.gbsm = true; break; }
	case 1057:
	  {
	    if ( cli.stype >= 0 )
	      {
		std::printf("Encountered --linear, but the spline type has already been set\n");
		std::exit(EXIT_FAILURE);
	      }
	    cli.stype = 0;
	    break;
	  }
	case 1058:
	  {
	    if ( cli.stype >= 0 )
	      {
		std::printf("Encountered --linear, but the spline type has already been set\n");
		std::exit(EXIT_FAILURE);
	      }
	    cli.stype = 1;
	    break;
	  }
        case 1059: { cli.pad = std::atoi(optarg); break; }
        case 1060: { cli.printfreq = std::atoi(optarg); break; }
	case 1061: { cli.sasm = true; break; } 
	  //case 1062: { cli.gxsm = true; break; } 
        case '?': // getopt_long already printed an error
          std::printf("%s: use -h for usage\n", argv[0]);
          std::exit(EXIT_FAILURE);
          break;

        default:
          std::printf("An error occured while reading command line options\n");
          std::exit(EXIT_FAILURE);
        }
    };



  for ( int iarg=0; optind<argc; ++optind, ++iarg )
    {
      if ( iarg > 0 )
	{
	  std::printf("Error: Unexpected command line argument: %s\n",argv[optind]);
	  std::exit(EXIT_FAILURE);
	}
      //cli.metafile = argv[optind];
    }

  if ( cli.ptol >= 1 )
    {
      std::printf("Error: --ptol must be less than 1.0\n");
      std::exit(EXIT_FAILURE);
    }
  
  if ( cli.ptol < 0 )
    {
      std::printf("Error: --ptol must be nonnegative\n");
      std::exit(EXIT_FAILURE);
    }
  
  std::size_t npath = (cli.ipath.size()>0) + (cli.opath.size()>0)
    + (cli.chk.size()>0);
  
  if ( npath == 0 )
    {
      if ( cli.boot )
	{
	  std::printf("Error: Cannot use --boot option when not using --ipath --opath --chk\n");
	  std::exit(EXIT_FAILURE);
	}
      if ( cli.rbf + cli.hist + cli.bspl + (cli.arbf>0) + (cli.wavg>1) > 1 )
	{
	  std::printf("Error: --rbf, --arbf, --hist, --wavg, and --bspl are mutually exclusive\n");
	  std::exit(EXIT_FAILURE);
	}
      
      if ( cli.disang.size() == 0 )
	{
	  std::printf("Error: --disang is a required argument\n");
	  std::exit(EXIT_FAILURE);
	}
      
      if ( cli.mdin.size() == 0 )
	{
	  std::printf("Error: --mdin is a required argument\n");
	  std::exit(EXIT_FAILURE);
	}
      
      if ( cli.curit < 0 )
	{
	  std::printf("Error: --curit is a required argument\n");
	  std::exit(EXIT_FAILURE);
	}

    }
  else if ( npath == 3 )
    {
      if ( cli.disang.size() > 0 )
	{
	  std::printf("Error: --disang cannot be used with --ipath --opath --chk\n");
	  std::exit(EXIT_FAILURE);
	}
      if ( cli.mdin.size() > 0 )
	{
	  std::printf("Error: --mdin cannot be used with --ipath --opath --chk\n");
	  std::exit(EXIT_FAILURE);
	}
      if ( cli.curit >= 0 )
	{
	  std::printf("Error: --curit cannot be used with --ipath --opath --chk\n");
	  std::exit(EXIT_FAILURE);
	}
      
    }
  else
    {
      std::printf("Error: Either all of these options must be used or none of them: --ipath --opath --chk\n");
      std::exit(EXIT_FAILURE);
    }


  // if ( cli.nquad < 1 )
  //   {
  //     std::printf("--nquad must be >= 1\n");
  //     std::exit(EXIT_FAILURE);
  //   }
  
  // cli.smooth_nmax = std::max(cli.smooth_nmax,cli.smooth_nmin);
  // if ( cli.smooth_nmin >= 3 )
  //   {
  //     if ( cli.smooth_nmin % 2 == 0 )
  // 	{
  // 	  std::printf("--smooth-nmin must be an odd integer\n");
  // 	  std::exit(EXIT_FAILURE);
  // 	}
  //     if ( cli.smooth_nmax % 2 == 0 )
  // 	{
  // 	  std::printf("--smooth-nmax must be an odd integer\n");
  // 	  std::exit(EXIT_FAILURE);
  // 	}
  //   }

  // cli.dampfactor = std::min(1.,std::max(0.,cli.dampfactor));

  if ( cli.acc_maxit < 1 )
    {
      cli.acc = false;
    }
  cli.acc_oobk = std::max( cli.acc_oobk, 0. );
  // if ( cli.acc and cli.dacc )
  //   {
  //     std::printf("--dacc and --acc cannot be used at the same time\n");
  //     std::exit(EXIT_FAILURE);
  //   }
  // if ( cli.acc and cli.tetris )
  //   {
  //     std::printf("--tetris and --acc cannot be used at the same time\n");
  //     std::exit(EXIT_FAILURE);
  //   }
  
  // cli.sim_layers = std::max(cli.sim_layers,cli.conv_layers);

  
  if ( cli.prefix.find('/') != std::string::npos )
    {
      std::printf("--prefix cannot contain a slash: '/'\n");
      std::exit(EXIT_FAILURE);
    }

  if ( cli.prefix.find(' ') != std::string::npos )
    {
      std::printf("--prefix cannot contain a space: ' '\n");
      std::exit(EXIT_FAILURE);
    }

  if ( cli.msm + cli.gbsm + cli.sasm > 1 )
    {
      std::printf("Can only use one of --msm --gbsm --sasm\n");
      std::exit(EXIT_FAILURE);
    }

  if ( cli.gbsm and cli.nlayers < 0 )
    {
      cli.nlayers = 2;
    }

  if ( cli.stype < 0 )
    {
      if ( cli.msm and cli.maxit < 1 )
	{
	  cli.stype = 0;
	}
      else
	{
	  cli.stype = 1;
	};
    }

  if ( cli.maxit < 0 )
    {
      if ( npath == 3 )
	{
	  cli.maxit = 0;
	}
      else
	{
	  cli.maxit = 300;
	}
    }

  cli.printfreq = std::max(1,cli.printfreq);

  cli.pad = std::max(1,cli.pad);


  if ( cli.npathpts < 1 )
    {
      if ( npath != 3 and (cli.gbsm or cli.sasm) )
	{
	  cli.npathpts = 100;
	}
    };
  
  // if ( npath == 3 )
  //   {
  //     cli.gbsm = false;
  //     cli.msm = true;
  //   }
  return cli;
}
