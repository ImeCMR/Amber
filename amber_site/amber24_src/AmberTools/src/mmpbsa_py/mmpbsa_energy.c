/* mmpbsa_energies.nab: NAB program to perform MM/PBSA and MM/GBSA analyses on 
 *                      a trajectory file akin to sander's imin=5 ability. The
 *                      command-line usage and output file format is designed 
 *                      to be close enough to sander's so that the interface and
 *                      parsers for MMPBSA.py do not have to change to accomodate
 *                      them
 *
 * By Jason Swails, Dwight McGee, Bill Miller III, Adrian Roitberg
 * 2/22/2011
 */
/* Conversion to nabc initiated by Dave Case, Sept. 2022  */

#include "nabc.h"
#include "AmberNetcdf.h"
#include "sff.h"

FILE* nabout;


// Some global variable declarations

struct CLOptions {
   char *mdin, *prmtop, *inpcrd, *traj, *mdout;
};

struct GBOptions {
   int igb, gbsa, rungb;
   double extdiel, saltcon, surften, rgbmax;
};

struct PBOptions {
   int runpb, inp, smoothopt, radiopt, npbopt, solvopt, maxitn;
   int nfocus, bcopt, eneopt, fscale, dbfopt;
   double epsin, epsout, istrng, dprob, iprob, accept, fillratio;
   double space, cutnb, sprob, cavity_surften, cavity_offset;
};


static struct GBOptions gbopt;
static struct PBOptions pbopt;
static struct CLOptions clopt;

// GB input file reading function
int ReadGBin(char* filename) {
   char option[50], *line, value[50];
   FILE* gbfile;

   gbfile = fopen(filename, "r");
   if ( gbfile == NULL ) {
      fprintf(stderr, "Error: Cannot open %s for reading!\n", filename);
      return -1;
   }

   // Make sure it's a GB input file
   line = NULL;
   size_t len = 0;
   getline(&line, &len, gbfile);
   if ( strncmp(line,"GB",2) )
      return 1;

   gbopt.rungb = 1;

   while (getline(&line, &len, gbfile) > -1) {
      sscanf(line, "%s = %s", option, value);
      if      ( ! strcmp(option,"igb")) gbopt.igb = atoi(value);
      else if ( ! strcmp(option,"gbsa")) gbopt.gbsa = atoi(value);
      else if ( ! strcmp(option,"extdiel")) gbopt.extdiel = atof(value);
      else if ( ! strcmp(option,"saltcon")) gbopt.saltcon = atof(value);
      else if ( ! strcmp(option,"surften")) gbopt.surften = atof(value);
      else if ( ! strcmp(option,"rgbmax")) gbopt.rgbmax = atof(value);
      else {
         fprintf(stderr, "Error: Unknown GB option %s!", option);
         fclose(gbfile);
         return -1;
      }
   }

   return fclose(gbfile);

};


// PB input file reading function
int ReadPBin(char* filename) {
   char option[50], *line, value[50];
   FILE* pbfile;
   int nread;

   pbfile = fopen(filename, "r");
   if ( pbfile == NULL ) {
      fprintf(stderr, "Error: Cannot open %s for reading!\n", filename);
      return -1;
   }

   // Check to make sure it's a PB input file
   line = NULL;
   size_t len = 0;
   getline(&line, &len, pbfile);
   if ( strncmp(line,"PB",2) ){
      fprintf(stderr, "Error: Cannot determine type of input file (%s)!",
              filename);
      return 1;
   }

   pbopt.runpb = 1;

   while (getline(&line, &len, pbfile) > -1) {
      sscanf(line, "%s = %s", option, value);
      if      ( ! strcmp(option,"inp") ) pbopt.inp = atoi(value);
      else if ( ! strcmp(option,"smoothopt") ) pbopt.smoothopt = atoi(value);
      else if ( ! strcmp(option,"radiopt") ) pbopt.radiopt = atoi(value);
      else if ( ! strcmp(option,"npbopt") ) pbopt.npbopt = atoi(value);
      else if ( ! strcmp(option,"solvopt") ) pbopt.solvopt = atoi(value);
      else if ( ! strcmp(option,"maxitn") ) pbopt.maxitn = atoi(value);
      else if ( ! strcmp(option,"nfocus") ) pbopt.nfocus = atoi(value);
      else if ( ! strcmp(option,"bcopt") ) pbopt.bcopt = atoi(value);
      else if ( ! strcmp(option,"eneopt") ) pbopt.eneopt = atoi(value);
      else if ( ! strcmp(option,"fscale") ) pbopt.fscale = atoi(value);
      else if ( ! strcmp(option,"dbfopt") ) pbopt.dbfopt = atoi(value);
      else if ( ! strcmp(option,"epsin") ) pbopt.epsin = atof(value);
      else if ( ! strcmp(option,"epsout") ) pbopt.epsout = atof(value);
      else if ( ! strcmp(option,"istrng") ) pbopt.istrng = atof(value);
      else if ( ! strcmp(option,"dprob") ) pbopt.dprob = atof(value);
      else if ( ! strcmp(option,"iprob") ) pbopt.iprob = atof(value);
      else if ( ! strcmp(option,"accept") ) pbopt.accept = atof(value);
      else if ( ! strcmp(option,"fillratio") ) pbopt.fillratio = atof(value);
      else if ( ! strcmp(option,"space") ) pbopt.space = atof(value);
      else if ( ! strcmp(option,"cutnb") ) pbopt.cutnb = atof(value);
      else if ( ! strcmp(option,"sprob") ) pbopt.sprob = atof(value);
      else if ( ! strcmp(option,"cavity_surften") ) pbopt.cavity_surften = atof(value);
      else if ( ! strcmp(option,"cavity_offset") ) pbopt.cavity_offset = atof(value);
      else {
         fprintf(stderr, "Error: Unknown PB option %s!\n", option);
         fclose(pbfile);
         return -1; 
      }

   }

   return (fclose(pbfile));

};

// Sets the default values for each set
int SetDefaults() {
   
   // First set filename defaults
   clopt.mdin = strdup("mdin");
   clopt.prmtop = strdup("prmtop");
   clopt.inpcrd = strdup("pdb");
   clopt.traj = strdup("mdcrd");
   clopt.mdout = strdup("mdout");

   // Next set default GB Options
   gbopt.rungb = 0;
   gbopt.igb = 1;
   gbopt.gbsa = 0;
   gbopt.extdiel = 78.5;
   gbopt.saltcon = 0.0;
   gbopt.surften = 0.0072;
   gbopt.rgbmax = 999.0;

   // Next set default PB Options
   pbopt.runpb = 0;
   pbopt.inp = 2;
   pbopt.smoothopt = 1;
   pbopt.radiopt = 1;
   pbopt.npbopt = 0;
   pbopt.solvopt = 1;
   pbopt.maxitn = 100;
   pbopt.nfocus = 2;
   pbopt.fscale = 8;
   pbopt.dbfopt = 1;
   pbopt.epsin = 1.0;
   pbopt.epsout = 80.0;
   pbopt.istrng = 0.0;
   pbopt.dprob = 1.4;
   pbopt.iprob = 2.0;
   pbopt.accept = 0.001;
   pbopt.fillratio = 2.0;
   pbopt.space = 0.5;
   pbopt.bcopt = 5;
   pbopt.eneopt = 2;
   pbopt.cutnb = 0.0;
   pbopt.sprob = 0.557;
   pbopt.cavity_surften = 0.0378;
   pbopt.cavity_offset = -0.5692;

   return 0;
};

// Prints usage statement of the program
int printusage(char* prog_name) {
   printf(" Usage: %s [-O] -i mdin -o mdout -p prmtop -c pdb -y mdcrd\n", prog_name);
   return 0;
};


// Checks to make sure that allowed values are used
int CheckValues() {
   int isinerr;
   isinerr = 0;
   if (gbopt.rungb) {
      if (gbopt.igb != 1 && gbopt.igb != 2 && gbopt.igb != 5 && gbopt.igb !=7 && gbopt.igb !=8) {
         fprintf(stderr, "Error: IGB must be 1, 2, 5, 7, or 8!\n");
         isinerr = 1;
    } if (gbopt.extdiel <= 0) {
         fprintf(stderr, "Error: EXTDIEL must be positive!\n");
         isinerr = 1;
    } if (gbopt.surften < 0) {
         fprintf(stderr, "Error: SURFTEN must be non-negative!\n");
         isinerr = 1;
    } if (gbopt.saltcon < 0) {
         fprintf(stderr, "Error: SALTCON must be non-negative!\n");
         isinerr = 1;
    } if (gbopt.rgbmax <= 0) {
         fprintf(stderr, "Error: RGBMAX must be positive!\n");
         isinerr = 1;
    } if (gbopt.rgbmax < 20 && gbopt.rgbmax > 0) {
         fprintf(stderr, "Warning: Low value for RGBMAX. Consider using default.\n");
      }
 } else {
      if (pbopt.inp != 0 && pbopt.inp != 1 && pbopt.inp != 2) {
         fprintf(stderr, "Error: INP must be 0, 1, or 2!\n");
         isinerr = 1;
    } if ( pbopt.inp == 1 ) {
         fprintf(stderr,"Warning: inp=1 was old default\n");
         if (floor(10000.0 * pbopt.cavity_surften) == floor(10000.0 * 0.0378)) {
            fprintf(stderr, "Warning: cavity_surften=0.0378 not recommended for inp=1, switching to inp=1 default value: 0.0050\n");
            pbopt.cavity_surften = 0.00500;
         }
         if (floor(10000.0 * pbopt.cavity_offset) == floor(-0.5692 * 10000.0)) {
            fprintf(stderr, "Warning: cavity_offset=-0.5692 not recommended for inp=1, switching to inp=1 default value: 0.000\n");
            pbopt.cavity_offset = 0.000;
         }
         if (floor(pbopt.sprob * 10000.0) == floor(0.557 * 10000.0)) {
            fprintf(stderr, "Warning: sprob=.557 not recommended for inp=1, switching to inp=1 default value: 1.400\n");
            pbopt.sprob = 1.400;
         } 
         if (pbopt.radiopt != 0 && pbopt.inp == 1) {
            fprintf(stderr, "Warning: radiopt should be set to 0 for inp=1\n");
            pbopt.radiopt = 0;
         }
    } if (pbopt.smoothopt != 0 && pbopt.smoothopt != 1 && pbopt.smoothopt != 2){
         fprintf(stderr, "Error: SMOOTHOPT must be 0, 1, or 2!\n");
         isinerr = 1;
    } if (pbopt.radiopt != 1 && pbopt.radiopt != 0) {
         fprintf(stderr, "Error: RADIOPT must be 0, 1, or 2!\n");
         isinerr = 1;
    } if (pbopt.npbopt != 0 && pbopt.npbopt != 1) {
         fprintf(stderr, "Error: NPBOPT must be 0 or 1!\n");
         isinerr = 1;
    } if (pbopt.solvopt < 1 || pbopt.solvopt == 7 || pbopt.solvopt > 8) {
         fprintf(stderr, "Error: SOLVOPT must be 1, 2, 3, 4, 5, 6, or 8!\n");
         isinerr = 1;
    } if (pbopt.maxitn < 1) {
         fprintf(stderr, "Error: MAXITN must be a positive integer!\n");
         isinerr = 1;
    } if (pbopt.nfocus != 1 && pbopt.nfocus != 2) {
         fprintf(stderr, "Error: NFOCUS must be 1 or 2!\n");
         isinerr = 1;
    } if (pbopt.fscale <= 0) {
         fprintf(stderr, "Error: NFOCUS must be non-negative!\n");
         isinerr = 1;
    } if (pbopt.epsin < 0 || pbopt.epsout < 0) {
         fprintf(stderr, "Error: EPSIN/OUT must be non-negative!\n");
         isinerr = 1;
    } if (pbopt.istrng < 0) {
         fprintf(stderr, "Error: ISTRNG must be non-negative!\n");
         isinerr = 1;
    } if (pbopt.dprob < 0) {
         fprintf(stderr, "Error: DPROB must be non-negative!\n");
         isinerr = 1;
    } if (pbopt.iprob < 0) {
         fprintf(stderr, "Error: IPROB must be non-negative!\n");
         isinerr = 1;
    } if (pbopt.accept <= 0) {
         fprintf(stderr, "Error: ACCEPT must be positive!\n");
         isinerr = 1;
    } if (pbopt.fillratio <= 0) {
         fprintf(stderr, "Error: FILLRATIO must be positive!\n");
         isinerr = 1;
    } if (pbopt.space <= 0) {
         fprintf(stderr, "Error: SPACE must be positive!\n");
         isinerr = 1;
    } if (pbopt.bcopt != 1 && pbopt.bcopt != 5 && pbopt.bcopt != 6 && pbopt.bcopt != 10) {
         fprintf(stderr, "Error: BCOPT must be 1, 5, 6, or 8!\n");
         isinerr = 1;
    } if (pbopt.eneopt != 1 && pbopt.eneopt != 2) {
         fprintf(stderr, "Error: ENEOPT must be 1 or 2!\n");
         isinerr = 1;
    } if (pbopt.cutnb < 0) {
         fprintf(stderr, "Error: CUTNB must be non-negative!\n");
         isinerr = 1;
    } if (pbopt.sprob < 0) {
         fprintf(stderr, "Error: SPROB must be non-negative!\n");
         isinerr = 1;
    } if (pbopt.cavity_surften < 0 ) {
         fprintf(stderr, "Error: CAVITY_SURFTEN must be non-negative!\n");
         isinerr = 1;
      }
   }
   return (isinerr);
};


/*
 *       BEGIN MAIN PROGRAM
 */

int main( int argc, char* argv[] )
{

nabout = stdout;

// Variable declaration
int i, j, fr_result, frame, trajDone;
double kappa, cut;
double *coords, *grad;
struct AmberNetcdf nctraj;
FILE *output, *asciitraj;
char *prog_name;

/* Variable Descriptions:
 * 
 * i, j        : counter variables
 * fr_result   : File Read RESULT (return code of ReadGB/PBin
 * frame       : keeps track of frame number we're doing for iteration counter
 * trajDone    : flag to indicate if we're done getting frames from ASCII traj.
 * coords      : x, y, z-coordinate array
 * grad        : x, y, z-gradient array for mme call
 * nctraj      : NetCDF trajectory struct
 * output      : File to redirect nabout to
 * asciitraj   : ASCII trajectory file (in case it's not netcdf)
 * line        : string buffer to hold a line pulled from asciitraj
 * prog_name   : name of the program
 * mol         : system molecule
 */

/* First parse the command-line; -r is a dummy flag to maintain
 * sander-like behavior
 */

SetDefaults();

prog_name = argv[0];
i = 1;
while (i < argc) {
   if ( ! strcmp(argv[i],"-i") )
      clopt.mdin = argv[++i];
   else if ( ! strcmp(argv[i],"-p") )
      clopt.prmtop = argv[++i];
   else if ( ! strcmp(argv[i],"-c") )
      clopt.inpcrd = argv[++i];
   else if ( ! strcmp(argv[i],"-o") )
      clopt.mdout = argv[++i];
   else if ( ! strcmp(argv[i],"-y") )
      clopt.traj = argv[++i];
   else if ( ! strcmp(argv[i],"-r") )
      i++;
   else if ( ! strcmp(argv[i],"-h" ) || ! strcmp(argv[i],"--help" ) || 
             ! strcmp(argv[i],"--h" ) ) {
      printusage(prog_name);
      exit(0);
  }else if ( strcmp( argv[i],"-O") ) {
      fprintf(stderr, "Error: Bad flag %s!\n", argv[i]);
      printusage(prog_name);
      exit(1);
   }

   i++;
}

/* First try to read the mdin file as a GB input file. If ReadGBin returns 1,
 * then we know it's a PB input file, so try to read it as a PB input file.
 * If that doesn't work and/or ReadGBin returns -1, we know we've reached an
 * error and should quit
 */
fr_result = ReadGBin(clopt.mdin);

if (fr_result == -1)
   exit(-1);
else if (fr_result == 1) {
   if (ReadPBin(clopt.mdin) != 0)
      exit(-1);
}

// Check to make sure that all values used are allowed
if (CheckValues() == 1) exit(1);

// Make nabout printing to go output file
output = fopen(clopt.mdout, "w");
nabout = output;

// Initialize the molecule and options
PARMSTRUCT_T *prm;
prm = rdparm( clopt.prmtop );
coords = malloc( prm->Nat3 * (sizeof(double)) );
grad = malloc( prm->Nat3 * (sizeof(double)) );

// Set the mm_options and call mm_init
char *q;
if (gbopt.rungb == 1) {
   kappa = sqrt(0.10806*gbopt.saltcon);
   // Set default cut at 999.0, but make sure it's >= rgbmax
   cut = 999.0;
   if (gbopt.rgbmax > cut) cut = gbopt.rgbmax;

   asprintf(&q, "e_debug=2, gb=%d, rgbmax=%lf, gbsa=%d, surften=%lf, cut=%lf", 
              gbopt.igb, gbopt.rgbmax, gbopt.gbsa, gbopt.surften, cut);
   mm_options( q );

   asprintf(&q, "epsext=%lf, kappa=%lf", gbopt.extdiel, kappa);
   mm_options( q );

} else {
   asprintf(&q, "e_debug=3, ipb=2, inp=%d, epsin=%lf, epsout=%lf, smoothopt=%d",
              pbopt.inp, pbopt.epsin, pbopt.epsout, pbopt.smoothopt);
   mm_options( q );

   asprintf(&q, "istrng=%lf, radiopt=%d, dprob=%lf, iprob=%lf, npbopt=%d",
         pbopt.istrng, pbopt.radiopt, pbopt.dprob, pbopt.iprob, pbopt.npbopt);
   mm_options( q );

   asprintf(&q, "solvopt=%d, accept=%lf, maxitn=%d, fillratio=%lf, space=%lf",
     pbopt.solvopt, pbopt.accept, pbopt.maxitn, pbopt.fillratio, pbopt.space);
   mm_options( q );

   asprintf(&q, "nfocus=%d, fscale=%d, bcopt=%d, eneopt=%d, cutnb=%lf, sprob=%lf",
            pbopt.nfocus, pbopt.fscale, pbopt.bcopt, pbopt.eneopt, pbopt.cutnb, 
            pbopt.sprob);
   mm_options( q );

   asprintf(&q, "cavity_surften=%lf, cavity_offset=%lf", 
            pbopt.cavity_surften, pbopt.cavity_offset);
   mm_options( q );

}

   // nothing frozen or constrained for now:
   int* frozen = parseMaskString( "@ZZZ", prm, coords, 2 );
   int* constrained = parseMaskString( "@ZZZ", prm, coords, 2 );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );


// Now begin to loop through trajectory snapshots
frame = 0;
if ( netcdfLoad(&nctraj, clopt.traj) == 0 ) {
   netcdfLoad(&nctraj, clopt.traj);
   fprintf(nabout, "Processing NetCDF trajectory (%s)\n\n", clopt.traj);
   while( netcdfGetNextFrame(&nctraj, coords, grad) ) {
      fprintf(nabout, "Processing frame %d\n", nctraj.currentFrame);
      mme(coords, grad, &frame);
      fprintf(nabout, "\n");
      frame++;
   }
   netcdfClose(&nctraj);
} else { // Now we assume ASCII
   fprintf(nabout, "Processing ASCII trajectory (%s)\n\n", clopt.traj);
   asciitraj = fopen(clopt.traj, "r");
   
   if (asciitraj == NULL) {
      fprintf(stderr, "Error: Could not open trajectory file (%s) for reading!", clopt.traj);
      exit(-1);
   }

   trajDone = 0;
   char *line = NULL;
   size_t len = 0;
   getline(&line, &len, asciitraj); // eat the title line

   for (i = 1;;i++) { // frame loop
      for (j = 0; j < prm->Nat3; j++) { // coord. loop
         if(fscanf(asciitraj, "%lf", &coords[j]) < 1) {
            trajDone = 1; // if we've reached here we can go no farther
            break;
         }
      }
      if (trajDone) break; // if we finished, quit out of the loop now!
  
      fprintf(nabout, "Processing frame %d\n", i);
      mme(coords, grad, &frame);
      fprintf(nabout, "\n");
      frame++;
   } 
}

if ( gbopt.rungb )
   fprintf(nabout, "MM/GBSA processing done!\n");
else
   fprintf(nabout, "MM/PBSA processing done!\n");


exit(fclose(output));

}
