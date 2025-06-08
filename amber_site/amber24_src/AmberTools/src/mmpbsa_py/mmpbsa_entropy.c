/* Nab program to calculate the entropy {program of MMPBSA.py}
   Last edit 3/7/11
   Written by Dwight McGee, Bill Miller III, Jason Swails */
/* Conversion to nabc initiated by Dave Case, Sept. 2022  */

#include <nabc.h>
#include <AmberNetcdf.h>

FILE* nabout;

//-----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{

nabout = stdout;

//Define All Variables

PARMSTRUCT_T *prm;
XMIN_OPT_T xo;
struct AmberNetcdf  nc_ascii;

int      natm;
int      framenum, traj_complete;
int      i, j;
double   energy, min_energy;           // energy returned, miminization energy
FILE     *asciitraj;                   //trajectory file (ASCII format)
//---------------------------------------------------------------------------------

/* Usage: mmpbsa_entropy.nab {pdb} {prmtop} {maxiter} {drms} {'string of MM options'} {nc/trj}
                           argv[1] argv[2]   argv[3]  argv[4]    argv[5]              argv[6]       
  Double Check to ensure MMPBSA.py is giving the correct amount of arguments */  
if(argc !=7){
   printf("    Bad command-line arguments\n");
   printf("Usage: %s {pdb} {prmtop} {maxiter} {drms} {'string of MM options'} {nc/trj}\n", argv[0]);
   exit(1);
}
//--------------------------------------------------------------------------------
// Define and Setup Necessary Parameters for MM/Min

prm = rdparm( argv[2] );         // read in prmtop from command line
natm=prm->Natom;                   // number of atoms
double* xyz = malloc( prm->Nat3 * (sizeof(double)) );
double* grad = malloc( prm->Nat3 * (sizeof(double)) );

//MM parameters
mm_options(argv[5]);               // set parameters for MM/Min calc

   // nothing frozen or constrained for now:
   int* frozen = parseMaskString( "@ZZZ", prm, xyz, 2 );
   int* constrained = parseMaskString( "@ZZZ", prm, xyz, 2 );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );

/*NOTE: mm_options() must be called before mme_init()
  NOTE: mme_init() must be called before mme() function
  NOTE: :-) means convergence 
  NOTE: :-( indicates CG iteration encountered negatitve curvature and had to abort */   

//XMIN Minimization parameters
xmin_opt_init(&xo);               // initialize xmin optimization
xo.maxiter     = atoi(argv[3]);   // max number of iterations allowed for XMIN
xo.grms_tol    = atof(argv[4]);   // convergence criterion
xo.method      = 3;               // Trunated Newton Conjugate gradient algorithm
xo.numdiff     = 1;               // method used to approximate the product of the Hessian matrix
xo.m_lbfgs     = 3;               // size of the L-BFG memory 
xo.ls_method   = 2;               // Do not change (line-search method)
xo.ls_maxiter  = 20;              // max number of line search steps per single minimization step
xo.ls_maxatmov = 0.15;            // Maximum (co-ordinate) movement per F.O.D allowed in line search
xo.print_level = 2;               // highest print verbose
//--------------------------------------------------------------------------------

//Begin processing frames from netcdf/ASCII trajectory
framenum = 0;
if (netcdfLoad(&nc_ascii, argv[6]) == 0){
   netcdfLoad(&nc_ascii, argv[6]);
   printf("\n Processing NETCDF traj\n\n");
   while(netcdfGetNextFrame(&nc_ascii, xyz, grad)){
      energy=mme(xyz, grad, &framenum);
      energy=xmin(mme, &natm, xyz, grad, &energy, &min_energy, &xo);
      if(min_energy < atof(argv[4])){    //check convergence
        printf("     ----Convergence Satisfied---- \n\n");
        nmode(xyz, 3*natm, mme2, 0, 1, 0.0, 0.0, 0); //entropy calc
        framenum++;
      }
      else{
         printf(" \n    |----Not minimized properly!----|\n ");
         printf("   |---- Entropy not Calculated---|\n\n");
      }
   }
  netcdfClose(&nc_ascii);
}
else{
   printf("\n Processing ASCII traj\n\n");
   asciitraj= fopen(argv[6], "r");
   if(asciitraj==NULL){ //Check for the existence of the mdcrd
      printf("\n Unable to open mdcrd (%s) !\n\n", argv[6]);
      exit(1);
   }
   traj_complete=0;
   char *line = NULL;
   size_t len = 0;
   getline(&line, &len, asciitraj);
   for(i = 1;;i++){ // frame loop
      for(j = 0; j< prm->Nat3; j++){
         if(fscanf(asciitraj, "%lf", &xyz[j]) < 1){
            traj_complete = 1;
            break;
         }
      }     
      if (traj_complete) break;
         energy=mme(xyz, grad, &framenum);
         energy=xmin(mme, &natm, xyz, grad, &energy, &min_energy, &xo);
         if (min_energy< atof(argv[4])){ //check convergence
            printf("     ----Convergence Satisfied---- \n\n");
            nmode(xyz, 3*natm, mme2, 0, 1, 0.0, 0.0, 0); //calc entropy
            framenum++;
         }
      else{
        printf(" \n    |----Not minimized properly!----|\n ");
        printf("   |---- Entropy not Calculated---|\n\n");
      }
   }
}

}
