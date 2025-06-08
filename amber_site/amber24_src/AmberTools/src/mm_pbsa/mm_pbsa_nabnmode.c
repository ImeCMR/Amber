//
// NAB code for nmode calculations, following XMIN minimization
//
// Holger Gohlke
//   Last change: 07.02.2010
//

#include <nabc.h>

FILE* nabout;

int main( int argc, char* argv[] )
{

nabout = stdout;

// Prepare
XMIN_OPT_T xo;
MOLECULE_T *m;
PARMSTRUCT_T *prm;
int ier, natm;
double energy, grms;

// Check input
//printf("%s\n", argv[1]);
//printf("%s\n", argv[2]);
//printf("%s\n", argv[3]);
//printf("%s\n", argv[4]);
//printf("%s\n", argv[5]);

// Read molecule and prmtop
m = getpdb(argv[1], "");
prm = rdparm( argv[2] );
natm=prm->Natom;
double* xyz = malloc( prm->Nat3 * (sizeof(double)) );
double* grad = malloc( prm->Nat3 * (sizeof(double)) );

setxyz_from_mol(&m, NULL, (POINT_T*)xyz);

// Init MM calculations
mm_options(argv[3]);

// nothing frozen or constrained for now:
int* frozen = parseMaskString( "@ZZZ", prm, xyz, 2 );
int* constrained = parseMaskString( "@ZZZ", prm, xyz, 2 );

mme_init_sff( prm, frozen, constrained, NULL, NULL );


// Do XMIN minimization
xmin_opt_init( &xo );
xo.maxiter = atoi(argv[4]);
xo.grms_tol = atof(argv[5]);
xo.method = 3;
xo.numdiff = 1;
xo.m_lbfgs = 3;
xo.ls_method = 2;
xo.ls_maxiter = 20;
xo.ls_maxatmov = 0.15;
xo.print_level = 1;
int zero = 0;

energy = mme( xyz, grad, &zero);
energy = xmin( mme, &natm, xyz, grad, &energy, &grms, &xo );
if(grms > xo.grms_tol) return -3;

// Calc normal modes:
ier = nmode( xyz, 3*natm, mme2, 0, 0, 0.0, 0.0, 0);
return ier;

}
