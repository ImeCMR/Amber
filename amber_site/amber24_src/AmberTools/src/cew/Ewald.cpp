#include "Ewald.hpp"

namespace hfdf
{
  // amber mm=consistent conversions
   double const AMBERELE = 18.2223;
   double const BOHRS_TO_A = 0.529177249;
   double const CODATA08_A_TO_BOHRS = 1. / hfdf::BOHRS_TO_A;
   double const CODATA08_AU_TO_KCAL = hfdf::AMBERELE * hfdf::AMBERELE * hfdf::CODATA08_A_TO_BOHRS;

  // use EXTERN &gau-consistent conversions
  // double const AMBERELE = 18.2223;
  // double const BOHRS_TO_A = 1. / 1.889726132873;
  // double const CODATA08_A_TO_BOHRS = 1. / hfdf::BOHRS_TO_A;
  // double const CODATA08_AU_TO_KCAL = 6.2750946943E+02;

  double const AU_PER_AMBER_CHARGE = 1. / hfdf::AMBERELE;
  double const AU_PER_AMBER_ANGSTROM = 1. / hfdf::BOHRS_TO_A;
  double const AU_PER_AMBER_KCAL_PER_MOL = hfdf::AU_PER_AMBER_CHARGE * hfdf::AU_PER_AMBER_CHARGE / hfdf::AU_PER_AMBER_ANGSTROM;
  double const AU_PER_AMBER_MMPOT = 1. / hfdf::AU_PER_AMBER_ANGSTROM;
  double const AU_PER_AMBER_FORCE = hfdf::AU_PER_AMBER_KCAL_PER_MOL / hfdf::AU_PER_AMBER_ANGSTROM;

}





//int MYRANK;




#include <vector>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cassert>




namespace ccdl
{
  double const PI = 3.141592653589793238462643383279502884197;
  double const TWO_PI = 2. * PI;
  double const FOUR_PI = 4. * PI;

  double anint( double const x );

  double wrap( double const x, double const l );

  void wrap( double * R, double const * ucell, double const * recip );
  //void wrap2( double * R, double const * ucell, double const * recip );

  void bspline_one_pass( double * c, double const w, int const n );

  void bspline_diff( int const order, double const * array, double * diff );

  void bspline_eval( double const w, int const order, double * array );

  void bspline_eval( double const w, int const order, double * array, double * darray );

  void bspline_dct( int const N, int const order, double const * in, double * out, 
		    bool opt_infl = true );

  void bspline_hdct( int const N, int const order, double const * in, double * out,
		     bool opt_infl = true );

  void bspline_periodic( double x, int const nx, int const order, double * w );

  void bspline_periodic( double x, int const nx, int const order,double * w, int * gidx );

  void bspline_periodic( double x, int const nx, int const order,double * w,double * dw,int * gidx );

  void localpts_periodic( double x, int const nx, int const nlocal, int * gidx );

  void cross3( double const * a, double const * b, double * axb );

  double dot3( double const * a, double const * b );

  double gamma_sum_ratio( double frac_crd, int order );
}

///////////////////////////////////////////////////////////////////////////
// round to nearest integer; return as double
///////////////////////////////////////////////////////////////////////////
double ccdl::anint( double const x )
{
  return ((x)>0? std::floor((x)+0.5) : std::ceil((x)-0.5));
}
///////////////////////////////////////////////////////////////////////////
// simple Rij wrap for orthogonal unit cells
///////////////////////////////////////////////////////////////////////////
double ccdl::wrap( double const x, double const l )
{
  return x - l * ccdl::anint( x/l );
}
///////////////////////////////////////////////////////////////////////////
// Rij wrap for general lattice vectors
///////////////////////////////////////////////////////////////////////////
void ccdl::wrap( double * R, double const * ucell, double const * recip )
{
  double frac[3] = { ccdl::anint(R[0] * recip[0] + R[1] * recip[1] + R[2] * recip[2]),
		     ccdl::anint(R[0] * recip[3] + R[1] * recip[4] + R[2] * recip[5]),
		     ccdl::anint(R[0] * recip[6] + R[1] * recip[7] + R[2] * recip[8]) };
  R[0] -= frac[0] * ucell[0] + frac[1] * ucell[3] + frac[2] * ucell[6];
  R[1] -= frac[0] * ucell[1] + frac[1] * ucell[4] + frac[2] * ucell[7];
  R[2] -= frac[0] * ucell[2] + frac[1] * ucell[5] + frac[2] * ucell[8];
}
// void ccdl::wrap2( double * R, double const * ucell, double const * recip )
// {
//   double frac[3] = { R[0] * recip[0] + R[1] * recip[1] + R[2] * recip[2],
// 		     R[0] * recip[3] + R[1] * recip[4] + R[2] * recip[5],
// 		     R[0] * recip[6] + R[1] * recip[7] + R[2] * recip[8]};
//   frac[0] -= ccdl::anint(frac[0]);
//   frac[1] -= ccdl::anint(frac[1]);
//   frac[2] -= ccdl::anint(frac[2]);
//   R[0] = frac[0] * ucell[0] + frac[1] * ucell[3] + frac[2] * ucell[6];
//   R[1] = frac[0] * ucell[1] + frac[1] * ucell[4] + frac[2] * ucell[7];
//   R[2] = frac[0] * ucell[2] + frac[1] * ucell[5] + frac[2] * ucell[8];
// }

///////////////////////////////////////////////////////////////////////////
// bspline utility function (internal helper function)
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_one_pass( double * c, double const w, int const n )
{
  int nm1 = n-1;
  double const div = 1. / nm1;
  c[nm1] = div * w * c[nm1 - 1];
  for ( int j=1; j<nm1; ++j )
    c[nm1 - j] = div * ((w + j) * c[nm1 - j - 1] + (n - j - w) * c[nm1 - j]);
  c[0] = div * (1 - w) * c[0];
}
///////////////////////////////////////////////////////////////////////////
// bspline derivative utility function (internal helper function)
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_diff( int const order, double const * array, double * diff )
{
  assert( order > 1 );
  int const nm1 = order-1;
  diff[0] = -array[0];
  for ( int j=1; j<nm1; ++j )
    diff[j] = array[j-1] - array[j];
  diff[nm1] = array[nm1-1];
}
///////////////////////////////////////////////////////////////////////////
// 1d bspline evaluation.  
// w is a number between 0 and 1.
// If order is even, then w has a value of 0.5 between 2 grid points.
// If order is odd, then w has a value of 0.5 on a grid point.
// The first value in array is the left-most-point
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_eval( double const w, int const order, double * array )
{
  array[0] = 1. - w;
  array[1] = w;
  if (order > 2)
    {
      // One pass to order 3:
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];
      if (order > 3)
        {
          // One pass to order 4:         
          double const div = 1./3.;
          array[3] = div * w * array[2];
          array[2] = div * ((w + 1.) * array[1] + (3. - w) * array[2]);
          array[1] = div * ((w + 2.) * array[0] + (2. - w) * array[1]);
          array[0] = div * (1. - w) * array[0];
          // and the rest
          for ( int k = 5; k < order+1; ++k )
            ccdl::bspline_one_pass(array, w, k);
        };
    };
}
///////////////////////////////////////////////////////////////////////////
// 1d bspline evaluation and first derivative
// w is a number between 0 and 1.
// If order is even, then w has a value of 0.5 between 2 grid points.
// If order is odd, then w has a value of 0.5 on a grid point.
// The first value in array is the left-most-point
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_eval( double const w, int const order, double * array, double * darray )
{
  assert( order > 2 );
  double const div = 1./3.;

  array[0] = 1. - w;
  array[1] = w;

  if (order == 4)
    {
      // One pass to order 3:
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];
      
      darray[0] = -array[0];
      darray[1] = array[0]-array[1];
      darray[2] = array[1]-array[2];
      darray[3] = array[2];

      // One pass to order 4:     
      array[3] = div * w * array[2];
      array[2] = div * ((w + 1.) * array[1] + (3. - w) * array[2]);
      array[1] = div * ((w + 2.) * array[0] + (2. - w) * array[1]);
      array[0] = div * (1. - w) * array[0];
      
    }
  else if ( order > 4 )
    {
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];

      array[3] = div * w * array[2];
      array[2] = div * ((w + 1.) * array[1] + (3. - w) * array[2]);
      array[1] = div * ((w + 2.) * array[0] + (2. - w) * array[1]);
      array[0] = div * (1. - w) * array[0];

      // and the rest
      for ( int k = 5; k < order; ++k ) // don't do k==order
        ccdl::bspline_one_pass(array, w, k);

      ccdl::bspline_diff(order,array,darray);

      // One more recursion: // do the k==order
      ccdl::bspline_one_pass(array, w, order);

    }
  else // order == 3
    {
      darray[0] = -array[0];
      darray[1] = array[0]-array[1];
      darray[2] = array[1];

      // One pass to order 3:
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];
    };
}
///////////////////////////////////////////////////////////////////////////
// Evaluates a 1d bspline weights for a function with fractional crd x
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_periodic
( double x, 
  int const nx,
  int const order,
  double * w )
{
  x = x*nx - (order%2)*0.5;
  ccdl::bspline_eval(x-std::floor(x),order,w);
}
///////////////////////////////////////////////////////////////////////////
// Evaluates a 1d bspline weights for a function with fractional crd x
// and also returns the grid indices, considering periodicity
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_periodic
( double x, 
  int const nx,
  int const order,
  double * w,
  int * gidx )
{
  x = x*nx - (order%2)*0.5;
  int ilo = std::floor(x);
  x -= ilo;
  // the first modulus takes us to [-nx+1,nx-1]
  // then we add nx, so we are in the range [0,2*nx-1]
  // the second modulus then puts us in [0,nx-1]
  // Each ilo in the loop uses a modulus, so we don't need
  // the second modulus on the first calculation of ilo
  //ilo = ((ilo+1-order/2)%nx+nx)%nx;
  ilo   =  (ilo+1-order/2)%nx+nx;//%nx;
  for ( int b=0; b<order; ++b, ++ilo )
    gidx[b] = ilo%nx;
  ///////////////////////////
  ccdl::bspline_eval(x,order,w);
}
///////////////////////////////////////////////////////////////////////////
// Evaluates a 1d bspline weights and first derivatives 
// for a function with fractional crd x
// and also returns the grid indices, considering periodicity
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_periodic
( double x, 
  int const nx,
  int const order,
  double * w,
  double * dw,
  int * gidx )
{
  x = x*nx - (order%2)*0.5;
  int ilo = std::floor(x);
  ccdl::bspline_eval(x-ilo,order,w,dw);
  // the first modulus takes us to [-nx+1,nx-1]
  // then we add nx, so we are in the range [0,2*nx-1]
  // the second modulus then puts us in [0,nx-1]
  // Each ilo in the loop uses a modulus, so we don't need
  // the second modulus on the first calculation of ilo
  //ilo = ((ilo+1-order/2)%nx+nx)%nx;
  ilo   =  (ilo+1-order/2)%nx+nx;//%nx;
  x = nx;
  for ( int b=0; b<order; ++b, ++ilo )
    {
      gidx[b] = ilo%nx;
      dw[b] *= x;
    };
  ///////////////////////////
}

///////////////////////////////////////////////////////////////////////////
// Returns an array that index the "order" grid points nearest
// to an atom with fractional crd x, considering periodicity
///////////////////////////////////////////////////////////////////////////
void ccdl::localpts_periodic
( double x, 
  int const nx,
  int const order,
  int * gidx )
{
  x = x*nx - (order%2)*0.5;
  int ilo = std::floor(x);
  //x -= ilo; // unused in this routine
  //ilo = ((ilo+1-order/2)%nx+nx)%nx; // dont need 2nd % because its in the loop
  ilo = ((ilo+1-order/2)%nx+nx); //%nx;
  for ( int b=0; b<order; ++b, ++ilo )
    gidx[b] = ilo%nx;
}
///////////////////////////////////////////////////////////////////////////
// computes the DFT coefficients of a 1d Bspline 
// (all frequencies - the outer loops in FFTW)
// If opt_infl = true (default,optional), then the dft coefficients
// are "monkey'd with" to make them amber-compatible.
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_dct( int const N, int const order, double const * in, double * out, 
			bool opt_infl )
{
  std::fill( out, out + N, 0. );
  
  int const o = (order-1)/2;
  int const Np = (order)/2+1;
  int const Nh = N-o;
  
  double const tpion = ccdl::TWO_PI/N;
  for ( int m=0; m<N/2+1; ++m )
    {
      double const k = m * tpion;
      for ( int n=0; n<Np; ++n )
	out[m] += in[n+o] * std::cos(k*n);
      for ( int i=0; i<o; ++i )
	out[m] += in[i] * std::cos(k*(Nh+i));
    };
  if ( opt_infl )
    for ( int m=0; m<N/2+1; ++m )
      out[m] *= ccdl::gamma_sum_ratio( (m < N/2+1 ? m : N-m)/(double)(N), order );
  for ( int m=N/2+1; m<N; ++m )
    out[m]=out[N-m];

}
///////////////////////////////////////////////////////////////////////////
// computes the DFT coefficients of a 1d Bspline 
// (only half the frequencies - the inner loop in FFTW)
// If opt_infl = true (default,optional), then the dft coefficients
// are "monkey'd with" to make them amber-compatible.
///////////////////////////////////////////////////////////////////////////
void ccdl::bspline_hdct( int const N, int const order, double const * in, double * out, 
			 bool opt_infl )
{
  std::fill( out, out + N/2+1, 0. );

  int const Nmax = N/2+1;
  int const o = (order-1)/2;
  int const Np = (order)/2+1;
  int const Nh = N-o;
  
  double const tpion = ccdl::TWO_PI/N;
  for ( int m=0; m<Nmax; ++m )
    {
      double const k = m * tpion;
      for ( int n=0; n<Np; ++n )
	out[m] += in[n+o] * std::cos(k*n);
      for ( int i=0; i<o; ++i )
	out[m] += in[i] * std::cos(k*(Nh+i));
    };
  if ( opt_infl )
    for ( int m=0; m<Nmax; ++m )
      out[m] *= ccdl::gamma_sum_ratio( m/(double)(N), order );
}
///////////////////////////////////////////////////////////////////////////
// computes the cross product of two vectors of length 3
///////////////////////////////////////////////////////////////////////////
void ccdl::cross3( double const * a, double const * b, double * axb )
{
  axb[0] = a[1]*b[2] - a[2]*b[1];
  axb[1] = b[0]*a[2] - a[0]*b[2];
  axb[2] = a[0]*b[1] - a[1]*b[0];
}
///////////////////////////////////////////////////////////////////////////
// computes the dot product of two vectors of length 3
///////////////////////////////////////////////////////////////////////////
double ccdl::dot3( double const * a, double const * b )
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
///////////////////////////////////////////////////////////////////////////
// computes a scale factor for the b-spline dft coefficient
// This scale factor isn't necessary, but some people use it.
// (The default in amber is to use it -- and the option to turn it off
//  is an undocumented secret).
// It is included here for compatibility purposes.
//   frac_crd is m/N, where m \in (-N/2 ... N/2)
///////////////////////////////////////////////////////////////////////////
double ccdl::gamma_sum_ratio( double frac_crd, int order )
{
  double gn  = 1.;
  double g2n = 1.;
  frac_crd *= ccdl::PI;
  for ( int i=1; i<51; ++i )
    {
      double c = std::pow( frac_crd / ( frac_crd + ccdl::PI*i ) , order );
      double s = std::pow( frac_crd / ( frac_crd - ccdl::PI*i ) , order );
      gn  += c + s;
      g2n += c*c + s*s;
    }
  return g2n/gn;
}










cew::Ewald::Ewald
#ifdef WITH_MPI
( int const comm,
  int const * p_natom,
  int const * p_nquant,
  int const * p_nquant_nlink,
  int const * p_iqmatoms,
  double const * p_charges,
  double const * p_qm_charges,
  double const * p_nucq,
  int const * p_nfft1,
  int const * p_nfft2,
  int const * p_nfft3,
  int const * p_order,
  double const * p_ew_coeff,
  double const * p_cut )
#else
( int const * p_natom,
  int const * p_nquant,
  int const * p_nquant_nlink,
  int const * p_iqmatoms,
  double const * p_charges,
  double const * p_qm_charges,
  double const * p_nucq,
  int const * p_nfft1,
  int const * p_nfft2,
  int const * p_nfft3,
  int const * p_order,
  double const * p_ew_coeff,
  double const * p_cut )
#endif
  : 
#ifdef WITH_MPI
  CommGroup(MPI_Comm_f2c(comm)),
  mpiSize(1), mpiRank(0),
#endif
  nat( *p_natom ),
  nqm( *p_nquant_nlink ),
  nquant( *p_nquant ),
  iqmatoms( p_iqmatoms, p_iqmatoms + *p_nquant_nlink ),
  charges( p_charges, p_charges + *p_natom ),
  qm_charges( p_qm_charges, p_qm_charges + *p_nquant_nlink ),
  qmq( p_qm_charges, p_qm_charges + *p_nquant_nlink ),
  nucq( p_nucq, p_nucq + *p_nquant_nlink ),
  //
  msgbuf( 1 + 9 + (3*nat+1) + (3*nat+1) + (3*nat), 0 ),
  ucell_changed( msgbuf.data() ),
  ucell( msgbuf.data() + 1 ),
  crd( ucell + 9 ),
  frc( crd + (3*nat) + 1 ),
  imgcrd( frc + (3*nat) + 1 ),
  //
  qmc( 3*nqm, 0. ),
  qmp( nqm, 0. ),
  //
  opt_infl( true ),
  recip(9,0.),
  n1( *p_nfft1 ), 
  n2( *p_nfft2 ),
  n3( *p_nfft3 ),
  nt( n1*n2*n3 ),
  nk3( n3/2+1 ),
  nk( n1*n2*nk3 ),
  order( *p_order ),
  beta( *p_ew_coeff ),
  zeta( beta*beta ),
  qmcut( *p_cut ),
  V(0.),
  w1(order,0.), w2(order,0.), w3(order,0.),
  bf1(n1,0.), bf2(n2,0.), bf3(nk3,0.),
  kernel(nk,0.),
  localfrc( 3*nat + 1, 0. ),
  localene( 0. ),
  ChargeCorrection( 0. ),
  internalene( 0. ),
  tmp_grid( nt, 0. )
{

  //MYRANK=0;
#ifdef WITH_MPI
  MPI_Comm_size( CommGroup, &mpiSize );
  MPI_Comm_rank( CommGroup, &mpiRank );
  //MYRANK=mpiRank;
#endif

  
  
  grid    = fftw_alloc_real( nt );

  dftcoef = reinterpret_cast< std::complex<double> * >( fftw_alloc_complex( nk ) );

  fplan   = new fftw_plan
    ( fftw_plan_dft_r2c_3d
      ( n1, n2, n3, 
        grid, 
        reinterpret_cast< fftw_complex *>( dftcoef ), 
        FFTW_ESTIMATE ) );

  rplan   = new fftw_plan
    ( fftw_plan_dft_c2r_3d
      ( n1, n2, n3, 
        reinterpret_cast< fftw_complex *>( dftcoef ), 
        grid, 
        FFTW_ESTIMATE ) );
}






cew::Ewald::~Ewald()
{
  fftw_destroy_plan( *fplan );
  delete fplan;
  fplan = NULL;

  fftw_destroy_plan( *rplan );
  delete rplan;
  rplan = NULL;

  fftw_free(reinterpret_cast< fftw_complex * &>(dftcoef));
  dftcoef = NULL;

  fftw_free( grid );
  grid = NULL;

  //fftw_cleanup();
}




void cew::Ewald::SetCrd
( double const * sander__coords_after_lnk_fix,
  double const * nblist__ucell )
{
  *ucell_changed = 0;

  for ( int i=0; i<9; ++i )
    {
      double t = nblist__ucell[i] * hfdf::AU_PER_AMBER_ANGSTROM;
      if ( std::abs(ucell[i]-t) > 1.e-6 ) 
	*ucell_changed = 1;
      ucell[i] = t;
    };
  int ncrd = 3*nat;
  for ( int i=0; i<ncrd; ++i )
    crd[i] = sander__coords_after_lnk_fix[i] * hfdf::AU_PER_AMBER_ANGSTROM;


#ifdef WITH_MPI
  MPI_Bcast( msgbuf.data(), 1 + 9 + 3*nat, MPI_DOUBLE, 0, CommGroup );
#endif

  if ( ( *ucell_changed > 0.5 ) or ( V < 0.001 ) )
    UpdateKernel();
  ImageCrdAroundQM();
}



#ifdef WITH_MPI
void cew::Ewald::SetCrd()
{
  MPI_Bcast( msgbuf.data(), 1 + 9 + 3*nat, MPI_DOUBLE, 0, CommGroup );
  if ( ( *ucell_changed > 0.5 ) or ( V < 0.001 ) )
    UpdateKernel();
  ImageCrdAroundQM();
}
#endif





#include <cstdio>



void cew::Ewald::prescf
( double const * sander__coords_after_lnk_fix,
  double const * nblist__ucell,
  double * p_qmcrd,
  int * p_nmm, double * p_mmcrd, double * p_mmcharges )
{

  std::fill( grid, grid+nt, 0. );
  std::fill( dftcoef, dftcoef+nk, std::complex<double>(0.,0.) );

  
  internalene = 0.;
  SetCrd( sander__coords_after_lnk_fix, nblist__ucell );
  
  for ( std::size_t i=0, n=qmc.size(); i<n; ++i )
    p_qmcrd[i] = qmc[i] / hfdf::AU_PER_AMBER_ANGSTROM;
  
  for ( std::size_t i=0, n=mmq.size(); i<n; ++i )
    p_mmcharges[i] = mmq[i];
  
  for ( std::size_t i=0, n=mmc.size(); i<n; ++i )
    p_mmcrd[i] = mmc[i] / hfdf::AU_PER_AMBER_ANGSTROM;
  
  *p_nmm = immatoms.size();
  
  int const ncrd = 3*nat;
  int const nmm  = immatoms.size();
  std::fill( frc, frc + ncrd + 1, 0. );
  std::fill( localfrc.data(), localfrc.data() + ncrd + 1, 0. );

  /*
  std::printf("qmcrd\n");
  for ( int i=0; i<nqm; ++i )
    {
      std::printf("%7i %7i ",i,iqmatoms[i]);
      for ( int k=0; k<3; ++k )
	std::printf("%12.5f",p_qmcrd[k+i*3]);
      std::printf("%12.5f",qm_charges[i]);
      std::printf("\n");
    };
  std::printf("mmcrd\n");
  for ( int i=0; i<*p_nmm; ++i )
    {
      std::printf("%7i %7i ",i,immatoms[i]);
      for ( int k=0; k<3; ++k )
	std::printf("%12.5f",p_mmcrd[k+i*3]);
      std::printf("%12.5f",p_mmcharges[i]);
      std::printf("\n");
    };
  */

  
  double localene = 0.;

  std::copy( qm_charges.data(), qm_charges.data() + nqm, qmq.data() );

  double g[3] = { 0.,0.,0. };

  
  // -0.5 <qm'|jn|qm'>

#ifdef WITH_MPI
  if ( ! mpiRank )
#endif
    {
      CptPlaneWavePotControlPts( nqm, qmc.data(), qmq.data() );
      for ( int i=0; i<nqm; ++i )
	{
	  int a = iqmatoms[i];
	  double p = CptPlaneWavePotAndGrd( qmc.data() + 3*i, g );
	  localene -= 0.5 * p * qmq[i];
	  localfrc[0+a*3] -= qmq[i] * g[0];
	  localfrc[1+a*3] -= qmq[i] * g[1];
	  localfrc[2+a*3] -= qmq[i] * g[2];
	};
    }

  for ( int i=0; i<nqm; ++i )
    {
      int a = iqmatoms[i];
      charges[a] = qm_charges[i];
    };

  /*
  std::printf("prescf localfrc %20.10e\n",localene);
  for ( int i=0; i<nqm; ++i )
    {
      std::printf("%7i %7i ",i,iqmatoms[i]);
      for ( int k=0; k<3; ++k )
	std::printf("%12.5f",localfrc[k+i*3]);
      std::printf("\n");
    }
  */
  
  // <Rt|jn|mm+qm'>

  CptPlaneWavePotControlPts( nat, imgcrd, charges.data() );

  

  
  std::fill( qmp.data(), qmp.data() + nqm, 0. );
  std::fill( tmp_grid.data(), tmp_grid.data() + nt, 0. );
  
  internalene = localene;
}




void cew::Ewald::postscf
( double const * qmgrd,
  double const * mmgrd,
  double * eamber,
  double * frcamber )
{
  int const ncrd = 3*nat;
  int const nmm  = immatoms.size();

  double ene = 0.;
  double localene = internalene;
  double g[3] = { 0.,0.,0. };


  // everything labeled frc is actually a gradient EXCEPT for frcamber
  
#ifdef WITH_MPI
  if ( ! mpiRank )
#endif
    {
      for ( int i=0; i<nqm; ++i )
	{
	  int a = iqmatoms[i];
	  //std::printf("qmfrc au %4i %12.5f %12.5f %12.5f\n",i+1,qmgrd[0+i*3] * hfdf::AU_PER_AMBER_FORCE,qmgrd[1+i*3] * hfdf::AU_PER_AMBER_FORCE,qmgrd[2+i*3] * hfdf::AU_PER_AMBER_FORCE);

	  localfrc[0+a*3] += qmgrd[0+i*3] * hfdf::AU_PER_AMBER_FORCE;
	  localfrc[1+a*3] += qmgrd[1+i*3] * hfdf::AU_PER_AMBER_FORCE;
	  localfrc[2+a*3] += qmgrd[2+i*3] * hfdf::AU_PER_AMBER_FORCE;
	};
      for ( int i=0; i<nmm; ++i )
	{
	  int a = immatoms[i];
	  //std::printf("mmfrc au %4i %12.5f %12.5f %12.5f\n",i+1,mmgrd[0+(i)*3] * hfdf::AU_PER_AMBER_FORCE,mmgrd[1+(i)*3] * hfdf::AU_PER_AMBER_FORCE,mmgrd[2+(i)*3] * hfdf::AU_PER_AMBER_FORCE);

	  localfrc[0+a*3] += mmgrd[0+(i)*3] * hfdf::AU_PER_AMBER_FORCE;
	  localfrc[1+a*3] += mmgrd[1+(i)*3] * hfdf::AU_PER_AMBER_FORCE;
	  localfrc[2+a*3] += mmgrd[2+(i)*3] * hfdf::AU_PER_AMBER_FORCE;
	};
    }

  
  std::copy( tmp_grid.data(), tmp_grid.data() + nt, grid );
  std::copy( nucq.data(), nucq.data() + nqm, qmq.data() );

		
#ifdef WITH_MPI
  for ( int i=0; i<nqm; ++i ) qmq[i] /= mpiSize;
#endif
  
  AccDensityControlPts( nqm, qmc.data(), qmq.data() );
  CptPlaneWavePotControlPts();

  for ( int i=0; i<nat; ++i )
    {
      CptPlaneWaveGrd( imgcrd + i*3, g );
      // if ( i == 0 )
      // 	{
      // 	  std::printf("grad %4i%15.10f%15.10f%15.10f%15.10f%15.10f%15.10f%15.10f\n",
      // 		      i,charges[i],imgcrd[0+i*3],imgcrd[1+i*3],imgcrd[2+i*3],
      // 		      charges[i] *g[0]/hfdf::AU_PER_AMBER_FORCE,
      // 		      charges[i] *g[1]/hfdf::AU_PER_AMBER_FORCE,
      // 		      charges[i] *g[2]/hfdf::AU_PER_AMBER_FORCE);
      // 	}
	
      frc[0+i*3] += charges[i] * g[0];
      frc[1+i*3] += charges[i] * g[1];
      frc[2+i*3] += charges[i] * g[2];
    }
  
#ifdef WITH_MPI
  std::copy( frc, frc + ncrd, crd );
  MPI_Reduce( crd, frc, ncrd, MPI_DOUBLE, MPI_SUM, 0, CommGroup );
  if ( mpiRank )
    {
      std::fill( frc, frc + ncrd, 0. );
      ene = 0.;
      localene = 0.;
    }
  else
#endif
    {
      for ( int a=0; a < ncrd; ++a )
	frc[a] += localfrc[a];
      ene += localene;
    }

  *eamber += ene / hfdf::AU_PER_AMBER_KCAL_PER_MOL;
  
  // change grd to a frc and accumulate
  for ( int a=0; a < ncrd; ++a )
    frcamber[a] -= frc[a] / hfdf::AU_PER_AMBER_FORCE;

  
  for ( int i=0; i<nqm; ++i )
    {
      int a = iqmatoms[i];
      charges[a] = 0.;
    };
  
}


  

double cew::Ewald::EnergyAndGradients()
{
  int const ncrd = 3*nat;
  int const nmm  = immatoms.size();
  std::fill( frc, frc + ncrd + 1, 0. );
  std::fill( localfrc.data(), localfrc.data() + ncrd + 1, 0. );
  double ene = 0.;
  double localene = 0.;
  double g[3] = { 0.,0.,0. };

  //std::printf("%zu %zu\n",qmq.size(),options->qm_charges.size());
  std::copy( qm_charges.data(), qm_charges.data() + nqm, qmq.data() );


  
  // -0.5 <qm'|jn|qm'>

#ifdef WITH_MPI
  if ( ! mpiRank )
#endif
    {
      CptPlaneWavePotControlPts( nqm, qmc.data(), qmq.data() );
      for ( int i=0; i<nqm; ++i )
	{
	  int a = iqmatoms[i];
	  double p = CptPlaneWavePotAndGrd( qmc.data() + 3*i, g );
	  localene -= 0.5 * p * qmq[i];
	  localfrc[0+a*3] -= qmq[i] * g[0];
	  localfrc[1+a*3] -= qmq[i] * g[1];
	  localfrc[2+a*3] -= qmq[i] * g[2];
	};
    }

  for ( int i=0; i<nqm; ++i )
    {
      int a = iqmatoms[i];
      charges[a] = qm_charges[i];
    };
  
  // <Rt|jn|mm+qm'>

  CptPlaneWavePotControlPts( nat, imgcrd, charges.data() );

  std::fill( qmp.data(), qmp.data() + nqm, 0. );

  //
  // end pre-scf
  //
  
  //scf->SetAtomCrd( qmc.data(), qmp.data(), nmm, mmc.data(), mmq.data() );
  //ene = scf->ScfProcedure();

  
  // nucq pointer is only valid after SetAtomCrd
  //double const * nucq = scf->GetNuclearCharges();



  //
  // begin post-scf
  // requires: grd, nucq
  //

  // this includes the left-hand b-spline derivatives
  //std::fill( crd, crd + 3*(nqm+nmm), 0. );
  //scf->CptGrd();

  
#ifdef WITH_MPI
  if ( ! mpiRank )
#endif
    {
      //scf->GetGrd( nqm+nmm, crd ); 
      
      for ( int i=0; i<nqm; ++i )
	{
	  int a = iqmatoms[i];
	  localfrc[0+a*3] += crd[0+i*3];
	  localfrc[1+a*3] += crd[1+i*3];
	  localfrc[2+a*3] += crd[2+i*3];
	};
      for ( int i=0; i<nmm; ++i )
	{
	  int a = immatoms[i];
	  localfrc[0+a*3] += crd[0+(i+nqm)*3];
	  localfrc[1+a*3] += crd[1+(i+nqm)*3];
	  localfrc[2+a*3] += crd[2+(i+nqm)*3];
	};
    }

 
  /////////////////////////////////////////////////////////////////////////////

  
  std::copy( tmp_grid.data(), tmp_grid.data() + nt, grid );
  std::copy( nucq.data(), nucq.data() + nqm, qmq.data() );
#ifdef WITH_MPI
  for ( int i=0; i<nqm; ++i ) qmq[i] /= mpiSize;
#endif
  AccDensityControlPts( nqm, qmc.data(), qmq.data() );
  CptPlaneWavePotControlPts();
  for ( int i=0; i<nat; ++i )
    {
      CptPlaneWaveGrd( imgcrd + i*3, g );
      frc[0+i*3] += charges[i] * g[0];
      frc[1+i*3] += charges[i] * g[1];
      frc[2+i*3] += charges[i] * g[2];
    }
#ifdef WITH_MPI
  std::copy( frc, frc + ncrd, crd );
  MPI_Reduce( crd, frc, ncrd, MPI_DOUBLE, MPI_SUM, 0, CommGroup );
  if ( mpiRank )
    {
      std::fill( frc, frc + ncrd, 0. );
      ene = 0.;
      localene = 0.;
    }
  else
#endif
    {
      for ( int a=0; a < ncrd; ++a )
	frc[a] += localfrc[a];
      ene += localene;
    }


  for ( int i=0; i<nqm; ++i )
    {
      int a = iqmatoms[i];
      charges[a] = 0.;
    };
  
  ///////////////////////////////////////////////////////////////////////////////

  return ene;
}



void cew::Ewald::UpdateKernel()
{
  // A column of ucell is a real space lattice vector, i.e.
  // a_1 = { ucell[0], ucell[1], ucell[2] };
  //
  // The reciprocal space lattice vectors are the columns of recip

  std::fill( recip.data(), recip.data() + 9, 0. );
  ccdl::cross3( ucell + 3, ucell + 6, recip.data() + 0 );
  ccdl::cross3( ucell + 6, ucell + 0, recip.data() + 3 );
  ccdl::cross3( ucell + 0, ucell + 3, recip.data() + 6 );
  V = ccdl::dot3( ucell + 0, recip.data() + 0 );
  for ( int i=0; i<9; ++i ) recip[i] /= V;


  ccdl::bspline_eval( ( order%2 == 0 ? 0.0 : 0.5 ), order, w1.data() );
  ccdl::bspline_dct(  n1, order, w1.data(), bf1.data(), opt_infl );
  if ( n2 == n1 ) 
    bf2=bf1;
  else
    ccdl::bspline_dct(  n2, order, w1.data(), bf2.data(), opt_infl );
  if ( n3 == n1 )
    std::copy( bf1.data(), bf1.data() + nk3, bf3.data() );
  else if ( n3 == n2 )
    std::copy( bf2.data(), bf2.data() + nk3, bf3.data() );
  else
    ccdl::bspline_hdct( n3, order, w1.data(), bf3.data(), opt_infl );

  /*
  std::printf("zeta %20.10e %3i %3i %3i %20.10e\n",zeta,n1,n2,n3,V);
  */
  
  double t = -0.25 / zeta;
  int ik = 0;
  for ( int i1 = 0; i1 < n1; ++i1 )
    {
      double m1 = ( i1 < n1/2+1 ? i1 : -(n1-i1) );
      double k1[3] = { m1 * recip[0],
                       m1 * recip[1],
                       m1 * recip[2] };
      for ( int i2 = 0; i2 < n2; ++i2 )
        {
          double m2 = ( i2 < n2/2+1 ? i2 : -(n2-i2) );
          double b12 = bf1[i1] * bf2[i2];
          double k12[3] = { k1[0] + m2 * recip[3],
                            k1[1] + m2 * recip[4],
                            k1[2] + m2 * recip[5] };
          for ( int i3 = 0; i3 < nk3; ++i3, ++ik )
            {
              // k = 2pi ( k1 * a1* + k2 * a2* + k3 * a3* )
              // where k1,k2,k3 are integers and ai* is the i'th reciprocal
              // space lattice vector
              double kv[3] = { ccdl::TWO_PI * ( k12[0] + i3 * recip[6] ),
                               ccdl::TWO_PI * ( k12[1] + i3 * recip[7] ),
                               ccdl::TWO_PI * ( k12[2] + i3 * recip[8] )};
              double k2 = kv[0]*kv[0] + kv[1]*kv[1] + kv[2]*kv[2];
              // bk is the DFT coefficient of the bspline weight
              double bk = b12 * bf3[i3];
              // gk is the Fourier coefficient of the gaussian <k|g> = e^{-k2/4zeta}
              double gk = std::exp(k2*t);
              // The electrostatic potential of <r|k> is 4pi/k2 <r|k>
              // The extra factor of 1/V is the <k|k'>^{-1} part of the fit
              kernel[ik] = (ccdl::FOUR_PI/k2)*gk/(bk*bk*V); 
            }
        }
    }
  // remove the k==0 term
  kernel[0] = 0.;
  
}





void qmmm_image_backend
( int const nat,
  double * crd,
  double const * ucell,
  double const * recip,
  double const qmcut,
  int const nquant,
  int const nquant_nlink,
  int const * iqmatoms, // 0 .. nquant-1 .. nquant_nlink-1
  std::vector<int> & mmatoms ) // 0 .. nat-1
{
  //double const A = 0.529177249;
  double const qmcut2 = qmcut*qmcut;

  mmatoms.resize(0);
  
  double avg[3] = { 0.,0.,0. };

  std::vector<double> Rsums(nquant,0.);
  for ( int i=0; i < nquant; ++i )
    {
      int const a = iqmatoms[i];
      for ( int j=i+1; j<nquant; ++j )
  	{
  	  int const b = iqmatoms[j];
  	  double R[3] = { crd[0+a*3] - crd[0+b*3],
  			  crd[1+a*3] - crd[1+b*3],
  			  crd[2+a*3] - crd[2+b*3] };
  	  ccdl::wrap( R, ucell, recip );
  	  double Rmag=std::sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
	  Rsums[i] += Rmag;
	  Rsums[j] += Rmag;
  	};
    };

  double rmin = 1.e+10;
  int imin = 0;
  for ( int i=0; i < nquant; ++i )
    {
      if ( Rsums[i] < rmin )
	{
	  imin = i;
	  rmin = Rsums[i];
	};
    }
  for ( int k=0; k<3; ++k )
    avg[k] = crd[k+iqmatoms[imin]*3];


  ccdl::wrap( avg, ucell, recip );

  
  
  double box[6] = { 1.e+6, -1.e+6,  1.e+6, -1.e+6,  1.e+6, -1.e+6 };


  for ( int i=0; i < nquant; ++i )
    {
      int const a = iqmatoms[i];
      double R[3] = { crd[0+a*3] - avg[0],
		      crd[1+a*3] - avg[1],
		      crd[2+a*3] - avg[2] };
      ccdl::wrap( R, ucell, recip );
      box[0] = std::min( box[0], R[0] );
      box[1] = std::max( box[1], R[0] );
      box[2] = std::min( box[2], R[1] );
      box[3] = std::max( box[3], R[1] );
      box[4] = std::min( box[4], R[2] );
      box[5] = std::max( box[5], R[2] );
    };


  box[0] += avg[0]-qmcut;
  box[1] += avg[0]+qmcut;
  box[2] += avg[1]-qmcut;
  box[3] += avg[1]+qmcut;
  box[4] += avg[2]-qmcut;
  box[5] += avg[2]+qmcut;
  
  avg[0] = 0.5 * ( box[0]+box[1] );
  avg[1] = 0.5 * ( box[2]+box[3] );
  avg[2] = 0.5 * ( box[4]+box[5] );



  box[0] -= avg[0];
  box[1] -= avg[0];
  box[2] -= avg[1];
  box[3] -= avg[1];
  box[4] -= avg[2];
  box[5] -= avg[2];
  

  int const * iqmatoms_end = iqmatoms + nquant_nlink;
  int const * b_end = iqmatoms + nquant;
  for ( int a = 0; a < nat; ++a )
    {
      crd[0+a*3] -= avg[0];
      crd[1+a*3] -= avg[1];
      crd[2+a*3] -= avg[2];
      ccdl::wrap( crd+a*3, ucell, recip );
    };
  for ( int a=0; a < nat; ++a )
    {
      if ( crd[0+a*3] > box[0] and crd[0+a*3] < box[1] and
      	   crd[1+a*3] > box[2] and crd[1+a*3] < box[3] and
      	   crd[2+a*3] > box[4] and crd[2+a*3] < box[5] )
	{
	  int const * p = std::find( iqmatoms, iqmatoms_end, a );
	  if ( p == iqmatoms_end ) // "a" is not a qm nor link atom
	    {
	      for ( int const * b=iqmatoms; b < b_end; ++b )
		{
		  double R[3] = 
		    { crd[0+a*3]-crd[0+*b*3],
		      crd[1+a*3]-crd[1+*b*3],
		      crd[2+a*3]-crd[2+*b*3] };

		  ccdl::wrap( R, ucell, recip );

		  double const r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

		  if ( r2 < qmcut2 )
		    {
		      mmatoms.push_back( a );
		      break;
		    };
		};
	    };
	};
      
    };



}







void cew::Ewald::ImageCrdAroundQM()
{
  int ncrd = 3 * nat;
  std::copy( crd, crd + ncrd, imgcrd );
  qmmm_image_backend
    ( nat, imgcrd, 
      ucell, recip.data(), qmcut, 
      nquant, nqm, // nqm = nquant_nlink
      iqmatoms.data(), immatoms );

  
  {
    int i = iqmatoms[0];
    double d[3] = { crd[0+i*3] - imgcrd[0+i*3],
		    crd[1+i*3] - imgcrd[1+i*3],
		    crd[2+i*3] - imgcrd[2+i*3] };
    
    for ( int a=0; a<nat; ++a )
      for ( int k=0; k<3; ++k )
	imgcrd[k+a*3] += d[k];
  }
  

  int nmm = immatoms.size();
  mmc.resize( 3*nmm );
  mmq.resize( nmm );

  for ( int i=0; i<nqm; ++i )
    {
      int a = iqmatoms[i];
      qmc[0+i*3] = imgcrd[0+a*3];
      qmc[1+i*3] = imgcrd[1+a*3];
      qmc[2+i*3] = imgcrd[2+a*3];
    }

  for ( int i=0; i<nmm; ++i )
    {
      int a = immatoms[i];
      mmc[0+i*3] = imgcrd[0+a*3];
      mmc[1+i*3] = imgcrd[1+a*3];
      mmc[2+i*3] = imgcrd[2+a*3];
      mmq[i]     = charges[a];
    }

}




void cew::Ewald::AccTmpDensityControlPts
( double const * coord, 
  double const chg,
  double * grd )
{
  int t1[12],t2[12],t3[12];

  double frac[3] = { coord[0] * recip[0] + coord[1] * recip[1] + coord[2] * recip[2],
		     coord[0] * recip[3] + coord[1] * recip[4] + coord[2] * recip[5],
		     coord[0] * recip[6] + coord[1] * recip[7] + coord[2] * recip[8] };
  
  ccdl::bspline_periodic( frac[0], n1, order, w1.data(), bf1.data(), t1 );
  ccdl::bspline_periodic( frac[1], n2, order, w2.data(), bf2.data(), t2 );
  ccdl::bspline_periodic( frac[2], n3, order, w3.data(), bf3.data(), t3 );
  double dedf[3] = { 0.,0.,0. };
  for ( int i=0; i<order; ++i )
    {
      for ( int j=0; j<order; ++j )
	{
	  int    const t12  = (t2[j]+(t1[i])*n2)*n3;
	  double const w1w2 = w1[i]*w2[j];
	  double const w1d2 = w1[i]*bf2[j];
	  double const d1w2 = bf1[i]*w2[j];

	  double const q12 = chg*w1w2;

	  for ( int k=0; k<order; ++k )
	    {
	      double const g = grid[ t3[k]+t12 ];
	      dedf[0] += g * d1w2*w3[k];
	      dedf[1] += g * w1d2*w3[k];
	      dedf[2] += g * w1w2*bf3[k];
	      tmp_grid[ t3[k]+t12 ] += q12*w3[k];
	    };
	};
    };
  grd[0] = chg * ( dedf[0]*recip[0] + dedf[1]*recip[3] + dedf[2]*recip[6] );
  grd[1] = chg * ( dedf[0]*recip[1] + dedf[1]*recip[4] + dedf[2]*recip[7] );
  grd[2] = chg * ( dedf[0]*recip[2] + dedf[1]*recip[5] + dedf[2]*recip[8] );
}



void cew::Ewald::AccDensityControlPts
( double const * coord, 
  double const chg )
{
  int t1[12],t2[12],t3[12];
  double frac[3] = { coord[0] * recip[0] + coord[1] * recip[1] + coord[2] * recip[2],
		     coord[0] * recip[3] + coord[1] * recip[4] + coord[2] * recip[5],
		     coord[0] * recip[6] + coord[1] * recip[7] + coord[2] * recip[8]};
  
  ccdl::bspline_periodic( frac[0], n1, order, w1.data(), t1 );
  ccdl::bspline_periodic( frac[1], n2, order, w2.data(), t2 );
  ccdl::bspline_periodic( frac[2], n3, order, w3.data(), t3 );
  for ( int i=0; i<order; ++i )
    {
      double q1 = w1[i]*chg;
      for ( int j=0; j<order; ++j )
	{
	  int t12 = (t2[j]+(t1[i])*n2)*n3;
	  double q12 = q1*w2[j];
	  for ( int k=0; k<order; ++k )
	    grid[ t3[k]+t12 ] += q12*w3[k];
	};
    };
}


void cew::Ewald::AccDensityControlPts
( int const natom, 
  double const * coord, 
  double const * chg )
{
  int t1[12],t2[12],t3[12];
  for ( int a = 0; a < natom; ++a )
    {
      double frac[3] = { coord[0+a*3] * recip[0] + coord[1+a*3] * recip[1] + coord[2+a*3] * recip[2],
                         coord[0+a*3] * recip[3] + coord[1+a*3] * recip[4] + coord[2+a*3] * recip[5],
                         coord[0+a*3] * recip[6] + coord[1+a*3] * recip[7] + coord[2+a*3] * recip[8]};

      ccdl::bspline_periodic( frac[0], n1, order, w1.data(), t1 );
      ccdl::bspline_periodic( frac[1], n2, order, w2.data(), t2 );
      ccdl::bspline_periodic( frac[2], n3, order, w3.data(), t3 );

      for ( int i=0; i<order; ++i )
        {
          double q1 = w1[i]*chg[a];
          for ( int j=0; j<order; ++j )
            {
              int t12 = (t2[j]+(t1[i])*n2)*n3;
              double q12 = q1*w2[j];
              for ( int k=0; k<order; ++k )
                grid[ t3[k]+t12 ] += q12*w3[k];
            };
        };
    };
}

#include <fstream>
#include <iomanip>

void cew::Ewald::CptPlaneWavePotControlPts()
{
  double Qtot = 0.;
  for ( int a=0; a<nt; ++a ) Qtot += grid[a];
  ChargeCorrection = - (ccdl::PI/zeta)*(Qtot/nt);

  // Compute the DFT coefficients of the (spreaded) Gaussian charge density
  fftw_execute( *fplan );

  // Convolute with the kernel
  for ( int k = 0; k < nk; ++k )
    dftcoef[k] *= kernel[k];
  
  // Evaluate the recip-space potential (control points) at the grid points
  fftw_execute( *rplan );


}



void cew::Ewald::CptPlaneWavePotControlPts
( int const natom, 
  double const * coord, 
  double const * chg )
{
  std::fill( grid, grid + nt, 0. );
  AccDensityControlPts( natom, coord, chg );
  CptPlaneWavePotControlPts();
}






double cew::Ewald::CptPlaneWavePot( double const * pt )
{
  double p = ChargeCorrection;

  int t1[12],t2[12],t3[12];

  double frac[3] = { pt[0] * recip[0] + pt[1] * recip[1] + pt[2] * recip[2],
		     pt[0] * recip[3] + pt[1] * recip[4] + pt[2] * recip[5],
		     pt[0] * recip[6] + pt[1] * recip[7] + pt[2] * recip[8] };
  
  ccdl::bspline_periodic( frac[0], n1, order, w1.data(), t1 );
  ccdl::bspline_periodic( frac[1], n2, order, w2.data(), t2 );
  ccdl::bspline_periodic( frac[2], n3, order, w3.data(), t3 );
  
  for ( int i=0; i<order; ++i )
    for ( int j=0; j<order; ++j )
      {
	int    const t12 = (t2[j]+(t1[i])*n2)*n3;
	double const w12 = w1[i]*w2[j];
	for ( int k=0; k<order; ++k )
	  p += grid[ t3[k]+t12 ] * w12*w3[k];
      };
  return p;
}


double cew::Ewald::CptPlaneWavePotAndGrd( double const * pt, double * grd )
{
  double p = ChargeCorrection;

  int t1[12],t2[12],t3[12];

  double frac[3] = { pt[0] * recip[0] + pt[1] * recip[1] + pt[2] * recip[2],
		     pt[0] * recip[3] + pt[1] * recip[4] + pt[2] * recip[5],
		     pt[0] * recip[6] + pt[1] * recip[7] + pt[2] * recip[8] };
  
  ccdl::bspline_periodic( frac[0], n1, order, w1.data(), bf1.data(), t1 );
  ccdl::bspline_periodic( frac[1], n2, order, w2.data(), bf2.data(), t2 );
  ccdl::bspline_periodic( frac[2], n3, order, w3.data(), bf3.data(), t3 );
  double dedf[3] = { 0.,0.,0. };
  for ( int i=0; i<order; ++i )
    {
      for ( int j=0; j<order; ++j )
	{
	  int    const t12  = (t2[j]+(t1[i])*n2)*n3;
	  double const w1w2 = w1[i]*w2[j];
	  double const w1d2 = w1[i]*bf2[j];
	  double const d1w2 = bf1[i]*w2[j];
	  double const w12  = w1[i]*w2[j];

	  for ( int k=0; k<order; ++k )
	    {
	      double const g = grid[ t3[k]+t12 ];
	      p      += g * w12*w3[k];
	      dedf[0] += g * d1w2*w3[k];
	      dedf[1] += g * w1d2*w3[k];
	      dedf[2] += g * w1w2*bf3[k];
	    };
	};
    };
  grd[0] = ( dedf[0]*recip[0] + dedf[1]*recip[3] + dedf[2]*recip[6] );
  grd[1] = ( dedf[0]*recip[1] + dedf[1]*recip[4] + dedf[2]*recip[7] );
  grd[2] = ( dedf[0]*recip[2] + dedf[1]*recip[5] + dedf[2]*recip[8] );

  return p;
}



void cew::Ewald::CptPlaneWaveGrd( double const * pt, double * grd )
{
  int t1[12],t2[12],t3[12];

  double frac[3] = { pt[0] * recip[0] + pt[1] * recip[1] + pt[2] * recip[2],
		     pt[0] * recip[3] + pt[1] * recip[4] + pt[2] * recip[5],
		     pt[0] * recip[6] + pt[1] * recip[7] + pt[2] * recip[8] };
  
  ccdl::bspline_periodic( frac[0], n1, order, w1.data(), bf1.data(), t1 );
  ccdl::bspline_periodic( frac[1], n2, order, w2.data(), bf2.data(), t2 );
  ccdl::bspline_periodic( frac[2], n3, order, w3.data(), bf3.data(), t3 );
  double dedf[3] = { 0.,0.,0. };
  for ( int i=0; i<order; ++i )
    {
      for ( int j=0; j<order; ++j )
	{
	  int    const t12  = (t2[j]+(t1[i])*n2)*n3;
	  double const w1w2 = w1[i]*w2[j];
	  double const w1d2 = w1[i]*bf2[j];
	  double const d1w2 = bf1[i]*w2[j];

	  for ( int k=0; k<order; ++k )
	    {
	      double const g = grid[ t3[k]+t12 ];
	      dedf[0] += g * d1w2*w3[k];
	      dedf[1] += g * w1d2*w3[k];
	      dedf[2] += g * w1w2*bf3[k];
	    };
	};
    };
  grd[0] = ( dedf[0]*recip[0] + dedf[1]*recip[3] + dedf[2]*recip[6] );
  grd[1] = ( dedf[0]*recip[1] + dedf[1]*recip[4] + dedf[2]*recip[7] );
  grd[2] = ( dedf[0]*recip[2] + dedf[1]*recip[5] + dedf[2]*recip[8] );
}
