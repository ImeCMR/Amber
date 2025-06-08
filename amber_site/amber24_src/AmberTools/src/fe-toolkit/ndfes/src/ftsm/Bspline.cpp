#include <cmath>
#include <iostream>
#include <cassert>
#include "Bspline.hpp"


namespace ndfes
{
  void bspline_one_pass( double * c, double const w, int const n );
  void bspline_diff( int const order, double const * array, double * diff );
  void bspline_eval( double const w, int const order, double * array ); 
  void bspline_eval( double const w, int const order, double * array, double * darray );
}


void ndfes::bspline_one_pass( double * c, double const w, int const n )
{
  int nm1 = n-1;
  double const div = 1. / nm1;
  c[nm1] = div * w * c[nm1 - 1];
  for ( int j=1; j<nm1; ++j )
    c[nm1 - j] = div * ((w + j) * c[nm1 - j - 1] + (n - j - w) * c[nm1 - j]);
  c[0] = div * (1 - w) * c[0];
}

void ndfes::bspline_diff( int const order, double const * array, double * diff )
{
  assert( order > 1 );
  int const nm1 = order-1;
  diff[0] = -array[0];
  for ( int j=1; j<nm1; ++j )
    diff[j] = array[j-1] - array[j];
  diff[nm1] = array[nm1-1];
}


void ndfes::bspline_eval( double const w, int const order, double * array )
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
	    {
	      ndfes::bspline_one_pass(array, w, k);
	    }
        };
    };
}

void ndfes::bspline_eval( double const w, int const order, double * array, double * darray )
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
	{
	  ndfes::bspline_one_pass(array, w, k);
	}

      ndfes::bspline_diff(order,array,darray);

      // One more recursion: // do the k==order
      ndfes::bspline_one_pass(array, w, order);

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


void ndfes::bspline_aperiodic
( double x, 
  double const lenx, 
  int const nx,
  int const order, 
  double * w,
  int * gidx )
{
  //double const ix = x;
  x = (x/lenx)*(nx-1) - (order%2)*0.5;
  int ilo = std::floor(x);
  //std::printf("%23.14e %20.14e ilo %i\n",x-ilo,x,ilo);
  
  // maybe this should really be
  // if ( ilo + 1 - order/2 == -1 )
  if ( ilo < nx/2-1 )
    {
      if ( x-ilo > 1. - 5.e-14 )
        {
          //std::printf("boundary case increase ilo %20.10e %i\n",x,ilo);
          ilo += 1;
          x = ilo;
        }
    }
  
  else if ( ilo + 1 - (order+1)/2 + (order-1) == nx )
    {
      if ( x-ilo < 5.e-14 )
        {
          //std::printf("boundary case decrease ilo %20.10e %i\n",x,ilo);
          ilo -= 1;
          x = ilo+1;
        }
    }
  
  //else if ( ilo > order/2 and x-ilo
  //std::printf("x,ilo: %26.17e %i\n",x,ilo);
  if ( order == 1 )
    {
      w[0] = 1.;
      gidx[0] = (x-ilo < 0.5) ? ilo : ilo+1;
      //std::printf("x %20.10e %20.10e lenx %20.10e nx %4i order %4i gidx[%i] %i (ilo=%i) = %20.10e\n",ix,x-ilo,lenx,nx,order,0,gidx[0],ilo,w[0]);
      //std::cout.flush();

    }
  else
    {
      //std::printf("bspline frac %20.10e %i\n",x-ilo,ilo);
      ndfes::bspline_eval(x-ilo,order,w);
      ilo += 1 - order/2;
      for ( int b=0; b<order; ++b, ++ilo )
        {
          gidx[b] = ilo;
          //std::printf("x %20.10e lenx %20.10e nx %4i order %4i gidx[%i] %i = %20.10e\n",ix,lenx,nx,order,b,ilo,w[b]);
          if ( ilo < 0 or ilo > nx-1 )
            {
              std::cerr << "ndfes::bspline_aperiodic node index out of bounds " << ilo << " " << nx << "\n";
              std::exit(1);
            }
      //gidx[b] = std::min(nx-1,std::max(0,ilo));
      //gidx[b] = ((ilo+nx)%nx+nx)%nx;
      //if ( ilo < 0 or ilo >= nx ) w[b]=0.;
      //std::printf("%4i(%4i) ",gidx[b],nx);
        }
    }
  //std::printf("\n");
}




void ndfes::bspline_aperiodic
( double x, 
  double const lenx, 
  int const nx,
  int const order, 
  double * w,
  double * dw,
  int * gidx )
{
  //double const ix = x;
  double dx = ((double)(nx-1.)/((double)lenx));
  x = (x/lenx)*(nx-1) - (order%2)*0.5;
  int ilo = std::floor(x);
  //std::printf("%23.14e %20.14e ilo %i\n",x-ilo,x,ilo);
  // maybe this should really be
  // if ( ilo + 1 - order/2 == -1 )
  if ( ilo < nx/2-1 )
    {
      if ( x-ilo > 1. - 5.e-14 )
        {
          //std::printf("boundary case increase ilo %20.10e %i\n",x,ilo);
          ilo += 1;
          x = ilo;
        }
    }
  
  else if ( ilo + 1 - (order+1)/2 + (order-1) == nx )
    {
      if ( x-ilo < 5.e-14 )
        {
          //std::printf("boundary case decrease ilo %20.10e %i\n",x,ilo);
          ilo -= 1;
          x = ilo+1;
        }
    }
  
  //else if ( ilo > order/2 and x-ilo
  //std::printf("x,ilo: %26.17e %i\n",x,ilo);
  if ( order == 1 )
    {
      w[0] = 1.;
      dw[0] = 0.;
      gidx[0] = (x-ilo < 0.5) ? ilo : ilo+1;
      //std::printf("x %20.10e %20.10e lenx %20.10e nx %4i order %4i gidx[%i] %i (ilo=%i) = %20.10e\n",ix,x-ilo,lenx,nx,order,0,gidx[0],ilo,w[0]);
      //std::cout.flush();

    }
  else
    {
      //std::printf("bspline frac %20.10e %i\n",x-ilo,ilo);
      ndfes::bspline_eval(x-ilo,order,w,dw);
      ilo += 1 - order/2;
      for ( int b=0; b<order; ++b, ++ilo )
        {
          gidx[b] = ilo;
	  dw[b] *= dx;
          //std::printf("x %20.10e lenx %20.10e nx %4i order %4i gidx[%i] %i = %20.10e\n",ix,lenx,nx,order,b,ilo,w[b]);
          if ( ilo < 0 or ilo > nx-1 )
            {
              std::cerr << "ndfes::bspline_aperiodic (deriv) node index out of bounds " << ilo << " " << nx << "\n";
              std::exit(1);
            }
      //gidx[b] = std::min(nx-1,std::max(0,ilo));
      //gidx[b] = ((ilo+nx)%nx+nx)%nx;
      //if ( ilo < 0 or ilo >= nx ) w[b]=0.;
      //std::printf("%4i(%4i) ",gidx[b],nx);
        }
    }
  //std::printf("\n");
}










void ndfes::bspline_periodic
( double x, 
  double const lenx, 
  int const nx,
  int const order, 
  double * w,
  int * gidx )
{
  x = (x/lenx)*nx - (order%2)*0.5;
  int ilo = std::floor(x);
  x -= ilo;
  if ( order == 1 )
    {
      w[0] = 1.;
      gidx[0] = (x < 0.5) ? ilo : ilo+1;
    }
  else
    {
      ndfes::bspline_eval(x,order,w);
      // the first modulus takes us to [-nx+1,nx-1]
      // then we add nx, so we are in the range [0,2*nx-1]
      // the second modulus then puts us in [0,nx-1]
      // Each ilo in the loop uses a modulus, so we don't need
      // the second modulus on the first calculation of ilo
      ilo   =  (ilo+1-order/2)%nx+nx;
      for ( int b=0; b<order; ++b, ++ilo )
        {
          gidx[b] = ilo%nx;
        }
    }
  //std::printf("\n");
}



void ndfes::bspline_periodic
( double x, 
  double const lenx, 
  int const nx,
  int const order, 
  double * w,
  double * dw,
  int * gidx )
{
  double dx = ((double)(nx-1.)/((double)lenx));
  x = (x/lenx)*nx - (order%2)*0.5;
  int ilo = std::floor(x);
  x -= ilo;
  if ( order == 1 )
    {
      w[0] = 1.;
      gidx[0] = (x < 0.5) ? ilo : ilo+1;
    }
  else
    {
      ndfes::bspline_eval(x,order,w,dw);
      // the first modulus takes us to [-nx+1,nx-1]
      // then we add nx, so we are in the range [0,2*nx-1]
      // the second modulus then puts us in [0,nx-1]
      // Each ilo in the loop uses a modulus, so we don't need
      // the second modulus on the first calculation of ilo
      ilo   =  (ilo+1-order/2)%nx+nx;
      for ( int b=0; b<order; ++b, ++ilo )
        {
          gidx[b] = ilo%nx;
	  dw[b] *= dx;
        }
    }
  //std::printf("\n");
}


