#ifndef _ndfes_Bspline_hpp_
#define _ndfes_Bspline_hpp_

namespace ndfes
{
  void bspline_aperiodic( double x, double const lenx, int const nx, int const order, double * w, int * gidx );
  
  void bspline_aperiodic( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx );

  void bspline_periodic( double x, double const lenx, int const nx, int const order, double * w, int * gidx );

  void bspline_periodic( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx );

}

#endif
