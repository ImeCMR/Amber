#include <algorithm>
#include "PeriodicUtils.hpp"


extern "C"
{
  extern void
  dger_(const int *m, const int *n, double const *alpha,
        double const *x, const int *incx,
        double const *y, const int *incy,
        double *a, const int *lda);
}



double wraptorange( double const x, double const lo, double const hi )
{
  double mid = 0.5*(lo+hi);
  mid = wrap(x-mid,hi-lo)+mid;
  if ( mid == hi ) mid = lo;
  return mid;
}



// void ndfes::LinearWtsToMeshWts
// ( std::vector< std::vector<double> > const & lsp,
//   std::vector<double> & Ws )
// {
//   int ndim = lsp.size();
  
//   Ws = lsp[0];

//   int inc=1;
//   double alpha=1.;
//   for ( int dim=1; dim < ndim; ++dim )
//     {
//       std::vector<double> Wfast(Ws);
//       std::vector<double> Wslow(lsp[dim]);
//       int nf=Wfast.size();
//       int ns=Wslow.size();
//       Ws.resize( nf*ns );
//       std::fill( Ws.data(), Ws.data() + Ws.size(), 0. );
//       dger_( &nf, &ns, &alpha, Wfast.data(), &inc,
// 	     Wslow.data(), &inc, Ws.data(), &nf );
//     }
// }
