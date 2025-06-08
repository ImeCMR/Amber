#include "State.hpp"
#include "GetBeta.hpp"
#include "PeriodicUtils.hpp"



ndfes::State::State()
  : ndim(0),
    iham(0),
    temperature(298.),
    beta(0.)
{}


ndfes::State::State
( std::size_t const indim,
  std::size_t const iiham,
  double const * icenter,
  double const * ifconst,
  double const itemperature )
  : ndim(indim),
    iham(iiham),
    center(icenter, icenter+indim),
    fconst(ifconst, ifconst+indim),
    temperature(itemperature),
    beta( ndfes::GetBeta(itemperature) ),
    statineff(1)
{
  //for ( int i=0; i<ndim; ++i )
  //fconst[i] = ConvertKcalToKT( fconst[i] );
}


double ndfes::State::CptBiasEnergy( int const * dimisper, double const * pt ) const
{
  double W = 0.;
  for ( std::size_t idim=0; idim<ndim; ++idim )
    {
      double dx = pt[idim] - center[idim];
      if ( dimisper[idim] > 0 )
	{
	  dx = wrap(dx,360.);
	}
      W += fconst[idim] * dx * dx;
    }
  return W;
}


bool ndfes::State::operator==( ndfes::State const & rhs ) const
{
   double const TOL = 1.e-10;
   double diff=0.;
   double x = 0.;
   x = iham - rhs.iham;
   diff += x*x;
   x = center.size() - rhs.center.size();
   diff += x*x;
   x = temperature - rhs.temperature;
   diff += x*x;
   if ( diff < TOL )
     {
         for ( std::size_t dim=0; dim<ndim; ++dim )
	 {
              x = center[dim] - rhs.center[dim];
	      diff += x*x;
	      x = fconst[dim] - rhs.fconst[dim];
	      diff += x*x;
	 }
     }
   return diff < TOL;
}

