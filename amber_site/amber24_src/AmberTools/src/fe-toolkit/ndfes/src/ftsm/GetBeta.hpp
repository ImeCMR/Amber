#ifndef _get_beta_hpp_
#define _get_beta_hpp_

namespace ndfes
{
  double GetBeta( double const T );
}

inline double ndfes::GetBeta(double const T)
{
  double const k_au = 3.16681534222374e-06;
  double const au_per_kcal = 1.59360145069066e-03;
  double const k_kcal = (k_au / au_per_kcal);
  return 1. / ( k_kcal * T );
}

#endif
