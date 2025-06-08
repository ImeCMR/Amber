#include "Sample.hpp"

ndfes::Sample::Sample()
  : stateidx(0),
    hamidx(0),
    sbinidx(0)
{
}

ndfes::Sample::Sample
( std::size_t const iistate,
  std::size_t const indim,
  double const * ipt,
  std::size_t const inham,
  std::size_t const iiham,
  double const * iunbiasedpotenes )
  : stateidx(iistate),
    hamidx(iiham),
    sbinidx(0),
    pt( ipt, ipt+indim ),
    potenes( iunbiasedpotenes, iunbiasedpotenes + inham )
{
}


void ndfes::Sample::AddToPotEnes( double const * c )
{
  for ( std::size_t i=0, n=potenes.size(); i<n; ++i )
    {
      potenes[i] += c[i];
    }
}
