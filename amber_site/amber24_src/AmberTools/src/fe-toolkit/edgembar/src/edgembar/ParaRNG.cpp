#include <vector>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "ParaRNG.hpp"



namespace ParaRNG
{
  class RNG
  {
  public:
    RNG();
    int GetInt( std::uniform_int_distribution<int> & dist );
    std::random_device device;
    std::vector<std::default_random_engine> rgen;
  };
  
}



ParaRNG::RNG::RNG()
{
  int N = 1;
#ifdef WITH_OPENMP
  N = omp_get_max_threads();
#endif
  for ( int i=0; i<N; ++i )
    {
      rgen.emplace_back( std::default_random_engine(device()) );
      //rgen.emplace_back( std::default_random_engine(0) );
    }
}


int ParaRNG::RNG::GetInt( std::uniform_int_distribution<int> & dist )
{
#ifdef WITH_OPENMP
  return dist(rgen[omp_get_thread_num()]);
#else
  return dist(rgen[0]);
#endif
}




static ParaRNG::RNG StaticRNG;


int ParaRNG::RandInt( std::uniform_int_distribution<int> & dist )
{
  return StaticRNG.GetInt(dist);
}
