#ifndef _RunAvg_hpp_
#define _RunAvg_hpp_

#include <cmath>



namespace ndfes
{


  class RunAvg
  {
  public:
    RunAvg();
    void reset();
    void push_back( double x );
    std::size_t size() const;
    double mean() const;
    double var() const;
    double stddev() const;
    double stderr( bool isbootstrap ) const;
  private:
    std::size_t numpts;
    double oldmu;
    double oldS;
    double newmu;
    double newS;
  };
  

  
}



inline ndfes::RunAvg::RunAvg()
  : numpts(0),
    oldmu(0.),
    oldS(0.),
    newmu(0.),
    newS(0.)
{}

inline void ndfes::RunAvg::reset()
{
  numpts = 0;
  oldmu = 0.;
  oldS = 0.;
  newmu = 0.;
  newS = 0.;
}

inline std::size_t ndfes::RunAvg::size() const
{
  return numpts;
}

inline double ndfes::RunAvg::mean() const
{
  return (numpts > 0) ? newmu : 0.0;
}

inline double ndfes::RunAvg::var() const
{
  return ( (numpts > 1) ? (newS/(numpts - 1)) : 0.0 );
}

inline double ndfes::RunAvg::stddev() const
{
  return ( (numpts > 1) ? std::sqrt( var() ) : 0.0 );
}

inline double ndfes::RunAvg::stderr( bool isbootstrap ) const
{
  return isbootstrap ? stddev() : ( (numpts > 1) ? std::sqrt(var() / numpts) : 0.0 );
}


#endif

