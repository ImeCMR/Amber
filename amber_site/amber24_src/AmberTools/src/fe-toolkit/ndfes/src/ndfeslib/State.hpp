#ifndef _ndfes_State_hpp_
#define _ndfes_State_hpp_

#include <vector>
#include <memory>

namespace ndfes
{
  class State
  {
  public:
    
    State();
    
    State( std::size_t const ndim,
	   std::size_t const iham,
	   double const * center,
	   double const * fconst,
	   double const temperature );

    std::size_t GetHamIdx() const { return iham; }

    std::size_t GetNumDims() const { return ndim; }

    double GetTemperature() const { return temperature; }

    double GetBeta() const { return beta; }
    
    double const * GetCenter() const { return center.data(); };

    double const * GetFConst() const { return fconst.data(); };

    double CptBiasEnergy( int const * dimisper, double const * pt ) const;

    double ConvertKcalToKT( double const E ) const { return beta*E; }

    double ConvertKTToKcal( double const E ) const { return E/beta; }

    bool operator==( ndfes::State const & rhs ) const;

    void SetStatIneff( std::size_t const g  ) { statineff = g; }

    std::size_t GetStatIneff() const { return statineff; }

  private:

    std::size_t ndim;
    std::size_t iham;
    std::vector<double> center;
    std::vector<double> fconst;
    double temperature;
    double beta;
    //double biased_energy;
    std::size_t statineff;
  };


  typedef std::shared_ptr< ndfes::State > pState;
}

#endif
