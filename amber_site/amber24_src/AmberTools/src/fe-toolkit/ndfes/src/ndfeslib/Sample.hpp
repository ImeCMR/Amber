#ifndef _ndfes_Sample_hpp_
#define _ndfes_Sample_hpp_

#include <memory>
#include <vector>

namespace ndfes
{
  class Sample
  {
  public:

    Sample();

    Sample( std::size_t const istate,
	    std::size_t const ndim,
	    double const * pt,
	    std::size_t const nham,
	    std::size_t const iham,
	    double const * unbiasedpotenes );

    std::size_t GetStateIdx() const { return stateidx; };
    
    std::size_t GetHamIdx() const { return hamidx; };
    
    std::size_t GetSpatialBinArrayIdx() const { return sbinidx; };

    void SetSpatialBinArrayIdx( std::size_t i ) { sbinidx=i; };
    
    double const * GetPt() const { return pt.data(); };

    double const * GetPotEnes() const { return potenes.data(); };

    double GetPotEne( std::size_t i ) const { return potenes[i]; };

    double GetDeltaEne( std::size_t i ) const { return potenes[i]-potenes[hamidx]; };

    double GetMaxBias() const { return maxbias; }

    void SetMaxBias( double w ) { maxbias = w; }

    void AddToPotEnes( double const * c );
    
  private:

    std::size_t stateidx;
    std::size_t hamidx;
    std::size_t sbinidx;
    std::vector<double> pt;
    std::vector<double> potenes;
    double maxbias;
    
  };
  
  typedef std::shared_ptr< ndfes::Sample > pSample;
}

#endif
