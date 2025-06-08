#ifndef _edgembar_Sim_hpp_
#define _edgembar_Sim_hpp_

#include <vector>
#include <string>
#include <iostream>
#include "CalcIdx.hpp"
#include "StateValues.hpp"



namespace edgembar
{
  class Sim
  {
  public:

    Sim();

    Sim( int const simidx );

    Sim( int const simidx,
	 std::vector<int> eneidxs,
	 std::vector<std::string> const & statelabels,
	 std::string const datadir,
	 double const beta,
	 int const autoeqmode );

    void ReadFiles( double const fstart, double const fstop, int const stride );
    void OnlyKeepRange( double const start, double const stop );
    void OnlyKeepProd();
    void Bootstrap();

    void WriteDebugInfo( std::ostream & cout ) const;

    int FindLocalIdx( int const eneidx ) const;
    int FindLocalIdx( edgembar::Sim const * other ) const;

    void PrecomputeExps( double const ptol );
    
  public:
    
    int SimIdx;
    std::vector<int> EneIdxs;
    std::vector<std::string> StateLabels;
    std::string DataDir;
    double Beta;
    int AutoEqMode;
    int NumSamples;
    int OrigNumSamples;

    int LocalSimIdx;
    std::vector<double> Emat;
    std::vector<double> Mdat;
    std::vector<double> MaxZ;
    std::vector<double> AvgEnes;
    std::vector<double> DVDL;
    
    int CurStride;
    int OrigStride;
    int ProdStart;
    int ProdStride;
    bool IsConverged;

  private:
    
    int CptStride() const;
    void CptAutoEquil( int & start, int & stride, bool & isconv, double const ptol ) const;

    
  };
}

#endif
