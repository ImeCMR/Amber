#ifndef _edgembar_Trial_hpp_
#define _edgembar_Trial_hpp_

#include <string>
#include <vector>
#include <ostream>

#include "ReadMode.hpp"
#include "Sim.hpp"
#include "CalcIdx.hpp"
#include "StateValues.hpp"


namespace edgembar
{
  class Edge;
  class Env;
  class Stage;
}

namespace edgembar
{

  class Trial
  {
  public:


    Trial();

    Trial( std::string const name,
	   std::vector<std::string> const & statelabels,
	   std::string const & datadir,
	   double const beta,
	   edgembar::ReadMode const readmode,
	   int const autoeqmode,
	   double const shift );

    void StoreAvgEnes();
    
    void SetLinkedList( edgembar::Edge * edge,
			edgembar::Env * env,
			edgembar::Stage * stage );

    edgembar::Trial ExtractRange( double const start, double const stop ) const;

    edgembar::Trial ExtractProd() const;
    
    edgembar::Trial ExtractBootstrap() const;


    void ExtractRangeInPlace( double const start, double const stop );

    void ExtractProdInPlace();
    
    void ExtractBootstrapInPlace();


    int GetNumStates() const { return mStateLabels.size(); }
    
    double CptObjective( double const * p, double * g ) const;
    void TestObjective() const;
    
    std::vector<edgembar::Sim *> GetSims();
    std::vector<edgembar::Sim const *> GetConstSims() const;

    void WriteDebugInfo( std::ostream & cout ) const;

    std::string GetName() const { return mName; }
    
    edgembar::Edge * GetEdge() { return mEdge; };
    edgembar::Env * GetEnv() { return mEnv; };
    edgembar::Stage * GetStage() { return mStage; };
    
    edgembar::Edge const * GetEdge() const { return mEdge; };
    edgembar::Env const * GetEnv() const { return mEnv; };
    edgembar::Stage const * GetStage() const { return mStage; };
    

    void SetParamOffset( int o ) { mParamOffset=o; }
    int GetParamOffset() const { return mParamOffset; }

    void MakeExpAvgGuess( double * p ) const;

    std::vector<double> GetOverlaps() const;

    std::vector<double> GetEntropies() const;
    
    std::string GetPython() const;

    void UseUniformWts();

    double GetAvgEne( int i ) const { return mAvgEnes[i]; }

    double GetShift() const { return mShift; }
    
  private:

    std::string mName;
    std::vector<std::string> mStateLabels;
    std::string mDataDir;
    edgembar::ReadMode mReadMode;
    double mBeta;
    double mShift;

    int mParamOffset;
    
    edgembar::Edge * mEdge;
    edgembar::Env * mEnv;
    edgembar::Stage * mStage;

    std::vector<edgembar::Sim> mSims;
    std::vector<double> mAvgEnes;
    bool mUniformWts;
    
  };
  
}

#endif

