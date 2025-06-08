#ifndef _edgembar_Stage_hpp_
#define _edgembar_Stage_hpp_

#include <string>
#include <vector>

#include "ReadMode.hpp"
#include "Trial.hpp"
#include "CalcIdx.hpp"
#include "StateValues.hpp"


namespace edgembar
{

  class Stage
  {
  public:

    Stage();
    
    Stage( std::string const name,
	   std::vector<edgembar::Trial> const & trials );

    void SetLinkedList( edgembar::Edge * edge,
                        edgembar::Env * env );


    void ExtractRangeInPlace( double const start, double const stop );

    void ExtractProdInPlace();
    
    void ExtractBootstrapInPlace();

    std::vector<edgembar::Sim *> GetSims();
    std::vector<edgembar::Sim const *> GetConstSims() const;
    
    std::vector<edgembar::Trial *> GetTrials();
    std::vector<edgembar::Trial const *> GetConstTrials() const;
    
    std::string GetName() const { return mName; }
    
    edgembar::Edge * GetEdge() { return mEdge; }
    edgembar::Env * GetEnv() { return mEnv; }
    
    edgembar::Edge const * GetEdge() const { return mEdge; }
    edgembar::Env const * GetEnv() const { return mEnv; }

    edgembar::Trial * GetTrialBegin() { return mTrials.data(); }
    edgembar::Trial * GetTrialEnd() { return mTrials.data()+mTrials.size(); }

    int GetNumTrials() const { return mTrials.size(); }
    
    std::string GetPython() const;
    
  private:

    std::string mName;
    std::vector<edgembar::Trial> mTrials;
    
    edgembar::Edge * mEdge;
    edgembar::Env * mEnv;

  };
  
}

#endif

