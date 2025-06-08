#ifndef _edgembar_Env_hpp_
#define _edgembar_Env_hpp_

#include <string>
#include <vector>
#include "Stage.hpp"
#include "CalcIdx.hpp"
#include "StateValues.hpp"

namespace edgembar
{
  
  class Env
  {
  public:

    Env();
    
    Env( std::string const name,
	 std::vector<edgembar::Stage> const & stages );

    void SetLinkedList( edgembar::Edge * edge );

    void ExtractRangeInPlace( double const start, double const stop );

    void ExtractProdInPlace();
    
    void ExtractBootstrapInPlace();

    std::vector<edgembar::Sim *> GetSims();
    std::vector<edgembar::Sim const *> GetConstSims() const;
    
    std::vector<edgembar::Trial *> GetTrials();
    std::vector<edgembar::Trial const *> GetConstTrials() const;

    std::string GetName() const { return mName; }
    
    edgembar::Edge * GetEdge() { return mEdge; }
    
    edgembar::Edge const * GetEdge() const { return mEdge; }

    edgembar::Stage * GetStageBegin() { return mStages.data(); }
    edgembar::Stage * GetStageEnd() { return mStages.data()+mStages.size(); }

    
    std::string GetPython() const;
    
    
  private:
    std::string mName;
    std::vector<edgembar::Stage> mStages;

    edgembar::Edge * mEdge;
  };
  
}

#endif
