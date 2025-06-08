#ifndef _edgembar_StateValues_hpp_
#define _edgembar_StateValues_hpp_

#include <vector>
#include <string>

namespace edgembar
{
  class StateValues
  {
  public:
    
    StateValues() { mResult.clear(); mBootstrap.clear();};

    void SetResult( std::vector<double> const & x );
    
    void PushBootstrap( std::vector<double> const & x );
    
    void Join( edgembar::StateValues const & x );
    
    std::vector<double> GetResult() const { return mResult; }
    
    std::vector<double> GetStderr() const;

    void reset();
    
  private:    
    std::vector<double> mResult;
    std::vector< std::vector<double> > mBootstrap;
    void Assert( std::string msg, std::vector<double> const & x ) const;
  };
}

#endif

