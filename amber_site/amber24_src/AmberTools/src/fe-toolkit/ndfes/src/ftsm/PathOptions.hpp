#ifndef _ndfes_PathOptions_hpp_
#define _ndfes_PathOptions_hpp_

#include <vector>

namespace ndfes
{

  struct PathOptions
  {
  public:
    
    PathOptions()
      : npathpts(-1),
	stype(1),
	smooth(true),
        fix0(false),
        fix1(false),
        scalefc(false),
        smoothfc(false),
        tdispfc(0.75),
	acc(true),
        acc_oobk(300.),
        acc_maxit(100),
        beta(1.6886561526804216)
    {}

    void ClearBounds()
    {
      minbounds.resize(0);
      maxbounds.resize(0);
    }
    
    int npathpts;
    
    int stype;
    bool smooth;
    
    bool fix0;
    bool fix1;

    bool scalefc;
    bool smoothfc;
    double tdispfc;
    
    std::vector<double> deffc;
    std::vector<double> minfc;
    std::vector<double> maxfc;

    bool acc;
    double acc_oobk;
    std::size_t acc_maxit;
    double beta;
    
    std::vector<double> minbounds;
    std::vector<double> maxbounds;
  };
  

}

#endif
