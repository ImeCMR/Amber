#ifndef _edgembar_CalcIdx_hpp_
#define _edgembar_CalcIdx_hpp_

#include <ostream>

namespace edgembar
{
  
  struct CalcIdx
  {
    CalcIdx()
      : BaseIdx(0),
	IsFwd(false),
	IsRev(false),
	IsCon(false),
	IsHalf1(false),
	IsHalf2(false),
	IsBootstrap(false),
	RangeLo(0),
	RangeHi(1),
	ConVal(0),
	Chisq(0)
    {}
    
    int BaseIdx;
    
    bool IsFwd;
    bool IsRev;
    bool IsCon;
    bool IsHalf1;
    bool IsHalf2;
    
    bool IsBootstrap;
    
    double RangeLo;
    double RangeHi;
    double ConVal;
    double Chisq;

    void Write( std::ostream & cout ) const;
  };

}

#endif
