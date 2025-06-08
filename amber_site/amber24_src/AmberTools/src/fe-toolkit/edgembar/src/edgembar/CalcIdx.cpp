#include <iomanip>

#include "CalcIdx.hpp"


void edgembar::CalcIdx::Write( std::ostream & cout ) const
{
  cout << "BaseIdx " << BaseIdx;
  
  if ( IsBootstrap )
    {
      cout << " IsBootstrap";
    }
  
  if ( IsFwd )
    {
      cout << " IsFwd";
    }
  if ( IsRev )
    {
      cout << " IsRev";
    }
  if ( IsHalf1 )
    {
      cout << " IsHalf1";
    }
  if ( IsHalf2 )
    {
      cout << " IsHalf2";
    }
  if ( IsFwd or IsRev or IsHalf1 or IsHalf2 )
    {
      cout << " RangeLo " << std::fixed << std::setw(5) << std::setprecision(2) << RangeLo;
      cout << " RangeHi " << std::fixed << std::setw(5) << std::setprecision(2) << RangeHi;
    }

  if ( IsCon )
    {
      cout << " IsCon";
      cout << " ConVal " << std::fixed << std::setw(10) << std::setprecision(4) << ConVal;
    }
  
  
}
