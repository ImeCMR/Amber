#ifndef _edgembar_ReadMode_hpp_
#define _edgembar_ReadMode_hpp_

#include <string>
#include <vector>

namespace edgembar
{
  enum ReadMode { INVALID=-1, AUTO=0,
		  MBAR=1, MBAREXP0=2, MBAREXP1=3, MBAREXP=4,
		  BAR=5, BAREXP0=6, BAREXP1=7, BAREXP=8 };

  edgembar::ReadMode GetMode( std::string mode );

  std::string GetModeString( edgembar::ReadMode mode );
  
  bool CanEvalMode( int const ns,
		    int const * avail,
		    edgembar::ReadMode const mode );

  edgembar::ReadMode AutoDetectMode( int const ns,
				     int const * avail );

  std::vector< std::vector<int> > GetSimIndexes
  ( int const ns,
    edgembar::ReadMode const mode );
  
}

#endif
