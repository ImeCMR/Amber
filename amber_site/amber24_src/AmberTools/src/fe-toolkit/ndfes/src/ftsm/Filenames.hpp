#ifndef _ndfes_Filenames_hpp_
#define _ndfes_Filenames_hpp_

#include <string>
#include <sstream>
#include <iomanip>

namespace ndfes
{  
  std::string GetDirName( int i, int pad );

  std::string GetImgName( int i, int pad );

  std::string GetInitName( int i, int pad );
}



inline std::string ndfes::GetDirName( int i, int pad )
{
  std::string dname;
  if ( i == 0 )
    {
      dname = "init";
    }
  else
    {
      std::stringstream sdir;
      sdir << "it" <<  std::setfill('0') << std::setw(pad) << i;
      dname = sdir.str();
    }
  return dname;
}
  
inline std::string ndfes::GetImgName( int i, int pad )
{
  std::stringstream sdir;
  sdir << "img" <<  std::setfill('0') << std::setw(pad) << i;
  return sdir.str();
}

inline std::string ndfes::GetInitName( int i, int pad )
{
  std::stringstream sdir;
  sdir << "init" <<  std::setfill('0') << std::setw(pad) << i;
  return sdir.str();
}



#endif
