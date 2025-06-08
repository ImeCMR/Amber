#ifndef _XMLIO_Tokenize_hpp_
#define _XMLIO_Tokenize_hpp_

#include <string>
#include <algorithm>

namespace XMLIO
{
  std::string tokenize( std::string str );
  std::string tokenize( char const *  car );
  std::string uppercase( std::string str );
}

inline std::string XMLIO::tokenize( std::string str )
{
  std::size_t const upper = str.find_last_not_of( " \t\f\v\n\r" ) + 1;
  std::size_t const lower = str.find_first_not_of( " \t\f\v\n\r" );
  str.erase( upper );
  str.erase( 0 , lower );
  return str;
}

inline std::string XMLIO::tokenize( char const *  car )
{
  return XMLIO::tokenize( std::string(car) );
}

inline std::string XMLIO::uppercase( std::string str )
{
  std::transform( str.begin(), str.end(), str.begin(), ::toupper );
  return str;
}

#endif

