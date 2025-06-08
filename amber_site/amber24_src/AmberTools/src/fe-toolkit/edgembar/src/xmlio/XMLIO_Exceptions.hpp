#ifndef _XMLIO_Exceptions_H_
#define _XMLIO_Exceptions_H_

#include <string>
#include <sstream>

namespace XMLIO
{
  
  class exception : public std::exception
  {
  public:
    exception( char const * msg ) : _what(msg) {};
    exception( std::string const msg ) : _what(msg) {};
    exception( std::stringstream const & msg ) : _what(msg.str()) {};
    virtual ~exception() throw() {};
    virtual void rethrow() const { throw *this; };
    char const * what() const throw() { return _what.c_str(); };
  protected:
    std::string _what;
  };


  class cant_construct : public XMLIO::exception
  {
  public:
    cant_construct( char const * msg ) : XMLIO::exception(msg) {};
    cant_construct( std::string const msg ) : XMLIO::exception(msg) {};
    cant_construct( std::stringstream const & msg ) : XMLIO::exception(msg) {};
    virtual ~cant_construct() throw() {};
    virtual void rethrow() const { throw *this; };
  };


  class cant_openfile : public XMLIO::exception
  {
  public:
    cant_openfile( char const * msg ) : XMLIO::exception(msg) {};
    cant_openfile( std::string const msg ) : XMLIO::exception(msg) {};
    cant_openfile( std::stringstream const & msg ) : XMLIO::exception(msg) {};
    virtual ~cant_openfile() throw() {};
    virtual void rethrow() const { throw *this; };
  };

}

#endif
