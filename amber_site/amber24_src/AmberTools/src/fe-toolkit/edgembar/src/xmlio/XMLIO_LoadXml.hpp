#ifndef _XMLIO_LoadXml_H_
#define _XMLIO_LoadXml_H_

#include "XMLIO_XmlUtilities.hpp"
#include <iostream>
#include <fstream>

namespace XMLIO
{
  
  inline void LoadXml( std::istream & in , XMLIO::xml_document & xml )
  {
    XMLIO::xml_parse_result xml_is_loaded = xml.load(in);
    if ( ! xml_is_loaded )
      {
	std::stringstream msg;
	msg << "GetXmlNode failed to load xml file.  Description: " 
	    << xml_is_loaded.description();
	throw XMLIO::exception(msg);
      };
  }
  
  
  inline void LoadXml( std::string const & file , XMLIO::xml_document & xml )
  {
    std::ifstream str(file.c_str());
    XMLIO::LoadXml(str,xml);
  }
  
  
}


#endif
