#include "XMLIO_XmlUtilities.hpp"
#include <cstring>



std::string XMLIO::FamilyLine
( XMLIO::xml_node const & node )
{
  assert( node );

  std::string s;
  for ( XMLIO::xml_node me = node; me; me = me.parent() )
    {
      if ( strlen(me.name()) == 0 ) break;
      s.insert(0,1,'>');
      s.insert(0,me.name());
      s.insert(0,1,'<');
    };
  return s;
}



void XMLIO::ExtractValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  std::string & v )
{
  assert( node );

  XMLIO::xml_node child;
  try
    {
      child = XMLIO::ExtractNode(node,f1);
    }
  catch ( XMLIO::exception const & e )
    {
      std::stringstream msg;
      msg << "XMLIO::ExtractValue couldn't extract non-optional field"
	  << " <" << f1 << "> at "
	  << XMLIO::FamilyLine(node) << "<" << f1 << ">" << "\n";
      msg << "\nbecause:\n" << e.what();
      throw XMLIO::exception(msg);
    };

  v = child.child_value();
  //std::cout << "string value = '" << v << "' for node '" << child.name() << "'\n";
}


void XMLIO::ExtractValue
( XMLIO::xml_node const & node, 
  std::string & v )
{
  assert( node );
  v = node.child_value();
}



void XMLIO::ExtractOptionalValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  std::string & v )
{
  XMLIO::xml_node child = node.child(f1);
  if ( child )
    v = child.child_value();
}



void XMLIO::ExtractValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  bool & v )
{
  assert( node );

  XMLIO::xml_node child;
  try
    {
      child = XMLIO::ExtractNode(node,f1);
    }
  catch ( XMLIO::exception const & e )
    {
      std::stringstream msg;
      msg << "XMLIO::ExtractValue couldn't extract non-optional field"
	  << " <" << f1 << "> at "
	  << XMLIO::FamilyLine(node) << "<" << f1 << ">" << "\n";
      msg << "\nbecause:\n" << e.what();
      throw XMLIO::exception(msg);
    };


  std::string strbool = XMLIO::uppercase( child.child_value() );
  strbool = XMLIO::tokenize(strbool);
  if ( strbool[0] == 'T' or strbool[0] == '1' ) 
    { v = true; }
  else if ( strbool[0] == 'F' or strbool[0] == '0' ) 
    { v = false; }
  else 
    { 
      std::stringstream msg;
      msg << "XMLIO::ExtractValue couldn't convert bool from"
	  << " '" << strbool << "' in field"
	  << " <" << f1 << "> at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    };
  
}




void XMLIO::ExtractValue
( XMLIO::xml_node const & node, 
  bool & v )
{
  assert( node );

  std::string strbool = XMLIO::uppercase( node.child_value() );
  strbool = XMLIO::tokenize(strbool);
  if ( strbool[0] == 'T' or strbool[0] == '1' ) 
    { v = true; }
  else if ( strbool[0] == 'F' or strbool[0] == '0' ) 
    { v = false; }
  else 
    { 
      std::stringstream msg;
      msg << "XMLIO::ExtractValue couldn't convert bool from"
	  << " '" << strbool << "' at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    };
  
}



void XMLIO::ExtractOptionalValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  bool & v )
{
  XMLIO::xml_node child = node.child(f1);
  if ( child )
    {
      std::string strbool = XMLIO::uppercase( child.child_value() );
      strbool = XMLIO::tokenize(strbool);
      if ( strbool[0] == 'T' or strbool[0] == '1' ) 
	{ v = true; }
      else if ( strbool[0] == 'F' or strbool[0] == '0' ) 
	{ v = false; };
    };
}






void XMLIO::ExtractAttribute
( XMLIO::xml_node const & node, 
  char const * f1, 
  std::string & v )
{
  assert( node );

  if ( ! node.attribute(f1) )
    {
      std::stringstream msg;
      msg << "XMLIO::ExtractAttribute couldn't extract non-optional attr"
	  << " '" << f1 << "' at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    }
  else
    {
      v = node.attribute(f1).value();
      v = XMLIO::tokenize(v);
    };
}



void XMLIO::ExtractAttribute
( XMLIO::xml_node const & node, 
  char const * f1, 
  bool & v )
{
  assert( node );

  if ( ! node.attribute(f1) )
    {
      std::stringstream msg;
      msg << "XMLIO::ExtractAttribute couldn't extract non-optional attr"
	  << " '" << f1 << "' at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    }
  else
    {
      std::string strbool = XMLIO::uppercase( node.attribute(f1).value() );
      strbool = XMLIO::tokenize(strbool);
      if ( strbool[0] == 'T' or strbool[0] == '1' ) 
	{ v = true; }
      else if ( strbool[0] == 'F' or strbool[0] == '0' ) 
	{ v = false; }
      else 
	{ 
	  std::stringstream msg;
	  msg << "XMLIO::ExtractAttribute couldn't convert bool from"
	      << " '" << strbool << "' in attr"
	      << " '" << f1 << "' at "
	  << XMLIO::FamilyLine(node);
	  throw XMLIO::exception(msg);
	};
      
    };
}
