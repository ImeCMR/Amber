#ifndef _XMLIO_XmlUtilities_H_
#define _XMLIO_XmlUtilities_H_

#include <string>
#include <sstream>
#include <cassert>
#include "pugixml.hpp"
#include "XMLIO_Exceptions.hpp"
#include "XMLIO_Tokenize.hpp"


namespace XMLIO
{
  using namespace pugi;
}


namespace XMLIO
{

  std::string FamilyLine( XMLIO::xml_node const & node );


  inline XMLIO::xml_node ExtractNode( XMLIO::xml_node const & node, char const * f1 );

  inline XMLIO::xml_node ExtractNode( XMLIO::xml_node const & node, char const * f1, char const * f2 );

  inline XMLIO::xml_node ExtractNode( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3 );


  inline XMLIO::xml_node ExtractOptionalNode( XMLIO::xml_node const & node, char const * f1 );

  inline XMLIO::xml_node ExtractOptionalNode( XMLIO::xml_node const & node, char const * f1, char const * f2 );

  inline XMLIO::xml_node ExtractOptionalNode( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3 );



  template <class T>
  void ExtractValue( XMLIO::xml_node const & node, T & v );

  void ExtractValue( XMLIO::xml_node const & node, std::string & v );

  void ExtractValue( XMLIO::xml_node const & node, bool & v );

  template <class T>
  void ExtractValue( XMLIO::xml_node const & node, char const * f1, T & v );

  void ExtractValue( XMLIO::xml_node const & node, char const * f1, std::string & v );

  void ExtractValue( XMLIO::xml_node const & node, char const * f1, bool & v );

  template <class T>
  void ExtractValue( XMLIO::xml_node const & node, char const * f1, char const * f2, T & v );

  template <class T>
  void ExtractValue( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3, T & v );



  template <class T>
  void ExtractOptionalValue( XMLIO::xml_node const & node, char const * f1, T & v );

  void ExtractOptionalValue( XMLIO::xml_node const & node, char const * f1, std::string & v );

  void ExtractOptionalValue( XMLIO::xml_node const & node, char const * f1, bool & v );

  template <class T>
  void ExtractOptionalValue( XMLIO::xml_node const & node, char const * f1, char const * f2, T & v );

  template <class T>
  void ExtractOptionalValue( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3, T & v );




  template <class T>
  inline void ExtractToken( XMLIO::xml_node const & node, char const * f1, T & v );

  inline void ExtractToken( XMLIO::xml_node const & node, char const * f1, std::string & v );

  inline void ExtractToken( XMLIO::xml_node const & node, char const * f1, bool & v );

  template <class T>
  inline void ExtractToken( XMLIO::xml_node const & node, char const * f1, char const * f2, T & v );

  template <class T>
  inline void ExtractToken( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3, T & v );





  template <class T>
  inline void ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, T & v );

  inline void ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, std::string & v );

  inline void ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, bool & v );

  template <class T>
  inline void ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, char const * f2, T & v );

  template <class T>
  inline void ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3, T & v );






  template <class T>
  void ExtractAttribute( XMLIO::xml_node const & node, char const * f1, T & v );

  void ExtractAttribute( XMLIO::xml_node const & node, char const * f1, std::string & v );

  void ExtractAttribute( XMLIO::xml_node const & node, char const * f1, bool & v );

  template <class T>
  void ExtractAttribute( XMLIO::xml_node const & node, char const * f1, char const * f2, T & v );

  template <class T>
  void ExtractAttribute( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3, T & v );

  
}




inline XMLIO::xml_node XMLIO::ExtractNode( XMLIO::xml_node const & node, char const * f1 )
{
  assert( node );
  if ( ! node.child(f1) )
    {
      std::stringstream msg;
      msg << "XMLIO::ExtractNode couldn't extract child"
	  << " <" << f1 << "> at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    };
  return node.child(f1);
}


inline XMLIO::xml_node XMLIO::ExtractOptionalNode( XMLIO::xml_node const & node, char const * f1 )
{
  return node.child(f1);
}


inline XMLIO::xml_node XMLIO::ExtractNode( XMLIO::xml_node const & node, char const * f1, char const * f2 )
{
  assert( node );
  if ( ! node.child(f1) )
    {
      std::stringstream msg;
      msg << "XMLIO::ExtractNode couldn't extract child"
	  << " <" << f1 << "> at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    };
  return XMLIO::ExtractNode(node.child(f1),f2);
}


inline XMLIO::xml_node XMLIO::ExtractOptionalNode( XMLIO::xml_node const & node, char const * f1, char const * f2 )
{
  return node.child(f1).child(f2);
}


inline XMLIO::xml_node XMLIO::ExtractNode( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3 )
{
  assert( node );
  if ( ! node.child(f1) )
    {
      std::stringstream msg;
      msg << "XMLIO::ExtractNode couldn't extract child"
	  << " <" << f1 << "> at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    };
  return XMLIO::ExtractNode(node.child(f1),f2,f3);
}


inline XMLIO::xml_node XMLIO::ExtractOptionalNode( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3 )
{
  return node.child(f1).child(f2).child(f3);
}


template <class T>
void XMLIO::ExtractValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  T & v )
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
    }

  std::stringstream ss;
  ss << child.child_value();
  ss >> v;
  if ( ss.fail() )
    {
      ss.clear();
      ss.seekg(0, std::ios::beg);
      std::stringstream msg;
      msg << "XMLIO::ExtractValue couldn't convert"
	  << " '" << ss.str() << "' in field"
	  << " <" << f1 << "> at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    };
}


template <class T>
void XMLIO::ExtractValue
( XMLIO::xml_node const & node, 
  T & v )
{
  assert( node );

  std::stringstream ss;
  ss << node.child_value();
  ss >> v;
  if ( ss.fail() )
    {
      ss.clear();
      ss.seekg(0, std::ios::beg);
      std::stringstream msg;
      msg << "XMLIO::ExtractValue couldn't convert"
	  << " '" << ss.str() << "' at "
	  << XMLIO::FamilyLine(node);
      throw XMLIO::exception(msg);
    };
}



template <class T>
void XMLIO::ExtractOptionalValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  T & v )
{
  XMLIO::xml_node child = node.child(f1);
  if ( child )
    {
      std::stringstream ss;
      ss << child.child_value();
      T tmp;
      ss >> tmp;
      if ( ! ss.fail() ) v=tmp;
    };
}



template <class T>
void XMLIO::ExtractValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  char const * f2,
  T & v )
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
	  << " <" << f2 << "> at "
	  << XMLIO::FamilyLine(node) << "<" << f1 << ">" << "\n";
      msg << "\nbecause:\n" << e.what();
      throw XMLIO::exception(msg);
    };

  XMLIO::ExtractValue(child,f2,v);
}

template <class T>
void XMLIO::ExtractOptionalValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  char const * f2,
  T & v )
{
  XMLIO::xml_node child = node.child(f1);
  XMLIO::ExtractOptionalValue(child,f2,v);
}



template <class T>
void XMLIO::ExtractValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  char const * f2,
  char const * f3,
  T & v )
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
	  << " <" << f3 << "> at "
	  << XMLIO::FamilyLine(node) << "<" << f1 << ">" << "<" << f2 << ">" << "\n";
      msg << "\nbecause:\n" << e.what();
      throw XMLIO::exception(msg);
    };

  XMLIO::ExtractValue(child,f2,f3,v);
}

template <class T>
void XMLIO::ExtractOptionalValue
( XMLIO::xml_node const & node, 
  char const * f1, 
  char const * f2,
  char const * f3,
  T & v )
{
  XMLIO::xml_node child = node.child(f1);
  XMLIO::ExtractOptionalValue(child,f2,f3,v);
}




template <class T>
inline void XMLIO::ExtractToken( XMLIO::xml_node const & node, char const * f1, T & v )
{
  XMLIO::ExtractValue(node,f1,v);
}


template <class T>
inline void XMLIO::ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, T & v )
{
  XMLIO::ExtractOptionalValue(node,f1,v);
}


inline void XMLIO::ExtractToken( XMLIO::xml_node const & node, char const * f1, std::string & v )
{
  XMLIO::ExtractValue(node,f1,v);
  XMLIO::tokenize(v);
}



inline void XMLIO::ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, std::string & v )
{
  XMLIO::ExtractOptionalValue(node,f1,v);
  XMLIO::tokenize(v);
}


inline void XMLIO::ExtractToken( XMLIO::xml_node const & node, char const * f1, bool & v )
{
  XMLIO::ExtractValue(node,f1,v);
}


inline void XMLIO::ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, bool & v )
{
  XMLIO::ExtractOptionalValue(node,f1,v);
}

template <class T>
inline void XMLIO::ExtractToken( XMLIO::xml_node const & node, char const * f1, char const * f2, T & v )
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
      msg << "XMLIO::ExtractToken couldn't extract non-optional field"
	  << " <" << f2 << "> at "
	  << XMLIO::FamilyLine(node) << "<" << f1 << ">" << "\n";
      msg << "\nbecause:\n" << e.what();
      throw XMLIO::exception(msg);
    };

  XMLIO::ExtractToken(child,f2,v);
}


template <class T>
inline void XMLIO::ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, char const * f2, T & v )
{
  XMLIO::xml_node child = node.child(f1);
  XMLIO::ExtractOptionalToken(child,f2,v);
}




template <class T>
inline void XMLIO::ExtractToken( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3, T & v )
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
      msg << "XMLIO::ExtractToken couldn't extract non-optional field"
	  << " <" << f3 << "> at "
	  << XMLIO::FamilyLine(node) << "<" << f1 << ">" << "<" << f2 << ">" << "\n";
      msg << "\nbecause:\n" << e.what();
      throw XMLIO::exception(msg);
    };

  XMLIO::ExtractToken(child,f2,f3,v);
}


template <class T>
inline void XMLIO::ExtractOptionalToken( XMLIO::xml_node const & node, char const * f1, char const * f2, char const * f3, T & v )
{
  XMLIO::xml_node child = node.child(f1);
  XMLIO::ExtractOptionalToken(child,f2,f3,v);
}









template <class T>
void XMLIO::ExtractAttribute
( XMLIO::xml_node const & node, 
  char const * f1, 
  T & v )
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
      std::stringstream ss;
      ss << node.attribute(f1).value();
      ss >> v;
      if ( ss.fail() )
	{
	  ss.clear();
	  ss.seekg(0, std::ios::beg);
	  std::stringstream msg;
	  msg << "XMLIO::ExtractAttribute couldn't convert"
	      << " '" << ss.str() << "' in attr"
	      << " '" << f1 << "' at "
	      << XMLIO::FamilyLine(node);
	  throw XMLIO::exception(msg);
	};
    };
}




template <class T>
void XMLIO::ExtractAttribute
( XMLIO::xml_node const & node, 
  char const * f1, 
  char const * f2,
  T & v )
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
      msg << "XMLIO::ExtractAttribute couldn't extract non-optional attr"
	  << " '" << f2 << "' at "
	  << XMLIO::FamilyLine(node) << "<" << f1 << ">" << "\n";
      msg << "\nbecause:\n" << e.what();
      throw XMLIO::exception(msg);
    };

  XMLIO::ExtractAttribute(child,f2,v);
}



template <class T>
void XMLIO::ExtractAttribute
( XMLIO::xml_node const & node, 
  char const * f1, 
  char const * f2,
  char const * f3,
  T & v )
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
      msg << "XMLIO::ExtractAttribute couldn't extract non-optional field"
	  << " <" << f3 << "> at "
	  << XMLIO::FamilyLine(node) << "<" << f1 << ">" << "<" << f2 << ">" << "\n";
      msg << "\nbecause:\n" << e.what();
      throw XMLIO::exception(msg);
    };

  XMLIO::ExtractAttribute(child,f2,f3,v);
}


#endif
