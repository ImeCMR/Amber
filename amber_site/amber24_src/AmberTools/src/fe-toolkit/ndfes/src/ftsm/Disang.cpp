#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <iomanip>
#include <vector>

#include "Disang.hpp"

//
// internal linkage utilities
//
namespace
{
  double const DEG_PER_RAD = 180. / M_PI;
  double const RAD_PER_DEG = M_PI / 180.;

  std::string TrimAfter( char const * c, std::string str );
  
  std::string Strip( std::string str );
  
  std::vector<std::string> Split( std::string line );
  
  std::vector<std::string> Split( std::string line, char const delim );
  
  std::vector< std::vector<std::string> >
  SplitWords
  ( std::string line, char const delim );

  void KeysAndValues
  ( std::string line, char const delim,
    std::vector<std::string> & tokens,
    std::vector< std::vector<std::string> > & values );
  
  void ProcessRestraintLine
  ( std::string line,
    std::vector<std::string> & keys,
    std::vector< std::vector<std::string> > & values );  
}



namespace
{

  std::string TrimAfter( char const * c, std::string str )
  {
    size_t pos = str.find( c );
    if ( pos != std::string::npos )
      {
	str.erase( pos );
      }
    return str;
  }
  
  std::string Strip( std::string str )
  {
    std::size_t const upper = str.find_last_not_of( " \t\f\v\n\r" ) + 1;
    std::size_t const lower = str.find_first_not_of( " \t\f\v\n\r" );
    str.erase( upper );
    str.erase( 0 , lower );
    return str;
  }
  
  
  
  std::vector<std::string> Split( std::string line )
  {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::copy( std::istream_iterator<std::string>(iss),
	       std::istream_iterator<std::string>(),
	       std::back_inserter(tokens) );
    return tokens;
  }
  
  std::vector<std::string> Split( std::string line, char const delim )
  {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    while ( std::getline(iss,line,delim) )
      {
	tokens.push_back( line );
      }
    return tokens;
  }
  
  std::vector< std::vector<std::string> > SplitWords( std::string line, char const delim )
  {
    std::vector< std::vector<std::string> > swords;
    std::vector<std::string> secs( Split(line,delim) );
    for ( std::vector<std::string>::iterator p=secs.begin(); p != secs.end(); ++p )
      {
	swords.push_back( Split(*p) );
      }
    return swords;
  }
  
  
  void KeysAndValues( std::string line, char const delim,
		      std::vector<std::string> & tokens,
		      std::vector< std::vector<std::string> > & values )
  {
    std::vector< std::vector<std::string> > swords( SplitWords(line,delim) );
    tokens.resize(0);
    values.resize(0);
    for ( std::size_t i=0; i<swords.size(); ++i )
      {
	if ( i+1 < swords.size() )
	  {
	    tokens.push_back( swords[i].back() );
	  }
	if ( i > 0 )
	  {
	    std::vector<std::string>::iterator pend = swords[i].end();
	    if ( i+1 < swords.size() )
	      {
		--pend;
	      }
	    values.push_back( std::vector<std::string>( swords[i].begin(), pend ) );
	  }
      }
  }
  
  void ProcessRestraintLine
  ( std::string line,
    std::vector<std::string> & keys,
    std::vector< std::vector<std::string> > & values )
  {
    std::replace( line.begin(), line.end(), '\n', ' ');
    std::replace( line.begin(), line.end(), ',',  ' ');
    std::replace( line.begin(), line.end(), '/',  ' ');
    
    std::vector<std::string> tokens( Split(line) );
    
    std::stringstream ss;
    for ( std::vector<std::string>::iterator
	    p=tokens.begin(); p != tokens.end(); )
      {
	if ( (p->compare("&rst") == 0)
	     or (p->compare("/") == 0)
	     or (p->compare("&end") == 0) )
	  {
	    p = tokens.erase( p );
	  }
	else
	  {
	    ss << " " << *p;
	    ++p;
	  }
      }
    line = ss.str();
    
    //std::vector<std::string> keys;
    //std::vector< std::vector<std::string> > values;
    
    keys.resize(0);
    values.resize(0);
    
    KeysAndValues( line, '=', keys, values );
    
    // for ( std::size_t i=0; i<keys.size(); ++i )
    //   {
    //     std::cout << keys[i] << " : ";
    //     for ( std::vector<std::string>::iterator
    // 	      p = values[i].begin(); p != values[i].end(); ++p )
    // 	{
    // 	  std::cout << *p << " ";
    // 	}
    //     std::cout << "\n";
    //   }
  }
  
}




void amber::Restraint::Write( std::ostream & cout ) const
{
  std::string t = "    ";
  cout << "&rst\n";

#define FMTF std::fixed << std::setw(14) << std::setprecision(8)
#define FMTS std::fixed << std::setw(5) << std::setprecision(1)
#define FMTW std::setw(14)

  cout << t << "iat = ";
  for ( std::vector<int>::const_iterator
	  p=iat.begin(); p!=iat.end(); ++p )
    {
      if ( p != iat.begin() )
	{
	  cout << ",";
	}
      cout << std::setw(5) << *p;
    }
  cout << "\n";

  if ( rstwt.size() > 0 )
    {
      cout << t << "rstwt = ";
      for ( std::vector<double>::const_iterator
	      p=rstwt.begin(); p!=rstwt.end(); ++p )
	{
	  if ( p != rstwt.begin() )
	    {
	      cout << ",";
	    }
	  cout << FMTS << *p;
	}
      cout << "\n";
    }
  
  cout << t << "r1  = " << FMTF << r1
       << ", r4  = " << FMTF << r4
       << "\n";

  if ( tidx > 0 )
    {
      std::ostringstream ss;
      ss << "RC" << tidx;
      cout << t << "r2  = "   << FMTW << ss.str()
	   << ", r3  = " << FMTF << ss.str()
	   << "\n";
    }
  else
    {
      cout << t << "r2  = "   << FMTF << r2
	   << ", r3  = " << FMTF << r3
	   << "\n";
    }
  
  double f = 1.;
  if ( isangle )
    {
      f = RAD_PER_DEG * RAD_PER_DEG;
    }

  
  cout << t << "rk2 = " << FMTF << rk2 * f
       << ", rk3 = " << FMTF << rk3 * f
       << "\n";

  for ( std::vector<std::string>::const_iterator
	  p = extra.begin(); p != extra.end(); ++p )
    {
      cout << t << *p << "\n";
    }
  
  cout << "/\n";

#undef FMTF
#undef FMTS
#undef FMTW
}


std::ostream & operator<<( std::ostream & cout, amber::Restraint const & res )
{
  res.Write(cout);
  return cout;
}



amber::Restraint::Restraint( std::string line )
  : r1(-1000.),
    r2(0.0),
    r3(0.0),
    r4(1000.),
    rk2(0.0),
    rk3(0.0),
    tidx(0),
    isr12(false),
    isangle(false),
    isdihed(false)
{
  std::vector<std::string> keys;
  std::vector< std::vector<std::string> > vals;
  ::ProcessRestraintLine( line, keys, vals );

  
  for ( std::size_t i=0; i<keys.size(); ++i )
    {
      if ( keys[i].compare("iat") == 0 )
	{
	  for ( std::vector<std::string>::iterator
		  p = vals[i].begin(); p != vals[i].end(); ++p )
	    {
	      iat.push_back( std::atoi( p->c_str() ) );
	    }
	}
      else if ( keys[i].compare("r1") == 0 )
	{
	  r1 = std::atof(vals[i][0].c_str());
	}
      else if ( keys[i].compare("r2") == 0 )
	{
	  if ( vals[i][0].compare("RC") > 0 )
	    {
	      int idx = std::atoi( vals[i][0].substr(2).c_str() );
	      if ( tidx > 0 )
		{
		  if ( idx != tidx )
		    {
		      std::cerr << "Inconsistent template variable "
				<< vals[i][0] << " when processing: "
				<< line << "\n";
		      std::exit(EXIT_FAILURE);
		    }
		}
	      tidx = idx;
	    }
	  else if ( tidx > 0 )
	    {
	      std::cerr << "Expected r2 template index RC" << tidx
			<< " when processing: "
			<< line << "\n";
	      std::exit(EXIT_FAILURE);
	    }
	  else
	    {
	      r2 = std::atof( vals[i][0].c_str() );
	    }
	}
      else if ( keys[i].compare("r3") == 0 )
	{
	  if ( vals[i][0].compare("RC") > 0 )
	    {
	      int idx = std::atoi( vals[i][0].substr(2).c_str() );
	      if ( tidx > 0 )
		{
		  if ( idx != tidx )
		    {
		      std::cerr << "Inconsistent template variable "
				<< vals[i][0] << " when processing: "
				<< line << "\n";
		      std::exit(EXIT_FAILURE);
		    }
		}
	      tidx = idx;
	    }
	  else if ( tidx > 0 )
	    {
	      std::cerr << "Expected r3 template index RC" << tidx
			<< " when processing: "
			<< line << "\n";
	      std::exit(EXIT_FAILURE);
	    }
	  else
	    {
	      r3 = std::atof( vals[i][0].c_str() );
	    }
	}
      else if ( keys[i].compare("r4") == 0 )
	{
	  r4 = std::atof(vals[i][0].c_str());
	}
      else if ( keys[i].compare("rk2") == 0 )
	{
	  rk2 = std::atof(vals[i][0].c_str());
	}
      else if ( keys[i].compare("rk3") == 0 )
	{
	  rk3 = std::atof(vals[i][0].c_str());
	}
      else if ( keys[i].compare("rstwt") == 0 )
	{
	  for ( std::vector<std::string>::iterator
		  p = vals[i].begin(); p != vals[i].end(); ++p )
	    {
	      rstwt.push_back( std::atof( p->c_str() ) );
	    }
	}
      else
	{
	  std::ostringstream ss;
	  ss << keys[i] << " =";
	  for ( std::vector<std::string>::iterator
		  p = vals[i].begin(); p != vals[i].end(); ++p )
	    {
	      ss << " " << *p;
	    }
	  extra.push_back( ss.str() );
	}
	  
    }

  if ( iat.size() > 2 )
    {
      if ( not (rstwt.size() > 1) )
	{
	  isangle = true;
	}
    }
  if ( iat.size() == 4 )
    {
      if ( rstwt.size() > 1 )
	{
	  isr12 = true;
	}
      else
	{
	  isdihed = true;
	}
    }

  if ( isangle )
    {
      double f = DEG_PER_RAD * DEG_PER_RAD;
      rk2 *= f;
      rk3 *= f;
    }

}



amber::Disang::Disang( std::string fname )
{
  Read(fname);
}

void amber::Disang::Read( std::string fname )
{
  std::ifstream cin;
  cin.open(fname.c_str());
  if ( ! cin )
    {
      std::cerr << "Could not open " << fname << std::endl;
      std::exit(EXIT_FAILURE);
    }
  Read(cin);
}
    
void amber::Disang::Read( std::istream & cin )
{
  std::string line;
  std::string cmd;
  bool reading = false;
  while ( std::getline(cin,line) )
    {
      line = Strip(TrimAfter("!",line));
      if ( line.size() > 0 )
	{
	  if ( line.compare(0,4,"&rst") == 0 )
	    {
	      if ( reading and cmd.size() > 0 )
		{
		  mRestraints.push_back( amber::Restraint( cmd ) );
		  cmd = "";
		}
	      bool hasend = false;
	      if ( line.find("&end") != std::string::npos )
		{
		  hasend = true;
		}
	      if ( line.find("/") != std::string::npos )
		{
		  hasend = true;
		}
	      if ( hasend )
		{
		  reading = false;
		  mRestraints.push_back( amber::Restraint( line ) );
		  cmd = "";
		}
	      else
		{
		  reading = true;
		  cmd = line;
		};
	    }
	  else if ( reading )
	    {
	      cmd = cmd + " " + line;
	      bool hasend = false;
	      if ( line.find("&end") != std::string::npos )
		{
		  hasend = true;
		}
	      if ( line.find("/") != std::string::npos )
		{
		  hasend = true;
		}
	      if ( hasend )
		{
		  reading = false;
		  mRestraints.push_back( amber::Restraint( cmd ) );
		  cmd = "";
		}
	    }
	}
    }
  if ( cmd.size() > 0 )
    {
      mRestraints.push_back( amber::Restraint( cmd ) );
    }
}
    
void amber::Disang::Write( std::string fname ) const
{
  std::ofstream cout;
  cout.open(fname.c_str());
  Write(cout);
}

void amber::Disang::Write( std::ostream & cout ) const
{
  for ( std::size_t i=0; i<mRestraints.size(); ++i )
    {
      mRestraints[i].Write(cout);
      cout << "\n";
    }
}

bool amber::Disang::ValidateTemplate() const
{
  bool ok = true;
  int nactive = 0;
  for ( std::size_t i=0; i<mRestraints.size(); ++i )
    {
      int tidx = mRestraints[i].tidx;
      if ( tidx > 0 )
	{
	  ++nactive;
	  if ( tidx != nactive )
	    {
	      ok = false;
	    }
	}
    }
  if ( nactive == 0 )
    {
      ok = false;
    }
  return ok;
}
    
std::vector<int> amber::Disang::GetTemplateIdxs() const
{
  std::vector<int> tidxs;
  for ( std::size_t i=0; i<mRestraints.size(); ++i )
    {
      int tidx = mRestraints[i].tidx;
      if ( tidx > 0 )
	{
	  tidxs.push_back( i );
	}
    };
  return tidxs;
}

std::vector<bool> amber::Disang::GetTemplateIdxIsAngle() const
{
  std::vector<bool> isang;
  for ( std::size_t i=0; i<mRestraints.size(); ++i )
    {
      int tidx = mRestraints[i].tidx;
      if ( tidx > 0 )
	{
	  isang.push_back( mRestraints[i].isangle );
	}
    }
  return isang;
}


std::vector<double> amber::Disang::GetTemplateForceConsts() const
{
  std::vector<double> fcs;
  for ( std::size_t i=0; i<mRestraints.size(); ++i )
    {
      int tidx = mRestraints[i].tidx;
      if ( tidx > 0 )
	{
	  fcs.push_back( mRestraints[i].rk2 );
	}
    }
  return fcs;
}


amber::Disang amber::Disang::FillTemplateValues
( double const * rcs,
  double const * fcs ) const
{
  amber::Disang odis;
  for ( std::size_t i=0; i<mRestraints.size(); ++i )
    {
      amber::Restraint res( mRestraints[i] );
      int tidx = mRestraints[i].tidx;
      if ( tidx > 0 )
	{
	  res.r2 = rcs[tidx-1];
	  res.r3 = rcs[tidx-1];
	  res.rk2 = fcs[tidx-1];
	  res.rk3 = fcs[tidx-1];
	  res.tidx = 0;
	}
      odis.mRestraints.push_back( res );
    }
  return odis;
}


std::vector<double> amber::Disang::GetCenters
( std::vector<int> const & idxs ) const
{
  std::vector<double> rcs( idxs.size(), 0 );
  for ( std::size_t i=0; i<idxs.size(); ++i )
    {
      if ( idxs[i] < 0 or (std::size_t)idxs[i] > mRestraints.size() )
	{
	  std::cerr << "amber::Disang::GetCenters could not read "
		    << "restraint " << idxs[i] << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      rcs[i] = mRestraints[idxs[i]].r2;
    }
  return rcs;
}

std::vector<double> amber::Disang::GetForceConsts
( std::vector<int> const & idxs ) const
{
  std::vector<double> fcs( idxs.size(), 0 );
  for ( std::size_t i=0; i<idxs.size(); ++i )
    {
      if ( idxs[i] < 0 or (std::size_t)idxs[i] > mRestraints.size() )
	{
	  std::cerr << "amber::Disang::GetForceConsts could not read "
		    << "restraint " << idxs[i] << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      fcs[i] = mRestraints[idxs[i]].rk2;
    }
  return fcs;
}


std::ostream & operator<<( std::ostream & cout, amber::Disang const & dis )
{
  dis.Write(cout);
  return cout;
}
