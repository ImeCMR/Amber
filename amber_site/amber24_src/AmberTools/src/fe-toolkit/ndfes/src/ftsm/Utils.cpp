#include <sys/types.h>
#include <sys/stat.h>
// # include <linux/limits.h>
#include <unistd.h>

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <vector>

#include "Utils.hpp"



namespace
{

  std::string FindAndReplace
  ( std::string original, 
    std::string const searchString, 
    std::string const replaceString )
  {
    std::string::size_type pos = 0;
    while ( (pos = original.find(searchString, pos)) != std::string::npos ) 
      {
	original.replace( pos, searchString.size(), replaceString );
	pos++;
      };
    
    return original;
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
  
  std::vector<std::string> SplitPath( std::string line )
  {
    //std::cout << "SplitPath inp " << line << "\n";
    //std::cout << "          out " << FindAndReplace(line,"//","/") << "\n";
    std::vector<std::string> parts(Split( FindAndReplace(line,"//","/") , '/'));
    while ( true )
      {
	if ( parts.size() > 1 )
	  {
	    if ( parts[0].compare(".") == 0 )
	      {
		parts.erase(parts.begin());
	      }
	    else
	      {
		break;
	      }
	  }
	else
	  {
	    break;
	  };
      };
    return parts;
  }

  std::string JoinPath( std::vector<std::string> paths )
  {
    std::stringstream ss;
    std::size_t npath = paths.size();
    for ( std::size_t i=0; i<npath; ++i )
      {
	if ( i > 0 )
	  {
	    ss << "/";
	  }
	ss << paths[i];
      }
    return ss.str();
  }
  
}





bool DirExists( std::string dname )
{
  bool exists = false;
  
  struct stat info;

  if( stat( dname.c_str(), &info ) != 0 )
    {
      exists = false;
      //printf( "cannot access %s\n", pathname );
    }
  else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
    {
      exists = true;
      //printf( "%s is a directory\n", pathname );
    }
  else
    {
      exists = false;
      //printf( "%s is no directory\n", pathname );
    }

  return exists;
}


bool FileExists( std::string fname )
{
  struct stat buffer;   
  return (stat (fname.c_str(), &buffer) == 0); 
}

bool FileIsLink( std::string fname )
{
  bool islink = false;
  struct stat buffer;
  int x = lstat (fname.c_str(), &buffer);
  if ( x == 0 )
    {
      if (S_ISLNK(buffer.st_mode))
	{
	  islink = true;
	}
    }
  return islink; 
}

void CopyFile( std::string sfile, std::string dfile )
{
  std::ifstream src(sfile.c_str(), std::ios::binary);
  std::ofstream dst(dfile.c_str(), std::ios::binary);
  dst << src.rdbuf();
}

std::string GetCwd()
{
  std::size_t const myPATH_MAX = 4096*4;
  char buff[myPATH_MAX];
  char * res = getcwd( buff, myPATH_MAX );
  if ( res == NULL )
    {
      std::cerr << "Failure calling getcwd"
		<< std::endl;
      std::exit(EXIT_FAILURE);
    }
  std::string cwd( buff );
  return cwd;
}


bool RelativeSymLink( std::string oldfile, std::string newlink )
{

  if ( not FileExists( oldfile ) )
    {
      std::cerr << "Cannot create symbolic link "
		<< newlink << " to file " << oldfile
		<< " because the file does not exist"
		<< std::endl;
      return false;
    }
  
  if ( FileExists( newlink ) )
    {
      std::cerr << "Cannot create symbolic link "
		<< newlink << " to file " << oldfile
		<< " because the link already exists as file";
      std::cerr << std::endl;
      return false;
    }
  else if ( FileIsLink( newlink ) )
    {
      std::cerr << "Cannot create symbolic link "
		<< newlink << " to file " << oldfile
		<< " because the link already exists as link";
    }


  std::vector<std::string> oldpath( SplitPath(oldfile) );
  std::vector<std::string> newpath( SplitPath(newlink) );
  std::string cwd( GetCwd() );

  std::vector<std::string> tpath( newpath.begin(), newpath.end()-1 );
  std::string newdir( JoinPath( tpath ) );

  //std::cout << "newdir " << newdir << " tpath.size " << tpath.size() << "\n";
  
  if ( newdir.size() > 0 )
    {
      if ( not DirExists( newdir ) )
	{
	  std::cerr << "Cannot create symbolic link "
		    << newlink << " to file " << oldfile
		    << " because directory " << newdir
		    << " does not exist"
		    << std::endl;
	  return false;
	}
    }

  std::vector<std::string> rpath( oldpath );
  std::string dotdot("..");
  rpath.insert(rpath.begin(),tpath.size(),dotdot);
  std::string rfile( JoinPath( rpath ) );

  int res = 0;
  if ( newdir.size() > 0 )
    {
      res = chdir( newdir.c_str() );
    }
  if ( res != 0 )
    {
      std::cerr << "Cannot create symbolic link "
		<< newlink << " to file " << oldfile
		<< " because chdir to '" << newdir
		<< "' failed"
		<< std::endl;
      return false;
    }

  //std::cout << "link " << rfile << " " << newpath.back();
  
  res = symlink( rfile.c_str(), newpath.back().c_str() );

  if ( res != 0 )
    {
      std::cerr << "Cannot create symbolic link "
		<< newlink << " to file " << oldfile
		<< " because 'symlink(\"" << rfile
		<< "\",\"" << newpath.back()
		<< "\") failed"
		<< std::endl;
      return false;
    }
  
  if ( newdir.size() > 0 )
    {
      res = chdir( cwd.c_str() );
    }
  if ( res != 0 )
    {
      std::cerr << "Cannot create symbolic link "
		<< newlink << " to file " << oldfile
		<< " because 'chdir \"" << cwd << "\" failed"
		<< std::endl;
      return false;
    }

  return true;
}



double SwitchOff( double const r , double const rlo , double const rhi )
{
  double s = 1.;
  if ( r >= rhi )      s = 0.;
  else if ( r <= rlo ) s = 1.;
  else
    {
      double const u = (rhi-r)/(rhi-rlo);
      double const u3 = u*u*u;
      double const u4 = u3*u;
      double const u5 = u4*u;
      s = 10.*u3 - 15.*u4 + 6.*u5;
    };
  return s;
}

double SwitchOn( double const r , double const rlo , double const rhi )
{
  return 1-SwitchOff(r,rlo,rhi);
}
