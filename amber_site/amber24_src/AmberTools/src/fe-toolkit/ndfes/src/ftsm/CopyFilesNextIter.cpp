#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>


#include "CopyFilesNextIter.hpp"
#include "Disang.hpp"
#include "Utils.hpp"
#include "Filenames.hpp"

// #include <unistd.h>
// #include <cstdio>
// int main()
// {
//   int chdir_res = chdir("foo");
//   std::printf("chdir_res %i\n",chdir_res);
//   int symlnk_res = symlink("../it01/img02.rst7","a.rst7");
//   std::printf("symlink_res %i\n",chdir_res);
//   chdir_res = chdir("..");
//   std::printf("chdir_res %i\n",chdir_res);
//   return 0;
// }






void ndfes::CopyFilesNextIter
( std::size_t const ndim,
  std::size_t const nsim,
  std::string const disang,
  std::string const mdin,
  int const curit,
  std::string const nextdir,
  std::string const prefix,
  std::vector<std::string> const extra_prefixes,
  bool cp_nearest,
  bool multidir,
  int firstit,
  int pad,
  double const * rcs,
  double const * fcs )
{

  
  amber::Disang tdisang(disang);
  if ( not tdisang.ValidateTemplate() )
    {
      std::cerr << "Disang " << disang << " is not a template" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  std::vector<int> tidxs( tdisang.GetTemplateIdxs() );
  if ( tidxs.size() != ndim )
    {
      std::cerr << "CopyFilesNextIter expected " << ndim
		<< " reaction coordinates, but " << disang
		<< " has " << tidxs.size() << " templated restraints"
		<< std::endl;
      std::exit(EXIT_FAILURE);
    }

  
  
  std::vector<std::string> tmdin;
  {
    std::ifstream cin;
    cin.open( mdin.c_str() );
    if ( not cin )
      {
	std::cerr << "Could not open " << mdin << std::endl;
	std::exit(EXIT_FAILURE);
      }
    std::string line;
    while ( std::getline( cin, line ) )
      {
	tmdin.push_back( line );
      }

    bool found_disang = false;
    bool found_dumpave = false;
    for ( std::vector<std::string>::iterator p=tmdin.begin(),
	    pend=tmdin.end(); p!=pend; ++p )
      {
	if ( p->find("DISANG") != std::string::npos )
	  {
	    found_disang = true;
	  }
	if ( p->find("DUMPAVE") != std::string::npos )
	  {
	    found_dumpave = true;
	  }
      }
    if ( not found_disang )
      {
	std::cerr << "Template mdin " << mdin << " did not contain DISANG"
		  << std::endl;
	std::exit(EXIT_FAILURE);
      }
    if ( not found_dumpave )
      {
	std::cerr << "Template mdin " << mdin << " did not contain DUMPAVE"
		  << std::endl;
	std::exit(EXIT_FAILURE);
      }
  }



  
  std::vector<std::string> searches;
  for ( std::size_t i=0; i<extra_prefixes.size(); ++i )
    {
      if ( extra_prefixes[i].size() > 0 )
	{
	  searches.push_back( extra_prefixes[i] );
	}
    }
  if ( prefix.size() > 0 )
    {
      searches.push_back( prefix );
    }

  {
    std::vector<std::string> us;
    for ( std::vector<std::string>::iterator p=searches.begin(),
	    pend = searches.end(); p != pend; ++p )
      {
	bool found = false;
	for ( std::vector<std::string>::iterator q=us.begin(),
		qend = us.end(); q != qend; ++q )
	  {
	    if ( p->compare(*q) == 0 )
	      {
		found = true;
	      };
	  }
	if ( not found )
	  {
	    us.push_back( *p );
	  }
      }
    searches = us;
  }

  if ( searches.size() == 0 )
    {
      searches.push_back( "" );
    }



  
  std::vector<std::string> old_files;
  std::vector<double> old_rcs;
  
  if ( multidir or cp_nearest )
    {
      for ( std::vector<std::string>::iterator sdir=searches.begin(),
	      sdirend = searches.end(); sdir != sdirend; ++sdir )
	{
	  for ( int it=firstit; it < curit+1; ++it )
	    {
	      if ( (not multidir) and cp_nearest and (it != curit) )
		{
		  continue;
		}
	      
	      std::stringstream idir;
	      idir << ndfes::GetDirName(it,pad);

	      
	      if ( sdir->size() > 0 )
		{
		  idir << "/" << (*sdir);
		}

	      for ( std::size_t img=0; img<nsim; ++img )
		{
		  std::stringstream base;
		  base << idir.str() << "/" << ndfes::GetImgName(img+1,pad);
		  //base << idir.str() << "/img"
		  //     << std::setfill('0') << std::setw(3) << img+1;
		  std::string rst = base.str() + ".rst7";
		  std::string dis = base.str() + ".disang";

		  //std::cout << "rst " << rst << " dis " << dis << std::endl;
		  //std::cout << FileExists(rst) << " " << FileExists( dis )
		  //<< std::endl;
		  if ( FileExists( rst ) and FileExists( dis ) )
		    {
		      amber::Disang idisang( dis );
		      std::vector<double> rcs( idisang.GetCenters( tidxs ) );
		      old_files.push_back( rst );
		      for ( std::size_t dim=0; dim<ndim; ++dim )
			{
			  old_rcs.push_back( rcs[dim] );
			}
		    }
		}
	      
	    }
	}
    }

  std::string odir_path;
  {
    std::stringstream ss;
    ss << nextdir;
    if ( prefix.size() > 0 )
      {
	ss << "/" << prefix;
      }
    odir_path = ss.str();
  }


  
  std::string cdir;
  {
    std::stringstream ss;
    ss << ndfes::GetDirName(curit,pad);
    // if ( curit == 0 )
    //   {
    // 	ss << "init";
    //   }
    // else
    //   {
    // 	ss << "it"
    // 	   << std::setfill('0') << std::setw(3) << curit;
    //   }
    if ( prefix.size() > 0 )
      {
	ss << "/" << prefix;
      }
    cdir = ss.str();
  }

  

  {
    std::ofstream cout;
    cout.open( odir_path + "/sims.txt" );
    for ( std::size_t i=0; i<nsim; ++i )
      {
	cout << std::setw(3) << i+1
	     << " "
	     << std::fixed << std::setw(8) << std::setprecision(4)
	     << i/(nsim-1.);
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    cout << " "
		 << std::fixed << std::setw(15) << std::setprecision(8)
		 << rcs[dim+i*ndim];
	  }
	cout << "\n";
      }
  }


  for ( std::size_t i=0; i<nsim; ++i )
    {

      amber::Disang odisang(tdisang.FillTemplateValues(rcs+i*ndim,fcs+i*ndim));

      std::string old_rst;
      if ( (cp_nearest or multidir) )
	{
	  if ( old_files.size() == 0 )
	    {
	      std::cerr << "Failed to locate old restart files" << std::endl;
	      std::exit(EXIT_FAILURE);
	    };
	  
	  std::size_t iold=0;
	  double dmin=1.e+30;
	  for ( std::size_t k=0; k<old_files.size(); ++k )
	    {
	      double d = 0;
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  double dx = rcs[dim+i*ndim] - old_rcs[dim+k*ndim];
		  d += dx*dx;
		}
	      d = std::sqrt(d);
	      if ( d < dmin )
		{
		  dmin = d;
		  iold = k;
		}
	    }
	  if ( iold < old_files.size() )
	    {
	      old_rst = old_files[iold];
	    }
	  else
	    {
	      std::cerr << "Failed to locate restart for img " << i+1 << std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	}
      else
	{
	  std::stringstream ss;
	  ss << cdir << "/" << ndfes::GetImgName(i+1,pad) << ".rst7";
	  // ss << cdir << "/img"
	  //    << std::setfill('0') << std::setw(3) << i+1
	  //    << ".rst7";
	  old_rst = ss.str();
	}


      {
	std::stringstream fname;
	fname << odir_path << "/" << ndfes::GetImgName(i+1,pad) << ".disang";
	// fname << odir_path << "/img"
	//       << std::setfill('0') << std::setw(3) << i+1
	//       << ".disang";
	
	std::cout << "   Writing " << fname.str() << "\n";
	
	std::ofstream cout;
	cout.open( fname.str().c_str() );
	cout << odisang;
      }



      {
	std::stringstream fname;
	fname << odir_path << "/" << ndfes::GetImgName(i+1,pad) << ".mdin";
	// fname << odir_path << "/img"
	//       << std::setfill('0') << std::setw(3) << i+1
	//       << ".mdin";

	std::cout << "   Writing " << fname.str() << "\n";

	std::ofstream cout;
	cout.open( fname.str().c_str() );
	for ( std::vector<std::string>::iterator p=tmdin.begin(),
		pend=tmdin.end(); p!=pend; ++p )
	  {
	    if ( p->find("DISANG") != std::string::npos )
	      {
		// cout << "DISANG=img"
		//      << std::setfill('0') << std::setw(3) << i+1
		//      << ".disang\n";
		cout << "DISANG=" << ndfes::GetImgName(i+1,pad) << ".disang\n";
	      }
	    else if ( p->find("DUMPAVE") != std::string::npos )
	      {
		// cout << "DUMPAVE=img"
		//      << std::setfill('0') << std::setw(3) << i+1
		//      << ".dumpave\n";
		cout << "DUMPAVE=" << ndfes::GetImgName(i+1,pad) << ".dumpave\n";
	      }
	    else
	      {
		cout << *p << "\n";
	      }
	  }
	
      }
      

      std::string orst;
      {
	std::stringstream ss;
	// ss << odir_path << "/init"
	//    << std::setfill('0') << std::setw(3) << i+1
	//    << ".rst7";

	ss <<  odir_path << "/" << ndfes::GetInitName(i+1,pad) << ".rst7";
	
	orst = ss.str();
      }
      
      std::cout << "   Copying file "
		<< old_rst << " -> " << orst << "\n\n";
	  
      if ( FileExists( orst ) )
	{
	  std::remove( orst.c_str() );
	}
      else if ( FileIsLink( orst ) )
	{
	  unlink( orst.c_str() );
	}
      
      // if ( FileExists( orst ) )
      // 	{
      // 	  if ( FileIsLink( orst ) )
      // 	    {
      // 	      unlink( orst.c_str() );
      // 	    }
      // 	  else
      // 	    {
      // 	      std::remove( orst.c_str() );
      // 	    }
      // 	}

      bool created = RelativeSymLink(old_rst,orst);
      if ( not created )
	{
	  std::cerr << "Failed to create symbolic link: " << orst
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      
    }
  
}
