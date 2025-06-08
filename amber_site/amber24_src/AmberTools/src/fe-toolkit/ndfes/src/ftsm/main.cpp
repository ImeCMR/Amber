
//
// Compile with
//
// g++ main.cpp $(python3-config --includes --libs)
//

#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <iterator>

#include "GetBeta.hpp"
#include "Options.hpp"
#include "FES.hpp"
#include "Utils.hpp"
#include "Smoothing.hpp"
#include "PCurve.hpp"
#include "PathData.hpp"
#include "Convergence.hpp"
#include "PeriodicUtils.hpp"
#include "MeshUtils.hpp"
#include "Tube.hpp"
#include "Disang.hpp"
#include "CopyFilesNextIter.hpp"
#include "Filenames.hpp"

#ifdef WITH_OPENMP
#include <omp.h>
#endif

bool IsFloat( std::string myString )
{
  std::istringstream iss(myString);
  double f;
  iss >> std::noskipws >> f; // noskipws considers leading whitespace invalid
  // Check the entire string was consumed and if either failbit or badbit is set
  return iss.eof() && !iss.fail(); 
}


template <typename T>
double CalcAvg( std::vector<T> X )
{
  double s = 0;
  for ( std::size_t i=0; i<X.size(); ++i )
    {
      s += X[i];
    }
  return s / X.size();
}


std::vector<double> CalcAvg( std::size_t const ndim, std::size_t nsim, std::vector<double> X )
{
  std::vector<double> avg(ndim,0);
  for ( std::size_t i=0; i<nsim; ++i )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  avg[dim] += X[dim+i*ndim] / nsim;
	};
    };
  return avg;
}


std::vector<double> ExpandAvg( std::size_t const ndim, std::size_t nsim, std::vector<double> avg )
{
  std::vector<double> data(ndim*nsim,0);
  for ( std::size_t i=0; i<nsim; ++i )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  data[dim+i*ndim] = avg[dim];
	};
    };
  return data;
}





void ReadPathFromDumpaves
( std::vector<int> idxs,
  std::vector< std::string > const & idumps,
  int const dumpmode,
  std::vector<double> & nobs,
  std::vector<double> & vals,
  std::vector<double> & errs )
{

  int ndumps = idumps.size();
  int ndim = idxs.size();
  int maxcol = ndim + 1;
  
  if ( dumpmode == 0 )
    {
      std::cerr << "dumpmode=0 : unclear how to interpret"
		<< " the columns in the dumpave." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  else if ( dumpmode == 1 )
    {
      for ( int i=0; i<ndim; ++i )
	{
	  idxs[i] = i;
	}
    }
  else
    {
      maxcol = 0;
      for ( int i=0; i<ndim; ++i )
	{
	  maxcol = std::max(maxcol,idxs[i]);
	}
      ++maxcol;
    }
  
  nobs.assign(ndumps,0);
  vals.assign(ndim*ndumps,0);
  errs.assign(ndim*ndumps,0);
  
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
  for ( int i=0; i<ndumps; ++i )
    {
      std::ifstream cin;
      cin.open( idumps[i].c_str() );
      std::vector<double> cols(maxcol,0);
      std::vector<double> data;
      std::string line;
      while ( std::getline(cin,line) )
	{
	  if ( line.size() > 0 )
	    {
	      std::stringstream sline(line);
	      for ( int c=0; c<maxcol; ++c )
		{
		  sline >> cols[c];
		}
	      for ( int k=0; k<ndim; ++k )
		{
		  // + 1 because column zero is the time value
		  data.push_back( cols[ idxs[k] + 1 ] );
		}
	    };
	}
      int nsamples = data.size() / ndim;
      nobs[i] = nsamples;
      std::vector<double> avgs(ndim,0);
      for ( int j=0; j<nsamples; ++j )
	{
	  for ( int k=0; k<ndim; ++k )
	    {
	      avgs[k] += data[k+j*ndim];
	    }
	}
      for ( int k=0; k<ndim; ++k )
	{
	  avgs[k] /= nsamples;
	}
      std::vector<double> vars(ndim,0);
      for ( int j=0; j<nsamples; ++j )
	{
	  for ( int k=0; k<ndim; ++k )
	    {
	      double dx = data[k+j*ndim] - avgs[k];
	      vars[k] += dx*dx;
	    }
	}
      for ( int k=0; k<ndim; ++k )
	{
	  vars[k] /= nsamples;
	}

      
      
      for ( int k=0; k<ndim; ++k )
	{
	  vals[k+i*ndim] = avgs[k];
	  errs[k+i*ndim] = std::sqrt(vars[k]);
	}
    }
}



void ReadMetafile( std::string fname,
		   std::size_t const ndim,
		   std::vector<double> & rcs,
		   std::vector<double> & fcs )
{
  rcs.resize(0);
  fcs.resize(0);

  std::ifstream cin;
  cin.open( fname.c_str() );
  std::string line;
  while ( std::getline(cin,line) )
    {
      std::vector<std::string> tokens;
      std::istringstream iss(line);
      std::copy( std::istream_iterator<std::string>(iss),
                 std::istream_iterator<std::string>(),
                 std::back_inserter(tokens) );
      if ( tokens.size() == 0 )
        {
          continue;
        }
      else
        {
          int ndim2 = (tokens.size()-3);
          if ( ndim2 % 2 )
            {
              //std::cerr << "Incorrect number of columns in file "
	      // << fname << " on line:\n" << line << "\n";
              //std::exit(EXIT_FAILURE);
	      throw 20;
            };
          
	  std::size_t myndim = ndim2/2;
	  if ( ndim != myndim )
	    {
	      throw 21;
	      //std::cerr << "Error: ndim mismatch detected in "
	      //<< fname << " : " << myndim
	      //<< " but expected " << ndim << std::endl;
	      //std::exit(EXIT_FAILURE);
	    }

	  if ( IsFloat(tokens[2]) )
	    {
	      throw 21;
	    }
	  
          std::vector<double> biasloc(ndim,0.);
          std::vector<double> biasfcs(ndim,0.);
          for ( std::size_t i=0; i<ndim; ++i )
            {
              biasloc[i] = std::atof(tokens[3+2*i].c_str());
              biasfcs[i] = std::atof(tokens[4+2*i].c_str());
            };

	  for ( std::size_t i=0; i<ndim; ++i )
	    {
	      rcs.push_back( biasloc[i] );
	      fcs.push_back( biasfcs[i] );
	    }
	};
    }
}



void ReadTxtPath( std::string fname,
		  std::size_t const ndim,
		  std::vector<double> & rcs )
{
  rcs.resize(0);
  std::ifstream cin;
  cin.open( fname.c_str() );
  std::string line;
  while ( std::getline(cin,line) )
    {
      std::vector<std::string> tokens;
      std::istringstream iss(line);
      std::copy( std::istream_iterator<std::string>(iss),
                 std::istream_iterator<std::string>(),
                 std::back_inserter(tokens) );
      if ( tokens.size() == 0 )
        {
          continue;
        }
      else
        {
	  
	  
	  // std::printf("line %s\n",line.c_str());
	  // std::printf("iss %s\n",iss.str().c_str());

	  // for ( std::size_t it=0; it<tokens.size(); ++it )
	  //   {
	  //     std::printf("token: %s\n",tokens[it].c_str());
	  //   }
	  
	  int ncol = (tokens.size()-2);
	  if ( ncol >= (int)ndim )
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  rcs.push_back( std::atof(tokens[2+dim].c_str()) );
		  //fcs.push_back( 100. );
		  //std::printf("read %12.5f\n",rcs.back());
		}
	    }
	  else
	    {
	      throw 22;
	    }
	}
    }
}




std::vector<double> ReadReservations
( std::string curdir )
{
  std::vector<double> vs;

  std::string fname = curdir + "/reservations.txt";
  if ( FileExists(fname) )
    {
      std::cout << "\nReading " << fname << std::endl;
      
      std::ifstream cin;
      cin.open( fname.c_str() );

      std::string line;
      while ( std::getline(cin,line) )
	{
	  std::stringstream sline(line);
	  double v;
	  while ( sline >> v )
	    {
	      vs.push_back(v);
	    }
	};
    }
  return vs;
}

void WriteReservations
( std::string curdir, std::size_t const ndim, std::vector<double> vs )
{
  std::string fname = curdir + "/reservations.txt";
  std::cout << "\nWriting " << fname << std::endl;
  std::ofstream cout;
  cout.open( fname.c_str() );
  std::size_t npts = vs.size()/ndim;
  for ( std::size_t i=0; i<npts; ++i )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  cout << std::scientific << std::setw(19)
	       << std::setprecision(10) << vs[dim+i*ndim];
	}
      cout << "\n";
    }
   cout.close();
}



std::vector<double> ReadPrevUnoccSims
( std::string prevdir,
  std::shared_ptr<ndfes::FES> fes )
{
  std::vector<double> prcs;

  std::string fname = prevdir + "/unoccsims.txt";
  if ( FileExists(fname) )
    {
      std::size_t ndim = fes->mDimInfo.GetNumDims();
      std::ifstream cin;
      cin.open( fname.c_str() );
      
      std::string line;
      while ( std::getline( cin, line ) )
	{
	  std::vector<std::string> tokens;
	  std::istringstream iss(line);
	  std::copy( std::istream_iterator<std::string>(iss),
		     std::istream_iterator<std::string>(),
		     std::back_inserter(tokens) );
	  if ( tokens.size() == 0 )
	    {
	      continue;
	    }
	  else if ( tokens.size() >= ndim )
	    {
	      std::vector<double> rcs;
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  rcs.push_back( std::atof( tokens[dim].c_str() ) );
		}

	      
	      bool isocc = fes->PointIsInRangeAndNonzeroOcc( rcs.data() );
	      
	      if ( ! isocc )
		{
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      prcs.push_back( rcs[dim] );
		    }
		}
	    }
	}	
    }
  
  return prcs;
}


std::vector<std::size_t> CheckUnoccSims
( std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions popts,
  std::vector<double> rcs,
  std::vector<double> & fcs,
  std::vector<double> & prcs )
{
  std::vector<std::size_t> imodified;
  
  std::size_t ndim = fes->mDimInfo.GetNumDims();
  std::size_t nsim = rcs.size() / ndim;
  std::size_t nprev = prcs.size() / ndim;
  for ( std::size_t i=0; i<nsim; ++i )
    {
      bool isocc = fes->PointIsInRangeAndNonzeroOcc( rcs.data() + i*ndim );

      if ( ! isocc )
	{
	  //std::cout << "PT " << i << " is unocc" << std::endl;
	  bool found = false;
	  for ( std::size_t j=0; j<nprev; ++j )
	    {
	      double d = 0;
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  double dx = rcs[dim+i*ndim] - prcs[dim+j*ndim];
		  d += dx*dx;
		}
	      d = std::sqrt(d);
	      if ( d < 1.e-5 )
		{
		  found = true;
		  break;
		}
	    }
	  if ( found )
	    {
	      imodified.push_back(i);
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  fcs[dim+i*ndim] = popts.maxfc[dim];
		}
	    }
	  else
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  prcs.push_back( rcs[dim+i*ndim] );
		}
	    }
	}
    }
  return imodified;
}


void WritePrevUnoccSims
( std::string curdir,
  std::size_t ndim,
  std::vector<double> prcs )
{
  std::size_t nsim = prcs.size() / ndim;
  std::string fname = curdir + "/unoccsims.txt";
  std::ofstream cout;
  cout.open( fname.c_str() );
  for ( std::size_t i=0; i<nsim; ++i )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  cout << std::fixed << std::setw(17) << std::setprecision(10)
	       << prcs[dim+i*ndim];
	}
      cout << "\n";
    }
  cout.close();
}


void PrintPrevUnoccSims( std::vector<std::size_t> modfcidxs, std::ostream & fh )
{
  if ( modfcidxs.size() > 0 )
    {
      fh << "\n";
      fh << "There are " << modfcidxs.size() << " proposed simulations in "
	 << "unoccupied bins that were previously simulated.\n";
      fh << "Setting the force constants to the maximum allowable values.\n";
      for ( std::size_t imod=0; imod<modfcidxs.size(); ++imod )
	{
	  fh << "isim " << std::setw(3) << modfcidxs[imod] + 1
	     << " FC SET TO MAX\n";
	}
    }
}



namespace ndfes
{

  class SimData
  {
  public:

    SimData( Options const & my_opts );

    //ndfes::PathIter MakeInitialGuess() const;

    std::shared_ptr< ndfes::FES > MakeInterp() const;
    
    bool mIsPrepareSims;
    
    Options mOpts;
    ndfes::PathOptions mPopts;

    double mBeta;
    std::size_t mNumDim;
    std::size_t mNumOutSim;
    std::size_t mNumPathPts;
    std::vector<int> mTidxs;
    std::vector<bool> mTidxIsAngle;
    std::vector<double> mDefFC;
    double mNumSamples;


    ndfes::PathIter mInitGuess;
    ndfes::PathSims mPrevSimAvgs;
    
    std::shared_ptr< ndfes::FES > mFES;
    //std::shared_ptr< ndfes::FES > mPenFES;
    ndfes::PathOpt mSimPaths;
  };
  
  
  class PathOptimize
  {
  public:
    
    PathOptimize( std::shared_ptr<ndfes::SimData> my_simdata,
		  std::shared_ptr<ndfes::FES> my_fes );

    void PredictPathSims
    ( std::size_t const nsim,
      std::vector<double> & rcs,
      std::vector<double> & fcs,
      std::vector<double> & ts );

    void PredictUniformPathSims
    ( std::size_t const nsim,
      std::vector<double> & rcs,
      std::vector<double> & fcs,
      std::vector<double> & ts );

    ndfes::PathOpt
    GetMicroIters( std::size_t const nmax ) const;
        
    std::shared_ptr<ndfes::PCurve> GetOutputSpline() const;

    std::shared_ptr<ndfes::SimData> mSimData;
    std::shared_ptr<ndfes::FES> mFES;
    ndfes::PathOpt mPaths;
    ndfes::PathOpt mSimPaths;
    
  };

}

ndfes::SimData::SimData
( Options const & my_opts )
  : mIsPrepareSims( my_opts.ipath.size() == 0 ),
    mOpts(my_opts),
    mBeta( ndfes::GetBeta( my_opts.temp ) ),
    mNumDim(0),
    mNumPathPts(0),
    mNumSamples(100.)
{

  if ( not mIsPrepareSims )
    {
      //////////////////////////////////////////////
      // CLI OPTIMIZATION
      //////////////////////////////////////////////
      
      // CHECK DIRS AND FILES

      if ( ! FileExists(mOpts.ipath) )
	{
	  std::cerr << "ipath " << mOpts.ipath << " does not exist"
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      
      if ( ! FileExists(mOpts.chk) )
	{
	  std::cerr << "chk " << mOpts.chk << " does not exist"
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      

      
      std::string chk = mOpts.chk;
      std::cout << "Reading model " << mOpts.model
		<< " from " << chk << "\n";
      //fes = ndfes::ReadFES( chk, mOpts.model );
      mFES.reset( new ndfes::FES( chk, mOpts.model ) );
      
      mNumDim = mFES->mDimInfo.GetNumDims();


      mTidxIsAngle.assign(mNumDim,false);
      for ( std::size_t dim=0; dim<mNumDim; ++dim )
	{
	  double w = mFES->mDimInfo.GetXmax()[dim] - mFES->mDimInfo.GetXmin()[dim];
	  if ( mFES->mDimInfo.IsPeriodic(dim) or w > 15. )
	    {
	      mTidxIsAngle[dim] = true;
	    }
	}
  
      std::vector<double> minfc( mNumDim, 0 );
      std::vector<double> maxfc( mNumDim, 0 );
      mDefFC.assign( mNumDim, 100. );
      for ( std::size_t dim=0; dim<mNumDim; ++dim )
	{
	  if ( mTidxIsAngle[dim] )
	    {
	      minfc[dim] = mOpts.minfc_ang;
	      maxfc[dim] = mOpts.maxfc_ang;
	      if ( mDefFC[dim] < minfc[dim] or mDefFC[dim] > maxfc[dim] )
		{
		  mDefFC[dim] = 0.5*(minfc[dim]+maxfc[dim]);
		}
	    }
	  else
	    {
	      minfc[dim] = mOpts.minfc;
	      maxfc[dim] = mOpts.maxfc;
	      mDefFC[dim] = 1.;
	      if ( mDefFC[dim] < minfc[dim] or mDefFC[dim] > maxfc[dim] )
		{
		  mDefFC[dim] = 0.5*(minfc[dim]+maxfc[dim]);
		}
	    }
	}


      
      try
	{
	  std::vector<double> rcs;
	  std::vector<double> fcs;
	  ReadMetafile( mOpts.ipath, mNumDim, rcs, fcs );
	  std::size_t nsim = rcs.size()/mNumDim;
	  //mDefFC = CalcAvg( mNumDim, nsim, fcs );
	  mInitGuess.mPath = ndfes::PathSpl( mOpts.stype, false, mNumDim, nsim, rcs.data() );
	  mInitGuess.mSims = ndfes::PathSims( mNumDim, nsim, rcs.data(), fcs.data() );
	}
      catch( int ex1 )
	{
	  try
	    {
	      std::vector<double> rcs;
	      ReadTxtPath( mOpts.ipath, mNumDim, rcs );
	      std::size_t nsim = rcs.size()/mNumDim;
	      //mDefFC = CalcAvg( mNumDim, nsim, fcs );
	      std::vector<double> fcs( ExpandAvg( mNumDim, nsim, mDefFC ) );
	      mInitGuess.mPath = ndfes::PathSpl( mOpts.stype, false, mNumDim, nsim, rcs.data() );
	      mInitGuess.mSims = ndfes::PathSims( mNumDim, nsim, rcs.data(), fcs.data() );
	    }
	  catch( int ex2 )
	    {
	      ndfes::PathOpt pathopt( mOpts.ipath );
	      mInitGuess = pathopt.mIters.back();
	      if ( mInitGuess.mSims.mNumSim == 0 )
		{
		  std::size_t nsim = mInitGuess.mPath.mNumPts;
		  std::vector<double> rcs( mInitGuess.mPath.GetUniformPts(nsim) );
		  std::vector<double> fcs( ExpandAvg( mNumDim, nsim, mDefFC ) );
		  mInitGuess.mSims = ndfes::PathSims( mNumDim, nsim, rcs.data(), fcs.data() );
		}
	    }
	}
      mPrevSimAvgs = mInitGuess.mSims;
      if ( mOpts.npathpts > 1 )
	{
	  mNumPathPts = mOpts.npathpts;
	}
      else
	{
	  mNumPathPts = mInitGuess.mPath.mNumPts;
	}
  

      mPopts.npathpts = mOpts.npathpts;
      mPopts.beta = mBeta;
      mPopts.stype = mOpts.stype;
      mPopts.fix0 = mOpts.fix0;
      mPopts.fix1 = mOpts.fix1;
      mPopts.smooth = mOpts.smooth;
      mPopts.smoothfc = mOpts.smoothfc;
      mPopts.tdispfc = mOpts.tdispfc;
      mPopts.minfc = minfc;
      mPopts.maxfc = maxfc;
      mPopts.deffc = mDefFC;
      // if ( mOpts.acc )
      // 	{
      // 	  mOpts.maxit = 0;
      // 	  mOpts.acc = true;
      // 	}
      mPopts.acc = true;
      mPopts.acc_oobk = mOpts.acc_oobk;
      mPopts.acc_maxit = mOpts.acc_maxit;
  
    }
  else
    {
      //////////////////////////////////////////////
      // FTSM SIMULATION
      //////////////////////////////////////////////

      // CHECK DIRS AND FILES
      
      std::string curdir = ndfes::GetDirName(mOpts.curit,mOpts.pad);
      std::string nextdir = ndfes::GetDirName(mOpts.curit+1,mOpts.pad);
      std::string prevdir = "";
      if ( mOpts.curit > 0 )
	{
	  prevdir = ndfes::GetDirName(mOpts.curit-1,mOpts.pad);
	}
      if ( mOpts.odir.size() > 0 )
	{
	  nextdir = mOpts.odir;
	}
      
      if ( ! DirExists(curdir) )
	{
	  std::cerr << "Current directory does not exist " << curdir
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      
      if ( ! DirExists(nextdir) )
	{
	  std::cerr << "Directory " << nextdir << " does not exist"
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}

      if ( mOpts.prefix.size() > 0 )
	{
	  std::string tdir =  curdir + "/" + mOpts.prefix;
	  if ( ! DirExists(tdir) )
	    {
	      std::cerr << "Directory does not exist " << tdir
			<< std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	  
	  tdir =  nextdir + "/" + mOpts.prefix;
	  if ( ! DirExists(tdir) )
	    {
	      std::cerr << "Directory does not exist " << tdir
			<< std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	}
      
      if ( ! FileExists(mOpts.mdin) )
	{
	  std::cerr << "mdin " << mOpts.mdin << " does not exist"
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      
      if ( ! FileExists(mOpts.disang) )
	{
	  std::cerr << "Template disang " << mOpts.disang << " does not exist"
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}


      // FIND DISANGS AND DUMPAVES
      
      std::vector< std::string > idisangs;
      std::vector< std::string > idumps;
      std::vector< std::string > irsts;
      int dumpmode = 0;
      
      for ( int i=0; i<200; ++i )
	{
	  std::string simg = ndfes::GetImgName(i+1,mOpts.pad);
	  std::stringstream base;
	  base << curdir << "/";
	  if ( mOpts.prefix.size() > 0 )
	    {
	      base << mOpts.prefix << "/";
	    }
	  base << simg;
	  std::string rst = base.str() + ".rst7";
	  std::string disang = base.str() + ".disang";
	  std::string cdump = base.str() + ".dumpave";
	  base.str("");
	  base.clear();
	  base << curdir << "/analysis/dumpaves/" << curdir << "/";
	  if ( mOpts.prefix.size() > 0 )
	    {
	      base << mOpts.prefix << "/";
	    }
	  base << simg << ".dumpave";
	  std::string idump = base.str();
	  
	  if ( FileExists(rst) )
	    {
	      if ( ! FileExists(disang) )
		{
		  std::cerr << "Restart " << rst << " does not have a "
			    << "corresponding disang " << disang << std::endl;
		  std::exit(EXIT_FAILURE);
		}
	      if ( FileExists(idump) )
		{
		  if ( dumpmode == 0 )
		    {
		      dumpmode = 1;
		    }
		  else if ( dumpmode == 2 )
		    {
		      std::cerr << "Not all of the dumpave files are "
				<< "present in " << curdir
				<< "/analysis/dumpaves/";
		      if ( mOpts.prefix.size() > 0 )
			{
			  std::cerr << mOpts.prefix << "/";
			}
		      std::cerr << curdir << std::endl;
		      std::exit(EXIT_FAILURE);
		    };
		  irsts.push_back(rst);
		  idisangs.push_back(disang);
		  idumps.push_back(idump);
		}
	      else if ( FileExists(cdump) )
		{
		  if ( dumpmode == 0 )
		    {
		      dumpmode = 2;
		    }
		  else if ( dumpmode == 1 )
		    {
		      std::cerr << "Some of the dumpaves are located in "
				<< curdir << " and others in "
				<< curdir << "/analysis/dumpaves/";
		      if ( mOpts.prefix.size() > 0 )
			{
			  std::cerr << mOpts.prefix << "/";
			}
		      std::cerr << curdir << std::endl;
		      std::exit(EXIT_FAILURE);
		    }
		  
		  irsts.push_back(rst);
		  idisangs.push_back(disang);
		  idumps.push_back(cdump);
		}
	      else
		{
		  std::cerr << "Restart " << rst << " does not have a "
			    << "corresponding dumpave " << idump << std::endl;
		  std::exit(EXIT_FAILURE);
		}
	      
	      std::cout << "IMG: \n"
			<< "     restart:  " << irsts.back() << "\n"
			<< "     dumpave:  " << idumps.back() << "\n"
			<< "     disang:   " << idisangs.back() << "\n";
	      
	    }
	  else
	    {
	      break;
	    };
	}
      std::cout << "\n";

      if ( mOpts.npathpts > 1 )
	{
	  mNumPathPts = mOpts.npathpts;
	}
      else
	{
	  mNumPathPts = irsts.size();
	}
      
      std::cout << "Reading template disang: " << mOpts.disang << "\n\n";
  
      //ndfes::ReadTemplateDisang( mOpts.disang, tidxs, tidxisangle, deffc );
      {
	amber::Disang tdisang( mOpts.disang );
	mTidxs = tdisang.GetTemplateIdxs();
	mTidxIsAngle = tdisang.GetTemplateIdxIsAngle();
	mDefFC = tdisang.GetTemplateForceConsts();
      }
      mNumDim = mTidxs.size();
      
      std::cout << "Num sims: " << irsts.size() << "\n\n";
      std::cout << "Num dims: " << mNumDim << "\n\n";
      std::cout << "Restraint indexes:";
      for ( std::size_t dim=0; dim<mNumDim; ++dim )
	{
	  std::cout << std::setw(3) << mTidxs[dim]+1;
	}
      std::cout << "\n\n";
      
      std::cout << "Is angle: [";
      for ( std::size_t dim=0; dim<mNumDim; ++dim )
	{
	  if ( mTidxIsAngle[dim] )
	    {
	      std::cout << "True";
	    }
	  else
	    {
	      std::cout << "False";
	    }
	  if ( dim+1 < mNumDim )
	    {
	      std::cout << ", ";
	    }
	  else
	    {
	      std::cout << "]\n\n";
	    }
	}
      
      
      // READ DUMPAVES / DISANGS
      {
	std::vector<double> nobs;
	std::vector<double> avgs;
	std::vector<double> stds;
	ReadPathFromDumpaves( mTidxs, idumps, dumpmode, nobs, avgs, stds );
	std::vector<double> rcs( mNumDim*idisangs.size(), 0 );
	std::vector<double> fcs( mNumDim*idisangs.size(), 0 );
	for ( std::size_t idis=0; idis<idisangs.size(); ++idis )
	  {
	    amber::Disang tdis( idisangs[idis] );
	    std::vector<double> trcs( tdis.GetCenters(mTidxs) );
	    std::vector<double> tfcs( tdis.GetForceConsts(mTidxs) );
	    for ( std::size_t dim=0; dim<mNumDim; ++dim )
	      {
		rcs[dim+idis*mNumDim] = trcs[dim];
		fcs[dim+idis*mNumDim] = tfcs[dim];
	      }
	  };
	mNumSamples = CalcAvg( nobs );
	mPrevSimAvgs = ndfes::PathSims( mNumDim, idisangs.size(), avgs.data(), fcs.data() );
	mInitGuess.mPath = ndfes::PathSpl( mOpts.stype, false, mNumDim, idisangs.size(), rcs.data() );
	mInitGuess.mSims = ndfes::PathSims( mNumDim, idisangs.size(), rcs.data(), fcs.data() );
      }
      
      std::vector<double> minfc( mNumDim, 0 );
      std::vector<double> maxfc( mNumDim, 0 );
      for ( std::size_t dim=0; dim<mNumDim; ++dim )
	{
	  if ( mTidxIsAngle[dim] )
	    {
	      minfc[dim] = mOpts.minfc_ang;
	      maxfc[dim] = mOpts.maxfc_ang;
	    }
	  else
	    {
	      minfc[dim] = mOpts.minfc;
	      maxfc[dim] = mOpts.maxfc;
	    }
	};
      
      
      mPopts.npathpts = mOpts.npathpts;
      mPopts.beta = mBeta;
      mPopts.stype = mOpts.stype;
      mPopts.fix0 = mOpts.fix0;
      mPopts.fix1 = mOpts.fix1;
      mPopts.smooth = mOpts.smooth;
      mPopts.smoothfc = mOpts.smoothfc;
      mPopts.tdispfc = mOpts.tdispfc;
      mPopts.minfc = minfc;
      mPopts.maxfc = maxfc;
      mPopts.deffc = mDefFC;
      mPopts.acc = true;
      if ( mOpts.msm and mOpts.maxit < 1 )
	{
	  mPopts.acc = false;
	}
      mPopts.acc_oobk = mOpts.acc_oobk;
      mPopts.acc_maxit = mOpts.acc_maxit;
  
      if ( mPopts.acc )
	{
	  std::stringstream chk;
	  chk << curdir << "/analysis/metafile.all.chk";
	  if ( ! FileExists( chk.str() ) )
	    {
	      std::cerr << "ndfes checkpoint file does not exist "
			<< chk.str() << std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	  else
	    {
	      std::cout << "Reading model " << mOpts.model
			<< " from " << chk.str() << "\n";
	    }
	  //fes = ndfes::ReadFES( chk.str(), mOpts.model );
	  mFES.reset( new ndfes::FES( chk.str(), mOpts.model ) );
	}
	
      if ( mOpts.curit > 0 )
	{
	  std::stringstream pathxml;
	      
	  pathxml << prevdir;
	  if ( mOpts.prefix.size() > 0 )
	    {
	      pathxml << "/" << mOpts.prefix;
	    }
	  pathxml << "/path.sims.xml";
	  if ( ! FileExists( pathxml.str() ) )
	    {
	      std::cerr << "File not found: " << pathxml.str() << std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	    
	  std::cout << "Reading previous path.sims from: " << pathxml.str() << "\n\n\n\n";
	    
	  //simpaths = ndfes::ReadPathsFromPkl( pathpkl.str() );
	  mSimPaths = ndfes::PathOpt( pathxml.str() );
	  mInitGuess = mSimPaths.mIters.back();
	}
      else
	{
	  mSimPaths.mIters.push_back(mInitGuess);
	};
    } // IsPrepareSims


  if ( mFES )
    {
      if ( (not mOpts.rbf) and (not mOpts.hist) and
	   (not mOpts.bspl) and (mOpts.arbf<1) and
	   (mOpts.wavg<3) )
	{
	  if ( mFES->mIsMBAR )
	    {
	      mOpts.rbf = true;
	    }
	  else
	    {
	      mOpts.bspl = true;
	    }
	}
    }
  
}
    

// ndfes::PathIter ndfes::SimData::MakeInitialGuess() const
// {
//   ndfes::PathIter initpath;
//   std::size_t NumSim = mNumSim;
//   std::size_t NumDim = mNumDim;
//   std::size_t NumPts = mNumPathPts;
//   if ( not IsPrepareSims )
//     {
//       std::vector<double> new_fcs( ExpandAvg( NumDim, NumSim, deffc ) );
	
//       initpath.reset( new ndfes::Path( NumDim, NumSim,
// 				       rcs, new_fcs, nobs,
// 				       avgs, stds,
// 				       mPopts ) );
      
//       if ( NumPts != NumSim )
// 	{
// 	  std::vector<double> myrcs( NumDim * NumPts, 0 );
// 	  std::vector<double> myfcs( NumDim * NumPts, 0 );
// 	  double mynobsval = nobs[0];
// 	  std::vector<double> mynobs( NumPts, mynobsval );
// 	  std::vector<double> mystds( NumDim * NumPts, 0 );
// 	  double nm1 = NumPts-1;
// 	  for ( std::size_t i=0; i<NumPts; ++i )
// 	    {
// 	      initpath->pspl->GetValue( i/nm1, myrcs.data()+i*NumDim );
// 	      for ( std::size_t dim=0; dim<NumDim; ++dim )
// 		{
// 		  myfcs[dim+i*NumDim] = deffc[dim];
// 		  mystds[dim+i*NumDim] = 1. / (2.*(popts.beta)*deffc[dim]);
// 		}
// 	    }
// 	  initpath.reset( new ndfes::Path( NumDim, NumPts,
// 					   myrcs, myfcs,
// 					   mynobs, myrcs, mystds,
// 					   popts ) );
// 	}
//     }
//   else
//     {
      
//       if ( fes )
// 	{
// 	  if ( mOpts.curit > 0 )
// 	    {
// 	      {
// 		std::cout << "Setting initial guess to " << GetDirName(opts.curit-1) << "\n\n";
// 		initpath.reset( new ndfes::Path( *(simpaths.back()) ) );
// 	      }
	      
// 	      {
// 		initpath = initpath->Next( popts );
// 		// This will need ->UpdateDensity once we have a fes object
// 	      }
// 	    } // curit > 0
// 	  else // curit == 0
// 	    {
// 	      if ( NumSim != NumPts )
// 		{
// 		  ndfes::PCurve pspl( NumDim, NumSim, rcs.data(), true, 100, 1000 );
		  
// 		  std::vector<double> ts( NumPts, 0 );
// 		  std::vector<double> myrcs( NumDim*NumPts, 0 );
// 		  double nm1=NumPts-1;
// 		  for ( std::size_t i=0; i<NumPts; ++i )
// 		    {
// 		      ts[i] = i / nm1;
// 		      pspl.GetValue( ts[i], myrcs.data() + i*NumDim );
// 		    }
		  
// 		  std::vector<double> myfcs( ExpandAvg( NumDim, NumPts, deffc ) );
// 		  std::vector<double> mystds( ExpandAvg( NumDim, NumPts, CalcAvg( NumDim, NumSim, stds ) ) );
// 		  std::vector<double> mynobs(NumPts,mNumSamples);
		  
// 		  initpath.reset
// 		    ( new ndfes::Path( NumDim, NumPts,
// 				       myrcs, myfcs,
// 				       mynobs, myrcs, mystds,
// 				       popts ) );
		  
// 		} // nsim != NumRealSim
// 	      else 
// 		{ // nsim == NumRealSim
		  
// 		  std::vector<double> myfcs( ExpandAvg( NumDim, NumPts, deffc ) );
		  
// 		  initpath.reset
// 		    ( new ndfes::Path( NumDim, NumPts,
// 				       rcs, myfcs,
// 				       nobs, avgs, stds,
// 				       popts ) );
	      
// 		} // nsim == NumRealSim
// 	    } // curit == 0
// 	} // if fes
//       else
// 	{ // not fes
// 	  if ( NumSim != NumPts )
// 	    {
// 	      std::cerr << "ERROR: Traditional FTSM requires --npathpts=NSIM or --nspathpts=0" << std::endl;
// 	      std::exit(EXIT_FAILURE);
// 	    }
// 	  else
// 	    {
// 	      std::vector<double> new_fcs( ExpandAvg( NumDim, NumPts, deffc ) );
// 	      initpath.reset( new ndfes::Path( NumDim, NumPts,
// 					       rcs, new_fcs,
// 					       nobs, avgs, stds,
// 					       popts ) );
// 	    } // nsim == NumRealSim
// 	} // not fes
      
//     } // not IsPrepareSims


  
//   return initpath;
// }


std::shared_ptr< ndfes::FES > ndfes::SimData::MakeInterp() const
{
  std::shared_ptr< ndfes::FES > ofes;

  if ( mFES )
    {
      ofes.reset( new ndfes::FES( *mFES ) );

      int minsize = mOpts.minsize;
      
      //ofes->AddOccPen( occpen );

      if ( mOpts.buffer and (mOpts.wavg<3) )
	{
	  ofes = ofes->AddPenaltyBuffer(minsize);
	  minsize=0;
	}
      else if ( mOpts.wavg > 2 )
	{
	  int nbuf = (mOpts.wavg + (mOpts.wavg%2))/2;
	  ofes = ofes->AddPenaltyBuffer(minsize,nbuf);
	}
      else if ( minsize > 0 and mOpts.arbf > 0 )
	{
	  ofes = ofes->DeleteBinsBySize(minsize);
	}

      if ( mOpts.arbf > 0 )
	{
	  std::cout << "Creating ARBF interpolator: " << mOpts.arbf << std::endl;
	  ofes->SetupARBF( mOpts.arbf, true );
	}
      else if ( mOpts.rbf )
	{
	  std::cout << "Creating RBF interpolator" << std::endl;
	  ofes->SetupRBF( mOpts.shape, minsize, mOpts.rbfmaxerr, true );
	}
      else if ( mOpts.hist )
	{
	  std::cout << "Creating Histogram interpolator" << std::endl;
	  ofes->SetupHist();
	}
      else if ( mOpts.bspl )
	{
	  std::cout << "Creating B-Spline interpolator" << std::endl;
	  ofes->SetupBspl();
	}
      else if ( mOpts.wavg > 2 )
	{
	  std::cout << "Creating WAVG interpolator" << std::endl;
	  ofes->SetupWAvg( mOpts.wavg, mOpts.wavg_niter );
	}
      std::cout << "Finished creating interpolator" << std::endl;

      ofes->SetupNoQuad();

    }
  
  return ofes;
}




ndfes::PathOptimize::PathOptimize
( std::shared_ptr<ndfes::SimData> my_simdata,
  std::shared_ptr<ndfes::FES> my_fes )
  : mSimData(my_simdata),
    mFES(my_fes),
    mSimPaths( my_simdata->mSimPaths )
{
  if ( (not mFES) or mSimData->mOpts.maxit < 1 )
    {
      // standard FTSM
      std::cout << "No path optimization\n";
      mPaths.mIters.push_back( mSimData->mInitGuess );
      if ( mSimData->mIsPrepareSims )
	{
	  ndfes::PathSpl pspl( mSimData->mPopts.stype,
			       mSimData->mPopts.smooth,
			       mSimData->mPrevSimAvgs.mNumDim,
			       mSimData->mPrevSimAvgs.mNumSim,
			       mSimData->mPrevSimAvgs.mRCs.data() );
	  mPaths.mIters.push_back( ndfes::PathIter( pspl ) );
	}	  
    }
  else if ( mFES )
    {
      Options opts(mSimData->mOpts);

      mPaths.mIters.push_back( mSimData->mInitGuess );
      for ( int it=0; it<opts.maxit; ++it )
	{
	  bool verbose = false;
	  if ( it == 0 )
	    {
	      verbose = true;
	    }
	  else if ( (it+1) % mSimData->mOpts.printfreq == 0 )
	    {
	      verbose = true;
	    }
	  if ( verbose )
	    {
	      std::cout << "Synthetic iteration " << it+1 << std::endl;
	    }
	  ndfes::PathIter npath( mPaths.mIters.back().Next( mFES, mSimData->mPopts ) );
	  
	  double pmin = 1;
	  bool samemeans = false;
	  ndfes::FTSMCheckSameMeans( npath, mPaths.mIters.back(),
				     mSimData->mPopts.deffc,
				     mSimData->mPopts.beta,
				     mSimData->mOpts.ptol,
				     samemeans, pmin, verbose );

	  if ( samemeans and not verbose )
	    {
	      std::cout << "Synthetic iteration " << it+1 << std::endl;
	      pmin = 1;
	      samemeans = false;
	      ndfes::FTSMCheckSameMeans( npath, mPaths.mIters.back(),
					 mSimData->mPopts.deffc,
					 mSimData->mPopts.beta,
					 mSimData->mOpts.ptol,
					 samemeans, pmin, true );
	    }
	  
	  mPaths.mIters.push_back( npath );
	  
	  if ( samemeans )
	    {
	      std::cout << "Converged because means are same\n";
	      break;
	    }
	  else
	    {
	      bool samepath = mPaths.CheckRepeatedPath(mFES);
	      if ( samepath )
		{
		  std::cout << "Terminated because a path was repeated\n";
		  break;
		}
	    }
	  
	  
	  if ( it == opts.maxit - 1 )
	    {
	      std::cout << "Terminated because maxit reached\n"
			<< "Convergence would be met if "
			<< "ptol <= " << std::fixed << std::setprecision(3)
			<< pmin << "\n";
	    }
	  
	} // main it loop
      
    } // has fes
  
  mSimPaths.mIters.push_back( mPaths.mIters.back() );
  
}


std::shared_ptr<ndfes::PCurve> ndfes::PathOptimize::GetOutputSpline() const
{
  return mPaths.mIters.back().mPath.mSpl;
}




void ndfes::PathOptimize::PredictPathSims
( std::size_t const nsim,
  std::vector<double> & rcs,
  std::vector<double> & fcs,
  std::vector<double> & ts )
{
  if ( (not mSimData->mOpts.predfc) or (not mFES) )
    {
      PredictUniformPathSims( nsim, rcs, fcs, ts );
    }
  else
    {
      ndfes::PathOptions topts( mSimData->mPopts );
      topts.ClearBounds();
      mPaths.mIters.back().PredictUniformCentroids( nsim, mFES, topts, ts );
      mSimPaths.mIters.back() = mPaths.mIters.back();
      rcs = mPaths.mIters.back().mSims.mRCs;
      fcs = mPaths.mIters.back().mSims.mFCs;
    }
}

void ndfes::PathOptimize::PredictUniformPathSims
( std::size_t const nsim,
  std::vector<double> & rcs,
  std::vector<double> & fcs,
  std::vector<double> & ts )
{
  ndfes::PathOptions topts( mSimData->mPopts );
  topts.ClearBounds();
  mPaths.mIters.back().PredictUniformSims( nsim, topts );
  mSimPaths.mIters.back() = mPaths.mIters.back();
  rcs = mPaths.mIters.back().mSims.mRCs;
  fcs = mPaths.mIters.back().mSims.mFCs;
  ts.assign( nsim, 0 );
  double const dt = 1./(nsim-1.);
  for ( std::size_t i=0; i<nsim; ++i )
    {
      ts[i] = i * dt;
    };
}



ndfes::PathOpt
ndfes::PathOptimize::GetMicroIters( std::size_t const nmax ) const
{
  ndfes::PathOpt micro;
  std::size_t npath = mPaths.mIters.size();
  double didx = (npath-1.) / (nmax-1.);
  int iold = -1;
  for ( std::size_t j=0; j<nmax; ++j )
    {
      int i = j * didx + 0.5;
      if ( i != iold )
	{
	  micro.mIters.push_back( mPaths.mIters[i] );
	  iold = i;
	}
    }
  return micro;
}

////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////


std::vector<double> CptDisplacedMovement
( std::size_t const ndim,
  std::size_t const nlev,
  double const * twidths,
  ndfes::PathSpl const & initpath,
  double const * pt )
{
  std::vector<double> dvec( initpath.mSpl->GetPointClosestTo( pt ) );
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      dvec[dim] = pt[dim] - dvec[dim];
    } 
  double norm = 0;
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      norm += dvec[dim]*dvec[dim];
    };
  norm = std::sqrt(norm);
  if ( norm < 1.e-8 )
    {
      std::fill(dvec.begin(),dvec.end(),0.);
    }
  else
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  dvec[dim] /= norm;
	}
	  
      double dmin = 1.e+30;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  //std::printf("dvec %20.10e\n",dvec[dim]);
	  if ( std::abs(dvec[dim]) > 1.e-8 )
	    {
	      //std::printf("%12.5f %12.5f %12.5f\n",dvec[dim],twidths[dim], twidths[dim]/std::abs(dvec[dim]));
	      dmin = std::min(dmin, nlev*twidths[dim]/std::abs(dvec[dim]));
	    }
	}
      if ( dmin > 10000 )
	{
	  dmin = 1;
	}
      //std::printf("scale %lu\n",i);
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  //std::printf(" %12.5f %12.5f %12.5f %12.5f\n",dvec[dim],dmin,pts1[dim+i*ndim],dmin*dvec[dim]+pts1[dim+i*ndim]);
	  dvec[dim] *= dmin;
	}
    }

  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      dvec[dim] += pt[dim];
    }
  
  return dvec;
}


void ExtendSimDisplacements
( ndfes::PathSpl const & initpath,
  ndfes::PathSpl const & optpath,
  std::shared_ptr<ndfes::FES> fes,
  std::size_t const nsim,
  double * rcs )
{
  std::size_t ndim = fes->mDimInfo.GetNumDims();
  double const * twidths = fes->mDimInfo.GetTargetWidths();

  //std::vector<double> pts0(initpath.GetUniformPts(nsim)); 
  std::vector<double> pts1(optpath.GetUniformPts(nsim));
  for ( std::size_t i=0; i<nsim; ++i )
    {
      std::vector<double> dvec( initpath.mSpl->GetPointClosestTo( pts1.data() + i *ndim ) );
      for ( std::size_t dim=0; dim<ndim; ++dim )
      {
         dvec[dim] = pts1[dim+i*ndim] - dvec[dim];
      } 
      double norm = 0;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  norm += dvec[dim]*dvec[dim];
	};
      norm = std::sqrt(norm);
      //std::printf("norm %20.10e\n",norm); 
      if ( norm < 1.e-8 )
	{
	  std::fill(dvec.begin(),dvec.end(),0.);
	}
      else
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      dvec[dim] /= norm;
	    }
	  
	  double dmin = 1.e+30;
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      //std::printf("dvec %20.10e\n",dvec[dim]);
	      if ( std::abs(dvec[dim]) > 1.e-8 )
		{
		  //std::printf("%12.5f %12.5f %12.5f\n",dvec[dim],twidths[dim], twidths[dim]/std::abs(dvec[dim]));
		  dmin = std::min(dmin, twidths[dim]/std::abs(dvec[dim]));
		}
	    }
	  if ( dmin > 10000 )
	    {
	      dmin = 1;
	    }
	  //std::printf("scale %lu\n",i);
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      //std::printf(" %12.5f %12.5f %12.5f %12.5f\n",dvec[dim],dmin,pts1[dim+i*ndim],dmin*dvec[dim]+pts1[dim+i*ndim]);
	      dvec[dim] *= dmin;
	    }
	}
      
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  rcs[dim+i*ndim] = dvec[dim] + pts1[dim+i*ndim];
	}
      
    }
}




void ExtendSimDisplacementsIfFewerSamples
( ndfes::PathSpl const & initpath,
  ndfes::PathSpl const & optpath,
  std::shared_ptr<ndfes::FES> fes,
  std::size_t const nsim,
  double * rcs )
{
  std::size_t ndim = fes->mDimInfo.GetNumDims();

  std::vector<double> pts(optpath.GetUniformPts(nsim));
  
  ExtendSimDisplacements(initpath,optpath,fes,nsim,rcs);

  for ( std::size_t i=0; i<nsim; ++i )
    {
      std::vector<double> cext(ndim,0);
      fes->mDimInfo.CptBinCenterFromPt( rcs+i*ndim, cext.data() );
      
      std::vector<double> corg(ndim,0);
      fes->mDimInfo.CptBinCenterFromPt( pts.data()+i*ndim, corg.data() );

      bool extIsReserved = false;
      for ( std::size_t j=0; j<i; ++j )
	{
	  double d = 0;
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      double e = cext[dim] - rcs[dim+j*ndim];
	      d += e*e;
	    }
	  d = std::sqrt(d);
	  //std::printf("%lu %lu ext %12.5f",i,j,d);
	  if ( d < 0.001 )
	    {
	      extIsReserved = true;
	    }
	}

      bool orgIsReserved = false;
      for ( std::size_t j=0; j<i; ++j )
	{
	  double d = 0;
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      double e = corg[dim] - rcs[dim+j*ndim];
	      d += e*e;
	    }
	  d = std::sqrt(d);
	  //std::printf("%lu %lu org %12.5f",i,j,d);
	  if ( d < 0.001 )
	    {
	      orgIsReserved = true;
	    }
	}

      //std::cout << extIsReserved << " " << orgIsReserved << "\n";
      if ( extIsReserved and orgIsReserved )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      rcs[dim+i*ndim] = pts[dim+i*ndim];
	    }
	}
      else if ( orgIsReserved )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      rcs[dim+i*ndim] = cext[dim];
	    }
	}
      else if ( extIsReserved )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      rcs[dim+i*ndim] = corg[dim];
	    }
	}
      else
	{
	  std::size_t norg = 0;
	  std::size_t next = 0;


	  if ( fes->PointIsInRange( corg.data() ) )
	    {
	      std::vector<std::size_t> bidx( fes->mDimInfo.CptBinIdxs( corg.data() ) );
	      std::size_t gidx = fes->mDimInfo.CptGlbBinIdx( bidx.data() );
	      std::unordered_map<std::size_t,std::size_t>::const_iterator
		p = fes->mGlbIdxMap.find(gidx);
	      if ( p != fes->mGlbIdxMap.end() )
		{
		  std::size_t ibin = p->second;
		  norg = fes->mBinSizes[ibin];
		}
	    }
	  
	  
	  if ( fes->PointIsInRange( cext.data() ) )
	    {
	      std::vector<std::size_t> bidx( fes->mDimInfo.CptBinIdxs( cext.data() ) );
	      std::size_t gidx = fes->mDimInfo.CptGlbBinIdx( bidx.data() );
	      std::unordered_map<std::size_t,std::size_t>::const_iterator
		p = fes->mGlbIdxMap.find(gidx);
	      if ( p != fes->mGlbIdxMap.end() )
		{
		  std::size_t ibin = p->second;
		  next = fes->mBinSizes[ibin];
		}
	    }
	  
	  if ( norg == 0 or norg < next )
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  rcs[dim+i*ndim] = corg[dim];
		}
	    }
	  else
	    {
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  rcs[dim+i*ndim] = cext[dim];
		}
	    }
	  
	}
    }
  
}





void SafelyBinPoints
( ndfes::PathSpl const & optpath,
  std::shared_ptr<ndfes::FES> fes,
  std::size_t const nsim,
  double * rcs )
{
  std::size_t ndim = fes->mDimInfo.GetNumDims();
  
  std::vector<double> pts(optpath.GetUniformPts(nsim));
  
  for ( std::size_t i=0; i<nsim; ++i )
    {
      std::vector<double> corg(ndim,0);
      fes->mDimInfo.CptBinCenterFromPt( pts.data()+i*ndim, corg.data() );
      
      bool orgIsReserved = false;
      for ( std::size_t j=0; j<i; ++j )
	{
	  double d = 0;
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      double e = corg[dim] - rcs[dim+j*ndim];
	      d += e*e;
	    }
	  d = std::sqrt(d);
	  //std::printf("%lu %lu org %12.5f\n",i,j,d);
	  if ( d < 0.001 )
	    {
	      orgIsReserved = true;
	    }
	}
      
      // std::cout << orgIsReserved << "\n";
      // for ( std::size_t dim=0; dim<ndim; ++dim )
      // 	{
      // 	  std::printf("%12.5f",pts[dim+i*ndim]);
      // 	}
      // std::printf(" ");
      
      // std::printf("\n");
      
      if ( orgIsReserved )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      rcs[dim+i*ndim] = pts[dim+i*ndim];
	    }
	}
      else
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      rcs[dim+i*ndim] = corg[dim];
	    }
	}
    }
}


void DXSMPlacement
( ndfes::PathSpl const & initpath,
  ndfes::PathSpl const & optpath,
  std::shared_ptr<ndfes::FES> fes,
  int const curit,
  std::size_t const nsim,
  double * rcs )
{
  std::size_t const ndim = fes->mDimInfo.GetNumDims();
  std::fill(rcs,rcs+ndim*nsim,0.);

  std::size_t const sched_levels[4] = { 0, 1, 0, 2 };
  double const sched_disp[3] = { 0., -1/3., 1/3. };
  std::size_t const maxlevel = sched_levels[ curit % 4 ];
  double const defdisp = sched_disp[ curit % 3 ];

  double const * twidths = fes->mDimInfo.GetTargetWidths();
  
  double const dt = 1./(nsim-1.);
  for ( std::size_t isim=0; isim<nsim; ++isim )
    {
      std::vector<double> ts;
      ts.push_back(isim*dt);
      if ( isim > 0 )
	{
	  ts.push_back((isim-1./3.)*dt);
	}
      if ( isim+1 < nsim )
	{
	  ts.push_back((isim+1./3.)*dt);
	  //ts.push_back((isim+1./2.)*dt);
	}

      bool found = false;
      for ( std::size_t lev=0; lev < maxlevel+1; ++lev )
	{
	  for ( std::size_t it=0; it<ts.size(); ++it )
	    {
	      double t = ts[it];

	      // std::printf("Test %3lu %3lu %3lu %12.5f\n",
	      // 		  isim,lev,it,t);
	      
	      std::vector<double> pt(ndim,0);
	      optpath.GetValue(t,pt.data());
	      std::vector<double> dvec(pt);
	      if ( lev > 0 )
		{
		  dvec = CptDisplacedMovement(ndim,lev,twidths,initpath,pt.data());
		}
	      // for ( std::size_t dim=0; dim<ndim; ++dim )
	      // 	{
	      // 	  std::printf("%12.5f",dvec[dim]);
	      // 	}
	      //std::printf("\n");
	      std::size_t nocc = fes->BinOccAtPt(dvec.data());
	      if ( nocc == 0 )
		{
		  //std::printf("found\n");
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      rcs[dim+isim*ndim] = dvec[dim];
		    }
		  found = true;
		  break;
		} 
	    }

	  if ( found )
	    {
	      break;
	    }
	  else
	    {
	      double t = std::max(0.,std::min(1.,(isim+defdisp)*dt));
	      //std::printf("Not found %12.5f %12.5f %3lu\n",defdisp,t,maxlevel);
	      std::vector<double> pt(ndim,0);
	      optpath.GetValue(t,pt.data());
	      std::vector<double> dvec(pt);
	      if ( maxlevel > 0 )
		{
		  dvec = CptDisplacedMovement(ndim,maxlevel,twidths,
					      initpath,pt.data());
		}
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  rcs[dim+isim*ndim] = dvec[dim];
		}
	    }
	}
      
    }
  
}

////////////////////////////////////////////////////////////////////


void OptPathNoUpdate( Options & opts )
{
  //double const beta = ndfes::GetBeta( opts.temp );

  std::shared_ptr< ndfes::SimData > simdata( new ndfes::SimData( opts ) );
  std::shared_ptr< ndfes::FES > fes( simdata->MakeInterp() );

  ndfes::PathOptimize pathopt( simdata, fes );
  
  std::cout << "Returned from optimization" << std::endl;

  if ( pathopt.mPaths.mIters.size() < 1 )
    {
      std::cout << "Optimization yielded no paths" << std::endl;
      return;
    }


  
 

  std::shared_ptr<ndfes::PCurve> ospl = pathopt.GetOutputSpline();
  

  std::size_t npathpts = simdata->mNumPathPts;
  std::size_t ndim = simdata->mNumDim;

  std::vector<double> rcs;
  std::vector<double> fcs;
  std::vector<double> ts;

  double mynobs = simdata->mNumSamples;
  std::vector<double> nobs( npathpts, mynobs );
  
  if ( not opts.predfc )
    {      
      if ( opts.maxit == 0 and opts.npathpts < 2 )
	{
	  ospl = simdata->mInitGuess.mPath.mSpl;
	  rcs = ospl->mX;
	  ts = ospl->mT;
	  npathpts = ts.size();
	  nobs.assign(npathpts,mynobs);
	  fcs = ExpandAvg( ndim, npathpts, simdata->mDefFC );
	}
      else
	{
	  //std::cout << nsim << " " << paths.back()->nsim << std::endl;
	  pathopt.PredictUniformPathSims( npathpts, rcs, fcs, ts );
	}
    }
  else
    {
      pathopt.PredictPathSims( npathpts, rcs, fcs, ts );
    }


  std::cout << "\nFinished creating the new string\n" << std::endl;

  bool isapath = true;
  if ( opts.gbsm )
    {
      isapath = false;

      std::vector<double> reservedpts;
      
      std::vector<ndfes::AreaSummary> simidxs
	( pathopt.mPaths.mIters.back().PredictGridBasedSims
	  ( npathpts, fes, simdata->mPopts, opts.nlayers, reservedpts ) );

      ndfes::SamplingConvPrint
	( simidxs, opts.nlayers,
	  100, std::cout );
      
    }
  

  {
    std::ofstream cout;
    cout.open( opts.opath );
#define FMTE std::scientific << std::setw(20) << std::setprecision(10)
    for ( std::size_t i=0; i<npathpts; ++i )
      {
	cout << std::setw(3) << i+1
	     << FMTE << ts[i];
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    cout << FMTE << rcs[dim+i*ndim];
	  }
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    cout << FMTE << fcs[dim+i*ndim];
	  }
	cout << "\n";
      }
    cout.close();
#undef FMTE
  }

  {
    std::stringstream dname;
    dname << opts.opath << ".";
    if ( fes->mInterpMode == 1 )
      {
	dname << "hist";
      }
    else if ( fes->mInterpMode == 2 )
      {
	dname << "rbf";
      }
    else if ( fes->mInterpMode == 3 )
      {
	dname << "bspl";
      }
    else if ( fes->mInterpMode == 4 )
      {
	dname << "arbf";
      }
    else if ( fes->mInterpMode == 5 )
      {
	dname << "wavg";
      }
    dname << "." << opts.model << ".dat";

    std::ofstream cout;
    cout.open( dname.str().c_str() );

#define FMTF std::fixed << std::setw(14) << std::setprecision(8)
#define FMTH std::fixed << std::setw(9) << std::setprecision(3)
#define FMTI std::fixed << std::setw(13) << std::setprecision(3)
#define FMTE std::scientific << std::setw(15) << std::setprecision(6)
#define FMTB std::scientific << std::setw(23) << std::setprecision(14)

    std::size_t nspl = npathpts;
    std::vector<double> splts( ts );
    std::vector<double> splpts( rcs );
    if ( opts.nsplpts > 1 and isapath )
      {
	nspl = opts.nsplpts;
	splts.resize( nspl );
	splpts.resize( ndim*nspl );
	double nm1 = nspl-1;
	std::shared_ptr< ndfes::PCurve > pspl( pathopt.GetOutputSpline() );
	for ( std::size_t i=0; i<nspl; ++i )
	  {
	    splts[i] = i/nm1;
	    pspl->GetValue(splts[i],splpts.data() + i*ndim);
	  }
      }

    std::vector<double> wsplpts(splpts);
    ndfes::WrapPath( fes->mDimInfo, wsplpts );


    
    for ( std::size_t i=0; i<nspl; ++i )
      {

	/////////////////////////
	// GRADIENT TEST
	/////////////////////////
	// {
	//   double const DEL=5.e-5;
	//   std::vector<double> anag(ndim,0);

	//   for ( std::size_t dim=0; dim<ndim; ++dim )
	//     {
	//       std::printf("%12.5f",splpts[dim+i*ndim]);
	//     }
	//   std::printf("\n");
	  
	//   double anav = 0;
	//   bool ok = fes->GetValue( splpts.data() + i*ndim,
	// 			   anav, anag.data() );
	
	//   if ( ok )
	//     {
	//       std::vector<double> numg(ndim,0);
	//       for ( std::size_t dim=0; dim<ndim; ++dim )
	// 	{
	// 	  splpts[dim+i*ndim] += DEL;
	// 	  double vhi = 0;
	// 	  fes->GetValue(splpts.data() + i*ndim, vhi);
	// 	  splpts[dim+i*ndim] -= 2*DEL;
	// 	  double vlo = 0;
	// 	  fes->GetValue(splpts.data() + i*ndim, vlo);
	// 	  splpts[dim+i*ndim] += DEL;
	// 	  numg[dim] = (vhi-vlo)/(2.*DEL);
	// 	}

	//       std::printf("%4lu",i);
	//       for ( std::size_t dim=0; dim<ndim; ++dim )
	// 	{
	// 	  std::printf(" %12.3e",anag[dim]);
	// 	}
	//       for ( std::size_t dim=0; dim<ndim; ++dim )
	// 	{
	// 	  std::printf(" %12.3e",numg[dim]);
	// 	}
	//       for ( std::size_t dim=0; dim<ndim; ++dim )
	// 	{
	// 	  std::printf(" %12.3e",anag[dim]-numg[dim]);
	// 	}
	//       std::printf("\n");
	      
	      
	//     }
	// }
	/////////////////////////
	/////////////////////////

	cout << std::setw(3) << i+1
	     << FMTB << splts[i];
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    cout << FMTB << splpts[dim+i*ndim];
	  }
	
	double val=0;
	double err=0;
	double re=0;
	std::size_t ns=0;
	//double avgocc = fes->CptAvgOccAtPt( wsplpts.data() + i*ndim, 3 );
	fes->GetValueErrorAndEntropy(wsplpts.data() + i*ndim,val,err,re,ns);
	cout << FMTE << val
	     << FMTI << err
	     << FMTH << re
	     << std::setw(7) << ns
	     << "\n";

      }
#undef FMTE
#undef FMTF
#undef FMTH
#undef FMTI
#undef FMTB
    cout.close();

  }

  if ( not opts.dry_run )
    {
      //ndfes::SavePathPickles( ndim, opts.opath + ".pkl", pathopt.GetMicroIters(10) );
      ndfes::PathOpt micro( pathopt.GetMicroIters(10) );
      micro.WriteXml( opts.opath + ".xml" );
    }

}





//////////////////////////////////////////////////////////////////////////////



int main_new( int argc, char * argv[] )
{

  Options opts( ReadOptions( argc, argv ) );


  //std::printf("opts.occpen %20.10e\n",opts.occpen); 

  if ( opts.ipath.size() > 0 )
    {
      OptPathNoUpdate( opts );
      return 0;
    }


  //double const beta = ndfes::GetBeta( opts.temp );
  std::string curdir = ndfes::GetDirName(opts.curit,opts.pad);
  std::string nextdir = ndfes::GetDirName(opts.curit+1,opts.pad);
  std::string prevdir = "";
  if ( opts.curit > 0 )
    {
      prevdir = ndfes::GetDirName(opts.curit-1,opts.pad);
    }
  if ( opts.odir.size() > 0 )
    {
      nextdir = opts.odir;
    }


  std::shared_ptr<ndfes::SimData> simdata( new ndfes::SimData( opts ) );
  std::shared_ptr<ndfes::FES> fes( simdata->MakeInterp() );

  //std::shared_ptr<ndfes::Path> initpath( simdata->MakeInitialGuess() );

  std::cout << "Creating new string for next iteration in:  " << nextdir << "\n";

  ndfes::PathOptimize pathopt( simdata, fes );

  std::size_t ndim = simdata->mNumDim;
  std::size_t nsim = simdata->mPrevSimAvgs.mNumSim;

  std::shared_ptr<ndfes::PCurve> ospl = pathopt.GetOutputSpline();
    
  //std::vector< std::shared_ptr<ndfes::Path> > paths = pathopt.paths;
  //std::shared_ptr<ndfes::Path> opath( new ndfes::Path( *(paths.back()) ) );
  //std::shared_ptr<ndfes::PCurve> ospl = opath->pspl; //pathopt.GetOutputSpline();
  //bool do_conv = ( opts.conv_layers > 0 and opts.conv_samples > 0 and fes );

  std::vector<double> new_rcs( ndim * nsim, 0 );
  std::vector<double> new_fcs( ndim * nsim, 0 );
  std::vector<double> new_ts( nsim, 0 );

  //std::shared_ptr< ndfes::PathOptimize > penopt;
  if ( not opts.gbsm )
    {
      if ( opts.sasm )
	{

	  pathopt.PredictUniformPathSims( nsim, new_rcs, new_fcs, new_ts );
	  
	  DXSMPlacement( pathopt.mPaths.mIters[0].mPath,
			 pathopt.mPaths.mIters.back().mPath,
			 fes, opts.curit, nsim, new_rcs.data() );
	  
	  pathopt.mPaths.mIters.back().mSims
	    = ndfes::PathSims(ndim,nsim,new_rcs.data(),new_fcs.data());
	  pathopt.mSimPaths.mIters.back() = pathopt.mPaths.mIters.back();
	}
      // else if ( opts.dxsm )
      // 	{
      // 	  pathopt.PredictUniformPathSims( nsim, new_rcs, new_fcs, new_ts );
      // 	  if ( opts.curit % 2 == 1 )
      // 	    {
      // 	      ExtendSimDisplacements
      // 		( pathopt.mPaths.mIters[0].mPath,
      // 		  pathopt.mPaths.mIters.back().mPath,
      // 		  fes, nsim, new_rcs.data() );
      // 	    }
      // 	  else
      // 	    {
      // 	      double dt = 1./(nsim-1.);
      // 	      double const shifts[4] = { 0., 1./3., 0., -1./3. };
      // 	      double shift = shifts[ (opts.curit/2) % 4 ];
      // 	      for ( std::size_t i=0; i<nsim; ++i )
      // 		{
      // 		  double t = std::max(0.,std::min(1.,(i+shift)*dt));
      // 		  pathopt.mPaths.mIters.back().mPath.GetValue( t, new_rcs.data() + i*ndim );
      // 		}
      // 	    }
	  
      // 	  pathopt.mPaths.mIters.back().mSims
      // 	    = ndfes::PathSims(ndim,nsim,new_rcs.data(),new_fcs.data());
      // 	  pathopt.mSimPaths.mIters.back() = pathopt.mPaths.mIters.back();

      // 	}
      else
	{
	  pathopt.PredictPathSims( nsim, new_rcs, new_fcs, new_ts );
	}
    }
  else
    {
      pathopt.PredictUniformPathSims( nsim, new_rcs, new_fcs, new_ts );
    }

  
  std::cout << "\nFinished creating the new string\n" << std::endl;

  

  ndfes::PathOptions popts( simdata->mPopts );
  
  double const RMSDTOL = 0.005;

  double distolgood = 0;
  double angtolgood = 0;
  double rmsdgood = 0;
  bool maxdispisbig = false;
  bool avgdispisbig = false;
  std::ofstream cfh;
  std::vector<std::ostream *> couts;
  {
    couts.push_back( &std::cout );
    if ( not opts.dry_run )
      {
	std::string fname = curdir + "/pathconv.txt";
	cfh.open( fname.c_str() );
	couts.push_back( &cfh );
      }

    ndfes::FTSMSlopeTest
      ( pathopt.mSimPaths,
	opts.mlen, opts.distol, opts.angtol, RMSDTOL,
	simdata->mTidxIsAngle,
	distolgood, angtolgood, rmsdgood, maxdispisbig, avgdispisbig,
	couts );
  }

  



  
  bool isapath = true;

  if ( opts.sasm  or opts.gbsm )
    {
      isapath = false;
    }
  
  if ( opts.gbsm )
    {
      isapath = false;
      
      std::vector<double> mynobs(nsim,simdata->mNumSamples);
      
      
      std::vector<double> prcs;
      prcs = ReadPrevUnoccSims( prevdir, fes );


      
      // std::vector<std::size_t> simidxs 
      // 	( ndfes::SamplingConvTest
      // 	  ( ndim, nsim, ospl, fes, opts.nlayers ) );

      // bool thisconv = ndfes::SamplingConvResult(simidxs,opts.conv_samples);
      
      bool isrefinement = ( distolgood < opts.distol ) and ( angtolgood < opts.angtol );

      //std::vector<ndfes::AreaSummary> simidxs;
      int reqlayers = opts.nlayers;
      {
	
	//bool myfound = false;
	for ( int layers=0; layers < std::min(opts.nlayers+1,2); ++layers )
	  {
	    std::vector<ndfes::AreaSummary> mysimidxs 
	      ( ndfes::SamplingConvTest
		( ndim, nsim, ospl, fes, layers ) );
	    
	    int mysamples = std::max( 1, 100/(layers+1) );
	    ndfes::SamplingConvPrint
	      ( mysimidxs, layers,
		mysamples, std::cout );
	      
	    // bool myconv = ndfes::SamplingConvResult(mysimidxs,100);
	    // if ( (not myconv) and (not myfound) )
	    //   {
	    // 	reqlayers = layers;
	    // 	myfound = true;
	    //   }
	    // simidxs = mysimidxs;
	  }
      }
      
      reqlayers = opts.curit % (opts.nlayers+1);
      


      std::vector<std::size_t> modfcidxs;
      if ( isrefinement and opts.predfc )
	{
	  pathopt.PredictPathSims( nsim, new_rcs, new_fcs, new_ts );
	}
      else
	{

	  ndfes::PathOptions tmpopts(popts);
	  if ( ! avgdispisbig )
	    {
	      for ( std::size_t icout=0; icout<couts.size(); ++icout )
		{
		  std::cout << "Setting minfc to ";
		  for ( std::size_t dim=0; dim<ndim; ++dim )
		    {
		      std::cout << std::fixed << std::setprecision(2) << std::setw(7) << simdata->mDefFC[dim];
		    }
		  std::cout << " because the average displacement is small\n";
		}
	      
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  tmpopts.minfc[dim] = simdata->mDefFC[dim];
		}
	    }
	  
	  
	  std::vector<double> reservedpts; //( ReadReservations(nextdir) );

	  pathopt.mPaths.mIters.back().PredictGridBasedSims
	    ( nsim, fes, tmpopts, reqlayers, reservedpts );
	  pathopt.mSimPaths.mIters.back() = pathopt.mPaths.mIters.back();
	  new_rcs = pathopt.mPaths.mIters.back().mSims.mRCs;
	  new_fcs = pathopt.mPaths.mIters.back().mSims.mFCs;
	    
	  modfcidxs = CheckUnoccSims( fes, simdata->mPopts, new_rcs, new_fcs, prcs );
	  
	  // if ( not opts.dry_run )
	  //   {
	  //     for ( std::size_t i=0; i<new_rcs.size(); ++i )
	  // 	{
	  // 	  reservedpts.push_back( new_rcs[i] );
	  // 	}
	  //     WriteReservations( nextdir, ndim, reservedpts );
	  //   }
	    
	}


      if ( not opts.dry_run )
	{
	  WritePrevUnoccSims( curdir, ndim, prcs );
	}
      
    }
  

  if ( not opts.dry_run )
    {
      std::cout << std::endl
		<< "Preparing simulations in: " << nextdir << "\n"
		<< std::endl;

      int firstit = 0;

      
      if ( opts.neqit >= 0 )
	{
	  if ( opts.curit > opts.neqit )
	    {
	      firstit = std::min( opts.neqit + 1, opts.curit - opts.neqit );
	    }
	}

      //firstit = std::max( firstit, std::max( 0, opts.curit - std::max(2,opts.sim_layers) ) );
      //firstit = std::min( firstit, std::max( 0, opts.curit - std::max(1,opts.nlayers/2) ) );
      firstit = std::max( firstit, std::max( 0, opts.curit - std::max(2,opts.nlayers) ) );

      
      //std::cout << "firstit = " << firstit << "\n";
      
      // ndfes::CopyFiles( ndim, simdata->mNumSim,
      // 			opts.disang, opts.mdin,
      // 			opts.curit, nextdir, opts.prefix,
      // 			opts.extra_prefixes,
      // 			opts.cp_nearest, not isapath, firstit,
      // 			new_rcs.data(), new_fcs.data() );

      ndfes::CopyFilesNextIter
	( ndim, nsim,
	  opts.disang, opts.mdin,
	  opts.curit, nextdir, opts.prefix,
	  opts.extra_prefixes,
	  opts.cp_nearest, not isapath, firstit, opts.pad,
	  new_rcs.data(), new_fcs.data() );

      
    }


  if ( not opts.dry_run )
    {
      std::stringstream sout;
      sout << curdir;
      if ( opts.prefix.size() > 0 )
	{
	  sout << "/" << opts.prefix;
	}
      
      //ndfes::SavePathPickles( ndim, sout.str() + "/path.micro.pkl", pathopt.GetMicroIters(10) );
      //ndfes::SavePathPickles( ndim, sout.str() + "/path.sims.pkl", pathopt.simpaths );

      pathopt.GetMicroIters(10).WriteXml( sout.str() + "/path.micro.xml" );
      pathopt.mSimPaths.WriteXml( sout.str() + "/path.sims.xml" );
      //if ( penopt )
      //{
      //penopt->GetMicroIters(10).WriteXml( sout.str() + "/path.pen.xml" );
      // }	 
    }

  
  if ( fes and not opts.dry_run )
    {
      std::stringstream fname;
      fname << curdir << "/analysis/";
      if ( opts.prefix.size() > 0 )
	{
	  fname << "path." << opts.prefix << ".dat";
	}
      else
	{
	  fname << "path.dat";
	}
      
      std::cout << std::endl
		<< "Writing 1D FES Projection to " << fname.str()
		<< std::endl;
          std::ofstream cout;
    cout.open( fname.str().c_str() );

#define FMTF std::fixed << std::setw(14) << std::setprecision(8)
#define FMTH std::fixed << std::setw(9) << std::setprecision(3)
#define FMTI std::fixed << std::setw(13) << std::setprecision(3)
#define FMTE std::scientific << std::setw(15) << std::setprecision(6)
#define FMTB std::scientific << std::setw(23) << std::setprecision(14)

    std::size_t nspl = ospl->mT.size();

    std::vector<double> splpts( ospl->mX );
    ndfes::WrapPath( fes->mDimInfo, splpts );
    
    for ( std::size_t i=0; i<nspl; ++i )
      {
	cout << std::setw(3) << i+1
	     << FMTB << ospl->mT[i];
	
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    cout << FMTB << ospl->mX[dim+i*ndim];
	  }
	
	double val=0;
	double err=0;
	double re=0;
	std::size_t ns=0;
	fes->GetValueErrorAndEntropy(splpts.data() + i*ndim,val,err,re,ns);
	cout << FMTE << val
	     << FMTI << err
	     << FMTH << re
	     << std::setw(7) << ns
	     << "\n";
      }
#undef FMTE
#undef FMTF
#undef FMTH
#undef FMTI
    cout.close();

    }
  


  return 0;
  
}

  



/////////////////////////////////////////////////////////////////////////////







int main( int argc, char * argv[] )
{
  return main_new( argc, argv );
}



