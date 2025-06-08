#define _USE_MATH_DEFINES
#include <cmath>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <random>

#include "Sim.hpp"
#include "Stats.hpp"
#include "StatsMod.hpp"

#include "ParaRNG.hpp"


int ReadDat( std::string fname,
	     std::vector<double> & es )
{
  es.resize(0);
  es.reserve(1000);
  int nskip = 0;

  std::ifstream cin;
  cin.open(fname);

  if ( cin.good() )
    {
      double xprev = -1.;
      /*
      std::string line;
      while ( std::getline(cin,line) )
	{
	  std::istringstream iss(line);
	  double xval,yval;
	  if ( iss >> xval >> yval )
	    {
	      if ( std::abs( xval-xprev ) < 0.0001 )
		{
		  ++nskip;
		}
	      else
		{
		  xprev=xval;
		  es.push_back( yval );
		}
	    }
	}
      */

      
      double xval,yval;
      while ( cin >> xval >> yval )
	{
	  if ( std::abs( xval-xprev ) < 0.0001 )
	    {
	      ++nskip;
	    }
	  else
	    {
	      xprev=xval;
	      es.push_back( yval );
	    }
	}
      

      /*
      cin.seekg(0, std::ios::end);
      std::size_t length = cin.tellg();
      std::string buffer;
      if ( (int)length > 0 )
	{
	  buffer.resize(static_cast<std::string::size_type>(length));
	  cin.seekg(0);
	  cin.read(&buffer.front(), length);
	  std::stringstream ss(buffer);
	  double xval,yval;
	  while ( ss >> xval >> yval )
	    {
	      if ( std::abs( xval-xprev ) < 0.0001 )
		{
		  ++nskip;
		}
	      else
		{
		  xprev=xval;
		  es.push_back( yval );
		}
	    }
	}
      */
      
    }
  else
    {
      std::cerr << "Failed to open " << fname << "\n";
      std::exit(EXIT_FAILURE);
    }

  /*
  if ( (int)es.size() == 0 )
    {
      std::cerr << "File contained no data " << fname << "\n";
      std::exit(EXIT_FAILURE);
    }
  */

  return nskip;
}


edgembar::Sim::Sim()
  : SimIdx(0),
    Beta(0),
    AutoEqMode(1),
    NumSamples(0),
    OrigNumSamples(0),
    LocalSimIdx(0),
    CurStride(0),
    OrigStride(0),
    ProdStart(0),
    ProdStride(0),
    IsConverged(true)
{
  
}

edgembar::Sim::Sim( int const simidx )
  : SimIdx(simidx),
    Beta(0),
    AutoEqMode(1),
    NumSamples(0),
    OrigNumSamples(0),
    LocalSimIdx(0),
    CurStride(0),
    OrigStride(0),
    ProdStart(0),
    ProdStride(0),
    IsConverged(true)
{
  
}


edgembar::Sim::Sim
( int const simidx,
  std::vector<int> eneidxs,
  std::vector<std::string> const & statelabels,
  std::string const datadir,
  double const beta,
  int const autoeqmode )
  : SimIdx(simidx),
    EneIdxs(eneidxs),
    StateLabels(statelabels),
    DataDir(datadir),
    Beta(beta),
    AutoEqMode(autoeqmode),
    NumSamples(0),
    OrigNumSamples(0),
    LocalSimIdx(0),
    AvgEnes(eneidxs.size(),0),
    CurStride(0),
    OrigStride(0),
    ProdStart(0),
    ProdStride(0),
    IsConverged(true)
{

  for ( int i=0, n=EneIdxs.size(); i<n; ++i )
    {
      if ( EneIdxs[i] == SimIdx )
	{
	  LocalSimIdx = i;
	}
    }
}


int edgembar::Sim::FindLocalIdx( int const eneidx ) const
{
  /*
  std::printf("Looking for %3i ",eneidx);
  std::printf("have:");
  for ( int i=0, n=EneIdxs.size(); i<n; ++i )
    {
      std::printf(" %3i",EneIdxs[i]);
    }
  std::printf("\n");
  */
  
  int idx = -1;
  for ( int i=0, n=EneIdxs.size(); i<n; ++i )
    {
      if ( eneidx == EneIdxs[i] )
	{
	  idx = i;
	  break;
	}
    };
  //std::printf("Returning %3i\n",idx);
  return idx;
}

int edgembar::Sim::FindLocalIdx( edgembar::Sim const * other ) const
{
  return FindLocalIdx( other->SimIdx );
}




void edgembar::Sim::ReadFiles( double const fstart, double const fstop, int const stride )
{

  std::vector< std::vector<double> > data;

  //std::cout << "EneIdxs.size() " << EneIdxs.size() << "\n";
  int ns = EneIdxs.size();

  for ( int i=0; i<ns; ++i )
    {

      std::stringstream sname;
      sname << DataDir << "/efep_" << StateLabels[SimIdx]
	    << "_" << StateLabels[EneIdxs[i]] << ".dat";

      std::vector<double> vals;
      ReadDat( sname.str(), vals );

      for ( int k=0, nk=vals.size(); k<nk; ++k )
	{
	  vals[k] *= Beta;
	}
      
      if ( NumSamples == 0 )
	{
	  NumSamples = vals.size();
	}
      else if ( (int)vals.size() != NumSamples )
	{
	  std::cerr << "File " << sname.str() << " has "
		    << vals.size() << " samples, but "
		    << DataDir << "/efep_" << StateLabels[SimIdx]
		    << "_" << StateLabels[0] << ".dat"
		    << " had " << NumSamples << " samples\n";
	  std::exit(EXIT_FAILURE);
	}
      data.push_back(vals);
    };

  int istart = std::min(NumSamples,std::max(0,(int)(NumSamples*fstart+0.5)));
  int istop = std::min(NumSamples,std::max(istart+1,(int)(NumSamples*fstop+0.5)));
  //std::printf("%5i %5i %5i\n",NumSamples,istart,istop);
  for ( std::size_t i=0, n=data.size(); i<n; ++i )
    {
      std::vector<double> tmp;
      for ( int ii=istart; ii<istop; ii += stride )
	{
	  tmp.push_back( data[i][ii] );
	}
      data[i] = tmp;
      //data[i].erase( data[i].begin() + istop, data[i].end() );
      //data[i].erase( data[i].begin(), data[i].begin() + istart );
    }
  
  NumSamples = 0;
  if ( ns > 0 )
    {
      NumSamples = data[0].size();
    }
  
  
  Emat.assign(ns*NumSamples,0.);
  
  for ( int k=0; k<NumSamples; ++k )
    {
      for ( int i=0; i<ns; ++i )
	{
	  Emat[i+k*ns] = data[i][k]; // - data[LocalSimIdx][k];
	}
    }
  OrigNumSamples = NumSamples;
  OrigStride = -1;
  CurStride = OrigStride;


  if ( ns > 0 )
    {
      std::stringstream sname;
      sname << DataDir << "/dvdl_" << StateLabels[SimIdx] << ".dat";
      
      std::ifstream cin;
      cin.open(sname.str().c_str());
      if ( cin.good() )
	{
	  ReadDat( sname.str(), DVDL );
	  
	  int n = DVDL.size();
	  if ( n == 0 )
	    {
	      std::cerr << "Ignoring DVDL data because " << sname.str()
			<< " has " << DVDL.size()
			<< " values" << "\n";
	      DVDL.resize(0);
	    }
	  istart = std::min(n,std::max(0,(int)(n*fstart+0.5)));
	  istop = std::min(n,std::max(istart+1,(int)(n*fstop+0.5)));
	  std::vector<double> tmp;
	  for ( int ii=istart; ii<istop; ii += stride )
	    {
	      tmp.push_back( DVDL[ii] );
	    }
	  DVDL = tmp;
	};
    }

  
  /*
  MaxZ.assign(NumSamples,0.);
  Mdat.assign( ns*NumSamples, 0. );

  for ( int k=0; k<NumSamples; ++k )
    {
      double maxz = -1.e+30;
      for ( int i=0; i<ns; ++i )
	{
	  double z = -Emat[i+k*ns];
	  maxz = std::max(maxz,z);
	}
      MaxZ[k] = maxz;
      for ( int i=0; i<ns; ++i )
	{
	  double z = -Emat[i+k*ns];
	  Mdat[i+k*ns] = std::exp(z-maxz);
	}
    };
  

  OrigNumSamples = NumSamples;
  OrigStride = CptStride();
  CptAutoEquil( ProdStart, ProdStride, IsConverged );
  CurStride = OrigStride;
  */
}


void edgembar::Sim::PrecomputeExps( double const ptol )
{
  int ns = EneIdxs.size();
  MaxZ.assign(NumSamples,0.);
  Mdat.assign( ns*NumSamples, 0. );
  
  std::vector<double> myavg(ns,0);
  for ( int k=0; k<NumSamples; ++k )
    {
      double maxz = -1.e+30;
      for ( int i=0; i<ns; ++i )
	{
	  Emat[i+k*ns] -= AvgEnes[i];
	  double z = -Emat[i+k*ns];
	  maxz = std::max(maxz,z);
	  myavg[i] += Emat[i+k*ns];
	}
      MaxZ[k] = maxz;
      //std::printf("maxz %5i %13.4e\n",k,MaxZ[k]);
      for ( int i=0; i<ns; ++i )
	{
	  double z = -Emat[i+k*ns];
	  Mdat[i+k*ns] = std::exp(z-maxz);
	}
    };
  // for ( int i=0; i<ns; ++i )
  //   {
  //     myavg[i] /= NumSamples;
  //     std::printf("%3i %13.4e %13.4e\n",i,myavg[i],AvgEnes[i]);
  //   }
  
  OrigNumSamples = NumSamples;
  OrigStride = CptStride();

  CptAutoEquil( ProdStart, ProdStride, IsConverged, ptol );
  CurStride = OrigStride;
}



int edgembar::Sim::CptStride() const
{
  std::vector<double> data(NumSamples);
  int const ns = EneIdxs.size();
  for ( int i=0; i<NumSamples; ++i )
    {
      data[i] = Emat[LocalSimIdx + i*ns];
    };
  
  //int s0 = edgembar::CptSampleStride(data.size(),data.data());
  int s0 = 1;
  int sp = s0;
  int sm = s0;
  if ( LocalSimIdx > 0 )
    {
      for ( int i=0; i<NumSamples; ++i )
	{
	  data[i] = Emat[(LocalSimIdx-1) + i*ns] - Emat[LocalSimIdx + i*ns];
	};
      sm = ccdl::CptSampleStride(data.size(),data.data());
    }
  if ( LocalSimIdx < ns-1 )
    {
      for ( int i=0; i<NumSamples; ++i )
	{
	  data[i] = Emat[(LocalSimIdx+1) + i*ns] - Emat[LocalSimIdx + i*ns];
	};
      sp = ccdl::CptSampleStride(data.size(),data.data());
    }
  return std::max(s0,std::max(sp,sm));
}


//void edgembar::Sim::StoreAutoEquil()
//{
//  CptAutoEquil( ProdStart, ProdStride, IsConverged );
//}


void edgembar::Sim::CptAutoEquil( int & start, int & stride, bool & isconv, double const ptol ) const
{

  double maxeq = 0.75;
  int const ns = EneIdxs.size();

  int pm = -1;
  int sm = -1;
  int pp = -1;
  int sp = -1;
  bool cm = true;
  bool cp = true;
  
  std::vector<double> datam(NumSamples);

  int nblocks = 20;
  
  if ( LocalSimIdx > 0 )
    {
      //std::ofstream cout;
      //cout.open("foo.dat");
      //std::cout << "REV " << LocalSimIdx << std::endl;
      double avgs = AvgEnes[LocalSimIdx-1] - AvgEnes[LocalSimIdx]; 
      for ( int i=0; i<NumSamples; ++i )
	{
	  datam[i] = Emat[(LocalSimIdx-1) + i*ns] - Emat[LocalSimIdx + i*ns] + avgs;
	  //cout << i << " " << std::scientific << std::setw(20)
	  //   << std::setprecision(10) << datam[i] << "\n";
	};
      //cout.close();
      
      //cm = edgembar::CptSampleStartAndStride(datam.size(),datam.data(),maxeq,pm,sm);
      //std::printf("datam.size(),nblocks %lu %i\n",datam.size(),nblocks);
      //std::printf("Rev\n");

      if ( AutoEqMode == 1 )
	{
	  cm = edgembar::CptSampleStartAndStrideByBlock_v3(datam.size(),datam.data(),maxeq,nblocks,pm,sm,ptol);
	}
      else
	{
	  std::cerr << "Invalid AutoEqMode " << AutoEqMode << " in edgembar::Sim::CptAutoEquil\n";
	}

      
      //std::printf("foo.dat start,stride %i %i\n",pm,sm);
      //std::exit(EXIT_SUCCESS);

    }

  std::vector<double> datap(NumSamples);

  if ( LocalSimIdx < ns-1 )
    {
      std::ofstream cout;
      //cout.open("foo.dat");

      //std::cout << "FWD " << LocalSimIdx << std::endl;
      double avgs = AvgEnes[LocalSimIdx+1] - AvgEnes[LocalSimIdx]; 
      for ( int i=0; i<NumSamples; ++i )
	{
	  datap[i] = Emat[(LocalSimIdx+1) + i*ns] - Emat[LocalSimIdx + i*ns] + avgs;
	  
	  //cout << i << " " << std::scientific << std::setw(20)
	  // << std::setprecision(10) << datap[i] << "\n";
	};
      //cout.close();
      //std::printf("Fwd\n");
      //cp = edgembar::CptSampleStartAndStride(datap.size(),datap.data(),maxeq,pp,sp);
      if ( AutoEqMode == 1 )
	{
	  cp = edgembar::CptSampleStartAndStrideByBlock_v3(datap.size(),datap.data(),maxeq,nblocks,pp,sp,ptol);
	}
      else
	{
	  std::cerr << "Invalid AutoEqMode " << AutoEqMode << " in edgembar::Sim::CptAutoEquil\n";
	}

    }

  start = std::max(pp,pm);

  if ( LocalSimIdx > 0 )
    {
      datam.erase( datam.begin(), datam.begin() + start );
      sm = ccdl::CptSampleStride(datam.size(),datam.data());
    }
  
  if ( LocalSimIdx < ns-1 )
    {
      datap.erase( datap.begin(), datap.begin() + start );
      sp = ccdl::CptSampleStride(datap.size(),datap.data());
    }
  stride = std::max(sp,sm);

  isconv = true;
  if ( ! cp )
    {
      isconv = false;
    }
  if ( ! cm )
    {
      isconv = false;
    };
}



void edgembar::Sim::OnlyKeepRange( double const start, double const stop )
{
  int istart = std::max(0, (int)(start * NumSamples+0.5));
  int istop = std::min(NumSamples, (int)(stop * NumSamples+0.5));
  int ns = EneIdxs.size();

  int ndvdl = DVDL.size();
  if ( ndvdl > 0 )
    {
      int dStart = std::max(0, (int)(start * ndvdl+0.5));
      int dStop  = std::min(ndvdl, (int)(stop * ndvdl+0.5));
      DVDL.erase( DVDL.begin()+dStop, DVDL.end() );
      DVDL.erase( DVDL.begin(), DVDL.begin()+dStart );
    }
  
  MaxZ.erase( MaxZ.begin()+istop, MaxZ.end() );
  MaxZ.erase( MaxZ.begin(), MaxZ.begin()+istart );
  
  Emat.erase( Emat.begin()+istop*ns, Emat.end() );
  Emat.erase( Emat.begin(), Emat.begin()+istart*ns );
  if ( (int)Emat.size() != (istop-istart)*ns )
    {
      std::cerr << "Programming error in edgembar::Sim::OnlyKeepRange : "
		<< "Emat Size inconsistency " << Emat.size() << " "
		<< istart << " " << istop << " " << ns << std::endl;
      std::exit(EXIT_FAILURE);
    }
  
  Mdat.erase( Mdat.begin()+istop*ns, Mdat.end() );
  Mdat.erase( Mdat.begin(), Mdat.begin()+istart*ns );
  if ( (int)Mdat.size() != (istop-istart)*ns )
    {
      std::cerr << "Programming error in edgembar::Sim::OnlyKeepRange : "
		<< "Mdat Size inconsistency " << Mdat.size() << " "
		<< istart << " " << istop << " " << ns << std::endl;
      std::exit(EXIT_FAILURE);
    }

  
  NumSamples = istop-istart;
  CurStride = CptStride();
}


namespace edgembar
{
  std::vector<double> StridedCopyVec
  ( int const start,
    int const stop,
    int const stride,
    double const * vin );
  
  std::vector<double> StridedCopyMat
  ( int const start,
    int const stop,
    int const stride,
    int const nfast,
    double const * vin );


  std::vector<double> CopyVec
  ( int const start,
    int const stop,
    double const * vin );
  
  std::vector<double> CopyMat
  ( int const start,
    int const stop,
    int const stride,
    double const * vin );

  
}


std::vector<double>
edgembar::StridedCopyVec
( int const start,
  int const stop,
  int const stride,
  double const * vin )
{
  std::vector<double> vout;
  int s = std::max(1,stride);
  for ( int i=start; i<stop; i += s )
    {
      vout.push_back( vin[i] );
    }
  return vout;
}

std::vector<double>
edgembar::CopyVec
( int const start,
  int const stop,
  double const * vin )
{
  std::vector<double> vout( vin+start, vin+stop );
  return vout;
}


std::vector<double>
edgembar::StridedCopyMat
( int const start,
  int const stop,
  int const stride,
  int const nfast,
  double const * vin )
{
  std::vector<double> vout;
  int s = std::max(1,stride);
  //std::cout << start << " " << stop << " " << nfast << " " << s << "\n";
  for ( int i=start; i<stop; i += s )
    {
      for ( int j=0; j<nfast; ++j )
	{
	  vout.push_back( vin[j+i*nfast] );
	};
    }
  return vout;
}


std::vector<double>
edgembar::CopyMat
( int const start,
  int const stop,
  int const nfast,
  double const * vin )
{
  std::vector<double> vout;
  //std::cout << start << " " << stop << " " << nfast << " " << s << "\n";
  for ( int i=start; i<stop; ++i )
    {
      for ( int j=0; j<nfast; ++j )
	{
	  vout.push_back( vin[j+i*nfast] );
	};
    }
  return vout;
}


  


void edgembar::Sim::OnlyKeepProd()
{  
  int istart = ProdStart;
  int istop = NumSamples;
  int ns = EneIdxs.size();

  if ( NumSamples > 0 and istart >= 0 )
    {
      //MaxZ = edgembar::StridedCopyVec(istart,istop,ProdStride,MaxZ.data());
      //Emat = edgembar::StridedCopyMat(istart,istop,ProdStride,ns,Emat.data());
      //Mdat = edgembar::StridedCopyMat(istart,istop,ProdStride,ns,Mdat.data());
      //CurStride = 1;
      int ndvdl = DVDL.size();
      if ( ndvdl > 0 )
	{
	  double start = ((double)istart)/((double)istop);
	  int dStart = std::max(0, (int)(start * ndvdl+0.5));
	  DVDL = edgembar::CopyVec(dStart,ndvdl,DVDL.data());
	}
      MaxZ = edgembar::CopyVec(istart,istop,MaxZ.data());
      Emat = edgembar::CopyMat(istart,istop,ns,Emat.data());
      Mdat = edgembar::CopyMat(istart,istop,ns,Mdat.data());
      CurStride = ProdStride;

      NumSamples = Emat.size() / ns;
    };
}



// void edgembar::Sim::PrepareMBAR()
// {
//   int ns = EneIdxs.size();
//   Mdat.assign( ns*NumSamples, 0. );
//   for ( int k=0; k<NumSamples; ++k )
//     {
//       double maxz = -1.e+30;
//       for ( int i=0; i<ns; ++i )
// 	{
// 	  double z = -Emat[i+k*ns];
// 	  maxz = std::max(maxz,z);
// 	}
//       for ( int i=0; i<ns; ++i )
// 	{
// 	  double z = -Emat[i+k*ns];
// 	  Mdat[i+k*ns] = std::exp(z-maxz);
// 	}
//     };
// }



void edgembar::Sim::Bootstrap()
{
  //std::vector<double> D(DVDL);
  std::vector<double> Z(MaxZ);
  std::vector<double> E(Emat);
  std::vector<double> M(Mdat);
  int ns = EneIdxs.size();
  std::uniform_int_distribution<int> dist(0,NumSamples-1);

  int k=0;
  while ( k < NumSamples )
    {
      int o = ParaRNG::RandInt(dist);
      for ( int u=0; u<CurStride; ++u )
	{
	  if ( k < NumSamples )
	    {
	      for ( int e=0; e<ns; ++e )
		{
		  M[e+k*ns] = Mdat[e+o*ns];
		  E[e+k*ns] = Emat[e+o*ns];
		}
	      //if ( DVDL.size() > 0 )
	      //{
	      //  D[k] = DVDL[o];
	      //}
	      Z[k] = MaxZ[o];
	      //std::printf("%5i / %5i : %5i\n",k,NumSamples,o);
	      o = (o+1) % NumSamples;
	      k += 1;
	    };
	}
    };
  //DVDL = D;
  Mdat = M;
  Emat = E;
  MaxZ = Z;
}



void edgembar::Sim::WriteDebugInfo( std::ostream & cout ) const
{
  if ( (int) StateLabels.size() > 0 )
    {
      cout << std::left << std::setw(20) << StateLabels[SimIdx]
       << " " << DataDir
       << "\n"
       << std::setw(6) << ""
       << "SimIdx: "
       << std::setw(3) << SimIdx
       << " LocalSimIdx: "
       << std::setw(3) << LocalSimIdx
       << " EneIdxs: ";
      for ( int i=0, n=EneIdxs.size(); i<n; ++i )
	{
	  cout << std::setw(4) << EneIdxs[i];
	};
      
      cout << "\n"
	   << std::setw(6) << ""
	   << "NumSamples: "
	   << std::setw(6) << NumSamples
	   << " CurStride: "
	   << std::setw(4) << CurStride
	   << " OrigStride: "
	   << std::setw(4) << OrigStride
	   << " ProdStart: "
	   << std::setw(5) << ProdStart
	   << " ProdStride: "
	   << std::setw(4) << ProdStride
	   << " IsConverged: "
	   << IsConverged
	   << "\n";
    };
}
