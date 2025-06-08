#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <algorithm>

#include "PathData.hpp"
#include "Smoothing.hpp"
#include "BiasedMin.hpp"
#include "LinearInterp.hpp"
#include "Tube.hpp"
#include "Utils.hpp"
#include "PredictCentroids.hpp"


///////////////////////////////////////////////////////////////////

namespace
{
  std::vector<double>
  GetScaledForceConstants
  ( std::size_t const ndim,
    std::size_t nsim,
    double const * minfcs,
    double const * maxfcs,
    double const * fc0,
    double const * ws )
  {
    std::vector<double> fcs(ndim*nsim,0);
    for ( std::size_t i=0; i<nsim; ++i )
      {
	double w = std::min(1.,std::max(-1.,ws[i]));
	if ( w > 0. )
	  {
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		fcs[dim+i*ndim] = (1.-w)*fc0[dim+i*ndim] + w*maxfcs[dim];
	      }
	  }
	else
	  {
	    w = std::abs(w);
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		fcs[dim+i*ndim] = (1.-w)*fc0[dim+i*ndim] + w*minfcs[dim];
	      }
	  }
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    fcs[dim+i*ndim] = std::min(maxfcs[dim],
				       std::max(minfcs[dim],
						fcs[dim+i*ndim]));
	  }
      }
    return fcs;
  }


  void CalcDisplacementPercentages
  ( std::size_t const ndim,
    std::size_t const onsim,
    double * percs,
    double * cosas,
    std::vector<bool> const & done,
    double const * test_rcs,
    double const * test_fcs,
    double const * dws,
    std::shared_ptr<ndfes::FES> fes,
    ndfes::PathOptions const & popts )
  {
    std::fill( percs, percs + onsim, popts.tdispfc );
    std::fill( cosas, cosas + onsim, 0. );

    std::vector<double> trcs( test_rcs, test_rcs+ndim*onsim );
    std::vector<double> tfcs( test_fcs, test_fcs+ndim*onsim );
    std::vector<double> eff_means
      ( ndfes::PredictCentroids
	( ndim, onsim, trcs, tfcs, done, fes, popts ) );

    
    std::vector<double> uw(ndim,0);
    std::vector<double> dx(ndim,0);
    for ( std::size_t i=0; i<onsim; ++i )
      {
	double lw = 0;
	double lx = 0;
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    lw += dws[dim+i*ndim] * dws[dim+i*ndim];
	    dx[dim] = eff_means[dim+i*ndim] - test_rcs[dim+i*ndim];
	    lx += dx[dim] * dx[dim];
	  }
    
	lw = std::sqrt(lw);
	lx = std::sqrt(lx);
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    uw[dim] = dws[dim+i*ndim] / lw;
	  }
	if ( lx > 0.001 )
	  {
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		cosas[i] += uw[dim] * dx[dim] / lx;
	      };
	    percs[i] = std::abs(cosas[i]) * (lx/lw);
	    // std::printf("perc %2lu %9.4f %9.4f [",i,percs[i],cosas[i]);
	    // for ( std::size_t dim=0; dim<ndim; ++dim )
	    //   {
	    //     std::printf("%9.4f",dx[dim]);
	    //   }
	    // std::printf("]%9.4f %9.4f\n",lx,lw);
	  }
      }
  }



  std::vector<double> GetUniformPathDeltas
  ( std::size_t const ndim,
    std::size_t const nsim,
    double const * rcs )
  {
    std::vector<double> dws(ndim*nsim,0);
    for ( std::size_t i=0; i<nsim; ++i )
      {
	if ( i == 0 )
	  {
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		dws[dim+i*ndim] = rcs[dim+(i+1)*ndim]-rcs[dim+(i)*ndim];
	      }
	  }
	else if ( i == nsim-1 )
	  {
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		dws[dim+i*ndim] = rcs[dim+(i)*ndim]-rcs[dim+(i-1)*ndim];
	      }
	  }
	else
	  {
	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		dws[dim+i*ndim] = 0.5*(rcs[dim+(i+1)*ndim]-rcs[dim+(i-1)*ndim]);
	      }
	  }
      }
    return dws;
  }





  // std::vector<double> ReadReservations
  // ( std::string curdir )
  // {
  //   std::vector<double> vs;
    
  //   std::string fname = curdir + "/reservations.txt";
  //   if ( FileExists(fname) )
  //     {
  // 	std::cout << "Reading " << fname << std::endl;
	
  // 	std::ifstream cin;
  // 	cin.open( fname.c_str() );
	
  // 	std::string line;
  // 	while ( std::getline(cin,line) )
  // 	  {
  // 	    std::stringstream sline(line);
  // 	    double v;
  // 	    while ( sline >> v )
  // 	      {
  // 		vs.push_back(v);
  // 	      }
  // 	  };
  //     }
  //   return vs;
  // }

  // void WriteReservations
  // ( std::string curdir, std::size_t const ndim, std::vector<double> vs )
  // {
  //   std::string fname = curdir + "/reservations.txt";
  //   std::cout << "Writing " << fname << std::endl;
  //   std::ofstream cout;
  //   cout.open( fname.c_str() );
  //   std::size_t npts = vs.size()/ndim;
  //   for ( std::size_t i=0; i<npts; ++i )
  //     {
  // 	for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	  {
  // 	    cout << std::scientific << std::setw(19)
  // 		 << std::setprecision(10) << vs[dim+i*ndim];
  // 	  }
  // 	cout << "\n";
  //     }
  //   cout.close();
  // }

  
}




///////////////////////////////////////////////////////////////////


ndfes::PathSpl::PathSpl()
  : mSType(0),
    mSmooth(true),
    mNumDim(0),
    mNumPts(0)
{}


ndfes::PathSpl::PathSpl
( int const stype, bool const smooth,
  std::size_t const ndim, std::size_t const npts, 
  double const * ipts )
  : mSType(stype),
    mSmooth(smooth),
    mNumDim(ndim),
    mNumPts(npts),
    mPts(ipts,ipts+npts*ndim)
{
  ResetSpl();
}

ndfes::PathSpl::PathSpl
( int const stype, bool const smooth,
  std::size_t const ndim, std::size_t const npts, 
  double const * ipts, double const * ts )
  : mSType(stype),
    mSmooth(smooth),
    mNumDim(ndim),
    mNumPts(npts),
    mPts(ipts,ipts+npts*ndim),
    mT(ts,ts+npts)
{
  ResetSpl();
}


void ndfes::PathSpl::ResetSpl()
{
  std::vector<double> cpts(mPts);
  if ( mSmooth )
    {
      int smooth_nmin = 3;
      int smooth_nmax = 11;
      
      cpts = ndfes::GetIterWinAvg(mNumDim,mNumPts,cpts.data(),
                                  smooth_nmin, smooth_nmax);

    }

  bool akima = mSType;
  if ( mT.size() == mNumPts )
    {
      mSpl.reset( new ndfes::PCurve
		  ( mNumDim, mNumPts, cpts.data(),
		    mT.data(), akima ) );
    }
  else
    {
      int maxit = 100;
      int nseg = 1000;
      
      mSpl.reset( new ndfes::PCurve
		  ( mNumDim, mNumPts, cpts.data(),
		    akima, maxit, nseg ) );
      mT = mSpl->mT;
    }
}

std::vector<double> ndfes::PathSpl::GetControlPts() const
{
  return mSpl->mX;
}

std::vector<double>
ndfes::PathSpl::GetUniformPts( std::size_t const n ) const
{
  double const dt = 1. / (n-1.);
  std::vector<double> pts(n*mNumDim,0.);
  for ( std::size_t i=0; i<n; ++i )
    {
      double t = i * dt;
      mSpl->GetValue(t, pts.data() + i*mNumDim );
    }
  return pts;
}

void ndfes::PathSpl::GetValue( double const t, double * x ) const
{
  return mSpl->GetValue(t,x);
}
  

void ndfes::PathSpl::WriteXml( std::ostream & cout, int const ilev ) const
{
#define IND1 std::setw(2*(ilev+0)) << ""
#define IND2 std::setw(2*(ilev+1)) << ""
#define IND3 std::setw(2*(ilev+2)) << ""

  std::vector<double> cpts( GetControlPts() );
  std::vector<double> ts( mSpl->mT );
  
  cout << IND1 << "<path>\n";
  cout << IND2 << "<type>" << mSType << "</type>\n";
  //cout << IND2 << "<smooth>" << (int)mSmooth << "</smooth>\n";
  cout << IND2 << "<smooth>" << 0 << "</smooth>\n";
  for ( std::size_t i=0; i<mNumPts; ++i )
    {
      cout << IND2 << "<pt idx=\"" << i << "\">\n";
      cout << IND3
	   << "<t>"
	   << std::scientific << std::setw(23) << std::setprecision(14)
	   << ts[i]
	   << "</t>\n";
      
      for ( std::size_t dim=0; dim<mNumDim; ++dim )
	{
	  cout << IND3
	       << "<rc dim=\"" << dim << "\">"
	       << std::scientific << std::setw(21) << std::setprecision(12)
	    //<< mPts[dim+i*mNumDim]
	       << cpts[dim+i*mNumDim]
	       << "</rc>"
	       << "\n";
	}
      cout << IND2 << "</pt>\n";
    }
  cout << IND1 << "</path>\n";

#undef IND1
#undef IND2
#undef IND3
}


ndfes::PathSpl::PathSpl( XMLIO::xml_node node )
  : mSType(0),
    mSmooth(true),
    mNumDim(0),
    mNumPts(0)
{

  if ( not node )
    {
      std::cerr << "Failed to read PathSpl from xml" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  //std::cout << "In PathSpl, extract type" << std::endl;
  XMLIO::ExtractValue(node,"type",mSType);
  //std::cout << "In PathSpl, extract smooth" << std::endl;
  XMLIO::ExtractValue(node,"smooth",mSmooth);

  std::size_t nt = 0;
  
  for ( XMLIO::xml_node pnode = node.child("pt");
	pnode; pnode = pnode.next_sibling("pt") )
    {
      //std::cout << "In PathSpl, loop pt" << std::endl;

      ++mNumPts;

      XMLIO::xml_node tnode = pnode.child("t");
      if ( tnode )
	{
	  ++nt;
	}
      
      if ( mNumPts == 1 )
	{
	  mNumDim = 0;
	  for ( XMLIO::xml_node qnode = pnode.child("rc");
		qnode; qnode = qnode.next_sibling("rc") )
	    {
	      //std::cout << "In PathSpl, loop rc" << std::endl;
	      ++mNumDim;
	    }
	}
    }


  mPts.assign( mNumDim*mNumPts, 0. );
  if ( nt == mNumPts )
    {
      mT.assign( mNumPts, 0. );
    }
  
  for ( XMLIO::xml_node pnode = node.child("pt");
	pnode; pnode = pnode.next_sibling("pt") )
    {
      int idx = 0;
      //std::cout << "In PathSpl, idx" << std::endl;
      XMLIO::ExtractAttribute(pnode,"idx",idx);

      if ( mT.size() > 0 )
	{
	  XMLIO::ExtractValue(pnode,"t",mT[idx]);
	};
      
      
      for ( XMLIO::xml_node qnode = pnode.child("rc");
	    qnode; qnode = qnode.next_sibling("rc") )
	{
	  int dim = 0;
	  //std::cout << "In PathSpl, dim" << std::endl;
	  XMLIO::ExtractAttribute(qnode,"dim",dim);
	  //std::cout << "In PathSpl, value " << mPts.size() << " " << dim+idx*mNumDim << std::endl;

	  XMLIO::ExtractValue(qnode,mPts[dim+idx*mNumDim]);
	  //std::cout << "In PathSpl, got" << mPts[dim+idx*mNumDim] << std::endl;

	}
    }

  ResetSpl();
}


///////////////////////////////////////////////////////////////////



ndfes::PathSims::PathSims()
  : mNumDim(0),
    mNumSim(0)
{}

ndfes::PathSims::PathSims
( std::size_t const ndim, std::size_t const nsim,
  double const * prcs, double const * pfcs )
  : mNumDim(ndim),
    mNumSim(nsim),
    mRCs( prcs, prcs + ndim*nsim ),
    mFCs( pfcs, pfcs + ndim*nsim )
{}


ndfes::PathSims::PathSims( XMLIO::xml_node node )
  : mNumDim(0),
    mNumSim(0)
{
  if ( not node )
    {
      return;
    };

  for ( XMLIO::xml_node pnode = node.child("pt");
	pnode; pnode = pnode.next_sibling("pt") )
    {
      ++mNumSim;
      if ( mNumSim == 1 )
	{
	  mNumDim = 0;
	  for ( XMLIO::xml_node qnode = pnode.child("rc");
		qnode; qnode = qnode.next_sibling("rc") )
	    {
	      ++mNumDim;
	    }
	}
    }

  mRCs.assign( mNumSim*mNumDim, 0. );
  mFCs.assign( mNumSim*mNumDim, 0. );

  
  for ( XMLIO::xml_node pnode = node.child("pt");
	pnode; pnode = pnode.next_sibling("pt") )
    {
      int idx = 0;
      XMLIO::ExtractAttribute(pnode,"idx",idx);
      
      for ( XMLIO::xml_node qnode = pnode.child("rc");
	    qnode; qnode = qnode.next_sibling("rc") )
	{
	  int dim = 0;
	  XMLIO::ExtractAttribute(qnode,"dim",dim);
	  XMLIO::ExtractValue(qnode,mRCs[dim+idx*mNumDim]);
	}

      for ( XMLIO::xml_node qnode = pnode.child("fc");
	    qnode; qnode = qnode.next_sibling("fc") )
	{
	  int dim = 0;
	  XMLIO::ExtractAttribute(qnode,"dim",dim);
	  XMLIO::ExtractValue(qnode,mFCs[dim+idx*mNumDim]);
	}
    }
  
}



void ndfes::PathSims::WriteXml( std::ostream & cout, int const ilev ) const
{
#define IND1 std::setw(2*(ilev+0)) << ""
#define IND2 std::setw(2*(ilev+1)) << ""
#define IND3 std::setw(2*(ilev+2)) << ""

  if ( mNumSim == 0 )
    {
      return;
    };
  
  cout << IND1 << "<sims>\n";
  for ( std::size_t i=0; i<mNumSim; ++i )
    {
      cout << IND2 << "<pt idx=\"" << i << "\">\n";
      for ( std::size_t dim=0; dim<mNumDim; ++dim )
	{
	  cout << IND3
	       << "<rc dim=\"" << dim << "\">"
	       << std::scientific << std::setw(21) << std::setprecision(12)
	       << mRCs[dim+i*mNumDim]
	       << "</rc>"
	       << "\n";
	}
      for ( std::size_t dim=0; dim<mNumDim; ++dim )
	{
	  cout << IND3
	       << "<fc dim=\"" << dim << "\">"
	       << std::scientific << std::setw(21) << std::setprecision(12)
	       << mFCs[dim+i*mNumDim]
	       << "</fc>"
	       << "\n";
	}
      cout << IND2 << "</pt>\n";
    }
  cout << IND1 << "</sims>\n";

#undef IND1
#undef IND2
#undef IND3
}


///////////////////////////////////////////////////////////////////


ndfes::PathIter::PathIter( XMLIO::xml_node node )
  : mPath( node.child("path") ),
    mSims( node.child("sims") )
{
}


void ndfes::PathIter::WriteXml( std::ostream & cout, int const ilev ) const
{
#define IND1 std::setw(2*(ilev+0)) << ""
  
  cout << IND1 << "<iter>\n";

  mPath.WriteXml( cout, ilev + 1 );
  mSims.WriteXml( cout, ilev + 1 );
  
  cout << IND1 << "</iter>\n";

#undef IND1
}


///////////////////////////////////////////////////////////////////


ndfes::PathOpt::PathOpt( std::string fname )
{
  XMLIO::xml_document xml;
  XMLIO::LoadXml( fname, xml );

  XMLIO::xml_node node( XMLIO::ExtractNode(xml,"opt") );

  for ( XMLIO::xml_node inode = node.child("iter");
	inode; inode = inode.next_sibling("iter") )
    {
      mIters.push_back( ndfes::PathIter( inode ) );
    }
}


void ndfes::PathOpt::WriteXml( std::ostream & cout, int const ilev ) const
{
#define IND1 std::setw(2*(ilev+0)) << ""
  
  cout << IND1 << "<opt>\n";

  for ( std::size_t i=0; i<mIters.size(); ++i )
    {
      mIters[i].WriteXml( cout, ilev + 1 );
    };
  
  cout << IND1 << "</opt>\n";

#undef IND1
}


void ndfes::PathOpt::WriteXml( std::string fname ) const
{
  std::ofstream cout;
  cout.open(fname.c_str());
  WriteXml(cout,0);
}





///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


std::vector<double> ndfes::PathIter::PredictCentroids
( std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts ) const
{
  std::vector<double> centroids;

  if ( popts.acc )
    {
      
      std::size_t npts = mPath.mNumPts;
      if ( popts.npathpts > 1 )
	{
	  npts = popts.npathpts;
	}
      std::size_t ndim = mPath.mNumDim;
      std::vector<double> rcs( mPath.GetUniformPts( npts ) );
      std::vector<double> fcs( rcs.size(), 0 );
      centroids.assign( rcs.size(), 0 );
      for ( std::size_t i=0; i<npts; ++i )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      fcs[dim+i*ndim] = popts.deffc[dim];
	    }
	};


      std::vector<double> wrcs(rcs);
      ndfes::WrapPath(fes->mDimInfo,wrcs);
      
      ndfes::BiasedMins( ndim, npts, wrcs.data(), fcs.data(),
			 centroids.data(), popts.acc_oobk,
			 fes, popts.acc_maxit, 1.e-13,
			 popts.minbounds, popts.maxbounds );

      ndfes::UnwrapCentroids(fes->mDimInfo,centroids,rcs);
    }

  return centroids;
}





ndfes::PathIter ndfes::PathIter::Next
( std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts ) const
{
  std::vector<double> rcs( PredictCentroids(fes,popts) );
  std::size_t ndim = mPath.mNumDim;
  std::size_t npts = rcs.size() / ndim;
  if ( popts.fix0 )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  rcs[dim] = mPath.mPts[dim];
	}
    }
  if ( popts.fix1 )
    {
      std::size_t o = (npts-1)*ndim;
      std::size_t p = (mPath.mNumPts-1)*ndim;
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  rcs[dim+o] = mPath.mPts[dim+p];
	}
    }


  // for ( std::size_t i=0; i<npts; ++i )
  //   {
  //     std::printf("%5lu",i);
  //     for ( std::size_t dim=0; dim<ndim; ++dim )
  // 	{
  // 	  std::printf(" %12.6f",rcs[dim+i*ndim]);
  // 	}
  //     std::printf("\n");
  //   }
  
  ndfes::PathSpl pspl( popts.stype, popts.smooth, ndim, npts, rcs.data() );
  return ndfes::PathIter( pspl );
}



void ndfes::PathIter::PredictUniformSims
( std::size_t const onsim,
  ndfes::PathOptions const & popts )
{
  std::size_t const ndim = mPath.mNumDim;
  std::vector<double> new_rcs( mPath.GetUniformPts( onsim ) );
  std::vector<double> new_fcs( onsim*ndim, 0 );
  for ( std::size_t i=0; i<onsim; ++i )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  new_fcs[dim+i*ndim] = popts.deffc[dim];
	}
    }
  mSims = ndfes::PathSims(ndim,onsim,new_rcs.data(),new_fcs.data());
}



void ndfes::PathIter::PredictUniformSims
( std::size_t const onsim,
  std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts )
{
#define FMTF std::fixed << std::setw(8) << std::setprecision(2)

  std::size_t const ndim = mPath.mNumDim;
  std::vector<double> new_rcs( mPath.GetUniformPts( onsim ) );
  std::vector<double> new_fcs(ndim*onsim,0);
  for ( std::size_t i=0; i<onsim; ++i )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  new_fcs[dim+i*ndim] = popts.deffc[dim];
	}
    }

  double const * minfc = popts.minfc.data();
  double const * maxfc = popts.maxfc.data();
  double tdisp = popts.tdispfc;
  std::vector<double> dws( GetUniformPathDeltas(ndim,onsim,new_rcs.data()) );

  std::vector<double> optfc(ndim*onsim,0);
  std::vector<double> optw(onsim,0);

  std::vector<double> percs(onsim,0);
  std::vector<double> cosas(onsim,0);
  
  std::vector< std::vector<double> > ws(onsim);
  std::vector< std::vector<double> > ps(onsim);
  std::vector< std::vector<double> > coss(onsim);
  
  std::vector<bool> done( onsim, false );
  std::vector<double> w( onsim, 0. );      
  std::vector<double> wp( onsim, 0. );      
  std::vector<double> fc;
  
  fc = GetScaledForceConstants
    (ndim,onsim,minfc,maxfc,
     new_fcs.data(),w.data());
  CalcDisplacementPercentages
    ( ndim,onsim, percs.data(), cosas.data(), done,
      new_rcs.data(), fc.data(), dws.data(),
      fes, popts );


  for ( std::size_t i=0; i<onsim; ++i )
    {
      ws[i].push_back( w[i] );
      ps[i].push_back( percs[i] );
      coss[i].push_back( cosas[i] );
    }
  
  std::fill( w.data(), w.data() + onsim, 1. );
  
  fc = GetScaledForceConstants
    (ndim,onsim,minfc,maxfc,
     new_fcs.data(),w.data());
  CalcDisplacementPercentages
    ( ndim, onsim, percs.data(), cosas.data(), done,
      new_rcs.data(), fc.data(), dws.data(),
      fes, popts );
  
  for ( std::size_t i=0; i<onsim; ++i )
    {
      // if w=1 yielded a larger percentage, then dont keep it
      // because that doesn't make sense
      if ( percs[i] < ps[i][0] )
	{
	  ws[i].push_back( w[i] );
	  ps[i].push_back( percs[i] );
	  coss[i].push_back( cosas[i] );
	};
    }
  
  std::fill( w.data(), w.data() + onsim, -1. );
      
  fc = GetScaledForceConstants
    (ndim,onsim,minfc,maxfc,
     new_fcs.data(),w.data());
  CalcDisplacementPercentages
    ( ndim,onsim,percs.data(), cosas.data(), done,
      new_rcs.data(), fc.data(), dws.data(),
      fes, popts );



  for ( std::size_t i=0; i<onsim; ++i )
    {
      // if w=0 yielded a smaller percentage, then dont keep it
      // because that doesn't make sense
      if ( percs[i] > ps[i][0] )
	{
	  ws[i].push_back( w[i] );
	  ps[i].push_back( percs[i] );
	  coss[i].push_back( cosas[i] );
	};
    }
  
  std::vector< ndfes::LinearInterp > interps( onsim );
  for ( std::size_t i=0; i<onsim; ++i )
    {
      interps[i].reset( ps[i].size(), ps[i].data(), ws[i].data() );
    }


  for ( std::size_t it=0; it<10; ++it )
    {
      for ( std::size_t i=0; i<onsim; ++i )
	{
	  w[i] = interps[i].GetValue( tdisp );
	  wp[i] = std::min(1.,std::max(-1.,w[i]));
	}
	  
      fc = GetScaledForceConstants
	(ndim,onsim,minfc,maxfc,
	 new_fcs.data(),wp.data());
	  
      CalcDisplacementPercentages
	( ndim,onsim,percs.data(), cosas.data(), done,
	  new_rcs.data(), fc.data(), dws.data(),
	  fes, popts );

      
      std::cout << "FC optimization iteration " << it+1 << "\n";
      for ( std::size_t i=0; i<onsim; ++i )
	{
	  std::string donestr = "F";
	  if ( done[i] )
	    {
	      donestr="T";
	    };
	  if ( not done[i] )
	    {
	      std::cout << std::setw(3) << i+1 << " " << donestr << " w:"
			<< std::fixed << std::setw(6) << std::setprecision(3)
			<< wp[i]
			<< " p:"
			<< std::fixed << std::setw(6) << std::setprecision(3)
			<< percs[i] 
			<< " fc:";
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  std::cout << FMTF << fc[dim+i*ndim];
		}
	      std::cout << "\n";
	    };
	}
      
      bool all_done = true;
      for ( std::size_t i=0; i<onsim; ++i )
	{
	  if ( not done[i] )
	    {
	      all_done = false;
	      optw[i] = wp[i];
	      
	      if ( interps[i].push_back( percs[i], wp[i] ) )
		{ // we found something new
		  
		  ws[i].push_back( wp[i] );
		  ps[i].push_back( percs[i] );
		  coss[i].push_back( cosas[i] );
		  
		  // is it close?
		  if ( std::abs( percs[i] - tdisp ) < 0.001 )
		    {
		      if ( popts.scalefc )
			{
			  optw[i] = wp[i] * std::abs(cosas[i]);
			}
		      else
			{
			  optw[i] = wp[i];
			}
		      done[i] = true;
		    }
		}
	      else 
		{ // we already saw this wp
		  done[i] = true;
		  
		  std::size_t w0idx = 0;
		  std::size_t idx = 0;
		  double mindp = 1.e+8;
		  double mindw = 1.e+8;
		  for ( std::size_t kk=0; kk<ws[i].size(); ++kk )
		    {
		      double p = ps[i][kk];
		      double w = ws[i][kk];
		      double dp = std::abs(p-tdisp);
		      if ( dp < mindp )
			{
			  idx = kk;
			  mindp = dp;
			}
		      if ( std::abs(w) < mindw )
			{
			  w0idx = kk;
			  mindw = std::abs(w);
			};
		    }
		  
		  double myw = ws[i][idx];
		  double myp = ps[i][idx];
		  
		  if ( myp < tdisp and
		       ps[i][w0idx] < myp and
		       myw > 0 )
		    {
		      wp[i] = 0;
		      percs[i] = ps[i][w0idx];
		      cosas[i] = coss[i][w0idx];
		    }
		  
		  if ( popts.scalefc )
		    {
		      optw[i] = wp[i] * std::abs(cosas[i]);
		    }
		  else
		    {
		      optw[i] = wp[i];
		    }
		}
	    };
	}
      
      if ( all_done )
	{	
	  break;
	};
      
    }
  
      
  fc = GetScaledForceConstants
    (ndim,onsim,minfc,maxfc,
     new_fcs.data(),optw.data());
  
  for ( std::size_t i=0; i<onsim; ++i )
    {
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  new_fcs[dim+i*ndim] = fc[dim+i*ndim];
	}
    }
    
  {
    std::cout << "Replacing end-point FCs with "
	      << "nearest neighbor values\n";
    
    std::size_t m1 = onsim-1;
    std::size_t m2 = onsim-2;
    for ( std::size_t dim=0; dim<ndim; ++dim )
      {
	new_fcs[dim+0*ndim] = new_fcs[dim+1*ndim];
	new_fcs[dim+m1*ndim] = new_fcs[dim+m2*ndim];
      }
    
    
    std::cout << "\nPredicted force constants from "
	      << "uniform simulations\n";
    for ( std::size_t i=0; i<onsim; ++i )
      {
	std::cout << std::setw(3) << i+1;
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    std::cout << " " << FMTF << new_fcs[dim+i*ndim];
	  }
	std::cout << "\n";
      }
    std::cout << "\n";
    
  }

  
  if ( popts.smoothfc )
    {
      std::vector<double> sfcs = 
	ndfes::GetIterWinAvg(ndim,onsim,new_fcs.data(),3,11);
      for ( std::size_t i=0; i<onsim; ++i )
	{
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      new_fcs[dim+i*ndim] =
		std::min(maxfc[dim],
			 std::max(minfc[dim],
				  sfcs[dim+i*ndim]));
	    }
	}
      
      std::cout << "\nPredicted force constants from "
		<< "uniform simulations after smoothing\n";
      for ( std::size_t i=0; i<onsim; ++i )
	{
	  std::cout << std::setw(3) << i+1;
	  for ( std::size_t dim=0; dim<ndim; ++dim )
	    {
	      std::cout << " " << FMTF << new_fcs[dim+i*ndim];
	    }
	  std::cout << "\n";
	}
      std::cout << "\n";
      
    }


#undef FMTF


  mSims = ndfes::PathSims( ndim, onsim, new_rcs.data(), new_fcs.data() );
  
}
    



void ndfes::PathIter::PredictUniformCentroids
( std::size_t const onsim,
  std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts,
  std::vector<double> & new_ts )
{
  new_ts.assign( onsim, 0. );
  
  std::size_t const ndim = mPath.mNumDim;
  PredictUniformSims( onsim, fes, popts );

  std::vector<double> new_rcs( mSims.mRCs );
  std::vector<double> new_fcs( mSims.mFCs );

  std::vector<double> cpts
    ( ndfes::PredictCentroids
      ( ndim, onsim, new_rcs, new_fcs, fes, popts ) );
  
  int maxit=100;
  int nseg=1000;
  bool isakima = false;
  if ( mPath.mSType > 0 )
    {
      isakima = true;
    }
  ndfes::PCurve cspl( ndim, onsim, cpts.data(), isakima, maxit, nseg );

  std::vector<double> uts(onsim,0);
  {
    double nm1 = onsim-1;
    for ( std::size_t i=0; i<onsim; ++i )
      {
	uts[i] = i/nm1;
      }
  }
  ndfes::LinearInterp ct2ut( onsim, cspl.mT.data(), uts.data() );

  std::vector<double> ucts(onsim,0);
  for ( std::size_t i=0; i<onsim; ++i )
    {
      ucts[i] = ct2ut.GetValue( uts[i] );
    }
  std::sort( ucts.begin(), ucts.end() );


  std::vector<ndfes::LinearInterp> fcspl(ndim);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      std::vector<double> tfc(onsim,0);
      for ( std::size_t i=0; i<onsim; ++i )
	{
	  tfc[i] = new_fcs[dim+i*ndim];
	}
      fcspl[dim].reset( onsim, uts.data(), tfc.data() );
    };

  for ( std::size_t i=0; i<onsim; ++i )
    {
      double t = ucts[i];
      new_ts[i] = t;
      mPath.GetValue( t, new_rcs.data() + i*ndim );
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  new_fcs[dim+i*ndim] = fcspl[dim].GetValue(t);
	}
    }
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      new_fcs[dim+0*ndim]         = new_fcs[dim+1*ndim];
      new_fcs[dim+(onsim-1)*ndim] = new_fcs[dim+(onsim-2)*ndim];
    }


  {
    std::vector<double> cpts
      ( ndfes::PredictCentroids
	( ndim, onsim, new_rcs, new_fcs, fes, popts ) );
    
    ndfes::PCurve fspl( ndim, onsim, cpts.data(), isakima, maxit, nseg );
    

    std::cout << "\nPredicted centroid distributions.\n"
	      << "Col 1: uniform discretization\n"
	      << "Col 2: predicted centroids from uniform window centers\n"
	      << "Col 3: predicted centroids from uniform centroids\n";

#define FMTF std::fixed << std::setw(9) << std::setprecision(5)
    
    for ( std::size_t i=0; i<onsim; ++i )
      {
	std::cout << FMTF << uts[i]
		  << " " << FMTF << cspl.mT[i]
		  << " " << FMTF << fspl.mT[i]
		  << "\n";
      }
    
#undef FMTF

    
    std::cout << "\n\nPredicted force constants from uniform centroids\n";

#define FMTF std::fixed << std::setw(9) << std::setprecision(2)
    
    for ( std::size_t i=0; i<onsim; ++i )
      {
	std::cout << std::setw(3) << i+1;
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    std::cout << FMTF << new_fcs[dim+i*ndim];
	  }
	std::cout << "\n";
      }
    std::cout << "\n";
  }
  
#undef FMTF

  mSims.mRCs = new_rcs;
  mSims.mFCs = new_fcs;
  
}







// void ndfes::PathIter::PredictMinOccSims
// ( std::size_t const onsim,
//   std::shared_ptr<ndfes::FES> fes,
//   ndfes::PathOptions const & popts )
// {
//   new_ts.assign( onsim, 0. );
  
//   std::size_t const ndim = mPath.mNumDim;
//   PredictUniformSims( onsim, popts );

//   std::vector<double> new_rcs( mSims.mRCs );
//   std::vector<double> new_fcs( mSims.mFCs );

//   std::size_t nscan=10;
//   double tfrac = 2./5.;
//   double dt = 1./(onsim-1.);
//   double tprev=0;
//   for ( std::size_t i=0; i<onsim; ++i )
//     {
//       double t = i*dt;
//       double tlo = std::max(0.,t-tfrac*dt);
//       double thi = std::min(1.,t+tfrac*dt);
//       if ( i > 0 )
// 	{
// 	  tlo = std::max(tlo,tprev+0.25*dt);
// 	  thi = std::max(thi,tlo);
// 	}
//       double tmin=tlo;
//       std::size_t nmin=100000000;
//       for ( std::size_t iscan=0; iscan<nscan; ++iscan )
// 	{
// 	  double myt = 0;
// 	  std::vector<double> c( ndim, 0 );
// 	  std::size_t nocc = 0;
	  
// 	  myt = t + iscan*((tlo-t)/(nscan-1.));
// 	  mPath.GetValue( myt, c.data() );
// 	  nocc = fes->BinOccAtPt( c.data() );
// 	  if ( nocc < nmin )
// 	    {
// 	      nmin = nocc;
// 	      tmin = myt;
// 	    }

// 	  myt = t + iscan*((thi-t)/(nscan-1.));
// 	  mPath.GetValue( myt, c.data() );
// 	  nocc = fes->BinOccAtPt( c.data() );
// 	  if ( nocc < nmin )
// 	    {
// 	      nmin = nocc;
// 	      tmin = myt;
// 	    }	  
// 	}
//       tprev=tmin;

//       mPath.GetValue( tmin, new_rcs.data() + i*ndim );
//     }
  
//   mSims.mRCs = new_rcs;
//   mSims.mFCs = new_fcs;
  
// }






std::vector<ndfes::AreaSummary> ndfes::PathIter::PredictGridBasedSims
( std::size_t const nsim,
  std::shared_ptr<ndfes::FES> fes,
  ndfes::PathOptions const & popts,
  int const nlayers,
  std::vector<double> & reservedpts )
{
  PredictUniformSims( nsim, popts );

  std::size_t const ndim = mPath.mNumDim;
  return ndfes::SamplingConvOpt
    ( ndim, nsim,
      mSims.mRCs.data(), mSims.mFCs.data(),
      mPath.mSpl, fes, popts,
      nlayers, reservedpts );
}



///////////////////////////////////////////////////////////////////


bool ndfes::PathOpt::CheckRepeatedPath( std::shared_ptr<ndfes::FES> fes ) const
{
  bool same = false;

  if ( mIters.size() > 2 )
    {
      ndfes::PathSpl const * last = &(mIters.back().mPath);
      std::size_t nsim = last->mNumPts;
      std::size_t ndim = last->mNumDim;
      
      bool allocc = true;
      for ( std::size_t isim=0; isim<nsim; ++isim )
	{
	  bool isocc = fes->PointIsOcc( last->mPts.data() + isim*ndim );
	  if ( not isocc )
	    {
	      allocc = false;
	    }
	}
      
      double TOL = 1.e-7;
      if ( not allocc )
	{
	  TOL = 1.e-3;
	}

      
      int istart = (int)mIters.size() - 7;
      if ( istart < 0 )
	{
	  istart = 0;
	}

      std::vector<double> disps;
      
      for ( std::size_t i=istart; i<mIters.size()-1; ++i )
	{
	  ndfes::PathSpl const * cur = &(mIters[i].mPath);
	  if ( cur->mNumPts != nsim )
	    {
	      continue;
	    }
	  double maxdisp = -1.e+30;
	  for ( std::size_t isim=0; isim<nsim; ++isim )
	    {
	      double disp = 0;
	      for ( std::size_t dim=0; dim<ndim; ++dim )
		{
		  double dx = last->mPts[dim+isim*ndim] - cur->mPts[dim+isim*ndim];
		  disp += dx*dx;
		};
	      disp = std::sqrt(disp);
	      if ( disp > maxdisp )
		{
		  maxdisp = disp;
		};
	    };
	  disps.push_back( maxdisp );

	}

      std::size_t ndisp = disps.size();
      std::size_t imin = 0;
      double maxdisp = 1.e+30;
      for ( std::size_t i=0; i<ndisp; ++i )
	{
	  if ( disps[i] < maxdisp )
	    {
	      maxdisp = disps[i];
	      imin = i;
	    };
	}
      if ( imin < ndisp - 1 and maxdisp < TOL )
	{
	  same = true;
	}
    }
  return same;
}
