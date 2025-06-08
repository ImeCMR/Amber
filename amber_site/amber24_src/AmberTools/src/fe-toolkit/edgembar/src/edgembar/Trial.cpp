#define _USE_MATH_DEFINES
#include <cmath>

#include "Trial.hpp"
#include "Edge.hpp"
#include "Env.hpp"
#include "Stage.hpp"
#include "StatsMod.hpp"

#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>





static bool file_exists (std::string const name)
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}





edgembar::Trial::Trial()
  : mReadMode(edgembar::AUTO),
    mBeta(0),
    mShift(0),
    mParamOffset(0),
    mEdge(NULL),
    mEnv(NULL),
    mStage(NULL),
    mUniformWts(false)
{

}


edgembar::Trial::Trial
( std::string const name,
  std::vector<std::string> const & statelabels,
  std::string const & datadir,
  double const beta,
  edgembar::ReadMode const readmode,
  int const autoeqmode,
  double const shift )
  : mName(name),
    mStateLabels(statelabels),
    mDataDir(datadir),
    mReadMode(readmode),
    mBeta(beta),
    mShift(shift),
    mParamOffset(0),
    mEdge(NULL),
    mEnv(NULL),
    mStage(NULL),
    mUniformWts(false)
{

  int ns = mStateLabels.size();
  std::vector<int> avail(ns*ns,0);
  for ( int t=0; t<ns; ++t )
    {
      for ( int e=0; e<ns; ++e )
	{
	  std::stringstream sname;
	  sname << mDataDir << "/" << "efep_"
		<< mStateLabels[t] << "_"
		<< mStateLabels[e] << ".dat";
	  if ( file_exists(sname.str()) )
	    {
	      avail[e+t*ns] = 1;
	    };
	}
    }

  if ( mReadMode == edgembar::AUTO )
    {
      mReadMode = edgembar::AutoDetectMode(ns,avail.data());
      if ( mReadMode == edgembar::INVALID )
	{
	  std::cerr << "The files in " << mDataDir << " are incomplete,"
		    << " so a calculation cannot be performed.\n";
	  std::exit(EXIT_FAILURE);
	};
    }

  if ( ! edgembar::CanEvalMode(ns,avail.data(),mReadMode) )
    {
      std::cerr << "The files in " << mDataDir << " are incomplete,"
		<< " so a calculation cannot be performed.\n";
      std::exit(EXIT_FAILURE);
    }
  
  std::vector< std::vector<int> > simidxs
    ( edgembar::GetSimIndexes(ns,mReadMode) );

  /*
  for ( std::size_t i=0; i<simidxs.size(); ++i )
    {
      std::cout << "simidxs[" << i << "] (" << ns << ") = ";
      for ( std::size_t k=0; k<simidxs[i].size(); ++k )
	{
	  std::cout << std::setw(3) << simidxs[i][k];
	}
      std::cout << "\n";
    }
  
  
  std::string mstr = edgembar::GetModeString(mReadMode);
  std::cout << std::setw(10) << mstr << " " << mDataDir << "\n";
  for ( int t=0; t<ns; ++t )
    {
      for ( int e=0; e<ns; ++e )
	{
	  std::printf("%3i",avail[e+t*ns]);
	}
      std::printf("\n");
    }
  std::printf("\n");

  for ( int t=0; t<ns; ++t )
    {
      std::printf("%3i : ",t);
      for ( int e=0, ne=simidxs[t].size(); e<ne; ++e )
	{
	  std::printf("%3i",simidxs[t][e]);
	}
      std::printf("\n");
    }
  std::printf("\n");
  */
  

  for ( int t=0; t<ns; ++t )
    {
      if ( (int)simidxs[t].size() > 0 )
	{
	  mSims.push_back( edgembar::Sim(t,simidxs[t],mStateLabels,
					 mDataDir,mBeta,autoeqmode) );
	}
      else
	{
	  mSims.push_back( edgembar::Sim(t) );
	}
    };
  
  
}



void edgembar::Trial::StoreAvgEnes()
{
  std::size_t const ns = mStateLabels.size();
  mAvgEnes.assign( ns, 0. );

  std::size_t istart=0;
  std::size_t istop=ns;
  bool do0 = false;
  bool do1 = false;
  switch ( mReadMode )
    {
    case edgembar::MBAR:
    case edgembar::BAR:
      {
	break;
      }
    case edgembar::MBAREXP0:
    case edgembar::BAREXP0:
      {
	istart=1;
	do0 = true;
	break;
      }
    case edgembar::MBAREXP1:
    case edgembar::BAREXP1:
      {
	istop=ns-1;
	do1 = true;
	break;
      }
    case edgembar::MBAREXP:
    case edgembar::BAREXP:
      {
	istart=1;
	istop=ns-1;
	do0 = true;
	do1 = true;
	break;
      }
    default:
      {
	std::cerr << "Invalid read mode in "
		  << "edgembar::Trial::ShiftEnesAndPrecomputeAvgs\n";
	std::exit(1);
      }
    };


  for ( std::size_t i=istart; i<istop; ++i )
    {
      double avg = 0.;
      std::size_t h = mSims[i].LocalSimIdx;
      std::size_t nh = mSims[i].EneIdxs.size();
      std::size_t nk = mSims[i].NumSamples;
      for ( std::size_t k=0; k<nk; ++k )
	{
	  avg += mSims[i].Emat[h+k*nh];
	}
      mAvgEnes[i] = avg / nk;
    }

  if ( do0 )
    {
      double avg = 0.;
      std::size_t i = 1;
      std::size_t h = mSims[i].FindLocalIdx( 0 );
      if ( h > ns-1 )
	{
	  std::cerr << "FindLocalIdx (do0) produced " << h
		    << " on state " << ns-1 << " for sim "
		    << i << std::endl;
	  std::exit(1);
	}
      std::size_t nh = mSims[i].EneIdxs.size();
      std::size_t nk = mSims[i].NumSamples;
      for ( std::size_t k=0; k<nk; ++k )
	{
	  avg += mSims[i].Emat[h+k*nh];
	}
      mAvgEnes[0] = avg / nk;
    }

  if ( do1 )
    {
      double avg = 0.;
      std::size_t i = ns-2;
      std::size_t h = mSims[i].FindLocalIdx( ns-1 );
      if ( h > ns-1 )
	{
	  std::cerr << "FindLocalIdx (do1) produced " << h
		    << " on state " << ns-1 << " for sim "
		    << i << std::endl;
	  std::exit(1);
	}
      std::size_t nh = mSims[i].EneIdxs.size();
      std::size_t nk = mSims[i].NumSamples;
      for ( std::size_t k=0; k<nk; ++k )
	{
	  avg += mSims[i].Emat[h+k*nh];
	}
      mAvgEnes[ns-1] = avg / nk;
    }
  
  
  for ( std::size_t isim=0; isim<ns; ++isim )
    {
      int nene = mSims[isim].EneIdxs.size();
      for ( int i=0; i<nene; ++i )
      	{
      	  int idx = mSims[isim].EneIdxs[i];
      	  mSims[isim].AvgEnes[i] = mAvgEnes[idx];
      	}
    }
}



void edgembar::Trial::UseUniformWts()
{
  mUniformWts=true;
}


void edgembar::Trial::SetLinkedList
( edgembar::Edge * edge,
  edgembar::Env * env,
  edgembar::Stage * stage )
{
  mEdge=edge;
  mEnv=env;
  mStage=stage;
}


edgembar::Trial edgembar::Trial::ExtractRange
( double const start,
  double const stop ) const
{
  edgembar::Trial T(*this);

  int ns = mSims.size();
  for ( int i=0; i<ns; ++i )
    {
      T.mSims[i].OnlyKeepRange(start,stop);
    }
  return T;
}


edgembar::Trial edgembar::Trial::ExtractProd() const
{
  edgembar::Trial T(*this);

  int ns = mSims.size();
  for ( int i=0; i<ns; ++i )
    {
      T.mSims[i].OnlyKeepProd();
    }
  return T;
}


edgembar::Trial edgembar::Trial::ExtractBootstrap() const
{
  edgembar::Trial T(*this);

  int ns = mSims.size();
  for ( int i=0; i<ns; ++i )
    {
      T.mSims[i].Bootstrap();
    }
  return T;
}






void edgembar::Trial::ExtractRangeInPlace
( double const start,
  double const stop )
{
  int ns = mSims.size();
  for ( int i=0; i<ns; ++i )
    {
      mSims[i].OnlyKeepRange(start,stop);
    }
}


void edgembar::Trial::ExtractProdInPlace()
{
  int ns = mSims.size();
  for ( int i=0; i<ns; ++i )
    {
      mSims[i].OnlyKeepProd();
    }
}


void edgembar::Trial::ExtractBootstrapInPlace()
{
  int ns = mSims.size();
  for ( int i=0; i<ns; ++i )
    {
      mSims[i].Bootstrap();
    }
}



std::vector<edgembar::Sim const *> edgembar::Trial::GetConstSims() const
{
  std::vector<edgembar::Sim const *> p;
  for ( int i=0, n=mSims.size(); i<n; ++i )
    {
      p.push_back( mSims.data() + i );
    }
  return p;
}

std::vector<edgembar::Sim *> edgembar::Trial::GetSims()
{
  std::vector<edgembar::Sim *> p;
  for ( int i=0, n=mSims.size(); i<n; ++i )
    {
      p.push_back( mSims.data() + i );
    }
  return p;
}



void edgembar::Trial::WriteDebugInfo( std::ostream & cout ) const
{
  cout << "Trial"
       << " Name: " << std::left << std::setw(14) << GetName()
       << " Stage: " << std::left << std::setw(12) << GetStage()->GetName()
       << " Env: " << std::left << std::setw(12) << GetEnv()->GetName()
       << " Edge: " << std::left << GetEdge()->GetName()
       << "\n"
       << std::setw(6) << ""
       << "Mode: " << std::left << std::setw(10)
       << edgembar::GetModeString(mReadMode)
       << " StateLabels:";
  for ( int i=0, n=mStateLabels.size(); i<n; ++i )
    {
      cout << " " << mStateLabels[i];
    }
  cout << "\n";
}


void edgembar::Trial::TestObjective() const
{
  double const DEL=5.e-4;
  
  int const Nstates = mStateLabels.size();
  
  std::vector<double> p(Nstates,0);
  std::vector<double> g(Nstates,0);
  for ( int i=0; i<Nstates; ++i )
    {
      p[i] = -i*2.5; //i+1;
    }

  double c = CptObjective(p.data(),g.data());
  std::printf("%20.10e\n",c);

  std::vector<double> tmp(Nstates,0);
  for ( int i=0; i<Nstates; ++i )
    {
      p[i] += DEL;
      double chi = CptObjective(p.data(),tmp.data());
      p[i] -= 2*DEL;
      double clo = CptObjective(p.data(),tmp.data());
      p[i] += DEL;

      double d = (chi-clo)/(2*DEL);
      std::printf("%3i ana: %13.4e num: %13.4e diff: %13.4e  %20.10e %20.10e\n",
		  i,g[i],d,g[i]-d,chi,clo);
    }
}



namespace edgembar
{
  double MBARObj( double const wt,
		  int const Nstates,
		  edgembar::Sim const * begin,
		  double const * p,
		  double * g );

  double BARObj( double const wt,
		 int const Nstates,
		 edgembar::Sim const * begin,
		 double const * p,
		 double * g );
  
  double FwdEXPObj( double const wt,
		    edgembar::Sim const * sim,
		    double const * psim,
		    double * gsim );
  
  double RevEXPObj( double const wt,
		    edgembar::Sim const * sim,
		    double const * pene,
		    double * gene );
}


double edgembar::MBARObj
( double const wt,
  int const Nstates,
  edgembar::Sim const * begin,
  double const * p,
  double * g )
{
  double chisq = 0;
  int Nsamples = 0;
  for ( int i=0; i<Nstates; ++i )
    {
      Nsamples += begin[i].NumSamples;
    }
  double const ooN = 1./Nsamples;
	
  std::vector<double> b( Nstates, 0. );
  //double maxz = -1.e+30;
  for ( int i=0; i<Nstates; ++i )
    {
      //std::printf("NumSamples[%i] = %i %12.3e\n",i,begin[i].NumSamples,p[i]);
      double NoverN = 1.;
      if ( begin[i].NumSamples > 0 )
	{
	  NoverN = ooN * begin[i].NumSamples;
	  b[i] = - p[i] - std::log( NoverN );
	}
      else
	{
	  b[i] = - p[i];
	}
      chisq += wt * NoverN * b[i];
      g[i] -= wt * NoverN;
      //maxz = std::max(maxz,-b[i]);
    }

  //std::printf("pre %20.10e %20.10e\n",chisq,ooN);
  
  for ( int i=0; i<Nstates; ++i )
    {
      //b[i] = std::exp(-b[i]-maxz);
      b[i] = std::exp(-b[i]);
    }

  

  for ( int i=0; i<Nstates; ++i )
    {
      int const joff = begin[i].LocalSimIdx - i;
      int const nene = begin[i].EneIdxs.size();
      double const * M = begin[i].Mdat.data();
      double const * Z = begin[i].MaxZ.data();
      for ( int k=0, nk=begin[i].NumSamples; k<nk; ++k )
	{
	  double den = 0.;
	  for ( int j=0; j<Nstates; ++j )
	    {
	      den += M[j+joff+k*nene] * b[j];
	    };
	  //chisq += ooN * ( std::log(den) + maxz + Z[k] );
	  //chisq += ooN * ( std::log(den) + maxz );
	  //chisq += wt * ooN * ( std::log(den) );
	  chisq += wt * ooN * ( std::log(den) + Z[k] );
	  den = ooN/den;
	  for ( int j=0; j<Nstates; ++j )
	    {
	      g[j] += wt * den * M[j+joff+k*nene] * b[j];
	    };
	};
    };

  //std::printf("post %20.10e\n",chisq);

  
  return chisq;
}



double edgembar::BARObj
( double const wt,
  int const Nstates,
  edgembar::Sim const * begin,
  double const * p,
  double * g )
{
  double chisq = 0;
  for ( int i=0; i<Nstates-1; ++i )
    {
      chisq += edgembar::MBARObj(wt,2,begin+i,p+i,g+i);
    }
  return chisq;
}


double edgembar::FwdEXPObj
( double const wt,
  edgembar::Sim const * sim,
  double const * psim,
  double * gsim )
{
  /*
    VERSION 1
    dF       = F1 - F0
    dU(x)    = U1(x) - U0(x)
    dF       = - ln (1/N0) \sum_{k=1}^{N0} exp[-dU(x)]
    exp[-dF] = (1/N0) \sum_{k=1}^{N0} exp[-dU(x)]
    1        = (1/N0) \sum_{k=1}^{N0} exp[dF-dU(x)]
    0        = (1/N0) \sum_{k=1}^{N0} exp[dF-dU(x)] - 1   [Eq 1]
    Eq. 1 is satisfied by minimizing:
    dX/d(dF) = 0
    X = (1/N0) \sum_{k=1}^{N0} exp[dF-dU(x)] - dF
    And the individual derivatives come from chain rule:
    dX/dF1 = dX/d(dF) * d(dF)/dF1 = + dX/d(dF)
    dX/dF0 = dX/d(dF) * d(dF)/dF0 = - dX/d(dF)

    VERSION 2
    Alternatively, let
    A        = (1/N0) \sum_{k=1}^{N0} exp[dF-dU(x)]
    then
    1        = A
    ln(1)    = ln( A )
    0        = ln( A )
    and this is satisified by minimizing
    X = A ( ln(A) - 1 )
    whose derivatives are
    dX/dF1 =   ln( A )
    dX/dF0 = - ln( A )
    The above expressions can be scaled by an arbitrary constant
    (because 0*c=0).  I choose the constant to by 1/N0 so that
    the objective function looks like the fast mbar function
    (which is (1/N) * ln( something )
    0 = ln(A)/N0
    X = A ( ln(A) - 1 )/N0
    dXdF1 =  ln(A)/N0
    dXdF0 = -ln(A)/N0

    VERSION 3
    No, no. That's not right. Let
    B = (1/N0) \sum_{k=1}^{N0} exp(-dU(x))
    Then Zwanzig's equation reads
    dF = - ln(B)
    Rearranging, (negate and exponentiate both sides)
    exp[-dF] = B
    Now multiply both sides by exp[dF]
    1 = exp[dF]*B
    Now take the natural log of both sides
    0 = ln( exp[dF]*B )
    If we let A = exp[dF]*B, then this could also be written as
    0 = ln(A)
    Now integrate wrt dF
    X = (1/2) [ ln(exp[dF]*B) ]^2
    but ln(uv) = ln(u) + ln(v), so
    X = (1/2) [ ln(exp[dF]) + ln(B) ]^2
    X = (1/2) [ dF + ln(B) ]^2
    Alternatively,
    X = (1/2) [ ln(A) ]^2
    Rhe derivatives are
    dX/dF1 =  ln(A)
    dX/dF0 = -ln(A)

    Note that X is literally a quadratic equation
    X = (1/2) ( dF - c )^2
    centered at c = -ln(B)

   */
  double const dF = psim[1]-psim[0];
  double const expdF = std::exp(dF);
  int const N0 = sim->NumSamples;
  int const i = sim->LocalSimIdx; // initial state
  int const j = i+1; // final state

  //std::printf("FwdExpOBJ N0 %i %i %i \n",N0,i,j);
  
  //std::cout << i << " " << j << " " << sim->EneIdxs[i] << " " << sim->EneIdxs[j] << "\n";
  
  int const nene = sim->EneIdxs.size();
  double const * M = sim->Mdat.data();
  double esum = 0.;
  for ( int k=0; k<N0; ++k )
    {
      esum += expdF * M[j+k*nene] / M[i+k*nene];
      //esum += M[j+k*nene] / M[i+k*nene];
      
      // std::printf("%5i %14.5e %14.5e %14.3f %14.5e %14.5e %14.5e\n",
      // 		  k,
      // 		  sim->Emat[j+k*nene],
      // 		  sim->Emat[i+k*nene],
      // 		  dF-(sim->Emat[j+k*nene]-sim->Emat[i+k*nene]),
      // 		  std::exp(dF - (sim->Emat[j+k*nene]-sim->Emat[i+k*nene])),
      // 		  expdF * M[j+k*nene] / M[i+k*nene],
      // 		  esum );
      
    };
  esum /= N0;

  // v1
  //double X = wt*(esum-dF);
  //double dXddF = wt*(esum - 1);

  // v2
  //double logsum = std::log(esum+1.e-300);
  //double X = wt*esum * (logsum-1.); // /N0;
  //double dXddF = wt*logsum; // /N0;

  // v3
  double logsum = std::log(esum+1.e-300);
  double X = wt*logsum*logsum/2.;
  double dXddF = wt*logsum;
  //double X = wt*( 0.5*dF*dF + dF*logsum );
  //double dXddF = wt*(dF+logsum);

  //std::printf("esum,logsum %20.10e %20.10e %20.10e\n",esum,logsum,X);
  
  gsim[1] += dXddF; // final state
  gsim[0] -= dXddF; // initial state
  return X;
}


double edgembar::RevEXPObj
( double const wt,
  edgembar::Sim const * sim,
  double const * pene,
  double * gene )
{
  /*
    In this notation, F1 is the final state and F0 is the initial state
    However, in the variables of this routine, the final state is 0
    and the initial state is 1 -- hence why I call it "reversed"
    exponential averaging.

    VERSION 1
    dF       = F1 - F0   (final - initial)
    dU(x)    = U1(x) - U0(x)
    dF       = - ln (1/N0) \sum_{k=1}^{N0} exp[-dU(x)]
    exp[-dF] = (1/N0) \sum_{k=1}^{N0} exp[-dU(x)]
    1        = (1/N0) \sum_{k=1}^{N0} exp[dF-dU(x)]
    0        = (1/N0) \sum_{k=1}^{N0} exp[dF-dU(x)] - 1   [Eq 1]
    Eq. 1 is satisfied by minimizing:
    dX/d(dF) = 0
    X = (1/N0) \sum_{k=1}^{N0} exp[dF-dU(x)] - dF
    And the individual derivatives come from chain rule:
    dX/dF1 = dX/d(dF) * d(dF)/dF1 = + dX/d(dF)
    dX/dF0 = dX/d(dF) * d(dF)/dF0 = - dX/d(dF)

    VERSION 2
    Alternatively, let
    A        = (1/N0) \sum_{k=1}^{N0} exp[dF-dU(x)]
    then
    1        = A
    ln(1)    = ln( A )
    0        = ln( A )
    and this is satisified by minimizing
    X = A ( ln(A) - 1 )
    whose derivatives are
    dX/dF1 =   ln( A )
    dX/dF0 = - ln( A )
    The above expressions can be scaled by an arbitrary constant
    (because 0*c=0).  I choose the constant to by 1/N0 so that
    the objective function looks like the fast mbar function
    (which is (1/N) * ln( something )
    0 = ln(A)/N0
    X = A ( ln(A) - 1 )/N0
    dXdF1 =  ln(A)/N0
    dXdF0 = -ln(A)/N0

    VERSION 3
    No, no. That's not right. Let
    B = (1/N0) \sum_{k=1}^{N0} exp(-dU(x))
    Then Zwanzig's equation reads
    dF = - ln(B)
    Rearranging, (negate and exponentiate both sides)
    exp[-dF] = B
    Now multiply both sides by exp[dF]
    1 = exp[dF]B
    Now take the natural log of both sides
    0 = ln( exp[dF] B )
    If we let A = exp[dF] B, then this could also be written as
    0 = ln(A)
    Now integrate wrt dF
    X = (1/2) [ ln(exp[dF]B) ]^2
    but ln(uv) = ln(u) + ln(v), so
    X = (1/2) [ ln(exp[dF]) + ln(B) ]^2
    X = (1/2) [ dF + ln(B) ]^2
    Alternatively,
    X = (1/2) [ ln(A) ]^2
    Rhe derivatives are
    dX/dF1 =  ln(A)
    dX/dF0 = -ln(A)

    Note that X is literally a quadratic equation
    X = (1/2) ( dF - c )^2
    centered at c = -ln(B)
   */
  
  double const dF = pene[0]-pene[1]; // final - initial
  double const expdF = std::exp(dF);
  int const N1 = sim->NumSamples;
  
  int const i = sim->LocalSimIdx; // i is the initial state
  int const j = i-1; // j is the final state
  
  int const nene = sim->EneIdxs.size();
  double const * M = sim->Mdat.data();
  double esum = 0.;
  for ( int k=0; k<N1; ++k )
    {
      esum += expdF * M[j+k*nene] / M[i+k*nene];
    };
  esum /= N1;

  // v1
  //double X = wt*(esum - dF);
  //double dXddF = wt*(esum - 1);

  // v2
  //double logsum = std::log(esum+1.e-300);
  //double X = wt*esum * (logsum-1.); // /N1;
  //double dXddF = wt*logsum; // /N1;

  // v3
  double logsum = std::log(esum+1.e-300);
  double X = wt*logsum*logsum/2.;
  double dXddF = wt*logsum; 
  
  gene[0] += dXddF; // final state
  gene[1] -= dXddF; // initial state
  return X;
}



/*
  fi = - ln (1/N) \sum_k exp(-ui) / [ sum_h (Nh/N) exp(fh-uh) ]
  exp(-fi) = (1/N) \sum_k exp(-ui) / [ sum_h (Nh/N) exp(fh-uh) ]
  0 = 1 - (1/N) \sum_k exp(fi-ui) / [ sum_h (Nh/N) exp(fh-uh) ]
  0 = Ni - \sum_k [ (Ni/N) exp(fi-ui) ] / [ sum_h (Nh/N) exp(fh-uh) ]
  0 = - \sum_k { [ exp( ln(Ni/N) + fi-ui) ] / [ sum_h exp( ln(Nh/N) + fh-uh) ] + (Ni/N) }

 */



double edgembar::Trial::CptObjective( double const * p, double * g ) const
{
  int const Nstates = mStateLabels.size();

  int ntrials = GetStage()->GetNumTrials();
  double const wt = 1. / ntrials;

  //double const mbarwt = wt * (1./Nstates);
  //double const barwt = mbarwt * (2./Nstates);
  //double const expwt = mbarwt * (1./Nstates);
  //double const barwt = mbarwt;
  //double const expwt = mbarwt;

  
  double chisq = 0;
  std::fill(g,g+Nstates,0);

  
  switch ( mReadMode )
    {
    case edgembar::MBAR:
      {
	double mbarwt = wt;
	if ( mUniformWts )
	  {
	    mbarwt = wt;
	  }
	//double const mbarwt = wt;
	chisq = edgembar::MBARObj(mbarwt,Nstates,mSims.data(),p,g);
	break;
      }
    case edgembar::MBAREXP0:
      {
	double mbarwt = wt;
	if ( Nstates > 2 )
	  {
	    if ( ! mUniformWts )
	      {
		mbarwt = wt;
	      }
	    chisq = edgembar::MBARObj(mbarwt,Nstates-1,mSims.data()+1,p+1,g+1);
	  };
	//double const mbarwt = wt;
	double const expwt = wt;
	chisq += edgembar::RevEXPObj(expwt,mSims.data()+1,p,g);
	break;
      }
    case edgembar::MBAREXP1:
      {
	double mbarwt = wt;
	if ( Nstates > 2 )
	  {
	    if ( ! mUniformWts )
	      {
		mbarwt = wt;
	      }
	    chisq = edgembar::MBARObj(mbarwt,Nstates-1,mSims.data(),p,g);
	  };
	//double const mbarwt = wt;
	double const expwt = wt;
	int const o = Nstates-2;
	chisq += edgembar::FwdEXPObj(expwt,mSims.data()+o,p+o,g+o);
	break;
      }
    case edgembar::MBAREXP:
      {
	double mbarwt = wt;
	if ( Nstates > 3 )
	  {
	    if ( ! mUniformWts )
	      {
		mbarwt = wt;
	      }
	    chisq = edgembar::MBARObj(mbarwt,Nstates-2,mSims.data()+1,p+1,g+1);
	  }
	//double const mbarwt = wt;
	double const expwt = wt;
	int const o = Nstates-2;
	chisq += edgembar::RevEXPObj(expwt,mSims.data()+1,p,g);
	chisq += edgembar::FwdEXPObj(expwt,mSims.data()+o,p+o,g+o);
	break;
      }
    case edgembar::BAR:
      {
	double barwt = wt ; // * (2./(Nstates));
	if ( mUniformWts )
	  {
	    barwt = wt;
	  }
	//double const barwt = wt;
	chisq = edgembar::BARObj(barwt,Nstates,mSims.data(),p,g);
	break;
      }
    case edgembar::BAREXP0:
      {
	double barwt = wt;
	if ( Nstates > 2 )
	  {
	    if ( ! mUniformWts )
	      {
		barwt = wt ; // * (2./(Nstates-1.));
	      }
	    chisq = edgembar::BARObj(barwt,Nstates-1,mSims.data()+1,p+1,g+1);
	  }
	//double const barwt = wt;
	double const expwt = wt;
	chisq += edgembar::RevEXPObj(expwt,mSims.data()+1,p,g);
	break;
      }
    case edgembar::BAREXP1:
      {
	double barwt = wt;
	if ( Nstates > 2 )
	  {
	    if ( ! mUniformWts )
	      {
		barwt = wt ; // * (2./(Nstates-1.));
	      }
	    chisq = edgembar::BARObj(barwt,Nstates-1,mSims.data(),p,g);
	  }
	//double const barwt = wt;
	double const expwt = wt;
	int const o = Nstates-2;
	chisq += edgembar::FwdEXPObj(expwt,mSims.data()+o,p+o,g+o);
	break;
      }
    case edgembar::BAREXP:
      {
	double barwt = wt;
	if ( Nstates > 3 )
	  {
	    if ( ! mUniformWts )
	      {
		barwt = wt ; // * (2./(Nstates-2.));
	      }
	    chisq = edgembar::BARObj(barwt,Nstates-2,mSims.data()+1,p+1,g+1);
	  }
	//double const barwt = wt;
	double const expwt = wt;
	int const o = Nstates-2;
	chisq += edgembar::RevEXPObj(expwt,mSims.data()+1,p,g);
	chisq += edgembar::FwdEXPObj(expwt,mSims.data()+o,p+o,g+o);
	break;
      }
    default:
      {
	std::cerr << "Programming error in edgembar::Trial::CptObjective :"
		  << "unimplemented readmode " << GetModeString(mReadMode)
		  << std::endl;
	std::exit(EXIT_FAILURE);
      }
    }
  return chisq;
}



namespace edgembar
{
  double FwdExp( edgembar::Sim const * sim );
  double RevExp( edgembar::Sim const * sim );
  double GapRE( edgembar::Sim const * sim, edgembar::Sim const * ene );
  double GapOverlap( edgembar::Sim const * A, edgembar::Sim const * B );
}

double edgembar::FwdExp( edgembar::Sim const * sim )
{
  int const N0 = sim->NumSamples;
  int const i = sim->LocalSimIdx; // initial state
  int const j = i+1; // final state
  
  int const nene = sim->EneIdxs.size();
  double const * E = sim->Emat.data();
  double maxz = -1.e+100;
  for ( int k=0; k<N0; ++k )
    {
      double z = E[i+k*nene]-E[j+k*nene];
      maxz = std::max(maxz,z);
      //esum += M[j+k*nene] / M[i+k*nene];
    };
  double esum = 0.;
  for ( int k=0; k<N0; ++k )
    {
      double z = (E[i+k*nene]-E[j+k*nene]);
      esum += std::exp(z-maxz);
      //esum += M[j+k*nene] / M[i+k*nene];
    };
  esum /= N0;
  return - ( std::log( esum + 1.e-300 ) + maxz );
}

double edgembar::RevExp( edgembar::Sim const * sim )
{
  int const N0 = sim->NumSamples;
  int const i = sim->LocalSimIdx; // initial state
  int const j = i-1; // final state
  
  int const nene = sim->EneIdxs.size();
  // double const * M = sim->Mdat.data();
  // double esum = 0.;
  // for ( int k=0; k<N0; ++k )
  //   {
  //     esum += M[j+k*nene] / M[i+k*nene];
  //   };
  double const * E = sim->Emat.data();
  double maxz = -1.e+100;
  for ( int k=0; k<N0; ++k )
    {
      double z = E[i+k*nene]-E[j+k*nene];
      maxz = std::max(maxz,z);
    };
  double esum = 0.;
  for ( int k=0; k<N0; ++k )
    {
      double z = (E[i+k*nene]-E[j+k*nene]);
      esum += std::exp(z-maxz);
    };
  esum /= N0;
  return std::log( esum + 1.e-300 ) + maxz;
}


double edgembar::GapRE( edgembar::Sim const * sim, edgembar::Sim const * ene )
{
  int const N0 = sim->NumSamples;
  int const i = sim->LocalSimIdx; // initial state
  int const j = sim->FindLocalIdx(ene); // final state
  
  //std::printf("ene_in_sim %4i %4i  N0 %4i\n",i,j,N0);

  if ( j < 0 or N0 <= 0 )
    {
      return -1;
    };

  std::vector<double> ws(N0,0.);
  
  int const nene = sim->EneIdxs.size();
  /*
  double const * M = sim->Mdat.data();
  double esum = 0.;
  for ( int k=0; k<N0; ++k )
    {
      ws[k] = M[j+k*nene] / M[i+k*nene];
      esum += ws[k];
    };
  for ( int k=0; k<N0; ++k )
    {
      ws[k] /= esum;
    };
  double num = 0.;
  double den = std::log( (double)N0 );
  for ( int k=0; k<N0; ++k )
    {
      num += ws[k] * std::log( ws[k] + 1.e-300 );
    }
  */

  double const * E = sim->Emat.data();
  double maxz = -1.e+100;
  for ( int k=0; k<N0; ++k )
    {
      double z = E[i+k*nene] - E[j+k*nene];
      maxz = std::max(maxz,z);
    }
  double esum = 0.;
  for ( int k=0; k<N0; ++k )
    {
      double z = E[i+k*nene] - E[j+k*nene];
      ws[k] = std::exp(z-maxz);
      esum += ws[k];
    }
    for ( int k=0; k<N0; ++k )
    {
      ws[k] /= esum;
    };
  double num = 0.;
  double den = std::log( (double)N0 );
  for ( int k=0; k<N0; ++k )
    {
      num += ws[k] * std::log( ws[k] + 1.e-300 );
    }

  return - num/den;
}



double edgembar::GapOverlap( edgembar::Sim const * A, edgembar::Sim const * B )
{
  int const B_in_A = A->FindLocalIdx( B );
  int const A_in_B = B->FindLocalIdx( A );

  //std::printf("B_in_A %4i  A_in_B %4i\n",B_in_A,A_in_B);
  
  if ( B_in_A < 0 or
       A_in_B < 0 or
       A->NumSamples == 0 or
       B->NumSamples == 0 )
    {
      return -1;
    }
  else if ( A == B )
    {
      return 1;
    }
  
  
  double mA = 0.;
  double vA = 0.;
  {
    int const N0 = A->NumSamples;
    int const i = A->LocalSimIdx; // initial state
    int const j = B_in_A; // final state
    
    int const nene = A->EneIdxs.size();
    double const * M = A->Emat.data();
    std::vector<double> dU( N0,0. );
    for ( int k=0; k<N0; ++k )
      {
	dU[k] = M[j+k*nene] - M[i+k*nene];
	mA += dU[k];
      };
    mA /= N0;
    for ( int k=0; k<N0; ++k )
      {
	double x = dU[k]-mA;
	vA += x*x;
      }
    vA /= N0;
  }

  double mB = 0.;
  double vB = 0.;
  {
    int const N0 = B->NumSamples;
    int const i = B->LocalSimIdx; // initial state
    int const j = A_in_B; // final state
    
    int const nene = B->EneIdxs.size();
    double const * M = B->Emat.data();
    std::vector<double> dU( N0,0. );
    for ( int k=0; k<N0; ++k )
      {
	dU[k] = - (M[j+k*nene] - M[i+k*nene]);
	mB += dU[k];
      };
    mB /= N0;
    for ( int k=0; k<N0; ++k )
      {
	double x = dU[k]-mB;
	vB += x*x;
      }
    vB /= N0;
  }

  double zA = 1.;
  if ( vA > 1.e-10 )
    {
      zA = 0.5 / vA;
    }
  double zB = 1.;
  if ( vB > 1.e-10 )
    {
      zB = 0.5 / vB;
    };
  double zab = zA*zB/(zA+zB);
  double sab = std::sqrt( zab/M_PI ) * std::exp( -zab * std::pow(mA-mB,2) );
  double zaa = zA*zA/(zA+zA);
  double saa = std::sqrt( zaa/M_PI );
  double zbb = zB*zB/(zB+zB);
  double sbb = std::sqrt( zbb/M_PI );

  //std::printf("gau %13.4e %13.4e   %13.4e %13.4e\n",vA,mA,vB,mB);
  return sab/std::max(saa,sbb);
}


std::vector<double> edgembar::Trial::GetOverlaps() const
{
  int const n = mSims.size();
  std::vector<double> S(n*n,0.);
  for ( int i=0; i<n; ++i )
    {
      for ( int j=0; j<n; ++j )
	{
	  S[j+i*n] = edgembar::GapOverlap( mSims.data() + i, mSims.data() + j );
	};
    };
  return S;
}


std::vector<double> edgembar::Trial::GetEntropies() const
{
  int n = mSims.size();
  std::vector<double> RE(n*n,0.);
  for ( int i=0; i<n; ++i )
    {
      for ( int j=0; j<n; ++j )
	{
	  //std::printf("Get GapRE %4i %4i\n",i,j);
	  RE[j+i*n] = edgembar::GapRE( mSims.data() + i, mSims.data() + j );
	};
    };
  return RE;
}


void edgembar::Trial::MakeExpAvgGuess( double * p ) const
{
  int const Nstates = mStateLabels.size();
  std::fill(p,p+Nstates,0.);


  
  switch ( mReadMode )
    {
    case edgembar::MBAR:
    case edgembar::BAR:
      {
	for ( int i=0; i<Nstates-1; ++i )
	  {
	    double fwd = edgembar::FwdExp( mSims.data()+i );
	    double rev = edgembar::RevExp( mSims.data()+i+1 );
	    p[i+1] = p[i] + 0.5*(fwd+rev);
	  }
	break;
      }
    case edgembar::MBAREXP0:
    case edgembar::BAREXP0:
      {
	p[1] = edgembar::RevExp( mSims.data()+1 );
	for ( int i=1; i<Nstates-1; ++i )
	  {
	    double fwd = edgembar::FwdExp( mSims.data()+i );
	    double rev = edgembar::RevExp( mSims.data()+i+1 );
	    p[i+1] = p[i] + 0.5*(fwd+rev);
	  }
	break;
      }
    case edgembar::MBAREXP1:
    case edgembar::BAREXP1:
      {
	for ( int i=0; i<Nstates-2; ++i )
	  {
	    double fwd = edgembar::FwdExp( mSims.data()+i );
	    double rev = edgembar::RevExp( mSims.data()+i+1 );
	    p[i+1] = p[i] + 0.5*(fwd+rev);
	  }
	p[Nstates-1] = p[Nstates-2] + edgembar::FwdExp( mSims.data()+Nstates-2 );
	break;
      }
    case edgembar::MBAREXP:
    case edgembar::BAREXP:
      {
	p[1] = edgembar::RevExp( mSims.data()+1 );
	for ( int i=1; i<Nstates-2; ++i )
	  {
	    double fwd = edgembar::FwdExp( mSims.data()+i );
	    double rev = edgembar::RevExp( mSims.data()+i+1 );
	    p[i+1] = p[i] + 0.5*(fwd+rev);
	  }
	p[Nstates-1] = p[Nstates-2] + edgembar::FwdExp( mSims.data()+Nstates-2 );
	break;
      }
    default:
      {
	std::cerr << "Programming error in edgembar::Trial::MakeExpAvgGuess :"
		  << "unimplemented readmode " << GetModeString(mReadMode)
		  << std::endl;
	std::exit(EXIT_FAILURE);
      }
    }

  for ( int i=0; i<Nstates; ++i )
    {
      if ( not std::isfinite(p[i]) )
	{
	  if ( i == 0 )
	    {
	      p[i] = 0;
	    }
	  else
	    {
	      p[i] = p[i-1];
	    }
	}
      if ( i < Nstates - 1 )
	{
	  double d = std::abs(p[i+1]-p[i]) / mBeta;
	  if ( d > 100 )
	    {
	      p[i+1] = p[i];
	    }
	};
    }
  
}


/*
void edgembar::Trial::GetResults_Production
( double const beta,
  std::vector<edgembar::CalcIdx> const & calcs,
  std::vector<edgembar::StateValues> const & vals,
  std::vector<double> & res ) const
{
  int const ns = mStateLabels.size();
  int const o = mParamOffset;
  res.assign( ns * 2, 0. );
  
  for ( int i=0, n=calcs.size(); i<n; ++i )
    {
      if ( ! calcs[i].IsCon and
	   ! calcs[i].IsFwd and
	   ! calcs[i].IsRev and
	   ! calcs[i].IsBootstrap )
	{
	  
	  int ibase = calcs[i].BaseIdx;
	  
	  
	  std::vector<double> v( vals[ibase].GetResult() );
	  std::vector<double> e( vals[ibase].GetStderr() );

	  for ( std::size_t k=0; k<ns; ++k )
	    {
	      res[0+k*2] = v[k+o];
	      res[1+k*2] = e[k+o];
	    }
	  break;
	}
    }

}
*/


std::string edgembar::Trial::GetPython() const
{
  int ns = mStateLabels.size();
  
  std::stringstream ss;
  ss << "edgembar.Trial("
     << "\"" << mName << "\","
     << "\"" << mDataDir << "\",[";
  for ( int i=0; i<ns; ++i )
    {
      ss << "\"" << mStateLabels[i] << "\"";
      if ( i < ns-1 )
	{
	  ss << ",";
	}
    }
  ss << "],"
     << "\"" << edgembar::GetModeString(mReadMode) << "\","
     << "offset=" << mParamOffset;

  std::vector<double> S( GetOverlaps() );
  std::vector<double> RE( GetEntropies() );

  bool hasdvdl = true;
  for ( int i=0; i<ns; ++i )
    {
      if ( (int)mSims[i].DVDL.size() == 0 )
	{
	  hasdvdl = false;
	  break;
	}
    }
  std::vector<double> dvdl_avgs(ns,0);
  std::vector<double> dvdl_errs(ns,0);
  if ( hasdvdl )
    {
      for ( int i=0; i<ns; ++i )
	{
	  int n = mSims[i].DVDL.size();
	  int g = ccdl::CptStatIneff(n,mSims[i].DVDL.data());
	  std::vector<double> tvec;
	  for ( int j=0; j<n; j += g )
	    {
	      tvec.push_back( mSims[i].DVDL[j] + mShift );
	      //std::cerr << std::scientific << std::setw(13) << std::setprecision(4) << tvec.back() << "\n";
	    }
	  //std::cerr << "tvec.size() " << i << " " << tvec.size() << std::endl;
	  ccdl::CptMeanAndStderr( tvec.size(),tvec.data(),
				  dvdl_avgs[i], dvdl_errs[i] );
	  //std::cerr << "dvdl_errs " << dvdl_errs[i] << std::endl;

	  dvdl_avgs[i] = 0;
	  for ( int j=0; j<n; ++j )
	    {
	      dvdl_avgs[i] += mSims[i].DVDL[j];
	    }
	  dvdl_avgs[i] /= n;
	  dvdl_avgs[i] += mShift/mBeta;
	};
    }
  
  std::stringstream res;
  res << ",results=[";
  for ( int i=0; i<ns; ++i )
    {
      std::stringstream m;
      if ( mSims[i].OrigNumSamples > 0 )
	{
	  m << "edgembar.SimProperty("
	    << mSims[i].OrigStride << ","
	    << mSims[i].OrigNumSamples << ","
	    << mSims[i].ProdStart << ","
	    << mSims[i].ProdStride << ","
	    << mSims[i].CurStride << ","
	    << mSims[i].NumSamples << ","
	    << mSims[i].AutoEqMode << ",";
	  
	  if ( mSims[i].IsConverged )
	    {
	      m << "True,";
	    }
	  else
	    {
	      m << "False,";
	    };
	  
	  m << "np.array([";
	  for ( int k=0; k<ns; ++k )
	    {
	      if ( k > 0 )
		{
		  m << ",";
		}
	      m << std::scientific << std::setprecision(3) << S[k+i*ns];
	    }
	  m << "]),np.array([";
	  for ( int k=0; k<ns; ++k )
	    {
	      if ( k > 0 )
		{
		  m << ",";
		}
	      m << std::scientific << std::setprecision(3) << RE[k+i*ns];
	    }
	  m << "])";

	  if ( hasdvdl )
	    {
	      m << ","
		<< std::scientific << std::setprecision(10) << dvdl_avgs[i]
		<< ","
		<< std::scientific << std::setprecision(10) << dvdl_errs[i];
	    }
	  else
	    {
	      m << ",None,None";
	    }

	  m << ")";
	}
      else
	{
	  m << "None";
	}
      if ( i > 0 )
	{
	  res << "," << m.str();
	}
      else
	{
	  res << m.str();
	}
    };
  res << "]";
	  
  ss << res.str() << ")";
  
  /*
  std::printf("\nOverlap\n");
  for ( int j=0; j<ns; ++j )
    {
      for ( int i=0; i<ns; ++i )
	{
	  std::printf("%8.3f",S[i+j*ns]);
	};
      std::printf("\n");
    };
  
  std::printf("Entropy\n");
  for ( int j=0; j<ns; ++j )
    {
      for ( int i=0; i<ns; ++i )
	{
	  std::printf("%8.3f",RE[i+j*ns]);
	};
      std::printf("\n");
    };
  */
  
  return ss.str();
}
