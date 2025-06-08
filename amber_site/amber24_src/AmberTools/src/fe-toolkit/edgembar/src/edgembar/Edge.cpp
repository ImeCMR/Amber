#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>    

#include "Edge.hpp"
#include "Opt.hpp"

edgembar::Edge::Edge
( std::string const name,
  edgembar::Env const & complex,
  edgembar::Env const & solvated )
  : mName(name),
    mNumParam(0),
    mComplex(complex),
    mSolvated(solvated)
{
  mComplex.SetLinkedList(this);
  mSolvated.SetLinkedList(this);


  for ( int ienv = 0; ienv < 2; ++ienv )
    {
      edgembar::Env * env = NULL;
      double coefsign = 1.;
      if ( ienv == 0 )
	{
	  env = &mComplex;
	  coefsign = 1;
	}
      else
	{
	  env = &mSolvated;
	  coefsign = -1;
	}
      
      for ( edgembar::Stage
	      * s = env->GetStageBegin(),
	      * send = env->GetStageEnd();
	    s != send; ++s )
	{
	  int ntrials = s->GetNumTrials();
	  for ( edgembar::Trial
		  * t = s->GetTrialBegin(),
		  * tend = s->GetTrialEnd();
		t != tend; ++t )
	    {
	      int nstates = t->GetNumStates();
	      
	      // offset to the lambda=0 state
	      int o = mNumParam;
	      t->SetParamOffset( o );
	      
	      // lambda=0 state
	      mFreeEneIdxs.push_back( o );
	      mFreeEneCoefs.push_back( -coefsign / ntrials );
	      mFreeShifts.push_back( 0. );
	      
	      // lambda=1 state
	      mFreeEneIdxs.push_back( nstates + o - 1 );
	      mFreeEneCoefs.push_back(  coefsign / ntrials );
	      mFreeShifts.push_back( t->GetShift() );
	      
	      mNumParam += nstates;
	    }
	}
    }
  RemoveConstraints();
}


void edgembar::Edge::StoreFreeAvgEnes()
{
  mFreeAvgEnes.resize(0);
  mAvgEnes.resize(0);
  
  for ( int ienv = 0; ienv < 2; ++ienv )
    {
      edgembar::Env * env = NULL;
      if ( ienv == 0 )
       	{
       	  env = &mComplex;
       	}
      else
       	{
       	  env = &mSolvated;
       	}
      
      for ( edgembar::Stage
	      * s = env->GetStageBegin(),
	      * send = env->GetStageEnd();
	    s != send; ++s )
	{
	  for ( edgembar::Trial
		  * t = s->GetTrialBegin(),
		  * tend = s->GetTrialEnd();
		t != tend; ++t )
	    {
	      int nstates = t->GetNumStates();
	      
	      // lambda=0 state
	      mFreeAvgEnes.push_back( t->GetAvgEne( 0 ) );
	      
	      // lambda=1 state
	      mFreeAvgEnes.push_back( t->GetAvgEne( nstates-1 ) );

	      for ( int i=0; i<nstates; ++i )
		{
		  mAvgEnes.push_back( t->GetAvgEne(i) - t->GetAvgEne(0) );
		  if ( i < nstates - 1 )
		    {
		      mShifts.push_back( 0. );
		    }
		  else
		    {
		      mShifts.push_back( t->GetShift() );
		    }
		}
	      
	    }
	}
    }
}




std::vector<edgembar::Sim const *> edgembar::Edge::GetConstSims() const
{
  std::vector<edgembar::Sim const *> p( mComplex.GetConstSims() );
  std::vector<edgembar::Sim const *> q( mSolvated.GetConstSims() );
  p.insert(p.end(),q.begin(),q.end());
  return p;
}

std::vector<edgembar::Sim *> edgembar::Edge::GetSims() 
{
  std::vector<edgembar::Sim *> p( mComplex.GetSims() );
  std::vector<edgembar::Sim *> q( mSolvated.GetSims() );
  p.insert(p.end(),q.begin(),q.end());
  return p;
}




std::vector<edgembar::Trial const *> edgembar::Edge::GetConstTrials() const
{
  std::vector<edgembar::Trial const *> p( mComplex.GetConstTrials() );
  std::vector<edgembar::Trial const *> q( mSolvated.GetConstTrials() );
  p.insert(p.end(),q.begin(),q.end());
  return p;
}

std::vector<edgembar::Trial *> edgembar::Edge::GetTrials() 
{
  std::vector<edgembar::Trial *> p( mComplex.GetTrials() );
  std::vector<edgembar::Trial *> q( mSolvated.GetTrials() );
  p.insert(p.end(),q.begin(),q.end());
  return p;
}




std::shared_ptr<edgembar::Edge>
edgembar::Edge::ExtractCopy() const
{
  std::shared_ptr<edgembar::Edge> edge( new edgembar::Edge(*this) );
  edge->mComplex.SetLinkedList(edge.get());
  edge->mSolvated.SetLinkedList(edge.get());
  return edge;
}
  
std::shared_ptr<edgembar::Edge>
edgembar::Edge::ExtractRange( double const start, double const stop ) const
{
  std::shared_ptr<edgembar::Edge> edge( new edgembar::Edge(*this) );
  edge->mComplex.SetLinkedList(edge.get());
  edge->mSolvated.SetLinkedList(edge.get());
  edge->mComplex.ExtractRangeInPlace(start,stop);
  edge->mSolvated.ExtractRangeInPlace(start,stop);
  return edge;
}
  
std::shared_ptr<edgembar::Edge>
edgembar::Edge::ExtractProd() const
{
  std::shared_ptr<edgembar::Edge> edge( new edgembar::Edge(*this) );
  edge->mComplex.SetLinkedList(edge.get());
  edge->mSolvated.SetLinkedList(edge.get());
  edge->mComplex.ExtractProdInPlace();
  edge->mSolvated.ExtractProdInPlace();
  return edge;
}
    
std::shared_ptr<edgembar::Edge>
edgembar::Edge::ExtractBootstrap() const
{
  std::shared_ptr<edgembar::Edge> edge( new edgembar::Edge(*this) );
  edge->mComplex.SetLinkedList(edge.get());
  edge->mSolvated.SetLinkedList(edge.get());
  edge->mComplex.ExtractBootstrapInPlace();
  edge->mSolvated.ExtractBootstrapInPlace();
  return edge;
}


void edgembar::Edge::RemoveConstraints()
{
  mConstraints.reset( mNumParam );
}

void edgembar::Edge::ApplyMinimalConstraints()
{
  mConstraints.reset( mNumParam );
  if ( (int)(mFreeEneIdxs.size()) % 2 == 1 )
    {
      std::cerr << "Error in edgembar::Constraint::ApplyMinimalConstraints() "
		<< "Expected mFreeEneIdxs to appear in pairs\n";
      std::exit(EXIT_FAILURE);
    }
  int const ncon = mFreeEneIdxs.size() / 2;
  std::vector<int> idxs(1,0);
  std::vector<double> const coefs(1,1.);
  for ( int icon=0; icon<ncon; ++icon )
    {
      idxs[0] = mFreeEneIdxs[0+icon*2];
      mConstraints.push_back(idxs,coefs,0.);
    }
  mConstraints.finalize();
}


void edgembar::Edge::ApplyConstraint( double const cval )
{
  mConstraints.reset( mNumParam );
  if ( (int)(mFreeEneIdxs.size()) % 2 == 1 )
    {
      std::cerr << "Error in edgembar::Constraint::ApplyConstraint() "
		<< "Expected mFreeEneIdxs to appear in pairs\n";
      std::exit(EXIT_FAILURE);
    }
  int const ncon = mFreeEneIdxs.size() / 2;
  std::vector<int> idxs(1,0);
  std::vector<double> const coefs(1,1.);
  for ( int icon=0; icon<ncon; ++icon )
    {
      idxs[0] = mFreeEneIdxs[0+icon*2];
      mConstraints.push_back(idxs,coefs,0.);
    }
  double coff = 0.;
  for ( std::size_t i=0,n=mFreeEneIdxs.size(); i<n; ++i )
    {
      coff += mFreeEneCoefs[i] * ( mFreeAvgEnes[i] + mFreeShifts[i] );
    }
  mConstraints.push_back(mFreeEneIdxs,mFreeEneCoefs,cval-coff);
  mConstraints.finalize();
}


std::vector<double> edgembar::Edge::GetFreeFromParams
( std::vector<double> const & p ) const
{
  return mConstraints.GetFreeFromParams(p);
}
    
std::vector<double> edgembar::Edge::GetParamsFromFree
( std::vector<double> const & q ) const
{
  return mConstraints.GetParamsFromFree(q);
}

// std::vector<double> edgembar::Edge::GetStateFreeEneFromParams
// ( std::vector<double> const & p ) const
// {
//   std::vector<double> fe(p);
//   for ( std::size_t i=0, n=mAvgEnes.size(); i<n; ++i )
//     {
//       fe[i] += mAvgEnes[i];
//     }
//   return fe;
// }



double edgembar::Edge::Optimize
( std::vector<double> const & pguess,
  std::vector<double> & popt,
  std::vector<double> & fopt,
  double const tol,
  int const verbosity ) const
{
  double chisq = edgembar::Optimize(this,pguess,popt,tol,verbosity);
  fopt = popt;
  for ( std::size_t i=0, n=popt.size(); i<n; ++i )
    {
      fopt[i] += mAvgEnes[i] + mShifts[i];
    }
  return chisq;
}


std::vector<double> edgembar::Edge::MakeExpAvgGuess() const
{
  int np = mNumParam;
  std::vector<double> p(np,0.);
  std::vector< edgembar::Trial const * > trials(GetConstTrials());
  int ntrial = trials.size();
  //bool biggap = false;
  for ( int i=0; i<ntrial; ++i )
    {
      int o = trials[i]->GetParamOffset();
      trials[i]->MakeExpAvgGuess( p.data() + o );
      //if ( mygap )
      //{
      //  biggap = true;
      //}
    }
  
  //if ( biggap )
  //{
  //  for ( int i=0; i<np; ++i )
  //	{
  //	  p[i] = 0;
  //	};
      //std::printf("note: reverting to --null guess\n");
  //}
  return p;
}


double edgembar::Edge::GetFreeEnergy( std::vector<double> const & fopt ) const
{
  if ( (int)fopt.size() != mNumParam )
    {
      std::cerr << "Error in edgembar::Edge::GetFreeEnergy input parameter "
		<< "array does not match expected number of parameters "
		<< fopt.size() << " " << mNumParam << "\n";
      std::exit( EXIT_FAILURE );
    };
  
  int nidx = mFreeEneIdxs.size();
  double ene = 0.;
  for ( int i=0; i<nidx; ++i )
    {
      ene += mFreeEneCoefs[i] * fopt[ mFreeEneIdxs[i] ];
      //ene += mFreeEneCoefs[i] * ( p[ mFreeEneIdxs[i] ] + mFreeAvgEnes[i] );
    }
  return ene;
}


std::string edgembar::Edge::GetPython
( int argc,
  char * argv[],
  double const beta,
  std::vector<edgembar::CalcIdx> const & calcs,
  std::vector<edgembar::StateValues> const & vals,
  double const ptol ) const
{
  typedef std::vector<edgembar::CalcIdx>::const_iterator calcit;

#define FMTF std::fixed << std::setprecision(8)
#define FMTE std::scientific << std::setprecision(12)

  std::vector<double> fe;
  for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
    {
      if ( ! p->IsCon and ! p->IsFwd and ! p->IsRev and ! p->IsBootstrap )
	{
	  int ibase = p->BaseIdx;
	  fe = vals[ibase].GetResult();
	  break;
	}
    }

  
  
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION undefined
#endif
  
  
  std::stringstream ss;
  ss << "edgembar.Edge("
     << "\"" << mName << "\","
     << mComplex.GetPython()
     << ","
     << mSolvated.GetPython()
     << ",";

  ss << "results=edgembar.Results(";
  
  
  ss << "version=\"" << PACKAGE_VERSION << "\",";
  ss << "command=\"\"\"";
  for ( int i=0; i<argc; ++i )
    {
      ss << argv[i];
      if ( i < argc-1 )
	{
	  ss << " ";
	}
    }
  ss << "\"\"\",";

   std::chrono::time_point<std::chrono::system_clock> cur_time
     = std::chrono::system_clock::now();
  std::time_t cur_time_value = std::chrono::system_clock::to_time_t(cur_time);
  std::string time_string = std::ctime(&cur_time_value);
  time_string.pop_back();
  ss << "date=\"" << time_string << "\",";

  ss << "ptol=" << FMTF << ptol << ",";
  
  ss << "prod=np.array([";
  {
    for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
      {
	if ( ! p->IsCon and ! p->IsFwd and ! p->IsRev and ! p->IsBootstrap )
	  {
	    int ibase = p->BaseIdx;
	    std::vector<double> v( vals[ibase].GetResult() );
	    std::vector<double> e( vals[ibase].GetStderr() );
	    for ( int i=0; i<mNumParam; ++i )
	      {
		ss << "[" << FMTF << v[i]/beta << "," << e[i]/beta << "]";
		if ( i < mNumParam-1 )
		  {
		    ss << ",";
		  };
	      }
	    break;
	  }
      }
  }
  ss << "])";

  bool hasCon=false;
  for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
    {
      if ( p->IsCon and ! p->IsBootstrap )
	{
	  hasCon = true;
	  break;
	}
    }
  if ( hasCon )
    {
      ss << ",con=[";

      bool isfirst=true;
      for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
	{
	  if ( p->IsCon and ! p->IsBootstrap )
	    {
	      if ( ! isfirst )
		{ ss << ","; }
	      else
		{ isfirst=false; }
	      
	      int ibase = p->BaseIdx;
	      ss << "edgembar.ConstraintData("
		 << FMTE << p->ConVal/beta << ","
		 << FMTE << p->Chisq/beta << ",";
	      
	      std::vector<double> v( vals[ibase].GetResult() );
	      std::vector<double> e( vals[ibase].GetStderr() );
	      ss << "np.array([";
	      for ( int i=0; i<mNumParam; ++i )
		{
		  ss << "[" << FMTF << v[i]/beta << "," << e[i]/beta << "]";
		  if ( i < mNumParam-1 )
		    {
		      ss << ",";
		    };
		}
	      ss << "])";

	      ss << ")";
	    }
	  
	}
      
      ss << "]";
    }



  int basecalcall = -1;
  bool hasFwd=false;
  for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
    {
      if ( p->IsFwd and ! p->IsBootstrap )
	{
	  hasFwd = true;
	  break;
	}
    }
  if ( hasFwd )
    {
      ss << ",fwd=[";

      bool isfirst=true;
      for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
	{
	  if ( p->IsFwd and ! p->IsBootstrap )
	    {
	      if ( ! isfirst )
		{ ss << ","; }
	      else
		{ isfirst=false; }

	      if ( p->RangeHi > 0.995 )
		{
		  basecalcall = p->BaseIdx;
		}
	      
	      int ibase = p->BaseIdx;
	      ss << "edgembar.TimeSeriesData("
		 << FMTF << p->RangeHi << ",";
	      
	      std::vector<double> v( vals[ibase].GetResult() );
	      std::vector<double> e( vals[ibase].GetStderr() );
	      ss << "np.array([";
	      for ( int i=0; i<mNumParam; ++i )
		{
		  ss << "[" << FMTF << v[i]/beta << "," << e[i]/beta << "]";
		  if ( i < mNumParam-1 )
		    {
		      ss << ",";
		    };
		}
	      ss << "])";

	      ss << ")";
	    }
	  
	}
      
      ss << "]";
    }



  bool hasRev=false;
  for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
    {
      if ( p->IsRev and ! p->IsBootstrap )
	{
	  hasRev = true;
	  break;
	}
    }
  if ( hasRev )
    {
      ss << ",rev=[";

      bool isfirst=true;
      for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
	{
	  if ( p->IsRev and ! p->IsBootstrap )
	    {
	      if ( ! isfirst )
		{ ss << ","; }
	      else
		{ isfirst=false; }
	      
	      int ibase = p->BaseIdx;
	      ss << "edgembar.TimeSeriesData("
		 << FMTF << 1. - p->RangeLo << ",";
	      
	      std::vector<double> v( vals[ibase].GetResult() );
	      std::vector<double> e( vals[ibase].GetStderr() );
	      ss << "np.array([";
	      for ( int i=0; i<mNumParam; ++i )
		{
		  ss << "[" << FMTF << v[i]/beta << "," << e[i]/beta << "]";
		  if ( i < mNumParam-1 )
		    {
		      ss << ",";
		    };
		}
	      ss << "])";

	      ss << ")";
	    }
	  
	}
      if ( basecalcall >= 0 )
	{
	  int ibase = basecalcall;
	  ss << ",edgembar.TimeSeriesData("
	     << FMTF << 1. << ",";
	      
	  std::vector<double> v( vals[ibase].GetResult() );
	  std::vector<double> e( vals[ibase].GetStderr() );
	  ss << "np.array([";
	  for ( int i=0; i<mNumParam; ++i )
	    {
	      ss << "[" << FMTF << v[i]/beta << "," << e[i]/beta << "]";
	      if ( i < mNumParam-1 )
		{
		  ss << ",";
		};
	    }
	  ss << "])";
	  ss << ")";
	}
      
      ss << "]";
    }

  

  bool hasHalf1=false;
  for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
    {
      if ( p->IsHalf1 and ! p->IsBootstrap )
	{
	  hasHalf1 = true;
	  break;
	}
    }

  if ( hasHalf1 )
    {
      ss << ",fhalf=[";

      bool isfirst=true;
      for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
	{
	  if ( p->IsHalf1 and ! p->IsBootstrap )
	    {
	      if ( ! isfirst )
		{ ss << ","; }
	      else
		{ isfirst=false; }
	      
	      int ibase = p->BaseIdx;
	      ss << "edgembar.TimeSeriesData("
		 << FMTF << p->RangeLo << ",";
	      
	      std::vector<double> v( vals[ibase].GetResult() );
	      std::vector<double> e( vals[ibase].GetStderr() );
	      ss << "np.array([";
	      for ( int i=0; i<mNumParam; ++i )
		{
		  ss << "[" << FMTF << v[i]/beta << "," << e[i]/beta << "]";
		  if ( i < mNumParam-1 )
		    {
		      ss << ",";
		    };
		}
	      ss << "])";

	      ss << ")";
	    }
	  
	}
      ss << "]";
    }
  
  bool hasHalf2=false;
  for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
    {
      if ( p->IsHalf2 and ! p->IsBootstrap )
	{
	  hasHalf2 = true;
	  break;
	}
    }

  if ( hasHalf2 )
    {
      ss << ",lhalf=[";

      bool isfirst=true;
      for ( calcit p=calcs.begin(),pend=calcs.end(); p!=pend; ++p )
	{
	  if ( p->IsHalf2 and ! p->IsBootstrap )
	    {
	      if ( ! isfirst )
		{ ss << ","; }
	      else
		{ isfirst=false; }
	      
	      int ibase = p->BaseIdx;
	      ss << "edgembar.TimeSeriesData("
		 << FMTF << p->RangeLo << ",";
	      
	      std::vector<double> v( vals[ibase].GetResult() );
	      std::vector<double> e( vals[ibase].GetStderr() );
	      ss << "np.array([";
	      for ( int i=0; i<mNumParam; ++i )
		{
		  ss << "[" << FMTF << v[i]/beta << "," << e[i]/beta << "]";
		  if ( i < mNumParam-1 )
		    {
		      ss << ",";
		    };
		}
	      ss << "])";

	      ss << ")";
	    }
	  
	}
      ss << "]";
    }
  
  
  ss << "))";
  return ss.str();
}


void edgembar::Edge::UseUniformWts()
{
  std::vector<edgembar::Trial *> p( GetTrials() );
  for ( std::size_t i=0, n=p.size(); i<n; ++i )
    {
      p[i]->UseUniformWts();
    }
}
