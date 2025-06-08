#include <sstream>

#include "Stage.hpp"

edgembar::Stage::Stage()
  : mEdge(NULL),
    mEnv(NULL)
{
  
}


edgembar::Stage::Stage
( std::string const name,
  std::vector<edgembar::Trial> const & trials )
  : mName(name),
    mTrials(trials),
    mEdge(NULL),
    mEnv(NULL)
{
  
}




void edgembar::Stage::SetLinkedList
( edgembar::Edge * edge,
  edgembar::Env * env )
{
  mEdge=edge;
  mEnv=env;
  for ( int i=0, n=mTrials.size(); i<n; ++i )
    {
      mTrials[i].SetLinkedList(edge,env,this);
    }
}


void edgembar::Stage::ExtractRangeInPlace
( double const start,
  double const stop )
{
  int ns = mTrials.size();
  for ( int i=0; i<ns; ++i )
    {
      mTrials[i].ExtractRangeInPlace(start,stop);
    }
}


void edgembar::Stage::ExtractProdInPlace()
{
  int ns = mTrials.size();
  for ( int i=0; i<ns; ++i )
    {
      mTrials[i].ExtractProdInPlace();
    }
}


void edgembar::Stage::ExtractBootstrapInPlace()
{
  int ns = mTrials.size();
  for ( int i=0; i<ns; ++i )
    {
      mTrials[i].ExtractBootstrapInPlace();
    }
}




std::vector<edgembar::Sim const *> edgembar::Stage::GetConstSims() const
{
  std::vector<edgembar::Sim const *> p;
  for ( int i=0, n=mTrials.size(); i<n; ++i )
    {
      std::vector<edgembar::Sim const *> q( mTrials[i].GetConstSims() );
      p.insert(p.end(),q.begin(),q.end());
    }
  return p;
}

std::vector<edgembar::Sim *> edgembar::Stage::GetSims()
{
  std::vector<edgembar::Sim *> p;
  for ( int i=0, n=mTrials.size(); i<n; ++i )
    {
      std::vector<edgembar::Sim *> q( mTrials[i].GetSims() );
      p.insert(p.end(),q.begin(),q.end());
    }
  return p;
}

std::vector<edgembar::Trial const *> edgembar::Stage::GetConstTrials() const
{
  std::vector<edgembar::Trial const *> p;
  for ( int i=0, n=mTrials.size(); i<n; ++i )
    {
      p.push_back( mTrials.data() + i );
    }
  return p;
}

std::vector<edgembar::Trial *> edgembar::Stage::GetTrials()
{
  std::vector<edgembar::Trial *> p;
  for ( int i=0, n=mTrials.size(); i<n; ++i )
    {
      p.push_back( mTrials.data() + i );
    }
  return p;
}


std::string edgembar::Stage::GetPython() const
{
  std::stringstream ss;
  ss << "edgembar.Stage("
     << "\"" << mName << "\",[";
  for ( int i=0, n=mTrials.size(); i<n; ++i )
    {
      ss << mTrials[i].GetPython();
      if ( i < n-1 )
        {
          ss << ",";
        };
    }
  ss << "])";

  return ss.str();
}
