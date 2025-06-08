#include <sstream>

#include "Env.hpp"

edgembar::Env::Env()
  : mEdge(NULL)
{

}


edgembar::Env::Env
( std::string const name,
  std::vector<edgembar::Stage> const & stages )
  : mName(name),
    mStages(stages),
    mEdge(NULL)
{

}



void edgembar::Env::SetLinkedList
( edgembar::Edge * edge )
{
  mEdge=edge;
  for ( int i=0, n=mStages.size(); i<n; ++i )
    {
      mStages[i].SetLinkedList(edge,this);
    }
}



void edgembar::Env::ExtractRangeInPlace
( double const start,
  double const stop )
{
  int ns = mStages.size();
  for ( int i=0; i<ns; ++i )
    {
      mStages[i].ExtractRangeInPlace(start,stop);
    }
}


void edgembar::Env::ExtractProdInPlace()
{
  int ns = mStages.size();
  for ( int i=0; i<ns; ++i )
    {
      mStages[i].ExtractProdInPlace();
    }
}


void edgembar::Env::ExtractBootstrapInPlace()
{
  int ns = mStages.size();
  for ( int i=0; i<ns; ++i )
    {
      mStages[i].ExtractBootstrapInPlace();
    }
}



std::vector<edgembar::Sim const *> edgembar::Env::GetConstSims() const
{
  std::vector<edgembar::Sim const *> p;
  for ( int i=0, n=mStages.size(); i<n; ++i )
    {
      std::vector<edgembar::Sim const *> q( mStages[i].GetConstSims() );
      p.insert(p.end(),q.begin(),q.end());
    }
  return p;
}

std::vector<edgembar::Sim *> edgembar::Env::GetSims()
{
  std::vector<edgembar::Sim *> p;
  for ( int i=0, n=mStages.size(); i<n; ++i )
    {
      std::vector<edgembar::Sim *> q( mStages[i].GetSims() );
      p.insert(p.end(),q.begin(),q.end());
    }
  return p;
}




std::vector<edgembar::Trial const *> edgembar::Env::GetConstTrials() const
{
  std::vector<edgembar::Trial const *> p;
  for ( int i=0, n=mStages.size(); i<n; ++i )
    {
      std::vector<edgembar::Trial const *> q( mStages[i].GetConstTrials() );
      p.insert(p.end(),q.begin(),q.end());
    }
  return p;
}

std::vector<edgembar::Trial *> edgembar::Env::GetTrials()
{
  std::vector<edgembar::Trial *> p;
  for ( int i=0, n=mStages.size(); i<n; ++i )
    {
      std::vector<edgembar::Trial *> q( mStages[i].GetTrials() );
      p.insert(p.end(),q.begin(),q.end());
    }
  return p;
}




std::string edgembar::Env::GetPython() const
{
  std::stringstream ss;
  if ( (int)mStages.size() > 0 )
    {
      ss << "edgembar.Env("
	 << "\"" << mName << "\",[";
      for ( int i=0, n=mStages.size(); i<n; ++i )
	{
	  ss << mStages[i].GetPython();
	  if ( i < n-1 )
	    {
	      ss << ",";
	    };
	}
      ss << "])";
    }
  else
    {
      ss << "None";
    };
  
  return ss.str();
}
