#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <iostream>
#include "ReadInput.hpp"
#include "ReadMode.hpp"
#include "Edge.hpp"
#include "../xmlio/xmlio.hpp"

#ifdef WITH_OPENMP
#include <omp.h>
#endif

std::shared_ptr<edgembar::Edge>
edgembar::ReadInput
( std::string fname,
  double const beta,
  edgembar::ReadMode mode,
  double const fstart,
  double const fstop,
  int const stride,
  double const ptol,
  int const autoeqmode )
{
  std::ifstream cin;
  cin.open(fname.c_str());
  if ( ! cin )
    {
      std::cerr << "File not found " << fname << std::endl;
      std::exit(EXIT_FAILURE);
    }

  std::string edgename;
  edgembar::Env cmpl;
  edgembar::Env solv;

  
  XMLIO::xml_document xml;
  XMLIO::LoadXml( fname, xml );

  XMLIO::xml_node edgenode( XMLIO::ExtractNode(xml,"edge") );
  XMLIO::ExtractAttribute(edgenode,"name",edgename);


  //
  // Loop over <env>
  //
  
  for ( XMLIO::xml_node envnode = edgenode.first_child();
	envnode; envnode = envnode.next_sibling() )
    {
      std::string nodename( envnode.name() );
      if ( nodename.compare("env") != 0 )
	{
	  std::cerr << "Expected <env> node, but found "
		    << nodename << "\n";
	  std::cerr << XMLIO::FamilyLine(envnode) << "\n";
	  std::exit(EXIT_FAILURE);
	};
      std::string envname;
      XMLIO::ExtractAttribute(envnode,"name",envname);

      if ( envname.compare("target") != 0 and
	   envname.compare("reference") != 0 )
	{
	  std::cerr << "Invalid environment encountered: " << envname << "\n";
	  std::exit(EXIT_FAILURE);
	}


      std::vector<edgembar::Stage> stages;
      
      //
      // Loop over <stage>
      //
      
      for ( XMLIO::xml_node stagenode = envnode.first_child();
	    stagenode; stagenode = stagenode.next_sibling() )
	{
	  nodename = stagenode.name();
	  if ( nodename.compare("stage") != 0 )
	    {
	      std::cerr << "Expected <stage> node, but found "
			<< nodename << "\n";
	      std::cerr << XMLIO::FamilyLine(stagenode) << "\n";
	      std::exit(EXIT_FAILURE);
	    };
	  std::string stagename;
	  XMLIO::ExtractAttribute(stagenode,"name",stagename);


	  
	  
	  //
	  // Loop over <trial>
	  //
	  std::vector<edgembar::Trial> trials;
	  
	  for ( XMLIO::xml_node trialnode = stagenode.first_child();
		trialnode; trialnode = trialnode.next_sibling() )
	    {
	      nodename = trialnode.name();
	      if ( nodename.compare("trial") != 0 )
		{
		  std::cerr << "Expected <trial> node, but found "
			    << nodename << "\n";
		  std::cerr << XMLIO::FamilyLine(trialnode) << "\n";
		  std::exit(EXIT_FAILURE);
		};
	      std::string trialname;
	      XMLIO::ExtractAttribute(trialnode,"name",trialname);
	      edgembar::ReadMode tmode = mode;
	      if ( trialnode.attribute("mode") )
		{
		  tmode = edgembar::GetMode( trialnode.attribute("mode").value() );
		}
	      
	      //
	      // Read <dir> and <ene> 
	      //
	      
	      std::string datadir;
	      std::vector<std::string> ene;
	      double shift = 0;
	  
	      for ( XMLIO::xml_node datanode = trialnode.first_child();
		    datanode; datanode = datanode.next_sibling() )
		{
		  nodename = datanode.name();
		  if ( nodename.compare("dir") == 0 )
		    {
		      datadir = datanode.child_value();
		    }
		  else if ( nodename.compare("ene") == 0 )
		    {
		      ene.push_back( datanode.child_value() );
		    }
		  else if ( nodename.compare("shift") == 0 )
		    {
		      std::istringstream tstr( datanode.child_value() );
		      if ( ! (tstr >> shift) )
			{
			  std::cerr << "Could not convert shift to float '"
				    << datanode.child_value() << "' while "
				    << "processing:\n";
			  std::cerr << XMLIO::FamilyLine(datanode) << "\n";
			  std::exit(EXIT_FAILURE);
			}
		    }
		  else
		    {
		      std::cerr << "Expected <dir> or <ene> nodes, but "
				<< "found " << nodename << "\n";
		      std::cerr << XMLIO::FamilyLine(datanode) << "\n";
		      std::exit(EXIT_FAILURE);
		    };
		  
		};

	      if ( datadir.size() == 0 )
		{
		  std::cerr << "Trial is missing a <dir> node\n";
		  std::cerr << XMLIO::FamilyLine(trialnode) << "\n";
		  std::exit(EXIT_FAILURE);
		};
	      
	      if ( ene.size() < 2 )
		{
		  std::cerr << "Trial is missing <ene> nodee\n";
		  std::cerr << XMLIO::FamilyLine(trialnode) << "\n";
		  std::exit(EXIT_FAILURE);
		};

	      trials.push_back( edgembar::Trial(trialname,ene,datadir,
						beta,tmode,autoeqmode,
						beta*shift) );
	      

	    };

	  stages.push_back( edgembar::Stage(stagename,trials) );

	};

      if ( envname.compare("target") == 0 )
	{
	  cmpl = edgembar::Env("target",stages);
	}
      else
	{
	  solv = edgembar::Env("reference",stages);
	}
      
    };

  std::shared_ptr<edgembar::Edge> edge( new edgembar::Edge(edgename,cmpl,solv) );

  //
  // Note: This could/should be placed inside the Edge constructor
  // I'd need to update the constructor arguments to accept fstart, fstop
  //
  {
    std::vector<edgembar::Sim *> sims( edge->GetSims() );
    int const n = sims.size();
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int i=0; i<n; ++i )
      {
	sims[i]->ReadFiles( fstart, fstop, stride );
      }
    
    std::vector<edgembar::Trial *> trials( edge->GetTrials() );
    int const nt = trials.size();
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int i=0; i<nt; ++i )
      {
	trials[i]->StoreAvgEnes();
      }
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int i=0; i<n; ++i )
      {
	sims[i]->PrecomputeExps( ptol );
      }
    edge->StoreFreeAvgEnes();
  }

  /*
  {
    //std::shared_ptr<edgembar::Edge> tedge( edge->ExtractRange(0,0.5) );
    std::vector<edgembar::Trial *> trials( edge->GetTrials() );
    int const n = trials.size();
    for ( int i=n-1; i<n; ++i )
      {
	trials[i]->WriteDebugInfo( std::cout );
	std::vector<edgembar::Sim *> osims( trials[i]->GetSims() );
	for ( int j=0, m=osims.size(); j<m; ++j )
	  {
	    osims[j]->WriteDebugInfo( std::cout );
	  };
      }

    trials[n-1]->TestObjective();
  }
  */
  
  
  return edge;
}
