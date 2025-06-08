#ifndef _PathData_hpp_
#define _PathData_hpp_

#include <vector>
#include <memory>
#include <ostream>
#include "../xmlio/xmlio.hpp"
#include "PCurve.hpp"
#include "FES.hpp"
#include "PathOptions.hpp"
#include "TubeUtils.hpp"

namespace ndfes
{
  
  
  class PathSpl
  {
  public:
    
    PathSpl();

    PathSpl( XMLIO::xml_node node );
    
    PathSpl( int const stype, bool const smooth,
	     std::size_t const ndim, std::size_t const npts, 
	     double const * ipts );

    PathSpl( int const stype, bool const smooth,
	     std::size_t const ndim, std::size_t const npts, 
	     double const * ipts, double const * ts );

    std::vector<double> GetUniformPts( std::size_t const npts ) const;

    std::vector<double> GetControlPts() const;
    
    void GetValue( double const t, double * x ) const;

    void WriteXml( std::ostream & cout, int const ilev ) const;

    
    
  public:
    
    int mSType;
    bool mSmooth;
    std::size_t mNumDim;
    std::size_t mNumPts;
    std::vector<double> mPts;
    std::vector<double> mT;
    std::shared_ptr< ndfes::PCurve > mSpl;

  private:

    void ResetSpl();
  };


  class PathSims
  {
  public:
    
    PathSims();

    PathSims( std::size_t const ndim, std::size_t const nsim,
	      double const * prcs, double const * pfcs );

    PathSims( XMLIO::xml_node node );

    void WriteXml( std::ostream & cout, int const ilev ) const;


    std::size_t mNumDim;
    std::size_t mNumSim;
    std::vector<double> mRCs;
    std::vector<double> mFCs;
    
  };


  class PathIter
  {
  public:

    PathIter() {};

    PathIter( ndfes::PathSpl const & pspl ) : mPath(pspl) {}
    
    PathIter( XMLIO::xml_node node );
    
    void WriteXml( std::ostream & cout, int const ilev ) const;
    
    std::vector<double> PredictCentroids
    ( std::shared_ptr<ndfes::FES> fes,
      ndfes::PathOptions const & popts ) const;
    
    ndfes::PathIter Next
    ( std::shared_ptr<ndfes::FES> fes,
      ndfes::PathOptions const & popts ) const;

    
    void PredictUniformSims
    ( std::size_t const nsim,
      ndfes::PathOptions const & popts );
    
    void PredictUniformSims
    ( std::size_t const nsim,
      std::shared_ptr<ndfes::FES> fes,
      ndfes::PathOptions const & popts );
    
    void PredictUniformCentroids
    ( std::size_t const nsim,
      std::shared_ptr<ndfes::FES> fes,
      ndfes::PathOptions const & popts,
      std::vector<double> & ts );

    std::vector<ndfes::AreaSummary> PredictGridBasedSims
    ( std::size_t const nsim,
      std::shared_ptr<ndfes::FES> fes,
      ndfes::PathOptions const & popts,
      int const nlayers,
      std::vector<double> & reservedpts );

    
    ndfes::PathSpl mPath;
    ndfes::PathSims mSims;
  };


  class PathOpt
  {
  public:

    PathOpt() {};
    
    PathOpt( std::string fname );

    void WriteXml( std::ostream & cout, int const ilev ) const;

    void WriteXml( std::string fname ) const;

    bool CheckRepeatedPath( std::shared_ptr<ndfes::FES> fes ) const;
    
    std::vector<ndfes::PathIter> mIters;
  };
  
}

#endif
