#ifndef _ndfes_vFEPData_hpp_
#define _ndfes_vFEPData_hpp_

#include <vector>

#include "DimInfo.hpp"
#include "SpatialBin.hpp"
#include "State.hpp"
#include "Sample.hpp"

namespace ndfes
{

  class vFEPData
  {
  public:

    vFEPData();
    
    vFEPData( ndfes::DimInfo const & info,
	      std::size_t const insubbins,
	      std::vector< ndfes::State > const & istates,
	      std::vector< ndfes::SpatialBin > const & isbins );

    void reset( ndfes::DimInfo const & info,
		std::size_t const insubbins,
		std::vector< ndfes::State > const & istates,
		std::vector< ndfes::SpatialBin > const & isbins );
    
    void CornerValsToSubCellVals( std::size_t const ibin,
				  double const * cornervals,
				  double * subcellvals ) const;
    
    void SubCellGrdsToCornerGrds( std::size_t const ibin,
				  double const * subcellvals,
				  double * cornervals ) const;
    
    std::vector<double> CptLinearTerm
    ( std::vector< ndfes::Sample > const & samples,
      std::vector< std::vector<std::size_t> > const & sbinsamples,
      std::vector< ndfes::SpatialBin > const & sbins ) const;
    
    double CptChisq( std::vector<double> const & x,
		     std::vector<double> & g,
		     std::vector<double> const & cornerhs ) const;

    void InterpvFEP( double const * c,
		     std::vector<double> const & p,
		     std::vector<double> const & dp,
		     double & val,
		     double & err ) const;
    
  public:
    
    std::size_t ndim;
    std::size_t bsplorder;
    std::size_t nsubbins;
    
    std::size_t nbspl;
    std::size_t nstates;
    std::size_t nbins;
    
    std::size_t ncorner;
    std::size_t subbinsize;

    ndfes::DimInfo diminfo;
    std::vector<std::size_t> glbcidxs;
    std::vector<std::size_t> ucidxs; // ncorner * nbins
    std::vector<double> c2sc; // ncorner * subbinsize
    std::vector<double> zwts; // nstates * subbinsize * nbins
  };
  
}

#endif
