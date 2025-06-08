#ifndef _ndfes_SystemInfo_hpp_
#define _ndfes_SystemInfo_hpp_

#include <cstddef>
#include <iostream>
#include <string>

#include "DimInfo.hpp"
#include "State.hpp"
#include "Sample.hpp"
#include "SpatialBin.hpp"
#include "RunAvg.hpp"
#include "vFEPData.hpp"

namespace ndfes
{
  class SystemInfo
  {
  public:
    
    SystemInfo();
    
    SystemInfo( ndfes::DimInfo const & dinfo,
		double const targetTemp,
		double const DeltaDoS,
		bool const sdosFullSampling,
		std::string const dosHist,
		std::size_t const numhams );


    void PrintSummary( std::ostream & cout ) const;

    void WriteChkpt_MBAR
    ( std::string fname,
      std::vector< std::vector<double> > const & Fbins,
      std::vector< std::vector<double> > const & Ebins,
      std::vector< std::vector<double> > const & Sbins ) const;

    void WriteChkpt_vFEP
    ( std::string fname,
      std::vector<double> const & Fs,
      std::vector<double> const & Es ) const;

    void WriteChkptXML_MBAR
    ( std::string fname,
      std::vector< std::vector<double> > const & Fbins,
      std::vector< std::vector<double> > const & Ebins,
      std::vector< std::vector<double> > const & Sbins ) const;

    void WriteChkptXML_vFEP
    ( std::string fname,
      std::vector<double> const & Fs,
      std::vector<double> const & Es ) const;

    
    std::size_t InsertState( std::size_t const iham,
		     double const * center,
		     double const * fconst,
		     double const temperature );
    
    std::size_t InsertState( ndfes::State const & s );


    void InsertSample( ndfes::Sample const & s ) { samples.push_back(s); }

    void SetStatIneff( std::size_t const istate, std::size_t const g );
    
    void CptSpatialHistogram();
    
    void ShiftPotEnesToAvg();

    
    ndfes::DimInfo const & GetDimInfo() const { return diminfo; }

    double GetBeta() const { return beta; }

    void SetDeltaDoS( double const d ) { deltaDoS = d; }

    double GetDeltaDoS() const { return deltaDoS; }
    
    std::string GetHistPrefix() const { return sdosHistPrefix; }
    
    void SetHistPrefix( std::string p ) { sdosHistPrefix=p; }
    
    void SetTargetTemp( double const T );
    
    double GetTargetTemp() const { return temperature; }

    double ConvertKcalToKT( double const E ) const { return beta*E; }

    double ConvertKTToKcal( double const E ) const { return E/beta; }


    
    std::size_t GetNumStates() const { return states.size(); }
    
    ndfes::State const * GetStates() const { return states.data(); }
        
    ndfes::State const & GetState( std::size_t i ) const { return states[i]; }

    
    std::size_t GetNumSamples() const { return samples.size(); }
    
    ndfes::Sample const * GetSamples() const { return samples.data(); }
        
    ndfes::Sample const & GetSample( std::size_t i ) const { return samples[i]; }


    std::size_t GetNumBins() const { return sbins.size(); }

    ndfes::SpatialBin const * GetBins() const { return sbins.data(); }

    ndfes::SpatialBin const & GetBin( std::size_t i ) const { return sbins[i]; }
    
    
    void CptBiasEnergies( double * W, double const * pt ) const;


    void ResampleStates( std::vector< std::vector<std::size_t> > & ssamples,
			 std::vector< std::vector<std::size_t> > & bsamples ) const;
    
    void ResampleBin( std::size_t ibin, std::vector<std::size_t> & bsamples ) const;

    void GetSampleIdxs( std::vector< std::vector<std::size_t> > & ssamples,
			std::vector< std::vector<std::size_t> > & bsamples ) const;
    
    double CptFastMBARObjective
    ( std::vector<double> const & bs,
      std::vector< std::vector<std::size_t> > const & idxs ) const;
    
    double CptFastMBARObjective
    ( std::vector<double> const & bs,
      std::vector<double> & gs,
      std::vector< std::vector<std::size_t> > const & idxs ) const;

    void CptWTP
    ( std::vector<double> const & Fin,
      std::vector< std::vector<std::size_t> > const & spts,
      std::vector< std::vector<std::size_t> > const & bpts,
      std::vector<double> & Fout,
      std::vector<double> & Fham,
      std::vector< std::vector<double> > & Fbins,
      std::vector< std::vector<double> > & Sbins ) const;

    std::size_t InitvFEP( std::size_t const nsubbins );

    std::vector<double> CptvFEPLinearTerm
    ( std::vector< std::vector<std::size_t> > const & bpts ) const;

    std::size_t GetParamCornerGlbIdx( std::size_t ic ) const;
    
    double CptvFEPObjective
    ( std::vector<double> const & bs,
      std::vector<double> & gs,
      std::vector<double> const & cornerh ) const;

    void InterpvFEP
    ( double const * c,
      std::vector<double> const & p,
      std::vector<double> const & dp,
      double & val,
      double & err ) const;

    void CheckForConsistentTemperatures( bool expert ) const;
    
  private:
    
    std::vector<ndfes::SpatialBin>::iterator FindBin( std::size_t const gidx );
    
    std::vector<ndfes::SpatialBin>::const_iterator FindBin( std::size_t const gidx ) const;
    
  private:
    
    ndfes::DimInfo diminfo;
    double temperature;
    double beta;
    double deltaDoS;
    bool sdosFullSampling;
    std::string sdosHistPrefix;
    std::size_t nham;
    std::vector< ndfes::State > states;
    std::vector< ndfes::Sample > samples;
    std::vector< ndfes::SpatialBin > sbins;
    std::vector<std::size_t> ucidxs;
    std::vector< std::vector<std::size_t> > statesamples;
    std::vector< std::vector<std::size_t> > sbinsamples;
    std::shared_ptr<ndfes::vFEPData> vfep;

};

typedef std::shared_ptr< ndfes::SystemInfo > pSystemInfo;

}

#endif

