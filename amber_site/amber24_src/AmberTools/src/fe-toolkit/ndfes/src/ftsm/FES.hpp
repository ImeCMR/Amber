#ifndef _ndfes_FES_hpp_
#define _ndfes_FES_hpp_

#include <vector>
#include <memory>
#include <unordered_map>
#include "DimInfo.hpp"
#include "SpatialBin.hpp"

namespace ndfes
{
  class FES
  {
  public:
    
    FES();

    FES( std::string xmlname, int const imodel );


    void SetupRBF( double const shape,
		   std::size_t const minsize,
		   double const maxerr,
		   bool const with_inverse );
    void SetupHist();
    void SetupBspl();
    void SetupARBF( int const delta, bool const with_inverse );
    void PreparevFEP();
    void SetupWAvg( int const order, int const niter );
    
    void EvalRBF( std::size_t const npt, double const * pts, double * vals ) const;
    void EvalRBF( std::size_t const npt, double const * pts, double * vals, double * grds ) const;

    void EvalRBF( std::vector<double> const & cpts, std::vector<double> const & cwts, std::size_t const npts, double const * pts, double * vals ) const;
    void EvalRBF( std::vector<double> const & cpts, std::vector<double> const & cwts, std::size_t const npts, double const * pts, double * vals, double * grds ) const;

    void EvalWAvgValue( double const * pt, double & ene ) const;
    void EvalWAvgValueAndError( double const * pt, double & ene, double & err ) const;
    void EvalWAvgValueAndGrad( double const * pt, double & ene, double * grd ) const;


    std::size_t BinOccAtPt( double const * pt ) const;
    bool PointIsOcc( double const * pt ) const;
    bool PointIsInRange( double const * pt ) const;
    bool PointIsInRangeAndNonzeroOcc( double const * pt ) const;

    void EvalQuadInBin( std::size_t ibin,
			std::vector<double> & work,
			double * vals );


    void SetupQuad( int nquad );
    void SetupNoQuad();

    void BootstrapFES( double const booterror );


    bool GetValue( double const * pt, double & val ) const;
    
    bool GetValue( double const * pt, double & val, double * grd ) const;
    
    bool GetPenalizedValue( double const oobk, double const * pt, double & val, double * grd ) const;
    
    bool GetBiasedValue( double const * rc, double const * fc,
			 double const oobk, double const * pt,
			 double & val, double * grd ) const;

    bool GetValueErrorAndEntropy( double const * pt,
				  double & val,
				  double & err,
				  double & re,
				  std::size_t & ns ) const;

    std::vector<double> GetModifiedBinValues() const;
    std::vector<double> GetModifiedBinValues( std::vector<double> const & binvals ) const;

    std::shared_ptr<ndfes::FES> DeleteBinsBySize( std::size_t const minsize ) const;

    
    std::shared_ptr<ndfes::FES> AddPenaltyBuffer( std::size_t const minsize ) const;
    std::shared_ptr<ndfes::FES> AddPenaltyBuffer_ORIG( std::size_t const minsize ) const;

    std::shared_ptr<ndfes::FES> AddPenaltyBuffer( std::size_t const minsize, std::size_t const nlayers ) const;

    
    int NumLayersOfOcc( double const * pt ) const;

    void AddOccPen( double const occpen );

    double CptAvgOccAtPt( double const * pt, std::size_t const order ) const;
    
    
    bool mIsMBAR;
    int mInterpMode;
    std::size_t mNumBins;
    std::size_t mOrder;
    std::size_t mNBspl;
    std::size_t mNumCorners;
    
    ndfes::DimInfo mDimInfo;
    std::vector<ndfes::SpatialBin> mBins;

    std::vector<std::size_t> mGlbCidxs; // total num corners
    std::vector<double> mGlbCvals; // total num corners
    std::vector<double> mGlbCerrs; // total num corners
    std::vector<std::size_t> mLocCidxs; // mNumBins
    
    std::vector<double> mRBFCpts; 
    std::vector<double> mRBFCwts;
    std::vector<double> mRBFAinv;
    std::vector<double> mRBFCerr;
    double mRBFShape;
    
    std::vector<double> mBinVals; // mNumBins
    std::vector<double> mBinErrs; // mNumBins
    std::vector<double> mBinREs;  // mNumBins
    std::vector<std::size_t> mBinSizes; // mNumBins


    // quadrature within a single bin centered at the origin
    std::size_t mNumQuad;
    std::size_t mSubBinSize;
    std::vector<double> mQuadPts;
    std::vector<double> mQuadWts;
    std::vector<double> mC2SC;

    std::vector< std::vector<double> > mMeshFE;

    int mARBFDelta;
    bool mWithInverse;
    std::vector< std::vector<std::size_t> > mARBFCidx;
    std::vector< std::vector<double> > mARBFCpts;
    std::vector< std::vector<double> > mARBFCwts;
    std::vector< std::vector<double> > mARBFAinv;

    std::unordered_map<std::size_t,std::size_t> mGlbIdxMap;
    
    // std::vector<double> mARBFA;
    // std::vector<double> mARBFB;
    // std::vector<int> mARBFdidx;
    
  };
}

#endif
