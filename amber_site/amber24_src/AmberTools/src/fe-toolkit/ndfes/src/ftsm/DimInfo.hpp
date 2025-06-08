#ifndef _ndfes_DimInfo_hpp_
#define _ndfes_DimInfo_hpp_

#include <vector>
#include <memory>
#include <ostream>
#include <cstddef>


namespace ndfes
{
  int ifloor( double const x );
  int iceil( double const x );
}


namespace ndfes
{
  class DimInfo
  {
  public:
    
    DimInfo();

    DimInfo( int const ndim,
	     int const bsplorder,
	     int const * isperiodic,
	     double const * target_widths );

    DimInfo( int const ndim,
	     int const bsplorder,
	     int const * isperiodic,
	     double const * xmins,
	     double const * xmaxs,
	     int const * dimsizes );

    int GetNumDims() const { return ndim; }

    int GetBsplOrder() const { return bsplorder; }

    int GetNumCorners() const { return ncorner; }

    int const * GetIsPeriodic() const { return isperiodic.data(); }
    
    bool IsPeriodic( int const idim ) const { return isperiodic[idim] > 0; }

    double const * GetTargetWidths() const { return target_widths.data(); }

    double const * GetXmin() const { return xmin.data(); }

    double const * GetXmax() const { return xmax.data(); }

    int const * GetDimSizes() const { return dimsizes.data(); }


    void CptBinCenterFromPt( double const * pt, double * c ) const;
    
    void ResetRange();
    
    void ModifyRange( double const * pt );

    void FinalizeRange();

    std::vector<std::size_t> CptBinIdxs( double const * pt ) const;
    
    //std::vector<int> CptBinIdxs( int const glbbinidx ) const;

    std::size_t CptGlbBinIdx( std::size_t const * binidxs ) const;

    std::size_t CptGlbBinIdx( double const * pt ) const;

    std::vector< std::vector<std::size_t> > CptCornerIdxs( std::size_t const * binidxs ) const;
    
    std::vector<std::size_t> CptCornerGlbIdxs( std::vector< std::vector<std::size_t> > const & binidxs ) const;

    std::vector<std::size_t> CptCornerGlbIdxs( std::size_t const * binidxs ) const;

    std::vector<double> CptBinCenter( std::size_t const * binidxs ) const;    

    void WriteChkpt( std::ostream & cout, int const indent ) const;


    // void Wrap( double * crd ) const;

    // void WrapDiffToOrigin( double * crd, double const * origin ) const;
    
    
  private:
    
    int ndim;
    int bsplorder;
    int nbspl;
    int ncorner;
    std::vector<int> isperiodic;
    std::vector<double> target_widths;
    std::vector<double> xmin;
    std::vector<double> xmax;
    std::vector<int> dimsizes;

  };


  typedef std::shared_ptr<ndfes::DimInfo> pDimInfo;


  void UnwrapPath( ndfes::DimInfo const & info, std::vector<double> & pts );
  
  void WrapPath( ndfes::DimInfo const & info, std::vector<double> & pts );

  void UnwrapCentroids( ndfes::DimInfo const & info,
			std::vector<double> & centroids,
			std::vector<double> const & rcs );
  
  void WrapPathToNearestOccBin
  ( ndfes::DimInfo const & info,
    std::vector< std::vector<std::size_t> > const & mbinidxs,
    std::vector<double> & pts );
  
  
}

#endif
