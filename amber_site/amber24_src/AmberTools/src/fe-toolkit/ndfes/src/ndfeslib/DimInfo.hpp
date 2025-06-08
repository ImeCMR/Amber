#ifndef _ndfes_DimInfo_hpp_
#define _ndfes_DimInfo_hpp_

#include <vector>
#include <memory>
#include <ostream>
#include <cstddef>

namespace ndfes
{
  class DimInfo
  {
  public:
    
    DimInfo();

    DimInfo( std::size_t const ndim,
	     std::size_t const bsplorder,
	     int const * isperiodic,
	     double const * target_widths );

    DimInfo( std::size_t const ndim,
	     std::size_t const bsplorder,
	     int const * isperiodic,
	     double const * xmins,
	     double const * xmaxs,
	     int const * dimsizes );

    std::size_t GetNumDims() const { return ndim; }

    std::size_t GetBsplOrder() const { return bsplorder; }

    std::size_t GetNumCorners() const { return ncorner; }

    int const * GetIsPeriodic() const { return isperiodic.data(); }
    
    bool IsPeriodic( std::size_t const idim ) const { return isperiodic[idim] > 0; }

    double const * GetTargetWidths() const { return target_widths.data(); }

    double const * GetXmin() const { return xmin.data(); }

    double const * GetXmax() const { return xmax.data(); }

    int const * GetDimSizes() const { return dimsizes.data(); }
    
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

    void WriteChkpt( std::ostream & cout, std::size_t const indent ) const;
    
    
  private:
    
    std::size_t ndim;
    std::size_t bsplorder;
    std::size_t nbspl;
    std::size_t ncorner;
    std::vector<int> isperiodic;
    std::vector<double> target_widths;
    std::vector<double> xmin;
    std::vector<double> xmax;
    std::vector<int> dimsizes;

  };


  typedef std::shared_ptr<ndfes::DimInfo> pDimInfo;
  
}

#endif
