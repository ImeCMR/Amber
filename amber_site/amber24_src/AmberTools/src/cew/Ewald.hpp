#ifndef _cew_ewald_hpp_
#define _cew_ewald_hpp_

#include <complex>
#include <vector>
#include <memory>
#include <fftw3.h>

// # include "ScfOptions.hpp"


#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace cew
{

  class Ewald
  {

  public:

#ifdef WITH_MPI
    Ewald( int const CommGroup,
	   int const * p_natom,
	   int const * p_nquant,
	   int const * p_nquant_nlink,
	   int const * p_iqmatoms,
	   double const * p_charges,
	   double const * p_qm_charges,
	   double const * p_nucq,
	   int const * p_nfft1,
	   int const * p_nfft2,
	   int const * p_nfft3,
	   int const * p_order,
	   double const * p_ew_coeff,
	   double const * p_cut  );
#else
    Ewald( int const * p_natom,
	   int const * p_nquant,
	   int const * p_nquant_nlink,
	   int const * p_iqmatoms,
	   double const * p_charges,
	   double const * p_qm_charges,
	   double const * p_nucq,
	   int const * p_nfft1,
	   int const * p_nfft2,
	   int const * p_nfft3,
	   int const * p_order,
	   double const * p_ew_coeff,
	   double const * p_cut );
#endif

    ~Ewald();

    void SetCrd
    ( double const * sander__coords_after_lnk_fix,
      double const * nblist__ucell );

#ifdef WITH_MPI
    void SetCrd();
#endif


    void prescf( double const * sander__coords_after_lnk_fix,
		 double const * nblist__ucell,
		 double * p_qmcrd,
		 int * nmm, double * mmcrd, double * mmcharges );

    void postscf( double const * qmgrd,
		  double const * mmgrd,
		  double * eamber,
		  double * frcamber );
    
    double EnergyAndGradients();

    void AccTmpDensityControlPts
    ( double const * coord, 
      double const chg,
      double * grd );

    void AccDensityControlPts
    ( double const * coord, 
      double const chg );

    void AccDensityControlPts
    ( int const natom, 
      double const * coord, 
      double const * chg );

    void CptPlaneWavePotControlPts();

    void CptPlaneWavePotControlPts
    ( int const natom, 
      double const * coord, 
      double const * chg );

    double CptPlaneWavePot( double const * pt );

    double CptPlaneWavePotAndGrd( double const * pt, double * grd );

    void CptPlaneWaveGrd( double const * pt, double * grd );

    //void SetFrameWrites( std::string fname, int freq );
    //void SetStep( int nstep );
    

#ifdef WITH_MPI
    MPI_Comm CommGroup;
    int mpiSize,mpiRank;
#endif

    //std::tr1::shared_ptr< cew::ScfSolver > scf;
    //cew::ScfOptions * options;
    // pointers to options
    int nat,nqm;
    int nquant;
    std::vector<int> iqmatoms;
    std::vector<double> charges;
    std::vector<double> qm_charges;
    std::vector<double> qmq;
    std::vector<double> nucq;

    std::vector<double> msgbuf;
    // pointers to msgbuf
    double * ucell_changed;
    double * ucell;
    double * crd;
    double * frc;
    double * imgcrd;
    double * qmcrd;

    // near-field arrays
    std::vector<double> qmc;
    std::vector<double> qmp;
    std::vector<int> immatoms;
    std::vector<double> mmc;
    std::vector<double> mmq;

    // fft
    bool opt_infl;
    std::vector<double> recip;
    int n1,n2,n3,nt,nk3,nk,order;
    double beta,zeta,qmcut,V;
    std::vector<double> w1,w2,w3,bf1,bf2,bf3;
    std::vector<double> kernel;
    std::vector<double> localfrc;
    double localene,ChargeCorrection;

    double internalene;
    
    // fftw
    std::vector<double> tmp_grid;
    double * grid;
    std::complex<double> * dftcoef;
    fftw_plan * fplan;
    fftw_plan * rplan;

    //bool write_qmmm_frame;
    //int qmmm_frame_freq;
    //std::string qmmm_frame_filename;
    //int mdstep;
    
  private:

    void UpdateKernel();
    void ImageCrdAroundQM();


  };

}




#endif
