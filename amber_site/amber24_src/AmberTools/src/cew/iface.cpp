#include "iface.hpp"
#include "Ewald.hpp"

#include <memory>
#include <cstdlib>
#include <iostream>

#ifdef WITH_MPI
#include <mpi.h>
#endif

std::shared_ptr< cew::Ewald > pCEW;



void cew_init_
(
#ifdef WITH_MPI
 int const * CommGroup,
#endif
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
 double const * p_cut  )
{
  pCEW.reset( new cew::Ewald
	      (
#ifdef WITH_MPI
	       *CommGroup,
#endif
	       p_natom,
	       p_nquant,
	       p_nquant_nlink,
	       p_iqmatoms,
	       p_charges,
	       p_qm_charges,
	       p_nucq,
	       p_nfft1,
	       p_nfft2,
	       p_nfft3,
	       p_order,
	       p_ew_coeff,
	       p_cut  ) );

}


void cew_prescf_
( double const * sander__coords_after_lnk_fix,
  double const * nblist__ucell,
  double * qmcrd,
  int * nmm, double * mmcharges, double * mmcrd )
{
  if ( pCEW )
    {
      pCEW->prescf( sander__coords_after_lnk_fix,
		    nblist__ucell,
		    qmcrd,
		    nmm, mmcrd, mmcharges );
    }
  else
    {
      std::cerr << "cew_prescf_ called without initializing cew" << std::endl;
      std::abort();
    }
}



void cew_postscf_
( double const * qmfrc,
  double const * mmfrc,
  double * eamber,
  double * frcamber )
{
  if ( pCEW )
    {
      pCEW->postscf( qmfrc, mmfrc, eamber, frcamber );
    }
  else
    {
      std::cerr << "cew_postscf_ called without initializing cew" << std::endl;
      std::abort();
    }
}



void cew_getpotatpt_
( double const * pt,
  double * pot )
{
  *pot = pCEW->CptPlaneWavePot(pt);
}

void cew_getgrdatpt_
( double const * pt,
  double * grd )
{
  pCEW->CptPlaneWaveGrd(pt,grd);
}
  
void cew_accdensatpt_
( double const * pt,
  double const * dens,
  double * grd )
{
  pCEW->AccTmpDensityControlPts(pt,*dens,grd);
}
