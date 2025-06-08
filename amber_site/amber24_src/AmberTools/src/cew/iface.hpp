#ifndef _cew_iface_hpp_
#define _cew_iface_hpp_

extern "C"
{
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
   double const * p_cut  );

  
  // iqmatoms (in): length nquant_nlink, and is zero-based indexing
  // charges (in): length natom in atomic units
  // qm_charges (in): length nquant_nlink in atomic units
  // nucq (in): length nquant_nlink in atomic units
  // ew_coeff (in): scalar in atomic units
  // cut (in): scalar in atomic units
  // other quantities are input scalar integers

  
  void cew_prescf_
  ( double const * sander__coords_after_lnk_fix,
    double const * nblist__ucell,
    double * qmcrd,
    int * nmm, double * mmcharges, double * mmcrd );

  // coords: length natom in angstrom
  // ucell: length 9 in angstrom
  // qmcrd (out): length 9 returned in angstrom
  // nmm (out): length 1 
  // mmcharges (out): length nmm in atomic units
  // mmcrd (out): length 3*nmm in angstrom
  
  
  void cew_postscf_
  ( double const * qmgrd,
    double const * mmgrd,
    double * eamber,
    double * frcamber );

  // qmgrd (in): length 3*nquant_nlink in kcal/mol/A
  // mmgrd (in): length 3*nmm in kcal/mol/A
  // eamber (inout): length 1 in kcal/mol
  // frc (inout): frc in kcal/mol/A

  
  void cew_getpotatpt_
  ( double const * pt,
    double * pot );
  
  // pt (in): length 3 in atomic units
  // pot (out): length 1 in atomic units


  void cew_getgrdatpt_
  ( double const * pt,
    double * grd );
  
  // pt (in): length 3 in atomic units
  // grd (out): length 3 in atomic units

  
  void cew_accdensatpt_
  ( double const * pt,
    double const * dens,
    double * grd );
  
  // pt (in): length 3 in atomic units
  // dens (out): length 1 in atomic units
  // grd (out): length 3 in atomic units


  
  /*
    int    const * sander__natom,
    int    const * qmmm_struct__nquant, 
    int    const * qmmm_struct__nlink, 
    int    const * qmmm_struct__iqm_atomic_numbers,
    int    const * qmmm_struct__iqmatoms,
    int    const * qmmm_struct__link_pairs,
    double const * qmmm_struct__scaled_mm_charges,
    double const * qmmm_struct__qm_resp_charges,
    double const * qmmm_struct__mm_link_pair_resp_charges,
    double const * ew_box__ew_coeff,
    int    const * ew_box__nfft1,
    int    const * ew_box__nfft2,
    int    const * ew_box__nfft3,
      int    const * ew_box__order,
      double const * qmmm_nml__qmcut,
      int   const * qmmm_mpi__commqmmm
  */
  
  
}

#endif
