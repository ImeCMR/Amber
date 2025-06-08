#include "../include/dprec.fh"
#include "../include/assert.fh"

module pol_gauss_mdin
  implicit none
  private

  integer :: ipgm
  integer :: pol_gauss_verbose = 0
  integer :: use_average = 0
  integer :: dipole_print = 0
  integer :: saved_nprint = 0
  integer :: dipole_scf_init = 3
  integer :: dipole_scf_init_order = 3
  integer :: dipole_scf_init_step = 2
  integer :: scf_solv_opt = 3
  integer :: scf_sor_niter = 100
  integer :: scf_cg_niter = 50
  integer :: scf_local_niter = 3
  _REAL_  :: dipole_scf_tol = 0.01d0
  _REAL_  :: scf_sor_coefficient = 0.65d0
  _REAL_  :: scf_local_cut=4.0d0
  _REAL_  :: ee_dsum_cut=9.0d0
  _REAL_  :: ee_damped_cut=4.5d0
  _REAL_  :: erfc_tol = 1.0e-5
  _REAL_  :: ee_gauss_cut = 3.12342d0 ! corresponding to max error of 10^-5 in erfc

  ! pack these into common blocks for easy bcast
  integer, parameter :: pgm_ctrl_int_cnt = 11
  integer, parameter :: pgm_ctrl_dbl_cnt = 7
  common / pgm_ctrl_int/ pol_gauss_verbose, use_average, dipole_print, saved_nprint, &
                         dipole_scf_init, dipole_scf_init_order, dipole_scf_init_step, &
                         scf_solv_opt, scf_sor_niter, scf_cg_niter, scf_local_niter
  save :: / pgm_ctrl_int /
  common / pgm_ctrl_dbl / dipole_scf_tol, scf_sor_coefficient, scf_local_cut, &
                          ee_dsum_cut, ee_damped_cut, erfc_tol, ee_gauss_cut
  save :: / pgm_ctrl_dbl / 

  public POL_GAUSS_read_mdin,ipgm,pol_gauss_verbose,dipole_print,saved_nprint, &
         dipole_scf_tol,use_average,dipole_scf_init,dipole_scf_init_order,dipole_scf_init_step,&
         scf_solv_opt,scf_cg_niter,scf_local_niter, &
         scf_sor_coefficient,scf_sor_niter, &
         ee_dsum_cut,ee_damped_cut,scf_local_cut, pgm_ctrl_int_cnt, pgm_ctrl_dbl_cnt
#ifdef MPI
  public pGM_mdin_bcast
#endif

  contains
!-------------------------------------------------------------------------------
#ifdef MPI
subroutine pGM_mdin_bcast()
  include 'mpif.h'
# include "extra.h"
# include "parallel.h"

  integer :: ibuf(pgm_ctrl_int_cnt)
  common / pgm_ctrl_int / ibuf
  save :: / pgm_ctrl_int /
  _REAL_ :: buf(pgm_ctrl_dbl_cnt)
  common / pgm_ctrl_dbl / buf
  save :: / pgm_ctrl_dbl /

  integer :: ier

  call mpi_bcast(ibuf, pgm_ctrl_int_cnt, MPI_INTEGER, 0, commsander, ier)
  call mpi_bcast(buf, pgm_ctrl_dbl_cnt, MPI_DOUBLE_PRECISION, 0, commsander, ier)
  !write(40+mytaskid, *) "MDIN: On processor", mytaskid
  !write(40+mytaskid, *) "ibuf", ibuf
  !write(40+mytaskid, *) "buf", buf
end subroutine pGM_mdin_bcast
#endif
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_read_mdin(nf)
  use file_io_dat

  integer,intent(in) :: nf

  namelist/pol_gauss/pol_gauss_verbose,use_average,dipole_print, &
                     dipole_scf_tol,dipole_scf_init,dipole_scf_init_order,dipole_scf_init_step, &
                     scf_solv_opt,scf_sor_coefficient,scf_sor_niter, &
                     scf_cg_niter,scf_local_niter,scf_local_cut, &
                     ee_dsum_cut,ee_damped_cut, &
                     erfc_tol,ee_gauss_cut

  read(nf,nml=pol_gauss)

  scf_cg_niter = max(2, scf_cg_niter)
end subroutine POL_GAUSS_read_mdin
!-------------------------------------------------------------------------------
end module pol_gauss_mdin
!-------------------------------------------------------------------------------
