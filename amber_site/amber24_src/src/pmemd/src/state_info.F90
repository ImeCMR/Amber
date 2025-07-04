#include "copyright.i"

!*******************************************************************************
!
! Module:  state_info_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module state_info_mod

  implicit none

  ! This is basically an enumeration of state info placed in state info arrays,
  ! in an effort to improve maintainability/clarity of the code.  This is
  ! unfortunately the best solution available, as there is no good mechanism
  ! for performing the same operation across all elements of a record (type),
  ! and we must be able to do total and running averages and rmsd on this
  ! data

  integer, parameter    :: si_tot_ene = 1
  integer, parameter    :: si_kin_ene = 2
  integer, parameter    :: si_solute_kin_ene = 3
  integer, parameter    :: si_solvent_kin_ene = 4
  integer, parameter    :: si_volume = 5
  integer, parameter    :: si_tot_press = 6
  integer, parameter    :: si_tot_ekcmt = 7
  integer, parameter    :: si_tot_virial = 8
  integer, parameter    :: si_pot_ene = 9
  integer, parameter    :: si_vdw_ene = 10
  integer, parameter    :: si_elect_ene = 11
  integer, parameter    :: si_hbond_ene = 12
  integer, parameter    :: si_bond_ene = 13
  integer, parameter    :: si_angle_ene = 14
  integer, parameter    :: si_dihedral_ene = 15
  integer, parameter    :: si_vdw_14_ene = 16
  integer, parameter    :: si_elect_14_ene = 17
  integer, parameter    :: si_restraint_ene = 18
  integer, parameter    :: si_dvdl = 19
  integer, parameter    :: si_density = 20
  integer, parameter    :: si_pme_err_est = 21
  integer, parameter    :: si_angle_ub_ene = 22
  integer, parameter    :: si_dihedral_imp_ene = 23
  integer, parameter    :: si_cmap_ene = 24
  integer, parameter    :: si_surf_ene = 25
  integer, parameter    :: si_gamma_ten = 26
  integer, parameter    :: si_amd_boost = 27
  integer, parameter    :: si_gamd_boost = 28
  integer, parameter    :: si_emap_ene = 29
  integer, parameter    :: si_efield_ene = 30
  integer, parameter    :: si_gamd_ppi   = 31
  integer, parameter    :: si_gamd_ppi2   = 32
  integer, parameter    :: si_phmd_ene    = 33
  integer, parameter    :: si_cnt = 34         ! update as needed

end module state_info_mod
