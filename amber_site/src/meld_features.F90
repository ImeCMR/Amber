!*******************************************************************************
!
! Module:   meld_features_mod
!
! Description: This module will contain the features related to MELD implementation
!
!*******************************************************************************

module meld_features_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  output_sorted_restraint_energies
!
! Description: <TBS>
!
!*******************************************************************************

  subroutine output_sorted_restraint_energies(nmrnum, nstep, nstlim, restraint_data)
    implicit none
    integer, intent(in) :: nmrnum, nstep, nstlim
    double precision, intent(in) :: restraint_data(2, nmrnum)
    double precision :: temp_num, temp_energy
    integer :: i, j, ierr, nsteps_nmr
    double precision, allocatable :: sorted_data(:,:)

    nsteps_nmr = nstep + 1

    if (nsteps_nmr == nstlim) then
      allocate(sorted_data(2, nmrnum))
      sorted_data = restraint_data

      ! Sort by energy (ascending)
      do i = 1, nmrnum - 1
        do j = i + 1, nmrnum
          if (sorted_data(2, i) > sorted_data(2, j)) then
            temp_energy = sorted_data(2, i)
            sorted_data(2, i) = sorted_data(2, j)
            sorted_data(2, j) = temp_energy

            temp_num = sorted_data(1, i)
            sorted_data(1, i) = sorted_data(1, j)
            sorted_data(1, j) = temp_num
          end if
        end do
      end do

      open(unit=10, file="restraint_energies.txt", status="unknown", position="append", action="write", iostat=ierr)
      if (ierr /= 0) then
        write(*, *) "Error opening file for restraint energies"
      else
        do i = 1, nmrnum
          write(10, '(F10.4, F10.4)') sorted_data(1, i), sorted_data(2, i)
        end do
        close(10)
      end if

      deallocate(sorted_data)
    end if
  end subroutine output_sorted_restraint_energies


!*******************************************************************************
!
! Subroutine:  select_active_restraints
!
! Description: Given energies for each restraint, select the top fraction as active.
!
!*******************************************************************************

subroutine select_active_restraints(nmrnum, restraint_energies, fraction, active)
  implicit none
  integer, intent(in) :: nmrnum
  double precision, intent(in) :: restraint_energies(nmrnum)
  double precision, intent(in) :: fraction ! e.g. 0.6 for 60%
  logical, intent(out) :: active(nmrnum)
  integer :: i, j, n_active
  integer, allocatable :: idx(:)
  double precision, allocatable :: energies_copy(:)

  allocate(idx(nmrnum))
  allocate(energies_copy(nmrnum))
  energies_copy = restraint_energies
  do i = 1, nmrnum
    idx(i) = i
  end do

  ! Simple selection sort for indices by energy (descending)
  do i = 1, nmrnum-1
    do j = i+1, nmrnum
      if (energies_copy(i) < energies_copy(j)) then
        call swap(energies_copy(i), energies_copy(j))
        call swap(idx(i), idx(j))
      end if
    end do
  end do

  !n_active = int(fraction * nmrnum + 0.5)
  n_active = ceiling(fraction * nmrnum)
  active(:) = .false.
  do i = 1, n_active
    active(idx(i)) = .true.
  end do

  deallocate(idx, energies_copy)
end subroutine select_active_restraints

! Helper swap subroutine
subroutine swap(a, b)
  implicit none
  double precision, intent(inout) :: a
  double precision :: tmp
  tmp = a
  a = b
  b = tmp
end subroutine swap

subroutine swap(a, b)
  implicit none
  integer, intent(inout) :: a
  integer :: tmp
  tmp = a
  a = b
  b = tmp
end subroutine swap

!*******************************************************************************
!
! Subroutine:  update_active_restraints
!
! Description: Store the logical array of active restraints in a module variable.
!
!*******************************************************************************

logical, allocatable, save :: meld_active_restraints(:)

subroutine update_active_restraints(nmrnum, active)
  implicit none
  integer, intent(in) :: nmrnum
  logical, intent(in) :: active(nmrnum)
  if (.not. allocated(meld_active_restraints)) allocate(meld_active_restraints(nmrnum))
  meld_active_restraints = active
end subroutine update_active_restraints

!*******************************************************************************
!
! Subroutine:  get_active_restraints
!
! Description: Return the current logical array of active restraints.
!
!*******************************************************************************

subroutine get_active_restraints(active)
  implicit none
  logical, intent(out), allocatable :: active(:)
  if (.not. allocated(meld_active_restraints)) then
    allocate(active(0))
  else
    allocate(active(size(meld_active_restraints)))
    active = meld_active_restraints
  end if
end subroutine get_active_restraints




end module meld_features_mod