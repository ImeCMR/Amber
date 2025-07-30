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



end module meld_features_mod