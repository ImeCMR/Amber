!<compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2010-2012 by 
! Andriy Kovalenko, Tyler Luchko and David A. Case.
! 
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
! 
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
! 
! You should have received a copy of the GNU General Public License in the
! ../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
! 
! Users of the 3D-RISM capability found here are requested to acknowledge
! use of the software in reports and publications.  Such acknowledgement
! should include the following citations:
! 
! 1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
! ibid. 112:10391-10417 (2000).   
! 
! 2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
! F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
! 
! 3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
! and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"

module rism3d_mrc
  use ISO_FORTRAN_ENV
  use safemem
  use rism_report_c
  implicit none

contains

  !> Write volumetric data to a file in MRC 2014 format.  see
  !!
  !! https://www.ccpem.ac.uk/mrc_format/mrc2014.php
  !!
  !! When writing in parallel, each process must call this function
  !! with its local data. Data transfer is handled internally.  We
  !! assume decomposition in the z-axis.
  !!
  !! MRC format is also compatible with the older CCP4
  !! format. However, some readers may ignore MRC-specific information
  !! such as the origin offset.
  !!
  !! @param[in] file File name to write to.
  !! @param[in] data Data to write in a n(1)*n(2)*n(3) linear array.
  !! @param[in] grid Grid object.
  !! @param[in] solute Solute object.
  !! @param[in] o_rank (optional) MPI process rank.
  !! @param[in] o_nproc (optional) MPI number of processes.
  !! @param[in] o_comm (optional) MPI communicator.
  subroutine rism3d_mrc_map_write (file, data, grid, solute, o_rank, o_nproc, o_comm)
    use iso_fortran_env
    use constants, only: PI
    use rism_util, only : freeUnit, rmExPrec
    use rism3d_grid_c
    use rism3d_solute_c
    implicit none
#if defined(MPI)
    include 'mpif.h'
    integer status(MPI_STATUS_SIZE)
#endif /*defined(MPI)*/
    character(len=*), intent(in) :: file
    type(rism3d_grid), intent(in) :: grid
    _REAL_, target, intent(in) :: data(grid%localDimsR(1), grid%localDimsR(2), grid%localDimsR(3)) !, centerOfMass(3)
    type(rism3d_solute), intent(in) :: solute
    integer, optional :: o_rank, o_nproc, o_comm
    
    _REAL_ :: nCellVec(1,3)
    
    integer :: rank = 0, nproc = 1, comm = 0
    integer :: i,j,k, irank, err
    integer(kind=int64) :: count
    integer, parameter :: dataperline = 3
    integer :: unit, iostat

#ifdef MPI
    _REAL_, pointer :: wrk_data(:, :, :) => NULL()
#endif /*MPI*/

    integer :: id
    _REAL_ :: minValue, maxValue, meanValue, rmsd!, totalValue
    logical, parameter :: bigEndian = ichar(transfer(1,'a')) == 0
    ! Up to 80-character long label describing file origin.
    character(len=*), parameter :: amberLabel = 'Amber 3D-RISM MRC map volumetric data.'
    _REAL_ :: temp(1, 3)
    
#ifdef RISM_DEBUG
    write(0, *) "writeMRC", rank
    call flush(6)
#endif /*RISM_DEBUG*/
    
    unit = freeUnit()
    if (present(o_rank)) rank = o_rank
    if (present(o_nproc)) nproc = o_nproc
    if (present(o_comm)) comm = o_comm
    if (rank == 0) then
       ! Unfortunately gfortran does not support form='BINARY' for
       ! compiler-independent binary output, but that probably won't
       ! be an issue for mrc readers.
       ! if (gfortran) then
       open(unit=unit, file=file, iostat=iostat, access='stream', status='replace', form='unformatted')
       ! else ! Intel Fortran
       !     open(unit=unit, file=file, iostat=iostat, access='sequential', status='replace', form='binary')
       ! end if
       if (iostat /= 0) then
          call rism_report_error("opening "//trim(file))
       end if
    end if
    ! Minimum, maximum, and mean density values.
#if defined(MPI)
    call MPI_REDUCE(minval(data), minValue, 1, MPI_DOUBLE, MPI_MIN, 0, comm, err)
    call MPI_REDUCE(maxval(data), maxValue, 1, MPI_DOUBLE, MPI_MAX, 0, comm, err)       
    call MPI_ALLREDUCE(sum(data), meanValue, 1, MPI_DOUBLE, MPI_SUM, comm, err)       
    call MPI_ALLREDUCE(size(data), count, 1, MPI_INTEGER8, MPI_SUM, comm, err)
    meanValue = meanValue/count
    call MPI_REDUCE(sum((meanValue - data)**2), rmsd, 1, MPI_DOUBLE, MPI_SUM, 0, comm, err)
    rmsd = sqrt(rmsd/count)
#else
    minValue = minval(data)
    maxValue = maxval(data)
    meanValue = sum(data) / size(data)
    rmsd = sqrt(sum((meanValue - data)**2) / size(data))
#endif /*defined(MPI)*/
    if (rank == 0) then
       ! Write header.

       ! Number of columns, rows, and sections (fastest to slowest changing).
       ! NX, NY, NZ
       write(unit) int(grid%globalDimsR, int32)

       ! Since values are stored as reals, mode == 2.
       ! MODE
       write(unit) int(2, int32)

       ! There is no offset for column, row, or section.
       ! NXSTART, NYSTART, NZSTART
       write(unit) int((/ 0, 0, 0 /), int32)
       
       ! Number of intervals along X, Y, Z.
       ! MX, MY, MZ
       write(unit) int(grid%globalDimsR, int32)

       ! Cell dimensions (Angstroms).
       ! CELLA
       write(unit) real(grid%boxLength, real32)

       ! Cell angles (degrees).
       ! CELLB
       write(unit) real(grid%unitCellAngles * 180 / PI, real32)

       ! Map column, rows, sects to X, Y, Z (1, 2, 3).
       ! MAPC, MAPR, MAPs
       write(unit) int((/ 1, 2, 3 /), int32)

       ! rmsd = sqrt(rmsd / size(data))
       ! meanValue = totalValue / size(data)
       ! DMIN, DMAX, DMEAN
       write(unit) real(minValue, 4), real(maxValue, 4), real(meanValue, real32)

       ! Space group number.  We assume P 1.
       ! ISPG
       write(unit) int(1, int32)

       ! Number of bytes used for storing symmetry operators.
       ! In our case, none.
       ! NSYMBT
       write(unit) int(0, int32)

       ! extra space used for anything - 0 by default
       ! The format is a bit confusing here.  It indicates 25 words
       ! but it is the third and fourth are EXTTYP AND NVERSION
       ! EXTRA
       do i=1, 2
          write(unit) int(i, int32)
       end do

       ! code for the type of extended header. One of CCP4, MCRO, SERI, AGAR, FEI1, HDF5
       ! We don't use the extended header so it shouldn't matter
       ! EXTTYP
       write(unit) 'MRCO'

       ! version of the MRC format
       ! The version of the MRC format that the file adheres to, specified as a 32-bit integer and calculated as:
       ! Year * 10 + version within the year (base 0)
       ! For the current format change, the value would be 20140. 
       ! NVERSION
       write(unit) int(20140, int32)
       do i=1, 21
          write(unit) int(i, int32)
       end do

       ! The rest of EXTRA, words 29-49
       
       ! phase origin (pixels) or origin of subvolume (A)
       ! This should be the same origin as for DX files
       ! ORIGIN
       if ( .not. grid%periodic ) then
          ! For non-periodic solutes, the solute center of mass (or
          ! center of geometry) was translated to the center of a grid
          ! with an origin at 0,0,0.  When writing out the grid,
          ! offset the origin so the original CoM is in the middle of
          ! the grid.
          write(unit) real(-solute%translation, real32)
       else
           ! For periodic solutes w/ centering=0, there is no translation, but the
          ! best view should have the CoM of the solute somewhere in
          ! the grid. Translate the origin to choose the periodic
          ! image that contains the CoM.
          if ( sum(abs(solute%translation)) == 0d0 ) then
             ! 1. Transform the CoM to fractional coordinates of the unit cell vectors (i.e., reciprocal space).
             nCellVec = transpose( matmul( grid%unitCellVectorsK,  reshape(solute%centerOfMass - solute%translation, [3,1]) ) )
             ! 2. Round down to the nearest whole unit cell vector
             nCellVec = floor(nCellVec)
             ! 3. Set the origin to whatever multiple of each unit cell vectors
             temp = matmul(nCellVec, grid%unitCellVectorsR)
          else
             temp = 0d0
          endif
          write(unit) real( temp(1,:) - solute%translation, real32)
       end if
       
       ! Character string 'MAP ' to identify file type.
       ! MAP
       write(unit) 'MAP '

       ! Machine stamp indicating endianness.
       ! MACHST
       if (bigEndian) then
          ! 0x11 0x11 0x00 0x00
          write(unit) int(z'11110000', int32)
       else
          ! 0x44 0x44 0x00 0x00 but 0x44 0x41 0x00 0x00 (z'00004144') is also ok
          write(unit) int(z'00004444', int32)
       end if
       
       ! RMS deviation of map from mean density.
       ! RMS
       write(unit) real(rmsd, real32)

       ! Number of labels being used.
       ! NLABL
       write(unit) int(1, int32)

       ! Ten 80-character labels.
       ! LABEL
       write(unit) amberLabel
       do id = 1, (9 * 80) + (80 - len(amberLabel))
          write(unit) int(0, int8)
       end do

       ! Symmetry records would go here, but we are not using any.
       
       ! Write volumetric data. This is in column-major format, so
       ! rank-0 always writes out first.
       write(unit) real(data, real32)
    end if
#if defined(MPI)
    ! only rank-0 needs temp data
    if(rank == 0) then
       wrk_data => safemem_realloc(wrk_data, ubound(data, 1), ubound(data, 2), ubound(data, 3), .false.)
    end if
    ! proceed through each process in order and transfer data to
    ! rank-0 for output.  This walks through the z-index.
    do i = 1, nproc-1
       if (rank == 0) then
          call mpi_recv(wrk_data, size(wrk_data), mpi_double, i, 0, comm, status, err)
          write(unit) real(wrk_data,real32)
       elseif(rank == i) then
          call mpi_send(data, size(data), mpi_double, 0, 0, comm, err)
       end if
    end do
    if (safemem_dealloc(wrk_data) /= 0) call rism_report_error("MRC_MAP_WRITE: wrk_data deallocation failed")
#endif /*defined(MPI)*/
  end subroutine rism3d_mrc_map_write
  
end module rism3d_mrc
