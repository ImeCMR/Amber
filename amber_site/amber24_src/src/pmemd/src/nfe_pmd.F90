! Pre-processing
#ifndef NFE_UTILS_H
#define NFE_UTILS_H

#ifndef NFE_DISABLE_ASSERT
#  define nfe_assert(stmt) if (.not.(stmt)) call afailed(__FILE__, __LINE__)
#  define nfe_assert_not_reached() call afailed(__FILE__, __LINE__)
#  define NFE_PURE_EXCEPT_ASSERT
#  define NFE_USE_AFAILED use nfe_lib_mod, only : afailed
#else
#  define nfe_assert(s) 
#  define nfe_assert_not_reached() 
#  define NFE_PURE_EXCEPT_ASSERT pure
#  define NFE_USE_AFAILED
#endif /* NFE_DISABLE_ASSERT */

#define NFE_OUT_OF_MEMORY call out_of_memory(__FILE__, __LINE__)

#ifdef MPI
#  define NFE_MASTER_ONLY_BEGIN if (mytaskid.eq.0) then
#  define NFE_MASTER_ONLY_END end if
#else
#  define NFE_MASTER_ONLY_BEGIN
#  define NFE_MASTER_ONLY_END
#endif /* MPI */

#define NFE_ERROR   ' ** NFE-Error ** : '
#define NFE_WARNING ' ** NFE-Warning ** : '
#define NFE_INFO    ' NFE : '

#endif /* NFE_UTILS_H */
!-------------------------------------------------------------------------

module nfe_pmd_mod

use nfe_lib_mod, only : SL => STRING_LENGTH, PMD_OUTPUT_UNIT, PMD_CV_UNIT
use nfe_colvar_mod

implicit none

private

! old REMD subroutines, not used now ------
#ifdef MPI
public :: on_delta
public :: on_exchange
#endif /* MPI */
public :: on_multipmemd_exit
! -----------------------------------------

public :: on_pmemd_init
public :: on_pmemd_exit

public :: on_force

!- - - - - - - - - - - - - - - - P R I V A T E - - - - - - - - - - - - - - - -

integer, private, parameter :: pmd_UNIT = PMD_OUTPUT_UNIT
integer, private, parameter :: CV_UNIT = PMD_CV_UNIT

character(*), private, parameter :: DEFAULT_OUTPUT_FILE = 'nfe-pmd'
character(*), private, parameter :: DEFAULT_CV_FILE = 'nfe-pmd-cv'

integer, private, parameter :: DEFAULT_OUTPUT_FREQ = 50

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer, private, save :: ncolvars = 0 ! .gt.0 means "active"
type(colvar_t), private, allocatable, save :: cv(:)

double precision, private, pointer, save :: anchor(:) => null() ! master
double precision, private, pointer, save :: a_position(:) => null() ! master
double precision, private, pointer, save :: a_strength(:) => null() ! master

double precision, private, allocatable, save :: cv_inst(:)
double precision, private, allocatable, save :: f_cv(:)

character(SL), private, save :: output_file = DEFAULT_OUTPUT_FILE
character(SL), private, save :: output_fmt
character(SL), private, save :: cv_file = DEFAULT_CV_FILE

integer, private, save :: output_freq = DEFAULT_OUTPUT_FREQ
integer, private, save :: mdstep ! = runmd.f::nstep + 1 (not zeroed on exchange)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef MPI
double precision, private, pointer, save :: o_anchor(:) => null() ! master
#endif /* MPI */

namelist / pmd /     output_file, output_freq, cv_file

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

contains

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef MPI
subroutine on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)

   use nfe_colvar_mod, only : colvar_difference
   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   integer, intent(in) :: o_masterrank
   logical, intent(in) :: need_U_xx

   double precision, intent(inout) :: U_mm, U_mo, U_om, U_oo

   double precision :: U_o(2), U_m(2)

   integer :: n, error

   if (ncolvars.eq.0) &
      return

   nfe_assert(multipmemd_rem().ne.0)
   nfe_assert(mytaskid.eq.0) ! master
   nfe_assert(pmemd_master_comm.ne.mpi_comm_null)

   ! exchange cv_inst(:) with the partner [store partner values in f_cv]
   call mpi_sendrecv &
      (cv_inst, ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
          f_cv, ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       pmemd_master_comm, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)

   ! evaluate 'my' values
   U_m(1) = ZERO ! U_mm = U_m(x_m)
   U_m(2) = ZERO ! U_mo = U_m(x_o)

   do n = 1, ncolvars
     if (cv_inst(n).le.a_position(4*n-3)) then
        U_m(1) = a_strength(2*n-1)*colvar_difference(cv(n),a_position(4*n-3),a_position(4*n-2))*cv_inst(n)
     else if (cv_inst(n).le.a_position(4*n-2)) then
        U_m(1) = a_strength(2*n-1)*colvar_difference(cv(n),cv_inst(n),a_position(4*n-2))**2/2
     else if (cv_inst(n).le.a_position(4*n-1)) then
        U_m(1) = 0.0
     else if (cv_inst(n).le.a_position(4*n)) then
        U_m(1) = a_strength(2*n)*colvar_difference(cv(n),cv_inst(n),a_position(4*n-1))**2/2
     else
        U_m(1) = a_strength(2*n)*colvar_difference(cv(n), a_position(4*n), a_position(4*n-1))*cv_inst(n)
     end if

     if (f_cv(n).le.a_position(4*n-3)) then
        U_m(2) = a_strength(2*n-1)*colvar_difference(cv(n),a_position(4*n-3),a_position(4*n-2))*f_cv(n)
     else if (f_cv(n).le.a_position(4*n-2)) then
        U_m(2) = a_strength(2*n-1)*colvar_difference(cv(n),f_cv(n),a_position(4*n-2))**2/2
     else if (f_cv(n).le.a_position(4*n-1)) then
        U_m(2) = 0.0
     else if (f_cv(n).le.a_position(4*n)) then
        U_m(2) = a_strength(2*n)*colvar_difference(cv(n),f_cv(n),a_position(4*n-1))**2/2
     else
        U_m(2) = a_strength(2*n)*colvar_difference(cv(n), a_position(4*n),a_position(4*n-1))*f_cv(n)
     end if
   end do

   ! get partner's U_m? (i.e., U_o? in this replica)
   call mpi_sendrecv &
      (U_m, 2, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       U_o, 2, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       pmemd_master_comm, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)

   if (need_U_xx) then
      U_mm = U_mm + U_m(1)
      U_mo = U_mo + U_m(2)
      U_om = U_om + U_o(2)
      U_oo = U_oo + U_o(1)
   end if

end subroutine on_delta

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_exchange(o_masterrank)

#  ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
   use parallel_dat_mod
#  endif /* NFE_DISABLE_ASSERT */

   implicit none

   integer, intent(in) :: o_masterrank

   character(SL) :: o_output_file
   integer :: o_output_freq, error

   if (ncolvars.eq.0) &
      return

   nfe_assert(multipmemd_rem().ne.0)
   nfe_assert(mytaskid.eq.0) ! master
   nfe_assert(pmemd_master_comm.ne.mpi_comm_null) ! master

   ! slow & naive

   call mpi_sendrecv(output_file, SL, MPI_CHARACTER, o_masterrank, 5, &
                   o_output_file, SL, MPI_CHARACTER, o_masterrank, 5, &
                     pmemd_master_comm, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)
   output_file = o_output_file
   call flush_UNIT(pmd_UNIT, output_file)

   call mpi_sendrecv(output_freq, 1, MPI_INTEGER, o_masterrank, 6, &
                   o_output_freq, 1, MPI_INTEGER, o_masterrank, 6, &
                     pmemd_master_comm, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)
   output_freq = o_output_freq

   nfe_assert(associated(anchor))
   nfe_assert(associated(o_anchor))

   call mpi_sendrecv &
      (anchor, 6*ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, 7, &
     o_anchor, 6*ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, 7, &
       pmemd_master_comm, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)

   anchor(1:6*ncolvars) = o_anchor(1:6*ncolvars)

end subroutine on_exchange
#endif /* MPI */

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_multipmemd_exit()

   use nfe_colvar_mod, only : colvar_cleanup
#  ifdef MPI
   use nfe_lib_mod, only : multipmemd_rem
#  endif /* MPI */
   use parallel_dat_mod

   implicit none

   integer :: n

   if (ncolvars.gt.0) then
      do n = 1, ncolvars
         call colvar_cleanup(cv(n))
      end do
      deallocate(cv, f_cv, cv_inst)
      NFE_MASTER_ONLY_BEGIN
      nullify(a_position, a_strength)
      deallocate(anchor)
#     ifdef MPI
      if (multipmemd_rem().ne.0) &
         deallocate(o_anchor)
#     endif /* MPI */
      NFE_MASTER_ONLY_END
   end if

   mdstep = 0
   ncolvars = 0

end subroutine on_multipmemd_exit

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_pmemd_init(mdin_unit, amass)

   use nfe_lib_mod
   use nfe_colvar_mod
   use file_io_mod
   use parallel_dat_mod

   implicit none

   integer, intent(in) :: mdin_unit

   double precision, intent(in) :: amass(*)

   integer :: n, error, ifind, i
   character(80) :: buf

   NFE_MASTER_ONLY_BEGIN

   nfe_assert(ncolvars.eq.0)
   
#  ifdef MPI
   if (multipmemd_numgroup().gt.1) then
      n = master_rank + 1
      write (unit = output_file, fmt = '(a,a,i3.3,a)') &
           DEFAULT_OUTPUT_FILE, '.', n, '.dat'
   else
#  endif /* MPI */
      output_file = DEFAULT_OUTPUT_FILE//'.dat'
#  ifdef MPI
   end if
#  endif /* MPI */

   rewind(mdin_unit)
   call nmlsrc('pmd', mdin_unit, ifind)

   ! no pmd section
   if (ifind.eq.0) goto 1

   if (pmemd_imin().ne.0) &
      call fatal('imin.ne.0 is not supported')

   rewind(mdin_unit)
   read(mdin_unit,nml=pmd,err=666)
   ncolvars = 0

   ! collective variables
   call amopen(CV_UNIT, cv_file, 'O', 'F', 'R')
   
   do
     call nmlsrc('colvar', CV_UNIT, ifind)
     if (ifind.eq.0) exit
     read(CV_UNIT,'(a80)') buf
     ncolvars = ncolvars + 1
   end do

   if (ncolvars.eq.0) &
      call fatal('no variable(s) in the CV file')

   allocate(cv(ncolvars), anchor(6*ncolvars), &
      f_cv(ncolvars), cv_inst(ncolvars), stat = error)
   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   a_position => anchor(0*ncolvars + 1:4*ncolvars)
   a_strength => anchor(4*ncolvars + 1:6*ncolvars)

   n = 1
   rewind(CV_UNIT)

   do while (n.le.ncolvars)  
         call colvar_nlread(CV_UNIT,cv(n))
         a_strength(2*n-1) = anchor_strength(1)
         a_strength(2*n)   = anchor_strength(2)
         a_position(4*n-3) = anchor_position(1)
         a_position(4*n-2) = anchor_position(2)
         a_position(4*n-1) = anchor_position(3)
         a_position(4*n)   = anchor_position(4)

         if ( a_position(4*n-3).gt.a_position(4*n-2)) &
            call fatal('anchor_position 1 cannot be larger than anchor_position 2')

         if ( a_position(4*n-2).gt.a_position(4*n-1)) &
            call fatal('anchor_position 2 cannot be larger than anchor_position 3')

         if ( a_position(4*n-1).gt.a_position(4*n)) &
            call fatal('anchor_position 3 cannot be larger than anchor_position 4')

         if ( a_strength(2*n-1).lt.0 .or. a_strength(2*n).lt.0) &
            call fatal('anchor_strength cannot be negative')
         
         if (colvar_has_refcrd(cv(n))) then
            refcrd_file = trim(refcrd_file)
            refcrd_len = len_trim(refcrd_file)
         end if
         if (colvar_is_quaternion(cv(n))) then
          allocate(cv(n)%q_index, stat = error)
           if (error.ne.0) &
            NFE_OUT_OF_MEMORY
            cv(n)%q_index = q_index
         end if
         if (colvar_has_axis(cv(n))) then
           allocate(cv(n)%axis(3), stat = error)
             if (error.ne.0) &
                NFE_OUT_OF_MEMORY
           i = 1
           do while (i.le.3)
             cv(n)%axis(i) = axis(i)
             i = i + 1
           end do
         end if

         nfe_atm_cnt = nfe_atm_cnt + cv_ni
         do i=1,cv_ni
            call AddToList(nfe_atm_lst,cv_i(i))
         end do

         n = n + 1
   end do

   output_freq = min(output_freq, pmemd_nstlim())
   output_freq = max(1, output_freq)

1  continue

   NFE_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(ncolvars, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(refcrd_len, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(refcrd_file, refcrd_len, MPI_CHARACTER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#  endif /* MPI */

   if (ncolvars.eq.0) &
      return

#  ifdef MPI
   if (mytaskid.ne.0) then
      allocate(cv(ncolvars), f_cv(ncolvars), cv_inst(ncolvars), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY
   end if
#  endif /* MPI */

   do n = 1, ncolvars
      call colvar_bootstrap(cv(n), n, amass)
   end do

   mdstep = 0

   NFE_MASTER_ONLY_BEGIN

   open (unit = pmd_UNIT, file = output_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'failed to open ''', trim(output_file), ''' for writing'
      call terminate()
   end if

   write (unit = pmd_UNIT, fmt = '(a,66(''=''))') '# = NFE%PMD '
   do n = 1, ncolvars
      write (unit = pmd_UNIT, &
         fmt = '(a,'//pfmt(n)//',                    &
               & a,'//pfmt(a_position(4*n-3), 6)//', &
               & a,'//pfmt(a_position(4*n-2), 6)//', &
               & a,'//pfmt(a_position(4*n-1), 6)//', &
               & a,'//pfmt(a_position(4*n),   6)//', &
               & /                                   &
               & a,'//pfmt(a_strength(2*n-1), 6)//', &
               & a,'//pfmt(a_strength(2*n),   6)//', &
               & a)')                                &
         '#   << anchor(', n, ') : position = ', a_position(4*n-3), ', ',a_position(4*n-2),', ',& 
           a_position(4*n-1),', ',a_position(4*n),&
         '#                  strength = ', a_strength(2*n-1),', ', a_strength(2*n), ' >>'
   end do
   write (unit = pmd_UNIT, &
      fmt = '(a,77(''-''),/a,'//pfmt(ncolvars)//',a,/a,77(''=''))') '# ', &
      '# MD time (ps), CV(1:', ncolvars, ')', '# '

   call flush_UNIT(pmd_UNIT, output_file)

   write (unit = output_fmt, fmt = '(a,'//pfmt(ncolvars)//',a)') &
      '(f12.4,', ncolvars, '(1x,f16.8))'

   ! print summary & we'r done

   write (unit = OUT_UNIT, fmt = '(a,a)') NFE_INFO, &
      '~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ P I N N E D  M.D. ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~'

   write (unit = OUT_UNIT, fmt = '(a,/a,a,a)') NFE_INFO, NFE_INFO, &
      'output_file = ', trim(output_file)
   write (unit = OUT_UNIT, &
      fmt = '(a,a,'//pfmt(output_freq)//',a,'//pfmt(output_freq*pmemd_timestep(), 4)//',a)') &
        NFE_INFO, 'output_freq = ', output_freq, ' (', &
        output_freq*pmemd_timestep(), ' ps)'

   write (unit = OUT_UNIT, fmt = '(a)') NFE_INFO

   do n = 1, ncolvars
      write (unit = OUT_UNIT, &
         fmt = '(a,a,'//pfmt(n)//',                  &
               & a,'//pfmt(a_position(4*n-3), 6)//', &
               & a,'//pfmt(a_position(4*n-2), 6)//', &
               & a,'//pfmt(a_position(4*n-1), 6)//', &
               & a, '//pfmt(a_position(4*n),  6)//', &
               & / a,                                &
               & a,'//pfmt(a_strength(2*n-1), 6)//', &
               & a,'//pfmt(a_strength(2*n),   6)//', &
               & a)')                                &
         NFE_INFO, 'CV #', n, ' << anchor : position = ', a_position(4*n-3),', ', a_position(4*n-2),', ',& 
         a_position(4*n-1),', ',a_position(4*n), NFE_INFO,&
         '                  strength = ', a_strength(2*n-1),', ', a_strength(2*n), ' >>'

      call colvar_print(cv(n), OUT_UNIT)
      write (unit = OUT_UNIT, fmt = '(a)') NFE_INFO
      if (colvar_is_quaternion(cv(n))) then
          write (unit = OUT_UNIT, fmt = '(a,a,I3)') NFE_INFO, &
          ' <> q_index = ', cv(n)%q_index
      end if
      if (colvar_has_axis(cv(n))) then
           write (unit = OUT_UNIT, fmt = '(a,a,f8.4,a,f8.4,a,f8.4,a)') NFE_INFO, &
           ' <> axis = [',cv(n)%axis(1),', ', cv(n)%axis(2),', ', cv(n)%axis(3),']'
      end if
   end do

   write (unit = OUT_UNIT, fmt = '(a,a/)') NFE_INFO, &
      '~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~'

   NFE_MASTER_ONLY_END
   return
   
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &pmd namelist!'
    call terminate()
end subroutine on_pmemd_init

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_pmemd_exit()

   use nfe_lib_mod, only : close_UNIT
   use nfe_colvar_mod
   use parallel_dat_mod

   implicit none

   integer :: n

   if (ncolvars.gt.0) then
      do n = 1, ncolvars
         call colvar_cleanup(cv(n))
      end do
      deallocate(cv, f_cv, cv_inst)
      NFE_MASTER_ONLY_BEGIN
      call close_UNIT(pmd_UNIT)
      nullify(a_position, a_strength)
      deallocate(anchor)
      NFE_MASTER_ONLY_END
   end if
end subroutine on_pmemd_exit

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_force(x, f, pot)

   use nfe_colvar_mod
   use nfe_lib_mod
   use parallel_dat_mod
   use gbl_constants_mod, only : PI

   implicit none

   double precision, intent(in) :: x(3,*)
   double precision, intent(inout) :: f(3,*)
   double precision, intent(inout) :: pot

   integer :: n, m 
   double precision :: fix_cv, apmean
   integer, DIMENSION(4) :: cv_q = (/COLVAR_QUATERNION0, COLVAR_QUATERNION1, &
                                     COLVAR_QUATERNION2, COLVAR_QUATERNION3/)
   double precision :: norm4(100), cv_N(4, 100)

#  ifdef MPI
   integer :: error
#  endif /* MPI */

   if (ncolvars.eq.0) &
      return

   do n = 1, ncolvars
      cv_inst(n) = colvar_value(cv(n), x)
   end do

   NFE_MASTER_ONLY_BEGIN
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n))) then
        do m = 1, 4
         if (cv(n)%type == cv_q(m)) then
            cv_N(m, cv(n)%q_index) = cv_inst(n)
         end if
        end do
      end if
   end do
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n))) then
       norm4(cv(n)%q_index) = sqrt(cv_N(1,cv(n)%q_index)**2 &
                                 + cv_N(2,cv(n)%q_index)**2 &
                                 + cv_N(3,cv(n)%q_index)**2 &
                                 + cv_N(4,cv(n)%q_index)**2)
      end if
   end do
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n))) then
         cv_inst(n) = cv_inst(n) / norm4(cv(n)%q_index)
      else
         cv_inst(n) = cv_inst(n)
      end if
   end do

   do n = 1, ncolvars

     fix_cv = cv_inst(n)
     if (colvar_is_periodic(cv(n))) then
        apmean = (a_position(4*n-2) + a_position(4*n-1))*0.5

10      if (fix_cv - apmean .gt. PI) then
            fix_cv = fix_cv - PI - PI
            goto 10
        else if (apmean - fix_cv .gt. PI) then
            fix_cv = fix_cv + PI + PI
            goto 10
        end if
     end if

     if (fix_cv.le.a_position(4*n-3)) then
        f_cv(n) = -a_strength(2*n-1)*colvar_difference(cv(n),a_position(4*n-3),a_position(4*n-2))
        pot = pot + a_strength(2*n-1)* &
              colvar_difference(cv(n),a_position(4*n-3),a_position(4*n-2))*fix_cv + &
              0.5*a_strength(2*n-1)*(a_position(4*n-2)**2-a_position(4*n-3)**2) 
     else if (fix_cv.le.a_position(4*n-2)) then
        f_cv(n) = -a_strength(2*n-1)*colvar_difference(cv(n),cv_inst(n),a_position(4*n-2))
        pot = pot + a_strength(2*n-1)*colvar_difference(cv(n),cv_inst(n),a_position(4*n-2))**2/2
     else if (fix_cv.le.a_position(4*n-1)) then
        f_cv(n) = 0.0
     else if (fix_cv.le.a_position(4*n)) then
        f_cv(n) = -a_strength(2*n)*colvar_difference(cv(n),cv_inst(n),a_position(4*n-1))
        pot = pot + a_strength(2*n)*colvar_difference(cv(n),cv_inst(n),a_position(4*n-1))**2/2
     else
        f_cv(n) = -a_strength(2*n)*colvar_difference(cv(n), a_position(4*n), a_position(4*n-1))
        pot = pot + a_strength(2*n)* &
              colvar_difference(cv(n), a_position(4*n), a_position(4*n-1))*fix_cv + &
              0.5*a_strength(2*n)*(a_position(4*n-1)**2-a_position(4*n)**2)
     end if

   end do
   NFE_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(f_cv, ncolvars, &
      MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(pot, 1, &
      MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#  endif /* MPI */

   ! FIXME: virial
   do n = 1, ncolvars
      call colvar_force(cv(n), x, f_cv(n), f)
   end do

   NFE_MASTER_ONLY_BEGIN

   if (nfe_real_mdstep) then

      if (mod(mdstep, output_freq).eq.0) then
         write (unit = pmd_UNIT, fmt = output_fmt) &
            pmemd_mdtime(), cv_inst(1:ncolvars)
         call flush_UNIT(pmd_UNIT, output_file)
      end if

      mdstep = mdstep + 1

   end if ! real_mdstep
   NFE_MASTER_ONLY_END

end subroutine on_force

end module nfe_pmd_mod
