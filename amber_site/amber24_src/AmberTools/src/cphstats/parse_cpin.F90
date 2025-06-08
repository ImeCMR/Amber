
#include "constants.h"

subroutine parse_limits(cpein_name, ierr, TITR_RES_C, TITR_STATES_C, ATOM_CHRG_C, MAX_H_COUNT)

   implicit none

   ! File unit
   integer, parameter  :: CPEIN_UNIT = 10

   ! The cpin, cein or cpein file name
   character(len=FN_LEN), intent(in) :: cpein_name

   integer, intent(out) :: ierr, TITR_RES_C, TITR_STATES_C, ATOM_CHRG_C, MAX_H_COUNT

   integer :: ntres, ntstates, natchrg, maxh, ifind

   namelist /cnstphe_limits/ ntres, ntstates, natchrg, maxh

   ! Set default values
   ntres = 50
   ntstates = 200
   natchrg = 1000
   maxh = 4

   ! Open the input file and read the cnstphe_limits namelist
   open(unit=CPEIN_UNIT, file=cpein_name, status='OLD', iostat=ierr)
   if (ierr .ne. 0) then
      write(0, '(a)') 'Failed opening the file at parse_limits'
      return
   end if

   ! Checking if namelist exists
   call nmlsrc('cnstphe_limits', CPEIN_UNIT, ifind)
   if (ifind .ne. 0) then        ! Namelist found. Read it:
      ! Read the namelist, bailing on error
      read(CPEIN_UNIT, nml=cnstphe_limits, iostat=ierr)
      if (ierr .ne. 0) then
         write(0, '(a)') 'Failed reading the file at parse_limits'
         return
      end if
   end if

   ! Close unit
   close (CPEIN_UNIT)

   ! Setting limits
   TITR_RES_C = ntres
   TITR_STATES_C = ntstates
   ATOM_CHRG_C = natchrg
   MAX_H_COUNT = maxh

end subroutine parse_limits

#ifdef REDOX
subroutine parse_cein(trescnt, eleccnt, stateinf, resname, cein_name, is_cpein, ierr, TITR_RES_C, TITR_STATES_C, ATOM_CHRG_C)
#else
subroutine parse_cpin(trescnt, protcnt, stateinf, resname, cpin_name, is_cpein, ierr, TITR_RES_C, TITR_STATES_C, ATOM_CHRG_C)
#endif

   implicit none

   ! The stateinf struct
#ifdef REDOX
   type :: const_e_info
#else
   type :: const_ph_info
#endif
      sequence
      integer :: num_states
      integer :: first_atom
      integer :: num_atoms
      integer :: first_state
      integer :: first_charge
#ifdef REDOX
   end type const_e_info
#else
   end type const_ph_info
#endif

   ! The namelist variables

   integer, intent(in) :: TITR_RES_C, TITR_STATES_C, ATOM_CHRG_C

   integer             :: trescnt
   integer             :: is_cpein
   integer             :: protcnt(0:TITR_STATES_C-1)
   integer             :: eleccnt(0:TITR_STATES_C-1)
   integer             :: resstate(0:TITR_RES_C-1)
   integer             :: ierr
   double precision    :: pka_corr(0:TITR_STATES_C-1)
   double precision    :: eo_corr(0:TITR_STATES_C-1)

#ifdef REDOX
   integer             :: cefirst_sol
   integer             :: ce_igb
   double precision    :: ce_intdiel
#else
   integer             :: cphfirst_sol
   integer             :: cph_igb
   double precision    :: cph_intdiel
#endif
   integer             :: cphefirst_sol
   integer             :: cphe_igb
   double precision    :: cphe_intdiel

   double precision    :: statene(0:TITR_STATES_C-1)
   double precision    :: chrgdat(0:ATOM_CHRG_C-1)

   character(len=40)   :: resname(0:TITR_RES_C)

#ifdef REDOX
   type(const_e_info) :: stateinf(0:TITR_RES_C-1)
   type(const_e_info) :: null_cnste_info = const_e_info(0,0,0,0,0)
#else
   type(const_ph_info) :: stateinf(0:TITR_RES_C-1)
   type(const_ph_info) :: null_cnstph_info = const_ph_info(0,0,0,0,0)
#endif

   ! Is our cpin file read yet?

   logical             :: is_read = .false.

#ifdef REDOX
   ! File unit
   integer, parameter  :: CEIN_UNIT = 10

   ! The cein name
   character(len=FN_LEN), intent(in) :: cein_name
#else
   ! File unit
   integer, parameter  :: CPIN_UNIT = 10

   ! The cpin name
   character(len=FN_LEN), intent(in) :: cpin_name
#endif

   ! The public functions

   ! We read it as a namelist
   namelist /cnstphe/ stateinf, resstate, protcnt, eleccnt, chrgdat, statene, &
                      pka_corr, eo_corr, trescnt, resname, cphefirst_sol, &
                      cphe_igb, cphe_intdiel

#ifdef REDOX
   ! We read it as a namelist
   namelist /cnste/ stateinf, resstate, eleccnt, chrgdat, statene, eo_corr, &
                     trescnt, resname, cefirst_sol, ce_igb, ce_intdiel

   ! Initialize the namelist variables
   stateinf(:) = null_cnste_info
   cefirst_sol = 0
   ce_igb = 0
   ce_intdiel = 0.d0
#else
   ! We read it as a namelist
   namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene, pka_corr, &
                     trescnt, resname, cphfirst_sol, cph_igb, cph_intdiel

   ! Initialize the namelist variables
   stateinf(:) = null_cnstph_info
   cphfirst_sol = 0
   cph_igb = 0
   cph_intdiel = 0.d0
#endif
   trescnt = 0
   resstate(:) = 0
   protcnt(:) = 0
   eleccnt(:) = 0
   pka_corr(:) = 1000.d0
   eo_corr(:) = 0.d0
   chrgdat(:) = 0.d0
   statene(:) = 0.d0
   resname(:) = ' '
   cphefirst_sol = 0
   cphe_igb = 0
   cphe_intdiel = 0.d0
   ierr = 0

#ifdef REDOX
   ! Open the unit, bailing on error
   open(unit=CEIN_UNIT, file=cein_name, status='OLD', iostat=ierr)
   if (ierr .ne. 0) then
      write(0, '(a)') 'Failed opening the file'
      return
   end if

   ! Read the namelist, bailing on error
   if (is_cpein .eq. 1) then
     read(CEIN_UNIT, nml=cnstphe, iostat=ierr)
   else
     read(CEIN_UNIT, nml=cnste, iostat=ierr)
   end if
   if (ierr .ne. 0) then
      write(0, '(a)') 'Failed reading the file'
      return
   end if
#else
   ! Open the unit, bailing on error
   open(unit=CPIN_UNIT, file=cpin_name, status='OLD', iostat=ierr)
   if (ierr .ne. 0) then
      write(0, '(a)') 'Failed opening the file'
      return
   end if

   ! Read the namelist, bailing on error
   if (is_cpein .eq. 1) then
     read(CPIN_UNIT, nml=cnstphe, iostat=ierr)
   else
     read(CPIN_UNIT, nml=cnstph, iostat=ierr)
   end if
   if (ierr .ne. 0) then
      write(0, '(a)') 'Failed reading the file'
      return
   end if
#endif

   ! If we got this far, then our file is read
   is_read = .true.

   return

#ifdef REDOX
end subroutine parse_cein
#else
end subroutine parse_cpin
#endif

subroutine nmlsrc(name,iun,ifind)
   implicit none
   integer:: i, ifind, ilen, iun, j
   ! Subroutine NaMeList SeaRCh

   ! Taken from $AMBEHROME/AmberTools/src/sander/nmlsrc.F90

   character(len=*) name
   character(len=len(name)) txt2
   character(len=8192) aline
   character(len=:), allocatable :: txt1
   ilen = len(name)

   do i = 1,999
      read(iun,1000,end=500,err=500) aline
      do j = 1,80
         if (aline(j:j) /= ' ') then
            if (aline(j:j) == '&' .or. aline(j:j) == '$') then
               if (80-j+1 >= ilen) then
                  ! This includes a trailing space to match the word length.
                  txt1 = aline(j+1:j+ilen+1)
                  call Upcase(aline(j+1:j+ilen+1),txt1)
                  call Upcase(name,txt2)
                  if (txt1 == txt2) goto 200
               else
                  exit
               end if
            else
               exit
            end if
         end if
      end do
   end do

   ! Namelist flag not found:
   500 ifind = 0
   rewind(iun)
   return

   ! Namelist flag found:
   200 ifind = 1
   backspace(iun)
   return

   ! Format statements:
   1000 format(a)
end subroutine nmlsrc

subroutine nmlsrc_from_file(iname,lenname,cpein_name,ifind,ierr)
   implicit none
   integer:: ifind, ierr, lenname
   character(len=FN_LEN) iname
   character(len=lenname) name

   ! File unit
   integer, parameter  :: CPEIN_UNIT = 10

   ! The cpin, cein or cpein file name
   character(len=FN_LEN), intent(in) :: cpein_name

   ! Converting variable
   name = iname

   ! Open the input file and read the cnstphe_limits namelist
   open(unit=CPEIN_UNIT, file=cpein_name, status='OLD', iostat=ierr)
   if (ierr .ne. 0) then
      write(0, '(a)') 'Failed opening the file at nmlsrc_from_file'
      return
   end if

   ! Checking if namelist exists
   call nmlsrc(name, CPEIN_UNIT, ifind)

   ! Close unit
   close (CPEIN_UNIT)
   
end subroutine nmlsrc_from_file

subroutine Upcase(string,upper)

   implicit none

   character(len=*), intent(in) :: string
   character(len=len(string)), intent(out) :: upper
   integer :: i, ic

   do i = 1, len(string)
      ic = iachar(string(i:i))
      if ( ic>96 .and. ic<123 ) then
         ic = ic - 32
      end if
      upper(i:i) = achar(ic)
   end do

end subroutine Upcase