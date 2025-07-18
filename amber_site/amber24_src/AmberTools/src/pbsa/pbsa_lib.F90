! <compile=optimized>
#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ A collection of useful (but difficult-to-place) subroutines.
module pbsa_lib

   implicit none

contains

! These subroutines should have few-to-no dependencies, or we will quickly
! reach a point where we get cyclic dependencies that kills compilation

#if defined SANDER || defined LIBPBSA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Makes a string completely upper-case
subroutine upper(instring)

  implicit none

! Passed variables

  character(*), intent(in out) :: instring

! Local variables

  integer :: i ! counter

  do i = 1, len_trim(instring)

    select case (instring(i:i))
      case('a')
        instring(i:i) = 'A'
      case('b')
        instring(i:i) = 'B'
      case('c')
        instring(i:i) = 'C'
      case('d')
        instring(i:i) = 'D'
      case('e')
        instring(i:i) = 'E'
      case('f')
        instring(i:i) = 'F'
      case('g')
        instring(i:i) = 'G'
      case('h')
        instring(i:i) = 'H'
      case('i')
        instring(i:i) = 'I'
      case('j')
        instring(i:i) = 'J'
      case('k')
        instring(i:i) = 'K'
      case('l')
        instring(i:i) = 'L'
      case('m')
        instring(i:i) = 'M'
      case('n')
        instring(i:i) = 'N'
      case('o')
        instring(i:i) = 'O'
      case('p')
        instring(i:i) = 'P'
      case('q')
        instring(i:i) = 'Q'
      case('r')
        instring(i:i) = 'R'
      case('s')
        instring(i:i) = 'S'
      case('t')
        instring(i:i) = 'T'
      case('u')
        instring(i:i) = 'U'
      case('v')
        instring(i:i) = 'V'
      case('w')
        instring(i:i) = 'W'
      case('x')
        instring(i:i) = 'X'
      case('y')
        instring(i:i) = 'Y'
      case('z')
        instring(i:i) = 'Z'
    end select
  end do

  return

end subroutine upper

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Strips leading whitespace
subroutine strip(instring)

  implicit none

! Passed variables

  character(*), intent(in out) :: instring

! Local variables

  integer :: i, begin

  do i = 1, len_trim(instring)

    if (instring(i:i) .gt. ' ') then
      begin = i
      exit
    end if

  end do

  instring = instring(begin:len_trim(instring))

  return

end subroutine strip
#endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Determines how many words are in a string
subroutine get_num_tokens(string, token_num)

  implicit none

! Passed arguments

  character(*), intent(in) :: string

  integer, intent(out)     :: token_num

! Local variables

  integer :: string_loc  ! our location in the string
  integer :: iend        ! last non-whitespace character location

  string_loc = 1
  iend = len_trim(string)
  token_num = 0

  do while (string_loc .le. iend)

    if ( string(string_loc:string_loc) .le. ' ' ) then
      string_loc = string_loc + 1
    else

      do while ( string(string_loc:string_loc) .gt. ' ' )
        string_loc = string_loc + 1
      end do

      token_num = token_num + 1
    end if
  end do

end subroutine get_num_tokens

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Gets a specific word from a string
subroutine get_token(string, num, token)

  implicit none

! Passed arguments

  character(*), intent(in)  :: string  ! The string to parse
  character(*), intent(out) :: token   ! The token to return

  integer, intent(in)       :: num     ! Which token to return

! Local variables

  integer   :: num_tokens
  integer   :: istart
  integer   :: iend
  integer   :: string_loc
  integer   :: token_count

  ! Uncomment the below chunk of code for a "safe" get_num_tokens at the
  ! expense of calling get_num_tokens() each time a specific token is
  ! pulled from the string. When it's commented out, token will just be
  ! a blank string upon return

! call get_num_tokens(string, num_tokens)

! if (num .gt. num_tokens)
!   write(mdout, *) ' Error in get_token: Looking for more tokens than &
!                     &there are in string'
!   call mexit(6,1)
! end if

  ! Now get the num'th token

  token_count = 0
  istart = 1
  iend = len_trim(string)
  token = ' '

  do while (istart .le. iend)

    if (string(istart:istart) .le. ' ') then

      istart = istart + 1

    else

      do string_loc = istart, iend
        if ( string(string_loc:string_loc) .le. ' ' ) exit
      end do

      token_count = token_count + 1

      ! If this is the token we want, store it and return
      if ( token_count .eq. num ) then
        token = string(istart:string_loc-1)
        return
      end if

      istart = string_loc ! Move to the next token

    end if

  end do

end subroutine get_token

end module pbsa_lib
