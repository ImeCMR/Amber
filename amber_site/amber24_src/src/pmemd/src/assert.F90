#include "../include/assert.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assertion failure reporter
subroutine Aass( condition, file, line )

   implicit none
   character(*) condition
   character(*) file
   integer      line

   write( 6,'(5A,I6,A)') 'ASSERTion ''', condition, ''' failed in ' &
         , file, ' at line ', line, '.'
   call exit(1)
   return
end subroutine Aass 
