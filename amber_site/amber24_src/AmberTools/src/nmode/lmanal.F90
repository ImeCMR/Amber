
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of program lmanal here]
program lmanal
   
   !    Program to compute correlation functions from langevin
   !    modes.
   !                           j. kottalam    22 Nov 88
   !    Modified to accept perturbation input from gas-phase
   !    modes.
   !                           d.a.case       Oct 89
   implicit double precision (a-h,o-z)
   logical bose, first
#  include "lmsizes.h"
#  include "files.h"
   dimension wr(6*ma), wi(6*ma), el(3*ma,mvect), x(3*ma), &
         iclass(6*ma),amass(ma)
   namelist /data/ kup, lup, nvect, ntrun, tf, np, bose, natom
   
   first = .true.
   call lmfil
   call amopen (5, lmdin, 'O', 'F', 'R')
   call amopen (6, lmdout, owrite, 'F', 'W')
   call amopen (54, inpcrd, 'O', 'F', 'R')
   write(6,*)
   write(6,'(a)') &
         '*****************************************************'
   write(6,'(a)') &
         '       Langevin mode analysis'
   write(6,'(a)') &
         '*****************************************************'
   inml = 0
   10 kup = 1
   lup = 2
   nvect = 24
   ntrun = 1
   tf = 2.0
   np = 1000
   bose = .false.
   inml = inml + 1
   read (5,data,end=99)
   write(6,*)
   write(6,'(a14,2i5)') 'ntrun, nvect: ',ntrun,nvect
   write(6,'(a19,3i5,f10.3)') 'kup, lup, np, tf:  ',kup,lup,np,tf
   if (bose) then
      write(6,'(a,f10.3)') 'bose: true'
   else
      write(6,'(a,f10.3)') 'bose: false'
   end if
   write(6,*)
   mn = 3*ma
   
   !  get masses from prmtop file:
   if( ntrun > 6 ) call rdparml( parm, amass )

   if (ntrun <= 3 .or. ntrun > 6) then
      
      !       ----- read minimized structure from file inpcrd
      
      if( first ) call getref (natom, 1, x)
      first = .false.
      if (natom > ma) then
         write (6,*) ' lmanal: Arrays not dimensioned to handle', &
               natom, ' atoms'
         call mexit(6, 1)
      end if
      
      if (nvect > mvect) then
         write (6,*) ' lmanal: Arrays not dimensioned to handle', &
               nvect, ' eigen vectors'
         call mexit(6, 1)
      end if
      
      nat6 = 6*natom
      nat3 = 3*natom
      write(6,'(a,a)') &
            ' Correlation functions will goto file: ',plot
      write(6,*)
      
      if( ntrun <=3 ) then
         call corfl (nat3, nat6, nvect, &
            kup, lup, ntrun, tf, np, bose, x)
      else
         call corfv (nat3, nat6, nvect, &
            kup, lup, ntrun, tf, np, bose, x, amass)
      endif
   else if (ntrun >= 4) then
      call distr (mn, nat6, wr, wi, ntrun)
   end if
   goto 10
   
   99 call mexit(0, 0)
end program lmanal 
