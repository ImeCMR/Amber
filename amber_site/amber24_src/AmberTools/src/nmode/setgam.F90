!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setgam here]
subroutine setgam (natom, a, amass, exp, eta, gamma, x, &
      ioseen, hrmax,igraph)
   
   !     ----- to set up the friction matrix gamma
   
   implicit double precision (a-h,o-z)
   parameter (maxat=3000)
   dimension expr(maxat)
   character(len=4) label
   dimension a(6*natom,6*natom), amass(natom), exp(natom), gamma(*)
   dimension x(3*natom)
   character(len=4) igraph(natom)
#  include "files.h"
   namelist /expos/ expr
   nr3 = 3*natom
   nr6 = 6*natom
   
   !     ----- set up hydrodynamic radii
   
   if (natom > maxat) then
      write(6,*) 'MAXAT is too small in subroutine setgam!'
      call mexit(6, 1)
   end if
   do i = 1, natom
      expr(i) = 0.0
   end do
   
   !     err=40 prevents amopen() from being used here
   open (10, file = expfil, status = 'old', err=40)
   read (10, expos)
   close (10)
   do i=1,natom
      exp(i) = expr(i)
   end do
   
   big = exp(1)
   do i = 2, natom
      if (exp(i) > big) big = exp(i)
   end do
   
   factor = hrmax / sqrt(big)
   do i = 1, natom
      exp(i) = sqrt(exp(i)) * factor
      write(label,'(a4)') igraph(i)
      if (label(1:1) == 'H') exp(i) = min(0.2d0,exp(i))
   end do
   goto 50
   
   40 continue
   write (6,*)
   write (6,'(t1,a)') ' expfile is not present. Hydrodynamic radius'
   write (6,'(t1,a,f6.4,a)') &
         '   of ',hrmax,' angs. is used for all atoms'
   write (6,'(t1,a)') '   except hydrogen, where max. value is 0.2'
   write (6,*)
   do i = 1, natom
      exp(i) = hrmax
      write(label,'(a4)') igraph(i)
      if (label(1:1) == 'H') exp(i) = min(0.2d0,exp(i))
   end do
   
   50 continue
   
   !     ----- set up friction matrix gamma
   
   if (ioseen == 0) then
      
      !       ----- diagonal gamma
      
      sxpita = 6.0 * 3.14159 * eta
      ij = 0
      do j = 1, 3 * natom
         do i = 1, j-1
            ij = ij + 1
            gamma(ij) =  0.0
         end do
         ij = ij + 1
         jat = (j-1)/3 + 1
         gamma(ij) = sxpita * exp(jat) / amass(jat)
      end do
      
   else
      
      !       ----- gamma with hydrodynamic interaction
      
      call oseen (ioseen, natom, nr3, nr6, a, eta, exp, x, &
            gamma, amass)
      do j = 1, nr6
         do i = 1, nr6
            a(i,j) = 0.0
         end do
      end do
      
   end if
   
   !     ----- copy gamma into matrix a
   
   do j = 1, nr3
      do i = 1, j
         index = j*(j-1)/2 + i
         a(i+nr3,j+nr3) = - gamma(index)
         a(j+nr3,i+nr3) = - gamma(index)
      end do
   end do
   
   return
end subroutine setgam 
