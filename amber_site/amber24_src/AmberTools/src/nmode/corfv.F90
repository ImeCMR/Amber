
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ much like corfl, but for velocity autocorrelation functions
subroutine corfv (nat3, nat6, nvect, kup, lup, ntrun, tf, np, bose, x, &
                  amass )
   
   !     ----- to calculate correlation functions from langevin modes
   
   implicit double precision (a-h,o-z)
   double precision kt
   logical first,bose
   complex lambda, delijc, zero
   
   !     automatic arrays:
   dimension wr(nat6),wi(nat6),el(nat3,nvect),iclass(nat6)
   dimension del(nat3,nat3), x(nat3), amass(*)

#  include "files.h"
   !                  ----CONSQ = hc/2kT in cm, with T=300K
   !                       (use for quantum, Bose statistics)
   parameter (kt = 0.6)
   parameter (consq = 2.39805d-3)
   data first /.true./
   
   zero = (0.0,0.0)
   rzero = 0.0
   nat = nat3/3
   totmass = 0.0
   do i=1,nat
      totmass = totmass + amass(i)
   end do
   write(6,*) 'totmass = ', totmass
   
   !     ----- read langevin modes from file lmode
   
   if (first) call amopen (18, plot, owrite, 'F', 'W')
   first = .false.
   call lmin (nat3, nat6, nvect, wr, wi, el, iclass)
   
   !     ----- set tf and tdel
   
   ti = 0.0
   tf = tf * 20.455
   tdel = (tf-ti)/np
       do nt=1,np+1
         t = ti + (nt-1)*tdel
   
#undef SINGLE_ATOM 
#ifdef SINGLE_ATOM
   i = kup
   j = lup    

#else
   do i=1,nat3
     do j=1,nat3
#endif

   !     ----- for each t, calculate the correlation matrix for vk, vl
   !        as in Eq. 7.1 of lamm and szabo j. chem. phys. 85, 7334 (1986)
   
   del(i,j) = 0.0
   n = 1  ! check this--looks likely to be wrong for velocities
   
   do while (n.le.nvect)
      lambda = cmplx ( wr(n), wi(n) )
      if (bose) then
         argq = sqrt(wr(n)**2 + wi(n)**2)*108.59*consq
         qcorr = argq/tanh(argq)
      end if
#ifdef debug
      if (i == 1 .and. j == 1 .and. wr(n) /= 0.0 .and. t == 0.0) &  
         write(6,'(i4,5f10.4)') n,wi(n)*108.59,1./(wr(n)*20.455), &
         sqrt(wr(n)**2 + wi(n)**2)*108.59,argq,qcorr
#endif
      if (iclass(n) == 5 .or. lambda == zero) then
         ! ---skip this mode
         n = n + 1
      else if (iclass(n) == 2 .or. iclass(n) == 4) then
         ! ---purely real mode
         ! if (real(lambda)*t .lt. -50.) cycle
         delijc = exp(lambda*t) * el(i,n) * el(j,n)
         if (bose) delijc = delijc*qcorr
         del(i,j) = del(i,j) + kt*real(delijc)
         write(6,'(4i4,4e14.5)') i,j, n,iclass(n),lambda,t,del(i,j)
         n = n + 1
      else if (iclass(n) == 3) then
         ! ---purely imaginary mode
         ! if (real(lambda)*t .lt. -50.) cycle
         delijc = exp(lambda*t) &
               * cmplx( rzero, el(i,n)) * cmplx( rzero, el(j,n))
         if (bose) delijc = delijc*qcorr
         del(i,j) = del(i,j) + kt*real(delijc)
         write(6,'(4i4,4e14.5)') i,j, n,iclass(n),lambda,t,del(i,j)
         n = n + 1
      else if (iclass(n) == 1) then
         ! ---complex mode
         ! if (real(lambda)*t .lt. -50.) cycle
         delijc = exp(lambda*t) &
               * cmplx( el(i,n), el(i,n+1) ) &
               * cmplx( el(j,n), el(j,n+1) )
         if (bose) delijc = delijc*qcorr
         del(i,j) = del(i,j) + 2.0 * kt*real(delijc)
         write(6,'(4i4,4e14.5)') i,j, n,iclass(n),lambda,t,del(i,j)
         n = n + 2
      end if  ! (iclass(n) == 5 .or. lambda == zero)

   end do  ! n
   
#ifdef SINGLE_ATOM
   write(18,'(f10.3,2x,f12.5)') t/20.455,del(i,j)
#else
   end do  ! j
   end do  ! i

   ! here we need to combine the info from the i,j correlation functions time t

   corfvx = 0.0
   corfvy = 0.0
   corfvz = 0.0
   corfvxy = 0
   corfvxz = 0
   corfvyz = 0

   if( ntrun == 7 ) then  !  translational motion: atomic velocities

      do i = 1,nat
         i3  = 3*(i-1)
         do j = 1,nat
            j3 = 3*(j-1)
            corfvx = corfvx + amass(i)*amass(j)*del(i3+1,j3+1)
            corfvy = corfvy + amass(i)*amass(j)*del(i3+2,j3+2)
            corfvz = corfvz + amass(i)*amass(j)*del(i3+3,j3+3)
            corfvxy = corfvxy + amass(i)*amass(j)*del(i3+1,j3+2)
            corfvxz = corfvxz + amass(i)*amass(j)*del(i3+1,j3+3)
            corfvyz = corfvyz + amass(i)*amass(j)*del(i3+2,j3+3)
         end do
      end do

   else if (ntrun == 8) then   !   rotational motion: angular velocities

      ! avoid divisions by zero (kludge? or workable?)
      do i=1,nat3
         if( abs(x(i)) < 1.d-5 ) x(i) = 1.d10
      end do

      do i = 1,nat
         i3  = 3*(i-1)
         do j = 1,nat
            j3 = 3*(j-1)
            corfvx = corfvx + amass(i)*amass(j)*  &
               (   del(i3+2,j3+2)/(x(i3+3)*x(j3+3)) &
                 - del(i3+3,j3+2)/(x(i3+2)*x(j3+3)) &
                 - del(i3+2,j3+3)/(x(i3+3)*x(j3+2)) &
                 + del(i3+3,j3+3)/(x(i3+2)*x(j3+2)) )
                 
            corfvy = corfvy + amass(i)*amass(j)*  &
               (   del(i3+1,j3+1)/(x(i3+3)*x(j3+3)) &
                 - del(i3+3,j3+1)/(x(i3+1)*x(j3+3)) &
                 - del(i3+1,j3+3)/(x(i3+3)*x(j3+1)) &
                 + del(i3+3,j3+3)/(x(i3+1)*x(j3+1)) )

            corfvz = corfvz + amass(i)*amass(j)*  &
               (   del(i3+1,j3+1)/(x(i3+2)*x(j3+2)) &
                 - del(i3+1,j3+2)/(x(i3+2)*x(j3+1)) &
                 - del(i3+2,j3+1)/(x(i3+1)*x(j3+2)) &
                 + del(i3+2,j3+2)/(x(i3+1)*x(j3+1)) )

            corfvxy = corfvxy + amass(i)*amass(j)*  &
               (   del(i3+2,j3+1)/(x(i3+3)*x(j3+3)) &
                 - del(i3+3,j3+1)/(x(i3+2)*x(j3+3)) &
                 - del(i3+2,j3+3)/(x(i3+3)*x(j3+1)) &
                 + del(i3+3,j3+3)/(x(i3+2)*x(j3+1)) )

            corfvxz = corfvxz + amass(i)*amass(j)*  &
               (   del(i3+2,j3+1)/(x(i3+3)*x(j3+2)) &
                 - del(i3+3,j3+2)/(x(i3+2)*x(j3+1)) &
                 - del(i3+2,j3+1)/(x(i3+3)*x(j3+2)) &
                 + del(i3+3,j3+2)/(x(i3+2)*x(j3+1)) )

            corfvyz = corfvyz + amass(i)*amass(j)*  &
               (   del(i3+1,j3+1)/(x(i3+3)*x(j3+2)) &
                 - del(i3+3,j3+2)/(x(i3+1)*x(j3+1)) &
                 - del(i3+1,j3+1)/(x(i3+3)*x(j3+2)) &
                 + del(i3+3,j3+2)/(x(i3+1)*x(j3+1)) )

         end do
      end do

   end if

   corfvx = corfvx/(totmass*totmass)
   corfvy = corfvy/(totmass*totmass)
   corfvz = corfvz/(totmass*totmass)
   corfvxy = corfvxy/(totmass*totmass)
   corfvxz = corfvxz/(totmass*totmass)
   corfvyz = corfvyz/(totmass*totmass)
   write(18,'(f8.3,2x,6f16.8)') t/20.455,corfvx,corfvy,corfvz,corfvxy, &
         corfvxz, corfvyz
#endif

   end do  ! nt
   return
end subroutine corfv 
