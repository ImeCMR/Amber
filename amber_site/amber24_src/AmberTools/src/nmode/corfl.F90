
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine corfl here]
subroutine corfl (nat3, nat6, nvect, kup, lup, ntrun, tf, np, bose, x)
   
   !     ----- to calculate correlation functions from langevin modes
   
   implicit double precision (a-h,o-z)
   double precision kt
   logical first,bose
   complex lambda, lkin, llin, lkjn, lljn, delijc, zero
   
   !     automatic arrays:
   dimension wr(nat6),wi(nat6),el(nat3,nvect),iclass(nat6)

   dimension e(3), del(3,3), del0(3,3), x(*)
#  include "files.h"
   !                  ----CONSQ = hc/2kT in cm, with T=300K
   !                       (use for quantum, Bose statistics)
   parameter (kt = 0.6)
   parameter (consq = 2.39805d-3)
   data first /.true./
   
   !     functions
   k(i) = 3*(kup-1) + i
   l(i) = 3*(lup-1) + i
   !     end functions
   
   zero = (0.0,0.0)
   rzero = 0.0
   
   !     ----- read langevin modes from file lmode
   
   if (first) call amopen (18, plot, owrite, 'F', 'W')
   first = .false.
   call lmin (nat3, nat6, nvect, wr, wi, el, iclass)
   
   !     ----- set tf and tdel
   
   ti = 0.0
   tf = tf * 20.455
   tdel = (tf-ti)/np
   
   !     ----- prepare plot file
   
   write(18,'(a51,i4,a1,i4)') &
         '#   time           C(t), normalized, for atom pair ',kup,'-',lup
   
   !    ---e(i) is the unit vector along the l-k bond:
   
   req = 0.0
   do i = 1, 3
      e(i) = x(k(i)) - x(l(i))
      req = req + e(i)*e(i)
   end do
   req = sqrt(req)
   do i = 1, 3
      e(i) = e(i) / req
   end do
   
   t = ti
   do nt = 1,np+1
      
      !     ----- for each t, calculate the correlation matrix for delta
      !        as in Eq. 7.16 of lamm and szabo j. chem. phys. 85, 7334 (1986)
      !        Note that the rhs of this eq. should be multiplied by kT
      
      do i = 1, 3
         ki = k(i)
         li = l(i)
         do j = 1, 3
            kj = k(j)
            lj = l(j)
            del(i,j) = 0.0
            n = 1   ! orig: n = 13
            
            do while (n.le.nvect)
               lambda = cmplx ( wr(n), wi(n) )
               if (bose) then
                  argq = sqrt(wr(n)**2 + wi(n)**2)*108.59*consq
                  qcorr = argq/tanh(argq)
#  ifdef debug
                  if (wr(n) /= 0.0 .and. t == 0.0) &
                        write(6,'(i4,5f10.4)') n,wi(n)*108.59,1./(wr(n)*20.455), &
                        sqrt(wr(n)**2 + wi(n)**2)*108.59,argq,qcorr
#  endif
               end if
               if (iclass(n) == 5 .or. lambda == zero) then
                  !                                                     ---skip this mode
                  n = n + 1
               else if (iclass(n) == 2 .or. iclass(n) == 4) then
                  !                       ---purely real mode
                  lkin = el(ki,n)
                  llin = el(li,n)
                  lkjn = el(kj,n)
                  lljn = el(lj,n)
                  !               write(6,'(i4,3e14.5)') iclass(n),lambda,t
                  !               if (real(lambda)*t .lt. -50.) go to 35
                  delijc = (exp(lambda*t)/(lambda*lambda)) &
                        * (lkin - llin) * (lkjn - lljn)
                  if (bose) delijc = delijc*qcorr
                  del(i,j) = del(i,j) - kt*real(delijc)
                  n = n + 1
               else if (iclass(n) == 3) then
                  !                        ---purely imaginary mode
                  lkin = cmplx (rzero, el(ki,n))
                  llin = cmplx (rzero, el(li,n))
                  lkjn = cmplx (rzero, el(kj,n))
                  lljn = cmplx (rzero, el(lj,n))
                  !               write(6,'(i4,3e14.5)') iclass(n),lambda,t
                  !               if (real(lambda)*t .lt. -50.) go to 35
                  delijc = (exp(lambda*t)/(lambda*lambda)) &
                        * (lkin - llin) * (lkjn - lljn)
                  if (bose) delijc = delijc*qcorr
                  del(i,j) = del(i,j) - kt*real(delijc)
                  n = n + 1
               else if (iclass(n) == 1) then
                  !                        ---complex mode
                  lkin = cmplx ( el(ki,n), el(ki,n+1) )
                  llin = cmplx ( el(li,n), el(li,n+1) )
                  lkjn = cmplx ( el(kj,n), el(kj,n+1) )
                  lljn = cmplx ( el(lj,n), el(lj,n+1) )
                  !               write(6,'(i4,3e14.5)') iclass(n),lambda,t
                  !               if (real(lambda)*t .lt. -50.) go to 35
                  delijc = (exp(lambda*t)/(lambda*lambda)) &
                        * (lkin - llin) * (lkjn - lljn)
                  if (bose) delijc = delijc*qcorr
                  del(i,j) = del(i,j) - 2.0 * kt*real(delijc)
                  n = n + 2
               end if  ! (iclass(n) == 5 .or. lambda == zero)
            end do
            if (t == 0.0) del0(i,j) = del(i,j)
         end do
      end do
      
      
      if (ntrun == 1) then
         
         !         ----- correlation in length, Eq. (10.2) of lamm and szabo:
         
         rtr0 = 0.0
         do i = 1, 3
            do j = 1,3
               rtr0 = rtr0 + e(i)*e(j)*del(i,j)
            end do
         end do
         if (t == ti) r0r0 = rtr0
         write(18,'(f10.3,2x,f12.5)') t/20.455,rtr0/r0r0
         
      else if (ntrun == 2) then
         
         !         ----- librational correlation function, using eq. 7.12
         !               of lamm and szabo (without beta on the lhs):
         
         sum = 0.0
         do i = 1, 3
            sum = sum - del0(i,i) + del(i,i)
            do j = 1,3
               sum = sum + (del0(i,j) - del(i,j)) * e(i)*e(j)
            end do
         end do
         p2 = 1.0 + (3.0/req**2) * sum
         write(18,'(f10.3,2x,f12.5)') t/20.455,p2
         
      else if (ntrun == 3) then
         
         !    ---- use eq. 45 of Henry & Szabo [JCP 82:4753(1985)]
         !         to get motional correction for dipole-dipole correlation:
         
         suma = 0.0
         do i = 1, 3
            suma = suma - 3.0*(del0(i,i) + del(i,i))
            do j = 1,3
               suma = suma + 15.0*((del0(i,j) - del(i,j)) * e(i)*e(j))
            end do
         end do
         gamma =  1.0 + (0.5/req**2) * suma
         gamma = gamma*gamma
         write(18,'(f10.3,2x,f12.5)') t/20.455,gamma
         
      end if
      
      t = t + tdel
   end do
   
   return
end subroutine corfl 
