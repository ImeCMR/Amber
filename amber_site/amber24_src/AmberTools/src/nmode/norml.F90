!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine norml here]
subroutine norml (neig, nm, n, z, wr, wi, gamma, iclass, ind)
   
   !     ----- to compute normalized Langevin modes, i.e., the L matrix
   
   implicit double precision (a-h,o-z)
   complex piki, pili, const, sum, elki, lambda
   complex  minusi, zero
   dimension z(nm,n), wr(n), wi(n), gamma(*), iclass(neig), ind(n)
   
   index (k,l) = max(k,l) * (max(k,l)-1) / 2 + min(k,l)
   
   n2 = n/2
   minusi = (0.0, -1.0)
   zero   = (0.0,  0.0)
   rzero  = 0.0
   
   !     ----- this can be easily figured out.  I will even give you a hint.
   !     ----- see eqs. (3.10) and (4.1) Lamm and Szabo, JCP, 85, 7336 (1986)
   !     ----- also see routine rg for the arrays wr, wi and z.
   
   ik = 1
   do while (ik <= neig)
      i = ind(ik)
      lambda = cmplx (wr(i), wi(i))
      
      if (wi(i) /= 0.0) then
         
         !         ----- a complex pair
         
         sum = 0.0
         do k = 1, n2
            piki = cmplx (z(k,i), z(k,i+1))
            do l = 1, n2
               pili = cmplx (z(l,i), z(l,i+1))
               kl = index(k,l)
               sum = sum + piki * gamma(kl) * pili
            end do
            sum = sum + 2.0 * (lambda * piki * piki)
         end do
         if (sum == zero) then
            const = cmplx (1.0, 0.0)
         else
            const = sqrt (lambda / sum)
         end if
         do k = 1, n2
            piki = cmplx (z(k,i), z(k,i+1))
            elki = piki * const
            z(k,i) = real (elki)
            z(k,i+1) = real (minusi*elki)
         end do
         iclass(ik) = 1
         iclass(ik+1) = 1
         ik = ik + 2
         
      else if (wr(i) > 0.0) then
         
         !         ----- lambda is  positive real
         
         sum = 0.0
         do k = 1, n2
            do l = 1, n2
               kl = index(k,l)
               sum = sum + z(k,i) * gamma(kl) * z(l,i)
            end do
            sum = sum + 2.0 * wr(i) * z(k,i) * z(k,i)
         end do
         const = sqrt (lambda/sum)
         do k = 1, n2
            z(k,i) = z(k,i) * const
         end do
         iclass(ik) = 2
         ik = ik + 1
         
      else if (wr(i) < 0.0) then
         
         !         ----- negative real
         
         sum = 0.0
         do k = 1, n2
            do l = 1, n2
               kl = index(k,l)
               sum = sum + z(k,i) * gamma(kl) * z(l,i)
            end do
            sum = sum + 2.0 * wr(i) * z(k,i)*z(k,i)
         end do
         if (real(sum) > 0.0) then
            const = sqrt (  - lambda / sum )
            iclass(ik) = 3
         else if (real(sum) < 0.0) then
            const = sqrt (lambda / sum)
            iclass(ik) = 4
         end if
         do k = 1, n2
            z(k,i) = const * z(k,i)
         end do
         ik = ik + 1
         
      else
         
         !         ----- must be exactly zero
         
         iclass(ik) = 5
         ik = ik + 1
         
      end if
      
   end do
   
   return
end subroutine norml 
