!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine spsort here]
subroutine spsort (n, dbl, wr, wi, ind, nreal)
   
   !     ----- a special sort routine called by lmode
   
   implicit double precision (a-h,o-z)
   dimension dbl(n), ind(n), wr(n), wi(n)
   
   !     ----- First make ind in natural order
   
   do i = 1, n
      ind(i) = i
      dbl(i) = abs(wi(i))
   end do
   
   !     ----- sort keeping track of which element goes where
   
   do i = 1, n-1
      
      !       ----- get pointer to minimum value of remaining entries
      
      ip = i
      do j = i+1, n
         if (dbl(j) < dbl(ip)) ip = j
      end do
      
      !       ----- interchange ip-th and i-th dbl
      
      dt = dbl(i)
      dbl(i) = dbl(ip)
      dbl(ip) = dt
      
      !       ----- interchange ip-th and i-th ind
      
      it = ind(i)
      ind(i) = ind(ip)
      ind(ip) = it
      
   end do
   
   !     ---- oops! Bring back complex conjugates in original order
   !     ---- i.e., (a + ib) followed by (a - ib)
   
   i = 1
   do while (i <= n)
      it = ind(i)
      if (wi(it) /= 0.0 .and. wr(it) == wr(ind(i+1))) then
         if (wi(it) < 0.0) then
            ind(i) = ind(i+1)
            ind(i+1) = it
         end if
         i = i + 2
      else
         i = i + 1
      end if
   end do
   
   !     ----- count reals
   
   nreal = 0
   do while (abs(wi(ind(nreal+1))) < 1.0e-3 .and. nreal < n)
      nreal = nreal + 1
      dbl(nreal) = abs(wr(ind(nreal)))
   end do
   
   !     ----- sort reals
   
   do i = 1, nreal-1
      ip = i
      do j = i+1, nreal
         if (dbl(j) < dbl(ip)) ip = j
      end do
      dt = dbl(i)
      dbl(i) = dbl(ip)
      dbl(ip) = dt
      it = ind(i)
      ind(i) = ind(ip)
      ind(ip) = it
   end do
   
   return
end subroutine spsort 
