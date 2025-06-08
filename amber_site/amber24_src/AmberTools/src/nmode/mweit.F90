!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mweit here]
subroutine mweit (n, neig, ind, z, amass, winv)
   implicit double precision (a-h,o-z)
   
   dimension ind(n), z(n,n), amass(n/6), winv(n/2)
   
   k = 1
   do i=1,n/6
      winv(k) = 1.0/sqrt(amass(i))
      winv(k+1) = winv(k)
      winv(k+2) = winv(k)
      k = k + 3
   end do
   
   do kj = 1, neig
      j = ind (kj)
      do i = 1, n/2
         z(i,j) = z(i,j) * winv(i)
      end do
   end do
   
   return
end subroutine mweit 
