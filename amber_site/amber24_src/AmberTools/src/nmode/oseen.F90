!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine oseen here]
subroutine oseen (ioseen, natom, nr3, nr6, a, eta, expos, x, &
      gamma, amass)
   
   implicit double precision (a-h,o-z)
   dimension a(nr6,nr6), expos(natom), x(nr3), amass(natom), det(2)
   dimension gamma(*)
   
   sxpita = 6.0 * 3.14159 * eta
   etpita = 8.0 * 3.14159 * eta
   ont    = 1.0 / 3.0
   
   i = 0
   do iat = 1, natom
      if (expos(iat) /= 0.0) then
         
         !         ----- diagonal blocks
         
         tat = 1.0 / (sxpita * expos(iat))
         do ix = 1, 3
            do jx = 1, 3
               if (ix == jx) then
                  a(i+ix,i+jx) = tat
               else
                  a(i+ix,i+jx) = 0.0
               end if
            end do
         end do
         
         !         ----- off-diagonal blocks
         
         j = i + 3
         iat3 = iat * 3
         xi = x(iat3-2)
         yi = x(iat3-1)
         zi = x(iat3)
         ix = i + 1
         iy = i + 2
         iz = i + 3
         do jat = iat+1, natom
            if (expos(jat) /= 0.0) then
               jat3 = jat * 3
               delx = xi - x(jat3-2)
               dely = yi - x(jat3-1)
               delz = zi - x(jat3)
               r2 = delx*delx + dely*dely + delz*delz
               r = sqrt(r2)
               jx = j + 1
               jy = j + 2
               jz = j + 3
               
               !             ----- oseen tensor
               
               const = 1.0 / (etpita*r)
               a(ix,jx) = const * ( 1.0 + delx*delx / r2 )
               a(iy,jy) = const * ( 1.0 + dely*dely / r2 )
               a(iz,jz) = const * ( 1.0 + delz*delz / r2 )
               a(ix,jy) = const * delx * dely / r2
               a(ix,jz) = const * delx * delz / r2
               a(iy,jz) = const * dely * delz / r2
               a(iy,jx) = a(ix,jy)
               a(iz,jx) = a(ix,jz)
               a(iz,jy) = a(iy,jz)
               
               !             ----- rotne prager correction
               
               if (ioseen > 1) then
                  const = const / r2 * (expos(iat)**2 + expos(jat)**2)
                  a(ix,jx) = a(ix,jx) + const * ( ont - delx*delx / r2 )
                  a(iy,jy) = a(iy,jy) + const * ( ont - dely*dely / r2 )
                  a(iz,jz) = a(iz,jz) + const * ( ont - delz*delz / r2 )
                  a(ix,jy) = a(ix,jy) - const * delx * dely / r2
                  a(ix,jz) = a(ix,jz) - const * delx * delz / r2
                  a(iy,jz) = a(iy,jz) - const * dely * delz / r2
                  a(iy,jx) = a(ix,jy)
                  a(iz,jx) = a(ix,jz)
                  a(iz,jy) = a(iy,jz)
               end if
               j = j + 3
            end if  ! (expos(jat) /= 0.0)
         end do
         i = i + 3
      end if  
   end do
   
   ntdim = i
   
   do i = 1, ntdim
      do j = 1, i-1
         a(i,j) = a(j,i)
      end do
   end do
   
#ifdef LINPACK
   job = 1
   call dgeco (a, nr6, ntdim, a(1,nr3+2), rcond, a(1,nr3+1))
   call dgedi (a, nr6, ntdim, a(1,nr3+2), det, a(1,nr3+1), job)
#else
   
   ! --- place matrix in contiguous locations beginning in the middle
   !       of the complete matrix a(nr6,nr6):
   
   k = 0
   ic = nr3 + 1
   do i=1,ntdim
      do j=1,ntdim
         k = k + 1
         if (k > nr6) then
            ic = ic + 1
            k = 1
         end if
         a(k,ic) = a(i,j)
      end do
   end do
   
   !  --- invert the matrix
   
   call matinv(a(1,nr3+1),ntdim,det,a(1,1),a(1,2))
   
   !  --- restore the inverted matrix to its rightful place:
   
   k = 0
   ic = nr3 + 1
   do i=1,ntdim
      do j=1,ntdim
         k = k + 1
         if (k > nr6) then
            ic = ic + 1
            k = 1
         end if
         a(i,j) = a(k,ic)
      end do
   end do
#endif
   
   ij = 0
   ja = 0
   do j = 1, 3*natom
      jat = (j+2)/3
      if (expos(jat) == 0.0) then
         do i = 1, j
            ij = ij + 1
            gamma(ij) = 0.0
         end do
      else
         ja = ja + 1
         ia = 0
         do i = 1, j
            iat = (i+2)/3
            ij = ij + 1
            if (expos(iat) == 0.0) then
               gamma(ij) = 0.0
            else
               ia = ia + 1
               gamma(ij) = a(ia,ja) / sqrt (amass(iat)*amass(jat))
            end if
         end do
      end if
   end do
   
   return
end subroutine oseen 
