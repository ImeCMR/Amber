#include "copyright.h"
#include "../include/dprec.fh"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)

subroutine bicg(xm,ym,zm,nbnd,nz_num,epsin,epsout,h,&
   c2,index,index2,bv,phi,xs,maxitn,itn,tol,inorm,norm)
   
!  bicg csr wrapper for both cpu and gpu calls
!
!  tol <- accept
!  bv <- f
!  phi <- u
!  xs  <- u0
   
   implicit none
   integer xm,ym,zm,xmymzm,nbnd,nz_num
   _REAL_ :: c2(nbnd,27)
   _REAL_ :: bv(1:xm*ym*zm)
   _REAL_ :: phi(1:xm*ym*zm),xs(1:xm*ym*zm)
   integer :: index(xm,ym,zm),index2(xm,ym,zm)
   integer maxitn, itn, ier
   _REAL_ epsin, epsout, h
   _REAL_ tol
   _REAL_ inorm, norm

   _REAL_, parameter :: FOURPI=4.0d0*3.14159265358979323846

   _REAL_,allocatable :: a(:),u(:),f(:)
   integer,allocatable :: ia(:),ja(:)

#ifdef nvCUDA
   real stol, snorm
   real, allocatable :: su(:), sf(:), sa(:)
#endif

   xmymzm = xm*ym*zm

   if ( allocated(a) ) then
      deallocate(a, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( a(1:nz_num), stat=ier); REQUIRE(ier==0)
   if ( allocated(ja) ) then
      deallocate(ja, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( ja(1:xmymzm+1), stat=ier); REQUIRE(ier==0) ! this is row index
   if ( allocated(ia) ) then
      deallocate(ia, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( ia(1:nz_num), stat=ier); REQUIRE(ier==0) ! this is column index
   if ( allocated(f) ) then
      deallocate(f, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( f(1:xmymzm), stat=ier); REQUIRE(ier==0)

   call setcsr(xm,ym,zm,nbnd,nz_num,index,index2,c2,bv,a,f,ia,ja)

#ifdef nvCUDA
   if ( allocated(sa) ) then
      deallocate(sa, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( sa(1:nz_num), stat=ier); REQUIRE(ier==0)
   if ( allocated(su) ) then
      deallocate(su, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( su(1:xmymzm), stat=ier); REQUIRE(ier==0)
   if ( allocated(sf) ) then
      deallocate(sf, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( sf(1:xmymzm), stat=ier); REQUIRE(ier==0)

   ja = ja - 1 ! convert to c style before passing
   ia = ia - 1 ! convert to c style before passing
   sa = a
   sf = f
   su = xs
   stol = tol
   inorm = sqrt(dot_product(sf, sf))
   ! ---- CUSPARSE version
   call cusparse_bicg_wrapper(su, sf, ja, ia, sa, xmymzm, nz_num, maxitn, stol, itn, snorm)
   ! ---- CUSP version
   !call cusp_bicg_wrapper(su, sf, ja, ia, sa, xmymzm, nz_num, maxitn, stol, itn, snorm)
   norm = snorm
   xs = su
   !! ---- double precision below ----
   !! ---- CUSPARSE only
   !!inorm = sqrt(dot_product(f, f))
   !!call cusparse_bicg_wrapper( xs,  f, ja, ia,  a, xmymzm, nz_num, maxitn,  tol, itn,  norm)
   !!if ( norm >= tol*inorm .or. itn >= maxitn ) then
   !!   xs(1:xmymzm) = 0.0d0
   !!   phi(1:xmymzm) = 0.0d0
   !!   write(6,*) 'PBSA WARNING: BICG failed to converge!'
   !!else
   !!   phi(1:xmymzm) = xs(1:xmymzm)
   !!end if
#else
   inorm = sqrt(dot_product(f, f))
   call cpu_bicg_wrapper( xs,  f, ja, ia,  a, xmymzm, nz_num, maxitn,  tol, itn,  norm)
#endif
   if ( norm >= tol*inorm .or. itn >= maxitn ) then
      xs(1:xmymzm) = 0.0d0
      phi(1:xmymzm) = 0.0d0
      write(6,*) 'PBSA WARNING: BICG failed to converge!'
   else
      phi(1:xmymzm) = xs(1:xmymzm)
   end if

   !write(6,*) 'BICG/CSR Initial Norm', inorm
   !write(6,*) 'BICG/CSR Final Norm', norm
   !write(6,*) 'BICG/CSR Iteration needed', itn

contains

subroutine setcsr(l,m,n,maxirr,nz_num,index,index2,c2,ff,a,f,ia,ja)

   implicit none
   integer l,m,n,maxirr,nz_num
   integer :: index(l,m,n), index2(l,m,n)
   _REAL_ :: c2(maxirr,27)
   _REAL_ :: ff(l,m,n)
   _REAL_ :: a(nz_num),f(l*m*n)
   integer :: ia(nz_num), ja(l*m*n+1)

   integer ncount,i,j,k,ne,ii,ir,nc1,i0,j0,k0,ndis,ne1
   integer :: position1

   ! ncount tells which ipb option this is, 0 for ipb6 and 1 for ipb7 and ipb8
   ncount = 0
   if ( nz_num == 7*l*m*n-2*(1+l+l*m) ) ncount = 1

   nz_num = 0
   do k = 1, n
   do j = 1, m
   do i = 1, l
 
      position1=npos(l,m,n,i,j,k)
      ja(position1)=nz_num+1
   
      if ( index(i,j,k) == 1 ) then
   
         if (k-1 >= 1) then
            nz_num=nz_num+1
            a(nz_num) = -epsin/h/h
            ia(nz_num) = npos(l,m,n,i,j,k-1)
         end if

         if (j-1 >= 1) then
            nz_num=nz_num+1
            a(nz_num) = -epsin/h/h
            ia(nz_num) = npos(l,m,n,i,j-1,k)
         end if
         if ((j-1 == 0).and.(k > 1)) then
            nz_num=nz_num+1
            a(nz_num) = 0.d0
            ia(nz_num) = npos(l,m,n,i,j-1,k)
         end if

         if (i-1 >= 1) then
            nz_num=nz_num+1
            a(nz_num) = -epsin/h/h
            ia(nz_num) = npos(l,m,n,i-1,j,k)
         end if
         if ((i-1 == 0).and.((k > 1).or.(j > 1))) then
            nz_num=nz_num+1
            a(nz_num) = 0.d0
            ia(nz_num) = npos(l,m,n,i-1,j,k)
         end if
   
         nz_num=nz_num+1
         a(nz_num) = 6.d0*epsin/h/h
         ia(nz_num) = position1
         f(position1) = ff(i,j,k)/h/h/h*FOURPI
   
         if (i+1 <= l) then
            nz_num=nz_num+1
            a(nz_num) = -epsin/h/h
            ia(nz_num) = npos(l,m,n,i+1,j,k)
         end if
         if ((i+1 == l+1).and.((k < n).or.(j < m))) then
            nz_num=nz_num+1
            a(nz_num) = 0.d0
            ia(nz_num) = npos(l,m,n,i+1,j,k)
         end if

         if (j+1 <= m) then
            nz_num=nz_num+1
            a(nz_num) = -epsin/h/h
            ia(nz_num) = npos(l,m,n,i,j+1,k)
         end if
         if ((j+1 == m+1).and.(k < n)) then
            nz_num=nz_num+1
            a(nz_num) = 0.d0
            ia(nz_num) = npos(l,m,n,i,j+1,k)
         end if

         if (k+1 <= n) then
            nz_num=nz_num+1
            a(nz_num) = -epsin/h/h
            ia(nz_num) = npos(l,m,n,i,j,k+1)
         end if
   
      else if (index(i,j,k) == 5 ) then
   
         if (k-1 >= 1) then
            nz_num=nz_num+1
            a(nz_num) = -epsout/h/h
            ia(nz_num) = npos(l,m,n,i,j,k-1)
         end if

         if (j-1 >= 1) then
            nz_num=nz_num+1
            a(nz_num) = -epsout/h/h
            ia(nz_num) = npos(l,m,n,i,j-1,k)
         end if
         if ((j-1 == 0).and.(k > 1)) then
            nz_num=nz_num+1
            a(nz_num) = 0.d0
            ia(nz_num) = npos(l,m,n,i,j-1,k)
         end if

         if (i-1 >= 1) then
            nz_num=nz_num+1
            a(nz_num) = -epsout/h/h
            ia(nz_num) = npos(l,m,n,i-1,j,k)
         end if
         if ((i-1 == 0).and.((k > 1).or.(j > 1))) then
            nz_num=nz_num+1
            a(nz_num) = 0.d0
            ia(nz_num) = npos(l,m,n,i-1,j,k)
         end if
   
         nz_num=nz_num+1
         a(nz_num) = 6.d0*epsout/h/h
         ia(nz_num) = position1
         f(position1) = ff(i,j,k)/h/h/h*FOURPI
   
         if (i+1 <= l) then
            nz_num=nz_num+1
            a(nz_num) = -epsout/h/h
            ia(nz_num) = npos(l,m,n,i+1,j,k)
         end if
         if ((i+1 == l+1).and.((k < n).or.(j < m))) then
            nz_num=nz_num+1
            a(nz_num) = 0.d0
            ia(nz_num) = npos(l,m,n,i+1,j,k)
         end if

         if (j+1 <= m) then
            nz_num=nz_num+1
            a(nz_num) = -epsout/h/h
            ia(nz_num) = npos(l,m,n,i,j+1,k)
         end if
         if ((j+1 == m+1).and.(k < n)) then
            nz_num=nz_num+1
            a(nz_num) = 0.d0
            ia(nz_num) = npos(l,m,n,i,j+1,k)
         end if

         if (k+1 <= n) then
            nz_num=nz_num+1
            a(nz_num) = -epsout/h/h
            ia(nz_num) = npos(l,m,n,i,j,k+1)
         end if
 
      else
   
         f(position1) = -ff(i,j,k)/h/h/h*FOURPI

         if (ncount == 0) then 

            do k0 = k-1, k+1
            do j0 = j-1, j+1
            do i0 = i-1, i+1
               ii=npos(3,3,3,k0-k+2,j0-j+2,i0-i+2)
               ir = index2(i,j,k)
               if ( abs(c2(ir,ii)) > 1.d-10 ) then
                  nz_num = nz_num + 1
                  a(nz_num) = -c2(ir,ii)
                  ia(nz_num) = npos(l,m,n,i0,j0,k0)
               end if
            end do
            end do
            end do

         else

            do k0 = k-1, k+1
            do j0 = j-1, j+1
            do i0 = i-1, i+1
               if (abs(i-i0)+abs(j-j0)+abs(k-k0) <= 1) then
                  if (i0-i == -1) ii = 1
                  if (j0-j == -1) ii = 2
                  if (k0-k == -1) ii = 3
                  if (i0 == i .and. j0==j .and. k0==k) ii = 4
                  if (k0-k == 1) ii = 5
                  if (j0-j == 1) ii = 6
                  if (i0-i == 1) ii = 7

                  ir = index2(i,j,k)
                  if ( abs(c2(ir,ii)) > 1.d-10 ) then
                     nz_num = nz_num + 1
                     a(nz_num) = -c2(ir,ii)
                     ia(nz_num) = npos(l,m,n,i0,j0,k0)
                  end if
               end if
            end do
            end do
            end do

         end if
   
      end if
  
   end do
   end do
   end do
   ja(l*m*n+1)=nz_num+1

end subroutine setcsr

integer function npos(l,m,n,i,j,k)

   implicit none
   integer l,m,n,i,j,k
   npos = i + (j-1)*l + (k-1)*l*m

end function npos

subroutine cpu_bicg_wrapper(x, b, row, col, val, xmymzm, nz_num, maxitn, accept, itn, residual)

   implicit none
   integer xmymzm, nz_num, maxitn, itn
   _REAL_ x(xmymzm), b(xmymzm), val(nz_num)
   integer row(xmymzm+1), col(nz_num)
   _REAL_ accept, residual

   _REAL_, allocatable :: p(:), q(:)
   _REAL_, allocatable :: r(:), r0(:), t(:)

   _REAL_ tol2
   _REAL_ nrmr, nrmr0, rho, rhop
   _REAL_ alpha, beta, omega
   integer k

   tol2 = accept**2 ! accept is 1-norm in pbsa

   if ( allocated(p) ) then
      deallocate(p, stat = ier); REQUIRE(ier==0)
   end if
   allocate(p(xmymzm), stat=ier); REQUIRE(ier==0)
   if ( allocated(q) ) then
      deallocate(q, stat = ier); REQUIRE(ier==0)
   end if
   allocate(q(xmymzm), stat=ier); REQUIRE(ier==0)
   if ( allocated(r) ) then
      deallocate(r, stat = ier); REQUIRE(ier==0)
   end if
   allocate(r(xmymzm), stat=ier); REQUIRE(ier==0)
   if ( allocated(r0) ) then
      deallocate(r0, stat = ier); REQUIRE(ier==0)
   end if
   allocate(r0(xmymzm), stat=ier); REQUIRE(ier==0)
   if ( allocated(t) ) then
      deallocate(t, stat = ier); REQUIRE(ier==0)
   end if
   allocate(t(xmymzm), stat=ier); REQUIRE(ier==0)

   ! initialize r, r0, p, norm
   r = b
   r0 = r
   p = r
   nrmr0 = dot_product(r,r) 
   nrmr = nrmr0
   rho = nrmr0

   ! BICG iteration
   do k = 0, maxitn
      ! compute q=Ap & alpha
      call csrmv(xmymzm,nz_num,val,row,col,p,q)
      alpha = rho/dot_product(r0,q)
      ! compute s = r - \alpha q, note s is just r
      r = r - alpha * q
      ! compute omega = (t^{T} s) / (t^{T} t), t = As
      call csrmv(xmymzm,nz_num,val,row,col,r,t)
      omega = dot_product(t,r)/dot_product(t,t)
      ! update x = x + alpha p + omega s
      ! update r = s - omega t, t = As
      x = x + alpha * p + omega * r
      r = r - omega * t
      ! compute new norm
      nrmr = dot_product(r,r)
      !print *, " itn ", k, " residual ", sqrt(nrmr), " tol2 ", tol2, " init norm ", nrmr0
      if (nrmr < tol2*nrmr0) exit
      ! compute beta = (rho_i/rho_i-1) (alpha/omega)
      rhop = rho;
      rho = dot_product(r0,r)
      beta = (rho/rhop)*(alpha/omega)
      ! compute p = r + beta (p - omega q), q = Ap
      p = r + beta * ( p - omega * q)
   end do

   itn = k
   residual = sqrt(nrmr) ! return sqrt'd 2-norm

   deallocate(p,q,r,r0,t)

end subroutine cpu_bicg_wrapper

subroutine csrmv(n,nz,val,row,col,x,y)

   implicit none
   integer :: n, nz
   integer :: row(n+1), col(nz)
   _REAL_ :: val(nz), x(n), y(n)

   integer i, j
   do i = 1, n
      y(i)  = 0.0d0
      do j = row(i), row(i+1) - 1
          y(i) = y(i) + val(j) * x(col(j))
      end do
   end do

end subroutine csrmv

end subroutine bicg

#ifdef nvCUDA
subroutine bicg_dia(xm,ym,zm,nbnd,epsin,epsout,h,&
   c2,index,index2,bv,phi,xs,maxitn,itn,tol,inorm,norm)

   ! bicg dia wrapper for gpu calls only
   !
   ! tol <- accept
   ! bv <- f
   ! phi <- u
   ! xs  <- u0

   implicit none
   integer xm, ym, zm, nbnd
   _REAL_ :: c2(nbnd,7)
   _REAL_ :: bv(1:xm*ym*zm)
   _REAL_ :: phi(1:xm*ym*zm),xs(1:xm*ym*zm)
   integer :: index(xm,ym,zm),index2(xm,ym,zm)
   integer maxitn, itn
   _REAL_ epsin, epsout, h
   _REAL_ tol
   _REAL_ inorm, norm

   _REAL_, parameter :: FOURPI=4.0d0*3.14159265358979323846

   integer :: ioffset(7)
   real stol, snorm
   real, allocatable :: su(:), sf(:), sa(:)

   integer xmymzm, ier

   xmymzm = xm*ym*zm

   if ( allocated(sa) ) then
      deallocate(sa, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( sa(1:7*xmymzm), stat=ier); REQUIRE(ier==0)
   if ( allocated(su) ) then
      deallocate(su, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( su(1:xmymzm), stat=ier); REQUIRE(ier==0)
   if ( allocated(sf) ) then
      deallocate(sf, stat = ier); REQUIRE(ier==0)
   end if
   allocate ( sf(1:xmymzm), stat=ier); REQUIRE(ier==0)

   call setdia(xm,ym,zm,nbnd,epsin,epsout,h,c2,index,index2,bv,ioffset,sa,sf)

   stol = tol
   su(1:xmymzm) = xs(1:xmymzm)
   inorm = sqrt(dot_product(sf, sf))
   ! ---- CUSPARSE version
   !call cusparse_bicg_wrapper(x, b, ioffset, sa, xm, ym, zm maxitn, saccept, l_itn, snorm)
   ! ---- CUSP version
   call cusp_bicg_wrapper(su, sf, ioffset, sa, xm, ym, zm, maxitn, stol, itn, snorm)
   norm = snorm

   if ( norm >= tol*inorm .or. itn >= maxitn ) then
      xs(1:xmymzm) = 0.0d0
      phi(1:xmymzm) = 0.0d0
      write(6,*) 'PBSA WARNING: BICG failed to converge!'
   else
      xs(1:xmymzm) = su(1:xmymzm)
      phi(1:xmymzm) = su(1:xmymzm)
   end if

   !write(6,*) 'BICG/DIA Initial Norm', inorm
   !write(6,*) 'BICG/DIA Final Norm', norm
   !write(6,*) 'BICG/DIA Iteration needed', itn

contains

subroutine setdia(xm,ym,zm,nbnd,epsin,epsout,h,c2,index,index2,bv,offset,dia,f)
!
!  all the arrays here have been properly padded with zeros, to form an n-element
!  array, the same as diagnol one. When putting all together, ordered as 1, 2, 3,
!  dia, 5, 6, 7, to form a 7n long array, the non-zero elements should be:
!  z-1 term: off_dia1: 1+xm*ym ~ n
!  y-1 term: off_dia2: 1+xm+n ~ 2n
!  x-1 term: off_dia3: 1+1+2n ~ 3n
!  dia term:      dia: 1+3n ~ 4n
!  x+1 term: off_dia5: 1+4n ~ 5n-1
!  y+1 term: off_dia6: 1+5n ~ 6n-xm
!  z+1 term: off_dia7: 1+6n ~ 7n-xm*ym
!

   implicit none
   integer xm,ym,zm,nbnd
   _REAL_ epsin, epsout, h
   _REAL_ :: c2(nbnd,7)
   integer :: index(xm,ym,zm),index2(xm,ym,zm)
   _REAL_ :: bv(xm,ym,zm)
   integer :: offset(7)
   real :: dia(1:7*xm*ym*zm)
   real :: f(1:xm*ym*zm)

   integer xmymzm
   integer :: i, j, k, ncount
   _REAL_, parameter :: FOURPI=4.0d0*3.14159265358979323846

   xmymzm = xm*ym*zm

   ncount = 0
   do k = 1, zm
   do j = 1, ym
   do i = 1, xm

      ncount = ncount + 1

      if ( index(i,j,k) == 1 ) then
         if ( k-1 >= 1   ) dia(ncount         ) = -epsin/h/h
         if ( k-1 == 0   ) dia(ncount         ) = 0d0
         if ( j-1 >= 1   ) dia(ncount+  xmymzm) = -epsin/h/h
         if ( j-1 == 0   ) dia(ncount+  xmymzm) = 0d0
         if ( i-1 >= 1   ) dia(ncount+2*xmymzm) = -epsin/h/h
         if ( i-1 == 0   ) dia(ncount+2*xmymzm) = 0d0
                           dia(ncount+3*xmymzm) = 6.d0*epsin/h/h
         if ( i+1 <= xm  ) dia(ncount+4*xmymzm) = -epsin/h/h
         if ( i+1 == xm+1) dia(ncount+4*xmymzm) = 0d0
         if ( j+1 <= ym  ) dia(ncount+5*xmymzm) = -epsin/h/h
         if ( j+1 == ym+1) dia(ncount+5*xmymzm) = 0d0
         if ( k+1 <= zm  ) dia(ncount+6*xmymzm) = -epsin/h/h
         if ( k+1 == zm+1) dia(ncount+6*xmymzm) = 0d0
                             f(ncount) = bv(i,j,k)/h/h/h*FOURPI
      else if ( index(i,j,k) == 5 ) then
         if ( k-1 >= 1   ) dia(ncount         ) = -epsout/h/h
         if ( k-1 == 0   ) dia(ncount         ) = 0d0
         if ( j-1 >= 1   ) dia(ncount+  xmymzm) = -epsout/h/h
         if ( j-1 == 0   ) dia(ncount+  xmymzm) = 0d0
         if ( i-1 >= 1   ) dia(ncount+2*xmymzm) = -epsout/h/h
         if ( i-1 == 0   ) dia(ncount+2*xmymzm) = 0d0
                           dia(ncount+3*xmymzm) = 6.d0*epsout/h/h
         if ( i+1 <= xm  ) dia(ncount+4*xmymzm) = -epsout/h/h
         if ( i+1 == xm+1) dia(ncount+4*xmymzm) = 0d0
         if ( j+1 <= ym  ) dia(ncount+5*xmymzm) = -epsout/h/h
         if ( j+1 == ym+1) dia(ncount+5*xmymzm) = 0d0
         if ( k+1 <= zm  ) dia(ncount+6*xmymzm) = -epsout/h/h
         if ( k+1 == zm+1) dia(ncount+6*xmymzm) = 0d0
                             f(ncount) = bv(i,j,k)/h/h/h*FOURPI
      else
         dia(ncount         ) = -c2(index2(i,j,k),3)
         dia(ncount+  xmymzm) = -c2(index2(i,j,k),2)
         dia(ncount+2*xmymzm) = -c2(index2(i,j,k),1)
         dia(ncount+3*xmymzm) = -c2(index2(i,j,k),4)
         dia(ncount+4*xmymzm) = -c2(index2(i,j,k),7)
         dia(ncount+5*xmymzm) = -c2(index2(i,j,k),6)
         dia(ncount+6*xmymzm) = -c2(index2(i,j,k),5)
           f(ncount) = -bv(i,j,k)/h/h/h*FOURPI
      end if

   end do
   end do
   end do

   offset = (/-xm*ym, -xm, -1, 0, 1, xm, xm*ym/)

end subroutine setdia

end subroutine bicg_dia
#endif
