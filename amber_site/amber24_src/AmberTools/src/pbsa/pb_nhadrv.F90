! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"
#include "timer.h"

module nha

#include "pb_constants.h"

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Algoritm interface routine
subroutine pb_nhadrv( npbstep,npbgrid,nstlim,atmfirst,atmlast,npbopt,solvopt,level,nfocus,bcopt,&
                      natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                      xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                      maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                      pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                      gcrd,acrg,&
                      nbnd,iepsav,insas,lvlset,&
                      chgrd,saltgrd,phi,&
                      bv,cphi,xs,&
                      ngrdcrg,grdcrg,qgrdcrg,ipb&
                    )

   use pbtimer_module
   implicit none

   ! all the driver variables are shared among all "contained" routines, so are
   ! not redeclared in the containted routines anymore, except there is a need
   ! to remap their dimensions and to copy to other variables.

   ! passed variables

   integer npbstep, npbgrid, nstlim, atmfirst, atmlast,natom
   integer npbopt, solvopt, level, nfocus, bcopt
   _REAL_  h, savh(nfocus), gox, goy, goz, savgox(nfocus), savgoy(nfocus), savgoz(nfocus)
   integer xm, ym, zm, xmym, xmymzm, savxm(nfocus), savym(nfocus), savzm(nfocus)
   integer maxitn, itn
   _REAL_  fmiccg, accept, laccept, wsor, lwsor, inorm, norm
   _REAL_  pbkappa, pbkb, pbtemp, ivalence, istrng, eps0, epsin, epsout, ionene
   _REAL_  gcrd(3,natom), acrg(natom)
   integer nbnd, iepsav(4,xmymzm)

   integer ngrdcrg,grdcrg(3,*),ipb
   _REAL_ qgrdcrg(*)

   integer insas(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)
   _REAL_  lvlset(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)

   _REAL_  chgrd(xm,ym,zm), saltgrd(xmymzm), phi(xmymzm)
   _REAL_  bv(xm,ym,zm), cphi(xmymzm), xs(xmymzm+2*xmym)
   integer ii,i,j,k

   ! local varialbes

   _REAL_ rh

   rh = ONE/h

   call pbtimer_start(PBTIME_PBBUILDSYS)

   if ( bcopt /= 2 .and. bcopt /= 6 ) then
      write(6,*) 'PB Bomb in pb_nhndrv(): bcopt can only be 2 or 6'
      call mexit(6,1)
   end if

   ! a. compute FD coulumbic potentials on interface grid points
   ! this is only needed for singular free PB
   ! bcopt == 2 is for conductor boundary
   ! bcopt == 6 is for free boundary
   ! bcopt == 11 is reserved for periodic boundary, not done yet

   cphi = ZERO
   call pb_dbcgrd( cphi(1), insas )

   ! b. convert grid boundary potentials into effective charges and store in bv()

   bv = ZERO
   call pb_bndcnd( bv(1,1,1), chgrd(1,1,1) )

   ! c. the interface grid points' jump conditions will also be stored in
   ! bv to be done inside iim()

   call pbtimer_stop(PBTIME_PBBUILDSYS)

   call iim(gox,goy,goz,xm,ym,zm,lvlset,insas,nbnd,iepsav, &
            epsin/eps0,epsout/eps0, &
            bv(1,1,1),phi(1),xs(1),h,atmfirst,atmlast,bcopt,&
            maxitn,itn,accept,norm,inorm,gcrd,natom,&
            ngrdcrg,grdcrg,qgrdcrg,chgrd,ipb,cphi&
           )

contains

subroutine pb_dbcgrd( cphi, insas )

   _REAL_ cphi(xm,ym,zm)
   integer insas(0:xm+1,0:ym+1,0:zm+1)

   integer i, j, k
   integer i0, j0, k0
   integer ip
   _REAL_ tmp

   tmp = ZERO

   do ip = 1, nbnd
      i0 = iepsav(1,ip); j0 = iepsav(2,ip); k0 = iepsav(3,ip)

      i = i0 - 1
      if ( insas(i ,j0,k0) > 0 ) then
         if ( cphi(i ,j0,k0) == ZERO ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
      end if

      i = i0 + 1
      if ( insas(i ,j0,k0) > 0 ) then
         if ( cphi(i ,j0,k0) == ZERO ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
      end if

      j = j0 - 1
      if ( insas(i0,j ,k0) > 0 ) then
         if ( cphi(i0,j ,k0) == ZERO ) then
            call get_coulpot(i0,j ,k0,tmp);
            cphi(i0,j ,k0) = tmp/epsin
         end if
      end if

      j = j0 + 1
      if ( insas(i0,j ,k0) > 0 ) then
         if ( cphi(i0,j ,k0) == ZERO ) then
            call get_coulpot(i0,j ,k0,tmp)
            cphi(i0,j ,k0) = tmp/epsin
         end if
      end if

      k = k0 - 1
      if ( insas(i0,j0,k ) > 0 ) then
         if ( cphi(i0,j0,k ) == ZERO ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
      end if

      k = k0 + 1
      if ( insas(i0,j0,k ) > 0 ) then
         if ( cphi(i0,j0,k ) == ZERO ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
      end if

      if ( insas(i0,j0,k0) > 0 ) then
         if ( cphi(i0,j0,k0) == ZERO ) then
            call get_coulpot(i0,j0,k0,tmp)
            cphi(i0,j0,k0) = tmp/epsin
         end if
      end if

   end do

end subroutine pb_dbcgrd

subroutine get_coulpot(i,j,k,pot)

   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green

   integer i,j,k
   _REAL_ pot

   integer iatm
   integer itmp,jtmp,ktmp
   integer idx,idy,idz

   _REAL_ factor,qtmp,rinv,xtmp,ytmp,ztmp,dx,dy,dz
   _REAL_ a,a1,b,b1,c,c1

   factor = ONE/(FOURPI)/h

   pot = ZERO
   do iatm = atmfirst, atmlast
      xtmp = gcrd(1,iatm); ytmp = gcrd(2,iatm); ztmp = gcrd(3,iatm)
      qtmp = factor*acrg(iatm)

      dx = abs(i-xtmp); dy = abs(j-ytmp); dz = abs(k-ztmp)
      if (dx < 40.d0 .and. dy < 40.d0 .and. dz < 40.d0) then
         idx = floor(dx); idy = floor(dy); idz = floor(dz)
         a=dx-idx;b=dy-idy;c=dz-idz
         a1 = 1 - a; b1 = 1 - b; c1 = 1 - c
         rinv = a1*b1*c1*green(idx  ,idy  ,idz  ) &
               +a *b1*c1*green(idx+1,idy  ,idz  ) &
               +a1*b *c1*green(idx  ,idy+1,idz  ) &
               +a *b *c1*green(idx+1,idy+1,idz  ) &
               +a1*b1*c *green(idx  ,idy  ,idz+1) &
               +a *b1*c *green(idx+1,idy  ,idz+1) &
               +a1*b *c *green(idx  ,idy+1,idz+1) &
               +a *b *c *green(idx+1,idy+1,idz+1)
      else
         rinv = ONE/sqrt(dble(dx**2 + dy**2 + dz**2))
      end if
      pot = pot + qtmp*rinv
   end do  !  iatm = atmfirst, atmlast

end subroutine get_coulpot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign Debye-Huckel potential for the boundary charge grid
subroutine pb_bndcnd( bv, chgrd )

   implicit none
   
   ! Common variables
    
   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green
     
   ! Passed variables
    
   _REAL_ bv(xm,ym,zm), chgrd(xm,ym,zm)
    
   ! Local variables
    
   integer i, j, k, iatm, ii
   integer xmtmp, ymtmp, zmtmp, ix, iy, iz, itmp, jtmp, ktmp, idx, idy, idz
   _REAL_ htmp, goxtmp, goytmp, goztmp
   _REAL_ qtmp, factor
   _REAL_ x, y, z, dx, dy, dz, xtmp, ytmp, ztmp
   _REAL_ xi, yi, zi, aa, bb, cc, aa1, bb1, cc1
   _REAL_ r, rinv

   ! part a: level = 1 cases ::::: 
   ! bcopt = 2
   ! zero potential in the singularty-free PB.
   ! the boundary will be all solvent
    
   if ( level == 1 .and. bcopt == 2 ) then
    
   ! bcopt = 4 .or. bcopt = 6
   ! sum of atom charge deby-huckel potentials in the singular (4) or singularity-free (6) PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. ( bcopt == 4 .or. bcopt == 6 ) ) then

      do iatm = atmfirst, atmlast
         xtmp = gcrd(1,iatm); ytmp = gcrd(2,iatm); ztmp = gcrd(3,iatm)
         qtmp = INV_FOURPI*acrg(iatm)
          
         ! k=0 and k=zm+1 faces
          
         do j = 1, ym; do i = 1, xm
            dx = i-xtmp; dy = j-ytmp; dz = ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h*r))*qtmp/r
             
            dz = zm+1-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do
          
         ! j=0 and ym+1 faces
          
         do k = 1, zm; do i = 1, xm
            dx = i-xtmp; dy  = ytmp; dz  = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h*r))*qtmp/r
             
            dy = ym+1-ytmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do
          
         ! i=0 and i=xm+1 faces
          
         do k = 1, zm; do j = 1, ym
            dx = xtmp; dy = j-ytmp; dz = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h*r))*qtmp/r
             
            dx = xm+1-xtmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do
      end do  !  iatm = atmfirst, atmlast
    
   ! part b: level > 1 case 
   ! electrostatic focusing
    
   else if ( level > 1 ) then
      xmtmp  = savxm(level-1) ; ymtmp  = savym(level-1) ; zmtmp  = savzm(level-1)
      htmp   = savh(level-1)
      goxtmp = savgox(level-1); goytmp = savgoy(level-1); goztmp = savgoz(level-1)
       
      ! k=0 and k=zm+1 faces
       
      do j = 1, ym; do i = 1, xm
          
         x  = gox + h*i        ; y  = goy + h*j        ; z  = goz
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa  = xi - dble( ix ); bb  = yi - dble( iy ); cc  = zi - dble( iz )
         aa1 = ONE - aa       ; bb1 = ONE - bb       ; cc1 = ONE - cc
         bv(i,j,1 ) = bv(i,j,1 ) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         z  = goz + h*(zm+1)
         zi = (z - goztmp)/htmp
         iz = int( zi )
         cc  = zi - dble( iz )
         cc1 = ONE - cc
         bv(i,j,zm) = bv(i,j,zm) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
   
      ! j=0 and j=ym+1 faces
   
      do k = 1, zm; do i = 1, xm
          
         x  = gox + h*i        ; y  = goy              ; z  = goz + h*k
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa = xi - dble( ix ); bb = yi - dble( iy ); cc = zi - dble( iz )
         aa1 = ONE - aa      ; bb1 = ONE - bb      ; cc1 = ONE - cc
         bv(i,1 ,k) = bv(i,1 ,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         y  = goy + h*(ym+1)
         yi = (y - goytmp)/htmp
         iy = int( yi )
         bb  = yi - dble( iy )
         bb1 = ONE - bb
         bv(i,ym,k) = bv(i,ym,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
        
      ! i=0 and i=xm+1 faces
        
      do k = 1, zm; do j = 1, ym
          
         x  = gox              ; y  = goy + h*j        ; z  = goz + h*k
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa  = xi - dble( ix ); bb  = yi - dble( iy ); cc  = zi - dble( iz )
         aa1 = ONE - aa       ; bb1 = ONE - bb       ; cc1 = ONE - cc
         bv(1 ,j,k) = bv(1 ,j,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         x  = gox + h * (xm+1)
         xi = (x - goxtmp)/htmp
         ix = int( xi )
         aa  = xi - dble( ix )
         aa1 = ONE - aa
         bv(xm,j,k) = bv(xm,j,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
       
   else
       
      ! unknown bcopt
       
      write(6, *) 'PB bomb in pb_nhadrv(): unknown BC option', bcopt
      call mexit(6, 1)
   end if 
 
end subroutine pb_bndcnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ phi interpretation, xm, ym, zm are for the previous phi map
_REAL_ function phintp( xmtmp,ymtmp,zmtmp,ix,iy,iz,aa,bb,cc,aa1,bb1,cc1 )
    
   ! Passed variables
    
   integer, intent(in) :: xmtmp, ymtmp, zmtmp, ix, iy, iz
   _REAL_, intent(in) :: aa, bb, cc, aa1, bb1, cc1
    
   ! Local Variables
    
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
    
   ! determine the position of the point w.r.t. the map
    
   bb1cc1 = bb1*cc1; bb_cc1 = bb *cc1; bb1cc  = bb1*cc ; bb_cc  = bb *cc

   ! triliner interpolation
    
   phintp = aa1*bb1cc1*phi( ix   + xmtmp*( iy-1 + ymtmp*( iz-1 ) ) ) + &
            aa *bb1cc1*phi( ix+1 + xmtmp*( iy-1 + ymtmp*( iz-1 ) ) ) + &
            aa1*bb_cc1*phi( ix   + xmtmp*( iy   + ymtmp*( iz-1 ) ) ) + &
            aa *bb_cc1*phi( ix+1 + xmtmp*( iy   + ymtmp*( iz-1 ) ) ) + &
            aa1*bb1cc *phi( ix   + xmtmp*( iy-1 + ymtmp*( iz   ) ) ) + &
            aa *bb1cc *phi( ix+1 + xmtmp*( iy-1 + ymtmp*( iz   ) ) ) + &
            aa1*bb_cc *phi( ix   + xmtmp*( iy   + ymtmp*( iz   ) ) ) + &
            aa *bb_cc *phi( ix+1 + xmtmp*( iy   + ymtmp*( iz   ) ) )
                
end function phintp

end subroutine pb_nhadrv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ This is the real IIM driver 
subroutine iim(gox,goy,goz,l,m,n,lvlset,insas,nbnd,iepsav,epsin,epsout,bv,u,u0,h,&
               atmfirst,atmlast,bcopt,maxitn,itn,accept,norm,inorm,gcrd,natom,&
               ngrdcrg,grdcrg,qgrdcrg,chgrd,ipb,cphi&
              )

   ! xs <- gox, ys <- goy, zs <- goz
   ! l  <- xm , m  <- ym , n  <- zm
   ! bi <- epsin, bo <- epsout
   ! u  <- phi, u0 <- xs

   use pbtimer_module
   use density_surface, only : index, index2, x, y, z, cirreg, bndatmlst, bndatmptr
   use iim_util, only : wp, qp, qyp, qzp, wyp, wzp, wyyp, wzzp, wyzp
   use solvent_accessibility, only : dprob,radip3
   use poisson_boltzmann, only : acrg, acrd
  
   implicit none

   ! passed variables
      
   _REAL_  gox,goy,goz,h,epsout,epsin
   integer l,m,n,nbnd,atmfirst,atmlast,bcopt,natom
   integer maxitn,itn
   _REAL_ accept,norm,inorm
   integer insas(0:l+1,0:m+1,0:n+1),iepsav(1:4,1:l*m*n)
   _REAL_  bv(l,m,n),u0(l,m,n),u(l,m,n)
   _REAL_  lvlset(0:l+1,0:m+1,0:n+1)
   _REAL_ gcrd(3,natom)

   integer ngrdcrg,grdcrg(3,*),ipb
   _REAL_ qgrdcrg(*)
   _REAL_ chgrd(l,m,n), cphi(l,m,n)

   ! local variables

   _REAL_, parameter :: eps0 = 8.8542D-12 / (1.6022D-19)**2 /  &
                               (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer, parameter :: nq=27
   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ bi,bo
   _REAL_ beta_max,rhs
   _REAL_ err
   _REAL_ :: intersect, jump_condition
   integer :: nbndx, nbndy, nbndz
   
   _REAL_, allocatable :: c2(:, :)
   !character(10) str
   !integer :: ndot
   !ndot=0
   !open(unit=58,file='interface.dot')
   !write (58, *) nbnd
   !write (58, '("DOTS")')

   ! interfacing variables 

   bi = epsin; bo = epsout

   ! allocating working arrays

   call pbtimer_start(PBTIME_PBBUILDSYS)

   allocate (c2(nbnd, 7)) ! A matrix coefficients

   beta_max = max(bi,bo)
   nz_num = 0
   c2 = ZERO
   rhs = ZERO

   if (ipb == 7) then

   ! new_jump_condition HA method, aka X factor method

   call X_factor_HA( c2, bv, bi, bo )

   else if (ipb == 8) then

   ! Second order HA method

   call second_order_HA( c2, bv )
   
   end if

   !close(58)

   ! 7-band nonzero conefficients in A

   nz_num = 7*l*m*n-2*(1+l+l*m)

   call pbtimer_stop(PBTIME_PBBUILDSYS)

   ! entering the linear system solver

   call pbtimer_start(PBTIME_PBSOLV)
#ifdef nvCUDA
   call bicg_dia(l,m,n,nbnd,epsin,epsout,h,c2,index,index2,bv,u,u0,maxitn,itn,accept,inorm,norm)
#else
   call bicg(l,m,n,nbnd,nz_num,epsin,epsout,h,c2,index,index2,bv,u,u0,maxitn,itn,accept,inorm,norm)
#endif
   call pbtimer_stop(PBTIME_PBSOLV)

   ! converting back to the PBSA unit for potential
 
   do k =1, n
   do j =1, m
   do i =1, l
      u(i,j,k) = u(i,j,k) * INV_FOURPI / eps0
   end do
   end do
   end do
   
   deallocate (c2)

contains

subroutine second_order_HA( c2, bv )

   implicit none
   _REAL_ :: c2(nbnd,7), bv(l,m,n)

   integer :: i, j, k, ir

   do ir = 1, nbnd
      i = iepsav(1,ir)
      j = iepsav(2,ir)
      k = iepsav(3,ir)

      ! x direction

      if ((lvlset(i,j,k)*lvlset(i-1,j,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i+1,j,k) > 0d0)) then
         intersect=find_intersect(i,j,k,i-1,j,k,ir)
         if (lvlset(i-1,j,k) < 0d0) then
            c2(ir,1) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,7) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else
            c2(ir,1) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,7) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if ((lvlset(i,j,k)*lvlset(i+1,j,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i-1,j,k) > 0d0)) then
         intersect=find_intersect(i,j,k,i+1,j,k,ir)
         if (lvlset(i+1,j,k) < 0d0) then
            c2(ir,1) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,7) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else
            c2(ir,1) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,7) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if ((lvlset(i,j,k)*lvlset(i+1,j,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i-1,j,k) < 0d0)) then
         if (lvlset(i,j,k) > 0d0) then
            intersect=find_intersect(i,j,k,i-1,j,k,ir)
            c2(ir,1) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i+1,j,k,ir)
            c2(ir,7) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            intersect=find_intersect(i,j,k,i-1,j,k,ir)
            c2(ir,1) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            intersect=find_intersect(i,j,k,i+1,j,k,ir)
            c2(ir,7) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k) == 0d0) then
         intersect=0d0
         if ((lvlset(i-1,j,k) > 0d0) .and. (lvlset(i+1,j,k) <= 0d0)) then
            c2(ir,1) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,7) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i-1,j,k) <= 0d0) .and. (lvlset(i+1,j,k) > 0d0)) then
            c2(ir,1) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,7) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i-1,j,k) > 0d0) .and. (lvlset(i+1,j,k) > 0d0)) then
            c2(ir,1) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,7) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i-1,j,k) == 0d0) then
         if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i+1,j,k) > 0d0)) then
            intersect=1d0
            c2(ir,1) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,7) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i+1,j,k) <= 0d0)) then
            intersect=1d0
            c2(ir,1) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i+1,j,k,ir)
            c2(ir,7) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i+1,j,k) > 0d0)) then
            intersect=find_intersect(i,j,k,i+1,j,k,ir)
            c2(ir,1) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,7) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if (lvlset(i+1,j,k) == 0d0) then
         if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i-1,j,k) > 0d0)) then
            intersect=find_intersect(i,j,k,i-1,j,k,ir)
            c2(ir,1) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,7) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i-1,j,k) > 0d0)) then
            intersect=1d0
            c2(ir,1) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,7) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i-1,j,k) < 0d0)) then
            intersect=find_intersect(i,j,k,i-1,j,k,ir)
            c2(ir,1) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=1d0
            c2(ir,7) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         end if
      else if (lvlset(i,j,k)<=0 .and. lvlset(i-1,j,k)<=0 .and. lvlset(i+1,j,k)<=0) then
         c2(ir,1) = bi
         c2(ir,7) = bi
      else if (lvlset(i,j,k)>0 .and. lvlset(i-1,j,k)>0 .and. lvlset(i+1,j,k)>0) then
         c2(ir,1) = bo
         c2(ir,7) = bo
      else
         write(6,*) 'PBSA BOMB in pb_nhadrv: 2nd-HA fails to set matrix in x-direction'
         call mexit(6,1)
      end if
      
      ! y direction

      if ((lvlset(i,j,k)*lvlset(i,j-1,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j+1,k) > 0d0)) then
         intersect=find_intersect(i,j,k,i,j-1,k,ir)
         if (lvlset(i,j-1,k) < 0d0) then
            c2(ir,2) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,6) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else
            c2(ir,2) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,6) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if ((lvlset(i,j,k)*lvlset(i,j+1,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j-1,k) > 0d0)) then
         intersect=find_intersect(i,j,k,i,j+1,k,ir)
         if (lvlset(i,j+1,k) < 0d0) then
            c2(ir,2) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,6) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else
            c2(ir,2) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,6) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if ((lvlset(i,j,k)*lvlset(i,j+1,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j-1,k) < 0d0)) then
         if (lvlset(i,j,k) > 0d0) then
            intersect=find_intersect(i,j,k,i,j-1,k,ir)
            c2(ir,2) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i,j+1,k,ir)
            c2(ir,6) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            intersect=find_intersect(i,j,k,i,j-1,k,ir)
            c2(ir,2) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            intersect=find_intersect(i,j,k,i,j+1,k,ir)
            c2(ir,6) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k) == 0d0) then
         intersect=0d0
         if ((lvlset(i,j-1,k) > 0d0) .and. (lvlset(i,j+1,k) <= 0d0)) then
            c2(ir,2) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,6) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i,j-1,k) <= 0d0) .and. (lvlset(i,j+1,k) > 0d0)) then
            c2(ir,2) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,6) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i,j-1,k) > 0d0) .and. (lvlset(i,j+1,k) > 0d0)) then
            c2(ir,2) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,6) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j-1,k) == 0d0) then
         if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j+1,k) > 0d0)) then
            intersect=1d0
            c2(ir,2) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,6) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j+1,k) <= 0d0)) then
            intersect=1d0
            c2(ir,2) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i,j+1,k,ir)
            c2(ir,6) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i,j+1,k) > 0d0)) then
            intersect=find_intersect(i,j,k,i,j+1,k,ir)
            c2(ir,2) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,6) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if (lvlset(i,j+1,k) == 0d0) then
         if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i,j-1,k) > 0d0)) then
            intersect=find_intersect(i,j,k,i,j-1,k,ir)
            c2(ir,2) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,6) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j-1,k) > 0d0)) then
            intersect=1d0
            c2(ir,2) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,6) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j-1,k) < 0d0)) then
            intersect=find_intersect(i,j,k,i,j-1,k,ir)
            c2(ir,2) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=1d0
            c2(ir,6) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         end if
      else if (lvlset(i,j,k)<=0 .and. lvlset(i,j-1,k)<=0 .and. lvlset(i,j+1,k)<=0) then
         c2(ir,2) = bi
         c2(ir,6) = bi
      else if (lvlset(i,j,k)>0 .and. lvlset(i,j-1,k)>0 .and. lvlset(i,j+1,k)>0) then
         c2(ir,2) = bo
         c2(ir,6) = bo
      else
         write(6,*) 'PBSA BOMB in pb_nhadrv: 2nd-HA fails to set matrix in y-direction'
         call mexit(6,1)
      end if
     
      ! z direction

      if ((lvlset(i,j,k)*lvlset(i,j,k-1) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j,k+1) > 0d0)) then
         intersect=find_intersect(i,j,k,i,j,k-1,ir)
         if (lvlset(i,j,k-1) < 0d0) then
            c2(ir,3) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,5) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else
            c2(ir,3) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,5) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if ((lvlset(i,j,k)*lvlset(i,j,k+1) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j,k-1) > 0d0)) then
         intersect=find_intersect(i,j,k,i,j,k+1,ir)
         if (lvlset(i,j,k+1) < 0d0) then
            c2(ir,3) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,5) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else
            c2(ir,3) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,5) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if ((lvlset(i,j,k)*lvlset(i,j,k+1) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j,k-1) < 0d0)) then
         if (lvlset(i,j,k) > 0d0) then
            intersect=find_intersect(i,j,k,i,j,k-1,ir)
            c2(ir,3) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i,j,k+1,ir)
            c2(ir,5) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            intersect=find_intersect(i,j,k,i,j,k-1,ir)
            c2(ir,3) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            intersect=find_intersect(i,j,k,i,j,k+1,ir)
            c2(ir,5) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k) == 0d0) then
         intersect=0d0
         if ((lvlset(i,j,k-1) > 0d0) .and. (lvlset(i,j,k+1) <= 0d0)) then
            c2(ir,3) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,5) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i,j,k-1) <= 0d0) .and. (lvlset(i,j,k+1) > 0d0)) then
            c2(ir,3) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,5) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i,j,k-1) > 0d0) .and. (lvlset(i,j,k+1) > 0d0)) then
            c2(ir,3) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,5) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k-1) == 0d0) then
         if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j,k+1) > 0d0)) then
            intersect=1d0
            c2(ir,3) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,5) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j,k+1) <= 0d0)) then
            intersect=1d0
            c2(ir,3) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i,j,k+1,ir)
            c2(ir,5) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i,j,k+1) > 0d0)) then
            intersect=find_intersect(i,j,k,i,j,k+1,ir)
            c2(ir,3) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,5) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         end if
      else if (lvlset(i,j,k+1) == 0d0) then
         if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i,j,k-1) > 0d0)) then
            intersect=find_intersect(i,j,k,i,j,k-1,ir)
            c2(ir,3) = 2d0*bi*1d0/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
            c2(ir,5) = 2d0*bi*(intersect+(1d0-intersect)*(bi/bo))/&
                       (intersect*(1d0+intersect)+(bi/bo)*(1d0-intersect)*(2d0+intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j,k-1) > 0d0)) then
            intersect=1d0
            c2(ir,3) = 2d0*bi*((bi/bo)*intersect+(1d0-intersect))/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
            c2(ir,5) = 2d0*bi*(bi/bo)/&
                       ((bi/bo)**2*(1d0+intersect)*intersect+(bi/bo)*(2d0+intersect)*(1d0-intersect))
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j,k-1) < 0d0)) then
            intersect=find_intersect(i,j,k,i,j,k-1,ir)
            c2(ir,3) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=1d0
            c2(ir,5) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         end if
      else if (lvlset(i,j,k)<=0 .and. lvlset(i,j,k-1)<=0 .and. lvlset(i,j,k+1)<=0) then
         c2(ir,3) = bi
         c2(ir,5) = bi
      else if (lvlset(i,j,k)>0 .and. lvlset(i,j,k-1)>0 .and. lvlset(i,j,k+1)>0) then
         c2(ir,3) = bo
         c2(ir,5) = bo
      else
         write(6,*) 'PBSA BOMB in pb_nhadrv: 2nd-HA fails to set matrix in z-direction'
         call mexit(6,1)
      end if
     
      c2(ir,4) = -(c2(ir,1)+c2(ir,2)+c2(ir,3)+c2(ir,5)+c2(ir,6)+c2(ir,7))
      
      !if (lvlset(i-1,j,k) <= 0d0) rhs = rhs - c2(ir,1)*coulomb_potential(dble(i)-1d0,dble(j),dble(k))
      !if (lvlset(i,j-1,k) <= 0d0) rhs = rhs - c2(ir,2)*coulomb_potential(dble(i),dble(j)-1d0,dble(k))
      !if (lvlset(i,j,k-1) <= 0d0) rhs = rhs - c2(ir,3)*coulomb_potential(dble(i),dble(j),dble(k)-1d0)
      !if (lvlset(i,j,k) <= 0d0) rhs = rhs - c2(ir,4)*coulomb_potential(dble(i),dble(j),dble(k))
      !if (lvlset(i,j,k+1) <= 0d0) rhs = rhs - c2(ir,5)*coulomb_potential(dble(i),dble(j),dble(k)+1d0)
      !if (lvlset(i,j+1,k) <= 0d0) rhs = rhs - c2(ir,6)*coulomb_potential(dble(i),dble(j)+1d0,dble(k))
      !if (lvlset(i+1,j,k) <= 0d0) rhs = rhs - c2(ir,7)*coulomb_potential(dble(i)+1d0,dble(j),dble(k))

      if (lvlset(i-1,j,k) <= 0d0) rhs = rhs - c2(ir,1)*cphi(i-1,j,k)*FOURPI*eps0
      if (lvlset(i,j-1,k) <= 0d0) rhs = rhs - c2(ir,2)*cphi(i,j-1,k)*FOURPI*eps0
      if (lvlset(i,j,k-1) <= 0d0) rhs = rhs - c2(ir,3)*cphi(i,j,k-1)*FOURPI*eps0
      if (lvlset(i,j,k) <= 0d0) rhs = rhs - c2(ir,4)*cphi(i,j,k)*FOURPI*eps0
      if (lvlset(i,j,k+1) <= 0d0) rhs = rhs - c2(ir,5)*cphi(i,j,k+1)*FOURPI*eps0
      if (lvlset(i,j+1,k) <= 0d0) rhs = rhs - c2(ir,6)*cphi(i,j+1,k)*FOURPI*eps0
      if (lvlset(i+1,j,k) <= 0d0) rhs = rhs - c2(ir,7)*cphi(i+1,j,k)*FOURPI*eps0

      rhs = rhs - chgrd(i,j,k)*FOURPI/h
      bv(i,j,k) = rhs*h*h*h/FOURPI
      rhs = ZERO
   end do

end subroutine

subroutine X_factor_HA( c2, bv, bi, bo )

   implicit none
   _REAL_ :: c2(nbnd,7), bv(l,m,n)
   _REAL_ :: bi, bo

   integer :: i, j, k, ir

   do ir = 1, nbnd
      i = iepsav(1,ir)
      j = iepsav(2,ir)
      k = iepsav(3,ir)

      ! x direction

      if ((lvlset(i,j,k)*lvlset(i-1,j,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i+1,j,k) > 0d0)) then
         intersect=find_intersect(i,j,k,i-1,j,k,ir)
         call new_jump_condition(dble(i)-intersect,dble(j),dble(k),ir,1,jump_condition)
         if (lvlset(i-1,j,k) < 0d0) then
            bi = bo * jump_condition
            c2(ir,1) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            c2(ir,7) = bo
         else
            bo = bi / jump_condition
            c2(ir,1) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,7) = bi
         end if
      else if ((lvlset(i,j,k)*lvlset(i+1,j,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i-1,j,k) > 0d0)) then
         intersect=find_intersect(i,j,k,i+1,j,k,ir)
         call new_jump_condition(dble(i)+intersect,dble(j),dble(k),ir,1,jump_condition)
         if (lvlset(i+1,j,k) < 0d0) then
            bi = bo * jump_condition
            c2(ir,1) = bo
            c2(ir,7) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            bo = bi / jump_condition
            c2(ir,1) = bi
            c2(ir,7) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if ((lvlset(i,j,k)*lvlset(i+1,j,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i-1,j,k) < 0d0)) then
         if (lvlset(i,j,k) > 0d0) then
            intersect=find_intersect(i,j,k,i-1,j,k,ir)
            call new_jump_condition(dble(i)-intersect,dble(j),dble(k),ir,1,jump_condition)
            bi = bo * jump_condition
            c2(ir,1) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i+1,j,k,ir)
            call new_jump_condition(dble(i)+intersect,dble(j),dble(k),ir,1,jump_condition)
            bi = bo * jump_condition
            c2(ir,7) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            intersect=find_intersect(i,j,k,i-1,j,k,ir)
            call new_jump_condition(dble(i)-intersect,dble(j),dble(k),ir,1,jump_condition)
            bo = bi / jump_condition
            c2(ir,1) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            intersect=find_intersect(i,j,k,i+1,j,k,ir)
            call new_jump_condition(dble(i)+intersect,dble(j),dble(k),ir,1,jump_condition)
            bo = bi / jump_condition
            c2(ir,7) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k) == 0d0) then
         intersect=0d0
         call new_jump_condition(dble(i),dble(j),dble(k),ir,1,jump_condition)
         bo = bi / jump_condition
         if ((lvlset(i-1,j,k) > 0d0) .and. (lvlset(i+1,j,k) <= 0d0)) then
            c2(ir,1) = 1d0/((1d0-intersect)/bo+intersect/bi)
            c2(ir,7) = bi
         else if ((lvlset(i-1,j,k) <= 0d0) .and. (lvlset(i+1,j,k) > 0d0)) then
            c2(ir,1) = bi
            c2(ir,7) = 1d0/((1d0-intersect)/bo+intersect/bi)
         else if ((lvlset(i-1,j,k) > 0d0) .and. (lvlset(i+1,j,k) > 0d0)) then
            c2(ir,1) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,7) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i-1,j,k) == 0d0) then
         if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i+1,j,k) > 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i)-intersect,dble(j),dble(k),ir,1,jump_condition)
            bi = bo * jump_condition
            c2(ir,1) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            c2(ir,7) = bo
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i+1,j,k) <= 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i)-intersect,dble(j),dble(k),ir,1,jump_condition)
            bi = bo * jump_condition
            c2(ir,1) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i+1,j,k,ir)
            call new_jump_condition(dble(i)+intersect,dble(j),dble(k),ir,1,jump_condition)
            bi = bo * jump_condition
            c2(ir,7) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i+1,j,k) > 0d0)) then
            intersect=find_intersect(i,j,k,i+1,j,k,ir)
            call new_jump_condition(dble(i)+intersect,dble(j),dble(k),ir,1,jump_condition)
            bo = bi / jump_condition
            c2(ir,1) = bi
            c2(ir,7) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i+1,j,k) == 0d0) then
         if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i-1,j,k) > 0d0)) then
            intersect=find_intersect(i,j,k,i-1,j,k,ir)
            call new_jump_condition(dble(i)-intersect,dble(j),dble(k),ir,1,jump_condition)
            bo = bi / jump_condition
            c2(ir,1) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,7) = bi
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i-1,j,k) > 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i)+intersect,dble(j),dble(k),ir,1,jump_condition)
            bi = bo * jump_condition
            c2(ir,1) = bo
            c2(ir,7) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i-1,j,k) < 0d0)) then
            intersect=find_intersect(i,j,k,i-1,j,k,ir)
            call new_jump_condition(dble(i)-intersect,dble(j),dble(k),ir,1,jump_condition)
            bi = bo * jump_condition
            c2(ir,1) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=1d0
            call new_jump_condition(dble(i)+intersect,dble(j),dble(k),ir,1,jump_condition)
            bi = bo * jump_condition
            c2(ir,7) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         end if
      else if (lvlset(i,j,k)<=0 .and. lvlset(i-1,j,k)<=0 .and. lvlset(i+1,j,k)<=0) then
         c2(ir,1) = bi
         c2(ir,7) = bi
      else if (lvlset(i,j,k)>0 .and. lvlset(i-1,j,k)>0 .and. lvlset(i+1,j,k)>0) then
         c2(ir,1) = bo
         c2(ir,7) = bo
      else
         write(6,*) 'PBSA BOMB in pb_nhadrv: X-HA fails to set matrix in x-direction'
         call mexit(6,1)
      end if
      bi = epsin
      bo = epsout
      
      ! y direction

      if ((lvlset(i,j,k)*lvlset(i,j-1,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j+1,k) > 0d0)) then
         intersect=find_intersect(i,j,k,i,j-1,k,ir)
         call new_jump_condition(dble(i),dble(j)-intersect,dble(k),ir,2,jump_condition)
         if (lvlset(i,j-1,k) < 0d0) then
            bi = bo * jump_condition
            c2(ir,2) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            c2(ir,6) = bo
         else
            bo = bi / jump_condition
            c2(ir,2) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,6) = bi
         end if
      else if ((lvlset(i,j,k)*lvlset(i,j+1,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j-1,k) > 0d0)) then
         intersect=find_intersect(i,j,k,i,j+1,k,ir)
         call new_jump_condition(dble(i),dble(j)+intersect,dble(k),ir,2,jump_condition)
         if (lvlset(i,j+1,k) < 0d0) then
            bi = bo * jump_condition
            c2(ir,2) = bo
            c2(ir,6) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            bo = bi / jump_condition
            c2(ir,2) = bi
            c2(ir,6) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if ((lvlset(i,j,k)*lvlset(i,j+1,k) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j-1,k) < 0d0)) then
         if (lvlset(i,j,k) > 0d0) then
            intersect=find_intersect(i,j,k,i,j-1,k,ir)
            call new_jump_condition(dble(i),dble(j)-intersect,dble(k),ir,2,jump_condition)
            bi = bo * jump_condition
            c2(ir,2) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i,j+1,k,ir)
            call new_jump_condition(dble(i),dble(j)+intersect,dble(k),ir,2,jump_condition)
            bi = bo * jump_condition
            c2(ir,6) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            intersect=find_intersect(i,j,k,i,j-1,k,ir)
            call new_jump_condition(dble(i),dble(j)-intersect,dble(k),ir,2,jump_condition)
            bo = bi / jump_condition
            c2(ir,2) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            intersect=find_intersect(i,j,k,i,j+1,k,ir)
            call new_jump_condition(dble(i),dble(j)+intersect,dble(k),ir,2,jump_condition)
            bo = bi / jump_condition
            c2(ir,6) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k) == 0d0) then
         intersect=0d0
         call new_jump_condition(dble(i),dble(j),dble(k),ir,2,jump_condition)
         bo = bi / jump_condition
         if ((lvlset(i,j-1,k) > 0d0) .and. (lvlset(i,j+1,k) <= 0d0)) then
            c2(ir,2) = 1d0/((1d0-intersect)/bo+intersect/bi)
            c2(ir,6) = bi
         else if ((lvlset(i,j-1,k) <= 0d0) .and. (lvlset(i,j+1,k) > 0d0)) then
            c2(ir,2) = bi
            c2(ir,6) = 1d0/((1d0-intersect)/bo+intersect/bi)
         else if ((lvlset(i,j-1,k) > 0d0) .and. (lvlset(i,j+1,k) > 0d0)) then
            c2(ir,2) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,6) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j-1,k) == 0d0) then
         if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j+1,k) > 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i),dble(j)-intersect,dble(k),ir,2,jump_condition)
            bi = bo * jump_condition
            c2(ir,2) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            c2(ir,6) = bo
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j+1,k) <= 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i),dble(j)-intersect,dble(k),ir,2,jump_condition)
            bi = bo * jump_condition
            c2(ir,2) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i,j+1,k,ir)
            call new_jump_condition(dble(i),dble(j)+intersect,dble(k),ir,2,jump_condition)
            bi = bo * jump_condition
            c2(ir,6) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i,j+1,k) > 0d0)) then
            intersect=find_intersect(i,j,k,i,j+1,k,ir)
            call new_jump_condition(dble(i),dble(j)+intersect,dble(k),ir,2,jump_condition)
            bo = bi / jump_condition
            c2(ir,2) = bi
            c2(ir,6) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j+1,k) == 0d0) then
         if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i,j-1,k) > 0d0)) then
            intersect=find_intersect(i,j,k,i,j-1,k,ir)
            call new_jump_condition(dble(i),dble(j)-intersect,dble(k),ir,2,jump_condition)
            bo = bi / jump_condition
            c2(ir,2) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,6) = bi
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j-1,k) > 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i),dble(j)+intersect,dble(k),ir,2,jump_condition)
            bi = bo * jump_condition
            c2(ir,2) = bo
            c2(ir,6) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j-1,k) < 0d0)) then
            intersect=find_intersect(i,j,k,i,j-1,k,ir)
            call new_jump_condition(dble(i),dble(j)-intersect,dble(k),ir,2,jump_condition)
            bi = bo * jump_condition
            c2(ir,2) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=1d0
            call new_jump_condition(dble(i),dble(j)+intersect,dble(k),ir,2,jump_condition)
            bi = bo * jump_condition
            c2(ir,6) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         end if
      else if (lvlset(i,j,k)<=0 .and. lvlset(i,j-1,k)<=0 .and. lvlset(i,j+1,k)<=0) then
         c2(ir,2) = bi
         c2(ir,6) = bi
      else if (lvlset(i,j,k)>0 .and. lvlset(i,j-1,k)>0 .and. lvlset(i,j+1,k)>0) then
         c2(ir,2) = bo
         c2(ir,6) = bo
      else
         write(6,*) 'PBSA BOMB in pb_nhadrv: X-HA fails to set matrix in y-direction'
         call mexit(6,1)
      end if
      bi = epsin
      bo = epsout
    
      ! z direction

      if ((lvlset(i,j,k)*lvlset(i,j,k-1) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j,k+1) > 0d0)) then
         intersect=find_intersect(i,j,k,i,j,k-1,ir)
         call new_jump_condition(dble(i),dble(j),dble(k)-intersect,ir,3,jump_condition)
         if (lvlset(i,j,k-1) < 0d0) then
            bi = bo * jump_condition
            c2(ir,3) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            c2(ir,5) = bo
         else
            bo = bi / jump_condition
            c2(ir,3) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,5) = bi
         end if
      else if ((lvlset(i,j,k)*lvlset(i,j,k+1) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j,k-1) > 0d0)) then
         intersect=find_intersect(i,j,k,i,j,k+1,ir)
         call new_jump_condition(dble(i),dble(j),dble(k)+intersect,ir,3,jump_condition)
         if (lvlset(i,j,k+1) < 0d0) then
            bi = bo * jump_condition
            c2(ir,3) = bo
            c2(ir,5) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            bo = bi / jump_condition
            c2(ir,3) = bi
            c2(ir,5) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if ((lvlset(i,j,k)*lvlset(i,j,k+1) < 0d0) .and. (lvlset(i,j,k)*lvlset(i,j,k-1) < 0d0)) then
         if (lvlset(i,j,k) > 0d0) then
            intersect=find_intersect(i,j,k,i,j,k-1,ir)
            call new_jump_condition(dble(i),dble(j),dble(k)-intersect,ir,3,jump_condition)
            bi = bo * jump_condition
            c2(ir,3) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i,j,k+1,ir)
            call new_jump_condition(dble(i),dble(j),dble(k)+intersect,ir,3,jump_condition)
            bi = bo * jump_condition
            c2(ir,5) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else
            intersect=find_intersect(i,j,k,i,j,k-1,ir)
            call new_jump_condition(dble(i),dble(j),dble(k)-intersect,ir,3,jump_condition)
            bo = bi / jump_condition
            c2(ir,3) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            intersect=find_intersect(i,j,k,i,j,k+1,ir)
            call new_jump_condition(dble(i),dble(j),dble(k)+intersect,ir,3,jump_condition)
            bo = bi / jump_condition
            c2(ir,5) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k) == 0d0) then
         intersect=0d0
         call new_jump_condition(dble(i),dble(j),dble(k),ir,3,jump_condition)
         bo = bi / jump_condition
         if ((lvlset(i,j,k-1) > 0d0) .and. (lvlset(i,j,k+1) <= 0d0)) then
            c2(ir,3) = 1d0/((1d0-intersect)/bo+intersect/bi)
            c2(ir,5) = bi
         else if ((lvlset(i,j,k-1) <= 0d0) .and. (lvlset(i,j,k+1) > 0d0)) then
            c2(ir,3) = bi
            c2(ir,5) = 1d0/((1d0-intersect)/bo+intersect/bi)
         else if ((lvlset(i,j,k-1) > 0d0) .and. (lvlset(i,j,k+1) > 0d0)) then
            c2(ir,3) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,5) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k-1) == 0d0) then
         if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j,k+1) > 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i),dble(j),dble(k)-intersect,ir,3,jump_condition)
            bi = bo * jump_condition
            c2(ir,3) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            c2(ir,5) = bo
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j,k+1) <= 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i),dble(j),dble(k)-intersect,ir,3,jump_condition)
            bi = bo * jump_condition
            c2(ir,3) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=find_intersect(i,j,k,i,j,k+1,ir)
            call new_jump_condition(dble(i),dble(j),dble(k)+intersect,ir,3,jump_condition)
            bi = bo * jump_condition
            c2(ir,5) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i,j,k+1) > 0d0)) then
            intersect=find_intersect(i,j,k,i,j,k+1,ir)
            call new_jump_condition(dble(i),dble(j),dble(k)+intersect,ir,3,jump_condition)
            bo = bi / jump_condition
            c2(ir,3) = bi
            c2(ir,5) = 1d0/(intersect/bi+(1d0-intersect)/bo)
         end if
      else if (lvlset(i,j,k+1) == 0d0) then
         if ((lvlset(i,j,k) < 0d0) .and. (lvlset(i,j,k-1) > 0d0)) then
            intersect=find_intersect(i,j,k,i,j,k-1,ir)
            call new_jump_condition(dble(i),dble(j),dble(k)-intersect,ir,3,jump_condition)
            bo = bi / jump_condition
            c2(ir,3) = 1d0/(intersect/bi+(1d0-intersect)/bo)
            c2(ir,5) = bi
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j,k-1) > 0d0)) then
            intersect=1d0
            call new_jump_condition(dble(i),dble(j),dble(k)+intersect,ir,3,jump_condition)
            bi = bo * jump_condition
            c2(ir,3) = bo
            c2(ir,5) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         else if ((lvlset(i,j,k) > 0d0) .and. (lvlset(i,j,k-1) < 0d0)) then
            intersect=find_intersect(i,j,k,i,j,k-1,ir)
            call new_jump_condition(dble(i),dble(j),dble(k)-intersect,ir,3,jump_condition)
            bi = bo * jump_condition
            c2(ir,3) = 1d0/(intersect/bo+(1d0-intersect)/bi)
            intersect=1d0
            call new_jump_condition(dble(i),dble(j),dble(k)+intersect,ir,3,jump_condition)
            bi = bo * jump_condition
            c2(ir,5) = 1d0/(intersect/bo+(1d0-intersect)/bi)
         end if
      else if (lvlset(i,j,k)<=0 .and. lvlset(i,j,k-1)<=0 .and. lvlset(i,j,k+1)<=0) then
         c2(ir,3) = bi
         c2(ir,5) = bi
      else if (lvlset(i,j,k)>0 .and. lvlset(i,j,k-1)>0 .and. lvlset(i,j,k+1)>0) then
         c2(ir,3) = bo
         c2(ir,5) = bo
      else
         write(6,*) 'PBSA BOMB in pb_nhadrv: X-HA fails to set matrix in z-direction'
         call mexit(6,1)
      end if
     
      c2(ir,4) = -(c2(ir,1)+c2(ir,2)+c2(ir,3)+c2(ir,5)+c2(ir,6)+c2(ir,7))
      
      !if (lvlset(i-1,j,k) <= 0d0) rhs = rhs - c2(ir,1)*coulomb_potential(dble(i)-1d0,dble(j),dble(k))
      !if (lvlset(i,j-1,k) <= 0d0) rhs = rhs - c2(ir,2)*coulomb_potential(dble(i),dble(j)-1d0,dble(k))
      !if (lvlset(i,j,k-1) <= 0d0) rhs = rhs - c2(ir,3)*coulomb_potential(dble(i),dble(j),dble(k)-1d0)
      !if (lvlset(i,j,k) <= 0d0) rhs = rhs - c2(ir,4)*coulomb_potential(dble(i),dble(j),dble(k))
      !if (lvlset(i,j,k+1) <= 0d0) rhs = rhs - c2(ir,5)*coulomb_potential(dble(i),dble(j),dble(k)+1d0)
      !if (lvlset(i,j+1,k) <= 0d0) rhs = rhs - c2(ir,6)*coulomb_potential(dble(i),dble(j)+1d0,dble(k))
      !if (lvlset(i+1,j,k) <= 0d0) rhs = rhs - c2(ir,7)*coulomb_potential(dble(i)+1d0,dble(j),dble(k))

      if (lvlset(i-1,j,k) <= 0d0) rhs = rhs - c2(ir,1)*cphi(i-1,j,k)*FOURPI*eps0
      if (lvlset(i,j-1,k) <= 0d0) rhs = rhs - c2(ir,2)*cphi(i,j-1,k)*FOURPI*eps0
      if (lvlset(i,j,k-1) <= 0d0) rhs = rhs - c2(ir,3)*cphi(i,j,k-1)*FOURPI*eps0
      if (lvlset(i,j,k) <= 0d0) rhs = rhs - c2(ir,4)*cphi(i,j,k)*FOURPI*eps0
      if (lvlset(i,j,k+1) <= 0d0) rhs = rhs - c2(ir,5)*cphi(i,j,k+1)*FOURPI*eps0
      if (lvlset(i,j+1,k) <= 0d0) rhs = rhs - c2(ir,6)*cphi(i,j+1,k)*FOURPI*eps0
      if (lvlset(i+1,j,k) <= 0d0) rhs = rhs - c2(ir,7)*cphi(i+1,j,k)*FOURPI*eps0

      rhs = rhs - chgrd(i,j,k)*FOURPI/h
      bv(i,j,k) = rhs*h*h*h/FOURPI

      bi = epsin
      bo = epsout
      rhs = ZERO
   end do

end subroutine

subroutine new_jump_condition(x, y, z, l, direction, jump_condition)

   implicit none

   _REAL_ x, y, z, jump_condition
   integer l, direction

   logical, external :: disnan

   integer :: k
   _REAL_ :: distance_reverse_cubic
   _REAL_, dimension(3) :: normal, coulomb_field
   integer :: jp, jatm
   integer :: first, last
   _REAL_ :: dx, dy, dz, disth,factor

   coulomb_field = 0d0
   normal = 0d0

   factor=0.5/(dprob/h)
   first = bndatmptr(l-1)+1
   last = bndatmptr(l)
   do jp = first, last
      jatm = bndatmlst(jp)
      dx = x-gcrd(1,jatm)
      dy = y-gcrd(2,jatm)
      dz = z-gcrd(3,jatm)
      disth = sqrt(dx**2+dy**2+dz**2)
      dx = dx/disth*factor
      dy = dy/disth*factor
      dz = dz/disth*factor
      normal(1:3) = normal(1:3) + density_1deriv_atom((disth-radip3(jatm)/h)*factor,dx,dy,dz)
   end do
   normal(1:3) = normal(1:3)/sqrt(sum(normal(1:3)**2))

   x = x*h + gox
   y = y*h + goy
   z = z*h + goz
   do k = atmfirst, atmlast
      distance_reverse_cubic = 1d0/sqrt((acrd(1,k)-x)**2+(acrd(2,k)-y)**2+(acrd(3,k)-z)**2)**3
      coulomb_field(1) = coulomb_field(1) + (x-acrd(1,k))*acrg(k)*distance_reverse_cubic
      coulomb_field(2) = coulomb_field(2) + (y-acrd(2,k))*acrg(k)*distance_reverse_cubic
      coulomb_field(3) = coulomb_field(3) + (z-acrd(3,k))*acrg(k)*distance_reverse_cubic
   end do

   jump_condition = ((epsin/epsout)*2d0*(dot_product(normal,coulomb_field))*normal(direction) +&
                     (epsin/epsout)*2d0*(coulomb_field(direction)-(dot_product(normal,coulomb_field))*normal(direction))&
                    )/&
                    (2d0*(dot_product(normal,coulomb_field))*normal(direction) +&
                     (epsin/epsout)*2d0*(coulomb_field(direction)-(dot_product(normal,coulomb_field))*normal(direction))&
                    )

   if ( disnan(jump_condition) ) jump_condition = ONE

   if ( jump_condition > ONE) jump_condition = ONE
   ! when there are atoms nearby, need additional safeguaring for the computed X values.
   if ( last-first > 0 ) then
      if ( jump_condition < epsin/epsout ) jump_condition = epsin/epsout
   else
      if ( jump_condition < accept ) jump_condition = accept
   end if

end subroutine new_jump_condition

function density_1deriv_atom(dist,dx,dy,dz)
   implicit none

   _REAL_, dimension(3) :: density_1deriv_atom
   _REAL_ :: dist, dx, dy, dz

   integer m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))

   ! All dist and components are in the unit of solvent probe diameter

   density_1deriv_atom(1:3) = ZERO
   if ( dist > ONE ) then
   else if ( dist <= ZERO ) then
      density_1deriv_atom(1) = -4.527143d0 * dx
      density_1deriv_atom(2) = -4.527143d0 * dy
      density_1deriv_atom(3) = -4.527143d0 * dz
   else
      do m = 1, 5
         if ( dist > dash(m) .and. dist <= dash(m+1) ) then
            density_1deriv_atom(1) = spcoef(m,2)                  *dx + &
                               TWO  *spcoef(m,3)*(dist-dash(m))   *dx + &
                               THREE*spcoef(m,4)*(dist-dash(m))**2*dx
            density_1deriv_atom(2) = spcoef(m,2)                  *dy + &
                               TWO  *spcoef(m,3)*(dist-dash(m))   *dy + &
                               THREE*spcoef(m,4)*(dist-dash(m))**2*dy
            density_1deriv_atom(3) = spcoef(m,2)                  *dz + &
                               TWO  *spcoef(m,3)*(dist-dash(m))   *dz + &
                               THREE*spcoef(m,4)*(dist-dash(m))**2*dz
         end if
      end do
   end if

end function density_1deriv_atom

function find_intersect(i,j,k,i_end,j_end,k_end,l)

   implicit none

   integer :: i,j,k,i_end,j_end,k_end,l
   _REAL_ :: find_intersect
   _REAL_ :: t

   if (lvlset(i,j,k) == 0d0) then
      find_intersect = 0d0
      return
   end if
   if (lvlset(i_end,j_end,k_end) == 0d0) then
      find_intersect = 1d0
      return
   end if
   if (i /= i_end) then
      if (lvlset(i,j,k) < 0) call root(dble(min(i,i_end)),dble(max(i,i_end)),dble(2*(i-i_end)+i_end),&
                                       lvlset(min(i,i_end),j,k),lvlset(max(i,i_end),j,k),lvlset(2*(i-i_end)+i_end,j,k),t)
      if (lvlset(i,j,k) > 0) call root(dble(min(i,i_end)),dble(max(i,i_end)),dble(2*(i_end-i)+i),&
                                       lvlset(min(i,i_end),j,k),lvlset(max(i,i_end),j,k),lvlset(2*(i_end-i)+i,j,k),t)
      find_intersect = abs(t-dble(i))
   end if
   if (j /= j_end) then
      if (lvlset(i,j,k) < 0) call root(dble(min(j,j_end)),dble(max(j,j_end)),dble(2*(j-j_end)+j_end),&
                                       lvlset(i,min(j,j_end),k),lvlset(i,max(j,j_end),k),lvlset(i,2*(j-j_end)+j_end,k),t)
      if (lvlset(i,j,k) > 0) call root(dble(min(j,j_end)),dble(max(j,j_end)),dble(2*(j_end-j)+j),&
                                       lvlset(i,min(j,j_end),k),lvlset(i,max(j,j_end),k),lvlset(i,2*(j_end-j)+j,k),t)
      find_intersect = abs(t-dble(j))
   end if
   if (k /= k_end) then
      if (lvlset(i,j,k) < 0) call root(dble(min(k,k_end)),dble(max(k,k_end)),dble(2*(k-k_end)+k_end),&
                                       lvlset(i,j,min(k,k_end)),lvlset(i,j,max(k,k_end)),lvlset(i,j,2*(k-k_end)+k_end),t)
      if (lvlset(i,j,k) > 0) call root(dble(min(k,k_end)),dble(max(k,k_end)),dble(2*(k_end-k)+k),&
                                       lvlset(i,j,min(k,k_end)),lvlset(i,j,max(k,k_end)),lvlset(i,j,2*(k_end-k)+k),t)
      find_intersect = abs(t-dble(k))
   end if

end function find_intersect

subroutine root(x0,x1,x2,f0,f1,f2,t0)

   implicit none

   ! passed variables

   _REAL_ x0,x1,x2,f0,f1,f2,t0

   ! local variables

   _REAL_ b,c,a0,b0,c0,t,r1,r2

   b = (f0-f1)/(x0-x1)
   c = f2 - f1 - b*(x2-x1)
   c = c/( (x2-x0)*(x2-x1))

   a0 = c
   b0 = b - c*(x0+x1)
   c0 = f1 -b*x1 + c*x0*x1

   if ( a0 == 0 ) then
      t0 = -c0/b0
      return
   end if

   t = b0*b0 - 4.0d0*a0*c0

   ! If t <=0, must be double root t is close to zero

   if ( t <= 0.0d0 ) then
      t0 = -b0/(2.0d0*a0)
      return
   end if

   t = sqrt(t)
   if ( b0 >= 0.0d0 ) then
      r1 = (-b0-t)/(2.0d0*a0)
   else
      r1 = (-b0+t)/(2.0d0*a0)
   end if

   r2 = -b0/a0-r1

   if ( x0 <= r1 + 1.0d-7 .and. r1 <= x1+1.0d-7 ) then
      t0 = r1
   else
      t0 = r2
   end if

   if ( x0 > t0 ) t0 = x0
   if ( x1 < t0 ) t0 = x1

end subroutine root

function find_lvlset(x,y,z,l)
   implicit none
   
   _REAL_ :: find_lvlset
   _REAL_ :: x, y, z
   integer :: l

   integer :: jp, jatm
   integer :: atmfirst, atmlast 
   _REAL_ :: dx, dy, dz, disth,factor

   factor=0.5/(dprob/h)   
   find_lvlset = -ONE
   atmfirst = bndatmptr(l-1)+1
   atmlast = bndatmptr(l)
   do jp = atmfirst, atmlast
      jatm = bndatmlst(jp)
      dx = x-gcrd(1,jatm)
      dy = y-gcrd(2,jatm)
      dz = z-gcrd(3,jatm)
      disth = sqrt(dx**2+dy**2+dz**2)
      find_lvlset = find_lvlset + density_atom((disth-radip3(jatm)/h)*factor)
   end do

end function find_lvlset

function density_atom(dist)
   implicit none

   _REAL_ :: density_atom
   _REAL_ :: dist

   integer m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))

   ! All dist and components are in the unit of solvent probe diameter

   density_atom = ZERO
   if ( dist > ONE ) then
   else if ( dist <= ZERO ) then
      density_atom = 1.0d0 - 4.527143d0*dist
   else
      do m = 1, 5
         if ( dist > dash(m) .and. dist <= dash(m+1) ) then
            density_atom = spcoef(m,1)                   + &
                           spcoef(m,2)*(dist-dash(m))    + &
                           spcoef(m,3)*(dist-dash(m))**2 + &
                           spcoef(m,4)*(dist-dash(m))**3
         end if
      end do
   end if

end function density_atom

end subroutine iim

end module nha
