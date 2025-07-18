#include "../include/dprec.fh"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine irre31 here]
subroutine irre31(l,m,n,h,hx,hy,hz,ifail,i0,j0,k0,info,bmax,b_in,b_out,x,y,z,phi,index, &
              nq,maxirr,index2,cirreg,xx,rhs)

   ! set up the linear equation for an irregular point (i0,j0,k0)
   ! determin the coefficients in matrix A and the rhs to satisfy
   ! the specified jump conditions

   use iim_util
   implicit none

   ! parameter settings

   integer im,ime,immax,ns,in,inmax,imnn,lwar,liwar
   parameter(im=36, ime=10, immax=im, ns=27, in=ns, &
         inmax=in, imnn=im+2*in, &
         lwar=3*inmax*inmax/2+10*inmax+2*immax+1, &
         liwar=in)

   ! external functions

   _REAL_, external :: fmin

   ! interface variables

   integer ifail,i0,j0,k0,info,nq,maxirr
   integer index(l,m,n),index2(l,m,n)
   _REAL_  bmax,b_out,b_in
   _REAL_  rhs
   _REAL_  x(0:l+1), y(0:m+1), z(0:n+1)
   _REAL_  phi(0:l+1,0:m+1, 0:n+1)
   _REAL_  cirreg(15, maxirr)
   _REAL_  xx(in)

   ! local variables

   integer i,j,nc,kr,kctr,ndis,l,m,n,nirreg,k,ii,ir,iout,iprint,ji
   integer iwar(liwar)

   _REAL_  wy,wz,wyy,wzz,wyz,q,qy,qz,fkk_in,fkk_out &
          ,hxx,hyy,hzz,h1,eps,rho,rho1,bl2jmp,bl3jmp &
          ,hxyz,tmp,tmp1,fkk,corr &
          ,fkjmp,fjmp,uump,tmp2,temp,tmp3,x1,y1,z1  &
          ,xyy,xzz,xyz,h,hx,hy,hz,hmax,ujmp
   _REAL_  t(3,3),a(20)
   _REAL_  xloc(3),tmpx(3)
   _REAL_  bl_out(3),bl_in(3)
   _REAL_  cc(inmax,inmax), ca(immax,inmax), caa(immax,inmax)
   _REAL_  cd(inmax),cb(immax),xl(in),xu(in),uu(imnn)
   _REAL_  war(lwar)
   _REAL_  tmpcoe(7)

   ! retrieving info for the irregular point

   ir = index2(i0, j0, k0)

   x1     = cirreg( 1, ir)
   y1     = cirreg( 2, ir)
   z1     = cirreg( 3, ir)

   xyy    = cirreg( 4, ir)
   xzz    = cirreg( 5, ir)
   xyz    = cirreg( 6, ir)

   t(1,1) = cirreg( 7, ir)
   t(1,2) = cirreg( 8, ir)
   t(1,3) = cirreg( 9, ir)
   t(2,1) = cirreg(10, ir)
   t(2,2) = cirreg(11, ir)
   t(2,3) = cirreg(12, ir)
   t(3,1) = cirreg(13, ir)
   t(3,2) = cirreg(14, ir)
   t(3,3) = cirreg(15, ir)

   ! this is the jump of the rhs ...

   call fjmps(x1,y1,z1, fjmp)

   ! this is the jump of the linear term ...

   call fkjmps(x1,y1,z1, fkjmp)

   ! retrieve the jump in u, its first derivatives,
   ! and second derivatives

   ujmp = wp(ir)
   wy   = wyp(ir)
   wz   = wzp(ir)
   wyy  = wyyp(ir)
   wzz  = wzzp(ir)
   wyz  = wyzp(ir)

   ! retrieve the jump in beta u_n, its first derivative,
   ! and second derivatives

   q =  qp(ir)
   qy = qyp(ir)
   qz = qzp(ir)

   ! compute the derivatives of beta of inside and out
   ! and transform them into local coordinate system

   call betas(hx,hy,hz,x1,y1,z1,t,bl_in,bl_out,b_in,b_out)

   ! the linear term inside and outside

   fkk_in = 0.0d0!fk_in(x1,y1,z1)
   fkk_out = 0.0d0!fk_out(x1,y1,z1)

   hxx=hx*hx
   hyy=hy*hy
   hzz=hz*hz
   h1 = min(hxx,hyy,hzz)

   eps=0.1e-11

   ! set up ca() and cb() for the nearest grid points
   ! Jun: what are these?
   ! ca: distance/distance products
   ! cb: beta/beta derivatives

   do i=1,immax
      do j=1,inmax
         ca(i,j) = 0.0
      end do
      cb(i) = 0.0
   end do

   nc = 0
   kr = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1
            kr = kr +1
            if (i == i0 .and. j == j0 .and. k == k0) then
               kr = kr -1
               kctr = nc
            else
               ca(kr+ime,nc) = 1.0d0
               cb(kr+ime) = -1.0d-20
            end if
            xl(nc) = 0.0
            xu(nc) = 2.0*bmax/h1
         end do
      end do
   end do

   xl(kctr)=-12.0*bmax/h1
   xu(kctr)=0.0

   rho = b_in/b_out
   rho1 = 1 - rho
   bl2jmp = bl_out(2)-bl_in(2)
   bl3jmp = bl_out(3)-bl_in(3)
   temp = fkjmp/b_out

   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1

            ! ----------- local coordinates transformation

            tmpx(1) = x(i) - x1
            tmpx(2) = y(j) - y1
            tmpx(3) = z(k) - z1
            call matvec(3,3,t, tmpx, xloc)

            ! ----------- form the coefficients of equality constraints

            if (index(i,j,k) <= 3) then
               ca(1,nc)  = 1.0
               ca(2,nc)  = xloc(1)
               ca(3,nc)  = xloc(2)
               ca(4,nc)  = xloc(3)
               ca(5,nc)  = 0.5*xloc(1)*xloc(1)
               ca(6,nc)  = 0.5*xloc(2)*xloc(2)
               ca(7,nc)  = 0.5*xloc(3)*xloc(3)
               ca(8,nc)  = xloc(1)*xloc(2)
               ca(9,nc)  = xloc(1)*xloc(3)
               ca(10,nc) = xloc(2)*xloc(3)
            else
               ca(1,nc)  = 1.0-0.5*temp*xloc(1)*xloc(1)
               ca(2,nc)  = rho*xloc(1) &
                     - 0.5*xloc(1)*xloc(1)*(rho1*(xyy+xzz) &
                     - bl_in(1)/b_out + rho*bl_out(1)/b_out) &
                     + 0.5*xloc(2)*xloc(2)*xyy*rho1 &
                     + 0.5*xloc(3)*xloc(3)*xzz*rho1 &
                     + xloc(1)*xloc(2)*(bl_in(2)/b_out &
                     - rho*bl_out(2)/b_out) &
                     + xloc(1)*xloc(3)*(bl_in(3)/b_out &
                     - rho*bl_out(3)/b_out) &
                     + xloc(2)*xloc(3)*xyz*rho1
               ca(3,nc)  = - 0.5*xloc(1)*xloc(1)*bl2jmp/b_out &
                     + xloc(1)*xloc(2)*xyy*rho1 &
                     + xloc(1)*xloc(3)*xyz*rho1 &
                     + xloc(2)
               ca(4,nc)  = - 0.5*xloc(1)*xloc(1)*bl3jmp/b_out &
                     + xloc(1)*xloc(3)*xzz*rho1 &
                     + xloc(1)*xloc(2)*xyz*rho1 &
                     + xloc(3)
               ca(5,nc)  = 0.5*xloc(1)*xloc(1)*rho
               ca(6,nc)  = 0.5*xloc(2)*xloc(2) &
                     + 0.5*(rho-1)*xloc(1)*xloc(1)
               ca(7,nc)  = 0.5*xloc(3)*xloc(3) &
                     + 0.5*(rho-1)*xloc(1)*xloc(1)
               ca(8,nc)  = xloc(1)*xloc(2)*rho
               ca(9,nc)  = xloc(1)*xloc(3)*rho
               ca(10,nc) = xloc(2)*xloc(3)
            end if  ! (index(i,j,k) <= 3)
         end do  !  k=k0-1,k0+1
      end do  !  j=j0-1,j0+1
   end do  !  i=i0-1,i0+1


   if (info > 3) then
      cb(1) = fkjmp
   end if

   cb(2)  = -bl_in(1)
   cb(3)  = -bl_in(2)
   cb(4)  = -bl_in(3)
   cb(5)  = -b_in
   cb(6)  = -b_in
   cb(7)  = -b_in

   do i=1,inmax
      do j=1,inmax
         cc(i,j) = 0.0
      end do
      cc(i,i) = 1.0
   end do
   hxyz=hx*hx+hy*hy+hz*hz
   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1
            ndis=abs(i-i0)+abs(j-j0)+abs(k-k0)
            if (ndis <= 1) then
               if (phi(i,j,k) <= 0.0) then
                  cd(nc)=-3.0*b_in!fb_in(b_in,b_out,x(i),y(j),z(k))/hxyz
               else
                  cd(nc)=-3.0*b_out!fb_out(b_in,b_out,x(i),y(j),z(k))/hxyz
               end if
            else
               cd(nc)=0.0
            end if

            if (ndis == 0) cd(nc)=-6.0*cd(nc)

         end do
      end do
   end do
   if (info > 3) then
      call cpymat(immax,inmax,ca,caa)
   end if

   iprint = 0
   iwar(1) = 1

   ! ----- solve quadratic programming

   call ql0001(im,ime,immax,in,inmax,imnn,cc,cd,ca,cb,xl,xu, &
      xx,uu,iout,ifail,iprint,war,lwar,iwar,liwar,eps)

   tmp=0.0
   do ii=1,10
      tmp1=0.0
      do ji=1,inmax
         tmp1=tmp1+ca(ii,ji)*xx(ji)
      end do
      if(abs(tmp1+cb(ii)) > tmp) tmp = abs(tmp1+cb(ii))
   end do

   if(ifail > 10) then
      write(*,*) "Warning!", i0,j0,k0, "Irre31 Inconsistent Constraints"
      return
      do i=1, in
         xx(i) = 0.0
      end do
      call regula(l,m,n,hx,hy,hz,b_in,b_out,i0,j0,k0,info,x,y,z, tmpcoe,rhs)
      xx(5)  = tmpcoe(2)
      xx(23) = tmpcoe(3)
      xx(11) = tmpcoe(4)
      xx(17) = tmpcoe(5)
      xx(13) = tmpcoe(6)
      xx(15) = tmpcoe(7)
      xx(14) = tmpcoe(1)
      return
   end if

   do i=1,20
      a(i) = 0.0
   end do

   !About the local coordinates:note that in the test point i dir and ksi are opposite

   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1

            nc = nc + 1
            tmpx(1) = x(i) - x1
            tmpx(2) = y(j) - y1
            tmpx(3) = z(k) - z1
            call matvec(3,3,t, tmpx, xloc)

            if (index(i,j,k) <= 3) then

               a(1) =a(1) +xx(nc)
               a(3) =a(3) +xx(nc)*xloc(1)
               a(5) =a(5) +xx(nc)*xloc(2)
               a(7) =a(7) +xx(nc)*xloc(3)
               a(9) =a(9) +xx(nc)*xloc(1)*xloc(1)*0.5
               a(11)=a(11)+xx(nc)*xloc(2)*xloc(2)*0.5
               a(13)=a(13)+xx(nc)*xloc(3)*xloc(3)*0.5
               a(15)=a(15)+xx(nc)*xloc(1)*xloc(2)
               a(17)=a(17)+xx(nc)*xloc(1)*xloc(3)
               a(19)=a(19)+xx(nc)*xloc(2)*xloc(3)
            else
               a(2) =a(2) +xx(nc)
               a(4) =a(4) +xx(nc)*xloc(1)
               a(6) =a(6) +xx(nc)*xloc(2)
               a(8) =a(8) +xx(nc)*xloc(3)
               a(10)=a(10)+xx(nc)*xloc(1)*xloc(1)*0.5
               a(12)=a(12)+xx(nc)*xloc(2)*xloc(2)*0.5
               a(14)=a(14)+xx(nc)*xloc(3)*xloc(3)*0.5
               a(16)=a(16)+xx(nc)*xloc(1)*xloc(2)
               a(18)=a(18)+xx(nc)*xloc(1)*xloc(3)
               a(20)=a(20)+xx(nc)*xloc(2)*xloc(3)
            end if
         end do
      end do  !  j=j0-1,j0+1
   end do  !  i=i0-1,i0+1
   tmp  = 1.0/b_out
   tmp1 =   a(4) + a(10)*(xyy+xzz-bl_out(1)*tmp) &
         - a(12)*xyy - a(14)*xzz - a(16)*bl_out(2)*tmp &
         - a(18)*bl_out(3)*tmp - a(20)*xyz
   tmp2 =   a(6) - a(10)*bl_out(2)*tmp &
         + a(16)*xyy + a(18)*xyz
   tmp3 =   a(8) - a(10)*bl_out(3)*tmp &
         + a(16)*xyz + a(18)*xzz
   corr =   a(10)*((fjmp-fkk_out*ujmp)*tmp-wyy-wzz) &
         + a(12)*wyy + a(14)*wzz + a(16)*qy*tmp &
         + a(18)*qz*tmp + a(20)*wyz + a(2)*ujmp &
         + tmp*tmp1*q + tmp2*wy + tmp3*wz

   if (info > 3) then
      corr = corr + fkk_out*ujmp - fjmp
   end if

   if (info <= 3) then
      fkk = 0.0d0!fk_in(x(i0),y(j0),z(k0))
      rhs = 0.0d0!ff_in(x(i0),y(j0),z(k0))
   else
      fkk = 0.0d0!fk_out(x(i0),y(j0),z(k0))
      rhs = 0.0d0!ff_out(x(i0),y(j0),z(k0))
   end if

   xx(kctr) = xx(kctr)+fkk
   rhs = rhs + corr

end subroutine irre31
