!!
!! Copyright (C) 2011 Andreas van Hameren. 
!!
!! This file is part of OneLOop-2.2.1.
!!
!! OneLOop-2.2.1 is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! OneLOop-2.2.1 is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with OneLOop-2.2.1.  If not, see <http://www.gnu.org/licenses/>.
!!


module avh_olo_func
  use avh_olo_kinds
  use avh_olo_units
!
  implicit none
!
  type :: qmplx_type
   complex(kindc2) :: c
   integer         :: p
  end type
!
  interface mysqrt
    module procedure mysqrt_0,mysqrt_r,mysqrt_i
  end interface
!
  interface qonv
    module procedure qonv_r,qonv_0,qonv_i
  end interface
!
  interface operator (*)
    module procedure prduct,prduct_r
  end interface
  interface operator (/)
    module procedure ratio,ratio_r
  end interface
!
  interface eta5
    module procedure eta5_0
  end interface
  interface eta3
    module procedure eta3_r,eta3_0
  end interface
  interface eta2
    module procedure eta2_r,eta2_0
  end interface
!
contains
!
!
   function mysqrt_0(xx) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! negative imaginary.
!*******************************************************************
   complex(kindc2) ,intent(in) :: xx
   complex(kindc2) :: rslt ,zz
   real(kindr2) :: xim,xre
   xim = aimag(xx)
   if (xim.eq.R0P0) then
     xre = real(xx)
     if (xre.ge.R0P0) then
       zz = cmplx(sqrt(xre),R0P0,kind=kindc2)
     else
       zz = cmplx(R0P0,-sqrt(-xre),kind=kindc2)
     endif
   else
     zz = sqrt(xx)
   endif
   rslt = zz
   end function

   function mysqrt_r(xx,sgn) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! imaginary and has the same sign as  sgn .
!*******************************************************************
   complex(kindc2) ,intent(in) :: xx
   real(kindr2)    ,intent(in) :: sgn
   complex(kindc2) :: rslt ,zz
   real(kindr2) :: xim,xre
   xim = aimag(xx)
   if (xim.eq.R0P0) then
     xre = real(xx)
     if (xre.ge.R0P0) then
       zz = cmplx(sqrt(xre),R0P0,kind=kindc2)
     else
       zz = cmplx(R0P0,sign(sqrt(-xre),sgn),kind=kindc2)
     endif
   else
     zz = sqrt(xx)
   endif
   rslt = zz
   end function

   function mysqrt_i(xx,sgn) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! imaginary and has the same sign as  sgn .
!*******************************************************************
   complex(kindc2) ,intent(in) :: xx
   integer         ,intent(in) :: sgn
   complex(kindc2) :: rslt ,zz
   real(kindr2) :: xim,xre
   xim = aimag(xx)
   if (xim.eq.R0P0) then
     xre = real(xx)
     if (xre.ge.R0P0) then
       zz = cmplx(sqrt(xre),R0P0,KIND=kindc2)
     else
       zz = cmplx(R0P0,sign(sqrt(-xre),real(sgn,KIND=kindr2)),KIND=kindc2)
     endif
   else
     zz = sqrt(xx)
   endif
   rslt = zz
   end function


   subroutine solabc( x1,x2 ,dd ,aa,bb,cc ,imode )
!*******************************************************************
! Returns the solutions  x1,x2  to the equation  aa*x^2+bb*x+cc=0
! Also returns  dd = aa*(x1-x2)
! If  imode=/=0  it uses  dd  as input as value of  sqrt(b^2-4*a*c)
!*******************************************************************
   complex(kindc2) ,intent(out)   :: x1,x2
   complex(kindc2) ,intent(inout) :: dd
   complex(kindc2) ,intent(in) :: aa,bb,cc
   integer         ,intent(in) :: imode
   complex(kindc2) :: qq,hh
   real(kindr2) :: r1,r2
!
   if (aa.eq.C0P0) then
     if (bb.eq.C0P0) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop solabc: ' &
         ,'no solutions, returning 0'
       x1 = C0P0
       x2 = C0P0
       dd = C0P0
     else
       x1 = -cc/bb
       x2 = x1
       dd = bb
     endif
   elseif (cc.eq.C0P0) then
     dd = -bb
     x1 = dd/aa
     x2 = C0P0
   else
     if (imode.eq.0) dd = sqrt(bb*bb - 4*aa*cc)
     qq = -bb+dd
     hh = -bb-dd
     r1 = abs(qq)
     r2 = abs(hh)
     if (r1.ge.r2) then
       x1 = qq/(2*aa)
       x2 = (2*cc)/qq
     else
       qq = hh
       x2 = qq/(2*aa)
       x1 = (2*cc)/qq
     endif
   endif
   end subroutine


   subroutine rfun(rr,dd ,qq)
!*******************************************************************
! Returns  rr  such that  qq = rr + 1/rr  and  Im(rr)  has the same
! sign as  Im(qq) .
! If  Im(qq)  is zero, then  Im(rr)  is negative or zero.
! If  Im(rr)  is zero, then  |rr| > 1/|rr| .
! Also returns  dd = rr - 1/rr .
!*******************************************************************
   complex(kindc2) ,intent(out) :: rr,dd
   complex(kindc2) ,intent(in)  :: qq
   complex(kindc2) :: r2
   real(kindr2) :: aa,bb
   integer :: ir,ik
   complex(kindc2) ,parameter :: two=2*C1P0,four=4*C1P0
   dd = sqrt(qq*qq-four)
   rr = qq+dd
   r2 = qq-dd
   aa = abs(rr)
   bb = abs(r2)
   if (bb.gt.aa) then
     rr = r2
     dd = -dd
   endif
   aa = aimag(qq)
   bb = aimag(rr)
   if (aa.eq.R0P0) then
     if (bb.le.R0P0) then
       rr = rr/two
     else
       rr = two/rr
       dd = -dd
     endif
   else
     ik = int(sign(R1P0,aa))
     ir = int(sign(R1P0,bb))
     if (ir.eq.ik) then
       rr = rr/two
     else
       rr = two/rr
       dd = -dd
     endif
   endif
   end subroutine

   subroutine rfun0(rr ,dd,qq)
!*******************************************************************
! Like rfun, but now  dd  is input, which may get a minus sign
!*******************************************************************
   complex(kindc2) ,intent(out)   :: rr
   complex(kindc2) ,intent(inout) :: dd
   complex(kindc2) ,intent(in)  :: qq
   complex(kindc2) :: r2
   real(kindr2) :: aa,bb
   integer :: ir,ik
   complex(kindc2) ,parameter :: two=2*C1P0
   rr = qq+dd
   r2 = qq-dd
   aa = abs(rr)
   bb = abs(r2)
   if (bb.gt.aa) then
     rr = r2
     dd = -dd
   endif
   aa = aimag(qq)
   bb = aimag(rr)
   if (aa.eq.R0P0) then
     if (bb.le.R0P0) then
       rr = rr/two
     else
       rr = two/rr
       dd = -dd
     endif
   else
     ik = int(sign(R1P0,aa))
     ir = int(sign(R1P0,bb))
     if (ir.eq.ik) then
       rr = rr/two
     else
       rr = two/rr
       dd = -dd
     endif
   endif
   end subroutine


   function qonv_r(xx,sgn) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz  becomes the
! sign of  sgn .
!*******************************************************************
   complex(kindc2) ,intent(in) :: xx
   real(kindr2)    ,intent(in) :: sgn
   type(qmplx_type) :: rslt
   real(kindr2) :: xre,xim
   xre = real(xx)
   if (xre.ge.R0P0) then
     rslt%c = xx
     rslt%p = 0
   else
     xim = aimag(xx)
     if (xim.eq.R0P0) then
       rslt%c = cmplx(-xre,R0P0,kind=kindc2)
       rslt%p = int(sign(R1P0,sgn))
     else
       rslt%c = -xx
       rslt%p = int(sign(R1P0,xim)) ! xim = -Im(rslt%c)
     endif
   endif
   end function
!
   function qonv_i(xx,sgn) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz  becomes the
! sign of  sgn .
!*******************************************************************
   complex(kindc2) ,intent(in) :: xx
   integer         ,intent(in) :: sgn
   type(qmplx_type) :: rslt
   real(kindr2) :: xre,xim
   xre = real(xx)
   if (xre.ge.R0P0) then
     rslt%c = xx
     rslt%p = 0
   else
     xim = aimag(xx)
     if (xim.eq.R0P0) then
       rslt%c = cmplx(-xre,R0P0,kind=kindc2)
       rslt%p = sign(1,sgn)
     else
       rslt%c = -xx
       rslt%p = int(sign(R1P0,xim)) ! xim = -Im(rslt%c)
     endif
   endif
   end function
!
   function qonv_0(xx) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz=1
!*******************************************************************
   complex(kindc2) ,intent(in) :: xx
   type(qmplx_type) :: rslt
   real(kindr2) :: xre,xim
   xre = real(xx)
   if (xre.ge.R0P0) then
     rslt%c = xx
     rslt%p = 0
   else
     xim = aimag(xx)
     if (xim.eq.R0P0) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop qonv: ' &
         ,'negative input with undefined sign for the imaginary part, ' &
         ,'putting +ieps'
       rslt%c = cmplx(-xre,R0P0,kind=kindc2)
       rslt%p = 1
     else
       rslt%c = -xx
       rslt%p = int(sign(R1P0,xim)) ! xim = -Im(rslt%c)
     endif
   endif
   end function
!
   function directly(xx,ix) result(rslt)
!*******************************************************************
!*******************************************************************
   complex(kindc2) ,intent(in) :: xx
   integer         ,intent(in) :: ix
   type(qmplx_type) :: rslt
   rslt%c = xx
   rslt%p = ix
   end function


   function sheet(xx) result(ii)
!*******************************************************************
! Returns the number of the Riemann-sheet (times 2) for the complex
! number  xx*exp(ix*imag*pi) . The real part of xx is assumed to be
! positive or zero. Examples:
! xx=1+imag, ix=-1 -> ii= 0 
! xx=1+imag, ix= 1 -> ii= 2 
! xx=1-imag, ix=-1 -> ii=-2 
! xx=1-imag, ix= 1 -> ii= 0 
! xx=1     , ix= 1 -> ii= 0  convention that log(-1)=pi on
! xx=1     , ix=-1 -> ii=-2  the principal Riemann-sheet
!*******************************************************************
   type(qmplx_type) ,intent(in) :: xx
   integer :: ii,jj
   real(kindr2) :: xim
   ii = xx%p/2*2
   jj = xx%p-ii
   xim = aimag(xx%c)
   if (xim.le.R0P0) then ! also xim=0 <==> log(-1)=pi, not -pi
     if (jj.eq.-1) ii = ii-2
   else
     if (jj.eq. 1) ii = ii+2
   endif
   end function


   function prduct(yy,xx) result(zz)
!*******************************************************************
! Return the product  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
   type(qmplx_type) ,intent(in) :: yy,xx
   type(qmplx_type) :: zz
   zz%c = yy%c*xx%c
   zz%p = yy%p+xx%p
   if (real(zz%c).lt.R0P0) then
     zz%p = zz%p + int(sign(R1P0,aimag(xx%c)))
     zz%c = -zz%c
   endif
   end function

   function prduct_r(yy,xx) result(zz)
!*******************************************************************
! Return the product  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
   type(qmplx_type) ,intent(in) :: yy
   real(kindr2)     ,intent(in) :: xx
   type(qmplx_type) :: zz
   zz%c = yy%c*abs(xx)
   zz%p = yy%p
   end function

   function ratio(yy,xx) result(zz)
!*******************************************************************
! Return the ratio  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
   type(qmplx_type) ,intent(in) :: yy,xx
   type(qmplx_type) :: zz
   zz%c = yy%c/xx%c
   zz%p = yy%p-xx%p
   if (real(zz%c).lt.R0P0) then
     zz%p = zz%p - int(sign(R1P0,aimag(xx%c)))
     zz%c = -zz%c
   endif
   end function
!
   function ratio_r(yy,xx) result(zz)
!*******************************************************************
!*******************************************************************
   type(qmplx_type) ,intent(in) :: yy
   real(kindr2)     ,intent(in) :: xx
   type(qmplx_type) :: zz
   zz%c = yy%c/abs(xx)
   zz%p = yy%p
   end function
!
!
   function eta3_r( aa,sa ,bb,sb ,cc,sc ) result(rslt)
!*******************************************************************
! 2*pi*imag times the result of
!     theta(-Im(a))*theta(-Im(b))*theta( Im(c))
!   - theta( Im(a))*theta( Im(b))*theta(-Im(c))
! where a,b,c are interpreted as a+i|eps|sa, b+i|eps|sb, c+i|eps|sc
!*******************************************************************
   complex(kindc2) ,intent(in) :: aa,bb,cc
   real(kindr2)    ,intent(in) :: sa,sb,sc
   complex(kindc2) :: rslt
   real(kindr2) :: ima,imb,imc
   ima = aimag(aa)
   imb = aimag(bb)
   imc = aimag(cc)
   if (ima.eq.R0P0) ima = sa
   if (imb.eq.R0P0) imb = sb
   if (imc.eq.R0P0) imc = sc
   ima = sign(R1P0,ima)
   imb = sign(R1P0,imb)
   imc = sign(R1P0,imc)
   if (ima.eq.imb.and.ima.ne.imc) then
     rslt = cmplx(R0P0,imc*TWOPI,kind=kindc2)
   else
     rslt = C0P0
   endif
   end function
!
   function eta3_0( aa ,bb ,cc ) result(rslt)
!*******************************************************************
! 2*pi*imag times the result of
!     theta(-Im(a))*theta(-Im(b))*theta( Im(c))
!   - theta( Im(a))*theta( Im(b))*theta(-Im(c))
! where a,b,c are interpreted as a+i|eps|sa, b+i|eps|sb, c+i|eps|sc
!*******************************************************************
   complex(kindc2) ,intent(in) :: aa,bb,cc
   complex(kindc2) :: rslt
   real(kindr2) :: ima,imb,imc
   ima = aimag(aa)
   imb = aimag(bb)
   imc = aimag(cc)
   ima = sign(R1P0,ima)
   imb = sign(R1P0,imb)
   imc = sign(R1P0,imc)
   if (ima.eq.imb.and.ima.ne.imc) then
     rslt = cmplx(R0P0,imc*TWOPI,kind=kindc2)
   else
     rslt = C0P0
   endif
   end function
!
   function eta5_0( aa ,b1,c1 ,b2,c2 ) result(rslt)
!*******************************************************************
! eta3(aa,b1,c1) - eta3(aa,b2,c2)
!*******************************************************************
   complex(kindc2) ,intent(in) :: aa,b1,c1 ,b2,c2
   complex(kindc2) :: rslt
   real(kindr2) :: imaa,imb1,imc1,imb2,imc2
   imaa = sign(R1P0,aimag(aa))
   imb1 = sign(R1P0,aimag(b1))
   imb2 = sign(R1P0,aimag(b2))
   imc1 = sign(R1P0,aimag(c1))
   imc2 = sign(R1P0,aimag(c2))
   if (imaa.eq.imb1) then
     if (imaa.eq.imb2) then
       if (imc1.eq.imc2) then
         rslt = C0P0
       elseif (imaa.ne.imc1) then
         rslt = cmplx(R0P0, imc1*TWOPI,kind=kindc2)
       else
         rslt = cmplx(R0P0,-imc2*TWOPI,kind=kindc2)
       endif
     elseif (imaa.ne.imc1) then
       rslt = cmplx(R0P0, imc1*TWOPI,kind=kindc2)
     else
       rslt = C0P0
     endif
   elseif (imaa.eq.imb2.and.imaa.ne.imc2) then
     rslt = cmplx(R0P0,-imc2*TWOPI,kind=kindc2)
   else
     rslt = C0P0
   endif
   end function
 
   function eta2_r( aa,sa ,bb,sb ) result(rslt)
!*******************************************************************
! The same as  eta3, but with  c=a*b, so that
!   eta(a,b) = log(a*b) - log(a) - log(b)
!*******************************************************************
   complex(kindc2) ,intent(in) :: aa,bb
   real(kindr2)    ,intent(in) :: sa,sb
   complex(kindc2) :: rslt
   real(kindr2) :: rea,reb,ima,imb,imab
   rea = real(aa)  ;ima = aimag(aa)
   reb = real(bb)  ;imb = aimag(bb)
   imab = rea*imb + reb*ima
   if (ima.eq.R0P0) ima = sa
   if (imb.eq.R0P0) imb = sb
   if (imab.eq.R0P0) imab = sign(rea,sb) + sign(reb,sa)
   ima  = sign(R1P0,ima)
   imb  = sign(R1P0,imb)
   imab = sign(R1P0,imab)
   if (ima.eq.imb.and.ima.ne.imab) then
     rslt = cmplx(R0P0,imab*TWOPI,kind=kindc2)
   else
     rslt = C0P0
   endif
   end function
! 
   function eta2_0( aa ,bb ) result(rslt)
!*******************************************************************
!*******************************************************************
   complex(kindc2) ,intent(in) :: aa,bb
   complex(kindc2) :: rslt
   real(kindr2) :: rea,reb,ima,imb,imab
   rea = real(aa)  ;ima = aimag(aa)
   reb = real(bb)  ;imb = aimag(bb)
   rea = rea*imb
   reb = reb*ima
   imab = rea+reb
   ima  = sign(R1P0,ima)
   imb  = sign(R1P0,imb)
   imab = sign(R1P0,imab)
   if (ima.eq.imb.and.ima.ne.imab) then
     rslt = cmplx(R0P0,imab*TWOPI,kind=kindc2)
   else
     rslt = C0P0
   endif
   end function 
!
end module


module avh_olo_loga
!*******************************************************************
! log( |xx|*exp(imag*pi*iph) ) = log|xx| + imag*pi*iph
!*******************************************************************
  use avh_olo_kinds
  use avh_olo_units
  use avh_olo_func
  implicit none
  private
  public :: loga
  real(kindr2) ,parameter :: pi=TWOPI/2
contains
!
  function loga(xx,iph) result(rslt)
  real(kindr2) ,intent(in) :: xx
  integer      ,intent(in) :: iph
  complex(kindc2) :: rslt
  real(kindr2) :: rr
!
  rr = abs(xx)
  if (rr.eq.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop loga: ' &
      ,'|xx|=',rr
  endif
  rslt = cmplx(log(rr),iph*pi,kind=kindc2)
  end function
!
end module


module avh_olo_bern
!*******************************************************************
! the first nn Bernoulli numbers
!*******************************************************************
  use avh_olo_kinds
  implicit none
  private
  public :: init_bern,rbern,cbern
  integer ,parameter :: nn=40
  real(kindr2)    :: rbern(nn) !PROTECTED
  complex(kindc2) :: cbern(nn) !PROTECTED
  integer :: ndigits=0

  ! protected :: rbern, cbern
contains
!  
  subroutine init_bern(ndig)
  integer ,intent(in) :: ndig
  integer :: jj
  integer ,parameter :: d=kindr2
  if (ndigits.eq.ndig) return ;ndigits=ndig
  rbern(1:nn) = R0P0
  rbern( 1) = -1._d/2._d
  rbern( 2) =  1._d/6._d
  rbern( 4) = -1._d/30._d
  rbern( 6) =  1._d/42._d
  rbern( 8) = -1._d/30._d
  rbern(10) =  5._d/66._d
  rbern(12) = -691._d/2730._d
  rbern(14) =  7._d/6._d
  rbern(16) = -3617._d/510._d
  rbern(18) =  43867._d/798._d
  rbern(20) = -174611._d/330._d
  rbern(22) =  854513._d/138._d
  rbern(24) = -236364091._d/2730._d
  rbern(26) =  8553103._d/6._d
  rbern(28) = -23749461029._d/870._d
  rbern(30) =  8615841276005._d/14322._d
  rbern(32) = -7709321041217._d/510._d
  rbern(34) =  2577687858367._d/6._d
  rbern(36) = -26315271553053477373._d/1919190._d
  rbern(38) =  2929993913841559._d/6._d
  rbern(40) = -261082718496449122051._d/13530._d
  do jj=1,nn
    cbern(jj) = cmplx(rbern(jj),kind=kindc2)
  enddo
  end subroutine
!
end module


module avh_olo_li2a
!*******************************************************************
!                  /1    ln(1-zz*t)
! avh_olo_li2a = - |  dt ---------- 
!                  /0        t
! with  zz = 1 - |xx|*exp(imag*pi*iph)
! Examples:
!   In order to get the dilog of  1.1  use  xx=1.1, iph=0
!   In order to get the dilog of -1.1  use  xx=1.1, iph=1
! Add multiples of  2  to  iph  in order to get the result on
! different Riemann-sheets.
!*******************************************************************
  use avh_olo_kinds
  use avh_olo_func
  use avh_olo_bern
  implicit none
  private
  public :: init_li2a,li2a
  real(kindr2)    ,parameter :: pi=TWOPI/2
  complex(kindc2) ,parameter :: pi2o6=C1P0*TWOPI*TWOPI/24
  integer :: nn=16
  integer :: ndigits=0
contains
!
  subroutine init_li2a(ndig)
  integer ,intent(in) :: ndig
  if (ndigits.eq.ndig) return ;ndigits=ndig
  call init_bern(ndigits)
  if     (ndigits.lt.24) then
    nn = 16
  else
    nn = 30
  endif
  end subroutine

  function li2a(xx,iph) result(rslt)
  real(kindr2) ,intent(in) :: xx
  integer      ,intent(in) :: iph
  complex(kindc2) :: rslt
  real(kindr2) :: rr,yy,lyy,loy,zz,z2,liox
  integer :: ii,ntwo,ione
  logical :: positive , r_gt_1 , y_lt_h
!
  rr = abs(xx)
  ntwo = iph/2*2
  ione = iph - ntwo
  positive = (ione.eq.0)
! 
  if     (rr.eq.R0P0) then
    rslt = pi2o6
  elseif (rr.eq.R1P0.and.positive) then
    rslt = C0P0
  else
    yy  = rr
    lyy = log(rr)
    if (.not.positive) yy = -yy
!
    r_gt_1 = (rr.gt.R1P0)
    if (r_gt_1) then
      yy   = R1P0/yy
      lyy  = -lyy
      ntwo = -ntwo
      ione = -ione
    endif
    loy = log(R1P0-yy) ! log(1-yy) is always real
!
    y_lt_h = (yy.lt.R5M1)
    if (y_lt_h) then
      zz = -loy ! log(1-yy) is real
    else
      zz = -lyy ! yy>0.5 => log(yy) is real
    endif
!
    z2 = zz*zz
    liox = rbern(nn)
    do ii=nn,4,-2
      liox = rbern(ii-2) + liox*z2/(ii*(ii+1))
    enddo
    liox = rbern(1) + liox*zz/3
    liox = zz + liox*z2/2
!
    rslt = cmplx(liox,kind=kindc2)
!
    if (y_lt_h) then
      rslt = pi2o6 - rslt - cmplx(loy*lyy,loy*pi*ione,kind=kindc2)
    endif
!
    rslt = rslt + cmplx( R0P0 , -loy*pi*ntwo ,kind=kindc2)
!
    if (r_gt_1) rslt = -rslt - cmplx(-lyy,iph*pi,kind=kindc2)**2/2
  endif
  end function
!
end module


module avh_olo_loga2
!*******************************************************************
! log(xx)/(1-xx)  with  xx = log|xx| + imag*pi*iph
!*******************************************************************
  use avh_olo_kinds
  use avh_olo_units
  use avh_olo_func
  implicit none
  private
  public :: init_loga2,loga2
  real(kindr2) :: thrs=epsilon(R1P0)
  integer :: ndigits=0
contains
!
  subroutine init_loga2(ndig)
  integer ,intent(in) :: ndig
  if (ndigits.eq.ndig) return ;ndigits=ndig
  thrs = 10*thrs
  end subroutine
!
  function loga2(xx,iph) result(rslt)
  use avh_olo_loga ,only : loga
  real(kindr2) ,intent(in) :: xx
  integer      ,intent(in) :: iph
  complex(kindc2) :: rslt
  real(kindr2) :: omx
!
  if (mod(iph,2).eq.0) then
    omx = R1P0-abs(xx)
  else
    omx = R1P0+abs(xx)
  endif
!
  if (iph.ne.0) then
    if (omx.eq.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop loga2: ' &
        ,'1-xx,iph=',omx,iph
      rslt = C0P0
    else
      rslt = loga(xx,iph)/cmplx(omx,kind=kindc2)
    endif
  else
    if (abs(omx).lt.thrs) then
      rslt = cmplx(-R1P0-omx/2,kind=kindc2)
    else
      rslt = loga(xx,iph)/cmplx(omx,kind=kindc2)
    endif
  endif
  end function
!
end module


module avh_olo_logc
!*******************************************************************
! Returns  log( |Re(xx)| + imag*Im(xx) ) + imag*pi*iph
!*******************************************************************
  use avh_olo_kinds
  use avh_olo_units
  use avh_olo_func
  implicit none
  private
  public :: logc
  complex(kindc2) ,parameter :: ipi=CiP0*TWOPI/2
contains
!
  function logc(xx) result(rslt)
  type(qmplx_type) ,intent(in) :: xx
  complex(kindc2) :: rslt
  if (xx%c.eq.C0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop logc: xx%c =',xx%c
    rslt = C0P0
  else
    rslt = log( cmplx(abs(real(xx%c)),aimag(xx%c),kind=kindc2) ) &
         + ipi*xx%p
  endif
  end function
!
end module


module avh_olo_li2c
!*******************************************************************
!                  /1    ln(1-zz*t)
! avh_olo_li2c = - |  dt ---------- 
!                  /0        t
! with  zz = 1 - ( |Re(xx)| + imag*Im(xx) )*exp(imag*pi*iph)
! Examples:
!   In order to get the dilog of  1+imag  use  xx=1+imag, iph= 0
!   In order to get the dilog of  1-imag  use  xx=1-imag, iph= 0
!   In order to get the dilog of -1+imag  use  xx=1-imag, iph= 1
!   In order to get the dilog of -1-imag  use  xx=1+imag, iph=-1
! Add multiples of  2  to  iph  in order to get the result on
! different Riemann-sheets.
!*******************************************************************
  use avh_olo_kinds
  use avh_olo_func
  use avh_olo_bern
  use avh_olo_li2a
  implicit none
  private
  public :: init_li2c,li2c
  complex(kindc2) ,parameter :: ipi=CiP0*TWOPI/2
  complex(kindc2) ,parameter :: pi2o6=C1P0*TWOPI*TWOPI/24
  integer :: nn=18
  integer :: ndigits=0
contains
!
  subroutine init_li2c(ndig)
  integer ,intent(in) :: ndig
  if (ndigits.eq.ndig) return ;ndigits=ndig
  call init_li2a(ndigits)
  call init_bern(ndigits)
  if     (ndigits.lt.24) then
    nn = 18
  else
    nn = 36
  endif
  end subroutine

  function li2c(xx) result(rslt)
  type(qmplx_type) :: xx
  complex(kindc2) :: rslt ,yy,lyy,loy,zz,z2
  real(kindr2) :: rex,imx
  integer :: ii,iyy
  logical :: x_gt_1 , y_lt_h
!
  rex = real(xx%c)
  imx = aimag(xx%c)
! 
  if (imx.eq.R0P0) then
    rslt = li2a(rex,xx%p)
  else
    rex = abs(rex)
!
    if (mod(xx%p,2).eq.0) then
      yy = cmplx(rex,imx,kind=kindc2)
      iyy = xx%p
    else
      yy = cmplx(-rex,-imx,kind=kindc2)
! Notice that  iyy=xx%p/2*2  does not deal correctly with the
! situation when  xx%p-xx%p/2*2 = sign(Im(xx%c)) . The following does:
      iyy = xx%p + nint(sign(R1P0,imx))
    endif
!
    x_gt_1 = (abs(xx%c).gt.R1P0)
    if (x_gt_1) then
      yy = C1P0/yy
      iyy = -iyy
    endif
    lyy = log(yy)
    loy = log(C1P0-yy)
!
    y_lt_h = (real(yy).lt.R5M1)
    if (y_lt_h) then
      zz = -loy
    else
      zz = -lyy
    endif
!
    z2 = zz*zz
    rslt = cbern(nn)
    do ii=nn,4,-2
      rslt = cbern(ii-2) + rslt*z2/(ii*(ii+1))
    enddo
    rslt = cbern(1) + rslt*zz/3
    rslt = zz + rslt*z2/2
!
    if (y_lt_h) rslt = pi2o6 - rslt - loy*lyy
!
    rslt = rslt - loy*ipi*iyy
!
    if (x_gt_1) rslt = -rslt - (lyy + ipi*iyy)**2/2
  endif
  end function
!
end module


module avh_olo_logc2
!*******************************************************************
! log(xx)/(1-xx)
! with  log(xx) = log( |Re(xx)| + imag*Im(xx) ) + imag*pi*iph
!*******************************************************************
  use avh_olo_kinds
  use avh_olo_units
  use avh_olo_func
  implicit none
  private
  public :: init_logc2,logc2
  real(kindr2) :: thrs=epsilon(R1P0)
  integer :: ndigits=0
contains
!
  subroutine init_logc2(ndig)
  integer ,intent(in) :: ndig
  if (ndigits.eq.ndig) return ;ndigits=ndig
  thrs = 10*thrs
  end subroutine
!
  function logc2(xx) result(rslt)
  use avh_olo_logc ,only : logc
  type(qmplx_type) ,intent(in) :: xx
  complex(kindc2) :: rslt ,omx
  if (mod(xx%p,2).eq.0) then
    omx = cmplx(1d0-abs(real(xx%c)),-aimag(xx%c),kind=kindc2)
  else
    omx = cmplx(1d0+abs(real(xx%c)), aimag(xx%c),kind=kindc2)
  endif
  if (xx%p.ne.0) then
    if (omx.eq.C0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop logc2: ' &
        ,'1-xx%c,xx%p=',omx,xx%p
      rslt = C0P0
    else
      rslt = logc(xx)/omx
    endif
  else
    if (abs(omx).lt.thrs) then
      rslt = -C1P0-omx/2
    else
      rslt = logc(xx)/omx
    endif
  endif
  end function
!
end module


module avh_olo_li2c2
!*******************************************************************
! avh_olo_li2c2 = ( li2(x1) - li2(x2) )/(x1%c-x2%c)
!
!                    /1    ln(1-zz*t)
! where  li2(x1) = - |  dt ----------
!                    /0        t
! with  zz = 1 - ( |Re(x1%c)| + imag*Im(x1%c) )*exp(imag*pi*x1%p)
! and similarly for li2(x2)
!*******************************************************************
  use avh_olo_kinds
  use avh_olo_units
  use avh_olo_func
  use avh_olo_li2c
  use avh_olo_logc2
  implicit none
  private
  public :: init_li2c2,li2c2
  complex(kindc2) ,parameter :: ipi=CiP0*TWOPI/2
  real(kindr2)    ,parameter :: thrs1=epsilon(R1P0)
  real(kindr2) :: thrs=0.11_kindr2
  integer      :: nmax=12
  integer :: ndigits=0
contains
!
  subroutine init_li2c2(ndig)
  integer ,intent(in) :: ndig
  if (ndigits.eq.ndig) return ;ndigits=ndig
  call init_logc2(ndigits)
  call init_li2c(ndigits)
  if     (ndigits.lt.16) then
    thrs = 0.11_kindr2 ! double precision
    nmax = 12
  elseif (ndigits.lt.24) then
    thrs = 0.02_kindr2 ! guess
    nmax = 12
  else
    thrs = 0.008_kindr2 ! quadruple precision
    nmax = 12
  endif
  end subroutine
!
  function li2c2(x1,x2) result(rslt)
  type(qmplx_type) ,intent(in) :: x1,x2
  complex(kindc2) :: rslt
  complex(kindc2) :: x1r,x2r,delta,xx,xr,omx,del,hh,ff(0:20),zz
  integer :: ih,ii
!
  if (mod(x1%p,2).eq.0) then
    x1r = cmplx( abs(real(x1%c)), aimag(x1%c),kind=kindc2)
  else
    x1r = cmplx(-abs(real(x1%c)),-aimag(x1%c),kind=kindc2)
  endif     
  if (mod(x2%p,2).eq.0) then
    x2r = cmplx( abs(real(x2%c)), aimag(x2%c),kind=kindc2)
  else
    x2r = cmplx(-abs(real(x2%c)),-aimag(x2%c),kind=kindc2)
  endif
  delta = x1r-x2r
!
  if (x1%p.ne.x2%p) then
    if (delta.eq.C0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop li2c2: ' &
        ,'x1%p,x2%p,delta=',x1%p,x2%p,delta
      rslt = C0P0
    else
      rslt = ( li2c(x1)-li2c(x2) )/delta
    endif
  else
    if (abs(delta/x1%c).gt.thrs) then
      rslt = ( li2c(x1)-li2c(x2) )/delta
    else
      xx  = x1%c
      xr  = x1r
      omx = C1P0-xr
      del = delta
      hh = C1P0-x2r
      if (abs(hh).gt.abs(omx)) then
        xx = x2%c
        xr = x2r
        omx = hh
        del = -delta
      endif
      if (abs(omx).lt.thrs1) then
        zz = -C1P0-omx/2-del/4
      else
        ih = x1%p - x1%p/2*2
        ff(0) = logc2(directly(xx,ih))
        hh = -C1P0
        do ii=1,nmax
          hh = -hh/xr
          ff(ii) = ( hh/ii + ff(ii-1) )/omx
        enddo
        zz = ff(nmax)/(nmax+1)
        do ii=nmax-1,0,-1
          zz = ff(ii)/(ii+1) - zz*del
        enddo
      endif
      ih = x1%p-ih
      if (ih.ne.0) then
        omx = C1P0-x1r
        zz = zz - ih*ipi*logc2(qonv((C1P0-x2r)/omx))/omx
      endif
      rslt = zz
    endif
  endif
  end function
!
end module
