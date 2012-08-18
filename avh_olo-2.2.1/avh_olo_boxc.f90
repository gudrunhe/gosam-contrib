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


module avh_olo_boxc
   use avh_olo_kinds
   use avh_olo_units
   use avh_olo_func
   implicit none
   private
   public :: boxc,init_boxc
   real(kindr2) :: thrss3fun=epsilon(R1P0)*1000
   integer :: ndigits=0
contains

   subroutine init_boxc(ndig)
   integer ,intent(in) :: ndig
   if (ndigits.eq.ndig) return ;ndigits=ndig
   if     (ndigits.lt.16) then ;thrss3fun = epsilon(R1P0)*1000
   elseif (ndigits.lt.24) then ;thrss3fun = epsilon(R1P0)*30000
   else                        ;thrss3fun = epsilon(R1P0)*1000000
   endif
   end subroutine
   

   subroutine boxc( rslt ,pp_in,mm_in ,ap_in )
!*******************************************************************
! Finite 1-loop scalar 4-point function for complex internal masses
! Based on the formulas from
!   Dao Thi Nhung and Le Duc Ninh, arXiv:0902.0325 [hep-ph]
!   G. 't Hooft and M.J.G. Veltman, Nucl.Phys.B153:365-401,1979 
!*******************************************************************
   use avh_olo_box ,only: base,casetable,ll=>permtable
   complex(kindc2) ,intent(out)   :: rslt(0:2)
   complex(kindc2) ,intent(inout) :: pp_in(6),mm_in(4)
   real(kindr2)    ,intent(in)  :: ap_in(6)
   complex(kindc2) :: pp(6),mm(4)
   real(kindr2) :: ap(6),aptmp(6),rem,imm,hh
   complex(kindc2) :: a,b,c,d,e,f,g,h,j,k,dpe,epk,x1,x2,sdnt,o1,j1,e1
   integer :: icase,jcase,ii,jj
   integer ,parameter :: lp(6,3)=&
            reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
   integer ,parameter :: lm(4,3)=&
            reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
   real(kindr2) ,parameter :: small=epsilon(R1P0)
!
!  write(*,*) 'MESSAGE from OneLOop boxc: you are calling me' !CALLINGME
!
   rslt = C0P0
!
   hh = R0P0
   do ii=1,6
     aptmp(ii) = ap_in(ii)
     if (aptmp(ii).gt.hh) hh = aptmp(ii)
   enddo
   hh = 100*small*hh
   do ii=1,6
     if (aptmp(ii).lt.hh) aptmp(ii) = R0P0
   enddo
!
   if (aptmp(5).eq.R0P0.or.aptmp(6).eq.R0P0) then
     if (aptmp(1).eq.R0P0.or.aptmp(3).eq.R0P0) then
       if (aptmp(2).eq.R0P0.or.aptmp(4).eq.R0P0) then
         if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
           ,'no choice with |s| and |t| large enough, putting them by hand'
         if (aptmp(5).eq.R0P0) then
           aptmp(5) = hh
           pp_in(5) = cmplx(sign(hh,real(pp_in(5))),kind=kindc2)
         endif
         if (aptmp(6).eq.R0P0) then
           aptmp(6) = hh
           pp_in(6) = cmplx(sign(hh,real(pp_in(6))),kind=kindc2)
         endif
         jj = 1
       else
         jj = 3
       endif
     else
       jj = 2
     endif
   else
     jj = 1
   endif
   do ii=1,6
     ap(ii) = aptmp(lp(ii,jj))
     if (ap(ii).gt.R0P0) then ;pp(ii) = pp_in(lp(ii,jj))
     else                     ;pp(ii) = C0P0
     endif
   enddo
   do ii=1,4
     rem =  real(mm_in(lm(ii,jj)))
     imm = aimag(mm_in(lm(ii,jj)))
     hh = small*abs(rem)
     if (abs(imm).lt.hh) imm = -hh
     mm(ii) = cmplx(rem,imm,kind=kindc2)
   enddo
!
   icase = 0
   do ii=1,4
     if (ap(ii).gt.R0P0) icase = icase + base(ii)
   enddo
!
   if (icase.lt.15) then
! at least one exernal mass equal zero
     jcase = casetable(icase)
     if (jcase.eq.0.or.jcase.eq.1.or.jcase.eq.5) then
! two opposite masses equal zero
       a = pp(ll(5,icase)) - pp(ll(1,icase))
       c = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(3,icase))
       g = pp(ll(2,icase))
       h = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
       d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
       e = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
       f = mm(ll(4,icase))
       j = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
       rslt(0) = t13fun( a,c,g,h ,d,e,f,j )
     else
       a = pp(ll(3,icase))
       b = pp(ll(2,icase))
       c = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
       h = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(6,icase)) + pp(ll(2,icase))
       j = pp(ll(5,icase)) - pp(ll(1,icase)) - pp(ll(2,icase))
       d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
       e = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
       k = (mm(ll(1,icase)) - mm(ll(2,icase))) + pp(ll(6,icase)) - pp(ll(4,icase))
       f = mm(ll(4,icase))
       epk = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
       rslt(0) = tfun( a,b  ,c  ,h,j ,d,e  ,f,k ) &
               - tfun( a,b+j,c+h,h,j ,d,epk,f,k )
     endif
   else
! no extenal mass equal zero
     if    (real((pp(5)-pp(1)-pp(2))**2-4*pp(1)*pp(2)).gt.R0P0)then ;icase=0 !12, no permutation
     elseif(real((pp(6)-pp(2)-pp(3))**2-4*pp(2)*pp(3)).gt.R0P0)then ;icase=8 !23, 1 cyclic permutation
     elseif(real((pp(4)-pp(5)-pp(3))**2-4*pp(5)*pp(3)).gt.R0P0)then ;icase=4 !34, 2 cyclic permutations
     elseif(real((pp(4)-pp(1)-pp(6))**2-4*pp(1)*pp(6)).gt.R0P0)then ;icase=2 !41, 3 cyclic permutations
     else
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
         ,'no positive lambda, returning 0'
       return
     endif
     a = pp(ll(3,icase))
     b = pp(ll(2,icase))
     g = pp(ll(1,icase))
     c = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
     h = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(6,icase)) + pp(ll(2,icase))
     j = pp(ll(5,icase)) - pp(ll(1,icase)) - pp(ll(2,icase))
     d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
     e = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
     k = (mm(ll(1,icase)) - mm(ll(2,icase))) + pp(ll(6,icase)) - pp(ll(4,icase))
     f = mm(ll(4,icase))
     dpe = (mm(ll(2,icase)) - mm(ll(4,icase))) - pp(ll(6,icase))
     epk = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
     call solabc( x1,x2 ,sdnt ,g,j,b ,0 )
     if (aimag(sdnt).ne.R0P0) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
         ,'no real solution for alpha, returning 0'
       return
     endif
!BAD        if (abs(real(x1)).gt.abs(real(x2))) then
     if (abs(real(x1)).lt.abs(real(x2))) then !BETTER
       sdnt = x1
       x1 = x2
       x2 = sdnt
     endif
     o1 = C1P0-x1
     j1 = j+2*g*x1
     e1 = e+k*x1
     rslt(0) =   -tfun( a+b+c,g    ,j+h,c+2*b+(h+j)*x1,j1    ,dpe,k  ,f,e1 ) &
             + o1*tfun( a    ,b+g+j,c+h,c+h*x1        ,o1*j1 ,d  ,epk,f,e1 ) &
             + x1*tfun( a    ,b    ,c  ,c+h*x1        ,-j1*x1,d  ,e  ,f,e1 )
   endif
   end subroutine


   function t13fun( aa,cc,gg,hh ,dd,ee,ff,jj ) result(rslt)
!*******************************************************************
! /1   /x                             y
! | dx |  dy -----------------------------------------------------
! /0   /0    (gy^2 + hxy + dx + jy + f)*(ax^2 + cxy + dx + ey + f)
!
! jj should have negative imaginary part
!*******************************************************************
   complex(kindc2) ,intent(in) :: aa,cc,gg,hh ,dd,ee,ff,jj
   complex(kindc2) :: rslt ,kk,ll,nn,y1,y2,sdnt,ieps
   real(kindr2) ,parameter :: small=epsilon(R1P0)**2
!
!  write(*,*) 'MESSAGE from OneLOop t13fun: you are calling me' !CALLINGME
!
   ieps = cmplx(R0P0,small*abs(real(ff)),kind=kindc2)
!
   kk = hh*aa - cc*gg
   ll = aa*dd + hh*ee - dd*gg - cc*jj
   nn = dd*(ee - jj) + (hh - cc)*(ff-ieps)
   call solabc( y1,y2 ,sdnt ,kk,ll,nn ,0 )
!
   rslt = - s3fun( y1,y2 ,C0P0,C1P0 ,aa   ,ee+cc,dd+ff ) &
          + s3fun( y1,y2 ,C0P0,C1P0 ,gg   ,jj+hh,dd+ff ) &
          - s3fun( y1,y2 ,C0P0,C1P0 ,gg+hh,dd+jj,ff    ) &
          + s3fun( y1,y2 ,C0P0,C1P0 ,aa+cc,ee+dd,ff    )
!
   rslt = rslt/kk
   end function


   function t1fun( aa,cc,gg,hh ,dd,ee,ff,jj ) result(rslt)
!*******************************************************************
! /1   /x                         1
! | dx |  dy ----------------------------------------------
! /0   /0    (g*x + h*x + j)*(a*x^2 + c*xy + d*x + e*y + f)
!
! jj should have negative imaginary part
!*******************************************************************
   complex(kindc2) ,intent(in) :: aa,cc,gg,hh ,dd,ee,ff,jj
   complex(kindc2) ::rslt ,kk,ll,nn,y1,y2,sdnt,ieps
   real(kindr2) ,parameter :: small=epsilon(R1P0)**2
!
!  write(*,*) 'MESSAGE from OneLOop t1fun: you are calling me' !CALLINGME
!
   ieps = cmplx(R0P0,small*abs(real(ff)),kind=kindc2)
!
   kk = hh*aa - cc*gg
   ll = hh*dd - cc*jj - ee*gg
   nn = hh*(ff-ieps) - ee*jj
   call solabc( y1,y2 ,sdnt ,kk,ll,nn ,0 )
!
   rslt = - s3fun( y1,y2 ,C0P0,C1P0 ,aa+cc,dd+ee,ff ) &
          + s3fun( y1,y2 ,C0P0,C1P0 ,C0P0 ,gg+hh,jj ) &
          - s3fun( y1,y2 ,C0P0,C1P0 ,C0P0 ,gg   ,jj ) &
          + s3fun( y1,y2 ,C0P0,C1P0 ,aa   ,dd   ,ff )
!
   rslt = rslt/kk
   end function


   function tfun( aa,bb,cc ,gin,hin ,dd,ee,ff ,jin ) result(rslt)
!*******************************************************************
! /1   /x                             1
! | dx |  dy ------------------------------------------------------
! /0   /0    (g*x + h*x + j)*(a*x^2 + b*y^2 + c*xy + d*x + e*y + f)
!*******************************************************************
   complex(kindc2) ,intent(in) :: aa,bb,cc ,gin,hin ,dd,ee,ff ,jin
   complex(kindc2) :: rslt ,gg,hh,jj,zz(2),beta,tmpa(2),tmpb(2) &
                  ,tmpc(2),kiz(2),ll,nn,kk,y1,y2,yy(2,2),sdnt,ieps
   real(kindr2) :: sj,ab1,ab2,ac1,ac2,abab,acac,abac,det,ap1,ap2 &
                  ,apab,apac,x1(2,2),x2(2,2),xmin
   integer :: iz,iy,izmin
   logical :: pp(2,2),p1,p2
   real(kindr2) ,parameter :: small=epsilon(R1P0)**2
!
!  write(*,*) 'MESSAGE from OneLOop tfun: you are calling me' !CALLINGME
!
   sj = aimag(jin)
   if (sj.eq.R0P0) then
     sj = -R1P0
   else
     sj = sign(R1P0,aimag(jin))
   endif
   gg = -sj*gin
   hh = -sj*hin
   jj = -sj*jin
!
   if     (bb.eq.C0P0) then
     rslt = -sj*t1fun( aa,cc,gg,hh ,dd,ee,ff,jj )
     return
   elseif (aa.eq.C0P0) then
     rslt = -sj*t1fun( bb+cc,-cc,-gg-hh,gg,-dd-ee-2*(bb+cc),dd+cc,dd+ee+bb+cc+ff,gg+hh+jj )
     return
   endif
!
   ieps = cmplx(R0P0,small*abs(real(ff)),kind=kindc2)
!
   call solabc( zz(1),zz(2) ,sdnt ,bb,cc,aa ,0 )
   if (abs(zz(1)).gt.abs(zz(2))) then
     beta = zz(1)
     zz(1) = zz(2)
     zz(2) = beta
   endif
!
   do iz=1,2
     beta = zz(iz)
     tmpa(iz) = gg + beta*hh
     tmpb(iz) = cc + 2*beta*bb
     tmpc(iz) = dd + beta*ee
     kiz(iz) =        bb*tmpa(iz)               - hh*tmpb(iz)
     ll      =        ee*tmpa(iz) - hh*tmpc(iz) - jj*tmpb(iz)
     nn      = (ff-ieps)*tmpa(iz) - jj*tmpc(iz)
     call solabc( yy(iz,1),yy(iz,2) ,sdnt ,kiz(iz),ll,nn ,0 )
     if (abs(aimag(beta)).ne.R0P0) then
       ab1 =  real(-beta)
       ab2 = aimag(-beta)
       ac1 = ab1+R1P0 !real(C1P0-beta)
       ac2 = ab2     !aimag(C1P0-beta)
       abab = ab1*ab1 + ab2*ab2
       acac = ac1*ac1 + ac2*ac2
       abac = ab1*ac1 + ab2*ac2
       det = abab*acac - abac*abac
       do iy=1,2
         ap1 =  real(yy(iz,iy))
         ap2 = aimag(yy(iz,iy))
         apab = ap1*ab1 + ap2*ab2
         apac = ap1*ac1 + ap2*ac2
         x1(iz,iy) = ( acac*apab - abac*apac )/det
         x2(iz,iy) = (-abac*apab + abab*apac )/det
       enddo
     else
       do iy=1,2
         x1(iz,iy) = -R1P0
         x2(iz,iy) = -R1P0
       enddo
     endif
   enddo
   xmin = R1P0
   izmin = 2
   do iz=1,2
   do iy=1,2
     if ( x1(iz,iy).ge.R0P0.and.x2(iz,iy).ge.R0P0 &
                 .and.x1(iz,iy)+x2(iz,iy).le.R1P0 ) then
       pp(iz,iy) = .true.
       if (x1(iz,iy).lt.xmin) then
         xmin = x1(iz,iy)
         izmin = iz
       endif
       if (x2(iz,iy).lt.xmin) then
         xmin = x2(iz,iy)
         izmin = iz
       endif
     else
       pp(iz,iy) = .false.
     endif
   enddo
   enddo
   iz = izmin+1
   if (iz.eq.3) iz = 1
!
   beta = zz(iz)
   kk = kiz(iz)
   y1 = yy(iz,1)
   y2 = yy(iz,2)
   p1 = pp(iz,1)
   p2 = pp(iz,2)
!
   rslt = + s3fun( y1,y2 ,beta ,C1P0      ,C0P0    ,hh   ,gg+jj    ) &
          - s3fun( y1,y2 ,C0P0 ,C1P0-beta ,C0P0    ,gg+hh,   jj    ) &
          + s3fun( y1,y2 ,C0P0 ,    -beta ,C0P0    ,gg   ,   jj    ) &
          - s3fun( y1,y2 ,beta ,C1P0      ,bb      ,cc+ee,aa+dd+ff ) &
          + s3fun( y1,y2 ,C0P0 ,C1P0-beta ,aa+bb+cc,dd+ee,ff       ) &
          - s3fun( y1,y2 ,C0P0 ,    -beta ,aa      ,dd   ,ff       )
!
   sdnt = plnr( y1,y2 ,p1,p2, tmpa(iz),tmpb(iz),tmpc(iz) )
   if (aimag(beta).le.R0P0) then ;rslt = rslt + sdnt
   else                          ;rslt = rslt - sdnt
   endif
!
   rslt = -sj*rslt/kk
   end function


   function s3fun( y1i,y2i ,dd,ee ,aa,bb,cin ) result(rslt)
!*******************************************************************
! Calculate
!            ( S3(y1i) - S3(y2i) )/( y1i - y2i )
! where
!               /1    ee * ln( aa*x^2 + bb*x + cc )
!       S3(y) = |  dx -----------------------------
!               /0           ee*x - y - dd
!
! y1i,y2i should have a non-zero imaginary part
!*******************************************************************
   use avh_olo_logc ,only: logc
   complex(kindc2) ,intent(in) ::  y1i,y2i ,dd,ee ,aa,bb,cin
   complex(kindc2) :: rslt ,y1,y2,fy1y2,z1,z2,tmp,cc
   real(kindr2) ::rea,reb,rez1,rez2,imz1,imz2,simc
   real(kindr2) ,parameter :: small=epsilon(R1P0)**2
!
   if (ee.eq.C0P0) then
     rslt = C0P0
     return
   endif
!
   cc = cin
   rea = abs(aa)
   reb = abs(bb)
   simc = abs(cc)
   if (simc.lt.thrss3fun*min(rea,reb)) cc = C0P0
!
   simc = aimag(cc)
   if (simc.eq.R0P0) then
     simc = aimag(bb)
     if (simc.eq.R0P0) simc = -R1P0
   endif
   simc = sign(R1P0,simc)
!
   y1 = (dd+y1i)/ee
   y2 = (dd+y2i)/ee
   if (aimag(y1).eq.R0P0) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'y1 has zero imaginary part'
   endif
   if (aimag(y2).eq.R0P0) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'y2 has zero imaginary part'
   endif
   fy1y2 = r0fun( y1,y2 )
!
   if     (aa.ne.C0P0) then
!
     call solabc( z1,z2 ,tmp ,aa,bb,cc ,0 )
     rea  = sign(R1P0,real(aa))
     rez1 = real(z1)
     rez2 = real(z2) 
     imz1 = aimag(z1) ! sign(Im(a*z1*z2)) = simc
     imz2 = aimag(z2)
     if (imz1.eq.R0P0) imz1 = simc*rea*sign(R1P0,rez2)*abs(small*rez1)
     if (imz2.eq.R0P0) imz2 = simc*rea*sign(R1P0,rez1)*abs(small*rez2)
     z1 = cmplx( rez1,imz1,kind=kindc2 )
     z2 = cmplx( rez2,imz2,kind=kindc2 )
     rslt = fy1y2 * ( logc(qonv(aa,simc)) &
                    + eta3( -z1,-imz1,-z2,-imz2,C0P0,simc*rea ) ) &
          + r1fun( z1,y1,y2,fy1y2 ) &
          + r1fun( z2,y1,y2,fy1y2 )
!
   elseif (bb.ne.C0P0) then
!
     z1 = -cc/bb ! - i|eps|Re(b)
     reb  = real(bb)
     rez1 = real(z1)
     imz1 = aimag(z1)
     if (abs(imz1).eq.R0P0) then
       imz1 = -simc*reb*abs(small*rez1/reb)
       z1 = cmplx( rez1,imz1,kind=kindc2 )
     endif
     rslt = fy1y2 * ( logc(qonv(bb,simc)) &
                    + eta3(bb,simc ,-z1,-imz1 ,cc,simc) ) &
          + r1fun( z1,y1,y2,fy1y2 )
!
   elseif (cc.ne.C0P0) then
!
     rslt = logc( qonv(cc,simc) )*fy1y2
!
   else!if (aa=bb=cc=0)
!
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'cc equal zero, returning 0'
     rslt = C0P0
!
   endif
!
   rslt = rslt/ee
   end function


   function r1fun( zz,y1,y2,fy1y2 ) result(rslt)
!*******************************************************************
! calculates  ( R1(y1,z) - R1(y2,z) )/( y1 - y2 )
! where
!                          /     / 1-y \       / 1-z \ \
!      R1(y,z) = ln(y-z) * | log |-----| - log |-----| |
!                          \     \ -y  /       \ -z  / / 
!
!                      /    y-z \       /    y-z \
!                - Li2 |1 - ----| + Li2 |1 - ----|
!                      \    -z  /       \    1-z /
!
!                                     / 1-y1 \       / 1-y2 \
!                                 log |------| - log |------| 
! input fy1y2 should be equal to      \  -y1 /       \  -y2 /
!                                 ---------------------------
!                                           y1 - y2
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_li2c ,only: li2c
   use avh_olo_logc2 ,only: logc2
   use avh_olo_li2c2 ,only: li2c2
   complex(kindc2) ,intent(in) :: y1,y2,zz,fy1y2
   complex(kindc2) :: rslt ,oz
   type(qmplx_type) :: q1z,q2z,qq
   real(kindr2) :: h12,hz1,hz2,hzz,hoz
   logical :: zzsmall,ozsmall
!
   oz = C1P0-zz
   h12 = abs(y1-y2)
   hz1 = abs(y1-zz)
   hz2 = abs(y2-zz)
   hzz = abs(zz)
   hoz = abs(oz)
   q1z = qonv(y1-zz)
   q2z = qonv(y2-zz)
!
   zzsmall = .false.
   ozsmall = .false.
   if     (hzz.lt.hz1.and.hzz.lt.hz2.and.hzz.lt.hoz) then ! |z| < |y1-z|,|y2-z|
     zzsmall = .true.
     rslt = fy1y2*logc( q1z ) &
          - ( logc(q1z*q2z)/2 + logc(qonv((y2-C1P0)/y2)) &
                                        - logc(qonv(oz)) )*logc2(q1z/q2z)/(y2-zz)
   elseif (hoz.lt.hz1.and.hoz.lt.hz2) then ! |1-z| < |y1-z|,|y2-z|
     ozsmall = .true.
     rslt = fy1y2*logc( q1z ) &
          - (-logc(q1z*q2z)/2 + logc(qonv((y2-C1P0)/y2)) &
                                       + logc(qonv(-zz)) )*logc2(q1z/q2z)/(y2-zz)
   elseif (h12.le.hz2.and.hz2.le.hz1) then ! |y1-y2| < |y2-z| < |y1-z|
     rslt = fy1y2*logc( q1z ) - r0fun( y2,zz )*logc2( q1z/q2z )        
   elseif (h12.le.hz1.and.hz1.le.hz2) then ! |y1-y2| < |y2-z| < |y1-z|
     rslt = fy1y2*logc( q2z ) - r0fun( y1,zz )*logc2( q2z/q1z )        
   else!if(hz1.lt.h12.or.hz2.lt.h12) then ! |y2-z|,|y1-z| < |y1-y2|
     rslt = C0P0
     if (hz1.ne.R0P0) rslt = rslt + (y1-zz)*logc( q1z )*r0fun( y1,zz )
     if (hz2.ne.R0P0) rslt = rslt - (y2-zz)*logc( q2z )*r0fun( y2,zz )
     rslt = rslt/(y1-y2)
   endif
!
   if (zzsmall) then ! |z| < |y1-z|,|y2-z|
     qq  = qonv(-zz)
     rslt = rslt + ( li2c( qq/q1z ) - li2c( qq/q2z ) )/(y1-y2)
   else
     qq  = qonv(-zz)
     rslt = rslt + li2c2( q1z/qq ,q2z/qq )/zz
   endif
!
   if (ozsmall) then ! |1-z| < |y1-z|,|y2-z|
     qq  = qonv(oz)
     rslt = rslt - ( li2c( qq/q1z ) - li2c( qq/q2z ) )/(y1-y2)
   else
     qq = qonv(oz)
     rslt = rslt + li2c2( q1z/qq ,q2z/qq )/oz
   endif
   end function


   function r0fun( y1,y2 ) result(rslt)
!*******************************************************************
!      / 1-y1 \       / 1-y2 \
!  log |------| - log |------| 
!      \  -y1 /       \  -y2 /
!  ---------------------------
!            y1 - y2
!
! y1,y2 should have non-zero imaginary parts
!*******************************************************************
   use avh_olo_logc2 ,only: logc2
   complex(kindc2) ,intent(in) :: y1,y2
   complex(kindc2) :: rslt ,oy1,oy2
   oy1 = C1P0-y1
   oy2 = C1P0-y2
   rslt = logc2( qonv(-y2)/qonv(-y1) )/y1 &
        + logc2( qonv(oy2)/qonv(oy1) )/oy1
   end function


   function plnr( y1,y2 ,p1,p2 ,aa,bb,cc ) result(rslt)
!*******************************************************************
!                   /   a    \          /   a    \
!            p1*log |--------| - p2*log |--------| 
!                   \ b*y1+c /          \ b*y2+c /
! 2*pi*imag* -------------------------------------
!                           y1 - y2
! 
! p1,p2 are logical, to be interpreted as 0,1 in the formula above 
!*******************************************************************
   use avh_olo_logc ,only: logc
   use avh_olo_logc2 ,only: logc2
   complex(kindc2) ,intent(in) :: y1,y2 ,aa,bb,cc
   logical         ,intent(in) :: p1,p2
   complex(kindc2) :: rslt ,x1,x2,xx
   type(qmplx_type) :: q1,q2
   complex(kindc2) ,parameter :: twopii=CiP0*TWOPI
!
   if (p1) then
     x1 = bb*y1 + cc
     xx = aa/x1
     if (aimag(xx).eq.R0P0) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop plnr: ' &
         ,'aa/x1 has zero imaginary part'
     endif
     q1 = qonv(xx)
   endif
   if (p2) then
     x2 = bb*y2 + cc
     xx = aa/x2
     if (aimag(xx).eq.R0P0) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop plnr: ' &
         ,'aa/x2 has zero imaginary part'
     endif
     q2 = qonv(xx)
   endif
   if (p1) then
     if (p2) then
       rslt = logc2( q2/q1 ) * twopii*bb/x2
     else
       rslt = logc( q1 ) * twopii/(y1-y2)
     endif
   elseif (p2) then
     rslt = logc( q2 ) * twopii/(y2-y1) ! minus sign
   else
     rslt = C0P0
   endif
   end function


end module
