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


module avh_olo_bub
  use avh_olo_kinds
  use avh_olo_units
  implicit none
  private
  public :: init_bub ,tadp ,bub0 ,bub11
  integer ,parameter :: d=kindr2
  integer ,parameter :: ntrmmax=20
  real(kindr2) ,parameter :: thrslistd(ntrmmax)=&
    (/5e-5_d,5e-3_d,0.05_d,0.10_d,0.15_d,0.20_d,0.30_d,0.40_d &
     ,0.50_d,0.60_d,0.65_d,0.68_d,0.72_d,0.74_d,0.76_d,0.78_d &
     ,0.80_d,0.82_d,0.83_d,0.84_d/)
  real(kindr2) ,parameter :: thrslisth(ntrmmax)=&
    (/7e-8_d,5e-4_d,2e-3_d,1e-2_d,3e-2_d,6e-2_d,0.11_d,0.17_d &
     ,0.22_d,0.28_d,0.33_d,0.37_d,0.42_d,0.47_d,0.51_d,0.54_d &
     ,0.58_d,0.60_d,0.62_d,0.65_d/)
  real(kindr2) ,parameter :: thrslistq(ntrmmax)=&
    (/1e-10_d,5e-5_d,1e-4_d,1e-3_d,7e-3_d,0.02_d,0.04_d,0.07_d &
      ,0.10_d,0.13_d,0.17_d,0.20_d,0.25_d,0.30_d,0.34_d,0.38_d &
      ,0.42_d,0.44_d,0.47_d,0.50_d/)
  real(kindr2) :: thrs=0.07_d
  real(kindr2) :: thrsexp=0.01_d
  real(kindr2) :: thrslist(1:ntrmmax)=thrslistd(1:ntrmmax)
  integer      :: ntrm=11,nnexp=7
  complex(kindc2) :: aaexp(8)=C0P0
  integer :: ndigits=0
contains
!
  subroutine init_bub(ndig)
  integer ,intent(in) :: ndig
  integer :: ii
  if (ndigits.eq.ndig) return ;ndigits=ndig
  if     (ndigits.lt.16) then
    thrs = 0.07_kindr2    ! double precision,
    ntrm = 11             ! tested to suit also b11
    thrsexp = 0.01_kindr2 !
    nnexp = 7             !
    thrslist = thrslistd  ! double precision
  elseif (ndigits.lt.24) then
    thrs = 0.02_kindr2     ! guess
    ntrm = 11              !
    thrsexp = 0.001_kindr2 !
    nnexp = 7              !
    thrslist = thrslisth   !
  else
    thrs = 0.005_kindr2     ! quadruple precision, not tested
    ntrm = 11               !
    thrsexp = 0.0001_kindr2 !
    nnexp = 7               !
    thrslist = thrslistq    ! quadruple precision
  endif
  do ii=1,nnexp
    aaexp(ii) = C1P0/(ii*(ii+1))
  enddo
  end subroutine


  subroutine tadp( rslt ,mm ,amm ,rmu2 )
!*******************************************************************
! The 1-loop scalar 1-point function.
!*******************************************************************
  use avh_olo_func
  use avh_olo_logc ,only : logc
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: mm
  real(kindr2)    ,intent(in)  :: amm,rmu2
!
!  write(*,*) 'MESSAGE from OneLOop tadp: you are calling me' !CALLINGME
!
  rslt(2) = C0P0
  if (amm.eq.R0P0) then
    rslt(1) = C0P0
    rslt(0) = C0P0
  else
    rslt(1) = mm
    rslt(0) = mm - mm*logc( qonv(mm/rmu2,-1) )
  endif
  end subroutine


  subroutine bub0( rslt ,pp,m1i,m2i ,app,am1i,am2i ,rmu2 )
!*******************************************************************
! The 1-loop scalar 2-point function. Based on the formulas from
! A. Denner, Fortsch.Phys.41:307-420,1993 arXiv:0709.1075 [hep-ph]
!*******************************************************************
  use avh_olo_func
  use avh_olo_logc ,only: logc
  use avh_olo_logc2 ,only: logc2
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: pp,m1i,m2i
  real(kindr2)    ,intent(in)  :: app,am1i,am2i,rmu2
  complex(kindc2) :: cc(0:ntrmmax),m1,m2,hh,aa,bb,rr,dd
  type(qmplx_type) :: qmm,qz1
  complex(kindc2) ,parameter :: two=C1P0*2
  real(kindr2) :: am1,am2,tt
  integer :: ii
!
!  write(*,*) 'MESSAGE from OneLOop bub0: you are calling me' !CALLINGME
!
  tt = max(am1i,am2i)
  if (am1i.lt.tt) then
    m1=m1i ;am1=am1i
    m2=m2i ;am2=am2i
  else
    m1=m2i ;am1=am2i
    m2=m1i ;am2=am1i
  endif
!
  rslt(2) = C0P0     
  rslt(1) = C1P0
!
  if (am2.eq.R0P0) then
    if (app.eq.R0P0) then
      rslt(1) = C0P0     
      rslt(0) = C0P0     
    else
      rslt(0) = two - logc(qonv(-pp/rmu2,-1))
    endif
  else!if(am2.ne.R0P0)
    tt = app/tt
    if (am1.eq.R0P0) then
      qmm = qonv(m2/rmu2,-1)
      if     (app.eq.R0P0) then
        rslt(0) = C1P0 - logc(qmm)
      elseif (pp.eq.m2) then
        rslt(0) = two - logc(qmm)
      elseif (tt.lt.R1P0) then
        hh = m2-pp
        rslt(0) = two + (hh/pp)*logc(qonv(hh/rmu2,-1)/qmm) - logc(qmm)
      else!if (tt.ge.R1P0) then
        hh = m2-pp
        rslt(0) = two - (m2/pp)*logc(qmm) + (hh/pp)*logc(qonv(hh/rmu2,-1))
      endif
    else!if(am1.ne.R0P0)
      if (app.eq.R0P0) then
         qz1 = qonv(m1/rmu2,-1)
         rslt(0) = C1P0 + logc2(qz1/qonv(m2/rmu2,-1)) - logc(qz1)
      else!if(pp.ne.C0P0)
        if     (tt.le.thrs) then
          call expans( cc ,m1,m2 ,am1,am2 ,rmu2 ,ntrm)
          rslt(0) = cc(ntrm)
          do ii=ntrm-1,0,-1
            rslt(0) = cc(ii) + pp*rslt(0)
          enddo
        elseif (tt.lt.R1P0) then
          hh = mysqrt(m1)
          bb = mysqrt(m2)
          aa = hh*bb ! sm1*sm2
          bb = hh/bb ! sm1/sm2
          hh = (m1+m2-pp)/aa
          dd = (m2-m1)**2 + ( pp - 2*(m1+m2) )*pp
          dd = mysqrt(dd)/aa
          call rfun0( rr ,dd ,hh )
          qz1 = qonv(bb,-1) ! sm1/sm2
          rslt(0) = two - logc(qonv(m2/rmu2,-1)) &
                        + logc(qz1)*two*m1/(aa*rr-m1) &
                        + logc2(qz1*qonv(rr,-1))*dd*aa/(aa*rr-m1+pp)
        else
          hh = mysqrt(m1)
          bb = mysqrt(m2)
          aa = hh*bb ! sm1*sm2
          bb = hh/bb ! sm1/sm2
          hh = (m1+m2-pp)/aa
          call rfun( rr,dd ,hh )
          rslt(0) = two - logc(qonv(aa/rmu2,-1)) &
                  + (logc(qonv(bb,-1))*(m2-m1) + logc(qonv(rr,-1))*dd*aa)/pp
        endif
!        call expans( cc ,m1,m2 ,am1,am2 ,rmu2 ,ntrm) !DEBUG
!        hh = cc(ntrm)                            !DEBUG
!        do ii=ntrm-1,0,-1                        !DEBUG
!          hh = cc(ii) + pp*hh                    !DEBUG
!        enddo                                    !DEBUG
!        write(*,'(a4,2d24.16)') 'exp:',hh        !DEBUG
      endif
    endif
  endif
  end subroutine


  subroutine bub11( b11,b00,b1,b0 ,pp,m0,m1 ,app,am0,am1 ,rmu2 )
!*******************************************************************
! Return the Passarino-Veltman functions b11,b00,b1,b0 , for
!
!      C   /      d^(Dim)q
!   ------ | -------------------- = b0
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!      C   /    d^(Dim)q q^mu
!   ------ | -------------------- = p^mu b1
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!      C   /  d^(Dim)q q^mu q^nu
!   ------ | -------------------- = g^{mu,nu} b00 + p^mu p^nu b11
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!*******************************************************************
  complex(kindc2) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindc2) ,intent(in)  :: pp,m0,m1
  real(kindr2)    ,intent(in)  :: app,am0,am1,rmu2
  complex(kindc2) :: a1(0:2),a0(0:2),ff,gg,c1,c2,cc(0:ntrmmax)
  real(kindr2) :: rr,maxm
  integer :: ii
!
  maxm = max(am0,am1)
  if (maxm.eq.R0P0) then
    if (app.eq.R0P0) then
      b0  = C0P0
      b1  = C0P0
      b00 = C0P0
      b11 = C0P0
      return
    endif
    rr = R1P0+thrs
  else
    rr = app/maxm
  endif
!
  ff = pp - m1 + m0
  gg = m0 + m1 - pp/3
  b0(2)  = C0P0
  b1(2)  = C0P0
  b00(2) = C0P0
  b11(2) = C0P0
  b0(1)  = C1P0
  b1(1)  = -C1P0/2
  b00(1) = gg/4
  b11(1) = C1P0/3
  call tadp( a1 ,m0 ,am0 ,rmu2 )
  call tadp( a0 ,m1 ,am1 ,rmu2 )
!
  if (rr.le.thrs) then
!    write(*,*) 'expansion' !DEBUG
    call expans( cc ,m0,m1 ,am0,am1 ,rmu2 ,ntrm )
    c2 = cc(ntrm)
    do ii=ntrm-1,2,-1
      c2 = cc(ii) + pp*c2
    enddo
    c1 = cc(1) + pp*c2
    b0(0)  = cc(0) + pp*c1
    b1(0)  = -( cc(0) + ff*c1 )/2
    b00(0) = ( a0(0) + ff*b1(0) + 2*m0*b0(0) + gg )/6
    b11(0) = cc(0) + (ff+m0-m1)*cc(1) + ff*ff*c2 - m0*c1
    b11(0) = ( b11(0) + C1P0/6 )/3
  else
    call bub0( b0 ,pp,m0,m1 ,app,am0,am1 ,rmu2 )
    b1(0)  = ( a1(0) - a0(0) - ff*b0(0) )/(2*pp)
    b00(0) = ( a0(0) + ff*b1(0) + 2*m0*b0(0) + gg )/6
    b11(0) = ( a0(0) - 2*ff*b1(0) - m0*b0(0) - gg/2 )/(3*pp)
  endif
!
  end subroutine


  subroutine expans( cc ,m1i,m2i ,am1i,am2i ,rmu2 ,ntrm )
!*******************************************************************
! Returns the first 1+ntrm coefficients of the expansion in p^2 of
! the finite part of the 1-loop scalar 2-point function 
!*******************************************************************
  use avh_olo_func
  use avh_olo_logc ,only: logc
  integer         ,intent(in)  :: ntrm
  complex(kindc2) ,intent(out) :: cc(0:ntrm)
  complex(kindc2) ,intent(in)  :: m1i,m2i
  real(kindr2)    ,intent(in)  :: am1i,am2i,rmu2
  complex(kindc2) :: m1,m2,zz,oz,xx,logz,tt(ntrm)
  type(qmplx_type) :: qm1,qm2,qzz
  real(kindr2) :: am1,am2
  integer :: ii
  real(kindr2) ::rr
!
!  write(*,*) 'MESSAGE from OneLOop bub expans: you are calling me' !CALLINGME
!
  if (am1i.lt.am2i) then
    m1=m1i ;am1=am1i
    m2=m2i ;am2=am2i
  else
    m1=m2i ;am1=am2i
    m2=m1i ;am2=am1i
  endif
!
  if (am2.eq.R0P0) then
!
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop bub expans: ' &
      ,'m1=m2=0, returning 0'
    do ii=0,ntrm
      cc(ii) = C0P0
    enddo
!
  else
!
    qm1 = qonv(m1/rmu2,-1) 
    qm2 = qonv(m2/rmu2,-1)
    qzz = qm1/qm2
    if (mod(qzz%p,2).eq.0) then
      zz = qzz%c
    else
      zz = -qzz%c
    endif
!
    if (m1.eq.C0P0) then
      cc(0) = C1P0 - logc(qm2)
    else
      oz = C1P0-zz
      rr = abs(oz)
      if (rr.lt.thrsexp) then
        xx = aaexp(nnexp)
        do ii=nnexp-1,1,-1
          xx = aaexp(ii) + oz*xx
        enddo
        xx = oz*xx
      else
        logz = logc( qzz )
        xx = zz*logz + oz
        xx = xx/oz
      endif
      cc(0) = xx - logc(qm2)
    endif
!
    zz = C1P0-zz
    xx = C1P0
    call expans1(tt ,ntrm,zz)
    do ii=1,ntrm
      xx = xx*m2
      cc(ii) = tt(ii)/(ii*xx)
    enddo
!
  endif
  end subroutine


  subroutine expans1(tt ,ntrm,zz)
!*******************************************************************
! Returns  tt(n) = int( ( x*(1-x)/(1-zz*x) )^n , x=0..1 )
! for  n=1...ntrm  and  |zz|=<1
!
! Gives at least 2 correct digits (4 at quad.) for tt(ntrm),
! and increasingly more digits for tt(i<ntrm)
!
! Uses recursion on integrals of the type
!    int( x^m * (1-x)^n / (1-z*x)^n , x=0..1 )
! and
!    int( x^m * (1-x)^n / (1+y*x)^(n+2) , x=0..1 )
! where  y = z/(1-z)
! The latter integrals are related to the original ones via the
! substitution  x <- 1-x  followed by  x <- (1-x)/(1+y*x)
!*******************************************************************
  integer         ,intent(in)  :: ntrm
  complex(kindc2) ,intent(in)  :: zz
  complex(kindc2) ,intent(out) :: tt(ntrm)
  complex(kindc2) :: tu(ntrm),tv(ntrm) ,tt0,tu0,tv0,yy,y2,oy
  real(kindr2) :: rr
  integer :: nn,ii,jj
!
  rr = real(zz)
  nn = ntrm
  if (nn.lt.1) nn = 1
  if (nn.gt.ntrmmax) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop bub expans1: ' &
      ,'ntrm =',nn,' > nmax =',ntrmmax,', using ntrm=nmax'
    nn = ntrmmax
    do ii=nn+1,ntrm
      tt(ii) = C0P0
    enddo
  endif
!
  if (zz.eq.C1P0) then
!    write(*,'(a16,i4)') 'simple expansion',nn !DEBUG
    do ii=1,nn
      tt(ii) = cmplx(R1P0/(ii+1),kind=kindc2)
    enddo
  elseif (rr.lt.thrslist(nn)) then
! Backward recursion, number of correct decimals constant, so need
! full precision from the start
!    write(*,'(a8,i4,d24.16)') 'Backward',nn,rr !DEBUG
    call expans2(tt(nn),tv(nn-1),tu(nn-1) ,nn,zz)
    do ii=nn-1,2,-1
      jj = ii+1
      tt(ii  ) = 2*tv(ii) - (zz*tt(jj)*ii)/jj
      tu(ii-1) = (2+R1P0/ii)*tt(ii) - zz*tu(ii)
      tv(ii-1) = (C1P0-zz)*tu(ii-1) + zz*( 2*tt(ii) - zz*tu(ii) )
    enddo
    tt(1) = 2*tv(1) - zz*tt(2)/2
  else
! Foreward recursion, number of correct decimals decreases
!    write(*,'(a8,i4,d24.16)') 'Foreward',nn,rr !DEBUG
    yy = zz/(C1P0-zz)
    y2 = yy*yy
    oy = C1P0+yy ! C1P0/(C1P0-zz)
    tt0 = C1P0-zz ! 1/(1+y)
    tu0 = ( oy*log(oy)-yy )/( y2*oy )
    tv0 = tt0/2
    tt(1) = ( tt0-2*tu0 )/( 2*yy )
    tv(1) = ( tv0 - 3*tt(1) )/( 3*yy )
    tu(1) = ( oy*tu0 - 2*yy*tt(1) - tv0 )/y2
    do ii=2,nn
      jj = ii-1
      tt(ii) = ii*( tt(jj)-2*tu(jj) )/( (ii+1)*yy )
      tv(ii) = ( ii*tv(jj) - (ii+ii+1)*tt(ii) )/( (ii+2)*yy )
      tu(ii) = ( oy*tu(jj) - 2*yy*tt(ii) - tv(jj) )/y2
    enddo
    yy = oy
    do ii=1,nn
      oy = oy*yy
      tt(ii) = oy*tt(ii)
    enddo
  endif
  end subroutine


  subroutine expans2(ff,fa,fb ,nn_in,zz)
!*******************************************************************
! ff = Beta(nn+1,nn+1) * 2F1(nn  ,nn+1;2*nn+2;zz)
! fa = Beta(nn+1,nn  ) * 2F1(nn-1,nn+1;2*nn+1;zz)
! fb = Beta(nn  ,nn+1) * 2F1(nn  ,nn  ;2*nn+1;zz)
!*******************************************************************
  complex(kindc2) ,intent(out) :: ff,fa,fb
  complex(kindc2) ,intent(in)  :: zz
  integer         ,intent(in)  :: nn_in
  integer ,parameter :: nmax=100
  integer :: aa,bb,cc,ii,ntrm
  complex(kindc2) ,save :: qq(0:nmax),qa(0:nmax),qb(0:nmax),gg,ga
  real(kindr2) ,save :: logprec=-36.0_kindr2
  integer ,save :: nn=0
  real(kindr2) :: ac0,bc0,ai,bi,ci,ac,bc
  if (nn.ne.nn_in) then
    nn = nn_in
    aa = nn-1
    bb = nn
    cc = nn+nn+1
    qq(0) = C1P0
    qa(0) = C1P0
    qb(0) = C1P0
    ac0 = real(aa,KIND=kindr2)/real(cc,KIND=kindr2)
    bc0 = real(bb,KIND=kindr2)/real(cc,KIND=kindr2)
    ntrm = nmax
    do ii=1,ntrm
      ai = real(aa+ii,KIND=kindr2)
      bi = real(bb+ii,KIND=kindr2)
      ci = real(cc+ii,KIND=kindr2)
      ac = ai/ci
      bc = bi/ci
      qq(ii) = qq(ii-1) * ai*bc  / ii
      qa(ii) = qa(ii-1) * ac0*bi / ii
      qb(ii) = qb(ii-1) * ai*bc0 / ii
      ac0 = ac
      bc0 = bc
    enddo
    ai = R1P0
    do ii=2,nn-1
      ai = ai*ii
    enddo
    ci = ai
    cc = nn+nn
    do ii=nn,cc
      ci = ci*ii
    enddo
    bi = ai*nn
    gg = bi*bi/(ci*(cc+1))
    ga = ai*bi/ci
    logprec = log(epsilon(R1P0))
  endif
!
  ai = abs(zz)
  if (ai.gt.R0P0) then
    ntrm = 1 + int(logprec/log(ai))
  else
    ntrm = 1
  endif
  if (ntrm.gt.nmax) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop bub expans2: ' &
      ,'ntrm =',ntrm,' > nmax =',nmax,', putting ntrm=nmax'
    ntrm = nmax
  endif
!
  ff = qq(ntrm)
  fa = qa(ntrm)
  fb = qb(ntrm)
  do ii=ntrm-1,0,-1
    ff = qq(ii) + ff*zz
    fa = qa(ii) + fa*zz
    fb = qb(ii) + fb*zz
  enddo
  ff = gg*ff
  fa = ga*fa
  fb = ga*fb
  end subroutine
!
end module
