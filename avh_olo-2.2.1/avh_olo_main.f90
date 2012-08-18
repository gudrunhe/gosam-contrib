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


module avh_olo
  use avh_olo_kinds
  use avh_olo_units
  use avh_olo_print
!
  implicit none
  private
  public :: olo_kind ,olo_unit ,olo_scale ,olo_onshell ,olo_setting &
           ,olo_a0 ,olo_b0 ,olo_b11 ,olo_c0 ,olo_d0
!
  integer ,parameter :: olo_kind = kindr2
!
  integer      :: ndigits = 0        ! corrected in subroutine hello
  real(kindr2) :: onshellthrs = R0P0 ! corrected in subroutine hello
  logical      :: nonzerothrs = .false.
!
  real(kindr2) :: muscale = R1P0
!
  character(99) ,parameter :: warnonshell=&
       'it seems you forgot to put some input explicitly on shell. ' &
     //'You may  call olo_onshell  to cure this.'
!
  logical :: intro=.true.
!
  interface olo_a0
    module procedure loc_a0r,a0rr,loc_a0c,a0cr
  end interface 
  interface olo_b0
    module procedure b0rr,b0rrr,b0rc,b0rcr,b0cc,b0ccr
  end interface 
  interface olo_b11
    module procedure b11rr,b11rrr,b11rc,b11rcr,b11cc,b11ccr
  end interface 
  interface olo_c0
    module procedure c0rr,c0rrr,c0rc,c0rcr,c0cc,c0ccr
  end interface 
  interface olo_d0
    module procedure d0rr,d0rrr,d0rc,d0rcr,d0cc,d0ccr
  end interface 

contains
 
 
  subroutine hello
!*******************************************************************
!*******************************************************************
  use avh_olo_loga2 ,only: init_loga2
  use avh_olo_li2c2 ,only: init_li2c2
  use avh_olo_bub   ,only: init_bub
  use avh_olo_boxc  ,only: init_boxc
!
  intro = .false.
!
  write(*,'(a72)') '########################################################################'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '#                     You are using OneLOop-2.2.1                      #'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '# for the evaluation of 1-loop scalar 1-, 2-, 3- and 4-point functions #'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '# author: Andreas van Hameren <hamerenREMOVETHIS@ifj.edu.pl>           #'
  write(*,'(a72)') '#   date: 07-09-2011                                                   #'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '# Please cite                                                          #'
  write(*,'(a72)') '#    A. van Hameren,                                                   #'
  write(*,'(a72)') '#      Comput.Phys.Commun. 182 (2011) 2427-2438, arXiv:1007.4716       #'
  write(*,'(a72)') '#    A. van Hameren, C.G. Papadopoulos and R. Pittau,                  #'
  write(*,'(a72)') '#      JHEP 0909:106,2009, arXiv:0903.4665                             #'
  write(*,'(a72)') '# in publications with results obtained with the help of this program. #'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '########################################################################'
!
  ndigits = int(digits(R1P0)*log(radix(R1P0)*R1P0)/log(R1P0*10))
      if (ndigits.lt.16) then ;onshellthrs = epsilon(R1P0)*100
  elseif (ndigits.lt.24) then ;onshellthrs = epsilon(R1P0)*1000
  elseif (ndigits.lt.32) then ;onshellthrs = epsilon(R1P0)*10000
  else                        ;onshellthrs = epsilon(R1P0)*1000000
  endif
!
  call init_print( ndigits )
  call init_loga2( ndigits )
  call init_li2c2( ndigits )
  call init_bub(   ndigits )
  call init_boxc(  ndigits )
!
  end subroutine
 
 
  subroutine olo_unit( val ,message )
!*******************************************************************
!*******************************************************************
  integer     ,intent(in) :: val
  character(*),intent(in),optional :: message
  if (intro) call hello
  if (present(message)) then ;call set_unit( message ,val )
  else                       ;call set_unit( 'all'   ,val )
  endif
  end subroutine
 
 
  subroutine olo_scale( val )
!*******************************************************************
!*******************************************************************
  real(kindr2) ,intent(in) :: val
  if (intro) call hello
  muscale = val
  end subroutine
 
 
  subroutine olo_onshell( thrs )
!*******************************************************************
!*******************************************************************
  real(kindr2) ,intent(in) :: thrs
  if (intro) call hello
  nonzerothrs = .true.
  onshellthrs = thrs
  end subroutine


  subroutine olo_setting( iunit )
!*******************************************************************
!*******************************************************************
  integer,optional,intent(in) :: iunit
  integer :: nunit
  if (intro) call hello
  nunit = munit
  if (present(iunit)) nunit = iunit
  if (nunit.le.0) return
!
  write(nunit,*) 'MESSAGE from OneLOop: real kind parameter =',trim(myprint(kindr2))
  write(nunit,*) 'MESSAGE from OneLOop: significant digits =',trim(myprint(ndigits))
!
  if (nonzerothrs) then
    write(nunit,*) 'MESSAGE from OneLOop: on-shell threshold =',trim(myprint(onshellthrs))
  else
    write(nunit,*) 'MESSAGE from OneLOop: on-shell threshold is not set'
  endif
!
  write(nunit,*) 'MESSAGE from OneLOop: scale (mu, not mu^2) =',trim(myprint(muscale))
!
  end subroutine
 
 
!*******************************************************************
!
!           C   / d^(Dim)q
! rslt = ------ | -------- 
!        i*pi^2 / (q^2-mm)
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  mm = mass squared
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************
!
  subroutine loc_a0r( rslt ,mm  )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: tadp
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: mm
!
  include 'avh_olo_a0_a.h90'
!
  mulocal = muscale
!
  include 'avh_olo_a0_b.h90'
!
  end subroutine
!
  subroutine a0rr( rslt ,mm ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: tadp
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: mm
  real(kindr2)    ,intent(in)  :: rmu
!
  include 'avh_olo_a0_a.h90'
!
  mulocal = rmu
!
  include 'avh_olo_a0_b.h90'
!
  end subroutine
!
  subroutine loc_a0c( rslt ,mm )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: tadp
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: mm
!
  include 'avh_olo_a0_a.h90'
!
  mulocal = muscale
!
  include 'avh_olo_a0_b.h90'
!
  end subroutine
!
  subroutine a0cr( rslt ,mm ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: tadp
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: mm
  real(kindr2)    ,intent(in)  :: rmu
!
  include 'avh_olo_a0_a.h90'
!
  mulocal = rmu
!
  include 'avh_olo_a0_b.h90'
!
  end subroutine


!*******************************************************************
!
!           C   /      d^(Dim)q
! rslt = ------ | --------------------
!        i*pi^2 / [q^2-m1][(q+k)^2-m2]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  pp = k^2, m1,m2 = mass squared
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************
!
  subroutine b0rr( rslt ,pp ,m1,m2 )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub0
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: pp,m1,m2
!
  include 'avh_olo_b0_a.h90'
!
  mulocal = muscale
!
  app = abs(pp)
  am1 = abs(m1)
  am2 = abs(m2)
!
  include 'avh_olo_b0_b.h90'
!
  end subroutine
!
  subroutine b0rrr( rslt ,pp ,m1,m2 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub0
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: pp,m1,m2,rmu
!
  include 'avh_olo_b0_a.h90'
!
  mulocal = rmu
!
  app = abs(pp)
  am1 = abs(m1)
  am2 = abs(m2)
!
  include 'avh_olo_b0_b.h90'
!
  end subroutine
!
  subroutine b0rc( rslt ,pp ,m1,m2 )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub0
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: pp
  complex(kindc2) ,intent(in)  :: m1,m2
  real(kindr2) :: hh
!
  include 'avh_olo_b0_a.h90'
!
  mulocal = muscale
!
  app = abs(pp)
!
  am1 = real(r1)
  hh  = aimag(r1)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = cmplx( am1 ,-hh ,kind=kindc2 )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = real(r2)
  hh  = aimag(r2)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = cmplx( am2 ,-hh ,kind=kindc2 )
  endif
  am2 = abs(am2) + abs(hh)
!
  include 'avh_olo_b0_b.h90'
!
  end subroutine
!
  subroutine b0rcr( rslt ,pp,m1,m2 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub0
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: pp ,rmu
  complex(kindc2) ,intent(in)  :: m1,m2
  real(kindr2) :: hh
!
  include 'avh_olo_b0_a.h90'
!
  mulocal = rmu
!
  app = abs(pp)
!
  am1 = real(r1)
  hh  = aimag(r1)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = cmplx( am1 ,-hh ,kind=kindc2 )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = real(r2)
  hh  = aimag(r2)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = cmplx( am2 ,-hh ,kind=kindc2 )
  endif
  am2 = abs(am2) + abs(hh)
!
  include 'avh_olo_b0_b.h90'
!
  end subroutine
!
  subroutine b0cc( rslt ,pp,m1,m2 )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub0
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: pp,m1,m2
  real(kindr2) :: hh
!
  include 'avh_olo_b0_a.h90'
!
  mulocal = muscale
!
  app = real(ss)
  if (aimag(ss).ne.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = cmplx( app ,kind=kindc2 )
  endif
  app = abs(app)
!
  am1 = real(r1)
  hh  = aimag(r1)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = cmplx( am1 ,-hh ,kind=kindc2 )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = real(r2)
  hh  = aimag(r2)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = cmplx( am2 ,-hh ,kind=kindc2 )
  endif
  am2 = abs(am2) + abs(hh)
!
  include 'avh_olo_b0_b.h90'
!
  end subroutine
!
  subroutine b0ccr( rslt ,pp,m1,m2 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub0
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: pp,m1,m2
  real(kindr2)    ,intent(in)  :: rmu
  real(kindr2) :: hh
!
  include 'avh_olo_b0_a.h90'
!
  mulocal = rmu
!
  app = real(ss)
  if (aimag(ss).ne.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = cmplx( app ,kind=kindc2 )
  endif
  app = abs(app)
!
  am1 = real(r1)
  hh  = aimag(r1)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = cmplx( am1 ,-hh ,kind=kindc2 )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = real(r2)
  hh  = aimag(r2)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = cmplx( am2 ,-hh ,kind=kindc2 )
  endif
  am2 = abs(am2) + abs(hh)
!
  include 'avh_olo_b0_b.h90'
!
  end subroutine


!*******************************************************************
! Return the Papparino-Veltman functions b11,b00,b1,b0 , for
!
!      C   /      d^(Dim)q
!   ------ | -------------------- = b0
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
!      C   /    d^(Dim)q q^mu
!   ------ | -------------------- = p^mu b1
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
!      C   /  d^(Dim)q q^mu q^nu
!   ------ | -------------------- = g^{mu,nu} b00 + p^mu p^nu b11
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************
!
  subroutine b11rr( b11,b00,b1,b0 ,pp,m1,m2 )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub11
  complex(kindc2) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2)    ,intent(in)  :: pp,m1,m2
!
  include 'avh_olo_b11_a.h90'
!
  mulocal = muscale
!
  app = abs(pp)
  am1 = abs(m1)
  am2 = abs(m2)
!
  include 'avh_olo_b11_b.h90'
!
  end subroutine
!
  subroutine b11rrr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub11
  complex(kindc2) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2)    ,intent(in)  :: pp,m1,m2,rmu
!
  include 'avh_olo_b11_a.h90'
!
  mulocal = rmu
!
  app = abs(pp)
  am1 = abs(m1)
  am2 = abs(m2)
!
  include 'avh_olo_b11_b.h90'
!
  end subroutine
!
  subroutine b11rc( b11,b00,b1,b0 ,pp,m1,m2 )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub11
  complex(kindc2) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2)    ,intent(in)  :: pp
  complex(kindc2) ,intent(in)  :: m1,m2
  real(kindr2) :: hh
!
  include 'avh_olo_b11_a.h90'
!
  mulocal = muscale
!
  app = abs(pp)
!
  am1 = real(r1)
  hh  = aimag(r1)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = cmplx( am1 ,-hh ,kind=kindc2 )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = real(r2)
  hh  = aimag(r2)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = cmplx( am2 ,-hh ,kind=kindc2 )
  endif
  am2 = abs(am2) + abs(hh)
!
  include 'avh_olo_b11_b.h90'
!
  end subroutine
!
  subroutine b11rcr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub11
  complex(kindc2) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2)    ,intent(in)  :: pp ,rmu
  complex(kindc2) ,intent(in)  :: m1,m2
  real(kindr2) :: hh
!
  include 'avh_olo_b11_a.h90'
!
  mulocal = rmu
!
  app = abs(pp)
!
  am1 = real(r1)
  hh  = aimag(r1)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = cmplx( am1 ,-hh ,kind=kindc2 )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = real(r2)
  hh  = aimag(r2)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = cmplx( am2 ,-hh ,kind=kindc2 )
  endif
  am2 = abs(am2) + abs(hh)
!
  include 'avh_olo_b11_b.h90'
!
  end subroutine
!
  subroutine b11cc( b11,b00,b1,b0 ,pp,m1,m2 )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub11
  complex(kindc2) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindc2) ,intent(in)  :: pp,m1,m2
  real(kindr2) :: hh
!
  include 'avh_olo_b11_a.h90'
!
  mulocal = muscale
!
  app = real(ss)
  if (aimag(ss).ne.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = cmplx( app ,kind=kindc2 )
  endif
  app = abs(app)
!
  am1 = real(r1)
  hh  = aimag(r1)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = cmplx( am1 ,-hh ,kind=kindc2 )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = real(r2)
  hh  = aimag(r2)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = cmplx( am2 ,-hh ,kind=kindc2 )
  endif
  am2 = abs(am2) + abs(hh)
!
  include 'avh_olo_b11_b.h90'
!
  end subroutine
!
  subroutine b11ccr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_bub ,only: bub11
  complex(kindc2) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindc2) ,intent(in)  :: pp,m1,m2
  real(kindr2)    ,intent(in)  :: rmu
  real(kindr2) :: hh
!
  include 'avh_olo_b11_a.h90'
!
  mulocal = rmu
!
  app = real(ss)
  if (aimag(ss).ne.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = cmplx( app ,kind=kindc2 )
  endif
  app = abs(app)
!
  am1 = real(r1)
  hh  = aimag(r1)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = cmplx( am1 ,-hh ,kind=kindc2 )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = real(r2)
  hh  = aimag(r2)
  if (hh.gt.R0P0) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = cmplx( am2 ,-hh ,kind=kindc2 )
  endif
  am2 = abs(am2) + abs(hh)
!
  include 'avh_olo_b11_b.h90'
!
  end subroutine


!*******************************************************************
! calculates
!               C   /               d^(Dim)q
!            ------ | ---------------------------------------
!            i*pi^2 / [q^2-m1] [(q+k1)^2-m2] [(q+k1+k2)^2-m3]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps)
!             * GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
!
! input:  p1=k1^2, p2=k2^2, p3=(k1+k2)^2,  m1,m2,m3=squared masses
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************
!
  subroutine c0rr( rslt ,p1,p2,p3 ,m1,m2,m3 )
!*******************************************************************
!*******************************************************************
  use avh_olo_tri
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2) ,intent(in)  :: p1,p2,p3 ,m1,m2,m3
  real(kindr2) :: pp(3),mm(3)
!
  include 'avh_olo_c0_a.h90'
!
  mulocal = muscale
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_c0_b.h90'
!
  end subroutine
!
  subroutine c0rrr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_tri
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2) ,intent(in)  :: p1,p2,p3 ,m1,m2,m3 ,rmu
  real(kindr2) :: pp(3),mm(3)
!
  include 'avh_olo_c0_a.h90'
!
  mulocal = rmu
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_c0_b.h90'
!
  end subroutine
!
  subroutine c0rc( rslt ,p1,p2,p3 ,m1,m2,m3 )
!*******************************************************************
!*******************************************************************
  use avh_olo_tri
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: p1,p2,p3
  complex(kindc2) ,intent(in)  :: m1,m2,m3
  real(kindr2)    :: pp(3)
  complex(kindc2) :: mm(3)
  real(kindr2) :: hh
!
  include 'avh_olo_c0_a.h90'
!
  mulocal = muscale
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = real(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = cmplx( am(ii) ,-hh ,kind=kindc2 )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_c0_b.h90'
!
  end subroutine
!
  subroutine c0rcr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_tri
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: p1,p2,p3 ,rmu
  complex(kindc2) ,intent(in)  :: m1,m2,m3
  real(kindr2)    :: pp(3)
  complex(kindc2) :: mm(3)
  real(kindr2) :: hh
!
  include 'avh_olo_c0_a.h90'
!
  mulocal = rmu
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = real(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = cmplx( am(ii) ,-hh ,kind=kindc2 )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_c0_b.h90'
!
  end subroutine
!
  subroutine c0cc( rslt ,p1,p2,p3 ,m1,m2,m3 )
!*******************************************************************
!*******************************************************************
  use avh_olo_tri
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: p1,p2,p3 ,m1,m2,m3
  complex(kindc2) :: pp(3),mm(3)
  real(kindr2) :: hh
!
  include 'avh_olo_c0_a.h90'
!
  mulocal = muscale
!
  do ii=1,3
    ap(ii) = real(pp(ii))
    if (aimag(pp(ii)).ne.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = cmplx( ap(ii) ,kind=kindc2 )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = real(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = cmplx( am(ii) ,-hh ,kind=kindc2 )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_c0_b.h90'
!
  end subroutine
!
  subroutine c0ccr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_tri
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: p1,p2,p3 ,m1,m2,m3
  real(kindr2)    ,intent(in)  :: rmu
  complex(kindc2) :: pp(3),mm(3)
  real(kindr2) :: hh
!
  include 'avh_olo_c0_a.h90'
!
  mulocal = rmu
!
  do ii=1,3
    ap(ii) = real(pp(ii))
    if (aimag(pp(ii)).ne.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = cmplx( ap(ii) ,kind=kindc2 )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = real(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = cmplx( am(ii) ,-hh ,kind=kindc2 )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_c0_b.h90'
!
  end subroutine


!*******************************************************************
! calculates
!
!    C   /                      d^(Dim)q
! ------ | --------------------------------------------------------
! i*pi^2 / [q^2-m1][(q+k1)^2-m2][(q+k1+k2)^2-m3][(q+k1+k2+k3)^2-m4]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps)
!             * GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
!
! input:  p1=k1^2, p2=k2^2, p3=k3^2, p4=(k1+k2+k3)^2, 
!         p12=(k1+k2)^2, p23=(k2+k3)^2, 
!         m1,m2,m3,m4=squared masses
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  avh_olo_onshell  to find out how this
! routines decides when to return IR-divergent cases.
!*******************************************************************
!
  subroutine d0rr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
!*******************************************************************
!*******************************************************************
  use avh_olo_box
  use avh_olo_boxc
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4
  real(kindr2) :: pp(6),mm(4)
!
  include 'avh_olo_d0_a.h90'
!
  mulocal = muscale
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_d0_b.h90'
!
  end subroutine
!
  subroutine d0rrr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_box
  use avh_olo_boxc
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu
  real(kindr2) :: pp(6),mm(4)
!
  include 'avh_olo_d0_a.h90'
!
  mulocal = rmu
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_d0_b.h90'
!
  end subroutine
!
  subroutine d0rc( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
!*******************************************************************
!*******************************************************************
  use avh_olo_box
  use avh_olo_boxc
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindc2) ,intent(in)  :: m1,m2,m3,m4
  real(kindr2)    :: pp(6)
  complex(kindc2) :: mm(4)
  real(kindr2) :: hh
!
  include 'avh_olo_d0_a.h90'
!
  mulocal = muscale
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = real(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = cmplx( am(ii) ,-hh ,kind=kindc2 )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_d0_b.h90'
!
  end subroutine
!
  subroutine d0rcr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_box
  use avh_olo_boxc
  complex(kindc2) ,intent(out) :: rslt(0:2)
  real(kindr2)    ,intent(in)  :: p1,p2,p3,p4,p12,p23 ,rmu
  complex(kindc2) ,intent(in)  :: m1,m2,m3,m4
  real(kindr2)    :: pp(6)
  complex(kindc2) :: mm(4)
  real(kindr2) :: hh
!
  include 'avh_olo_d0_a.h90'
!
  mulocal = rmu
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = real(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = cmplx( am(ii) ,-hh ,kind=kindc2 )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_d0_b.h90'
!
  end subroutine
!
  subroutine d0cc( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
!*******************************************************************
!*******************************************************************
  use avh_olo_box
  use avh_olo_boxc
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4
  complex(kindc2) :: pp(6),mm(4)
  real(kindr2) :: hh
!
  include 'avh_olo_d0_a.h90'
!
  mulocal = muscale
!
  do ii=1,6
    ap(ii) = real(pp(ii))
    if (aimag(pp(ii)).ne.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = cmplx( ap(ii) ,R0P0 ,kind=kindc2 )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = real(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = cmplx( am(ii) ,-hh ,kind=kindc2 )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_d0_b.h90'
!
  end subroutine
!
  subroutine d0ccr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
!*******************************************************************
!*******************************************************************
  use avh_olo_box
  use avh_olo_boxc
  complex(kindc2) ,intent(out) :: rslt(0:2)
  complex(kindc2) ,intent(in)  :: p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4
  real(kindr2)    ,intent(in)  :: rmu
  complex(kindc2) :: pp(6),mm(4)
  real(kindr2) :: hh
!
  include 'avh_olo_d0_a.h90'
!
  mulocal = rmu
!
  do ii=1,6
    ap(ii) = real(pp(ii))
    if (aimag(pp(ii)).ne.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = cmplx( ap(ii) ,R0P0 ,kind=kindc2 )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = real(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.R0P0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = cmplx( am(ii) ,-hh ,kind=kindc2 )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  include 'avh_olo_d0_b.h90'
!
  end subroutine

end module
