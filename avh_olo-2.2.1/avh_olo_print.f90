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


module avh_olo_print
  use avh_olo_kinds
  implicit none
  private
  public :: myprint,init_print
!
  integer ,parameter :: noverh=10 !maximally 6 decimals for exponent
  integer :: ndigits=19
  integer :: nefrmt=19+noverh
!
  interface myprint
    module procedure printr,printc,printi
  end interface
!
contains
!
  subroutine init_print( ndig )
  integer ,intent(in) :: ndig
  ndigits = ndig+ndig/4+1
  nefrmt  = ndigits+noverh
  end subroutine
!
  function printc( zz ) result(rslt)
  complex(kindc2) ,intent(in) :: zz
  character(nefrmt*2+3) :: rslt
  rslt = '('//trim(printr(real(zz))) &
       //','//trim(printr(aimag(zz)           )) &
       //')'
  rslt = adjustl(rslt)
  end function
!
  function printr( xx ) result(rslt)
  real(kindr2) ,intent(in) :: xx
  character(nefrmt  ) :: rslt
  character(nefrmt+1) :: cc
  character(10) :: aa,bb
  write(aa,'(i10)') nefrmt+1  ;aa=adjustl(aa)
  write(bb,'(i10)') ndigits   ;bb=adjustl(bb)
  aa = '(e'//trim(aa)//'.'//trim(bb)//')'
  write(cc,aa) xx  ;cc=adjustl(cc)
  if (cc(1:2).eq.'-0') then ;rslt = '-'//cc(3:ndigits*2)
  else                      ;rslt = ' '//cc(2:ndigits*2)
  endif
  end function
!
  function printi( ii ) result(rslt)
  integer ,intent(in) :: ii
  character(ndigits) :: rslt
  character(ndigits) :: cc
  character(10) :: aa
  write(aa,'(i10)') ndigits ;aa=adjustl(aa)
  aa = '(i'//trim(aa)//')'
  write(cc,aa) ii ;cc=adjustl(cc)
  if (cc(1:1).ne.'-') then ;rslt=' '//cc
  else                     ;rslt=cc 
  endif
  end function
!
end module
