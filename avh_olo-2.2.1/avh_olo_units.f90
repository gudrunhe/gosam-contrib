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


module avh_olo_units
  implicit none
  integer :: eunit=6 !PROTECTED
  integer :: wunit=6 !PROTECTED
  integer :: munit=6 !PROTECTED
  integer :: punit=0 !PROTECTED ! print all

  ! protected :: eunit, wunit, munit, punit
contains
  subroutine set_unit( message ,val )
!***********************************************************************
! message is intended to be one of the following:
! 'printall', 'message' ,'warning' ,'error'
!***********************************************************************
  character(*) ,intent(in) :: message
  integer      ,intent(in) :: val
  if (.false.) then
  elseif (message(1:8).eq.'printall') then ;punit=val
  elseif (message(1:7).eq.'message' ) then ;munit=val
  elseif (message(1:7).eq.'warning' ) then ;wunit=val
  elseif (message(1:5).eq.'error'   ) then ;eunit=val
  else
    eunit=val
    wunit=val
    munit=val
    punit=0
  endif
  end subroutine
end module
