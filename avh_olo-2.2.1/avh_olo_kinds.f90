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


module avh_olo_kinds
  use avh_olo_xkind
!
  implicit none
  private
  public :: kindr2,kindc2 &
           ,R0P0,R1P0,R5M1,TWOPI,SQRT2,C0P0,C1P0,CiP0
!
  integer ,parameter :: kindr2 = olo_xkind
  integer ,parameter :: kindc2 = kindr2
!
  real(kindr2) ,parameter :: R0P0=0._kindr2
  real(kindr2) ,parameter :: R1P0=1._kindr2
  real(kindr2) ,parameter :: R5M1=0.5_kindr2
!                                  1 2345678901234567890123456789012 
  real(kindr2) ,parameter :: TWOPI=6.2831853071795864769252867665590_kindr2
  real(kindr2) ,parameter :: SQRT2=1.4142135623730950488016887242097_kindr2
  complex(kindc2) ,parameter :: C0P0 = (0._kindr2,0._kindr2)
  complex(kindc2) ,parameter :: C1P0 = (1._kindr2,0._kindr2)
  complex(kindc2) ,parameter :: CiP0 = (0._kindr2,1._kindr2)
!  
end module
