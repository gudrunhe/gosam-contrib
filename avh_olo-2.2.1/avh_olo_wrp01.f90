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


      subroutine avh_olo_mu_set(mu)
      use avh_olo
      implicit none
      real(olo_kind) ,intent(in) :: mu
      call olo_scale( mu )
      end subroutine

      subroutine avh_olo_onshell(thrs)
      use avh_olo
      implicit none
      real(olo_kind) ,intent(in) :: thrs
      call olo_onshell( thrs )
      end subroutine

      subroutine avh_olo_unit( unit_in )
      use avh_olo
      implicit none
      integer ,intent(in) :: unit_in
      call olo_unit( unit_in ,'all' )
      end subroutine

      subroutine avh_olo_printall( unit_in )
      use avh_olo
      implicit none
      integer ,intent(in) :: unit_in
      call olo_unit( unit_in ,'printall' )
      end subroutine

      subroutine avh_olo_a0c( rslt ,mm )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: rslt(0:2)
      complex(olo_kind) ,intent(in)  :: mm
      call olo_a0( rslt ,mm )
      end subroutine

      subroutine avh_olo_b0c( rslt ,pp,m1,m2 )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: rslt(0:2)
      complex(olo_kind) ,intent(in)  :: pp,m1,m2
      call olo_b0( rslt ,pp,m1,m2 )
      end subroutine

      subroutine avh_olo_b11c( b11,b00,b1,b0 ,pp,m1,m2 )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
      complex(olo_kind) ,intent(in)  :: pp,m1,m2
      call olo_b11( b11,b00,b1,b0 ,pp,m1,m2 )
      end subroutine

      subroutine avh_olo_c0c( rslt ,p1,p2,p3 ,m1,m2,m3 )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: rslt(0:2)
      complex(olo_kind) ,intent(in)  :: p1,p2,p3 ,m1,m2,m3
      call olo_c0( rslt ,p1,p2,p3 ,m1,m2,m3 )
      end subroutine

      subroutine avh_olo_d0c( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: rslt(0:2)
      complex(olo_kind) ,intent(in)  :: p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4
      call olo_d0( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
      end subroutine

      subroutine avh_olo_a0m( rslt ,mm )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: rslt(0:2)
      real(olo_kind)    ,intent(in)  :: mm
      call olo_a0( rslt ,mm )
      end subroutine

      subroutine avh_olo_b0m( rslt ,pp,m1,m2 )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: rslt(0:2)
      real(olo_kind)    ,intent(in)  :: pp,m1,m2
      call olo_b0( rslt ,pp,m1,m2 )
      end subroutine

      subroutine avh_olo_b11m( b11,b00,b1,b0 ,pp,m1,m2 )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
      real(olo_kind)    ,intent(in)  :: pp,m1,m2
      call olo_b11( b11,b00,b1,b0 ,pp,m1,m2 )
      end subroutine

      subroutine avh_olo_c0m( rslt ,p1,p2,p3 ,m1,m2,m3 )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: rslt(0:2)
      real(olo_kind)    ,intent(in)  :: p1,p2,p3 ,m1,m2,m3
      call olo_c0( rslt ,p1,p2,p3 ,m1,m2,m3 )
      end subroutine

      subroutine avh_olo_d0m( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
      use avh_olo
      implicit none
      complex(olo_kind) ,intent(out) :: rslt(0:2)
      real(olo_kind)    ,intent(in)  :: p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4
      call olo_d0( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
      end subroutine
