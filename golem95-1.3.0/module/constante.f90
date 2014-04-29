!****h* src/module/constante
! NAME
!
!  Module constante
!
! USAGE
!
!  use constante
!
! DESCRIPTION
!
!  This module is used to get the values of different constants
!
! OUTPUT
!
!  This module exports 13 parameters
!  * pi -- a real (type ki), pi
!  * gammae -- a real (type ki), the Euler constant
!  * i_ -- a complex (type ki), the square root of -1
!  * pi3 -- a real (type ki), pi**2/3
!  * pi6 -- a real (type ki), pi**2/6
!  * pi12 -- a real (type ki), pi**2/12
!  * un -- a real (type ki), 1
!  * zero -- a real (type ki), 0
!  * czero -- a complex (type ki), (0,0)
!  * cun -- a complex (type ki), (1,0)
!  * b_null -- an integer 0
!  * s_null -- an array with shape 0
!  * nullarray -- the same as s_null !!!!!!
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!
!*****
module constante
  !
  use precision_golem
  implicit none
  !
  private :: ki
  !
  real(ki), parameter :: pi = 3.1415926535897932384626433832795028841&
  &971693993751_ki
  real(ki), parameter :: gammae = 0.577215664901532860606512090082402&
  &4310421593359399_ki
  real(ki), parameter :: pi3 = pi**2/3.0_ki
  real(ki), parameter :: pi6 = pi**2/6.0_ki
  real(ki), parameter :: pi12 = pi**2/12.0_ki
  real(ki), parameter :: un = 1._ki
  real(ki), parameter :: zero = 0._ki
  complex(ki), parameter :: i_ = (0._ki,1._ki)
  complex(ki), parameter :: czero = (0.0_ki,0.0_ki)
  complex(ki), parameter :: cun = (1._ki,0._ki)
  integer, parameter :: b_null = 0
  integer, dimension(0), parameter :: s_null = 0
  integer, dimension(0), parameter :: nullarray = 0
  !
end module constante
