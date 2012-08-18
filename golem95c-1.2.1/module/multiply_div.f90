! 
!****h* src/module/multiply_div
! NAME
!
!  Module multiply_div
!
! USAGE
!
!  use multiply_div
!
! DESCRIPTION
!
!  This module contains the function mult_div, This function computes 
!  numericaly (1+alpha*epsilon)*(A/epsilon+B). The type of the output array is
!  identical to the type of the input array.
!
! OUTPUT
!
!  This module exports the function mult_div
!
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!
!*****
!
module multiply_div
  !
  use precision_golem
  !
  implicit none
  !
  private 
  !
  interface mult_div
     !
     module procedure mult_div_r
     module procedure mult_div_c
     !
  end interface
  !  
  public :: mult_div
  !
  contains
  !
  !****f* src/module/multiply_div/mult_div_r
  ! NAME
  !
  !  Function mult_div_r
  !
  ! USAGE
  !
  !  real_dim4 = mult_div_r(alpha,array)
  !
  ! DESCRIPTION
  !
  !  This function computes numericaly (1+alpha*epsilon)*(A/epsilon+B)
  !  with A = a1 + i*a2 and B = b1 + i*b2. The returned result is put
  !  into an array t (rank 1, shape 4) where t(1) = a1, t(2) = a2,
  !  t(3) = b1+alpha*a1, t(4) = b2+alpha*a2. 
  !
  ! INPUTS
  !
  !  * alpha -- a real (type ki)
  !  * array -- a real (type ki) array of rank 1, shape 4
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  a real (type ki) array of rank 1 and shape 4
  !
  ! NOTES
  !
  ! The return value of this function is a real array of shape 4,
  ! contrary to the complex array returned by mult_div_c.
  !
  ! EXAMPLE
  !
  !  resu = multipy_div_r(alpha,array)
  !  resu(1) = array(1)
  !  resu(2) = array(2)
  !  resu(3) = array(3) + alpha*array(1)
  !  resu(4) = array(4) + alpha*array(2)
  !
  !*****
  function mult_div_r(alpha,array)
    !
    real(ki), intent(in) :: alpha
    real(ki), intent(in), dimension(4) :: array
    real(ki), dimension(4) :: mult_div_r
    !
    mult_div_r = array
    mult_div_r(3) = mult_div_r(3) + alpha*array(1)
    mult_div_r(4) = mult_div_r(4) + alpha*array(2)
    !
  end function mult_div_r
  !
  !
  !****f* src/module/multiply_div/mult_div_c
  ! NAME
  !
  !  Function mult_div_c
  !
  ! USAGE
  !
  !  cmplx_dim2 = mult_div_c(alpha,array)
  !
  ! DESCRIPTION
  !
  !  This function computes numericaly (1+alpha*epsilon)*(A/epsilon+B)
  !  with A and B complex. The returned result is put
  !  into an complex array t (rank 1, shape 2) where t(1) = A,
  !  t(2) = B + alpha*A. 
  !
  ! INPUTS
  !
  !  * alpha -- a real (type ki)
  !  * array -- a complex (type ki) array of rank 1, shape 2
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  a complex (type ki) array of rank 1 and shape 2
  !
  ! NOTES
  !
  ! The return value of this function is a complex array of shape 2,
  ! contrary to the real array returned by mult_div_r.
  !
  ! EXAMPLE
  !
  !  resu = multipy_div_c(alpha,array)
  !  resu(1) = array(1)
  !  resu(2) = array(2) + alpha*array(1)
  !
  !*****
  function mult_div_c(alpha,array)
    !
    real(ki), intent(in) :: alpha
    complex(ki), intent(in), dimension(2) :: array
    complex(ki), dimension(2) :: mult_div_c
    !
    mult_div_c = array
    mult_div_c(2) = mult_div_c(2) + alpha*array(1)
    !
  end function mult_div_c
  !
end module multiply_div
