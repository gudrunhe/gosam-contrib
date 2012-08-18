!
!****h* src/form_factor/form_factor_1p
! NAME
!
!  Module form_factor_1p
!
! USAGE
!
!  use form_factor_1p
!
! DESCRIPTION
!
!  This module contains the form factor for tadpoles.
!
! OUTPUT
!
!  It exports the functions:
!  * a10 -- a function to compute A^{1,0}
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * generic_function_1p (src/integrals/one_point/generic_function_1p.f90)
!  * matrice_s (src/kinematic/matrice_s.f90)
!  * array (src/module/array.f90)
!  * form_factor_type (src/module/form_factor_type.f90)
!  * constante (src/module/constante.f90)
!
!*****
!
module form_factor_1p
  !
  use precision_golem
  use generic_function_1p
  use matrice_s
  use array
  use form_factor_type
  use constante, only: czero
  implicit none
  !
  private
  !
  interface a10
    !
    module procedure a10_b, a10_s
    !
  end interface
  !
  !
  public :: a10
  !
  contains
    !
    !****f* src/form_factor/form_factor_1p/a10_b
    ! NAME
    !
    !  Function a10_b
    !
    ! USAGE
    !
    !  type(form_factor) = a10_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{1,0}.
    ! 
    ! INPUTS
    !
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a10_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: a10_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp
      !
      b_pro = pminus(b_ref,b_pin)
      !
      temp(:) = czero
      temp(2:3) = f1p(s_mat_p,b_pro)
      a10_b = temp
      !
    end function a10_b
    !
    !****f* src/form_factor/form_factor_1p/a10_s
    ! NAME
    !
    !  Function a10_s
    !
    ! USAGE
    !
    !  type(form_factor) = a10_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{1,0}.
    ! 
    ! INPUTS
    !
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a10_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a10_s
      !
      a10_s = a10_b(packb(set))
      !
    end function a10_s
    !
    !
end module form_factor_1p
