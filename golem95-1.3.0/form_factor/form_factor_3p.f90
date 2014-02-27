!
!****h* src/form_factor/form_factor_3p
! NAME
!
!  Module form_factor_3p
!
! USAGE
!
!  use form_factor_3p
!
! DESCRIPTION
!
!  This module contains the different form factors for three point amplitudes.
!
! OUTPUT
!
!  It exports six functions:
!  * a30 -- a function to compute A^{3,0}
!  * a31 -- a function to compute A^{3,1}
!  * a32 -- a function to compute A^{3,2}
!  * a33 -- a function to compute A^{3,3}
!  * b32 -- a function to compute B^{3,2}
!  * b33 -- a function to compute B^{3,3}
!
!  Note that a3xx and b3xx are generic functions which can be called either with a
!  set of integers or with an integer whose bits represents the set
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * generic_function_3p (src/integrals/three_point/generic_function_3p.f90)
!  * matrice_s (src/kinematic/matrice_s.f90)
!  * array (src/module/array.f90)
!  * form_factor_type (src/module/form_factor_type.f90)
!  * constante (src/module/constante.f90)
!
!*****
module form_factor_3p
  !
  use precision_golem
  use generic_function_3p
  use matrice_s
  use array
  use form_factor_type
  use constante, only: czero
  !
  implicit none
  !
  private 
  !
  interface a30
    !
    module procedure a30_b, a30_s
    !
  end interface
  !
  interface a31
    !
    module procedure a31_b, a31_s
    !
  end interface
  !
  interface a32
    !
    module procedure a32_b, a32_s
    !
  end interface
  !
  interface a33
    !
    module procedure a33_b, a33_s
    !
  end interface
  !
  interface b32
    !
    module procedure b32_b, b32_s
    !
  end interface
  !
  interface b33
    !
    module procedure b33_b, b33_s
    !
  end interface
  !
  public :: a30,a31,a32,a33,b32,b33
  !
  contains
    !
    !****f* src/form_factor/form_factor_3p/a30_b
    ! NAME
    !
    !  Function a30_b
    !
    ! USAGE
    !
    !  type(form_factor) = a30_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{3,0}.
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
    !  and s_mat_p
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
    function a30_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: a30_b
      !
      integer :: b_pro
      !
      b_pro = pminus(b_ref,b_pin)
      !
      a30_b = f3p(s_mat_p,b_pro)
      !
    end function a30_b
    !
    !****f* src/form_factor/form_factor_3p/a30_s
    ! NAME
    !
    !  Function a30_s
    !
    ! USAGE
    !
    !  type(form_factor) = a30_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{3,0}.
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
    function a30_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a30_s
      !
      a30_s = a30_b(packb(set))
      !
    end function a30_s
    !
    !****f* src/form_factor/form_factor_3p/a31_b
    ! NAME
    !
    !  Function a31_b
    !
    ! USAGE
    !
    !  type(form_factor) = a31_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{3,1}(l_1).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
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
    function a31_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: a31_b
      !
      integer :: b_pro
      !
      b_pro = pminus(b_ref,b_pin)
      !
      a31_b = -f3p(s_mat_p,b_pro,l1)
      !
    end function a31_b
    !
    !****f* src/form_factor/form_factor_3p/a31_s
    ! NAME
    !
    !  Function a31_s
    !
    ! USAGE
    !
    !  type(form_factor) = a31_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{3,1}(l_1).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
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
    function a31_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a31_s
      !
      a31_s = a31_b(l1,packb(set))
      !
    end function a31_s
    !
    !****f* src/form_factor/form_factor_3p/a32_b
    ! NAME
    !
    !  Function a32_b
    !
    ! USAGE
    !
    !  type(form_factor) = a32_b(l1,l2,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{3,2}(l1,l2).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
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
    function a32_b(l1,l2,b_pin)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in) :: b_pin
      type(form_factor) :: a32_b
      !
      integer :: b_pro
      !
      b_pro = pminus(b_ref,b_pin)
      !
      a32_b = f3p(s_mat_p,b_pro,l1,l2)
      !
    end function a32_b
    !
    !****f* src/form_factor/form_factor_3p/a32_s
    ! NAME
    !
    !  Function a32_s
    !
    ! USAGE
    !
    !  type(form_factor) = a32_s(l1,l2,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{3,2}(l1,l2).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
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
    function a32_s(l1,l2,set)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a32_s
      !
      a32_s = a32_b(l1,l2,packb(set))
      !
    end function a32_s
    !
    !****f* src/form_factor/form_factor_3p/a33_b
    ! NAME
    !
    !  Function a33_b
    !
    ! USAGE
    !
    !  type(form_factor) = a33_b(l1,l2,l3,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{3,3}(l1,l2,l3).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
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
    function a33_b(l1,l2,l3,b_pin)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in) :: b_pin
      type(form_factor) :: a33_b
      !
      integer :: b_pro
      !
      b_pro = pminus(b_ref,b_pin)
      !
      a33_b = -f3p(s_mat_p,b_pro,l1,l2,l3)
      !
    end function a33_b
    !
    !****f* src/form_factor/form_factor_3p/a33_s
    ! NAME
    !
    !  Function a33_s
    !
    ! USAGE
    !
    !  type(form_factor) = a33_s(l1,l2,l3,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{3,3}(l1,l2,l3).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
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
    function a33_s(l1,l2,l3,set)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a33_s
      !
      a33_s = a33_b(l1,l2,l3,packb(set))
      !
    end function a33_s
    !
    !****f* src/form_factor/form_factor_3p/b32_b
    ! NAME
    !
    !  Function b32_b
    !
    ! USAGE
    !
    !  type(form_factor) = b32_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{3,2}.
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
    !  and s_mat_p
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
    function b32_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: b32_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp
      !
      b_pro = pminus(b_ref,b_pin)
      !
      temp(:) = czero
      temp(2:3) = -f3p_np2(s_mat_p,b_pro)/2._ki
      b32_b = temp
      !
    end function b32_b
    !
    !****f* src/form_factor/form_factor_3p/b32_s
    ! NAME
    !
    !  Function b32_s
    !
    ! USAGE
    !
    !  type(form_factor) = b32_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{3,2}.
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
    function b32_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b32_s
      !
      b32_s = b32_b(packb(set))
      !
    end function b32_s
    !
    !****f* src/form_factor/form_factor_3p/b33_b
    ! NAME
    !
    !  Function b33_b
    !
    ! USAGE
    !
    !  type(form_factor) = b33_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{3,3}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
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
    function b33_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: b33_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp
      !
      b_pro = pminus(b_ref,b_pin)
      !
      temp(:) = czero
      temp(2:3) = f3p_np2(s_mat_p,b_pro,l1)/2._ki
      b33_b = temp
      !
    end function b33_b
    !
    !****f* src/form_factor/form_factor_3p/b33_s
    ! NAME
    !
    !  Function b33_s
    !
    ! USAGE
    !
    !  type(form_factor) = b33_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{3,3}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
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
    function b33_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b33_s
      !
      b33_s = b33_b(l1,packb(set))
      !
    end function b33_s
    !
end module form_factor_3p
