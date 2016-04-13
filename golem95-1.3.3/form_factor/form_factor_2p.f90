!
!****h* src/form_factor/form_factor_2p
! NAME
!
!  Module form_factor_2p
!
! USAGE
!
!  use form_factor_2p
!
! DESCRIPTION
!
!  This module contains the different form factors for two point amplitudes.
!
! OUTPUT
!
!  It exports four functions:
!  * a20 -- a function to compute A^{2,0}
!  * a21 -- a function to compute A^{2,1}
!  * a22 -- a function to compute A^{2,2}
!  * b22 -- a function to compute B^{2,2}
!
!  Note that a2xx and b2xx are generic functions which can be called either with a
!  set of integers or with an integer whose bits represents the set
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * generic_function_2p (src/integrals/two_point/generic_function_2p.f90)
!  * matrice_s (src/kinematic/matrice_s.f90)
!  * array (src/module/array.f90)
!  * form_factor_type (src/module/form_factor_type.f90)
!  * constante (src/module/constante.f90)
!
!*****
module form_factor_2p
  !
  use precision_golem
  use matrice_s
  use form_factor_type
  use generic_function_2p
  use array
  use constante, only: czero
  !
  implicit none
  !
  private
  !
  interface a20
    !
    module procedure a20_b, a20_s
    !
  end interface
  !
  interface a21
    !
    module procedure a21_b, a21_s
    !
  end interface
  !
  interface a22
    !
    module procedure a22_b, a22_s
    !
  end interface
  !
  interface b22
    !
    module procedure b22_b, b22_s
    !
  end interface
  !
  public :: a20,a21,a22,b22
  !
  contains
    !
    !****f* src/form_factor/form_factor_2p/a20_b
    ! NAME
    !
    !  Function a20_b
    !
    ! USAGE
    !
    !  type(form_factor) = a20_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{2,0}.
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
    function a20_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: a20_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp
      !
      b_pro = pminus(b_ref,b_pin)
      !
      temp(:) = czero
      temp(2:3) = f2p(s_mat_p,b_pro)
      a20_b = temp
      !
    end function a20_b
    !
    !****f* src/form_factor/form_factor_2p/a20_s
    ! NAME
    !
    !  Function a20_s
    !
    ! USAGE
    !
    !  type(form_factor) = a20_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{2,0}.
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
    function a20_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a20_s
      !
      a20_s = a20_b(packb(set))
      !
    end function a20_s
    !
    !****f* src/form_factor/form_factor_2p/a21_b
    ! NAME
    !
    !  Function a21_b
    !
    ! USAGE
    !
    !  type(form_factor) = a21_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{2,1}(l_1).
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
    function a21_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: a21_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp
      !
      b_pro = pminus(b_ref,b_pin)
      !
      temp(:) = czero
      temp(2:3) = -f2p(s_mat_p,b_pro,l1)
      a21_b = temp
      !
    end function a21_b
    !
    !****f* src/form_factor/form_factor_2p/a21_s
    ! NAME
    !
    !  Function a21_s
    !
    ! USAGE
    !
    !  type(form_factor) = a21_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{2,1}(l_1).
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
    function a21_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a21_s
      !
      a21_s = a21_b(l1,packb(set))
      !
    end function a21_s
    !
    !****f* src/form_factor/form_factor_2p/a22_b
    ! NAME
    !
    !  Function a22_b
    !
    ! USAGE
    !
    !  type(form_factor) = a22_b(l1,l2,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{2,2}(l1,l2).
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
    function a22_b(l1,l2,b_pin)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in) :: b_pin
      type(form_factor) :: a22_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp
      !
      b_pro = pminus(b_ref,b_pin)
      !
      temp(:) = czero
      temp(2:3) = f2p(s_mat_p,b_pro,l1,l2)
      a22_b = temp
      !
    end function a22_b
    !
    !****f* src/form_factor/form_factor_2p/a22_s
    ! NAME
    !
    !  Function a22_s
    !
    ! USAGE
    !
    !  type(form_factor) = a22_s(l1,l2,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{2,2}(l1,l2).
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
    function a22_s(l1,l2,set)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a22_s
      !
      a22_s = a22_b(l1,l2,packb(set))
      !
    end function a22_s
    !
    !****f* src/form_factor/form_factor_2p/b22_b
    ! NAME
    !
    !  Function b22_b
    !
    ! USAGE
    !
    !  type(form_factor) = b22_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{2,2}.
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
    function b22_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: b22_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp
      !
      b_pro = pminus(b_ref,b_pin)
      !
      temp(:) = czero
      temp(2:3) = -f2p_np2(s_mat_p,b_pro)/2._ki
      b22_b = temp
      !
    end function b22_b
    !
    !****f* src/form_factor/form_factor_2p/b22_s
    ! NAME
    !
    !  Function b22_s
    !
    ! USAGE
    !
    !  type(form_factor) = b22_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{2,2}.
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
    !
    !*****
    function b22_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b22_s
      !
      b22_s = b22_b(packb(set))
      !
    end function b22_s
    !
end module form_factor_2p
