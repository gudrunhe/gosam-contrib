! NAME
! SYNOPSIS
!****h* src/module/form_factor_type
! NAME
!
!  Module form_factor_type
!
! USAGE
!
!  use form_factor_type
!
! DESCRIPTION
!
!  This module contains two type definitions : the form factors and
!  epsilon type. This module overloads the *, /, +, -, = and ** operators
!
! OUTPUT
!
!  This module exports two types:
!  * form_factor -- define the type of the fom factors
!  * epsilon_type -- define the type for object having an epsilon expansion
!
!  five operators:
!  * * -- overload of the multiplication operator for form_factor and epsilon_type object
!  * / -- overload of the division operator for form_factor and epsilon_type object
!  * + -- overload of the addition operator for form_factor and epsilon_type object
!  * - -- overload of the subtraction operator for form_factor and epsilon_type object
!  * = -- overload of the assignment operator for form_factor and epsilon_type object
!  * ** -- overload of the power operator for form_factor and epsilon_type object
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!  * constante (src/module/constante.f90)
!
!
! AUTHOR
!   Thomas Reiter
!
! CREATION DATE
!   Oct 19, 2007
!
!*****
module form_factor_type
   !
   use precision_golem, only: ki
   use constante
   !
   implicit none
   !
   private :: ki
   !
   !****t* src/module/form_factor_type/form_factor
   ! NAME
   !   form_factor -- represents the result of a form factor
   !
   ! SYNOPSIS
   !   type form_factor
   !
   ! SOURCE
   type form_factor
      complex(ki) :: a
      complex(ki) :: b
      complex(ki) :: c
   end type form_factor
   !
   ! NOTES
   !   * a is the coefficient of the 1/epsilon^2 pole
   !   * b is the coefficient of the 1/epsilon pole
   !   * c is the coefficient of the finite term
   !****
   !
   !****t* src/module/form_factor_type/epsilon_type
   ! NAME
   !   epsilon_type -- a type that represents positive
   !     powers of epsilon
   !
   ! SYNOPSIS
   !   type epsilon_type
   !
   ! SOURCE
   type epsilon_type
      complex(ki) :: coefficient
      integer :: power
   end type epsilon_type
   !
   !****
   !
   !****t* src/module/form_factor_type/eps
   ! NAME
   !   eps -- singleton object of the epsilon-type.
   !
   ! SYNOPSIS
   !   type(epsilon_type), parameter :: eps
   !
   ! SOURCE
   type(epsilon_type), parameter :: eps = epsilon_type(1.0, 1)
   !
   ! EXAMPLE
   !   type(form_factor) :: ff1 = a20((/.../))
   !   type(form_factor) :: ff2 = eps * ff1
   !****

   !****** src/module/form_factor_type/multiplication
   ! NAME
   !   c * ff -- Multiplication of a form_factor with a scalar
   !
   ! SYNOPSIS
   interface operator(*)
   !****
      module procedure mul_complex_ff
      module procedure mul_ff_complex
      module procedure mul_real_ff
      module procedure mul_ff_real
      module procedure mul_integer_ff
      module procedure mul_ff_integer
      module procedure mul_ff_eps
      module procedure mul_eps_ff
      module procedure mul_eps_eps
      module procedure mul_eps_real
      module procedure mul_eps_complex
      module procedure mul_eps_integer
      module procedure mul_real_eps
      module procedure mul_complex_eps
      module procedure mul_integer_eps
   end interface
   !
   private :: mul_complex_ff, mul_ff_complex
   private :: mul_real_ff, mul_ff_real
   private :: mul_integer_ff, mul_ff_integer
   private :: mul_ff_eps, mul_eps_ff, mul_eps_eps
   private :: mul_eps_real, mul_real_eps
   private :: mul_eps_complex, mul_complex_eps
   private :: mul_eps_integer, mul_integer_eps
   !
   !****** src/module/form_factor_type/division
   ! NAME
   !   ff / c -- Division of a form_factor by a scalar
   !
   ! SYNOPSIS
   interface operator(/)
   !****
      module procedure div_ff_complex
      module procedure div_ff_real
      module procedure div_ff_integer
   end interface
   !
   private :: div_ff_complex, div_ff_real, div_ff_integer
   !
   !****** src/module/form_factor_type/sum
   ! NAME
   !   ff + x, + ff -- Sums involving form_factor(s)
   !
   ! SYNOPSIS
   interface operator(+)
   !****
      module procedure add_complex_ff
      module procedure add_ff_complex
      module procedure add_real_ff
      module procedure add_ff_real
      module procedure add_integer_ff
      module procedure add_ff_integer
      module procedure add_ff_ff
      module procedure plus_ff
   end interface

   private :: add_complex_ff, add_ff_complex
   private :: add_real_ff, add_ff_real
   private :: add_integer_ff, add_ff_integer
   private :: add_ff_ff, plus_ff
   
   !****** src/module/form_factor_type/subtraction
   ! NAME
   !   ff - x, - ff -- Subtractions involving form_factor(s)
   !
   ! SYNOPSIS
   interface operator(-)
   !****
      module procedure sub_complex_ff
      module procedure sub_ff_complex
      module procedure sub_real_ff
      module procedure sub_ff_real
      module procedure sub_integer_ff
      module procedure sub_ff_integer
      module procedure sub_ff_ff
      module procedure minus_ff
   end interface

   private :: sub_complex_ff, sub_ff_complex
   private :: sub_real_ff, sub_ff_real
   private :: sub_integer_ff, sub_ff_integer
   private :: sub_ff_ff, minus_ff
  

   !****** src/module/form_factor_type/assignment
   ! NAME
   !   ff = x -- Assignment to a form_factor
   !
   ! SYNOPSIS
   interface assignment(=)
   !
   ! NOTES
   !   In the assignment of a complex array the RHS must have the form
   !      (/ a, b, c /).
   !   In the assignment of a real array the RHS must have the form
   !      (/ real(a), aimag(a), real(b), aimag(b), real(c), aimag(c) /).
   !   The later form reflects the convention in the form factors of Golem90.
   !****
      module procedure assign_ff_complex
      module procedure assign_ff_real
      module procedure assign_ff_integer
      module procedure assign_ff_complex_array3
      module procedure assign_ff_real_array6
   end interface

   private :: assign_ff_complex, assign_ff_real, assign_ff_integer 
   private :: assign_ff_complex_array3, assign_ff_real_array6

   interface operator(**)
      module procedure pow_eps_int
   end interface

   private :: pow_eps_int
contains

   pure elemental function mul_complex_ff(x, ff) result(r)
      !
      complex(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r%a = x * ff%a
      r%b = x * ff%b
      r%c = x * ff%c
      !
   end function mul_complex_ff
   !
   pure elemental function mul_ff_complex(ff, x) result(r)
      !
      complex(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r%a = x * ff%a
      r%b = x * ff%b
      r%c = x * ff%c
      !
   end function mul_ff_complex
   !
   pure elemental function mul_real_ff(x, ff) result(r)
      !
      real(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0.0_ki, ki)
      r%a = z * ff%a
      r%b = z * ff%b
      r%c = z * ff%c
      !
   end function mul_real_ff
   !
   pure elemental function mul_ff_real(ff, x) result(r)
      !
      real(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0.0_ki, ki)
      r%a = z * ff%a
      r%b = z * ff%b
      r%c = z * ff%c
      !
   end function mul_ff_real
   !
   pure elemental function mul_integer_ff(x, ff) result(r)
      implicit none
      integer, intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0, ki)
      r%a = z * ff%a
      r%b = z * ff%b
      r%c = z * ff%c
      !
   end function mul_integer_ff
   !
   pure elemental function mul_ff_integer(ff, x) result(r)
      !
      integer, intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0, ki)
      r%a = z * ff%a
      r%b = z * ff%b
      r%c = z * ff%c
      !
   end function mul_ff_integer
   !
   pure elemental function pow_eps_int(eps, power) result(r)
      !
      type(epsilon_type), intent(in) :: eps
      integer, intent(in) :: power
      type(epsilon_type) :: r
      !
      r%coefficient = eps%coefficient ** power
      r%power = power * eps%power
      !
   end function pow_eps_int

   pure elemental function mul_eps_ff(x, ff) result(r)
      !
      type(epsilon_type), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      if (x%power >= 3) then
         !
         r%a = 0.0_ki
         r%b = 0.0_ki
         r%c = 0.0_ki
         !
      elseif (x%power == 2) then
         !
         r%a = 0.0_ki
         r%b = 0.0_ki 
         r%c = x%coefficient * ff%a
         !
      else
         !
         r%a = 0.0_ki
         r%b = x%coefficient * ff%a
         r%c = x%coefficient * ff%b
         !
      end if
      !
   end function mul_eps_ff
   !
   pure elemental function mul_ff_eps(ff, x) result(r)
      !
      type(epsilon_type), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      if (x%power >= 3) then
         !
         r = 0.0_ki
         !
      elseif (X%power == 2) then
         !
         r%a = 0.0_ki 
         r%b = 0.0_ki 
         r%c = x%coefficient * r%a
         !
      else
         !
         r%a = 0.0_ki
         r%b = x%coefficient * ff%a
         r%c = x%coefficient * ff%b
         !
      end if
      !
   end function mul_ff_eps
   !
   pure elemental function mul_eps_eps(eps1, eps2) result(r)
      !
      type(epsilon_type), intent(in) :: eps1, eps2
      type(epsilon_type) :: r
      !
      r%coefficient = eps1%coefficient * eps2%coefficient
      r%power = eps1%power + eps2%power
      !
   end function  mul_eps_eps
   !
   pure elemental function mul_eps_complex(eps, x) result(r)
      !
      type(epsilon_type), intent(in) :: eps
      complex(ki), intent(in) :: x
      type(epsilon_type) :: r
      !
      r%power = eps%power
      r%coefficient = x * eps%coefficient
      !
   end function  mul_eps_complex
   !
   pure elemental function mul_complex_eps(x, eps) result(r)
     !
      complex(ki), intent(in) :: x
      type(epsilon_type), intent(in) :: eps
      type(epsilon_type) :: r
      !
      r%power = eps%power
      r%coefficient = x * eps%coefficient
      !
   end function  mul_complex_eps
   !
   pure elemental function mul_eps_real(eps, x) result(r)
      !
      type(epsilon_type), intent(in) :: eps
      real(ki), intent(in) :: x
      type(epsilon_type) :: r
      !
      r%power = eps%power
      r%coefficient = x * eps%coefficient
      !
   end function  mul_eps_real
   !
   pure elemental function mul_real_eps(x, eps) result(r)
      !
      real(ki), intent(in) :: x
      type(epsilon_type), intent(in) :: eps
      type(epsilon_type) :: r
      !
      r%power = eps%power
      r%coefficient = x * eps%coefficient
      !
   end function  mul_real_eps
   !
   pure elemental function mul_eps_integer(eps, x) result(r)
      !
      type(epsilon_type), intent(in) :: eps
      integer, intent(in) :: x
      type(epsilon_type) :: r
      !
      r%power = eps%power
      r%coefficient = cmplx(x, 0, ki) * eps%coefficient
      !
   end function  mul_eps_integer
   !
   pure elemental function mul_integer_eps(x, eps) result(r)
      !
      integer, intent(in) :: x
      type(epsilon_type), intent(in) :: eps
      type(epsilon_type) :: r
      !
      r%power = eps%power
      r%coefficient = cmplx(x, 0, ki) * eps%coefficient
      !
   end function  mul_integer_eps

   pure elemental function div_ff_complex(ff, x) result(r)
      !
      complex(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r%a = ff%a / x
      r%b = ff%b / x
      r%c = ff%c / x
      !
   end function div_ff_complex
   !
   pure elemental function div_ff_real(ff, x) result(r)
      !
      real(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0.0_ki, ki)
      r%a = ff%a / z
      r%b = ff%b / z
      r%c = ff%c / z
      !
   end function div_ff_real
   !
   pure elemental function div_ff_integer(ff, x) result(r)
      !
      integer, intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0, ki)
      r%a = ff%a / z
      r%b = ff%b / z
      r%c = ff%c / z
      !
   end function div_ff_integer
   !
   pure elemental function add_complex_ff(x, ff) result(r)
      !
      complex(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c + x
      !
   end function add_complex_ff
   !
   pure elemental function add_ff_complex(ff, x) result(r)
      !
      complex(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c + x
      !
   end function add_ff_complex
   !
   pure elemental function add_real_ff(x, ff) result(r)
      !
      real(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0.0_ki, ki)
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c + z
      !
   end function add_real_ff
   !
   pure elemental function add_ff_real(ff, x) result(r)
      !
      real(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0.0_ki, ki)
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c + z
      !
   end function add_ff_real
   !
   pure elemental function add_integer_ff(x, ff) result(r)
      !
      integer, intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0, ki)
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c + z
      !
   end function add_integer_ff
   !
   pure elemental function add_ff_integer(ff, x) result(r)
      !
      integer, intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0, ki)
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c + z
      !
   end function add_ff_integer
   !
   pure elemental function add_ff_ff(x, ff) result(r)
      !
      type(form_factor), intent(in) :: ff, x
      type(form_factor) :: r
      !
      r%a = ff%a + x%a
      r%b = ff%b + x%b
      r%c = ff%c + x%c
      !
   end function add_ff_ff
   !
   pure elemental function plus_ff(ff) result(r)
      !
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r = ff
      !
   end function plus_ff
   !
   pure elemental function sub_complex_ff(x, ff) result(r)
      !
      complex(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r%a = ff%a
      r%b = ff%b
      r%c = x - ff%c
      !
   end function sub_complex_ff
   !
   pure elemental function sub_ff_complex(ff, x) result(r)
      !
      complex(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c - x
      !
   end function sub_ff_complex
   !
   pure elemental function sub_real_ff(x, ff) result(r)
      !
      real(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0.0_ki, ki)
      r%a = ff%a
      r%b = ff%b
      r%c = z - ff%c
      !
   end function sub_real_ff
   !
   pure elemental function sub_ff_real(ff, x) result(r)
      !
      real(ki), intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0.0_ki, ki)
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c - z
      !
   end function sub_ff_real
   !
   pure elemental function sub_integer_ff(x, ff) result(r)
      !
      integer, intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0, ki)
      r%a = ff%a
      r%b = ff%b
      r%c = z - ff%c
      !
   end function sub_integer_ff
   !
   pure elemental function sub_ff_integer(ff, x) result(r)
      !
      integer, intent(in) :: x
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      complex(ki) :: z
      !
      z = cmplx(x, 0, ki)
      r%a = ff%a
      r%b = ff%b
      r%c = ff%c - z
      !
   end function sub_ff_integer
   !
   pure elemental function sub_ff_ff(ff, x) result(r)
      !
      type(form_factor), intent(in) :: ff, x
      type(form_factor) :: r
      !
      r%a = ff%a - x%a
      r%b = ff%b - x%b
      r%c = ff%c - x%c
      !
   end function sub_ff_ff
   !
   pure elemental function minus_ff(ff) result(r)
      !
      type(form_factor), intent(in) :: ff
      type(form_factor) :: r
      !
      r%a = -ff%a
      r%b = -ff%b
      r%c = -ff%c
      !
   end function minus_ff
   !
   !~ pure elemental subroutine assign_ff_complex(ff, x)
   pure subroutine assign_ff_complex(ff, x)
      !
      type(form_factor), intent(out) :: ff
      complex(ki), intent(in) :: x
      !
      ff%a = (0.0_ki, 0.0_ki)
      ff%b = (0.0_ki, 0.0_ki)
      ff%c = x
      !
   end subroutine assign_ff_complex
   !
   !~ pure elemental subroutine assign_ff_real(ff, x)
   pure subroutine assign_ff_real(ff, x)
      !
      type(form_factor), intent(out) :: ff
      real(ki), intent(in) :: x
      !
      ff%a = (0.0_ki, 0.0_ki)
      ff%b = (0.0_ki, 0.0_ki)
      ff%c = cmplx(x, 0.0_ki, ki)
      !
   end subroutine assign_ff_real
   !
   !~ pure elemental subroutine assign_ff_integer(ff, x)
   pure subroutine assign_ff_integer(ff, x)
      !
      type(form_factor), intent(out) :: ff
      integer, intent(in) :: x
      !
      ff%a = (0.0_ki, 0.0_ki)
      ff%b = (0.0_ki, 0.0_ki)
      ff%c = cmplx(x, 0, ki)
      !
   end subroutine assign_ff_integer
   !
   pure subroutine assign_ff_complex_array3(ff, x)
      !
      type(form_factor), intent(out) :: ff
      complex(ki), dimension(1:3), intent(in) :: x
      !
      ff%a = x(1)
      ff%b = x(2)
      ff%c = x(3)
      !
   end subroutine assign_ff_complex_array3
   !
   pure subroutine assign_ff_real_array6(ff, x)
      !
      type(form_factor), intent(out) :: ff
      real(ki), dimension(1:6), intent(in) :: x
      !
      ff%a = x(1) + i_ * x(2)
      ff%b = x(3) + i_ * x(4)
      ff%c = x(5) + i_ * x(6)
      !
   end subroutine assign_ff_real_array6
   !
end module form_factor_type
