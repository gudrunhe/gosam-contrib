module     maccu
   !
   ! Synopsis:    routines for floating-point accumulation
   !              of sums with small relative error
   ! Author:      Thomas Reiter <thomasr@nikhef.nl>
   ! Date:        28 Jul. 2010
   ! Language:    Fortran 95
   ! Description: Implementation of the algorithm presented in
   !              ``AN ALGORITHM FOR FLOATING-POINT ACCUMULATION OF SUMS
   !                WITH SMALL RELATIVE ERROR''
   !              Michael Malcolm, STAN-CS-70-163, June 1970
   !
   ! Example:
   !              arr = (/ 1.2345E+20_ki, 1.0_ki, -1.2345E+20_ki/)
   !              print*, sum(arr), sorted_sum(arr)
   !              ! output:
   !              !   0.0      1.0
   !
   use precision, only: ki
   implicit none

   private

   integer, parameter, private :: min_ex_ki = minexponent(1.0_ki)
   integer, parameter, private :: max_ex_ki = maxexponent(1.0_ki)

   type     accumulator_type
     real(ki), dimension(min_ex_ki:max_ex_ki) :: a = 0.0_ki
   end type accumulator_type

   interface add_accu
      module procedure add_accu_c
      module procedure add_accu_r
   end interface

   interface reduce_accu
      module procedure reduce_accu_c
      module procedure reduce_accu_r
   end interface

   public :: accumulator_type, add_accu, reduce_accu

contains
   pure elemental subroutine add_accu_r(acc, t)
      implicit none
      type(accumulator_type), intent(inout) :: acc
      real(ki), intent(in) :: t

      real(ki) :: r, d
      integer :: e, i
      real(ki) :: radix_ki
      
      radix_ki = scale(1.0_ki, 1)

      r = fraction(t)
      e = exponent(t)
      i = e

      do while (r .ne. 0.0_ki .and. i .gt. min_ex_ki)
         r = scale(r, 1)
         i = i - 1
         ! The following two lines extract the first
         ! digit from the number. Hence, d plays the
         ! role of a_{ij} in the original publication.
         d = aint(r)
         r = r - d

         ! This is step 3 of the original algorithm:
         ! Add the digit to the according accummulator.
         acc%a(i) = acc%a(i) + d
      end do
   end  subroutine add_accu_r

   pure elemental subroutine add_accu_c(acc_re, acc_im, t)
      implicit none
      type(accumulator_type), intent(inout) :: acc_re, acc_im
      complex(ki), intent(in) :: t
         
      call add_accu(acc_re, real(t, ki))
      call add_accu(acc_im, aimag(t))
   end  subroutine add_accu_c

   pure elemental function reduce_accu_r(acc) result(f)
      ! This routine is step 4 of the original algorithm:
      ! Sum in decreasing order.
      implicit none
      type(accumulator_type), intent(in) :: acc
      real(ki) :: f

      integer :: e

      f = 0.0_ki

      do e = max_ex_ki, min_ex_ki, -1
         if (acc%a(e) .ne. 0.0_ki) f = f + scale(acc%a(e), e)
      end do
   end  function reduce_accu_r

   pure elemental function reduce_accu_c(acc_re, acc_im) result(f)
      implicit none
      type(accumulator_type), intent(in) :: acc_re, acc_im
      complex(ki) :: f

      f = cmplx(reduce_accu_r(acc_re), reduce_accu_r(acc_im), ki)
   end  function reduce_accu_c

end module maccu
