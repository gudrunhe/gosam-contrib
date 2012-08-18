! 
!****h* src/module/s_matrix_type
! NAME
!
!  Module s_matrix_type
!
! USAGE
!
!  use s_matrix_type
!
! DESCRIPTION
!
!  This module contains a type definition for the kinematic s_matrix, 
!  intended to mimic a run-time polymorphism.
!
! OUTPUT
!
!  This module exports the derived type:
!  * s_matrix_poly 
!
!  One function:
!  * assign_s_matrix -- associates the pointers in the s_matrix to a given real or complex matrix
!   
!  Subroutines:
!  * nullify_s_matrix -- nullifies the pointers
!  * set_s_matrix_bits -- sets an integer which describes the positions of complex masses in the s matrix.
!                         sets an integer which describes the positions of zero mass entries.
!  * fill_s_matrix -- fills the real array associated with a complex array.
!
! USES
!
!  * precision_golem (src/module/preci_double.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * equal (src/module/equal.f90)
!  * constante (src/module/constante.f90)
!
!*****
module s_matrix_type
  !
  use precision_golem, only: ki
  use sortie_erreur
  use equal
  use constante, only : zero
  !
  implicit none
  !
  !
  !****t* src/module/s_matrix_type/s_matrix_poly
  ! NAME
  !   s_matrix_poly
  !
  ! SYNOPSIS
  !   type s_matrix_poly
  !
  ! SOURCE
  type s_matrix_poly
     !
     real(ki), dimension(:,:), pointer :: pt_real
     complex(ki), dimension(:,:), pointer :: pt_cmplx
     integer :: b_cmplx, b_zero
     !
  end type s_matrix_poly
  !
  ! NOTES
  !   * pt_real points to a real array (s_mat_r) if associated.
  !   * pt_cmplx points to a complex array (s_mat_c) if associated
  !   * b_cmplx is a bit-integer encoding the positions of 
  !     complex mass entries in the S matrix.
  !   * b_zero is a bit-integer encoding the positions of
  !     vanishing masses.
  !
  !****
  !
  interface assign_s_matrix
     !
     module procedure assign_s_matrix_r
     module procedure assign_s_matrix_c
     !
  end interface
  !
  !
  private
  !
  !
  public :: assign_s_matrix, set_s_matrix_bits, nullify_s_matrix, fill_s_matrix,s_matrix_poly 
  !
  !
contains
  !
  !
  !
  !****f* src/module/s_matrix_type/assign_s_matrix
  ! NAME
  !
  !  assign_s_matrix
  !
  ! USAGE
  !
  !  assign_s_matrix(s_mat_r)
  !  assign_s_matrix(s_mat_c,s_mat_r)
  !
  ! DESCRIPTION
  !
  !  This function associates the global (type s_matrix_poly) s_mat_p
  !  with the given real or complex input matrix. In the case a complex
  !  matrix is given, a real matrix, which will contain the real part of the 
  !  complex matrix, has also be given as an argument.
  !
  ! INPUTS
  !
  !  A real and a complex matrix. Or just a real matrix.
  !
  !
  ! RETURN VALUE
  !
  !  a type (s_matrix_poly) is returned.
  !
  !
  !*****
  !
  function assign_s_matrix_r(s_mat_r) result (s_mat_p)
    real(ki), dimension(:,:), target, intent(in) :: s_mat_r
    type(s_matrix_poly) :: s_mat_p
    !
    s_mat_p%pt_real => s_mat_r
    nullify(s_mat_p%pt_cmplx)
    s_mat_p%b_cmplx = 0
    s_mat_p%b_zero = -1
    !
  end function assign_s_matrix_r
  !
  !
  function assign_s_matrix_c(s_mat_c, s_mat_r) result (s_mat_p)
    complex(ki), dimension(:,:), target, intent(in) :: s_mat_c
    real(ki), dimension(:,:), target, intent(in) :: s_mat_r
    type(s_matrix_poly) :: s_mat_p
    !
    s_mat_p%pt_cmplx => s_mat_c
    s_mat_p%pt_real => s_mat_r
    s_mat_p%b_cmplx = -1
    s_mat_p%b_zero = -1
    !
  end function assign_s_matrix_c 
  !
  !****f* src/module/s_matrix_type/fill_s_matrix
  ! NAME
  !
  !  Subroutine fill_s_matrix
  !
  ! USAGE
  !
  !  call fill_s_matrix(s_mat_p)
  !
  ! DESCRIPTION
  !
  !  This procedure fills the associated real array with the
  !  real entries of a complex s_matrix if the corresponding
  !  pointer to the complex array is associated.
  !
  ! INPUTS
  !
  !  *  a type (s_matrix_poly) object
  !
  ! RETURN VALUE
  !
  !  none
  !
  !
  !*****
  !
  subroutine fill_s_matrix(s_mat_p)
    type(s_matrix_poly) :: s_mat_p
    !
    if (associated(s_mat_p%pt_cmplx) ) then
       !
       s_mat_p%pt_real = real(s_mat_p%pt_cmplx,ki)
       !
    end if
    !
  end subroutine fill_s_matrix    
  !
  !****f* src/module/s_matrix_type/set_s_matrix_bits
  ! NAME
  !
  !  Subroutine set_s_matrix_bits
  !
  ! USAGE
  !
  !  call set_s_matrix_bits(s_mat_p,set_ref)
  !
  ! DESCRIPTION
  !
  !  This procedure checks the diagonal of the given matrix for complex
  !  entries with non-vanishing imaginary part and vanishing masses as well.
  !  The results are encoded in bit-integers included in the derived type s_mat_p. 
  !
  ! INPUTS
  !
  !  *  a type (s_matrix_poly) object
  !  *  a reference set of numbers associated with the diagonal entries.
  !     (typically an array (/1,2,...,n/) ).
  !
  ! RETURN VALUE
  !
  !  none
  !
  !
  !*****
  !
  subroutine set_s_matrix_bits(s_mat_poly,set_ref)
    type(s_matrix_poly), intent (inout) :: s_mat_poly
    integer, intent(in), dimension(:) :: set_ref
    integer, dimension(size(set_ref)) :: position_c, position_r
    complex(ki), dimension(size(set_ref)) :: diagonal
    real(ki), dimension(size(set_ref)) :: diag_imag, diag_real
    integer :: i,n
    !
    n = size(set_ref)
    position_c = 0
    position_r = 0
    !
    if (associated(s_mat_poly%pt_cmplx)) then
       !
       do i = 1, n
          diagonal(i) = s_mat_poly%pt_cmplx(i,i)
       end do
       !
       diag_imag = aimag(diagonal)
       !
       if (minval(diag_imag) .lt. 0._ki) then
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In subroutine set_s_matrix_bits:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'The S matrix contains masses with positive imaginary part!'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'This might lead to wrong results!'
          call catch_exception(1)
          !
       end if
       !
       where (diag_imag .gt. 2.0_ki*epsilon(1.0_ki)) position_c = set_ref
              ! twice epsilon to give consistency with for rounding-error size imaginary parts
       where (position_c .ne. 0) position_c = ibset(0,pos=position_c)
       s_mat_poly%b_cmplx = sum (position_c)
    end if
    !
    do i = 1, n
       diag_real(i) = s_mat_poly%pt_real(i,i)
    end do
    !
    where ( (position_c == 0) .and. (equal_real(diag_real,zero) ) ) position_r = set_ref
    !
    where (position_r .ne. 0) position_r = ibset(0,pos=position_r)
    s_mat_poly%b_zero = sum (position_r)
    !
  end subroutine set_s_matrix_bits
  !
  !
  !****f* src/module/s_matrix_type/nullify_s_matrix
  ! NAME
  !
  !  Subroutine nullify_s_matrix
  !
  ! USAGE
  !
  !  nullify_s_matrix(s_mat_p)
  !
  ! DESCRIPTION
  !
  !  This procedure nullifies the pointers in the input object.
  !
  ! INPUTS
  !
  !  *  type (s_matrix_poly) object
  !
  ! RETURN VALUE
  !
  !  none
  !
  !
  !*****
  !
  subroutine nullify_s_matrix(s_mat_p)
    type(s_matrix_poly) :: s_mat_p
    !
    if (associated(s_mat_p%pt_real)) then
       nullify(s_mat_p%pt_real)
    end if
    !
    if (associated(s_mat_p%pt_cmplx)) then
       nullify(s_mat_p%pt_cmplx)
    end if
    !
  end subroutine nullify_s_matrix
  !
end module s_matrix_type
