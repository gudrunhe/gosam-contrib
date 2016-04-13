!****h* src/integrals/two_point/generic_function_1p
! NAME
!
!  Module generic_function_1p
!
! USAGE
!
!  use generic_function_1p
!
! DESCRIPTION
!
!  This module contains the generic routines to compute 
!  one point functions in n  dimensions
!
! OUTPUT
!
!  It exports one public routine:
!  * f1p -- a function to compute the one point function in n dimensions
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * array (src/module/array.f90)
!  * logarithme (src/module/z_log.f90)
!  * constante (src/module/constante.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * parametre (src/module/parametre.f90)
!  * s_matrix_type (src/module/s_matrix_type.f90)
!  * equal (src/module/equal.f90)
!
!*****
module generic_function_1p
  !
  use precision_golem
  use array
  use logarithme
  use constante, only:czero, zero
  use sortie_erreur
  use parametre
  use s_matrix_type
  use equal
  !
  implicit none
  !
  private
  !
  interface f1p
     !
     module procedure f1p_r, f1p_c
     module procedure f1p_p
     !
  end interface

  public :: f1p
  !
contains
  !
  !****f* src/integrals/one_point/generic_function_1p/f1p
  ! NAME
  !
  !  Function f1p
  !
  ! USAGE
  !
  !  real_dim4 = f1p(s_mat_p,b_pro,parf1)
  !
  ! DESCRIPTION
  !
  !  This function computes the generic two point function in n dimensions, 
  !  with or without Feynman parameters in the numerator
  !
  ! INPUTS
  !
  !  * s_mat_(r/c/p) -- a real/complex (type ki)/type(s_matrix_poly) array of rank 2, the S matrix
  !  * b_pro -- an integer which represents the set of the four unpinched
  !    propagators
  !  * parf1 -- an integer (optional), the label of the one Feynman parameter
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki) array of rank 1 and shape 2
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  !
  function f1p_p(s_mat_p,b_pro,parf1)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent(in) :: b_pro
    integer, intent(in), optional :: parf1
    complex(ki), dimension(2) :: f1p_p
    !
    if (iand(s_mat_p%b_cmplx, b_pro) .eq. 0 ) then
       !
       f1p_p = f1p_r(s_mat_p%pt_real, b_pro, parf1=parf1)
       !
    else
       !
       f1p_p = f1p_c(s_mat_p%pt_cmplx, b_pro, parf1=parf1)
       !
    end if
    !
  end function f1p_p
  !
  function f1p_r(s_mat_r,b_pro,parf1)
    !
    real(ki), intent (in), dimension(:,:) :: s_mat_r
    integer, intent (in) :: b_pro
    integer, intent (in), optional :: parf1
    complex(ki), dimension(2) :: f1p_r
    !
    integer :: par1
    real(ki) :: mass1
    integer :: m1
    integer, dimension(1) :: s
    !
    if (present(parf1)) then
       par1 = parf1
    else 
       par1 = 0
    end if
    !
    if (par1 /= 0) par1 = locateb(par1, b_pro)
    !
    if (par1 == -1) then
       !
       f1p_r(:) = czero
       !
    else
       !
       s = unpackb(b_pro,countb(b_pro))
       !
       m1 = s(1)
       !
       mass1 = -s_mat_r(m1,m1)/2._ki
       !
       if ( equal_real(mass1,zero) ) then
          !
          f1p_r(:) = czero
          !
       else
          !
          if  (par1 == 0) then
             ! 
             f1p_r(1) = cmplx(mass1,0._ki,ki)
             !
             if (rat_or_tot_par%tot_selected) then
                !
                f1p_r(2) = mass1*(1._ki - z_log(mass1/mu2_scale_par,-1._ki) )
                !
             else if (rat_or_tot_par%rat_selected) then
                !
                f1p_r(2) = cmplx(mass1,0._ki,ki)
                !
             end if
             !
          else if  (par1 /= 0) then
             !
             f1p_r(:) = czero
             !
          end if
          !
       end if
       !
    end if
    !
  end function f1p_r
  !
  function f1p_c(s_mat_c,b_pro,parf1)
    !
    complex(ki), intent (in), dimension(:,:) :: s_mat_c
    integer, intent (in) :: b_pro
    integer, intent (in), optional :: parf1
    complex(ki), dimension(2) :: f1p_c
    !
    integer :: par1
    complex(ki) :: mass1
    integer :: m1
    integer, dimension(1) :: s
    !
    if (present(parf1)) then
       par1 = parf1
    else 
       par1 = 0
    end if
    !
    if (par1 /= 0) par1 = locateb(par1, b_pro)
    !
    if (par1 == -1) then
       !
       f1p_c(:) = czero
       !
    else
       !
       s = unpackb(b_pro,countb(b_pro))
       !
       m1 = s(1)
       !
       mass1 = -s_mat_c(m1,m1)/2._ki
       !
       ! This function is only called with non_vanishing mass1
       !
       if  (par1 == 0) then
          ! 
          f1p_c(1) = mass1
          !
          if (rat_or_tot_par%tot_selected) then
             !
             f1p_c(2) = mass1*(1._ki - z_log(mass1/mu2_scale_par,-1._ki) )
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f1p_c(2) = mass1
             !
          end if
          !
       else if  (par1 /= 0) then
          !
          f1p_c(:) = czero
          !
       end if
       !
    end if
    !
  end function f1p_c
  !
end module generic_function_1p
