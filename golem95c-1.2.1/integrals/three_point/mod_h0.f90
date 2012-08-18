!
!****h* src/integrals/three_point/func_h0
! NAME
!
!  Module func_h0
!
! USAGE
!
!  use func_h0
!
! DESCRIPTION
!
!  This module is specific for the function h0 defined
!  by h0(x,alpha) = (-x-i lambda)^(alpha)/x
!  with alpha << 1. The three functions h0d, h0e and h0f
!  are defined as:
!  h0(x,alpha) = h0d(x) + alpha h0e(x) + alpha^2 h0f(x)
!  
!
! OUTPUT
!
!  This module exports three functions:
!  * h0d -- function
!  * h0e -- function
!  * h0f -- function
!  
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * parametre (src/module/parametre.f90)
!  * logarithme (src/module/z_log.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!
!*****
module func_h0
  !
  use precision_golem
  use parametre
  use logarithme
  use sortie_erreur
  implicit none
  !
  private 
  public :: h0d,h0e,h0f
  !
  contains
    !
    !****f* src/integrals/three_point/func_h0/h0d
    ! NAME
    !
    !  Function h0d
    !
    ! USAGE
    !
    !  real_dim2 = h0d(x)
    !
    ! DESCRIPTION
    !
    !  Compute the function 1/x
    !
    ! INPUTS
    !
    !  * x -- a real (type ki)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  It returns a real (type ki) array of rank 1 and shape 2
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function h0d(x)
      !
      real(ki), intent(in) :: x
      real(ki), dimension(2) :: h0d
      !
      h0d(1) = 1._ki/x
      h0d(2) = 0._ki
      !
    end function h0d
    !
    !****f* src/integrals/three_point/func_h0/h0e
    ! NAME
    !
    !  Function h0e
    !
    ! USAGE
    !
    !  real_dim2 = h0e(x)
    !
    ! DESCRIPTION
    !
    !  Compute the function ln(-x)/x
    !
    ! INPUTS
    !
    !  * x -- a real (type ki)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, the returned value depends on the global variable rat_or_tot_par
    !
    ! RETURN VALUE
    !
    !  It returns a real (type ki) array of rank 1 and shape 2
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function h0e(x)
      !
      real(ki), intent(in) :: x
      real(ki), dimension(2) :: h0e
      !
      if (rat_or_tot_par%tot_selected) then
        !
        h0e(1) = real(z_log(-x/mu2_scale_par,-1._ki)/x,ki)
        h0e(2) = aimag(z_log(-x/mu2_scale_par,-1._ki)/x)
        !
      else if (rat_or_tot_par%rat_selected) then
        !
        h0e = 0._ki
        !
      end if
      !
    end function h0e
    !
    !****f* src/integrals/three_point/func_h0/h0f
    ! NAME
    !
    !  Function h0f
    !
    ! USAGE
    !
    !  real_dim2 = h0f(x)
    !
    ! DESCRIPTION
    !
    !  Compute the function 1/2 ln(-x)^2/x
    !
    ! INPUTS
    !
    !  * x -- a real (type ki)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, the returned value depends on the global variable rat_or_tot_par
    !
    ! RETURN VALUE
    !
    !  It returns a real (type ki) array of rank 1 and shape 2
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function h0f(x)
      !
      real(ki), intent(in) :: x
      real(ki), dimension(2) :: h0f
      !
      if (rat_or_tot_par%tot_selected) then
        !
        h0f(1) = real(1._ki/2._ki*z_log2(-x/mu2_scale_par,-1._ki)/x,ki)
        h0f(2) = aimag(1._ki/2._ki*z_log2(-x/mu2_scale_par,-1._ki)/x)
        !
      else if (rat_or_tot_par%rat_selected) then
        !
        h0f = 0._ki
        !
      end if
      !
    end function h0f
    !
end module func_h0
