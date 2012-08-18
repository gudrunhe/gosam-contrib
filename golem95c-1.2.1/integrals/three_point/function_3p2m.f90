!
!****h* src/integral/three_point/function_3p2m
! NAME
!
!  Module function_3p2m
!
! USAGE
!
!  use function_3p2m
!
! DESCRIPTION
!
!  This module is used to compute the two off-shell external leg three point function
!  with no internal mass with/without Feynman parameters in n, n+2 dimensions
!
! OUTPUT
!
!  This module exports two functions:
!  * f3p2m -- a function for the computation of the two off-shell external leg three 
!    point function with/without Feynman parameters in n dimensions
!  * f3p2m_np2 -- a function for the computation of the two off-shell external leg three 
!    point function with/without Feynman parameters in n+2 dimensions
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * logarithme (src/module/z_log.f90)
!  * func_h0 (src/integrals/three_point/mod_h0.f90)
!  * func_he (src/integrals/three_point/mod_he.f90)
!  * func_hf (src/integrals/three_point/mod_hf.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!
!*****
module function_3p2m
  !
  use precision_golem
  use logarithme
  use func_h0
  use func_he
  use func_hf
  use sortie_erreur
  implicit none
  !
  private 
  !
  public :: f3p2m, f3p2m_np2
  !
  contains
    !
    !
    !****f* src/integral/three_point/function_3p2m/f3p2m
    ! NAME
    !
    !  Function f3p2m
    !
    ! USAGE
    !
    !  real_dim6 = f3p2m(s13,s23,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the two off-shell external leg three point function in n dimensions
    !  with up to three Feynman parameters in the numerator.
    !  It retuns an array of 6 reals corresponding to the real/imaginary
    !  part of the coefficient of the 1/epsilon^2 term, real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s13 -- real (type ki), the value of the S matrix element corresponding to the first external off-shell leg
    !  * s23 -- real (type ki), the value of the S matrix element corresponding to the second external off-shell leg
    !  * par1 -- an integer, the label of the third Feynman parameter
    !  * par2 -- an integer, the label of the second Feynman parameter
    !  * par3 -- an integer, the label of the first Feynman parameter
    !  Note that par1,par2 and par3 are supposed to be ordered, i.e.
    !  par1 <= par2 <= par3, note also that put zero for par1, par2 or par3
    !  if this Feynman parameter does not exist.
    !  Use the routine tri_int(t_in,t_out) to order the labels in the module 
    !  tri_croissant (src/module/tri.f90)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  An real (type ki) array of rank 1 and shape 6 corresponding to 
    !  the real/imaginary part of the coefficient of the 1/epsilon^2 term,
    !  real/imaginary part of the coefficient of the 1/epsilon term
    !  and the real/imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    ! two mass three point function without Feynman parameters 
    ! f3p2m(s13,s23,0,0,0) 
    ! with one Feynman parameter at the numerator z_1 
    ! f3p2m(s13,s23,0,0,1)
    ! with three Feynman parameters at the numerator z_2^2 z_3
    ! f3p2m(s13,s23,2,2,3) 
    !
    !*****
    function f3p2m(s23,s13,par1,par2,par3)
      !
      real(ki), intent (in) :: s23,s13
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(6) :: f3p2m
      !
      f3p2m = 0._ki
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        f3p2m(3:4) = he(1,s13,s23)
        f3p2m(5:6) = hf(1,s13,s23)
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        if (par3 == 1) then
          !
          f3p2m(3:4) = he(2,s13,s23)
          f3p2m(5:6) = hf(2,s13,s23)-f3p2m(3:4)
          !
        else if (par3 == 2) then
          !
          f3p2m(3:4) = he(2,s23,s13)
          f3p2m(5:6) = hf(2,s23,s13)-f3p2m(3:4)
          !
        else if (par3 == 3) then
          !
          f3p2m(5:6) = he(1,s13,s23)
          !
        end if
        !
      ! cas avec deux parametres de feynman au numerateur
      else if ( (par1==0) ) then
        !
        if ( (par2 == 1) .and. (par3 == 1) ) then
          !
          f3p2m(3:4) = he(3,s13,s23)
          f3p2m(5:6) = hf(3,s13,s23)-3._ki/2._ki*f3p2m(3:4)
          !
        else if ( (par2 == 2) .and. (par3 == 2) ) then
          !
          f3p2m(3:4) = he(3,s23,s13)
          f3p2m(5:6) = hf(3,s23,s13)-3._ki/2._ki*f3p2m(3:4)
          !
        else if ( (par2 == 3) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(1,s13,s23)/2._ki
          !
       else if ( (par2 == 1) .and. (par3 == 2) ) then
          !
          f3p2m(3:4) = he(2,s13,s23)-he(3,s13,s23)
          f3p2m(5:6) = hf(2,s13,s23)-hf(3,s13,s23)-3._ki/2._ki*f3p2m(3:4)
          !
        else if ( (par2 == 1) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(2,s13,s23)/2._ki
          !
       else if ( (par2 == 2) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(2,s23,s13)/2._ki
          !
        end if
        !
      ! cas avec trois parametres de feynman au numerateur
      else
        !
        if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 1) ) then
          !
          f3p2m(3:4) = he(4,s13,s23)
          f3p2m(5:6) = hf(4,s13,s23)-11._ki/6._ki*f3p2m(3:4)
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          f3p2m(3:4) = he(4,s23,s13)
          f3p2m(5:6) = hf(4,s23,s13)-11._ki/6._ki*f3p2m(3:4)
          !
        else if ( (par1 == 3) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(1,s13,s23)/3._ki
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 2) ) then
          !
          f3p2m(3:4) = he(3,s13,s23)-he(4,s13,s23)
          f3p2m(5:6) = hf(3,s13,s23)-hf(4,s13,s23)-11._ki/6._ki*f3p2m(3:4)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          f3p2m(3:4) = he(3,s23,s13)-he(4,s23,s13)
          f3p2m(5:6) = hf(3,s23,s13)-hf(4,s23,s13)-11._ki/6._ki*f3p2m(3:4)
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(3,s13,s23)/3._ki
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(3,s23,s13)/3._ki
          !
        else if ( (par1 == 1) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(2,s13,s23)/6._ki
          !
        else if ( (par1 == 2) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(2,s23,s13)/6._ki
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          f3p2m(5:6) = he(2,s13,s23)/3._ki-he(3,s13,s23)/3._ki
          !
        end if
        !
      end if
      !
      ! On change \epsilon_{ir} en -\epsilon_{uv}
      !
      f3p2m(3:4) = -f3p2m(3:4)
      !
    end function f3p2m
    !
    !
    !****f* src/integral/three_point/function_3p2m/f3p2m_np2
    ! NAME
    !
    !  Function f3p2m_np2
    !
    ! USAGE
    !
    !  real_dim4 = f3p2m_np2(s13,s23,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the two off-shell external leg three point function in n+2 dimensions. 
    !  with up to one Feynman parameter in the numerator.
    !  It retuns an array of 4 reals corresponding to the real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s13 -- real (type ki), the value of the S matrix element corresponding to the first external off-shell leg
    !  * s23 -- real (type ki), the value of the S matrix element corresponding to the second external off-shell leg
    !  * par1 -- an integer, the label of the third Feynman parameter = 0
    !  * par2 -- an integer, the label of the second Feynman parameter = 0
    !  * par3 -- an integer, the label of the first Feynman parameter
    !  Note that par1,par2 and par3 are supposed to be ordered, i.e.
    !  par1 <= par2 <= par3, note also that put zero for par1, par2 or par3
    !  if this Feynman parameter does not exist.
    !  Use the routine tri_int(t_in,t_out) to order the labels in the module 
    !  tri_croissant (src/module/tri.f90)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  An real (type ki) array of rank 1 and shape 4 corresponding to 
    !  the real/imaginary part of the coefficient of the 1/epsilon term
    !  and the real/imaginary part of the constant term. If par1 and/or par2
    !  are different from zero, an error is returned.
    !
    ! EXAMPLE
    !
    ! two mass three point function without Feynman parameters 
    ! f3p2m_np2(s13,s23,0,0,0) 
    ! with one Feynman parameter at the numerator z_1 
    ! f3p2m_np2(s13,s23,0,0,1)
    !
    !*****
    function f3p2m_np2(s23,s13,par1,par2,par3)
      !
      real(ki), intent (in) :: s23,s13
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(4) :: f3p2m_np2
      !
      f3p2m_np2 = 0._ki
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        f3p2m_np2(1) = -1._ki/2._ki
        f3p2m_np2(2) = 0._ki
        f3p2m_np2(3:4) = (s13*h0e(s13)+s23*he(1,s13,s23))/2._ki
        f3p2m_np2(3) = f3p2m_np2(3)-3._ki/2._ki
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        f3p2m_np2(1) = -1._ki/6._ki
        f3p2m_np2(2) = 0._ki
        !
        if (par3 == 1) then
          !
          f3p2m_np2(3:4) = (s13*h0e(s13)+s23*he(2,s13,s23))/6._ki
          !
        else if (par3 == 2) then
          !
          f3p2m_np2(3:4) = (s23*h0e(s23)+s13*he(2,s23,s13))/6._ki
          !
        else if (par3 == 3) then
          !
          f3p2m_np2(3:4) = (s13*h0e(s13)+s23*he(1,s13,s23))/6._ki
          !
        end if
        !
        f3p2m_np2(3) = f3p2m_np2(3)-8._ki/18._ki
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'error in function f3p2m_np2'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'no need of two mass six dimensional 3-point function &
                          &with more than one Feynman parameter in the numerator'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'Feynman param 1: %d1'
        tab_erreur_par(3)%arg_int = par1
        tab_erreur_par(4)%a_imprimer = .true.
        tab_erreur_par(4)%chaine = 'Feynman param 2: %d1'
        tab_erreur_par(4)%arg_int = par2
        tab_erreur_par(5)%a_imprimer = .true.
        tab_erreur_par(5)%chaine = 'Feynman param 3: %d1'
        tab_erreur_par(5)%arg_int = par3
        call catch_exception(0)
        !
      end if
      !
    end function f3p2m_np2
    !
end module function_3p2m
