!
!****h* src/integral/three_point/function_3p1m_2mi
! NAME
!
!  Module function_3p1m_2mi
!
! USAGE
!
!  use function_3p1m_2mi
!
! DESCRIPTION
!
!  This module is used to compute the one off-shell external leg two internal mass three point function
!  with/without Feynman parameters in n, n+2 dimensions
!
! OUTPUT
!
!  This module exports two functions:
!  * f3p1m_2mi -- a function for the computation of the one off-shell external leg two internal mass three 
!    point function with/without Feynman parameters in n dimensions
!  * f3p1m_2mi_np2 -- a function for the computation of the one off-shell external leg two internal mass three 
!    point function with/without Feynman parameters in n+2 dimensions
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * func_gn (src/integrals/three_point/mod_gn.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90) 
!         only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
!  * constante (src/module/constante.f90) only : un,i_,pi6
!  * parametre (src/module/parametre.f90) only : rat_or_tot_par,mu2_scale_par
!  * array (src/module/array.f90) only : packb
!
!*****
module function_3p1m_2mi
  !
  use precision_golem
  use func_gn
  use sortie_erreur, only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
  use constante, only : un,i_,pi6
  use parametre, only : rat_or_tot_par,mu2_scale_par
  use array, only : packb
  implicit none
  !
  private 
  !
  public :: f3p1m_2mi, f3p1m_2mi_np2
  !
  contains
    !
    !
    !****f* src/integral/three_point/function_3p1m_2mi/f3p1m_2mi
    ! NAME
    !
    !  Function f3p1m_2mi
    !
    ! USAGE
    !
    !  real_dim6 = f3p1m_2mi(s13,m1_sq,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the one off-shell external leg two internal mass three point function in n dimensions
    !  with up to three Feynman parameters in the numerator.
    !  It retuns an array of 6 reals corresponding to the real/imaginary
    !  part of the coefficient of the 1/epsilon^2 term, real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s13 -- real (type ki), the value of the S matrix element corresponding to the external off-shell leg
    !  * m1_sq -- real (type ki), the value of the first internal mass squared
    !  * m3_sq -- real (type ki), the value of the second internal mass squared
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
    !
    !*****
    function f3p1m_2mi(s13,m1_sq,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s13,m1_sq,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(6) :: f3p1m_2mi
      !
      real(ki) :: a,b,c
      real(ki) :: lmu2
      real(ki) :: true_thresh,false_thresh
      logical :: dist
      !
      f3p1m_2mi = 0._ki
      !
      a = s13+m1_sq+m3_sq
      b = ( m1_sq-m3_sq-(s13+m1_sq+m3_sq) )
      c = m3_sq
      ! one tests if we are closer from the real threshold or from the false threshold
      true_thresh = (sqrt(m1_sq)+sqrt(m3_sq))**2
      false_thresh = (sqrt(m1_sq)-sqrt(m3_sq))**2
      dist = abs(a-true_thresh) <= abs(a-false_thresh)
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        f3p1m_2mi(3:4) = -0.5_ki*ge(1,a,b,c,dist)
        f3p1m_2mi(5:6) = -0.5_ki*gf(1,a,b,c,dist)
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        if (par3 == 1) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(2,a,b,c,dist)
          !
        else if (par3 == 2) then
          !
          f3p1m_2mi(3:4) = -0.5_ki*ge(1,a,b,c,dist)
          f3p1m_2mi(5:6) = -0.5_ki*gf(1,a,b,c,dist)+ge(1,a,b,c,dist)
          !
        else if (par3 == 3) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(2,a,-b-2*a,a+b+c,dist)
          !
        end if
        !
      ! cas avec deux parametres de feynman au numerateur
      else if ( (par1==0) ) then
        !
        if ( (par2 == 1) .and. (par3 == 1) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(3,a,b,c,dist)/2._ki
          !
        else if ( (par2 == 2) .and. (par3 == 2) ) then
          !
          f3p1m_2mi(3:4) = -0.5_ki*ge(1,a,b,c,dist)
          f3p1m_2mi(5:6) = -0.5_ki*gf(1,a,b,c,dist)+3._ki/2._ki*ge(1,a,b,c,dist)
          !
        else if ( (par2 == 3) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(3,a,-b-2*a,a+b+c,dist)/2._ki
          !
        else if ( (par2 == 1) .and. (par3 == 2) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(2,a,b,c,dist)/2._ki
          !
        else if ( (par2 == 1) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -(ge(2,a,b,c,dist)-ge(3,a,b,c,dist))/2._ki
          !
        else if ( (par2 == 2) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(2,a,-b-2*a,a+b+c,dist)/2._ki
          !
        end if
        !
      ! cas avec trois parametres de feynman au numerateur
      else
        !
        if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 1) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(4,a,b,c,dist)/3._ki
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          f3p1m_2mi(3:4) = -0.5_ki*ge(1,a,b,c,dist)
          f3p1m_2mi(5:6) = -0.5_ki*gf(1,a,b,c,dist)+11._ki/6._ki*ge(1,a,b,c,dist)
          !
        else if ( (par1 == 3) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(4,a,-b-2*a,a+b+c,dist)/3._ki
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 2) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(3,a,b,c,dist)/6._ki
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(2,a,b,c,dist)/3._ki
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -(ge(3,a,b,c,dist)-ge(4,a,b,c,dist))/3._ki
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(2,a,-b-2*a,a+b+c,dist)/3._ki
          !
        else if ( (par1 == 1) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -(ge(3,a,-b-2*a,a+b+c,dist)-ge(4,a,-b-2*a,a+b+c,dist))/3._ki
          !
        else if ( (par1 == 2) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -ge(3,a,-b-2*a,a+b+c,dist)/6._ki
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          f3p1m_2mi(3:4) = 0._ki
          f3p1m_2mi(5:6) = -(ge(2,a,b,c,dist)-ge(3,a,b,c,dist))/6._ki
          !
        end if
        !
      end if
      !
      ! On change \epsilon_{ir} en -\epsilon_{uv}
      !
      f3p1m_2mi(3:4) = -f3p1m_2mi(3:4)
      !
      ! On factorise r_{\gamma}
      !
      f3p1m_2mi(5:6) = f3p1m_2mi(5:6)+pi6*f3p1m_2mi(1:2)
      !
      ! on ajoute la dependence en mu^2
      !
      lmu2 = log(mu2_scale_par)
      f3p1m_2mi(5:6) = f3p1m_2mi(5:6) + f3p1m_2mi(3:4)*lmu2 + f3p1m_2mi(1:2)*lmu2**2/2._ki      
      f3p1m_2mi(3:4) = f3p1m_2mi(3:4) + f3p1m_2mi(1:2)*lmu2
      !
    end function f3p1m_2mi
    !
    !****f* src/integral/three_point/function_3p1m_2mi/f3p1m_2mi_np2
    ! NAME
    !
    !  Function f3p1m_2mi_np2
    !
    ! USAGE
    !
    !  real_dim4 = f3p1m_2mi_np2(s13,m1_sq,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the one off-shell external leg two internal mass three point function in n+2 dimensions. 
    !  with up to one Feynman parameter in the numerator.
    !  It retuns an array of 4 reals corresponding to the real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s13 -- real (type ki), the value of the S matrix element corresponding to the external off-shell leg
    !  * m1_sq -- real (type ki), the value of the first internal mass squared
    !  * m3_sq -- real (type ki), the value of the second internal mass squared
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
    !
    !*****
    function f3p1m_2mi_np2(s13,m1_sq,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s13,m1_sq,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(4) :: f3p1m_2mi_np2
      !
      real(ki) :: a,b,c
      real(ki) :: lmu2
      !
      f3p1m_2mi_np2 = 0._ki
      !
      a = s13+m1_sq+m3_sq
      b = ( m1_sq-m3_sq-(s13+m1_sq+m3_sq) )
      c = m3_sq
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        f3p1m_2mi_np2(1) = -1._ki/2._ki
        f3p1m_2mi_np2(2) = 0._ki
        f3p1m_2mi_np2(3:4) = gl(1,a,b,c)/2._ki
        f3p1m_2mi_np2(3) = f3p1m_2mi_np2(3)-2._ki/4._ki
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        f3p1m_2mi_np2(1) = -1._ki/6._ki
        f3p1m_2mi_np2(2) = 0._ki
        !
        if (par3 == 1) then
          !
          f3p1m_2mi_np2(3:4) = gl(2,a,b,c)/3._ki
          f3p1m_2mi_np2(3) = f3p1m_2mi_np2(3)-1._ki/9._ki
          !
        else if (par3 == 2) then
          !
          f3p1m_2mi_np2(3:4) = gl(1,a,b,c)/6._ki
          f3p1m_2mi_np2(3) = f3p1m_2mi_np2(3)-5._ki/18._ki
          !
        else if (par3 == 3) then
          !
          f3p1m_2mi_np2(3:4) = gl(2,a,-b-2*a,a+b+c)/3._ki
          f3p1m_2mi_np2(3) = f3p1m_2mi_np2(3)-1._ki/9._ki
          !
        end if
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f3p1m_2mi_np2:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'no need of 3-point integrals in 6 dimension &
                          &with more than one Feynman parameter in the numerator'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'The value of Feynman parameters in argument: %d1'
        tab_erreur_par(3)%arg_int_tab = (/packb((/par1,par2,par3/)),4/)
        call catch_exception(0)
        !
      end if
      !
      ! on ajoute la dependence en mu^2
      !
      lmu2 = log(mu2_scale_par)
      f3p1m_2mi_np2(3:4) = f3p1m_2mi_np2(3:4) + f3p1m_2mi_np2(1:2)*lmu2
      !
    end function f3p1m_2mi_np2
    !
end module function_3p1m_2mi
