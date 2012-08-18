!
!****h* src/integrals/three_point/func_he
! NAME
!
!  Module func_he
!
! USAGE
!
!  use func_he
!
! DESCRIPTION
!
!  This module contains several functions for the computation of
!  int^1_0 dy y^(n-1)/(y*z1+(1-y)*z3) where z1 and z3 are complex numbers
!
! OUTPUT
!
!  This modules exports three functions:
!  * he -- a function
!  * he_gen -- a function
!  * he_c -- a function
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * numerical_evaluation (src/numerical/mod_numeric.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * parametre (src/module/parametre.f90)
!  * logarithme (src/module/z_log.f90)
!  * constante (src/module/constante.f90)
!
!*****
module func_he
  !
  use precision_golem
  use numerical_evaluation
  use sortie_erreur, only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
  use parametre, only : coupure_3p2m,rat_or_tot_par,tolerance,alpha_par,&
       & beta_par,lambda_par
  use logarithme, only : z_log
  use constante, only : i_,un,czero
  implicit none
  !
  real(ki) :: a1_glob,a3_glob,eps_glob
  complex(ki) :: a1_glob_c, a3_glob_c
  real(ki) :: plus_grand_glob
  integer :: expo_glob
  !
  private
  !
  interface he
     !
     module procedure he_rarg
     module procedure he_carg
     !
  end interface
  !
  public :: he,he_gen,he_c
  !
contains
  !
  !****f* src/integrals/three_point/func_he/he
  ! NAME
  !
  !  Function he
  !  Note that this function is an interface for two other functions
  !  he_rarg and he_carg
  !
  ! USAGE
  !
  !  real_dim2 = he(n,a1,a3)
  !
  ! DESCRIPTION
  !
  !  This function computes:
  !  - int^1_0 dy y^(n-1)/(y*z1+(1-y)*z3)
  !  where z1 = -a1 -i lambda and z3 = -a3 - i lambda
  !  For n=1, it is equal to: - (ln(z1)-ln(z3))/(z1-z3)
  !  compatible with the definition of HnE
  !  It switches to numerical evaluation if 
  !  |a1-a3|/max(|a1|,|a3|) < coupure_3p2m
  !
  ! INPUTS
  !
  !  * n -- an integer, the power of y in the integrand
  !  * a1 -- a real/complex (type ki), z1 (time -1)
  !  * a3 -- a real/complex (type ki), z3 (time -1)
  !  or
  !  * n -- an integer, the power of y in the integrand
  !  * a1 -- a complex (type ki), z1 (time -1)
  !  * a3 -- a complex (type ki), z3 (time -1)
  !
  ! SIDE EFFECTS
  !
  !  No side effect, the returned value depends on the global variables
  !  rat_or_tot_par, coupure_3p2m
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
  function he_rarg(n,a1,a3)
    !
    integer, intent(in) :: n
    real(ki), intent(in) :: a1,a3
    real(ki), dimension(2) :: he_rarg
    !
    complex(ki) :: rest
    complex(ki) :: abserr
    complex(ki), dimension(4) :: ver
    real(ki) :: g1,g3
    !
    plus_grand_glob = max(abs(a1),abs(a3))
    g1 = a1/plus_grand_glob
    g3 = a3/plus_grand_glob
    ! les variables a1_glob, a3_glob, expo_glob et eps_glob sont globales 
    a1_glob = -g1
    a3_glob = -g3
    expo_glob = n
    ! on choisit eps_glob de telle facon que le pole soit hors du contour
    eps_glob = -sign(un,a1-a3)
    !
    ! mettre une coupure d'ordre 1 !!!!!
    if (abs(g1-g3) > coupure_3p2m) then
       !
       ver = 0._ki
       !
       if (n >= 1) then
          !
          if (rat_or_tot_par%tot_selected) then
             !
             ver(1) = (z_log(-g1,-1._ki)-z_log(-g3,-1._ki))/(g1-g3)
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             ver(1) = 0._ki
             !
          end if
          !
       end if
       !
       if (n >= 2) then
          !
          ver(2) = (-g3*ver(1)+1._ki)/(g1-g3)
          !
       end if
       !
       if (n >= 3) then
          !
          ver(3) = (-g3*ver(2)+1._ki/2._ki)/(g1-g3)
          !
       end if
       !
       if (n >= 4) then
          !
          ver(4) = (-g3*ver(3)+1._ki/3._ki )/(g1-g3)
          !
       end if
       !
       he_rarg(1) = real(ver(n),ki)/plus_grand_glob
       he_rarg(2) = aimag(ver(n))/plus_grand_glob
       !
    else if ( (abs(g1-g3) <= coupure_3p2m) .and. &
         (rat_or_tot_par%tot_selected) ) then
       !
       origine_info_par = "he_arg"
       num_grand_b_info_par = n
       denom_grand_b_info_par = abs(a1-a3)
       !
       call generic_eval_numer(eval_numer_he,0._ki,1._ki,tolerance,rest,abserr)
       !
       tab_erreur_par(1)%a_imprimer = .true.
       tab_erreur_par(1)%chaine = 'In function he_rarg:'
       tab_erreur_par(2)%a_imprimer = .true.
       tab_erreur_par(2)%chaine = 'the error returned by adapt_gauss1 is: %z0'
       tab_erreur_par(2)%arg_comp = abserr
       call catch_exception(1)
       !
       he_rarg(1) = real(rest,ki)/plus_grand_glob
       he_rarg(2) = aimag(rest)/plus_grand_glob
       !
    else
       !
       tab_erreur_par(1)%a_imprimer = .true.
       tab_erreur_par(1)%chaine = 'In function he_rarg (file mod_he.f90)'
       tab_erreur_par(2)%a_imprimer = .true.
       tab_erreur_par(2)%chaine = 'the choice rat has been made, it is singular when a1=a3 %d0'
       tab_erreur_par(2)%arg_real=abs(g1-g3)
       call catch_exception(0)
       !
    end if
    !
  end function he_rarg
  !
  function he_carg(n,a1,a3)
    !
    integer, intent(in) :: n
    complex(ki), intent(in) :: a1,a3
    real(ki), dimension(2) :: he_carg
    !
    complex(ki) :: rest
    complex(ki) :: abserr
    complex(ki), dimension(4) :: ver
    complex(ki) :: g1,g3
    !
    !~ plus_grand_glob = max(abs(real(a1,ki)), abs(aimag(a1)), abs(real(a3,ki)), abs(aimag(a3)) )
    plus_grand_glob = max(abs(a1),abs(a3))
    g1 = a1/plus_grand_glob
    g3 = a3/plus_grand_glob
    ! les variables a1_glob, a3_glob, expo_glob et eps_glob sont globales 
    a1_glob_c = -g1
    a3_glob_c = -g3
    expo_glob = n
    ! mettre une coupure d'ordre 1 !!!!!
    if (abs(g1-g3) > coupure_3p2m) then
       !
       ver(:) = czero
       !
       if (n >= 1) then
          !
          if (rat_or_tot_par%tot_selected) then
             !
             ver(1) = (log(-g1)-log(-g3))/(g1-g3)
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             ver(1) = 0._ki
             !
          end if
          !
       end if
       !
       if (n >= 2) then
          !
          ver(2) = (-g3*ver(1)+1._ki)/(g1-g3)
          !
       end if
       !
       if (n >= 3) then
          !
          ver(3) = (-g3*ver(2)+1._ki/2._ki)/(g1-g3)
          !
       end if
       !
       if (n >= 4) then
          !
          ver(4) = (-g3*ver(3)+1._ki/3._ki )/(g1-g3)
          !
       end if
       !
       he_carg(1) = real(ver(n),ki)/plus_grand_glob
       he_carg(2) = aimag(ver(n))/plus_grand_glob
       !
    else if ( (abs(g1-g3) <= coupure_3p2m) .and. &
         (rat_or_tot_par%tot_selected) ) then
       !
       ! we choose eps_glob in such a way that the pole is outside the contour
       ! we are in the case that sign(Im(g3)) = sign(Im(g1)) and sign(Re(g3)) = sign(Re(g1))
       ! in this case, we choose eps_glob such that eps_glob*(Re(g1)-Re(g3)) < 0 if Im(g1) > 0
       ! and eps_glob*(Re(g1)-Re(g3)) > 0 if Im(g1) < 0
       if ( sign(un,aimag(g1)) == sign(un,aimag(g3)) ) then
         eps_glob = -sign(un,aimag(g1))*sign(un,real(g1,ki)-real(g3,ki))
       else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function he_carg (file mod_he.f90) Im(g1) and Im(g3) do not the same sign'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Im(g1): %z0'
          tab_erreur_par(2)%arg_comp = aimag(g1)
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'Im(g3): %z0'
          tab_erreur_par(3)%arg_comp = aimag(g3)
          call catch_exception(0)
       end if
       !
       origine_info_par = "he_carg"
       num_grand_b_info_par = n
       denom_grand_b_info_par = abs(a1-a3)
       !
       call generic_eval_numer(eval_numer_he_c,0._ki,1._ki,tolerance,rest,abserr)
       !
       tab_erreur_par(1)%a_imprimer = .true.
       tab_erreur_par(1)%chaine = 'In function he_carg:'
       tab_erreur_par(2)%a_imprimer = .true.
       tab_erreur_par(2)%chaine = 'the error returned by adapt_gauss1 is: %z0'
       tab_erreur_par(2)%arg_comp = abserr
       call catch_exception(1)
       !
       he_carg(1) = real(rest,ki)/plus_grand_glob
       he_carg(2) = aimag(rest)/plus_grand_glob
       !
    else
       !
       tab_erreur_par(1)%a_imprimer = .true.
       tab_erreur_par(1)%chaine = 'In function he_carg (file mod_he.f90)'
       tab_erreur_par(2)%a_imprimer = .true.
       tab_erreur_par(2)%chaine = 'the choice rat has been made, it is singular when a1=a3 %d0'
       tab_erreur_par(2)%arg_real=abs(g1-g3)
       call catch_exception(0)
       !
    end if
    !
  end function he_carg
  !
  !
  !****if* src/integrals/three_point/func_he/eval_numer_he
  ! NAME
  !
  !  Function eval_numer_he
  !
  ! USAGE
  !
  !  complex = eval_numer_he(u)
  !
  ! DESCRIPTION
  !
  !  This is the integrand for the numerical evaluation of he
  !
  ! INPUTS
  !
  !  * u -- a real (type ki), the integration variable
  !
  ! SIDE EFFECTS
  !
  !  No side effect. The variables a1_glob, a3_glob, expo_glob and eps_glob
  !  are global in this module whereas variables lambda_par,beta_par,
  !  alpha_par are given by the module parametre
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki)
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  function eval_numer_he(u)
    !
    real(ki), intent (in) :: u
    complex(ki) :: eval_numer_he
    !
    real(ki) :: x,y
    complex(ki) :: z,jacob
    !
    x = u
    y = -lambda_par*u**alpha_par*(1._ki-u)**beta_par
    jacob = 1._ki - eps_glob*i_*lambda_par*u**(alpha_par-1._ki)&
         *(1._ki-u)**(beta_par-1._ki)*(alpha_par*(1._ki-u)-beta_par*u)
    z = x + eps_glob*i_*y
    eval_numer_he = z**(expo_glob-1)/(z*a1_glob+(1._ki-z)*a3_glob)
    eval_numer_he = -eval_numer_he*jacob
    !
  end function eval_numer_he
  !
  !****if* src/integrals/three_point/func_he/eval_numer_he
  ! NAME
  !
  !  Function eval_numer_he_c
  !
  ! USAGE
  !
  !  complex = eval_numer_he_c(u)
  !
  ! DESCRIPTION
  !
  !  This is the integrand for the numerical evaluation of he_carg
  !
  ! INPUTS
  !
  !  * u -- a real (type ki), the integration variable
  !
  ! SIDE EFFECTS
  !
  !  No side effect. The variables a1_glob_c, a3_glob_c, expo_glob and eps_glob
  !  are global in this module whereas variables lambda_par,beta_par,
  !  alpha_par are given by the module parametre
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki)
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  function eval_numer_he_c(u)
    !
    real(ki), intent (in) :: u
    complex(ki) :: eval_numer_he_c
    !
    real(ki) :: x,y
    complex(ki) :: z,jacob
    !
    x = u
    y = -lambda_par*u**alpha_par*(1._ki-u)**beta_par
    jacob = 1._ki - eps_glob*i_*lambda_par*u**(alpha_par-1._ki)&
         *(1._ki-u)**(beta_par-1._ki)*(alpha_par*(1._ki-u)-beta_par*u)
    z = x + eps_glob*i_*y
    eval_numer_he_c = z**(expo_glob-1)/(z*a1_glob_c+(1._ki-z)*a3_glob_c)
    eval_numer_he_c = -eval_numer_he_c*jacob
    !
  end function eval_numer_he_c
  !
  !****f* src/integrals/three_point/func_he/he_c
  ! NAME
  !
  !  Function he_c
  !
  ! USAGE
  !
  !  complex = he_c(n,a1,a3)
  !
  ! DESCRIPTION
  !
  ! This function computes the same thing as he
  ! but it returns a complex instead of a real array of rank 1 and shape 2
  !
  ! INPUTS
  !
  !  * n -- an integer, the power of y in the integrand
  !  * a1 -- a real (type ki), the real part of z1 (time -1)
  !  * a3 -- a real (type ki), the real part of z3 (time -1)
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki)
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  function he_c(n,a1,a3)
    !
    integer, intent(in) :: n
    real(ki), intent(in) :: a1,a3
    complex(ki) :: he_c
    !
    real(ki), dimension(2) :: temp
    !
    temp = he(n,a1,a3)
    he_c = cmplx(temp(1),temp(2),ki)
    !
  end function he_c
  !
  !****f* src/integrals/three_point/func_he/he_gen
  ! NAME
  !
  !  Function he_gen
  !
  ! USAGE
  !
  !  real_dim2 = he_gen(n,a1,b1,a3,b3)
  !
  ! DESCRIPTION
  !
  !  This function computes:
  !  int^1_0 dy y^n/(y*z1+(1-y)*z3)
  !  where z1 = a1 + i b1 and z3 = a3 + i b3
  !  For n=1, it is equal to: (ln(z1)-ln(z3))/(z1-z3)
  !  It switches to numerical evaluation if 
  !  |a1-a3|/max(|a1|,|a3|) < coupure_3p2m
  !
  ! INPUTS
  !
  !  * n -- an integer, the power of y in the integrand
  !  * a1 -- a real (type ki), the real part of z1 
  !  * b1 -- a real (type ki), the imaginary part of z1 
  !  * a3 -- a real (type ki), the real part of z3 
  !  * b3 -- a real (type ki), the imaginary part of z3 
  !
  ! SIDE EFFECTS
  !
  !  No side effect, the returned value depends on the global variables
  !  rat_or_tot_par, coupure_3p2m
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
  function he_gen(n,a1,b1,a3,b3)
    !
    integer, intent(in) :: n
    real(ki), intent(in) :: a1,b1,a3,b3
    real(ki), dimension(2) :: he_gen
    !
    complex(ki) :: rest
    complex(ki) :: abserr
    complex(ki), dimension(4) :: ver
    !
    plus_grand_glob = max(abs(a1),abs(a3))
    ! les variables a1_gen_glob, a3_gen_glob, expo_gen_glob et eps_gen_glob
    ! sont globales
    a1_glob = a1/plus_grand_glob
    a3_glob = a3/plus_grand_glob
    expo_glob = n
    ! on choisit eps de telle facon que le pole soit hors du contour
    eps_glob = sign(un,b1*a3-b3*a1)
    !
    ! mettre une coupure d'ordre 1 !!!!!
    if (abs(a1_glob-a3_glob) > coupure_3p2m) then
       !
       ver(1) = (z_log(a1,b1)-z_log(a3,b3))/(a1-a3)
       ver(2) = (-a3*ver(1)+1._ki)/(a1-a3)
       ver(3) = (-a3*ver(2)+1._ki/2._ki)/(a1-a3)
       ver(4) = (-a3*ver(3)+1._ki/3._ki )/(a1-a3)
       !
       rest = ver(n)
       !
    else
       !
       origine_info_par = "he_gen"
       num_grand_b_info_par = n
       denom_grand_b_info_par = abs(a1-a3)
       !
       call generic_eval_numer(eval_numer_he_gen,0._ki,1._ki,tolerance,rest,abserr)
       !
       tab_erreur_par(1)%a_imprimer = .true.
       tab_erreur_par(1)%chaine = 'In function he_gen:'
       tab_erreur_par(2)%a_imprimer = .true.
       tab_erreur_par(2)%chaine = 'the error returned by adapt_gauss1 is: %z0'
       tab_erreur_par(2)%arg_comp = abserr
       call catch_exception(1)
       !
    end if
    !
    he_gen(1) = real(rest,ki)/plus_grand_glob
    he_gen(2) = aimag(rest)/plus_grand_glob
    !
  end function he_gen
  !
  !****if* src/integrals/three_point/func_he/eval_numer_he_gen
  ! NAME
  !
  !  Function eval_numer_he_gen
  !
  ! USAGE
  !
  !  complex = eval_numer_he_gen(u)
  !
  ! DESCRIPTION
  !
  !  This is the integrand for the numerical evaluation of he_gen
  !
  ! INPUTS
  !
  !  * u -- a real (type ki), the integration variable
  !
  ! SIDE EFFECTS
  !
  !  No side effect. The variables a1_glob, a3_glob, expo_glob and eps_glob
  !  are global in this module whereas variables lambda_par,beta_par,
  !  alpha_par are given by the module parametre
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki)
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  function eval_numer_he_gen(u)
    !
    real(ki), intent (in) :: u
    complex(ki) :: eval_numer_he_gen
    !
    real(ki) :: x,y
    complex(ki) :: z,jacob
    !
    x = u
    y = -lambda_par*u**alpha_par*(1._ki-u)**beta_par
    jacob = 1._ki - eps_glob*i_*lambda_par*u**(alpha_par-1._ki)&
         *(1._ki-u)**(beta_par-1._ki)*(alpha_par*(1._ki-u)-beta_par*u)
    z = x + eps_glob*i_*y
    eval_numer_he_gen = z**(expo_glob-1)/(z*a1_glob+(1._ki-z)*a3_glob)
    eval_numer_he_gen = eval_numer_he_gen*jacob
    !
  end function eval_numer_he_gen
  !
end module func_he
