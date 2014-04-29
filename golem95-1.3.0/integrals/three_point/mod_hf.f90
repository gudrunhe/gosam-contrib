!
!****h* src/integrals/three_point/func_hf
! NAME
!
!  Module func_hf
!
! USAGE
!
!  use func_hf
!
! DESCRIPTION
!
!  This module contains several functions for the computation of
!  - int^1_0 dy y^n*ln(y*z1+(1-y)*z3)/(y*z1+(1-y)*z3) where z1 and 
!  z3 are complex numbers
!
! OUTPUT
!
!  This modules exports three functions:
!  * hf -- a function
!  * hf_gen -- a function
!  * hf_c -- a function
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
module func_hf
  use precision_golem
  use numerical_evaluation
  use sortie_erreur, only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
  use parametre, only : coupure_3p2m,rat_or_tot_par,tolerance,alpha_par,beta_par,lambda_par,mu2_scale_par
  use logarithme, only : z_log,z_log2
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
  interface hf
     !
     module procedure hf_rarg
     module procedure hf_carg
     !
  end interface
  !
  public :: hf,hf_gen,hf_c
  !
contains
  !
  !****f* src/integrals/three_point/func_hf/hf
  ! NAME
  !
  !  Function hf
  !  Note that this function is an interface for two other functions
  !  hf_rarg and hf_carg
  !
  ! USAGE
  !
  !  real_dim2 = hf(n,a1,a3)
  !
  ! DESCRIPTION
  !
  !  This function computes:
  !  - int^1_0 dy y^(n-1)*ln(y*z1+(1-y)*z3)/(y*z1+(1-y)*z3)
  !  where z1 = a1 + i b1 and z3 = a3 + i b3
  !  For n=1, it is equal to: -(ln^2(z1)-ln^2(z3))/2/(z1-z3)
  !  compatible with the definition of HnF
  !  It switches to numerical evaluation if 
  !  |a1-a3|/max(|a1|,|a3|) < coupure_3p2m
  !
  ! INPUTS
  !
  !  * n -- an integer, the power of y in the integrand
  !  * a1 -- a real (type ki), the real part of z1 (time -1)
  !  * a3 -- a real (type ki), the real part of z3 (time -1)
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
  function hf_rarg(n,a1,a3)
    !
    integer, intent(in) :: n
    real(ki), intent(in) :: a1,a3
    real(ki), dimension(2) :: hf_rarg
    !
    complex(ki) :: rest
    complex(ki) :: abserr
    complex(ki), dimension(4) :: ver,verm,vert
    real(ki) :: g1,g3,lm
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
    if (rat_or_tot_par%tot_selected) then
       !
       if (abs(g1-g3) > coupure_3p2m) then
          !
          lm = log(plus_grand_glob/mu2_scale_par)
          ver = 0._ki
          verm = 0._ki
          vert = 0._ki
          !
          if (n >= 1) then
             !
             ver(1) = ( z_log2(-g1,-1._ki) - z_log2(-g3,-1._ki) )/(g1-g3)/2._ki
             verm(1) = (z_log(-g1,-1._ki)-z_log(-g3,-1._ki))/(g1-g3)
             vert(1) = ver(1)+lm*verm(1)
             !
          end if
          !
          if (n >= 2) then
             !
             ver(2) = ( -g3*ver(1) - 1._ki + ( g1*z_log(-g1,-1._ki) &
                  - g3*z_log(-g3,-1._ki) )/(g1-g3) )/(g1-g3)
             verm(2) = (-g3*verm(1)+1._ki)/(g1-g3)
             vert(2) = ver(2)+lm*verm(2)
             !
          end if
          !
          if (n >= 3) then
             !
             ver(3) = ( -g3*ver(2) + 1._ki/4._ki*( -3._ki*g3**2 &
                  + 2*z_log(-g3,-1._ki)*g3**2 - g1**2 + 4._ki*g1*g3 &
                  + 2._ki*z_log(-g1,-1._ki)*g1**2  &
                  - 4*z_log(-g1,-1._ki)*g1*g3 )/(g1-g3)**2 )/(g1-g3)
             verm(3) = (-g3*verm(2)+1._ki/2._ki)/(g1-g3)
             vert(3) = ver(3)+lm*verm(3)
             !
          end if
          !
          if (n >= 4) then
             !
             ver(4) = ( -g3*ver(3) + 1._ki/18._ki*( -6._ki*z_log(-g3,-1._ki)*g3**3 &
                  + 11._ki*g3**3 - 18._ki*g3**2*g1 - 2._ki*g1**3  &
                  + 9._ki*g1**2*g3 - 18._ki*z_log(-g1,-1._ki)*g1**2*g3 &
                  + 18._ki*z_log(-g1,-1._ki)*g1*g3**2 &
                  + 6._ki*z_log(-g1,-1._ki)*g1**3 )/(g1-g3)**3 )/(g1-g3)
             verm(4) = (-g3*verm(3)+1._ki/3._ki )/(g1-g3)
             vert(4) = ver(4)+lm*verm(4)
             !
          end if
          !
          rest = vert(n)
          !
       else if ( (abs(g1-g3) <= coupure_3p2m) .and. (abs(g1-g3) >= tiny(g1)) ) then
          !
          origine_info_par = "hf"
          num_grand_b_info_par = n
          denom_grand_b_info_par = abs(a1-a3)
          !
          call generic_eval_numer(eval_numer_hf,0._ki,1._ki,tolerance,rest,abserr)
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function hf_rarg:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the error returned by adapt_gauss1 is: %z0'
          tab_erreur_par(2)%arg_comp = abserr
          call catch_exception(1)
          !
       else if ( (abs(g1-g3) <= tiny(g1)) ) then
          !
          rest = (z_log(-g1,-1._ki) + log(plus_grand_glob/mu2_scale_par))/(g1*real(n,ki))
          !
       end if
       !
       hf_rarg(1) = real(rest,ki)/plus_grand_glob
       hf_rarg(2) = aimag(rest)/plus_grand_glob
       !
    else if (rat_or_tot_par%rat_selected) then
       !
       if (abs(g1-g3) > coupure_3p2m) then
          !
          ver = 0._ki
          vert = 0._ki
          !
          if (n >= 1) then
             ver(1) = 0._ki
             vert(1) = 0._ki
          end if
          !
          if (n >= 2) then
             !
             ver(2) = ( -g3*ver(1) - 1._ki )/(g1-g3)
             vert(2) = ver(2)
             !
          end if
          !
          if (n >= 3) then
             !
             ver(3) = ( -g3*ver(2) + 1._ki/4._ki*( -3._ki*g3**2 &
                  - g1**2 + 4._ki*g1*g3 )/(g1-g3)**2 )/(g1-g3)
             vert(3) = ver(3)
             !
          end if
          !
          if (n >= 4) then
             !
             ver(4) = ( -g3*ver(3) + 1._ki/18._ki*(  &
                  + 11._ki*g3**3 - 18._ki*g3**2*g1 - 2._ki*g1**3  &
                  + 9._ki*g1**2*g3 )/(g1-g3)**3 )/(g1-g3)
             vert(4) = ver(4)
             !
          end if
          !
          rest = vert(n)
          !
          hf_rarg(1) = real(rest,ki)/plus_grand_glob
          hf_rarg(2) = aimag(rest)/plus_grand_glob
          !
       else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function hf_rarg (file mod_hf.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the choice rat has been made, it is singular when a1=a3 %d0'
          tab_erreur_par(2)%arg_real=abs(g1-g3)
          call catch_exception(0)
          !
       end if
       !
    end if
    !
  end function hf_rarg
  !
  function hf_carg(n,a1,a3)
    !
    integer, intent(in) :: n
    complex(ki), intent(in) :: a1,a3
    real(ki), dimension(2) :: hf_carg
    !
    complex(ki) :: rest
    complex(ki) :: abserr
    complex(ki), dimension(4) :: ver,verm,vert
    complex(ki) :: g1,g3
    real(ki) :: lm
    !
    !~ plus_grand_glob = max(abs(real(a1,ki)), abs(aimag(a1)), abs(real(a3,ki)), abs(aimag(a3)) )
    plus_grand_glob = max(abs(a1), abs(a3))
    g1 = a1/plus_grand_glob
    g3 = a3/plus_grand_glob
    ! les variables a1_glob, a3_glob, expo_glob et eps_glob sont globales 
    a1_glob_c = -g1
    a3_glob_c = -g3
    expo_glob = n
    !
    ! mettre une coupure d'ordre 1 !!!!!
    if (rat_or_tot_par%tot_selected) then
       !
       if (abs(g1-g3) > coupure_3p2m) then
          !
          lm = log(plus_grand_glob/mu2_scale_par)
          ver(:) = czero
          verm(:) = czero
          vert(:) = czero
          !
          if (n >= 1) then
             !
             ver(1) = ( log(-g1)**2 - log(-g3)**2 )/(g1-g3)/2._ki
             verm(1) = (log(-g1)-log(-g3))/(g1-g3)
             vert(1) = ver(1)+lm*verm(1)
             !
          end if
          !
          if (n >= 2) then
             !
             ver(2) = ( -g3*ver(1) - 1._ki + ( g1*log(-g1) &
                  - g3*log(-g3) )/(g1-g3) )/(g1-g3)
             verm(2) = (-g3*verm(1)+1._ki)/(g1-g3)
             vert(2) = ver(2)+lm*verm(2)
             !
          end if
          !
          if (n >= 3) then
             !
             ver(3) = ( -g3*ver(2) + 1._ki/4._ki*( -3._ki*g3**2 &
                  + 2*log(-g3)*g3**2 - g1**2 + 4._ki*g1*g3 &
                  + 2._ki*log(-g1)*g1**2  &
                  - 4*log(-g1)*g1*g3 )/(g1-g3)**2 )/(g1-g3)
             verm(3) = (-g3*verm(2)+1._ki/2._ki)/(g1-g3)
             vert(3) = ver(3)+lm*verm(3)
             !
          end if
          !
          if (n >= 4) then
             !
             ver(4) = ( -g3*ver(3) + 1._ki/18._ki*( -6._ki*log(-g3)*g3**3 &
                  + 11._ki*g3**3 - 18._ki*g3**2*g1 - 2._ki*g1**3  &
                  + 9._ki*g1**2*g3 - 18._ki*log(-g1)*g1**2*g3 &
                  + 18._ki*log(-g1)*g1*g3**2 &
                  + 6._ki*log(-g1)*g1**3 )/(g1-g3)**3 )/(g1-g3)
             verm(4) = (-g3*verm(3)+1._ki/3._ki )/(g1-g3)
             vert(4) = ver(4)+lm*verm(4)
             !
          end if
          !
          rest = vert(n)
          !
       else if ( (abs(g1-g3) <= coupure_3p2m) .and. (abs(g1-g3) >= tiny(real(g1,ki))) ) then
          !
       ! we choose eps_glob in such a way that the pole is outside the contour
       ! we are in the case that sign(Im(g3)) = sign(Im(g1)) and sign(Re(g3)) = sign(Re(g1))
       ! in this case, we choose eps_glob such that eps_glob*(Re(g1)-Re(g3)) < 0 if Im(g1) > 0
       ! and eps_glob*(Re(g1)-Re(g3)) > 0 if Im(g1) < 0
          if ( sign(un,aimag(g1)) == sign(un,aimag(g3)) ) then
            eps_glob = -sign(un,aimag(g1))*sign(un,real(g1,ki)-real(g3,ki))
          else
             tab_erreur_par(1)%a_imprimer = .true.
             tab_erreur_par(1)%chaine = 'In function hf_carg (file mod_hf.f90) Im(g1) and Im(g3) do not the same sign'
             tab_erreur_par(2)%a_imprimer = .true.
             tab_erreur_par(2)%chaine = 'Im(g1): %z0'
             tab_erreur_par(2)%arg_comp = aimag(g1)
             tab_erreur_par(3)%a_imprimer = .true.
             tab_erreur_par(3)%chaine = 'Im(g3): %z0'
             tab_erreur_par(3)%arg_comp = aimag(g3)
             call catch_exception(0)
          end if
          !
          origine_info_par = "hf_carg"
          num_grand_b_info_par = n
          denom_grand_b_info_par = abs(a1-a3)
          !
          call generic_eval_numer(eval_numer_hf_c,0._ki,1._ki,tolerance,rest,abserr)
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function hf:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the error returned by adapt_gauss1 is: %z0'
          tab_erreur_par(2)%arg_comp = abserr
          call catch_exception(1)
          !
       else if ( (abs(g1-g3) <= tiny(real(g1,ki))) ) then
          !
          rest = (log(-g1) + log(plus_grand_glob/mu2_scale_par))/(g1*real(n,ki))
          !
       end if
       !
       hf_carg(1) = real(rest,ki)/plus_grand_glob
       hf_carg(2) = aimag(rest)/plus_grand_glob
       !
    else if (rat_or_tot_par%rat_selected) then
       !
       if (abs(g1-g3) > coupure_3p2m) then
          !
          ver = 0._ki
          vert = 0._ki
          !
          if (n >= 1) then
             ver(1) = 0._ki
             vert(1) = 0._ki
          end if
          !
          if (n >= 2) then
             !
             ver(2) = ( -g3*ver(1) - 1._ki )/(g1-g3)
             vert(2) = ver(2)
             !
          end if
          !
          if (n >= 3) then
             !
             ver(3) = ( -g3*ver(2) + 1._ki/4._ki*( -3._ki*g3**2 &
                  - g1**2 + 4._ki*g1*g3 )/(g1-g3)**2 )/(g1-g3)
             vert(3) = ver(3)
             !
          end if
          !
          if (n >= 4) then
             !
             ver(4) = ( -g3*ver(3) + 1._ki/18._ki*(  &
                  + 11._ki*g3**3 - 18._ki*g3**2*g1 - 2._ki*g1**3  &
                  + 9._ki*g1**2*g3 )/(g1-g3)**3 )/(g1-g3)
             vert(4) = ver(4)
             !
          end if
          !
          rest = vert(n)
          !
          hf_carg(1) = real(rest,ki)/plus_grand_glob
          hf_carg(2) = aimag(rest)/plus_grand_glob
          !
       else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function hf (file mod_hf.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the choice rat has been made, it is singular when a1=a3 %d0'
          tab_erreur_par(2)%arg_real=abs(g1-g3)
          call catch_exception(0)
          !
       end if
       !
    end if
    !
  end function hf_carg
  !
  ! variables a1_glob, a3_glob, expo_glob and eps_glob are global in this 
  ! module whereas variables lambda_par,beta_par,alpha_par are given by the 
  ! module parametre
  !****if* src/integrals/three_point/func_hf/eval_numer_hf
  ! NAME
  !
  !  Function eval_numer_hf
  !
  ! USAGE
  !
  !  complex = eval_numer_hf(u)
  !
  ! DESCRIPTION
  !
  !  This is the integrand for the numerical evaluation of hf_rarg
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
  function eval_numer_hf(u)
    !
    real(ki), intent (in) :: u
    complex(ki) :: eval_numer_hf
    !
    real(ki) :: x,y
    complex(ki) :: z,jacob
    !
    x = u
    y = -lambda_par*u**alpha_par*(1._ki-u)**beta_par
    jacob = 1._ki - eps_glob*i_*lambda_par*u**(alpha_par-1._ki)&
         *(1._ki-u)**(beta_par-1._ki)*(alpha_par*(1._ki-u)-beta_par*u)
    z = x + eps_glob*i_*y
    eval_numer_hf = z**(expo_glob-1)*( log(z*a1_glob+(1._ki-z)*a3_glob) &
         + log(plus_grand_glob/mu2_scale_par) )/(z*a1_glob+(1._ki-z)*a3_glob)
    eval_numer_hf = -eval_numer_hf*jacob
    !
  end function eval_numer_hf
  !
  ! variables a1_glob, a3_glob, expo_glob and eps_glob are global in this 
  ! module whereas variables lambda_par,beta_par,alpha_par are given by the 
  ! module parametre
  !****if* src/integrals/three_point/func_hf/eval_numer_hf
  ! NAME
  !
  !  Function eval_numer_hf_c
  !
  ! USAGE
  !
  !  complex = eval_numer_hf_c(u)
  !
  ! DESCRIPTION
  !
  !  This is the integrand for the numerical evaluation of hf_carg
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
  function eval_numer_hf_c(u)
    !
    real(ki), intent (in) :: u
    complex(ki) :: eval_numer_hf_c
    !
    real(ki) :: x,y
    complex(ki) :: z,jacob
    !
    x = u
    y = -lambda_par*u**alpha_par*(1._ki-u)**beta_par
    jacob = 1._ki - eps_glob*i_*lambda_par*u**(alpha_par-1._ki)&
         *(1._ki-u)**(beta_par-1._ki)*(alpha_par*(1._ki-u)-beta_par*u)
    z = x + eps_glob*i_*y
    eval_numer_hf_c = z**(expo_glob-1)*( log(z*a1_glob_c+(1._ki-z)*a3_glob_c) &
         + log(plus_grand_glob/mu2_scale_par) )/(z*a1_glob_c+(1._ki-z)*a3_glob_c)
    eval_numer_hf_c = -eval_numer_hf_c*jacob
    !
  end function eval_numer_hf_c
  !
  ! This function computes the same thing as he
  ! but it returns a complex instead of two dim array
  !****f* src/integrals/three_point/func_hf/hf_c
  ! NAME
  !
  !  Function hf_c
  !
  ! USAGE
  !
  !  complex = hf_c(n,a1,a3)
  !
  ! DESCRIPTION
  !
  ! This function computes the same thing as hf
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
  function hf_c(n,a1,a3)
    !
    integer, intent(in) :: n
    real(ki), intent(in) :: a1,a3
    complex(ki) :: hf_c
    !
    real(ki), dimension(2) :: temp
    !
    temp = hf(n,a1,a3)
    hf_c = cmplx(temp(1),temp(2),ki)
    !
  end function hf_c
  !
  !****f* src/integrals/three_point/func_hf/hf_gen
  ! NAME
  !
  !  Function hf_gen
  !
  ! USAGE
  !
  !  real_dim2 = hf_gen(n,a1,b1,a3,b3)
  !
  ! DESCRIPTION
  !
  !  This function computes:
  !  int^1_0 dy y^n*ln(y*z1+(1-y)*z3)/(y*z1+(1-y)*z3)
  !  where z1 = a1 + i b1 and z3 = a3 + i b3
  !  For n=1, it is equal to: (ln^2(z1)-ln^2(z3))/(z1-z3)
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
  function hf_gen(n,a1,b1,a3,b3)
    !
    integer, intent(in) :: n
    real(ki), intent(in) :: a1,b1,a3,b3
    real(ki), dimension(2) :: hf_gen
    !
    complex(ki) :: rest
    complex(ki) :: abserr
    complex(ki), dimension(4) :: ver,verm
    real(ki) :: lm
    !
    plus_grand_glob = max(abs(a1),abs(a3))
    ! les variables a1_glob, a3_glob, expo_glob et eps_glob
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
       lm = log(plus_grand_glob/mu2_scale_par)
       !
       ver(1) = (z_log2(a1,b1)-z_log2(a3,b3))/(a1-a3)/2._ki
       ver(2) = (-a3*ver(1)-1._ki+(a1*z_log(a1,b1)-a3*z_log(a3,b3))/(a1-a3))&
            /(a1-a3)
       ver(3) = ( -a3*ver(2)+1._ki/4._ki*(-3._ki*a3**2+2._ki*z_log(a3,b3)*a3**2 &
            -a1**2+4._ki*a1*a3+2._ki*z_log(a1,b1)*a1**2-4._ki*z_log(a1,b1) &
            *a1*a3)/(a1-a3)**2 )/(a1-a3)
       ver(4) = ( -a3*ver(3)+1._ki/18._ki*(-6._ki*z_log(a3,b3)*a3**3 &
            +11._ki*a3**3-18._ki*a3**2*a1-18._ki*z_log(a1,b1)*a1**2*a3 &
            +18._ki*z_log(a1,b1)*a1*a3**2-2._ki*a1**3+9._ki*a1**2*a3 &
            +6._ki*z_log(a1,b1)*a1**3)/(a1-a3)**3 )/(a1-a3)
       !
       verm(1) = (z_log(a1,b1)-z_log(a3,b3))/(a1-a3)
       verm(2) = (-a3*ver(1)+1._ki)/(a1-a3)
       verm(3) = (-a3*ver(2)+1._ki/2._ki)/(a1-a3)
       verm(4) = (-a3*ver(3)+1._ki/3._ki )/(a1-a3)
       !
       rest = ver(n)+lm*verm(n)
       !
    else if ( (abs(a1_glob-a3_glob) <= coupure_3p2m) .and. &
         (abs(a1_glob-a3_glob) >= tiny(a1_glob)) ) then
       !
       origine_info_par = "hf_gen"
       num_grand_b_info_par = n
       denom_grand_b_info_par = abs(a1-a3)
       !
       call generic_eval_numer(eval_numer_hf_gen,0._ki,1._ki,tolerance,rest,abserr)
       !
       tab_erreur_par(1)%a_imprimer = .true.
       tab_erreur_par(1)%chaine = 'In function hf_gen:'
       tab_erreur_par(2)%a_imprimer = .true.
       tab_erreur_par(2)%chaine = 'the error returned by adapt_gauss1 is: %z0'
       tab_erreur_par(2)%arg_comp = abserr
       call catch_exception(1)
       !
    else if ( (abs(a1_glob-a3_glob) <= tiny(a1_glob)) ) then
       !
       rest = (z_log(-a1_glob,-1._ki) + log(plus_grand_glob/mu2_scale_par))/a1_glob/real(n,ki)
       !
    end if
    !
    hf_gen(1) = real(rest,ki)/plus_grand_glob
    hf_gen(2) = aimag(rest)/plus_grand_glob
    !
  end function hf_gen
  !
  !****if* src/integrals/three_point/func_hf/eval_numer_hf_gen
  ! NAME
  !
  !  Function eval_numer_hf_gen
  !
  ! USAGE
  !
  !  complex = eval_numer_hf_gen(u)
  !
  ! DESCRIPTION
  !
  !  This is the integrand for the numerical evaluation of hf_gen
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
  function eval_numer_hf_gen(u)
    !
    real(ki), intent (in) :: u
    complex(ki) :: eval_numer_hf_gen
    !
    real(ki) :: x,y
    complex(ki) :: z,jacob
    !
    x = u
    y = -lambda_par*u**alpha_par*(1._ki-u)**beta_par
    jacob = 1._ki - eps_glob*i_*lambda_par*u**(alpha_par-1._ki)&
         *(1._ki-u)**(beta_par-1._ki)*(alpha_par*(1._ki-u)-beta_par*u)
    z = x + eps_glob*i_*y
    eval_numer_hf_gen = z**(expo_glob-1)*( log(z*a1_glob+(1._ki-z)*a3_glob) &
         + log(plus_grand_glob/mu2_scale_par) )&
         /(z*a1_glob+(1._ki-z)*a3_glob)
    eval_numer_hf_gen = eval_numer_hf_gen*jacob
    !
  end function eval_numer_hf_gen
  !
end module func_hf
