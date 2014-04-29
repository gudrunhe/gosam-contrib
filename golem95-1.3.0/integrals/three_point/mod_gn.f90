!
!****h* src/integrals/three_point/func_gn
! NAME
!
!  Module func_gn
!
! USAGE
!
!  use func_gn
!
! DESCRIPTION
!
!  This module contains several functions for the computation of
!  int^1_0 dx x^(n-1)*ln(a*x^2+b*x+c-i*lambda)/(a*x^2+b*x+c-i*lambda) 
!  where a, b and c are real numbers
!
! OUTPUT
!
!  This modules exports three functions:
!  * ge -- a function
!  * gl -- a function
!  * gf -- a function
!
! USES
!
!  * precision_golem (src/module/precision.f90)
!  * numerical_evaluation (src/numerical/mod_numeric.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!         only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
!  * parametre (src/module/parametre.f90)
!  * logarithme (src/module/z_log.f90)
!  * dilogarithme (src/module/zdilog.f90)
!  * constante (src/module/constante.f90) only : i_,un,pi
!
!*****
module func_gn
  !
  use precision_golem
  use numerical_evaluation
  use sortie_erreur, only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
  use parametre
  use logarithme
  use dilogarithme
  use constante, only : i_,un,pi
  implicit none
  !
  private
  real(ki) :: a_glob,b_glob,c_glob
  real(ki) :: lm_glob
  real(ki) :: lambda_glob
  integer :: expo_glob
  logical :: dist_glob
  public :: gf,ge,gl
  !
  contains
    !
    !****f* src/integrals/three_point/func_gn/ge
    ! NAME
    !
    !  Function ge
    !
    ! USAGE
    !
    !  real_dim2 = ge(n,a,b,c,dist)
    !
    ! DESCRIPTION
    !
    !  This function computes:
    !  int^1_0 dx x^(n-1)/(a*x^2+b*x+c-i*lambda) 
    !  where a, b and c are reals
    !  It switches to numerical evaluation if 
    !  (b^2-4*a*c) < coupure_3p1m_2mi
    !  Around the Landau pole, the divergent part is extracted analytically,
    !  only the rest is computed numerically
    !
    ! INPUTS
    !
    !  * n -- an integer, the power of x in the integrand
    !  * a -- a real (type ki), coefficient of x^2
    !  * b -- a real (type ki), coefficient of x^1
    !  * c -- a real (type ki), coefficient of x^0
    !  * dist -- a logical, true if we are close to the real threshold
    !
    ! SIDE EFFECTS
    !
    !  No side effect, the returned value depends on the global variables
    !  rat_or_tot_par, coupure_3p1m_2mi
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
    function ge(n,ax,bx,cx,dist)
      !
      integer, intent(in) :: n
      real(ki), intent(in) :: ax,bx,cx
      logical, intent(in) :: dist
      real(ki), dimension(2) :: ge
      !
      real(ki) :: a,b,c
      complex(ki) :: veri,veri_rat,rest,abserr
      real(ki) :: x1,x2,deltax,denom
      complex(ki) :: cx1,cx2,cdeltax
      real(ki) :: delta
      complex(ki) :: div_part
      real(ki) :: coeff
      complex(ki) :: extra_part,integ1
      real(ki) :: plus_grand
      logical :: not_small_a
      real(ki) :: sb,sa,aa,x2t
      real(ki) :: cut_e = 1.e-3_ki
      !
      plus_grand = max(abs(ax),abs(bx),abs(cx))
      a = ax/plus_grand
      b = bx/plus_grand
      c = cx/plus_grand
      !
      expo_glob = n
      sb = sign(un,b)
      delta = b*b-4._ki*a*c
      a_glob = a
      b_glob = b
      c_glob = c
      lambda_glob = lambda_par
      dist_glob = dist
      !
      origine_info_par = "ge"
      num_grand_b_info_par = n
      denom_grand_b_info_par = abs(delta)
      !
      if (dist) then
        !
        if (delta > 0._ki) then
          !
          if (a > 0._ki) then
            !
            div_part = 2._ki*i_*pi/sqrt(delta)
            !
          else if (a < 0._ki) then
            !
            div_part = 2._ki*i_*pi/sqrt(delta)
            !
          end if
          !
        else if (delta <= 0._ki) then
          !
          if (a > 0._ki) then
            !
            div_part = 2._ki*pi/sqrt(-delta)
            !
          else if (a < 0._ki) then
            !
            div_part = -2._ki*pi/sqrt(-delta)
            !
          end if
          !
        end if
        !
        integ1 = 1._ki/a*(z_log(a+b+c,-1._ki)-z_log(c,-1._ki))
        !
        select case(n)
          !
          case(1)
            !
            coeff = 1._ki
            extra_part = 0._ki
            !
          case(2)
            !
            coeff = -b/(2._ki*a)
            extra_part = 0.5_ki*integ1
            !
          case(3)
            !
            coeff = (b**2-2._ki*a*c)/(2._ki*a**2)
            extra_part = (1._ki-b/2._ki*integ1)/a
            !
          case(4)
            !
            coeff = (3._ki*a*b*c-b**3)/(2._ki*a**3)
            extra_part = ( a/2._ki-b + (b*b-c*a)/2._ki*integ1 )/a/a
            !
        end select
        !
      else
        !
        div_part = 0._ki
        extra_part = 0._ki
        integ1 = 0._ki
        !
      end if
      !
      if (delta >= 0._ki) then
        !
        aa = abs(a)
        not_small_a = aa > cut_e
        !
        if (a == 0._ki) then
          !
          sa = 1._ki ! arbitrary value
          !
        else
          !
          sa = sign(un,a)
          !
        end if
        !
        x2t = -(b+sb*sqrt(delta))/2._ki
        !
        if (not_small_a) then
          !
          x1 = (-b + sb*sqrt(delta))/(2._ki*a)
          x2 = (-b - sb*sqrt(delta))/(2._ki*a)
          denom = (x1-x2)*a
          !
        else ! no need of x2 and deltax in this case
          !
          x1 = c/x2t
          !x2 = x2t/a
          denom = sb*sqrt(delta)
          !
        end if
        !
        !
        if (sqrt(abs(delta)) >  coupure_3p1m_2mi) then
          !
          if (n == 1) then
            !
            if (not_small_a) then
              !
              veri_rat = 0._ki
              veri =  ( z_log((x1-1._ki)/x1,sb) - z_log((x2-1._ki)/x2,-sb) ) &
                      /denom
              !
            else
              !
              veri_rat = 0._ki ! dummy value
              veri = ( z_log((x1-1._ki)/x1,sb) - z_log((x2t-a)/x2t,-sb) )/denom
              !
            end if
            !
          else if (n == 2) then
            !
            if (not_small_a) then
              !
              veri_rat = 0._ki
              veri =  ( x1*z_log((x1-1._ki)/x1,sb) - x2*z_log((x2-1._ki)/x2,-sb) ) &
                      /denom
              !
            else
              !
              veri_rat = 0._ki ! dummy value
              veri = ( x1*z_log((x1-1._ki)/x1,sb) - q(1,a/x2t,-sb) )/denom
              !
            end if
            !
          else if (n == 3) then
            !
            if (not_small_a) then
              !
              veri_rat =  1._ki/a
              veri = ( x1**2*z_log((x1-1._ki)/x1,sb) - x2**2*z_log((x2-1._ki)/x2,-sb) ) &
                      /denom
              !
            else
              !
              veri_rat = 0._ki ! dummy value
              veri = ( x1**2*z_log((x1-1._ki)/x1,sb)+x1 - q(2,a/x2t,-sb) )/denom
              !
            end if
            !
          else if (n == 4) then
            !
            if (not_small_a) then
              !
              veri_rat =  ( 1._ki/2._ki + x1 + x2 )/a
              veri = ( x1**3*z_log((x1-1._ki)/x1,sb) - x2**3*z_log((x2-1._ki)/x2,-sb) ) &
                      /denom
              !
            else
              !
              veri_rat = 0._ki ! dummy value
              veri = ( x1**3*z_log((x1-1._ki)/x1,sb)+x1**2+x1*0.5_ki - q(3,a/x2t,-sb) )/denom
              !
            end if
            !
          else
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In function ge (file mod_gn.f90)'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = 'n should be 1,2,3,4 but is %d0'
            tab_erreur_par(3)%arg_int = n
            call catch_exception(0)
            !
          end if
          !
          if ( rat_or_tot_par%tot_selected  ) then
            !
            rest = veri + veri_rat
            !
          else if ( ( rat_or_tot_par%rat_selected ) .and. not_small_a ) then
            !
            rest =  veri_rat
            !
          else !if ( rat_or_tot_par%rat_selected  ) then
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In function ge (file mod_gn.f90)'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = 'the value of a is too small to extract the rational part : a= %f0'
            tab_erreur_par(3)%arg_real = a
            call catch_exception(1)
            !
          end if
          !
        else if ( (sqrt(abs(delta)) <=  coupure_3p1m_2mi) .and. &
           & (rat_or_tot_par%tot_selected) ) then
          !
          call generic_eval_numer(eval_numer_ge,0._ki,1._ki,tolerance,rest,abserr)
          !
          rest = rest + coeff*div_part + extra_part
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function ge (file mod_gn.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the choice rat has been made, it is'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'singular when abs(delta) is too small: %f0'
          tab_erreur_par(3)%arg_real = abs(delta)
          call catch_exception(0)
          !
        end if
        !
      else
        !
        cx1 = -b/(2._ki*a) + i_*sqrt(-delta)/(2._ki*abs(a))
        cx2 = -b/(2._ki*a) - i_*sqrt(-delta)/(2._ki*abs(a))
        cdeltax = cx1-cx2
        !
        if (sqrt(abs(delta)) >  coupure_3p1m_2mi) then
          !
          if (n == 1) then
            !
            veri_rat = 0._ki
            veri =  ( log((cx1-1._ki)/cx1) - log((cx2-1._ki)/cx2) ) &
                    /cdeltax/a
            !
          else if (n == 2) then
            !
            veri_rat = 0._ki
            veri =  ( cx1*log((cx1-1._ki)/cx1) - cx2*log((cx2-1._ki)/cx2) ) &
                    /cdeltax/a
            !
          else if (n == 3) then
            !
            veri_rat = 1._ki/a
            veri = ( cx1**2*log((cx1-1._ki)/cx1) - cx2**2*log((cx2-1._ki)/cx2) ) &
                    /cdeltax/a
            !
          else if (n == 4) then
            !
            veri_rat = ( 1._ki/2._ki + cx1 + cx2 )/a
            veri = ( cx1**3*log((cx1-1._ki)/cx1) - cx2**3*log((cx2-1._ki)/cx2) ) &
                    /cdeltax/a
            !
          else
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In function ge (file mod_gn.f90)'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = 'n should be 1,2,3,4 but is %d0'
            tab_erreur_par(3)%arg_int = n
            call catch_exception(0)
            !
          end if
          !
          if ( rat_or_tot_par%tot_selected  ) then
            !
            rest = veri + veri_rat
            !
          else !if ( rat_or_tot_par%rat_selected ) then
            !
            rest =  veri_rat
            !
          end if
          !
        else if ( (sqrt(abs(delta)) <=  coupure_3p1m_2mi) .and. &
           & (rat_or_tot_par%tot_selected) ) then
          !
          call generic_eval_numer(eval_numer_ge,0._ki,1._ki,tolerance,rest,abserr)
          !
          rest = rest + coeff*div_part + extra_part
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function ge (file mod_gn.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the choice rat has been made, it is'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'singular when abs(delta) is too small: %f0'
          tab_erreur_par(3)%arg_real = abs(delta)
          call catch_exception(0)
          !
        end if
        !
      end if
      !
      ge(1) = real(rest,ki)/plus_grand
      ge(2) = aimag(rest)/plus_grand
      !
      !
    end function ge
    !
    !****f* src/integrals/three_point/func_gn/gl
    ! NAME
    !
    !  Function gl
    !
    ! USAGE
    !
    !  real_dim2 = gl(n,a,b,c)
    !
    ! DESCRIPTION
    !
    !  This function computes:
    !  int^1_0 dx x^(n-1)*ln(a*x^2+b*x+c-i*lambda) 
    !  where a, b and c are reals
    !  here, no need to switch to numerical evaluation
    !  no numerical problems when (b^2-4*a*c) = 0
    !
    ! INPUTS
    !
    !  * n -- an integer, the power of x in the integrand
    !  * a -- a real (type ki), coefficient of x^2
    !  * b -- a real (type ki), coefficient of x^1
    !  * c -- a real (type ki), coefficient of x^0
    !
    ! SIDE EFFECTS
    !
    !  No side effect, the returned value depends on the global variable
    !  rat_or_tot_par
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
    function gl(n,ax,bx,cx)
      !
      integer, intent(in) :: n
      real(ki), intent(in) :: ax,bx,cx
      real(ki), dimension(2) :: gl
      !
      real(ki) :: a,b,c
      complex(ki) :: veri,veri_rat,rest
      real(ki) :: x1,x2
      complex(ki) :: cx1,cx2
      real(ki) :: delta,plus_grand,lm
      logical :: not_small_a
      real(ki) :: sb,sa,aa,x2t
      real(ki) :: cut_l = 1.e-5_ki
      !
      plus_grand = max(abs(ax),abs(bx),abs(cx))
      lm = log(plus_grand)
      a = ax/plus_grand
      b = bx/plus_grand
      c = cx/plus_grand
      !
      sb = sign(un,b)
      delta = b*b-4._ki*a*c
      !
      if (delta >= 0._ki) then
        !
        aa = abs(a)
        not_small_a = aa > cut_l
        !
        if (a == 0._ki) then
          !
          sa = 1._ki ! arbitrary value
          !
        else
          !
          sa = sign(un,a)
          !
        end if
        !
        x2t = -(b+sb*sqrt(delta))/2._ki
        !
        if (not_small_a) then
          !
          x1 = (-b + sb*sqrt(delta))/(2._ki*a)
          x2 = (-b - sb*sqrt(delta))/(2._ki*a)
          !
        else
          !
          x1 = c/x2t
          !x2 = x2t/a
          !
        end if
        !
        if (n == 1) then
          !
          if (not_small_a) then
            !
            veri_rat = -2._ki
            !veri =  ( z_log(a,-1._ki) + (1._ki-x1)*z_log(1._ki-x1,-1._ki) + x1*z_log(-x1,-1._ki) &
                    !+ (1._ki-x2)*z_log(1._ki-x2,1._ki) + x2*z_log(-x2,1._ki) )
            veri =  ( z_log(a,-1._ki) + (1._ki-x1)*z_log(1._ki-x1,-sb) + x1*z_log(-x1,-sb) &
                    + (1._ki-x2)*z_log(1._ki-x2,sb) + x2*z_log(-x2,sb) )
            !
          else
            !
            veri_rat = 0._ki
            !veri = ( -2._ki + (1._ki-x1)*z_log(1._ki-x1,-sb) + x1*z_log(-x1,-sb) &
                    !+ z_log(sa,-1._ki) + z_log(sa*(a-x2t),sb) - q(1,a/x2t,-sb) )
            veri = ( -2._ki - x1*( z_log(1._ki-x1,-sb) - z_log(-x1,-sb) ) &
                    + z_log(a+b+c,-1._ki) - q(1,a/x2t,-sb) )
            !
          end if
          !
        else if (n == 2) then
          !
          if (not_small_a) then
            !
            veri_rat = -(1._ki + x1 + x2)/2._ki
            veri =  ( z_log(a,-1._ki) + (1._ki-x1**2)*z_log(1._ki-x1,-sb) + x1**2*z_log(-x1,-sb) &
                    + (1._ki-x2**2)*z_log(1._ki-x2,sb) + x2**2*z_log(-x2,sb) )/2._ki
            !
          else
            !
            veri_rat = 0._ki
            !veri = ( -(1._ki+x1) + (1._ki-x1**2)*z_log(1._ki-x1,-sb) + x1**2*z_log(-x1,-sb) &
                    !+ z_log(sa,-1._ki) + z_log(sa*(a-x2t),sb) - q(2,a/x2t,-sb) )/2._ki
            veri = ( -(1._ki+x1) - x1**2*( z_log(1._ki-x1,-sb) - z_log(-x1,-sb) ) &
                    + z_log(a+b+c,-1._ki) - q(2,a/x2t,-sb) )/2._ki
            !
          end if
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function gl (file mod_gn.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'value of n not implemented'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'n= %d0'
          tab_erreur_par(3)%arg_int = n
          call catch_exception(0)
          !
        end if
        !
        veri =  veri + veri_rat
        rest = veri
        !
      else
        !
        !
        cx1 = -b/(2._ki*a) + i_*sqrt(-delta)/(2._ki*abs(a))
        cx2 = -b/(2._ki*a) - i_*sqrt(-delta)/(2._ki*abs(a))
        !
        if (n == 1) then
          !
          veri_rat = -2._ki
          veri =  ( z_log(a,-1._ki) + (1._ki-cx1)*log(1._ki-cx1) + cx1*log(-cx1) &
                  + (1._ki-cx2)*log(1._ki-cx2) + cx2*log(-cx2) )
          !
        else if (n == 2) then
          !
          veri_rat = -(1._ki + cx1 + cx2)/2._ki
          veri =  ( z_log(a,-1._ki) + (1._ki-cx1**2)*log(1._ki-cx1) + cx1**2*log(-cx1) &
                  + (1._ki-cx2**2)*log(1._ki-cx2) + cx2**2*log(-cx2) )/2._ki
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function gl (file mod_gn.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'value of n not implemented'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'n= %d0'
          tab_erreur_par(3)%arg_int = n
          call catch_exception(0)
          !
        end if
        !
        veri =  veri + veri_rat
        rest = veri
        !
      end if
      !
      gl(1) = real(rest,ki) + lm/real(n,ki)
      gl(2) = aimag(rest)
      !
    end function gl
    !
    !****f* src/integrals/three_point/func_gn/gf
    ! NAME
    !
    !  Function gf
    !
    ! USAGE
    !
    !  real_dim2 = gf(n,a,b,c,dist)
    !
    ! DESCRIPTION
    !
    !  This function computes:
    !  int^1_0 dx x^(n-1)*ln(a*x^2+b*x+c-i*lambda)/(a*x^2+b*x+c-i*lambda) 
    !  where a, b and c are reals
    !  It switches to numerical evaluation if 
    !  (b^2-4*a*c) < coupure_3p1m_2mi
    !  Around the Landau pole, the divergent part is extracted analytically,
    !  only the rest is computed numerically
    !
    ! INPUTS
    !
    !  * n -- an integer, the power of x in the integrand
    !  * a -- a real (type ki), coefficient of x^2
    !  * b -- a real (type ki), coefficient of x^1
    !  * c -- a real (type ki), coefficient of x^0
    !  * dist -- a logical, true if we are close to the real threshold
    !
    ! SIDE EFFECTS
    !
    !  No side effect, the returned value depends on the global variables
    !  rat_or_tot_par, coupure_3p1m_2mi
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
    function gf(n,ax,bx,cx,dist)
      !
      integer, intent(in) :: n
      real(ki), intent(in) :: ax,bx,cx
      logical, intent(in) :: dist
      real(ki), dimension(2) :: gf
      !
      real(ki) :: a,b,c
      complex(ki) :: veri,verim,rest,abserr
      real(ki) :: x1,x2,deltax,denom
      complex(ki) :: cx1,cx2,cdeltax
      real(ki) :: delta,plus_grand
      complex(ki) :: div_part
      logical :: not_small_a
      real(ki) :: sb,sa,aa,x2t
      real(ki) :: cut_f = 1.e-6_ki
      !
      plus_grand = max(abs(ax),abs(bx),abs(cx))
      lm_glob = log(plus_grand)
      a = ax/plus_grand
      b = bx/plus_grand
      c = cx/plus_grand
      !
      expo_glob = n
      sb = sign(un,b)
      delta = b*b-4._ki*a*c
      a_glob = a
      b_glob = b
      c_glob = c
      lambda_glob = lambda_par
      dist_glob = dist
      !
      origine_info_par = "gf"
      num_grand_b_info_par = n
      denom_grand_b_info_par = abs(delta)
      !
      ! divergent part for Landau singularities
      if (dist) then
        !
        if (delta > 0._ki) then
          !
          if (a > 0._ki) then
            !
            div_part = 2._ki*i_*pi/sqrt(delta)*(log(delta/a) - i_*pi + lm_glob )
            !
          else if (a < 0._ki) then
            !
            div_part = 2._ki*i_*pi/sqrt(delta)*(log(-delta/a) + lm_glob)
            !
          end if
          !
        else if (delta <= 0._ki) then
          !
          if (a > 0._ki) then
            !
            div_part = 2._ki*pi/sqrt(-delta)*(log(-delta/a) + lm_glob)
            !
          else if (a < 0._ki) then
            !
            div_part = -2._ki*pi/sqrt(-delta)*(log(delta/a) - i_*pi + lm_glob)
            !
          end if
          !
        end if
        !
      else
        !
        div_part = 0._ki
        !
      end if
      !
      if (delta >= 0._ki) then
        !
        aa = abs(a)
        not_small_a = aa > cut_f
        !
        if (a == 0._ki) then
          !
          sa = 1._ki ! arbitrary value
          !
        else
          !
          sa = sign(un,a)
          !
        end if
        !
        x2t = -(b+sb*sqrt(delta))/2._ki
        !
        if (not_small_a) then
          ! we do not care which is x1 and which is x2, in this way
          ! whatever th value of b is, x2 is the root which is sent 
          ! to infinity when a --> 0
          x1 = (-b + sb*sqrt(delta))/(2._ki*a)
          x2 = (-b - sb*sqrt(delta))/(2._ki*a)
          deltax = (x1-x2)
          denom = deltax*a
        else ! no need of x2 and deltax in this case
          !
          x1 = c/x2t
          !x2 = x2t/a
          denom = sb*sqrt(delta)
          !
        end if
        !
        if (sqrt(abs(delta)) >  coupure_3p1m_2mi) then
          !
          if (not_small_a) then
            !
          !veri =  ( -zdilog((1._ki-x2)/(1._ki-x1),sign(1._ki,2._ki-x1-x2)) + zdilog(x2/x1,sign(1._ki,-x1-x2)) &
                  !+ zdilog((1._ki-x1)/(1._ki-x2),sign(1._ki,-2._ki+x1+x2)) - zdilog(x1/x2,sign(1._ki,x1+x2)) &
                  !+ ( 2._ki*z_log(deltax,1._ki) - i_*pi + z_log(a,-1._ki) ) &
                  !*( z_log((x1-1._ki)/x1,1._ki) - z_log((x2-1._ki)/x2,-1._ki) ) &
                  !)/deltax/a
          !verim =  ( z_log((x1-1._ki)/x1,1._ki) - z_log((x2-1._ki)/x2,-1._ki) ) &
                  !/deltax/a
          veri =  ( -zdilog((1._ki-x2)/(1._ki-x1),sign(sb,2._ki-x1-x2)) + zdilog(x2/x1,sign(sb,-x1-x2)) &
                  + zdilog((1._ki-x1)/(1._ki-x2),sign(sb,-2._ki+x1+x2)) - zdilog(x1/x2,sign(sb,x1+x2)) &
                  + ( 2._ki*z_log(deltax,sb) - i_*pi + z_log(a,-1._ki) ) &
                  *( z_log((x1-1._ki)/x1,sb) - z_log((x2-1._ki)/x2,-sb) ) &
                  )/denom
          verim =  ( z_log((x1-1._ki)/x1,sb) - z_log((x2-1._ki)/x2,-sb) ) &
                  /denom
            !
          else
            !
            ! in order that the small imaginary part get the right sign, we must have x2t > a
            !
            veri = ( 2._ki*( zdilog(a*(x2t-c)/(x2t*(a-x2t)),-sa) - zdilog(a*c/x2t/x2t,-sa) ) &
                    + ( z_log((x1-1._ki)/x1,sb) - z_log((x2t-a)/x2t,-sb) )* &
                    ( 0.5_ki*(z_log(1._ki-x1,-sb) + z_log(x1,sb) - z_log(sa*x2t-aa,-sb) - z_log(-sa*x2t,sb)) &
                     + 2._ki*z_log(sb*sa*sqrt(delta),sb) + z_log(sa,-1._ki)-i_*pi) &
                  )/denom
            verim = ( z_log((x1-1._ki)/x1,sb) - z_log((x2t-a)/x2t,-sb) )/denom
            !
          end if
          !
          rest = ( veri + verim*lm_glob )/plus_grand
          !
        else
          !
          call generic_eval_numer(eval_numer_gf,0._ki,1._ki,tolerance,rest,abserr)
          !
          rest = (rest + div_part)/plus_grand
          !
        end if
        !
      else
        !
        cx1 = -b/(2._ki*a) + i_*sqrt(-delta)/(2._ki*abs(a))
        cx2 = -b/(2._ki*a) - i_*sqrt(-delta)/(2._ki*abs(a))
        cdeltax = cx1-cx2
        !
        if (sqrt(abs(delta)) >  coupure_3p1m_2mi) then
          !
          veri =  ( -cdilog((1._ki-cx2)/(1._ki-cx1)) + cdilog(cx2/cx1) &
                  + cdilog((1._ki-cx1)/(1._ki-cx2)) - cdilog(cx1/cx2) &
                  + ( 2._ki*log(cdeltax) - i_*pi + z_log(a,-1._ki) ) &
                  *( log((cx1-1._ki)/cx1) - log((cx2-1._ki)/cx2) ) &
                  )/cdeltax/a
          verim =  ( log((cx1-1._ki)/cx1) - log((cx2-1._ki)/cx2) ) &
                  /cdeltax/a
          !
          rest = ( veri + verim*lm_glob )/plus_grand
          !
        else
          !
          call generic_eval_numer(eval_numer_gf,0._ki,1._ki,tolerance,rest,abserr)
          !
          rest = (  rest + div_part )/plus_grand
          !
        end if
        !
      end if
      !
      gf(1) = real(rest,ki)
      gf(2) = aimag(rest)
      !
      !
    end function gf
    !
    !****if* src/integrals/three_point/func_gn/eval_numer_ge
    ! NAME
    !
    !  Function eval_numer_ge
    !
    ! USAGE
    !
    !  complex = eval_numer_ge(u)
    !
    ! DESCRIPTION
    !
    !  This is the integrand for the numerical evaluation for ge.
    !  Depending if 
    !  part 1/( (z-x_1)*(z-x_2) )
    !
    ! INPUTS
    !
    !  * u -- a real (type ki), the integration variable
    !
    ! SIDE EFFECTS
    !
    !  No side effect. The variables of type xxx_glob
    !  are global in this module
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
    function eval_numer_ge(u)
      !
      real(ki), intent (in) :: u
      complex(ki) :: eval_numer_ge
      !
      real(ki) :: x,y
      complex(ki) :: z,jacob
      real(ki) :: sigma,delta
      complex(ki) :: integ2
      !
      if (dist_glob) then
        !
        sigma = 1._ki+b_glob/a_glob/2._ki
        !
      else
        !
        sigma = -b_glob/a_glob/2._ki
        !
      end if
      !
      delta = b_glob*b_glob-4._ki*a_glob*c_glob
      !
      x = u
      !
      if ( (sigma <= 1._ki) .and. (sigma >= 0._ki) ) then
        !
        y = lambda_glob*sign(un,a_glob)*u*(u-1._ki)*(u-sigma)
        z = x + i_*y
        jacob = 1._ki + i_*lambda_glob*sign(un,a_glob)*( (u-1._ki)*(u-sigma) + u*(u-1._ki) + u*(u-sigma) )
        !
      else
        !
        y = lambda_glob*sign(un,a_glob*sigma)*u*(u-1._ki)
        z = x - i_*y
        jacob = 1._ki - i_*lambda_glob*sign(un,a_glob*sigma)*( (u-1._ki) + u )
        !
      end if
      !
      if (dist_glob) then
        !
        integ2 = -b_glob/(2._ki*a_glob+b_glob) &
          &/(-delta*a_glob/(2._ki*a_glob+b_glob)**2*z*z+delta/(2._ki*a_glob+b_glob)*z+c_glob)
        !
        select case(expo_glob)
          !
          case(1)
            !
            eval_numer_ge = -integ2
            !
          case(2)
            !
            eval_numer_ge = b_glob/2._ki/a_glob*integ2
            !
          case(3)
            !
            eval_numer_ge = (c_glob-b_glob*b_glob/2._ki/a_glob)*integ2/a_glob 
            !
          case(4)
            !
            eval_numer_ge = -b_glob*(3._ki*a_glob*c_glob-b_glob*b_glob)/2._ki/a_glob**3&
              &*integ2
            !
          end select
          !
      else
        !
        eval_numer_ge = z**(expo_glob-1)/(a_glob*z*z+b_glob*z+c_glob)
        !
      end if
      !
      eval_numer_ge = eval_numer_ge*jacob
      !
    end function eval_numer_ge
    !
    !****if* src/integrals/three_point/func_gn/eval_numer_gf
    ! NAME
    !
    !  Function eval_numer_gf
    !
    ! USAGE
    !
    !  complex = eval_numer_gf(u)
    !
    ! DESCRIPTION
    !
    !  This is the integrand for the numerical evaluation of gf,
    !  part ln(z-x_1)/( (z-x_1)*(z-x_2) )
    !
    ! INPUTS
    !
    !  * u -- a real (type ki), the integration variable
    !
    ! SIDE EFFECTS
    !
    !  No side effect. The variables of type xxx_glob
    !  are global in this module
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
    function eval_numer_gf(u)
      !
      real(ki), intent (in) :: u
      complex(ki) :: eval_numer_gf
      !
      real(ki) :: x,y
      complex(ki) :: z,jacob
      real(ki) :: sigma,delta,a_coeff
      !
      !
      delta = b_glob*b_glob-4._ki*a_glob*c_glob
      !
      if (dist_glob) then
        !
        sigma = 1._ki+b_glob/a_glob/2._ki
        a_coeff = -delta
        !
      else
        !
        sigma = -b_glob/a_glob/2._ki
        a_coeff = a_glob
        !
      end if
      !
      x = u
      !
      if ( (sigma <= 1._ki) .and. (sigma >= 0._ki) ) then
        !
        y = lambda_glob*sign(un,a_coeff)*u*(u-1._ki)*(u-sigma)
        z = x + i_*y
        jacob = 1._ki + i_*lambda_glob*sign(un,a_coeff)*( (u-1._ki)*(u-sigma) + u*(u-1._ki) + u*(u-sigma) )
        !
      else
        !
        y = lambda_glob*sign(un,a_coeff*sigma)*u*(u-1._ki)
        z = x + i_*y
        jacob = 1._ki + i_*lambda_glob*sign(un,a_coeff*sigma)*( (u-1._ki) + u )
        !
      end if
      !
      if (dist_glob) then
        !
        eval_numer_gf = -b_glob/(2._ki*a_glob+b_glob) &
          &*log( (1._ki-2._ki*a_glob/(2._ki*a_glob+b_glob)*z)**2 /&
          & (-delta*a_glob/(2._ki*a_glob+b_glob)**2*z*z+delta/(2._ki*a_glob+b_glob)*z+c_glob) )&
          &/(-delta*a_glob/(2._ki*a_glob+b_glob)**2*z*z+delta/(2._ki*a_glob+b_glob)*z+c_glob)
        eval_numer_gf = eval_numer_gf  + lm_glob*( b_glob/(2._ki*a_glob+b_glob) &
          &/(-delta*a_glob/(2._ki*a_glob+b_glob)**2*z*z+delta/(2._ki*a_glob+b_glob)*z+c_glob) )
        !
      else
        !
        eval_numer_gf = (log(a_glob*z*z+b_glob*z+c_glob)+lm_glob)/(a_glob*z*z+b_glob*z+c_glob)
        !
      end if
      !
      eval_numer_gf = eval_numer_gf*jacob
      !
    end function eval_numer_gf
    !
end module func_gn

