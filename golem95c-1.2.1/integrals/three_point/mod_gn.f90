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
!  This modules exports four functions:
!  * ge -- a function
!  * ge_c -- a function
!  * gf -- a function
!  * gf_c -- a function
!
! USES
!
!  * precision_golem (src/module/precision.f90)
!  * numerical_evaluation (src/numerical/mod_numeric.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90) only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
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
  real(ki) :: a_glob,b_glob,c_glob,eps_glob
  complex(ki) :: x1_glob,x2_glob
  real(ki) :: plus_grand_glob
  real(ki) :: lambda_glob,alpha_glob,beta_glob
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
    function ge(n,a,b,c,dist)
    !~ function ge(n,a,b,c)
      !
      integer, intent(in) :: n
      real(ki), intent(in) :: a,b,c
      logical, intent(in) :: dist
      real(ki), dimension(2) :: ge
      !
      complex(ki) :: veri,veri_rat,rest,abserr
      complex(ki) :: residue
      complex(ki) :: extra_imag
      real(ki) :: pole1
      real(ki) :: x1,x2,deltax
      complex(ki) :: cx1,cx2,cdeltax
      real(ki) :: delta
      logical :: l1
      complex(ki) :: div_part
      real(ki) :: coeff
      !
      plus_grand_glob = max(abs(a),abs(b),abs(c))
      expo_glob = n
      delta = b*b-4._ki*a*c
      a_glob = a
      b_glob = b
      c_glob = c
      lambda_glob = lambda_par
      alpha_glob = alpha_par
      beta_glob = beta_par
      dist_glob = dist
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
      else
        !
        div_part = 0._ki
        !
      end if
      !
      select case(n)
        !
        case(1)
          !
          coeff = 1._ki
          !
        case(2)
          !
          coeff = -b/(2._ki*a)
          !
        case(3)
          !
          coeff = (b**2-2._ki*a*c)/(2._ki*a**2)
          !
        case(4)
          !
          coeff = (3._ki*a*b*c-b**3)/(2._ki*a**3)
          !
      end select
      !
      if (delta >= 0._ki) then
        !
        x1 = (-b + sqrt(delta))/(2._ki*a)
        x2 = (-b - sqrt(delta))/(2._ki*a)
        x1_glob = cmplx(x1,0._ki)
        x2_glob = cmplx(x2,0._ki)
        deltax = x1-x2
        !
        if (sqrt(abs(delta)) >  coupure_3p1m_2mi) then
        !~ if (abs(delta) >  0._ki) then
          !
          if (n == 1) then
            !
            veri_rat = 0._ki
            veri =  ( z_log((x1-1._ki)/x1,1._ki) - z_log((x2-1._ki)/x2,-1._ki) ) &
                    /deltax/a
            !
          else if (n == 2) then
            !
            veri_rat = 0._ki
            veri =  ( x1*z_log((x1-1._ki)/x1,1._ki) - x2*z_log((x2-1._ki)/x2,-1._ki) ) &
                    /deltax/a
            !
          else if (n == 3) then
            !
            veri_rat =  1._ki/a
            veri = ( x1**2*z_log((x1-1._ki)/x1,1._ki) - x2**2*z_log((x2-1._ki)/x2,-1._ki) ) &
                    /deltax/a
            !
          else if (n == 4) then
            !
            veri_rat =  ( 1._ki/2._ki + x1 + x2 )/a
            veri = ( x1**3*z_log((x1-1._ki)/x1,1._ki) - x2**3*z_log((x2-1._ki)/x2,-1._ki) ) &
                    /deltax/a
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
            ! to please the compiler
            stop
            !
          end if
          !
          if ( rat_or_tot_par%tot_selected  ) then
            !
            rest = veri + veri_rat
            !
          else !if ( rat_or_tot_par%rat_selected  ) then
            !
            rest =  veri_rat
            !
          end if
          !
        else if ( (sqrt(abs(delta)) <=  coupure_3p1m_2mi) .and. &
           & (rat_or_tot_par%tot_selected) ) then
          !
          !~ pole1 = x2
          !~ eps_glob = 1._ki
          !
          !~ call generic_eval_numer(eval_numer_ge,0._ki,1._ki,1.0E-8_ki,rest,abserr)
          call generic_eval_numer(eval_numer_ge,0._ki,1._ki,tolerance,rest,abserr)
          !
          !~ if ( (pole1 >= 0._ki) .and. (pole1 <= 1._ki) ) then
            !~ !
            !~ residue = x2**(n-1)/(x2-x1)/a_glob
            !~ extra_imag = -2._ki*i_*pi*residue
            !~ !
          !~ else
            !~ !
            !~ extra_imag = 0._ki
            !~ !
          !~ end if
          !~ !
          !~ rest = rest + extra_imag
          !
          !~ write(*,*) 'exact result is:',coeff*div_part
          !~ write(*,*) 'err numer result is:',abserr,rest
          !~ write(*,*) 'n is:',n
          !~ select case(n)
          !~ case(1)
          !~ write(*,*) 'analytical result is:', &
            !~ ( z_log((x1-1._ki)/x1,1._ki) - z_log((x2-1._ki)/x2,-1._ki) )/deltax/a
          !~ case(2)
          !~ write(*,*) 'analytical result is:', &
            !~ ( x1*z_log((x1-1._ki)/x1,1._ki) - x2*z_log((x2-1._ki)/x2,-1._ki) )/deltax/a
          !~ case(3)
          !~ write(*,*) 'analytical result is:', &
            !~ ( x1**2*z_log((x1-1._ki)/x1,1._ki) - x2**2*z_log((x2-1._ki)/x2,-1._ki) )/deltax/a + 1._ki/a
          !~ case(4)
          !~ write(*,*) 'analytical result is:', &
            !~ ( x1**3*z_log((x1-1._ki)/x1,1._ki) - x2**3*z_log((x2-1._ki)/x2,-1._ki) )/deltax/a + ( 1._ki/2._ki + x1 + x2 )/a
            !~ end select
          !~ if (dist) then
          rest = rest + coeff*div_part
          !~ write(*,*) 'numer result is:',rest
          !~ end if
        !else if ( (abs(delta) <=  coupure_3p1m_2mi) .and. &
        !   & (rat_or_tot_par%rat_selected) ) then
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
          ! to please the compiler
          stop
          !
        end if
        !
      else
        !
        origine_info_par = "ge"
        num_grand_b_info_par = n
        !! Thomas Reiter, 01/March/2011:
        !! original code:
        ! denom_grand_b_info_par = abs(deltax)
        !! cannot be right, deltax was never initialized
        !! probably: delta instead of deltax.
        denom_grand_b_info_par = abs(delta)
        !
        !
        cx1 = -b/(2._ki*a) + i_*sqrt(-delta)/(2._ki*abs(a))
        cx2 = -b/(2._ki*a) - i_*sqrt(-delta)/(2._ki*abs(a))
        x1_glob = cx1
        x2_glob = cx2
        cdeltax = cx1-cx2
        !
        if (sqrt(abs(delta)) >  coupure_3p1m_2mi) then
        !~ if (abs(delta) >  0._ki) then
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
            ! to please the compiler
            stop
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
          !~ eps_glob = 1._ki
          !~ !
          !~ call inside_contour(cx2,l1)
          !
          !~ call generic_eval_numer(eval_numer_ge,0._ki,1._ki,1.0E-8_ki,rest,abserr)
          call generic_eval_numer(eval_numer_ge,0._ki,1._ki,tolerance,rest,abserr)
          !
          !~ if ( l1) then
            !~ !
            !~ residue = cx2**(n-1)/(cx2-cx1)/a_glob
            !~ extra_imag = -2._ki*i_*pi*residue
            !~ !
          !~ else
            !~ !
            !~ extra_imag = 0._ki
            !~ !
          !~ end if
          !~ !
          !~ rest = rest + extra_imag
          !
        !else if ( (abs(delta) <=  coupure_3p1m_2mi) .and. &
        !   & (rat_or_tot_par%rat_selected) ) then
          !~ write(*,*) 'exact result is:',coeff*div_part
          !~ write(*,*) 'err numer result is:',abserr,rest
          !~ write(*,*) 'n is:',n
          !~ select case(n)
          !~ case(1)
          !~ write(*,*) 'analytical result is:', &
          !~ ( log((cx1-1._ki)/cx1) - log((cx2-1._ki)/cx2) )/cdeltax/a
          !~ case(2)
          !~ write(*,*) 'analytical result is:', &
          !~ ( cx1*log((cx1-1._ki)/cx1) - cx2*log((cx2-1._ki)/cx2) )/cdeltax/a
          !~ case(3)
          !~ write(*,*) 'analytical result is:', &
          !~ ( cx1**2*log((cx1-1._ki)/cx1) - cx2**2*log((cx2-1._ki)/cx2) )/cdeltax/a + 1._ki/a
          !~ case(4)
          !~ write(*,*) 'analytical result is:', &
          !~ ( cx1**3*log((cx1-1._ki)/cx1) - cx2**3*log((cx2-1._ki)/cx2) )/cdeltax/a +  ( 1._ki/2._ki + cx1 + cx2 )/a
          !~ end select
          rest = rest + coeff*div_part
          !~ write(*,*) 'numer result is:',rest
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
          ! to please the compiler
          stop
          !
        end if
        !
      end if
      !
      ge(1) = real(rest,ki)
      ge(2) = aimag(rest)
      !~ write(*,*) 'test ge :',rest,abserr
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
    function gl(n,a,b,c)
      !
      integer, intent(in) :: n
      real(ki), intent(in) :: a,b,c
      real(ki), dimension(2) :: gl
      !
      complex(ki) :: veri,veri_rat,rest
      real(ki) :: x1,x2
      complex(ki) :: cx1,cx2
      real(ki) :: delta
      !
      !~ plus_grand_glob = max(abs(a),abs(b),abs(c))
      delta = b*b-4._ki*a*c
      !
      if (delta >= 0._ki) then
        !
        x1 = (-b + sqrt(delta))/(2._ki*a)
        x2 = (-b - sqrt(delta))/(2._ki*a)
        x1_glob = cmplx(x1,0._ki)
        x2_glob = cmplx(x2,0._ki)
        !
        if (n == 1) then
          !
          veri_rat = -2._ki
          veri =  ( z_log(a,-1._ki) + (1._ki-x1)*z_log(1._ki-x1,-1._ki) + x1*z_log(-x1,-1._ki) &
                  + (1._ki-x2)*z_log(1._ki-x2,1._ki) + x2*z_log(-x2,1._ki) )
          !
        else if (n == 2) then
          !
          veri_rat = -(1._ki + x1 + x2)/2._ki
          veri =  ( z_log(a,-1._ki) + (1._ki-x1**2)*z_log(1._ki-x1,-1._ki) + x1**2*z_log(-x1,-1._ki) &
                  + (1._ki-x2**2)*z_log(1._ki-x2,1._ki) + x2**2*z_log(-x2,1._ki) )/2._ki
          !
        else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function gl (file mod_gn.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'value of n not implemented'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'n= %d0'
          tab_erreur_par(3)%arg_int = n
          call catch_exception(0)
          !
          ! to please the compiler
          stop
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
        x1_glob = cx1
        x2_glob = cx2
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
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function gl (file mod_gn.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'value of n not implemented'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'n= %d0'
          tab_erreur_par(3)%arg_int = n
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
        end if
        !
        veri =  veri + veri_rat
        rest = veri
        !
      end if
      !
      gl(1) = real(rest,ki)
      gl(2) = aimag(rest)
      !
      !
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
    !~ function gf(n,a,b,c)
    function gf(n,a,b,c,dist)
      !
      integer, intent(in) :: n
      real(ki), intent(in) :: a,b,c
      logical, intent(in) :: dist
      real(ki), dimension(2) :: gf
      !
      complex(ki) :: veri,rest,abserr
      complex(ki) :: rest1,abserr1,rest2,abserr2,rest3,abserr3
      complex(ki) :: residue1,residue2
      complex(ki) :: extra_imag1,extra_imag2
      real(ki) :: pole1,pole2
      real(ki) :: x1,x2,deltax
      complex(ki) :: cx1,cx2,cdeltax
      real(ki) :: delta
      logical :: l1,l2
      complex(ki) :: div_part
      !
      plus_grand_glob = max(abs(a),abs(b),abs(c))
      expo_glob = n
      delta = b*b-4._ki*a*c
      a_glob = a
      b_glob = b
      c_glob = c
      lambda_glob = lambda_par
      alpha_glob = alpha_par
      beta_glob = beta_par
      dist_glob = dist
      !write(*,*) 'delta et a :',delta,a,dist
      ! divergent part for Landau singularities
      if (dist) then
        !
        if (delta > 0._ki) then
          !
          if (a > 0._ki) then
            !
            div_part = 2._ki*i_*pi/sqrt(delta)*(log(delta/a) - i_*pi)
            !
          else if (a < 0._ki) then
            !
            div_part = 2._ki*i_*pi/sqrt(delta)*log(-delta/a)
            !
          end if
          !
        else if (delta <= 0._ki) then
          !
          if (a > 0._ki) then
            !
            div_part = 2._ki*pi/sqrt(-delta)*log(-delta/a)
            !
          else if (a < 0._ki) then
            !
            div_part = -2._ki*pi/sqrt(-delta)*(log(delta/a) - i_*pi)
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
        x1 = (-b + sqrt(delta))/(2._ki*a)
        x2 = (-b - sqrt(delta))/(2._ki*a)
        x1_glob = cmplx(x1,0._ki)
        x2_glob = cmplx(x2,0._ki)
        deltax = x1-x2
        !
        if (sqrt(abs(delta)) >  coupure_3p1m_2mi) then
        !~ if (abs(delta) >  0._ki) then
          !
          veri =  ( -zdilog((1._ki-x2)/(1._ki-x1),sign(1._ki,2._ki-x1-x2)) + zdilog(x2/x1,sign(1._ki,-x1-x2)) &
                  + zdilog((1._ki-x1)/(1._ki-x2),sign(1._ki,-2._ki+x1+x2)) - zdilog(x1/x2,sign(1._ki,x1+x2)) &
                  + ( 2._ki*z_log(deltax,1._ki) - i_*pi + z_log(a,-1._ki) ) &
                  *( z_log((x1-1._ki)/x1,1._ki) - z_log((x2-1._ki)/x2,-1._ki) ) &
                  )/deltax/a
          !
          rest = veri
          !
        else
          !
          !~ pole1 = x2
          !~ eps_glob = 1._ki
          !
          !~ call generic_eval_numer(eval_numer_gf1,0._ki,1._ki,1.0E-8_ki,rest1,abserr1)
          !~ call generic_eval_numer(eval_numer_gf1,0._ki,1._ki,tolerance,rest1,abserr1)
          call generic_eval_numer(eval_numer_gf,0._ki,1._ki,tolerance,rest,abserr)
          !
          !~ if ( (pole1 >= 0._ki) .and. (pole1 <= 1._ki) ) then
            !~ !
            !~ residue1 = ( z_log(a_glob,-1._ki) + z_log(x2-x1,-1._ki) )/(x2-x1)/a_glob
            !~ extra_imag1 = -2._ki*i_*pi*residue1
            !~ !
          !~ else
            !~ !
            !~ extra_imag1 = 0._ki
            !~ !
          !~ end if
          !~ !
          !~ rest1 = rest1 + extra_imag1
          !~ !
          !~ pole2 = x1
          !~ eps_glob = -1._ki
          !~ !
          !~ call generic_eval_numer(eval_numer_gf2,0._ki,1._ki,1.0E-8_ki,rest2,abserr2)
          !~ call generic_eval_numer(eval_numer_gf2,0._ki,1._ki,tolerance,rest2,abserr2)
          !~ !
          !~ if ( (pole2 >= 0._ki) .and. (pole2 <= 1._ki) ) then
            !~ !
            !~ residue2 = z_log(x1-x2,1._ki)/(x1-x2)/a_glob
            !~ extra_imag2 = 2._ki*i_*pi*residue2
            !~ !
          !~ else
            !~ !
            !~ extra_imag2 = 0._ki
            !~ !
          !~ end if
          !~ !
          !~ rest2 = rest2 + extra_imag2
          !
          !~ rest = rest1+rest2
          !~ abserr = abserr1 + abserr2
          !~ write(*,*) 'exact result is:',2._ki*i_*pi/sqrt(delta)*(log(delta/a)-i_*pi)
          !~ write(*,*) 'err numer result is:',abserr
          !~ write(*,*) 'analytical result is:', &
          !~ ( -zdilog((1._ki-x2)/(1._ki-x1),sign(1._ki,2._ki-x1-x2)) + zdilog(x2/x1,sign(1._ki,-x1-x2)) &
                  !~ + zdilog((1._ki-x1)/(1._ki-x2),sign(1._ki,-2._ki+x1+x2)) - zdilog(x1/x2,sign(1._ki,x1+x2)) &
                  !~ + ( 2._ki*z_log(deltax,1._ki) - i_*pi + z_log(a,-1._ki) ) &
                  !~ *( z_log((x1-1._ki)/x1,1._ki) - z_log((x2-1._ki)/x2,-1._ki) ) &
                  !~ )/deltax/a
          !~ write(*,*) 'numer result is:',rest
          !~ if (dist) then
          !~ rest = rest + 2._ki*i_*pi/sqrt(delta)*(log(delta/a)-i_*pi)
          rest = rest + div_part
          !~ end if
          !~ write(*,*) 'numer result is:',rest
          !
        end if
        !
      else
        origine_info_par = "gf"
        num_grand_b_info_par = n
        denom_grand_b_info_par = abs(deltax)
        !
        !
        cx1 = -b/(2._ki*a) + i_*sqrt(-delta)/(2._ki*abs(a))
        cx2 = -b/(2._ki*a) - i_*sqrt(-delta)/(2._ki*abs(a))
        x1_glob = cx1
        x2_glob = cx2
        cdeltax = cx1-cx2
        !
        if (sqrt(abs(delta)) >  coupure_3p1m_2mi) then
        !~ if (abs(delta) >  0._ki) then
          !
          veri =  ( -cdilog((1._ki-cx2)/(1._ki-cx1)) + cdilog(cx2/cx1) &
                  + cdilog((1._ki-cx1)/(1._ki-cx2)) - cdilog(cx1/cx2) &
                  + ( 2._ki*log(cdeltax) - i_*pi + z_log(a,-1._ki) ) &
                  *( log((cx1-1._ki)/cx1) - log((cx2-1._ki)/cx2) ) &
                  )/cdeltax/a
          !
          rest = veri
          !
        else
          !
          !~ eps_glob = 1._ki
          !~ call inside_contour(cx2,l1)
          !~ !
          !~ call generic_eval_numer(eval_numer_gf1,0._ki,1._ki,tolerance,rest1,abserr1)
          call generic_eval_numer(eval_numer_gf,0._ki,1._ki,tolerance,rest,abserr)
          !
          !~ if ( l1) then
            !~ !
            !~ residue1 = ( z_log(a_glob,-1._ki) + log(cx2-cx1) )/(cx2-cx1)/a_glob
            !~ extra_imag1 = -2._ki*i_*pi*residue1
            !~ !
          !~ else
            !~ !
            !~ extra_imag1 = 0._ki
            !~ !
          !~ end if
          !~ !
          !~ rest1 = rest1 + extra_imag1
          !~ !
          !~ eps_glob = -1._ki
          !~ call inside_contour(cx1,l2)
          !~ !
          !~ call generic_eval_numer(eval_numer_gf2,0._ki,1._ki,1.0E-8_ki,rest2,abserr2)
          !~ !
          !~ if ( l2 ) then
            !~ !
            !~ residue2 = log(cx1-cx2)/(cx1-cx2)/a_glob
            !~ extra_imag2 = 2._ki*i_*pi*residue2
            !~ !
          !~ else
            !~ !
            !~ extra_imag2 = 0._ki
            !~ !
          !~ end if
          !~ !
          !~ rest2 = rest2 + extra_imag2
          !~ !
          !~ !
          !~ rest = rest2 + rest1
          !~ abserr = abserr1 + abserr2
          !~ if (dist) then
          !~ rest = rest + 2._ki*i_*pi/sqrt(delta)*(log(delta/a)-i_*pi)
          !~ write(*,*) 'analytical result is:', &
          !~ ( -cdilog((1._ki-cx2)/(1._ki-cx1)) + cdilog(cx2/cx1) &
                  !~ + cdilog((1._ki-cx1)/(1._ki-cx2)) - cdilog(cx1/cx2) &
                  !~ + ( 2._ki*log(cdeltax) - i_*pi + z_log(a,-1._ki) ) &
                  !~ *( log((cx1-1._ki)/cx1) - log((cx2-1._ki)/cx2) ) &
                  !~ )/cdeltax/a
          rest = rest + div_part
          !~ write(*,*) 'numer result is:',rest
          !~ end if
          !
        end if
        !
      end if
      !
      gf(1) = real(rest,ki)
      gf(2) = aimag(rest)
      !~ write(*,*) 'test gf :',rest,abserr
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
    !  This is the integrand for the numerical evaluation of ge,
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
    !~ function eval_numer_ge(u)
      !~ !
      !~ real(ki), intent (in) :: u
      !~ complex(ki) :: eval_numer_ge
      !~ !
      !~ real(ki) :: x,y
      !~ real(ki) :: eps
      !~ complex(ki) :: z,jacob
      !~ !
      !~ eps = eps_glob
      !~ x = u
      !~ y = lambda_glob*u**alpha_glob*(1._ki-u)**beta_glob
      !~ jacob = 1._ki - eps*i_*lambda_glob*u**(alpha_glob-1._ki)&
              !~ *(1._ki-u)**(beta_glob-1._ki)*(alpha_glob*(1._ki-u)-beta_glob*u)
      !~ z = x - eps*i_*y
      !~ eval_numer_ge = z**(expo_glob-1)/(z-x1_glob)/(z-x2_glob)/a_glob
      !~ eval_numer_ge = eval_numer_ge*jacob
      !~ !
    !~ end function eval_numer_ge
    function eval_numer_ge(u)
      !
      real(ki), intent (in) :: u
      complex(ki) :: eval_numer_ge
      !
      real(ki) :: x,y
      real(ki) :: eps
      complex(ki) :: z,jacob
      real(ki) :: sigma
      !
      !~ eps = eps_glob
      sigma = -b_glob/a_glob/2._ki
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
      select case(expo_glob)
      !
      case(1)
      !
      eval_numer_ge = -( &
        &1._ki/(a_glob-b_glob*z+c_glob*z*z) + &
        &1._ki/(a_glob+b_glob*z+c_glob*z*z) + &
        &1._ki/(a_glob*z*z-b_glob*z+c_glob) )
      !
      case(2)
      !
      eval_numer_ge = -( &
        &(-2._ki*b_glob) &
        &/( (a_glob-b_glob*z+c_glob*z*z)*(a_glob+b_glob*z+c_glob*z*z) ) + &
        &(-z)/(a_glob*z*z-b_glob*z+c_glob) )
      !
      case(3)
      !
      eval_numer_ge = ( &
        &(-2._ki*b_glob**2 + 2._ki*c_glob*a_glob + 2._ki*c_glob*c_glob*z*z) &
        &/( (a_glob-b_glob*z+c_glob*z*z)*(a_glob+b_glob*z+c_glob*z*z) ) + &
        &(-b_glob*z+c_glob)/(a_glob*z*z-b_glob*z+c_glob) + 1._ki )/a_glob
      !
      case(4)
      !
      eval_numer_ge = -( &
        &2._ki*b_glob*(-b_glob**2 + 2._ki*c_glob*a_glob + 2._ki*c_glob*c_glob*z*z) &
        &/( (a_glob-b_glob*z+c_glob*z*z)*(a_glob+b_glob*z+c_glob*z*z) ) + &
        &((c_glob*a_glob-b_glob**2)*z+b_glob*c_glob)/(a_glob*z*z-b_glob*z+c_glob) &
        & - (a_glob*z-b_glob) )/a_glob**2
      !
      end select
      else
      eval_numer_ge = z**(expo_glob-1)/(a_glob*z*z+b_glob*z+c_glob)
      end if
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
      real(ki) :: eps
      complex(ki) :: z,jacob
      real(ki) :: sigma
      !
      !~ eps = eps_glob
      !
      !~ sigma = -real(b,ki)/a/2._ki
      sigma = -b_glob/a_glob/2._ki
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
      !~ eval_numer_gf = (z_log(a_glob,-1._ki) +log(z-x1_glob) )/(z-x1_glob)/(z-x2_glob)/a_glob
      if (dist_glob) then
      eval_numer_gf = -( &
        &log((a_glob-b_glob*z+c_glob*z*z)/(z*z))/(a_glob-b_glob*z+c_glob*z*z) + &
        &log((a_glob+b_glob*z+c_glob*z*z)/(z*z))/(a_glob+b_glob*z+c_glob*z*z) + &
        &log(a_glob*z*z-b_glob*z+c_glob)/(a_glob*z*z-b_glob*z+c_glob) )
      else
      eval_numer_gf = log(a_glob*z*z+b_glob*z+c_glob)/(a_glob*z*z+b_glob*z+c_glob)
      end if
      eval_numer_gf = eval_numer_gf*jacob
      !
    end function eval_numer_gf
    !
    !****if* src/integrals/three_point/func_gn/eval_numer_gf1
    ! NAME
    !
    !  Function eval_numer_gf1
    !
    ! USAGE
    !
    !  complex = eval_numer_gf1(u)
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
    function eval_numer_gf1(u)
      !
      real(ki), intent (in) :: u
      complex(ki) :: eval_numer_gf1
      !
      real(ki) :: x,y
      real(ki) :: eps
      complex(ki) :: z,jacob
      !
      eps = eps_glob
      x = u
      y = lambda_glob*u**alpha_glob*(1._ki-u)**beta_glob
      jacob = 1._ki - eps*i_*lambda_glob*u**(alpha_glob-1._ki)&
              *(1._ki-u)**(beta_glob-1._ki)*(alpha_glob*(1._ki-u)-beta_glob*u)
      z = x - eps*i_*y
      eval_numer_gf1 = (z_log(a_glob,-1._ki) +log(z-x1_glob) )/(z-x1_glob)/(z-x2_glob)/a_glob
      eval_numer_gf1 = eval_numer_gf1*jacob
      !
    end function eval_numer_gf1
    !
    !****if* src/integrals/three_point/func_gn/eval_numer_gf2
    ! NAME
    !
    !  Function eval_numer_gf2
    !
    ! USAGE
    !
    !  complex = eval_numer_gf2(u)
    !
    ! DESCRIPTION
    !
    !  This is the integrand for the numerical evaluation of gf,
    !  part ln(z-x_2)/( (z-x_1)*(z-x_2) )
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
    function eval_numer_gf2(u)
      !
      real(ki), intent (in) :: u
      complex(ki) :: eval_numer_gf2
      !
      real(ki) :: x,y
      real(ki) :: eps
      complex(ki) :: z,jacob
      !
      eps = eps_glob
      x = u
      y = lambda_glob*u**alpha_glob*(1._ki-u)**beta_glob
      jacob = 1._ki - eps*i_*lambda_glob*u**(alpha_glob-1._ki)&
              *(1._ki-u)**(beta_glob-1._ki)*(alpha_glob*(1._ki-u)-beta_glob*u)
      z = x - eps*i_*y
      eval_numer_gf2 = log(z-x2_glob)/(z-x1_glob)/(z-x2_glob)/a_glob
      eval_numer_gf2 = eval_numer_gf2*jacob
      !
    end function eval_numer_gf2
    !
    !****if* src/integrals/three_point/func_gn/inside_contour
    ! NAME
    !
    !  Subroutine inside_contour
    !
    ! USAGE
    !
    !  call inside_contour(pole,yes_or_no)
    !
    ! DESCRIPTION
    !
    !  This subroutine tests if the pole is inside the contour or not
    !
    ! INPUTS
    !
    !  * pole -- a complex (type ki), the pole
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  It returns a logical, true if the pole is inside, false otherwise
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine inside_contour(pole,yes_or_no)
      !
      complex(ki), INTENT(IN) :: pole
      logical, INTENT(OUT) :: yes_or_no
      !
      real(ki) :: distance,x,y
      !
      yes_or_no = .false.
      x = real(pole,ki)
      y = aimag(pole)
      distance = lambda_glob*x**alpha_glob*(1._ki-x)**beta_glob
      if ( abs(distance-y) <= 0.1_ki ) then
        !
        lambda_glob = 2._ki*lambda_glob
        distance = lambda_glob*x**alpha_glob*(1._ki-x)**beta_glob
        !
      end if
      !
      if ( (x >= 0._ki) .and. (x <= 1._ki) .and. (abs(y) <= distance) .and. (sign(1._ki,y) == sign(1._ki,-eps_glob)) ) then
        !
        yes_or_no = .true.
        !
      end if
      !
    end subroutine inside_contour
    !
end module func_gn

