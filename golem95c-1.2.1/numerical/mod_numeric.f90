!****h* src/numerical/numerical_evaluation
! NAME
!
! Module numerical_evaluation
!
! USAGE
!
!  use numerical_evaluation
!
! DESCRIPTION
!
!  This module contains a generic routine for a one dimensional integration.
!  Up to now, the routine used is adapt_gauss1 (in file mod_adapt_gauss.f90).
!  To add a new integration routine, wrap it in a module, load this module 
!  in this file using use association and add a new if case in the routine generic_eval_numer
!  also modify the value of choix accordingly. Of course, do not forget to modify the
!  Makefile (or better the script configure,pl) in such a way that this new module is compiled
!
! OUTPUT
!
!  With this module, one can access to the routine generic_eval_numer
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * parametre (src/module/parametre.f90)
!  * adapt_gauss (src/numerical/mod_adapt_gauss.f90)
!
!*****
!
module numerical_evaluation
  !
  use precision_golem
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use adapt_gauss
  implicit none
  !
  integer, parameter :: choix = 1
  private :: ki,choix
  public :: generic_eval_numer
  !
  contains
    !
    !****f* src/numerical/numerical_evaluation/generic_eval_numer
    ! NAME
    !
    !  Subroutine generic_eval_numer
    !
    ! USAGE
    !
    !  call generic_eval_numer(func,b_inf,b_sup,tol,rest,abserr)
    !
    ! DESCRIPTION
    !
    !  Generic routine for the one dimensional integration.
    !
    ! INPUTS
    !
    !  * func -- an external function as declared by the interface block
    !  * b_inf -- a real (type ki), the lower bound of the integration range
    !  * b_sup -- a real (type ki), the upper bound of the integration range
    !  * tol -- a real (type ki), the tolerance asked by the user
    !
    ! SIDE EFFECTS
    !
    !  no side effects
    !
    ! RETURN VALUE
    !
    !  * rest -- a complex (type ki), the result of the integration
    !  * abserr -- a complex (type ki), the absolute value of the estimated error
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine generic_eval_numer(func,b_inf,b_sup,tol,rest,abserr)
      !
      real(ki), intent (in) :: b_inf,b_sup,tol
      complex(ki), intent (out) :: rest
      complex(ki), intent (out) :: abserr
      !
      interface
        !
        function func(x)
          use precision_golem
          real(ki), intent (in) :: x
          complex(ki) :: func
        end function func
        !
      end interface
      !
      if (choix == 1) then
        !
        call adapt_gauss1(func,b_inf,b_sup,tol,rest,abserr)
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine generic_eval_numer'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the value of the variable choix is not correct : choix=%d0'
        tab_erreur_par(2)%arg_int = choix
        call catch_exception(0)
        !
      end if
      !
    end subroutine generic_eval_numer
    !
end module numerical_evaluation
