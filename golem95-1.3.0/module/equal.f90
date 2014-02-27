! 
!****h* src/module/equal
! NAME
!
!  Module equal
!
! USAGE
!
!  use equal
!
! DESCRIPTION
!
!  This module is used to compare two objects
!
! OUTPUT
!
!  This module exports two functions:
!  * equal_real -- to compare two reals
!  * cut_s -- sets the momentum s equal to zero, if the ratio s/M or the absoulute value
!             is smaller than some global parameters.
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!  * parametre (src/module/parametre.f90)
!  * constante (src/module/constante.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!
!*****
module equal
  !
  use precision_golem
  use parametre, only: cut_s_over_m, cut_s_abs
  use constante, only: zero, czero
  use sortie_erreur, only : tab_erreur_par,catch_exception
  implicit none
  !
  private
  !
  interface cut_s
     !
     module procedure cut_s_s
     module procedure cut_s_sm_r, cut_s_sm_c
     module procedure cut_s_smm_r, cut_s_smm_c
     !
  end interface
  !
  public :: equal_real, cut_s
  !
  contains
    !
    !****f* src/module/equal/equal_real
    ! NAME
    !
    !  Function equal_real
    !
    ! USAGE
    !
    !  logical = equal_real(xa,xb,echelle)
    !
    ! DESCRIPTION
    !
    !  This function compares two real (of same kind) by computing their 
    !  difference and compares it to epsilon (of the same kind)
    !  true if | xa - xb | <= echelle*epsilon, false otherwise
    !
    ! INPUTS
    !
    !  * xa -- a real of type ki
    !  * xb -- a real of type ki
    !  * echelle -- a real of type ki (optional)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  It returns a logical
    !
    ! EXAMPLE
    !
    !
    !*****
    pure elemental function equal_real(xa,xb,echelle)
      implicit none
      !
      real(ki), intent(in) :: xa,xb
      real(ki), intent(in), optional :: echelle
      logical equal_real
      !
      real(ki) :: my_epsilon
      !
      my_epsilon = epsilon(1._ki)
      if (present(echelle)) my_epsilon=echelle*my_epsilon
      equal_real = abs(xa-xb) <= my_epsilon
      !
    end function equal_real
    !
    !
    !****f* src/module/equal/cut_s
    ! NAME
    !
    !  subroutine cut_s
    !
    ! USAGE
    !
    !  call cut_s(s,m1,m2)
    !  call cut_s(s,m)
    !  call cut_s(s)
    !
    ! DESCRIPTION
    !
    !  This function sets s to zero, if either of the following conditions is fulfilled.
    !  abs(s) is smaller than cut_s_abs.
    !  abs(s/Sum(m)) is smaller than cut_s_over_m.
    !  Calling this routine improves stability in the form factor calculations.
    !
    ! INPUTS
    !
    !  * s -- a real of type ki
    !  * m, m1, m2 -- real/complex type
    !
    ! SIDE EFFECTS
    !
    !  Possible change of s to zero.
    !
    ! RETURN VALUE
    !
    ! EXAMPLE
    !
    !
    !*****
    subroutine cut_s_smm_r(s,mass1,mass2)
      implicit none
      real(ki), intent(inout) :: s
      real(ki), intent(in) :: mass1, mass2
      !
      real(ki) :: sum_mass
      !
      sum_mass = mass1 + mass2
      call cut_s(s, sum_mass)
      !
    end subroutine cut_s_smm_r
    !
    subroutine cut_s_sm_r(s,mass)
      implicit none
      real(ki), intent(inout) :: s
      real(ki), intent(in) :: mass
      !
      if (equal_real(mass,zero)) then
         !
         call cut_s(s)
         !
      else if ( abs(s/mass) .le. cut_s_over_m) then
         !
         s = zero
         !
      else if ( abs(s) .le. cut_s_abs ) then
         !
         tab_erreur_par(1)%a_imprimer = .true.
         tab_erreur_par(1)%chaine = 'in cut_s: s is set to zero because its absolute value is lower than cut_s_abs!'
         tab_erreur_par(2)%a_imprimer = .true.
         tab_erreur_par(2)%chaine = 'The ratio s/Mass is not lower than the parameter cut_s_over_m!'
         tab_erreur_par(3)%a_imprimer = .true.
         tab_erreur_par(3)%chaine = 's= %f0'
         tab_erreur_par(3)%arg_real = s
         tab_erreur_par(4)%a_imprimer = .true.
         tab_erreur_par(4)%chaine = 'mass= %f0'
         tab_erreur_par(4)%arg_real = mass
         !
         call catch_exception(1)
         !
         s = zero
         !
      end if
      !
    end subroutine cut_s_sm_r
    !
    subroutine cut_s_s(s)
      implicit none
      real(ki), intent(inout) :: s
      !
      if ( s /= zero) then
         !
         if ( abs(s) .le. cut_s_abs ) then
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'in cut_s: s is set to zero because its absolute value is lower than cut_s_abs!'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = 's= %f0'
            tab_erreur_par(2)%arg_real = s
            !
            call catch_exception(1)
            !
            s = zero
            !
         end if
         !
      end if
      !
    end subroutine cut_s_s
    !
    subroutine cut_s_smm_c(s,mass1,mass2)
      implicit none
      real(ki), intent(inout) :: s
      complex(ki), intent(in) :: mass1, mass2
      !
      real(ki) :: sum_mass
      !
      sum_mass = real(mass1 + mass2,ki)
      call cut_s(s, sum_mass)
      !
    end subroutine cut_s_smm_c
    !
    subroutine cut_s_sm_c(s,mass_c)
      implicit none
      real(ki), intent(inout) :: s
      complex(ki), intent(in) :: mass_c
      !
      real(ki) :: mass
      !
      mass = real(mass_c,ki)
      !
      if (equal_real(mass,zero)) then
         !
         call cut_s(s)
         !
      else if ( abs(s/mass) .le. cut_s_over_m) then
         !
         s = zero
         !
      else if ( abs(s) .le. cut_s_abs ) then
         !
         tab_erreur_par(1)%a_imprimer = .true.
         tab_erreur_par(1)%chaine = 'in cut_s: s is set to zero because its absolute value is lower than cut_s_abs!'
         tab_erreur_par(2)%a_imprimer = .true.
         tab_erreur_par(2)%chaine = 'The ratio s/Mass is not lower than the parameter cut_s_over_m!'
         tab_erreur_par(3)%a_imprimer = .true.
         tab_erreur_par(3)%chaine = 's= %f0'
         tab_erreur_par(3)%arg_real = s
         tab_erreur_par(4)%a_imprimer = .true.
         tab_erreur_par(4)%chaine = 'mass= %f0'
         tab_erreur_par(4)%arg_real = mass
         !
         call catch_exception(1)
         !
         s = zero
         !
      end if
      !
    end subroutine cut_s_sm_c
    !
end module equal
