!
!****h* src/module/parametre
! NAME
!
!  Module parametre
!
! USAGE
!
!  use parametre
!
! DESCRIPTION
!
!  This module is used to pass some variables used by many functions of 
!  the GOLEM program. Note that these variables can be rewritten.
!  It contains also a routine to print the parameters (it prints on unit 6)
!
! OUTPUT
!
!  It exports the variables:
!  * tolerance -- a real (type ki), the tolerance for the numerical integration
!  * lambda_par -- a real (type ki), a parameter of the contour deformation
!  * alpha_par -- a real (type ki), a parameter of the contour deformation
!  * beta_par -- a real (type ki), a parameter of the contour deformation
!  * coupure_3p2m -- a real (type ki), a cut between numerical and analytical 
!                    computation for two mass three point functions
!  * coupure_3p3m -- a real (type ki), a cut between numerical and analytical 
!                    computation for three mass three point functions
!  * coupure_4p1m -- a real (type ki), a cut between numerical and analytical 
!                    computation for one mass four point functions
!  * coupure_4p2m_opp -- a real (type ki), a cut between numerical and analytical 
!                    computation for two mass opposite four point functions
!  * coupure_4p2m_adj -- a real (type ki), a cut between numerical and analytical 
!                    computation for two mass adjacent four point functions
!  * coupure_4p3m -- a real (type ki), a cut between numerical and analytical 
!                    computation for three mass four point functions
!  * coupure_4p4m -- a real (type ki), a cut between numerical and analytical 
!                    computation for four mass four point functions (not active)
!  * coupure_3p2m_1mi --  a real (type ki), a cut between numerical and analytical 
!                    computation for one internal mass two external mass three point functions
!  * rat_or_tot_par -- a character (len=3)
!  * rmass_or_cmass_par -- a character (len=5)
!  * if_print_info_par -- a logical, if true it prints some informations concerning the numerical
!                    integration
!  * if_print_warn_par -- a logical, if true it prints some informations concerning the warning 
!                    about numerical precision
!  * accuracy_par -- the accuracy for the matrix inversion and the numerical integration
!  * not_enough_accuracy_par -- a flag to ring the bell if the accuracy is not reached in
!                               in the matrix inversion and the numerical integration
!  * mu2_scale_par -- the square of the renormalisation scale
!  * subroutine print_parameter -- to print these variables
!  * withlt  -- flag to use LoopTools instead of avh_olo for finite D0,C0 
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!
!*****
module parametre
  !
  use precision_golem
  implicit none
  !
  ! here everything is public except what has been defined explicitly as private
  private :: ki
  !
  type rat_or_tot_string
     character(len=3) :: image
     logical :: rat_selected
     logical :: tot_selected
  end type rat_or_tot_string
  !
  type rmass_or_cmass_string
     character(len=5) :: image
     logical :: rmass_selected
     logical :: cmass_selected
  end type rmass_or_cmass_string
  !
  !
  real(ki),save :: tolerance = 1.e-8_ki  ! precision for Gaussian integration
  !
  real(ki),save :: lambda_par = 1._ki  ! parameters of the contour
  real(ki),save :: alpha_par = 1._ki   ! deformation
  real(ki),save :: beta_par = 1._ki    !
  !
  real(ki),save :: coupure_3p1m_2mi = 5.e-3_ki ! value to switch between analytical and numerical
                                           ! evaluation for one mass three point functions (mod_gn)
  real(ki),save :: coupure_3p2m = 5.e-3_ki ! value to switch between analytical and numerical
                                           ! evaluation for two mass three point functions
  real(ki),save :: coupure_3p3m = 5.e-3_ki ! value to switch between analytical and numerical
                                           ! evaluation for three mass three point functions
  real(ki),save :: coupure_4p1m = 5.e-3_ki ! value to switch between analytical and numerical
                                           ! evaluation for one/zero mass four point functions
  real(ki),save :: coupure_4p2m_opp = 5.e-3_ki ! value to switch between analytical and numerical
                                               ! evaluation for two opposite mass four point functions
  real(ki),save :: coupure_4p2m_adj = 5.e-3_ki! value to switch between analytical and numerical
                                              ! evaluation for two adjacent mass four point functions
  real(ki),save :: coupure_4p3m = 5.e-3_ki ! value to switch between analytical and numerical
                                           ! evaluation for three mass four point functions
  real(ki),save :: coupure_4p4m = 0._ki ! value to switch between analytical and numerical
                                        ! evaluation for four mass four point functions
  real(ki),save :: coupure_3p2m_1mi = 5.e-3_ki ! value to switch between analytical and numerical
                                        ! evaluation for one internal mass two external mass three point functions
  real(ki),save :: cut_s_abs = 10._ki*epsilon(1._ki)
  real(ki),save :: cut_s_over_m = 1000000._ki*epsilon(1._ki)

  type(rat_or_tot_string),parameter :: tot = rat_or_tot_string('tot', .false., .true.)
  type(rat_or_tot_string),parameter :: rat = rat_or_tot_string('rat', .true., .false.)

  type(rat_or_tot_string),save :: rat_or_tot_par = tot
  !type(rat_or_tot_string),save :: rat_or_tot_par = rat


  type(rmass_or_cmass_string),parameter :: rmass = rmass_or_cmass_string('rmass', .true., .false.)
  type(rmass_or_cmass_string),parameter :: cmass = rmass_or_cmass_string('cmass', .false., .true.)

  type(rmass_or_cmass_string),save :: rmass_or_cmass_par = cmass

  logical, save :: if_print_info_par = .false. ! if true print information
  logical, save :: if_print_warn_par = .false. ! if true print information for warning
  real(ki),save :: accuracy_par = 1.e-10_ki ! the accuracy for the matrix inversion and the numerical integration
  logical, save :: not_enough_accuracy_par = .false. !
  real(ki), save :: mu2_scale_par = 1._ki ! the square of the renormalisation scale
  logical, save :: olo = .false.  ! flag set to true as soon as avh_olo has been called once
  ! added to include LT option Jan2011
  logical, save :: withlt = .false.  

 interface assignment(=)
    module procedure assign_rat_or_tot_string
 end interface

 interface operator(==)
    module procedure equals_rat_or_tot_string
    module procedure equals_rat_or_tot_string_revd
 end interface

 private :: assign_rat_or_tot_string
 private :: equals_rat_or_tot_string
 private :: equals_rat_or_tot_string_revd

 interface assignment(=)
    module procedure assign_rmass_or_cmass_string
 end interface

 interface operator(==)
    module procedure equals_rmass_or_cmass_string
    module procedure equals_rmass_or_cmass_string_r
 end interface

 private :: assign_rmass_or_cmass_string
 private :: equals_rmass_or_cmass_string
 private :: equals_rmass_or_cmass_string_r

 contains
   !
    pure subroutine assign_rat_or_tot_string(rot, ch)
       implicit none
       type(rat_or_tot_string), intent(out) :: rot
       character(len=3), intent(in) :: ch

       rot%image = ch
       rot%rat_selected = ch .eq. 'rat'
       rot%tot_selected = ch .eq. 'tot'
    end  subroutine assign_rat_or_tot_string
    !
    pure function equals_rat_or_tot_string(rot, ch) result(test)
       implicit none
       type(rat_or_tot_string), intent(in) :: rot
       character(len=3), intent(in) :: ch
       logical :: test

       test = rot%image .eq. ch
    end  function equals_rat_or_tot_string
    !
    pure function equals_rat_or_tot_string_revd(ch, rot) result(test)
       implicit none
       character(len=3), intent(in) :: ch
       type(rat_or_tot_string), intent(in) :: rot
       logical :: test

       test = rot%image .eq. ch
    end  function equals_rat_or_tot_string_revd
    !
    !
    pure subroutine assign_rmass_or_cmass_string(roc, ch)
      implicit none
      type(rmass_or_cmass_string), intent(out) :: roc
      character(len=5), intent(in) :: ch
      
      roc%image = ch
      roc%rmass_selected = ch .eq. 'rmass'
      roc%cmass_selected = ch .eq. 'cmass'
    end  subroutine assign_rmass_or_cmass_string
    !
    pure function equals_rmass_or_cmass_string(roc, ch) result(test)
      implicit none
      type(rmass_or_cmass_string), intent(in) :: roc
      character(len=5), intent(in) :: ch
      logical :: test
      
      test = roc%image .eq. ch
    end  function equals_rmass_or_cmass_string
    !
    pure function equals_rmass_or_cmass_string_r(ch, roc) result(test)
      implicit none
      character(len=5), intent(in) :: ch
      type(rmass_or_cmass_string), intent(in) :: roc
      logical :: test
      
      test = roc%image .eq. ch
    end  function equals_rmass_or_cmass_string_r
    !
    !
    !****f* src/module/print_parameter
    ! NAME
    !
    !  Subroutine print_parameter
    !
    ! USAGE
    !
    !  call print_parameter()
    !
    ! DESCRIPTION
    !
    !  This routine print the variables defined in the module parametre
    !
    ! INPUTS
    !
    !  No inputs
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  It prints on the unit 6
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine print_parameter()
      !
      integer :: unit
      !
      unit = 6
      !
      write(unit,*) 'tolerance :',tolerance
      write(unit,*) 'lambda_par :',lambda_par
      write(unit,*) 'alpha_par :',alpha_par
      write(unit,*) 'beta_par :',beta_par
      write(unit,*) 'coupure_3p2m :',coupure_3p2m
      write(unit,*) 'coupure_3p3m :',coupure_3p3m
      write(unit,*) 'coupure_4p1m :',coupure_4p1m
      write(unit,*) 'coupure_4p2m_opp :',coupure_4p2m_opp
      write(unit,*) 'coupure_4p2m_adj :',coupure_4p2m_adj
      write(unit,*) 'coupure_4p3m :',coupure_4p3m
      write(unit,*) 'coupure_4p4m :',coupure_4p4m
      write(unit,*) 'rat_or_tot_par : ',rat_or_tot_par
      write(unit,*) 'rmass_or_cmass : ',rmass_or_cmass_par
      write(unit,*) 'if_print_info_par : ',if_print_info_par
      write(unit,*) 'if_print_warn_par : ',if_print_warn_par
      write(unit,*) 'accuracy_par : ',accuracy_par
      write(unit,*) 'not_enough_accuracy_par : ',not_enough_accuracy_par
      write(unit,*) 'mu2_scale_par : ',mu2_scale_par
      write(unit,*) 'accuracy_par : ',accuracy_par
      write(unit,*) 'not_enough_accuracy_par : ',not_enough_accuracy_par
      write(unit,*) 'mu2_scale_par : ',mu2_scale_par
      write(unit,*) 'olo : ',olo
      !
    end subroutine print_parameter
    !
end module parametre
