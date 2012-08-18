!****h* src/module/sortie_erreur
! NAME
!
!  Module sortie_erreur
!
! USAGE
!
!  use sortie_erreur
!
! DESCRIPTION
!
!  This module is used to generate error exception or to print some information from
!  a function/subroutine
!
! OUTPUT
!
!  This module exports:
!  * erreur -- derived type
!  * tab_erreur_par -- an array of 7 derived type erreur
!  * catch_exception -- a subroutine to perform an action depending on the level
!  * print_type -- a subroutine to print the type erreur
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!  * parametre (src/module/parametre.f90)
!  * array (src/module/array.f90)
!
!*****
module sortie_erreur
  !
  use precision_golem
  use parametre, only : if_print_info_par,if_print_warn_par,not_enough_accuracy_par
  use array, only : unpackb
  implicit none
  !
  private
  !
  !****t* src/module/sortie_erreur/erreur
  ! NAME
  !
  !   erreur -- derived type, to print error/info 
  !
  ! SYNOPSIS
  !
  !   type erreur
  !
  ! SOURCE
  !
  type erreur
    !
    character(len=256) :: chaine
    logical :: a_imprimer = .false.
    integer :: arg_int
    real(ki) :: arg_real
    complex(ki) :: arg_comp
    character(len=32) :: arg_char
    integer, dimension(2) :: arg_int_tab
    !
  end type erreur
  !
  ! NOTES
  !
  !   * set erreur%a_imprimer = .true. to print it
  !   * arg_in_tab(1) : packb(tab), arg_int_tab(2) : size(tab)
  !
  !****
  !
  integer :: max_err = 7
  !
  ! an array of 7 derived type is reserved
  !
  type(erreur), dimension(7), save :: tab_erreur_par
  character (len=132),save :: origine_info_par = "" ! the type of function which are integrated
  real(ki),save :: num_grand_b_info_par = 0._ki  ! the numerator of B
  real(ki),save :: denom_grand_b_info_par = 0._ki  ! the denominator of B
  character (len=22),save :: origine_inv_info_par = "" ! the size of the matrix to inverse
  !
  public :: erreur,tab_erreur_par,catch_exception
  public :: origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
  public :: origine_inv_info_par
  !
  contains
    !
    !****f* src/module/sortie_erreur/catch_exception
    ! NAME
    !
    !  Subroutine catch_exception
    !
    ! USAGE
    !
    !  call catch_exception(level)
    !
    ! DESCRIPTION
    !
    !  For the error exception
    !  This routine prints on the unit 0 (stderr for fortran)  
    !  the array tab_erreur_par 
    !
    ! INPUTS
    !
    !  * level -- a integer : 0 the program stops, 1 warning, 2 info
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  write on the unit 0 (stderr for fortran) : level=0,1
    !  or write on the unit 12 : level=2 
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine catch_exception(level)
      !
      integer, intent(in) :: level
      !
      integer :: unit
      integer :: i
      !
      select case(level)
        !
        case(0)
          !
          unit = 0
          !
          write(unit,*) '+++++++++++++++ERROR+++++++++++++++++++++++'
          write(unit,*) 'The program stops because'
          !
          do i=1,max_err
            !
            if (tab_erreur_par(i)%a_imprimer) call print_type(unit,tab_erreur_par(i))
            !
          end do
          !
          stop
          !
        case(1)
          !
          unit = 0
          !
          if (if_print_warn_par) then
            !
            write(unit,*) '+++++++++++++++WARNING+++++++++++++++++++++++'
            !
            do i=1,max_err
              !
              if (tab_erreur_par(i)%a_imprimer) then
                 call print_type(unit,tab_erreur_par(i))
                 tab_erreur_par(i)%a_imprimer = .false.
              end if
              !
            end do
            write(unit,*) 'Type of Feynman integrals :',trim(origine_info_par)
            write(unit,*) 'Numerator of B :',num_grand_b_info_par
            write(unit,*) 'Denominator of B :',denom_grand_b_info_par
            write(unit,*) 'Type of matrix :',trim(origine_inv_info_par)
            !
          end if
          !
          not_enough_accuracy_par = .true.
          !
        case(2)
          !
          if (if_print_info_par) then
            !
            unit = 12
            !
            write(unit,*) '+++++++++++++++++INFO++++++++++++++++++++++'
            !
            do i=1,max_err
              !
              if (tab_erreur_par(i)%a_imprimer) then 
                 call print_type(unit,tab_erreur_par(i))
                 tab_erreur_par(i)%a_imprimer = .false.
              end if
              !
            end do
            !
         else
            !
            do i=1,max_err
               !
               tab_erreur_par(i)%a_imprimer = .false.
               !
            end do
            !
         end if
          !
        case default
          !
          unit = 0
          !
          write(unit,*) 'The level argument of the routine catch_exception must be less or equal than 2'
          write(unit,*) 'this argument is :',level
          !
          stop
          !
      end select
      !
    end subroutine catch_exception
    !
    !****if* src/module/sortie_erreur/print_type
    ! NAME
    !
    !  Subroutine print_type
    !
    ! USAGE
    !
    !  call print_type(unit,type_err)
    !
    ! DESCRIPTION
    !
    !  For the error exception
    !  This routine prints on the unit 0 (stderr for fortran)  
    !  a string of characters and an integer/real/complex/string of 
    !  characters or an array of integers/reals
    !
    ! INPUTS
    !
    !  * unit -- an integer, the unit where to print
    !  * type_err -- type(erreur), the derived type which will printed
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  write on the unit unit
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine print_type(unit,type_err)
      !
      integer, intent(in) :: unit
      type(erreur), intent(in) :: type_err
      !
      integer :: i
      character(len=3), dimension(5) :: car =(/'%d0','%f0','%z0','%c0','%d1'/)
      integer, dimension(5) :: l
      !
      do i=1,size(l)
        !
        l(i) = index(trim(type_err%chaine),car(i))
        !
      end do
      !
      if (maxval(l) == 0) then
        !
        write(unit,*) trim(type_err%chaine)
        !
      else if (l(1) /= 0) then
        !
        write(unit,*) type_err%chaine(1:l(1)-1),type_err%arg_int
        !
      else if (l(2) /= 0) then
        !
        write(unit,*) type_err%chaine(1:l(2)-1),type_err%arg_real
        !
      else if (l(3) /= 0) then
        !
        write(unit,*) type_err%chaine(1:l(3)-1),type_err%arg_comp
        !
      else if (l(4) /= 0) then
        !
        write(unit,*) type_err%chaine(1:l(4)-1),type_err%arg_char
        !
      else if (l(5) /= 0) then
        !
        write(unit,*) type_err%chaine(1:l(5)-1),unpackb(type_err%arg_int_tab(1),type_err%arg_int_tab(2))
        !
      end if
      !
    end subroutine print_type
    !
end module  sortie_erreur
