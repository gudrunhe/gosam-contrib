! 
!****h* src/module/translate
! NAME
!
!  Module translate
!
! USAGE
!
!  use translate
!
! DESCRIPTION
!
!  This module is used to translate an array of n (=2m) reals into an array
!  of m complexs
!
! OUTPUT
!
!  It exports:
!  * to_complex -- a subroutine to translate an array of n (=2m) reals into an array
!                  of m complexs
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!
!*****
module translate
  !
  use precision_golem
  use sortie_erreur
  implicit none
  !
  private 
  public :: to_complex
  contains
    !
    !****f* src/module/translate/to_complex
    ! NAME
    !
    !  Subroutine to_complex
    !
    ! USAGE
    !
    !  call to_complex(t,z)
    !
    ! DESCRIPTION
    !
    !  This subroutine transforms an array of reals of rank 1 and shape 2*m 
    !  t in an array of complexs of size m z, it returns z(i) = t(i) + i_*t(i+1). 
    !  If size of t is odd, the subroutine to_complex returns an error
    !
    ! INPUTS
    !
    !  * t -- a real array (type ki) of rank 1
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * z -- a complex array (type ki) of rank 1 and shape size(t)/2
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine to_complex(t,z)
      !
      real(ki), intent(in), dimension(:) :: t
      complex(ki), intent(out) ,dimension(:) :: z
      !
      integer :: dim_t,i,j
      !
      dim_t = size(t)
      !
      if (mod(dim_t,2) == 0) then
        !
        do i = 1,dim_t,2
          !
          j = (i+1)/2
          z(j) = cmplx(t(i),t(i+1),ki)
          !
        end do
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'error in subroutine to_complex'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'The size of the first argument array is odd %d0'
        tab_erreur_par(2)%arg_int = dim_t
        call catch_exception(0)
        !
      end if
      !
    end subroutine to_complex
  !
end module translate
