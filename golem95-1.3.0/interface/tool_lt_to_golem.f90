! 
!****h* src/interface/tool_lt_to_golem
! NAME
!
!  Module tool_lt_to_golem
!
! USAGE
!
!  use tool_lt_to_golem
!
! DESCRIPTION
!
!  This module contains one function to build the interface between LoopTools
!  and Golem
!
! OUTPUT
!
!  This module exports one function:
!  * extract -- extract the numbers contained in a string
!
! USES
!
!
!
!*****
module tool_lt_to_golem
  !
  implicit none
  !
  private 
  !
  public :: extract
  !
  contains
    !
    !****f* src/interface/tool_lt_to_golem/extract
    ! NAME
    !
    !  Subroutine extract
    !
    ! USAGE
    !
    !  call extract(chaine,tab_int)
    !
    ! DESCRIPTION
    !
    !  This routine takes a string of characters, extracts the numbers and puts them into an array
    !
    ! INPUTS
    !
    !  * chaine -- a character of unknown length
    !  * tab_int -- an integer array of rank 1 whose extend is length(chaine), 
    !               it is filled with -1
    !
    ! SIDE EFFECTS
    !
    !  no side effect
    !
    ! RETURN VALUE
    !
    !  it returns an integer array of rank 1
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine extract(chaine,tab_int)
      !
      character (len=*), intent (in) :: chaine
      integer, dimension(len(chaine)) :: tab_int
      !
      integer :: i
      character (len=10) :: chiffre = '0123456789'
      character (len=1) :: c
      !
      tab_int = -1
      !
      do i=1, len(chaine)
        !
        c = chaine(i:i+1)
        !
        if (verify(c,chiffre) == 0) then
          !
          tab_int(i) = iachar(c) - 48
          !
        end if
        !
      end do
      !
    end subroutine extract
    !
end module tool_lt_to_golem
