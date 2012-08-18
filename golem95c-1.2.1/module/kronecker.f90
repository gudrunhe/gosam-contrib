! 
!****h* src/module/kronecker
! NAME
!
!  Module kronecker
!
! USAGE
!
!  use kronecker
!
! DESCRIPTION
!
!  This module contains two functions delta and deltab which correspond 
!  respectively to the Kronecker symbol \delta_{ij}
!  and 1-\delta_{ij}. These two functions have two integer arguments and 
!  return an integer
!
! OUTPUT
!
!  This module exports two functions:
!  * delta -- the Kronecker symbol
!  * deltab -- 1-delta
!
! USES
!
!  No uses
!
!*****
module kronecker
  !
  implicit none
  !
  contains
    !
    !****f* src/module/kronecker/delta
    ! NAME
    !
    !  Function delta
    !
    ! USAGE
    !
    !  integer = delta(i,j)
    !
    ! DESCRIPTION
    !
    !  This is the Kronecker symbol \delta_{ij}
    !
    ! INPUTS
    !
    !  * i -- an integer
    !  * j -- an integer
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  It returns an integer 0 or 1
    !
    ! EXAMPLE
    !
    !  k = delta(2,3) --> k = 0
    !  k = delta(3,3) --> k = 1
    !
    !*****
    function delta(i,j)
      !
      integer, intent (in) :: i
      integer, intent (in) :: j
      integer :: delta
      !
      if (i == j) then
        !
        delta = 1
        !
      else
        !
        delta = 0
        !
      end if
      !
    end function delta
    !
    !****f* src/module/kronecker/deltab
    ! NAME
    !
    !  Function deltab
    !
    ! USAGE
    !
    !  integer = deltab(i,j)
    !
    ! DESCRIPTION
    !
    !  This is one minus the Kronecker symbol, 1-\delta_{ij}
    !
    ! INPUTS
    !
    !  * i -- an integer
    !  * j -- an integer
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  It returns an integer 0 or 1
    !
    ! EXAMPLE
    !
    !  k = deltab(2,3) --> k = 1
    !  k = deltab(3,3) --> k = 0
    !
    !*****
    function deltab(i,j)
      !
      integer, intent (in) :: i
      integer, intent (in) :: j
      integer :: deltab
      !
      if (i == j) then
        !
        deltab = 0
        !
      else
        !
        deltab = 1
        !
      end if
      !
    end function deltab
    !
end module kronecker
