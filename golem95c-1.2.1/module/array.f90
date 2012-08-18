! 
!****h* src/module/array
! NAME
!
!  Module array
!
! USAGE
!
!  use array
!
! DESCRIPTION
!
!  This module contains six functions which enable set manipulations knowing that
!  a set of integers is represented with the digits of an integer. The six functions
!  are : packb, unpackb, pminus, punion, countb and locateb
!
! OUTPUT
!
!  This module exports six functions:
!
!  * packb -- to transform a set of integers into an integer (unique transformation)
!  * unpackb -- to perform the inverse operation as packb do
!  * pminus -- to subtract two sets
!  * punion -- to add two sets
!  * countb -- to count the number of element of the set
!  * locateb -- to give the location of an element in a set
!
! USES
!
!  none
!
!*****
!
module array
  !
  implicit none
  !
  integer, dimension(0:255), parameter :: bit_count = (/0,1,1,2,1,2,2,3,1,2,2,3&
       &,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2&
       &,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4&
       &,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5&
       &,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6&
       &,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3&
       &,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5&
       &,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8/)
  !
  integer, dimension(0:2047), parameter :: bit_sets = (/&
       &-1,-1,-1,-1,-1,-1,-1,-1,&
       &0,-1,-1,-1,-1,-1,-1,-1,&
       &1,-1,-1,-1,-1,-1,-1,-1,&
       &0,1,-1,-1,-1,-1,-1,-1,&
       &2,-1,-1,-1,-1,-1,-1,-1,&
       &0,2,-1,-1,-1,-1,-1,-1,&
       &1,2,-1,-1,-1,-1,-1,-1,&
       &0,1,2,-1,-1,-1,-1,-1,&
       &3,-1,-1,-1,-1,-1,-1,-1,&
       &0,3,-1,-1,-1,-1,-1,-1,&
       &1,3,-1,-1,-1,-1,-1,-1,&
       &0,1,3,-1,-1,-1,-1,-1,&
       &2,3,-1,-1,-1,-1,-1,-1,&
       &0,2,3,-1,-1,-1,-1,-1,&
       &1,2,3,-1,-1,-1,-1,-1,&
       &0,1,2,3,-1,-1,-1,-1,&
       &4,-1,-1,-1,-1,-1,-1,-1,&
       &0,4,-1,-1,-1,-1,-1,-1,&
       &1,4,-1,-1,-1,-1,-1,-1,&
       &0,1,4,-1,-1,-1,-1,-1,&
       &2,4,-1,-1,-1,-1,-1,-1,&
       &0,2,4,-1,-1,-1,-1,-1,&
       &1,2,4,-1,-1,-1,-1,-1,&
       &0,1,2,4,-1,-1,-1,-1,&
       &3,4,-1,-1,-1,-1,-1,-1,&
       &0,3,4,-1,-1,-1,-1,-1,&
       &1,3,4,-1,-1,-1,-1,-1,&
       &0,1,3,4,-1,-1,-1,-1,&
       &2,3,4,-1,-1,-1,-1,-1,&
       &0,2,3,4,-1,-1,-1,-1,&
       &1,2,3,4,-1,-1,-1,-1,&
       &0,1,2,3,4,-1,-1,-1,&
       &5,-1,-1,-1,-1,-1,-1,-1,&
       &0,5,-1,-1,-1,-1,-1,-1,&
       &1,5,-1,-1,-1,-1,-1,-1,&
       &0,1,5,-1,-1,-1,-1,-1,&
       &2,5,-1,-1,-1,-1,-1,-1,&
       &0,2,5,-1,-1,-1,-1,-1,&
       &1,2,5,-1,-1,-1,-1,-1,&
       &0,1,2,5,-1,-1,-1,-1,&
       &3,5,-1,-1,-1,-1,-1,-1,&
       &0,3,5,-1,-1,-1,-1,-1,&
       &1,3,5,-1,-1,-1,-1,-1,&
       &0,1,3,5,-1,-1,-1,-1,&
       &2,3,5,-1,-1,-1,-1,-1,&
       &0,2,3,5,-1,-1,-1,-1,&
       &1,2,3,5,-1,-1,-1,-1,&
       &0,1,2,3,5,-1,-1,-1,&
       &4,5,-1,-1,-1,-1,-1,-1,&
       &0,4,5,-1,-1,-1,-1,-1,&
       &1,4,5,-1,-1,-1,-1,-1,&
       &0,1,4,5,-1,-1,-1,-1,&
       &2,4,5,-1,-1,-1,-1,-1,&
       &0,2,4,5,-1,-1,-1,-1,&
       &1,2,4,5,-1,-1,-1,-1,&
       &0,1,2,4,5,-1,-1,-1,&
       &3,4,5,-1,-1,-1,-1,-1,&
       &0,3,4,5,-1,-1,-1,-1,&
       &1,3,4,5,-1,-1,-1,-1,&
       &0,1,3,4,5,-1,-1,-1,&
       &2,3,4,5,-1,-1,-1,-1,&
       &0,2,3,4,5,-1,-1,-1,&
       &1,2,3,4,5,-1,-1,-1,&
       &0,1,2,3,4,5,-1,-1,&
       &6,-1,-1,-1,-1,-1,-1,-1,&
       &0,6,-1,-1,-1,-1,-1,-1,&
       &1,6,-1,-1,-1,-1,-1,-1,&
       &0,1,6,-1,-1,-1,-1,-1,&
       &2,6,-1,-1,-1,-1,-1,-1,&
       &0,2,6,-1,-1,-1,-1,-1,&
       &1,2,6,-1,-1,-1,-1,-1,&
       &0,1,2,6,-1,-1,-1,-1,&
       &3,6,-1,-1,-1,-1,-1,-1,&
       &0,3,6,-1,-1,-1,-1,-1,&
       &1,3,6,-1,-1,-1,-1,-1,&
       &0,1,3,6,-1,-1,-1,-1,&
       &2,3,6,-1,-1,-1,-1,-1,&
       &0,2,3,6,-1,-1,-1,-1,&
       &1,2,3,6,-1,-1,-1,-1,&
       &0,1,2,3,6,-1,-1,-1,&
       &4,6,-1,-1,-1,-1,-1,-1,&
       &0,4,6,-1,-1,-1,-1,-1,&
       &1,4,6,-1,-1,-1,-1,-1,&
       &0,1,4,6,-1,-1,-1,-1,&
       &2,4,6,-1,-1,-1,-1,-1,&
       &0,2,4,6,-1,-1,-1,-1,&
       &1,2,4,6,-1,-1,-1,-1,&
       &0,1,2,4,6,-1,-1,-1,&
       &3,4,6,-1,-1,-1,-1,-1,&
       &0,3,4,6,-1,-1,-1,-1,&
       &1,3,4,6,-1,-1,-1,-1,&
       &0,1,3,4,6,-1,-1,-1,&
       &2,3,4,6,-1,-1,-1,-1,&
       &0,2,3,4,6,-1,-1,-1,&
       &1,2,3,4,6,-1,-1,-1,&
       &0,1,2,3,4,6,-1,-1,&
       &5,6,-1,-1,-1,-1,-1,-1,&
       &0,5,6,-1,-1,-1,-1,-1,&
       &1,5,6,-1,-1,-1,-1,-1,&
       &0,1,5,6,-1,-1,-1,-1,&
       &2,5,6,-1,-1,-1,-1,-1,&
       &0,2,5,6,-1,-1,-1,-1,&
       &1,2,5,6,-1,-1,-1,-1,&
       &0,1,2,5,6,-1,-1,-1,&
       &3,5,6,-1,-1,-1,-1,-1,&
       &0,3,5,6,-1,-1,-1,-1,&
       &1,3,5,6,-1,-1,-1,-1,&
       &0,1,3,5,6,-1,-1,-1,&
       &2,3,5,6,-1,-1,-1,-1,&
       &0,2,3,5,6,-1,-1,-1,&
       &1,2,3,5,6,-1,-1,-1,&
       &0,1,2,3,5,6,-1,-1,&
       &4,5,6,-1,-1,-1,-1,-1,&
       &0,4,5,6,-1,-1,-1,-1,&
       &1,4,5,6,-1,-1,-1,-1,&
       &0,1,4,5,6,-1,-1,-1,&
       &2,4,5,6,-1,-1,-1,-1,&
       &0,2,4,5,6,-1,-1,-1,&
       &1,2,4,5,6,-1,-1,-1,&
       &0,1,2,4,5,6,-1,-1,&
       &3,4,5,6,-1,-1,-1,-1,&
       &0,3,4,5,6,-1,-1,-1,&
       &1,3,4,5,6,-1,-1,-1,&
       &0,1,3,4,5,6,-1,-1,&
       &2,3,4,5,6,-1,-1,-1,&
       &0,2,3,4,5,6,-1,-1,&
       &1,2,3,4,5,6,-1,-1,&
       &0,1,2,3,4,5,6,-1,&
       &7,-1,-1,-1,-1,-1,-1,-1,&
       &0,7,-1,-1,-1,-1,-1,-1,&
       &1,7,-1,-1,-1,-1,-1,-1,&
       &0,1,7,-1,-1,-1,-1,-1,&
       &2,7,-1,-1,-1,-1,-1,-1,&
       &0,2,7,-1,-1,-1,-1,-1,&
       &1,2,7,-1,-1,-1,-1,-1,&
       &0,1,2,7,-1,-1,-1,-1,&
       &3,7,-1,-1,-1,-1,-1,-1,&
       &0,3,7,-1,-1,-1,-1,-1,&
       &1,3,7,-1,-1,-1,-1,-1,&
       &0,1,3,7,-1,-1,-1,-1,&
       &2,3,7,-1,-1,-1,-1,-1,&
       &0,2,3,7,-1,-1,-1,-1,&
       &1,2,3,7,-1,-1,-1,-1,&
       &0,1,2,3,7,-1,-1,-1,&
       &4,7,-1,-1,-1,-1,-1,-1,&
       &0,4,7,-1,-1,-1,-1,-1,&
       &1,4,7,-1,-1,-1,-1,-1,&
       &0,1,4,7,-1,-1,-1,-1,&
       &2,4,7,-1,-1,-1,-1,-1,&
       &0,2,4,7,-1,-1,-1,-1,&
       &1,2,4,7,-1,-1,-1,-1,&
       &0,1,2,4,7,-1,-1,-1,&
       &3,4,7,-1,-1,-1,-1,-1,&
       &0,3,4,7,-1,-1,-1,-1,&
       &1,3,4,7,-1,-1,-1,-1,&
       &0,1,3,4,7,-1,-1,-1,&
       &2,3,4,7,-1,-1,-1,-1,&
       &0,2,3,4,7,-1,-1,-1,&
       &1,2,3,4,7,-1,-1,-1,&
       &0,1,2,3,4,7,-1,-1,&
       &5,7,-1,-1,-1,-1,-1,-1,&
       &0,5,7,-1,-1,-1,-1,-1,&
       &1,5,7,-1,-1,-1,-1,-1,&
       &0,1,5,7,-1,-1,-1,-1,&
       &2,5,7,-1,-1,-1,-1,-1,&
       &0,2,5,7,-1,-1,-1,-1,&
       &1,2,5,7,-1,-1,-1,-1,&
       &0,1,2,5,7,-1,-1,-1,&
       &3,5,7,-1,-1,-1,-1,-1,&
       &0,3,5,7,-1,-1,-1,-1,&
       &1,3,5,7,-1,-1,-1,-1,&
       &0,1,3,5,7,-1,-1,-1,&
       &2,3,5,7,-1,-1,-1,-1,&
       &0,2,3,5,7,-1,-1,-1,&
       &1,2,3,5,7,-1,-1,-1,&
       &0,1,2,3,5,7,-1,-1,&
       &4,5,7,-1,-1,-1,-1,-1,&
       &0,4,5,7,-1,-1,-1,-1,&
       &1,4,5,7,-1,-1,-1,-1,&
       &0,1,4,5,7,-1,-1,-1,&
       &2,4,5,7,-1,-1,-1,-1,&
       &0,2,4,5,7,-1,-1,-1,&
       &1,2,4,5,7,-1,-1,-1,&
       &0,1,2,4,5,7,-1,-1,&
       &3,4,5,7,-1,-1,-1,-1,&
       &0,3,4,5,7,-1,-1,-1,&
       &1,3,4,5,7,-1,-1,-1,&
       &0,1,3,4,5,7,-1,-1,&
       &2,3,4,5,7,-1,-1,-1,&
       &0,2,3,4,5,7,-1,-1,&
       &1,2,3,4,5,7,-1,-1,&
       &0,1,2,3,4,5,7,-1,&
       &6,7,-1,-1,-1,-1,-1,-1,&
       &0,6,7,-1,-1,-1,-1,-1,&
       &1,6,7,-1,-1,-1,-1,-1,&
       &0,1,6,7,-1,-1,-1,-1,&
       &2,6,7,-1,-1,-1,-1,-1,&
       &0,2,6,7,-1,-1,-1,-1,&
       &1,2,6,7,-1,-1,-1,-1,&
       &0,1,2,6,7,-1,-1,-1,&
       &3,6,7,-1,-1,-1,-1,-1,&
       &0,3,6,7,-1,-1,-1,-1,&
       &1,3,6,7,-1,-1,-1,-1,&
       &0,1,3,6,7,-1,-1,-1,&
       &2,3,6,7,-1,-1,-1,-1,&
       &0,2,3,6,7,-1,-1,-1,&
       &1,2,3,6,7,-1,-1,-1,&
       &0,1,2,3,6,7,-1,-1,&
       &4,6,7,-1,-1,-1,-1,-1,&
       &0,4,6,7,-1,-1,-1,-1,&
       &1,4,6,7,-1,-1,-1,-1,&
       &0,1,4,6,7,-1,-1,-1,&
       &2,4,6,7,-1,-1,-1,-1,&
       &0,2,4,6,7,-1,-1,-1,&
       &1,2,4,6,7,-1,-1,-1,&
       &0,1,2,4,6,7,-1,-1,&
       &3,4,6,7,-1,-1,-1,-1,&
       &0,3,4,6,7,-1,-1,-1,&
       &1,3,4,6,7,-1,-1,-1,&
       &0,1,3,4,6,7,-1,-1,&
       &2,3,4,6,7,-1,-1,-1,&
       &0,2,3,4,6,7,-1,-1,&
       &1,2,3,4,6,7,-1,-1,&
       &0,1,2,3,4,6,7,-1,&
       &5,6,7,-1,-1,-1,-1,-1,&
       &0,5,6,7,-1,-1,-1,-1,&
       &1,5,6,7,-1,-1,-1,-1,&
       &0,1,5,6,7,-1,-1,-1,&
       &2,5,6,7,-1,-1,-1,-1,&
       &0,2,5,6,7,-1,-1,-1,&
       &1,2,5,6,7,-1,-1,-1,&
       &0,1,2,5,6,7,-1,-1,&
       &3,5,6,7,-1,-1,-1,-1,&
       &0,3,5,6,7,-1,-1,-1,&
       &1,3,5,6,7,-1,-1,-1,&
       &0,1,3,5,6,7,-1,-1,&
       &2,3,5,6,7,-1,-1,-1,&
       &0,2,3,5,6,7,-1,-1,&
       &1,2,3,5,6,7,-1,-1,&
       &0,1,2,3,5,6,7,-1,&
       &4,5,6,7,-1,-1,-1,-1,&
       &0,4,5,6,7,-1,-1,-1,&
       &1,4,5,6,7,-1,-1,-1,&
       &0,1,4,5,6,7,-1,-1,&
       &2,4,5,6,7,-1,-1,-1,&
       &0,2,4,5,6,7,-1,-1,&
       &1,2,4,5,6,7,-1,-1,&
       &0,1,2,4,5,6,7,-1,&
       &3,4,5,6,7,-1,-1,-1,&
       &0,3,4,5,6,7,-1,-1,&
       &1,3,4,5,6,7,-1,-1,&
       &0,1,3,4,5,6,7,-1,&
       &2,3,4,5,6,7,-1,-1,&
       &0,2,3,4,5,6,7,-1,&
       &1,2,3,4,5,6,7,-1,&
       &0,1,2,3,4,5,6,7/)
  !
contains
  !
  !****f* src/module/packb
  ! NAME
  !
  !  Function packb
  !
  ! USAGE
  !
  !  integer = packb(set)
  !
  ! DESCRIPTION
  !
  !  This function transforms a set of integers into
  !  an integer, this integer is unique
  ! Apparently Fortran allows to use arrays
  ! for the second argument which saves us a loop.
  !
  ! The elements in set have to be <= 31 which should
  ! not be a problem for realistic applications.
  !
  ! INPUTS
  !
  !  * set -- a set of integer
  !
  ! SIDE EFFECTS
  !
  !  No side effect (pure function)
  !
  ! RETURN VALUE
  !
  !  an integer 
  !
  ! EXAMPLE
  !
  !  i = packb( (/1,2,3/) ) 
  !  i is 14 which is in binary base 1110
  !
  !*****
  !
  pure function packb(set) result(bits)
    !
    integer, intent(in), dimension(:) :: set
    integer :: bits
    !
    bits = sum( ibset(0,pos=set) )
    !
  end function packb
  !
  !****f* src/module/unpackb
  ! NAME
  !
  !  Function unpackb
  !
  ! USAGE
  !
  !  integer_set = unpackb(bits,dim)
  !
  ! DESCRIPTION
  !
  !  This function performs the inverse operation
  !  as packb does : from an integer, it reconstructs the
  !  set of integers
  !
  ! INPUTS
  !
  !  * bits -- an integer
  !  * dim -- an integer, the dimension of the set obtained
  !
  ! SIDE EFFECTS
  !
  !  No side effect (pure function)
  !
  ! RETURN VALUE
  !
  !  an integer array of rank 1 and shape dim
  !
  ! EXAMPLE
  !
  !  set = unpackb( 14 ) 
  !  set is (/1,2,3/) because the binary representation of 14 is 1110
  !
  !*****
  !
  pure function unpackb(bits,dim)
    !
    integer, intent(in) :: bits
    integer, intent(in) :: dim
    integer, dimension(dim) :: unpackb
    !
    integer :: i,k,n
    !
    if (bits < 256) then
       !
       n = bits * 8
       unpackb = bit_sets(n:n + (dim-1))
       !
    else
       !
       i = bits
       k = 0
       n = 1
       !
       do while (i /= 0)
          !
          if (modulo(i,2) == 1)  then
             !
             unpackb(n) = k
             n = n + 1
             !
          end if
          !
          k = k+1
          i = ishft(i,-1)
          !
       end do
       !
    end if
    !
  end function unpackb
  !
  !****f* src/module/pminus
  ! NAME
  !
  !  Function pminus
  !
  ! USAGE
  !
  !  integer = pminus(bits1,bits2)
  !
  ! DESCRIPTION
  !
  !  This function subtracts the set which is 
  !  represented by bits2 to the one that is 
  !  represented by bits1. If the two sets set1 and set2
  !  are defined by set1=unpackb(bits1,dim1)
  !  and set2=unpackb(bits2,dim2), then ib = pminus(bits1,bits2)
  !  gives an integer such that unpackb(ib,dim_ib) is the set
  !  of integers of shape dim1-dim2 (dim1 > dim2) which contains
  !  the elements of set1 which do not belong to set2
  !  Note that if dim1 < dim2, the result returns is pminus(bits2,bits1)
  !  If none of the elements of set2 belongs to set1, then 
  !  pminus(bits1,bits2) = bits1
  !
  ! INPUTS
  !
  !  * bits1 -- an integer
  !  * bits2 -- an integer
  !
  ! SIDE EFFECTS
  !
  !  No side effect (pure function)
  !
  ! RETURN VALUE
  !
  !  an integer 
  !
  ! EXAMPLE
  !
  !  i1 = packb( (/1,2,3/) ) 
  !  i2 = packb( (/2/) )
  !  i3 = pminus(i1,i2)
  !  unpackb(i3) is the set (/1,3/)
  !
  !*****
  !
  pure function pminus(bits1, bits2) result(bits)
    !
    integer, intent(in) :: bits1, bits2
    integer :: bits
    !
    integer :: cits1, cits2
    !
    if ( bits1 >= bits2 ) then
       !
       cits1 = bits1
       cits2 = bits2
       !
    else
       !
       cits1 = bits2
       cits2 = bits1
       !
    end if
    !
    bits = iand(cits1,not(cits2))
    !
  end function pminus
  !
  !****f* src/module/punion
  ! NAME
  !
  !  Function punion
  !
  ! USAGE
  !
  !  integer = punion(bits1,bits2)
  !
  ! DESCRIPTION
  !
  !  This function adds the set which is 
  !  represented by bits2 to the one that is 
  !  represented by bits1. If the two sets set1 and set2
  !  are defined by set1=unpackb(bits1,dim1)
  !  and set2=unpackb(bits2,dim2), then ib = punion(bits1,bits2)
  !  gives an integer such that unpackb(ib,dim_ib) is the set
  !  of integers of shape dim1+dim2 which contains
  !  the elements of set1 and those of set2
  !  Note that if some elements of set2 belong to set1, they do not
  !  appear twice
  !
  ! INPUTS
  !
  !  * bits1 -- an integer
  !  * bits2 -- an integer
  !
  ! SIDE EFFECTS
  !
  !  No side effect (pure function)
  !
  ! RETURN VALUE
  !
  !  an integer 
  !
  ! EXAMPLE
  !
  !  i1 = packb( (/1,3,4/) ) 
  !  i2 = packb( (/2/) )
  !  i3 = punion(i1,i2)
  !  unpackb(i3) is the set (/1,2,3,4/)
  !
  !*****
  !
  pure function punion(bits1, bits2) result(bits)
    !
    integer, intent(in) :: bits1, bits2
    integer :: bits
    !
    bits = ior(bits1,bits2)
    !
  end function punion
  !
  !
  !****f* src/module/countb
  ! NAME
  !
  !  Function countb
  !
  ! USAGE
  !
  !  integer = countb(bits)
  !
  ! DESCRIPTION
  !
  !  This function computes the shape of the rank 1 integer 
  !  set given by unpackb(bits,dim)
  !
  ! INPUTS
  !
  !  * bits -- an integer
  !
  ! SIDE EFFECTS
  !
  !  No side effect (pure function)
  !
  ! RETURN VALUE
  !
  !  an integer 
  !
  ! EXAMPLE
  !
  !  i1 = packb( (/1,2,3/) )
  !  i2 = countb(i1)
  !  i2 is 3
  !
  !*****
  !
  pure function countb(bits)
    !
    integer, intent(in) :: bits
    integer :: countb
    !
    integer :: i
    !
    if (bits < 256) then
       !
       countb = bit_count(bits)
       !
    else
       !
       countb = 0
       i = bits
       do while (i /= 0)
          !
          countb = countb + bit_count(iand(i, 255))
          i = ishft(i,-8)
          !
       end do
       !
    end if
    !
  end function countb
  !
  !
  !****f* src/module/locateb
  ! NAME
  !
  !  Function locateb
  !
  ! USAGE
  !
  !  integer = locateb(i,bits)
  !
  ! DESCRIPTION
  !
  !  The function locateb returns the location of the element i  
  !  in the set given by unpackb(bits,countb(bits)). 
  !  If i does not belong to bits, the function locateb 
  !  returns -1
  !  
  !
  ! INPUTS
  !
  !  * i -- an integer, the element of a
  !  * bits -- an integer
  !
  ! SIDE EFFECTS
  !
  !  No side effect (elemental)
  !
  ! RETURN VALUE
  !
  !  It returns an integer, the location of i in the array a
  !
  ! EXAMPLE
  !
  !  bits = packb( (/3,5,6,7/) )
  !  j = locateb(5,bits) --> j is equal to 2
  !  j = locateb(6,bits) --> j is equal to 3
  !  j = locateb(4,bits) --> j is equal to -1
  !  Note that if the set is not ordered, the packing
  !  orders it.
  !  Note also that this function has the attribute elemental
  !  that means that, the argument can be a set of integers:
  !  locateb( (/3,7/) , bits) will return (/1,4/)
  !  
  !
  !*****
  !
  elemental function locateb(i,bits)
    !
    integer, intent(in) :: bits
    integer, intent(in) :: i
    integer :: locateb
    !
    integer :: ib
    !
    if (btest(bits,i)) then
       !
       ib = ibits(bits,0,i)
       if (ib < 256) then
          !
          locateb = bit_count(ib)+1
          !
       else
          !
          locateb = countb(ib)+1
          !
       end if
       !
    else
       !
       locateb = -1
       !
    end if
    !
  end function locateb
  !
end module array
