! 
!****h* src/module/tri_croissant
! NAME
!
!  Module tri_croissant
!
! USAGE
!
!  use tri_croissant
!
! DESCRIPTION
!
!  This module is used to sort an integer array or shift its elements by certain amount
!
! OUTPUT
!
!  This module exports:
!  * tri_int -- a subroutine to sort out an integer array
!  * shift_param -- a subroutine to shift (modulo n) the elements of an integer array
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!
!*****
module tri_croissant
  !
  use precision_golem
  implicit none
  !
  private
  !
  !
  public :: tri_int2, tri_int3, tri_int4, shift_param, exchange_param
  !
  !
  contains
    !
    !****f* src/module/tri_croissant/tri_int
    ! NAME
    !
    !  Subroutine tri_int
    !
    ! USAGE
    !
    !  call tri_int(t_in,t_out)
    !
    ! DESCRIPTION
    !
    !  This routine sorts in increasing order an integer array t_int and put the 
    !  result in the integer array t_out
    !
    ! INPUTS
    !
    !  * t_int -- an integer array of rank 1, the array to sort
    !
    ! SIDE EFFECTS
    !
    !     NONE
    !
    ! RETURN VALUE
    !
    !  * t_out -- an integer array of rank 1, the sorted array
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    ! on trie par ordre croissant le tableau t_in et 
    ! on met le resultat dans t_out. On utilise les procedures 
    ! compare et exchange

    pure subroutine tri_int2(t_in,t_out)
       implicit none
       integer, intent(in), dimension(2) :: t_in
       integer, intent(out), dimension(2) :: t_out
       integer, dimension(2) :: tmp

       if (t_in(1) < t_in(2)) then
          tmp(1) = t_in(1)
          tmp(2) = t_in(2)
       else
          tmp(1) = t_in(2)
          tmp(2) = t_in(1)
       end if
       t_out(:) = tmp(:)
    end  subroutine tri_int2

    pure subroutine tri_int3(t_in,t_out)
       implicit none
       integer, intent(in), dimension(3) :: t_in
       integer, intent(out), dimension(3) :: t_out
       integer, dimension(3) :: tmp

       tmp(:) = t_in(:)

       if (tmp(1) > tmp(3)) then
          if (tmp(1) < tmp(2)) then
             t_out(:) = tmp((/3,1,2/))
          else
             if (tmp(2) > tmp(3)) then
                t_out(:) = tmp((/3,2,1/))
             else
                t_out(:) = tmp((/2,3,1/))
             end if
          end if
       else
          if (tmp(2) < tmp(1)) then
             t_out(:) = tmp((/2,1,3/))
          else
             if (tmp(2) < tmp(3)) then
                t_out(:) = tmp((/1,2,3/))
             else
                t_out(:) = tmp((/1,3,2/))
             end if
          end if
       end if
    end  subroutine tri_int3

    subroutine tri_int4(t_in,t_out)
      implicit none
      !
      integer, intent(in), dimension(4) :: t_in
      integer, intent(out), dimension(4) :: t_out
      !
      integer :: j,i,alpha,beta
      !

      t_out(1) = t_in(1)
      do i=2,4
         alpha = t_in(i)

         do j=1,i-1
            if (alpha < t_out(j)) then
               beta = t_out(j)
               t_out(j) = alpha
               alpha = beta
            end if
         end do
         t_out(i) = alpha
      end do
    end subroutine tri_int4
    !
    !****f* src/module/tri_croissant/shift_param
    ! NAME
    !
    !  Subroutine shift_param
    !
    ! USAGE
    !
    !  call shift_param(z_param_ini,shift,modd,z_param_out)
    !
    ! DESCRIPTION
    !
    !  This routine shifts the array z_param_ini of Feynman parameters
    !  the shift is done as following: 
    !  z --> z+shift if z+shift <= modd 
    !  else mod(z+shift-1,modd)+1 if z+shift > modd
    !  the result is put into the array z_param_out
    !
    ! INPUTS
    !
    !  * z_param_ini -- an integer array of rank 1, the array to shift
    !  * shift -- an integer, the value of the shift
    !  * modd -- an integer, the shift is made modulo modd
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * z_param_out -- an integer array of rank 1, the shifted array
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine shift_param(z_param_ini,shift,modd,z_param_out)
      !
      integer, intent(in) :: shift,modd
      integer, intent(in), dimension(modd) :: z_param_ini
      integer, intent(out), dimension(modd) :: z_param_out
      !
      integer, dimension(modd) :: temp
      !
      where (z_param_ini .ne. 0)
         temp(:) = modulo(z_param_ini(:)+shift-1, modd) + 1
      elsewhere
         temp(:) = 0
      endwhere
      !
      select case(modd)
         case(1)
            z_param_out = temp
         case(2)
            call tri_int2(temp,z_param_out)
         case(3)
            call tri_int3(temp,z_param_out)
         case(4)
            call tri_int4(temp,z_param_out)
         case default
            print*, "shift_param: unimplemented value of modd: ", modd
            stop
      end select
      !
    end subroutine shift_param
    !
    !
    !****f* src/module/tri_croissant/exchange_param
    ! NAME
    !
    !  Subroutine exchange_param
    !
    ! USAGE
    !
    !  call exchange_param(z_param_ini,tab,modd,z_param_out)
    !
    ! DESCRIPTION
    !
    !  This routine exchanges in the array z_param_ini of Feynman parameters
    !  the label tab(1) with the label tab(2)
    !  the result is put into the array z_param_out
    !
    ! INPUTS
    !
    !  * z_param_ini -- an integer array of rank 1, the array to shift
    !  * tab -- an integer array of rank 1, the two labels to exchange
    !  * modd -- an integer, the shift is made modulo modd
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * z_param_out -- an integer array of rank 1, the shifted array
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine exchange_param(z_param_ini,tab,modd,z_param_out)
      !
      integer, intent(in) :: modd
      integer, dimension(2), intent(in) :: tab
      integer, intent(in), dimension(modd) :: z_param_ini
      integer, intent(out), dimension(modd) :: z_param_out
      !
      integer, dimension(modd) :: temp
      !
      where (z_param_ini == tab(1))
         temp = tab(2)
      elsewhere (z_param_ini == tab(2))
         temp = tab(1)
      elsewhere
         temp = z_param_ini
      end where
      !
      select case(modd)
         case(1)
            z_param_out = temp
         case(2)
            call tri_int2(temp,z_param_out)
         case(3)
            call tri_int3(temp,z_param_out)
         case(4)
            call tri_int4(temp,z_param_out)
         case default
            print*, "shift_param: modd too large: ", modd
            stop
      end select
      !
    end subroutine exchange_param
    !
end module tri_croissant
!
