! 
!****h* src/module/cache_generic
! NAME
!
!  Module cache_generic
!
! USAGE
!
!  use cache_generic
!
! DESCRIPTION
!
!  This module is used to reserve some memory to store already computed results
!
! OUTPUT
!
!  This module exports three routines:
!  * allocate_cache -- to reserve the memory
!  * reset_cache -- to force the re-computation of the cache arrays
!  * clear_cache -- to clear the reserved memory
!
!
!
!*****
module cache_generic
  !
  use precision_golem
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use form_factor_type
  implicit none
  !
  private :: ki
  !

  integer(kind=8), dimension(10240):: cache_generic_tag;
  integer,dimension(10240) :: cache_generic_next;
  integer,dimension(1024) :: cache_table;
  type(form_factor),dimension(10240) ::  cache_generic_val;
  integer :: cache_generic_count

  public :: allocate_cache_generic, clear_cache_generic, reset_cache_generic
  public :: cache_generic_get_value, cache_generic_put_value

  contains
    subroutine allocate_cache_generic()
     !
     call reset_cache_generic()
     !
    end subroutine allocate_cache_generic

    subroutine reset_cache_generic
      cache_generic_count = 0
      cache_table = 0
    end subroutine reset_cache_generic

    subroutine clear_cache_generic()
          call reset_cache_generic()
          !
    end subroutine clear_cache_generic

    pure integer(kind=8) function get_hash_value(leg_count,dim_nplus,b_pin,l_count,l) result(hash)
       implicit none
       integer, intent(in) :: leg_count
       integer, intent(in) :: dim_nplus
       integer, intent (in) :: b_pin
       integer, intent (in) :: l_count
       integer, intent (in),dimension(:) :: l
       integer :: i

       hash = 0 

       do i = 1, l_count
          hash = hash*7 + l(i)
       end do
       hash = hash*128 + b_pin
       hash = hash*30 + dim_nplus
       hash = hash*7 + leg_count

    end function get_hash_value



    logical function cache_generic_get_value(leg_count,dim_nplus,b_pin,l_count,l,val) result(res)
       implicit none
       integer, intent(in) :: leg_count
       integer, intent(in) :: dim_nplus
       integer, intent (in) :: b_pin
       integer, intent (in) :: l_count
       integer, intent (in),dimension(:) :: l


       type(form_factor), intent(out) :: val

       integer(kind=8) :: hash_val
       integer :: hash_index, i, j, k;
       logical :: found

       hash_val=get_hash_value(leg_count,dim_nplus,b_pin,l_count,l)
       hash_index=mod(hash_val, 1021) + 1
       i = cache_table(hash_index)
       found = .false.
       do while ((i/=0) .and. (.not. found))
               found = cache_generic_tag(i)==hash_val
               j=i
               i = cache_generic_next(i)
       end do
       if (found) then
               val = cache_generic_val(j)
               res=.true.
               return
       end if
       res=.false.
       return
    end function cache_generic_get_value

    subroutine cache_generic_put_value(leg_count,dim_nplus,b_pin,l_count,l,val)
       implicit none
       integer, intent(in) :: leg_count
       integer, intent(in) :: dim_nplus
       integer, intent (in) :: b_pin
       integer, intent (in) :: l_count
       integer, intent (in),dimension(:) :: l

       type(form_factor), intent(in) :: val

       integer(kind=8) :: hash_val
       integer :: hash_index, i, j, k;
       logical :: found

       hash_val=get_hash_value(leg_count,dim_nplus,b_pin,l_count,l);
       hash_index=mod(hash_val, 1021) + 1
       i = cache_table(hash_index)
       found = .false.
       do while ((i/=0) .and. (.not. found))
               found = cache_generic_tag(i)==hash_val
               j=i
               i=cache_generic_next(i)
       enddo
       if (found) then
               cache_generic_val(j)=val
       else
               cache_generic_count = cache_generic_count+1
               !if (cache_table(hash_index) /= 0) then
               !        write (*, *) "collision", cache_generic_tag(cache_table(hash_index))
               !        write (*, *) "current", hash_val, leg_count,dim_nplus,b_pin,l_count,l,val
               !end if
               if (cache_generic_count>10240) then
                       return
               end if
               cache_generic_val(cache_generic_count) = val
               cache_generic_tag(cache_generic_count) = hash_val
               cache_generic_next(cache_generic_count) = cache_table(hash_index)
               cache_table(hash_index) = cache_generic_count
       end if
       !write (*, *) cache_generic_count
       
    end subroutine cache_generic_put_value

end module cache_generic
