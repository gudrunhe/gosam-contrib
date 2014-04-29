! 
!****h* src/module/cache
! NAME
!
!  Module cache
!
! USAGE
!
!  use cache
!
! DESCRIPTION
!
!  This module is used to reserve some memory to store already computed four/three
!  point functions
!
! OUTPUT
!
!  This module exports three routines:
!  * allocate_cache -- to reserve the memory
!  * reset_cache -- to force the re-computation of the cache arrays
!  * clear_cache -- to clear the reserved memory
!
! USES
!
!  * sortie_erreur (src/module/sortie_erreur.f90)
!
!
!*****
module cache
  !
  use precision_golem
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use form_factor_type
  use cache_generic
  implicit none
  !
  private :: ki
  !
  integer, private :: err
  logical, dimension(:,:,:,:,:), allocatable :: computed_f4p_np2
  complex(ki),dimension(:,:,:,:,:), allocatable :: results_f4p_np2
  logical, dimension(:,:,:), allocatable :: computed_f4p_np4
  complex(ki),dimension(:,:,:,:), allocatable :: results_f4p_np4
  logical, dimension(:,:,:,:,:,:), allocatable :: computed_f3p
  real(ki),dimension(:,:,:,:,:,:,:), allocatable :: results_f3p
  logical, dimension(:,:,:,:,:,:), allocatable :: computed_f3p_np2
  real(ki),dimension(:,:,:,:,:,:,:), allocatable :: results_f3p_np2


  !
  ! everything public except err, ki
  !public :: allocate_cache, reset_cache, clear_cache
  !
  contains
    !
    !****f* src/module/cache/allocate_cache
    ! NAME
    !
    !  Subroutine allocate_cache
    !
    ! USAGE
    !
    !  call allocate_cache(dim_s)
    !
    ! DESCRIPTION
    !
    !  This subroutine allocates the necessary memory to store
    !  the n+2/n+4 four point functions and the n/n+2 three point functions
    !
    ! INPUTS
    !
    !  * dim_s -- an integer, the dimension of the S matrix
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !*****
    subroutine allocate_cache(dim_s)
      !
      integer, intent(in) :: dim_s
      !
      allocate(computed_f4p_np2(0:dim_s,0:dim_s,0:dim_s,0:dim_s,0:dim_s),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocate_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for computed_f4p_np2'
        call catch_exception(0)
        !
      end if
      !
      allocate(computed_f4p_np4(0:dim_s,0:dim_s,0:dim_s),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocate_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for computed_f4p_np4'
        call catch_exception(0)
        !
      end if
      !
      allocate(computed_f3p(0:dim_s,0:dim_s,0:dim_s,0:3,0:3,0:3),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocate_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for computed_f3p'
        call catch_exception(0)
        !
      end if
      !
      allocate(computed_f3p_np2(0:dim_s,0:dim_s,0:dim_s,0:3,0:3,0:3),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocate_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for computed_f3p_np2'
        call catch_exception(0)
        !
      end if
      computed_f4p_np2 = .false.
      computed_f4p_np4 = .false.
      computed_f3p = .false.
      computed_f3p_np2 = .false.
      !
      allocate(results_f4p_np2(0:dim_s,0:dim_s,0:dim_s,0:dim_s,0:dim_s),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocate_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for results_f4p_np2'
        call catch_exception(0)
        !
      end if
      !
      allocate(results_f4p_np4(0:dim_s,0:dim_s,0:dim_s,2),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocate_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for computed_f4p_np4'
        call catch_exception(0)
        !
      end if
      !
      allocate(results_f3p(0:dim_s,0:dim_s,0:dim_s,0:3,0:3,0:3,6),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocate_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for results_f3p'
        call catch_exception(0)
        !
      end if
      !
      allocate(results_f3p_np2(0:dim_s,0:dim_s,0:dim_s,0:3,0:3,0:3,4),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocate_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for results_f3p_np2'
        call catch_exception(0)
        !
      end if
      !
      call allocate_cache_generic()
      !
    end subroutine allocate_cache
    !
    !
    !****f* src/module/cache/reset_cache
    ! NAME
    !
    !  Subroutine reset_cache
    !
    ! USAGE
    !
    !  call reset_cache()
    !
    ! DESCRIPTION
    !
    !  This subroutine forces the cache arrays to be computed again 
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
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !*****
    subroutine reset_cache()
      !
      computed_f4p_np2 = .false.
      computed_f4p_np4 = .false.
      computed_f3p = .false.
      computed_f3p_np2 = .false.

      call reset_cache_generic()
      !
    end subroutine reset_cache
    !
    !
    !
    !****f* src/module/cache/clear_cache
    ! NAME
    !
    !  Subroutine clear_cache
    !
    ! USAGE
    !
    !  call clear_cache()
    !
    ! DESCRIPTION
    !
    !  This subroutine deallocates the reserved memory to store
    !  the n+2/n+4 four point functions and the n/n+2 three point functions.
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
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !*****
    subroutine clear_cache()
      !
      cache_generic_count=0
      !
      deallocate(computed_f4p_np2,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine clear_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for computed_f4p_np2'
        call catch_exception(0)
        !
      end if
      !
      deallocate(results_f4p_np2,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine clear_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for results_f4p_np2'
        call catch_exception(0)
        !
      end if
      !
      deallocate(computed_f4p_np4,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine clear_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for computed_f4p_np4'
        call catch_exception(0)
        !
      end if
      !
      deallocate(results_f4p_np4,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine clear_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for results_f4p_np4'
        call catch_exception(0)
        !
      end if
      !
      deallocate(computed_f3p,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine clear_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for computed_f3p'
        call catch_exception(0)
        !
      end if
      !
      deallocate(results_f3p,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine clear_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for results_f3p'
        call catch_exception(0)
        !
      end if
      !
      deallocate(computed_f3p_np2,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine clear_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for computed_f3p_np2'
        call catch_exception(0)
        !
      end if
      !
      deallocate(results_f3p_np2,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine clear_cache'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for results_f3p_np2'
        call catch_exception(0)
        !
      end if
      !
      call clear_cache_generic()
      !
    end subroutine clear_cache
end module cache
