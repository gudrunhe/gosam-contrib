!
!****h* src/kinematic/matrice_s
! NAME
!
!  Module matrice_s
!
! USAGE
!
!  use matrice_s
!
! DESCRIPTION
!
!  This module is used : to reserve some memory in order to pass the S matrix, its
!  shape, the set of propagator labels; to compute the inverse
!  of S matrix and the related quantities : the b's and sumb, also for all
!  possible reduced matrices. The S matrix is allocated here and also its dimension
!  and it returns the result through the three functions inv_s, b and sumb.
!
!
! OUTPUT
!
!  This module exports five variables:
!  * dim_s -- an integer, the shape of the initial S matrix
!  * set_ref -- an integer array, the set of initial propagators
!  * s_mat_c -- a complex (type ki) array of rank 2, the S matrix.
!  * s_mat_p -- a derived type, including the S matrix for either real or complex masses
!               and integer-bits encoding the positions of masses with non-vanishing Im-part,
!               and vanishing masses.
!  * s_mat -- A pointer associated with s_mat_c. The user can fill s_mat or s_mat_c with complex or
!             real values.
!
!  and also eleven functions:
!  * initgolem95 -- calls allocation_s, initializes the cache, associates s_mat.
!  * allocation_s -- to allocate the required memory
!  * deallocation_s -- to deallocate the used memory
!  * preparesmatrix -- fill s_mat_r with the real part of s_mat_c, sets the bit integers in s_mat_p
!                        calls init_invs.
!  * init_invs -- to fill all the array for the inverse of the S matrix
!                 and the inverse of the reduce S matrix
!  * inv_s -- it contains the inverse of the S matrix
!  * hj -- it contains H matrix (pseudo-inverse of G)
!  * b -- it contains the b coefficients
!  * sumb -- it contains the B coeficient
!  * norma_sumb -- it contains the normalised B coefficient
!  * exitgolem95 -- deallocates memory, clear the cache.
!
!  Only dim_s and set_ref take a value in this module, not the other variables
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * cache (src/module/cache.f90)
!  * inverse_matrice (src/kinematic/inverse_matrice.f90)
!  * tri_croissant (src/module/tri_croissant.f90)
!  * array (src/module/array.f90)
!  * parametre (src/module/parametre.f90)
!  * s_matrix_type (src/module/s_matrix_type.f90)
!
!*****
!
module matrice_s
  !
  use precision_golem
  use sortie_erreur
  use cache
  use inverse_matrice
  use tri_croissant
  use array
  use parametre
  use s_matrix_type
  use constante, only:czero
  implicit none
  !
  !
  private
  !
  integer :: dim_s
  integer, dimension(:), allocatable :: set_ref
  integer, dimension(6) :: ref_vector = (/ 1, 2, 3, 4, 5, 6 /)
  integer :: b_ref
  real(ki), dimension(:,:), allocatable :: s_mat_r
  complex(ki), dimension(:,:), allocatable, target :: s_mat_c
  type(s_matrix_poly) :: s_mat_p
  complex(ki), dimension(:,:), pointer :: s_mat
  !
  public :: dim_s, set_ref, b_ref, s_mat_c, s_mat_p, s_mat
  !
  ! The first index of the following arrays is (b_pin/2)+1, the remaining indices
  ! are for the indices of the representing rank-n tensor.
  !
  real(ki), dimension(:,:,:), allocatable :: hjj_r
  complex(ki), dimension(:,:,:), allocatable :: hjj_c
  !
  real(ki), dimension(:,:,:), allocatable :: invs_n_r
  complex(ki), dimension(:,:,:), allocatable :: invs_n_c
  !
  real(ki), dimension(:,:), allocatable :: b_n_r
  complex(ki), dimension(:,:), allocatable :: b_n_c
  !
  real(ki), dimension(:), allocatable :: sumb_n_r
  complex(ki), dimension(:), allocatable :: sumb_n_c
  !
  real(ki), dimension(:), allocatable :: norma_sumb_n_r
  complex(ki), dimension(:), allocatable :: norma_sumb_n_c
  !
  integer :: err
  !
  public :: allocation_s, deallocation_s, init_invs, inv_s, hj, b, sumb, norma_sumb
  public :: initgolem95, preparesmatrix, prepare_s_matrix_local, exitgolem95
  !
  interface put_to_zero
     !
     module procedure put_to_zero_r, put_to_zero_c
     !
  end interface
  !
  contains
    !
    !****f* src/kinematic/matrice_s/initgolem95
    ! NAME
    !
    !  Subroutine initgolem95
    !
    ! USAGE
    !
    !  call initgolem95(dim, opt_set)
    !
    ! DESCRIPTION
    !
    !  This subroutine is the first of three macro functions which needs to be called by the user.
    !  It allocates memory for all internal matrices needed in subsequent calculations.
    !  The caching system is initialized.
    !  A pointer s_mat is associated with a complex matrix s_mat_c.
    !  This is the s matrix which has to be filled after initgolem95() is called.
    !  The argument 'dim' sets the maximal number of external legs.
    !  An optional array for the numbering of propagators can be given.
    !  The default is set to (/ 1, ... , dim /)
    !
    ! INPUTS
    !
    !  * dim -- an integer, the maximal number of external legs
    !  * opt_set -- an optional integer array for the numbering of propagators
    !
    ! SIDE EFFECTS
    !
    ! A call to allocation_s is made, implying all side effects given there.
    ! The caching system is initialized.
    ! A pointer 's_mat' is associated with the global matrix s_mat_c.
    ! The internal parameter rmass_or_cmass_par is set to cmass. If a purely real
    ! s matrix is given by the user it will be set to rmass in the call of
    ! preparesmatrix.
    !
    ! RETURN VALUE
    !
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine initgolem95(dim, opt_set)
      !
      integer, intent(in) :: dim
      integer, dimension(:), optional, intent(in) :: opt_set
      !
      rmass_or_cmass_par = cmass
      !
      call allocation_s(dim)
      !
      if (present(opt_set) ) then
         !
         set_ref(:) = opt_set(:)
         !
      else
         !
         set_ref(:) = ref_vector(1:dim)
         !
      end if
      !
      b_ref = packb(set_ref)
      !
      call allocate_cache(dim)
      !
      s_mat => s_mat_c
      !
    end subroutine initgolem95
    !
    !****f* src/kinematic/matrice_s/allocation_s
    ! NAME
    !
    !  Subroutine allocation_s
    !
    ! USAGE
    !
    !  call allocation_s(dim)
    !
    ! DESCRIPTION
    !
    !  This subroutine reserves the memory for several internal objects.
    !  In case of rmass_or_cmass_par==cmass, there complex copies of each preceding array
    !  are also allocated.
    !
    !  After memory allocation, s_mat_p is then assigned the matrix s_mat_c or s_mat_r,
    !  respectively. The corresponding pointers in s_mat_p are associated or nullified.
    !  In case a complex matrix is assigned, there will be also a pointer associated with
    !  a real matrix, which has entries according to the real part of the complex matrix.
    !
    ! INPUTS
    !
    !  * dim -- an integer, the maximal number of external legs
    !
    ! SIDE EFFECTS
    !
    !  This routine modify the value of the variable dim_s
    !  It initialises invs_n, hjj, b_n, sumb_n, norma_sumb_n to zero
    !  It associates the global objects s_mat_p with s_mat_r or s_mat_c.
    !
    ! RETURN VALUE
    !
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine allocation_s(dim)
      !
      integer, intent(in) :: dim
      !
      if (rmass_or_cmass_par%rmass_selected) then
         !
         call allocation_s_r(dim)
         s_mat_p = assign_s_matrix(s_mat_r)
         !
      else if (rmass_or_cmass_par%cmass_selected) then
         !
         call allocation_s_r(dim)
         call allocation_s_c(dim)
         s_mat_p = assign_s_matrix(s_mat_c,s_mat_r)
         !
      else
         !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'rmass_or_cmass_par has wrong value'
        call catch_exception(0)
         !
      end if
      !
    end subroutine allocation_s
    !
    subroutine allocation_s_r(dim)
      !
      integer, intent(in) :: dim
      integer :: err
      !
      dim_s = dim
      !
      allocate(s_mat_c(dim_s,dim_s),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for s_mat'
        call catch_exception(0)
        !
      end if
      !
      allocate(s_mat_r(dim_s,dim_s),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for s_mat'
        call catch_exception(0)
        !
      end if
      !
      allocate(set_ref(dim_s),stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot allocate memory for set_ref'
        call catch_exception(0)
        !
      end if
      if (dim>= 1 .and. dim <= 6 ) then
              allocate( invs_n_r(2**dim,dim,dim), &
                hjj_r(dim,dim,dim),b_n_r(2**dim,dim), sumb_n_r(2**dim), &
                norma_sumb_n_r(2**dim), stat=err)
              if (err /= 0) then
                      !
                      tab_erreur_par(1)%a_imprimer = .true.
                      tab_erreur_par(1)%chaine = 'In subroutine allocation_s'
                      tab_erreur_par(2)%a_imprimer = .true.
                      tab_erreur_par(2)%chaine = 'cannot allocate memory for invs_n_r or ...'
                      call catch_exception(0)
                      !
              end if
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine allocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'dimension %d0 not supported.'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
      invs_n_r = 0._ki
      b_n_r = 0._ki
      sumb_n_r = 0._ki
      norma_sumb_n_r = 0._ki
      !
    end subroutine allocation_s_r
    !
    subroutine allocation_s_c(dim)
      !
      integer, intent(in) :: dim
      integer :: err
      !
      dim_s = dim
      !
      !
      !
      if (dim>= 1 .and. dim <= 6 ) then
              allocate( invs_n_c(2**dim,dim,dim), &
                hjj_c(dim,dim,dim),b_n_c(2**dim,dim), sumb_n_c(2**dim), &
                norma_sumb_n_c(2**dim), stat=err)
              if (err /= 0) then
                      !
                      tab_erreur_par(1)%a_imprimer = .true.
                      tab_erreur_par(1)%chaine = 'In subroutine allocation_s'
                      tab_erreur_par(2)%a_imprimer = .true.
                      tab_erreur_par(2)%chaine = 'cannot allocate memory for invs_n_r or ...'
                      call catch_exception(0)
                      !
              end if
      end if
      !
      invs_n_c = czero
      b_n_c = czero
      sumb_n_c = czero
      norma_sumb_n_c = czero
      !
    end subroutine allocation_s_c
    !
    !****f* src/kinematic/matrice_s/preparesmatrix
    ! NAME
    !
    !  Subroutine preparesmatrix
    !
    ! USAGE
    !
    !  call preparesmatrix()
    !  call prepare_s_matrix_local(s_mat_p_loc,set_ref_loc)
    !
    ! DESCRIPTION
    !
    !  This subroutine prepares the global or local s_mat_p object, consisting
    !  of pointers to s_mat_c and s_mat_r and integer bits b_cmplx and b_zero.
    !  A call to init_invs is made to fill the inverse matrices needed
    !  in the form factor calculations.
    !  If the user has defined a purely real s matrix, the internal parameter
    !  rmass_or_cmass_par is set to rmass and only the real branch of the library
    !  is used.
    !  In the complex case, form factors which are not affected by complex
    !  masses will be called with a sub matrix of s_mat_r, the real part of s_mat_c.
    !  The routine also sets the bits for complex mass and zero mass-
    !  entries.
    !  The subroutine prepare_s_matrix_local is used internally to prepare local type
    !  s_matrix_poly objects. This subroutine does not interact with the inverse matrices
    !  and the caching system.
    !
    ! INPUTS
    !
    !  For prepare_s_matrix_local, s_mat_p and set_ref need to be given.
    !
    ! RETURN VALUE
    !
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine preparesmatrix()
      !
      rmass_or_cmass_par = cmass
      if (.not. associated(s_mat_p%pt_cmplx) ) s_mat_p%pt_cmplx => s_mat_c
      !
      call fill_s_matrix(s_mat_p)
      call set_s_matrix_bits(s_mat_p,set_ref)
      !
      if ( s_mat_p%b_cmplx == 0 ) then
         !
         rmass_or_cmass_par = rmass
         nullify(s_mat_p%pt_cmplx)
         !
      end if
      !
      call reset_cache()
      call init_invs()
      !
    end subroutine preparesmatrix
    !
    subroutine prepare_s_matrix_local(s_mat_poly,set_ref_loc)
      type(s_matrix_poly),intent (inout) :: s_mat_poly
      integer, dimension(:) :: set_ref_loc
      !
      call fill_s_matrix(s_mat_poly)
      call set_s_matrix_bits(s_mat_poly,set_ref_loc)
      !
    end subroutine prepare_s_matrix_local
    !
    !
    !****f* src/kinematic/matrice_s/deallocation_s
    ! NAME
    !
    !  Subroutine deallocation_s
    !
    ! USAGE
    !
    !  call deallocation_s()
    !
    ! DESCRIPTION
    !
    !  This subroutine deallocates the memory reserved by the preceeding
    !  subroutine.
    !  The pointers in s_mat_p are nullified.
    !
    ! INPUTS
    !
    !  No input
    !
    ! SIDE EFFECTS
    !
    !  This routine destroys all the variables initialised in the
    !  preceeding subroutine as well as any associations in s_mat_p.
    !
    ! RETURN VALUE
    !
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine deallocation_s()
      !
      call nullify_s_matrix(s_mat_p)
      !
      if (rmass_or_cmass_par%rmass_selected) then
         !
         call deallocation_s_r()
         !
      else if (rmass_or_cmass_par%cmass_selected) then
         !
         call deallocation_s_c()
         call deallocation_s_r()
         !
      else
         !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine deallocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'rmass_or_cmass has wrong value'
        call catch_exception(0)
         !
      end if
      !
    end subroutine deallocation_s
    !
    subroutine deallocation_s_r()
      !
      integer :: err
      !
      deallocate(s_mat_r,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine deallocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for s_mat'
        call catch_exception(0)
      end if
      !
      deallocate(set_ref,stat=err)
      !
      if (err /= 0) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine deallocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for s_ref'
        call catch_exception(0)
      end if
      !
      deallocate(invs_n_r,hjj_r,b_n_r, sumb_n_r, norma_sumb_n_r,stat=err)
      !
      if (err /= 0) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine deallocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for invs_n_r or hjj or b or sumb_n_r or norma_sumb_n_r'
        call catch_exception(0)
      end if
      !
    end subroutine deallocation_s_r
    !
    subroutine deallocation_s_c()
      !
      integer :: err
      !
      deallocate(s_mat_c,stat=err)
      !
      if (err /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine deallocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for s_mat'
        call catch_exception(0)
      end if
      !
      deallocate(invs_n_c,b_n_c,sumb_n_c,norma_sumb_n_c,hjj_c,stat=err)
      !
      if (err /= 0) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine deallocation_s'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'cannot deallocate memory for invs_n or b_n or sumb_n or norma_sumb_n or hjj_c'
        call catch_exception(0)
      end if
      !
    end subroutine deallocation_s_c
    !
    !****f* src/kinematic/inversion/init_invs
    ! NAME
    !
    !  Subroutine init_invs
    !
    ! USAGE
    !
    !  call init_invs()
    !
    ! DESCRIPTION
    !
    ! This function comes in two copies for real masses and complex masses.
    ! The respective arrays are filled.
    !
    ! This routine fills the arrays:
    ! invs_n, hjj, b_n, sumb_n, norma_sumb_n
    !
    ! One can print a typical error due to the numerical inversion
    !
    ! INPUTS
    !
    !  No input
    !
    ! SIDE EFFECTS
    !
    !  This routine modifies the values of the real or complex arrays
    !  invs_n, hjj, b_n, sumb_n, norma_sumb_n
    !
    ! RETURN VALUE
    !
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    subroutine init_invs()
      !
      !
      if (rmass_or_cmass_par%cmass_selected) then
         !
         call init_invs_c()
         call init_invs_r()
         !
      else if (rmass_or_cmass_par%rmass_selected) then
         !
      call init_invs_r()
      !
      else
         !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In subroutine init_invs case()'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'rmass_or_cmass_par has wrong value'
          call catch_exception(0)
          !
      end if
      !
    end subroutine init_invs
    !
    !
    subroutine init_invs_r()
      !
      integer :: i1,i2,i,j,k,pin_count
      real(ki), dimension(:,:), allocatable :: temp_mat_r,temp1_mat_r

      real(ki) :: error,tmp_error
      real(ki) :: plus_grand
      integer, dimension(6) :: pinch
      !
      plus_grand = maxval(array=abs(s_mat_r))
      b_ref = packb(set_ref)
      allocate(temp_mat_r(dim_s,dim_s),temp1_mat_r(dim_s,dim_s),stat=err)
      !
      if (err /= 0) then
         !
         tab_erreur_par(1)%a_imprimer = .true.
         tab_erreur_par(1)%chaine = 'In subroutine init_invs'
         tab_erreur_par(2)%a_imprimer = .true.
         tab_erreur_par(2)%chaine = 'cannot allocate memory for temp_mat and temp1_mat'
         call catch_exception(0)
         !
      end if
      !
      error=0
      do i = 0, 2**dim_s-2 ! iterate over all possible pinches
        pinch=-1
        pinch=unpackb(i*2,size(pinch))
        pin_count=countb(i)
        origine_inv_info_par = achar(dim_s+48)//'x'//achar(dim_s+48)//' matrix'
        !
        if (pin_count>0) then
                origine_inv_info_par = trim(origine_inv_info_par)//'pinch'
                origine_inv_info_par = trim(origine_inv_info_par)//' '//achar(pinch(j)+48)
                call put_to_zero(pinch(1),s_mat_r,temp_mat_r)
        else
                temp_mat_r=s_mat_r
        end if
        !
        do j = 2, pin_count
               origine_inv_info_par = trim(origine_inv_info_par)//' '//achar(pinch(j)+48)
               call put_to_zero(pinch(j),temp_mat_r,temp1_mat_r)
               temp_mat_r=temp1_mat_r
        end do
        invs_n_r(i+1,:,:)=0._ki
        call inverse(temp_mat_r,invs_n_r(i+1,:,:),tmp_error,pinch(1),pinch(2),pinch(3),pinch(4),pinch(5))
        b_n_r(i+1,:) = sum(invs_n_r(i+1,:,:),dim=1)
        sumb_n_r(i+1) = sum(b_n_r(i+1,:))
        norma_sumb_n_r(i+1) = sumb_n_r(i+1)*plus_grand
        if(pin_count<=2 .and. tmp_error>error) then ! for compatibility with old error variable
                error=tmp_error
        end if
      end do
      !
      if (dim_s==6) then
         do i=1,6
          do i1=1,6
            do i2=1,6
              hjj_r(i1,i2,i) = -2._ki*(  invs_n_r(1,i1,i2) &
                                 - invs_n_r(1,i,i1)*b_n_r(1,i2)/b_n_r(1,i) &
                                 - invs_n_r(1,i,i2)*b_n_r(1,i1)/b_n_r(1,i) &
                                 + invs_n_r(1,i,i)*b_n_r(1,i1)*b_n_r(1,i2) &
                                 /b_n_r(1,i)**2 )
            end do
          end do
        end do
      end if
      !
      deallocate(temp_mat_r,temp1_mat_r,stat=err)
      !
      if (err /= 0) then
         !
         tab_erreur_par(1)%a_imprimer = .true.
         tab_erreur_par(1)%chaine = 'In subroutine init_invs'
         tab_erreur_par(2)%a_imprimer = .true.
         tab_erreur_par(2)%chaine = 'cannot deallocate memory for temp_mat and temp1_mat'
         call catch_exception(0)
         !
      end if
    end subroutine init_invs_r
    !
    subroutine init_invs_c()
      !
      integer :: i1,i2,i,j,pin_count
      complex(ki), dimension(:,:), allocatable :: temp_mat_c,temp1_mat_c
      real(ki) :: error,tmp_error
      real(ki) :: plus_grand
      integer, dimension(6) :: pinch
      !
      plus_grand = maxval(array=abs(s_mat_c))
      b_ref = packb(set_ref)
      allocate(temp_mat_c(dim_s,dim_s),temp1_mat_c(dim_s,dim_s), stat=err)
      !
      if (err /= 0) then
         !
         tab_erreur_par(1)%a_imprimer = .true.
         tab_erreur_par(1)%chaine = 'In subroutine init_invs'
         tab_erreur_par(2)%a_imprimer = .true.
         tab_erreur_par(2)%chaine = 'cannot allocate memory for temp_mat and temp1_mat'
         call catch_exception(0)
         !
      end if
      !
      error=0
      do i = 0, 2**dim_s-2 ! iterate over all possible pinches
        pinch=-1
        pinch=unpackb(i*2,6)
        pin_count=countb(i)
        origine_inv_info_par = achar(dim_s+48)//'x'//achar(dim_s+48)//' matrix'
        if (pin_count>0) then
                origine_inv_info_par = trim(origine_inv_info_par)//'pinch'
                origine_inv_info_par = trim(origine_inv_info_par)//' '//achar(pinch(j)+48)
                call put_to_zero(pinch(1),s_mat_c,temp_mat_c)
        else
                temp_mat_c=s_mat_c
        end if
        !
        do j = 2, pin_count
               origine_inv_info_par = trim(origine_inv_info_par)//' '//achar(pinch(j)+48)
               call put_to_zero(pinch(j),temp_mat_c,temp1_mat_c)
               temp_mat_c=temp1_mat_c
        end do
        call inverse(temp_mat_c,invs_n_c(i+1,:,:),tmp_error,pinch(1),pinch(2),pinch(3),pinch(4),pinch(5))
        b_n_c(i+1,:) = sum(invs_n_c(i+1,:,:),dim=1)
        sumb_n_c(i+1) = sum(b_n_c(i+1,:))
        norma_sumb_n_c(i+1) = sumb_n_c(i+1)*plus_grand
        if(pin_count<=2 .and. tmp_error>error) then ! to be compatible with old error variable
                error=tmp_error
        end if
      end do
      !
      if (dim_s==6) then
         do i=1,6
          do i1=1,6
            do i2=1,6
              hjj_c(i1,i2,i) = -2._ki*(  invs_n_c(1,i1,i2) &
                                 - invs_n_c(1,i,i1)*b_n_c(1,i2)/b_n_c(1,i) &
                                 - invs_n_c(1,i,i2)*b_n_c(1,i1)/b_n_c(1,i) &
                                 + invs_n_c(1,i,i)*b_n_c(1,i1)*b_n_c(1,i2) &
                                 /b_n_c(1,i)**2 )
            end do
          end do
        end do
      end if
      !
      deallocate(temp_mat_c,temp1_mat_c,stat=err)
      !
      if (err /= 0) then
         !
         tab_erreur_par(1)%a_imprimer = .true.
         tab_erreur_par(1)%chaine = 'In subroutine init_invs'
         tab_erreur_par(2)%a_imprimer = .true.
         tab_erreur_par(2)%chaine = 'cannot deallocate memory for temp_mat and temp1_mat_c'
         call catch_exception(0)
         !
      end if
    end subroutine init_invs_c
    !
    !****if* src/kinematic/inversion/put_to_zero
    ! NAME
    !
    !  Subroutine put_to_zero
    !
    ! USAGE
    !
    !  call put_to_zero(i,mati,matf)
    !
    ! DESCRIPTION
    !
    ! This routine put to 0 the line and the column i of the square matrix mati
    ! It returns a square matrix matf of dim n x n (n being the
    ! dimension of mati). It is overloaded with a real and complex version.
    !
    ! INPUTS
    !
    !  * i -- an integer, the value of the line/column to be put to zero
    !  * mati -- an real/complex (type ki) array of rank 2
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * matf -- an real/copmplex (type ki) array of rank 2
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine put_to_zero_r(i,mati,matf)
      !
      integer, intent(in) :: i
      real(ki), dimension(:,:), intent(in) :: mati
      real(ki), dimension(size(mati,1),size(mati,2)), intent(out) :: matf
      !
      integer :: n
      !
      n = size(mati,1)          ! la matrice mati est carree
      matf = mati
      matf(i,:) = 0._ki
      matf(:,i) = 0._ki
      !
    end subroutine put_to_zero_r
    !
    subroutine put_to_zero_c(i,mati,matf)
      !
      integer, intent(in) :: i
      complex(ki), dimension(:,:), intent(in) :: mati
      complex(ki), dimension(size(mati,1),size(mati,2)), intent(out) :: matf
      !
      integer :: n
      !
      n = size(mati,1)          ! la matrice mati est carree
      matf = mati
      matf(i,:) = czero
      matf(:,i) = czero
      !
    end subroutine put_to_zero_c
    !
    !****f* src/kinematic/inversion/inv_s
    ! NAME
    !
    !  Function inv_s
    !
    ! USAGE
    !
    !  complex = inv_s(i,j,set)
    !
    ! DESCRIPTION
    !
    !  This function gives the generic inverse of the S matrix whatever
    !  its dimension (<=6)
    !
    ! INPUTS
    !
    !  * i -- an integer, line number
    !  * j -- an integer, row number
    !  * set -- an integer array of rank 1, the set of pinch propagators
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  Warning: Now a complex (type ki) is returned! [TK Sep10]
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    function inv_s(i,j,b_pin)
      !
      integer, intent (in) :: i,j
      integer, intent (in) :: b_pin
      complex(ki) :: inv_s
      !
      call check_pin(b_pin,"inv_s")
      if (ior(s_mat_p%b_cmplx,b_pin) .eq. b_pin) then
         !
         inv_s = cmplx(inv_s_r(i,j,b_pin),0._ki,ki)
         !
      else
         !
         inv_s = inv_s_c(i,j,b_pin)
         !
      end if
      !write (*, *) "INV_S(",i,j,b_pin,")=", inv_s
      !
    end function inv_s
    !
    function inv_s_r(i,j,b_pin)
      integer, intent (in) :: i,j
      integer, intent (in) :: b_pin
      real(ki) :: inv_s_r
      !
      inv_s_r =  invs_n_r((b_pin/2)+1,i,j)
    end function inv_s_r
    !
    function inv_s_c(i,j,b_pin)
      integer, intent (in) :: i,j
      integer, intent (in) :: b_pin
      complex(ki) :: inv_s_c
      !
      inv_s_c =  invs_n_c((b_pin/2)+1,i,j)
    end function inv_s_c
    !
    !****f* src/kinematic/inversion/hj
    ! NAME
    !
    !  Function hj
    !
    ! USAGE
    !
    !  complex = hj(i,j,set)
    !
    ! DESCRIPTION
    !
    ! This function gives the H matrix (pseudo-inverse of G) (dim=6)
    !
    ! INPUTS
    !
    !  * i -- an integer, line number
    !  * j -- an integer, row number
    !  * set -- an integer array of rank 1, the set of pinch propagators
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  Warning: Now a complex (type ki) is returned! [TK Sep10]
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    function hj(i,j,b_pin)
      !
      integer, intent (in) :: i,j
      integer, intent (in) :: b_pin
      complex(ki) :: hj
      !
      call check_pin(b_pin,'hj')
      if (ior(s_mat_p%b_cmplx,b_pin) .eq. b_pin) then
         !
         hj = cmplx(hj_r(i,j,b_pin),0._ki,ki)
         !
      else
         !
         hj = hj_c(i,j,b_pin)
         !
      end if
      !
    end function hj
    !
    function hj_r(i,j,b_pin)
      !
      integer, intent (in) :: i,j
      integer, intent (in) :: b_pin
      real(ki) :: hj_r
      !
      integer :: k
      integer, dimension(1) :: set
      integer :: dim_set
      !
      if (b_pin < 256) then
         dim_set = bit_count(b_pin)
         if (dim_set /= 0) then
           !
           !allocate(set(1:dim_set))
           k = bit_sets(b_pin*8)
           !
         else
           k = 0
         end if
      else
         dim_set = countb(b_pin)
         if (dim_set /= 0) then
           !
           !allocate(set(1:dim_set))
           set = unpackb(b_pin,1)
           k = set(1)
           !
         else
           k = 0
         end if
      end if
      !
      select case(dim_s)
      !
      case(6)                   ! case where we start with a 6-point amplitude
        !
        if (dim_set == 1) then
          if ( (i == k) .or. (j == k) ) then
            !
            hj_r = 0._ki
            !
          else
            !
            hj_r = hjj_r(i,j,k)
            !
          end if
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function hj, for 6-point'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the array set has not the right dimension: %d0'
          tab_erreur_par(2)%arg_int = dim_set
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
        end if
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function hj'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of the S matrix is not&
                          & correct %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      end select
    end function hj_r
    !
    function hj_c(i,j,b_pin)
      !
      integer, intent (in) :: i,j
      integer, intent (in) :: b_pin
      complex(ki) :: hj_c
      !
      integer :: k
      integer, dimension(1) :: set
      integer :: dim_set
      !
      if (b_pin < 256) then
         dim_set = bit_count(b_pin)
         if (dim_set /= 0) then
           !
           !allocate(set(1:dim_set))
           k = bit_sets(b_pin*8)
           !
         else
           k = 0
         end if
      else
         dim_set = countb(b_pin)
         if (dim_set /= 0) then
           !
           !allocate(set(1:dim_set))
           set = unpackb(b_pin,1)
           k = set(1)
           !
         else
           k = 0
         end if
      end if
      !
      select case(dim_s)
      !
      case(6)                   ! case where we start with a 6-point amplitude
        !
        if (dim_set == 1) then
          if ( (i == k) .or. (j == k) ) then
            !
            hj_c = czero
            !
          else
            !
            hj_c = hjj_c(i,j,k)
            !
          end if
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function hj, for 6-point'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the array set has not the right dimension: %d0'
          tab_erreur_par(2)%arg_int = dim_set
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
        end if
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function hj'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of the S matrix is not&
                          & correct %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      end select
    end function hj_c
    !
    !****f* src/kinematic/inversion/b
    ! NAME
    !
    !  Function b
    !
    ! USAGE
    !
    !  complex = b(i,set)
    !
    ! DESCRIPTION
    !
    ! This function gives the b coefficients whatever the S matrix dimension (<=6)
    !
    ! INPUTS
    !
    !  * i -- an integer, label of the b coefficients
    !  * set -- an integer array of rank 1, the set of pinch propagators
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  Warning: Now a complex (type ki) is returned! [TK Sep10]
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    function b(i,b_pin)
      !
      integer, intent (in) :: i
      integer, intent (in) :: b_pin
      complex(ki) :: b
      !
      call check_pin(b_pin,'b')
      if (ior(s_mat_p%b_cmplx,b_pin) .eq. b_pin) then
         !
         b = cmplx(b_r(i,b_pin),0._ki,ki)
         !
      else
         !
         b = b_c(i,b_pin)
         !
      end if
      ! write (*, *) "B(",i,b_pin,")=", b
      !
    end function b
    !
    function b_r(i,b_pin)
      integer, intent (in) :: i
      integer, intent (in) :: b_pin
      real(ki) :: b_r
      b_r= b_n_r((b_pin)/2+1,i)
    end function b_r
    function b_c(i,b_pin)
      integer, intent (in) :: i
      integer, intent (in) :: b_pin
      complex(ki) :: b_c
      b_c= b_n_c((b_pin/2)+1,i)
    end function b_c
    !
    !****f* src/kinematic/inversion/sumb
    ! NAME
    !
    !  Function sumb
    !
    ! USAGE
    !
    !  complex = sumb(set)
    !
    ! DESCRIPTION
    !
    ! This function gives the B coefficient whatever the S matrix dimension (<=6)
    !
    ! INPUTS
    !
    !  * set -- an integer array of rank 1, the set of pinch propagators
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  Warning: Now a complex (type ki) is returned! [TK Sep10]
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    function sumb(b_pin)
      !
      integer, intent (in) :: b_pin
      complex(ki) :: sumb
      !
      if (ior(s_mat_p%b_cmplx,b_pin) .eq. b_pin) then
         !
         sumb = cmplx(sumb_r(b_pin),0._ki,ki)
         !
      else
         !
         sumb = sumb_c(b_pin)
         !
      end if
      !
    end function sumb
    !
    function sumb_r(b_pin)
      !
      integer, intent (in) :: b_pin
      real(ki) :: sumb_r
      !
      sumb_r=sumb_n_r(b_pin/2+1)
    end function sumb_r
    function sumb_c(b_pin)
      !
      integer, intent (in) :: b_pin
      complex(ki) :: sumb_c
      sumb_c=sumb_n_c(b_pin/2+1)
    end function sumb_c
    !
    !****f* src/kinematic/inversion/norma_sumb
    ! NAME
    !
    !  Function norma_sumb
    !
    ! USAGE
    !
    !  complex = norma_sumb(set)
    !
    ! DESCRIPTION
    !
    ! This function gives the B coefficient whatever the S matrix dimension (<=6)
    ! divided by the greatest (in absolute value) element of the S matrix
    !
    ! INPUTS
    !
    !  * set -- an integer array of rank 1, the set of pinch propagators
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  Warning: Now a complex (type ki) is returned! [TK Sep10]
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    function norma_sumb(b_pin)
      !
      integer, intent(in) :: b_pin
      complex(ki) :: norma_sumb
      !
      call check_pin(b_pin,'norma_sumb')
      if (ior(s_mat_p%b_cmplx,b_pin) .eq. b_pin) then
         !
         norma_sumb = cmplx(norma_sumb_r(b_pin),0._ki,ki)
         !
      else
         !
         norma_sumb = norma_sumb_c(b_pin)
         !
      end if
      !
    end function norma_sumb
    !
    function norma_sumb_r(b_pin)
      !
      integer, intent (in) :: b_pin
      real(ki) :: norma_sumb_r
      !
      norma_sumb_r=norma_sumb_n_r(b_pin/2+1)
    end function
    function norma_sumb_c(b_pin)
      !
      integer, intent (in) :: b_pin
      complex(ki) :: norma_sumb_c
      !
      norma_sumb_c=norma_sumb_n_c(b_pin/2+1)
    end function
    !
    !
    !****f* src/kinematic/matrice_s/exitgolem95
    ! NAME
    !
    !  Subroutine exitgolem95
    !
    ! USAGE
    !
    !  call exitgolem95()
    !
    ! DESCRIPTION
    !
    !  This subroutine should be called at the end of the form factor calculation.
    !  It frees all memory previously allocated, it clears the cache and nullifies pointers.
    !
    ! INPUTS
    !
    ! SIDE EFFECTS
    !
    ! RETURN VALUE
    !
    !  No return value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine exitgolem95()
      !
      rmass_or_cmass_par = cmass
      !
      nullify(s_mat)
      !
      call deallocation_s()
      !
      call clear_cache()
      !
    end subroutine exitgolem95
    !
    subroutine check_pin(b_pin,func)
      integer, intent(in) :: b_pin
      character(*), intent(in) :: func
      !
      if (dim_s<= 0 .or. dim_s>=7) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine '//trim(func(:))
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'dimension %d0 not supported.'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
      end if
      !
      if (b_pin>=2**(dim_s+1)-1 .or. iand(b_pin,1)==1) then ! do not allow complete pinch or "0" pinched
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine '//func(:)
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'parameter b_pin=%d0 invalid.'
        tab_erreur_par(2)%arg_int = b_pin
        call catch_exception(0)
      end if
    end subroutine check_pin
end module matrice_s
