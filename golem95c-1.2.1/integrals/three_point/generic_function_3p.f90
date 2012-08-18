!****h* src/integrals/three_point/generic_function_3p
! NAME
!
!  Module generic_function_3p
!
! USAGE
!
!  use generic_function_3p
!
! DESCRIPTION
!
!  This module contains the generic routines to compute the
!  three point functions in n and n+2 dimensions
!
! OUTPUT
!
!  It exports two public routines:
!  * f3p(_sc)     -- a function to compute the three point function in n dimensions
!  * f3p_np2(_sc) -- a function to compute the three point function in n+2 dimensions
!    Calling the functions with _sc returns a real array. These calls will not be cached.
!   
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * matrice_s (src/kinematic/matrice_s.f90) ( only : dim_s, b_ref )
!  * s_matrix_type (src/module/s_matrix_type.f90)
!  * array (src/module/array.f90)
!  * tri_croissant (src/module/tri.f90)
!  * constante (src/module/constante.f90)
!  * function_3p1m (src/integrals/three_point/function_3p1m.f90)
!  * function_3p2m (src/integrals/three_point/function_3p2m.f90)
!  * function_3p3m (src/integrals/three_point/function_3p3m.f90)
!  * cache (src/module/cache.f90)
!  * equal (src/module/equal.f90)
!
!*****
module generic_function_3p
  !
  use precision_golem
  use matrice_s, only : dim_s,b_ref
  use s_matrix_type
  use array
  use tri_croissant
  use constante
  use function_3p0m_1mi
  use function_3p1m
  use function_3p1m_1mi
  use function_3p1m_2mi
  use function_3p2m
  use function_3p2m_1mi
  use function_3p3m
  use parametre
  use function_3p_finite
  use cache
  use equal
  use sortie_erreur, only : tab_erreur_par,catch_exception
  implicit none
  !
  private
  !
  integer, dimension(:), allocatable :: set
  integer, dimension(3) :: set_tot
  !
  interface f3p_sc
     module procedure f3p_sc_r, f3p_sc_c
     module procedure f3p_sc_p
  end interface
  !
  interface f3p_np2_sc
     module procedure f3p_np2_sc_r, f3p_np2_sc_c
     module procedure f3p_np2_sc_p
  end interface
  !
  public :: f3p, f3p_np2
!  public :: f3p_ra, f3p_np2_ra !!! return real arrays. not needed in current implementation
  public :: f3p_sc, f3p_np2_sc !!! (non-cached) return real arrays. 
  !can be called with real/complex arrays in addition.
  !
contains
  !
  !****f* src/integrals/three_point/generic_function_3p/f3p
  ! NAME
  !
  !  Function f3p
  !
  ! USAGE
  !
  !  complex_dim3 = f3p(s_mat_p, b_pro, parf1, parf2, parf3)
  !
  ! DESCRIPTION
  !
  !  This function computes the generic three point function in n dimensions, 
  !  with or without Feynman parameters in the numerator.
  !
  ! INPUTS
  !
  !  * s_mat_p -- a type s_matrix_poly object, the S matrix
  !  * b_pro -- an integer whose digits represents the set of the three unpinched
  !    propagators
  !  * parf1 -- an integer (optional), the label of the one Feynman parameter
  !  * parf2 -- an integer (optional), the label of the second Feynman parameter
  !  * parf3 -- an integer (optional), the label of the third Feynman parameter
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki) array of rank 1 and shape 3
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  !
  function f3p(s_mat_p, b_pro,parf1,parf2,parf3)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent (in) :: b_pro
    integer, intent (in), optional :: parf1,parf2,parf3
    complex(ki), dimension(3) :: f3p
    real(ki), dimension(6) :: f3p_real
    !
    f3p_real = f3p_ra(s_mat_p,b_pro,parf1=parf1,parf2=parf2,parf3=parf3)
    f3p(1) = f3p_real(1) + i_ * f3p_real(2)
    f3p(2) = f3p_real(3) + i_ * f3p_real(4)
    f3p(3) = f3p_real(5) + i_ * f3p_real(6)
    !
  end function f3p
  !
  function f3p_ra(s_mat_p,b_pro,parf1,parf2,parf3)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent (in) :: b_pro
    integer, intent (in), optional :: parf1,parf2,parf3
    real(ki),dimension(6) :: f3p_ra
    !
    integer :: par1,par2,par3
    integer :: par_cache1,par_cache2,par_cache3
    integer, dimension(3) :: z_param_ini,z_param_out
    integer :: taille
    integer :: b_pin
    integer, dimension(3) :: s
    !
    par1 = 0
    par2 = 0
    par3 = 0
    !
    if (present(parf1)) par1 = parf1 
    if (present(parf2)) par2 = parf2
    if (present(parf3)) par3 = parf3
    !
    z_param_ini = (/ par1,par2,par3 /)
    !
    where (z_param_ini /= 0)
       !
       z_param_ini = locateb(z_param_ini,b_pro)
       !
    elsewhere
       !
       z_param_ini = 0
       !
    end where
    !
    if ( minval(z_param_ini) == -1 ) then
       !
       f3p_ra = 0.0_ki
       !
    else
       !
       s = unpackb(b_pro,countb(b_pro))
       taille = dim_s - size(s)
       !
       select case(taille)
          !
       case(0)
          !
          set_tot = 0
          !
       case(1)
          !
          allocate(set(1:taille))
          b_pin = pminus(b_ref,b_pro)
          set = unpackb(b_pin,countb(b_pin))
          set_tot(1:2) = 0
          set_tot(3) = set(1)
          !
       case(2)
          !
          allocate(set(1:taille))
          b_pin = pminus(b_ref,b_pro)
          set = unpackb(b_pin,countb(b_pin))
          set_tot(1) = 0
          set_tot(2:3) = set
          !
       case(3)
          !
          allocate(set(1:taille))
          b_pin = pminus(b_ref,b_pro)
          set = unpackb(b_pin,countb(b_pin))
          set_tot = set
          !
       case default
          !
          set_tot = 0
          taille = 0
          !
       end select
       !
       call tri_int3(z_param_ini,z_param_out)
       par_cache1 = z_param_out(1)
       par_cache2 = z_param_out(2)
       par_cache3 = z_param_out(3)
       !
       cache : if ( computed_f3p(set_tot(1),set_tot(2),set_tot(3),&
            &par_cache1,par_cache2,par_cache3) ) then
          !
          f3p_ra = results_f3p(set_tot(1),set_tot(2),set_tot(3),&
               &par_cache1,par_cache2,par_cache3,:)
          !
       else cache
          !
          f3p_ra = f3p_sc(s_mat_p,s,par_cache1,par_cache2,par_cache3)
          !
          computed_f3p(set_tot(1),set_tot(2),set_tot(3),&
               &par_cache1,par_cache2,par_cache3) = .true.
          results_f3p(set_tot(1),set_tot(2),set_tot(3),&
               &par_cache1,par_cache2,par_cache3,:) = f3p_ra
          !
       end if cache
       !
       if (taille /= 0) deallocate(set)
       !
    end if
    !
  end function f3p_ra
  !
  !****f* src/integrals/three_point/generic_function_3p/f3p_sc
  ! NAME
  !
  !  Function f3p_sc
  !
  ! USAGE
  !
  !  real_dim6 = f3p_sc(s_mat,s,parf1,parf2,parf3)
  !
  ! DESCRIPTION
  !
  !  This function computes the generic three point function in n dimensions, 
  !  with or without Feynman parameters in the numerator without using a cache
  !
  ! INPUTS
  !
  !  * s_mat -- a real/complex (type ki)/s_matrix_poly array of rank 2, the S matrix
  !  * s -- an integer array of rank 1 and shape 3, the set of the three unpinched
  !    propagators
  !  * parf1 -- an integer (optional), the label of the one Feynman parameter
  !  * parf2 -- an integer (optional), the label of the second Feynman parameter
  !  * parf3 -- an integer (optional), the label of the third Feynman parameter
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  It returns a real (type ki) array of rank 1 and shape 6
  !
  ! EXAMPLE
  !
  !
  !
  !*****

  function f3p_sc_p(s_mat_p,s,parf1,parf2,parf3)
    type(s_matrix_poly) :: s_mat_p
    integer, intent (in), dimension(3) :: s
    integer, intent(in), optional :: parf1, parf2, parf3
    real(ki), dimension(6) :: f3p_sc_p
    !
    if (iand(s_mat_p%b_cmplx, packb(s)) .eq. 0 ) then
       !
       f3p_sc_p = f3p_sc_r(s_mat_p%pt_real, s, parf1=parf1,parf2=parf2,parf3=parf3)
       !
    else
       !
       f3p_sc_p = f3p_sc_c(s_mat_p%pt_cmplx, s, parf1=parf1,parf2=parf2,parf3=parf3)
       !
    end if
    !
  end function f3p_sc_p
  !
  function f3p_sc_r(s_mat_r,s,parf1,parf2,parf3)
    !
    real(ki), intent (in), dimension(:,:) :: s_mat_r
    integer, intent (in), dimension(3) :: s
    integer, intent (in), optional :: parf1,parf2,parf3
    real(ki),dimension(6) :: f3p_sc_r
    !
    integer :: par1,par2,par3
    integer, dimension(3) :: z_param_ini,z_param_out
    real(ki) :: arg1,arg2,arg3,s1,s2,s3
    real(ki) :: mass1,mass2,mass3
    integer :: m1,m2,m3
    logical, dimension(3) :: argz, mz,sz
    !
    par1 = 0
    par2 = 0
    par3 = 0
    !
    if (present(parf1)) par1 = parf1 
    if (present(parf2)) par2 = parf2
    if (present(parf3)) par3 = parf3
    !
    if ( (par1 == -1) .or. (par2 == -1) .or. (par3 == -1) ) then
       !
       f3p_sc_r(:) = 0._ki
       !
    else
       ! symetrie: la place de z1,z2,z3 n'a pas d'importance, on les met 
       ! dans l'ordre croissant
       z_param_ini(1) = par1
       z_param_ini(2) = par2
       z_param_ini(3) = par3
       !
       m1 = s(1)
       m2 = s(2)
       m3 = s(3)
       !
       arg1 = s_mat_r(m1,m2)
       arg2 = s_mat_r(m2,m3)
       arg3 = s_mat_r(m1,m3)
       !
       argz(1) = equal_real(arg1,zero)
       argz(2) = equal_real(arg2,zero)
       argz(3) = equal_real(arg3,zero)
       !
       ! internal masses
       mass1 = -s_mat_r(m1,m1)/2._ki
       mass2 = -s_mat_r(m2,m2)/2._ki
       mass3 = -s_mat_r(m3,m3)/2._ki
       !
       mz(1) = equal_real(mass1,zero)
       mz(2) = equal_real(mass2,zero)
       mz(3) = equal_real(mass3,zero)
       !
       ! external p_i^2
       s1=arg3+mass1+mass3
       s2=arg1+mass1+mass2
       s3=arg2+mass2+mass3
       !
       sz(1) = equal_real(s1,zero)
       sz(2) = equal_real(s2,zero)
       sz(3) = equal_real(s3,zero)
       !
       call cut_s(s1,mass1,mass3)
       call cut_s(s2,mass1,mass2)
       call cut_s(s3,mass2,mass3)
       !
       ! initialize all components
       f3p_sc_r(:) = 0._ki
       !
       ! the integrals are classified by the off-shellness of the external legs
       !
       !~ case with one light-like, two massive on-shell legs: QL tri5
       !
       if ( ( argz(1) ) .and. ( argz(2) ) .and. ( argz(3) ) ) then
          !
          ! comment 11.08.10: single out call with s_matrix being zero:
          !
          if  (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) )  ) then
             tab_erreur_par(1)%a_imprimer = .true.
             tab_erreur_par(1)%chaine = 'In function generic_function_3p.f90:'
             tab_erreur_par(2)%a_imprimer = .true.
             tab_erreur_par(2)%chaine = 'function call with all arguments zero!'
             call catch_exception(1)
             !    
             !
             ! case with one internal mass, two on-shell massive legs
             !~ QL tri5, two on-shell massive legs
             !
          else if (  ( .not.( mz(3) ) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             call tri_int3(z_param_ini,z_param_out)
             !
             f3p_sc_r = f3p0m_1mi(mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( .not.( mz(1) ) ) .and. ( mz(2) ) .and. ( mz(3) ) ) then
             !
             call shift_param(z_param_ini,2,3,z_param_out)
             !
             f3p_sc_r = f3p0m_1mi(mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( .not.( mz(2) ) ) .and. ( mz(1) ) .and. ( mz(3) ) ) then
             !
             call shift_param(z_param_ini,1,3,z_param_out)
             !
             f3p_sc_r = f3p0m_1mi(mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else
             ! comment 11.08.10: only one internal mass possible kinematically!
             tab_erreur_par(1)%a_imprimer = .true.
             tab_erreur_par(1)%chaine = 'In function generic_function_3p.f90:'
             tab_erreur_par(2)%a_imprimer = .true.
             tab_erreur_par(2)%chaine = 'function call with one lightlike, two massive external legs &
                  & and more than one internal mass!'
             tab_erreur_par(3)%a_imprimer = .true.
             tab_erreur_par(3)%chaine = 'This should not be allowed kinematically!'
             call catch_exception(1)
             !
             call tri_int3(z_param_ini,z_param_out)
             ! finite triangle
             f3p_sc_r(3:6)=f3p_finite("ndi",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if ! end scan internal mass
          !
          !
          ! cases with one off-shell leg (one argi nonzero): 
          !~ QL tri1,tri4,tri6
          !
       else if ( ( argz(1) ) .and. ( argz(2) ) ) then
          !
          call tri_int3(z_param_ini,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             !~ QL tri1 
             f3p_sc_r = f3p1m(arg3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             !~ case with one internal mass: QL tri4
             !
          else if (  ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_sc_r = f3p1m_1mi(arg3,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(2) ) .and. ( mz(3) ) ) then
             !
             call exchange_param(z_param_ini,(/1,3/),3,z_param_out)
             !
             f3p_sc_r = f3p1m_1mi(arg3,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(1) ) .and. ( mz(3) ) ) then
             ! no, this one is finite, corrected June 3, 2010
             f3p_sc_r(3:6)=f3p_finite("ndi",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             !~ case with two internal masses: QL tri6
             !
          else if ( mz(2) ) then
             !
             f3p_sc_r = f3p1m_2mi(arg3,mass1,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else 
             ! finite triangle
             f3p_sc_r(3:6)=f3p_finite("ndi",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          ! cyclic permutation of arguments
          !
       else if ( (argz(2)) .and. (argz(3)) ) then
          !
          call shift_param(z_param_ini,2,3,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_sc_r = f3p1m(arg1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass
             !
          else if (  ( mz(2) ) .and. ( mz(3) ) ) then
             !
             f3p_sc_r = f3p1m_1mi(arg1,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(3) ) .and. ( mz(1) ) ) then
             !
             ! labels (/1,2/) -> (/2,3/) corrected 19.7.
             call exchange_param(z_param_ini,(/2,3/),3,z_param_out)
             !
             f3p_sc_r = f3p1m_1mi(arg1,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(2) ) .and. ( mz(1) ) ) then
             !
             ! finite
             f3p_sc_r(3:6)=f3p_finite("ndi",s2,s3,s1,mass2,mass3,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with two internal masses
             !
          else if ( mz(3) ) then
             !
             ! mass labels corrected 19.7.10
             f3p_sc_r = f3p1m_2mi(arg1,mass2,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else 
             ! finite triangle
             f3p_sc_r(3:6)=f3p_finite("ndi",s2,s3,s1,mass2,mass3,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
       else if ( (argz(1)) .and. (argz(3)) ) then
          !
          call shift_param(z_param_ini,1,3,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_sc_r = f3p1m(arg2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass
             !
          else if (  ( mz(1) ) .and. ( mz(3) ) ) then
             !
             f3p_sc_r = f3p1m_1mi(arg2,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(2) ) .and. ( mz(1) ) ) then
             !
             call exchange_param(z_param_ini,(/1,2/),3,z_param_out)
             !
             f3p_sc_r = f3p1m_1mi(arg2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(2) ) .and. ( mz(3) ) ) then
             !
             f3p_sc_r(3:6) = f3p_finite("ndi",s3,s1,s2,mass3,mass1,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with two internal masses
             !
          else if ( mz(1) ) then
             !
             ! mass labels corrected 19.7.10
             f3p_sc_r = f3p1m_2mi(arg2,mass3,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else 
             ! finite triangle
             f3p_sc_r(3:6)=f3p_finite("ndi",s3,s1,s2,mass3,mass1,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          !
          ! cases with two off-shell legs:  QL tri2,tri3
          !
       else if ( (argz(1)) .and. (.not.(argz(2))) .and. (.not.(argz(3))) ) then
          !
          call tri_int3(z_param_ini,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             !~ QL tri2
             f3p_sc_r = f3p2m(arg2,arg3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass: QL tri3
             !
          else if (  ( mz(1) ) .and. ( mz(2) ) .and.  &
                                ! corrected 3.6.10
               &         (.not. equal_real(mass3,s1) ) .and. (.not. equal_real(mass3,s3) )  ) then
             !
             f3p_sc_r = f3p2m_1mi(arg2,arg3,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else 
             ! finite triangle
             f3p_sc_r(3:6)=f3p_finite("ndi",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          ! permute arguments
       else if ( (argz(2)) .and. (.not.(argz(1))) .and. (.not.(argz(3))) ) then
          !
          call shift_param(z_param_ini,2,3,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_sc_r = f3p2m(arg3,arg1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass: QL tri3 
             !
          else if (  ( mz(2) ) .and. ( mz(3) ) .and. &
                                ! corrected 3.6.10                          
               &         (.not. equal_real(mass1,s2) ) .and. (.not. equal_real(mass1,s1) ) ) then
             !
             f3p_sc_r = f3p2m_1mi(arg3,arg1,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else 
             ! finite triangle
             f3p_sc_r(3:6)=f3p_finite("ndi",s2,s3,s1,mass2,mass3,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !	    
          end if
          !
       else if ( (argz(3)) .and. (.not.(argz(1))) .and. (.not.(argz(2))) ) then
          !
          call shift_param(z_param_ini,1,3,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_sc_r = f3p2m(arg1,arg2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass: QL tri3
             ! 
          else if (  ( mz(1) ) .and. ( mz(3) ) .and.  &
               &         (.not. equal_real(mass2,s3) ) .and. (.not. equal_real(mass2,s2) ) ) then
             !
             f3p_sc_r = f3p2m_1mi(arg1,arg2,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else 
             ! finite triangle
             f3p_sc_r(3:6)=f3p_finite("ndi",s3,s1,s2,mass3,mass1,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
       else  if ( .not.(argz(3)) .and. (.not.(argz(1))) .and. (.not.(argz(2))) ) then
          ! finite
          !
          if (  ( mz(1) ) .and. ( mz(2) )  .and. ( mz(3) ) ) then
             !
             if ( (abs(arg3) >= abs(arg1)) .and. (abs(arg3) >= abs(arg2)) ) then
                !
                call tri_int3(z_param_ini,z_param_out)
                f3p_sc_r(3:6) = f3p3m("ndi",arg3,arg1,arg2, &
                     z_param_out(1),z_param_out(2),z_param_out(3))
                !
             else if ( (abs(arg1) >= abs(arg3)) .and. (abs(arg1) >= abs(arg2)) ) then
                !
                call shift_param(z_param_ini,2,3,z_param_out)
                f3p_sc_r(3:6) = f3p3m("ndi",arg1,arg2,arg3, &
                     z_param_out(1),z_param_out(2),z_param_out(3))
                !
             else if ( (abs(arg2) >= abs(arg3)) .and. (abs(arg2) >= abs(arg1)) ) then
                !
                call shift_param(z_param_ini,1,3,z_param_out)
                f3p_sc_r(3:6) = f3p3m("ndi",arg2,arg3,arg1, &
                     z_param_out(1),z_param_out(2),z_param_out(3))
                !
             end if ! end if scan abs(argi)
             !
          else  ! internal masses present
             !
             ! finite triangle with internal masses
             call tri_int3(z_param_ini,z_param_out)
             f3p_sc_r(3:6)=f3p_finite("ndi",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
            !
          end if ! end if 3-off-shell legs triangle
          !
       else  ! other values of arg should not occur
          !
          ! finite triangle with internal masses
          f3p_sc_r(3:6)=f3p_finite("ndi",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
          ! 
       end if  ! end if arg1, arg2,arg3 nonzero
       !
    end if ! end if par1==-1 ...
    !
  end function f3p_sc_r
  !
  function f3p_sc_c(s_mat_c,s,parf1,parf2,parf3)
    !
    complex(ki), intent (in), dimension(:,:) :: s_mat_c
    integer, intent (in), dimension(3) :: s
    integer, intent (in), optional :: parf1,parf2,parf3
    real(ki),dimension(6) :: f3p_sc_c
    !
    integer :: par1,par2,par3
    integer, dimension(3) :: z_param_ini,z_param_out
    complex(ki) :: arg1,arg2,arg3
    real(ki) ::s1,s2,s3
    complex(ki) :: mass1,mass2,mass3
    integer :: m1,m2,m3
    logical, dimension(3) :: argz, mz
    logical :: finite = .true.
    !
    !
    par1 = 0
    par2 = 0
    par3 = 0
    !
    if (present(parf1)) par1 = parf1 
    if (present(parf2)) par2 = parf2
    if (present(parf3)) par3 = parf3
    !
    !
    if ( (par1 == -1) .or. (par2 == -1) .or. (par3 == -1) ) then
       !
       f3p_sc_c(:) = 0._ki
       !
    else
       ! symetrie: la place de z1,z2,z3 n'a pas d'importance, on les met 
       ! dans l'ordre croissant
       z_param_ini(1) = par1
       z_param_ini(2) = par2
       z_param_ini(3) = par3
       !
       m1 = s(1)
       m2 = s(2)
       m3 = s(3)
       !
       arg1 = s_mat_c(m1,m2)
       arg2 = s_mat_c(m2,m3)
       arg3 = s_mat_c(m1,m3)
       !
       argz(1) = equal_real(real(arg1,ki),zero) .and. equal_real(aimag(arg1),zero) 
       argz(2) = equal_real(real(arg2,ki),zero) .and. equal_real(aimag(arg2),zero) 
       argz(3) = equal_real(real(arg3,ki),zero) .and. equal_real(aimag(arg3),zero) 
       !
       ! internal masses
       mass1 = -s_mat_c(m1,m1)/2._ki
       mass2 = -s_mat_c(m2,m2)/2._ki
       mass3 = -s_mat_c(m3,m3)/2._ki
       !
       mz(1) = equal_real(real(mass1,ki),zero) .and. equal_real(aimag(mass1),zero)
       mz(2) = equal_real(real(mass2,ki),zero) .and. equal_real(aimag(mass2),zero)
       mz(3) = equal_real(real(mass3,ki),zero) .and. equal_real(aimag(mass3),zero)
       !
       ! external p_i^2
       s1 = real(arg3+mass1+mass3,ki)
       s2 = real(arg1+mass1+mass2,ki)
       s3 = real(arg2+mass2+mass3,ki)
       !
       call cut_s(s1,mass1,mass3)
       call cut_s(s2,mass1,mass2)
       call cut_s(s3,mass2,mass3)
       !
       ! initialize all components
       !
       f3p_sc_c(:) = 0._ki
       !
       ! the integrals are classified by the off-shellness of the external legs
       !
       ! In complex case, there is only one divergent triangle, QL tri3
       !
       finite = .true.
       !
       if ( (argz(1)) .and. (.not.(argz(2))) .and. (.not.(argz(3))) ) then
          !
          call tri_int3(z_param_ini,z_param_out)
          !
          if ( ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_sc_c = f3p2m_1mi(arg2,arg3,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             finite = .false.
             !
          end if
          !
          ! permute arguments
       else if ( (argz(2)) .and. (.not.(argz(1))) .and. (.not.(argz(3))) ) then
          !
          call shift_param(z_param_ini,2,3,z_param_out)
          !
          if ( ( mz(2) ) .and. ( mz(3) ) ) then
             !
             f3p_sc_c = f3p2m_1mi(arg3,arg1,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             finite = .false.
             !    
          end if
          !
       else if ( (argz(3)) .and. (.not.(argz(1))) .and. (.not.(argz(2))) ) then
          !
          call shift_param(z_param_ini,1,3,z_param_out)
          !
          if ( ( mz(1) ) .and. ( mz(3) )  ) then
	  ! .and. (.not. equal_real(s2,zero) )
             !
             f3p_sc_c = f3p2m_1mi(arg1,arg2,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             finite = .false.
             !
          end if
          !
       end if !argz
       !
       if (finite) then !finite triangle
          !
          call tri_int3(z_param_ini,z_param_out)
          !
          f3p_sc_c(3:6)=f3p_finite("ndi",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
          ! 
       end if  ! end if call finite
       !
    end if ! end if par1==-1 ...
    !
  end function f3p_sc_c
  !
  !****f* src/integrals/three_point/generic_function_3p/f3p_np2
  ! NAME
  !
  !  Function f3p_np2
  !
  ! USAGE
  !
  !  complex_dim2 = f3p_np2_ca(s_mat,b_pro,parf1,parf2,parf3)
  ! 
  ! DESCRIPTION
  !
  !  This function computes the generic three point function in n+2 dimensions, 
  !  with or without Feynman parameters in the numerator
  !
  ! INPUTS
  !
  !  * s_mat -- a s_matrix_poly type object, the S matrix
  !  * b_pro -- an integer whose digits represents the set of the three unpinched
  !    propagators
  !  * parf1 -- an integer (optional), the label of the one Feynman parameter
  !  * parf2 -- an integer (optional), the label of the second Feynman parameter
  !  * parf3 -- an integer (optional), the label of the third Feynman parameter
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki) array of rank 1 and shape 2
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  !
  function f3p_np2(s_mat_p, b_pro,parf1,parf2,parf3)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent (in) :: b_pro
    integer, intent (in), optional :: parf1,parf2,parf3
    complex(ki), dimension(2) :: f3p_np2
    real(ki), dimension(4) :: f3p_np2_real
    !
    f3p_np2_real = f3p_np2_ra(s_mat_p,b_pro,parf1=parf1,parf2=parf2,parf3=parf3)
    !
    f3p_np2(1) = f3p_np2_real(1) + i_ * f3p_np2_real(2)
    f3p_np2(2) = f3p_np2_real(3) + i_ * f3p_np2_real(4)
    !
  end function f3p_np2
  !
  function f3p_np2_ra(s_mat_p,b_pro,parf1,parf2,parf3)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent (in) :: b_pro
    integer, intent (in), optional :: parf1,parf2,parf3
    real(ki),dimension(4) :: f3p_np2_ra
    !
    integer :: par1,par2,par3
    integer :: par_cache1,par_cache2,par_cache3
    integer, dimension(3) :: z_param_ini,z_param_out
    integer :: taille
    integer :: b_pin
    integer, dimension(3) :: s
    !
    par1 = 0
    par2 = 0
    par3 = 0
    !
    if (present(parf1)) par1 = parf1 
    if (present(parf2)) par2 = parf2
    if (present(parf3)) par3 = parf3
    !
    !
    z_param_ini = (/ par1,par2,par3 /)
    !
    where (z_param_ini /= 0)
       !
       z_param_ini = locateb(z_param_ini,b_pro)
       !
    elsewhere
       !
       z_param_ini = 0
       !
    end where
    !
    if ( minval(z_param_ini) == -1 ) then
       !
       f3p_np2_ra = 0._ki
       !
    else
       !
       s = unpackb(b_pro,countb(b_pro))
       taille = dim_s - size(s)
       !
       select case(taille)
          !
       case(0)
          !
          set_tot = 0
          !
       case(1)
          !
          allocate(set(1:taille))
          b_pin = pminus(b_ref,b_pro)
          set = unpackb(b_pin,countb(b_pin))
          set_tot(1:2) = 0
          set_tot(3) = set(1)
          !
       case(2)
          !
          allocate(set(1:taille))
          b_pin = pminus(b_ref,b_pro)
          set = unpackb(b_pin,countb(b_pin))
          set_tot(1) = 0
          set_tot(2:3) = set
          !
       case(3)
          !
          allocate(set(1:taille))
          b_pin = pminus(b_ref,b_pro)
          set = unpackb(b_pin,countb(b_pin))
          set_tot = set
          !
       case default
          !
          set_tot = 0
          taille = 0
          !
          !
       end select
       !
       call tri_int3(z_param_ini,z_param_out)
       par_cache1 = z_param_out(1)
       par_cache2 = z_param_out(2)
       par_cache3 = z_param_out(3)
       !
       cache : if ( computed_f3p_np2(set_tot(1),set_tot(2),set_tot(3),&
            &par_cache1,par_cache2,par_cache3) ) then
          !
          f3p_np2_ra = results_f3p_np2(set_tot(1),set_tot(2),set_tot(3),&
               &par_cache1,par_cache2,par_cache3,:)
          !
       else cache
          !
          f3p_np2_ra = f3p_np2_sc(s_mat_p,s,par_cache1,par_cache2,par_cache3)
          !
          computed_f3p_np2(set_tot(1),set_tot(2),set_tot(3),&
               &par_cache1,par_cache2,par_cache3) = .true.
          results_f3p_np2(set_tot(1),set_tot(2),set_tot(3),&
               &par_cache1,par_cache2,par_cache3,:) = f3p_np2_ra
          !
       end if cache
       !
       if (taille /= 0) deallocate(set)
       !
    end if
    !
  end function f3p_np2_ra
  !
  !****f* src/integrals/three_point/generic_function_3p/f3p_np2_sc
  ! NAME
  !
  !  Function f3p_np2_sc
  !
  ! USAGE
  !
  !  real_dim4 = f3p_np2_sc(s_mat,s,parf1,parf2,parf3)
  !
  ! DESCRIPTION
  !
  !  This function computes the generic three point function in n+2 dimensions, 
  !  with or without Feynman parameters in the numerator
  !
  ! INPUTS
  !
  !  * s_mat -- a real/complex (type ki)/s_matrix_poly array of rank 2, the S matrix
  !  * s -- an integer array of rank 1 and shape 3, the set of the three unpinched
  !    propagators
  !  * parf1 -- an integer (optional), the label of the one Feynman parameter
  !  * parf2 -- an integer (optional), the label of the second Feynman parameter
  !  * parf3 -- an integer (optional), the label of the third Feynman parameter
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  It returns a real (type ki) array of rank 1 and shape 4
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  function f3p_np2_sc_p(s_mat_p,s,parf1,parf2,parf3)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent (in), dimension(3) :: s
    integer, intent (in), optional :: parf1,parf2,parf3
    real(ki),dimension(4) :: f3p_np2_sc_p
    !
    if (iand(s_mat_p%b_cmplx, packb(s)) .eq. 0 ) then
       !
       f3p_np2_sc_p = f3p_np2_sc_r(s_mat_p%pt_real, s, parf1=parf1,parf2=parf2,parf3=parf3)
       !
    else
       !
       f3p_np2_sc_p = f3p_np2_sc_c(s_mat_p%pt_cmplx, s, parf1=parf1,parf2=parf2,parf3=parf3)
       !
    end if
    !
  end function f3p_np2_sc_p
  !
  function f3p_np2_sc_r(s_mat_r,s,parf1,parf2,parf3)
    !
    real(ki), intent (in), dimension(:,:) :: s_mat_r
    integer, intent (in), dimension(3) :: s
    integer, intent (in), optional :: parf1,parf2,parf3
    real(ki),dimension(4) :: f3p_np2_sc_r
    !
    integer :: par1,par2,par3
    integer, dimension(3) :: z_param_ini,z_param_out
    real(ki) :: arg1,arg2,arg3
    real(ki) :: mass1,mass2,mass3,s1,s2,s3
    integer :: m1,m2,m3
    logical, dimension(3) :: argz,mz
    !
!    calls_f3p_np2_sc = calls_f3p_np2_sc + 1
    !
    par1 = 0
    par2 = 0
    par3 = 0
    !
    if (present(parf1)) par1 = parf1 
    if (present(parf2)) par2 = parf2
    if (present(parf3)) par3 = parf3
    !
    !
    if ( (par1 == -1) .or. (par2 == -1) .or. (par3 == -1) ) then
       !
       f3p_np2_sc_r(:) = 0._ki
       !
    else
       ! symetrie: la place de z1,z2,z3 n'a pas d'importance, on les met 
       ! dans l'ordre croissant
       z_param_ini(1) = par1
       z_param_ini(2) = par2
       z_param_ini(3) = par3
       !
       m1 = s(1)
       m2 = s(2)
       m3 = s(3)
       !
       !
       arg1 = s_mat_r(m1,m2)
       arg2 = s_mat_r(m2,m3)
       arg3 = s_mat_r(m1,m3)
       !
       argz(1) = equal_real(arg1,zero)
       argz(2) = equal_real(arg2,zero)
       argz(3) = equal_real(arg3,zero)
       !
       ! internal masses
       mass1 = -s_mat_r(m1,m1)/2._ki
       mass2 = -s_mat_r(m2,m2)/2._ki
       mass3 = -s_mat_r(m3,m3)/2._ki
       !
       mz(1) = equal_real(mass1,zero)
       mz(2) = equal_real(mass2,zero)
       mz(3) = equal_real(mass3,zero)
       !
       !
       s1=arg3+mass1+mass3
       s2=arg1+mass1+mass2
       s3=arg2+mass2+mass3
       !
       call cut_s(s1,mass1,mass3)
       call cut_s(s2,mass1,mass2)
       call cut_s(s3,mass2,mass3)
       !
       ! the integrals are classified by the off-shellness of the external legs
       !
       ! case with all external legs on shell
       !
       ! initialize all components
       f3p_np2_sc_r(:) = 0._ki
       !
       !
       if ( ( argz(1) ) .and. ( argz(2) ) .and. ( argz(3) ) ) then
          !
          ! zero internal mass
          !
          if  (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) )  ) then
             tab_erreur_par(1)%a_imprimer = .true.
             tab_erreur_par(1)%chaine = 'In function generic_function_3p.f90:'
             tab_erreur_par(2)%a_imprimer = .true.
             tab_erreur_par(2)%chaine = 'function call with all arguments zero!'
             call catch_exception(1)
             !
             ! case with one internal mass
             !
          else if (  ( .not.( mz(3) ) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             call tri_int3(z_param_ini,z_param_out)
             !
             f3p_np2_sc_r = f3p0m_1mi_np2(mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( .not.( mz(1) ) ) .and. ( mz(2) ) .and. ( mz(3) ) ) then
             !
             call shift_param(z_param_ini,2,3,z_param_out)
             !
             f3p_np2_sc_r = f3p0m_1mi_np2(mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( .not.( mz(2) ) ) .and. ( mz(1) ) .and. ( mz(3) ) ) then
             !
             call shift_param(z_param_ini,1,3,z_param_out)
             !
             f3p_np2_sc_r = f3p0m_1mi_np2(mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else
             !
             ! comment 11.08.10: only one internal mass possible kinematically!
             tab_erreur_par(1)%a_imprimer = .true.
             tab_erreur_par(1)%chaine = 'In function generic_function_3p.f90:'
             tab_erreur_par(2)%a_imprimer = .true.
             tab_erreur_par(2)%chaine = 'function call with one lightlike, two massive external legs &
                  & and more than one internal mass!'
             tab_erreur_par(3)%a_imprimer = .true.
             tab_erreur_par(3)%chaine = 'This should not be allowed kinematically!'
             call catch_exception(1)
             ! 
             !
             call tri_int3(z_param_ini,z_param_out)
             f3p_np2_sc_r=f3p_finite("n+2",s1,s2,s3,mass1,mass2,mass3, &
                  & z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          !
          ! case with one external leg off shell
          !
       else if ( ( argz(1) ) .and. ( argz(2) ) ) then
          !
          !
          call tri_int3(z_param_ini,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_r = f3p1m_np2(arg3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass
             !
          else if (  ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_r = f3p1m_1mi_np2(arg3,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(2) ) .and. ( mz(3) ) ) then
             !
             call exchange_param(z_param_ini,(/1,3/),3,z_param_out)
             !
             f3p_np2_sc_r = f3p1m_1mi_np2(arg3,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(1) ) .and. ( mz(3) ) ) then
             !
             ! comment 11.08.10: this triangle is 'finite' (in 4dim).
             f3p_np2_sc_r = f3p_finite("n+2",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with two internal masses
             !
          else if ( mz(2) ) then
             !
             f3p_np2_sc_r = f3p1m_2mi_np2(arg3,mass1,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else
             f3p_np2_sc_r=f3p_finite("n+2",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          !
       else if ( (argz(2)) .and. (argz(3)) ) then
          !
          call shift_param(z_param_ini,2,3,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_r = f3p1m_np2(arg1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass
             !
          else if (  ( mz(2) ) .and. ( mz(3) ) ) then
             !
             f3p_np2_sc_r = f3p1m_1mi_np2(arg1,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(3) ) .and. ( mz(1) ) ) then
             !
             ! changed 11.08.10: (/1,2/) -> (/2,3/)
             !
             call exchange_param(z_param_ini,(/2,3/),3,z_param_out)
             !
             f3p_np2_sc_r = f3p1m_1mi_np2(arg1,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(2) ) .and. ( mz(1) ) ) then
             !
             ! comment 11.08.10: this triangle is 'finite' (in 4dim)
             f3p_np2_sc_r=f3p_finite("n+2",s2,s3,s1,mass2,mass3,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with two internal masses
             !
          else if ( mz(3) ) then
             ! comment 11.08.10: masses swapped!
             f3p_np2_sc_r = f3p1m_2mi_np2(arg1,mass2,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else
             f3p_np2_sc_r=f3p_finite("n+2",s2,s3,s1,mass2,mass3,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          !
       else if ( (argz(1)) .and. (argz(3)) ) then
          !
          call shift_param(z_param_ini,1,3,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_r = f3p1m_np2(arg2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal masse
             !
          else if (  ( mz(1) ) .and. ( mz(3) ) ) then
             !
             f3p_np2_sc_r = f3p1m_1mi_np2(arg2,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(2) ) .and. ( mz(1) ) ) then
             !
             call exchange_param(z_param_ini,(/1,2/),3,z_param_out)
             !
             f3p_np2_sc_r = f3p1m_1mi_np2(arg2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else if (  ( mz(2) ) .and. ( mz(3) ) ) then
             !
             ! comment 11.08.10: this triangle is 'finite' (in 4dim)
             f3p_np2_sc_r=f3p_finite("n+2",s3,s1,s2,mass3,mass1,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with two internal masses
             !
          else if ( mz(1) ) then
             !
             ! comment 11.08.10: masses swapped!
             f3p_np2_sc_r = f3p1m_2mi_np2(arg2,mass3,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else
             f3p_np2_sc_r=f3p_finite("n+2",s3,s1,s2,mass3,mass1,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          !
          ! case with two external legs off shell
          !
       else if ( (argz(1)) .and. (.not.(argz(2))) .and. (.not.(argz(3))) ) then
          !
          call tri_int3(z_param_ini,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_r = f3p2m_np2(arg2,arg3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass
             !
          else if (  ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_r = f3p2m_1mi_np2(arg2,arg3,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else
             f3p_np2_sc_r=f3p_finite("n+2",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
       else if ( (argz(2)) .and. (.not.(argz(1))) .and. (.not.(argz(3))) ) then
          !
          call shift_param(z_param_ini,2,3,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_r = f3p2m_np2(arg3,arg1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass
             !
          else if (  ( mz(2) ) .and. ( mz(3) ) ) then
             !
             f3p_np2_sc_r = f3p2m_1mi_np2(arg3,arg1,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else
             f3p_np2_sc_r=f3p_finite("n+2",s2,s3,s1,mass2,mass3,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          !
       else if ( (argz(3)) .and. (.not.(argz(1))) .and. (.not.(argz(2))) ) then
          !
          call shift_param(z_param_ini,1,3,z_param_out)
          !
          if (  ( mz(3) ) .and. ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_r = f3p2m_np2(arg1,arg2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
             ! case with one internal mass
             !
          else if (  ( mz(1) ) .and. ( mz(3) ) ) then
             !
             f3p_np2_sc_r = f3p2m_1mi_np2(arg1,arg2,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          else
             f3p_np2_sc_r=f3p_finite("n+2",s3,s1,s2,mass3,mass1,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if
          !
          !
       else if ( .not.(argz(1)) .and. (.not.(argz(2))) .and. (.not.(argz(3))) ) then
          !
          if (  ( mz(1) ) .and. ( mz(2) )  .and. ( mz(3) ) ) then

             if ( (abs(arg3) >= abs(arg1)) .and. (abs(arg3) >= abs(arg2)) ) then
                !
                call tri_int3(z_param_ini,z_param_out)
                f3p_np2_sc_r = f3p3m("n+2",arg3,arg1,arg2, &
                     z_param_out(1),z_param_out(2),z_param_out(3))
                !
             else if ( (abs(arg1) >= abs(arg3)) .and. (abs(arg1) >= abs(arg2)) ) then
                !
                call shift_param(z_param_ini,2,3,z_param_out)
                f3p_np2_sc_r = f3p3m("n+2",arg1,arg2,arg3, &
                     z_param_out(1),z_param_out(2),z_param_out(3))
                !
             else if ( (abs(arg2) >= abs(arg3)) .and. (abs(arg2) >= abs(arg1)) ) then
                !
                call shift_param(z_param_ini,1,3,z_param_out)
                f3p_np2_sc_r = f3p3m("n+2",arg2,arg3,arg1, &
                     z_param_out(1),z_param_out(2),z_param_out(3))
                !
             end if
             !
          else
             call tri_int3(z_param_ini,z_param_out)
             f3p_np2_sc_r=f3p_finite("n+2",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             !
          end if ! end scan internal masses
          !
       end if ! end scan argi
       !
    end if  ! end if par1=-1...
    !
  end function f3p_np2_sc_r
  !
  function f3p_np2_sc_c(s_mat_c,s,parf1,parf2,parf3)
    !
    complex(ki), intent (in), dimension(:,:) :: s_mat_c
    integer, intent (in), dimension(3) :: s
    integer, intent (in), optional :: parf1,parf2,parf3
    real(ki),dimension(4) :: f3p_np2_sc_c
    !
    integer :: par1,par2,par3
    integer, dimension(3) :: z_param_ini,z_param_out
    complex(ki) :: arg1,arg2,arg3
    complex(ki) :: mass1, mass2, mass3
    real(ki) :: s1,s2,s3
    integer :: m1,m2,m3
    logical, dimension(3) :: argz,mz
    logical :: finite = .true.
    !
    !
    par1 = 0
    par2 = 0
    par3 = 0
    !
    if (present(parf1)) par1 = parf1 
    if (present(parf2)) par2 = parf2
    if (present(parf3)) par3 = parf3
    !
    !
    if ( (par1 == -1) .or. (par2 == -1) .or. (par3 == -1) ) then
       !
       f3p_np2_sc_c(:) = 0._ki
       !
    else
       ! symetrie: la place de z1,z2,z3 n'a pas d'importance, on les met 
       ! dans l'ordre croissant
       z_param_ini(1) = par1
       z_param_ini(2) = par2
       z_param_ini(3) = par3
       !
       m1 = s(1)
       m2 = s(2)
       m3 = s(3)
       !
      !
       arg1 = s_mat_c(m1,m2)
       arg2 = s_mat_c(m2,m3)
       arg3 = s_mat_c(m1,m3)
       !
       argz(1) = equal_real(real(arg1,ki),zero) .and. equal_real(aimag(arg1),zero) 
       argz(2) = equal_real(real(arg2,ki),zero) .and. equal_real(aimag(arg2),zero) 
       argz(3) = equal_real(real(arg3,ki),zero) .and. equal_real(aimag(arg3),zero) 
       !
       ! internal masses
       mass1 = -s_mat_c(m1,m1)/2._ki
       mass2 = -s_mat_c(m2,m2)/2._ki
       mass3 = -s_mat_c(m3,m3)/2._ki
       !
       mz(1) = equal_real(real(mass1,ki),zero) .and. equal_real(aimag(mass1),zero)
       mz(2) = equal_real(real(mass2,ki),zero) .and. equal_real(aimag(mass2),zero)
       mz(3) = equal_real(real(mass3,ki),zero) .and. equal_real(aimag(mass3),zero)
       !
       ! external p_i^2
       s1 = real(arg3+mass1+mass3,ki)
       s2 = real(arg1+mass1+mass2,ki)
       s3 = real(arg2+mass2+mass3,ki)
       !
       call cut_s(s1,mass1,mass3)
       call cut_s(s2,mass1,mass2)
       call cut_s(s3,mass2,mass3)
       !
       ! initialize all components
       f3p_np2_sc_c(:) = 0._ki
       !
       !
       finite = .true.
       ! Similar to 4dim case, there is only one triangle with singular S-matrix.
       ! case with two external legs off shell
       !
       if ( (argz(1)) .and. (.not.(argz(2))) .and. (.not.(argz(3))) ) then
          !
          call tri_int3(z_param_ini,z_param_out)
          !
          ! case with one internal mass
          !
          if (  ( mz(1) ) .and. ( mz(2) ) ) then
             !
             f3p_np2_sc_c = f3p2m_1mi_np2(arg2,arg3,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
             finite = .false.
             !
          end if
          !
       else if ( (argz(2)) .and. (.not.(argz(1))) .and. (.not.(argz(3))) ) then
          !
          call shift_param(z_param_ini,2,3,z_param_out)
          !
          ! case with one internal mass
          !
          if (  ( mz(2) ) .and. ( mz(3) ) ) then
             !
             f3p_np2_sc_c = f3p2m_1mi_np2(arg3,arg1,mass1,z_param_out(1),z_param_out(2),z_param_out(3))
             finite = .false.
             !
          end if
          !
          !
       else if ( (argz(3)) .and. (.not.(argz(1))) .and. (.not.(argz(2))) ) then
          !
          call shift_param(z_param_ini,1,3,z_param_out)
          !
          ! case with one internal mass
          !
          if (  ( mz(1) ) .and. ( mz(3) ) ) then
             !
             f3p_np2_sc_c = f3p2m_1mi_np2(arg1,arg2,mass2,z_param_out(1),z_param_out(2),z_param_out(3))
             finite = .false.
             !
          end if
          !
       end if !args
       !
       if (finite) then ! call to finite triangle 
          !
          call tri_int3(z_param_ini,z_param_out)
          f3p_np2_sc_c=f3p_finite("n+2",s1,s2,s3,mass1,mass2,mass3,z_param_out(1),z_param_out(2),z_param_out(3))
          !
       end if ! end finite
       !
    end if  ! end if par1=-1...
    !
  end function f3p_np2_sc_c
  !
end module generic_function_3p
