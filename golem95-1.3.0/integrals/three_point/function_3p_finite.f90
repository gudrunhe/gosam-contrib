! 
!~ 13.5.2011: no need for LT, avh_olo any longer, 
!~ use own 3-point functions by JPhi,Eric
!
!****h* src/integral/three_point/function_3pC0i
! NAME
!
!  Module function_3p_finite
!
! USAGE
!
!  use function_3p_finite
!
! DESCRIPTION
!
!  This module is used to compute IR finite three point functions
!  with/without Feynman parameters in n, n+2 dimensions
!
! OUTPUT
!
!  This module exports the functions:
!  * f3p_finite, C0  -- functions for the computation of IR finite 
!  three-point functions with/without Feynman parameters in n, n+2 dimensions
!  * f3p_finite_c -- a function which computes the same thing as f3p_finite, only 
!    the format of the return values is different
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * numerical_evaluation (src/numerical/mod_numeric.f90)
!  * dilogarithme (src/module/zdilog.f90)
!  * logarithme (src/module/z_log.f90)
!  * constante (src/module/constante.f90)
!  * parametre (src/module/parametre.f90)
!  * array (src/module/array.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * generic_function_2p (src/integrals/two_point/generic_function_2p.f90)
!  * multiply_div (src/module/multiply_div.f90)
!  * s_matrix_type (src/module/s_matrix_type.f90)
!
!*****
module function_3p_finite
  !
  use precision_golem
  use numerical_evaluation
  use dilogarithme
  use logarithme
  use constante
  use parametre
  use array
  use equal
  use sortie_erreur
  use generic_function_2p
  use multiply_div
  use s_matrix_type
  use matrice_s, only : prepare_s_matrix_local
  implicit none
  !
  private 
  !
  real(ki) :: eps_glob
  integer :: par1_glob,par2_glob,par3_glob
  character (len=3) :: dim_glob
  !
  real(ki), dimension(3) :: b_real
  real(ki) :: sumb_real
  real(ki), dimension(3,3) :: invs_real, s_mat_real
  !
  complex(ki), dimension(3) :: b_complex
  complex(ki) :: sumb_complex
  complex(ki), dimension(3,3) :: invs_complex, s_mat_complex
  !
  type (s_matrix_poly) :: s_mat_p_loc
  !
  integer, dimension(3) :: par
  integer, dimension(3) :: s = (/1,2,3/)
  !
  logical, dimension(:), allocatable :: deja_calcule
  real(ki),dimension(:,:), allocatable :: resultat
  logical, dimension(:), allocatable :: deja_calcule_c
  complex(ki), dimension(:,:), allocatable :: resultat_c
  !
  logical, dimension(:,:), allocatable :: deja_calcule2
  real(ki),dimension(:,:,:), allocatable :: resultat2
  logical, dimension(:,:), allocatable :: deja_calcule2_c
  complex(ki),dimension(:,:,:), allocatable :: resultat2_c
  !
  logical, dimension(:), allocatable :: deja_calcule_np2
  real(ki),dimension(:,:), allocatable :: resultat_np2
  logical, dimension(:), allocatable :: deja_calcule_np2_c
  complex(ki),dimension(:,:), allocatable :: resultat_np2_c
  !
  logical, dimension(:,:,:), allocatable :: deja_calcule22
  real(ki),dimension(:,:,:,:), allocatable :: resultat22
  logical, dimension(:,:,:), allocatable :: deja_calcule22_c
  complex(ki),dimension(:,:,:,:), allocatable :: resultat22_c
  !
  ! added by jpg 30/04/2011
  real(ki) :: delta_glob                !small imaginary part
  complex(ki) :: det_s,num_b1,num_b2,num_b3,num_sumb
  real(ki) :: det_g
  complex(ki), dimension(3,3) :: g_glob
  complex(ki) :: a,b,c,d,e,f
  logical, dimension(3) :: tab_test_glob
  logical, dimension(3) :: tab_test_g2_glob
  !
  !
  interface f3p_finite
     !
     module procedure f3p_finite_rarg
     module procedure f3p_finite_carg
     !
  end interface
  !
  interface C0
     !
     module procedure C0_rarg
     module procedure C0_carg
     !
  end interface
  !
  public :: f3p_finite,f3p_finite_c,C0
  !
  !
 contains
  !
  !****f* src/integral/three_point/function_3p_finite/f3p_finite
  ! NAME
  !
  !  Function f3p_finite
  !
  ! USAGE
  !
  !  real_dim4 = f3p_finite(dim,s1,s2,s3,m1,m2,m3,par1,par2,par3)
  !
  ! DESCRIPTION
  !
  !  This function computes the IR finite three-point 
  !  function in n and n+2 dimensions. 
  !
  ! INPUTS
  !
  !  * s1 -- a real (type ki), p1^2
  !  * s2 -- a real (type ki), p2^2
  !  * s3 -- a real (type ki), p3^2
  !  * m1 -- a real/complex (type ki), the first mass squared
  !  * m2 -- a real/complex (type ki), the second mass squared
  !  * m3 -- a real/complex (type ki), the third mass squared
  !  * par1 -- an integer, the label of the third Feynman parameter
  !  * par2 -- an integer, the label of the second Feynman parameter
  !  * par3 -- an integer, the label of the first Feynman parameter
  !  Note that par1,par2 and par3 are supposed to be ordered, i.e.
  !  par1 <= par2 <= par3, note also that put zero for par1, par2 or par3
  !  if this Feynman parameter does not exist.
  !  Use the routine tri_int(t_in,t_out) to order the labels in the module 
  !  tri_croissant (src/module/tri.f90)
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  An real (type ki) array of rank 1 and shape 6,
  !  where the last two entries corresponding to 
  !  the real/imaginary part of the constant term. 
  !  the first 4 entries are always zero, but the shape should be 
  !  uniform for all triangles called in generic_function_3p
  !
  !*****
  function f3p_finite_rarg(dim,s1,s2,s3,m1,m2,m3,par1,par2,par3)
    implicit none
    !
    character (len=3), intent (in) :: dim
    real(ki), intent (in) :: s1,s2,s3,m1,m2,m3
    integer, intent (in) :: par1,par2,par3
    real(ki), dimension(4) :: f3p_finite_rarg
    !
    integer :: nb_par
    real(ki) :: lamb, detS3
    real(ki) :: plus_grand
    real(ki) :: norma
    complex(ki) :: resto,abserro
    complex(ki), dimension(3,3) :: s_mat_loc
    integer :: i
    real(ki) :: s1p,s2p,s3p,m1p,m2p,m3p
    !
    !write(6,*) " MR: entering f3p_finite_rarg  ",s1," s2 ",s2," s3 ",s3," m1 ",m1," m2 ",m2," m3 ",m3
    par = (/par1,par2,par3/)
    !
    s_mat_real(1,:) = (/-m1*2._ki,s2-m1-m2,s1-m1-m3/)
    s_mat_real(2,:) = (/s2-m1-m2,-m2*2._ki,s3-m2-m3/)
    s_mat_real(3,:) = (/s1-m1-m3,s3-m2-m3,-m3*2._ki/)
    !
    ! on redefinit la matrice S de telle facon a ce que ces elements
    ! soient entre -1 et 1
    ! this can't be done if partly invariants instead of 
    ! s_mat_real are used as arguments
    !
    if (rat_or_tot_par%tot_selected) then
    !
      plus_grand = maxval(array=abs(s_mat_real))
    !
    else if (rat_or_tot_par%rat_selected) then
    !
      plus_grand = 1._ki
    !
    end if
    !
    !
    s_mat_real = s_mat_real/plus_grand
    !
    s_mat_p_loc = assign_s_matrix(s_mat_real)
    call prepare_s_matrix_local(s_mat_p_loc,s)
    !
    detS3 = -(s_mat_real(1, 3)**2*s_mat_real(2, 2)) + 2*s_mat_real(1, 2)*s_mat_real(1, 3)*s_mat_real(2, 3) -  &
         &     s_mat_real(1, 2)**2*s_mat_real(3, 3) + s_mat_real(1, 1)*(-s_mat_real(2, 3)**2 +   &
         &     s_mat_real(2, 2)*s_mat_real(3, 3))

    b_real(1) = (-s_mat_real(2, 3)**2 + s_mat_real(1, 3)*(-s_mat_real(2, 2) + s_mat_real(2, 3)) +     &
         &      s_mat_real(1, 2)*(s_mat_real(2, 3) - s_mat_real(3, 3)) + s_mat_real(2, 2)*s_mat_real(3, 3))/detS3
    b_real(2) = (-s_mat_real(1, 3)**2 + s_mat_real(1, 3)*s_mat_real(2, 3) +  &
         &    s_mat_real(1, 2)*(s_mat_real(1, 3) - s_mat_real(3, 3)) +   &
         &    s_mat_real(1, 1)*(-s_mat_real(2, 3) + s_mat_real(3, 3)))/detS3
    b_real(3) = (-s_mat_real(1, 2)**2 - s_mat_real(1, 3)*s_mat_real(2, 2) +   &
         &    s_mat_real(1, 1)*(s_mat_real(2, 2) - s_mat_real(2, 3)) +      &
         &    s_mat_real(1, 2)*(s_mat_real(1, 3) + s_mat_real(2, 3)))/detS3
    !
    sumb_real = (-s_mat_real(1, 2)**2 - s_mat_real(1, 3)**2 + s_mat_real(1, 1)*s_mat_real(2, 2) - &
         &   2*s_mat_real(1, 3)*(s_mat_real(2, 2) - s_mat_real(2, 3)) - 2*s_mat_real(1, 1)*s_mat_real(2, 3) - &
         &   s_mat_real(2, 3)**2 + 2*s_mat_real(1, 2)*(s_mat_real(1, 3) + s_mat_real(2, 3) - s_mat_real(3, 3)) + &
         &   s_mat_real(1, 1)*s_mat_real(3, 3) + s_mat_real(2, 2)*s_mat_real(3, 3))/detS3
    !
    invs_real(1,1) = (-s_mat_real(2, 3)**2 + s_mat_real(2, 2)*s_mat_real(3, 3))/detS3
    invs_real(1,2) = (s_mat_real(1, 3)*s_mat_real(2, 3) - s_mat_real(1, 2)*s_mat_real(3, 3))/detS3
    invs_real(1,3) = (-(s_mat_real(1, 3)*s_mat_real(2, 2)) + s_mat_real(1, 2)*s_mat_real(2, 3))/detS3
    invs_real(2,1) = invs_real(1,2)
    invs_real(2,2) = (-s_mat_real(1, 3)**2 + s_mat_real(1, 1)*s_mat_real(3, 3))/detS3
    invs_real(2,3) = (s_mat_real(1, 2)*s_mat_real(1, 3) - s_mat_real(1, 1)*s_mat_real(2, 3))/detS3
    invs_real(3,1) = invs_real(1,3)
    invs_real(3,2) = invs_real(2,3)
    invs_real(3,3) = (-s_mat_real(1, 2)**2 + s_mat_real(1, 1)*s_mat_real(2, 2))/detS3
    !
    !lamb = 2._ki*s_mat_real(1,3)*s_mat_real(2,3)+2._ki*s_mat_real(1,2)*s_mat_real(2,3)&
         !+2._ki*s_mat_real(1,3)*s_mat_real(1,2)-s_mat_real(1,3)*s_mat_real(1,3)-s_mat_real(1,2)*s_mat_real(1,2)&
         !-s_mat_real(2,3)*s_mat_real(2,3)
    lamb = (-s_mat_real(1, 2)**2 - s_mat_real(1, 3)**2 + s_mat_real(1, 1)*s_mat_real(2, 2) - &
         &   2*s_mat_real(1, 3)*(s_mat_real(2, 2) - s_mat_real(2, 3)) - 2*s_mat_real(1, 1)*s_mat_real(2, 3) - &
         &   s_mat_real(2, 3)**2 + 2*s_mat_real(1, 2)*(s_mat_real(1, 3) + s_mat_real(2, 3) - s_mat_real(3, 3)) + &
         &   s_mat_real(1, 1)*s_mat_real(3, 3) + s_mat_real(2, 2)*s_mat_real(3, 3))
    !
    nb_par = count(mask=par/=0)
    !
    s_mat_loc = cmplx(s_mat_real,0._ki,ki)
    !
    m1p = -s_mat_loc(1,1)/2._ki
    m2p = -s_mat_loc(2,2)/2._ki
    m3p = -s_mat_loc(3,3)/2._ki
    s1p = s_mat_loc(1,3) - (s_mat_loc(1,1)+s_mat_loc(3,3))/2._ki
    s2p = s_mat_loc(1,2) - (s_mat_loc(1,1)+s_mat_loc(2,2))/2._ki
    s3p = s_mat_loc(2,3) - (s_mat_loc(2,2)+s_mat_loc(3,3))/2._ki
    !
    a = s_mat_loc(1,3)-1._ki/2._ki*s_mat_loc(1,1)-1._ki/2._ki*s_mat_loc(3,3)
    b = s_mat_loc(2,3)-1._ki/2._ki*s_mat_loc(2,2)-1._ki/2._ki*s_mat_loc(3,3)
    c = s_mat_loc(1,2)-s_mat_loc(1,3)-s_mat_loc(2,3)+s_mat_loc(3,3)
    d = -s_mat_loc(1,3)+s_mat_loc(1,1)
    e = s_mat_loc(1,3)-s_mat_loc(1,2)
    f = -1._ki/2._ki*s_mat_loc(1,1)
    !
      !write(*,*) 'test int_log 0:',b,s3
    !
    ! les coefficients des differentes formes quadratiques
    g_glob(1,1) = b
    g_glob(1,2) = c+e
    g_glob(1,3) = a+d+f
    g_glob(2,1) = a
    g_glob(2,2) = d
    g_glob(2,3) = f
    g_glob(3,1) = a+b+c
    g_glob(3,2) = d+e
    g_glob(3,3) = f
    !
    do i=1,3
      tab_test_glob(i) =  equal_real(real(g_glob(i,1),ki),0._ki,1.e+1_ki) .and. equal_real(abs(g_glob(i,2)),0._ki,1.e5_ki)
    end do
    do i=1,3
      !write(*,*) 'test tab_test :',real(g_glob(i,1),ki),abs(g_glob(i,2))
      !write(*,*) 'test tab_test :',equal_real(real(g_glob(i,1),ki),0._ki,1.e+1_ki),equal_real(abs(g_glob(i,2)),0._ki,1.e5_ki)
      tab_test_g2_glob(i) =  equal_real(real(g_glob(i,1),ki),0._ki,1.e+1_ki) .and. .not.(equal_real(abs(g_glob(i,2)),0._ki,1.e5_ki))
    end do
    !~ delta_glob = tiny(1._ki)
    delta_glob = epsilon(1._ki)
    !
    call compute_deter(s_mat_loc(1,3),s_mat_loc(1,2),s_mat_loc(2,3),s_mat_loc(1,1),s_mat_loc(2,2),s_mat_loc(3,3))
    !
    if (nb_par == 0) then
       !
       norma = -1._ki/2._ki
       !
    else if (nb_par == 1) then
       !
       norma = -1._ki/6._ki
       !
    else if (nb_par == 2) then
       !
       ! Different normalisations: depends on whether the feynman parameters are equal or not 
       ! Use the fact that the two parameters will be par2 and par3
       if (par2==par3) then
          norma = -1._ki/12._ki
       else   
          norma = -1._ki/24._ki
       endif
       !
    else
       !
       norma = 0._ki
       !
    end if
    !
    ! memory allocation to save time in the recursion
    !
    allocate(deja_calcule(4))
    allocate(resultat(4,2))
    allocate(deja_calcule2(3,4))
    allocate(resultat2(3,4,4))
    allocate(deja_calcule_np2(4))
    allocate(resultat_np2(4,4))
    allocate(deja_calcule22(3,4,4))
    allocate(resultat22(3,4,4,4))
    !
    ! initialisation
    !
    deja_calcule = .false.
    resultat = 0._ki
    deja_calcule2 = .false.
    resultat2 = 0._ki
    deja_calcule_np2 = .false.
    resultat_np2 = 0._ki
    deja_calcule22 = .false.
    resultat22 = 0._ki
    !
    !
    f3p_finite_rarg = 0._ki
    !
       !write(*,*) 'test f3p_finite_rarg 0 :',count(mask=tab_test_glob),abs(sumb_real),lamb
       !write(*,*) 'test f3p_finite_rarg 0 :',&
       !&detS3,b_real(1)*detS3,b_real(2)*detS3,b_real(3)*detS3
       !write(*,*) 'test f3p_finite_rarg 0 :',sumb_real*detS3
     !write(*,*) 'test f3p_finite :',abs(sumb_real),coupure_3p3m
    !~ write(*,*) 'test f3p_finite :',abs(sumb_real),coupure_3p3m
    !if ( (abs(sumb_real) > coupure_3p3m) .or. (count(mask=tab_test_glob) == 2) ) then
    !if ( (abs(lamb) > coupure_3p3m) .or. (count(mask=tab_test_glob) == 2) ) then
    if ( (abs(sumb_real) > coupure_3p3m) .or. (count(mask=tab_test_glob) == 2) ) then
    !write(6,*) " MR: f3p_finite_rarg in coup3p3m if "
     !if (abs(sumb_real) < 0._ki) then
       !
       !write(*,*) 'test f3p_finite_rarg 1 :',count(mask=tab_test_glob)
       ! always use analytic computation until massive numerical is implemented
        ! branching removed completely Feb 4, 2011 because if detS3=0, 
       ! sumb_real=NAN so numerical branch was entered, although scalar integral 
       ! still finite if detS3=0 without any z_i=0 (LLS)
       !
      ! when invariants as input are used, no division by plus_grand!
       !
       if (dim == "ndi") then
          !
     !write(6,*) " MR: f3p_finite_rarg start a3pC0i_rarg "
         !f3p_finite_rarg(3:4)= a3pC0i_rarg(s1,s2,s3,m1,m2,m3,par1,par2,par3)&
         !   &/plus_grand
        f3p_finite_rarg(3:4)= a3pC0i_rarg(s1p,s2p,s3p,m1p,m2p,m3p,par1,par2,par3)&
           &/plus_grand
      !write(6,*) " MR: f3p_finite_rarg a3pC0i_rarg done ",f3p_finite_rarg(3:4)
         !
       else if (dim == "n+2") then
          !
          f3p_finite_rarg = a3pC0i_np2_rarg(s1p,s2p,s3p,m1p,m2p,m3p,par1,par2,par3)
          !f3p_finite_rarg = a3pC0i_np2_rarg(s1,s2,s3,m1,m2,m3,par1,par2,par3)
          f3p_finite_rarg(3) = f3p_finite_rarg(3)-log(plus_grand)*norma
	  ! mu2_scale_par is already contained in the bubbles, 
	  ! scaling of s_mat_real  does not enter here
          !
       end if
       !
    else
     !
     !write(6,*) " MR: f3p_finite_rarg in coup3p3m else, num "
     ! numerical computation
     !
     dim_glob = dim
     par1_glob = par1
     par2_glob = par2
     par3_glob = par3
     !
       resto = 0._ki
       abserro = 0._ki
       !
       origine_info_par = "f3p_finite, dimension "//dim
       num_grand_b_info_par = lamb
       denom_grand_b_info_par = detS3
      !
      call generic_eval_numer(eval_numer_g0,0._ki,1._ki,tolerance,resto,abserro)
      !
       if (dim == "ndi") then      
          !
           resto = resto/plus_grand
          !
        else if (dim == "n+2") then
          !
          f3p_finite_rarg(1) = norma
          f3p_finite_rarg(2) = 0._ki
          resto = resto-log(plus_grand/mu2_scale_par)*norma
          !
        end if
       !
        f3p_finite_rarg(3) = real(resto,ki)
        f3p_finite_rarg(4) = aimag(resto)
       !
    end if  ! end if analytic or numeric
    !
    deallocate(deja_calcule)
    deallocate(resultat)
    deallocate(deja_calcule2)
    deallocate(resultat2)
    deallocate(deja_calcule_np2)
    deallocate(resultat_np2)
    deallocate(deja_calcule22)
    deallocate(resultat22)
    !
    call nullify_s_matrix(s_mat_p_loc)
    !
     !write(6,*) " MR: leaving f3p_finite_rarg ",f3p_finite_rarg

  end function f3p_finite_rarg
  !
  function f3p_finite_carg(dim,s1r,s2r,s3r,m1,m2,m3,par1,par2,par3)
    implicit none
    !
    character (len=3), intent (in) :: dim
    real(ki), intent (in) :: s1r, s2r, s3r
    complex(ki), intent (in) :: m1,m2,m3
    integer, intent (in) :: par1,par2,par3
    real(ki), dimension(4) :: f3p_finite_carg
    !
    complex(ki) :: s1, s2, s3
    integer :: nb_par
    complex(ki) :: lamb, detS3
    real(ki) :: plus_grand
    complex(ki) :: norma
    complex(ki) :: resto,abserro
    complex(ki) :: temp0
    complex(ki), dimension(2) :: temp
    complex(ki), dimension(3,3) :: s_mat_loc
    real(ki) :: s1rp,s2rp,s3rp
    complex(ki) :: m1p,m2p,m3p
    integer :: i
    !
     !write(6,*) " MR: entering f3p_finite_carg s1r ",s1r," s2r ",s2r," s3r ",s3r," m1 ",m1," m2 ",m2," m3 ",m3

    par = (/par1,par2,par3/)
    !
    s1 = cmplx(s1r, 0._ki,ki)
    s2 = cmplx(s2r, 0._ki,ki)
    s3 = cmplx(s3r, 0._ki,ki)
    !
    s_mat_complex(1,:) = (/-m1*2._ki,s2-m1-m2,s1-m1-m3/)
    s_mat_complex(2,:) = (/s2-m1-m2,-m2*2._ki,s3-m2-m3/)
    s_mat_complex(3,:) = (/s1-m1-m3,s3-m2-m3,-m3*2._ki/)
    !
    ! on redefinit la matrice S de telle facon a ce que ces elements
    ! soient entre -1 et 1
    ! this can't be done if partly invariants instead of 
    ! s_mat_complex are used as arguments
    !
    !~ plus_grand = 1._ki
    if (rat_or_tot_par%tot_selected) then
    !
     plus_grand = maxval(array=abs(s_mat_complex))
    !
    else if (rat_or_tot_par%rat_selected) then
    !
      plus_grand = 1._ki
    !
    end if
    !
    s_mat_complex = s_mat_complex/plus_grand
    !
    s_mat_real = real(s_mat_complex,ki)
    !
    s_mat_p_loc = assign_s_matrix(s_mat_complex,s_mat_real)
    call prepare_s_matrix_local(s_mat_p_loc, s)
    !
    detS3 = -(s_mat_complex(1, 3)**2*s_mat_complex(2, 2)) + 2*s_mat_complex(1, 2)*s_mat_complex(1, 3)*s_mat_complex(2, 3) -  &
         &     s_mat_complex(1, 2)**2*s_mat_complex(3, 3) + s_mat_complex(1, 1)*(-s_mat_complex(2, 3)**2 +   &
         &     s_mat_complex(2, 2)*s_mat_complex(3, 3))

    b_complex(1) = (-s_mat_complex(2, 3)**2 + s_mat_complex(1, 3)*(-s_mat_complex(2, 2) + s_mat_complex(2, 3)) +     &
         &      s_mat_complex(1, 2)*(s_mat_complex(2, 3) - s_mat_complex(3, 3)) + s_mat_complex(2, 2)*s_mat_complex(3, 3))/detS3
    b_complex(2) = (-s_mat_complex(1, 3)**2 + s_mat_complex(1, 3)*s_mat_complex(2, 3) +  &
         &    s_mat_complex(1, 2)*(s_mat_complex(1, 3) - s_mat_complex(3, 3)) +   &
         &    s_mat_complex(1, 1)*(-s_mat_complex(2, 3) + s_mat_complex(3, 3)))/detS3
    b_complex(3) = (-s_mat_complex(1, 2)**2 - s_mat_complex(1, 3)*s_mat_complex(2, 2) +   &
         &    s_mat_complex(1, 1)*(s_mat_complex(2, 2) - s_mat_complex(2, 3)) +      &
         &    s_mat_complex(1, 2)*(s_mat_complex(1, 3) + s_mat_complex(2, 3)))/detS3
    !
    sumb_complex = (-s_mat_complex(1, 2)**2 - s_mat_complex(1, 3)**2 + s_mat_complex(1, 1)*s_mat_complex(2, 2) - &
         &   2*s_mat_complex(1, 3)*(s_mat_complex(2, 2) - s_mat_complex(2, 3)) - 2*s_mat_complex(1, 1)*s_mat_complex(2, 3) - &
         &   s_mat_complex(2, 3)**2 + 2*s_mat_complex(1, 2)*(s_mat_complex(1, 3) + s_mat_complex(2, 3) - s_mat_complex(3, 3)) + &
         &   s_mat_complex(1, 1)*s_mat_complex(3, 3) + s_mat_complex(2, 2)*s_mat_complex(3, 3))/detS3
    !
    invs_complex(1,1) = (-s_mat_complex(2, 3)**2 + s_mat_complex(2, 2)*s_mat_complex(3, 3))/detS3
    invs_complex(1,2) = (s_mat_complex(1, 3)*s_mat_complex(2, 3) - s_mat_complex(1, 2)*s_mat_complex(3, 3))/detS3
    invs_complex(1,3) = (-(s_mat_complex(1, 3)*s_mat_complex(2, 2)) + s_mat_complex(1, 2)*s_mat_complex(2, 3))/detS3
    invs_complex(2,1) = invs_complex(1,2)
    invs_complex(2,2) = (-s_mat_complex(1, 3)**2 + s_mat_complex(1, 1)*s_mat_complex(3, 3))/detS3
    invs_complex(2,3) = (s_mat_complex(1, 2)*s_mat_complex(1, 3) - s_mat_complex(1, 1)*s_mat_complex(2, 3))/detS3
    invs_complex(3,1) = invs_complex(1,3)
    invs_complex(3,2) = invs_complex(2,3)
    invs_complex(3,3) = (-s_mat_complex(1, 2)**2 + s_mat_complex(1, 1)*s_mat_complex(2, 2))/detS3
    !
    lamb = (-s_mat_complex(1, 2)**2 - s_mat_complex(1, 3)**2 + s_mat_complex(1, 1)*s_mat_complex(2, 2) - &
         &   2*s_mat_complex(1, 3)*(s_mat_complex(2, 2) - s_mat_complex(2, 3)) - 2*s_mat_complex(1, 1)*s_mat_complex(2, 3) - &
         &   s_mat_complex(2, 3)**2 + 2*s_mat_complex(1, 2)*(s_mat_complex(1, 3) + s_mat_complex(2, 3) - s_mat_complex(3, 3)) + &
         &   s_mat_complex(1, 1)*s_mat_complex(3, 3) + s_mat_complex(2, 2)*s_mat_complex(3, 3))
!     lamb = 2._ki*s_mat_complex(1,3)*s_mat_complex(2,3)+2._ki*s_mat_complex(1,2)*s_mat_complex(2,3)&
!         +2._ki*s_mat_complex(1,3)*s_mat_complex(1,2)-s_mat_complex(1,3)*s_mat_complex(1,3)-s_mat_complex(1,2)*s_mat_complex(1,2)&
!         -s_mat_complex(2,3)*s_mat_complex(2,3)
    !
    nb_par = count(mask=par/=0)
    !
    s_mat_loc = s_mat_complex
    !
    m1p = m1/plus_grand
    m2p = m2/plus_grand
    m3p = m3/plus_grand
    s1rp = s1r/plus_grand
    s2rp = s2r/plus_grand
    s3rp = s3r/plus_grand
    !
    a = s_mat_loc(1,3)-1._ki/2._ki*s_mat_loc(1,1)-1._ki/2._ki*s_mat_loc(3,3)
    b = s_mat_loc(2,3)-1._ki/2._ki*s_mat_loc(2,2)-1._ki/2._ki*s_mat_loc(3,3)
    c = s_mat_loc(1,2)-s_mat_loc(1,3)-s_mat_loc(2,3)+s_mat_loc(3,3)
    d = -s_mat_loc(1,3)+s_mat_loc(1,1)
    e = s_mat_loc(1,3)-s_mat_loc(1,2)
    f = -1._ki/2._ki*s_mat_loc(1,1)
    !
    !
    ! les coefficients des differentes formes quadratiques
    g_glob(1,1) = b
    g_glob(1,2) = c+e
    g_glob(1,3) = a+d+f
    g_glob(2,1) = a
    g_glob(2,2) = d
    g_glob(2,3) = f
    g_glob(3,1) = a+b+c
    g_glob(3,2) = d+e
    g_glob(3,3) = f
    !
    do i=1,3
      tab_test_glob(i) =  equal_real(real(g_glob(i,1),ki),0._ki,1.e+1_ki) .and. equal_real(abs(g_glob(i,2)),0._ki,1.e5_ki)
    end do
    do i=1,3
      tab_test_g2_glob(i) =  equal_real(real(g_glob(i,1),ki),0._ki,1.e+1_ki) .and. .not.(equal_real(abs(g_glob(i,2)),0._ki,1.e5_ki))
    end do
    delta_glob = epsilon(1._ki)
    !
    !~ b_inf = 0._ki
    !~ b_sup = 1._ki
    !
    call compute_deter(s_mat_loc(1,3),s_mat_loc(1,2),s_mat_loc(2,3),s_mat_loc(1,1),s_mat_loc(2,2),s_mat_loc(3,3))
    !
    !
    if (nb_par == 0) then
       !
       norma = -1._ki/2._ki
       !
    else if (nb_par == 1) then
       !
       norma = -1._ki/6._ki
       !
    else if (nb_par == 2) then
       !
       ! Different normalisations: depends on whether the feynman parameters are equal or not 
       ! Use the fact that the two parameters will be par2 and par3
       if (par2==par3) then
          norma = -1._ki/12._ki
       else   
          norma = -1._ki/24._ki
       endif
       !
    else
       !
       norma = 0._ki
       !
    end if
    !
    ! memory allocation to save time in the recursion
    !
    allocate(deja_calcule_c(4))
    allocate(resultat_c(4,1) )
    allocate(deja_calcule2_c(3,4))
    allocate(resultat2_c(3,4,2))
    allocate(deja_calcule_np2_c(4))
    allocate(resultat_np2_c(4,2))
    allocate(deja_calcule22_c(3,4,4))
    allocate(resultat22_c(3,4,4,2))
    !
    ! initialisation
    !
    deja_calcule_c = .false.
    resultat_c = czero
    deja_calcule2_c = .false.
    resultat2_c = czero
    deja_calcule_np2_c = .false.
    resultat_np2_c = czero
    deja_calcule22_c = .false.
    resultat22_c = czero
    !
    !
    f3p_finite_carg(:) = 0._ki
    !
     !write(*,*) 'test f3p_finite_carg :',abs(sumb_complex),coupure_3p3m
    if (abs(sumb_complex) > coupure_3p3m) then
    ! if (abs(sumb_complex) > 0._ki) then
       !
    !write(6,*) " MR: f3p_finite_carg in coup3p3m if "
     ! always use analytic computation until massive numerical is implemented
       ! 
       !
       if (dim == "ndi") then
          !
          !temp0 = a3pC0i_carg(s1r,s2r,s3r,m1,m2,m3,par1,par2,par3)
          temp0 = a3pC0i_carg(s1rp,s2rp,s3rp,m1p,m2p,m3p,par1,par2,par3)
          f3p_finite_carg(3) = real(temp0,ki)/plus_grand
          f3p_finite_carg(4) = aimag(temp0)/plus_grand
          !
       else if (dim == "n+2") then
          !
          !temp = a3pC0i_np2_carg(s1r,s2r,s3r,m1,m2,m3,par1,par2,par3)
          temp = a3pC0i_np2_carg(s1rp,s2rp,s3rp,m1p,m2p,m3p,par1,par2,par3)
          f3p_finite_carg(1) = real(temp(1),ki)
          f3p_finite_carg(2) = aimag(temp(1))
          f3p_finite_carg(3) = real(temp(2),ki)
          f3p_finite_carg(4) = aimag(temp(2))
          f3p_finite_carg(3) = f3p_finite_carg(3)-log(plus_grand)*norma
	  ! mu2_scale_par is already contained in the bubbles, 
	  ! scaling of s_mat_complex  does not enter here
          !
       end if
       !
    else
     !
    !write(6,*) " MR: f3p_finite_carg in coup3p3m else, i.e num "
    ! numerical computation
     !
     dim_glob = dim
     par1_glob = par1
     par2_glob = par2
     par3_glob = par3
     !
       resto = 0._ki
       abserro = 0._ki
       !
       origine_info_par = "f3p_finite, dimension "//dim
       num_grand_b_info_par = real(lamb,ki)
       denom_grand_b_info_par = real(detS3,ki)
      !
      call generic_eval_numer(eval_numer_g0,0._ki,1._ki,tolerance,resto,abserro)
      !
       if (dim == "ndi") then      
          !
           resto = resto/plus_grand
           !~ resto = resto
          !
        else if (dim == "n+2") then
          !
          f3p_finite_carg(1) = norma
          f3p_finite_carg(2) = 0._ki
           resto = resto-log(plus_grand/mu2_scale_par)*norma
           !~ resto = resto
          !
        end if
       !
        f3p_finite_carg(3) = real(resto,ki)
        f3p_finite_carg(4) = aimag(resto)
       !
    end if  ! end if analytic or numeric
    !
       ! numerical computation disabled
       !
      deallocate(deja_calcule_c)
      deallocate(resultat_c)
      deallocate(deja_calcule2_c)
      deallocate(resultat2_c)
      deallocate(deja_calcule_np2_c)
      deallocate(resultat_np2_c)
      deallocate(deja_calcule22_c)
      deallocate(resultat22_c)
    !
    call nullify_s_matrix(s_mat_p_loc)
    !
  end function f3p_finite_carg
  !
  !****f* src/integral/three_point/function_3pC0i/f3p_finite_c
  ! NAME
  !
  !  Function f3p_finite_c
  !
  ! USAGE
  !
  !  complex_dim2 = f3p_finite_c(s1,s2,s3,m1,m2,m3,par1,par2,par3)
  !
  ! DESCRIPTION
  !
  !  It computes the same as the function f3p_finite, but the returned
  !  value is a complex (type ki) array of rank 1 and shape 2
  !
  ! INPUTS
  !
  !  * dim -- a character (length 3), to compute in n or n+2 dimensions
  !  * s1 -- a real (type ki), p1^2
  !  * s2 -- a real (type ki), p2^2
  !  * s3 -- a real (type ki), p3^2
  !  * m1 -- a real (type ki), the first mass^2
  !  * m2 -- a real (type ki), the second mass^2
  !  * m3 -- a real (type ki), the third mass^2
  !  * par1 -- an integer, the label of the third Feynman parameter
  !  * par2 -- an integer, the label of the second Feynman parameter
  !  * par3 -- an integer, the label of the first Feynman parameter
  !  Note that par1,par2 and par3 are supposed to be ordered, i.e.
  !  par1 <= par2 <= par3, note also that put zero for par1, par2 or par3
  !  if this Feynman parameter does not exist.
  !  Use the routine tri_int(t_in,t_out) to order the labels in the module 
  !  tri_croissant (src/module/tri.f90)
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  An complex (type ki) array of rank 1 and shape 2 corresponding to 
  !  the (real part,imaginary part) of the coefficient of the 1/epsilon term
  !  and the (real part,imaginary part) of the constant term. if par1 and/or par2
  !  are different from zero for dim="n+2", an error is returned.
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  function f3p_finite_c(dim,s1,s2,s3,m1,m2,m3,par1,par2,par3)
    !
    use translate
    implicit none
    !
    character (len=3), intent (in) :: dim
    real(ki), intent (in) :: s1,s2,s3,m1,m2,m3
    integer, intent (in) :: par1,par2,par3
    complex(ki), dimension(2) :: f3p_finite_c
    !
    real(ki), dimension(4) :: res4
    !
    res4 = f3p_finite(dim,s1,s2,s3,m1,m2,m3,par1,par2,par3)
    call to_complex(res4,f3p_finite_c)
    !
  end function f3p_finite_c
  !
  !****if* src/integral/three_point/function_3pC0i/a3pC0i
  ! NAME
  !
  !  Function a3pC0i
  !
  ! USAGE
  !
  !  real_dim2 = a3pC0i_rarg(s1,s2,s3,m1r,m2r,m3r,par1,par2,par3)
  !  complex_dim1 = a3pC0i_carg(s1,s2,s3,m1c,m2c,m3c,par1,par2,par3)
  !
  ! DESCRIPTION
  !
  !  This recursive function implements the formula of ref 1
  !
  ! INPUTS
  !
  !  * s1 -- a real (type ki), p1^2
  !  * s2 -- a real (type ki), p2^2
  !  * s3 -- a real (type ki), p3^2
  !  * m1(r/c) -- a real/complex (type ki), the first internal mass squared
  !  * m2(r/c) -- a real/complex (type ki), the second internal mass squared
  !  * m3(r/c) -- a real/complex (type ki), the third internal mass squared
  !  * par1 -- an integer, the label of the third Feynman parameter
  !  * par2 -- an integer, the label of the second Feynman parameter
  !  * par3 -- an integer, the label of the first Feynman parameter
  !
  ! SIDE EFFECTS
  !
  !  This function modifies the value of the local (for the module) variables:
  !  * deja_calcule(_c), deja_calcule2, deja_calcule_np2(_c) and deja_calcule22
  !  * resultat(_c), resultat2, resultat_np2(_c) and resultat22
  !
  ! RETURN VALUE
  !
  !  It returns a real/complex (type ki) array of rank 1 and shape 2/1
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  recursive function a3pC0i_rarg(s1,s2,s3,m1,m2,m3,par1,par2,par3) result(res_3pC0i_rarg)
    !
    implicit none
    real(ki), intent (in) :: s1,s2,s3,m1,m2,m3
    integer, intent (in) :: par1,par2,par3
    real(ki), dimension(2) :: res_3pC0i_rarg
    !
    integer :: j
    integer :: nb_par_loc
    integer, dimension(3) :: par_loc,par_plus
    real(ki), dimension(4) :: truc1
    real(ki), dimension(2) :: temp0
    real(ki), dimension(4) :: temp1,temp2,temp3
    real(ki), dimension(4) :: temp10,temp11,temp12
    complex(ki) :: ctemp
    integer :: ib,b_pro,b_pro_mj
    !
    b_pro = packb(s)
    !
    par_loc = (/par1,par2,par3/)
    par_plus = par_loc+1
    nb_par_loc = count(mask=par_loc/=0)
    !
    ! cas sans parametre de feynman au numerateur
    !
    if (nb_par_loc == 0) then
       ctemp = C0(s1,s2,s3,m1,m2,m3)
       res_3pC0i_rarg(1) = real(ctemp,ki)
       res_3pC0i_rarg(2) = aimag(ctemp)
       !
       ! cas avec un parametre de feynman au numerateur
       !
    else if (nb_par_loc == 1) then
       !
       if (deja_calcule(1)) then
          !
          temp0 = resultat(1,:)
          !
       else
          !
          temp0 = a3pC0i_rarg(s1,s2,s3,m1,m2,m3,0,0,0)
          resultat(1,:) = temp0
          deja_calcule(1) = .true.
          !
       end if
       !
       temp1 = 0._ki
       temp2 = 0._ki
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if (deja_calcule2(j,1)) then
                !
                truc1 = resultat2(j,1,:)
                !
             else
                !
                truc1 = f2p_ra(s_mat_p_loc,b_pro_mj) !returns real array!
                resultat2(j,1,:) = truc1
                deja_calcule2(j,1) = .true.
                !
             end if
             !
             temp1 = temp1 + b_real(j)*truc1
             temp2 = temp2 + invs_real(j,par3)*truc1
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_rarg(1) = b_real(par3)*(temp0(1) - temp1(3))/sumb_real + temp2(3)
       res_3pC0i_rarg(2) = b_real(par3)*(temp0(2) - temp1(4))/sumb_real + temp2(4)
       !
       ! cas avec deux parametres de feynman au numerateur
       !
    else if (nb_par_loc == 2) then
       !
       if (deja_calcule_np2(par_plus(3))) then
          !
          temp11 = resultat_np2(par_plus(3),:)
          !
       else
          !
          temp11 = a3pC0i_np2_rarg(s1,s2,s3,m1,m2,m3,0,0,par3)
          resultat_np2(par_plus(3),:) = temp11
          deja_calcule_np2(par_plus(3)) = .true.
          !
       end if
       !
       temp10 = resultat_np2(1,:)
       temp3 = invs_real(par2,par3)*temp10
       temp1 = b_real(par2)*temp11
       temp1 = mult_div(-2._ki/3._ki,temp1)*3._ki
       temp2 = 0._ki
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if (j /= par3) then
                !
                if (deja_calcule2(j,par_plus(3))) then
                   !
                   truc1 = resultat2(j,par_plus(3),:)
                   !
                else
                   !
                   truc1 = f2p_ra(s_mat_p_loc,b_pro_mj,par3) !returns real array!
                   resultat2(j,par_plus(3),:) = truc1
                   deja_calcule2(j,par_plus(3)) = .true.
                   !
                end if
                !
                temp2 = temp2 + invs_real(j,par2)*truc1
                !
             end if
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_rarg(1) = -temp3(3) + temp1(3) + temp2(3)
       res_3pC0i_rarg(2) = -temp3(4) + temp1(4) + temp2(4)
       !
       ! cas avec trois parametres de feynman au numerateur
       !
    else
       !
       temp12 = a3pC0i_np2_rarg(s1,s2,s3,m1,m2,m3,0,par2,par3)
       temp10 = resultat_np2(par_plus(3),:)
       temp11 = resultat_np2(par_plus(2),:)
       temp3 = invs_real(par1,par2)*temp10 &
            + invs_real(par1,par3)*temp11
       temp1 = b_real(par1)*temp12
       temp1 = mult_div(-1._ki/2._ki,temp1)*4._ki
       temp2 = 0._ki
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if ( (j /= par3) .and. (j /= par2) ) then
                !
                truc1 = resultat22(j,par_plus(2),par_plus(3),:)
                temp2 = temp2 + invs_real(j,par1)*truc1
                !
             end if
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_rarg(1) = -temp3(3) + temp1(3) + temp2(3)
       res_3pC0i_rarg(2) = -temp3(4) + temp1(4) + temp2(4)
       !
    end if
    !
  end function a3pC0i_rarg
  !
  recursive function a3pC0i_carg(s1,s2,s3,m1,m2,m3,par1,par2,par3) result(res_3pC0i_carg)
    !
    implicit none
    real(ki), intent (in) :: s1,s2,s3
    complex(ki), intent (in) :: m1,m2,m3
    integer, intent (in) :: par1,par2,par3
    complex(ki) :: res_3pC0i_carg
    !
    integer :: j
    integer :: nb_par_loc
    integer, dimension(3) :: par_loc,par_plus
    complex(ki), dimension(2) :: truc1_c
    complex(ki) :: temp0
    complex(ki), dimension(2) :: temp1, temp2, temp3
    complex(ki), dimension(2) :: temp10, temp11, temp12
    complex(ki) :: ctemp
    integer :: ib,b_pro,b_pro_mj
    !
    b_pro = packb(s)
    !
    par_loc = (/par1,par2,par3/)
    par_plus = par_loc+1
    nb_par_loc = count(mask=par_loc/=0)
    !
    ! cas sans parametre de feynman au numerateur
    !
    if (nb_par_loc == 0) then
       ctemp = C0(s1,s2,s3,m1,m2,m3)
       res_3pC0i_carg = ctemp
       !
       ! cas avec un parametre de feynman au numerateur
       !
    else if (nb_par_loc == 1) then
       !
       if (deja_calcule_c(1)) then
          !
          temp0 = resultat_c(1,1)
          !
       else
          !
          temp0 = a3pC0i_carg(s1,s2,s3,m1,m2,m3,0,0,0)
          resultat_c(1,1) = temp0
          deja_calcule_c(1) = .true.
          !
       end if
       !
       temp1(:) = czero
       temp2(:) = czero
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if (deja_calcule2_c(j,1)) then
                !
                truc1_c = resultat2_c(j,1,:)
                !
             else
                !
                truc1_c = f2p(s_mat_p_loc,b_pro_mj) !returns complex array!
                resultat2_c(j,1,:) = truc1_c
                deja_calcule2_c(j,1) = .true.
                !
             end if
             !
             temp1 = temp1 + b_complex(j)*truc1_c
             temp2 = temp2 + invs_complex(j,par3)*truc1_c
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_carg = b_complex(par3)*(temp0 - temp1(2))/sumb_complex + temp2(2)
       !
       ! cas avec deux parametres de feynman au numerateur
       !
    else if (nb_par_loc == 2) then
       !
       if (deja_calcule_np2_c(par_plus(3))) then
          !
          temp11 = resultat_np2_c(par_plus(3),:)
          !
       else
          !
          temp11 = a3pC0i_np2_carg(s1,s2,s3,m1,m2,m3,0,0,par3)
          resultat_np2_c(par_plus(3),:) = temp11
          deja_calcule_np2_c(par_plus(3)) = .true.
          !
       end if
       !
       temp10 = resultat_np2_c(1,:)
       temp3 = invs_complex(par2,par3)*temp10
       temp1 = b_complex(par2)*temp11
       temp1 = mult_div(-2._ki/3._ki,temp1)*3._ki
       temp2(:) = czero
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if (j /= par3) then
                !
                if (deja_calcule2_c(j,par_plus(3))) then
                   !
                   truc1_c = resultat2_c(j,par_plus(3),:)
                   !
                else
                   !
                   truc1_c = f2p(s_mat_p_loc,b_pro_mj,par3) !returns complex array!
                   resultat2_c(j,par_plus(3),:) = truc1_c
                   deja_calcule2_c(j,par_plus(3)) = .true.
                   !
                end if
                !
                temp2 = temp2 + invs_complex(j,par2)*truc1_c
                !
             end if
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_carg = -temp3(2) + temp1(2) + temp2(2)
       !
       ! cas avec trois parametres de feynman au numerateur
       !
    else
       !
       temp12 = a3pC0i_np2_carg(s1,s2,s3,m1,m2,m3,0,par2,par3)
       temp10 = resultat_np2_c(par_plus(3),:)
       temp11 = resultat_np2_c(par_plus(2),:)
       temp3 = invs_complex(par1,par2)*temp10 &
            + invs_complex(par1,par3)*temp11
       temp1 = b_complex(par1)*temp12
       temp1 = mult_div(-1._ki/2._ki,temp1)*4._ki
       temp2(:) = czero
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if ( (j /= par3) .and. (j /= par2) ) then
                !
                truc1_c = resultat22_c(j,par_plus(2),par_plus(3),:) !!returns complex array!
                temp2 = temp2 + invs_complex(j,par1)*truc1_c
                !
             end if
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_carg = -temp3(2) + temp1(2) + temp2(2)
       !
    end if
    !
  end function a3pC0i_carg
  !
  !****if* src/integral/three_point/function_3pC0i/a3pC0i_np2
  ! NAME
  !
  !  Function a3pC0i_np2
  !
  ! USAGE
  !
  !  real_dim4 = a3pC0i_np2_rarg(s1,s2,s3,m1r,m2r,m3r,par1,par2,par3)
  !  complex_dim2 = a3pC0i_np2_carg(s1,s2,s3,m1c,m2c,m3c,par1,par2,par3)
  !
  ! DESCRIPTION
  !
  !  This recursive function implements the formula of ref 1
  !
  ! INPUTS
  !
  !  * s1 -- a real (type ki), p1^2
  !  * s2 -- a real (type ki), p2^2
  !  * s3 -- a real (type ki), p3^2
  !  * m1(r/c) -- a real/complex (type ki), the first internal mass squared
  !  * m2(r/c) -- a real/complex (type ki), the second internal mass squared
  !  * m3(r/c) -- a real/complex (type ki), the third internal mass squared
  !  * par1 -- an integer, the label of the third Feynman parameter
  !  * par2 -- an integer, the label of the second Feynman parameter
  !  * par3 -- an integer, the label of the first Feynman parameter
  !
  ! SIDE EFFECTS
  !
  !  This function modifies the value of the local (for the module) variables:
  !  * deja_calcule(_c), deja_calcule2, deja_calcule_np2(_c) and deja_calcule22
  !  * resultat(_c), resultat2, resultat_np2(_c) and resultat22
  !
  ! RETURN VALUE
  !
  !  It returns a real/complex (type ki) array of rank 1 and shape 4/2
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  !
  recursive function a3pC0i_np2_rarg(s1,s2,s3,m1,m2,m3,par1,par2,par3) result(res_3pC0i_np2_rarg)
    !
    implicit none

    real(ki), intent (in) :: s1,s2,s3,m1,m2,m3
    integer, intent (in) :: par1,par2,par3
    real(ki), dimension(4) :: res_3pC0i_np2_rarg
    !
    integer :: j
    integer :: nb_par_loc
    integer, dimension(3) :: par_loc,par_plus
    real(ki), dimension(4) :: truc1,truc2
    real(ki), dimension(2) :: temp0
    real(ki), dimension(4) :: temp1,temp2,temp3
    real(ki), dimension(4) :: temp10,temp11
    integer :: ib,b_pro,b_pro_mj
    !
    b_pro = packb(s)
    !
    par_loc = (/par1,par2,par3/)
    par_plus = par_loc+1
    nb_par_loc = count(mask=par_loc/=0)
    !
    ! cas sans parametre de feynman au numerateur
    !
    if (nb_par_loc == 0) then
       !
       if (deja_calcule(1)) then
          !
          temp0 = resultat(1,:)
          !
       else
          !
          temp0 = a3pC0i_rarg(s1,s2,s3,m1,m2,m3,0,0,0)
          resultat(1,:) = temp0
          deja_calcule(1) = .true.
          !
       end if
       !
       temp1 = 0._ki
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if (deja_calcule2(j,1)) then
                !
                truc1 = resultat2(j,1,:)
                !
             else
                !
                truc1 = f2p_ra(s_mat_p_loc,b_pro_mj) !returns real array!
                resultat2(j,1,:) = truc1
                deja_calcule2(j,1) = .true.
                !
             end if
             !
             temp1 = temp1 + b_real(j)*truc1
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_np2_rarg(1) = (- temp1(1))/sumb_real
       res_3pC0i_np2_rarg(2) = (- temp1(2))/sumb_real
       res_3pC0i_np2_rarg(3) = (temp0(1) - temp1(3))/sumb_real
       res_3pC0i_np2_rarg(4) = (temp0(2) - temp1(4))/sumb_real
       res_3pC0i_np2_rarg = mult_div(1._ki,res_3pC0i_np2_rarg)/2._ki
       !
       ! cas avec un parametre de feynman au numerateur
       !
    else if (nb_par_loc == 1) then
       !
       if (deja_calcule_np2(1)) then
          !
          temp10 = resultat_np2(1,:)
          !
       else
          !
          temp10 = a3pC0i_np2_rarg(s1,s2,s3,m1,m2,m3,0,0,0)
          resultat_np2(1,:) = temp10
          deja_calcule_np2(1) = .true.
          !
       end if
       !
       temp3 = b_real(par3)*temp10
       temp1 = 0._ki
       temp2 = 0._ki
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if (deja_calcule2(j,1)) then
                !
                truc1 = resultat2(j,1,:)
                !
             else
                !
                truc1 = f2p_ra(s_mat_p_loc,b_pro_mj) !returns real array!
                resultat2(j,1,:) = truc1
                deja_calcule2(j,1) = .true.
                !
             end if
             !
             temp1 = temp1 + invs_real(j,par3)*truc1
             !
             if (j /= par3) then
                !
                if (deja_calcule2(j,par_plus(3))) then
                   !
                   truc2 = resultat2(j,par_plus(3),:)
                   !
                else
                   !
                   truc2 = f2p_ra(s_mat_p_loc,b_pro_mj,par3) !returns real array!
                   resultat2(j,par_plus(3),:) = truc2
                   deja_calcule2(j,par_plus(3)) = .true.
                   !
                end if
                !
                temp2 = temp2 + b_real(j)*truc2
                !
             end if
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       temp1 = mult_div(2._ki/3._ki,temp1)/3._ki
       temp2 = mult_div(2._ki/3._ki,temp2)/3._ki
       res_3pC0i_np2_rarg(1) = (temp3(1) + temp1(1) - temp2(1))/sumb_real
       res_3pC0i_np2_rarg(2) = (temp3(2) + temp1(2) - temp2(2))/sumb_real
       res_3pC0i_np2_rarg(3) = (temp3(3) + temp1(3) - temp2(3))/sumb_real
       res_3pC0i_np2_rarg(4) = (temp3(4) + temp1(4) - temp2(4))/sumb_real
       !
       ! cas avec deux parametres de feynman au numerateur
       !
    else if (nb_par_loc == 2) then
       !
       temp0 = a3pC0i_rarg(s1,s2,s3,m1,m2,m3,0,par2,par3)
       !
       if (deja_calcule_np2(par_plus(2))) then
          !
          temp10 = resultat_np2(par_plus(2),:)
          !
       else
          !
          temp10 = a3pC0i_np2_rarg(s1,s2,s3,m1,m2,m3,0,0,par2)
          resultat_np2(par_plus(2),:) = temp10
          deja_calcule_np2(par_plus(2)) = .true.
          !
       end if
       !
       !
       if (deja_calcule_np2(par_plus(3))) then
          !
          temp11 = resultat_np2(par_plus(3),:)
          !
       else
          !
          temp11 = a3pC0i_np2_rarg(s1,s2,s3,m1,m2,m3,0,0,par3)
          resultat_np2(par_plus(3),:) = temp11
          deja_calcule_np2(par_plus(3)) = .true.
          !
       end if
       !
       temp3 =  b_real(par3)*temp10 + b_real(par2)*temp11
       temp1 = 0._ki
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if ( (j /= par2) .and. (j /= par3) ) then
                !
                if (deja_calcule22(j,par_plus(2),par_plus(3))) then
                   !
                   truc1 = resultat22(j,par_plus(2),par_plus(3),:)
                   !
                else
                   !
                   truc1 = f2p_ra(s_mat_p_loc,b_pro_mj,par2,par3) !returns real array!
                   resultat22(j,par_plus(2),par_plus(3),:) = truc1
                   deja_calcule22(j,par_plus(2),par_plus(3)) = .true.
                   !
                end if
                !
                temp1 = temp1 + b_real(j)*truc1
                !
             end if
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_np2_rarg(1) = (temp3(1) - temp1(1))/sumb_real
       res_3pC0i_np2_rarg(2) = (temp3(2) - temp1(2))/sumb_real
       res_3pC0i_np2_rarg(3) = (temp0(1) + temp3(3) - temp1(3))/sumb_real
       res_3pC0i_np2_rarg(4) = (temp0(2) + temp3(4) - temp1(4))/sumb_real
       res_3pC0i_np2_rarg = mult_div(1._ki/2._ki,res_3pC0i_np2_rarg)/4._ki
       !
       ! cas avec trois parametres de feynman au numerateur
       !
    else
       !
       tab_erreur_par(1)%a_imprimer = .true.
       tab_erreur_par(1)%chaine = 'Error in f3p_finite:'
       tab_erreur_par(2)%a_imprimer = .true.
       tab_erreur_par(2)%chaine = 'rank 3 6-dim 3-point function should not be needed'
       call catch_exception(0)
       !
    end if
    !
  end function a3pC0i_np2_rarg
  !
  recursive function a3pC0i_np2_carg(s1,s2,s3,m1,m2,m3,par1,par2,par3) result(res_3pC0i_np2_carg)
    !
    implicit none
    !
    real(ki), intent (in) :: s1,s2,s3
    complex(ki), intent (in) :: m1,m2,m3
    integer, intent (in) :: par1,par2,par3
    complex(ki), dimension(2) :: res_3pC0i_np2_carg
    !
    integer :: j
    integer :: nb_par_loc
    integer, dimension(3) :: par_loc,par_plus
    complex(ki), dimension(2) :: truc1_c, truc2_c
    complex(ki) :: temp0
    complex(ki), dimension(2) :: temp1, temp2, temp3
    complex(ki), dimension(2) :: temp10, temp11
    integer :: ib,b_pro,b_pro_mj
    !
    b_pro = packb(s)
    !
    par_loc = (/par1,par2,par3/)
    par_plus = par_loc+1
    nb_par_loc = count(mask=par_loc/=0)
    !
    ! cas sans parametre de feynman au numerateur
    !
    if (nb_par_loc == 0) then
       !
       if (deja_calcule_c(1)) then
          !
          temp0 = resultat_c(1,1)
          !
       else
          !
          temp0 = a3pC0i_carg(s1,s2,s3,m1,m2,m3,0,0,0)
          resultat_c(1,1) = temp0
          deja_calcule_c(1) = .true.
          !
       end if
       !
       temp1(:) = czero
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if (deja_calcule2_c(j,1)) then
                !
                truc1_c = resultat2_c(j,1,:)
                !
             else
                !
                truc1_c = f2p(s_mat_p_loc,b_pro_mj) !!! returns complex array!
                resultat2_c(j,1,:) = truc1_c
                deja_calcule2_c(j,1) = .true.
                !
             end if
             !
             temp1 = temp1 + b_complex(j)*truc1_c
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_np2_carg(1) = (- temp1(1))/sumb_complex
       res_3pC0i_np2_carg(2) = (temp0 - temp1(2))/sumb_complex
       res_3pC0i_np2_carg = mult_div(1._ki,res_3pC0i_np2_carg)/2._ki
       !
       ! cas avec un parametre de feynman au numerateur
       !
    else if (nb_par_loc == 1) then
       !
       if (deja_calcule_np2_c(1)) then
          !
          temp10 = resultat_np2_c(1,:)
          !
       else
          !
          temp10 = a3pC0i_np2_carg(s1,s2,s3,m1,m2,m3,0,0,0)
          resultat_np2_c(1,:) = temp10
          deja_calcule_np2_c(1) = .true.
          !
       end if
       !
       temp3 = b_complex(par3)*temp10
       temp1(:) = czero
       temp2(:) = czero
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if (deja_calcule2_c(j,1)) then
                !
                truc1_c = resultat2_c(j,1,:)
                !
             else
                !
                truc1_c = f2p(s_mat_p_loc,b_pro_mj) !!! returns complex array!
                resultat2_c(j,1,:) = truc1_c
                deja_calcule2_c(j,1) = .true.
                !
             end if
             !
             temp1 = temp1 + invs_complex(j,par3)*truc1_c
             !
             if (j /= par3) then
                !
                if (deja_calcule2_c(j,par_plus(3))) then
                   !
                   truc2_c = resultat2_c(j,par_plus(3),:)
                   !
                else
                   !
                   truc2_c = f2p(s_mat_p_loc,b_pro_mj,par3) !!! returns complex array!
                   resultat2_c(j,par_plus(3),:) = truc2_c
                   deja_calcule2_c(j,par_plus(3)) = .true.
                   !
                end if
                !
                temp2 = temp2 + b_complex(j)*truc2_c
                !
             end if
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       temp1 = mult_div(2._ki/3._ki,temp1)/3._ki
       temp2 = mult_div(2._ki/3._ki,temp2)/3._ki
       res_3pC0i_np2_carg(1) = (temp3(1) + temp1(1) - temp2(1))/sumb_complex
       res_3pC0i_np2_carg(2) = (temp3(2) + temp1(2) - temp2(2))/sumb_complex
       !
       ! cas avec deux parametres de feynman au numerateur
       !
    else if (nb_par_loc == 2) then
       !
       temp0 = a3pC0i_carg(s1,s2,s3,m1,m2,m3,0,par2,par3)
       !
       if (deja_calcule_np2_c(par_plus(2))) then
          !
          temp10 = resultat_np2_c(par_plus(2),:)
          !
       else
          !
          temp10 = a3pC0i_np2_carg(s1,s2,s3,m1,m2,m3,0,0,par2)
          resultat_np2_c(par_plus(2),:) = temp10
          deja_calcule_np2_c(par_plus(2)) = .true.
          !
       end if
       !
       if (deja_calcule_np2_c(par_plus(3))) then
          !
          temp11 = resultat_np2_c(par_plus(3),:)
          !
       else
          !
          temp11 = a3pC0i_np2_carg(s1,s2,s3,m1,m2,m3,0,0,par3)
          resultat_np2_c(par_plus(3),:) = temp11
          deja_calcule_np2_c(par_plus(3)) = .true.
          !
       end if
       !
       temp3 =  b_complex(par3)*temp10 + b_complex(par2)*temp11
       temp1(:) = czero
       !
       ib = b_pro
       j = 0
       !
       do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
             !
             b_pro_mj = ibclr(b_pro,j)
             !
             if ( (j /= par2) .and. (j /= par3) ) then
                !
                if (deja_calcule22_c(j,par_plus(2),par_plus(3))) then
                   !
                   truc1_c = resultat22_c(j,par_plus(2),par_plus(3),:)
                   !
                else
                   !
                   truc1_c = f2p(s_mat_p_loc,b_pro_mj,par2,par3) !!!returns complex array
                   resultat22_c(j,par_plus(2),par_plus(3),:) = truc1_c
                   deja_calcule22_c(j,par_plus(2),par_plus(3)) = .true.
                   !
                end if
                !
                temp1 = temp1 + b_complex(j)*truc1_c
                !
             end if
             !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
       end do
       !
       res_3pC0i_np2_carg(1) = (temp3(1) - temp1(1))/sumb_complex
       res_3pC0i_np2_carg(2) = (temp0 + temp3(2) - temp1(2))/sumb_complex
       res_3pC0i_np2_carg = mult_div(1._ki/2._ki,res_3pC0i_np2_carg)/4._ki
       !
       ! cas avec trois parametres de feynman au numerateur
       !
    else
       !
       tab_erreur_par(1)%a_imprimer = .true.
       tab_erreur_par(1)%chaine = 'Error in f3p_finite:'
       tab_erreur_par(2)%a_imprimer = .true.
       tab_erreur_par(2)%chaine = 'rank 3 6-dim 3-point function should not be needed'
       call catch_exception(0)
       !
    end if
    !
  end function a3pC0i_np2_carg
  !
  !****f* src/integral/three_point/function_3pC0i/C0
  ! NAME
  !
  ! Function C0
  !
  ! USAGE
  !
  !  complex = C0(s1,s2,s3,m1,m2,m3)
  !
  ! DESCRIPTION
  !
  !  This function computes finite scalar three point functions
  !  with internal masses in 4 dimensions
  !
  ! INPUTS
  !
  !  * s1 -- a real/complex (type ki), p1^2
  !  * s2 -- a real/complex (type ki), p2^2
  !  * s3 -- a real/complex (type ki), p3^2
  !  * m1 -- a real/complex (type ki), the first internal mass squared
  !  * m2 -- a real/complex (type ki), the second internal mass squared
  !  * m3 -- a real/complex (type ki), the third internal mass squared
  !
  ! SIDE EFFECTS
  !
  !  No side effect, it uses the value of rat_or_tot_par 
  !  (in src/module/parametre.f90)
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki)
  !
  !  
  !*****
  function C0_rarg(s1,s2,s3,m1sq,m2sq,m3sq)
    implicit none
    !
    !  s1,s2,s3 = SQUARED external momenta
    !  m1sq,m2sq,m3sq = SQUARED internal masses
    real(ki), intent(in) :: s1,s2,s3,m1sq,m2sq,m3sq
    real(ki) :: s1r,s2r,s3r
    ! real(ki) :: del
    complex(ki) :: C0_rarg
    complex(ki) :: c0_golem
    !  del = epsilon(1._ki)
    !
    s1r = s1
    s2r = s2
    s3r = s3
    !
!    if (equal_real(s1r,zero) ) s1r = 0._ki
!    if (equal_real(s2r,zero) ) s2r = 0._ki
!    if (equal_real(s3r,zero) ) s3r = 0._ki   
    !
    if (rat_or_tot_par%tot_selected) then
       !
       c0_golem =  i3p3m_3mi()
       !
       C0_rarg=c0_golem
       !
    else !if (rat_or_tot_par%rat_selected) then
       !
       C0_rarg = czero
       !
    end if
    !
  end function C0_rarg
  !
  function C0_carg(s1,s2,s3,m1sq,m2sq,m3sq)
    implicit none
    !
    !  s1,s2,s3 = SQUARED external momenta
    !  m1sq,m2sq,m3sq = SQUARED internal complex masses
    real(ki), intent(in) :: s1,s2,s3
    complex(ki), intent(in) :: m1sq,m2sq,m3sq
    complex(ki_avh) :: cp1,cp2,cp3,cm1,cm2,cm3
    real(ki) :: s1r,s2r,s3r
    ! real(ki) :: del
    complex(ki) :: C0_carg
    complex(ki) :: c0_golem
    !integer :: i,j,k
    !  del = epsilon(1._ki)
    !
    s1r = s1
    s2r = s2
    s3r = s3
    !
!    if (equal_real(s1r,zero) ) s1r = 0._ki
!    if (equal_real(s2r,zero) ) s2r = 0._ki
!    if (equal_real(s3r,zero) ) s3r = 0._ki
    !    
    if (rat_or_tot_par%tot_selected) then
       !
       !
       cp1 = cmplx(s1r,0._ki_avh,kind=ki_avh)
       cp2 = cmplx(s2r,0._ki_avh,kind=ki_avh)
       cp3 = cmplx(s3r,0._ki_avh,kind=ki_avh)
       !
       cm1 = cmplx(m1sq,kind=ki_avh)
       cm2 = cmplx(m2sq,kind=ki_avh)
       cm3 = cmplx(m3sq,kind=ki_avh)
       !
       !
       ! we use now the golem implementation of C0 (jpg)
       c0_golem =  i3p3m_3mi()
       !
       C0_carg=c0_golem
       !
       !
    else !if (rat_or_tot_par%rat_selected) then
       !
       C0_carg = czero
       !
    end if
    !
  end function C0_carg
  !
  ! *************************************************
   !added by jpg l2 30/04/2011
   subroutine compute_deter(s13,s12,s23,s11,s22,s33)
      !
      complex(ki), intent(in) :: s13,s12,s23,s11,s22,s33
      !
      complex(ki) :: temp
      !
      det_s = s11*(s22*s33-s23**2)-s12*(s12*s33-s13*s23)+s13*(s12*s23-s13*s22)
      !
      num_b1 = s22*s33-s12*s33-s23**2+s13*s23+s12*s23-s13*s22
      num_b2 = -s12*s33+s11*s33+s13*s23-s11*s23-s13**2+s12*s13
      num_b3 = s12*s23-s11*s23-s13*s22+s11*s22+s12*s13-s12**2
      num_sumb = s22*s33-2*s12*s33+s11*s33-s23**2+2._ki*s13*s23+2._ki*s12*s23-2._ki&
                          *s11*s23-2._ki*s13*s22+s11*s22-s13**2+2._ki*s12*s13-s12**2
      !
      temp = (-s33+2._ki*s13-s11)*(-s33+2._ki*s23-s22)-(-s33+s23+s13-s12)**2
      det_g = real(temp,ki) ! det(G) is real by construction
      !write(*,*) 'test compute_deter :',s13,s12,s23,s11,s22,s33
      !
   end subroutine compute_deter
   !
   function i3p3m_3mi()
      !
      complex(ki) :: i3p3m_3mi
      real(ki), dimension(3) :: s
      real(ki) :: s_alpha
      !
     !write(6,*) "MR:  in  i3p3m_3mi oneone ",g_glob(1,1)," onetwo ",g_glob(1,2), " onethree ", g_glob(1,3),&
!&" twoone ",g_glob(2,1)," twotwo ",g_glob(2,2), " twothree ", g_glob(2,3) ,&
!&" threeone ",g_glob(3,1)," threetwo ",g_glob(3,2), " threethree ", g_glob(3,3) 
     !s_alpha = -sign(un,real(g_glob(2,1)-g_glob(3,1),ki))
      !s_alpha = sign(un,real(g_glob(2,1)-g_glob(3,1),ki))
      !if (sign(un,real(num_b1)) == s_alpha) then
        !s(1) = 1._ki
      !else
        !s(1) = -1._ki
      !end if
      !if (sign(un,real(num_b2,ki)) == s(1)) then
        !s(2) = s(1)
      !else
        !s(2) = -s(1)
      !end if
      !if (sign(un,real(num_b3,ki)) == s(1)) then
        !s(3) = -s(1)
      !else
        !s(3) = s(1)
      !end if
      !write(*,*) 'test i3p3m_3mi :',tab_test_g2_glob
      if (tab_test_g2_glob(1) ) then
        s(1) = -sign(un,real(g_glob(1,2),ki)/real(num_b1,ki))
        s(2) = s(1)
        s(3) = -s(1)
      else if (tab_test_g2_glob(2) ) then
        s(2) = -sign(un,real(g_glob(2,2),ki)/real(num_b2,ki))
        s(1) = s(2)
        s(3) = -s(2)
      else if (tab_test_g2_glob(3) ) then
      !write(*,*) 'test i3p3m_3mi prime :',g_glob(3,2),num_b3
        s(3) = -sign(un,real(g_glob(3,2),ki)/real(num_b3,ki))
        s(1) = -s(3)
        s(2) = -s(3)
      else
        s(1) = 1._ki
        s(2) = 1._ki
        s(3) = -1._ki
      end if
      !write(*,*) 'differents s :',s(1),s(2),s(3)
      !i3p3m_3mi = -(j3p3m_3mi(real(g_glob(1,1),ki),g_glob(1,2),g_glob(1,3),1)&
                          !&+j3p3m_3mi(real(g_glob(2,1),ki),g_glob(2,2),g_glob(2,3),2)&
                          !&+j3p3m_3mi(real(g_glob(3,1),ki),g_glob(3,2),g_glob(3,3),3) )
      i3p3m_3mi = -(j3p3m_3mi(real(g_glob(1,1),ki),g_glob(1,2),g_glob(1,3),1,s(1))&
                          &+j3p3m_3mi(real(g_glob(2,1),ki),g_glob(2,2),g_glob(2,3),2,s(2))&
                          &+j3p3m_3mi(real(g_glob(3,1),ki),g_glob(3,2),g_glob(3,3),3,s(3)) )
      !
    !write(6,*) "MR:  i3p3m_3mi donee ",i3p3m_3mi 


   end function i3p3m_3mi
   !
   !function j3p3m_3mi(g2,g1,g0,i)
   function j3p3m_3mi(g2,g1,g0,i,s)
    !
     real(ki), intent(in) :: g2
     complex(ki), intent(in) :: g1,g0
     integer, intent(in) :: i
     real(ki), intent(in) :: s
     complex(ki) :: j3p3m_3mi
     !
     complex(ki) :: y0,num_b,coeff
     !
       !write(6,*) "Enter j3p3m_3mi(g2 ",g2," g1 ",g1," g0 ",g0," i ",i
       !write(6,*) "with det_g  ",det_g,"  det_s  ",det_s

    select case(i)
     !
     case(1)
       !
       num_b = num_b1
       !
     case(2)
       !
       num_b = num_b2
       !
     case(3)
       !
       num_b = num_b3
       !
     end select
     !
     ! we can have the case where g2=g1=0, in this case, the b_i corresponding is zero
     !
     if (tab_test_glob(i)) then
       !
       j3p3m_3mi = 0._ki
       !
     else
       !
       if (det_g < 0._ki) then
         !
         ! case where g2=0, this case appear only for det(g) < 0
         !
         if (abs(g2) < 1.e-10_ki) then
           !
         !write(6,*) "y0 1 num_b",num_b," s ",s
         y0 = 2._ki*(g0+det_s/2._ki/det_g)/(-g1+s*num_b/sqrt(-det_g))
          !write(6,*) "y0 1 ",y0
           !
         else
           !
         !write(6,*) "y0 2 "
           y0 = -g1/(2._ki*g2)-s*num_b/(2._ki*g2*sqrt(-det_g))
            !write(6,*) "y0 2 ",y0
         !
         end if
         !
         coeff = 1._ki/(2._ki*sqrt(-det_g))
         !
       else if (det_g > 0._ki) then
         !
         !write(6,*) "y0 2 "
         y0 = -g1/(2._ki*g2)-i_*s*num_b/(2._ki*abs(g2)*sqrt(det_g))
           !write(6,*) "y0 2 ",y0
        coeff = i_*sign(un,g2)/(2._ki*sqrt(det_g))
         !
       end if
       !
     !write(6,*) "Call j3p3m_3mi coeff ",coeff," g2 ",g2," g1 ",g1," g0 ",g0," y0 ",y0, " last ", -det_s/det_g/2._ki-i_*delta_glob
       j3p3m_3mi = 2._ki*s*coeff*int_log(g2,g1,g0,y0,-det_s/det_g/2._ki-i_*delta_glob)
       !
     end if
     !
     !write(6,*) "Exit j3p3m_3mi(g2 ",g2," g1 ",g1," g0 ",g0," i ",i, " with ", j3p3m_3mi 
      !
   end function j3p3m_3mi
   !
   ! compute the integral
   ! \int^1_0 ds (ln(a*s^2 + b*s + c - i*\lambda)-ln(arg_c))/(s-z_1)
   ! with arg_c = (a*z_1^2 + b*z_1 + c - i*\lambda)
   ! analytically. In the case of complex masses, b and c become complex
   ! but Im(b)*s+Im(c) are of the type (l1-l2)*s+l2 where l1,l2 are real 
   ! and have the same sign, so (l1-l2)*s+l2 has always the same sign when
   ! s runs between 0 and 1 : the sign of l1 (or l2)
   !
   !
   function int_log(a,b,c,z1,arg_c)
      !
      real(ki), intent(in) :: a
      complex(ki), intent(in) :: b,c
      complex(ki), intent(in) :: z1,arg_c
      complex(ki) :: int_log
      !
      complex(ki) :: s_plus,s_moins
      complex(ki) :: arg_sqrt
      real(ki) :: small
      real(ki) :: epsilon_gv,delta_gv
      complex(ki) :: rest,eta_part,delta
      complex(ki) :: eta_temp
      real(ki) :: z1_r,s_r_plus,s_r_moins,b_r,c_r
      real(ki) :: grand
      real(ki) :: eps_b
      !
      grand = 1.e+13_ki
      arg_sqrt = b*b - 4._ki*a*c
      ! on suppose qu'il y a une petite patie imaginaire
      ! negative : c --> c - i lambda
      small = delta_glob
      !delta_glob = 1.e-12_ki
      !small = delta_glob*1.e+0_ki
      !write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      !write(*,*) 'test int_log :',a,b,c,z1,arg_c
      !faire le cas ou a =0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      if (equal_real(aimag(b),zero) .and. equal_real(aimag(c),zero) ) then
        !
        b_r = real(b,ki)
        c_r = real(c,ki)
        !
        if ( equal_real(aimag(z1),zero) ) then ! z1 is real, no eta function in this case
          !
          z1_r = real(z1,ki)
          !
          if (real(arg_sqrt,ki) <= zero) then ! negative discriminant
            !
            delta = i_*sqrt(abs(arg_sqrt))*sign(un,a)
            !
            s_plus= (-b + delta)/(2._ki*a)
            s_moins = (-b - delta)/(2._ki*a)
            !
            int_log = cdilog(z1/(z1-s_plus)) - cdilog((z1-1._ki)/(z1-s_plus)) &
            &+ cdilog(z1/(z1-s_moins)) - cdilog((z1-1._ki)/(z1-s_moins))
            !
          else ! positive discriminant
            !
            ! careful treatement of the case where a --> 0, one root is -c/b,
            ! the other goes to infinity but is frozen if greater than grand
            ! The fact that one root is infinite does not disturb the
            ! computation of dilog for this root appear at the denominator.
            ! Obiouvsly, if a-->0, the discriminant is positif, so we always
            ! ended to this case
            delta = sqrt(arg_sqrt)
            !
            if (b_r >= zero) then
              !
              s_plus = -2._ki*c/(b+delta)
              !
              if (abs(a) <= b_r/(grand+sign(un,a)*c_r/b_r) ) then
                !
                s_moins = -sign(un,a)*grand
                !
              else
                !
                s_moins = (-b - delta)/(2._ki*a)
                !
              end if
              !
            else if (b_r < zero) then
              !
              s_moins = -2._ki*c/(b-delta)
              !
              if (abs(a) <= -b_r/(grand-sign(un,a)*c_r/b_r) ) then
                !
                s_plus = sign(un,a)*grand
                !
              else
                !
                s_plus = (-b + delta)/(2._ki*a)
                !
              end if
              !
            end if
            !s_plus= (-b + delta)/(2._ki*a)
            !s_moins = (-b - delta)/(2._ki*a)
            s_r_plus = real(s_plus,ki)
            s_r_moins = real(s_moins,ki)
            !
            !int_log = zdilog(z1/(z1-s_plus),1._ki) - zdilog((z1-1._ki)/(z1-s_plus),1._ki) &
            !&+ zdilog(z1/(z1-s_moins),-1._ki) - zdilog((z1-1._ki)/(z1-s_moins),-1._ki) 
            int_log = zdilog(z1_r/(z1_r-s_r_plus),sign(un,z1_r)) - zdilog((z1_r-1._ki)/(z1_r-s_r_plus),sign(un,z1_r-1._ki)) &
            &+ zdilog(z1_r/(z1_r-s_r_moins),-sign(un,z1_r)) - zdilog((z1_r-1._ki)/(z1_r-s_r_moins),-sign(un,z1_r-1._ki)) 
            !
          end if
          !
        else ! z1 is complex, s_plus and s_moins are still complex conjugate
          !
          epsilon_gv = un
          delta_gv = -sign(un,aimag(arg_c))
          !
          if (real(arg_sqrt,ki) <= zero) then ! negative discriminant
            !
            delta = i_*sqrt(abs(arg_sqrt))*sign(un,a)
            !
            ! be careful that for the case where z1 has a real part which is -b/2/a
            ! we have to keep a small real part
            s_plus = (-b + delta)/(2._ki*a) + small
            s_moins = (-b - delta)/(2._ki*a) - small
            !s_plus = (-b + delta)/(2._ki*a)
            !s_moins = (-b - delta)/(2._ki*a)
            !
            ! we treat explicitely the case where Re(z1) = -b/2/a
            !
            if (equal_real(real(z1,ki),real(-b/2._ki/a,ki))) then
              if ( (sign(un,aimag(z1-s_plus)) > 0._ki) .and. (sign(un,aimag(z1-s_moins)) > 0._ki) ) then
                eta_temp = -2._ki*i_*pi
              else
                eta_temp = 0._ki
              end if
            else
              eta_temp = eta(z1-s_plus,z1-s_moins)
            end if
            !
            ! s_plus and s_moins are complex conjugate, so eta(-s_plus,-s_moins) = 0
            eta_part = -eta_temp &
            & -eta(a-i_*small*epsilon_gv,1._ki/(a-i_*small*delta_gv))
            !
            rest = cdilog(z1/(z1-s_plus)) - cdilog((z1-1._ki)/(z1-s_plus)) &
            &+ eta(-s_plus,1._ki/(z1-s_plus))*log(z1/(z1-s_plus)) &
            &- eta(1._ki-s_plus,1._ki/(z1-s_plus))*log((z1-1._ki)/(z1-s_plus)) &
            &+ cdilog(z1/(z1-s_moins)) - cdilog((z1-1._ki)/(z1-s_moins)) &
            &+ eta(-s_moins,1._ki/(z1-s_moins))*log(z1/(z1-s_moins)) &
            &- eta(1._ki-s_moins,1._ki/(z1-s_moins))*log((z1-1._ki)/(z1-s_moins)) 
            !
          else ! positive discriminant
            !
            delta = sqrt(arg_sqrt)
            ! we need the small imaginary part for z1 is complex when
            ! Re(z1) = -b/2/a
            if (b_r >= zero) then
              s_plus = -2._ki*c/(b+delta)
              if (abs(a) <= b_r/(grand+sign(un,a)*c_r/b_r) ) then
                s_moins = -sign(un,a)*grand
              else
                s_moins = (-b - delta)/(2._ki*a)
              end if
            else if (b_r < zero) then
              s_moins = -2._ki*c/(b-delta)
              if (abs(a) <= -b_r/(grand-sign(un,a)*c_r/b_r) ) then
                s_plus = sign(un,a)*grand
              else
                s_plus = (-b + delta)/(2._ki*a)
              end if
            end if
            s_plus = s_plus + i_*small
            s_moins = s_moins - i_*small
            !s_plus= (-b + delta)/(2._ki*a)+i_*small
            !s_moins = (-b - delta)/(2._ki*a)-i_*small
            !
            ! we treat explicitely the case where Re(z1) = -b/2/a
            !
            !write(*,*) 'test equal :',equal_real(real(z1,ki),real(-b/2._ki/a,ki))
            if (equal_real(real(z1,ki),real(-b/2._ki/a,ki))) then
              if (sign(un,a) == sign(un,aimag(z1))) then
                eta_temp = -2._ki*i_*pi*sign(un,a)
              else
                eta_temp = 0._ki
              end if
            else
              eta_temp = eta(z1-s_plus,z1-s_moins)
            end if
            !
            ! s_plus and s_moins are complex conjugate, so eta(-s_plus,-s_moins) = 0
            eta_part = -eta_temp &
            & -eta(a-i_*small*epsilon_gv,1._ki/(a-i_*small*delta_gv))
            !
            rest = cdilog(z1/(z1-s_plus)) - cdilog((z1-1._ki)/(z1-s_plus)) &
            &+ eta(-s_plus,1._ki/(z1-s_plus))*log(z1/(z1-s_plus)) &
            &- eta(1._ki-s_plus,1._ki/(z1-s_plus))*log((z1-1._ki)/(z1-s_plus)) &
            &+ cdilog(z1/(z1-s_moins)) - cdilog((z1-1._ki)/(z1-s_moins)) &
            &+ eta(-s_moins,1._ki/(z1-s_moins))*log(z1/(z1-s_moins)) &
            &- eta(1._ki-s_moins,1._ki/(z1-s_moins))*log((z1-1._ki)/(z1-s_moins)) 
            !
          end if
          !
          eta_part = eta_part*log((z1-1._ki)/z1)
          !
          int_log = rest + eta_part
           !
        end if
          !
      else ! b,c and z1 are complexs
        !
        delta = sqrt(arg_sqrt)
        eps_b = real(exp(eta(b,b)/2._ki),ki)
        delta = delta*eps_b
        ! treatement of the case where a --> 0
        if ( equal_real(eps_b,un) ) then
          s_plus = -2._ki*c/(b+delta)
          !if (abs(-b/a+c/b) >= grand) then
          if (abs(a) <= abs(-b)/grand) then
            !s_moins = grand*exp(i_*atan2(aimag(-b/a+c/b),real(-b/a+c/b,ki)))
            s_moins = -grand*sign(un,a)*b
          else
            s_moins = (-b - delta)/(2._ki*a)
          end if
        else if ( equal_real(eps_b,-un) ) then
          s_moins = -2._ki*c/(b-delta)
          !if (abs(-b/a+c/b) >= grand) then
          if (abs(a) <= abs(-b)/grand) then
            !s_plus = grand*exp(i_*atan2(real(-b/a+c/b,ki),aimag(-b/a+c/b)))
            s_plus = -grand*sign(un,a)*b
          else
            s_plus = (-b + delta)/(2._ki*a)
          end if
        end if
        !s_plus= (-b + delta)/(2._ki*a)
        !s_moins = (-b - delta)/(2._ki*a)
        !
        epsilon_gv = un
        delta_gv = -sign(un,aimag(arg_c))
        !
        eta_part = eta(-s_plus,-s_moins)-eta(z1-s_plus,z1-s_moins) &
        & -eta(a-i_*small*epsilon_gv,1._ki/(a-i_*small*delta_gv))
        !
        eta_part = eta_part*log((z1-1._ki)/z1)
        !
        rest = cdilog(z1/(z1-s_plus)) - cdilog((z1-1._ki)/(z1-s_plus)) &
        &+ eta(-s_plus,1._ki/(z1-s_plus))*log(z1/(z1-s_plus)) &
        &- eta(1._ki-s_plus,1._ki/(z1-s_plus))*log((z1-1._ki)/(z1-s_plus)) &
        &+ cdilog(z1/(z1-s_moins)) - cdilog((z1-1._ki)/(z1-s_moins)) &
        &+ eta(-s_moins,1._ki/(z1-s_moins))*log(z1/(z1-s_moins)) &
        &- eta(1._ki-s_moins,1._ki/(z1-s_moins))*log((z1-1._ki)/(z1-s_moins)) 
        !
        int_log = rest + eta_part
        !
      end if
      !
    end function int_log
    !
   !!
   !! compute the integral
   !! \int^1_0 ds (ln(a*s^2 + b*s + c - i*\lambda)-ln(arg_c))/(s-z_1)
   !! with arg_c = (a*z_1^2 + b*z_1 + c - i*\lambda)
   !! analytically. In the case of complex masses, b and c become complex
   !! but Im(b)*s+Im(c) are of the type (l1-l2)*s+l2 where l1,l2 are real 
   !! and have the same sign, so (l1-l2)*s+l2 has always the same sign when
   !! s runs between 0 and 1 : the sign of l1 (or l2)
   !!
   !!
   !function int_log(a,b,c,z1,arg_c)
      !!
      !real(ki), intent(in) :: a
      !complex(ki), intent(in) :: b,c
      !complex(ki), intent(in) :: z1,arg_c
      !complex(ki) :: int_log
      !!
      !complex(ki) :: s_plus,s_moins
      !complex(ki) :: arg_sqrt
      !real(ki) :: small
      !real(ki) :: epsilon_gv,delta_gv
      !complex(ki) :: rest,eta_part,delta
      !complex(ki) :: eta_temp
      !!
      !arg_sqrt = b*b - 4._ki*a*c
      !! on suppose qu'il y a une petite patie imaginaire
      !! negative : c --> c - i lambda
      !!small = delta_glob
      !small = delta_glob*1.e+0_ki
      !write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      !write(*,*) 'test int_log :',a,b,c,z1,arg_c
      !!faire le cas ou a =0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !if (equal_real(aimag(b),zero) .and. equal_real(aimag(c),zero) ) then
        !!
        !if ( equal_real(aimag(z1),zero) ) then ! z1 is real, no eta function in this case
          !!
          !if (real(arg_sqrt,ki) <= zero) then
            !!
            !delta = i_*sqrt(abs(arg_sqrt))*sign(un,a)
            !!
            !s_plus= (-b + delta)/(2._ki*a)
            !s_moins = (-b - delta)/(2._ki*a)
            !!
            !int_log = cdilog(z1/(z1-s_plus)) - cdilog((z1-1._ki)/(z1-s_plus)) &
            !&+ cdilog(z1/(z1-s_moins)) - cdilog((z1-1._ki)/(z1-s_moins))
            !!write(*,*) 'tist 1:',cdilog(z1/(z1-s_plus)),- cdilog((z1-1._ki)/(z1-s_plus))
            !!write(*,*) 'tist 1:',cdilog(z1/(z1-s_moins)), - cdilog((z1-1._ki)/(z1-s_moins))
            !!
          !else
            !!
            !delta = sqrt(arg_sqrt)
            !s_plus= (-b + delta)/(2._ki*a)
            !s_moins = (-b - delta)/(2._ki*a)
            !!
            !!int_log = zdilog(z1/(z1-s_plus),1._ki) - zdilog((z1-1._ki)/(z1-s_plus),1._ki) &
            !!&+ zdilog(z1/(z1-s_moins),-1._ki) - zdilog((z1-1._ki)/(z1-s_moins),-1._ki) 
            !int_log = zdilog(z1/(z1-s_plus),sign(un,real(z1,ki))) - zdilog((z1-1._ki)/(z1-s_plus),sign(un,real(z1-1._ki))) &
            !&+ zdilog(z1/(z1-s_moins),-sign(un,real(z1,ki))) - zdilog((z1-1._ki)/(z1-s_moins),-sign(un,real(z1-1._ki,ki))) 
            !!write(*,*) 'tist 2:',s_plus,s_moins
            !!write(*,*) 'tist 2:',int_log
            !!
          !end if
          !!
        !else ! z1 is complex, s_plus and s_moins are still complex conjugate
          !!
          !epsilon_gv = un
          !delta_gv = -sign(un,aimag(arg_c))
          !!
          !if (real(arg_sqrt,ki) <= zero) then
            !!
            !delta = i_*sqrt(abs(arg_sqrt))*sign(un,a)
            !!
            !! be careful that for the case where z1 has a real part which is -b/2/a
            !! we have to keep a small real part
            !s_plus = (-b + delta)/(2._ki*a) + small
            !s_moins = (-b - delta)/(2._ki*a) - small
            !!
            !! we treat explicitely the case where Re(z1) = -b/2/a
            !!
            !if (equal_real(real(z1,ki),real(-b/2._ki/a,ki))) then
              !if ( (sign(un,aimag(z1-s_plus)) > 0._ki) .and. (sign(un,aimag(z1-s_moins)) > 0._ki) ) then
                !eta_temp = -2._ki*i_*pi
              !else
                !eta_temp = 0._ki
              !end if
            !else
              !eta_temp = eta(z1-s_plus,z1-s_moins)
            !end if
            !!
            !! s_plus and s_moins are complex conjugate, so eta(-s_plus,-s_moins) = 0
            !eta_part = -eta_temp &
            !& -eta(a-i_*small*epsilon_gv,1._ki/(a-i_*small*delta_gv))
            !!
      !!write(*,*) 'test int_log 0 :',z1-s_plus,z1-s_moins
      !!write(*,*) 'test int_log 1 :',real(z1,ki),real(-b/2._ki/a,ki),eta(z1-s_plus,z1-s_moins)
            !rest = cdilog(z1/(z1-s_plus)) - cdilog((z1-1._ki)/(z1-s_plus)) &
            !&+ eta(-s_plus,1._ki/(z1-s_plus))*log(z1/(z1-s_plus)) &
            !&- eta(1._ki-s_plus,1._ki/(z1-s_plus))*log((z1-1._ki)/(z1-s_plus)) &
            !&+ cdilog(z1/(z1-s_moins)) - cdilog((z1-1._ki)/(z1-s_moins)) &
            !&+ eta(-s_moins,1._ki/(z1-s_moins))*log(z1/(z1-s_moins)) &
            !&- eta(1._ki-s_moins,1._ki/(z1-s_moins))*log((z1-1._ki)/(z1-s_moins)) 
      !!write(*,*) 'test int_log res :',eta_part,rest
            !!
          !else
            !!
            !delta = sqrt(arg_sqrt)
            !! we need the small imaginary part for z1 is complex when
            !! Re(z1) = -b/2/a
            !s_plus= (-b + delta)/(2._ki*a)+i_*small
            !s_moins = (-b - delta)/(2._ki*a)-i_*small
            !!
            !! we treat explicitely the case where Re(z1) = -b/2/a
            !!
            !if (equal_real(real(z1,ki),real(-b/2._ki/a,ki))) then
              !if (sign(un,a) == sign(un,aimag(z1))) then
                !eta_temp = -2._ki*i_*pi*sign(un,a)
              !else
                !eta_temp = 0._ki
              !end if
            !else
              !eta_temp = eta(z1-s_plus,z1-s_moins)
            !end if
            !!
            !! s_plus and s_moins are complex conjugate, so eta(-s_plus,-s_moins) = 0
            !eta_part = -eta_temp &
            !& -eta(a-i_*small*epsilon_gv,1._ki/(a-i_*small*delta_gv))
      !!write(*,*) 'test int_log 0 :',z1-s_plus,z1-s_moins
      !!write(*,*) 'test int_log 1 :',real(z1,ki),real(-b/2._ki/a,ki),eta_temp
            !rest = cdilog(z1/(z1-s_plus)) - cdilog((z1-1._ki)/(z1-s_plus)) &
            !&+ eta(-s_plus,1._ki/(z1-s_plus))*log(z1/(z1-s_plus)) &
            !&- eta(1._ki-s_plus,1._ki/(z1-s_plus))*log((z1-1._ki)/(z1-s_plus)) &
            !&+ cdilog(z1/(z1-s_moins)) - cdilog((z1-1._ki)/(z1-s_moins)) &
            !&+ eta(-s_moins,1._ki/(z1-s_moins))*log(z1/(z1-s_moins)) &
            !&- eta(1._ki-s_moins,1._ki/(z1-s_moins))*log((z1-1._ki)/(z1-s_moins)) 
      !!write(*,*) 'test int_log res :',eta_part,rest
            !!
          !end if
          !!
          !eta_part = eta_part*log((z1-1._ki)/z1)
          !!
          !int_log = rest + eta_part
           !!
        !end if
          !!
      !else
        !!
        !delta = sqrt(arg_sqrt)
        !s_plus= (-b + delta)/(2._ki*a)
        !s_moins = (-b - delta)/(2._ki*a)
        !!
        !epsilon_gv = un
        !delta_gv = -sign(un,aimag(arg_c))
        !!
        !eta_part = eta(-s_plus,-s_moins)-eta(z1-s_plus,z1-s_moins) &
        !& -eta(a-i_*small*epsilon_gv,1._ki/(a-i_*small*delta_gv))
        !!
        !eta_part = eta_part*log((z1-1._ki)/z1)
        !!
        !rest = cdilog(z1/(z1-s_plus)) - cdilog((z1-1._ki)/(z1-s_plus)) &
        !&+ eta(-s_plus,1._ki/(z1-s_plus))*log(z1/(z1-s_plus)) &
        !&- eta(1._ki-s_plus,1._ki/(z1-s_plus))*log((z1-1._ki)/(z1-s_plus)) &
        !&+ cdilog(z1/(z1-s_moins)) - cdilog((z1-1._ki)/(z1-s_moins)) &
        !&+ eta(-s_moins,1._ki/(z1-s_moins))*log(z1/(z1-s_moins)) &
        !&- eta(1._ki-s_moins,1._ki/(z1-s_moins))*log((z1-1._ki)/(z1-s_moins)) 
        !!
        !int_log = rest + eta_part
        !!
      !end if
      !!
    !end function int_log
    !
   function eval_numer_g0(u)
      !
      real(ki), intent (in) :: u
      complex(ki) :: eval_numer_g0
      !
      real(ki) :: x,y
      complex(ki) :: z1,jacob1
      complex(ki) :: z2,jacob2
      !
      eval_numer_g0 = ( & 
      &+log_quad(real(g_glob(1,1),ki),g_glob(1,2),g_glob(1,3),1,par1_glob,par2_glob,par3_glob,dim_glob,u) &
      &+log_quad(real(g_glob(2,1),ki),g_glob(2,2),g_glob(2,3),2,par1_glob,par2_glob,par3_glob,dim_glob,u) &
      &+log_quad(real(g_glob(3,1),ki),g_glob(3,2),g_glob(3,3),3,par1_glob,par2_glob,par3_glob,dim_glob,u) &
      & +rat_part(par1_glob,par2_glob,par3_glob,dim_glob,u))
      !
    end function eval_numer_g0
    !
    ! integrale du type :
    ! g(u)*log(a*u^2+b*u+c-i*lamb)/D(u)
    ! lorsque z = u + i*f(u)*eps, la partie imaginaire de l'argument du log
    ! est : (2*a*u+b)*eps*u*f(u) - lamb
    ! 2*a*u+b > 0 si u > \Sigma et a > 0 ou u < \Sigma et a < 0
    ! 2*a*u+b < 0 si u > \Sigma et a < 0 ou u < \Sigma et a > 0
    ! il faut donc que eps*f(u) soit de signe oppose a 2*a*u+b
    ! Or le Re(pole) est dehors de [0,1], donc seule la coupure du log importe
    ! Ca marche aussi si nous avons des masses imaginaires : dans ce cas
    ! la partie imaginaire de l'argument du log devient:
    ! (2*a*u+Re(b))*eps*u*f(u) + Im(b)*u+Im(c)
    ! Mais les masses internes ont toutes des parties imaginaires de meme signe (<0)
    ! donc Im(b)*u+Im(c) est < 0 quand u varie entre 0 et 1, on est donc
    ! ramene au cas precedent
    !
      function log_quad(a,b,c,flag,par1,par2,par3,dim,u)
      !
      real(ki), intent(in) :: a
      complex(ki), intent(in) :: b,c
      real(ki), intent(in) :: u
      integer, intent(in) :: flag
      integer, intent (in) :: par1,par2,par3
      character (len=3), intent (in) :: dim
      complex(ki) :: log_quad
      !
      real(ki) :: small
      real(ki) :: x,y
      complex(ki) :: z,jacob
      real(ki) :: sigma
      !
      !
      small = delta_glob
      !
      sigma = -real(b,ki)/a/2._ki
      !
      x = u
      !
      if ( (sigma <= 1._ki) .and. (sigma >= 0._ki) ) then
        !
        y = lambda_par*sign(un,a)*u*(u-1._ki)*(u-sigma)
        z = x + i_*y
        jacob = 1._ki + i_*lambda_par*sign(un,a)*( (u-1._ki)*(u-sigma) + u*(u-1._ki) + u*(u-sigma) )
        !
      else
        !
        y = lambda_par*sign(un,a*sigma)*u*(u-1._ki)
        z = x - i_*y
        jacob = 1._ki - i_*lambda_par*sign(un,a*sigma)*( (u-1._ki) + u )
        !
      end if
      !
      log_quad = ( log(a*z*z+b*z+c)*fg(z,flag,par1,par2,par3,dim) )*jacob
      !
      !
    end function log_quad
    !
    function rat_part(par1,par2,par3,dim,u)
      !
      real(ki), intent(in) :: u
      integer, intent (in) :: par1,par2,par3
      character (len=3), intent (in) :: dim
      complex(ki) :: rat_part
      !
      real(ki) :: x,y
      complex(ki) :: z,jacob
      !
      ! for this part, the sign of the contour does matter
      ! because there is no pole in it
      !
      x = u
      y = lambda_par*u*(u-1._ki)
      z = x - i_*y
      jacob = 1._ki - i_*lambda_par*( (u-1._ki) + u )
      !
      rat_part = fg(z,4,par1,par2,par3,dim)*jacob
      !
    end function rat_part
    !
    function fg(z,flag,par1,par2,par3,dim)
      !
      complex(ki), intent(in) :: z
      integer, intent(in) :: flag
      integer, intent (in) :: par1,par2,par3
      character (len=3), intent (in) :: dim
      complex(ki) :: fg
      !
      integer, dimension(3) :: par
      integer :: nb_par
      complex(ki) :: g1,g2,g3
      complex(ki) :: den1,den2,den3
      !
      complex(ki) ::&
            &tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,&
            &tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20,&
            &tmp21,tmp22,tmp23,tmp24,tmp25,tmp26,tmp27,tmp28,tmp29,tmp30,&
            &tmp31,tmp32,tmp33,tmp34,tmp35,tmp36,tmp37,tmp38,tmp39,tmp40,&
            &tmp41,tmp42,tmp43,tmp44,tmp45,tmp46,tmp47,tmp48,tmp49,tmp50,&
            &tmp51,tmp52,tmp53,tmp54,tmp55,tmp56,tmp57,tmp58,tmp59,tmp60,&
            &tmp61,tmp62
      !
      par = (/par1,par2,par3/)
      nb_par = count(mask=par/=0)
      !
      g1 = b*z*z+(c+e)*z+a+d+f
      g2 = a*z*z+d*z+f
      g3 = (a+b+c)*z*z+(d+e)*z+f
      !
      den1 = 2._ki*det_g*g1+det_s
      den2 = 2._ki*det_g*g2+det_s
      den3 = 2._ki*det_g*g3+det_s
      !
      if (dim == "ndi") then
        !
        if (nb_par == 0) then
          !
          select case (flag)
            !
            case(1)
            !
            fg = -num_b1/(2._ki*det_g*g1+det_s)
            !
            case(2)
            !
            fg = -num_b2/(2._ki*det_g*g2+det_s)
            !
            case(3)
            !
            fg = -num_b3/(2._ki*det_g*g3+det_s)
            !
            case(4)
            !
            fg = 0._ki
            !
          end select
          !
        else if (nb_par == 1) then
          !
          select case(par3)
          !
          case(1)
            !
            select case(flag)
            !
            case(1)
              !
              fg=2._ki*b*g1/den1-2._ki*num_b1*g1*(-c**2+4._ki*a*b+2._ki*b*d-c*e&
                &)/den1**2
              !
            case(2)
              !
              fg=(-num_b2+c*g2+num_b2*z)/den2-2._ki*num_b2*g2*(-det_g+4._ki*a*b&
                &-c**2+det_g*z-c*e+2._ki*b*d)/den2**2
              !
            case(3)
              !
              fg=(-num_b3-g3*c-2._ki*g3*b+num_b3*z)/den3-2._ki*num_b3*g3*(-det_&
                &g+4._ki*a*b-c**2+det_g*z-c*e+2._ki*b*d)/den3**2
              !
            case(4)
              !
              fg=num_b3*(-1._ki+z)/den3+num_b2*(-1._ki+z)/den2
              !
            end select
            !
          case(2)
            !
            select case(flag)
            !
            case(1)
              !
              fg=(c*g1-num_b1*z)/den1+2._ki*num_b1*g1*(det_g*z-c*d+2._ki*a*e)/d&
                &en1**2
              !
            case(2)
              !
              fg=2._ki*a*g2/den2+2._ki*num_b2*(-c*d+2._ki*a*e)*g2/den2**2
              !
            case(3)
              !
              fg=(-g3*c-2._ki*g3*a-num_b3*z)/den3+2._ki*num_b3*g3*(det_g*z-c*d+&
                &2._ki*a*e)/den3**2
              !
            case(4)
              !
              fg=-num_b3*z/den3-num_b1*z/den1
              !
            end select
            !
          case(3)
            !
            select case(flag)
            !
            case(1)
              !
              fg=(-2._ki*b*g1-num_b1-c*g1+num_b1*z)/den1-2._ki*num_b1*g1*(-det_&
                &g+c*e-2._ki*b*d+det_g*z-c*d+2._ki*a*e)/den1**2
              !
            case(2)
              !
              fg=(-c*g2-num_b2*z-2._ki*a*g2)/den2-2._ki*num_b2*g2*(-det_g*z+c*e&
                &-2._ki*b*d-c*d+2._ki*a*e)/den2**2
              !
            case(3)
              !
              fg=2._ki*g3*(c+a+b)/den3-2._ki*num_b3*g3*(-c*d-2._ki*b*d+c*e+2._k&
                &i*a*e)/den3**2
              !
            case(4)
              !
              fg=-num_b2*z/den2+num_b1*(-1._ki+z)/den1
              !
            end select
            !
          end select
          !
        else if (nb_par == 2) then
        !
        select case(par2)
        !
        case(1)
            !
            select case(par3)
            !
            case(1)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=-2._ki*num_b1*g1*(-4._ki*b*g1+4._ki*b*d+4._ki*b*f+4._ki*b*a-2.&
                      &_ki*c*e-c**2-e**2)/den1**2
                    !
                case(2)
                    !
                    fg=(-num_b2+2._ki*c*g2+2._ki*e*g2+2._ki*num_b2*z-num_b2*z**2)/den&
                      &2-2._ki*num_b2*g2*(-4._ki*b*g2+2._ki*det_g*z-2._ki*c*e+4._ki*b*&
                      &d-det_g+4._ki*b*a-c**2-det_g*z**2+4._ki*b*f-e**2)/den2**2
                    !
                case(3)
                    !
                    fg=(-num_b3-2._ki*c*g3-4._ki*b*g3-2._ki*e*g3+2._ki*num_b3*z-num_b&
                      &3*z**2)/den3-2._ki*num_b3*g3*(-4._ki*b*g3+2._ki*det_g*z-2._ki*c&
                      &*e+4._ki*b*d-det_g+4._ki*b*a-c**2-det_g*z**2+4._ki*b*f-e**2)/de&
                      &n3**2
                    !
                case(4)
                    !
                    fg=-num_b3*(z-1._ki)**2/den3-num_b2*(z-1._ki)**2/den2
                    !
                end select
                !
            case(2)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=-g1*(e+c)/den1+2._ki*num_b1*g1*(2._ki*c*g1-c*d+2._ki*e*a+e*d-2&
                      &._ki*c*f)/den1**2
                    !
                case(2)
                    !
                    fg=g2*(2._ki*a+d)/den2+2._ki*num_b2*g2*(2._ki*c*g2-c*d+2._ki*e*a+&
                      &e*d-2._ki*c*f)/den2**2
                    !
                case(3)
                    !
                    fg=(-2._ki*a*g3-c*g3-num_b3*z-d*g3+e*g3+num_b3*z**2)/den3+2._ki*n&
                      &um_b3*g3*(2._ki*c*g3+det_g*z-c*d+2._ki*e*a-det_g*z**2+e*d-2._ki&
                      &*c*f)/den3**2
                    !
                case(4)
                    !
                    fg=num_b3*z*(z-1._ki)/den3
                    !
                end select
                !
            case(3)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=g1*(c+e+2._ki*b)/den1-2._ki*num_b1*g1*(2._ki*c*g1+4._ki*b*g1-c&
                      &*d+2._ki*e*a+c*e-2._ki*b*d+e*d-4._ki*b*f+e**2-2._ki*c*f)/den1**&
                      &2
                    !
                case(2)
                    !
                    fg=(-c*g2-2._ki*e*g2-num_b2*z+num_b2*z**2-2._ki*a*g2-d*g2)/den2-2&
                      &._ki*num_b2*g2*(2._ki*c*g2+4._ki*b*g2-det_g*z+c*e-2._ki*b*d-c*d&
                      &+det_g*z**2-4._ki*b*f+e**2-2._ki*c*f+2._ki*e*a+e*d)/den2**2
                    !
                case(3)
                    !
                    fg=g3*(2._ki*c+2._ki*b+e+2._ki*a+d)/den3-2._ki*num_b3*g3*(2._ki*c&
                      &*g3+4._ki*b*g3-c*d+2._ki*e*a+c*e-2._ki*b*d+e*d-4._ki*b*f+e**2-2&
                      &._ki*c*f)/den3**2
                    !
                case(4)
                    !
                    fg=num_b2*z*(z-1._ki)/den2
                    !
                end select
                !
            end select
            !
        case(2)
            !
            select case(par3)
            !
            case(1)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=-g1*(e+c)/den1+2._ki*num_b1*g1*(2._ki*c*g1-c*d+2._ki*e*a+e*d-2&
                      &._ki*c*f)/den1**2
                    !
                case(2)
                    !
                    fg=g2*(2._ki*a+d)/den2+2._ki*num_b2*g2*(2._ki*c*g2-c*d+2._ki*e*a+&
                      &e*d-2._ki*c*f)/den2**2
                    !
                case(3)
                    !
                    fg=(-2._ki*a*g3-c*g3-num_b3*z-d*g3+e*g3+num_b3*z**2)/den3+2._ki*n&
                      &um_b3*g3*(2._ki*c*g3+det_g*z-c*d+2._ki*e*a-det_g*z**2+e*d-2._ki&
                      &*c*f)/den3**2
                    !
                case(4)
                    !
                    fg=num_b3*z*(z-1._ki)/den3
                    !
                end select
                !
            case(2)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=(-2._ki*g1*d-4._ki*a*g1-num_b1*z**2)/den1-2._ki*num_b1*g1*(-4.&
                      &_ki*a*g1-det_g*z**2+4._ki*a*f-d**2)/den1**2
                    !
                case(2)
                    !
                    fg=-2._ki*num_b2*g2*(-4._ki*a*g2+4._ki*a*f-d**2)/den2**2
                    !
                case(3)
                    !
                    fg=(2._ki*d*g3-num_b3*z**2)/den3-2._ki*num_b3*g3*(-4._ki*a*g3-det&
                      &_g*z**2+4._ki*a*f-d**2)/den3**2
                    !
                case(4)
                    !
                    fg=-num_b3*z**2/den3-num_b1*z**2/den1
                    !
                end select
                !
            case(3)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=(g1*e+2._ki*c*g1-num_b1*z+2._ki*g1*d+4._ki*a*g1+num_b1*z**2)/d&
                      &en1+2._ki*num_b1*g1*(-2._ki*c*g1-4._ki*a*g1+det_g*z-e*d+2._ki*c&
                      &*f-det_g*z**2+4._ki*a*f-d**2)/den1**2
                    !
                case(2)
                    !
                    fg=-d*g2/den2+2._ki*num_b2*g2*(-2._ki*c*g2-4._ki*a*g2+4._ki*a*f-e&
                      &*d+2._ki*c*f-d**2)/den2**2
                    !
                case(3)
                    !
                    fg=-g3*(e+d)/den3+2._ki*num_b3*g3*(-2._ki*c*g3-4._ki*a*g3+4._ki*a&
                      &*f-e*d+2._ki*c*f-d**2)/den3**2
                    !
                case(4)
                    !
                    fg=num_b1*z*(z-1._ki)/den1
                    !
                end select
                !
            end select
            !
        case(3)
            !
            select case(par3)
            !
            case(1)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=g1*(c+e+2._ki*b)/den1-2._ki*num_b1*g1*(2._ki*c*g1+4._ki*b*g1-c&
                      &*d+2._ki*e*a+c*e-2._ki*b*d+e*d-4._ki*b*f+e**2-2._ki*c*f)/den1**&
                      &2
                    !
                case(2)
                    !
                    fg=(-c*g2-2._ki*e*g2-num_b2*z+num_b2*z**2-2._ki*a*g2-d*g2)/den2-2&
                      &._ki*num_b2*g2*(2._ki*c*g2+4._ki*b*g2-det_g*z+c*e-2._ki*b*d-c*d&
                      &+det_g*z**2-4._ki*b*f+e**2-2._ki*c*f+2._ki*e*a+e*d)/den2**2
                    !
                case(3)
                    !
                    fg=g3*(2._ki*c+2._ki*b+e+2._ki*a+d)/den3-2._ki*num_b3*g3*(2._ki*c&
                      &*g3+4._ki*b*g3-c*d+2._ki*e*a+c*e-2._ki*b*d+e*d-4._ki*b*f+e**2-2&
                      &._ki*c*f)/den3**2
                    !
                case(4)
                    !
                    fg=num_b2*z*(z-1._ki)/den2
                    !
                end select
                !
            case(2)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=(g1*e+2._ki*c*g1-num_b1*z+2._ki*g1*d+4._ki*a*g1+num_b1*z**2)/d&
                      &en1+2._ki*num_b1*g1*(-2._ki*c*g1-4._ki*a*g1+det_g*z-e*d+2._ki*c&
                      &*f-det_g*z**2+4._ki*a*f-d**2)/den1**2
                    !
                case(2)
                    !
                    fg=-d*g2/den2+2._ki*num_b2*g2*(-2._ki*c*g2-4._ki*a*g2+4._ki*a*f-e&
                      &*d+2._ki*c*f-d**2)/den2**2
                    !
                case(3)
                    !
                    fg=-g3*(e+d)/den3+2._ki*num_b3*g3*(-2._ki*c*g3-4._ki*a*g3+4._ki*a&
                      &*f-e*d+2._ki*c*f-d**2)/den3**2
                    !
                case(4)
                    !
                    fg=num_b1*z*(z-1._ki)/den1
                    !
                end select
                !
            case(3)
                !
                select case(flag)
                !
                case(1)
                    !
                    fg=(-2._ki*g1*e-4._ki*c*g1-4._ki*b*g1-2._ki*g1*d-4._ki*a*g1-num_b&
                      &1*z**2+2._ki*num_b1*z-num_b1)/den1-2._ki*num_b1*g1*(-4._ki*b*g1&
                      &-4._ki*c*g1-4._ki*a*g1-det_g+4._ki*b*f-e**2+2._ki*det_g*z-2._ki&
                      &*e*d+4._ki*c*f-det_g*z**2+4._ki*a*f-d**2)/den1**2
                    !
                case(2)
                    !
                    fg=(2._ki*e*g2-num_b2*z**2+2._ki*d*g2)/den2-2._ki*num_b2*g2*(-4._&
                      &ki*b*g2-4._ki*c*g2-4._ki*a*g2+4._ki*a*f-det_g*z**2+4._ki*b*f-e*&
                      &*2-2._ki*e*d+4._ki*c*f-d**2)/den2**2
                    !
                case(3)
                    !
                    fg=-2._ki*num_b3*g3*(-4._ki*b*g3-4._ki*c*g3-4._ki*a*g3-2._ki*e*d+&
                      &4._ki*a*f+4._ki*c*f+4._ki*b*f-d**2-e**2)/den3**2
                    !
                case(4)
                    !
                    fg=-num_b2*z**2/den2-num_b1*(z-1._ki)**2/den1
                    !
                end select
                !
            end select
            !
        end select
        !
        else if (nb_par == 3) then
        !
         select case(par1)
         !
         case(1)
            !
            select case(par2)
            !
            case(1)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*d
                     tmp2 = b*f
                     tmp3 = b*g1
                     tmp4 = b*b
                     tmp5 = c*e
                     tmp6 = e*e
                     tmp7 = den1*den1
                     fg = (2.0_ki*(-(2.0_ki*(4.0_ki*(tmp2-tmp3)-(det_g+tmp6))*(2.0_ki*tmp1+det_g-tmp5&
                     &)/(tmp7*den1)*g1*num_b1+(3.0_ki*(b*det_g*g1-2.0_ki*tmp3*num_b1)+2.0_ki*(-(b*c*e&
                     &*g1+tmp6*num_b1))+4.0_ki*(2.0_ki*(tmp1*num_b1+tmp2*num_b1)+tmp3*tmp3+tmp4*d*g1+&
                     &a*b*num_b1-(tmp4*f*g1+tmp5*num_b1))+tmp3*tmp6+det_g*num_b1-c*c*num_b1)/tmp7))*g&
                     &1)
                     !
                  case(2)
                     !
                     tmp1 = b*f
                     tmp2 = tmp1*num_b2
                     tmp3 = b*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = b*d
                     tmp6 = tmp5*num_b2
                     tmp7 = e*e
                     tmp8 = c*e
                     tmp9 = tmp8*num_b2
                     tmp10 = det_g*num_b2
                     tmp11 = tmp10*z
                     tmp12 = z*z
                     tmp13 = b*c
                     tmp14 = c*g2
                     tmp15 = num_b2*z
                     tmp16 = den2*den2
                     tmp17 = e*g2
                     fg = ((-(4.0_ki*(4.0_ki*(tmp1-tmp3)-(tmp12*det_g+tmp7))*(2.0_ki*tmp5+det_g*z-tmp&
                     &8)/(tmp16*den2)*g2*num_b2+(24.0_ki*(tmp2-tmp4)+12.0_ki*(tmp6+tmp4*z)+8.0_ki*(a*&
                     &b*num_b2-tmp2*z)+6.0_ki*(tmp11-(tmp7*num_b2+tmp10*tmp12+tmp9))+2.0_ki*(tmp15*tm&
                     &p7-(det_g*e*g2*z+c*c*num_b2+tmp12*tmp9+tmp10))+4.0_ki*(tmp11*tmp12+tmp12*tmp6+t&
                     &mp13*g2*g2-(tmp5*e*g2+tmp13*f*g2))+3.0_ki*tmp14*tmp7+tmp12*c*det_g*g2)/tmp16))*&
                     &g2+(3.0_ki*(tmp14+tmp15+2.0_ki*tmp17-tmp12*num_b2)+tmp12*tmp14+tmp12*tmp15-(2.0&
                     &_ki*tmp17*z+num_b2))/den2)
                     !
                  case(3)
                     !
                     tmp1 = b*f
                     tmp2 = tmp1*num_b3
                     tmp3 = b*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = b*d
                     tmp6 = tmp5*num_b3
                     tmp7 = e*e
                     tmp8 = c*e
                     tmp9 = tmp8*num_b3
                     tmp10 = det_g*num_b3
                     tmp11 = tmp10*z
                     tmp12 = z*z
                     tmp13 = b*c
                     tmp14 = c*g3
                     tmp15 = num_b3*z
                     tmp16 = den3*den3
                     tmp17 = e*g3
                     fg = ((-(4.0_ki*(4.0_ki*(tmp1-tmp3)-(tmp12*det_g+tmp7))*(2.0_ki*tmp5+det_g*z-tmp&
                     &8)/(tmp16*den3)*g3*num_b3+(24.0_ki*(tmp2-tmp4)+12.0_ki*(tmp6+tmp4*z)+6.0_ki*(tm&
                     &p11-(tmp7*num_b3+tmp10*tmp12+tmp9))+8.0_ki*(a*b*num_b3+b*b*f*g3-(tmp2*z+tmp3*tm&
                     &p3))+4.0_ki*(tmp11*tmp12+tmp12*tmp6+tmp13*f*g3+tmp5*e*g3-tmp13*g3*g3)+2.0_ki*(t&
                     &mp15*tmp7+det_g*e*g3*z-(tmp12*b*det_g*g3+c*c*num_b3+tmp3*tmp7+tmp12*tmp9+tmp10)&
                     &)-(tmp12*c*det_g*g3+3.0_ki*tmp14*tmp7))/tmp16))*g3+(6.0_ki*(-(tmp3+tmp17))+2.0_&
                     &ki*(tmp17*z-tmp12*tmp3)+3.0_ki*(tmp15-(tmp12*num_b3+tmp14))+tmp12*tmp15-(tmp12*&
                     &tmp14+num_b3))/den3)
                     !
                  case(4)
                     !
                     tmp1 = b*g1
                     tmp2 = num_b2*z
                     tmp3 = z*z
                     tmp4 = num_b3*z
                     tmp5 = b*g3
                     tmp6 = b*f
                     tmp7 = b*d
                     tmp8 = e*e
                     tmp9 = c*e
                     tmp10 = 2.0_ki*(tmp3*det_g+tmp7*z)
                     tmp11 = tmp9*z
                     fg = ((3.0_ki*tmp1+1.0_ki/2.0_ki*num_b1)/den1+(3.0_ki*(tmp2-tmp3*num_b2)+3.0_ki/&
                     &2.0_ki*tmp2*tmp3+1.0_ki/2.0_ki*tmp3*c*g2-(e*g2*z+num_b2))/den2+(3.0_ki*(tmp4-tm&
                     &p3*num_b3)+3.0_ki/2.0_ki*tmp3*tmp4+e*g3*z-(1.0_ki/2.0_ki*tmp3*c*g3+tmp3*tmp5+nu&
                     &m_b3))/den3-((4.0_ki*(b*g2-tmp6)+tmp10+tmp8-tmp11)/(den2*den2)*g2*num_b2*z+(4.0&
                     &_ki*(tmp5-tmp6)+tmp10+tmp8-tmp11)/(den3*den3)*g3*num_b3*z+(4.0_ki*(tmp1-tmp6)+2&
                     &.0_ki*(tmp7+det_g)+tmp8-tmp9)/(den1*den1)*g1*num_b1))
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*c
                     tmp2 = c*f
                     tmp3 = b*d
                     tmp4 = c*g1
                     tmp5 = e*g1
                     tmp6 = d*e
                     tmp7 = c*e
                     tmp8 = den1*den1
                     fg = (((4.0_ki*(tmp1*f*g1+tmp3*num_b1*z+a*e*num_b1-(tmp3*e*g1+tmp1*d*g1+tmp2*num&
                     &_b1+tmp1*g1*g1))+2.0_ki*(tmp4*num_b1+tmp5*c*c+tmp6*num_b1+det_g*num_b1*z+b*det_&
                     &g*g1*z-(c*det_g*g1+c*d*num_b1+tmp7*num_b1*z))+tmp4*e*e-(det_g*e*g1+tmp5*num_b1)&
                     &)/tmp8-(4.0_ki*(2.0_ki*(tmp2-tmp4)+det_g*z-tmp6)*(2.0_ki*tmp3+det_g-tmp7)/(tmp8&
                     &*den1)*g1*num_b1+(c+e+2.0_ki*b*z)/den1))*g1)
                     !
                  case(2)
                     !
                     tmp1 = c*f
                     tmp2 = tmp1*num_b2
                     tmp3 = c*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = d*e
                     tmp6 = tmp5*num_b2
                     tmp7 = c*d
                     tmp8 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(4.0_ki*(a+d)-(2.0_ki*d*z+g2))/den2+(4.0_ki*(2.0_ki*(tmp4-t&
                     &mp2)+tmp6+tmp2*z+a*e*num_b2)+2.0_ki*(b*d*d*g2+c*c*f*g2-(tmp7*e*g2+tmp7*num_b2+t&
                     &mp6*z+3.0_ki*tmp4*z+tmp3*tmp3))+d*det_g*g2*z-e*g2*num_b2)/tmp8-4.0_ki*(2.0_ki*(&
                     &tmp1-tmp3)-tmp5)*(2.0_ki*b*d+det_g*z-c*e)/(tmp8*den2)*g2*num_b2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = c*f
                     tmp2 = tmp1*num_b3
                     tmp3 = c*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = b*c
                     tmp6 = g3*g3
                     tmp7 = d*e
                     tmp8 = tmp7*num_b3
                     tmp9 = det_g*num_b3
                     tmp10 = z*z
                     tmp11 = tmp9*z
                     tmp12 = b*d
                     tmp13 = b*g3
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = e*g3
                     tmp17 = den3*den3
                     tmp18 = d*g3
                     tmp19 = num_b3*z
                     fg = (((8.0_ki*(tmp4-tmp2)+4.0_ki*(tmp8+tmp10*tmp11+tmp2*z+tmp5*tmp6+tmp10*tmp12&
                     &*num_b3+tmp12*e*g3+a*e*num_b3-(tmp5*f*g3+tmp10*tmp9))+2.0_ki*(tmp11+tmp3*tmp3+t&
                     &mp14*e*g3-(tmp10*b*det_g*g3+c*c*f*g3+tmp10*tmp15*num_b3+tmp8*z+3.0_ki*tmp4*z+tm&
                     &p14*num_b3+tmp13*d*d))+det_g*e*g3*z-(d*det_g*g3*z+tmp10*c*det_g*g3+tmp3*e*e+tmp&
                     &16*num_b3))/tmp17-4.0_ki*(2.0_ki*(tmp1-tmp3)+tmp10*det_g-tmp7)*(2.0_ki*tmp12+de&
                     &t_g*z-tmp15)/(tmp17*den3)*g3*num_b3)*g3+(2.0_ki*(tmp16+tmp10*tmp13+tmp10*num_b3&
                     &-(a*g3+tmp18))+1.0_ki/2.0_ki*tmp6+tmp10*tmp3+tmp18*z-(tmp16*z+tmp10*tmp19+tmp3+&
                     &tmp19))/den3)
                     !
                  case(4)
                     !
                     tmp1 = c*g1
                     tmp2 = c*f
                     tmp3 = b*d*z
                     tmp4 = d*e
                     tmp5 = c*e*z
                     tmp6 = c*g3
                     tmp7 = z*z
                     fg = (((2.0_ki*(tmp2-c*g2)-tmp4)/(den2*den2)*num_b2*z-1.0_ki/2.0_ki*d/den2*z)*g2&
                     &+(tmp1+1.0_ki/2.0_ki*e*g1-(b*g1*z+1.0_ki/2.0_ki*num_b1*z))/den1+1.0_ki/2.0_ki*(&
                     &2.0_ki*(2.0_ki*num_b3*z+b*g3*z-num_b3)+tmp6*z+d*g3-(e*g3+3.0_ki*tmp7*num_b3))/d&
                     &en3*z+(2.0_ki*(tmp2+tmp3+det_g*z-tmp1)-(tmp5+tmp4))/(den1*den1)*g1*num_b1+(2.0_&
                     &ki*(tmp2+tmp3+tmp7*det_g-tmp6)-(tmp5+tmp4))/(den3*den3)*g3*num_b3*z)
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*g1
                     tmp2 = b*d
                     tmp3 = tmp2*num_b1
                     tmp4 = b*f
                     tmp5 = b*b
                     tmp6 = b*c
                     tmp7 = c*e
                     tmp8 = tmp7*num_b1
                     tmp9 = c*f
                     tmp10 = det_g*num_b1
                     tmp11 = e*e
                     tmp12 = b*det_g*g1
                     tmp13 = c*g1
                     tmp14 = e*g1
                     tmp15 = d*e
                     tmp16 = den1*den1
                     fg = (((c+e+2.0_ki*b*z)/den1+4.0_ki*(2.0_ki*tmp2+det_g-tmp7)*(4.0_ki*(tmp4-tmp1)&
                     &+2.0_ki*(tmp9-tmp13)+det_g*z-(det_g+tmp15+tmp11))/(tmp16*den1)*g1*num_b1-(8.0_k&
                     &i*(tmp5*f*g1-(tmp5*d*g1+tmp4*num_b1+tmp3+tmp1*tmp1))+4.0_ki*(tmp8+tmp1*num_b1+t&
                     &mp3*z+tmp6*e*g1+tmp6*f*g1+a*e*num_b1-(tmp6*d*g1+tmp2*e*g1+tmp9*num_b1+tmp6*g1*g&
                     &1))+2.0_ki*(tmp10*z+tmp11*num_b1+tmp12*z+tmp13*num_b1+tmp14*c*c+tmp15*num_b1-(c&
                     &*det_g*g1+c*d*num_b1+tmp8*z+tmp1*tmp11+3.0_ki*tmp12+tmp10))+tmp11*tmp13-(det_g*&
                     &e*g1+tmp14*num_b1))/tmp16)*g1)
                     !
                  case(2)
                     !
                     tmp1 = b*f
                     tmp2 = tmp1*num_b2
                     tmp3 = b*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = c*f
                     tmp6 = tmp5*num_b2
                     tmp7 = c*g2
                     tmp8 = tmp7*num_b2
                     tmp9 = e*e
                     tmp10 = b*c
                     tmp11 = g2*g2
                     tmp12 = b*d
                     tmp13 = tmp12*num_b2
                     tmp14 = d*e
                     tmp15 = tmp14*num_b2
                     tmp16 = det_g*num_b2
                     tmp17 = z*z
                     tmp18 = tmp16*z
                     tmp19 = c*d
                     tmp20 = c*e
                     tmp21 = tmp20*num_b2
                     tmp22 = num_b2*z
                     tmp23 = e*g2
                     tmp24 = den2*den2
                     tmp25 = d*g2
                     fg = ((4.0_ki*(2.0_ki*tmp12+det_g*z-tmp20)*(4.0_ki*(tmp1-tmp3)+2.0_ki*(tmp5-tmp7&
                     &)-(tmp17*det_g+tmp9+tmp14))/(tmp24*den2)*g2*num_b2-(16.0_ki*(tmp4-tmp2)+6.0_ki*&
                     &(-(tmp8*z+2.0_ki*tmp4*z))+8.0_ki*(tmp8+tmp2*z-tmp6)+4.0_ki*(tmp15+tmp16*tmp17+t&
                     &mp6*z+tmp9*num_b2+tmp10*f*g2+tmp12*e*g2+a*e*num_b2-(tmp17*tmp18+tmp13*tmp17+tmp&
                     &10*tmp11+tmp13))+2.0_ki*(tmp21+tmp17*tmp21+tmp3*d*d+c*c*f*g2+det_g*e*g2*z-(tmp1&
                     &9*e*g2+tmp22*tmp9+tmp19*num_b2+tmp15*z+tmp7*tmp7+tmp18))+d*det_g*g2*z-(tmp17*c*&
                     &det_g*g2+3.0_ki*tmp7*tmp9+tmp23*num_b2))/tmp24)*g2+(2.0_ki*(tmp17*num_b2+tmp23*&
                     &z-(a*g2+tmp25+2.0_ki*tmp23))+1.0_ki/2.0_ki*tmp11+tmp25*z-(tmp17*tmp7+tmp17*tmp2&
                     &2+tmp7+tmp22))/den2)
                     !
                  case(3)
                     !
                     tmp1 = b*f
                     tmp2 = tmp1*num_b3
                     tmp3 = b*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = b*c
                     tmp6 = f*g3
                     tmp7 = c*f
                     tmp8 = tmp7*num_b3
                     tmp9 = c*g3
                     tmp10 = tmp9*num_b3
                     tmp11 = e*e
                     tmp12 = b*d
                     tmp13 = d*e
                     tmp14 = tmp13*num_b3
                     tmp15 = c*d
                     tmp16 = c*e
                     tmp17 = den3*den3
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(-(e*z+d*z))+4.0_ki*(a+b+c+d+e)-g3)/den3+4.0_ki*(2.&
                     &0_ki*(2.0_ki*(tmp1-tmp3)+tmp7-tmp9)-(tmp13+tmp11))*(2.0_ki*tmp12+det_g*z-tmp16)&
                     &/(tmp17*den3)*g3*num_b3-(4.0_ki*(tmp14+tmp11*num_b3+tmp8*z+a*e*num_b3-tmp12*num&
                     &_b3)+8.0_ki*(2.0_ki*(tmp4-tmp2)+tmp10+tmp3*tmp3+tmp2*z+tmp5*g3*g3-(tmp5*f*g3+tm&
                     &p6*b*b+tmp8))+2.0_ki*(3.0_ki*(-(2.0_ki*tmp4*z+tmp10*z))+tmp9*tmp9+tmp11*tmp3+tm&
                     &p11*tmp9+tmp16*num_b3+tmp15*e*g3-(tmp11*num_b3*z+tmp6*c*c+tmp3*d*d+tmp15*num_b3&
                     &+tmp14*z))-(det_g*e*g3*z+d*det_g*g3*z+e*g3*num_b3))/tmp17)*g3)
                     !
                  case(4)
                     !
                     tmp1 = b*g1
                     tmp2 = c*g1
                     tmp3 = b*f
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = tmp4*z
                     tmp7 = e*e
                     tmp8 = c*e
                     tmp9 = d*e
                     tmp10 = tmp8*z
                     tmp11 = c*g2
                     tmp12 = z*z
                     fg = ((-((4.0_ki*(tmp3-b*g3)+2.0_ki*(tmp5-c*g3)-(tmp9+tmp7))/(den3*den3)*num_b3*&
                     &z+1.0_ki/2.0_ki*(d+e)/den3*z))*g3+(tmp1*z+1.0_ki/2.0_ki*num_b1*z-(1.0_ki/2.0_ki&
                     &*e*g1+1.0_ki/2.0_ki*num_b1+tmp2+3.0_ki*tmp1))/den1+(2.0_ki*(2.0_ki*(b*g2-tmp3)+&
                     &tmp11+tmp6+tmp12*det_g-tmp5)+tmp7+tmp9-tmp10)/(den2*den2)*g2*num_b2*z-((2.0_ki*&
                     &(2.0_ki*(tmp3-tmp1)+tmp5+tmp6+det_g*z-(det_g+tmp4+tmp2))+tmp8-(tmp9+tmp7+tmp10)&
                     &)/(den1*den1)*g1*num_b1+1.0_ki/2.0_ki*(2.0_ki*(num_b2-(2.0_ki*num_b2*z+e*g2))+t&
                     &mp11*z+3.0_ki*tmp12*num_b2-d*g2)/den2*z))
                     !
                  end select
                  !
               end select
               !
            case(2)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*c
                     tmp2 = c*f
                     tmp3 = b*d
                     tmp4 = c*g1
                     tmp5 = e*g1
                     tmp6 = d*e
                     tmp7 = c*e
                     tmp8 = den1*den1
                     fg = (((4.0_ki*(tmp1*f*g1+tmp3*num_b1*z+a*e*num_b1-(tmp3*e*g1+tmp1*d*g1+tmp2*num&
                     &_b1+tmp1*g1*g1))+2.0_ki*(tmp4*num_b1+tmp5*c*c+tmp6*num_b1+det_g*num_b1*z+b*det_&
                     &g*g1*z-(c*det_g*g1+c*d*num_b1+tmp7*num_b1*z))+tmp4*e*e-(det_g*e*g1+tmp5*num_b1)&
                     &)/tmp8-(4.0_ki*(2.0_ki*(tmp2-tmp4)+det_g*z-tmp6)*(2.0_ki*tmp3+det_g-tmp7)/(tmp8&
                     &*den1)*g1*num_b1+(c+e+2.0_ki*b*z)/den1))*g1)
                     !
                  case(2)
                     !
                     tmp1 = c*f
                     tmp2 = tmp1*num_b2
                     tmp3 = c*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = d*e
                     tmp6 = tmp5*num_b2
                     tmp7 = c*d
                     tmp8 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(4.0_ki*(a+d)-(2.0_ki*d*z+g2))/den2+(4.0_ki*(2.0_ki*(tmp4-t&
                     &mp2)+tmp6+tmp2*z+a*e*num_b2)+2.0_ki*(b*d*d*g2+c*c*f*g2-(tmp7*e*g2+tmp7*num_b2+t&
                     &mp6*z+3.0_ki*tmp4*z+tmp3*tmp3))+d*det_g*g2*z-e*g2*num_b2)/tmp8-4.0_ki*(2.0_ki*(&
                     &tmp1-tmp3)-tmp5)*(2.0_ki*b*d+det_g*z-c*e)/(tmp8*den2)*g2*num_b2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = c*f
                     tmp2 = tmp1*num_b3
                     tmp3 = c*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = b*c
                     tmp6 = g3*g3
                     tmp7 = d*e
                     tmp8 = tmp7*num_b3
                     tmp9 = det_g*num_b3
                     tmp10 = z*z
                     tmp11 = tmp9*z
                     tmp12 = b*d
                     tmp13 = b*g3
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = e*g3
                     tmp17 = den3*den3
                     tmp18 = d*g3
                     tmp19 = num_b3*z
                     fg = (((8.0_ki*(tmp4-tmp2)+4.0_ki*(tmp8+tmp10*tmp11+tmp2*z+tmp5*tmp6+tmp10*tmp12&
                     &*num_b3+tmp12*e*g3+a*e*num_b3-(tmp5*f*g3+tmp10*tmp9))+2.0_ki*(tmp11+tmp3*tmp3+t&
                     &mp14*e*g3-(tmp10*b*det_g*g3+c*c*f*g3+tmp10*tmp15*num_b3+tmp8*z+3.0_ki*tmp4*z+tm&
                     &p14*num_b3+tmp13*d*d))+det_g*e*g3*z-(d*det_g*g3*z+tmp10*c*det_g*g3+tmp3*e*e+tmp&
                     &16*num_b3))/tmp17-4.0_ki*(2.0_ki*(tmp1-tmp3)+tmp10*det_g-tmp7)*(2.0_ki*tmp12+de&
                     &t_g*z-tmp15)/(tmp17*den3)*g3*num_b3)*g3+(2.0_ki*(tmp16+tmp10*tmp13+tmp10*num_b3&
                     &-(a*g3+tmp18))+1.0_ki/2.0_ki*tmp6+tmp10*tmp3+tmp18*z-(tmp16*z+tmp10*tmp19+tmp3+&
                     &tmp19))/den3)
                     !
                  case(4)
                     !
                     tmp1 = c*g1
                     tmp2 = c*f
                     tmp3 = b*d*z
                     tmp4 = d*e
                     tmp5 = c*e*z
                     tmp6 = c*g3
                     tmp7 = z*z
                     fg = (((2.0_ki*(tmp2-c*g2)-tmp4)/(den2*den2)*num_b2*z-1.0_ki/2.0_ki*d/den2*z)*g2&
                     &+(tmp1+1.0_ki/2.0_ki*e*g1-(b*g1*z+1.0_ki/2.0_ki*num_b1*z))/den1+1.0_ki/2.0_ki*(&
                     &2.0_ki*(2.0_ki*num_b3*z+b*g3*z-num_b3)+tmp6*z+d*g3-(e*g3+3.0_ki*tmp7*num_b3))/d&
                     &en3*z+(2.0_ki*(tmp2+tmp3+det_g*z-tmp1)-(tmp5+tmp4))/(den1*den1)*g1*num_b1+(2.0_&
                     &ki*(tmp2+tmp3+tmp7*det_g-tmp6)-(tmp5+tmp4))/(den3*den3)*g3*num_b3*z)
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = a*e
                     tmp3 = c*f
                     tmp4 = c*g1
                     tmp5 = d*g1
                     tmp6 = c*c
                     tmp7 = c*d
                     tmp8 = d*e
                     tmp9 = den1*den1
                     tmp10 = det_g*z
                     fg = (((3.0_ki*(2.0_ki*tmp4*num_b1*z+c*det_g*g1*z)+4.0_ki*(2.0_ki*(tmp1*num_b1-a&
                     &*f*num_b1)+a*c*e*g1-(tmp3*num_b1*z+tmp2*num_b1*z))+2.0_ki*(tmp1*e*e+d*d*num_b1+&
                     &tmp6*f*g1+tmp7*num_b1*z+tmp8*num_b1*z-(det_g*num_b1*z*z+tmp7*e*g1+tmp5*tmp6+tmp&
                     &4*tmp4))+tmp5*num_b1+det_g*e*g1*z)/tmp9+4.0_ki*(2.0_ki*(tmp3-tmp4)+tmp10-tmp8)*&
                     &(tmp10+2.0_ki*tmp2-tmp7)/(tmp9*den1)*g1*num_b1-1.0_ki/2.0_ki*(4.0_ki*(2.0_ki*a+&
                     &d)+2.0_ki*(3.0_ki*c*z+e*z)+g1)/den1)*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = d*d
                     tmp3 = c*g2
                     tmp4 = den2*den2
                     fg = (((8.0_ki*(a*g2*num_b2-a*f*num_b2)+4.0_ki*(tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2))&
                     &+tmp2*tmp3+2.0_ki*tmp2*num_b2+d*g2*num_b2)/tmp4+4.0_ki*(2.0_ki*(c*f-tmp3)-d*e)*&
                     &(2.0_ki*a*e-c*d)/(tmp4*den2)*g2*num_b2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*g3
                     tmp2 = a*c
                     tmp3 = g3*g3
                     tmp4 = det_g*num_b3
                     tmp5 = z*z
                     tmp6 = a*e
                     tmp7 = c*f
                     tmp8 = c*g3
                     tmp9 = d*d
                     tmp10 = c*d
                     tmp11 = d*e
                     tmp12 = d*g3
                     tmp13 = den3*den3
                     fg = ((4.0_ki*(2.0_ki*(tmp7-tmp8)+tmp5*det_g-tmp11)*(2.0_ki*tmp6+det_g*z-tmp10)/&
                     &(tmp13*den3)*g3*num_b3-(8.0_ki*(a*f*num_b3-tmp1*num_b3)+4.0_ki*(tmp2*f*g3+tmp4*&
                     &tmp5*z+tmp5*tmp6*num_b3+tmp7*num_b3*z-(a*d*e*g3+tmp2*tmp3))+2.0_ki*(tmp1*e*e+c*&
                     &c*f*g3+tmp5*a*det_g*g3-(3.0_ki*tmp8*num_b3*z+tmp11*num_b3*z+tmp10*e*g3+tmp10*tm&
                     &p5*num_b3+tmp9*num_b3+tmp4*tmp5+tmp8*tmp8))+tmp8*tmp9+tmp5*c*det_g*g3+det_g*e*g&
                     &3*z-(d*det_g*g3*z+tmp12*num_b3))/tmp13)*g3+(2.0_ki*(tmp12+tmp1*tmp5)+1.0_ki/2.0&
                     &_ki*tmp3+tmp5*tmp8+tmp5*num_b3*z+e*g3*z-(tmp5*num_b3+tmp12*z))/den3)
                     !
                  case(4)
                     !
                     tmp1 = a*e
                     tmp2 = c*f
                     tmp3 = c*g1
                     tmp4 = c*d
                     tmp5 = d*e
                     tmp6 = c*g3
                     tmp7 = z*z
                     fg = (1.0_ki/2.0_ki*(2.0_ki*(a*g3*z-num_b3*z)+tmp6*z+3.0_ki*tmp7*num_b3+e*g3-d*g&
                     &3)/den3*z-((2.0_ki*(tmp2+tmp1*z+tmp7*det_g-tmp6)-(tmp4*z+tmp5))/(den3*den3)*g3*&
                     &num_b3*z+(2.0_ki*(tmp1+tmp2+det_g*z-tmp3)-(tmp5+tmp4))/(den1*den1)*g1*num_b1*z+&
                     &1.0_ki/2.0_ki*(3.0_ki*tmp3+e*g1-num_b1*z)/den1*z))
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = b*c
                     tmp3 = a*e
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = c*g1
                     tmp7 = e*e
                     tmp8 = c*det_g*g1
                     tmp9 = tmp6*num_b1
                     tmp10 = d*g1
                     tmp11 = c*c
                     tmp12 = e*g1
                     tmp13 = det_g*num_b1
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = d*e
                     tmp17 = det_g*e*g1
                     tmp18 = den1*den1
                     tmp19 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(3.0_ki*c*z+e*z)+4.0_ki*(2.0_ki*a+d+b*z)+g1)/den1-(&
                     &4.0_ki*(2.0_ki*(tmp5-tmp6)+tmp19-tmp16)*(2.0_ki*(tmp3-tmp4)+tmp15+tmp19-(det_g+&
                     &tmp14))/(tmp18*den1)*g1*num_b1+(3.0_ki*(tmp8*z+2.0_ki*tmp9*z)+4.0_ki*(2.0_ki*(t&
                     &mp1*num_b1-a*f*num_b1)+tmp2*f*g1+tmp4*num_b1*z+a*c*e*g1-(tmp5*num_b1*z+tmp4*e*g&
                     &1+tmp3*num_b1*z+tmp2*d*g1+tmp2*g1*g1))+2.0_ki*(tmp1*tmp7+tmp11*tmp12+tmp13*z+d*&
                     &d*num_b1+tmp11*f*g1+tmp14*num_b1*z+tmp16*num_b1*z+b*det_g*g1*z-(tmp15*num_b1*z+&
                     &tmp14*e*g1+tmp13*z*z+tmp10*tmp11+tmp9+tmp8+tmp6*tmp6))+tmp10*num_b1+tmp17*z+tmp&
                     &6*tmp7-(tmp12*num_b1+tmp17))/tmp18))*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = c*f
                     tmp3 = tmp2*num_b2
                     tmp4 = c*g2
                     tmp5 = tmp4*num_b2
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = tmp7*num_b2
                     tmp9 = c*d
                     tmp10 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z-d)+g2)/den2-(4.0_ki*(2.0_ki*(tmp2-tmp4)-tmp7)*&
                     &(2.0_ki*(a*e-b*d)+c*e-(det_g*z+tmp9))/(tmp10*den2)*g2*num_b2+(4.0_ki*(2.0_ki*(a&
                     &*g2*num_b2-a*f*num_b2)+tmp5+tmp3*z+tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2+tmp3))+2.0_ki&
                     &*(tmp8+tmp6*num_b2+tmp6*b*g2+c*c*f*g2-(tmp9*e*g2+tmp8*z+3.0_ki*tmp5*z+tmp4*tmp4&
                     &))+tmp4*tmp6+d*g2*num_b2+d*det_g*g2*z-e*g2*num_b2)/tmp10))*g2)
                     !
                  case(3)
                     !
                     tmp1 = z*z
                     tmp2 = a*g3
                     tmp3 = c*g3
                     tmp4 = a*c
                     tmp5 = g3*g3
                     tmp6 = b*c
                     tmp7 = c*f
                     tmp8 = a*e
                     tmp9 = b*d
                     tmp10 = c*d
                     tmp11 = d*d
                     tmp12 = e*e
                     tmp13 = d*e
                     tmp14 = c*e
                     tmp15 = den3*den3
                     fg = (((8.0_ki*(a*f*num_b3-tmp2*num_b3)+2.0_ki*(tmp12*tmp2+tmp1*tmp14*num_b3+tmp&
                     &11*b*g3+tmp1*a*det_g*g3+tmp1*b*det_g*g3+tmp1*c*det_g*g3-(tmp1*tmp10*num_b3+tmp1&
                     &3*num_b3+tmp11*num_b3))+4.0_ki*(tmp7*num_b3+tmp1*tmp8*num_b3+tmp4*f*g3+tmp6*f*g&
                     &3+c*c*f*g3-(a*d*e*g3+tmp9*e*g3+tmp10*e*g3+tmp1*tmp9*num_b3+tmp5*tmp6+tmp4*tmp5+&
                     &tmp3*num_b3+tmp3*tmp3))+tmp11*tmp3+tmp12*tmp3+e*g3*num_b3-d*g3*num_b3)/tmp15-(4&
                     &.0_ki*(2.0_ki*(tmp7-tmp3)+tmp1*det_g-tmp13)*(2.0_ki*(tmp8-tmp9)+tmp14-tmp10)/(t&
                     &mp15*den3)*g3*num_b3+(2.0_ki*(tmp1*a+tmp1*b+tmp1*c)+d+e+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = z*z
                     tmp2 = a*e
                     tmp3 = b*d
                     tmp4 = c*d
                     tmp5 = c*e
                     tmp6 = c*f
                     tmp7 = d*e
                     tmp8 = c*g1
                     tmp9 = e*g1
                     tmp10 = tmp8*z
                     fg = (((2.0_ki*(tmp2-tmp3)+tmp5-tmp4)/(den3*den3)*tmp1*num_b3-(a+b+c)*tmp1/den3)&
                     &*g3+(1.0_ki/2.0_ki*d/den2*z-(2.0_ki*(tmp6-c*g2)-tmp7)/(den2*den2)*num_b2*z)*g2+&
                     &(3.0_ki/2.0_ki*tmp10+1.0_ki/2.0_ki*tmp9*z+1.0_ki/2.0_ki*num_b1*z+b*g1*z-(1.0_ki&
                     &/2.0_ki*tmp1*num_b1+1.0_ki/2.0_ki*tmp9+tmp8))/den1+(2.0_ki*(tmp8+tmp1*det_g+tmp&
                     &2*z+tmp6*z-(det_g*z+tmp3*z+tmp6+tmp10))+tmp7+tmp5*z-(tmp7*z+tmp4*z))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               end select
               !
            case(3)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*g1
                     tmp2 = b*d
                     tmp3 = tmp2*num_b1
                     tmp4 = b*f
                     tmp5 = b*b
                     tmp6 = b*c
                     tmp7 = c*e
                     tmp8 = tmp7*num_b1
                     tmp9 = c*f
                     tmp10 = det_g*num_b1
                     tmp11 = e*e
                     tmp12 = b*det_g*g1
                     tmp13 = c*g1
                     tmp14 = e*g1
                     tmp15 = d*e
                     tmp16 = den1*den1
                     fg = (((c+e+2.0_ki*b*z)/den1+4.0_ki*(2.0_ki*tmp2+det_g-tmp7)*(4.0_ki*(tmp4-tmp1)&
                     &+2.0_ki*(tmp9-tmp13)+det_g*z-(det_g+tmp15+tmp11))/(tmp16*den1)*g1*num_b1-(8.0_k&
                     &i*(tmp5*f*g1-(tmp5*d*g1+tmp4*num_b1+tmp3+tmp1*tmp1))+4.0_ki*(tmp8+tmp1*num_b1+t&
                     &mp3*z+tmp6*e*g1+tmp6*f*g1+a*e*num_b1-(tmp6*d*g1+tmp2*e*g1+tmp9*num_b1+tmp6*g1*g&
                     &1))+2.0_ki*(tmp10*z+tmp11*num_b1+tmp12*z+tmp13*num_b1+tmp14*c*c+tmp15*num_b1-(c&
                     &*det_g*g1+c*d*num_b1+tmp8*z+tmp1*tmp11+3.0_ki*tmp12+tmp10))+tmp11*tmp13-(det_g*&
                     &e*g1+tmp14*num_b1))/tmp16)*g1)
                     !
                  case(2)
                     !
                     tmp1 = b*f
                     tmp2 = tmp1*num_b2
                     tmp3 = b*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = c*f
                     tmp6 = tmp5*num_b2
                     tmp7 = c*g2
                     tmp8 = tmp7*num_b2
                     tmp9 = e*e
                     tmp10 = b*c
                     tmp11 = g2*g2
                     tmp12 = b*d
                     tmp13 = tmp12*num_b2
                     tmp14 = d*e
                     tmp15 = tmp14*num_b2
                     tmp16 = det_g*num_b2
                     tmp17 = z*z
                     tmp18 = tmp16*z
                     tmp19 = c*d
                     tmp20 = c*e
                     tmp21 = tmp20*num_b2
                     tmp22 = num_b2*z
                     tmp23 = e*g2
                     tmp24 = den2*den2
                     tmp25 = d*g2
                     fg = ((4.0_ki*(2.0_ki*tmp12+det_g*z-tmp20)*(4.0_ki*(tmp1-tmp3)+2.0_ki*(tmp5-tmp7&
                     &)-(tmp17*det_g+tmp9+tmp14))/(tmp24*den2)*g2*num_b2-(16.0_ki*(tmp4-tmp2)+6.0_ki*&
                     &(-(tmp8*z+2.0_ki*tmp4*z))+8.0_ki*(tmp8+tmp2*z-tmp6)+4.0_ki*(tmp15+tmp16*tmp17+t&
                     &mp6*z+tmp9*num_b2+tmp10*f*g2+tmp12*e*g2+a*e*num_b2-(tmp17*tmp18+tmp13*tmp17+tmp&
                     &10*tmp11+tmp13))+2.0_ki*(tmp21+tmp17*tmp21+tmp3*d*d+c*c*f*g2+det_g*e*g2*z-(tmp1&
                     &9*e*g2+tmp22*tmp9+tmp19*num_b2+tmp15*z+tmp7*tmp7+tmp18))+d*det_g*g2*z-(tmp17*c*&
                     &det_g*g2+3.0_ki*tmp7*tmp9+tmp23*num_b2))/tmp24)*g2+(2.0_ki*(tmp17*num_b2+tmp23*&
                     &z-(a*g2+tmp25+2.0_ki*tmp23))+1.0_ki/2.0_ki*tmp11+tmp25*z-(tmp17*tmp7+tmp17*tmp2&
                     &2+tmp7+tmp22))/den2)
                     !
                  case(3)
                     !
                     tmp1 = b*f
                     tmp2 = tmp1*num_b3
                     tmp3 = b*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = b*c
                     tmp6 = f*g3
                     tmp7 = c*f
                     tmp8 = tmp7*num_b3
                     tmp9 = c*g3
                     tmp10 = tmp9*num_b3
                     tmp11 = e*e
                     tmp12 = b*d
                     tmp13 = d*e
                     tmp14 = tmp13*num_b3
                     tmp15 = c*d
                     tmp16 = c*e
                     tmp17 = den3*den3
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(-(e*z+d*z))+4.0_ki*(a+b+c+d+e)-g3)/den3+4.0_ki*(2.&
                     &0_ki*(2.0_ki*(tmp1-tmp3)+tmp7-tmp9)-(tmp13+tmp11))*(2.0_ki*tmp12+det_g*z-tmp16)&
                     &/(tmp17*den3)*g3*num_b3-(4.0_ki*(tmp14+tmp11*num_b3+tmp8*z+a*e*num_b3-tmp12*num&
                     &_b3)+8.0_ki*(2.0_ki*(tmp4-tmp2)+tmp10+tmp3*tmp3+tmp2*z+tmp5*g3*g3-(tmp5*f*g3+tm&
                     &p6*b*b+tmp8))+2.0_ki*(3.0_ki*(-(2.0_ki*tmp4*z+tmp10*z))+tmp9*tmp9+tmp11*tmp3+tm&
                     &p11*tmp9+tmp16*num_b3+tmp15*e*g3-(tmp11*num_b3*z+tmp6*c*c+tmp3*d*d+tmp15*num_b3&
                     &+tmp14*z))-(det_g*e*g3*z+d*det_g*g3*z+e*g3*num_b3))/tmp17)*g3)
                     !
                  case(4)
                     !
                     tmp1 = b*g1
                     tmp2 = c*g1
                     tmp3 = b*f
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = tmp4*z
                     tmp7 = e*e
                     tmp8 = c*e
                     tmp9 = d*e
                     tmp10 = tmp8*z
                     tmp11 = c*g2
                     tmp12 = z*z
                     fg = ((-((4.0_ki*(tmp3-b*g3)+2.0_ki*(tmp5-c*g3)-(tmp9+tmp7))/(den3*den3)*num_b3*&
                     &z+1.0_ki/2.0_ki*(d+e)/den3*z))*g3+(tmp1*z+1.0_ki/2.0_ki*num_b1*z-(1.0_ki/2.0_ki&
                     &*e*g1+1.0_ki/2.0_ki*num_b1+tmp2+3.0_ki*tmp1))/den1+(2.0_ki*(2.0_ki*(b*g2-tmp3)+&
                     &tmp11+tmp6+tmp12*det_g-tmp5)+tmp7+tmp9-tmp10)/(den2*den2)*g2*num_b2*z-((2.0_ki*&
                     &(2.0_ki*(tmp3-tmp1)+tmp5+tmp6+det_g*z-(det_g+tmp4+tmp2))+tmp8-(tmp9+tmp7+tmp10)&
                     &)/(den1*den1)*g1*num_b1+1.0_ki/2.0_ki*(2.0_ki*(num_b2-(2.0_ki*num_b2*z+e*g2))+t&
                     &mp11*z+3.0_ki*tmp12*num_b2-d*g2)/den2*z))
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = b*c
                     tmp3 = a*e
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = c*g1
                     tmp7 = e*e
                     tmp8 = c*det_g*g1
                     tmp9 = tmp6*num_b1
                     tmp10 = d*g1
                     tmp11 = c*c
                     tmp12 = e*g1
                     tmp13 = det_g*num_b1
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = d*e
                     tmp17 = det_g*e*g1
                     tmp18 = den1*den1
                     tmp19 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(3.0_ki*c*z+e*z)+4.0_ki*(2.0_ki*a+d+b*z)+g1)/den1-(&
                     &4.0_ki*(2.0_ki*(tmp5-tmp6)+tmp19-tmp16)*(2.0_ki*(tmp3-tmp4)+tmp15+tmp19-(det_g+&
                     &tmp14))/(tmp18*den1)*g1*num_b1+(3.0_ki*(tmp8*z+2.0_ki*tmp9*z)+4.0_ki*(2.0_ki*(t&
                     &mp1*num_b1-a*f*num_b1)+tmp2*f*g1+tmp4*num_b1*z+a*c*e*g1-(tmp5*num_b1*z+tmp4*e*g&
                     &1+tmp3*num_b1*z+tmp2*d*g1+tmp2*g1*g1))+2.0_ki*(tmp1*tmp7+tmp11*tmp12+tmp13*z+d*&
                     &d*num_b1+tmp11*f*g1+tmp14*num_b1*z+tmp16*num_b1*z+b*det_g*g1*z-(tmp15*num_b1*z+&
                     &tmp14*e*g1+tmp13*z*z+tmp10*tmp11+tmp9+tmp8+tmp6*tmp6))+tmp10*num_b1+tmp17*z+tmp&
                     &6*tmp7-(tmp12*num_b1+tmp17))/tmp18))*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = c*f
                     tmp3 = tmp2*num_b2
                     tmp4 = c*g2
                     tmp5 = tmp4*num_b2
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = tmp7*num_b2
                     tmp9 = c*d
                     tmp10 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z-d)+g2)/den2-(4.0_ki*(2.0_ki*(tmp2-tmp4)-tmp7)*&
                     &(2.0_ki*(a*e-b*d)+c*e-(det_g*z+tmp9))/(tmp10*den2)*g2*num_b2+(4.0_ki*(2.0_ki*(a&
                     &*g2*num_b2-a*f*num_b2)+tmp5+tmp3*z+tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2+tmp3))+2.0_ki&
                     &*(tmp8+tmp6*num_b2+tmp6*b*g2+c*c*f*g2-(tmp9*e*g2+tmp8*z+3.0_ki*tmp5*z+tmp4*tmp4&
                     &))+tmp4*tmp6+d*g2*num_b2+d*det_g*g2*z-e*g2*num_b2)/tmp10))*g2)
                     !
                  case(3)
                     !
                     tmp1 = z*z
                     tmp2 = a*g3
                     tmp3 = c*g3
                     tmp4 = a*c
                     tmp5 = g3*g3
                     tmp6 = b*c
                     tmp7 = c*f
                     tmp8 = a*e
                     tmp9 = b*d
                     tmp10 = c*d
                     tmp11 = d*d
                     tmp12 = e*e
                     tmp13 = d*e
                     tmp14 = c*e
                     tmp15 = den3*den3
                     fg = (((8.0_ki*(a*f*num_b3-tmp2*num_b3)+2.0_ki*(tmp12*tmp2+tmp1*tmp14*num_b3+tmp&
                     &11*b*g3+tmp1*a*det_g*g3+tmp1*b*det_g*g3+tmp1*c*det_g*g3-(tmp1*tmp10*num_b3+tmp1&
                     &3*num_b3+tmp11*num_b3))+4.0_ki*(tmp7*num_b3+tmp1*tmp8*num_b3+tmp4*f*g3+tmp6*f*g&
                     &3+c*c*f*g3-(a*d*e*g3+tmp9*e*g3+tmp10*e*g3+tmp1*tmp9*num_b3+tmp5*tmp6+tmp4*tmp5+&
                     &tmp3*num_b3+tmp3*tmp3))+tmp11*tmp3+tmp12*tmp3+e*g3*num_b3-d*g3*num_b3)/tmp15-(4&
                     &.0_ki*(2.0_ki*(tmp7-tmp3)+tmp1*det_g-tmp13)*(2.0_ki*(tmp8-tmp9)+tmp14-tmp10)/(t&
                     &mp15*den3)*g3*num_b3+(2.0_ki*(tmp1*a+tmp1*b+tmp1*c)+d+e+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = z*z
                     tmp2 = a*e
                     tmp3 = b*d
                     tmp4 = c*d
                     tmp5 = c*e
                     tmp6 = c*f
                     tmp7 = d*e
                     tmp8 = c*g1
                     tmp9 = e*g1
                     tmp10 = tmp8*z
                     fg = (((2.0_ki*(tmp2-tmp3)+tmp5-tmp4)/(den3*den3)*tmp1*num_b3-(a+b+c)*tmp1/den3)&
                     &*g3+(1.0_ki/2.0_ki*d/den2*z-(2.0_ki*(tmp6-c*g2)-tmp7)/(den2*den2)*num_b2*z)*g2+&
                     &(3.0_ki/2.0_ki*tmp10+1.0_ki/2.0_ki*tmp9*z+1.0_ki/2.0_ki*num_b1*z+b*g1*z-(1.0_ki&
                     &/2.0_ki*tmp1*num_b1+1.0_ki/2.0_ki*tmp9+tmp8))/den1+(2.0_ki*(tmp8+tmp1*det_g+tmp&
                     &2*z+tmp6*z-(det_g*z+tmp3*z+tmp6+tmp10))+tmp7+tmp5*z-(tmp7*z+tmp4*z))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*g1
                     tmp2 = a*g1
                     tmp3 = b*c
                     tmp4 = d*g1
                     tmp5 = b*b
                     tmp6 = tmp4*tmp5
                     tmp7 = f*g1
                     tmp8 = tmp3*d
                     tmp9 = tmp8*g1
                     tmp10 = b*d
                     tmp11 = tmp10*num_b1
                     tmp12 = c*det_g
                     tmp13 = tmp12*g1
                     tmp14 = c*g1
                     tmp15 = tmp14*num_b1
                     tmp16 = e*g1
                     tmp17 = c*c
                     tmp18 = tmp16*tmp17
                     tmp19 = det_g*num_b1
                     tmp20 = a*c*e
                     tmp21 = tmp20*g1
                     tmp22 = tmp3*e
                     tmp23 = tmp22*g1
                     tmp24 = b*det_g
                     tmp25 = tmp24*g1
                     tmp26 = c*e
                     tmp27 = tmp26*num_b1
                     tmp28 = tmp13*z
                     tmp29 = d*d
                     tmp30 = e*e
                     tmp31 = tmp17*tmp4
                     tmp32 = det_g*e*g1
                     tmp33 = c*d
                     tmp34 = den1*den1
                     tmp35 = d*f
                     tmp36 = tmp12*f
                     tmp37 = tmp10*det_g
                     tmp38 = det_g*det_g
                     tmp39 = d*det_g*e
                     tmp12 = tmp12*e
                     tmp40 = det_g*z
                     fg = (((3.0_ki*(2.0_ki*(tmp15*z-tmp25)+tmp28)+4.0_ki*(2.0_ki*(tmp11*z+tmp2*num_b&
                     &1+tmp5*tmp7+tmp3*f*g1-(a*f*num_b1+tmp10*e*g1+tmp3*g1*g1+tmp9+tmp6+tmp1*tmp1))+t&
                     &mp18+tmp21+tmp23+tmp19*z+tmp25*z-(c*f*num_b1*z+a*e*num_b1*z+tmp27*z+tmp1*num_b1&
                     &+tmp15+tmp13+tmp11))+2.0_ki*(tmp27+tmp14*tmp30+tmp17*tmp7+tmp2*tmp30+tmp29*num_&
                     &b1+tmp33*num_b1*z+d*e*num_b1*z-(tmp33*e*g1+tmp19*z*z+tmp16*num_b1+tmp1*tmp30+tm&
                     &p32+tmp31+tmp19+tmp14*tmp14))+tmp32*z+tmp4*num_b1)/tmp34+4.0_ki*(2.0_ki*(tmp31+&
                     &tmp37+tmp39+tmp10*tmp30+tmp12*z+tmp36*z+a*det_g*e*z-(tmp30*a*d+tmp38*z+tmp30*tm&
                     &p33+tmp17*tmp35+tmp28))+4.0_ki*(2.0_ki*(tmp6+tmp9-(tmp8*f+tmp35*tmp5))+tmp13+tm&
                     &p25+tmp20*f+tmp22*f+tmp17*e*f+tmp29*b*e-(tmp37*z+tmp24*f+tmp36+tmp23+tmp21+tmp1&
                     &8))+tmp38+tmp40*tmp40+tmp26*tmp29+tmp30*det_g-(tmp33*det_g*z+tmp39*z+tmp26*tmp3&
                     &0+tmp12))/(tmp34*den1)*g1*num_b1-1.0_ki/2.0_ki*(8.0_ki*(a+b*z)+4.0_ki*(d-b)+2.0&
                     &_ki*(3.0_ki*c*z+e*z)+g1)/den1)*g1)
                     !
                  case(2)
                     !
                     tmp1 = b*g2
                     tmp2 = tmp1*num_b2
                     tmp3 = c*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = b*f*num_b2
                     tmp6 = c*f*num_b2
                     tmp7 = a*c
                     tmp8 = g2*g2
                     tmp9 = b*c
                     tmp10 = d*d
                     tmp11 = c*c
                     tmp12 = d*e*num_b2
                     tmp13 = det_g*num_b2
                     tmp14 = z*z
                     tmp15 = a*d
                     tmp16 = b*d
                     tmp17 = c*d
                     tmp18 = e*e
                     tmp19 = e*g2
                     tmp20 = num_b2*z
                     tmp21 = c*e
                     tmp22 = d*det_g
                     tmp23 = d*g2
                     tmp24 = c*det_g
                     tmp25 = tmp24*g2
                     tmp26 = den2*den2
                     tmp27 = d*f
                     tmp28 = b*b
                     tmp29 = tmp9*d
                     tmp30 = tmp7*e
                     tmp31 = tmp9*e
                     tmp32 = b*det_g
                     tmp33 = det_g*z
                     fg = (((12.0_ki*(-(tmp4*z+tmp2*z))+8.0_ki*(tmp2+tmp4+tmp5*z+tmp6*z+a*g2*num_b2-(&
                     &a*f*num_b2+tmp6+tmp5))+2.0_ki*(tmp10*num_b2+tmp13*tmp14+tmp18*num_b2+tmp14*tmp2&
                     &1*num_b2+tmp22*g2*z+det_g*e*g2*z-(tmp19*num_b2+tmp18*tmp20))+4.0_ki*(tmp12+tmp1&
                     &*tmp10+tmp11*f*g2+tmp16*e*g2+tmp7*f*g2+tmp9*f*g2-(tmp17*e*g2+tmp15*e*g2+tmp14*t&
                     &mp16*num_b2+tmp13*tmp14*z+tmp8*tmp9+tmp7*tmp8+tmp12*z+tmp3*tmp3))+tmp10*tmp3+tm&
                     &p23*num_b2-(3.0_ki*tmp18*tmp3+tmp14*tmp25))/tmp26+4.0_ki*(8.0_ki*(tmp23*tmp28+t&
                     &mp29*g2-(tmp29*f+tmp27*tmp28))+2.0_ki*(tmp11*tmp23+tmp16*tmp18+tmp14*tmp16*det_&
                     &g+tmp22*e*z-(tmp17*tmp18+tmp15*tmp18+tmp11*tmp27))+4.0_ki*(tmp25*z+tmp30*f+tmp3&
                     &1*f+tmp10*b*e+tmp11*e*f+tmp32*g2*z-(tmp32*f*z+tmp24*f*z+tmp31*g2+tmp30*g2+tmp11&
                     &*tmp19))+tmp10*tmp21+tmp18*tmp33+tmp33*tmp33*z-(tmp14*tmp24*e+tmp18*tmp21))/(tm&
                     &p26*den2)*g2*num_b2)*g2+(2.0_ki*(tmp19+tmp23-(tmp23*z+tmp19*z))+tmp14*tmp20+tmp&
                     &14*tmp3-(tmp14*num_b2+tmp8))/den2)
                     !
                  case(3)
                     !
                     tmp1 = b*c
                     tmp2 = g3*g3
                     tmp3 = b*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = c*g3
                     tmp6 = f*g3
                     tmp7 = c*c
                     tmp8 = c*d
                     tmp9 = tmp5*num_b3
                     tmp10 = a*g3
                     tmp11 = b*f*num_b3
                     tmp12 = b*b
                     tmp13 = c*f*num_b3
                     tmp14 = a*c
                     tmp15 = d*d
                     tmp16 = d*e*num_b3
                     tmp17 = a*d
                     tmp18 = z*z
                     tmp19 = b*d
                     tmp20 = e*e
                     tmp21 = e*g3
                     tmp22 = a*det_g
                     tmp23 = b*det_g
                     tmp24 = tmp23*g3
                     tmp25 = c*det_g
                     tmp26 = tmp25*g3
                     tmp27 = c*e
                     tmp28 = d*g3
                     tmp29 = d*det_g
                     tmp30 = den3*den3
                     tmp31 = d*f
                     tmp32 = tmp1*d
                     tmp33 = tmp14*e
                     tmp34 = tmp1*e
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+4.0_ki*(tmp18*a+tmp18*b+tmp18*c)+3.0_ki*g&
                     &3)/den3+4.0_ki*(2.0_ki*(tmp19*tmp20+tmp26*z+tmp28*tmp7+tmp18*tmp22*e-(tmp25*f*z&
                     &+tmp18*tmp19*det_g+tmp31*tmp7+tmp20*tmp8+tmp17*tmp20))+4.0_ki*(2.0_ki*(tmp12*tm&
                     &p28+tmp32*g3-(tmp32*f+tmp12*tmp31))+tmp24*z+tmp33*f+tmp34*f+tmp15*b*e+tmp7*e*f-&
                     &(tmp23*f*z+tmp34*g3+tmp33*g3+tmp21*tmp7))+tmp15*tmp27+tmp18*tmp25*e+tmp20*det_g&
                     &*z+tmp29*e*z-(tmp18*tmp8*det_g+tmp20*tmp27))/(tmp30*den3)*g3*num_b3-(6.0_ki*(2.&
                     &0_ki*(tmp4*z+tmp1*f*g3-tmp1*tmp2)+tmp6*tmp7+tmp9*z-(tmp8*e*g3+tmp5*tmp5))+4.0_k&
                     &i*(2.0_ki*(tmp11+tmp13+tmp12*tmp6+a*f*num_b3-(tmp11*z+tmp10*num_b3+tmp9+tmp4+tm&
                     &p3*tmp3))+tmp15*tmp3+tmp14*f*g3+tmp18*a*e*num_b3-(tmp19*e*g3+tmp18*tmp19*num_b3&
                     &+tmp17*e*g3+tmp14*tmp2+tmp13*z+tmp16))+2.0_ki*(tmp10*tmp20+tmp16*z+tmp18*tmp24+&
                     &tmp18*tmp26+tmp21*num_b3+tmp18*tmp22*g3+tmp18*tmp27*num_b3+tmp20*num_b3*z-(tmp1&
                     &8*tmp8*num_b3+tmp20*num_b3+tmp20*tmp3+tmp15*num_b3))+tmp15*tmp5+tmp29*g3*z+det_&
                     &g*e*g3*z-(tmp28*num_b3+tmp20*tmp5))/tmp30)*g3)
                     !
                  case(4)
                     !
                     tmp1 = c*g1
                     tmp2 = b*g1
                     tmp3 = e*g1
                     tmp4 = tmp1*z
                     tmp5 = z*z
                     tmp6 = b*f
                     tmp7 = c*f
                     tmp8 = b*d
                     tmp9 = tmp8*z
                     tmp10 = d*e
                     tmp11 = tmp5*det_g
                     tmp12 = a*e*z
                     tmp13 = c*e
                     tmp14 = tmp13*z
                     tmp15 = e*e
                     tmp16 = c*d*z
                     tmp17 = c*g2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(a*z+b*z+c*z)+d+e)/den3*z-(4.0_ki*(b*g3-tmp6)+2.0_k&
                     &i*(tmp12+c*g3-(tmp9+tmp7))+tmp10+tmp14+tmp15-tmp16)/(den3*den3)*num_b3*z)*g3+(2&
                     &.0_ki*(tmp1-tmp2*z)+3.0_ki/2.0_ki*(2.0_ki*tmp2-tmp4)+tmp3+1.0_ki/2.0_ki*num_b1+&
                     &1.0_ki/2.0_ki*tmp5*num_b1-(num_b1*z+1.0_ki/2.0_ki*tmp3*z))/den1+1.0_ki/2.0_ki*(&
                     &2.0_ki*(-(num_b2*z+e*g2+d*g2))+tmp17*z+3.0_ki*tmp5*num_b2)/den2*z-((2.0_ki*(2.0&
                     &_ki*(tmp17+b*g2-(tmp7+tmp6))+tmp10+tmp11+tmp9)+tmp15-tmp14)/(den2*den2)*g2*num_&
                     &b2*z+(2.0_ki*(2.0_ki*(tmp1+tmp2-(det_g*z+tmp9+tmp7+tmp6))+tmp10+tmp11+tmp12+tmp&
                     &14+tmp8+det_g+tmp7*z-tmp4)+tmp15-(tmp10*z+tmp16+tmp13))/(den1*den1)*g1*num_b1))
                     !
                  end select
                  !
               end select
               !
            end select
            !
         case(2)
            !
            select case(par2)
            !
            case(1)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*c
                     tmp2 = c*f
                     tmp3 = b*d
                     tmp4 = c*g1
                     tmp5 = e*g1
                     tmp6 = d*e
                     tmp7 = c*e
                     tmp8 = den1*den1
                     fg = (((4.0_ki*(tmp1*f*g1+tmp3*num_b1*z+a*e*num_b1-(tmp3*e*g1+tmp1*d*g1+tmp2*num&
                     &_b1+tmp1*g1*g1))+2.0_ki*(tmp4*num_b1+tmp5*c*c+tmp6*num_b1+det_g*num_b1*z+b*det_&
                     &g*g1*z-(c*det_g*g1+c*d*num_b1+tmp7*num_b1*z))+tmp4*e*e-(det_g*e*g1+tmp5*num_b1)&
                     &)/tmp8-(4.0_ki*(2.0_ki*(tmp2-tmp4)+det_g*z-tmp6)*(2.0_ki*tmp3+det_g-tmp7)/(tmp8&
                     &*den1)*g1*num_b1+(c+e+2.0_ki*b*z)/den1))*g1)
                     !
                  case(2)
                     !
                     tmp1 = c*f
                     tmp2 = tmp1*num_b2
                     tmp3 = c*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = d*e
                     tmp6 = tmp5*num_b2
                     tmp7 = c*d
                     tmp8 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(4.0_ki*(a+d)-(2.0_ki*d*z+g2))/den2+(4.0_ki*(2.0_ki*(tmp4-t&
                     &mp2)+tmp6+tmp2*z+a*e*num_b2)+2.0_ki*(b*d*d*g2+c*c*f*g2-(tmp7*e*g2+tmp7*num_b2+t&
                     &mp6*z+3.0_ki*tmp4*z+tmp3*tmp3))+d*det_g*g2*z-e*g2*num_b2)/tmp8-4.0_ki*(2.0_ki*(&
                     &tmp1-tmp3)-tmp5)*(2.0_ki*b*d+det_g*z-c*e)/(tmp8*den2)*g2*num_b2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = c*f
                     tmp2 = tmp1*num_b3
                     tmp3 = c*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = b*c
                     tmp6 = g3*g3
                     tmp7 = d*e
                     tmp8 = tmp7*num_b3
                     tmp9 = det_g*num_b3
                     tmp10 = z*z
                     tmp11 = tmp9*z
                     tmp12 = b*d
                     tmp13 = b*g3
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = e*g3
                     tmp17 = den3*den3
                     tmp18 = d*g3
                     tmp19 = num_b3*z
                     fg = (((8.0_ki*(tmp4-tmp2)+4.0_ki*(tmp8+tmp10*tmp11+tmp2*z+tmp5*tmp6+tmp10*tmp12&
                     &*num_b3+tmp12*e*g3+a*e*num_b3-(tmp5*f*g3+tmp10*tmp9))+2.0_ki*(tmp11+tmp3*tmp3+t&
                     &mp14*e*g3-(tmp10*b*det_g*g3+c*c*f*g3+tmp10*tmp15*num_b3+tmp8*z+3.0_ki*tmp4*z+tm&
                     &p14*num_b3+tmp13*d*d))+det_g*e*g3*z-(d*det_g*g3*z+tmp10*c*det_g*g3+tmp3*e*e+tmp&
                     &16*num_b3))/tmp17-4.0_ki*(2.0_ki*(tmp1-tmp3)+tmp10*det_g-tmp7)*(2.0_ki*tmp12+de&
                     &t_g*z-tmp15)/(tmp17*den3)*g3*num_b3)*g3+(2.0_ki*(tmp16+tmp10*tmp13+tmp10*num_b3&
                     &-(a*g3+tmp18))+1.0_ki/2.0_ki*tmp6+tmp10*tmp3+tmp18*z-(tmp16*z+tmp10*tmp19+tmp3+&
                     &tmp19))/den3)
                     !
                  case(4)
                     !
                     tmp1 = c*g1
                     tmp2 = c*f
                     tmp3 = b*d*z
                     tmp4 = d*e
                     tmp5 = c*e*z
                     tmp6 = c*g3
                     tmp7 = z*z
                     fg = (((2.0_ki*(tmp2-c*g2)-tmp4)/(den2*den2)*num_b2*z-1.0_ki/2.0_ki*d/den2*z)*g2&
                     &+(tmp1+1.0_ki/2.0_ki*e*g1-(b*g1*z+1.0_ki/2.0_ki*num_b1*z))/den1+1.0_ki/2.0_ki*(&
                     &2.0_ki*(2.0_ki*num_b3*z+b*g3*z-num_b3)+tmp6*z+d*g3-(e*g3+3.0_ki*tmp7*num_b3))/d&
                     &en3*z+(2.0_ki*(tmp2+tmp3+det_g*z-tmp1)-(tmp5+tmp4))/(den1*den1)*g1*num_b1+(2.0_&
                     &ki*(tmp2+tmp3+tmp7*det_g-tmp6)-(tmp5+tmp4))/(den3*den3)*g3*num_b3*z)
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = a*e
                     tmp3 = c*f
                     tmp4 = c*g1
                     tmp5 = d*g1
                     tmp6 = c*c
                     tmp7 = c*d
                     tmp8 = d*e
                     tmp9 = den1*den1
                     tmp10 = det_g*z
                     fg = (((3.0_ki*(2.0_ki*tmp4*num_b1*z+c*det_g*g1*z)+4.0_ki*(2.0_ki*(tmp1*num_b1-a&
                     &*f*num_b1)+a*c*e*g1-(tmp3*num_b1*z+tmp2*num_b1*z))+2.0_ki*(tmp1*e*e+d*d*num_b1+&
                     &tmp6*f*g1+tmp7*num_b1*z+tmp8*num_b1*z-(det_g*num_b1*z*z+tmp7*e*g1+tmp5*tmp6+tmp&
                     &4*tmp4))+tmp5*num_b1+det_g*e*g1*z)/tmp9+4.0_ki*(2.0_ki*(tmp3-tmp4)+tmp10-tmp8)*&
                     &(tmp10+2.0_ki*tmp2-tmp7)/(tmp9*den1)*g1*num_b1-1.0_ki/2.0_ki*(4.0_ki*(2.0_ki*a+&
                     &d)+2.0_ki*(3.0_ki*c*z+e*z)+g1)/den1)*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = d*d
                     tmp3 = c*g2
                     tmp4 = den2*den2
                     fg = (((8.0_ki*(a*g2*num_b2-a*f*num_b2)+4.0_ki*(tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2))&
                     &+tmp2*tmp3+2.0_ki*tmp2*num_b2+d*g2*num_b2)/tmp4+4.0_ki*(2.0_ki*(c*f-tmp3)-d*e)*&
                     &(2.0_ki*a*e-c*d)/(tmp4*den2)*g2*num_b2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*g3
                     tmp2 = a*c
                     tmp3 = g3*g3
                     tmp4 = det_g*num_b3
                     tmp5 = z*z
                     tmp6 = a*e
                     tmp7 = c*f
                     tmp8 = c*g3
                     tmp9 = d*d
                     tmp10 = c*d
                     tmp11 = d*e
                     tmp12 = d*g3
                     tmp13 = den3*den3
                     fg = ((4.0_ki*(2.0_ki*(tmp7-tmp8)+tmp5*det_g-tmp11)*(2.0_ki*tmp6+det_g*z-tmp10)/&
                     &(tmp13*den3)*g3*num_b3-(8.0_ki*(a*f*num_b3-tmp1*num_b3)+4.0_ki*(tmp2*f*g3+tmp4*&
                     &tmp5*z+tmp5*tmp6*num_b3+tmp7*num_b3*z-(a*d*e*g3+tmp2*tmp3))+2.0_ki*(tmp1*e*e+c*&
                     &c*f*g3+tmp5*a*det_g*g3-(3.0_ki*tmp8*num_b3*z+tmp11*num_b3*z+tmp10*e*g3+tmp10*tm&
                     &p5*num_b3+tmp9*num_b3+tmp4*tmp5+tmp8*tmp8))+tmp8*tmp9+tmp5*c*det_g*g3+det_g*e*g&
                     &3*z-(d*det_g*g3*z+tmp12*num_b3))/tmp13)*g3+(2.0_ki*(tmp12+tmp1*tmp5)+1.0_ki/2.0&
                     &_ki*tmp3+tmp5*tmp8+tmp5*num_b3*z+e*g3*z-(tmp5*num_b3+tmp12*z))/den3)
                     !
                  case(4)
                     !
                     tmp1 = a*e
                     tmp2 = c*f
                     tmp3 = c*g1
                     tmp4 = c*d
                     tmp5 = d*e
                     tmp6 = c*g3
                     tmp7 = z*z
                     fg = (1.0_ki/2.0_ki*(2.0_ki*(a*g3*z-num_b3*z)+tmp6*z+3.0_ki*tmp7*num_b3+e*g3-d*g&
                     &3)/den3*z-((2.0_ki*(tmp2+tmp1*z+tmp7*det_g-tmp6)-(tmp4*z+tmp5))/(den3*den3)*g3*&
                     &num_b3*z+(2.0_ki*(tmp1+tmp2+det_g*z-tmp3)-(tmp5+tmp4))/(den1*den1)*g1*num_b1*z+&
                     &1.0_ki/2.0_ki*(3.0_ki*tmp3+e*g1-num_b1*z)/den1*z))
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = b*c
                     tmp3 = a*e
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = c*g1
                     tmp7 = e*e
                     tmp8 = c*det_g*g1
                     tmp9 = tmp6*num_b1
                     tmp10 = d*g1
                     tmp11 = c*c
                     tmp12 = e*g1
                     tmp13 = det_g*num_b1
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = d*e
                     tmp17 = det_g*e*g1
                     tmp18 = den1*den1
                     tmp19 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(3.0_ki*c*z+e*z)+4.0_ki*(2.0_ki*a+d+b*z)+g1)/den1-(&
                     &4.0_ki*(2.0_ki*(tmp5-tmp6)+tmp19-tmp16)*(2.0_ki*(tmp3-tmp4)+tmp15+tmp19-(det_g+&
                     &tmp14))/(tmp18*den1)*g1*num_b1+(3.0_ki*(tmp8*z+2.0_ki*tmp9*z)+4.0_ki*(2.0_ki*(t&
                     &mp1*num_b1-a*f*num_b1)+tmp2*f*g1+tmp4*num_b1*z+a*c*e*g1-(tmp5*num_b1*z+tmp4*e*g&
                     &1+tmp3*num_b1*z+tmp2*d*g1+tmp2*g1*g1))+2.0_ki*(tmp1*tmp7+tmp11*tmp12+tmp13*z+d*&
                     &d*num_b1+tmp11*f*g1+tmp14*num_b1*z+tmp16*num_b1*z+b*det_g*g1*z-(tmp15*num_b1*z+&
                     &tmp14*e*g1+tmp13*z*z+tmp10*tmp11+tmp9+tmp8+tmp6*tmp6))+tmp10*num_b1+tmp17*z+tmp&
                     &6*tmp7-(tmp12*num_b1+tmp17))/tmp18))*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = c*f
                     tmp3 = tmp2*num_b2
                     tmp4 = c*g2
                     tmp5 = tmp4*num_b2
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = tmp7*num_b2
                     tmp9 = c*d
                     tmp10 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z-d)+g2)/den2-(4.0_ki*(2.0_ki*(tmp2-tmp4)-tmp7)*&
                     &(2.0_ki*(a*e-b*d)+c*e-(det_g*z+tmp9))/(tmp10*den2)*g2*num_b2+(4.0_ki*(2.0_ki*(a&
                     &*g2*num_b2-a*f*num_b2)+tmp5+tmp3*z+tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2+tmp3))+2.0_ki&
                     &*(tmp8+tmp6*num_b2+tmp6*b*g2+c*c*f*g2-(tmp9*e*g2+tmp8*z+3.0_ki*tmp5*z+tmp4*tmp4&
                     &))+tmp4*tmp6+d*g2*num_b2+d*det_g*g2*z-e*g2*num_b2)/tmp10))*g2)
                     !
                  case(3)
                     !
                     tmp1 = z*z
                     tmp2 = a*g3
                     tmp3 = c*g3
                     tmp4 = a*c
                     tmp5 = g3*g3
                     tmp6 = b*c
                     tmp7 = c*f
                     tmp8 = a*e
                     tmp9 = b*d
                     tmp10 = c*d
                     tmp11 = d*d
                     tmp12 = e*e
                     tmp13 = d*e
                     tmp14 = c*e
                     tmp15 = den3*den3
                     fg = (((8.0_ki*(a*f*num_b3-tmp2*num_b3)+2.0_ki*(tmp12*tmp2+tmp1*tmp14*num_b3+tmp&
                     &11*b*g3+tmp1*a*det_g*g3+tmp1*b*det_g*g3+tmp1*c*det_g*g3-(tmp1*tmp10*num_b3+tmp1&
                     &3*num_b3+tmp11*num_b3))+4.0_ki*(tmp7*num_b3+tmp1*tmp8*num_b3+tmp4*f*g3+tmp6*f*g&
                     &3+c*c*f*g3-(a*d*e*g3+tmp9*e*g3+tmp10*e*g3+tmp1*tmp9*num_b3+tmp5*tmp6+tmp4*tmp5+&
                     &tmp3*num_b3+tmp3*tmp3))+tmp11*tmp3+tmp12*tmp3+e*g3*num_b3-d*g3*num_b3)/tmp15-(4&
                     &.0_ki*(2.0_ki*(tmp7-tmp3)+tmp1*det_g-tmp13)*(2.0_ki*(tmp8-tmp9)+tmp14-tmp10)/(t&
                     &mp15*den3)*g3*num_b3+(2.0_ki*(tmp1*a+tmp1*b+tmp1*c)+d+e+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = z*z
                     tmp2 = a*e
                     tmp3 = b*d
                     tmp4 = c*d
                     tmp5 = c*e
                     tmp6 = c*f
                     tmp7 = d*e
                     tmp8 = c*g1
                     tmp9 = e*g1
                     tmp10 = tmp8*z
                     fg = (((2.0_ki*(tmp2-tmp3)+tmp5-tmp4)/(den3*den3)*tmp1*num_b3-(a+b+c)*tmp1/den3)&
                     &*g3+(1.0_ki/2.0_ki*d/den2*z-(2.0_ki*(tmp6-c*g2)-tmp7)/(den2*den2)*num_b2*z)*g2+&
                     &(3.0_ki/2.0_ki*tmp10+1.0_ki/2.0_ki*tmp9*z+1.0_ki/2.0_ki*num_b1*z+b*g1*z-(1.0_ki&
                     &/2.0_ki*tmp1*num_b1+1.0_ki/2.0_ki*tmp9+tmp8))/den1+(2.0_ki*(tmp8+tmp1*det_g+tmp&
                     &2*z+tmp6*z-(det_g*z+tmp3*z+tmp6+tmp10))+tmp7+tmp5*z-(tmp7*z+tmp4*z))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               end select
               !
            case(2)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = a*e
                     tmp3 = c*f
                     tmp4 = c*g1
                     tmp5 = d*g1
                     tmp6 = c*c
                     tmp7 = c*d
                     tmp8 = d*e
                     tmp9 = den1*den1
                     tmp10 = det_g*z
                     fg = (((3.0_ki*(2.0_ki*tmp4*num_b1*z+c*det_g*g1*z)+4.0_ki*(2.0_ki*(tmp1*num_b1-a&
                     &*f*num_b1)+a*c*e*g1-(tmp3*num_b1*z+tmp2*num_b1*z))+2.0_ki*(tmp1*e*e+d*d*num_b1+&
                     &tmp6*f*g1+tmp7*num_b1*z+tmp8*num_b1*z-(det_g*num_b1*z*z+tmp7*e*g1+tmp5*tmp6+tmp&
                     &4*tmp4))+tmp5*num_b1+det_g*e*g1*z)/tmp9+4.0_ki*(2.0_ki*(tmp3-tmp4)+tmp10-tmp8)*&
                     &(tmp10+2.0_ki*tmp2-tmp7)/(tmp9*den1)*g1*num_b1-1.0_ki/2.0_ki*(4.0_ki*(2.0_ki*a+&
                     &d)+2.0_ki*(3.0_ki*c*z+e*z)+g1)/den1)*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = d*d
                     tmp3 = c*g2
                     tmp4 = den2*den2
                     fg = (((8.0_ki*(a*g2*num_b2-a*f*num_b2)+4.0_ki*(tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2))&
                     &+tmp2*tmp3+2.0_ki*tmp2*num_b2+d*g2*num_b2)/tmp4+4.0_ki*(2.0_ki*(c*f-tmp3)-d*e)*&
                     &(2.0_ki*a*e-c*d)/(tmp4*den2)*g2*num_b2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*g3
                     tmp2 = a*c
                     tmp3 = g3*g3
                     tmp4 = det_g*num_b3
                     tmp5 = z*z
                     tmp6 = a*e
                     tmp7 = c*f
                     tmp8 = c*g3
                     tmp9 = d*d
                     tmp10 = c*d
                     tmp11 = d*e
                     tmp12 = d*g3
                     tmp13 = den3*den3
                     fg = ((4.0_ki*(2.0_ki*(tmp7-tmp8)+tmp5*det_g-tmp11)*(2.0_ki*tmp6+det_g*z-tmp10)/&
                     &(tmp13*den3)*g3*num_b3-(8.0_ki*(a*f*num_b3-tmp1*num_b3)+4.0_ki*(tmp2*f*g3+tmp4*&
                     &tmp5*z+tmp5*tmp6*num_b3+tmp7*num_b3*z-(a*d*e*g3+tmp2*tmp3))+2.0_ki*(tmp1*e*e+c*&
                     &c*f*g3+tmp5*a*det_g*g3-(3.0_ki*tmp8*num_b3*z+tmp11*num_b3*z+tmp10*e*g3+tmp10*tm&
                     &p5*num_b3+tmp9*num_b3+tmp4*tmp5+tmp8*tmp8))+tmp8*tmp9+tmp5*c*det_g*g3+det_g*e*g&
                     &3*z-(d*det_g*g3*z+tmp12*num_b3))/tmp13)*g3+(2.0_ki*(tmp12+tmp1*tmp5)+1.0_ki/2.0&
                     &_ki*tmp3+tmp5*tmp8+tmp5*num_b3*z+e*g3*z-(tmp5*num_b3+tmp12*z))/den3)
                     !
                  case(4)
                     !
                     tmp1 = a*e
                     tmp2 = c*f
                     tmp3 = c*g1
                     tmp4 = c*d
                     tmp5 = d*e
                     tmp6 = c*g3
                     tmp7 = z*z
                     fg = (1.0_ki/2.0_ki*(2.0_ki*(a*g3*z-num_b3*z)+tmp6*z+3.0_ki*tmp7*num_b3+e*g3-d*g&
                     &3)/den3*z-((2.0_ki*(tmp2+tmp1*z+tmp7*det_g-tmp6)-(tmp4*z+tmp5))/(den3*den3)*g3*&
                     &num_b3*z+(2.0_ki*(tmp1+tmp2+det_g*z-tmp3)-(tmp5+tmp4))/(den1*den1)*g1*num_b1*z+&
                     &1.0_ki/2.0_ki*(3.0_ki*tmp3+e*g1-num_b1*z)/den1*z))
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*f
                     tmp2 = a*c
                     tmp3 = z*z
                     tmp4 = a*e
                     tmp5 = a*g1
                     tmp6 = c*g1
                     tmp7 = d*d
                     tmp8 = c*d
                     tmp9 = den1*den1
                     fg = -(((2.0_ki*(2.0_ki*tmp5+d*g1)+tmp3*num_b1-tmp6*z)/den1*z-((8.0_ki*(a*a*e*g1&
                     &-tmp1*num_b1*z)+2.0_ki*(tmp7*num_b1*z+d*det_g*g1*z-tmp3*tmp8*num_b1)+4.0_ki*(tm&
                     &p2*f*g1+tmp3*tmp4*num_b1+3.0_ki*tmp5*num_b1*z+tmp3*det_g*num_b1*z+a*d*e*g1+a*de&
                     &t_g*g1*z-(tmp2*d*g1+tmp2*g1*g1))-(tmp3*c*det_g*g1+3.0_ki*tmp6*tmp7))/tmp9+4.0_k&
                     &i*(4.0_ki*(tmp1-tmp5)-(tmp3*det_g+tmp7))*(2.0_ki*tmp4+det_g*z-tmp8)/(tmp9*den1)&
                     &*g1*num_b1)*g1))
                     !
                  case(2)
                     !
                     tmp1 = 4.0_ki*(a*f-a*g2)-d*d
                     tmp2 = den2*den2
                     fg = (2.0_ki*(tmp1/tmp2*a+2.0_ki*(2.0_ki*a*e-c*d)/(tmp2*den2)*tmp1*num_b2)*g2*g2&
                     &)
                     !
                  case(3)
                     !
                     tmp1 = a*g3
                     tmp2 = a*f
                     tmp3 = a*c
                     tmp4 = z*z
                     tmp5 = a*e
                     tmp6 = c*g3
                     tmp7 = d*d
                     tmp8 = c*d
                     tmp9 = den3*den3
                     fg = -(((2.0_ki*(tmp1*z-d*g3)+tmp4*num_b3+tmp6*z)/den3*z-(4.0_ki*(4.0_ki*(tmp2-t&
                     &mp1)-(tmp4*det_g+tmp7))*(2.0_ki*tmp5+det_g*z-tmp8)/(tmp9*den3)*g3*num_b3-(8.0_k&
                     &i*(tmp2*num_b3*z+a*a*f*g3-tmp1*tmp1)+2.0_ki*(tmp4*tmp8*num_b3+d*det_g*g3*z-(tmp&
                     &4*a*det_g*g3+tmp7*num_b3*z+tmp1*tmp7))+4.0_ki*(tmp3*f*g3+a*d*e*g3-(tmp4*det_g*n&
                     &um_b3*z+tmp4*tmp5*num_b3+3.0_ki*tmp1*num_b3*z+tmp3*g3*g3))-(tmp4*c*det_g*g3+3.0&
                     &_ki*tmp6*tmp7))/tmp9)*g3))
                     !
                  case(4)
                     !
                     tmp1 = a*f
                     tmp2 = a*g1
                     tmp3 = z*z
                     tmp4 = tmp3*det_g+a*e*z
                     tmp5 = d*d
                     tmp6 = c*d*z
                     tmp7 = a*g3
                     fg = ((2.0_ki*(2.0_ki*(tmp2-tmp1)+tmp4)+tmp5-tmp6)/(den1*den1)*g1*num_b1*z+(2.0_&
                     &ki*(2.0_ki*(tmp7-tmp1)+tmp4)+tmp5-tmp6)/(den3*den3)*g3*num_b3*z-(1.0_ki/2.0_ki*&
                     &(2.0_ki*(tmp7*z-d*g3)+3.0_ki*tmp3*num_b3+c*g3*z)/den3*z+1.0_ki/2.0_ki*(2.0_ki*(&
                     &2.0_ki*tmp2+d*g1)+3.0_ki*tmp3*num_b1-c*g1*z)/den1*z))
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = e*g1
                     tmp2 = a*f
                     tmp3 = a*g1
                     tmp4 = c*g1
                     tmp5 = a*c
                     tmp6 = g1*g1
                     tmp7 = det_g*num_b1
                     tmp8 = z*z
                     tmp9 = a*e
                     tmp10 = tmp9*num_b1
                     tmp11 = c*f
                     tmp12 = d*d
                     tmp13 = c*det_g*g1
                     tmp14 = d*g1
                     tmp15 = c*c
                     tmp16 = num_b1*z
                     tmp17 = c*d
                     tmp18 = tmp17*num_b1
                     tmp19 = d*e
                     tmp20 = den1*den1
                     tmp21 = det_g*z
                     fg = ((-(4.0_ki*(tmp21+2.0_ki*tmp9-tmp17)*(2.0_ki*(tmp11-tmp4)+4.0_ki*(tmp2-tmp3&
                     &)+tmp21-(tmp8*det_g+tmp19+tmp12))/(tmp20*den1)*g1*num_b1+(8.0_ki*(tmp1*a*a-tmp2&
                     &*num_b1*z)+3.0_ki*(tmp13*z-tmp12*tmp4)+6.0_ki*(2.0_ki*tmp3*num_b1*z+tmp4*num_b1&
                     &*z)+2.0_ki*(tmp12*tmp16+tmp18*z+tmp3*e*e+tmp15*f*g1+tmp19*num_b1*z+d*det_g*g1*z&
                     &-(tmp17*e*g1+tmp18*tmp8+tmp14*tmp15+tmp4*tmp4))+4.0_ki*(tmp10*tmp8+tmp5*e*g1+tm&
                     &p5*f*g1+tmp7*tmp8*z+a*d*e*g1+a*det_g*g1*z-(tmp5*d*g1+tmp11*num_b1*z+tmp7*tmp8+t&
                     &mp5*tmp6+tmp10*z))+tmp14*num_b1+det_g*e*g1*z-tmp13*tmp8)/tmp20))*g1+(2.0_ki*(tm&
                     &p14*z+2.0_ki*tmp3*z)+1.0_ki/2.0_ki*tmp6+tmp1*z+tmp16*tmp8+3.0_ki*tmp4*z-(tmp8*n&
                     &um_b1+tmp4*tmp8))/den1)
                     !
                  case(2)
                     !
                     tmp1 = a*a
                     tmp2 = a*c
                     tmp3 = d*d
                     tmp4 = den2*den2
                     fg = ((-(4.0_ki*(2.0_ki*a*e-c*d)*(4.0_ki*(a*f-a*g2)+2.0_ki*(c*f-c*g2)-(d*e+tmp3)&
                     &)/(tmp4*den2)*num_b2+(8.0_ki*(tmp1*f-tmp1*g2)+4.0_ki*(tmp2*f-(a*d*e+tmp2*g2))+t&
                     &mp3*c+d*num_b2-2.0_ki*tmp3*a)/tmp4))*g2*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*g3
                     tmp2 = a*c
                     tmp3 = f*g3
                     tmp4 = a*f
                     tmp5 = c*f
                     tmp6 = c*g3
                     tmp7 = d*d
                     tmp8 = c*d
                     tmp9 = d*e
                     tmp10 = den3*den3
                     fg = (((4.0_ki*(2.0_ki*(tmp3*a*a+tmp2*f*g3+tmp4*num_b3*z-(tmp2*g3*g3+tmp1*tmp1))&
                     &+tmp5*num_b3*z)+2.0_ki*(3.0_ki*(-(tmp6*num_b3*z+2.0_ki*tmp1*num_b3*z))+tmp1*e*e&
                     &+tmp3*c*c-(tmp9*num_b3*z+tmp8*e*g3+tmp7*num_b3*z+tmp6*tmp7+tmp1*tmp7+tmp6*tmp6)&
                     &)+d*det_g*g3*z+det_g*e*g3*z-d*g3*num_b3)/tmp10-(4.0_ki*(2.0_ki*(2.0_ki*(tmp4-tm&
                     &p1)+tmp5-tmp6)-(tmp9+tmp7))*(2.0_ki*a*e+det_g*z-tmp8)/(tmp10*den3)*g3*num_b3+1.&
                     &0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = a*f
                     tmp2 = c*f
                     tmp3 = d*d
                     tmp4 = d*e
                     tmp5 = c*g1
                     tmp6 = z*z
                     tmp7 = a*g1
                     tmp8 = a*e
                     tmp9 = c*d
                     fg = (((4.0_ki*(tmp1-a*g3)+2.0_ki*(tmp2-c*g3)-(tmp4+tmp3))/(den3*den3)*num_b3*z-&
                     &1.0_ki/2.0_ki*(d+e)/den3*z)*g3+1.0_ki/2.0_ki*(2.0_ki*(2.0_ki*tmp7+d*g1)+3.0_ki*&
                     &(tmp5+tmp6*num_b1-num_b1*z)+e*g1-tmp5*z)/den1*z-(4.0_ki*(tmp7-tmp1)+2.0_ki*(tmp&
                     &5+tmp6*det_g+tmp8*z-(det_g*z+tmp8+tmp2))+tmp3+tmp4+tmp9-tmp9*z)/(den1*den1)*g1*&
                     &num_b1*z)
                     !
                  end select
                  !
               end select
               !
            case(3)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = b*c
                     tmp3 = a*e
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = c*g1
                     tmp7 = e*e
                     tmp8 = c*det_g*g1
                     tmp9 = tmp6*num_b1
                     tmp10 = d*g1
                     tmp11 = c*c
                     tmp12 = e*g1
                     tmp13 = det_g*num_b1
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = d*e
                     tmp17 = det_g*e*g1
                     tmp18 = den1*den1
                     tmp19 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(3.0_ki*c*z+e*z)+4.0_ki*(2.0_ki*a+d+b*z)+g1)/den1-(&
                     &4.0_ki*(2.0_ki*(tmp5-tmp6)+tmp19-tmp16)*(2.0_ki*(tmp3-tmp4)+tmp15+tmp19-(det_g+&
                     &tmp14))/(tmp18*den1)*g1*num_b1+(3.0_ki*(tmp8*z+2.0_ki*tmp9*z)+4.0_ki*(2.0_ki*(t&
                     &mp1*num_b1-a*f*num_b1)+tmp2*f*g1+tmp4*num_b1*z+a*c*e*g1-(tmp5*num_b1*z+tmp4*e*g&
                     &1+tmp3*num_b1*z+tmp2*d*g1+tmp2*g1*g1))+2.0_ki*(tmp1*tmp7+tmp11*tmp12+tmp13*z+d*&
                     &d*num_b1+tmp11*f*g1+tmp14*num_b1*z+tmp16*num_b1*z+b*det_g*g1*z-(tmp15*num_b1*z+&
                     &tmp14*e*g1+tmp13*z*z+tmp10*tmp11+tmp9+tmp8+tmp6*tmp6))+tmp10*num_b1+tmp17*z+tmp&
                     &6*tmp7-(tmp12*num_b1+tmp17))/tmp18))*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = c*f
                     tmp3 = tmp2*num_b2
                     tmp4 = c*g2
                     tmp5 = tmp4*num_b2
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = tmp7*num_b2
                     tmp9 = c*d
                     tmp10 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z-d)+g2)/den2-(4.0_ki*(2.0_ki*(tmp2-tmp4)-tmp7)*&
                     &(2.0_ki*(a*e-b*d)+c*e-(det_g*z+tmp9))/(tmp10*den2)*g2*num_b2+(4.0_ki*(2.0_ki*(a&
                     &*g2*num_b2-a*f*num_b2)+tmp5+tmp3*z+tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2+tmp3))+2.0_ki&
                     &*(tmp8+tmp6*num_b2+tmp6*b*g2+c*c*f*g2-(tmp9*e*g2+tmp8*z+3.0_ki*tmp5*z+tmp4*tmp4&
                     &))+tmp4*tmp6+d*g2*num_b2+d*det_g*g2*z-e*g2*num_b2)/tmp10))*g2)
                     !
                  case(3)
                     !
                     tmp1 = z*z
                     tmp2 = a*g3
                     tmp3 = c*g3
                     tmp4 = a*c
                     tmp5 = g3*g3
                     tmp6 = b*c
                     tmp7 = c*f
                     tmp8 = a*e
                     tmp9 = b*d
                     tmp10 = c*d
                     tmp11 = d*d
                     tmp12 = e*e
                     tmp13 = d*e
                     tmp14 = c*e
                     tmp15 = den3*den3
                     fg = (((8.0_ki*(a*f*num_b3-tmp2*num_b3)+2.0_ki*(tmp12*tmp2+tmp1*tmp14*num_b3+tmp&
                     &11*b*g3+tmp1*a*det_g*g3+tmp1*b*det_g*g3+tmp1*c*det_g*g3-(tmp1*tmp10*num_b3+tmp1&
                     &3*num_b3+tmp11*num_b3))+4.0_ki*(tmp7*num_b3+tmp1*tmp8*num_b3+tmp4*f*g3+tmp6*f*g&
                     &3+c*c*f*g3-(a*d*e*g3+tmp9*e*g3+tmp10*e*g3+tmp1*tmp9*num_b3+tmp5*tmp6+tmp4*tmp5+&
                     &tmp3*num_b3+tmp3*tmp3))+tmp11*tmp3+tmp12*tmp3+e*g3*num_b3-d*g3*num_b3)/tmp15-(4&
                     &.0_ki*(2.0_ki*(tmp7-tmp3)+tmp1*det_g-tmp13)*(2.0_ki*(tmp8-tmp9)+tmp14-tmp10)/(t&
                     &mp15*den3)*g3*num_b3+(2.0_ki*(tmp1*a+tmp1*b+tmp1*c)+d+e+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = z*z
                     tmp2 = a*e
                     tmp3 = b*d
                     tmp4 = c*d
                     tmp5 = c*e
                     tmp6 = c*f
                     tmp7 = d*e
                     tmp8 = c*g1
                     tmp9 = e*g1
                     tmp10 = tmp8*z
                     fg = (((2.0_ki*(tmp2-tmp3)+tmp5-tmp4)/(den3*den3)*tmp1*num_b3-(a+b+c)*tmp1/den3)&
                     &*g3+(1.0_ki/2.0_ki*d/den2*z-(2.0_ki*(tmp6-c*g2)-tmp7)/(den2*den2)*num_b2*z)*g2+&
                     &(3.0_ki/2.0_ki*tmp10+1.0_ki/2.0_ki*tmp9*z+1.0_ki/2.0_ki*num_b1*z+b*g1*z-(1.0_ki&
                     &/2.0_ki*tmp1*num_b1+1.0_ki/2.0_ki*tmp9+tmp8))/den1+(2.0_ki*(tmp8+tmp1*det_g+tmp&
                     &2*z+tmp6*z-(det_g*z+tmp3*z+tmp6+tmp10))+tmp7+tmp5*z-(tmp7*z+tmp4*z))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = e*g1
                     tmp2 = a*f
                     tmp3 = a*g1
                     tmp4 = c*g1
                     tmp5 = a*c
                     tmp6 = g1*g1
                     tmp7 = det_g*num_b1
                     tmp8 = z*z
                     tmp9 = a*e
                     tmp10 = tmp9*num_b1
                     tmp11 = c*f
                     tmp12 = d*d
                     tmp13 = c*det_g*g1
                     tmp14 = d*g1
                     tmp15 = c*c
                     tmp16 = num_b1*z
                     tmp17 = c*d
                     tmp18 = tmp17*num_b1
                     tmp19 = d*e
                     tmp20 = den1*den1
                     tmp21 = det_g*z
                     fg = ((-(4.0_ki*(tmp21+2.0_ki*tmp9-tmp17)*(2.0_ki*(tmp11-tmp4)+4.0_ki*(tmp2-tmp3&
                     &)+tmp21-(tmp8*det_g+tmp19+tmp12))/(tmp20*den1)*g1*num_b1+(8.0_ki*(tmp1*a*a-tmp2&
                     &*num_b1*z)+3.0_ki*(tmp13*z-tmp12*tmp4)+6.0_ki*(2.0_ki*tmp3*num_b1*z+tmp4*num_b1&
                     &*z)+2.0_ki*(tmp12*tmp16+tmp18*z+tmp3*e*e+tmp15*f*g1+tmp19*num_b1*z+d*det_g*g1*z&
                     &-(tmp17*e*g1+tmp18*tmp8+tmp14*tmp15+tmp4*tmp4))+4.0_ki*(tmp10*tmp8+tmp5*e*g1+tm&
                     &p5*f*g1+tmp7*tmp8*z+a*d*e*g1+a*det_g*g1*z-(tmp5*d*g1+tmp11*num_b1*z+tmp7*tmp8+t&
                     &mp5*tmp6+tmp10*z))+tmp14*num_b1+det_g*e*g1*z-tmp13*tmp8)/tmp20))*g1+(2.0_ki*(tm&
                     &p14*z+2.0_ki*tmp3*z)+1.0_ki/2.0_ki*tmp6+tmp1*z+tmp16*tmp8+3.0_ki*tmp4*z-(tmp8*n&
                     &um_b1+tmp4*tmp8))/den1)
                     !
                  case(2)
                     !
                     tmp1 = a*a
                     tmp2 = a*c
                     tmp3 = d*d
                     tmp4 = den2*den2
                     fg = ((-(4.0_ki*(2.0_ki*a*e-c*d)*(4.0_ki*(a*f-a*g2)+2.0_ki*(c*f-c*g2)-(d*e+tmp3)&
                     &)/(tmp4*den2)*num_b2+(8.0_ki*(tmp1*f-tmp1*g2)+4.0_ki*(tmp2*f-(a*d*e+tmp2*g2))+t&
                     &mp3*c+d*num_b2-2.0_ki*tmp3*a)/tmp4))*g2*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*g3
                     tmp2 = a*c
                     tmp3 = f*g3
                     tmp4 = a*f
                     tmp5 = c*f
                     tmp6 = c*g3
                     tmp7 = d*d
                     tmp8 = c*d
                     tmp9 = d*e
                     tmp10 = den3*den3
                     fg = (((4.0_ki*(2.0_ki*(tmp3*a*a+tmp2*f*g3+tmp4*num_b3*z-(tmp2*g3*g3+tmp1*tmp1))&
                     &+tmp5*num_b3*z)+2.0_ki*(3.0_ki*(-(tmp6*num_b3*z+2.0_ki*tmp1*num_b3*z))+tmp1*e*e&
                     &+tmp3*c*c-(tmp9*num_b3*z+tmp8*e*g3+tmp7*num_b3*z+tmp6*tmp7+tmp1*tmp7+tmp6*tmp6)&
                     &)+d*det_g*g3*z+det_g*e*g3*z-d*g3*num_b3)/tmp10-(4.0_ki*(2.0_ki*(2.0_ki*(tmp4-tm&
                     &p1)+tmp5-tmp6)-(tmp9+tmp7))*(2.0_ki*a*e+det_g*z-tmp8)/(tmp10*den3)*g3*num_b3+1.&
                     &0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = a*f
                     tmp2 = c*f
                     tmp3 = d*d
                     tmp4 = d*e
                     tmp5 = c*g1
                     tmp6 = z*z
                     tmp7 = a*g1
                     tmp8 = a*e
                     tmp9 = c*d
                     fg = (((4.0_ki*(tmp1-a*g3)+2.0_ki*(tmp2-c*g3)-(tmp4+tmp3))/(den3*den3)*num_b3*z-&
                     &1.0_ki/2.0_ki*(d+e)/den3*z)*g3+1.0_ki/2.0_ki*(2.0_ki*(2.0_ki*tmp7+d*g1)+3.0_ki*&
                     &(tmp5+tmp6*num_b1-num_b1*z)+e*g1-tmp5*z)/den1*z-(4.0_ki*(tmp7-tmp1)+2.0_ki*(tmp&
                     &5+tmp6*det_g+tmp8*z-(det_g*z+tmp8+tmp2))+tmp3+tmp4+tmp9-tmp9*z)/(den1*den1)*g1*&
                     &num_b1*z)
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = c*g1
                     tmp3 = tmp2*num_b1
                     tmp4 = e*g1
                     tmp5 = a*a
                     tmp6 = tmp4*tmp5
                     tmp7 = det_g*num_b1
                     tmp8 = z*z
                     tmp9 = a*c
                     tmp10 = tmp9*e
                     tmp11 = tmp10*g1
                     tmp12 = a*e
                     tmp13 = tmp12*num_b1
                     tmp14 = c*f*num_b1
                     tmp15 = c*det_g
                     tmp16 = tmp15*g1
                     tmp17 = tmp16*z
                     tmp18 = g1*g1
                     tmp19 = e*e
                     tmp20 = b*c
                     tmp21 = d*g1
                     tmp22 = c*c
                     tmp23 = tmp21*tmp22
                     tmp24 = tmp7*z
                     tmp25 = tmp9*d
                     tmp26 = tmp25*g1
                     tmp27 = a*d
                     tmp28 = a*det_g
                     tmp29 = tmp28*g1*z
                     tmp30 = tmp20*d
                     tmp31 = tmp30*g1
                     tmp32 = b*d
                     tmp33 = c*d
                     tmp34 = tmp33*num_b1
                     tmp35 = d*e*num_b1
                     tmp36 = d*d
                     tmp37 = tmp22*tmp4
                     tmp38 = num_b1*z
                     tmp39 = c*e
                     tmp40 = d*det_g
                     tmp41 = det_g*e*g1
                     tmp42 = den1*den1
                     tmp43 = e*f
                     tmp44 = tmp28*e
                     tmp45 = tmp15*f
                     tmp46 = det_g*z
                     tmp47 = tmp46*tmp46
                     tmp48 = tmp33*det_g
                     tmp49 = tmp40*e
                     fg = (((6.0_ki*(tmp17-tmp3)+12.0_ki*(tmp3*z+tmp1*num_b1*z)+8.0_ki*(tmp11+tmp6-(a&
                     &*f*num_b1*z+tmp7*tmp8+tmp14*z+tmp13*z))+2.0_ki*(tmp37+tmp21*num_b1+tmp36*tmp38+&
                     &tmp41*z+tmp40*g1*z+b*det_g*g1*z-(tmp39*num_b1*z+tmp34*tmp8+tmp35+tmp16))+4.0_ki&
                     &*(tmp14+tmp24+tmp29+tmp1*tmp19+tmp13*tmp8+tmp24*tmp8+tmp34*z+tmp35*z+tmp20*f*g1&
                     &+tmp22*f*g1+tmp27*e*g1+tmp32*num_b1*z+tmp9*f*g1-(tmp33*e*g1+tmp32*e*g1+tmp18*tm&
                     &p9+tmp18*tmp20+tmp31+tmp26+tmp23+tmp2*tmp2))+tmp19*tmp2-(tmp4*num_b1+3.0_ki*tmp&
                     &2*tmp36+tmp16*tmp8+tmp41))/tmp42+4.0_ki*(8.0_ki*(tmp10*f+tmp43*tmp5-(tmp6+tmp11&
                     &))+2.0_ki*(tmp16+tmp47+tmp22*tmp43+tmp36*tmp39+tmp36*b*e-(tmp32*det_g*z+tmp49*z&
                     &+tmp48*z+tmp44*tmp8+tmp12*tmp36+tmp45+tmp37))+4.0_ki*(tmp23+tmp26+tmp31+tmp44*z&
                     &+tmp45*z+tmp28*f*z-(tmp22*d*f+tmp30*f+tmp25*f+tmp19*tmp27+tmp29+tmp17))+tmp49+t&
                     &mp33*tmp36+tmp48*tmp8+tmp15*e*z-(det_g*det_g*z+tmp47*z+tmp36*tmp46+tmp19*tmp33)&
                     &)/(tmp42*den1)*g1*num_b1)*g1+(2.0_ki*(tmp2+tmp8*num_b1-(b*g1*z+tmp4*z+tmp21*z+3&
                     &.0_ki*tmp2*z+2.0_ki*tmp1*z))+tmp4+tmp2*tmp8-(tmp38*tmp8+tmp38+tmp18))/den1)
                     !
                  case(2)
                     !
                     tmp1 = a*g2
                     tmp2 = a*c
                     tmp3 = f*g2
                     tmp4 = a*a
                     tmp5 = a*d
                     tmp6 = c*g2
                     tmp7 = d*d
                     tmp8 = c*c
                     tmp9 = d*g2
                     tmp10 = c*d
                     tmp11 = e*g2
                     tmp12 = d*det_g
                     tmp13 = den2*den2
                     tmp14 = e*f
                     tmp15 = tmp2*e
                     tmp16 = e*e
                     tmp17 = tmp2*d
                     tmp18 = b*c*d
                     tmp19 = c*det_g
                     fg = (((4.0_ki*(2.0_ki*(tmp3*tmp4+tmp2*f*g2-(tmp5*e*g2+tmp2*g2*g2+tmp1*tmp1))+c*&
                     &f*num_b2*z)+2.0_ki*(tmp3*tmp8+tmp6*tmp7+tmp9*num_b2+tmp7*b*g2-(d*e*num_b2*z+3.0&
                     &_ki*tmp6*num_b2*z+tmp10*e*g2+tmp1*tmp7+tmp6*tmp6))+tmp12*g2*z-tmp11*num_b2)/tmp&
                     &13+4.0_ki*(2.0_ki*(tmp14*tmp8+tmp19*g2*z+tmp7*b*e+tmp7*c*e-(tmp7*a*e+tmp19*f*z+&
                     &tmp11*tmp8))+4.0_ki*(2.0_ki*(tmp14*tmp4+tmp15*f-(tmp15*g2+tmp11*tmp4))+tmp17*g2&
                     &+tmp18*g2+tmp8*tmp9-(tmp8*d*f+tmp18*f+tmp17*f+tmp16*tmp5))+tmp10*tmp7+tmp12*e*z&
                     &-tmp10*tmp16)/(tmp13*den2)*g2*num_b2-1.0_ki/2.0_ki*(g2+2.0_ki*d*z)/den2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*c
                     tmp2 = g3*g3
                     tmp3 = a*g3
                     tmp4 = c*g3
                     tmp5 = f*g3
                     tmp6 = c*c
                     tmp7 = c*d
                     tmp8 = a*a
                     tmp9 = e*e
                     tmp10 = b*c
                     tmp11 = a*d
                     tmp12 = a*e
                     tmp13 = z*z
                     tmp14 = b*d
                     tmp15 = d*d
                     tmp16 = d*g3
                     tmp17 = a*det_g
                     tmp18 = tmp17*g3
                     tmp19 = c*det_g
                     tmp20 = tmp19*g3
                     tmp21 = c*e
                     tmp22 = e*g3
                     tmp23 = d*det_g
                     tmp24 = den3*den3
                     tmp25 = e*f
                     tmp26 = tmp1*e
                     tmp27 = tmp1*d
                     tmp28 = tmp10*d
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+4.0_ki*(tmp13*a+tmp13*b+tmp13*c)+3.0_ki*g&
                     &3)/den3+4.0_ki*(2.0_ki*(tmp15*tmp21+tmp25*tmp6+tmp13*tmp17*e+tmp15*b*e+tmp19*f*&
                     &z-(tmp13*tmp14*det_g+tmp22*tmp6+tmp20*z+tmp12*tmp15))+4.0_ki*(2.0_ki*(tmp25*tmp&
                     &8+tmp26*f-(tmp26*g3+tmp22*tmp8))+tmp16*tmp6+tmp27*g3+tmp28*g3+tmp17*f*z-(tmp6*d&
                     &*f+tmp28*f+tmp27*f+tmp18*z+tmp11*tmp9))+tmp15*tmp7+tmp13*tmp19*e-(tmp23*e*z+tmp&
                     &15*det_g*z+tmp13*tmp7*det_g+tmp7*tmp9))/(tmp24*den3)*g3*num_b3-(6.0_ki*(2.0_ki*&
                     &(tmp1*f*g3-(tmp3*num_b3*z+tmp1*tmp2))+tmp5*tmp6-(tmp7*e*g3+tmp4*num_b3*z+tmp4*t&
                     &mp4))+4.0_ki*(2.0_ki*(tmp5*tmp8+a*f*num_b3*z-tmp3*tmp3)+tmp3*tmp9+tmp10*f*g3+tm&
                     &p12*tmp13*num_b3+c*f*num_b3*z-(tmp14*e*g3+tmp13*tmp14*num_b3+tmp11*e*g3+tmp10*t&
                     &mp2))+2.0_ki*(tmp13*tmp18+tmp13*tmp20+tmp13*tmp21*num_b3+tmp15*b*g3+tmp13*b*det&
                     &_g*g3-(d*e*num_b3*z+tmp15*num_b3*z+tmp13*tmp7*num_b3+tmp16*num_b3+tmp15*tmp3))+&
                     &tmp22*num_b3+tmp4*tmp9+tmp23*g3*z+det_g*e*g3*z-tmp15*tmp4)/tmp24)*g3)
                     !
                  case(4)
                     !
                     tmp1 = a*f
                     tmp2 = c*f
                     tmp3 = a*e
                     tmp4 = tmp3*z
                     tmp5 = b*d*z
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = c*d
                     tmp9 = tmp8*z
                     tmp10 = c*e*z
                     tmp11 = z*z
                     tmp12 = c*g1
                     tmp13 = tmp12*z
                     tmp14 = num_b1*z
                     tmp15 = a*g1*z
                     tmp16 = e*g1
                     tmp17 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(a*z+b*z+c*z)+d+e)/den3*z-(4.0_ki*(tmp1-a*g3)+2.0_k&
                     &i*(tmp2+tmp4-(c*g3+tmp5))+tmp10-(tmp9+tmp7+tmp6))/(den3*den3)*num_b3*z)*g3+((2.&
                     &0_ki*(tmp2-c*g2)-tmp7)/(den2*den2)*num_b2*z-1.0_ki/2.0_ki*d/den2*z)*g2+(3.0_ki*&
                     &(tmp11*num_b1-tmp13)+3.0_ki/2.0_ki*(-(tmp11*tmp14+tmp14))+tmp12+1.0_ki/2.0_ki*t&
                     &mp16+1.0_ki/2.0_ki*tmp11*tmp12-(d*g1*z+b*g1*z+tmp16*z+2.0_ki*tmp15))/den1+(4.0_&
                     &ki*(tmp13+tmp15-(tmp2*z+tmp11*det_g+tmp1*z+tmp4))+2.0_ki*(tmp17+tmp2+tmp5+tmp9+&
                     &tmp11*tmp17+tmp11*tmp3+tmp7*z-tmp12)+tmp6*z-(tmp11*tmp8+tmp7+tmp10))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               end select
               !
            end select
            !
         case(3)
            !
            select case(par2)
            !
            case(1)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*g1
                     tmp2 = b*d
                     tmp3 = tmp2*num_b1
                     tmp4 = b*f
                     tmp5 = b*b
                     tmp6 = b*c
                     tmp7 = c*e
                     tmp8 = tmp7*num_b1
                     tmp9 = c*f
                     tmp10 = det_g*num_b1
                     tmp11 = e*e
                     tmp12 = b*det_g*g1
                     tmp13 = c*g1
                     tmp14 = e*g1
                     tmp15 = d*e
                     tmp16 = den1*den1
                     fg = (((c+e+2.0_ki*b*z)/den1+4.0_ki*(2.0_ki*tmp2+det_g-tmp7)*(4.0_ki*(tmp4-tmp1)&
                     &+2.0_ki*(tmp9-tmp13)+det_g*z-(det_g+tmp15+tmp11))/(tmp16*den1)*g1*num_b1-(8.0_k&
                     &i*(tmp5*f*g1-(tmp5*d*g1+tmp4*num_b1+tmp3+tmp1*tmp1))+4.0_ki*(tmp8+tmp1*num_b1+t&
                     &mp3*z+tmp6*e*g1+tmp6*f*g1+a*e*num_b1-(tmp6*d*g1+tmp2*e*g1+tmp9*num_b1+tmp6*g1*g&
                     &1))+2.0_ki*(tmp10*z+tmp11*num_b1+tmp12*z+tmp13*num_b1+tmp14*c*c+tmp15*num_b1-(c&
                     &*det_g*g1+c*d*num_b1+tmp8*z+tmp1*tmp11+3.0_ki*tmp12+tmp10))+tmp11*tmp13-(det_g*&
                     &e*g1+tmp14*num_b1))/tmp16)*g1)
                     !
                  case(2)
                     !
                     tmp1 = b*f
                     tmp2 = tmp1*num_b2
                     tmp3 = b*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = c*f
                     tmp6 = tmp5*num_b2
                     tmp7 = c*g2
                     tmp8 = tmp7*num_b2
                     tmp9 = e*e
                     tmp10 = b*c
                     tmp11 = g2*g2
                     tmp12 = b*d
                     tmp13 = tmp12*num_b2
                     tmp14 = d*e
                     tmp15 = tmp14*num_b2
                     tmp16 = det_g*num_b2
                     tmp17 = z*z
                     tmp18 = tmp16*z
                     tmp19 = c*d
                     tmp20 = c*e
                     tmp21 = tmp20*num_b2
                     tmp22 = num_b2*z
                     tmp23 = e*g2
                     tmp24 = den2*den2
                     tmp25 = d*g2
                     fg = ((4.0_ki*(2.0_ki*tmp12+det_g*z-tmp20)*(4.0_ki*(tmp1-tmp3)+2.0_ki*(tmp5-tmp7&
                     &)-(tmp17*det_g+tmp9+tmp14))/(tmp24*den2)*g2*num_b2-(16.0_ki*(tmp4-tmp2)+6.0_ki*&
                     &(-(tmp8*z+2.0_ki*tmp4*z))+8.0_ki*(tmp8+tmp2*z-tmp6)+4.0_ki*(tmp15+tmp16*tmp17+t&
                     &mp6*z+tmp9*num_b2+tmp10*f*g2+tmp12*e*g2+a*e*num_b2-(tmp17*tmp18+tmp13*tmp17+tmp&
                     &10*tmp11+tmp13))+2.0_ki*(tmp21+tmp17*tmp21+tmp3*d*d+c*c*f*g2+det_g*e*g2*z-(tmp1&
                     &9*e*g2+tmp22*tmp9+tmp19*num_b2+tmp15*z+tmp7*tmp7+tmp18))+d*det_g*g2*z-(tmp17*c*&
                     &det_g*g2+3.0_ki*tmp7*tmp9+tmp23*num_b2))/tmp24)*g2+(2.0_ki*(tmp17*num_b2+tmp23*&
                     &z-(a*g2+tmp25+2.0_ki*tmp23))+1.0_ki/2.0_ki*tmp11+tmp25*z-(tmp17*tmp7+tmp17*tmp2&
                     &2+tmp7+tmp22))/den2)
                     !
                  case(3)
                     !
                     tmp1 = b*f
                     tmp2 = tmp1*num_b3
                     tmp3 = b*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = b*c
                     tmp6 = f*g3
                     tmp7 = c*f
                     tmp8 = tmp7*num_b3
                     tmp9 = c*g3
                     tmp10 = tmp9*num_b3
                     tmp11 = e*e
                     tmp12 = b*d
                     tmp13 = d*e
                     tmp14 = tmp13*num_b3
                     tmp15 = c*d
                     tmp16 = c*e
                     tmp17 = den3*den3
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(-(e*z+d*z))+4.0_ki*(a+b+c+d+e)-g3)/den3+4.0_ki*(2.&
                     &0_ki*(2.0_ki*(tmp1-tmp3)+tmp7-tmp9)-(tmp13+tmp11))*(2.0_ki*tmp12+det_g*z-tmp16)&
                     &/(tmp17*den3)*g3*num_b3-(4.0_ki*(tmp14+tmp11*num_b3+tmp8*z+a*e*num_b3-tmp12*num&
                     &_b3)+8.0_ki*(2.0_ki*(tmp4-tmp2)+tmp10+tmp3*tmp3+tmp2*z+tmp5*g3*g3-(tmp5*f*g3+tm&
                     &p6*b*b+tmp8))+2.0_ki*(3.0_ki*(-(2.0_ki*tmp4*z+tmp10*z))+tmp9*tmp9+tmp11*tmp3+tm&
                     &p11*tmp9+tmp16*num_b3+tmp15*e*g3-(tmp11*num_b3*z+tmp6*c*c+tmp3*d*d+tmp15*num_b3&
                     &+tmp14*z))-(det_g*e*g3*z+d*det_g*g3*z+e*g3*num_b3))/tmp17)*g3)
                     !
                  case(4)
                     !
                     tmp1 = b*g1
                     tmp2 = c*g1
                     tmp3 = b*f
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = tmp4*z
                     tmp7 = e*e
                     tmp8 = c*e
                     tmp9 = d*e
                     tmp10 = tmp8*z
                     tmp11 = c*g2
                     tmp12 = z*z
                     fg = ((-((4.0_ki*(tmp3-b*g3)+2.0_ki*(tmp5-c*g3)-(tmp9+tmp7))/(den3*den3)*num_b3*&
                     &z+1.0_ki/2.0_ki*(d+e)/den3*z))*g3+(tmp1*z+1.0_ki/2.0_ki*num_b1*z-(1.0_ki/2.0_ki&
                     &*e*g1+1.0_ki/2.0_ki*num_b1+tmp2+3.0_ki*tmp1))/den1+(2.0_ki*(2.0_ki*(b*g2-tmp3)+&
                     &tmp11+tmp6+tmp12*det_g-tmp5)+tmp7+tmp9-tmp10)/(den2*den2)*g2*num_b2*z-((2.0_ki*&
                     &(2.0_ki*(tmp3-tmp1)+tmp5+tmp6+det_g*z-(det_g+tmp4+tmp2))+tmp8-(tmp9+tmp7+tmp10)&
                     &)/(den1*den1)*g1*num_b1+1.0_ki/2.0_ki*(2.0_ki*(num_b2-(2.0_ki*num_b2*z+e*g2))+t&
                     &mp11*z+3.0_ki*tmp12*num_b2-d*g2)/den2*z))
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = b*c
                     tmp3 = a*e
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = c*g1
                     tmp7 = e*e
                     tmp8 = c*det_g*g1
                     tmp9 = tmp6*num_b1
                     tmp10 = d*g1
                     tmp11 = c*c
                     tmp12 = e*g1
                     tmp13 = det_g*num_b1
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = d*e
                     tmp17 = det_g*e*g1
                     tmp18 = den1*den1
                     tmp19 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(3.0_ki*c*z+e*z)+4.0_ki*(2.0_ki*a+d+b*z)+g1)/den1-(&
                     &4.0_ki*(2.0_ki*(tmp5-tmp6)+tmp19-tmp16)*(2.0_ki*(tmp3-tmp4)+tmp15+tmp19-(det_g+&
                     &tmp14))/(tmp18*den1)*g1*num_b1+(3.0_ki*(tmp8*z+2.0_ki*tmp9*z)+4.0_ki*(2.0_ki*(t&
                     &mp1*num_b1-a*f*num_b1)+tmp2*f*g1+tmp4*num_b1*z+a*c*e*g1-(tmp5*num_b1*z+tmp4*e*g&
                     &1+tmp3*num_b1*z+tmp2*d*g1+tmp2*g1*g1))+2.0_ki*(tmp1*tmp7+tmp11*tmp12+tmp13*z+d*&
                     &d*num_b1+tmp11*f*g1+tmp14*num_b1*z+tmp16*num_b1*z+b*det_g*g1*z-(tmp15*num_b1*z+&
                     &tmp14*e*g1+tmp13*z*z+tmp10*tmp11+tmp9+tmp8+tmp6*tmp6))+tmp10*num_b1+tmp17*z+tmp&
                     &6*tmp7-(tmp12*num_b1+tmp17))/tmp18))*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = c*f
                     tmp3 = tmp2*num_b2
                     tmp4 = c*g2
                     tmp5 = tmp4*num_b2
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = tmp7*num_b2
                     tmp9 = c*d
                     tmp10 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z-d)+g2)/den2-(4.0_ki*(2.0_ki*(tmp2-tmp4)-tmp7)*&
                     &(2.0_ki*(a*e-b*d)+c*e-(det_g*z+tmp9))/(tmp10*den2)*g2*num_b2+(4.0_ki*(2.0_ki*(a&
                     &*g2*num_b2-a*f*num_b2)+tmp5+tmp3*z+tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2+tmp3))+2.0_ki&
                     &*(tmp8+tmp6*num_b2+tmp6*b*g2+c*c*f*g2-(tmp9*e*g2+tmp8*z+3.0_ki*tmp5*z+tmp4*tmp4&
                     &))+tmp4*tmp6+d*g2*num_b2+d*det_g*g2*z-e*g2*num_b2)/tmp10))*g2)
                     !
                  case(3)
                     !
                     tmp1 = z*z
                     tmp2 = a*g3
                     tmp3 = c*g3
                     tmp4 = a*c
                     tmp5 = g3*g3
                     tmp6 = b*c
                     tmp7 = c*f
                     tmp8 = a*e
                     tmp9 = b*d
                     tmp10 = c*d
                     tmp11 = d*d
                     tmp12 = e*e
                     tmp13 = d*e
                     tmp14 = c*e
                     tmp15 = den3*den3
                     fg = (((8.0_ki*(a*f*num_b3-tmp2*num_b3)+2.0_ki*(tmp12*tmp2+tmp1*tmp14*num_b3+tmp&
                     &11*b*g3+tmp1*a*det_g*g3+tmp1*b*det_g*g3+tmp1*c*det_g*g3-(tmp1*tmp10*num_b3+tmp1&
                     &3*num_b3+tmp11*num_b3))+4.0_ki*(tmp7*num_b3+tmp1*tmp8*num_b3+tmp4*f*g3+tmp6*f*g&
                     &3+c*c*f*g3-(a*d*e*g3+tmp9*e*g3+tmp10*e*g3+tmp1*tmp9*num_b3+tmp5*tmp6+tmp4*tmp5+&
                     &tmp3*num_b3+tmp3*tmp3))+tmp11*tmp3+tmp12*tmp3+e*g3*num_b3-d*g3*num_b3)/tmp15-(4&
                     &.0_ki*(2.0_ki*(tmp7-tmp3)+tmp1*det_g-tmp13)*(2.0_ki*(tmp8-tmp9)+tmp14-tmp10)/(t&
                     &mp15*den3)*g3*num_b3+(2.0_ki*(tmp1*a+tmp1*b+tmp1*c)+d+e+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = z*z
                     tmp2 = a*e
                     tmp3 = b*d
                     tmp4 = c*d
                     tmp5 = c*e
                     tmp6 = c*f
                     tmp7 = d*e
                     tmp8 = c*g1
                     tmp9 = e*g1
                     tmp10 = tmp8*z
                     fg = (((2.0_ki*(tmp2-tmp3)+tmp5-tmp4)/(den3*den3)*tmp1*num_b3-(a+b+c)*tmp1/den3)&
                     &*g3+(1.0_ki/2.0_ki*d/den2*z-(2.0_ki*(tmp6-c*g2)-tmp7)/(den2*den2)*num_b2*z)*g2+&
                     &(3.0_ki/2.0_ki*tmp10+1.0_ki/2.0_ki*tmp9*z+1.0_ki/2.0_ki*num_b1*z+b*g1*z-(1.0_ki&
                     &/2.0_ki*tmp1*num_b1+1.0_ki/2.0_ki*tmp9+tmp8))/den1+(2.0_ki*(tmp8+tmp1*det_g+tmp&
                     &2*z+tmp6*z-(det_g*z+tmp3*z+tmp6+tmp10))+tmp7+tmp5*z-(tmp7*z+tmp4*z))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*g1
                     tmp2 = a*g1
                     tmp3 = b*c
                     tmp4 = d*g1
                     tmp5 = b*b
                     tmp6 = tmp4*tmp5
                     tmp7 = f*g1
                     tmp8 = tmp3*d
                     tmp9 = tmp8*g1
                     tmp10 = b*d
                     tmp11 = tmp10*num_b1
                     tmp12 = c*det_g
                     tmp13 = tmp12*g1
                     tmp14 = c*g1
                     tmp15 = tmp14*num_b1
                     tmp16 = e*g1
                     tmp17 = c*c
                     tmp18 = tmp16*tmp17
                     tmp19 = det_g*num_b1
                     tmp20 = a*c*e
                     tmp21 = tmp20*g1
                     tmp22 = tmp3*e
                     tmp23 = tmp22*g1
                     tmp24 = b*det_g
                     tmp25 = tmp24*g1
                     tmp26 = c*e
                     tmp27 = tmp26*num_b1
                     tmp28 = tmp13*z
                     tmp29 = d*d
                     tmp30 = e*e
                     tmp31 = tmp17*tmp4
                     tmp32 = det_g*e*g1
                     tmp33 = c*d
                     tmp34 = den1*den1
                     tmp35 = d*f
                     tmp36 = tmp12*f
                     tmp37 = tmp10*det_g
                     tmp38 = det_g*det_g
                     tmp39 = d*det_g*e
                     tmp12 = tmp12*e
                     tmp40 = det_g*z
                     fg = (((3.0_ki*(2.0_ki*(tmp15*z-tmp25)+tmp28)+4.0_ki*(2.0_ki*(tmp11*z+tmp2*num_b&
                     &1+tmp5*tmp7+tmp3*f*g1-(a*f*num_b1+tmp10*e*g1+tmp3*g1*g1+tmp9+tmp6+tmp1*tmp1))+t&
                     &mp18+tmp21+tmp23+tmp19*z+tmp25*z-(c*f*num_b1*z+a*e*num_b1*z+tmp27*z+tmp1*num_b1&
                     &+tmp15+tmp13+tmp11))+2.0_ki*(tmp27+tmp14*tmp30+tmp17*tmp7+tmp2*tmp30+tmp29*num_&
                     &b1+tmp33*num_b1*z+d*e*num_b1*z-(tmp33*e*g1+tmp19*z*z+tmp16*num_b1+tmp1*tmp30+tm&
                     &p32+tmp31+tmp19+tmp14*tmp14))+tmp32*z+tmp4*num_b1)/tmp34+4.0_ki*(2.0_ki*(tmp31+&
                     &tmp37+tmp39+tmp10*tmp30+tmp12*z+tmp36*z+a*det_g*e*z-(tmp30*a*d+tmp38*z+tmp30*tm&
                     &p33+tmp17*tmp35+tmp28))+4.0_ki*(2.0_ki*(tmp6+tmp9-(tmp8*f+tmp35*tmp5))+tmp13+tm&
                     &p25+tmp20*f+tmp22*f+tmp17*e*f+tmp29*b*e-(tmp37*z+tmp24*f+tmp36+tmp23+tmp21+tmp1&
                     &8))+tmp38+tmp40*tmp40+tmp26*tmp29+tmp30*det_g-(tmp33*det_g*z+tmp39*z+tmp26*tmp3&
                     &0+tmp12))/(tmp34*den1)*g1*num_b1-1.0_ki/2.0_ki*(8.0_ki*(a+b*z)+4.0_ki*(d-b)+2.0&
                     &_ki*(3.0_ki*c*z+e*z)+g1)/den1)*g1)
                     !
                  case(2)
                     !
                     tmp1 = b*g2
                     tmp2 = tmp1*num_b2
                     tmp3 = c*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = b*f*num_b2
                     tmp6 = c*f*num_b2
                     tmp7 = a*c
                     tmp8 = g2*g2
                     tmp9 = b*c
                     tmp10 = d*d
                     tmp11 = c*c
                     tmp12 = d*e*num_b2
                     tmp13 = det_g*num_b2
                     tmp14 = z*z
                     tmp15 = a*d
                     tmp16 = b*d
                     tmp17 = c*d
                     tmp18 = e*e
                     tmp19 = e*g2
                     tmp20 = num_b2*z
                     tmp21 = c*e
                     tmp22 = d*det_g
                     tmp23 = d*g2
                     tmp24 = c*det_g
                     tmp25 = tmp24*g2
                     tmp26 = den2*den2
                     tmp27 = d*f
                     tmp28 = b*b
                     tmp29 = tmp9*d
                     tmp30 = tmp7*e
                     tmp31 = tmp9*e
                     tmp32 = b*det_g
                     tmp33 = det_g*z
                     fg = (((12.0_ki*(-(tmp4*z+tmp2*z))+8.0_ki*(tmp2+tmp4+tmp5*z+tmp6*z+a*g2*num_b2-(&
                     &a*f*num_b2+tmp6+tmp5))+2.0_ki*(tmp10*num_b2+tmp13*tmp14+tmp18*num_b2+tmp14*tmp2&
                     &1*num_b2+tmp22*g2*z+det_g*e*g2*z-(tmp19*num_b2+tmp18*tmp20))+4.0_ki*(tmp12+tmp1&
                     &*tmp10+tmp11*f*g2+tmp16*e*g2+tmp7*f*g2+tmp9*f*g2-(tmp17*e*g2+tmp15*e*g2+tmp14*t&
                     &mp16*num_b2+tmp13*tmp14*z+tmp8*tmp9+tmp7*tmp8+tmp12*z+tmp3*tmp3))+tmp10*tmp3+tm&
                     &p23*num_b2-(3.0_ki*tmp18*tmp3+tmp14*tmp25))/tmp26+4.0_ki*(8.0_ki*(tmp23*tmp28+t&
                     &mp29*g2-(tmp29*f+tmp27*tmp28))+2.0_ki*(tmp11*tmp23+tmp16*tmp18+tmp14*tmp16*det_&
                     &g+tmp22*e*z-(tmp17*tmp18+tmp15*tmp18+tmp11*tmp27))+4.0_ki*(tmp25*z+tmp30*f+tmp3&
                     &1*f+tmp10*b*e+tmp11*e*f+tmp32*g2*z-(tmp32*f*z+tmp24*f*z+tmp31*g2+tmp30*g2+tmp11&
                     &*tmp19))+tmp10*tmp21+tmp18*tmp33+tmp33*tmp33*z-(tmp14*tmp24*e+tmp18*tmp21))/(tm&
                     &p26*den2)*g2*num_b2)*g2+(2.0_ki*(tmp19+tmp23-(tmp23*z+tmp19*z))+tmp14*tmp20+tmp&
                     &14*tmp3-(tmp14*num_b2+tmp8))/den2)
                     !
                  case(3)
                     !
                     tmp1 = b*c
                     tmp2 = g3*g3
                     tmp3 = b*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = c*g3
                     tmp6 = f*g3
                     tmp7 = c*c
                     tmp8 = c*d
                     tmp9 = tmp5*num_b3
                     tmp10 = a*g3
                     tmp11 = b*f*num_b3
                     tmp12 = b*b
                     tmp13 = c*f*num_b3
                     tmp14 = a*c
                     tmp15 = d*d
                     tmp16 = d*e*num_b3
                     tmp17 = a*d
                     tmp18 = z*z
                     tmp19 = b*d
                     tmp20 = e*e
                     tmp21 = e*g3
                     tmp22 = a*det_g
                     tmp23 = b*det_g
                     tmp24 = tmp23*g3
                     tmp25 = c*det_g
                     tmp26 = tmp25*g3
                     tmp27 = c*e
                     tmp28 = d*g3
                     tmp29 = d*det_g
                     tmp30 = den3*den3
                     tmp31 = d*f
                     tmp32 = tmp1*d
                     tmp33 = tmp14*e
                     tmp34 = tmp1*e
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+4.0_ki*(tmp18*a+tmp18*b+tmp18*c)+3.0_ki*g&
                     &3)/den3+4.0_ki*(2.0_ki*(tmp19*tmp20+tmp26*z+tmp28*tmp7+tmp18*tmp22*e-(tmp25*f*z&
                     &+tmp18*tmp19*det_g+tmp31*tmp7+tmp20*tmp8+tmp17*tmp20))+4.0_ki*(2.0_ki*(tmp12*tm&
                     &p28+tmp32*g3-(tmp32*f+tmp12*tmp31))+tmp24*z+tmp33*f+tmp34*f+tmp15*b*e+tmp7*e*f-&
                     &(tmp23*f*z+tmp34*g3+tmp33*g3+tmp21*tmp7))+tmp15*tmp27+tmp18*tmp25*e+tmp20*det_g&
                     &*z+tmp29*e*z-(tmp18*tmp8*det_g+tmp20*tmp27))/(tmp30*den3)*g3*num_b3-(6.0_ki*(2.&
                     &0_ki*(tmp4*z+tmp1*f*g3-tmp1*tmp2)+tmp6*tmp7+tmp9*z-(tmp8*e*g3+tmp5*tmp5))+4.0_k&
                     &i*(2.0_ki*(tmp11+tmp13+tmp12*tmp6+a*f*num_b3-(tmp11*z+tmp10*num_b3+tmp9+tmp4+tm&
                     &p3*tmp3))+tmp15*tmp3+tmp14*f*g3+tmp18*a*e*num_b3-(tmp19*e*g3+tmp18*tmp19*num_b3&
                     &+tmp17*e*g3+tmp14*tmp2+tmp13*z+tmp16))+2.0_ki*(tmp10*tmp20+tmp16*z+tmp18*tmp24+&
                     &tmp18*tmp26+tmp21*num_b3+tmp18*tmp22*g3+tmp18*tmp27*num_b3+tmp20*num_b3*z-(tmp1&
                     &8*tmp8*num_b3+tmp20*num_b3+tmp20*tmp3+tmp15*num_b3))+tmp15*tmp5+tmp29*g3*z+det_&
                     &g*e*g3*z-(tmp28*num_b3+tmp20*tmp5))/tmp30)*g3)
                     !
                  case(4)
                     !
                     tmp1 = c*g1
                     tmp2 = b*g1
                     tmp3 = e*g1
                     tmp4 = tmp1*z
                     tmp5 = z*z
                     tmp6 = b*f
                     tmp7 = c*f
                     tmp8 = b*d
                     tmp9 = tmp8*z
                     tmp10 = d*e
                     tmp11 = tmp5*det_g
                     tmp12 = a*e*z
                     tmp13 = c*e
                     tmp14 = tmp13*z
                     tmp15 = e*e
                     tmp16 = c*d*z
                     tmp17 = c*g2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(a*z+b*z+c*z)+d+e)/den3*z-(4.0_ki*(b*g3-tmp6)+2.0_k&
                     &i*(tmp12+c*g3-(tmp9+tmp7))+tmp10+tmp14+tmp15-tmp16)/(den3*den3)*num_b3*z)*g3+(2&
                     &.0_ki*(tmp1-tmp2*z)+3.0_ki/2.0_ki*(2.0_ki*tmp2-tmp4)+tmp3+1.0_ki/2.0_ki*num_b1+&
                     &1.0_ki/2.0_ki*tmp5*num_b1-(num_b1*z+1.0_ki/2.0_ki*tmp3*z))/den1+1.0_ki/2.0_ki*(&
                     &2.0_ki*(-(num_b2*z+e*g2+d*g2))+tmp17*z+3.0_ki*tmp5*num_b2)/den2*z-((2.0_ki*(2.0&
                     &_ki*(tmp17+b*g2-(tmp7+tmp6))+tmp10+tmp11+tmp9)+tmp15-tmp14)/(den2*den2)*g2*num_&
                     &b2*z+(2.0_ki*(2.0_ki*(tmp1+tmp2-(det_g*z+tmp9+tmp7+tmp6))+tmp10+tmp11+tmp12+tmp&
                     &14+tmp8+det_g+tmp7*z-tmp4)+tmp15-(tmp10*z+tmp16+tmp13))/(den1*den1)*g1*num_b1))
                     !
                  end select
                  !
               end select
               !
            case(2)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = b*c
                     tmp3 = a*e
                     tmp4 = b*d
                     tmp5 = c*f
                     tmp6 = c*g1
                     tmp7 = e*e
                     tmp8 = c*det_g*g1
                     tmp9 = tmp6*num_b1
                     tmp10 = d*g1
                     tmp11 = c*c
                     tmp12 = e*g1
                     tmp13 = det_g*num_b1
                     tmp14 = c*d
                     tmp15 = c*e
                     tmp16 = d*e
                     tmp17 = det_g*e*g1
                     tmp18 = den1*den1
                     tmp19 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(3.0_ki*c*z+e*z)+4.0_ki*(2.0_ki*a+d+b*z)+g1)/den1-(&
                     &4.0_ki*(2.0_ki*(tmp5-tmp6)+tmp19-tmp16)*(2.0_ki*(tmp3-tmp4)+tmp15+tmp19-(det_g+&
                     &tmp14))/(tmp18*den1)*g1*num_b1+(3.0_ki*(tmp8*z+2.0_ki*tmp9*z)+4.0_ki*(2.0_ki*(t&
                     &mp1*num_b1-a*f*num_b1)+tmp2*f*g1+tmp4*num_b1*z+a*c*e*g1-(tmp5*num_b1*z+tmp4*e*g&
                     &1+tmp3*num_b1*z+tmp2*d*g1+tmp2*g1*g1))+2.0_ki*(tmp1*tmp7+tmp11*tmp12+tmp13*z+d*&
                     &d*num_b1+tmp11*f*g1+tmp14*num_b1*z+tmp16*num_b1*z+b*det_g*g1*z-(tmp15*num_b1*z+&
                     &tmp14*e*g1+tmp13*z*z+tmp10*tmp11+tmp9+tmp8+tmp6*tmp6))+tmp10*num_b1+tmp17*z+tmp&
                     &6*tmp7-(tmp12*num_b1+tmp17))/tmp18))*g1)
                     !
                  case(2)
                     !
                     tmp1 = a*c
                     tmp2 = c*f
                     tmp3 = tmp2*num_b2
                     tmp4 = c*g2
                     tmp5 = tmp4*num_b2
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = tmp7*num_b2
                     tmp9 = c*d
                     tmp10 = den2*den2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z-d)+g2)/den2-(4.0_ki*(2.0_ki*(tmp2-tmp4)-tmp7)*&
                     &(2.0_ki*(a*e-b*d)+c*e-(det_g*z+tmp9))/(tmp10*den2)*g2*num_b2+(4.0_ki*(2.0_ki*(a&
                     &*g2*num_b2-a*f*num_b2)+tmp5+tmp3*z+tmp1*f*g2-(a*d*e*g2+tmp1*g2*g2+tmp3))+2.0_ki&
                     &*(tmp8+tmp6*num_b2+tmp6*b*g2+c*c*f*g2-(tmp9*e*g2+tmp8*z+3.0_ki*tmp5*z+tmp4*tmp4&
                     &))+tmp4*tmp6+d*g2*num_b2+d*det_g*g2*z-e*g2*num_b2)/tmp10))*g2)
                     !
                  case(3)
                     !
                     tmp1 = z*z
                     tmp2 = a*g3
                     tmp3 = c*g3
                     tmp4 = a*c
                     tmp5 = g3*g3
                     tmp6 = b*c
                     tmp7 = c*f
                     tmp8 = a*e
                     tmp9 = b*d
                     tmp10 = c*d
                     tmp11 = d*d
                     tmp12 = e*e
                     tmp13 = d*e
                     tmp14 = c*e
                     tmp15 = den3*den3
                     fg = (((8.0_ki*(a*f*num_b3-tmp2*num_b3)+2.0_ki*(tmp12*tmp2+tmp1*tmp14*num_b3+tmp&
                     &11*b*g3+tmp1*a*det_g*g3+tmp1*b*det_g*g3+tmp1*c*det_g*g3-(tmp1*tmp10*num_b3+tmp1&
                     &3*num_b3+tmp11*num_b3))+4.0_ki*(tmp7*num_b3+tmp1*tmp8*num_b3+tmp4*f*g3+tmp6*f*g&
                     &3+c*c*f*g3-(a*d*e*g3+tmp9*e*g3+tmp10*e*g3+tmp1*tmp9*num_b3+tmp5*tmp6+tmp4*tmp5+&
                     &tmp3*num_b3+tmp3*tmp3))+tmp11*tmp3+tmp12*tmp3+e*g3*num_b3-d*g3*num_b3)/tmp15-(4&
                     &.0_ki*(2.0_ki*(tmp7-tmp3)+tmp1*det_g-tmp13)*(2.0_ki*(tmp8-tmp9)+tmp14-tmp10)/(t&
                     &mp15*den3)*g3*num_b3+(2.0_ki*(tmp1*a+tmp1*b+tmp1*c)+d+e+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = z*z
                     tmp2 = a*e
                     tmp3 = b*d
                     tmp4 = c*d
                     tmp5 = c*e
                     tmp6 = c*f
                     tmp7 = d*e
                     tmp8 = c*g1
                     tmp9 = e*g1
                     tmp10 = tmp8*z
                     fg = (((2.0_ki*(tmp2-tmp3)+tmp5-tmp4)/(den3*den3)*tmp1*num_b3-(a+b+c)*tmp1/den3)&
                     &*g3+(1.0_ki/2.0_ki*d/den2*z-(2.0_ki*(tmp6-c*g2)-tmp7)/(den2*den2)*num_b2*z)*g2+&
                     &(3.0_ki/2.0_ki*tmp10+1.0_ki/2.0_ki*tmp9*z+1.0_ki/2.0_ki*num_b1*z+b*g1*z-(1.0_ki&
                     &/2.0_ki*tmp1*num_b1+1.0_ki/2.0_ki*tmp9+tmp8))/den1+(2.0_ki*(tmp8+tmp1*det_g+tmp&
                     &2*z+tmp6*z-(det_g*z+tmp3*z+tmp6+tmp10))+tmp7+tmp5*z-(tmp7*z+tmp4*z))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = e*g1
                     tmp2 = a*f
                     tmp3 = a*g1
                     tmp4 = c*g1
                     tmp5 = a*c
                     tmp6 = g1*g1
                     tmp7 = det_g*num_b1
                     tmp8 = z*z
                     tmp9 = a*e
                     tmp10 = tmp9*num_b1
                     tmp11 = c*f
                     tmp12 = d*d
                     tmp13 = c*det_g*g1
                     tmp14 = d*g1
                     tmp15 = c*c
                     tmp16 = num_b1*z
                     tmp17 = c*d
                     tmp18 = tmp17*num_b1
                     tmp19 = d*e
                     tmp20 = den1*den1
                     tmp21 = det_g*z
                     fg = ((-(4.0_ki*(tmp21+2.0_ki*tmp9-tmp17)*(2.0_ki*(tmp11-tmp4)+4.0_ki*(tmp2-tmp3&
                     &)+tmp21-(tmp8*det_g+tmp19+tmp12))/(tmp20*den1)*g1*num_b1+(8.0_ki*(tmp1*a*a-tmp2&
                     &*num_b1*z)+3.0_ki*(tmp13*z-tmp12*tmp4)+6.0_ki*(2.0_ki*tmp3*num_b1*z+tmp4*num_b1&
                     &*z)+2.0_ki*(tmp12*tmp16+tmp18*z+tmp3*e*e+tmp15*f*g1+tmp19*num_b1*z+d*det_g*g1*z&
                     &-(tmp17*e*g1+tmp18*tmp8+tmp14*tmp15+tmp4*tmp4))+4.0_ki*(tmp10*tmp8+tmp5*e*g1+tm&
                     &p5*f*g1+tmp7*tmp8*z+a*d*e*g1+a*det_g*g1*z-(tmp5*d*g1+tmp11*num_b1*z+tmp7*tmp8+t&
                     &mp5*tmp6+tmp10*z))+tmp14*num_b1+det_g*e*g1*z-tmp13*tmp8)/tmp20))*g1+(2.0_ki*(tm&
                     &p14*z+2.0_ki*tmp3*z)+1.0_ki/2.0_ki*tmp6+tmp1*z+tmp16*tmp8+3.0_ki*tmp4*z-(tmp8*n&
                     &um_b1+tmp4*tmp8))/den1)
                     !
                  case(2)
                     !
                     tmp1 = a*a
                     tmp2 = a*c
                     tmp3 = d*d
                     tmp4 = den2*den2
                     fg = ((-(4.0_ki*(2.0_ki*a*e-c*d)*(4.0_ki*(a*f-a*g2)+2.0_ki*(c*f-c*g2)-(d*e+tmp3)&
                     &)/(tmp4*den2)*num_b2+(8.0_ki*(tmp1*f-tmp1*g2)+4.0_ki*(tmp2*f-(a*d*e+tmp2*g2))+t&
                     &mp3*c+d*num_b2-2.0_ki*tmp3*a)/tmp4))*g2*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*g3
                     tmp2 = a*c
                     tmp3 = f*g3
                     tmp4 = a*f
                     tmp5 = c*f
                     tmp6 = c*g3
                     tmp7 = d*d
                     tmp8 = c*d
                     tmp9 = d*e
                     tmp10 = den3*den3
                     fg = (((4.0_ki*(2.0_ki*(tmp3*a*a+tmp2*f*g3+tmp4*num_b3*z-(tmp2*g3*g3+tmp1*tmp1))&
                     &+tmp5*num_b3*z)+2.0_ki*(3.0_ki*(-(tmp6*num_b3*z+2.0_ki*tmp1*num_b3*z))+tmp1*e*e&
                     &+tmp3*c*c-(tmp9*num_b3*z+tmp8*e*g3+tmp7*num_b3*z+tmp6*tmp7+tmp1*tmp7+tmp6*tmp6)&
                     &)+d*det_g*g3*z+det_g*e*g3*z-d*g3*num_b3)/tmp10-(4.0_ki*(2.0_ki*(2.0_ki*(tmp4-tm&
                     &p1)+tmp5-tmp6)-(tmp9+tmp7))*(2.0_ki*a*e+det_g*z-tmp8)/(tmp10*den3)*g3*num_b3+1.&
                     &0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = a*f
                     tmp2 = c*f
                     tmp3 = d*d
                     tmp4 = d*e
                     tmp5 = c*g1
                     tmp6 = z*z
                     tmp7 = a*g1
                     tmp8 = a*e
                     tmp9 = c*d
                     fg = (((4.0_ki*(tmp1-a*g3)+2.0_ki*(tmp2-c*g3)-(tmp4+tmp3))/(den3*den3)*num_b3*z-&
                     &1.0_ki/2.0_ki*(d+e)/den3*z)*g3+1.0_ki/2.0_ki*(2.0_ki*(2.0_ki*tmp7+d*g1)+3.0_ki*&
                     &(tmp5+tmp6*num_b1-num_b1*z)+e*g1-tmp5*z)/den1*z-(4.0_ki*(tmp7-tmp1)+2.0_ki*(tmp&
                     &5+tmp6*det_g+tmp8*z-(det_g*z+tmp8+tmp2))+tmp3+tmp4+tmp9-tmp9*z)/(den1*den1)*g1*&
                     &num_b1*z)
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = c*g1
                     tmp3 = tmp2*num_b1
                     tmp4 = e*g1
                     tmp5 = a*a
                     tmp6 = tmp4*tmp5
                     tmp7 = det_g*num_b1
                     tmp8 = z*z
                     tmp9 = a*c
                     tmp10 = tmp9*e
                     tmp11 = tmp10*g1
                     tmp12 = a*e
                     tmp13 = tmp12*num_b1
                     tmp14 = c*f*num_b1
                     tmp15 = c*det_g
                     tmp16 = tmp15*g1
                     tmp17 = tmp16*z
                     tmp18 = g1*g1
                     tmp19 = e*e
                     tmp20 = b*c
                     tmp21 = d*g1
                     tmp22 = c*c
                     tmp23 = tmp21*tmp22
                     tmp24 = tmp7*z
                     tmp25 = tmp9*d
                     tmp26 = tmp25*g1
                     tmp27 = a*d
                     tmp28 = a*det_g
                     tmp29 = tmp28*g1*z
                     tmp30 = tmp20*d
                     tmp31 = tmp30*g1
                     tmp32 = b*d
                     tmp33 = c*d
                     tmp34 = tmp33*num_b1
                     tmp35 = d*e*num_b1
                     tmp36 = d*d
                     tmp37 = tmp22*tmp4
                     tmp38 = num_b1*z
                     tmp39 = c*e
                     tmp40 = d*det_g
                     tmp41 = det_g*e*g1
                     tmp42 = den1*den1
                     tmp43 = e*f
                     tmp44 = tmp28*e
                     tmp45 = tmp15*f
                     tmp46 = det_g*z
                     tmp47 = tmp46*tmp46
                     tmp48 = tmp33*det_g
                     tmp49 = tmp40*e
                     fg = (((6.0_ki*(tmp17-tmp3)+12.0_ki*(tmp3*z+tmp1*num_b1*z)+8.0_ki*(tmp11+tmp6-(a&
                     &*f*num_b1*z+tmp7*tmp8+tmp14*z+tmp13*z))+2.0_ki*(tmp37+tmp21*num_b1+tmp36*tmp38+&
                     &tmp41*z+tmp40*g1*z+b*det_g*g1*z-(tmp39*num_b1*z+tmp34*tmp8+tmp35+tmp16))+4.0_ki&
                     &*(tmp14+tmp24+tmp29+tmp1*tmp19+tmp13*tmp8+tmp24*tmp8+tmp34*z+tmp35*z+tmp20*f*g1&
                     &+tmp22*f*g1+tmp27*e*g1+tmp32*num_b1*z+tmp9*f*g1-(tmp33*e*g1+tmp32*e*g1+tmp18*tm&
                     &p9+tmp18*tmp20+tmp31+tmp26+tmp23+tmp2*tmp2))+tmp19*tmp2-(tmp4*num_b1+3.0_ki*tmp&
                     &2*tmp36+tmp16*tmp8+tmp41))/tmp42+4.0_ki*(8.0_ki*(tmp10*f+tmp43*tmp5-(tmp6+tmp11&
                     &))+2.0_ki*(tmp16+tmp47+tmp22*tmp43+tmp36*tmp39+tmp36*b*e-(tmp32*det_g*z+tmp49*z&
                     &+tmp48*z+tmp44*tmp8+tmp12*tmp36+tmp45+tmp37))+4.0_ki*(tmp23+tmp26+tmp31+tmp44*z&
                     &+tmp45*z+tmp28*f*z-(tmp22*d*f+tmp30*f+tmp25*f+tmp19*tmp27+tmp29+tmp17))+tmp49+t&
                     &mp33*tmp36+tmp48*tmp8+tmp15*e*z-(det_g*det_g*z+tmp47*z+tmp36*tmp46+tmp19*tmp33)&
                     &)/(tmp42*den1)*g1*num_b1)*g1+(2.0_ki*(tmp2+tmp8*num_b1-(b*g1*z+tmp4*z+tmp21*z+3&
                     &.0_ki*tmp2*z+2.0_ki*tmp1*z))+tmp4+tmp2*tmp8-(tmp38*tmp8+tmp38+tmp18))/den1)
                     !
                  case(2)
                     !
                     tmp1 = a*g2
                     tmp2 = a*c
                     tmp3 = f*g2
                     tmp4 = a*a
                     tmp5 = a*d
                     tmp6 = c*g2
                     tmp7 = d*d
                     tmp8 = c*c
                     tmp9 = d*g2
                     tmp10 = c*d
                     tmp11 = e*g2
                     tmp12 = d*det_g
                     tmp13 = den2*den2
                     tmp14 = e*f
                     tmp15 = tmp2*e
                     tmp16 = e*e
                     tmp17 = tmp2*d
                     tmp18 = b*c*d
                     tmp19 = c*det_g
                     fg = (((4.0_ki*(2.0_ki*(tmp3*tmp4+tmp2*f*g2-(tmp5*e*g2+tmp2*g2*g2+tmp1*tmp1))+c*&
                     &f*num_b2*z)+2.0_ki*(tmp3*tmp8+tmp6*tmp7+tmp9*num_b2+tmp7*b*g2-(d*e*num_b2*z+3.0&
                     &_ki*tmp6*num_b2*z+tmp10*e*g2+tmp1*tmp7+tmp6*tmp6))+tmp12*g2*z-tmp11*num_b2)/tmp&
                     &13+4.0_ki*(2.0_ki*(tmp14*tmp8+tmp19*g2*z+tmp7*b*e+tmp7*c*e-(tmp7*a*e+tmp19*f*z+&
                     &tmp11*tmp8))+4.0_ki*(2.0_ki*(tmp14*tmp4+tmp15*f-(tmp15*g2+tmp11*tmp4))+tmp17*g2&
                     &+tmp18*g2+tmp8*tmp9-(tmp8*d*f+tmp18*f+tmp17*f+tmp16*tmp5))+tmp10*tmp7+tmp12*e*z&
                     &-tmp10*tmp16)/(tmp13*den2)*g2*num_b2-1.0_ki/2.0_ki*(g2+2.0_ki*d*z)/den2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*c
                     tmp2 = g3*g3
                     tmp3 = a*g3
                     tmp4 = c*g3
                     tmp5 = f*g3
                     tmp6 = c*c
                     tmp7 = c*d
                     tmp8 = a*a
                     tmp9 = e*e
                     tmp10 = b*c
                     tmp11 = a*d
                     tmp12 = a*e
                     tmp13 = z*z
                     tmp14 = b*d
                     tmp15 = d*d
                     tmp16 = d*g3
                     tmp17 = a*det_g
                     tmp18 = tmp17*g3
                     tmp19 = c*det_g
                     tmp20 = tmp19*g3
                     tmp21 = c*e
                     tmp22 = e*g3
                     tmp23 = d*det_g
                     tmp24 = den3*den3
                     tmp25 = e*f
                     tmp26 = tmp1*e
                     tmp27 = tmp1*d
                     tmp28 = tmp10*d
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+4.0_ki*(tmp13*a+tmp13*b+tmp13*c)+3.0_ki*g&
                     &3)/den3+4.0_ki*(2.0_ki*(tmp15*tmp21+tmp25*tmp6+tmp13*tmp17*e+tmp15*b*e+tmp19*f*&
                     &z-(tmp13*tmp14*det_g+tmp22*tmp6+tmp20*z+tmp12*tmp15))+4.0_ki*(2.0_ki*(tmp25*tmp&
                     &8+tmp26*f-(tmp26*g3+tmp22*tmp8))+tmp16*tmp6+tmp27*g3+tmp28*g3+tmp17*f*z-(tmp6*d&
                     &*f+tmp28*f+tmp27*f+tmp18*z+tmp11*tmp9))+tmp15*tmp7+tmp13*tmp19*e-(tmp23*e*z+tmp&
                     &15*det_g*z+tmp13*tmp7*det_g+tmp7*tmp9))/(tmp24*den3)*g3*num_b3-(6.0_ki*(2.0_ki*&
                     &(tmp1*f*g3-(tmp3*num_b3*z+tmp1*tmp2))+tmp5*tmp6-(tmp7*e*g3+tmp4*num_b3*z+tmp4*t&
                     &mp4))+4.0_ki*(2.0_ki*(tmp5*tmp8+a*f*num_b3*z-tmp3*tmp3)+tmp3*tmp9+tmp10*f*g3+tm&
                     &p12*tmp13*num_b3+c*f*num_b3*z-(tmp14*e*g3+tmp13*tmp14*num_b3+tmp11*e*g3+tmp10*t&
                     &mp2))+2.0_ki*(tmp13*tmp18+tmp13*tmp20+tmp13*tmp21*num_b3+tmp15*b*g3+tmp13*b*det&
                     &_g*g3-(d*e*num_b3*z+tmp15*num_b3*z+tmp13*tmp7*num_b3+tmp16*num_b3+tmp15*tmp3))+&
                     &tmp22*num_b3+tmp4*tmp9+tmp23*g3*z+det_g*e*g3*z-tmp15*tmp4)/tmp24)*g3)
                     !
                  case(4)
                     !
                     tmp1 = a*f
                     tmp2 = c*f
                     tmp3 = a*e
                     tmp4 = tmp3*z
                     tmp5 = b*d*z
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = c*d
                     tmp9 = tmp8*z
                     tmp10 = c*e*z
                     tmp11 = z*z
                     tmp12 = c*g1
                     tmp13 = tmp12*z
                     tmp14 = num_b1*z
                     tmp15 = a*g1*z
                     tmp16 = e*g1
                     tmp17 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(a*z+b*z+c*z)+d+e)/den3*z-(4.0_ki*(tmp1-a*g3)+2.0_k&
                     &i*(tmp2+tmp4-(c*g3+tmp5))+tmp10-(tmp9+tmp7+tmp6))/(den3*den3)*num_b3*z)*g3+((2.&
                     &0_ki*(tmp2-c*g2)-tmp7)/(den2*den2)*num_b2*z-1.0_ki/2.0_ki*d/den2*z)*g2+(3.0_ki*&
                     &(tmp11*num_b1-tmp13)+3.0_ki/2.0_ki*(-(tmp11*tmp14+tmp14))+tmp12+1.0_ki/2.0_ki*t&
                     &mp16+1.0_ki/2.0_ki*tmp11*tmp12-(d*g1*z+b*g1*z+tmp16*z+2.0_ki*tmp15))/den1+(4.0_&
                     &ki*(tmp13+tmp15-(tmp2*z+tmp11*det_g+tmp1*z+tmp4))+2.0_ki*(tmp17+tmp2+tmp5+tmp9+&
                     &tmp11*tmp17+tmp11*tmp3+tmp7*z-tmp12)+tmp6*z-(tmp11*tmp8+tmp7+tmp10))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               end select
               !
            case(3)
               !
               select case(par3)
               !
               case(1)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = b*g1
                     tmp2 = a*g1
                     tmp3 = b*c
                     tmp4 = d*g1
                     tmp5 = b*b
                     tmp6 = tmp4*tmp5
                     tmp7 = f*g1
                     tmp8 = tmp3*d
                     tmp9 = tmp8*g1
                     tmp10 = b*d
                     tmp11 = tmp10*num_b1
                     tmp12 = c*det_g
                     tmp13 = tmp12*g1
                     tmp14 = c*g1
                     tmp15 = tmp14*num_b1
                     tmp16 = e*g1
                     tmp17 = c*c
                     tmp18 = tmp16*tmp17
                     tmp19 = det_g*num_b1
                     tmp20 = a*c*e
                     tmp21 = tmp20*g1
                     tmp22 = tmp3*e
                     tmp23 = tmp22*g1
                     tmp24 = b*det_g
                     tmp25 = tmp24*g1
                     tmp26 = c*e
                     tmp27 = tmp26*num_b1
                     tmp28 = tmp13*z
                     tmp29 = d*d
                     tmp30 = e*e
                     tmp31 = tmp17*tmp4
                     tmp32 = det_g*e*g1
                     tmp33 = c*d
                     tmp34 = den1*den1
                     tmp35 = d*f
                     tmp36 = tmp12*f
                     tmp37 = tmp10*det_g
                     tmp38 = det_g*det_g
                     tmp39 = d*det_g*e
                     tmp12 = tmp12*e
                     tmp40 = det_g*z
                     fg = (((3.0_ki*(2.0_ki*(tmp15*z-tmp25)+tmp28)+4.0_ki*(2.0_ki*(tmp11*z+tmp2*num_b&
                     &1+tmp5*tmp7+tmp3*f*g1-(a*f*num_b1+tmp10*e*g1+tmp3*g1*g1+tmp9+tmp6+tmp1*tmp1))+t&
                     &mp18+tmp21+tmp23+tmp19*z+tmp25*z-(c*f*num_b1*z+a*e*num_b1*z+tmp27*z+tmp1*num_b1&
                     &+tmp15+tmp13+tmp11))+2.0_ki*(tmp27+tmp14*tmp30+tmp17*tmp7+tmp2*tmp30+tmp29*num_&
                     &b1+tmp33*num_b1*z+d*e*num_b1*z-(tmp33*e*g1+tmp19*z*z+tmp16*num_b1+tmp1*tmp30+tm&
                     &p32+tmp31+tmp19+tmp14*tmp14))+tmp32*z+tmp4*num_b1)/tmp34+4.0_ki*(2.0_ki*(tmp31+&
                     &tmp37+tmp39+tmp10*tmp30+tmp12*z+tmp36*z+a*det_g*e*z-(tmp30*a*d+tmp38*z+tmp30*tm&
                     &p33+tmp17*tmp35+tmp28))+4.0_ki*(2.0_ki*(tmp6+tmp9-(tmp8*f+tmp35*tmp5))+tmp13+tm&
                     &p25+tmp20*f+tmp22*f+tmp17*e*f+tmp29*b*e-(tmp37*z+tmp24*f+tmp36+tmp23+tmp21+tmp1&
                     &8))+tmp38+tmp40*tmp40+tmp26*tmp29+tmp30*det_g-(tmp33*det_g*z+tmp39*z+tmp26*tmp3&
                     &0+tmp12))/(tmp34*den1)*g1*num_b1-1.0_ki/2.0_ki*(8.0_ki*(a+b*z)+4.0_ki*(d-b)+2.0&
                     &_ki*(3.0_ki*c*z+e*z)+g1)/den1)*g1)
                     !
                  case(2)
                     !
                     tmp1 = b*g2
                     tmp2 = tmp1*num_b2
                     tmp3 = c*g2
                     tmp4 = tmp3*num_b2
                     tmp5 = b*f*num_b2
                     tmp6 = c*f*num_b2
                     tmp7 = a*c
                     tmp8 = g2*g2
                     tmp9 = b*c
                     tmp10 = d*d
                     tmp11 = c*c
                     tmp12 = d*e*num_b2
                     tmp13 = det_g*num_b2
                     tmp14 = z*z
                     tmp15 = a*d
                     tmp16 = b*d
                     tmp17 = c*d
                     tmp18 = e*e
                     tmp19 = e*g2
                     tmp20 = num_b2*z
                     tmp21 = c*e
                     tmp22 = d*det_g
                     tmp23 = d*g2
                     tmp24 = c*det_g
                     tmp25 = tmp24*g2
                     tmp26 = den2*den2
                     tmp27 = d*f
                     tmp28 = b*b
                     tmp29 = tmp9*d
                     tmp30 = tmp7*e
                     tmp31 = tmp9*e
                     tmp32 = b*det_g
                     tmp33 = det_g*z
                     fg = (((12.0_ki*(-(tmp4*z+tmp2*z))+8.0_ki*(tmp2+tmp4+tmp5*z+tmp6*z+a*g2*num_b2-(&
                     &a*f*num_b2+tmp6+tmp5))+2.0_ki*(tmp10*num_b2+tmp13*tmp14+tmp18*num_b2+tmp14*tmp2&
                     &1*num_b2+tmp22*g2*z+det_g*e*g2*z-(tmp19*num_b2+tmp18*tmp20))+4.0_ki*(tmp12+tmp1&
                     &*tmp10+tmp11*f*g2+tmp16*e*g2+tmp7*f*g2+tmp9*f*g2-(tmp17*e*g2+tmp15*e*g2+tmp14*t&
                     &mp16*num_b2+tmp13*tmp14*z+tmp8*tmp9+tmp7*tmp8+tmp12*z+tmp3*tmp3))+tmp10*tmp3+tm&
                     &p23*num_b2-(3.0_ki*tmp18*tmp3+tmp14*tmp25))/tmp26+4.0_ki*(8.0_ki*(tmp23*tmp28+t&
                     &mp29*g2-(tmp29*f+tmp27*tmp28))+2.0_ki*(tmp11*tmp23+tmp16*tmp18+tmp14*tmp16*det_&
                     &g+tmp22*e*z-(tmp17*tmp18+tmp15*tmp18+tmp11*tmp27))+4.0_ki*(tmp25*z+tmp30*f+tmp3&
                     &1*f+tmp10*b*e+tmp11*e*f+tmp32*g2*z-(tmp32*f*z+tmp24*f*z+tmp31*g2+tmp30*g2+tmp11&
                     &*tmp19))+tmp10*tmp21+tmp18*tmp33+tmp33*tmp33*z-(tmp14*tmp24*e+tmp18*tmp21))/(tm&
                     &p26*den2)*g2*num_b2)*g2+(2.0_ki*(tmp19+tmp23-(tmp23*z+tmp19*z))+tmp14*tmp20+tmp&
                     &14*tmp3-(tmp14*num_b2+tmp8))/den2)
                     !
                  case(3)
                     !
                     tmp1 = b*c
                     tmp2 = g3*g3
                     tmp3 = b*g3
                     tmp4 = tmp3*num_b3
                     tmp5 = c*g3
                     tmp6 = f*g3
                     tmp7 = c*c
                     tmp8 = c*d
                     tmp9 = tmp5*num_b3
                     tmp10 = a*g3
                     tmp11 = b*f*num_b3
                     tmp12 = b*b
                     tmp13 = c*f*num_b3
                     tmp14 = a*c
                     tmp15 = d*d
                     tmp16 = d*e*num_b3
                     tmp17 = a*d
                     tmp18 = z*z
                     tmp19 = b*d
                     tmp20 = e*e
                     tmp21 = e*g3
                     tmp22 = a*det_g
                     tmp23 = b*det_g
                     tmp24 = tmp23*g3
                     tmp25 = c*det_g
                     tmp26 = tmp25*g3
                     tmp27 = c*e
                     tmp28 = d*g3
                     tmp29 = d*det_g
                     tmp30 = den3*den3
                     tmp31 = d*f
                     tmp32 = tmp1*d
                     tmp33 = tmp14*e
                     tmp34 = tmp1*e
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+4.0_ki*(tmp18*a+tmp18*b+tmp18*c)+3.0_ki*g&
                     &3)/den3+4.0_ki*(2.0_ki*(tmp19*tmp20+tmp26*z+tmp28*tmp7+tmp18*tmp22*e-(tmp25*f*z&
                     &+tmp18*tmp19*det_g+tmp31*tmp7+tmp20*tmp8+tmp17*tmp20))+4.0_ki*(2.0_ki*(tmp12*tm&
                     &p28+tmp32*g3-(tmp32*f+tmp12*tmp31))+tmp24*z+tmp33*f+tmp34*f+tmp15*b*e+tmp7*e*f-&
                     &(tmp23*f*z+tmp34*g3+tmp33*g3+tmp21*tmp7))+tmp15*tmp27+tmp18*tmp25*e+tmp20*det_g&
                     &*z+tmp29*e*z-(tmp18*tmp8*det_g+tmp20*tmp27))/(tmp30*den3)*g3*num_b3-(6.0_ki*(2.&
                     &0_ki*(tmp4*z+tmp1*f*g3-tmp1*tmp2)+tmp6*tmp7+tmp9*z-(tmp8*e*g3+tmp5*tmp5))+4.0_k&
                     &i*(2.0_ki*(tmp11+tmp13+tmp12*tmp6+a*f*num_b3-(tmp11*z+tmp10*num_b3+tmp9+tmp4+tm&
                     &p3*tmp3))+tmp15*tmp3+tmp14*f*g3+tmp18*a*e*num_b3-(tmp19*e*g3+tmp18*tmp19*num_b3&
                     &+tmp17*e*g3+tmp14*tmp2+tmp13*z+tmp16))+2.0_ki*(tmp10*tmp20+tmp16*z+tmp18*tmp24+&
                     &tmp18*tmp26+tmp21*num_b3+tmp18*tmp22*g3+tmp18*tmp27*num_b3+tmp20*num_b3*z-(tmp1&
                     &8*tmp8*num_b3+tmp20*num_b3+tmp20*tmp3+tmp15*num_b3))+tmp15*tmp5+tmp29*g3*z+det_&
                     &g*e*g3*z-(tmp28*num_b3+tmp20*tmp5))/tmp30)*g3)
                     !
                  case(4)
                     !
                     tmp1 = c*g1
                     tmp2 = b*g1
                     tmp3 = e*g1
                     tmp4 = tmp1*z
                     tmp5 = z*z
                     tmp6 = b*f
                     tmp7 = c*f
                     tmp8 = b*d
                     tmp9 = tmp8*z
                     tmp10 = d*e
                     tmp11 = tmp5*det_g
                     tmp12 = a*e*z
                     tmp13 = c*e
                     tmp14 = tmp13*z
                     tmp15 = e*e
                     tmp16 = c*d*z
                     tmp17 = c*g2
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(a*z+b*z+c*z)+d+e)/den3*z-(4.0_ki*(b*g3-tmp6)+2.0_k&
                     &i*(tmp12+c*g3-(tmp9+tmp7))+tmp10+tmp14+tmp15-tmp16)/(den3*den3)*num_b3*z)*g3+(2&
                     &.0_ki*(tmp1-tmp2*z)+3.0_ki/2.0_ki*(2.0_ki*tmp2-tmp4)+tmp3+1.0_ki/2.0_ki*num_b1+&
                     &1.0_ki/2.0_ki*tmp5*num_b1-(num_b1*z+1.0_ki/2.0_ki*tmp3*z))/den1+1.0_ki/2.0_ki*(&
                     &2.0_ki*(-(num_b2*z+e*g2+d*g2))+tmp17*z+3.0_ki*tmp5*num_b2)/den2*z-((2.0_ki*(2.0&
                     &_ki*(tmp17+b*g2-(tmp7+tmp6))+tmp10+tmp11+tmp9)+tmp15-tmp14)/(den2*den2)*g2*num_&
                     &b2*z+(2.0_ki*(2.0_ki*(tmp1+tmp2-(det_g*z+tmp9+tmp7+tmp6))+tmp10+tmp11+tmp12+tmp&
                     &14+tmp8+det_g+tmp7*z-tmp4)+tmp15-(tmp10*z+tmp16+tmp13))/(den1*den1)*g1*num_b1))
                     !
                  end select
                  !
               case(2)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = a*g1
                     tmp2 = c*g1
                     tmp3 = tmp2*num_b1
                     tmp4 = e*g1
                     tmp5 = a*a
                     tmp6 = tmp4*tmp5
                     tmp7 = det_g*num_b1
                     tmp8 = z*z
                     tmp9 = a*c
                     tmp10 = tmp9*e
                     tmp11 = tmp10*g1
                     tmp12 = a*e
                     tmp13 = tmp12*num_b1
                     tmp14 = c*f*num_b1
                     tmp15 = c*det_g
                     tmp16 = tmp15*g1
                     tmp17 = tmp16*z
                     tmp18 = g1*g1
                     tmp19 = e*e
                     tmp20 = b*c
                     tmp21 = d*g1
                     tmp22 = c*c
                     tmp23 = tmp21*tmp22
                     tmp24 = tmp7*z
                     tmp25 = tmp9*d
                     tmp26 = tmp25*g1
                     tmp27 = a*d
                     tmp28 = a*det_g
                     tmp29 = tmp28*g1*z
                     tmp30 = tmp20*d
                     tmp31 = tmp30*g1
                     tmp32 = b*d
                     tmp33 = c*d
                     tmp34 = tmp33*num_b1
                     tmp35 = d*e*num_b1
                     tmp36 = d*d
                     tmp37 = tmp22*tmp4
                     tmp38 = num_b1*z
                     tmp39 = c*e
                     tmp40 = d*det_g
                     tmp41 = det_g*e*g1
                     tmp42 = den1*den1
                     tmp43 = e*f
                     tmp44 = tmp28*e
                     tmp45 = tmp15*f
                     tmp46 = det_g*z
                     tmp47 = tmp46*tmp46
                     tmp48 = tmp33*det_g
                     tmp49 = tmp40*e
                     fg = (((6.0_ki*(tmp17-tmp3)+12.0_ki*(tmp3*z+tmp1*num_b1*z)+8.0_ki*(tmp11+tmp6-(a&
                     &*f*num_b1*z+tmp7*tmp8+tmp14*z+tmp13*z))+2.0_ki*(tmp37+tmp21*num_b1+tmp36*tmp38+&
                     &tmp41*z+tmp40*g1*z+b*det_g*g1*z-(tmp39*num_b1*z+tmp34*tmp8+tmp35+tmp16))+4.0_ki&
                     &*(tmp14+tmp24+tmp29+tmp1*tmp19+tmp13*tmp8+tmp24*tmp8+tmp34*z+tmp35*z+tmp20*f*g1&
                     &+tmp22*f*g1+tmp27*e*g1+tmp32*num_b1*z+tmp9*f*g1-(tmp33*e*g1+tmp32*e*g1+tmp18*tm&
                     &p9+tmp18*tmp20+tmp31+tmp26+tmp23+tmp2*tmp2))+tmp19*tmp2-(tmp4*num_b1+3.0_ki*tmp&
                     &2*tmp36+tmp16*tmp8+tmp41))/tmp42+4.0_ki*(8.0_ki*(tmp10*f+tmp43*tmp5-(tmp6+tmp11&
                     &))+2.0_ki*(tmp16+tmp47+tmp22*tmp43+tmp36*tmp39+tmp36*b*e-(tmp32*det_g*z+tmp49*z&
                     &+tmp48*z+tmp44*tmp8+tmp12*tmp36+tmp45+tmp37))+4.0_ki*(tmp23+tmp26+tmp31+tmp44*z&
                     &+tmp45*z+tmp28*f*z-(tmp22*d*f+tmp30*f+tmp25*f+tmp19*tmp27+tmp29+tmp17))+tmp49+t&
                     &mp33*tmp36+tmp48*tmp8+tmp15*e*z-(det_g*det_g*z+tmp47*z+tmp36*tmp46+tmp19*tmp33)&
                     &)/(tmp42*den1)*g1*num_b1)*g1+(2.0_ki*(tmp2+tmp8*num_b1-(b*g1*z+tmp4*z+tmp21*z+3&
                     &.0_ki*tmp2*z+2.0_ki*tmp1*z))+tmp4+tmp2*tmp8-(tmp38*tmp8+tmp38+tmp18))/den1)
                     !
                  case(2)
                     !
                     tmp1 = a*g2
                     tmp2 = a*c
                     tmp3 = f*g2
                     tmp4 = a*a
                     tmp5 = a*d
                     tmp6 = c*g2
                     tmp7 = d*d
                     tmp8 = c*c
                     tmp9 = d*g2
                     tmp10 = c*d
                     tmp11 = e*g2
                     tmp12 = d*det_g
                     tmp13 = den2*den2
                     tmp14 = e*f
                     tmp15 = tmp2*e
                     tmp16 = e*e
                     tmp17 = tmp2*d
                     tmp18 = b*c*d
                     tmp19 = c*det_g
                     fg = (((4.0_ki*(2.0_ki*(tmp3*tmp4+tmp2*f*g2-(tmp5*e*g2+tmp2*g2*g2+tmp1*tmp1))+c*&
                     &f*num_b2*z)+2.0_ki*(tmp3*tmp8+tmp6*tmp7+tmp9*num_b2+tmp7*b*g2-(d*e*num_b2*z+3.0&
                     &_ki*tmp6*num_b2*z+tmp10*e*g2+tmp1*tmp7+tmp6*tmp6))+tmp12*g2*z-tmp11*num_b2)/tmp&
                     &13+4.0_ki*(2.0_ki*(tmp14*tmp8+tmp19*g2*z+tmp7*b*e+tmp7*c*e-(tmp7*a*e+tmp19*f*z+&
                     &tmp11*tmp8))+4.0_ki*(2.0_ki*(tmp14*tmp4+tmp15*f-(tmp15*g2+tmp11*tmp4))+tmp17*g2&
                     &+tmp18*g2+tmp8*tmp9-(tmp8*d*f+tmp18*f+tmp17*f+tmp16*tmp5))+tmp10*tmp7+tmp12*e*z&
                     &-tmp10*tmp16)/(tmp13*den2)*g2*num_b2-1.0_ki/2.0_ki*(g2+2.0_ki*d*z)/den2)*g2)
                     !
                  case(3)
                     !
                     tmp1 = a*c
                     tmp2 = g3*g3
                     tmp3 = a*g3
                     tmp4 = c*g3
                     tmp5 = f*g3
                     tmp6 = c*c
                     tmp7 = c*d
                     tmp8 = a*a
                     tmp9 = e*e
                     tmp10 = b*c
                     tmp11 = a*d
                     tmp12 = a*e
                     tmp13 = z*z
                     tmp14 = b*d
                     tmp15 = d*d
                     tmp16 = d*g3
                     tmp17 = a*det_g
                     tmp18 = tmp17*g3
                     tmp19 = c*det_g
                     tmp20 = tmp19*g3
                     tmp21 = c*e
                     tmp22 = e*g3
                     tmp23 = d*det_g
                     tmp24 = den3*den3
                     tmp25 = e*f
                     tmp26 = tmp1*e
                     tmp27 = tmp1*d
                     tmp28 = tmp10*d
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(d*z+e*z)+4.0_ki*(tmp13*a+tmp13*b+tmp13*c)+3.0_ki*g&
                     &3)/den3+4.0_ki*(2.0_ki*(tmp15*tmp21+tmp25*tmp6+tmp13*tmp17*e+tmp15*b*e+tmp19*f*&
                     &z-(tmp13*tmp14*det_g+tmp22*tmp6+tmp20*z+tmp12*tmp15))+4.0_ki*(2.0_ki*(tmp25*tmp&
                     &8+tmp26*f-(tmp26*g3+tmp22*tmp8))+tmp16*tmp6+tmp27*g3+tmp28*g3+tmp17*f*z-(tmp6*d&
                     &*f+tmp28*f+tmp27*f+tmp18*z+tmp11*tmp9))+tmp15*tmp7+tmp13*tmp19*e-(tmp23*e*z+tmp&
                     &15*det_g*z+tmp13*tmp7*det_g+tmp7*tmp9))/(tmp24*den3)*g3*num_b3-(6.0_ki*(2.0_ki*&
                     &(tmp1*f*g3-(tmp3*num_b3*z+tmp1*tmp2))+tmp5*tmp6-(tmp7*e*g3+tmp4*num_b3*z+tmp4*t&
                     &mp4))+4.0_ki*(2.0_ki*(tmp5*tmp8+a*f*num_b3*z-tmp3*tmp3)+tmp3*tmp9+tmp10*f*g3+tm&
                     &p12*tmp13*num_b3+c*f*num_b3*z-(tmp14*e*g3+tmp13*tmp14*num_b3+tmp11*e*g3+tmp10*t&
                     &mp2))+2.0_ki*(tmp13*tmp18+tmp13*tmp20+tmp13*tmp21*num_b3+tmp15*b*g3+tmp13*b*det&
                     &_g*g3-(d*e*num_b3*z+tmp15*num_b3*z+tmp13*tmp7*num_b3+tmp16*num_b3+tmp15*tmp3))+&
                     &tmp22*num_b3+tmp4*tmp9+tmp23*g3*z+det_g*e*g3*z-tmp15*tmp4)/tmp24)*g3)
                     !
                  case(4)
                     !
                     tmp1 = a*f
                     tmp2 = c*f
                     tmp3 = a*e
                     tmp4 = tmp3*z
                     tmp5 = b*d*z
                     tmp6 = d*d
                     tmp7 = d*e
                     tmp8 = c*d
                     tmp9 = tmp8*z
                     tmp10 = c*e*z
                     tmp11 = z*z
                     tmp12 = c*g1
                     tmp13 = tmp12*z
                     tmp14 = num_b1*z
                     tmp15 = a*g1*z
                     tmp16 = e*g1
                     tmp17 = det_g*z
                     fg = ((1.0_ki/2.0_ki*(2.0_ki*(a*z+b*z+c*z)+d+e)/den3*z-(4.0_ki*(tmp1-a*g3)+2.0_k&
                     &i*(tmp2+tmp4-(c*g3+tmp5))+tmp10-(tmp9+tmp7+tmp6))/(den3*den3)*num_b3*z)*g3+((2.&
                     &0_ki*(tmp2-c*g2)-tmp7)/(den2*den2)*num_b2*z-1.0_ki/2.0_ki*d/den2*z)*g2+(3.0_ki*&
                     &(tmp11*num_b1-tmp13)+3.0_ki/2.0_ki*(-(tmp11*tmp14+tmp14))+tmp12+1.0_ki/2.0_ki*t&
                     &mp16+1.0_ki/2.0_ki*tmp11*tmp12-(d*g1*z+b*g1*z+tmp16*z+2.0_ki*tmp15))/den1+(4.0_&
                     &ki*(tmp13+tmp15-(tmp2*z+tmp11*det_g+tmp1*z+tmp4))+2.0_ki*(tmp17+tmp2+tmp5+tmp9+&
                     &tmp11*tmp17+tmp11*tmp3+tmp7*z-tmp12)+tmp6*z-(tmp11*tmp8+tmp7+tmp10))/(den1*den1&
                     &)*g1*num_b1)
                     !
                  end select
                  !
               case(3)
                  !
                  select case(flag)
                  !
                  case(1)
                     !
                     tmp1 = c*g1
                     tmp2 = tmp1*num_b1
                     tmp3 = b*c
                     tmp4 = g1*g1
                     tmp5 = b*g1
                     tmp6 = c*f*num_b1
                     tmp7 = det_g*num_b1
                     tmp8 = tmp7*z
                     tmp9 = z*z
                     tmp10 = a*c
                     tmp11 = tmp10*e
                     tmp12 = tmp11*g1
                     tmp13 = a*e
                     tmp14 = tmp13*num_b1
                     tmp15 = a*g1
                     tmp16 = tmp3*d
                     tmp17 = tmp16*g1
                     tmp18 = b*d
                     tmp19 = tmp18*num_b1
                     tmp20 = c*det_g
                     tmp21 = tmp20*g1
                     tmp22 = tmp21*z
                     tmp23 = e*g1
                     tmp24 = a*a
                     tmp25 = tmp23*tmp24
                     tmp26 = d*g1
                     tmp27 = b*b
                     tmp28 = tmp26*tmp27
                     tmp29 = f*g1
                     tmp30 = e*e
                     tmp31 = b*det_g
                     tmp32 = tmp31*g1
                     tmp33 = c*c
                     tmp34 = tmp26*tmp33
                     tmp35 = tmp23*tmp33
                     tmp36 = d*e*num_b1
                     tmp37 = c*d
                     tmp38 = tmp37*num_b1
                     tmp39 = c*e
                     tmp40 = tmp39*num_b1
                     tmp41 = tmp10*d
                     tmp42 = tmp41*g1
                     tmp43 = a*d
                     tmp44 = a*det_g
                     tmp45 = tmp44*g1*z
                     tmp46 = tmp3*e
                     tmp47 = tmp46*g1
                     tmp48 = d*d
                     tmp49 = det_g*e*g1
                     tmp50 = num_b1*z
                     tmp51 = d*det_g
                     tmp52 = den1*den1
                     tmp53 = e*f
                     tmp54 = d*f
                     tmp55 = tmp20*f
                     tmp56 = tmp44*e
                     tmp57 = tmp18*det_g
                     tmp58 = det_g*det_g
                     tmp59 = det_g*z
                     tmp60 = tmp59*tmp59
                     tmp61 = tmp51*e
                     tmp62 = tmp37*det_g
                     tmp20 = tmp20*e
                     fg = ((-(4.0_ki*(12.0_ki*(tmp17+tmp11*f-(tmp16*f+tmp12))+8.0_ki*(tmp28+tmp24*tmp&
                     &53-(tmp27*tmp54+tmp25))+2.0_ki*(tmp57+tmp18*tmp30-(tmp56*tmp9+tmp13*tmp48))+4.0&
                     &_ki*(tmp32+tmp42+tmp46*f+tmp44*f*z-(tmp41*f+tmp31*f+tmp47+tmp45))+3.0_ki*(tmp60&
                     &+tmp61+tmp20*z+tmp39*tmp48-(tmp62*z+tmp61*z+tmp58*z+tmp30*tmp37))+6.0_ki*(tmp21&
                     &+tmp34+tmp33*tmp53+tmp55*z+tmp56*z+tmp48*b*e-(tmp57*z+tmp33*tmp54+tmp30*tmp43+t&
                     &mp55+tmp35+tmp22))+tmp58+tmp30*det_g+tmp37*tmp48+tmp62*tmp9-(tmp60*z+tmp48*tmp5&
                     &9+tmp30*tmp39+tmp20))/(tmp52*den1)*g1*num_b1+(18.0_ki*(tmp2*z-tmp2)+8.0_ki*(tmp&
                     &25+tmp27*tmp29+b*f*num_b1-(a*f*num_b1*z+tmp5*tmp5+tmp28))+2.0_ki*(tmp40+tmp48*t&
                     &mp50+tmp51*g1*z-(tmp38*tmp9+tmp30*num_b1+tmp30*tmp5))+3.0_ki*(tmp1*tmp30+tmp26*&
                     &num_b1+tmp49*z-(tmp23*num_b1+tmp1*tmp48+tmp49))+4.0_ki*(tmp45+tmp47+tmp14*tmp9+&
                     &tmp8*tmp9+tmp10*f*g1+tmp43*e*g1-(tmp10*tmp4+tmp7+tmp42+tmp19))+6.0_ki*(tmp35+tm&
                     &p15*tmp30+tmp29*tmp33+tmp32*z+tmp36*z+tmp38*z-(tmp37*e*g1+tmp40*z+tmp36+tmp34+t&
                     &mp32+tmp21+tmp1*tmp1))+12.0_ki*(tmp12+tmp6+tmp8+tmp19*z+tmp15*num_b1*z+tmp3*f*g&
                     &1-(tmp18*e*g1+tmp7*tmp9+tmp6*z+tmp5*num_b1+tmp3*tmp4+tmp14*z+tmp17))+9.0_ki*tmp&
                     &22-tmp21*tmp9)/tmp52))*g1+(2.0_ki*(2.0_ki*tmp15*z+tmp26*z)+6.0_ki*(tmp5*z-(tmp5&
                     &+tmp1))+3.0_ki*(tmp50+3.0_ki*tmp1*z+tmp23*z-(tmp9*num_b1+tmp23))+3.0_ki/2.0_ki*&
                     &tmp4+tmp50*tmp9-(tmp1*tmp9+num_b1))/den1)
                     !
                  case(2)
                     !
                     tmp1 = g2*g2
                     tmp2 = d*g2
                     tmp3 = num_b2*z
                     tmp4 = z*z
                     tmp5 = c*g2
                     tmp6 = e*g2
                     tmp7 = a*c
                     tmp8 = a*d
                     tmp9 = b*g2
                     tmp10 = a*g2
                     tmp11 = f*g2
                     tmp12 = a*a
                     tmp13 = d*d
                     tmp14 = c*c
                     tmp15 = c*d
                     tmp16 = b*c
                     tmp17 = b*d
                     tmp18 = e*e
                     tmp19 = d*det_g
                     tmp20 = c*e
                     tmp21 = c*det_g
                     tmp22 = tmp21*g2
                     tmp23 = den2*den2
                     tmp24 = tmp7*e
                     tmp25 = tmp16*d
                     tmp26 = e*f
                     tmp27 = d*f
                     tmp28 = b*b
                     tmp29 = tmp7*d
                     tmp30 = tmp16*e
                     tmp31 = b*det_g
                     tmp32 = det_g*z
                     fg = ((-(4.0_ki*(3.0_ki*(tmp13*tmp20+tmp19*e*z-tmp15*tmp18)+2.0_ki*(tmp17*tmp18+&
                     &tmp17*tmp4*det_g-tmp13*a*e)+8.0_ki*(tmp12*tmp26+tmp2*tmp28-(tmp27*tmp28+tmp12*t&
                     &mp6))+12.0_ki*(tmp24*f+tmp25*g2-(tmp25*f+tmp24*g2))+4.0_ki*(tmp29*g2+tmp30*f+tm&
                     &p31*g2*z-(tmp31*f*z+tmp30*g2+tmp29*f))+6.0_ki*(tmp14*tmp2+tmp14*tmp26+tmp22*z+t&
                     &mp13*b*e-(tmp21*f*z+tmp18*tmp8+tmp14*tmp6+tmp14*tmp27))+tmp13*tmp15+tmp18*tmp32&
                     &+tmp32*tmp32*z-(tmp21*tmp4*e+tmp18*tmp20))/(tmp23*den2)*g2*num_b2+(8.0_ki*(tmp1&
                     &1*tmp12+b*f*num_b2*z-tmp10*tmp10)+2.0_ki*(tmp20*tmp4*num_b2+det_g*e*g2*z-(tmp18&
                     &*tmp3+tmp10*tmp13))+4.0_ki*(tmp16*f*g2+tmp17*e*g2-(tmp4*det_g*num_b2*z+tmp17*tm&
                     &p4*num_b2+tmp1*tmp16))+12.0_ki*(tmp7*f*g2+c*f*num_b2*z-(tmp9*num_b2*z+tmp8*e*g2&
                     &+tmp1*tmp7))+3.0_ki*(tmp13*tmp5+tmp2*num_b2+tmp19*g2*z-(tmp6*num_b2+tmp18*tmp5)&
                     &)+6.0_ki*(tmp11*tmp14+tmp13*tmp9-(d*e*num_b2*z+3.0_ki*tmp5*num_b2*z+tmp15*e*g2+&
                     &tmp5*tmp5))-tmp22*tmp4)/tmp23))*g2+(3.0_ki/2.0_ki*(tmp1+2.0_ki*tmp2*z)+2.0_ki*t&
                     &mp6*z-(tmp4*tmp5+tmp3*tmp4))/den2)
                     !
                  case(3)
                     !
                     tmp1 = z*z
                     tmp2 = a*c
                     tmp3 = g3*g3
                     tmp4 = b*c
                     tmp5 = c*g3
                     tmp6 = f*g3
                     tmp7 = c*c
                     tmp8 = a*g3
                     tmp9 = b*g3
                     tmp10 = c*d
                     tmp11 = a*a
                     tmp12 = b*b
                     tmp13 = a*d
                     tmp14 = a*e
                     tmp15 = b*d
                     tmp16 = e*e
                     tmp17 = d*d
                     tmp18 = a*det_g
                     tmp19 = tmp18*g3
                     tmp20 = b*det_g
                     tmp21 = tmp20*g3
                     tmp22 = c*det_g
                     tmp23 = c*e
                     tmp24 = d*g3
                     tmp25 = e*g3
                     tmp26 = num_b3*z
                     tmp27 = den3*den3
                     tmp28 = tmp2*e
                     tmp29 = tmp4*d
                     tmp30 = e*f
                     tmp31 = d*f
                     tmp32 = tmp2*d
                     tmp33 = tmp4*e
                     tmp34 = det_g*z
                     fg = (((6.0_ki*(tmp16*tmp8+tmp17*tmp9)+3.0_ki*(tmp25*num_b3-tmp24*num_b3)+16.0_k&
                     &i*(tmp2*f*g3+tmp4*f*g3-(tmp3*tmp4+tmp2*tmp3))+12.0_ki*(tmp6*tmp7+tmp9*num_b3*z-&
                     &(tmp8*num_b3*z+tmp10*e*g3+tmp5*tmp5))+4.0_ki*(tmp1*tmp19+tmp1*tmp21+tmp1*tmp22*&
                     &g3+tmp1*tmp23*num_b3-tmp1*tmp10*num_b3)+2.0_ki*(tmp16*tmp26+d*det_g*g3*z+det_g*&
                     &e*g3*z-(tmp17*tmp8+tmp17*tmp26+tmp16*tmp9))+8.0_ki*(tmp11*tmp6+tmp12*tmp6+tmp1*&
                     &tmp14*num_b3+a*f*num_b3*z-(b*f*num_b3*z+tmp15*e*g3+tmp13*e*g3+tmp1*tmp15*num_b3&
                     &+tmp9*tmp9+tmp8*tmp8)))/tmp27-(4.0_ki*(3.0_ki*(tmp17*tmp23-tmp10*tmp16)+8.0_ki*&
                     &(tmp11*tmp30+tmp12*tmp24-(tmp12*tmp31+tmp11*tmp25))+2.0_ki*(tmp15*tmp16+tmp1*tm&
                     &p22*e-(tmp1*tmp10*det_g+tmp14*tmp17))+12.0_ki*(tmp28*f+tmp29*g3-(tmp29*f+tmp28*&
                     &g3))+6.0_ki*(tmp24*tmp7+tmp30*tmp7+tmp17*b*e-(tmp31*tmp7+tmp25*tmp7+tmp13*tmp16&
                     &))+4.0_ki*(tmp21*z+tmp32*g3+tmp33*f+tmp1*tmp18*e+tmp18*f*z-(tmp20*f*z+tmp1*tmp1&
                     &5*det_g+tmp33*g3+tmp32*f+tmp19*z))+tmp10*tmp17+tmp16*tmp34-(tmp17*tmp34+tmp16*t&
                     &mp23))/(tmp27*den3)*g3*num_b3+(2.0_ki*(d*z+e*z)+4.0_ki*(tmp1*a+tmp1*b+tmp1*c)+3&
                     &.0_ki*g3)/den3))*g3)
                     !
                  case(4)
                     !
                     tmp1 = num_b1*z
                     tmp2 = z*z
                     tmp3 = c*g1
                     tmp4 = tmp3*z
                     tmp5 = b*g1
                     tmp6 = e*g1
                     tmp7 = a*g1*z
                     tmp8 = c*f
                     tmp9 = det_g*z
                     tmp10 = tmp2*det_g
                     tmp11 = a*e
                     tmp12 = tmp11*z
                     tmp13 = b*d
                     tmp14 = tmp13*z
                     tmp15 = d*e
                     tmp16 = c*d
                     tmp17 = tmp16*z
                     tmp18 = c*e
                     tmp19 = tmp18*z
                     tmp20 = b*f
                     tmp21 = a*f
                     tmp22 = e*e
                     tmp23 = d*d
                     tmp24 = c*g2
                     fg = (((2.0_ki*(tmp19-tmp17)+4.0_ki*(tmp12+tmp21+b*g3-(a*g3+tmp20+tmp14))+tmp22-&
                     &tmp23)/(den3*den3)*num_b3*z-(2.0_ki*(a*z+b*z+c*z)+d+e)/den3*z)*g3+(3.0_ki/2.0_k&
                     &i*(3.0_ki*(tmp1+tmp4-tmp2*num_b1)+2.0_ki*(tmp5*z-(tmp5+tmp3))+tmp1*tmp2+tmp6*z-&
                     &(num_b1+tmp6))+2.0_ki*tmp7+d*g1*z-1.0_ki/2.0_ki*tmp2*tmp3)/den1+(3.0_ki*(2.0_ki&
                     &*(tmp24-tmp8)+tmp15)+2.0_ki*(2.0_ki*(b*g2-tmp20)+tmp10+tmp14)+tmp22-tmp19)/(den&
                     &2*den2)*g2*num_b2*z-((2.0_ki*(2.0_ki*(tmp20+tmp7-(tmp21*z+tmp5))+tmp11*tmp2+tmp&
                     &2*tmp9-(det_g+tmp13))+3.0_ki*(2.0_ki*(tmp14+tmp4+tmp8+tmp9-(tmp8*z+tmp3+tmp12+t&
                     &mp10))+tmp17+tmp15*z-(tmp19+tmp15))+tmp18+tmp23*z-(tmp16*tmp2+tmp22))/(den1*den&
                     &1)*g1*num_b1+1.0_ki/2.0_ki*(3.0_ki*(tmp2*num_b2-d*g2)+tmp24*z-2.0_ki*e*g2)/den2&
                     &*z))
                     !
                  end select
                  !
               end select
               !
            end select
            !
         end select
         !
        else
         !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
          & 'Unexpected value for nb_par = %d0'
          tab_erreur_par(2)%arg_int = nb_par
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      else if (dim == "n+2") then
        !
        if (nb_par == 0) then
          !
          select case (flag)
            !
            case(1)
            !
            fg = num_b1*g1/(2._ki*det_g*g1+det_s)
            !
            case(2)
            !
            fg = num_b2*g2/(2._ki*det_g*g2+det_s)
            !
            case(3)
            !
            fg = num_b3*g3/(2._ki*det_g*g3+det_s)
            !
            case(4)
            !
            fg = -1._ki/2._ki
            !
          end select
        else if (nb_par == 1) then
          !
          select case(par3)
          !
          case(1)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-b*g1**2/den1+num_b1*g1**2*(-c**2+4._ki*a*b-e*c+2._ki*b*d)/den&
                &1**2
              !
            case(2)
              !
              fg=-1._ki/2._ki*g2*(c*g2-2._ki*num_b2+2._ki*num_b2*z)/den2+num_b2&
                &*g2**2*(-det_g+4._ki*a*b-c**2+det_g*z-e*c+2._ki*b*d)/den2**2
              !
            case(3)
              !
              fg=1._ki/2._ki*g3*(g3*c+2._ki*g3*b+2._ki*num_b3-2._ki*z*num_b3)/d&
                &en3+num_b3*g3**2*(-det_g+4._ki*a*b-c**2+det_g*z-e*c+2._ki*b*d)/&
                &den3**2
              !
            case(4)
              !
              fg=-1._ki/4._ki-1._ki/2._ki*num_b3*g3*(z-1._ki)/den3-1._ki/2._ki*&
                &num_b2*g2*(z-1._ki)/den2
              !
            end select
            !
          case(2)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-1._ki/2._ki*g1*(c*g1-2._ki*num_b1*z)/den1-num_b1*g1**2*(det_g&
                &*z-c*d+2._ki*e*a)/den1**2
              !
            case(2)
              !
              fg=-a*g2**2/den2-num_b2*g2**2*(2._ki*e*a-c*d)/den2**2
              !
            case(3)
              !
              fg=1._ki/2._ki*g3*(2._ki*g3*a+g3*c+2._ki*z*num_b3)/den3-num_b3*g3&
                &**2*(det_g*z-c*d+2._ki*e*a)/den3**2
              !
            case(4)
              !
              fg=1._ki/2._ki*num_b3*z*g3/den3-1._ki/4._ki+1._ki/2._ki*num_b1*z*&
                &g1/den1
              !
            end select
            !
          case(3)
            !
            select case(flag)
            !
            case(1)
              !
              fg=1._ki/2._ki*g1*(c*g1+2._ki*g1*b+2._ki*num_b1-2._ki*num_b1*z)/d&
                &en1+num_b1*g1**2*(-det_g+e*c-2._ki*b*d+det_g*z-c*d+2._ki*e*a)/d&
                &en1**2
              !
            case(2)
              !
              fg=1._ki/2._ki*g2*(2._ki*g2*a+c*g2+2._ki*num_b2*z)/den2+num_b2*g2&
                &**2*(-det_g*z+e*c-2._ki*b*d-c*d+2._ki*e*a)/den2**2
              !
            case(3)
              !
              fg=-g3**2*(c+a+b)/den3+num_b3*g3**2*(2._ki*e*a+e*c-c*d-2._ki*b*d)&
                &/den3**2
              !
            case(4)
              !
              fg=1._ki/2._ki*num_b2*z*g2/den2-1._ki/4._ki-1._ki/2._ki*num_b1*g1&
                &*(z-1._ki)/den1
              !
            end select
            !
          end select
          !
        else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
          & 'Unexpected value for nb_par = %d0'
          tab_erreur_par(2)%arg_int = nb_par
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      else 
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = &
        & 'Unexpected value for dim = %c0'
        tab_erreur_par(2)%arg_char = dim
        call catch_exception(0)
        !
        stop
        !
      end if
      !
    end function fg
    !
end module function_3p_finite
