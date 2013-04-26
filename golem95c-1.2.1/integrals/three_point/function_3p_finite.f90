! 
!~ 30.4.2010: adapted from function_3p3m.f90 for internal masses
! 
!~ 24.6.2010: uses Andre van Hameren's OneLOop for finite C0
!~ 14.1.2011: include LT option in addition
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
!  * f3p_finite_c -- a function which computes the same thing as f3p_finite,  
!    only the format of the return values is different
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
  private :: ki
  real(ki) :: s12_glob,s23_glob,s13_glob
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
  private :: eps_glob,s12_glob,s23_glob,s13_glob,par1_glob,par2_glob,par3_glob,dim_glob
  private :: b_real, b_complex, sumb_real, sumb_complex, invs_real, invs_complex, s_mat_real, s_mat_complex, par, s
  private :: deja_calcule,resultat,deja_calcule2,resultat2,deja_calcule_np2,resultat_np2,deja_calcule22,resultat22
  private :: deja_calcule_c,resultat_c,deja_calcule2_c,resultat2_c,deja_calcule_np2_c,resultat_np2_c,deja_calcule22_c,resultat22_c
  private :: s_mat_p_loc
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
  interface
     subroutine avh_olo_c0m(rslt,p1,p2,p3,m1,m2,m3)
       use precision_golem, only: ki_avh
       implicit none
       complex(ki_avh), intent(out) :: rslt(0:2)
       real(ki_avh), intent(in) :: p1,p2,p3,m1,m2,m3
     end subroutine avh_olo_c0m
  end interface
  !
  !
  interface
     subroutine avh_olo_c0c(rslt,p1,p2,p3,m1,m2,m3)
       use precision_golem, only: ki_avh
       implicit none
       complex(ki_avh), intent(out) :: rslt(0:2)
       complex(ki_avh), intent(in) :: p1,p2,p3,m1,m2,m3
     end subroutine avh_olo_c0c
  end interface
  !
  ! added Jan2011 for LT option
!AC!        interface
!AC!         function C0(p1, p2, p1p2, m1, m2, m3)
!AC!            use precision_golem, only: ki_lt
!AC!            implicit none
!AC!            real(ki_lt), intent(in) :: p1, p2, p1p2, m1, m2, m3
!AC!           complex(ki_lt) :: C0
!AC!         end function
!AC!      end interface
!AC!      interface
!AC!         function C0C(p1, p2, p1p2, m1, m2, m3)
!AC!            use precision_golem, only: ki_lt
!AC!            implicit none
!AC!            real(ki_lt), intent(in) :: p1, p2, p1p2
!AC!            complex(ki_lt), intent(in) :: m1, m2, m3
!AC!            complex(ki_lt) :: C0C
!AC!         end function
!AC!      end interface
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
    integer :: nb_par
    real(ki) :: lamb, detS3
    real(ki) :: plus_grand
    real(ki) :: norma
    real(ki) :: s1p,s2p,s3p,m1p,m2p,m3p
    !complex(ki) :: resto,abserro
    !
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
    plus_grand = 1._ki
      if (rat_or_tot_par%tot_selected) then
    !
      plus_grand = maxval(array=abs(s_mat_real))
!      plus_grand = 1._ki
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
    !
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
!    lamb = 2._ki*s_mat_real(1,3)*s_mat_real(2,3)+2._ki*s_mat_real(1,2)*s_mat_real(2,3)&
!         +2._ki*s_mat_real(1,3)*s_mat_real(1,2)-s_mat_real(1,3)*s_mat_real(1,3)-s_mat_real(1,2)*s_mat_real(1,2)&
!         -s_mat_real(2,3)*s_mat_real(2,3)
    !
      lamb = (-s_mat_real(1, 2)**2 - s_mat_real(1, 3)**2 + s_mat_real(1, 1)*s_mat_real(2, 2) - &
         &   2*s_mat_real(1, 3)*(s_mat_real(2, 2) - s_mat_real(2, 3)) - 2*s_mat_real(1, 1)*s_mat_real(2, 3) - &
         &   s_mat_real(2, 3)**2 + 2*s_mat_real(1, 2)*(s_mat_real(1, 3) + s_mat_real(2, 3) - s_mat_real(3, 3)) + &
         &   s_mat_real(1, 1)*s_mat_real(3, 3) + s_mat_real(2, 2)*s_mat_real(3, 3))
    !
    nb_par = count(mask=par/=0)
    !
    m1p = -s_mat_real(1,1)/2._ki
    m2p = -s_mat_real(2,2)/2._ki
    m3p = -s_mat_real(3,3)/2._ki
    s1p = s_mat_real(1,3) - (s_mat_real(1,1)+s_mat_real(3,3))/2._ki
    s2p = s_mat_real(1,2) - (s_mat_real(1,1)+s_mat_real(2,2))/2._ki
    s3p = s_mat_real(2,3) - (s_mat_real(2,2)+s_mat_real(3,3))/2._ki
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
       norma =  -1._ki/24._ki
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
    !      if (abs(sumb_real) > coupure_3p3m) then
    ! if (abs(sumb_real) > 0._ki) then
       !
       ! always use analytic computation until massive numerical is implemented
        ! branching removed completely Feb 4, 2011 because if detS3=0, 
       ! sumb_real=NAN so numerical branch was entered, although scalar integral 
       ! still finite if detS3=0 without any z_i=0 (LLS)
       !
      ! when invariants as input are used, no division by plus_grand!
       !
       if (dim == "ndi") then
          !
 !         f3p_finite_rarg(3:4)= a3pC0i_rarg(s1,s2,s3,m1,m2,m3,par1,par2,par3)
	  f3p_finite_rarg(3:4)= a3pC0i_rarg(s1p,s2p,s3p,m1p,m2p,m3p,par1,par2,par3)&
           &/plus_grand
          !
       else if (dim == "n+2") then
          !
!          f3p_finite_rarg = a3pC0i_np2_rarg(s1,s2,s3,m1,m2,m3,par1,par2,par3)
         f3p_finite_rarg = a3pC0i_np2_rarg(s1p,s2p,s3p,m1p,m2p,m3p,par1,par2,par3)
         f3p_finite_rarg(3) = f3p_finite_rarg(3)-log(plus_grand)*norma
!          
	  ! mu2_scale_par is already contained in the bubbles, 
          !
       end if
       !
    ! else
       !
       ! numerical computation
       !
       !dim_glob = dim
       !par1_glob = par1
      ! par2_glob = par2
       !par3_glob = par3
       !
       !s13_glob = s_mat_real(1,3)
       !s12_glob = s_mat_real(1,2)
       !s23_glob = s_mat_real(2,3)
       !
       !resto = 0._ki
       !abserro = 0._ki
       !
       ! eps_glob = sign(1._ki,s23_glob-s13_glob)
       !
       !origine_info_par = "f3p_finite, dimension "//dim
      ! num_grand_b_info_par = lamb
      ! denom_grand_b_info_par = (2._ki*s_mat_real(1,3)*s_mat_real(1,2)*s_mat_real(2,3))
       !
       ! call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,resto,abserro)
      !
      ! if (dim == "ndi") then      
          !
      !     resto = resto/plus_grand
          !
      !  else if (dim == "n+2") then
          !
      !    f3p_finite_rarg(1) = norma
      !    f3p_finite_rarg(2) = 0._ki
      !     resto = resto-log(plus_grand/mu2_scale_par)*norma
          !
      !  end if
       !
      !  f3p_finite_rarg(3) = real(resto,ki)
      !  f3p_finite_rarg(4) = aimag(resto)
       !
    !end if  ! end if analytic or numeric
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
    !complex(ki) :: resto,abserro
    complex(ki) :: temp0
    complex(ki), dimension(2) :: temp
    real(ki) :: s1rp,s2rp,s3rp
    complex(ki) :: m1p,m2p,m3p
    !
!    write(6,*) '******** debug   f3p_finite_carg *************** '
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
    !plus_grand = 1._ki
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
    m1p = m1/plus_grand
    m2p = m2/plus_grand
    m3p = m3/plus_grand
    s1rp = s1r/plus_grand
    s2rp = s2r/plus_grand
    s3rp = s3r/plus_grand
    !
    s_mat_p_loc = assign_s_matrix(s_mat_complex,s_mat_real)
    !
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
    !
    nb_par = count(mask=par/=0)
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
       norma = -1._ki/24._ki
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
    !      if (abs(sumb_complex) > coupure_3p3m) then
    ! if (abs(sumb_complex) > 0._ki) then
       !
       ! always use analytic computation until massive numerical is implemented
       ! 
       !
       if (dim == "ndi") then
          !
!          temp0 = a3pC0i_carg(s1r,s2r,s3r,m1,m2,m3,par1,par2,par3)
          temp0 = a3pC0i_carg(s1rp,s2rp,s3rp,m1p,m2p,m3p,par1,par2,par3)
          f3p_finite_carg(3) = real(temp0,ki)/plus_grand
          f3p_finite_carg(4) = aimag(temp0)/plus_grand
          !
       else if (dim == "n+2") then
          !
!          temp = a3pC0i_np2_carg(s1r,s2r,s3r,m1,m2,m3,par1,par2,par3)
          temp = a3pC0i_np2_carg(s1rp,s2rp,s3rp,m1p,m2p,m3p,par1,par2,par3)
          f3p_finite_carg(1) = real(temp(1),ki)
          f3p_finite_carg(2) = aimag(temp(1))
          f3p_finite_carg(3) = real(temp(2),ki)
          f3p_finite_carg(4) = aimag(temp(2))
          f3p_finite_carg(3) = f3p_finite_carg(3)-log(plus_grand)*norma
	  ! mu2_scale_par is already contained in the bubbles, 
          !
       end if
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
 !
  !
  !*****
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
 !    write(6,*) 'calls a3pC0i_rarg'
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
    real(ki) :: del
    complex(ki) :: C0_rarg
    complex(ki_avh), dimension(0:2) :: C0olo
 !AC!   complex(ki_lt) :: C0
      del = epsilon(1._ki)
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
       ! changed to include LT option Jan2011
          ! use avh_olo or LT in finite case
!AC!    if (withlt) then
       !
!AC!    C0_rarg = C0(real(s1,ki_lt),real(s2,ki_lt),real(s3,ki_lt),&
!AC!        & real(m3sq,ki_lt),real(m1sq,ki_lt),real(m2sq,ki_lt))
       !
!AC!    else 
    ! use avh_olo
    ! on-shell cutoff decreased 26.4.2013 GH
         if (.not. olo) then
              call avh_olo_onshell(100._ki*del)
          !call avh_olo_onshell(1.e-10_ki)
          call avh_olo_mu_set(sqrt(mu2_scale_par))
          olo=.true.
         end if
       !
       !
        call avh_olo_c0m(C0olo,real(s1r,ki_avh),real(s2r,ki_avh),real(s3r,ki_avh), &
            &    real(m3sq,ki_avh),real(m1sq,ki_avh),real(m2sq,ki_avh))
       !
       !
       C0_rarg=C0olo(0)
       !
!AC!       end if ! end if withlt
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
    real(ki) :: del
    complex(ki) :: C0_carg
    complex(ki_avh), dimension(0:2) :: C0olo
   !AC!   complex(ki_lt) :: C0C
    !integer :: i,j,k
      del = epsilon(1._ki)
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
       ! changed to include LT option Jan2011
       ! use avh_olo or LT in finite case
!AC!     if (withlt) then
       !
!AC!     C0_carg = C0C(real(s1,ki_lt),real(s2,ki_lt),real(s3,ki_lt),&
!AC!         & cmplx(m3sq,kind=ki_lt),cmplx(m1sq,kind=ki_lt),cmplx(m2sq,kind=ki_lt))
       !
!AC!    else 
       ! use avh_olo
    ! on-shell cutoff decreased 26.4.2013 GH
       if (.not. olo) then
          call avh_olo_onshell(100._ki*del)
   !       call avh_olo_onshell(1.e-10_ki)
          call avh_olo_mu_set(sqrt(mu2_scale_par))
          olo=.true.
        end if
       !
       cp1 = cmplx(s1r,0._ki_avh,kind=ki_avh)
       cp2 = cmplx(s2r,0._ki_avh,kind=ki_avh)
       cp3 = cmplx(s3r,0._ki_avh,kind=ki_avh)
       !
       cm1 = cmplx(m1sq,kind=ki_avh)
       cm2 = cmplx(m2sq,kind=ki_avh)
       cm3 = cmplx(m3sq,kind=ki_avh)
       !
       call avh_olo_c0c(C0olo,cp1,cp2,cp3,cm3,cm1,cm2)
       !
       C0_carg=C0olo(0)
       !
!AC!       end if !end if withlt
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
  !
end module function_3p_finite
