! 
!****h* src/interface/tensor_integrals
! NAME
!
!  Module tensor_integrals
!
! USAGE
!
!  use tensor_integrals
!
! DESCRIPTION
!
!  This module provides an interface which allows to compute
!  tensor integrals rather than form factors.
!
! OUTPUT
!
!  This module exports the functions:
!  * init_smat   -- initialize the s_smat from vectors and masses
!  * ti1         -- tensor tadpoles
!  * ti2         -- tensor bubbles
!  * ti3         -- tensor triangles
!  * ti4         -- tensor boxes
!  * ti5         -- tensor pentagons
!  * ti6         -- tensor hexagons
!
! USES
!
!  precision_golem
!  form_factor_type
!  form_factor_1p
!  form_factor_2p
!  form_factor_3p
!  form_factor_4p
!  form_factor_5p
!  form_factor_6p
!  cache
!  matrice_s
!  spinor
!  array
!
!*****
module     tensor_integrals
use precision_golem, only: ki
use form_factor_type, only: form_factor, operator(*), operator(+)
use form_factor_1p, only: a10
use form_factor_2p, only: a20, a21, a22, b22
use form_factor_3p, only: a30, a31, a32, a33, b32, b33
use form_factor_4p, only: a40, a41, a42, a43, a44, b42, b43, b44, c44
use form_factor_5p, only: a50, a51, a52, a53, a54, a55, b52, b53, b54, b55, &
                        & c54, c55
use form_factor_6p, only: a60, a61, a62, a63, a64, a65, a66
use form_factor_higher_ranks, only: a12, b12, a23, b23, a34, b34, a45, b45, c45, &
                & a56, b56, c56, d56, a67
use cache, only: allocate_cache, clear_cache
use matrice_s, only: set_ref, s_mat, allocation_s, deallocation_s, init_invs, &
    &                b_ref
use spinor, only: scalar
use array, only: packb, unpackb, pminus
implicit none
private

private :: a10, a20, a21, a22, b22
private :: a30, a31, a32, a33, b32, b33
private :: a40, a41, a42, a43, a44, b42, b43, b44, c44 
private :: a50, a51, a52, a53, a54, a55, b52, b53, b54, b55, c54, c55
private :: a60, a61, a62, a63, a64, a65, a66

private :: ki, form_factor, allocate_cache, clear_cache, scalar
private :: packb, unpackb, pminus
private :: set_ref, s_mat, allocation_s, deallocation_s, init_invs, b_ref
integer, dimension(0), target, private :: loc_s_null = 0

private :: symmetric_A_coeff1
private :: symmetric_A_coeff2
private :: symmetric_A_coeff3
private :: symmetric_A_coeff4
private :: symmetric_A_coeff5
private :: symmetric_A_coeff6
private :: symmetric_B_coeff2
private :: symmetric_B_coeff3
private :: symmetric_B_coeff4
private :: symmetric_B_coeff5
private :: symmetric_C_coeff4
private :: symmetric_C_coeff5
private :: symmetric_D_coeff6

interface     symmetric_A_coeff
   module procedure symmetric_A_coeff1
   module procedure symmetric_A_coeff2
   module procedure symmetric_A_coeff3
   module procedure symmetric_A_coeff4
   module procedure symmetric_A_coeff5
   module procedure symmetric_A_coeff6
   module procedure symmetric_A_coeff7
end interface symmetric_A_coeff

interface     symmetric_B_coeff
   module procedure symmetric_B_coeff2
   module procedure symmetric_B_coeff3
   module procedure symmetric_B_coeff4
   module procedure symmetric_B_coeff5
   module procedure symmetric_B_coeff6
end interface symmetric_B_coeff

interface     symmetric_C_coeff
   module procedure symmetric_C_coeff4
   module procedure symmetric_C_coeff5
   module procedure symmetric_C_coeff6
end interface symmetric_C_coeff

interface     symmetric_D_coeff
   module procedure symmetric_D_coeff6
end interface symmetric_D_coeff

interface     init_smat
   module procedure init_smat1
   module procedure init_smat2
   module procedure init_smat3
   module procedure init_smat4
   module procedure init_smat5
   module procedure init_smat6
end interface init_smat

interface     ti1
   module procedure ti1r0
   module procedure ti1r1
   module procedure ti1r2
end interface ti1

interface     ti2
   module procedure ti2r0
   module procedure ti2r1
   module procedure ti2r2
   module procedure ti2r3
end interface ti2

interface     ti3
   module procedure ti3r0
   module procedure ti3r1
   module procedure ti3r2
   module procedure ti3r3
   module procedure ti3r4
end interface ti3

interface     ti4
   module procedure ti4r0
   module procedure ti4r1
   module procedure ti4r2
   module procedure ti4r3
   module procedure ti4r4
end interface ti4

interface     ti5
   module procedure ti5r0
   module procedure ti5r1
   module procedure ti5r2
   module procedure ti5r3
   module procedure ti5r4
   module procedure ti5r5
   module procedure ti5r6
end interface ti5

interface     ti6
   module procedure ti6r0
   module procedure ti6r1
   module procedure ti6r2
   module procedure ti6r3
   module procedure ti6r4
   module procedure ti6r5
   module procedure ti6r6
   module procedure ti6r7
end interface ti6

private :: ti1r0, ti1r1, ti1r2
private :: ti2r0, ti2r1, ti2r2, ti2r3
private :: ti3r0, ti3r1, ti3r2, ti3r3, ti3r4
private :: ti4r0, ti4r1, ti4r2, ti4r3, ti4r4, ti4r5
private :: ti5r0, ti5r1, ti5r2, ti5r3, ti5r4, ti5r5, ti5r6
private :: ti6r0, ti6r1, ti6r2, ti6r3, ti6r4, ti6r5, ti6r6

integer, parameter, public :: use_existing_smat = 1
integer, parameter, public :: keep_smat_on_exit = 2

public :: ti1, ti2, ti3, ti4, ti5, ti6, init_smat, done_smat

contains

pure elemental function chop(val,prec)
   implicit none
   real(ki), intent(in) :: val
   real(ki), optional, intent(in) :: prec
   real(ki) :: chop

   if (present(prec)) then
      if (abs(val) .gt. prec) then
         chop = val
      else
         chop = 0.0_ki
      end if
   else
      if (abs(val) .gt. 1.0E+04_ki*epsilon(1.0_ki)) then
         chop = val
      else
         chop = 0.0_ki
      end if
   end if
end  function chop

!---#[ init_smat :
subroutine     init_smat1(m1sq)
   implicit none
   real(ki), intent(in) :: m1sq

   call allocation_s(1)
   set_ref = (/1/)
   b_ref = packb(set_ref)
   s_mat(1,1) = -m1sq-m1sq

   call allocate_cache(1)
   call init_invs()

end subroutine init_smat1

subroutine     init_smat2(r1,m1sq,m2sq)
   implicit none
   real(ki), dimension(4), intent(in) :: r1
   real(ki), intent(in) :: m1sq,m2sq

   call allocation_s(2)
   set_ref = (/1,2/)
   b_ref = packb(set_ref)
   s_mat(1,1) = -m1sq-m1sq
   s_mat(1,2) = chop(scalar(r1,r1)) - m1sq - m2sq
   s_mat(2,1) = s_mat(1,2)
   s_mat(2,2) = -m2sq-m2sq

   call allocate_cache(2)
   call init_invs()
end subroutine init_smat2

subroutine     init_smat3(r1,r2,m1sq,m2sq,m3sq)
   implicit none
   real(ki), dimension(4), intent(in) :: r1,r2
   real(ki), intent(in) :: m1sq,m2sq,m3sq

   real(ki), dimension(4) :: delt

   call allocation_s(3)
   set_ref = (/1,2,3/)
   b_ref = packb(set_ref)
   s_mat(1,1) = -m1sq-m1sq
   delt = r1-r2
   s_mat(1,2) = chop(scalar(delt,delt)) - m1sq - m2sq
   delt = r1
   s_mat(1,3) = chop(scalar(delt,delt)) - m1sq - m3sq
   s_mat(2,1) = s_mat(1,2)
   s_mat(2,2) = -m2sq-m2sq
   delt = r2
   s_mat(2,3) = chop(scalar(delt,delt)) - m2sq - m3sq
   s_mat(3,1) = s_mat(1,3)
   s_mat(3,2) = s_mat(2,3)
   s_mat(3,3) = -m3sq-m3sq

   call allocate_cache(3)
   call init_invs()
end subroutine init_smat3

subroutine     init_smat4(r1,r2,r3,m1sq,m2sq,m3sq,m4sq)
   implicit none
   real(ki), dimension(4), intent(in) :: r1,r2,r3
   real(ki), intent(in) :: m1sq,m2sq,m3sq,m4sq

   real(ki), dimension(4) :: delt

   call allocation_s(4)
   set_ref = (/1,2,3,4/)
   b_ref = packb(set_ref)
   s_mat(1,1) = -m1sq-m1sq
   delt = r1-r2
   s_mat(1,2) = chop(scalar(delt,delt)) - m1sq - m2sq
   delt = r1-r3
   s_mat(1,3) = chop(scalar(delt,delt)) - m1sq - m3sq
   delt = r1
   s_mat(1,4) = chop(scalar(delt,delt)) - m1sq - m4sq
   s_mat(2,1) = s_mat(1,2)
   s_mat(2,2) = -m2sq-m2sq
   delt = r2-r3
   s_mat(2,3) = chop(scalar(delt,delt)) - m2sq - m3sq
   delt = r2
   s_mat(2,4) = chop(scalar(delt,delt)) - m2sq - m4sq
   s_mat(3,1) = s_mat(1,3)
   s_mat(3,2) = s_mat(2,3)
   s_mat(3,3) = -m3sq-m3sq
   delt = r3
   s_mat(3,4) = chop(scalar(delt,delt)) - m3sq - m4sq
   s_mat(4,1) = s_mat(1,4)
   s_mat(4,2) = s_mat(2,4)
   s_mat(4,3) = s_mat(3,4)
   s_mat(4,4) = -m4sq-m4sq

   call allocate_cache(4)
   call init_invs()
end subroutine init_smat4

subroutine     init_smat5(r1,r2,r3,r4,m1sq,m2sq,m3sq,m4sq,m5sq)
   implicit none
   real(ki), dimension(4), intent(in) :: r1,r2,r3,r4
   real(ki), intent(in) :: m1sq,m2sq,m3sq,m4sq,m5sq

   real(ki), dimension(4) :: delt

   call allocation_s(5)
   set_ref = (/1,2,3,4,5/)
   b_ref = packb(set_ref)
   s_mat(1,1) = -m1sq-m1sq
   delt = r1-r2
   s_mat(1,2) = chop(scalar(delt,delt)) - m1sq - m2sq
   delt = r1-r3
   s_mat(1,3) = chop(scalar(delt,delt)) - m1sq - m3sq
   delt = r1-r4
   s_mat(1,4) = chop(scalar(delt,delt)) - m1sq - m4sq
   delt = r1
   s_mat(1,5) = chop(scalar(delt,delt)) - m1sq - m5sq
   s_mat(2,1) = s_mat(1,2)
   s_mat(2,2) = -m2sq-m2sq
   delt = r2-r3
   s_mat(2,3) = chop(scalar(delt,delt)) - m2sq - m3sq
   delt = r2-r4
   s_mat(2,4) = chop(scalar(delt,delt)) - m2sq - m4sq
   delt = r2
   s_mat(2,5) = chop(scalar(delt,delt)) - m2sq - m5sq
   s_mat(3,1) = s_mat(1,3)
   s_mat(3,2) = s_mat(2,3)
   s_mat(3,3) = -m3sq-m3sq
   delt = r3-r4
   s_mat(3,4) = chop(scalar(delt,delt)) - m3sq - m4sq
   delt = r3
   s_mat(3,5) = chop(scalar(delt,delt)) - m3sq - m5sq
   s_mat(4,1) = s_mat(1,4)
   s_mat(4,2) = s_mat(2,4)
   s_mat(4,3) = s_mat(3,4)
   s_mat(4,4) = -m4sq-m4sq
   delt = r4
   s_mat(4,5) = chop(scalar(delt,delt)) - m4sq - m5sq
   s_mat(5,1) = s_mat(1,5)
   s_mat(5,2) = s_mat(2,5)
   s_mat(5,3) = s_mat(3,5)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = -m5sq-m5sq

   call allocate_cache(5)
   call init_invs()
end subroutine init_smat5

subroutine     init_smat6(r1,r2,r3,r4,r5,m1sq,m2sq,m3sq,m4sq,m5sq,m6sq)
   implicit none
   real(ki), dimension(4), intent(in) :: r1,r2,r3,r4,r5
   real(ki), intent(in) :: m1sq,m2sq,m3sq,m4sq,m5sq,m6sq

   real(ki), dimension(4) :: delt

   call allocation_s(6)
   set_ref = (/1,2,3,4,5,6/)
   b_ref = packb(set_ref)
   s_mat(1,1) = -m1sq-m1sq
   delt = r1-r2
   s_mat(1,2) = chop(scalar(delt,delt)) - m1sq - m2sq
   delt = r1-r3
   s_mat(1,3) = chop(scalar(delt,delt)) - m1sq - m3sq
   delt = r1-r4
   s_mat(1,4) = chop(scalar(delt,delt)) - m1sq - m4sq
   delt = r1-r5
   s_mat(1,5) = chop(scalar(delt,delt)) - m1sq - m5sq
   delt = r1
   s_mat(1,6) = chop(scalar(delt,delt)) - m1sq - m6sq
   s_mat(2,1) = s_mat(1,2)
   s_mat(2,2) = -m2sq-m2sq
   delt = r2-r3
   s_mat(2,3) = chop(scalar(delt,delt)) - m2sq - m3sq
   delt = r2-r4
   s_mat(2,4) = chop(scalar(delt,delt)) - m2sq - m4sq
   delt = r2-r5
   s_mat(2,5) = chop(scalar(delt,delt)) - m2sq - m5sq
   delt = r2
   s_mat(2,6) = chop(scalar(delt,delt)) - m2sq - m6sq
   s_mat(3,1) = s_mat(1,3)
   s_mat(3,2) = s_mat(2,3)
   s_mat(3,3) = -m3sq-m3sq
   delt = r3-r4
   s_mat(3,4) = chop(scalar(delt,delt)) - m3sq - m4sq
   delt = r3-r5
   s_mat(3,5) = chop(scalar(delt,delt)) - m3sq - m5sq
   delt = r3
   s_mat(3,6) = chop(scalar(delt,delt)) - m3sq - m6sq
   s_mat(4,1) = s_mat(1,4)
   s_mat(4,2) = s_mat(2,4)
   s_mat(4,3) = s_mat(3,4)
   s_mat(4,4) = -m4sq-m4sq
   delt = r4-r5
   s_mat(4,5) = chop(scalar(delt,delt)) - m4sq - m5sq
   delt = r4
   s_mat(4,6) = chop(scalar(delt,delt)) - m4sq - m6sq
   s_mat(5,1) = s_mat(1,5)
   s_mat(5,2) = s_mat(2,5)
   s_mat(5,3) = s_mat(3,5)
   s_mat(5,4) = s_mat(4,5)
   s_mat(5,5) = -m5sq-m5sq
   delt = r5
   s_mat(5,6) = chop(scalar(delt,delt)) - m5sq - m6sq
   s_mat(6,1) = s_mat(1,6)
   s_mat(6,2) = s_mat(2,6)
   s_mat(6,3) = s_mat(3,6)
   s_mat(6,4) = s_mat(4,6)
   s_mat(6,5) = s_mat(5,6)
   s_mat(6,6) = -m6sq-m6sq

   call allocate_cache(6)
   call init_invs()
end subroutine init_smat6

subroutine     done_smat()
   implicit none
   call clear_cache()
   call deallocation_s()
end subroutine done_smat
!---#] init_smat :
!---#[ One point tensor integrals :
subroutine ti1r0(tens,m1,flag,pinches)
   implicit none
   type(form_factor), intent(out) :: tens
   real(ki), intent(in) :: m1
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   logical :: f_init_smat, f_dispose_smat

   integer, dimension(:), pointer :: lpinches

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(m1*m1)

   tens = a10(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti1r0

subroutine ti1r1(tens,m1,flag,pinches)
   implicit none
   type(form_factor), dimension(4), intent(out) :: tens
   real(ki), intent(in) :: m1
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   logical :: f_init_smat, f_dispose_smat

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if(f_init_smat) call init_smat(m1*m1)

   tens(:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)

   if(f_dispose_smat) call done_smat()
end subroutine ti1r1
subroutine ti1r2(tens,m1,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4), intent(out) :: tens
   real(ki), intent(in) :: m1
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   real(ki), dimension(4,4) :: term
   logical :: f_init_smat, f_dispose_smat

   integer, dimension(:), pointer :: lpinches

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(m1*m1)

   call symmetric_B_coeff(term)
   tens(:,:) = term(:,:) * B12(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti1r2

!---#] One point tensor integrals :
!---#[ Two point tensor integrals :
subroutine ti2r0(tens,r1,m1,m2,flag,pinches)
   implicit none
   type(form_factor), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1
   real(ki), intent(in) :: m1, m2
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(r1,m1*m1,m2*m2)

   tens = a20(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti2r0

subroutine ti2r1(tens,r1,m1,m2,flag,pinches)
   implicit none
   type(form_factor), dimension(4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1
   real(ki), intent(in) :: m1, m2
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4) :: term
   integer, dimension(2) :: unpinched
   unpinched = unpackb(pminus(b_ref,packb(pinches)),2)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(r1,m1*m1,m2*m2)

   call symmetric_A_coeff(term,r1)
   tens(:) = term * A21(unpinched(1),lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti2r1

subroutine ti2r2(tens,r1,m1,m2,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1
   real(ki), intent(in) :: m1, m2
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4) :: term
   integer, dimension(2) :: unpinched
   unpinched = unpackb(pminus(b_ref,packb(pinches)),2)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(r1,m1*m1,m2*m2)

   call symmetric_A_coeff(term,r1,r1)
   tens(:,:) = term(:,:) * A22(unpinched(1),unpinched(1),lpinches)
   call symmetric_B_coeff(term)
   tens(:,:) = tens(:,:) + term(:,:) * B22(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti2r2
subroutine ti2r3(tens,r1,m1,m2,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1
   real(ki), intent(in) :: m1, m2
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4) :: term
   real(ki), dimension(1,4) :: rarr
   integer :: j1, j2, j3
   integer, dimension(2) :: unpinched

   rarr(1,:) = r1(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),2)

   if(f_init_smat) call init_smat(r1,m1*m1,m2*m2)

   tens(:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,2
   do j2=1,2
   do j3=1,2
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:))
   tens(:,:,:) = tens(:,:,:) &
             & + term(:,:,:) * A23(unpinched(j1),unpinched(j2),unpinched(j3),&
             & lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:))
   tens(:,:,:) = tens(:,:,:) &
             & + term(:,:,:) * B23(unpinched(j1),lpinches)
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti2r3
!---#] Two point tensor integrals :
!---#[ Three point tensor integrals :
subroutine ti3r0(tens,r1,r2,m1,m2,m3,flag,pinches)
   implicit none
   type(form_factor), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2
   real(ki), intent(in) :: m1, m2, m3
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(r1,r2,m1*m1,m2*m2,m3*m3)

   tens = A30(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti3r0

subroutine ti3r1(tens,r1,r2,m1,m2,m3,flag,pinches)
   implicit none
   type(form_factor), dimension(4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2
   real(ki), intent(in) :: m1, m2, m3
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   integer, dimension(3) :: unpinched

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),2)

   if(f_init_smat) call init_smat(r1,r2,m1*m1,m2*m2,m3*m3)

   tens(:) = r1(:) * A31(unpinched(1),lpinches) &
         & + r2(:) * A31(unpinched(2),lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti3r1

subroutine ti3r2(tens,r1,r2,m1,m2,m3,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2
   real(ki), intent(in) :: m1, m2, m3
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4) :: term
   real(ki), dimension(2,4) :: rarr
   integer :: j1, j2
   integer, dimension(3) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),3)

   if(f_init_smat) call init_smat(r1,r2,m1*m1,m2*m2,m3*m3)

   tens(:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,2
   do j2=1,2
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:))
   tens(:,:) = tens(:,:) &
           & + term(:,:) * A32(unpinched(j1),unpinched(j2),lpinches)
   end do
   end do
   call symmetric_B_coeff(term)
   tens(:,:) = tens(:,:) + term(:,:) * B32(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti3r2

subroutine ti3r3(tens,r1,r2,m1,m2,m3,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2
   real(ki), intent(in) :: m1, m2, m3
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4) :: term
   real(ki), dimension(2,4) :: rarr
   integer :: j1, j2, j3
   integer, dimension(3) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),3)

   if(f_init_smat) call init_smat(r1,r2,m1*m1,m2*m2,m3*m3)

   tens(:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,2
   do j2=1,2
   do j3=1,2
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:))
   tens(:,:,:) = tens(:,:,:) &
             & + term(:,:,:) * A33(unpinched(j1),unpinched(j2),unpinched(j3),&
             & lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:))
   tens(:,:,:) = tens(:,:,:) &
             & + term(:,:,:) * B33(unpinched(j1),lpinches)
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti3r3
subroutine ti3r4(tens,r1,r2,m1,m2,m3,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2
   real(ki), intent(in) :: m1, m2, m3
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4) :: term
   real(ki), dimension(2,4) :: rarr
   integer :: j1, j2, j3, j4
   integer, dimension(4) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),3)

   if(f_init_smat) call init_smat(r1,r2,m1*m1,m2*m2,m3*m3)

   tens(:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,3
   do j2=1,3
   do j3=1,3
   do j4=1,3
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:))
   tens(:,:,:,:) = tens(:,:,:,:) &
               & + term(:,:,:,:) * A44(unpinched(j1),unpinched(j2), &
               & unpinched(j3),unpinched(j4),lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:),rarr(j2,:))
   tens(:,:,:,:) = tens(:,:,:,:) &
               & + term(:,:,:,:) * B44(unpinched(j1),unpinched(j2),lpinches)
   end do
   end do
   call symmetric_C_coeff(term)
   tens(:,:,:,:) = tens(:,:,:,:) + term(:,:,:,:) * C44(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti3r4
!---#] Three point tensor integrals :
!---#[ Four point tensor integrals :
subroutine ti4r0(tens,r1,r2,r3,m1,m2,m3,m4,flag,pinches)
   implicit none
   type(form_factor), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3
   real(ki), intent(in) :: m1, m2, m3, m4
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(r1,r2,r3,m1*m1,m2*m2,m3*m3,m4*m4)

   tens = A40(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti4r0

subroutine ti4r1(tens,r1,r2,r3,m1,m2,m3,m4,flag,pinches)
   implicit none
   type(form_factor), dimension(4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3
   real(ki), intent(in) :: m1, m2, m3, m4
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4) :: term
   real(ki), dimension(3,4) :: rarr
   integer :: j1
   integer, dimension(4) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),4)

   if(f_init_smat) call init_smat(r1,r2,r3,m1*m1,m2*m2,m3*m3,m4*m4)

   tens(:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,3
   call symmetric_A_coeff(term,rarr(j1,:))
   tens(:) = tens(:) + term(:) * A41(unpinched(j1),lpinches)
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti4r1
subroutine ti4r2(tens,r1,r2,r3,m1,m2,m3,m4,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3
   real(ki), intent(in) :: m1, m2, m3, m4
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4) :: term
   real(ki), dimension(3,4) :: rarr
   integer :: j1, j2
   integer, dimension(4) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),4)

   if(f_init_smat) call init_smat(r1,r2,r3,m1*m1,m2*m2,m3*m3,m4*m4)

   tens(:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,3
   do j2=1,3
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:))
   tens(:,:) = tens(:,:) &
           & + term(:,:) * A42(unpinched(j1),unpinched(j2),lpinches)
   end do
   end do
   call symmetric_B_coeff(term)
   tens(:,:) = tens(:,:) + term(:,:) * B42(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti4r2

subroutine ti4r3(tens,r1,r2,r3,m1,m2,m3,m4,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3
   real(ki), intent(in) :: m1, m2, m3, m4
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4) :: term
   real(ki), dimension(3,4) :: rarr
   integer :: j1, j2, j3
   integer, dimension(4) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),4)

   if(f_init_smat) call init_smat(r1,r2,r3,m1*m1,m2*m2,m3*m3,m4*m4)

   tens(:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,3
   do j2=1,3
   do j3=1,3
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:))
   tens(:,:,:) = tens(:,:,:) &
             & + term(:,:,:) * A43(unpinched(j1),unpinched(j2),unpinched(j3),&
             & lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(unpinched(j1),:))
   tens(:,:,:) = tens(:,:,:) + term(:,:,:) * B43(j1,lpinches)
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti4r3

subroutine ti4r4(tens,r1,r2,r3,m1,m2,m3,m4,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3
   real(ki), intent(in) :: m1, m2, m3, m4
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4) :: term
   real(ki), dimension(3,4) :: rarr
   integer :: j1, j2, j3, j4
   integer, dimension(4) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),4)

   if(f_init_smat) call init_smat(r1,r2,r3,m1*m1,m2*m2,m3*m3,m4*m4)

   tens(:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,3
   do j2=1,3
   do j3=1,3
   do j4=1,3
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:))
   tens(:,:,:,:) = tens(:,:,:,:) &
               & + term(:,:,:,:) * A44(unpinched(j1),unpinched(j2), &
               & unpinched(j3),unpinched(j4),lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:),rarr(j2,:))
   tens(:,:,:,:) = tens(:,:,:,:) &
               & + term(:,:,:,:) * B44(unpinched(j1),unpinched(j2),lpinches)
   end do
   end do
   call symmetric_C_coeff(term)
   tens(:,:,:,:) = tens(:,:,:,:) + term(:,:,:,:) * C44(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti4r4
subroutine ti4r5(tens,r1,r2,r3,m1,m2,m3,m4,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3
   real(ki), intent(in) :: m1, m2, m3, m4
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4,4) :: term
   real(ki), dimension(3,4) :: rarr
   integer :: j1, j2, j3, j4, j5
   integer, dimension(4) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),4)

   if(f_init_smat) call init_smat(r1,r2,r3,m1*m1,m2*m2,m3*m3,m4*m4)

   tens(:,:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,4
   do j2=1,4
   do j3=1,4
   do j4=1,4
   do j5=1,4
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:),&
   &                           rarr(j5,:))
   tens(:,:,:,:,:) = tens(:,:,:,:,:) + term(:,:,:,:,:) * &
   & A45(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4),&
   & unpinched(j5),lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:))
   tens(:,:,:,:,:) = tens(:,:,:,:,:) + term(:,:,:,:,:) * &
   & B45(unpinched(j1),unpinched(j2),unpinched(j3),lpinches)
   end do
   end do
   call symmetric_C_coeff(term,rarr(j1,:))
   tens(:,:,:,:,:) = tens(:,:,:,:,:) + term(:,:,:,:,:) * &
   & C45(unpinched(j1),lpinches)
   end do
   if(f_dispose_smat) call done_smat()
end subroutine ti4r5
!---#] Four point tensor integrals :
!---#[ Five point tensor integrals :
subroutine ti5r0(tens,r1,r2,r3,r4,m1,m2,m3,m4,m5,flag,pinches)
   implicit none
   type(form_factor), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4
   real(ki), intent(in) :: m1, m2, m3, m4, m5
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(r1,r2,r3,r4,m1*m1,m2*m2,m3*m3,m4*m4,m5*m5)

   tens = A50(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti5r0
subroutine ti5r1(tens,r1,r2,r3,r4,m1,m2,m3,m4,m5,flag,pinches)
   implicit none
   type(form_factor), dimension(4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4
   real(ki), intent(in) :: m1, m2, m3, m4, m5
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4) :: term
   real(ki), dimension(4,4) :: rarr
   integer :: j1
   integer, dimension(5) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),5)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,m1*m1,m2*m2,m3*m3,m4*m4,m5*m5)

   tens(:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,4
   call symmetric_A_coeff(term,rarr(j1,:))
   tens(:) = tens(:) + term(:) * A51(unpinched(j1),lpinches)
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti5r1
subroutine ti5r2(tens,r1,r2,r3,r4,m1,m2,m3,m4,m5,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4
   real(ki), intent(in) :: m1, m2, m3, m4, m5
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4) :: term
   real(ki), dimension(4,4) :: rarr
   integer :: j1, j2
   integer, dimension(5) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),5)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,m1*m1,m2*m2,m3*m3,m4*m4,m5*m5)

   tens(:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,4
   do j2=1,4
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:))
   tens(:,:) = tens(:,:) &
           & + term(:,:) * A52(unpinched(j1),unpinched(j2),lpinches)
   end do
   end do
   call symmetric_B_coeff(term)
   tens(:,:) = tens(:,:) + term(:,:) * B52(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti5r2
subroutine ti5r3(tens,r1,r2,r3,r4,m1,m2,m3,m4,m5,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4
   real(ki), intent(in) :: m1, m2, m3, m4, m5
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4) :: term
   real(ki), dimension(4,4) :: rarr
   integer :: j1, j2, j3
   integer, dimension(5) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),5)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,m1*m1,m2*m2,m3*m3,m4*m4,m5*m5)

   tens(:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,4
   do j2=1,4
   do j3=1,4
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:))
   tens(:,:,:) = tens(:,:,:) &
             & + term(:,:,:) * A53(unpinched(j1),unpinched(j2), &
             & unpinched(j3),lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:))
   tens(:,:,:) = tens(:,:,:) + term(:,:,:) * B53(unpinched(j1),lpinches)
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti5r3
subroutine ti5r4(tens,r1,r2,r3,r4,m1,m2,m3,m4,m5,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4
   real(ki), intent(in) :: m1, m2, m3, m4, m5
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4) :: term
   real(ki), dimension(4,4) :: rarr
   integer :: j1, j2, j3, j4
   integer, dimension(5) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),5)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,m1*m1,m2*m2,m3*m3,m4*m4,m5*m5)

   tens(:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,4
   do j2=1,4
   do j3=1,4
   do j4=1,4
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:))
   tens(:,:,:,:) = tens(:,:,:,:) + term(:,:,:,:) * &
   & A54(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4),lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:),rarr(j2,:))
   tens(:,:,:,:) = tens(:,:,:,:) &
               & + term(:,:,:,:) * B54(unpinched(j1),unpinched(j2),lpinches)
   end do
   end do
   call symmetric_C_coeff(term)
   tens(:,:,:,:) = tens(:,:,:,:) + term(:,:,:,:) * C54(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti5r4
subroutine ti5r5(tens,r1,r2,r3,r4,m1,m2,m3,m4,m5,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4
   real(ki), intent(in) :: m1, m2, m3, m4, m5
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4,4) :: term
   real(ki), dimension(4,4) :: rarr
   integer :: j1, j2, j3, j4, j5
   integer, dimension(5) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),5)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,m1*m1,m2*m2,m3*m3,m4*m4,m5*m5)

   tens(:,:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,4
   do j2=1,4
   do j3=1,4
   do j4=1,4
   do j5=1,4
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:),&
   &                           rarr(j5,:))
   tens(:,:,:,:,:) = tens(:,:,:,:,:) + term(:,:,:,:,:) * &
   & A55(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4),&
   & unpinched(j5),lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:))
   tens(:,:,:,:,:) = tens(:,:,:,:,:) + term(:,:,:,:,:) * &
   & B55(unpinched(j1),unpinched(j2),unpinched(j3),lpinches)
   end do
   end do
   call symmetric_C_coeff(term,rarr(j1,:))
   tens(:,:,:,:,:) = tens(:,:,:,:,:) + term(:,:,:,:,:) * &
   & C55(unpinched(j1),lpinches)
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti5r5
subroutine ti5r6(tens,r1,r2,r3,r4,m1,m2,m3,m4,m5,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4
   real(ki), intent(in) :: m1, m2, m3, m4, m5
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4,4,4) :: term
   real(ki), dimension(5,4) :: rarr
   integer :: j1, j2, j3, j4, j5, j6
   integer, dimension(5) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),5)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5)

   tens(:,:,:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,5
   do j2=1,5
   do j3=1,5
   do j4=1,5
   do j5=1,5
   do j6=1,5
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:),&
   &                           rarr(j5,:),rarr(j6,:))
   tens(:,:,:,:,:,:) = tens(:,:,:,:,:,:) + term(:,:,:,:,:,:) * &
   & A66(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4),&
   & unpinched(j5),unpinched(j6),lpinches)
   end do
   end do
   call symmetric_B_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:))
   tens(:,:,:,:,:,:) = tens(:,:,:,:,:,:) + term(:,:,:,:,:,:) * &
   & B56(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4),lpinches)
   end do
   end do
   call symmetric_C_coeff(term,rarr(j1,:),rarr(j2,:))
   tens(:,:,:,:,:,:) = tens(:,:,:,:,:,:) + term(:,:,:,:,:,:) * &
   & C56(unpinched(j1),unpinched(j2),lpinches)
   end do
   end do

   call symmetric_D_coeff(term)
   tens(:,:,:,:,:,:) = tens(:,:,:,:,:,:) + term(:,:,:,:,:,:) * &
   & D56(lpinches)


   if(f_dispose_smat) call done_smat()
end subroutine ti5r6


!---#] Five point tensor integrals :
!---#[ Six point tensor integrals :
subroutine ti6r0(tens,r1,r2,r3,r4,r5,m1,m2,m3,m4,m5,m6,flag,pinches)
   implicit none
   type(form_factor), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4, r5
   real(ki), intent(in) :: m1, m2, m3, m4, m5, m6
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if

   if(f_init_smat) call init_smat(r1,r2,r3,r4,r5,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5,m6*m6)

   tens = A60(lpinches)

   if(f_dispose_smat) call done_smat()
end subroutine ti6r0

subroutine ti6r1(tens,r1,r2,r3,r4,r5,m1,m2,m3,m4,m5,m6,flag,pinches)
   implicit none
   type(form_factor), dimension(4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4, r5
   real(ki), intent(in) :: m1, m2, m3, m4, m5, m6
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4) :: term
   real(ki), dimension(5,4) :: rarr
   integer :: j1
   integer, dimension(6) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)
   rarr(5,:) = r5(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),6)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,r5,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5,m6*m6)

   tens(:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,5
   call symmetric_A_coeff(term,rarr(j1,:))
   tens(:) = tens(:) + term(:) * A61(unpinched(j1),lpinches)
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti6r1
subroutine ti6r2(tens,r1,r2,r3,r4,r5,m1,m2,m3,m4,m5,m6,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4, r5
   real(ki), intent(in) :: m1, m2, m3, m4, m5, m6
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4) :: term
   real(ki), dimension(5,4) :: rarr
   integer :: j1, j2
   integer, dimension(6) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)
   rarr(5,:) = r5(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),6)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,r5,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5,m6*m6)

   tens(:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,5
   do j2=1,5
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:))
   tens(:,:) = tens(:,:) &
           & + term(:,:) * A62(unpinched(j1),unpinched(j2),lpinches)
   end do
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti6r2
subroutine ti6r3(tens,r1,r2,r3,r4,r5,m1,m2,m3,m4,m5,m6,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4, r5
   real(ki), intent(in) :: m1, m2, m3, m4, m5, m6
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4) :: term
   real(ki), dimension(5,4) :: rarr
   integer :: j1, j2, j3
   integer, dimension(6) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)
   rarr(5,:) = r5(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),6)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,r5,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5,m6*m6)

   tens(:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,5
   do j2=1,5
   do j3=1,5
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:))
   tens(:,:,:) = tens(:,:,:) &
             & + term(:,:,:) * A63(unpinched(j1),unpinched(j2),unpinched(j3),&
             & lpinches)
   end do
   end do
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti6r3
subroutine ti6r4(tens,r1,r2,r3,r4,r5,m1,m2,m3,m4,m5,m6,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4, r5
   real(ki), intent(in) :: m1, m2, m3, m4, m5, m6
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4) :: term
   real(ki), dimension(5,4) :: rarr
   integer :: j1, j2, j3, j4
   integer, dimension(6) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)
   rarr(5,:) = r5(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),6)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,r5,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5,m6*m6)

   tens(:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,5
   do j2=1,5
   do j3=1,5
   do j4=1,5
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:))
   tens(:,:,:,:) = tens(:,:,:,:) + term(:,:,:,:) * &
   & A64(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4),lpinches)
   end do
   end do
   end do
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti6r4
subroutine ti6r5(tens,r1,r2,r3,r4,r5,m1,m2,m3,m4,m5,m6,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4, r5
   real(ki), intent(in) :: m1, m2, m3, m4, m5, m6
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4,4) :: term
   real(ki), dimension(5,4) :: rarr
   integer :: j1, j2, j3, j4, j5
   integer, dimension(6) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)
   rarr(5,:) = r5(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),6)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,r5,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5,m6*m6)

   tens(:,:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,5
   do j2=1,5
   do j3=1,5
   do j4=1,5
   do j5=1,5
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:),&
   &                           rarr(j5,:))
   tens(:,:,:,:,:) = tens(:,:,:,:,:) + term(:,:,:,:,:) * &
   & A65(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4), &
   & unpinched(j5),lpinches)
   end do
   end do
   end do
   end do
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti6r5
subroutine ti6r6(tens,r1,r2,r3,r4,r5,m1,m2,m3,m4,m5,m6,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4, r5
   real(ki), intent(in) :: m1, m2, m3, m4, m5, m6
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4,4,4) :: term
   real(ki), dimension(5,4) :: rarr
   integer :: j1, j2, j3, j4, j5, j6
   integer, dimension(6) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)
   rarr(5,:) = r5(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),6)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,r5,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5,m6*m6)

   tens(:,:,:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,5
   do j2=1,5
   do j3=1,5
   do j4=1,5
   do j5=1,5
   do j6=1,5
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:),&
   &                           rarr(j5,:),rarr(j6,:))
   tens(:,:,:,:,:,:) = tens(:,:,:,:,:,:) + term(:,:,:,:,:,:) * &
   & A66(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4),&
   & unpinched(j5),unpinched(j6),lpinches)
   end do
   end do
   end do
   end do
   end do
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti6r6
subroutine ti6r7(tens,r1,r2,r3,r4,r5,m1,m2,m3,m4,m5,m6,flag,pinches)
   implicit none
   type(form_factor), dimension(4,4,4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3, r4, r5
   real(ki), intent(in) :: m1, m2, m3, m4, m5, m6
   integer, optional, intent(in) :: flag
   integer, dimension(:), target, optional, intent(in) :: pinches
   integer, dimension(:), pointer :: lpinches
   logical :: f_init_smat, f_dispose_smat
   real(ki), dimension(4,4,4,4,4,4,4) :: term
   real(ki), dimension(5,4) :: rarr
   integer :: j1, j2, j3, j4, j5, j6, j7
   integer, dimension(6) :: unpinched

   rarr(1,:) = r1(:)
   rarr(2,:) = r2(:)
   rarr(3,:) = r3(:)
   rarr(4,:) = r4(:)
   rarr(5,:) = r5(:)

   if(present(flag)) then
      f_init_smat    = iand(flag, use_existing_smat) .eq. 0
      f_dispose_smat = iand(flag, keep_smat_on_exit) .eq. 0
   else
      f_init_smat    = .true.
      f_dispose_smat = .true.
   end if

   if (present(pinches)) then
      lpinches => pinches
   else
      lpinches => loc_s_null
   end if
   unpinched = unpackb(pminus(b_ref,packb(lpinches)),6)

   if(f_init_smat) call init_smat(r1,r2,r3,r4,r5,m1*m1,m2*m2,m3*m3,m4*m4,&
   &                              m5*m5,m6*m6)

   tens(:,:,:,:,:,:,:) = form_factor(0.0_ki, 0.0_ki, 0.0_ki)
   do j1=1,5
   do j2=1,5
   do j3=1,5
   do j4=1,5
   do j5=1,5
   do j6=1,5
   do j7=1,5
   call symmetric_A_coeff(term,rarr(j1,:),rarr(j2,:),rarr(j3,:),rarr(j4,:),&
   &                           rarr(j5,:),rarr(j6,:),rarr(j7,:))
   tens(:,:,:,:,:,:,:) = tens(:,:,:,:,:,:,:) + term(:,:,:,:,:,:,:) * &
   & A67(unpinched(j1),unpinched(j2),unpinched(j3),unpinched(j4),&
   & unpinched(j5),unpinched(j6),unpinched(j7),lpinches)
   end do
   end do
   end do
   end do
   end do
   end do
   end do

   if(f_dispose_smat) call done_smat()
end subroutine ti6r7

!---#] Six point tensor integrals :
!---#[ Symmetric A coefficients :
pure subroutine symmetric_A_coeff1(tens,r1)
   implicit none
   real(ki), dimension(4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1

   tens(:) = r1(:)
end subroutine symmetric_A_coeff1

pure subroutine symmetric_A_coeff2(tens,r1,r2)
   implicit none
   real(ki), dimension(4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1,r2

   integer :: i1,i2

   tens(:,:) = 0.0_ki

   !$omp parallel do
   do i1 = 1,4
   do i2 = 1,4
      tens(i1,i2) = r1(i1)*r2(i2)
   end do     
   end do
   !$omp end parallel do
end  subroutine symmetric_A_coeff2

pure subroutine symmetric_A_coeff3(tens,r1,r2,r3)
   implicit none
   real(ki), dimension(4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1,r2,r3

   integer :: i1,i2,i3

   !$omp parallel do
   do i1 = 1,4
   do i2 = 1,4
   do i3 = 1,4
      tens(i1,i2,i3) = r1(i1)*r2(i2)*r3(i3)
   end do     
   end do     
   end do
   !$omp end parallel do
end  subroutine symmetric_A_coeff3

pure subroutine symmetric_A_coeff4(tens,r1,r2,r3,r4)
   implicit none
   real(ki), dimension(4), intent(in) :: r1,r2,r3,r4
   real(ki), dimension(4,4,4,4), intent(out) :: tens

   integer :: i1,i2,i3,i4

   !$omp parallel do
   do i1 = 1,4
   do i2 = 1,4
   do i3 = 1,4
   do i4 = 1,4
      tens(i1,i2,i3,i4) = r1(i1)*r2(i2)*r3(i3)*r4(i4)
   end do     
   end do     
   end do     
   end do
   !$omp end parallel do
end  subroutine symmetric_A_coeff4

pure subroutine symmetric_A_coeff5(tens,r1,r2,r3,r4,r5)
   implicit none
   real(ki), dimension(4), intent(in) :: r1,r2,r3,r4,r5
   real(ki), dimension(4,4,4,4,4), intent(out) :: tens

   integer :: i1,i2,i3,i4,i5

   !$omp parallel do
   do i1 = 1,4
   do i2 = 1,4
   do i3 = 1,4
   do i4 = 1,4
   do i5 = 1,4
      tens(i1,i2,i3,i4,i5) = r1(i1)*r2(i2)*r3(i3)*r4(i4)*r5(i5)
   end do     
   end do     
   end do
   end do
   end do
   !$omp end parallel do
end  subroutine symmetric_A_coeff5

pure subroutine symmetric_A_coeff6(tens,r1,r2,r3,r4,r5,r6)
   implicit none
   real(ki), dimension(4), intent(in) :: r1,r2,r3,r4,r5,r6
   real(ki), dimension(4,4,4,4,4,4), intent(out) :: tens

   integer :: i1,i2,i3,i4,i5,i6

   !$omp parallel do
   do i1 = 1,4
   do i2 = 1,4
   do i3 = 1,4
   do i4 = 1,4
   do i5 = 1,4
   do i6 = 1,4
      tens(i1,i2,i3,i4,i5,i6) = r1(i1)*r2(i2)*r3(i3)*r4(i4)*r5(i5)*r6(i6)
   end do     
   end do     
   end do     
   end do
   end do
   end do
   !$omp end parallel do
end  subroutine symmetric_A_coeff6
pure subroutine symmetric_A_coeff7(tens,r1,r2,r3,r4,r5,r6,r7)
   implicit none
   real(ki), dimension(4), intent(in) :: r1,r2,r3,r4,r5,r6,r7
   real(ki), dimension(4,4,4,4,4,4,4), intent(out) :: tens

   integer :: i1,i2,i3,i4,i5,i6,i7

   !$omp parallel do
   do i1 = 1,4
   do i2 = 1,4
   do i3 = 1,4
   do i4 = 1,4
   do i5 = 1,4
   do i6 = 1,4
   do i7 = 1,4
      tens(i1,i2,i3,i4,i5,i6,i7) = r1(i1)*r2(i2)*r3(i3)*r4(i4)*r5(i5)&
      &                            *r6(i6)*r7(i7)
   end do     
   end do     
   end do     
   end do     
   end do
   end do
   end do
   !$omp end parallel do
end  subroutine symmetric_A_coeff7

!---#] Symmetric A coefficients :
!---#[ Symmetric B coefficients :
pure subroutine symmetric_B_coeff2(tens)
   implicit none
   real(ki), dimension(4,4), intent(out) :: tens

   tens(:,:) =  0.0_ki
   tens(1,1) =  1.0_ki
   tens(2,2) = -1.0_ki
   tens(3,3) = -1.0_ki
   tens(4,4) = -1.0_ki
end  subroutine symmetric_B_coeff2

pure subroutine symmetric_B_coeff3(tens,r1)
   implicit none
   real(ki), dimension(4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1

   real(ki), dimension(4) :: atens
   real(ki) :: term
   integer :: i1,i2,i3
   integer :: s1, s2

   call symmetric_A_coeff(atens,r1)

   !$omp parallel do
   do i1=1,4
   s1 = -sign(1, 2*i1-3)
   do i2=1,4
   s2 = -sign(1, 2*i2-3)
   do i3=1,4
      term = 0.0_ki
      if(i2==i3) term = term + s2*atens(i1)
      if(i1==i3) term = term + s1*atens(i2)
      if(i1==i2) term = term + s1*atens(i3)
      tens(i1,i2,i3) = term
   enddo
   enddo
   enddo
   !$omp end parallel do
end  subroutine symmetric_B_coeff3

pure subroutine symmetric_B_coeff4(tens,r1,r2)
   implicit none
   real(ki), dimension(4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2

   real(ki) :: term
   real(ki), dimension(4,4) :: atens
   integer :: i1,i2,i3,i4
   integer :: s1,s2,s3

   call symmetric_A_coeff(atens,r1,r2)
 
   !$omp parallel do
   do i1=1,4
   s1 = -sign(1, 2*i1-3)
   do i2=1,4
   s2 = -sign(1, 2*i2-3)
   do i3=1,4
   s3 = -sign(1, 2*i3-3)
   do i4=1,4
      term = 0.0_ki
      if(i3==i4) term = term + s3*atens(i1,i2)
      if(i2==i4) term = term + s2*atens(i1,i3)
      if(i2==i3) term = term + s2*atens(i1,i4)
      if(i1==i4) term = term + s1*atens(i2,i3)
      if(i1==i3) term = term + s1*atens(i2,i4)
      if(i1==i2) term = term + s1*atens(i3,i4)
      tens(i1,i2,i3,i4) = term
   enddo
   enddo
   enddo
   enddo
   !$omp end parallel do
end  subroutine symmetric_B_coeff4

pure subroutine symmetric_B_coeff5(tens,r1,r2,r3)
   implicit none
   real(ki), dimension(4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3

   real(ki) :: term
   real(ki), dimension(4,4,4) :: atens
   integer :: i1,i2,i3,i4,i5
   integer :: s1,s2,s3,s4

   call symmetric_A_coeff(atens,r1,r2,r3)
 
   !$omp parallel do
   do i1=1,4
   s1 = -sign(1, 2*i1-3)
   do i2=1,4
   s2 = -sign(1, 2*i2-3)
   do i3=1,4
   s3 = -sign(1, 2*i3-3)
   do i4=1,4
   s4 = -sign(1, 2*i4-3)
   do i5=1,4
      term = 0.0_ki
      if(i4==i5) term = term + s4*atens(i1,i2,i3)
      if(i3==i5) term = term + s3*atens(i1,i2,i4)
      if(i3==i4) term = term + s3*atens(i1,i2,i5)
      if(i2==i5) term = term + s2*atens(i1,i3,i4)
      if(i2==i4) term = term + s2*atens(i1,i3,i5)
      if(i2==i3) term = term + s2*atens(i1,i4,i5)
      if(i1==i5) term = term + s1*atens(i2,i3,i4)
      if(i1==i4) term = term + s1*atens(i2,i3,i5)
      if(i1==i3) term = term + s1*atens(i2,i4,i5)
      if(i1==i2) term = term + s1*atens(i3,i4,i5)
      tens(i1,i2,i3,i4,i5) = term
   enddo
   enddo
   enddo
   enddo
   enddo
   !$omp end parallel do
end  subroutine symmetric_B_coeff5
!---#] Symmetric B coefficients :
pure subroutine symmetric_B_coeff6(tens,r1,r2,r3,r4)
   implicit none
   real(ki), dimension(4,4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2, r3,r4

   real(ki) :: term
   real(ki), dimension(4,4,4,4) :: atens
   integer :: i1,i2,i3,i4,i5,i6
   integer :: s1,s2,s3,s4,s5

   call symmetric_A_coeff(atens,r1,r2,r3,r4)
 
   !$omp parallel do
   do i1=1,4
   s1 = -sign(1, 2*i1-3)
   do i2=1,4
   s2 = -sign(1, 2*i2-3)
   do i3=1,4
   s3 = -sign(1, 2*i3-3)
   do i4=1,4
   s4 = -sign(1, 2*i4-3)
   do i5=1,4
   s5 = -sign(1, 2*i5-3)
   do i6=1,4
      term = 0.0_ki
      if(i1==i2) term = term + s1*atens(i3,i4,i5,i6)
      if(i1==i3) term = term + s1*atens(i2,i4,i5,i6)
      if(i1==i4) term = term + s1*atens(i2,i3,i5,i6)
      if(i1==i5) term = term + s1*atens(i2,i3,i4,i6)
      if(i1==i6) term = term + s1*atens(i2,i3,i4,i5)
      if(i2==i3) term = term + s2*atens(i1,i4,i5,i6)
      if(i2==i4) term = term + s2*atens(i1,i3,i5,i6)
      if(i2==i5) term = term + s2*atens(i1,i3,i4,i6)
      if(i2==i6) term = term + s2*atens(i1,i3,i4,i5)
      if(i3==i4) term = term + s3*atens(i1,i2,i5,i6)
      if(i3==i5) term = term + s3*atens(i1,i2,i4,i6)
      if(i3==i6) term = term + s3*atens(i1,i2,i4,i5)
      if(i4==i5) term = term + s4*atens(i1,i2,i3,i6)
      if(i4==i6) term = term + s4*atens(i1,i2,i3,i5)
      if(i5==i6) term = term + s5*atens(i1,i2,i3,i4)
      tens(i1,i2,i3,i4,i5,i6) = term
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   !$omp end parallel do
end  subroutine symmetric_B_coeff6

!---#[ Symmetric C coefficients :
pure subroutine symmetric_C_coeff4(tens)
   implicit none
   real(ki), dimension(4,4,4,4), intent(out) :: tens

   real(ki) :: term
   integer :: i1,i2,i3,i4
   integer :: s1, s2, s3

   !$omp parallel do
   do i1=1,4
   s1 = -sign(1, 2*i1-3)
   do i2=1,4
   s2 = -sign(1, 2*i2-3)
   do i3=1,4
   s3 = -sign(1, 2*i3-3)
   do i4=1,4
      term = 0.0_ki
      if((i1==i2).and.(i3==i4)) term = term + s1*s3
      if((i1==i3).and.(i2==i4)) term = term + s1*s2
      if((i1==i4).and.(i2==i3)) term = term + s1*s2
      tens(i1,i2,i3,i4) = term
   enddo
   enddo
   enddo
   enddo
   !$omp end parallel do
end  subroutine symmetric_C_coeff4

pure subroutine symmetric_C_coeff5(tens,r1)
   implicit none
   real(ki), dimension(4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1

   real(ki), dimension(4) :: atens
   real(ki) :: term
   integer :: i1,i2,i3,i4,i5
   integer :: s1, s2, s3, s4

   call symmetric_A_coeff(atens,r1)

   !$omp parallel do
   do i1=1,4
   s1 = -sign(1, 2*i1-3)
   do i2=1,4
   s2 = -sign(1, 2*i2-3)
   do i3=1,4
   s3 = -sign(1, 2*i3-3)
   do i4=1,4
   s4 = -sign(1, 2*i4-3)
   do i5=1,4
      term = 0.0_ki
      if((i5==i2).and.(i3==i4)) term = term + s2*s3*atens(i1)
      if((i5==i3).and.(i2==i4)) term = term + s3*s2*atens(i1)
      if((i5==i4).and.(i2==i3)) term = term + s4*s2*atens(i1)
      if((i1==i5).and.(i3==i4)) term = term + s1*s3*atens(i2)
      if((i1==i3).and.(i5==i4)) term = term + s1*s4*atens(i2)
      if((i1==i4).and.(i5==i3)) term = term + s1*s3*atens(i2)
      if((i1==i2).and.(i5==i4)) term = term + s1*s4*atens(i3)
      if((i1==i5).and.(i2==i4)) term = term + s1*s2*atens(i3)
      if((i1==i4).and.(i2==i5)) term = term + s1*s2*atens(i3)
      if((i1==i2).and.(i3==i5)) term = term + s1*s3*atens(i4)
      if((i1==i3).and.(i2==i5)) term = term + s1*s2*atens(i4)
      if((i1==i5).and.(i2==i3)) term = term + s1*s2*atens(i4)
      if((i1==i2).and.(i3==i4)) term = term + s1*s3*atens(i5)
      if((i1==i3).and.(i2==i4)) term = term + s1*s2*atens(i5)
      if((i1==i4).and.(i2==i3)) term = term + s1*s2*atens(i5)
      tens(i1,i2,i3,i4,i5) = term
   enddo
   enddo
   enddo
   enddo
   enddo
   !$omp end parallel do
end  subroutine symmetric_C_coeff5
pure subroutine symmetric_C_coeff6(tens,r1,r2)
   implicit none
   real(ki), dimension(4,4,4,4,4,4), intent(out) :: tens
   real(ki), dimension(4), intent(in) :: r1, r2

   real(ki), dimension(4,4) :: atens
   real(ki) :: term
   integer :: i1,i2,i3,i4,i5, i6
   integer :: s1, s2, s3, s4, s5

   call symmetric_A_coeff(atens,r1,r2)

   !$omp parallel do
   do i1=1,4
   s1 = -sign(1, 2*i1-3)
   do i2=1,4
   s2 = -sign(1, 2*i2-3)
   do i3=1,4
   s3 = -sign(1, 2*i3-3)
   do i4=1,4
   s4 = -sign(1, 2*i4-3)
   do i5=1,4
   s5 = -sign(1, 2*i5-3)
   do i6=1,4
      term = 0.0_ki
      if((i1==i2).and.(i3==i4)) term = term + s1*s3*atens(i5,i6)
      if((i1==i2).and.(i3==i5)) term = term + s1*s3*atens(i4,i6)
      if((i1==i2).and.(i3==i6)) term = term + s1*s3*atens(i4,i5)
      if((i1==i2).and.(i4==i5)) term = term + s1*s4*atens(i3,i6)
      if((i1==i2).and.(i4==i6)) term = term + s1*s4*atens(i3,i5)
      if((i1==i2).and.(i5==i6)) term = term + s1*s5*atens(i3,i4)
      if((i1==i3).and.(i2==i4)) term = term + s1*s2*atens(i5,i6)
      if((i1==i3).and.(i2==i5)) term = term + s1*s2*atens(i4,i6)
      if((i1==i3).and.(i2==i6)) term = term + s1*s2*atens(i4,i5)
      if((i1==i3).and.(i4==i5)) term = term + s1*s4*atens(i2,i6)
      if((i1==i3).and.(i4==i6)) term = term + s1*s4*atens(i2,i5)
      if((i1==i3).and.(i5==i6)) term = term + s1*s5*atens(i2,i4)
      if((i1==i4).and.(i2==i3)) term = term + s1*s2*atens(i5,i6)
      if((i1==i4).and.(i2==i5)) term = term + s1*s2*atens(i3,i6)
      if((i1==i4).and.(i2==i6)) term = term + s1*s2*atens(i3,i5)
      if((i1==i4).and.(i3==i5)) term = term + s1*s3*atens(i2,i6)
      if((i1==i4).and.(i3==i6)) term = term + s1*s3*atens(i2,i5)
      if((i1==i4).and.(i5==i6)) term = term + s1*s5*atens(i2,i3)
      if((i1==i5).and.(i2==i3)) term = term + s1*s2*atens(i4,i6)
      if((i1==i5).and.(i2==i4)) term = term + s1*s2*atens(i3,i6)
      if((i1==i5).and.(i2==i6)) term = term + s1*s2*atens(i3,i4)
      if((i1==i5).and.(i3==i4)) term = term + s1*s3*atens(i2,i6)
      if((i1==i5).and.(i3==i6)) term = term + s1*s3*atens(i2,i4)
      if((i1==i5).and.(i4==i6)) term = term + s1*s4*atens(i2,i3)
      if((i1==i6).and.(i2==i3)) term = term + s1*s2*atens(i4,i5)
      if((i1==i6).and.(i2==i4)) term = term + s1*s2*atens(i3,i5)
      if((i1==i6).and.(i2==i5)) term = term + s1*s2*atens(i3,i4)
      if((i1==i6).and.(i3==i4)) term = term + s1*s3*atens(i2,i5)
      if((i1==i6).and.(i3==i5)) term = term + s1*s3*atens(i2,i4)
      if((i1==i6).and.(i4==i5)) term = term + s1*s4*atens(i2,i3)
      if((i2==i3).and.(i4==i5)) term = term + s2*s4*atens(i1,i6)
      if((i2==i3).and.(i4==i6)) term = term + s2*s4*atens(i1,i5)
      if((i2==i3).and.(i5==i6)) term = term + s2*s5*atens(i1,i4)
      if((i2==i4).and.(i3==i5)) term = term + s2*s3*atens(i1,i6)
      if((i2==i4).and.(i3==i6)) term = term + s2*s3*atens(i1,i5)
      if((i2==i4).and.(i5==i6)) term = term + s2*s5*atens(i1,i3)
      if((i2==i5).and.(i3==i4)) term = term + s2*s3*atens(i1,i6)
      if((i2==i5).and.(i3==i6)) term = term + s2*s3*atens(i1,i4)
      if((i2==i5).and.(i4==i6)) term = term + s2*s4*atens(i1,i3)
      if((i2==i6).and.(i3==i4)) term = term + s2*s3*atens(i1,i5)
      if((i2==i6).and.(i3==i5)) term = term + s2*s3*atens(i1,i4)
      if((i2==i6).and.(i4==i5)) term = term + s2*s4*atens(i1,i3)
      if((i3==i4).and.(i5==i6)) term = term + s3*s5*atens(i1,i2)
      if((i3==i5).and.(i4==i6)) term = term + s3*s4*atens(i1,i2)
      if((i3==i6).and.(i4==i5)) term = term + s3*s4*atens(i1,i2)
      tens(i1,i2,i3,i4,i5,i6) = term
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   !$omp end parallel do
end  subroutine symmetric_C_coeff6
!---#] Symmetric C coefficients :
!---#[ Symmetric D coefficients :
pure subroutine symmetric_D_coeff6(tens)
   implicit none
   real(ki), dimension(4,4,4,4,4,4), intent(out) :: tens

   real(ki) :: term
   integer :: i1,i2,i3,i4,i5,i6
   integer :: s1, s2, s3, s4, s5

   !$omp parallel do
   do i1=1,4
   s1 = -sign(1, 2*i1-3)
   do i2=1,4
   s2 = -sign(1, 2*i2-3)
   do i3=1,4
   s3 = -sign(1, 2*i3-3)
   do i4=1,4
   s4 = -sign(1, 2*i4-3)
   do i5=1,4
   s5 = -sign(1, 2*i5-3)
   do i6=1,4
      term = 0.0_ki
      if((i1==i2).and.(i3==i4) .and. (i5==i6)) term = term + s1*s3*s5
      if((i1==i2).and.(i3==i5) .and. (i4==i6)) term = term + s1*s3*s4
      if((i1==i2).and.(i3==i6) .and. (i4==i5)) term = term + s1*s3*s4
      if((i1==i3).and.(i2==i4) .and. (i5==i6)) term = term + s1*s2*s5
      if((i1==i3).and.(i2==i5) .and. (i4==i6)) term = term + s1*s2*s4
      if((i1==i3).and.(i2==i6) .and. (i4==i5)) term = term + s1*s2*s4
      if((i1==i4).and.(i2==i3) .and. (i5==i6)) term = term + s1*s2*s5
      if((i1==i4).and.(i2==i5) .and. (i3==i6)) term = term + s1*s2*s3
      if((i1==i4).and.(i2==i6) .and. (i3==i5)) term = term + s1*s2*s3
      if((i1==i5).and.(i2==i3) .and. (i4==i6)) term = term + s1*s2*s4
      if((i1==i5).and.(i2==i4) .and. (i3==i6)) term = term + s1*s2*s3
      if((i1==i5).and.(i2==i6) .and. (i3==i4)) term = term + s1*s2*s3
      if((i1==i6).and.(i2==i3) .and. (i4==i5)) term = term + s1*s2*s4
      if((i1==i6).and.(i2==i4) .and. (i3==i5)) term = term + s1*s2*s3
      if((i1==i6).and.(i2==i5) .and. (i3==i4)) term = term + s1*s2*s3
      tens(i1,i2,i3,i4,i5,i6) = term
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   !$omp end parallel do
end  subroutine symmetric_D_coeff6
!---#] Symmetric D coefficients :
end module tensor_integrals
