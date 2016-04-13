!
!****f* src/interface/gd0i
! NAME
!
!  Function gd0i
!
! USAGE
!
!  complex = gd0i(idt,s1,s2,s3,s4,s,t,m1,m2,m3,m4,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function is the LoopTools D0i function. 
!  The first argument is a character of length <= 6
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * idt -- a character of length <= 6, the type of form factors
!  * s1 -- a real (type ki), p1^2
!  * s2 -- a real (type ki), p2^2
!  * s3 -- a real (type ki), p3^2
!  * s4 -- a real (type ki), p4^2
!  * s -- a real (type ki), (p1+p2)^2
!  * t -- a real (type ki), (p2+p3)^2
!  * m1 -- a real (type ki), mass^2 of propagator 4
!  * m2 -- a real (type ki), mass^2 of propagator 1
!  * m3 -- a real (type ki), mass^2 of propagator 2
!  * m4 -- a real (type ki), mass^2 of propagator 3
!  * mu2 -- a real (type ki), renormalisation scale squared
!  * eps_flag -- an integer, a flag to select the coefficient in front the power of epsilon
!
! SIDE EFFECTS
!
!  No side effect
!
! RETURN VALUE
!
!  It returns a complex (type ki) corresponding
!  to the real part, imaginary part of the coefficient in front 1/epsilon^2 (eps_flag=-2),
!  the real part, imaginary part of the 1/epsilon term (eps_flag=-1) and the real part,
!  imaginary part of the constant term (eps_flag=0).
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * matrice_s, only: set_ref, s_mat, allocation_s, deallocation_s, init_invs (src/kinematic/matrice_s.f90)
!  * form_factor_type, only: form_factor (src/module/form_factor_type.f90)
!  * form_factor_4p (src/form_factor/form_factor_4p.f90)
!  * form_factor_higher_ranks (src/form_factor/form_factor_higher_ranks.f90)
!  * cache, only: allocate_cache, clear_cache (src/module/cache.f90)
!  * constante, only : b_null (src/module/constante.f90)
!  * tool_lt_to_golem, only : extract (src/interface/tool_lt_to_golem.f90)
!  * sortie_erreur, only : tab_erreur_par,catch_exception (src/module/sortie_erreur.f90)
!
!
! EXAMPLE
!
!
!
!*****
function gd0i(idt,s1,s2,s3,s4,s,t,m1,m2,m3,m4,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_4p ! module containing the four-point form factors (export all)
  use form_factor_higher_ranks
  use cache, only: allocate_cache, clear_cache
  use constante, only: b_null
  use tool_lt_to_golem, only : extract
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  character (len=*), intent (in) :: idt
  real(ki), intent (in) :: s1,s2,s3,s4,s,t,m1,m2,m3,m4,mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gd0i
  !
  integer, dimension(6) :: tab
  integer :: n1,n2
  integer, dimension(6) :: temp
  integer :: j1,j2,j3,j4,j5
  type(form_factor) :: ff
  ! real(ki) :: lmu2
  real(ki) :: mu2store
  character (len=6) :: idt_temp
  !
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
  !
 call initgolem95(4)
  !
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = t-m2-m4
  s_mat(1,4) = s1-m2-m1
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m4
  s_mat(2,4) = s-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m4
  s_mat(3,4) = s4-m4-m1
  !
  s_mat(4,1) = s_mat(1,4)
  s_mat(4,2) = s_mat(2,4)
  s_mat(4,3) = s_mat(3,4)
  s_mat(4,4) = -2._ki*m1 
  !
  call preparesmatrix()
  !
  idt_temp = idt ! on complete la chaine entrante par des blancs
  call extract(idt_temp,tab)
  n1= count(tab/=-1)
  n2= count(tab/=-1 .and. tab/=0)
  !
  temp = -2
  temp(1:n1) = pack(tab,tab/=-1)
  !
  select case(n1)
    !
    case(1)
      !
      j1 = temp(1)
      !
      select case(n2)
      !
      case(0)
        !
        ff = a40(b_null)
        !
      case(1)
        !
        ff = a41(j1,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gd0i (lt_to_golem.f90)'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'The value of n2 %d0'
        tab_erreur_par(2)%arg_int = n2
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'is not compatible with n1 %d0'
        tab_erreur_par(3)%arg_int = n1
        call catch_exception(0)
        !
        stop
        !
      end select
      !
    case(2)
      !
      j1 = temp(1)
      j2 = temp(2)
      !
      select case(n2)
      !
      case(0)
        !
        ff = b42(b_null)
        !
      case(2)
        !
        ff = a42(j1,j2,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gd0i (lt_to_golem.f90)'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'The value of n2 %d0'
        tab_erreur_par(2)%arg_int = n2
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'is not compatible with n1 %d0'
        tab_erreur_par(3)%arg_int = n1
        call catch_exception(0)
        !
        stop
        !
      end select
      !
    case(3)
      !
      j1 = temp(1)
      j2 = temp(2)
      j3 = temp(3)
      !
      select case(n2)
      !
      case(1)
        !
        ff = b43(j3,b_null)
        !
      case(3)
        !
        ff = a43(j1,j2,j3,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gd0i (lt_to_golem.f90)'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'The value of n2 %d0'
        tab_erreur_par(2)%arg_int = n2
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'is not compatible with n1 %d0'
        tab_erreur_par(3)%arg_int = n1
        call catch_exception(0)
        !
        stop
        !
      end select
      !
    case(4)
      !
      j1 = temp(1)
      j2 = temp(2)
      j3 = temp(3)
      j4 = temp(4)
      !
      select case(n2)
      !
      case(0)
        !
        ff = c44(b_null)
        !
      case(2)
        !
        ff = b44(j3,j4,b_null)
        !
      case(4)
        !
        ff = a44(j1,j2,j3,j4,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gd0i (lt_to_golem.f90)'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'The value of n2 %d0'
        tab_erreur_par(2)%arg_int = n2
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'is not compatible with n1 %d0'
        tab_erreur_par(3)%arg_int = n1
        call catch_exception(0)
        !
        stop
        !
      end select
      !
    case(5)
      !
      j1 = temp(1)
      j2 = temp(2)
      j3 = temp(3)
      j4 = temp(4)
      j5 = temp(5)
      !
      select case(n2)
      !
      case(1)
        !
        ff = c45(j5,b_null)
        !
      case(3)
        !
        ff = b45(j3,j4,j5,b_null)
      case(5)
        !
        ff = a45(j1,j2,j3,j4,j5,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gd0i (lt_to_golem.f90)'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'The value of n2 %d0'
        tab_erreur_par(2)%arg_int = n2
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'is not compatible with n1 %d0'
        tab_erreur_par(3)%arg_int = n1
        call catch_exception(0)
        !
        stop
        !
      end select
      !
    case default
      !
      tab_erreur_par(1)%a_imprimer = .true.
      tab_erreur_par(1)%chaine = 'In function gd0i (lt_to_golem.f90)'
      tab_erreur_par(2)%a_imprimer = .true.
      tab_erreur_par(2)%chaine = 'The value of n1 is not correct %d0'
      tab_erreur_par(2)%arg_int = n1
      call catch_exception(0)
      !
      stop
      !
  end select
  !
  ! lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
    ! gd0i = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
    gd0i = ff%c 
    !
  else if (eps_flag == -1) then
    !
    ! gd0i = ff%b + lmu2*ff%a 
    gd0i = ff%b 
    !
  else if (eps_flag == -2) then
    !
    gd0i = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gd0i (gd0.f90)'
    tab_erreur_par(2)%a_imprimer = .true.
    tab_erreur_par(2)%chaine = 'eps_flag should be -2, -1 or 0 but is %d0'
    tab_erreur_par(2)%arg_int = eps_flag
    call catch_exception(0)
    !
    stop
    !
  end if
  !
  mu2_scale_par=mu2store
  !
  call exitgolem95()
  !  
 end function gd0i
!
!
!****f* src/interface/gd0
! NAME
!
!  Function gd0
!
! USAGE
!
!  complex = gd0(s1,s2,s3,s4,s,t,m1,m2,m3,m4,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function is the scalar LoopTools D0 function. 
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * s1 -- a real (type ki), p1^2
!  * s2 -- a real (type ki), p2^2
!  * s3 -- a real (type ki), p3^2
!  * s4 -- a real (type ki), p4^2
!  * s -- a real (type ki), (p1+p2)^2
!  * t -- a real (type ki), (p2+p3)^2
!  * m1 -- a real (type ki), mass^2 of propagator 4
!  * m2 -- a real (type ki), mass^2 of propagator 1
!  * m3 -- a real (type ki), mass^2 of propagator 2
!  * m4 -- a real (type ki), mass^2 of propagator 3
!  * mu2 -- a real (type ki), renormalisation scale squared
!  * eps_flag -- an integer, a flag to select the coefficient in front the power of epsilon
!
! SIDE EFFECTS
!
!  No side effect
!
! RETURN VALUE
!
!  It returns a complex (type ki) corresponding
!  to the real part, imaginary part of the coefficient in front 1/epsilon^2 (eps_flag=-2),
!  the real part, imaginary part of the 1/epsilon term (eps_flag=-1) and the real part,
!  imaginary part of the constant term (eps_flag=0).
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * matrice_s, only: set_ref, s_mat, allocation_s, deallocation_s, init_invs (src/kinematic/matrice_s.f90)
!  * form_factor_type, only: form_factor (src/module/form_factor_type.f90)
!  * form_factor_4p (src/form_factor/form_factor_4p.f90)
!  * cache, only: allocate_cache, clear_cache (src/module/cache.f90)
!  * constante, only : b_null (src/module/constante.f90)
!  * sortie_erreur, only : tab_erreur_par,catch_exception (src/module/sortie_erreur.f90)
!
!
!*****
function gd0(s1,s2,s3,s4,s,t,m1,m2,m3,m4,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_4p ! module containing the four-point form factors (export all)
  use cache, only: allocate_cache, clear_cache
  use constante, only: b_null
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  real(ki), intent (in) :: s1,s2,s3,s4,s,t,m1,m2,m3,m4,mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gd0
  !
  type(form_factor) :: ff
  ! real(ki) :: lmu2
  real(ki) :: mu2store
  !
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
  !
  call initgolem95(4)
  !
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = t-m2-m4
  s_mat(1,4) = s1-m2-m1
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m4
  s_mat(2,4) = s-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m4
  s_mat(3,4) = s4-m4-m1
  !
  s_mat(4,1) = s_mat(1,4)
  s_mat(4,2) = s_mat(2,4)
  s_mat(4,3) = s_mat(3,4)
  s_mat(4,4) = -2._ki*m1 
  !
  call preparesmatrix()
  !
  ff = a40(b_null)
  !
  ! lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
    ! gd0 = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
    gd0 = ff%c
    !
  else if (eps_flag == -1) then
    !
    ! gd0 = ff%b + lmu2*ff%a 
    gd0 = ff%b 
    !
  else if (eps_flag == -2) then
    !
    gd0 = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gd0 (gd0.f90)'
    tab_erreur_par(2)%a_imprimer = .true.
    tab_erreur_par(2)%chaine = 'eps_flag should be -2, -1 or 0 but is %d0'
    tab_erreur_par(2)%arg_int = eps_flag
    call catch_exception(0)
    !
    stop
    !
  end if
  !
  mu2_scale_par=mu2store
  !
  call exitgolem95()
  !  
end function gd0
!
!****f* src/interface/gd0c
! NAME
!
!  Function gd0c
!
! USAGE
!
!  complex = gd0c(s1,s2,s3,s4,s,t,m1,m2,m3,m4,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function is the scalar LoopTools D0 function. 
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * s1 -- a complex (type ki), p1^2
!  * s2 -- a complex (type ki), p2^2
!  * s3 -- a complex (type ki), p3^2
!  * s4 -- a complex (type ki), p4^2
!  * s -- a complex (type ki), (p1+p2)^2
!  * t -- a complex (type ki), (p2+p3)^2
!  * m1 -- a complex (type ki), mass^2 of propagator 4
!  * m2 -- a complex (type ki), mass^2 of propagator 1
!  * m3 -- a complex (type ki), mass^2 of propagator 2
!  * m4 -- a complex (type ki), mass^2 of propagator 3
!  * mu2 -- a real (type ki), renormalisation scale squared
!  * eps_flag -- an integer, a flag to select the coefficient in front the power of epsilon
!
! SIDE EFFECTS
!
!  No side effect
!
! RETURN VALUE
!
!  It returns a complex (type ki) corresponding
!  to the real part, imaginary part of the coefficient in front 1/epsilon^2 (eps_flag=-2),
!  the real part, imaginary part of the 1/epsilon term (eps_flag=-1) and the real part,
!  imaginary part of the constant term (eps_flag=0).
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * matrice_s, only: set_ref, s_mat, allocation_s, deallocation_s, init_invs (src/kinematic/matrice_s.f90)
!  * form_factor_type, only: form_factor (src/module/form_factor_type.f90)
!  * form_factor_4p (src/form_factor/form_factor_4p.f90)
!  * cache, only: allocate_cache, clear_cache (src/module/cache.f90)
!  * constante, only : b_null (src/module/constante.f90)
!  * sortie_erreur, only : tab_erreur_par,catch_exception (src/module/sortie_erreur.f90)
!
!
!*****
function gd0c(s1,s2,s3,s4,s,t,m1,m2,m3,m4,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_4p ! module containing the four-point form factors (export all)
  use cache, only: allocate_cache, clear_cache
  use constante, only: b_null
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  complex(ki), intent (in) :: s1,s2,s3,s4,s,t,m1,m2,m3,m4
  real(ki), intent (in) :: mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gd0c
  !
  type(form_factor) :: ff
  ! real(ki) :: lmu2
  real(ki) :: mu2store
  !
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
  !
  call initgolem95(4)
  !
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = t-m2-m4
  s_mat(1,4) = s1-m2-m1
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m4
  s_mat(2,4) = s-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m4
  s_mat(3,4) = s4-m4-m1
  !
  s_mat(4,1) = s_mat(1,4)
  s_mat(4,2) = s_mat(2,4)
  s_mat(4,3) = s_mat(3,4)
  s_mat(4,4) = -2._ki*m1 
  !
  call preparesmatrix()
  !
  ff = a40(b_null)
  !
  ! lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
    ! gd0c = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
    gd0c = ff%c
    !
  else if (eps_flag == -1) then
    !
    ! gd0c = ff%b + lmu2*ff%a 
    gd0c = ff%b 
    !
  else if (eps_flag == -2) then
    !
    gd0c = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gd0c (gd0.f90)'
    tab_erreur_par(2)%a_imprimer = .true.
    tab_erreur_par(2)%chaine = 'eps_flag should be -2, -1 or 0 but is %d0'
    tab_erreur_par(2)%arg_int = eps_flag
    call catch_exception(0)
    !
    stop
    !
  end if
  !
  mu2_scale_par=mu2store
  !
  call exitgolem95()
  !  
end function gd0c
