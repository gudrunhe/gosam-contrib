!
!****f* src/interface/gc0i
! NAME
!
!  Function gc0i
!
! USAGE
!
!  complex = gc0i(idt,s1,s2,s3,m1,m2,m3,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function corresponds to the LoopTools C0i function. 
!  The first argument is a character of length <= 5
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * idt -- a character of length <= 5, the type of form factors
!  * s1 -- a real (type ki), p1^2
!  * s2 -- a real (type ki), p2^2
!  * s3 -- a real (type ki), p3^2
!  * m1 -- a real (type ki), mass of propagator 3
!  * m2 -- a real (type ki), mass of propagator 1
!  * m3 -- a real (type ki), mass of propagator 2
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
!  * form_factor_3p (src/form_factor/form_factor_3p.f90)
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
function gc0i(idt,s1,s2,s3,m1,m2,m3,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_3p 
  use form_factor_higher_ranks
  use cache, only: allocate_cache, clear_cache
  use constante, only: b_null
  use tool_lt_to_golem, only : extract
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use array, only: packb
  use parametre, only: mu2_scale_par
  implicit none
  !
  character (len=*), intent (in) :: idt
  real(ki), intent (in) :: s1,s2,s3,m1,m2,m3,mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gc0i
  !
  integer, dimension(5) :: tab
  integer :: n1,n2
  integer, dimension(5) :: temp
  integer :: j1,j2,j3,j4
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
  character (len=5) :: idt_temp
  !
  ! to avoid confusion if mu2_scale_par is not=1 in main program: 
  !   
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
  call initgolem95(3)
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m1
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
        ff = a30(b_null)
        !
      case(1)
        !
        ff = a31(j1,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gc0i (lt_to_golem.f90)'
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
        ff = b32(b_null)
        !
      case(2)
        !
        ff = a32(j1,j2,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gc0i (lt_to_golem.f90)'
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
        ff = b33(j3,b_null)
        !
      case(3)
        !
        ff = a33(j1,j2,j3,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gc0i (lt_to_golem.f90)'
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
        ff = c34(b_null)
        !
      case(2)
        !
        ff = b34(j3,j4,b_null)
        !
      case(4)
        !
        ff = a34(j1,j2,j3,j4,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gc0i (lt_to_golem.f90)'
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
      tab_erreur_par(1)%chaine = 'In function gc0i (lt_to_golem.f90)'
      tab_erreur_par(2)%a_imprimer = .true.
      tab_erreur_par(2)%chaine = 'The value of n1 is not correct %d0'
      tab_erreur_par(2)%arg_int = n1
      call catch_exception(0)
      !
      stop
      !
  end select
  !
  lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
   ! mu2_scale_par already contained in basic triangle functions
    ! note that local mu2 is used here
   ! gc0i = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
    gc0i = ff%c
    !
  else if (eps_flag == -1) then
    !
    ! gc0i = ff%b + lmu2*ff%a 
    gc0i = ff%b
    !
  else if (eps_flag == -2) then
    !
    gc0i = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gc0i (gc0.f90)'
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
 end function gc0i
!
!****f* src/interface/gc0iC
! NAME
!
!  Function gc0iC
!
! USAGE
!
!  complex = gc0iC(idt,s1,s2,s3,m1,m2,m3,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function corresponds to the LoopTools C0i function. 
!  The first argument is a character of length <= 5
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * idt -- a character of length <= 5, the type of form factors
!  * s1 -- a real (type ki), p1^2
!  * s2 -- a real (type ki), p2^2
!  * s3 -- a real (type ki), p3^2
!  * m1 -- a real (type ki), mass of propagator 3
!  * m2 -- a real (type ki), mass of propagator 1
!  * m3 -- a real (type ki), mass of propagator 2
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
!  * form_factor_3p (src/form_factor/form_factor_3p.f90)
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
function gc0iC(idt,s1,s2,s3,m1,m2,m3,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_3p 
  use form_factor_higher_ranks
  use cache, only: allocate_cache, clear_cache
  use constante, only: b_null
  use tool_lt_to_golem, only : extract
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use array, only: packb
  use parametre, only: mu2_scale_par
  implicit none
  !
  character (len=*), intent (in) :: idt
  complex(ki), intent (in) :: s1,s2,s3,m1,m2,m3
  real(ki), intent (in) :: mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gc0iC
  !
  integer, dimension(5) :: tab
  integer :: n1,n2
  integer, dimension(5) :: temp
  integer :: j1,j2,j3,j4
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
  character (len=5) :: idt_temp
  !
  ! to avoid confusion if mu2_scale_par is not=1 in main program: 
  !   
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
  call initgolem95(3)
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m1
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
        ff = a30(b_null)
        !
      case(1)
        !
        ff = a31(j1,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gc0iC (lt_to_golem.f90)'
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
        ff = b32(b_null)
        !
      case(2)
        !
        ff = a32(j1,j2,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gc0iC (lt_to_golem.f90)'
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
        ff = b33(j3,b_null)
        !
      case(3)
        !
        ff = a33(j1,j2,j3,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gc0iC (lt_to_golem.f90)'
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
        ff = c34(b_null)
        !
      case(2)
        !
        ff = b34(j3,j4,b_null)
        !
      case(4)
        !
        ff = a34(j1,j2,j3,j4,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gc0iC (lt_to_golem.f90)'
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
      tab_erreur_par(1)%chaine = 'In function gc0iC (lt_to_golem.f90)'
      tab_erreur_par(2)%a_imprimer = .true.
      tab_erreur_par(2)%chaine = 'The value of n1 is not correct %d0'
      tab_erreur_par(2)%arg_int = n1
      call catch_exception(0)
      !
      stop
      !
  end select
  !
  lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
   ! mu2_scale_par already contained in basic triangle functions
    ! note that local mu2 is used here
   ! gc0iC = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
    gc0iC = ff%c
    !
  else if (eps_flag == -1) then
    !
    ! gc0iC = ff%b + lmu2*ff%a 
    gc0iC = ff%b
    !
  else if (eps_flag == -2) then
    !
    gc0iC = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gc0iC (gc0.f90)'
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
 end function gc0iC
!
!****f* src/interface/gc0
! NAME
!
!  Function gc0
!
! USAGE
!
!  complex = gc0(s1,s2,s3,m1,m2,m3,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function corresponds to the scalar C0 function. 
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * s1 -- a real (type ki), p1^2
!  * s2 -- a real (type ki), p2^2
!  * s3 -- a real (type ki), p3^2
!  * m1 -- a real (type ki), mass^2 of propagator 3
!  * m2 -- a real (type ki), mass^2 of propagator 1
!  * m3 -- a real (type ki), mass^2 of propagator 2
!  * mu2 -- a real (type ki), renormalisation scale squared
!  * eps_flag -- an integer, a flag to select the pole coefficient 
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
!  * form_factor_3p (src/form_factor/form_factor_3p.f90)
!  * cache, only: allocate_cache, clear_cache (src/module/cache.f90)
!  * constante, only : b_null (src/module/constante.f90)
!  * sortie_erreur, only : tab_erreur_par,catch_exception (src/module/sortie_erreur.f90)
!
!
! EXAMPLE
!
!
!
!*****
function gc0(s1,s2,s3,m1,m2,m3,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type
  use form_factor_3p ! module containing the three-point form factors (export all)
  use cache, only: allocate_cache, clear_cache
  use constante
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use array, only: packb  
  use parametre, only: mu2_scale_par
  implicit none
  !
  real(ki), intent (in) :: s1,s2,s3,m1,m2,m3,mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gc0
  !
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
  !
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
 call initgolem95(3)
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m1
  !
  call preparesmatrix()
  !
   ff = a30(s_null)
  !
  lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
    !gc0 = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
   ! mu scale already contained in basic triangle functions
   ! but here use  local mu2
    gc0 = ff%c 
    !
  else if (eps_flag == -1) then
    !
   ! gc0 = ff%b + lmu2*ff%a 
    gc0 = ff%b
    !
  else if (eps_flag == -2) then
    !
    gc0 = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gc0 (gc0.f90)'
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
 end function gc0
!
!****f* src/interface/gc0c
! NAME
!
!  Function gc0c
!
! USAGE
!
!  complex = gc0c(s1,s2,s3,m1,m2,m3,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function corresponds to the scalar C0 function. 
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * s1 -- a complex (type ki), p1^2
!  * s2 -- a complex (type ki), p2^2
!  * s3 -- a complex (type ki), p3^2
!  * m1 -- a complex (type ki), mass^2 of propagator 3
!  * m2 -- a complex (type ki), mass^2 of propagator 1
!  * m3 -- a complex (type ki), mass^2 of propagator 2
!  * mu2 -- a real (type ki), renormalisation scale squared
!  * eps_flag -- an integer, a flag to select the pole coefficient 
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
!  * form_factor_3p (src/form_factor/form_factor_3p.f90)
!  * cache, only: allocate_cache, clear_cache (src/module/cache.f90)
!  * constante, only : b_null (src/module/constante.f90)
!  * sortie_erreur, only : tab_erreur_par,catch_exception (src/module/sortie_erreur.f90)
!
!
! EXAMPLE
!
!
!
!*****
function gc0c(s1,s2,s3,m1,m2,m3,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type
  use form_factor_3p ! module containing the three-point form factors (export all)
  use cache, only: allocate_cache, clear_cache
  use constante
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use array, only: packb  
  use parametre, only: mu2_scale_par
  implicit none
  !
  complex(ki), intent (in) :: s1,s2,s3,m1,m2,m3
  real(ki), intent (in) :: mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gc0c
  !
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
  !
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
 call initgolem95(3)
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m1
  !
  call preparesmatrix()
  !
   ff = a30(s_null)
  !
  lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
    !gc0c = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
   ! mu scale already contained in basic triangle functions
   ! but here use  local mu2
    gc0c = ff%c 
    !
  else if (eps_flag == -1) then
    !
   ! gc0c = ff%b + lmu2*ff%a 
    gc0c = ff%b
    !
  else if (eps_flag == -2) then
    !
    gc0c = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gc0c (gc0.f90)'
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
 end function gc0c
