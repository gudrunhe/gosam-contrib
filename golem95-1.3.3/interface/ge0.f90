!
!****f* src/interface/ge0i
! NAME
!
!  Function ge0i
!
! USAGE
!
!  complex = ge0i(idt,s1,s2,s3,s4,s5,s12,s23,s34,s45,s51,m1,m2,m3,m4,m5,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function is the LoopTools E0i function. 
!  The first argument is a character of length <= 7
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * idt -- a character of length <= 7, the type of form factors
!  * s1 -- a real (type ki), p1^2
!  * s2 -- a real (type ki), p2^2
!  * s3 -- a real (type ki), p3^2
!  * s4 -- a real (type ki), p4^2
!  * s5 -- a real (type ki), p5^2
!  * s12 -- a real (type ki), (p1+p2)^2
!  * s23 -- a real (type ki), (p2+p3)^2
!  * s34 -- a real (type ki), (p3+p4)^2
!  * s45 -- a real (type ki), (p4+p5)^2
!  * s51 -- a real (type ki), (p5+p1)^2
!  * m1 -- a real (type ki), mass^2 of propagator 5
!  * m2 -- a real (type ki), mass^2 of propagator 1
!  * m3 -- a real (type ki), mass^2 of propagator 2
!  * m4 -- a real (type ki), mass^2 of propagator 3
!  * m5 -- a real (type ki), mass^2 of propagator 4
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
!  * form_factor_5p (src/form_factor/form_factor_5p.f90)
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
function ge0i(idt,s1,s2,s3,s4,s5,s12,s23,s34,s45,s51,m1,m2,m3,m4,m5,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_5p ! module containing the four-point form factors (export all)
  use cache, only: allocate_cache, clear_cache
  use constante, only: b_null
  use tool_lt_to_golem, only : extract
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  character (len=*), intent (in) :: idt
  real(ki), intent (in) :: s1,s2,s3,s4,s5,s12,s23,s34,s45,s51,m1,m2,m3,m4,m5,mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: ge0i
  !
  integer, dimension(7) :: tab
  integer :: n1,n2
  integer, dimension(5) :: temp
  integer :: j1,j2,j3,j4,j5
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
  character (len=7) :: idt_temp
  !
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
  call initgolem95(5)
  !
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = s23-m2-m4
  s_mat(1,4) = s51-m2-m5
  s_mat(1,5) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m4
  s_mat(2,4) = s34-m3-m5
  s_mat(2,5) = s12-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m4
  s_mat(3,4) = s4-m4-m5
  s_mat(3,5) = s45-m4-m1
  !
  s_mat(4,1) = s_mat(1,4)
  s_mat(4,2) = s_mat(2,4)
  s_mat(4,3) = s_mat(3,4)
  s_mat(4,4) = -2._ki*m5
  s_mat(4,5) = s5-m1-m5
  !
  s_mat(5,1) = s_mat(1,5)
  s_mat(5,2) = s_mat(2,5)
  s_mat(5,3) = s_mat(3,5)
  s_mat(5,4) = s_mat(4,5)
  s_mat(5,5) = -2._ki*m1
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
        ff = a50(b_null)
        !
      case(1)
        !
        ff = a51(j1,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function ge0i (lt_to_golem.f90)'
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
        ff = b52(b_null)
        !
      case(2)
        !
        ff = a52(j1,j2,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function ge0i (lt_to_golem.f90)'
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
        ff = b53(j3,b_null)
        !
      case(3)
        !
        ff = a53(j1,j2,j3,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function ge0i (lt_to_golem.f90)'
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
        ff = c54(b_null)
        !
      case(2)
        !
        ff = b54(j3,j4,b_null)
        !
      case(4)
        !
        ff = a54(j1,j2,j3,j4,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function ge0i (lt_to_golem.f90)'
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
        ff = c55(j5,b_null)
        !
      case(3)
        !
        ff = b55(j3,j4,j5,b_null)
        !
      case(5)
        !
        ff = a55(j1,j2,j3,j4,j5,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function ge0i (lt_to_golem.f90)'
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
      tab_erreur_par(1)%chaine = 'In function ge0i (lt_to_golem.f90)'
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
    !  mu2_scale_par contained in basis integrals
      !ge0i = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
      ge0i = ff%c
    !
  else if (eps_flag == -1) then
    !
    ! ge0i = ff%b + lmu2*ff%a 
    ge0i = ff%b
    !
  else if (eps_flag == -2) then
    !
    ge0i = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function ge0i (ge0.f90)'
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
 end function ge0i
!
!****f* src/interface/ge0
! NAME
!
!  Function ge0
!
! USAGE
!
!  complex = ge0(s1,s2,s3,s4,s5,s12,s23,s34,s45,s51,m1,m2,m3,m4,m5,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function is the scalar E0 function. 
!  The last two arguments are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * s1 -- a real (type ki), p1^2
!  * s2 -- a real (type ki), p2^2
!  * s3 -- a real (type ki), p3^2
!  * s4 -- a real (type ki), p4^2
!  * s5 -- a real (type ki), p5^2
!  * s12 -- a real (type ki), (p1+p2)^2
!  * s23 -- a real (type ki), (p2+p3)^2
!  * s34 -- a real (type ki), (p3+p4)^2
!  * s45 -- a real (type ki), (p4+p5)^2
!  * s51 -- a real (type ki), (p5+p1)^2
!  * m1 -- a real (type ki), mass^2 of propagator 5
!  * m2 -- a real (type ki), mass^2 of propagator 1
!  * m3 -- a real (type ki), mass^2 of propagator 2
!  * m4 -- a real (type ki), mass^2 of propagator 3
!  * m5 -- a real (type ki), mass^2 of propagator 4
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
!  * form_factor_5p (src/form_factor/form_factor_5p.f90)
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
function ge0(s1,s2,s3,s4,s5,s12,s23,s34,s45,s51,m1,m2,m3,m4,m5,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_5p ! module containing the four-point form factors (export all)
  use cache, only: allocate_cache, clear_cache
  use constante, only: b_null
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  real(ki), intent (in) :: s1,s2,s3,s4,s5,s12,s23,s34,s45,s51,m1,m2,m3,m4,m5,mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: ge0
  !
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
  !
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
  !
  call initgolem95(5)
  !
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s2-m2-m3
  s_mat(1,3) = s23-m2-m4
  s_mat(1,4) = s51-m2-m5
  s_mat(1,5) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m3
  s_mat(2,3) = s3-m3-m4
  s_mat(2,4) = s34-m3-m5
  s_mat(2,5) = s12-m3-m1
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2._ki*m4
  s_mat(3,4) = s4-m4-m5
  s_mat(3,5) = s45-m4-m1
  !
  s_mat(4,1) = s_mat(1,4)
  s_mat(4,2) = s_mat(2,4)
  s_mat(4,3) = s_mat(3,4)
  s_mat(4,4) = -2._ki*m5
  s_mat(4,5) = s5-m1-m5
  !
  s_mat(5,1) = s_mat(1,5)
  s_mat(5,2) = s_mat(2,5)
  s_mat(5,3) = s_mat(3,5)
  s_mat(5,4) = s_mat(4,5)
  s_mat(5,5) = -2._ki*m1
  !
  call preparesmatrix()
  !
        !
    ff = a50(b_null)
        !
   lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
    ! ge0 = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
    ge0 = ff%c
    !
  else if (eps_flag == -1) then
    !
   ! ge0 = ff%b + lmu2*ff%a 
   ge0 = ff%b
    !
  else if (eps_flag == -2) then
    !
    ge0 = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function ge0 (ge0.f90)'
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
  !
end function ge0
