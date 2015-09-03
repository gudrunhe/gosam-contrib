! 
!****f* src/interface/gb0i
! NAME
!
!  Function gb0i
!
! USAGE
!
!  complex = gb0i(idt,s1,m1,m2,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function corresponds to the LoopTools b0i function.
!  The first argument is a character of length <= 4
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * idt -- a character of length <= 4, the type of form factors
!  * s1 -- a real (type ki), p1^2
!  * m1 -- a real (type ki), mass of propagator 2
!  * m2 -- a real (type ki), mass of propagator 1
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
!  * form_factor_2p (src/form_factor/form_factor_2p.f90)
!  * form_factor_higher_ranks (src/form_factor/form_factor_higher_ranks.f90)
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
function gb0i(idt,s1,m1,m2,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_2p 
  use form_factor_higher_ranks
  use constante
  use tool_lt_to_golem, only : extract
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  character (len=*), intent (in) :: idt
  real(ki), intent (in) :: s1,m1,m2,mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gb0i
  !
  integer, dimension(4) :: tab
  integer :: n1,n2
  integer, dimension(4) :: temp
  integer :: j1,j2,j3
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
  character (len=4) :: idt_temp
  !
  ! to avoid confusion if mu2_scale_par is not=1 in main program: 
  !   
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
  call initgolem95(2)
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m1
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
        ff = a20(b_null)
        !
      case(1)
        !
        ff = a21(j1,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gb0i (lt_to_golem.f90)'
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
        ff = b22(b_null)
        !
      case(2)
        !
        ff = a22(j1,j2,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gb0i (lt_to_golem.f90)'
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
        ff = b23(j3,b_null)
        !
      case(3)
        !
        ff = a23(j1,j2,j3,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gb0i (lt_to_golem.f90)'
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
      tab_erreur_par(1)%chaine = 'In function gb0i (lt_to_golem.f90)'
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
    !gb0i = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
    ! mu2_scale_par scale contained already in basic function,
    ! but local mu2 is used here
      gb0i = ff%c
   !
  else if (eps_flag == -1) then
    !
    !  gb0i = ff%b + lmu2*ff%a 
    gb0i = ff%b
    !
  else if (eps_flag == -2) then
    !
    gb0i = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gb0i (gb0.f90)'
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
end function gb0i
! 
!****f* src/interface/gb0iC
! NAME
!
!  Function gb0iC
!
! USAGE
!
!  complex = gb0iC(idt,s1,m1,m2,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function corresponds to the LoopTools b0i function.
!  The first argument is a character of length <= 4
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * idt -- a character of length <= 4, the type of form factors
!  * s1 -- a real (type ki), p1^2
!  * m1 -- a real (type ki), mass of propagator 2
!  * m2 -- a real (type ki), mass of propagator 1
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
!  * form_factor_2p (src/form_factor/form_factor_2p.f90)
!  * form_factor_higher_ranks (src/form_factor/form_factor_higher_ranks.f90)
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
function gb0iC(idt,s1,m1,m2,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_2p 
  use form_factor_higher_ranks
  use constante
  use tool_lt_to_golem, only : extract
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  character (len=*), intent (in) :: idt
  complex(ki), intent (in) :: s1,m1,m2
  real(ki), intent (in) :: mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gb0iC
  !
  integer, dimension(4) :: tab
  integer :: n1,n2
  integer, dimension(4) :: temp
  integer :: j1,j2,j3
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
  character (len=4) :: idt_temp
  !
  ! to avoid confusion if mu2_scale_par is not=1 in main program: 
  !   
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
  call initgolem95(2)
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m1
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
        ff = a20(b_null)
        !
      case(1)
        !
        ff = a21(j1,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gb0iC (lt_to_golem.f90)'
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
        ff = b22(b_null)
        !
      case(2)
        !
        ff = a22(j1,j2,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gb0iC (lt_to_golem.f90)'
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
        ff = b23(j3,b_null)
        !
      case(3)
        !
        ff = a23(j1,j2,j3,b_null)
        !
      case default
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function gb0iC (lt_to_golem.f90)'
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
      tab_erreur_par(1)%chaine = 'In function gb0iC (lt_to_golem.f90)'
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
    !gb0iC = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
    ! mu2_scale_par scale contained already in basic function,
    ! but local mu2 is used here
      gb0iC = ff%c
   !
  else if (eps_flag == -1) then
    !
    !  gb0iC = ff%b + lmu2*ff%a 
    gb0iC = ff%b
    !
  else if (eps_flag == -2) then
    !
    gb0iC = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gb0iC (gb0.f90)'
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
end function gb0iC









! 
!****f* src/interface/gb0
! NAME
!
!  Function gb0
!
! USAGE
!
!  complex = gb0(s1,m1,m2,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function corresponds to the scalar b0 function.
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * s1 -- a real (type ki), p1^2
!  * m1 -- a real (type ki), mass^2 of propagator 2
!  * m2 -- a real (type ki), mass^2 of propagator 1
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
!  * form_factor_2p (src/form_factor/form_factor_2p.f90)
!  * form_factor_higher_ranks (src/form_factor/form_factor_higher_ranks.f90)
!  * constante, only : b_null (src/module/constante.f90)
!  * sortie_erreur, only : tab_erreur_par,catch_exception (src/module/sortie_erreur.f90)
!
!
! EXAMPLE
!
!
!
!*****
function gb0(s1,m1,m2,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_2p 
  use form_factor_higher_ranks
  use constante
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  real(ki), intent (in) :: s1,m1,m2,mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gb0
  !
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
   !   
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
  call initgolem95(2)
  !
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m1
      !
  call preparesmatrix()
  !
  ff = a20(b_null)
      !
  lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
   ! gb0 = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
   ! mu2_scale_par contained already in basic function
   ! use local mu2 !
    gb0 = ff%c
    !
  else if (eps_flag == -1) then
    !
    !gb0 = ff%b + lmu2*ff%a 
    gb0 = ff%b
    !
  else if (eps_flag == -2) then
    !
    gb0 = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gb0 (gb0.f90)'
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
end function gb0
! 
!****f* src/interface/gb0c
! NAME
!
!  Function gb0c
!
! USAGE
!
!  complex = gb0c(s1,m1,m2,mu2,eps_flag)
!
! DESCRIPTION
!
!  This function corresponds to the scalar b0 function.
!  There are two arguments more which are the renormalisation 
!  scale squared and a flag which selects the coefficient of 
!  the Laurent series in epsilon
!
! INPUTS
!
!  * s1 -- a complex (type ki), p1^2
!  * m1 -- a complex (type ki), mass^2 of propagator 2
!  * m2 -- a complex (type ki), mass^2 of propagator 1
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
!  * form_factor_2p (src/form_factor/form_factor_2p.f90)
!  * form_factor_higher_ranks (src/form_factor/form_factor_higher_ranks.f90)
!  * constante, only : b_null (src/module/constante.f90)
!  * sortie_erreur, only : tab_erreur_par,catch_exception (src/module/sortie_erreur.f90)
!
!
! EXAMPLE
!
!
!
!*****
function gb0c(s1,m1,m2,mu2,eps_flag)
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor
  use form_factor_2p 
  use form_factor_higher_ranks
  use constante
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: mu2_scale_par
  implicit none
  !
  complex(ki), intent (in) :: s1,m1,m2
  real(ki), intent(in) :: mu2
  integer, intent(in) :: eps_flag
  complex(ki) :: gb0c
  !
  type(form_factor) :: ff
  real(ki) :: lmu2,mu2store
   !   
  mu2store=mu2_scale_par
  mu2_scale_par=mu2
 !
  call initgolem95(2)
  !
  !
  s_mat(1,1) = -2._ki*m2
  s_mat(1,2) = s1-m1-m2
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2._ki*m1
      !
  call preparesmatrix()
  !
  ff = a20(b_null)
      !
  lmu2 = log(mu2)
  !
  if (eps_flag == 0) then
    !
   ! gb0c = ff%c + lmu2*ff%b + lmu2**2*ff%a/2._ki
   ! mu2_scale_par contained already in basic function
   ! use local mu2 !
    gb0c = ff%c
    !
  else if (eps_flag == -1) then
    !
    !gb0c = ff%b + lmu2*ff%a 
    gb0c = ff%b
    !
  else if (eps_flag == -2) then
    !
    gb0c = ff%a
    !
  else
    !
    tab_erreur_par(1)%a_imprimer = .true.
    tab_erreur_par(1)%chaine = 'In function gb0c (gb0.f90)'
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
end function gb0c
