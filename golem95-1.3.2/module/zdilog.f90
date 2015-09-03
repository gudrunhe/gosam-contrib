!****h* src/module/dilogarithme
! NAME
!
!  Module dilogarithme
!
! USAGE
!
!  use dilogarithme
!
! DESCRIPTION
!
!  This module provides two public routines to compute the dilogarithm with
!  real and complex argument
!
! OUTPUT
!
!  It exports:
!  * zdilog -- a function which returns the dilogarithm with real (or complex) argument
!  * cdilog -- a function which returns the dilogarithm with complex argument
!
! NOTES
! 
!  zdilog can be called with complex argument and parameter s. If the imaginary part vanishes,
!  zdilog with real argument is called. Here, s becomes important.
!
! USES
!
!  * precision_golem (src/module/precision.f90)
!  * constante (src/module/constante.f90)
!  * logarithme (src/module/z_log.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * equal (src/module/equal.f90)
!
!*****
module dilogarithme
  !
  use precision_golem
  use constante, only : un,pi,pi6,pi12,zero,czero,cun
  use logarithme, only : z_log,z_log2
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only: if_print_warn_par
  use equal
  !
  implicit none
  !
  private
  real(ki), parameter :: zeta2 = 1.6449340668482264364724151666460252_ki
  ! bern_glob contient B_k(0)/(k+1)!
  real(ki), dimension (0:20), parameter :: bern_glob = (/ &
    1.00000000000000000000000000000000_ki,&
    -0.25000000000000000000000000000000_ki,&
    0.02777777777777777777777777777778_ki,&
    -0.00027777777777777777777777777778_ki,&
    0.00000472411186696900982615268330_ki,&
    -0.00000009185773074661963550852440_ki,&
    0.00000000189788699889709990720092_ki,&
    -0.00000000004064761645144225526806_ki,&
    0.00000000000089216910204564525552_ki,&
    -0.00000000000001993929586072107569_ki,&
    0.00000000000000045189800296199182_ki,&
    -0.00000000000000001035651761218125_ki,&
    0.00000000000000000023952186210262_ki,&
    -0.00000000000000000000558178587433_ki,&
    0.00000000000000000000013091507554_ki,&
    -0.00000000000000000000000308741980_ki,&
    0.00000000000000000000000007315976_ki,&
    -0.00000000000000000000000000174085_ki,&
    0.00000000000000000000000000004158_ki,&
    -0.00000000000000000000000000000100_ki,&
    0.00000000000000000000000000000002_ki&
  /)
  integer :: imax_glob = ki+ki/4 ! determine how many elements of the array 
                                 ! bern_glob have to be taken as a function of ki
  
  interface zdilog
     !
     module procedure zdilog_r
     module procedure zdilog_c
     !
  end interface
  !
  public :: zdilog, cdilog
  !
  contains
    !
    !****f* src/module/dilogarithme/zdilog
    ! NAME
    !
    !  Function zdilog
    !
    ! USAGE
    !
    !  complex = zdilog(a,s)
    !
    ! DESCRIPTION
    !
    !  This function returns the dilogarithm of a complex z, this complex number 
    !  has the specific form: z = a + i lambda s where lambda << 1.
    !  a can now be complex. If the imaginary part vanishes, the sign of s is relevant.
    !
    ! INPUTS
    !
    !  * a -- a real/complex (type ki), the argument
    !  * s -- a real (type ki), s = +/- 1, it gives the sign of the small imaginary part
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function zdilog_r(a,s)
      !
      real(ki), intent(in) :: a,s
      complex(ki) :: zdilog_r
      !
      if (abs(s) == un) then
        !
        if (a <= un) then
          !
          zdilog_r = dilog(a)
          !
        else
          !
          zdilog_r = -dilog(1._ki/a)-pi6-0.5_ki*z_log2(-a,-s)
          !
        endif
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'error in zdilog :'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the second argument must be 1. or -1. %f0'
        tab_erreur_par(2)%arg_real = s
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      endif
      !
    end function zdilog_r
    !
    function zdilog_c(a,s)
      !
      complex(ki), intent(in) :: a
      real(ki), intent(in) :: s
      complex(ki) :: zdilog_c
      !
      if (equal_real(aimag(a),zero) ) then
         !
         zdilog_c = zdilog_r(real(a,ki),s)
         !
      else
         !
         zdilog_c = cdilog(a)
         !
      end if
      !
    end function zdilog_c
    !****if* src/module/dilogarithme/dilog
    ! NAME
    !
    !  Function dilog
    !
    ! USAGE
    !
    !  real = dilog(x)
    !
    ! DESCRIPTION
    !
    !  This function return the dilogarithm of a real x for x < 1
    !
    ! INPUTS
    !
    !  * x -- a real (type ki), the argument of the dilogarithm
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a real (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function dilog(x)
      !
      real(ki), intent (in) :: x
      real(ki) :: dilog
      !
      real(ki), parameter :: un_demi = 0.5_ki
      real(ki) :: arg, temp_ln, add_on, s, temp, ln_arg
      integer :: i
      !
      if ( equal_real(x,zero) ) then
        !
        dilog = zero
        !
      else if ( equal_real(x,un) ) then
        !
        dilog = pi6
        !
      else if ( equal_real(x,-un) ) then
        !
        dilog = -pi12
        !
      else if ( equal_real(x,un_demi) ) then
        !
        dilog = pi12 - 0.5_ki*log(2._ki)**2
        !
      else
        !
        if (x < -un) then
          !
          arg = 1._ki/(1._ki-x)
          s = un
          temp_ln = log(1._ki-x)
          add_on = -pi6 + temp_ln*( 0.5_ki*temp_ln - log(-x) )
          !
        else if ( x < zero) then
          !
          arg = x/(x - 1._ki)
          s = -un
          add_on = -0.5_ki*log(1._ki-x)*log(1._ki-x)
          !
        else if ( x < un_demi) then
          !
          arg = x
          s = un
          add_on = 0._ki
          !
        else if ( x < un ) then
          !
          arg = 1._ki-x
          s = -un
          add_on = pi6 - log(x)*log(1._ki-x)
          !
        else
          !

          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'error in dilog :'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'dilog(x) called for x >= 1: x=%f0'
          tab_erreur_par(2)%arg_real = x
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
        end if
        !
        ln_arg = -log(1._ki-arg)
        temp = 1._ki
        dilog = bern_glob(0) + bern_glob(1)*ln_arg
        !
        do i = 2,imax_glob
          !
          temp = temp*ln_arg*ln_arg
          dilog = dilog + bern_glob(i)*temp
          !
        end do
        !
        dilog = dilog*ln_arg
        dilog = s*dilog + add_on
        !
      end if
      !
    end function dilog
    !
    !****f* src/module/dilogarithme/cdilog
    ! NAME
    !
    !  Function cdilog
    !
    ! USAGE
    !
    !  complex = cdilog(z)
    !
    ! DESCRIPTION
    !
    !  This function return the dilogarithm of a complex z, taken from T. Binoth
    !
    ! INPUTS
    !
    !  * z -- a complex (type ki), the argument of the dilogarithm
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    ! adapted from Binoth's program
    ! modified 10.5.2010: set Im part to zero if < delta to avoid 
    ! errors at z=(1._ki,delta)
    function cdilog(z)
      !
      complex(ki), intent (in) :: z
      complex(ki) :: cdilog
      real(ki) :: delta
      logical :: old_if_print_warn_par
      !
      delta=10._ki*epsilon(1._ki)
    !  
      if ( (real(z)<1._ki+delta).and.(real(z)>1._ki-delta).and.(abs(aimag(z))< 0.000000045 ) ) then 
              cdilog = zeta2
      !
      else if (z == czero) then
        !
        cdilog = czero
        !
      else if (z == cun) then
        !
        cdilog = zeta2
        !
      else if(abs(z) <= 0.5_ki) then
        !
        cdilog = cdilog6(z)
        !
      else if(real(z) < 0._ki) then
        !
       cdilog = cdilog2(z) 
        !
        ! > 0 changed to >=0 Feb 14, 2011
      else if(real(z) >= 0._ki) then
        !
        cdilog = cdilog3(z) 
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'error in function cdilog :'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the argument z is not in the good range :  %z0'
        tab_erreur_par(2)%arg_comp = z
        old_if_print_warn_par=if_print_warn_par
        if_print_warn_par=.true.
        call catch_exception(1)
        if_print_warn_par=old_if_print_warn_par
        !
        ! to please the compiler
        !stop
        cdilog = z
        !
      end if
      !
    end function cdilog
    !
    !****if* src/module/dilogarithme/cdilog2
    ! NAME
    !
    !  Function cdilog2
    !
    ! USAGE
    !
    !  complex = cdilog2(z)
    !
    ! DESCRIPTION
    !
    ! Transform the dilog function with Re(z) < 0 to a dilog with Re(z) >= 0
    !
    ! INPUTS
    !
    !  * z -- a complex (type ki), the argument of the dilogarithm
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function cdilog2(z)
      !
      complex(ki), intent(in) :: z
      complex(ki) :: cdilog2
      !
      if(real(z) >= 0._ki) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'spence function error in cdilog2 at z =  %z0'
        tab_erreur_par(1)%arg_comp = z
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      else !if(real(z) < 0._ki) then
        !
        cdilog2 = - cdilog3(1._ki - z) &
                 - (log(z))*(log(1._ki - z)) + zeta2
        !
      endif
      !
    end function cdilog2
    !
    !****if* src/module/dilogarithme/cdilog3
    ! NAME
    !
    !  Function cdilog3
    !
    ! USAGE
    !
    !  complex = cdilog3(z)
    !
    ! DESCRIPTION
    !
    ! Transform the dilog function with Re(z) >= 0 to a dilog with Re(z)>=0 & |z|<=1
    !
    ! INPUTS
    !
    !  * z -- a complex (type ki), the argument of the dilogarithm
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function cdilog3(z)
      !
      complex(ki), intent(in) :: z
      complex(ki) :: cdilog3
      !
      if(real(z) < 0._ki) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'spence function error in cdilog3 at z =  %z0'
        tab_erreur_par(1)%arg_comp = z
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      else if(abs(z) <= 1._ki) then
        !
        cdilog3 = cdilog4(z)
        !
      else ! if(abs(z) > 1._ki) then
        !
        cdilog3 = - cdilog4(1._ki /z) & 
                 - 1._ki /2._ki * (log(-z))**2 - zeta2
        !
      endif
      !
    end function cdilog3
    !     
    !****if* src/module/dilogarithme/cdilog4
    ! NAME
    !
    !  Function cdilog4
    !
    ! USAGE
    !
    !  complex = cdilog4(z)
    !
    ! DESCRIPTION
    !
    ! Separate the case |z| < 1/2 and |z| >= 1/2
    !
    ! INPUTS
    !
    !  * z -- a complex (type ki), the argument of the dilogarithm
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function cdilog4(z)
      !
      complex(ki), intent(in) :: z
      complex(ki) :: cdilog4
      complex(ki) :: z1, z2, z3, z4
      !
      if(real(z) < 0._ki) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'spence function error 1 in cdilog4 at z = %z0'
        tab_erreur_par(1)%arg_comp = z
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      else if(abs(z) > 1._ki) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'spence function error 2 in cdilog4 at z = %z0'
        tab_erreur_par(1)%arg_comp = z
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      else if(abs(z) <= 0.5_ki) then
        !
        cdilog4 = cdilog6(z) 
        !
      else ! if(abs(z) > 0.5_ki) then
        !
        z1 = sqrt(sqrt(z))
        z2 = sqrt(1._ki + sqrt(z))
        z3 = sqrt(z)
        z4 = 1._ki + sqrt(z)
        !
        cdilog4 = 2*( &
            2*(   cdilog5(z1) + cdilog5(1._ki /(1._ki + z1)) &
             + cdilog5(1._ki /z2) + cdilog5(1._ki /(1._ki + 1._ki/z2)) &
             - (log(1._ki + z1))*(log(z1)) &
             - (log(1._ki + 1._ki /z2))*(log(1._ki /z2)) &
             + 1._ki /2._ki * (log(1._ki + z1))**2 &
             + 1._ki /2._ki * (log(1._ki + 1._ki /z2))**2 &
             - 2*zeta2 &
            ) &
             - (log(z4))*(log(z3)) &
             + 1._ki /2._ki * (log(z4))**2 &
             - zeta2 &
                  )
      endif
      !
    end function cdilog4
    !
    !****if* src/module/dilogarithme/cdilog5
    ! NAME
    !
    !  Function cdilog5
    !
    ! USAGE
    !
    !  complex = cdilog5(z)
    !
    ! DESCRIPTION
    !
    ! compute the case |z| >= 1/2 and arg(z) <= Pi/8
    !
    ! INPUTS
    !
    !  * z -- a complex (type ki), the argument of the dilogarithm
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function cdilog5(z)
      !
      complex(ki), intent(in) :: z
      complex(ki) :: cdilog5
      !
      if(abs(z) > (1._ki + 1.e-3_ki)) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'spence function error in cdilog5 at z = %z0'
        tab_erreur_par(1)%arg_comp = z
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      else if(abs(z) <= 0.5_ki) then
        !
        ! this had been a place giving an error, as we didn't expect it ever to
        ! be entered.  However, for z near to (1,0), but just outside the range
        ! we give (very rare), we can still enter here. Neverthelesss, the
        ! function is well-behaved here so this shouldn't cause any trouble.
         cdilog5 = cdilog6(z) 

       !
      else if(abs(aimag(log(z))) > pi/8._ki .and. .not.equal_real(abs(aimag(log(z))),pi/8._ki) ) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'spence function error in cdilog5 at z = %z0'
        tab_erreur_par(1)%arg_comp = z
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      else
        !
        cdilog5 = - cdilog6(1._ki - z) & 
                 - log(z)*log(1._ki - z) + zeta2
        !
      end if
      !
    end function cdilog5
    !
    !****if* src/module/dilogarithme/cdilog6
    ! NAME
    !
    !  Function cdilog6
    !
    ! USAGE
    !
    !  complex = cdilog6(z)
    !
    ! DESCRIPTION
    !
    ! compute the case |z| <= 1/2
    !
    ! INPUTS
    !
    !  * z -- a complex (type ki), the argument of the dilogarithm
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function cdilog6(z)
      !
      complex(ki), intent (in) :: z
      complex(ki) :: cdilog6
      !
      complex(ki) :: lnz,temp
      integer :: i
      !
      lnz = -log(1._ki-z)
      temp = 1._ki
      cdilog6 = bern_glob(0)+bern_glob(1)*lnz
      !
      do i = 2,imax_glob
        !
        temp = temp*lnz*lnz
        cdilog6 = cdilog6 + bern_glob(i)*temp
        !
      end do
      !
      cdilog6 = cdilog6*lnz
      !
    end function cdilog6
  !
end module dilogarithme
