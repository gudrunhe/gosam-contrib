!****h* src/module/logarithme
! NAME
!
!  Module logarithme
!
! USAGE
!
!  use logarithme
!
! DESCRIPTION
!
!  This module provides three public routines to compute the logarithm,
!  the logarithm squared and a special function (generalisation of ln(1-z)/z) 
!  assuming that the argument is of the type z = a + i lambda s, where
!  lambda > 0 and << 1 and s = +/- 1. a can be a complex type. If its
!  imaginary part vanishes, the sign of s becomes relevant.
!  
!
! OUTPUT
!
!  It exports:
!  * z_log -- a function which returns the logarithm with a complex argument
!  * z_log2 -- a function which returns the logarithm squared with a complex argument
!  * q -- a recursive function which a generalisation of ln(1-z)/z
!
! NOTES
! 
!  z_log (and z_log2) can be called with complex argument and parameter s. If the imaginary part vanishes,
!  z_log (or z_log2) with real argument is called. Here, s becomes important.
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!  * constante (src/module/constante.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * equal (src/module/equal.f90)
!
!*****
module logarithme
  !
  use precision_golem
  use constante, only : pi,i_,un, zero
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use equal
  !
  implicit none
  !
  private
  !
  real(ki), parameter :: small_glob = 5.e-1_ki
  !
  interface z_log
     !
     module procedure z_log_r, z_log_c
     !
  end interface
  !
  interface z_log2
     !
     module procedure z_log2_r, z_log2_c
     !
  end interface
  !
  interface q
     !
     module procedure q_r, q_c
     !
  end interface
  !
  public :: z_log,z_log2,q,eta
  !
  contains
    !
    !****f* src/module/logarithme/z_log
    ! NAME
    !
    !  Function z_log
    !
    ! USAGE
    !
    !  complex = z_log(a,s)
    !
    ! DESCRIPTION
    !
    !  Compute the ln(z) with z = a + i lambda s
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
    ! NOTES
    !
    !  If the imaginary part of the argument vanishes, the sign of s becomes relevant.
    !
    !
    !*****
    !
    function z_log_r(a,s)
      !
      real(ki), intent(in) :: a,s
      complex(ki) :: z_log_r
      !
      if (abs(s) == 1._ki) then
        !
        if (a > 0._ki) then
          !
          z_log_r = cmplx(log(a), 0.0_ki, ki)
          !
        else
          !
          ! z_log_r = log(-a)+i_*pi*s
          z_log_r = cmplx(log(-a), pi*s, ki)
          !
        endif
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'error in z_log:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the second argument must be 1. or -1. %f0'
        tab_erreur_par(2)%arg_real = s
        call catch_exception(0)
        !
        ! to please the compiler:
        stop
        !
      endif
      !
    end function z_log_r
    !
    !
    function z_log_c(a,s)
      !
      complex(ki), intent(in) :: a
      real(ki), intent(in) :: s
      complex(ki) :: z_log_c
      !

      if (equal_real(aimag(a),zero) ) then
         !
         z_log_c = z_log_r(real(a,ki),s)
         !
      else
         !
         z_log_c = log(a)
         !
      end if
      !
    end function z_log_c
    !
    !****f* src/module/logarithme/z_log2
    ! NAME
    !
    !  Function z_log2
    !
    ! USAGE
    !
    !  complex = z_log2(a,s)
    !
    ! DESCRIPTION
    !
    !  Compute the ln(z)^2 with z = a + i lambda s
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
    ! NOTES
    !
    !  If the imaginary part of the argument vanishes, the sign of s becomes relevant.
    !
    !
    !
    !*****
    !
    function z_log2_r(a,s)
      !
      real(ki), intent(in) :: a,s
      complex(ki) :: z_log2_r

      real(ki) :: lga
      !
      if (abs(s) == 1._ki) then
        !
        if (a > 0._ki) then
          !
          lga = log(a)
          z_log2_r = cmplx(lga*lga, 0.0_ki, ki)
          !
        else
          !
          ! z_log2_r = log(-a)**2-pi**2+2._ki*pi*i_*s*log(-a)
          lga = log(-a)
          z_log2_r = cmplx((lga+pi)*(lga-pi), 2.0_ki*pi*s*lga, ki)
          !
        endif
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'error in z_log2:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the second argument must be 1. or -1. %f0'
        tab_erreur_par(2)%arg_real = s
        call catch_exception(0)
        !
        ! to please the compiler:
        stop
        !
      endif
      !
    end function z_log2_r
    !
    function z_log2_c(a,s)
      !
      complex(ki), intent(in) :: a
      real(ki), intent(in) ::s
      complex(ki) :: z_log2_c
      !
      if (equal_real(aimag(a),zero) ) then
         !
         z_log2_c = z_log2_r(real(a,ki),s)
         !
      else
         !
         z_log2_c = log(a)**2
         !
      end if
      !
    end function z_log2_c
    !
    !****f* src/module/logarithme/q
    ! NAME
    !
    !  Function q
    !
    ! USAGE
    !
    !  complex = q(n,x,s)
    !
    ! DESCRIPTION
    !
    !  It computes the function q defined recusively by
    !  q_n(X) = (q_{n-1}(X)+1/(n-1))/X 
    !  with q_1(X) = ln(1-X)/X 
    !  assuming that X = x + i*s*lambda and s=+/- 1, 
    !  Care is taken for small values of x.
    !  For x < small_glob
    !  q_n(x) = -( 1/n + \sum_{j=n+1}^\infinity x^{j-n}/j )
    !  Note that in this case there is no imaginary part.
    !
    ! INPUTS
    !
    !  * n -- an integer, the order of q
    !  * x -- a real/complex (type ki), the real part
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
    ! NOTES
    !  
    !  This function can now be called with complex x. If the imaginary part vanishes,
    !  the evaluation is switched to the version with real x and the original entry for 
    !  s becomes relevant.
    !
    !
    !
    !*****
    !
    recursive function q_r(n,x,s) result(resq)
      !
      integer, intent(in) :: n
      real(ki), intent(in) :: x,s
      complex(ki) :: resq
      !
      integer :: nm1
      real(ki) :: tt,temp,expo_x,denom
      !
      if (abs(x) > small_glob) then
        !
        if (n == 1) then
          resq = z_log(un-x,-s)/x
        else
          nm1 = n - 1
          resq = (q_r(nm1,x,s) + 1._ki/real(nm1,ki))/x
        end if
        !
      else ! no imaginary part in this case
        !
        denom = real(n,ki)
        tt = 1._ki/denom
        expo_x = un
        temp = 10._ki !artificial value to enter into the loop
        !
        do while(abs(tt-temp) >= epsilon(x))
          !
          temp = tt
          expo_x = x*expo_x
          denom = denom + un
          tt = tt + expo_x/denom
          !
        end do
        !
        resq = cmplx(-tt,0._ki,ki)
        !
      end if
      !
    end function q_r
    !
    recursive function q_c(n,x,s) result(resq)
      !
      integer, intent(in) :: n
      complex(ki), intent(in) :: x
      real(ki), intent (in) :: s
      complex(ki) :: resq
      !
      integer :: nm1
      complex(ki) :: tt,temp,denom,expo_x
      !
      if (equal_real(aimag(x), zero)) then
         !
         resq = q_r(n,real(x,ki),s)
         !
      else if (abs(x) > small_glob) then
        !
        if (n == 1) then
           !
           resq = z_log((cmplx(un,0._ki,ki) - x),-s)/x
           !
        else
           !
          nm1 = n - 1
          resq = (q_c(nm1,x,s) + 1._ki/real(nm1,ki))/x
          !
        end if
        !
      else 
        !
        denom = cmplx(real(n,ki),0._ki,ki)
        tt = 1._ki/denom
        expo_x = cmplx(un,0._ki,ki)
        temp = cmplx(10._ki,0._ki,ki) !artificial value to enter into the loop
        !
        do while(abs(tt-temp) >= epsilon(real(x,ki)))
          !
          temp = tt
          expo_x = x*expo_x
          denom = denom + cmplx(un,0._ki,ki)
          tt = tt + expo_x/denom
          !
        end do
        !
        resq = -tt
        !
      end if
      !
    end function q_c
    !
    !****f* src/module/logarithme/eta
    ! NAME
    !
    !  Function eta
    !
    ! USAGE
    !
    !  complex = eta(x,y,z)
    !
    ! DESCRIPTION
    !
    !  It computes the function eta defined by
    !  eta(x,y) - ln(x*y) - ln(x) - ln(y)
    !
    ! INPUTS
    !
    !  * x -- a complex (type ki)
    !  * y -- a complex (type ki)
    !  * z -- optional, a complex (type ki) : product x*y
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! NOTES
    !  
    !
    !
    !
    !*****
    !
    !function eta(z1,z2)
      !!
      !complex(ki), intent(in) :: z1,z2
      !complex(ki) :: eta
      !!
      !real(ki) :: im1,im2,imt
      !complex(ki) :: z1n,z2n
      !real(ki) :: r_max
      !!
      !r_max = max(abs(z1),abs(z2))
      !z1n = z1/r_max
      !z2n = z2/r_max
      !im1 = aimag(z1n)
      !im2 = aimag(z2n)
      !imt = aimag(z1n*z2n)
      !!
      !write(*,*) 'in eta :',imt
      !!write(*,*) 'in eta 1 :',z1n,z2n
      !!mettre imt = 0 si trop petit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!if (abs(imt) < epsilon(1._ki)) imt = 0._ki
      !!
      !eta = eta_mod(real(z1n,ki),im1,real(z2n,ki),im2,imt)
      !!
    !end function eta
    function eta(z1,z2,z1_z2)
      !
      complex(ki), intent(in) :: z1,z2
      complex(ki), intent(in),optional :: z1_z2
      complex(ki) :: eta
      !
      real(ki) :: im1,im2,imt
      real(ki) :: re1,re2
      !
      im1 = aimag(z1)
      im2 = aimag(z2)
      if (present(z1_z2)) then
        !
        imt = aimag(z1_z2)
        !
      else
        !
        imt = aimag(z1*z2)
        !
      end if
      re1 = real(z1,ki)
      re2 = real(z2,ki)
      !
      !write(*,*) 'in eta :',imt
      !write(*,*) 'in eta 1 :',z1n,z2n
      !mettre imt = 0 si trop petit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !if (abs(imt) < epsilon(1._ki)) imt = 0._ki
      !
      if ( (im1 >= 0._ki) .and. (im2 >= 0._ki) .and. (imt < 0._ki) ) then
        !
        eta = -2._ki*i_*pi
        !
      else if ( (im1 < 0._ki) .and. (im2 < 0._ki) .and. (imt >= 0._ki) ) then
        !
        eta = 2._ki*i_*pi
        !
      else if ( (im1 == 0._ki) .and. (im2 == 0._ki) &
                  &.and. (re1 > 0._ki) .and. (re2 > 0._ki) ) then
        !
        eta = -2._ki*i_*pi
        !
      else
        !
        eta = 0._ki
        !
      end if
      !
    end function eta
    !
    !
    !****f* src/module/logarithme/eta_mod
    ! NAME
    !
    !  Function eta_mod
    !
    ! USAGE
    !
    !  complex = eta_mod(im1,im2,imt)
    !
    ! DESCRIPTION
    !
    !  It computes the function eta defined by
    !  eta(x,y) - ln(x*y) - ln(x) - ln(y)
    !  the argument are Re(x), Im(x), Re(y), Im(y) and Im(x*y)
    !
    ! INPUTS
    !
    !  * re1 -- a real (type ki)
    !  * im2 -- a complex (type ki)
    !  * re2 -- a real (type ki)
    !  * im2 -- a complex (type ki)
    !  * imt -- a complex (type ki)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  This function returns a complex (type ki)
    !
    ! NOTES
    !  
    !
    !
    !
    !*****
    !
    !function eta_mod(re1,im1,re2,im2,imt)
      !!
      !real(ki), intent(in) :: re1,im1,re2,im2,imt
      !complex(ki) :: eta_mod
      !!
      !!write(*,*) 'eta_mod :',re1,im1,re2,im2,imt
      !if ( (im1 >= 0._ki) .and. (im2 >= 0._ki) .and. (imt < 0._ki) ) then
        !!
        !eta_mod = -2._ki*i_*pi
        !!
      !else if ( (im1 < 0._ki) .and. (im2 < 0._ki) .and. (imt >= 0._ki) ) then
        !!
        !eta_mod = 2._ki*i_*pi
        !!
      !else if ( (im1 == 0._ki) .and. (im2 == 0._ki) &
                  !&.and. (re1 > 0._ki) .and. (re2 > 0._ki) ) then
        !!
        !eta_mod = -2._ki*i_*pi
        !!
      !else
        !!
        !eta_mod = 0._ki
        !!
      !end if
      !!
    !end function eta_mod
    !
end module logarithme
