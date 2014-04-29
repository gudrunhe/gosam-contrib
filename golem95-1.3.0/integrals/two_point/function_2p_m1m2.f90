! changed 13.5.2010 to include scale (mu^2)^eps
! the default scale is 1, defined in parametre.f90 
!****h* src/integral/two_point/function_2p_m1m2
! NAME
!
!  Module function_2p_m1m2
!
! USAGE
!
!  use function_2p_m1m2
!
! DESCRIPTION
!
!  This module is used to compute the two-point function
!  I_2(s,m1^2,m2^2)
!  with/without Feynman parameters in n dimensions
!
! OUTPUT
!
!  This module exports the functions:
!  * f2p_m1m2 -- a function for the computation of 
!  two-point integrals
!  with non-zero momentum and two masses: I2^n({zj})(s,m1^2,m2^2)
!  with/without Feynman parameters, in n dimensions
!  one of the masses can be zero
!  massless case is already contained in generic_function_2p
!
!  
!  i2sm1m2: computes the scalar two point function
!  where both propagators have nonzero mass: 
!  I_2^n(s,m1^2,m2^2)
!
!  i2sm1: computes the scalar two point function
!  where only one propagator has nonzero mass: 
!  I_2^n(s,m^2,0)
!
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * logarithme (src/module/z_log.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * function_2p0m_1mi (src/integrals/two_point/function_2p0m_1mi.f90)
!
!*****
!
module function_2p_m1m2
  !
  use precision_golem
  use logarithme
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use constante
  use equal
  use function_2p0m_1mi
  use parametre, only : rat_or_tot_par,mu2_scale_par
  !
  implicit none
  !
  private 
  !
  interface f2p_m1m2
     !
     module procedure f2p_m1m2_r, f2p_m1m2_c
     !
  end interface
  !
  interface i2sm1m2
     !
     module procedure i2sm1m2_r, i2sm1m2_c
     !
  end interface
  !
  interface i2sm1
     !
     module procedure i2sm1_r, i2sm1_c
     !
  end interface
  !
  real(ki) :: grand_glob = huge(1._ki)
  !
  public :: f2p_m1m2, i2sm1m2, i2sm1
  !
contains
  !
  !
  !****f* src/integral/two_point/function_2p_m1m2/f2p_m1m2
  ! NAME
  !
  !  Function f2p_m1m2
  !
  ! USAGE
  !
  !  real_dim4 = f2p_m1m2(s,msq1_r,msq2_r,par1,par2)
  !  complex_dim2 =  f2p_m1m2(s,msq1_c,msq2_c,par1,par2)
  !
  ! DESCRIPTION
  !
  !  This function computes the 
  !  two point function in n dimensions
  !  with non-zero momentum and two massive propagators
  !  with up to two Feynman parameters in the numerator.
  !  It retuns an array of 4 reals / 2 complex corresponding to the real/imaginary
  !  part of the coefficient of the 
  !  1/epsilon term and the real/imaginary part of the 
  !  constant term.
  !  corresponds to eqs.(A.5),(A.7) in hep-ph/0504267
  !  note that for rank one A_j^{2,1}=MINUS I_2(j,...)
  !
  ! INPUTS
  !
  !  * m1,m2 -- real/complex (type ki), the value of the masses
  !  * par1 -- an integer, the label of one Feynman parameter
  !  * par2 -- an integer, the label of the second Feynman parameter
  !  Note that par1,par2 are ordered internally, i.e.
  !  par1 <= par2, note also to use zero for par1, par2 
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
  !  An real/complex (type ki) array of rank 1 and shape 4/2 corresponding to 
  !  the real/imaginary part of the coefficient of the 1/epsilon term
  !  and the real/imaginary part of the constant term.
  !
  ! EXAMPLE
  !
  ! light-like-momentum two point function without Feynman parameters 
  ! f2p_m1m2(s,m1sq,m2sq,0,0) 
  ! with one Feynman parameter in the numerator z_1 
  ! f2p_m1m2(s,m1sq,m2sq,0,1)
  ! with two Feynman parameters in the numerator z_2^2
  ! f2p_m1m2(s,m1sq,m2sq,2,2) 
  ! with two Feynman parameters in the numerator z1*z_2
  ! f2p_m1m2(s,m1sq,m2sq,1,2) 
  !
  !***** 
  function f2p_m1m2_r(s,m1,m2,par1,par2)
    ! m1 and m2 are the squared masses
    ! should only be called if s, m1, m2 are nonzero
    !
    real(ki), intent (in) :: s,m1,m2
    integer, intent (in) :: par1,par2
    real(ki), dimension(4) :: f2p_m1m2_r, i2sca
    !
    f2p_m1m2_r(:) = 0._ki
    i2sca =  i2sm1m2(s,m1,m2)
    !
    ! scalar case
    if ( (par1 == 0) .and. (par2 == 0) ) then
       ! 
       f2p_m1m2_r = i2sca
       !
       ! rank one: note that rat or tot is distinguished in i2sm1m2
    else if ( (par1 == 0) .and. (par2 == 1) ) then
       !
       f2p_m1m2_r =  i2sca/2._ki -   &
            &                 (m1-m2)/s/2._ki*( i2sca - i20m1m2(m1,m2) )
       !
    else if ( (par1 == 0) .and. (par2 == 2) ) then
       !
       f2p_m1m2_r =  i2sca/2._ki + &
            &                  (m1-m2)/s/2._ki*( i2sca - i20m1m2(m1,m2) )
       !
       ! rank two: rat singled out explicitly here
    else if ( (par1 == 1) .and. (par2 == 1) ) then
       !
       f2p_m1m2_r(1) = 1._ki/3._ki
       f2p_m1m2_r(2) = 0._ki
       !	
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_m1m2_r(3) = (-6._ki*m1**2 + 12._ki*m1*m2 - 6._ki*m2**2 + 9._ki*m1*s - 9._ki*m2*s + s**2 + &
               &                     6._ki*(m1**2 + m2**2 + m2*s + s**2 - 2._ki*m1*(m2 + s))* &
               &                     i2sca(3) +                                            &
               &                    6._ki*m1*(m1 - m2 - 2._ki*s)*real(z_log(m1/mu2_scale_par,-1._ki)) - &
               &                    6._ki*m2*(m1 - m2 - s)*real(z_log(m2/mu2_scale_par,-1._ki)))/(18._ki*s**2)
          !
          f2p_m1m2_r(4) = ( 6._ki*(m1**2 + m2**2 + m2*s + s**2 - 2._ki*m1*(m2 + s))*   &
               &                   i2sca(4) +                                  &
               &                   6._ki*m1*(m1 - m2 - 2*s)*aimag(z_log(m1/mu2_scale_par,-1._ki)) -           &
               &                   6._ki*m2*(m1 - m2 - s)*aimag(z_log(m2/mu2_scale_par,-1._ki)))/(18._ki*s**2)
          !
       else !if (rat_or_tot_par%rat_selected) then
          !
          f2p_m1m2_r(3) = (-6._ki*m1**2 + 12._ki*m1*m2 - 6._ki*m2**2 + 9._ki*m1*s - 9._ki*m2*s + s**2 + &
               &                     6._ki*(m1**2 + m2**2 + m2*s + s**2 - 2._ki*m1*(m2 + s))*     &
               &                     i2sca(3) )/(18._ki*s**2)
          !
          f2p_m1m2_r(4) = ( 6._ki*(m1**2 + m2**2 + m2*s + s**2 - 2._ki*m1*(m2 + s))*   &
               &                   i2sca(4) )/(18._ki*s**2)
          !
       end if ! end if rat or tot
       !   
    else if ( (par1 == 1) .and. (par2 == 2) ) then
       !
       f2p_m1m2_r(1) = 1._ki/6._ki
       f2p_m1m2_r(2) = 0._ki
       !	
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_m1m2_r(3) = (6._ki*m1**2 - 12._ki*m1*m2 + 6._ki*m2**2 - s**2 +   &
               &                    3._ki*(-2._ki*m1**2 - 2._ki*m2**2 + m2*s + s**2 + m1*(4._ki*m2 + s))*  &
               &                    i2sca(3) +      &
               &                    3._ki*m1*(-2._ki*m1 + 2._ki*m2 + s)*real(z_log(m1/mu2_scale_par,-1._ki)) +  &
               &                    3._ki*m2*( 2._ki*m1 - 2._ki*m2 + s)*&
               &                    real(z_log(m2/mu2_scale_par,-1._ki)))/(18._ki*s**2)
          !
          f2p_m1m2_r(4) = (3._ki*(-2._ki*m1**2 - 2._ki*m2**2 + m2*s + s**2 + m1*(4._ki*m2 + s))*  &
               &                    i2sca(4) +      &
               &                    3._ki*m1*(-2._ki*m1 + 2._ki*m2 + s)*aimag(z_log(m1/mu2_scale_par,-1._ki)) +  &
               &                    3._ki*m2*( 2._ki*m1 - 2._ki*m2 + s)*aimag(z_log(m2/mu2_scale_par,-1._ki)))/ &
               &                    (18._ki*s**2)
          !
       else !if (rat_or_tot_par%rat_selected) then
          !
          f2p_m1m2_r(3) = (6._ki*m1**2 - 12._ki*m1*m2 + 6._ki*m2**2 - s**2 +   &
               &                    3._ki*(-2._ki*m1**2 - 2._ki*m2**2 + m2*s + s**2 + m1*(4._ki*m2 + s))*  &
               &                    i2sca(3) )/(18._ki*s**2)
          !
          f2p_m1m2_r(4) = (3._ki*(-2._ki*m1**2 - 2._ki*m2**2 + m2*s + s**2 + m1*(4._ki*m2 + s))*  &
               &                    i2sca(4) )/(18._ki*s**2)
          !
       end if ! end if rat or tot
       !   
    else if ( (par1 == 2) .and. (par2 == 2) ) then
       !
       f2p_m1m2_r(1) = 1._ki/3._ki
       f2p_m1m2_r(2) = 0._ki
       !
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_m1m2_r(3) = (-6._ki*m1**2 + 12._ki*m1*m2 - 6._ki*m2**2 - 9._ki*m1*s + 9._ki*m2*s + s**2 +  &
               &                      6._ki*(m1**2 + (m2 - s)**2 + m1*(-2._ki*m2 + s))*i2sca(3) + &
               &                      6._ki*m1*(m1 - m2 + s)*real(z_log(m1/mu2_scale_par,-1._ki)) -  &
               &                      6._ki*m2*(m1 - m2 + 2._ki*s)*real(z_log(m2/mu2_scale_par,-1._ki)))/(18._ki*s**2)
          !
          f2p_m1m2_r(4) = ( 6._ki*(m1**2 + (m2 - s)**2 + m1*(-2._ki*m2 + s))*i2sca(4) + &
               &                      6._ki*m1*(m1 - m2 + s)*aimag(z_log(m1/mu2_scale_par,-1._ki)) -  &
               &                      6._ki*m2*(m1 - m2 + 2._ki*s)*aimag(z_log(m2/mu2_scale_par,-1._ki)))/(18._ki*s**2)
          !
       else !if (rat_or_tot_par%rat_selected) then
          !
          f2p_m1m2_r(3) = ( -6._ki*m1**2 + 12._ki*m1*m2 - 6._ki*m2**2 - 9._ki*m1*s + 9._ki*m2*s + s**2 +  &
               &                      6._ki*(m1**2 + (m2 - s)**2 + m1*(-2._ki*m2 + s))*  &
               &                     i2sca(3) )/(18._ki*s**2)
          !
          f2p_m1m2_r(4) = ( 6._ki*(m1**2 + (m2 - s)**2 + m1*(-2._ki*m2 + s))*  & 
               &                      i2sca(4) )/(18._ki*s**2)
          !
       end if ! end if rat or tot
       !   
    end if
    !
    !
  end function f2p_m1m2_r
  !
  !
  function f2p_m1m2_c(s,m1,m2,par1,par2)
    ! m1 and m2 are the squared masses
    ! should only be called if s, m1, m2 are nonzero
    !
    real(ki), intent (in) :: s
    complex(ki), intent(in) :: m1,m2
    integer, intent (in) :: par1,par2
    complex(ki), dimension(2) :: f2p_m1m2_c, i2sca
    complex(ki) :: ratpart
    !
    f2p_m1m2_c(:) = czero
    i2sca=i2sm1m2(s,m1,m2)
    !
    ! scalar case
    if ( (par1 == 0) .and. (par2 == 0) ) then
       ! 
       f2p_m1m2_c = i2sca
       !
       ! rank one: note that rat or tot is distinguished in i2sm1m2  
    else if ( (par1 == 0) .and. (par2 == 1) ) then
       !
       f2p_m1m2_c =  i2sca/2._ki -   &
            &                 (m1-m2)/s/2._ki*( i2sca - i20m1m2(m1,m2) )
       !
    else if ( (par1 == 0) .and. (par2 == 2) ) then
       !
       f2p_m1m2_c =  i2sca/2._ki + &
            &                  (m1-m2)/s/2._ki*( i2sca - i20m1m2(m1,m2) )
       !
       ! rank two: rat singled out explicitly here
    else if ( (par1 == 1) .and. (par2 == 1) ) then
       !
       f2p_m1m2_c(1) = 1._ki/3._ki
       !	
       ratpart = (-6._ki*(m1-m2)**2 + 9._ki*s*(m1-m2) + s**2 +  &
            &  6._ki*(m1**2 + m2**2 + m2*s + s**2 - 2._ki*m1*(m2 + s) )*i2sca(2) )/(18._ki*s**2)
       !
       if (rat_or_tot_par%rat_selected) then
          !
          f2p_m1m2_c(2) = ratpart
          !
       else !if (rat_or_tot_par%tot_selected) then 
          !
          f2p_m1m2_c(2) = ratpart + ( m1*(m1 - m2 - 2._ki*s)*z_log(m1/mu2_scale_par,-1._ki) + &
               &          m2*(m2 - m1 + s)*z_log(m2/mu2_scale_par,-1._ki) )/(3._ki*s**2)
          !
       end if ! end if rat or tot
       !   
    else if ( (par1 == 1) .and. (par2 == 2) ) then
       !
       f2p_m1m2_c(1) = 1._ki/6._ki
       !
       ratpart = ( 6._ki*(m1-m2)**2 - s**2 + 3._ki*(-2._ki*(m1-m2)**2 + s*(m1+m2) +s**2)*i2sca(2) )/(18._ki*s**2)
       !
       if (rat_or_tot_par%rat_selected) then
          !
          f2p_m1m2_c(2) = ratpart
          !
       else !if (rat_or_tot_par%tot_selected) then
          !
          f2p_m1m2_c(2) = ratpart + (m1*(s - 2._ki*(m1-m2))*z_log(m1/mu2_scale_par,-1._ki) + &
               &          m2*(s + 2._ki*(m1-m2))*z_log(m2/mu2_scale_par,-1._ki))/(6._ki*s**2)
          !
       end if ! end if rat or tot
       !   
    else if ( (par1 == 2) .and. (par2 == 2) ) then
       !
       f2p_m1m2_c(1) = 1._ki/3._ki
       !
       ratpart = (-6._ki*(m1-m2)**2-9._ki*(m1-m2)*s+s**2+6._ki*((m1-m2)**2 + s*(m1-2*m2) + &
            &     s**2)*i2sca(2))/(18._ki*s**2)
       !
       if (rat_or_tot_par%rat_selected) then
          !
          f2p_m1m2_c(2) = ratpart
          !
       else !if (rat_or_tot_par%tot_selected) then
          !
          f2p_m1m2_c(2) = ratpart + (m1*(m1 - m2 + s)*z_log(m1/mu2_scale_par,-1._ki) + &
               &          m2*(m2 - m1 -2._ki*s)*z_log(m2/mu2_scale_par,-1._ki))/(3._ki*s**2)
          !
       end if ! end if rat or tot
       !   
    end if
    !
    !
  end function f2p_m1m2_c
  !
  ! ****f* src/integral/two_point/i2sm1m2
  ! NAME
  !
  !  Function i2sm1m2
  !
  ! USAGE
  !
  !  real_dim4 = i2sm1m2_r(s,m1,m2)
  !  complex_dim2 = i2sm1m2_c(s,m1,m2)
  !
  ! DESCRIPTION
  !
  !  This function computes the scalar two point function
  !  with non-zero momentum and two nonzero masses: I_2(s,m1^2,m2^2)
  !  in n dimensions
  !
  ! INPUTS
  !
  !  * m1,m2 -- real/complex (type ki), the masses squared
  !
  ! SIDE EFFECTS
  !
  !  No side effect, it uses the value of rat_or_tot_par 
  !  (in src/module/parametre.f90)
  !
  ! RETURN VALUE
  !
  !  It returns a real/complex (type ki) array of rank 1 and shape 4/2
  !
  !*****
  !
  ! Direct computation of I2
  ! We have to compute \int^1_0 dx \ln(s x^2+x (-s+m1-m2)+m2-i \lambda)
  ! \ln(s x^2+x (-s+m1-m2)+m2-i \lambda) = \ln(s - i \lamda) + \ln(x - x_1) + x - x_2) + \eta(-x_1,-x_2)
  ! where x_1 and x_2 are the roots of the quadratic form
  ! Note that in the case of real masses, the eta function is 0 because x_1 and x_2 are complex conjugate
  ! (this is not the case when m1 and m2 are complex!!!!)
  ! \int^1_0 dx \ln(x-x_i) = (1-x_i) \ln(1-x_i) + x_i \ln(-x_i)
  ! When a mass is zero, x_1 (x_2) goes to 1 (0) but this formula is well behaved numerically
  ! because x*log(x) --> 0 when x-->0 even in Fortran .....
  !
  function i2sm1m2_r(s,m1,m2)
    !
    real(ki), intent(in) :: s,m1,m2
    real(ki), dimension(4) :: i2sm1m2_r
    !
    real(ki) :: delta, x1_r,x2_r
    complex(ki) :: x1_c,x2_c
    complex(ki) :: i2fin,part
    real(ki) :: sm,rap1,rap2,rapm
    real(ki) :: sq_ab,delta_t
    !
    delta = s*s+m1*m1+m2*m2-2._ki*s*m1-2._ki*s*m2-2._ki*m1*m2
    sm = sign(un,m1-m2)
    !
    if (equal_real(m1,0._ki) ) then
      !
      rap1 = grand_glob*sign(un,s)
      !
    else
      !
      rap1 = s/m1
      !
    end if
    !
    if (equal_real(m2,0._ki) ) then
      !
      rap2 = grand_glob*sign(un,s)
      !
    else
      !
      rap2 = s/m2
      !
    end if
    !
    if (equal_real(s,0._ki) .and. .not.(equal_real(abs(m1-m2),0._ki)) ) then
      !
      rapm = grand_glob*sign(un,s)
      !
    else
      !
      rapm = (m1-m2)/s
      !
    end if
    !
    delta_t = abs(s)*(1._ki+rapm*rapm)-2._ki*sign(un,s)*(m1+m2)
    sq_ab = sqrt(abs(s))
    !
    !
    ! no eta function in this case for x1 and x2 are complex conjugate
    !
    if (delta >= 0._ki) then
      !
      ! careful treatment of s -->0 (only in the case delta > 0)
      ! we have to distinguish the case where |s| << m1,|m1-m2| (or |s| << m2,|m1-m2|)
      ! and the case |s|,|m1-m2| << m1,m2
      !
      if ( (abs(rap1) <= 1.e-4_ki) .or. (abs(rap2) <= 1.e-4_ki) ) then
        !
        ! Case |s|,|m1-m2| << m1,m2 
        ! In this case both x1 and x2 --> infinity when s-->0
        !
        if ( abs(rapm) <= 1.e-4_ki) then
          !
          if (sq_ab > grand_glob) then
            !
            x1_r = sign(un,s)*grand_glob
            x2_r = -sign(un,s)*grand_glob
            !
          else
            !
            x1_r = (1._ki-rapm)/2._ki + sign(un,s)*sqrt(delta_t)/2._ki/sq_ab
            x2_r = (1._ki-rapm)/2._ki - sign(un,s)*sqrt(delta_t)/2._ki/sq_ab
            !
          end if
          !
          if (sm > 0._ki) then
            !
            part = z_log(m1/mu2_scale_par,-1._ki)-q(1,1._ki/x1_r,-1._ki)-q(1,1._ki/x2_r,1._ki)
            !
          else if (sm <= 0._ki) then
            !
            part = z_log(m2/mu2_scale_par,-1._ki)-q(1,1._ki/(1._ki-x1_r),-1._ki)-q(1,1._ki/(1._ki-x2_r),1._ki)
            !
          end if
          !
        else ! case where |s| << m1,|m1-m2| (or |s| << m2,|m1-m2|)
          !
          if (sm > 0._ki) then
            !
            x1_r = -2._ki*m2/(-s+m1-m2+sqrt(delta))
            !
            if (abs((m1-m2)/s) > grand_glob) then
              !
              x2_r = -grand_glob
              !
            else
              !
              x2_r = ( -(-s+m1-m2) - sqrt(delta) )/2._ki/s ! -i_*lambda
              !
            end if
            !
            part = z_log(m1/mu2_scale_par,-1._ki)-x1_r*z_log((1._ki-x1_r)/(-x1_r),1._ki)-q(1,1._ki/x2_r,1._ki)
            !
          else if (sm <= 0._ki) then
            !
            x2_r = 2._ki*m2/(s-m1+m2+sqrt(delta))
            !
            if (abs((m1-m2)/s) > grand_glob) then
              !
              x1_r = grand_glob
              !
            else
              !
              x1_r = ( -(-s+m1-m2) + sqrt(delta) )/2._ki/s ! i_*lambda
              !
            end if
            !
            part = z_log(m2/mu2_scale_par,-1._ki)-q(1,1._ki/(1._ki-x1_r),-1._ki)+(1._ki-x2_r)*z_log((1._ki-x2_r)/(-x2_r),-1._ki)
            !
          end if
          !
        end if
        !
      else ! Case where |s| is large
        !
        x1_r = ( -(-s+m1-m2) + sqrt(delta) )/2._ki/s ! +i_*lambda
        x2_r = ( -(-s+m1-m2) - sqrt(delta) )/2._ki/s ! -i_*lambda
        !
        part = z_log(s/mu2_scale_par,-1._ki) &
          & + (1._ki-x1_r)*z_log(1._ki-x1_r,-1._ki)+x1_r*z_log(-x1_r,-1._ki) &
          & + (1._ki-x2_r)*z_log(1._ki-x2_r,1._ki)+x2_r*z_log(-x2_r,1._ki)
        !
      end if
      !
    else
      !
      !
      ! Case |s|,|m1-m2| << m1,m2 
      ! In this case both x1 and x2 --> infinity when s-->0
      !
      if ( (abs(rapm) <= 1.e-4_ki) .and. ( (abs(rap1) <= 1.e-4_ki) .or. (abs(rap2) <= 1.e-4_ki) ) ) then
        !
        if (sq_ab > grand_glob) then
          !
          x1_c = i_*grand_glob
          x2_c = -i_*grand_glob
          !
        else
          !
          x1_c = (1._ki-rapm)/2._ki + i_*sqrt(abs(delta_t))/2._ki/sq_ab
          x2_c = (1._ki-rapm)/2._ki - i_*sqrt(abs(delta_t))/2._ki/sq_ab
          !
        end if
        !
        if (sm > 0._ki) then
          !
          part = z_log(m1/mu2_scale_par,-1._ki)-q(1,1._ki/x1_c,-1._ki)-q(1,1._ki/x2_c,1._ki)
          !
        else if (sm <= 0._ki) then
          !
          part = z_log(m2/mu2_scale_par,-1._ki)-q(1,1._ki/(1._ki-x1_c),-1._ki)-q(1,1._ki/(1._ki-x2_c),1._ki)
          !
        end if
        !
      else
        ! we don't care which is x_1 and x_2, they play a symetric role in the formula
        !
        x1_c = ( -(-s+m1-m2) + i_*sqrt(abs(delta)) )/2._ki/s
        x2_c = ( -(-s+m1-m2) - i_*sqrt(abs(delta)) )/2._ki/s
        !
        part = z_log(s/mu2_scale_par,-1._ki) &
          & + (1._ki-x1_c)*log(1._ki-x1_c)+x1_c*log(-x1_c) &
          & + (1._ki-x2_c)*log(1._ki-x2_c)+x2_c*log(-x2_c)
        !
      end if
      !
    end if
    !
    i2fin = 2._ki - part
    !
    ! *************************************************
    !
    i2sm1m2_r(1) = 1._ki
    i2sm1m2_r(2) = 0._ki
    !
    if (rat_or_tot_par%rat_selected) then
       !
       i2sm1m2_r(3) = 2._ki
       i2sm1m2_r(4) = 0._ki
       !
    else !if (rat_or_tot_par%tot_selected) then
       !
       i2sm1m2_r(3) = real(i2fin,ki)
       i2sm1m2_r(4) = aimag(i2fin)
       !
    end if
    !
  end function i2sm1m2_r
  !
  function i2sm1m2_c(s,m1,m2)
    !
    real(ki), intent(in) :: s
    complex(ki), intent(in) :: m1,m2
    complex(ki), dimension(2) :: i2sm1m2_c
    !
    complex(ki) :: delta,x1,x2
    complex(ki) :: i2fin,part
    real(ki) :: sm,rap1,rap2,rapm
    real(ki) :: m1_r,m2_r,sq_ab
    complex(ki) :: delta_t
    real(ki) :: lambda
    !
    delta = s*s+m1*m1+m2*m2-2._ki*s*m1-2._ki*s*m2-2._ki*m1*m2
    !
    m1_r = abs(m1)
    m2_r = abs(m2)
    sm = sign(un,m1_r-m2_r)
    lambda = epsilon(1._ki)
    !
    if (equal_real(m1_r,0._ki) ) then
      !
      rap1 = grand_glob*sign(un,s)
      !
    else
      !
      rap1 = s/m1_r
      !
    end if
    !
    if (equal_real(m2_r,0._ki) ) then
      !
      rap2 = grand_glob*sign(un,s)
      !
    else
      !
      rap2 = s/m2_r
      !
    end if
    !
    if (equal_real(s,0._ki) .and. .not.(equal_real(abs(m1-m2),0._ki)) ) then
      !
      rapm = grand_glob*sign(un,s)*(m1-m2)/abs(m1-m2)
      !
    else
      !
      rapm = (m1-m2)/s
      !
    end if
    !
    delta_t = abs(s)*(1._ki+rapm*rapm)-2._ki*sign(un,s)*(m1+m2)
    sq_ab = sqrt(abs(s))
    !
      !
      ! careful treatment of s -->0 (only in the case delta > 0)
      ! we have to distinguish the case where |s| << m1,|m1-m2| (or |s| << m2,|m1-m2|)
      ! and the case |s|,|m1-m2| << m1,m2
      !
      if ( (abs(rap1) <= 1.e-4_ki) .or. (abs(rap2) <= 1.e-4_ki) ) then
        !
        ! Case |s|,|m1-m2| << m1,m2 
        ! In this case both x1 and x2 --> infinity when s-->0
        !
        if ( abs(rapm) <= 1.e-4_ki) then
          !
          if (sq_ab > grand_glob) then
            !
            x1 = sign(un,s)*grand_glob*sqrt(delta_t)/abs(sqrt(delta_t))
            x2 = -sign(un,s)*grand_glob*sqrt(delta_t)/abs(sqrt(delta_t))
            !
          else
            !
            x1 = (1._ki-rapm)/2._ki + sign(un,s)*sqrt(delta_t)/2._ki/sq_ab
            x2 = (1._ki-rapm)/2._ki - sign(un,s)*sqrt(delta_t)/2._ki/sq_ab
            !
          end if
          !
          if (sm > 0._ki) then
            !
            part = z_log(m1/mu2_scale_par,-1._ki)-q(1,1._ki/x1,-1._ki)-q(1,1._ki/x2,1._ki) &
              & - eta(1._ki-x1,1._ki-x2,(m1-i_*lambda)/s) + eta(-x1,-x2,(m2-i_*lambda)/s)
            !
          else if (sm <= 0._ki) then
            !
            part = z_log(m2/mu2_scale_par,-1._ki)-q(1,1._ki/(1._ki-x1),-1._ki)-q(1,1._ki/(1._ki-x2),1._ki)
            !
          end if
          !
        else ! case where |s| << m1,|m1-m2| (or |s| << m2,|m1-m2|)
          !
          if (sm > 0._ki) then
            !
            x1 = -2._ki*m2/(-s+m1-m2+sqrt(delta))
            !
            if (abs((m1-m2)/s) > grand_glob) then
              !
              x2 = -grand_glob*(-s+m1-m2)/abs(-s+m1-m2)
              !
            else
              !
              x2 = ( -(-s+m1-m2) - sqrt(delta) )/2._ki/s ! -i_*lambda
              !
            end if
            !
            part = log(m1/mu2_scale_par) - x1*log((1._ki-x1)/(-x1)) - q(1,1._ki/x2,1._ki) &
              & - eta(1._ki-x1,1._ki-x2,(m1-i_*lambda)/s) + eta(-x1,-x2,(m2-i_*lambda)/s)
            !
          else if (sm <= 0._ki) then
            !
            x2 = 2._ki*m2/(s-m1+m2+sqrt(delta))
            !
            if (abs((m1-m2)/s) > grand_glob) then
              !
              x1 = grand_glob*(-s+m1-m2)/abs(-s+m1-m2)
              !
            else
              !
              x1 = ( -(-s+m1-m2) + sqrt(delta) )/2._ki/s ! i_*lambda
              !
            end if
            !
            part = log(m2/mu2_scale_par) - q(1,1._ki/(1._ki-x1),-1._ki) + (1._ki-x2)*log((1._ki-x2)/(-x2))
            !
          end if
          !
        end if
        !
      else ! Case where |s| is large
        !
        x1 = ( -(-s+m1-m2) + sqrt(delta) )/2._ki/s ! +i_*lambda
        x2 = ( -(-s+m1-m2) - sqrt(delta) )/2._ki/s ! -i_*lambda
        !
        part = z_log(s/mu2_scale_par,-1._ki) &
          & + (1._ki-x1)*log(1._ki-x1) + x1*log(-x1) &
          & + (1._ki-x2)*log(1._ki-x2) + x2*log(-x2) &
          & + eta(-x1,-x2,(m2-i_*lambda)/s)
        !
      end if
      !
    !
    i2fin = 2._ki - part
    !
    ! *************************************************
    !
    i2sm1m2_c(1) = 1._ki
    !

    if (rat_or_tot_par%rat_selected) then
       !
       i2sm1m2_c(2) = 2._ki
       !
    else !if (rat_or_tot_par%tot_selected) then
       !
       i2sm1m2_c(2) = i2fin
       !
    end if
    !
  end function i2sm1m2_c
  !
  !
  !
  ! ****f* src/integral/two_point/i2sm1
  ! NAME
  !
  !  Function i2sm1
  !
  ! USAGE
  !
  !  real_dim4 = i2sm1(s,msq_r)
  !  complex_dim2 = i2sm1(s,msq)
  !
  ! DESCRIPTION
  !
  !  This function computes the scalar two point function
  !  with non-zero momentum and m2=0: I_2(s,m^2,0)
  !  in n dimensions
  !
  ! INPUTS
  !
  !  * msq -- a real/complex (type ki), the mass squared
  !
  ! SIDE EFFECTS
  !
  !  No side effect, it uses the value of rat_or_tot_par 
  !  (in src/module/parametre.f90)
  !
  ! RETURN VALUE
  !
  !  It returns a real/complex (type ki) array of rank 1 and shape 4/2
  !
  !*****
  function i2sm1_r(s1,m1s)
    !
    real(ki), intent(in) :: s1,m1s
    real(ki), dimension(4) :: i2sm1_r
    real(ki) :: delta
    !
    delta=1000*epsilon(1._ki)
    !
    i2sm1_r(1) = 1._ki
    i2sm1_r(2) = 0._ki
    !
    if (rat_or_tot_par%tot_selected) then
       !
       if ( abs(m1s-s1) > delta) then
          !
          i2sm1_r(3) = real( -z_log(m1s/mu2_scale_par,-1._ki)+2._ki+ &
               &              (m1s-s1)/s1*( z_log((m1s-s1)/mu2_scale_par,-1._ki)  - &
               &              z_log(m1s/mu2_scale_par,-1._ki) ) )
          i2sm1_r(4) = aimag( -z_log(m1s/mu2_scale_par,-1._ki)+2._ki+ &
               &              (m1s-s1)/s1*( z_log((m1s-s1)/mu2_scale_par,-1._ki)  - &
               &               z_log(m1s/mu2_scale_par,-1._ki) ) )
          !
       else
          ! 
          i2sm1_r(3) = real( -z_log(m1s/mu2_scale_par,-1._ki)+2._ki )
          i2sm1_r(4) = aimag( -z_log(m1s/mu2_scale_par,-1._ki)+2._ki )
          !
       end if
       !
    else !if (rat_or_tot_par%rat_selected) then
       !
       i2sm1_r(3) = 2._ki
       i2sm1_r(4) = 0._ki
       !
    end if
    !
  end function i2sm1_r
  !
  function i2sm1_c(s1,m1s)
    !
    real(ki), intent(in) :: s1
    complex(ki), intent(in) :: m1s
    complex(ki), dimension(2) :: i2sm1_c
    real(ki) :: delta
    !
    delta=1000*epsilon(1._ki)
    !
    i2sm1_c(1) = cmplx(1._ki,0._ki,ki)
    !
    if (rat_or_tot_par%rat_selected) then
       !
       i2sm1_c(2) = cmplx(2._ki,0._ki,ki)
       !
    else !if (rat_or_tot_par%tot_selected) then
       !
       if ( abs( m1s - cmplx(s1,0._ki,ki) ) > delta) then
          !
          i2sm1_c(2) = 2._ki + ( (m1s-s1)*z_log((m1s-s1)/mu2_scale_par,-1._ki)  - &
               &       m1s*z_log(m1s/mu2_scale_par,-1._ki) )/s1
          !
       else
          !
          i2sm1_c(2) = 2._ki - z_log(m1s/mu2_scale_par,-1._ki)
          !
       end if
       !
    end if
    !
  end function i2sm1_c
  !
end module function_2p_m1m2
