!~ changed 13.5.2010 to include scale (mu^2)^eps
!~ the default scale is 1, defined in parametre.f90 
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
  public :: f2p_m1m2, i2sm1m2, i2sm1, i2sm1m2_old
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
  !  real_dim4 = i2sm1m2(s,m1sq_r,m2sq_r)
  !  complex_dim2 = i2sm1m2(s,m1sq_r,m2sq_r)
  !
  ! DESCRIPTION
  !
  !  This function computes the scalar two point function
  !  with non-zero momentum and two nonzero masses: I_2(s,m1^2,m2^2)
  !  in n dimensions
  !
  ! INPUTS
  !
  !  * msq1,m2sq -- real/complex (type ki), the masses squared
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
  function i2sm1m2_old(s,m1,m2)
    !
    real(ki), intent(in) :: s,m1,m2
    real(ki), dimension(4) :: i2sm1m2_old
    real(ki) :: delta,sig
    complex(ki) :: i2fin,tlog1,tlog2,x1,x2,t1,t2,r1,r2,rlog1,rlog2
    !
    delta = s**2+m2**2+m1**2-2._ki*s*m2-2._ki*s*m1-2._ki*m2*m1    
    sig = sign(1._ki,s)
    !
    if (delta >= 0._ki) then
       !
       x1 = (s+m1-m2+sqrt(delta))/(2._ki*s)  ! Im x1 =sign(s)*i*eps 
       x2 = (s+m1-m2-sqrt(delta))/(2._ki*s)   ! Im x2 =-sign(s)*i*eps 
       !
    else !if (delta < 0._ki) then
       !
       x1 = (s+m1-m2+(-sig*i_)*sqrt(-delta))/(2._ki*s)
       x2 = (s+m1-m2-(-sig*i_)*sqrt(-delta))/(2._ki*s)
       !
    end if
    !
    t1=(x1-1._ki)
    t2=(x2-1._ki)
    r1=(x1-1._ki)/x1
    r2=(x2-1._ki)/x2
    !
    if ( .not.(equal_real(aimag(t1),zero))) then
       tlog1= log(t1)
    else  
       tlog1=z_log(real(t1),1._ki*sig)
    endif
    !
    if ( .not.(equal_real(aimag(t2),zero))) then
       tlog2= log(t2)
    else  
       tlog2=z_log(real(t2),-1._ki*sig)
    endif
    !
    if ( .not.(equal_real(aimag(x1),zero))) then
       rlog1= log(x1)
    else  
       rlog1=z_log(real(x1),1._ki*sig)
    endif
    !
    if ( .not.(equal_real(aimag(x2),zero))) then
       rlog2= log(x2)
    else  
       rlog2=z_log(real(x2),-1._ki*sig)
    endif
    !
    !
    ! ***** to be checked ! **************************
    !
    i2fin = 2._ki-z_log(s/mu2_scale_par,1._ki) + t1*tlog1-x1*rlog1 + t2*tlog2-x2*rlog2
    !
    ! *************************************************
    !
    i2sm1m2_old(1) = 1._ki
    i2sm1m2_old(2) = 0._ki
    !
    if (rat_or_tot_par%tot_selected) then
       !
       i2sm1m2_old(3) = real(i2fin)
       i2sm1m2_old(4) = aimag(i2fin)
       !
    else !if (rat_or_tot_par%rat_selected) then
       !
       i2sm1m2_old(3) = 2._ki
       i2sm1m2_old(4) = 0._ki
       !
    end if
    !
  end function i2sm1m2_old
  !
  function i2sm1m2_r(s,m1,m2)
    !
    real(ki), intent(in) :: s,m1,m2
    real(ki), dimension(4) :: i2sm1m2_r
    real(ki) :: delta, smm
    complex(ki) :: xlog1, xlog2, x1, x2, sqrtd, lm1, lm2
    complex(ki) :: i2fin
    !
    delta = s**2+m2**2+m1**2-2._ki*s*m2-2._ki*s*m1-2._ki*m2*m1    
    !    
    smm = -s+m1+m2
    lm1 = cmplx(log(m1/mu2_scale_par),0._ki,ki)
    lm2 = cmplx(log(m2/mu2_scale_par),0._ki,ki)
    !
    if (equal_real(delta,zero) ) then 
       !
       sqrtd = czero
       xlog1 = czero !!! this is set to zero to allow a faster evaluation
       xlog2 = czero
       !
    else if (delta .gt. 0._ki) then
       !
       sqrtd = cmplx(sqrt(delta),0._ki,ki)
       x1 = cmplx(smm,0._ki,ki) + sqrtd   ! Im r1 = +i*eps 
       x2 = cmplx(smm,0._ki,ki) - sqrtd   ! Im r2 = -i*eps 
       xlog1 = z_log(x1,1._ki)
       xlog2 = z_log(x2,-1._ki)
       !
    else !if (delta .lt. 0._ki) then
       !
       sqrtd = cmplx(0._ki,sqrt(-delta),ki)
       x1 = cmplx(smm,0._ki,ki) + sqrtd
       x2 = cmplx(smm,0._ki,ki) - sqrtd
       xlog1 = z_log(x1,1._ki)
       xlog2 = z_log(x2,-1._ki)
       !
    end if
    !
    ! ***** to be checked ! **************************
    !
    i2fin = 2._ki+( (smm - cmplx(2._ki*m1,0._ki,ki) )*lm1 + &
         &          (smm - cmplx(2._ki*m2,0._ki,ki) )*lm2 + &
         &                               sqrtd*(xlog1-xlog2) )/2._ki/s
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
    complex(ki) , intent(in) :: m1, m2
    complex(ki), dimension(2) :: i2sm1m2_c
    real(ki) :: sig
    complex(ki) :: delta, smm, sc
    complex(ki) :: i2fin, xlog1, xlog2, x1, x2, sqrtd, lm1, lm2
    !
    sc = cmplx(s,zero,ki)
    !
    smm = -sc+m1+m2
    ! delta = sc**2+m2**2+m1**2-2._ki*sc*m2-2._ki*sc*m1-2._ki*m2*m1
    delta = smm*smm - 4._ki * m1*m2
    !
    lm1 = z_log(m1/mu2_scale_par,-1._ki)
    lm2 = z_log(m2/mu2_scale_par,-1._ki)
    !
    if ( (equal_real(aimag(delta),zero) .and. (real(delta,ki) .lt. zero) ) ) then 
       !
       sig = sign(un,s)
       sqrtd = (sig*i_)*sqrt(-delta)
       !
    else
       !
       sqrtd = sqrt(delta)
       !
    end if
    !
    x1 = (smm + sqrtd)/m1 
    x2 = (smm - sqrtd)/m1
    !
    !
    xlog1 = z_log(x1, 1._ki) !!! If im-part vanishes, the real part is always positive! (for real s).
    xlog2 = z_log(x2,-1._ki) !!! eps-prescription checked for maximal complex s,
                             !!! (si=-(Sqrt[m1gam1]+Sqrt[m2gam2])^2!
    !                        !!! in case of real masses, the right branch is taken!!
    !
    ! *** This needs to be checked ***
    !
    i2fin = 2._ki + ( (smm-2._ki*m1)*lm1 + &
         &  (smm-2._ki*m2)*lm2 + sqrtd*(xlog1-xlog2) )/2._ki/sc
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
  !
end module function_2p_m1m2
