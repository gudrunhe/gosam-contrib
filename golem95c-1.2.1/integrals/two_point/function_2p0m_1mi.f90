!~ changed 13.5.2010 to include scale (mu^2)^eps
!~ the default scale is 1, defined in parametre.f90 
!****h* src/integral/two_point/function_2p0m
! NAME
!
!  Module function_2p0m
!
! USAGE
!
!  use function_2p0m
!
! DESCRIPTION
!
!  This module is used to compute the two-point function
!  with zero momentum and two equal masses: I_2(0,m^2,m^2)
!  and the two-point function
!  with zero momentum and two different masses: I_2(0,m1^2,m2^2)
!  with/without Feynman parameters in n dimensions
!
! OUTPUT
!
!  This module exports the functions:
!  * f2p0m_1mi -- a function for the computation of the 
!  two-point integrals
!  with zero momentum and two equal masses: I2({j})(0,m^2,m^2)
!  with/without Feynman parameters, in n dimensions
!
!  * f2p0m_m1m2 -- a function for the computation of the 
!  two-point integrals
!  with zero momentum and two different masses: I2({j})(0,m1^2,m2^2)
!  with/without Feynman parameters, in n dimensions
!
!  scalar functions:
!
!  i20m1: computes the scalar two point function
!  with zero momentum and one propagator having nonzero mass: 
!  I_2^n(0,0,m^2)
!
!  i20mm: computes the scalar two point function
!  with zero momentum and two massive propagators 
!  with equal masses: I_2^n(0,m^2,m^2)
!
!  i20m1m2: computes the scalar two point function
!  with zero momentum and two massive propagators 
!  with different masses: I_2^n(0,m1^2,m2^2)
!
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * logarithme (src/module/z_log.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * constante (src/module/constante.f90)
!
!*****
module function_2p0m_1mi
  !
  use precision_golem
  use logarithme
  use sortie_erreur
  use equal
  use parametre
  use constante, only : zero, czero
  !
  implicit none
  !
  private  
  !
  interface i20m1
     !
     module procedure i20m1_r, i20m1_c
     !
  end interface
  !
  interface i20mm
     !
     module procedure i20mm_r, i20mm_c
     !
  end interface
  !
  interface f2p0m_1mi
     !
     module procedure f2p0m_1mi_r, f2p0m_1mi_c
     !
  end interface
  !
  interface f2p0m_m1m2
     !
     module procedure f2p0m_m1m2_r, f2p0m_m1m2_c
     !
  end interface
  !
  interface i20m1m2
     !
     module procedure i20m1m2_r, i20m1m2_c
     !
  end interface
  !
  public :: f2p0m_1mi, f2p0m_m1m2, i20m1, i20mm,i20m1m2
  !
contains
  !
  !
  !****f* src/integral/two_point/function_2p0m/f2p0m_1mi
  ! NAME
  !
  !  Function f2p0m_1mi
  !
  ! USAGE
  !
  !  real_dim4 = f2p0m_1mi(msq_r,par1,par2)
  !  complex_dim2 = f2p0m_1mi(msq_c,par1,par2)
  !
  ! DESCRIPTION
  !
  !  This function computes the two point function in n dimensions
  !  with zero momentum and two massive propagators with m1=m2
  !  with up to two Feynman parameters in the numerator.
  !  It retuns an array of (4 reals / 2 complex) corresponding to the
  !   real/imaginary part of the coefficient of the 1/epsilon term 
  !   and the real/imaginary part of the constant term.
  ! 
  ! corresponds to eqs.(A.9),(A.10) in hep-ph/0504267
  ! note overall minus sign has to be corrected in first line of (A.10)
  ! note also that for rank one A_j^{2,1}=MINUS I_2(j,...)
  !
  ! INPUTS
  !
  !  * m1_sq -- real/complex (type ki), the value of the mass
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
  !  the real/imaginary part of the coefficient of the coefficient 
  !  of the 1/epsilon term
  !  and the real/imaginary part of the constant term.
  !
  ! EXAMPLE
  !
  ! light-like-momentum two point function without Feynman parameters 
  ! f2p0m_1mi(msq,0,0) 
  ! with one Feynman parameter in the numerator z_1 
  ! f2p0m_1mi(msq,0,1)
  ! with two Feynman parameters in the numerator z_2**2
  ! f2p0m_1mi(msq,2,2) 
  ! with two Feynman parameters in the numerator z1*z_2
  ! f2p0m_1mi(msq,1,2) 
  !
  !***** 
  !
  function f2p0m_1mi_r(m1_sq,par1,par2)
    !
    real(ki), intent (in) :: m1_sq
    integer, intent (in) :: par1,par2
    real(ki), dimension(4) :: f2p0m_1mi_r
    !
    f2p0m_1mi_r(:) = 0._ki
    !
    ! scalar case
    if ( (par1 == 0) .and. (par2 == 0) ) then
       !
       f2p0m_1mi_r = i20mm(m1_sq)
       !
       ! rank one
    else if ( (par1 == 0) .and. (par2 == 1) ) then
       !
       f2p0m_1mi_r =  i20mm(m1_sq)/2._ki
       !
    else if ( (par1 == 0) .and. (par2 == 2) ) then
       !
       f2p0m_1mi_r =  i20mm(m1_sq)/2._ki
       !
       ! rank two
    else if ( (par1 == 1) .and. (par2 == 1) ) then
       !
       f2p0m_1mi_r = i20mm(m1_sq)/3._ki
       !
    else if ( (par1 == 1) .and. (par2 == 2) ) then
       !
       f2p0m_1mi_r = i20mm(m1_sq)/6._ki
       !
    else if ( (par1 == 2) .and. (par2 == 2) ) then
       !
       f2p0m_1mi_r = i20mm(m1_sq)/3._ki
       !
    end if
    !
    !
  end function f2p0m_1mi_r
  !
  function f2p0m_1mi_c(m1_sq,par1,par2)
    !
    complex(ki), intent (in) :: m1_sq
    integer, intent (in) :: par1,par2
    complex(ki), dimension(2) :: f2p0m_1mi_c
    !
    f2p0m_1mi_c(:) = czero
    !
    ! scalar case
    if ( (par1 == 0) .and. (par2 == 0) ) then
       !
       f2p0m_1mi_c = i20mm(m1_sq)
       !
       ! rank one
    else if ( (par1 == 0) .and. (par2 == 1) ) then
       !
       f2p0m_1mi_c =  i20mm(m1_sq)/2._ki
       !
    else if ( (par1 == 0) .and. (par2 == 2) ) then 
       !
       f2p0m_1mi_c =  i20mm(m1_sq)/2._ki
       !
       ! rank two
    else if ( (par1 == 1) .and. (par2 == 1) ) then
       !
       f2p0m_1mi_c = i20mm(m1_sq)/3._ki
       !
    else if ( (par1 == 1) .and. (par2 == 2) ) then 
       !
       f2p0m_1mi_c = i20mm(m1_sq)/6._ki
       !
    else if ( (par1 == 2) .and. (par2 == 2) ) then 
       !
       f2p0m_1mi_c = i20mm(m1_sq)/3._ki
       !
    end if
    !
    !
  end function f2p0m_1mi_c
  !
  !****f* src/integral/two_point/function_2p0m/f2p0m_m1m2
  ! NAME
  !
  !  Function f2p0m_m1m2
  !
  ! USAGE
  !
  !  real_dim6 = f2p0m_m1m2(m1sq_r,m2sq_r,par1,par2)
  !  complex_dim3 = f2p0m_m1m2(m1sq_c,m2sq_c,par1,par2)
  !
  ! DESCRIPTION
  !
  !  This function computes the two point function in n dimensions
  !  with zero momentum and two massive propagators with m1 not= m2
  !  with up to two Feynman parameters in the numerator.
  !  It retuns an array of (6 reals / 3 complex) corresponding to the real/imaginary
  !  part of the coefficient of the 1/epsilon**2 term, real/imaginary part of the
  !  coefficient of the 1/epsilon term and the real/imaginary part of the 
  !  constant term.
  ! corresponds to eqs.(A.8) in hep-ph/0504267
  !  note that for rank one A_j^{2,1}=MINUS I_2(j,...)
  !
  ! INPUTS
  !
  !  * m1_sq,m2_sq -- real/complex (type ki), the values of the masses
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
  !  An real/complex (type ki) array of rank 1 and shape 6/3 corresponding to 
  !  the real/imaginary part of the coefficient of the 1/epsilon**2 term,
  !  real/imaginary part of the coefficient of the 1/epsilon term
  !  and the real/imaginary part of the constant term.
  !
  ! EXAMPLE
  !
  ! light-like-momentum two point function without Feynman parameters 
  ! f2p0m_m1m2(m1sq,m2sq,0,0) 
  ! with one Feynman parameter in the numerator z_1 
  ! f2p0m_m1m2(m1sq,m2sq,0,1)
  ! with two Feynman parameters in the numerator z_2**2
  ! f2p0m_m1m2(m1sq,m2sq,2,2) 
  ! with two Feynman parameters in the numerator z1*z_2
  ! f2p0m_m1m2(m1sq,m2sq,1,2) 
  !
  !***** 
  function f2p0m_m1m2_r(m1_sq,m2_sq,par1,par2)
    !
    real(ki), intent (in) :: m1_sq,m2_sq
    integer, intent (in) :: par1,par2
    real(ki), dimension(4) :: f2p0m_m1m2_r
    real(ki) :: small,diffrm
    !
    f2p0m_m1m2_r(:) = 0._ki
    diffrm=sqrt(m1_sq)-sqrt(m2_sq)
    small=1.e-6_ki
    !
    ! scalar case
    if ( (par1 == 0) .and. (par2 == 0) ) then
       !
       f2p0m_m1m2_r = i20m1m2(m1_sq,m2_sq)
       !
       ! rank one, z1
    else if ( (par1 == 0) .and. (par2 == 1) ) then
       !
       f2p0m_m1m2_r(1) = 1._ki/2._ki
       f2p0m_m1m2_r(2) = 0._ki
       !
       !
       if (abs(diffrm) > small ) then 
          ! 
          if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_r(3) = -(-m1_sq**2 + 4._ki*m1_sq*m2_sq - 3._ki*m2_sq**2 +    &
                  &                             2._ki*m1_sq*(m1_sq - 2._ki*m2_sq)*real(z_log(m1_sq/mu2_scale_par,-1._ki)) +        &
                  &                             2._ki*m2_sq**2*real(z_log(m2_sq/mu2_scale_par,-1._ki)))/(4._ki*(m1_sq - m2_sq)**2)
             f2p0m_m1m2_r(4) =-(2._ki*m1_sq*(m1_sq - 2._ki*m2_sq)*aimag(z_log(m1_sq/mu2_scale_par,-1._ki)) +  &
                  &                           2._ki*m2_sq**2*aimag(z_log(m2_sq/mu2_scale_par,-1._ki)))/(4._ki*(m1_sq - m2_sq)**2)
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq - 3._ki*m2_sq)/(4._ki*(m1_sq - m2_sq))
             f2p0m_m1m2_r(4) = 0._ki
             !
          end if  ! end if rat or tot
          !
       else  ! use expansion in (m2sq-m1sq) up to order 3
          !    
          if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq) * ( 27._ki*m1_sq**2 - 9._ki*m1_sq*m2_sq + 2._ki*m2_sq**2 )/(120._ki*m1_sq**3) - &
                  &                   real(z_log(m1_sq/mu2_scale_par,-1._ki))/2._ki
             !
             f2p0m_m1m2_r(4) = - aimag(z_log(m1_sq/mu2_scale_par,-1._ki))/2._ki
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq) * ( 27._ki*m1_sq**2 - 9._ki*m1_sq*m2_sq + 2._ki*m2_sq**2 )/(120._ki*m1_sq**3)
             f2p0m_m1m2_r(4) = 0._ki
             !
          end if  ! end if rat or tot
          !    
       end if ! end if abs(diffrm) > small
       !	   
       ! rank one, z2
    else if ( (par1 == 0) .and. (par2 == 2) ) then
       !
       f2p0m_m1m2_r(1) = 1._ki/2._ki
       f2p0m_m1m2_r(2) = 0._ki
       !
       if (abs(diffrm) > small ) then 
          ! 
          if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_r(3) = -(-m2_sq**2 + 4._ki*m2_sq*m1_sq - 3._ki*m1_sq**2 +    &
                  &          2._ki*m2_sq*(m2_sq - 2._ki*m1_sq)*real(z_log(m2_sq/mu2_scale_par,-1._ki)) +        &
                  &          2._ki*m1_sq**2*real(z_log(m1_sq/mu2_scale_par,-1._ki)))/(4._ki*(m2_sq - m1_sq)**2)
             f2p0m_m1m2_r(4) =-(2._ki*m2_sq*(m2_sq - 2._ki*m1_sq)*aimag(z_log(m2_sq/mu2_scale_par,-1._ki)) +  &
                  &          2._ki*m1_sq**2*aimag(z_log(m1_sq/mu2_scale_par,-1._ki)))/(4._ki*(m2_sq - m1_sq)**2)
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_r(3) = (m2_sq - 3._ki*m1_sq)/(4._ki*(m2_sq - m1_sq))
             f2p0m_m1m2_r(4) = 0._ki
             !
          end if  ! end if rat or tot
          !
       else  ! use expansion in (m2sq-m1sq) up to order 3
          !    
          if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq)*(63._ki*m1_sq**2 - 31._ki*m1_sq*m2_sq + 8._ki*m2_sq**2 )/(120._ki*m1_sq**3) - &
                  &                  real(z_log(m1_sq/mu2_scale_par,-1._ki),ki)/2._ki
             !
             f2p0m_m1m2_r(4) = - aimag(z_log(m1_sq/mu2_scale_par,-1._ki))/2._ki
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq)*(63._ki*m1_sq**2 - 31._ki*m1_sq*m2_sq + 8._ki*m2_sq**2 )/(120._ki*m1_sq**3)
             f2p0m_m1m2_r(4) = 0._ki
             !
          end if  ! end if rat or tot
          !    
       end if ! end if abs(diffrm) > small
       !	   
       ! rank two
    else if ( (par1 == 1) .and. (par2 == 1) ) then
       !
       f2p0m_m1m2_r(1) = 1._ki/3._ki
       f2p0m_m1m2_r(2) = 0._ki
       !
       if (rat_or_tot_par%tot_selected) then
          !
          if (abs(diffrm) > small ) then 
             !   	write(6,*) 'using unexpanded expression, diffrm=',diffrm
             f2p0m_m1m2_r(3) = 1._ki/9._ki-m2_sq/(m1_sq-m2_sq)/6._ki +  &
                  &                          m2_sq**2/(m1_sq-m2_sq)**2/3._ki -  &
                  &                          real(z_log(m1_sq/mu2_scale_par,-1._ki))/3._ki  -   &
                  &    m2_sq**3*(real(z_log(m1_sq/mu2_scale_par,-1._ki))- &
                  &    real(z_log(m2_sq/mu2_scale_par,-1._ki)))/(m1_sq-m2_sq)**3/3._ki
             !
             f2p0m_m1m2_r(4) = -  aimag(z_log(m1_sq/mu2_scale_par,-1._ki))/3._ki  -   &
                  &  m2_sq**3*(aimag(z_log(m1_sq/mu2_scale_par,-1._ki))- &
                  &  aimag(z_log(m2_sq/mu2_scale_par,-1._ki)))/(m1_sq-m2_sq)**3/3._ki
             !
          else  ! use expansion in (m2sq-m1sq) up to order 3
             !
             f2p0m_m1m2_r(3) =  (m1_sq-m2_sq)*( 19._ki*m1_sq**2 - 5._ki*m1_sq*m2_sq + m2_sq**2 )/(180._ki*m1_sq**3) -&
                  &                         real(z_log(m1_sq/mu2_scale_par,-1._ki))/3._ki   
             !
	     ! mu2 dependence corrected 18.7.2012 GH
             f2p0m_m1m2_r(4) = - aimag(z_log(m1_sq/mu2_scale_par,-1._ki))/3._ki
             !
          end if ! end if abs(diffrm) > small
          !
       else if (rat_or_tot_par%rat_selected) then
          !
    	  if (abs(diffrm) > small ) then 
             f2p0m_m1m2_r(3) = 1._ki/9._ki-m2_sq/(m1_sq-m2_sq)/6._ki +  &
                  &                          m2_sq**2/(m1_sq-m2_sq)**2/3._ki
    	  else  ! use expansion in (m2sq-m1sq) up to order 3
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq)*( 19._ki*m1_sq**2 - 5._ki*m1_sq*m2_sq + m2_sq**2 )/(180._ki*m1_sq**3)
          end if ! end if abs(diffrm) > small
          f2p0m_m1m2_r(4) = 0._ki
          !
       end if   ! end if rat or tot
       !
       !
    else if ( (par1 == 1) .and. (par2 == 2) ) then
       !
       f2p0m_m1m2_r(1) = 1._ki/6._ki
       f2p0m_m1m2_r(2) = 0._ki
       !
       if (abs(diffrm) > small ) then 
          !
          if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_r(3) = (5._ki*m1_sq**3 - 27._ki*m1_sq**2*m2_sq + 27._ki*m1_sq*m2_sq**2 -  &
                  &             5._ki*m2_sq**3 - 6._ki*m1_sq**2*(m1_sq - 3._ki*m2_sq)*real(z_log(m1_sq/mu2_scale_par,-1._ki),ki) + &
                  &             6._ki*m2_sq**2*(-3._ki*m1_sq + m2_sq)*real(z_log(m2_sq/mu2_scale_par,-1._ki),ki))/  &
                  &             (36._ki*(m1_sq - m2_sq)**3)
             f2p0m_m1m2_r(4) =( - 6._ki*m1_sq**2*(m1_sq - 3._ki*m2_sq)*aimag(z_log(m1_sq/mu2_scale_par,-1._ki)) +   &
                  &                             6._ki*m2_sq**2*(-3._ki*m1_sq + m2_sq)*aimag(z_log(m2_sq/mu2_scale_par,-1._ki)))/  &
                  &                            (36._ki*(m1_sq - m2_sq)**3)
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_r(3) = (5._ki*m1_sq**2 - 22._ki*m1_sq*m2_sq + 5._ki*m2_sq**2)/(36._ki*(m1_sq - m2_sq)**2)
             f2p0m_m1m2_r(4) = 0._ki
             !
          end if ! end if rat or tot
          !
       else  !  use expansion in (m2sq-m1sq) up to order 3
          !
          if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq)*(43._ki*m1_sq**2 - 17._ki*m1_sq*m2_sq + 4._ki*m2_sq**2 )/(360._ki*m1_sq**3) - &
                  &                  real(z_log(m1_sq/mu2_scale_par,-1._ki))/6._ki
             !
             f2p0m_m1m2_r(4) = - aimag(z_log(m1_sq/mu2_scale_par,-1._ki))/6._ki
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq)*(43._ki*m1_sq**2 - 17._ki*m1_sq*m2_sq + 4._ki*m2_sq**2 )/(360._ki*m1_sq**3)
             f2p0m_m1m2_r(4) = 0._ki
             !
          end if ! end if rat or tot
          !
       end if ! end if abs(diffrm) > small
       !	
    else if ( (par1 == 2) .and. (par2 == 2) ) then
       !
       f2p0m_m1m2_r(1) = 1._ki/3._ki
       f2p0m_m1m2_r(2) = 0._ki
       !
       if (abs(diffrm) > small ) then 
          !
          if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_r(3) = (11._ki*m1_sq**3 - 18._ki*m1_sq**2*m2_sq +  &
                  &                             9._ki*m1_sq*m2_sq**2 - 2._ki*m2_sq**3 -     &
                  &                             6._ki*m1_sq**3*real(z_log(m1_sq/mu2_scale_par,-1._ki)) +  &
                  &			       6._ki*m2_sq*(3._ki*m1_sq**2 - 3._ki*m1_sq*m2_sq + &
                  &                                m2_sq**2)*real(z_log(m2_sq/mu2_scale_par,-1._ki)))/ &
                  &                                (18._ki*(m1_sq - m2_sq)**3)
             f2p0m_m1m2_r(4) =( -6._ki*m1_sq**3*aimag(z_log(m1_sq/mu2_scale_par,-1._ki)) +  &
                  &                            6._ki*m2_sq*(3._ki*m1_sq**2 - 3._ki*m1_sq*m2_sq + &
                  &                               m2_sq**2)*aimag(z_log(m2_sq/mu2_scale_par,-1._ki)))/  &
                  &                               (18._ki*(m1_sq - m2_sq)**3)
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_r(3) = (11._ki*m1_sq**2 - 7._ki*m1_sq*m2_sq + 2._ki*m2_sq**2 )/(18._ki*(m1_sq - m2_sq)**2)
             f2p0m_m1m2_r(4) = 0._ki
             !
          end if ! end if rat or tot
          !	
       else  !  use expansion in (m2_sq-m1sq) up to order 3
          !
          if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq)*(73._ki*m1_sq**2 - 38._ki*m1_sq*m2_sq + 10._ki*m2_sq**2 )/(180._ki*m1_sq**3) - &
                  &                   real(z_log(m1_sq/mu2_scale_par,-1._ki))/3._ki
             !
             f2p0m_m1m2_r(4) = - aimag(z_log(m1_sq/mu2_scale_par,-1._ki))/3._ki
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_r(3) = (m1_sq-m2_sq)*(73._ki*m1_sq**2 - 38._ki*m1_sq*m2_sq + 10._ki*m2_sq**2 )/(180._ki*m1_sq**3)
             f2p0m_m1m2_r(4) = 0._ki
             !
          end if ! end if rat or tot
          !
       end if ! end if abs(diffrm) > small
       !	
    end if ! end test values of par1,par2
    !
    !
  end function f2p0m_m1m2_r
  !
  !
  function f2p0m_m1m2_c(m1_sq,m2_sq,par1,par2)
    !
    complex(ki), intent (in) :: m1_sq,m2_sq
    integer, intent (in) :: par1,par2
    complex(ki), dimension(2) :: f2p0m_m1m2_c
    complex(ki) :: ratpart
    real(ki) :: small,diffrm
    !
    f2p0m_m1m2_c(:) = czero
    diffrm = sqrt(abs(m1_sq-m2_sq))
    small = 1.e-6_ki
    !
    ! scalar case
    if ( (par1 == 0) .and. (par2 == 0) ) then
       !
       f2p0m_m1m2_c = i20m1m2(m1_sq,m2_sq)
       !
       ! rank one, z1
    else if ( (par1 == 0) .and. (par2 == 1) ) then
       !
       f2p0m_m1m2_c(1) = 1._ki/2._ki
       !
       !
       if (diffrm > small ) then 
          ! 
          ratpart =  (m1_sq - 3._ki*m2_sq)/((m1_sq - m2_sq)*4._ki)
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) =  ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) =  ratpart - ( m1_sq*(m1_sq - 2._ki*m2_sq)*z_log(m1_sq/mu2_scale_par,-1._ki) +        &
                  &                         m2_sq**2*z_log(m2_sq/mu2_scale_par,-1._ki) )/(2._ki*(m1_sq - m2_sq)**2)
             !
          end if  ! end if rat or tot
          !
       else  ! use expansion in (m2sq-m1sq) up to order 3
          !
          ratpart = (m1_sq-m2_sq)*(27._ki*m1_sq**2 - 9._ki*m1_sq*m2_sq + 2._ki*m2_sq**2)/(120._ki*m1_sq**3)
          !
          if (rat_or_tot_par%rat_selected) then 
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - z_log(m1_sq/mu2_scale_par,-1._ki)/2._ki
             !
          end if  ! end if rat or tot
          !    
       end if ! end if abs(diffrm) > small
       !	   
       ! rank one, z2
    else if ( (par1 == 0) .and. (par2 == 2) ) then
       !
       f2p0m_m1m2_c(1) = 1._ki/2._ki
       !
       if (diffrm > small ) then 
          ! 
          ratpart = (3._ki*m1_sq - m2_sq)/(4._ki*(m1_sq-m2_sq))
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - ( m2_sq*(m2_sq - 2._ki*m1_sq)*z_log(m2_sq/mu2_scale_par,-1._ki)  +     &
                  &                m1_sq**2*z_log(m1_sq/mu2_scale_par,-1._ki) )/(2._ki*(m2_sq - m1_sq)**2)
             !
          end if  ! end if rat or tot
          !
       else  ! use expansion in (m2sq-m1sq) up to order 3
          !
          ratpart = (m1_sq - m2_sq)*(63._ki*m1_sq**2 - 31._ki*m1_sq*m2_sq + 8._ki*m2_sq**2)/(120._ki*m1_sq**3) 
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - z_log(m1_sq/mu2_scale_par,-1._ki)/2._ki
             !
          end if  ! end if rat or tot
          !    
       end if ! end if abs(diffrm) > small
       !	   
       ! rank two
    else if ( (par1 == 1) .and. (par2 == 1) ) then
       !
       f2p0m_m1m2_c(1) = 1._ki/3._ki
       !        
       if (diffrm > small ) then 
          !
          ratpart = (2._ki*m1_sq**2 - 7._ki*m1_sq*m2_sq + 11._ki*m2_sq**2)/(18._ki*(m1_sq - m2_sq)**2)
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - ( m1_sq*(m1_sq**2-3._ki*m1_sq*m2_sq+3._ki*m2_sq**2)*z_log(m1_sq/mu2_scale_par,-1._ki) - &
                  &          m2_sq**3*z_log(m2_sq/mu2_scale_par,-1._ki) )/(3._ki*(m1_sq - m2_sq)**3)
             !         
          end if  ! end if rat or tot
          !
       else  ! use expansion in (m2sq-m1sq) up to order 3
          !
          write(6,*) 'using expansion in m1s-m2s, diffrm=',diffrm
          !
          ratpart = (m1_sq - m2_sq)*(19._ki*m1_sq**2 - 5._ki*m1_sq*m2_sq + m2_sq**2)/(180._ki*m1_sq**3)
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - z_log(m1_sq/mu2_scale_par,-1._ki)/3._ki
             !
          end if ! end if rat or tot
          !
       end if   ! end if  abs(diffrm) > small
       !
       !
    else if ( (par1 == 1) .and. (par2 == 2) ) then
       !
       f2p0m_m1m2_c(1) = 1._ki/6._ki
       !
       if (diffrm > small ) then 
          !
          ratpart = (5._ki*m1_sq**2 - 22._ki*m1_sq*m2_sq + 5._ki*m2_sq**2)/(36._ki*(m1_sq - m2_sq)**2)
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - ( m1_sq**2*(m1_sq - 3._ki*m2_sq)*z_log(m1_sq/mu2_scale_par,-1._ki) - &
                  &          m2_sq**2*(m2_sq - 3._ki*m1_sq)*z_log(m2_sq/mu2_scale_par,-1._ki) )/(6._ki*(m1_sq - m2_sq)**3)
             !
          end if ! end if rat or tot
          !
       else  !  use expansion in (m2sq-m1sq) up to order 3
          !
          ratpart = (m1_sq - m2_sq)*(43._ki*m1_sq**2 - 17._ki*m1_sq*m2_sq + 4._ki*m2_sq**2)/(360._ki*m1_sq**3)
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - z_log(m1_sq/mu2_scale_par,-1._ki)/6._ki
             !
          end if ! end if rat or tot
          !
       end if ! end if abs(diffrm) > small
       !
    else if ( (par1 == 2) .and. (par2 == 2) ) then
       !
       f2p0m_m1m2_c(1) = 1._ki/3._ki
       !
       if (diffrm > small ) then 
          !
          ratpart = (11._ki*m1_sq**2 - 7._ki*m1_sq*m2_sq + 2._ki*m2_sq**2)/(18._ki*(m1_sq - m2_sq)**2)
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - ( m1_sq**3*z_log(m1_sq/mu2_scale_par,-1._ki) - &
                  &           m2_sq*(3._ki*m1_sq**2-3._ki*m1_sq*m2_sq+m2_sq**2)*z_log(m2_sq/mu2_scale_par,-1._ki) )/ &
                  &           (3._ki*(m1_sq-m2_sq)**3)
             !
          end if ! end if rat or tot
          !	
       else  !  use expansion in (m2_sq-m1sq) up to order 3
          !
          ratpart = (m1_sq - m2_sq)*(73._ki*m1_sq**2 - 38._ki*m1_sq*m2_sq + 10._ki*m2_sq**2)/(180._ki*m1_sq**3)
          !
          if (rat_or_tot_par%rat_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p0m_m1m2_c(2) = ratpart - z_log(m1_sq/mu2_scale_par,-1._ki)/3._ki
             !
          end if ! end if rat or tot
          !
       end if ! end if abs(diffrm) > small
       !	
    end if ! end test values of par1,par2
    !
    !
  end function f2p0m_m1m2_c
  !
  !
  !
  ! ************* scalar functions *****************
  !
  ! ****f* src/integral/two_point/i20m1
  ! NAME
  !
  !  Function i20m1
  !
  ! USAGE
  !
  !  real_dim4 = i20m1(msq_r)
  !  complex_dim2 = i20m1(msq_c)
  !
  ! DESCRIPTION
  !
  !  This function computes the scalar two point function
  !  with zero momentum and one nonzero mass: I_2(0,0,m**2)
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
  function i20m1_r(msq)
    !
    real(ki), intent(in) :: msq
    real(ki), dimension(4) :: i20m1_r
    !
    i20m1_r(1) = 1._ki
    i20m1_r(2) = 0._ki
    !
    if (rat_or_tot_par%tot_selected) then
       !
       i20m1_r(3) = 1._ki-real(z_log(msq/mu2_scale_par,-1._ki))
       i20m1_r(4) = -aimag(z_log(msq/mu2_scale_par,-1._ki))
       !
    else if (rat_or_tot_par%rat_selected) then
       !
       i20m1_r(3) = 1._ki
       i20m1_r(4) = 0._ki
       !
    end if
    !
  end function i20m1_r
  !
  function i20m1_c(msq)
    !
    complex(ki), intent(in) :: msq
    complex(ki), dimension(2) :: i20m1_c
    !
    i20m1_c(1) = cmplx(1._ki,0._ki,ki)
    !
    if (rat_or_tot_par%rat_selected) then
       !
       i20m1_c(2) = cmplx(1._ki,0._ki,ki)
       !
    else if (rat_or_tot_par%tot_selected) then
       !
       i20m1_c(2) = 1._ki - z_log(msq/mu2_scale_par,-1._ki)
       !
    end if
    !
  end function i20m1_c
  !
  ! ****f* src/integral/two_point/i20mm
  ! NAME
  !
  !  Function i20mm
  !
  ! USAGE
  !
  !  real_dim4 = i20mm(msq_r)
  !  complex_dim2 = i20mm(msq_c)
  !
  ! DESCRIPTION
  !
  !  This function computes the scalar two point function
  !  with zero momentum and two equal nonzero masses: I_2(0,m^2,m^2)
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
  function i20mm_r(msq)
    !
    real(ki), intent(in) :: msq
    real(ki), dimension(4) :: i20mm_r
    !
    i20mm_r(1) = 1._ki
    i20mm_r(2) = 0._ki
    !
    if (rat_or_tot_par%tot_selected) then
       !
       i20mm_r(3) = -real(z_log(msq/mu2_scale_par,-1._ki))
       i20mm_r(4) = -aimag(z_log(msq/mu2_scale_par,-1._ki))
       !
    else if (rat_or_tot_par%rat_selected) then
       !
       i20mm_r(3) = 0._ki
       i20mm_r(4) = 0._ki
       !
    end if
    !
  end function i20mm_r
  !
  function i20mm_c(msq)
    !
    complex(ki), intent(in) :: msq
    complex(ki), dimension(2) :: i20mm_c
    !
    i20mm_c(1) = cmplx(1._ki,0._ki,ki)
    !
    if (rat_or_tot_par%rat_selected) then
       !
       i20mm_c(2) = czero
       !
    else if (rat_or_tot_par%tot_selected) then
       !
       ! scale dependence corrected 18.7.2012 GH
       !
       i20mm_c(2) = -z_log(msq/mu2_scale_par,-1._ki)
       !
    end if
    !
  end function i20mm_c
  !
  ! ****f* src/integral/two_point/i20m1m2
  ! NAME
  !
  !  Function i20m1m2
  !
  ! USAGE
  !
  !  real_dim4 = i20m1m2(msq1_r,msq2_r)
  !  complex_dim2 = i20m1m2(msq1_c,msq2_c)
  !
  ! DESCRIPTION
  !
  !  This function computes the scalar two point function
  !  with zero momentum and two equal nonzero masses: I_2(0,m1**2,m2**2)
  !  in n dimensions
  !
  ! INPUTS
  !
  !  * m1sq,m2sq -- real/complex (type ki), the masses squared
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
  function i20m1m2_r(m1sq,m2sq)
    !
    real(ki), intent(in) :: m1sq,m2sq
    real(ki), dimension(4) :: i20m1m2_r
    !
    if (equal_real(m1sq,m2sq)) then
       !
       i20m1m2_r = i20mm(m1sq)
       !   
    else
       ! 
       if (rat_or_tot_par%tot_selected) then
          !
          i20m1m2_r = (m2sq*i20m1(m2sq)-m1sq*i20m1(m1sq))/(m2sq-m1sq)
          !
       else if (rat_or_tot_par%rat_selected) then
          !
          i20m1m2_r(1) = 1._ki
          i20m1m2_r(2) = 0._ki
          i20m1m2_r(3) = 1._ki
          i20m1m2_r(4) = 0._ki
          !
       end if  ! end if rat or tot
       !
    end if  ! end test if m1=m2
    ! 
  end function i20m1m2_r
  !
  function i20m1m2_c(m1sq,m2sq)
    !
    complex(ki), intent(in) :: m1sq,m2sq
    complex(ki), dimension(2) :: i20m1m2_c
    real(ki) :: diffm
    !
    diffm = abs(m1sq-m2sq)
    !
    if (equal_real(diffm,zero)) then
       !
       i20m1m2_c = i20mm(m1sq)
       !   
    else
       ! 
       if (rat_or_tot_par%rat_selected) then
          !
          i20m1m2_c(1) = cmplx(1._ki,0._ki,ki)
          i20m1m2_c(2) = cmplx(1._ki,0._ki,ki)
          !
       else if (rat_or_tot_par%tot_selected) then
          !
          i20m1m2_c = (m2sq*i20m1(m2sq)-m1sq*i20m1(m1sq))/(m2sq-m1sq)
          !
       end if  ! end if rat or tot
       !
    end if  ! end test if m1=m2
    ! 
  end function i20m1m2_c
  !
  !
  !
end module function_2p0m_1mi
