!
!****h* src/integral/three_point/function_3p1m_1mi
! NAME
!
!  Module function_3p1m_1mi
!
! USAGE
!
!  use function_3p1m_1mi
!
! DESCRIPTION
!
!  This module is used to compute the one off-shell external leg one internal mass three point function
!  with/without Feynman parameters in n, n+2 dimensions
!
! OUTPUT
!
!  This module exports two functions:
!  * f3p1m_1mi -- a function for the computation of the one off-shell external leg one internal mass three 
!    point function with/without Feynman parameters in n dimensions
!  * f3p1m_1mi_np2 -- a function for the computation of the one off-shell external leg one internal mass three 
!    point function with/without Feynman parameters in n+2 dimensions
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * logarithme (src/module/z_log.f90)
!  * dilogarithme (src/module/zdilog.f90)
!  * func_he (src/integrals/three_point/mod_he.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90) only : tab_erreur_par,catch_exception
!  * constante (src/module/constante.f90) only : un,pi6
!  * parametre (src/module/parametre.f90) only : rat_or_tot_par,mu2_scale_par
!  * array (src/module/array.f90) only : packb
!
!*****
module function_3p1m_1mi
  !
  use precision_golem
  use logarithme
  use dilogarithme
  use func_he
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use constante, only : un,pi6
  use parametre, only : rat_or_tot_par,mu2_scale_par
  use array, only : packb
  implicit none
  !
  private 
  !
  public :: f3p1m_1mi, f3p1m_1mi_np2
  !
  contains
    !
    !
    !****f* src/integral/three_point/function_3p1m_1mi/f3p1m_1mi
    ! NAME
    !
    !  Function f3p1m_1mi
    !
    ! USAGE
    !
    !  real_dim6 = f3p1m_1mi(s13,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the one off-shell external leg one internal mass three point function in n dimensions
    !  with up to three Feynman parameters in the numerator.
    !  It retuns an array of 6 reals corresponding to the real/imaginary
    !  part of the coefficient of the 1/epsilon^2 term, real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s13 -- real (type ki), the value of the S matrix element corresponding to the external off-shell leg
    !  * m3_sq -- real (type ki), the value of the internal mass squared
    !  * par1 -- an integer, the label of the third Feynman parameter
    !  * par2 -- an integer, the label of the second Feynman parameter
    !  * par3 -- an integer, the label of the first Feynman parameter
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  An real (type ki) array of rank 1 and shape 6 corresponding to 
    !  the real/imaginary part of the coefficient of the 1/epsilon^2 term,
    !  real/imaginary part of the coefficient of the 1/epsilon term
    !  and the real/imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !*****
    function f3p1m_1mi(s13,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(6) :: f3p1m_1mi
      !
      complex(ki) :: c_temp_d2,c_temp_d2_rat
      complex(ki) :: c_temp_d1,c_temp_d1_rat
      complex(ki) :: c_temp,c_temp_rat
      real(ki) :: lmu2
      !
      f3p1m_1mi = 0._ki
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        c_temp_d2=1._ki/2._ki/s13
        !
        c_temp_d2_rat=1._ki/2._ki/s13
        !
        c_temp_d1=z_log(-s13,-1._ki)/s13-1._ki/2._ki*z_log(m3_sq,-1._ki)/&
          &s13
        !
        c_temp_d1_rat=0._ki
        !
        c_temp=-zdilog((m3_sq+s13)/s13,-1._ki)/s13+1._ki/2._ki*z_log2(-s1&
          &3,-1._ki)/s13-1._ki/4._ki*z_log2(m3_sq,-1._ki)/s13
        !
        c_temp_rat=0._ki
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        if (par3 == 1) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/s13
          !
          c_temp_d1_rat=1._ki/s13
          !
          c_temp=(he_c(1,s13,-m3_sq)*s13-2._ki+z_log(m3_sq,-1._ki)+2._ki*he&
            &_c(1,s13,-m3_sq)*m3_sq)/s13
          !
          c_temp_rat=-2._ki/s13
          !
        else if (par3 == 2) then
          !
          c_temp_d2=1._ki/2._ki/s13
          !
          c_temp_d2_rat=1._ki/2._ki/s13
          !
          c_temp_d1=z_log(-s13,-1._ki)/s13-1._ki/s13-1._ki/2._ki*z_log(m3_s&
            &q,-1._ki)/s13
          !
          c_temp_d1_rat=-1._ki/s13
          !
          c_temp=-1._ki/4._ki*(4._ki*zdilog((m3_sq+s13)/s13,-1._ki)-2._ki*z&
            &_log2(-s13,-1._ki)+8._ki*z_log(-s13,-1._ki)+z_log2(m3_sq,-1._ki&
            &)-8._ki-4._ki*z_log(m3_sq,-1._ki))/s13
          !
          c_temp_rat=2._ki/s13
          !
        else if (par3 == 3) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=he_c(1,s13,-m3_sq)
          !
          c_temp_rat=0._ki
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p1m_1mi:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Unimplemented combination of parameters.'
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
        end if
        !
      ! cas avec deux parametres de feynman au numerateur
      else if ( (par1==0) ) then
        !
        if ( (par2 == 1) .and. (par3 == 1) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/2._ki/s13
          !
          c_temp_d1_rat=1._ki/2._ki/s13
          !
          c_temp=1._ki/2._ki*(he_c(2,s13,-m3_sq)*s13**2-s13+2._ki*he_c(2,s1&
            &3,-m3_sq)*m3_sq*s13-3._ki*m3_sq+2._ki*he_c(2,s13,-m3_sq)*m3_sq*&
            &*2+z_log(m3_sq,-1._ki)*m3_sq)/s13/m3_sq
          !
          c_temp_rat=-1._ki/2._ki*(2._ki*s13+m3_sq)/(s13+m3_sq)/s13
          !
        else if ( (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp_d2=1._ki/2._ki/s13
          !
          c_temp_d2_rat=1._ki/2._ki/s13
          !
          c_temp_d1=z_log(-s13,-1._ki)/s13-3._ki/2._ki/s13-1._ki/2._ki*z_lo&
            &g(m3_sq,-1._ki)/s13
          !
          c_temp_d1_rat=-3._ki/2._ki/s13
          !
          c_temp=-1._ki/4._ki*(4._ki*zdilog((s13+m3_sq)/s13,-1._ki)-2._ki*z&
            &_log2(-s13,-1._ki)+12._ki*z_log(-s13,-1._ki)-6._ki*z_log(m3_sq,&
            &-1._ki)-14._ki+z_log2(m3_sq,-1._ki))/s13
          !
          c_temp_rat=7._ki/2._ki/s13
          !
        else if ( (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/2._ki*(he_c(2,s13,-m3_sq)*s13-1._ki)/m3_sq
          !
          c_temp_rat=-1._ki/2._ki/(s13+m3_sq)
          !
        else if ( (par2 == 1) .and. (par3 == 2) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/2._ki/s13
          !
          c_temp_d1_rat=1._ki/2._ki/s13
          !
          c_temp=1._ki/2._ki*(he_c(1,s13,-m3_sq)*s13-3._ki+z_log(m3_sq,-1._&
            &ki)+2._ki*he_c(1,s13,-m3_sq)*m3_sq)/s13
          !
          c_temp_rat=-3._ki/2._ki/s13
          !
        else if ( (par2 == 1) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/2._ki*he_c(2,s13,-m3_sq)
          !
          c_temp_rat=1._ki/2._ki/(s13+m3_sq)
          !
        else if ( (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/2._ki*he_c(1,s13,-m3_sq)
          !
          c_temp_rat=0._ki
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p1m_1mi:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Unimplemented combination of parameters.'
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
        end if
        !
      ! cas avec trois parametres de feynman au numerateur
      else
        !
        if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 1) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/3._ki/s13
          !
          c_temp_d1_rat=1._ki/3._ki/s13
          !
          c_temp=1._ki/18._ki*(6._ki*he_c(3,s13,-m3_sq)*s13**3-3._ki*s13**2&
            &+18._ki*he_c(3,s13,-m3_sq)*m3_sq*s13**2+18._ki*he_c(3,s13,-m3_s&
            &q)*m3_sq**2*s13-12._ki*m3_sq*s13+12._ki*he_c(3,s13,-m3_sq)*m3_s&
            &q**3+6._ki*z_log(m3_sq,-1._ki)*m3_sq**2-22._ki*m3_sq**2)/s13/m3&
            &_sq**2
          !
          c_temp_rat=-1._ki/18._ki/m3_sq*(6._ki*s13**2+25._ki*m3_sq*s13+16.&
            &_ki*m3_sq**2)/(s13+m3_sq)/s13
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp_d2=1._ki/2._ki/s13
          !
          c_temp_d2_rat=1._ki/2._ki/s13
          !
          c_temp_d1=z_log(-s13,-1._ki)/s13-11._ki/6._ki/s13-1._ki/2._ki*z_l&
            &og(m3_sq,-1._ki)/s13
          !
          c_temp_d1_rat=-11._ki/6._ki/s13
          !
          c_temp=1._ki/36._ki*(-36._ki*zdilog((s13+m3_sq)/s13,-1._ki)+18._k&
            &i*z_log2(-s13,-1._ki)-132._ki*z_log(-s13,-1._ki)+66._ki*z_log(m&
            &3_sq,-1._ki)+170._ki-9._ki*z_log2(m3_sq,-1._ki))/s13
          !
          c_temp_rat=85._ki/18._ki/s13
          !
        else if ( (par1 == 3) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/6._ki*(2._ki*he_c(3,s13,-m3_sq)*s13**2-s13-m3_sq)/m3&
            &_sq**2
          !
          c_temp_rat=-1._ki/6._ki/m3_sq*(2._ki*s13+m3_sq)/(s13+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 2) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/6._ki/s13
          !
          c_temp_d1_rat=1._ki/6._ki/s13
          !
          c_temp=1._ki/18._ki*(3._ki*he_c(2,s13,-m3_sq)*s13**2-3._ki*s13+6.&
            &_ki*he_c(2,s13,-m3_sq)*m3_sq*s13-11._ki*m3_sq+6._ki*he_c(2,s13,&
            &-m3_sq)*m3_sq**2+3._ki*z_log(m3_sq,-1._ki)*m3_sq)/m3_sq/s13
          !
          c_temp_rat=-1._ki/18._ki*(8._ki*s13+5._ki*m3_sq)/s13/(s13+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/3._ki/s13
          !
          c_temp_d1_rat=1._ki/3._ki/s13
          !
          c_temp=1._ki/9._ki*(3._ki*he_c(1,s13,-m3_sq)*s13-11._ki+3._ki*z_l&
            &og(m3_sq,-1._ki)+6._ki*he_c(1,s13,-m3_sq)*m3_sq)/s13
          !
          c_temp_rat=-11._ki/9._ki/s13
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/3._ki*he_c(3,s13,-m3_sq)
          !
          c_temp_rat=1._ki/6._ki/(s13+m3_sq)
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/3._ki*he_c(1,s13,-m3_sq)
          !
          c_temp_rat=0._ki
          !
        else if ( (par1 == 1) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/6._ki*(2._ki*he_c(3,s13,-m3_sq)*s13-1)/m3_sq
          !
          c_temp_rat=-1._ki/6._ki/(s13+m3_sq)
          !
        else if ( (par1 == 2) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/6._ki*(he_c(2,s13,-m3_sq)*s13-1)/m3_sq
          !
          c_temp_rat=-1._ki/6._ki/(s13+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=0._ki
          !
          c_temp_d1_rat=0._ki
          !
          c_temp=1._ki/6._ki*he_c(2,s13,-m3_sq)
          !
          c_temp_rat=1._ki/6._ki/(s13+m3_sq)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p1m_1mi:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Unimplemented combination of parameters.'
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
        end if
        !
      end if
      !
      if ( (rat_or_tot_par%tot_selected)  ) then
        !
        f3p1m_1mi(1:2) = (/real(c_temp_d2,ki),aimag(c_temp_d2)/)
        f3p1m_1mi(3:4) = (/real(c_temp_d1,ki),aimag(c_temp_d1)/)
        f3p1m_1mi(5:6) = (/real(c_temp,ki),aimag(c_temp)/)
        !
      else !if ( (rat_or_tot_par%rat_selected)  ) then
        !
        f3p1m_1mi(1:2) = (/real(c_temp_d2_rat,ki),aimag(c_temp_d2_rat)/)
        f3p1m_1mi(3:4) = (/real(c_temp_d1_rat,ki),aimag(c_temp_d1_rat)/)
        f3p1m_1mi(5:6) = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
        !
      end if
      !
      ! On change \epsilon_{ir} en -\epsilon_{uv}
      !
      f3p1m_1mi(3:4) = -f3p1m_1mi(3:4)
      !
      ! On factorise r_{\gamma}
      !
      f3p1m_1mi(5:6) = f3p1m_1mi(5:6)+pi6*f3p1m_1mi(1:2)
      !
      ! on ajoute la dependence en mu^2
      !
      lmu2 = log(mu2_scale_par)
      f3p1m_1mi(5:6) = f3p1m_1mi(5:6) + f3p1m_1mi(3:4)*lmu2 + f3p1m_1mi(1:2)*lmu2**2/2._ki      
      f3p1m_1mi(3:4) = f3p1m_1mi(3:4) + f3p1m_1mi(1:2)*lmu2
      !
    end function f3p1m_1mi
    !
    !
    !****f* src/integral/three_point/function_3p1m_1mi/f3p1m_1mi_np2
    ! NAME
    !
    !  Function f3p1m_1mi_np2
    !
    ! USAGE
    !
    !  real_dim4 = f3p1m_1mi_np2(s13,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the one off-shell external leg one internal mass three point function in n+2 dimensions. 
    !  with up to one Feynman parameter in the numerator.
    !  It retuns an array of 4 reals corresponding to the real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s13 -- real (type ki), the value of the S matrix element corresponding to the external off-shell leg
    !  * m3_sq -- real (type ki), the value of the internal mass squared
    !  * par1 -- an integer, the label of the third Feynman parameter = 0
    !  * par2 -- an integer, the label of the second Feynman parameter = 0
    !  * par3 -- an integer, the label of the first Feynman parameter
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  An real (type ki) array of rank 1 and shape 4 corresponding to 
    !  the real/imaginary part of the coefficient of the 1/epsilon term
    !  and the real/imaginary part of the constant term. If par1 and/or par2
    !  are different from zero, an error is returned.
    !
    ! EXAMPLE
    !
    !
    !*****
    function f3p1m_1mi_np2(s13,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(4) :: f3p1m_1mi_np2
      !
      complex(ki) :: c_temp,c_temp_rat
      real(ki) :: lmu2
      !
      f3p1m_1mi_np2 = 0._ki
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        f3p1m_1mi_np2(1) = -1._ki/2._ki
        f3p1m_1mi_np2(2) = 0._ki
        !
        c_temp=1._ki/2._ki*z_log(m3_sq,-1._ki)-3._ki/2._ki+1._ki/2._ki*he&
          &_c(1,s13,-m3_sq)*s13
        !
        c_temp_rat=-3._ki/2._ki
        !
        if ( (rat_or_tot_par%tot_selected)  ) then
          !
          f3p1m_1mi_np2(3:4) = (/real(c_temp,ki),aimag(c_temp)/)
          !
        else !if ( (rat_or_tot_par%rat_selected)  ) then
          !
          f3p1m_1mi_np2(3:4) = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
          !
        end if
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        f3p1m_1mi_np2(1) = -1._ki/6._ki
        f3p1m_1mi_np2(2) = 0._ki
        !
        if (par3 == 1) then
          !
          c_temp=1._ki/6._ki*z_log(m3_sq,-1._ki)-1._ki/18._ki*(3._ki*s13**2&
            &*he_c(2,s13,-m3_sq)-6._ki*he_c(1,s13,-m3_sq)*s13*m3_sq-3._ki*s1&
            &3+11._ki*m3_sq)/m3_sq
          !
          c_temp_rat=-1._ki/18._ki*(8._ki*s13+11._ki*m3_sq)/(s13+m3_sq)
          !
        else if (par3 == 2) then
          !
          c_temp=1._ki/6._ki*z_log(m3_sq,-1._ki)-11._ki/18._ki+1._ki/6._ki*&
            &he_c(1,s13,-m3_sq)*s13
          !
          c_temp_rat=-11._ki/18._ki
          !
        else if (par3 == 3) then
          !
          c_temp=1._ki/6._ki*z_log(m3_sq,-1._ki)+1._ki/18._ki*(3._ki*s13**2&
            &*he_c(2,s13,-m3_sq)-3._ki*s13-5._ki*m3_sq)/m3_sq
          !
          c_temp_rat=-1._ki/18._ki*(8._ki*s13+5._ki*m3_sq)/(s13+m3_sq)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p1m_1mi_np2:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Unimplemented combination of parameters.'
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
        end if
        !
        if ( (rat_or_tot_par%tot_selected)  ) then
          !
          f3p1m_1mi_np2(3:4) = (/real(c_temp,ki),aimag(c_temp)/)
          !
        else !if ( (rat_or_tot_par%rat_selected)  ) then
          !
          f3p1m_1mi_np2(3:4) = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
          !
        end if
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f3p1m_1mi_np2:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'no need of 3-point integrals in 6 dimension &
                          &with more than one Feynman parameter in the numerator'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'The value of Feynman parameters in argument: %d1'
        tab_erreur_par(3)%arg_int_tab = (/packb((/par1,par2,par3/)),4/)
        call catch_exception(0)
        !
        ! to please the compiler
        stop
        !
      end if
      !
      ! on ajoute la dependence en mu^2
      !
      lmu2 = log(mu2_scale_par)
      f3p1m_1mi_np2(3:4) = f3p1m_1mi_np2(3:4) + f3p1m_1mi_np2(1:2)*lmu2
      !
    end function f3p1m_1mi_np2
    !
end module function_3p1m_1mi
