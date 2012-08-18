!
!****h* src/integral/three_point/function_3p0m_1mi
! NAME
!
!  Module function_3p0m_1mi
!
! USAGE
!
!  use function_3p0m_1mi
!
! DESCRIPTION
!
!  This module is used to compute the zero off-shell external leg one internal mass three point function
!  with/without Feynman parameters in n, n+2 dimensions
!
! OUTPUT
!
!  This module exports two functions:
!  * f3p0m_1mi -- a function for the computation of the zero off-shell external leg one internal mass three 
!    point function with/without Feynman parameters in n dimensions
!  * f3p0m_1mi_np2 -- a function for the computation of the zero off-shell external leg one internal mass three 
!    point function with/without Feynman parameters in n+2 dimensions
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * logarithme (src/module/z_log.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90) only : tab_erreur_par,catch_exception
!  * parametre (src/module/parametre.f90) only : rat_or_tot_par,mu2_scale_par
!  * array (src/module/array.f90) only : packb
!
!*****
module function_3p0m_1mi
  !
  use precision_golem
  use logarithme
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre, only : rat_or_tot_par,mu2_scale_par
  use array, only : packb
  implicit none
  !
  private 
  !
  public :: f3p0m_1mi, f3p0m_1mi_np2
  !
  contains
    !
    !
    !****f* src/integral/three_point/function_3p0m_1mi/f3p0m_1mi
    ! NAME
    !
    !  Function f3p0m_1mi
    !
    ! USAGE
    !
    !  real_dim6 = f3p0m_1mi(m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the zero off-shell external leg one internal mass three point function in n dimensions
    !  with up to three Feynman parameters in the numerator.
    !  It returns an array of 6 reals corresponding to the real/imaginary
    !  part of the coefficient of the 1/epsilon^2 term, real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * m3_sq -- real (type ki), the value of the internal mass squared
    !  * par1 -- an integer, the label of the third Feynman parameter
    !  * par2 -- an integer, the label of the second Feynman parameter
    !  * par3 -- an integer, the label of the first Feynman parameter
    !  Note that par1,par2 and par3 are supposed to be ordered, i.e.
    !  par1 <= par2 <= par3, note also that put zero for par1, par2 or par3
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
    !  An real (type ki) array of rank 1 and shape 6 corresponding to 
    !  the real/imaginary part of the coefficient of the 1/epsilon^2 term,
    !  real/imaginary part of the coefficient of the 1/epsilon term
    !  and the real/imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !*****
    function f3p0m_1mi(m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(6) :: f3p0m_1mi
      !
      complex(ki) :: c_temp_d2,c_temp_d2_rat
      complex(ki) :: c_temp_d1,c_temp_d1_rat
      complex(ki) :: c_temp,c_temp_rat
      real(ki) :: lmu2
      !
      f3p0m_1mi = 0._ki
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        c_temp_d2=0._ki
        !
        c_temp_d2_rat=0._ki
        !
        c_temp_d1=1._ki/2._ki/m3_sq
        !
        c_temp_d1_rat=1._ki/2._ki/m3_sq
        !
        c_temp=(1._ki+1._ki/2._ki*z_log(m3_sq,-1._ki))/m3_sq
        !
        c_temp_rat=1._ki/m3_sq
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
          c_temp_d1=1._ki/2._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/2._ki/m3_sq
          !
          c_temp=1._ki/2._ki*z_log(m3_sq,-1._ki)/m3_sq
          !
          c_temp_rat=0._ki
          !
        else if (par3 == 2) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/2._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/2._ki/m3_sq
          !
          c_temp=1._ki/2._ki*z_log(m3_sq,-1._ki)/m3_sq
          !
          c_temp_rat=0._ki
          !
        else if (par3 == 3) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=-1._ki/2._ki/m3_sq
          !
          c_temp_d1_rat=-1._ki/2._ki/m3_sq
          !
          c_temp=(1._ki-1._ki/2._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=1._ki/m3_sq
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p0m_1mi:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Unimplemented combination of parameters.'
          tab_erreur_par(3)%chaine = 'par1 = %d0, par2,par3 = %d1'
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
          c_temp_d1=1._ki/2._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/2._ki/m3_sq
          !
          c_temp=(-1._ki/2._ki+1._ki/2._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=-1._ki/2._ki/m3_sq
          !
        else if ( (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/2._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/2._ki/m3_sq
          !
          c_temp=(-1._ki/2._ki+1._ki/2._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=-1._ki/2._ki/m3_sq
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
          c_temp=-1._ki/2._ki/m3_sq
          !
          c_temp_rat=-1._ki/2._ki/m3_sq
          !
        else if ( (par2 == 1) .and. (par3 == 2) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/4._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/4._ki/m3_sq
          !
          c_temp=(-1._ki/4._ki+1._ki/4._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=-1._ki/4._ki/m3_sq
          !
        else if ( (par2 == 1) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=-1._ki/4._ki/m3_sq
          !
          c_temp_d1_rat=-1._ki/4._ki/m3_sq
          !
          c_temp=(3._ki/4._ki-1._ki/4._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=3._ki/4._ki/m3_sq
          !
        else if ( (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=-1._ki/4._ki/m3_sq
          !
          c_temp_d1_rat=-1._ki/4._ki/m3_sq
          !
          c_temp=(3._ki/4._ki-1._ki/4._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=3._ki/4._ki/m3_sq
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p0m_1mi:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Unimplemented combination of parameters.'
          tab_erreur_par(3)%chaine = 'par1 = %d0, par2,par3 = %d1'
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
          c_temp_d1=1._ki/2._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/2._ki/m3_sq
          !
          c_temp=(-5._ki/6._ki+1._ki/2._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=-5._ki/6._ki/m3_sq
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/2._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/2._ki/m3_sq
          !
          c_temp=(-5._ki/6._ki+1._ki/2._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=-5._ki/6._ki/m3_sq
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
          c_temp=-1._ki/6._ki/m3_sq
          !
          c_temp_rat=-1._ki/6._ki/m3_sq
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 2) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/6._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/6._ki/m3_sq
          !
          c_temp=(-5._ki/18._ki+1._ki/6._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=-5._ki/18._ki/m3_sq
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=1._ki/6._ki/m3_sq
          !
          c_temp_d1_rat=1._ki/6._ki/m3_sq
          !
          c_temp=(-5._ki/18._ki+1._ki/6._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=-5._ki/18._ki/m3_sq
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=-1._ki/6._ki/m3_sq
          !
          c_temp_d1_rat=-1._ki/6._ki/m3_sq
          !
          c_temp=(11._ki/18._ki-1._ki/6._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=11._ki/18._ki/m3_sq
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=-1._ki/6._ki/m3_sq
          !
          c_temp_d1_rat=-1._ki/6._ki/m3_sq
          !
          c_temp=(11._ki/18._ki-1._ki/6._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=11._ki/18._ki/m3_sq
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
          c_temp=-1._ki/6._ki/m3_sq
          !
          c_temp_rat=-1._ki/6._ki/m3_sq
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
          c_temp=-1._ki/6._ki/m3_sq
          !
          c_temp_rat=-1._ki/6._ki/m3_sq
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp_d2=0._ki
          !
          c_temp_d2_rat=0._ki
          !
          c_temp_d1=-1._ki/12._ki/m3_sq
          !
          c_temp_d1_rat=-1._ki/12._ki/m3_sq
          !
          c_temp=(11._ki/36._ki-1._ki/12._ki*z_log(m3_sq,-1._ki))/m3_sq
          !
          c_temp_rat=11._ki/36._ki/m3_sq
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p0m_1mi:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Unimplemented combination of parameters.'
          tab_erreur_par(3)%chaine = 'par1 = %d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
          !
        end if
        !
      end if
      !
      if ( (rat_or_tot_par%tot_selected)  ) then
        !
        f3p0m_1mi(1:2) = (/real(c_temp_d2,ki),aimag(c_temp_d2)/)
        f3p0m_1mi(3:4) = (/real(c_temp_d1,ki),aimag(c_temp_d1)/)
        f3p0m_1mi(5:6) = (/real(c_temp,ki),aimag(c_temp)/)
        !
      else !if ( (rat_or_tot_par%rat_selected)  ) then
        !
        f3p0m_1mi(1:2) = (/real(c_temp_d2_rat,ki),aimag(c_temp_d2_rat)/)
        f3p0m_1mi(3:4) = (/real(c_temp_d1_rat,ki),aimag(c_temp_d1_rat)/)
        f3p0m_1mi(5:6) = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
        !
      end if
      !
      ! On change \epsilon_{ir} en -\epsilon_{uv}
      !
      f3p0m_1mi(3:4) = -f3p0m_1mi(3:4)
      !
      ! on ajoute la dependence en mu^2
      !
      lmu2 = log(mu2_scale_par)
      f3p0m_1mi(5:6) = f3p0m_1mi(5:6) + f3p0m_1mi(3:4)*lmu2 + f3p0m_1mi(1:2)*lmu2**2/2._ki      
      f3p0m_1mi(3:4) = f3p0m_1mi(3:4) + f3p0m_1mi(1:2)*lmu2
      !
    end function f3p0m_1mi
    !
    !
    !****f* src/integral/three_point/function_3p0m_1mi/f3p0m_1mi_np2
    ! NAME
    !
    !  Function f3p0m_1mi_np2
    !
    ! USAGE
    !
    !  real_dim4 = f3p0m_1mi_np2(m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the zero off-shell external leg one internal mass three point function in n+2 dimensions. 
    !  with up to one Feynman parameter in the numerator.
    !  It retuns an array of 4 reals corresponding to the real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * m3_sq -- real (type ki), the value of the internal mass squared
    !  * par1 -- an integer, the label of the third Feynman parameter = 0
    !  * par2 -- an integer, the label of the second Feynman parameter = 0
    !  * par3 -- an integer, the label of the first Feynman parameter
    !  Note that par1,par2 and par3 are supposed to be ordered, i.e.
    !  par1 <= par2 <= par3, note also that put zero for par1, par2 or par3
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
    !  An real (type ki) array of rank 1 and shape 4 corresponding to 
    !  the real/imaginary part of the coefficient of the 1/epsilon term
    !  and the real/imaginary part of the constant term. If par1 and/or par2
    !  are different from zero, an error is returned.
    !
    ! EXAMPLE
    !
    !
    !*****
    function f3p0m_1mi_np2(m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(4) :: f3p0m_1mi_np2
      !
      complex(ki) :: c_temp,c_temp_rat
      real(ki) :: lmu2
      !
      f3p0m_1mi_np2 = 0._ki
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        f3p0m_1mi_np2(1) = -1._ki/2._ki
        f3p0m_1mi_np2(2) = 0._ki
        !
        c_temp=-3._ki/2._ki+1._ki/2._ki*z_log(m3_sq,-1._ki)
        !
        c_temp_rat=-3._ki/2._ki
        !
        if ( (rat_or_tot_par%tot_selected)  ) then
          !
          f3p0m_1mi_np2(3:4) = (/real(c_temp,ki),aimag(c_temp)/)
          !
        else !if ( (rat_or_tot_par%rat_selected)  ) then
          !
          f3p0m_1mi_np2(3:4) = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
          !
        end if
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        f3p0m_1mi_np2(1) = -1._ki/6._ki
        f3p0m_1mi_np2(2) = 0._ki
        !
        if (par3 == 1) then !changed: 11.08.10
          !
          c_temp=-11._ki/18._ki+1._ki/6._ki*z_log(m3_sq,-1._ki)
          !
          c_temp_rat=-11._ki/18._ki
          !
        else if (par3 == 2) then !changed: 11.08.10
          !
          c_temp=-11._ki/18._ki+1._ki/6._ki*z_log(m3_sq,-1._ki)
          !
          c_temp_rat=-11._ki/18._ki
          !
        else if (par3 == 3) then
          !
          c_temp=-5._ki/18._ki+1._ki/6._ki*z_log(m3_sq,-1._ki)
          !
          c_temp_rat=-5._ki/18._ki
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p0m_1mi_np2:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'Unimplemented combination of parameters.'
          tab_erreur_par(3)%chaine = 'par1 = %d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          ! to please the compiler
          stop
          !
          !
        end if
        !
        if ( (rat_or_tot_par%tot_selected)  ) then
          !
          f3p0m_1mi_np2(3:4) = (/real(c_temp,ki),aimag(c_temp)/)
          !
        else !if ( (rat_or_tot_par%rat_selected)  ) then
          !
          f3p0m_1mi_np2(3:4) = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
          !
        end if
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f3p0m_1mi_np2:'
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
      f3p0m_1mi_np2(3:4) = f3p0m_1mi_np2(3:4) + f3p0m_1mi_np2(1:2)*lmu2
      !
    end function f3p0m_1mi_np2
    !
end module function_3p0m_1mi
