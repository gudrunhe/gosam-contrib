!
!****h* src/integral/three_point/function_3p2m_1mi
! NAME
!
!  Module function_3p2m_1mi
!
! USAGE
!
!  use function_3p2m_1mi
!
! DESCRIPTION
!
!  This module is used to compute the two off-shell external leg one internal mass three point function
!  with/without Feynman parameters in n, n+2 dimensions
!
! OUTPUT
!
!  This module exports two functions:
!  * f3p2m_1mi -- a function for the computation of the two off-shell external leg one internal mass three 
!    point function with/without Feynman parameters in n dimensions
!  * f3p2m_1mi_np2 -- a function for the computation of the two off-shell external leg one internal mass three 
!    point function with/without Feynman parameters in n+2 dimensions
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * logarithme (src/module/z_log.f90)
!  * dilogarithme (src/module/zdilog.f90)
!  * func_he (src/integrals/three_point/mod_he.f90)
!  * func_hf (src/integrals/three_point/mod_hf.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!         only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
!  * constante (src/module/constante.f90) only : i_,un
!  * parametre (src/module/parametre.f90)
!         only : coupure_3p2m_1mi,rat_or_tot_par,tolerance,alpha_par,beta_par,lambda_par,mu2_scale_par
!  * array (src/module/array.f90) only : packb
!  * numerical_evaluation (src/numerical/mod_numeric.f90) only : generic_eval_numer
!  * matrice_s (src/kinematics/matrice_s.f90)
!
!*****
module function_3p2m_1mi
  !
  use precision_golem
  use logarithme
  use dilogarithme
  use func_he
  use func_hf
  use sortie_erreur, only : tab_erreur_par,catch_exception,origine_info_par,num_grand_b_info_par,denom_grand_b_info_par
  use constante, only : i_, un, czero
  use parametre, only : coupure_3p2m_1mi,rat_or_tot_par,tolerance,alpha_par,beta_par,lambda_par,mu2_scale_par
  use array, only : packb
  use numerical_evaluation, only : generic_eval_numer
  use matrice_s, only: find_plus_grand
  implicit none
  !
  private
  complex(ki) :: s13_glob,m3_sq_glob,s23_glob
  real(ki) :: eps_glob
  integer :: par1_glob,par2_glob,par3_glob
  character (len=3) :: dim_glob
  integer, dimension(3) :: par
  !
  interface f3p2m_1mi
     module procedure f3p2m_1mi_r
     module procedure f3p2m_1mi_c
  end interface
  !
  interface f3p2m_1mi_np2
     module procedure f3p2m_1mi_np2_r
     module procedure f3p2m_1mi_np2_c
  end interface
  !
!
  public :: f3p2m_1mi, f3p2m_1mi_np2
  !
  contains
    !
    !****f* src/integral/three_point/function_3p2m_1mi/f3p2m_1mi
    ! NAME
    !
    !  Function f3p2m_1mi
    !
    ! USAGE
    !
    !  real_dim6 = f3p2m_1mi(s23,s13,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the two off-shell external leg one internal mass three point function in n dimensions
    !  with up to three Feynman parameters in the numerator.
    !  It retuns an array of 6 reals corresponding to the real/imaginary
    !  part of the coefficient of the 1/epsilon^2 term, real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s23 -- real/complex (type ki), the value of the S matrix element corresponding to the first external off-shell leg
    !  * s13 -- real/complex (type ki), the value of the S matrix element corresponding to the second external off-shell leg
    !  * m3_sq -- real/complex (type ki), the value of the internal mass squared
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
    !
    function f3p2m_1mi_r(s23,s13,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(6) :: f3p2m_1mi_r
      !
      real(ki) :: lamb
      real(ki) :: plus_grand
      complex(ki) :: resto,abserro
      real(ki) :: as23,as13,am3_sq
      !
      !
      ! on redefinit la matrice S de telle facon a ce que ces elements
      ! soient entre -1 et 1
      !
      if (rat_or_tot_par%tot_selected) then
        !
        plus_grand = find_plus_grand((/ abs(s13),abs(s23),abs(m3_sq)/))
        !
      else !if (rat_or_tot_par%rat_selected) then
        !
        plus_grand = 1._ki
        !
      end if
      !
      as13 = s13/plus_grand
      as23 = s23/plus_grand
      am3_sq = m3_sq/plus_grand
      !
      lamb = as13-as23
      !
      f3p2m_1mi_r(:) = 0._ki
      !
      ! the correction for plus_grand are taken into account in he and hf
      !
      f3p2m_1mi_r = a3p2m_1mi_div_r(s23,s13,m3_sq,par1,par2,par3)
      !
      if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_3p2m_1mi) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f3p2m_1mi (in file function_3p2m_1mi.f90):&
        &the flag rat to compute the rational part is on &
        &and the program reaches a region of phase space in &
        &which det(G) = 0 . Be careful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to go on, he has to &
        &reduce the value of the parameter coupure_3p2m_1mi'
        call catch_exception(0)
        !
        stop
        !
      end if
      !
      if (abs(lamb) > coupure_3p2m_1mi) then
        !
        ! analytic computation
        !
        f3p2m_1mi_r(5:6) = f3p2m_1mi_r(5:6) + a3p2m_1mi_r(as23,as13,am3_sq,par1,par2,par3)&
                    &/plus_grand
        !
      else
        !
        ! numerical computation
        !
        dim_glob = "ndi"
        par1_glob = par1
        par2_glob = par2
        par3_glob = par3
        !
        s23_glob = cmplx(as23,0._ki,ki)
        s13_glob = cmplx(as13,0._ki,ki)
        m3_sq_glob = cmplx(am3_sq,0._ki,ki)
        !
        resto = 0._ki
        abserro = 0._ki
        !
        ! on pose z = x - i*eps*y (avec x et y > 0)
        ! z*s23+(1-z)*s23 = s23+x*(s23-s23)-i*eps*y*(s23-s23)
        ! on veut la partie imaginaire du meme signe que i*lambda
        ! => eps*(s23-s23) < 0
        !
        ! faire attention que suivant le signe de eps_glob, on tourne dans le
        ! sens des aiguilles d'une montre ou inversement
        ! eps_glob = 1, on ferme le contour vers le bas --> -2 i Pi residu
        ! eps_glob = -1, on ferme le contour vers le haut --> +2 i Pi residu
        !
        eps_glob = sign(1._ki,as23-as13)
        !
        origine_info_par = "f3p2m_1mi, dimension "
        num_grand_b_info_par = lamb
        denom_grand_b_info_par = 1._ki
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,resto,abserro)
        !
        resto = resto/plus_grand
        !
        f3p2m_1mi_r(5) = f3p2m_1mi_r(5) + real(resto,ki)
        f3p2m_1mi_r(6) = f3p2m_1mi_r(6) + aimag(resto)
        !
      end if
      !
      ! la dependance en mu2 se fait a travers les fonctions he,hf
      ! inutile de l'ajouter
      !
    end function f3p2m_1mi_r
    !
   function f3p2m_1mi_c(s23,s13,m3_sq,par1,par2,par3)
      !
      complex(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(6) :: f3p2m_1mi_c
      !
      complex(ki) :: lamb
      real(ki) :: plus_grand
      complex(ki) :: as23,as13,am3_sq
      complex(ki) :: resto,abserro
      !
      !
      ! We divide by the maximal real or imaginary value to get the real and 
      ! imaginary entries between -1 and 1
      !
      if (rat_or_tot_par%tot_selected) then
         !
         plus_grand = find_plus_grand(abs( (/ real( (/ s13, s23, m3_sq /), ki ) , &
              &                              aimag( (/ s13, s23, m3_sq /) ) /))) 
         !
      else !if (rat_or_tot_par%rat_selected) then
        !
        plus_grand = 1._ki
        !
      end if
      !
      as13 = s13/plus_grand
      as23 = s23/plus_grand
      am3_sq = m3_sq/plus_grand
      !
      lamb = as13-as23
      !
      f3p2m_1mi_c(:) = 0._ki
      !
      ! the correction for plus_grand are taken into account in he and hf
      !
      f3p2m_1mi_c = a3p2m_1mi_div_c(s23,s13,m3_sq,par1,par2,par3)
      !
      if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_3p2m_1mi) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = &
        &'In function f3p2m_1mi (in file function_3p2m_1mi.f90):&
        &the flag rat to compute the rational part is on &
        &and the program reaches a region of phase space in &
        &which det(G) = 0 . Be careful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to go on, he has to &
        &reduce the value of the parameter coupure_3p2m_1mi'
        call catch_exception(0)
        !
        stop
        !
      end if
      !
      if (abs(lamb) > coupure_3p2m_1mi) then
        !
        ! analytic computation
        !
        f3p2m_1mi_c(5:6) = f3p2m_1mi_c(5:6) + a3p2m_1mi_c(as23,as13,am3_sq,par1,par2,par3)&
                   &/plus_grand
        !
      else
        !
        ! numerical computation
        !
        dim_glob = "ndi"
        par1_glob = par1
        par2_glob = par2
        par3_glob = par3
        !
        s23_glob = as23
        s13_glob = as13
        m3_sq_glob = am3_sq
        !
        resto = 0._ki
        abserro = 0._ki
        !
        ! on pose z = x - i*eps*y (avec x et y > 0)
        ! z*s13+(1-z)*s23 = s23+x*(s13-s23)-i*eps*y*(s13-s23)
        ! now s13 and s23 are complex BUT HAVE THE SAME IMAGINARY PART
        ! i.e. s13-s23 is real. 
        ! We want the the argument of the log never cross the cut, that is to say that the
        ! sign(Im(arg_log)) is constant along the contour
        ! => eps = -sign(Im(s23))*sign(s13-s23)
        ! Note that with this prescription we avoid the pole when z*s13+(1-z)*s23=0
        !
        ! faire attention que suivant le signe de eps_glob, on tourne dans le
        ! sens des aiguilles d'une montre ou inversement
        ! eps_glob = 1, on ferme le contour vers le bas --> -2 i Pi residu
        ! eps_glob = -1, on ferme le contour vers le haut --> +2 i Pi residu
        !
        eps_glob = -sign(1._ki,aimag(s23_glob))*sign(1._ki,real(s13_glob-s23_glob,ki))
        !
        origine_info_par = "f3p2m_1mi_c, dimension "
        num_grand_b_info_par = lamb
        denom_grand_b_info_par = 1._ki
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,resto,abserro)
        !
        resto = resto/plus_grand
        !
        f3p2m_1mi_c(5) = f3p2m_1mi_c(5) + real(resto,ki)
        f3p2m_1mi_c(6) = f3p2m_1mi_c(6) + aimag(resto)
        !
      end if
      !
      ! la dependance en mu2 se fait a travers les fonctions he,hf
      ! inutile de l'ajouter
      !
      !
    end function f3p2m_1mi_c
    !
    !****f* src/integral/three_point/function_3p2m_1mi/f3p2m_1mi_np2
    ! NAME
    !
    !  Function f3p2m_1mi_np2
    !
    ! USAGE
    !
    !  real_dim4 = f3p2m_1mi_np2(s23,s13,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the two off-shell external leg one internal mass three point function in n+2 dimensions
    !  with up to three Feynman parameters in the numerator.
    !  It retuns an array of 4 reals corresponding to the real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s23 -- real/complex (type ki), the value of the S matrix element corresponding to the first external off-shell leg
    !  * s13 -- real/complex (type ki), the value of the S matrix element corresponding to the second external off-shell leg
    !  * m3_sq -- real/complex (type ki), the value of the internal mass squared
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
    !  An real (type ki) array of rank 1 and shape 4 corresponding to 
    !  the real/imaginary part of the coefficient of the 1/epsilon term
    !  and the real/imaginary part of the constant term. If par1 and/or par2
    !  are different from zero, an error is returned.
    !
    ! EXAMPLE
    !
    !
    !*****
    !
    function f3p2m_1mi_np2_r(s23,s13,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(4) :: f3p2m_1mi_np2_r
      !
      integer :: nb_par
      real(ki) :: lamb
      real(ki) :: plus_grand
      real(ki) :: norma
      complex(ki) :: resto,abserro
      real(ki) :: as23,as13,am3_sq
      real(ki) :: lmu2
      !
      par = (/par1,par2,par3/)
      !
      !
      ! on redefinit la matrice S de telle facon a ce que ces elements
      ! soient entre -1 et 1
      !
      if (rat_or_tot_par%tot_selected) then
        !
        plus_grand = max(abs(s13),abs(s23),abs(m3_sq))
        !
      else !if (rat_or_tot_par%rat_selected) then
        !
        plus_grand = 1._ki
        !
      end if
      !
      as13 = s13/plus_grand
      as23 = s23/plus_grand
      am3_sq = m3_sq/plus_grand
      !
      lamb = as13-as23
      !
      nb_par = count(mask=par/=0)
      !
      if (nb_par == 0) then
        !
        norma = -1._ki/2._ki
        !
      else if (nb_par == 1) then
        !
        norma = -1._ki/6._ki
        !
      else
        !
        norma = 0._ki
        !
      end if
      !
      !
      f3p2m_1mi_np2_r(:) = 0._ki
      !
      f3p2m_1mi_np2_r(1:2) = (/ norma, 0._ki /)
      !
      if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_3p2m_1mi) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f3p2m_1mi (in file function_3p2m_1mi.f90):&
        &the flag rat to compute the rational part is on &
        &and the program reaches a region of phase space in &
        &which det(G) = 0 . Be careful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to go on, he has to &
        &reduce the value of the parameter coupure_3p2m_1mi'
        call catch_exception(0)
        !
        stop
        !
      end if
      !
      if (abs(lamb) > coupure_3p2m_1mi) then
        !
        ! analytic computation
        !
        f3p2m_1mi_np2_r(3:4) = a3p2m_1mi_np2_r(as23,as13,am3_sq,par1,par2,par3)
        f3p2m_1mi_np2_r(3) = f3p2m_1mi_np2_r(3)-log(plus_grand)*norma
        !
      else
        !
        ! numerical computation
        !
        dim_glob = "n+2"
        par1_glob = par1
        par2_glob = par2
        par3_glob = par3
        !
        s23_glob = cmplx(as23,0._ki,ki)
        s13_glob = cmplx(as13,0._ki,ki)
        m3_sq_glob = cmplx(am3_sq,0._ki,ki)
        !
        resto = 0._ki
        abserro = 0._ki
        !
        ! on pose z = x - i*eps*y (avec x et y > 0)
        ! z*s23+(1-z)*s23 = s23+x*(s23-s23)-i*eps*y*(s23-s23)
        ! on veut la partie imaginaire du meme signe que i*lambda
        ! => eps*(s23-s23) < 0
        !
        ! faire attention que suivant le signe de eps_glob, on tourne dans le
        ! sens des aiguilles d'une montre ou inversement
        ! eps_glob = 1, on ferme le contour vers le bas --> -2 i Pi residu
        ! eps_glob = -1, on ferme le contour vers le haut --> +2 i Pi residu
        !
        eps_glob = sign(1._ki,as23-as13)
        !
        origine_info_par = "f3p2m_1mi_np2, dimension "
        num_grand_b_info_par = lamb
        denom_grand_b_info_par = 1._ki
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,resto,abserro)
        !
        resto = resto-log(plus_grand)*norma
        !
        f3p2m_1mi_np2_r(3) = real(resto,ki)
        f3p2m_1mi_np2_r(4) = aimag(resto)
        !
      end if
      !
      ! on ajoute la dependence en mu^2
      !
      lmu2 = log(mu2_scale_par)
      f3p2m_1mi_np2_r(3:4) = f3p2m_1mi_np2_r(3:4) + f3p2m_1mi_np2_r(1:2)*lmu2
      !
    end function f3p2m_1mi_np2_r
    !
    !
    function f3p2m_1mi_np2_c(s23,s13,m3_sq,par1,par2,par3)
      !
      complex(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(4) :: f3p2m_1mi_np2_c
      !
      integer :: nb_par
      complex(ki) :: lamb
      complex(ki) :: resto, abserro
      real(ki) :: plus_grand
      real(ki) :: norma
      complex(ki) :: as23,as13,am3_sq
      real(ki) :: lmu2
      !
      par = (/par1,par2,par3/)
      !
      !
      ! We divide by the maximal real or imaginary value to get the real and 
      ! imaginary entries between -1 and 1
      !
      if (rat_or_tot_par%tot_selected) then
         !
         plus_grand = max ( maxval ( abs( real( (/ s13, s23, m3_sq /), ki ) ) ), &
              &              maxval ( abs( aimag( (/ s13, s23, m3_sq /) ) ) ) )
         !
      else !if (rat_or_tot_par%rat_selected) then
        !
        plus_grand = 1._ki
        !
      end if
      !
      as13 = s13/plus_grand
      as23 = s23/plus_grand
      am3_sq = m3_sq/plus_grand
      !
      lamb = as13-as23
      !
      nb_par = count(mask=par/=0)
      !
      if (nb_par == 0) then
        !
        norma = -1._ki/2._ki
        !
      else if (nb_par == 1) then
        !
        norma = -1._ki/6._ki
        !
      else
        !
        norma = 0._ki
        !
      end if
      !
      !
      f3p2m_1mi_np2_c(:) = 0._ki
      !
      f3p2m_1mi_np2_c(1:2) = (/ norma, 0._ki /)
      !
      if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_3p2m_1mi) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f3p2m_1mi (in file function_3p2m_1mi.f90):&
        &the flag rat to compute the rational part is on &
        &and the program reaches a region of phase space in &
        &which det(G) = 0 . Be careful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to go on, he has to &
        &reduce the value of the parameter coupure_3p2m_1mi'
        call catch_exception(0)
        !
        stop
        !
      end if
      !
       if (abs(lamb) > coupure_3p2m_1mi) then 
        !
        ! analytic computation
        !
        f3p2m_1mi_np2_c(3:4) = a3p2m_1mi_np2_c(as23,as13,am3_sq,par1,par2,par3)
        f3p2m_1mi_np2_c(3) = f3p2m_1mi_np2_c(3)-log(plus_grand)*norma
        !
      else
        !
        ! numerical computation
        !
        dim_glob = "n+2"
        par1_glob = par1
        par2_glob = par2
        par3_glob = par3
        !
        s23_glob = as23
        s13_glob = as13
        m3_sq_glob = am3_sq
        !
        resto = 0._ki
        abserro = 0._ki
        !
        ! on pose z = x - i*eps*y (avec x et y > 0)
        ! z*s13+(1-z)*s23 = s23+x*(s13-s23)-i*eps*y*(s13-s23)
        ! now s13 and s23 are complex BUT HAVE THE SAME IMAGINARY PART
        ! i.e. s13-s23 is real. 
        ! We want the the argument of the log never cross the cut, that is to say that the
        ! sign(Im(arg_log)) is constant along the contour
        ! => eps = -sign(Im(s23))*sign(s13-s23)
        ! Note that with this prescription we avoid the pole when z*s13+(1-z)*s23=0
        !
        ! faire attention que suivant le signe de eps_glob, on tourne dans le
        ! sens des aiguilles d'une montre ou inversement
        ! eps_glob = 1, on ferme le contour vers le bas --> -2 i Pi residu
        ! eps_glob = -1, on ferme le contour vers le haut --> +2 i Pi residu
        !
        eps_glob = -sign(1._ki,aimag(s23_glob))*sign(1._ki,real(s13_glob-s23_glob,ki))
        !
        origine_info_par = "f3p2m_1mi_np2_c, dimension "
        num_grand_b_info_par = lamb
        denom_grand_b_info_par = 1._ki
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,resto,abserro)
        !
        resto = resto-log(plus_grand)*norma
        !
        f3p2m_1mi_np2_c(3) = real(resto,ki)
        f3p2m_1mi_np2_c(4) = aimag(resto)
        !
      end if
      !
      ! on ajoute la dependence en mu^2
      !
      lmu2 = log(mu2_scale_par)
      f3p2m_1mi_np2_c(3:4) = f3p2m_1mi_np2_c(3:4) + f3p2m_1mi_np2_c(1:2)*lmu2
      !
    end function f3p2m_1mi_np2_c

    !****if* src/integral/three_point/function_3p2m_1mi/a3p2m_1mi_div
    ! NAME
    !
    !  Function a3p2m_1mi_div
    !
    ! USAGE
    !
    !  real_dim6 = a3p2m_1mi_div(s23,s13,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the divergent part of the two off-shell external leg three point function in n dimensions
    !  with up to three Feynman parameters in the numerator.
    !  It retuns an array of 6 reals corresponding to the real/imaginary
    !  part of the coefficient of the 1/epsilon^2 term, real/imaginary part of the
    !  coefficient of the 1/epsilon term and the real/imaginary part of the 
    !  constant term.
    !
    ! INPUTS
    !
    !  * s23 -- real/complex (type ki), the value of the S matrix element corresponding to the first external off-shell leg
    !  * s13 -- real/complex (type ki), the value of the S matrix element corresponding to the second external off-shell leg
    !  * m3_sq -- real/complex (type ki), the value of the internal mass squared
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
    !
    function a3p2m_1mi_div_r(s23,s13,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(6) :: a3p2m_1mi_div_r
      !
      a3p2m_1mi_div_r(:) = 0._ki
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        a3p2m_1mi_div_r(3:4)=he(1,s13,s23)
        !
        a3p2m_1mi_div_r(5:6)=hf(1,s13,s23)
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        if (par3 == 1) then
          !
          a3p2m_1mi_div_r(3:4)=he(2,s13,s23)
          !
          a3p2m_1mi_div_r(5:6)=hf(2,s13,s23)
          !
        else if (par3 == 2) then
          !
          a3p2m_1mi_div_r(3:4)=he(1,s13,s23)-he(2,s13,s23)
          !
          a3p2m_1mi_div_r(5:6)=hf(1,s13,s23)-hf(2,s13,s23)
          !
        else if (par3 == 3) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        end if
        !
      ! cas avec deux parametres de feynman au numerateur
      else if ( (par1==0) ) then
        !
        if ( (par2 == 1) .and. (par3 == 1) ) then
          !
          a3p2m_1mi_div_r(3:4)=he(3,s13,s23)
          !
          a3p2m_1mi_div_r(5:6)=hf(3,s13,s23)
          !
        else if ( (par2 == 2) .and. (par3 == 2) ) then
          !
          a3p2m_1mi_div_r(3:4)=he(1,s13,s23)-2._ki*he(2,s13,s23)+he(3,s13,s23&
            &)
          !
          a3p2m_1mi_div_r(5:6)=hf(1,s13,s23)-2._ki*hf(2,s13,s23)+hf(3,s13,s23&
            &)
          !
        else if ( (par2 == 3) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        else if ( (par2 == 1) .and. (par3 == 2) ) then
          !
          a3p2m_1mi_div_r(3:4)=he(2,s13,s23)-he(3,s13,s23)
          !
          a3p2m_1mi_div_r(5:6)=hf(2,s13,s23)-hf(3,s13,s23)
          !
        else if ( (par2 == 1) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        else if ( (par2 == 2) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        end if
        !
      ! cas avec trois parametres de feynman au numerateur
      else
        !
        if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 1) ) then
          !
          a3p2m_1mi_div_r(3:4)=he(4,s13,s23)
          !
          a3p2m_1mi_div_r(5:6)=hf(4,s13,s23)
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          a3p2m_1mi_div_r(3:4)=he(1,s13,s23)-3._ki*he(2,s13,s23)-he(4,s13,s23&
            &)+3._ki*he(3,s13,s23)
          !
          a3p2m_1mi_div_r(5:6)=hf(1,s13,s23)-3._ki*hf(2,s13,s23)-hf(4,s13,s23&
            &)+3._ki*hf(3,s13,s23)
          !
        else if ( (par1 == 3) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 2) ) then
          !
          a3p2m_1mi_div_r(3:4)=-he(4,s13,s23)+he(3,s13,s23)
          !
          a3p2m_1mi_div_r(5:6)=-hf(4,s13,s23)+hf(3,s13,s23)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          a3p2m_1mi_div_r(3:4)=he(2,s13,s23)+he(4,s13,s23)-2._ki*he(3,s13,s23&
            &)
          !
          a3p2m_1mi_div_r(5:6)=hf(2,s13,s23)+hf(4,s13,s23)-2._ki*hf(3,s13,s23&
            &)
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        else if ( (par1 == 1) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        else if ( (par1 == 2) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          a3p2m_1mi_div_r(3:4)=0._ki
          !
          a3p2m_1mi_div_r(5:6)=0._ki
          !
        end if
        !
      end if
      !
      ! On change \epsilon_{ir} en -\epsilon_{uv}
      !
      a3p2m_1mi_div_r(3:4) = -a3p2m_1mi_div_r(3:4)
      !
    end function a3p2m_1mi_div_r
    !
    !
    function a3p2m_1mi_div_c(s23,s13,m3_sq,par1,par2,par3)
      !
      complex(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(6) :: a3p2m_1mi_div_c
      !
      a3p2m_1mi_div_c(:) = 0._ki
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then  
        !
        a3p2m_1mi_div_c(3:4)=he(1,s13,s23)
        !
        a3p2m_1mi_div_c(5:6)=hf(1,s13,s23)
        !
      ! cas avec un parametre de feynman au numerateur 
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        if (par3 == 1) then  
          !
          a3p2m_1mi_div_c(3:4)=he(2,s13,s23)
          !
          a3p2m_1mi_div_c(5:6)=hf(2,s13,s23)
          !
        else if (par3 == 2) then  
          !
          a3p2m_1mi_div_c(3:4)=he(1,s13,s23)-he(2,s13,s23)
          !
          a3p2m_1mi_div_c(5:6)=hf(1,s13,s23)-hf(2,s13,s23)
          !
        else if (par3 == 3) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        end if
        !
      ! cas avec deux parametres de feynman au numerateur
      else if ( (par1==0) ) then
        !
        if ( (par2 == 1) .and. (par3 == 1) ) then  
          !
          a3p2m_1mi_div_c(3:4)=he(3,s13,s23)
          !
          a3p2m_1mi_div_c(5:6)=hf(3,s13,s23)
          !
        else if ( (par2 == 2) .and. (par3 == 2) ) then  
          !
          a3p2m_1mi_div_c(3:4)=he(1,s13,s23)-2._ki*he(2,s13,s23)+he(3,s13,s23&
            &)
          !
          a3p2m_1mi_div_c(5:6)=hf(1,s13,s23)-2._ki*hf(2,s13,s23)+hf(3,s13,s23&
            &)
          !
        else if ( (par2 == 3) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        else if ( (par2 == 1) .and. (par3 == 2) ) then  
          !
          a3p2m_1mi_div_c(3:4)=he(2,s13,s23)-he(3,s13,s23)
          !
          a3p2m_1mi_div_c(5:6)=hf(2,s13,s23)-hf(3,s13,s23)
          !
        else if ( (par2 == 1) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        else if ( (par2 == 2) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        end if
        !
      ! cas avec trois parametres de feynman au numerateur
      else
        !
        if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 1) ) then  
          !
          a3p2m_1mi_div_c(3:4)=he(4,s13,s23)
          !
          a3p2m_1mi_div_c(5:6)=hf(4,s13,s23)
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 2) ) then  
          !
          a3p2m_1mi_div_c(3:4)=he(1,s13,s23)-3._ki*he(2,s13,s23)-he(4,s13,s23&
            &)+3._ki*he(3,s13,s23)
          !
          a3p2m_1mi_div_c(5:6)=hf(1,s13,s23)-3._ki*hf(2,s13,s23)-hf(4,s13,s23&
            &)+3._ki*hf(3,s13,s23)
          !
        else if ( (par1 == 3) .and. (par2 == 3) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 2) ) then  
          !
          a3p2m_1mi_div_c(3:4)=-he(4,s13,s23)+he(3,s13,s23)
          !
          a3p2m_1mi_div_c(5:6)=-hf(4,s13,s23)+hf(3,s13,s23)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 2) ) then  
          !
          a3p2m_1mi_div_c(3:4)=he(2,s13,s23)+he(4,s13,s23)-2._ki*he(3,s13,s23&
            &)
          !
          a3p2m_1mi_div_c(5:6)=hf(2,s13,s23)+hf(4,s13,s23)-2._ki*hf(3,s13,s23&
            &)
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        else if ( (par1 == 1) .and. (par2 == 3) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        else if ( (par1 == 2) .and. (par2 == 3) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 3) ) then  
          !
          a3p2m_1mi_div_c(3:4)=0._ki
          !
          a3p2m_1mi_div_c(5:6)=0._ki
          !
        end if
        !
      end if
      !
      ! On change \epsilon_{ir} en -\epsilon_{uv}
      !
      a3p2m_1mi_div_c(3:4) = -a3p2m_1mi_div_c(3:4)
      !
    end function a3p2m_1mi_div_c
    !
    !****if* src/integral/three_point/function_3p2m_1mi/a3p2m_1mi
    ! NAME
    !
    !  Function a3p2m_1mi
    !
    ! USAGE
    !
    !  real_dim2 = a3p2m_1mi(s23,s13,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the non divergent part two off-shell external leg three point function in n dimensions
    !  with up to three Feynman parameters in the numerator.
    !  It retuns an array of 2 reals corresponding to the real/imaginary
    !  part of the constant term.
    !
    ! INPUTS
    !
    !  * s23 -- real/complex (type ki), the value of the S matrix element corresponding to the first external off-shell leg
    !  * s13 -- real/complex (type ki), the value of the S matrix element corresponding to the second external off-shell leg
    !  * m3_sq -- real/complex (type ki), the value of the internal mass squared
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
    !  An real (type ki) array of rank 1 and shape 2 corresponding to 
    !  the real/imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !*****
    function a3p2m_1mi_r(s23,s13,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(2) :: a3p2m_1mi_r
      !
      complex(ki) :: c_temp,c_temp_rat
      real(ki) :: sc13,sc23
      !
      a3p2m_1mi_r(:) = 0._ki
      !
      sc13=sign(un,s13+m3_sq)
      !
      sc23=sign(un,s23+m3_sq)
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        c_temp=-1._ki/(s13-s23)*zdilog((s13+m3_sq)/s13,-1._ki)+zdilog((s2&
          &3+m3_sq)/s23,-1._ki)/(s13-s23)
        !
        c_temp_rat=czero
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        if (par3 == 1) then
          !
          c_temp=s23/(s13-s23)**2*zdilog((s13+m3_sq)/s13,-1._ki)-s23/(s13-s&
            &23)**2*zdilog((s23+m3_sq)/s23,-1._ki)-m3_sq*z_log(m3_sq,-1._ki)&
            &/(s13-s23)**2-s23/(s13-s23)**2*z_log(-s23,-1._ki)+(s23+m3_sq)/(&
            &s13-s23)**2*z_log(-s13,-1._ki)-1._ki/(s13-s23)+1._ki/s13*q(1,1.&
            &_ki+m3_sq/s13,-sc13)*(s23+m3_sq)*m3_sq/(s13-s23)**2
          !
          c_temp_rat=-1._ki/(s13-s23)
          !
        else if (par3 == 2) then
          !
          c_temp=-s13/(s13-s23)**2*zdilog((s13+m3_sq)/s13,-1._ki)+s13/(s13-&
            &s23)**2*zdilog((s23+m3_sq)/s23,-1._ki)+m3_sq*z_log(m3_sq,-1._ki&
            &)/(s13-s23)**2+s13/(s13-s23)**2*z_log(-s23,-1._ki)-(s13+m3_sq)/&
            &(s13-s23)**2*z_log(-s13,-1._ki)+1._ki/s23*q(1,1._ki+m3_sq/s23,-&
            &sc23)*m3_sq/(s13-s23)-m3_sq*(s13+m3_sq)/s13/(s13-s23)**2*q(1,1.&
            &_ki+m3_sq/s13,-sc13)+1._ki/(s13-s23)
          !
          c_temp_rat=1._ki/(s13-s23)
          !
        else if (par3 == 3) then
          !
          c_temp=-1._ki/(s13-s23)*z_log(-s23,-1._ki)+1._ki/s13*q(1,1._ki+m3&
            &_sq/s13,-sc13)*m3_sq/(s13-s23)-1._ki/s23*q(1,1._ki+m3_sq/s23,-s&
            &c23)*m3_sq/(s13-s23)+1._ki/(s13-s23)*z_log(-s13,-1._ki)
          !
          c_temp_rat=0._ki
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function a3p2m_1mi_r:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "Unimplemented combination of Feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      ! cas avec deux parametres de feynman au numerateur
      else if ( (par1==0) ) then
        !
        if ( (par2 == 1) .and. (par3 == 1) ) then
          !
          c_temp=-s23**2/(s13-s23)**3*zdilog((s13+m3_sq)/s13,-1._ki)+s23**2&
            &/(s13-s23)**3*zdilog((s23+m3_sq)/s23,-1._ki)+1._ki/2._ki*m3_sq*&
            &(-m3_sq+2._ki*s23)/(s13-s23)**3*z_log(m3_sq,-1._ki)-1._ki/2._ki&
            &*(s23+m3_sq)*(3._ki*s23-m3_sq)/(s13-s23)**3*z_log(-s13,-1._ki)+&
            &3._ki/2._ki*s23**2/(s13-s23)**3*z_log(-s23,-1._ki)-1._ki/2._ki/&
            &s13**2*m3_sq**2/(s13-s23)**3*(s23+m3_sq)**2*q(2,1._ki+m3_sq/s13&
            &,-sc13)-m3_sq*(s23+m3_sq)*(-m3_sq+s23)/s13/(s13-s23)**3*q(1,1._&
            &ki+m3_sq/s13,-sc13)-1._ki/4._ki*(3._ki*s13**3-12._ki*s13**2*s23&
            &-2._ki*m3_sq*s13**2+9._ki*s13*s23**2+4._ki*m3_sq*s13*s23+2._ki*&
            &s13*m3_sq**2-2._ki*m3_sq*s23**2-4._ki*s23*m3_sq**2-2._ki*m3_sq*&
            &*3)/s13/(s13-s23)**3
          !
          c_temp_rat=-1._ki/4._ki*(3._ki*s13**2+m3_sq*s13-9._ki*s23*s13-7._&
            &ki*m3_sq*s23)/(s13+m3_sq)/(s13-s23)**2
          !
        else if ( (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp=-s13**2/(s13-s23)**3*zdilog((s13+m3_sq)/s13,-1._ki)+s13**2&
            &/(s13-s23)**3*zdilog((s23+m3_sq)/s23,-1._ki)+1._ki/2._ki*m3_sq*&
            &(4._ki*s13-2._ki*s23-m3_sq)/(s13-s23)**3*z_log(m3_sq,-1._ki)-1.&
            &_ki/2._ki*(s13+m3_sq)*(3._ki*s13-m3_sq)/(s13-s23)**3*z_log(-s13&
            &,-1._ki)+1._ki/2._ki*(2._ki*s23*m3_sq+3._ki*s13**2-2._ki*m3_sq*&
            &s13)/(s13-s23)**3*z_log(-s23,-1._ki)-1._ki/2._ki*m3_sq**2*(s13+&
            &m3_sq)**2/s13**2/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)+m3_sq*&
            &(s13-s23-m3_sq)/s23/(s13-s23)**2*q(1,1._ki+m3_sq/s23,-sc23)+1._&
            &ki/2._ki/s23**2/(s13-s23)*m3_sq**2*q(2,1._ki+m3_sq/s23,-sc23)-m&
            &3_sq*(s13-m3_sq)*(s13+m3_sq)/s13/(s13-s23)**3*q(1,1._ki+m3_sq/s&
            &13,-sc13)+1._ki/4._ki/s13/s23*(3._ki*s13*s23**3-2._ki*m3_sq*s13&
            &*s23**2-12._ki*s13**2*s23**2+9._ki*s13**3*s23+4._ki*m3_sq*s13**&
            &2*s23+2._ki*m3_sq**3*s23+2._ki*m3_sq**2*s13*s23-2._ki*s13**3*m3&
            &_sq)/(s13-s23)**3
          !
          c_temp_rat=1._ki/4._ki*(7._ki*m3_sq*s13+9._ki*s23*s13-m3_sq*s23-3&
            &._ki*s23**2)/(s13-s23)**2/(s23+m3_sq)
          !
        else if ( (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/2._ki/(s13-s23)*z_log(-s13,-1._ki)-1._ki/2._ki/(s13-&
            &s23)*z_log(-s23,-1._ki)-1._ki/2._ki/s13**2*m3_sq**2/(s13-s23)*q&
            &(2,1._ki+m3_sq/s13,-sc13)-1._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)*&
            &m3_sq/(s13-s23)+1._ki/2._ki/s23**2/(s13-s23)*m3_sq**2*q(2,1._ki&
            &+m3_sq/s23,-sc23)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m3_sq/(s&
            &13-s23)-1._ki/2._ki*m3_sq/s13/s23
          !
          c_temp_rat=-1._ki/2._ki*m3_sq/(s13+m3_sq)/(s23+m3_sq)
          !
        else if ( (par2 == 1) .and. (par3 == 2) ) then
          !
          c_temp=s23*s13/(s13-s23)**3*zdilog((s13+m3_sq)/s13,-1._ki)-s23*s1&
            &3/(s13-s23)**3*zdilog((s23+m3_sq)/s23,-1._ki)-1._ki/2._ki*m3_sq&
            &*(-m3_sq+2._ki*s13)/(s13-s23)**3*z_log(m3_sq,-1._ki)+1._ki/2._k&
            &i*(s23*m3_sq-m3_sq**2+m3_sq*s13+3._ki*s23*s13)/(s13-s23)**3*z_l&
            &og(-s13,-1._ki)-1._ki/2._ki*(s23*m3_sq-m3_sq*s13+3._ki*s23*s13)&
            &/(s13-s23)**3*z_log(-s23,-1._ki)+1._ki/2._ki*(s23+m3_sq)*m3_sq*&
            &*2*(s13+m3_sq)/s13**2/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)+1&
            &._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)**2*m3_sq**2&
            &+m3_sq*(-m3_sq**2+s23*s13)/s13/(s13-s23)**3*q(1,1._ki+m3_sq/s13&
            &,-sc13)-1._ki/4._ki*(3._ki*s13**3-3._ki*s13*s23**2+2._ki*s23*m3&
            &_sq**2+2._ki*m3_sq**3)/s13/(s13-s23)**3
          !
          c_temp_rat=-3._ki/4._ki*(s13+s23)/(s13-s23)**2
          !
        else if ( (par2 == 1) .and. (par3 == 3) ) then
          !
          c_temp=-1._ki/2._ki*(-m3_sq+s23)/(s13-s23)**2*z_log(-s13,-1._ki)+&
            &1._ki/2._ki*(-m3_sq+s23)/(s13-s23)**2*z_log(-s23,-1._ki)-1._ki/&
            &2._ki/s13**2*(s23+m3_sq)/(s13-s23)**2*m3_sq**2*q(2,1._ki+m3_sq/&
            &s13,-sc13)-1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)&
            &**2*m3_sq**2+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)/(s13-s23)**2*&
            &m3_sq**2+1._ki/2._ki/s13*(s13**2-s23*s13-m3_sq*s13+s23*m3_sq+m3&
            &_sq**2)/(s13-s23)**2
          !
          c_temp_rat=1._ki/2._ki*s13/(s13-s23)/(s13+m3_sq)
          !
        else if ( (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/2._ki*(s13-m3_sq)/(s13-s23)**2*z_log(-s13,-1._ki)-1.&
            &_ki/2._ki*(s13-m3_sq)/(s13-s23)**2*z_log(-s23,-1._ki)+1._ki/2._&
            &ki*m3_sq**2*(s13+m3_sq)/s13**2/(s13-s23)**2*q(2,1._ki+m3_sq/s13&
            &,-sc13)+1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)**2&
            &*m3_sq**2-1._ki/2._ki/s23**2/(s13-s23)*m3_sq**2*q(2,1._ki+m3_sq&
            &/s23,-sc23)-1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)/(s13-s23)**2*m&
            &3_sq**2-1._ki/2._ki/s13/s23*(-s13*s23**2+s13**2*s23+s23*m3_sq**&
            &2+m3_sq*s13*s23-m3_sq*s13**2)/(s13-s23)**2
          !
          c_temp_rat=-1._ki/2._ki*s23/(s13-s23)/(s23+m3_sq)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function a3p2m_1mi_r:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "Unimplemented combination of Feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      ! cas avec trois parametres de feynman au numerateur
      else
        !
        if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 1) ) then
          !
          c_temp=s23**3/(s13-s23)**4*zdilog((s13+m3_sq)/s13,-1._ki)-s23**3/&
            &(s13-s23)**4*zdilog((s23+m3_sq)/s23,-1._ki)-1._ki/6._ki*m3_sq*(&
            &-3._ki*s23*m3_sq+2._ki*m3_sq**2+6._ki*s23**2)/(s13-s23)**4*z_lo&
            &g(m3_sq,-1._ki)+1._ki/6._ki*(s23+m3_sq)*(2._ki*m3_sq**2-5._ki*s&
            &23*m3_sq+11._ki*s23**2)/(s13-s23)**4*z_log(-s13,-1._ki)-11._ki/&
            &6._ki*s23**3/(s13-s23)**4*z_log(-s23,-1._ki)+m3_sq*(s23+m3_sq)*&
            &(m3_sq**2+s23**2-s23*m3_sq)/s13/(s13-s23)**4*q(1,1._ki+m3_sq/s1&
            &3,-sc13)+1._ki/3._ki/s13**3*(s23+m3_sq)**3/(s13-s23)**4*m3_sq**&
            &3*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/2._ki*m3_sq**2*(s23+m3_sq)**&
            &2*(-2._ki*m3_sq+s23)/s13**2/(s13-s23)**4*q(2,1._ki+m3_sq/s13,-s&
            &c13)-1._ki/36._ki*(22._ki*s13**5-6._ki*m3_sq*s13**4-99._ki*s23*&
            &s13**4-6._ki*m3_sq**2*s13**3+198._ki*s23**2*s13**3+36._ki*s23*m&
            &3_sq*s13**3-54._ki*m3_sq*s23**2*s13**2-121._ki*s23**3*s13**2+18&
            &._ki*m3_sq**3*s13**2+18._ki*m3_sq**2*s23**2*s13+24._ki*s23**3*m&
            &3_sq*s13-36._ki*s23*m3_sq**3*s13-30._ki*m3_sq**4*s13+6._ki*s23*&
            &*3*m3_sq**2+18._ki*s23*m3_sq**4+6._ki*m3_sq**5+18._ki*m3_sq**3*&
            &s23**2)/s13**2/(s13-s23)**4
          !
          c_temp_rat=-1._ki/36._ki*(22._ki*s13**4+38._ki*m3_sq*s13**3-77._k&
            &i*s23*s13**3-124._ki*s23*m3_sq*s13**2+121._ki*s23**2*s13**2+4._&
            &ki*m3_sq**2*s13**2-23._ki*s23*m3_sq**2*s13+218._ki*s23**2*s13*m&
            &3_sq+85._ki*s23**2*m3_sq**2)/(s13-s23)**3/(s13+m3_sq)**2
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp=-s13**3/(s13-s23)**4*zdilog((s13+m3_sq)/s13,-1._ki)+s13**3&
            &/(s13-s23)**4*zdilog((s23+m3_sq)/s23,-1._ki)+1._ki/6._ki*m3_sq*&
            &(2._ki*m3_sq**2+6._ki*s23*m3_sq+6._ki*s23**2+18._ki*s13**2-18._&
            &ki*s23*s13-9._ki*m3_sq*s13)/(s13-s23)**4*z_log(m3_sq,-1._ki)-1.&
            &_ki/6._ki*(s13+m3_sq)*(2._ki*m3_sq**2-5._ki*m3_sq*s13+11._ki*s1&
            &3**2)/(s13-s23)**4*z_log(-s13,-1._ki)+1._ki/6._ki*(18._ki*s23*m&
            &3_sq*s13+11._ki*s13**3-12._ki*m3_sq*s13**2+6._ki*m3_sq**2*s13-6&
            &._ki*m3_sq**2*s23-6._ki*s23**2*m3_sq)/(s13-s23)**4*z_log(-s23,-&
            &1._ki)+1._ki/3._ki/s23**3*m3_sq**3/(s13-s23)*q(3,1._ki+m3_sq/s2&
            &3,-sc23)-m3_sq*(s13+m3_sq)*(s13**2-m3_sq*s13+m3_sq**2)/s13/(s13&
            &-s23)**4*q(1,1._ki+m3_sq/s13,-sc13)-1._ki/3._ki*m3_sq**3*(s13+m&
            &3_sq)**3/s13**3/(s13-s23)**4*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/2&
            &._ki*m3_sq*(2._ki*m3_sq**2+2._ki*s13**2-4._ki*s23*s13+2._ki*s23&
            &**2-3._ki*m3_sq*s13+3._ki*s23*m3_sq)/s23/(s13-s23)**3*q(1,1._ki&
            &+m3_sq/s23,-sc23)-1._ki/2._ki*m3_sq**2*(s13-2._ki*m3_sq)*(s13+m&
            &3_sq)**2/s13**2/(s13-s23)**4*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/2&
            &._ki*m3_sq**2*(-m3_sq+s13-s23)/s23**2/(s13-s23)**2*q(2,1._ki+m3&
            &_sq/s23,-sc23)+1._ki/36._ki/s13**2/s23**2*(-6._ki*m3_sq**2*s13*&
            &*5+121._ki*s23**2*s13**5-24._ki*s23*m3_sq*s13**5+54._ki*s23**2*&
            &s13**4*m3_sq+36._ki*s23*m3_sq**2*s13**4-198._ki*s23**3*s13**4-3&
            &6._ki*s23**3*m3_sq*s13**3+99._ki*s23**4*s13**3-36._ki*s23**2*m3&
            &_sq**2*s13**3+24._ki*s23**3*m3_sq**2*s13**2+6._ki*s23**4*m3_sq*&
            &s13**2-22._ki*s23**5*s13**2-12._ki*s23**2*m3_sq**4*s13+6._ki*s2&
            &3**2*m3_sq**5)/(s13-s23)**4
          !
          c_temp_rat=1._ki/36._ki*(218._ki*s23*m3_sq*s13**2+121._ki*s23**2*&
            &s13**2+85._ki*m3_sq**2*s13**2-124._ki*s23**2*s13*m3_sq-23._ki*s&
            &23*m3_sq**2*s13-77._ki*s23**3*s13+38._ki*s23**3*m3_sq+4._ki*s23&
            &**2*m3_sq**2+22._ki*s23**4)/(s23+m3_sq)**2/(s13-s23)**3
          !
        else if ( (par1 == 3) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/3._ki/(s13-s23)*z_log(-s13,-1._ki)-1._ki/3._ki/(s13-&
            &s23)*z_log(-s23,-1._ki)-1._ki/3._ki/s23**3*m3_sq**3/(s13-s23)*q&
            &(3,1._ki+m3_sq/s23,-sc23)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*&
            &m3_sq/(s13-s23)+1._ki/3._ki/s13**3/(s13-s23)*m3_sq**3*q(3,1._ki&
            &+m3_sq/s13,-sc13)-1._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)*m3_sq/(s&
            &13-s23)-1._ki/s13**2/(s13-s23)*m3_sq**2*q(2,1._ki+m3_sq/s13,-sc&
            &13)+1._ki/s23**2/(s13-s23)*m3_sq**2*q(2,1._ki+m3_sq/s23,-sc23)-&
            &1._ki/6._ki*(5._ki*s23*s13-s23*m3_sq-m3_sq*s13)/s23**2*m3_sq/s1&
            &3**2
          !
          c_temp_rat=-1._ki/6._ki*(5._ki*s23*s13+3._ki*m3_sq*s13+m3_sq**2+3&
            &._ki*s23*m3_sq)*m3_sq/(s13+m3_sq)**2/(s23+m3_sq)**2
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 2) ) then
          !
          c_temp=-s23**2*s13/(s13-s23)**4*zdilog((s13+m3_sq)/s13,-1._ki)+s2&
            &3**2*s13/(s13-s23)**4*zdilog((s23+m3_sq)/s23,-1._ki)+1._ki/6._k&
            &i*m3_sq*(2._ki*m3_sq**2+6._ki*s23*s13-3._ki*m3_sq*s13)/(s13-s23&
            &)**4*z_log(m3_sq,-1._ki)-1._ki/6._ki*(2._ki*m3_sq**3+2._ki*s23*&
            &*2*m3_sq-2._ki*m3_sq**2*s23+4._ki*s23*m3_sq*s13-m3_sq**2*s13+11&
            &._ki*s23**2*s13)/(s13-s23)**4*z_log(-s13,-1._ki)+1._ki/6._ki*(2&
            &._ki*s23**2*m3_sq-2._ki*m3_sq**2*s23-2._ki*s23*m3_sq*s13+2._ki*&
            &m3_sq**2*s13+11._ki*s23**2*s13)/(s13-s23)**4*z_log(-s23,-1._ki)&
            &-m3_sq*(m3_sq**3+s23**2*s13)/s13/(s13-s23)**4*q(1,1._ki+m3_sq/s&
            &13,-sc13)-1._ki/3._ki*m3_sq**3*(s23+m3_sq)**2*(s13+m3_sq)/s13**&
            &3/(s13-s23)**4*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/3._ki/s23*q(1,1&
            &._ki+m3_sq/s23,-sc23)*m3_sq**3/(s13-s23)**3-1._ki/2._ki*m3_sq**&
            &2*(s23+m3_sq)*(-2._ki*m3_sq**2-m3_sq*s13+s23*s13)/s13**2/(s13-s&
            &23)**4*q(2,1._ki+m3_sq/s13,-sc13)-1._ki/36._ki*(11._ki*s13**5-6&
            &6._ki*s23*s13**4-6._ki*m3_sq*s13**4+12._ki*s23*m3_sq*s13**3+33.&
            &_ki*s23**2*s13**3+6._ki*m3_sq**2*s13**3-6._ki*m3_sq*s23**2*s13*&
            &*2+22._ki*s23**3*s13**2-12._ki*s23*m3_sq**2*s13**2-6._ki*m3_sq*&
            &*3*s13**2-12._ki*m3_sq**2*s23**2*s13+24._ki*m3_sq**4*s13+12._ki&
            &*s23*m3_sq**3*s13-6._ki*m3_sq**3*s23**2-6._ki*m3_sq**5-12._ki*s&
            &23*m3_sq**4)/s13**2/(s13-s23)**4
          !
          c_temp_rat=-1._ki/36._ki*(11._ki*s13**3-55._ki*s23*s13**2+5._ki*m&
            &3_sq*s13**2-49._ki*s23*m3_sq*s13-22._ki*s23**2*s13-22._ki*s23**&
            &2*m3_sq)/(s13-s23)**3/(s13+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp=s23*s13**2/(s13-s23)**4*zdilog((s13+m3_sq)/s13,-1._ki)-s23&
            &*s13**2/(s13-s23)**4*zdilog((s23+m3_sq)/s23,-1._ki)-1._ki/6._ki&
            &*m3_sq*(2._ki*m3_sq**2+6._ki*s13**2+3._ki*s23*m3_sq-6._ki*m3_sq&
            &*s13)/(s13-s23)**4*z_log(m3_sq,-1._ki)+1._ki/6._ki*(4._ki*s23*m&
            &3_sq*s13+11._ki*s23*s13**2+2._ki*m3_sq**3+2._ki*m3_sq*s13**2-2.&
            &_ki*m3_sq**2*s13-m3_sq**2*s23)/(s13-s23)**4*z_log(-s13,-1._ki)-&
            &1._ki/6._ki*(4._ki*s23*m3_sq*s13+11._ki*s23*s13**2-4._ki*m3_sq*&
            &s13**2+4._ki*m3_sq**2*s13-4._ki*m3_sq**2*s23)/(s13-s23)**4*z_lo&
            &g(-s23,-1._ki)+m3_sq*(s23*s13**2+m3_sq**3)/s13/(s13-s23)**4*q(1&
            &,1._ki+m3_sq/s13,-sc13)+1._ki/3._ki*(s23+m3_sq)*m3_sq**3*(s13+m&
            &3_sq)**2/s13**3/(s13-s23)**4*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/6&
            &._ki*m3_sq**2*(3._ki*s13-3._ki*s23-4._ki*m3_sq)/s23/(s13-s23)**&
            &3*q(1,1._ki+m3_sq/s23,-sc23)+1._ki/2._ki*m3_sq**2*(s13+m3_sq)*(&
            &s23*s13-2._ki*m3_sq**2-s23*m3_sq)/s13**2/(s13-s23)**4*q(2,1._ki&
            &+m3_sq/s13,-sc13)+1._ki/6._ki/s23**2*m3_sq**3/(s13-s23)**2*q(2,&
            &1._ki+m3_sq/s23,-sc23)-1._ki/36._ki/s13**2/s23*(11._ki*s23**4*s&
            &13**2-66._ki*s23**3*s13**3-6._ki*s23**3*s13**2*m3_sq+33._ki*s23&
            &**2*s13**4+6._ki*s23**2*m3_sq**4+12._ki*s23**2*s13**3*m3_sq+24.&
            &_ki*s23**2*s13**2*m3_sq**2+22._ki*s13**5*s23-6._ki*s23*m3_sq*s1&
            &3**4-18._ki*s23*m3_sq**4*s13+6._ki*s23*m3_sq**5-12._ki*s23*m3_s&
            &q**2*s13**3+6._ki*m3_sq**2*s13**4)/(s13-s23)**4
          !
          c_temp_rat=-1._ki/36._ki*(22._ki*s23*s13**2+22._ki*m3_sq*s13**2+4&
            &9._ki*s23*m3_sq*s13+55._ki*s23**2*s13-5._ki*s23**2*m3_sq-11._ki&
            &*s23**3)/(s13-s23)**3/(s23+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/3._ki*(m3_sq**2+s23**2-s23*m3_sq)/(s13-s23)**3*z_log&
            &(-s13,-1._ki)-1._ki/3._ki*(m3_sq**2+s23**2-s23*m3_sq)/(s13-s23)&
            &**3*z_log(-s23,-1._ki)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m3_&
            &sq**3/(s13-s23)**3+1._ki/3._ki/s13**3*(s23+m3_sq)**2*m3_sq**3/(&
            &s13-s23)**3*q(3,1._ki+m3_sq/s13,-sc13)-1._ki/3._ki/s23*q(1,1._k&
            &i+m3_sq/s23,-sc23)*m3_sq**3/(s13-s23)**3-1._ki/s13**2*(s23+m3_s&
            &q)*m3_sq**3/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/6._ki&
            &/s13**2*(s13**4+m3_sq*s13**3-4._ki*s23*s13**3+3._ki*s23**2*s13*&
            &*2-3._ki*m3_sq**2*s13**2-s23**2*s13*m3_sq+4._ki*s23*m3_sq**2*s1&
            &3+5._ki*m3_sq**3*s13-m3_sq**4-2._ki*s23*m3_sq**3-s23**2*m3_sq**&
            &2)/(s13-s23)**3
          !
          c_temp_rat=1._ki/6._ki*s13*(s13**2+3._ki*m3_sq*s13-3._ki*s23*s13-&
            &5._ki*s23*m3_sq)/(s13-s23)**2/(s13+m3_sq)**2
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/3._ki*(s13**2-m3_sq*s13+m3_sq**2)/(s13-s23)**3*z_log&
            &(-s13,-1._ki)-1._ki/3._ki*(s13**2-m3_sq*s13+m3_sq**2)/(s13-s23)&
            &**3*z_log(-s23,-1._ki)-1._ki/3._ki/s23**3*m3_sq**3/(s13-s23)*q(&
            &3,1._ki+m3_sq/s23,-sc23)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m&
            &3_sq**3/(s13-s23)**3+1._ki/3._ki*m3_sq**3*(s13+m3_sq)**2/s13**3&
            &/(s13-s23)**3*q(3,1._ki+m3_sq/s13,-sc13)-1._ki/3._ki/s23*q(1,1.&
            &_ki+m3_sq/s23,-sc23)*m3_sq**3/(s13-s23)**3-m3_sq**3*(s13+m3_sq)&
            &/s13**2/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/3._ki/s23&
            &**2*m3_sq**3/(s13-s23)**2*q(2,1._ki+m3_sq/s23,-sc23)-1._ki/6._k&
            &i/s13**2/s23**2*(-s23*m3_sq*s13**4-m3_sq**2*s13**4+3._ki*s23**2&
            &*s13**4-4._ki*s23**3*s13**3+4._ki*s23*m3_sq**2*s13**3+s23**3*s1&
            &3**2*m3_sq+s23**4*s13**2-3._ki*s23**2*s13**2*m3_sq**2-3._ki*m3_&
            &sq**3*s23**2*s13+s23**2*m3_sq**4)/(s13-s23)**3
          !
          c_temp_rat=-1._ki/6._ki*s23*(3._ki*s23*s13+5._ki*m3_sq*s13-s23**2&
            &-3._ki*s23*m3_sq)/(s13-s23)**2/(s23+m3_sq)**2
          !
        else if ( (par1 == 1) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp=-1._ki/6._ki*(-2._ki*m3_sq+s23)/(s13-s23)**2*z_log(-s13,-1&
            &._ki)+1._ki/6._ki*(-2._ki*m3_sq+s23)/(s13-s23)**2*z_log(-s23,-1&
            &._ki)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)/(s13-s23)**2*m3_sq**&
            &2+1._ki/3._ki/s13**3*(s23+m3_sq)/(s13-s23)**2*m3_sq**3*q(3,1._k&
            &i+m3_sq/s13,-sc13)-1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(&
            &s13-s23)**2*m3_sq**2-1._ki/2._ki*m3_sq**2*(s23+2._ki*m3_sq)/s13&
            &**2/(s13-s23)**2*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/6._ki/s23**2*&
            &m3_sq**3/(s13-s23)**2*q(2,1._ki+m3_sq/s23,-sc23)+1._ki/6._ki/s1&
            &3**2/s23*(-s23**2*m3_sq**2-s23**2*s13**2+2._ki*s23**2*s13*m3_sq&
            &+s23*s13**3+5._ki*s23*m3_sq**2*s13-s23*m3_sq**3-2._ki*s23*m3_sq&
            &*s13**2-m3_sq**2*s13**2)/(s13-s23)**2
          !
          c_temp_rat=1._ki/6._ki*(s23*s13**2+m3_sq*s13**2-m3_sq**2*s13+m3_s&
            &q**2*s23)/(s13-s23)/(s13+m3_sq)**2/(s23+m3_sq)
          !
        else if ( (par1 == 2) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/6._ki*(s13-2._ki*m3_sq)/(s13-s23)**2*z_log(-s13,-1._&
            &ki)-1._ki/6._ki*(s13-2._ki*m3_sq)/(s13-s23)**2*z_log(-s23,-1._k&
            &i)+1._ki/3._ki/s23**3*m3_sq**3/(s13-s23)*q(3,1._ki+m3_sq/s23,-s&
            &c23)-1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)/(s13-s23)**2*m3_sq**2&
            &-1._ki/3._ki*m3_sq**3*(s13+m3_sq)/s13**3/(s13-s23)**2*q(3,1._ki&
            &+m3_sq/s13,-sc13)+1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s&
            &13-s23)**2*m3_sq**2+1._ki/2._ki*m3_sq**2*(s13+2._ki*m3_sq)/s13*&
            &*2/(s13-s23)**2*q(2,1._ki+m3_sq/s13,-sc13)-1._ki/6._ki*m3_sq**2&
            &*(m3_sq+3._ki*s13-3._ki*s23)/s23**2/(s13-s23)**2*q(2,1._ki+m3_s&
            &q/s23,-sc23)-1._ki/6._ki/s13**2/s23**2*(-s23**3*s13**2+s23**2*s&
            &13**3+4._ki*m3_sq**2*s23**2*s13-m3_sq**3*s23**2+2._ki*m3_sq*s23&
            &**2*s13**2-2._ki*s23*m3_sq*s13**3-2._ki*s23*m3_sq**2*s13**2+m3_&
            &sq**2*s13**3)/(s13-s23)**2
          !
          c_temp_rat=-1._ki/6._ki*(s23**2*s13+m3_sq**2*s13-m3_sq**2*s23+s23&
            &**2*m3_sq)/(s23+m3_sq)**2/(s13-s23)/(s13+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp=-1._ki/6._ki*(2._ki*s23*s13-m3_sq*s13-s23*m3_sq+2._ki*m3_s&
            &q**2)/(s13-s23)**3*z_log(-s13,-1._ki)+1._ki/6._ki*(2._ki*s23*s1&
            &3-m3_sq*s13-s23*m3_sq+2._ki*m3_sq**2)/(s13-s23)**3*z_log(-s23,-&
            &1._ki)-1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m3_sq**3/(s13-s23)*&
            &*3-1._ki/3._ki*(s13+m3_sq)*m3_sq**3*(s23+m3_sq)/s13**3/(s13-s23&
            &)**3*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/3._ki/s23*q(1,1._ki+m3_sq&
            &/s23,-sc23)*m3_sq**3/(s13-s23)**3+1._ki/2._ki*m3_sq**3*(s13+s23&
            &+2._ki*m3_sq)/s13**2/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)-1.&
            &_ki/6._ki/s23**2*m3_sq**3/(s13-s23)**2*q(2,1._ki+m3_sq/s23,-sc2&
            &3)+1._ki/6._ki/s23/s13**2*(-s23**3*s13**2-m3_sq**2*s23**2*s13+2&
            &._ki*m3_sq*s23**2*s13**2+m3_sq**3*s23**2-4._ki*s23*m3_sq**3*s13&
            &-2._ki*s23*m3_sq*s13**3+s23*s13**4+s23*m3_sq**4+m3_sq**2*s13**3&
            &)/(s13-s23)**3
          !
          c_temp_rat=1._ki/6._ki*(s23*s13**2+m3_sq*s13**2+s23**2*s13+s23**2&
            &*m3_sq)/(s13-s23)**2/(s13+m3_sq)/(s23+m3_sq)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function a3p2m_1mi_r:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "Unimplemented combination of Feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      end if
      !
      if ( (rat_or_tot_par%tot_selected)  ) then
        !
        a3p2m_1mi_r=(/real(c_temp,ki),aimag(c_temp)/)
        !
      else !if ( (rat_or_tot_par%rat_selected)  ) then
        !
        a3p2m_1mi_r=(/real(c_temp_rat,ki),aimag(c_temp_rat)/)
        !
      end if
      !
    end function a3p2m_1mi_r
    !
    function a3p2m_1mi_c(s23,s13,m3_sq,par1,par2,par3)
      !
      complex(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(2) :: a3p2m_1mi_c
      !
      complex(ki) :: c_temp,c_temp_rat
      real(ki) :: sc13,sc23
      !
      a3p2m_1mi_c(:) = 0._ki
      !
      sc13=sign(un,real(s13+m3_sq,ki))
      !
      sc23=sign(un,real(s23+m3_sq,ki))
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        c_temp=-1._ki/(s13-s23)*zdilog((s13+m3_sq)/s13,-1._ki)+zdilog((s2&
          &3+m3_sq)/s23,-1._ki)/(s13-s23)
        !
        c_temp_rat = czero
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        if (par3 == 1) then
          !
          c_temp=s23/(s13-s23)**2*zdilog((s13+m3_sq)/s13,-1._ki)-s23/(s13-s&
            &23)**2*zdilog((s23+m3_sq)/s23,-1._ki)-m3_sq*z_log(m3_sq,-1._ki)&
            &/(s13-s23)**2-s23/(s13-s23)**2*z_log(-s23,-1._ki)+(s23+m3_sq)/(&
            &s13-s23)**2*z_log(-s13,-1._ki)-1._ki/(s13-s23)+1._ki/s13*q(1,1.&
            &_ki+m3_sq/s13,-sc13)*(s23+m3_sq)*m3_sq/(s13-s23)**2
          !
          c_temp_rat=-1._ki/(s13-s23)
          !
        else if (par3 == 2) then
          !
          c_temp=-s13/(s13-s23)**2*zdilog((s13+m3_sq)/s13,-1._ki)+s13/(s13-&
            &s23)**2*zdilog((s23+m3_sq)/s23,-1._ki)+m3_sq*z_log(m3_sq,-1._ki&
            &)/(s13-s23)**2+s13/(s13-s23)**2*z_log(-s23,-1._ki)-(s13+m3_sq)/&
            &(s13-s23)**2*z_log(-s13,-1._ki)+1._ki/s23*q(1,1._ki+m3_sq/s23,-&
            &sc23)*m3_sq/(s13-s23)-m3_sq*(s13+m3_sq)/s13/(s13-s23)**2*q(1,1.&
            &_ki+m3_sq/s13,-sc13)+1._ki/(s13-s23)
          !
          c_temp_rat=1._ki/(s13-s23)
          !
        else if (par3 == 3) then
          !
          c_temp=-1._ki/(s13-s23)*z_log(-s23,-1._ki)+1._ki/s13*q(1,1._ki+m3&
            &_sq/s13,-sc13)*m3_sq/(s13-s23)-1._ki/s23*q(1,1._ki+m3_sq/s23,-s&
            &c23)*m3_sq/(s13-s23)+1._ki/(s13-s23)*z_log(-s13,-1._ki)
          !
          c_temp_rat=0._ki
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'in function a3p2m_1mi_c:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "unimplemented combination of feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      ! cas avec deux parametres de feynman au numerateur
      else if ( (par1==0) ) then
        !
        if ( (par2 == 1) .and. (par3 == 1) ) then
          !
          c_temp=-s23**2/(s13-s23)**3*zdilog((s13+m3_sq)/s13,-1._ki)+s23**2&
            &/(s13-s23)**3*zdilog((s23+m3_sq)/s23,-1._ki)+1._ki/2._ki*m3_sq*&
            &(-m3_sq+2._ki*s23)/(s13-s23)**3*z_log(m3_sq,-1._ki)-1._ki/2._ki&
            &*(s23+m3_sq)*(3._ki*s23-m3_sq)/(s13-s23)**3*z_log(-s13,-1._ki)+&
            &3._ki/2._ki*s23**2/(s13-s23)**3*z_log(-s23,-1._ki)-1._ki/2._ki/&
            &s13**2*m3_sq**2/(s13-s23)**3*(s23+m3_sq)**2*q(2,1._ki+m3_sq/s13&
            &,-sc13)-m3_sq*(s23+m3_sq)*(-m3_sq+s23)/s13/(s13-s23)**3*q(1,1._&
            &ki+m3_sq/s13,-sc13)-1._ki/4._ki*(3._ki*s13**3-12._ki*s13**2*s23&
            &-2._ki*m3_sq*s13**2+9._ki*s13*s23**2+4._ki*m3_sq*s13*s23+2._ki*&
            &s13*m3_sq**2-2._ki*m3_sq*s23**2-4._ki*s23*m3_sq**2-2._ki*m3_sq*&
            &*3)/s13/(s13-s23)**3
          !
          c_temp_rat=-1._ki/4._ki*(3._ki*s13**2+m3_sq*s13-9._ki*s23*s13-7._&
            &ki*m3_sq*s23)/(s13+m3_sq)/(s13-s23)**2
          !
        else if ( (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp=-s13**2/(s13-s23)**3*zdilog((s13+m3_sq)/s13,-1._ki)+s13**2&
            &/(s13-s23)**3*zdilog((s23+m3_sq)/s23,-1._ki)+1._ki/2._ki*m3_sq*&
            &(4._ki*s13-2._ki*s23-m3_sq)/(s13-s23)**3*z_log(m3_sq,-1._ki)-1.&
            &_ki/2._ki*(s13+m3_sq)*(3._ki*s13-m3_sq)/(s13-s23)**3*z_log(-s13&
            &,-1._ki)+1._ki/2._ki*(2._ki*s23*m3_sq+3._ki*s13**2-2._ki*m3_sq*&
            &s13)/(s13-s23)**3*z_log(-s23,-1._ki)-1._ki/2._ki*m3_sq**2*(s13+&
            &m3_sq)**2/s13**2/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)+m3_sq*&
            &(s13-s23-m3_sq)/s23/(s13-s23)**2*q(1,1._ki+m3_sq/s23,-sc23)+1._&
            &ki/2._ki/s23**2/(s13-s23)*m3_sq**2*q(2,1._ki+m3_sq/s23,-sc23)-m&
            &3_sq*(s13-m3_sq)*(s13+m3_sq)/s13/(s13-s23)**3*q(1,1._ki+m3_sq/s&
            &13,-sc13)+1._ki/4._ki/s13/s23*(3._ki*s13*s23**3-2._ki*m3_sq*s13&
            &*s23**2-12._ki*s13**2*s23**2+9._ki*s13**3*s23+4._ki*m3_sq*s13**&
            &2*s23+2._ki*m3_sq**3*s23+2._ki*m3_sq**2*s13*s23-2._ki*s13**3*m3&
            &_sq)/(s13-s23)**3
          !
          c_temp_rat=1._ki/4._ki*(7._ki*m3_sq*s13+9._ki*s23*s13-m3_sq*s23-3&
            &._ki*s23**2)/(s13-s23)**2/(s23+m3_sq)
          !
        else if ( (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/2._ki/(s13-s23)*z_log(-s13,-1._ki)-1._ki/2._ki/(s13-&
            &s23)*z_log(-s23,-1._ki)-1._ki/2._ki/s13**2*m3_sq**2/(s13-s23)*q&
            &(2,1._ki+m3_sq/s13,-sc13)-1._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)*&
            &m3_sq/(s13-s23)+1._ki/2._ki/s23**2/(s13-s23)*m3_sq**2*q(2,1._ki&
            &+m3_sq/s23,-sc23)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m3_sq/(s&
            &13-s23)-1._ki/2._ki*m3_sq/s13/s23
          !
          c_temp_rat=-1._ki/2._ki*m3_sq/(s13+m3_sq)/(s23+m3_sq)
          !
        else if ( (par2 == 1) .and. (par3 == 2) ) then
          !
          c_temp=s23*s13/(s13-s23)**3*zdilog((s13+m3_sq)/s13,-1._ki)-s23*s1&
            &3/(s13-s23)**3*zdilog((s23+m3_sq)/s23,-1._ki)-1._ki/2._ki*m3_sq&
            &*(-m3_sq+2._ki*s13)/(s13-s23)**3*z_log(m3_sq,-1._ki)+1._ki/2._k&
            &i*(s23*m3_sq-m3_sq**2+m3_sq*s13+3._ki*s23*s13)/(s13-s23)**3*z_l&
            &og(-s13,-1._ki)-1._ki/2._ki*(s23*m3_sq-m3_sq*s13+3._ki*s23*s13)&
            &/(s13-s23)**3*z_log(-s23,-1._ki)+1._ki/2._ki*(s23+m3_sq)*m3_sq*&
            &*2*(s13+m3_sq)/s13**2/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)+1&
            &._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)**2*m3_sq**2&
            &+m3_sq*(-m3_sq**2+s23*s13)/s13/(s13-s23)**3*q(1,1._ki+m3_sq/s13&
            &,-sc13)-1._ki/4._ki*(3._ki*s13**3-3._ki*s13*s23**2+2._ki*s23*m3&
            &_sq**2+2._ki*m3_sq**3)/s13/(s13-s23)**3
          !
          c_temp_rat=-3._ki/4._ki*(s13+s23)/(s13-s23)**2
          !
        else if ( (par2 == 1) .and. (par3 == 3) ) then
          !
          c_temp=-1._ki/2._ki*(-m3_sq+s23)/(s13-s23)**2*z_log(-s13,-1._ki)+&
            &1._ki/2._ki*(-m3_sq+s23)/(s13-s23)**2*z_log(-s23,-1._ki)-1._ki/&
            &2._ki/s13**2*(s23+m3_sq)/(s13-s23)**2*m3_sq**2*q(2,1._ki+m3_sq/&
            &s13,-sc13)-1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)&
            &**2*m3_sq**2+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)/(s13-s23)**2*&
            &m3_sq**2+1._ki/2._ki/s13*(s13**2-s23*s13-m3_sq*s13+s23*m3_sq+m3&
            &_sq**2)/(s13-s23)**2
          !
          c_temp_rat=1._ki/2._ki*s13/(s13-s23)/(s13+m3_sq)
          !
        else if ( (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/2._ki*(s13-m3_sq)/(s13-s23)**2*z_log(-s13,-1._ki)-1.&
            &_ki/2._ki*(s13-m3_sq)/(s13-s23)**2*z_log(-s23,-1._ki)+1._ki/2._&
            &ki*m3_sq**2*(s13+m3_sq)/s13**2/(s13-s23)**2*q(2,1._ki+m3_sq/s13&
            &,-sc13)+1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)**2&
            &*m3_sq**2-1._ki/2._ki/s23**2/(s13-s23)*m3_sq**2*q(2,1._ki+m3_sq&
            &/s23,-sc23)-1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)/(s13-s23)**2*m&
            &3_sq**2-1._ki/2._ki/s13/s23*(-s13*s23**2+s13**2*s23+s23*m3_sq**&
            &2+m3_sq*s13*s23-m3_sq*s13**2)/(s13-s23)**2
          !
          c_temp_rat=-1._ki/2._ki*s23/(s13-s23)/(s23+m3_sq)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'in function a3p2m_1mi_c:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "unimplemented combination of feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      ! cas avec trois parametres de feynman au numerateur
      else
        !
        if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 1) ) then
          !
          c_temp=s23**3/(s13-s23)**4*zdilog((s13+m3_sq)/s13,-1._ki)-s23**3/&
            &(s13-s23)**4*zdilog((s23+m3_sq)/s23,-1._ki)-1._ki/6._ki*m3_sq*(&
            &-3._ki*s23*m3_sq+2._ki*m3_sq**2+6._ki*s23**2)/(s13-s23)**4*z_lo&
            &g(m3_sq,-1._ki)+1._ki/6._ki*(s23+m3_sq)*(2._ki*m3_sq**2-5._ki*s&
            &23*m3_sq+11._ki*s23**2)/(s13-s23)**4*z_log(-s13,-1._ki)-11._ki/&
            &6._ki*s23**3/(s13-s23)**4*z_log(-s23,-1._ki)+m3_sq*(s23+m3_sq)*&
            &(m3_sq**2+s23**2-s23*m3_sq)/s13/(s13-s23)**4*q(1,1._ki+m3_sq/s1&
            &3,-sc13)+1._ki/3._ki/s13**3*(s23+m3_sq)**3/(s13-s23)**4*m3_sq**&
            &3*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/2._ki*m3_sq**2*(s23+m3_sq)**&
            &2*(-2._ki*m3_sq+s23)/s13**2/(s13-s23)**4*q(2,1._ki+m3_sq/s13,-s&
            &c13)-1._ki/36._ki*(22._ki*s13**5-6._ki*m3_sq*s13**4-99._ki*s23*&
            &s13**4-6._ki*m3_sq**2*s13**3+198._ki*s23**2*s13**3+36._ki*s23*m&
            &3_sq*s13**3-54._ki*m3_sq*s23**2*s13**2-121._ki*s23**3*s13**2+18&
            &._ki*m3_sq**3*s13**2+18._ki*m3_sq**2*s23**2*s13+24._ki*s23**3*m&
            &3_sq*s13-36._ki*s23*m3_sq**3*s13-30._ki*m3_sq**4*s13+6._ki*s23*&
            &*3*m3_sq**2+18._ki*s23*m3_sq**4+6._ki*m3_sq**5+18._ki*m3_sq**3*&
            &s23**2)/s13**2/(s13-s23)**4
          !
          c_temp_rat=-1._ki/36._ki*(22._ki*s13**4+38._ki*m3_sq*s13**3-77._k&
            &i*s23*s13**3-124._ki*s23*m3_sq*s13**2+121._ki*s23**2*s13**2+4._&
            &ki*m3_sq**2*s13**2-23._ki*s23*m3_sq**2*s13+218._ki*s23**2*s13*m&
            &3_sq+85._ki*s23**2*m3_sq**2)/(s13-s23)**3/(s13+m3_sq)**2
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp=-s13**3/(s13-s23)**4*zdilog((s13+m3_sq)/s13,-1._ki)+s13**3&
            &/(s13-s23)**4*zdilog((s23+m3_sq)/s23,-1._ki)+1._ki/6._ki*m3_sq*&
            &(2._ki*m3_sq**2+6._ki*s23*m3_sq+6._ki*s23**2+18._ki*s13**2-18._&
            &ki*s23*s13-9._ki*m3_sq*s13)/(s13-s23)**4*z_log(m3_sq,-1._ki)-1.&
            &_ki/6._ki*(s13+m3_sq)*(2._ki*m3_sq**2-5._ki*m3_sq*s13+11._ki*s1&
            &3**2)/(s13-s23)**4*z_log(-s13,-1._ki)+1._ki/6._ki*(18._ki*s23*m&
            &3_sq*s13+11._ki*s13**3-12._ki*m3_sq*s13**2+6._ki*m3_sq**2*s13-6&
            &._ki*m3_sq**2*s23-6._ki*s23**2*m3_sq)/(s13-s23)**4*z_log(-s23,-&
            &1._ki)+1._ki/3._ki/s23**3*m3_sq**3/(s13-s23)*q(3,1._ki+m3_sq/s2&
            &3,-sc23)-m3_sq*(s13+m3_sq)*(s13**2-m3_sq*s13+m3_sq**2)/s13/(s13&
            &-s23)**4*q(1,1._ki+m3_sq/s13,-sc13)-1._ki/3._ki*m3_sq**3*(s13+m&
            &3_sq)**3/s13**3/(s13-s23)**4*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/2&
            &._ki*m3_sq*(2._ki*m3_sq**2+2._ki*s13**2-4._ki*s23*s13+2._ki*s23&
            &**2-3._ki*m3_sq*s13+3._ki*s23*m3_sq)/s23/(s13-s23)**3*q(1,1._ki&
            &+m3_sq/s23,-sc23)-1._ki/2._ki*m3_sq**2*(s13-2._ki*m3_sq)*(s13+m&
            &3_sq)**2/s13**2/(s13-s23)**4*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/2&
            &._ki*m3_sq**2*(-m3_sq+s13-s23)/s23**2/(s13-s23)**2*q(2,1._ki+m3&
            &_sq/s23,-sc23)+1._ki/36._ki/s13**2/s23**2*(-6._ki*m3_sq**2*s13*&
            &*5+121._ki*s23**2*s13**5-24._ki*s23*m3_sq*s13**5+54._ki*s23**2*&
            &s13**4*m3_sq+36._ki*s23*m3_sq**2*s13**4-198._ki*s23**3*s13**4-3&
            &6._ki*s23**3*m3_sq*s13**3+99._ki*s23**4*s13**3-36._ki*s23**2*m3&
            &_sq**2*s13**3+24._ki*s23**3*m3_sq**2*s13**2+6._ki*s23**4*m3_sq*&
            &s13**2-22._ki*s23**5*s13**2-12._ki*s23**2*m3_sq**4*s13+6._ki*s2&
            &3**2*m3_sq**5)/(s13-s23)**4
          !
          c_temp_rat=1._ki/36._ki*(218._ki*s23*m3_sq*s13**2+121._ki*s23**2*&
            &s13**2+85._ki*m3_sq**2*s13**2-124._ki*s23**2*s13*m3_sq-23._ki*s&
            &23*m3_sq**2*s13-77._ki*s23**3*s13+38._ki*s23**3*m3_sq+4._ki*s23&
            &**2*m3_sq**2+22._ki*s23**4)/(s23+m3_sq)**2/(s13-s23)**3
          !
        else if ( (par1 == 3) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/3._ki/(s13-s23)*z_log(-s13,-1._ki)-1._ki/3._ki/(s13-&
            &s23)*z_log(-s23,-1._ki)-1._ki/3._ki/s23**3*m3_sq**3/(s13-s23)*q&
            &(3,1._ki+m3_sq/s23,-sc23)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*&
            &m3_sq/(s13-s23)+1._ki/3._ki/s13**3/(s13-s23)*m3_sq**3*q(3,1._ki&
            &+m3_sq/s13,-sc13)-1._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)*m3_sq/(s&
            &13-s23)-1._ki/s13**2/(s13-s23)*m3_sq**2*q(2,1._ki+m3_sq/s13,-sc&
            &13)+1._ki/s23**2/(s13-s23)*m3_sq**2*q(2,1._ki+m3_sq/s23,-sc23)-&
            &1._ki/6._ki*(5._ki*s23*s13-s23*m3_sq-m3_sq*s13)/s23**2*m3_sq/s1&
            &3**2
          !
          c_temp_rat=-1._ki/6._ki*(5._ki*s23*s13+3._ki*m3_sq*s13+m3_sq**2+3&
            &._ki*s23*m3_sq)*m3_sq/(s13+m3_sq)**2/(s23+m3_sq)**2
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 2) ) then
          !
          c_temp=-s23**2*s13/(s13-s23)**4*zdilog((s13+m3_sq)/s13,-1._ki)+s2&
            &3**2*s13/(s13-s23)**4*zdilog((s23+m3_sq)/s23,-1._ki)+1._ki/6._k&
            &i*m3_sq*(2._ki*m3_sq**2+6._ki*s23*s13-3._ki*m3_sq*s13)/(s13-s23&
            &)**4*z_log(m3_sq,-1._ki)-1._ki/6._ki*(2._ki*m3_sq**3+2._ki*s23*&
            &*2*m3_sq-2._ki*m3_sq**2*s23+4._ki*s23*m3_sq*s13-m3_sq**2*s13+11&
            &._ki*s23**2*s13)/(s13-s23)**4*z_log(-s13,-1._ki)+1._ki/6._ki*(2&
            &._ki*s23**2*m3_sq-2._ki*m3_sq**2*s23-2._ki*s23*m3_sq*s13+2._ki*&
            &m3_sq**2*s13+11._ki*s23**2*s13)/(s13-s23)**4*z_log(-s23,-1._ki)&
            &-m3_sq*(m3_sq**3+s23**2*s13)/s13/(s13-s23)**4*q(1,1._ki+m3_sq/s&
            &13,-sc13)-1._ki/3._ki*m3_sq**3*(s23+m3_sq)**2*(s13+m3_sq)/s13**&
            &3/(s13-s23)**4*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/3._ki/s23*q(1,1&
            &._ki+m3_sq/s23,-sc23)*m3_sq**3/(s13-s23)**3-1._ki/2._ki*m3_sq**&
            &2*(s23+m3_sq)*(-2._ki*m3_sq**2-m3_sq*s13+s23*s13)/s13**2/(s13-s&
            &23)**4*q(2,1._ki+m3_sq/s13,-sc13)-1._ki/36._ki*(11._ki*s13**5-6&
            &6._ki*s23*s13**4-6._ki*m3_sq*s13**4+12._ki*s23*m3_sq*s13**3+33.&
            &_ki*s23**2*s13**3+6._ki*m3_sq**2*s13**3-6._ki*m3_sq*s23**2*s13*&
            &*2+22._ki*s23**3*s13**2-12._ki*s23*m3_sq**2*s13**2-6._ki*m3_sq*&
            &*3*s13**2-12._ki*m3_sq**2*s23**2*s13+24._ki*m3_sq**4*s13+12._ki&
            &*s23*m3_sq**3*s13-6._ki*m3_sq**3*s23**2-6._ki*m3_sq**5-12._ki*s&
            &23*m3_sq**4)/s13**2/(s13-s23)**4
          !
          c_temp_rat=-1._ki/36._ki*(11._ki*s13**3-55._ki*s23*s13**2+5._ki*m&
            &3_sq*s13**2-49._ki*s23*m3_sq*s13-22._ki*s23**2*s13-22._ki*s23**&
            &2*m3_sq)/(s13-s23)**3/(s13+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 2) ) then
          !
          c_temp=s23*s13**2/(s13-s23)**4*zdilog((s13+m3_sq)/s13,-1._ki)-s23&
            &*s13**2/(s13-s23)**4*zdilog((s23+m3_sq)/s23,-1._ki)-1._ki/6._ki&
            &*m3_sq*(2._ki*m3_sq**2+6._ki*s13**2+3._ki*s23*m3_sq-6._ki*m3_sq&
            &*s13)/(s13-s23)**4*z_log(m3_sq,-1._ki)+1._ki/6._ki*(4._ki*s23*m&
            &3_sq*s13+11._ki*s23*s13**2+2._ki*m3_sq**3+2._ki*m3_sq*s13**2-2.&
            &_ki*m3_sq**2*s13-m3_sq**2*s23)/(s13-s23)**4*z_log(-s13,-1._ki)-&
            &1._ki/6._ki*(4._ki*s23*m3_sq*s13+11._ki*s23*s13**2-4._ki*m3_sq*&
            &s13**2+4._ki*m3_sq**2*s13-4._ki*m3_sq**2*s23)/(s13-s23)**4*z_lo&
            &g(-s23,-1._ki)+m3_sq*(s23*s13**2+m3_sq**3)/s13/(s13-s23)**4*q(1&
            &,1._ki+m3_sq/s13,-sc13)+1._ki/3._ki*(s23+m3_sq)*m3_sq**3*(s13+m&
            &3_sq)**2/s13**3/(s13-s23)**4*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/6&
            &._ki*m3_sq**2*(3._ki*s13-3._ki*s23-4._ki*m3_sq)/s23/(s13-s23)**&
            &3*q(1,1._ki+m3_sq/s23,-sc23)+1._ki/2._ki*m3_sq**2*(s13+m3_sq)*(&
            &s23*s13-2._ki*m3_sq**2-s23*m3_sq)/s13**2/(s13-s23)**4*q(2,1._ki&
            &+m3_sq/s13,-sc13)+1._ki/6._ki/s23**2*m3_sq**3/(s13-s23)**2*q(2,&
            &1._ki+m3_sq/s23,-sc23)-1._ki/36._ki/s13**2/s23*(11._ki*s23**4*s&
            &13**2-66._ki*s23**3*s13**3-6._ki*s23**3*s13**2*m3_sq+33._ki*s23&
            &**2*s13**4+6._ki*s23**2*m3_sq**4+12._ki*s23**2*s13**3*m3_sq+24.&
            &_ki*s23**2*s13**2*m3_sq**2+22._ki*s13**5*s23-6._ki*s23*m3_sq*s1&
            &3**4-18._ki*s23*m3_sq**4*s13+6._ki*s23*m3_sq**5-12._ki*s23*m3_s&
            &q**2*s13**3+6._ki*m3_sq**2*s13**4)/(s13-s23)**4
          !
          c_temp_rat=-1._ki/36._ki*(22._ki*s23*s13**2+22._ki*m3_sq*s13**2+4&
            &9._ki*s23*m3_sq*s13+55._ki*s23**2*s13-5._ki*s23**2*m3_sq-11._ki&
            &*s23**3)/(s13-s23)**3/(s23+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 1) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/3._ki*(m3_sq**2+s23**2-s23*m3_sq)/(s13-s23)**3*z_log&
            &(-s13,-1._ki)-1._ki/3._ki*(m3_sq**2+s23**2-s23*m3_sq)/(s13-s23)&
            &**3*z_log(-s23,-1._ki)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m3_&
            &sq**3/(s13-s23)**3+1._ki/3._ki/s13**3*(s23+m3_sq)**2*m3_sq**3/(&
            &s13-s23)**3*q(3,1._ki+m3_sq/s13,-sc13)-1._ki/3._ki/s23*q(1,1._k&
            &i+m3_sq/s23,-sc23)*m3_sq**3/(s13-s23)**3-1._ki/s13**2*(s23+m3_s&
            &q)*m3_sq**3/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/6._ki&
            &/s13**2*(s13**4+m3_sq*s13**3-4._ki*s23*s13**3+3._ki*s23**2*s13*&
            &*2-3._ki*m3_sq**2*s13**2-s23**2*s13*m3_sq+4._ki*s23*m3_sq**2*s1&
            &3+5._ki*m3_sq**3*s13-m3_sq**4-2._ki*s23*m3_sq**3-s23**2*m3_sq**&
            &2)/(s13-s23)**3
          !
          c_temp_rat=1._ki/6._ki*s13*(s13**2+3._ki*m3_sq*s13-3._ki*s23*s13-&
            &5._ki*s23*m3_sq)/(s13-s23)**2/(s13+m3_sq)**2
          !
        else if ( (par1 == 2) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/3._ki*(s13**2-m3_sq*s13+m3_sq**2)/(s13-s23)**3*z_log&
            &(-s13,-1._ki)-1._ki/3._ki*(s13**2-m3_sq*s13+m3_sq**2)/(s13-s23)&
            &**3*z_log(-s23,-1._ki)-1._ki/3._ki/s23**3*m3_sq**3/(s13-s23)*q(&
            &3,1._ki+m3_sq/s23,-sc23)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m&
            &3_sq**3/(s13-s23)**3+1._ki/3._ki*m3_sq**3*(s13+m3_sq)**2/s13**3&
            &/(s13-s23)**3*q(3,1._ki+m3_sq/s13,-sc13)-1._ki/3._ki/s23*q(1,1.&
            &_ki+m3_sq/s23,-sc23)*m3_sq**3/(s13-s23)**3-m3_sq**3*(s13+m3_sq)&
            &/s13**2/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/3._ki/s23&
            &**2*m3_sq**3/(s13-s23)**2*q(2,1._ki+m3_sq/s23,-sc23)-1._ki/6._k&
            &i/s13**2/s23**2*(-s23*m3_sq*s13**4-m3_sq**2*s13**4+3._ki*s23**2&
            &*s13**4-4._ki*s23**3*s13**3+4._ki*s23*m3_sq**2*s13**3+s23**3*s1&
            &3**2*m3_sq+s23**4*s13**2-3._ki*s23**2*s13**2*m3_sq**2-3._ki*m3_&
            &sq**3*s23**2*s13+s23**2*m3_sq**4)/(s13-s23)**3
          !
          c_temp_rat=-1._ki/6._ki*s23*(3._ki*s23*s13+5._ki*m3_sq*s13-s23**2&
            &-3._ki*s23*m3_sq)/(s13-s23)**2/(s23+m3_sq)**2
          !
        else if ( (par1 == 1) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp=-1._ki/6._ki*(-2._ki*m3_sq+s23)/(s13-s23)**2*z_log(-s13,-1&
            &._ki)+1._ki/6._ki*(-2._ki*m3_sq+s23)/(s13-s23)**2*z_log(-s23,-1&
            &._ki)+1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)/(s13-s23)**2*m3_sq**&
            &2+1._ki/3._ki/s13**3*(s23+m3_sq)/(s13-s23)**2*m3_sq**3*q(3,1._k&
            &i+m3_sq/s13,-sc13)-1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(&
            &s13-s23)**2*m3_sq**2-1._ki/2._ki*m3_sq**2*(s23+2._ki*m3_sq)/s13&
            &**2/(s13-s23)**2*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/6._ki/s23**2*&
            &m3_sq**3/(s13-s23)**2*q(2,1._ki+m3_sq/s23,-sc23)+1._ki/6._ki/s1&
            &3**2/s23*(-s23**2*m3_sq**2-s23**2*s13**2+2._ki*s23**2*s13*m3_sq&
            &+s23*s13**3+5._ki*s23*m3_sq**2*s13-s23*m3_sq**3-2._ki*s23*m3_sq&
            &*s13**2-m3_sq**2*s13**2)/(s13-s23)**2
          !
          c_temp_rat=1._ki/6._ki*(s23*s13**2+m3_sq*s13**2-m3_sq**2*s13+m3_s&
            &q**2*s23)/(s13-s23)/(s13+m3_sq)**2/(s23+m3_sq)
          !
        else if ( (par1 == 2) .and. (par2 == 3) .and. (par3 == 3) ) then
          !
          c_temp=1._ki/6._ki*(s13-2._ki*m3_sq)/(s13-s23)**2*z_log(-s13,-1._&
            &ki)-1._ki/6._ki*(s13-2._ki*m3_sq)/(s13-s23)**2*z_log(-s23,-1._k&
            &i)+1._ki/3._ki/s23**3*m3_sq**3/(s13-s23)*q(3,1._ki+m3_sq/s23,-s&
            &c23)-1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)/(s13-s23)**2*m3_sq**2&
            &-1._ki/3._ki*m3_sq**3*(s13+m3_sq)/s13**3/(s13-s23)**2*q(3,1._ki&
            &+m3_sq/s13,-sc13)+1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s&
            &13-s23)**2*m3_sq**2+1._ki/2._ki*m3_sq**2*(s13+2._ki*m3_sq)/s13*&
            &*2/(s13-s23)**2*q(2,1._ki+m3_sq/s13,-sc13)-1._ki/6._ki*m3_sq**2&
            &*(m3_sq+3._ki*s13-3._ki*s23)/s23**2/(s13-s23)**2*q(2,1._ki+m3_s&
            &q/s23,-sc23)-1._ki/6._ki/s13**2/s23**2*(-s23**3*s13**2+s23**2*s&
            &13**3+4._ki*m3_sq**2*s23**2*s13-m3_sq**3*s23**2+2._ki*m3_sq*s23&
            &**2*s13**2-2._ki*s23*m3_sq*s13**3-2._ki*s23*m3_sq**2*s13**2+m3_&
            &sq**2*s13**3)/(s13-s23)**2
          !
          c_temp_rat=-1._ki/6._ki*(s23**2*s13+m3_sq**2*s13-m3_sq**2*s23+s23&
            &**2*m3_sq)/(s23+m3_sq)**2/(s13-s23)/(s13+m3_sq)
          !
        else if ( (par1 == 1) .and. (par2 == 2) .and. (par3 == 3) ) then
          !
          c_temp=-1._ki/6._ki*(2._ki*s23*s13-m3_sq*s13-s23*m3_sq+2._ki*m3_s&
            &q**2)/(s13-s23)**3*z_log(-s13,-1._ki)+1._ki/6._ki*(2._ki*s23*s1&
            &3-m3_sq*s13-s23*m3_sq+2._ki*m3_sq**2)/(s13-s23)**3*z_log(-s23,-&
            &1._ki)-1._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m3_sq**3/(s13-s23)*&
            &*3-1._ki/3._ki*(s13+m3_sq)*m3_sq**3*(s23+m3_sq)/s13**3/(s13-s23&
            &)**3*q(3,1._ki+m3_sq/s13,-sc13)+1._ki/3._ki/s23*q(1,1._ki+m3_sq&
            &/s23,-sc23)*m3_sq**3/(s13-s23)**3+1._ki/2._ki*m3_sq**3*(s13+s23&
            &+2._ki*m3_sq)/s13**2/(s13-s23)**3*q(2,1._ki+m3_sq/s13,-sc13)-1.&
            &_ki/6._ki/s23**2*m3_sq**3/(s13-s23)**2*q(2,1._ki+m3_sq/s23,-sc2&
            &3)+1._ki/6._ki/s23/s13**2*(-s23**3*s13**2-m3_sq**2*s23**2*s13+2&
            &._ki*m3_sq*s23**2*s13**2+m3_sq**3*s23**2-4._ki*s23*m3_sq**3*s13&
            &-2._ki*s23*m3_sq*s13**3+s23*s13**4+s23*m3_sq**4+m3_sq**2*s13**3&
            &)/(s13-s23)**3
          !
          c_temp_rat=1._ki/6._ki*(s23*s13**2+m3_sq*s13**2+s23**2*s13+s23**2&
            &*m3_sq)/(s13-s23)**2/(s13+m3_sq)/(s23+m3_sq)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'in function a3p2m_1mi_c:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "unimplemented combination of feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      end if
      !
      if ( (rat_or_tot_par%tot_selected)  ) then
        !
        a3p2m_1mi_c=(/real(c_temp,ki),aimag(c_temp)/)
        !
      else !if ( (rat_or_tot_par%rat_selected)  ) then
        !
        a3p2m_1mi_c=(/real(c_temp_rat,ki),aimag(c_temp_rat)/)
        !
      end if
      !
    end function a3p2m_1mi_c
    !
    !
    !****if* src/integral/three_point/function_3p2m_1mi/a3p2m_1mi_np2
    ! NAME
    !
    !  Function a3p2m_1mi_np2
    !
    ! USAGE
    !
    !  real_dim2 = a3p2m_1mi_np2(s23,s13,m3_sq,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the non divergent two off-shell external leg three point function in n+2 dimensions. 
    !  with up to one Feynman parameter in the numerator.
    !  It retuns an array of 2 reals corresponding to the real/imaginary part of the
    !  constant term.
    !
    ! INPUTS
    !
    !  * s23 -- real/complex (type ki), the value of the S matrix element corresponding to the first external off-shell leg
    !  * s13 -- real/complex (type ki), the value of the S matrix element corresponding to the second external off-shell leg
    !  * m3_sq -- real/complex (type ki), the value of the internal mass squared
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
    !  An real (type ki) array of rank 1 and shape 2 corresponding to 
    !  the real/imaginary part of the coefficient of the constant term. If par1 and/or par2
    !  are different from zero, an error is returned.
    !
    ! EXAMPLE
    !
    !
    !*****
    function a3p2m_1mi_np2_r(s23,s13,m3_sq,par1,par2,par3)
      !
      real(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(2) :: a3p2m_1mi_np2_r
      !
      complex(ki) :: c_temp,c_temp_rat
      real(ki) :: sc13,sc23
      !
      a3p2m_1mi_np2_r(:) = 0._ki
      !
      sc13=sign(un,s13+m3_sq)
      !
      sc23=sign(un,s23+m3_sq)
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        c_temp=1._ki/2._ki*(-s13+m3_sq)/(-s13+s23)*z_log(-s13,-1._ki)+1._&
          &ki/2._ki*(-m3_sq+s23)/(-s13+s23)*z_log(-s23,-1._ki)+1._ki/2._ki&
          &*m3_sq**2/s13/(-s13+s23)*q(1,1._ki+m3_sq/s13,-sc13)-3._ki/2._ki&
          &-1._ki/2._ki*m3_sq**2/s23/(-s13+s23)*q(1,1._ki+m3_sq/s23,-sc23)
        !
        c_temp_rat=-3._ki/2._ki
        !
        if ( (rat_or_tot_par%tot_selected)  ) then
          !
          a3p2m_1mi_np2_r = (/real(c_temp,ki),aimag(c_temp)/)
          !
        else !if ( (rat_or_tot_par%rat_selected)  ) then
          !
          a3p2m_1mi_np2_r = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
          !
        end if
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        !
        if (par3 == 1) then
          !
          c_temp=-1._ki/6._ki*(-m3_sq*s23+m3_sq**2+2._ki*s13*s23-s13**2)/(s&
            &13-s23)**2*z_log(-s13,-1._ki)+1._ki/6._ki*(-m3_sq*s23+m3_sq**2+&
            &s23**2)/(s13-s23)**2*z_log(-s23,-1._ki)-1._ki/3._ki/s13*q(1,1._&
            &ki+m3_sq/s13,-sc13)/(s13-s23)**2*m3_sq**3+1._ki/6._ki/s13**2*(s&
            &23+m3_sq)/(s13-s23)**2*m3_sq**3*q(2,1._ki+m3_sq/s13,-sc13)+1._k&
            &i/6._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)**2*m3_sq**3-1.&
            &_ki/18._ki*(3._ki*m3_sq**3+3._ki*s23*m3_sq**2-3._ki*s13*m3_sq**&
            &2-3._ki*s23*s13*m3_sq+3._ki*m3_sq*s13**2-19._ki*s23*s13**2+11._&
            &ki*s23**2*s13+8._ki*s13**3)/s13/(s13-s23)**2
          !
          c_temp_rat=-1._ki/18._ki*(8._ki*s13**2-11._ki*s13*s23+11._ki*s13*&
            &m3_sq-11._ki*m3_sq*s23)/(s13+m3_sq)/(s13-s23)
          !
        else if (par3 == 2) then
          !
          c_temp=1._ki/6._ki*(-s13*m3_sq+m3_sq**2+s13**2)/(s13-s23)**2*z_lo&
            &g(-s13,-1._ki)-1._ki/6._ki*(2._ki*s13*s23-s23**2+m3_sq**2-s13*m&
            &3_sq)/(s13-s23)**2*z_log(-s23,-1._ki)+1._ki/3._ki/s13*q(1,1._ki&
            &+m3_sq/s13,-sc13)/(s13-s23)**2*m3_sq**3-1._ki/6._ki*(s13+m3_sq)&
            &*m3_sq**3/s13**2/(s13-s23)**2*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/&
            &6._ki/s23**2/(s13-s23)*m3_sq**3*q(2,1._ki+m3_sq/s23,-sc23)-1._k&
            &i/6._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)**2*m3_sq**3+1.&
            &_ki/18._ki/s23/s13*(3._ki*m3_sq**3*s23+3._ki*s23*s13*m3_sq**2-3&
            &._ki*m3_sq**2*s13**2+3._ki*m3_sq*s23*s13**2-3._ki*s23**2*s13*m3&
            &_sq-8._ki*s23**3*s13-11._ki*s13**3*s23+19._ki*s23**2*s13**2)/(s&
            &13-s23)**2
          !
          c_temp_rat=-1._ki/18._ki*(11._ki*s13*m3_sq-11._ki*m3_sq*s23-8._ki&
            &*s23**2+11._ki*s13*s23)/(s13-s23)/(s23+m3_sq)
          !
        else if (par3 == 3) then
          !
          c_temp=-1._ki/6._ki*(-s13+2._ki*m3_sq)/(s13-s23)*z_log(-s13,-1._k&
            &i)+1._ki/6._ki*(2._ki*m3_sq-s23)/(s13-s23)*z_log(-s23,-1._ki)-1&
            &._ki/2._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m3_sq**2/(s13-s23)+1.&
            &_ki/6._ki/s13**2/(s13-s23)*m3_sq**3*q(2,1._ki+m3_sq/s13,-sc13)-&
            &1._ki/6._ki/s23**2/(s13-s23)*m3_sq**3*q(2,1._ki+m3_sq/s23,-sc23&
            &)+1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)*m3_sq**2/(s13-s23)&
            &+1._ki/18._ki*(3._ki*m3_sq**2-8._ki*s13*s23)/s13/s23
          !
          c_temp_rat=-1._ki/18._ki*(8._ki*s13*m3_sq+8._ki*s13*s23+5._ki*m3_&
            &sq**2+8._ki*m3_sq*s23)/(s23+m3_sq)/(s13+m3_sq)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'in function a3p2m_1mi_np2_r:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "unimplemented combination of feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
        if ( (rat_or_tot_par%tot_selected)  ) then
          !
          a3p2m_1mi_np2_r = (/real(c_temp,ki),aimag(c_temp)/)
          !
        else !if ( (rat_or_tot_par%rat_selected)  ) then
          !
          a3p2m_1mi_np2_r = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
          !
        end if
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a3p3m_np2:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = &
           &'no need of 3-point integrals in 6 dimension &
           &with more than one Feynman parameter in the numerator'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = &
           &'The value of Feynman parameters in argument: %d1'
        tab_erreur_par(3)%arg_int_tab = (/packb(par),3/)
        call catch_exception(0)
        !
        stop
        !
      end if
      !
    end function a3p2m_1mi_np2_r
    !
    function a3p2m_1mi_np2_c(s23,s13,m3_sq,par1,par2,par3)
      !
      complex(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(2) :: a3p2m_1mi_np2_c
      !
      complex(ki) :: c_temp,c_temp_rat
      real(ki) :: sc13,sc23
      !
      a3p2m_1mi_np2_c(:) = 0._ki
      !
      sc13=sign(un,real(s13+m3_sq,ki))
      !
      sc23=sign(un,real(s23+m3_sq,ki))
      !
      ! cas sans parametre de feynman au numerateur
      if ( (par1 == 0) .and. (par2 == 0) .and. (par3 == 0) ) then
        !
        c_temp=1._ki/2._ki*(-s13+m3_sq)/(-s13+s23)*z_log(-s13,-1._ki)+1._&
          &ki/2._ki*(-m3_sq+s23)/(-s13+s23)*z_log(-s23,-1._ki)+1._ki/2._ki&
          &*m3_sq**2/s13/(-s13+s23)*q(1,1._ki+m3_sq/s13,-sc13)-3._ki/2._ki&
          &-1._ki/2._ki*m3_sq**2/s23/(-s13+s23)*q(1,1._ki+m3_sq/s23,-sc23)
        !
        c_temp_rat=-3._ki/2._ki
        !
        if ( (rat_or_tot_par%tot_selected)  ) then
          !
          a3p2m_1mi_np2_c = (/real(c_temp,ki),aimag(c_temp)/)
          !
        else if ( (rat_or_tot_par%rat_selected)  ) then
          !
          a3p2m_1mi_np2_c = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'in function a3p2m_1mi_np2_c:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "unimplemented combination of feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      ! cas avec un parametre de feynman au numerateur
      else if ( (par1 == 0) .and. (par2 == 0) ) then
        !
        !
        if (par3 == 1) then
          !
          c_temp=-1._ki/6._ki*(-m3_sq*s23+m3_sq**2+2._ki*s13*s23-s13**2)/(s&
            &13-s23)**2*z_log(-s13,-1._ki)+1._ki/6._ki*(-m3_sq*s23+m3_sq**2+&
            &s23**2)/(s13-s23)**2*z_log(-s23,-1._ki)-1._ki/3._ki/s13*q(1,1._&
            &ki+m3_sq/s13,-sc13)/(s13-s23)**2*m3_sq**3+1._ki/6._ki/s13**2*(s&
            &23+m3_sq)/(s13-s23)**2*m3_sq**3*q(2,1._ki+m3_sq/s13,-sc13)+1._k&
            &i/6._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)**2*m3_sq**3-1.&
            &_ki/18._ki*(3._ki*m3_sq**3+3._ki*s23*m3_sq**2-3._ki*s13*m3_sq**&
            &2-3._ki*s23*s13*m3_sq+3._ki*m3_sq*s13**2-19._ki*s23*s13**2+11._&
            &ki*s23**2*s13+8._ki*s13**3)/s13/(s13-s23)**2
          !
          c_temp_rat=-1._ki/18._ki*(8._ki*s13**2-11._ki*s13*s23+11._ki*s13*&
            &m3_sq-11._ki*m3_sq*s23)/(s13+m3_sq)/(s13-s23)
          !
        else if (par3 == 2) then
          !
          c_temp=1._ki/6._ki*(-s13*m3_sq+m3_sq**2+s13**2)/(s13-s23)**2*z_lo&
            &g(-s13,-1._ki)-1._ki/6._ki*(2._ki*s13*s23-s23**2+m3_sq**2-s13*m&
            &3_sq)/(s13-s23)**2*z_log(-s23,-1._ki)+1._ki/3._ki/s13*q(1,1._ki&
            &+m3_sq/s13,-sc13)/(s13-s23)**2*m3_sq**3-1._ki/6._ki*(s13+m3_sq)&
            &*m3_sq**3/s13**2/(s13-s23)**2*q(2,1._ki+m3_sq/s13,-sc13)+1._ki/&
            &6._ki/s23**2/(s13-s23)*m3_sq**3*q(2,1._ki+m3_sq/s23,-sc23)-1._k&
            &i/6._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)/(s13-s23)**2*m3_sq**3+1.&
            &_ki/18._ki/s23/s13*(3._ki*m3_sq**3*s23+3._ki*s23*s13*m3_sq**2-3&
            &._ki*m3_sq**2*s13**2+3._ki*m3_sq*s23*s13**2-3._ki*s23**2*s13*m3&
            &_sq-8._ki*s23**3*s13-11._ki*s13**3*s23+19._ki*s23**2*s13**2)/(s&
            &13-s23)**2
          !
          c_temp_rat=-1._ki/18._ki*(11._ki*s13*m3_sq-11._ki*m3_sq*s23-8._ki&
            &*s23**2+11._ki*s13*s23)/(s13-s23)/(s23+m3_sq)
          !
        else if (par3 == 3) then
          !
          c_temp=-1._ki/6._ki*(-s13+2._ki*m3_sq)/(s13-s23)*z_log(-s13,-1._k&
            &i)+1._ki/6._ki*(2._ki*m3_sq-s23)/(s13-s23)*z_log(-s23,-1._ki)-1&
            &._ki/2._ki/s13*q(1,1._ki+m3_sq/s13,-sc13)*m3_sq**2/(s13-s23)+1.&
            &_ki/6._ki/s13**2/(s13-s23)*m3_sq**3*q(2,1._ki+m3_sq/s13,-sc13)-&
            &1._ki/6._ki/s23**2/(s13-s23)*m3_sq**3*q(2,1._ki+m3_sq/s23,-sc23&
            &)+1._ki/2._ki/s23*q(1,1._ki+m3_sq/s23,-sc23)*m3_sq**2/(s13-s23)&
            &+1._ki/18._ki*(3._ki*m3_sq**2-8._ki*s13*s23)/s13/s23
          !
          c_temp_rat=-1._ki/18._ki*(8._ki*s13*m3_sq+8._ki*s13*s23+5._ki*m3_&
            &sq**2+8._ki*m3_sq*s23)/(s23+m3_sq)/(s13+m3_sq)
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'in function a3p2m_1mi_np2_c:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "unimplemented combination of feynman parameters"
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'par1=%d0, par2,par3 = %d1'
          tab_erreur_par(3)%arg_int = par1
          tab_erreur_par(3)%arg_int_tab = (/par2,par3/)
          call catch_exception(0)
          !
          stop
          !
        end if
        !
        if ( (rat_or_tot_par%tot_selected)  ) then
          !
          a3p2m_1mi_np2_c = (/real(c_temp,ki),aimag(c_temp)/)
          !
        else !if ( (rat_or_tot_par%rat_selected)  ) then
          !
          a3p2m_1mi_np2_c = (/real(c_temp_rat,ki),aimag(c_temp_rat)/)
          !
        end if
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a3p3m_np2:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = &
             &'no need of 3-point integrals in 6 dimension &
             &with more than one Feynman parameter in the numerator'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = &
             &'The value of Feynman parameters in argument: %d1'
        tab_erreur_par(3)%arg_int_tab = (/packb(par),3/)
        call catch_exception(0)
        !
        stop
        !
      end if
      !
    end function a3p2m_1mi_np2_c
    !
    !
    !****if* src/integral/three_point/function_3p2m_1mi/eval_numer_gi
    ! NAME
    !
    !  Function eval_numer_gi
    !
    ! USAGE
    !
    !  complex = eval_numer_gi(u)
    !
    ! DESCRIPTION
    !
    !  This function is the integrand that will be computed numerically
    !
    ! INPUTS
    !
    !  * u -- a real (type ki), the integral variable
    !
    ! SIDE EFFECTS
    !
    !  No side effect, use the values of the local (for this module) variables
    !  eps_glob,s23_glob,s13_glob,s23_glob,par1_glob,par2_glob,par3_glob,dim_glob
    !  and also the global variables alpha_par,beta_par and lambda_par given
    !  by the module parametre (src/module/parametre.f90)
    !
    ! RETURN VALUE
    !
    !  It returns a complex (type ki) which is the value of the
    !  integrand at the value u
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function eval_numer_gi(u)
      !
      real(ki), intent (in) :: u
      complex(ki) :: eval_numer_gi
      !
      real(ki) :: x,y
      complex(ki) :: z,jacob
      !
      x = u
      y = lambda_par*u**alpha_par*(1._ki-u)**beta_par
      z = x - eps_glob*i_*y
      jacob = 1._ki - eps_glob*i_*lambda_par*u**(alpha_par-1._ki)&
              *(1._ki-u)**(beta_par-1._ki)*(alpha_par*(1._ki-u)-beta_par*u)
      !
      eval_numer_gi = fg(z,s23_glob,s13_glob,m3_sq_glob,&
                      &  par1_glob,par2_glob,par3_glob,&
                      &  dim_glob)
      eval_numer_gi = eval_numer_gi*jacob
      !
    end function eval_numer_gi
    !
    !****if* src/integral/three_point/function_3p2m_1mi/fg
    ! NAME
    !
    !  Function fg
    !
    ! USAGE
    !
    !  complex = fg(z,s23,s13,m3_sq,par1,par2,par3,dim)
    !
    ! DESCRIPTION
    !
    !  This function gives the structure of the integrand for the different cases
    !
    ! INPUTS
    !
    !  * z -- a complex (type ki), the integral variable
    !  * s23 -- complex (type ki), the value of the S matrix element corresponding to the first external off-shell leg
    !  * s13 -- complex (type ki), the value of the S matrix element corresponding to the second external off-shell leg
    !  * m3_sq -- complex (type ki), the value of the internal mass squared
    !  * par1 -- an integer, the label of the third Feynman parameter
    !  * par2 -- an integer, the label of the second Feynman parameter
    !  * par3 -- an integer, the label of the first Feynman parameter
    !  * dim -- a character (length 3), to compute in n or n+2 dimensions, 
    !    the values are "ndi", "n+2"
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  It returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function fg(z,s23,s13,m3_sq,par1,par2,par3,dim)
      !
      complex(ki), intent (in) :: z
      complex(ki), intent (in) :: s23,s13,m3_sq
      integer, intent (in) :: par1,par2,par3
      character (len=3) :: dim
      complex(ki) :: fg
      !
      integer, dimension(3) :: par
      integer :: nb_par
      complex(ki) :: d1_var,d2_var
      !
      par = (/par1,par2,par3/)
      nb_par = count(mask=par/=0)
      !
      d1_var=-z*(s13-s23)-s23
      !
      d2_var=-z*(s13-s23)-s23-m3_sq
      !
      if (dim == "ndi") then
        !
        if (nb_par == 0) then
          !
          fg=log(d1_var)*m3_sq/d1_var/d2_var-log(m3_sq)*m3_sq/d1_var/d2_var
          !
        else if (nb_par == 1) then
          !
          select case(par3)
          !
          case(1)
            !
            fg=m3_sq/d1_var*z/d2_var-log(d1_var)*m3_sq**2/d1_var*z/d2_var**2+&
              &log(m3_sq)*m3_sq**2/d1_var*z/d2_var**2+z/d1_var
            !
          case(2)
            !
            fg=m3_sq/d1_var/d2_var-log(d1_var)*m3_sq**2/d1_var/d2_var**2+log(&
              &m3_sq)*m3_sq**2/d1_var/d2_var**2-m3_sq/d1_var*z/d2_var+log(d1_v&
              &ar)*m3_sq**2/d1_var*z/d2_var**2-log(m3_sq)*m3_sq**2/d1_var*z/d2&
              &_var**2-z/d1_var+1._ki/d1_var
            !
          case(3)
            !
            fg=-1._ki/d2_var+m3_sq*log(d1_var)/d2_var**2-m3_sq*log(m3_sq)/d2_&
              &var**2
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'in function fg:'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
               & "par3 should be 1, 2 or 3 but is %d0"
            tab_erreur_par(2)%arg_int = par3
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else if (nb_par == 2) then
          !
          select case(par2)
          !
          case(1)
            !
            select case(par3)
            !
            case(1)
              !
              fg=1._ki/2._ki*m3_sq/d1_var*z**2/d2_var-m3_sq**2/d1_var*z**2/d2_v&
                &ar**2+log(d1_var)*m3_sq**3/d1_var*z**2/d2_var**3-log(m3_sq)*m3_&
                &sq**3/d1_var*z**2/d2_var**3+3._ki/2._ki*z**2/d1_var
              !
            case(2)
              !
              fg=log(d1_var)*m3_sq**3/d1_var*z/d2_var**3-log(m3_sq)*m3_sq**3/d1&
                &_var*z/d2_var**3-log(d1_var)*m3_sq**3/d1_var*z**2/d2_var**3+log&
                &(m3_sq)*m3_sq**3/d1_var*z**2/d2_var**3+1._ki/2._ki*m3_sq/d1_var&
                &*z/d2_var-m3_sq**2/d1_var*z/d2_var**2-1._ki/2._ki*m3_sq/d1_var*&
                &z**2/d2_var+m3_sq**2/d1_var*z**2/d2_var**2+3._ki/2._ki*z/d1_var&
                &-3._ki/2._ki*z**2/d1_var
              !
            case(3)
              !
              fg=-1._ki/2._ki*z/d2_var+m3_sq*z/d2_var**2-m3_sq**2*log(d1_var)*z&
                &/d2_var**3+m3_sq**2*log(m3_sq)*z/d2_var**3
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'in function fg:'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
                 & "par3 should be 1, 2 or 3 but is %d0"
              tab_erreur_par(2)%arg_int = par3
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(2)
            !
            select case(par3)
            !
            case(2)
              !
              fg=1._ki/2._ki*m3_sq/d1_var/d2_var-log(m3_sq)*m3_sq**3/d1_var/d2_&
                &var**3-m3_sq/d1_var*z/d2_var+2._ki*m3_sq**2/d1_var*z/d2_var**2+&
                &1._ki/2._ki*m3_sq/d1_var*z**2/d2_var-m3_sq**2/d1_var*z**2/d2_va&
                &r**2+log(d1_var)*m3_sq**3/d1_var/d2_var**3-m3_sq**2/d1_var/d2_v&
                &ar**2-2._ki*log(d1_var)*m3_sq**3/d1_var*z/d2_var**3+2._ki*log(m&
                &3_sq)*m3_sq**3/d1_var*z/d2_var**3+log(d1_var)*m3_sq**3/d1_var*z&
                &**2/d2_var**3-log(m3_sq)*m3_sq**3/d1_var*z**2/d2_var**3-3._ki*z&
                &/d1_var+3._ki/2._ki*z**2/d1_var+3._ki/2._ki/d1_var
              !
            case(3)
              !
              fg=-1._ki/2._ki/d2_var+m3_sq/d2_var**2-m3_sq**2*log(d1_var)/d2_va&
                &r**3+m3_sq**2*log(m3_sq)/d2_var**3+1._ki/2._ki*z/d2_var-m3_sq*z&
                &/d2_var**2+m3_sq**2*log(d1_var)*z/d2_var**3-m3_sq**2*log(m3_sq)&
                &*z/d2_var**3
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'in function fg:'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
                 & "par3 should be 2 or 3 but is %d0"
              tab_erreur_par(2)%arg_int = par3
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(3)
            !
            select case(par3)
            !
            case(3)
              !
              fg=-1._ki/2._ki/d2_var+m3_sq*log(d1_var)/d2_var**2-m3_sq*log(m3_s&
                &q)/d2_var**2-m3_sq/d2_var**2+m3_sq**2*log(d1_var)/d2_var**3-m3_&
                &sq**2*log(m3_sq)/d2_var**3
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'in function fg:'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
                 & "par3 should be 3 but is %d0"
              tab_erreur_par(2)%arg_int = par3
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'in function fg:'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
               & "par2 should be 1, 2 or 3 but is %d0"
            tab_erreur_par(2)%arg_int = par2
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else if (nb_par == 3) then
          !
          select case(par1)
          !
          case(1)
            !
            select case(par2)
            !
            case(1)
              !
              select case(par3)
              !
              case(1)
                !
                fg=1._ki/3._ki*m3_sq/d1_var*z**3/d2_var-1._ki/2._ki*m3_sq**2/d1_v&
                  &ar*z**3/d2_var**2+m3_sq**3/d1_var*z**3/d2_var**3-log(d1_var)*m3&
                  &_sq**4/d1_var*z**3/d2_var**4+log(m3_sq)*m3_sq**4/d1_var*z**3/d2&
                  &_var**4+11._ki/6._ki*z**3/d1_var
                !
              case(2)
                !
                fg=1._ki/3._ki*m3_sq/d1_var*z**2/d2_var-1._ki/2._ki*m3_sq**2/d1_v&
                  &ar*z**2/d2_var**2+m3_sq**3/d1_var*z**2/d2_var**3-1._ki/3._ki*m3&
                  &_sq/d1_var*z**3/d2_var+1._ki/2._ki*m3_sq**2/d1_var*z**3/d2_var*&
                  &*2-m3_sq**3/d1_var*z**3/d2_var**3-log(d1_var)*m3_sq**4/d1_var*z&
                  &**2/d2_var**4+log(m3_sq)*m3_sq**4/d1_var*z**2/d2_var**4+log(d1_&
                  &var)*m3_sq**4/d1_var*z**3/d2_var**4-log(m3_sq)*m3_sq**4/d1_var*&
                  &z**3/d2_var**4+11._ki/6._ki*z**2/d1_var-11._ki/6._ki*z**3/d1_va&
                  &r
                !
              case(3)
                !
                fg=-1._ki/3._ki*z**2/d2_var+1._ki/2._ki*m3_sq*z**2/d2_var**2-m3_s&
                  &q**2*z**2/d2_var**3+m3_sq**3*log(d1_var)*z**2/d2_var**4-m3_sq**&
                  &3*log(m3_sq)*z**2/d2_var**4
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'in function fg:'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                   & "par3 should be 1, 2 or 3 but is %d0"
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(2)
              !
              select case(par3)
              !
              case(2)
                !
                fg=1._ki/3._ki*m3_sq/d1_var*z/d2_var-1._ki/2._ki*m3_sq**2/d1_var*&
                  &z/d2_var**2+m3_sq**3/d1_var*z/d2_var**3-2._ki/3._ki*m3_sq/d1_va&
                  &r*z**2/d2_var+m3_sq**2/d1_var*z**2/d2_var**2-2._ki*m3_sq**3/d1_&
                  &var*z**2/d2_var**3+1._ki/3._ki*m3_sq/d1_var*z**3/d2_var-1._ki/2&
                  &._ki*m3_sq**2/d1_var*z**3/d2_var**2+m3_sq**3/d1_var*z**3/d2_var&
                  &**3-log(d1_var)*m3_sq**4/d1_var*z/d2_var**4+log(m3_sq)*m3_sq**4&
                  &/d1_var*z/d2_var**4+2._ki*log(d1_var)*m3_sq**4/d1_var*z**2/d2_v&
                  &ar**4-2._ki*log(m3_sq)*m3_sq**4/d1_var*z**2/d2_var**4-log(d1_va&
                  &r)*m3_sq**4/d1_var*z**3/d2_var**4+log(m3_sq)*m3_sq**4/d1_var*z*&
                  &*3/d2_var**4+11._ki/6._ki*z/d1_var-11._ki/3._ki*z**2/d1_var+11.&
                  &_ki/6._ki*z**3/d1_var
                !
              case(3)
                !
                fg=-1._ki/3._ki*z/d2_var+1._ki/2._ki*m3_sq*z/d2_var**2-m3_sq**2*z&
                  &/d2_var**3+m3_sq**3*log(d1_var)*z/d2_var**4-m3_sq**3*log(m3_sq)&
                  &*z/d2_var**4+1._ki/3._ki*z**2/d2_var-1._ki/2._ki*m3_sq*z**2/d2_&
                  &var**2+m3_sq**2*z**2/d2_var**3-m3_sq**3*log(d1_var)*z**2/d2_var&
                  &**4+m3_sq**3*log(m3_sq)*z**2/d2_var**4
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'in function fg:'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                   & "par3 should be 2 or 3 but is %d0"
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(3)
              !
              select case(par3)
              !
              case(3)
                !
                fg=-1._ki/6._ki*z/d2_var+1._ki/2._ki*m3_sq*z/d2_var**2-m3_sq**2*l&
                  &og(d1_var)*z/d2_var**3+m3_sq**2*log(m3_sq)*z/d2_var**3+m3_sq**2&
                  &*z/d2_var**3-m3_sq**3*log(d1_var)*z/d2_var**4+m3_sq**3*log(m3_s&
                  &q)*z/d2_var**4
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'in function fg:'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                   & "par3 should be 3 but is %d0"
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'in function fg:'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
                 & "par2 should be 1, 2 or 3 but is %d0"
              tab_erreur_par(2)%arg_int = par2
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(2)
            !
            select case(par2)
            !
            case(2)
              !
              select case(par3)
              !
              case(2)
                !
                fg=-3._ki*log(d1_var)*m3_sq**4/d1_var*z**2/d2_var**4+3._ki*log(m3&
                  &_sq)*m3_sq**4/d1_var*z**2/d2_var**4+3._ki*log(d1_var)*m3_sq**4/&
                  &d1_var*z/d2_var**4-3._ki*log(m3_sq)*m3_sq**4/d1_var*z/d2_var**4&
                  &+m3_sq/d1_var*z**2/d2_var-1._ki/3._ki*m3_sq/d1_var*z**3/d2_var-&
                  &m3_sq/d1_var*z/d2_var-3._ki/2._ki*m3_sq**2/d1_var*z**2/d2_var**&
                  &2+1._ki/2._ki*m3_sq**2/d1_var*z**3/d2_var**2-m3_sq**3/d1_var*z*&
                  &*3/d2_var**3+log(d1_var)*m3_sq**4/d1_var*z**3/d2_var**4-log(m3_&
                  &sq)*m3_sq**4/d1_var*z**3/d2_var**4+3._ki*m3_sq**3/d1_var*z**2/d&
                  &2_var**3-3._ki*m3_sq**3/d1_var*z/d2_var**3+1._ki/3._ki*m3_sq/d1&
                  &_var/d2_var-1._ki/2._ki*m3_sq**2/d1_var/d2_var**2+m3_sq**3/d1_v&
                  &ar/d2_var**3+3._ki/2._ki*m3_sq**2/d1_var*z/d2_var**2-log(d1_var&
                  &)*m3_sq**4/d1_var/d2_var**4+log(m3_sq)*m3_sq**4/d1_var/d2_var**&
                  &4-11._ki/2._ki*z/d1_var+11._ki/2._ki*z**2/d1_var-11._ki/6._ki*z&
                  &**3/d1_var+11._ki/6._ki/d1_var
                !
              case(3)
                !
                fg=-1._ki/3._ki/d2_var+1._ki/2._ki*m3_sq/d2_var**2-m3_sq**2/d2_va&
                  &r**3+m3_sq**3*log(d1_var)/d2_var**4-m3_sq**3*log(m3_sq)/d2_var*&
                  &*4+2._ki/3._ki*z/d2_var-m3_sq*z/d2_var**2+2._ki*m3_sq**2*z/d2_v&
                  &ar**3-2._ki*m3_sq**3*log(d1_var)*z/d2_var**4+2._ki*m3_sq**3*log&
                  &(m3_sq)*z/d2_var**4-1._ki/3._ki*z**2/d2_var+1._ki/2._ki*m3_sq*z&
                  &**2/d2_var**2-m3_sq**2*z**2/d2_var**3+m3_sq**3*log(d1_var)*z**2&
                  &/d2_var**4-m3_sq**3*log(m3_sq)*z**2/d2_var**4
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'in function fg:'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                   & "par3 should be 2 or 3 but is %d0"
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(3)
              !
              select case(par3)
              !
              case(3)
                !
                fg=-1._ki/6._ki/d2_var+1._ki/2._ki*m3_sq/d2_var**2-m3_sq**2*log(d&
                  &1_var)/d2_var**3+m3_sq**2*log(m3_sq)/d2_var**3+m3_sq**2/d2_var*&
                  &*3-m3_sq**3*log(d1_var)/d2_var**4+m3_sq**3*log(m3_sq)/d2_var**4&
                  &+1._ki/6._ki*z/d2_var-1._ki/2._ki*m3_sq*z/d2_var**2+m3_sq**2*lo&
                  &g(d1_var)*z/d2_var**3-m3_sq**2*log(m3_sq)*z/d2_var**3-m3_sq**2*&
                  &z/d2_var**3+m3_sq**3*log(d1_var)*z/d2_var**4-m3_sq**3*log(m3_sq&
                  &)*z/d2_var**4
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'in function fg:'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                   & "par3 should be 3 but is %d0"
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'in function fg:'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
                 & "par2 should be 2 or 3 but is %d0"
              tab_erreur_par(2)%arg_int = par2
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(3)
            !
            select case(par2)
            !
            case(3)
              !
              select case(par3)
              !
              case(3)
                !
                fg=-1._ki/3._ki/d2_var+m3_sq*log(d1_var)/d2_var**2-m3_sq*log(m3_s&
                  &q)/d2_var**2-3._ki/2._ki*m3_sq/d2_var**2+2._ki*m3_sq**2*log(d1_&
                  &var)/d2_var**3-2._ki*m3_sq**2*log(m3_sq)/d2_var**3-m3_sq**2/d2_&
                  &var**3+m3_sq**3*log(d1_var)/d2_var**4-m3_sq**3*log(m3_sq)/d2_va&
                  &r**4
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'in function fg:'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                   & "par3 should be 3 but is %d0"
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'in function fg:'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
                 & "par2 should be 3 but is %d0"
              tab_erreur_par(2)%arg_int = par2
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'in function fg:'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
               & "par1 should be 1, 2 or 3 but is %d0"
            tab_erreur_par(2)%arg_int = par1
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'in function fg:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "Unexpected value for nb_par = %d0"
          tab_erreur_par(2)%arg_int = nb_par
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      else if (dim == "n+2") then
        !
        if (nb_par == 0) then
          !
          fg=1._ki/2._ki*m3_sq**2*log(m3_sq)/d2_var**2-1._ki/2._ki*(m3_sq-d&
            &2_var)*(m3_sq+d2_var)/d2_var**2*log(d1_var)+1._ki/2._ki*(m3_sq-&
            &2._ki*d2_var)/d2_var
          !
        else if (nb_par == 1) then
          !
          select case(par3)
          !
          case(1)
            !
            fg=-1._ki/3._ki*m3_sq**3*log(m3_sq)*z/d2_var**3+1._ki/3._ki*log(d&
              &1_var)/d2_var**3*(m3_sq+d2_var)*(m3_sq**2-m3_sq*d2_var+d2_var**&
              &2)*z-1._ki/18._ki/d2_var**2*(-3._ki*m3_sq*d2_var+6._ki*m3_sq**2&
              &+13._ki*d2_var**2)*z
            !
          case(2)
            !
            fg=-1._ki/3._ki*m3_sq**3*log(m3_sq)/d2_var**3*(1._ki-z)+1._ki/3._&
              &ki*log(d1_var)/d2_var**3*(m3_sq+d2_var)*(m3_sq**2-m3_sq*d2_var+&
              &d2_var**2)*(1._ki-z)-1._ki/18._ki/d2_var**2*(-3._ki*m3_sq*d2_va&
              &r+6._ki*m3_sq**2+13._ki*d2_var**2)*(1._ki-z)
            !
          case(3)
            !
            fg=1._ki/6._ki*m3_sq**2*(3._ki*d2_var+2._ki*m3_sq)/d2_var**3*log(&
              &m3_sq)-1._ki/6._ki*(2._ki*m3_sq-d2_var)*(m3_sq+d2_var)**2/d2_va&
              &r**3*log(d1_var)+1._ki/18._ki*(6._ki*m3_sq*d2_var-5._ki*d2_var*&
              &*2+6._ki*m3_sq**2)/d2_var**2
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'in function fg:'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
               & "par3 should be 1, 2 or 3 but is %d0"
            tab_erreur_par(2)%arg_int = par3
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'in function fg:'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
             & "Unexpected value for nb_par = %d0"
          tab_erreur_par(2)%arg_int = nb_par
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'in function fg:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = &
           & "dim is %c0"
        tab_erreur_par(2)%arg_char = dim
        call catch_exception(0)
        !
        stop
        !
      end if
      !
    end function fg
    !
end module function_3p2m_1mi
