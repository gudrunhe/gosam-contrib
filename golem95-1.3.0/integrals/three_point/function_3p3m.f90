! 
!****h* src/integral/three_point/function_3p3m
! NAME
!
!  Module function_3p3m
!
! USAGE
!
!  use function_3p3m
!
! DESCRIPTION
!
!  This module is used to compute the three off-shell external leg three point function
!  with no internal leg with/without Feynman parameters in n, n+2 dimensions
!
! OUTPUT
!
!  This module exports three functions:
!  * f3p3m -- a function for the computation of the three mass three 
!    point function with/without Feynman parameters in n, n+2 dimensions
!  * f3p3m_c -- a function which computes the same thing as f3p3m, only 
!    the format of the return values is different
!  * i3_3mass -- a function for the computation of the scalar three mass three 
!    point function in n dimensions
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * numerical_evaluation (src/numerical/mod_numeric.f90)
!  * dilogarithme (src/module/zdilog.f90)
!  * logarithme (src/module/z_log.f90)
!  * constante (src/module/constante.f90)
!  * parametre (src/module/parametre.f90)
!  * array (src/module/array.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * generic_function_2p (src/integrals/two_point/generic_function_2p.f90)
!  * multiply_div (src/module/multiply_div.f90)
!  * s_matrix_type (src/module/s_matrix_type.f90)
!  * matrice_s (src/kinematic/matrice_s.f90)
!
!*****
module function_3p3m
  !
  use precision_golem
  use numerical_evaluation
  use dilogarithme
  use logarithme
  use constante
  use parametre
  use array
  use sortie_erreur
  use generic_function_2p
  use multiply_div
  use s_matrix_type
  use matrice_s, only : prepare_s_matrix_local, find_plus_grand
  implicit none
  !
  private
  !
  real(ki) :: s12_glob,s23_glob,s13_glob
  real(ki) :: eps_glob
  integer :: par1_glob,par2_glob,par3_glob
  character (len=3) :: dim_glob
  !
  real(ki), dimension(3) :: b
  real(ki) :: sumb
  real(ki), dimension(3,3) :: invs,s_mat_loc
  integer, dimension(3) :: par
  integer, dimension(3) :: s = (/1,2,3/)
  type (s_matrix_poly) :: s_mat_p_loc
  !
  logical, dimension(:), allocatable :: deja_calcule
  real(ki),dimension(:,:), allocatable :: resultat
  logical, dimension(:,:), allocatable :: deja_calcule2
  real(ki),dimension(:,:,:), allocatable :: resultat2
  logical, dimension(:), allocatable :: deja_calcule_np2
  real(ki),dimension(:,:), allocatable :: resultat_np2
  logical, dimension(:,:,:), allocatable :: deja_calcule22
  real(ki),dimension(:,:,:,:), allocatable :: resultat22
  !
 public :: f3p3m,i3_3mass,f3p3m_c
  !
  contains
    !
    !****f* src/integral/three_point/function_3p3m/f3p3m
    ! NAME
    !
    !  Function f3p3m
    !
    ! USAGE
    !
    !  real_dim4 = f3p3m(dim,m1,m2,m3,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This function computes the three off-shell external leg three point function in n 
    !  and n+2 dimension. It uses the formula of ref. 
    !  It switches to numerical evaluation if the Gram determinant is smaller than
    !  coupure_3p3m (in src/module/parametre.f90)
    !
    ! INPUTS
    !
    !  * dim -- a character (length 3), to compute in n or n+2 dimensions, 
    !    the values are "ndi", "n+2"
    !  * m1 -- a real (type ki), the first mass squared
    !  * m2 -- a real (type ki), the second mass squared
    !  * m3 -- a real (type ki), the third mass squared
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
    !  are different from zero for dim="n+2", an error is returned.
    !
    ! EXAMPLE
    !
    ! three mass three point function without Feynman parameters in n dimensions
    ! f3p3m("ndi",m1,m2,m3,0,0,0) 
    ! with one Feynman parameter at the numerator z_1 in n dimensions 
    ! f3p3m("ndi",m1,m2,m3,0,0,1) 
    ! with three Feynman parameters at the numerator z_2^2 z_3 in n dimensions 
    ! f3p3m("ndi",m1,m2,m3,2,2,3) 
    ! three mass three point function without Feynman parameters in n+2 dimensions 
    ! f3p3m("n+2",m1,m2,m3,0,0,0) 
    ! with one Feynman parameter at the numerator z_1 in n+2 dimensions 
    ! f3p3m("n+2",m1,m2,m3,0,0,1) 
    !
    !*****
    function f3p3m(dim,m1,m2,m3,par1,par2,par3)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: m1,m2,m3
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(4) :: f3p3m
      !
      integer :: nb_par
      real(ki) :: lamb
      real(ki) :: plus_grand
      real(ki) :: norma
      complex(ki) :: resto,abserro
      !
      par = (/par1,par2,par3/)
      !
      s_mat_loc(1,:) = (/0._ki,m2,m1/)
      s_mat_loc(2,:) = (/m2,0._ki,m3/)
      s_mat_loc(3,:) = (/m1,m3,0._ki/)
      !
      ! on redefinit la matrice S de telle facon a ce que ces elements
      ! soient entre -1 et 1
      !
      if (rat_or_tot_par%tot_selected) then
        !
        plus_grand = find_plus_grand(array=abs(s_mat_loc))
        !
      else !if (rat_or_tot_par%rat_selected) then
        !
        plus_grand = 1._ki
        !
      end if
      !
      s_mat_loc = s_mat_loc/plus_grand
      !
      s_mat_p_loc = assign_s_matrix(s_mat_loc)
      call prepare_s_matrix_local(s_mat_p_loc,s)
      !
      b(1) = (s_mat_loc(1,3)+s_mat_loc(1,2)-s_mat_loc(2,3))/(2._ki*s_mat_loc(1,3)*s_mat_loc(1,2))
      b(2) = (s_mat_loc(1,2)+s_mat_loc(2,3)-s_mat_loc(1,3))/(2._ki*s_mat_loc(1,2)*s_mat_loc(2,3))
      b(3) = (s_mat_loc(1,3)+s_mat_loc(2,3)-s_mat_loc(1,2))/(2._ki*s_mat_loc(1,3)*s_mat_loc(2,3))
      !
      sumb = (2._ki*s_mat_loc(1,3)*s_mat_loc(2,3)+2._ki*s_mat_loc(1,2)*s_mat_loc(2,3)&
      &+2._ki*s_mat_loc(1,3)*s_mat_loc(1,2)-s_mat_loc(1,3)*s_mat_loc(1,3)-s_mat_loc(1,2)*s_mat_loc(1,2)&
      &-s_mat_loc(2,3)*s_mat_loc(2,3))/(2._ki*s_mat_loc(1,3)*s_mat_loc(1,2)*s_mat_loc(2,3))
      !
      invs(1,1) = -s_mat_loc(2,3)/s_mat_loc(1,2)/s_mat_loc(1,3)/2._ki
      invs(1,2) = 1._ki/s_mat_loc(1,2)/2._ki
      invs(1,3) = 1._ki/s_mat_loc(1,3)/2._ki
      invs(2,1) = 1._ki/s_mat_loc(1,2)/2._ki
      invs(2,2) = -s_mat_loc(1,3)/s_mat_loc(1,2)/s_mat_loc(2,3)/2._ki
      invs(2,3) = 1._ki/s_mat_loc(2,3)/2._ki
      invs(3,1) = 1._ki/s_mat_loc(1,3)/2._ki
      invs(3,2) = 1._ki/s_mat_loc(2,3)/2._ki
      invs(3,3) = -s_mat_loc(1,2)/s_mat_loc(1,3)/s_mat_loc(2,3)/2._ki
      !
      lamb = 2._ki*s_mat_loc(1,3)*s_mat_loc(2,3)+2._ki*s_mat_loc(1,2)*s_mat_loc(2,3)&
      +2._ki*s_mat_loc(1,3)*s_mat_loc(1,2)-s_mat_loc(1,3)*s_mat_loc(1,3)-s_mat_loc(1,2)*s_mat_loc(1,2)&
      -s_mat_loc(2,3)*s_mat_loc(2,3)
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
      ! memory allocation to save time in the recursion
      !
      allocate(deja_calcule(4))
      allocate(resultat(4,2))
      allocate(deja_calcule2(3,4))
      allocate(resultat2(3,4,4))
      allocate(deja_calcule_np2(4))
      allocate(resultat_np2(4,4))
      allocate(deja_calcule22(3,4,4))
      allocate(resultat22(3,4,4,4))
      !
      ! initialisation
      !
      deja_calcule = .false.
      resultat = 0._ki
      deja_calcule2 = .false.
      resultat2 = 0._ki
      deja_calcule_np2 = .false.
      resultat_np2 = 0._ki
      deja_calcule22 = .false.
      resultat22 = 0._ki
      !
      f3p3m = 0._ki
      !
      if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_4p1m) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f3p3m (in file function_3p3m.f90): &
        &the flag rat to compute the rational part is on &
        &and the program reachs a region of phase space in &
        &which det(G) = 0  Be careful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to go on, he has to &
        &reduce the value of the parameter coupure_3p3m'
        call catch_exception(0)
        !
        stop
        !
      end if
      !
      if (abs(sumb) > coupure_3p3m) then
        !
        ! analytic computation
        !
        if (dim == "ndi") then
          !
          f3p3m(3:4)= a3p3m(s_mat_loc(1,3),s_mat_loc(1,2),s_mat_loc(2,3),par1,par2,par3)&
                      &/plus_grand
          !
        else if (dim == "n+2") then
          !
          f3p3m = a3p3m_np2(s_mat_loc(1,3),s_mat_loc(1,2),s_mat_loc(2,3),par1,par2,par3)
          f3p3m(3) = f3p3m(3)-log(plus_grand)*norma
          ! mu2_scale_par is already contained in the bubbles, 
          ! but scaling of s_mat_loc still needs to be undone
          !
        else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p3m (function_3p3m.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'dim = %c0'
          tab_erreur_par(2)%arg_char = dim
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      else
        !
        ! numerical computation
        !
        dim_glob = dim
        par1_glob = par1
        par2_glob = par2
        par3_glob = par3
        !
        s13_glob = s_mat_loc(1,3)
        s12_glob = s_mat_loc(1,2)
        s23_glob = s_mat_loc(2,3)
        !
        resto = 0._ki
        abserro = 0._ki
        !
        ! on pose z = x - i*eps*y (avec x et y > 0)
        ! z*s13+(1-z)*s23 = s23+x*(s13-s23)-i*eps*y*(s13-s23)
        ! on veut la partie imaginaire du meme signe que i*lambda
        ! => eps*(s13-s23) < 0
        !
        ! faire attention que suivant le signe de eps_glob, on tourne dans le
        ! sens des aiguilles d'une montre ou inversement
        ! eps_glob = 1, on ferme le contour vers le bas --> -2 i Pi residu
        ! eps_glob = -1, on ferme le contour vers le haut --> +2 i Pi residu
        !
        eps_glob = sign(1._ki,s23_glob-s13_glob)
        !
        origine_info_par = "f3p3m, dimension "//dim
        num_grand_b_info_par = lamb
        denom_grand_b_info_par = (2._ki*s_mat_loc(1,3)*s_mat_loc(1,2)*s_mat_loc(2,3))
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,resto,abserro)
        !
        if (dim == "ndi") then      
          !
          resto = resto/plus_grand
          !
        else if (dim == "n+2") then
          !
          f3p3m(1) = norma
          f3p3m(2) = 0._ki
          resto = resto-log(plus_grand/mu2_scale_par)*norma
          !
        else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function f3p3m (function_3p3m.f90)'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'dim = %c0'
          tab_erreur_par(2)%arg_char = dim
          call catch_exception(0)
          !
          stop
          !
        end if
        !
        f3p3m(3) = real(resto,ki)
        f3p3m(4) = aimag(resto)
        !
      end if
      !
      ! on libere la memoire
      !
      deallocate(deja_calcule)
      deallocate(resultat)
      deallocate(deja_calcule2)
      deallocate(resultat2)
      deallocate(deja_calcule_np2)
      deallocate(resultat_np2)
      deallocate(deja_calcule22)
      deallocate(resultat22)
      !
    end function f3p3m
    !
    !****f* src/integral/three_point/function_3p3m/f3p3m_c
    ! NAME
    !
    !  Function f3p3m_c
    !
    ! USAGE
    !
    !  complex_dim3 = f3p3m_c(dim,m1,m2,m3,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  It computes the same thing that the function f3p3m, but the returned
    !  value is a complex (type ki) array of rank 1 and shape 2
    !
    ! INPUTS
    !
    !  * dim -- a character (length 3), to compute in n or n+2 dimensions, 
    !    the values are "ndi", "n+2"
    !  * m1 -- a real (type ki), the first mass squared
    !  * m2 -- a real (type ki), the second mass squared
    !  * m3 -- a real (type ki), the third mass squared
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
    !  An complex (type ki) array of rank 1 and shape 2 corresponding to 
    !  the (real part,imaginary part) of the coefficient of the 1/epsilon term
    !  and the (real part,imaginary part) of the constant term. If par1 and/or par2
    !  are different from zero for dim="n+2", an error is returned.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function f3p3m_c(dim,m1,m2,m3,par1,par2,par3)
      !
      use translate
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: m1,m2,m3
      integer, intent (in) :: par1,par2,par3
      complex(ki), dimension(2) :: f3p3m_c
      !
      real(ki), dimension(4) :: res4
      !
      res4 = f3p3m(dim,m1,m2,m3,par1,par2,par3)
      call to_complex(res4,f3p3m_c)
      !
    end function f3p3m_c
    !
    !****if* src/integral/three_point/function_3p3m/a3p3m
    ! NAME
    !
    !  Function a3p3m
    !
    ! USAGE
    !
    !  real_dim2 = a3p3m(m1,m2,m3,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This recursive function implements the formula of ref 1
    !
    ! INPUTS
    !
    !  * m1 -- a real (type ki), the first mass squared
    !  * m2 -- a real (type ki), the second mass squared
    !  * m3 -- a real (type ki), the third mass squared
    !  * par1 -- an integer, the label of the third Feynman parameter
    !  * par2 -- an integer, the label of the second Feynman parameter
    !  * par3 -- an integer, the label of the first Feynman parameter
    !
    ! SIDE EFFECTS
    !
    !  This function modify the value of the local (for the module) variables:
    !  * deja_calcule, deja_calcule2, deja_calcule_np2 and deja_calcule22
    !  * resultat, resultat2, resultat_np2 and resultat22
    !
    ! RETURN VALUE
    !
    !  It returns a real (type ki) array of rank 1 and shape 2
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    recursive function a3p3m(m1,m2,m3,par1,par2,par3) result(res_3p3m)
      !
      real(ki), intent (in) :: m1,m2,m3
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(2) :: res_3p3m
      !
      integer :: j
      integer :: nb_par_loc
      integer, dimension(3) :: par_loc,par_plus
      real(ki), dimension(4) :: truc1
      real(ki), dimension(2) :: temp0
      real(ki), dimension(4) :: temp1,temp2,temp3
      real(ki), dimension(4) :: temp10,temp11,temp12
      complex(ki) :: ctemp
      integer :: ib,b_pro,b_pro_mj
      !
      b_pro = packb(s)
      !
      par_loc = (/par1,par2,par3/)
      par_plus = par_loc+1
      nb_par_loc = count(mask=par_loc/=0)
      !
      ! cas sans parametre de feynman au numerateur
      !
      if (nb_par_loc == 0) then
        ctemp = i3_3mass(m1,m2,m3)
        res_3p3m(1) = real(ctemp,ki)
        res_3p3m(2) = aimag(ctemp)
      !
      ! cas avec un parametre de feynman au numerateur
      !
      else if (nb_par_loc == 1) then
        !
        if (deja_calcule(1)) then
          !
          temp0 = resultat(1,:)
          !
        else
          !
          temp0 = a3p3m(m1,m2,m3,0,0,0)
          resultat(1,:) = temp0
          deja_calcule(1) = .true.
          !
        end if
        !
        temp1 = 0._ki
        temp2 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            !
            if (deja_calcule2(j,1)) then
              !
              truc1 = resultat2(j,1,:)
              !
            else
              !
              truc1 = f2p_ra(s_mat_p_loc,b_pro_mj) !returns real array!
              resultat2(j,1,:) = truc1
              deja_calcule2(j,1) = .true.
              !
            end if
            !
            temp1 = temp1 + b(j)*truc1
            temp2 = temp2 + invs(j,par3)*truc1
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do
        !
        res_3p3m(1) = b(par3)*(temp0(1) - temp1(3))/sumb + temp2(3)
        res_3p3m(2) = b(par3)*(temp0(2) - temp1(4))/sumb + temp2(4)
      !
      ! cas avec deux parametres de feynman au numerateur
      !
      else if (nb_par_loc == 2) then
        !
        if (deja_calcule_np2(par_plus(3))) then
          !
          temp11 = resultat_np2(par_plus(3),:)
          !
        else
          !
          temp11 = a3p3m_np2(m1,m2,m3,0,0,par3)
          resultat_np2(par_plus(3),:) = temp11
          deja_calcule_np2(par_plus(3)) = .true.
          !
        end if
        !
        temp10 = resultat_np2(1,:)
        temp3 = invs(par2,par3)*temp10
        temp1 = b(par2)*temp11
        temp1 = mult_div(-2._ki/3._ki,temp1)*3._ki
        temp2 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            !
            if (j /= par3) then
              !
              if (deja_calcule2(j,par_plus(3))) then
                !
                truc1 = resultat2(j,par_plus(3),:)
                !
              else
                !
                truc1 = f2p_ra(s_mat_p_loc,b_pro_mj,par3) !returns real array!
                resultat2(j,par_plus(3),:) = truc1
                deja_calcule2(j,par_plus(3)) = .true.
                !
              end if
              !
              temp2 = temp2 + invs(j,par2)*truc1
              !
            end if
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do
        !
        res_3p3m(1) = -temp3(3) + temp1(3) + temp2(3)
        res_3p3m(2) = -temp3(4) + temp1(4) + temp2(4)
      !
      ! cas avec trois parametres de feynman au numerateur
      !
      else
        !
        temp12 = a3p3m_np2(m1,m2,m3,0,par2,par3)
        temp10 = resultat_np2(par_plus(3),:)
        temp11 = resultat_np2(par_plus(2),:)
        temp3 = invs(par1,par2)*temp10 &
                + invs(par1,par3)*temp11
        temp1 = b(par1)*temp12
        temp1 = mult_div(-1._ki/2._ki,temp1)*4._ki
        temp2 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            !
            if ( (j /= par3) .and. (j /= par2) ) then
              !
              truc1 = resultat22(j,par_plus(2),par_plus(3),:)
              temp2 = temp2 + invs(j,par1)*truc1
              !
            end if
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do
        !
        res_3p3m(1) = -temp3(3) + temp1(3) + temp2(3)
        res_3p3m(2) = -temp3(4) + temp1(4) + temp2(4)
        !
      end if
      !
    end function a3p3m
    !
    !****if* src/integral/three_point/function_3p3m/a3p3m_np2
    ! NAME
    !
    !  Function a3p3m_np2
    !
    ! USAGE
    !
    !  real_dim4 = a3p3m_np2(m1,m2,m3,par1,par2,par3)
    !
    ! DESCRIPTION
    !
    !  This recursive function implements the formula of ref 1
    !
    ! INPUTS
    !
    !  * m1 -- a real (type ki), the first mass squared
    !  * m2 -- a real (type ki), the second mass squared
    !  * m3 -- a real (type ki), the third mass squared
    !  * par1 -- an integer, the label of the third Feynman parameter
    !  * par2 -- an integer, the label of the second Feynman parameter
    !  * par3 -- an integer, the label of the first Feynman parameter
    !
    ! SIDE EFFECTS
    !
    !  This function modify the value of the local (for the module) variables:
    !  * deja_calcule, deja_calcule2, deja_calcule_np2 and deja_calcule22
    !  * resultat, resultat2, resultat_np2 and resultat22
    !
    ! RETURN VALUE
    !
    !  It returns a real (type ki) array of rank 1 and shape 4
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    recursive function a3p3m_np2(m1,m2,m3,par1,par2,par3) result(res_3p3m_np2)
      !
      real(ki), intent (in) :: m1,m2,m3
      integer, intent (in) :: par1,par2,par3
      real(ki), dimension(4) :: res_3p3m_np2
      !
      integer :: j
      integer :: nb_par_loc
      integer, dimension(3) :: par_loc,par_plus
      real(ki), dimension(4) :: truc1,truc2
      real(ki), dimension(2) :: temp0
      real(ki), dimension(4) :: temp1,temp2,temp3
      real(ki), dimension(4) :: temp10,temp11
      integer :: ib,b_pro,b_pro_mj
      !
      b_pro = packb(s)
      !
      par_loc = (/par1,par2,par3/)
      par_plus = par_loc+1
      nb_par_loc = count(mask=par_loc/=0)
      !
      ! cas sans parametre de feynman au numerateur
      !
      if (nb_par_loc == 0) then
        !
        if (deja_calcule(1)) then
          !
          temp0 = resultat(1,:)
          !
        else
          !
          temp0 = a3p3m(m1,m2,m3,0,0,0)
          resultat(1,:) = temp0
          deja_calcule(1) = .true.
          !
        end if
        !
        temp1 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            !
            if (deja_calcule2(j,1)) then
              !
              truc1 = resultat2(j,1,:)
              !
            else
              !
              truc1 = f2p_ra(s_mat_p_loc,b_pro_mj) !returns real array!
              resultat2(j,1,:) = truc1
              deja_calcule2(j,1) = .true.
              !
            end if
            !
            temp1 = temp1 + b(j)*truc1
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do
        !
        res_3p3m_np2(1) = (- temp1(1))/sumb
        res_3p3m_np2(2) = (- temp1(2))/sumb
        res_3p3m_np2(3) = (temp0(1) - temp1(3))/sumb
        res_3p3m_np2(4) = (temp0(2) - temp1(4))/sumb
        res_3p3m_np2 = mult_div(1._ki,res_3p3m_np2)/2._ki
      !
      ! cas avec un parametre de feynman au numerateur
      !
      else if (nb_par_loc == 1) then
        !
        if (deja_calcule_np2(1)) then
          !
          temp10 = resultat_np2(1,:)
          !
        else
          !
          temp10 = a3p3m_np2(m1,m2,m3,0,0,0)
          resultat_np2(1,:) = temp10
          deja_calcule_np2(1) = .true.
          !
        end if
        !
        temp3 = b(par3)*temp10
        temp1 = 0._ki
        temp2 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            !
            if (deja_calcule2(j,1)) then
              !
              truc1 = resultat2(j,1,:)
              !
            else
              !
              truc1 = f2p_ra(s_mat_p_loc,b_pro_mj) !returns real array!
              resultat2(j,1,:) = truc1
              deja_calcule2(j,1) = .true.
              !
            end if
            !
            temp1 = temp1 + invs(j,par3)*truc1
            !
            if (j /= par3) then
              !
              if (deja_calcule2(j,par_plus(3))) then
                !
                truc2 = resultat2(j,par_plus(3),:)
                !
              else
                !
                truc2 = f2p_ra(s_mat_p_loc,b_pro_mj,par3) !returns real array!
                resultat2(j,par_plus(3),:) = truc2
                deja_calcule2(j,par_plus(3)) = .true.
                !
              end if
              !
              temp2 = temp2 + b(j)*truc2
              !
            end if
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do
        !
        temp1 = mult_div(2._ki/3._ki,temp1)/3._ki
        temp2 = mult_div(2._ki/3._ki,temp2)/3._ki
        res_3p3m_np2(1) = (temp3(1) + temp1(1) - temp2(1))/sumb
        res_3p3m_np2(2) = (temp3(2) + temp1(2) - temp2(2))/sumb
        res_3p3m_np2(3) = (temp3(3) + temp1(3) - temp2(3))/sumb
        res_3p3m_np2(4) = (temp3(4) + temp1(4) - temp2(4))/sumb
      !
      ! cas avec deux parametres de feynman au numerateur
      !
      else if (nb_par_loc == 2) then
        !
        temp0 = a3p3m(m1,m2,m3,0,par2,par3)
        !
        if (deja_calcule_np2(par_plus(2))) then
          !
          temp10 = resultat_np2(par_plus(2),:)
          !
        else
          !
          temp10 = a3p3m_np2(m1,m2,m3,0,0,par2)
          resultat_np2(par_plus(2),:) = temp10
          deja_calcule_np2(par_plus(2)) = .true.
          !
        end if
        !
        if (deja_calcule_np2(par_plus(3))) then
          !
          temp11 = resultat_np2(par_plus(3),:)
          !
        else
          !
          temp11 = a3p3m_np2(m1,m2,m3,0,0,par3)
          resultat_np2(par_plus(3),:) = temp11
          deja_calcule_np2(par_plus(3)) = .true.
          !
        end if
        !
        temp3 =  b(par3)*temp10 + b(par2)*temp11
        temp1 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            !
            if ( (j /= par2) .and. (j /= par3) ) then
              !
              if (deja_calcule22(j,par_plus(2),par_plus(3))) then
                !
                truc1 = resultat22(j,par_plus(2),par_plus(3),:)
                !
              else
                !
                truc1 = f2p_ra(s_mat_p_loc,b_pro_mj,par2,par3) !returns real array!
                resultat22(j,par_plus(2),par_plus(3),:) = truc1
                deja_calcule22(j,par_plus(2),par_plus(3)) = .true.
                !
              end if
              !
              temp1 = temp1 + b(j)*truc1
              !
            end if
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do
        !
        res_3p3m_np2(1) = (temp3(1) - temp1(1))/sumb
        res_3p3m_np2(2) = (temp3(2) - temp1(2))/sumb
        res_3p3m_np2(3) = (temp0(1) + temp3(3) - temp1(3))/sumb
        res_3p3m_np2(4) = (temp0(2) + temp3(4) - temp1(4))/sumb
        res_3p3m_np2 = mult_div(1._ki/2._ki,res_3p3m_np2)/4._ki
      !
      ! cas avec trois parametres de feynman au numerateur
      !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a3p3m_np2:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'no need of 3-point integrals in 6 dimension &
                          &with more than one Feynman parameter in the numerator'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'The value of Feynman parameters in argument: %d1'
        tab_erreur_par(3)%arg_int_tab = (/packb(par),4/)
        call catch_exception(0)
        !
      end if
      !
    end function a3p3m_np2
    !
    !****f* src/integral/three_point/function_3p3m/i3_3mass
    ! NAME
    !
    !  Function i3_3mass
    !
    ! USAGE
    !
    !  complex = i3_3mass(m1,m2,m3)
    !
    ! DESCRIPTION
    !
    !  This function computes the scalar three off-shell external leg three point function
    !  in n dimension
    !
    ! INPUTS
    !
    !  * m1 -- a real (type ki), the first mass squared
    !  * m2 -- a real (type ki), the second mass squared
    !  * m3 -- a real (type ki), the third mass squared
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of rat_or_tot_par 
    !  (in src/module/parametre.f90)
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
    function i3_3mass(m1,m2,m3)
      !
      real(ki), intent(in) :: m1,m2,m3
      complex(ki) :: i3_3mass
      !
      complex(ki) :: cx1,cx2
      real(ki) :: delta,sig,nsig
      real(ki) :: x1,x2
      !
      if (rat_or_tot_par%tot_selected) then
        !
        delta = m1**2+m2**2+m3**2-2._ki*m1*m2-2._ki*m1*m3-2._ki*m2*m3
        sig = sign(1._ki,m1)
        !
        if (delta >= 0._ki) then
          !
          x1 = (m1+m3-m2+sqrt(delta))/(2._ki*m1)
          x2 = (m1+m3-m2-sqrt(delta))/(2._ki*m1)
          !
          ! pour avoir une fonction symetrique en fonction des trois arguments
          ! il faut multiplier la partie imaginaire par sig
          !
          nsig = sig*sign(1._ki,m2-m3)
          sig = sig*sig
          !
          i3_3mass = 1._ki/sqrt(delta)*( 2._ki*zdilog(1._ki-1._ki/x1,-sig) &
                    + 2._ki*zdilog(1._ki-1._ki/(1._ki-x2),-sig) + pi**2/3._ki &
                    + 1._ki/2._ki*( &
                      z_log2((1._ki-x1)/x1,sig) + z_log2((1._ki-x2)/x2,-sig) &
                    - z_log2(x2/(1._ki-x1),nsig) + z_log2(x1/(1._ki-x2),-nsig) ) )
          !
        else !if (delta < 0._ki) then
          !
          cx1 = (m1+m3-m2+(-sig*i_)*sqrt(-delta))/(2._ki*m1)
          cx2 = (m1+m3-m2-(-sig*i_)*sqrt(-delta))/(2._ki*m1)
          i3_3mass = (sig*i_)/sqrt(-delta)*( 2._ki*cdilog(1._ki-1._ki/cx1) &
                    + 2._ki*cdilog(1._ki-1._ki/(1._ki-cx2)) + pi**2/3._ki &
                    + 1._ki/2._ki*( &
                      (log((1._ki-cx1)/cx1))**2 + (log((1._ki-cx2)/cx2))**2 &
                    - (log(cx2/(1._ki-cx1)))**2 + (log(cx1/(1._ki-cx2)))**2 ) )
          !
        end if
        !
      else !if (rat_or_tot_par%rat_selected) then
        !
        i3_3mass = 0._ki
        !
      end if
      !
   end function i3_3mass
    !
    !****if* src/integral/three_point/function_3p3m/eval_numer_gi
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
    !  eps_glob,s13_glob,s12_glob,s23_glob,par1_glob,par2_glob,par3_glob,dim_glob
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
      eval_numer_gi = fg(z,s13_glob,s12_glob,s23_glob,&
                      &  par1_glob,par2_glob,par3_glob,&
                      &  dim_glob)
      eval_numer_gi = eval_numer_gi*jacob
      !
    end function eval_numer_gi
    !
    !****if* src/integral/three_point/function_3p3m/fg
    ! NAME
    !
    !  Function fg
    !
    ! USAGE
    !
    !  complex = fg(z,s13,s12,s23,par1,par2,par3,dim)
    !
    ! DESCRIPTION
    !
    !  This function gives the structure of the integrand for the different cases
    !
    ! INPUTS
    !
    !  * z -- a complex (type ki), the integral variable
    !  * s13 -- a real (type ki), the first mass squared
    !  * s12 -- a real (type ki), the second mass squared
    !  * s23 -- a real (type ki), the third mass squared
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
    function fg(z,s13,s12,s23,par1,par2,par3,dim)
      !
      complex(ki), intent (in) :: z
      real(ki), intent (in) :: s13,s12,s23
      integer, intent (in) :: par1,par2,par3
      character (len=3) :: dim
      complex(ki) :: fg
      !
      integer, dimension(3) :: par
      integer :: nb_par
      complex(ki) :: c_var,d_var,h_var
      !
      par = (/par1,par2,par3/)
      nb_par = count(mask=par/=0)
      !
      c_var=z*s13+(1._ki-z)*s23
      !
      d_var=z*(1._ki-z)*s12-z*s13-(1._ki-z)*s23
      !
      h_var=z*(1._ki-z)*s12
      !
      if (dim == "ndi") then      
        if (nb_par == 0) then
          !
          fg=(log(z)+log(1._ki-z)+z_log(s12,1._ki)-log(c_var))/d_var
          !
        else if (nb_par == 1) then
          !
          select case(par3)
          !
          case(1)
            !
            fg=-z*(-d_var+c_var*log(z)+c_var*log(1._ki-z)+c_var*z_log(s12,1._ki&
              &)-c_var*log(c_var))/d_var**2
            !
          case(2)
            !
            fg=(-d_var+c_var*log(z)+c_var*log(1._ki-z)+c_var*z_log(s12,1._ki)-c&
              &_var*log(c_var))*(-1._ki+z)/d_var**2
            !
          case(3)
            !
            fg=(-log(c_var)*d_var-c_var*log(c_var)+d_var*log(z)+d_var*log(1._k&
              &i-z)+d_var*z_log(s12,1._ki)+c_var*log(z)+c_var*log(1._ki-z)+c_var&
              &*z_log(s12,1._ki)-d_var)/d_var**2
            !
          case default
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
            & 'par3 should be 1, 2 or 3 but is %d0'
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
              fg=1._ki/2._ki*z**2*(d_var**2-2._ki*c_var*d_var+2._ki*c_var**2*log(z)&
                &+2._ki*c_var**2*log(1._ki-z)+2._ki*c_var**2*z_log(s12,1._ki)-2._ki*c&
                &_var**2*log(c_var))/d_var**3
              !
            case(2)
              !
              fg=-1._ki/2._ki*z*(d_var**2-2._ki*c_var*d_var+2._ki*c_var**2*log(z)+2&
                &._ki*c_var**2*log(1._ki-z)+2._ki*c_var**2*z_log(s12,1._ki)-2._ki*c_v&
                &ar**2*log(c_var))*(-1._ki+z)/d_var**3
              !
            case(3)
              !
              fg=-1._ki/2._ki*z*(-2._ki*c_var*log(c_var)*d_var-2._ki*c_var**2*log(c&
                &_var)+2._ki*c_var*d_var*log(z)+2._ki*c_var*d_var*log(1._ki-z)+2._ki&
                &*c_var*d_var*z_log(s12,1._ki)+2._ki*c_var**2*log(z)+2._ki*c_var**2&
                &*log(1._ki-z)+2._ki*c_var**2*z_log(s12,1._ki)-d_var**2-2._ki*c_var*&
                &d_var)/d_var**3
              !
            case default
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              & 'par3 should be 1, 2 or 3 but is %d0'
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
              fg=1._ki/2._ki*(-1._ki+z)**2*(d_var**2-2._ki*c_var*d_var+2._ki*c_var**&
                &2*log(z)+2._ki*c_var**2*log(1._ki-z)+2._ki*c_var**2*z_log(s12,1._ki&
                &)-2._ki*c_var**2*log(c_var))/d_var**3
              !
            case(3)
              !
              fg=1._ki/2._ki*(-2._ki*c_var*log(c_var)*d_var-2._ki*c_var**2*log(c_va&
                &r)+2._ki*c_var*d_var*log(z)+2._ki*c_var*d_var*log(1._ki-z)+2._ki*c_&
                &var*d_var*z_log(s12,1._ki)+2._ki*c_var**2*log(z)+2._ki*c_var**2*lo&
                &g(1._ki-z)+2._ki*c_var**2*z_log(s12,1._ki)-d_var**2-2._ki*c_var*d_v&
                &ar)*(-1._ki+z)/d_var**3
              !
            case default
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              & 'par3 should be 2 or 3 but is %d0'
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
              fg=1._ki/2._ki*(-2._ki*log(c_var)*d_var**2-4._ki*c_var*log(c_var)*d_v&
                &ar-2._ki*c_var**2*log(c_var)+2._ki*d_var**2*log(z)+2._ki*d_var**2*&
                &log(1._ki-z)+2._ki*d_var**2*z_log(s12,1._ki)+4._ki*c_var*d_var*log(&
                &z)+4._ki*c_var*d_var*log(1._ki-z)+4._ki*c_var*d_var*z_log(s12,1._ki&
                &)+2._ki*c_var**2*log(z)+2._ki*c_var**2*log(1._ki-z)+2._ki*c_var**2*&
                &z_log(s12,1._ki)-3._ki*d_var**2-2._ki*c_var*d_var)/d_var**3
              !
            case default
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              & 'par3 should be 3 but is %d0'
              tab_erreur_par(2)%arg_int = par3
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
              & 'par2 should be 1, 2 or 3 but is %d0'
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
                fg=-1._ki/6._ki*z**3*(-2._ki*d_var**3+3._ki*c_var*d_var**2-6._ki*c_var&
                  &**2*d_var+6._ki*c_var**3*log(z)+6._ki*c_var**3*log(1._ki-z)+6._ki*c&
                  &_var**3*z_log(s12,1._ki)-6._ki*c_var**3*log(c_var))/d_var**4
                !
              case(2)
                !
                fg=1._ki/6._ki*z**2*(-2._ki*d_var**3+3._ki*c_var*d_var**2-6._ki*c_var*&
                  &*2*d_var+6._ki*c_var**3*log(z)+6._ki*c_var**3*log(1._ki-z)+6._ki*c_&
                  &var**3*z_log(s12,1._ki)-6._ki*c_var**3*log(c_var))*(-1._ki+z)/d_va&
                  &r**4
                !
              case(3)
                !
                fg=1._ki/6._ki*z**2*(-6._ki*c_var**2*log(c_var)*d_var-6._ki*c_var**3*&
                  &log(c_var)+6._ki*c_var**2*d_var*log(z)+6._ki*c_var**2*d_var*log(1&
                  &._ki-z)+6._ki*c_var**2*d_var*z_log(s12,1._ki)+6._ki*c_var**3*log(z)&
                  &+6._ki*c_var**3*log(1._ki-z)+6._ki*c_var**3*z_log(s12,1._ki)+d_var*&
                  &*3-3._ki*c_var*d_var**2-6._ki*c_var**2*d_var)/d_var**4
                !
              case default
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                & 'par3 should be 1, 2 or 3 but is %d0'
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
                fg=-1._ki/6._ki*z*(-1._ki+z)**2*(-2._ki*d_var**3+3._ki*c_var*d_var**2-&
                  &6._ki*c_var**2*d_var+6._ki*c_var**3*log(z)+6._ki*c_var**3*log(1._ki&
                  &-z)+6._ki*c_var**3*z_log(s12,1._ki)-6._ki*c_var**3*log(c_var))/d_v&
                  &ar**4
                !
              case(3)
              !
                fg=-1._ki/6._ki*z*(-6._ki*c_var**2*log(c_var)*d_var-6._ki*c_var**3*lo&
                  &g(c_var)+6._ki*c_var**2*d_var*log(z)+6._ki*c_var**2*d_var*log(1._k&
                  &i-z)+6._ki*c_var**2*d_var*z_log(s12,1._ki)+6._ki*c_var**3*log(z)+6&
                  &._ki*c_var**3*log(1._ki-z)+6._ki*c_var**3*z_log(s12,1._ki)+d_var**3&
                  &-3._ki*c_var*d_var**2-6._ki*c_var**2*d_var)*(-1._ki+z)/d_var**4
                !
              case default
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                & 'par3 should be 2 or 3 but is %d0'
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
                fg=-1._ki/6._ki*z*(-6._ki*c_var*log(c_var)*d_var**2-12._ki*c_var**2*l&
                  &og(c_var)*d_var-6._ki*c_var**3*log(c_var)+6._ki*c_var*d_var**2*lo&
                  &g(z)+6._ki*c_var*d_var**2*log(1._ki-z)+6._ki*c_var*d_var**2*z_log(&
                  &s12,1._ki)+12._ki*c_var**2*d_var*log(z)+12._ki*c_var**2*d_var*log(&
                  &1._ki-z)+12._ki*c_var**2*d_var*z_log(s12,1._ki)+6._ki*c_var**3*log(&
                  &z)+6._ki*c_var**3*log(1._ki-z)+6._ki*c_var**3*z_log(s12,1._ki)-2._ki&
                  &*d_var**3-9._ki*c_var*d_var**2-6._ki*c_var**2*d_var)/d_var**4
                !
              case default
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                & 'par3 should be 3 but is %d0'
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              & 'par2 should be 3 but is %d0'
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
                fg=1._ki/6._ki*(-1._ki+z)**3*(-2._ki*d_var**3+3._ki*c_var*d_var**2-6._k&
                  &i*c_var**2*d_var+6._ki*c_var**3*log(z)+6._ki*c_var**3*log(1._ki-z)&
                  &+6._ki*c_var**3*z_log(s12,1._ki)-6._ki*c_var**3*log(c_var))/d_var*&
                  &*4
                !
              case(3)
                !
                fg=1._ki/6._ki*(-1._ki+z)**2*(-6._ki*c_var**2*log(c_var)*d_var-6._ki*c&
                  &_var**3*log(c_var)+6._ki*c_var**2*d_var*log(z)+6._ki*c_var**2*d_v&
                  &ar*log(1._ki-z)+6._ki*c_var**2*d_var*z_log(s12,1._ki)+6._ki*c_var**&
                  &3*log(z)+6._ki*c_var**3*log(1._ki-z)+6._ki*c_var**3*z_log(s12,1._ki&
                  &)+d_var**3-3._ki*c_var*d_var**2-6._ki*c_var**2*d_var)/d_var**4
                !
              case default
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                & 'par3 should be 2 or 3 but is %d0'
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
                fg=1._ki/6._ki*(-6._ki*c_var*log(c_var)*d_var**2-12._ki*c_var**2*log(&
                  &c_var)*d_var-6._ki*c_var**3*log(c_var)+6._ki*c_var*d_var**2*log(z&
                  &)+6._ki*c_var*d_var**2*log(1._ki-z)+6._ki*c_var*d_var**2*z_log(s12&
                  &,1._ki)+12._ki*c_var**2*d_var*log(z)+12._ki*c_var**2*d_var*log(1._k&
                  &i-z)+12._ki*c_var**2*d_var*z_log(s12,1._ki)+6._ki*c_var**3*log(z)+&
                  &6._ki*c_var**3*log(1._ki-z)+6._ki*c_var**3*z_log(s12,1._ki)-2._ki*d_&
                  &var**3-9._ki*c_var*d_var**2-6._ki*c_var**2*d_var)*(-1._ki+z)/d_var&
                  &**4
                !
              case default
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                & 'par3 should be 3 but is %d0'
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              & 'par2 should be 2 or 3 but is %d0'
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
                fg=1._ki/6._ki*(-6._ki*log(c_var)*d_var**3-18._ki*c_var*log(c_var)*d_&
                  &var**2-18._ki*c_var**2*log(c_var)*d_var-6._ki*c_var**3*log(c_var)&
                  &+6._ki*d_var**3*log(z)+6._ki*d_var**3*log(1._ki-z)+6._ki*d_var**3*z&
                  &_log(s12,1._ki)+18._ki*c_var*d_var**2*log(z)+18._ki*c_var*d_var**2&
                  &*log(1._ki-z)+18._ki*c_var*d_var**2*z_log(s12,1._ki)+18._ki*c_var**&
                  &2*d_var*log(z)+18._ki*c_var**2*d_var*log(1._ki-z)+18._ki*c_var**2*&
                  &d_var*z_log(s12,1._ki)+6._ki*c_var**3*log(z)+6._ki*c_var**3*log(1.&
                  &_ki-z)+6._ki*c_var**3*z_log(s12,1._ki)-11._ki*d_var**3-15._ki*c_var*&
                  &d_var**2-6._ki*c_var**2*d_var)/d_var**4
                !
              case default
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                & 'par3 should be 3 but is %d0'
                tab_erreur_par(2)%arg_int = par3
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              & 'par2 should be 3 but is %d0'
              tab_erreur_par(2)%arg_int = par2
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
            & 'par1 should be 1, 2 or 3 but is %d0'
            tab_erreur_par(2)%arg_int = par1
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
          & 'Unexpected value for nb_par = %d0'
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
          fg=1._ki/2._ki*(c_var**2*log(-c_var)+d_var**2*log(z)+d_var**2*log(1&
            &._ki-z)+d_var**2*z_log(-s12,-1._ki)-c_var**2*log(z)-c_var**2*log(&
            &1._ki-z)-c_var**2*z_log(-s12,-1._ki)+c_var*d_var-d_var**2)/d_var*&
            &*2
          !
        else if (nb_par == 1) then
          !
          select case(par3)
          !
          case(1)
            !
            fg=-1._ki/18._ki*z*(6._ki*c_var**3*log(-c_var)-6._ki*d_var**3*log(z)-&
              &6._ki*d_var**3*log(1._ki-z)-6._ki*d_var**3*z_log(-s12,-1._ki)-6._ki*&
              &c_var**3*log(z)-6._ki*c_var**3*log(1._ki-z)-6._ki*c_var**3*z_log(-&
              &s12,-1._ki)+4._ki*d_var**3+6._ki*c_var**2*d_var-3._ki*c_var*d_var**&
              &2)/d_var**3
            !
          case(2)
            !
            fg=1._ki/18._ki*(6._ki*c_var**3*log(-c_var)-6._ki*d_var**3*log(z)-6._k&
              &i*d_var**3*log(1._ki-z)-6._ki*d_var**3*z_log(-s12,-1._ki)-6._ki*c_v&
              &ar**3*log(z)-6._ki*c_var**3*log(1._ki-z)-6._ki*c_var**3*z_log(-s12&
              &,-1._ki)+4._ki*d_var**3+6._ki*c_var**2*d_var-3._ki*c_var*d_var**2)*&
              &(-1._ki+z)/d_var**3
            !
          case(3)
            !
            fg=1._ki/18._ki*(9._ki*c_var**2*log(-c_var)*d_var+6._ki*c_var**3*log(&
              &-c_var)-9._ki*c_var**2*d_var*log(z)-9._ki*c_var**2*d_var*log(1._ki&
              &-z)-9._ki*c_var**2*d_var*z_log(-s12,-1._ki)-6._ki*c_var**3*log(z)-&
              &6._ki*c_var**3*log(1._ki-z)-6._ki*c_var**3*z_log(-s12,-1._ki)+3._ki*&
              &d_var**3*log(z)+3._ki*d_var**3*log(1._ki-z)+3._ki*d_var**3*z_log(-&
              &s12,-1._ki)-5._ki*d_var**3+6._ki*c_var**2*d_var+6._ki*c_var*d_var**&
              &2)/d_var**3
            !
          case default
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
            & 'par3 should be 1, 2 or 3 but is %d0'
            tab_erreur_par(2)%arg_int = par3
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
          & 'Unexpected value for nb_par = %d0'
          tab_erreur_par(2)%arg_int = nb_par
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      else 
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function fg (function_3p3m.f90):'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = &
        & 'Unexpected value for dim = %c0'
        tab_erreur_par(2)%arg_char = dim
        call catch_exception(0)
        !
        stop
        !
      end if
      !
    end function fg
    !
end module function_3p3m
!
