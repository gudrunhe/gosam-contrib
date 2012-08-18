! 
!****h* src/integrals/four_point/function_4p2m_opp
! NAME
!
!  Module function_4p2m_opp
!
! USAGE
!
!  use function_4p2m_opp
!
! DESCRIPTION
!
!  This module computes the six-dimensional and eight dimensional 
!  two opposite mass four point function with or without Feynman parameters
!  in the numerator.
!
! OUTPUT
!
!  This module exports three functions f4p2m_opp, f4p2m_opp_c and f2b
!  all the other subroutines/functions of this module are private 
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
!  * generic_function_3p (src/integrals/three_point/generic_function_3p.f90)
!  * translate (src/module/translate.f90)
!
!*****
module function_4p2m_opp
  !
  use precision_golem
  use numerical_evaluation
  use dilogarithme
  use logarithme
  use constante
  use parametre
  use array
  use sortie_erreur
  use generic_function_3p
  use translate
  implicit none
  !
  private 
  !
  real(ki) :: s24_glob,s34_glob,s12_glob,s13_glob
  real(ki) :: eps_glob
  integer :: par1_glob,par2_glob,par3_glob,par4_glob
  integer :: flag_glob
  character (len=3) :: dim_glob
  !
  real(ki), dimension(4) :: b
  real(ki) :: sumb
  real(ki), dimension(4,4) :: invs,s_mat
  integer, dimension(4) :: par
  integer, dimension(4) :: s = (/1,2,3,4/)
  !
  logical, dimension(:), allocatable :: deja_calcule
  real(ki),dimension(:,:), allocatable :: resultat
  logical, dimension(:,:), allocatable :: deja_calcule3
  real(ki),dimension(:,:,:), allocatable :: resultat3
  logical, dimension(:,:), allocatable :: deja_calcule3_np2
  real(ki),dimension(:,:,:), allocatable :: resultat3_np2
  logical, dimension(:,:,:), allocatable :: deja_calcule33
  real(ki),dimension(:,:,:,:), allocatable :: resultat33
  !
  !
  public :: f4p2m_opp,f2b,f4p2m_opp_c
  !
  contains
    !
    !****f* src/integrals/four_point/function_4p2m_opp/f4p2m_opp
    ! NAME
    !
    !  Function f4p2m_opp
    !
    ! USAGE
    !
    !  real_dim_4 = f4p2m_opp(dim,s24,s13,s12,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function computes the six dimensional/eight dimensional
    !  two opposit mass four point function with or without Feynman parameters 
    !  in the numerator.
    !
    ! INPUTS
    !
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           two opposit mass four point function, dim="n+4" eight dimensional
    !           two opposit mass four point function
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s34 -- a real (type ki), the S matrix element 3,4
    !  * par1 -- an integer, the label of the fourth Feynman parameter, if none, put 0
    !  * par2 -- an integer, the label of the third Feynman parameter, if none, put 0
    !  * par3 -- an integer, the label of the second Feynman parameter, if none, put 0
    !  * par4 -- an integer, the label of the first Feynman parameter, if none, put 0
    !
    !  Be careful that, in this function, the arguments par1, par2, par3 and par4
    !  are mandatory, otherwise use the generic four point function f4p_np2 (f4p_np4).
    !  
    !
    ! SIDE EFFECTS
    !
    !  No side effects
    !
    ! RETURN VALUE
    !
    !  this function returns an array of four reals (type ki) corresponding to the 
    !  real imaginary part of 1/epsilon coefficient, real, imaginary part of the 
    !  finite part (as epsilon --> 0)
    !
    ! EXAMPLE
    !
    !  If the user wants to compute:
    !  * a six dimensional two opposit mass four point function 
    !    with no Feynman parameters in the numerator:
    !    real_dim_4 = f4p2m_opp("n+2",s24,s13,s12,s34,0,0,0,0)
    !  * a eight dimensional two opposit mass four point function 
    !    with no Feynman parameters in the numerator:
    !    real_dim_4 = f4p2m_opp("n+4",s24,s13,s12,s34,0,0,0,0)
    !  * a six dimensional two opposit mass four point function 
    !    with the Feynman parameter z1 in the numerator:
    !    real_dim_4 = f4p2m_opp("n+2",s24,s13,s12,s34,0,0,0,1)
    !  * a six dimensional two opposit mass four point function 
    !    with the Feynman parameters z1^2*z2 in the numerator:
    !    real_dim_4 = f4p2m_opp("n+2",s24,s13,s12,s34,0,2,1,1)
    !
    !*****
    function f4p2m_opp(dim,s24,s13,s12,s34,par1,par2,par3,par4)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: s24,s13,s12,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(4) :: f4p2m_opp
      !
      integer :: nb_par
      real(ki) :: lamb
      real(ki) :: plus_grand
      real(ki) :: norma
      complex(ki) :: rest1,abserr1
      complex(ki) :: rest2,abserr2
      complex(ki) :: resto,abserro
      complex(ki) :: extra_imag1,extra_imag2
      real(ki) :: pole1,pole2
      complex(ki) :: residue1,residue2
      real(ki) :: t2,t3,t4,t5
      !
      par = (/par1,par2,par3,par4/)
      !
      s_mat(1,:) = (/0._ki,s12,s13,0._ki/)
      s_mat(2,:) = (/s12,0._ki,0._ki,s24/)
      s_mat(3,:) = (/s13,0._ki,0._ki,s34/)
      s_mat(4,:) = (/0._ki,s24,s34,0._ki/)
      !
      ! on redefinit la matrice S de telle facon a ce que ces elements
      ! soient entre -1 et 1
      !
      plus_grand = maxval(array=abs(s_mat))
      s_mat = s_mat/plus_grand
      !
      b(1) = (s_mat(3,4)-s_mat(2,4))/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      b(2) = (s_mat(3,4)-s_mat(1,3))/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      b(3) = (s_mat(1,2)-s_mat(2,4))/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      b(4) = (s_mat(1,2)-s_mat(1,3))/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      !
      sumb = 2._ki*(s_mat(3,4)+s_mat(1,2)-s_mat(1,3)-s_mat(2,4))&
                  &/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      !
      invs(1,1) = 0._ki
      invs(1,2) = s_mat(3,4)/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      invs(1,3) = -s_mat(2,4)/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      invs(1,4) = 0._ki
      invs(2,1) = s_mat(3,4)/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      invs(2,2) = 0._ki
      invs(2,3) = 0._ki
      invs(2,4) = -s_mat(1,3)/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      invs(3,1) = -s_mat(2,4)/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      invs(3,2) = 0._ki
      invs(3,3) = 0._ki
      invs(3,4) = s_mat(1,2)/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      invs(4,1) = 0._ki
      invs(4,2) = -s_mat(1,3)/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      invs(4,3) = s_mat(1,2)/(s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
      invs(4,4) = 0._ki
      !
      lamb = s_mat(3,4)+s_mat(1,2)-s_mat(1,3)-s_mat(2,4)
      !
      nb_par = count(mask=par/=0)
      !
      if (nb_par == 0) then
        norma = 1._ki/6._ki
      else if (nb_par == 1) then
        norma = 1._ki/24._ki
      else
        norma = 0._ki
      end if
      !
      ! memory allocation to save time in the recursion
      !
      allocate(deja_calcule(5))
      allocate(resultat(5,2))
      allocate(deja_calcule3(4,5))
      allocate(resultat3(4,5,6))
      allocate(deja_calcule3_np2(4,5))
      allocate(resultat3_np2(4,5,4))
      allocate(deja_calcule33(4,5,5))
      allocate(resultat33(4,5,5,6))
      !
      ! initialisation
      !
      deja_calcule = .false.
      resultat = 0._ki
      deja_calcule3 = .false.
      resultat3 = 0._ki
      deja_calcule3_np2 = .false.
      resultat3_np2 = 0._ki
      deja_calcule33 = .false.
      resultat33 = 0._ki
      !
      f4p2m_opp = 0._ki
      !
      if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_4p2m_opp) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f4p2m_opp (in file f4p2m_opp.f90): &
        &the flag rat to compute the rational part is on &
        &and the program reachs a region of phase space in &
        &which det(G) = 0  Becareful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to  go on, he has to &
        &reduce the value of the parameter coupure_4p2m_opp'
        call catch_exception(0)
      end if
      !
      if (abs(sumb) > coupure_4p2m_opp) then
        !
        !---#[ analytic computation:
        !
        if (dim == "n+2") then
          !
          f4p2m_opp(3:4)= a4p2m_opp_np2(s_mat(2,4),s_mat(1,3),&
                         &s_mat(1,2),s_mat(3,4),&
                         &par1,par2,par3,par4)/plus_grand
          !
        else if (dim == "n+4") then
          !
          f4p2m_opp = a4p2m_opp_np4(s_mat(2,4),s_mat(1,3),&
                      &s_mat(1,2),s_mat(3,4),&
                      &par1,par2,par3,par4)
          f4p2m_opp(3) = f4p2m_opp(3)-log(plus_grand)*norma
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = &
          & "In f4p2m_opp (function_4p2m_opp.f90)"
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
          & "Unimplemented choice: dim = %c0"
          tab_erreur_par(2)%arg_char = dim
          call catch_exception(0)
          !
          stop
          !
        end if
        !
        !---#] analytic computation:
        !
      else
        !
        !---#[ numerical computation:
        !
        dim_glob = dim
        par1_glob = par1
        par2_glob = par2
        par3_glob = par3
        par4_glob = par4
        !
        s12_glob = s_mat(1,2)
        s13_glob = s_mat(1,3)
        s24_glob = s_mat(2,4)
        s34_glob = s_mat(3,4)
        !
        t2 = (s24_glob+s13_glob-s12_glob-s34_glob)
        t3 = (s13_glob*s24_glob-s12_glob*s34_glob)
        t4 = s13_glob-s34_glob
        t5 = s12_glob-s13_glob
        !
        resto = 0._ki
        abserro = 0._ki
        !
        ! on pose z = x - i*eps*y (avec x et y > 0)
        ! z*s24+(1-z)*s34 = s34+x*(s24-s34)-i*eps*y*(s24-s34)
        ! on veut la partie imaginaire du meme signe que i*lambda
        ! => eps*(s24-s34) < 0
        !
        ! faire attention que suivant le signe de eps_glob, on tourne dans le
        ! sens des aiguilles d'une montre ou inversement
        ! eps_glob = 1, on ferme le contour vers le bas --> -2 i Pi residu
        ! eps_glob = -1, on ferme le contour vers le haut --> +2 i Pi residu
        !
        !
        num_grand_b_info_par = lamb
        denom_grand_b_info_par = (s_mat(1,2)*s_mat(3,4)-s_mat(1,3)*s_mat(2,4))
        !
        eps_glob = sign(1._ki,s34_glob-s24_glob)
        flag_glob = 1
        !
        origine_info_par = "f4p2m_opp part 1, dimension "//dim
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,rest1,abserr1)
        !
        residue1 = compute_residue(t2,t3,t4,t5)
        !
        pole1 = (s13_glob-s34_glob)/t2
        !
        if ( (pole1 >= 0._ki) .and. (pole1 <= 1._ki) &
            & .and. (eps_glob == sign(1._ki,t2)) ) then
          !
          extra_imag1 = -2._ki*i_*pi*residue1*eps_glob
          !
        else
          !
          extra_imag1 = 0._ki
          !
        end if
        !
        resto = resto + rest1 + extra_imag1
        abserro = abserro + abserr1
        !
        eps_glob = sign(1._ki,s13_glob-s12_glob)
        flag_glob = 2
        !
        origine_info_par = "f4p2m_opp part 2, dimension "//dim
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,rest2,abserr2)
        !
        ! le residue au pole pour la somme des deux parties est nul
        !
        residue2 = -residue1
        pole2 = pole1
        !
        if ( (pole2 >= 0._ki) .and. (pole2 <= 1._ki) &
            & .and. (eps_glob == sign(1._ki,t2)) ) then
          !
          extra_imag2 = -2._ki*i_*pi*residue2*eps_glob
          !
        else
          !
          extra_imag2 = 0._ki
          !
        end if
        !
        resto = resto + rest2 + extra_imag2
        abserro = abserro + abserr2
        !
        if (dim == "n+2") then      
          !
          resto = resto/plus_grand
          !
        else if (dim == "n+4") then
          !
          f4p2m_opp(1) = norma
          f4p2m_opp(2) = 0._ki
          resto = resto-log(plus_grand/mu2_scale_par)*norma
          !
        else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = &
          & "In f4p2m_opp (function_4p2m_opp.f90)"
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = &
          & "Illegal value for dim = %c0"
          tab_erreur_par(2)%arg_char = dim
          call catch_exception(0)
          !
          stop
          !
        end if
        !
        f4p2m_opp(3) = real(resto,ki)
        f4p2m_opp(4) = aimag(resto)
        !
        !---#] numerical computation:
      end if
      !
      ! on libere la memoire
      !
      deallocate(deja_calcule)
      deallocate(resultat)
      deallocate(deja_calcule3)
      deallocate(resultat3)
      deallocate(deja_calcule3_np2)
      deallocate(resultat3_np2)
      deallocate(deja_calcule33)
      deallocate(resultat33)
      !
    end function f4p2m_opp
    !
    !****f* src/integrals/four_point/function_4p2m_opp/f4p2m_opp_c
    ! NAME
    !
    !  Function f4p2m_opp_c
    !
    ! USAGE
    !
    !  complex_dim_2 = f4p2m_opp_c(dim,s24,s13,s12,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function computes the same thing that the fucntion f4p2m_opp
    !
    ! INPUTS
    !
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           two opposit mass four point function, dim="n+4" eight dimensional
    !           two opposit mass four point function
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 2,3
    !  * s34 -- a real (type ki), the S matrix element 3,4
    !  * par1 -- an integer, the label of the fourth Feynman parameter, if none, put 0
    !  * par2 -- an integer, the label of the third Feynman parameter, if none, put 0
    !  * par3 -- an integer, the label of the second Feynman parameter, if none, put 0
    !  * par4 -- an integer, the label of the first Feynman parameter, if none, put 0
    !  
    !
    ! SIDE EFFECTS
    !
    !  No side effects
    !
    ! RETURN VALUE
    !
    !  this function returns an array of two complexs (type ki) corresponding to the 
    !   1/epsilon coefficient and the finite part (as epsilon --> 0)
    !
    ! EXAMPLE
    !
    !  see function f4p2m_opp
    !
    !*****
    function f4p2m_opp_c(dim,s24,s13,s12,s34,par1,par2,par3,par4)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: s24,s13,s12,s34
      integer, intent (in) :: par1,par2,par3,par4
      complex(ki), dimension(2) :: f4p2m_opp_c
      !
      real(ki), dimension(4) :: res4
      !
      res4 = f4p2m_opp(dim,s24,s13,s12,s34,par1,par2,par3,par4)
      call to_complex(res4,f4p2m_opp_c)
      !
    end function f4p2m_opp_c
    !
    !****if* src/integrals/four_point/function_4p2m_opp/a4p2m_opp_np2
    ! NAME
    !
    !  recursive function a4p2m_opp_np2
    !
    ! USAGE
    !
    !  real_dim_2 = a4p2m_opp_np2(s24,s13,s12,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function is the core for the analytic computation of the six dimensional
    !  two opposit mass four point function. It is recursive and implement the formulae
    !  of JHEP 10 (2005) 015.
    !
    !
    ! INPUTS
    !
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s34 -- a real (type ki), the S matrix element 3,4
    !  * par1 -- an integer, the label of the fourth Feynman parameter, if none, put 0
    !  * par2 -- an integer, the label of the third Feynman parameter, if none, put 0
    !  * par3 -- an integer, the label of the second Feynman parameter, if none, put 0
    !  * par4 -- an integer, the label of the first Feynman parameter, if none, put 0
    !  
    !
    ! SIDE EFFECTS
    !
    !  No side effects
    !
    ! RETURN VALUE
    !
    !  this function returns an array of two reals (type ki) corresponding to the 
    !  real and imaginary part of the finite part (as epsilon --> 0)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    recursive function a4p2m_opp_np2(s24,s13,s12,s34,par1,par2,par3,par4) result(res_4p2m_opp_np2)
      !
      real(ki), intent (in) :: s24,s13,s12,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(2) :: res_4p2m_opp_np2
      !
      integer, dimension(3) :: smj
      integer :: j
      integer :: nb_par_loc
      integer, dimension(4) :: par_loc,par_plus
      real(ki), dimension(6) :: truc1,truc2,truc3
      real(ki), dimension(2) :: temp0
      real(ki), dimension(6) :: temp1,temp2,temp3,temp4
      real(ki), dimension(2) :: temp10,temp11,temp12,temp13,temp14,temp15
      complex(ki) :: ctemp
      integer :: ib,b_pro,b_pro_mj
      !
      b_pro = packb(s)
      !
      par_loc = (/par1,par2,par3,par4/)
      par_plus = par_loc+1
      nb_par_loc = count(mask=par_loc/=0)
      !
      ! cas sans parametre de feynman au numerateur
      !
      if (nb_par_loc == 0) then
        !
        ctemp = -1._ki*f2b(s24,s13,s12,s34)/(s24+s13-s12-s34)
        res_4p2m_opp_np2(1) = real(ctemp,ki)
        res_4p2m_opp_np2(2) = aimag(ctemp)
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
          temp0 = a4p2m_opp_np2(s24,s13,s12,s34,0,0,0,0)
          resultat(1,:) = temp0
          deja_calcule(1) = .true.
          !
        end if
        !
        temp0 = b(par4)*temp0
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
            smj = unpackb(b_pro_mj,countb(b_pro_mj))
            !
            if (deja_calcule3(j,1)) then
              truc1 = resultat3(j,1,:)
            else
              truc1 = f3p_sc(s_mat,smj)
              resultat3(j,1,:) = truc1
              deja_calcule3(j,1) = .true.
            end if
            !
            temp1 = temp1 + invs(j,par4)*truc1/2._ki
            !
            if (j /= par4) then
              !
              if (deja_calcule3(j,par_plus(4))) then
                !
                truc2 = resultat3(j,par_plus(4),:)
                !
              else
                !
                truc2 = f3p_sc(s_mat,smj,locateb(par4,b_pro_mj))
                resultat3(j,par_plus(4),:) = truc2
                deja_calcule3(j,par_plus(4)) = .true.
                !
              end if
              !
              temp2 = temp2 - b(j)*truc2/2._ki
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
        res_4p2m_opp_np2(1) = (temp0(1) + temp1(5) + temp2(5))/sumb
        res_4p2m_opp_np2(2) = (temp0(2) + temp1(6) + temp2(6))/sumb
      !
      ! cas avec deux parametres de feynman au numerateur
      !
      else if (nb_par_loc == 2) then
        !
        if (deja_calcule(par_plus(4))) then
          !
          temp10 = resultat(par_plus(4),:)
          !
        else
          !
          temp10 = a4p2m_opp_np2(s24,s13,s12,s34,0,0,0,par4)
          resultat(par_plus(4),:) = temp10
          deja_calcule(par_plus(4)) = .true.
          !
        end if
        !
        if (deja_calcule(par_plus(3))) then
          !
          temp11 = resultat(par_plus(3),:)
          !
        else
          !
          temp11 = a4p2m_opp_np2(s24,s13,s12,s34,0,0,0,par3)
          resultat(par_plus(3),:) = temp11
          deja_calcule(par_plus(3)) = .true.
          !
        end if
        !
        temp12 = resultat(1,:)
        temp0 = b(par3)*temp10+b(par4)*temp11 - invs(par3,par4)*temp12/2._ki
        !
        temp1 = 0._ki
        temp2 = 0._ki
        temp3 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            smj = unpackb(b_pro_mj,countb(b_pro_mj))
            !
            if (j /= par3) then
              !
              if (deja_calcule3(j,par_plus(3))) then
                truc1 = resultat3(j,par_plus(3),:)
              else
                truc1 = f3p_sc(s_mat,smj,locateb(par3,b_pro_mj))
                resultat3(j,par_plus(3),:) = truc1
                deja_calcule3(j,par_plus(3)) = .true.
              end if
              !
              temp1 = temp1 + invs(j,par4)*truc1/4._ki
              !
            end if
            !
            if (j /= par4) then
              !
              if (deja_calcule3(j,par_plus(4))) then
                !
                truc2 = resultat3(j,par_plus(4),:)
                !
              else
                !
                truc2 = f3p_sc(s_mat,smj,locateb(par4,b_pro_mj))
                resultat3(j,par_plus(4),:) = truc2
                deja_calcule3(j,par_plus(4)) = .true.
                !
              end if
              !
              temp2 = temp2 + invs(j,par3)*truc2/4._ki
              !
            end if
            !
            if ( (j /= par3) .and. (j /= par4) ) then
              !
              if (deja_calcule33(j,par_plus(3),par_plus(4))) then
                !
                truc3 = resultat33(j,par_plus(3),par_plus(4),:)
                !
               else
                !
                truc3 = f3p_sc(s_mat,smj,locateb(par3,b_pro_mj),locateb(par4,b_pro_mj))
                resultat33(j,par_plus(3),par_plus(4),:) = truc3
                deja_calcule33(j,par_plus(3),par_plus(4)) = .true.
                !
              end if
              !
              temp3 = temp3 - b(j)*truc3/2._ki
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
        res_4p2m_opp_np2(1) = (temp0(1) + temp1(5) + temp2(5) + temp3(5)) &
                                *2._ki/3._ki/sumb
        res_4p2m_opp_np2(2) = (temp0(2) + temp1(6) + temp2(6) + temp3(6)) &
                                *2._ki/3._ki/sumb
      !
      ! cas avec trois parametres de feynman au numerateur
      !
      else
        !
        temp10 = a4p2m_opp_np2(s24,s13,s12,s34,0,0,par2,par3)
        temp11 = a4p2m_opp_np2(s24,s13,s12,s34,0,0,par2,par4)
        temp12 = a4p2m_opp_np2(s24,s13,s12,s34,0,0,par3,par4)
        !
        temp13 = resultat(par_plus(4),:)
        temp14 = resultat(par_plus(3),:)
        temp15 = resultat(par_plus(2),:)
        !
        temp0 = b(par4)*temp10+b(par3)*temp11+b(par2)*temp12 &
               - ( invs(par2,par3)*temp13+invs(par2,par4)*temp14&
                  +invs(par3,par4)*temp15 )/3._ki
        !
        temp1 = 0._ki
        temp2 = 0._ki
        temp3 = 0._ki
        temp4 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            smj = unpackb(b_pro_mj,countb(b_pro_mj))
            !
            if ( (j /= par2) .and. (j /= par3) ) then
              !
              truc1 = resultat33(j,par_plus(2),par_plus(3),:)
              temp1 = temp1 + invs(j,par4)*truc1/6._ki
              !
            end if
            !
            if ( (j /= par2) .and. (j /= par4) ) then
              !
              truc2 = resultat33(j,par_plus(2),par_plus(4),:)
              temp2 = temp2 + invs(j,par3)*truc2/6._ki
              !
            end if
            !
            if ( (j /= par3) .and. (j /= par4) ) then
              !
              truc3 = resultat33(j,par_plus(3),par_plus(4),:)
              temp3 = temp3 + invs(j,par2)*truc3/6._ki
              !
            end if
            !
            if ( (j /= par2) .and. (j /= par3) .and. (j /= par4) ) then
              !
              temp4 = temp4 - b(j)*f3p_sc(s_mat,smj,locateb(par2,b_pro_mj), &
                                       locateb(par3,b_pro_mj),locateb(par4,b_pro_mj))/2._ki
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
        res_4p2m_opp_np2(1) = ( temp0(1) + temp1(5) + temp2(5) + temp3(5) &
                       + temp4(5) )/2._ki/sumb
        res_4p2m_opp_np2(2) = ( temp0(2) + temp1(6) + temp2(6) + temp3(6) &
                       + temp4(6) )/2._ki/sumb
        !
      end if
      !
    end function a4p2m_opp_np2
    !
    !****if* src/integrals/four_point/function_4p2m_opp/a4p2m_opp_np4
    ! NAME
    !
    !  recursive function a4p2m_opp_np4
    !
    ! USAGE
    !
    !  real_dim_4 = a4p2m_opp_np4(s24,s13,s12,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function is the core for the analytic computation of the eight dimensional
    !  two opposit mass four point function. It is recursive and implement the formulae
    !  of JHEP 10 (2005) 015.
    !
    !
    ! INPUTS
    !
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s34 -- a real (type ki), the S matrix element 3,4
    !  * par1 -- an integer, the label of the fourth Feynman parameter, if none, put 0
    !  * par2 -- an integer, the label of the third Feynman parameter, if none, put 0
    !  * par3 -- an integer, the label of the second Feynman parameter, if none, put 0
    !  * par4 -- an integer, the label of the first Feynman parameter, if none, put 0
    !  
    !
    ! SIDE EFFECTS
    !
    !  No side effects
    !
    ! RETURN VALUE
    !
    !  this function returns an array of four reals (type ki) corresponding to the 
    !  real imaginary part of 1/epsilon coefficient, real, imaginary part of the 
    !  finite part (as epsilon --> 0)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    recursive function a4p2m_opp_np4(s24,s13,s12,s34,par1,par2,par3,par4) result(res_4p2m_opp_np4)
      !
      real(ki), intent (in) :: s24,s13,s12,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(4) :: res_4p2m_opp_np4
      !
      integer, dimension(3) :: smj
      integer :: j
      integer :: nb_par_loc
      integer, dimension(4) :: par_loc,par_plus
      real(ki), dimension(4) :: truc1
      real(ki), dimension(2) :: temp0
      real(ki), dimension(4) :: temp1,temp2,temp3
      integer :: ib,b_pro,b_pro_mj
      !
      b_pro = packb(s)
      !
      par_loc = (/par1,par2,par3,par4/)
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
          temp0 = a4p2m_opp_np2(s24,s13,s12,s34,0,0,0,0)
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
            smj = unpackb(b_pro_mj,countb(b_pro_mj))
            !
            if (deja_calcule3_np2(j,1)) then
              !
              truc1 = resultat3_np2(j,1,:)
              !
            else
              !
              truc1 = f3p_np2_sc(s_mat,smj)
              resultat3_np2(j,1,:) = truc1
              deja_calcule3_np2(j,1) = .true.
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
        res_4p2m_opp_np4(1) = (-temp1(1))/(3._ki*sumb)
        res_4p2m_opp_np4(2) = (-temp1(2))/(3._ki*sumb)
        res_4p2m_opp_np4(3) = (temp0(1)-temp1(3)-2._ki/3._ki*temp1(1))/(3._ki*sumb)
        res_4p2m_opp_np4(4) = (temp0(2)-temp1(4)-2._ki/3._ki*temp1(2))/(3._ki*sumb)
      !
      ! cas avec un parametre de feynman au numerateur
      !
      else if (nb_par_loc == 1) then
        !
        temp0 = a4p2m_opp_np2(s24,s13,s12,s34,0,0,0,par4)/3._ki
        temp1 = b(par4)*a4p2m_opp_np4(s24,s13,s12,s34,0,0,0,0)
        temp2 = 0._ki
        temp3 = 0._ki
        !
        ib = b_pro
        j = 0
        !
        do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            smj = unpackb(b_pro_mj,countb(b_pro_mj))
            !
            truc1 = resultat3_np2(j,1,:)
            temp2 = temp2 + invs(j,par4)*truc1/6._ki
            !
            if (j /= par4) then
              !
              temp3 = temp3 - b(j)*f3p_np2_sc(s_mat,smj,locateb(par4,b_pro_mj))/2._ki
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
        res_4p2m_opp_np4(1) = ( temp1(1)+temp2(1)+temp3(1) )/(2._ki*sumb)
        res_4p2m_opp_np4(2) = ( temp1(2)+temp2(2)+temp3(2) )/(2._ki*sumb)
        res_4p2m_opp_np4(3) = ( temp1(3)+temp1(1)/6._ki+temp2(3)+temp2(1)/2._ki &
                           +temp3(3)+temp3(1)/2._ki+temp0(1) )/(2._ki*sumb)
        res_4p2m_opp_np4(4) = ( temp1(4)+temp1(2)/6._ki+temp2(4)+temp2(2)/2._ki &
                           +temp3(4)+temp3(2)/2._ki+temp0(2) )/(2._ki*sumb)
      !
      ! cas avec plus de un parametre de feynman au numerateur
      !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a4p2m_opp_np4:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'no need of four-point integrals in n+4 dimension &
                          &with more than one Feynman parameter in the numerator'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'The value of Feynman parameters in argument: %d1'
        tab_erreur_par(3)%arg_int_tab = (/packb(par),4/)
        call catch_exception(0)
        !
      end if
      !
    end function a4p2m_opp_np4
    !
    !****f* src/integrals/four_point/function_4p2m_opp/f2b
    ! NAME
    !
    !  function f2b
    !
    ! USAGE
    !
    !  complex = f2b(a,b,c,d)
    !
    ! DESCRIPTION
    !
    !  This function is the "finite part" of the scalar four dimensional two  
    !  opposit mass four point function. The expression has been taken in 
    !  Nucl. Phys. {\bf B615} (2001) , 385
    !
    !
    ! INPUTS
    !
    !  * a -- a real (type ki), (p1+p2)^2
    !  * b -- a real (type ki), (p2+p3)^2
    !  * c -- a real (type ki), p2^2
    !  * d -- a real (type ki), p4^2
    !  
    !
    ! SIDE EFFECTS
    !
    !  No side effects
    !  Affected by the variable rat_or_tot_par (in src/module/parametre.f90)
    !
    ! RETURN VALUE
    !
    !  this function returns a complex (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function f2b(a,b,c,d)
      !
      real(ki), intent(in) :: a,b,c,d
      complex(ki) :: f2b
      !
      if (rat_or_tot_par%tot_selected) then
        !
        f2b = -( zdilog(1._ki-c*d/(a*b),sign(un,c*d*(a+b)-a*b*(c+d))) &
                +( z_log(c*d/(a*b),sign(un,-c*d*(a+b)+a*b*(c+d))) &
                  -z_log(c/a,sign(un,a-c))-z_log(d/b,sign(un,b-d)) ) &
                  *z_log(1._ki-c*d/(a*b),sign(un,c*d*(a+b)-a*b*(c+d))) ) &
                +zdilog(1._ki-c/a,sign(un,c-a)) &
                +zdilog(1._ki-c/b,sign(un,c-b)) &
                +zdilog(1._ki-d/a,sign(un,d-a)) &
                +zdilog(1._ki-d/b,sign(un,d-b)) &
                +z_log2(a/b,sign(un,b-a))/2._ki
        !
      else !if (rat_or_tot_par%rat_selected) then
        !
        f2b = 0._ki
        !
      end if
      !
   end function f2b
    !
    !****if* src/integrals/four_point/function_4p2m_opp/eval_numer_gi
    ! NAME
    !
    !  function eval_numer_gi
    !
    ! USAGE
    !
    !  complex = eval_numer_gi(u)
    !
    ! DESCRIPTION
    !
    !  This function contains the integrand for the numerical computation in phase 
    !  space region where det(G) ~ 0
    !
    !
    ! INPUTS
    !
    !  * u -- a real (type ki), between 0 and 1
    !  
    !
    ! SIDE EFFECTS
    !
    !  No side effects
    !
    ! RETURN VALUE
    !
    !  this function returns a complex (type ki). It is called by 
    !  the routine adapt_gauss1 in the function f4p2m_opp
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
      eval_numer_gi = fg(z,s24_glob,s13_glob,s12_glob,s34_glob,&
                      &  par1_glob,par2_glob,par3_glob,par4_glob,flag_glob,&
                      &  dim_glob)
      eval_numer_gi = eval_numer_gi*jacob
      !
    end function eval_numer_gi
    !
    !****if* src/integrals/four_point/function_4p2m_opp/compute_residue
    ! NAME
    !
    !  Function compute_residue
    !
    ! USAGE
    !
    !  complex = compute_residue(t2,t3,t4,t5)
    !
    ! DESCRIPTION
    !
    !  This function computes the residue of the pole in the case where the pole
    !  is inside the contour
    !
    ! INPUTS
    !
    !  * t2 -- a real (type ki)
    !  * t3 -- a real (type ki)
    !  * t4 -- a real (type ki)
    !  * t5 -- a real (type ki)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the global variable (for this module) dim_glob
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
    function compute_residue(t2,t3,t4,t5)
      !
      real(ki), intent (in) :: t2,t3,t4,t5
      complex(ki) :: compute_residue
      !
      complex(ki) :: temp0
      integer, dimension(4) :: par
      integer :: nb_par
      !
      par = (/par1_glob,par2_glob,par3_glob,par4_glob/)
      nb_par = count(mask=par/=0)
      !
      if (dim_glob == "n+2") then
        !---#[ dim_glob == "n+2":
        if (nb_par == 0) then
          !---#[ nb_par == 0:
          !
          temp0=-z_log(1._ki/t2*t3,1._ki)/t2
          !
          !---#] nb_par == 0:
        else if (nb_par == 1) then
          !---#[ nb_par == 1:
          select case(par4_glob)
          !
          case(1)
            !
            temp0=(-1._ki/2._ki*(t5+t2)/t2*z_log(1._ki/t2*t3,1._ki)-1._ki/2._ki*t5/&
              &t2)/t2
            !
          case(2)
            !
            temp0=-1._ki/2._ki*t4/t2**2*z_log(1._ki/t2*t3,1._ki)
            !
          case(3)
            !
            temp0=-1._ki/2._ki*(t2-t4)*z_log(1._ki/t2*t3,1._ki)/t2**2
            !
          case(4)
            !
            temp0=(1._ki/2._ki*t5/t2*z_log(1._ki/t2*t3,1._ki)+1._ki/2._ki*t5/t2)/t2
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
            tab_erreur_par(2)%arg_int = par4_glob
            call catch_exception(0)
            !
            stop
            !
          end select
          !
          !---#] nb_par == 1:
        else if (nb_par == 2) then
          !---#[ nb_par == 2:
          !
          select case(par3_glob)
          !
          case(1)
            !
            select case(par4_glob)
            !
            case(1)
              !
              temp0=(-(t2+t5)**2/t2**2*z_log(1._ki/t2*t3,1._ki)/3._ki-t5*(3._ki*t5+&
                &4._ki*t2)/t2**2/6._ki)/t2
              !
            case(2)
              !
              temp0=(-(t3+t4*t5+t4*t2)/t2**2*z_log(1._ki/t2*t3,1._ki)/6._ki-t4*t5/&
                &t2**2/6._ki)/t2
              !
            case(3)
              !
              temp0=(-(-t3+t5*t2+t2**2-t4*t5-t4*t2)/t2**2*z_log(1._ki/t2*t3,1._ki&
                &)/6._ki-t5*(t2-t4)/t2**2/6._ki)/t2
              !
            case(4)
              !
              temp0=(t5*(t2+t5)/t2**2*z_log(1._ki/t2*t3,1._ki)/3._ki+t5*(3._ki*t5+2&
                &._ki*t2)/t2**2/6._ki)/t2
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
              tab_erreur_par(2)%arg_int = par4_glob
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(2)
            !
            select case(par4_glob)
            !
            case(2)
              !
              temp0=-t4**2/t2**3*z_log(1._ki/t2*t3,1._ki)/3._ki
              !
            case(3)
              !
              temp0=-t4*(t2-t4)*z_log(1._ki/t2*t3,1._ki)/t2**3/3._ki
              !
            case(4)
              !
              temp0=((t3+t4*t5)/t2**2*z_log(1._ki/t2*t3,1._ki)/6._ki+t4*t5/t2**2/6&
                &._ki)/t2
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
              tab_erreur_par(2)%arg_int = par4_glob
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(3)
            !
            select case(par4_glob)
            !
            case(3)
              !
              temp0=-(t2-t4)**2*z_log(1._ki/t2*t3,1._ki)/t2**3/3._ki
              !
            case(4)
              !
              temp0=((t5*t2-t4*t5-t3)/t2**2*z_log(1._ki/t2*t3,1._ki)/6._ki+t5*(t2-&
                &t4)/t2**2/6._ki)/t2
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
              tab_erreur_par(2)%arg_int = par4_glob
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(4)
            !
            select case(par4_glob)
            !
            case(4)
              !
              temp0=(-t5**2/t2**2*z_log(1._ki/t2*t3,1._ki)/3._ki-t5**2/t2**2/2._ki)&
                &/t2
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
              tab_erreur_par(2)%arg_int = par4_glob
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "Unexpected value for par3 = %d0"
            tab_erreur_par(2)%arg_int = par3_glob
            call catch_exception(0)
            !
            stop
            !
          end select
          !
          !---#] nb_par == 2:
        else if (nb_par == 3) then
          !---#[ nb_par == 3:
          !
          select case(par2_glob)
          !
          case(1)
            !
            select case(par3_glob)
            !
            case(1)
              !
              select case(par4_glob)
              !
              case(1)
                !
                temp0=(-(t5+t2)**3/t2**3*z_log(1._ki/t2*t3,1._ki)/4._ki-t5*(27._ki*t5&
                  &*t2+11._ki*t5**2+18._ki*t2**2)/t2**3/24._ki)/t2
                !
              case(2)
                !
                temp0=(-(t5+t2)*(t4*t5+t4*t2+2._ki*t3)/t2**3*z_log(1._ki/t2*t3,1._ki&
                  &)/12._ki-t5*(2._ki*t3+3._ki*t4*t5+4._ki*t4*t2)/t2**3/24._ki)/t2
                !
              case(3)
                !
                temp0=((t5+t2)*(t4*t5-t5*t2+t4*t2+2._ki*t3-t2**2)/t2**3*z_log(1._ki&
                  &/t2*t3,1._ki)/12._ki+t5*(4._ki*t4*t2-4._ki*t2**2+2._ki*t3+3._ki*t4*t5&
                  &-3._ki*t5*t2)/t2**3/24._ki)/t2
                !
              case(4)
                !
                temp0=(t5*(t5+t2)**2/t2**3*z_log(1._ki/t2*t3,1._ki)/4._ki+t5*(18._ki*&
                  &t5*t2+11._ki*t5**2+6._ki*t2**2)/t2**3/24._ki)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(2)
              !
              select case(par4_glob)
              !
              case(2)
                !
                temp0=(-t4*(t4*t5+t4*t2+2._ki*t3)/t2**3*z_log(1._ki/t2*t3,1._ki)/12.&
                  &_ki-t4**2*t5/t2**3/12._ki)/t2
                !
              case(3)
                !
                temp0=((-t2*t3+2._ki*t4*t3-t4*t5*t2-t4*t2**2+t4**2*t5+t4**2*t2)/t2&
                  &**3*z_log(1._ki/t2*t3,1._ki)/12._ki+t4*t5*(t4-t2)/t2**3/12._ki)/t2
                !
              case(4)
                !
                temp0=((2._ki*t3*t5+t2*t3+t4*t5**2+t4*t5*t2)/t2**3*z_log(1._ki/t2*t&
                  &3,1._ki)/12._ki+t5*(2._ki*t3+3._ki*t4*t5+2._ki*t4*t2)/t2**3/24._ki)/&
                  &t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(3)
              !
              select case(par4_glob)
              !
              case(3)
                !
                temp0=(-(t4-t2)*(t4*t5-t5*t2+t4*t2+2._ki*t3-t2**2)/t2**3*z_log(1._k&
                  &i/t2*t3,1._ki)/12._ki-(t4-t2)**2*t5/t2**3/12._ki)/t2
                !
              case(4)
                !
                temp0=(-(2._ki*t3*t5+t2*t3-t5**2*t2-t5*t2**2+t4*t5**2+t4*t5*t2)/t2&
                  &**3*z_log(1._ki/t2*t3,1._ki)/12._ki-t5*(2._ki*t3-3._ki*t5*t2-2._ki*t2&
                  &**2+3._ki*t4*t5+2._ki*t4*t2)/t2**3/24._ki)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(par4_glob)
              !
              case(4)
                !
                temp0=(-t5**2*(t5+t2)/t2**3*z_log(1._ki/t2*t3,1._ki)/4._ki-t5**2*(11&
                  &._ki*t5+9._ki*t2)/t2**3/24._ki)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "Unexpected value for par3 = %d0"
              tab_erreur_par(2)%arg_int = par3_glob
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(2)
            !
            select case(par3_glob)
            !
            case(2)
              !
              select case(par4_glob)
              !
              case(2)
                !
                temp0=-t4**3/t2**4*z_log(1._ki/t2*t3,1._ki)/4._ki
                !
              case(3)
                !
                temp0=t4**2*(t4-t2)*z_log(1._ki/t2*t3,1._ki)/t2**4/4._ki
                !
              case(4)
                !
                temp0=(t4*(2._ki*t3+t4*t5)/t2**3*z_log(1._ki/t2*t3,1._ki)/12._ki+t4**&
                  &2*t5/t2**3/12._ki)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(3)
              !
              select case(par4_glob)
              !
              case(3)
                !
                temp0=-t4*(t4-t2)**2*z_log(1._ki/t2*t3,1._ki)/t2**4/4._ki
                !
              case(4)
                !
                temp0=(-(-t2*t3+2._ki*t4*t3-t4*t5*t2+t4**2*t5)/t2**3*z_log(1._ki/t2&
                  &*t3,1._ki)/12._ki-t4*t5*(t4-t2)/t2**3/12._ki)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(par4_glob)
              !
              case(4)
                !
                temp0=(-t5*(2._ki*t3+t4*t5)/t2**3*z_log(1._ki/t2*t3,1._ki)/12._ki-t5*&
                  &(3._ki*t4*t5+2._ki*t3)/t2**3/24._ki)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "Unexpected value for par3 = %d0"
              tab_erreur_par(2)%arg_int = par3_glob
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(3)
            !
            select case(par3_glob)
            !
            case(3)
              !
              select case(par4_glob)
              !
              case(3)
                !
                temp0=(t4-t2)**3*z_log(1._ki/t2*t3,1._ki)/t2**4/4._ki
                !
              case(4)
                !
                temp0=((t4-t2)*(t4*t5+2._ki*t3-t5*t2)/t2**3*z_log(1._ki/t2*t3,1._ki)&
                  &/12._ki+(t4-t2)**2*t5/t2**3/12._ki)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(par4_glob)
              !
              case(4)
                !
                temp0=(t5*(t4*t5+2._ki*t3-t5*t2)/t2**3*z_log(1._ki/t2*t3,1._ki)/12._k&
                  &i+t5*(-3._ki*t5*t2+3._ki*t4*t5+2._ki*t3)/t2**3/24._ki)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "Unexpected value for par3 = %d0"
              tab_erreur_par(2)%arg_int = par3_glob
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(4)
            !
            select case(par3_glob)
            !
            case(4)
              !
              select case(par4_glob)
              !
              case(4)
                !
                temp0=(t5**3/t2**3*z_log(1._ki/t2*t3,1._ki)/4._ki+11._ki/24._ki*t5**3/&
                  &t2**3)/t2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
                tab_erreur_par(2)%arg_int = par4_glob
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "Unexpected value for par3 = %d0"
              tab_erreur_par(2)%arg_int = par3_glob
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "Unexpected value for par2 = %d0"
            tab_erreur_par(2)%arg_int = par2_glob
            call catch_exception(0)
            !
            stop
            !
          end select
          !
          !---#] nb_par == 3:
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = "Unexpected value for nb_par = %d0"
          tab_erreur_par(2)%arg_int = nb_par
          call catch_exception(0)
          !
          stop
          !
        end if
        !
        !---#] dim_glob == "n+2":
      else if (dim_glob == "n+4") then
        !---#[ dim_glob == "n+4":
        !
        if (nb_par == 0) then
          !---#[ nb_par == 0:
          temp0=-1._ki/6._ki/t2**2*t3*z_log(-1._ki/t2*t3,-1._ki)
          !---#] nb_par == 0:
        else if (nb_par == 1) then
          !---#[ nb_par == 1:
          select case(par4_glob)
          !
          case(1)
            !
            temp0=(-t3*(t5+t2)/t2**2*z_log(-1._ki/t2*t3,-1._ki)/12._ki-t3*t5/t2*&
              &*2/24._ki)/t2
            !
          case(2)
            !
            temp0=-t4/t2**3*t3*z_log(-1._ki/t2*t3,-1._ki)/12._ki
            !
          case(3)
            !
            temp0=-t3*(t2-t4)*z_log(-1._ki/t2*t3,-1._ki)/t2**3/12._ki
            !
          case(4)
            !
            temp0=(t3*t5/t2**2*z_log(-1._ki/t2*t3,-1._ki)/12._ki+t3*t5/t2**2/24.&
              &_ki)/t2
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "Unexpected value for par4 = %d0"
            tab_erreur_par(2)%arg_int = par4_glob
            call catch_exception(0)
            !
            stop
            !
          end select
          !
          !---#] nb_par == 1:
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = "Unexpected value for nb_par = %d0"
          tab_erreur_par(2)%arg_int = nb_par
          call catch_exception(0)
          !
          stop
          !
        end if
        !
        !---#] dim_glob == "n+4":
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = "In compute_residue (function_4p2m_opp.f90)"
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = "Unexpected value for dim = %c0"
        tab_erreur_par(2)%arg_char = dim_glob
        call catch_exception(0)
        !
        stop
        !
      end if
      !
      compute_residue = temp0
      !
    end function compute_residue
    !
    !****if* src/integrals/four_point/function_4p2m_opp/fg
    ! NAME
    !
    !  function fg
    !
    ! USAGE
    !
    !  complex = fg(z,s24,s13,s12,s34,par1,par2,par3,par4,flag,dim)
    !
    ! DESCRIPTION
    !
    !  This function contains the one dimensional integral representation of 
    !  the six/eight dimensional two opposit mass four point function
    !
    !
    ! INPUTS
    !
    !  * z -- a real (type ki), integration variable
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s34 -- a real (type ki), the S matrix element 3,4
    !  * par1 -- an integer, the label of the fourth Feynman parameter, if none, put 0
    !  * par2 -- an integer, the label of the third Feynman parameter, if none, put 0
    !  * par3 -- an integer, the label of the second Feynman parameter, if none, put 0
    !  * par4 -- an integer, the label of the first Feynman parameter, if none, put 0
    !  * flag -- TODO undocumented parameter
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           two opposit mass four point function, dim="n+4" eight dimensional
    !           two opposit mass four point function
    !  
    !
    ! SIDE EFFECTS
    !
    !  No side effects
    !
    ! RETURN VALUE
    !
    !  this function returns a complex (type ki) corresponding to the 
    !  one dimensional integral representation of the six/eight dimensional
    !  one/zero mass four point function
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function fg(z,s24,s13,s12,s34,par1,par2,par3,par4,flag,dim)
      !
      complex(ki), intent (in) :: z
      real(ki), intent (in) :: s24,s13,s12,s34
      integer, intent (in) :: par1,par2,par3,par4,flag
      character (len=3) :: dim
      complex(ki) :: fg
      !
      integer, dimension(4) :: par
      integer :: nb_par
      complex(ki) :: c_var,e_var,f_var
      !
      par = (/par1,par2,par3,par4/)
      nb_par = count(mask=par/=0)
      !
      c_var = z*s12+(1._ki-z)*s13
      !
      f_var = z*(s24-s12)+(1._ki-z)*(s34-s13)
      !
      e_var = z*s24+(1._ki-z)*s34
      !
      if (dim == "n+2") then      
        if (nb_par == 0) then
          !
          select case(flag)
          !
          case(1)
            !
            fg=-log(e_var)/f_var
            !
          case(2)
            !
            fg=1._ki/f_var*log(c_var)
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In function fb (function_4p2m_opp.f90)"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
            &"Parameter flag should be 1 or 2 but is %d0"
            tab_erreur_par(2)%arg_int = flag
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else if (nb_par == 1) then
          !
          select case(par4)
          !
          case(1)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-1._ki/2._ki*e_var*log(e_var)/f_var**2+1._ki/2._ki/f_var
              !
            case(2)
              !
              fg=1._ki/2._ki*log(c_var)/f_var**2*e_var
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = &
              &"In function fb (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              &"Parameter flag should be 1 or 2 but is %d0"
              tab_erreur_par(2)%arg_int = flag
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(2)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-1._ki/2._ki*z*log(e_var)/f_var
              !
            case(2)
              !
              fg=1._ki/2._ki*z/f_var*log(c_var)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = &
              &"In function fb (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              &"Parameter flag should be 1 or 2 but is %d0"
              tab_erreur_par(2)%arg_int = flag
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(3)
            !
            select case(flag)
            !
            case(1)
              !
              fg=1._ki/2._ki*(-1._ki+z)*log(e_var)/f_var
              !
            case(2)
              !
              fg=-1._ki/2._ki*(-1._ki+z)/f_var*log(c_var)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = &
              &"In function fb (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              &"Parameter flag should be 1 or 2 but is %d0"
              tab_erreur_par(2)%arg_int = flag
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(4)
            !
            select case(flag)
            !
            case(1)
              !
              fg=1._ki/2._ki*log(e_var)/f_var**2*c_var-1._ki/2._ki/f_var
              !
            case(2)
              !
              fg=-1._ki/2._ki*c_var/f_var**2*log(c_var)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = &
              &"In function fb (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              &"Parameter flag should be 1 or 2 but is %d0"
              tab_erreur_par(2)%arg_int = flag
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "par4 should be 1, 2, 3 or 4 but is %d0"
            tab_erreur_par(2)%arg_int = par4
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else if (nb_par == 2) then
          !
          select case(par3)
          !
          case(1)
          !
            select case(par4)
            !
            case(1)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/3._ki*e_var**2*log(e_var)/f_var**3-1._ki/6._ki*(c_var-3._ki*&
                  &e_var)/f_var**2
                !
              case(2)
                !
                fg=1._ki/3._ki*log(c_var)/f_var**3*e_var**2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(2)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/6._ki*z*e_var*log(e_var)/f_var**2+1._ki/6._ki*z/f_var
                !
              case(2)
                !
                fg=1._ki/6._ki*z*log(c_var)/f_var**2*e_var
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(3)
              !
              select case(flag)
              !
              case(1)
                !
                fg=1._ki/6._ki*(-1._ki+z)*e_var*log(e_var)/f_var**2-1._ki/6._ki*(-1._ki&
                  &+z)/f_var
                !
              case(2)
                !
                fg=-1._ki/6._ki*(-1._ki+z)*log(c_var)/f_var**2*e_var
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(flag)
              !
              case(1)
                !
                fg=1._ki/3._ki*e_var*log(e_var)/f_var**3*c_var-1._ki/6._ki*(c_var+e_v&
                  &ar)/f_var**2
                !
              case(2)
                !
                fg=-1._ki/3._ki*log(c_var)*c_var/f_var**3*e_var
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par4 should be 1, 2, 3 or 4 but is %d0"
              tab_erreur_par(2)%arg_int = par4
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(2)
            !
            select case(par4)
            !
            case(2)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/3._ki*z**2*log(e_var)/f_var
                !
              case(2)
                !
                fg=1._ki/3._ki*z**2/f_var*log(c_var)
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(3)
              !
              select case(flag)
              !
              case(1)
                !
                fg=1._ki/3._ki*z*(-1._ki+z)*log(e_var)/f_var
                !
              case(2)
                !
                fg=-1._ki/3._ki*z*(-1._ki+z)/f_var*log(c_var)
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(flag)
              !
              case(1)
                !
                fg=1._ki/6._ki*z*log(e_var)/f_var**2*c_var-1._ki/6._ki*z/f_var
                !
              case(2)
                !
                fg=-1._ki/6._ki*z*c_var/f_var**2*log(c_var)
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par4 should be 2, 3 or 4 but is %d0"
              tab_erreur_par(2)%arg_int = par4
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(3)
            !
            select case(par4)
            !
            case(3)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/3._ki*(-1._ki+z)**2*log(e_var)/f_var
                !
              case(2)
                !
                fg=1._ki/3._ki*(-1._ki+z)**2/f_var*log(c_var)
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/6._ki*(-1._ki+z)*log(e_var)/f_var**2*c_var+1._ki/6._ki*(-1._k&
                  &i+z)/f_var
                !
              case(2)
                !
                fg=1._ki/6._ki*c_var*(-1._ki+z)/f_var**2*log(c_var)
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par4 should be 3 or 4 but is %d0"
              tab_erreur_par(2)%arg_int = par4
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(4)
            !
            select case(par4)
            !
            case(4)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/3._ki/e_var**2*log(e_var)/f_var**3*(e_var*c_var**2*f_var+&
                  &c_var**2*f_var**2+e_var**2*c_var**2+2._ki*c_var*f_var**3-f_var**&
                  &2*e_var*c_var-e_var**2*c_var*f_var-2._ki*f_var**3*e_var+f_var**4&
                  &+f_var**2*e_var**2)+1._ki/6._ki*(-e_var+3._ki*c_var)/f_var**2
                !
              case(2)
                !
                fg=1._ki/3._ki*c_var**2/f_var**3*log(c_var)
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = &
                &"In function fb (function_4p2m_opp.f90)"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = &
                &"Parameter flag should be 1 or 2 but is %d0"
                tab_erreur_par(2)%arg_int = flag
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par4 should be 4 but is %d0"
              tab_erreur_par(2)%arg_int = par4
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "par3 should be 1, 2, 3 or 4 but is %d0"
            tab_erreur_par(2)%arg_int = par3
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else if (nb_par == 3) then
          !
          select case(par2)
          !
          case(1)
            !
            select case(par3)
            !
            case(1)
              !
              select case(par4)
              !
              case(1)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/4._ki*e_var**3*log(e_var)/f_var**4+1._ki/24._ki*(2._ki*c_var&
                    &**2-7._ki*c_var*e_var+11._ki*e_var**2)/f_var**3
                  !
                case(2)
                  !
                  fg=1._ki/4._ki*log(c_var)/f_var**4*e_var**3
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(2)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*z*e_var**2*log(e_var)/f_var**3-1._ki/24._ki*z*(c_var&
                    &-3._ki*e_var)/f_var**2
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*z*log(c_var)/f_var**3*e_var**2
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*(-1._ki+z)*e_var**2*log(e_var)/f_var**3+1._ki/24._ki*(&
                    &c_var-3._ki*e_var)*(-1._ki+z)/f_var**2
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*(-1._ki+z)*log(c_var)/f_var**3*e_var**2
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/4._ki*e_var**2*log(e_var)/f_var**4*c_var+1._ki/24._ki*(c_var&
                    &**2-5._ki*c_var*e_var-2._ki*e_var**2)/f_var**3
                  !
                case(2)
                  !
                  fg=-1._ki/4._ki*log(c_var)*c_var/f_var**4*e_var**2
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 1, 2, 3 or 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(2)
              !
              select case(par4)
              !
              case(2)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*z**2*e_var*log(e_var)/f_var**2+1._ki/12._ki*z**2/f_v&
                    &ar
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*z**2*log(c_var)/f_var**2*e_var
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*z*(-1._ki+z)*e_var*log(e_var)/f_var**2-1._ki/12._ki*z*&
                    &(-1._ki+z)/f_var
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*z*(-1._ki+z)*log(c_var)/f_var**2*e_var
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*z*e_var*log(e_var)/f_var**3*c_var-1._ki/24._ki*z*(c_v&
                    &ar+e_var)/f_var**2
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*z*log(c_var)*c_var/f_var**3*e_var
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 2, 3 or 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(3)
              !
              select case(par4)
              !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*(-1._ki+z)**2*e_var*log(e_var)/f_var**2+1._ki/12._ki*&
                    &(-1._ki+z)**2/f_var
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*(-1._ki+z)**2*log(c_var)/f_var**2*e_var
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*(-1._ki+z)*e_var*log(e_var)/f_var**3*c_var+1._ki/24.&
                    &_ki*(c_var+e_var)*(-1._ki+z)/f_var**2
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*(-1._ki+z)*log(c_var)*c_var/f_var**3*e_var
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 3 or 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(par4)
              !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki/e_var*log(e_var)/f_var**4*(c_var**2*f_var**2+3._ki*&
                    &c_var**2*e_var**2+2._ki*c_var**2*f_var*e_var+2._ki*c_var*f_var**3&
                    &-2._ki*c_var*f_var*e_var**2+f_var**4+f_var**2*e_var**2-2._ki*f_va&
                    &r**3*e_var)+1._ki/24._ki*(-e_var**2+2._ki*c_var**2+5._ki*c_var*e_va&
                    &r)/f_var**3
                  !
                case(2)
                  !
                  fg=1._ki/4._ki*log(c_var)*c_var**2/f_var**4*e_var
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par3 should be 1, 2, 3 or 4 but is %d0"
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
              select case(par4)
              !
              case(2)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/4._ki*z**3*log(e_var)/f_var
                  !
                case(2)
                  !
                  fg=1._ki/4._ki*z**3/f_var*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/4._ki*z**2*(-1._ki+z)*log(e_var)/f_var
                  !
                case(2)
                  !
                  fg=-1._ki/4._ki*z**2*(-1._ki+z)/f_var*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*z**2*log(e_var)/f_var**2*c_var-1._ki/12._ki*z**2/f_va&
                    &r
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*z**2*c_var/f_var**2*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 2, 3 or 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(3)
              !
              select case(par4)
              !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/4._ki*z*(-1._ki+z)**2*log(e_var)/f_var
                  !
                case(2)
                  !
                  fg=1._ki/4._ki*z*(-1._ki+z)**2/f_var*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*z*(-1._ki+z)*log(e_var)/f_var**2*c_var+1._ki/12._ki*z&
                    &*(-1._ki+z)/f_var
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*z*c_var*(-1._ki+z)/f_var**2*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 3 or 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(par4)
              !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*z/e_var**2*log(e_var)/f_var**3*(c_var**2*e_var**2+&
                    &2._ki*c_var**2*f_var*e_var+3._ki*c_var**2*f_var**2+6._ki*c_var*f_v&
                    &ar**3-2._ki*c_var*f_var*e_var**2-4._ki*c_var*f_var**2*e_var+3._ki*&
                    &f_var**4-6._ki*f_var**3*e_var+3._ki*f_var**2*e_var**2)+1._ki/24._ki&
                    &*z*(3._ki*c_var-e_var)/f_var**2
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*z*c_var**2/f_var**3*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par3 should be 2, 3 or 4 but is %d0"
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
              select case(par4)
              !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/4._ki*(-1._ki+z)**3*log(e_var)/f_var
                  !
                case(2)
                  !
                  fg=-1._ki/4._ki*(-1._ki+z)**3/f_var*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*(-1._ki+z)**2*log(e_var)/f_var**2*c_var-1._ki/12._ki*(&
                    &-1._ki+z)**2/f_var
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*c_var*(-1._ki+z)**2/f_var**2*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 3 or 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case(4)
              !
              select case(par4)
              !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*(-1._ki+z)/e_var**2*log(e_var)/f_var**3*(c_var**2*e_&
                    &var**2+2._ki*c_var**2*f_var*e_var+3._ki*c_var**2*f_var**2+6._ki*c_&
                    &var*f_var**3-2._ki*c_var*f_var*e_var**2-4._ki*c_var*f_var**2*e_va&
                    &r+3._ki*f_var**4-6._ki*f_var**3*e_var+3._ki*f_var**2*e_var**2)-1._k&
                    &i/24._ki*(3._ki*c_var-e_var)*(-1._ki+z)/f_var**2
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*c_var**2*(-1._ki+z)/f_var**3*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par3 should be 4 but is %d0"
              tab_erreur_par(2)%arg_int = par3
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(4)
            !
            select case(par3)
            !
            case(4)
              !
              select case(par4)
              !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/4._ki/e_var**2*log(e_var)/f_var**4*c_var*(-f_var**2*e_var*&
                    &*2+f_var**4+c_var**2*f_var**2+2._ki*c_var*f_var**3+c_var**2*e_va&
                    &r**2)-1._ki/24._ki*(2._ki*e_var**2+11._ki*c_var**2-7._ki*c_var*e_var&
                    &)/f_var**3
                  !
                case(2)
                  !
                  fg=-1._ki/4._ki*c_var**3/f_var**4*log(c_var)
                  !
                case default
                  !
                  tab_erreur_par(1)%a_imprimer = .true.
                  tab_erreur_par(1)%chaine = &
                  &"In function fb (function_4p2m_opp.f90)"
                  tab_erreur_par(2)%a_imprimer = .true.
                  tab_erreur_par(2)%chaine = &
                  &"Parameter flag should be 1 or 2 but is %d0"
                  tab_erreur_par(2)%arg_int = flag
                  call catch_exception(0)
                  !
                  stop
                  !
                end select
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
                tab_erreur_par(2)%a_imprimer = .true.
                tab_erreur_par(2)%chaine = "par4 should be 4 but is %d0"
                tab_erreur_par(2)%arg_int = par4
                call catch_exception(0)
                !
                stop
                !
              end select
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par3 should be 4 but is %d0"
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
            tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "par2 should be 4 but is %d0"
            tab_erreur_par(2)%arg_int = par2
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = "Unexpected value for nb_bar = %d0"
          tab_erreur_par(2)%arg_int = nb_par
          call catch_exception(0)
          !
          stop
          !
        end if
        !
      else if (dim == "n+4") then
        !
        if (nb_par == 0) then
          !
          select case(flag)
          !
          case(1)
            !
            fg=-1._ki/6._ki*e_var*log(-e_var)/f_var+4._ki/9._ki
            !
          case(2)
            !
            fg=1._ki/6._ki*c_var/f_var*log(-c_var)
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = &
            &"In function fb (function_4p2m_opp.f90)"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
            &"Parameter flag should be 1 or 2 but is %d0"
            tab_erreur_par(2)%arg_int = flag
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else if (nb_par == 1) then
          !
          select case(par4)
          !
          case(1)
          !
          select case(flag)
            !
            case(1)
              !
              fg=-1._ki/24._ki*e_var**2*log(-e_var)/f_var**2-1._ki/144._ki*(13._ki*c&
                &_var-19._ki*e_var)/f_var
              !
            case(2)
              !
              fg=1._ki/24._ki*c_var*(c_var+2._ki*f_var)/f_var**2*log(-c_var)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = &
              &"In function fb (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              &"Parameter flag should be 1 or 2 but is %d0"
              tab_erreur_par(2)%arg_int = flag
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(2)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-1._ki/12._ki*z*e_var*log(-e_var)/f_var+2._ki/9._ki*z
              !
            case(2)
              !
              fg=1._ki/12._ki*z*c_var/f_var*log(-c_var)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = &
              &"In function fb (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              &"Parameter flag should be 1 or 2 but is %d0"
              tab_erreur_par(2)%arg_int = flag
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(3)
            !
            select case(flag)
            !
            case(1)
              !
              fg=1._ki/12._ki*(-1._ki+z)*e_var*log(-e_var)/f_var+2._ki/9._ki-2._ki/9.&
                &_ki*z
              !
            case(2)
              !
              fg=-1._ki/12._ki*c_var*(-1._ki+z)/f_var*log(-c_var)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = &
              &"In function fb (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              &"Parameter flag should be 1 or 2 but is %d0"
              tab_erreur_par(2)%arg_int = flag
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case(4)
            !
            select case(flag)
            !
            case(1)
              !
              fg=1._ki/24._ki*log(-e_var)/f_var**2*(-f_var*e_var+c_var*e_var)-1._k&
                &i/144._ki*(-13._ki*e_var+19._ki*c_var)/f_var
              !
            case(2)
              !
              fg=-1._ki/24._ki*c_var**2/f_var**2*log(-c_var)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = &
              &"In function fb (function_4p2m_opp.f90)"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = &
              &"Parameter flag should be 1 or 2 but is %d0"
              tab_erreur_par(2)%arg_int = flag
              call catch_exception(0)
              !
              stop
              !
            end select
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "par4 should be 1, 2, 3 or 4 but is %d0"
            tab_erreur_par(2)%arg_int = par4
            call catch_exception(0)
            !
            stop
            !
          end select
          !
       else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = "Unexpected value for nb_bar = %d0"
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
        tab_erreur_par(1)%chaine = "In fg (function_4p2m_opp.f90):"
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = "Unexpected value for dim = %c0"
        tab_erreur_par(2)%arg_char = dim
        call catch_exception(0)
        !
        stop
        !
      end if
    end function fg
    !
end module function_4p2m_opp
!
