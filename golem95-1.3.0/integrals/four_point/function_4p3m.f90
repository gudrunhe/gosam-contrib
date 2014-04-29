! 
!****h* src/integrals/four_point/function_4p3m
! NAME
!
!  Module function_4p3m
!
! USAGE
!
!  use function_4p3m
!
! DESCRIPTION
!
!  This module computes the six-dimensional and eight dimensional 
!  three mass four point function with or without Feynman parameters
!  in the numerator.
!
! OUTPUT
!
!  This module exports three functions f4p3m, f4p3m_c and f3
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
module function_4p3m
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
  real(ki) :: s23_glob,s24_glob,s34_glob,s12_glob,s13_glob
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
  public :: f4p3m,f3,f4p3m_c
  !
  contains
    !
    !****f* src/integrals/four_point/function_4p3m/f4p3m
    ! NAME
    !
    !  Function f4p3m
    !
    ! USAGE
    !
    !  real_dim_4 = f4p3m(dim,s24,s13,s12,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function computes the six dimensional/eight dimensional
    !  three mass four point function with or without Feynman parameters 
    !  in the numerator.
    !
    ! INPUTS
    !
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           three mass four point function, dim="n+4" eight dimensional
    !           three mass four point function
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s23 -- a real (type ki), the S matrix element 2,3
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
    !  * a six dimensional three mass four point function 
    !    with no Feynman parameters in the numerator:
    !    real_dim_4 = f4p3m("n+2",s24,s13,s12,s23,s34,0,0,0,0)
    !  * a eight dimensional three mass four point function 
    !    with no Feynman parameters in the numerator:
    !    real_dim_4 = f4p3m("n+4",s24,s13,s12,s23,s34,0,0,0,0)
    !  * a six dimensional three mass four point function 
    !    with the Feynman parameter z1 in the numerator:
    !    real_dim_4 = f4p3m("n+2",s24,s13,s12,s23,s34,0,0,0,1)
    !  * a six dimensional three mass four point function 
    !    with the Feynman parameters z1^2*z2 in the numerator:
    !    real_dim_4 = f4p3m("n+2",s24,s13,s12,s23,s34,0,2,1,1)
    !
    !*****
    function f4p3m(dim,s24,s13,s12,s23,s34,par1,par2,par3,par4)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: s24,s13,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(4) :: f4p3m
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
      real(ki) :: t1,t2,t3,t4,t5,sign_arg
      !
      par = (/par1,par2,par3,par4/)
      !
      s_mat(1,:) = (/0._ki,s12,s13,0._ki/)
      s_mat(2,:) = (/s12,0._ki,s23,s24/)
      s_mat(3,:) = (/s13,s23,0._ki,s34/)
      s_mat(4,:) = (/0._ki,s24,s34,0._ki/)
      ! on redefinit la matrice S de telle facon a ce que ses elements
      ! soient entre -1 et 1
      plus_grand = maxval(array=abs(s_mat))
      s_mat = s_mat/plus_grand
      !
      b(1)=-(-s_mat(1,2)*s_mat(3,4)**2+s_mat(3,4)*s_mat(1,3)*s_m&
        &at(2,4)+s_mat(2,3)*s_mat(1,3)*s_mat(2,4)-2._ki*s_mat(2,4)&
        &*s_mat(3,4)*s_mat(2,3)+s_mat(1,2)*s_mat(3,4)*s_mat(2,3)-&
        &s_mat(2,4)**2*s_mat(1,3)+s_mat(1,2)*s_mat(3,4)*s_mat(2,4&
        &))/(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3))**2
      b(2)=-(s_mat(1,3)-s_mat(3,4))/(s_mat(1,2)*s_mat(3,4)-s_mat&
        &(2,4)*s_mat(1,3))
      b(3)=(s_mat(1,2)-s_mat(2,4))/(s_mat(1,2)*s_mat(3,4)-s_mat(&
        &2,4)*s_mat(1,3))
      b(4)=(2._ki*s_mat(2,3)*s_mat(1,3)*s_mat(1,2)-s_mat(2,3)*s_m&
        &at(1,3)*s_mat(2,4)+s_mat(1,2)**2*s_mat(3,4)-s_mat(1,2)*s&
        &_mat(3,4)*s_mat(1,3)+s_mat(2,4)*s_mat(1,3)**2-s_mat(1,2)&
        &*s_mat(3,4)*s_mat(2,3)-s_mat(2,4)*s_mat(1,3)*s_mat(1,2))&
        &/(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3))**2
      !
      sumb=2._ki*(s_mat(1,2)*s_mat(3,4)**2-s_mat(3,4)*s_mat(1,3)*&
        &s_mat(2,4)-s_mat(2,3)*s_mat(1,3)*s_mat(2,4)+s_mat(2,4)*s&
        &_mat(3,4)*s_mat(2,3)-s_mat(1,2)*s_mat(3,4)*s_mat(2,3)+s_&
        &mat(2,4)**2*s_mat(1,3)-s_mat(1,2)*s_mat(3,4)*s_mat(2,4)+&
        &s_mat(2,4)*s_mat(1,3)**2-s_mat(1,2)*s_mat(3,4)*s_mat(1,3&
        &)+s_mat(1,2)**2*s_mat(3,4)-s_mat(2,4)*s_mat(1,3)*s_mat(1&
        &,2)+s_mat(2,3)*s_mat(1,3)*s_mat(1,2))/(s_mat(1,2)*s_mat(&
        &3,4)-s_mat(2,4)*s_mat(1,3))**2
      !
      invs(1,1)=2._ki*s_mat(2,4)*s_mat(3,4)*s_mat(2,3)/(s_mat(1,2&
        &)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3))**2
      invs(1,2)=s_mat(3,4)/(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_m&
        &at(1,3))
      invs(1,3)=-s_mat(2,4)/(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_&
        &mat(1,3))
      invs(1,4)=-s_mat(2,3)*(s_mat(1,2)*s_mat(3,4)+s_mat(2,4)*s_&
        &mat(1,3))/(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3))*&
        &*2
      invs(2,1) = invs(1,2)
      invs(2,2)=0._ki
      invs(2,3)=0._ki
      invs(2,4)=-s_mat(1,3)/(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_&
        &mat(1,3))
      invs(3,1) = invs(1,3)
      invs(3,2) = invs(2,3)
      invs(3,3)=0._ki
      invs(3,4)=s_mat(1,2)/(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_m&
        &at(1,3))
      invs(4,1) = invs(1,4)
      invs(4,2) = invs(2,4)
      invs(4,3) = invs(3,4)
      invs(4,4)=2._ki*s_mat(2,3)*s_mat(1,3)*s_mat(1,2)/(s_mat(1,2&
        &)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3))**2
      !
      lamb = (s_mat(1,2)*s_mat(3,4)**2-s_mat(3,4)*s_mat(1,3)*&
        &s_mat(2,4)-s_mat(2,3)*s_mat(1,3)*s_mat(2,4)+s_mat(2,4)*s&
        &_mat(3,4)*s_mat(2,3)-s_mat(1,2)*s_mat(3,4)*s_mat(2,3)+s_&
        &mat(2,4)**2*s_mat(1,3)-s_mat(1,2)*s_mat(3,4)*s_mat(2,4)+&
        &s_mat(2,4)*s_mat(1,3)**2-s_mat(1,2)*s_mat(3,4)*s_mat(1,3&
        &)+s_mat(1,2)**2*s_mat(3,4)-s_mat(2,4)*s_mat(1,3)*s_mat(1&
        &,2)+s_mat(2,3)*s_mat(1,3)*s_mat(1,2))
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
      f4p3m = 0._ki
      !
      if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_4p3m) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f4p3m (in file f4p3m.f90): &
        &the flag rat to compute the rational part is on &
        &and the program reachs a region of phase space in &
        &which det(G) = 0  Becareful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to  go on, he has to &
        &reduce the value of the parameter coupure_4p3m'
        call catch_exception(0)
        !
        stop
        !
      end if
      !
      if (abs(sumb) > coupure_4p3m) then
        !
        ! analytic computation
        !
        if (dim == "n+2") then
          !
          f4p3m(3:4)= a4p3m_np2(s_mat(2,4),s_mat(1,3),&
                         &s_mat(1,2),s_mat(2,3),s_mat(3,4),&
                         &par1,par2,par3,par4)/plus_grand
          !
        else if (dim == "n+4") then
          !
          f4p3m = a4p3m_np4(s_mat(2,4),s_mat(1,3),&
                      &s_mat(1,2),s_mat(2,3),s_mat(3,4),&
                      &par1,par2,par3,par4)
          f4p3m(3) = f4p3m(3)-log(plus_grand)*norma
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
        par4_glob = par4
        !
        s12_glob = s_mat(1,2)
        s13_glob = s_mat(1,3)
        s23_glob = s_mat(2,3)
        s24_glob = s_mat(2,4)
        s34_glob = s_mat(3,4)
        !
        t1 = (s13_glob-s34_glob)*(s24_glob-s12_glob)
        t2 = (s24_glob+s13_glob-s12_glob-s34_glob)
        t3 = (s13_glob*s24_glob-s12_glob*s34_glob)
        t4 = s13_glob-s34_glob
        t5 = s12_glob-s13_glob
        !
        sign_arg = sign(un,(t1*s23_glob-t2*t3))
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
        denom_grand_b_info_par = (s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3))**2
        !
        eps_glob = sign(1._ki,s34_glob-s24_glob)
        flag_glob = 1
        !
        origine_info_par = "f4p3m part 1, dimension "//dim
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,rest1,abserr1)
        !
        residue1 = compute_residue(t1,t2,t3,t4,t5,s23_glob,sign_arg)
        !
        pole1 = (s13_glob-s34_glob)/t2
        !
        if ( (pole1 >= 0._ki) .and. (pole1 <= 1._ki) &
            & .and. (eps_glob == sign(1._ki,t2)) ) then
          extra_imag1 = -2._ki*i_*pi*residue1*eps_glob
        else
          extra_imag1 = 0._ki
        end if
        !
        resto = resto + rest1 + extra_imag1
        abserro = abserro + abserr1
        !
        eps_glob = sign(1._ki,s13_glob-s12_glob)
        flag_glob = 2
        !
        origine_info_par = "f4p3m part 2, dimension "//dim
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,rest2,abserr2)
        ! le residue au pole pour la somme des deux parties est nul
        residue2 = -residue1
        pole2 = pole1
        !
        if ( (pole2 >= 0._ki) .and. (pole2 <= 1._ki) &
            & .and. (eps_glob == sign(1._ki,t2)) ) then
          extra_imag2 = -2._ki*i_*pi*residue2*eps_glob
        else
          extra_imag2 = 0._ki
        end if
        !
        resto = resto + rest2 + extra_imag2
        abserro = abserro + abserr2
        !
        if (dim == "n+2") then      
          resto = resto/plus_grand
        else if (dim == "n+4") then
          f4p3m(1) = norma
          f4p3m(2) = 0._ki
          resto = resto-log(plus_grand/mu2_scale_par)*norma
        end if
        !
        f4p3m(3) = real(resto,ki)
        f4p3m(4) = aimag(resto)
        !
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
    end function f4p3m
    !
    !****f* src/integrals/four_point/function_4p3m/f4p3m_c
    ! NAME
    !
    !  Function f4p3m_c
    !
    ! USAGE
    !
    !  complex_dim_2 = f4p3m_c(dim,s24,s13,s12,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function computes the same thing that the function f4p3m
    !
    ! INPUTS
    !
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           two adjacent mass four point function, dim="n+4" eight dimensional
    !           two adjacent mass four point function
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s23 -- a real (type ki), the S matrix element 2,3
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
    !  see function f4p3m
    !
    !*****
    function f4p3m_c(dim,s24,s13,s12,s23,s34,par1,par2,par3,par4)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: s24,s13,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      complex(ki), dimension(2) :: f4p3m_c
      !
      real(ki), dimension(4) :: res4
      !
      res4 = f4p3m(dim,s24,s13,s12,s23,s34,par1,par2,par3,par4)
      call to_complex(res4,f4p3m_c)
      !
    end function f4p3m_c
    !
    !****if* src/integrals/four_point/function_4p3m/a4p3m_np2
    ! NAME
    !
    !  recursive function a4p3m_np2
    !
    ! USAGE
    !
    !  real_dim_2 = a4p3m_np2(s24,s13,s12,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function is the core for the analytic computation of the six dimensional
    !  three mass four point function. It is recursive and implement the formulae
    !  of JHEP 10 (2005) 015.
    !
    !
    ! INPUTS
    !
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s23 -- a real (type ki), the S matrix element 2,3
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
    recursive function a4p3m_np2(s24,s13,s12,s23,s34,par1,par2,par3,par4) result(res_4p3m_np2)
      !
      real(ki), intent (in) :: s24,s13,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(2) :: res_4p3m_np2
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
        if (deja_calcule3(1,1)) then
          !
          truc1 = resultat3(1,1,:)
          !
        else
          !
          truc1 = f3p_sc(s_mat,unpackb(ibclr(b_pro,1),3))
          resultat3(1,1,:) = truc1
          deja_calcule3(1,1) = .true.
          !
        end if
        !
        if (deja_calcule3(4,1)) then
          !
          truc2 = resultat3(4,1,:)
          !
        else
          !
          truc2 = f3p_sc(s_mat,unpackb(ibclr(b_pro,4),3))
          resultat3(4,1,:) = truc2
          deja_calcule3(4,1) = .true.
          !
        end if
        !
        ctemp = f3(s24,s13,s12,s23,s34)
        res_4p3m_np2(1) = (real(ctemp,ki) - b(1)*truc1(5) - b(4)*truc2(5))/sumb
        res_4p3m_np2(2) = (aimag(ctemp) - b(1)*truc1(6) - b(4)*truc2(6))/sumb
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
          temp0 = a4p3m_np2(s24,s13,s12,s23,s34,0,0,0,0)
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
              !
              truc1 = resultat3(j,1,:)
              !
            else
              !
              truc1 = f3p_sc(s_mat,smj)
              resultat3(j,1,:) = truc1
              deja_calcule3(j,1) = .true.
              !
            end if
            !
            temp1 = temp1 + invs(j,par4)*truc1/2._ki
            !
            if (j /= par4) then
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
        res_4p3m_np2(1) = (temp0(1) + temp1(5) + temp2(5))/sumb
        res_4p3m_np2(2) = (temp0(2) + temp1(6) + temp2(6))/sumb
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
          temp10 = a4p3m_np2(s24,s13,s12,s23,s34,0,0,0,par4)
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
          temp11 = a4p3m_np2(s24,s13,s12,s23,s34,0,0,0,par3)
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
                !
                truc1 = resultat3(j,par_plus(3),:)
                !
              else
                !
                truc1 = f3p_sc(s_mat,smj,locateb(par3,b_pro_mj))
                resultat3(j,par_plus(3),:) = truc1
                deja_calcule3(j,par_plus(3)) = .true.
                !
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
        res_4p3m_np2(1) = (temp0(1) + temp1(5) + temp2(5) + temp3(5)) &
                                *2._ki/3._ki/sumb
        res_4p3m_np2(2) = (temp0(2) + temp1(6) + temp2(6) + temp3(6)) &
                                *2._ki/3._ki/sumb
      !
      ! cas avec trois parametres de feynman au numerateur
      !
      else
        !
        temp10 = a4p3m_np2(s24,s13,s12,s23,s34,0,0,par2,par3)
        temp11 = a4p3m_np2(s24,s13,s12,s23,s34,0,0,par2,par4)
        temp12 = a4p3m_np2(s24,s13,s12,s23,s34,0,0,par3,par4)
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
        res_4p3m_np2(1) = ( temp0(1) + temp1(5) + temp2(5) + temp3(5) &
                       + temp4(5) )/2._ki/sumb
        res_4p3m_np2(2) = ( temp0(2) + temp1(6) + temp2(6) + temp3(6) &
                       + temp4(6) )/2._ki/sumb
      end if
      !
    end function a4p3m_np2
    !
    !****if* src/integrals/four_point/function_4p3m/a4p3m_np4
    ! NAME
    !
    !  recursive function a4p3m_np4
    !
    ! USAGE
    !
    !  real_dim_4 = a4p3m_np4(s24,s13,s12,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function is the core for the analytic computation of the eight dimensional
    !  three mass four point function. It is recursive and implement the formulae
    !  of JHEP 10 (2005) 015.
    !
    !
    ! INPUTS
    !
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s23 -- a real (type ki), the S matrix element 2,3
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
    recursive function a4p3m_np4(s24,s13,s12,s23,s34,par1,par2,par3,par4) result(res_4p3m_np4)
      !
      real(ki), intent (in) :: s24,s13,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(4) :: res_4p3m_np4
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
          temp0 = a4p3m_np2(s24,s13,s12,s23,s34,0,0,0,0)
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
        res_4p3m_np4(1) = (-temp1(1))/(3._ki*sumb)
        res_4p3m_np4(2) = (-temp1(2))/(3._ki*sumb)
        res_4p3m_np4(3) = (temp0(1)-temp1(3)-2._ki/3._ki*temp1(1))/(3._ki*sumb)
        res_4p3m_np4(4) = (temp0(2)-temp1(4)-2._ki/3._ki*temp1(2))/(3._ki*sumb)
      !
      ! cas avec un parametre de feynman au numerateur
      !
      else if (nb_par_loc == 1) then
        !
        temp0 = a4p3m_np2(s24,s13,s12,s23,s34,0,0,0,par4)/3._ki
        temp1 = b(par4)*a4p3m_np4(s24,s13,s12,s23,s34,0,0,0,0)
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
        res_4p3m_np4(1) = ( temp1(1)+temp2(1)+temp3(1) )/(2._ki*sumb)
        res_4p3m_np4(2) = ( temp1(2)+temp2(2)+temp3(2) )/(2._ki*sumb)
        res_4p3m_np4(3) = ( temp1(3)+temp1(1)/6._ki+temp2(3)+temp2(1)/2._ki &
                           +temp3(3)+temp3(1)/2._ki+temp0(1) )/(2._ki*sumb)
        res_4p3m_np4(4) = ( temp1(4)+temp1(2)/6._ki+temp2(4)+temp2(2)/2._ki &
                           +temp3(4)+temp3(2)/2._ki+temp0(2) )/(2._ki*sumb)
      !
      ! cas avec plus de un parametre de feynman au numerateur
      !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a4p3m_np4:'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'no need of four-point integrals in n+4 dimension &
                          &with more than one Feynman parameter in the numerator'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'The value of Feynman parameters in argument: %d1'
        tab_erreur_par(3)%arg_int_tab = (/packb(par),4/)
        call catch_exception(0)
        !
        stop
        !
      end if
      !
    end function a4p3m_np4
    !
    !****f* src/integrals/four_point/function_4p3m/f3
    ! NAME
    !
    !  function f3
    !
    ! USAGE
    !
    !  complex = f3(s,t,m2,m3,m4)
    !
    ! DESCRIPTION
    !
    !  This function is the "finite part" of the scalar four dimensional three  
    !  mass four point function.
    !
    !
    ! INPUTS
    !
    !  * s -- a real (type ki), (p1+p2)^2
    !  * t -- a real (type ki), (p2+p3)^2
    !  * m2 -- a real (type ki), p2^2
    !  * m3 -- a real (type ki), p3^2
    !  * m4 -- a real (type ki), p4^2
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
    function f3(s,t,m2,m3,m4)
      !
      real(ki), intent(in) :: s,t,m2,m3,m4
      complex(ki) :: f3
      !
      real(ki) :: lamb1
      !
      lamb1 = s*t-m2*m4
      !
      if (rat_or_tot_par%tot_selected) then
        !
        f3 = (  z_log(m3/s,sign(un,s-m3))*z_log(m4/s,sign(un,s-m4)) &
             &- z_log(t/s,sign(un,s-t))*z_log(m2/s,sign(un,s-m2)) &
             &- z_log(t/s,sign(un,s-t))*z_log(m3/s,sign(un,s-m3)) &
             &+ z_log(m2/s,sign(un,s-m2))*z_log(m3/s,sign(un,s-m3)) &
             &- 2._ki*zdilog(1._ki-m2/s,sign(un,m2-s)) &
             &- 2._ki*zdilog(1._ki-m4/t,sign(un,m4-t)) &
             &+ 2._ki &
             & *( zdilog(1._ki-m2*m4/(s*t),sign(un,m2*m4*(s+t)-s*t*(m2+m4))) &
                  +( z_log(m2*m4/(s*t),sign(un,-m2*m4*(s+t)+s*t*(m2+m4))) &
                    -z_log(m2/s,sign(un,s-m2))-z_log(m4/t,sign(un,t-m4)) ) &
                  *z_log(1._ki-m2*m4/(s*t),sign(un,m2*m4*(s+t)-s*t*(m2+m4))) ) &
             & )/lamb1
        !
      else !if (rat_or_tot_par%rat_selected) then
        !
        f3 = 0._ki
        !
      end if
      !
    end function f3
    !
    !****if* src/integrals/four_point/function_4p3m/eval_numer_gi
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
    !  the routine adapt_gauss1 in the function f4p3m
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
      eval_numer_gi = fg(z,s24_glob,s13_glob,s12_glob,s23_glob,s34_glob,&
                      &  par1_glob,par2_glob,par3_glob,par4_glob,flag_glob,&
                      &  dim_glob)
      eval_numer_gi = eval_numer_gi*jacob
      !
    end function eval_numer_gi
    !
    !****if* src/integrals/four_point/function_4p3m/compute_residue
    ! NAME
    !
    !  Function compute_residue
    !
    ! USAGE
    !
    !  complex = compute_residue(t1,t2,t3,t4,t5,t6,sign_arg)
    !
    ! DESCRIPTION
    !
    !  This function computes the residue of the pole in the case where the pole
    !  is inside the contour
    !
    ! INPUTS
    !
    !  * t1 -- a real (type ki)
    !  * t2 -- a real (type ki)
    !  * t3 -- a real (type ki)
    !  * t4 -- a real (type ki)
    !  * t5 -- a real (type ki)
    !  * t6 -- a real (type ki)
    !  * sign_arg -- a real (type ki)
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
    function compute_residue(t1,t2,t3,t4,t5,t6,sign_arg)
      !
      real(ki), intent (in) :: t1,t2,t3,t4,t5,t6,sign_arg
      complex(ki) :: compute_residue
      !
      complex(ki) :: temp0,stemp1,stemp2,stemp3,stemp4,&
                      &stemp5,stemp6,stemp7,stemp8,stemp9,&
                      &stemp10,stemp11
      integer, dimension(4) :: par
      integer :: nb_par
      !
      par = (/par1_glob,par2_glob,par3_glob,par4_glob/)
      nb_par = count(mask=par/=0)
      !
      if (dim_glob == "n+2") then      
        if (nb_par == 0) then
          !
          temp0=(-z_log(t1*t6/t2**2,1._ki)+q(1,(-t1*t6+t2*t3)/t2/t3,sign_arg&
            &))/t2
          !
        else if (nb_par == 1) then
          !
          select case(par4_glob)
          !
          case(1)
            !
            temp0=(-1._ki/2._ki*(t2**2+t6*t2+t5*t2-2._ki*t6*t4)/t2**2*z_log(t1*t&
              &6/t2**2,1._ki)+1._ki/2._ki/t3*(t2**2*t3+t5*t2*t3+t2*t3*t6-2._ki*t5*&
              &t1*t6-2._ki*t2*t1*t6-2._ki*t3*t6*t4)/t2**2*q(2,(t2*t3-t1*t6)/t2/t&
              &3,sign_arg)-1._ki/2._ki*(2._ki*t2**2+2._ki*t5*t2+t6*t2-2._ki*t6*t4)/&
              &t2**2)/t2
            !
          case(2)
            !
            temp0=(-1._ki/2._ki*(t4*t2**3+2._ki*t2*t3*t6+2._ki*t6*t4*t5*t2-6._ki*t&
              &3*t6*t4+2._ki*t5*t1*t6-4._ki*t4**2*t5*t6)/t2**4*z_log(t1*t6/t2**2&
              &,1._ki)+1._ki/2._ki*t4/t2*q(2,(t2*t3-t1*t6)/t2/t3,sign_arg)-1._ki/2&
              &._ki*(t4*t2**2*t3-4._ki*t4**2*t2*t3+2._ki*t2*t1*t3+2._ki*t2*t4*t5*t&
              &1-6._ki*t4*t1*t3+4._ki*t4**3*t3-4._ki*t4**2*t5*t1)*t6/t1/t2**4)/t2
            !
          case(3)
            !
            temp0=(-1._ki/2._ki*(t2**4-t4*t2**3+2._ki*t6*t2**2*t5-4._ki*t2*t3*t6-&
              &6._ki*t6*t4*t5*t2+6._ki*t3*t6*t4-2._ki*t5*t1*t6+4._ki*t4**2*t5*t6)/&
              &t2**4*z_log(t1*t6/t2**2,1._ki)+1._ki/2._ki/t2*(t2-t4)*q(2,(t2*t3-t&
              &1*t6)/t2/t3,sign_arg)-1._ki/2._ki*(t3*t2**3-5._ki*t4*t2**2*t3+2._ki&
              &*t2**2*t1*t5+8._ki*t4**2*t2*t3-4._ki*t2*t1*t3-6._ki*t2*t4*t5*t1+6.&
              &_ki*t4*t1*t3-4._ki*t4**3*t3+4._ki*t4**2*t5*t1)*t6/t1/t2**4)/t2
            !
          case(4)
            !
            temp0=(1._ki/2._ki*(t5*t2+t6*t2-2._ki*t6*t4)/t2**2*z_log(t1*t6/t2**2&
              &,1._ki)-1._ki/2._ki/t3*(t5*t2*t3+t2*t3*t6-2._ki*t5*t1*t6-2._ki*t3*t6&
              &*t4)/t2**2*q(2,(t2*t3-t1*t6)/t2/t3,sign_arg)+1._ki/2._ki*(t6*t2+2&
              &._ki*t5*t2-2._ki*t6*t4)/t2**2)/t2
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = &
            & 'In function compute_residue (function_4p3m.f90):'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = &
            & 'par4 = %d0'
            tab_erreur_par(2)%arg_int = par4_glob
            call catch_exception(0)
            !
            stop
            !
          end select
          !
        else if (nb_par == 2) then
          !
          select case(par3_glob)
          !
          case(1)
            !
            select case(par4_glob)
            !
            case(1)
              !
              stemp2=-(-2._ki*t1*t6**2+4._ki*t6**2*t4**2-4._ki*t6**2*t2*t4+t6**2*t&
                &2**2+t2**3*t6+t6*t2**2*t5-2._ki*t6*t4*t5*t2-2._ki*t4*t6*t2**2-t6*&
                &t3*t2+t2**4+t5**2*t2**2+2._ki*t5*t2**3)/t2**4*z_log(t1*t6/t2**2,&
                &1._ki)/3._ki
              !
              stemp3=(-6._ki*t1*t6*t2**2*t3*t5+t3**2*t2**4-3._ki*t1*t6*t2*t3*t5**&
                &2+2._ki*t2**3*t3**2*t5+t2**2*t3**2*t5**2+t2**2*t3**2*t6*t5-t6*t2&
                &*t3**3-2._ki*t2*t3**2*t6*t4*t5+3._ki*t2**2*t1**2*t6**2-3._ki*t3*t2&
                &**3*t1*t6+3._ki*t1**2*t6**2*t5**2+t1*t6**2*t3**2+6._ki*t1**2*t6**&
                &2*t5*t2-3._ki*t1*t6**2*t2*t3*t5+6._ki*t1*t6**2*t3*t4*t5-3._ki*t1*t&
                &6**2*t2**2*t3+6._ki*t1*t6**2*t2*t3*t4-2._ki*t2**2*t3**2*t6*t4+t2*&
                &*3*t3**2*t6-4._ki*t2*t3**2*t6**2*t4+t2**2*t3**2*t6**2+4._ki*t3**2&
                &*t6**2*t4**2)/t3**2/t2**4*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/3._k&
                &i+(2._ki*t1*t6**2*t3-12._ki*t6**2*t3*t4**2-3._ki*t6**2*t2**2*t3+12&
                &._ki*t6**2*t2*t3*t4+3._ki*t6*t1*t2*t5**2+3._ki*t1*t6*t2**3+6._ki*t6&
                &*t1*t2**2*t5+3._ki*t2*t3**2*t6-5._ki*t6*t3*t2**3+10._ki*t6*t2**2*t&
                &3*t4-5._ki*t6*t5*t2**2*t3+10._ki*t6*t2*t3*t4*t5-6._ki*t2**4*t3-6._k&
                &i*t2**2*t3*t5**2-17._ki*t2**3*t3*t5)/t2**4/t3/6._ki
              !
              stemp1=stemp2+stemp3
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            case(2)
              !
              stemp2=-(-36._ki*t6**2*t1*t5*t4+12._ki*t6**2*t1*t5*t2-12._ki*t1*t6**&
                &2*t3-36._ki*t6**2*t2*t3*t4-24._ki*t6**2*t4**2*t5*t2+6._ki*t6**2*t2&
                &**2*t4*t5+24._ki*t6**2*t5*t4**3+48._ki*t6**2*t3*t4**2+6._ki*t6**2*&
                &t2**2*t3+2._ki*t1*t6*t2**3+2._ki*t6*t2**4*t4-4._ki*t6*t2**3*t4**2+&
                &t4*t2**5+t4*t2**4*t5+t2**4*t3)/t2**6*z_log(t1*t6/t2**2,1._ki)/6.&
                &_ki
              !
              stemp3=-(-2._ki*t6*t4*t2*t3+4._ki*t4**2*t6*t3+3._ki*t4*t2*t1*t6-t2**&
                &2*t3*t4-t2*t3*t4*t5+3._ki*t4*t5*t1*t6-t2*t3**2+t1*t6*t3)/t2**3/t&
                &3*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/6._ki-(-36._ki*t1**2*t6**2*t4&
                &*t5+12._ki*t1**2*t6**2*t5*t2-12._ki*t1**2*t6**2*t3+18._ki*t6**2*t1&
                &*t2**2*t4*t5-72._ki*t6**2*t1*t4**2*t5*t2+18._ki*t1*t6**2*t2**2*t3&
                &+144._ki*t1*t6**2*t4**2*t3-108._ki*t1*t6**2*t2*t3*t4+72._ki*t6**2*&
                &t1*t5*t4**3+4._ki*t6**2*t2**3*t3*t4+48._ki*t6**2*t4**3*t2*t3-24._k&
                &i*t6**2*t2**2*t3*t4**2-32._ki*t6**2*t4**4*t3-8._ki*t6*t1*t2**3*t4&
                &**2+4._ki*t6*t1*t4*t2**4+3._ki*t1*t4*t2**5+3._ki*t1*t4*t2**4*t5+t1&
                &*t2**4*t3)/t2**6/t1/12._ki
              !
              stemp1=stemp2+stemp3
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            case(3)
              !
              stemp2=(24._ki*t6**2*t1*t5*t2-36._ki*t6**2*t1*t5*t4-12._ki*t1*t6**2*&
                &t3+18._ki*t6**2*t2**2*t3-60._ki*t6**2*t2*t3*t4-48._ki*t6**2*t4**2*&
                &t5*t2-6._ki*t6**2*t2**3*t5+48._ki*t6**2*t3*t4**2+24._ki*t6**2*t5*t&
                &4**3+30._ki*t6**2*t2**2*t4*t5+2._ki*t1*t6*t2**3-2._ki*t6*t2**5+6._k&
                &i*t6*t2**4*t4-4._ki*t6*t2**3*t4**2+t2**4*t3+t4*t2**4*t5+t4*t2**5&
                &-t5*t2**5-t2**6)/t2**6*z_log(t1*t6/t2**2,1._ki)/6._ki
              !
              stemp3=(2._ki*t6*t2**2*t3-6._ki*t6*t4*t2*t3+t2**2*t3*t5+3._ki*t4*t2*&
                &t1*t6-t2**2*t3*t4+4._ki*t4**2*t6*t3-t2*t3*t4*t5+3._ki*t4*t5*t1*t6&
                &-3._ki*t5*t2*t1*t6+t2**3*t3-t2*t3**2-3._ki*t2**2*t1*t6+t1*t6*t3)/&
                &t2**3/t3*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/6._ki+(24._ki*t1**2*t6&
                &**2*t5*t2-36._ki*t1**2*t6**2*t4*t5-12._ki*t1**2*t6**2*t3-180._ki*t&
                &1*t6**2*t2*t3*t4-18._ki*t6**2*t1*t5*t2**3+72._ki*t6**2*t1*t5*t4**&
                &3-144._ki*t6**2*t1*t4**2*t5*t2+90._ki*t6**2*t1*t2**2*t4*t5+144._ki&
                &*t1*t6**2*t4**2*t3+54._ki*t1*t6**2*t2**2*t3-32._ki*t6**2*t4**4*t3&
                &+80._ki*t6**2*t4**3*t2*t3-4._ki*t6**2*t2**4*t3+28._ki*t6**2*t2**3*&
                &t3*t4-72._ki*t6**2*t2**2*t3*t4**2-8._ki*t6*t1*t2**3*t4**2+12._ki*t&
                &6*t1*t4*t2**4-4._ki*t6*t1*t2**5+3._ki*t1*t4*t2**4*t5+t1*t2**4*t3+&
                &3._ki*t1*t4*t2**5-3._ki*t1*t5*t2**5-3._ki*t2**6*t1)/t2**6/t1/12._ki
              !
              !
              stemp1=stemp2+stemp3
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            case(4)
              !
              stemp2=(-4._ki*t1*t6**2-8._ki*t6**2*t2*t4+8._ki*t6**2*t4**2+2._ki*t6*&
                &*2*t2**2-2._ki*t4*t6*t2**2+2._ki*t6*t2**2*t5-2._ki*t6*t3*t2-4._ki*t&
                &6*t4*t5*t2+t2**3*t6+2._ki*t5*t2**3+2._ki*t5**2*t2**2)/t2**4*z_log&
                &(t1*t6/t2**2,1._ki)/6._ki
              !
              stemp3=-(2._ki*t2**3*t3**2*t5+6._ki*t1**2*t6**2*t5*t2-6._ki*t1*t6*t2&
                &**2*t3*t5+12._ki*t1*t6**2*t3*t4*t5+2._ki*t2**2*t3**2*t6*t5-4._ki*t&
                &2*t3**2*t6*t4*t5-6._ki*t1*t6**2*t2*t3*t5+2._ki*t2**2*t3**2*t6**2-&
                &8._ki*t2*t3**2*t6**2*t4+8._ki*t3**2*t6**2*t4**2+2._ki*t1*t6**2*t3*&
                &*2-2._ki*t6*t2*t3**3+2._ki*t2**2*t3**2*t5**2+6._ki*t1**2*t6**2*t5*&
                &*2-6._ki*t1*t6*t2*t3*t5**2+t2**3*t3**2*t6-2._ki*t2**2*t3**2*t6*t4&
                &-3._ki*t1*t6**2*t2**2*t3+6._ki*t1*t6**2*t2*t3*t4)/t3**2/t2**4*q(3&
                &,(t2*t3-t1*t6)/t2/t3,sign_arg)/6._ki-(4._ki*t1*t6**2*t3-24._ki*t6*&
                &*2*t3*t4**2+24._ki*t6**2*t2*t3*t4-6._ki*t6**2*t2**2*t3+6._ki*t6*t1&
                &*t2**2*t5+6._ki*t6*t1*t2*t5**2+6._ki*t2*t3**2*t6+20._ki*t6*t2*t3*t&
                &4*t5-10._ki*t6*t5*t2**2*t3-5._ki*t6*t3*t2**3+10._ki*t6*t2**2*t3*t4&
                &-12._ki*t2**2*t3*t5**2-12._ki*t2**3*t3*t5)/t2**4/t3/12._ki
              !
              stemp1=stemp2+stemp3
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            end select
            !
          case(2)
            !
            select case(par4_glob)
            !
            case(2)
              !
              stemp4=2._ki*t1*t6*t2**2*t3*t5+4._ki*t6**2*t1*t2**2*t4*t5+6._ki*t1*t&
                &6*t2*t3*t5**2-t4**2*t2**6/3._ki+t1**2*t6**2*t5*t2+16._ki*t2*t3**2&
                &*t6**2*t4+20._ki*t6**2*t1*t5**2*t4**2-4._ki*t6**2*t4**3*t5*t2**2-&
                &4._ki/3._ki*t6**2*t1*t4*t2**3+4._ki*t6**2*t4**4*t5*t2-48._ki*t6**2*&
                &t4**3*t3*t5+8._ki*t6**2*t4**3*t5**2*t2-8._ki*t1*t6**2*t2*t3*t4-t2&
                &**2*t1**2*t6**2/3._ki+4._ki*t1*t6**2*t3**2-8._ki*t1*t6**2*t2*t3*t5&
                &+2._ki*t6**2*t2**3*t3*t4+12._ki*t6**2*t4**3*t2*t3-10._ki*t6**2*t2*&
                &*2*t3*t4**2+4._ki/3._ki*t6*t1*t4*t2**4+6._ki*t2**2*t3**2*t6*t5-t3*&
                &t2**3*t1*t6-2._ki*t6**2*t4**2*t5**2*t2**2+t6**2*t4**2*t5*t2**3+1&
                &0._ki/3._ki*t6**2*t1*t4**2*t2**2+5._ki*t6*t4**2*t3*t2**3-4._ki*t6*t&
                &4**3*t5**3*t2-8._ki*t6**2*t4**4*t5**2-2._ki*t6*t2*t3**3
              !
              stemp3=-2._ki*t2**2*t3**2*t6**2-4._ki/3._ki*t6**2*t4**4*t2**2-4._ki/3&
                &._ki*t6*t4**3*t2**4+4._ki/3._ki*t6**2*t4**3*t2**3+2._ki/3._ki*t6*t4*&
                &*2*t2**5-2._ki*t1**2*t6**2*t5**2-26._ki*t3**2*t6**2*t4**2+32._ki*t&
                &1*t6**2*t3*t4*t5+4._ki*t6*t2**3*t3*t4*t5-10._ki*t6*t2**2*t3*t4**2&
                &*t5+2._ki*t6*t4**2*t5**3*t2**2-2._ki*t6*t4**3*t5**2*t2**2-2._ki*t6&
                &*t2**4*t3*t4+t6*t4**2*t5**2*t2**3-t6*t4**2*t5*t2**4+2._ki*t6*t4*&
                &*3*t5*t2**3+2._ki*t1*t6**2*t2**2*t3-4._ki*t2**2*t3**2*t6*t4+stemp&
                &4+4._ki*t6*t1*t4*t5**3*t2-2._ki*t6*t1*t4*t5*t2**3+2._ki*t6*t1*t4*t&
                &5**2*t2**2-24._ki*t2*t3**2*t6*t4*t5-10._ki*t6**2*t1*t4**2*t5*t2-8&
                &._ki*t6**2*t1*t4*t5**2*t2-30._ki*t6*t4**2*t3*t5**2*t2+12._ki*t6*t4&
                &*t3*t5**2*t2**2-t6**2*t4**2*t2**4/3._ki+t2**3*t3**2*t6+40._ki*t6*&
                &*2*t4**2*t3*t5*t2-8._ki*t6**2*t4*t3*t5*t2**2
              !
              stemp4=1._ki/t2**8*z_log(t1*t6/t2**2,1._ki)
              !
              stemp2=stemp3*stemp4
              !
              stemp4=t4**2/t2**2*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/3._ki
              !
              stemp8=-3._ki*t6*t4**2*t5**2*t1**3*t2**2+12._ki*t6*t4**3*t5**2*t1**&
                &3*t2-4._ki/3._ki*t6*t1**2*t4**2*t5*t3*t2**3-16._ki*t6*t1**2*t4**4*&
                &t5*t3*t2+8._ki*t6*t1**2*t4**3*t5*t3*t2**2-2._ki*t6*t4**3*t3*t2**3&
                &*t1**2+t6*t4**2*t2**4*t3*t1**2/3._ki-8._ki/3._ki*t6*t1**2*t4**5*t3&
                &*t2-4._ki/3._ki*t4*t3**2*t2**3*t1**2*t6+t4**2*t3**2*t2**4*t1*t6/6&
                &._ki+3._ki*t4**2*t5**2*t3*t2**3*t1**2+12._ki*t4*t3*t5**2*t1**3*t2*&
                &*2+3._ki/2._ki*t6*t4**2*t5*t2**3*t1**3+6._ki*t4**3*t3**2*t2**2*t1*&
                &*2-5._ki*t4**2*t3**2*t2**3*t1**2-2._ki*t4*t3*t2**4*t1**3+4._ki*t4*&
                &*4*t3**3*t2**3+t2**3*t3**2*t1**3-4._ki/3._ki*t4**3*t2**4*t1**3-2.&
                &_ki*t2*t3**3*t1**3-16._ki/3._ki*t4**5*t3**3*t2**2+t4**2*t3**3*t2**&
                &5/6._ki-16._ki/3._ki*t6*t1*t4**5*t3**2*t2-15._ki*t6*t4**2*t3*t2**2*&
                &t1**3+t6*t2**2*t3*t1**4
              !
              stemp7=-3._ki*t6*t3**2*t1**3*t2**2-12._ki*t6*t4**4*t5**2*t1**3+10._k&
                &i*t6*t4**2*t5**2*t1**4+5._ki/3._ki*t6*t4**2*t2**2*t1**4-2._ki*t6*t&
                &4**4*t2**2*t1**3-2._ki/3._ki*t6*t4*t2**3*t1**4+56._ki/3._ki*t6*t1**&
                &2*t4**4*t3**2+8._ki/3._ki*t6*t1*t4**6*t3**2-39._ki*t6*t4**2*t1**3*&
                &t3**2-t6*t1**3*t4**2*t2**4/2._ki+2._ki*t6*t1**3*t4**3*t2**3-12._ki&
                &*t1*t4**3*t3**3*t2**2-2._ki*t1*t4**4*t3**2*t2**3+28._ki/3._ki*t1*t&
                &4**4*t3**3*t2+4._ki/3._ki*t1*t4**5*t3**2*t2**2-12._ki*t1*t4**4*t3*&
                &*2*t2**2*t5+6._ki*t1*t4**3*t3**2*t2**3*t5-12._ki*t6*t1**3*t2**2*t&
                &3*t4*t5-24._ki*t4*t3**2*t1**3*t5*t2-30._ki*t4**2*t3**2*t5*t1**2*t&
                &2**2-30._ki*t4**2*t3*t5**2*t1**3*t2-6._ki*t6*t4**3*t5*t2**2*t1**3&
                &+18._ki*t6*t4**3*t3*t2*t1**3-4._ki*t6*t3*t5*t1**4*t2+60._ki*t6*t4*&
                &*2*t3*t5*t1**3*t2+stemp8
              !
              stemp8=4._ki*t6*t4**4*t3**2*t2**2*t1+3._ki*t6*t4*t3*t2**3*t1**3+16.&
                &_ki*t6*t3*t5*t1**4*t4-72._ki*t6*t4**3*t3*t5*t1**3-4._ki*t6*t2*t3*t&
                &1**4*t4-t1*t4**2*t3**2*t2**4*t5+8._ki*t1*t4**5*t3**2*t5*t2+4._ki*&
                &t6*t1**2*t4**4*t3*t2**2+36._ki*t4**3*t3**2*t5*t1**2*t2+10._ki*t6*&
                &t1**2*t2**2*t3**2*t4**2+32._ki/3._ki*t6*t1**2*t4**5*t5*t3-4._ki/3.&
                &_ki*t4**3*t3**2*t2**3*t1*t6+4._ki*t4**4*t5*t3*t2**2*t1**2+t4**2*t&
                &5*t3*t2**4*t1**2+6._ki*t1**2*t2**3*t3**2*t4*t5-2._ki*t4**4*t3*t2*&
                &*3*t1**2-t4**2*t3*t2**5*t1**2/2._ki-4._ki*t2**2*t3**2*t1**3*t4+t4&
                &*t3**2*t2**4*t1**2-2._ki*t4**3*t5**2*t2**2*t1**3+2._ki*t4**2*t5**&
                &3*t1**3*t2**2-2._ki/3._ki*t1*t4*t3**3*t2**4-t4**2*t3**2*t2**5*t1/&
                &6._ki+2._ki*t4**3*t5*t2**3*t1**3-t4**2*t5*t2**4*t1**3
              !
              stemp6=-4._ki*t4**3*t5**3*t1**3*t2+5._ki*t4**2*t3**3*t1*t2**3+5._ki*&
                &t4**2*t2**3*t1**3*t3+t4**3*t3**2*t2**4*t1+6._ki*t3**2*t5*t2**2*t&
                &1**3-8._ki*t3**3*t2**2*t1**2*t4+13._ki*t3**3*t2*t1**2*t4**2+24._ki&
                &*t6*t4*t3**2*t1**3*t2+6._ki*t6*t4**4*t5*t2*t1**3+stemp7+stemp8-2&
                &4._ki*t6*t1**2*t4**3*t3**2*t2+2._ki/3._ki*t4**2*t2**5*t1**3+t1**2*&
                &t2**3*t3**3+2._ki*t6*t3**2*t1**4-4._ki/3._ki*t4**3*t3**3*t2**4+8._k&
                &i/3._ki*t4**6*t3**3*t2-10._ki*t4**2*t5*t3*t2**2*t1**3-4._ki*t4**3*&
                &t5*t3*t2**3*t1**2+t4**2*t5**2*t2**3*t1**3-4._ki*t6*t4*t5**2*t1**&
                &4*t2+2._ki*t4**3*t3*t2**4*t1**2+4._ki*t4*t5*t3*t2**3*t1**3-12._ki*&
                &t4**3*t5**2*t3*t2**2*t1**2+12._ki*t4**4*t5**2*t3*t2*t1**2-5._ki*t&
                &6*t4**2*t5*t2*t1**4+2._ki*t6*t4*t5*t2**2*t1**4
              !
              stemp7=t6/t1**3/t2**8
              !
              stemp5=stemp6*stemp7
              !
              stemp3=stemp4+stemp5
              !
              stemp1=stemp2+stemp3
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            case(3)
              !
              stemp4=-72._ki*t6**2*t4**2*t3*t5*t2-t6*t1*t2**4*t5+2._ki*t6*t2**4*t&
                &3*t5+32._ki*t6**2*t4*t3*t5*t2**2+2._ki*t1**2*t6**2*t5**2+2._ki/3._k&
                &i*t6*t1*t2**5-2._ki/3._ki*t6**2*t1*t2**4+2._ki/3._ki*t6*t2**6*t4-t2&
                &**7*t4/3._ki+2._ki*t6**2*t1*t5*t2**3+t6**2*t4*t2**4*t5-2._ki*t6**2&
                &*t4*t5**2*t2**3-26._ki*t2*t3**2*t6**2*t4-t2**5*t6*t4*t5-20._ki*t6&
                &**2*t1*t5**2*t4**2+8._ki*t6**2*t4**3*t5*t2**2+10._ki/3._ki*t6**2*t&
                &1*t4*t2**3-4._ki*t6**2*t4**4*t5*t2+48._ki*t6**2*t4**3*t3*t5-16._ki&
                &*t6**2*t4**3*t5**2*t2-4._ki*t6**2*t2**3*t3*t5+2._ki*t6*t4*t5**3*t&
                &2**3+t6*t4*t5**2*t2**4-4._ki*t6**2*t1*t5**2*t2**2+6._ki*t6*t5**2*&
                &t2**3*t3+t6*t1*t5**2*t2**3+2._ki*t6*t1*t5**3*t2**2+30._ki*t6*t4**&
                &2*t3*t5**2*t2-30._ki*t6*t4*t3*t5**2*t2**2-4._ki*t6*t1*t4*t5**3*t2&
                &+2._ki*t6*t1*t4*t5*t2**3-2._ki*t6*t1*t4*t5**2*t2**2-6._ki*t1*t6*t2&
                &*t3*t5**2+2._ki*t6*t2*t3**3+6._ki*t2**2*t3**2*t6**2+5._ki/3._ki*t6*&
                &*2*t4**2*t2**4+4._ki/3._ki*t6**2*t4**4*t2**2-t4*t2**5*t6**2/3._ki+&
                &4._ki/3._ki*t6*t4**3*t2**4
              !
              stemp5=stemp4-8._ki/3._ki*t6**2*t4**3*t2**3-2._ki*t6*t4**2*t2**5-t6*&
                &t3*t2**5+t4**2*t2**6/3._ki+24._ki*t2*t3**2*t6*t4*t5-10._ki*t6*t2**&
                &3*t3*t4*t5+10._ki*t6*t2**2*t3*t4**2*t5+16._ki*t1*t6**2*t2*t3*t5-3&
                &2._ki*t1*t6**2*t3*t4*t5+26._ki*t3**2*t6**2*t4**2+t2**2*t1**2*t6**&
                &2/3._ki-2._ki*t2**3*t3**2*t6+t6**2*t2**4*t3-4._ki*t1*t6**2*t3**2+8&
                &._ki*t6**2*t4**4*t5**2+8._ki*t1*t6**2*t2*t3*t4-2._ki*t1*t6*t2**2*t&
                &3*t5-10._ki*t6**2*t1*t2**2*t4*t5+10._ki*t6**2*t1*t4**2*t5*t2
              !
              stemp3=stemp5+20._ki*t6**2*t1*t4*t5**2*t2-8._ki*t6**2*t2**3*t3*t4-1&
                &2._ki*t6**2*t4**3*t2*t3+18._ki*t6**2*t2**2*t3*t4**2-4._ki/3._ki*t6*&
                &t1*t4*t2**4+t3*t2**3*t1*t6+10._ki*t6**2*t4**2*t5**2*t2**2-5._ki*t&
                &6**2*t4**2*t5*t2**3-10._ki/3._ki*t6**2*t1*t4**2*t2**2-5._ki*t6*t4*&
                &*2*t3*t2**3+4._ki*t6*t4**3*t5**3*t2-6._ki*t6*t4**2*t5**3*t2**2+2.&
                &_ki*t6*t4**3*t5**2*t2**2+5._ki*t6*t2**4*t3*t4-3._ki*t6*t4**2*t5**2&
                &*t2**3+3._ki*t6*t4**2*t5*t2**4-2._ki*t6*t4**3*t5*t2**3-4._ki*t1*t6&
                &**2*t2**2*t3+4._ki*t2**2*t3**2*t6*t4-t1**2*t6**2*t5*t2-12._ki*t2*&
                &*2*t3**2*t6*t5
              !
              stemp4=1._ki/t2**8*z_log(t1*t6/t2**2,1._ki)
              !
              stemp2=stemp3*stemp4
              !
              stemp4=-(-t2+t4)*t4/t2**2*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/3._ki
              !
              stemp8=12._ki*t6*t4**3*t5*t2**2*t1**3-18._ki*t6*t4**3*t3*t2*t1**3-1&
                &5._ki/2._ki*t6*t4**2*t5*t2**3*t1**3+10._ki*t6*t4*t5**2*t1**4*t2+5.&
                &_ki*t6*t4**2*t5*t2*t1**4-8._ki*t1*t4**5*t3**2*t5*t2+20._ki*t1*t4**&
                &4*t3**2*t2**2*t5-18._ki*t1*t4**3*t3**2*t2**3*t5+30._ki*t4**2*t3*t&
                &5**2*t1**3*t2+10._ki*t4**2*t5*t3*t2**2*t1**3+8._ki*t4**3*t5*t3*t2&
                &**3*t1**2+2._ki/3._ki*t1**3*t2**6*t4-8._ki/3._ki*t4**6*t3**3*t2-t1*&
                &*3*t2**5*t3-2._ki*t2**3*t3**2*t1**3-3._ki*t1**2*t2**3*t3**3+2._ki*&
                &t2*t3**3*t1**3-6._ki*t6*t1**3*t2**3*t3*t5-4._ki/3._ki*t6*t1**2*t2*&
                &*4*t3*t4*t5+t6*t1**2*t2**5*t3*t4/3._ki+24._ki*t4*t3**2*t1**3*t5*t&
                &2-36._ki*t4**3*t3**2*t5*t1**2*t2-6._ki*t4**3*t3**2*t2**2*t1**2+9.&
                &_ki*t4**2*t3**2*t2**3*t1**2+5._ki*t4*t3*t2**4*t1**3-5._ki*t6*t4*t5&
                &*t2**2*t1**4+15._ki*t6*t4**2*t5**2*t1**3*t2**2-24._ki*t6*t4**3*t5&
                &**2*t1**3*t2+7._ki*t1*t4**2*t3**2*t2**4*t5-4._ki*t4**4*t5*t3*t2**&
                &2*t1**2-5._ki*t4**2*t5*t3*t2**4*t1**2-24._ki*t1**2*t2**3*t3**2*t4&
                &*t5
              !
              stemp7=-28._ki/3._ki*t4**4*t3**3*t2**3-2._ki*t6*t3**2*t1**4-3._ki/2._k&
                &i*t4**2*t3**3*t2**5-t6*t1**4*t2**4/3._ki+16._ki/3._ki*t4**3*t3**3*&
                &t2**4-10._ki*t4*t5*t3*t2**3*t1**3+24._ki*t4**3*t5**2*t3*t2**2*t1*&
                &*2-12._ki*t4**4*t5**2*t3*t2*t1**2-15._ki*t4**2*t5**2*t3*t2**3*t1*&
                &*2-30._ki*t4*t3*t5**2*t1**3*t2**2+54._ki*t4**2*t3**2*t5*t1**2*t2*&
                &*2+27._ki*t6*t4**2*t3*t2**2*t1**3+8._ki*t6*t3*t5*t1**4*t2-108._ki*&
                &t6*t4**2*t3*t5*t1**3*t2-28._ki/3._ki*t6*t4**4*t3**2*t2**2*t1-12._k&
                &i*t6*t4*t3*t2**3*t1**3-16._ki*t6*t3*t5*t1**4*t4-3._ki*t4**2*t5**2&
                &*t2**3*t1**3+2._ki*t4**4*t3*t2**3*t1**2+5._ki/2._ki*t4**2*t3*t2**5&
                &*t1**2+4._ki*t2**2*t3**2*t1**3*t4-4._ki*t4*t3**2*t2**4*t1**2+2._ki&
                &*t4**3*t5**2*t2**2*t1**3-6._ki*t4**2*t5**3*t1**3*t2**2-4._ki*t4**&
                &3*t3*t2**4*t1**2+11._ki/3._ki*t1*t4*t3**3*t2**4+7._ki/6._ki*t4**2*t&
                &3**2*t2**5*t1-2._ki*t4**3*t5*t2**3*t1**3+3._ki*t4**2*t5*t2**4*t1*&
                &*3+4._ki*t4**3*t5**3*t1**3*t2-13._ki*t4**2*t3**3*t1*t2**3-5._ki*t4&
                &**2*t2**3*t1**3*t3+stemp8
              !
              stemp8=stemp7-3._ki*t4**3*t3**2*t2**4*t1-12._ki*t3**2*t5*t2**2*t1**&
                &3+13._ki*t3**3*t2**2*t1**2*t4-13._ki*t3**3*t2*t1**2*t4**2+6._ki*t1&
                &**3*t2**3*t3*t5**2+3._ki*t1**2*t2**4*t3**2*t5+t4*t5**2*t2**4*t1*&
                &*3+2._ki*t4*t5**3*t2**3*t1**3-t1*t2**6*t3**2*t4/6._ki+56._ki/3._ki*&
                &t1*t4**3*t3**3*t2**2+10._ki/3._ki*t1*t4**4*t3**2*t2**3-28._ki/3._ki&
                &*t1*t4**4*t3**3*t2-4._ki/3._ki*t1*t4**5*t3**2*t2**2+2._ki*t2**4*t3&
                &*t5*t1**3-t4*t2**5*t1**3*t5-t4*t2**6*t1**2*t3/2._ki-t6*t1**3*t4*&
                &t2**5/2._ki+3._ki/2._ki*t6*t1**3*t2**4*t3-2._ki*t6*t2**2*t3*t1**4-2&
                &._ki*t4**2*t2**5*t1**3-t1*t4*t3**2*t2**5*t5+t2**5*t3*t4*t1**2*t5&
                &+3._ki*t4*t5**2*t2**4*t1**2*t3+16._ki/3._ki*t4**3*t3**2*t2**3*t1*t&
                &6+22._ki/3._ki*t4*t3**2*t2**3*t1**2*t6-3._ki/2._ki*t4**2*t3**2*t2**&
                &4*t1*t6+t4*t3**2*t2**5*t1*t6/6._ki+28._ki/3._ki*t6*t1**2*t4**2*t5*&
                &t3*t2**3+80._ki/3._ki*t6*t1**2*t4**4*t5*t3*t2-24._ki*t6*t1**2*t4**&
                &3*t5*t3*t2**2+6._ki*t6*t4**3*t3*t2**3*t1**2
              !
              stemp6=stemp8-7._ki/3._ki*t6*t4**2*t2**4*t3*t1**2+8._ki/3._ki*t6*t1**&
                &2*t4**5*t3*t2-t1*t3**3*t2**5/3._ki+8._ki*t6*t1*t4**5*t3**2*t2+48.&
                &_ki*t6*t1**3*t2**2*t3*t4*t5+t4*t3**3*t2**6/6._ki+4._ki/3._ki*t4**3*&
                &t2**4*t1**3+t1**2*t2**5*t3**2/2._ki+8._ki*t4**5*t3**3*t2**2+9._ki*&
                &t6*t3**2*t1**3*t2**2+12._ki*t6*t4**4*t5**2*t1**3-10._ki*t6*t4**2*&
                &t5**2*t1**4-5._ki/3._ki*t6*t4**2*t2**2*t1**4+2._ki*t6*t4**4*t2**2*&
                &t1**3+5._ki/3._ki*t6*t4*t2**3*t1**4-56._ki/3._ki*t6*t1**2*t4**4*t3*&
                &*2-8._ki/3._ki*t6*t1*t4**6*t3**2+39._ki*t6*t4**2*t1**3*t3**2+5._ki/&
                &2._ki*t6*t1**3*t4**2*t2**4+t6*t5*t2**3*t1**4-2._ki*t6*t1**4*t5**2&
                &*t2**2-4._ki*t6*t1**3*t4**3*t2**3-2._ki/3._ki*t3**2*t2**4*t1**2*t6&
                &+72._ki*t6*t4**3*t3*t5*t1**3+4._ki*t6*t2*t3*t1**4*t4-39._ki*t6*t4*&
                &t3**2*t1**3*t2-6._ki*t6*t4**4*t5*t2*t1**3-20._ki/3._ki*t6*t1**2*t4&
                &**4*t3*t2**2-26._ki*t6*t1**2*t2**2*t3**2*t4**2-32._ki/3._ki*t6*t1*&
                &*2*t4**5*t5*t3+112._ki/3._ki*t6*t1**2*t4**3*t3**2*t2-3._ki*t6*t4*t&
                &5**2*t2**3*t1**3+3._ki/2._ki*t6*t4*t5*t2**4*t1**3
              !
              stemp7=t6/t1**3/t2**8
              !
              stemp5=stemp6*stemp7
              !
              stemp3=stemp4+stemp5
              !
              stemp1=stemp2+stemp3
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            case(4)
              !
              stemp1=(-6._ki*t1*t6**2*t3+6._ki*t6**2*t1*t5*t2-18._ki*t6**2*t1*t5*t&
                &4-18._ki*t6**2*t2*t3*t4+3._ki*t6**2*t2**2*t4*t5+12._ki*t6**2*t5*t4&
                &**3-12._ki*t6**2*t4**2*t5*t2+3._ki*t6**2*t2**2*t3+24._ki*t6**2*t3*&
                &t4**2+2._ki*t1*t6*t2**3+2._ki*t6*t2**4*t4-4._ki*t6*t2**3*t4**2+t2*&
                &*4*t3+t4*t2**4*t5)/t2**6*z_log(t1*t6/t2**2,1._ki)/6._ki+(3._ki*t4*&
                &t5*t1*t6-t2*t3*t4*t5+t1*t6*t3-t2*t3**2-2._ki*t6*t4*t2*t3+4._ki*t4&
                &**2*t6*t3)/t2**3/t3*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/6._ki+(-6.&
                &_ki*t1**2*t6**2*t3+6._ki*t1**2*t6**2*t5*t2-18._ki*t1**2*t6**2*t4*t&
                &5-54._ki*t1*t6**2*t2*t3*t4+36._ki*t6**2*t1*t5*t4**3+72._ki*t1*t6**&
                &2*t4**2*t3-36._ki*t6**2*t1*t4**2*t5*t2+9._ki*t6**2*t1*t2**2*t4*t5&
                &+9._ki*t1*t6**2*t2**2*t3-16._ki*t6**2*t4**4*t3+24._ki*t6**2*t4**3*&
                &t2*t3-12._ki*t6**2*t2**2*t3*t4**2+2._ki*t6**2*t2**3*t3*t4+4._ki*t6&
                &*t1*t4*t2**4-8._ki*t6*t1*t2**3*t4**2+t1*t2**4*t3+3._ki*t1*t4*t2**&
                &4*t5)/t2**6/t1/12._ki
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            end select
            !
          case(3)
            !
            select case(par4_glob)
            !
            case(3)
              !
              stemp5=-t2**8/3._ki-t2**6*t6**2/3._ki-4._ki/3._ki*t6*t4**3*t2**4-32._k&
                &i*t6**2*t1*t4*t5**2*t2-30._ki*t6*t4**2*t3*t5**2*t2+48._ki*t6*t4*t&
                &3*t5**2*t2**2+4._ki*t6*t1*t4*t5**3*t2+16._ki*t6*t2**3*t3*t4*t5-10&
                &._ki*t6*t2**2*t3*t4**2*t5-t3*t2**3*t1*t6+2._ki/3._ki*t6*t2**7+3._ki&
                &*t2**3*t3**2*t6+18._ki*t6**2*t2**3*t3*t4+12._ki*t6**2*t4**3*t2*t3&
                &-26._ki*t6**2*t2**2*t3*t4**2+10._ki*t6*t4**2*t5**3*t2**2-2._ki*t6*&
                &t4**3*t5**2*t2**2-8._ki*t6*t2**4*t3*t4+5._ki*t6*t4**2*t5**2*t2**3&
                &-5._ki*t6*t4**2*t5*t2**4+2._ki*t6*t4**3*t5*t2**3
              !
              stemp4=stemp5+6._ki*t1*t6**2*t2**2*t3-4._ki*t2**2*t3**2*t6*t4+t1**2&
                &*t6**2*t5*t2-6._ki*t6**2*t1*t5*t2**3-6._ki*t6**2*t4*t2**4*t5+12._k&
                &i*t6**2*t4*t5**2*t2**3+36._ki*t2*t3**2*t6**2*t4+4._ki*t2**5*t6*t4&
                &*t5+20._ki*t6**2*t1*t5**2*t4**2-16._ki/3._ki*t6**2*t1*t4*t2**3+4._k&
                &i*t6**2*t4**4*t5*t2-48._ki*t6**2*t4**3*t3*t5+24._ki*t6**2*t4**3*t&
                &5**2*t2+16._ki*t6**2*t2**3*t3*t5-8._ki*t6*t4*t5**3*t2**3-4._ki*t6*&
                &t4*t5**2*t2**4+12._ki*t6**2*t1*t5**2*t2**2-18._ki*t6*t5**2*t2**3*&
                &t3-2._ki*t6*t1*t5**2*t2**3-4._ki*t6*t1*t5**3*t2**2+2._ki*t6*t1*t2*&
                &*4*t5-6._ki*t6*t2**4*t3*t5
              !
              stemp5=4._ki/3._ki*t6*t1*t4*t2**4+18._ki*t2**2*t3**2*t6*t5+104._ki*t6&
                &**2*t4**2*t3*t5*t2-72._ki*t6**2*t4*t3*t5*t2**2-12._ki*t6**2*t4**3&
                &*t5*t2**2-t4**2*t2**6/3._ki-2._ki*t6*t1*t4*t5*t2**3+2._ki*t6*t1*t4&
                &*t5**2*t2**2+6._ki*t1*t6*t2*t3*t5**2-24._ki*t2*t3**2*t6*t4*t5-24.&
                &_ki*t1*t6**2*t2*t3*t5+32._ki*t1*t6**2*t3*t4*t5+stemp4-26._ki*t6**2&
                &*t4**2*t5**2*t2**2+13._ki*t6**2*t4**2*t5*t2**3+10._ki/3._ki*t6**2*&
                &t1*t4**2*t2**2+5._ki*t6*t4**2*t3*t2**3-4._ki*t6*t4**3*t5**3*t2-8.&
                &_ki*t1*t6**2*t2*t3*t4+2._ki*t1*t6*t2**2*t3*t5+16._ki*t6**2*t1*t2**&
                &2*t4*t5-10._ki*t6**2*t1*t4**2*t5*t2
              !
              stemp3=stemp5+2._ki/3._ki*t2**7*t4-4._ki*t6**2*t2**4*t3+3._ki*t6*t3*t&
                &2**5+t6**2*t2**5*t5+t6*t2**5*t5**2-t2**2*t1**2*t6**2/3._ki+4._ki*&
                &t1*t6**2*t3**2-8._ki*t6**2*t4**4*t5**2-2._ki*t6*t2*t3**3-12._ki*t2&
                &**2*t3**2*t6**2-13._ki/3._ki*t6**2*t4**2*t2**4-4._ki/3._ki*t6**2*t4&
                &**4*t2**2+4._ki*t6**2*t4**3*t2**3+10._ki/3._ki*t6*t4**2*t2**5+2._ki&
                &*t4*t2**5*t6**2-2._ki*t1**2*t6**2*t5**2-26._ki*t3**2*t6**2*t4**2-&
                &4._ki/3._ki*t6*t1*t2**5+2._ki*t6**2*t1*t2**4-8._ki/3._ki*t6*t2**6*t4&
                &-t5*t2**6*t6-2._ki*t6**2*t2**4*t5**2+2._ki*t6*t2**4*t5**3
              !
              stemp4=1._ki/t2**8*z_log(t1*t6/t2**2,1._ki)
              !
              stemp2=stemp3*stemp4
              !
              stemp4=(-t2+t4)**2/t2**2*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/3._ki
              !
              stemp8=54._ki*t6*t4*t3**2*t1**3*t2+6._ki*t6*t4**4*t5*t2*t1**3+39._ki&
                &/2._ki*t6*t4**2*t5*t2**3*t1**3+t3**3*t2**7/6._ki-10._ki*t4**2*t5*t&
                &3*t2**2*t1**3-12._ki*t4**3*t5*t3*t2**3*t1**2+4._ki*t4**4*t5*t3*t2&
                &**2*t1**2+2._ki/3._ki*t1**3*t2**7-112._ki/3._ki*t6*t1**2*t4**4*t5*t&
                &3*t2+152._ki/3._ki*t6*t1**2*t4**3*t5*t3*t2**2-38._ki/3._ki*t6*t4**3&
                &*t3*t2**3*t1**2-2._ki*t1**2*t2**5*t3**2+3._ki*t2**3*t3**2*t1**3+5&
                &2._ki/3._ki*t4**4*t3**3*t2**3-6._ki*t2**5*t3*t4*t1**2*t5-18._ki*t4*&
                &t5**2*t2**4*t1**2*t3-14._ki/3._ki*t1*t4**4*t3**2*t2**3+28._ki/3._ki&
                &*t1*t4**4*t3**3*t2+4._ki/3._ki*t1*t4**5*t3**2*t2**2-6._ki*t2**4*t3&
                &*t5*t1**3+4._ki*t4*t2**5*t1**3*t5+3._ki*t4*t2**6*t1**2*t3+t3**2*t&
                &2**6*t1*t6/6._ki-6._ki*t6*t1**3*t2**4*t3-3._ki*t6*t5**2*t2**4*t1**&
                &3-76._ki/3._ki*t1*t4**3*t3**3*t2**2+t6*t1**4*t2**4-t1**2*t2**7*t3&
                &/2._ki+10._ki/3._ki*t4**2*t2**5*t1**3+5._ki/3._ki*t1*t3**3*t2**5-4._k&
                &i/3._ki*t4**3*t2**4*t1**3-44._ki/3._ki*t4**3*t3**3*t2**4+6._ki*t1**&
                &2*t2**3*t3**3-t2**6*t1**3*t5+2._ki*t5**3*t2**4*t1**3+41._ki/6._ki*&
                &t4**2*t3**3*t2**5
              !
              stemp9=25._ki/3._ki*t6*t4**2*t2**4*t3*t1**2-8._ki/3._ki*t6*t1**2*t4**&
                &5*t3*t2+28._ki/3._ki*t6*t1**2*t4**4*t3*t2**2+50._ki*t6*t1**2*t2**2&
                &*t3**2*t4**2+3._ki*t6*t1**3*t4*t2**5+41._ki/6._ki*t4**2*t3**2*t2**&
                &4*t1*t6-5._ki/3._ki*t4*t3**2*t2**5*t1*t6-39._ki*t6*t4**2*t3*t2**2*&
                &t1**3-18._ki*t6*t4**3*t5*t2**2*t1**3+18._ki*t6*t4**3*t3*t2*t1**3+&
                &48._ki*t4*t3*t5**2*t1**3*t2**2-24._ki*t4*t3**2*t1**3*t5*t2+36._ki*&
                &t4**3*t3**2*t5*t1**2*t2-t1*t2**7*t3**2/6._ki+stemp8+13._ki*t4**2*&
                &t5*t3*t2**4*t1**2+54._ki*t1**2*t2**3*t3**2*t4*t5+8._ki*t1*t4*t3**&
                &2*t2**5*t5
              !
              stemp7=stemp9-25._ki*t1*t4**2*t3**2*t2**4*t5+8._ki*t1*t4**5*t3**2*t&
                &5*t2-28._ki*t1*t4**4*t3**2*t2**2*t5+3._ki/2._ki*t6*t5*t2**5*t1**3+&
                &t6*t1**2*t2**6*t3/3._ki+3._ki*t6*t2**2*t3*t1**4-18._ki*t6*t3**2*t1&
                &**3*t2**2-12._ki*t6*t4**4*t5**2*t1**3+10._ki*t6*t4**2*t5**2*t1**4&
                &+5._ki/3._ki*t6*t4**2*t2**2*t1**4-2._ki*t6*t4**4*t2**2*t1**3-8._ki/&
                &3._ki*t6*t4*t2**3*t1**4+56._ki/3._ki*t6*t1**2*t4**4*t3**2+8._ki/3._k&
                &i*t6*t1*t4**6*t3**2-39._ki*t6*t4**2*t1**3*t3**2-13._ki/2._ki*t6*t1&
                &**3*t4**2*t2**4-3._ki*t6*t5*t2**3*t1**4+6._ki*t6*t1**4*t5**2*t2**&
                &2+6._ki*t6*t1**3*t4**3*t2**3
              !
              stemp9=10._ki/3._ki*t3**2*t2**4*t1**2*t6+32._ki/3._ki*t6*t1**2*t4**5*&
                &t5*t3-152._ki/3._ki*t6*t1**2*t4**3*t3**2*t2+18._ki*t6*t4*t5**2*t2*&
                &*3*t1**3-9._ki*t6*t4*t5*t2**4*t1**3+24._ki*t6*t1**3*t2**3*t3*t5+3&
                &2._ki/3._ki*t6*t1**2*t2**4*t3*t4*t5-8._ki/3._ki*t6*t1**2*t2**5*t3*t&
                &4-44._ki/3._ki*t4**3*t3**2*t2**3*t1*t6-64._ki/3._ki*t4*t3**2*t2**3*&
                &t1**2*t6+38._ki*t1*t4**3*t3**2*t2**3*t5-78._ki*t4**2*t3**2*t5*t1*&
                &*2*t2**2-30._ki*t4**2*t3*t5**2*t1**3*t2+16._ki*t4*t5*t3*t2**3*t1*&
                &*3-36._ki*t4**3*t5**2*t3*t2**2*t1**2+12._ki*t4**4*t5**2*t3*t2*t1*&
                &*2+39._ki*t4**2*t5**2*t3*t2**3*t1**2-4._ki/3._ki*t6*t1**2*t2**5*t3&
                &*t5
              !
              stemp8=stemp9-32._ki/3._ki*t6*t1*t4**5*t3**2*t2-108._ki*t6*t1**3*t2*&
                &*2*t3*t4*t5-12._ki*t6*t3*t5*t1**4*t2+156._ki*t6*t4**2*t3*t5*t1**3&
                &*t2+52._ki/3._ki*t6*t4**4*t3**2*t2**2*t1+27._ki*t6*t4*t3*t2**3*t1*&
                &*3+16._ki*t6*t3*t5*t1**4*t4-72._ki*t6*t4**3*t3*t5*t1**3-4._ki*t6*t&
                &2*t3*t1**4*t4+6._ki*t4**3*t3**2*t2**2*t1**2-13._ki*t4**2*t3**2*t2&
                &**3*t1**2-8._ki*t4*t3*t2**4*t1**3+5._ki*t4**2*t5**2*t2**3*t1**3-2&
                &._ki*t4**4*t3*t2**3*t1**2-13._ki/2._ki*t4**2*t3*t2**5*t1**2-4._ki*t&
                &2**2*t3**2*t1**3*t4+9._ki*t4*t3**2*t2**4*t1**2-2._ki*t4**3*t5**2*&
                &t2**2*t1**3
              !
              stemp6=stemp8+10._ki*t4**2*t5**3*t1**3*t2**2+6._ki*t4**3*t3*t2**4*t&
                &1**2-32._ki/3._ki*t1*t4*t3**3*t2**4-25._ki/6._ki*t4**2*t3**2*t2**5*&
                &t1+2._ki*t4**3*t5*t2**3*t1**3-5._ki*t4**2*t5*t2**4*t1**3-4._ki*t4*&
                &*3*t5**3*t1**3*t2+25._ki*t4**2*t3**3*t1*t2**3+5._ki*t4**2*t2**3*t&
                &1**3*t3+19._ki/3._ki*t4**3*t3**2*t2**4*t1+18._ki*t3**2*t5*t2**2*t1&
                &**3-18._ki*t3**3*t2**2*t1**2*t4+13._ki*t3**3*t2*t1**2*t4**2+t2**6&
                &*t3*t1**2*t5+3._ki*t2**5*t3*t5**2*t1**2-18._ki*t1**3*t2**3*t3*t5*&
                &*2-12._ki*t1**2*t2**4*t3**2*t5-4._ki*t4*t5**2*t2**4*t1**3-8._ki*t4&
                &*t5**3*t2**3*t1**3-t1*t2**6*t3**2*t5+4._ki/3._ki*t1*t2**6*t3**2*t&
                &4-32._ki/3._ki*t4**5*t3**3*t2**2-t6*t1**3*t2**6/2._ki+8._ki/3._ki*t4&
                &**6*t3**3*t2-2._ki*t2*t3**3*t1**3-5._ki/3._ki*t4*t3**3*t2**6+2._ki*&
                &t6*t3**2*t1**4+3._ki*t1**3*t2**5*t3-8._ki/3._ki*t1**3*t2**6*t4+t5*&
                &*2*t2**5*t1**3-16._ki*t6*t4*t5**2*t1**4*t2-5._ki*t6*t4**2*t5*t2*t&
                &1**4+8._ki*t6*t4*t5*t2**2*t1**4-39._ki*t6*t4**2*t5**2*t1**3*t2**2&
                &+36._ki*t6*t4**3*t5**2*t1**3*t2-100._ki/3._ki*t6*t1**2*t4**2*t5*t3&
                &*t2**3+stemp7
              !
              stemp7=t6/t1**3/t2**8
              !
              stemp5=stemp6*stemp7
              !
              stemp3=stemp4+stemp5
              !
              stemp1=stemp2+stemp3
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            case(4)
              !
              stemp2=-(12._ki*t6**2*t1*t5*t2-18._ki*t6**2*t1*t5*t4-6._ki*t1*t6**2*&
                &t3-30._ki*t6**2*t2*t3*t4+12._ki*t6**2*t5*t4**3+15._ki*t6**2*t2**2*&
                &t4*t5-24._ki*t6**2*t4**2*t5*t2+24._ki*t6**2*t3*t4**2+9._ki*t6**2*t&
                &2**2*t3-3._ki*t6**2*t2**3*t5+2._ki*t1*t6*t2**3-2._ki*t6*t2**5+6._ki&
                &*t6*t2**4*t4-4._ki*t6*t2**3*t4**2+t4*t2**4*t5-t5*t2**5+t2**4*t3)&
                &/t2**6*z_log(t1*t6/t2**2,1._ki)/6._ki
              !
              stemp3=-(-3._ki*t5*t2*t1*t6+2._ki*t6*t2**2*t3+t2**2*t3*t5-6._ki*t6*t&
                &4*t2*t3+4._ki*t4**2*t6*t3-t2*t3*t4*t5+3._ki*t4*t5*t1*t6+t1*t6*t3-&
                &t2*t3**2)/t2**3/t3*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/6._ki-(12._k&
                &i*t1**2*t6**2*t5*t2-18._ki*t1**2*t6**2*t4*t5-6._ki*t1**2*t6**2*t3&
                &-9._ki*t6**2*t1*t5*t2**3+72._ki*t1*t6**2*t4**2*t3-90._ki*t1*t6**2*&
                &t2*t3*t4+27._ki*t1*t6**2*t2**2*t3-72._ki*t6**2*t1*t4**2*t5*t2+45.&
                &_ki*t6**2*t1*t2**2*t4*t5+36._ki*t6**2*t1*t5*t4**3-36._ki*t6**2*t2*&
                &*2*t3*t4**2+40._ki*t6**2*t4**3*t2*t3+14._ki*t6**2*t2**3*t3*t4-2._k&
                &i*t6**2*t2**4*t3-16._ki*t6**2*t4**4*t3-8._ki*t6*t1*t2**3*t4**2+12&
                &._ki*t6*t1*t4*t2**4-4._ki*t6*t1*t2**5+3._ki*t1*t4*t2**4*t5+t1*t2**&
                &4*t3-3._ki*t1*t5*t2**5)/t2**6/t1/12._ki
              !
              stemp1=stemp2+stemp3
              !
              stemp2=1._ki/t2
              !
              temp0=stemp1*stemp2
              !
            end select
            !
          case(4)
            !
            select case(par4_glob)
            !
            case(4)
              !
              temp0=(-(-2._ki*t1*t6**2+4._ki*t6**2*t4**2-4._ki*t6**2*t2*t4+t6**2*t&
                &2**2-2._ki*t6*t4*t5*t2-t6*t3*t2+t6*t2**2*t5+t5**2*t2**2)/t2**4*z&
                &_log(t1*t6/t2**2,1._ki)/3._ki+(6._ki*t1*t6**2*t3*t4*t5-3._ki*t1*t6*&
                &t2*t3*t5**2+t1*t6**2*t3**2-3._ki*t1*t6**2*t2*t3*t5+3._ki*t1**2*t6&
                &**2*t5**2+t2**2*t3**2*t6**2+t2**2*t3**2*t5**2-t6*t2*t3**3-4._ki*&
                &t2*t3**2*t6**2*t4+t2**2*t3**2*t6*t5+4._ki*t3**2*t6**2*t4**2-2._ki&
                &*t2*t3**2*t6*t4*t5)/t3**2/t2**4*q(3,(t2*t3-t1*t6)/t2/t3,sign_ar&
                &g)/3._ki+(2._ki*t1*t6**2*t3-12._ki*t6**2*t3*t4**2+12._ki*t6**2*t2*t&
                &3*t4-3._ki*t6**2*t2**2*t3+3._ki*t6*t1*t2*t5**2-5._ki*t6*t5*t2**2*t&
                &3+10._ki*t6*t2*t3*t4*t5+3._ki*t2*t3**2*t6-6._ki*t2**2*t3*t5**2)/t2&
                &**4/t3/6)/t2
              !
            end select
            !
          end select
          !
        else if (nb_par == 3) then
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
                stemp2=-(4._ki*t6**2*t2*t3*t4-2._ki*t6*t4*t5**2*t2**2-4._ki*t6*t2**3&
                  &*t4*t5-2._ki*t6*t2**2*t3*t5-2._ki*t6**2*t1*t5*t2-2._ki*t6*t3*t2**3&
                  &+t6**2*t5*t2**3-2._ki*t6**2*t2**2*t3+3._ki*t5*t2**5+t2**5*t6+t6**&
                  &3*t2**3-8._ki*t6**3*t4**3-4._ki*t6**2*t5*t2**2*t4+4._ki*t6**2*t5*t&
                  &2*t4**2+12._ki*t6**3*t4**2*t2+t6*t2**3*t5**2-2._ki*t6*t2**4*t4+t5&
                  &**3*t2**3+t2**6+3._ki*t5**2*t2**4+t6**2*t2**4+2._ki*t6*t2**4*t5+1&
                  &2._ki*t6**3*t1*t4-6._ki*t6**3*t2**2*t4-4._ki*t6**2*t2**3*t4+4._ki*t&
                  &6**2*t2**2*t4**2-6._ki*t6**3*t1*t2-2._ki*t6**2*t1*t2**2)/t2**6*z_&
                  &log(t1*t6/t2**2,1._ki)/4._ki
                !
                stemp4=-(2._ki*t1**2*t6**2*t2**2-2._ki*t1*t6**2*t2**2*t3+t3**2*t6**&
                  &2*t2**2+4._ki*t1**2*t6**2*t5*t2+4._ki*t1*t6**2*t2*t3*t4-2._ki*t3*t&
                  &1*t6**2*t5*t2-4._ki*t3**2*t6**2*t4*t2+2._ki*t6**2*t1**2*t5**2+4._k&
                  &i*t3*t6**2*t4*t1*t5+2._ki*t3**2*t1*t6**2+4._ki*t3**2*t6**2*t4**2-&
                  &2._ki*t6*t1*t3*t2**3-4._ki*t1*t6*t2**2*t3*t5-2._ki*t6*t1*t3*t5**2*&
                  &t2-2._ki*t2*t3**3*t6+t3**2*t2**4+2._ki*t2**3*t3**2*t5+t3**2*t5**2&
                  &*t2**2)*(2._ki*t6*t3*t4+2._ki*t1*t6*t5-t3*t6*t2+2._ki*t1*t6*t2-t3*&
                  &t2*t5-t2**2*t3)/t3**3/t2**6*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/4&
                  &._ki
                !
                stemp7=2._ki*t6*t1*t3*t2**3*t5**2-11._ki/12._ki*t3**2*t2**3*t5**3+2.&
                  &_ki*t6*t1*t5*t2**4*t3+2._ki/3._ki*t6*t1*t3*t2**2*t5**3+7._ki/6._ki*t&
                  &6*t3**2*t2**3*t4*t5-11._ki/24._ki*t6**3*t2**3*t3**2-17._ki/6._ki*t6&
                  &**2*t2*t3**2*t4**2*t5+17._ki/6._ki*t6**2*t3**2*t4*t5*t2**2+11._ki/&
                  &3._ki*t6**3*t3**2*t4**3+t3**3*t2**3*t6/12._ki+5._ki/3._ki*t6*t3**2*&
                  &t2**2*t4*t5**2-11._ki/12._ki*t3**2*t2**6+t6**2*t1*t2**4*t3/2._ki+1&
                  &7._ki/6._ki*t6**2*t2**3*t3**2*t4-11._ki/6._ki*t6**2*t2*t3**3*t4-17.&
                  &_ki/24._ki*t6**2*t2**3*t3**2*t5-17._ki/6._ki*t6**2*t3**2*t4**2*t2**&
                  &2+2._ki/3._ki*t6*t1*t3*t2**5+5._ki/4._ki*t6**3*t1*t3**2*t2-5._ki/2._k&
                  &i*t6**3*t3**2*t1*t4-11._ki/2._ki*t6**3*t3**2*t2*t4**2
                !
                stemp6=stemp7+11._ki/4._ki*t6**3*t2**2*t3**2*t4-t6**2*t1**2*t5**2*t&
                  &2**2-t3**2*t2**2*t1*t6**2/12._ki+7._ki/6._ki*t6*t3**3*t2**2*t5-5._k&
                  &i/6._ki*t6*t3**2*t2**3*t5**2+5._ki/3._ki*t6*t2**4*t3**2*t4-7._ki/12&
                  &._ki*t6*t2**4*t3**2*t5-t6**2*t1**2*t5*t2**3-t6**2*t1**2*t2*t5**3&
                  &/3._ki-67._ki/12._ki*t3**2*t2**4*t5**2-17._ki/24._ki*t6**2*t3**2*t2*&
                  &*4-5._ki/6._ki*t6*t2**5*t3**2+11._ki/12._ki*t6**2*t2**2*t3**3+t6**2&
                  &*t1*t2**3*t3*t5-59._ki/12._ki*t2**5*t3**2*t5-t1**2*t2**4*t6**2/3.&
                  &_ki-t6**2*t1*t2*t3*t4*t5**2-t6**2*t1*t2*t5*t3**2/12._ki-2._ki*t6**&
                  &2*t1*t3*t4*t5*t2**2-t6**2*t1*t2**3*t3*t4+t6**2*t1*t5**2*t3*t2**&
                  &2/2._ki
                !
                stemp7=1._ki/t3**2/t2**6
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(2)
                !
                stemp2=-(2._ki*t3*t5*t2**5+t4*t2**7-36._ki*t6**3*t1**2*t5+2._ki*t6*t&
                  &1*t2**5-216._ki*t1*t6**3*t4*t5*t2-6._ki*t3*t4*t6*t2**4-96._ki*t6**&
                  &3*t4**4*t5-12._ki*t6**2*t2**4*t4**2-240._ki*t6**3*t3*t4**3+12._ki*&
                  &t6**3*t3*t2**3+2._ki*t6*t2**5*t3+12._ki*t6**2*t2**3*t4**3-18._ki*t&
                  &6**2*t1*t2**3*t4+288._ki*t6**3*t1*t5*t4**2+36._ki*t1*t6**3*t2**2*&
                  &t5+12._ki*t6**3*t2**3*t5*t4-72._ki*t6**3*t2**2*t5*t4**2+2._ki*t6*t&
                  &2**5*t4*t5-4._ki*t6*t2**4*t4**2*t5+180._ki*t6**3*t4*t1*t3+288._ki*&
                  &t6**3*t3*t2*t4**2-108._ki*t6**3*t3*t2**2*t4+2._ki*t6*t1*t2**4*t5-&
                  &72._ki*t6**3*t3*t2*t1+144._ki*t6**3*t4**3*t5*t2+2._ki*t6*t2**6*t4-&
                  &4._ki*t6*t2**5*t4**2+3._ki*t6**2*t2**5*t4+2._ki*t2**6*t4*t5+2._ki*t&
                  &3*t2**6+t4*t5**2*t2**5+6._ki*t6**2*t1*t2**4)/t2**8*z_log(t1*t6/t&
                  &2**2,1._ki)/12._ki
                !
                stemp4=(-4._ki*t3*t4*t5**2*t2*t1*t6+2._ki*t3**2*t2**3*t4*t5+4._ki*t1&
                  &**2*t6**2*t3*t2-6._ki*t1*t6*t3**2*t2**2-6._ki*t2*t6*t3**3*t4-2._ki&
                  &*t2*t1*t6**2*t3**2+6._ki*t4*t2**2*t1**2*t6**2+3._ki*t2**2*t4*t3**&
                  &2*t6**2+2._ki*t2**3*t4*t6*t3**2-4._ki*t2**2*t4**2*t6*t3**2+t3**2*&
                  &t4*t5**2*t2**2+12._ki*t4**3*t3**2*t6**2+2._ki*t4*t5*t3**2*t6*t2**&
                  &2-8._ki*t2**2*t4*t1*t6**2*t3+16._ki*t2*t4**2*t1*t6**2*t3+2._ki*t2*&
                  &*3*t3**3+t4*t3**2*t2**4+2._ki*t2**2*t5*t3**3+2._ki*t6*t3**3*t2**2&
                  &+4._ki*t5*t1**2*t6**2*t3+6._ki*t1*t6**2*t3**2*t4-12._ki*t2*t4**2*t&
                  &3**2*t6**2+6._ki*t4*t5**2*t6**2*t1**2+12._ki*t4*t5*t2*t6**2*t1**2&
                  &-4._ki*t4*t3*t2**3*t1*t6-6._ki*t2*t5*t1*t6*t3**2-8._ki*t3*t4*t5*t2&
                  &*t6**2*t1+16._ki*t3*t4**2*t5*t6**2*t1-4._ki*t3**2*t4**2*t5*t2*t6-&
                  &8._ki*t4*t3*t5*t2**2*t1*t6)/t2**5/t3**2*q(4,(t2*t3-t1*t6)/t2/t3,&
                  &sign_arg)/12._ki
                !
                stemp7=-5._ki/36._ki*t1*t2**6*t3*t4*t5-7._ki/36._ki*t1*t2**5*t3*t4*t5&
                  &**2-7._ki/18._ki*t6*t1*t2**6*t4*t3+t6*t1**2*t2**4*t4*t5**2/6._ki+2&
                  &._ki/3._ki*t6*t1*t4*t3**2*t2**4-2._ki/9._ki*t6*t1*t2**5*t3**2+t6*t1&
                  &**2*t2**5*t3/9._ki+t6**3*t3*t5*t1**3-25._ki/2._ki*t6**3*t1**2*t4*t&
                  &3**2+5._ki*t6**3*t1**2*t2*t3**2+110._ki/3._ki*t6**3*t1*t4**3*t3**2&
                  &+t1*t3**2*t2**6/36._ki+7._ki/9._ki*t6*t1*t2**5*t4**2*t3+15._ki*t6**&
                  &3*t1**2*t3*t4*t5*t2+3._ki/2._ki*t6**2*t1*t3*t2**4*t4**2-3._ki/2._ki&
                  &*t6**2*t1*t3*t2**3*t4**3-11._ki/6._ki*t6**3*t1*t3*t2**3*t5*t4+3._k&
                  &i/4._ki*t6**2*t1**2*t3*t2**3*t4-3._ki/8._ki*t6**2*t1*t3*t2**5*t4+t&
                  &6*t1**2*t2**6*t4/6._ki
                !
                stemp6=stemp7-t6**3*t4*t3**2*t2**4/4._ki-11._ki/6._ki*t6**3*t1*t3**2&
                  &*t2**3-6._ki*t6**3*t3**2*t2**2*t4**3-t6**2*t1**2*t3*t2**4/4._ki-4&
                  &4._ki*t6**3*t1*t2*t4**2*t3**2+11._ki*t6**3*t1*t3*t2**2*t5*t4**2+2&
                  &._ki*t6**3*t3**2*t2**3*t4**2+8._ki*t6**3*t3**2*t4**4*t2-2._ki/9._ki&
                  &*t1*t2**5*t3**2*t5-7._ki/36._ki*t1*t2**7*t3*t4+t6*t1**2*t2**5*t4*&
                  &t5/3._ki+44._ki/3._ki*t6**3*t1*t3*t4**4*t5+33._ki/2._ki*t6**3*t1*t2*&
                  &*2*t4*t3**2-20._ki*t6**3*t1**2*t3*t5*t4**2-5._ki/2._ki*t6**3*t1**2&
                  &*t3*t2**2*t5+t6*t1**2*t2**4*t5*t3/9._ki-22._ki*t6**3*t1*t3*t4**3*&
                  &t5*t2+7._ki/9._ki*t6*t1*t2**4*t3*t4**2*t5-7._ki/18._ki*t6*t1*t2**5*&
                  &t3*t4*t5-4._ki*t6**3*t3**2*t4**5
                !
                stemp7=1._ki/t1/t2**8/t3
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(3)
                !
                stemp2=-(-2._ki*t2**6*t4*t5+360._ki*t1*t6**3*t4*t5*t2+6._ki*t3*t4*t6&
                  &*t2**4+18._ki*t6**2*t1*t2**3*t4-6._ki*t6*t2**5*t4*t5+4._ki*t6*t2**&
                  &4*t4**2*t5+96._ki*t6**3*t4**4*t5+24._ki*t6**2*t2**4*t4**2-12._ki*t&
                  &6**2*t2**3*t4**3-6._ki*t6*t2**6*t4+4._ki*t6*t2**5*t4**2-t4*t2**7-&
                  &15._ki*t6**2*t2**5*t4+240._ki*t6**3*t3*t4**3-48._ki*t6**3*t3*t2**3&
                  &-2._ki*t3*t5*t2**5+2._ki*t5*t2**7-2._ki*t3*t2**6+t5**2*t2**6+2._ki*&
                  &t6*t2**7-288._ki*t6**3*t1*t5*t4**2-108._ki*t1*t6**3*t2**2*t5-84._k&
                  &i*t6**3*t2**3*t5*t4+216._ki*t6**3*t2**2*t5*t4**2-240._ki*t6**3*t4&
                  &**3*t5*t2+108._ki*t6**3*t3*t2*t1-180._ki*t6**3*t4*t1*t3-432._ki*t6&
                  &**3*t3*t2*t4**2+252._ki*t6**3*t3*t2**2*t4-2._ki*t6*t1*t2**4*t5-t4&
                  &*t5**2*t2**5-12._ki*t6**2*t1*t2**4-4._ki*t6*t2**5*t3+36._ki*t6**3*&
                  &t1**2*t5+12._ki*t6**3*t5*t2**4+2._ki*t6*t2**6*t5-2._ki*t6*t1*t2**5&
                  &+t2**8+3._ki*t6**2*t2**6)/t2**8*z_log(t1*t6/t2**2,1._ki)/12._ki
                !
                stemp6=t3**2*t5**2*t2**3/12._ki+t2**2*t1**2*t6**2*t5-t1**2*t6**2*t&
                  &3*t2/3._ki+t1*t6*t3**2*t2**2/2._ki+t5*t3**2*t6*t2**3/6._ki-t1*t6**&
                  &2*t3**2*t4/2._ki+t2*t6*t3**3*t4/2._ki+t2*t1*t6**2*t3**2/3._ki-t4*t&
                  &2**2*t1**2*t6**2/2._ki-t5*t1**2*t6**2*t3/3._ki-t4*t5**2*t6**2*t1*&
                  &*2/2._ki+2._ki*t2*t4**2*t3**2*t6**2-5._ki/4._ki*t2**2*t4*t3**2*t6**&
                  &2-t2**3*t4*t6*t3**2/2._ki+t2**2*t4**2*t6*t3**2/3._ki-t3**2*t2**3*&
                  &t4*t5/6._ki+t3**2*t2**3*t6**2/4._ki-4._ki/3._ki*t3*t4**2*t5*t6**2*t&
                  &1-t3*t2**2*t5**2*t1*t6/3._ki-2._ki/3._ki*t3*t5*t2**3*t1*t6-2._ki/3.&
                  &_ki*t3*t2**2*t6**2*t5*t1+t3**2*t2**5/12._ki
                !
                stemp5=stemp6-t2**3*t3**3/6._ki-2._ki/3._ki*t3*t2**3*t6**2*t1-t3*t2*&
                  &*4*t1*t6/3._ki+t2*t5**2*t1**2*t6**2/2._ki+t3**2*t4**2*t5*t2*t6/3.&
                  &_ki+t3*t4*t5**2*t2*t1*t6/3._ki+2._ki*t2**2*t4*t1*t6**2*t3+t3**2*t2&
                  &**4*t6/6._ki-t6*t3**3*t2**2/3._ki-t4**3*t3**2*t6**2-t4*t3**2*t2**&
                  &4/12._ki-t2**2*t5*t3**3/6._ki+t2**3*t1**2*t6**2/2._ki+t3**2*t2**4*&
                  &t5/6._ki-4._ki/3._ki*t2*t4**2*t1*t6**2*t3+2._ki/3._ki*t4*t3*t5*t2**2&
                  &*t1*t6-t4*t5*t2*t6**2*t1**2+t4*t3*t2**3*t1*t6/3._ki+t2*t5*t1*t6*&
                  &t3**2/2._ki+2._ki*t3*t4*t5*t2*t6**2*t1-t4*t5*t3**2*t6*t2**2/2._ki-&
                  &t3**2*t4*t5**2*t2**2/12._ki
                !
                stemp6=1._ki/t3**2/t2**5*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)
                !
                stemp4=stemp5*stemp6
                !
                stemp7=-t6**3*t3**2*t2**5/4._ki+3._ki/2._ki*t6**2*t1*t3*t2**3*t4**3-&
                  &3._ki*t6**2*t1*t3*t2**4*t4**2-t6*t1**2*t2**4*t5*t3/9._ki-44._ki/3.&
                  &_ki*t6**3*t1*t3*t4**4*t5-77._ki/2._ki*t6**3*t1*t2**2*t4*t3**2+20._k&
                  &i*t6**3*t1**2*t3*t5*t4**2+15._ki/2._ki*t6**3*t1**2*t3*t2**2*t5-3.&
                  &_ki/4._ki*t6**2*t1**2*t3*t2**3*t4+15._ki/8._ki*t6**2*t1*t3*t2**5*t4&
                  &-t6*t1**2*t2**5*t4*t5/3._ki-t6*t1**2*t2**4*t4*t5**2/6._ki+77._ki/6&
                  &._ki*t6**3*t1*t3*t2**3*t5*t4+110._ki/3._ki*t6**3*t1*t3*t4**3*t5*t2&
                  &+66._ki*t6**3*t1*t2*t4**2*t3**2-7._ki/18._ki*t6*t1*t2**6*t3*t5-33.&
                  &_ki*t6**3*t1*t3*t2**2*t5*t4**2-t1*t3**2*t2**6/36._ki+t6*t1**2*t2*&
                  &*7/6._ki-25._ki*t6**3*t1**2*t3*t4*t5*t2-11._ki/6._ki*t6**3*t1*t3*t5&
                  &*t2**4+7._ki/6._ki*t6*t1*t2**6*t4*t3-2._ki/3._ki*t6*t1*t4*t3**2*t2*&
                  &*4-7._ki/9._ki*t6*t1*t2**5*t4**2*t3+4._ki*t6**3*t3**2*t4**5
                !
                stemp6=stemp7-7._ki/36._ki*t1*t2**8*t3-7._ki/9._ki*t6*t1*t2**4*t3*t4*&
                  &*2*t5+7._ki/6._ki*t6*t1*t2**5*t3*t4*t5+5._ki/36._ki*t1*t2**6*t3*t4*&
                  &t5+7._ki/36._ki*t1*t2**5*t3*t4*t5**2-t6*t1**2*t2**5*t3/9._ki-3._ki/&
                  &8._ki*t6**2*t2**6*t1*t3+t6*t1**2*t5*t2**6/3._ki+t6*t1**2*t2**5*t5&
                  &**2/6._ki-7._ki/18._ki*t6*t1*t3*t2**7-t6**3*t3*t5*t1**3+25._ki/2._ki&
                  &*t6**3*t1**2*t4*t3**2-15._ki/2._ki*t6**3*t1**2*t2*t3**2-110._ki/3.&
                  &_ki*t6**3*t1*t4**3*t3**2+9._ki/4._ki*t6**3*t4*t3**2*t2**4+22._ki/3.&
                  &_ki*t6**3*t1*t3**2*t2**3+14._ki*t6**3*t3**2*t2**2*t4**3+t6**2*t1*&
                  &*2*t3*t2**4/2._ki-8._ki*t6**3*t3**2*t2**3*t4**2-12._ki*t6**3*t3**2&
                  &*t4**4*t2-t6*t1**2*t2**6*t4/6._ki+2._ki/9._ki*t1*t2**5*t3**2*t5+7.&
                  &_ki/36._ki*t1*t2**7*t3*t4+4._ki/9._ki*t6*t1*t2**5*t3**2-7._ki/36._ki*&
                  &t1*t2**6*t3*t5**2-5._ki/36._ki*t1*t3*t2**7*t5
                !
                stemp7=1._ki/t1/t2**8/t3
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(4)
                !
                stemp2=(-18._ki*t6**3*t1*t2+36._ki*t6**3*t1*t4-18._ki*t6**3*t2**2*t4&
                  &+3._ki*t5**3*t2**3+6._ki*t5**2*t2**4+3._ki*t5*t2**5+2._ki*t6**2*t2*&
                  &*4+t2**5*t6+3._ki*t6**3*t2**3-24._ki*t6**3*t4**3+36._ki*t6**3*t4**&
                  &2*t2-4._ki*t6**2*t1*t2**2-6._ki*t6**2*t2**2*t3+8._ki*t6**2*t2**2*t&
                  &4**2-8._ki*t6**2*t2**3*t4+3._ki*t6**2*t5*t2**3+3._ki*t6*t2**3*t5**&
                  &2+4._ki*t6*t2**4*t5-2._ki*t6*t2**4*t4-4._ki*t6*t3*t2**3-6._ki*t6**2&
                  &*t1*t5*t2-12._ki*t6**2*t5*t2**2*t4+12._ki*t6**2*t5*t2*t4**2+12._ki&
                  &*t6**2*t2*t3*t4-6._ki*t6*t4*t5**2*t2**2-8._ki*t6*t2**3*t4*t5-6._ki&
                  &*t6*t2**2*t3*t5)/t2**6*z_log(t1*t6/t2**2,1._ki)/12._ki
                !
                stemp6=t1**3*t6**3*t5**3+2._ki*t3**3*t6**3*t4**3-t3**3*t6*t2**5/12&
                  &._ki-3._ki/2._ki*t6**2*t5*t2*t1*t3**3-3._ki/2._ki*t2*t3*t1**2*t6**3*&
                  &t5**2+t2**2*t3**2*t1*t6**2*t5**2+4._ki*t3**2*t1*t6**3*t5*t4**2+t&
                  &2**2*t3**2*t1*t6*t5**3-3._ki/2._ki*t2*t3*t1**2*t6**2*t5**3+3._ki*t&
                  &3*t1**2*t6**3*t4*t5**2+t2**2*t3**2*t1*t6**3*t5-4._ki*t2*t3**2*t1&
                  &*t6**3*t5*t4-2._ki*t2*t3**2*t1*t6**2*t4*t5**2+t2**2*t3**3*t6**2*&
                  &t5*t4-t2*t3**3*t6**2*t5*t4**2+t2**2*t3**3*t6*t4*t5**2/2._ki+4._ki&
                  &*t6**3*t1**2*t3*t4*t5*t2+2._ki/3._ki*t3**3*t2**3*t6*t4*t5+t3*t2**&
                  &2*t6**3*t4*t1**2-2._ki/3._ki*t3**2*t2**3*t1*t6**2*t4-3._ki/2._ki*t3&
                  &*t2**3*t1**2*t6**2*t5+t3**2*t2**4*t1*t6*t5+2._ki*t5**2*t2*t1**3*&
                  &t6**3+2._ki/3._ki*t3**3*t2**3*t6**2*t4-t3**3*t2**4*t6*t5/3._ki-2._k&
                  &i/3._ki*t3**3*t2**2*t6**2*t4**2-t3*t2**3*t6**3*t1**2/2._ki+t3**2*&
                  &t2**4*t1*t6**2/3._ki
                !
                stemp5=stemp6+t3**3*t2**4*t6*t4/6._ki+t6*t5*t2**2*t3**4/2._ki-3._ki*&
                  &t2*t3**3*t6**3*t4**2+3._ki/2._ki*t2**2*t3**3*t6**3*t4+t3**3*t6**3&
                  &*t1*t4+t2**2*t1**3*t6**3*t5+t6**3*t5*t1**2*t3**2-t2**3*t3**3*t6&
                  &*t5**2/4._ki-t2**3*t3**3*t6**2*t5/4._ki+2._ki/3._ki*t6**3*t1**2*t2*&
                  &t3**2+2._ki/3._ki*t6**3*t1*t3**2*t2**3-t3**3*t2*t6**3*t1/2._ki-t3*&
                  &*4*t2*t6**2*t4-t6**2*t2**2*t1*t3**3-t3**3*t5**2*t2**4/2._ki-8._ki&
                  &/3._ki*t6**3*t1*t2**2*t4*t3**2-2._ki*t6**3*t1**2*t3*t2**2*t5+8._ki&
                  &/3._ki*t6**3*t1*t2*t4**2*t3**2-t3**3*t2**5*t5/4._ki-t2**3*t3**3*t&
                  &5**3/4._ki-t3**3*t2**4*t6**2/6._ki+4._ki/3._ki*t3**2*t2**3*t6**2*t5&
                  &*t1-8._ki/3._ki*t3**2*t2**2*t6**2*t4*t5*t1-3._ki*t3*t5**2*t2**2*t1&
                  &**2*t6**2+2._ki*t3**2*t5**2*t2**3*t1*t6+t6*t2**3*t3**4/3._ki+t3**&
                  &4*t2**2*t6**2/2._ki-t2**3*t3**3*t6**3/4._ki
                !
                stemp6=1._ki/t2**6/t3**3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)
                !
                stemp4=stemp5*stemp6
                !
                stemp7=11._ki/24._ki*t6**3*t2**3*t3**2-11._ki/3._ki*t6**3*t3**2*t4**3&
                  &-11._ki/12._ki*t6**2*t2**2*t3**3+17._ki/36._ki*t6**2*t3**2*t2**4+5.&
                  &_ki/18._ki*t6*t2**5*t3**2+t6**2*t1*t2*t5*t3**2/12._ki+4._ki/3._ki*t6&
                  &**2*t1*t3*t4*t5*t2**2+t6**2*t1*t2**3*t3*t4/3._ki-t6**2*t1*t5**2*&
                  &t3*t2**2/2._ki+t6**2*t1**2*t5*t2**3/3._ki+t6**2*t1**2*t2*t5**3/3.&
                  &_ki-t6**2*t1*t2**4*t3/6._ki-17._ki/9._ki*t6**2*t2**3*t3**2*t4+11._ki&
                  &/6._ki*t6**2*t2*t3**3*t4+17._ki/24._ki*t6**2*t2**3*t3**2*t5+17._ki/&
                  &9._ki*t6**2*t3**2*t4**2*t2**2-5._ki/4._ki*t6**3*t1*t3**2*t2+5._ki/2&
                  &._ki*t6**3*t3**2*t1*t4+11._ki/2._ki*t6**3*t3**2*t2*t4**2
                !
                stemp6=stemp7-11._ki/4._ki*t6**3*t2**2*t3**2*t4+2._ki/3._ki*t6**2*t1*&
                  &*2*t5**2*t2**2+t3**2*t2**2*t1*t6**2/18._ki-7._ki/6._ki*t6*t3**3*t2&
                  &**2*t5+5._ki/6._ki*t6*t3**2*t2**3*t5**2-5._ki/9._ki*t6*t2**4*t3**2*&
                  &t4+10._ki/9._ki*t6*t2**4*t3**2*t5+17._ki/8._ki*t3**2*t2**4*t5**2+17&
                  &._ki/6._ki*t6**2*t2*t3**2*t4**2*t5-17._ki/6._ki*t6**2*t3**2*t4*t5*t&
                  &2**2-4._ki/3._ki*t6*t1*t3*t2**3*t5**2-2._ki/3._ki*t6*t1*t5*t2**4*t3&
                  &-2._ki/3._ki*t6*t1*t3*t2**2*t5**3-20._ki/9._ki*t6*t3**2*t2**3*t4*t5&
                  &-5._ki/3._ki*t6*t3**2*t2**2*t4*t5**2-2._ki/3._ki*t6**2*t1*t2**3*t3*&
                  &t5+t6**2*t1*t2*t3*t4*t5**2+11._ki/12._ki*t2**5*t3**2*t5-7._ki/9._ki&
                  &*t3**3*t2**3*t6+11._ki/12._ki*t3**2*t2**3*t5**3
                !
                stemp7=1._ki/t3**2/t2**6
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
              !
            case(2)
              !
              select case(par4_glob)
              !
              case(2)
                !
                stemp5=95._ki*t6**3*t3**2*t4**3+7._ki/12._ki*t6**2*t4**2*t2**6+6._ki*&
                  &t6**2*t4*t5*t3*t2**4-30._ki*t6**2*t4**2*t5*t3*t2**3+36._ki*t6**2*&
                  &t4**3*t5*t3*t2**2+15._ki*t6**2*t4*t5**2*t3*t2**3+10._ki*t6**2*t1*&
                  &t4*t5**3*t2**2+6._ki*t6**2*t1*t4*t5**2*t2**3-25._ki*t6**2*t1*t4**&
                  &2*t5**3*t2-t6*t4**2*t2**7/4._ki-t4**2*t5*t2**7/12._ki-5._ki*t6**2*&
                  &t2**2*t3**3+20._ki*t6**3*t4**5*t5**2-t4*t3*t2**7/6._ki-t6**3*t4**&
                  &2*t2**5/4._ki+3._ki/2._ki*t6**3*t4**3*t2**4-3._ki/4._ki*t6**3*t1**2*&
                  &t2**3+2._ki*t6**3*t4**5*t2**2+7._ki/3._ki*t6**2*t4**4*t2**4-60._ki*&
                  &t6**2*t3**2*t4*t5*t2**2+3._ki/2._ki*t6**2*t3**2*t2**4-30._ki*t6**3&
                  &*t1*t4**2*t5*t2**2+6._ki*t6**3*t1*t4*t5*t2**3+75._ki*t6**3*t4**2*&
                  &t3*t5*t2**2+10._ki*t6**2*t1*t2**3*t3*t4+15._ki*t6**2*t1*t5**2*t3*&
                  &t2**2
                !
                stemp4=195._ki/2._ki*t6**2*t2*t3**2*t4**2*t5+90._ki*t6**2*t4**3*t5**&
                  &2*t3*t2-15._ki*t6**2*t1*t4**2*t5**2*t2**2-75._ki*t6**2*t4**2*t5**&
                  &2*t3*t2**2+6._ki*t6**2*t1*t2**3*t3*t5-15._ki*t6**2*t1*t2*t5*t3**2&
                  &-24._ki*t6**2*t1*t3*t4*t5*t2**2-t6*t4*t2**6*t1/2._ki+30._ki*t6**3*&
                  &t2**2*t3**2*t4+3._ki/2._ki*t6**2*t1**2*t5**2*t2**2-3._ki*t3**2*t2*&
                  &*2*t1*t6**2+15._ki/2._ki*t6**3*t1*t4**2*t2**3-60._ki*t6**2*t1*t2*t&
                  &3*t4*t5**2-10._ki*t6**3*t4*t3*t5*t2**3-15._ki*t6**3*t1*t3*t5*t2**&
                  &2-12._ki*t6**3*t1**2*t5*t2*t4+25._ki/2._ki*t6**2*t1*t4**2*t5*t2**3&
                  &+stemp5+15._ki*t6**2*t2*t3**3*t4+15._ki/2._ki*t6**2*t2**3*t3**2*t5&
                  &+39._ki/2._ki*t6**2*t3**2*t4**2*t2**2+15._ki*t6**3*t1*t3**2*t2-45.&
                  &_ki*t6**3*t3**2*t1*t4-5._ki/4._ki*t6**2*t4**2*t5*t2**5+5._ki*t6**2*&
                  &t4**3*t5*t2**4-5._ki*t6**2*t4**4*t5*t2**3-15._ki*t6**2*t4**3*t3*t&
                  &2**3
                !
                stemp5=stemp4-5._ki/2._ki*t6**2*t4*t3*t2**5-6._ki*t6**2*t4**3*t5**2*&
                  &t2**3+6._ki*t6**2*t4**4*t5**2*t2**2+3._ki/2._ki*t6**2*t4**2*t5**2*&
                  &t2**4-10._ki*t6**2*t4**3*t5**3*t2**2+10._ki*t6**2*t4**4*t5**3*t2+&
                  &5._ki/2._ki*t6**2*t4**2*t5**3*t2**3-35._ki/6._ki*t6**2*t1*t4**2*t2*&
                  &*4+7._ki/3._ki*t6**2*t1*t4*t2**5-5._ki/4._ki*t6**2*t1**2*t5*t2**3+5&
                  &._ki/2._ki*t6**2*t1**2*t2*t5**3-5._ki/2._ki*t6**2*t1*t2**4*t3-12._ki&
                  &*t6**2*t2**3*t3**2*t4-195._ki/2._ki*t6**3*t3**2*t2*t4**2-5._ki*t6*&
                  &*2*t1*t4*t5*t2**4-180._ki*t6**3*t4**3*t3*t5*t2+120._ki*t6**3*t1*t&
                  &3*t5*t2*t4-15._ki*t6**3*t1*t4*t5**2*t2**2+36._ki*t6**3*t1*t4**3*t&
                  &5*t2+39._ki*t6**3*t1*t2*t3*t4**2-195._ki*t6**3*t1*t4**2*t3*t5+75.&
                  &_ki*t6**3*t1*t4**2*t5**2*t2-24._ki*t6**3*t1*t2**2*t3*t4-t4**2*t2*&
                  &*8/12._ki-3._ki/2._ki*t6**3*t1*t4*t2**4
                !
                stemp3=stemp5-90._ki*t6**3*t1*t4**3*t5**2+3._ki*t6**3*t1*t2**3*t3+1&
                  &5._ki*t6**3*t1**2*t3*t5+3._ki*t6**3*t1**2*t2**2*t4-28._ki*t6**3*t4&
                  &**4*t3*t2-15._ki*t6**3*t4**2*t3*t2**3-30._ki*t6**3*t4**4*t5**2*t2&
                  &+36._ki*t6**3*t4**3*t3*t2**2+2._ki*t6**3*t4*t3*t2**4+15._ki*t6**3*&
                  &t4**3*t5**2*t2**2-9._ki*t6**3*t1*t4**3*t2**2-8._ki*t6**3*t4**5*t5&
                  &*t2-5._ki/2._ki*t6**3*t4**2*t5**2*t2**3-6._ki*t6**3*t4**3*t5*t2**3&
                  &+12._ki*t6**3*t4**4*t5*t2**2+t6**3*t4**2*t5*t2**4+140._ki*t6**3*t&
                  &4**4*t3*t5+30._ki*t6**3*t1**2*t5**2*t4-15._ki/2._ki*t6**3*t1**2*t5&
                  &**2*t2-3._ki*t6**3*t1**2*t2*t3+3._ki*t6**3*t1**2*t5*t2**2+25._ki/2&
                  &._ki*t6**2*t4**2*t3*t2**4-3._ki*t6**3*t4**4*t2**3+7._ki/12._ki*t1**&
                  &2*t2**4*t6**2-7._ki/3._ki*t6**2*t4**3*t2**5+t6*t4**3*t2**6/2._ki-5&
                  &._ki/2._ki*t6**3*t2**3*t3**2
                !
                stemp4=1._ki/t2**10*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=-(-t2**2*t3*t4+4._ki*t4*t2*t1*t6-3._ki*t6*t4*t2*t3+6._ki*t4**&
                  &2*t3*t6-2._ki*t3**2*t2+2._ki*t3*t1*t6+4._ki*t4*t1*t6*t5-t3*t5*t2*t&
                  &4)*t4/t2**4/t3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki
                !
                stemp9=-25._ki/2._ki*t6**3*t1**4*t4*t5**2*t2**2+22._ki*t6**3*t1**3*t&
                  &4**4*t5*t2**2-t1**3*t4**2*t2**8/9._ki+t6**2*t1**2*t4*t3**2*t2**5&
                  &-14._ki*t6**2*t1**2*t4**4*t3**2*t2**2-55._ki/3._ki*t6**3*t1**3*t4*&
                  &t3*t5*t2**3-45._ki*t6**3*t1**2*t4**3*t3**2*t2**2-8._ki*t6**3*t1**&
                  &2*t4**5*t3*t2**2-45._ki*t6**2*t1**3*t4**2*t5*t3*t2**3+54._ki*t6**&
                  &2*t1**3*t4**3*t5*t3*t2**2+15._ki*t6**2*t1**3*t4**4*t5**3*t2-5._ki&
                  &/4._ki*t6**3*t1**2*t4**2*t5*t2**4*t3-5._ki/12._ki*t6**2*t1*t2**5*t&
                  &3**3*t4-5._ki/8._ki*t6**2*t1*t4**2*t3**2*t2**5*t5-25._ki/2._ki*t6**&
                  &3*t1**4*t3*t5*t2**2+70._ki/3._ki*t6**2*t1*t4**4*t3**3*t2**2+25._ki&
                  &/6._ki*t6**2*t1*t2**4*t3**3*t4**2-15._ki/2._ki*t6**2*t1**4*t4**2*t&
                  &5**2*t2**2-5._ki/2._ki*t6**2*t1**4*t4*t5*t2**4-25._ki/2._ki*t6**2*t&
                  &1**4*t4**2*t5**3*t2
                !
                stemp8=stemp9+25._ki/4._ki*t6**2*t1**4*t4**2*t5*t2**3+5._ki*t6**2*t1&
                  &**4*t4*t5**3*t2**2+3._ki*t6**2*t1**4*t2**3*t3*t5-2._ki*t6**2*t1*t&
                  &4**6*t3**2*t2**2+4._ki*t6**2*t1*t4**5*t3**2*t2**3+t6**2*t1*t3**2&
                  &*t4**3*t2**5-10._ki*t6**3*t1*t4**5*t3**2*t2**2-4._ki*t6**3*t1**5*&
                  &t5*t2*t4+10._ki*t6**3*t1*t4**6*t3**2*t2+5._ki*t6**3*t1*t4**4*t3**&
                  &2*t2**3-t6**3*t1**5*t2**3/4._ki+15._ki/4._ki*t6**2*t1**3*t4**2*t5*&
                  &*3*t2**3+45._ki/2._ki*t6**2*t1**3*t4*t5**2*t3*t2**3-325._ki/2._ki*t&
                  &6**3*t1**4*t4**2*t3*t5+100._ki*t6**3*t1**4*t3*t5*t2*t4+15._ki/2._k&
                  &i*t6**2*t1**3*t4**3*t5*t2**4-15._ki/2._ki*t6**2*t1**3*t4**4*t5*t2&
                  &**3+5._ki*t6**2*t1*t4**3*t3**2*t5*t2**4-3._ki*t6**2*t1*t4**4*t3**&
                  &2*t2**4-95._ki/3._ki*t6**2*t1**2*t3**3*t4**3*t2+65._ki/2._ki*t6**2*&
                  &t1**2*t3**3*t4**2*t2**2
                !
                stemp9=20._ki*t6**2*t1*t4**5*t3**2*t5*t2**2-15._ki/8._ki*t6**2*t1**3&
                  &*t4**2*t5*t2**5-t6**2*t1*t4**2*t3**2*t2**6/8._ki+117._ki/4._ki*t6*&
                  &*2*t1**3*t3**2*t2**2*t4**2-18._ki*t6**2*t1**3*t3**2*t2**3*t4-15.&
                  &_ki*t6**2*t1**3*t4**3*t5**3*t2**2+125._ki/2._ki*t6**3*t1**4*t4**2*&
                  &t5**2*t2-6._ki*t6**2*t1**2*t4**3*t5*t3*t2**4-5._ki/4._ki*t6**3*t1*&
                  &t4**3*t3**2*t2**4-20._ki*t6**2*t1**2*t4**5*t3*t5**2*t2-70._ki*t6*&
                  &*2*t1**2*t3**2*t4**4*t5*t2+5._ki*t6**2*t1**4*t2**3*t3*t4-12._ki*t&
                  &6**2*t1**4*t2**2*t3*t4*t5-30._ki*t6**2*t1**4*t5**2*t3*t2*t4+55._k&
                  &i*t6**3*t1**3*t4*t3**2*t2**2+11._ki/3._ki*t6**3*t1**3*t4*t3*t2**4&
                  &+770._ki/3._ki*t6**3*t1**3*t4**4*t3*t5+11._ki/6._ki*t6**3*t1**3*t4*&
                  &*2*t5*t2**4-330._ki*t6**3*t1**3*t4**3*t3*t5*t2-715._ki/4._ki*t6**3&
                  &*t1**3*t3**2*t4**2*t2+9._ki*t6**2*t1**3*t4*t5*t3*t2**4
                !
                stemp7=stemp9+9._ki/4._ki*t6**2*t1**3*t4**2*t5**2*t2**4-44._ki/3._ki*&
                  &t6**3*t1**3*t4**5*t5*t2-15._ki*t6**2*t1**2*t4**3*t3*t5**2*t2**3+&
                  &10._ki/3._ki*t6**2*t1**2*t4**5*t3*t2**3-15._ki/2._ki*t6**2*t1**2*t4&
                  &**2*t3**2*t2**4+12._ki*t6**2*t1**2*t4**4*t5*t3*t2**3+90._ki*t6**2&
                  &*t1**2*t3**2*t4**3*t5*t2**2+3._ki*t6**2*t1**4*t4*t5**2*t2**3-15.&
                  &_ki/2._ki*t6**2*t1**4*t3**2*t5*t2+15._ki/2._ki*t6**2*t1**4*t5**2*t3&
                  &*t2**2+5._ki*t6**3*t1**5*t3*t5+10._ki*t6**3*t1**5*t5**2*t4-5._ki/2&
                  &._ki*t6**3*t1**5*t5**2*t2+t6**3*t1**5*t5*t2**2-t6**3*t1**5*t2*t3&
                  &-15._ki/2._ki*t6**3*t1**4*t4**3*t2**2+25._ki/4._ki*t6**3*t1**4*t4**&
                  &2*t2**3-5._ki/4._ki*t6**3*t1**4*t4*t2**4-75._ki/2._ki*t6**3*t1**4*t&
                  &3**2*t4-75._ki*t6**3*t1**4*t4**3*t5**2+stemp8
                !
                stemp9=stemp7+5._ki/2._ki*t6**3*t1**4*t2**3*t3+25._ki/2._ki*t6**3*t1*&
                  &*4*t3**2*t2+11._ki/4._ki*t6**3*t1**3*t4**3*t2**4+11._ki/3._ki*t6**3&
                  &*t1**3*t4**5*t2**2-55._ki/12._ki*t6**3*t1**3*t3**2*t2**3+110._ki/3&
                  &._ki*t6**3*t1**3*t4**5*t5**2+1045._ki/6._ki*t6**3*t1**3*t3**2*t4**&
                  &3-11._ki/24._ki*t6**3*t1**3*t4**2*t2**5-11._ki/2._ki*t6**3*t1**3*t4&
                  &**4*t2**3-15._ki/2._ki*t6**2*t1**3*t3**3*t2**2-35._ki/12._ki*t6**2*&
                  &t1**4*t4**2*t2**4-5._ki/4._ki*t6**2*t1**4*t2**4*t3-3._ki/2._ki*t6**&
                  &2*t1**4*t3**2*t2**2+5._ki/6._ki*t6**2*t1**2*t3**3*t2**4+10._ki/3._k&
                  &i*t6**2*t4**4*t3**3*t2**4-20._ki/3._ki*t6**2*t4**5*t3**3*t2**3+20&
                  &._ki/3._ki*t6**2*t4**6*t3**3*t2**2-8._ki/3._ki*t6**2*t4**7*t3**3*t2&
                  &+7._ki/6._ki*t6**2*t1**4*t4*t2**5+7._ki/8._ki*t6**2*t1**3*t4**2*t2*&
                  &*6
                !
                stemp8=stemp9+7._ki/2._ki*t6**2*t1**3*t4**4*t2**4-7._ki/2._ki*t6**2*t&
                  &1**3*t4**3*t2**5+9._ki/4._ki*t6**2*t1**3*t3**2*t2**4+t6*t1**3*t4*&
                  &*3*t2**6/2._ki-t6*t1**3*t4**2*t2**7/4._ki+t6**2*t2**6*t3**3*t4**2&
                  &/12._ki+t6**3*t1**5*t2**2*t4-10._ki*t6**2*t1*t4**6*t3**2*t5*t2-15&
                  &._ki*t6**2*t1*t4**4*t3**2*t5*t2**3-15._ki*t6**2*t1*t2**3*t3**3*t4&
                  &**3+5._ki/2._ki*t6**2*t1**2*t4**2*t3*t5**2*t2**4-8._ki*t6**2*t1**2&
                  &*t4**5*t5*t3*t2**2+5._ki*t6**2*t1**2*t4*t3**2*t5*t2**4-10._ki*t6*&
                  &*2*t1**2*t4*t3**3*t2**3-5._ki/12._ki*t6**2*t1**2*t4**2*t3*t2**6-5&
                  &._ki*t6**2*t1**2*t4**4*t3*t2**4+18._ki*t6**2*t1**2*t4**3*t3**2*t2&
                  &**3+5._ki/2._ki*t6**2*t1**2*t4**3*t3*t2**5-55._ki*t6**3*t1**3*t4**&
                  &4*t5**2*t2+55._ki/2._ki*t6**3*t1**3*t4**3*t5**2*t2**2+275._ki/2._ki&
                  &*t6**3*t1**3*t4**2*t3*t5*t2**2
                !
                stemp9=stemp8+40._ki*t6**3*t1**2*t4**5*t3*t5*t2+t6**3*t1**2*t4**2*&
                  &t3*t2**5/4._ki-5._ki/4._ki*t6**3*t1**2*t4*t3**2*t2**4-225._ki/2._ki*&
                  &t6**2*t1**3*t4**2*t5**2*t3*t2**2+135._ki*t6**2*t1**3*t4**3*t5**2&
                  &*t3*t2+5._ki*t6**3*t1**4*t4*t5*t2**3+45._ki/4._ki*t6**2*t1**3*t2**&
                  &3*t3**2*t5+75._ki/4._ki*t6**2*t1**3*t4**2*t3*t2**4-55._ki/12._ki*t6&
                  &**3*t1**3*t4**2*t5**2*t2**3-55._ki/2._ki*t6**3*t1**3*t4**2*t3*t2*&
                  &*3+66._ki*t6**3*t1**3*t4**3*t3*t2**2-154._ki/3._ki*t6**3*t1**3*t4*&
                  &*4*t3*t2+t6**3*t1*t4**2*t3**2*t2**5/8._ki+30._ki*t6**3*t1**4*t4**&
                  &3*t5*t2+65._ki/2._ki*t6**3*t1**4*t2*t3*t4**2-11._ki*t6**3*t1**3*t4&
                  &**3*t5*t2**3-5._ki/6._ki*t6**2*t2**5*t3**3*t4**3-t1**3*t4**2*t5*t&
                  &2**7/9._ki+585._ki/4._ki*t6**2*t1**3*t3**2*t5*t2*t4**2-90._ki*t6**2&
                  &*t1**3*t3**2*t5*t2**2*t4
                !
                stemp6=stemp9-45._ki/2._ki*t6**2*t1**3*t4**3*t3*t2**3-15._ki/4._ki*t6&
                  &**2*t1**3*t4*t3*t2**5+45._ki/2._ki*t6**2*t1**3*t3**3*t2*t4-9._ki*t&
                  &6**2*t1**3*t4**3*t5**2*t2**3-75._ki/2._ki*t6**2*t1**2*t3**2*t4**2&
                  &*t5*t2**3+30._ki*t6**2*t1**2*t4**4*t3*t5**2*t2**2-40._ki/3._ki*t6*&
                  &*2*t1*t4**5*t3**3*t2-2._ki*t6**3*t1**2*t4**3*t3*t2**4+25._ki/2._ki&
                  &*t6**3*t1**2*t4**2*t3**2*t2**3-20._ki*t6**3*t1**2*t4**6*t3*t5+4.&
                  &_ki*t6**3*t1**2*t4**6*t3*t2+6._ki*t6**3*t1**2*t4**4*t3*t2**3+t6**&
                  &2*t1**2*t4**2*t5*t3*t2**5-30._ki*t6**3*t1**2*t4**4*t3*t5*t2**2+1&
                  &0._ki*t6**3*t1**2*t4**3*t3*t5*t2**3-20._ki*t6**3*t1**4*t2**2*t3*t&
                  &4-25._ki*t6**3*t1**4*t4**2*t5*t2**2-t1**3*t4*t3*t2**7/18._ki-40._k&
                  &i*t6**3*t1**2*t4**5*t3**2-4._ki*t6**3*t1*t4**7*t3**2+9._ki*t6**2*&
                  &t1**3*t4**4*t5**2*t2**2+70._ki*t6**3*t1**2*t4**4*t3**2*t2
                !
                stemp7=1._ki/t1**3/t2**10
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(3)
                !
                stemp5=54._ki*t6**3*t1*t4**2*t5*t2**2-24._ki*t6**3*t1*t4*t5*t2**3-t&
                  &6*t1*t2**7/4._ki+7._ki/12._ki*t6**2*t4*t2**7+3._ki/4._ki*t6*t4**2*t2&
                  &**7+t6**3*t2**5*t3-3._ki/4._ki*t6**3*t1*t2**5-t6**3*t4*t2**6/4._ki&
                  &+t4**2*t2**8/12._ki+10._ki*t6**2*t4*t3*t2**5+12._ki*t6**2*t4**3*t5&
                  &**2*t2**3+3._ki*t6**3*t1*t5*t2**4-15._ki/2._ki*t6**3*t1*t5**2*t2**&
                  &3-5._ki/4._ki*t6**2*t2**6*t4*t5+3._ki/2._ki*t6**2*t4*t5**2*t2**5+5.&
                  &_ki/2._ki*t6**2*t4*t5**3*t2**4+10._ki*t6**3*t2**3*t3**2+3._ki/2._ki*&
                  &t6**3*t1**2*t2**3-9._ki/2._ki*t6**3*t4**3*t2**4-7._ki/3._ki*t6**2*t&
                  &4**4*t2**4-2._ki*t6**3*t4**5*t2**2+45._ki*t6**3*t3**2*t1*t4-135._k&
                  &i/2._ki*t6**3*t2**2*t3**2*t4+3._ki*t3**2*t2**2*t1*t6**2-3._ki/2._ki&
                  &*t6**2*t1**2*t5**2*t2**2-30._ki*t6**2*t1*t5**2*t3*t2**2+280._ki*t&
                  &6**3*t4**3*t3*t5*t2-195._ki*t6**3*t4**2*t3*t5*t2**2-195._ki/2._ki*&
                  &t6**2*t2*t3**2*t4**2*t5-195._ki*t6**3*t1*t3*t5*t2*t4+195._ki/2._ki&
                  &*t6**2*t3**2*t4*t5*t2**2+60._ki*t6**3*t1*t4*t5**2*t2**2
                !
                stemp4=-25._ki/2._ki*t6**2*t1*t4**2*t5*t2**3-20._ki*t6**3*t4**5*t5**&
                  &2-95._ki*t6**3*t3**2*t4**3+stemp5-36._ki*t6**3*t1*t4**3*t5*t2-39.&
                  &_ki*t6**3*t1*t2*t3*t4**2-15._ki*t6**2*t1*t4*t5**2*t2**3-25._ki*t6*&
                  &*2*t1*t4*t5**3*t2**2-t4*t2**9/12._ki+15._ki*t6**2*t1*t4**2*t5**2*&
                  &t2**2-24._ki*t6**2*t4*t5*t3*t2**4-15._ki/2._ki*t6**2*t4**2*t5**2*t&
                  &2**4+20._ki*t6**2*t4**3*t5**3*t2**2-10._ki*t6**2*t4**4*t5**3*t2-2&
                  &5._ki/2._ki*t6**2*t4**2*t5**3*t2**3+t6**3*t4*t5*t2**5+35._ki/6._ki*&
                  &t6**2*t1*t4**2*t2**4+t6*t4*t2**6*t1/2._ki+3._ki*t6**2*t3*t5*t2**5&
                  &+6._ki*t6**3*t1*t4*t2**4+90._ki*t6**3*t1*t4**3*t5**2-9._ki*t6**3*t&
                  &1*t2**3*t3-15._ki*t6**3*t1**2*t3*t5-3._ki*t6**3*t1**2*t2**2*t4-27&
                  &._ki/2._ki*t6**3*t1*t4**2*t2**3+135._ki*t6**2*t4**2*t5**2*t3*t2**2&
                  &+25._ki/2._ki*t6**2*t1*t4*t5*t2**4+12._ki*t6**3*t1**2*t5*t2*t4-90.&
                  &_ki*t6**2*t4**3*t5**2*t3*t2-140._ki*t6**3*t4**4*t3*t5-30._ki*t6**3&
                  &*t1**2*t5**2*t4+15._ki*t6**3*t1**2*t5**2*t2+3._ki*t6**3*t1**2*t2*&
                  &t3
                !
                stemp5=stemp4-5._ki*t6**3*t3*t2**4*t5-5._ki/2._ki*t6**3*t4*t5**2*t2*&
                  &*4+25._ki*t6**2*t1*t4**2*t5**3*t2-6._ki*t6**3*t1**2*t5*t2**2+25._k&
                  &i/4._ki*t6**2*t4**2*t5*t2**5-10._ki*t6**2*t4**3*t5*t2**4+5._ki*t6*&
                  &*2*t4**4*t5*t2**3+15._ki*t6**2*t4**3*t3*t2**3-45._ki/2._ki*t6**2*t&
                  &4**2*t3*t2**4+3._ki*t6**2*t1*t5**2*t2**4-5._ki/2._ki*t6**2*t1*t5*t&
                  &2**5+5._ki*t6**2*t1*t5**3*t2**3+15._ki/2._ki*t6**2*t3*t5**2*t2**4-&
                  &6._ki*t6**2*t4**4*t5**2*t2**2+39._ki/2._ki*t6**2*t2**3*t3**2*t4+28&
                  &5._ki/2._ki*t6**3*t3**2*t2*t4**2-15._ki*t6**2*t2*t3**3*t4-45._ki/2.&
                  &_ki*t6**2*t2**3*t3**2*t5-39._ki/2._ki*t6**2*t3**2*t4**2*t2**2-45._k&
                  &i/2._ki*t6**3*t1*t3**2*t2+54._ki*t6**2*t4**2*t5*t3*t2**3+60._ki*t6&
                  &**2*t1*t2*t3*t4*t5**2-36._ki*t6**2*t4**3*t5*t3*t2**2+28._ki*t6**3&
                  &*t4**4*t3*t2+50._ki*t6**3*t4**4*t5**2*t2-56._ki*t6**3*t4**3*t3*t2&
                  &**2-11._ki*t6**3*t4*t3*t2**4-45._ki*t6**3*t4**3*t5**2*t2**2+39._ki&
                  &*t6**3*t4**2*t3*t2**3-60._ki*t6**2*t4*t5**2*t3*t2**3+t4**2*t5*t2&
                  &**7/12._ki
                !
                stemp3=stemp5+15._ki*t6**2*t1*t2*t5*t3**2+55._ki*t6**3*t4*t3*t5*t2*&
                  &*3+15._ki/2._ki*t6**2*t2**2*t3**3+t4*t3*t2**7/6._ki+45._ki*t6**3*t1&
                  &*t3*t5*t2**2-5._ki/4._ki*t6**2*t3*t2**6+7._ki/6._ki*t6**2*t1*t2**6+&
                  &9._ki*t6**3*t1*t4**3*t2**2+8._ki*t6**3*t4**5*t5*t2+35._ki/2._ki*t6*&
                  &*3*t4**2*t5**2*t2**3+18._ki*t6**3*t4**3*t5*t2**3-20._ki*t6**3*t4*&
                  &*4*t5*t2**2-7._ki*t6**3*t4**2*t5*t2**4+195._ki*t6**3*t1*t4**2*t3*&
                  &t5+24._ki*t6**2*t1*t3*t4*t5*t2**2-t6*t4**3*t2**6/2._ki-10._ki*t6**&
                  &2*t1*t2**3*t3*t4-35._ki/6._ki*t6**2*t1*t4*t2**5-135._ki*t6**3*t1*t&
                  &4**2*t5**2*t2+39._ki*t6**3*t1*t2**2*t3*t4-t4*t2**8*t5/12._ki-9._ki&
                  &/2._ki*t6**2*t3**2*t2**4+5._ki/4._ki*t6**2*t1**2*t5*t2**3-5._ki/2._k&
                  &i*t6**2*t1**2*t2*t5**3+5._ki*t6**2*t1*t2**4*t3-35._ki/12._ki*t6**2&
                  &*t4**2*t2**6+5._ki*t6**3*t4**4*t2**3-t3*t2**8/12._ki+7._ki/4._ki*t6&
                  &**3*t4**2*t2**5+14._ki/3._ki*t6**2*t4**3*t2**5-12._ki*t6**2*t1*t2*&
                  &*3*t3*t5-7._ki/12._ki*t1**2*t2**4*t6**2-t6*t4*t2**8/4._ki
                !
                stemp4=1._ki/t2**10*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=-(t4**2*t3*t2**2-4._ki*t4**2*t2*t1*t6+9._ki*t4**2*t3*t6*t2+t&
                  &3*t4**2*t2*t5+4._ki*t4*t2**2*t1*t6+2._ki*t4*t3**2*t2-t4*t3*t2**3-&
                  &t3**2*t2**2+4._ki*t4*t2*t1*t6*t5-4._ki*t4**2*t5*t1*t6-2._ki*t4*t3*&
                  &t1*t6-6._ki*t3*t4**3*t6+t3*t2*t1*t6-t4*t3*t5*t2**2-3._ki*t4*t3*t6&
                  &*t2**2)/t2**4/t3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki
                !
                stemp9=143._ki/2._ki*t6**3*t1**3*t4**2*t3*t2**3+385._ki/12._ki*t6**3*&
                  &t1**3*t4**2*t5**2*t2**3-308._ki/3._ki*t6**3*t1**3*t4**3*t3*t2**2-&
                  &5._ki/4._ki*t6**2*t1**4*t5*t2**5-45._ki/4._ki*t6**3*t1**4*t4**2*t2*&
                  &*3+5._ki*t6**3*t1**4*t4*t2**4+75._ki/2._ki*t6**3*t1**4*t3**2*t4+75&
                  &._ki*t6**3*t1**4*t4**3*t5**2-15._ki/2._ki*t6**3*t1**4*t2**3*t3-585&
                  &._ki/4._ki*t6**2*t1**3*t3**2*t5*t2*t4**2+45._ki/2._ki*t6**2*t1**3*t&
                  &4**3*t3*t2**3+39._ki/2._ki*t6**2*t1**2*t4**2*t3**2*t2**4-20._ki*t6&
                  &**2*t1**2*t4**4*t5*t3*t2**3+15._ki/2._ki*t6**2*t1**4*t3**2*t5*t2+&
                  &45._ki/4._ki*t6**2*t1**3*t3*t5**2*t2**4+585._ki/4._ki*t6**2*t1**3*t&
                  &3**2*t5*t2**2*t4-15._ki*t6**2*t1**3*t4**4*t5**3*t2+15._ki*t6**2*t&
                  &1**3*t4*t3*t2**5-45._ki/2._ki*t6**2*t1**3*t3**3*t2*t4-5._ki/12._ki*&
                  &t6**2*t1**2*t4*t3*t2**7+154._ki/3._ki*t6**3*t1**3*t4**4*t3*t2-15.&
                  &_ki*t6**2*t1**4*t5**2*t3*t2**2+18._ki*t6**2*t1**3*t4**3*t5**2*t2*&
                  &*3-55._ki/12._ki*t6**3*t1**3*t4*t5**2*t2**4+5._ki/2._ki*t6**2*t1**2&
                  &*t2**5*t3**2*t5
                !
                stemp8=-135._ki/4._ki*t6**2*t1**3*t4**2*t3*t2**4-135._ki/4._ki*t6**2*&
                  &t1**3*t2**3*t3**2*t5-60._ki*t6**3*t1**2*t4**5*t3*t5*t2+75._ki/8._k&
                  &i*t6**2*t1**3*t4**2*t5*t2**5-140._ki*t6**2*t1**2*t3**2*t4**3*t5*&
                  &t2**2-15._ki*t6**2*t1**3*t4**3*t5*t2**4+15._ki/2._ki*t6**2*t1**3*t&
                  &4**4*t5*t2**3-9._ki*t6**2*t1**3*t4**4*t5**2*t2**2+195._ki/2._ki*t6&
                  &**2*t1**2*t3**2*t4**2*t5*t2**3+70._ki*t6**3*t1**2*t4**4*t3*t5*t2&
                  &**2+9._ki/2._ki*t6**2*t1**3*t3*t5*t2**5+20._ki*t6**3*t1*t4**5*t3**&
                  &2*t2**2-14._ki*t6**3*t1**2*t4**4*t3*t2**3-100._ki*t6**3*t1**2*t4*&
                  &*4*t3**2*t2-5._ki/4._ki*t6**3*t1**2*t4*t5*t2**5*t3+t6**3*t1**2*t4&
                  &*t3*t2**6/4._ki-40._ki*t6**3*t1**2*t4**3*t3*t5*t2**3-11._ki/2._ki*t&
                  &6**2*t1**2*t4*t3**2*t2**5-6._ki*t6**2*t1**4*t2**3*t3*t5-325._ki/2&
                  &._ki*t6**3*t1**4*t3*t5*t2*t4-85._ki/6._ki*t6**2*t1*t2**4*t3**3*t4*&
                  &*2-20._ki*t6**2*t1*t4**3*t3**2*t5*t2**4-5._ki*t6**2*t1**4*t2**3*t&
                  &3*t4+35._ki/12._ki*t6**2*t1*t2**5*t3**3*t4-75._ki/4._ki*t6**2*t1**3&
                  &*t4**2*t5**3*t2**3+stemp9
                !
                stemp9=-50._ki*t6**2*t1**2*t4**4*t3*t5**2*t2**2+14._ki*t6**2*t1**2*&
                  &t4**4*t3**2*t2**2-90._ki*t6**2*t1**3*t4*t5**2*t3*t2**3+40._ki/3._k&
                  &i*t6**2*t1*t4**5*t3**3*t2+12._ki*t6**2*t1**4*t2**2*t3*t4*t5+9._ki&
                  &/4._ki*t6**2*t1**3*t4*t5**2*t2**5+7._ki*t6**2*t1**3*t4**3*t2**5-2&
                  &7._ki/4._ki*t6**2*t1**3*t3**2*t2**4-t6*t1**3*t4**3*t2**6/2._ki+3._k&
                  &i/4._ki*t6*t1**3*t4**2*t2**7-25._ki/4._ki*t6**3*t1**4*t5**2*t2**3-&
                  &15._ki/8._ki*t6**2*t1**3*t3*t2**6+3._ki/2._ki*t6**2*t1**4*t5**2*t2*&
                  &*4+5._ki/2._ki*t6**2*t1**4*t5**3*t2**3+45._ki*t6**2*t1**2*t4**3*t3&
                  &*t5**2*t2**3-75._ki/4._ki*t6**3*t1**4*t3**2*t2-t6**3*t1**5*t2**2*&
                  &t4-10._ki*t6**3*t1**5*t5**2*t4+5._ki*t6**3*t1**5*t5**2*t2-2._ki*t6&
                  &**3*t1**5*t5*t2**2+t6**3*t1**5*t2*t3+15._ki/2._ki*t6**3*t1**4*t4*&
                  &*3*t2**2+40._ki*t6**3*t1**2*t4**5*t3**2+30._ki*t6**2*t1**4*t5**2*&
                  &t3*t2*t4-11._ki/12._ki*t6**2*t2**6*t3**3*t4**2
                !
                stemp7=stemp9-5._ki/24._ki*t6**2*t1*t2**6*t3**3+25._ki/6._ki*t6**2*t2&
                  &**5*t3**3*t4**3+11._ki/6._ki*t6**3*t1**3*t2**5*t3+95._ki/3._ki*t6**&
                  &2*t1**2*t3**3*t4**3*t2+7._ki*t6**2*t1*t4**4*t3**2*t2**4+45._ki/4.&
                  &_ki*t6**3*t1**2*t4**2*t5*t2**4*t3+t1**3*t4*t3*t2**7/18._ki-t1**3*&
                  &t4*t2**8*t5/9._ki+5._ki/2._ki*t6**3*t1**4*t5*t2**4-10._ki/3._ki*t6**&
                  &2*t1**2*t3**3*t2**4-5._ki/8._ki*t6**2*t1*t4*t3**2*t2**6*t5+2._ki*t&
                  &6**2*t1*t4**6*t3**2*t2**2+4._ki*t6**3*t1**5*t5*t2*t4-35._ki/2._ki*&
                  &t6**2*t1**2*t4**2*t3*t5**2*t2**4-5._ki/8._ki*t6**3*t1**4*t2**5+8.&
                  &_ki*t6**2*t1**2*t4**5*t5*t3*t2**2-11._ki/8._ki*t6**3*t1*t4**2*t3**&
                  &2*t2**5-6._ki*t6**2*t1*t4**5*t3**2*t2**3-55._ki/6._ki*t6**3*t1**3*&
                  &t3*t2**4*t5-15._ki/2._ki*t6**2*t1**2*t4**3*t3*t2**5-95._ki/2._ki*t6&
                  &**2*t1**2*t3**3*t4**2*t2**2-t6**2*t1*t4*t3**2*t2**7/8._ki-14._ki*&
                  &t6**3*t1*t4**6*t3**2*t2+18._ki*t6**2*t1**2*t4**3*t5*t3*t2**4-7._k&
                  &i*t6**2*t1**2*t4**2*t5*t3*t2**5+stemp8
                !
                stemp9=stemp7+4._ki*t6**3*t1*t4**7*t3**2-5._ki*t6**3*t1**5*t3*t5+7.&
                  &_ki/8._ki*t6**2*t1**3*t4*t2**7-t6*t1**3*t4*t2**8/4._ki-20._ki*t6**3&
                  &*t1**4*t4*t5*t2**3-30._ki*t6**2*t1*t4**5*t3**2*t5*t2**2-30._ki*t6&
                  &**3*t1**4*t4**3*t5*t2+605._ki/6._ki*t6**3*t1**3*t4*t3*t5*t2**3+45&
                  &._ki/8._ki*t6**2*t1*t4**2*t3**2*t2**5*t5+9._ki/8._ki*t6**2*t1*t4**2&
                  &*t3**2*t2**6+33._ki*t6**3*t1**3*t4**3*t5*t2**3-65._ki/2._ki*t6**3*&
                  &t1**4*t2*t3*t4**2+t6**3*t1*t4*t3**2*t2**6/8._ki+50._ki*t6**3*t1**&
                  &4*t4*t5**2*t2**2+25._ki/4._ki*t6**3*t1*t4**3*t3**2*t2**4+95._ki*t6&
                  &**3*t1**2*t4**3*t3**2*t2**2+10._ki*t6**2*t1*t4**6*t3**2*t5*t2+20&
                  &._ki*t6**2*t1**2*t4**5*t3*t5**2*t2-15._ki*t6**3*t1*t4**4*t3**2*t2&
                  &**3-t1**3*t4*t2**9/9._ki-10._ki*t6**2*t4**4*t3**3*t2**4+40._ki/3._k&
                  &i*t6**2*t4**5*t3**3*t2**3-28._ki/3._ki*t6**2*t4**6*t3**3*t2**2+8.&
                  &_ki/3._ki*t6**2*t4**7*t3**3*t2
                !
                stemp8=stemp9-35._ki/12._ki*t6**2*t1**4*t4*t2**5-35._ki/8._ki*t6**2*t&
                  &1**3*t4**2*t2**6-7._ki/2._ki*t6**2*t1**3*t4**4*t2**4-55._ki/2._ki*t&
                  &6**2*t1**2*t4*t3**2*t5*t2**4-4._ki*t6**2*t1*t3**2*t4**3*t2**5+27&
                  &5._ki/3._ki*t6**3*t1**3*t4**4*t5**2*t2-100._ki/3._ki*t6**2*t1*t4**4&
                  &*t3**3*t2**2+44._ki/3._ki*t6**3*t1**3*t4**5*t5*t2-110._ki/3._ki*t6*&
                  &*3*t1**3*t4**4*t5*t2**2-495._ki/4._ki*t6**3*t1**3*t4*t3**2*t2**2+&
                  &35._ki*t6**2*t1*t4**4*t3**2*t5*t2**3+15._ki/2._ki*t6**2*t1**4*t4**&
                  &2*t5**2*t2**2-36._ki*t6**2*t1**3*t4*t5*t3*t2**4+95._ki/3._ki*t6**2&
                  &*t1*t2**3*t3**3*t4**3+75._ki/2._ki*t6**3*t1**4*t3*t5*t2**2+25._ki/&
                  &4._ki*t6**2*t1**4*t4*t5*t2**4-165._ki/2._ki*t6**3*t1**3*t4**3*t5**&
                  &2*t2**2+35._ki/12._ki*t6**2*t1**2*t4**2*t3*t2**6+45._ki/2._ki*t6**2&
                  &*t1**2*t4*t3**3*t2**3+405._ki/2._ki*t6**2*t1**3*t4**2*t5**2*t3*t2&
                  &**2-121._ki/6._ki*t6**3*t1**3*t4*t3*t2**4-135._ki*t6**2*t1**3*t4**&
                  &3*t5**2*t3*t2-770._ki/3._ki*t6**3*t1**3*t4**4*t3*t5+25._ki/2._ki*t6&
                  &**2*t1**4*t4**2*t5**3*t2-9._ki/4._ki*t6**3*t1**2*t4**2*t3*t2**5+1&
                  &2._ki*t6**3*t1**2*t4**5*t3*t2**2
                !
                stemp9=stemp8+25._ki/3._ki*t6**2*t1**2*t4**4*t3*t2**4-225._ki/2._ki*t&
                  &6**3*t1**4*t4**2*t5**2*t2+65._ki/2._ki*t6**3*t1**4*t2**2*t3*t4+55&
                  &._ki/6._ki*t6**3*t1**3*t4**4*t2**3+35._ki/12._ki*t6**2*t1**4*t4**2*&
                  &t2**4+5._ki/2._ki*t6**2*t1**4*t2**4*t3+3._ki/2._ki*t6**2*t1**4*t3**&
                  &2*t2**2-5._ki/8._ki*t6**3*t1**2*t3**2*t2**5+t6**2*t1**2*t3**2*t2*&
                  &*6/2._ki+t6**2*t2**7*t3**3*t4/12._ki+t1**3*t4**2*t5*t2**7/9._ki+45&
                  &._ki/4._ki*t6**2*t1**3*t3**3*t2**2-25._ki/4._ki*t6**2*t1**4*t4**2*t&
                  &5*t2**3+7._ki/12._ki*t6**2*t1**4*t2**6+20._ki*t6**3*t1**2*t4**6*t3&
                  &*t5-45._ki/4._ki*t6**2*t1**3*t4**2*t5**2*t2**4-t1**3*t3*t2**8/36.&
                  &_ki-117._ki/4._ki*t6**2*t1**3*t3**2*t2**2*t4**2-28._ki*t6**2*t1**2*&
                  &t4**3*t3**2*t2**3+35._ki/4._ki*t6**3*t1**2*t4*t3**2*t2**4-25._ki/2&
                  &._ki*t6**2*t1**4*t4*t5**3*t2**2+45._ki*t6**3*t1**4*t4**2*t5*t2**2&
                  &-77._ki/6._ki*t6**3*t1**3*t4**2*t5*t2**4+1540._ki/3._ki*t6**3*t1**3&
                  &*t4**3*t3*t5*t2+t6**3*t1**5*t2**3/2._ki
                !
                stemp6=stemp9+117._ki/4._ki*t6**2*t1**3*t3**2*t2**3*t4+30._ki*t6**2*&
                  &t1**3*t4**3*t5**3*t2**2+t1**3*t4**2*t2**8/9._ki+8._ki*t6**3*t1**2&
                  &*t4**3*t3*t2**4-85._ki/2._ki*t6**3*t1**2*t4**2*t3**2*t2**3+70._ki*&
                  &t6**2*t1**2*t3**2*t4**4*t5*t2-10._ki/3._ki*t6**2*t1**2*t4**5*t3*t&
                  &2**3-715._ki/2._ki*t6**3*t1**3*t4**2*t3*t5*t2**2-15._ki/2._ki*t6**2&
                  &*t1**4*t4*t5**2*t2**3+11._ki/6._ki*t6**3*t1**3*t4*t5*t2**5+t6**2*&
                  &t1**2*t2**6*t3*t4*t5+325._ki/2._ki*t6**3*t1**4*t4**2*t3*t5+81._ki*&
                  &t6**2*t1**3*t4**2*t5*t3*t2**3-4._ki*t6**3*t1**2*t4**6*t3*t2-54._k&
                  &i*t6**2*t1**3*t4**3*t5*t3*t2**2-15._ki/8._ki*t6**2*t1**3*t2**6*t4&
                  &*t5+1045._ki/4._ki*t6**3*t1**3*t3**2*t4**2*t2-11._ki/24._ki*t6**3*t&
                  &1**3*t4*t2**6-33._ki/4._ki*t6**3*t1**3*t4**3*t2**4-11._ki/3._ki*t6*&
                  &*3*t1**3*t4**5*t2**2+55._ki/3._ki*t6**3*t1**3*t3**2*t2**3-110._ki/&
                  &3._ki*t6**3*t1**3*t4**5*t5**2-1045._ki/6._ki*t6**3*t1**3*t3**2*t4*&
                  &*3+77._ki/24._ki*t6**3*t1**3*t4**2*t2**5+5._ki/2._ki*t6**2*t1**2*t3&
                  &*t5**2*t2**5*t4+15._ki/4._ki*t6**2*t1**3*t4*t5**3*t2**4
                !
                stemp7=1._ki/t1**3/t2**10
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(4)
                !
                stemp2=(-12._ki*t6**2*t2**4*t4**2+12._ki*t6**2*t2**3*t4**3+t6*t2**6&
                  &*t4+3._ki*t6**2*t2**5*t4-120._ki*t6**3*t3*t4**3+6._ki*t6**3*t3*t2*&
                  &*3+2._ki*t6*t2**5*t3+2._ki*t3*t5*t2**5-18._ki*t6**2*t1*t2**3*t4+2.&
                  &_ki*t6*t2**5*t4*t5-4._ki*t6*t2**4*t4**2*t5-6._ki*t3*t4*t6*t2**4+14&
                  &4._ki*t6**3*t1*t5*t4**2+18._ki*t1*t6**3*t2**2*t5+6._ki*t6**3*t2**3&
                  &*t5*t4-36._ki*t6**3*t2**2*t5*t4**2+72._ki*t6**3*t4**3*t5*t2-36._ki&
                  &*t6**3*t3*t2*t1+90._ki*t6**3*t4*t1*t3-54._ki*t6**3*t3*t2**2*t4+2.&
                  &_ki*t6*t1*t2**4*t5+144._ki*t6**3*t3*t2*t4**2-108._ki*t1*t6**3*t4*t&
                  &5*t2+t3*t2**6+t2**6*t4*t5-2._ki*t6*t2**5*t4**2+t4*t5**2*t2**5+6.&
                  &_ki*t6**2*t1*t2**4-18._ki*t6**3*t1**2*t5+t6*t1*t2**5-48._ki*t6**3*&
                  &t4**4*t5)/t2**8*z_log(t1*t6/t2**2,1._ki)/12._ki
                !
                stemp4=(-t3**2*t4*t5**2*t2**2-2._ki*t6*t3**3*t2**2-2._ki*t2**2*t5*t&
                  &3**3-t3**2*t2**3*t4*t5-2._ki*t1**2*t6**2*t3*t2+3._ki*t1*t6*t3**2*&
                  &t2**2-6._ki*t1*t6**2*t3**2*t4+6._ki*t2*t6*t3**3*t4+2._ki*t2*t1*t6*&
                  &*2*t3**2-4._ki*t5*t1**2*t6**2*t3-6._ki*t4*t5**2*t6**2*t1**2+12._ki&
                  &*t2*t4**2*t3**2*t6**2-3._ki*t2**2*t4*t3**2*t6**2-t2**3*t4*t6*t3*&
                  &*2+2._ki*t2**2*t4**2*t6*t3**2+4._ki*t4*t3*t5*t2**2*t1*t6-6._ki*t4*&
                  &t5*t2*t6**2*t1**2+6._ki*t2*t5*t1*t6*t3**2+8._ki*t3*t4*t5*t2*t6**2&
                  &*t1-12._ki*t4**3*t3**2*t6**2-16._ki*t3*t4**2*t5*t6**2*t1+4._ki*t3*&
                  &*2*t4**2*t5*t2*t6+4._ki*t3*t4*t5**2*t2*t1*t6+4._ki*t2**2*t4*t1*t6&
                  &**2*t3-8._ki*t2*t4**2*t1*t6**2*t3-2._ki*t4*t5*t3**2*t6*t2**2-t2**&
                  &3*t3**3)/t3**2/t2**5*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki
                !
                stemp7=3._ki/2._ki*t6**2*t1*t3*t2**3*t4**3-3._ki/2._ki*t6**2*t1*t3*t2&
                  &**4*t4**2-3._ki/4._ki*t6**2*t1**2*t3*t2**3*t4+t1*t3**2*t2**6/9._ki&
                  &+2._ki*t6**3*t3**2*t4**5+7._ki/36._ki*t6*t1*t2**6*t4*t3+7._ki/36._ki&
                  &*t1*t2**6*t3*t4*t5+7._ki/36._ki*t1*t2**5*t3*t4*t5**2-2._ki/3._ki*t6&
                  &*t1*t4*t3**2*t2**4-7._ki/18._ki*t6*t1*t2**5*t4**2*t3-15._ki/2._ki*t&
                  &6**3*t1**2*t3*t4*t5*t2+11._ki/12._ki*t6**3*t1*t3*t2**3*t5*t4-22._k&
                  &i/3._ki*t6**3*t1*t3*t4**4*t5-33._ki/4._ki*t6**3*t1*t2**2*t4*t3**2-&
                  &t6*t1**2*t2**5*t4*t5/6._ki-t6*t1**2*t2**4*t4*t5**2/6._ki-t6*t1**2&
                  &*t2**4*t5*t3/9._ki-7._ki/9._ki*t6*t1*t2**4*t3*t4**2*t5+7._ki/18._ki*&
                  &t6*t1*t2**5*t3*t4*t5
                !
                stemp6=stemp7+t6**2*t1**2*t3*t2**4/4._ki-t6**3*t3**2*t2**3*t4**2-4&
                  &._ki*t6**3*t3**2*t4**4*t2+2._ki/9._ki*t1*t2**5*t3**2*t5+2._ki/9._ki*&
                  &t6*t1*t2**5*t3**2-t6*t1**2*t2**5*t3/18._ki-t6**3*t3*t5*t1**3/2._k&
                  &i+25._ki/4._ki*t6**3*t1**2*t4*t3**2-5._ki/2._ki*t6**3*t1**2*t2*t3**&
                  &2-55._ki/3._ki*t6**3*t1*t4**3*t3**2+t6**3*t4*t3**2*t2**4/8._ki+11.&
                  &_ki/12._ki*t6**3*t1*t3**2*t2**3+3._ki*t6**3*t3**2*t2**2*t4**3+10._k&
                  &i*t6**3*t1**2*t3*t5*t4**2+5._ki/4._ki*t6**3*t1**2*t3*t2**2*t5+11.&
                  &_ki*t6**3*t1*t3*t4**3*t5*t2+22._ki*t6**3*t1*t2*t4**2*t3**2-11._ki/&
                  &2._ki*t6**3*t1*t3*t2**2*t5*t4**2+3._ki/8._ki*t6**2*t1*t3*t2**5*t4
                !
                stemp7=1._ki/t1/t2**8/t3
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
              !
            case(3)
              !
              select case(par4_glob)
              !
              case(3)
                !
                stemp5=195._ki/2._ki*t6**2*t2*t3**2*t4**2*t5-135._ki*t6**2*t3**2*t4*&
                  &t5*t2**2+39._ki/2._ki*t6**3*t1*t4**2*t2**3-27._ki/2._ki*t6**3*t1*t4&
                  &*t2**4-45._ki/2._ki*t6**2*t4*t3*t2**5-18._ki*t6**2*t4**3*t5**2*t2*&
                  &*3-12._ki*t6**3*t1*t5*t2**4+30._ki*t6**3*t1*t5**2*t2**3-15._ki*t6*&
                  &*2*t4*t5**3*t2**4-9._ki*t6**2*t1*t5**2*t2**4+15._ki/2._ki*t6**2*t1&
                  &*t5*t2**5-15._ki*t6**2*t1*t5**3*t2**3-30._ki*t6**2*t3*t5**2*t2**4&
                  &+6._ki*t6**2*t4**4*t5**2*t2**2+90._ki*t6**2*t4**3*t5**2*t3*t2-195&
                  &._ki*t6**2*t4**2*t5**2*t3*t2**2-160._ki*t6**3*t4*t3*t5*t2**3-t2**&
                  &9*t5/12._ki+t4*t2**9/6._ki-90._ki*t6**3*t1*t3*t5*t2**2-195._ki*t6**&
                  &3*t1*t4**2*t3*t5+7._ki/12._ki*t6**2*t2**8-5._ki/4._ki*t6**2*t5*t2**&
                  &7+95._ki*t6**3*t3**2*t4**3-t4**2*t5*t2**7/12._ki+t4*t2**8*t5/6._ki&
                  &-10._ki*t6**2*t2**2*t3**3+3._ki/2._ki*t6**2*t5**2*t2**6-45._ki*t6**&
                  &3*t3**2*t1*t4-375._ki/2._ki*t6**3*t3**2*t2*t4**2+120._ki*t6**3*t2*&
                  &*2*t3**2*t4+3._ki/2._ki*t6**2*t1**2*t5**2*t2**2+30._ki*t6**3*t1*t3&
                  &**2*t2-3._ki*t3**2*t2**2*t1*t6**2
                !
                stemp4=-12._ki*t6**2*t3*t5*t2**5-25._ki/4._ki*t6**3*t4**2*t2**5-25._k&
                  &i*t6**3*t2**3*t3**2+7._ki/3._ki*t6**2*t4**4*t2**4+195._ki*t6**3*t1&
                  &*t4**2*t5**2*t2-7._ki*t6**3*t4**4*t2**3+54._ki*t6**2*t4*t5*t3*t2*&
                  &*4-9._ki/4._ki*t6**3*t1**2*t2**3-5._ki*t6**3*t2**5*t3+2._ki*t6**3*t&
                  &4**5*t2**2+19._ki/2._ki*t6**3*t4**3*t2**4+t6*t4*t2**8+20._ki*t6**3&
                  &*t4**5*t5**2+39._ki/2._ki*t6**2*t4**2*t5**2*t2**4-30._ki*t6**2*t4*&
                  &*3*t5**3*t2**2+10._ki*t6**2*t4**4*t5**3*t2+65._ki/2._ki*t6**2*t4**&
                  &2*t5**3*t2**3-8._ki*t6**3*t4*t5*t2**5-35._ki/6._ki*t6**2*t1*t4**2*&
                  &t2**4+9._ki*t6**2*t3**2*t2**4-15._ki*t6**2*t1*t4**2*t5**2*t2**2+2&
                  &5._ki*t6**3*t4**2*t5*t2**4+140._ki*t6**3*t4**4*t3*t5+32._ki*t6**3*&
                  &t4*t3*t2**4+95._ki*t6**3*t4**3*t5**2*t2**2+76._ki*t6**3*t4**3*t3*&
                  &t2**2-90._ki*t6**3*t1*t4**3*t5**2+18._ki*t6**3*t1*t2**3*t3-7._ki*t&
                  &6**2*t4**3*t2**5+7._ki/12._ki*t1**2*t2**4*t6**2-54._ki*t6**3*t1*t2&
                  &**2*t3*t4+stemp5+t3*t2**8/6._ki+91._ki/12._ki*t6**2*t4**2*t2**6-20&
                  &._ki*t6**2*t1*t4*t5*t2**4-5._ki/4._ki*t6*t4**2*t2**7
                !
                stemp5=54._ki*t6**3*t1*t4*t5*t2**3+t6*t1*t2**7/2._ki-78._ki*t6**3*t1&
                  &*t4**2*t5*t2**2-65._ki/4._ki*t6**2*t4**2*t5*t2**5-5._ki*t6**2*t4**&
                  &4*t5*t2**3-15._ki*t6**2*t4**3*t3*t2**3+30._ki*t6**3*t1**2*t5**2*t&
                  &4-45._ki/2._ki*t6**3*t1**2*t5**2*t2-3._ki*t6**3*t1**2*t2*t3+25._ki*&
                  &t6**3*t3*t2**4*t5+15._ki*t6**2*t4**3*t5*t2**4-7._ki/2._ki*t6**2*t4&
                  &*t2**7-7._ki/2._ki*t6**2*t1*t2**6+t6*t4**3*t2**6/2._ki+5._ki*t6**2*&
                  &t3*t2**6+3._ki*t6**3*t1*t2**5-t6**3*t2**7/4._ki-5._ki/2._ki*t6**3*t&
                  &5**2*t2**5+40._ki*t6**2*t1*t4*t5**3*t2**2+t6**3*t5*t2**6+2._ki*t6&
                  &**3*t4*t2**6+5._ki/2._ki*t6**2*t5**3*t2**5-5._ki/4._ki*t6**2*t1**2*&
                  &t5*t2**3+24._ki*t6**2*t1*t4*t5**2*t2**3-t4*t3*t2**7/6._ki+375._ki*&
                  &t6**3*t4**2*t3*t5*t2**2+15._ki*t6**3*t1**2*t3*t5-28._ki*t6**3*t4*&
                  &*4*t3*t2-75._ki*t6**3*t4**2*t3*t2**3-70._ki*t6**3*t4**4*t5**2*t2+&
                  &15._ki/2._ki*t6**2*t2**6*t4*t5-9._ki*t6**2*t4*t5**2*t2**5+3._ki*t6*&
                  &*3*t1**2*t2**2*t4-60._ki*t6**2*t1*t2*t3*t4*t5**2-380._ki*t6**3*t4&
                  &**3*t3*t5*t2
                !
                stemp3=stemp5-t6*t2**9/4._ki+stemp4+5._ki/2._ki*t6**2*t1**2*t2*t5**3&
                  &-15._ki/2._ki*t6**2*t1*t2**4*t3-27._ki*t6**2*t2**3*t3**2*t4+15._ki*&
                  &t6**2*t2*t3**3*t4+45._ki*t6**2*t2**3*t3**2*t5+39._ki/2._ki*t6**2*t&
                  &3**2*t4**2*t2**2-15._ki*t6**2*t1*t2*t5*t3**2-25._ki*t6**2*t1*t4**&
                  &2*t5**3*t2-78._ki*t6**2*t4**2*t5*t3*t2**3+270._ki*t6**3*t1*t3*t5*&
                  &t2*t4-135._ki*t6**3*t1*t4*t5**2*t2**2+36._ki*t6**2*t4**3*t5*t3*t2&
                  &**2+135._ki*t6**2*t4*t5**2*t3*t2**3-12._ki*t6**3*t1**2*t5*t2*t4+3&
                  &6._ki*t6**3*t1*t4**3*t5*t2+39._ki*t6**3*t1*t2*t3*t4**2-24._ki*t6**&
                  &2*t1*t3*t4*t5*t2**2+25._ki/2._ki*t6**2*t1*t4**2*t5*t2**3+10._ki*t6&
                  &**2*t1*t2**3*t3*t4+45._ki*t6**2*t1*t5**2*t3*t2**2-t4**2*t2**8/12&
                  &._ki+28._ki/3._ki*t6**2*t1*t4*t2**5-t6*t4*t2**6*t1/2._ki+18._ki*t6**&
                  &2*t1*t2**3*t3*t5-t2**10/12._ki+20._ki*t6**3*t4*t5**2*t2**4+65._ki/&
                  &2._ki*t6**2*t4**2*t3*t2**4-9._ki*t6**3*t1*t4**3*t2**2-8._ki*t6**3*&
                  &t4**5*t5*t2-125._ki/2._ki*t6**3*t4**2*t5**2*t2**3-38._ki*t6**3*t4*&
                  &*3*t5*t2**3+28._ki*t6**3*t4**4*t5*t2**2+9._ki*t6**3*t1**2*t5*t2**&
                  &2
                !
                stemp4=1._ki/t2**10*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=-(-t3*t2**3+4._ki*t2**2*t1*t6+t2**2*t3*t4-3._ki*t3*t6*t2**2-&
                  &t3*t5*t2**2+9._ki*t6*t4*t2*t3+2._ki*t3**2*t2+4._ki*t5*t2*t1*t6+t3*&
                  &t5*t2*t4-4._ki*t4*t2*t1*t6-6._ki*t4**2*t3*t6-4._ki*t4*t1*t6*t5-2._k&
                  &i*t3*t1*t6)*(-t4+t2)/t2**4/t3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)&
                  &/12._ki
                !
                stemp9=-5._ki/4._ki*t6**3*t1**2*t2**6*t3*t5-41._ki/8._ki*t6**2*t1*t4*&
                  &*2*t3**2*t2**6+30._ki*t6**3*t1**4*t4**3*t5*t2+65._ki/2._ki*t6**3*t&
                  &1**4*t2*t3*t4**2-81._ki/2._ki*t6**2*t1**3*t3**2*t2**3*t4-t6**2*t1&
                  &*t2**8*t3**2/8._ki-55._ki/6._ki*t6**3*t1**3*t2**5*t3+11._ki/3._ki*t6&
                  &**3*t1**3*t4*t2**6+t6**3*t1**2*t2**7*t3/4._ki+9._ki/4._ki*t6**2*t1&
                  &**3*t5**2*t2**6-55._ki/12._ki*t6**3*t1**3*t5**2*t2**5+11._ki/6._ki*&
                  &t6**3*t1**3*t5*t2**6-45._ki*t6**2*t1**3*t4**3*t5**3*t2**2+16._ki*&
                  &t6**2*t1**2*t4*t3**2*t2**5-40._ki*t6**2*t1**2*t4*t3**3*t2**3-7._k&
                  &i/4._ki*t6**2*t1**4*t2**6+15._ki*t6**2*t1**3*t4**4*t5**3*t2-165._k&
                  &i*t6**3*t1**2*t4**3*t3**2*t2**2+55._ki*t6**2*t1*t4**3*t3**2*t5*t&
                  &2**4-45._ki*t6**2*t1**3*t3*t5**2*t2**4-45._ki/2._ki*t6**2*t1**3*t4&
                  &*t5**3*t2**4+195._ki/4._ki*t6**2*t1**3*t4**2*t5**3*t2**3+405._ki/2&
                  &._ki*t6**2*t1**3*t4*t5**2*t3*t2**3+45._ki/2._ki*t6**2*t1**3*t4**3*&
                  &t5*t2**4+54._ki*t6**2*t1**3*t4**3*t5*t3*t2**2-15._ki/2._ki*t6**2*t&
                  &1**3*t4**4*t5*t2**3+91._ki/8._ki*t6**2*t1**3*t4**2*t2**6
                !
                stemp8=7._ki/2._ki*t6**2*t1**3*t4**4*t2**4-21._ki/2._ki*t6**2*t1**3*t&
                  &4**3*t2**5+27._ki/2._ki*t6**2*t1**3*t3**2*t2**4+t6*t1**3*t4**3*t2&
                  &**6/2._ki-5._ki/4._ki*t6*t1**3*t4**2*t2**7+9._ki*t6**2*t1**4*t2**3*&
                  &t3*t5-18._ki*t6**2*t1**3*t3*t5*t2**5-125._ki/12._ki*t6**2*t1*t2**5&
                  &*t3**3*t4+117._ki/4._ki*t6**2*t1**3*t3**2*t2**2*t4**2-40._ki/3._ki*&
                  &t6**2*t1*t4**5*t3**3*t2-13._ki*t6**2*t1*t4**4*t3**2*t2**4+117._ki&
                  &/4._ki*t6**2*t1**3*t4**2*t5**2*t2**4-205._ki/8._ki*t6**2*t1*t4**2*&
                  &t3**2*t2**5*t5+80._ki*t6**2*t1**2*t4*t3**2*t5*t2**4-95._ki/3._ki*t&
                  &6**2*t1**2*t3**3*t4**3*t2+275._ki/6._ki*t6**3*t1**3*t4**2*t5*t2**&
                  &4+110._ki/3._ki*t6**3*t1**3*t4*t5**2*t2**4+125._ki/2._ki*t6**2*t1**&
                  &2*t3**3*t4**2*t2**2-1375._ki/4._ki*t6**3*t1**3*t3**2*t4**2*t2+41.&
                  &_ki/4._ki*t6**3*t1**2*t4**2*t3*t2**5+5._ki*t6**2*t1**4*t2**3*t3*t4&
                  &+45._ki/4._ki*t6**2*t1**3*t2**6*t4*t5-125._ki/4._ki*t6**3*t1**2*t4*&
                  &t3**2*t2**4-125._ki/12._ki*t6**2*t1**2*t4**2*t3*t2**6-35._ki/3._ki*&
                  &t6**2*t1**2*t4**4*t3*t2**4-t1**3*t4*t3*t2**7/18._ki-117._ki*t6**2&
                  &*t1**3*t4**2*t5*t3*t2**3+5._ki/2._ki*t6**3*t1**4*t2**5+stemp9
                !
                stemp9=-t1**3*t2**9*t5/9._ki-t1**3*t4**2*t2**8/9._ki+325._ki/2._ki*t6&
                  &**3*t1**4*t4**2*t5**2*t2-75._ki*t6**3*t1**4*t3*t5*t2**2+t1**3*t3&
                  &*t2**8/18._ki-8._ki*t6**2*t1**2*t2**6*t3*t4*t5+38._ki*t6**2*t1**2*&
                  &t4**3*t3**2*t2**3+40._ki*t6**2*t1*t4**5*t3**2*t5*t2**2-45._ki*t6*&
                  &*3*t1**4*t2**2*t3*t4-22._ki*t6**3*t1**2*t4**3*t3*t2**4-16._ki*t6*&
                  &*3*t1**2*t4**5*t3*t2**2+95._ki/6._ki*t6**2*t1**2*t4**3*t3*t2**5+2&
                  &5._ki*t6**2*t1**2*t4**2*t5*t3*t2**5+2._ki/9._ki*t1**3*t4*t2**9-2._k&
                  &i*t6**2*t1*t4**6*t3**2*t2**2+8._ki*t6**2*t1*t4**5*t3**2*t2**3-20&
                  &._ki*t6**2*t1**2*t3*t5**2*t2**5*t4+25._ki/2._ki*t6**3*t1**2*t4*t5*&
                  &t2**5*t3-205._ki/4._ki*t6**3*t1**2*t4**2*t5*t2**4*t3-27._ki/2._ki*t&
                  &6**2*t1**3*t4*t5**2*t2**5-12._ki*t6**2*t1**4*t2**2*t3*t4*t5-30._k&
                  &i*t6**2*t1**4*t5**2*t3*t2*t4-1375._ki/12._ki*t6**3*t1**3*t4**2*t5&
                  &**2*t2**3-38._ki*t6**2*t1**2*t4**3*t5*t3*t2**4-20._ki*t6**2*t1**2&
                  &*t4**5*t3*t5**2*t2+t6**2*t2**8*t3**3/12._ki-65._ki*t6**3*t1**4*t4&
                  &**2*t5*t2**2+stemp8
                !
                stemp7=205._ki/2._ki*t6**3*t1**2*t4**2*t3**2*t2**3-20._ki*t6**3*t1**&
                  &2*t4**6*t3*t5-10._ki*t6**2*t1*t4**6*t3**2*t5*t2-15._ki/2._ki*t6**2&
                  &*t1**4*t5**3*t2**3+15._ki/4._ki*t6**2*t1**4*t5*t2**5-21._ki/4._ki*t&
                  &6**2*t1**3*t4*t2**7+t6*t1**3*t4*t2**8+61._ki/12._ki*t6**2*t2**6*t&
                  &3**3*t4**2+5._ki/4._ki*t6**2*t1*t2**6*t3**3-85._ki/6._ki*t6**2*t2**&
                  &5*t3**3*t4**3+15._ki/4._ki*t6**2*t1**3*t5**3*t2**5-15._ki/8._ki*t6*&
                  &*2*t1**3*t5*t2**7-325._ki/2._ki*t6**3*t1**4*t4**2*t3*t5+7._ki/8._ki&
                  &*t6**2*t1**3*t2**8+11._ki*t6**2*t1*t3**2*t4**3*t2**5-275._ki/2._ki&
                  &*t6**3*t1**3*t4**2*t3*t2**3-5._ki/2._ki*t6**3*t1**2*t4*t3*t2**6+1&
                  &30._ki/3._ki*t6**2*t1*t4**4*t3**3*t2**2-44._ki/3._ki*t6**3*t1**3*t4&
                  &*t5*t2**5-3._ki/2._ki*t6**3*t1*t4*t3**2*t2**6+stemp9-3._ki/4._ki*t6&
                  &**3*t1**5*t2**3+4._ki*t6**3*t1**2*t4**6*t3*t2-t6**2*t2**7*t3**3*&
                  &t4-10._ki*t6**3*t1**4*t5*t2**4+25._ki/3._ki*t6**2*t1**2*t3**3*t2**&
                  &4-t1**3*t4**2*t5*t2**7/9._ki+225._ki*t6**3*t1**4*t3*t5*t2*t4-65._k&
                  &i*t6**2*t1*t4**4*t3**2*t5*t2**3
                !
                stemp9=stemp7-585._ki/2._ki*t6**2*t1**3*t4**2*t5**2*t3*t2**2-55._ki*&
                  &t6**2*t1*t2**3*t3**3*t4**3-11._ki/24._ki*t6**3*t1**3*t2**7-85._ki/&
                  &4._ki*t6**3*t1*t4**3*t3**2*t2**4+12._ki*t6**2*t1**4*t4*t5**2*t2**&
                  &3-15._ki/2._ki*t6**2*t1**4*t3**2*t5*t2+45._ki/2._ki*t6**2*t1**4*t5*&
                  &*2*t3*t2**2+81._ki*t6**2*t1**3*t4*t5*t3*t2**4-44._ki/3._ki*t6**3*t&
                  &1**3*t4**5*t5*t2-t6*t1**3*t2**9/4._ki+11._ki/3._ki*t6**3*t1**3*t4*&
                  &*5*t2**2+209._ki/12._ki*t6**3*t1**3*t4**3*t2**4-275._ki/6._ki*t6**3&
                  &*t1**3*t3**2*t2**3+110._ki/3._ki*t6**3*t1**3*t4**5*t5**2+1045._ki/&
                  &6._ki*t6**3*t1**3*t3**2*t4**3-275._ki/24._ki*t6**3*t1**3*t4**2*t2*&
                  &*5-77._ki/6._ki*t6**3*t1**3*t4**4*t2**3-75._ki/2._ki*t6**2*t1**2*t4&
                  &**2*t3**2*t2**4-385._ki/3._ki*t6**3*t1**3*t4**4*t5**2*t2+1045._ki/&
                  &6._ki*t6**3*t1**3*t4**3*t5**2*t2**2+t6**3*t1*t3**2*t2**7/8._ki+t6&
                  &**2*t1**2*t5*t2**7*t3+28._ki*t6**2*t1**2*t4**4*t5*t3*t2**3+26._ki&
                  &*t6**3*t1**2*t4**4*t3*t2**3+205._ki/6._ki*t6**2*t1*t2**4*t3**3*t4&
                  &**2+190._ki*t6**2*t1**2*t3**2*t4**3*t5*t2**2
                !
                stemp8=stemp9+5._ki/2._ki*t6**2*t1**2*t3*t5**2*t2**6-t1**3*t2**10/9&
                  &._ki+5._ki/4._ki*t6**2*t1*t4*t3**2*t2**7-15._ki*t6**2*t1**3*t3**3*t&
                  &2**2-35._ki/12._ki*t6**2*t1**4*t4**2*t2**4-15._ki/4._ki*t6**2*t1**4&
                  &*t2**4*t3-3._ki/2._ki*t6**2*t1**4*t3**2*t2**2+15._ki/4._ki*t6**3*t1&
                  &**2*t3**2*t2**5-5._ki/2._ki*t6**2*t1**2*t3**2*t2**6-405._ki/2._ki*t&
                  &6**2*t1**3*t3**2*t5*t2**2*t4+135._ki*t6**2*t1**3*t4**3*t5**2*t3*&
                  &t2+585._ki/4._ki*t6**2*t1**3*t3**2*t5*t2*t4**2-34._ki*t6**3*t1*t4*&
                  &*5*t3**2*t2**2+25._ki/4._ki*t6**2*t1*t4*t3**2*t2**6*t5-40._ki*t6**&
                  &3*t1**2*t4**5*t3**2-4._ki*t6**3*t1*t4**7*t3**2+5._ki*t6**3*t1**5*&
                  &t3*t5+t6**3*t1**5*t2**2*t4+10._ki*t6**3*t1**5*t5**2*t4-15._ki/2._k&
                  &i*t6**3*t1**5*t5**2*t2+3._ki*t6**3*t1**5*t5*t2**2-t6**3*t1**5*t2&
                  &*t3-45._ki/2._ki*t6**2*t1**3*t4**3*t3*t2**3+61._ki/8._ki*t6**3*t1*t&
                  &4**2*t3**2*t2**5+1375._ki/2._ki*t6**3*t1**3*t4**2*t3*t5*t2**2+80.&
                  &_ki*t6**3*t1**2*t4**5*t3*t5*t2-4._ki*t6**3*t1**5*t5*t2*t4+18._ki*t&
                  &6**3*t1*t4**6*t3**2*t2
                !
                stemp9=stemp8-135._ki/4._ki*t6**2*t1**3*t4*t3*t2**5+275._ki/6._ki*t6*&
                  &*3*t1**3*t3*t2**4*t5+35._ki*t6**3*t1*t4**4*t3**2*t2**3+10._ki/3._k&
                  &i*t6**2*t1**2*t4*t3*t2**7-14._ki*t6**2*t1**2*t4**4*t3**2*t2**2-5&
                  &._ki/8._ki*t6**2*t1*t2**7*t3**2*t5-130._ki*t6**3*t1**2*t4**4*t3*t5&
                  &*t2**2-225._ki/2._ki*t6**3*t1**4*t4*t5**2*t2**2+154._ki/3._ki*t6**3&
                  &*t1**3*t4**4*t5*t2**2-209._ki/3._ki*t6**3*t1**3*t4**3*t5*t2**3-70&
                  &._ki*t6**2*t1**2*t3**2*t4**4*t5*t2-15._ki/2._ki*t6**2*t1**4*t4**2*&
                  &t5**2*t2**2+25._ki*t6**3*t1**4*t5**2*t2**3-5._ki/12._ki*t6**2*t1**&
                  &2*t3*t2**8+15._ki/2._ki*t6**2*t1**3*t3*t2**6-9._ki/2._ki*t6**2*t1**&
                  &4*t5**2*t2**4-25._ki/2._ki*t6**2*t1**2*t2**5*t3**2*t5+418._ki/3._ki&
                  &*t6**3*t1**3*t4**3*t3*t2**2+110._ki*t6**3*t1**2*t4**3*t3*t5*t2**&
                  &3-375._ki/2._ki*t6**2*t1**2*t3**2*t4**2*t5*t2**3+45._ki/2._ki*t6**2&
                  &*t1**3*t3**3*t2*t4+770._ki/3._ki*t6**3*t1**3*t4**4*t3*t5-10._ki*t6&
                  &**2*t1**4*t4*t5*t2**4-25._ki/2._ki*t6**2*t1**4*t4**2*t5**3*t2+10.&
                  &_ki/3._ki*t6**2*t1**2*t4**5*t3*t2**3+220._ki*t6**3*t1**3*t4*t3**2*&
                  &t2**2-154._ki/3._ki*t6**3*t1**3*t4**4*t3*t2
                !
                stemp6=stemp9+130._ki*t6**3*t1**2*t4**4*t3**2*t2+70._ki*t6**2*t1**2&
                  &*t4**4*t3*t5**2*t2**2-15._ki/2._ki*t6**3*t1**4*t4**3*t2**2+65._ki/&
                  &4._ki*t6**3*t1**4*t4**2*t2**3-45._ki/4._ki*t6**3*t1**4*t4*t2**4-75&
                  &._ki/2._ki*t6**3*t1**4*t3**2*t4-75._ki*t6**3*t1**4*t4**3*t5**2+15.&
                  &_ki*t6**3*t1**4*t2**3*t3+25._ki*t6**3*t1**4*t3**2*t2+176._ki/3._ki*&
                  &t6**3*t1**3*t4*t3*t2**4-27._ki*t6**2*t1**3*t4**3*t5**2*t2**3+25.&
                  &_ki/4._ki*t6**2*t1**4*t4**2*t5*t2**3+135._ki/2._ki*t6**2*t1**3*t2**&
                  &3*t3**2*t5+195._ki/4._ki*t6**2*t1**3*t4**2*t3*t2**4-195._ki/8._ki*t&
                  &6**2*t1**3*t4**2*t5*t2**5+2._ki/9._ki*t1**3*t4*t2**8*t5+70._ki/3._k&
                  &i*t6**2*t4**4*t3**3*t2**4-68._ki/3._ki*t6**2*t4**5*t3**3*t2**3+12&
                  &._ki*t6**2*t4**6*t3**3*t2**2-8._ki/3._ki*t6**2*t4**7*t3**3*t2+14._k&
                  &i/3._ki*t6**2*t1**4*t4*t2**5-880._ki/3._ki*t6**3*t1**3*t4*t3*t5*t2&
                  &**3-95._ki*t6**2*t1**2*t4**3*t3*t5**2*t2**3-2090._ki/3._ki*t6**3*t&
                  &1**3*t4**3*t3*t5*t2+9._ki*t6**2*t1**3*t4**4*t5**2*t2**2+20._ki*t6&
                  &**2*t1**4*t4*t5**3*t2**2-8._ki*t6**2*t1**2*t4**5*t5*t3*t2**2+45.&
                  &_ki*t6**3*t1**4*t4*t5*t2**3+125._ki/2._ki*t6**2*t1**2*t4**2*t3*t5*&
                  &*2*t2**4
                !
                stemp7=1._ki/t1**3/t2**10
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(4)
                !
                stemp2=(-6._ki*t6*t2**5*t4*t5+4._ki*t6*t2**4*t4**2*t5+6._ki*t3*t4*t6&
                  &*t2**4+18._ki*t6**2*t1*t2**3*t4-144._ki*t6**3*t1*t5*t4**2-t3*t2**&
                  &6+t5**2*t2**6+3._ki*t6**2*t2**6+t6*t2**7+t5*t2**7-t2**6*t4*t5-t4&
                  &*t5**2*t2**5-12._ki*t6**2*t1*t2**4+18._ki*t6**3*t1**2*t5+6._ki*t6*&
                  &*3*t5*t2**4+2._ki*t6*t2**6*t5-t6*t1*t2**5+48._ki*t6**3*t4**4*t5+2&
                  &4._ki*t6**2*t2**4*t4**2-12._ki*t6**2*t2**3*t4**3-3._ki*t6*t2**6*t4&
                  &+2._ki*t6*t2**5*t4**2-15._ki*t6**2*t2**5*t4+120._ki*t6**3*t3*t4**3&
                  &-90._ki*t6**3*t4*t1*t3-216._ki*t6**3*t3*t2*t4**2+126._ki*t6**3*t3*&
                  &t2**2*t4-2._ki*t6*t1*t2**4*t5-54._ki*t1*t6**3*t2**2*t5-42._ki*t6**&
                  &3*t2**3*t5*t4+108._ki*t6**3*t2**2*t5*t4**2-120._ki*t6**3*t4**3*t5&
                  &*t2+54._ki*t6**3*t3*t2*t1-24._ki*t6**3*t3*t2**3-4._ki*t6*t2**5*t3-&
                  &2._ki*t3*t5*t2**5+180._ki*t1*t6**3*t4*t5*t2)/t2**8*z_log(t1*t6/t2&
                  &**2,1._ki)/12._ki
                !
                stemp6=t4*t5*t2*t6**2*t1**2/2._ki+t3**2*t4*t5**2*t2**2/12._ki+t3*t2&
                  &**2*t5**2*t1*t6/3._ki+t3*t5*t2**3*t1*t6/3._ki+t6*t3**3*t2**2/3._ki&
                  &+2._ki/3._ki*t3*t2**2*t6**2*t5*t1+t2**3*t3**3/12._ki-t3**2*t2**3*t&
                  &6**2/4._ki-t2**2*t4*t1*t6**2*t3+2._ki/3._ki*t2*t4**2*t1*t6**2*t3+t&
                  &4*t5*t3**2*t6*t2**2/2._ki+t3**2*t2**3*t4*t5/12._ki+t3*t2**3*t6**2&
                  &*t1/3._ki-t2*t5**2*t1**2*t6**2/2._ki-t2**2*t1**2*t6**2*t5/2._ki-t3&
                  &**2*t5**2*t2**3/12._ki-t3**2*t2**4*t6/12._ki+t4**3*t3**2*t6**2+t2&
                  &**2*t5*t3**3/6._ki
                !
                stemp5=stemp6-t3**2*t2**4*t5/12._ki+t1**2*t6**2*t3*t2/6._ki-t1*t6*t&
                  &3**2*t2**2/4._ki-t5*t3**2*t6*t2**3/6._ki+t1*t6**2*t3**2*t4/2._ki-t&
                  &2*t6*t3**3*t4/2._ki-t2*t1*t6**2*t3**2/3._ki+t5*t1**2*t6**2*t3/3._k&
                  &i+t4*t5**2*t6**2*t1**2/2._ki-2._ki*t2*t4**2*t3**2*t6**2+5._ki/4._ki&
                  &*t2**2*t4*t3**2*t6**2+t2**3*t4*t6*t3**2/4._ki-t2**2*t4**2*t6*t3*&
                  &*2/6._ki-t4*t3*t5*t2**2*t1*t6/3._ki-t2*t5*t1*t6*t3**2/2._ki-2._ki*t&
                  &3*t4*t5*t2*t6**2*t1+4._ki/3._ki*t3*t4**2*t5*t6**2*t1-t3**2*t4**2*&
                  &t5*t2*t6/3._ki-t3*t4*t5**2*t2*t1*t6/3._ki
                !
                stemp6=1._ki/t3**2/t2**5*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)
                !
                stemp4=stemp5*stemp6
                !
                stemp7=7._ki/36._ki*t6*t1*t3*t2**7-25._ki/4._ki*t6**3*t1**2*t4*t3**2+&
                  &15._ki/4._ki*t6**3*t1**2*t2*t3**2+55._ki/3._ki*t6**3*t1*t4**3*t3**2&
                  &-9._ki/8._ki*t6**3*t4*t3**2*t2**4-11._ki/3._ki*t6**3*t1*t3**2*t2**3&
                  &-7._ki*t6**3*t3**2*t2**2*t4**3-t6**2*t1**2*t3*t2**4/2._ki+4._ki*t6&
                  &**3*t3**2*t2**3*t4**2+6._ki*t6**3*t3**2*t4**4*t2+t6**3*t3*t5*t1*&
                  &*3/2._ki-4._ki/9._ki*t6*t1*t2**5*t3**2+t6*t1**2*t2**5*t3/18._ki+3._k&
                  &i/8._ki*t6**2*t2**6*t1*t3-t6*t1**2*t5*t2**6/6._ki-2._ki/9._ki*t1*t2&
                  &**5*t3**2*t5+7._ki/36._ki*t1*t2**6*t3*t5**2+7._ki/36._ki*t1*t3*t2**&
                  &7*t5-t6*t1**2*t2**5*t5**2/6._ki+22._ki/3._ki*t6**3*t1*t3*t4**4*t5+&
                  &77._ki/4._ki*t6**3*t1*t2**2*t4*t3**2-10._ki*t6**3*t1**2*t3*t5*t4**&
                  &2-15._ki/4._ki*t6**3*t1**2*t3*t2**2*t5
                !
                stemp6=-55._ki/3._ki*t6**3*t1*t3*t4**3*t5*t2+stemp7+t6*t1**2*t2**4*&
                  &t4*t5**2/6._ki+t6*t1**2*t2**4*t5*t3/9._ki+11._ki/12._ki*t6**3*t1*t3&
                  &*t5*t2**4-7._ki/12._ki*t6*t1*t2**6*t4*t3+2._ki/3._ki*t6*t1*t4*t3**2&
                  &*t2**4+7._ki/18._ki*t6*t1*t2**5*t4**2*t3+t6**3*t3**2*t2**5/8._ki-2&
                  &._ki*t6**3*t3**2*t4**5-t1*t3**2*t2**6/9._ki-33._ki*t6**3*t1*t2*t4*&
                  &*2*t3**2+33._ki/2._ki*t6**3*t1*t3*t2**2*t5*t4**2+7._ki/18._ki*t6*t1&
                  &*t2**6*t3*t5-7._ki/6._ki*t6*t1*t2**5*t3*t4*t5+7._ki/9._ki*t6*t1*t2*&
                  &*4*t3*t4**2*t5+3._ki*t6**2*t1*t3*t2**4*t4**2+3._ki/4._ki*t6**2*t1*&
                  &*2*t3*t2**3*t4-15._ki/8._ki*t6**2*t1*t3*t2**5*t4+t6*t1**2*t2**5*t&
                  &4*t5/6._ki-7._ki/36._ki*t1*t2**6*t3*t4*t5-7._ki/36._ki*t1*t2**5*t3*t&
                  &4*t5**2-3._ki/2._ki*t6**2*t1*t3*t2**3*t4**3+25._ki/2._ki*t6**3*t1**&
                  &2*t3*t4*t5*t2-77._ki/12._ki*t6**3*t1*t3*t2**3*t5*t4
                !
                stemp7=1._ki/t1/t2**8/t3
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
              !
            case(4)
              !
              select case(par4_glob)
              !
              case(4)
                !
                stemp2=-(-18._ki*t6**3*t1*t2+36._ki*t6**3*t1*t4-24._ki*t6**3*t4**3+3&
                  &._ki*t6**3*t2**3-18._ki*t6**3*t2**2*t4+36._ki*t6**3*t4**2*t2-2._ki*&
                  &t6**2*t1*t2**2-6._ki*t6**2*t1*t5*t2+t6**2*t2**4+3._ki*t6**2*t5*t2&
                  &**3-6._ki*t6**2*t2**2*t3-4._ki*t6**2*t2**3*t4+12._ki*t6**2*t5*t2*t&
                  &4**2+4._ki*t6**2*t2**2*t4**2-12._ki*t6**2*t5*t2**2*t4+12._ki*t6**2&
                  &*t2*t3*t4-6._ki*t6*t4*t5**2*t2**2+3._ki*t6*t2**3*t5**2+2._ki*t6*t2&
                  &**4*t5-4._ki*t6*t2**3*t4*t5-2._ki*t6*t3*t2**3-6._ki*t6*t2**2*t3*t5&
                  &+3._ki*t5**3*t2**3+3._ki*t5**2*t2**4)/t2**6*z_log(t1*t6/t2**2,1._k&
                  &i)/12._ki
                !
                stemp6=-2._ki/3._ki*t3**2*t2**3*t6**2*t5*t1-t3**2*t5**2*t2**3*t1*t6&
                  &+3._ki/2._ki*t3*t5**2*t2**2*t1**2*t6**2-2._ki*t6**3*t1**2*t3*t4*t5&
                  &*t2+3._ki/2._ki*t2*t3*t1**2*t6**3*t5**2+3._ki/2._ki*t6**2*t5*t2*t1*&
                  &t3**3-t3**3*t2**3*t6*t4*t5/3._ki-t2**2*t3**2*t1*t6**2*t5**2-4._ki&
                  &*t3**2*t1*t6**3*t5*t4**2-t2**2*t3**3*t6*t4*t5**2/2._ki+t2*t3**3*&
                  &t6**2*t5*t4**2-t2**2*t3**3*t6**2*t5*t4-t2**2*t3**2*t1*t6**3*t5-&
                  &3._ki*t3*t1**2*t6**3*t4*t5**2+3._ki/2._ki*t2*t3*t1**2*t6**2*t5**3-&
                  &t2**2*t3**2*t1*t6*t5**3-2._ki*t3**3*t6**3*t4**3-t1**3*t6**3*t5**&
                  &3+t2**3*t3**3*t6**3/4._ki-t3**4*t2**2*t6**2/2._ki-t6*t2**3*t3**4/&
                  &6._ki+t3**3*t2**4*t6**2/12._ki+t2**3*t3**3*t5**3/4._ki
                !
                stemp5=stemp6+t3**3*t5**2*t2**4/4._ki+4._ki/3._ki*t3**2*t2**2*t6**2*&
                  &t4*t5*t1-t6**3*t1*t3**2*t2**3/3._ki-t6**3*t1**2*t2*t3**2/3._ki-3.&
                  &_ki/2._ki*t2**2*t3**3*t6**3*t4+3._ki*t2*t3**3*t6**3*t4**2-t6*t5*t2&
                  &**2*t3**4/2._ki+t3**3*t2**2*t6**2*t4**2/3._ki+t3**3*t2**4*t6*t5/6&
                  &._ki-t3**3*t2**3*t6**2*t4/3._ki-t5**2*t2*t1**3*t6**3+t6**2*t2**2*&
                  &t1*t3**3/2._ki+t3**4*t2*t6**2*t4+t3**3*t2*t6**3*t1/2._ki+t2**3*t3&
                  &**3*t6**2*t5/4._ki+t2**3*t3**3*t6*t5**2/4._ki-t6**3*t5*t1**2*t3**&
                  &2-t3**3*t6**3*t1*t4+4._ki*t2*t3**2*t1*t6**3*t5*t4+2._ki*t2*t3**2*&
                  &t1*t6**2*t4*t5**2+4._ki/3._ki*t6**3*t1*t2**2*t4*t3**2+t6**3*t1**2&
                  &*t3*t2**2*t5-4._ki/3._ki*t6**3*t1*t2*t4**2*t3**2
                !
                stemp6=1._ki/t2**6/t3**3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)
                !
                stemp4=stemp5*stemp6
                !
                stemp5=-(66._ki*t3**2*t2**3*t5**3+48._ki*t6**2*t1*t3*t4*t5*t2**2+72&
                  &._ki*t6**2*t1*t2*t3*t4*t5**2-28._ki*t3**3*t2**3*t6+33._ki*t6**3*t2&
                  &**3*t3**2+66._ki*t3**2*t2**4*t5**2-264._ki*t6**3*t3**2*t4**3-66._k&
                  &i*t6**2*t2**2*t3**3+17._ki*t6**2*t3**2*t2**4-120._ki*t6*t3**2*t2*&
                  &*2*t4*t5**2-68._ki*t6**2*t2**3*t3**2*t4+24._ki*t6**2*t1**2*t2*t5*&
                  &*3+180._ki*t6**3*t3**2*t1*t4-90._ki*t6**3*t1*t3**2*t2+68._ki*t6**2&
                  &*t3**2*t4**2*t2**2+51._ki*t6**2*t2**3*t3**2*t5+132._ki*t6**2*t2*t&
                  &3**3*t4+24._ki*t6**2*t1**2*t5**2*t2**2-198._ki*t6**3*t2**2*t3**2*&
                  &t4+396._ki*t6**3*t3**2*t2*t4**2+40._ki*t6*t2**4*t3**2*t5+60._ki*t6&
                  &*t3**2*t2**3*t5**2-84._ki*t6*t3**3*t2**2*t5+2._ki*t3**2*t2**2*t1*&
                  &t6**2-24._ki*t6**2*t1*t2**3*t3*t5+6._ki*t6**2*t1*t2*t5*t3**2-36._k&
                  &i*t6**2*t1*t5**2*t3*t2**2+204._ki*t6**2*t2*t3**2*t4**2*t5-204._ki&
                  &*t6**2*t3**2*t4*t5*t2**2-48._ki*t6*t1*t3*t2**3*t5**2-48._ki*t6*t1&
                  &*t3*t2**2*t5**3-80._ki*t6*t3**2*t2**3*t4*t5)/t3**2/t2**6/72._ki
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
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
                stemp6=30._ki*t3*t5**2*t2**2*t1**2*t6**2-60._ki*t3**2*t5**2*t2**3*t&
                  &1*t6+60._ki*t3**2*t2**3*t6**2*t5*t1+480._ki*t6**2*t4**4*t3*t5**2*&
                  &t2**2-80._ki*t6**2*t4**4*t5*t2**3*t3+90._ki*t6**2*t4**2*t3*t5**2*&
                  &t2**4+9._ki*t6**2*t1*t3*t2**4*t4**2-75._ki*t6*t1*t4*t5**4*t2**2*t&
                  &3-3._ki*t6**2*t1*t3*t2**5*t4-30._ki*t6**2*t1*t4*t5*t2**4*t3-540._k&
                  &i*t6**2*t1*t4**2*t3*t5**2*t2**2-720._ki*t6**2*t1*t4**2*t5**3*t2*&
                  &t3+240._ki*t6**2*t1*t4*t5**3*t3*t2**2+180._ki*t6**2*t4*t3**2*t5**&
                  &2*t2**3-120._ki*t6*t1*t4*t5**3*t3*t2**3-120._ki*t6*t4**2*t3*t5**3&
                  &*t2**4-75._ki*t6*t4**2*t5**4*t3*t2**3-36._ki*t6*t4*t5*t2**5*t3**2&
                  &+108._ki*t6*t4**2*t3**2*t5*t2**4-36._ki*t6*t4**2*t5**2*t2**5*t3+9&
                  &._ki*t6**3*t1*t4**2*t5*t2**4-45._ki*t6**3*t4**2*t3*t5**2*t2**3
                !
                stemp5=stemp6-20._ki*t3**4*t2**2*t6**2-5._ki*t2**3*t3**3*t6**3-5._ki&
                  &*t1**3*t6**3*t5**3-4._ki*t3**3*t6*t2**5+40._ki*t6**3*t4**6*t5**3-&
                  &t6**2*t4**4*t2**6-50._ki*t2**2*t3**2*t1*t6*t5**3+40._ki*t2*t3*t1*&
                  &*2*t6**2*t5**3+225._ki*t3*t1**2*t6**3*t4*t5**2-45._ki*t2**2*t3**2&
                  &*t1*t6**3*t5+450._ki*t2*t3**2*t1*t6**3*t5*t4-600._ki*t2*t3**2*t1*&
                  &t6**2*t4*t5**2+640._ki*t6**2*t4**4*t3*t5**3*t2-15._ki*t6**2*t4**2&
                  &*t5*t2**5*t3+70._ki*t6**2*t4**3*t5*t2**4*t3-36._ki*t6*t1*t4*t3*t2&
                  &**4*t5**2-150._ki*t6*t4*t5**3*t3**2*t2**3-300._ki*t3**2*t2**2*t6*&
                  &*2*t4*t5*t1+200._ki*t3**3*t2**3*t6*t4*t5-15._ki*t3*t2**2*t6**3*t4&
                  &*t1**2+25._ki*t3**2*t2**3*t1*t6**2*t4-5._ki*t3*t2**3*t1**2*t6**2*&
                  &t5+6._ki*t6*t1*t2**5*t3*t4*t5
                !
                stemp6=-3._ki*t6*t1*t2**6*t4*t3-80._ki*t6**2*t5*t2*t1*t3**3-1125._ki&
                  &*t6**3*t4**3*t3**2*t5*t2+405._ki*t6**3*t4**2*t3**2*t5*t2**2-400.&
                  &_ki*t2**2*t3**3*t6**2*t5*t4+60._ki*t6**2*t1*t4**2*t5**4*t2**2+7._k&
                  &i*t6**2*t1*t4**3*t5*t2**4-140._ki*t6**2*t1*t4**3*t5**4*t2-3._ki*t&
                  &6**2*t1*t4**2*t5*t2**5-15._ki*t6**2*t1*t4**2*t5**2*t2**4+60._ki*t&
                  &6**2*t1*t4**2*t5**3*t2**3-140._ki*t6**2*t1*t4**3*t5**3*t2**2-3._k&
                  &i/2._ki*t6**2*t1**2*t4*t5*t2**4-t4**3*t2**9/4._ki+stemp5+30._ki*t6&
                  &**2*t1**2*t4*t5**4*t2-15._ki/2._ki*t6**2*t1**2*t4*t5**2*t2**3+30.&
                  &_ki*t6**2*t1**2*t4*t5**3*t2**2+35._ki*t6**2*t1*t4**3*t5**2*t2**3-&
                  &12._ki*t3**2*t2**4*t1*t6*t5-45._ki*t2*t3*t1**2*t6**3*t5**2+120._ki&
                  &*t2**2*t3**2*t1*t6**2*t5**2-855._ki*t3**2*t1*t6**3*t5*t4**2
                !
                stemp4=stemp6-14._ki*t6*t4**3*t5*t3*t2**5+6._ki*t6*t4**2*t5*t3*t2**&
                  &6-180._ki*t6*t4*t3**2*t5**2*t2**4+450._ki*t6*t4**2*t5**3*t3**2*t2&
                  &**2-12._ki*t6*t1*t4**2*t5**3*t2**4+t6**2*t4**3*t2**7/4._ki+315._ki&
                  &*t3**3*t6**3*t4**3+t6**2*t4**5*t2**5+10._ki*t6*t2**3*t3**4+90._ki&
                  &*t6**2*t1*t4**2*t5*t2**3*t3+180._ki*t6**2*t1*t4*t3*t5**2*t2**3+1&
                  &500._ki*t6**2*t4**3*t3**2*t5**2*t2+90._ki*t6**2*t4*t3**2*t5*t2**4&
                  &-540._ki*t6**2*t4**2*t3**2*t5*t2**3-560._ki*t6**2*t4**3*t3*t5**3*&
                  &t2**2+10._ki*t3**3*t2**4*t6**2-1080._ki*t6**2*t4**2*t3**2*t5**2*t&
                  &2**2+750._ki*t6**2*t4**3*t3**2*t5*t2**2+760._ki*t2*t3**3*t6**2*t5&
                  &*t4**2+250._ki*t2**2*t3**3*t6*t4*t5**2+9._ki*t6**3*t1**2*t4*t5*t2&
                  &**3-27._ki*t6**3*t1**2*t4**2*t5*t2**2-45._ki*t6**3*t1**2*t4*t5**3&
                  &*t2+75._ki*t6**3*t1*t2**2*t3*t4**3
                !
                stemp6=stemp4-45._ki*t6**3*t3**2*t2**3*t4*t5-135._ki*t6**3*t1*t4*t3&
                  &*t5**2*t2**2+810._ki*t6**3*t1*t4**2*t3*t5**2*t2+315._ki*t6**3*t4*&
                  &*3*t3*t5**2*t2**2-720._ki*t6**3*t4**4*t3*t5**2*t2+t6*t1*t2**5*t3&
                  &**2-100._ki*t3**3*t2**3*t6**2*t4-40._ki*t3**3*t2**4*t6*t5+190._ki*&
                  &t3**3*t2**2*t6**2*t4**2+3._ki*t3*t2**3*t6**3*t1**2-5._ki*t3**2*t2&
                  &**4*t1*t6**2+20._ki*t3**3*t2**4*t6*t4+25._ki*t6*t5*t2**2*t3**4-28&
                  &5._ki*t2*t3**3*t6**3*t4**2+75._ki*t2**2*t3**3*t6**3*t4-105._ki*t3*&
                  &*3*t6**3*t1*t4+84._ki*t6*t4**3*t5**2*t2**4*t3+540._ki*t6*t4**2*t3&
                  &**2*t5**2*t2**3+175._ki*t6*t4**3*t3*t5**4*t2**2+280._ki*t6*t4**3*&
                  &t3*t5**3*t2**3+48._ki*t6**3*t1*t4**4*t5*t2**2-1125._ki*t6**3*t1*t&
                  &4**3*t3*t5**2
                !
                stemp5=stemp6+9._ki*t6**3*t1*t2**4*t3*t4-45._ki*t6**3*t1*t4**2*t5**&
                  &3*t2**2-54._ki*t6**3*t1*t2**3*t3*t4**2+210._ki*t6**3*t1*t4**3*t5*&
                  &*3*t2-42._ki*t6**3*t1*t4**3*t5*t2**3-3._ki*t6*t1*t4**2*t2**6*t5+3&
                  &._ki*t6*t1*t4**2*t5**2*t2**5-15._ki*t6*t1*t4**2*t5**5*t2**2-30._ki&
                  &*t6*t1*t4**2*t5**4*t2**3+2._ki*t6**2*t4**4*t5*t2**5-2._ki*t6**2*t&
                  &4**5*t5*t2**4-10._ki*t6**2*t4**5*t5**2*t2**3+10._ki*t6**2*t4**4*t&
                  &5**2*t2**4-5._ki/2._ki*t6**2*t4**3*t5**2*t2**5-36._ki*t6**3*t4**5*&
                  &t3*t2**2-10._ki*t6*t4**3*t5**4*t2**4+20._ki*t6*t4**4*t5**4*t2**3-&
                  &t6*t4**3*t5*t2**7+2._ki*t6*t4**4*t5*t2**6-5._ki*t6*t4**3*t5**5*t2&
                  &**3+10._ki*t6*t4**4*t5**5*t2**2+8._ki*t6*t4**4*t5**3*t2**4-3._ki*t&
                  &6*t4**2*t3*t2**7
                !
                stemp6=stemp5-4._ki*t6*t4**3*t5**3*t2**5+7._ki*t6*t4**3*t3*t2**6+3.&
                  &_ki*t6*t4*t3**2*t2**6-9._ki*t6*t4**2*t3**2*t2**5+t6*t4**3*t5**2*t&
                  &2**6-2._ki*t6*t4**4*t5**2*t2**5+120._ki*t6**2*t4**2*t5**3*t3*t2**&
                  &3-420._ki*t6**2*t4**3*t3*t5**2*t2**3+t2**2*t1**3*t6**3*t5+45._ki*&
                  &t6**3*t5*t1**2*t3**2-50._ki*t2**3*t3**3*t6*t5**2+40._ki*t2**3*t3*&
                  &*3*t6**2*t5-t6**2*t1**2*t3*t2**4/2._ki+30._ki*t3**3*t2*t6**3*t1+7&
                  &0._ki*t3**4*t2*t6**2*t4-20._ki*t6**2*t2**2*t1*t3**3+990._ki*t6**3*&
                  &t4**4*t3**2*t5+t6**3*t4**3*t5*t2**5-8._ki*t6**3*t4**6*t5*t2**2-5&
                  &._ki*t6**3*t4**3*t5**3*t2**3+3._ki*t6**3*t4**2*t3*t2**5-60._ki*t6*&
                  &*3*t4**5*t5**3*t2
                !
                stemp3=stemp6+30._ki*t6**3*t4**4*t5**3*t2**2-6._ki*t6**3*t4**4*t5*t&
                  &2**4+48._ki*t6**3*t4**4*t3*t2**3-21._ki*t6**3*t4**3*t3*t2**4+540.&
                  &_ki*t6**3*t4**5*t3*t5**2+12._ki*t6**3*t4**5*t5*t2**3+135._ki*t6**3&
                  &*t1**2*t4**2*t5**3-240._ki*t6**3*t1*t4**4*t5**3+3._ki/4._ki*t6**2*&
                  &t1**2*t4*t2**5+3._ki/2._ki*t6**2*t1*t4**2*t2**6-7._ki/2._ki*t6**2*t&
                  &1*t4**3*t2**5-15._ki/2._ki*t6**2*t4*t3**2*t2**5+45._ki*t6**2*t4**2&
                  &*t3**2*t2**4-3._ki/2._ki*t6**2*t4**2*t3*t2**6+7._ki*t6**2*t4**3*t3&
                  &*t2**5+10._ki*t6**2*t4**3*t5**4*t2**3-40._ki*t6**2*t4**4*t5**4*t2&
                  &**2+40._ki*t6**2*t4**5*t5**4*t2-125._ki/2._ki*t6**2*t4**3*t3**2*t2&
                  &**3-8._ki*t6**2*t4**4*t3*t2**4+10._ki*t6**2*t4**3*t5**3*t2**4-40.&
                  &_ki*t6**2*t4**4*t5**3*t2**3+40._ki*t6**2*t4**5*t5**3*t2**2-t6**2*&
                  &t4**3*t5*t2**6/2._ki
                !
                stemp4=1._ki/t2**12*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=t4**3/t2**3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/4._ki
                !
                stemp10=20._ki*t6*t1**3*t4**4*t3**2*t5*t2**5+80._ki*t6*t1**3*t4**6*&
                  &t3**2*t5*t2**3-420._ki*t6*t1**4*t4**3*t3**2*t5**2*t2**3-600._ki*t&
                  &6*t1**5*t3**3*t5*t2**2*t4+105._ki*t6*t1**5*t4**3*t5*t2**4*t3+180&
                  &._ki*t6*t1**5*t4**2*t5**3*t3*t2**3+960._ki*t6*t1**5*t4**4*t3*t5**&
                  &3*t2-300._ki*t6*t1**6*t4*t3**2*t5**2*t2+135._ki*t6*t1**5*t4*t3**2&
                  &*t5*t2**4+270._ki*t6*t1**5*t4*t3**2*t5**2*t2**3-810._ki*t6*t1**5*&
                  &t4**2*t3**2*t5*t2**3+32._ki*t2**3*t3**5*t4**8-96._ki*t1**4*t4**4*&
                  &t5*t2**4*t3**2+80._ki/3._ki*t1**2*t4**4*t3**3*t5*t2**6+175._ki/3._k&
                  &i*t1**2*t4**3*t3**4*t5*t2**5+400._ki/3._ki*t1**2*t4**6*t3**3*t5**&
                  &2*t2**3+450._ki*t1**4*t4**2*t3**3*t2**3*t5**2-100._ki*t1**2*t4**5&
                  &*t3**3*t5**2*t2**4+100._ki/3._ki*t1**2*t4**4*t3**3*t5**2*t2**5+16&
                  &0._ki/3._ki*t6*t1**2*t2**4*t4**4*t3**4+84._ki*t1**4*t4**3*t5*t2**5&
                  &*t3**2+320._ki*t1**3*t4**4*t3**3*t2**4*t5-240._ki*t1**3*t4**5*t3*&
                  &*3*t2**3*t5+10._ki*t1**3*t4**3*t3**2*t5**2*t2**6-60._ki*t1**3*t4*&
                  &*4*t3**2*t5**2*t2**5-32._ki/3._ki*t6*t1**2*t4**8*t3**3*t2**2
                !
                stemp9=40._ki/3._ki*t6*t1**2*t4**5*t3**3*t2**5-176._ki/3._ki*t6*t1**2&
                  &*t4**7*t3**4*t2-1620._ki*t6*t1**5*t4**2*t3**2*t5**2*t2**2+8._ki/3&
                  &._ki*t1**2*t4**4*t3**3*t2**7+135._ki/2._ki*t6*t1**5*t4**2*t3**2*t2&
                  &**4-15._ki/4._ki*t6*t1**5*t4**3*t5**2*t2**5-150._ki*t6*t1**6*t4*t5&
                  &*t3**2*t2**2+40._ki/3._ki*t6*t1**4*t4**3*t5**3*t3*t2**4-80._ki*t6*&
                  &t1**4*t4**4*t5**3*t3*t2**3+60._ki*t6*t1**4*t4**2*t3**2*t5**2*t2*&
                  &*4+40._ki/3._ki*t6*t1**4*t4**6*t3*t2**3*t5+160._ki*t6*t1**4*t4**5*&
                  &t5**3*t3*t2**2-720._ki*t6*t1**4*t4**5*t3**2*t5**2*t2-200._ki/3._ki&
                  &*t1**2*t4**7*t3**3*t5**2*t2**2-80._ki*t1**2*t4**5*t3**3*t5*t2**5&
                  &+2250._ki*t6*t1**5*t4**3*t3**2*t5**2*t2-630._ki*t6*t1**5*t4**3*t3&
                  &*t5**2*t2**3+720._ki*t6*t1**5*t4**4*t3*t5**2*t2**2+15._ki*t6*t1**&
                  &5*t4**3*t5**3*t2**4+t6*t1**2*t4**3*t3**3*t2**7/3._ki-10._ki/3._ki*&
                  &t6*t1**2*t4**4*t3**3*t2**6+80._ki/3._ki*t6*t1**2*t4**7*t3**3*t2**&
                  &3+400._ki/3._ki*t6*t1**2*t4**6*t3**4*t2**2+32._ki*t6*t1*t2**2*t3**&
                  &4*t4**8+4._ki/3._ki*t6*t1**2*t4**3*t5*t2**6*t3**3-40._ki*t6*t1**4*&
                  &t4**4*t3**2*t2**4-90._ki*t6*t1**4*t2**4*t3**3*t4**2+stemp10
                !
                stemp10=stemp9-50._ki*t6*t1**4*t3**4*t2**3*t4+120._ki*t6*t1**6*t4*t&
                  &5**3*t3*t2**2+90._ki*t6*t1**6*t4*t3*t5**2*t2**3+5._ki/4._ki*t1*t4*&
                  &*3*t5*t2**7*t3**4-25._ki/2._ki*t1*t4**4*t5*t2**6*t3**4+40._ki*t1*t&
                  &4**7*t3**4*t2**4-44._ki/3._ki*t6**2*t1**5*t4**6*t5*t2**2-11._ki*t6&
                  &**2*t1**5*t4**4*t5*t2**4-110._ki*t6**2*t1**5*t4**5*t5**3*t2+75._k&
                  &i*t6**2*t1**7*t4*t3*t5**2-15._ki*t6**2*t1**7*t3*t5**2*t2+400._ki*&
                  &t1**3*t4**4*t3**3*t5**2*t2**3+25._ki/3._ki*t1**3*t4**3*t3**2*t5**&
                  &3*t2**5-475._ki/2._ki*t1**4*t4**2*t3**4*t5*t2**2-90._ki*t1**4*t4**&
                  &2*t3**2*t5**2*t2**5-20._ki*t1**4*t4**3*t5**3*t3*t2**5+80._ki*t1**&
                  &4*t4**4*t5**3*t3*t2**4+22._ki*t6**2*t1**5*t4**5*t5*t2**3-66._ki*t&
                  &6**2*t1**5*t4**5*t3*t2**2-300._ki*t6**2*t1**4*t4**6*t3**2*t5+t6*&
                  &*2*t1**4*t4**3*t3*t2**6/4._ki-375._ki/2._ki*t6**2*t1**4*t4**3*t3**&
                  &3*t2**2+320._ki/3._ki*t6*t1**2*t4**7*t3**3*t5*t2**2-140._ki*t1**3*&
                  &t4**3*t3**3*t2**5*t5+20._ki*t1**3*t4**2*t3**3*t2**6*t5+120._ki*t1&
                  &**3*t4**5*t3**2*t5**2*t2**4
                !
                stemp8=2._ki*t1**3*t3**2*t4**3*t5*t2**7-16._ki*t1**3*t3**2*t4**6*t5&
                  &*t2**4-80._ki*t1**3*t4**6*t3**2*t5**2*t2**3+160._ki/3._ki*t6*t1**2&
                  &*t4**5*t3**3*t5*t2**4-225._ki/2._ki*t1**3*t4**2*t3**4*t2**4*t5-16&
                  &5._ki/2._ki*t6**2*t1**5*t4**2*t3*t5**2*t2**3+3._ki*t6**2*t1**7*t4*&
                  &t5*t2**3-15._ki*t6**2*t1**7*t4*t5**3*t2+675._ki*t6**2*t1**6*t4**2&
                  &*t3*t5**2*t2+55._ki*t6**2*t1**5*t4**4*t5**3*t2**2+990._ki*t6**2*t&
                  &1**5*t4**5*t3*t5**2-40._ki*t2**4*t3**5*t4**7+5._ki/6._ki*t1**3*t3*&
                  &*5*t2**5+2._ki*t4**4*t3**5*t2**7-5._ki*t1**4*t3**4*t2**5+80._ki/3.&
                  &_ki*t2**5*t3**5*t4**6-t4**3*t3**5*t2**8/6._ki-10._ki*t4**5*t3**5*t&
                  &2**6-4._ki*t1**5*t3**3*t2**5+350._ki*t1**4*t4**3*t3**2*t5**3*t2**&
                  &3-80._ki*t1**4*t4**5*t5**3*t3*t2**3-50._ki*t1**4*t4**5*t5**4*t3*t&
                  &2**2-25._ki/2._ki*t1**4*t4**3*t5**4*t3*t2**4-880._ki*t6*t1**4*t4**&
                  &4*t3**3*t5*t2+420._ki*t1**4*t4**3*t3**2*t5**2*t2**4+24._ki*t1**4*&
                  &t4**4*t3*t2**5*t5**2-270._ki*t6*t1**6*t4**2*t3*t5**2*t2**2-15._ki&
                  &*t6*t1**6*t4*t5*t2**4*t3+stemp10
                !
                stemp10=stemp8-1425._ki/2._ki*t6**2*t1**6*t4**2*t3**2*t5-45._ki*t6**&
                  &2*t1**6*t2**3*t3*t4**2+70._ki/3._ki*t6*t1**3*t4**3*t3**3*t2**5+t6&
                  &*t1**4*t4**4*t3*t2**6-210._ki*t6*t1**4*t4**3*t3**4*t2-20._ki/3._ki&
                  &*t6*t1**3*t4**6*t3**2*t2**4-1045._ki/2._ki*t6**2*t1**5*t3**3*t4**&
                  &2*t2-5._ki/2._ki*t6*t1**3*t4**2*t3**3*t2**6+30._ki*t6*t1**3*t4**2*&
                  &t3**4*t2**4+5._ki*t6*t1**3*t4**5*t3**2*t2**5-5._ki/3._ki*t6*t1**3*&
                  &t4**4*t3**2*t2**6+5._ki/24._ki*t6*t1**3*t4**3*t3**2*t2**7+21._ki/2&
                  &._ki*t6*t1**5*t4**3*t3*t2**5+15._ki*t6*t1**5*t4**3*t5**4*t2**3+60&
                  &._ki*t6*t1**5*t4**5*t5**4*t2+5._ki*t1**4*t3**5*t2**3+10._ki*t1**5*&
                  &t2**3*t3**4-5._ki*t1*t4**4*t3**4*t2**7+20._ki*t1*t4**5*t3**4*t2**&
                  &6+t1*t4**3*t3**4*t2**8/2._ki+100._ki*t1*t2**3*t3**5*t4**6+2._ki*t1&
                  &**4*t4**4*t3*t2**7-t1**4*t4**3*t3*t2**8/2._ki-35._ki/2._ki*t1**4*t&
                  &3**5*t2**2*t4-2._ki*t1**4*t4**5*t3*t2**6-7._ki*t1**4*t4**3*t3**2*&
                  &t2**6
                !
                stemp9=t1**3*t4**4*t3**2*t2**7+5._ki*t1**3*t4*t3**4*t2**6-25._ki/2.&
                  &_ki*t1**3*t3**5*t2**4*t4-40._ki*t1*t4**6*t3**4*t2**5-44._ki*t1*t2*&
                  &*2*t3**5*t4**7+3._ki/4._ki*t1*t2**7*t3**5*t4**2-16._ki*t1*t4**8*t3&
                  &**4*t2**3+7._ki*t1**5*t4**3*t3*t2**6+20._ki*t1**5*t3**3*t2**4*t4-&
                  &4._ki*t1**5*t4**3*t5**3*t2**5+50._ki*t1**4*t4*t3**4*t2**4-25._ki/2&
                  &._ki*t1**4*t3**4*t5*t2**4-9._ki*t1**5*t4**2*t3**2*t2**5+36._ki*t1*&
                  &*4*t4**2*t3**3*t2**5-50._ki*t1**4*t3**3*t4**3*t2**4-2._ki*t1**5*t&
                  &4**4*t5**2*t2**5+95._ki/2._ki*t1**3*t3**5*t2**3*t4**2+125._ki*t1**&
                  &3*t4**3*t3**4*t2**4+4._ki/3._ki*t1**3*t4**6*t3**2*t2**5-110._ki*t1&
                  &**3*t4**4*t3**4*t2**3+2._ki*t1**3*t4**2*t3**3*t2**7-105._ki/2._ki*&
                  &t1**3*t3**5*t2**2*t4**3-t1**3*t4**3*t3**2*t2**8/6._ki-2._ki*t1**3&
                  &*t4**5*t3**2*t2**6-45._ki*t1**3*t4**2*t3**4*t2**5-14._ki*t1**3*t4&
                  &**3*t3**3*t2**6+32._ki*t1**3*t4**4*t3**3*t2**5+stemp10
                !
                stemp10=stemp9-24._ki*t1**3*t4**5*t3**3*t2**4+3._ki/4._ki*t6*t1**6*t&
                  &4**2*t2**6+10._ki/3._ki*t6*t1**4*t3**4*t2**4-5._ki/2._ki*t6*t1**6*t&
                  &3**2*t2**4+t6**2*t1**7*t3*t2**3+45._ki*t6**2*t1**7*t4**2*t5**3-1&
                  &75._ki/2._ki*t6**2*t1**6*t3**3*t4+1155._ki/2._ki*t6**2*t1**5*t4**3*&
                  &t3**3+25._ki*t6**2*t1**6*t3**3*t2+15._ki*t6*t1**5*t4**4*t5**2*t2*&
                  &*4-15._ki*t6*t1**5*t4**5*t5**2*t2**3-120._ki*t6*t1**5*t4**4*t5*t2&
                  &**3*t3+135._ki*t6*t1**5*t4**2*t3*t5**2*t2**4+1125._ki*t6*t1**5*t4&
                  &**3*t3**2*t5*t2**2-165._ki/2._ki*t6**2*t1**5*t3**2*t2**3*t4*t5-77&
                  &._ki/2._ki*t6**2*t1**5*t4**3*t3*t2**4+1815._ki*t6**2*t1**5*t4**4*t&
                  &3**2*t5+11._ki/2._ki*t6**2*t1**5*t4**2*t3*t2**5+25._ki*t1**3*t4**2&
                  &*t3**3*t5**2*t2**5-500._ki*t1**4*t4**3*t3**3*t5*t2**3+625._ki/2._k&
                  &i*t1**3*t4**3*t3**4*t2**3*t5-625._ki*t1**4*t4**3*t3**3*t5**2*t2*&
                  &*2+50._ki*t1**4*t4**4*t5**4*t3*t2**3-400._ki*t1**4*t4**4*t3**2*t5&
                  &**3*t2**2-18._ki*t1**4*t4**2*t5*t2**6*t3**2-480._ki*t1**4*t4**4*t&
                  &3**2*t5**2*t2**3
                !
                stemp7=-500._ki/3._ki*t1**2*t4**6*t3**4*t5*t2**2-160._ki/3._ki*t1**2*&
                  &t4**7*t3**3*t5*t2**3+8._ki*t1**4*t4**4*t3**2*t2**5-40._ki*t1*t4**&
                  &8*t3**4*t5*t2**2+100._ki*t1*t4**7*t3**4*t5*t2**3+275._ki/2._ki*t6*&
                  &*2*t1**5*t3**3*t2**2*t4+88._ki*t6**2*t1**5*t4**4*t3*t2**3+11._ki/&
                  &6._ki*t6**2*t1**5*t4**3*t5*t2**5-55._ki/6._ki*t6**2*t1**5*t4**3*t5&
                  &**3*t2**3-320._ki/3._ki*t6*t1**4*t4**6*t5**3*t3*t2-6._ki*t1**4*t4*&
                  &*3*t3*t2**6*t5**2-24._ki*t1**4*t4**5*t3*t2**4*t5**2-60._ki*t1**4*&
                  &t4*t3**3*t5*t2**5-120._ki*t6*t1**2*t4**5*t3**4*t2**3+t6*t1**2*t4&
                  &**2*t3**4*t2**6+375._ki*t6**2*t1**6*t3**2*t5*t2*t4-45._ki/4._ki*t6&
                  &**2*t1**4*t4**2*t3**2*t5*t2**4+60._ki*t6*t1**5*t4**5*t5**3*t2**2&
                  &-60._ki*t6*t1**5*t4**4*t5**3*t2**3-3._ki*t6*t1**5*t4**5*t5*t2**4-&
                  &12._ki*t6*t1**5*t4**4*t3*t2**4-9._ki/4._ki*t6*t1**5*t4**2*t3*t2**6&
                  &+60._ki*t6*t1**5*t3**3*t5*t2**3-375._ki/4._ki*t6*t1**5*t4**3*t3**2&
                  &*t2**3-40._ki*t6*t1**3*t4**7*t3**2*t5*t2**2-45._ki/2._ki*t6*t1**5*&
                  &t4**2*t5*t2**5*t3-840._ki*t6*t1**5*t4**3*t3*t5**3*t2**2+3._ki*t6*&
                  &t1**5*t4**4*t5*t2**5+stemp10
                !
                stemp10=stemp7-3._ki/4._ki*t6*t1**5*t4**3*t5*t2**6-128._ki/3._ki*t6*t&
                  &1**2*t4**8*t3**3*t5*t2-40._ki/3._ki*t6*t1**2*t4**4*t5*t2**5*t3**3&
                  &+2._ki*t6*t1*t4**4*t3**4*t2**6+80._ki/3._ki*t6*t1*t2**4*t3**4*t4**&
                  &6+1155._ki/2._ki*t6**2*t1**5*t4**3*t3*t5**2*t2**2-1320._ki*t6**2*t&
                  &1**5*t4**4*t3*t5**2*t2+1485._ki/2._ki*t6**2*t1**5*t4**2*t3**2*t5*&
                  &t2**2-t6*t1*t4**3*t3**4*t2**7/6._ki-120._ki*t6*t1**3*t4**5*t3**2*&
                  &t5**2*t2**3+160._ki*t6*t1**3*t4**6*t3**2*t5**2*t2**2-80._ki*t6*t1&
                  &**3*t4**7*t3**2*t5**2*t2+480._ki*t6*t1**3*t4**5*t3**3*t5*t2**2-8&
                  &00._ki/3._ki*t6*t1**3*t4**6*t3**3*t5*t2+280._ki/3._ki*t6*t1**3*t4**&
                  &3*t3**3*t5*t2**4+220._ki*t6*t1**3*t4**4*t3**4*t2**2-200._ki/3._ki*&
                  &t6*t1**3*t4**6*t3**3*t2**2-125._ki*t6*t1**3*t4**3*t3**4*t2**3-14&
                  &0._ki*t6*t1**3*t4**5*t3**4*t2-5._ki/2._ki*t6*t1**3*t3**4*t4*t2**5+&
                  &10._ki/3._ki*t6*t1**3*t4**7*t3**2*t2**3+120._ki*t6*t1**3*t4**5*t3*&
                  &*3*t2**3-80._ki*t6*t1**3*t4**4*t3**3*t2**4-275._ki*t1**3*t4**4*t3&
                  &**4*t2**2*t5-300._ki*t1**3*t4**5*t3**3*t5**2*t2**2+108._ki*t1**5*&
                  &t4**2*t3**2*t5*t2**4
                !
                stemp9=10._ki*t6*t1**4*t4**3*t3*t5**2*t2**5-360._ki*t6*t1**4*t4**2*&
                  &t3**3*t2**3*t5-320._ki/3._ki*t6*t1**2*t4**6*t3**3*t5*t2**3-10._ki*&
                  &t6*t1*t4**5*t3**4*t2**5+1140._ki*t6*t1**5*t4**2*t3**3*t5*t2-100.&
                  &_ki*t1*t4**6*t3**4*t5*t2**4-75._ki*t1**4*t4**2*t3**2*t5**3*t2**4+&
                  &stemp10+175._ki*t1**5*t4**3*t3*t5**4*t2**2-320._ki*t6*t1**3*t4**4&
                  &*t3**3*t5*t2**3-35._ki/3._ki*t6*t1**2*t3**4*t4**3*t2**5-t6*t1**4*&
                  &t4**3*t3*t2**7/6._ki-220._ki*t6*t1**4*t4**4*t3**3*t2**2+540._ki*t1&
                  &**5*t4**2*t3**2*t5**2*t2**3-36._ki*t1**5*t4**2*t5**2*t2**5*t3-36&
                  &._ki*t1**5*t4*t5*t2**5*t3**2-75._ki*t1**5*t4**2*t5**4*t3*t2**3-5.&
                  &_ki/3._ki*t6*t1**4*t4**3*t3*t2**6*t5+1000._ki*t6*t1**4*t4**3*t3**3&
                  &*t2**2*t5-60._ki*t6*t1**4*t4**4*t3*t5**2*t2**4+30._ki*t6*t1**4*t4&
                  &**2*t3**2*t5*t2**5+120._ki*t6*t1**4*t4**5*t3*t5**2*t2**3-360._ki*&
                  &t6*t1**6*t4**2*t5**3*t2*t3+45._ki*t6*t1**6*t4**2*t5*t2**3*t3+60.&
                  &_ki*t6*t1**6*t3**2*t5**2*t2**2-10._ki/3._ki*t1**2*t4**3*t5*t3**3*t&
                  &2**7-25._ki/4._ki*t1**2*t4**2*t3**4*t5*t2**6+45._ki*t6**2*t1**4*t4&
                  &**2*t3**3*t2**3
                !
                stemp10=stemp9+6._ki*t6**2*t1**4*t4**5*t3*t2**4+4._ki*t6**2*t1**4*t&
                  &4**7*t3*t2**2-60._ki*t6*t1**5*t4**4*t5**4*t2**2+960._ki*t6*t1**4*&
                  &t4**4*t3**2*t5**2*t2**2-35._ki*t6**2*t1**6*t4**3*t5*t2**3-20._ki*&
                  &t6*t1**4*t4**5*t3*t2**4*t5+10._ki*t6*t1**4*t4**4*t3*t2**5*t5-2._k&
                  &i*t6*t1**4*t4**5*t3*t2**5-150._ki*t6*t1**5*t4*t3**3*t2**3+190._ki&
                  &*t6*t1**4*t3**4*t2**2*t4**2+30._ki*t6*t1**4*t4**5*t3**2*t2**3+35&
                  &._ki/2._ki*t6*t1**4*t4**3*t3**2*t2**5-5._ki/2._ki*t6*t1**4*t4**2*t3&
                  &**2*t2**6+10._ki*t6*t1**4*t2**5*t3**3*t4-1875._ki/2._ki*t6**2*t1**&
                  &6*t4**3*t3*t5**2-75._ki/2._ki*t6**2*t1**6*t3**2*t5*t2**2+15._ki/2.&
                  &_ki*t6**2*t1**6*t4**2*t5*t2**4+15._ki/2._ki*t6**2*t1**6*t2**4*t3*t&
                  &4+40._ki*t6**2*t1**6*t4**4*t5*t2**2-225._ki/2._ki*t6**2*t1**6*t4*t&
                  &3*t5**2*t2**2-5._ki/2._ki*t6*t1**3*t4**3*t5*t3**2*t2**6-5._ki*t6*t&
                  &1**3*t4**3*t5**2*t2**5*t3**2-4125._ki/2._ki*t6**2*t1**5*t4**3*t3*&
                  &*2*t5*t2-32._ki/3._ki*t2**2*t3**5*t4**9-75._ki/2._ki*t6**2*t1**6*t4&
                  &**2*t5**3*t2**2+175._ki*t6**2*t1**6*t4**3*t5**3*t2
                !
                stemp8=stemp10+125._ki/2._ki*t6**2*t1**6*t2**2*t3*t4**3+4._ki/3._ki*t&
                  &6*t1**4*t4**6*t3*t2**4+30._ki*t6*t1**6*t4**2*t5**4*t2**2+7._ki/2.&
                  &_ki*t6*t1**6*t4**3*t5*t2**4-40._ki*t6*t1**6*t3**3*t5*t2-3._ki/2._ki&
                  &*t6*t1**6*t4**2*t5*t2**5-55._ki/6._ki*t6**2*t1**5*t2**3*t3**3+220&
                  &._ki/3._ki*t6**2*t1**5*t4**6*t5**3-210._ki*t6**2*t1**4*t4**5*t3**3&
                  &-44._ki*t6**2*t1**3*t4**7*t3**3-16._ki/3._ki*t6**2*t1**2*t4**9*t3*&
                  &*3-10._ki*t6*t1**6*t3**3*t2**2+3._ki/2._ki*t6*t1**5*t4**5*t2**5-3.&
                  &_ki/2._ki*t6*t1**5*t4**4*t2**6+3._ki/8._ki*t6*t1**5*t4**3*t2**7-30.&
                  &_ki*t6*t1**5*t3**4*t2**2+15._ki*t6*t1**5*t3**3*t2**4-7._ki/4._ki*t6&
                  &*t1**6*t4**3*t2**5+15._ki*t6**2*t1**7*t3**2*t5-200._ki*t6**2*t1**&
                  &6*t4**4*t5**3-95._ki*t1**4*t4**2*t3**4*t2**3-15._ki/4._ki*t6**2*t1&
                  &**4*t3**3*t2**4*t4-8._ki*t6**2*t1**4*t4**6*t3*t2**3+40._ki*t6*t1*&
                  &*4*t4*t3**3*t5*t2**4-80._ki*t6*t1**4*t4**6*t3*t5**2*t2**2+250._ki&
                  &*t6*t1**4*t2**3*t3**3*t4**3+250._ki*t1**5*t4*t3**3*t5**2*t2**2+4&
                  &50._ki*t1**5*t4**2*t5**3*t3**2*t2**2
                !
                stemp10=stemp8+200._ki*t1**5*t4*t5*t3**3*t2**3-150._ki*t1**5*t4*t5*&
                  &*3*t3**2*t2**3-40._ki*t6*t1*t2**3*t3**4*t4**7-32._ki/3._ki*t6*t1*t&
                  &2*t3**4*t4**9-2._ki*t6**2*t1**4*t4**4*t3*t2**5-60._ki*t6**2*t1**4&
                  &*t4**7*t3*t5**2+330._ki*t6**2*t1**4*t4**4*t3**3*t2+40._ki*t6**2*t&
                  &1**3*t4**4*t3**3*t2**3+3._ki/4._ki*t6**2*t1**3*t4**2*t3**3*t2**5+&
                  &100._ki*t6**2*t1**3*t4**6*t3**3*t2-14._ki*t1**5*t4**3*t5*t3*t2**5&
                  &+6._ki*t1**5*t4**2*t5*t3*t2**6-120._ki*t1**5*t4**2*t3*t5**3*t2**4&
                  &+50._ki*t1*t4**5*t3**4*t5*t2**5+280._ki*t1**5*t4**3*t3*t5**3*t2**&
                  &3+100._ki*t1**3*t4**5*t3**2*t5**3*t2**3-80._ki/3._ki*t6*t1**2*t4**&
                  &6*t3**3*t2**4-200._ki/3._ki*t1**3*t4**6*t3**2*t5**3*t2**2-50._ki*t&
                  &1**3*t4**4*t3**2*t5**3*t2**4+30._ki*t6*t1**6*t4**2*t5**3*t2**3-7&
                  &0._ki*t6*t1**6*t4**3*t5**4*t2+35._ki/2._ki*t6*t1**6*t4**3*t5**2*t2&
                  &**3-70._ki*t6*t1**6*t4**3*t5**3*t2**2+285._ki*t6*t1**5*t3**3*t2**&
                  &2*t4**2+105._ki*t6*t1**5*t3**4*t2*t4-45._ki/4._ki*t6*t1**5*t4*t3**&
                  &2*t2**5
                !
                stemp9=stemp10-5._ki/2._ki*t1**2*t4**2*t3**4*t2**7+70._ki/3._ki*t1**2&
                  &*t4**3*t3**4*t2**6-200._ki/3._ki*t1**2*t4**6*t3**4*t2**3-t1**2*t4&
                  &**3*t3**3*t2**8/3._ki-16._ki/3._ki*t1**2*t4**7*t3**3*t2**4+32._ki/3&
                  &._ki*t1**2*t4**6*t3**3*t2**5-8._ki*t1**2*t4**5*t3**3*t2**6-80._ki*&
                  &t1**2*t4**4*t3**4*t2**5+120._ki*t1**2*t4**5*t3**4*t2**4-70._ki*t1&
                  &**2*t2**2*t3**5*t4**5+15._ki*t1**2*t2**5*t3**5*t4**2+110._ki*t1**&
                  &2*t2**3*t3**5*t4**4-5._ki/4._ki*t1**2*t4*t3**5*t2**6-125._ki/2._ki*&
                  &t1**2*t2**4*t3**5*t4**3-35._ki/4._ki*t1*t2**6*t3**5*t4**3-50._ki*t&
                  &1**5*t3**3*t5**2*t2**3-3._ki*t1**5*t4**2*t3*t2**7+3._ki*t1**5*t4*&
                  &t3**2*t2**6-10._ki*t1**5*t4**3*t5**4*t2**4+t1**5*t4**3*t5**2*t2*&
                  &*6-40._ki*t1**5*t3**3*t5*t2**4+20._ki*t1**5*t4**4*t5**4*t2**3+2._k&
                  &i*t1**5*t4**4*t5*t2**6+25._ki*t1**5*t3**4*t5*t2**2-5._ki*t1**5*t4&
                  &**3*t5**5*t2**3+10._ki*t1**5*t4**4*t5**5*t2**2+3._ki/2._ki*t1**4*t&
                  &4**2*t3**2*t2**7
                !
                stemp10=stemp9-6._ki*t1**4*t4*t3**3*t2**6-t1**5*t4**3*t5*t2**7-90.&
                  &_ki*t1*t2**4*t3**5*t4**5+40._ki*t1*t2**5*t3**5*t4**4+540._ki*t6**2&
                  &*t1**4*t4**5*t3**2*t5*t2-360._ki*t6**2*t1**4*t4**4*t3**2*t5*t2**&
                  &2+8._ki*t1**5*t4**4*t5**3*t2**4-10._ki*t6*t1**3*t4**2*t3**3*t5*t2&
                  &**5+40._ki*t6*t1**3*t4**4*t3**2*t5**2*t2**4-60._ki*t6*t1**3*t4**5&
                  &*t3**2*t5*t2**4+120._ki*t6**2*t1**4*t4**6*t3*t5**2*t2-15._ki/2._ki&
                  &*t6**2*t1**3*t4**4*t5*t2**4*t3**2-60._ki*t6**2*t1**3*t4**6*t3**2&
                  &*t5*t2**2+3._ki/4._ki*t6**2*t1**3*t4**3*t5*t2**5*t3**2-24._ki*t6**&
                  &2*t1**3*t4**8*t3**2*t5-180._ki*t1**5*t4*t3**2*t5**2*t2**4+360._ki&
                  &*t1**4*t4**2*t3**3*t5*t2**4+30._ki*t6**2*t1**3*t4**5*t3**2*t5*t2&
                  &**3+60._ki*t6**2*t1**3*t4**7*t3**2*t5*t2-5._ki*t6**2*t1**7*t2**2*&
                  &t3*t4-9._ki*t6**2*t1**7*t4**2*t5*t2**2-35._ki/4._ki*t6**2*t1**3*t3&
                  &**3*t4**3*t2**4+105._ki*t6**2*t1**4*t4**3*t5*t3**2*t2**3-15._ki/4&
                  &._ki*t6**2*t1**4*t4**3*t5**2*t2**4*t3+30._ki*t6**2*t1**4*t4**4*t5&
                  &**2*t2**3*t3-90._ki*t6**2*t1**4*t4**5*t3*t5**2*t2**2+84._ki*t1**5&
                  &*t4**3*t5**2*t2**4*t3
                !
                stemp6=stemp10-4._ki*t1**4*t4**4*t3*t2**6*t5+125._ki*t1**4*t3**4*t5&
                  &*t2**3*t4-75._ki*t1**4*t4*t3**3*t2**4*t5**2+t1**4*t4**3*t3*t2**7&
                  &*t5+4._ki*t1**4*t4**5*t3*t2**5*t5+30._ki*t6*t1**6*t2**3*t3**2*t5+&
                  &25._ki/2._ki*t6*t1**6*t2**3*t4*t3**2+9._ki/2._ki*t6*t1**6*t4**2*t3*&
                  &t2**4-15._ki/2._ki*t6*t1**6*t4**2*t5**2*t2**4-3._ki/2._ki*t6*t1**6*&
                  &t4*t3*t2**5-210._ki*t6*t1**4*t4**3*t3**2*t5*t2**4+480._ki*t6*t1**&
                  &4*t4**4*t3**2*t5*t2**3-360._ki*t6*t1**4*t4**5*t3**2*t5*t2**2+25.&
                  &_ki/2._ki*t1**3*t4*t3**4*t2**5*t5-12._ki*t1**3*t3**2*t4**4*t5*t2**&
                  &6+24._ki*t1**3*t3**2*t4**5*t5*t2**5-175._ki*t1**3*t4**3*t3**3*t5*&
                  &*2*t2**4-200._ki*t1**2*t4**4*t3**4*t5*t2**4-25._ki/6._ki*t1**2*t4*&
                  &*3*t5**2*t2**6*t3**3+320._ki/3._ki*t1**2*t4**6*t3**3*t5*t2**4+300&
                  &._ki*t1**2*t4**5*t3**4*t5*t2**3+t6**2*t1**2*t4**4*t3**3*t2**5+40&
                  &._ki/3._ki*t6**2*t1**2*t4**6*t3**3*t2**3+16._ki*t6**2*t1**2*t4**8*&
                  &t3**3*t2-5._ki*t6**2*t1**2*t4**5*t3**3*t2**4-90._ki*t6**2*t1**3*t&
                  &4**5*t3**3*t2**2-20._ki*t6**2*t1**2*t4**7*t3**3*t2**2-t6**2*t1**&
                  &2*t4**3*t3**3*t2**6/12._ki
                !
                stemp7=t6/t2**12/t1**5
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(3)
                !
                stemp6=-30._ki*t6**3*t1*t4*t5**3*t2**3-250._ki*t2**2*t3**3*t6*t4*t5&
                  &**2-5._ki/4._ki*t6**2*t4**3*t2**7-45._ki*t6**3*t1*t3*t5**2*t2**3-t&
                  &6*t1*t2**5*t3**2-190._ki*t3**3*t2**2*t6**2*t4**2-6._ki*t3*t2**3*t&
                  &6**3*t1**2+10._ki*t3**2*t2**4*t1*t6**2-20._ki*t3**3*t2**4*t6*t4-2&
                  &5._ki*t6*t5*t2**2*t3**4+410._ki*t2*t3**3*t6**3*t4**2-165._ki*t2**2&
                  &*t3**3*t6**3*t4+70._ki*t6**2*t1*t4*t5*t2**4*t3-30._ki*t6**3*t4*t3&
                  &*t5**2*t2**4+2._ki*t6*t1*t4*t5**2*t2**6-9._ki*t6**2*t1*t3*t2**4*t&
                  &4**2+7._ki*t6**2*t1*t3*t2**5*t4-1395._ki*t6**3*t1*t4**2*t3*t5**2*&
                  &t2-24._ki*t6*t4*t5**2*t2**6*t3-25._ki*t6*t1*t5**4*t3*t2**3+6._ki*t&
                  &6**3*t1*t4*t5*t2**5+80._ki*t6**2*t1*t3*t5**3*t2**3+60._ki*t6**2*t&
                  &4*t3*t5**2*t2**5-10._ki*t6**2*t1*t4*t5**2*t2**5-40._ki*t6*t1*t5**&
                  &3*t3*t2**4-30._ki*t3*t5**2*t2**2*t1**2*t6**2+60._ki*t3**2*t5**2*t&
                  &2**3*t1*t6-10._ki*t6**2*t1*t3*t5*t2**5+2._ki*t6*t1*t2**6*t3*t5+80&
                  &._ki*t3**3*t2**4*t6*t5
                !
                stemp5=80._ki*t6**2*t4**4*t5**3*t2**3-40._ki*t6**2*t4**5*t5**3*t2**&
                  &2+5._ki/2._ki*t6**2*t4**3*t5*t2**6-4._ki*t6**2*t4**4*t5*t2**5+2._ki&
                  &*t6**2*t4**5*t5*t2**4+10._ki*t6**2*t4**5*t5**2*t2**3-20._ki*t6**2&
                  &*t4**4*t5**2*t2**4+5._ki*t1**3*t6**3*t5**3-8._ki*t6*t1*t4*t5**3*t&
                  &2**5-20._ki*t6*t1*t4*t5**4*t2**4+40._ki*t6**2*t1*t4*t5**4*t2**3-2&
                  &1._ki*t6**3*t1**2*t4*t5*t2**3-4._ki*t6*t4**2*t5**3*t2**6+20._ki*t6&
                  &**2*t2**2*t1*t3**3-90._ki*t6**3*t4**4*t5**3*t2**2+18._ki*t6**3*t4&
                  &**4*t5*t2**4-76._ki*t6**3*t4**4*t3*t2**3+57._ki*t6**3*t4**3*t3*t2&
                  &**4-10._ki*t6*t1*t4*t5**5*t2**3-40._ki*t6**3*t4**6*t5**3+27._ki*t6&
                  &**3*t1**2*t4**2*t5*t2**2+105._ki*t6**3*t1**2*t4*t5**3*t2-t6**2*t&
                  &4**5*t2**5+540._ki*t6**2*t1*t4**2*t3*t5**2*t2**2+720._ki*t6**2*t1&
                  &*t4**2*t5**3*t2*t3+105._ki*t3**3*t6**3*t1*t4-t2**2*t1**3*t6**3*t&
                  &5-45._ki*t6**3*t5*t1**2*t3**2-8._ki*t6*t4**4*t5**3*t2**4+8._ki*t6*&
                  &t4**2*t3*t2**7+12._ki*t6*t4**3*t5**3*t2**5+stemp6
                !
                stemp6=-t6*t1*t3*t2**7+90._ki*t2*t3*t1**2*t6**3*t5**2-560._ki*t6**2&
                  &*t1*t4*t5**3*t3*t2**2-660._ki*t6**2*t4*t3**2*t5**2*t2**3-48._ki*t&
                  &6**3*t1*t4**4*t5*t2**2-75._ki*t6**3*t1*t2**2*t3*t4**3+36._ki*t6*t&
                  &1*t4*t3*t2**4*t5**2+1125._ki*t6**3*t1*t4**3*t3*t5**2+320._ki*t6*t&
                  &4**2*t3*t5**3*t2**4+350._ki*t6*t4*t5**3*t3**2*t2**3-240._ki*t2**2&
                  &*t3**2*t1*t6**2*t5**2-t4**2*t2**10/4._ki+50._ki*t2**2*t3**2*t1*t6&
                  &*t5**3-120._ki*t2**3*t3**3*t6**2*t5-t6*t4**2*t2**8*t5-60._ki*t6*t&
                  &3**2*t2**5*t5**2+t6*t4**2*t5**2*t2**7-33._ki*t6**3*t1*t2**4*t3*t&
                  &4+200._ki*t6*t4**2*t5**4*t3*t2**3+84._ki*t6*t4*t5*t2**5*t3**2+855&
                  &._ki*t3**2*t1*t6**3*t5*t4**2+8._ki*t3**3*t6*t2**5-315._ki*t3**3*t6&
                  &**3*t4**3-540._ki*t6**3*t4**5*t3*t5**2-20._ki*t6**3*t4**5*t5*t2**&
                  &3-135._ki*t6**3*t1**2*t4**2*t5**3+240._ki*t6**3*t1*t4**4*t5**3-3.&
                  &_ki/4._ki*t6**2*t1**2*t4*t2**5-4._ki*t6**2*t1*t4**2*t2**6+7._ki/2._k&
                  &i*t6**2*t1*t4**3*t2**5-108._ki*t6*t4**2*t3**2*t5*t2**4
                !
                stemp4=-5._ki/2._ki*t6**2*t3**2*t2**6-10._ki*t6*t2**3*t3**4+t6**2*t1&
                  &**2*t2**6/4._ki+96._ki*t6*t4**2*t5**2*t2**5*t3-10._ki*t6**2*t4*t5*&
                  &t2**6*t3-520._ki*t6**2*t4**2*t5**3*t3*t2**3-7._ki*t6*t4**3*t3*t2*&
                  &*6-7._ki*t6*t4*t3**2*t2**6+9._ki*t6*t4**2*t3**2*t2**5-3._ki*t6*t4*&
                  &*3*t5**2*t2**6-84._ki*t6*t4**3*t5**2*t2**4*t3-540._ki*t6*t4**2*t3&
                  &**2*t5**2*t2**3-5._ki*t6*t4**2*t5**5*t2**4-990._ki*t6**3*t4**4*t3&
                  &**2*t5-7._ki*t6**3*t4**3*t5*t2**5+8._ki*t6**3*t4**6*t5*t2**2+35._k&
                  &i*t6**3*t4**3*t5**3*t2**3-5._ki/2._ki*t6**2*t4**2*t5**2*t2**6-175&
                  &._ki*t6*t4**3*t3*t5**4*t2**2+195._ki*t6**3*t1*t4**2*t5**3*t2**2+9&
                  &3._ki*t6**3*t1*t2**3*t3*t4**2-450._ki*t6*t4**2*t5**3*t3**2*t2**2-&
                  &750._ki*t6**2*t4**3*t3**2*t5*t2**2-390._ki*t6**3*t1*t4**3*t5**3*t&
                  &2+78._ki*t6**3*t1*t4**3*t5*t2**3-3._ki*t6*t1*t4**2*t5**2*t2**5+12&
                  &._ki*t6*t1*t4**2*t5**3*t2**4-50._ki*t6*t3**2*t5**3*t2**4-2._ki*t6*&
                  &t4*t3*t2**8-10._ki*t6*t4**2*t5**4*t2**5+stemp6+stemp5
                !
                stemp6=10._ki*t6**2*t4**2*t5**3*t2**5-t6**2*t4**2*t5*t2**7/2._ki-t6&
                  &**2*t4*t3*t2**7+25._ki/2._ki*t6**2*t4**3*t5**2*t2**5+30._ki*t6*t4*&
                  &*3*t5**4*t2**4-20._ki*t6*t4**4*t5**4*t2**3+3._ki*t6*t4**3*t5*t2**&
                  &7-2._ki*t6*t4**4*t5*t2**6+15._ki*t6*t4**3*t5**5*t2**3-10._ki*t6*t4&
                  &**4*t5**5*t2**2+100._ki*t2**3*t3**3*t6*t5**2+20._ki*t2**3*t3**3*t&
                  &6**3-640._ki*t6**2*t4**4*t3*t5**3*t2+55._ki/2._ki*t6**2*t4*t3**2*t&
                  &2**5-155._ki/2._ki*t6**2*t4**2*t3**2*t2**4+13._ki/2._ki*t6**2*t4**2&
                  &*t3*t2**6+2._ki*t6**3*t4*t3*t2**6+60._ki*t6**2*t3**2*t2**4*t5**2+&
                  &10._ki*t6**2*t4**2*t5**4*t2**4-t6**2*t1**2*t5*t2**5/2._ki+30._ki*t&
                  &6**2*t2**5*t3**2*t5+t6**2*t1*t4*t2**7-15._ki*t6**3*t3**2*t2**4*t&
                  &5-15._ki*t6**3*t1**2*t5**3*t2**2+15._ki*t6*t1*t4**2*t5**5*t2**2+3&
                  &0._ki*t6*t1*t4**2*t5**4*t2**3-30._ki*t3**3*t2**4*t6**2+30._ki*t3**&
                  &4*t2**2*t6**2+65._ki*t6**2*t4**2*t5*t2**5*t3-39._ki*t6**3*t1*t4**&
                  &2*t5*t2**4
                !
                stemp5=-280._ki*t6*t4**3*t3*t5**3*t2**3+270._ki*t6**3*t4**2*t3*t5**&
                  &2*t2**3-855._ki*t6**3*t4**3*t3*t5**2*t2**2+1140._ki*t6**3*t4**4*t&
                  &3*t5**2*t2+1695._ki*t6**3*t4**3*t3**2*t5*t2-70._ki*t3**4*t2*t6**2&
                  &*t4-45._ki*t3**3*t2*t6**3*t1+3._ki*t6**3*t1*t3*t2**5+3._ki*t6**3*t&
                  &1**2*t5*t2**4+t6**2*t1**2*t3*t2**4/2._ki-12._ki*t6*t5*t3**2*t2**6&
                  &+2._ki*t6**2*t4**4*t2**6-130._ki*t6**2*t4**3*t5*t2**4*t3+stemp6-2&
                  &._ki*t6*t1*t4*t5*t2**7-990._ki*t6**3*t4**2*t3**2*t5*t2**2-1500._ki&
                  &*t6**2*t4**3*t3**2*t5**2*t2+120._ki*t6*t1*t4*t5**3*t3*t2**3+14._k&
                  &i*t6*t4**3*t5*t3*t2**5+780._ki*t6**2*t4**3*t3*t5**2*t2**3+t6*t3*&
                  &*2*t2**7+930._ki*t6**2*t4**2*t3**2*t5*t2**3-330._ki*t6**2*t4*t3**&
                  &2*t5*t2**4-16._ki*t6*t4**2*t5*t3*t2**6+420._ki*t6*t4*t3**2*t5**2*&
                  &t2**4+t4**3*t2**9/4._ki-40._ki*t2*t3*t1**2*t6**2*t5**3-30._ki*t6**&
                  &2*t1**2*t4*t5**4*t2+5._ki*t3*t2**3*t1**2*t6**2*t5-480._ki*t6**2*t&
                  &4**4*t3*t5**2*t2**2-30._ki*t6**2*t1**2*t4*t5**3*t2**2+15._ki/2._ki&
                  &*t6**2*t1**2*t4*t5**2*t2**3
                !
                stemp6=-35._ki*t6**2*t1*t4**3*t5**2*t2**3-6._ki*t6*t1*t2**5*t3*t4*t&
                  &5-225._ki*t3*t1**2*t6**3*t4*t5**2-160._ki*t6**2*t1*t4**2*t5**4*t2&
                  &**2-7._ki*t6**2*t1*t4**3*t5*t2**4+3._ki*t6*t1*t4**2*t2**6*t5+140.&
                  &_ki*t6**2*t1*t4**3*t5**4*t2+80._ki*t6**2*t4*t5**3*t2**4*t3+stemp4&
                  &+40._ki*t6**2*t1*t4*t5**3*t2**4+3._ki*t6*t1*t2**6*t4*t3+12._ki*t3*&
                  &*2*t2**4*t1*t6*t5+135._ki*t2**2*t3**2*t1*t6**3*t5+80._ki*t6**2*t4&
                  &**4*t5*t2**3*t3+80._ki*t6**2*t5*t2*t1*t3**3+8._ki*t6**2*t1*t4**2*&
                  &t5*t2**5+40._ki*t6**2*t1*t4**2*t5**2*t2**4-390._ki*t6**2*t4**2*t3&
                  &*t5**2*t2**4-13._ki*t6**2*t4**3*t3*t2**5-50._ki*t6**2*t4**3*t5**4&
                  &*t2**3+80._ki*t6**2*t4**4*t5**4*t2**2-40._ki*t6**2*t4**5*t5**4*t2&
                  &+125._ki/2._ki*t6**2*t4**3*t3**2*t2**3+8._ki*t6**2*t4**4*t3*t2**4-&
                  &50._ki*t6**2*t4**3*t5**3*t2**4-120._ki*t3**2*t2**3*t6**2*t5*t1+75&
                  &._ki*t6*t1*t4*t5**4*t2**2*t3+stemp5-2._ki*t6**2*t1*t2**6*t4*t5-16&
                  &0._ki*t6**2*t1*t4**2*t5**3*t2**3+140._ki*t6**2*t1*t4**3*t5**3*t2*&
                  &*2
                !
                stemp3=stemp6+1040._ki*t6**2*t4**3*t3*t5**3*t2**2+3._ki/2._ki*t6**2*&
                  &t1**2*t4*t5*t2**4+300._ki*t3**2*t2**2*t6**2*t4*t5*t1-200._ki*t3**&
                  &3*t2**3*t6*t4*t5-18._ki*t6**3*t4**2*t3*t2**5+36._ki*t6**3*t4**5*t&
                  &3*t2**2+100._ki*t6**3*t4**5*t5**3*t2+60._ki*t6**2*t1*t3*t5**2*t2*&
                  &*4-90._ki*t6**2*t1*t4**2*t5*t2**3*t3-420._ki*t6**2*t1*t4*t3*t5**2&
                  &*t2**3-720._ki*t2*t3**2*t1*t6**3*t5*t4+1860._ki*t6**2*t4**2*t3**2&
                  &*t5**2*t2**2+600._ki*t2*t3**2*t1*t6**2*t4*t5**2-50._ki*t6*t4*t5**&
                  &4*t2**4*t3-5._ki/2._ki*t6**2*t1**2*t5**2*t2**4+10._ki*t6**2*t1**2*&
                  &t5**3*t2**3+10._ki*t6**2*t1**2*t5**4*t2**2-5._ki*t6**3*t4**2*t5**&
                  &3*t2**4+t6**3*t4**2*t5*t2**6+160._ki*t3**3*t2**3*t6**2*t4-760._ki&
                  &*t2*t3**3*t6**2*t5*t4**2+640._ki*t2**2*t3**3*t6**2*t5*t4+225._ki*&
                  &t6**3*t3**2*t2**3*t4*t5-t6**2*t2**6*t1*t3+2._ki*t6*t4**4*t5**2*t&
                  &2**5+4._ki*t6*t4*t3*t2**7*t5-80._ki*t6*t4*t5**3*t2**5*t3-25._ki*t3&
                  &**2*t2**3*t1*t6**2*t4+15._ki*t3*t2**2*t6**3*t4*t1**2+495._ki*t6**&
                  &3*t1*t4*t3*t5**2*t2**2-12._ki*t6*t1*t3*t5**2*t2**5+t6**2*t4**2*t&
                  &2**8/4._ki
                !
                stemp4=1._ki/t2**12*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=(-t4+t2)*t4**2/t2**3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/4._k&
                  &i
                !
                stemp11=-5._ki/3._ki*t6*t4*t2**7*t3**2*t1**4-2._ki*t6*t1**6*t4**2*t2&
                  &**6-40._ki/3._ki*t6*t1**4*t3**4*t2**4+5._ki*t6*t1**6*t3**2*t2**4+1&
                  &75._ki/2._ki*t6**2*t1**6*t3**3*t4-1155._ki/2._ki*t6**2*t1**5*t4**3*&
                  &t3**3-77._ki/6._ki*t6**2*t1**5*t4**3*t5*t2**5-418._ki/3._ki*t6**2*t&
                  &1**5*t4**4*t3*t2**3-520._ki/3._ki*t1**2*t4**5*t3**4*t2**4+70._ki*t&
                  &1**2*t2**2*t3**5*t4**5-95._ki/2._ki*t1**2*t2**5*t3**5*t4**2-460._k&
                  &i/3._ki*t1**2*t2**3*t3**5*t4**4+95._ki/12._ki*t1**2*t4*t3**5*t2**6&
                  &+755._ki/6._ki*t1**2*t2**4*t3**5*t4**3+155._ki/4._ki*t1*t2**6*t3**5&
                  &*t4**3-1300._ki/3._ki*t1**2*t4**5*t3**4*t5*t2**3
                !
                stemp10=stemp11-320._ki/3._ki*t1**2*t4**4*t3**3*t5*t2**6+825._ki/2._k&
                  &i*t6**2*t1**6*t4*t3*t5**2*t2**2+200._ki/3._ki*t1**3*t4**6*t3**2*t&
                  &5**3*t2**2+150._ki*t1**3*t4**4*t3**2*t5**3*t2**4-200._ki*t6*t1**4&
                  &*t4*t3**3*t5*t2**4-250._ki*t1**5*t4*t3**3*t5**2*t2**2-625._ki/3._k&
                  &i*t1**2*t4**3*t3**4*t5*t2**5-375._ki*t6**2*t1**4*t4**3*t5*t3**2*&
                  &t2**3-200._ki*t1**2*t4**6*t3**3*t5**2*t2**3+120._ki*t6**2*t1**3*t&
                  &4**6*t3**2*t5*t2**2+16._ki/3._ki*t6**2*t1**2*t4**9*t3**3-3._ki/2._k&
                  &i*t6*t1**5*t4**5*t2**5+3._ki*t6*t1**5*t4**4*t2**6-15._ki/8._ki*t6*&
                  &t1**5*t4**3*t2**7+80._ki*t1*t4**6*t3**4*t2**5+44._ki*t1*t2**2*t3*&
                  &*5*t4**7-7._ki*t1*t2**7*t3**5*t4**2
                !
                stemp11=16._ki*t1*t4**8*t3**4*t2**3-2325._ki/2._ki*t6**2*t1**6*t4**2&
                  &*t3*t5**2*t2+220._ki*t6*t1**4*t4**4*t3**3*t2**2+13._ki/6._ki*t6*t1&
                  &*t4**3*t3**4*t2**7+280._ki*t6*t1**3*t4**5*t3**2*t5**2*t2**3+40._k&
                  &i/3._ki*t6*t4**2*t5**3*t2**5*t1**4*t3+95._ki/4._ki*t6**2*t1**4*t3*&
                  &*3*t2**4*t4-465._ki/4._ki*t6*t1**5*t4**2*t3**2*t2**4+75._ki/4._ki*t&
                  &6*t1**5*t4**3*t5**2*t2**5-156._ki*t1**4*t4**3*t5*t2**5*t3**2+96.&
                  &_ki*t1**4*t4**4*t5*t2**4*t3**2+120._ki*t6*t1**5*t4**4*t5**4*t2**2&
                  &+240._ki*t6*t1**5*t4*t3**3*t2**3+100._ki*t1**4*t4**3*t5**3*t3*t2*&
                  &*5-450._ki*t1**5*t4**2*t5**3*t3**2*t2**2-160._ki*t6*t1**2*t4**5*t&
                  &3**3*t5*t2**4-448._ki/3._ki*t6*t1**2*t4**7*t3**3*t5*t2**2
                !
                stemp9=stemp11-820._ki/3._ki*t6*t1**4*t3**4*t2**2*t4**2-160._ki*t1**&
                  &4*t4**4*t5**3*t3*t2**4+35._ki*t6*t1**6*t4*t5*t2**4*t3+2790._ki*t6&
                  &*t1**5*t4**2*t3**2*t5**2*t2**2-30._ki*t6*t1**4*t4**5*t3**2*t2**3&
                  &-95._ki/2._ki*t6*t1**4*t4**3*t3**2*t2**5+t4**2*t3**4*t2**9*t1/2._k&
                  &i+t6*t1**6*t4*t2**7/2._ki-t6*t1**6*t3*t2**6/2._ki+100._ki*t1**5*t3&
                  &**3*t5**2*t2**3+8._ki*t1**5*t4**2*t3*t2**7-7._ki*t1**5*t4*t3**2*t&
                  &2**6+30._ki*t1**5*t4**3*t5**4*t2**4-3._ki*t1**5*t4**3*t5**2*t2**6&
                  &-70._ki*t1**3*t4**3*t3**2*t5**2*t2**6+570._ki*t6*t1**4*t4**3*t3**&
                  &2*t5*t2**4+stemp10
                !
                stemp11=180._ki*t6*t1**4*t4**4*t3*t5**2*t2**4+270._ki*t6*t1**6*t4**&
                  &2*t3*t5**2*t2**2+6215._ki/2._ki*t6**2*t1**5*t4**3*t3**2*t5*t2+30.&
                  &_ki*t1**4*t4**3*t3*t2**6*t5**2-4._ki/3._ki*t6*t1**4*t4**6*t3*t2**4&
                  &-75._ki/2._ki*t6**2*t1**6*t3*t5**2*t2**3-t6*t4**2*t3**4*t2**8*t1/&
                  &6._ki+78._ki*t1**4*t4**2*t5*t2**6*t3**2-2250._ki*t6*t1**5*t4**3*t3&
                  &**2*t5**2*t2-50._ki*t1**5*t4*t5**4*t2**4*t3+24._ki*t1**4*t4**5*t3&
                  &*t2**4*t5**2+220._ki*t1**4*t4*t3**3*t5*t2**5-25._ki*t6**2*t1**6*t&
                  &4*t5**3*t2**3+11._ki/6._ki*t6**2*t4**2*t2**6*t1**5*t5+t1**4*t4**2&
                  &*t3*t2**8*t5+40._ki*t6*t1**6*t3*t5**3*t2**3+140._ki*t6*t1**3*t4**&
                  &5*t3**4*t2
                !
                stemp10=stemp11-15._ki/4._ki*t6*t1**5*t3**2*t2**6-5._ki/6._ki*t6*t3**&
                  &4*t2**6*t1**3+10._ki/3._ki*t6*t3**3*t2**6*t1**4+5._ki/2._ki*t6**2*t&
                  &1**6*t3*t2**5-565._ki/3._ki*t1**3*t4**3*t3**4*t2**4-4._ki/3._ki*t1*&
                  &*3*t4**6*t3**2*t2**5+110._ki*t1**3*t4**4*t3**4*t2**3-12._ki*t1**3&
                  &*t4**2*t3**3*t2**7+105._ki/2._ki*t1**3*t3**5*t2**2*t4**3+7._ki/6._k&
                  &i*t1**3*t4**3*t3**2*t2**8-8._ki*t1**4*t4**4*t3**2*t2**5-13._ki/2.&
                  &_ki*t1**4*t4**2*t3**2*t2**7+22._ki*t1**4*t4*t3**3*t2**6+3._ki*t1**&
                  &5*t4**3*t5*t2**7+170._ki*t1*t2**4*t3**5*t4**5-110._ki*t1*t2**5*t3&
                  &**5*t4**4-56._ki*t1*t4**7*t3**4*t2**4
                !
                stemp11=25._ki*t1*t4**4*t3**4*t2**7+4._ki*t1**5*t4*t3*t2**7*t5-33._k&
                  &i/4._ki*t6**2*t1**3*t4**3*t5*t2**5*t3**2-120._ki*t1**3*t4**2*t3**&
                  &3*t2**6*t5-180._ki*t6*t1**5*t3**3*t5*t2**3+240._ki*t6*t1**4*t4**4&
                  &*t5**3*t3*t2**3-60._ki*t1*t4**5*t3**4*t2**6-11._ki/2._ki*t1*t4**3*&
                  &t3**4*t2**8-136._ki*t1*t2**3*t3**5*t4**6-4._ki*t1**4*t4**4*t3*t2*&
                  &*7+5._ki/2._ki*t1**4*t4**3*t3*t2**8+35._ki/2._ki*t1**4*t3**5*t2**2*&
                  &t4+2._ki*t1**4*t4**5*t3*t2**6+13._ki*t1**4*t4**3*t3**2*t2**6-3._ki&
                  &*t1**3*t4**4*t3**2*t2**7-24._ki*t1**5*t4*t5**2*t2**6*t3-60._ki*t1&
                  &**4*t4*t3**2*t5**2*t2**6
                !
                stemp8=stemp11+500._ki*t1**4*t4**3*t3**3*t5*t2**3-2825._ki/6._ki*t1*&
                  &*3*t4**3*t3**4*t2**3*t5+50._ki/3._ki*t6*t1**2*t4**4*t3**3*t2**6-5&
                  &._ki*t6*t1**6*t3*t5*t2**5+95._ki/6._ki*t6*t1**3*t3**4*t4*t2**5+24.&
                  &_ki*t6**2*t1**3*t4**8*t3**2*t5-80._ki*t1**5*t4*t5**3*t2**5*t3-200&
                  &._ki*t1**5*t4*t5*t3**3*t2**3-108._ki*t1**5*t4**2*t3**2*t5*t2**4+1&
                  &00._ki/3._ki*t6*t1**4*t4**5*t3*t2**4*t5-30._ki*t6*t1**4*t4**4*t3*t&
                  &2**5*t5+375._ki/4._ki*t6*t1**5*t4**3*t3**2*t2**3+180._ki*t1**3*t4*&
                  &*4*t3**2*t5**2*t2**5-200._ki*t1**3*t4**5*t3**2*t5**2*t2**4+t1**5&
                  &*t3**2*t2**7+72._ki*t6*t1*t2**3*t3**4*t4**7+stemp10+stemp9
                !
                stemp11=128._ki/3._ki*t6*t1**2*t4**8*t3**3*t5*t2+170._ki*t6**2*t1**3&
                  &*t4**5*t3**3*t2**2+5._ki*t6**2*t1**6*t4*t5*t2**5-55._ki/2._ki*t6**&
                  &2*t5*t2**4*t3**2*t1**5-35._ki/3._ki*t6*t1**3*t4**5*t3**2*t2**5-10&
                  &._ki/3._ki*t6*t1**3*t4**7*t3**2*t2**3+35._ki/3._ki*t6*t1**4*t4**3*t&
                  &3*t2**6*t5-2080._ki/3._ki*t6*t1**3*t4**5*t3**3*t5*t2**2+640._ki/3.&
                  &_ki*t6*t1**2*t4**6*t3**3*t5*t2**3+350._ki*t1**5*t4*t5**3*t3**2*t2&
                  &**3-200._ki/3._ki*t2**5*t3**5*t4**6+5._ki/3._ki*t3**4*t2**7*t1**3+t&
                  &6**2*t4*t3**3*t2**6*t1**3/2._ki-360._ki*t6*t1**4*t4**2*t3**2*t5**&
                  &2*t2**4-40._ki/3._ki*t6*t1**4*t4**6*t3*t2**3*t5+40._ki/3._ki*t4*t3*&
                  &*3*t2**7*t1**3*t5
                !
                stemp10=stemp11-60._ki*t6*t1**6*t2**3*t3**2*t5-5._ki/2._ki*t6*t4**2*&
                  &t5*t3**2*t2**7*t1**3+80._ki*t6*t1**4*t4**6*t3*t5**2*t2**2+275._ki&
                  &*t1**3*t4**4*t3**4*t2**2*t5-14._ki*t1**3*t3**2*t4**3*t5*t2**7-78&
                  &0._ki*t6**2*t1**4*t4**5*t3**2*t5*t2-112._ki/3._ki*t6*t1**2*t4**7*t&
                  &3**3*t2**3+200._ki/3._ki*t6*t1**2*t4**4*t5*t2**5*t3**3+110._ki/3._k&
                  &i*t6*t1*t4**5*t3**4*t2**5-8._ki*t1**5*t4**4*t5**3*t2**4+4._ki/3._k&
                  &i*t4*t2**8*t3**3*t1**3+95._ki*t1**4*t4**2*t3**4*t2**3-50._ki*t1**&
                  &5*t3**2*t5**3*t2**4-2._ki*t1**5*t4*t3*t2**8-10._ki*t1**5*t4**2*t5&
                  &**4*t2**5-210._ki*t6*t1**6*t4*t3*t5**2*t2**3+15._ki*t6*t1**4*t4**&
                  &2*t3**2*t2**6
                !
                stemp11=1170._ki*t6*t1**5*t4**3*t3*t5**2*t2**3-1900._ki/3._ki*t1**3*&
                  &t4**4*t3**3*t5**2*t2**3+300._ki*t1**3*t4**5*t3**3*t5**2*t2**2-15&
                  &._ki/2._ki*t1**4*t3**5*t2**3+20._ki/3._ki*t6*t1**3*t4**4*t3**2*t2**&
                  &6-15._ki/8._ki*t6*t1**3*t4**3*t3**2*t2**7+36._ki*t6**2*t1**2*t4**7&
                  &*t3**3*t2**2-5._ki*t6*t4**2*t5**2*t2**6*t1**3*t3**2+13._ki/12._ki*&
                  &t6**2*t1**2*t4**3*t3**3*t2**6-6._ki*t6**2*t1**2*t4**4*t3**3*t2**&
                  &5-544._ki/3._ki*t6*t1**2*t4**6*t3**4*t2**2-25._ki/2._ki*t6*t1**6*t2&
                  &**3*t4*t3**2+385._ki/6._ki*t6**2*t1**5*t4**3*t5**3*t2**3+84._ki*t1&
                  &**5*t4*t5*t2**5*t3**2+1140._ki*t6*t1**4*t4**3*t3**2*t5**2*t2**3+&
                  &16._ki*t1**3*t3**2*t4**6*t5*t2**4-25._ki/6._ki*t4*t5*t2**7*t3**4*t&
                  &1**2
                !
                stemp9=stemp11+80._ki*t1**3*t4**6*t3**2*t5**2*t2**3-720._ki*t6*t1**&
                  &5*t4**4*t3*t5**2*t2**2+780._ki*t6**2*t1**4*t4**4*t3**2*t5*t2**2-&
                  &39._ki/2._ki*t6*t1**5*t4**3*t3*t2**5-75._ki*t6*t1**5*t4**3*t5**3*t&
                  &2**4-35._ki/2._ki*t6*t1**6*t4**3*t5**2*t2**3+70._ki*t6*t1**6*t4**3&
                  &*t5**3*t2**2+200._ki*t1*t4**6*t3**4*t5*t2**4+80._ki*t1**5*t3**3*t&
                  &5*t2**4-2._ki*t1**5*t4**4*t5*t2**6-25._ki*t1**5*t3**4*t5*t2**2+15&
                  &._ki*t1**5*t4**3*t5**5*t2**3-10._ki*t1**5*t4**4*t5**5*t2**2-75._ki&
                  &/2._ki*t6**2*t1**6*t3**3*t2-45._ki*t6**2*t1**7*t4**2*t5**3-2._ki*t&
                  &6**2*t1**7*t3*t2**3-20._ki*t1**5*t4**4*t5**4*t2**3+stemp10
                !
                stemp11=325._ki/2._ki*t6**2*t1**6*t4**2*t5**3*t2**2+800._ki/3._ki*t6*&
                  &t1**3*t4**6*t3**3*t5*t2-75._ki*t6*t1**5*t4**3*t5**4*t2**3+40._ki/&
                  &3._ki*t6*t5*t2**5*t3**3*t1**4-10._ki/3._ki*t1**3*t3**5*t2**5+90._ki&
                  &*t6*t1**5*t3**2*t2**4*t5**2+15._ki*t6*t1**5*t4**2*t5**4*t2**4-12&
                  &5._ki/2._ki*t1**3*t4*t3**4*t2**5*t5+10._ki/3._ki*t6*t1**4*t4**5*t3*&
                  &t2**5+40._ki*t6*t1**3*t4**7*t3**2*t5*t2**2+200._ki/3._ki*t1**2*t4*&
                  &*7*t3**3*t5**2*t2**2-620._ki*t1**4*t4**2*t3**3*t5*t2**4+560._ki/3&
                  &._ki*t1**2*t4**5*t3**3*t5*t2**5-32._ki/3._ki*t1**2*t4**4*t3**3*t2*&
                  &*7-80._ki*t6*t1**6*t4**2*t5**4*t2**2-285._ki*t6*t1**5*t3**3*t2**2&
                  &*t4**2+320._ki*t1**5*t4**2*t3*t5**3*t2**4
                !
                stemp10=stemp11-5._ki/3._ki*t6*t4*t3**3*t2**7*t1**3-75._ki*t6**2*t1*&
                  &*7*t4*t3*t5**2+20._ki*t6*t1**6*t4*t5**3*t2**4+40._ki*t6*t4*t3**2*&
                  &t2**5*t1**4*t5**2+825._ki/2._ki*t6**2*t1**5*t3**2*t2**3*t4*t5-30.&
                  &_ki*t6*t1**5*t4**4*t5**2*t2**4+14._ki*t1**5*t4**3*t5*t3*t2**5+110&
                  &._ki*t6*t1**4*t3**4*t2**3*t4-105._ki*t6*t1**5*t3**4*t2*t4-128._ki/&
                  &3._ki*t6*t1*t2**2*t3**4*t4**8-44._ki/3._ki*t6*t1**2*t4**3*t5*t2**6&
                  &*t3**3-12._ki*t4**4*t3**5*t2**7-t6*t1**6*t2**6*t4*t5+30._ki*t6*t1&
                  &**6*t3*t5**2*t2**4-40._ki*t6*t1**2*t4**5*t3**3*t2**5-25._ki/6._ki*&
                  &t4**2*t5**2*t2**7*t1**2*t3**3+480._ki*t1**4*t4**4*t3**2*t5**2*t2&
                  &**3
                !
                stemp7=-16._ki*t1**5*t4**2*t5*t3*t2**6+stemp8-250._ki/3._ki*t6*t1**3&
                  &*t4**3*t3**3*t2**5-325._ki*t6**2*t1**6*t4**3*t5**3*t2-125._ki/2._k&
                  &i*t6**2*t1**6*t2**2*t3*t4**3-5._ki*t6*t1**6*t4*t5**2*t2**5-540._k&
                  &i*t1**5*t4**2*t3**2*t5**2*t2**3+495._ki*t6**2*t1**5*t4**2*t3*t5*&
                  &*2*t2**3+700._ki/3._ki*t1**2*t4**5*t3**3*t5**2*t2**4-400._ki/3._ki*&
                  &t1**2*t4**4*t3**3*t5**2*t2**5+stemp10-180._ki*t6*t1**4*t4**2*t3*&
                  &*2*t5*t2**5-200._ki*t6*t1**4*t4**5*t3*t5**2*t2**3+33._ki*t6**2*t1&
                  &**5*t4**4*t5*t2**4+110._ki/3._ki*t4**5*t3**5*t2**6-2._ki*t1**4*t3*&
                  &*3*t2**7+15._ki*t6*t1**5*t4**5*t5**2*t2**3-1000._ki/3._ki*t6*t1**3&
                  &*t4**3*t3**3*t5*t2**4+176._ki/3._ki*t6*t1**2*t4**7*t3**4*t2-7._ki*&
                  &t6**2*t1**3*t4**2*t3**3*t2**5+90._ki*t6*t1**5*t4*t3*t5**2*t2**5+&
                  &96._ki*t1**5*t4**2*t5**2*t2**5*t3-150._ki*t1**3*t4**2*t3**3*t5**2&
                  &*t2**5+20._ki*t6*t1**6*t4*t5**4*t2**3+36._ki*t1**3*t3**2*t4**4*t5&
                  &*t2**6-128._ki/3._ki*t2**3*t3**5*t4**8+72._ki*t2**4*t3**5*t4**7+15&
                  &._ki*t1**4*t3**4*t2**5-t4**2*t3**5*t2**9/6._ki+8._ki*t1**5*t3**3*t&
                  &2**5-175._ki/3._ki*t1**3*t4**3*t3**2*t5**3*t2**5-136._ki*t6**2*t1*&
                  &*3*t4**6*t3**3*t2+165._ki/4._ki*t6*t1**5*t4*t3**2*t2**5-55._ki/6._k&
                  &i*t6**2*t4**2*t5**3*t2**4*t1**5+stemp9
                !
                stemp10=-100._ki/3._ki*t6**2*t1**2*t4**6*t3**3*t2**3-64._ki/3._ki*t6*&
                  &*2*t1**2*t4**8*t3**3*t2-5._ki/3._ki*t6*t4**2*t2**7*t1**4*t3*t5-11&
                  &40._ki*t6*t1**5*t4**2*t3**3*t5*t2-12._ki*t1**4*t3**2*t5*t4*t2**7-&
                  &50._ki*t1**4*t4*t5**3*t2**5*t3**2-20._ki*t1**4*t4**2*t3*t5**3*t2*&
                  &*6+220._ki*t6*t1**4*t2**4*t3**3*t4**2-800._ki/3._ki*t6*t1**4*t4**5&
                  &*t5**3*t3*t2**2-50._ki*t6*t1**4*t2**5*t3**3*t4+190._ki/3._ki*t6*t1&
                  &**4*t4**4*t3**2*t2**4+10._ki*t6*t4**2*t3*t2**6*t1**4*t5**2-5._ki*&
                  &t6**2*t1**7*t5**3*t2**2+t4*t3**5*t2**8*t1/2._ki+25._ki/6._ki*t5*t2&
                  &**6*t3**4*t1**3-5._ki*t1**5*t4**2*t5**5*t2**4-t1**5*t4**2*t2**8*&
                  &t5-60._ki*t1**5*t3**2*t2**5*t5**2+t1**5*t4**2*t5**2*t2**7-10._ki*&
                  &t1**5*t2**3*t3**4-7._ki*t1**5*t4**3*t3*t2**6-20._ki*t1**5*t3**3*t&
                  &2**4*t4+10._ki*t6*t1**6*t3**3*t2**2-5._ki/3._ki*t4*t3**4*t2**8*t1*&
                  &*2+12._ki*t1**5*t4**3*t5**3*t2**5-80._ki*t1**4*t4*t3**4*t2**4+75.&
                  &_ki/2._ki*t1**4*t3**4*t5*t2**4+9._ki*t1**5*t4**2*t3**2*t2**5-62._ki&
                  &*t1**4*t4**2*t3**3*t2**5+360._ki*t6*t1**6*t4**2*t5**3*t2*t3-3._ki&
                  &/4._ki*t6*t1**5*t4**2*t5*t2**7-6._ki*t1**4*t4**2*t5**2*t2**7*t3+1&
                  &50._ki*t6*t1**6*t4*t5*t3**2*t2**2
                !
                stemp11=stemp10-25._ki/2._ki*t1**4*t4**2*t3*t5**4*t2**5+155._ki/4._ki&
                  &*t6**2*t1**3*t3**3*t4**3*t2**4-285._ki/2._ki*t6**2*t1**4*t4**2*t3&
                  &**3*t2**3+120._ki*t6*t1**5*t4**4*t5*t2**3*t3-14._ki*t6**2*t1**4*t&
                  &4**5*t3*t2**4+30._ki*t6**2*t1**7*t3*t5**2*t2-7._ki*t6**2*t1**7*t4&
                  &*t5*t2**3-40._ki*t1**3*t3**2*t4**5*t5*t2**5+475._ki*t1**3*t4**3*t&
                  &3**3*t5**2*t2**4+380._ki*t1**3*t4**3*t3**3*t2**5*t5+135._ki/4._ki*&
                  &t6**2*t1**4*t4**3*t5**2*t2**4*t3+55._ki/3._ki*t6**2*t1**2*t4**5*t&
                  &3**3*t2**4-80._ki*t6*t1**3*t4**4*t3**2*t5*t2**5-120._ki*t6*t1**3*&
                  &t4**6*t3**2*t5*t2**3-25._ki*t1**3*t4*t3**4*t2**6+55._ki/2._ki*t1**&
                  &3*t3**5*t2**4*t4
                !
                stemp9=stemp11+110._ki/3._ki*t6**2*t1**5*t2**3*t3**3-220._ki/3._ki*t6&
                  &**2*t1**5*t4**6*t5**3+210._ki*t6**2*t1**4*t4**5*t3**3+44._ki*t6**&
                  &2*t1**3*t4**7*t3**3-4._ki*t6**2*t1**4*t4**7*t3*t2**2+10._ki/3._ki*&
                  &t1**3*t4**5*t3**2*t2**6+110._ki*t1**3*t4**2*t3**4*t2**5+38._ki*t1&
                  &**3*t4**3*t3**3*t2**6-152._ki/3._ki*t1**3*t4**4*t3**3*t2**5+24._ki&
                  &*t1**3*t4**5*t3**3*t2**4+3._ki/8._ki*t6*t1**5*t4**2*t2**8+35._ki*t&
                  &6**2*t1**7*t4*t5**3*t2-165._ki*t6**2*t1**5*t4**4*t5**3*t2**2+25.&
                  &_ki/3._ki*t4**2*t5**3*t2**6*t1**3*t3**2+230._ki/3._ki*t6*t1**3*t4**&
                  &2*t3**3*t5*t2**5-160._ki*t6*t1**3*t4**4*t3**2*t5**2*t2**4+200._ki&
                  &*t1**5*t4**2*t5**4*t3*t2**3
                !
                stemp11=550._ki/3._ki*t6**2*t1**5*t4**5*t5**3*t2+680._ki/3._ki*t6*t1*&
                  &*2*t4**5*t3**4*t2**3-15._ki*t6*t1**5*t4*t5*t2**6*t3-28._ki/3._ki*t&
                  &6*t1**2*t4**2*t3**4*t2**6+880._ki*t6*t1**4*t4**4*t3**3*t5*t2-650&
                  &._ki*t1**4*t4**3*t3**2*t5**3*t2**3+300._ki*t6*t1**6*t4*t3**2*t5**&
                  &2*t2-90._ki*t6**2*t1**3*t4**5*t3**2*t5*t2**3-1125._ki*t6*t1**5*t4&
                  &**3*t3**2*t5*t2**2-495._ki*t6*t1**5*t4*t3**2*t5*t2**4-150._ki*t1*&
                  &t4**5*t3**4*t5*t2**5+40._ki*t1*t4**8*t3**4*t5*t2**2+80._ki*t1**4*&
                  &t4**5*t5**3*t3*t2**3-240._ki*t6*t1**3*t4**6*t3**2*t5**2*t2**2-84&
                  &._ki*t6**2*t1**3*t4**7*t3**2*t5*t2+44._ki/3._ki*t6**2*t1**5*t4**6*&
                  &t5*t2**2-110._ki/3._ki*t6**2*t1**5*t4**5*t5*t2**3
                !
                stemp10=stemp11+195._ki/2._ki*t6*t1**5*t4**2*t5*t2**5*t3-3135._ki/2.&
                  &_ki*t6**2*t1**5*t4**3*t3*t5**2*t2**2-5._ki/12._ki*t3**5*t2**7*t1**&
                  &2-3._ki*t6*t1**4*t4**4*t3*t2**6-440._ki/3._ki*t6*t1**2*t2**4*t4**4&
                  &*t3**4-200._ki/3._ki*t6*t1*t2**4*t3**4*t4**6-280._ki*t6*t1**6*t4*t&
                  &5**3*t3*t2**2+13._ki/6._ki*t4**3*t3**5*t2**8+66._ki*t6**2*t1**5*t4&
                  &**5*t3*t2**2+1560._ki*t6*t1**5*t4**3*t3*t5**3*t2**2+50._ki*t1**4*&
                  &t3**3*t4**3*t2**4+2._ki*t1**5*t4**4*t5**2*t2**5-205._ki/3._ki*t1**&
                  &3*t3**5*t2**3*t4**2+45._ki*t6*t1**5*t3**4*t2**2-45._ki*t6*t1**5*t&
                  &3**3*t2**4+7._ki/4._ki*t6*t1**6*t4**3*t2**5-15._ki*t6**2*t1**7*t3*&
                  &*2*t5
                !
                stemp11=stemp10+200._ki*t6**2*t1**6*t4**4*t5**3-3._ki/2._ki*t6*t1**5&
                  &*t4*t3*t2**7-280._ki*t1**5*t4**3*t3*t5**3*t2**3+420._ki*t1**5*t4*&
                  &t3**2*t5**2*t2**4+210._ki*t6*t1**4*t4**3*t3**4*t2+20._ki*t6*t4*t3&
                  &**2*t2**6*t1**4*t5-175._ki*t1**5*t4**3*t3*t5**4*t2**2-990._ki*t6*&
                  &t1**5*t4*t3**2*t5**2*t2**3+300._ki*t6**2*t1**4*t4**6*t3**2*t5+50&
                  &0._ki/3._ki*t1**2*t4**6*t3**4*t5*t2**2-6._ki*t6*t1**5*t4**4*t5*t2*&
                  &*5+15._ki/4._ki*t6*t1**5*t4**3*t5*t2**6-9._ki/4._ki*t6**2*t1**4*t4*&
                  &*3*t3*t2**6+15._ki*t6*t1**5*t4**2*t5**3*t2**5+50._ki*t1**4*t4**5*&
                  &t5**4*t3*t2**2+125._ki/2._ki*t1**4*t4**3*t5**4*t3*t2**4
                !
                stemp8=stemp11-11._ki/3._ki*t6*t1**2*t4**3*t3**3*t2**7+140._ki*t6*t1&
                  &**3*t4**5*t3**2*t5*t2**4-60._ki*t6*t1**5*t4**5*t5**3*t2**2+225._k&
                  &i/2._ki*t6**2*t1**6*t3**2*t5*t2**2+120._ki*t6*t1**5*t4**4*t5**3*t&
                  &2**3-70._ki*t6*t1**4*t4**3*t3*t5**2*t2**5-140._ki*t1*t4**7*t3**4*&
                  &t5*t2**3+45._ki*t6*t1**5*t2**5*t3**2*t5+5._ki*t6**2*t1**7*t2**2*t&
                  &3*t4+9._ki*t6**2*t1**7*t4**2*t5*t2**2+12._ki*t6**2*t1**4*t4**6*t3&
                  &*t2**3-1130._ki/3._ki*t6*t1**4*t2**3*t3**3*t4**3-60._ki*t6*t1**5*t&
                  &4**5*t5**4*t2+50._ki/3._ki*t4*t3**3*t2**6*t1**3*t5**2-t6**2*t4**2&
                  &*t3**3*t2**7*t1**2/12._ki+65._ki*t6**2*t1**6*t4**3*t5*t2**3+10._ki&
                  &*t4**2*t3**2*t2**7*t1**3*t5**2+stemp9
                !
                stemp11=stemp8+209._ki/2._ki*t6**2*t1**5*t4**3*t3*t2**4+2090._ki*t6*&
                  &*2*t1**5*t4**4*t3*t5**2*t2-120._ki*t6**2*t1**4*t4**4*t5**2*t2**3&
                  &*t3+2._ki*t4**2*t3**2*t2**8*t1**3*t5+320._ki/3._ki*t6*t1**4*t4**6*&
                  &t5**3*t3*t2+210._ki*t6**2*t1**4*t4**5*t3*t5**2*t2**2+160._ki/3._ki&
                  &*t1**2*t4**7*t3**3*t5*t2**3+880._ki*t6*t1**4*t4**2*t3**3*t2**3*t&
                  &5-33._ki*t6**2*t1**5*t4**2*t3*t2**5-780._ki*t1**4*t4**3*t3**2*t5*&
                  &*2*t2**4-1815._ki*t6**2*t1**5*t4**4*t3**2*t5+115._ki/6._ki*t6*t1**&
                  &3*t4**2*t3**3*t2**6-95._ki*t6*t1**3*t4**2*t3**4*t2**4-180._ki*t6*&
                  &*2*t1**4*t4**6*t3*t5**2*t2+8._ki*t6**2*t1**4*t4**4*t3*t2**5
                !
                stemp10=stemp11+8._ki*t1**4*t4**4*t3*t2**6*t5-48._ki*t1**4*t4**4*t3&
                  &*t2**5*t5**2-t6*t1**4*t4**2*t3*t2**8/6._ki-280._ki/3._ki*t6*t1**4*&
                  &t4**3*t5**3*t3*t2**4+60._ki*t6**2*t1**4*t4**7*t3*t5**2-55._ki*t6*&
                  &*2*t4*t3*t2**4*t1**5*t5**2-9._ki/2._ki*t6*t1**6*t4**2*t3*t2**4-7.&
                  &_ki/2._ki*t6*t1**6*t4**3*t5*t2**4+1395._ki*t6*t1**5*t4**2*t3**2*t5&
                  &*t2**3+30._ki*t1**2*t4**3*t5*t3**3*t2**7+625._ki*t1**4*t4**3*t3**&
                  &3*t5**2*t2**2-20._ki/3._ki*t6*t4*t5*t2**6*t3**3*t1**3+4._ki/3._ki*t&
                  &6*t4**2*t5*t2**7*t1**2*t3**3-1815._ki*t6**2*t1**5*t4**2*t3**2*t5&
                  &*t2**2+345._ki/4._ki*t6**2*t1**4*t4**2*t3**2*t5*t2**4+t6**2*t4**2&
                  &*t3*t2**7*t1**4/4._ki-585._ki*t6*t1**5*t4**2*t3*t5**2*t2**4
                !
                stemp11=stemp10+960._ki*t6*t1**5*t3**3*t5*t2**2*t4+32._ki/3._ki*t2**&
                  &2*t3**5*t4**9+2255._ki/3._ki*t6**2*t1**5*t3**3*t4**2*t2-460._ki*t6&
                  &**2*t1**4*t4**4*t3**3*t2-12._ki*t6*t1*t4**4*t3**4*t2**6-10._ki/3.&
                  &_ki*t4**2*t5*t3**3*t2**8*t1**2-775._ki*t1**4*t4**2*t3**3*t2**3*t5&
                  &**2-1520._ki*t6*t1**4*t4**4*t3**2*t5**2*t2**2+575._ki/12._ki*t1**2&
                  &*t4**2*t3**4*t5*t2**6-200._ki*t1**4*t3**4*t5*t2**3*t4+1425._ki/2.&
                  &_ki*t6**2*t1**6*t4**2*t3**2*t5+155._ki/2._ki*t6**2*t1**6*t2**3*t3*&
                  &t4**2+720._ki*t6*t1**4*t4**5*t3**2*t5**2*t2-500._ki/3._ki*t1**3*t4&
                  &**5*t3**2*t5**3*t2**3-15._ki/4._ki*t6*t1**5*t4**2*t5**2*t2**6+2._k&
                  &i/3._ki*t6*t4*t3**4*t2**7*t1**2
                !
                stemp9=stemp11+120._ki*t6*t1**5*t4*t5**3*t2**4*t3-600._ki*t6**2*t1*&
                  &*6*t3**2*t5*t2*t4+11._ki/3._ki*t6**2*t4*t2**6*t3*t1**5-15._ki/4._ki&
                  &*t6**2*t4**2*t5**2*t2**5*t1**4*t3-55._ki/4._ki*t1*t4**3*t5*t2**7*&
                  &t3**4+1300._ki/3._ki*t1**2*t4**4*t3**4*t5*t2**4+10._ki*t6*t1**3*t4&
                  &**6*t3**2*t2**4+stemp7+5._ki/24._ki*t6*t4**2*t3**2*t2**8*t1**3+3.&
                  &_ki/4._ki*t6**2*t4**2*t5*t2**6*t1**3*t3**2+75._ki/2._ki*t1**2*t4**3&
                  &*t5**2*t2**6*t3**3-15._ki/2._ki*t6**2*t4*t5*t2**5*t3**2*t1**4-195&
                  &._ki*t6*t1**5*t4**3*t5*t2**4*t3-920._ki/3._ki*t6*t1**3*t4**4*t3**4&
                  &*t2**2-5._ki/4._ki*t6**2*t3**3*t2**5*t1**4+3._ki*t6*t1**5*t4**5*t5&
                  &*t2**4-760._ki*t6*t1**4*t4**4*t3**2*t5*t2**3+1875._ki/2._ki*t6**2*&
                  &t1**6*t4**3*t3*t5**2
                !
                stemp11=stemp9-100._ki*t1**4*t4**4*t5**4*t3*t2**3+400._ki*t1**4*t4*&
                  &*4*t3**2*t5**3*t2**2+275._ki*t1**4*t4*t3**3*t2**4*t5**2+40._ki*t6&
                  &*t1**6*t3**3*t5*t2+4._ki*t6*t1**6*t4**2*t5*t2**5+32._ki/3._ki*t6*t&
                  &1*t2*t3**4*t4**9-4520._ki/3._ki*t6*t1**4*t4**3*t3**3*t2**2*t5+5._k&
                  &i/4._ki*t4**2*t5*t2**8*t1*t3**4+125._ki/2._ki*t1*t4**4*t5*t2**6*t3&
                  &**4+160._ki/3._ki*t6*t1**2*t4**6*t3**3*t2**4+32._ki/3._ki*t6*t1**2*&
                  &t4**8*t3**3*t2**2-65._ki/2._ki*t6**2*t1**6*t4**2*t5*t2**4+45._ki/2&
                  &._ki*t6*t1**3*t4**3*t5*t3**2*t2**6+360._ki*t6*t1**4*t4**5*t3**2*t&
                  &5*t2**2+20._ki*t6*t1**6*t4**2*t5**2*t2**4+7._ki/2._ki*t6*t1**6*t4*&
                  &t3*t2**5
                !
                stemp10=stemp11+80._ki*t6*t1**3*t4**7*t3**2*t5**2*t2-12._ki*t1**5*t&
                  &5*t3**2*t2**6-4._ki*t1**5*t4**2*t5**3*t2**6-25._ki*t1**4*t3**3*t5&
                  &**2*t2**5+t6**2*t1**7*t5*t2**4-t4**2*t3**3*t2**9*t1**2/3._ki-20.&
                  &_ki*t1**4*t3**3*t5*t2**6+t1**4*t4*t3**2*t2**8-t1**4*t4**2*t3*t2*&
                  &*9/2._ki-t1**3*t4**2*t3**2*t2**9/6._ki-55._ki/2._ki*t6**2*t1**6*t2*&
                  &*4*t3*t4-40._ki*t6**2*t1**6*t4**4*t5*t2**2-80._ki*t6*t1**6*t4**2*&
                  &t5**3*t2**3-780._ki*t6*t1**5*t4**2*t5**3*t3*t2**3+755._ki/3._ki*t6&
                  &*t1**3*t4**3*t3**4*t2**3-110._ki*t6**2*t1**3*t4**4*t3**3*t2**3-1&
                  &20._ki*t6*t1**6*t3**2*t5**2*t2**2
                !
                stemp11=stemp10+200._ki/3._ki*t6*t1**3*t4**6*t3**3*t2**2+520._ki/3._k&
                  &i*t6*t1**3*t4**4*t3**3*t2**4-520._ki/3._ki*t6*t1**3*t4**5*t3**3*t&
                  &2**3-960._ki*t6*t1**5*t4**4*t3*t5**3*t2+t6*t4**2*t3**3*t2**8*t1*&
                  &*2/3._ki+325._ki*t1**4*t4**2*t3**2*t5**3*t2**4-84._ki*t1**5*t4**3*&
                  &t5**2*t2**4*t3+45._ki*t6*t1**3*t4**3*t5**2*t2**5*t3**2+115._ki/6.&
                  &_ki*t1**2*t4**2*t3**4*t2**7-250._ki/3._ki*t1**2*t4**3*t3**4*t2**6+&
                  &200._ki/3._ki*t1**2*t4**6*t3**4*t2**3+3._ki*t1**2*t4**3*t3**3*t2**&
                  &8+16._ki/3._ki*t1**2*t4**7*t3**3*t2**4-45._ki*t6*t1**6*t4**2*t5*t2&
                  &**3*t3+390._ki*t1**4*t4**2*t3**2*t5**2*t2**5-605._ki/2._ki*t6**2*t&
                  &1**5*t3**3*t2**2*t4+39._ki/4._ki*t6*t1**5*t4**2*t3*t2**6
                !
                stemp6=stemp11-160._ki*t1**2*t4**6*t3**3*t5*t2**4+7._ki/6._ki*t6*t1*&
                  &*4*t4**3*t3*t2**7-16._ki*t1**2*t4**6*t3**3*t2**5+56._ki/3._ki*t1**&
                  &2*t4**5*t3**3*t2**6+520._ki/3._ki*t1**2*t4**4*t3**4*t2**5-5._ki*t1&
                  &**4*t4**3*t3*t2**7*t5-4._ki*t1**4*t4**5*t3*t2**5*t5+275._ki*t1**3&
                  &*t4**2*t3**4*t2**4*t5-1520._ki/3._ki*t1**3*t4**4*t3**3*t2**4*t5+7&
                  &55._ki/2._ki*t6**2*t1**4*t4**3*t3**3*t2**2+2080._ki/3._ki*t6*t1**3*&
                  &t4**4*t3**3*t5*t2**3+75._ki/2._ki*t6**2*t1**3*t4**4*t5*t2**4*t3**&
                  &2+70._ki*t6*t1**6*t4**3*t5**4*t2+240._ki*t1**3*t4**5*t3**3*t2**3*&
                  &t5+475._ki/2._ki*t1**4*t4**2*t3**4*t5*t2**2+155._ki/3._ki*t6*t1**2*&
                  &t3**4*t4**3*t2**5-990._ki*t6**2*t1**5*t4**5*t3*t5**2+12._ki*t6*t1&
                  &**5*t4**4*t3*t2**4
                !
                stemp7=t6/t1**5/t2**12
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(4)
                !
                stemp5=20._ki*t6**2*t1*t2*t3*t4*t5**2+5._ki*t6**2*t1*t2*t5*t3**2+5.&
                  &_ki*t6**3*t1*t3*t5*t2**2+65._ki*t6**3*t1*t4**2*t3*t5-25._ki*t6**3*&
                  &t1*t4**2*t5**2*t2-14._ki/3._ki*t6**2*t1*t2**3*t3*t4+16._ki*t6**3*t&
                  &1*t2**2*t3*t4+7._ki/6._ki*t6**2*t4*t3*t2**5+10._ki/3._ki*t6**2*t4**&
                  &3*t5**3*t2**2-10._ki/3._ki*t6**2*t4**4*t5**3*t2-5._ki/6._ki*t6**2*t&
                  &4**2*t5**3*t2**3+5._ki*t6**2*t1*t4**2*t2**4-2._ki*t6**2*t1*t4*t2*&
                  &*5+t6*t4*t2**6*t1/2._ki+20._ki*t6**3*t1*t4**2*t5*t2**2+56._ki/3._ki&
                  &*t6**3*t4**4*t3*t2+10._ki*t6**3*t4**2*t3*t2**3+10._ki*t6**3*t4**4&
                  &*t5**2*t2-24._ki*t6**3*t4**3*t3*t2**2-4._ki/3._ki*t6**3*t4*t3*t2**&
                  &4-5._ki*t6**3*t4**3*t5**2*t2**2+9._ki*t6**3*t1*t4**3*t2**2
                !
                stemp4=stemp5+16._ki/3._ki*t6**3*t4**5*t5*t2+5._ki/6._ki*t6**3*t4**2*&
                  &t5**2*t2**3+4._ki*t6**3*t4**3*t5*t2**3-8._ki*t6**3*t4**4*t5*t2**2&
                  &-2._ki/3._ki*t6**3*t4**2*t5*t2**4-4._ki*t6**3*t1*t4*t5*t2**3-25._ki&
                  &*t6**3*t4**2*t3*t5*t2**2-2._ki*t6**2*t4**4*t2**4-35._ki/6._ki*t6**&
                  &2*t4**2*t3*t2**4+7._ki/12._ki*t6**2*t4**2*t5*t2**5-7._ki/3._ki*t6**&
                  &2*t4**3*t5*t2**4+7._ki/3._ki*t6**2*t4**4*t5*t2**3+7._ki*t6**2*t4**&
                  &3*t3*t2**3-140._ki/3._ki*t6**3*t4**4*t3*t5-10._ki*t6**3*t1**2*t5**&
                  &2*t4+5._ki/2._ki*t6**3*t1**2*t5**2*t2+2._ki*t6**3*t1**2*t2*t3-2._ki&
                  &*t6**3*t1**2*t5*t2**2+60._ki*t6**3*t4**3*t3*t5*t2-40._ki*t6**3*t1&
                  &*t3*t5*t2*t4+5._ki*t6**3*t1*t4*t5**2*t2**2-24._ki*t6**3*t1*t4**3*&
                  &t5*t2
                !
                stemp5=-26._ki*t6**3*t1*t2*t3*t4**2-5._ki*t6**2*t1*t5**2*t3*t2**2-6&
                  &5._ki/2._ki*t6**2*t2*t3**2*t4**2*t5+20._ki*t6**2*t3**2*t4*t5*t2**2&
                  &+8._ki*t6**3*t1**2*t5*t2*t4-35._ki/6._ki*t6**2*t1*t4**2*t5*t2**3+5&
                  &._ki/3._ki*t6**2*t2**2*t3**3-t1**2*t2**4*t6**2/2._ki-20._ki/3._ki*t6&
                  &**3*t4**5*t5**2+t6**3*t4**2*t2**5/4._ki-3._ki/2._ki*t6**3*t4**3*t2&
                  &**4+3._ki/4._ki*t6**3*t1**2*t2**3-2._ki*t6**3*t4**5*t2**2+t4**2*t5&
                  &*t2**7/12._ki+t4*t3*t2**7/6._ki+2._ki*t6**2*t4**3*t2**5+5._ki/6._ki*&
                  &t6**3*t2**3*t3**2-95._ki/3._ki*t6**3*t3**2*t4**3+7._ki/3._ki*t6**2*&
                  &t1*t4*t5*t2**4+t6*t4**2*t2**7/4._ki+3._ki*t6**3*t4**4*t2**3-t6**2&
                  &*t4**2*t2**6/2._ki
                !
                stemp3=stemp5-t6*t4**3*t2**6/2._ki+25._ki*t6**2*t4**2*t5**2*t3*t2**&
                  &2+stemp4-30._ki*t6**2*t4**3*t5**2*t3*t2-5._ki*t6**2*t4*t5**2*t3*t&
                  &2**3-10._ki/3._ki*t6**2*t1*t4*t5**3*t2**2+25._ki/3._ki*t6**2*t1*t4*&
                  &*2*t5**3*t2+10._ki/3._ki*t6**3*t4*t3*t5*t2**3+7._ki/12._ki*t6**2*t1&
                  &**2*t5*t2**3-5._ki/6._ki*t6**2*t1**2*t2*t5**3+7._ki/6._ki*t6**2*t1*&
                  &t2**4*t3-5._ki*t6**2*t2*t3**3*t4-5._ki/2._ki*t6**2*t2**3*t3**2*t5-&
                  &5._ki*t6**3*t1*t3**2*t2+15._ki*t6**3*t3**2*t1*t4+65._ki/2._ki*t6**3&
                  &*t3**2*t2*t4**2-10._ki*t6**3*t2**2*t3**2*t4-15._ki/2._ki*t6**3*t1*&
                  &t4**2*t2**3+3._ki/2._ki*t6**3*t1*t4*t2**4+30._ki*t6**3*t1*t4**3*t5&
                  &**2-2._ki*t6**3*t1*t2**3*t3-5._ki*t6**3*t1**2*t3*t5-3._ki*t6**3*t1&
                  &**2*t2**2*t4
                !
                stemp4=1._ki/t2**10*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=-(-2._ki*t3*t1*t6+2._ki*t3**2*t2+3._ki*t6*t4*t2*t3-6._ki*t4**2&
                  &*t3*t6+t3*t5*t2*t4-4._ki*t4*t1*t6*t5)*t4/t2**4/t3*q(4,(t2*t3-t1*&
                  &t6)/t2/t3,sign_arg)/12._ki
                !
                stemp9=22._ki/3._ki*t6**3*t1**3*t4**3*t5*t2**3+25._ki/6._ki*t6**3*t1*&
                  &*4*t4*t5**2*t2**2+7._ki/3._ki*t6**2*t1**2*t4**4*t3*t2**4-7._ki/6._k&
                  &i*t6**2*t1**2*t4**3*t3*t2**5+2._ki/3._ki*t6**3*t1**5*t2*t3+15._ki/&
                  &2._ki*t6**3*t1**4*t4**3*t2**2-25._ki/4._ki*t6**3*t1**4*t4**2*t2**3&
                  &+t1**3*t4*t3*t2**7/18._ki+5._ki/4._ki*t6**3*t1**4*t4*t2**4+25._ki*t&
                  &6**3*t1**4*t4**3*t5**2-5._ki/3._ki*t6**3*t1**4*t2**3*t3-11._ki/4._k&
                  &i*t6**3*t1**3*t4**3*t2**4-11._ki/3._ki*t6**3*t1**3*t4**5*t2**2+55&
                  &._ki/36._ki*t6**3*t1**3*t3**2*t2**3-110._ki/9._ki*t6**3*t1**3*t4**5&
                  &*t5**2-1045._ki/18._ki*t6**3*t1**3*t3**2*t4**3+11._ki/24._ki*t6**3*&
                  &t1**3*t4**2*t2**5
                !
                stemp8=stemp9+11._ki/2._ki*t6**3*t1**3*t4**4*t2**3+5._ki/2._ki*t6**2*&
                  &t1**3*t3**3*t2**2+5._ki/2._ki*t6**2*t1**4*t4**2*t2**4+7._ki/12._ki*&
                  &t6**2*t1**4*t2**4*t3+5._ki/12._ki*t6**3*t1*t4**3*t3**2*t2**4-70._k&
                  &i/9._ki*t6**2*t1*t4**4*t3**3*t2**2-25._ki/18._ki*t6**2*t1*t2**4*t3&
                  &**3*t4**2+20._ki/3._ki*t6**2*t1**2*t4**5*t3*t5**2*t2+70._ki/3._ki*t&
                  &6**2*t1**2*t3**2*t4**4*t5*t2-14._ki/9._ki*t6**2*t1**2*t4**5*t3*t2&
                  &**3-30._ki*t6**2*t1**2*t3**2*t4**3*t5*t2**2-5._ki/18._ki*t6**2*t1*&
                  &*2*t3**3*t2**4-10._ki/9._ki*t6**2*t4**4*t3**3*t2**4+20._ki/9._ki*t6&
                  &**2*t4**5*t3**3*t2**3-20._ki/9._ki*t6**2*t4**6*t3**3*t2**2+25._ki/&
                  &2._ki*t6**3*t1**4*t3**2*t4-25._ki/6._ki*t6**3*t1**4*t3**2*t2
                !
                stemp9=8._ki/9._ki*t6**2*t4**7*t3**3*t2-t6**2*t1**4*t4*t2**5-3._ki/4&
                  &._ki*t6**2*t1**3*t4**2*t2**6-3._ki*t6**2*t1**3*t4**4*t2**4+3._ki*t&
                  &6**2*t1**3*t4**3*t2**5-t6*t1**3*t4**3*t2**6/2._ki+t6*t1**3*t4**2&
                  &*t2**7/4._ki-t6**2*t2**6*t3**3*t4**2/36._ki+5._ki/18._ki*t6**2*t2**&
                  &5*t3**3*t4**3+t1**3*t4**2*t5*t2**7/9._ki-70._ki/3._ki*t6**3*t1**2*&
                  &t4**4*t3**2*t2+55._ki/9._ki*t6**3*t1**3*t4*t3*t5*t2**3+15._ki*t6**&
                  &3*t1**2*t4**3*t3**2*t2**2+16._ki/3._ki*t6**3*t1**2*t4**5*t3*t2**2&
                  &-44._ki/3._ki*t6**3*t1**3*t4**4*t5*t2**2-55._ki/3._ki*t6**3*t1**3*t&
                  &4*t3**2*t2**2-22._ki/9._ki*t6**3*t1**3*t4*t3*t2**4
                !
                stemp7=stemp9+stemp8-40._ki/3._ki*t6**3*t1**2*t4**5*t3*t5*t2+10._ki*&
                  &t6**3*t1**2*t4**4*t3*t5*t2**2-25._ki/6._ki*t6**3*t1**2*t4**2*t3**&
                  &2*t2**3+4._ki/3._ki*t6**3*t1**2*t4**3*t3*t2**4-770._ki/9._ki*t6**3*&
                  &t1**3*t4**4*t3*t5+55._ki/3._ki*t6**3*t1**3*t4**4*t5**2*t2-55._ki/6&
                  &._ki*t6**3*t1**3*t4**3*t5**2*t2**2+25._ki/2._ki*t6**2*t1**2*t3**2*&
                  &t4**2*t5*t2**3-10._ki*t6**2*t1**2*t4**4*t3*t5**2*t2**2-11._ki/9._k&
                  &i*t6**3*t1**3*t4**2*t5*t2**4+7._ki/6._ki*t6**2*t1**4*t4*t5*t2**4+&
                  &25._ki/6._ki*t6**2*t1**4*t4**2*t5**3*t2+10._ki/3._ki*t6**3*t1*t4**5&
                  &*t3**2*t2**2+8._ki/3._ki*t6**3*t1**5*t5*t2*t4+20._ki/3._ki*t6**3*t1&
                  &**2*t4**6*t3*t5-8._ki/3._ki*t6**3*t1**2*t4**6*t3*t2-4._ki*t6**3*t1&
                  &**2*t4**4*t3*t2**3
                !
                stemp9=stemp7-35._ki/12._ki*t6**2*t1**4*t4**2*t5*t2**3+110._ki*t6**3&
                  &*t1**3*t4**3*t3*t5*t2+715._ki/12._ki*t6**3*t1**3*t3**2*t4**2*t2+5&
                  &5._ki/36._ki*t6**3*t1**3*t4**2*t5**2*t2**3+40._ki/9._ki*t6**2*t1*t4&
                  &**5*t3**3*t2+88._ki/9._ki*t6**3*t1**3*t4**5*t5*t2-10._ki/3._ki*t6**&
                  &3*t1**2*t4**3*t3*t5*t2**3-5._ki/3._ki*t6**2*t1**4*t4*t5**3*t2**2-&
                  &7._ki/3._ki*t6**2*t1**4*t2**3*t3*t4+55._ki/3._ki*t6**3*t1**3*t4**2*&
                  &t3*t2**3+40._ki/3._ki*t6**3*t1**2*t4**5*t3**2+4._ki/3._ki*t6**3*t1*&
                  &t4**7*t3**2-5._ki/3._ki*t6**3*t1**5*t3*t5-t6**3*t1**5*t2**2*t4-10&
                  &._ki/3._ki*t6**3*t1**5*t5**2*t4+5._ki/6._ki*t6**3*t1**5*t5**2*t2
                !
                stemp8=stemp9-2._ki/3._ki*t6**3*t1**5*t5*t2**2-44._ki*t6**3*t1**3*t4&
                  &**3*t3*t2**2-10._ki/3._ki*t6**3*t1*t4**6*t3**2*t2+5._ki*t6**2*t1*t&
                  &4**4*t3**2*t5*t2**3+5._ki*t6**2*t1*t2**3*t3**3*t4**3+75._ki/2._ki*&
                  &t6**2*t1**3*t4**2*t5**2*t3*t2**2+308._ki/9._ki*t6**3*t1**3*t4**4*&
                  &t3*t2-15._ki/4._ki*t6**2*t1**3*t2**3*t3**2*t5-45._ki*t6**2*t1**3*t&
                  &4**3*t5**2*t3*t2+10._ki*t6**2*t1**4*t5**2*t3*t2*t4+5._ki/2._ki*t6*&
                  &*2*t1**4*t3**2*t5*t2-195._ki/4._ki*t6**2*t1**3*t3**2*t5*t2*t4**2-&
                  &35._ki/4._ki*t6**2*t1**3*t4**2*t3*t2**4-275._ki/6._ki*t6**3*t1**3*t&
                  &4**2*t3*t5*t2**2+7._ki/4._ki*t6**2*t1**3*t4*t3*t2**5-5._ki/3._ki*t6&
                  &**3*t1*t4**4*t3**2*t2**3+25._ki/6._ki*t6**3*t1**4*t3*t5*t2**2-5._k&
                  &i/2._ki*t6**2*t1**4*t5**2*t3*t2**2
                !
                stemp9=stemp8-t6**3*t1**2*t4**2*t3*t2**5/6._ki+30._ki*t6**2*t1**3*t&
                  &3**2*t5*t2**2*t4+21._ki/2._ki*t6**2*t1**3*t4**3*t3*t2**3-125._ki/6&
                  &._ki*t6**3*t1**4*t4**2*t5**2*t2+7._ki/8._ki*t6**2*t1**3*t4**2*t5*t&
                  &2**5+95._ki/9._ki*t6**2*t1**2*t3**3*t4**3*t2-65._ki/6._ki*t6**2*t1*&
                  &*2*t3**3*t4**2*t2**2-20._ki/3._ki*t6**2*t1*t4**5*t3**2*t5*t2**2+1&
                  &0._ki/3._ki*t6**2*t1*t4**6*t3**2*t5*t2-15._ki/2._ki*t6**2*t1**3*t3*&
                  &*3*t2*t4-7._ki/2._ki*t6**2*t1**3*t4**3*t5*t2**4+7._ki/2._ki*t6**2*t&
                  &1**3*t4**4*t5*t2**3+40._ki/3._ki*t6**3*t1**4*t2**2*t3*t4+50._ki/3.&
                  &_ki*t6**3*t1**4*t4**2*t5*t2**2+325._ki/6._ki*t6**3*t1**4*t4**2*t3*&
                  &t5-100._ki/3._ki*t6**3*t1**4*t3*t5*t2*t4-t6**3*t1*t4**2*t3**2*t2*&
                  &*5/24._ki
                !
                stemp6=stemp9+5._ki/12._ki*t6**3*t1**2*t4*t3**2*t2**4+5._ki/12._ki*t6&
                  &**3*t1**2*t4**2*t5*t2**4*t3+5._ki/36._ki*t6**2*t1*t2**5*t3**3*t4+&
                  &5._ki*t6**2*t1**2*t4**3*t3*t5**2*t2**3+10._ki/3._ki*t6**2*t1**2*t4&
                  &*t3**3*t2**3+t6**3*t1**5*t2**3/4._ki-5._ki/3._ki*t6**2*t1*t4**3*t3&
                  &**2*t5*t2**4+5._ki*t6**2*t1**3*t4**3*t5**3*t2**2-5._ki*t6**2*t1**&
                  &3*t4**4*t5**3*t2-10._ki/3._ki*t6**3*t1**4*t4*t5*t2**3-20._ki*t6**3&
                  &*t1**4*t4**3*t5*t2+5._ki/24._ki*t6**2*t1*t4**2*t3**2*t2**5*t5-5._k&
                  &i/6._ki*t6**2*t1**2*t4**2*t3*t5**2*t2**4-5._ki/3._ki*t6**2*t1**2*t&
                  &4*t3**2*t5*t2**4+7._ki/36._ki*t6**2*t1**2*t4**2*t3*t2**6-65._ki/3.&
                  &_ki*t6**3*t1**4*t2*t3*t4**2-5._ki/4._ki*t6**2*t1**3*t4**2*t5**3*t2&
                  &**3-15._ki/2._ki*t6**2*t1**3*t4*t5**2*t3*t2**3
                !
                stemp7=1._ki/t1**3/t2**10
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
              !
            case(3)
              !
              select case(par4_glob)
              !
              case(3)
                !
                stemp7=-660._ki*t6*t4*t3**2*t5**2*t2**4+450._ki*t6*t4**2*t5**3*t3**&
                  &2*t2**2-12._ki*t6*t1*t4**2*t5**3*t2**4+3._ki*t6*t1*t4**2*t5**2*t2&
                  &**5-15._ki*t6*t1*t4**2*t5**5*t2**2-30._ki*t6*t1*t4**2*t5**4*t2**3&
                  &-3._ki*t6*t1*t4**2*t2**6*t5-140._ki*t6**2*t1*t4*t5**3*t2**4-400._k&
                  &i*t6**2*t4*t5**3*t2**4*t3+30._ki*t6**2*t1**2*t4*t5**4*t2-15._ki/2&
                  &._ki*t6**2*t1**2*t4*t5**2*t2**3+30._ki*t6**2*t1**2*t4*t5**3*t2**2&
                  &+35._ki*t6**2*t1*t4**3*t5**2*t2**3+260._ki*t6**2*t1*t4**2*t5**4*t&
                  &2**2+7._ki*t6**2*t1*t4**3*t5*t2**4-140._ki*t6**2*t1*t4**3*t5**4*t&
                  &2-13._ki*t6**2*t1*t4**2*t5*t2**5
                !
                stemp6=stemp7-65._ki*t6**2*t1*t4**2*t5**2*t2**4+260._ki*t6**2*t1*t4&
                  &**2*t5**3*t2**3-140._ki*t6**2*t1*t4**3*t5**3*t2**2-3._ki/2._ki*t6*&
                  &*2*t1**2*t4*t5*t2**4+90._ki*t6**2*t1*t4**2*t5*t2**3*t3+660._ki*t6&
                  &**2*t1*t4*t3*t5**2*t2**3-110._ki*t6**2*t1*t4*t5*t2**4*t3-540._ki*&
                  &t6**2*t1*t4**2*t3*t5**2*t2**2-720._ki*t6**2*t1*t4**2*t5**3*t2*t3&
                  &+880._ki*t6**2*t1*t4*t5**3*t3*t2**2-120._ki*t6**2*t4**4*t5**4*t2*&
                  &*2+40._ki*t6**2*t4**5*t5**4*t2-125._ki/2._ki*t6**2*t4**3*t3**2*t2*&
                  &*3-8._ki*t6**2*t4**4*t3*t2**4+130._ki*t6**2*t4**3*t5**3*t2**4-120&
                  &._ki*t6**2*t4**4*t5**3*t2**3-13._ki/2._ki*t6**2*t4**3*t5*t2**6+6._k&
                  &i*t6**2*t4**4*t5*t2**5
                !
                stemp5=-2._ki*t6**2*t4**5*t5*t2**4-10._ki*t6**2*t4**5*t5**2*t2**3+3&
                  &0._ki*t6**2*t4**4*t5**2*t2**4-65._ki/2._ki*t6**2*t4**3*t5**2*t2**5&
                  &-50._ki*t6*t4**3*t5**4*t2**4+20._ki*t6*t4**4*t5**4*t2**3-5._ki*t6*&
                  &t4**3*t5*t2**7+2._ki*t6*t4**4*t5*t2**6-25._ki*t6*t4**3*t5**5*t2**&
                  &3+10._ki*t6*t4**4*t5**5*t2**2+8._ki*t6*t4**4*t5**3*t2**4-13._ki*t6&
                  &*t4**2*t3*t2**7-20._ki*t6*t4**3*t5**3*t2**5+7._ki*t6*t4**3*t3*t2*&
                  &*6+11._ki*t6*t4*t3**2*t2**6-9._ki*t6*t4**2*t3**2*t2**5+5._ki*t6*t4&
                  &**3*t5**2*t2**6-2._ki*t6*t4**4*t5**2*t2**5-13._ki*t6**3*t4*t3*t2*&
                  &*6-240._ki*t6**2*t3**2*t2**4*t5**2-60._ki*t6**2*t4**2*t5**4*t2**4&
                  &+t6**2*t1**2*t5*t2**5+5._ki*t6**2*t1**2*t5**2*t2**4-20._ki*t6**2*&
                  &t1**2*t5**3*t2**3-20._ki*t6**2*t1**2*t5**4*t2**2+40._ki*t6**3*t4*&
                  &*2*t5**3*t2**4-8._ki*t6**3*t4**2*t5*t2**6+t6*t1*t2**5*t3**2-220.&
                  &_ki*t3**3*t2**3*t6**2*t4-120._ki*t3**3*t2**4*t6*t5+190._ki*t3**3*t&
                  &2**2*t6**2*t4**2+9._ki*t3*t2**3*t6**3*t1**2-15._ki*t3**2*t2**4*t1&
                  &*t6**2+20._ki*t3**3*t2**4*t6*t4+25._ki*t6*t5*t2**2*t3**4+stemp6
                !
                stemp6=-535._ki*t2*t3**3*t6**3*t4**2+290._ki*t2**2*t3**3*t6**3*t4-1&
                  &05._ki*t3**3*t6**3*t1*t4+3._ki*t6**2*t2**6*t1*t3+2._ki*t6*t1*t3*t2&
                  &**7+t2**2*t1**3*t6**3*t5+45._ki*t6**3*t5*t1**2*t3**2-150._ki*t2**&
                  &3*t3**3*t6*t5**2+240._ki*t2**3*t3**3*t6**2*t5+4._ki*t6*t4**2*t2**&
                  &8*t5+180._ki*t3**2*t2**3*t6**2*t5*t1-300._ki*t3**2*t2**2*t6**2*t4&
                  &*t5*t1+200._ki*t3**3*t2**3*t6*t4*t5-15._ki*t3*t2**2*t6**3*t4*t1**&
                  &2+25._ki*t3**2*t2**3*t1*t6**2*t4-5._ki*t3*t2**3*t1**2*t6**2*t5+6.&
                  &_ki*t6*t1*t2**5*t3*t4*t5-3._ki*t6*t1*t2**6*t4*t3-12._ki*t3**2*t2**&
                  &4*t1*t6*t5-80._ki*t6**2*t5*t2*t1*t3**3-135._ki*t2*t3*t1**2*t6**3*&
                  &t5**2+360._ki*t2**2*t3**2*t1*t6**2*t5**2-855._ki*t3**2*t1*t6**3*t&
                  &5*t4**2+t4**2*t2**10/2._ki-t4**3*t2**9/4._ki+t6**2*t4*t2**9/4._ki-&
                  &3._ki/2._ki*t6**2*t4**2*t2**8+60._ki*t3**3*t2**4*t6**2-5._ki*t1**3*&
                  &t6**3*t5**3+10._ki*t6*t2**3*t3**4-40._ki*t3**4*t2**2*t6**2+40._ki*&
                  &t6**3*t4**6*t5**3-3._ki*t6*t3**2*t2**7+stemp5-630._ki*t6**3*t3**2&
                  &*t2**3*t4*t5
                !
                stemp7=stemp6-1080._ki*t6**3*t1*t4*t3*t5**2*t2**2+1980._ki*t6**3*t1&
                  &*t4**2*t3*t5**2*t2+33._ki*t6**3*t1**2*t4*t5*t2**3-27._ki*t6**3*t1&
                  &**2*t4**2*t5*t2**2-165._ki*t6**3*t1**2*t4*t5**3*t2+75._ki*t6**3*t&
                  &1*t2**2*t3*t4**3+48._ki*t6**3*t1*t4**4*t5*t2**2-1125._ki*t6**3*t1&
                  &*t4**3*t3*t5**2+72._ki*t6**3*t1*t2**4*t3*t4-465._ki*t6**3*t1*t4**&
                  &2*t5**3*t2**2-132._ki*t6**3*t1*t2**3*t3*t4**2+570._ki*t6**3*t1*t4&
                  &**3*t5**3*t2-114._ki*t6**3*t1*t4**3*t5*t2**3+93._ki*t6**3*t1*t4**&
                  &2*t5*t2**4-855._ki*t6**3*t4**2*t3*t5**2*t2**3+1695._ki*t6**3*t4**&
                  &3*t3*t5**2*t2**2-1560._ki*t6**3*t4**4*t3*t5**2*t2
                !
                stemp4=stemp7-2265._ki*t6**3*t4**3*t3**2*t5*t2+1845._ki*t6**3*t4**2&
                  &*t3**2*t5*t2**2+40._ki*t6**2*t4**5*t5**3*t2**2+40._ki*t2*t3*t1**2&
                  &*t6**2*t5**3+225._ki*t3*t1**2*t6**3*t4*t5**2-270._ki*t2**2*t3**2*&
                  &t1*t6**3*t5+990._ki*t2*t3**2*t1*t6**3*t5*t4-600._ki*t2*t3**2*t1*t&
                  &6**2*t4*t5**2-880._ki*t2**2*t3**3*t6**2*t5*t4+760._ki*t2*t3**3*t6&
                  &**2*t5*t4**2+175._ki*t6*t4*t5**4*t2**4*t3-14._ki*t6*t4*t3*t2**7*t&
                  &5+280._ki*t6*t4*t5**3*t2**5*t3+250._ki*t2**2*t3**3*t6*t4*t5**2+15&
                  &0._ki*t6**3*t1*t4*t5**3*t2**3+180._ki*t6**3*t1*t3*t5**2*t2**3+84.&
                  &_ki*t6*t4*t5**2*t2**6*t3+50._ki*t6*t1*t5**4*t3*t2**3+80._ki*t6*t1*&
                  &t5**3*t3*t2**4
                !
                stemp6=stemp4-300._ki*t6**2*t4*t3*t5**2*t2**5+35._ki*t6**2*t1*t4*t5&
                  &**2*t2**5-140._ki*t6**2*t1*t4*t5**4*t2**3+1440._ki*t6**2*t4*t3**2&
                  &*t5**2*t2**3+1500._ki*t6**2*t4**3*t3**2*t5**2*t2+720._ki*t6**2*t4&
                  &*t3**2*t5*t2**4-1320._ki*t6**2*t4**2*t3**2*t5*t2**3-1520._ki*t6**&
                  &2*t4**3*t3*t5**3*t2**2-2640._ki*t6**2*t4**2*t3**2*t5**2*t2**2+75&
                  &0._ki*t6**2*t4**3*t3**2*t5*t2**2+640._ki*t6**2*t4**4*t3*t5**3*t2-&
                  &155._ki*t6**2*t4**2*t5*t2**5*t3+190._ki*t6**2*t4**3*t5*t2**4*t3+4&
                  &._ki*t6*t1*t4*t5*t2**7+24._ki*t6*t1*t3*t5**2*t2**5-4._ki*t6*t1*t4*&
                  &t5**2*t2**6+195._ki*t6**3*t4*t3*t5**2*t2**4+180._ki*t6*t3**2*t2**&
                  &5*t5**2-4._ki*t6*t4**2*t5**2*t2**7+150._ki*t6*t3**2*t5**3*t2**4+7&
                  &._ki*t6*t4*t3*t2**8+40._ki*t6*t4**2*t5**4*t2**5-60._ki*t6**2*t4**2&
                  &*t5**3*t2**5+3._ki*t6**2*t4**2*t5*t2**7+5._ki*t6**2*t4*t3*t2**7-1&
                  &20._ki*t6**2*t2**5*t3**2*t5-7._ki/2._ki*t6**2*t1*t4*t2**7+75._ki*t6&
                  &**3*t3**2*t2**4*t5+45._ki*t6**3*t1**2*t5**3*t2**2-9._ki*t6**3*t1*&
                  &*2*t5*t2**4-12._ki*t6**3*t1*t3*t2**5-t6**2*t1**2*t3*t2**4/2._ki+6&
                  &0._ki*t3**3*t2*t6**3*t1+70._ki*t3**4*t2*t6**2*t4
                !
                stemp5=stemp6-20._ki*t6**2*t2**2*t1*t3**3+36._ki*t6*t5*t3**2*t2**6+&
                  &16._ki*t6*t4**2*t5**3*t2**6+15._ki*t6**2*t4**2*t5**2*t2**6+20._ki*&
                  &t6*t4**2*t5**5*t2**4+t6**3*t4*t5*t2**7-5._ki*t6**3*t4*t5**3*t2**&
                  &5-15._ki*t6**3*t3*t5**2*t2**5-15._ki*t6**3*t1*t5**3*t2**4+3._ki*t6&
                  &**3*t1*t5*t2**6+10._ki*t6**2*t4*t5**3*t2**6-5._ki/2._ki*t6**2*t4*t&
                  &5**2*t2**7-5._ki*t6**2*t5*t2**7*t3+10._ki*t6**2*t4*t5**4*t2**5+30&
                  &._ki*t6**2*t3*t5**2*t2**6+40._ki*t6**2*t5**3*t2**5*t3-t6*t1*t5*t2&
                  &**8-40._ki*t6*t3*t5**3*t2**6-25._ki*t6*t5**4*t2**5*t3+2._ki*t6*t5*&
                  &t3*t2**8-t6**2*t4*t2**8*t5/2._ki-5._ki*t6**2*t1*t5**2*t2**6-t6**2&
                  &*t1*t5*t2**7+20._ki*t6**2*t1*t5**4*t2**4+20._ki*t6**2*t1*t5**3*t2&
                  &**5-10._ki*t6*t4*t5**4*t2**6-t6*t4*t2**9*t5-4._ki*t6*t4*t5**3*t2*&
                  &*7+t6*t2**8*t4*t5**2-5._ki*t6*t4*t5**5*t2**5-12._ki*t6*t3*t5**2*t&
                  &2**7-10._ki*t6*t1*t5**4*t2**5-5._ki*t6*t1*t5**5*t2**4-4._ki*t6*t1*&
                  &t5**3*t2**6+t6*t1*t5**2*t2**7
                !
                stemp6=stemp5+990._ki*t6**3*t4**4*t3**2*t5+25._ki*t6**3*t4**3*t5*t2&
                  &**5-8._ki*t6**3*t4**6*t5*t2**2-125._ki*t6**3*t4**3*t5**3*t2**3+57&
                  &._ki*t6**3*t4**2*t3*t2**5-36._ki*t6**3*t4**5*t3*t2**2-140._ki*t6**&
                  &3*t4**5*t5**3*t2+190._ki*t6**3*t4**4*t5**3*t2**2-38._ki*t6**3*t4*&
                  &*4*t5*t2**4+104._ki*t6**3*t4**4*t3*t2**3-113._ki*t6**3*t4**3*t3*t&
                  &2**4+540._ki*t6**3*t4**5*t3*t5**2+28._ki*t6**3*t4**5*t5*t2**3+135&
                  &._ki*t6**3*t1**2*t4**2*t5**3-240._ki*t6**3*t1*t4**4*t5**3+3._ki/4.&
                  &_ki*t6**2*t1**2*t4*t2**5+13._ki/2._ki*t6**2*t1*t4**2*t2**6-7._ki/2.&
                  &_ki*t6**2*t1*t4**3*t2**5-60._ki*t6**2*t4*t3**2*t2**5+110._ki*t6**2&
                  &*t4**2*t3**2*t2**4-31._ki/2._ki*t6**2*t4**2*t3*t2**6+19._ki*t6**2*&
                  &t4**3*t3*t2**5+130._ki*t6**2*t4**3*t5**4*t2**3+50._ki*t6**2*t4*t5&
                  &*t2**6*t3+1240._ki*t6**2*t4**2*t5**3*t3*t2**3-1140._ki*t6**2*t4**&
                  &3*t3*t5**2*t2**3+480._ki*t6**2*t4**4*t3*t5**2*t2**2-80._ki*t6**2*&
                  &t4**4*t5*t2**3*t3+930._ki*t6**2*t4**2*t3*t5**2*t2**4-75._ki*t6*t1&
                  &*t4*t5**4*t2**2*t3+7._ki*t6**2*t1*t2**6*t4*t5-180._ki*t6**2*t1*t3&
                  &*t5**2*t2**4+9._ki*t6**2*t1*t3*t2**4*t4**2-11._ki*t6**2*t1*t3*t2*&
                  &*5*t4-50._ki*t2**3*t3**3*t6**3
                !
                stemp7=stemp6+315._ki*t3**3*t6**3*t4**3+10._ki*t6**2*t3**2*t2**6-t6&
                  &**2*t1**2*t2**6/2._ki-12._ki*t3**3*t6*t2**5+t6**3*t3*t2**7-t6**2*&
                  &t3*t2**8/2._ki+t6**2*t1*t2**8/2._ki-t6*t3*t2**9+13._ki/4._ki*t6**2*&
                  &t4**3*t2**7-3._ki*t6**2*t4**4*t2**6+t6**2*t4**5*t2**5+30._ki*t3*t&
                  &5**2*t2**2*t1**2*t6**2-60._ki*t3**2*t5**2*t2**3*t1*t6-t4*t2**11/&
                  &4._ki-50._ki*t2**2*t3**2*t1*t6*t5**3-30._ki*t6**3*t1*t4*t5*t2**5-2&
                  &40._ki*t6**2*t1*t3*t5**3*t2**3
                !
                stemp3=stemp7+30._ki*t6**2*t1*t3*t5*t2**5-4._ki*t6*t1*t2**6*t3*t5+1&
                  &6._ki*t6*t1*t4*t5**3*t2**5+40._ki*t6*t1*t4*t5**4*t2**4+20._ki*t6*t&
                  &1*t4*t5**5*t2**3-36._ki*t6*t1*t4*t3*t2**4*t5**2-550._ki*t6*t4*t5*&
                  &*3*t3**2*t2**3-520._ki*t6*t4**2*t3*t5**3*t2**4-325._ki*t6*t4**2*t&
                  &5**4*t3*t2**3-132._ki*t6*t4*t5*t2**5*t3**2+108._ki*t6*t4**2*t3**2&
                  &*t5*t2**4-156._ki*t6*t4**2*t5**2*t2**5*t3+84._ki*t6*t4**3*t5**2*t&
                  &2**4*t3+540._ki*t6*t4**2*t3**2*t5**2*t2**3+175._ki*t6*t4**3*t3*t5&
                  &**4*t2**2+280._ki*t6*t4**3*t3*t5**3*t2**3-120._ki*t6*t1*t4*t5**3*&
                  &t3*t2**3-14._ki*t6*t4**3*t5*t3*t2**5+26._ki*t6*t4**2*t5*t3*t2**6
                !
                stemp4=1._ki/t2**12*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=(-t4+t2)**2*t4/t2**3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/4._k&
                  &i
                !
                stemp11=20._ki*t6*t3**2*t5**2*t2**6*t1**4-5._ki/2._ki*t6*t4*t5*t3**2&
                  &*t2**8*t1**3-5._ki*t6*t3**2*t4*t2**7*t1**3*t5**2+4._ki/3._ki*t6*t3&
                  &**3*t4*t2**8*t1**2*t5+15._ki*t6*t1**5*t4*t5**4*t2**5+45._ki*t6*t1&
                  &**5*t3*t5**2*t2**6+60._ki*t6*t1**5*t5**3*t2**5*t3+10._ki*t6*t5*t3&
                  &**2*t2**7*t1**4-10._ki/3._ki*t6*t2**7*t3**3*t1**3*t5+4._ki/3._ki*t6&
                  &*t1**4*t4**2*t3*t2**8-19._ki/3._ki*t6*t4*t3**4*t2**7*t1**2-600._ki&
                  &*t6*t1**5*t4*t5**3*t2**4*t3-25._ki/12._ki*t6*t4**2*t3**2*t2**8*t1&
                  &**3-450._ki*t6*t1**5*t4*t3*t5**2*t2**5+9._ki/2._ki*t6*t1**5*t4**2*&
                  &t5*t2**7+15._ki/2._ki*t6*t1**5*t4*t3*t2**7-180._ki*t6*t1**5*t2**5*&
                  &t3**2*t5+45._ki/2._ki*t6*t1**5*t4**2*t5**2*t2**6+240._ki*t6**2*t1*&
                  &*4*t4**6*t3*t5**2*t2
                !
                stemp10=stemp11-255._ki/2._ki*t6**2*t1**3*t4**4*t5*t2**4*t3**2-204.&
                  &_ki*t6**2*t1**3*t4**6*t3**2*t5*t2**2+183._ki/4._ki*t6**2*t1**3*t4*&
                  &*3*t5*t2**5*t3**2-24._ki*t6**2*t1**3*t4**8*t3**2*t5-278._ki*t6**2&
                  &*t1**3*t4**5*t3**3*t2**2-172._ki/3._ki*t6**2*t1**2*t4**7*t3**3*t2&
                  &**2-85._ki/12._ki*t6**2*t1**2*t4**3*t3**3*t2**6+73._ki/3._ki*t6**2*&
                  &t1**2*t4**4*t3**3*t2**5+208._ki/3._ki*t6**2*t1**2*t4**6*t3**3*t2*&
                  &*3+80._ki/3._ki*t6**2*t1**2*t4**8*t3**3*t2-155._ki/3._ki*t6**2*t1**&
                  &2*t4**5*t3**3*t2**4+108._ki*t6**2*t1**3*t4**7*t3**2*t5*t2-5._ki*t&
                  &6**2*t1**7*t2**2*t3*t4-9._ki*t6**2*t1**7*t4**2*t5*t2**2-1130._ki*&
                  &t6*t1**4*t4**3*t3**2*t5*t2**4+1040._ki*t6*t1**4*t4**4*t3**2*t5*t&
                  &2**3-360._ki*t6*t1**4*t4**5*t3**2*t5*t2**2-205._ki*t6*t1**3*t4**3&
                  &*t5**2*t2**5*t3**2-8305._ki/2._ki*t6**2*t1**5*t4**3*t3**2*t5*t2
                !
                stemp11=-910._ki/3._ki*t6*t1**3*t4**2*t3**3*t5*t2**5+440._ki*t6*t1**&
                  &3*t4**4*t3**2*t5**2*t2**4-260._ki*t6*t1**3*t4**5*t3**2*t5*t2**4-&
                  &680._ki/3._ki*t6*t1**2*t4**4*t5*t2**5*t3**3-310._ki/3._ki*t6*t1*t4*&
                  &*5*t3**4*t2**5+146._ki/3._ki*t6*t1*t4**4*t3**4*t2**6+25._ki/3._ki*t&
                  &4*t5**3*t2**7*t1**3*t3**2+125._ki/3._ki*t4**2*t5**2*t2**7*t1**2*t&
                  &3**3+108._ki*t1**5*t4**2*t3**2*t5*t2**4+84._ki*t1**5*t4**3*t5**2*&
                  &t2**4*t3-12._ki*t1**4*t4**4*t3*t2**6*t5+275._ki*t1**4*t3**4*t5*t2&
                  &**3*t4-600._ki*t1**4*t4*t3**3*t2**4*t5**2-195._ki/4._ki*t6*t1**5*t&
                  &4**3*t5**2*t2**5-150._ki*t6*t1**6*t4*t5*t3**2*t2**2-270._ki*t6*t1&
                  &**6*t4**2*t3*t5**2*t2**2-90._ki*t6*t1**5*t4**2*t5**3*t2**5-320._k&
                  &i/3._ki*t6*t1**4*t4**6*t5**3*t3*t2+1000._ki/3._ki*t6*t1**4*t4**3*t&
                  &5**3*t3*t2**4
                !
                stemp9=stemp11-775._ki*t1**4*t4**2*t3**2*t5**3*t2**4-340._ki*t1*t4*&
                  &*6*t3**4*t5*t2**4-1980._ki*t6*t1**5*t4**2*t3**2*t5*t2**3+175._ki*&
                  &t1**5*t4**3*t3*t5**4*t2**2+880._ki*t1**4*t4**2*t3**3*t5*t2**4-16&
                  &._ki*t4**2*t3**2*t2**8*t1**3*t5-96._ki*t1**4*t4**4*t5*t2**4*t3**2&
                  &+228._ki*t1**4*t4**3*t5*t2**5*t3**2+100._ki/3._ki*t4*t5*t2**7*t3**&
                  &4*t1**2+13._ki*t1**4*t4**3*t3*t2**7*t5+4._ki*t1**4*t4**5*t3*t2**5&
                  &*t5-500._ki*t1**4*t4**3*t3**3*t5*t2**3+3775._ki/6._ki*t1**3*t4**3*&
                  &t3**4*t2**3*t5+175._ki*t1**3*t4*t3**4*t2**5*t5-76._ki*t1**3*t3**2&
                  &*t4**4*t5*t2**6+56._ki*t1**3*t3**2*t4**5*t5*t2**5-2825._ki/3._ki*t&
                  &1**3*t4**3*t3**3*t5**2*t2**4+540._ki*t1**5*t4**2*t3**2*t5**2*t2*&
                  &*3-132._ki*t1**5*t4*t5*t2**5*t3**2+stemp10
                !
                stemp11=-325._ki*t1**5*t4**2*t5**4*t3*t2**3+280._ki*t1**5*t4**3*t3*&
                  &t5**3*t2**3-200._ki/3._ki*t1**2*t4**7*t3**3*t5**2*t2**2+9._ki/2._ki&
                  &*t6*t1**6*t4**2*t3*t2**4+825._ki*t6**2*t1**6*t3**2*t5*t2*t4-1040&
                  &._ki/3._ki*t1**2*t4**5*t3**3*t5*t2**5+1120._ki/3._ki*t6*t1**2*t4**5&
                  &*t3**3*t5*t2**4-55._ki*t6*t1**6*t4*t5*t2**4*t3+2250._ki*t6*t1**5*&
                  &t4**3*t3**2*t5**2*t2-1710._ki*t6*t1**5*t4**3*t3*t5**2*t2**3+720.&
                  &_ki*t6*t1**5*t4**4*t3*t5**2*t2**2+250._ki*t1**5*t4*t3**3*t5**2*t2&
                  &**2+450._ki*t1**5*t4**2*t5**3*t3**2*t2**2-550._ki*t1**5*t4*t5**3*&
                  &t3**2*t2**3-1025._ki/2._ki*t1**3*t4**2*t3**4*t2**4*t5+2080._ki/3._k&
                  &i*t1**3*t4**4*t3**3*t2**4*t5-240._ki*t1**3*t4**5*t3**3*t2**3*t5+&
                  &250._ki*t1**3*t4**3*t3**2*t5**2*t2**6-2260._ki/3._ki*t1**3*t4**3*t&
                  &3**3*t2**5*t5
                !
                stemp10=stemp11+380._ki*t1**3*t4**2*t3**3*t2**6*t5-380._ki*t1**3*t4&
                  &**4*t3**2*t5**2*t2**5+280._ki*t1**3*t4**5*t3**2*t5**2*t2**4+50._k&
                  &i*t1**3*t3**2*t4**3*t5*t2**7-16._ki*t1**3*t3**2*t4**6*t5*t2**4-9&
                  &3._ki/4._ki*t6*t1**5*t4**2*t3*t2**6+360._ki*t6*t1**5*t3**3*t5*t2**&
                  &3-375._ki/4._ki*t6*t1**5*t4**3*t3**2*t2**3-40._ki*t6*t1**3*t4**7*t&
                  &3**2*t5*t2**2+220._ki*t6*t1**3*t4**4*t3**2*t5*t2**5+160._ki*t6*t1&
                  &**3*t4**6*t3**2*t5*t2**3+210._ki*t6*t1**3*t4**3*t3**3*t2**5+19._k&
                  &i/3._ki*t6*t1**4*t4**4*t3*t2**6-210._ki*t6*t1**4*t4**3*t3**4*t2-4&
                  &0._ki/3._ki*t6*t1**3*t4**6*t3**2*t2**4-200._ki/3._ki*t6*t1**3*t4**6&
                  &*t3**3*t2**2-425._ki*t6*t1**3*t4**3*t3**4*t2**3-140._ki*t6*t1**3*&
                  &t4**5*t3**4*t2-325._ki/6._ki*t6*t1**3*t3**4*t4*t2**5
                !
                stemp11=10._ki/3._ki*t6*t1**3*t4**7*t3**2*t2**3+680._ki/3._ki*t6*t1**&
                  &3*t4**5*t3**3*t2**3-2275._ki/12._ki*t1**2*t4**2*t3**4*t5*t2**6+33&
                  &0._ki*t6*t1**6*t4*t3*t5**2*t2**3-1640._ki*t6*t1**4*t4**2*t3**3*t2&
                  &**3*t5-770._ki/3._ki*t6**2*t1**5*t4**5*t5**3*t2+75._ki*t6**2*t1**7&
                  &*t4*t3*t5**2-45._ki*t6**2*t1**7*t3*t5**2*t2+11._ki*t6**2*t1**7*t4&
                  &*t5*t2**3-55._ki*t6**2*t1**7*t4*t5**3*t2+1650._ki*t6**2*t1**6*t4*&
                  &*2*t3*t5**2*t2+1595._ki/3._ki*t6**2*t1**5*t3**3*t2**2*t4+572._ki/3&
                  &._ki*t6**2*t1**5*t4**4*t3*t2**3+275._ki/6._ki*t6**2*t1**5*t4**3*t5&
                  &*t2**5-130._ki*t6*t4*t3**2*t2**6*t1**4*t5-320._ki/3._ki*t6*t4**2*t&
                  &5**3*t2**5*t1**4*t3+25._ki*t6*t4**2*t5*t3**2*t2**7*t1**3-325._ki/&
                  &2._ki*t1**4*t4**3*t5**4*t3*t2**4+88._ki/3._ki*t1**2*t4**4*t3**3*t2&
                  &**7-455._ki/6._ki*t1**2*t4**2*t3**4*t2**7
                !
                stemp8=stemp11+210._ki*t1**2*t4**3*t3**4*t2**6-200._ki/3._ki*t1**2*t&
                  &4**6*t3**4*t2**3-41._ki/3._ki*t1**2*t4**3*t3**3*t2**8-16._ki/3._ki*&
                  &t1**2*t4**7*t3**3*t2**4+64._ki/3._ki*t1**2*t4**6*t3**3*t2**5-104.&
                  &_ki/3._ki*t1**2*t4**5*t3**3*t2**6-920._ki/3._ki*t1**2*t4**4*t3**4*t&
                  &2**5+680._ki/3._ki*t1**2*t4**5*t3**4*t2**4-70._ki*t1**2*t2**2*t3**&
                  &5*t4**5+110._ki*t1**2*t2**5*t3**5*t4**2+590._ki/3._ki*t1**2*t2**3*&
                  &t3**5*t4**4-325._ki/12._ki*t1**2*t4*t3**5*t2**6-425._ki/2._ki*t1**2&
                  &*t2**4*t3**5*t4**3-475._ki/4._ki*t1*t2**6*t3**5*t4**3-6._ki*t4**2*&
                  &t3**4*t2**9*t1-150._ki*t1**5*t3**3*t5**2*t2**3-13._ki*t1**5*t4**2&
                  &*t3*t2**7+11._ki*t1**5*t4*t3**2*t2**6+stemp9+stemp10
                !
                stemp11=-50._ki*t1**5*t4**3*t5**4*t2**4+5._ki*t1**5*t4**3*t5**2*t2*&
                  &*6-120._ki*t1**5*t3**3*t5*t2**4+20._ki*t1**5*t4**4*t5**4*t2**3+2.&
                  &_ki*t1**5*t4**4*t5*t2**6+25._ki*t1**5*t3**4*t5*t2**2-25._ki*t1**5*&
                  &t4**3*t5**5*t2**3+416._ki/3._ki*t6*t1*t2**4*t3**4*t4**6-85._ki/6._k&
                  &i*t6*t1*t4**3*t3**4*t2**7-44._ki/3._ki*t6**2*t4**2*t2**6*t1**5*t5&
                  &+220._ki/3._ki*t6**2*t4**2*t5**3*t2**4*t1**5+715._ki/2._ki*t6**2*t4&
                  &*t3*t2**4*t1**5*t5**2-143._ki/6._ki*t6**2*t4*t2**6*t3*t1**5+75._ki&
                  &/2._ki*t6**2*t4**2*t5**2*t2**5*t1**4*t3-9._ki*t6**2*t4**2*t5*t2**&
                  &6*t1**3*t3**2+320._ki*t6*t1**3*t4**6*t3**2*t5**2*t2**2-80._ki*t6*&
                  &t1**3*t4**7*t3**2*t5**2*t2+2720._ki/3._ki*t6*t1**3*t4**5*t3**3*t5&
                  &*t2**2-800._ki/3._ki*t6*t1**3*t4**6*t3**3*t5*t2
                !
                stemp10=stemp11+840._ki*t6*t1**3*t4**3*t3**3*t5*t2**4+880._ki/3._ki*&
                  &t1**2*t4**4*t3**3*t5*t2**6-19._ki/4._ki*t6**2*t4*t3**3*t2**6*t1**&
                  &3-t6*t3**4*t4*t2**9*t1/6._ki-t6**2*t3**3*t4*t2**8*t1**2/12._ki-6.&
                  &_ki*t1**4*t4**2*t3*t2**8*t5-325._ki/3._ki*t4*t3**3*t2**6*t1**3*t5*&
                  &*2+10._ki*t4*t3**2*t5**2*t2**8*t1**3-25._ki/2._ki*t1**4*t4*t5**4*t&
                  &3*t2**6-6._ki*t1**4*t3*t5**2*t4*t2**8+t1**4*t2**9*t4*t5*t3+10._ki&
                  &*t1**5*t4**4*t5**5*t2**2-6._ki*t1**4*t3**2*t2**8*t5-25._ki*t1**4*&
                  &t5**3*t2**6*t3**2-t1**4*t3*t4*t2**10/2._ki-30._ki*t1**4*t3**2*t2*&
                  &*7*t5**2-t1**3*t3**2*t4*t2**10/6._ki+8._ki*t1**4*t4**4*t3**2*t2**&
                  &5+31._ki/2._ki*t1**4*t4**2*t3**2*t2**7
                !
                stemp11=-48._ki*t1**4*t4*t3**3*t2**6-5._ki*t1**5*t4**3*t5*t2**7-278&
                  &._ki*t1*t2**4*t3**5*t4**5+240._ki*t1*t2**5*t3**5*t4**4+72._ki*t1*t&
                  &4**7*t3**4*t2**4-85._ki*t1*t4**4*t3**4*t2**7+140._ki*t1*t4**5*t3*&
                  &*4*t2**6+61._ki/2._ki*t1*t4**3*t3**4*t2**8+172._ki*t1*t2**3*t3**5*&
                  &t4**6+6._ki*t1**4*t4**4*t3*t2**7-13._ki/2._ki*t1**4*t4**3*t3*t2**8&
                  &-35._ki/2._ki*t1**4*t3**5*t2**2*t4-2._ki*t1**4*t4**5*t3*t2**6-19._k&
                  &i*t1**4*t4**3*t3**2*t2**6+19._ki/3._ki*t1**3*t4**4*t3**2*t2**7+70&
                  &._ki*t1**3*t4*t3**4*t2**6-145._ki/3._ki*t1**3*t3**5*t2**4*t4-136._k&
                  &i*t1*t4**6*t3**4*t2**5-44._ki*t1*t2**2*t3**5*t4**7
                !
                stemp9=stemp11+133._ki/4._ki*t1*t2**7*t3**5*t4**2-16._ki*t1*t4**8*t3&
                  &**4*t2**3+7._ki*t1**5*t4**3*t3*t2**6+20._ki*t1**5*t3**3*t2**4*t4+&
                  &40._ki/3._ki*t4*t3**4*t2**8*t1**2-20._ki*t1**5*t4**3*t5**3*t2**5+1&
                  &10._ki*t1**4*t4*t3**4*t2**4-75._ki*t1**4*t3**4*t5*t2**4-9._ki*t1**&
                  &5*t4**2*t3**2*t2**5+88._ki*t1**4*t4**2*t3**3*t2**5-50._ki*t1**4*t&
                  &3**3*t4**3*t2**4-2._ki*t1**5*t4**4*t5**2*t2**5+535._ki/6._ki*t1**3&
                  &*t3**5*t2**3*t4**2+755._ki/3._ki*t1**3*t4**3*t3**4*t2**4+4._ki/3._k&
                  &i*t1**3*t4**6*t3**2*t2**5-110._ki*t1**3*t4**4*t3**4*t2**3+38._ki*&
                  &t1**3*t4**2*t3**3*t2**7-105._ki/2._ki*t1**3*t3**5*t2**2*t4**3-25.&
                  &_ki/6._ki*t1**3*t4**3*t3**2*t2**8+stemp10
                !
                stemp11=-14._ki/3._ki*t1**3*t4**5*t3**2*t2**6-205._ki*t1**3*t4**2*t3&
                  &**4*t2**5-226._ki/3._ki*t1**3*t4**3*t3**3*t2**6+208._ki/3._ki*t1**3&
                  &*t4**4*t3**3*t2**5-24._ki*t1**3*t4**5*t3**3*t2**4+25._ki/3._ki*t3*&
                  &*3*t5**2*t2**7*t1**3-t4*t2**10*t1**2*t3**3/3._ki+8._ki*t1**5*t4**&
                  &4*t5**3*t2**4-26._ki/3._ki*t4*t2**8*t3**3*t1**3+t4*t3**4*t2**10*t&
                  &1/2._ki-95._ki*t1**4*t4**2*t3**4*t2**3+150._ki*t1**5*t3**2*t5**3*t&
                  &2**4+7._ki*t1**5*t4*t3*t2**8+40._ki*t1**5*t4**2*t5**4*t2**5+36._ki&
                  &*t1**5*t5*t3**2*t2**6+16._ki*t1**5*t4**2*t5**3*t2**6+100._ki*t1**&
                  &4*t3**3*t5**2*t2**5+10._ki/3._ki*t4**2*t3**3*t2**9*t1**2+80._ki*t1&
                  &**4*t3**3*t5*t2**6
                !
                stemp10=stemp11-5._ki*t1**4*t4*t3**2*t2**8+3._ki*t1**4*t4**2*t3*t2*&
                  &*9+950._ki*t1**4*t4**3*t3**2*t5**3*t2**3-80._ki*t1**4*t4**5*t5**3&
                  &*t3*t2**3-80._ki*t1**3*t4**6*t3**2*t5**2*t2**3-455._ki/6._ki*t6*t1&
                  &**3*t4**2*t3**3*t2**6+220._ki*t6*t1**3*t4**2*t3**4*t2**4+65._ki/3&
                  &._ki*t6*t1**3*t4**5*t3**2*t2**5-55._ki/3._ki*t6*t1**3*t4**4*t3**2*&
                  &t2**6+205._ki/24._ki*t6*t1**3*t4**3*t3**2*t2**7+57._ki/2._ki*t6*t1*&
                  &*5*t4**3*t3*t2**5+195._ki*t6*t1**5*t4**3*t5**4*t2**3+60._ki*t6*t1&
                  &**5*t4**5*t5**4*t2-1155._ki*t6**2*t1**5*t3**2*t2**3*t4*t5-1243._k&
                  &i/6._ki*t6**2*t1**5*t4**3*t3*t2**4+1815._ki*t6**2*t1**5*t4**4*t3*&
                  &*2*t5+209._ki/2._ki*t6**2*t1**5*t4**2*t3*t2**5+stemp8+stemp9+195.&
                  &_ki*t6*t1**5*t4**3*t5**3*t2**4
                !
                stemp11=stemp10-3._ki/4._ki*t6*t1**5*t4*t2**8*t5-t6*t1**4*t4*t3*t2*&
                  &*9/6._ki+45._ki*t6*t1**5*t4**4*t5**2*t2**4-15._ki*t6*t1**5*t4**5*t&
                  &5**2*t2**3+1140._ki*t6*t1**5*t4**2*t3**3*t5*t2-120._ki*t6*t1**5*t&
                  &4**4*t5*t2**3*t3+1395._ki*t6*t1**5*t4**2*t3*t5**2*t2**4+1125._ki*&
                  &t6*t1**5*t4**3*t3**2*t5*t2**2-465._ki/2._ki*t6*t1**5*t4**2*t5*t2*&
                  &*5*t3-2280._ki*t6*t1**5*t4**3*t3*t5**3*t2**2+9._ki*t6*t1**5*t4**4&
                  &*t5*t2**5-39._ki/4._ki*t6*t1**5*t4**3*t5*t2**6+60._ki*t6*t1**5*t4*&
                  &*5*t5**3*t2**2-180._ki*t6*t1**5*t4**4*t5**3*t2**3-3._ki*t6*t1**5*&
                  &t4**5*t5*t2**4-12._ki*t6*t1**5*t4**4*t3*t2**4-200._ki/3._ki*t1**3*&
                  &t4**6*t3**2*t5**3*t2**2-950._ki/3._ki*t1**3*t4**4*t3**2*t5**3*t2*&
                  &*4-275._ki*t1**3*t4**4*t3**4*t2**2*t5
                !
                stemp7=stemp11+2600._ki/3._ki*t1**3*t4**4*t3**3*t5**2*t2**3-300._ki*&
                  &t1**3*t4**5*t3**3*t5**2*t2**2+625._ki/3._ki*t1**3*t4**3*t3**2*t5*&
                  &*3*t2**5+475._ki*t1**3*t4**2*t3**3*t5**2*t2**5-t1**5*t3*t2**9+14&
                  &6._ki/3._ki*t4**4*t3**5*t2**7+8._ki*t1**4*t3**3*t2**7+160._ki/3._ki*&
                  &t2**3*t3**5*t4**8-12._ki*t1**5*t3**3*t2**5-30._ki*t1**4*t3**4*t2*&
                  &*5+25._ki/3._ki*t1**3*t3**5*t2**5-t3**5*t4*t2**10/6._ki+10._ki*t1**&
                  &5*t2**3*t3**4+t1**4*t3**2*t2**9/2._ki+7._ki/3._ki*t4**2*t3**5*t2**&
                  &9+200._ki*t1**5*t4*t5*t3**3*t2**3-128._ki/3._ki*t6*t1**2*t4**8*t3*&
                  &*3*t5*t2+1045._ki/3._ki*t6**2*t1**5*t4**4*t5**3*t2**2+210._ki*t6**&
                  &2*t1**3*t4**5*t3**2*t5*t2**3-615._ki/4._ki*t6**2*t1**4*t4**3*t5**&
                  &2*t2**4*t3
                !
                stemp11=stemp7+990._ki*t6**2*t1**5*t4**5*t3*t5**2-1520._ki/3._ki*t6*&
                  &t1**4*t4**4*t5**3*t3*t2**3+1140._ki*t6*t1**4*t4**2*t3**2*t5**2*t&
                  &2**4+40._ki/3._ki*t6*t1**4*t4**6*t3*t2**3*t5+1120._ki/3._ki*t6*t1**&
                  &4*t4**5*t5**3*t3*t2**2-880._ki*t6*t1**4*t4**4*t3**3*t5*t2+1180._k&
                  &i/3._ki*t6*t1**3*t4**4*t3**4*t2**2-272._ki/3._ki*t6*t1**2*t4**6*t3&
                  &**3*t2**4-32._ki/3._ki*t6*t1**2*t4**8*t3**3*t2**2+3._ki/4._ki*t6**2&
                  &*t3**2*t4*t2**7*t1**3*t5-15._ki/4._ki*t6**2*t3*t4*t2**6*t1**4*t5*&
                  &*2+7._ki/6._ki*t6**2*t4**2*t3**3*t2**7*t1**2+10._ki*t6*t5**2*t4*t2&
                  &**7*t1**4*t3-5._ki/3._ki*t6*t4*t5*t3*t2**8*t1**4+t6*t4*t3**3*t2**&
                  &9*t1**2/3._ki+5._ki/24._ki*t6*t4*t2**9*t1**3*t3**2-78._ki*t1**4*t4*&
                  &*3*t3*t2**6*t5**2-24._ki*t1**4*t4**5*t3*t2**4*t5**2
                !
                stemp10=stemp11-480._ki*t1**4*t4*t3**3*t5*t2**5-170._ki/3._ki*t6*t1*&
                  &*2*t4**4*t3**3*t2**6-14._ki*t1**5*t4**3*t5*t3*t2**5-1088._ki/3._ki&
                  &*t6*t1**2*t4**6*t3**3*t5*t2**3+192._ki*t6*t1**2*t4**7*t3**3*t5*t&
                  &2**2+26._ki*t1**5*t4**2*t5*t3*t2**6-520._ki*t1**5*t4**2*t3*t5**3*&
                  &t2**4-900._ki*t6**2*t1**6*t4*t3*t5**2*t2**2+350._ki*t1*t4**5*t3**&
                  &4*t5*t2**5-40._ki*t1*t4**8*t3**4*t5*t2**2+180._ki*t1*t4**7*t3**4*&
                  &t5*t2**3+305._ki/4._ki*t1*t4**3*t5*t2**7*t3**4-425._ki/2._ki*t1*t4*&
                  &*4*t5*t2**6*t3**4+50._ki*t6*t4**2*t5**2*t2**6*t1**3*t3**2-55._ki/&
                  &6._ki*t6**2*t4*t5**3*t2**5*t1**5+t6**2*t4*t2**8*t1**4*t3/4._ki-55&
                  &._ki/2._ki*t6**2*t3*t5**2*t2**5*t1**5+11._ki/6._ki*t6**2*t4*t5*t2**&
                  &7*t1**5-360._ki*t6*t1**5*t3**2*t2**4*t5**2
                !
                stemp11=stemp10-90._ki*t6*t1**5*t4**2*t5**4*t2**4-200._ki/3._ki*t6*t&
                  &5*t2**5*t3**3*t1**4+40._ki/3._ki*t6*t4*t3**3*t2**7*t1**3+7._ki/2._k&
                  &i*t6*t1**6*t2**6*t4*t5-90._ki*t6*t1**6*t3*t5**2*t2**4+35._ki/2._ki&
                  &*t6*t1**6*t4*t5**2*t2**5-70._ki*t6*t1**6*t4*t5**4*t2**3-120._ki*t&
                  &6*t1**6*t3*t5**3*t2**3+15._ki*t6*t1**6*t3*t5*t2**5-1025._ki/6._ki*&
                  &t1**2*t4**3*t5**2*t2**6*t3**3+700._ki/3._ki*t1**3*t4**5*t3**2*t5*&
                  &*3*t2**3+440._ki*t6*t1**6*t4*t5**3*t3*t2**2-1375._ki/6._ki*t6**2*t&
                  &1**5*t4**3*t5**3*t2**3-3960._ki*t6*t1**5*t4**2*t3**2*t5**2*t2**2&
                  &-410._ki*t6*t1**4*t2**4*t3**3*t4**2-1300._ki/3._ki*t1**2*t4**5*t3*&
                  &*3*t5**2*t2**4-16._ki*t6*t4**2*t5*t2**7*t1**2*t3**3+280._ki/3._ki*&
                  &t6*t1**2*t4**5*t3**3*t2**5
                !
                stemp9=stemp11-176._ki/3._ki*t6*t1**2*t4**7*t3**4*t2-1112._ki/3._ki*t&
                  &6*t1**2*t4**5*t3**4*t2**3+133._ki/3._ki*t6*t1**2*t4**2*t3**4*t2**&
                  &6+320._ki*t6*t1**2*t2**4*t4**4*t3**4+61._ki/3._ki*t6*t1**2*t4**3*t&
                  &3**3*t2**7+48._ki*t6*t1**2*t4**7*t3**3*t2**3+688._ki/3._ki*t6*t1**&
                  &2*t4**6*t3**4*t2**2+160._ki/3._ki*t6*t1*t2**2*t3**4*t4**8+244._ki/&
                  &3._ki*t6*t1**2*t4**3*t5*t2**6*t3**3+7._ki/3._ki*t6*t4**2*t3**4*t2*&
                  &*8*t1+15._ki*t6*t1**5*t4*t5**3*t2**6-15._ki/4._ki*t6*t1**5*t4*t5**&
                  &2*t2**7-15._ki/2._ki*t6*t1**5*t5*t2**7*t3-15._ki/4._ki*t6**2*t2**6*&
                  &t3**2*t1**4*t5+150._ki*t6**2*t1**6*t3*t5**2*t2**3+125._ki*t6**2*t&
                  &1**6*t4*t5**3*t2**3-25._ki*t6**2*t1**6*t4*t5*t2**5+275._ki/2._ki*t&
                  &6**2*t5*t2**4*t3**2*t1**5+2._ki*t4*t3**2*t2**9*t1**3*t5+1155._ki/&
                  &2._ki*t6**2*t1**5*t4**3*t3**3
                !
                stemp11=stemp9+50._ki*t6**2*t1**6*t3**3*t2-275._ki/3._ki*t6**2*t1**5&
                  &*t2**3*t3**3+220._ki/3._ki*t6**2*t1**5*t4**6*t5**3-210._ki*t6**2*t&
                  &1**4*t4**5*t3**3-44._ki*t6**2*t1**3*t4**7*t3**3-16._ki/3._ki*t6**2&
                  &*t1**2*t4**9*t3**3-10._ki*t6*t1**6*t3**3*t2**2+3._ki/2._ki*t6*t1**&
                  &5*t4**5*t2**5-9._ki/2._ki*t6*t1**5*t4**4*t2**6+39._ki/8._ki*t6*t1**&
                  &5*t4**3*t2**7-60._ki*t6*t1**5*t3**4*t2**2+90._ki*t6*t1**5*t3**3*t&
                  &2**4-7._ki/4._ki*t6*t1**6*t4**3*t2**5+15._ki*t6**2*t1**7*t3**2*t5-&
                  &200._ki*t6**2*t1**6*t4**4*t5**3-9._ki/4._ki*t6*t1**5*t4**2*t2**8+1&
                  &5._ki*t6*t1**5*t3**2*t2**6+5._ki*t6*t3**4*t2**6*t1**3
                !
                stemp10=stemp11-50._ki/3._ki*t6*t3**3*t2**6*t1**4+5._ki/2._ki*t6**2*t&
                  &1**6*t5*t2**6+t6*t2**8*t1**2*t3**4/3._ki-5._ki/6._ki*t6*t2**8*t3**&
                  &2*t1**4-5._ki/2._ki*t6*t1**6*t5**2*t2**6-t6*t1**6*t5*t2**7/2._ki+1&
                  &0._ki*t6*t1**6*t5**4*t2**4+10._ki*t6*t1**6*t5**3*t2**5+3._ki/8._ki*&
                  &t6*t1**5*t4*t2**9-3._ki/4._ki*t6*t2**8*t3*t1**5-5._ki/6._ki*t6*t3**&
                  &3*t2**8*t1**3+15._ki*t6**2*t1**7*t5**3*t2**2-3._ki*t6**2*t1**7*t5&
                  &*t2**4-10._ki*t6**2*t1**6*t3*t2**5+15._ki/2._ki*t6**2*t3**3*t2**5*&
                  &t1**4+3._ki/2._ki*t6*t1**6*t3*t2**6-7._ki/4._ki*t6*t1**6*t4*t2**7+t&
                  &6**2*t2**7*t3**3*t1**3/4._ki+11._ki/6._ki*t6**2*t2**7*t3*t1**5-25.&
                  &_ki/2._ki*t6**2*t1**6*t5**3*t2**4
                !
                stemp11=stemp10+45._ki*t6**2*t1**7*t4**2*t5**3-10._ki*t1**5*t4*t5**&
                  &4*t2**6-t1**5*t4*t2**9*t5-4._ki*t1**5*t4*t5**3*t2**7+t1**5*t2**8&
                  &*t4*t5**2-5._ki*t1**5*t4*t5**5*t2**5-12._ki*t1**5*t3*t5**2*t2**7-&
                  &40._ki*t1**5*t3*t5**3*t2**6-25._ki*t1**5*t5**4*t2**5*t3+2._ki*t1**&
                  &5*t5*t3*t2**8-25._ki/12._ki*t2**8*t3**4*t1**2*t5-19._ki/4._ki*t4*t3&
                  &**5*t2**8*t1-125._ki/6._ki*t5*t2**6*t3**4*t1**3+20._ki*t1**5*t4**2&
                  &*t5**5*t2**4+4._ki*t1**5*t4**2*t2**8*t5+180._ki*t1**5*t3**2*t2**5&
                  &*t5**2-4._ki*t1**5*t4**2*t5**2*t2**7-110._ki*t6**2*t1**6*t2**3*t3&
                  &*t4**2+1100._ki/3._ki*t1**2*t4**4*t3**3*t5**2*t2**5
                !
                stemp8=stemp11-1380._ki*t6**2*t1**4*t4**4*t3**2*t5*t2**2+330._ki*t6&
                  &**2*t1**4*t4**2*t3**3*t2**3+26._ki*t6**2*t1**4*t4**5*t3*t2**4+4.&
                  &_ki*t6**2*t1**4*t4**7*t3*t2**2-180._ki*t6*t1**5*t4**4*t5**4*t2**2&
                  &-330._ki*t6*t1**5*t4*t3**3*t2**3+1070._ki/3._ki*t6*t1**4*t3**4*t2*&
                  &*2*t4**2+30._ki*t6*t1**4*t4**5*t3**2*t2**3+565._ki/6._ki*t6*t1**4*&
                  &t4**3*t3**2*t2**5-95._ki/2._ki*t6*t1**4*t4**2*t3**2*t2**6+140._ki*&
                  &t6*t1**4*t2**5*t3**3*t4-260._ki/3._ki*t6*t1**4*t4**4*t3**2*t2**4-&
                  &580._ki/3._ki*t6*t1**4*t3**4*t2**3*t4-300._ki*t6*t1**6*t4*t3**2*t5&
                  &**2*t2+1080._ki*t6*t1**5*t4*t3**2*t5*t2**4-3680._ki/3._ki*t6*t1**3&
                  &*t4**4*t3**3*t5*t2**3-475._ki/3._ki*t6*t1**2*t3**4*t4**3*t2**5-25&
                  &._ki/6._ki*t6*t1**4*t4**3*t3*t2**7-220._ki*t6*t1**4*t4**4*t3**3*t2&
                  &**2+4._ki/3._ki*t6*t1**4*t4**6*t3*t2**4
                !
                stemp11=stemp8+130._ki*t6*t1**6*t4**2*t5**4*t2**2+7._ki/2._ki*t6*t1*&
                  &*6*t4**3*t5*t2**4-40._ki*t6*t1**6*t3**3*t5*t2-13._ki/2._ki*t6*t1**&
                  &6*t4**2*t5*t2**5+130._ki*t6*t1**6*t4**2*t5**3*t2**3-70._ki*t6*t1*&
                  &*6*t4**3*t5**4*t2+35._ki/2._ki*t6*t1**6*t4**3*t5**2*t2**3-70._ki*t&
                  &6*t1**6*t4**3*t5**3*t2**2+285._ki*t6*t1**5*t3**3*t2**2*t4**2+105&
                  &._ki*t6*t1**5*t3**4*t2*t4-90._ki*t6*t1**5*t4*t3**2*t2**5+165._ki*t&
                  &6*t1**5*t4**2*t3**2*t2**4-5._ki/6._ki*t3**4*t2**9*t1**2-14._ki*t1*&
                  &*5*t4*t3*t2**7*t5+280._ki*t1**5*t4*t5**3*t2**5*t3+84._ki*t1**5*t4&
                  &*t5**2*t2**6*t3+300._ki*t1**4*t4*t3**2*t5**2*t2**6+250._ki*t1**4*&
                  &t4*t5**3*t2**5*t3**2
                !
                stemp10=stemp11+120._ki*t1**4*t4**2*t3*t5**3*t2**6+60._ki*t1**4*t3*&
                  &*2*t5*t4*t2**7+36._ki*t1**4*t4**2*t5**2*t2**7*t3+75._ki*t1**4*t4*&
                  &*2*t3*t5**4*t2**5+4._ki/3._ki*t1**3*t4**2*t3**2*t2**9+20._ki/3._ki*&
                  &t5*t3**3*t2**8*t1**3+13._ki/4._ki*t6*t1**6*t4**2*t2**6+100._ki/3._k&
                  &i*t6*t1**4*t3**4*t2**4-15._ki/2._ki*t6*t1**6*t3**2*t2**4+3._ki*t6*&
                  &*2*t1**7*t3*t2**3-175._ki/2._ki*t6**2*t1**6*t3**3*t4-140._ki/3._ki*&
                  &t6*t1**4*t4**5*t3*t2**4*t5+190._ki/3._ki*t6*t1**4*t4**4*t3*t2**5*&
                  &t5-2260._ki*t6*t1**4*t4**3*t3**2*t5**2*t2**3-14._ki/3._ki*t6*t1**4&
                  &*t4**5*t3*t2**5-325._ki/4._ki*t6**2*t1**4*t3**3*t2**4*t4-16._ki*t6&
                  &**2*t1**4*t4**6*t3*t2**3-22._ki*t6**2*t1**4*t4**4*t3*t2**5+590._k&
                  &i*t6**2*t1**4*t4**4*t3**3*t2
                !
                stemp11=stemp10+240._ki*t6**2*t1**3*t4**4*t3**3*t2**3+133._ki/4._ki*&
                  &t6**2*t1**3*t4**2*t3**3*t2**5+172._ki*t6**2*t1**3*t4**6*t3**3*t2&
                  &-475._ki/4._ki*t6**2*t1**3*t3**3*t4**3*t2**4+945._ki*t6**2*t1**4*t&
                  &4**3*t5*t3**2*t2**3+330._ki*t6**2*t1**4*t4**4*t5**2*t2**3*t3-390&
                  &._ki*t6**2*t1**4*t4**5*t3*t5**2*t2**2-95._ki*t6**2*t1**6*t4**3*t5&
                  &*t2**3-1425._ki/2._ki*t6**2*t1**6*t4**2*t3**2*t5-1875._ki/2._ki*t6*&
                  &*2*t1**6*t4**3*t3*t5**2-225._ki*t6**2*t1**6*t3**2*t5*t2**2+155._k&
                  &i/2._ki*t6**2*t1**6*t4**2*t5*t2**4+60._ki*t6**2*t1**6*t2**4*t3*t4&
                  &+40._ki*t6**2*t1**6*t4**4*t5*t2**2-775._ki/2._ki*t6**2*t1**6*t4**2&
                  &*t5**3*t2**2+475._ki*t6**2*t1**6*t4**3*t5**3*t2+125._ki/2._ki*t6**&
                  &2*t1**6*t2**2*t3*t4**3-209._ki/3._ki*t6**2*t1**5*t4**4*t5*t2**4+4&
                  &0._ki/3._ki*t6*t4*t5**3*t2**6*t1**4*t3
                !
                stemp9=stemp11+2160._ki*t6*t1**5*t4*t3**2*t5**2*t2**3-1320._ki*t6*t&
                  &1**5*t3**3*t5*t2**2*t4+285._ki*t6*t1**5*t4**3*t5*t2**4*t3+1860._k&
                  &i*t6*t1**5*t4**2*t5**3*t3*t2**3+960._ki*t6*t1**5*t4**4*t3*t5**3*&
                  &t2+560._ki*t6*t1**4*t4*t3**3*t5*t2**4-80._ki*t6*t1**4*t4**6*t3*t5&
                  &**2*t2**2+250._ki*t6*t1**4*t4**3*t3*t5**2*t2**5+1510._ki/3._ki*t6*&
                  &t1**4*t2**3*t3**3*t4**3-20._ki*t1**4*t4*t5**3*t3*t2**7+2080._ki*t&
                  &6*t1**4*t4**4*t3**2*t5**2*t2**2-720._ki*t6*t1**4*t4**5*t3**2*t5*&
                  &*2*t2-520._ki*t6*t1**3*t4**5*t3**2*t5**2*t2**3-660._ki*t1**5*t4*t&
                  &3**2*t5**2*t2**4+t2**9*t3**5*t1/4._ki+5._ki/2._ki*t3**5*t2**7*t1**&
                  &2+t6*t1**6*t2**8/4._ki-32._ki/3._ki*t2**2*t3**5*t4**9-3._ki*t1**5*t&
                  &3**2*t2**7+10._ki*t1**4*t3**5*t2**3
                !
                stemp11=stemp9-85._ki/6._ki*t4**3*t3**5*t2**8-344._ki/3._ki*t2**4*t3*&
                  &*5*t4**7-310._ki/3._ki*t4**5*t3**5*t2**6-25._ki/3._ki*t3**4*t2**7*t&
                  &1**3-70._ki*t6*t1**6*t4*t5**3*t2**4+40._ki/3._ki*t6*t4**2*t2**7*t1&
                  &**4*t3*t5-260._ki*t6*t4*t3**2*t2**5*t1**4*t5**2-80._ki*t6*t4**2*t&
                  &3*t2**6*t1**4*t5**2+160._ki/3._ki*t6*t4*t5*t2**6*t3**3*t1**3-5._ki&
                  &/2._ki*t6**2*t4**2*t3*t2**7*t1**4+60._ki*t6**2*t4*t5*t2**5*t3**2*&
                  &t1**4-4._ki*t6*t4**2*t3**3*t2**8*t1**2+65._ki/6._ki*t6*t4*t2**7*t3&
                  &**2*t1**4+75._ki*t6*t1**5*t4*t5*t2**6*t3-10._ki/3._ki*t4*t5*t3**3*&
                  &t2**9*t1**2-80._ki*t4**2*t3**2*t2**7*t1**3*t5**2-25._ki/6._ki*t3**&
                  &3*t4*t2**8*t1**2*t5**2+5._ki/4._ki*t3**4*t4*t2**9*t1*t5
                !
                stemp10=stemp11+175._ki*t1**5*t4*t5**4*t2**4*t3-5885._ki/6._ki*t6**2&
                  &*t1**5*t3**3*t4**2*t2-44._ki/3._ki*t6**2*t1**5*t4**6*t5*t2**2+154&
                  &._ki/3._ki*t6**2*t1**5*t4**5*t5*t2**3-66._ki*t6**2*t1**5*t4**5*t3*&
                  &t2**2-300._ki*t6**2*t1**4*t4**6*t3**2*t5+41._ki/4._ki*t6**2*t1**4*&
                  &t4**3*t3*t2**6-1275._ki/2._ki*t6**2*t1**4*t4**3*t3**3*t2**2+6215.&
                  &_ki/2._ki*t6**2*t1**5*t4**3*t3*t5**2*t2**2-2860._ki*t6**2*t1**5*t4&
                  &**4*t3*t5**2*t2-3135._ki/2._ki*t6**2*t1**5*t4**2*t3*t5**2*t2**3+6&
                  &765._ki/2._ki*t6**2*t1**5*t4**2*t3**2*t5*t2**2-1365._ki/4._ki*t6**2&
                  &*t1**4*t4**2*t3**2*t5*t2**4+1020._ki*t6**2*t1**4*t4**5*t3**2*t5*&
                  &t2-156._ki*t1**5*t4**2*t5**2*t2**5*t3-920._ki/3._ki*t6*t1**3*t4**4&
                  &*t3**3*t2**4-344._ki/3._ki*t6*t1*t2**3*t3**4*t4**7-32._ki/3._ki*t6*&
                  &t1*t2*t3**4*t4**9-125._ki/3._ki*t6*t1**4*t4**3*t3*t2**6*t5+6040._k&
                  &i/3._ki*t6*t1**4*t4**3*t3**3*t2**2*t5
                !
                stemp11=stemp10-380._ki*t6*t1**4*t4**4*t3*t5**2*t2**4+570._ki*t6*t1&
                  &**4*t4**2*t3**2*t5*t2**5+280._ki*t6*t1**4*t4**5*t3*t5**2*t2**3-3&
                  &60._ki*t6*t1**6*t4**2*t5**3*t2*t3+45._ki*t6*t1**6*t4**2*t5*t2**3*&
                  &t3+180._ki*t6*t1**6*t3**2*t5**2*t2**2+90._ki*t6*t1**6*t2**3*t3**2&
                  &*t5+25._ki/2._ki*t6*t1**6*t2**3*t4*t3**2-65._ki/2._ki*t6*t1**6*t4**&
                  &2*t5**2*t2**4-11._ki/2._ki*t6*t1**6*t4*t3*t2**5-50._ki*t1**4*t4**5&
                  &*t5**4*t3*t2**2+1140._ki*t1**4*t4**3*t3**2*t5**2*t2**4+72._ki*t1*&
                  &*4*t4**4*t3*t2**5*t5**2-625._ki*t1**4*t4**3*t3**3*t5**2*t2**2+15&
                  &0._ki*t1**4*t4**4*t5**4*t3*t2**3-400._ki*t1**4*t4**4*t3**2*t5**3*&
                  &t2**2-186._ki*t1**4*t4**2*t5*t2**6*t3**2-480._ki*t1**4*t4**4*t3**&
                  &2*t5**2*t2**3-200._ki/3._ki*t4**2*t5**3*t2**6*t1**3*t3**2
                !
                stemp6=stemp11-500._ki/3._ki*t1**2*t4**6*t3**4*t5*t2**2-410._ki/3._ki&
                  &*t1**2*t4**3*t5*t3**3*t2**7+2._ki/3._ki*t3**3*t2**9*t1**3-2300._ki&
                  &/3._ki*t1**2*t4**4*t3**4*t5*t2**4+640._ki/3._ki*t1**2*t4**6*t3**3*&
                  &t5*t2**4-205._ki/2._ki*t6*t1**3*t4**3*t5*t3**2*t2**6+1700._ki/3._ki&
                  &*t1**2*t4**5*t3**4*t5*t2**3+416._ki/3._ki*t2**5*t3**5*t4**6+525._k&
                  &i*t1**2*t4**3*t3**4*t5*t2**5+800._ki/3._ki*t1**2*t4**6*t3**3*t5**&
                  &2*t2**3-260._ki/3._ki*t4*t3**3*t2**7*t1**3*t5-15._ki*t4**2*t5*t2**&
                  &8*t1*t3**4+100._ki/3._ki*t4**2*t5*t3**3*t2**8*t1**2+1100._ki*t1**4&
                  &*t4**2*t3**3*t2**3*t5**2-475._ki/2._ki*t1**4*t4**2*t3**4*t5*t2**2&
                  &-930._ki*t1**4*t4**2*t3**2*t5**2*t2**5-260._ki*t1**4*t4**3*t5**3*&
                  &t3*t2**5+240._ki*t1**4*t4**4*t5**3*t3*t2**4-160._ki/3._ki*t1**2*t4&
                  &**7*t3**3*t5*t2**3-60._ki*t6**2*t1**4*t4**7*t3*t5**2
                !
                stemp7=t6/t1**5/t2**12
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(4)
                !
                stemp5=95._ki/3._ki*t6**3*t3**2*t4**3-5._ki/2._ki*t6**2*t2**2*t3**3+t&
                  &6*t4*t2**8/4._ki+20._ki/3._ki*t6**3*t4**5*t5**2-7._ki/4._ki*t6**3*t4&
                  &**2*t2**5-3._ki/2._ki*t6**3*t1**2*t2**3+2._ki*t6**3*t4**5*t2**2+2.&
                  &_ki*t6**2*t4**4*t2**4-5._ki*t6**3*t4**4*t2**3-3._ki/4._ki*t6*t4**2*&
                  &t2**7+t6*t4**3*t2**6/2._ki+7._ki/12._ki*t6**2*t3*t2**6-t6**2*t4*t2&
                  &**7/2._ki-t6**2*t1*t2**6+t6*t1*t2**7/4._ki+t4*t2**8*t5/12._ki-t4*t&
                  &3*t2**7/6._ki-7._ki/12._ki*t6**2*t1**2*t5*t2**3+5._ki/6._ki*t6**2*t1&
                  &**2*t2*t5**3-7._ki/3._ki*t6**2*t1*t2**4*t3+5._ki*t6**2*t2*t3**3*t4&
                  &+15._ki/2._ki*t6**2*t2**3*t3**2*t5+15._ki/2._ki*t6**3*t1*t3**2*t2-1&
                  &5._ki*t6**3*t3**2*t1*t4-95._ki/2._ki*t6**3*t3**2*t2*t4**2+45._ki/2.&
                  &_ki*t6**3*t2**2*t3**2*t4+27._ki/2._ki*t6**3*t1*t4**2*t2**3
                !
                stemp4=-6._ki*t6**3*t1*t4*t2**4-4._ki*t6**2*t4**3*t2**5+t6**3*t4*t2&
                  &**6/4._ki-10._ki/3._ki*t6**3*t2**3*t3**2+16._ki*t6**3*t1*t4*t5*t2**&
                  &3+65._ki*t6**3*t4**2*t3*t5*t2**2-280._ki/3._ki*t6**3*t4**3*t3*t5*t&
                  &2-20._ki*t6**3*t1*t4*t5**2*t2**2+65._ki*t6**3*t1*t3*t5*t2*t4+24._k&
                  &i*t6**3*t1*t4**3*t5*t2+26._ki*t6**3*t1*t2*t3*t4**2-20._ki*t6**2*t&
                  &1*t2*t3*t4*t5**2-5._ki*t6**2*t1*t2*t5*t3**2+14._ki/3._ki*t6**2*t1*&
                  &t2**3*t3*t4+10._ki*t6**2*t1*t5**2*t3*t2**2+65._ki/2._ki*t6**2*t2*t&
                  &3**2*t4**2*t5-8._ki*t6**3*t1**2*t5*t2*t4+35._ki/6._ki*t6**2*t1*t4*&
                  &*2*t5*t2**3-35._ki/6._ki*t6**2*t1*t4*t5*t2**4-45._ki*t6**2*t4**2*t&
                  &5**2*t3*t2**2+20._ki*t6**2*t4*t5**2*t3*t2**3+25._ki/3._ki*t6**2*t1&
                  &*t4*t5**3*t2**2-25._ki/3._ki*t6**2*t1*t4**2*t5**3*t2-65._ki/2._ki*t&
                  &6**2*t3**2*t4*t5*t2**2-55._ki/3._ki*t6**3*t4*t3*t5*t2**3-15._ki*t6&
                  &**3*t1*t3*t5*t2**2-65._ki*t6**3*t1*t4**2*t3*t5+stemp5
                !
                stemp5=-5._ki/2._ki*t6**2*t3*t5**2*t2**4-20._ki/3._ki*t6**2*t4**3*t5*&
                  &*3*t2**2+10._ki/3._ki*t6**2*t4**4*t5**3*t2+25._ki/6._ki*t6**2*t4**2&
                  &*t5**3*t2**3-2._ki/3._ki*t6**3*t4*t5*t2**5-5._ki*t6**2*t1*t4**2*t2&
                  &**4+5._ki*t6**2*t1*t4*t2**5-t6*t4*t2**6*t1/2._ki-30._ki*t6**3*t1*t&
                  &4**3*t5**2+6._ki*t6**3*t1*t2**3*t3+5._ki*t6**3*t1**2*t3*t5+3._ki*t&
                  &6**3*t1**2*t2**2*t4-56._ki/3._ki*t6**3*t4**4*t3*t2-26._ki*t6**3*t4&
                  &**2*t3*t2**3-50._ki/3._ki*t6**3*t4**4*t5**2*t2+112._ki/3._ki*t6**3*&
                  &t4**3*t3*t2**2+22._ki/3._ki*t6**3*t4*t3*t2**4+15._ki*t6**3*t4**3*t&
                  &5**2*t2**2-9._ki*t6**3*t1*t4**3*t2**2-16._ki/3._ki*t6**3*t4**5*t5*&
                  &t2-35._ki/6._ki*t6**3*t4**2*t5**2*t2**3-12._ki*t6**3*t4**3*t5*t2**&
                  &3+40._ki/3._ki*t6**3*t4**4*t5*t2**2+14._ki/3._ki*t6**3*t4**2*t5*t2*&
                  &*4+140._ki/3._ki*t6**3*t4**4*t3*t5+10._ki*t6**3*t1**2*t5**2*t4-5._k&
                  &i*t6**3*t1**2*t5**2*t2
                !
                stemp3=-2._ki*t6**3*t1**2*t2*t3+5._ki/3._ki*t6**3*t3*t2**4*t5+5._ki/6&
                  &._ki*t6**3*t4*t5**2*t2**4+4._ki*t6**3*t1**2*t5*t2**2+5._ki/2._ki*t6&
                  &**3*t1*t5**2*t2**3+7._ki/12._ki*t6**2*t2**6*t4*t5-t4**2*t5*t2**7/&
                  &12._ki+3._ki/4._ki*t6**3*t1*t2**5+9._ki/2._ki*t6**3*t4**3*t2**4+5._ki&
                  &/2._ki*t6**2*t4**2*t2**6+t1**2*t2**4*t6**2/2._ki-2._ki/3._ki*t6**3*&
                  &t2**5*t3+30._ki*t6**2*t4**3*t5**2*t3*t2-5._ki/6._ki*t6**2*t4*t5**3&
                  &*t2**4-5._ki/3._ki*t6**2*t1*t5**3*t2**3+7._ki/6._ki*t6**2*t1*t5*t2*&
                  &*5+stemp4+stemp5+t3*t2**8/12._ki+45._ki*t6**3*t1*t4**2*t5**2*t2-2&
                  &6._ki*t6**3*t1*t2**2*t3*t4-36._ki*t6**3*t1*t4**2*t5*t2**2+21._ki/2&
                  &._ki*t6**2*t4**2*t3*t2**4-35._ki/12._ki*t6**2*t4**2*t5*t2**5+14._ki&
                  &/3._ki*t6**2*t4**3*t5*t2**4-7._ki/3._ki*t6**2*t4**4*t5*t2**3-7._ki*&
                  &t6**2*t4**3*t3*t2**3-14._ki/3._ki*t6**2*t4*t3*t2**5-2._ki*t6**3*t1&
                  &*t5*t2**4
                !
                stemp4=1._ki/t2**10*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=-(t4*t3*t5*t2**2-4._ki*t4*t2*t1*t6*t5-t3*t2*t1*t6+3._ki*t4*t&
                  &3*t6*t2**2-9._ki*t4**2*t3*t6*t2+4._ki*t4**2*t5*t1*t6-2._ki*t4*t3**&
                  &2*t2+2._ki*t4*t3*t1*t6+t3**2*t2**2+6._ki*t3*t4**3*t6-t3*t4**2*t2*&
                  &t5)/t3/t2**4*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki
                !
                stemp9=195._ki/4._ki*t6**2*t1**3*t3**2*t5*t2*t4**2-195._ki/4._ki*t6**&
                  &2*t1**3*t3**2*t5*t2**2*t4-21._ki/2._ki*t6**2*t1**3*t4**3*t3*t2**3&
                  &-7._ki*t6**2*t1**3*t4*t3*t2**5+25._ki/6._ki*t6**2*t1**4*t4*t5**3*t&
                  &2**2+7._ki/3._ki*t6**2*t1**4*t2**3*t3*t4-10._ki*t6**2*t1**4*t5**2*&
                  &t3*t2*t4-5._ki/2._ki*t6**2*t1**4*t3**2*t5*t2-275._ki/9._ki*t6**3*t1&
                  &**3*t4**4*t5**2*t2+55._ki/2._ki*t6**3*t1**3*t4**3*t5**2*t2**2+25.&
                  &_ki/4._ki*t6**3*t1**4*t3**2*t2+10._ki*t6**2*t1*t4**5*t3**2*t5*t2**&
                  &2-10._ki/3._ki*t6**2*t1*t4**6*t3**2*t5*t2-35._ki/3._ki*t6**2*t1*t4*&
                  &*4*t3**2*t5*t2**3-95._ki/9._ki*t6**2*t1*t2**3*t3**3*t4**3+5._ki*t6&
                  &**2*t1**4*t5**2*t3*t2**2+11._ki/24._ki*t6**3*t1*t4**2*t3**2*t2**5&
                  &+55._ki/18._ki*t6**3*t1**3*t3*t2**4*t5+40._ki/3._ki*t6**3*t1**4*t4*&
                  &t5*t2**3-308._ki/9._ki*t6**3*t1**3*t4**4*t3*t2
                !
                stemp8=stemp9-t6**3*t1**5*t2**3/2._ki+5._ki/8._ki*t6**3*t1**4*t2**5+&
                  &715._ki/6._ki*t6**3*t1**3*t4**2*t3*t5*t2**2+45._ki/4._ki*t6**2*t1**&
                  &3*t2**3*t3**2*t5+63._ki/4._ki*t6**2*t1**3*t4**2*t3*t2**4-35._ki/8.&
                  &_ki*t6**2*t1**3*t4**2*t5*t2**5+7._ki*t6**2*t1**3*t4**3*t5*t2**4-7&
                  &._ki/2._ki*t6**2*t1**3*t4**4*t5*t2**3+20._ki/3._ki*t6**2*t1*t4**3*t&
                  &3**2*t5*t2**4+t1**3*t3*t2**8/36._ki+15._ki/2._ki*t6**2*t1**3*t3**3&
                  &*t2*t4-15._ki*t6**2*t1**2*t4**3*t3*t5**2*t2**3+20._ki*t6**3*t1**4&
                  &*t4**3*t5*t2+65._ki/3._ki*t6**3*t1**4*t2*t3*t4**2-22._ki*t6**3*t1*&
                  &*3*t4**3*t5*t2**3-50._ki/3._ki*t6**3*t1**4*t4*t5**2*t2**2+5._ki/24&
                  &._ki*t6**2*t1*t4*t3**2*t2**6*t5+7._ki/36._ki*t6**2*t1**2*t4*t3*t2*&
                  &*7-5._ki/6._ki*t6**2*t1**2*t2**5*t3**2*t5+3._ki/2._ki*t6**3*t1**2*t&
                  &4**2*t3*t2**5+55._ki/6._ki*t6**2*t1**2*t4*t3**2*t5*t2**4
                !
                stemp9=stemp8-15._ki/2._ki*t6**2*t1**2*t4*t3**3*t2**3-5._ki/6._ki*t6*&
                  &*2*t1**2*t3*t5**2*t2**5*t4+5._ki/12._ki*t6**3*t1**2*t4*t5*t2**5*t&
                  &3-t6**3*t1**2*t4*t3*t2**6/6._ki-t6**3*t1*t4*t3**2*t2**6/24._ki-49&
                  &._ki/36._ki*t6**2*t1**2*t4**2*t3*t2**6-35._ki/9._ki*t6**2*t1**2*t4*&
                  &*4*t3*t2**4+7._ki/2._ki*t6**2*t1**2*t4**3*t3*t2**5-25._ki/12._ki*t6&
                  &**3*t1*t4**3*t3**2*t2**4-20._ki/3._ki*t6**2*t1**2*t4**5*t3*t5**2*&
                  &t2-15._ki/8._ki*t6**2*t1*t4**2*t3**2*t2**5*t5-88._ki/9._ki*t6**3*t1&
                  &**3*t4**5*t5*t2-70._ki/3._ki*t6**2*t1**2*t3**2*t4**4*t5*t2+14._ki/&
                  &9._ki*t6**2*t1**2*t4**5*t3*t2**3+140._ki/3._ki*t6**2*t1**2*t3**2*t&
                  &4**3*t5*t2**2-65._ki/2._ki*t6**2*t1**2*t3**2*t4**2*t5*t2**3+20._ki&
                  &*t6**3*t1**2*t4**5*t3*t5*t2-70._ki/3._ki*t6**3*t1**2*t4**4*t3*t5*&
                  &t2**2-16._ki/3._ki*t6**3*t1**2*t4**3*t3*t2**4+85._ki/6._ki*t6**3*t1&
                  &**2*t4**2*t3**2*t2**3
                !
                stemp7=stemp9-20._ki/3._ki*t6**3*t1**2*t4**6*t3*t5+8._ki/3._ki*t6**3*&
                  &t1**2*t4**6*t3*t2+28._ki/3._ki*t6**3*t1**2*t4**4*t3*t2**3-325._ki/&
                  &6._ki*t6**3*t1**4*t4**2*t3*t5+325._ki/6._ki*t6**3*t1**4*t3*t5*t2*t&
                  &4+40._ki/3._ki*t6**3*t1**2*t4**3*t3*t5*t2**3+100._ki/3._ki*t6**3*t1&
                  &**2*t4**4*t3**2*t2-605._ki/18._ki*t6**3*t1**3*t4*t3*t5*t2**3-95._k&
                  &i/3._ki*t6**3*t1**2*t4**3*t3**2*t2**2-8._ki*t6**3*t1**2*t4**5*t3*&
                  &t2**2-35._ki/12._ki*t6**3*t1**2*t4*t3**2*t2**4+7._ki/8._ki*t6**2*t1&
                  &**3*t2**6*t4*t5-15._ki/4._ki*t6**3*t1**2*t4**2*t5*t2**4*t3-5._ki/4&
                  &._ki*t6**2*t1**3*t4*t5**3*t2**4-15._ki/4._ki*t6**2*t1**3*t3*t5**2*&
                  &t2**4-135._ki/2._ki*t6**2*t1**3*t4**2*t5**2*t3*t2**2+45._ki*t6**2*&
                  &t1**3*t4**3*t5**2*t3*t2-t1**3*t4**2*t5*t2**7/9._ki-40._ki/3._ki*t6&
                  &**3*t1**2*t4**5*t3**2-4._ki/3._ki*t6**3*t1*t4**7*t3**2+5._ki/3._ki*&
                  &t6**3*t1**5*t3*t5+t6**3*t1**5*t2**2*t4
                !
                stemp9=10._ki/3._ki*t6**3*t1**5*t5**2*t4-5._ki/3._ki*t6**3*t1**5*t5**&
                  &2*t2+4._ki/3._ki*t6**3*t1**5*t5*t2**2-2._ki/3._ki*t6**3*t1**5*t2*t3&
                  &-15._ki/2._ki*t6**3*t1**4*t4**3*t2**2+45._ki/4._ki*t6**3*t1**4*t4**&
                  &2*t2**3-5._ki*t6**3*t1**4*t4*t2**4-25._ki/2._ki*t6**3*t1**4*t3**2*&
                  &t4-25._ki*t6**3*t1**4*t4**3*t5**2+5._ki*t6**3*t1**4*t2**3*t3+33._k&
                  &i/4._ki*t6**3*t1**3*t4**3*t2**4+11._ki/3._ki*t6**3*t1**3*t4**5*t2*&
                  &*2-55._ki/9._ki*t6**3*t1**3*t3**2*t2**3-t1**3*t4*t3*t2**7/18._ki+t&
                  &1**3*t4*t2**8*t5/9._ki+110._ki/9._ki*t6**3*t1**3*t4**5*t5**2-77._ki&
                  &/24._ki*t6**3*t1**3*t4**2*t2**5-55._ki/6._ki*t6**3*t1**3*t4**4*t2*&
                  &*3-15._ki/4._ki*t6**2*t1**3*t3**3*t2**2-5._ki/2._ki*t6**2*t1**4*t4*&
                  &*2*t2**4-7._ki/6._ki*t6**2*t1**4*t2**4*t3
                !
                stemp8=stemp9+5._ki/24._ki*t6**3*t1**2*t3**2*t2**5-t6**2*t2**7*t3**&
                  &3*t4/36._ki-5._ki/3._ki*t6**3*t1**4*t5*t2**4+10._ki/9._ki*t6**2*t1**&
                  &2*t3**3*t2**4+10._ki/3._ki*t6**2*t4**4*t3**3*t2**4-40._ki/9._ki*t6*&
                  &*2*t4**5*t3**3*t2**3+28._ki/9._ki*t6**2*t4**6*t3**3*t2**2-8._ki/9.&
                  &_ki*t6**2*t4**7*t3**3*t2+5._ki/2._ki*t6**2*t1**4*t4*t2**5+15._ki/4.&
                  &_ki*t6**2*t1**3*t4**2*t2**6+3._ki*t6**2*t1**3*t4**4*t2**4-6._ki*t6&
                  &**2*t1**3*t4**3*t2**5+t6*t1**3*t4**3*t2**6/2._ki-3._ki/4._ki*t6*t1&
                  &**3*t4**2*t2**7+25._ki/12._ki*t6**3*t1**4*t5**2*t2**3+7._ki/8._ki*t&
                  &6**2*t1**3*t3*t2**6-5._ki/6._ki*t6**2*t1**4*t5**3*t2**3+7._ki/12._k&
                  &i*t6**2*t1**4*t5*t2**5-3._ki/4._ki*t6**2*t1**3*t4*t2**7+t6*t1**3*&
                  &t4*t2**8/4._ki+11._ki/36._ki*t6**2*t2**6*t3**3*t4**2
                !
                stemp9=stemp8+5._ki/72._ki*t6**2*t1*t2**6*t3**3-25._ki/18._ki*t6**2*t&
                  &2**5*t3**3*t4**3-11._ki/9._ki*t6**3*t1**3*t2**5*t3+11._ki/24._ki*t6&
                  &**3*t1**3*t4*t2**6+1045._ki/18._ki*t6**3*t1**3*t3**2*t4**3+50._ki/&
                  &3._ki*t6**2*t1**2*t4**4*t3*t5**2*t2**2-40._ki/9._ki*t6**2*t1*t4**5&
                  &*t3**3*t2+100._ki/9._ki*t6**2*t1*t4**4*t3**3*t2**2-11._ki/9._ki*t6*&
                  &*3*t1**3*t4*t5*t2**5+55._ki/36._ki*t6**3*t1**3*t4*t5**2*t2**4+85.&
                  &_ki/18._ki*t6**2*t1*t2**4*t3**3*t4**2-35._ki/12._ki*t6**2*t1**4*t4*&
                  &t5*t2**4-20._ki/3._ki*t6**3*t1*t4**5*t3**2*t2**2-8._ki/3._ki*t6**3*&
                  &t1**5*t5*t2*t4+14._ki/3._ki*t6**3*t1*t4**6*t3**2*t2-t6**2*t1**4*t&
                  &2**6/2._ki-25._ki/6._ki*t6**2*t1**4*t4**2*t5**3*t2+35._ki/12._ki*t6*&
                  &*2*t1**4*t4**2*t5*t2**3+stemp7+35._ki/6._ki*t6**2*t1**2*t4**2*t3*&
                  &t5**2*t2**4
                !
                stemp6=stemp9+220._ki/9._ki*t6**3*t1**3*t4**4*t5*t2**2+165._ki/4._ki*&
                  &t6**3*t1**3*t4*t3**2*t2**2+121._ki/9._ki*t6**3*t1**3*t4*t3*t2**4+&
                  &770._ki/9._ki*t6**3*t1**3*t4**4*t3*t5+77._ki/9._ki*t6**3*t1**3*t4**&
                  &2*t5*t2**4-1540._ki/9._ki*t6**3*t1**3*t4**3*t3*t5*t2-1045._ki/12._k&
                  &i*t6**3*t1**3*t3**2*t4**2*t2-385._ki/36._ki*t6**3*t1**3*t4**2*t5*&
                  &*2*t2**3-143._ki/3._ki*t6**3*t1**3*t4**2*t3*t2**3+616._ki/9._ki*t6*&
                  &*3*t1**3*t4**3*t3*t2**2+5._ki*t6**3*t1*t4**4*t3**2*t2**3-25._ki/2&
                  &._ki*t6**3*t1**4*t3*t5*t2**2+75._ki/2._ki*t6**3*t1**4*t4**2*t5**2*&
                  &t2-65._ki/3._ki*t6**3*t1**4*t2**2*t3*t4-30._ki*t6**3*t1**4*t4**2*t&
                  &5*t2**2-10._ki*t6**2*t1**3*t4**3*t5**3*t2**2+5._ki*t6**2*t1**3*t4&
                  &**4*t5**3*t2-35._ki/36._ki*t6**2*t1*t2**5*t3**3*t4+25._ki/4._ki*t6*&
                  &*2*t1**3*t4**2*t5**3*t2**3+30._ki*t6**2*t1**3*t4*t5**2*t3*t2**3-&
                  &95._ki/9._ki*t6**2*t1**2*t3**3*t4**3*t2+95._ki/6._ki*t6**2*t1**2*t3&
                  &**3*t4**2*t2**2
                !
                stemp7=1._ki/t1**3/t2**10
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
              !
            case(4)
              !
              select case(par4_glob)
              !
              case(4)
                !
                stemp2=-(-24._ki*t6**3*t2**2*t5*t4**2+48._ki*t6**3*t4**3*t5*t2-24._k&
                  &i*t6**3*t3*t2*t1+60._ki*t6**3*t4*t1*t3+96._ki*t6**3*t3*t2*t4**2-3&
                  &6._ki*t6**3*t3*t2**2*t4+2._ki*t6*t1*t2**4*t5+2._ki*t6*t2**5*t4*t5-&
                  &4._ki*t6*t2**4*t4**2*t5-6._ki*t3*t4*t6*t2**4-18._ki*t6**2*t1*t2**3&
                  &*t4+2._ki*t3*t5*t2**5+6._ki*t6**2*t1*t2**4-12._ki*t6**3*t1**2*t5-3&
                  &2._ki*t6**3*t4**4*t5-12._ki*t6**2*t2**4*t4**2+96._ki*t6**3*t1*t5*t&
                  &4**2+12._ki*t1*t6**3*t2**2*t5+4._ki*t6**3*t2**3*t5*t4-72._ki*t1*t6&
                  &**3*t4*t5*t2+t4*t5**2*t2**5+12._ki*t6**2*t2**3*t4**3+3._ki*t6**2*&
                  &t2**5*t4-80._ki*t6**3*t3*t4**3+4._ki*t6**3*t3*t2**3+2._ki*t6*t2**5&
                  &*t3)/t2**8*z_log(t1*t6/t2**2,1._ki)/12._ki
                !
                stemp4=(-6._ki*t2*t6*t3**3*t4+2._ki*t2**2*t5*t3**3-4._ki*t3**2*t4**2&
                  &*t5*t2*t6+2._ki*t4*t5*t3**2*t6*t2**2+6._ki*t1*t6**2*t3**2*t4-8._ki&
                  &*t3*t4*t5*t2*t6**2*t1+4._ki*t5*t1**2*t6**2*t3-6._ki*t2*t5*t1*t6*t&
                  &3**2-4._ki*t3*t4*t5**2*t2*t1*t6+t3**2*t4*t5**2*t2**2+16._ki*t3*t4&
                  &**2*t5*t6**2*t1+2._ki*t6*t3**3*t2**2-12._ki*t2*t4**2*t3**2*t6**2+&
                  &12._ki*t4**3*t3**2*t6**2-2._ki*t2*t1*t6**2*t3**2+3._ki*t2**2*t4*t3&
                  &**2*t6**2+6._ki*t4*t5**2*t6**2*t1**2)/t3**2/t2**5*q(4,(t2*t3-t1*&
                  &t6)/t2/t3,sign_arg)/12._ki
                !
                stemp5=-(-8._ki*t6*t1**2*t2**4*t5*t3-56._ki*t6*t1*t2**4*t3*t4**2*t5&
                  &+28._ki*t6*t1*t2**5*t3*t4*t5-48._ki*t6*t1*t4*t3**2*t2**4-360._ki*t&
                  &6**3*t1**2*t3*t4*t5*t2+44._ki*t6**3*t1*t3*t2**3*t5*t4+96._ki*t6**&
                  &3*t3**2*t4**5-352._ki*t6**3*t1*t3*t4**4*t5-396._ki*t6**3*t1*t2**2&
                  &*t4*t3**2+480._ki*t6**3*t1**2*t3*t5*t4**2+60._ki*t6**3*t1**2*t3*t&
                  &2**2*t5+6._ki*t6**3*t4*t3**2*t2**4+44._ki*t6**3*t1*t3**2*t2**3+14&
                  &4._ki*t6**3*t3**2*t2**2*t4**3+18._ki*t6**2*t1**2*t3*t2**4+16._ki*t&
                  &1*t2**5*t3**2*t5-48._ki*t6**3*t3**2*t2**3*t4**2-192._ki*t6**3*t3*&
                  &*2*t4**4*t2+14._ki*t1*t2**5*t3*t4*t5**2+528._ki*t6**3*t1*t3*t4**3&
                  &*t5*t2+108._ki*t6**2*t1*t3*t2**3*t4**3-108._ki*t6**2*t1*t3*t2**4*&
                  &t4**2-54._ki*t6**2*t1**2*t3*t2**3*t4+27._ki*t6**2*t1*t3*t2**5*t4-&
                  &12._ki*t6*t1**2*t2**4*t4*t5**2+1056._ki*t6**3*t1*t2*t4**2*t3**2-2&
                  &64._ki*t6**3*t1*t3*t2**2*t5*t4**2+16._ki*t6*t1*t2**5*t3**2-24._ki*&
                  &t6**3*t3*t5*t1**3+300._ki*t6**3*t1**2*t4*t3**2-120._ki*t6**3*t1**&
                  &2*t2*t3**2-880._ki*t6**3*t1*t4**3*t3**2)/t1/t2**8/t3/72._ki
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
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
                stemp6=-3._ki/4._ki*t4**2*t2**10+t4**3*t2**9/4._ki-132._ki*t6**3*t4**&
                  &2*t3*t2**5+36._ki*t6**3*t4**5*t3*t2**2+180._ki*t6**3*t4**5*t5**3*&
                  &t2-330._ki*t6**3*t4**4*t5**3*t2**2+66._ki*t6**3*t4**4*t5*t2**4-13&
                  &2._ki*t6**3*t4**4*t3*t2**3+189._ki*t6**3*t4**3*t3*t2**4-540._ki*t6&
                  &**3*t4**5*t3*t5**2-36._ki*t6**3*t4**5*t5*t2**3-135._ki*t6**3*t1**&
                  &2*t4**2*t5**3+240._ki*t6**3*t1*t4**4*t5**3-3._ki/4._ki*t6**2*t1**2&
                  &*t4*t2**5-9._ki*t6**2*t1*t4**2*t2**6+7._ki/2._ki*t6**2*t1*t4**3*t2&
                  &**5+105._ki*t6**2*t4*t3**2*t2**5-285._ki/2._ki*t6**2*t4**2*t3**2*t&
                  &2**4+57._ki/2._ki*t6**2*t4**2*t3*t2**6-25._ki*t6**2*t4**3*t3*t2**5&
                  &-250._ki*t6**2*t4**3*t5**4*t2**3+160._ki*t6**2*t4**4*t5**4*t2**2-&
                  &40._ki*t6**2*t4**5*t5**4*t2+125._ki/2._ki*t6**2*t4**3*t3**2*t2**3+&
                  &8._ki*t6**2*t4**4*t3*t2**4-250._ki*t6**2*t4**3*t5**3*t2**4+160._ki&
                  &*t6**2*t4**4*t5**3*t2**3-40._ki*t6**2*t4**5*t5**3*t2**2+25._ki/2.&
                  &_ki*t6**2*t4**3*t5*t2**6-8._ki*t6**2*t4**4*t5*t2**5+2._ki*t6**2*t4&
                  &**5*t5*t2**4+10._ki*t6**2*t4**5*t5**2*t2**3-40._ki*t6**2*t4**4*t5&
                  &**2*t2**4+125._ki/2._ki*t6**2*t4**3*t5**2*t2**5+70._ki*t6*t4**3*t5&
                  &**4*t2**4-20._ki*t6*t4**4*t5**4*t2**3
                !
                stemp5=7._ki*t6*t4**3*t5*t2**7-2._ki*t6*t4**4*t5*t2**6+35._ki*t6*t4*&
                  &*3*t5**5*t2**3-10._ki*t6*t4**4*t5**5*t2**2-8._ki*t6*t4**4*t5**3*t&
                  &2**4+18._ki*t6*t4**2*t3*t2**7+28._ki*t6*t4**3*t5**3*t2**5-7._ki*t6&
                  &*t4**3*t3*t2**6-15._ki*t6*t4*t3**2*t2**6+9._ki*t6*t4**2*t3**2*t2*&
                  &*5-7._ki*t6*t4**3*t5**2*t2**6+2._ki*t6*t4**4*t5**2*t2**5-9._ki*t6*&
                  &*3*t4*t5*t2**7+45._ki*t6**3*t4*t5**3*t2**5+90._ki*t6**3*t3*t5**2*&
                  &t2**5+75._ki*t6**3*t1*t5**3*t2**4-15._ki*t6**3*t1*t5*t2**6-70._ki*&
                  &t6**2*t4*t5**3*t2**6+35._ki/2._ki*t6**2*t4*t5**2*t2**7+25._ki*t6**&
                  &2*t5*t2**7*t3-70._ki*t6**2*t4*t5**4*t2**5-150._ki*t6**2*t3*t5**2*&
                  &t2**6-200._ki*t6**2*t5**3*t2**5*t3+stemp6+45._ki*t6**3*t4*t3*t2**&
                  &6+600._ki*t6**2*t3**2*t2**4*t5**2+190._ki*t6**2*t4**2*t5**4*t2**4&
                  &-3._ki/2._ki*t6**2*t1**2*t5*t2**5-15._ki/2._ki*t6**2*t1**2*t5**2*t2&
                  &**4+30._ki*t6**2*t1**2*t5**3*t2**3+30._ki*t6**2*t1**2*t5**4*t2**2&
                  &-165._ki*t6**3*t4**2*t5**3*t2**4+33._ki*t6**3*t4**2*t5*t2**6-t6*t&
                  &1*t2**5*t3**2+280._ki*t3**3*t2**3*t6**2*t4+160._ki*t3**3*t2**4*t6&
                  &*t5-190._ki*t3**3*t2**2*t6**2*t4**2-12._ki*t3*t2**3*t6**3*t1**2
                !
                stemp6=stemp5+20._ki*t3**2*t2**4*t1*t6**2-20._ki*t3**3*t2**4*t6*t4-&
                  &25._ki*t6*t5*t2**2*t3**4+660._ki*t2*t3**3*t6**3*t4**2-450._ki*t2**&
                  &2*t3**3*t6**3*t4+105._ki*t3**3*t6**3*t1*t4-6._ki*t6**2*t2**6*t1*t&
                  &3-3._ki*t6*t1*t3*t2**7-t2**2*t1**3*t6**3*t5-45._ki*t6**3*t5*t1**2&
                  &*t3**2+200._ki*t2**3*t3**3*t6*t5**2-400._ki*t2**3*t3**3*t6**2*t5-&
                  &9._ki*t6*t4**2*t2**8*t5-360._ki*t6*t3**2*t2**5*t5**2+9._ki*t6*t4**&
                  &2*t5**2*t2**7-300._ki*t6*t3**2*t5**3*t2**4-15._ki*t6*t4*t3*t2**8-&
                  &90._ki*t6*t4**2*t5**4*t2**5+190._ki*t6**2*t4**2*t5**3*t2**5-19._ki&
                  &/2._ki*t6**2*t4**2*t5*t2**7-14._ki*t6**2*t4*t3*t2**7+300._ki*t6**2&
                  &*t2**5*t3**2*t5+15._ki/2._ki*t6**2*t1*t4*t2**7-225._ki*t6**3*t3**2&
                  &*t2**4*t5-90._ki*t6**3*t1**2*t5**3*t2**2+18._ki*t6**3*t1**2*t5*t2&
                  &**4+30._ki*t6**3*t1*t3*t2**5+t6**2*t1**2*t3*t2**4/2._ki-75._ki*t3*&
                  &*3*t2*t6**3*t1-70._ki*t3**4*t2*t6**2*t4+20._ki*t6**2*t2**2*t1*t3*&
                  &*3-72._ki*t6*t5*t3**2*t2**6-36._ki*t6*t4**2*t5**3*t2**6-95._ki/2._k&
                  &i*t6**2*t4**2*t5**2*t2**6+1120._ki*t2**2*t3**3*t6**2*t5*t4-760._k&
                  &i*t2*t3**3*t6**2*t5*t4**2
                !
                stemp7=stemp6-375._ki*t6*t4*t5**4*t2**4*t3+30._ki*t6*t4*t3*t2**7*t5&
                  &-600._ki*t6*t4*t5**3*t2**5*t3-250._ki*t2**2*t3**3*t6*t4*t5**2-420&
                  &._ki*t6**3*t1*t4*t5**3*t2**3-450._ki*t6**3*t1*t3*t5**2*t2**3-180.&
                  &_ki*t6*t4*t5**2*t2**6*t3-75._ki*t6*t1*t5**4*t3*t2**3+1350._ki*t6**&
                  &3*t3**2*t2**3*t4*t5+1890._ki*t6**3*t1*t4*t3*t5**2*t2**2-2565._ki*&
                  &t6**3*t1*t4**2*t3*t5**2*t2+1125._ki*t6**3*t1*t4**3*t3*t5**2-126.&
                  &_ki*t6**3*t1*t2**4*t3*t4+855._ki*t6**3*t1*t4**2*t5**3*t2**2+171._k&
                  &i*t6**3*t1*t2**3*t3*t4**2-750._ki*t6**3*t1*t4**3*t5**3*t2+150._ki&
                  &*t6**3*t1*t4**3*t5*t2**3-171._ki*t6**3*t1*t4**2*t5*t2**4
                !
                stemp4=stemp7+36._ki*t6*t1*t4*t3*t2**4*t5**2+750._ki*t6*t4*t5**3*t3&
                  &**2*t2**3+720._ki*t6*t4**2*t3*t5**3*t2**4+450._ki*t6*t4**2*t5**4*&
                  &t3*t2**3+180._ki*t6*t4*t5*t2**5*t3**2-108._ki*t6*t4**2*t3**2*t5*t&
                  &2**4+216._ki*t6*t4**2*t5**2*t2**5*t3-84._ki*t6*t4**3*t5**2*t2**4*&
                  &t3-540._ki*t6*t4**2*t3**2*t5**2*t2**3-175._ki*t6*t4**3*t3*t5**4*t&
                  &2**2-280._ki*t6*t4**3*t3*t5**3*t2**3+120._ki*t6*t1*t4*t5**3*t3*t2&
                  &**3-45._ki*t6*t4**2*t5**5*t2**4+14._ki*t6*t4**3*t5*t3*t2**5-36._ki&
                  &*t6*t4**2*t5*t3*t2**6+900._ki*t6*t4*t3**2*t5**2*t2**4-450._ki*t6*&
                  &t4**2*t5**3*t3**2*t2**2+12._ki*t6*t1*t4**2*t5**3*t2**4-3._ki*t6*t&
                  &1*t4**2*t5**2*t2**5
                !
                stemp7=15._ki*t6*t1*t4**2*t5**5*t2**2+30._ki*t6*t1*t4**2*t5**4*t2**&
                  &3+3._ki*t6*t1*t4**2*t2**6*t5+300._ki*t6**2*t1*t4*t5**3*t2**4+1120&
                  &._ki*t6**2*t4*t5**3*t2**4*t3+1980._ki*t6**3*t4**2*t3*t5**2*t2**3-&
                  &2835._ki*t6**3*t4**3*t3*t5**2*t2**2+1980._ki*t6**3*t4**4*t3*t5**2&
                  &*t2+3._ki/4._ki*t4*t2**11+2835._ki*t6**3*t4**3*t3**2*t5*t2-2970._ki&
                  &*t6**3*t4**2*t3**2*t5*t2**2-30._ki*t6**2*t1**2*t4*t5**4*t2+15._ki&
                  &/2._ki*t6**2*t1**2*t4*t5**2*t2**3-30._ki*t6**2*t1**2*t4*t5**3*t2*&
                  &*2-35._ki*t6**2*t1*t4**3*t5**2*t2**3-360._ki*t6**2*t1*t4**2*t5**4&
                  &*t2**2-7._ki*t6**2*t1*t4**3*t5*t2**4+140._ki*t6**2*t1*t4**3*t5**4&
                  &*t2
                !
                stemp6=stemp7+18._ki*t6**2*t1*t4**2*t5*t2**5+3._ki*t6*t1*t5*t2**8+1&
                  &60._ki*t6*t3*t5**3*t2**6+100._ki*t6*t5**4*t2**5*t3-8._ki*t6*t5*t3*&
                  &t2**8+7._ki/2._ki*t6**2*t4*t2**8*t5+20._ki*t6**2*t1*t5**2*t2**6+4.&
                  &_ki*t6**2*t1*t5*t2**7-80._ki*t6**2*t1*t5**4*t2**4-80._ki*t6**2*t1*&
                  &t5**3*t2**5+50._ki*t6*t4*t5**4*t2**6+5._ki*t6*t4*t2**9*t5+20._ki*t&
                  &6*t4*t5**3*t2**7-5._ki*t6*t2**8*t4*t5**2+25._ki*t6*t4*t5**5*t2**5&
                  &+48._ki*t6*t3*t5**2*t2**7+30._ki*t6*t1*t5**4*t2**5+15._ki*t6*t1*t5&
                  &**5*t2**4+12._ki*t6*t1*t5**3*t2**6
                !
                stemp5=stemp6-3._ki*t6*t1*t5**2*t2**7-990._ki*t6**3*t4**4*t3**2*t5-&
                  &63._ki*t6**3*t4**3*t5*t2**5+8._ki*t6**3*t4**6*t5*t2**2+315._ki*t6*&
                  &*3*t4**3*t5**3*t2**3-10._ki*t6*t2**3*t3**4+90._ki*t6**2*t1*t4**2*&
                  &t5**2*t2**4-360._ki*t6**2*t1*t4**2*t5**3*t2**3+140._ki*t6**2*t1*t&
                  &4**3*t5**3*t2**2+3._ki/2._ki*t6**2*t1**2*t4*t5*t2**4-90._ki*t6**2*&
                  &t1*t4**2*t5*t2**3*t3-900._ki*t6**2*t1*t4*t3*t5**2*t2**3+150._ki*t&
                  &6**2*t1*t4*t5*t2**4*t3+540._ki*t6**2*t1*t4**2*t3*t5**2*t2**2+720&
                  &._ki*t6**2*t1*t4**2*t5**3*t2*t3-t2**12/4._ki-6._ki*t6**3*t3*t2**7+&
                  &4._ki*t6**2*t4**4*t2**6-10._ki*t6*t5**4*t2**7+50._ki*t3**4*t2**2*t&
                  &6**2-25._ki*t6**2*t3**2*t2**6+10._ki*t6**2*t5**4*t2**6+4._ki*t6*t3&
                  &*t2**9+3._ki/4._ki*t6**2*t1**2*t2**6-t6**2*t4**5*t2**5+16._ki*t3**&
                  &3*t6*t2**5-5._ki*t6*t5**5*t2**6-t6*t2**10*t5+5._ki/2._ki*t6**2*t3*&
                  &t2**8-25._ki/4._ki*t6**2*t4**3*t2**7+19._ki/4._ki*t6**2*t4**2*t2**8&
                  &-7._ki/4._ki*t6**2*t4*t2**9-100._ki*t3**3*t2**4*t6**2-4._ki*t6*t5**&
                  &3*t2**8+100._ki*t2**3*t3**3*t6**3-2._ki*t6**2*t1*t2**8-315._ki*t3*&
                  &*3*t6**3*t4**3
                !
                stemp7=stemp5+10._ki*t6**2*t5**3*t2**7-5._ki/2._ki*t6**2*t5**2*t2**8&
                  &-40._ki*t6**3*t4**6*t5**3+t6*t5**2*t2**9+t6**3*t5*t2**8-5._ki*t6*&
                  &*3*t5**3*t2**6+5._ki*t1**3*t6**3*t5**3-t6**2*t2**9*t5/2._ki+6._ki*&
                  &t6*t3**2*t2**7+360._ki*t6**2*t1*t3*t5**2*t2**4-9._ki*t6**2*t1*t3*&
                  &t2**4*t4**2+15._ki*t6**2*t1*t3*t2**5*t4-30._ki*t3*t5**2*t2**2*t1*&
                  &*2*t6**2+60._ki*t3**2*t5**2*t2**3*t1*t6-240._ki*t3**2*t2**3*t6**2&
                  &*t5*t1+300._ki*t3**2*t2**2*t6**2*t4*t5*t1-200._ki*t3**3*t2**3*t6*&
                  &t4*t5
                !
                stemp6=stemp7+15._ki*t3*t2**2*t6**3*t4*t1**2-25._ki*t3**2*t2**3*t1*&
                  &t6**2*t4+5._ki*t3*t2**3*t1**2*t6**2*t5-6._ki*t6*t1*t2**5*t3*t4*t5&
                  &-1200._ki*t6**2*t1*t4*t5**3*t3*t2**2-2520._ki*t6**2*t4*t3**2*t5**&
                  &2*t2**3-1500._ki*t6**2*t4**3*t3**2*t5**2*t2-1260._ki*t6**2*t4*t3*&
                  &*2*t5*t2**4+1710._ki*t6**2*t4**2*t3**2*t5*t2**3+2000._ki*t6**2*t4&
                  &**3*t3*t5**3*t2**2+3420._ki*t6**2*t4**2*t3**2*t5**2*t2**2-750._ki&
                  &*t6**2*t4**3*t3**2*t5*t2**2-640._ki*t6**2*t4**4*t3*t5**3*t2+285.&
                  &_ki*t6**2*t4**2*t5*t2**5*t3-2280._ki*t6**2*t4**2*t5**3*t3*t2**3+1&
                  &500._ki*t6**2*t4**3*t3*t5**2*t2**3-480._ki*t6**2*t4**4*t3*t5**2*t&
                  &2**2+80._ki*t6**2*t4**4*t5*t2**3*t3-1710._ki*t6**2*t4**2*t3*t5**2&
                  &*t2**4
                !
                stemp7=stemp6+75._ki*t6*t1*t4*t5**4*t2**2*t3-15._ki*t6**2*t1*t2**6*&
                  &t4*t5+3._ki*t6*t1*t2**6*t4*t3+12._ki*t3**2*t2**4*t1*t6*t5+80._ki*t&
                  &6**2*t5*t2*t1*t3**3+180._ki*t2*t3*t1**2*t6**3*t5**2-480._ki*t2**2&
                  &*t3**2*t1*t6**2*t5**2+855._ki*t3**2*t1*t6**3*t5*t4**2+50._ki*t2**&
                  &2*t3**2*t1*t6*t5**3-40._ki*t2*t3*t1**2*t6**2*t5**3-225._ki*t3*t1*&
                  &*2*t6**3*t4*t5**2+450._ki*t2**2*t3**2*t1*t6**3*t5-1260._ki*t2*t3*&
                  &*2*t1*t6**3*t5*t4+600._ki*t2*t3**2*t1*t6**2*t4*t5**2-120._ki*t6*t&
                  &1*t5**3*t3*t2**4+840._ki*t6**2*t4*t3*t5**2*t2**5-75._ki*t6**2*t1*&
                  &t4*t5**2*t2**5+300._ki*t6**2*t1*t4*t5**4*t2**3
                !
                stemp3=stemp7-45._ki*t6**3*t1**2*t4*t5*t2**3+27._ki*t6**3*t1**2*t4*&
                  &*2*t5*t2**2+225._ki*t6**3*t1**2*t4*t5**3*t2-75._ki*t6**3*t1*t2**2&
                  &*t3*t4**3-48._ki*t6**3*t1*t4**4*t5*t2**2-140._ki*t6**2*t4*t5*t2**&
                  &6*t3+t6**2*t2**10/4._ki+stemp4-250._ki*t6**2*t4**3*t5*t2**4*t3-6.&
                  &_ki*t6*t1*t4*t5*t2**7-36._ki*t6*t1*t3*t5**2*t2**5+6._ki*t6*t1*t4*t&
                  &5**2*t2**6-675._ki*t6**3*t4*t3*t5**2*t2**4+84._ki*t6**3*t1*t4*t5*&
                  &t2**5+480._ki*t6**2*t1*t3*t5**3*t2**3-60._ki*t6**2*t1*t3*t5*t2**5&
                  &+6._ki*t6*t1*t2**6*t3*t5-24._ki*t6*t1*t4*t5**3*t2**5-60._ki*t6*t1*&
                  &t4*t5**4*t2**4-30._ki*t6*t1*t4*t5**5*t2**3
                !
                stemp4=1._ki/t2**12*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=(-t4+t2)**3/t2**3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/4._ki
                !
                stemp11=-300._ki*t6*t1**4*t2**5*t3**3*t4+110._ki*t6*t1**4*t4**4*t3*&
                  &*2*t2**4+660._ki*t6*t1**4*t2**4*t3**3*t4**2+300._ki*t6*t1**4*t3**&
                  &4*t2**3*t4-600._ki*t6*t1**6*t4*t5**3*t3*t2**2-450._ki*t6*t1**6*t4&
                  &*t3*t5**2*t2**3+300._ki*t6*t1**6*t4*t3**2*t5**2*t2-1890._ki*t6*t1&
                  &**5*t4*t3**2*t5*t2**4-3780._ki*t6*t1**5*t4*t3**2*t5**2*t2**3+256&
                  &5._ki*t6*t1**5*t4**2*t3**2*t5*t2**3+1680._ki*t6*t1**5*t3**3*t5*t2&
                  &**2*t4-375._ki*t6*t1**5*t4**3*t5*t2**4*t3-3420._ki*t6*t1**5*t4**2&
                  &*t5**3*t3*t2**3-32._ki*t6**2*t1**2*t4**8*t3**3*t2+121._ki*t6**2*t&
                  &1**2*t4**5*t3**3*t2**4-414._ki*t6**2*t1**3*t4**5*t3**2*t5*t2**3+&
                  &5._ki*t6**2*t1**7*t2**2*t3*t4+9._ki*t6**2*t1**7*t4**2*t5*t2**2+12&
                  &5._ki*t6**2*t1**6*t4**3*t5*t2**3+1425._ki/2._ki*t6**2*t1**6*t4**2*&
                  &t3**2*t5
                !
                stemp10=stemp11+285._ki/2._ki*t6**2*t1**6*t2**3*t3*t4**2+1875._ki/2.&
                  &_ki*t6**2*t1**6*t4**3*t3*t5**2-285._ki/2._ki*t6**2*t1**6*t4**2*t5*&
                  &t2**4-55._ki*t6*t4**2*t2**7*t1**4*t3*t5-840._ki*t6*t1**4*t4**3*t5&
                  &**3*t3*t2**4-1050._ki*t6**2*t1**6*t3**2*t5*t2*t4-600._ki*t1**5*t4&
                  &*t5**3*t2**5*t3-180._ki*t1**5*t4*t5**2*t2**6*t3-840._ki*t1**4*t4*&
                  &t3**2*t5**2*t2**6+70._ki*t6*t1**6*t4**3*t5**3*t2**2-285._ki*t6*t1&
                  &**5*t3**3*t2**2*t4**2-105._ki*t6*t1**5*t3**4*t2*t4+315._ki/2._ki*t&
                  &6*t1**5*t4*t3**2*t2**5-855._ki/4._ki*t6*t1**5*t4**2*t3**2*t2**4+1&
                  &50._ki*t6*t1**6*t4*t5*t3**2*t2**2+2160._ki*t6**2*t1**4*t4**4*t3**&
                  &2*t5*t2**2-4._ki*t6**2*t1**4*t4**7*t3*t2**2+240._ki*t6*t1**5*t4**&
                  &4*t5**4*t2**2+420._ki*t6*t1**5*t4*t3**3*t2**3-440._ki*t6*t1**4*t3&
                  &**4*t2**2*t4**2
                !
                stemp11=-30._ki*t6*t1**4*t4**5*t3**2*t2**3-315._ki/2._ki*t6*t1**4*t4&
                  &**3*t3**2*t2**5+110._ki*t6*t1**4*t4**2*t3**2*t2**6+165._ki*t6**2*&
                  &t3*t5**2*t2**5*t1**5-33._ki/2._ki*t6**2*t4*t5*t2**7*t1**5+900._ki*&
                  &t6*t1**5*t3**2*t2**4*t5**2+285._ki*t6*t1**5*t4**2*t5**4*t2**4+20&
                  &0._ki*t6*t5*t2**5*t3**3*t1**4-55._ki*t6*t4*t3**3*t2**7*t1**3-15._k&
                  &i/2._ki*t6*t1**6*t2**6*t4*t5+180._ki*t6*t1**6*t3*t5**2*t2**4-75._k&
                  &i/2._ki*t6*t1**6*t4*t5**2*t2**5+150._ki*t6*t1**6*t4*t5**4*t2**3+2&
                  &40._ki*t6*t1**6*t3*t5**3*t2**3-30._ki*t6*t1**6*t3*t5*t2**5+150._ki&
                  &*t6*t1**6*t4*t5**3*t2**4+900._ki*t6*t4*t3**2*t2**5*t1**4*t5**2+3&
                  &30._ki*t6*t4**2*t3*t2**6*t1**4*t5**2+stemp10-1720._ki*t6*t1**3*t4&
                  &**3*t3**3*t5*t2**4-64._ki*t1**2*t4**4*t3**3*t2**7
                !
                stemp9=stemp11+425._ki/2._ki*t1**2*t4**2*t3**4*t2**7-430._ki*t1**2*t&
                  &4**3*t3**4*t2**6+200._ki/3._ki*t1**2*t4**6*t3**4*t2**3+43._ki*t1**&
                  &2*t4**3*t3**3*t2**8+16._ki/3._ki*t1**2*t4**7*t3**3*t2**4-80._ki/3.&
                  &_ki*t1**2*t4**6*t3**3*t2**5+56._ki*t1**2*t4**5*t3**3*t2**6+480._ki&
                  &*t1**2*t4**4*t3**4*t2**5-280._ki*t1**2*t4**5*t3**4*t2**4+70._ki*t&
                  &1**2*t2**2*t3**5*t4**5-425._ki/2._ki*t1**2*t2**5*t3**5*t4**2-240.&
                  &_ki*t1**2*t2**3*t3**5*t4**4+275._ki/4._ki*t1**2*t4*t3**5*t2**6+645&
                  &._ki/2._ki*t1**2*t2**4*t3**5*t4**3+1155._ki/4._ki*t1*t2**6*t3**5*t4&
                  &**3+73._ki/2._ki*t4**2*t3**4*t2**9*t1+200._ki*t1**5*t3**3*t5**2*t2&
                  &**3+18._ki*t1**5*t4**2*t3*t2**7-15._ki*t1**5*t4*t3**2*t2**6+70._ki&
                  &*t1**5*t4**3*t5**4*t2**4-7._ki*t1**5*t4**3*t5**2*t2**6
                !
                stemp11=stemp9+15._ki*t6*t1**5*t4**5*t5**2*t2**3-1140._ki*t6*t1**5*&
                  &t4**2*t3**3*t5*t2+120._ki*t6*t1**5*t4**4*t5*t2**3*t3-1125._ki*t6*&
                  &t1**5*t4**3*t3**2*t5*t2**2+855._ki/2._ki*t6*t1**5*t4**2*t5*t2**5*&
                  &t3+3000._ki*t6*t1**5*t4**3*t3*t5**3*t2**2-12._ki*t6*t1**5*t4**4*t&
                  &5*t2**5-60._ki*t6*t1**5*t4**5*t5**3*t2**2+240._ki*t6*t1**5*t4**4*&
                  &t5**3*t2**3+3._ki*t6*t1**5*t4**5*t5*t2**4+12._ki*t6*t1**5*t4**4*t&
                  &3*t2**4+171._ki/4._ki*t6*t1**5*t4**2*t3*t2**6+110._ki/3._ki*t4*t5*t&
                  &3**3*t2**9*t1**2+720._ki*t1**5*t4**2*t3*t5**3*t2**4-690._ki*t1*t4&
                  &**5*t3**4*t5*t2**5+40._ki*t1*t4**8*t3**4*t5*t2**2-1155._ki/4._ki*t&
                  &1*t4**3*t5*t2**7*t3**4+1125._ki/2._ki*t1*t4**4*t5*t2**6*t3**4-4._k&
                  &i*t3**3*t2**9*t1**3
                !
                stemp10=stemp11-5._ki/2._ki*t6*t3**2*t5*t2**9*t1**3+73._ki/3._ki*t6*t&
                  &4**2*t3**3*t2**8*t1**2+160._ki*t1**5*t3**3*t5*t2**4-20._ki*t1**5*&
                  &t4**4*t5**4*t2**3-2._ki*t1**5*t4**4*t5*t2**6-25._ki*t1**5*t3**4*t&
                  &5*t2**2+35._ki*t1**5*t4**3*t5**5*t2**3-10._ki*t1**5*t4**4*t5**5*t&
                  &2**2+30._ki*t1**4*t3**2*t2**8*t5+125._ki*t1**4*t5**3*t2**6*t3**2+&
                  &7._ki/2._ki*t1**4*t3*t4*t2**10+150._ki*t1**4*t3**2*t2**7*t5**2-50.&
                  &_ki/3._ki*t1**3*t3**5*t2**5-18._ki*t4*t3**2*t2**9*t1**3*t5+32._ki/3&
                  &._ki*t6*t1**2*t4**8*t3**3*t2**2-39._ki/4._ki*t6**2*t3**2*t4*t2**7*&
                  &t1**3*t5+165._ki/4._ki*t6**2*t3*t4*t2**6*t1**4*t5**2+44._ki/3._ki*t&
                  &6**2*t1**5*t4**6*t5*t2**2-33._ki/4._ki*t6**2*t4**2*t3**3*t2**7*t1&
                  &**2-90._ki*t6*t5**2*t4*t2**7*t1**4*t3-120._ki*t6*t4*t5**3*t2**6*t&
                  &1**4*t3
                !
                stemp11=stemp10+15._ki*t6*t4*t5*t3*t2**8*t1**4-13._ki/3._ki*t6*t4*t3&
                  &**3*t2**9*t1**2-55._ki/24._ki*t6*t4*t2**9*t1**3*t3**2-120._ki*t6*t&
                  &3**2*t5**2*t2**6*t1**4-5._ki*t1**5*t5**5*t2**6+5._ki/2._ki*t6*t3**&
                  &4*t4*t2**9*t1+5._ki/4._ki*t6**2*t3**3*t4*t2**8*t1**2+3._ki/4._ki*t6&
                  &**2*t3**2*t5*t2**8*t1**3+275._ki*t4**2*t5**3*t2**6*t1**3*t3**2+5&
                  &00._ki/3._ki*t1**2*t4**6*t3**4*t5*t2**2+160._ki/3._ki*t1**2*t4**7*t&
                  &3**3*t5*t2**3-7._ki/2._ki*t6*t1**6*t4**3*t5*t2**4+t3**4*t2**11*t1&
                  &/2._ki+430._ki*t1**2*t4**3*t5*t3**3*t2**7+2125._ki/4._ki*t1**2*t4**&
                  &2*t3**4*t5*t2**6+1200._ki*t1**2*t4**4*t3**4*t5*t2**4+1075._ki/2._k&
                  &i*t1**2*t4**3*t5**2*t2**6*t3**3-1120._ki*t6*t1**3*t4**5*t3**3*t5&
                  &*t2**2+800._ki/3._ki*t6*t1**3*t4**6*t3**3*t5*t2-168._ki*t1**4*t3**&
                  &2*t5*t4*t2**7
                !
                stemp8=stemp11-t3**5*t2**11/6._ki+60._ki*t6**2*t1**7*t3*t5**2*t2-13&
                  &2._ki*t6**2*t1**3*t4**7*t3**2*t5*t2-15._ki/4._ki*t6**2*t5**2*t2**7&
                  &*t1**4*t3+375._ki/4._ki*t6*t1**5*t4**3*t5**2*t2**5+55._ki/2._ki*t6*&
                  &t4*t5*t3**2*t2**8*t1**3+55._ki*t6*t3**2*t4*t2**7*t1**3*t5**2-52.&
                  &_ki/3._ki*t6*t3**3*t4*t2**8*t1**2*t5-105._ki*t6*t1**5*t4*t5**4*t2*&
                  &*5-225._ki*t6*t1**5*t3*t5**2*t2**6-495._ki/2._ki*t6**2*t4*t5*t2**5&
                  &*t3**2*t1**4-300._ki*t6*t1**5*t5**3*t2**5*t3-60._ki*t6*t5*t3**2*t&
                  &2**7*t1**4-75._ki/2._ki*t6*t4*t2**7*t3**2*t1**4-210._ki*t6*t1**5*t&
                  &4*t5*t2**6*t3+285._ki*t6*t1**5*t4**2*t5**3*t2**5+320._ki/3._ki*t6*&
                  &t1**4*t4**6*t5**3*t3*t2+24._ki*t6**2*t1**3*t4**8*t3**2*t5+880._ki&
                  &*t6*t1**4*t4**4*t5**3*t3*t2**3-2640._ki*t6*t1**4*t4**2*t3**2*t5*&
                  &*2*t2**4-960._ki*t6*t1**5*t4**4*t3*t5**3*t2
                !
                stemp11=stemp8+80._ki*t6*t1**4*t4**6*t3*t5**2*t2**2-630._ki*t6*t1**&
                  &4*t4**3*t3*t5**2*t2**5-630._ki*t6*t1**4*t2**3*t3**3*t4**3+2640._k&
                  &i*t6*t1**4*t4**2*t3**3*t2**3*t5-2640._ki*t6*t1**4*t4**4*t3**2*t5&
                  &**2*t2**2+720._ki*t6*t1**4*t4**5*t3**2*t5**2*t2+60._ki*t6*t1**4*t&
                  &4**5*t3*t2**4*t5-110._ki*t6*t1**4*t4**4*t3*t2**5*t5+3780._ki*t6*t&
                  &1**4*t4**3*t3**2*t5**2*t2**3+6._ki*t6*t1**4*t4**5*t3*t2**5+19._ki&
                  &*t1**4*t4**2*t3*t2**8*t5+375._ki*t4*t3**3*t2**6*t1**3*t5**2-736.&
                  &_ki*t6*t1**2*t4**5*t3**3*t5*t2**4+600._ki*t6*t1**2*t4**4*t5*t2**5&
                  &*t3**3+242._ki*t6*t1*t4**5*t3**4*t2**5-275._ki/2._ki*t4*t5*t2**7*t&
                  &3**4*t1**2+1920._ki*t6*t1**3*t4**4*t3**3*t5*t2**3+385._ki*t6*t1**&
                  &2*t3**4*t4**3*t2**5+21._ki/2._ki*t6*t1**4*t4**3*t3*t2**7
                !
                stemp10=stemp11+220._ki*t6*t1**4*t4**4*t3**3*t2**2-4._ki/3._ki*t6*t1&
                  &**4*t4**6*t3*t2**4-180._ki*t6*t1**6*t4**2*t5**4*t2**2+1425._ki*t1&
                  &**4*t4**2*t3**2*t5**3*t2**4-240._ki*t6*t1**6*t3**2*t5**2*t2**2-1&
                  &20._ki*t6*t1**6*t2**3*t3**2*t5-25._ki/2._ki*t6*t1**6*t2**3*t4*t3**&
                  &2-9._ki/2._ki*t6*t1**6*t4**2*t3*t2**4+45._ki*t6*t1**6*t4**2*t5**2*&
                  &t2**4+15._ki/2._ki*t6*t1**6*t4*t3*t2**5+1890._ki*t6*t1**4*t4**3*t3&
                  &**2*t5*t2**4-1320._ki*t6*t1**4*t4**4*t3**2*t5*t2**3+360._ki*t6*t1&
                  &**4*t4**5*t3**2*t5*t2**2+3._ki/8._ki*t6*t1**5*t2**10-90._ki*t4*t3*&
                  &*2*t5**2*t2**8*t1**3+175._ki/2._ki*t1**4*t4*t5**4*t3*t2**6+140._ki&
                  &*t1**4*t4*t5**3*t3*t2**7+42._ki*t1**4*t3*t5**2*t4*t2**8-7._ki*t1*&
                  &*4*t2**9*t4*t5*t3+3630._ki*t6**2*t1**5*t4**2*t3*t5**2*t2**3-5445&
                  &._ki*t6**2*t1**5*t4**2*t3**2*t5*t2**2
                !
                stemp11=stemp10+3825._ki/4._ki*t6**2*t1**4*t4**2*t3**2*t5*t2**4+105&
                  &._ki*t6*t1**4*t4**3*t3*t2**6*t5+300._ki*t4*t3**3*t2**7*t1**3*t5+8&
                  &0._ki*t1**3*t4**6*t3**2*t5**2*t2**3-146._ki*t6*t1**2*t4**2*t3**4*&
                  &t2**6-600._ki*t6*t1**2*t2**4*t4**4*t3**4-77._ki*t6*t1**2*t4**3*t3&
                  &**3*t2**7+150._ki*t6*t1**2*t4**4*t3**3*t2**6-176._ki/3._ki*t6*t1**&
                  &2*t4**7*t3**3*t2**3-704._ki/3._ki*t6*t1**2*t4**7*t3**3*t5*t2**2+1&
                  &28._ki/3._ki*t6*t1**2*t4**8*t3**3*t5*t2-64._ki*t6*t1*t2**2*t3**4*t&
                  &4**8-308._ki*t6*t1**2*t4**3*t5*t2**6*t3**3-33._ki/2._ki*t6*t4**2*t&
                  &3**4*t2**8*t1+645._ki*t6*t1**3*t4**3*t5**2*t2**5*t3**2+10395._ki/&
                  &2._ki*t6**2*t1**5*t4**3*t3**2*t5*t2+425._ki/2._ki*t6*t1**3*t4**2*t&
                  &3**3*t2**6-425._ki*t6*t1**3*t4**2*t3**4*t2**4-35._ki*t6*t1**3*t4*&
                  &*5*t3**2*t2**5+40._ki*t6*t1**3*t4**4*t3**2*t2**6
                !
                stemp9=stemp11-215._ki/8._ki*t6*t1**3*t4**3*t3**2*t2**7+630._ki*t6**&
                  &2*t1**4*t4**5*t3*t5**2*t2**2-300._ki*t6**2*t1**4*t4**6*t3*t5**2*&
                  &t2+675._ki/2._ki*t6**2*t1**3*t4**4*t5*t2**4*t3**2+312._ki*t6**2*t1&
                  &**3*t4**6*t3**2*t5*t2**2-693._ki/4._ki*t6**2*t1**3*t4**3*t5*t2**5&
                  &*t3**2+414._ki*t6**2*t1**3*t4**5*t3**3*t2**2+84._ki*t6**2*t1**2*t&
                  &4**7*t3**3*t2**2-76._ki*t6**2*t1**2*t4**4*t3**3*t2**5-380._ki/3._k&
                  &i*t6**2*t1**2*t4**6*t3**3*t2**3-25._ki/2._ki*t1**4*t3**5*t2**3+4.&
                  &_ki*t1**5*t3*t2**9+40._ki*t6*t1**6*t3**3*t5*t2+9._ki*t6*t1**6*t4**&
                  &2*t5*t2**5-180._ki*t6*t1**6*t4**2*t5**3*t2**3+70._ki*t6*t1**6*t4*&
                  &*3*t5**4*t2-35._ki/2._ki*t6*t1**6*t4**3*t5**2*t2**3+330._ki*t4**2*&
                  &t3**2*t2**7*t1**3*t5**2+275._ki/6._ki*t3**3*t4*t2**8*t1**2*t5**2-&
                  &65._ki/4._ki*t3**4*t4*t2**9*t1*t5-375._ki*t1**5*t4*t5**4*t2**4*t3
                !
                stemp11=stemp9+30._ki*t1**5*t4*t3*t2**7*t5+625._ki*t1**4*t4**3*t3**&
                  &3*t5**2*t2**2-200._ki*t1**4*t4**4*t5**4*t3*t2**3-960._ki*t6*t1**3&
                  &*t4**4*t3**2*t5**2*t2**4+400._ki*t1**4*t4**4*t3**2*t5**3*t2**2+3&
                  &42._ki*t1**4*t4**2*t5*t2**6*t3**2+480._ki*t1**4*t4**4*t3**2*t5**2&
                  &*t2**3+1664._ki/3._ki*t6*t1**2*t4**6*t3**3*t5*t2**3+3._ki/2._ki*t1*&
                  &*3*t3**2*t4*t2**10-8._ki*t1**4*t4**4*t3**2*t2**5-57._ki/2._ki*t1**&
                  &4*t4**2*t3**2*t2**7+84._ki*t1**4*t4*t3**3*t2**6+7._ki*t1**5*t4**3&
                  &*t5*t2**7+414._ki*t1*t2**4*t3**5*t4**5-450._ki*t1*t2**5*t3**5*t4*&
                  &*4-88._ki*t1*t4**7*t3**4*t2**4+225._ki*t1*t4**4*t3**4*t2**7-276._k&
                  &i*t1*t4**5*t3**4*t2**6-231._ki/2._ki*t1*t4**3*t3**4*t2**8
                !
                stemp10=stemp11-208._ki*t1*t2**3*t3**5*t4**6-8._ki*t1**4*t4**4*t3*t&
                  &2**7+25._ki/2._ki*t1**4*t4**3*t3*t2**8+35._ki/2._ki*t1**4*t3**5*t2*&
                  &*2*t4+2._ki*t1**4*t4**5*t3*t2**6+25._ki*t1**4*t4**3*t3**2*t2**6-1&
                  &1._ki*t1**3*t4**4*t3**2*t2**7-150._ki*t1**3*t4*t3**4*t2**6+75._ki*&
                  &t1**3*t3**5*t2**4*t4+5._ki/4._ki*t3**4*t5*t2**10*t1-25._ki/6._ki*t3&
                  &**3*t5**2*t2**9*t1**2+t1**5*t5**2*t2**9-220._ki*t6*t4*t5*t2**6*t&
                  &3**3*t1**3+51._ki/4._ki*t6**2*t4**2*t3*t2**7*t1**4-10._ki*t1**5*t2&
                  &**3*t3**4-2._ki*t2**9*t3**5*t1-t1**5*t2**10*t5-20._ki*t1**4*t3**3&
                  &*t2**7+32._ki/3._ki*t2**2*t3**5*t4**9+6._ki*t1**5*t3**2*t2**7-152.&
                  &_ki*t4**4*t3**5*t2**7
                !
                stemp11=stemp10-10._ki*t1**5*t5**4*t2**7+25._ki*t3**4*t2**7*t1**3+1&
                  &68._ki*t2**4*t3**5*t4**7-35._ki/4._ki*t3**5*t2**7*t1**2-t1**3*t3**&
                  &2*t2**11/6._ki-5._ki/2._ki*t1**4*t3**2*t2**9-4._ki*t1**5*t5**3*t2**&
                  &8+50._ki*t1**4*t3**4*t2**5+35._ki/6._ki*t3**4*t2**9*t1**2+16._ki*t1&
                  &**5*t3**3*t2**5-64._ki*t2**3*t3**5*t4**8-t2**11*t1**2*t3**3/3._ki&
                  &-760._ki/3._ki*t2**5*t3**5*t4**6+5._ki/2._ki*t3**5*t4*t2**10-t1**4*&
                  &t3*t2**11/2._ki+242._ki*t4**5*t3**5*t2**6-33._ki/2._ki*t4**2*t3**5*&
                  &t2**9-152._ki*t6*t1*t4**4*t3**4*t2**6-760._ki/3._ki*t6*t1*t2**4*t3&
                  &**4*t4**6+377._ki/6._ki*t6*t1*t4**3*t3**4*t2**7
                !
                stemp7=stemp11+121._ki/2._ki*t6**2*t4**2*t2**6*t1**5*t5-605._ki/2._ki&
                  &*t6**2*t4**2*t5**3*t2**4*t1**5+165._ki/2._ki*t6**2*t4*t2**6*t3*t1&
                  &**5-765._ki/4._ki*t6**2*t4**2*t5**2*t2**5*t1**4*t3+219._ki/4._ki*t6&
                  &**2*t4**2*t5*t2**6*t1**3*t3**2-605._ki*t6**2*t1**5*t4**4*t5**3*t&
                  &2**2-400._ki*t6*t1**3*t4**6*t3**2*t5**2*t2**2+80._ki*t6*t1**3*t4*&
                  &*7*t3**2*t5**2*t2+825._ki/4._ki*t6**2*t1**4*t3**3*t2**4*t4+20._ki*&
                  &t6**2*t1**4*t4**6*t3*t2**3+48._ki*t6**2*t1**4*t4**4*t3*t2**5+60.&
                  &_ki*t6**2*t1**4*t4**7*t3*t5**2-720._ki*t6**2*t1**4*t4**4*t3**3*t2&
                  &-219._ki/2._ki*t6**2*t1**3*t4**2*t3**3*t2**5-208._ki*t6**2*t1**3*t&
                  &4**6*t3**3*t2+1155._ki/4._ki*t6**2*t1**3*t3**3*t4**3*t2**4-1935._k&
                  &i*t6**2*t1**4*t4**3*t5*t3**2*t2**3+1935._ki/4._ki*t6**2*t1**4*t4*&
                  &*3*t5**2*t2**4*t3-720._ki*t6**2*t1**4*t4**4*t5**2*t2**3*t3+275._k&
                  &i*t1**3*t4**4*t3**4*t2**2*t5-1100._ki*t1**3*t4**4*t3**3*t5**2*t2&
                  &**3+300._ki*t1**3*t4**5*t3**3*t5**2*t2**2
                !
                stemp11=stemp7-525._ki*t1**3*t4**3*t3**2*t5**3*t2**5+440._ki*t6*t4*&
                  &*2*t5**3*t2**5*t1**4*t3-255._ki/2._ki*t6*t4**2*t5*t3**2*t2**7*t1*&
                  &*3-1100._ki*t1**3*t4**2*t3**3*t5**2*t2**5+377._ki/6._ki*t4**3*t3**&
                  &5*t2**8-630._ki*t1**3*t4**3*t3**2*t5**2*t2**6+1260._ki*t1**3*t4**&
                  &3*t3**3*t2**5*t5-880._ki*t1**3*t4**2*t3**3*t2**6*t5+660._ki*t1**3&
                  &*t4**4*t3**2*t5**2*t2**5-360._ki*t1**3*t4**5*t3**2*t5**2*t2**4-1&
                  &26._ki*t1**3*t3**2*t4**3*t5*t2**7+16._ki*t1**3*t3**2*t4**6*t5*t2*&
                  &*4-1260._ki*t6**2*t1**4*t4**5*t3**2*t5*t2-42._ki*t6**2*t1**4*t4**&
                  &5*t3*t2**4-66._ki*t6**2*t1**5*t4**5*t5*t2**3+66._ki*t6**2*t1**5*t&
                  &4**5*t3*t2**2+300._ki*t6**2*t1**4*t4**6*t3**2*t5-129._ki/4._ki*t6*&
                  &*2*t1**4*t4**3*t3*t2**6+1935._ki/2._ki*t6**2*t1**4*t4**3*t3**3*t2&
                  &**2
                !
                stemp10=stemp11-5._ki*t6*t3**2*t5**2*t2**8*t1**3+4._ki/3._ki*t6*t3**&
                  &3*t5*t2**9*t1**2+40._ki/3._ki*t6*t3*t5**3*t2**7*t1**4+10._ki*t6*t3&
                  &*t2**8*t1**4*t5**2-10395._ki/2._ki*t6**2*t1**5*t4**3*t3*t5**2*t2*&
                  &*2+3630._ki*t6**2*t1**5*t4**4*t3*t5**2*t2-105._ki*t6**2*t1**6*t2*&
                  &*4*t3*t4-40._ki*t6**2*t1**6*t4**4*t5*t2**2+1575._ki*t6**2*t1**6*t&
                  &4*t3*t5**2*t2**2+1425._ki/2._ki*t6**2*t1**6*t4**2*t5**3*t2**2-625&
                  &._ki*t6**2*t1**6*t4**3*t5**3*t2-125._ki/2._ki*t6**2*t1**6*t2**2*t3&
                  &*t4**3+121._ki*t6**2*t1**5*t4**4*t5*t2**4+330._ki*t6**2*t1**5*t4*&
                  &*5*t5**3*t2-75._ki*t6**2*t1**7*t4*t3*t5**2-800._ki/3._ki*t1**2*t4*&
                  &*6*t3**3*t5*t2**4-700._ki*t1**2*t4**5*t3**4*t5*t2**3-640._ki*t1**&
                  &2*t4**4*t3**3*t5*t2**6-1075._ki*t1**2*t4**3*t3**4*t5*t2**5-1000.&
                  &_ki/3._ki*t1**2*t4**6*t3**3*t5**2*t2**3
                !
                stemp11=stemp10+365._ki/4._ki*t4**2*t5*t2**8*t1*t3**4-170._ki*t4**2*&
                  &t5*t3**3*t2**8*t1**2-1425._ki*t1**4*t4**2*t3**3*t2**3*t5**2+16._k&
                  &i*t1**4*t4**4*t3*t2**6*t5-350._ki*t1**4*t3**4*t5*t2**3*t4+1050._k&
                  &i*t1**4*t4*t3**3*t2**4*t5**2-25._ki*t1**4*t4**3*t3*t2**7*t5-4._ki&
                  &*t1**4*t4**5*t3*t2**5*t5+500._ki*t1**4*t4**3*t3**3*t5*t2**3-1575&
                  &._ki/2._ki*t1**3*t4**3*t3**4*t2**3*t5-375._ki*t1**3*t4*t3**4*t2**5&
                  &*t5-72._ki*t1**3*t3**2*t4**5*t5*t2**5+1575._ki*t1**3*t4**3*t3**3*&
                  &t5**2*t2**4-540._ki*t1**5*t4**2*t3**2*t5**2*t2**3-2475._ki/2._ki*t&
                  &6**2*t4*t3*t2**4*t1**5*t5**2-832._ki/3._ki*t6*t1**2*t4**6*t3**4*t&
                  &2**2-2565._ki*t6*t1**5*t4**2*t3*t5**2*t2**4+475._ki/2._ki*t1**4*t4&
                  &**2*t3**4*t5*t2**2+1710._ki*t1**4*t4**2*t3**2*t5**2*t2**5+500._ki&
                  &*t1**4*t4**3*t5**3*t3*t2**5
                !
                stemp9=stemp11-320._ki*t1**4*t4**4*t5**3*t3*t2**4+150._ki*t1**4*t4*&
                  &*3*t3*t2**6*t5**2+24._ki*t1**4*t4**5*t3*t2**4*t5**2+420._ki*t6*t1&
                  &**3*t4**5*t3**2*t5*t2**4+840._ki*t1**4*t4*t3**3*t5*t2**5+75._ki/4&
                  &._ki*t6*t1**5*t4**3*t5*t2**6+14._ki*t1**5*t4**3*t5*t3*t2**5-36._ki&
                  &*t1**5*t4**2*t5*t3*t2**6+840._ki*t6*t1**3*t4**5*t3**2*t5**2*t2**&
                  &3+645._ki/2._ki*t6*t1**3*t4**3*t5*t3**2*t2**6-75._ki/2._ki*t6*t1**5&
                  &*t4**3*t3*t2**5-375._ki*t6*t1**5*t4**3*t5**4*t2**3-60._ki*t6*t1**&
                  &5*t4**5*t5**4*t2+2475._ki*t6**2*t1**5*t3**2*t2**3*t4*t5+693._ki/2&
                  &._ki*t6**2*t1**5*t4**3*t3*t2**4-1815._ki*t6**2*t1**5*t4**4*t3**2*&
                  &t5-242._ki*t6**2*t1**5*t4**2*t3*t2**5+1210._ki*t6**2*t1**5*t3**3*&
                  &t4**2*t2-800._ki*t1**2*t4**4*t3**3*t5**2*t2**5-1200._ki*t6*t1**4*&
                  &t4*t3**3*t5*t2**4+208._ki*t1*t4**6*t3**4*t2**5
                !
                stemp11=stemp9+216._ki*t1**5*t4**2*t5**2*t2**5*t3+180._ki*t1**5*t4*&
                  &t5*t2**5*t3**2+450._ki*t1**5*t4**2*t5**4*t3*t2**3-280._ki*t1**5*t&
                  &4**3*t3*t5**3*t2**3-300._ki*t1**3*t4**5*t3**2*t5**3*t2**3+200._ki&
                  &/3._ki*t1**3*t4**6*t3**2*t5**3*t2**2+550._ki*t1**3*t4**4*t3**2*t5&
                  &**3*t2**4+176._ki/3._ki*t6*t1**2*t4**7*t3**4*t2+520._ki*t1*t4**6*t&
                  &3**4*t5*t2**4+44._ki*t1*t2**2*t3**5*t4**7-219._ki/2._ki*t1*t2**7*t&
                  &3**5*t4**2+16._ki*t1*t4**8*t3**4*t2**3-7._ki*t1**5*t4**3*t3*t2**6&
                  &-20._ki*t1**5*t3**3*t2**4*t4-55._ki*t4*t3**4*t2**8*t1**2+28._ki*t1&
                  &**5*t4**3*t5**3*t2**5-140._ki*t1**4*t4*t3**4*t2**4+125._ki*t1**4*&
                  &t3**4*t5*t2**4+9._ki*t1**5*t4**2*t3**2*t2**5
                !
                stemp10=stemp11-114._ki*t1**4*t4**2*t3**3*t2**5+50._ki*t1**4*t3**3*&
                  &t4**3*t2**4+2._ki*t1**5*t4**4*t5**2*t2**5-110._ki*t1**3*t3**5*t2*&
                  &*3*t4**2-10._ki/3._ki*t3**3*t5*t2**10*t1**2+25._ki/3._ki*t3**2*t5**&
                  &3*t2**8*t1**3-315._ki*t1**3*t4**3*t3**4*t2**4-4._ki/3._ki*t1**3*t4&
                  &**6*t3**2*t2**5+110._ki*t1**3*t4**4*t3**4*t2**3-88._ki*t1**3*t4**&
                  &2*t3**3*t2**7+105._ki/2._ki*t1**3*t3**5*t2**2*t4**3+21._ki/2._ki*t1&
                  &**3*t4**3*t3**2*t2**8+6._ki*t1**3*t4**5*t3**2*t2**6+330._ki*t1**3&
                  &*t4**2*t3**4*t2**5+126._ki*t1**3*t4**3*t3**3*t2**6-88._ki*t1**3*t&
                  &4**4*t3**3*t2**5+24._ki*t1**3*t4**5*t3**3*t2**4+11._ki/6._ki*t6**2&
                  &*t2**8*t1**5*t5-55._ki/6._ki*t6**2*t5**3*t2**6*t1**5-25._ki/2._ki*t&
                  &6**2*t1**6*t5*t2**6-8._ki/3._ki*t6*t2**8*t1**2*t3**4
                !
                stemp11=stemp10+5._ki*t6*t2**8*t3**2*t1**4+10._ki*t6*t1**6*t5**2*t2&
                  &**6+2._ki*t6*t1**6*t5*t2**7-40._ki*t6*t1**6*t5**4*t2**4-40._ki*t6*&
                  &t1**6*t5**3*t2**5-21._ki/8._ki*t6*t1**5*t4*t2**9+15._ki/4._ki*t6*t2&
                  &**8*t3*t1**5+35._ki/6._ki*t6*t3**3*t2**8*t1**3-30._ki*t6**2*t1**7*&
                  &t5**3*t2**2+6._ki*t6**2*t1**7*t5*t2**4+25._ki*t6**2*t1**6*t3*t2**&
                  &5-105._ki/4._ki*t6**2*t3**3*t2**5*t1**4-3._ki*t6*t1**6*t3*t2**6+15&
                  &._ki/4._ki*t6*t1**6*t4*t2**7-2._ki*t6**2*t2**7*t3**3*t1**3-t6**2*t&
                  &3**3*t2**9*t1**2/12._ki-11._ki*t6**2*t2**7*t3*t1**5+125._ki/2._ki*t&
                  &6**2*t1**6*t5**3*t2**4-50._ki*t3**3*t5**2*t2**7*t1**3+11._ki/3._ki&
                  &*t4*t2**10*t1**2*t3**3
                !
                stemp8=stemp11-8._ki*t1**5*t4**4*t5**3*t2**4+30._ki*t4*t2**8*t3**3*&
                  &t1**3-13._ki/2._ki*t4*t3**4*t2**10*t1-25._ki/2._ki*t1**4*t5**4*t2**&
                  &7*t3-6._ki*t1**4*t3*t2**9*t5**2-20._ki*t1**4*t5**3*t3*t2**8+95._ki&
                  &*t1**4*t4**2*t3**4*t2**3-9._ki/2._ki*t6*t1**6*t4**2*t2**6-200._ki/&
                  &3._ki*t6*t1**4*t3**4*t2**4+10._ki*t6*t1**6*t3**2*t2**4-4._ki*t6**2&
                  &*t1**7*t3*t2**3-45._ki*t6**2*t1**7*t4**2*t5**3+175._ki/2._ki*t6**2&
                  &*t1**6*t3**3*t4-1155._ki/2._ki*t6**2*t1**5*t4**3*t3**3-125._ki/2._k&
                  &i*t6**2*t1**6*t3**3*t2+550._ki/3._ki*t6**2*t1**5*t2**3*t3**3-220.&
                  &_ki/3._ki*t6**2*t1**5*t4**6*t5**3+210._ki*t6**2*t1**4*t4**5*t3**3+&
                  &44._ki*t6**2*t1**3*t4**7*t3**3+16._ki/3._ki*t6**2*t1**2*t4**9*t3**&
                  &3+132._ki*t1**3*t3**2*t4**4*t5*t2**6-600._ki*t6*t1**5*t3**3*t5*t2&
                  &**3
                !
                stemp11=stemp8+375._ki/4._ki*t6*t1**5*t4**3*t3**2*t2**3+40._ki*t6*t1&
                  &**3*t4**7*t3**2*t5*t2**2-480._ki*t6*t1**3*t4**4*t3**2*t5*t2**5-2&
                  &00._ki*t6*t1**3*t4**6*t3**2*t5*t2**3-430._ki*t6*t1**3*t4**3*t3**3&
                  &*t2**5-11._ki*t6*t1**4*t4**4*t3*t2**6+210._ki*t6*t1**4*t4**3*t3**&
                  &4*t2+50._ki/3._ki*t6*t1**3*t4**6*t3**2*t2**4+200._ki/3._ki*t6*t1**3&
                  &*t4**6*t3**3*t2**2-40._ki/3._ki*t6*t1**4*t4**6*t3*t2**3*t5-480._ki&
                  &*t6*t1**4*t4**5*t5**3*t3*t2**2+880._ki*t6*t1**4*t4**4*t3**3*t5*t&
                  &2-480._ki*t6*t1**3*t4**4*t3**4*t2**2+416._ki/3._ki*t6*t1**2*t4**6*&
                  &t3**3*t2**4+377._ki/12._ki*t6**2*t1**2*t4**3*t3**3*t2**6+66._ki*t4&
                  &**2*t3**2*t2**8*t1**3*t5+96._ki*t1**4*t4**4*t5*t2**4*t3**2-300._k&
                  &i*t1**4*t4**3*t5*t2**5*t3**2+270._ki*t6*t1**6*t4**2*t3*t5**2*t2*&
                  &*2
                !
                stemp10=stemp11+75._ki*t6*t1**6*t4*t5*t2**4*t3+5130._ki*t6*t1**5*t4&
                  &**2*t3**2*t5**2*t2**2-2250._ki*t6*t1**5*t4**3*t3**2*t5**2*t2+225&
                  &0._ki*t6*t1**5*t4**3*t3*t5**2*t2**3-720._ki*t6*t1**5*t4**4*t3*t5*&
                  &*2*t2**2-375._ki*t6*t1**5*t4**3*t5**3*t2**4+21._ki/4._ki*t6*t1**5*&
                  &t4*t2**8*t5+3._ki/2._ki*t6*t1**4*t4*t3*t2**9-60._ki*t6*t1**5*t4**4&
                  &*t5**2*t2**4+292._ki/3._ki*t6*t4**2*t5*t2**7*t1**2*t3**3+552._ki*t&
                  &6*t1**2*t4**5*t3**4*t2**3+660._ki*t6*t1**4*t4**4*t3*t5**2*t2**4-&
                  &15._ki*t6**2*t1**7*t4*t5*t2**3+75._ki*t6**2*t1**7*t4*t5**3*t2-427&
                  &5._ki/2._ki*t6**2*t1**6*t4**2*t3*t5**2*t2-990._ki*t6**2*t1**5*t4**&
                  &5*t3*t5**2-825._ki*t6**2*t1**5*t3**3*t2**2*t4-242._ki*t6**2*t1**5&
                  &*t4**4*t3*t2**3-231._ki/2._ki*t6**2*t1**5*t4**3*t5*t2**5+450._ki*t&
                  &6*t4*t3**2*t2**6*t1**4*t5-255._ki*t6*t4**2*t5**2*t2**6*t1**3*t3*&
                  &*2
                !
                stemp11=stemp10-700._ki*t1**4*t4*t5**3*t2**5*t3**2-380._ki*t1**4*t4&
                  &**2*t3*t5**3*t2**6-114._ki*t1**4*t4**2*t5**2*t2**7*t3-475._ki/2._k&
                  &i*t1**4*t4**2*t3*t5**4*t2**5-175._ki*t1**5*t4**3*t3*t5**4*t2**2+&
                  &900._ki*t1**5*t4*t3**2*t5**2*t2**4-1140._ki*t1**4*t4**2*t3**3*t5*&
                  &t2**4+1155._ki/2._ki*t6**2*t1**5*t4**3*t5**3*t2**3+165._ki/2._ki*t6&
                  &**2*t4*t5**3*t2**5*t1**5-11._ki/4._ki*t6**2*t4*t2**8*t1**4*t3+850&
                  &._ki*t6*t1**3*t4**2*t3**3*t5*t2**5-220._ki*t1*t4**7*t3**4*t5*t2**&
                  &3+200._ki/3._ki*t1**2*t4**7*t3**3*t5**2*t2**2+560._ki*t1**2*t4**5*&
                  &t3**3*t5*t2**5+700._ki*t1**2*t4**5*t3**3*t5**2*t2**4-75._ki*t4*t5&
                  &**3*t2**7*t1**3*t3**2-425._ki/2._ki*t4**2*t5**2*t2**7*t1**2*t3**3&
                  &+375._ki*t6**2*t1**6*t3**2*t5*t2**2-108._ki*t1**5*t4**2*t3**2*t5*&
                  &t2**4-450._ki*t6**2*t1**3*t4**4*t3**3*t2**3
                !
                stemp9=stemp11-84._ki*t1**5*t4**3*t5**2*t2**4*t3-250._ki*t1**5*t4*t&
                  &3**3*t5**2*t2**2-450._ki*t1**5*t4**2*t5**3*t3**2*t2**2-200._ki*t1&
                  &**5*t4*t5*t3**3*t2**3+750._ki*t1**5*t4*t5**3*t3**2*t2**3+825._ki*&
                  &t1**3*t4**2*t3**4*t2**4*t5-1275._ki/2._ki*t6**2*t1**4*t4**2*t3**3&
                  &*t2**3-t6*t1**6*t2**8-880._ki*t1**3*t4**4*t3**3*t2**4*t5+240._ki*&
                  &t1**3*t4**5*t3**3*t2**3*t5+70._ki/3._ki*t6*t2**7*t3**3*t1**3*t5-1&
                  &1._ki/2._ki*t6*t1**4*t4**2*t3*t2**8+91._ki/3._ki*t6*t4*t3**4*t2**7*&
                  &t1**2+1680._ki*t6*t1**5*t4*t5**3*t2**4*t3+85._ki/8._ki*t6*t4**2*t3&
                  &**2*t2**8*t1**3+1260._ki*t6*t1**5*t4*t3*t5**2*t2**5-57._ki/4._ki*t&
                  &6*t1**5*t4**2*t5*t2**7-21._ki*t6*t1**5*t4*t3*t2**7+450._ki*t6*t1*&
                  &*5*t2**5*t3**2*t5-285._ki/4._ki*t6*t1**5*t4**2*t5**2*t2**6-5._ki/3&
                  &._ki*t6*t2**9*t1**4*t3*t5
                !
                stemp11=stemp9-184._ki*t6*t1**2*t4**5*t3**3*t2**5-1250._ki*t1**4*t4&
                  &**3*t3**2*t5**3*t2**3+80._ki*t1**4*t4**5*t5**3*t3*t2**3+50._ki*t1&
                  &**4*t4**5*t5**4*t3*t2**2+625._ki/2._ki*t1**4*t4**3*t5**4*t3*t2**4&
                  &-1500._ki*t1**4*t4**3*t3**2*t5**2*t2**4-96._ki*t1**4*t4**4*t3*t2*&
                  &*5*t5**2+10._ki*t6*t1**6*t3**3*t2**2-3._ki/2._ki*t6*t1**5*t4**5*t2&
                  &**5+6._ki*t6*t1**5*t4**4*t2**6-75._ki/8._ki*t6*t1**5*t4**3*t2**7+7&
                  &5._ki*t6*t1**5*t3**4*t2**2-150._ki*t6*t1**5*t3**3*t2**4+7._ki/4._ki&
                  &*t6*t1**6*t4**3*t2**5-15._ki*t6**2*t1**7*t3**2*t5+200._ki*t6**2*t&
                  &1**6*t4**4*t5**3-t6*t3**4*t2**10*t1/6._ki+57._ki/8._ki*t6*t1**5*t4&
                  &**2*t2**8-75._ki/2._ki*t6*t1**5*t3**2*t2**6
                !
                stemp10=stemp11-35._ki/2._ki*t6*t3**4*t2**6*t1**3+50._ki*t6*t3**3*t2&
                  &**6*t1**4-15._ki/4._ki*t6*t1**5*t5**2*t2**8+15._ki*t6*t1**5*t5**4*&
                  &t2**6+15._ki*t6*t1**5*t5**3*t2**7-3._ki/4._ki*t6*t1**5*t2**9*t5-t6&
                  &*t1**4*t3*t2**10/6._ki+5._ki/24._ki*t6*t2**10*t1**3*t3**2+t6*t3**3&
                  &*t2**10*t1**2/3._ki+t6**2*t2**9*t1**4*t3/4._ki-300._ki*t1**5*t3**2&
                  &*t5**3*t2**4-15._ki*t1**5*t4*t3*t2**8-90._ki*t1**5*t4**2*t5**4*t2&
                  &**5-72._ki*t1**5*t5*t3**2*t2**6-36._ki*t1**5*t4**2*t5**3*t2**6-25&
                  &0._ki*t1**4*t3**3*t5**2*t2**5+2._ki*t2**10*t1**3*t3**2*t5+10._ki*t&
                  &5**2*t3**2*t2**9*t1**3-17._ki*t4**2*t3**3*t2**9*t1**2-200._ki*t1*&
                  &*4*t3**3*t5*t2**6+14._ki*t1**4*t4*t3**2*t2**8
                !
                stemp11=stemp10-19._ki/2._ki*t1**4*t4**2*t3*t2**9-11._ki/2._ki*t1**3*&
                  &t4**2*t3**2*t2**9-40._ki*t5*t3**3*t2**8*t1**3+50._ki*t1**5*t4*t5*&
                  &*4*t2**6+5._ki*t1**5*t4*t2**9*t5+20._ki*t1**5*t4*t5**3*t2**7-5._ki&
                  &*t1**5*t2**8*t4*t5**2+25._ki*t1**5*t4*t5**5*t2**5+48._ki*t1**5*t3&
                  &*t5**2*t2**7+160._ki*t1**5*t3*t5**3*t2**6+100._ki*t1**5*t5**4*t2*&
                  &*5*t3-8._ki*t1**5*t5*t3*t2**8+175._ki/12._ki*t2**8*t3**4*t1**2*t5+&
                  &91._ki/4._ki*t4*t3**5*t2**8*t1+125._ki/2._ki*t5*t2**6*t3**4*t1**3-4&
                  &5._ki*t1**5*t4**2*t5**5*t2**4-9._ki*t1**5*t4**2*t2**8*t5-360._ki*t&
                  &1**5*t3**2*t2**5*t5**2+9._ki*t1**5*t4**2*t5**2*t2**7+t1**4*t3*t2&
                  &**10*t5
                !
                stemp6=stemp11+645._ki*t6*t1**3*t4**3*t3**4*t2**3+140._ki*t6*t1**3*&
                  &t4**5*t3**4*t2+275._ki/2._ki*t6*t1**3*t3**4*t4*t2**5-10._ki/3._ki*t&
                  &6*t1**3*t4**7*t3**2*t2**3-280._ki*t6*t1**3*t4**5*t3**3*t2**3+480&
                  &._ki*t6*t1**3*t4**4*t3**3*t2**4+168._ki*t6*t1*t2**3*t3**4*t4**7+3&
                  &2._ki/3._ki*t6*t1*t2*t3**4*t4**9-2520._ki*t6*t1**4*t4**3*t3**3*t2*&
                  &*2*t5-1320._ki*t6*t1**4*t4**2*t3**2*t5*t2**5-360._ki*t6*t1**4*t4*&
                  &*5*t3*t5**2*t2**3+360._ki*t6*t1**6*t4**2*t5**3*t2*t3-45._ki*t6*t1&
                  &**6*t4**2*t5*t2**3*t3-105._ki*t6*t1**5*t4*t5**3*t2**6+105._ki/4._k&
                  &i*t6*t1**5*t4*t5**2*t2**7+75._ki/2._ki*t6*t1**5*t5*t2**7*t3+105._k&
                  &i/4._ki*t6**2*t2**6*t3**2*t1**4*t5-375._ki*t6**2*t1**6*t3*t5**2*t&
                  &2**3-350._ki*t6**2*t1**6*t4*t5**3*t2**3+70._ki*t6**2*t1**6*t4*t5*&
                  &t2**5-825._ki/2._ki*t6**2*t5*t2**4*t3**2*t1**5+91._ki/4._ki*t6**2*t&
                  &4*t3**3*t2**6*t1**3
                !
                stemp7=t6/t1**5/t2**12
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              case(4)
                !
                stemp5=-t4*t2**8*t5/6._ki+t4**2*t5*t2**7/12._ki+25._ki/4._ki*t6**3*t4&
                  &**2*t2**5+t4*t3*t2**7/6._ki-13._ki/2._ki*t6**2*t4**2*t2**6-95._ki/3&
                  &._ki*t6**3*t3**2*t4**3+7._ki/12._ki*t6**2*t5*t2**7-2._ki*t6**3*t4*t&
                  &2**6-t6*t4**3*t2**6/2._ki-2._ki/3._ki*t6**3*t5*t2**6+25._ki/3._ki*t6&
                  &**3*t2**3*t3**2+10._ki/3._ki*t6**3*t2**5*t3+6._ki*t6**2*t4**3*t2**&
                  &5-24._ki*t6**3*t1*t4**3*t5*t2-26._ki*t6**3*t1*t2*t3*t4**2+8._ki*t6&
                  &**3*t1**2*t5*t2*t4-35._ki/6._ki*t6**2*t1*t4**2*t5*t2**3+28._ki/3._k&
                  &i*t6**2*t1*t4*t5*t2**4+65._ki*t6**2*t4**2*t5**2*t3*t2**2-30._ki*t&
                  &6**2*t4**3*t5**2*t3*t2+7._ki/12._ki*t6**2*t1**2*t5*t2**3-5._ki/6._k&
                  &i*t6**2*t1**2*t2*t5**3+7._ki/2._ki*t6**2*t1*t2**4*t3-5._ki*t6**2*t&
                  &2*t3**3*t4-15._ki*t6**2*t2**3*t3**2*t5-10._ki*t6**3*t1*t3**2*t2+1&
                  &5._ki*t6**3*t3**2*t1*t4+125._ki/2._ki*t6**3*t3**2*t2*t4**2-40._ki*t&
                  &6**3*t2**2*t3**2*t4
                !
                stemp4=stemp5-39._ki/2._ki*t6**3*t1*t4**2*t2**3+27._ki/2._ki*t6**3*t1&
                  &*t4*t2**4+30._ki*t6**3*t1*t4**3*t5**2-12._ki*t6**3*t1*t2**3*t3+7.&
                  &_ki*t6**2*t4**3*t3*t2**3+21._ki/2._ki*t6**2*t4*t3*t2**5+8._ki*t6**3&
                  &*t1*t5*t2**4-10._ki*t6**3*t1*t5**2*t2**3-7._ki/2._ki*t6**2*t2**6*t&
                  &4*t5+5._ki*t6**2*t4*t5**3*t2**4-7._ki/2._ki*t6**2*t1*t5*t2**5+5._ki&
                  &*t6**2*t1*t5**3*t2**3+10._ki*t6**2*t3*t5**2*t2**4+10._ki*t6**2*t4&
                  &**3*t5**3*t2**2-10._ki/3._ki*t6**2*t4**4*t5**3*t2-65._ki/6._ki*t6**&
                  &2*t4**2*t5**3*t2**3+16._ki/3._ki*t6**3*t4*t5*t2**5+5._ki*t6**2*t1*&
                  &t4**2*t2**4-8._ki*t6**2*t1*t4*t2**5+t6*t4*t2**6*t1/2._ki+10._ki/3.&
                  &_ki*t6**2*t2**2*t3**3-t6*t4*t2**8-20._ki/3._ki*t6**3*t4**5*t5**2-1&
                  &9._ki/2._ki*t6**3*t4**3*t2**4+9._ki/4._ki*t6**3*t1**2*t2**3-2._ki*t6&
                  &**3*t4**5*t2**2-2._ki*t6**2*t4**4*t2**4+7._ki*t6**3*t4**4*t2**3-t&
                  &1**2*t2**4*t6**2/2._ki
                !
                stemp5=stemp4+5._ki/4._ki*t6*t4**2*t2**7-7._ki/3._ki*t6**2*t3*t2**6+3&
                  &._ki*t6**2*t4*t2**7+3._ki*t6**2*t1*t2**6-t6*t1*t2**7/2._ki-3._ki*t6&
                  &**3*t1*t2**5+5._ki/6._ki*t6**3*t5**2*t2**5-5._ki/6._ki*t6**2*t5**3*&
                  &t2**5+20._ki*t6**2*t1*t2*t3*t4*t5**2+5._ki*t6**2*t1*t2*t5*t3**2-1&
                  &4._ki/3._ki*t6**2*t1*t2**3*t3*t4-15._ki*t6**2*t1*t5**2*t3*t2**2-65&
                  &._ki/2._ki*t6**2*t2*t3**2*t4**2*t5+45._ki*t6**2*t3**2*t4*t5*t2**2+&
                  &160._ki/3._ki*t6**3*t4*t3*t5*t2**3+30._ki*t6**3*t1*t3*t5*t2**2+65.&
                  &_ki*t6**3*t1*t4**2*t3*t5-t3*t2**8/6._ki+t2**9*t5/12._ki-t6**2*t2**&
                  &8/2._ki+t6*t2**9/4._ki+t6**3*t2**7/4._ki-45._ki*t6**2*t4*t5**2*t3*t&
                  &2**3-40._ki/3._ki*t6**2*t1*t4*t5**3*t2**2+25._ki/3._ki*t6**2*t1*t4*&
                  &*2*t5**3*t2-65._ki*t6**3*t1*t4**2*t5**2*t2+36._ki*t6**3*t1*t2**2*&
                  &t3*t4+52._ki*t6**3*t1*t4**2*t5*t2**2
                !
                stemp3=stemp5-36._ki*t6**3*t1*t4*t5*t2**3-125._ki*t6**3*t4**2*t3*t5&
                  &*t2**2+380._ki/3._ki*t6**3*t4**3*t3*t5*t2-90._ki*t6**3*t1*t3*t5*t2&
                  &*t4+45._ki*t6**3*t1*t4*t5**2*t2**2-5._ki*t6**3*t1**2*t3*t5-3._ki*t&
                  &6**3*t1**2*t2**2*t4+56._ki/3._ki*t6**3*t4**4*t3*t2+50._ki*t6**3*t4&
                  &**2*t3*t2**3+70._ki/3._ki*t6**3*t4**4*t5**2*t2-152._ki/3._ki*t6**3*&
                  &t4**3*t3*t2**2-64._ki/3._ki*t6**3*t4*t3*t2**4-95._ki/3._ki*t6**3*t4&
                  &**3*t5**2*t2**2+9._ki*t6**3*t1*t4**3*t2**2+16._ki/3._ki*t6**3*t4**&
                  &5*t5*t2+125._ki/6._ki*t6**3*t4**2*t5**2*t2**3+76._ki/3._ki*t6**3*t4&
                  &**3*t5*t2**3-56._ki/3._ki*t6**3*t4**4*t5*t2**2-50._ki/3._ki*t6**3*t&
                  &4**2*t5*t2**4-140._ki/3._ki*t6**3*t4**4*t3*t5-10._ki*t6**3*t1**2*t&
                  &5**2*t4+15._ki/2._ki*t6**3*t1**2*t5**2*t2+2._ki*t6**3*t1**2*t2*t3-&
                  &25._ki/3._ki*t6**3*t3*t2**4*t5-20._ki/3._ki*t6**3*t4*t5**2*t2**4-6.&
                  &_ki*t6**3*t1**2*t5*t2**2-91._ki/6._ki*t6**2*t4**2*t3*t2**4+91._ki/1&
                  &2._ki*t6**2*t4**2*t5*t2**5-7._ki*t6**2*t4**3*t5*t2**4+7._ki/3._ki*t&
                  &6**2*t4**4*t5*t2**3
                !
                stemp4=1._ki/t2**10*z_log(t1*t6/t2**2,1._ki)
                !
                stemp2=stemp3*stemp4
                !
                stemp4=-(-t4+t2)*(t3*t5*t2**2+3._ki*t3*t6*t2**2-t3*t5*t2*t4-4._ki*t&
                  &5*t2*t1*t6-9._ki*t6*t4*t2*t3-2._ki*t3**2*t2+6._ki*t4**2*t3*t6+4._ki&
                  &*t4*t1*t6*t5+2._ki*t3*t1*t6)/t3/t2**4*q(4,(t2*t3-t1*t6)/t2/t3,si&
                  &gn_arg)/12._ki
                !
                stemp9=-125._ki/6._ki*t6**2*t1**2*t3**3*t4**2*t2**2-40._ki/3._ki*t6**&
                  &2*t1*t4**5*t3**2*t5*t2**2+10._ki/3._ki*t6**2*t1*t4**6*t3**2*t5*t2&
                  &+65._ki/3._ki*t6**2*t1*t4**4*t3**2*t5*t2**3+55._ki/3._ki*t6**2*t1*t&
                  &2**3*t3**3*t4**3+195._ki/2._ki*t6**2*t1**3*t4**2*t5**2*t3*t2**2-4&
                  &5._ki*t6**2*t1**3*t4**3*t5**2*t3*t2-195._ki/4._ki*t6**2*t1**3*t3**&
                  &2*t5*t2*t4**2+135._ki/2._ki*t6**2*t1**3*t3**2*t5*t2**2*t4+21._ki/2&
                  &._ki*t6**2*t1**3*t4**3*t3*t2**3-35._ki/12._ki*t6**2*t1**4*t4**2*t5&
                  &*t2**3-20._ki/3._ki*t6**2*t1**4*t4*t5**3*t2**2-7._ki/3._ki*t6**2*t1&
                  &**4*t2**3*t3*t4+10._ki*t6**2*t1**4*t5**2*t3*t2*t4+5._ki/2._ki*t6**&
                  &2*t1**4*t3**2*t5*t2-15._ki/2._ki*t6**2*t1**4*t5**2*t3*t2**2-25._ki&
                  &/12._ki*t6**2*t1*t4*t3**2*t2**6*t5-14._ki/9._ki*t6**2*t1**2*t4*t3*&
                  &t2**7+25._ki/6._ki*t6**2*t1**2*t2**5*t3**2*t5+125._ki/12._ki*t6**3*&
                  &t1**2*t4*t3**2*t2**4-21._ki/4._ki*t6**2*t1**3*t2**6*t4*t5+205._ki/&
                  &12._ki*t6**3*t1**2*t4**2*t5*t2**4*t3
                !
                stemp8=stemp9-25._ki/3._ki*t6**3*t1**4*t5**2*t2**3-7._ki/2._ki*t6**2*&
                  &t1**3*t3*t2**6+5._ki/2._ki*t6**2*t1**4*t5**3*t2**3-7._ki/4._ki*t6**&
                  &2*t1**4*t5*t2**5+9._ki/2._ki*t6**2*t1**3*t4*t2**7-t6*t1**3*t4*t2*&
                  &*8-61._ki/36._ki*t6**2*t2**6*t3**3*t4**2-5._ki/12._ki*t6**2*t1*t2**&
                  &6*t3**3+85._ki/18._ki*t6**2*t2**5*t3**3*t4**3+55._ki/36._ki*t6**3*t&
                  &1**3*t5**2*t2**5-11._ki/9._ki*t6**3*t1**3*t5*t2**6+55._ki/9._ki*t6*&
                  &*3*t1**3*t2**5*t3-11._ki/3._ki*t6**3*t1**3*t4*t2**6-t6**3*t1**2*t&
                  &2**7*t3/6._ki-5._ki/4._ki*t6**2*t1**3*t5**3*t2**5+7._ki/8._ki*t6**2*&
                  &t1**3*t5*t2**7+63._ki/4._ki*t6**2*t1**3*t4*t3*t2**5-15._ki/2._ki*t6&
                  &**2*t1**3*t3**3*t2*t4+95._ki/3._ki*t6**2*t1**2*t4**3*t3*t5**2*t2*&
                  &*3-125._ki/6._ki*t6**2*t1**2*t4**2*t3*t5**2*t2**4-80._ki/3._ki*t6**&
                  &2*t1**2*t4*t3**2*t5*t2**4+40._ki/3._ki*t6**2*t1**2*t4*t3**3*t2**3&
                  &+175._ki/36._ki*t6**2*t1**2*t4**2*t3*t2**6
                !
                stemp9=stemp8+49._ki/9._ki*t6**2*t1**2*t4**4*t3*t2**4+3._ki/2._ki*t6*&
                  &*2*t1**4*t2**6+3._ki/4._ki*t6**3*t1**5*t2**3-t1**3*t3*t2**8/18._ki&
                  &+t6*t1**3*t2**9/4._ki-5._ki/2._ki*t6**3*t1**4*t2**5-t6**2*t2**8*t3&
                  &**3/36._ki+t1**3*t2**9*t5/9._ki+11._ki/24._ki*t6**3*t1**3*t2**7+275&
                  &._ki/18._ki*t6**3*t1**3*t3**2*t2**3-110._ki/9._ki*t6**3*t1**3*t4**5&
                  &*t5**2-1045._ki/18._ki*t6**3*t1**3*t3**2*t4**3+275._ki/24._ki*t6**3&
                  &*t1**3*t4**2*t2**5+77._ki/6._ki*t6**3*t1**3*t4**4*t2**3+5._ki*t6**&
                  &2*t1**3*t3**3*t2**2+5._ki/2._ki*t6**2*t1**4*t4**2*t2**4+7._ki/4._ki&
                  &*t6**2*t1**4*t2**4*t3-5._ki/4._ki*t6**3*t1**2*t3**2*t2**5+t6**2*t&
                  &2**7*t3**3*t4/3._ki+20._ki/3._ki*t6**3*t1**4*t5*t2**4-25._ki/9._ki*t&
                  &6**2*t1**2*t3**3*t2**4-t6**3*t1*t3**2*t2**7/24._ki
                !
                stemp7=stemp9-70._ki/9._ki*t6**2*t4**4*t3**3*t2**4+68._ki/9._ki*t6**2&
                  &*t4**5*t3**3*t2**3-4._ki*t6**2*t4**6*t3**3*t2**2+8._ki/9._ki*t6**2&
                  &*t4**7*t3**3*t2-4._ki*t6**2*t1**4*t4*t2**5-39._ki/4._ki*t6**2*t1**&
                  &3*t4**2*t2**6-3._ki*t6**2*t1**3*t4**4*t2**4+9._ki*t6**2*t1**3*t4*&
                  &*3*t2**5-t6*t1**3*t4**3*t2**6/2._ki+5._ki/4._ki*t6*t1**3*t4**2*t2*&
                  &*7+2090._ki/9._ki*t6**3*t1**3*t4**3*t3*t5*t2+1375._ki/12._ki*t6**3*&
                  &t1**3*t3**2*t4**2*t2+1375._ki/36._ki*t6**3*t1**3*t4**2*t5**2*t2**&
                  &3+275._ki/3._ki*t6**3*t1**3*t4**2*t3*t2**3-836._ki/9._ki*t6**3*t1**&
                  &3*t4**3*t3*t2**2+308._ki/9._ki*t6**3*t1**3*t4**4*t3*t2-45._ki/2._ki&
                  &*t6**2*t1**3*t2**3*t3**2*t5-20._ki*t6**3*t1**4*t4**3*t5*t2-65._ki&
                  &/3._ki*t6**3*t1**4*t2*t3*t4**2+418._ki/9._ki*t6**3*t1**3*t4**3*t5*&
                  &t2**3+75._ki/2._ki*t6**3*t1**4*t4*t5**2*t2**2-308._ki/9._ki*t6**3*t&
                  &1**3*t4**4*t5*t2**2-220._ki/3._ki*t6**3*t1**3*t4*t3**2*t2**2
                !
                stemp9=stemp7-352._ki/9._ki*t6**3*t1**3*t4*t3*t2**4-770._ki/9._ki*t6*&
                  &*3*t1**3*t4**4*t3*t5-275._ki/9._ki*t6**3*t1**3*t4**2*t5*t2**4+5._k&
                  &i/24._ki*t6**2*t1*t2**7*t3**2*t5+20._ki/3._ki*t6**2*t1**2*t3*t5**2&
                  &*t2**5*t4-25._ki/6._ki*t6**3*t1**2*t4*t5*t2**5*t3+5._ki/3._ki*t6**3&
                  &*t1**2*t4*t3*t2**6+t6**3*t1*t4*t3**2*t2**6/2._ki+88._ki/9._ki*t6**&
                  &3*t1**3*t4**5*t5*t2+385._ki/9._ki*t6**3*t1**3*t4**4*t5**2*t2-1045&
                  &._ki/18._ki*t6**3*t1**3*t4**3*t5**2*t2**2-1375._ki/6._ki*t6**3*t1**&
                  &3*t4**2*t3*t5*t2**2-80._ki/3._ki*t6**3*t1**2*t4**5*t3*t5*t2+130._k&
                  &i/3._ki*t6**3*t1**2*t4**4*t3*t5*t2**2-110._ki/3._ki*t6**3*t1**2*t4&
                  &**3*t3*t5*t2**3-130._ki/3._ki*t6**3*t1**2*t4**4*t3**2*t2+880._ki/9&
                  &._ki*t6**3*t1**3*t4*t3*t5*t2**3+55._ki*t6**3*t1**2*t4**3*t3**2*t2&
                  &**2-133._ki/18._ki*t6**2*t1**2*t4**3*t3*t2**5+85._ki/12._ki*t6**3*t&
                  &1*t4**3*t3**2*t2**4+20._ki/3._ki*t6**2*t1**2*t4**5*t3*t5**2*t2+70&
                  &._ki/3._ki*t6**2*t1**2*t3**2*t4**4*t5*t2
                !
                stemp8=stemp9-14._ki/9._ki*t6**2*t1**2*t4**5*t3*t2**3-190._ki/3._ki*t&
                  &6**2*t1**2*t3**2*t4**3*t5*t2**2+125._ki/2._ki*t6**2*t1**2*t3**2*t&
                  &4**2*t5*t2**3+t1**3*t4*t3*t2**7/18._ki-2._ki/9._ki*t1**3*t4*t2**8*&
                  &t5+t1**3*t4**2*t5*t2**7/9._ki+7._ki/36._ki*t6**2*t1**2*t3*t2**8+40&
                  &._ki/3._ki*t6**3*t1**2*t4**5*t3**2+4._ki/3._ki*t6**3*t1*t4**7*t3**2&
                  &-5._ki/3._ki*t6**3*t1**5*t3*t5-t6**3*t1**5*t2**2*t4-10._ki/3._ki*t6&
                  &**3*t1**5*t5**2*t4+5._ki/2._ki*t6**3*t1**5*t5**2*t2-2._ki*t6**3*t1&
                  &**5*t5*t2**2+2._ki/3._ki*t6**3*t1**5*t2*t3+15._ki/2._ki*t6**3*t1**4&
                  &*t4**3*t2**2-65._ki/4._ki*t6**3*t1**4*t4**2*t2**3+45._ki/4._ki*t6**&
                  &3*t1**4*t4*t2**4+25._ki/2._ki*t6**3*t1**4*t3**2*t4+25._ki*t6**3*t1&
                  &**4*t4**3*t5**2-10._ki*t6**3*t1**4*t2**3*t3-25._ki/3._ki*t6**3*t1*&
                  &*4*t3**2*t2-209._ki/12._ki*t6**3*t1**3*t4**3*t2**4
                !
                stemp9=stemp8-11._ki/3._ki*t6**3*t1**3*t4**5*t2**2+32._ki/3._ki*t6**3&
                  &*t1**2*t4**5*t3*t2**2-3._ki/4._ki*t6**2*t1**3*t2**8-70._ki/3._ki*t6&
                  &**2*t1**2*t4**4*t3*t5**2*t2**2+40._ki/9._ki*t6**2*t1*t4**5*t3**3*&
                  &t2-130._ki/9._ki*t6**2*t1*t4**4*t3**3*t2**2+88._ki/9._ki*t6**3*t1**&
                  &3*t4*t5*t2**5-110._ki/9._ki*t6**3*t1**3*t4*t5**2*t2**4-205._ki/18.&
                  &_ki*t6**2*t1*t2**4*t3**3*t4**2+14._ki/3._ki*t6**2*t1**4*t4*t5*t2**&
                  &4+25._ki/6._ki*t6**2*t1**4*t4**2*t5**3*t2+15._ki/2._ki*t6**2*t1**3*&
                  &t4*t5**3*t2**4+15._ki*t6**2*t1**3*t3*t5**2*t2**4+125._ki/36._ki*t6&
                  &**2*t1*t2**5*t3**3*t4+205._ki/24._ki*t6**2*t1*t4**2*t3**2*t2**5*t&
                  &5+5._ki/12._ki*t6**3*t1**2*t2**6*t3*t5-5._ki/6._ki*t6**2*t1**2*t3*t&
                  &5**2*t2**6-91._ki/4._ki*t6**2*t1**3*t4**2*t3*t2**4+91._ki/8._ki*t6*&
                  &*2*t1**3*t4**2*t5*t2**5-21._ki/2._ki*t6**2*t1**3*t4**3*t5*t2**4+7&
                  &._ki/2._ki*t6**2*t1**3*t4**4*t5*t2**3-55._ki/3._ki*t6**2*t1*t4**3*t&
                  &3**2*t5*t2**4
                !
                stemp6=stemp9+15._ki*t6**2*t1**3*t4**3*t5**3*t2**2-5._ki*t6**2*t1**&
                  &3*t4**4*t5**3*t2-65._ki/4._ki*t6**2*t1**3*t4**2*t5**3*t2**3-135._k&
                  &i/2._ki*t6**2*t1**3*t4*t5**2*t3*t2**3+95._ki/9._ki*t6**2*t1**2*t3*&
                  &*3*t4**3*t2-41._ki/6._ki*t6**3*t1**2*t4**2*t3*t2**5-325._ki/6._ki*t&
                  &6**3*t1**4*t4**2*t5**2*t2+30._ki*t6**3*t1**4*t2**2*t3*t4+130._ki/&
                  &3._ki*t6**3*t1**4*t4**2*t5*t2**2+325._ki/6._ki*t6**3*t1**4*t4**2*t&
                  &3*t5-75._ki*t6**3*t1**4*t3*t5*t2*t4-61._ki/24._ki*t6**3*t1*t4**2*t&
                  &3**2*t2**5-275._ki/18._ki*t6**3*t1**3*t3*t2**4*t5-30._ki*t6**3*t1*&
                  &*4*t4*t5*t2**3+44._ki/3._ki*t6**3*t1**2*t4**3*t3*t2**4-205._ki/6._k&
                  &i*t6**3*t1**2*t4**2*t3**2*t2**3+20._ki/3._ki*t6**3*t1**2*t4**6*t3&
                  &*t5-8._ki/3._ki*t6**3*t1**2*t4**6*t3*t2-52._ki/3._ki*t6**3*t1**2*t4&
                  &**4*t3*t2**3+34._ki/3._ki*t6**3*t1*t4**5*t3**2*t2**2+8._ki/3._ki*t6&
                  &**3*t1**5*t5*t2*t4-6._ki*t6**3*t1*t4**6*t3**2*t2-35._ki/3._ki*t6**&
                  &3*t1*t4**4*t3**2*t2**3+25._ki*t6**3*t1**4*t3*t5*t2**2
                !
                stemp7=1._ki/t1**3/t2**10
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
              !
            case(4)
              !
              select case(par4_glob)
              !
              case(4)
                !
                stemp2=-(-6._ki*t6*t2**5*t4*t5+4._ki*t6*t2**4*t4**2*t5+6._ki*t3*t4*t&
                  &6*t2**4+18._ki*t6**2*t1*t2**3*t4-96._ki*t6**3*t1*t5*t4**2-36._ki*t&
                  &1*t6**3*t2**2*t5+36._ki*t6**3*t3*t2*t1-60._ki*t6**3*t4*t1*t3-144.&
                  &_ki*t6**3*t3*t2*t4**2+84._ki*t6**3*t3*t2**2*t4-12._ki*t6**2*t1*t2*&
                  &*4+72._ki*t6**3*t2**2*t5*t4**2-80._ki*t6**3*t4**3*t5*t2-28._ki*t6*&
                  &*3*t2**3*t5*t4+12._ki*t6**3*t1**2*t5+4._ki*t6**3*t5*t2**4+2._ki*t6&
                  &*t2**6*t5+32._ki*t6**3*t4**4*t5-t4*t5**2*t2**5+120._ki*t1*t6**3*t&
                  &4*t5*t2-2._ki*t6*t1*t2**4*t5-12._ki*t6**2*t2**3*t4**3+80._ki*t6**3&
                  &*t3*t4**3-16._ki*t6**3*t3*t2**3+3._ki*t6**2*t2**6-4._ki*t6*t2**5*t&
                  &3-15._ki*t6**2*t2**5*t4+24._ki*t6**2*t2**4*t4**2-2._ki*t3*t5*t2**5&
                  &+t5**2*t2**6)/t2**8*z_log(t1*t6/t2**2,1._ki)/12._ki
                !
                stemp4=(-16._ki*t3*t4**2*t5*t6**2*t1+4._ki*t3**2*t4**2*t5*t2*t6-12.&
                  &_ki*t4**3*t3**2*t6**2+24._ki*t3*t4*t5*t2*t6**2*t1+4._ki*t2*t1*t6**&
                  &2*t3**2+4._ki*t3*t4*t5**2*t2*t1*t6-t3**2*t4*t5**2*t2**2+6._ki*t2*&
                  &t6*t3**3*t4-6._ki*t4*t5*t3**2*t6*t2**2-15._ki*t2**2*t4*t3**2*t6**&
                  &2+24._ki*t2*t4**2*t3**2*t6**2-8._ki*t3*t2**2*t6**2*t5*t1-6._ki*t1*&
                  &t6**2*t3**2*t4-4._ki*t6*t3**3*t2**2+6._ki*t2*t5*t1*t6*t3**2-2._ki*&
                  &t2**2*t5*t3**3-6._ki*t4*t5**2*t6**2*t1**2-4._ki*t5*t1**2*t6**2*t3&
                  &+3._ki*t3**2*t2**3*t6**2+2._ki*t5*t3**2*t6*t2**3+6._ki*t2*t5**2*t1&
                  &**2*t6**2-4._ki*t3*t2**2*t5**2*t1*t6+t3**2*t5**2*t2**3)/t3**2/t2&
                  &**5*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki
                !
                stemp7=-44._ki/9._ki*t6**3*t1*t3*t4**4*t5-t6*t1**2*t2**4*t5*t3/9._ki&
                  &-7._ki/9._ki*t6*t1*t2**4*t3*t4**2*t5+7._ki/6._ki*t6*t1*t2**5*t3*t4*&
                  &t5-11._ki/18._ki*t6**3*t1*t3*t5*t2**4-5._ki/2._ki*t6**3*t1**2*t2*t3&
                  &**2+4._ki/9._ki*t6*t1*t2**5*t3**2-3._ki/8._ki*t6**2*t2**6*t1*t3+t6*&
                  &t1**2*t2**5*t5**2/6._ki-t6**3*t3*t5*t1**3/3._ki+25._ki/6._ki*t6**3*&
                  &t1**2*t4*t3**2-77._ki/6._ki*t6**3*t1*t2**2*t4*t3**2-t6*t1**2*t2**&
                  &4*t4*t5**2/6._ki+5._ki/2._ki*t6**3*t1**2*t3*t2**2*t5+20._ki/3._ki*t6&
                  &**3*t1**2*t3*t5*t4**2-2._ki/3._ki*t6*t1*t4*t3**2*t2**4-25._ki/3._ki&
                  &*t6**3*t1**2*t3*t4*t5*t2+110._ki/9._ki*t6**3*t1*t3*t4**3*t5*t2+77&
                  &._ki/18._ki*t6**3*t1*t3*t2**3*t5*t4
                !
                stemp6=stemp7+22._ki*t6**3*t1*t2*t4**2*t3**2-110._ki/9._ki*t6**3*t1*&
                  &t4**3*t3**2+3._ki/4._ki*t6**3*t4*t3**2*t2**4+22._ki/9._ki*t6**3*t1*&
                  &t3**2*t2**3+t6**2*t1**2*t3*t2**4/2._ki-8._ki/3._ki*t6**3*t3**2*t2*&
                  &*3*t4**2+14._ki/3._ki*t6**3*t3**2*t2**2*t4**3-t6**3*t3**2*t2**5/1&
                  &2._ki-11._ki*t6**3*t1*t3*t2**2*t5*t4**2+4._ki/3._ki*t6**3*t3**2*t4*&
                  &*5+3._ki/2._ki*t6**2*t1*t3*t2**3*t4**3-7._ki/18._ki*t6*t1*t2**6*t3*&
                  &t5+15._ki/8._ki*t6**2*t1*t3*t2**5*t4+7._ki/36._ki*t1*t2**5*t3*t4*t5&
                  &**2-3._ki/4._ki*t6**2*t1**2*t3*t2**3*t4-3._ki*t6**2*t1*t3*t2**4*t4&
                  &**2-4._ki*t6**3*t3**2*t4**4*t2+2._ki/9._ki*t1*t2**5*t3**2*t5-7._ki/&
                  &36._ki*t1*t2**6*t3*t5**2
                !
                stemp7=1._ki/t1/t2**8/t3
                !
                stemp5=stemp6*stemp7
                !
                stemp3=stemp4+stemp5
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
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
                stemp2=(-6._ki*t6**3*t1*t2+12._ki*t6**3*t1*t4+t6**3*t2**3-8._ki*t6**&
                  &3*t4**3+12._ki*t6**3*t4**2*t2-6._ki*t6**3*t2**2*t4-2._ki*t6**2*t1*&
                  &t5*t2+4._ki*t6**2*t5*t2*t4**2+4._ki*t6**2*t2*t3*t4+t6**2*t5*t2**3&
                  &-4._ki*t6**2*t5*t2**2*t4-2._ki*t6**2*t2**2*t3+t6*t2**3*t5**2-2._ki&
                  &*t6*t4*t5**2*t2**2-2._ki*t6*t2**2*t3*t5+t5**3*t2**3)/t2**6*z_log&
                  &(t1*t6/t2**2,1._ki)/4._ki
                !
                stemp3=-(-2._ki*t1*t6*t5+t3*t6*t2+t3*t2*t5-2._ki*t6*t3*t4)*(t3**2*t&
                  &6**2*t2**2-4._ki*t3**2*t6**2*t4*t2-2._ki*t3*t1*t6**2*t5*t2+2._ki*t&
                  &3**2*t1*t6**2+4._ki*t3**2*t6**2*t4**2+4._ki*t3*t6**2*t4*t1*t5+2._k&
                  &i*t6**2*t1**2*t5**2-2._ki*t2*t3**3*t6-2._ki*t6*t1*t3*t5**2*t2+t3*&
                  &*2*t5**2*t2**2)/t2**6/t3**3*q(4,(t2*t3-t1*t6)/t2/t3,sign_arg)/4&
                  &._ki+(60._ki*t6**3*t3**2*t1*t4-30._ki*t6**3*t1*t3**2*t2+11._ki*t6**&
                  &3*t2**3*t3**2-88._ki*t6**3*t3**2*t4**3-66._ki*t6**3*t2**2*t3**2*t&
                  &4+132._ki*t6**3*t3**2*t2*t4**2+8._ki*t6**2*t1**2*t2*t5**3-12._ki*t&
                  &6**2*t1*t5**2*t3*t2**2+24._ki*t6**2*t1*t2*t3*t4*t5**2+2._ki*t6**2&
                  &*t1*t2*t5*t3**2-68._ki*t6**2*t3**2*t4*t5*t2**2+44._ki*t6**2*t2*t3&
                  &**3*t4+68._ki*t6**2*t2*t3**2*t4**2*t5-22._ki*t6**2*t2**2*t3**3+17&
                  &._ki*t6**2*t2**3*t3**2*t5-16._ki*t6*t1*t3*t2**2*t5**3-40._ki*t6*t3&
                  &**2*t2**2*t4*t5**2+20._ki*t6*t3**2*t2**3*t5**2-28._ki*t6*t3**3*t2&
                  &**2*t5+22._ki*t3**2*t2**3*t5**3)/t3**2/t2**6/24._ki
                !
                stemp1=stemp2+stemp3
                !
                stemp2=1._ki/t2
                !
                temp0=stemp1*stemp2
                !
              end select
              !
            end select
            !
          end select
          !
        end if
        !
      else if (dim_glob == "n+4") then
        if (nb_par == 0) then
          !
          temp0=(-1._ki/6._ki*(-4._ki*t6**2*t1*t3+4._ki*t6**2*t1*t5*t2-8._ki*t5*&
            &t4*t6**2*t1-2._ki*t6**2*t1*t2**2+4._ki*t6**2*t1*t2*t4+2._ki*t6**2*&
            &t3*t2**2-8._ki*t6**2*t3*t2*t4+8._ki*t6**2*t3*t4**2+2._ki*t2**3*t1*&
            &t6+t2**4*t3)/t2**5*z_log(-t1*t6/t2**2,-1._ki)+1._ki/6._ki*t3/t2*q(&
            &2,(t2*t3-t1*t6)/t2/t3,sign_arg)-1._ki/6._ki*(-2._ki*t1*t3+2._ki*t1*&
            &t2*t4+2._ki*t1*t5*t2-4._ki*t1*t5*t4-t1*t2**2+12._ki*t3*t4**2-12._ki&
            &*t3*t2*t4+3._ki*t3*t2**2)*t6**2/t2**5)/t2
          !
        else if (nb_par == 1) then
          !
          select case(par4_glob)
          !
          case(1)
            !
            stemp2=(9._ki*t1**2*t5*t6**3-3._ki*t1**2*t6**3*t2+3._ki*t1*t6**3*t2*&
              &*3-36._ki*t1*t5*t6**3*t4**2+12._ki*t1*t6**3*t2*t4**2-12._ki*t1*t6*&
              &*3*t2**2*t4+36._ki*t1*t5*t6**3*t2*t4+18._ki*t1*t3*t6**3*t2-36._ki*&
              &t1*t3*t6**3*t4-9._ki*t1*t5*t6**3*t2**2-3._ki*t1*t2**4*t6**2+6._ki*&
              &t1*t2**3*t6**2*t4-t1*t5*t2**4*t6-t1*t2**5*t6-3._ki*t3*t6**3*t2**&
              &3+24._ki*t3*t6**3*t4**3+18._ki*t3*t6**3*t2**2*t4-36._ki*t3*t6**3*t&
              &2*t4**2-t3*t2**5*t6+2._ki*t3*t2**4*t6*t4-t3*t2**6-t3*t5*t2**5)/t&
              &2**7*z_log(-t1*t6/t2**2,-1._ki)/12._ki
            !
            stemp3=-(2._ki*t5*t1*t6-t3*t5*t2+2._ki*t2*t1*t6-t3*t2*t6+2._ki*t3*t6&
              &*t4-t3*t2**2)/t2**3*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki+(-2&
              &._ki*t1**2*t6**3*t2+6._ki*t1**2*t5*t6**3-15._ki*t1*t5*t6**3*t2**2-&
              &20._ki*t1*t6**3*t2**2*t4+60._ki*t1*t5*t6**3*t2*t4-60._ki*t1*t5*t6*&
              &*3*t4**2+30._ki*t1*t3*t6**3*t2-60._ki*t1*t3*t6**3*t4+20._ki*t1*t6*&
              &*3*t2*t4**2+5._ki*t1*t6**3*t2**3-3._ki*t1*t2**4*t6**2+6._ki*t1*t2*&
              &*3*t6**2*t4-132._ki*t3*t6**3*t2*t4**2+66._ki*t3*t6**3*t2**2*t4+88&
              &._ki*t3*t6**3*t4**3-11._ki*t3*t6**3*t2**3-2._ki*t3*t2**5*t6+4._ki*t&
              &3*t2**4*t6*t4-2._ki*t3*t2**6-2._ki*t3*t5*t2**5)/t2**7/24._ki
            !
            stemp1=stemp2+stemp3
            !
            stemp2=1._ki/t2
            !
            temp0=stemp1*stemp2
            !
          case(2)
            !
            stemp4=-6._ki*t6**2*t4**2*t5**2*t3*t2**2-9._ki*t6**2*t4*t3**2*t2**2&
              &*t5-9._ki/2._ki*t1**2*t6**3*t4*t5*t2+6._ki*t6**3*t4**2*t2**2*t3*t5&
              &+12._ki*t1*t3*t6**3*t2*t4**2-5._ki/6._ki*t1*t6**2*t4*t2**4*t5+3._ki&
              &/2._ki*t1*t6**3*t4*t5*t2**3+12._ki*t6**2*t4**2*t3**2*t5*t2+6._ki*t&
              &6**2*t4**3*t5**2*t3*t2+3._ki/2._ki*t6**2*t4*t5**2*t3*t2**3-3._ki/2&
              &._ki*t1*t6**3*t5**2*t2**2*t4-24._ki*t1*t6**3*t3*t5*t4**2-2._ki*t1*&
              &t6**2*t4**2*t5**3*t2+10._ki*t6**3*t4**3*t3**2+t1**2*t6**2*t2**4/&
              &2._ki-t6**2*t3**3*t2**2-3._ki/4._ki*t1**2*t2**3*t6**3+3._ki*t1*t6**&
              &2*t5**2*t2**2*t3-3._ki*t1*t6**2*t3**2*t5*t2+6._ki*t1*t6**3*t4**2*&
              &t5**2*t2-3._ki*t1*t6**3*t3*t5*t2**2+6._ki*t1*t6**3*t4**3*t5*t2-6.&
              &_ki*t1*t6**3*t4**2*t5*t2**2+18._ki*t1*t6**3*t3*t5*t2*t4+9._ki/2._ki&
              &*t1**2*t6**3*t5**2*t4+9._ki/4._ki*t1**2*t6**3*t2**2*t4+3._ki*t1**2&
              &*t6**3*t5*t3-3._ki/2._ki*t1**2*t2*t6**3*t5**2-5._ki/12._ki*t1**2*t2&
              &**3*t6**2*t5+t1**2*t6**2*t5**3*t2/2._ki-3._ki/4._ki*t1*t6**3*t4*t2&
              &**4-5._ki/6._ki*t1*t6**2*t3*t2**4
            !
            stemp3=stemp4+3._ki*t1*t6**3*t3**2*t2-6._ki*t1*t6**3*t4**3*t5**2-15&
              &._ki/2._ki*t1*t6**3*t4*t3**2-3._ki*t1*t6**3*t4**3*t2**2+9._ki/2._ki*&
              &t6**3*t2**2*t3**2*t4+3._ki/2._ki*t1**2*t5*t6**3*t2**2-3._ki/2._ki*t&
              &1**2*t6**3*t2*t3+t1*t6**2*t2**5*t4-2._ki*t1*t6**2*t4**2*t2**4+8.&
              &_ki*t6**3*t4**4*t5*t3+t6**3*t4*t2**4*t3/2._ki-12._ki*t6**3*t4**2*t&
              &3**2*t2-3._ki*t6**3*t4**2*t2**3*t3+6._ki*t6**3*t4**3*t3*t2**2-4._k&
              &i*t6**3*t4**4*t3*t2+3._ki/2._ki*t6**2*t2**3*t3**2*t5+5._ki/2._ki*t6&
              &**2*t4*t3**3*t2-5._ki/3._ki*t6**2*t4**3*t2**3*t3+5._ki/3._ki*t6**2*&
              &t4**2*t2**4*t3-5._ki/12._ki*t6**2*t4*t2**5*t3+3._ki/2._ki*t1*t6**3*&
              &t3*t2**3+3._ki*t1*t6**3*t4**2*t2**3-t1*t6*t4*t2**6/4._ki-9._ki*t1*&
              &t6**2*t5**2*t2*t3*t4+5._ki/3._ki*t1*t6**2*t4**2*t2**3*t5+5._ki/2._k&
              &i*t1*t6**2*t4*t3*t2**3-12._ki*t6**3*t4**3*t5*t3*t2-t6**3*t4*t2**&
              &3*t3*t5-9._ki*t4*t1*t6**3*t3*t2**2+t1*t6**2*t4*t5**3*t2**2-t4*t3&
              &*t2**7/12._ki-t3**2*t6**3*t2**3/2._ki
            !
            stemp4=1._ki/t2**9*z_log(-t1*t6/t2**2,-1._ki)
            !
            stemp2=stemp3*stemp4
            !
            stemp4=t4*t3/t2**2*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki
            !
            stemp8=11._ki*t1**2*t6*t4**2*t2**2*t3*t5-t1**3*t4**2*t5**3*t2+5._ki&
              &*t1**3*t6*t4**3*t5*t2-5._ki*t1**3*t6*t4**2*t5*t2**2+5._ki*t1**3*t&
              &6*t4**2*t5**2*t2-20._ki*t1**3*t6*t3*t5*t4**2-5._ki/2._ki*t1**3*t6*&
              &t3*t5*t2**2-11._ki/2._ki*t1**2*t6*t4**2*t2**3*t3-22._ki/3._ki*t1**2&
              &*t6*t4**4*t3*t2+4._ki*t1*t6*t4**4*t3**2*t2-3._ki/2._ki*t1**2*t3**3&
              &*t2**2-2._ki/3._ki*t4**5*t3**3*t2-t4**3*t3**3*t2**3+4._ki/3._ki*t4*&
              &*4*t3**3*t2**2-t4*t3**3*t2**5/24._ki+10._ki*t1**3*t6*t4**2*t3*t2-&
              &2._ki*t1*t6*t4**5*t3**2+3._ki/2._ki*t1**4*t6*t5**2*t4-t1**4*t3*t2*&
              &t6/2._ki
            !
            stemp7=stemp8+t1**4*t5*t6*t2**2/2._ki-3._ki/2._ki*t1*t2**3*t3**3*t4-&
              &5._ki/12._ki*t1**3*t4*t5*t2**4+5._ki/6._ki*t1**3*t4**2*t5*t2**3+33.&
              &_ki/4._ki*t1**2*t6*t4*t3**2*t2**2-22._ki*t1**2*t6*t4**3*t5*t3*t2-2&
              &2._ki*t1**2*t6*t4**2*t3**2*t2-11._ki/6._ki*t1**2*t6*t4*t2**3*t3*t5&
              &+11._ki/12._ki*t1**2*t3*t2**4*t6*t4+44._ki/3._ki*t1**2*t6*t4**4*t5*&
              &t3-3._ki/2._ki*t1**4*t6*t4*t5*t2+15._ki*t1**3*t6*t3*t5*t2*t4+11._ki&
              &*t1**2*t6*t4**3*t3*t2**2-t1*t3**2*t2**4*t6*t4/8._ki-9._ki/2._ki*t1&
              &**3*t5**2*t2*t3*t4-9._ki*t1**2*t4**2*t5**2*t3*t2**2+9._ki*t1**2*t&
              &4**3*t5**2*t3*t2+18._ki*t1**2*t4**2*t3**2*t5*t2+9._ki/4._ki*t1**2*&
              &t4*t5**2*t3*t2**3
            !
            stemp8=stemp7-27._ki/2._ki*t1**2*t4*t3**2*t2**2*t5-3._ki*t1*t4**2*t5&
              &*t3**2*t2**3-4._ki*t1*t4**4*t5*t3**2*t2+6._ki*t1*t4**3*t5*t3**2*t&
              &2**2+t1*t4*t5*t3**2*t2**4/2._ki-5._ki/4._ki*t1**3*t6*t4*t5**2*t2**&
              &2+5._ki/4._ki*t1**3*t6*t4*t5*t2**3-t1**4*t2**3*t6/4._ki+t4**2*t3**&
              &3*t2**4/3._ki+t1**3*t4*t2**5/2._ki-t1**3*t4**2*t2**4-5._ki/12._ki*t&
              &1**3*t2**4*t3-11._ki/12._ki*t1**2*t2**3*t3**2*t6+t1**4*t6*t5*t3-t&
              &1**4*t2*t6*t5**2/2._ki+3._ki/4._ki*t1**4*t4*t6*t2**2-5._ki/2._ki*t1*&
              &*3*t6*t4**3*t2**2-25._ki/4._ki*t1**3*t6*t4*t3**2+5._ki/2._ki*t1**3*&
              &t6*t4**2*t2**3
            !
            stemp6=stemp8-5._ki/8._ki*t1**3*t6*t4*t2**4+5._ki/4._ki*t1**3*t6*t3*t&
              &2**3+5._ki/2._ki*t1**3*t6*t3**2*t2-5._ki*t1**3*t6*t4**3*t5**2+55._k&
              &i/3._ki*t1**2*t6*t4**3*t3**2+5._ki/4._ki*t1**3*t4*t3*t2**3+3._ki/2.&
              &_ki*t1**3*t5**2*t2**2*t3-3._ki/2._ki*t1**3*t3**2*t5*t2+t1**3*t4*t5&
              &**3*t2**2/2._ki-5._ki/8._ki*t1**2*t4*t2**5*t3+5._ki/2._ki*t1**2*t4**&
              &2*t2**4*t3-5._ki/2._ki*t1**2*t4**3*t2**3*t3+9._ki/4._ki*t1**2*t2**3&
              &*t3**2*t5+15._ki/4._ki*t1**2*t4*t3**3*t2-10._ki/3._ki*t1*t4**3*t3**&
              &3*t2+4._ki*t1*t4**2*t3**3*t2**2-3._ki*t1*t6*t4**3*t2**2*t3**2+t1*&
              &t3**2*t2**3*t6*t4**2-15._ki/2._ki*t1**3*t6*t4*t3*t2**2+t1*t2**4*t&
              &3**3/6._ki
            !
            stemp7=t6**2/t2**9/t1**2
            !
            stemp5=stemp6*stemp7
            !
            stemp3=stemp4+stemp5
            !
            stemp1=stemp2+stemp3
            !
            stemp2=1._ki/t2
            !
            temp0=stemp1*stemp2
            !
          case(3)
            !
            stemp5=5._ki/2._ki*t1*t6**2*t4*t2**4*t5+15._ki/2._ki*t1*t6**3*t5**2*t&
              &2**2*t4+t4*t3*t2**7/12._ki-15._ki/2._ki*t1*t6**3*t4*t5*t2**3-18._ki&
              &*t6**3*t4**2*t2**2*t3*t5+24._ki*t1*t6**3*t3*t5*t4**2-3._ki*t1*t6*&
              &*2*t4*t5**3*t2**2-5._ki/6._ki*t1*t6**2*t2**5*t5+3._ki/2._ki*t1*t6**&
              &3*t5*t2**4-3._ki/2._ki*t1*t6**3*t5**2*t2**3+3._ki/2._ki*t6**2*t5**2&
              &*t2**4*t3+25._ki/12._ki*t6**2*t4*t2**5*t3+5._ki/3._ki*t6**2*t4**3*t&
              &2**3*t3-10._ki/3._ki*t6**2*t4**2*t2**4*t3+3._ki/2._ki*t1**2*t2**3*t&
              &6**3+t6**3*t3*t2**5/2._ki-5._ki/12._ki*t6**2*t3*t2**6-10._ki*t6**3*&
              &t4**3*t3**2+12._ki*t6**2*t4**2*t5**2*t3*t2**2
            !
            stemp4=stemp5+2._ki*t1*t6**2*t4**2*t5**3*t2-6._ki*t1*t6**2*t5**2*t2&
              &**2*t3+3._ki*t1*t6**2*t3**2*t5*t2+3._ki/2._ki*t6**2*t3**3*t2**2-t1&
              &**2*t6**2*t2**4/2._ki-3._ki/4._ki*t1*t6**3*t2**5+2._ki*t3**2*t6**3*&
              &t2**3-t1*t6*t2**7/4._ki-12._ki*t1*t6**3*t4**2*t5**2*t2+9._ki*t1*t6&
              &**3*t3*t5*t2**2-12._ki*t6**2*t4**2*t3**2*t5*t2-6._ki*t6**2*t4**3*&
              &t5**2*t3*t2+9._ki/2._ki*t1**2*t6**3*t4*t5*t2-6._ki*t1*t6**3*t4**3*&
              &t5*t2-9._ki/2._ki*t1**2*t6**3*t5**2*t4-9._ki/4._ki*t1**2*t6**3*t2**&
              &2*t4-3._ki*t1**2*t6**3*t5*t3+3._ki*t1**2*t2*t6**3*t5**2+5._ki/12._k&
              &i*t1**2*t2**3*t6**2*t5
            !
            stemp5=stemp4-t1**2*t6**2*t5**3*t2/2._ki-15._ki/2._ki*t6**2*t4*t5**2&
              &*t3*t2**3+15._ki*t4*t1*t6**3*t3*t2**2+15._ki/4._ki*t1*t6**3*t4*t2*&
              &*4+5._ki/3._ki*t1*t6**2*t3*t2**4-9._ki/2._ki*t1*t6**3*t3**2*t2+6._ki&
              &*t1*t6**3*t4**3*t5**2+15._ki/2._ki*t1*t6**3*t4*t3**2-6._ki*t1*t6**&
              &3*t4**2*t2**3+20._ki*t6**3*t4**3*t5*t3*t2-30._ki*t1*t6**3*t3*t5*t&
              &2*t4+12._ki*t1*t6**3*t4**2*t5*t2**2+9._ki*t1*t6**2*t5**2*t2*t3*t4&
              &+15._ki*t6**2*t4*t3**2*t2**2*t5+3._ki*t1*t6**3*t4**3*t2**2-21._ki/&
              &2._ki*t6**3*t2**2*t3**2*t4-3._ki*t1**2*t5*t6**3*t2**2+3._ki/2._ki*t&
              &1**2*t6**3*t2*t3
            !
            stemp3=stemp5-3._ki*t1*t6**2*t2**5*t4+2._ki*t1*t6**2*t4**2*t2**4-t6&
              &**3*t5*t3*t2**4+t1*t6*t4*t2**6/4._ki+7._ki*t6**3*t4*t2**3*t3*t5-t&
              &3*t2**8/12._ki-5._ki/3._ki*t1*t6**2*t4**2*t2**3*t5+t1*t6**2*t5**3*&
              &t2**3-9._ki/2._ki*t1*t6**3*t3*t2**3-8._ki*t6**3*t4**4*t5*t3-7._ki/2&
              &._ki*t6**3*t4*t2**4*t3+18._ki*t6**3*t4**2*t3**2*t2+9._ki*t6**3*t4*&
              &*2*t2**3*t3-10._ki*t6**3*t4**3*t3*t2**2+4._ki*t6**3*t4**4*t3*t2-9&
              &._ki/2._ki*t6**2*t2**3*t3**2*t5-5._ki/2._ki*t6**2*t4*t3**3*t2+t1*t6&
              &**2*t2**6-12._ki*t1*t3*t6**3*t2*t4**2-5._ki/2._ki*t1*t6**2*t4*t3*t&
              &2**3
            !
            stemp4=1._ki/t2**9*z_log(-t1*t6/t2**2,-1._ki)
            !
            stemp2=stemp3*stemp4
            !
            stemp4=-(-t2+t4)*t3/t2**2*q(3,(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki
            !
            !
            stemp8=33._ki/2._ki*t1**2*t6*t4**2*t2**3*t3-44._ki/3._ki*t1**2*t6*t4*&
              &*4*t5*t3-2._ki*t4**4*t3**3*t2**2+9._ki/2._ki*t1**3*t5**2*t2*t3*t4+&
              &18._ki*t1**2*t4**2*t5**2*t3*t2**2-9._ki*t1**2*t4**3*t5**2*t3*t2+5&
              &._ki/4._ki*t1**3*t6*t5*t2**4-2._ki/3._ki*t1*t2**4*t3**3-5._ki/8._ki*t&
              &1**2*t3*t2**6+3._ki/8._ki*t4*t3**3*t2**5+9._ki/8._ki*t1*t3**2*t2**4&
              &*t6*t4+9._ki/4._ki*t1**2*t3**3*t2**2+7._ki/3._ki*t4**3*t3**3*t2**3+&
              &2._ki/3._ki*t4**5*t3**3*t2-18._ki*t1**2*t4**2*t3**2*t5*t2-45._ki/4.&
              &_ki*t1**2*t4*t5**2*t3*t2**3+11._ki/12._ki*t1**2*t3*t2**5*t6-t1*t6*&
              &t3**2*t2**5/8._ki-t1**4*t6*t5*t3+t1**4*t2*t6*t5**2-3._ki/4._ki*t1*&
              &*4*t4*t6*t2**2+5._ki/2._ki*t1**3*t6*t4**3*t2**2
            !
            stemp7=stemp8+11._ki/3._ki*t1**2*t2**3*t3**2*t6-4._ki*t1*t3**2*t2**3&
              &*t6*t4**2-11._ki/6._ki*t1**2*t6*t5*t3*t2**4-55._ki/3._ki*t1**2*t6*t&
              &4**3*t3*t2**2-33._ki*t1**2*t6*t4**2*t2**2*t3*t5+25._ki/4._ki*t1**3&
              &*t6*t4*t5**2*t2**2+25._ki/2._ki*t1**3*t6*t4*t3*t2**2-5._ki*t1**3*t&
              &6*t4**3*t5*t2-25._ki/4._ki*t1**3*t6*t4*t5*t2**3+9._ki*t1*t4**2*t5*&
              &t3**2*t2**3+t1**3*t4**2*t5**3*t2-3._ki/2._ki*t1**3*t4*t5**3*t2**2&
              &+25._ki/8._ki*t1**2*t4*t2**5*t3-5._ki*t1**2*t4**2*t2**4*t3+5._ki/2.&
              &_ki*t1**2*t4**3*t2**3*t3-27._ki/4._ki*t1**2*t2**3*t3**2*t5-15._ki/4&
              &._ki*t1**2*t4*t3**3*t2+10._ki/3._ki*t1*t4**3*t3**3*t2-6._ki*t1*t4**&
              &2*t3**3*t2**2+10._ki*t1**3*t6*t4**2*t5*t2**2+45._ki/2._ki*t1**2*t4&
              &*t3**2*t2**2*t5+t1**3*t5**3*t2**3/2._ki+4._ki*t1*t4**4*t5*t3**2*t&
              &2
            !
            stemp8=stemp7+3._ki/2._ki*t1**4*t6*t4*t5*t2-10._ki*t1**3*t6*t4**2*t3&
              &*t2-10._ki*t1**3*t6*t4**2*t5**2*t2-3._ki/2._ki*t1**3*t4*t2**5+t1**&
              &3*t4**2*t2**4-5._ki/8._ki*t1**3*t6*t2**5+t1**4*t2**3*t6/2._ki+20._k&
              &i*t1**3*t6*t3*t5*t4**2+5._ki/6._ki*t1**3*t2**4*t3+22._ki/3._ki*t1**&
              &2*t6*t4**4*t3*t2+7._ki*t1*t6*t4**3*t2**2*t3**2+15._ki/2._ki*t1**3*&
              &t6*t3*t5*t2**2-77._ki/4._ki*t1**2*t6*t4*t3**2*t2**2-25._ki*t1**3*t&
              &6*t3*t5*t2*t4-5._ki/12._ki*t1**3*t2**5*t5-5._ki/4._ki*t1**3*t6*t5**&
              &2*t2**3+25._ki/4._ki*t1**3*t6*t4*t3**2+25._ki/8._ki*t1**3*t6*t4*t2*&
              &*4-15._ki/4._ki*t1**3*t6*t3*t2**3+110._ki/3._ki*t1**2*t6*t4**3*t5*t&
              &3*t2+33._ki*t1**2*t6*t4**2*t3**2*t2-6._ki*t1*t6*t4**4*t3**2*t2
            !
            stemp6=stemp8-10._ki*t1*t4**3*t5*t3**2*t2**2-7._ki/2._ki*t1*t4*t5*t3&
              &**2*t2**4+77._ki/6._ki*t1**2*t6*t4*t2**3*t3*t5-3._ki/2._ki*t1**4*t6&
              &*t5**2*t4-15._ki/4._ki*t1**3*t6*t3**2*t2+5._ki*t1**3*t6*t4**3*t5**&
              &2-55._ki/3._ki*t1**2*t6*t4**3*t3**2+2._ki*t1*t6*t4**5*t3**2-5._ki*t&
              &1**3*t6*t4**2*t2**3+9._ki/4._ki*t1**2*t5**2*t2**4*t3+t1*t3**2*t5*&
              &t2**5/2._ki-t3**3*t2**6/24._ki+7._ki/2._ki*t1*t2**3*t3**3*t4+5._ki/4&
              &._ki*t1**3*t4*t5*t2**4-5._ki/6._ki*t1**3*t4**2*t5*t2**3+t1**4*t3*t&
              &2*t6/2._ki-t1**4*t5*t6*t2**2-5._ki/4._ki*t1**3*t4*t3*t2**3-3._ki*t1&
              &**3*t5**2*t2**2*t3+3._ki/2._ki*t1**3*t3**2*t5*t2-4._ki/3._ki*t4**2*&
              &t3**3*t2**4-77._ki/12._ki*t1**2*t3*t2**4*t6*t4+t1**3*t2**6/2._ki
            !
            stemp7=t6**2/t2**9/t1**2
            !
            stemp5=stemp6*stemp7
            !
            stemp3=stemp4+stemp5
            !
            stemp1=stemp2+stemp3
            !
            stemp2=1._ki/t2
            !
            temp0=stemp1*stemp2
            !
          case(4)
            !
            stemp2=-(-3._ki*t1**2*t6**3*t2+9._ki*t1**2*t5*t6**3+36._ki*t1*t5*t6*&
              &*3*t2*t4+18._ki*t1*t3*t6**3*t2+3._ki*t1*t6**3*t2**3-9._ki*t1*t5*t6&
              &**3*t2**2-36._ki*t1*t3*t6**3*t4+12._ki*t1*t6**3*t2*t4**2-36._ki*t1&
              &*t5*t6**3*t4**2-12._ki*t1*t6**3*t2**2*t4-6._ki*t1*t2**4*t6**2+12.&
              &_ki*t1*t2**3*t6**2*t4-2._ki*t1*t5*t2**4*t6+24._ki*t3*t6**3*t4**3-3&
              &._ki*t3*t6**3*t2**3+18._ki*t3*t6**3*t2**2*t4-36._ki*t3*t6**3*t2*t4&
              &**2+4._ki*t3*t2**4*t6*t4-2._ki*t3*t2**5*t6-2._ki*t3*t5*t2**5)/t2**&
              &7*z_log(-t1*t6/t2**2,-1._ki)/24._ki
            !
            stemp3=(2._ki*t5*t1*t6-t3*t2*t6-t3*t5*t2+2._ki*t3*t6*t4)/t2**3*q(3,&
              &(t2*t3-t1*t6)/t2/t3,sign_arg)/12._ki-(6._ki*t1**2*t5*t6**3-2._ki*t&
              &1**2*t6**3*t2+30._ki*t1*t3*t6**3*t2-15._ki*t1*t5*t6**3*t2**2-60._k&
              &i*t1*t5*t6**3*t4**2+20._ki*t1*t6**3*t2*t4**2-20._ki*t1*t6**3*t2**&
              &2*t4+60._ki*t1*t5*t6**3*t2*t4-60._ki*t1*t3*t6**3*t4+5._ki*t1*t6**3&
              &*t2**3-6._ki*t1*t2**4*t6**2+12._ki*t1*t2**3*t6**2*t4+88._ki*t3*t6*&
              &*3*t4**3-132._ki*t3*t6**3*t2*t4**2+66._ki*t3*t6**3*t2**2*t4-11._ki&
              &*t3*t6**3*t2**3+8._ki*t3*t2**4*t6*t4-4._ki*t3*t2**5*t6-4._ki*t3*t5&
              &*t2**5)/t2**7/48._ki
            !
            stemp1=stemp2+stemp3
            !
            stemp2=1._ki/t2
            !
            temp0=stemp1*stemp2
            !
          end select
          !
        end if
        !
      end if
      !
      compute_residue = temp0
      !
    end function compute_residue
    !
    !****if* src/integrals/four_point/function_4p3m/fg
    ! NAME
    !
    !  function fg
    !
    ! USAGE
    !
    !  complex = fg(z,s24,s13,s12,s23,s34,par1,par2,par3,par4,flag,dim)
    !
    ! DESCRIPTION
    !
    !  This function contains the one dimensional integral representation of 
    !  the six/eight dimensional three mass four point function
    !
    !
    ! INPUTS
    !
    !  * z -- a real (type ki), integration variable
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s12 -- a real (type ki), the S matrix element 1,2
    !  * s23 -- a real (type ki), the S matrix element 2,3
    !  * s34 -- a real (type ki), the S matrix element 3,4
    !  * par1 -- an integer, the label of the fourth Feynman parameter, if none, put 0
    !  * par2 -- an integer, the label of the third Feynman parameter, if none, put 0
    !  * par3 -- an integer, the label of the second Feynman parameter, if none, put 0
    !  * par4 -- an integer, the label of the first Feynman parameter, if none, put 0
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           three mass four point function, dim="n+4" eight dimensional
    !           three mass four point function
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
    !  three mass four point function
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function fg(z,s24,s13,s12,s23,s34,par1,par2,par3,par4,flag,dim)
      !
      complex(ki), intent (in) :: z
      real(ki), intent (in) :: s24,s13,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4,flag
      character (len=3) :: dim
      complex(ki) :: fg
      !
      integer, dimension(4) :: par
      integer :: nb_par
      complex(ki) :: c_var,d_var,e_var,f_var,g_var,h_var
      !
      par = (/par1,par2,par3,par4/)
      nb_par = count(mask=par/=0)
      !
      c_var = z*s12+(1._ki-z)*s13
      !
      f_var = z*(s24-s12)+(1._ki-z)*(s34-s13)
      !
      g_var = z*(1._ki-z)*s23-z*s24-(1._ki-z)*s34
      !
      d_var = z*(1._ki-z)*s23-z*s12-(1._ki-z)*s13
      !
      e_var = z*s24+(1._ki-z)*s34
      !
      h_var = z*(1._ki-z)*s23
      !
      if (dim == "n+2") then      
        if (nb_par == 0) then
          select case(flag)
          !
          case(1)
            !
            fg=(-(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var+log(e_var)&
              &*e_var)/f_var/g_var
            !
          case(2)
            !
            fg=((log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var-c_var*log(c&
              &_var))/f_var/d_var
            !
          end select
        else if (nb_par == 1) then
          select case(par4)
          !
          case(1)
            !
            select case(flag)
            !
            case(1)
              !
              fg=1._ki/2._ki*(-(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var*&
                &*2+log(e_var)*e_var**2)/f_var**2/g_var
              !
            case(2)
              !
              fg=-(1._ki/2._ki*c_var**2/f_var*log(c_var)-1._ki/2._ki*(log(z)&
                &+log(1._ki-z)+z_log(s23,1._ki))*h_var**2/f_var)/d_var**2-(&
                &1._ki/2._ki*c_var/f_var**2*log(c_var)*(2._ki*f_var+c_var)+1&
                &._ki/2._ki*(c_var*f_var-(log(z)+log(1._ki-z)+z_log(s23,1._ki&
                &))*h_var**2)/f_var**2)/d_var
              !
            end select
            !
          case(2)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-(-1._ki/2._ki*z*(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_v&
                &ar**2/f_var**2*d_var+1._ki/2._ki*z*log(e_var)/f_var*e_var*&
                &*2)/g_var**2-(-1._ki/2._ki*z*(log(z)+log(1._ki-z)+z_log(s23&
                &,1._ki))*h_var/f_var**2*(c_var-d_var)+z*(log(z)+log(1._ki-&
                &z)+z_log(s23,1._ki))*h_var**2/f_var**3*d_var+1._ki/2._ki*z*&
                &(c_var+f_var)/f_var)/g_var-z*(log(z)+log(1._ki-z)+z_log(s&
                &23,1._ki))*h_var*(-d_var-h_var+g_var+c_var)/f_var**3
              !
            case(2)
              !
              fg=-(-1._ki/2._ki*z*(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_v&
                &ar/f_var**2*c_var*g_var-1._ki/2._ki*z*c_var**2/f_var*log(c&
                &_var))/d_var**2-(-1._ki/2._ki*z*(log(z)+log(1._ki-z)+z_log(&
                &s23,1._ki))*h_var/f_var**2*(h_var-g_var)-z*(log(z)+log(1.&
                &_ki-z)+z_log(s23,1._ki))*h_var/f_var**3*c_var*g_var-1._ki/2&
                &._ki*z*c_var/f_var)/d_var
              !
            end select
            !
          case(3)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-(1._ki/2._ki*(-1._ki+z)*(log(z)+log(1._ki-z)+z_log(s23,1._k&
                &i))*h_var**2/f_var**2*d_var-1._ki/2._ki*(-1._ki+z)*log(e_va&
                &r)/f_var*e_var**2)/g_var**2-(1._ki/2._ki*(-1._ki+z)*(log(z)&
                &+log(1._ki-z)+z_log(s23,1._ki))*h_var/f_var**2*(c_var-d_va&
                &r)-(-1._ki+z)*(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var*&
                &*2/f_var**3*d_var-1._ki/2._ki*(-c_var+z*c_var-f_var+z*f_va&
                &r)/f_var)/g_var+(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_v&
                &ar*(-d_var-h_var+g_var+c_var)*(-1._ki+z)/f_var**3
              !
            case(2)
              !
              fg=-(1._ki/2._ki*(-1._ki+z)*(log(z)+log(1._ki-z)+z_log(s23,1._k&
                &i))*h_var/f_var**2*c_var*g_var+1._ki/2._ki*c_var**2*(-1._ki&
                &+z)/f_var*log(c_var))/d_var**2-(1._ki/2._ki*(-1._ki+z)*(log&
                &(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var/f_var**2*(h_var-g&
                &_var)+(-1._ki+z)*(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_v&
                &ar/f_var**3*c_var*g_var+1._ki/2._ki*(-1._ki+z)*c_var/f_var)&
                &/d_var
              !
            end select
            !
          case(4)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-(-1._ki/2._ki*log(e_var)/f_var*e_var**2+1._ki/2._ki*(log(z&
                &)+log(1._ki-z)+z_log(s23,1._ki))*h_var**2/f_var)/g_var**2-&
                &(-1._ki/2._ki*log(e_var)/f_var**2*e_var*(-c_var+f_var)-1._k&
                &i/2._ki*h_var*((log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var&
                &+f_var)/f_var**2)/g_var-1._ki/2._ki/f_var
              !
            case(2)
              !
              fg=1._ki/2._ki*(-(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var*&
                &*2+c_var**2*log(c_var))/f_var**2/d_var
              !
            end select
            !
          end select
        else if (nb_par == 2) then
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
                fg=-1._ki/3._ki*h_var**3/f_var**3/g_var*(log(z)+log(1._ki-z)+&
                  &z_log(s23,1._ki))+1._ki/3._ki/g_var*log(e_var)/f_var**3*e_v&
                  &ar**3-5._ki/6._ki*c_var/f_var**2
                !
              case(2)
                !
                fg=-1._ki/3._ki/f_var**3*log(c_var)*c_var/d_var**3*(3._ki*d_v&
                  &ar*f_var**2*h_var+3._ki*c_var*d_var**2*f_var+c_var**2*d_v&
                  &ar**2+c_var**2*d_var*f_var+c_var**2*f_var**2)+1._ki/3._ki*&
                  &h_var**3*(d_var**2+d_var*f_var+f_var**2)/d_var**3/f_var*&
                  &*3*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/6._ki*c_var*&
                  &(2._ki*c_var*d_var+2._ki*c_var*f_var-5._ki*d_var**2+5._ki*d_&
                  &var*f_var)/d_var**2/f_var**2
                !
              end select
              !
            case(2)
              !
              select case(flag)
              !
              case(1)
                !
                fg=1._ki/6._ki*z*h_var**2*(-2._ki*f_var*h_var*g_var+f_var**2*&
                  &h_var+2._ki*c_var*f_var*g_var-g_var**2*f_var-f_var**2*g_v&
                  &ar-3._ki*h_var*g_var**2)/f_var**4/g_var**2*(log(z)+log(1.&
                  &_ki-z)+z_log(s23,1._ki))-1._ki/6._ki*z/g_var**2*log(e_var)/f&
                  &_var**2*e_var**3-1._ki/6._ki*z*(g_var*f_var+f_var**2+c_var&
                  &**2+2._ki*c_var*f_var)/f_var**2/g_var-1._ki/2._ki*z*h_var**&
                  &2*(-d_var-h_var+g_var+2._ki*c_var)/f_var**4*(log(z)+log(1&
                  &._ki-z)+z_log(s23,1._ki))+1._ki/6._ki*z/f_var
                !
              case(2)
                !
                fg=1._ki/6._ki*c_var**2*z*(3._ki*d_var*f_var+c_var*d_var+2._ki&
                  &*c_var*f_var)/d_var**3/f_var**2*log(c_var)+1._ki/6._ki*z*h&
                  &_var**2*(-2._ki*c_var*d_var*f_var**2-2._ki*c_var*f_var**3+&
                  &d_var*f_var**2*h_var+d_var**2*f_var**2+d_var*f_var**3-2.&
                  &_ki*c_var*d_var**2*f_var+2._ki*d_var**2*f_var*h_var-2._ki*f&
                  &_var*d_var**3+6._ki*c_var*d_var**3)/d_var**3/f_var**4*(lo&
                  &g(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/6._ki*c_var*z*(2._k&
                  &i*c_var*f_var+c_var*d_var+2._ki*d_var*f_var)/d_var**2/f_v&
                  &ar**2
                !
              end select
              !
            case(3)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/6._ki*(z-1)*h_var**2*(-2._ki*f_var*h_var*g_var+f_va&
                  &r**2*h_var+2._ki*c_var*f_var*g_var-g_var**2*f_var-f_var**&
                  &2*g_var-3._ki*h_var*g_var**2)/f_var**4/g_var**2*(log(z)+l&
                  &og(1._ki-z)+z_log(s23,1._ki))+1._ki/6._ki*(z-1)/g_var**2*log&
                  &(e_var)/f_var**2*e_var**3+1._ki/6._ki*(z-1)*(g_var*f_var+f&
                  &_var**2+c_var**2+2._ki*c_var*f_var)/f_var**2/g_var+1._ki/2&
                  &._ki*(z-1)*h_var**2*(-d_var-h_var+g_var+2._ki*c_var)/f_var&
                  &**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/6._ki*(z-1)&
                  &/f_var
                !
              case(2)
                !
                fg=-1._ki/6._ki*c_var**2*(z-1)*(3._ki*d_var*f_var+c_var*d_var&
                  &+2._ki*c_var*f_var)/d_var**3/f_var**2*log(c_var)-1._ki/6._k&
                  &i*(z-1)*h_var**2*(-2._ki*c_var*d_var*f_var**2-2._ki*c_var*&
                  &f_var**3+d_var*f_var**2*h_var+d_var**2*f_var**2+d_var*f_&
                  &var**3-2._ki*c_var*d_var**2*f_var+2._ki*d_var**2*f_var*h_v&
                  &ar-2._ki*f_var*d_var**3+6._ki*c_var*d_var**3)/d_var**3/f_v&
                  &ar**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/6._ki*(z-&
                  &1)*c_var*(2._ki*c_var*f_var+c_var*d_var+2._ki*d_var*f_var)&
                  &/d_var**2/f_var**2
                !
              end select
              !
            case(4)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/6._ki*h_var**3*(-2._ki*g_var+f_var)/g_var**2/f_var*&
                  &*3*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/6._ki*e_var*&
                  &*2*(-f_var*e_var+2._ki*c_var*g_var-g_var*f_var)/g_var**2/&
                  &f_var**3*log(e_var)+1._ki/6._ki*(g_var+f_var+c_var)**2/f_v&
                  &ar**2/g_var+1._ki/6._ki/f_var**2*(g_var-2._ki*h_var)
                !
              case(2)
                !
                fg=1._ki/6._ki*c_var**2*(2._ki*c_var*d_var+c_var*f_var+3._ki*d&
                  &_var*f_var)/d_var**2/f_var**3*log(c_var)-1._ki/6._ki*h_var&
                  &**3*(2._ki*d_var+f_var)/d_var**2/f_var**3*(log(z)+log(1._k&
                  &i-z)+z_log(s23,1._ki))+1._ki/6._ki/d_var/f_var**2*c_var**2
                !
              end select
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
                fg=-1._ki/3._ki*z**2*h_var*(f_var**6+27._ki*g_var**5*f_var+f_&
                  &var**4*c_var**2-f_var**5*g_var+f_var**4*g_var**2+24._ki*f&
                  &_var**3*g_var**3+42._ki*f_var**2*g_var**4+2._ki*f_var**5*c&
                  &_var+24._ki*f_var**2*g_var**3*c_var-f_var**4*c_var*g_var+&
                  &33._ki*g_var**4*f_var*c_var+6._ki*g_var**3*f_var*c_var**2+&
                  &6._ki*g_var**6+6._ki*g_var**4*c_var**2+12._ki*g_var**5*c_va&
                  &r)/f_var**5/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._ki)&
                  &)+1._ki/3._ki*z**2/g_var**3*log(e_var)/f_var*e_var**3+1._ki&
                  &/6._ki*z**2*(c_var+f_var)*(2._ki*c_var+2._ki*f_var-g_var)/g&
                  &_var**2/f_var+1._ki/3._ki*z**2*(log(z)+log(1._ki-z)+z_log(s&
                  &23,1._ki))*h_var*(18._ki*c_var*d_var**2+6._ki*d_var**2*f_va&
                  &r+6._ki*f_var*c_var**2+3._ki*c_var*d_var*f_var+18._ki*d_var&
                  &**3-18._ki*g_var*d_var**2+6._ki*g_var**2*d_var-6._ki*c_var*&
                  &g_var**2-6._ki*g_var*d_var*f_var+3._ki*g_var**2*f_var+d_va&
                  &r*f_var**2-2._ki*f_var**2*g_var-3._ki*c_var*f_var*g_var+2.&
                  &_ki*c_var*f_var**2)/f_var**5
                !
              case(2)
                !
                fg=-1._ki/3._ki*c_var**3*z**2/d_var**3/f_var*log(c_var)+1._ki&
                  &/3._ki*z**2*h_var*(f_var**4*c_var**2+4._ki*c_var*d_var**3*&
                  &f_var**2-6._ki*c_var**2*d_var**3*f_var-c_var*f_var**4*d_v&
                  &ar+d_var**4*f_var**2-2._ki*d_var**3*f_var**3+f_var**4*d_v&
                  &ar**2-3._ki*d_var**4*c_var*f_var+6._ki*d_var**4*c_var**2)/&
                  &d_var**3/f_var**5*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1&
                  &._ki/6._ki*c_var*z**2*(2._ki*c_var-d_var)/d_var**2/f_var
                !
              end select
              !
            case(3)
              !
              select case(flag)
              !
              case(1)
                !
                fg=1._ki/3._ki*z*(z-1)*h_var*(f_var**6+27._ki*g_var**5*f_var+&
                  &f_var**4*c_var**2-f_var**5*g_var+f_var**4*g_var**2+24._ki&
                  &*f_var**3*g_var**3+42._ki*f_var**2*g_var**4+2._ki*f_var**5&
                  &*c_var+24._ki*f_var**2*g_var**3*c_var-f_var**4*c_var*g_va&
                  &r+33._ki*g_var**4*f_var*c_var+6._ki*g_var**3*f_var*c_var**&
                  &2+6._ki*g_var**6+6._ki*g_var**4*c_var**2+12._ki*g_var**5*c_&
                  &var)/f_var**5/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._k&
                  &i))-1._ki/3._ki*z*(z-1)/g_var**3*log(e_var)/f_var*e_var**3&
                  &-1._ki/6._ki*z*(z-1)*(c_var+f_var)*(2._ki*c_var+2._ki*f_var-&
                  &g_var)/g_var**2/f_var-1._ki/3._ki*z*(log(z)+log(1._ki-z)+z_&
                  &log(s23,1._ki))*h_var*(18._ki*c_var*d_var**2+6._ki*d_var**2&
                  &*f_var+6._ki*f_var*c_var**2+3._ki*c_var*d_var*f_var+18._ki*&
                  &d_var**3-18._ki*g_var*d_var**2+6._ki*g_var**2*d_var-6._ki*c&
                  &_var*g_var**2-6._ki*g_var*d_var*f_var+3._ki*g_var**2*f_var&
                  &+d_var*f_var**2-2._ki*f_var**2*g_var-3._ki*c_var*f_var*g_v&
                  &ar+2._ki*c_var*f_var**2)*(z-1)/f_var**5
                !
              case(2)
                !
                fg=1._ki/3._ki*z*c_var**3*(z-1)/d_var**3/f_var*log(c_var)-1.&
                  &_ki/3._ki*z*(z-1)*h_var*(f_var**4*c_var**2+4._ki*c_var*d_va&
                  &r**3*f_var**2-6._ki*c_var**2*d_var**3*f_var-c_var*f_var**&
                  &4*d_var+d_var**4*f_var**2-2._ki*d_var**3*f_var**3+f_var**&
                  &4*d_var**2-3._ki*d_var**4*c_var*f_var+6._ki*d_var**4*c_var&
                  &**2)/d_var**3/f_var**5*(log(z)+log(1._ki-z)+z_log(s23,1._k&
                  &i))+1._ki/6._ki*z*(z-1)*c_var*(2._ki*c_var-d_var)/d_var**2/&
                  &f_var
                !
              end select
              !
            case(4)
              !
              select case(flag)
              !
              case(1)
                !
                fg=1._ki/6._ki*z*h_var**2*(-2._ki*f_var**2*h_var*g_var+2._ki*f&
                  &_var**3*h_var+c_var*g_var*f_var**2+2._ki*f_var**2*g_var**&
                  &2-2._ki*f_var**3*g_var+2._ki*g_var**2*f_var*h_var-2._ki*g_v&
                  &ar**2*c_var*f_var+4._ki*g_var**3*f_var+6._ki*g_var**3*h_va&
                  &r)/g_var**3/f_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki)&
                  &)+1._ki/6._ki*z*e_var**2*(-2._ki*f_var*e_var+c_var*g_var-2.&
                  &_ki*g_var*f_var)/g_var**3/f_var**2*log(e_var)-1._ki/6._ki*z&
                  &*(4._ki*c_var*f_var**2+2._ki*f_var*c_var**2-g_var**2*f_var&
                  &+f_var**2*g_var+2._ki*f_var**3-c_var**2*g_var)/f_var**2/g&
                  &_var**2+1._ki/2._ki*z*h_var**2*(-2._ki*d_var+c_var-2._ki*h_v&
                  &ar+2._ki*g_var)/f_var**4*(log(z)+log(1._ki-z)+z_log(s23,1.&
                  &_ki))-1._ki/6._ki*z/f_var
                !
              case(2)
                !
                fg=-1._ki/6._ki*z*c_var**3/d_var**2/f_var**2*log(c_var)-1._ki&
                  &/6._ki*z*h_var**2*(-2._ki*c_var*d_var*f_var-c_var*f_var**2&
                  &+2._ki*d_var*f_var*h_var-2._ki*d_var**2*f_var+2._ki*d_var*f&
                  &_var**2+3._ki*c_var*d_var**2)/f_var**4/d_var**2*(log(z)+l&
                  &og(1._ki-z)+z_log(s23,1._ki))-1._ki/6._ki/f_var**2*z*c_var**&
                  &2/d_var
                !
              end select
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
                fg=-1._ki/3._ki*(z-1)**2*h_var*(f_var**6+27._ki*g_var**5*f_va&
                  &r+f_var**4*c_var**2-f_var**5*g_var+f_var**4*g_var**2+24.&
                  &_ki*f_var**3*g_var**3+42._ki*f_var**2*g_var**4+2._ki*f_var*&
                  &*5*c_var+24._ki*f_var**2*g_var**3*c_var-f_var**4*c_var*g_&
                  &var+33._ki*g_var**4*f_var*c_var+6._ki*g_var**3*f_var*c_var&
                  &**2+6._ki*g_var**6+6._ki*g_var**4*c_var**2+12._ki*g_var**5*&
                  &c_var)/f_var**5/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1&
                  &._ki))+1._ki/3._ki*(z-1)**2/g_var**3*log(e_var)/f_var*e_var&
                  &**3+1._ki/6._ki*(z-1)**2*(c_var+f_var)*(2._ki*c_var+2._ki*f_&
                  &var-g_var)/g_var**2/f_var+1._ki/3._ki*(log(z)+log(1._ki-z)+&
                  &z_log(s23,1._ki))*h_var*(z-1)**2*(18._ki*c_var*d_var**2+6.&
                  &_ki*d_var**2*f_var+6._ki*f_var*c_var**2+3._ki*c_var*d_var*f&
                  &_var+18._ki*d_var**3-18._ki*g_var*d_var**2+6._ki*g_var**2*d&
                  &_var-6._ki*c_var*g_var**2-6._ki*g_var*d_var*f_var+3._ki*g_v&
                  &ar**2*f_var+d_var*f_var**2-2._ki*f_var**2*g_var-3._ki*c_va&
                  &r*f_var*g_var+2._ki*c_var*f_var**2)/f_var**5
                !
              case(2)
                !
                fg=-1._ki/3._ki*c_var**3*(z-1)**2/f_var/d_var**3*log(c_var)+&
                  &1._ki/3._ki*(z-1)**2*h_var*(f_var**4*c_var**2+4._ki*c_var*d&
                  &_var**3*f_var**2-6._ki*c_var**2*d_var**3*f_var-c_var*f_va&
                  &r**4*d_var+d_var**4*f_var**2-2._ki*d_var**3*f_var**3+f_va&
                  &r**4*d_var**2-3._ki*d_var**4*c_var*f_var+6._ki*d_var**4*c_&
                  &var**2)/d_var**3/f_var**5*(log(z)+log(1._ki-z)+z_log(s23,&
                  &1._ki))-1._ki/6._ki*(z-1)**2*c_var*(2._ki*c_var-d_var)/d_var&
                  &**2/f_var
                !
              end select
              !
            case(4)
              !
              select case(flag)
              !
              case(1)
                !
                fg=-1._ki/6._ki*(z-1)*h_var**2*(-2._ki*f_var**2*h_var*g_var+2&
                  &._ki*f_var**3*h_var+c_var*g_var*f_var**2+2._ki*f_var**2*g_&
                  &var**2-2._ki*f_var**3*g_var+2._ki*g_var**2*f_var*h_var-2._k&
                  &i*g_var**2*c_var*f_var+4._ki*g_var**3*f_var+6._ki*g_var**3&
                  &*h_var)/f_var**4/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,&
                  &1._ki))-1._ki/6._ki*(z-1)*e_var**2*(-2._ki*f_var*e_var+c_var&
                  &*g_var-2._ki*g_var*f_var)/f_var**2/g_var**3*log(e_var)+1.&
                  &_ki/6._ki*(z-1)*(4._ki*c_var*f_var**2+2._ki*f_var*c_var**2-g&
                  &_var**2*f_var+f_var**2*g_var+2._ki*f_var**3-c_var**2*g_va&
                  &r)/f_var**2/g_var**2-1._ki/2._ki*(z-1)*h_var**2*(-2._ki*d_v&
                  &ar+c_var-2._ki*h_var+2._ki*g_var)/f_var**4*(log(z)+log(1._k&
                  &i-z)+z_log(s23,1._ki))+1._ki/6._ki*(z-1)/f_var
                !
              case(2)
                !
                fg=1._ki/6._ki*c_var**3/d_var**2*(z-1)/f_var**2*log(c_var)+1&
                  &._ki/6._ki*(z-1)*h_var**2*(-2._ki*c_var*d_var*f_var-c_var*f&
                  &_var**2+2._ki*d_var*f_var*h_var-2._ki*d_var**2*f_var+2._ki*&
                  &d_var*f_var**2+3._ki*c_var*d_var**2)/f_var**4/d_var**2*(l&
                  &og(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/6._ki*(z-1)*c_var&
                  &**2/f_var**2/d_var
                !
              end select
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
                fg=-1._ki/3._ki*h_var**3*(-g_var*f_var+f_var**2+g_var**2)/f_&
                  &var**3/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._k&
                  &i/3._ki/g_var**3*log(e_var)/f_var**3*e_var*(-c_var**2*g_v&
                  &ar*f_var+c_var**2*f_var**2+c_var**2*g_var**2+2._ki*c_var*&
                  &f_var**3+c_var*g_var*f_var**2-g_var**2*c_var*f_var+2._ki*&
                  &f_var**3*g_var+f_var**4+f_var**2*g_var**2)+1._ki/6._ki*(g_&
                  &var+f_var+c_var)*(2._ki*f_var**2+g_var*f_var+2._ki*c_var*f&
                  &_var-2._ki*c_var*g_var-g_var**2)/f_var**2/g_var**2+1._ki/6&
                  &._ki/f_var**2*(3._ki*c_var+g_var)
                !
              case(2)
                !
                fg=-1._ki/3._ki*c_var**3/d_var/f_var**3*log(c_var)+1._ki/3._ki&
                  &*(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var**3/f_var**3/&
                  &d_var
                !
              end select
              !
            end select
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
                  fg=-1._ki/4._ki*h_var**4/f_var**4/g_var*(log(z)+log(1._ki-z)+&
                    &z_log(s23,1._ki))+1._ki/4._ki/g_var*log(e_var)/f_var**4*e_v&
                    &ar**4-1._ki/12._ki*c_var/f_var**3*(-13._ki*g_var+13._ki*f_va&
                    &r+21._ki*c_var)
                  !
                case(2)
                  !
                  fg=-1._ki/4._ki/f_var**4*log(c_var)*c_var/d_var**4*(2._ki*d_v&
                    &ar*f_var+c_var*d_var+c_var*f_var)*(2._ki*d_var*f_var**2*h&
                    &_var+c_var**2*d_var**2+2._ki*c_var*d_var**2*f_var+c_var**&
                    &2*f_var**2)+1._ki/4._ki*h_var**4*(d_var+f_var)*(d_var**2+f&
                    &_var**2)/f_var**4/d_var**4*(log(z)+log(1._ki-z)+z_log(s23&
                    &,1._ki))-1._ki/24._ki*c_var*(6._ki*c_var**2*d_var**2+6._ki*c_&
                    &var**2*d_var*f_var+6._ki*c_var**2*f_var**2-42._ki*c_var*d_&
                    &var**3+21._ki*c_var*d_var**2*f_var+21._ki*c_var*f_var**2*d&
                    &_var+26._ki*d_var**4-52._ki*d_var**3*f_var+26._ki*f_var**2*&
                    &d_var**2)/f_var**3/d_var**3
                  !
                end select
                !
              case(2)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*z*h_var**3*(-3._ki*h_var*f_var*g_var+f_var**2&
                    &*h_var+3._ki*c_var*f_var*g_var+3._ki*g_var**2*f_var-g_var*&
                    &f_var**2-12._ki*c_var*g_var**2)/f_var**5/g_var**2*(log(z)&
                    &+log(1._ki-z)+z_log(s23,1._ki))-1._ki/12._ki*z/g_var**2*log(&
                    &e_var)/f_var**3*e_var**4-1._ki/12._ki*z/f_var**3*(f_var**3&
                    &+c_var**3+3._ki*c_var*h_var*f_var-6._ki*c_var*f_var*g_var)&
                    &/g_var
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*z*c_var**2*(4._ki*c_var*d_var**2*f_var+8._ki*c&
                    &_var*f_var**2*d_var+6._ki*f_var**2*d_var**2+c_var**2*d_va&
                    &r**2+2._ki*c_var**2*d_var*f_var+3._ki*c_var**2*f_var**2)/d&
                    &_var**4/f_var**3*log(c_var)+1._ki/12._ki*z*h_var**3*(-3._ki&
                    &*c_var*f_var**3*d_var-3._ki*c_var*f_var**4+f_var**3*d_var&
                    &*h_var+f_var**3*d_var**2+f_var**4*d_var-3._ki*c_var*d_var&
                    &**2*f_var**2+2._ki*d_var**2*f_var**2*h_var+f_var**2*d_var&
                    &**3-3._ki*d_var**3*f_var*c_var+3._ki*f_var*d_var**3*h_var-&
                    &3._ki*f_var*d_var**4+12._ki*c_var*d_var**4)/d_var**4/f_var&
                    &**5*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/24._ki*z*c_&
                    &var*(4._ki*c_var**2*d_var*f_var+6._ki*c_var**2*f_var**2+6.&
                    &_ki*c_var*d_var**2*f_var+13._ki*c_var*f_var**2*d_var+2._ki*&
                    &c_var**2*d_var**2-6._ki*d_var**3*f_var+6._ki*f_var**2*d_va&
                    &r**2)/f_var**3/d_var**3
                  !
                end select
                !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*(-1._ki+z)*h_var**3*(-3._ki*h_var*f_var*g_var&
                    &+f_var**2*h_var+3._ki*c_var*f_var*g_var+3._ki*g_var**2*f_v&
                    &ar-g_var*f_var**2-12._ki*c_var*g_var**2)/f_var**5/g_var**&
                    &2*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/12._ki*(-1._ki&
                    &+z)/g_var**2*log(e_var)/f_var**3*e_var**4+1._ki/12._ki*(-1&
                    &._ki+z)/f_var**3*(f_var**3+c_var**3+3._ki*c_var*h_var*f_va&
                    &r-6._ki*c_var*f_var*g_var)/g_var
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*c_var**2*(4._ki*c_var*d_var**2*f_var+8._ki*c_&
                    &var*f_var**2*d_var+6._ki*f_var**2*d_var**2+c_var**2*d_var&
                    &**2+2._ki*c_var**2*d_var*f_var+3._ki*c_var**2*f_var**2)*(-&
                    &1._ki+z)/d_var**4/f_var**3*log(c_var)-1._ki/12._ki*(-1._ki+z&
                    &)*h_var**3*(-3._ki*c_var*f_var**3*d_var-3._ki*c_var*f_var*&
                    &*4+f_var**3*d_var*h_var+f_var**3*d_var**2+f_var**4*d_var&
                    &-3._ki*c_var*d_var**2*f_var**2+2._ki*d_var**2*f_var**2*h_v&
                    &ar+f_var**2*d_var**3-3._ki*d_var**3*f_var*c_var+3._ki*f_va&
                    &r*d_var**3*h_var-3._ki*f_var*d_var**4+12._ki*c_var*d_var**&
                    &4)/d_var**4/f_var**5*(log(z)+log(1._ki-z)+z_log(s23,1._ki)&
                    &)-1._ki/24._ki*(-1._ki+z)*c_var*(4._ki*c_var**2*d_var*f_var+&
                    &6._ki*c_var**2*f_var**2+6._ki*c_var*d_var**2*f_var+13._ki*c&
                    &_var*f_var**2*d_var+2._ki*c_var**2*d_var**2-6._ki*d_var**3&
                    &*f_var+6._ki*f_var**2*d_var**2)/d_var**3/f_var**3
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*h_var**4*(-3._ki*g_var+f_var)/g_var**2/f_var&
                    &**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/12._ki*e_va&
                    &r**3*(-f_var*e_var+3._ki*c_var*g_var-g_var*f_var)/g_var**&
                    &2/f_var**4*log(e_var)+1._ki/24._ki*(-6._ki*g_var**3-6._ki*g_&
                    &var**2*f_var+2._ki*f_var**3+6._ki*c_var*h_var*g_var+6._ki*c&
                    &_var*h_var*f_var+2._ki*c_var**3+c_var**2*g_var-12._ki*c_va&
                    &r*g_var**2-12._ki*c_var*f_var*g_var+6._ki*g_var**2*h_var)/&
                    &f_var**3/g_var
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*c_var**2*(3._ki*c_var**2*d_var**2+2._ki*c_var*&
                    &*2*d_var*f_var+c_var**2*f_var**2+8._ki*c_var*d_var**2*f_v&
                    &ar+4._ki*c_var*f_var**2*d_var+6._ki*f_var**2*d_var**2)/d_v&
                    &ar**3/f_var**4*log(c_var)-1._ki/12._ki*h_var**4*(3._ki*d_va&
                    &r**2+2._ki*d_var*f_var+f_var**2)/d_var**3/f_var**4*(log(z&
                    &)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/24._ki*c_var**2*(4._ki&
                    &*c_var*d_var+2._ki*c_var*f_var-7._ki*d_var**2+7._ki*d_var*f&
                    &_var)/f_var**3/d_var**2
                  !
                end select
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
                  fg=-1._ki/12._ki*z**2*h_var**2*(-4._ki*f_var**3*g_var**3+f_va&
                    &r**4*c_var**2+2._ki*c_var*f_var**5+6._ki*f_var**2*g_var**3&
                    &*c_var-2._ki*f_var**4*g_var*c_var-12._ki*g_var**4*f_var*c_&
                    &var+f_var**6+3._ki*f_var**2*g_var**4+3._ki*f_var**4*g_var*&
                    &*2-2._ki*f_var**5*g_var+30._ki*g_var**4*c_var**2)/f_var**6&
                    &/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/12._k&
                    &i*z**2/g_var**3*log(e_var)/f_var**2*e_var**4+1._ki/24._ki*&
                    &z**2*(-g_var*f_var**2+2._ki*f_var**3+6._ki*c_var*h_var*f_v&
                    &ar+2._ki*c_var**3-c_var**2*g_var-8._ki*c_var*f_var*g_var)/&
                    &f_var**2/g_var**2
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*z**2*c_var**3*(c_var*d_var+4._ki*d_var*f_var&
                    &+3._ki*c_var*f_var)/d_var**4/f_var**2*log(c_var)+1._ki/12.&
                    &_ki*z**2*h_var**2*(18._ki*c_var*d_var**4*f_var**2-30._ki*d_&
                    &var**4*c_var**2*f_var-12._ki*d_var**5*c_var*f_var-2._ki*c_&
                    &var*f_var**5*d_var+f_var**5*d_var**2+f_var**4*c_var**2*d&
                    &_var-2._ki*f_var**4*c_var*d_var**2-7._ki*f_var**3*d_var**4&
                    &+3._ki*f_var**4*d_var**3+3._ki*d_var**5*f_var**2+3._ki*c_va&
                    &r**2*f_var**5+30._ki*d_var**5*c_var**2)/f_var**6/d_var**4&
                    &*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/24._ki*c_var*z&
                    &**2*(-2._ki*d_var**2*f_var-c_var*d_var**2+6._ki*f_var*c_va&
                    &r**2+5._ki*d_var*c_var*f_var+2._ki*d_var*c_var**2)/f_var**&
                    &2/d_var**3
                  !
                end select
                !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*z*(-1._ki+z)*h_var**2*(-4._ki*f_var**3*g_var**&
                    &3+f_var**4*c_var**2+2._ki*c_var*f_var**5+6._ki*f_var**2*g_&
                    &var**3*c_var-2._ki*f_var**4*g_var*c_var-12._ki*g_var**4*f_&
                    &var*c_var+f_var**6+3._ki*f_var**2*g_var**4+3._ki*f_var**4*&
                    &g_var**2-2._ki*f_var**5*g_var+30._ki*g_var**4*c_var**2)/f_&
                    &var**6/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._k&
                    &i/12._ki*z*(-1._ki+z)/g_var**3*log(e_var)/f_var**2*e_var**&
                    &4-1._ki/24._ki*z*(-1._ki+z)*(-g_var*f_var**2+2._ki*f_var**3+&
                    &6._ki*c_var*h_var*f_var+2._ki*c_var**3-c_var**2*g_var-8._ki&
                    &*c_var*f_var*g_var)/f_var**2/g_var**2
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*z*c_var**3*(c_var*d_var+4._ki*d_var*f_var+3._k&
                    &i*c_var*f_var)*(-1._ki+z)/d_var**4/f_var**2*log(c_var)-1.&
                    &_ki/12._ki*z*(-1._ki+z)*h_var**2*(18._ki*c_var*d_var**4*f_va&
                    &r**2-30._ki*d_var**4*c_var**2*f_var-12._ki*d_var**5*c_var*&
                    &f_var-2._ki*c_var*f_var**5*d_var+f_var**5*d_var**2+f_var*&
                    &*4*c_var**2*d_var-2._ki*f_var**4*c_var*d_var**2-7._ki*f_va&
                    &r**3*d_var**4+3._ki*f_var**4*d_var**3+3._ki*d_var**5*f_var&
                    &**2+3._ki*c_var**2*f_var**5+30._ki*d_var**5*c_var**2)/f_va&
                    &r**6/d_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/&
                    &24._ki*z*(-1._ki+z)*c_var*(-2._ki*d_var**2*f_var-c_var*d_va&
                    &r**2+6._ki*f_var*c_var**2+5._ki*d_var*c_var*f_var+2._ki*d_v&
                    &ar*c_var**2)/f_var**2/d_var**3
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*z*h_var**3*(-2._ki*f_var**2*h_var*g_var+f_var&
                    &**3*h_var+g_var*f_var**2*c_var+2._ki*f_var**2*g_var**2-f_&
                    &var**3*g_var+3._ki*g_var**2*f_var*h_var-3._ki*g_var**2*c_v&
                    &ar*f_var-3._ki*g_var**3*f_var+6._ki*g_var**3*c_var)/f_var*&
                    &*5/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/12&
                    &._ki*z*e_var**3*(-f_var*e_var+c_var*g_var-g_var*f_var)/f_&
                    &var**3/g_var**3*log(e_var)-1._ki/24._ki*z/f_var**3*(-6._ki*&
                    &g_var**2*c_var*f_var+6._ki*f_var*c_var*h_var*g_var+6._ki*c&
                    &_var*h_var*f_var**2+2._ki*f_var*c_var**3+f_var**3*g_var+2&
                    &._ki*f_var**4-2._ki*g_var*c_var**3-9._ki*g_var*f_var*c_var*&
                    &*2-12._ki*g_var*f_var**2*c_var)/g_var**2
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*z*c_var**3*(2._ki*d_var*f_var+c_var*d_var+c_&
                    &var*f_var)/d_var**3/f_var**3*log(c_var)-1._ki/12._ki*z*h_v&
                    &ar**3*(-2._ki*c_var*f_var**2*d_var-c_var*f_var**3+d_var*f&
                    &_var**2*h_var+2._ki*f_var**2*d_var**2+f_var**3*d_var-3._ki&
                    &*c_var*d_var**2*f_var+3._ki*d_var**2*f_var*h_var-3._ki*d_v&
                    &ar**3*f_var+6._ki*c_var*d_var**3)/d_var**3/f_var**5*(log(&
                    &z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/24._ki*z*c_var**2*(2&
                    &._ki*c_var*f_var+2._ki*c_var*d_var+3._ki*d_var*f_var)/f_var&
                    &**3/d_var**2
                  !
                end select
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
                  fg=-1._ki/12._ki*(-1._ki+z)**2*h_var**2*(-4._ki*f_var**3*g_var&
                    &**3+f_var**4*c_var**2+2._ki*c_var*f_var**5+6._ki*f_var**2*&
                    &g_var**3*c_var-2._ki*f_var**4*g_var*c_var-12._ki*g_var**4*&
                    &f_var*c_var+f_var**6+3._ki*f_var**2*g_var**4+3._ki*f_var**&
                    &4*g_var**2-2._ki*f_var**5*g_var+30._ki*g_var**4*c_var**2)/&
                    &f_var**6/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1&
                    &._ki/12._ki*(-1._ki+z)**2/g_var**3*log(e_var)/f_var**2*e_va&
                    &r**4+1._ki/24._ki*(-1._ki+z)**2*(-g_var*f_var**2+2._ki*f_var&
                    &**3+6._ki*c_var*h_var*f_var+2._ki*c_var**3-c_var**2*g_var-&
                    &8._ki*c_var*f_var*g_var)/f_var**2/g_var**2
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*c_var**3*(-1._ki+z)**2*(c_var*d_var+4._ki*d_v&
                    &ar*f_var+3._ki*c_var*f_var)/f_var**2/d_var**4*log(c_var)+&
                    &1._ki/12._ki*(-1._ki+z)**2*h_var**2*(18._ki*c_var*d_var**4*f&
                    &_var**2-30._ki*d_var**4*c_var**2*f_var-12._ki*d_var**5*c_v&
                    &ar*f_var-2._ki*c_var*f_var**5*d_var+f_var**5*d_var**2+f_v&
                    &ar**4*c_var**2*d_var-2._ki*f_var**4*c_var*d_var**2-7._ki*f&
                    &_var**3*d_var**4+3._ki*f_var**4*d_var**3+3._ki*d_var**5*f_&
                    &var**2+3._ki*c_var**2*f_var**5+30._ki*d_var**5*c_var**2)/f&
                    &_var**6/d_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1.&
                    &_ki/24._ki*(-1._ki+z)**2*c_var*(-2._ki*d_var**2*f_var-c_var*&
                    &d_var**2+6._ki*f_var*c_var**2+5._ki*d_var*c_var*f_var+2._ki&
                    &*d_var*c_var**2)/f_var**2/d_var**3
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*(-1._ki+z)*h_var**3*(-2._ki*f_var**2*h_var*g_&
                    &var+f_var**3*h_var+g_var*f_var**2*c_var+2._ki*f_var**2*g_&
                    &var**2-f_var**3*g_var+3._ki*g_var**2*f_var*h_var-3._ki*g_v&
                    &ar**2*c_var*f_var-3._ki*g_var**3*f_var+6._ki*g_var**3*c_va&
                    &r)/f_var**5/g_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._ki)&
                    &)-1._ki/12._ki*(-1._ki+z)*e_var**3*(-f_var*e_var+c_var*g_va&
                    &r-g_var*f_var)/f_var**3/g_var**3*log(e_var)+1._ki/24._ki*(&
                    &-1._ki+z)/f_var**3*(-6._ki*g_var**2*c_var*f_var+6._ki*f_var&
                    &*c_var*h_var*g_var+6._ki*c_var*h_var*f_var**2+2._ki*f_var*&
                    &c_var**3+f_var**3*g_var+2._ki*f_var**4-2._ki*g_var*c_var**&
                    &3-9._ki*g_var*f_var*c_var**2-12._ki*g_var*f_var**2*c_var)/&
                    &g_var**2
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*c_var**3*(2._ki*d_var*f_var+c_var*d_var+c_var&
                    &*f_var)*(-1._ki+z)/d_var**3/f_var**3*log(c_var)+1._ki/12._k&
                    &i*(-1._ki+z)*h_var**3*(-2._ki*c_var*f_var**2*d_var-c_var*f&
                    &_var**3+d_var*f_var**2*h_var+2._ki*f_var**2*d_var**2+f_va&
                    &r**3*d_var-3._ki*c_var*d_var**2*f_var+3._ki*d_var**2*f_var&
                    &*h_var-3._ki*d_var**3*f_var+6._ki*c_var*d_var**3)/f_var**5&
                    &/d_var**3*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/24._k&
                    &i*(-1._ki+z)*c_var**2*(2._ki*c_var*f_var+2._ki*c_var*d_var+&
                    &3._ki*d_var*f_var)/f_var**3/d_var**2
                  !
                end select
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
                  fg=-1._ki/12._ki*h_var**4*(f_var**2-2._ki*g_var*f_var+3._ki*g_&
                    &var**2)/g_var**3/f_var**4*(log(z)+log(1._ki-z)+z_log(s23,&
                    &1._ki))+1._ki/12._ki/g_var**3*log(e_var)/f_var**4*e_var**2*&
                    &(c_var**2*f_var**2-2._ki*g_var*f_var*c_var**2+3._ki*c_var*&
                    &*2*g_var**2-2._ki*g_var**2*c_var*f_var+2._ki*c_var*f_var**&
                    &3+f_var**4+2._ki*f_var**3*g_var+f_var**2*g_var**2)+1._ki/2&
                    &4._ki*(-6._ki*c_var**2*g_var**2-16._ki*g_var*f_var**2*c_var&
                    &+6._ki*c_var*h_var*f_var**2-17._ki*g_var*f_var*c_var**2-18&
                    &._ki*g_var**2*c_var*f_var+2._ki*f_var*c_var**3+12._ki*f_var&
                    &*c_var*h_var*g_var-4._ki*g_var*c_var**3-6._ki*g_var**3*c_v&
                    &ar+6._ki*c_var*h_var*g_var**2+3._ki*f_var**3*g_var+2._ki*f_&
                    &var**4)/f_var**3/g_var**2
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*c_var**3*(3._ki*c_var*d_var+c_var*f_var+4._ki&
                    &*d_var*f_var)/f_var**4/d_var**2*log(c_var)+1._ki/12._ki*h_&
                    &var**4*(3._ki*d_var+f_var)/f_var**4/d_var**2*(log(z)+log(&
                    &1._ki-z)+z_log(s23,1._ki))-1._ki/12._ki/d_var/f_var**3*c_var&
                    &**3
                  !
                end select
                !
              end select
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
                  fg=1._ki/4._ki*z**3*h_var*(2._ki*f_var**7*h_var*c_var+4._ki*f_&
                    &var**4*h_var*g_var**3*c_var+f_var**6*h_var*c_var**2+2._ki&
                    &*g_var**4*f_var**3*h_var**2-2._ki*g_var*f_var**6*h_var**2&
                    &+2._ki*g_var**2*f_var**5*h_var**2-2._ki*g_var**3*f_var**4*&
                    &h_var**2-3._ki*f_var**3*h_var*g_var**5+3._ki*f_var**4*h_va&
                    &r*g_var**4-3._ki*f_var**5*h_var*g_var**3+3._ki*f_var**6*h_&
                    &var*g_var**2+f_var**7*h_var*g_var+f_var**8*h_var-4._ki*f_&
                    &var**3*h_var*g_var**3*c_var**2+2._ki*f_var**4*h_var*c_var&
                    &**2*g_var**2-f_var**5*h_var*g_var*c_var**2-2._ki*f_var**5&
                    &*h_var*c_var*g_var**2+6._ki*g_var**4*f_var**2*h_var*c_var&
                    &**2+8._ki*g_var**5*f_var**2*h_var*c_var-10._ki*g_var**5*f_&
                    &var*h_var*c_var**2-6._ki*f_var**3*h_var*c_var*g_var**4+g_&
                    &var**3*f_var**5*c_var-g_var**2*f_var**5*c_var**2+10._ki*g&
                    &_var**6*c_var**2*f_var+10._ki*g_var**5*c_var**3*f_var-4._k&
                    &i*g_var**6*f_var**2*c_var+3._ki*f_var**3*g_var**5*c_var-8&
                    &._ki*g_var**4*f_var**2*c_var**3-6._ki*g_var**5*f_var**2*c_&
                    &var**2-g_var*f_var**7*c_var-2._ki*g_var**2*f_var**4*c_var&
                    &**3-2._ki*g_var**4*f_var**4*c_var+g_var*f_var**5*c_var**3&
                    &+g_var*f_var**6*c_var**2+2._ki*g_var**4*f_var**3*c_var**2&
                    &+4._ki*g_var**3*f_var**3*c_var**3-g_var**5*f_var**4-g_var&
                    &*f_var**8-20._ki*g_var**6*c_var**3-g_var**3*f_var**6-g_va&
                    &r**2*f_var**7+g_var**4*f_var**5+g_var**6*f_var**3)/f_var&
                    &**7/g_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/4&
                    &._ki*z**3/g_var**4*log(e_var)/f_var*e_var**4-1._ki/24._ki*z&
                    &**3*(2._ki*g_var**2*f_var-3._ki*g_var*f_var**2+6._ki*f_var*&
                    &*3+6._ki*c_var**3+18._ki*c_var*h_var*f_var-24._ki*c_var*f_v&
                    &ar*g_var-3._ki*c_var**2*g_var+2._ki*c_var*g_var**2)/g_var*&
                    &*3/f_var
                  !
                case(2)
                  !
                  fg=1._ki/4._ki*z**3*c_var**4/d_var**4/f_var*log(c_var)+1._ki/&
                    &4._ki*z**3*h_var*(2._ki*d_var**3*f_var**4*h_var**2+3._ki*d_&
                    &var**3*f_var**5*h_var-2._ki*d_var**4*f_var**3*h_var**2-6.&
                    &_ki*d_var**4*f_var**4*h_var+3._ki*f_var**3*d_var**5*h_var-&
                    &4._ki*c_var*d_var**3*f_var**4*h_var+10._ki*d_var**5*c_var*&
                    &*2*f_var*h_var+4._ki*f_var**3*d_var**3*c_var**2*h_var+2._k&
                    &i*f_var**4*d_var**2*c_var**2*h_var-8._ki*d_var**5*f_var**&
                    &2*h_var*c_var-f_var**6*c_var**3+20._ki*d_var**6*c_var**3-&
                    &3._ki*d_var**4*f_var**5+3._ki*d_var**5*f_var**4-d_var**6*f&
                    &_var**3+d_var**3*f_var**6+f_var**5*d_var*c_var**2*h_var-&
                    &2._ki*c_var*d_var**2*f_var**5*h_var-16._ki*c_var**2*d_var*&
                    &*4*f_var**2*h_var+14._ki*c_var*d_var**4*f_var**3*h_var+4.&
                    &_ki*d_var**6*f_var**2*c_var-c_var*d_var**2*f_var**6-50._ki&
                    &*c_var**3*d_var**5*f_var-c_var*d_var**3*f_var**5+9._ki*c_&
                    &var*d_var**4*f_var**4-11._ki*c_var*d_var**5*f_var**3+26._k&
                    &i*c_var**2*d_var**5*f_var**2-10._ki*d_var**6*c_var**2*f_v&
                    &ar-4._ki*f_var**3*c_var**3*d_var**3-2._ki*f_var**4*c_var**&
                    &3*d_var**2-f_var**5*c_var**3*d_var-18._ki*f_var**3*d_var*&
                    &*4*c_var**2+f_var**5*d_var**2*c_var**2+f_var**6*d_var*c_&
                    &var**2+38._ki*f_var**2*d_var**4*c_var**3)/f_var**7/d_var*&
                    &*4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/24._ki*c_var&
                    &*z**3*(6._ki*c_var**2-3._ki*c_var*d_var+2._ki*d_var**2)/d_v&
                    &ar**3/f_var
                  !
                end select
                !
              case(3)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/4._ki*z**2*(-1._ki+z)*h_var*(2._ki*f_var**7*h_var*c_&
                    &var+4._ki*f_var**4*h_var*g_var**3*c_var+f_var**6*h_var*c_&
                    &var**2+2._ki*g_var**4*f_var**3*h_var**2-2._ki*g_var*f_var*&
                    &*6*h_var**2+2._ki*g_var**2*f_var**5*h_var**2-2._ki*g_var**&
                    &3*f_var**4*h_var**2-3._ki*f_var**3*h_var*g_var**5+3._ki*f_&
                    &var**4*h_var*g_var**4-3._ki*f_var**5*h_var*g_var**3+3._ki*&
                    &f_var**6*h_var*g_var**2+f_var**7*h_var*g_var+f_var**8*h_&
                    &var-4._ki*f_var**3*h_var*g_var**3*c_var**2+2._ki*f_var**4*&
                    &h_var*c_var**2*g_var**2-f_var**5*h_var*g_var*c_var**2-2.&
                    &_ki*f_var**5*h_var*c_var*g_var**2+6._ki*g_var**4*f_var**2*&
                    &h_var*c_var**2+8._ki*g_var**5*f_var**2*h_var*c_var-10._ki*&
                    &g_var**5*f_var*h_var*c_var**2-6._ki*f_var**3*h_var*c_var*&
                    &g_var**4+g_var**3*f_var**5*c_var-g_var**2*f_var**5*c_var&
                    &**2+10._ki*g_var**6*c_var**2*f_var+10._ki*g_var**5*c_var**&
                    &3*f_var-4._ki*g_var**6*f_var**2*c_var+3._ki*f_var**3*g_var&
                    &**5*c_var-8._ki*g_var**4*f_var**2*c_var**3-6._ki*g_var**5*&
                    &f_var**2*c_var**2-g_var*f_var**7*c_var-2._ki*g_var**2*f_v&
                    &ar**4*c_var**3-2._ki*g_var**4*f_var**4*c_var+g_var*f_var*&
                    &*5*c_var**3+g_var*f_var**6*c_var**2+2._ki*g_var**4*f_var*&
                    &*3*c_var**2+4._ki*g_var**3*f_var**3*c_var**3-g_var**5*f_v&
                    &ar**4-g_var*f_var**8-20._ki*g_var**6*c_var**3-g_var**3*f_&
                    &var**6-g_var**2*f_var**7+g_var**4*f_var**5+g_var**6*f_va&
                    &r**3)/f_var**7/g_var**4*(log(z)+log(1._ki-z)+z_log(s23,1.&
                    &_ki))+1._ki/4._ki*z**2*(-1._ki+z)/g_var**4*log(e_var)/f_var*&
                    &e_var**4+1._ki/24._ki*z**2*(2._ki*g_var**2*f_var-3._ki*g_var&
                    &*f_var**2+6._ki*f_var**3+6._ki*c_var**3+18._ki*c_var*h_var*&
                    &f_var-24._ki*c_var*f_var*g_var-3._ki*c_var**2*g_var+2._ki*c&
                    &_var*g_var**2)*(-1._ki+z)/g_var**3/f_var
                  !
                case(2)
                  !
                  fg=-1._ki/4._ki*z**2*c_var**4*(-1._ki+z)/f_var/d_var**4*log(c&
                    &_var)-1._ki/4._ki*z**2*(-1._ki+z)*h_var*(2._ki*d_var**3*f_va&
                    &r**4*h_var**2+3._ki*d_var**3*f_var**5*h_var-2._ki*d_var**4&
                    &*f_var**3*h_var**2-6._ki*d_var**4*f_var**4*h_var+3._ki*f_v&
                    &ar**3*d_var**5*h_var-4._ki*c_var*d_var**3*f_var**4*h_var+&
                    &10._ki*d_var**5*c_var**2*f_var*h_var+4._ki*f_var**3*d_var*&
                    &*3*c_var**2*h_var+2._ki*f_var**4*d_var**2*c_var**2*h_var-&
                    &8._ki*d_var**5*f_var**2*h_var*c_var-f_var**6*c_var**3+20.&
                    &_ki*d_var**6*c_var**3-3._ki*d_var**4*f_var**5+3._ki*d_var**&
                    &5*f_var**4-d_var**6*f_var**3+d_var**3*f_var**6+f_var**5*&
                    &d_var*c_var**2*h_var-2._ki*c_var*d_var**2*f_var**5*h_var-&
                    &16._ki*c_var**2*d_var**4*f_var**2*h_var+14._ki*c_var*d_var&
                    &**4*f_var**3*h_var+4._ki*d_var**6*f_var**2*c_var-c_var*d_&
                    &var**2*f_var**6-50._ki*c_var**3*d_var**5*f_var-c_var*d_va&
                    &r**3*f_var**5+9._ki*c_var*d_var**4*f_var**4-11._ki*c_var*d&
                    &_var**5*f_var**3+26._ki*c_var**2*d_var**5*f_var**2-10._ki*&
                    &d_var**6*c_var**2*f_var-4._ki*f_var**3*c_var**3*d_var**3-&
                    &2._ki*f_var**4*c_var**3*d_var**2-f_var**5*c_var**3*d_var-&
                    &18._ki*f_var**3*d_var**4*c_var**2+f_var**5*d_var**2*c_var&
                    &**2+f_var**6*d_var*c_var**2+38._ki*f_var**2*d_var**4*c_va&
                    &r**3)/f_var**7/d_var**4*(log(z)+log(1._ki-z)+z_log(s23,1.&
                    &_ki))-1._ki/24._ki*z**2*(-1._ki+z)*c_var*(6._ki*c_var**2-3._ki&
                    &*c_var*d_var+2._ki*d_var**2)/d_var**3/f_var
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*z**2*h_var**2*(3._ki*f_var**7-3._ki*f_var**6*&
                    &g_var+3._ki*f_var**5*g_var**2-3._ki*f_var**4*g_var**3+3._ki&
                    &*f_var**3*g_var**4+6._ki*f_var**6*c_var+2._ki*f_var**4*g_v&
                    &ar**2*c_var-f_var**4*g_var*c_var**2-4._ki*f_var**5*g_var*&
                    &c_var-3._ki*g_var**5*f_var**2-10._ki*g_var**5*c_var**2+3._k&
                    &i*c_var**2*f_var**5-2._ki*g_var**4*f_var**2*c_var-2._ki*g_&
                    &var**4*f_var*c_var**2+8._ki*g_var**5*f_var*c_var)/f_var**&
                    &6/g_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/12.&
                    &_ki*z**2*e_var**3*(-3._ki*f_var*e_var+c_var*g_var-3._ki*g_v&
                    &ar*f_var)/f_var**2/g_var**4*log(e_var)+1._ki/24._ki*z**2*(&
                    &-18._ki*g_var**2*c_var*f_var+c_var**2*g_var**2-f_var**2*g&
                    &_var**2+18._ki*f_var*c_var*h_var*g_var+18._ki*c_var*h_var*&
                    &f_var**2+6._ki*f_var*c_var**3+3._ki*f_var**3*g_var+6._ki*f_&
                    &var**4-32._ki*g_var*f_var**2*c_var-19._ki*g_var*f_var*c_va&
                    &r**2-2._ki*g_var*c_var**3)/f_var**2/g_var**3
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*c_var**4*z**2/d_var**3/f_var**2*log(c_var)-1&
                    &._ki/12._ki*z**2*h_var**2*(f_var**4*c_var**2+10._ki*d_var**&
                    &3*c_var*f_var**2-8._ki*d_var**3*c_var**2*f_var-2._ki*d_var&
                    &*c_var*f_var**4+3._ki*f_var**2*d_var**4-6._ki*d_var**3*f_v&
                    &ar**3+3._ki*f_var**4*d_var**2-8._ki*d_var**4*c_var*f_var+1&
                    &0._ki*d_var**4*c_var**2)/d_var**3/f_var**6*(log(z)+log(1.&
                    &_ki-z)+z_log(s23,1._ki))+1._ki/24._ki*c_var**2*z**2*(2._ki*c_&
                    &var-d_var)/f_var**2/d_var**2
                  !
                end select
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
                  fg=1._ki/4._ki*z*(-1._ki+z)**2*h_var*(2._ki*f_var**7*h_var*c_v&
                    &ar+4._ki*f_var**4*h_var*g_var**3*c_var+f_var**6*h_var*c_v&
                    &ar**2+2._ki*g_var**4*f_var**3*h_var**2-2._ki*g_var*f_var**&
                    &6*h_var**2+2._ki*g_var**2*f_var**5*h_var**2-2._ki*g_var**3&
                    &*f_var**4*h_var**2-3._ki*f_var**3*h_var*g_var**5+3._ki*f_v&
                    &ar**4*h_var*g_var**4-3._ki*f_var**5*h_var*g_var**3+3._ki*f&
                    &_var**6*h_var*g_var**2+f_var**7*h_var*g_var+f_var**8*h_v&
                    &ar-4._ki*f_var**3*h_var*g_var**3*c_var**2+2._ki*f_var**4*h&
                    &_var*c_var**2*g_var**2-f_var**5*h_var*g_var*c_var**2-2._k&
                    &i*f_var**5*h_var*c_var*g_var**2+6._ki*g_var**4*f_var**2*h&
                    &_var*c_var**2+8._ki*g_var**5*f_var**2*h_var*c_var-10._ki*g&
                    &_var**5*f_var*h_var*c_var**2-6._ki*f_var**3*h_var*c_var*g&
                    &_var**4+g_var**3*f_var**5*c_var-g_var**2*f_var**5*c_var*&
                    &*2+10._ki*g_var**6*c_var**2*f_var+10._ki*g_var**5*c_var**3&
                    &*f_var-4._ki*g_var**6*f_var**2*c_var+3._ki*f_var**3*g_var*&
                    &*5*c_var-8._ki*g_var**4*f_var**2*c_var**3-6._ki*g_var**5*f&
                    &_var**2*c_var**2-g_var*f_var**7*c_var-2._ki*g_var**2*f_va&
                    &r**4*c_var**3-2._ki*g_var**4*f_var**4*c_var+g_var*f_var**&
                    &5*c_var**3+g_var*f_var**6*c_var**2+2._ki*g_var**4*f_var**&
                    &3*c_var**2+4._ki*g_var**3*f_var**3*c_var**3-g_var**5*f_va&
                    &r**4-g_var*f_var**8-20._ki*g_var**6*c_var**3-g_var**3*f_v&
                    &ar**6-g_var**2*f_var**7+g_var**4*f_var**5+g_var**6*f_var&
                    &**3)/f_var**7/g_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._k&
                    &i))-1._ki/4._ki*z*(-1._ki+z)**2/g_var**4*log(e_var)/f_var*e&
                    &_var**4-1._ki/24._ki*z*(-1._ki+z)**2*(2._ki*g_var**2*f_var-3&
                    &._ki*g_var*f_var**2+6._ki*f_var**3+6._ki*c_var**3+18._ki*c_v&
                    &ar*h_var*f_var-24._ki*c_var*f_var*g_var-3._ki*c_var**2*g_v&
                    &ar+2._ki*c_var*g_var**2)/g_var**3/f_var
                  !
                case(2)
                  !
                  fg=1._ki/4._ki*z*c_var**4*(-1._ki+z)**2/f_var/d_var**4*log(c_&
                    &var)+1._ki/4._ki*z*(-1._ki+z)**2*h_var*(2._ki*d_var**3*f_var&
                    &**4*h_var**2+3._ki*d_var**3*f_var**5*h_var-2._ki*d_var**4*&
                    &f_var**3*h_var**2-6._ki*d_var**4*f_var**4*h_var+3._ki*f_va&
                    &r**3*d_var**5*h_var-4._ki*c_var*d_var**3*f_var**4*h_var+1&
                    &0._ki*d_var**5*c_var**2*f_var*h_var+4._ki*f_var**3*d_var**&
                    &3*c_var**2*h_var+2._ki*f_var**4*d_var**2*c_var**2*h_var-8&
                    &._ki*d_var**5*f_var**2*h_var*c_var-f_var**6*c_var**3+20._k&
                    &i*d_var**6*c_var**3-3._ki*d_var**4*f_var**5+3._ki*d_var**5&
                    &*f_var**4-d_var**6*f_var**3+d_var**3*f_var**6+f_var**5*d&
                    &_var*c_var**2*h_var-2._ki*c_var*d_var**2*f_var**5*h_var-1&
                    &6._ki*c_var**2*d_var**4*f_var**2*h_var+14._ki*c_var*d_var*&
                    &*4*f_var**3*h_var+4._ki*d_var**6*f_var**2*c_var-c_var*d_v&
                    &ar**2*f_var**6-50._ki*c_var**3*d_var**5*f_var-c_var*d_var&
                    &**3*f_var**5+9._ki*c_var*d_var**4*f_var**4-11._ki*c_var*d_&
                    &var**5*f_var**3+26._ki*c_var**2*d_var**5*f_var**2-10._ki*d&
                    &_var**6*c_var**2*f_var-4._ki*f_var**3*c_var**3*d_var**3-2&
                    &._ki*f_var**4*c_var**3*d_var**2-f_var**5*c_var**3*d_var-1&
                    &8._ki*f_var**3*d_var**4*c_var**2+f_var**5*d_var**2*c_var*&
                    &*2+f_var**6*d_var*c_var**2+38._ki*f_var**2*d_var**4*c_var&
                    &**3)/f_var**7/d_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._k&
                    &i))+1._ki/24._ki*z*(-1._ki+z)**2*c_var*(6._ki*c_var**2-3._ki*&
                    &c_var*d_var+2._ki*d_var**2)/d_var**3/f_var
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=1._ki/12._ki*z*(-1._ki+z)*h_var**2*(3._ki*f_var**7-3._ki*f_v&
                    &ar**6*g_var+3._ki*f_var**5*g_var**2-3._ki*f_var**4*g_var**&
                    &3+3._ki*f_var**3*g_var**4+6._ki*f_var**6*c_var+2._ki*f_var*&
                    &*4*g_var**2*c_var-f_var**4*g_var*c_var**2-4._ki*f_var**5*&
                    &g_var*c_var-3._ki*g_var**5*f_var**2-10._ki*g_var**5*c_var*&
                    &*2+3._ki*c_var**2*f_var**5-2._ki*g_var**4*f_var**2*c_var-2&
                    &._ki*g_var**4*f_var*c_var**2+8._ki*g_var**5*f_var*c_var)/f&
                    &_var**6/g_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1.&
                    &_ki/12._ki*z*(-1._ki+z)*e_var**3*(-3._ki*f_var*e_var+c_var*g&
                    &_var-3._ki*g_var*f_var)/f_var**2/g_var**4*log(e_var)-1._ki&
                    &/24._ki*z*(-1._ki+z)*(-18._ki*g_var**2*c_var*f_var+c_var**2&
                    &*g_var**2-f_var**2*g_var**2+18._ki*f_var*c_var*h_var*g_va&
                    &r+18._ki*c_var*h_var*f_var**2+6._ki*f_var*c_var**3+3._ki*f_&
                    &var**3*g_var+6._ki*f_var**4-32._ki*g_var*f_var**2*c_var-19&
                    &._ki*g_var*f_var*c_var**2-2._ki*g_var*c_var**3)/f_var**2/g&
                    &_var**3
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*z*c_var**4/d_var**3*(-1._ki+z)/f_var**2*log(&
                    &c_var)+1._ki/12._ki*z*(-1._ki+z)*h_var**2*(f_var**4*c_var**&
                    &2+10._ki*d_var**3*c_var*f_var**2-8._ki*d_var**3*c_var**2*f&
                    &_var-2._ki*d_var*c_var*f_var**4+3._ki*f_var**2*d_var**4-6.&
                    &_ki*d_var**3*f_var**3+3._ki*f_var**4*d_var**2-8._ki*d_var**&
                    &4*c_var*f_var+10._ki*d_var**4*c_var**2)/d_var**3/f_var**6&
                    &*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/24._ki*z*(-1._k&
                    &i+z)*c_var**2*(2._ki*c_var-d_var)/f_var**2/d_var**2
                  !
                end select
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
                  fg=1._ki/12._ki*z*h_var**3*(-3._ki*f_var**3*h_var*g_var+3._ki*&
                    &f_var**4*h_var+f_var**3*g_var*c_var+3._ki*f_var**3*g_var*&
                    &*2-3._ki*f_var**4*g_var+3._ki*h_var*f_var**2*g_var**2-2._ki&
                    &*c_var*g_var**2*f_var**2-3._ki*f_var**2*g_var**3-3._ki*g_v&
                    &ar**3*f_var*h_var+3._ki*f_var*g_var**3*c_var+3._ki*g_var**&
                    &4*f_var-4._ki*g_var**4*c_var)/g_var**4/f_var**5*(log(z)+l&
                    &og(1._ki-z)+z_log(s23,1._ki))-1._ki/12._ki*z/g_var**4*log(e_&
                    &var)/f_var**3*e_var**2*(-2._ki*g_var*f_var*c_var**2+c_var&
                    &**2*g_var**2+3._ki*c_var**2*f_var**2+6._ki*c_var*f_var**3+&
                    &4._ki*g_var*f_var**2*c_var-2._ki*g_var**2*c_var*f_var+3._ki&
                    &*f_var**4+6._ki*f_var**3*g_var+3._ki*f_var**2*g_var**2)-1.&
                    &_ki/24._ki*z/f_var**3*(9._ki*f_var**4*g_var+6._ki*f_var**5+2&
                    &._ki*g_var**2*c_var**3-18._ki*f_var*g_var**3*c_var+18._ki*h&
                    &_var*f_var*g_var**2*c_var-4._ki*f_var*g_var*c_var**3-35._k&
                    &i*f_var**2*g_var*c_var**2+18._ki*h_var*f_var**3*c_var-40.&
                    &_ki*f_var**3*g_var*c_var+6._ki*f_var**2*c_var**3+2._ki*f_va&
                    &r**3*g_var**2-18._ki*c_var**2*g_var**2*f_var-54._ki*c_var*&
                    &g_var**2*f_var**2+36._ki*f_var**2*c_var*h_var*g_var)/g_va&
                    &r**3
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*z*c_var**4/d_var**2/f_var**3*log(c_var)+1._ki&
                    &/12._ki*z*h_var**3*(-3._ki*d_var*c_var*f_var-f_var**2*c_va&
                    &r+3._ki*d_var*f_var*h_var-3._ki*d_var**2*f_var+3._ki*f_var*&
                    &*2*d_var+4._ki*c_var*d_var**2)/f_var**5/d_var**2*(log(z)+&
                    &log(1._ki-z)+z_log(s23,1._ki))+1._ki/12._ki*z/f_var**3*c_var&
                    &**3/d_var
                  !
                end select
                !
              end select
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
                  fg=-1._ki/4._ki*(-1._ki+z)**3*h_var*(2._ki*f_var**7*h_var*c_va&
                    &r+4._ki*f_var**4*h_var*g_var**3*c_var+f_var**6*h_var*c_va&
                    &r**2+2._ki*g_var**4*f_var**3*h_var**2-2._ki*g_var*f_var**6&
                    &*h_var**2+2._ki*g_var**2*f_var**5*h_var**2-2._ki*g_var**3*&
                    &f_var**4*h_var**2-3._ki*f_var**3*h_var*g_var**5+3._ki*f_va&
                    &r**4*h_var*g_var**4-3._ki*f_var**5*h_var*g_var**3+3._ki*f_&
                    &var**6*h_var*g_var**2+f_var**7*h_var*g_var+f_var**8*h_va&
                    &r-4._ki*f_var**3*h_var*g_var**3*c_var**2+2._ki*f_var**4*h_&
                    &var*c_var**2*g_var**2-f_var**5*h_var*g_var*c_var**2-2._ki&
                    &*f_var**5*h_var*c_var*g_var**2+6._ki*g_var**4*f_var**2*h_&
                    &var*c_var**2+8._ki*g_var**5*f_var**2*h_var*c_var-10._ki*g_&
                    &var**5*f_var*h_var*c_var**2-6._ki*f_var**3*h_var*c_var*g_&
                    &var**4+g_var**3*f_var**5*c_var-g_var**2*f_var**5*c_var**&
                    &2+10._ki*g_var**6*c_var**2*f_var+10._ki*g_var**5*c_var**3*&
                    &f_var-4._ki*g_var**6*f_var**2*c_var+3._ki*f_var**3*g_var**&
                    &5*c_var-8._ki*g_var**4*f_var**2*c_var**3-6._ki*g_var**5*f_&
                    &var**2*c_var**2-g_var*f_var**7*c_var-2._ki*g_var**2*f_var&
                    &**4*c_var**3-2._ki*g_var**4*f_var**4*c_var+g_var*f_var**5&
                    &*c_var**3+g_var*f_var**6*c_var**2+2._ki*g_var**4*f_var**3&
                    &*c_var**2+4._ki*g_var**3*f_var**3*c_var**3-g_var**5*f_var&
                    &**4-g_var*f_var**8-20._ki*g_var**6*c_var**3-g_var**3*f_va&
                    &r**6-g_var**2*f_var**7+g_var**4*f_var**5+g_var**6*f_var*&
                    &*3)/f_var**7/g_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki&
                    &))+1._ki/4._ki*(-1._ki+z)**3/g_var**4*log(e_var)/f_var*e_va&
                    &r**4+1._ki/24._ki*(-1._ki+z)**3*(2._ki*g_var**2*f_var-3._ki*g&
                    &_var*f_var**2+6._ki*f_var**3+6._ki*c_var**3+18._ki*c_var*h_&
                    &var*f_var-24._ki*c_var*f_var*g_var-3._ki*c_var**2*g_var+2.&
                    &_ki*c_var*g_var**2)/g_var**3/f_var
                  !
                case(2)
                  !
                  fg=-1._ki/4._ki*c_var**4*(-1._ki+z)**3/f_var/d_var**4*log(c_v&
                    &ar)-1._ki/4._ki*(-1._ki+z)**3*h_var*(2._ki*d_var**3*f_var**4&
                    &*h_var**2+3._ki*d_var**3*f_var**5*h_var-2._ki*d_var**4*f_v&
                    &ar**3*h_var**2-6._ki*d_var**4*f_var**4*h_var+3._ki*f_var**&
                    &3*d_var**5*h_var-4._ki*c_var*d_var**3*f_var**4*h_var+10._k&
                    &i*d_var**5*c_var**2*f_var*h_var+4._ki*f_var**3*d_var**3*c&
                    &_var**2*h_var+2._ki*f_var**4*d_var**2*c_var**2*h_var-8._ki&
                    &*d_var**5*f_var**2*h_var*c_var-f_var**6*c_var**3+20._ki*d&
                    &_var**6*c_var**3-3._ki*d_var**4*f_var**5+3._ki*d_var**5*f_&
                    &var**4-d_var**6*f_var**3+d_var**3*f_var**6+f_var**5*d_va&
                    &r*c_var**2*h_var-2._ki*c_var*d_var**2*f_var**5*h_var-16._k&
                    &i*c_var**2*d_var**4*f_var**2*h_var+14._ki*c_var*d_var**4*&
                    &f_var**3*h_var+4._ki*d_var**6*f_var**2*c_var-c_var*d_var*&
                    &*2*f_var**6-50._ki*c_var**3*d_var**5*f_var-c_var*d_var**3&
                    &*f_var**5+9._ki*c_var*d_var**4*f_var**4-11._ki*c_var*d_var&
                    &**5*f_var**3+26._ki*c_var**2*d_var**5*f_var**2-10._ki*d_va&
                    &r**6*c_var**2*f_var-4._ki*f_var**3*c_var**3*d_var**3-2._ki&
                    &*f_var**4*c_var**3*d_var**2-f_var**5*c_var**3*d_var-18._k&
                    &i*f_var**3*d_var**4*c_var**2+f_var**5*d_var**2*c_var**2+&
                    &f_var**6*d_var*c_var**2+38._ki*f_var**2*d_var**4*c_var**3&
                    &)/f_var**7/d_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))&
                    &-1._ki/24._ki*(-1._ki+z)**3*c_var*(6._ki*c_var**2-3._ki*c_var&
                    &*d_var+2._ki*d_var**2)/d_var**3/f_var
                  !
                end select
                !
              case(4)
                !
                select case(flag)
                !
                case(1)
                  !
                  fg=-1._ki/12._ki*(-1._ki+z)**2*h_var**2*(3._ki*f_var**7-3._ki*f&
                    &_var**6*g_var+3._ki*f_var**5*g_var**2-3._ki*f_var**4*g_var&
                    &**3+3._ki*f_var**3*g_var**4+6._ki*f_var**6*c_var+2._ki*f_va&
                    &r**4*g_var**2*c_var-f_var**4*g_var*c_var**2-4._ki*f_var**&
                    &5*g_var*c_var-3._ki*g_var**5*f_var**2-10._ki*g_var**5*c_va&
                    &r**2+3._ki*c_var**2*f_var**5-2._ki*g_var**4*f_var**2*c_var&
                    &-2._ki*g_var**4*f_var*c_var**2+8._ki*g_var**5*f_var*c_var)&
                    &/f_var**6/g_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-&
                    &1._ki/12._ki*(-1._ki+z)**2*e_var**3*(-3._ki*f_var*e_var+c_va&
                    &r*g_var-3._ki*g_var*f_var)/f_var**2/g_var**4*log(e_var)+1&
                    &._ki/24._ki*(-1._ki+z)**2*(-18._ki*g_var**2*c_var*f_var+c_va&
                    &r**2*g_var**2-f_var**2*g_var**2+18._ki*f_var*c_var*h_var*&
                    &g_var+18._ki*c_var*h_var*f_var**2+6._ki*f_var*c_var**3+3._k&
                    &i*f_var**3*g_var+6._ki*f_var**4-32._ki*g_var*f_var**2*c_va&
                    &r-19._ki*g_var*f_var*c_var**2-2._ki*g_var*c_var**3)/f_var*&
                    &*2/g_var**3
                  !
                case(2)
                  !
                  fg=1._ki/12._ki*c_var**4*(-1._ki+z)**2/f_var**2/d_var**3*log(&
                    &c_var)-1._ki/12._ki*(-1._ki+z)**2*h_var**2*(f_var**4*c_var*&
                    &*2+10._ki*d_var**3*c_var*f_var**2-8._ki*d_var**3*c_var**2*&
                    &f_var-2._ki*d_var*c_var*f_var**4+3._ki*f_var**2*d_var**4-6&
                    &._ki*d_var**3*f_var**3+3._ki*f_var**4*d_var**2-8._ki*d_var*&
                    &*4*c_var*f_var+10._ki*d_var**4*c_var**2)/d_var**3/f_var**&
                    &6*(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/24._ki*(-1._ki&
                    &+z)**2*c_var**2*(2._ki*c_var-d_var)/f_var**2/d_var**2
                  !
                end select
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
                  fg=-1._ki/12._ki*(-1._ki+z)*h_var**3*(-3._ki*f_var**3*h_var*g_&
                    &var+3._ki*f_var**4*h_var+f_var**3*g_var*c_var+3._ki*f_var*&
                    &*3*g_var**2-3._ki*f_var**4*g_var+3._ki*h_var*f_var**2*g_va&
                    &r**2-2._ki*c_var*g_var**2*f_var**2-3._ki*f_var**2*g_var**3&
                    &-3._ki*g_var**3*f_var*h_var+3._ki*f_var*g_var**3*c_var+3._k&
                    &i*g_var**4*f_var-4._ki*g_var**4*c_var)/g_var**4/f_var**5*&
                    &(log(z)+log(1._ki-z)+z_log(s23,1._ki))+1._ki/12._ki*(-1._ki+z&
                    &)/g_var**4*log(e_var)/f_var**3*e_var**2*(-2._ki*g_var*f_v&
                    &ar*c_var**2+c_var**2*g_var**2+3._ki*c_var**2*f_var**2+6._k&
                    &i*c_var*f_var**3+4._ki*g_var*f_var**2*c_var-2._ki*g_var**2&
                    &*c_var*f_var+3._ki*f_var**4+6._ki*f_var**3*g_var+3._ki*f_va&
                    &r**2*g_var**2)+1._ki/24._ki*(-1._ki+z)/f_var**3*(9._ki*f_var&
                    &**4*g_var+6._ki*f_var**5+2._ki*g_var**2*c_var**3-18._ki*f_v&
                    &ar*g_var**3*c_var+18._ki*h_var*f_var*g_var**2*c_var-4._ki*&
                    &f_var*g_var*c_var**3-35._ki*f_var**2*g_var*c_var**2+18._ki&
                    &*h_var*f_var**3*c_var-40._ki*f_var**3*g_var*c_var+6._ki*f_&
                    &var**2*c_var**3+2._ki*f_var**3*g_var**2-18._ki*c_var**2*g_&
                    &var**2*f_var-54._ki*c_var*g_var**2*f_var**2+36._ki*f_var**&
                    &2*c_var*h_var*g_var)/g_var**3
                  !
                case(2)
                  !
                  fg=-1._ki/12._ki*c_var**4/d_var**2*(-1._ki+z)/f_var**3*log(c_&
                    &var)-1._ki/12._ki*(-1._ki+z)*h_var**3*(-3._ki*d_var*c_var*f_&
                    &var-f_var**2*c_var+3._ki*d_var*f_var*h_var-3._ki*d_var**2*&
                    &f_var+3._ki*f_var**2*d_var+4._ki*c_var*d_var**2)/f_var**5/&
                    &d_var**2*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1._ki/12._ki&
                    &*(-1._ki+z)/f_var**3*c_var**3/d_var
                  !
                end select
                !
              end select
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
                  fg=-1._ki/4._ki*h_var**4*(-g_var+f_var)*(f_var**2+g_var**2)/&
                    &f_var**4/g_var**4*(log(z)+log(1._ki-z)+z_log(s23,1._ki))-1&
                    &._ki/4._ki*e_var*(2._ki*g_var*f_var**2*e_var+c_var**2*g_var&
                    &**2+f_var**2*g_var**2+c_var**2*f_var**2+2._ki*c_var*f_var&
                    &**3+f_var**4)*(-f_var*e_var+c_var*g_var-g_var*f_var)/f_v&
                    &ar**4/g_var**4*log(e_var)+1._ki/24._ki*(15._ki*f_var**4*g_v&
                    &ar-18._ki*g_var**4*c_var+6._ki*f_var**5+6._ki*g_var**2*c_va&
                    &r**3-18._ki*g_var**3*c_var**2-72._ki*f_var*g_var**3*c_var+&
                    &54._ki*h_var*f_var*g_var**2*c_var+18._ki*g_var**3*c_var*h_&
                    &var-6._ki*f_var*g_var*c_var**3-51._ki*f_var**2*g_var*c_var&
                    &**2+18._ki*h_var*f_var**3*c_var-48._ki*f_var**3*g_var*c_va&
                    &r+6._ki*f_var**2*c_var**3+11._ki*f_var**3*g_var**2-57._ki*c&
                    &_var**2*g_var**2*f_var-106._ki*c_var*g_var**2*f_var**2+54&
                    &._ki*f_var**2*c_var*h_var*g_var)/f_var**3/g_var**3
                  !
                case(2)
                  !
                  fg=1._ki/4._ki*c_var**4/d_var/f_var**4*log(c_var)-1._ki/4._ki*&
                    &(log(z)+log(1._ki-z)+z_log(s23,1._ki))*h_var**4/f_var**4/d&
                    &_var
                  !
                end select
                !
              end select
              !
            end select
            !
          end select
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
          fg=-1._ki/6._ki*h_var**2*(f_var*h_var*g_var-f_var**2*h_var-g&
            &_var*f_var*c_var-2._ki*f_var*g_var**2+2._ki*g_var*f_var**2&
            &+2._ki*c_var*g_var**2)/f_var**3/g_var**2*(log(z)+log(1._ki&
            &-z)+z_log(-s23,-1._ki))-1._ki/6._ki/g_var**2*log(-e_var)/f_&
            &var*e_var**3-1._ki/18._ki*(6._ki*f_var*c_var+3._ki*c_var**2-&
            &5._ki*f_var*g_var+3._ki*f_var**2)/f_var/g_var
          !
          case(2)
          !
          fg=1._ki/6._ki*c_var**3/d_var**2/f_var*log(-c_var)-1._ki/6._ki&
            &*h_var**2*(c_var*f_var*d_var+f_var**2*c_var-d_var*f_var*&
            &h_var+2._ki*d_var**2*f_var-2._ki*d_var*f_var**2-2._ki*c_var&
            &*d_var**2)/f_var**3/d_var**2*(log(z)+log(1._ki-z)+z_log(-&
            &s23,-1._ki))+1._ki/6._ki*c_var**2/f_var/d_var
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
              fg=1._ki/24._ki*h_var**3*(-2._ki*f_var*h_var*g_var+f_var**2*h&
                &_var+2._ki*c_var*g_var*f_var+4._ki*f_var*g_var**2-2._ki*f_v&
                &ar**2*g_var-6._ki*c_var*g_var**2)/f_var**4/g_var**2*(log(&
                &z)+log(1._ki-z)+z_log(-s23,-1._ki))-1._ki/24._ki/g_var**2*lo&
                &g(-e_var)/f_var**2*e_var**4-1._ki/144._ki/f_var**2*(-13._ki&
                &*f_var**2*g_var+6._ki*f_var**3+6._ki*c_var**3+18._ki*c_var*&
                &h_var*f_var-18._ki*c_var*g_var*f_var)/g_var
              !
            case(2)
              !
              fg=1._ki/24._ki*c_var**3*(c_var*d_var+2._ki*c_var*f_var+4._ki*&
                &d_var*f_var)/d_var**3/f_var**2*log(-c_var)+1._ki/24._ki*h_&
                &var**3*(-2._ki*f_var**2*c_var*d_var-2._ki*c_var*f_var**3+d&
                &_var*f_var**2*h_var+2._ki*d_var**2*f_var**2+2._ki*f_var**3&
                &*d_var-2._ki*c_var*d_var**2*f_var+2._ki*f_var*d_var**2*h_v&
                &ar-4._ki*d_var**3*f_var+6._ki*c_var*d_var**3)/f_var**4/d_v&
                &ar**3*(log(z)+log(1._ki-z)+z_log(-s23,-1._ki))+1._ki/24._ki*&
                &c_var**2*(2._ki*c_var*f_var+c_var*d_var+3._ki*d_var*f_var)&
                &/d_var**2/f_var**2
              !
            end select
            !
          case(2)
            !
            select case(flag)
            !
            case(1)
              !
              fg=-1._ki/12._ki*z*h_var**2*(-2._ki*f_var**5*g_var+3._ki*f_var&
                &**4*g_var**2-3._ki*f_var**3*g_var**3+f_var**6+3._ki*f_var*&
                &*2*g_var**4+6._ki*g_var**4*c_var**2+2._ki*f_var**2*c_var*g&
                &_var**3-2._ki*f_var**4*c_var*g_var-6._ki*f_var*g_var**4*c_&
                &var+c_var**2*f_var**4+2._ki*c_var*f_var**5)/f_var**5/g_va&
                &r**3*(log(z)+log(1._ki-z)+z_log(-s23,-1._ki))+1._ki/12._ki*z&
                &/g_var**3*log(-e_var)/f_var*e_var**4+1._ki/72._ki*z*(7._ki*&
                &f_var*g_var**2-3._ki*f_var**2*g_var+6._ki*f_var**3+6._ki*c_&
                &var**3+18._ki*c_var*h_var*f_var-24._ki*c_var*g_var*f_var-3&
                &._ki*c_var**2*g_var)/f_var/g_var**2
              !
            case(2)
              !
              fg=-1._ki/12._ki*z*c_var**4/d_var**3/f_var*log(-c_var)+1._ki/&
                &12._ki*z*h_var**2*(c_var**2*f_var**4+8._ki*c_var*f_var**2*&
                &d_var**3-6._ki*c_var**2*d_var**3*f_var-2._ki*c_var*f_var**&
                &4*d_var+3._ki*f_var**2*d_var**4-6._ki*f_var**3*d_var**3+3.&
                &_ki*f_var**4*d_var**2-6._ki*d_var**4*c_var*f_var+6._ki*c_va&
                &r**2*d_var**4)/f_var**5/d_var**3*(log(z)+log(1._ki-z)+z_l&
                &og(-s23,-1._ki))-1._ki/24._ki*z*c_var**2*(2._ki*c_var-d_var)&
                &/f_var/d_var**2
              !
            end select
            !
          case(3)
            !
            select case(flag)
            !
            case(1)
              !
              fg=1._ki/12._ki*(z-1)*h_var**2*(-2._ki*f_var**5*g_var+3._ki*f_&
                &var**4*g_var**2-3._ki*f_var**3*g_var**3+f_var**6+3._ki*f_v&
                &ar**2*g_var**4+6._ki*g_var**4*c_var**2+2._ki*f_var**2*c_va&
                &r*g_var**3-2._ki*f_var**4*c_var*g_var-6._ki*f_var*g_var**4&
                &*c_var+c_var**2*f_var**4+2._ki*c_var*f_var**5)/f_var**5/g&
                &_var**3*(log(z)+log(1._ki-z)+z_log(-s23,-1._ki))-1._ki/12._k&
                &i*(z-1)/g_var**3*log(-e_var)/f_var*e_var**4-1._ki/72._ki*(&
                &z-1)*(7._ki*f_var*g_var**2-3._ki*f_var**2*g_var+6._ki*f_var&
                &**3+6._ki*c_var**3+18._ki*c_var*h_var*f_var-24._ki*c_var*g_&
                &var*f_var-3._ki*c_var**2*g_var)/f_var/g_var**2
              !
            case(2)
              !
              fg=1._ki/12._ki*c_var**4/d_var**3*(z-1)/f_var*log(-c_var)-1.&
                &_ki/12._ki*(z-1)*h_var**2*(c_var**2*f_var**4+8._ki*c_var*f_&
                &var**2*d_var**3-6._ki*c_var**2*d_var**3*f_var-2._ki*c_var*&
                &f_var**4*d_var+3._ki*f_var**2*d_var**4-6._ki*f_var**3*d_va&
                &r**3+3._ki*f_var**4*d_var**2-6._ki*d_var**4*c_var*f_var+6.&
                &_ki*c_var**2*d_var**4)/f_var**5/d_var**3*(log(z)+log(1._ki&
                &-z)+z_log(-s23,-1._ki))+1._ki/24._ki*(z-1)*c_var**2*(2._ki*c&
                &_var-d_var)/f_var/d_var**2
              !
            end select
            !
          case(4)
            !
            select case(flag)
            !
            case(1)
              !
              fg=1._ki/24._ki*h_var**3*(-2._ki*f_var**2*h_var*g_var+2._ki*f_&
                &var**3*h_var+c_var*f_var**2*g_var+3._ki*g_var**2*f_var**2&
                &-3._ki*g_var*f_var**3+2._ki*f_var*g_var**2*h_var-2._ki*c_va&
                &r*g_var**2*f_var-3._ki*g_var**3*f_var+3._ki*g_var**3*c_var&
                &)/g_var**3/f_var**4*(log(z)+log(1._ki-z)+z_log(-s23,-1._ki&
                &))+1._ki/24._ki*e_var**3*(-2._ki*f_var*e_var+c_var*g_var-3.&
                &_ki*f_var*g_var)/g_var**3/f_var**2*log(-e_var)-1._ki/144._k&
                &i/f_var**2*(-13._ki*g_var**2*f_var**2-36._ki*c_var*g_var**&
                &2*f_var+36._ki*f_var*c_var*h_var*g_var+36._ki*c_var*h_var*&
                &f_var**2+12._ki*f_var*c_var**3+12._ki*g_var*f_var**3+12._ki&
                &*f_var**4-36._ki*c_var**2*g_var*f_var-6._ki*c_var**3*g_var&
                &-54._ki*c_var*f_var**2*g_var)/g_var**2
              !
            case(2)
              !
              fg=-1._ki/24._ki*c_var**4/d_var**2/f_var**2*log(-c_var)-1._ki&
                &/24._ki*h_var**3*(-2._ki*c_var*f_var*d_var-c_var*f_var**2+&
                &2._ki*f_var*d_var*h_var-3._ki*f_var*d_var**2+3._ki*f_var**2&
                &*d_var+3._ki*d_var**2*c_var)/f_var**4/d_var**2*(log(z)+lo&
                &g(1._ki-z)+z_log(-s23,-1._ki))-1._ki/24._ki/f_var**2*c_var**&
                &3/d_var
              !
            end select
            !
          end select
          !
       end if
        !
      end if
    end function fg
    !
end module function_4p3m
!
