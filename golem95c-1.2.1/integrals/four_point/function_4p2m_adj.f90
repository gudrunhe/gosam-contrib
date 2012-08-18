! 
!****h* src/integrals/four_point/function_4p2m_adj
! NAME
!
!  Module function_4p2m_adj
!
! USAGE
!
!  use function_4p2m_adj
!
! DESCRIPTION
!
!  This module computes the six-dimensional and eight dimensional 
!  two adjacent mass four point function with or without Feynman parameters
!  in the numerator.
!
! OUTPUT
!
!  This module exports three functions f4p2m_adj, f4p2m_adj_c and f2a
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
module function_4p2m_adj
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
  real(ki) :: s23_glob,s24_glob,s34_glob,s13_glob
  real(ki) :: eps_glob
  integer :: par1_glob,par2_glob,par3_glob,par4_glob
  character (len=3) :: dim_glob
  !
  real(ki), dimension(4) :: b
  real(ki) :: sumb
  real(ki), dimension(4,4) :: invs,s_mat
  integer, dimension(4) :: par
  integer, dimension(4) :: s = (/1,2,3,4/)
  real(ki) :: lamb
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
 public :: f4p2m_adj,f2a,f4p2m_adj_c
  !
  contains
    !
    !****f* src/integrals/four_point/function_4p2m_adj/f4p2m_adj
    ! NAME
    !
    !  Function f4p2m_adj
    !
    ! USAGE
    !
    !  real_dim_4 = f4p2m_adj(dim,s24,s13,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function computes the six dimensional/eight dimensional
    !  two adjacent mass four point function with or without Feynman parameters 
    !  in the numerator.
    !
    ! INPUTS
    !
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           two adjacent mass four point function, dim="n+4" eight dimensional
    !           two adjacent mass four point function
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
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
    !  * a six dimensional two adjacent mass four point function 
    !    with no Feynman parameters in the numerator:
    !    real_dim_4 = f4p2m_adj("n+2",s24,s13,s23,s34,0,0,0,0)
    !  * a eight dimensional two adjacent mass four point function 
    !    with no Feynman parameters in the numerator:
    !    real_dim_4 = f4p2m_adj("n+4",s24,s13,s23,s34,0,0,0,0)
    !  * a six dimensional two adjacent mass four point function 
    !    with the Feynman parameter z1 in the numerator:
    !    real_dim_4 = f4p2m_adj("n+2",s24,s13,s23,s34,0,0,0,1)
    !  * a six dimensional two adjacent mass four point function 
    !    with the Feynman parameters z1^2*z2 in the numerator:
    !    real_dim_4 = f4p2m_adj("n+2",s24,s13,s23,s34,0,2,1,1)
    !
    !*****
    function f4p2m_adj(dim,s24,s13,s23,s34,par1,par2,par3,par4)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: s24,s13,s23,s34
      integer, intent (in) :: par1,par2,par3,par4 
      real(ki), dimension(4) :: f4p2m_adj
      !
      integer :: nb_par
      real(ki) :: plus_grand
      real(ki) :: norma
      complex(ki) :: resto,abserro
      !
      par = (/par1,par2,par3,par4/)
      !
      s_mat(1,:) = (/0._ki,0._ki,s13,0._ki/)
      s_mat(2,:) = (/0._ki,0._ki,s23,s24/)
      s_mat(3,:) = (/s13,s23,0._ki,s34/)
      s_mat(4,:) = (/0._ki,s24,s34,0._ki/)
      ! on redefinit la matrice S de telle facon a ce que ses elements
      ! soient entre -1 et 1
      plus_grand = maxval(array=abs(s_mat))
      s_mat = s_mat/plus_grand
      !
      b(1) = (s_mat(1,3)*s_mat(2,4)-s_mat(1,3)*s_mat(2,3)&
            &-s_mat(1,3)*s_mat(3,4)+2._ki*s_mat(2,3)*s_mat(3,4))&
            &/(s_mat(1,3)**2*s_mat(2,4))
      b(2) = (s_mat(1,3)-s_mat(3,4))/(s_mat(1,3)*s_mat(2,4))
      b(3) = 1._ki/s_mat(1,3)
      b(4) = (s_mat(1,3)-s_mat(2,3))/(s_mat(1,3)*s_mat(2,4))
      !
      sumb = 2._ki*(s_mat(1,3)**2-s_mat(1,3)*s_mat(2,3)&
            &-s_mat(1,3)*s_mat(3,4)+s_mat(1,3)*s_mat(2,4)&
            &+s_mat(2,3)*s_mat(3,4))/(s_mat(1,3)**2*s_mat(2,4))
      !
      invs(1,1) = 2._ki*s_mat(2,3)/s_mat(2,4)*s_mat(3,4)/s_mat(1,3)**2
      invs(1,2) = -1._ki/s_mat(1,3)/s_mat(2,4)*s_mat(3,4)
      invs(1,3) = 1._ki/s_mat(1,3)
      invs(1,4) = -1._ki/s_mat(1,3)*s_mat(2,3)/s_mat(2,4)
      invs(2,1) = -1._ki/s_mat(1,3)/s_mat(2,4)*s_mat(3,4)
      invs(2,2) = 0._ki
      invs(2,3) = 0._ki
      invs(2,4) = 1._ki/s_mat(2,4)
      invs(3,1) = 1._ki/s_mat(1,3)
      invs(3,2) = 0._ki
      invs(3,3) = 0._ki
      invs(3,4) = 0._ki
      invs(4,1) = -1._ki/s_mat(1,3)*s_mat(2,3)/s_mat(2,4)
      invs(4,2) = 1._ki/s_mat(2,4)
      invs(4,3) = 0._ki
      invs(4,4) = 0._ki
      !
      lamb = s_mat(1,3)**2-s_mat(1,3)*s_mat(2,3)&
            &-s_mat(1,3)*s_mat(3,4)+s_mat(1,3)*s_mat(2,4)&
            &+s_mat(2,3)*s_mat(3,4)
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
      f4p2m_adj = 0._ki
      !
      if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_4p2m_adj) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f4p2m_adj (in file f4p2m_adj.f90): &
        &the flag rat to compute the rational part is on &
        &and the program reachs a region of phase space in &
        &which det(G) = 0  Becareful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to  go on, he has to &
        &reduce the value of the parameter coupure_4p2m_adj'
        call catch_exception(0)
      end if
      !
      if (abs(sumb) > coupure_4p2m_adj) then
        !
        ! analytic computation
        !
        if (dim == "n+2") then
          !
          f4p2m_adj(3:4)= a4p2m_adj_np2(s_mat(2,4),s_mat(1,3),&
                            &s_mat(2,3),s_mat(3,4),&
                            &par1,par2,par3,par4)/plus_grand
          !
        else if (dim == "n+4") then
          !
          f4p2m_adj = a4p2m_adj_np4(s_mat(2,4),s_mat(1,3),&
                      &s_mat(2,3),s_mat(3,4),&
                      &par1,par2,par3,par4)
          f4p2m_adj(3) = f4p2m_adj(3)-log(plus_grand)*norma
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
        s13_glob = s_mat(1,3)
        s23_glob = s_mat(2,3)
        s24_glob = s_mat(2,4)
        s34_glob = s_mat(3,4)
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
        eps_glob = sign(1._ki,s34_glob-s24_glob)
        !
        origine_info_par = "f4p2m_adj, dimension "//dim
        num_grand_b_info_par = lamb
        denom_grand_b_info_par = (s_mat(1,3)**2*s_mat(2,4))
        !
        call generic_eval_numer(eval_numer_gi,0._ki,1._ki,tolerance,resto,abserro)
        !
        if (dim == "n+2") then      
          resto = resto/plus_grand
        else if (dim == "n+4") then
          f4p2m_adj(1) = norma
          f4p2m_adj(2) = 0._ki
          resto = resto-log(plus_grand/mu2_scale_par)*norma
        end if
        !
        f4p2m_adj(3) = real(resto,ki)
        f4p2m_adj(4) = aimag(resto)
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
    end function f4p2m_adj
    !
    !****f* src/integrals/four_point/function_4p2m_adj/f4p2m_adj_c
    ! NAME
    !
    !  Function f4p2m_adj_c
    !
    ! USAGE
    !
    !  complex_dim_2 = f4p2m_adj_c(dim,s24,s13,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function computes the same thing that the function f4p2m_adj
    !
    ! INPUTS
    !
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           two adjacent mass four point function, dim="n+4" eight dimensional
    !           two adjacent mass four point function
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
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
    !  see function f4p2m_adj
    !
    !*****
    function f4p2m_adj_c(dim,s24,s13,s23,s34,par1,par2,par3,par4)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: s24,s13,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      complex(ki), dimension(2) :: f4p2m_adj_c
      !
      real(ki), dimension(4) :: res4
      !
      res4 = f4p2m_adj(dim,s24,s13,s23,s34,par1,par2,par3,par4)
      call to_complex(res4,f4p2m_adj_c)
      !
    end function f4p2m_adj_c
    !
    !****if* src/integrals/four_point/function_4p2m_adj/a4p2m_adj_np2
    ! NAME
    !
    !  recursive function a4p2m_adj_np2
    !
    ! USAGE
    !
    !  real_dim_2 = a4p2m_adj_np2(s24,s13,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function is the core for the analytic computation of the six dimensional
    !  two adjacent mass four point function. It is recursive and implement the formulae
    !  of JHEP 10 (2005) 015.
    !
    !
    ! INPUTS
    !
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
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
    recursive function a4p2m_adj_np2(s24,s13,s23,s34,par1,par2,par3,par4) result(res_4p2m_adj_np2)
      !
      real(ki), intent (in) :: s24,s13,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(2) :: res_4p2m_adj_np2
      !
      integer, dimension(3) :: smj,sm1
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
        !~ sm1 = s .minus. (/1/)
        sm1 = unpackb(ibclr(b_pro,1),3)
        !
        if (deja_calcule3(1,1)) then
          !
          truc1 = resultat3(1,1,:)
          !
        else
          !
          truc1 = f3p_sc(s_mat,sm1)
          resultat3(1,1,:) = truc1
          deja_calcule3(1,1) = .true.
          !
        end if
        !
        ctemp = f2a(s24,s13,s23,s34)
        res_4p2m_adj_np2(1) = -s13*( 2._ki*real(ctemp,ki) &
                          + s13*s24*b(1)*truc1(5) )/(2._ki*lamb)
        res_4p2m_adj_np2(2) = -s13*( 2._ki*aimag(ctemp) &
                          + s13*s24*b(1)*truc1(6) )/(2._ki*lamb)
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
          temp0 = a4p2m_adj_np2(s24,s13,s23,s34,0,0,0,0)
          resultat(1,:) = temp0
          deja_calcule(1) = .true.
          !
        end if
        !
        temp0 = b(par4)*temp0
        !
        temp1 = 0._ki
        temp2 = 0._ki
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
        res_4p2m_adj_np2(1) = (temp0(1) + temp1(5) + temp2(5))/sumb
        res_4p2m_adj_np2(2) = (temp0(2) + temp1(6) + temp2(6))/sumb
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
          temp10 = a4p2m_adj_np2(s24,s13,s23,s34,0,0,0,par4)
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
          temp11 = a4p2m_adj_np2(s24,s13,s23,s34,0,0,0,par3)
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
        res_4p2m_adj_np2(1) = (temp0(1) + temp1(5) + temp2(5) + temp3(5)) &
                          *2._ki/3._ki/sumb
        res_4p2m_adj_np2(2) = (temp0(2) + temp1(6) + temp2(6) + temp3(6)) &
                          *2._ki/3._ki/sumb
      !
      ! cas avec trois parametres de feynman au numerateur
      !
      else
        !
        temp10 = a4p2m_adj_np2(s24,s13,s23,s34,0,0,par2,par3)
        temp11 = a4p2m_adj_np2(s24,s13,s23,s34,0,0,par2,par4)
        temp12 = a4p2m_adj_np2(s24,s13,s23,s34,0,0,par3,par4)
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
        res_4p2m_adj_np2(1) = ( temp0(1) + temp1(5) + temp2(5) + temp3(5) &
                       + temp4(5) )/2._ki/sumb
        res_4p2m_adj_np2(2) = ( temp0(2) + temp1(6) + temp2(6) + temp3(6) &
                       + temp4(6) )/2._ki/sumb
        !
      end if
      !
    end function a4p2m_adj_np2
    !
    !****if* src/integrals/four_point/function_4p2m_adj/a4p2m_adj_np4
    ! NAME
    !
    !  recursive function a4p2m_adj_np4
    !
    ! USAGE
    !
    !  real_dim_4 = a4p2m_adj_np4(s24,s13,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function is the core for the analytic computation of the eight dimensional
    !  two adjacent mass four point function. It is recursive and implement the formulae
    !  of JHEP 10 (2005) 015.
    !
    !
    ! INPUTS
    !
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
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
    recursive function a4p2m_adj_np4(s24,s13,s23,s34,par1,par2,par3,par4) result(res_4p2m_adj_np4)
      !
      real(ki), intent (in) :: s24,s13,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(4) :: res_4p2m_adj_np4
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
          temp0 = a4p2m_adj_np2(s24,s13,s23,s34,0,0,0,0)
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
        res_4p2m_adj_np4(1) = (-temp1(1))/(3._ki*sumb)
        res_4p2m_adj_np4(2) = (-temp1(2))/(3._ki*sumb)
        res_4p2m_adj_np4(3) = (temp0(1)-temp1(3)-2._ki/3._ki*temp1(1))/(3._ki*sumb)
        res_4p2m_adj_np4(4) = (temp0(2)-temp1(4)-2._ki/3._ki*temp1(2))/(3._ki*sumb)
      !
      ! cas avec un parametre de feynman au numerateur
      !
      else if (nb_par_loc == 1) then
        !
        temp0 = a4p2m_adj_np2(s24,s13,s23,s34,0,0,0,par4)/3._ki
        temp1 = b(par4)*a4p2m_adj_np4(s24,s13,s23,s34,0,0,0,0)
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
        res_4p2m_adj_np4(1) = ( temp1(1)+temp2(1)+temp3(1) )/(2._ki*sumb)
        res_4p2m_adj_np4(2) = ( temp1(2)+temp2(2)+temp3(2) )/(2._ki*sumb)
        res_4p2m_adj_np4(3) = ( temp1(3)+temp1(1)/6._ki+temp2(3)+temp2(1)/2._ki &
                           +temp3(3)+temp3(1)/2._ki+temp0(1) )/(2._ki*sumb)
        res_4p2m_adj_np4(4) = ( temp1(4)+temp1(2)/6._ki+temp2(4)+temp2(2)/2._ki &
                           +temp3(4)+temp3(2)/2._ki+temp0(2) )/(2._ki*sumb)
      !
      ! cas avec plus de un parametre de feynman au numerateur
      !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a4p2m_adj_np4:'
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
    end function a4p2m_adj_np4
    !
    !****f* src/integrals/four_point/function_4p2m_adj/f2a
    ! NAME
    !
    !  function f2a
    !
    ! USAGE
    !
    !  complex = f2a(u,v,w,x)
    !
    ! DESCRIPTION
    !
    !  This function is the "finite part" of the scalar four dimensional two  
    !  adjacent mass four point function. The expression has been taken in 
    !  Nucl. Phys. {\bf B615} (2001) , 385
    !
    !
    ! INPUTS
    !
    !  * u -- a real (type ki), (p1+p2)^2
    !  * v -- a real (type ki), (p2+p3)^2
    !  * w -- a real (type ki), p3^2
    !  * x -- a real (type ki), p4^2
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
    function f2a(u,v,w,x)
      !
      real(ki), intent(in) :: u,v,w,x
      complex(ki) :: f2a
      !
      f2a =   zdilog(1._ki-w/v,sign(un,w-v)) &
            + zdilog(1._ki-x/v,sign(un,x-v)) &
            + z_log(u/v,sign(un,v-u))*z_log(x/v,sign(un,v-x))/2._ki &
            + z_log(w/v,sign(un,v-w))*z_log(u/x,sign(un,x-u))/2._ki
      !
    end function f2a
    !
    !****if* src/integrals/four_point/function_4p2m_adj/eval_numer_gi
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
    !  the routine adapt_gauss1 in the function f4p2m_adj
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
      eval_numer_gi = fg(z,s24_glob,s13_glob,s23_glob,s34_glob,&
                      &  par1_glob,par2_glob,par3_glob,par4_glob,&
                      &  dim_glob)
      eval_numer_gi = eval_numer_gi*jacob
      !
    end function eval_numer_gi
    !
    !****if* src/integrals/four_point/function_4p2m_adj/fg
    ! NAME
    !
    !  function fg
    !
    ! USAGE
    !
    !  complex = fg(z,s24,s13,s23,s34,par1,par2,par3,par4,dim)
    !
    ! DESCRIPTION
    !
    !  This function contains the one dimensional integral representation of 
    !  the six/eight dimensional two adjacent mass four point function
    !
    !
    ! INPUTS
    !
    !  * z -- a real (type ki), integration variable
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s23 -- a real (type ki), the S matrix element 2,3
    !  * s34 -- a real (type ki), the S matrix element 3,4
    !  * par1 -- an integer, the label of the fourth Feynman parameter, if none, put 0
    !  * par2 -- an integer, the label of the third Feynman parameter, if none, put 0
    !  * par3 -- an integer, the label of the second Feynman parameter, if none, put 0
    !  * par4 -- an integer, the label of the first Feynman parameter, if none, put 0
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           two adjacent mass four point function, dim="n+4" eight dimensional
    !           two adjacent mass four point function
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
    !  two adjacent mass four point function
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function fg(z,s24,s13,s23,s34,par1,par2,par3,par4,dim)
      !
      complex(ki), intent (in) :: z
      real(ki), intent (in) :: s24,s13,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      character (len=3) :: dim
      complex(ki) :: fg
      !
      integer, dimension(4) :: par
      integer :: nb_par
      complex(ki) :: c_var,d_var,e_var,f_var,g_var,h_var
      complex(ki) :: umz
      !
      par = (/par1,par2,par3,par4/)
      nb_par = count(mask=par/=0)
      umz = 1._ki - z
      !
      c_var = s13
      !
      f_var = z*s24+(1._ki-z)*(s34-s13)
      !
      g_var = z*(1._ki-z)*s23-z*s24-(1._ki-z)*s34
      !
      d_var = z*s23-s13
      !
      e_var = z*s24+(1._ki-z)*s34
      !
      h_var = z*s23
      !
      if (dim == "n+2") then      
        if (nb_par == 0) then
          !
          fg=log(e_var)/g_var/f_var*e_var-c_var/d_var/f_var*(log(1._ki-z)+z_&
            &log(s13,1._ki))-1._ki/g_var*(log(z)+log(1._ki-z)+z_log(s23,1._ki))/&
            &d_var*h_var
          !
        else if (nb_par == 1) then
          !
          select case(par4)
          !
          case(1)
            !
            fg=-1._ki/2._ki*(2._ki*d_var*f_var+c_var*d_var-c_var*d_var*z+c_var*f&
              &_var)*c_var/f_var**2/d_var**2*(log(1._ki-z)+z_log(s13,1._ki))+1._k&
              &i/2._ki/g_var*log(e_var)/f_var**2*e_var**2-1._ki/2._ki*c_var/f_var&
              &/d_var-1._ki/2._ki/g_var*(log(z)+log(1._ki-z)+z_log(s23,1._ki))/d_v&
              &ar**2*h_var**2
            !
          case(2)
            !
            fg=1._ki/2._ki*z/g_var**2*(log(z)+log(1._ki-z)+z_log(s23,1._ki))/d_va&
              &r**2*h_var*((1._ki-z)*d_var*h_var+c_var*g_var-d_var*g_var)-1._ki/&
              &2._ki/g_var**2*log(e_var)*z/f_var*e_var**2+1._ki/2._ki*z*c_var**2/&
              &d_var**2/f_var*(log(1._ki-z)+z_log(s13,1._ki))-1._ki/2._ki/g_var*z/&
              &d_var*c_var-1._ki/2._ki*z/g_var
            !
          case(3)
            !
            fg=-1._ki/2._ki*(-1._ki+z)/g_var**2*(log(z)+log(1._ki-z)+z_log(s23,1.&
              &_ki))/d_var**2*h_var*((1._ki-z)*d_var*h_var+c_var*g_var-d_var*g_v&
              &ar)+1._ki/2._ki/g_var**2*log(e_var)/f_var*e_var**2*(-1._ki+z)-1._ki&
              &/2._ki*c_var**2*(-1._ki+z)/d_var**2/f_var*(log(1._ki-z)+z_log(s13,&
              &1._ki))-1._ki/2._ki/g_var/d_var*c_var+1._ki/2._ki/g_var*z/d_var*c_va&
              &r+1._ki/2._ki*z/g_var-1._ki/2._ki/g_var
            !
          case(4)
            !
            fg=1._ki/2._ki*(-1._ki+z)/g_var**2*(log(z)+log(1._ki-z)+z_log(s23,1._k&
              &i))/d_var*h_var**2+1._ki/2._ki/g_var**2*log(e_var)/f_var**2*e_var&
              &*(f_var*e_var+g_var*f_var-c_var*g_var+c_var*g_var*z)-1._ki/2._ki*&
              &(-1._ki+z)/g_var/f_var*h_var-1._ki/2._ki*(-1._ki+z)*c_var**2/d_var/&
              &f_var**2*(log(1._ki-z)+z_log(s13,1._ki))-1._ki/2._ki/f_var
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
              fg=-1._ki/3._ki*c_var/f_var**3/d_var**3*(log(umz)+z_log(s13,1._ki))*&
                &(c_var**2*d_var**2*umz**2+3._ki*d_var*f_var**2*h_var+c_var**2*f_&
                &var**2+c_var**2*d_var*f_var*umz+3._ki*c_var*d_var**2*f_var*umz)+&
                &1._ki/6._ki*c_var/d_var**2/f_var**2*(-5._ki*d_var**2*umz-4._ki*c_va&
                &r*d_var*umz+5._ki*d_var*g_var+2._ki*c_var*g_var)-1._ki/3._ki/g_var*&
                &(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d_var**3*h_var**3+1._ki&
                &/3._ki/g_var*log(e_var)/f_var**3*e_var**3
              !
            case(2)
              !
              fg=1._ki/6._ki*c_var**2/f_var**2/d_var**3*(log(umz)+z_log(s13,1._ki)&
                &)*(1._ki-umz)*(3._ki*d_var*f_var+c_var*d_var*umz+2._ki*c_var*f_var&
                &)+1._ki/6._ki*h_var**2/g_var**2/d_var**3*(log(1._ki-umz)+log(umz)+&
                &z_log(s23,1._ki))*(1._ki-umz)*(-d_var*g_var+d_var*h_var*umz+2._ki*&
                &c_var*g_var)+1._ki/6._ki/g_var/d_var**2/f_var*(1._ki-umz)*(d_var**&
                &2*g_var-c_var**2*d_var*umz-2._ki*c_var*d_var**2*umz+2._ki*c_var**&
                &2*g_var+2._ki*c_var*d_var*g_var-d_var**3*umz)+1._ki/6._ki*log(e_va&
                &r)*e_var**3*(-1._ki+umz)/g_var**2/f_var**2
              !
            case(3)
              !
              fg=1._ki/6._ki*c_var**2/f_var**2/d_var**3*(log(umz)+z_log(s13,1._ki)&
                &)*(3._ki*d_var*f_var+c_var*d_var*umz+2._ki*c_var*f_var)*umz+1._ki/&
                &6._ki*h_var**2/g_var**2/d_var**3*(log(1._ki-umz)+log(umz)+z_log(s&
                &23,1._ki))*umz*(-d_var*g_var+d_var*h_var*umz+2._ki*c_var*g_var)-1&
                &._ki/6._ki/g_var**2*log(e_var)/f_var**2*e_var**3*umz+1._ki/6._ki/g_&
                &var/d_var**2/f_var*(d_var**2*g_var-c_var**2*d_var*umz-2._ki*c_va&
                &r*d_var**2*umz+2._ki*c_var**2*g_var+2._ki*c_var*d_var*g_var-d_var&
                &**3*umz)*umz
              !
            case(4)
              !
              fg=1._ki/6._ki*c_var**2/f_var**3/d_var**2*(log(umz)+z_log(s13,1._ki)&
                &)*umz*(3._ki*d_var*f_var+2._ki*c_var*d_var*umz+c_var*f_var)-1._ki/&
                &6._ki/g_var**2*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d_var**2&
                &*h_var**3*umz+1._ki/6._ki/g_var**2*log(e_var)/f_var**3*e_var**2*(&
                &f_var*e_var+g_var*f_var-2._ki*c_var*g_var*umz)+1._ki/6._ki/d_var/g&
                &_var/f_var**2*(-2._ki*h_var*d_var*g_var*umz+c_var**2*d_var*umz**&
                &2+d_var**3*umz**2+2._ki*c_var*d_var**2*umz**2+d_var*g_var**2+c_v&
                &ar**2*g_var*umz)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
              fg=-1._ki/3._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d&
                &_var**3*h_var*(-1._ki+umz)**2*(d_var**4*umz**2+2._ki*c_var*d_var*&
                &*3*umz**2+c_var**2*d_var**2*umz**2+c_var**2*g_var*d_var*umz-c_v&
                &ar*d_var*g_var**2-2._ki*g_var*d_var**3*umz+d_var**2*g_var**2-c_v&
                &ar*g_var*d_var**2*umz+c_var**2*g_var**2)+1._ki/6._ki*(c_var+d_var&
                &)/d_var**2/g_var**2*(-1._ki+umz)**2*(2._ki*c_var*g_var+2._ki*c_var&
                &*d_var*umz+2._ki*d_var**2*umz-3._ki*d_var*g_var)-1._ki/3._ki*(-1._ki&
                &+umz)**2*(-log(e_var)*e_var**3*d_var**3+c_var**3*g_var**3*log(u&
                &mz)+c_var**3*g_var**3*z_log(s13,1._ki))/g_var**3/f_var/d_var**3
              !
            case(3)
              !
              fg=-1._ki/3._ki*c_var**3/d_var**3/f_var*(log(umz)+z_log(s13,1._ki))*&
                &(1._ki-umz)*umz-1._ki/3._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log&
                &(s23,1._ki))/d_var**3*h_var*(1._ki-umz)*umz*(d_var**4*umz**2+2._ki&
                &*c_var*d_var**3*umz**2+c_var**2*d_var**2*umz**2+c_var**2*g_var*&
                &d_var*umz-c_var*d_var*g_var**2-2._ki*g_var*d_var**3*umz+d_var**2&
                &*g_var**2-c_var*g_var*d_var**2*umz+c_var**2*g_var**2)+1._ki/3._ki&
                &/g_var**3*log(e_var)/f_var*e_var**3*(1._ki-umz)*umz+1._ki/6._ki*(c&
                &_var+d_var)/d_var**2/g_var**2*(1._ki-umz)*(2._ki*c_var*g_var+2._ki&
                &*c_var*d_var*umz+2._ki*d_var**2*umz-3._ki*d_var*g_var)*umz
              !
            case(4)
              !
              fg=-1._ki/6._ki*c_var**3/d_var**2/f_var**2*(log(umz)+z_log(s13,1._ki&
                &))*(1._ki-umz)*umz+1._ki/6._ki*h_var**2/g_var**3/d_var**2*(log(1._k&
                &i-umz)+log(umz)+z_log(s23,1._ki))*(1._ki-umz)*umz*(-2._ki*d_var*g_&
                &var+2._ki*d_var*h_var*umz+c_var*g_var)-1._ki/6._ki/g_var**3*log(e_&
                &var)/f_var**2*e_var**2*(1._ki-umz)*(2._ki*f_var*e_var+2._ki*g_var*&
                &f_var-c_var*g_var*umz)-1._ki/6._ki/g_var**2/f_var/d_var*(1._ki-umz&
                &)*(2._ki*d_var**3*umz**2+4._ki*c_var*d_var**2*umz**2+2._ki*c_var**&
                &2*d_var*umz**2-4._ki*c_var*d_var*g_var*umz+d_var*g_var**2-3._ki*d&
                &_var**2*g_var*umz-c_var**2*g_var*umz)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
              fg=-1._ki/3._ki*c_var**3/d_var**3/f_var*(log(umz)+z_log(s13,1._ki))*&
                &umz**2-1._ki/3._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(s23,1._k&
                &i))/d_var**3*h_var*umz**2*(d_var**4*umz**2+2._ki*c_var*d_var**3*&
                &umz**2+c_var**2*d_var**2*umz**2+c_var**2*g_var*d_var*umz-c_var*&
                &d_var*g_var**2-2._ki*g_var*d_var**3*umz+d_var**2*g_var**2-c_var*&
                &g_var*d_var**2*umz+c_var**2*g_var**2)+1._ki/3._ki/g_var**3*log(e_&
                &var)/f_var*e_var**3*umz**2+1._ki/6._ki*(c_var+d_var)/d_var**2/g_v&
                &ar**2*umz**2*(2._ki*c_var*g_var+2._ki*c_var*d_var*umz+2._ki*d_var*&
                &*2*umz-3._ki*d_var*g_var)
              !
            case(4)
              !
              fg=-1._ki/6._ki*c_var**3/d_var**2/f_var**2*(log(umz)+z_log(s13,1._ki&
                &))*umz**2+1._ki/6._ki*h_var**2/g_var**3/d_var**2*(log(1._ki-umz)+l&
                &og(umz)+z_log(s23,1._ki))*umz**2*(-2._ki*d_var*g_var+2._ki*d_var*h&
                &_var*umz+c_var*g_var)-1._ki/6._ki/g_var**3*log(e_var)/f_var**2*e_&
                &var**2*umz*(2._ki*f_var*e_var+2._ki*g_var*f_var-c_var*g_var*umz)-&
                &1._ki/6._ki/g_var**2/f_var/d_var*(2._ki*d_var**3*umz**2+4._ki*c_var&
                &*d_var**2*umz**2+2._ki*c_var**2*d_var*umz**2-4._ki*c_var*d_var*g_&
                &var*umz+d_var*g_var**2-3._ki*d_var**2*g_var*umz-c_var**2*g_var*u&
                &mz)*umz
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
              fg=-1._ki/3._ki*c_var**3/d_var/f_var**3*(log(umz)+z_log(s13,1._ki))*&
                &umz**2-1._ki/3._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(s23,1._k&
                &i))/d_var*h_var**3*umz**2+1._ki/3._ki/g_var**3*log(e_var)/f_var**&
                &3*e_var*(-c_var*g_var**2*f_var*umz+c_var**2*f_var**2*umz**2+c_v&
                &ar**2*g_var**2*umz**2-c_var**2*g_var*f_var*umz**2+g_var**2*f_va&
                &r**2+f_var**4+2._ki*c_var*f_var**3*umz+c_var*g_var*f_var**2*umz+&
                &2._ki*g_var*f_var**3)+1._ki/6._ki/g_var**2/f_var**2*(g_var-c_var*u&
                &mz-d_var*umz)*(-2._ki*c_var*d_var*umz**2-2._ki*d_var**2*umz**2+g_&
                &var**2+d_var*g_var*umz+4._ki*c_var*g_var*umz)
              !
            case default
              !
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
            tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=1._ki/4._ki/f_var**4*(log(umz)+z_log(s13,1._ki))*c_var/d_var**4*(&
                  &c_var**2*d_var**2*umz**2+2._ki*c_var*d_var**2*f_var*umz+2._ki*d_v&
                  &ar*f_var**2*h_var+c_var**2*f_var**2)*(-2._ki*d_var*f_var-umz*c_v&
                  &ar*d_var-c_var*f_var)+1._ki/24._ki*c_var*(-1._ki+umz)*(-18._ki*c_va&
                  &r**2*d_var-42._ki*c_var*d_var**2-26._ki*d_var**3-18._ki*d_var*umz*&
                  &c_var**2-42._ki*umz*c_var*d_var**2-26._ki*umz*d_var**3+63._ki*c_va&
                  &r*g_var*d_var+52._ki*g_var*d_var**2+18._ki*c_var**2*g_var)/d_var*&
                  &*2/f_var**3+1._ki/4._ki/g_var*log(e_var)/f_var**4*e_var**4-7._ki/8&
                  &._ki*c_var**2/d_var**2/f_var**3*g_var**2-7._ki/4._ki*c_var**2/f_va&
                  &r**3-3._ki/4._ki*c_var**3/d_var/f_var**3+3._ki/4._ki*c_var**3/d_var&
                  &**2/f_var**3*g_var+21._ki/8._ki*c_var**2/d_var/f_var**3*g_var-13.&
                  &_ki/12._ki*c_var*d_var/f_var**3+13._ki/6._ki*c_var/f_var**3*g_var-1&
                  &3._ki/12._ki*c_var/d_var/f_var**3*g_var**2-1._ki/4._ki*c_var**3/d_v&
                  &ar**3/f_var**3*g_var**2-1._ki/4._ki/g_var*(log(1._ki-umz)+log(umz)&
                  &+z_log(s23,1._ki))/d_var**4*h_var**4
                !
              case(2)
                !
                fg=1._ki/12._ki*c_var**2/f_var**3/d_var**4*(log(umz)+z_log(s13,1._ki&
                  &))*(1._ki-umz)*(3._ki*c_var**2*f_var**2+c_var**2*d_var**2*umz**2+&
                  &6._ki*d_var**2*f_var**2+8._ki*c_var*d_var*f_var**2+2._ki*c_var**2*&
                  &d_var*f_var*umz+4._ki*c_var*d_var**2*f_var*umz)-1._ki/12._ki/g_var&
                  &**2*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d_var**4*h_var**3*&
                  &(1._ki-umz)*(-d_var*h_var*umz-3._ki*c_var*g_var+d_var*g_var)-1._ki&
                  &/4._ki/g_var/f_var**2*c_var*h_var*(1._ki-umz)*umz**2+1._ki/2._ki/f_&
                  &var**2*c_var*(1._ki-umz)*umz+1._ki/6._ki/f_var**2*d_var*(1._ki-umz)&
                  &*umz+19._ki/24._ki/f_var**2*c_var**2/d_var*(1._ki-umz)*umz+5._ki/12&
                  &._ki/f_var**2*c_var**3/d_var**2*(1._ki-umz)*umz-1._ki/12._ki/g_var/&
                  &f_var**2*d_var**2*(1._ki-umz)*umz**2-1._ki/12._ki/g_var/f_var**2*c&
                  &_var**3/d_var*(1._ki-umz)*umz**2+1._ki/24._ki*(-1._ki+umz)*(2._ki*g_&
                  &var**3*f_var*d_var**3+2._ki*log(e_var)*e_var**4*d_var**3+6._ki*c_&
                  &var**3*g_var**3*f_var+13._ki*c_var**2*g_var**3*f_var*d_var+6._ki*&
                  &c_var*g_var**3*f_var*d_var**2)/f_var**3/g_var**2/d_var**3
                !
              case(3)
                !
                fg=-1._ki/12._ki/g_var**2*log(e_var)/f_var**3*e_var**4*umz-1._ki/12.&
                  &_ki/g_var**2*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d_var**4*h&
                  &_var**3*umz*(-d_var*h_var*umz-3._ki*c_var*g_var+d_var*g_var)+1._k&
                  &i/12._ki*c_var**2/f_var**3/d_var**4*(log(umz)+z_log(s13,1._ki))*(&
                  &3._ki*c_var**2*f_var**2+c_var**2*d_var**2*umz**2+6._ki*d_var**2*f&
                  &_var**2+8._ki*c_var*d_var*f_var**2+2._ki*c_var**2*d_var*f_var*umz&
                  &+4._ki*c_var*d_var**2*f_var*umz)*umz-1._ki/4._ki/g_var/f_var**2*c_&
                  &var*h_var*umz**3+1._ki/2._ki/f_var**2*c_var*umz**2+1._ki/6._ki/f_va&
                  &r**2*d_var*umz**2+19._ki/24._ki/f_var**2*c_var**2/d_var*umz**2+5.&
                  &_ki/12._ki/f_var**2*c_var**3/d_var**2*umz**2-1._ki/4._ki*g_var/d_va&
                  &r/f_var**2*c_var*umz-1._ki/12._ki/g_var/f_var**2*c_var**3/d_var*u&
                  &mz**3-1._ki/12._ki/g_var/f_var**2*d_var**2*umz**3-13._ki/24._ki*g_v&
                  &ar/d_var**2/f_var**2*c_var**2*umz-1._ki/4._ki*g_var/d_var**3/f_va&
                  &r**2*c_var**3*umz-1._ki/12._ki*g_var*(-1._ki+umz)/f_var**2-1._ki/12&
                  &._ki*g_var/f_var**2
                !
              case(4)
                !
                fg=-1._ki/12._ki/g_var**2*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/&
                  &d_var**3*h_var**4*umz-1._ki/12._ki/g_var**2*log(e_var)/f_var**4*e&
                  &_var**3*(-f_var*e_var-g_var*f_var+3._ki*c_var*g_var*umz)+1._ki/12&
                  &._ki*c_var**2/d_var**3/f_var**4*(log(umz)+z_log(s13,1._ki))*(c_va&
                  &r**2*f_var**2+3._ki*c_var**2*d_var**2*umz**2+6._ki*d_var**2*f_var&
                  &**2+4._ki*c_var*d_var*f_var**2+2._ki*c_var**2*d_var*f_var*umz+8._k&
                  &i*c_var*d_var**2*f_var*umz)*umz+1._ki/4._ki/g_var/f_var**3*c_var*&
                  &d_var*h_var*umz**3+1._ki/4._ki*g_var/f_var**3*h_var*umz+1._ki/4._ki&
                  &*c_var**3/d_var/f_var**3*umz**2-1._ki/2._ki*c_var*d_var/f_var**3*&
                  &umz**2-1._ki/4._ki/f_var**3*d_var**2*umz**2+1._ki/24._ki*c_var**2/f&
                  &_var**3*umz**2+1._ki/12._ki/g_var/f_var**3*d_var**3*umz**3-1._ki/1&
                  &2._ki*c_var**3/d_var**2/f_var**3*g_var*umz-7._ki/24._ki*c_var**2/d&
                  &_var/f_var**3*g_var*umz+1._ki/12._ki/g_var/f_var**3*c_var**3*umz*&
                  &*3-1._ki/12._ki*g_var**2/f_var**3
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=-1._ki/12._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/&
                  &d_var**4*h_var**2*(-1._ki+umz)**2*(c_var**2*d_var**2*umz**2+2._ki&
                  &*c_var*d_var**3*umz**2+d_var**4*umz**2+3._ki*c_var**2*g_var**2+d&
                  &_var**2*g_var**2-2._ki*c_var*d_var*g_var**2+2._ki*c_var**2*d_var*&
                  &g_var*umz-2._ki*d_var**3*g_var*umz)+1._ki/12._ki*c_var**3/d_var**4&
                  &/f_var**2*(log(umz)+z_log(s13,1._ki))*(-1._ki+umz)**2*(-4._ki*d_va&
                  &r*f_var-umz*c_var*d_var-3._ki*c_var*f_var)+1._ki/4._ki/g_var**2/f_&
                  &var*c_var*h_var*(-1._ki+umz)**2*umz**2-1._ki/3._ki/g_var/f_var*c_v&
                  &ar*(-1._ki+umz)**2*umz-5._ki/24._ki/g_var/f_var*d_var*(-1._ki+umz)*&
                  &*2*umz+1._ki/12._ki/g_var**2/f_var*d_var**2*(-1._ki+umz)**2*umz**2&
                  &-1._ki/24._ki/g_var/f_var*c_var**2/d_var*(-1._ki+umz)**2*umz+1._ki/&
                  &12._ki/g_var**2/f_var*c_var**3/d_var*(-1._ki+umz)**2*umz**2+1._ki/&
                  &12._ki/g_var/f_var*c_var**3/d_var**2*(-1._ki+umz)**2*umz+1._ki/24.&
                  &_ki*(-1._ki+umz)**2*(3._ki*g_var**3*f_var*d_var**3+2._ki*log(e_var)&
                  &*e_var**4*d_var**3-5._ki*c_var**2*g_var**3*f_var*d_var-6._ki*c_va&
                  &r**3*g_var**3*f_var+2._ki*c_var*g_var**3*f_var*d_var**2)/f_var**&
                  &2/g_var**3/d_var**3
                !
              case(3)
                !
                fg=1._ki/12._ki/g_var**3*log(e_var)/f_var**2*e_var**4*(1._ki-umz)*um&
                  &z+1._ki/12._ki*c_var**3/d_var**4/f_var**2*(log(umz)+z_log(s13,1._k&
                  &i))*(1._ki-umz)*(-4._ki*d_var*f_var-umz*c_var*d_var-3._ki*c_var*f_&
                  &var)*umz-1._ki/12._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(s23,&
                  &1._ki))/d_var**4*h_var**2*(1._ki-umz)*umz*(c_var**2*d_var**2*umz*&
                  &*2+2._ki*c_var*d_var**3*umz**2+d_var**4*umz**2+3._ki*c_var**2*g_v&
                  &ar**2+d_var**2*g_var**2-2._ki*c_var*d_var*g_var**2+2._ki*c_var**2&
                  &*d_var*g_var*umz-2._ki*d_var**3*g_var*umz)+1._ki/4._ki/g_var**2/f_&
                  &var*c_var*h_var*(1._ki-umz)*umz**3-1._ki/3._ki/g_var/f_var*c_var*(&
                  &1._ki-umz)*umz**2-5._ki/24._ki/g_var/f_var*d_var*(1._ki-umz)*umz**2&
                  &+1._ki/12._ki/f_var*c_var/d_var*(1._ki-umz)*umz-5._ki/24._ki/f_var*c&
                  &_var**2/d_var**2*(1._ki-umz)*umz+1._ki/12._ki/g_var**2/f_var*c_var&
                  &**3/d_var*(1._ki-umz)*umz**3-1._ki/4._ki/f_var*c_var**3/d_var**3*(&
                  &1._ki-umz)*umz-1._ki/24._ki/g_var/f_var*c_var**2/d_var*(1._ki-umz)*&
                  &umz**2+1._ki/12._ki/g_var/f_var*c_var**3/d_var**2*(1._ki-umz)*umz*&
                  &*2+1._ki/12._ki/g_var**2/f_var*d_var**2*(1._ki-umz)*umz**3-1._ki/8.&
                  &_ki*umz*(-1._ki+umz)/f_var
                !
              case(4)
                !
                fg=-1._ki/12._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/&
                  &d_var**3*h_var**3*(1._ki-umz)*umz*(-d_var*h_var*umz-c_var*g_var+&
                  &d_var*g_var)+1._ki/12._ki*c_var**3/f_var**3/d_var**3*(log(umz)+z_&
                  &log(s13,1._ki))*(1._ki-umz)*umz*(-2._ki*d_var*f_var-umz*c_var*d_va&
                  &r-c_var*f_var)+1._ki/12._ki/g_var**3*log(e_var)/f_var**3*e_var**3&
                  &*(1._ki-umz)*(-f_var*e_var-g_var*f_var+c_var*g_var*umz)-1._ki/4._k&
                  &i/g_var**2/f_var**2*c_var*d_var*h_var*(1._ki-umz)*umz**3-1._ki/4.&
                  &_ki/f_var**2*c_var*(1._ki-umz)*umz-1._ki/6._ki/f_var**2*d_var*(1._ki&
                  &-umz)*umz+3._ki/8._ki/f_var**2/g_var*c_var**2*(1._ki-umz)*umz**2+5&
                  &._ki/24._ki/g_var/f_var**2*d_var**2*(1._ki-umz)*umz**2+1._ki/12._ki/&
                  &g_var/f_var**2*c_var**3/d_var*(1._ki-umz)*umz**2-1._ki/12._ki/f_va&
                  &r**2*c_var**3/d_var**2*(1._ki-umz)*umz-1._ki/8._ki/f_var**2*c_var*&
                  &*2/d_var*(1._ki-umz)*umz-1._ki/12._ki/g_var**2/f_var**2*d_var**3*(&
                  &1._ki-umz)*umz**3+1._ki/2._ki/f_var**2/g_var*c_var*d_var*(1._ki-umz&
                  &)*umz**2-1._ki/12._ki/g_var**2/f_var**2*c_var**3*(1._ki-umz)*umz**&
                  &3-1._ki/24._ki*g_var*(-1._ki+umz)/f_var**2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=1._ki/12._ki/g_var**3*log(e_var)/f_var**2*e_var**4*umz**2-1._ki/1&
                  &2._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d_var**4&
                  &*h_var**2*umz**2*(c_var**2*d_var**2*umz**2+2._ki*c_var*d_var**3*&
                  &umz**2+d_var**4*umz**2+3._ki*c_var**2*g_var**2+d_var**2*g_var**2&
                  &-2._ki*c_var*d_var*g_var**2+2._ki*c_var**2*d_var*g_var*umz-2._ki*d&
                  &_var**3*g_var*umz)+1._ki/12._ki*c_var**3/d_var**4/f_var**2*(log(u&
                  &mz)+z_log(s13,1._ki))*(-4._ki*d_var*f_var-umz*c_var*d_var-3._ki*c_&
                  &var*f_var)*umz**2+1._ki/4._ki/g_var**2/f_var*c_var*h_var*umz**4-1&
                  &._ki/3._ki/g_var/f_var*c_var*umz**3-5._ki/24._ki/g_var/f_var*d_var*&
                  &umz**3+1._ki/12._ki/g_var**2/f_var*d_var**2*umz**4-1._ki/4._ki/f_va&
                  &r*c_var**3/d_var**3*umz**2+1._ki/12._ki/f_var*c_var/d_var*umz**2-&
                  &5._ki/24._ki/f_var*c_var**2/d_var**2*umz**2+1._ki/12._ki/g_var**2/f&
                  &_var*c_var**3/d_var*umz**4+1._ki/12._ki/g_var/f_var*c_var**3/d_va&
                  &r**2*umz**3-1._ki/24._ki/g_var/f_var*c_var**2/d_var*umz**3+1._ki/8&
                  &._ki*(-1._ki+umz**2)/f_var+1._ki/8._ki/f_var
                !
              case(4)
                !
                fg=1._ki/12._ki*c_var**3/f_var**3/d_var**3*(log(umz)+z_log(s13,1._ki&
                  &))*umz**2*(-2._ki*d_var*f_var-umz*c_var*d_var-c_var*f_var)+1._ki/&
                  &12._ki/g_var**3*log(e_var)/f_var**3*e_var**3*(-f_var*e_var-g_var&
                  &*f_var+c_var*g_var*umz)*umz-1._ki/12._ki/g_var**3*(log(1._ki-umz)+&
                  &log(umz)+z_log(s23,1._ki))/d_var**3*h_var**3*umz**2*(-d_var*h_va&
                  &r*umz-c_var*g_var+d_var*g_var)-1._ki/4._ki/g_var**2/f_var**2*c_va&
                  &r*d_var*h_var*umz**4-1._ki/4._ki/f_var**2*c_var*umz**2-1._ki/6._ki/&
                  &f_var**2*d_var*umz**2-1._ki/12._ki/f_var**2*c_var**3/d_var**2*umz&
                  &**2-1._ki/8._ki/f_var**2*c_var**2/d_var*umz**2+1._ki/12._ki/g_var/f&
                  &_var**2*c_var**3/d_var*umz**3+1._ki/2._ki/f_var**2/g_var*c_var*d_&
                  &var*umz**3-1._ki/12._ki/g_var**2/f_var**2*c_var**3*umz**4+5._ki/24&
                  &._ki/g_var/f_var**2*d_var**2*umz**3+3._ki/8._ki/f_var**2/g_var*c_v&
                  &ar**2*umz**3-1._ki/12._ki/g_var**2/f_var**2*d_var**3*umz**4+1._ki/&
                  &24._ki*g_var*(-1._ki+umz)/f_var**2+1._ki/24._ki*g_var/f_var**2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=-1._ki/12._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/&
                  &d_var**2*h_var**4*umz**2+1._ki/12._ki/g_var**3*log(e_var)/f_var**&
                  &4*e_var**2*(f_var**4+3._ki*c_var**2*g_var**2*umz**2+c_var**2*f_v&
                  &ar**2*umz**2-2._ki*c_var**2*g_var*f_var*umz**2+2._ki*g_var*f_var*&
                  &*3+g_var**2*f_var**2+2._ki*f_var**3*c_var*umz-2._ki*c_var*g_var**&
                  &2*f_var*umz)+1._ki/12._ki*c_var**3/d_var**2/f_var**4*(log(umz)+z_&
                  &log(s13,1._ki))*umz**2*(-4._ki*d_var*f_var-3._ki*umz*c_var*d_var-c&
                  &_var*f_var)+1._ki/4._ki/g_var**2/f_var**3*c_var*d_var**2*h_var*um&
                  &z**4-1._ki/6._ki*c_var/f_var**3*g_var*umz+1._ki/24._ki*g_var/f_var*&
                  &*3*d_var*umz-1._ki/4._ki/g_var/f_var**3*c_var**3*umz**3+11._ki/24.&
                  &_ki*c_var**2/f_var**3*umz**2+1._ki/12._ki/g_var**2/f_var**3*d_var*&
                  &*4*umz**4-5._ki/24._ki/g_var/f_var**3*d_var**3*umz**3+7._ki/12._ki*&
                  &c_var*d_var/f_var**3*umz**2+1._ki/12._ki/g_var**2/f_var**3*c_var*&
                  &*3*d_var*umz**4-2._ki/3._ki/g_var/f_var**3*c_var*d_var**2*umz**3-&
                  &17._ki/24._ki/g_var/f_var**3*c_var**2*d_var*umz**3+1._ki/8._ki/f_va&
                  &r**3*d_var**2*umz**2-1._ki/12._ki*c_var**3/d_var/f_var**3*umz**2-&
                  &1._ki/24._ki*g_var**2/f_var**3
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=1._ki/4._ki/g_var**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d_&
                  &var**4*h_var*(-1._ki+umz)**3*(-d_var*h_var*umz-c_var*g_var+d_var&
                  &*g_var)*(-2._ki*d_var**2*h_var*g_var*umz+c_var**2*d_var**2*umz**&
                  &2+2._ki*c_var*d_var**3*umz**2+d_var**4*umz**2+c_var**2*g_var**2+&
                  &d_var**2*g_var**2)+3._ki/4._ki/g_var**3*c_var*h_var*(-1._ki+umz)**&
                  &3*umz**2-1._ki/g_var**2*c_var*(-1._ki+umz)**3*umz-5._ki/8._ki/g_var&
                  &**2*d_var*(-1._ki+umz)**3*umz+1._ki/4._ki/g_var**2*c_var**3/d_var*&
                  &*2*(-1._ki+umz)**3*umz+1._ki/4._ki/g_var**3*c_var**3/d_var*(-1._ki+&
                  &umz)**3*umz**2-1._ki/8._ki/g_var**2*c_var**2/d_var*(-1._ki+umz)**3&
                  &*umz+1._ki/4._ki/g_var**3*d_var**2*(-1._ki+umz)**3*umz**2-1._ki/24.&
                  &_ki*(-1._ki+umz)**3*(6._ki*c_var**4*g_var**4*log(umz)+6._ki*c_var**&
                  &4*g_var**4*z_log(s13,1._ki)-6._ki*log(e_var)*e_var**4*d_var**4+3.&
                  &_ki*c_var**2*f_var*d_var**2*g_var**3-11._ki*f_var*d_var**4*g_var*&
                  &*3-2._ki*c_var*f_var*d_var**3*g_var**3-6._ki*c_var**3*f_var*d_var&
                  &*g_var**3)/f_var/d_var**4/g_var**4
                !
              case(3)
                !
                fg=-1._ki/4._ki/g_var**4*log(e_var)/f_var*e_var**4*(-1._ki+umz)**2*u&
                  &mz-1._ki/4._ki/g_var**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/&
                  &d_var**4*h_var*(-1._ki+umz)**2*umz*(-d_var*h_var*umz-c_var*g_var&
                  &+d_var*g_var)*(-2._ki*d_var**2*h_var*g_var*umz+c_var**2*d_var**2&
                  &*umz**2+2._ki*c_var*d_var**3*umz**2+d_var**4*umz**2+c_var**2*g_v&
                  &ar**2+d_var**2*g_var**2)-3._ki/4._ki/g_var**3*c_var*h_var*(-1._ki+&
                  &umz)**2*umz**3+1._ki/g_var**2*c_var*(-1._ki+umz)**2*umz**2+5._ki/8&
                  &._ki/g_var**2*d_var*(-1._ki+umz)**2*umz**2+1._ki/4._ki*c_var**4/f_v&
                  &ar/d_var**4*(log(umz)+z_log(s13,1._ki))*(-1._ki+umz)**2*umz-1._ki/&
                  &4._ki/g_var**2*c_var**3/d_var**2*(-1._ki+umz)**2*umz**2-1._ki/4._ki&
                  &/g_var**3*c_var**3/d_var*(-1._ki+umz)**2*umz**3+1._ki/8._ki/d_var*&
                  &*2/g_var*c_var**2*(-1._ki+umz)**2*umz-1._ki/4._ki/d_var**3/g_var*c&
                  &_var**3*(-1._ki+umz)**2*umz-1._ki/12._ki/d_var/g_var*c_var*(-1._ki+&
                  &umz)**2*umz+1._ki/8._ki/g_var**2*c_var**2/d_var*(-1._ki+umz)**2*um&
                  &z**2-1._ki/4._ki/g_var**3*d_var**2*(-1._ki+umz)**2*umz**3-11._ki/24&
                  &._ki*(-1._ki+umz)**2*umz/g_var
                !
              case(4)
                !
                fg=-1._ki/12._ki/g_var**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/&
                  &d_var**3*h_var**2*(-1._ki+umz)**2*umz*(3._ki*d_var**4*umz**2+3._ki&
                  &*c_var**2*d_var**2*umz**2+6._ki*c_var*d_var**3*umz**2-2._ki*c_var&
                  &*d_var*g_var**2+3._ki*d_var**2*g_var**2-4._ki*c_var*g_var*d_var**&
                  &2*umz+2._ki*c_var**2*d_var*g_var*umz+c_var**2*g_var**2-6._ki*d_va&
                  &r**3*g_var*umz)-1._ki/12._ki/g_var**4*log(e_var)/f_var**2*e_var**&
                  &3*(-1._ki+umz)**2*(-3._ki*f_var*e_var-3._ki*g_var*f_var+c_var*g_va&
                  &r*umz)+3._ki/4._ki/f_var/g_var**3*c_var*d_var*h_var*(-1._ki+umz)**&
                  &2*umz**3+7._ki/12._ki/g_var/f_var*c_var*(-1._ki+umz)**2*umz+11._ki/&
                  &24._ki/g_var/f_var*d_var*(-1._ki+umz)**2*umz+1._ki/12._ki*c_var**4/&
                  &d_var**3/f_var**2*(log(umz)+z_log(s13,1._ki))*(-1._ki+umz)**2*umz&
                  &+1._ki/4._ki/f_var/g_var**3*d_var**3*(-1._ki+umz)**2*umz**3-19._ki/&
                  &24._ki/f_var/g_var**2*c_var**2*(-1._ki+umz)**2*umz**2+1._ki/4._ki/f&
                  &_var/g_var**3*c_var**3*(-1._ki+umz)**2*umz**3-4._ki/3._ki/f_var/g_&
                  &var**2*c_var*d_var*(-1._ki+umz)**2*umz**2-5._ki/8._ki/g_var**2/f_v&
                  &ar*d_var**2*(-1._ki+umz)**2*umz**2-1._ki/12._ki/g_var/f_var*c_var*&
                  &*3/d_var**2*(-1._ki+umz)**2*umz+1._ki/24._ki/g_var/f_var*c_var**2/&
                  &d_var*(-1._ki+umz)**2*umz-1._ki/12._ki/g_var**2/f_var*c_var**3/d_v&
                  &ar*(-1._ki+umz)**2*umz**2-1._ki/12._ki*(-1._ki+umz)**2/f_var
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=-1._ki/4._ki/g_var**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d&
                  &_var**4*h_var*(1._ki-umz)*umz**2*(-d_var*h_var*umz-c_var*g_var+d&
                  &_var*g_var)*(-2._ki*d_var**2*h_var*g_var*umz+c_var**2*d_var**2*u&
                  &mz**2+2._ki*c_var*d_var**3*umz**2+d_var**4*umz**2+c_var**2*g_var&
                  &**2+d_var**2*g_var**2)-1._ki/4._ki/g_var**4*log(e_var)/f_var*e_va&
                  &r**4*(1._ki-umz)*umz**2-3._ki/4._ki/g_var**3*c_var*h_var*(1._ki-umz&
                  &)*umz**4+1._ki/g_var**2*c_var*(1._ki-umz)*umz**3+5._ki/8._ki/g_var*&
                  &*2*d_var*(1._ki-umz)*umz**3-1._ki/4._ki/g_var**3*d_var**2*(1._ki-um&
                  &z)*umz**4+1._ki/8._ki/g_var**2*c_var**2/d_var*(1._ki-umz)*umz**3+1&
                  &._ki/4._ki*c_var**4/f_var/d_var**4*(log(umz)+z_log(s13,1._ki))*(1.&
                  &_ki-umz)*umz**2-1._ki/12._ki/d_var/g_var*c_var*(1._ki-umz)*umz**2-1&
                  &._ki/4._ki/d_var**3/g_var*c_var**3*(1._ki-umz)*umz**2+1._ki/8._ki/d_&
                  &var**2/g_var*c_var**2*(1._ki-umz)*umz**2-1._ki/4._ki/g_var**3*c_va&
                  &r**3/d_var*(1._ki-umz)*umz**4-1._ki/4._ki/g_var**2*c_var**3/d_var*&
                  &*2*(1._ki-umz)*umz**3+11._ki/24._ki*umz**2*(-1._ki+umz)/g_var
                !
              case(4)
                !
                fg=-1._ki/12._ki/g_var**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/&
                  &d_var**3*h_var**2*(1._ki-umz)*umz**2*(3._ki*d_var**4*umz**2+3._ki*&
                  &c_var**2*d_var**2*umz**2+6._ki*c_var*d_var**3*umz**2-2._ki*c_var*&
                  &d_var*g_var**2+3._ki*d_var**2*g_var**2-4._ki*c_var*g_var*d_var**2&
                  &*umz+2._ki*c_var**2*d_var*g_var*umz+c_var**2*g_var**2-6._ki*d_var&
                  &**3*g_var*umz)-1._ki/12._ki/g_var**4*log(e_var)/f_var**2*e_var**3&
                  &*(1._ki-umz)*(-3._ki*f_var*e_var-3._ki*g_var*f_var+c_var*g_var*umz&
                  &)*umz+3._ki/4._ki/f_var/g_var**3*c_var*d_var*h_var*(1._ki-umz)*umz&
                  &**4+7._ki/12._ki/g_var/f_var*c_var*(1._ki-umz)*umz**2+11._ki/24._ki/&
                  &g_var/f_var*d_var*(1._ki-umz)*umz**2-5._ki/8._ki/g_var**2/f_var*d_&
                  &var**2*(1._ki-umz)*umz**3+1._ki/12._ki*c_var**4/d_var**3/f_var**2*&
                  &(log(umz)+z_log(s13,1._ki))*(1._ki-umz)*umz**2-1._ki/12._ki/g_var**&
                  &2/f_var*c_var**3/d_var*(1._ki-umz)*umz**3-1._ki/12._ki/g_var/f_var&
                  &*c_var**3/d_var**2*(1._ki-umz)*umz**2+1._ki/24._ki/g_var/f_var*c_v&
                  &ar**2/d_var*(1._ki-umz)*umz**2+1._ki/4._ki/f_var/g_var**3*d_var**3&
                  &*(1._ki-umz)*umz**4-19._ki/24._ki/f_var/g_var**2*c_var**2*(1._ki-um&
                  &z)*umz**3+1._ki/4._ki/f_var/g_var**3*c_var**3*(1._ki-umz)*umz**4-4&
                  &._ki/3._ki/f_var/g_var**2*c_var*d_var*(1._ki-umz)*umz**3+1._ki/12._k&
                  &i*umz*(-1._ki+umz)/f_var
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=-1._ki/12._ki/g_var**4*log(e_var)/f_var**3*e_var**2*(1._ki-umz)*(&
                  &6._ki*f_var**3*c_var*umz+c_var**2*g_var**2*umz**2-2._ki*c_var**2*&
                  &g_var*f_var*umz**2+3._ki*c_var**2*f_var**2*umz**2+4._ki*c_var*f_v&
                  &ar**2*g_var*umz+6._ki*g_var*f_var**3+3._ki*g_var**2*f_var**2-2._ki&
                  &*c_var*g_var**2*f_var*umz+3._ki*f_var**4)-1._ki/12._ki/g_var**4*(l&
                  &og(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d_var**2*h_var**3*(1._ki-&
                  &umz)*umz**2*(-3._ki*d_var*h_var*umz-c_var*g_var+3._ki*d_var*g_var&
                  &)-3._ki/4._ki/g_var**3/f_var**2*c_var*d_var**2*h_var*(1._ki-umz)*u&
                  &mz**4+1._ki/6._ki/f_var**2*c_var*(1._ki-umz)*umz+1._ki/24._ki/f_var*&
                  &*2*d_var*(1._ki-umz)*umz-1._ki/4._ki/g_var**3/f_var**2*c_var**3*d_&
                  &var*(1._ki-umz)*umz**4-1._ki/12._ki/g_var/f_var**2*c_var**3/d_var*&
                  &(1._ki-umz)*umz**2+1._ki/12._ki*c_var**4/d_var**2/f_var**3*(log(um&
                  &z)+z_log(s13,1._ki))*(1._ki-umz)*umz**2-11._ki/24._ki/g_var/f_var**&
                  &2*d_var**2*(1._ki-umz)*umz**2+5._ki/8._ki/g_var**2/f_var**2*d_var*&
                  &*3*(1._ki-umz)*umz**3+5._ki/12._ki/g_var**2/f_var**2*c_var**3*(1._k&
                  &i-umz)*umz**3-17._ki/24._ki/f_var**2/g_var*c_var**2*(1._ki-umz)*um&
                  &z**2-1._ki/4._ki/g_var**3/f_var**2*d_var**4*(1._ki-umz)*umz**4+35.&
                  &_ki/24._ki/g_var**2/f_var**2*c_var**2*d_var*(1._ki-umz)*umz**3-13.&
                  &_ki/12._ki/f_var**2/g_var*c_var*d_var*(1._ki-umz)*umz**2+5._ki/3._ki&
                  &/g_var**2/f_var**2*c_var*d_var**2*(1._ki-umz)*umz**3-1._ki/24._ki*&
                  &g_var*(-1._ki+umz)/f_var**2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=-1._ki/4._ki/g_var**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d&
                  &_var**4*h_var*umz**3*(-d_var*h_var*umz-c_var*g_var+d_var*g_var)&
                  &*(-2._ki*d_var**2*h_var*g_var*umz+c_var**2*d_var**2*umz**2+2._ki*&
                  &c_var*d_var**3*umz**2+d_var**4*umz**2+c_var**2*g_var**2+d_var**&
                  &2*g_var**2)-1._ki/4._ki/g_var**4*log(e_var)/f_var*e_var**4*umz**3&
                  &-3._ki/4._ki/g_var**3*c_var*h_var*umz**5+1._ki/g_var**2*c_var*umz*&
                  &*4+5._ki/8._ki/g_var**2*d_var*umz**4-1._ki/4._ki/g_var**3*d_var**2*&
                  &umz**5-1._ki/4._ki/g_var**2*c_var**3/d_var**2*umz**4+1._ki/4._ki*c_&
                  &var**4/f_var/d_var**4*(log(umz)+z_log(s13,1._ki))*umz**3-1._ki/4.&
                  &_ki/g_var**3*c_var**3/d_var*umz**5-1._ki/12._ki/g_var*c_var/d_var*&
                  &umz**3-1._ki/4._ki/g_var*c_var**3/d_var**3*umz**3+1._ki/8._ki/g_var&
                  &*c_var**2/d_var**2*umz**3+1._ki/8._ki/g_var**2*c_var**2/d_var*umz&
                  &**4-11._ki/24._ki*(-1._ki+umz)*(1._ki+umz+umz**2)/g_var-11._ki/24._ki&
                  &/g_var
                !
              case(4)
                !
                fg=-1._ki/12._ki/g_var**4*log(e_var)/f_var**2*e_var**3*umz**2*(-3._k&
                  &i*f_var*e_var-3._ki*g_var*f_var+c_var*g_var*umz)-1._ki/12._ki/g_va&
                  &r**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d_var**3*h_var**2&
                  &*umz**3*(c_var**2*g_var**2+3._ki*c_var**2*d_var**2*umz**2+6._ki*c&
                  &_var*d_var**3*umz**2+3._ki*d_var**4*umz**2-2._ki*c_var*d_var*g_va&
                  &r**2-6._ki*d_var**3*g_var*umz+3._ki*d_var**2*g_var**2+2._ki*c_var*&
                  &*2*d_var*g_var*umz-4._ki*c_var*g_var*d_var**2*umz)+3._ki/4._ki/g_v&
                  &ar**3/f_var*c_var*d_var*h_var*umz**5+7._ki/12._ki/g_var/f_var*c_v&
                  &ar*umz**3+11._ki/24._ki/g_var/f_var*d_var*umz**3-5._ki/8._ki/g_var*&
                  &*2/f_var*d_var**2*umz**4+1._ki/24._ki/g_var/f_var*c_var**2/d_var*&
                  &umz**3+1._ki/4._ki/g_var**3/f_var*c_var**3*umz**5-1._ki/12._ki/g_va&
                  &r/f_var*c_var**3/d_var**2*umz**3-1._ki/12._ki/g_var**2/f_var*c_va&
                  &r**3/d_var*umz**4+1._ki/4._ki/g_var**3/f_var*d_var**3*umz**5-4._ki&
                  &/3._ki/g_var**2/f_var*c_var*d_var*umz**4+1._ki/12._ki*c_var**4/d_v&
                  &ar**3/f_var**2*(log(umz)+z_log(s13,1._ki))*umz**3-19._ki/24._ki/g_&
                  &var**2/f_var*c_var**2*umz**4-1._ki/12._ki*(-1._ki+umz**2)/f_var-1.&
                  &_ki/12._ki/f_var
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
                fg=-1._ki/12._ki/g_var**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/&
                  &d_var**2*h_var**3*umz**3*(-3._ki*d_var*h_var*umz-c_var*g_var+3._k&
                  &i*d_var*g_var)-1._ki/12._ki/g_var**4*log(e_var)/f_var**3*e_var**2&
                  &*(6._ki*f_var**3*c_var*umz+c_var**2*g_var**2*umz**2-2._ki*c_var**&
                  &2*g_var*f_var*umz**2+3._ki*c_var**2*f_var**2*umz**2+4._ki*c_var*f&
                  &_var**2*g_var*umz+6._ki*g_var*f_var**3+3._ki*g_var**2*f_var**2-2.&
                  &_ki*c_var*g_var**2*f_var*umz+3._ki*f_var**4)*umz-3._ki/4._ki/g_var*&
                  &*3/f_var**2*c_var*d_var**2*h_var*umz**5+1._ki/6._ki/f_var**2*c_va&
                  &r*umz**2+1._ki/24._ki/f_var**2*d_var*umz**2-17._ki/24._ki/f_var**2/&
                  &g_var*c_var**2*umz**3+5._ki/8._ki/g_var**2/f_var**2*d_var**3*umz*&
                  &*4+1._ki/12._ki*c_var**4/d_var**2/f_var**3*(log(umz)+z_log(s13,1.&
                  &_ki))*umz**3+5._ki/12._ki/g_var**2/f_var**2*c_var**3*umz**4-1._ki/4&
                  &._ki/g_var**3/f_var**2*d_var**4*umz**5-13._ki/12._ki/f_var**2/g_va&
                  &r*c_var*d_var*umz**3+5._ki/3._ki/g_var**2/f_var**2*c_var*d_var**2&
                  &*umz**4+35._ki/24._ki/g_var**2/f_var**2*c_var**2*d_var*umz**4-1._k&
                  &i/4._ki/g_var**3/f_var**2*c_var**3*d_var*umz**5-1._ki/12._ki/g_var&
                  &/f_var**2*c_var**3/d_var*umz**3-11._ki/24._ki/g_var/f_var**2*d_va&
                  &r**2*umz**3+1._ki/24._ki*g_var*(-1._ki+umz)/f_var**2+1._ki/24._ki*g_&
                  &var/f_var**2
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = "par3 should be 3 or 4 but is %d0"
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
                fg=-1._ki/4._ki/g_var**4*(log(1._ki-umz)+log(umz)+z_log(s23,1._ki))/d&
                  &_var*h_var**4*umz**3-1._ki/4._ki/g_var**4*log(e_var)/f_var**4*e_v&
                  &ar*(-f_var*e_var-g_var*f_var+c_var*g_var*umz)*(c_var**2*f_var**&
                  &2*umz**2+c_var**2*g_var**2*umz**2+g_var**2*f_var**2+2._ki*f_var*&
                  &*2*g_var*e_var+f_var**4+2._ki*f_var**3*c_var*umz)+3._ki/4._ki/g_va&
                  &r**3/f_var**3*c_var*d_var**3*h_var*umz**5-1._ki/6._ki*g_var/f_var&
                  &**3*c_var*umz+1._ki/8._ki*g_var/f_var**3*d_var*umz+1._ki/4._ki/g_va&
                  &r**3/f_var**3*c_var**3*d_var**2*umz**5+3._ki/4._ki/g_var/f_var**3&
                  &*c_var**3*umz**3-1._ki/2._ki/f_var**3*c_var**2*umz**2+11._ki/24._ki&
                  &/g_var/f_var**3*d_var**3*umz**3-17._ki/8._ki/g_var**2/f_var**3*c_&
                  &var**2*d_var**2*umz**4+15._ki/8._ki/g_var/f_var**3*c_var**2*d_var&
                  &*umz**3-2._ki/g_var**2/f_var**3*c_var*d_var**3*umz**4-5._ki/8._ki/&
                  &g_var**2/f_var**3*d_var**4*umz**4-1._ki/6._ki/f_var**3*c_var*d_va&
                  &r*umz**2+1._ki/4._ki/g_var**3/f_var**3*d_var**5*umz**5+1._ki/4._ki*&
                  &c_var**4/d_var/f_var**4*(log(umz)+z_log(s13,1._ki))*umz**3-3._ki/&
                  &4._ki/g_var**2/f_var**3*c_var**3*d_var*umz**4+19._ki/12._ki/g_var/&
                  &f_var**3*c_var*d_var**2*umz**3-1._ki/8._ki/f_var**3*d_var**2*umz*&
                  &*2-1._ki/12._ki*g_var**2/f_var**3
                !
              case default
                !
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
              tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
            tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = "par2 should be 1, 2, 3 or 4 but is %d0"
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
          tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
          fg=-1._ki/6._ki/g_var**2*log(-e_var)/f_var*e_var**3+1._ki/6._ki*(-1._k&
            &i+z)/g_var**2*(log(z)+log(1._ki-z)+z_log(-s23,-1._ki))/d_var**2*h&
            &_var**2*((-1._ki+z)*d_var*h_var-c_var*g_var+2._ki*d_var*g_var)-1.&
            &_ki/6._ki*(-1._ki+z)*c_var**3/d_var**2/f_var*(log(1._ki-z)+z_log(-s&
            &13,-1._ki))+1._ki/3._ki/g_var*z*c_var+1._ki/6._ki*d_var/g_var*z-1._ki&
            &/3._ki/g_var*c_var-1._ki/6._ki/d_var/g_var*c_var**2+1._ki/6._ki/d_va&
            &r/g_var*z*c_var**2+4._ki/9._ki-1._ki/6._ki*d_var/g_var
          !
        else if (nb_par == 1) then
          !
          select case(par4)
          !
          case(1)
            !
            fg=1._ki/24._ki*c_var**3/d_var**3/f_var**2*(log(umz)+z_log(-s13,-1.&
              &_ki))*umz*(4._ki*d_var*f_var+c_var*d_var*umz+2._ki*c_var*f_var)-1.&
              &_ki/24._ki*h_var**3/d_var**3/g_var**2*(log(1._ki-umz)+log(umz)+z_l&
              &og(-s23,-1._ki))*umz*(-2._ki*c_var*g_var-d_var*h_var*umz+2._ki*d_v&
              &ar*g_var)-1._ki/144._ki/d_var**2/g_var/f_var*(-25._ki*g_var*d_var*&
              &*3*umz+18._ki*c_var*h_var*d_var**2*umz**2+6._ki*c_var**3*d_var*um&
              &z**2+6._ki*d_var**4*umz**2+19._ki*g_var**2*d_var**2-18._ki*c_var*d&
              &_var**2*g_var*umz-18._ki*c_var**2*d_var*g_var*umz-12._ki*c_var**3&
              &*g_var*umz)-1._ki/24._ki/g_var**2*log(-e_var)/f_var**2*e_var**4
            !
          case(2)
            !
            fg=-1._ki/12._ki*c_var**4/d_var**3/f_var*(log(umz)+z_log(-s13,-1._ki&
              &))*(1._ki-umz)*umz-1._ki/12._ki/g_var**3*(log(1._ki-umz)+log(umz)+z&
              &_log(-s23,-1._ki))/d_var**3*h_var**2*(1._ki-umz)*umz*(c_var**2*d_&
              &var**2*umz**2+2._ki*c_var*d_var**3*umz**2+d_var**4*umz**2+g_var*&
              &*2*c_var**2-2._ki*d_var*c_var*g_var**2-2._ki*c_var*d_var**2*g_var&
              &*umz+3._ki*g_var**2*d_var**2-3._ki*g_var*d_var**3*umz+c_var**2*d_&
              &var*g_var*umz)+1._ki/72._ki/g_var**2/d_var**2*(1._ki-umz)*(-15._ki*&
              &g_var*d_var**3*umz+18._ki*c_var*h_var*d_var**2*umz**2+6._ki*c_var&
              &**3*d_var*umz**2+6._ki*d_var**4*umz**2+16._ki*g_var**2*d_var**2-2&
              &4._ki*c_var*d_var**2*g_var*umz-3._ki*c_var**2*d_var*g_var*umz+6._k&
              &i*c_var**3*g_var*umz)-1._ki/12._ki*log(-e_var)*e_var**4*(-1._ki+um&
              &z)/g_var**3/f_var
            !
          case(3)
            !
            fg=-1._ki/12._ki*c_var**4/d_var**3/f_var*(log(umz)+z_log(-s13,-1._ki&
              &))*umz**2-1._ki/12._ki/g_var**3*(log(1._ki-umz)+log(umz)+z_log(-s2&
              &3,-1._ki))/d_var**3*h_var**2*umz**2*(c_var**2*d_var**2*umz**2+2.&
              &_ki*c_var*d_var**3*umz**2+d_var**4*umz**2+g_var**2*c_var**2-2._ki&
              &*d_var*c_var*g_var**2-2._ki*c_var*d_var**2*g_var*umz+3._ki*g_var*&
              &*2*d_var**2-3._ki*g_var*d_var**3*umz+c_var**2*d_var*g_var*umz)+1&
              &._ki/12._ki/g_var**3*log(-e_var)/f_var*e_var**4*umz+1._ki/72._ki/g_&
              &var**2/d_var**2*umz*(-15._ki*g_var*d_var**3*umz+18._ki*c_var*h_va&
              &r*d_var**2*umz**2+6._ki*c_var**3*d_var*umz**2+6._ki*d_var**4*umz*&
              &*2+16._ki*g_var**2*d_var**2-24._ki*c_var*d_var**2*g_var*umz-3._ki*&
              &c_var**2*d_var*g_var*umz+6._ki*c_var**3*g_var*umz)
            !
          case(4)
            !
            fg=-1._ki/24._ki*c_var**4/d_var**2/f_var**2*(log(umz)+z_log(-s13,-1&
              &._ki))*umz**2-1._ki/24._ki*h_var**3/d_var**2/g_var**3*(log(1._ki-um&
              &z)+log(umz)+z_log(-s23,-1._ki))*umz**2*(3._ki*d_var*g_var-2._ki*d_&
              &var*h_var*umz-c_var*g_var)-1._ki/24._ki/g_var**3*log(-e_var)/f_va&
              &r**2*e_var**3*(2._ki*f_var*e_var+3._ki*f_var*g_var-c_var*g_var*um&
              &z)-1._ki/144._ki/d_var/g_var**2/f_var*(12._ki*d_var**4*umz**3+36._k&
              &i*c_var*h_var*d_var**2*umz**3+13._ki*d_var*g_var**3-36._ki*c_var*&
              &*2*d_var*g_var*umz**2-g_var**2*d_var**2*umz+18._ki*d_var*c_var*g&
              &_var**2*umz-54._ki*c_var*d_var**2*g_var*umz**2+12._ki*d_var*c_var&
              &**3*umz**3-6._ki*c_var**3*g_var*umz**2-24._ki*g_var*d_var**3*umz*&
              &*2)
            !
          case default
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
          tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
        tab_erreur_par(1)%chaine = "In fg (function_4p2m_adj.f90):"
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
end module function_4p2m_adj
!
