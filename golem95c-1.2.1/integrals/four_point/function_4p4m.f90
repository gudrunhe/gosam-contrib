! 
!****h* src/integrals/four_point/function_4p4m
! NAME
!
!  Module function_4p4m
!
! USAGE
!
!  use function_4p4m
!
! DESCRIPTION
!
!  This module computes the six-dimensional and eight dimensional 
!  three mass four point function with or without Feynman parameters
!  in the numerator.
!
! OUTPUT
!
!  This module exports three functions f4p4m, f4p4m_c and f4
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
module function_4p4m
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
  public :: f4p4m,f4,f4p4m_c
  !
  contains
    !
    !****f* src/integrals/four_point/function_4p4m/f4p4m
    ! NAME
    !
    !  Function f4p4m
    !
    ! USAGE
    !
    !  real_dim_4 = f4p4m(dim,s24,s13,s12,s23,s34,par1,par2,par3,par4)
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
    !    real_dim_4 = f4p4m("n+2",s24,s13,s12,s23,s34,0,0,0,0)
    !  * a eight dimensional three mass four point function 
    !    with no Feynman parameters in the numerator:
    !    real_dim_4 = f4p4m("n+4",s24,s13,s12,s23,s34,0,0,0,0)
    !  * a six dimensional three mass four point function 
    !    with the Feynman parameter z1 in the numerator:
    !    real_dim_4 = f4p4m("n+2",s24,s13,s12,s23,s34,0,0,0,1)
    !  * a six dimensional three mass four point function 
    !    with the Feynman parameters z1^2*z2 in the numerator:
    !    real_dim_4 = f4p4m("n+2",s24,s13,s12,s23,s34,0,2,1,1)
    !
    !*****
    function f4p4m(dim,s24,s13,s14,s12,s23,s34,par1,par2,par3,par4)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: s24,s13,s14,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(4) :: f4p4m
      !
      integer :: nb_par
      real(ki) :: lamb,det_s
      real(ki) :: plus_grand
      real(ki) :: norma
      !
      par = (/par1,par2,par3,par4/)
      !
      s_mat(1,:) = (/0._ki,s12,s13,s14/)
      s_mat(2,:) = (/s12,0._ki,s23,s24/)
      s_mat(3,:) = (/s13,s23,0._ki,s34/)
      s_mat(4,:) = (/s14,s24,s34,0._ki/)
      ! on redefinit la matrice S de telle facon a ce que ses elements
      ! soient entre -1 et 1
      plus_grand = maxval(array=abs(s_mat))
      s_mat = s_mat/plus_grand
      !
      det_s = (-2._ki*s_mat(2,3)*s_mat(1,3)*s_mat(2,4)*s_mat(1,4)+s_mat(1,2&
        &)**2*s_mat(3,4)**2-2._ki*s_mat(1,2)*s_mat(3,4)*s_mat(2,3)*s_mat(&
        &1,4)-2._ki*s_mat(1,2)*s_mat(3,4)*s_mat(1,3)*s_mat(2,4)+s_mat(2,4&
        &)**2*s_mat(1,3)**2+s_mat(1,4)**2*s_mat(2,3)**2)
      !
      b(1)=-(-s_mat(1,4)*s_mat(2,3)**2-s_mat(2,4)**2*s_mat(1,3)+s_mat(1&
        &,2)*s_mat(3,4)*s_mat(2,3)-s_mat(1,2)*s_mat(3,4)**2+s_mat(3,4)*s&
        &_mat(2,3)*s_mat(1,4)+s_mat(3,4)*s_mat(1,3)*s_mat(2,4)+s_mat(1,2&
        &)*s_mat(3,4)*s_mat(2,4)-2._ki*s_mat(2,4)*s_mat(3,4)*s_mat(2,3)+s&
        &_mat(2,3)*s_mat(1,3)*s_mat(2,4)+s_mat(2,3)*s_mat(2,4)*s_mat(1,4&
        &))/det_s
      !
      b(2)=-(-s_mat(1,4)**2*s_mat(2,3)-s_mat(2,4)*s_mat(1,3)**2+s_mat(1&
        &,2)*s_mat(3,4)*s_mat(1,4)+s_mat(1,2)*s_mat(3,4)*s_mat(1,3)-s_ma&
        &t(1,2)*s_mat(3,4)**2+s_mat(3,4)*s_mat(2,3)*s_mat(1,4)+s_mat(3,4&
        &)*s_mat(1,3)*s_mat(2,4)+s_mat(2,3)*s_mat(1,3)*s_mat(1,4)-2._ki*s&
        &_mat(1,3)*s_mat(3,4)*s_mat(1,4)+s_mat(2,4)*s_mat(1,3)*s_mat(1,4&
        &))/det_s
      !
      b(3)=(s_mat(1,4)**2*s_mat(2,3)+s_mat(2,4)**2*s_mat(1,3)-s_mat(1,2&
        &)*s_mat(3,4)*s_mat(2,4)+s_mat(1,2)**2*s_mat(3,4)-s_mat(2,4)*s_m&
        &at(1,3)*s_mat(1,2)-s_mat(2,4)*s_mat(1,3)*s_mat(1,4)-s_mat(2,3)*&
        &s_mat(2,4)*s_mat(1,4)-s_mat(1,4)*s_mat(2,3)*s_mat(1,2)-s_mat(1,&
        &2)*s_mat(3,4)*s_mat(1,4)+2._ki*s_mat(2,4)*s_mat(1,4)*s_mat(1,2))&
        &/det_s
      !
      b(4)=(2._ki*s_mat(2,3)*s_mat(1,3)*s_mat(1,2)-s_mat(1,4)*s_mat(2,3)&
        &*s_mat(1,2)-s_mat(2,3)*s_mat(1,3)*s_mat(2,4)-s_mat(2,3)*s_mat(1&
        &,3)*s_mat(1,4)+s_mat(1,2)**2*s_mat(3,4)-s_mat(1,2)*s_mat(3,4)*s&
        &_mat(1,3)+s_mat(2,4)*s_mat(1,3)**2-s_mat(1,2)*s_mat(3,4)*s_mat(&
        &2,3)-s_mat(2,4)*s_mat(1,3)*s_mat(1,2)+s_mat(1,4)*s_mat(2,3)**2)&
        &/det_s
      !
      lamb=2._ki*(-s_mat(1,4)*s_mat(2,3)*s_mat(1,2)+s_mat(1,4)**2*s_mat(&
        &2,3)-s_mat(3,4)*s_mat(2,3)*s_mat(1,4)-s_mat(2,3)*s_mat(1,3)*s_m&
        &at(2,4)-s_mat(1,2)*s_mat(3,4)*s_mat(1,3)+s_mat(1,2)**2*s_mat(3,&
        &4)+s_mat(2,4)*s_mat(1,3)**2+s_mat(1,4)*s_mat(2,3)**2+s_mat(2,3)&
        &*s_mat(1,3)*s_mat(1,2)-s_mat(2,3)*s_mat(1,3)*s_mat(1,4)-s_mat(1&
        &,2)*s_mat(3,4)*s_mat(2,3)-s_mat(2,4)*s_mat(1,3)*s_mat(1,2)-s_ma&
        &t(3,4)*s_mat(1,3)*s_mat(2,4)-s_mat(2,3)*s_mat(2,4)*s_mat(1,4)+s&
        &_mat(2,4)*s_mat(3,4)*s_mat(2,3)+s_mat(2,4)**2*s_mat(1,3)-s_mat(&
        &2,4)*s_mat(1,3)*s_mat(1,4)-s_mat(1,2)*s_mat(3,4)*s_mat(2,4)-s_m&
        &at(1,2)*s_mat(3,4)*s_mat(1,4)+s_mat(2,4)*s_mat(1,4)*s_mat(1,2)+&
        &s_mat(1,3)*s_mat(3,4)*s_mat(1,4)+s_mat(1,2)*s_mat(3,4)**2)
        !
      sumb=lamb/det_s
      !
      invs(1,1)=2._ki*s_mat(2,4)*s_mat(3,4)*s_mat(2,3)/det_s
      !
      invs(1,2)=s_mat(3,4)*(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3)&
        &-s_mat(1,4)*s_mat(2,3))/det_s
      !
      invs(1,3)=-s_mat(2,4)*(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3&
        &)+s_mat(1,4)*s_mat(2,3))/det_s
      !
      invs(1,4)=-s_mat(2,3)*(s_mat(2,4)*s_mat(1,3)+s_mat(1,2)*s_mat(3,4&
        &)-s_mat(1,4)*s_mat(2,3))/det_s
      !
      invs(2,1) = invs(1,2)
      !
      invs(2,2)=2._ki*s_mat(1,3)*s_mat(3,4)*s_mat(1,4)/det_s
      !
      invs(2,3)=-s_mat(1,4)*(s_mat(2,4)*s_mat(1,3)+s_mat(1,2)*s_mat(3,4&
        &)-s_mat(1,4)*s_mat(2,3))/det_s
      !
      invs(2,4)=-s_mat(1,3)*(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3&
        &)+s_mat(1,4)*s_mat(2,3))/det_s
      !
      invs(3,1) = invs(1,3)
      !
      invs(3,2) = invs(2,3)
      !
      invs(3,3)=2._ki*s_mat(2,4)*s_mat(1,4)*s_mat(1,2)/det_s
      !
      invs(3,4)=s_mat(1,2)*(s_mat(1,2)*s_mat(3,4)-s_mat(2,4)*s_mat(1,3)&
        &-s_mat(1,4)*s_mat(2,3))/det_s
      !
      invs(4,1) = invs(1,4)
      !
      invs(4,2) = invs(2,4)
      !
      invs(4,3) = invs(3,4)
      !
      invs(4,4)=2._ki*s_mat(2,3)*s_mat(1,3)*s_mat(1,2)/det_s
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
      f4p4m = 0._ki
      !
      !~ if ( (rat_or_tot_par%rat_selected) .and. (abs(lamb) <= coupure_4p4m) ) then
      if ( (rat_or_tot_par%rat_selected) .and. (4._ki <= 2._ki) ) then
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function f4p4m (in file f4p4m.f90): &
        &the flag rat to compute the rational part is on &
        &and the program reachs a region of phase space in &
        &which det(G) = 0  Becareful that the rational part &
        &is not well behaved in this region&
        &Nevertheless if the user wants to  go on, he has to &
        &reduce the value of the parameter coupure_4p4m'
        call catch_exception(0)
      end if
      !
      !~ if (abs(sumb) > coupure_4p4m) then
      if (4._ki > 2._ki) then
        !
        ! analytic computation
        !
        if (dim == "n+2") then
          !
          f4p4m(3:4)= a4p4m_np2(s_mat(2,4),s_mat(1,3),&
                         &s_mat(1,4),s_mat(1,2),s_mat(2,3),s_mat(3,4),&
                         &par1,par2,par3,par4)/plus_grand
          !
        else if (dim == "n+4") then
          !
          f4p4m = a4p4m_np4(s_mat(2,4),s_mat(1,3),&
                      &s_mat(1,4),s_mat(1,2),s_mat(2,3),s_mat(3,4),&
                      &par1,par2,par3,par4)
          f4p4m(3) = f4p4m(3)-log(plus_grand)*norma
          !
        end if
        !
      else
        !
        ! numerical computation
        !
        ! not implemented
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
    end function f4p4m
    !
    !****f* src/integrals/four_point/function_4p4m/f4p4m_c
    ! NAME
    !
    !  Function f4p4m_c
    !
    ! USAGE
    !
    !  complex_dim_2 = f4p4m_c(dim,s24,s13,s14,s12,s23,s34,par1,par2,par3,par4)
    !
    ! DESCRIPTION
    !
    !  This function computes the same thing that the function f4p4m
    !
    ! INPUTS
    !
    !  * dim -- a character (dimension 3), dim="n+2" six dimensional 
    !           two adjacent mass four point function, dim="n+4" eight dimensional
    !           two adjacent mass four point function
    !  * s24 -- a real (type ki), the S matrix element 2,4
    !  * s13 -- a real (type ki), the S matrix element 1,3
    !  * s14 -- a real (type ki), the S matrix element 1,4
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
    !  see function f4p4m
    !
    !*****
    function f4p4m_c(dim,s24,s13,s14,s12,s23,s34,par1,par2,par3,par4)
      !
      character (len=3), intent (in) :: dim
      real(ki), intent (in) :: s24,s13,s14,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      complex(ki), dimension(2) :: f4p4m_c
      !
      real(ki), dimension(4) :: res4
      !
      res4 = f4p4m(dim,s24,s13,s14,s12,s23,s34,par1,par2,par3,par4)
      call to_complex(res4,f4p4m_c)
      !
    end function f4p4m_c
    !
    !****if* src/integrals/four_point/function_4p4m/a4p4m_np2
    ! NAME
    !
    !  recursive function a4p4m_np2
    !
    ! USAGE
    !
    !  real_dim_2 = a4p4m_np2(s24,s13,s14,s12,s23,s34,par1,par2,par3,par4)
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
    !  * s14 -- a real (type ki), the S matrix element 1,4
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
    recursive function a4p4m_np2(s24,s13,s14,s12,s23,s34,par1,par2,par3,par4) result(res_4p4m_np2)
      !
      real(ki), intent (in) :: s24,s13,s14,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(2) :: res_4p4m_np2
      !      
      integer, dimension(3) :: smj
      integer :: j
      integer :: nb_par_loc
      integer, dimension(4) :: par_loc,par_plus
      real(ki), dimension(6) :: truc1,truc2,truc3,truc4
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
          truc4 = resultat3(4,1,:)
          !
        else
          !
          truc4 = f3p_sc(s_mat,unpackb(ibclr(b_pro,4),3))
          resultat3(4,1,:) = truc4
          deja_calcule3(4,1) = .true.
          !
        end if
        !
        if (deja_calcule3(3,1)) then
          !
          truc3 = resultat3(3,1,:)
          !
        else
          !
          truc3 = f3p_sc(s_mat,unpackb(ibclr(b_pro,3),3))
          resultat3(3,1,:) = truc3
          deja_calcule3(3,1) = .true.
          !
        end if
        !
        if (deja_calcule3(2,1)) then
          !
          truc2 = resultat3(2,1,:)
          !
        else
          !
          truc2 = f3p_sc(s_mat,unpackb(ibclr(b_pro,2),3))
          resultat3(2,1,:) = truc2
          deja_calcule3(2,1) = .true.
          !
        end if
        !
        ctemp = f4(s24,s13,s14,s12,s23,s34)
        res_4p4m_np2(1) = (real(ctemp,ki) - b(1)*truc1(5) - b(2)*truc2(5) - b(3)*truc3(5) - b(4)*truc4(5))/sumb
        res_4p4m_np2(2) = (aimag(ctemp) - b(1)*truc1(6) - b(2)*truc2(6) - b(3)*truc3(6) - b(4)*truc4(6))/sumb
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
          temp0 = a4p4m_np2(s24,s13,s14,s12,s23,s34,0,0,0,0)
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
        res_4p4m_np2(1) = (temp0(1) + temp1(5) + temp2(5))/sumb
        res_4p4m_np2(2) = (temp0(2) + temp1(6) + temp2(6))/sumb
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
          temp10 = a4p4m_np2(s24,s13,s14,s12,s23,s34,0,0,0,par4)
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
          temp11 = a4p4m_np2(s24,s13,s14,s12,s23,s34,0,0,0,par3)
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
        res_4p4m_np2(1) = (temp0(1) + temp1(5) + temp2(5) + temp3(5)) &
                                *2._ki/3._ki/sumb
        res_4p4m_np2(2) = (temp0(2) + temp1(6) + temp2(6) + temp3(6)) &
                                *2._ki/3._ki/sumb
      !
      ! cas avec trois parametres de feynman au numerateur
      !
      else
        !
        temp10 = a4p4m_np2(s24,s13,s14,s12,s23,s34,0,0,par2,par3)
        temp11 = a4p4m_np2(s24,s13,s14,s12,s23,s34,0,0,par2,par4)
        temp12 = a4p4m_np2(s24,s13,s14,s12,s23,s34,0,0,par3,par4)
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
        res_4p4m_np2(1) = ( temp0(1) + temp1(5) + temp2(5) + temp3(5) &
                       + temp4(5) )/2._ki/sumb
        res_4p4m_np2(2) = ( temp0(2) + temp1(6) + temp2(6) + temp3(6) &
                       + temp4(6) )/2._ki/sumb
      end if
      !
    end function a4p4m_np2
    !
    !****if* src/integrals/four_point/function_4p4m/a4p4m_np4
    ! NAME
    !
    !  recursive function a4p4m_np4
    !
    ! USAGE
    !
    !  real_dim_4 = a4p4m_np4(s24,s13,s14,s12,s23,s34,par1,par2,par3,par4)
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
    !  * s14 -- a real (type ki), the S matrix element 1,4
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
    recursive function a4p4m_np4(s24,s13,s14,s12,s23,s34,par1,par2,par3,par4) result(res_4p4m_np4)
      !
      real(ki), intent (in) :: s24,s13,s14,s12,s23,s34
      integer, intent (in) :: par1,par2,par3,par4
      real(ki), dimension(4) :: res_4p4m_np4
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
          temp0 = a4p4m_np2(s24,s13,s14,s12,s23,s34,0,0,0,0)
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
        res_4p4m_np4(1) = (-temp1(1))/(3._ki*sumb)
        res_4p4m_np4(2) = (-temp1(2))/(3._ki*sumb)
        res_4p4m_np4(3) = (temp0(1)-temp1(3)-2._ki/3._ki*temp1(1))/(3._ki*sumb)
        res_4p4m_np4(4) = (temp0(2)-temp1(4)-2._ki/3._ki*temp1(2))/(3._ki*sumb)
      !
      ! cas avec un parametre de feynman au numerateur
      !
      else if (nb_par_loc == 1) then
        !
        temp0 = a4p4m_np2(s24,s13,s14,s12,s23,s34,0,0,0,par4)/3._ki
        temp1 = b(par4)*a4p4m_np4(s24,s13,s14,s12,s23,s34,0,0,0,0)
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
        res_4p4m_np4(1) = ( temp1(1)+temp2(1)+temp3(1) )/(2._ki*sumb)
        res_4p4m_np4(2) = ( temp1(2)+temp2(2)+temp3(2) )/(2._ki*sumb)
        res_4p4m_np4(3) = ( temp1(3)+temp1(1)/6._ki+temp2(3)+temp2(1)/2._ki &
                           +temp3(3)+temp3(1)/2._ki+temp0(1) )/(2._ki*sumb)
        res_4p4m_np4(4) = ( temp1(4)+temp1(2)/6._ki+temp2(4)+temp2(2)/2._ki &
                           +temp3(4)+temp3(2)/2._ki+temp0(2) )/(2._ki*sumb)
      !
      ! cas avec plus de un parametre de feynman au numerateur
      !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a4p4m_np4:'
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
    end function a4p4m_np4
    !
    !****f* src/integrals/four_point/function_4p4m/f4
    ! NAME
    !
    !  function f4
    !
    ! USAGE
    !
    !  complex = f4(s,t,s1,s2,s3,s4)
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
    !  * s1 -- a real (type ki), p1^2
    !  * s2 -- a real (type ki), p2^2
    !  * s3 -- a real (type ki), p3^2
    !  * s4 -- a real (type ki), p4^2
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
    !**********************************************************************
    !
    !     the _ki massless box function  !!
    !
    !     note: parameters adapted for double precision
    !           from Binoth et al.
    !
    function f4(s,t,s1,s2,s3,s4)
      !
      real(ki), intent(in) :: s,t,s1,s2,s3,s4
      complex(ki) :: f4
      !
      complex(ki) :: k12,k13,k14,k23,k24,k34
      complex(ki) :: a, b, c, d, e, f, discr, srdelta
      complex(ki) :: x1, x2
      complex(ki) :: eta
      !
      if (rat_or_tot_par%tot_selected) then
        !
        ! parameters:
        !
        eta= (1.e-25_ki,0._ki)
        !
        ! scaled momenta:
        ! kij = - ( p_{i} + ... + p_{j-1} )^2
        !
        k12 = -s1
        k13 = -s
        k14 = -s4
        k23 = -s2
        k24 = -t
        k34 = -s3 
        !
        ! short hands:
        !
        a = k24*k34
        b = k13*k24+k12*k34-k14*k23 
        c = k12*k13
        d = k23
        !
        e = (k34-i_*eta)/(k13-i_*eta) 
        f = (k24-i_*eta)/(k12-i_*eta)        
        !
        ! the discriminant:
        !
        discr   = (b/a)**2  - 4._ki*( c + i_*eta*d*1._ki )/( a )
        srdelta = sqrt( discr ) 
        !
        ! the roots:
        !
        x1 = -( b/a - srdelta )/2._ki 
        x2 = -( b/a + srdelta )/2._ki 
        !
        ! the massless box function:
        !
        f4 = 1._ki /a / srdelta*&
        ( + log(-x1)**2/2._ki &
          + cdilog(1._ki+x1*e) + pheta(-x1,e)*log(1._ki+x1*e)&
          + cdilog(1._ki+x1*f) + pheta(-x1,f)*log(1._ki+x1*f)&
          - log(-x1)*( log(k12-i_*eta)+log(k13-i_*eta)&
                          -log(k14-i_*eta)-log(k23-i_*eta) )&
          - log(-x2)**2/2._ki &
          - cdilog(1._ki+x2*e) - pheta(-x2,e)*log(1._ki+x2*e)&
          - cdilog(1._ki+x2*f) - pheta(-x2,f)*log(1._ki+x2*f)&
          + log(-x2)*( log(k12-i_*eta)+log(k13-i_*eta)&
          - log(k14-i_*eta)-log(k23-i_*eta) )   )
          !
      else if (rat_or_tot_par%rat_selected) then
        !
        f4 = 0._ki
        !
      end if
        !
      end function f4
!**********************************************************************
      function pheta(u,v)
        !
        complex(ki) :: u,v
        complex(ki) :: pheta
        !
        complex(ki) :: w
        !
        w = u*v
        !      pheta = i*2._ki*pi*( 
        !     &       +theta(-dimag(u))*theta(-dimag(v))*theta( dimag(w))
        !     &       -theta( dimag(u))*theta( dimag(v))*theta(-dimag(w)) )
        !     if (dimag(w).eq.0._ki.and.dble(w).le.0._ki) write (*,*) 'pheta !'
        pheta = log( w ) - log( u ) - log( v )
        !
      end function pheta
    !
    !
end module function_4p4m
!
