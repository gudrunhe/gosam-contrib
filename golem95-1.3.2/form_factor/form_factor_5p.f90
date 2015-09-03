!
!****h* src/form_factor/form_factor_5p
! NAME
!
!  Module form_factor_5p
!
! USAGE
!
!  use form_factor_5p
!
! DESCRIPTION
!
!  This module contains the different form factors for five point amplitudes.
!
! OUTPUT
!
!  It exports twelve functions:
!  * a50 -- a function to compute A^{5,0}
!  * a51 -- a function to compute A^{5,1}
!  * a52 -- a function to compute A^{5,2}
!  * a53 -- a function to compute A^{5,3}
!  * a54 -- a function to compute A^{5,4}
!  * a55 -- a function to compute A^{5,5}
!  * b52 -- a function to compute B^{5,2}
!  * b53 -- a function to compute B^{5,3}
!  * b54 -- a function to compute B^{5,4}
!  * b55 -- a function to compute B^{5,5}
!  * c54 -- a function to compute C^{5,4}
!  * c55 -- a function to compute C^{5,5}
!
!  Note that a5xx, b5xx and c5xx are generic functions which can be called either with a
!  set of integers or with an integer whose bits represents the set
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * generic_function_4p (src/integrals/four_point/generic_function_4p.f90)
!  * generic_function_3p (src/integrals/three_point/generic_function_3p.f90)
!  * matrice_s (src/kinematic/matrice_s.f90)
!  * array (src/module/array.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * multiply_div (src/module/multiply_div.f90)
!  * form_factor_type (src/module/form_factor_type.f90)
!
!*****
module form_factor_5p
  !
  use precision_golem
  use generic_function_4p
  use generic_function_3p
  use array
  use matrice_s
  use sortie_erreur
  use multiply_div
  use form_factor_type
  use constante, only: czero
  implicit none
  !
  private
  integer :: b_pin_glob,b_pro_glob
  !
  interface a50
    !
    module procedure a50_b, a50_s
    !
  end interface
  !
  interface a51
    !
    module procedure a51_b, a51_s
    !
  end interface
  !
  interface a52
    !
    module procedure a52_b, a52_s
    !
  end interface
  !
  interface a53
    !
    module procedure a53_b, a53_s
    !
  end interface
  !
  interface a54
    !
    module procedure a54_b, a54_s
    !
  end interface
  !
  interface a55
    !
    module procedure a55_b, a55_s
    !
  end interface
  !
  interface b52
    !
    module procedure b52_b, b52_s
    !
  end interface
  !
  interface b53
    !
    module procedure b53_b, b53_s
    !
  end interface
  !
  interface b54
    !
    module procedure b54_b, b54_s
    !
  end interface
  !
  interface b55
    !
    module procedure b55_b, b55_s
    !
  end interface
  !
  interface c54
    !
    module procedure c54_b, c54_s
    !
  end interface
  !
  interface c55
    !
    module procedure c55_b, c55_s
    !
  end interface
  !
  public :: a50,a51,a52,a53,a54,a55,b52,b53,b54,b55,c54,c55
  !
  contains
    !
    !****f* src/form_factor/form_factor_5p/a50_b
    ! NAME
    !
    !  Function a50_b
    !
    ! USAGE
    !
    !  type(form_factor) = a50_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,0}.
    ! 
    ! INPUTS
    !
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a50_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: a50_b
      !
      integer :: j,k
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj
      integer :: b_pro,b_pro_mj,b_pro_mjk
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp1 = czero
        temp2(:) = czero
        !
        ib = b_pro
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            b_pin_pj = punion( b_pin,ibset(0,j) )
            !
            temp1 = temp1 + b(j,b_pin)*sumb(b_pin_pj)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj)
            !
            ibj = b_pro_mj
            k = 0
            !
            second_pinch: do while (ibj /= 0)
              !
              if (modulo(ibj,2) == 1) then
                !
                b_pro_mjk = ibclr(b_pro_mj,k)
                !
                temp2 = temp2 + b(j,b_pin)*b(k,b_pin_pj)*f3p(s_mat_p,b_pro_mjk)
                !
              end if
              !
              k = k+1
              ibj = ishft(ibj,-1)
              !
            end do second_pinch
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp2(3) = temp2(3) + temp1
        a50_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a50'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a50_b
    !
    !****f* src/form_factor/form_factor_5p/a50_s
    ! NAME
    !
    !  Function a50_s
    !
    ! USAGE
    !
    !  type(form_factor) = a50_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,0}.
    ! 
    ! INPUTS
    !
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a50_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a50_s
      !
      a50_s = a50_b(packb(set))
      !
    end function a50_s
    !
    !****f* src/form_factor/form_factor_5p/a51_b
    ! NAME
    !
    !  Function a51_b
    !
    ! USAGE
    !
    !  type(form_factor) = a51_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,1}(l_1).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a51_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: a51_b
      !
      integer :: j,k
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj
      integer :: b_pro,b_pro_mj,b_pro_mjk
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp1 = czero
        temp2(:) = czero
        !
        ib = b_pro
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            b_pin_pj = punion( b_pin,ibset(0,j) )
            !
            temp1 = temp1 - inv_s(j,l1,b_pin)*sumb(b_pin_pj)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj)
            !
            ibj = b_pro_mj
            k = 0
            !
            second_pinch: do while (ibj /= 0)
              !
              if (modulo(ibj,2) == 1) then
                !
                b_pro_mjk = ibclr(b_pro_mj,k)
                !
                temp2 = temp2 - inv_s(j,l1,b_pin)*b(k,b_pin_pj)*f3p(s_mat_p,b_pro_mjk)
                !
              end if
              !
              k = k+1
              ibj = ishft(ibj,-1)
              !
            end do second_pinch
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp2(3) = temp2(3) + temp1
        a51_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a51'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a51_b
    !
    !****f* src/form_factor/form_factor_5p/a51_s
    ! NAME
    !
    !  Function a51_s
    !
    ! USAGE
    !
    !  type(form_factor) = a51_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,1}(l_1).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a51_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a51_s
      !
      a51_s = a51_b(l1,packb(set))
      !
    end function a51_s
    !
    !****f* src/form_factor/form_factor_5p/a52_b
    ! NAME
    !
    !  Function a52_b
    !
    ! USAGE
    !
    !  type(form_factor) = a52_b(l1,l2,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,2}(l1,l2).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a52_b(l1,l2,b_pin)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in) :: b_pin
      type(form_factor) :: a52_b
      !
      integer :: j,k
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj
      integer :: b_pro,b_pro_mj,b_pro_mjk
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp1 = czero
        temp2(:) = czero
        !
        ib = b_pro
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            b_pin_pj = punion( b_pin,ibset(0,j) )
            !
            temp1 = temp1 + ( inv_s(j,l1,b_pin)*b(l2,b_pin) + inv_s(j,l2,b_pin)*b(l1,b_pin) &
                             - 2._ki*inv_s(l1,l2,b_pin)*b(j,b_pin) &
                             + b(j,b_pin)*inv_s(l1,l2,b_pin_pj) )*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj)
            !
            ibj = b_pro_mj
            k = 0
            !
            second_pinch: do while (ibj /= 0)
              !
              if (modulo(ibj,2) == 1) then
                !
                b_pro_mjk = ibclr(b_pro_mj,k)
                !
                temp2 = temp2 + ( inv_s(j,l2,b_pin)*inv_s(k,l1,b_pin_pj) &
                                 + inv_s(j,l1,b_pin)*inv_s(k,l2,b_pin_pj) ) &
                                *f3p(s_mat_p,b_pro_mjk)/2._ki
                !
              end if
              !
              k = k+1
              ibj = ishft(ibj,-1)
              !
            end do second_pinch
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp2(3) = temp2(3) + temp1
        a52_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a52'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a52_b
    !
    !****f* src/form_factor/form_factor_5p/a52_s
    ! NAME
    !
    !  Function a52_s
    !
    ! USAGE
    !
    !  type(form_factor) = a52_s(l1,l2,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,2}(l1,l2).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a52_s(l1,l2,set)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a52_s
      !
      a52_s = a52_b(l1,l2,packb(set))
      !
    end function a52_s
    !
    !****f* src/form_factor/form_factor_5p/a53_b
    ! NAME
    !
    !  Function a53_b
    !
    ! USAGE
    !
    !  type(form_factor) = a53_b(l1,l2,l3,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,3}(l1,l2,l3).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a53_b(l1,l2,l3,b_pin)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in) :: b_pin
      type(form_factor) :: a53_b
      !
      complex(ki), dimension(3) :: t53
      !
      if (dim_s >= 5) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        t53 = f53(l1,l2,l3) + f53(l1,l3,l2) + f53(l3,l2,l1)
        a53_b = t53
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a53'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a53_b
    !
    !****f* src/form_factor/form_factor_5p/a53_s
    ! NAME
    !
    !  Function a53_s
    !
    ! USAGE
    !
    !  type(form_factor) = a53_s(l1,l2,l3,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,3}(l1,l2,l3).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a53_s(l1,l2,l3,set)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a53_s
      !
      a53_s = a53_b(l1,l2,l3,packb(set))
      !
    end function a53_s
    !
    !****if* src/form_factor/form_factor_5p/f53
    ! NAME
    !
    !  Function f53
    !
    ! USAGE
    !
    !  real_dim6 = f53(k1,k2,k3)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a53
    !
    ! INPUTS
    !
    !  * k1 -- an integer
    !  * k2 -- an integer
    !  * k3 -- an integer
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the global (to this module) variable b_pin_glob,b_pro_glob
    !  defined in a53
    !
    ! RETURN VALUE
    !
    !  It returns an array of six reals (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function f53(k1,k2,k3)
      !
      integer, intent(in) :: k1,k2,k3
      complex(ki), dimension(3) :: f53
      !
      integer :: j,k
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      complex(ki) :: bk1_cache, bk2_cache, invsk1k2_cache
      complex(ki) :: invsjk1_cache, invsjk2_cache
      integer :: ib,ibj
      integer :: b_pro_mj,b_pro_mjk
      integer :: b_pin_pj
      !
      temp1 = czero
      temp2(:) = czero
      !
      ib = b_pro_glob
      j = 0
      !

      bk1_cache = b(k1,b_pin_glob)
      bk2_cache = b(k2,b_pin_glob)
      invsk1k2_cache = inv_s(k1,k2,b_pin_glob)
      
      first_pinch: do while (ib /= 0)
         !
         if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro_glob,j)
            b_pin_pj = ibset(b_pin_glob,j)
            !
            invsjk1_cache = inv_s(j,k1,b_pin_glob)
            invsjk2_cache = inv_s(j,k2,b_pin_glob)

            temp1 = temp1 - ( invsjk1_cache*bk2_cache &
                 + invsjk2_cache*bk1_cache &
                 - 2._ki*invsk1k2_cache*b(j,b_pin_glob) &
                 + b(j,b_pin_glob)*inv_s(k1,k2,b_pin_pj) ) &
                 *f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,k3)*2._ki/3._ki &
                 + inv_s(j,k3,b_pin_glob)*inv_s(k1,k2,b_pin_pj) &
                 *f4p_np2(s_mat_p,b_pro_mj,b_pin_pj)/3._ki
            !
            ibj = b_pro_mj
            k = 0
            !
            
            second_pinch: do while (ibj /= 0)
               !
               if (modulo(ibj,2) == 1) then
                  !
                  b_pro_mjk = ibclr(b_pro_mj,k)
                  !
                  temp2 = temp2 - ( invsjk1_cache*inv_s(k,k2,b_pin_pj) &
                       + invsjk2_cache*inv_s(k,k1,b_pin_pj) ) &
                       *f3p(s_mat_p,b_pro_mjk,k3)/6._ki
                  !
               end if
               !
               k = k+1
               ibj = ishft(ibj,-1)
               !
            end do second_pinch
            !
         end if
         !
         j = j+1
         ib= ishft(ib,-1)
         !
      end do first_pinch
      !
      f53 = temp2
      f53(3) = f53(3) + temp1
      !
    end function f53
    !
    !****f* src/form_factor/form_factor_5p/a54_b
    ! NAME
    !
    !  Function a54_b
    !
    ! USAGE
    !
    !  type(form_factor) = a54_b(l1,l2,l3,l4,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,4}(l1,l2,l3,l4).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l4 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a54_b(l1,l2,l3,l4,b_pin)
      !
      integer, intent (in) :: l1,l2,l3,l4
      integer, intent (in) :: b_pin
      type(form_factor) :: a54_b
      !
      complex(ki), dimension(3) :: t54
      !
      if (dim_s >= 5) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        t54 = (  f54(l1,l2,l3,l4) + f54(l1,l3,l2,l4) + f54(l1,l4,l3,l2) &
               + f54(l3,l2,l1,l4) + f54(l4,l2,l3,l1) + f54(l3,l4,l1,l2) &
               + g54(l1,l2,l3,l4) + g54(l2,l1,l3,l4) + g54(l3,l2,l1,l4) &
               + g54(l4,l2,l3,l1) )/4._ki
        a54_b = t54
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a54'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a54_b
    !
    !****f* src/form_factor/form_factor_5p/a54_s
    ! NAME
    !
    !  Function a54_s
    !
    ! USAGE
    !
    !  type(form_factor) = a54_s(l1,l2,l3,l4,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,4}(l1,l2,l3,l4).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l4 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a54_s(l1,l2,l3,l4,set)
      !
      integer, intent (in) :: l1,l2,l3,l4
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a54_s
      !
      a54_s = a54_b(l1,l2,l3,l4,packb(set))
      !
    end function a54_s
    !
    !****if* src/form_factor/form_factor_5p/f54
    ! NAME
    !
    !  Function f54
    !
    ! USAGE
    !
    !  real_dim6 = f54(k1,k2,k3,k4)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a54
    !
    ! INPUTS
    !
    !  * k1 -- an integer
    !  * k2 -- an integer
    !  * k3 -- an integer
    !  * k4 -- an integer
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the global (to this module) variable b_pin_glob,b_pro_glob
    !  defined in a54
    !
    ! RETURN VALUE
    !
    !  It returns an array of six reals (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function f54(k1,k2,k3,k4)
      !
      integer, intent(in) :: k1,k2,k3,k4
      complex(ki), dimension(3) :: f54
      !
      integer :: j,k
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj
      integer :: b_pro_mj,b_pro_mjk
      integer :: b_pin_pj
      !
      temp1 = czero
      temp2(:) = czero
      !
      ib = b_pro_glob
      j = 0
      !
      first_pinch: do while (ib /= 0)
        !
        if (modulo(ib,2) == 1)  then
          !
          b_pro_mj = ibclr(b_pro_glob,j)
          b_pin_pj = punion( b_pin_glob,ibset(0,j) )
          !
          temp1 = temp1 - 2._ki*( 2._ki*inv_s(k3,k4,b_pin_glob)*b(j,b_pin_glob) &
                                - inv_s(j,k4,b_pin_glob)*b(k3,b_pin_glob) &
                                - inv_s(j,k3,b_pin_glob)*b(k4,b_pin_glob) &
                                - b(j,b_pin_glob)*inv_s(k3,k4,b_pin_pj) ) &
                          *f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,k1,k2)
          !
          ibj = b_pro_mj
          k = 0
          !
          second_pinch: do while (ibj /= 0)
            !
            if (modulo(ibj,2) == 1) then
              !
              b_pro_mjk = ibclr(b_pro_mj,k)
              !
              temp2 = temp2 + ( inv_s(j,k3,b_pin_glob)*inv_s(k,k4,b_pin_pj) &
                               + inv_s(j,k4,b_pin_glob)*inv_s(k,k3,b_pin_pj) ) &
                              *f3p(s_mat_p,b_pro_mjk,k1,k2)/3._ki
              !
            end if
            !
            k = k+1
            ibj = ishft(ibj,-1)
            !
          end do second_pinch
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      f54 = temp2
      f54(3) = f54(3) + temp1
      !
    end function f54
    !
    !****if* src/form_factor/form_factor_5p/g54
    ! NAME
    !
    !  Function g54
    !
    ! USAGE
    !
    !  real_dim6 = g54(k1,k2,k3,k4)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a54
    !
    ! INPUTS
    !
    !  * k1 -- an integer
    !  * k2 -- an integer
    !  * k3 -- an integer
    !  * k4 -- an integer
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the global (to this module) variable b_pin_glob,b_pro_glob
    !  defined in a54
    !
    ! RETURN VALUE
    !
    !  It returns an array of six reals (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function g54(k1,k2,k3,k4)
      !
      integer, intent(in) :: k1,k2,k3,k4
      complex(ki), dimension(3) :: g54
      !
      integer :: j
      complex(ki) :: temp1
      integer :: ib
      integer :: b_pro_mj
      integer :: b_pin_pj
      !
      temp1 = czero
      !
      ib = b_pro_glob
      j = 0
      !
      first_pinch: do while (ib /= 0)
        !
        if (modulo(ib,2) == 1)  then
          !
          b_pro_mj = ibclr(b_pro_glob,j)
          b_pin_pj = punion( b_pin_glob,ibset(0,j) )
          !
          temp1 = temp1 - 2._ki*( inv_s(j,k4,b_pin_glob)*inv_s(k2,k3,b_pin_pj) &
                                + inv_s(j,k3,b_pin_glob)*inv_s(k2,k4,b_pin_pj) &
                                + inv_s(j,k2,b_pin_glob)*inv_s(k3,k4,b_pin_pj) ) &
                          *f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,k1)/3._ki
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      g54(:) = czero
      g54(3) = g54(3) + temp1
      !
    end function g54
    !
    !****f* src/form_factor/form_factor_5p/a55_b
    ! NAME
    !
    !  Function a55_b
    !
    ! USAGE
    !
    !  type(form_factor) = a55_b(l1,l2,l3,l4,l5,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,5}(l1,l2,l3,l4,l5).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l4 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l5 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a55_b(l1,l2,l3,l4,l5,b_pin)
      !
      integer, intent (in) :: l1,l2,l3,l4,l5
      integer, intent (in) :: b_pin
      type(form_factor) :: a55_b
      !
      complex(ki), dimension(3) :: t55
      !
      if (dim_s >= 5) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        t55 = (  f55(l1,l2,l3,l4,l5) + f55(l1,l2,l4,l3,l5) &
               + f55(l1,l2,l5,l4,l3) + f55(l1,l4,l3,l2,l5) &
               + f55(l1,l5,l3,l4,l2) + f55(l4,l2,l3,l1,l5) &
               + f55(l5,l2,l3,l4,l1) + f55(l1,l4,l5,l2,l3) &
               + f55(l4,l2,l5,l1,l3) + f55(l4,l5,l3,l1,l2) & 
               + g55(l1,l2,l3,l4,l5) + g55(l1,l3,l2,l4,l5) &
               + g55(l3,l2,l1,l4,l5) + g55(l4,l2,l3,l1,l5) &
               + g55(l1,l4,l3,l2,l5) + g55(l5,l2,l3,l4,l1) &
               + g55(l1,l5,l3,l4,l2) + g55(l3,l4,l1,l2,l5) &
               + g55(l3,l5,l1,l4,l2) + g55(l4,l5,l3,l1,l2) )/5._ki
        a55_b = t55
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a55'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a55_b
    !
    !****f* src/form_factor/form_factor_5p/a55_s
    ! NAME
    !
    !  Function a55_s
    !
    ! USAGE
    !
    !  type(form_factor) = a55_s(l1,l2,l3,l4,l5,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{5,5}(l1,l2,l3,l4,l5).
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l4 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l5 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function a55_s(l1,l2,l3,l4,l5,set)
      !
      integer, intent (in) :: l1,l2,l3,l4,l5
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a55_s
      !
      a55_s = a55_b(l1,l2,l3,l4,l5,packb(set))
      !
    end function a55_s
    !
    !****if* src/form_factor/form_factor_5p/f55
    ! NAME
    !
    !  Function f55
    !
    ! USAGE
    !
    !  real_dim6 = f55(k1,k2,k3,k4,k5)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a55
    !
    ! INPUTS
    !
    !  * k1 -- an integer
    !  * k2 -- an integer
    !  * k3 -- an integer
    !  * k4 -- an integer
    !  * k5 -- an integer
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the global (to this module) variable b_pin_glob,b_pro_glob
    !  defined in a55
    !
    ! RETURN VALUE
    !
    !  It returns an array of six reals (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function f55(k1,k2,k3,k4,k5)
      !
      integer, intent(in) :: k1,k2,k3,k4,k5
      complex(ki), dimension(3) :: f55
      !
      integer :: j,k
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj
      integer :: b_pro_mj,b_pro_mjk
      integer :: b_pin_pj
      !
      temp1 = czero
      temp2(:) = czero
      !
      ib = b_pro_glob
      j = 0
      !
      first_pinch: do while (ib /= 0)
        !
        if (modulo(ib,2) == 1)  then
          !
          b_pro_mj = ibclr(b_pro_glob,j)
          b_pin_pj = punion( b_pin_glob,ibset(0,j) )
          !
          temp1 = temp1 + 2._ki*( 2._ki*inv_s(k4,k5,b_pin_glob)*b(j,b_pin_glob) &
                                - inv_s(j,k4,b_pin_glob)*b(k5,b_pin_glob) &
                                - inv_s(j,k5,b_pin_glob)*b(k4,b_pin_glob) &
                                - b(j,b_pin_glob)*inv_s(k4,k5,b_pin_pj) ) &
                          *f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,parf1=k1,parf2=k2,parf3=k3)
          !
          ibj = b_pro_mj
          k = 0
          !
          second_pinch: do while (ibj /= 0)
            !
            if (modulo(ibj,2) == 1) then
              !
              b_pro_mjk = ibclr(b_pro_mj,k)
              !
              temp2 = temp2 - ( inv_s(j,k5,b_pin_glob)*inv_s(k,k4,b_pin_pj) &
                               + inv_s(j,k4,b_pin_glob)*inv_s(k,k5,b_pin_pj) ) &
                              *f3p(s_mat_p,b_pro_mjk,parf1=k1,parf2=k2,parf3=k3)/4._ki       
              !
            end if
            !
            k = k+1
            ibj = ishft(ibj,-1)
            !
          end do second_pinch
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      f55 = temp2
      f55(3) = f55(3) + temp1
      !
    end function f55
    !
    !****if* src/form_factor/form_factor_5p/g55
    ! NAME
    !
    !  Function g55
    !
    ! USAGE
    !
    !  real_dim6 = g55(k1,k2,k3,k4,k5)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a55
    !
    ! INPUTS
    !
    !  * k1 -- an integer
    !  * k2 -- an integer
    !  * k3 -- an integer
    !  * k4 -- an integer
    !  * k5 -- an integer
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the global (to this module) variable b_pin_glob,b_pro_glob
    !  defined in a55
    !
    ! RETURN VALUE
    !
    !  It returns an array of six reals (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function g55(k1,k2,k3,k4,k5)
      !
      integer, intent(in) :: k1,k2,k3,k4,k5
      complex(ki), dimension(3) :: g55
      !
      integer :: j
      complex(ki) :: temp1
      integer :: ib
      integer :: b_pro_mj
      integer :: b_pin_pj
      !
      temp1 = czero
      !
      ib = b_pro_glob
      j = 0
      !
      first_pinch: do while (ib /= 0)
        !
        if (modulo(ib,2) == 1)  then
          !
          b_pro_mj = ibclr(b_pro_glob,j)
          b_pin_pj = punion( b_pin_glob,ibset(0,j) )
          !
          temp1 = temp1 + ( inv_s(j,k4,b_pin_glob)*inv_s(k3,k5,b_pin_pj) &
                           + inv_s(j,k3,b_pin_glob)*inv_s(k4,k5,b_pin_pj) &
                           + inv_s(j,k5,b_pin_glob)*inv_s(k3,k4,b_pin_pj) ) &
                          *f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,parf1=k1,parf2=k2)/2._ki
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      g55(:) = czero
      g55(3) = g55(3) + temp1
      !
    end function g55
    !
    !****f* src/form_factor/form_factor_5p/b52_b
    ! NAME
    !
    !  Function b52_b
    !
    ! USAGE
    !
    !  type(form_factor) = b52_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{5,2}.
    ! 
    ! INPUTS
    !
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function b52_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: b52_b
      !
      integer :: j
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib
      integer :: b_pro,b_pro_mj
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp1 = czero
        !
        ib = b_pro
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            b_pin_pj = punion( b_pin,ibset(0,j) )
            !
            temp1 = temp1 - b(j,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj)/2._ki
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp2(:) = czero
        temp2(3) = temp1
        b52_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function b52'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function b52_b
    !
    !****f* src/form_factor/form_factor_5p/b52_s
    ! NAME
    !
    !  Function b52_s
    !
    ! USAGE
    !
    !  type(form_factor) = b52_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{5,2}.
    ! 
    ! INPUTS
    !
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function b52_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b52_s
      !
      b52_s = b52_b(packb(set))
      !
    end function b52_s
    !
    !****f* src/form_factor/form_factor_5p/b53_b
    ! NAME
    !
    !  Function b53_b
    !
    ! USAGE
    !
    !  type(form_factor) = b53_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{5,3}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function b53_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: b53_b
      !
      integer :: j
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib
      integer :: b_pro,b_pro_mj
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp1 = czero
        !
        ib = b_pro
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            b_pin_pj = ibset(b_pin,j)
            !
            temp1 = temp1 + ( b(j,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,l1) &
                             + inv_s(j,l1,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj)/2._ki )/3._ki
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp2(:) = czero
        temp2(3) = temp1
        b53_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function b53'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function b53_b
    !
    !****f* src/form_factor/form_factor_5p/b53_s
    ! NAME
    !
    !  Function b53_s
    !
    ! USAGE
    !
    !  type(form_factor) = b53_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{5,3}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function b53_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b53_s
      !
      b53_s = b53_b(l1,packb(set))
      !
    end function b53_s
    !
    !****f* src/form_factor/form_factor_5p/b54_b
    ! NAME
    !
    !  Function b54_b
    !
    ! USAGE
    !
    !  type(form_factor) = b54_b(l1,l2,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{5,4}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function b54_b(l1,l2,b_pin)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in) :: b_pin
      type(form_factor) :: b54_b
      !
      integer :: j,k
      complex(ki) :: temp1
      complex(ki), dimension(2) :: temp2
      complex(ki), dimension(3) :: temp3
      integer :: ib,ibj
      integer :: b_pro,b_pro_mj,b_pro_mjk
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp1 = czero
        temp2(:) = czero
        !
        ib = b_pro
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            b_pin_pj = ibset(b_pin,j)
            !
            temp2 = temp2 + mult_div(-2._ki/3._ki,f4p_np4(s_mat_p,b_pro_mj,b_pin_pj)) &
                    * (inv_s(l1,l2,b_pin)*b(j,b_pin) - 0.5_ki*inv_s(j,l1,b_pin)*b(l2,b_pin) &
                                                   - 0.5_ki*inv_s(j,l2,b_pin)*b(l1,b_pin))
            !
            temp1 = temp1 - b(j,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,l1,l2)
            temp1 = temp1 - 0.5_ki*(inv_s(j,l1,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,l2) &
                  + inv_s(j,l2,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,l1))
            !
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp3(:) = czero
        temp3(2:3) = temp2
        temp3(3) = temp3(3) + temp1
        temp3 = temp3/4._ki
        b54_b = temp3
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function b54'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function b54_b
    !
    !****f* src/form_factor/form_factor_5p/b54_s
    ! NAME
    !
    !  Function b54_s
    !
    ! USAGE
    !
    !  type(form_factor) = b54_s(l1,l2,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{5,4}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function b54_s(l1,l2,set)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b54_s
      !
      b54_s = b54_b(l1,l2,packb(set))
      !
    end function b54_s
    !
    !****f* src/form_factor/form_factor_5p/b55_b
    ! NAME
    !
    !  Function b55_b
    !
    ! USAGE
    !
    !  type(form_factor) = b55_b(l1,l2.l3,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{5,5}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function b55_b(l1,l2,l3,b_pin)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in) :: b_pin
      type(form_factor) :: b55_b
      !
      integer :: j
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      complex(ki), dimension(2) :: temp3
      integer :: ib
      integer :: b_pro_mj
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        temp1 = czero
        temp3 = czero
        !
        ib = b_pro_glob
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro_glob,j)
            b_pin_pj = ibset(b_pin_glob,j)

            temp3 = temp3 + mult_div(-1._ki/2._ki,f4p_np4(s_mat_p,b_pro_mj,b_pin_pj,l1))* &
             (0.5_ki*(inv_s(j,l3,b_pin)*b(l2,b_pin) + inv_s(j,l2,b_pin)*b(l3,b_pin)) - inv_s(l2,l3,b_pin)*b(j,b_pin))
            temp3 = temp3 + mult_div(-1._ki/2._ki,f4p_np4(s_mat_p,b_pro_mj,b_pin_pj,l2))* &
             (0.5_ki*(inv_s(j,l3,b_pin)*b(l1,b_pin) + inv_s(j,l1,b_pin)*b(l3,b_pin)) - inv_s(l1,l3,b_pin)*b(j,b_pin))
            temp3 = temp3 + mult_div(-1._ki/2._ki,f4p_np4(s_mat_p,b_pro_mj,b_pin_pj,l3))* &
             (0.5_ki*(inv_s(j,l1,b_pin)*b(l2,b_pin) + inv_s(j,l2,b_pin)*b(l1,b_pin)) - inv_s(l2,l1,b_pin)*b(j,b_pin))

            !
            temp1 = temp1 + b(j,b_pin_glob)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,l1,l2,l3)
            temp1 = temp1 + 0.5_ki* inv_s(j,l3,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,l1,l2)
            temp1 = temp1 + 0.5_ki* inv_s(j,l2,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,l1,l3)
            temp1 = temp1 + 0.5_ki* inv_s(j,l1,b_pin)*f4p_np2(s_mat_p,b_pro_mj,b_pin_pj,l2,l3)
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp2(:) = czero
        temp2(2:3) =  temp3(1:2)
        temp2(3) = temp2(3) + temp1
        temp2 = temp2/5._ki
        b55_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function b55'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function b55_b
    !
    !****f* src/form_factor/form_factor_5p/b55_s
    ! NAME
    !
    !  Function b55_s
    !
    ! USAGE
    !
    !  type(form_factor) = b55_s(l1,l2.l3,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{5,5}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l2 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * l3 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function b55_s(l1,l2,l3,set)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b55_s
      !
      b55_s = b55_b(l1,l2,l3,packb(set))
      !
    end function b55_s
    !
    !****f* src/form_factor/form_factor_5p/c54_b
    ! NAME
    !
    !  Function c54_b
    !
    ! USAGE
    !
    !  type(form_factor) = c54_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor C^{5,4}.
    ! 
    ! INPUTS
    !
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat_p
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function c54_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: c54_b
      !
      integer :: j
      complex(ki), dimension(2) :: temp2
      complex(ki), dimension(3) :: temp3
      integer :: ib
      integer :: b_pro,b_pro_mj
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp2(:) = czero
        !
        ib = b_pro
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            b_pin_pj = punion( b_pin,ibset(0,j) )
            !
            temp2 = temp2 + mult_div(-2._ki/3._ki,f4p_np4(s_mat_p,b_pro_mj,b_pin_pj)) &
                           *b(j,b_pin)
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp3(:) = czero
        temp3(2:3) = temp2/4._ki
        c54_b = temp3
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function c54'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function c54_b
    !
    !****f* src/form_factor/form_factor_5p/c54_s
    ! NAME
    !
    !  Function c54_s
    !
    ! USAGE
    !
    !  type(form_factor) = c54_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor C^{5,4}.
    ! 
    ! INPUTS
    !
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function c54_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: c54_s
      !
      c54_s = c54_b(packb(set))
      !
    end function c54_s
    !
    !****f* src/form_factor/form_factor_5p/c55_b
    ! NAME
    !
    !  Function c55_b
    !
    ! USAGE
    !
    !  type(form_factor) = c55_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor C^{5,5}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables b_ref
    !  and s_mat
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function c55_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: c55_b
      !
      integer :: j
      complex(ki), dimension(2) :: temp2
      complex(ki), dimension(3) :: temp3
      integer :: ib
      integer :: b_pro,b_pro_mj
      integer :: b_pin_pj
      !
      if (dim_s >= 5) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp2(:) = czero
        !
        ib = b_pro
        j = 0
        !
        first_pinch: do while (ib /= 0)
          !
          if (modulo(ib,2) == 1)  then
            !
            b_pro_mj = ibclr(b_pro,j)
            b_pin_pj = punion( b_pin,ibset(0,j) )
            !
            temp2 = temp2 - mult_div(-1._ki/2._ki,f4p_np4(s_mat_p,b_pro_mj,b_pin_pj,l1)) &
                           *b(j,b_pin) - inv_s(j,l1,b_pin)*f4p_np4(s_mat_p,b_pro_mj,b_pin_pj)/4._ki
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp3(:) = czero
        temp3(2:3) = temp2/5._ki
        c55_b = temp3
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function c55'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 5: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function c55_b
    !
    !****f* src/form_factor/form_factor_5p/c55_s
    ! NAME
    !
    !  Function c55_s
    !
    ! USAGE
    !
    !  type(form_factor) = c55_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor C^{5,5}.
    ! 
    ! INPUTS
    !
    ! * l1 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  The result returned is of the type form_factor
    !  It returns an array of three complex (type ki) corresponding
    !  to the real part, imaginary part of the coefficient in front 1/epsilon^2,
    !  the real part, imaginary part of the 1/epsilon term and the real part,
    !  imaginary part of the constant term.
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    function c55_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: c55_s
      !
      c55_s = c55_b(l1,packb(set))
      !
    end function c55_s
    !
end module form_factor_5p
