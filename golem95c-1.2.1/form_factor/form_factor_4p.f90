!
!****h* src/form_factor/form_factor_4p
! NAME
!
!  Module form_factor_4p
!
! USAGE
!
!  use form_factor_4p
!
! DESCRIPTION
!
!  This module contains the different form factors for four point amplitudes.
!
! OUTPUT
!
!  It exports nine functions:
!  * a40 -- a function to compute A^{4,0}
!  * a41 -- a function to compute A^{4,1}
!  * a42 -- a function to compute A^{4,2}
!  * a43 -- a function to compute A^{4,3}
!  * a44 -- a function to compute A^{4,4}
!  * b42 -- a function to compute B^{4,2}
!  * b43 -- a function to compute B^{4,3}
!  * b44 -- a function to compute B^{4,4}
!  * c44 -- a function to compute C^{4,4}
!
!  Note that a4xx, b4xx and c4xx are generic functions which can be called either with a
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
!  * form_factor_type (src/module/form_factor_type.f90)
!  * constante (src/module/constante.f90)
!
!*****
!
module form_factor_4p
  !
  use precision_golem
  use generic_function_4p
  use generic_function_3p
  use matrice_s
  use array
  use sortie_erreur
  use form_factor_type
  use constante, only: czero
  !
  implicit none
  !
  private
  !
  integer :: b_pin_glob,b_pro_glob
  !
  interface a40
    !
    module procedure a40_b, a40_s
    !
  end interface
  !
  interface a41
    !
    module procedure a41_b, a41_s
    !
  end interface
  !
  interface a42
    !
    module procedure a42_b, a42_s
    !
  end interface
  !
  interface a43
    !
    module procedure a43_b, a43_s
    !
  end interface
  !
  interface a44
    !
    module procedure a44_b, a44_s
    !
  end interface
  !
  interface b42
    !
    module procedure b42_b, b42_s
    !
  end interface
  !
  interface b43
    !
    module procedure b43_b, b43_s
    !
  end interface
  !
  interface b44
    !
    module procedure b44_b, b44_s
    !
  end interface
  !
  interface c44
    !
    module procedure c44_b, c44_s
    !
  end interface
  !
  !
  public :: a40,a41,a42,a43,a44,b42,b43,b44,c44
  !
  !
  contains
    !
    !
    !****f* src/form_factor/form_factor_4p/a40_b
    ! NAME
    !
    !  Function a40_b
    !
    ! USAGE
    !
    !  type(form_factor) = a40_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,0}.
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
    function a40_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: a40_b
      !
      integer :: j
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      !real(ki) :: m1s,m2s,m3s,m4s
      integer :: ib
      integer :: b_pro,b_pro_mj
      !integer :: m1,m2,m3,m4
      integer, dimension(4) :: s
      !
      if (dim_s >= 4) then
        !
        b_pro = pminus(b_ref,b_pin)
        s = unpackb(b_pro,countb(b_pro))
        !
        !
        ! test if no internal masses are present
        !
        no_masses: if (iand(s_mat_p%b_zero, b_pro) .eq. b_pro ) then !case no-internal masses 
           !
           temp1 = sumb(b_pin)*f4p_np2(s_mat_p,b_pro,b_pin)
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
                 !
                 temp2 = temp2 + b(j,b_pin)*f3p(s_mat_p,b_pro_mj)
                 !
              end if
              !
              j = j+1
              ib= ishft(ib,-1)
              !
           end do first_pinch
           !
           temp2(3) = temp2(3) + temp1
           a40_b = temp2
           !
        else  ! internal masses are present: use 4-dim boxes 
           !
           temp2(1:3) = f4p(s_mat_p,b_pro,b_pin)
           a40_b = temp2
           !
	end if no_masses ! end test if internal masses are present
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a40'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a40_b
    !
    !****f* src/form_factor/form_factor_4p/a40_s
    ! NAME
    !
    !  Function a40_s
    !
    ! USAGE
    !
    !  type(form_factor) = a40_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,0}.
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
    function a40_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a40_s
      !
      a40_s = a40_b(packb(set))
      !
    end function a40_s
    !
    !****f* src/form_factor/form_factor_4p/a41_b
    ! NAME
    !
    !  Function a41_b
    !
    ! USAGE
    !
    !  type(form_factor) = a41_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,1}(l_1).
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
    function a41_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: a41_b
      !
      integer :: j
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib
      integer :: b_pro,b_pro_mj
      !
      if (dim_s >= 4) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp1 = -b(l1,b_pin)*f4p_np2(s_mat_p,b_pro,b_pin)
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
            !
            temp2 = temp2 - inv_s(j,l1,b_pin)*f3p(s_mat_p,b_pro_mj)
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp2(3) = temp2(3) + temp1
        a41_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a41'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a41_b
    !
    !****f* src/form_factor/form_factor_4p/a41_s
    ! NAME
    !
    !  Function a41_s
    !
    ! USAGE
    !
    !  type(form_factor) = a41_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,1}(l_1).
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
    function a41_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a41_s
      !
      a41_s = a41_b(l1,packb(set))
      !
    end function a41_s
    !
    !****f* src/form_factor/form_factor_4p/a42_b
    ! NAME
    !
    !  Function a42_b
    !
    ! USAGE
    !
    !  type(form_factor) = a42_b(l1,l2,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,2}(l1,l2).
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
    function a42_b(l1,l2,b_pin)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in) :: b_pin
      type(form_factor) :: a42_b
      !
      integer :: j
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib
      integer :: b_pro,b_pro_mj
      !
      if (dim_s >= 4) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp1 =  b(l1,b_pin)*f4p_np2(s_mat_p,b_pro,b_pin,l2) &
               + b(l2,b_pin)*f4p_np2(s_mat_p,b_pro,b_pin,l1) &
               - inv_s(l1,l2,b_pin)*f4p_np2(s_mat_p,b_pro,b_pin)
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
            !
            temp2 = temp2 + ( inv_s(j,l1,b_pin)*f3p(s_mat_p,b_pro_mj,l2) &
                            &+ inv_s(j,l2,b_pin)*f3p(s_mat_p,b_pro_mj,l1) )/2._ki
            !
          end if
          !
          j = j+1
          ib= ishft(ib,-1)
          !
        end do first_pinch
        !
        temp2(3) = temp2(3) + temp1
        a42_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a42'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a42_b
    !
    !****f* src/form_factor/form_factor_4p/a42_s
    ! NAME
    !
    !  Function a42_s
    !
    ! USAGE
    !
    !  type(form_factor) = a42_s(l1,l2,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,2}(l1,l2).
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
    function a42_s(l1,l2,set)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a42_s
      !
      a42_s = a42_b(l1,l2,packb(set))
      !
    end function a42_s
    !
    !****f* src/form_factor/form_factor_4p/a43_b
    ! NAME
    !
    !  Function a43_b
    !
    ! USAGE
    !
    !  type(form_factor) = a43_b(l1,l2,l3,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,3}(l1,l2,l3).
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
    function a43_b(l1,l2,l3,b_pin)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in) :: b_pin
      type(form_factor) :: a43_b
      !
      complex(ki), dimension(3) :: t43
      !
      if (dim_s >= 4) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        t43 = f43(l1,l2,l3) + f43(l2,l1,l3) + f43(l3,l2,l1)
        a43_b = t43
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a43'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a43_b
     !
    !****f* src/form_factor/form_factor_4p/a43_s
    ! NAME
    !
    !  Function a43_s
    !
    ! USAGE
    !
    !  type(form_factor) = a43_s(l1,l2,l3,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,3}(l1,l2,l3).
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
    function a43_s(l1,l2,l3,set)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a43_s
      !
      a43_s = a43_b(l1,l2,l3,packb(set))
      !
    end function a43_s
    !
    !****if* src/form_factor/form_factor_4p/f43
    ! NAME
    !
    !  Function f43
    !
    ! USAGE
    !
    !  real_dim6 = f43(k1,k2,k3)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a43
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
    !  defined in a43
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
    function f43(k1,k2,k3)
      !
      integer, intent(in) :: k1,k2,k3
      complex(ki), dimension(3) :: f43
      !
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: j
      integer :: ib
      integer :: b_pro_mj
      !
      temp1 = 2._ki/3._ki*inv_s(k2,k3,b_pin_glob)*f4p_np2(s_mat_p,b_pro_glob,b_pin_glob,k1) &
             - b(k1,b_pin_glob)*f4p_np2(s_mat_p,b_pro_glob,b_pin_glob,k2,k3)
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
          !
          temp2 = temp2 - inv_s(j,k1,b_pin_glob) &
                         *f3p(s_mat_p,b_pro_mj,k2,k3)/3._ki
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      f43 = temp2
      f43(3) = f43(3) + temp1
      !
    end function f43
    !
    !****f* src/form_factor/form_factor_4p/a44_b
    ! NAME
    !
    !  Function a44_b
    !
    ! USAGE
    !
    !  type(form_factor) = a44_b(l1,l2,l3,l4,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,4}(l1,l2,l3,l4).
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
    function a44_b(l1,l2,l3,l4,b_pin)
      !
      integer, intent (in) :: l1,l2,l3,l4
      integer, intent (in) :: b_pin
      type(form_factor) :: a44_b
      !
      complex(ki), dimension(3) :: t44
      !
      if (dim_s >= 4) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        t44 =  f44(l1,l2,l3,l4) + f44(l1,l3,l2,l4) + f44(l1,l4,l3,l2) &
             + f44(l3,l2,l1,l4) + f44(l4,l2,l3,l1) + f44(l3,l4,l1,l2) &
             + g44(l1,l2,l3,l4) + g44(l2,l1,l3,l4) + g44(l3,l1,l2,l4) &
             + g44(l4,l1,l2,l3)
        a44_b =  t44
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a44'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a44_b
    !
    !****f* src/form_factor/form_factor_4p/a44_s
    ! NAME
    !
    !  Function a44_s
    !
    ! USAGE
    !
    !  type(form_factor) = a44_s(l1,l2,l3,l4,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{4,4}(l1,l2,l3,l4).
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
    function a44_s(l1,l2,l3,l4,set)
      !
      integer, intent (in) :: l1,l2,l3,l4
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a44_s
      !
      a44_s = a44_b(l1,l2,l3,l4,packb(set))
      !
    end function a44_s
    !
    !****if* src/form_factor/form_factor_4p/f44
    ! NAME
    !
    !  Function f44
    !
    ! USAGE
    !
    !  real_dim6 = f44(k1,k2,k3,k4)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a44
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
    !  defined in a44
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
    function f44(k1,k2,k3,k4)
      !
      integer, intent(in) :: k1,k2,k3,k4
      complex(ki), dimension(3) :: f44
      !
      complex(ki) :: temp1
      !
      temp1 = -1._ki/2._ki*inv_s(k1,k2,b_pin_glob) &
               *f4p_np2(s_mat_p,b_pro_glob,b_pin_glob,k3,k4)
      f44(:) = czero
      f44(3) = temp1
      !
    end function f44
    !
    !****if* src/form_factor/form_factor_4p/g44
    ! NAME
    !
    !  Function g44
    !
    ! USAGE
    !
    !  real_dim6 = g44(k1,k2,k3,k4)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a44
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
    !  defined in a44
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
    function g44(k1,k2,k3,k4)
      !
      integer, intent(in) :: k1,k2,k3,k4
      complex(ki), dimension(3) :: g44
      !
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: j
      integer :: ib
      integer :: b_pro_mj
      !
      temp1 = b(k1,b_pin_glob)*f4p_np2(s_mat_p,b_pro_glob,b_pin_glob,k2,k3,k4)
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
          !
          if ( (j /= k2) .and. (j /= k3)  .and. (j /= k4) ) then
            !
            temp2 = temp2 + inv_s(j,k1,b_pin_glob) &
                           *f3p(s_mat_p,b_pro_mj,k2,k3,k4)/4._ki
            !
          end if
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      g44 = temp2
      g44(3) = g44(3) + temp1
      !
    end function g44
    !
    !****f* src/form_factor/form_factor_4p/b42_b
    ! NAME
    !
    !  Function b42_b
    !
    ! USAGE
    !
    !  type(form_factor) = b42_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{4,2}.
    ! 
    ! INPUTS
    !
    ! * set -- an array of integers of rank 1 corresponding to the label 
    !          of the propagators pinched (removed from the original set
    !          which is in the global variable set_ref)
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
    function b42_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: b42_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp2
      !
      if (dim_s >= 4) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp2(:) = czero
        temp2(3) = -f4p_np2(s_mat_p,b_pro,b_pin)/2._ki
        b42_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function b42'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function b42_b
    !
    !****f* src/form_factor/form_factor_4p/b42_s
    ! NAME
    !
    !  Function b42_s
    !
    ! USAGE
    !
    !  type(form_factor) = b42_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{4,2}.
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
    function b42_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b42_s
      !
     b42_s = b42_b(packb(set))
      !
    end function b42_s
    !
    !****f* src/form_factor/form_factor_4p/b43_b
    ! NAME
    !
    !  Function b43_b
    !
    ! USAGE
    !
    !  type(form_factor) = b43_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{4,3}.
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
    function b43_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: b43_b
      !
      integer :: b_pro
      !~ real(ki), dimension(2) :: temp1
      complex(ki), dimension(3) :: temp2
      !
      if (dim_s >= 4) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp2(:) = czero
        temp2(3) = f4p_np2(s_mat_p,b_pro,b_pin,l1)/2._ki
        b43_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function b43'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function b43_b
    !
    !****f* src/form_factor/form_factor_4p/b43_s
    ! NAME
    !
    !  Function b43_s
    !
    ! USAGE
    !
    !  type(form_factor) = b43_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{4,3}.
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
    function b43_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b43_s
      !
      b43_s = b43_b(l1,packb(set))
      !
    end function b43_s
    !
    !****f* src/form_factor/form_factor_4p/b44_b
    ! NAME
    !
    !  Function b44_b
    !
    ! USAGE
    !
    !  type(form_factor) = b44_b(l1,l2,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{4,4}.
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
    function b44_b(l1,l2,b_pin)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in) :: b_pin
      type(form_factor) :: b44_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp2
      !
      if (dim_s >= 4) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp2(:) = czero
        temp2(3) = -f4p_np2(s_mat_p,b_pro,b_pin,l1,l2)/2._ki
        b44_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function b44'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function b44_b
    !
    !****f* src/form_factor/form_factor_4p/b44_s
    ! NAME
    !
    !  Function b44_s
    !
    ! USAGE
    !
    !  type(form_factor) = b44_s(l1,l2,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor B^{4,4}.
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
    function b44_s(l1,l2,set)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: b44_s
      !
      b44_s = b44_b(l1,l2,packb(set))
      !
    end function b44_s
    !
    !****f* src/form_factor/form_factor_4p/c44_b
    ! NAME
    !
    !  Function c44_b
    !
    ! USAGE
    !
    !  type(form_factor) = c44_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor C^{4,4}.
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
    function c44_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: c44_b
      !
      integer :: b_pro
      complex(ki), dimension(3) :: temp3
      !
      if (dim_s >= 4) then
        !
        b_pro = pminus(b_ref,b_pin)
        !
        temp3(:) = czero
        !
        temp3(2:3) = f4p_np4(s_mat_p,b_pro,b_pin)/4._ki
        c44_b = temp3
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function c44'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 4: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function c44_b
    !
    !****f* src/form_factor/form_factor_4p/c44_s
    ! NAME
    !
    !  Function c44_s
    !
    ! USAGE
    !
    !  type(form_factor) = c44_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor C^{4,4}.
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
    function c44_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: c44_s
      !
     c44_s = c44_b(packb(set))
      !
    end function c44_s
    !
end module form_factor_4p
