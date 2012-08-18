!
!****h* src/form_factor/form_factor_6p
! NAME
!
!  Module form_factor_6p
!
! USAGE
!
!  use form_factor_6p
!
! DESCRIPTION
!
!  This module contains the different form factors for six point amplitudes.
!
! OUTPUT
!
!  It exports seven functions:
!  * a60 -- a function to compute A^{6,0}
!  * a61 -- a function to compute A^{6,1}
!  * a62 -- a function to compute A^{6,2}
!  * a63 -- a function to compute A^{6,3}
!  * a64 -- a function to compute A^{6,4}
!  * a65 -- a function to compute A^{6,5}
!  * a66 -- a function to compute A^{6,6}
!
!  Note that a6xx are generic functions which can be called either with a
!  set of integers or with an integer whose bits represents the set
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * generic_function_4p (src/integrals/four_point/generic_function_4p.f90)
!  * generic_function_3p (src/integrals/three_point/generic_function_3p.f90)
!  * form_factor_5p (src/form_factor/form_factor_5p.f90)
!  * array (src/module/array.f90)
!  * matrice_s (src/kinematic/matrice_s.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * form_factor_type (src/module/form_factor_type.f90)
!  * constante (src/module/constante.f90)
!
!*****
module form_factor_6p
  !
  use precision_golem
  use generic_function_4p
  use generic_function_3p
  use form_factor_5p
  use array
  use matrice_s
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
  interface a60
    !
    module procedure a60_b, a60_s
    !
  end interface
  !
  interface a61
    !
    module procedure a61_b, a61_s
    !
  end interface
  !
  interface a62
    !
    module procedure a62_b, a62_s
    !
  end interface
  !
  interface a63
    !
    module procedure a63_b, a63_s
    !
  end interface
  !
  interface a64
    !
    module procedure a64_b, a64_s
    !
  end interface
  !
  interface a65
    !
    module procedure a65_b, a65_s
    !
  end interface
  !
  interface a66
    !
    module procedure a66_b, a66_s
    !
  end interface
  !
  !
  public :: a60,a61,a62,a63,a64,a65,a66
  !
  !
  contains
    !
    !
    !****f* src/form_factor/form_factor_6p/a60_b
    ! NAME
    !
    !  Function a60_b
    !
    ! USAGE
    !
    !  type(form_factor) = a60_b(b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,0}.
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
    function a60_b(b_pin)
      !
      integer, intent (in) :: b_pin
      type(form_factor) :: a60_b
      !
      integer :: j,k,l
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj,ibjk
      integer :: b_pro,b_pro_mj,b_pro_mjk,b_pro_mjkl
      integer :: b_pin_pj,b_pin_pjk
      !
      if (dim_s >= 6) then
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
            ibj = b_pro_mj
            k = 0
            !
            second_pinch: do while (ibj /= 0)
              !
              if (modulo(ibj,2) == 1) then
                !
                b_pro_mjk = ibclr(b_pro_mj,k)
                b_pin_pjk = punion( b_pin_pj,ibset(0,k) )
                temp1 = temp1 + b(j,b_pin)*b(k,b_pin_pj)*sumb(b_pin_pjk)*f4p_np2(s_mat_p,b_pro_mjk,b_pin_pjk)
                !
                ibjk = b_pro_mjk
                l = 0
                !
                third_pinch: do while (ibjk /= 0)
                  !
                  if (modulo(ibjk,2) == 1) then
                    !
                    b_pro_mjkl = ibclr(b_pro_mjk,l)
                    temp2 = temp2 + b(j,b_pin)*b(k,b_pin_pj)*b(l,b_pin_pjk) &
                                    *f3p(s_mat_p,b_pro_mjkl)
                    !
                  end if
                  !
                  l = l + 1
                  ibjk = ishft(ibjk,-1)
                  !
                end do third_pinch
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
        a60_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a60'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 6: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a60_b
     !
    !****f* src/form_factor/form_factor_6p/a60_s
    ! NAME
    !
    !  Function a60_s
    !
    ! USAGE
    !
    !  type(form_factor) = a60_s(set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,0}.
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
    function a60_s(set)
      !
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a60_s
      !
      a60_s = a60_b(packb(set))
      !
    end function a60_s
    !
    !****f* src/form_factor/form_factor_5p/a61_b
    ! NAME
    !
    !  Function a61_b
    !
    ! USAGE
    !
    !  type(form_factor) = a61_b(l1,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,1}(l_1).
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
    function a61_b(l1,b_pin)
      !
      integer, intent (in) :: l1
      integer, intent (in) :: b_pin
      type(form_factor) :: a61_b
      !
      integer :: j,k,l
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj,ibjk
      integer :: b_pro,b_pro_mj,b_pro_mjk,b_pro_mjkl
      integer :: b_pin_pj,b_pin_pjk
      !
      if (dim_s >= 6) then
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
            ibj = b_pro_mj
            k = 0
            !
            second_pinch: do while (ibj /= 0)
              !
              if (modulo(ibj,2) == 1) then
                !
                b_pro_mjk = ibclr(b_pro_mj,k)
                b_pin_pjk = punion( b_pin_pj,ibset(0,k) )
                temp1 = temp1 - inv_s(j,l1,b_pin)*b(k,b_pin_pj)*sumb(b_pin_pjk) &
                                *f4p_np2(s_mat_p,b_pro_mjk,b_pin_pjk)
                !
                ibjk = b_pro_mjk
                l = 0
                !
                third_pinch: do while (ibjk /= 0)
                  !
                  if (modulo(ibjk,2) == 1) then
                    !
                    b_pro_mjkl = ibclr(b_pro_mjk,l)
                    temp2 = temp2 - inv_s(j,l1,b_pin)*b(k,b_pin_pj)*b(l,b_pin_pjk) &
                                    *f3p(s_mat_p,b_pro_mjkl)
                    !
                  end if
                  !
                  l = l + 1
                  ibjk = ishft(ibjk,-1)
                  !
                end do third_pinch
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
        a61_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a61'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 6: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a61_b
    !
    !****f* src/form_factor/form_factor_5p/a61_s
    ! NAME
    !
    !  Function a61_s
    !
    ! USAGE
    !
    !  type(form_factor) = a61_s(l1,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,1}(l_1).
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
    function a61_s(l1,set)
      !
      integer, intent (in) :: l1
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a61_s
      !
      a61_s = a61_b(l1,packb(set))
      !
    end function a61_s
    !
    !****f* src/form_factor/form_factor_5p/a62_b
    ! NAME
    !
    !  Function a62_b
    !
    ! USAGE
    !
    !  type(form_factor) = a62_b(l1,l2,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,2}(l1,l2).
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
    function a62_b(l1,l2,b_pin)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in) :: b_pin
      type(form_factor) :: a62_b
      !
      integer :: j,k,l
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj,ibjk
      integer :: b_pro,b_pro_mj,b_pro_mjk,b_pro_mjkl
      integer :: b_pin_pj,b_pin_pjk
      !
      if (dim_s >= 6) then
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
            ibj = b_pro_mj
            k = 0
            !
            second_pinch: do while (ibj /= 0)
              !
              if (modulo(ibj,2) == 1) then
                !
                b_pro_mjk = ibclr(b_pro_mj,k)
                b_pin_pjk = punion( b_pin_pj,ibset(0,k) )
                temp1 = temp1 + ( inv_s(j,l1,b_pin)*inv_s(k,l2,b_pin_pj) &
                                 + inv_s(j,l2,b_pin)*inv_s(k,l1,b_pin_pj) ) &
                                *sumb(b_pin_pjk)*f4p_np2(s_mat_p,b_pro_mjk,b_pin_pjk)/2._ki
                !
                ibjk = b_pro_mjk
                l = 0
                !
                third_pinch: do while (ibjk /= 0)
                  !
                  if (modulo(ibjk,2) == 1) then
                    !
                     b_pro_mjkl = ibclr(b_pro_mjk,l)
                     temp2 = temp2 + ( inv_s(j,l1,b_pin)*inv_s(k,l2,b_pin_pj) &
                                       + inv_s(j,l2,b_pin)*inv_s(k,l1,b_pin_pj) ) &
                                      *b(l,b_pin_pjk)*f3p(s_mat_p,b_pro_mjkl)/2._ki
                    !
                  end if
                  !
                  l = l + 1
                  ibjk = ishft(ibjk,-1)
                  !
                end do third_pinch
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
        a62_b = temp2
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a62'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 6: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a62_b
    !
    !****f* src/form_factor/form_factor_5p/a62_s
    ! NAME
    !
    !  Function a62_s
    !
    ! USAGE
    !
    !  type(form_factor) = a62_s(l1,l2,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,2}(l1,l2).
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
    function a62_s(l1,l2,set)
      !
      integer, intent (in) :: l1,l2
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a62_s
      !
      a62_s = a62_b(l1,l2,packb(set))
      !
   end function a62_s
    !
    !****f* src/form_factor/form_factor_5p/a63_b
    ! NAME
    !
    !  Function a63_b
    !
    ! USAGE
    !
    !  type(form_factor) = a63_b(l1,l2,l3,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,3}(l1,l2,l3).
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
    function a63_b(l1,l2,l3,b_pin)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in) :: b_pin
      type(form_factor) :: a63_b
      !
      complex(ki), dimension(3) :: t63
      !
      if (dim_s >= 6) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        t63 = f63(l1,l2,l3) + f63(l1,l3,l2) + f63(l3,l2,l1)
        a63_b = t63
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a63'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 6: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a63_b
    !
    !****f* src/form_factor/form_factor_5p/a63_s
    ! NAME
    !
    !  Function a63_s
    !
    ! USAGE
    !
    !  type(form_factor) = a63_s(l1,l2,l3,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,3}(l1,l2,l3).
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
    function a63_s(l1,l2,l3,set)
      !
      integer, intent (in) :: l1,l2,l3
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a63_s
      !
      a63_s = a63_b(l1,l2,l3,packb(set))
      !
    end function a63_s
    !
    !****if* src/form_factor/form_factor_6p/f63
    ! NAME
    !
    !  Function f63
    !
    ! USAGE
    !
    !  real_dim6 = f63(k1,k2,k3)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a63
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
    !  defined in a63
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
    function f63(k1,k2,k3)
      !
      integer, intent(in) :: k1,k2,k3
      complex(ki), dimension(3) :: f63
      !
      integer :: j,k,l
      complex(ki) :: temp1
      complex(ki), dimension(3) :: temp2
      integer :: ib,ibj,ibjk
      integer :: b_pro_mj,b_pro_mjk,b_pro_mjkl
      integer :: b_pin_pj,b_pin_pjk
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
          ibj = b_pro_mj
          k = 0
          !
          second_pinch: do while (ibj /= 0)
            !
            if (modulo(ibj,2) == 1) then
              !
              b_pro_mjk = ibclr(b_pro_mj,k)
              b_pin_pjk = punion( b_pin_pj,ibset(0,k) )
              temp1 = temp1 - inv_s(k1,k2,b_pin_glob)*inv_s(j,k3,b_pin_glob)*b(k,b_pin_pj) &
                              *f4p_np2(s_mat_p,b_pro_mjk,b_pin_pjk)/3._ki &
                            - inv_s(j,k3,b_pin_glob)*(  inv_s(k,k1,b_pin_pj)*b(k2,b_pin_pj) &
                                               + inv_s(k,k2,b_pin_pj)*b(k1,b_pin_pj) &
                                               - 2._ki*inv_s(k1,k2,b_pin_pj)*b(k,b_pin_pj) &
                                               + b(k,b_pin_pj)*inv_s(k1,k2,b_pin_pjk) ) &
                                              *f4p_np2(s_mat_p,b_pro_mjk,b_pin_pjk)/3._ki
              !
              ibjk = b_pro_mjk
              l = 0
              !
              third_pinch: do while (ibjk /= 0)
                !
                if (modulo(ibjk,2) == 1) then
                  !
                  b_pro_mjkl = ibclr(b_pro_mjk,l)
                  temp2 = temp2 - inv_s(j,k3,b_pin_glob)*( &
                                     inv_s(k,k2,b_pin_pj)*inv_s(l,k1,b_pin_pjk) &
                                   + inv_s(k,k1,b_pin_pj)*inv_s(l,k2,b_pin_pjk) ) &
                                  *f3p(s_mat_p,b_pro_mjkl)/6._ki
                  !
                end if
                !
                l = l + 1
                ibjk = ishft(ibjk,-1)
                !
              end do third_pinch
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
      f63 = temp2
      f63(3) = f63(3) + temp1
      !
    end function f63
    !
    !****f* src/form_factor/form_factor_6p/a64_b
    ! NAME
    !
    !  Function a64_b
    !
    ! USAGE
    !
    !  type(form_factor) = a64_b(l1,l2,l3,l4,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,4}(l1,l2,l3,l4).
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
    function a64_b(l1,l2,l3,l4,b_pin)
      !
      integer, intent (in) :: l1,l2,l3,l4
      integer, intent (in) :: b_pin
      type(form_factor) :: a64_b
      !
      if (dim_s >= 6) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        a64_b=  - (  f64(l1,l2,l3,l4) + f64(l2,l1,l3,l4) + f64(l3,l2,l1,l4) &
                  + f64(l4,l2,l3,l1) )/4._ki
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a64'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 6: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a64_b
    !
    !****f* src/form_factor/form_factor_6p/a64_s
    ! NAME
    !
    !  Function a64_s
    !
    ! USAGE
    !
    !  type(form_factor) = a64_s(l1,l2,l3,l4,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,4}(l1,l2,l3,l4).
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
    function a64_s(l1,l2,l3,l4,set)
      !
      integer, intent (in) :: l1,l2,l3,l4
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a64_s
      !
      a64_s = a64_b(l1,l2,l3,l4,packb(set))
      !
    end function a64_s
    !
    !****if* src/form_factor/form_factor_6p/f64
    ! NAME
    !
    !  Function f64
    !
    ! USAGE
    !
    !  real_dim6 = f64(k1,k2,k3,k4)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a64
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
    !  No side effect, it uses the global (to this module) variable b_pro_glob, b_pin_glob
    !  defined in a64
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
    function f64(k1,k2,k3,k4)
      !
      integer, intent(in) :: k1,k2,k3,k4
      type(form_factor) :: f64
      !
      integer :: j
      type(form_factor) :: temp1,temp2
      integer :: ib
      integer :: b_pin_pj
      !
      temp1 = czero
      temp2 = czero
      !
      ib = b_pro_glob
      j = 0
      !
      first_pinch: do while (ib /= 0)
        !
        if (modulo(ib,2) == 1)  then
          !
          b_pin_pj = punion( b_pin_glob,ibset(0,j) )
          !
          temp1 = temp1 - 2._ki*inv_s(j,k1,b_pin_glob)*(  inv_s(k2,k3,b_pin_glob)*b53(k4,b_pin_pj) &
                                                 + inv_s(k2,k4,b_pin_glob)*b53(k3,b_pin_pj) &
                                                 + inv_s(k3,k4,b_pin_glob)*b53(k2,b_pin_pj) )
          temp2 = temp2 + inv_s(j,k1,b_pin_glob)*a53(k2,k3,k4,b_pin_pj)
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      f64 = temp1 + temp2
      !
    end function f64
    !
    !****f* src/form_factor/form_factor_6p/a65_b
    ! NAME
    !
    !  Function a65_b
    !
    ! USAGE
    !
    !  type(form_factor) = a65_b(l1,l2,l3,l4,l5,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,5}(l1,l2,l3,l4,l5).
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
    function a65_b(l1,l2,l3,l4,l5,b_pin)
      !
      integer, intent (in) :: l1,l2,l3,l4,l5
      integer, intent (in) :: b_pin
      type(form_factor) :: a65_b
      !
      if (dim_s >= 6) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        a65_b =  - ( f65(l1,l2,l3,l4,l5) + f65(l2,l1,l3,l4,l5) + f65(l3,l2,l1,l4,l5) &
                 + f65(l4,l2,l3,l1,l5) + f65(l5,l2,l3,l4,l1) )/5._ki
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a65'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 6: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a65_b
    !
    !****f* src/form_factor/form_factor_6p/a65_s
    ! NAME
    !
    !  Function a65_s
    !
    ! USAGE
    !
    !  type(form_factor) = a65_s(l1,l2,l3,l4,l5,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,5}(l1,l2,l3,l4,l5).
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
    function a65_s(l1,l2,l3,l4,l5,set)
      !
      integer, intent (in) :: l1,l2,l3,l4,l5
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a65_s
      !
      a65_s = a65_b(l1,l2,l3,l4,l5,packb(set))
      !
    end function a65_s
    !
    !****if* src/form_factor/form_factor_6p/f65
    ! NAME
    !
    !  Function f65
    !
    ! USAGE
    !
    !  real_dim6 = f65(k1,k2,k3,k4,k5)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a65
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
    !  No side effect, it uses the global (to this module) variable b_pro_glob,b_pin_glob
    !  defined in a65
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
    function f65(k1,k2,k3,k4,k5)
      !
      integer, intent(in) :: k1,k2,k3,k4,k5
      type(form_factor) :: f65
      !
      integer :: j
      type(form_factor) :: temp1,temp2
      integer :: ib
      integer :: b_pin_pj
      !
      temp1 = czero
      temp2 = czero
      !
      ib = b_pro_glob
      j = 0
      !
      first_pinch: do while (ib /= 0)
         !
        if (modulo(ib,2) == 1)  then
          !
          b_pin_pj = punion( b_pin_glob,ibset(0,j) )
          !
          temp1 = temp1 + inv_s(j,k1,b_pin_glob)*( &
                             -2._ki*inv_s(k2,k3,b_pin_glob)*b54(k4,k5,b_pin_pj) &
                             -2._ki*inv_s(k2,k4,b_pin_glob)*b54(k3,k5,b_pin_pj) &
                             -2._ki*inv_s(k2,k5,b_pin_glob)*b54(k3,k4,b_pin_pj) &
                             -2._ki*inv_s(k3,k4,b_pin_glob)*b54(k2,k5,b_pin_pj) &
                             -2._ki*inv_s(k3,k5,b_pin_glob)*b54(k2,k4,b_pin_pj) &
                             -2._ki*inv_s(k4,k5,b_pin_glob)*b54(k2,k3,b_pin_pj) )
          temp1 = temp1 + inv_s(j,k1,b_pin_glob)*( &			   
                             +4._ki*inv_s(k2,k3,b_pin_glob)*inv_s(k4,k5,b_pin_glob) &
                             +4._ki*inv_s(k2,k4,b_pin_glob)*inv_s(k3,k5,b_pin_glob) &
                             +4._ki*inv_s(k2,k5,b_pin_glob)*inv_s(k3,k4,b_pin_glob) ) &
                             *c54(b_pin_pj)
          temp2 = temp2 + inv_s(j,k1,b_pin_glob)*a54(k2,k3,k4,k5,b_pin_pj)
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      f65 = temp1 + temp2
      !
    end function f65
    !
    !****f* src/form_factor/form_factor_6p/a66_b
    ! NAME
    !
    !  Function a66_b
    !
    ! USAGE
    !
    !  type(form_factor) = a66_b(l1,l2,l3,l4,l5,l6,b_pin)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,6}(l1,l2,l3,l4,l5,l6).
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
    ! * l6 -- an integer corresponding to a label in the set of the three 
    !         remaining propagators
    ! * b_pin -- an integer whose bits represent an array of integers of rank 1 corresponding 
    !            to the label of the propagators pinched (removed from the original set
    !            which is in the global variable b_ref)
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the value of the global variables set_ref
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
    function a66_b(l1,l2,l3,l4,l5,l6,b_pin)
      !
      integer, intent (in) :: l1,l2,l3,l4,l5,l6
      integer, intent (in) :: b_pin
      type(form_factor) :: a66_b
      !
      if (dim_s >= 6) then
        !
        b_pro_glob = pminus(b_ref,b_pin)
        b_pin_glob = b_pin
        !
        a66_b =  - (  f66(l1,l2,l3,l4,l5,l6) + f66(l2,l1,l3,l4,l5,l6) &
                  + f66(l3,l2,l1,l4,l5,l6) + f66(l4,l2,l3,l1,l5,l6) &
                  + f66(l5,l2,l3,l4,l1,l6) + f66(l6,l2,l3,l4,l5,l1) )/6._ki
        !
      else
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In function a66'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the dimension of dim_s is less than 6: %d0'
        tab_erreur_par(2)%arg_int = dim_s
        call catch_exception(0)
        !
      end if
      !
    end function a66_b
    !
    !****f* src/form_factor/form_factor_6p/a66_s
    ! NAME
    !
    !  Function a66_s
    !
    ! USAGE
    !
    !  type(form_factor) = a66_s(l1,l2,l3,l4,l5,l6,set)
    !
    ! DESCRIPTION
    !
    ! This function defines the form factor A^{6,6}(l1,l2,l3,l4,l5,l6).
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
    ! * l6 -- an integer corresponding to a label in the set of the three 
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
    function a66_s(l1,l2,l3,l4,l5,l6,set)
      !
      integer, intent (in) :: l1,l2,l3,l4,l5,l6
      integer, intent (in), dimension(:) :: set
      type(form_factor) :: a66_s
      !
      a66_s = a66_b(l1,l2,l3,l4,l5,l6,packb(set))
      !
    end function a66_s
    !
    !****if* src/form_factor/form_factor_6p/f66
    ! NAME
    !
    !  Function f66
    !
    ! USAGE
    !
    !  real_dim6 = f66(k1,k2,k3,k4,k5,k6)
    !
    ! DESCRIPTION
    !
    !  A function to simplify the writting of the function a66
    !
    ! INPUTS
    !
    !  * k1 -- an integer
    !  * k2 -- an integer
    !  * k3 -- an integer
    !  * k4 -- an integer
    !  * k5 -- an integer
    !  * k6 -- an integer
    !
    ! SIDE EFFECTS
    !
    !  No side effect, it uses the global (to this module) variable b_pro_glob,b_pin_glob
    !  defined in a66
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
    function f66(k1,k2,k3,k4,k5,k6)
      !
      integer, intent(in) :: k1,k2,k3,k4,k5,k6
      type(form_factor) :: f66
      !
      integer :: j
      type(form_factor) :: temp1,temp2
      integer :: ib
      integer :: b_pin_pj
      !
      temp1 = czero
      temp2 = czero
      !
      ib = b_pro_glob
      j = 0
      !
      first_pinch: do while (ib /= 0)
        !
        if (modulo(ib,2) == 1)  then
          !
          b_pin_pj = punion( b_pin_glob,ibset(0,j) )
          !
          temp1 = temp1 + inv_s(j,k1,b_pin_glob)*( &
                             -2._ki*inv_s(k2,k3,b_pin_glob)*b55(k4,k5,k6,b_pin_pj) &
                             -2._ki*inv_s(k2,k4,b_pin_glob)*b55(k3,k5,k6,b_pin_pj) &
                             -2._ki*inv_s(k2,k5,b_pin_glob)*b55(k3,k4,k6,b_pin_pj) &
                             -2._ki*inv_s(k2,k6,b_pin_glob)*b55(k3,k4,k5,b_pin_pj) &
                             -2._ki*inv_s(k3,k4,b_pin_glob)*b55(k2,k5,k6,b_pin_pj) &
                             -2._ki*inv_s(k3,k5,b_pin_glob)*b55(k2,k4,k6,b_pin_pj) &
                             -2._ki*inv_s(k3,k6,b_pin_glob)*b55(k2,k4,k5,b_pin_pj) &
                             -2._ki*inv_s(k4,k5,b_pin_glob)*b55(k2,k3,k6,b_pin_pj) &
                             -2._ki*inv_s(k4,k6,b_pin_glob)*b55(k2,k3,k5,b_pin_pj) &
                             -2._ki*inv_s(k5,k6,b_pin_glob)*b55(k2,k3,k4,b_pin_pj) )
          !
          temp1 = temp1 + inv_s(j,k1,b_pin_glob)*( &			   
                             ( 4._ki*inv_s(k2,k3,b_pin_glob)*inv_s(k4,k5,b_pin_glob) &
                              +4._ki*inv_s(k2,k4,b_pin_glob)*inv_s(k3,k5,b_pin_glob) &
                              +4._ki*inv_s(k2,k5,b_pin_glob)*inv_s(k3,k4,b_pin_glob) ) &
                             *c55(k6,b_pin_pj) &
                             + ( 4._ki*inv_s(k2,k3,b_pin_glob)*inv_s(k4,k6,b_pin_glob) &
                                +4._ki*inv_s(k2,k4,b_pin_glob)*inv_s(k3,k6,b_pin_glob) &
                                +4._ki*inv_s(k2,k6,b_pin_glob)*inv_s(k3,k4,b_pin_glob) ) &
                             *c55(k5,b_pin_pj) &
                             + ( 4._ki*inv_s(k2,k3,b_pin_glob)*inv_s(k5,k6,b_pin_glob) &
                                +4._ki*inv_s(k2,k6,b_pin_glob)*inv_s(k3,k5,b_pin_glob) &
                                +4._ki*inv_s(k2,k5,b_pin_glob)*inv_s(k3,k6,b_pin_glob) ) &
                             *c55(k4,b_pin_pj) &
                             + ( 4._ki*inv_s(k2,k6,b_pin_glob)*inv_s(k4,k5,b_pin_glob) &
                                +4._ki*inv_s(k2,k4,b_pin_glob)*inv_s(k6,k5,b_pin_glob) &
                                +4._ki*inv_s(k2,k5,b_pin_glob)*inv_s(k6,k4,b_pin_glob) ) &
                             *c55(k3,b_pin_pj) &
                             + ( 4._ki*inv_s(k6,k3,b_pin_glob)*inv_s(k4,k5,b_pin_glob) &
                                +4._ki*inv_s(k6,k4,b_pin_glob)*inv_s(k3,k5,b_pin_glob) &
                                +4._ki*inv_s(k6,k5,b_pin_glob)*inv_s(k3,k4,b_pin_glob) ) &
                             *c55(k2,b_pin_pj) )
          !
          temp2 = temp2 + inv_s(j,k1,b_pin_glob)*a55(k2,k3,k4,k5,k6,b_pin_pj)
          !
        end if
        !
        j = j+1
        ib= ishft(ib,-1)
        !
      end do first_pinch
      !
      f66 = temp1 + temp2
      !
    end function f66
    !
end module form_factor_6p
