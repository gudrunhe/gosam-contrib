!****h* src/numerical/adapt_gauss
! NAME
!
!  Module adapt_gauss (file src/numerical/adapt_gauss.f90)
!
! USAGE
!
!  use adapt_gauss
!
! DESCRIPTION
!
!  This module contains several routines for a one dimensional
!  integration using Gauss Kronrod method
!
! OUTPUT
!
!  The only subroutine which can be used by use association in adapt_gaus1,
!  all the other subroutines/functions of this module are private
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * array (src/module/array.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!
!*****
!
module adapt_gauss
  !
  use precision_golem
  use array, only : packb
  use sortie_erreur, only : tab_erreur_par,catch_exception
  implicit none
  !
  private :: ki
  real(ki) :: tol_glob
  complex(ki) :: err_glob,res_glob
  logical :: encore_glob
  integer :: compt_call_glob,compt_cell_glob
  character (len=3), parameter :: name_par = 'kro'
  !~ integer, parameter :: n0_par = 10, n1_par = 21
  integer, parameter :: n0_par = 7, n1_par = 15
  integer, parameter :: compt_max_par = 40 ! maximal number of iteration
  real(ki), parameter :: tol_max_par = 1.e-12_ki
  integer, parameter :: nb_cell_max_par = 100000 ! maximal number of cells 
                                                 !(limited by the memory requirement)
  public :: adapt_gauss1
  !
  type intervalle
    real(ki),dimension(1) :: point
    real(ki) :: taille
    complex(ki) :: resultat
    complex(ki) :: erreur
    logical :: divise
    type(intervalle), pointer :: suivant
  end type intervalle
  !
  contains
    !
    !****f* src/numerical/adapt_gauss/adapt_gauss1
    ! NAME
    !
    !  Subroutine adapt_gauss1 (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call adapt_gauss1(func,b_inf,b_sup,tol,rest,abserr)
    !
    ! DESCRIPTION
    !
    ! This subroutine performs a one dimensional adaptative integration of the 
    ! function func between b_inf and b_sup using a Gaussian quadrature
    ! with Kronrod polynomial. The integrand is assumed to be complex.  For a
    ! certain criterium (function test_error), the range of integration is split 
    ! into two. All the cells are put into a chained list whose element are of 
    ! type intervalle 
    !
    ! INPUTS
    !
    !  this subroutine takes as inputs:
    !  * func -- an external function as declared by the interface block
    !  * b_inf -- a real (type ki), the lower bound of the integration range
    !  * b_sup -- a real (type ki), the upper bound of the integration range
    !  * tol -- a real (type ki), the tolerance asked by the user
    !
    ! SIDE EFFECTS
    !
    !  no side effects
    !
    ! RETURN VALUE
    !
    !  it returns:
    !  * rest -- a complex (type ki), the result of the integration
    !  * abserr -- a complex (type ki), the absolute value of the estimated error
    !
    ! EXAMPLE
    !
    !  to integrate a function f
    !  between 0 and 1 with a tolerance of 0.0001
    !  the result is put in result
    !  and the relative error returned in error 
    !
    !  call adapt_gauss1(f,0._ki,1._ki,1.e-4_ki,result,error)
    !
    !*****
    !
    subroutine adapt_gauss1(func,b_inf,b_sup,tol,rest,abserr)
      !
      real(ki), intent (in) :: b_inf,b_sup,tol
      complex(ki), intent (out) :: rest
      complex(ki), intent (out) :: abserr
      interface
        function func(x)
          use precision_golem
          real(ki), intent (in) :: x
          complex(ki) :: func
        end function func
      end interface
      !
      type(intervalle), pointer :: new
      complex(ki) :: rest1,abserr1
      integer :: compt
      !
      ! initialisation
      !
      rest = 0._ki
      compt = 0
      compt_call_glob = 0
      compt_cell_glob = 0
      tol_glob = 1._ki*tol
      res_glob = 0._ki
      err_glob = 0._ki
      encore_glob = .true.
      !
      ! first evaluation into the entire range
      !
      call gauss1(func,b_inf,b_sup,rest1,abserr1)
      !
      compt_call_glob = compt_call_glob + n1_par
      err_glob = abserr1
      compt = compt + 1
      compt_cell_glob = compt_cell_glob + 1
      res_glob = rest1
      !
      if ( test_error(err_glob,tol_glob) ) then
        !
        rest = res_glob
        abserr = err_glob
        !
      else ! one divides
        !
        call creation(new,b_inf,b_sup) ! creation of the nested list
        res_glob = 0._ki
        err_glob = 0._ki
        ! call imprime(new)
        !
        do while ( encore_glob .and. (compt < compt_max_par) )
          !
          encore_glob = .false.
          !
          if (tol_glob >= tol_max_par) then
            !
            tol_glob = tol_glob/2._ki
            !
          else
            !
            tol_glob = tol_glob
            !
          end if
          !
          call decoupe(new,func) ! each cells marked TRUE are split in two
          ! call imprime(new)
          !
          if (compt == (compt_max_par-1) ) then
            !
            call recupere_total(new) ! result of all remaining cells is collected
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'finished because the number of iterations reachs its maximum'
            call catch_exception(1)
            !
          else if (compt_cell_glob > nb_cell_max_par) then
            !
            call recupere_total(new) ! result of all remaining cells is collected
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'finished because the number of cells reachs its maximum'
            call catch_exception(1)
            !
          else
            !
            ! the results of FALSE cells are collected and the cells destroyed
            call recupere_partiel(new) 
            !
          end if
          ! call imprime(new)
          !
          compt = compt + 1
          !
        end do
        !
        rest = res_glob
        abserr = err_glob
        !
        call libere(new)
        !
      end if
      !
      tab_erreur_par(1)%a_imprimer = .true.
      tab_erreur_par(1)%chaine = 'Statistic in subroutine adapt_gauss1'
      tab_erreur_par(2)%a_imprimer = .true.
      tab_erreur_par(2)%chaine = 'number of function calls: %d0'
      tab_erreur_par(2)%arg_int = compt_call_glob
      tab_erreur_par(3)%a_imprimer = .true.
      tab_erreur_par(3)%chaine = 'number of cells: %d0'
      tab_erreur_par(3)%arg_int = compt_cell_glob
      tab_erreur_par(4)%a_imprimer = .true.
      tab_erreur_par(4)%chaine = 'number of iteration: %d0'
      tab_erreur_par(4)%arg_int = compt
      tab_erreur_par(5)%a_imprimer = .true.
      tab_erreur_par(5)%chaine = 'number of Gauss points: %d1'
      tab_erreur_par(5)%arg_int_tab = packb( (/n0_par,n1_par/) )
      tab_erreur_par(6)%a_imprimer = .true.
      tab_erreur_par(6)%chaine = 'Type of polynom: %c0'
      tab_erreur_par(6)%arg_char = name_par
      tab_erreur_par(7)%a_imprimer = .true.
      tab_erreur_par(7)%chaine = 'Tolerance: %f0'
      tab_erreur_par(7)%arg_real = tol
      call catch_exception(2)
      !
    end subroutine adapt_gauss1
    !
    !****if* src/numerical/adapt_gauss/relative_error
    ! NAME
    !
    !  function test_error (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  logical = test_error(error,tol)
    !
    ! DESCRIPTION
    !
    ! This function tests the error compared to the tolerance required, it
    ! returns a value of type logical
    !
    ! INPUTS
    !
    !  this function take two arguments:
    !  * error -- a complex (type ki), the absolute error
    !  * tol -- the required tolerance
    !
    ! SIDE EFFECTS
    !
    !  no side effects
    !
    ! RETURN VALUE
    !
    !  test_error is a logical
    !
    ! EXAMPLE
    !
    !*****
    !
    function test_error(error,tol)
      !
      complex(ki), intent(in) :: error
      real(ki), intent(in) :: tol
      !
      logical :: test_error
      real(ki) :: tempa,tempb
      !
      tempa = real(error,ki)
      tempb = aimag(error)
      test_error = ( (abs(tempa) <= tol) .and. (abs(tempb) <= tol) ) .or. &
                  & ( (abs(tempa) <= tol_max_par) .and. (abs(tempb) <= tol_max_par) )
      !
    end function test_error
    !
    !****if* src/numerical/adapt_gauss/creation
    ! NAME
    !
    !  subroutine creation (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call creation(new,b_inf,b_sup)
    !
    ! DESCRIPTION
    !
    !  This routine creates a chained list new whose elments are of type
    !  intervalle. The pointer of the last element of the list must be null
    !  (i.e. points on nothing). Norma f95 is used. Variables ending by _glob are
    !  global for this module
    !
    !
    ! INPUTS
    !
    !  * new -- an intervalle type, a pointer on the chained list
    !  * b_inf -- a real (type ki), lower value of the list
    !  * b_sup -- a real (type ki), upper value of the list
    !
    ! SIDE EFFECTS
    !
    !  no side effects
    !
    ! RETURN VALUE
    !
    ! after the call, new is pointer on a chained list 
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine creation(new,b_inf,b_sup)
      !
      type(intervalle), pointer :: new
      real(ki), intent (in) :: b_inf,b_sup
      !
      type(intervalle), pointer :: init,fin
      integer :: res
      !
      allocate(fin,stat=res)
      !
      if (res /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine creation (module numerical_evaluation)'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'enable to allocate fin %d0'
        tab_erreur_par(2)%arg_int = res
        call catch_exception(0)
        !
      end if
      !
      allocate(init,stat=res)
      !
      if (res /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine creation (module numerical_evaluation)'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'enable to allocate init %d0'
        tab_erreur_par(2)%arg_int = res
        call catch_exception(0)
        !
      end if
      !
      fin%point = (/b_sup/)
      fin%taille = 0._ki
      fin%resultat = (0._ki,0._ki)
      fin%erreur = 0._ki
      fin%divise = .false.
      fin%suivant => null()
      new => fin
      !
      init%point = (/b_inf/)
      init%taille = b_sup - b_inf
      init%resultat = res_glob
      init%erreur = err_glob
      init%divise = .true.
      init%suivant => fin
      new => init
      !
    end subroutine creation
    ! 
    !****if* src/numerical/adapt_gauss/decoupe
    ! NAME
    !
    !  recursive subroutine decoupe (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call decoupe(new,func)
    !
    ! DESCRIPTION
    !
    !  For each cell of the chained list, this routine splits a 1-dimension cell 
    !  into 2 subcells if the cell is marked true For the two sub-cells, it computes 
    !  the integral in the sub-cell and marks it true of false depending the error 
    !  returned. Note that since this subroutine is recursive, it acts
    !  globaly on the whole chained list containing the cells
    !
    !
    ! INPUTS
    !
    !  * new -- a pointer (type intervalle) on the chained list
    !  * func -- an external function R --> C, the integrand
    !
    ! SIDE EFFECTS
    !
    !  this routine modify the chained list whose the first pointer is new
    !
    ! RETURN VALUE
    !
    !  no return value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    recursive subroutine decoupe(new,func)
      !
      type(intervalle), pointer :: new
      interface
        function func(x)
          use precision_golem
          real(ki), intent (in) :: x
          complex(ki) :: func
        end function func
      end interface
      !
      type(intervalle), pointer :: nouveau => null()
      real(ki), dimension(1) :: vx
      real(ki) :: n_taille
      complex(ki) :: rest1_loc,abserr1_loc
      logical :: div
      integer :: res
      !
      if (associated(new%suivant)) then ! if this is not the last cell
        !
        if (new%divise) then ! if one has to split
          !
          n_taille = new%taille/2._ki
          vx = (/n_taille/)
          new%taille = n_taille
          !
          call gauss1(func,new%point(1),new%point(1)+new%taille,&
                      &rest1_loc,abserr1_loc)
          !
          compt_call_glob = compt_call_glob + n1_par
          new%erreur = abserr1_loc
          compt_cell_glob = compt_cell_glob + 1
          new%resultat = rest1_loc
          !
          if ( test_error(new%erreur,tol_glob) ) then
            !
            div = .false.
            !
          else
            !
            div = .true.
            !
          end if
          !
          new%divise = div
          allocate(nouveau,stat=res)
          !
          if (res /= 0) then
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In subroutine decoupe (module numerical_evaluation)'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = 'the allocation runs into trouble %d0'
            tab_erreur_par(2)%arg_int = res
            call catch_exception(0)
            !
          end if
          !
          nouveau%point = new%point + vx
          nouveau%taille = n_taille
          call gauss1(func,nouveau%point(1),nouveau%point(1)+nouveau%taille,&
                      &rest1_loc,abserr1_loc)
          compt_call_glob = compt_call_glob + n1_par
          nouveau%erreur = abserr1_loc
          nouveau%resultat = rest1_loc
          !
          if ( test_error(nouveau%erreur,tol_glob) ) then
            !
            div = .false.
            !
          else
            !
            div = .true.
            !
          end if
          !
          nouveau%divise = div
          nouveau%suivant => new%suivant
          new%suivant => nouveau
          compt_cell_glob = compt_cell_glob + 1
          call decoupe(nouveau%suivant,func)
          !
        else
          !
          call decoupe(new%suivant,func)
          !
        end if
        !
      end if
      !
    end subroutine decoupe
    !
    !****if* src/numerical/adapt_gauss/recupere_total
    ! NAME
    !
    !  recursive subroutine recupere_total (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call recupere_total(new)
    !
    ! DESCRIPTION
    !
    !  This routine collects the value of the integrant in each cell 
    !  (as well as the error). It acts globaly on the chained list
    !
    !
    ! INPUTS
    !
    !  * new -- a pointer on the first entry of the chained list
    !
    ! SIDE EFFECTS
    !
    !  this routine modifies the global (for this module) variables:
    !  * res_glob -- sum of the result of the integral for each cell
    !  * err_glob -- sum of the error for each cell
    !
    ! RETURN VALUE
    !
    !  no returned value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    recursive subroutine recupere_total(new)
      !
      type(intervalle), pointer :: new
      !
      if (associated(new%suivant)) then
        !
        res_glob = res_glob + new%resultat
        err_glob = err_glob + new%erreur
        !
        call recupere_total(new%suivant)
        !
      end if
      !
    end subroutine recupere_total 
    !
    !****if* src/numerical/adapt_gauss/recupere_partiel
    ! NAME
    !
    !  recursive subroutine recupere_partiel (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call recupere_partiel(new)
    !
    ! DESCRIPTION
    !
    !  This routine collects the value of the integrant only in cell marked FALSE 
    !  (as well as the error) and removes these cells to save space
    !
    !
    ! INPUTS
    !
    !  * new -- a pointer (type inetrvalle) pointing on the first entry of the chained list
    !
    ! SIDE EFFECTS
    !
    !  this routine modifies the global (for this module) variables:
    !  * res_glob -- sum of the result of the integral for each cell
    !  * err_glob -- sum of the error for each cell
    !  * encore_glob -- a logical, to know if there are cells marked true
    !
    ! RETURN VALUE
    !
    !  no returned value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    recursive subroutine recupere_partiel(new)
      !
      type(intervalle), pointer :: new
      !
      type(intervalle), pointer :: temp
      type(intervalle), pointer :: temp1
      integer :: res
      !
      if (associated(new%suivant)) then
        !
        call recupere_partiel(new%suivant)
        encore_glob = encore_glob .or. new%divise
        !
        if (.not.new%divise) then
          !
          res_glob = res_glob + new%resultat
          err_glob = err_glob + new%erreur
          temp => new
          temp1 => new%suivant
          new => temp1
          !
          deallocate(temp,stat=res)
          !
          if (res /= 0) then
            !
            tab_erreur_par(1)%a_imprimer = .true.
            tab_erreur_par(1)%chaine = 'In subroutine recupere_partiel (module numerical_evaluation)'
            tab_erreur_par(2)%a_imprimer = .true.
            tab_erreur_par(2)%chaine = 'the deallocation runs into trouble %d0'
            tab_erreur_par(2)%arg_int = res
            call catch_exception(0)
            !
          end if
          !
        end if
        !
      end if
      !
    end subroutine recupere_partiel 
    !
    !****if* src/numerical/adapt_gauss/imprime
    ! NAME
    !
    !  recursive subroutine imprime (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call imprime(new)
    !
    ! DESCRIPTION
    !
    !  This routine prints the value of the element of the chained list 
    !  i.e. the structure of each cell
    !
    !
    ! INPUTS
    !
    !  * new -- a pointer (type inetrvalle) pointing on the first entry of the chained list
    !
    ! SIDE EFFECTS
    !
    !  no side effects
    !
    ! RETURN VALUE
    !
    !  no returned value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    recursive subroutine imprime(new)
      !
      type(intervalle), pointer :: new
      !
      if (associated(new%suivant)) then
        !
        write (*,'("borne inf",1f16.8,2x,"taille ",f16.8,2x,&
          &" on divise? ",l2)')&
          & new%point,new%taille,new%divise
        write (*,'("resultat",2f16.8,2x,"erreur ",2e16.8)')&
          & new%resultat,new%erreur
        !
        call imprime(new%suivant)
        !
      end if
      !
    end subroutine imprime 
    !
    !****if* src/numerical/adapt_gauss/libere
    ! NAME
    !
    !  recursive subroutine libere (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call libere(new)
    !
    ! DESCRIPTION
    !
    !  This routine deallocates all the elements of the chained list
    !
    ! INPUTS
    !
    !  * new -- a pointer (type inetrvalle) pointing on the first entry of the chained list
    !
    ! SIDE EFFECTS
    !
    !  destroy the chained list
    !
    ! RETURN VALUE
    !
    !  no returned value
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    recursive subroutine libere(new)
      !
      type(intervalle), pointer :: new
      !
      integer :: res
      !
      if (associated(new%suivant)) call libere(new%suivant)
      !
      deallocate(new,stat=res)
      !
      if (res /= 0) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'In subroutine libere (module numerical_evaluation)'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the deallocation runs into trouble %d0'
        tab_erreur_par(2)%arg_int = res
        call catch_exception(0)
        !
      end if
      !
    end subroutine libere
    !
    !****if* src/numerical/adapt_gauss/give_me_the_weight
    ! NAME
    !
    !  subroutine give_me_the_weight (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call give_me_the_weight(nomb,name,weight,zero)
    !
    ! DESCRIPTION
    !
    !  This routine gives the weight and zero for Legendre polynomials of
    !  degrees : n=7,10 (created with maple file poid_gauss.m)
    !  and also the weight and zero for Kronrod polynomials of degrees 15,21
    !  in such a way the zeros of the Legendre polynomial of degree 7 (resp. 10)
    !  is between the zeros of the Kronrod polynomial of degree 15 (resp. 21).
    !
    ! INPUTS
    !
    !  * nomb -- an integer, the degree of the Legendre/Kronrod polynomials
    !  * name -- a character (dimension 3), the name of the polynom
    !
    ! SIDE EFFECTS
    !
    !  no side effects
    !
    ! RETURN VALUE
    !
    !  * weight -- an array of real (type ki) containing the weights
    !  * zero -- an array of real (type ki) containing the zero of the polynomials
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine give_me_the_weight(nomb,name,weight,zero)
      !
      integer, intent (in) :: nomb
      character (len=3), intent (in) :: name
      real(ki), intent (out), dimension(:) :: weight,zero
      !
      ! poids et abscisses pour les polynomes de Legendre 
      ! (calcule avec poid_gauss.m)
      !
      real(ki), dimension(7) :: zero_leg_7=(/&
                                  &-0.9491079123427585_ki,-0.7415311855993944_ki,&
                                  &-0.4058451513773972_ki,0.0000000000000000_ki,&
                                  &0.4058451513773972_ki,0.7415311855993944_ki,&
                                  &0.9491079123427585_ki&
                                  &/)
      real(ki), dimension(7) :: weight_leg_7=(/&
                                  &0.1294849661688697_ki,0.2797053914892767_ki,&
                                  &0.3818300505051189_ki,0.4179591836734694_ki,&
                                  &0.3818300505051189_ki,0.2797053914892767_ki,&
                                  &0.1294849661688697_ki&
                                  &/)
      real(ki), dimension(10) :: zero_leg_10=(/&
                                  &-0.9739065285171717_ki,-0.8650633666889845_ki,&
                                  &-0.6794095682990244_ki,-0.4333953941292472_ki,&
                                  &-0.1488743389816312_ki,0.1488743389816312_ki,&
                                  &0.4333953941292472_ki,0.6794095682990244_ki,&
                                  &0.8650633666889845_ki,0.9739065285171717_ki&
                                  &/)
      real(ki), dimension(10) :: weight_leg_10=(/&
                                  &0.0666713443086881_ki,0.1494513491505806_ki,&
                                  &0.2190863625159820_ki,0.2692667193099964_ki,&
                                  &0.2955242247147529_ki,0.2955242247147529_ki,&
                                  &0.2692667193099964_ki,0.2190863625159820_ki,&
                                  &0.1494513491505806_ki,0.0666713443086881_ki&
                                  &/)
      !
      ! poids et abscisses pour les polynomes de Kronrod (pris sur le WEB)
      !
      real(ki), dimension(15) :: zero_kro_15=(/&
                                  &-0.9914553711208126_ki,-0.9491079123427585_ki,&
                                  &-0.8648644233597691_ki,-0.7415311855993944_ki,&
                                  &-0.5860872354676911_ki,-0.4058451513773972_ki,&
                                  &-0.2077849550789850_ki,0.0e+00_ki,&
                                  &0.2077849550789850_ki,0.4058451513773972_ki,&
                                  &0.5860872354676911_ki,0.7415311855993944_ki,&
                                  &0.8648644233597691_ki,0.9491079123427585_ki,&
                                  &0.9914553711208126_ki&
                                  &/)
      real(ki), dimension(15) :: weight_kro_15=(/&
                                  &0.2293532201052922e-01_ki,0.6309209262997855e-01_ki,&
                                  &0.1047900103222502_ki,0.1406532597155259_ki,&
                                  &0.1690047266392679_ki,0.1903505780647854_ki,&
                                  &0.2044329400752989_ki,0.2094821410847278_ki,&
                                  &0.2044329400752989_ki,0.1903505780647854_ki,&
                                  &0.1690047266392679_ki,0.1406532597155259_ki,&
                                  &0.1047900103222502_ki,0.6309209262997855e-01_ki,&
                                  &0.2293532201052922e-01_ki&
                                  &/)
      !
      real(ki), dimension(21) :: zero_kro_21=(/&
                                    &-0.9956571630258081e+00_ki,-0.9739065285171717e+00_ki,&
                                    &-0.9301574913557082e+00_ki,-0.8650633666889845e+00_ki,&
                                    &-0.7808177265864169e+00_ki,-0.6794095682990244e+00_ki,&
                                    &-0.5627571346686047e+00_ki,-0.4333953941292472e+00_ki,&
                                    &-0.2943928627014602e+00_ki,-0.1488743389816312e+00_ki,&
                                    &0.0e+00_ki,0.1488743389816312e+00_ki,&
                                    &0.2943928627014602e+00_ki,0.4333953941292472e+00_ki,&
                                    &0.5627571346686047e+00_ki,0.6794095682990244e+00_ki,&
                                    &0.7808177265864169e+00_ki,0.8650633666889845e+00_ki,&
                                    &0.9301574913557082e+00_ki,0.9739065285171717e+00_ki,&
                                    &0.9956571630258081e+00_ki&
                                    &/)
      !
      real(ki), dimension(21) :: weight_kro_21=(/&
                                    &0.1169463886737187e-01_ki,0.3255816230796473e-01_ki,&
                                    &0.5475589657435200e-01_ki,0.7503967481091995e-01_ki,&
                                    &0.9312545458369761e-01_ki,0.1093871588022976e+00_ki,&
                                    &0.1234919762620659e+00_ki,0.1347092173114733e+00_ki,&
                                    &0.1427759385770601e+00_ki,0.1477391049013385e+00_ki,&
                                    &0.1494455540029169e+00_ki,0.1477391049013385e+00_ki,&
                                    &0.1427759385770601e+00_ki,0.1347092173114733e+00_ki,&
                                    &0.1234919762620659e+00_ki,0.1093871588022976e+00_ki,&
                                    &0.9312545458369761e-01_ki,0.7503967481091995e-01_ki,&
                                    &0.5475589657435200e-01_ki,0.3255816230796473e-01_ki,&
                                    &0.1169463886737187e-01_ki&
                                    &/)
      !
      if (name == 'leg') then
        !
        if (nomb == 7) then
          !
          weight = weight_leg_7
          zero = zero_leg_7
          !
        else if (nomb == 10) then
          !
          weight = weight_leg_10
          zero = zero_leg_10
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In subroutine give_me_the_weight'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the number of point is incorrect for & 
                          &Legendre weights: %d0'
          tab_erreur_par(2)%arg_int = nomb
          call catch_exception(0)
          !
        end if
        !
      else if (name == 'kro') then
        !
        if (nomb == 15) then
          !
          weight = weight_kro_15
          zero = zero_kro_15
          !
        else if (nomb == 21) then
          !
          weight = weight_kro_21
          zero = zero_kro_21
          !
        else
          !
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'In subroutine give_me_the_weight'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'the number of point is incorrect for & 
                          &Kronrod weights: %d0'
          tab_erreur_par(2)%arg_int = nomb
          call catch_exception(0)
          !
        end if
        !
      end if
      !
    end subroutine give_me_the_weight
    !
    !
    !****if* src/numerical/adapt_gauss/gauss1
    ! NAME
    !
    !  subroutine gauss1 (file src/numerical/adapt_gauss.f90)
    !
    ! USAGE
    !
    !  call gauss1(func,b_inf,b_sup,rest,err)
    !
    ! DESCRIPTION
    !
    !   This routine computes the one dimensional integration using 
    !   the Gauss-Kronrod method of the function func in the range
    !   b_inf to b_sup
    !
    !
    ! INPUTS
    !
    !  * func -- an external function from R to C
    !  * b_inf -- a real (type ki), lower bound of the integration
    !  * b_sup -- a real (type ki), upper bound of the integration
    !
    ! SIDE EFFECTS
    !
    !  no side effects
    !
    ! RETURN VALUE
    !
    !  * rest -- a complex (type ki), the result of the integration
    !  * err -- a complex (type ki), the estimate of the absolute error returned
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine gauss1(func,b_inf,b_sup,rest,err)
      !
      real(ki), intent (in) :: b_inf,b_sup
      complex(ki), intent (out) :: rest
      complex(ki), intent (out) :: err
      !
      interface
        !
        function func(x)
          use precision_golem
          real(ki), intent (in) :: x
          complex(ki) :: func
        end function func
        !
      end interface
      !
      real(ki), dimension(n1_par) :: weight,zero
      real(ki), dimension(:), allocatable :: weight_leg,zero_leg
      integer :: i,j
      real(ki) :: hm,hp
      real(ki) :: argg,argk
      complex(ki) :: tg,tk
      complex(ki) :: restg,restk,restg_leg
      !
      ! call for the Kronrod zeros and weights
      !
      call give_me_the_weight(n1_par,'kro',weight,zero)
      rest = 0._ki
      err = 1._ki
      allocate(weight_leg(n0_par),zero_leg(n0_par))
      !
      ! call for the Legendre zeros and weights
      !
      call give_me_the_weight(n0_par,'leg',weight_leg,zero_leg)
      restg_leg = rest
      restg = rest
      restk = rest
      hm = (b_sup - b_inf)/2._ki
      hp = (b_sup + b_inf)/2._ki
      !
      do i=1,n1_par-1,2
        !
        argk = hm*zero(i)+hp
        argg = hm*zero(i+1)+hp
        tk = func(argk)
        tg = func(argg)
        restk = restk + weight(i)*tk
        restg = restg + weight(i+1)*tg
        j = (i+1)/2
        restg_leg = restg_leg + weight_leg(j)*tg
        !
      end do
      !
      argk = hm*zero(n1_par)+hp
      restk = restk + weight(n1_par)*func(argk)
      rest = (restg+restk)*hm
      err = rest-restg_leg*hm
      !
      deallocate(weight_leg,zero_leg)
      !
    end subroutine gauss1
    !
end module adapt_gauss
