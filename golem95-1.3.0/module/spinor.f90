! 
!****h* src/module/spinor
! NAME
!
!  Module spinor
!
! USAGE
!
!  use spinor
!
! DESCRIPTION
!
!  This module contains all the function to compute the spinorial products,
!  scalar products and epsilon_tensor
!
! OUTPUT
!
!  It exports:
!  * ket -- a function to compute the ket spinor
!  * bra -- a function to compute the bra spinor
!  * pslash -- a function to compute p^{\mu} \gamma_{\mu}
!  * bra_ket -- a function to compute the spinorial product
!  * eps_prod_sca -- a function to compute the scalar product e_i.p_j
!  * eps_prod_eps -- a function to compute the scalar product e_i.e_j
!  * scalar -- a function to compute the scalar product
!  * e_ -- a function to compute the epsilon tensor
!
! USES
!
!  * precision_golem (src/module/precision_golem.f90)
!
!*****
module spinor
  !
  use precision_golem
  use constante, only : i_
  implicit none
  !
  private 
  !
  complex(ki), dimension(16), parameter :: g0_col = (/0._ki,0._ki,1._ki,0._ki,&
                                                    &0._ki,0._ki,0._ki,1._ki,&
                                                    &1._ki,0._ki,0._ki,0._ki,&
                                                    &0._ki,1._ki,0._ki,0._ki/)
  complex(ki), dimension(16), parameter :: g1_col = (/0._ki,0._ki,0._ki,1._ki,&
                                                    &0._ki,0._ki,1._ki,0._ki,&
                                                    &0._ki,-1._ki,0._ki,0._ki,&
                                                    &-1._ki,0._ki,0._ki,0._ki/)
  complex(ki), dimension(16), parameter :: g2_col = (/0._ki,0._ki,0._ki,1._ki,&
                                                    &0._ki,0._ki,-1._ki,0._ki,&
                                                    &0._ki,-1._ki,0._ki,0._ki,&
                                                    &1._ki,0._ki,0._ki,0._ki/)
  complex(ki), dimension(16), parameter :: g3_col = (/0._ki,0._ki,1._ki,0._ki,&
                                                    &0._ki,0._ki,0._ki,-1._ki,&
                                                    &-1._ki,0._ki,0._ki,0._ki,&
                                                    &0._ki,1._ki,0._ki,0._ki/)
  !complex(ki), dimension(4,4), parameter :: gamma0 = reshape(g0_col,(/4,4/))
  !complex(ki), dimension(4,4), parameter :: gamma1 = reshape(g1_col,(/4,4/))
  !complex(ki), dimension(4,4), parameter :: gamma2 = i_*reshape(g2_col,(/4,4/))
  !complex(ki), dimension(4,4), parameter :: gamma3 = reshape(g3_col,(/4,4/))
  !
  public ket,bra,pslash,bra_ket,eps_prod_sca,eps_prod_eps,scalar,e_
  !
  contains
    !
    !****f* src/module/spinor/ket
    ! NAME
    !
    !  Function ket
    !
    ! USAGE
    !
    !  complex_dim4_1 = ket(p,i)
    !
    ! DESCRIPTION
    !
    !  This function computes the spinor using the chinese's paper
    !  Nucl. Phys. B291 (1987) 392-428 equation A.16
    !  modified for non physical configuration E < 0.
    !  The functions bra and ket verify the conditions:
    !  <-p-|q+> = <p-|-q+> = i <p-|q+>
    !  <-p+|q-> = <p+|-q-> = i <p+|q->
    !  <-p-|-q+> = - <p-|q+>
    !  <-p+|-q-> = - <p+|q->
    !
    ! INPUTS
    !
    !  * p -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !  * i -- an integer, the value of the helicity (= 1,-1)
    !
    ! SIDE EFFECTS
    !
    !  No side effect (pure function)
    !
    ! RETURN VALUE
    !
    !  It returns a complex array (type ki) of rank 2 and shape 4,1
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    pure function ket(p,i)
      !
      real(ki), dimension(4), intent(in) :: p
      integer, intent(in) :: i
      complex(ki), dimension(4,1) :: ket
      !
      real(ki) :: p_plus,p_moins,t1
      complex(ki) :: p_perp,c1,extra_phase
      !
      p_plus = p(1)+p(4)
      p_moins = p(1)-p(4)
      p_perp = sign(1._ki,p(1))*(p(2)+i_*p(3))
      extra_phase = (1._ki,0._ki)
      !
      ! for non physical configuration
      !
      if (p(1) < 0._ki) then
        !
        extra_phase = i_
        !
      end if
      !
      t1 = sqrt(abs(p_plus))
      !
      if (p_plus == 0._ki) then
        !
        c1 = sqrt(p_moins)
        !
      else
        !
        c1 = p_perp/t1
        !
      end if
      !
      if (i == 1) then
        !
        ket(1,:) = extra_phase*t1
        ket(2,:) = extra_phase*c1
        ket(3,:) = 0._ki
        ket(4,:) = 0._ki
        !
      else if (i == -1) then
        !
        ket(1,:) = 0._ki
        ket(2,:) = 0._ki
        ket(3,:) = extra_phase*conjg(c1)
        ket(4,:) = -extra_phase*t1
        !
      end if
      !
    end function ket
    !
    !****f* src/module/spinor/bra
    ! NAME
    !
    !  Function bra
    !
    ! USAGE
    !
    !  complex_dim1_4 = bra(p,i)
    !
    ! DESCRIPTION
    !
    !  This function computes the spinor bra using 
    !  bra(p) = ket(p)^{\dagger} \gamma_0
    !  modified for non physical configuration p_0 < 0.
    !  The functions bra and ket verify the conditions:
    !  <-p-|q+> = <p-|-q+> = i <p-|q+>
    !  <-p+|q-> = <p+|-q-> = i <p+|q->
    !  <-p-|-q+> = - <p-|q+>
    !  <-p+|-q-> = - <p+|q->
    !
    ! INPUTS
    !
    !  * p -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !  * i -- an integer, the value of the helicity (= 1,-1)
    !
    ! SIDE EFFECTS
    !
    !  No side effect (pure function)
    !
    ! RETURN VALUE
    !
    !  It returns a complex array (type ki) of rank 2 and shape 1,4
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    pure function bra(p,i)
      !
      real(ki), dimension(4), intent(in) :: p
      integer, intent(in) :: i
      complex(ki), dimension(1,4) :: bra
      !
      complex(ki), dimension(1,4) :: tra
      complex(ki), dimension(4,4) :: gamma0
      
      gamma0 = reshape(g0_col,(/4,4/))
      !
      tra = transpose(ket(p,i))
      bra = sign(1._ki,p(1))*conjg(tra)
      bra = matmul(bra,gamma0)
      !
    end function bra
    !
    !****f* src/module/spinor/pslash
    ! NAME
    !
    !  Function pslash
    !
    ! USAGE
    !
    !  complex_dim4_4 = pslash(p)
    !
    ! DESCRIPTION
    !
    !  This function computes p_{\mu} \gamma^{\mu}, i.e.
    !  p0 gamma0 - p1 gamma1 - p2 gamma2 - p3 gamma3
    !  taking the Chinese convention for the gamma matrices in
    !  Weyl representation
    !
    ! INPUTS
    !
    !  * p -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !
    ! SIDE EFFECTS
    !
    !  No side effect (pure function)
    !
    ! RETURN VALUE
    !
    !  It returns a complex array (type ki) of rank 2 and shape 4,4
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    pure function pslash(p)
      !
      real(ki), dimension(4), intent(in) :: p
      complex(ki), dimension(4,4) :: pslash
      !
      real(ki) :: p_plus,p_moins
      complex(ki) :: p_perp
      !
      p_plus = p(1) + p(4)
      p_moins = p(1) - p(4)
      p_perp = p(2) + i_*p(3)
      !
      pslash = 0._ki
      pslash(1,3) = p_plus
      pslash(1,4) = conjg(p_perp)
      pslash(2,3) = p_perp
      pslash(2,4) = p_moins
      pslash(3,1) = p_moins
      pslash(3,2) = -conjg(p_perp)
      pslash(4,1) = -p_perp
      pslash(4,2) = p_plus
      !
    end function pslash
    !
    !****f* src/module/spinor/bra_ket
    ! NAME
    !
    !  Function bra_ket
    !
    ! USAGE
    !
    !  complex = bra_ket(p1,i1,p2,i2,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)
    !
    ! DESCRIPTION
    !
    !  This function computes <p1 i1|k1slash*k2slash*...*k10slash|p2 i2>
    !  i1 and i2 = +/- 1 are the helicities
    !  where the inner argument k1slash,k2slash,...,k10slash are optional
    !
    ! INPUTS
    !
    !  * p1 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !  * i1 -- an integer, the value of the helicity (= 1,-1)
    !  * p2 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !  * i2 -- an integer, the value of the helicity (= 1,-1)
    !  * k1 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k2 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k3 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k4 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k5 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k6 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k7 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k8 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k9 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !  * k10 -- a real array (type ki) of rank 1, shape 4; a 4-momentum optional
    !
    ! SIDE EFFECTS
    !
    !  No side effect (pure function)
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
    pure function bra_ket(p1,i1,p2,i2,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)
      !
      real(ki), dimension(4), intent(in) :: p1,p2
      real(ki), dimension(4), intent(in), optional :: k1,k2,k3,k4,k5,&
                                                     &k6,k7,k8,k9,k10
      integer, intent(in) :: i1,i2
      complex(ki) :: bra_ket
      !
      complex(ki), dimension(1,1) :: temp
      complex(ki), dimension(4,4) :: c_mat
      complex(ki), dimension(4,1) :: c_col
      integer :: nb_arg
      logical :: test
      !
      nb_arg = 2
      ! calcul du nombre d'arguments optionnels
      if (present(k1)) nb_arg = nb_arg + 1
      if (present(k2)) nb_arg = nb_arg + 1
      if (present(k3)) nb_arg = nb_arg + 1
      if (present(k4)) nb_arg = nb_arg + 1
      if (present(k5)) nb_arg = nb_arg + 1
      if (present(k6)) nb_arg = nb_arg + 1
      if (present(k7)) nb_arg = nb_arg + 1
      if (present(k8)) nb_arg = nb_arg + 1
      if (present(k9)) nb_arg = nb_arg + 1
      if (present(k10)) nb_arg = nb_arg + 1
      !
      test = ( (modulo(nb_arg,2) == 0) .and. (i1*i2 == -1) ) .or. &
             ( (modulo(nb_arg,2) == 1) .and. (i1*i2 == 1) )
      !
      if ( test ) then
        !
        if (present(k1)) then
          !
          c_mat = pslash(k1)
          !
          if (present(k2)) c_mat = matmul(c_mat,pslash(k2))
          if (present(k3)) c_mat = matmul(c_mat,pslash(k3))
          if (present(k4)) c_mat = matmul(c_mat,pslash(k4))
          if (present(k5)) c_mat = matmul(c_mat,pslash(k5))
          if (present(k6)) c_mat = matmul(c_mat,pslash(k6))
          if (present(k7)) c_mat = matmul(c_mat,pslash(k7))
          if (present(k8)) c_mat = matmul(c_mat,pslash(k8))
          if (present(k9)) c_mat = matmul(c_mat,pslash(k9))
          if (present(k10)) c_mat = matmul(c_mat,pslash(k10))
          !
          c_col = matmul(c_mat,ket(p2,i2))
          temp = matmul(bra(p1,i1),c_col)
          !
        else
          !
          temp = matmul(bra(p1,i1),ket(p2,i2))
          !
        end if
        !
        bra_ket = temp(1,1)
        !
      else
        !
        bra_ket = 0._ki
        !
      end if
      !
    end function bra_ket
    !
    !****f* src/module/spinor/eps_prod_sca
    ! NAME
    !
    !  Function eps_prod_sca
    !
    ! USAGE
    !
    !  complex = eps_prod_sca(i,r1,p1,p2)
    !
    ! DESCRIPTION
    !
    !  This function computes e^i(p1).p2 where r1 is the reference momentum
    !  be careful that p2 is assumed to be a lightlike vector
    !
    ! INPUTS
    !
    !  * i -- an integer, the value of the helicity (= 1,-1)
    !  * r1 -- a real array (type ki) of rank 1, shape 4; the refrence momentum
    !  * p1 -- a real array (type ki) of rank 1, shape 4; the momentum of the spin 1
    !  * p2 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !
    ! SIDE EFFECTS
    !
    !  No side effect (pure function)
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
    pure function eps_prod_sca(i,r1,p1,p2)
      !
      integer, intent(in) :: i
      real(ki), dimension(4), intent(in) :: r1,p1,p2
      complex(ki) :: eps_prod_sca
      !
      complex(ki) :: ctemp,cjtemp
      real(ki) :: denom
      !
      ctemp = bra_ket(r1,-1,p1,1)
      cjtemp = conjg(ctemp)
      denom = ctemp*cjtemp*sqrt(2._ki)
      !
      if (i == 1) then
        !
        eps_prod_sca = bra_ket(r1,-1,p1,-1,p2)*cjtemp/denom
        !
      else if (i == -1) then
        !
        eps_prod_sca = bra_ket(r1,1,p1,1,p2)*ctemp/denom
      else
        !
        eps_prod_sca = huge(1.0_ki)
        !
      end if
       !
   end function eps_prod_sca
    !
    !****f* src/module/spinor/eps_prod_eps
    ! NAME
    !
    !  Function eps_prod_eps
    !
    ! USAGE
    !
    !  complex = eps_prod_eps(i1,r1,p1,i2,r2,p2)
    !
    ! DESCRIPTION
    !
    !  This function computes e^i(p1).e^j(p) where r1 is the reference momemtum 
    !  for e(p1) and  r2 is the reference momemtum for e(p2)
    !
    ! INPUTS
    !
    !  * i1 -- an integer, the value of the helicity (= 1,-1)
    !  * r1 -- a real array (type ki) of rank 1, shape 4; the refrence momentum
    !  * p1 -- a real array (type ki) of rank 1, shape 4; the momentum of the spin 1
    !  * i2 -- an integer, the value of the helicity (= 1,-1)
    !  * r2 -- a real array (type ki) of rank 1, shape 4; the refrence momentum
    !  * p2 -- a real array (type ki) of rank 1, shape 4; the momentum of the spin 1
    !
    ! SIDE EFFECTS
    !
    !  No side effect (pure function)
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
    pure function eps_prod_eps(i1,r1,p1,i2,r2,p2)
      !
      integer, intent(in) :: i1,i2
      real(ki), dimension(4), intent(in) :: r1,p1,r2,p2
      complex(ki) :: eps_prod_eps
      !
      complex(ki) :: c1temp,c1jtemp,c2temp,c2jtemp
      real(ki) :: denom
      !
      c1temp = bra_ket(r1,-1,p1,1)
      c1jtemp = conjg(c1temp)
      c2temp = bra_ket(r2,-1,p2,1)
      c2jtemp = conjg(c2temp)
      denom = c1temp*c1jtemp*c2temp*c2jtemp
      !
      if ( (i1 == 1) .and. (i2 == 1) ) then
        !
        eps_prod_eps = bra_ket(r2,-1,r1,1)*bra_ket(p1,1,p2,-1) &
                        *c1jtemp*c2jtemp/denom
        !
      else if ( (i1 == 1) .and. (i2 == -1) ) then
        !
        eps_prod_eps = bra_ket(r2,1,p1,-1)*bra_ket(r1,-1,p2,1) &
                        *c1jtemp*c2temp/denom
        !
      else if ( (i1 == -1) .and. (i2 == 1) ) then
        !
        eps_prod_eps = bra_ket(r1,1,p2,-1)*bra_ket(r2,-1,p1,1) &
                        *c2jtemp*c1temp/denom
        !
      else if ( (i1 == -1) .and. (i2 == -1) ) then
        !
        eps_prod_eps = bra_ket(p2,-1,p1,1)*bra_ket(r1,1,r2,-1) &
                        *c1temp*c2temp/denom
      else
        eps_prod_eps = huge(1.0_ki)
        !
      end if
      !
    end function eps_prod_eps
    !
    !****f* src/module/spinor/scalar
    ! NAME
    !
    !  Function scalar
    !
    ! USAGE
    !
    !  real = scalar(p1,p2)
    !
    ! DESCRIPTION
    !
    !  This function compute the scalar product of two 4 momentum
    !
    ! INPUTS
    !
    !  * p1 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !  * p2 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !
    ! SIDE EFFECTS
    !
    !  No side effect (pure function)
    !
    ! RETURN VALUE
    !
    !  It returns a real (type ki)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    pure function scalar(p1,p2)
      !
      real(ki), intent (in), dimension(4) :: p1,p2
      real(ki) :: scalar
      !
      scalar = p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3) - p1(4)*p2(4)
      !
    end function scalar
    !
    !****f* src/module/spinor/e_
    ! NAME
    !
    !  Function e_
    !
    ! USAGE
    !
    !  complex = e_(k1,k2,k3,k4)
    !
    ! DESCRIPTION
    !
    !  This function gives the antisymetric tensor epsilon
    !  From Thomas Reiter
    !
    ! INPUTS
    !
    !  * k1 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !  * k2 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !  * k3 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !  * k4 -- a real array (type ki) of rank 1, shape 4; a 4-momentum
    !
    ! SIDE EFFECTS
    !
    !  No side effect (pure function)
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
    pure function  e_(k1,k2,k3,k4)
      !
      real(ki), intent (in), dimension(4) :: k1, k2, k3, k4
      complex(ki) :: e_
      !
      real(ki) :: res
      real(ki) :: k12, k23, k34, k13, k14, k24
      !
      k12 = k3(1)*k4(2)-k3(2)*k4(1)
      k23 = k3(2)*k4(3)-k3(3)*k4(2)
      k34 = k3(3)*k4(4)-k3(4)*k4(3)
      k13 = k3(1)*k4(3)-k3(3)*k4(1)
      k14 = k3(1)*k4(4)-k3(4)*k4(1)
      k24 = k3(2)*k4(4)-k3(4)*k4(2)
      !   
      res =  k1(1)*(k2(2)*k34 - k2(3)*k24 + k2(4)*k23)&
      &    + k1(2)*(k2(3)*k14 - k2(1)*k34 - k2(4)*k13)&
      &    + k1(3)*(k2(1)*k24 - k2(2)*k14 + k2(4)*k12)&
      &    + k1(4)*(k2(2)*k13 - k2(1)*k23 - k2(3)*k12)
      !
      e_ =  i_*res
      !
    end function e_
    !
end module spinor
!
