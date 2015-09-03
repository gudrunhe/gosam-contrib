!****h* src/integrals/two_point/generic_function_2p
!~ changed 13.5.2010 to include scale (mu^2)^eps
!~ the default scale is 1, defined in parametre.f90 
! NAME
!
!  Module generic_function_2p
!
! USAGE
!
!  use generic_function_2p
!
! DESCRIPTION
!
!  This module contains the generic routines to compute the
!  two point functions in n and n+2 dimensions
!
! OUTPUT
!
!  It exports two public routine:
!  * f2p -- a function to compute the two point function in n dimensions
!  * f2p_np2 -- a function to compute the two point function in n+2 dimensions
!
! USES
!
!  * precision (src/module/precision.f90)
!  * array (src/module/array.f90)
!  * logarithme (src/module/z_log.f90)
!  * tri_croissant (src/module/tri.f90)
!  * constante (src/module/constante.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * parametre (src/module/parametre.f90)
!  * function_2p0m_1mi (src/integrals/two_point/function_2p0m_1mi.f90)
!  * function_2p_m1m2 (src/integrals/two_point/function_2p_m1m2.f90)
!  * s_matrix_type (src/module/s_matrix_type.f90)
!
!*****
module generic_function_2p
  !
  use precision_golem
  use array
  use logarithme
  use tri_croissant
  use constante
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use parametre
  use equal
  use s_matrix_type
  use function_2p0m_1mi
  use function_2p_m1m2
  !
  implicit none
  !
  private
  !
  interface f2p
    !
    module procedure f2p_p
    !
  end interface
  !
  interface f2p_np2
    !
    module procedure f2p_np2_p
    !
  end interface
  !
  public :: f2p, f2p_np2
  public :: f2p_ra, f2p_np2_ra
  !
contains
  !
  !****f* src/integrals/two_point/generic_function_2p/f2p
  ! NAME
  !
  !  Function f2p
  !
  ! USAGE
  !
  !  cmplx_dim2 = f2p(s_mat_p,b_pro,parf1,parf2)
  !
  ! DESCRIPTION
  !
  !  This function computes the generic two point function in n dimensions, 
  !  with or without Feynman parameters in the numerator
  !
  ! INPUTS
  !
  !  * s_mat -- a derived type s_matrix_poly, the S matrix
  !  * b_pro -- an integer which represents the set of the four unpinched
  !    propagators
  !  * parf1 -- an integer (optional), the label of the one Feynman parameter
  !  * parf2 -- an integer (optional), the label of the second Feynman parameter
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki) array of rank 1 and shape 2
  !
  !
  !*****
  function f2p_ra(s_mat_p,b_pro,parf1,parf2)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent(in) :: b_pro
    integer, intent(in), optional :: parf1, parf2
    real(ki), dimension(4) :: f2p_ra
    complex(ki), dimension(2) :: f2p_ca
    !
    f2p_ca = f2p_p(s_mat_p,b_pro,parf1=parf1,parf2=parf2)
    !
    f2p_ra(1) = real(f2p_ca(1),ki)
    f2p_ra(2) = aimag(f2p_ca(1))
    f2p_ra(3) = real(f2p_ca(2),ki)
    f2p_ra(4) = aimag(f2p_ca(2))
    !
  end function f2p_ra
  !
  function f2p_p(s_mat_p,b_pro,parf1,parf2)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent(in) :: b_pro
    integer, intent(in), optional :: parf1, parf2
    complex(ki), dimension(2) :: f2p_p
    !
    if (iand(s_mat_p%b_cmplx, b_pro) .eq. 0 ) then
      !
      f2p_p = f2p_r(s_mat_p%pt_real, b_pro, parf1=parf1, parf2=parf2)
      !
    else
      !
      f2p_p = f2p_c(s_mat_p%pt_cmplx, b_pro, parf1=parf1, parf2=parf2)
      !
    end if
    !
  end function f2p_p
  !
  function f2p_r(s_mat_r,b_pro,parf1,parf2)
    !
    real(ki), intent (in), dimension(:,:) :: s_mat_r
    integer, intent (in) :: b_pro
    integer, intent (in), optional :: parf1,parf2
    complex(ki), dimension(2) :: f2p_r
    !
    real(ki), dimension(4) :: i2sonem,f2p_rr
    integer :: par1,par2
    integer, dimension(2) :: z_param_ini,z_param_out
    real(ki) :: arg1, s12, mass1, mass2, diffm
    integer :: m1,m2,dim_pro
    integer, dimension(2) :: s
    logical :: sz, mz1, mz2
    !
    par1 = 0
    par2 = 0
    !
    if (present(parf1)) par1 = parf1 
    if (present(parf2)) par2 = parf2
    !
    z_param_ini = (/ par1,par2 /)
    !
    where (z_param_ini /= 0)
      !
      z_param_ini = locateb(z_param_ini,b_pro)
      !
    elsewhere
      !
      z_param_ini = 0
      !
    end where
    !
    if ( minval(z_param_ini) == -1 ) then
      !
      f2p_rr(:) = 0._ki
      !
    else
      !
      call tri_int2(z_param_ini,z_param_out)
      !
      if (b_pro<256) then
        dim_pro = bit_count(b_pro)
        s = bit_sets(8*b_pro:8*b_pro+dim_pro-1)
      else
        dim_pro = countb(b_pro)
        s = unpackb(b_pro,dim_pro)
      end if
      !
      m1 = s(1)
      m2 = s(2)
      !
      arg1 = s_mat_r(m1,m2)
      !
      ! internal masses	
      mass1 = -s_mat_r(m1,m1)/2._ki
      mass2 = -s_mat_r(m2,m2)/2._ki
      s12 = arg1+mass1+mass2
      !
      call cut_s(s12,mass1,mass2)
      !
      mz1 = equal_real(mass1, zero,1000._ki)   ! 1000 added by MR 10.11.11
      mz2 = equal_real(mass2, zero,1000._ki)   ! 1000 added by MR 10.11.11
      sz = equal_real(s12,zero,1000._ki)   ! 1000 added by MR 10.11.11
      !
      diffm=mass1-mass2
      !	
      if ( (sz) .and. (mz1) .and. (mz2) ) then
        !
        f2p_rr(:) = 0._ki
        !  (scaleless two-point function is zero)
        !
      else if ( .not.(sz) .and. (mz1) .and. (mz2) ) then 
	      ! massless case
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_rr(1) = 1._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = 2._ki-real(z_log(-s12/mu2_scale_par,-1._ki))
            f2p_rr(4) = -aimag(z_log(-s12/mu2_scale_par,-1._ki))
            !
          else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 2._ki
            f2p_rr(4) = 0._ki
            !
          end if
          !
        else if ( (z_param_out(1) == 0) .and. (z_param_out(2) /= 0) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = 1._ki-real(z_log(-s12/mu2_scale_par,-1._ki))/2._ki
            f2p_rr(4) = -aimag(z_log(-s12/mu2_scale_par,-1._ki))/2._ki
            !
          else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 1._ki
            f2p_rr(4) = 0._ki
            !
          end if
          !
        else if ( (z_param_out(1) /= 0) .and. (z_param_out(2) /= 0) ) then
          !
          if (z_param_out(1) == z_param_out(2)) then
            !
            f2p_rr(1) = 1._ki/3._ki
            f2p_rr(2) = 0._ki
            !
            if (rat_or_tot_par%tot_selected) then
              !
              f2p_rr(3) = 13._ki/18._ki-real(z_log(-s12/mu2_scale_par,-1._ki))/3._ki
              f2p_rr(4) = -aimag(z_log(-s12/mu2_scale_par,-1._ki))/3._ki
              !
            else if (rat_or_tot_par%rat_selected) then
              !
              f2p_rr(3) = 13._ki/18._ki
              f2p_rr(4) = 0._ki
              !
            end if
            !
          else
            !
            f2p_rr(1) = 1._ki/6._ki
            f2p_rr(2) = 0._ki
            !
            if (rat_or_tot_par%tot_selected) then
              !
              f2p_rr(3) = 5._ki/18._ki-real(z_log(-s12/mu2_scale_par,-1._ki),ki)/6._ki
              f2p_rr(4) = -aimag(z_log(-s12/mu2_scale_par,-1._ki))/6._ki
              !
            else if (rat_or_tot_par%rat_selected) then
              !
              f2p_rr(3) = 5._ki/18._ki
              f2p_rr(4) = 0._ki
              !
            end if  ! end if rat or tot
            !
          end if ! end if z1==z2
          !
        end if  ! end test value of z1,z2
        !
        !
        !*************** massive cases *******************************
        ! added 07.08.09
        ! assumes real masses in numerator
        ! ************************************************************
      else if (  (sz) .and. (.not.(mz1)) .and. (mz2) ) then 
        ! case p^2=0, m1 nonzero, m2=0
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_rr = i20m1(mass1)
          !   
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = -(-1._ki + 2._ki*real(z_log(mass1/mu2_scale_par,-1._ki)))/4._ki
            f2p_rr(4) = -aimag(z_log(mass1/mu2_scale_par,-1._ki))/2._ki
            !
          else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 1._ki/4._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = -(-3._ki + 2._ki*real(z_log(mass1/mu2_scale_par,-1._ki)))/4._ki
            f2p_rr(4) = -aimag(z_log(mass1/mu2_scale_par,-1._ki))/2._ki
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 3._ki/4._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr(1) = 1._ki/3._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (1._ki - 3._ki*real(z_log(mass1/mu2_scale_par,-1._ki)))/9._ki
            f2p_rr(4) =  - aimag(z_log(mass1/mu2_scale_par,-1._ki))/3._ki
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 1._ki/9._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/6._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (5._ki - 6._ki*real(z_log(mass1/mu2_scale_par,-1._ki)))/36._ki
            f2p_rr(4) =  - aimag(z_log(mass1/mu2_scale_par,-1._ki))/6._ki
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 5._ki/36._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/3._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (11 - 6*real(z_log(mass1/mu2_scale_par,-1._ki)))/18._ki
            f2p_rr(4) =  - aimag(z_log(mass1/mu2_scale_par,-1._ki))/3._ki
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 11._ki/18._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        end if  ! end test value of z1,z2
        !
        ! ******************
      else if (  (sz) .and. (mz1) .and. (.not.(mz2)) ) then 
        ! case p^2=0, m2 nonzero, m1=0
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_rr = i20m1(mass2)
          !   
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = -(-3._ki + 2._ki*real(z_log(mass2/mu2_scale_par,-1._ki)))/4._ki
            f2p_rr(4) = -aimag(z_log(mass2/mu2_scale_par,-1._ki))/2._ki
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 3._ki/4._ki
            f2p_rr(4) = 0._ki
            !
	        end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = -(-1._ki + 2._ki*real(z_log(mass2/mu2_scale_par,-1._ki)))/4._ki
            f2p_rr(4) = -aimag(z_log(mass2/mu2_scale_par,-1._ki))/2._ki
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 1._ki/4._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr(1) = 1._ki/3._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (11._ki - 6._ki*real(z_log(mass2/mu2_scale_par,-1._ki)))/18._ki
            f2p_rr(4) =  - aimag(z_log(mass2/mu2_scale_par,-1._ki))/3._ki
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 11._ki/18._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/6._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (5._ki - 6._ki*real(z_log(mass2/mu2_scale_par,-1._ki)))/36._ki
            f2p_rr(4) =  - aimag(z_log(mass2/mu2_scale_par,-1._ki))/6._ki
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 5._ki/36._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/3._ki
          f2p_rr(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (1._ki - 3._ki*real(z_log(mass2/mu2_scale_par,-1._ki)))/9._ki
            f2p_rr(4) =  - aimag(z_log(mass2/mu2_scale_par,-1._ki))/3._ki
            !
          else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = 1._ki/9._ki
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        end if  ! end test value of z1,z2
        !
        ! ******************
        ! ** eq. (A.10) ****
      else if ( (sz) .and. (.not.(mz1)) .and. (equal_real(diffm,zero)) ) then 
        ! case p^2=0, m1 nonzero, m2=m1
        ! write(6,*) 'case (2c): s12 =0, m1 nonzero, m2=m1'
        !	  
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_rr = f2p0m_1mi(mass1,0,0)
          !   
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr = f2p0m_1mi(mass1,0,1)
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p0m_1mi(mass1,0,2)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr = f2p0m_1mi(mass1,1,1)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p0m_1mi(mass1,1,2)
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p0m_1mi(mass1,2,2)
          !
        end if  ! end test value of z1,z2
        ! 
        ! ******************
        ! ** eq. (A.8) ****
      else if ( (sz) .and. (.not.(mz1)) .and. (.not.(mz2)) .and. .not.(equal_real(diffm,zero)) ) then 
        ! case p^2=0, m1 nonzero, m2 nonzero, m2 NOT=m1
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_rr = f2p0m_m1m2(mass1,mass2,0,0)
          !   
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr = f2p0m_m1m2(mass1,mass2,0,1)
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p0m_m1m2(mass1,mass2,0,2)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr = f2p0m_m1m2(mass1,mass2,1,1)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p0m_m1m2(mass1,mass2,1,2)
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p0m_m1m2(mass1,mass2,2,2)
          !
        end if  ! end test value of z1,z2
        !
        ! ************ now case s12 nonzero **********************    
      else if (  (.not.(sz)) .and. (mz1) .and. (.not.(mz2)) ) then 
        ! case  p^2 nonzero, m1=0, m2 nonzero
        !
	      i2sonem=i2sm1(s12,mass2)
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_rr = i2sonem
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = -(mass2 - (mass2 + s12)*i2sonem(3) - & 
                 &                  mass2*real(z_log(mass2/mu2_scale_par,-1._ki)))/(2._ki*s12)
            f2p_rr(4) = -(           - (mass2 + s12)*i2sonem(4) - & 
                 &	                mass2*aimag(z_log(mass2/mu2_scale_par,-1._ki)))/(2._ki*s12)
            !
          else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = -( mass2 - (mass2 + s12)*i2sonem(3) )/(2._ki*s12)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = -( (mass2 - s12)*i2sonem(3) + mass2*(-1._ki + &
                 &        real(z_log(mass2/mu2_scale_par,-1._ki))) )/(2._ki*s12)
            f2p_rr(4) = -( (mass2 - s12)*i2sonem(4) + &
                 &        mass2*aimag(z_log(mass2/mu2_scale_par,-1._ki)) )/(2._ki*s12)
            !
          else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = -((mass2 - s12)*i2sonem(3) - mass2)/(2._ki*s12)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !	       
          f2p_rr(1) = 1._ki/3._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (-6._ki*mass2**2 - 9._ki*mass2*s12 + s12**2 + & 
                  &       6._ki*(mass2**2 + mass2*s12 + s12**2)*i2sonem(3) + & 
                  &       6._ki*mass2*(mass2 + s12)* &
                  &       real(z_log(mass2/mu2_scale_par,-1._ki)) )/(18._ki*s12**2)
            f2p_rr(4) = ( 6._ki*(mass2**2 + mass2*s12 + s12**2)*i2sonem(4) + & 
                  &       6._ki*mass2*(mass2 + s12)* &
                  &       aimag(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = ( -6._ki*mass2**2 - 9._ki*mass2*s12 + s12**2 + &
                 &       6._ki*(mass2**2 + mass2*s12 + s12**2)*i2sonem(3) )/(18._ki*s12**2)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot	       
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/6._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (6._ki*mass2**2 - s12**2 +   &
                 &      3._ki*(-2._ki*mass2**2 + mass2*s12 + s12**2)*i2sonem(3) + &
                 &      3._ki*mass2*(-2._ki*mass2 + s12)* &
                 &      real(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            !
            f2p_rr(4) = ( 3._ki*(-2._ki*mass2**2 + mass2*s12 + s12**2)*i2sonem(4) + &
                 &       3._ki*mass2*(-2._ki*mass2 + s12)* &
                 &       aimag(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = (6._ki*mass2**2 - s12**2 +   &
                 &      3._ki*(-2*mass2**2 + mass2*s12 + s12**2)*i2sonem(3) )/(18._ki*s12**2)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot	        
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/3._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (-6._ki*mass2**2 + 9._ki*mass2*s12 + s12**2 + &
                 &       6._ki*(mass2 - s12)**2*i2sonem(3) +          &
                 &       6._ki*mass2*(mass2 - 2._ki*s12)* &
                 &       real(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            f2p_rr(4) = ( 6._ki*(mass2 - s12)**2*i2sonem(4) + &
                 &       6._ki*mass2*(mass2 - 2._ki*s12)* &
                 &       aimag(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = (-6._ki*mass2**2 + 9._ki*mass2*s12 + s12**2 + &
                 &       6._ki*(mass2 - s12)**2*i2sonem(3) )/(18._ki*s12**2)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot        
          !
        end if  ! end test value of z1,z2
        ! ******************
      else if ( (.not.(sz)) .and. (.not.(mz1)) .and. (mz2) ) then 
        ! case p^2 nonzero, m1 nonzero, m2=0
        !
	      i2sonem=i2sm1(s12,mass1)
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_rr = i2sonem
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = -( (mass1 - s12)*i2sonem(3) + &
                &        mass1*(-1._ki + real(z_log(mass1/mu2_scale_par,-1._ki))))/(2._ki*s12)
            f2p_rr(4) = -( (mass1 - s12)*i2sonem(4) + &
                &        mass1*aimag(z_log(mass1/mu2_scale_par,-1._ki)))/(2._ki*s12)
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = -((mass1 - s12)*i2sonem(3) - mass1)/(2._ki*s12)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot     
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/2._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = -(mass1 - (mass1 + s12)*i2sonem(3) - &
                 &	      mass1*real(z_log(mass1/mu2_scale_par,-1._ki)))/(2._ki*s12)
            f2p_rr(4) = -( - (mass1 + s12)*i2sonem(4) - &
                 &        mass1*aimag(z_log(mass1/mu2_scale_par,-1._ki)))/(2._ki*s12)
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = -( mass1 - (mass1 + s12)*i2sonem(3) )/(2._ki*s12)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr(1) = 1._ki/3._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (-6._ki*mass1**2 + 9._ki*mass1*s12 + s12**2 + & 
                 &	     6._ki*(mass1 - s12)**2*i2sonem(3) +  &
                 &       6._ki*mass1*(mass1 - 2*s12)* &
                 &       real(z_log(mass1/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            f2p_rr(4) =( 6._ki*(mass1 - s12)**2*i2sonem(4) + &
                 &      6._ki*mass1*(mass1 - 2._ki*s12)* &
                 &      aimag(z_log(mass1/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = (-6._ki*mass1**2 + 9._ki*mass1*s12 + s12**2 + & 
                 &       6._ki*(mass1 - s12)**2*i2sonem(3) )/(18._ki*s12**2)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot        
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/6._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (6._ki*mass1**2 - s12**2 +   &
                 &      3._ki*(-2._ki*mass1**2 + mass1*s12 + s12**2)*i2sonem(3) + &
                 &      3._ki*mass1*(-2._ki*mass1 + s12)* &
                 &      real(z_log(mass1/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            f2p_rr(4) = ( 3._ki*(-2._ki*mass1**2 + mass1*s12 + s12**2)*i2sonem(4) + &
                 &       3._ki*mass1*(-2._ki*mass1 + s12)* &
                 &       aimag(z_log(mass1/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = (6._ki*mass1**2 - s12**2 +   &
                 &      3._ki*(-2._ki*mass1**2 + mass1*s12 + s12**2)*i2sonem(3) )/(18._ki*s12**2)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot	        
          !
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr(1) = 1._ki/3._ki
          f2p_rr(2) = 0._ki
          !     
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_rr(3) = (-6._ki*mass1**2 - 9._ki*mass1*s12 + s12**2 +  &
                 &	     6._ki*(mass1**2 + mass1*s12 + s12**2)*i2sonem(3) + &
                 &       6._ki*mass1*(mass1 + s12)* &
                 &       real(z_log(mass1/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            f2p_rr(4) = ( 6._ki*(mass1**2 + mass1*s12 + s12**2)*i2sonem(4) + &
                 &       6._ki*mass1*(mass1 + s12)* &
                 &       aimag(z_log(mass1/mu2_scale_par,-1._ki)))/(18._ki*s12**2)
            !
	        else if (rat_or_tot_par%rat_selected) then
            !
            f2p_rr(3) = (-6._ki*mass1**2 - 9._ki*mass1*s12 + s12**2 + &
                 &       6._ki*(mass1**2 + mass1*s12 + s12**2)*i2sonem(3))/(18._ki*s12**2)
            f2p_rr(4) = 0._ki
            !
          end if  ! end if rat or tot	       
          !
        end if  ! end test value of z1,z2
	      !
        ! ******************
      else if ( (.not.(sz)) .and. (.not.(mz1)) .and. (.not.(mz2)) ) then 
        ! case p^2 nonzero, m1 nonzero, m2 nonzero, eq.(A.5)
        ! includes case m1=m2
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_rr = f2p_m1m2(s12,mass1,mass2,0,0)
          !   
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr = f2p_m1m2(s12,mass1,mass2,0,1)
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p_m1m2(s12,mass1,mass2,0,2)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_rr = f2p_m1m2(s12,mass1,mass2,1,1)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p_m1m2(s12,mass1,mass2,1,2)
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_rr = f2p_m1m2(s12,mass1,mass2,2,2)
          !
        end if  ! end test value of z1,z2
        ! ******************************************************************
	      !
      else 
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'something wrong with arguments of two-point function f2p'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 's12= %f0'
        tab_erreur_par(2)%arg_real = s12
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'm1s= %f0'
        tab_erreur_par(3)%arg_real = mass1
        tab_erreur_par(4)%a_imprimer = .true.
        tab_erreur_par(4)%chaine = 'm2s= %f0'
        tab_erreur_par(4)%arg_real = mass2
        !
        call catch_exception(0)
      end if  ! end if s12,m1,m2 eq. zero
      !
    end if  ! end if ( minval(z_param_ini) == -1 )
    !
    f2p_r(1) = f2p_rr(1) + i_ * f2p_rr(2)
    f2p_r(2) = f2p_rr(3) + i_ * f2p_rr(4)
    !
  end function f2p_r
  !
  function f2p_c(s_mat_c,b_pro,parf1,parf2)
    !
    complex(ki), intent (in), dimension(:,:) :: s_mat_c
    integer, intent (in) :: b_pro
    integer, intent (in), optional :: parf1,parf2
    complex(ki), dimension(2) :: f2p_c
    !
    real(ki), dimension(4) :: f2p_cc
    complex(ki), dimension(2) :: i2sonem
    integer :: par1, par2
    integer, dimension(2) :: z_param_ini, z_param_out
    complex(ki) :: mass1, mass2, ratpart, diffm
    real(ki) :: s12
    integer :: m1,m2,dim_pro
    integer, dimension(2) :: s
    logical :: sz, mz1, mz2, diffz
    !
    par1 = 0
    par2 = 0
    !
    if (present(parf1)) par1 = parf1 
    if (present(parf2)) par2 = parf2
    !
    z_param_ini = (/ par1,par2 /)
    !
    where (z_param_ini /= 0)
      !
      z_param_ini = locateb(z_param_ini,b_pro)
      !
    elsewhere
      !
      z_param_ini = 0
      !
    end where
    !
    if ( minval(z_param_ini) == -1 ) then
      !
      f2p_c(:) = czero
      !
    else
      !
      call tri_int2(z_param_ini,z_param_out)
      !
      if (b_pro<256) then
        dim_pro = bit_count(b_pro)
         s = bit_sets(8*b_pro:8*b_pro+dim_pro-1)
      else
         dim_pro = countb(b_pro)
         s = unpackb(b_pro,dim_pro)
      end if
      !
      m1 = s(1)
      m2 = s(2)
      !
      ! internal masses	
      mass1 = -s_mat_c(m1,m1)/2._ki
      mass2 = -s_mat_c(m2,m2)/2._ki
      s12 = real(s_mat_c(m1,m2)+mass1+mass2,ki)
      !
      call cut_s(s12,mass1,mass2)
      !
      diffm = mass1-mass2
      !	
      mz1 = ( equal_real(real(mass1,ki), zero,1000._ki) .and. equal_real(aimag(mass1), zero,1000._ki) )  ! 1000 added by MR 10.11.11
      mz2 = ( equal_real(real(mass2,ki), zero,1000._ki) .and. equal_real(aimag(mass2), zero,1000._ki) )  ! 1000 added by MR 10.11.11
      sz = equal_real(s12,zero)  ! 1000 added by MR 10.11.11
      !
      diffz = ( equal_real(real(diffm,ki), zero) .and. equal_real(aimag(diffm), zero) ) ! -- not altered by MR
      !
      ! implement the massless case, for any complex masses of order 10^-14 that come this far.
      if ( (sz) .and. (mz1) .and. (mz2) ) then
        !
        f2p_cc(:) = 0._ki
        !  (scaleless two-point function is zero)
        !
      elseif ( .not.(sz) .and. (mz1) .and. (mz2) ) then 
	      ! massless case
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_cc(1) = 1._ki
          f2p_cc(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_cc(3) = 2._ki-real(z_log(-s12/mu2_scale_par,-1._ki))
            f2p_cc(4) = -aimag(z_log(-s12/mu2_scale_par,-1._ki))
            !
          else if (rat_or_tot_par%rat_selected) then
            !
            f2p_cc(3) = 2._ki
            f2p_cc(4) = 0._ki
            !
          end if
          !
        else if ( (z_param_out(1) == 0) .and. (z_param_out(2) /= 0) ) then
          !
          f2p_cc(1) = 1._ki/2._ki
          f2p_cc(2) = 0._ki
          !
          if (rat_or_tot_par%tot_selected) then
            !
            f2p_cc(3) = 1._ki-real(z_log(-s12/mu2_scale_par,-1._ki))/2._ki
            f2p_cc(4) = -aimag(z_log(-s12/mu2_scale_par,-1._ki))/2._ki
            !
          else if (rat_or_tot_par%rat_selected) then
            !
            f2p_cc(3) = 1._ki
            f2p_cc(4) = 0._ki
            !
          end if
          !
        else if ( (z_param_out(1) /= 0) .and. (z_param_out(2) /= 0) ) then
          !
          if (z_param_out(1) == z_param_out(2)) then
            !
            f2p_cc(1) = 1._ki/3._ki
            f2p_cc(2) = 0._ki
            !
            if (rat_or_tot_par%tot_selected) then
              !
              f2p_cc(3) = 13._ki/18._ki-real(z_log(-s12/mu2_scale_par,-1._ki))/3._ki
              f2p_cc(4) = -aimag(z_log(-s12/mu2_scale_par,-1._ki))/3._ki
              !
            else if (rat_or_tot_par%rat_selected) then
              !
              f2p_cc(3) = 13._ki/18._ki
              f2p_cc(4) = 0._ki
              !
            end if
            !
          else
            !
            f2p_cc(1) = 1._ki/6._ki
            f2p_cc(2) = 0._ki
            !
            if (rat_or_tot_par%tot_selected) then
              !
              f2p_cc(3) = 5._ki/18._ki-real(z_log(-s12/mu2_scale_par,-1._ki),ki)/6._ki
              f2p_cc(4) = -aimag(z_log(-s12/mu2_scale_par,-1._ki))/6._ki
              !
            else if (rat_or_tot_par%rat_selected) then
              !
              f2p_cc(3) = 5._ki/18._ki
              f2p_cc(4) = 0._ki
              !
            end if  ! end if rat or tot
            !
          end if ! end if z1==z2
          !
        end if  ! end test value of z1,z2
        ! 
        f2p_c(1) = f2p_cc(1) + i_ * f2p_cc(2)
        f2p_c(2) = f2p_cc(3) + i_ * f2p_cc(4)
        !
        !
        ! *************** massive cases, complex *******************************
        ! **  (this function is only called with at least one non-zero mass)  **
        ! **********************************************************************
        !
      else if (  sz .and. (.not. mz1) .and. mz2 ) then
        ! case p^2=0, m1 nonzero, m2=0
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_c = i20m1(mass1)
          !   
        else if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c(1) = 1._ki/2._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 1._ki/4._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) =  -(-1._ki + 2._ki*z_log(mass1/mu2_scale_par,-1._ki))/4._ki
            !
          end if  ! end if rat or tot
          !
        else if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/2._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 3._ki/4._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = -(-3._ki + 2._ki*z_log(mass1/mu2_scale_par,-1._ki))/4._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c(1) = 1._ki/3._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 1._ki/9._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = (1._ki - 3._ki*z_log(mass1/mu2_scale_par,-1._ki))/9._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/6._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 5._ki/36._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) =  (5._ki - 6._ki*z_log(mass1/mu2_scale_par,-1._ki))/36._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/3._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 11._ki/18._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = (11._ki - 6._ki*z_log(mass1/mu2_scale_par,-1._ki))/18._ki
            !
          end if  ! end if rat or tot
          !
        end if  ! end test value of z1,z2
        !
        ! ******************
      else if (  sz .and. mz1 .and. (.not. mz2) ) then
        ! case p^2=0, m2 nonzero, m1=0
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_c = i20m1(mass2)
          !   
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then 
          !
          f2p_c(1) = 1._ki/2._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 3._ki/4._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = -(-3._ki + 2._ki*z_log(mass2/mu2_scale_par,-1._ki))/4._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/2._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 1._ki/4._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = -(-1._ki + 2._ki*z_log(mass2/mu2_scale_par,-1._ki))/4._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c(1) = 1._ki/3._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 11._ki/18._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) =  (11._ki - 6._ki*z_log(mass2/mu2_scale_par,-1._ki))/18._ki
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/6._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 5._ki/36._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = (5._ki - 6._ki*z_log(mass2/mu2_scale_par,-1._ki))/36._ki
            !
          end if  ! end if rat or tot
          !      
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/3._ki
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = 1._ki/9._ki
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) =  (1 - 3*z_log(mass2/mu2_scale_par,-1._ki))/9._ki
            !
          end if  ! end if rat or tot
          !
        end if  ! end test value of z1,z2
        !
        ! ******************
        ! ** eq. (A.10) ****
        !
      else if (  sz .and. (.not. mz1) .and. diffz ) then 
        ! case p^2=0, m1 nonzero, m2=m1
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_c = f2p0m_1mi(mass1,0,0)
          !   
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c = f2p0m_1mi(mass1,0,1)
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p0m_1mi(mass1,0,2)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c = f2p0m_1mi(mass1,1,1)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p0m_1mi(mass1,1,2)
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p0m_1mi(mass1,2,2)
          !
        end if  ! end test value of z1,z2
        !    
        ! ******************
        ! ** eq. (A.8) ****
      else if (  sz .and. (.not. mz1) .and. (.not. mz2)  .and. (.not. diffz) ) then
        ! case p^2=0, m1 nonzero, m2 nonzero, m2 NOT=m1
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_c = f2p0m_m1m2(mass1,mass2,0,0)
          !   
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c = f2p0m_m1m2(mass1,mass2,0,1)
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p0m_m1m2(mass1,mass2,0,2)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c = f2p0m_m1m2(mass1,mass2,1,1)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p0m_m1m2(mass1,mass2,1,2)
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p0m_m1m2(mass1,mass2,2,2)
          !
        end if  ! end test value of z1,z2
        !
        ! ************ now case s12 nonzero **********************    
        !
      else if (  (.not. sz) .and. mz1 .and. (.not.mz2) ) then
        ! case  p^2 nonzero, m1=0, m2 nonzero
        !
        i2sonem=i2sm1(s12,mass2)
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_c = i2sonem
          !
        else if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c(1) = 1._ki/2._ki
          !     
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) =  -( mass2 - (mass2 + s12)*i2sonem(2) )/(2._ki*s12)
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = -(mass2 - (mass2 + s12)*i2sonem(2) - & 
                 &           mass2*z_log(mass2/mu2_scale_par,-1._ki))/(2._ki*s12)
            !
          end if  ! end if rat or tot
          !
        else if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/2._ki
          !     
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) =  -((mass2 - s12)*i2sonem(2) - mass2)/(2._ki*s12)
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = -((mass2 - s12)*i2sonem(2) + & 
               &           mass2*(-1._ki + z_log(mass2/mu2_scale_par,-1._ki)))/(2._ki*s12)
            !
          end if  ! end if rat or tot
          !
        else if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !	
          f2p_c(1) = 1._ki/3._ki
          !
          ratpart =  (-6._ki*mass2**2 - 9._ki*mass2*s12 + s12**2 + &
              &       6._ki*(mass2**2 + mass2*s12 + s12**2)*i2sonem(2))/(18._ki*s12**2)
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = ratpart
            !
          else if (rat_or_tot_par%tot_selected) then !!!checked!!!
            !
            f2p_c(2) = ratpart + mass2*(mass2 + s12)* &
              &         z_log(mass2/mu2_scale_par,-1._ki)/(3._ki*s12**2)
            !
          end if  ! end if rat or tot	       
          !
        else if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/6._ki
          !    
          ratpart =  (6._ki*mass2**2 - s12**2 +   &
               &      3._ki*(-2._ki*mass2**2 + mass2*s12 + s12**2)*i2sonem(2) )/(18._ki*s12**2)
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = ratpart
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = ratpart + mass2*(-2._ki*mass2 + s12)* &
              &         z_log(mass2/mu2_scale_par,-1._ki)/(6._ki*s12**2)
            !
          end if  ! end if rat or tot	               
          !   
        else if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/3._ki
          !     
          ratpart =  (-6._ki*mass2**2 + 9._ki*mass2*s12 + s12**2 + &
               &       6._ki*(mass2 - s12)**2*i2sonem(2) )/(18._ki*s12**2)
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = ratpart
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = ratpart + mass2*(mass2 - 2._ki*s12)* & 
             &         z_log(mass2/mu2_scale_par,-1._ki)/(3._ki*s12**2)
            !
          end if  ! end if rat or tot        
          !
        end if  ! end test value of z1,z2
        !
        ! ******************
        !
      else if ( (.not. sz) .and. (.not. mz1) .and. mz2 ) then
        ! case p^2 nonzero, m1 nonzero, m2=0
        !
        i2sonem = i2sm1(s12,mass1)
        !
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_c = i2sonem
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c(1) = 1._ki/2._ki
          !     
          ratpart =  -((mass1 - s12)*i2sonem(2) - mass1)/(2._ki*s12)
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = ratpart
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = ratpart - mass1*z_log(mass1/mu2_scale_par,-1._ki)/(2._ki*s12)
            !
          end if  ! end if rat or tot     
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/2._ki
          !     
          ratpart =  -( mass1 - (mass1 + s12)*i2sonem(2) )/(2._ki*s12)
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = ratpart
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = ratpart + mass1*z_log(mass1/mu2_scale_par,-1._ki)/(2._ki*s12)
            !
          end if  ! end if rat or tot
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c(1) = 1._ki/3._ki
          !     
          ratpart = (-6._ki*mass1**2 + 9._ki*mass1*s12 + s12**2 + & 
               &      6._ki*(mass1 - s12)**2*i2sonem(2) )/(18._ki*s12**2)
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = ratpart
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = ratpart + mass1*(mass1 - 2._ki*s12)* &
              &         z_log(mass1/mu2_scale_par,-1._ki)/(3._ki*s12**2)
            !
          end if  ! end if rat or tot        
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/6._ki
          !   
          ratpart =  (6._ki*mass1**2 - s12**2 +   &
               &      3._ki*(-2._ki*mass1**2 + mass1*s12 + s12**2)*i2sonem(2) )/(18._ki*s12**2)
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c = ratpart
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = ratpart + mass1*(-2._ki*mass1 + s12)* &
             &         z_log(mass1/mu2_scale_par,-1._ki)/(6._ki*s12**2)
            !
          end if  ! end if rat or tot	        
          !        
          !   
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c(1) = 1._ki/3._ki
          !     
          ratpart = (-6._ki*mass1**2 - 9._ki*mass1*s12 + s12**2 + &
               &      6._ki*(mass1**2 + mass1*s12 + s12**2)*i2sonem(2))/(18._ki*s12**2)
          !
          if (rat_or_tot_par%rat_selected) then
            !
            f2p_c(2) = ratpart
            !
          else if (rat_or_tot_par%tot_selected) then
            !
            f2p_c(2) = ratpart + mass1*(mass1 + s12)* &
              &         z_log(mass1/mu2_scale_par,-1._ki)/(3._ki*s12**2)
            !
          end if  ! end if rat or tot	       
          !
        end if  ! end test value of z1,z2
        !
        ! ******************
        !
      else if (  (.not. sz) .and. (.not. mz1) .and. (.not. mz2) ) then
        ! case p^2 nonzero, m1 nonzero, m2 nonzero, eq.(A.5)
        ! includes case m1=m2
        !   
        if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 0) ) then
          !
          f2p_c = f2p_m1m2(s12,mass1,mass2,0,0)
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c = f2p_m1m2(s12,mass1,mass2,0,1)
          !
        else  if ( (z_param_out(1) == 0) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p_m1m2(s12,mass1,mass2,0,2)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 1) ) then
          !
          f2p_c = f2p_m1m2(s12,mass1,mass2,1,1)
          !
        else  if ( (z_param_out(1) == 1) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p_m1m2(s12,mass1,mass2,1,2)
          !
        else  if ( (z_param_out(1) == 2) .and. (z_param_out(2) == 2) ) then
          !
          f2p_c = f2p_m1m2(s12,mass1,mass2,2,2)
          !
        end if  ! end test value of z1,z2
        ! ******************************************************************
        !
      else 
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'something wrong with arguments of two-point function f2p'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 's12= %f0'
        tab_erreur_par(2)%arg_real = s12
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'm1s= %f0'
        tab_erreur_par(3)%arg_real = mass1
        tab_erreur_par(4)%a_imprimer = .true.
        tab_erreur_par(4)%chaine = 'm2s= %f0'
        tab_erreur_par(4)%arg_real = mass2
        !
        call catch_exception(0)
        !
      end if  ! end if s12,m1,m2 eq. zero
      !
    end if  ! end if ( minval(z_param_ini) == -1 )
    !
  end function f2p_c
  !
  !****f* src/integrals/two_point/generic_function_2p/f2p_np2
  ! NAME
  !
  !  Function f2p_np2
  !
  ! USAGE
  !
  !  cmplx_dim2 = f2p_np2(s_mat_p,b_pro)
  !
  ! DESCRIPTION
  !
  !  This function computes the generic two point function in n+2 dimensions, 
  !  without Feynman parameters in the numerator
  !
  ! INPUTS
  !
  !  * s_mat_p -- a s_matrix_poly type object
  !  * b_pro -- an integer which represents the set of the four unpinched
  !    propagators
  !
  ! SIDE EFFECTS
  !
  !  No side effect
  !
  ! RETURN VALUE
  !
  !  It returns a complex (type ki) array of rank 1 and shape 2
  !
  ! EXAMPLE
  !
  !
  !
  !*****
  !
  function f2p_np2_ra(s_mat_p,b_pro)
    !
    type(s_matrix_poly) :: s_mat_p
    integer, intent(in) :: b_pro
    real(ki), dimension(4) :: f2p_np2_ra
    complex(ki), dimension(2) :: f2p_np2_ca
    !
    f2p_np2_ca = f2p_np2_p(s_mat_p,b_pro)
    !
    f2p_np2_ra(1) = real(f2p_np2_ca(1),ki)
    f2p_np2_ra(2) = aimag(f2p_np2_ca(1))
    f2p_np2_ra(3) = real(f2p_np2_ca(2),ki)
    f2p_np2_ra(4) = aimag(f2p_np2_ca(2))
    !
  end function f2p_np2_ra
  !
  function f2p_np2_p(s_mat_p,b_pro)
    type(s_matrix_poly) :: s_mat_p
    integer, intent(in) :: b_pro
    complex(ki), dimension(2) :: f2p_np2_p
    !
    if (iand(s_mat_p%b_cmplx, b_pro) .eq. 0 ) then
       !
       f2p_np2_p = f2p_np2_r(s_mat_p%pt_real, b_pro)
       !
    else
       !
       f2p_np2_p = f2p_np2_c(s_mat_p%pt_cmplx, b_pro)
       !
    end if
    !
  end function f2p_np2_p
  !
  function f2p_np2_r(s_mat_r,b_pro)
    !
    real(ki), intent (in), dimension(:,:) :: s_mat_r
    integer, intent (in) :: b_pro
    complex(ki),dimension(2) :: f2p_np2_r
    ! 
    real(ki), dimension(4) :: f2p_np2_rr,i2sonem,i2sca
    real(ki) :: arg1,s12,mass1,mass2,diffrm,small
    real(ki) :: arg_log
    integer :: m1,m2
    integer, dimension(2) :: s
    logical :: sz, mz1, mz2
    !
    small=1.e-6_ki
    !
    s = unpackb(b_pro,countb(b_pro))
    !
    m1 = s(1)
    m2 = s(2)
    !
    arg1 = s_mat_r(m1,m2)
    !
    ! internal masses	
    mass1 = -s_mat_r(m1,m1)/2._ki
    mass2 = -s_mat_r(m2,m2)/2._ki
    diffrm = sqrt(mass1)-sqrt(mass2)
    !
    s12 = arg1+mass1+mass2
    !
    call cut_s(s12, mass1, mass2)
    !
     mz1 = equal_real(mass1, zero,1000._ki)   ! 1000 added by MR 10.11.11
     mz2 = equal_real(mass2, zero,1000._ki)   ! 1000 added by MR 10.11.11
     sz = equal_real(s12, zero,1000._ki)   ! 1000 added by MR 10.11.11
    !
    !
    if ( (sz) .and. (mz1) .and. (mz2) ) then
       !
       f2p_np2_rr(:) = 0._ki
       ! (f2p_np2_rr with no scale is zero)
       !
    else if ( (.not.(sz)) .and. (mz1) .and. (mz2) ) then 
       ! massless case
       !
       arg_log = arg1/mu2_scale_par
       !
       f2p_np2_rr(1) = 1._ki/6._ki
       f2p_np2_rr(2) = 0._ki
       !
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_rr(3) = 4._ki/9._ki-real(z_log(-arg_log,-1._ki),ki)/6._ki
          f2p_np2_rr(4) = -aimag(z_log(-arg_log,-1._ki))/6._ki
          !
       else if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_rr(3) = 4._ki/9._ki
          f2p_np2_rr(4) = 0._ki
          !
       end if
       !
       f2p_np2_rr = s12*f2p_np2_rr
       !
       !*************** massive cases *******************************
       ! added 07.08.09
       ! ************************************************************
    else if ( (sz) .and. (.not.(mz1)) .and. (mz2) ) then 
       ! case p^2=0, m1 nonzero, m2=0
       !
       f2p_np2_rr(1) = -1._ki/2._ki
       f2p_np2_rr(2) = 0._ki
       !     
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_rr(3) = -2._ki*(3._ki - 2*real(z_log(mass1/mu2_scale_par,-1._ki)))/8._ki
          f2p_np2_rr(4) = aimag(z_log(mass1/mu2_scale_par,-1._ki))/2._ki
          !
       else if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_rr(3) = -3._ki/4._ki
          f2p_np2_rr(4) = 0._ki
          !
       end if  ! end if rat or tot
       !   
       f2p_np2_rr = mass1*f2p_np2_rr
       !
       ! ******************
    else if ( (sz) .and. (mz1) .and. (.not.(mz2)) ) then 
       ! case p^2=0, m2 nonzero, m1=0
       !
       f2p_np2_rr(1) = -1._ki/2._ki
       f2p_np2_rr(2) = 0._ki
       !     
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_rr(3) = - (3._ki - 2*real(z_log(mass2/mu2_scale_par,-1._ki)))/4._ki
          f2p_np2_rr(4) = aimag(z_log(mass2/mu2_scale_par,-1._ki))/2._ki
          !
       else if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_rr(3) = -3._ki/4._ki
          f2p_np2_rr(4) = 0._ki
          !
       end if  ! end if rat or tot
       !   
       f2p_np2_rr = mass2*f2p_np2_rr
       !
       ! ******************
       ! ** eq. (A.10), uses f2p_np2_rr = -2*B22 ****
    else if ( (sz).and.(.not.(mz1)).and. &
         &    (equal_real(diffrm,zero)) ) then 
       ! case p^2=0, m1 nonzero, m2=m1
       !
       f2p_np2_rr = -mass1*i20m1(mass1)
       !   
       ! ******************
    else if ( (sz) .and. (.not.(mz1)) .and. (.not.(mz2)) .and. .not.(equal_real(diffrm,zero)) ) then 
       ! case p^2=0, m1 nonzero, m2 nonzero, m1 not= m2
       !
       f2p_np2_rr(1) = -(mass1+mass2)/2._ki
       f2p_np2_rr(2) = 0._ki
       !     
       if (abs(diffrm) > small ) then 
          ! 
          if (rat_or_tot_par%tot_selected) then
             !
             f2p_np2_rr(3) = - (3._ki*mass1**2 - 3._ki*mass2**2 - & 
                  &            2._ki*mass1**2*real(z_log(mass1/mu2_scale_par,-1._ki)) + &
                  &            2._ki*mass2**2*real(z_log(mass2/mu2_scale_par,-1._ki)))/ &
                  &            (4._ki*(mass1 - mass2))
             f2p_np2_rr(4) = - ( - 2*mass1**2*aimag(z_log(mass1/mu2_scale_par,-1._ki)) + &
                  &            2*mass2**2*aimag(z_log(mass2/mu2_scale_par,-1._ki)))/ &
                  &            (4._ki*(mass1 - mass2))
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p_np2_rr(3) = - (3*mass1**2 - 3*mass2**2)/(4._ki*(mass1 - mass2))
             f2p_np2_rr(4) = 0._ki
             !
          end if  ! end if rat or tot
          !
       else  ! use expansion in (m2sq-m1sq) up to order 3
          !    
          if (rat_or_tot_par%tot_selected) then
             !
             f2p_np2_rr(3) = (mass1 + mass2)*( -19._ki*mass1**2 + 8._ki*mass1*mass2 - mass2**2 + &
                  &         12._ki*mass1**2*real(z_log(mass1/mu2_scale_par,-1._ki)))/(24._ki*mass1**2)
             f2p_np2_rr(4) = (mass1 + mass2)*( 12._ki*mass1**2* &
                  &         aimag(z_log(mass1/mu2_scale_par,-1._ki)))/(24._ki*mass1**2)
             !
          else if (rat_or_tot_par%rat_selected) then
             !
             f2p_np2_rr(3) = (mass1 + mass2)*( -19._ki*mass1**2 + & 
                  &          8._ki*mass1*mass2 - mass2**2 )/(24._ki*mass1**2)
             f2p_np2_rr(4) = 0._ki
             !
          end if  ! end if rat or tot
          !    
       end if ! end if abs(diffrm) > small
       !	   
       !   
       ! ******************
       !  s12 nonzero
       !********************
    else if ( (.not.(sz)) .and. (mz1) .and. (.not.(mz2)) ) then 
       ! case p^2 nonzero, m1=0, m2 nonzero
       !
       i2sonem=i2sm1(s12,mass2)
       !
       f2p_np2_rr(1) = -(3._ki*mass2 - s12)/6._ki
       f2p_np2_rr(2) = 0._ki
       !
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_rr(3) =  - (3._ki*mass2**2 + 9._ki*mass2*s12 - 2._ki*s12**2 - &
               &	     3._ki*(mass2 - s12)**2*i2sonem(3) - 3._ki*mass2*(mass2 + s12)* &
               &             real(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12)
          !
          f2p_np2_rr(4) =  - ( - 3._ki*(mass2 - s12)**2*i2sonem(4) - 3._ki*mass2*(mass2 + s12)* &
               &              aimag(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12)
          !
       else if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_rr(3) = - (3._ki*mass2**2 + 9._ki*mass2*s12 - 2._ki*s12**2 - & 
               &            3._ki*(mass2 - s12)**2*i2sonem(3) )/(18._ki*s12)
          !
          f2p_np2_rr(4) =  0._ki
          !
       end if  ! end if rat or tot
       !   
       ! ******************
    else if ( (.not.(sz)) .and. (.not.(mz1)) .and. (mz2) ) then 
       ! case p^2 nonzero, m1 nonzero, m2=0
       !
       i2sonem=i2sm1(s12,mass1)
       !
       f2p_np2_rr(1) = -(3._ki*mass1 - s12)/6._ki
       f2p_np2_rr(2) = 0._ki
       !
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_rr(3) = - (3._ki*mass1**2 + 9._ki*mass1*s12 - 2._ki*s12**2 - &
               &            3._ki*(mass1 - s12)**2*i2sonem(3) -  &
               &            3._ki*mass1*(mass1 + s12)* &
               &            real(z_log(mass1/mu2_scale_par,-1._ki)))/(18._ki*s12)
          !
          f2p_np2_rr(4) = - ( - 3._ki*(mass1 - s12)**2*i2sonem(4) - &
               &             3._ki*mass1*(mass1 + s12)* &
               &             aimag(z_log(mass1/mu2_scale_par,-1._ki)))/(18._ki*s12)
          !
       else if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_rr(3) = - (3*mass1**2 + 9*mass1*s12 - 2*s12**2 - &
               &            3*(mass1 - s12)**2*i2sonem(3) )/(18._ki*s12)
          !
          f2p_np2_rr(4) =  0._ki
          !
       end if  ! end if rat or tot
       !   
       ! ******************
    else if ( (.not.(sz)) .and. (.not.(mz1)) .and. (.not.(mz2)) ) then 
       ! case p^2 nonzero, m1 nonzero, m2 nonzero, -2*B22 of eq.(A.5)
       !
       i2sca=i2sm1m2(s12,mass1,mass2)
       !
       f2p_np2_rr(1) = - (3._ki*mass1 + 3._ki*mass2 - s12)/6._ki
       f2p_np2_rr(2) = 0._ki
       !
       if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_rr(3) = - (3._ki*mass1**2 - 6._ki*mass1*mass2 + 3._ki*mass2**2 + & 
               &            9._ki*mass1*s12 + 9._ki*mass2*s12 - 2._ki*s12**2 - &
               &            3._ki*(mass1**2 + (mass2 - s12)**2 - &
               &            2._ki*mass1*(mass2 + s12))*i2sca(3) - &
               &            3._ki*mass1*(mass1 - mass2 + s12)* &
               &            real(z_log(mass1/mu2_scale_par,-1._ki)) + &
               &            3._ki*mass2*(mass1 - mass2 - s12)* &
               &            real(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12)
          !
          f2p_np2_rr(4) = - (-3._ki*( mass1**2 + (mass2 - s12)**2  - &
               &             2._ki*mass1*(mass2 + s12) )*i2sca(4) - &
               &             3._ki*mass1*(mass1 - mass2 + s12)* &
               &             aimag(z_log(mass1/mu2_scale_par,-1._ki)) + &
               &             3._ki*mass2*(mass1 - mass2 - s12)* &
               &             aimag(z_log(mass2/mu2_scale_par,-1._ki)))/(18._ki*s12)
          !
       else if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_rr(3) = - (3._ki*mass1**2 - 6._ki*mass1*mass2 + 3._ki*mass2**2 + &
               &            9._ki*mass1*s12 + 9._ki*mass2*s12 - 2._ki*s12**2 - &
               &            3._ki*(mass1**2 + (mass2 - s12)**2 - &
               &            2._ki*mass1*(mass2 + s12))*i2sca(3) )/(18._ki*s12)
          !
          f2p_np2_rr(4) =  0._ki
          !
       end if  ! end if rat or tot
       !   
       ! ******************************************************************
    end if  ! end test if s12,m1,m2 zero
    !
    f2p_np2_r(1) = f2p_np2_rr(1) + i_ * f2p_np2_rr(2)
    f2p_np2_r(2) = f2p_np2_rr(3) + i_ * f2p_np2_rr(4)    
    !
  end function f2p_np2_r
  !
  !
  function f2p_np2_c(s_mat_c,b_pro)
    !
    complex(ki), intent (in), dimension(:,:) :: s_mat_c
    integer, intent (in) :: b_pro
    complex(ki),dimension(2) :: f2p_np2_c
    !
    complex(ki), dimension(2) :: i2sonem_c, i2sca_c
    complex(ki) :: mass1, mass2, lambda
    real(ki) :: s12, diffrm, small
    integer :: m1,m2
    integer, dimension(2) :: s
    logical :: sz, m1z, m2z, diffz
    !
    small=1.e-6_ki
    !
    s = unpackb(b_pro,countb(b_pro))
    !
    m1 = s(1)
    m2 = s(2)
    !
    ! internal masses	
    mass1 = -s_mat_c(m1,m1)/2._ki
    mass2 = -s_mat_c(m2,m2)/2._ki
    diffrm = sqrt(abs(mass1-mass2))
    !
    s12 = real(s_mat_c(m1,m2)+mass1+mass2,ki)
    !
    call cut_s(s12, mass1, mass2)
    !
     m1z = ( equal_real(real(mass1,ki), zero,1000._ki) .and. equal_real(aimag(mass1), zero,1000._ki) )   ! 1000 added by MR 10.11.11
     m2z = ( equal_real(real(mass2,ki), zero,1000._ki) .and. equal_real(aimag(mass2), zero,1000._ki) )   ! 1000 added by MR 10.11.11
    !
     sz = equal_real(s12,zero,1000._ki)   ! 1000 added by MR 10.11.11
    !
    diffz = equal_real(diffrm, zero)  ! -- no 1000 put in
    !
    ! *************** massive cases, complex ************************
    ! ** (this function is called with at least one non-zero mass) **
    ! ***************************************************************
    !
    if ( sz .and. (.not. m1z) .and. m2z ) then
       ! case p^2=0, m1 nonzero, m2=0
       !
       f2p_np2_c(1) = -1._ki/2._ki
       !     
       if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_c(2) = -3._ki/4._ki
          !
       else if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_c(2) = - (3._ki - 2._ki*z_log(mass1/mu2_scale_par,-1._ki))/4._ki
          !
       end if  ! end if rat or tot
       !   
       f2p_np2_c = mass1*f2p_np2_c
       !
       ! ******************
    else if ( sz .and. m1z .and. (.not. m2z) ) then
       ! case p^2=0, m2 nonzero, m1=0
       !
       f2p_np2_c(1) = -1._ki/2._ki
       !     
       if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_c(2) = -3._ki/4._ki
          !
       else if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_c(2) = - (3._ki - 2._ki*z_log(mass2/mu2_scale_par,-1._ki))/4._ki
          !
       end if  ! end if rat or tot
       !   
       f2p_np2_c = mass2*f2p_np2_c
       !
       ! ******************
       ! ** eq. (A.10), uses f2p_np2_c = -2*B22 ****
    else if ( sz .and. (.not. m1z) .and. diffz ) then
       ! case p^2=0, m1 nonzero, m2=m1
       !
       f2p_np2_c = -mass1*i20m1(mass1)
       !   
       ! ******************
    else if ( sz .and. (.not. m1z) .and. (.not. m2z) .and. (.not. diffz) ) then
       ! case p^2=0, m1 nonzero, m2 nonzero, m1 not= m2
       !
       f2p_np2_c(1) = -(mass1+mass2)/2._ki
       !     
       if (diffrm > small ) then
          ! 
          if (rat_or_tot_par%rat_selected) then
             !
             f2p_np2_c(2) = -3._ki*(mass1 + mass2)/4._ki
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p_np2_c(2) = (-3*mass1 - 3*mass2  + & 
                  &             (2*mass1**2*z_log(mass1/mu2_scale_par,-1._ki) - &
                  &             2*mass2**2*z_log(mass2/mu2_scale_par,-1._ki) )/(mass1 - mass2) )/4._ki
             !
          end if  ! end if rat or tot
          !
       else  ! use expansion in (m2sq-m1sq) up to order 3
          !    
          if (rat_or_tot_par%rat_selected) then
             !
             f2p_np2_c(2) =  (mass1 + mass2)*(  &
                  &	         -19._ki*mass1**2 + 8._ki*mass1*mass2 - mass2**2 )/(24._ki*mass1**2)
             !
          else if (rat_or_tot_par%tot_selected) then
             !
             f2p_np2_c(2) =  (mass1 + mass2)*(  &
                  &	        -19._ki*mass1**2 + 8._ki*mass1*mass2 - mass2**2 + &
                  &              12._ki*mass1**2*z_log(mass1/mu2_scale_par,-1._ki) )/(24._ki*mass1**2)
             !
          end if  ! end if rat or tot
          !    
       end if ! end if abs(diffrm) > small
       !	   
       !   
       ! ******************
       !  s12 nonzero
       !********************
       !
    else if ( (.not. sz) .and. m1z .and. (.not. m2z) ) then
       ! case p^2 nonzero, m1=0, m2 nonzero
       !
       i2sonem_c = i2sm1(s12,mass2)
       !
       f2p_np2_c(1) = (s12 - 3._ki*mass2)/6._ki
       !
       if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_c(2) = - (3._ki*mass2**2 + 9._ki*mass2*s12 - 2._ki*s12**2 - & 
               &                3._ki*(mass2 - s12)**2*i2sonem_c(2) )/(18._ki*s12)
          !
       else if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_c(2) =  - (3._ki*mass2**2 + 9._ki*mass2*s12 - 2._ki*s12**2 - &
               &	         3._ki*(mass2 - s12)**2*i2sonem_c(2) -  &
               &                 3._ki*mass2*(mass2 + s12)* &
               &                 z_log(mass2/mu2_scale_par,-1._ki))/(18._ki*s12)
          !
       end if  ! end if rat or tot
       !   
       ! ******************
    else if ( (.not. sz) .and. (.not. m1z) .and. m2z ) then 
       ! case p^2 nonzero, m1 nonzero, m2=0
       !
       i2sonem_c=i2sm1(s12,mass1)
       !
       f2p_np2_c(1) = (s12 - 3._ki*mass1)/6._ki
       !
       if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_c(2) = - (3._ki*mass1**2 + 9._ki*mass1*s12 - 2._ki*s12**2 - &
               &                3._ki*(mass1 - s12)**2*i2sonem_c(2) )/(18._ki*s12) 
          !
       else if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_c(2) = - (3._ki*mass1**2 + 9._ki*mass1*s12 - 2._ki*s12**2 - &
               &                3._ki*(mass1 - s12)**2*i2sonem_c(2) -  &
               &                3._ki*mass1*(mass1 + s12)* &
               &                z_log(mass1/mu2_scale_par,-1._ki))/(18._ki*s12)
          !
       end if  ! end if rat or tot
       !   
       ! ******************
    else if ( (.not. sz) .and. (.not. m1z) .and. (.not.m2z) ) then
       ! case p^2 nonzero, m1 nonzero, m2 nonzero, -2*B22 of eq.(A.5)
       !
       i2sca_c=i2sm1m2(s12,mass1,mass2)
       !
       f2p_np2_c(1) = (s12 - 3._ki*mass1 - 3._ki*mass2)/6._ki
       !
       lambda = s12**2 + mass1**2 + mass2**2 - 2._ki*s12*(mass1 + mass2) - 2._ki*mass1*mass2
       !
       if (rat_or_tot_par%rat_selected) then
          !
          f2p_np2_c(2) = - (3._ki*(mass1 - mass2)**2 + 9._ki*s12*(mass1 + mass2) - &
               &                2._ki*s12**2 - 3._ki*lambda*i2sca_c(2) )/(18._ki*s12)
          !
       else if (rat_or_tot_par%tot_selected) then
          !
          f2p_np2_c(2) = - (3._ki*(mass1 - mass2)**2 + 9._ki*s12*(mass1 + mass2) - &
               &                2._ki*s12**2 - 3._ki*lambda*i2sca_c(2) - &
               &                3._ki*mass1*(mass1 - mass2 + s12)* &
               &                z_log(mass1/mu2_scale_par,-1._ki) + &
               &                3._ki*mass2*(mass1 - mass2 - s12)* &
               &                z_log(mass2/mu2_scale_par,-1._ki) )/(18._ki*s12)
          !
       end if  ! end if rat or tot
       !   
       ! ******************************************************************
    end if  ! end test if s12,m1,m2 zero
    !
  end function f2p_np2_c
  !
end module generic_function_2p
