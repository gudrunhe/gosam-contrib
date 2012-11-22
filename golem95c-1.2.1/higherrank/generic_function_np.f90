! vim:ts=3:sw=3
!****h* src/higherrank/generic_function_np.f90
!
!
! NAME
!
!  Module generic_function_np
!
! USAGE
!
!  use generic_function_np
!
! DESCRIPTION
!
!  This module contains the generic routines to compute the
!  n<=4 point functions in m dimensions with arbitrary rank
!
! OUTPUT
!
!  It exports the public routine: fnp_generic
!
! USES
!
!  * TODO
!  * array (src/module/array.f90)
!  * cache (src/module/cache.f90)
!  * constante (src/module/constante.f90)
!  * equal (src/module/equal.f90)
!  * form_factor_type (src/module/form_factor_type.f90)
!  * generic_function_3p (src/integrals/generic_function_3p.f90)
!  * matrice_s (src/kinematic/matrice_s.f90)
!  * precision (src/module/precision_golem.f90)
!  * s_matrix_type (src/module/s_matrix_type.f90)
!  * sortie_erreur (src/module/sortie_erreur.f90)
!  * tri_croissant (src/module/tri.f90)
!
!*****
module generic_function_np
  !
  use precision_golem
  use form_factor_4p
  use form_factor_3p
  use form_factor_2p
  use form_factor_1p
  use generic_function_4p
  use generic_function_3p
  use generic_function_2p
  use generic_function_1p
  use matrice_s
  use array
  use inverse_matrice
  use form_factor_type
  use s_matrix_type
  use constante, only: czero
  use logarithme
  use sortie_erreur
  use parametre
  use equal
  use cache
  !

  implicit none
  !
  private
  !
public :: fnp_generic
private :: reduce_generic, reduce_pave_generic
!
!

private ::  f2p_ndim_0p_generic;

!
  !
  contains
    !
    !

    !  fnp_generic()
    !
    !
    ! calculates form factor integrals in a generic way by
    ! calling reduce_generic() or reduce_pave_generic()
    !
    ! Parameters:
    !        leg_count
    !        dim_nplus => dimension = 4+dim_nplus - 2*epsilon
    !        b_pin : pinned legs
    !        l_count : number of Feynman parameters in the numerator
    !        l : array with Feynman parameters
    !        depth (optional,internal only): current recursion depth

    recursive function fnp_generic(leg_count,dim_nplus,b_pin,l_count,l,depth) result(return_val)
       implicit none
       integer, intent(in) :: leg_count
       integer, intent(in) :: dim_nplus
       integer, intent (in) :: b_pin
       integer, intent (in) :: l_count
       integer, intent (in),dimension(:) :: l
       integer, intent (in), optional :: depth
       type(form_factor) :: return_val, control
       complex(ki) :: detS
       integer :: cur_depth,i,b_used,hash

       real(ki) ::  limit_small_detS, abs_sumb
       integer, save :: miss

       limit_small_detS = 1E-15_ki


       if (present(depth)) then
          cur_depth=depth
       else
          cur_depth=0
       end if

        b_used=pminus(b_ref,b_pin)




       ! simple Cache system
       hash=0
      if (leg_count<=4 .and. cur_depth <= 100) then ! do not use cache for control
          ! max 6 legs supported
          do i = 1, l_count
             hash = hash*7 + l(i)
          end do
          hash = hash*128 + b_used
          hash = hash*20 + dim_nplus
          hash = hash*7 + leg_count
          ! search backwards in cache, last entries are more probable
          do i = cache_generic_count, 1, -1
             if (cache_generic_tag(i)==hash) then
                  return_val=cache_generic(i)
                  return
             end if
             miss = miss + 1
          end do
       end if

       do i  = 1, l_count
           if (locateb(l(i),b_used) > leg_count .and. leg_count>1 ) then
              write (*, *) "WRONG CALL:", leg_count, dim_nplus,l_count, b_used
              write (*, *) "depth:",cur_depth,"", "legs:", leg_count,"", "dim-n:",dim_nplus,"", "l_count",l_count, "l",l
              write (*, *) l(i), b_pin, locateb(l(i),b_used)

              return_val = 0
              return

              ! stop
           end if
        end do

        ! test if S-matrix contains only zeros

        if(test_mat_zero(s_mat_p,leg_count,b_pin)) then
           return_val=0
           !if (cur_depth==0) then
           !   write (*,*) "Zero S-matrix."
           !end if
           return
        end if


      if ((leg_count>=2)  &
             .and. (cur_depth<=100) .and. ( dim_nplus>=2 .or. (l_count>leg_count) &
                                           .or. (leg_count==4 .and. l_count>0) ) ) then
          detS = calc_determinant(s_mat_p,leg_count,b_pin)

          abs_sumb = abs(sumb(b_pin))

          ! test if detS is below threshold or if sumb is NaN
          if (( abs(detS) <=limit_small_detS) .or. (.not. abs_sumb>=0._ki)) then

             return_val=reduce_pave_generic(leg_count,dim_nplus,b_pin,l_count,l,cur_depth)


             if (hash /= 0 .and. cache_generic_count < 1000) then
                cache_generic_count = cache_generic_count + 1
                cache_generic(cache_generic_count) = return_val
                cache_generic_tag(cache_generic_count) = hash
             end if

             return
          end if
       end if

       return_val =  reduce_generic(leg_count,dim_nplus,b_pin,l_count,l,cur_depth)

       if ( (hash /= 0) .and. (cache_generic_count < 1000)) then
          cache_generic_count = cache_generic_count + 1
          cache_generic(cache_generic_count) = return_val
          cache_generic_tag(cache_generic_count) = hash
       end if

    end function fnp_generic

    !
    !
    ! reduce_generic
    !
    !
    !
    ! uses already implemented integrals or reduce assuming detS /= 0
    !
    recursive function reduce_generic(leg_count,dim_nplus,b_pin,l_count,l,depth) result(return_val)
        implicit none
        integer, intent(in) :: leg_count
        integer, intent(in) :: dim_nplus
        integer, intent (in) :: b_pin
        integer, intent (in) :: l_count
        integer, intent (in),dimension(:) :: l
        integer, intent (in), optional :: depth
        integer :: i,j,k,b_used,b_tmp,ib,cur_depth
        type(form_factor) :: return_val, temp1, temp2
        complex(ki), dimension(3) :: ret_temp, tmp4,tmp5
        complex(ki)  :: mass1, tmp
        real(ki) :: r_tmp1,r_tmp2
        integer :: n ! dimension = n + dim_nplus =  4-2epsilon + dim_nplus TODO: float
        integer :: y
        integer,dimension(size(l)) :: l_tmp1
        integer,dimension(1) :: s1
        logical :: usable

        n = 4
        if (present(depth)) then
                cur_depth=depth
        else
                cur_depth=0
        end if



        b_used=pminus(b_ref,b_pin)
        ret_temp=0._ki




        ret_temp = 0

        ! ===
        ! some checks if everything is consistent
        ! ===
        if ((countb(b_used) /= leg_count ) .or. (size(l)/= l_count)) then
                tab_erreur_par(1)%a_imprimer = .true.
                tab_erreur_par(1)%chaine = 'Assert failed: internal error in file generic_function_np.f90'
                call catch_exception(0)
                stop
        end if

        if (leg_count<=0 .or. leg_count>=5) then
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'Internal error: not implemented, in file generic_function_np.f90.'
              tab_erreur_par(2)%a_imprimer = .true.
              tab_erreur_par(2)%chaine = '  Only leg_count = 1,2,3,4 supported.'
              call catch_exception(0)
        end if

        ! =================================================
        ! use existing (master) integrals if available
        !
        !
        if (l_count == 0 .and. dim_nplus==0) then
           select case (leg_count)
           !case (5)
           !   ret_temp(2:3) = a50(b_pin)
           case (4)
              return_val = a40(b_pin)
           case (3)
              ret_temp(1:3) = f3p(s_mat_p,b_used)
              return_val= ret_temp
           case (2)
              ret_temp(2:3) = f2p(s_mat_p,b_used)
              return_val= ret_temp
           case (1)
              ret_temp(2:3) = f1p(s_mat_p,b_used)
              return_val= ret_temp
           case default
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'Internal error: not yet implemented, in file generic_function_np.f90'
              call catch_exception(0)
           end select
           return
        end if
        if (l_count == 2 .and. dim_nplus==0 .and. leg_count==2) then
           ret_temp(:) = czero
           ret_temp(2:3) = f2p(s_mat_p, b_used,l(1),l(2))
           return_val= ret_temp
           return
        end if


        !if (.false. .and. l_count == 2 .and. dim_nplus==2 .and. leg_count==3) then
        !!if (l_count == 2 .and. dim_nplus==2 .and. leg_count==3 ) then
        !    ! TODO Trap:
        !      ! The program stops because error in function f3p2m_np2 no need of
        !      ! two mass six dimensional 3-point function with more
        !      ! than one Feynman parameter in the numerator
        !   ret_temp(:) = czero
        !   ret_temp(2:3) = f3p_np2(s_mat_p, b_used,l(1),l(2))
        !   return_val= ret_temp
        !   return
        !end if


        if (l_count == 3 .and. dim_nplus==2 .and. leg_count==4) then
           ret_temp(:) = czero
           ret_temp(3) = f4p_np2(s_mat_p, b_used,b_pin,l(1),l(2),l(3))
           return_val= ret_temp
           return
        end if
        !if (l_count == 4 .and. dim_nplus==2 .and. leg_count==4) then
        !   ret_temp(:) = czero
        !   ret_temp(3) = f4p_np2(s_mat_p, b_used,b_pin,l(1),l(2),l(3),l(4))
        !   return_val= ret_temp
        !   return
        !end if
        if (l_count == 4 .and. dim_nplus==0 .and. leg_count==4) then
           return_val= a44(l(1),l(2),l(3),l(4),b_pin)
           !ret_temp(1:3) = f4p(s_mat_p,b_used,b_pin,(1),l(2),l(3),l(4))
           return
        end if


        if (l_count == 3 .and. dim_nplus==0 .and. leg_count==4) then
           return_val= -a43(l(1),l(2),l(3),b_pin)
           return
        end if

        if (l_count == 2 .and. dim_nplus==0 .and. leg_count==4) then
           return_val= a42(l(1),l(2),b_pin)
           return
        end if

        if (leg_count==3 .and. dim_nplus==0) then
           select case (l_count)
           case(3)
               return_val = f3p(s_mat_p, b_used,l(1),l(2),l(3))
               return
           case(2)
               return_val = f3p(s_mat_p, b_used,l(1),l(2))
               return
           end select
       end if

       if (l_count == 1 .and. dim_nplus==0) then
           select case (leg_count)
           case(4)
              return_val=-a41(l(1),b_pin)
              return
           case (3)
              !     return_val = -a31(l(1),b_pin)
              return_val = f3p(s_mat_p, b_used,l(1))
              return
           case (2)
              ! return_val = -a21(l(1),b_pin)
              ret_temp(2:3) = f2p(s_mat_p,b_used,l(1))
              return_val=ret_temp
              return
          end select
       end if
       if (l_count == 1 .and. dim_nplus==2) then
                select case (leg_count)
                case (4)
                        ret_temp(:) = czero
                        ret_temp(3) = f4p_np2(s_mat_p,b_used,b_pin,l(1))
                        return_val=ret_temp
                        return
               ! case (3)
               !         ret_temp(:) = czero
               !         ret_temp(2:3) = f3p_np2(s_mat_p,b_used,l(1))
               !         return_val=ret_temp
               !         return
               end select
        end if


        if (l_count == 0 .and. dim_nplus==2 ) then
           select case (leg_count)
           case (4)
              ret_temp(:) = czero
              ret_temp(3) = f4p_np2(s_mat_p,b_used,b_pin)
              return_val=ret_temp
              return
           case (3)
              ret_temp(:) = czero
              ret_temp(2:3) = f3p_np2(s_mat_p,b_used)
              return_val=ret_temp
              return
           case (2)
              ret_temp(:) = czero
              ret_temp(2:3) = f2p_np2(s_mat_p,b_used)
              return_val=ret_temp
              return
           end select
        end if

        if (l_count == 2 .and. dim_nplus==2 ) then
             select case (leg_count)
                    case (4)
                       return_val= f4p_np2(s_mat_p,b_used,b_pin,l(1),l(2))
                       return
             end select
       end if

        if (l_count == 0 .and. dim_nplus==4 ) then
             select case (leg_count)
                    case (41)
                       return_val= f4p_np4(s_mat_p,b_used,b_pin)
                       return
             end select
       end if




       if (leg_count == 1 .and. dim_nplus==0) then
          return_val = a10(b_pin)
          return
       end if


        if (leg_count==1 .and. dim_nplus /= 0) then
                ret_temp(:) = czero
                s1=unpackb(b_used,countb(b_used))
                k=s1(1)

                if (iand(s_mat_p%b_cmplx, b_used) .eq. 0 ) then
                   mass1=-s_mat_p%pt_real(k,k)/2
                else
                   mass1=-s_mat_p%pt_cmplx(k,k)/2
                end if

                return_val = f1p_ndim_generic(mass1,dim_nplus)

                return
        end if


        ! ==================================================
        !     Reduction of higher dimensional scalar integrals
        !
        ! !     cf. eq. (68) and (83) from  hep-ph/0504267
        !
        if (l_count == 0) then ! dim_nplus>0
               !
               if (dim_nplus<=0) then
                  write (*, *) "internal error"
                  stop
               end if

               temp1 = fnp_generic(leg_count,dim_nplus-2,b_pin,l_count,l,cur_depth+1)
               ib=b_used
               j=0
               temp2=0
               do while (ib /= 0)
                       if (iand(ib,1) == 1 ) then
                               b_tmp=ibset(b_pin,j)
                               temp2=temp2 + b(j,b_pin)*fnp_generic(leg_count-1, dim_nplus-2,b_tmp,0,l,cur_depth+1)
                       end if
                       ib=ishft(ib,-1)
                       j=j+1
               end do

               tmp4(1) =temp1%a-temp2%a
               tmp4(2) =temp1%b-temp2%b
               tmp4(3) =temp1%c-temp2%c

               tmp=n+dim_nplus - leg_count - 1

               if (tmp/=0) then
                 ret_temp(3)= 1._ki/tmp * (tmp4(3) + 2._ki/tmp * tmp4(2)  &
                                                   + 4._ki/tmp/tmp * tmp4(1) )
                 ret_temp(2)= 1._ki/tmp * (tmp4(2) + 2._ki/tmp * tmp4(1))
                 ret_temp(1)= 1._ki/tmp * tmp4(1)
               else ! tmp==0
                        tab_erreur_par(1)%a_imprimer = .true.
                        tab_erreur_par(1)%chaine = 'Internal error:' &
                                //'This case is not yet implemented, in file generic_function_np.f90'
                        call catch_exception(0)
               end if
               if (sumb(b_pin)==0) then
                       return_val = 0
                         tab_erreur_par(1)%a_imprimer = .true.
                        tab_erreur_par(1)%chaine = 'Internal error:' &
                                //'This case (sumb=0) is not yet implemented, in file generic_function_np.f90'
                        call catch_exception(0)
                       return
               end if

               return_val = ret_temp/sumb(b_pin)
               return
        end if


        ! =======================================================
        !     Reduction of integrals with Feynman parameter in the
        !     numerator
        !
        !Reduce, l_count >= 1
        !  sec. 4.2 in hep-ph/0504267
        !  general formular: eq. (58) in hep-ph/9911342 (different convention)
        !



        ! == First terms
        ! ==  - invS * I^(n+2)

        temp1=0

        if (leg_count < 5) then

           l_tmp1=0
           do k = 2, l_count
                   l_tmp1(1:k-2) =l(2:k-1)
                   l_tmp1(k-1:l_count-2) =l(k+1:l_count)
                   i=max(l_count-2,0)
                   temp1=temp1-inv_s(l(1),l(k),b_pin)*fnp_generic(leg_count,dim_nplus+2,b_pin,i,l_tmp1(1:l_count-2),cur_depth+1)
           end do



           ! == Second term
           ! ==  - b * I^(n+2)

           temp2 = fnp_generic(leg_count,dim_nplus+2,b_pin,l_count-1,l(2:l_count),cur_depth+1)
           tmp=(leg_count-n-dim_nplus -(l_count-1._ki)-1._ki)
           tmp4(:)=(/ temp1%a, temp1%b, temp1%c /)
           tmp5(:)=(/ temp2%a, temp2%b, temp2%c /)

           if (tmp/= 0) then
              tmp4(3)=tmp4(3)-(b(l(1),b_pin) * ( tmp * tmp5(3) + 2._ki*tmp5(2) ))
              tmp4(2)=tmp4(2)-(b(l(1),b_pin) * ( tmp * tmp5(2) + 2._ki*tmp5(1) ))
              tmp4(1)=tmp4(1)-(b(l(1),b_pin) * ( tmp * tmp5(1) ))
           else
           end if
           temp1%a=tmp4(1)
           temp1%b=tmp4(2)
           temp1%c=tmp4(3)




        end if


        ! == Third term
        ! ==  - c * invS * I^(n)(... S \ {j} )


        ib = b_used
        j=0

        do while (ib /= 0)
           if (iand(ib,1)==1 ) then
              b_tmp=ibset(b_pin,j) ! equal to: b_tmp = punion( b_pin,ibset(0,j) )
              do k = 1,1 ! not symmetrised yet
                 usable = .true.
                 ! exclude zero integrals
                 do i = 1,l_count
                    if (l(i) == j .and. k /= i ) then
                       usable =.false.
                       exit
                    end if
                 end do
                 if ( usable .and. (leg_count-1>0) ) then
                    l_tmp1(1:k-1) =l(1:k-1)
                    l_tmp1(k:l_count-1) =l(k+1:l_count)
                    temp2 = fnp_generic(leg_count-1,dim_nplus,b_tmp,l_count-1,l_tmp1(1:l_count-1),cur_depth+1)
                    temp1=temp1 + inv_s(j,l(k),b_pin)*temp2
                 end if
              end do
           end if
           ib=ishft(ib,-1)
           j=j+1
        enddo

        return_val = temp1
    end function reduce_generic


    ! Test if matrix is zero

    function test_mat_zero (mat_p,used_size,b_pin)
       implicit none
       type(s_matrix_poly), intent(in) :: mat_p
       integer, intent(in) :: used_size
       integer, intent(in) :: b_pin
       integer,dimension(used_size) :: s
       logical :: test_mat_zero
       integer :: b_used
       integer :: i,j

       b_used=pminus(b_ref,b_pin)
       s = unpackb(b_used,countb(b_used))

       if (iand(mat_p%b_cmplx, b_used) .eq. 0 )  then
          do i = 1, used_size
             do j = i, used_size
                if ( .not. equal_real(mat_p%pt_real(s(i),s(j)),zero)) then
                   test_mat_zero = .false.
                   return
                end if
             end do
          end do
       else
          do i = 1, used_size
             do j = i, used_size
                if ( .not. equal_zero(mat_p%pt_cmplx(s(i),s(j)))) then
                   test_mat_zero = .false.
                   return
                end if
             end do
          end do
       end if
       test_mat_zero = .true.
       return
    end function test_mat_zero


    ! Calculate determinant assume: mat_p symmetric
recursive  function calc_determinant(mat_p,used_size,b_pin) result(detS)
       implicit none
       type(s_matrix_poly), intent(in) :: mat_p
       integer, intent(in) :: used_size
       integer, intent(in) :: b_pin
       complex(ki) :: detS

       integer,dimension(used_size) :: s
       integer :: b_used

       ! check consistency
       b_used=pminus(b_ref,b_pin)

       if (countb(b_used) /= used_size ) then
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'Assert failed: internal error in file generic_function_np.f90 (inconsistency)'
          call catch_exception(0)
          stop
       end if

       s = unpackb(b_used,countb(b_used))

       if (used_size == 1) then
          if (iand(mat_p%b_cmplx, b_used) .eq. 0 ) then
             detS= mat_p%pt_real(s(1),s(1))
          else
             detS= mat_p%pt_cmplx(s(1),s(1))
          end if
       else if (used_size == 2) then
          ! calculate det(S), use s_ij=s_ji
          if (iand(mat_p%b_cmplx, b_used) .eq. 0 ) then
             detS = mat_p%pt_real(s(1),s(1))*mat_p%pt_real(s(2),s(2)) - &
                mat_p%pt_real(s(1),s(2))*mat_p%pt_real(s(1),s(2))
          else
             detS = mat_p%pt_cmplx(s(1),s(1))*mat_p%pt_cmplx(s(2),s(2)) - &
                mat_p%pt_cmplx(s(1),s(2))*mat_p%pt_cmplx(s(1),s(2))
          end if
       else if (used_size == 3) then
          if (iand(mat_p%b_cmplx, b_used) .eq. 0 ) then
             ! calculate det(S), use s_ij=s_ji
             detS= mat_p%pt_real(s(1),s(3))**2 * mat_p%pt_real(s(2),s(2)) &
                + 2*mat_p%pt_real(s(1),s(2)) * mat_p%pt_real(s(1),s(3))*mat_p%pt_real(s(2),s(3)) &
                + mat_p%pt_real(s(1),s(2))**2 * mat_p%pt_real(s(3),s(3)) &
                + mat_p%pt_real(s(1),s(1))*( mat_p%pt_real(s(2),s(3))**2 + &
                mat_p%pt_real(s(2),s(2))*mat_p%pt_real(s(3),s(3)) )
          else
             ! calculate det(S), use s_ij=s_ji
             detS= mat_p%pt_cmplx(s(1),s(3))**2 * mat_p%pt_cmplx(s(2),s(2)) &
                + 2*mat_p%pt_cmplx(s(1),s(2)) * mat_p%pt_cmplx(s(1),s(3))*mat_p%pt_cmplx(s(2),s(3)) &
                + mat_p%pt_cmplx(s(1),s(2))**2 * mat_p%pt_cmplx(s(3),s(3)) &
                + mat_p%pt_cmplx(s(1),s(1))*( mat_p%pt_cmplx(s(2),s(3))**2 + &
                mat_p%pt_cmplx(s(2),s(2))*mat_p%pt_cmplx(s(3),s(3)) )
          end if
       else if (used_size == 4) then
          ! calculate det(S), use s_ij=s_ji
          if (iand(mat_p%b_cmplx, b_used) .eq. 1 ) then
             detS= mat_p%pt_cmplx(s(1),s(2))*mat_p%pt_cmplx(s(2),s(4))*mat_p%pt_cmplx(s(3),s(3))*mat_p%pt_cmplx(s(4),s(1)) &
               - mat_p%pt_cmplx(s(1),s(2))*mat_p%pt_cmplx(s(2),s(3))*mat_p%pt_cmplx(s(3),s(4))*mat_p%pt_cmplx(s(4),s(1)) &
               - mat_p%pt_cmplx(s(1),s(1))*mat_p%pt_cmplx(s(2),s(4))*mat_p%pt_cmplx(s(3),s(3))*mat_p%pt_cmplx(s(4),s(2)) &
               + mat_p%pt_cmplx(s(1),s(1))*mat_p%pt_cmplx(s(2),s(3))*mat_p%pt_cmplx(s(3),s(4))*mat_p%pt_cmplx(s(4),s(2)) &
               - mat_p%pt_cmplx(s(1),s(2))*mat_p%pt_cmplx(s(2),s(4))*mat_p%pt_cmplx(s(3),s(1))*mat_p%pt_cmplx(s(4),s(3)) &
               - mat_p%pt_cmplx(s(1),s(1))*mat_p%pt_cmplx(s(2),s(4))*mat_p%pt_cmplx(s(3),s(2))*mat_p%pt_cmplx(s(4),s(3)) &
               + mat_p%pt_cmplx(s(1),s(2))*mat_p%pt_cmplx(s(2),s(1))*mat_p%pt_cmplx(s(3),s(4))*mat_p%pt_cmplx(s(4),s(3)) &
               - mat_p%pt_cmplx(s(1),s(1))*mat_p%pt_cmplx(s(2),s(2))*mat_p%pt_cmplx(s(3),s(4))*mat_p%pt_cmplx(s(4),s(3)) &
               - mat_p%pt_cmplx(s(1),s(4))*(mat_p%pt_cmplx(s(2),s(3))*(mat_p%pt_cmplx(s(3),s(2))*mat_p%pt_cmplx(s(4),s(1)) &
               - mat_p%pt_cmplx(s(3),s(1))*mat_p%pt_cmplx(s(4),s(2))) + &
               mat_p%pt_cmplx(s(2),s(2))*(-(mat_p%pt_cmplx(s(3),s(3))*mat_p%pt_cmplx(s(4),s(1))) &
               + mat_p%pt_cmplx(s(3),s(1))*mat_p%pt_cmplx(s(4),s(3))) + &
               mat_p%pt_cmplx(s(2),s(1))*(mat_p%pt_cmplx(s(3),s(3))*mat_p%pt_cmplx(s(4),s(2)) &
               - mat_p%pt_cmplx(s(3),s(2))*mat_p%pt_cmplx(s(4),s(3)))) - &
               (mat_p%pt_cmplx(s(1),s(2))*(mat_p%pt_cmplx(s(2),s(3))*mat_p%pt_cmplx(s(3),s(1)) &
               - mat_p%pt_cmplx(s(2),s(1))*mat_p%pt_cmplx(s(3),s(3))) + &
               mat_p%pt_cmplx(s(1),s(1))*(-(mat_p%pt_cmplx(s(2),s(3))*mat_p%pt_cmplx(s(3),s(2))) &
               + mat_p%pt_cmplx(s(2),s(2))*mat_p%pt_cmplx(s(3),s(3))))*mat_p%pt_cmplx(s(4),s(4)) &
               - mat_p%pt_cmplx(s(1),s(3))*(mat_p%pt_cmplx(s(2),s(4))*(-(mat_p%pt_cmplx(s(3),s(2))*mat_p%pt_cmplx(s(4),s(1))) &
               + mat_p%pt_cmplx(s(3),s(1))*mat_p%pt_cmplx(s(4),s(2))) + &
               mat_p%pt_cmplx(s(2),s(2))*(mat_p%pt_cmplx(s(3),s(4))*mat_p%pt_cmplx(s(4),s(1)) &
               - mat_p%pt_cmplx(s(3),s(1))*mat_p%pt_cmplx(s(4),s(4))) + &
               mat_p%pt_cmplx(s(2),s(1))*(-(mat_p%pt_cmplx(s(3),s(4))*mat_p%pt_cmplx(s(4),s(2))) &
               + mat_p%pt_cmplx(s(3),s(2))*mat_p%pt_cmplx(s(4),s(4))))
          else
             detS= mat_p%pt_real(s(1),s(2))*mat_p%pt_real(s(2),s(4))*mat_p%pt_real(s(3),s(3))*mat_p%pt_real(s(4),s(1)) &
               - mat_p%pt_real(s(1),s(2))*mat_p%pt_real(s(2),s(3))*mat_p%pt_real(s(3),s(4))*mat_p%pt_real(s(4),s(1)) &
               - mat_p%pt_real(s(1),s(1))*mat_p%pt_real(s(2),s(4))*mat_p%pt_real(s(3),s(3))*mat_p%pt_real(s(4),s(2)) &
               + mat_p%pt_real(s(1),s(1))*mat_p%pt_real(s(2),s(3))*mat_p%pt_real(s(3),s(4))*mat_p%pt_real(s(4),s(2)) &
               - mat_p%pt_real(s(1),s(2))*mat_p%pt_real(s(2),s(4))*mat_p%pt_real(s(3),s(1))*mat_p%pt_real(s(4),s(3)) &
               - mat_p%pt_real(s(1),s(1))*mat_p%pt_real(s(2),s(4))*mat_p%pt_real(s(3),s(2))*mat_p%pt_real(s(4),s(3)) &
               + mat_p%pt_real(s(1),s(2))*mat_p%pt_real(s(2),s(1))*mat_p%pt_real(s(3),s(4))*mat_p%pt_real(s(4),s(3)) &
               - mat_p%pt_real(s(1),s(1))*mat_p%pt_real(s(2),s(2))*mat_p%pt_real(s(3),s(4))*mat_p%pt_real(s(4),s(3)) &
               - mat_p%pt_real(s(1),s(4))*(mat_p%pt_real(s(2),s(3))*(mat_p%pt_real(s(3),s(2))*mat_p%pt_real(s(4),s(1)) &
               - mat_p%pt_real(s(3),s(1))*mat_p%pt_real(s(4),s(2))) + &
               mat_p%pt_real(s(2),s(2))*(-(mat_p%pt_real(s(3),s(3))*mat_p%pt_real(s(4),s(1))) &
               + mat_p%pt_real(s(3),s(1))*mat_p%pt_real(s(4),s(3))) + &
               mat_p%pt_real(s(2),s(1))*(mat_p%pt_real(s(3),s(3))*mat_p%pt_real(s(4),s(2)) &
               - mat_p%pt_real(s(3),s(2))*mat_p%pt_real(s(4),s(3)))) - &
               (mat_p%pt_real(s(1),s(2))*(mat_p%pt_real(s(2),s(3))*mat_p%pt_real(s(3),s(1)) &
               - mat_p%pt_real(s(2),s(1))*mat_p%pt_real(s(3),s(3))) + &
               mat_p%pt_real(s(1),s(1))*(-(mat_p%pt_real(s(2),s(3))*mat_p%pt_real(s(3),s(2))) &
               + mat_p%pt_real(s(2),s(2))*mat_p%pt_real(s(3),s(3))))*mat_p%pt_real(s(4),s(4)) &
               - mat_p%pt_real(s(1),s(3))*(mat_p%pt_real(s(2),s(4))*(-(mat_p%pt_real(s(3),s(2))*mat_p%pt_real(s(4),s(1))) &
               + mat_p%pt_real(s(3),s(1))*mat_p%pt_real(s(4),s(2))) + &
               mat_p%pt_real(s(2),s(2))*(mat_p%pt_real(s(3),s(4))*mat_p%pt_real(s(4),s(1)) &
               - mat_p%pt_real(s(3),s(1))*mat_p%pt_real(s(4),s(4))) + &
               mat_p%pt_real(s(2),s(1))*(-(mat_p%pt_real(s(3),s(4))*mat_p%pt_real(s(4),s(2))) &
               + mat_p%pt_real(s(3),s(2))*mat_p%pt_real(s(4),s(4))))
          end if
      else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'Determinant calulcation for this case not implemented yet.'
          call catch_exception(0)
          stop
       end if
    end function calc_determinant


    elemental function equal_zero(z)
        implicit none
        complex(ki), intent(in) :: z
        logical :: equal_zero

        equal_zero=equal_real(real(z,ki),zero) .and. equal_real(aimag(z),zero)
     end function  equal_zero




     !
     ! Generic integral reduction by Passarino-Veltmann method
     !
     !

     recursive function reduce_pave_generic(leg_count,dim_nplus,b_pin,l_count,l,depth) result(return_val)
        implicit none
        integer, intent(in) :: leg_count
        integer, intent(in) :: dim_nplus
        integer, intent (in) :: b_pin
        integer, intent (in) :: l_count
        integer, intent (in),dimension(:) :: l
        integer, intent (in), optional :: depth
        integer :: i,j,k,ii,b_used,b_tmp,cur_depth,b_pro,mass_zero_count
        integer :: onshell_count
        type(form_factor) :: return_val, tmp, temp2
        complex(ki), dimension(3) :: ret_temp
        complex(ki), dimension(leg_count)  :: m,r,f
        integer :: n ! dimension n = 4-2epsilon TODO: float
        integer,dimension(size(l)) :: lpos
        integer,dimension(size(l)+1) :: l_tmp1
        integer,dimension(leg_count) :: s

        real(ki), dimension(leg_count-1,leg_count-1) :: z_mat_r, inv_z_mat_r
        complex(ki), dimension(leg_count-1,leg_count-1) :: z_mat_c, inv_z_mat_c
        real(ki) :: inv_error
        logical :: use_complex_z_mat,usable

        n = 4


        if (present(depth)) then
           cur_depth=depth
        else
           cur_depth=0
        end if

        b_used=pminus(b_ref,b_pin)
        b_pro=b_used

        lpos = locateb(l,b_used)


        ret_temp=0._ki



        ! some checks if everything is consistent
        if ((countb(b_used) /= leg_count ) .or. (size(l)/= l_count)) then
           write (*, *) countb(b_used),"PPPPPXXXX", size(l), l_count, b_used
           tab_erreur_par(1)%a_imprimer = .true.
           tab_erreur_par(1)%chaine = 'Assert failed: internal error in file generic_function_np.f90'
           call catch_exception(0)
           stop
        end if

        s = unpackb(b_used,countb(b_used))

         if ((l_count > 0)) then
            usable = .true.
            do j = 1, l_count
               if (lpos(j)==leg_count) then
                  i=j
                  usable = .false.
                  exit
               end if
            end do
            if ( .not. usable) then
               l_tmp1(1:i-1) = l(1:i-1)
               l_tmp1(i:l_count-1) = l(i+1:l_count)
               tmp = fnp_generic(leg_count,dim_nplus, b_pin, l_count-1,l_tmp1(1:l_count-1),cur_depth+1)
               l_tmp1(1:l_count) = l
               do j = 1, leg_count-1
                  l_tmp1(i) = s(j)
                  tmp = tmp - fnp_generic(leg_count,dim_nplus, b_pin, l_count,l_tmp1(1:l_count),cur_depth+1)
               end do
               return_val = tmp
               return
            end if
         end if

        !
        ! extract masses and momenta
        !

        mass_zero_count=0
        onshell_count=0



        if (iand(s_mat_p%b_cmplx, b_pro) .eq. 0 ) then
           do i = 1, leg_count
              m(i) = -s_mat_p%pt_real(s(i),s(i))/2._ki
              if (equal_zero(m(i))) then
                 mass_zero_count= mass_zero_count+1
              end if
           end do
           do i = 1, leg_count-1
              r(i) =  s_mat_p%pt_real(s(i),s(leg_count)) + m(i) + m(leg_count)
              f(i) = r(i) - m(i) + m(leg_count)
              if (equal_zero(r(i))) then
                 onshell_count=onshell_count+1
              end if
           end do
           r(leg_count) = 0._ki
           do i = 1, leg_count-1
              do j = 1, leg_count-1
                 if (i==j) then
                    z_mat_r(i,j) = 2._ki*real(r(i),ki)
                 else
                    z_mat_r(i,j) = -s_mat_p%pt_real(s(i),s(j)) + real(r(i) + r(j) - m(i) - m(j),ki)
                 end if
              end do
           end do
           call inverse(z_mat_r,inv_z_mat_r,inv_error)
           use_complex_z_mat = .false.
        else
           do i = 1, leg_count
              m(i) = -s_mat_p%pt_cmplx(s(i),s(i))/2._ki
              if (equal_zero(m(i))) then
                 mass_zero_count= mass_zero_count+1
              end if
           end do
           do i = 1, leg_count-1
              r(i) =  s_mat_p%pt_cmplx(s(i),s(leg_count)) + m(i) + m(leg_count)
              f(i) = r(i) - m(i) + m(leg_count)
              if (equal_zero(r(i))) then
                 onshell_count=onshell_count+1
              end if
           end do
           r(leg_count) = 0._ki
           do i = 1, leg_count-1
              do j = 1, leg_count-1
                 if (i==j) then
                    z_mat_c(i,j) = 2._ki*r(i)
                 else
                    z_mat_c(i,j) = ( -s_mat_p%pt_real(s(i),s(j)) + r(i) + r(j) - m(i) - m(j))
                 end if
              end do
           end do
           call inverse(z_mat_c,inv_z_mat_c,inv_error)
           use_complex_z_mat = .true.
        end if

        if ( (leg_count==2) .and. (mass_zero_count /= 2) .and. ( onshell_count==1 )  &
             .and. ( l_count==0 .or.  mass_zero_count == 1 .or. equal_zero(m(1)-m(2)) )) then
           return_val = f2p_ndim_0p_generic(m(1),m(2),dim_nplus,l_count,lpos);
           return
        end if

        !
        ! Reduction case 1
        !

!        if ((dim_nplus >= 2) .and. ( l_count== 0 ) ) then
        if ((dim_nplus >= 2) ) then

           ! only leg_count == 2,3,4 implemented yet
           if ( leg_count /= 3 .and. leg_count /= 2 .and. leg_count /= 4) then
              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'Assert failed: wrong arguments (not implemented) in file generic_function_np.f90'
              call catch_exception(0)
              stop
           end if

           ! check for special cases

           if (leg_count==3) then
              if ( (use_complex_z_mat .and.  equal_zero(-z_mat_c(1,2)+(z_mat_c(1,1)+z_mat_c(2,2))/2._ki )) &
                .or. (.not. use_complex_z_mat .and. equal_real(-z_mat_r(1,2)+(z_mat_r(1,1)+z_mat_r(2,2))/2._ki,0._ki))) then
                 onshell_count=onshell_count+1
              end if
           end if

           if (leg_count==3 .and. mass_zero_count==leg_count .and. onshell_count>0 .and. dim_nplus ==4 .and. l_count==0 ) then
              if (use_complex_z_mat) then
                  return_val =  f3p_np4_m0(r(1),r(2), -z_mat_c(1,2)+(z_mat_c(1,1)+z_mat_c(2,2))/2._ki)
              else
                  return_val = f3p_np4_m0(r(1),r(2),cmplx(-z_mat_r(1,2)+(z_mat_r(1,1)+z_mat_r(2,2))/2._ki,0._ki,ki))
              end if
              return
           end if


           ! hep-ph/0509141 5.10 or arXiv:1105.0920 A.8
           ! used factor 1/(2(D+P-1-N))=  1/(4*(2-eps)) = 1/8 + eps/16 + eps^2/32
           !                                                           + O(eps^3)
           !                            for N=3, P=3

           b_tmp = ibset(b_pin,s(leg_count))

           usable = .true.
            do j = 1, l_count
               if ( lpos(j) > leg_count-1 ) then
                  usable = .false.
                  exit
               end if
            end do

           if (usable) then
              tmp = (-1._ki)**(l_count)*(-0.5_ki)**((dim_nplus-2)/2) * &
                 fnp_generic(leg_count-1,dim_nplus-2,b_tmp,l_count,l,cur_depth+1)
           else
              tmp = 0
           end if

           temp2 = tmp
           tmp = tmp + 2._ki*m(leg_count)*(-1._ki)**(l_count)*(-.5_ki)**((dim_nplus-2)/2) * &
              fnp_generic(leg_count,dim_nplus-2,b_pin,l_count,l,cur_depth+1)



           l_tmp1 = 0
           l_tmp1(2:) = l(:)
           do j = 1, leg_count-1
              l_tmp1(1)=s(j)
              tmp = tmp + f(j)* (-1._ki)**(l_count+1)*(-.5_ki)**((dim_nplus-2)/2)* &
               fnp_generic(leg_count,dim_nplus-2,b_pin,l_count+1,l_tmp1(1:l_count+1),cur_depth+1)

           end do

           k = (n+dim_nplus-2)+(l_count+2)-1-leg_count ! (D+P-1-N) without - 2*eps
           if (k==0) then

              tab_erreur_par(1)%a_imprimer = .true.
              tab_erreur_par(1)%chaine = 'Assert failed: unknown case in generic_function_np.f90'
              call catch_exception(0)
              stop
           end if

           ret_temp = (/ tmp%a, tmp%b, tmp%c /)



           ret_temp(3) = ret_temp(3)/(2*k) + ret_temp(2)/k/k + ret_temp(1)*2._ki/k/k/k
           ret_temp(2) = ret_temp(2)/(2*k) + ret_temp(1)/k/k
           ret_temp(1) = ret_temp(1)/(2*k)

           return_val = ret_temp*(-2._ki)**(dim_nplus/2)*(-1)**(l_count) ! factor because of g_mu,nu encoded by dim_nplus


        !
        !  Reduction case 2
        !
        else if (l_count > 0) then
           ! Formula 5.8 and 5.11 from 0509141v2

           temp2=0._ki


           do i = 1, leg_count-1
              usable = .true.
              do j = 2, l_count
                 if (l(j) == s(i)) then
                    usable = .false.
                    exit
                 end if
              end do

              if (usable) then
                 b_tmp = ibset(b_pin,s(i))
                 l_tmp1(1:l_count-1) = l(2:l_count)
                 tmp  = (-1._ki)**(l_count-1) * (-0.5_ki)**((dim_nplus)/2) * &
                    fnp_generic(leg_count-1,dim_nplus,b_tmp,l_count-1, l_tmp1(1:l_count-1),cur_depth+1)
              else
                 tmp = 0._ki
              end if


              b_tmp = ibset(b_pin,s(leg_count))

              k=0
              usable = .true.
              do j= 2, l_count
                if (lpos(j) == leg_count) then
                  usable = .false.
                  exit
                end if
              end do

              if (usable) then
                 tmp = tmp - (-1._ki)**(l_count-1) * (-.5_ki)**((dim_nplus)/2) &
                    * fnp_generic(leg_count-1,dim_nplus,b_tmp,l_count-1,l(2:l_count),cur_depth+1)
              end if

              tmp = tmp - f(i) *  (-1._ki)**(l_count-1) *  (-.5_ki)**((dim_nplus)/2) &
                     * fnp_generic(leg_count,dim_nplus,b_pin,l_count-1,l(2:),cur_depth+1)

              ! tmp contains now S_ki2..iP

              l_tmp1=0
              do ii = 2, l_count
                if (i==lpos(ii)) then
                    !k=0
                    !do j = 2, l_count
                    !   if(l(j) /= l(ii)) then
                    !     k=k+1
                    !     l_tmp1(k) = l(j)
                    !  end if
                    !end do

                    l_tmp1(1:ii-2) = l(2:ii-1)
                    l_tmp1(ii-1:l_count-2) = l(ii+1:l_count)
                    k = l_count-2
                    tmp = tmp - 2._ki*(-1._ki)**(k)*(-.5_ki)**((dim_nplus+2)/2) &
                       *fnp_generic(leg_count,dim_nplus+2,b_pin,k,l_tmp1(1:k),cur_depth+1)
                 end if
              end do

              if (use_complex_z_mat) then
                 temp2 = temp2 + inv_z_mat_c(lpos(1),i)*tmp
              else
                 temp2 = temp2 + inv_z_mat_r(lpos(1),i)*tmp
              end if
           end do
           return_val=temp2*(-1._ki)**(l_count) * (-2._ki)**(dim_nplus/2)

       else
           tab_erreur_par(1)%a_imprimer = .true.
           tab_erreur_par(1)%chaine = 'Assert failed: Not implemented yet. generic_function_np.f90'
           call catch_exception(0)
           stop
        end if

     end function reduce_pave_generic




    function f3p_np4_m0(r1s,r2s,r1mr2s)
       complex(ki),intent (in) :: r1s,r2s,r1mr2s
       type(form_factor) :: f3p_np4_m0
       complex(ki), dimension(3) :: ret_temp

       complex(ki) :: Q,P
       logical :: zr1s,zr2s,zr1mr2s

       ! case: one scale or two scales

       Q = 0._ki
       P = 0._ki




       zr1s=equal_zero(r1s)
       zr2s=equal_zero(r2s)
       zr1mr2s=equal_zero(r1mr2s)

       if ( zr1s .and. zr2s .and. (.not. zr1mr2s) ) then
          Q = r1mr2s
       else if( zr2s .and. zr1mr2s .and. (.not. zr1s) ) then
          Q = r1s
       else if ( zr1s .and. zr1mr2s .and. (.not. zr2s) ) then
          Q = r2s
       else if ( (.not. zr1s) .and. (.not. zr2s) .and. zr1mr2s ) then
          Q = r1s
          P = r2s
       else if ( zr1s .and. (.not. zr2s) .and. (.not. zr1mr2s) ) then
          Q = r2s
          P = r1mr2s
       else if ( (.not. zr1s) .and. zr2s .and. (.not. zr1mr2s) ) then
          Q = r1s
          P = r1mr2s
       else if ( zr1s .and. zr1mr2s .and. zr2s ) then
          ! no scale, return zero
          ret_temp=0
          f3p_np4_m0=ret_temp
          return
       else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'Assert failed: f3p_np4_m0: three scales,  not implented'
          call catch_exception(0)
          stop
       end if
       ret_temp = 0._ki


       if ( P == 0 .and. Q /= 0 ) then
          if (rat_or_tot_par%rat_selected) then
             ret_temp(3) =  1._ki/144._ki*Q*( - 19._ki)
             ret_temp(2) =  Q / (-24._ki)
          else if (rat_or_tot_par%tot_selected) then
             ret_temp(3) =  1._ki/144._ki*Q*(6._ki*z_log(-Q/mu2_scale_par,-1._ki) - 19._ki)
             ret_temp(2) =  Q / (-24._ki)
          end if

       else if ( P /= 0 .and. Q /= 0 ) then
          if (rat_or_tot_par%rat_selected) then
             ret_temp(2) = ( Q + P ) / (-24._ki)
             ret_temp(3) = (-19._ki* Q**2 +  19._ki*P**2 ) / 144._ki / (Q-P)
          else if (rat_or_tot_par%tot_selected) then
             ret_temp(2) = ( Q + P ) / (-24._ki)
             ret_temp(3) = (-19._ki* Q**2 + 6._ki * Q**2 *z_log(-Q/mu2_scale_par,-1._ki)+ &
                19._ki*P**2 - 6._ki*P**2*z_log(-P/mu2_scale_par,-1._ki) ) / 144._ki / (Q-P)
          end if
       else
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'Assert failed: f3p_np4_m0: internal error.'
          call catch_exception(0)
          stop
       end if



       f3p_np4_m0 = ret_temp

    end function f3p_np4_m0

    function f1p_ndim_generic(mass1, dim_nplus) result (return_val)
        implicit none
        complex(ki), intent(in) :: mass1;
        integer, intent(in) :: dim_nplus
        type(form_factor) :: return_val
        complex(ki), dimension(3) :: ret_temp
        integer :: i,j,y
        real (ki) :: r_tmp1, r_tmp2
        ! Formula (3.1) from hep-ph/0509141


        ret_temp = 0

        if ( (equal_real(real(mass1,ki),zero)) .and.  (equal_real(aimag(mass1),zero))) then ! all massless tadpoles are zero
           return_val=ret_temp;
           return
        end if

        !
        ! !  explicit formulars

        !   if (l_count==0 .and. dim_nplus == 2) then
        !      ret_temp(2) = -0.5_ki*mass1*mass1
        !      if (rat_or_tot_par%tot_selected) then
        !         ret_temp(3) = mass1*mass1 * (-3._ki + 2._ki*z_log(mass1/mu2_scale_par,-1._ki) ) / 4._ki
        !      else
        !         ret_temp(3) = mass1*mass1 * (-3._ki ) / 4._ki
        !      end if

        !      return_val=ret_temp*(1._ki)
        !      return
        !   end if

        !   if (l_count==0 .and. dim_nplus == 4) then
        !      ret_temp(2) =  mass1**3._ki / 6._ki
        !      if (rat_or_tot_par%tot_selected) then
        !         ret_temp(3) = mass1**3._ki * (11._ki - 6._ki *z_log(mass1/mu2_scale_par,-1._ki)) / 36._ki
        !      else
        !         ret_temp(3) = mass1**3._ki * (11._ki) / 36._ki
        !      end if
        !      return_val=ret_temp*(1._ki)
        !      return
        !   end if

        !   if (l_count==0 .and. dim_nplus == 6) then
        !      ret_temp(2) =- (mass1**4) / 24._ki
        !      if (rat_or_tot_par%tot_selected) then
        !         ret_temp(3) = mass1**4 * (-25._ki + 12._ki *z_log(mass1/mu2_scale_par,-1._ki)) / 288._ki
        !      else
        !         ret_temp(3) = mass1**4 * (-25._ki) / 288._ki
        !      end if
        !      return_val=ret_temp*(1._ki)
        !      return
        !   end if

        !   if (l_count==0 .and. dim_nplus == 8) then
        !      ret_temp(2) = (mass1**5) / 120._ki
        !      if (rat_or_tot_par%tot_selected) then
        !         ret_temp(3) = mass1**5 * (137._ki - 60._ki *z_log(mass1/mu2_scale_par,-1._ki)) / 7200._ki
        !      else
        !         ret_temp(3) = mass1**5 * (137._ki) / 7200._ki
        !      end if
        !      return_val=ret_temp*(1._ki)
        !      return
        !   end if


         ! Formula (3.4) from hep-ph/0509141

        ! we need  A0 (equal to f1p)
        ret_temp(2) = mass1
        if (rat_or_tot_par%tot_selected) then
           ret_temp(3) = mass1*(1._ki - z_log(mass1/mu2_scale_par,-1._ki))
        else if (rat_or_tot_par%rat_selected) then
           ret_temp(3) = mass1* (1._ki)
        end if


        y = dim_nplus

        if (y > 0) then ! To be tested, y is assumed even
                 i = y/2
                 r_tmp1=1._ki
                 r_tmp2=0._ki
                 do j = 2, i+1
                         r_tmp1=r_tmp1*j
                 end do
                 do j = 1, i
                         r_tmp2 = r_tmp2+1._ki/(1._ki+j)
                 end do
                 ret_temp(3) = ret_temp(3) + mass1*r_tmp2
                ! ret_temp = mass1**i / (2._ki**i * r_tmp1)  * (ret_temp)  * (-2._ki)**(i)
                 ret_temp = mass1**i / (r_tmp1)  * (ret_temp) * (-1._ki)**(i) ! equal to line above

                ! Factor (-2._ki)**(i) needed as formular is for A_(00...)

         end if
         return_val = ret_temp

    end function f1p_ndim_generic


    pure function harmonic_number(n) result (return_val)
       implicit none
       integer, intent(in) :: n;
       integer :: i
       real(ki) :: return_val
       return_val = 0
       do i = 1, n
          return_val = return_val + 1._ki/i
       end do
    end function


    ! l_pos here only 1 or 2 allowed, corresponding to m1,m2
    function  f2p_ndim_0p_generic(m1,m2,dim_nplus,l_count,l_pos) result (return_val)
        implicit none
        complex(ki) :: m1, m2, diffm;

        integer, intent(in) :: dim_nplus
        integer, intent (in) :: l_count
        integer, intent (in),dimension(:) :: l_pos
        type(form_factor) :: return_val
        integer :: leq1_count ,k1,k2;
        real(ki) :: k1fak,k2fak,ksump1fak, harmonic1,harmonic2
        logical :: zm1,zm2,zdiffm
        integer :: i
        complex(ki), dimension(3) :: ret_temp


        diffm = m1-m2
        zm1 = equal_zero(m1)
        zm2 = equal_zero(m2)
        zdiffm = equal_zero(diffm)

        return_val = 0
        ret_temp = 0


        if (zm1 .and. zm2) then
           return
        end if

        if (l_count==0) then
           if ( .not. zm1 .and. zm2) then
              ! Formula (A.9) hep/0504267
              return_val = f1p_ndim_generic(m1,dim_nplus)/m1
              return
          else if ( zm1 .and. .not. zm2) then
              return_val = f1p_ndim_generic(m2,dim_nplus)/m2
              return
          else if ( zdiffm )  then
             if ( dim_nplus>=2) then
                return_val = -1*f1p_ndim_generic(m1,dim_nplus-2) ! A.10
             else
              ret_temp(2) = 1
              if (rat_or_tot_par%tot_selected) then
                 ret_temp(3) = (z_log(m1/mu2_scale_par,-1._ki))
              else if (rat_or_tot_par%rat_selected) then
                 ret_temp(3) = 0
              end if
              return_val = ret_temp
             end if
             return
          else
             ! Formula (A.9) hep/0504267
             return_val = f1p_ndim_generic(m2,dim_nplus-2) - f1p_ndim_generic(m1,dim_nplus-2) / (m2-m1)  ! A.9
             return
          end if

       else ! l_count > 0

          leq1_count=0
          do i = 1, l_count
             if (l_pos(i) == 1) then
                leq1_count = leq1_count+1
             endif
          end do

          if (zm1 .and. .not. zm2) then
             ! swap m1 and m2 and in l_pos 1<->2
             zm1=zm2
             zm2=.true.
             m1=m2
             m2=0
             leq1_count = l_count-leq1_count
          end if

          k1= leq1_count
          k2= l_count- leq1_count

          k1fak = 1;
          do i = 1, leq1_count
             k1fak = k1fak*i
          end do
          k2fak = 1
          do i = 2, l_count-leq1_count
             k2fak = k2fak*i
          end do
          ksump1fak = 1
          do i = 2,l_count+1
             ksump1fak=ksump1fak*i
          end do


         if ( zdiffm )  then
            if (dim_nplus==0) then
               ret_temp(2) = k1fak * k2fak / ksump1fak
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = ret_temp(2) * ( -z_log(m1/mu2_scale_par,-1._ki))
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = 0
               end if
            else if (dim_nplus==2) then
               ret_temp(2) = - k1fak * k2fak / ksump1fak * m1
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = -ret_temp(2) * (  z_log(m1/mu2_scale_par,-1._ki) - 1._ki)
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -1._ki)
               end if
            else if (dim_nplus==4) then
               ret_temp(2) = k1fak * k2fak / ksump1fak * m1 * m1  / 2._ki
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = ret_temp(2) * ( 3._ki - 2._ki * z_log(m1/mu2_scale_par,-1._ki)) / 2._ki
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = ret_temp(2) * ( 3._ki ) / 2._ki
               end if
            else if (dim_nplus==6) then
               ret_temp(2) = - k1fak * k2fak / ksump1fak * m1**3  / 6._ki
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -11._ki + 6._ki*z_log(m1/mu2_scale_par,-1._ki)) / 6._ki
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -11._ki ) / 36._ki
               end if
            else if (dim_nplus==8) then
               ret_temp(2) = k1fak * k2fak / ksump1fak * m1**4  / 24._ki
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = ret_temp(2) * ( 25._ki - 12._ki*z_log(m1/mu2_scale_par,-1._ki)) / 12._ki
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = ret_temp(2) * ( 25._ki ) / 12._ki
               end if
            else if (dim_nplus==10) then
               ret_temp(2) = - k1fak * k2fak / ksump1fak * m1**5  / 120._ki
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -137._ki + 60._ki*z_log(m1/mu2_scale_par,-1._ki)) / 60._ki
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -137._ki ) / 60._ki
               end if
            else
               tab_erreur_par(1)%a_imprimer = .true.
               tab_erreur_par(1)%chaine = 'Assert failed: f2p_ndim_0p_generic case (m1=m2,dim>14,l_count>0) not implemented.'
               call catch_exception(0)
               stop
            end if
            return_val = ret_temp
         else if ( (.not. zm1 .and. zm2)) then

            do i = l_count+2, l_count+1 + dim_nplus/2
               ksump1fak = ksump1fak * i
            end do
            do i = leq1_count+1, leq1_count + dim_nplus/2
               k1fak = k1fak * i
            end do


            if (dim_nplus==0) then
               ret_temp(2) = k1fak * k2fak / ksump1fak
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = - ret_temp(2) * ( harmonic_number(k1) - harmonic_number(1+k1+k2) + &
                     z_log(m1/mu2_scale_par,-1._ki))
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = - ret_temp(2) * ( harmonic_number(k1) - harmonic_number(1+k1+k2))
               end if
            else if (dim_nplus==2) then
               ret_temp(2) = - k1fak * k2fak / ksump1fak * m1

               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = -ret_temp(2) * (-1._ki + harmonic_number(k1+1) - harmonic_number(2+k1+k2) + &
                     &  z_log(m1/mu2_scale_par,-1._ki))
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = -ret_temp(2) * (-1._ki + harmonic_number(k1+1) - harmonic_number(2+k1+k2))
               end if
            else if (dim_nplus==4) then
               ret_temp(2) = k1fak * k2fak / ksump1fak * m1 * m1  / 2._ki
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = ret_temp(2) * ( 3._ki - 2._ki*(harmonic_number(2+k1) -  harmonic_number(3+k1+k2) + &
                    & z_log(m1/mu2_scale_par,-1._ki))) / 2._ki
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = ret_temp(2) * ( 3._ki - 2._ki*(harmonic_number(2+k1) - harmonic_number(3+k1+k2))) / 2._ki
               end if
            else if (dim_nplus==6) then
               ret_temp(2) = - k1fak * k2fak / ksump1fak * m1**3  / 6._ki
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -11._ki + 6._ki*(harmonic_number(3+k1) - harmonic_number(4+k1+k2) + &
                   &  z_log(m1/mu2_scale_par,-1._ki))) / 6._ki
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -11._ki + 6._ki*(harmonic_number(3+k1) - harmonic_number(4+k1+k2))) / 6._ki
               end if
            else if (dim_nplus==8) then
               ret_temp(2) = k1fak * k2fak / ksump1fak * m1**4  / 24._ki
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = ret_temp(2) * ( 25._ki - 12._ki*(harmonic_number(4+k1) - harmonic_number(5+k1+k2) + &
                   &  z_log(m1/mu2_scale_par,-1._ki))) / 12._ki
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = ret_temp(2) * ( 25._ki - 12._ki*(harmonic_number(4+k1) - harmonic_number(5+k1+k2))) / 12._ki
               end if
            else if (dim_nplus==10) then
               ret_temp(2) = - k1fak * k2fak / ksump1fak * m1**5  / 120._ki
               if (rat_or_tot_par%tot_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -137._ki + 60._ki*(harmonic_number(5+k1) - harmonic_number(6+k1+k2) &
                    & + z_log(m1/mu2_scale_par,-1._ki))) / 60._ki
               else if (rat_or_tot_par%rat_selected) then
                  ret_temp(3) = -ret_temp(2) * ( -137._ki + 60._ki*(harmonic_number(5+k1) - harmonic_number(6+k1+k2))) / 60._ki
               end if
            else
               tab_erreur_par(1)%a_imprimer = .true.
               tab_erreur_par(1)%chaine = 'Assert failed: f2p_ndim_0p_generic case (m2=0 or m1=0,dim>14,l_count>0) not implemented.'
               call catch_exception(0)
               stop
            end if

            return_val = ret_temp

          else
               tab_erreur_par(1)%a_imprimer = .true.
               tab_erreur_par(1)%chaine = 'Assert failed: f2p_ndim_0p_generic case (m1 != m2,l_count>0) not implemented.'
               call catch_exception(0)
               stop
         end if

      end if
    end function f2p_ndim_0p_generic




end module generic_function_np
