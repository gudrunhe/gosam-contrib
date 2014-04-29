module mtens
   use precision, only: ki
   use constants, only: czip, cone, im, two, half, zip,ctwo,one,six,three, &
                        chaf
   use options, only: verbosity
   use mfunctions
   use mcgs
   implicit none
   private
   save

   complex(ki), dimension(0:209,4) :: qg

   logical :: qg_list_is_initialized = .false.

   integer :: myrank, mynleg

   public :: tensor_reconstruction
   public :: numetens

   private :: ki, czip, cone, im, two, half, zip, verbosity

contains

   subroutine tensor_reconstruction(numeval,nleg,rank)

      integer, intent(in) :: nleg, rank
      complex(ki), dimension(0:209) :: xneval,known
      complex(ki), dimension(0:9) :: dxneval,dknown
      integer :: ig

      interface
         function     numeval(ncut, Q, mu2)
           use precision
           implicit none
           integer, intent(in) :: ncut
           complex(ki), dimension(4), intent(in) :: Q
           complex(ki), intent(in) :: mu2
           complex(ki) :: numeval
         end function numeval
      end interface

      mynleg=nleg
      myrank=rank

      if (.not. qg_list_is_initialized) call init_qg_list

      call init_cg(numeval,rank,xneval,dxneval)
 
      if (rank.eq.1) then 
         do ig=1,4
            cg(ig) =  xneval(ig)
         enddo
         if (verbosity .ge. 3) then
            do ig=0,4
               print*, ig, cg(ig)
            enddo
         end if

      elseif (rank.eq.2) then 
         call solve_system_rk2(xneval,known) 
         if (nleg.le.3) call solve_system_extra(dxneval,dknown,nleg,rank)
       elseif (rank.eq.3) then
         call solve_system_rk3(xneval,known)
         if (nleg.le.3) call solve_system_extra(dxneval,dknown,nleg,rank)
      elseif (rank.eq.4) then
         call solve_system_rk4(xneval,known) 
         if (nleg.le.4) call solve_system_extra(dxneval,dknown,nleg,rank)
      elseif (rank.eq.5) then
         call solve_system_rk5(xneval,known)
      elseif (rank.ge.6) then
         call solve_system_rk6(xneval,known)
      elseif (rank.ge.7) then
         print*, 'rank',rank,' not implemented'
         stop
      endif

!    print*, 'allora=',- (&
!            &  +  cg(35) + cg(36) + cg(37) + cg(38)&
!            &  + (cg(51) - cg(56) + cg(52) - cg(53)&
!            &  + +cg(54) - cg(55))/3.0_ki ) *3.0_ki/8.0_ki&
!
!            &  + ((cg(35) + cg(36) + cg(37) + cg(38))&
!            &  +  (cg(51) + cg(56) + cg(52) + cg(53)&
!            &  +  +cg(54) + cg(55))/3.0_ki ) *3.0_ki/4.0_ki



!    print*, 'allora=', (&
!            &  +  cg(35) + cg(36) + cg(37) + cg(38)&
!            &  + (cg(51) - cg(56) + cg(52) - cg(53)&
!            &  + +cg(54) - cg(55))/3.0_ki ) *1.0_ki/8.0_ki,&
!            abs((&
!            &  +  cg(35) + cg(36) + cg(37) + cg(38)&
!            &  + (cg(51) - cg(56) + cg(52) - cg(53)&
!            &  + +cg(54) - cg(55))/3.0_ki ) *1.0_ki/8.0_ki)

!    print*, 'allora=',+ (&
!            &  +  cg(35) + cg(36) + cg(37) + cg(38)&
!            &  + (cg(51) + cg(56) + cg(52) + cg(53)&
!            &  + +cg(54) + cg(55)) ) /four!&
!            &         - (&
!            &  +  cg(35) + cg(36) + cg(37) + cg(38)&
!            &  + (cg(51) + cg(56) + cg(52) + cg(53)&
!            &  + +cg(54) + cg(55))) /eight

!,&
!            abs((&
!            &  +  cg(35) + cg(36) + cg(37) + cg(38)&
!            &  + (cg(51) + cg(56) + cg(52) + cg(53)&
!            &  + +cg(54) + cg(55)) ) /four)


      
    end subroutine tensor_reconstruction


    subroutine solve_system_extra(dxneval,dknown,nleg, rank)
      implicit none
      complex(ki), dimension(0:9) :: dxneval,dknown
      integer :: ig
      !complex(ki) :: mu2
      integer :: nleg, rank, diff
     
      diff = nleg-rank
      
      if ((nleg.eq.2).and.(rank.eq.2)) then
         cgx(1) = dxneval(0) ! termine mu2
      endif

      if ((nleg.eq.3).and.(rank.ge.2)) then
         cgx(1) = dxneval(0)! termine mu2 
         if (rank.eq.3) then
            cgx(3)=dxneval(1)-cgx(1)-sub1(1)-sub2(1)-sub3(1) ! termini mu2*q1
            cgx(4)=dxneval(2)-cgx(1)-sub1(2)-sub2(2)-sub3(2) ! termini mu2*q2
            cgx(5)=dxneval(3)-cgx(1)-sub1(3)-sub2(3)-sub3(3) ! termini mu2*q3
            cgx(6)=dxneval(4)-cgx(1)-sub1(4)-sub2(4)-sub3(4) ! termini mu2*q4
         endif
      endif

      if ((nleg.eq.4).and.(rank.eq.4)) then
! mu2 and mu4
         cgx(2)= (dxneval(0)+ dxneval(5))/two ! termine mu2**2 
         cgx(1)= (dxneval(0)- dxneval(5))/two ! termine mu2 
         do ig =1,4
            !MU2=1
            dknown(ig)=dxneval(ig)-cgx(2)-cgx(1)-sub1(ig)-sub2(ig)-sub3(ig)-sub4(ig) 
            dknown(ig+5)=dxneval(ig+5)-cgx(2)-cgx(1)-sub1(ig+4)-sub2(ig+4)-sub3(ig+4)-sub4(ig+4)
         enddo
! mu2*q
         cgx(3) = (dknown(1) - dknown(6))/two
         cgx(4) = (dknown(2) - dknown(7))/two
         cgx(5) = (dknown(3) - dknown(8))/two
         cgx(6) = (dknown(4) - dknown(9))/two
! mu2*q*q
         cgx(7) = (dknown(1) + dknown(6))/two
         cgx(8) = (dknown(2) + dknown(7))/two
         cgx(9) = (dknown(3) + dknown(8))/two
         cgx(10)= (dknown(4) + dknown(9))/two
      endif

!      do ig =1,10
!         print*, 'cgx',ig,'=',cgx(ig)
!      enddo
!      stop     
 end subroutine solve_system_extra

subroutine solve_system_rk2(xneval,known)
   implicit none
   complex(ki), dimension(0:209),intent(inout) :: xneval,known
   integer :: ig

   do ig=0,3
      cg(1+ig) = (xneval(1+ig)-xneval(5+ig))/two
      cg(5+ig) = (xneval(1+ig)+xneval(5+ig))/two
   enddo
   cg(9)  = xneval(9)  - sub1(9)- sub2(9)
   cg(10) = xneval(10) - sub1(10)- sub2(10)
   cg(11) = xneval(11) - sub1(11)- sub2(11)
   cg(12) = xneval(12) - sub1(12)- sub2(12)
   cg(13) = xneval(13) - sub1(13)- sub2(13)
   cg(14) = xneval(14) - sub1(14)- sub2(14)           
         
   if (verbosity .ge. 1) then
      do ig=0,14
         print*, ig, cg(ig)
      enddo
   end if
end subroutine solve_system_rk2

subroutine solve_system_rk3(xneval,known)
   implicit none
   complex(ki), dimension(0:209),intent(inout) :: xneval,known
   integer :: ig, igs
   integer, dimension(0:5,3) :: is2

   do ig=0,3
      cg(1+ig)  = ( - 1.0_ki/2.0_ki*xneval(ig+1) - 1.0_ki/6.0_ki*xneval(ig+5) &
           + 8.0_ki/3.0_ki*xneval(ig+15))
      cg(5+ig)  = (1.0_ki/2.0_ki*xneval(ig+1) + 1.0_ki/2.0_ki*xneval(ig+5))
      cg(15+ig) = (xneval(ig+1) - 1.0_ki/3.0_ki*xneval(ig+5) - &
           8.0_ki/3.0_ki*xneval(ig+15))
   enddo 

   !1-2:   9,19,22
   !1-3:  10,20,23
   !1-4:  11,21,24
   !2-3:  12,25,27
   !2-4:  13,26,28
   !3-4:  14,29,30

   is2(0,:) =  (/ 9 ,19,22 /)
   is2(1,:) =  (/ 10,20,23 /)
   is2(2,:) =  (/ 11,21,24 /)
   is2(3,:) =  (/ 12,25,27 /)
   is2(4,:) =  (/ 13,26,28 /)
   is2(5,:) =  (/ 14,29,30 /)

   do igs=9,30
      known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)
   enddo
  
   do ig=0,5
      cg(is2(ig,1))= - 3.0_ki/2.0_ki*known(9+ig) - 1.0_ki/2.0_ki*known(19+ig) +&
           4.0_ki*known(25+ig)
      cg(is2(ig,2)) = 2.0_ki*known(9+ig) - 4.0_ki*known(25+ig)
      cg(is2(ig,3)) = 1.0_ki/2.0_ki*known(9+ig) + 1.0_ki/2.0_ki*known(19+ig)
   enddo
! finally k=3
   do igs=31,34
      known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)
      cg(igs) =  known(igs)
   enddo

   if (verbosity .ge. 1) then
      do ig=0,34
         print*, ig, cg(ig)
      enddo
   end if

end subroutine solve_system_rk3

subroutine solve_system_rk4(xneval,known)
   implicit none
   complex(ki), dimension(0:209),intent(inout) :: xneval,known
   integer :: ig, igs
   integer, dimension(0:5,6) :: is2
   integer, dimension(0:3,4) :: is3


! one q
   do ig=0,3 
      cg(1+ig)  = ( - 1.0_ki/6.0_ki*xneval(ig+1) + 1.0_ki/6.0_ki*xneval(ig+5) + 4.0_ki/&
           & 3.0_ki*xneval(ig+15) - 4.0_ki/3.0_ki*xneval(ig+35))
      
      cg(5+ig)  = ( - 1.0_ki/6.0_ki*xneval(ig+1) - 1.0_ki/6.0_ki*xneval(ig+5) + 8.0_ki/&
           & 3.0_ki*xneval(ig+15) + 8.0_ki/3.0_ki*xneval(ig+35))
      
      cg(15+ig) = (2.0_ki/3.0_ki*xneval(ig+1) - 2.0_ki/3.0_ki*xneval(ig+5) - 4.0_ki/3.0_&
           &ki*xneval(ig+15) + 4.0_ki/3.0_ki*xneval(ig+35))
      
      cg(35+ig) = (2.0_ki/3.0_ki*xneval(ig+1) + 2.0_ki/3.0_ki*xneval(ig+5) - 8.0_ki/3.0_&
           &ki*xneval(ig+15) - 8.0_ki/3.0_ki*xneval(ig+35))
   enddo

   ! two q
   do igs=1,56
      known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-sub4(igs)
   enddo

   !1-2:   9,19,22,39,42,51
   !1-3:  10,20,23,40,43,52
   !1-4:  11,21,24,41,44,53
   !2-3:  12,25,27,45,47,54
   !2-4:  13,26,28,46,48,55
   !3-4:  14,29,30,49,50,56

   is2(0,:) =  (/ 9 ,19,22,39,42,51 /)
   is2(1,:) =  (/ 10,20,23,40,43,52 /)
   is2(2,:) =  (/ 11,21,24,41,44,53 /)
   is2(3,:) =  (/ 12,25,27,45,47,54 /)
   is2(4,:) =  (/ 13,26,28,46,48,55 /)
   is2(5,:) =  (/ 14,29,30,49,50,56 /)


    do ig=0,5
       cg(is2(ig,1))  = 1.0_ki/6.0_ki*known(ig+9) + 1.0_ki/2.0_ki*known(ig+19) - 2.0_ki/3.0&
            &_ki*known(ig+25) - 2.0_ki/3.0_ki*known(ig+39) - 16.0_ki/3.0_ki*known(ig+&
            & 45)

       cg(is2(ig,2)) = - 1.0_ki/4.0_ki*known(ig+19) + known(ig+25) + known(ig+39) - 1.0_ki/&
            & 4.0_ki*known(ig+51)
       
       cg(is2(ig,3)) = 1.0_ki/2.0_ki*known(ig+9) + 1.0_ki/4.0_ki*known(ig+19) - known(ig+25)&
            &  - known(ig+39) - 1.0_ki/4.0_ki*known(ig+51)

       cg(is2(ig,4)) = 4.0_ki/3.0_ki*known(ig+9) - 4.0_ki*known(ig+25) - 4.0_ki/3.0_ki*&
            & known(ig+39)

       cg(is2(ig,5)) = - known(ig+9) - 3.0_ki/4.0_ki*known(ig+19) + 11.0_ki/3.0_ki*&
            & known(ig+25) + known(ig+39) + 16.0_ki/3.0_ki*known(ig+45) + 1.0_ki/4.0_ki&
            & *known(ig+51)
       
       cg(is2(ig,6)) = 1.0_ki/4.0_ki*known(ig+19) + known(ig+25) + known(ig+39) + 1.0_ki/4.0_&
            &ki*known(ig+51)
    enddo

    ! three q
    do igs = 57,68
       known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-sub4(igs)
    enddo
    do igs = 31,34
       known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-sub4(igs)
    enddo

    !1-2-3: 31,57,60,62
    !1-2-4: 32,58,61,63
    !1-3-4: 33,59,64,65
    !2-3-4: 34,66,67,68
   is3(0,:) =  (/ 31,57,60,62 /)
   is3(1,:) =  (/ 32,58,61,63 /)
   is3(2,:) =  (/ 33,59,64,65 /)
   is3(3,:) =  (/ 34,66,67,68 /)

   
   do ig=0,3
      cg(is3(ig,1))  =  - 1.0_ki/2.0_ki*known(ig+31) - 1.0_ki/2.0_ki*known(ig+57) - 1.0_&
           &ki/2.0_ki*known(ig+61) - 1.0_ki/2.0_ki*known(ig+65)
      
      cg(is3(ig,2))  = 1.0_ki/2.0_ki*known(ig+31) + 1.0_ki/2.0_ki*known(ig+57)
      
      cg(is3(ig,3))  = 1.0_ki/2.0_ki*known(ig+31) + 1.0_ki/2.0_ki*known(ig+61)
      
      cg(is3(ig,4))  = 1.0_ki/2.0_ki*known(ig+31) + 1.0_ki/2.0_ki*known(ig+65)
   enddo

   !  four q 
   cg(69) = xneval(69)-sub1(69)-sub2(69)-sub3(69)-sub4(69)

   if (verbosity .ge. 1) then
      do ig=0,69
         print*, ig, cg(ig)
      enddo
   end if
 end subroutine solve_system_rk4


subroutine solve_system_rk5(xneval,known)
   implicit none
   complex(ki), dimension(0:209),intent(inout) :: xneval,known
   integer :: ig, igs
   integer, dimension(0:5,10) :: is2
   integer, dimension(0:3,10) :: is3

  ! one q
   do ig=0,3 
      cg(1+ig)  =  - 1.0_ki/3.0_ki*xneval(ig+1) + 1.0_ki/9.0_ki*xneval(ig+5) + 16.0_ki/&
           & 9.0_ki*xneval(ig+15) - 16.0_ki/15.0_ki*xneval(ig+35) + 1.0_ki/90.0_ki*&
           & xneval(ig+70)
      
      cg(5+ig)  =  - 1.0_ki/6.0_ki*xneval(ig+1) - 1.0_ki/6.0_ki*xneval(ig+5) + 8.0_ki/&
           & 3.0_ki*xneval(ig+15) + 8.0_ki/3.0_ki*xneval(ig+35)
      
      cg(15+ig) = 3.0_ki/2.0_ki*xneval(ig+1) - 7.0_ki/18.0_ki*xneval(ig+5) - 32.0_ki/9.&
           &0_ki*xneval(ig+15) - 1.0_ki/18.0_ki*xneval(ig+70)
      
      cg(35+ig) = 2.0_ki/3.0_ki*xneval(ig+1) + 2.0_ki/3.0_ki*xneval(ig+5) - 8.0_ki/3.0_&
           &ki*xneval(ig+15) - 8.0_ki/3.0_ki*xneval(ig+35)
      
      cg(70+ig) = - 2.0_ki/3.0_ki*xneval(ig+1) - 2.0_ki/9.0_ki*xneval(ig+5) + 16.0_ki/&
           & 9.0_ki*xneval(ig+15) + 16.0_ki/15.0_ki*xneval(ig+35) + 2.0_ki/45.0_ki*&
           & xneval(ig+70)
   enddo
  
   ! two q
   do igs=9,97
      known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-sub4(igs)-sub5(igs)
   enddo

   !1-2:   9,19,51,22,39,42,74,77,86,89
   !1-3:  10,20,52,23,40,43,75,78,87,90
   !1-4:  11,21,53,24,41,44,76,79,88,91
   !2-3:  12,25,54,27,45,47,80,82,92,94
   !2-4:  13,26,55,28,46,48,81,83,93,95
   !3-4:  14,29,56,30,49,50,84,85,96,97

   is2(0,:) =  (/ 9 ,19,22,39,42,51,74,77,86,89 /)
   is2(1,:) =  (/ 10,20,23,40,43,52,75,78,87,90 /)
   is2(2,:) =  (/ 11,21,24,41,44,53,76,79,88,91 /)
   is2(3,:) =  (/ 12,25,27,45,47,54,80,82,92,94 /)
   is2(4,:) =  (/ 13,26,28,46,48,55,81,83,93,95 /)
   is2(5,:) =  (/ 14,29,30,49,50,56,84,85,96,97 /)

   do ig=0,5
      cg(is2(ig,1))  = - 1.0_ki/36.0_ki*known(ig+9) + 5.0_ki/36.0_ki*known(ig+19) - 4.0_k&
           &i/9.0_ki*known(ig+25) - 4.0_ki/3.0_ki*known(ig+39) - 32.0_ki/9.0_ki*&
           & known(ig+45) - 1.0_ki/36.0_ki*known(ig+51) + 7.0_ki/12.0_ki*known(ig+74)&
           &  + 4.0_ki/9.0_ki*known(ig+80) + 4.0_ki/9.0_ki*known(ig+86) - 8.0_ki/9.0&
           &_ki*known(ig+92)
   
      cg(is2(ig,2)) = - 11.0_ki/12.0_ki*known(ig+9) - 1.0_ki/12.0_ki*known(ig+19) + 8.0_&
           &ki/3.0_ki*known(ig+25) + 8.0_ki/3.0_ki*known(ig+39) - 1.0_ki/12.0_ki*&
           & known(ig+51) - 11.0_ki/12.0_ki*known(ig+74) + 4.0_ki/3.0_ki*known(ig+80)&
           &  + 4.0_ki/3.0_ki*known(ig+92)
   
      cg(is2(ig,3)) =  - 35.0_ki/36.0_ki*known(ig+9) - 11.0_ki/36.0_ki*known(ig+19) + 16.&
        &0_ki/9.0_ki*known(ig+25) + 32.0_ki/9.0_ki*known(ig+45) - 5.0_ki/36.0_ki&
        & *known(ig+51) + 1.0_ki/12.0_ki*known(ig+74) + 20.0_ki/9.0_ki*known(ig+80)&
        &  + 8.0_ki/9.0_ki*known(ig+86) - 4.0_ki/9.0_ki*known(ig+92)
   
      cg(is2(ig,4)) = 1.0_ki/9.0_ki*known(ig+9) + 1.0_ki/9.0_ki*known(ig+19) + 4.0_ki/9.0&
        &_ki*known(ig+25) + 4.0_ki/3.0_ki*known(ig+39) + 32.0_ki/9.0_ki*known(ig+&
        & 45) + 1.0_ki/9.0_ki*known(ig+51) - 1.0_ki/3.0_ki*known(ig+74) - 4.0_ki/&
        & 9.0_ki*known(ig+80) - 16.0_ki/9.0_ki*known(ig+86) - 4.0_ki/9.0_ki*&
        & known(ig+92)
   
      cg(is2(ig,5)) = 1.0_ki/6.0_ki*known(ig+9) - 1.0_ki/2.0_ki*known(ig+19) + 1.0_ki/6.0&
           &_ki*known(ig+51) - 1.0_ki/2.0_ki*known(ig+74) + 4.0_ki/3.0_ki*known(ig+86&
           & ) + 4.0_ki/3.0_ki*known(ig+92)
   
      cg(is2(ig,6)) = 1.0_ki/4.0_ki*known(ig+9) + 1.0_ki/4.0_ki*known(ig+19) + 1.0_ki/4.0&
           &_ki*known(ig+51) + 1.0_ki/4.0_ki*known(ig+74)
   
      cg(is2(ig,7)) = 2.0_ki/3.0_ki*known(ig+9) - 8.0_ki/3.0_ki*known(ig+25) - 8.0_ki/3.0&
           &_ki*known(ig+39) + 2.0_ki/3.0_ki*known(ig+74)
   
      cg(is2(ig,8)) = 2.0_ki/3.0_ki*known(ig+9) + 2.0_ki/3.0_ki*known(ig+19) - 8.0_ki/3.0&
           &_ki*known(ig+80) - 8.0_ki/3.0_ki*known(ig+86)
      
      cg(is2(ig,9)) = 5.0_ki/9.0_ki*known(ig+9) - 1.0_ki/9.0_ki*known(ig+19) - 16.0_ki/9.&
           &0_ki*known(ig+25) - 32.0_ki/9.0_ki*known(ig+45) - 1.0_ki/9.0_ki*known(ig+&
           & 51) - 1.0_ki/3.0_ki*known(ig+74) + 4.0_ki/9.0_ki*known(ig+80) + 16.0_ki&
           &/9.0_ki*known(ig+86) + 4.0_ki/9.0_ki*known(ig+92)
   
      cg(is2(ig,10)) = 1.0_ki/2.0_ki*known(ig+9) - 1.0_ki/6.0_ki*known(ig+19) - 1.0_ki/6.0&
           &_ki*known(ig+51) + 1.0_ki/2.0_ki*known(ig+74) - 4.0_ki/3.0_ki*known(ig+80&
           & ) - 4.0_ki/3.0_ki*known(ig+92)
   enddo
 
! three q ------------------------------------------------------------
   
   !1-2-3: 31,57,60,62, 98,101,103,110,112,116
   !1-2-4: 32,58,61,63, 99,102,104,111,113,117
   !1-3-4: 33,59,64,65,100,105,106,114,115,118
   !2-3-4: 34,66,67,68,107,108,109,119,120,121

   is3(0,:) =  (/ 31,57,60,62, 98,101,103,110,112,116 /)
   is3(1,:) =  (/ 32,58,61,63, 99,102,104,111,113,117 /)
   is3(2,:) =  (/ 33,59,64,65,100,105,106,114,115,118 /)
   is3(3,:) =  (/ 34,66,67,68,107,108,109,119,120,121 /)

   do igs = 31,34
      known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-sub4(igs)-sub5(igs)
   enddo
   do igs = 57,121
      known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-sub4(igs)-sub5(igs)
   enddo
   
   do ig=0,3
      cg(is3(ig,1))  =  - 11.0_ki/4.0_ki*known(ig+31) + 1.0_ki/3.0_ki*known(ig+57) + 1.0&
        &_ki/3.0_ki*known(ig+61) + 1.0_ki/3.0_ki*known(ig+65) + 1.0_ki/4.0_ki*&
        & known(ig+98) + 1.0_ki/4.0_ki*known(ig+102) + 1.0_ki/4.0_ki*known(ig+106)&
        &  + 8.0_ki/3.0_ki*known(ig+110) + 8.0_ki/3.0_ki*known(ig+114) + 8.0_ki/&
        & 3.0_ki*known(ig+118)

      cg(is3(ig,2))  =  - 1.0_ki/4.0_ki*known(ig+61) - 1.0_ki/4.0_ki*known(ig+65) - 1.0_&
           &ki/4.0_ki*known(ig+98) - 1.0_ki/4.0_ki*known(ig+102)
   
      cg(is3(ig,3))  =  - 1.0_ki/4.0_ki*known(ig+57) - 1.0_ki/4.0_ki*known(ig+65) - 1.0_&
           &ki/4.0_ki*known(ig+98) - 1.0_ki/4.0_ki*known(ig+106)
   
      cg(is3(ig,4))  = - 1.0_ki/4.0_ki*known(ig+57) - 1.0_ki/4.0_ki*known(ig+61) - 1.0_&
        &ki/4.0_ki*known(ig+102) - 1.0_ki/4.0_ki*known(ig+106)
      
      cg(is3(ig,5))  = known(ig+31) - 1.0_ki/3.0_ki*known(ig+57) - 8.0_ki/3.0_ki*known(ig+110)
      
      cg(is3(ig,6))  = known(ig+31) - 1.0_ki/3.0_ki*known(ig+61) - 8.0_ki/3.0_ki*known(ig+114)
   
      cg(is3(ig,7))  = known(ig+31) - 1.0_ki/3.0_ki*known(ig+65) - 8.0_ki/3.0_ki*known(ig+118)
   
      cg(is3(ig,8))  = 1.0_ki/4.0_ki*known(ig+31) + 1.0_ki/4.0_ki*known(ig+57) + 1.0_ki/&
           & 4.0_ki*known(ig+61) + 1.0_ki/4.0_ki*known(ig+98)
   
      cg(is3(ig,9))  = 1.0_ki/4.0_ki*known(ig+31) + 1.0_ki/4.0_ki*known(ig+57) + 1.0_ki/&
           & 4.0_ki*known(ig+65) + 1.0_ki/4.0_ki*known(ig+102)
   
      cg(is3(ig,10))  = 1.0_ki/4.0_ki*known(ig+31) + 1.0_ki/4.0_ki*known(ig+61) + 1.0_ki/&
           & 4.0_ki*known(ig+65) + 1.0_ki/4.0_ki*known(ig+106)
   enddo

! four q ------------------------------------------------------------- !   
   do igs = 122,125
      known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-sub4(igs)-sub5(igs)
   enddo
   known(69) = xneval(69)-sub1(69)-sub2(69)-sub3(69)-sub4(69)-sub5(69)

!!$   cg(69)  = - known(69) - 1.0_ki/2.0_ki*known(122) - 1.0_ki/2.0_ki*&
!!$        & known(123) - 1.0_ki/2.0_ki*known(124) - 1.0_ki/2.0_ki*known(125)   
!!$   cg(122)  = 1.0_ki/2.0_ki*known(69) + 1.0_ki/2.0_ki*known(122)
!!$   cg(123)  = 1.0_ki/2.0_ki*known(69) + 1.0_ki/2.0_ki*known(123)
!!$   cg(124)  = 1.0_ki/2.0_ki*known(69) + 1.0_ki/2.0_ki*known(124)
!!$   cg(125)  = 1.0_ki/2.0_ki*known(69) + 1.0_ki/2.0_ki*known(125)
  
   cg(69) = (-two*known(69) - known(122) - known (123) - known (124) - known (125))/two
   cg(122) =   (known(69) + known(122))/two
   cg(123) =   (known(69) + known(123))/two
   cg(124) =   (known(69) + known(124))/two
   cg(125) =   (known(69) + known(125))/two

   if (verbosity .ge. 1) then
      do ig=0,125
         print*, ig, cg(ig)
      enddo
   end if
!   stop
 end subroutine solve_system_rk5

subroutine solve_system_rk6(xneval,known)
   implicit none
   complex(ki), dimension(0:209),intent(inout) :: xneval,known
   integer :: ig, igs 
   integer, dimension(0:5,15) :: is2
   integer, dimension(0:3,20) :: is3
   ! one q
   do ig=0,3              
      cg(1+ig)  =  - 2.0_ki/9.0_ki*xneval(ig+1) + 2.0_ki/9.0_ki*xneval(ig+5) + 64.0_ki/&
           & 45.0_ki*xneval(ig+15) - 64.0_ki/45.0_ki*xneval(ig+35) + 1.0_ki/180.0_ki*&
           & xneval(ig+70) - 1.0_ki/180.0_ki*xneval(ig+126)
      
      cg(5+ig)  =  - 2.0_ki/9.0_ki*xneval(ig+1) - 2.0_ki/9.0_ki*xneval(ig+5) + 128.0_ki&
           &/45.0_ki*xneval(ig+15) + 128.0_ki/45.0_ki*xneval(ig+35) + 1.0_ki/360.0_ki&
           & *xneval(ig+70) + 1.0_ki/360.0_ki*xneval(ig+126)

      cg(15+ig) = 17.0_ki/18.0_ki*xneval(ig+1) - 17.0_ki/18.0_ki*xneval(ig+5) - 16.0_ki&
           &/9.0_ki*xneval(ig+15) + 16.0_ki/9.0_ki*xneval(ig+35) - 1.0_ki/36.0_ki*&
           & xneval(ig+70) + 1.0_ki/36.0_ki*xneval(ig+126)

      cg(35+ig) = 17.0_ki/18.0_ki*xneval(ig+1) + 17.0_ki/18.0_ki*xneval(ig+5) - 32.0_ki&
           &/9.0_ki*xneval(ig+15) - 32.0_ki/9.0_ki*xneval(ig+35) - 1.0_ki/72.0_ki*&
           & xneval(ig+70) - 1.0_ki/72.0_ki*xneval(ig+126)

      cg(70+ig) =  - 2.0_ki/9.0_ki*xneval(ig+1) + 2.0_ki/9.0_ki*xneval(ig+5) + 16.0_ki/&
           & 45.0_ki*xneval(ig+15) - 16.0_ki/45.0_ki*xneval(ig+35) + 1.0_ki/45.0_ki*&
           & xneval(ig+70) - 1.0_ki/45.0_ki*xneval(ig+126)

      cg(126+ig)=  - 2.0_ki/9.0_ki*xneval(ig+1) - 2.0_ki/9.0_ki*xneval(ig+5) + 32.0_ki/&
           & 45.0_ki*xneval(ig+15) + 32.0_ki/45.0_ki*xneval(ig+35) + 1.0_ki/90.0_ki*&
           & xneval(ig+70) + 1.0_ki/90.0_ki*xneval(ig+126)
   enddo

   ! two q
   do igs=9,159
      known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-&
           sub4(igs)-sub5(igs)-sub6(igs)
   enddo
   !1-2:   9,19,51,22,39,42,74,77,86,89,130,133,142,145,166
   !1-3:  10,20,52,23,40,43,75,78,87,90,131,134,143,146,167
   !1-4:  11,21,53,24,41,44,76,79,88,91,132,135,144,147,168
   !2-3:  12,25,54,27,45,47,80,82,92,94,136,138,148,150,169
   !2-4:  13,26,55,28,46,48,81,83,93,95,137,139,149,151,170
   !3-4:  14,29,56,30,49,50,84,85,96,97,140,141,152,153,171

   is2(0,:) =  (/ 9 ,19,22,39,42,51,74,77,86,89,130,133,142,145,166 /)
   is2(1,:) =  (/ 10,20,23,40,43,52,75,78,87,90,131,134,143,146,167 /)
   is2(2,:) =  (/ 11,21,24,41,44,53,76,79,88,91,132,135,144,147,168 /)
   is2(3,:) =  (/ 12,25,27,45,47,54,80,82,92,94,136,138,148,150,169 /)
   is2(4,:) =  (/ 13,26,28,46,48,55,81,83,93,95,137,139,149,151,170 /)
   is2(5,:) =  (/ 14,29,30,49,50,56,84,85,96,97,140,141,152,153,171 /)
   
   do ig=0,5
      cg(is2(ig,11))  =  - 1.0_ki/4.0_ki*known(ig+9) - 11.0_ki/36.0_ki*known(ig+19) + 1.0_k&
           &i/3.0_ki*known(ig+25) + 2.0_ki/3.0_ki*known(ig+39) - 24.0_ki/5.0_ki*&
           & known(ig+45) - 1.0_ki/36.0_ki*known(ig+51) - 1.0_ki/12.0_ki*known(ig+74)&
           &  + 5.0_ki/9.0_ki*known(ig+80) + 4.0_ki/3.0_ki*known(ig+86) + 1.0_ki/3.0&
           &_ki*known(ig+92) + known(ig+130) + 2.0_ki/9.0_ki*known(ig+136) - 136.0_ki/&
           & 45.0_ki*known(ig+142) + 1.0_ki/45.0_ki*known(ig+148) + 1.0_ki/45.0_ki*&
           & known(ig+154)

      cg(is2(ig,2)) = - 5.0_ki/12.0_ki*known(ig+9) + 13.0_ki/12.0_ki*known(ig+19) + 2.0_&
           &ki*known(ig+25) + 10.0_ki/3.0_ki*known(ig+39) + 16.0_ki/3.0_ki*known(ig+&
           & 45) + 5.0_ki/12.0_ki*known(ig+51) - 13.0_ki/12.0_ki*known(ig+74) - 8.0_&
           &ki/3.0_ki*known(ig+86) + 8.0_ki/3.0_ki*known(ig+92) - 10.0_ki/3.0_ki*&
           & known(ig+130) - 2.0_ki*known(ig+136) - 16.0_ki/3.0_ki*known(ig+142)
      
      cg(is2(ig,3)) = - 5.0_ki/12.0_ki*known(ig+9) + 11.0_ki/12.0_ki*known(ig+19) + 2.0_&
           &ki*known(ig+25) + 10.0_ki/3.0_ki*known(ig+39) + 32.0_ki/3.0_ki*known(ig+&
           & 45) + 5.0_ki/12.0_ki*known(ig+51) - 11.0_ki/12.0_ki*known(ig+74) - 8.0_&
           &ki/3.0_ki*known(ig+86) + 8.0_ki/3.0_ki*known(ig+92) - 10.0_ki/3.0_ki*&
           & known(ig+130) - 2.0_ki*known(ig+136) - 32.0_ki/3.0_ki*known(ig+142)
      
      cg(is2(ig,4)) = 5.0_ki/12.0_ki*known(ig+9) - 13.0_ki/36.0_ki*known(ig+19) - 7.0_ki/&
           & 3.0_ki*known(ig+25) - 11.0_ki/3.0_ki*known(ig+39) - 5.0_ki/36.0_ki*&
           & known(ig+51) + 13.0_ki/12.0_ki*known(ig+74) + 13.0_ki/9.0_ki*known(ig+80)&
           &  - 29.0_ki/9.0_ki*known(ig+92) + 11.0_ki/9.0_ki*known(ig+130) + 7.0_ki/&
           & 9.0_ki*known(ig+136) + 80.0_ki/9.0_ki*known(ig+142) - 1.0_ki/9.0_ki*&
           & known(ig+154)
      
      cg(is2(ig,5)) = 5.0_ki/12.0_ki*known(ig+9) + 11.0_ki/36.0_ki*known(ig+19) + 2.0_ki/&
           & 3.0_ki*known(ig+25) - 7.0_ki/3.0_ki*known(ig+39) + 8.0_ki/3.0_ki*&
           & known(ig+45) - 5.0_ki/36.0_ki*known(ig+51) + 5.0_ki/12.0_ki*known(ig+74)&
           &  - 14.0_ki/9.0_ki*known(ig+80) - 4.0_ki/3.0_ki*known(ig+86) - 10.0_ki/&
           & 9.0_ki*known(ig+92) - 8.0_ki/9.0_ki*known(ig+130) + 7.0_ki/9.0_ki*&
           & known(ig+136) + 56.0_ki/9.0_ki*known(ig+142) - 1.0_ki/9.0_ki*known(ig+148)
      
      cg(is2(ig,6)) = - 5.0_ki/12.0_ki*known(ig+9) - 7.0_ki/4.0_ki*known(ig+19) - 8.0_ki&
           &/3.0_ki*known(ig+39) - 32.0_ki/3.0_ki*known(ig+45) - 5.0_ki/12.0_ki*&
           & known(ig+51) + 11.0_ki/12.0_ki*known(ig+74) + 8.0_ki/3.0_ki*known(ig+80)&
           &  + 16.0_ki/3.0_ki*known(ig+86) - 8.0_ki/3.0_ki*known(ig+92) + 16.0_ki/&
           & 3.0_ki*known(ig+130) + 8.0_ki/3.0_ki*known(ig+136) + 32.0_ki/3.0_ki*&
           & known(ig+142)
      
      cg(is2(ig,7)) = 1.0_ki/3.0_ki*known(ig+9) - 1.0_ki/3.0_ki*known(ig+19) - 4.0_ki/3.0&
           &_ki*known(ig+25) - 4.0_ki/3.0_ki*known(ig+39) - 1.0_ki/3.0_ki*known(ig+51&
           & ) + 1.0_ki/3.0_ki*known(ig+74) + 4.0_ki/3.0_ki*known(ig+130) + 4.0_ki/&
           & 3.0_ki*known(ig+136)
      
      cg(is2(ig,8)) = 1.0_ki/3.0_ki*known(ig+9) - known(ig+19) - 4.0_ki/3.0_ki*known(ig+25)&
           &  - 4.0_ki*known(ig+39) - 32.0_ki/3.0_ki*known(ig+45) - 1.0_ki/3.0_ki*&
           & known(ig+51) + known(ig+74) + 8.0_ki/3.0_ki*known(ig+86) - 8.0_ki/3.0_ki*&
           & known(ig+92) + 4.0_ki*known(ig+130) + 4.0_ki/3.0_ki*known(ig+136) + 32.0_k&
           &i/3.0_ki*known(ig+142)
      
      cg(is2(ig,9)) = 1.0_ki/3.0_ki*known(ig+9) + 1.0_ki/3.0_ki*known(ig+19) - 2.0_ki/3.0&
           &_ki*known(ig+25) + 2.0_ki/3.0_ki*known(ig+39) - 1.0_ki/3.0_ki*known(ig+51&
           & ) - 1.0_ki/3.0_ki*known(ig+74) - 2.0_ki/3.0_ki*known(ig+130) + 2.0_ki/&
           & 3.0_ki*known(ig+136)
      
      cg(is2(ig,10)) = 1.0_ki/3.0_ki*known(ig+9) - known(ig+19) - 2.0_ki/3.0_ki*known(ig+25)&
           &  - 2.0_ki*known(ig+39) - 16.0_ki/3.0_ki*known(ig+45) - 1.0_ki/3.0_ki*&
           & known(ig+51) + known(ig+74) + 8.0_ki/3.0_ki*known(ig+86) - 8.0_ki/3.0_ki*&
           & known(ig+92) + 2.0_ki*known(ig+130) + 2.0_ki/3.0_ki*known(ig+136) + 16.0_k&
           &i/3.0_ki*known(ig+142)
      
      cg(is2(ig,11))=  - 1.0_ki/3.0_ki*known(ig+9) + 1.0_ki/9.0_ki*known(ig+19) + 4.0_ki/&
           & 3.0_ki*known(ig+25) + 4.0_ki/3.0_ki*known(ig+39) + 1.0_ki/9.0_ki*&
           & known(ig+51) - 1.0_ki/3.0_ki*known(ig+74) - 4.0_ki/9.0_ki*known(ig+80) + &
           & 4.0_ki/9.0_ki*known(ig+92) - 4.0_ki/9.0_ki*known(ig+130) - 4.0_ki/9.0_k&
           &i*known(ig+136) - 64.0_ki/45.0_ki*known(ig+142) + 4.0_ki/45.0_ki*&
           & known(ig+154)

      cg(is2(ig,12))= - 1.0_ki/3.0_ki*known(ig+9) + 1.0_ki/9.0_ki*known(ig+19) + 4.0_ki/&
           & 3.0_ki*known(ig+39) + 32.0_ki/15.0_ki*known(ig+45) + 1.0_ki/9.0_ki*&
           & known(ig+51) - 1.0_ki/3.0_ki*known(ig+74) + 8.0_ki/9.0_ki*known(ig+80) + &
           & 8.0_ki/9.0_ki*known(ig+92) - 8.0_ki/9.0_ki*known(ig+130) - 4.0_ki/9.0_k&
           &i*known(ig+136) - 32.0_ki/9.0_ki*known(ig+142) + 4.0_ki/45.0_ki*known(ig+148)

      cg(is2(ig,13))= 1.0_ki/3.0_ki*known(ig+9) + 1.0_ki/3.0_ki*known(ig+19) - 4.0_ki/3.0&
           &_ki*known(ig+25) - 4.0_ki/3.0_ki*known(ig+39) + 1.0_ki/3.0_ki*known(ig+51&
           & ) + 1.0_ki/3.0_ki*known(ig+74) - 4.0_ki/3.0_ki*known(ig+130) - 4.0_ki/&
           & 3.0_ki*known(ig+136)

      cg(is2(ig,14))= 1.0_ki/3.0_ki*known(ig+9) + 5.0_ki/3.0_ki*known(ig+19) + 4.0_ki/3.0&
           &_ki*known(ig+25) + 4.0_ki*known(ig+39) + 32.0_ki/3.0_ki*known(ig+45) + 1.0&
           &_ki/3.0_ki*known(ig+51) - known(ig+74) - 8.0_ki/3.0_ki*known(ig+80) - 16.0&
           &_ki/3.0_ki*known(ig+86) + 8.0_ki/3.0_ki*known(ig+92) - 4.0_ki*known(ig+&
           & 130) - 4.0_ki/3.0_ki*known(ig+136) - 32.0_ki/3.0_ki*known(ig+142)

      cg(is2(ig,15))= 1.0_ki/3.0_ki*known(ig+9) - 1.0_ki/9.0_ki*known(ig+19) + 8.0_ki/3.0&
           &_ki*known(ig+39) + 1.0_ki/3.0_ki*known(ig+51) - known(ig+74) - 8.0_ki/9.0_&
           &ki*known(ig+80) + 8.0_ki/3.0_ki*known(ig+92) - 8.0_ki/9.0_ki*known(ig+136&
           & ) - 64.0_ki/9.0_ki*known(ig+142)
   enddo
 
!--------------------------------------------------------------------
      do igs=31,199
         known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-&
              sub4(igs)-sub5(igs)-sub6(igs)
      enddo
! three q
!1-2-3: 31,57,60,62, 98,101,103,110,112,116,154,157,159,172,174,178,180,184,186,200
!1-2-4: 32,58,61,63, 99,102,104,111,113,117,155,158,160,173,175,179,181,185,187,201
!1-3-4: 33,59,64,65,100,105,106,114,115,118,156,161,162,176,177,182,183,188,189,202
!2-3-4: 34,66,67,68,107,108,109,119,120,121,163,164,165,190,191,192,193,194,195,203      
   is3(0,:) =  (/ 31,57,60,62, 98,101,103,110,112,116,154,157,159,172,174,178,180,184,186,200 /)
   is3(1,:) =  (/ 32,58,61,63, 99,102,104,111,113,117,155,158,160,173,175,179,181,185,187,201 /)
   is3(2,:) =  (/ 33,59,64,65,100,105,106,114,115,118,156,161,162,176,177,182,183,188,189,202 /)
   is3(3,:) =  (/ 34,66,67,68,107,108,109,119,120,121,163,164,165,190,191,192,193,194,195,203 /)
   
   do ig=0,3
      cg(is3(ig,1))  = - 2687.0_ki/4680.0_ki*known(ig+31) - 293.0_ki/4680.0_ki*&
           & known(ig+57) + 37.0_ki/180.0_ki*known(ig+61) + 2269.0_ki/4680.0_ki*&
           & known(ig+65) + 23.0_ki/234.0_ki*known(ig+98) - 337.0_ki/936.0_ki*&
           & known(ig+102) - 547.0_ki/2340.0_ki*known(ig+106) - 121.0_ki/195.0_ki*&
           & known(ig+110) + 892.0_ki/585.0_ki*known(ig+114) + 4.0_ki/3.0_ki*known(ig+&
           & 118) - 1711.0_ki/1755.0_ki*known(ig+160) - 734.0_ki/1755.0_ki*&
           & known(ig+164) + 17824.0_ki/8775.0_ki*known(ig+168) - 83.0_ki/1755.0_ki&
           & *known(ig+172) - 1858.0_ki/2925.0_ki*known(ig+176) + 8588.0_ki/1755.0_k&
           &i*known(ig+180) + 692.0_ki/351.0_ki*known(ig+184) + 36544.0_ki/8775.0_k&
           &i*known(ig+188) - 224.0_ki/1755.0_ki*known(ig+192) + 1264.0_ki/975.0_ki&
           & *known(ig+196)

      cg(is3(ig,2))  = - 47.0_ki/39.0_ki*known(ig+31) + 1.0_ki/39.0_ki*known(ig+57) + &
           & 1.0_ki/12.0_ki*known(ig+61) - 95.0_ki/156.0_ki*known(ig+65) - 9.0_ki/&
           & 52.0_ki*known(ig+98) + 25.0_ki/156.0_ki*known(ig+102) + 4.0_ki/39.0_ki&
           & *known(ig+106) + 12.0_ki/13.0_ki*known(ig+110) + 80.0_ki/39.0_ki*&
           & known(ig+114) + 76.0_ki/117.0_ki*known(ig+160) + 200.0_ki/117.0_ki*&
           & known(ig+164) + 1856.0_ki/585.0_ki*known(ig+168) + 80.0_ki/117.0_ki*&
           & known(ig+172) + 88.0_ki/195.0_ki*known(ig+176) - 128.0_ki/117.0_ki*&
           & known(ig+180) - 448.0_ki/117.0_ki*known(ig+184) - 8128.0_ki/585.0_ki*&
           & known(ig+188) - 160.0_ki/117.0_ki*known(ig+192) - 48.0_ki/65.0_ki*&
           & known(ig+196)

      cg(is3(ig,3))  = - 269.0_ki/468.0_ki*known(ig+31) + 29.0_ki/234.0_ki*known(ig+57&
           & ) - 1.0_ki/4.0_ki*known(ig+61) - 305.0_ki/468.0_ki*known(ig+65) - 43.0_&
           &ki/234.0_ki*known(ig+98) + 37.0_ki/117.0_ki*known(ig+102) + 115.0_ki/&
           & 468.0_ki*known(ig+106) + 32.0_ki/13.0_ki*known(ig+110) - 140.0_ki/117.0&
           &_ki*known(ig+114) + 868.0_ki/351.0_ki*known(ig+160) + 560.0_ki/351.0_ki&
           & *known(ig+164) - 4912.0_ki/1755.0_ki*known(ig+168) - 10.0_ki/351.0_ki*&
           & known(ig+172) + 574.0_ki/585.0_ki*known(ig+176) - 3104.0_ki/351.0_ki*&
           & known(ig+180) - 1504.0_ki/351.0_ki*known(ig+184) + 5696.0_ki/1755.0_ki&
           & *known(ig+188) + 176.0_ki/351.0_ki*known(ig+192) - 128.0_ki/65.0_ki*&
           & known(ig+196)

      cg(is3(ig,4))  = - 101.0_ki/195.0_ki*known(ig+31) - 7.0_ki/260.0_ki*known(ig+57)&
           &  - 4.0_ki/45.0_ki*known(ig+61) + 241.0_ki/260.0_ki*known(ig+65) + 29.0_&
           &ki/156.0_ki*known(ig+98) - 5.0_ki/78.0_ki*known(ig+102) - 253.0_ki/780.&
           &0_ki*known(ig+106) + 58.0_ki/195.0_ki*known(ig+110) - 116.0_ki/585.0_ki&
           & *known(ig+114) + 8.0_ki/3.0_ki*known(ig+118) - 74.0_ki/585.0_ki*known(ig+&
           & 160) - 556.0_ki/585.0_ki*known(ig+164) - 5584.0_ki/2925.0_ki*known(ig+&
           & 168) - 472.0_ki/585.0_ki*known(ig+172) - 722.0_ki/975.0_ki*known(ig+&
           & 176) + 1192.0_ki/585.0_ki*known(ig+180) + 136.0_ki/39.0_ki*known(ig+&
           & 184) + 1344.0_ki/325.0_ki*known(ig+188) + 944.0_ki/585.0_ki*known(ig+&
           & 192) + 1328.0_ki/975.0_ki*known(ig+196)

      cg(is3(ig,5))  = 41.0_ki/65.0_ki*known(ig+31) - 43.0_ki/195.0_ki*known(ig+57) + 1.&
           &0_ki/10.0_ki*known(ig+61) - 19.0_ki/130.0_ki*known(ig+65) + 1.0_ki/26.0&
           &_ki*known(ig+98) + 1.0_ki/26.0_ki*known(ig+102) + 166.0_ki/195.0_ki*&
           & known(ig+106) - 404.0_ki/195.0_ki*known(ig+110) + 32.0_ki/65.0_ki*&
           & known(ig+114) - 656.0_ki/195.0_ki*known(ig+160) - 284.0_ki/195.0_ki*&
           & known(ig+164) - 256.0_ki/975.0_ki*known(ig+168) + 32.0_ki/195.0_ki*&
           & known(ig+172) - 48.0_ki/325.0_ki*known(ig+176) + 448.0_ki/195.0_ki*&
           & known(ig+180) + 64.0_ki/39.0_ki*known(ig+184) - 256.0_ki/975.0_ki*&
           & known(ig+188) - 64.0_ki/195.0_ki*known(ig+192) + 192.0_ki/325.0_ki*&
           & known(ig+196)

      cg(is3(ig,6))  = - 1771.0_ki/4680.0_ki*known(ig+31) + 431.0_ki/4680.0_ki*&
           & known(ig+57) - 143.0_ki/360.0_ki*known(ig+61) - 989.0_ki/2340.0_ki*&
           & known(ig+65) - 173.0_ki/936.0_ki*known(ig+98) + 103.0_ki/234.0_ki*&
           & known(ig+102) - 191.0_ki/2340.0_ki*known(ig+106) + 144.0_ki/65.0_ki*&
           & known(ig+110) - 994.0_ki/585.0_ki*known(ig+114) + 6142.0_ki/1755.0_ki*&
           & known(ig+160) + 4808.0_ki/1755.0_ki*known(ig+164) - 12328.0_ki/8775.0_k&
           &i*known(ig+168) - 19.0_ki/1755.0_ki*known(ig+172) + 2401.0_ki/2925.0_ki&
           & *known(ig+176) - 14696.0_ki/1755.0_ki*known(ig+180) - 2768.0_ki/351.0_k&
           &i*known(ig+184) - 45088.0_ki/8775.0_ki*known(ig+188) + 428.0_ki/1755.0_&
           &ki*known(ig+192) - 576.0_ki/325.0_ki*known(ig+196)

      cg(is3(ig,7))  = 187.0_ki/390.0_ki*known(ig+31) + 41.0_ki/780.0_ki*known(ig+57)&
           &  - 1.0_ki/120.0_ki*known(ig+61) - 121.0_ki/1560.0_ki*known(ig+65) + 19.&
           &0_ki/104.0_ki*known(ig+98) + 5.0_ki/312.0_ki*known(ig+102) - 187.0_ki/&
           & 780.0_ki*known(ig+106) + 19.0_ki/65.0_ki*known(ig+110) - 38.0_ki/195.0_&
           &ki*known(ig+114) - 4.0_ki/3.0_ki*known(ig+118) - 1.0_ki/585.0_ki*&
           & known(ig+160) - 914.0_ki/585.0_ki*known(ig+164) - 1256.0_ki/2925.0_ki*&
           & known(ig+168) - 38.0_ki/585.0_ki*known(ig+172) - 73.0_ki/975.0_ki*&
           & known(ig+176) + 1028.0_ki/585.0_ki*known(ig+180) + 548.0_ki/117.0_ki*&
           & known(ig+184) + 3424.0_ki/2925.0_ki*known(ig+188) + 76.0_ki/585.0_ki*&
           & known(ig+192) + 32.0_ki/975.0_ki*known(ig+196)

      cg(is3(ig,8))  = 41.0_ki/260.0_ki*known(ig+31) + 29.0_ki/260.0_ki*known(ig+57) + &
           & 3.0_ki/20.0_ki*known(ig+61) - 21.0_ki/130.0_ki*known(ig+65) + 7.0_ki/&
           & 52.0_ki*known(ig+98) - 3.0_ki/26.0_ki*known(ig+102) + 3.0_ki/65.0_ki*&
           & known(ig+106) - 12.0_ki/65.0_ki*known(ig+110) + 8.0_ki/65.0_ki*known(ig+&
           & 114) - 164.0_ki/195.0_ki*known(ig+160) - 136.0_ki/195.0_ki*known(ig+&
           & 164) - 64.0_ki/975.0_ki*known(ig+168) + 8.0_ki/195.0_ki*known(ig+172)&
           &  - 12.0_ki/325.0_ki*known(ig+176) + 112.0_ki/195.0_ki*known(ig+180) + &
           & 16.0_ki/39.0_ki*known(ig+184) - 64.0_ki/975.0_ki*known(ig+188) - 16.0_k&
           &i/195.0_ki*known(ig+192) + 48.0_ki/325.0_ki*known(ig+196)

      cg(is3(ig,9))  = 41.0_ki/260.0_ki*known(ig+31) + 29.0_ki/260.0_ki*known(ig+57) - &
           & 1.0_ki/10.0_ki*known(ig+61) + 23.0_ki/260.0_ki*known(ig+65) - 3.0_ki/&
           & 26.0_ki*known(ig+98) + 7.0_ki/52.0_ki*known(ig+102) + 3.0_ki/65.0_ki*&
           & known(ig+106) - 12.0_ki/65.0_ki*known(ig+110) + 8.0_ki/65.0_ki*known(ig+&
           & 114) - 164.0_ki/195.0_ki*known(ig+160) - 136.0_ki/195.0_ki*known(ig+&
           & 164) - 64.0_ki/975.0_ki*known(ig+168) + 8.0_ki/195.0_ki*known(ig+172)&
           &  - 12.0_ki/325.0_ki*known(ig+176) + 112.0_ki/195.0_ki*known(ig+180) + &
           & 16.0_ki/39.0_ki*known(ig+184) - 64.0_ki/975.0_ki*known(ig+188) - 16.0_k&
           &i/195.0_ki*known(ig+192) + 48.0_ki/325.0_ki*known(ig+196)
      
      cg(is3(ig,10))  = 41.0_ki/260.0_ki*known(ig+31) - 9.0_ki/65.0_ki*known(ig+57) + 3.0&
           &_ki/20.0_ki*known(ig+61) + 23.0_ki/260.0_ki*known(ig+65) - 3.0_ki/26.0_&
           &ki*known(ig+98) - 3.0_ki/26.0_ki*known(ig+102) + 77.0_ki/260.0_ki*&
           & known(ig+106) - 12.0_ki/65.0_ki*known(ig+110) + 8.0_ki/65.0_ki*known(ig+&
           & 114) - 164.0_ki/195.0_ki*known(ig+160) - 136.0_ki/195.0_ki*known(ig+&
           & 164) - 64.0_ki/975.0_ki*known(ig+168) + 8.0_ki/195.0_ki*known(ig+172)&
           &  - 12.0_ki/325.0_ki*known(ig+176) + 112.0_ki/195.0_ki*known(ig+180) + &
           & 16.0_ki/39.0_ki*known(ig+184) - 64.0_ki/975.0_ki*known(ig+188) - 16.0_k&
           &i/195.0_ki*known(ig+192) + 48.0_ki/325.0_ki*known(ig+196)

      cg(is3(ig,11))  = - 34.0_ki/195.0_ki*known(ig+31) + 14.0_ki/195.0_ki*known(ig+57)&
           &  - 2.0_ki/15.0_ki*known(ig+61) + 38.0_ki/195.0_ki*known(ig+65) - 2.0_ki&
           &/39.0_ki*known(ig+98) - 2.0_ki/39.0_ki*known(ig+102) - 16.0_ki/65.0_ki&
           & *known(ig+106) + 64.0_ki/65.0_ki*known(ig+110) - 128.0_ki/195.0_ki*&
           & known(ig+114) + 1064.0_ki/585.0_ki*known(ig+160) + 616.0_ki/585.0_ki*&
           & known(ig+164) + 1024.0_ki/2925.0_ki*known(ig+168) - 128.0_ki/585.0_ki*&
           & known(ig+172) + 64.0_ki/325.0_ki*known(ig+176) - 1792.0_ki/585.0_ki*&
           & known(ig+180) - 256.0_ki/117.0_ki*known(ig+184) + 1024.0_ki/2925.0_ki*&
           & known(ig+188) + 256.0_ki/585.0_ki*known(ig+192) - 256.0_ki/325.0_ki*&
           & known(ig+196)

      cg(is3(ig,12))  = 11.0_ki/468.0_ki*known(ig+31) - 55.0_ki/468.0_ki*known(ig+57) + &
           & 7.0_ki/36.0_ki*known(ig+61) + 73.0_ki/234.0_ki*known(ig+65) + 17.0_ki/&
           & 468.0_ki*known(ig+98) - 25.0_ki/117.0_ki*known(ig+102) + 7.0_ki/234.0_k&
           &i*known(ig+106) - 16.0_ki/13.0_ki*known(ig+110) + 148.0_ki/117.0_ki*&
           & known(ig+114) - 460.0_ki/351.0_ki*known(ig+160) - 176.0_ki/351.0_ki*&
           & known(ig+164) + 4432.0_ki/1755.0_ki*known(ig+168) + 70.0_ki/351.0_ki*&
           & known(ig+172) - 274.0_ki/585.0_ki*known(ig+176) + 1136.0_ki/351.0_ki*&
           & known(ig+180) + 544.0_ki/351.0_ki*known(ig+184) - 4928.0_ki/1755.0_ki*&
           & known(ig+188) - 296.0_ki/351.0_ki*known(ig+192) + 64.0_ki/65.0_ki*&
           & known(ig+196)

      cg(is3(ig,13))  = 149.0_ki/195.0_ki*known(ig+31) + 19.0_ki/130.0_ki*known(ig+57)&
           &  + 19.0_ki/180.0_ki*known(ig+61) - 129.0_ki/260.0_ki*known(ig+65) + 19.&
           &0_ki/156.0_ki*known(ig+98) + 19.0_ki/156.0_ki*known(ig+102) - 19.0_ki/&
           & 390.0_ki*known(ig+106) + 38.0_ki/195.0_ki*known(ig+110) - 76.0_ki/585.0&
           &_ki*known(ig+114) - 8.0_ki/3.0_ki*known(ig+118) - 694.0_ki/585.0_ki*&
           & known(ig+160) - 956.0_ki/585.0_ki*known(ig+164) - 10544.0_ki/2925.0_ki&
           & *known(ig+168) + 148.0_ki/585.0_ki*known(ig+172) + 298.0_ki/975.0_ki*&
           & known(ig+176) + 2072.0_ki/585.0_ki*known(ig+180) + 56.0_ki/13.0_ki*&
           & known(ig+184) + 7232.0_ki/975.0_ki*known(ig+188) - 296.0_ki/585.0_ki*&
           & known(ig+192) - 224.0_ki/325.0_ki*known(ig+196)

      cg(is3(ig,14))  = 107.0_ki/195.0_ki*known(ig+31) - 2.0_ki/195.0_ki*known(ig+57) + &
           & 1.0_ki/15.0_ki*known(ig+61) + 41.0_ki/195.0_ki*known(ig+65) + 4.0_ki/&
           & 39.0_ki*known(ig+98) + 4.0_ki/39.0_ki*known(ig+102) - 151.0_ki/195.0_ki&
           & *known(ig+106) - 176.0_ki/195.0_ki*known(ig+110) - 56.0_ki/195.0_ki*&
           & known(ig+114) + 368.0_ki/585.0_ki*known(ig+160) - 608.0_ki/585.0_ki*&
           & known(ig+164) - 5792.0_ki/2925.0_ki*known(ig+168) - 56.0_ki/585.0_ki*&
           & known(ig+172) - 436.0_ki/975.0_ki*known(ig+176) + 2336.0_ki/585.0_ki*&
           & known(ig+180) + 512.0_ki/117.0_ki*known(ig+184) + 12928.0_ki/2925.0_ki&
           & *known(ig+188) + 112.0_ki/585.0_ki*known(ig+192) + 704.0_ki/975.0_ki*&
           & known(ig+196)

      cg(is3(ig,15))  =  - 6.0_ki/65.0_ki*known(ig+31) - 9.0_ki/65.0_ki*known(ig+57) - 1.&
           &0_ki/10.0_ki*known(ig+61) - 21.0_ki/130.0_ki*known(ig+65) - 3.0_ki/26.0&
           &_ki*known(ig+98) - 3.0_ki/26.0_ki*known(ig+102) + 3.0_ki/65.0_ki*&
           & known(ig+106) - 12.0_ki/65.0_ki*known(ig+110) + 8.0_ki/65.0_ki*known(ig+&
           & 114) + 356.0_ki/195.0_ki*known(ig+160) + 128.0_ki/65.0_ki*known(ig+164&
           & ) + 672.0_ki/325.0_ki*known(ig+168) + 8.0_ki/195.0_ki*known(ig+172) + &
           & 484.0_ki/975.0_ki*known(ig+176) - 928.0_ki/195.0_ki*known(ig+180) - 64.&
           &0_ki/13.0_ki*known(ig+184) - 1408.0_ki/325.0_ki*known(ig+188) - 16.0_ki&
           &/195.0_ki*known(ig+192) - 896.0_ki/975.0_ki*known(ig+196)

      cg(is3(ig,16))  = 557.0_ki/390.0_ki*known(ig+31) + 23.0_ki/390.0_ki*known(ig+57)&
           &  + 1.0_ki/30.0_ki*known(ig+61) + 73.0_ki/195.0_ki*known(ig+65) + 19.0_k&
           &i/78.0_ki*known(ig+98) - 10.0_ki/39.0_ki*known(ig+102) + 7.0_ki/195.0_k&
           &i*known(ig+106) - 96.0_ki/65.0_ki*known(ig+110) - 328.0_ki/195.0_ki*&
           & known(ig+114) - 1856.0_ki/585.0_ki*known(ig+160) - 2224.0_ki/585.0_ki*&
           & known(ig+164) - 9856.0_ki/2925.0_ki*known(ig+168) - 328.0_ki/585.0_ki*&
           & known(ig+172) - 548.0_ki/975.0_ki*known(ig+176) + 4768.0_ki/585.0_ki*&
           & known(ig+180) + 1216.0_ki/117.0_ki*known(ig+184) + 40064.0_ki/2925.0_ki&
           & *known(ig+188) + 656.0_ki/585.0_ki*known(ig+192) + 384.0_ki/325.0_ki*&
           & known(ig+196)

      cg(is3(ig,17))  = 17.0_ki/390.0_ki*known(ig+31) - 7.0_ki/390.0_ki*known(ig+57) - 2.&
           &0_ki/15.0_ki*known(ig+61) - 19.0_ki/390.0_ki*known(ig+65) - 2.0_ki/13.0&
           &_ki*known(ig+98) + 1.0_ki/78.0_ki*known(ig+102) + 4.0_ki/65.0_ki*&
           & known(ig+106) - 16.0_ki/65.0_ki*known(ig+110) + 32.0_ki/195.0_ki*&
           & known(ig+114) + 904.0_ki/585.0_ki*known(ig+160) + 1016.0_ki/585.0_ki*&
           & known(ig+164) - 256.0_ki/2925.0_ki*known(ig+168) + 32.0_ki/585.0_ki*&
           & known(ig+172) - 16.0_ki/325.0_ki*known(ig+176) - 2672.0_ki/585.0_ki*&
           & known(ig+180) - 560.0_ki/117.0_ki*known(ig+184) - 256.0_ki/2925.0_ki*&
           & known(ig+188) - 64.0_ki/585.0_ki*known(ig+192) + 64.0_ki/325.0_ki*&
           & known(ig+196)

      cg(is3(ig,18))  = - 4.0_ki/65.0_ki*known(ig+31) - 6.0_ki/65.0_ki*known(ig+57) - 1.&
           &0_ki/15.0_ki*known(ig+61) - 7.0_ki/65.0_ki*known(ig+65) - 1.0_ki/13.0_k&
           &i*known(ig+98) - 1.0_ki/13.0_ki*known(ig+102) + 2.0_ki/65.0_ki*known(ig+&
           & 106) - 8.0_ki/65.0_ki*known(ig+110) + 16.0_ki/195.0_ki*known(ig+114)&
           &  + 64.0_ki/195.0_ki*known(ig+160) + 256.0_ki/195.0_ki*known(ig+164) + &
           & 3424.0_ki/975.0_ki*known(ig+168) + 92.0_ki/195.0_ki*known(ig+172) - 8.0&
           &_ki/325.0_ki*known(ig+176) - 272.0_ki/195.0_ki*known(ig+180) - 128.0_ki&
           &/39.0_ki*known(ig+184) - 6976.0_ki/975.0_ki*known(ig+188) - 184.0_ki/&
           & 195.0_ki*known(ig+192) + 32.0_ki/325.0_ki*known(ig+196)

      cg(is3(ig,19))  = 37.0_ki/390.0_ki*known(ig+31) - 7.0_ki/65.0_ki*known(ig+57) + 4.0&
           &_ki/45.0_ki*known(ig+61) + 8.0_ki/195.0_ki*known(ig+65) - 7.0_ki/78.0_k&
           &i*known(ig+98) - 7.0_ki/78.0_ki*known(ig+102) + 79.0_ki/390.0_ki*&
           & known(ig+106) - 28.0_ki/195.0_ki*known(ig+110) + 56.0_ki/585.0_ki*&
           & known(ig+114) - 556.0_ki/585.0_ki*known(ig+160) + 376.0_ki/585.0_ki*&
           & known(ig+164) + 6784.0_ki/2925.0_ki*known(ig+168) - 68.0_ki/585.0_ki*&
           & known(ig+172) - 28.0_ki/975.0_ki*known(ig+176) + 608.0_ki/585.0_ki*&
           & known(ig+180) - 80.0_ki/39.0_ki*known(ig+184) - 4672.0_ki/975.0_ki*&
           & known(ig+188) + 136.0_ki/585.0_ki*known(ig+192) + 112.0_ki/975.0_ki*&
           & known(ig+196)

      cg(is3(ig,20))  = 6.0_ki/65.0_ki*known(ig+31) + 9.0_ki/65.0_ki*known(ig+57) + 1.0_k&
           &i/10.0_ki*known(ig+61) + 21.0_ki/130.0_ki*known(ig+65) + 3.0_ki/26.0_ki&
           & *known(ig+98) + 3.0_ki/26.0_ki*known(ig+102) - 3.0_ki/65.0_ki*known(ig+&
           & 106) + 12.0_ki/65.0_ki*known(ig+110) - 8.0_ki/65.0_ki*known(ig+114) + &
           & 164.0_ki/195.0_ki*known(ig+160) + 136.0_ki/195.0_ki*known(ig+164) + 64.&
           &0_ki/975.0_ki*known(ig+168) - 8.0_ki/195.0_ki*known(ig+172) + 12.0_ki/&
           & 325.0_ki*known(ig+176) - 112.0_ki/195.0_ki*known(ig+180) - 16.0_ki/39.0&
           &_ki*known(ig+184) + 64.0_ki/975.0_ki*known(ig+188) + 16.0_ki/195.0_ki*&
           & known(ig+192) - 48.0_ki/325.0_ki*known(ig+196)
      
     enddo
    

!---------------------------------------------------------------------------
      do igs=69,209
         known(igs) = xneval(igs)-sub1(igs)-sub2(igs)-sub3(igs)-&
                                  sub4(igs)-sub5(igs)-sub6(igs)
      enddo
!  four q 
      cg(69)  = (1082.0_ki/369.0_ki*known(69) + 215.0_ki/246.0_ki*known(122)&
           &  + 2617.0_ki/246.0_ki*known(123) - 1751.0_ki/246.0_ki*known(124)&
           &  - 73.0_ki/82.0_ki*known(125) - 472.0_ki/41.0_ki*known(200) - &
           & 616.0_ki/123.0_ki*known(201) - 9352.0_ki/369.0_ki*known(202) + &
           & 7184.0_ki/369.0_ki*known(203) + 48.0_ki/41.0_ki*known(204) - 728.&
           &0_ki/369.0_ki*known(205) + 48.0_ki/41.0_ki*known(206) - 1360.0_ki&
           &/369.0_ki*known(207) - 18.0_ki/41.0_ki*known(208) + 52.0_ki/123.0&
           &_ki*known(209))
      
      cg(122)  = ( - 505.0_ki/246.0_ki*known(69) + 429.0_ki/82.0_ki*known(&
           & 122) + 20.0_ki/41.0_ki*known(123) + 200.0_ki/41.0_ki*known(124)&
           &  + 40.0_ki/41.0_ki*known(125) - 1280.0_ki/41.0_ki*known(200) - &
           & 744.0_ki/41.0_ki*known(201) - 680.0_ki/123.0_ki*known(202) - &
           & 1832.0_ki/123.0_ki*known(203) - 120.0_ki/41.0_ki*known(204) + &
           & 224.0_ki/123.0_ki*known(205) - 120.0_ki/41.0_ki*known(206) + 40.0&
           &_ki/123.0_ki*known(207) + 4.0_ki/41.0_ki*known(208) - 16.0_ki/41.&
           &0_ki*known(209))
      
      cg(123)  = ( - 355.0_ki/246.0_ki*known(69) - 26.0_ki/41.0_ki*known(122&
           & ) - 21.0_ki/82.0_ki*known(123) - 23.0_ki/41.0_ki*known(124) - 21.&
           &0_ki/41.0_ki*known(125) + 344.0_ki/41.0_ki*known(200) + 120.0_ki/&
           & 41.0_ki*known(201) - 176.0_ki/123.0_ki*known(202) + 232.0_ki/123.&
           &0_ki*known(203) + 104.0_ki/41.0_ki*known(204) - 52.0_ki/123.0_ki&
           & *known(205) + 104.0_ki/41.0_ki*known(206) + 184.0_ki/123.0_ki*&
           & known(207) + 2.0_ki/41.0_ki*known(208) - 8.0_ki/41.0_ki*known(&
           & 209))
      
      cg(124)  = ( - 2905.0_ki/738.0_ki*known(69) - 1223.0_ki/246.0_ki*&
           & known(122) - 1583.0_ki/123.0_ki*known(123) + 1021.0_ki/123.0_ki*&
           & known(124) + 38.0_ki/41.0_ki*known(125) + 1408.0_ki/41.0_ki*&
           & known(200) + 802.0_ki/41.0_ki*known(201) + 12904.0_ki/369.0_ki*&
           & known(202) - 8780.0_ki/369.0_ki*known(203) - 424.0_ki/123.0_ki*&
           & known(204) + 827.0_ki/369.0_ki*known(205) - 32.0_ki/41.0_ki*&
           & known(206) + 1180.0_ki/369.0_ki*known(207) + 12.0_ki/41.0_ki*&
           & known(208) + 20.0_ki/123.0_ki*known(209))
      
      cg(125)  = ( - 811.0_ki/738.0_ki*known(69) + 31.0_ki/246.0_ki*known(&
           & 122) - 878.0_ki/123.0_ki*known(123) + 767.0_ki/246.0_ki*known(&
           & 124) + 155.0_ki/82.0_ki*known(125) + 472.0_ki/41.0_ki*known(200)&
           &  + 14.0_ki/41.0_ki*known(201) + 6400.0_ki/369.0_ki*known(202) - &
           & 2756.0_ki/369.0_ki*known(203) - 472.0_ki/123.0_ki*known(204) + &
           & 113.0_ki/369.0_ki*known(205) - 48.0_ki/41.0_ki*known(206) + 868.0&
           &_ki/369.0_ki*known(207) + 18.0_ki/41.0_ki*known(208) - 52.0_ki/&
           & 123.0_ki*known(209))

      cg(196)  = (1.0_ki/3.0_ki*known(69) - known(122) + 8.0_ki/3.0_ki*&
           & known(201))
      
      cg(197)  = (1.0_ki/3.0_ki*known(69) - known(123) + 8.0_ki/3.0_ki*&
           & known(202))
      
      cg(198)  = (1.0_ki/3.0_ki*known(69) - known(124) + 8.0_ki/3.0_ki*&
           & known(203))
      
      cg(199)  = (1.0_ki/3.0_ki*known(69) - known(125) + 8.0_ki/3.0_ki*&
           & known(204))
      
      cg(204)  = ( - 2.0_ki*known(122) - 2.0_ki*known(123) + 16.0_ki*known(&
           & 200) + 8.0_ki*known(201) + 8.0_ki*known(202))
      
      cg(205)  = (4.0_ki*known(69) + known(122) - 5.0_ki*known(124) - 4.0_ki&
           & *known(201) + 16.0_ki*known(203) - 2.0_ki*known(205))
      
      cg(206)  = ( - 178.0_ki/123.0_ki*known(69) - 153.0_ki/41.0_ki*known(&
           & 122) + 62.0_ki/41.0_ki*known(123) + 5.0_ki/41.0_ki*known(124) - &
           & 40.0_ki/41.0_ki*known(125) + 624.0_ki/41.0_ki*known(200) + 580.0_&
           &ki/41.0_ki*known(201) - 304.0_ki/123.0_ki*known(202) - 136.0_ki/&
           & 123.0_ki*known(203) + 120.0_ki/41.0_ki*known(204) + 22.0_ki/123.0&
           &_ki*known(205) + 120.0_ki/41.0_ki*known(206) - 40.0_ki/123.0_ki*&
           & known(207) - 4.0_ki/41.0_ki*known(208) + 16.0_ki/41.0_ki*known(&
           & 209))
      
      cg(207)  = ( - 1.0_ki/3.0_ki*known(69) + 3.0_ki/2.0_ki*known(122) + 5.0&
           &_ki*known(123) + 1.0_ki/2.0_ki*known(124) - 16.0_ki*known(200)&
           &  - 6.0_ki*known(201) - 40.0_ki/3.0_ki*known(202) - 4.0_ki/3.0_ki&
           & *known(203) + 1.0_ki/3.0_ki*known(205) - 4.0_ki/3.0_ki*known(207&
           & ))
      
      cg(208)  = (280.0_ki/123.0_ki*known(69) + 93.0_ki/82.0_ki*known(122)&
           &  - 92.0_ki/41.0_ki*known(123) + 5.0_ki/82.0_ki*known(124) + 21.0_&
           &ki/41.0_ki*known(125) - 344.0_ki/41.0_ki*known(200) - 202.0_ki/&
           & 41.0_ki*known(201) + 832.0_ki/123.0_ki*known(202) - 68.0_ki/123.0&
           &_ki*known(203) - 104.0_ki/41.0_ki*known(204) + 11.0_ki/123.0_ki*&
           & known(205) - 104.0_ki/41.0_ki*known(206) - 20.0_ki/123.0_ki*&
           & known(207) - 2.0_ki/41.0_ki*known(208) + 8.0_ki/41.0_ki*known(&
           & 209))

      cg(209)  = (284.0_ki/369.0_ki*known(69) + 304.0_ki/123.0_ki*known(122)&
           &  + 968.0_ki/123.0_ki*known(123) - 406.0_ki/123.0_ki*known(124)&
           &  - 38.0_ki/41.0_ki*known(125) - 752.0_ki/41.0_ki*known(200) - &
           & 392.0_ki/41.0_ki*known(201) - 7984.0_ki/369.0_ki*known(202) + &
           & 3368.0_ki/369.0_ki*known(203) + 424.0_ki/123.0_ki*known(204) - &
           & 212.0_ki/369.0_ki*known(205) + 32.0_ki/41.0_ki*known(206) - 688.0&
           &_ki/369.0_ki*known(207) - 12.0_ki/41.0_ki*known(208) - 20.0_ki/&
           & 123.0_ki*known(209))
      
      
      if (verbosity .ge. 1) then
         do ig=0,209
            print*, ig, cg(ig)
         enddo
      end if
 end subroutine solve_system_rk6


subroutine     init_qg_list
   implicit none
   

  
      qg(0,:) = (/ czip, czip, czip, czip /)

      qg(1,:) = (/ cone, czip, czip, czip /)
      qg(2,:) = (/ czip, cone, czip, czip /)
      qg(3,:) = (/ czip, czip, cone, czip /)
      qg(4,:) = (/ czip, czip, czip, cone /)

      qg(5,:) = (/ -cone, czip, czip, czip /)
      qg(6,:) = (/  czip,-cone, czip, czip /)
      qg(7,:) = (/  czip, czip,-cone, czip /)
      qg(8,:) = (/  czip, czip, czip,-cone /)

      qg(9,:) = (/ cone, cone, czip, czip /)
      qg(10,:) =(/ cone, czip, cone, czip /) 
      qg(11,:) =(/ cone, czip, czip, cone /) 
      qg(12,:) =(/ czip, cone, cone, czip /) 
      qg(13,:) =(/ czip, cone, czip, cone /) 
      qg(14,:) =(/ czip, czip, cone, cone /) 

      qg(15,:) = (/ chaf, czip, czip, czip /)
      qg(16,:) = (/ czip, chaf, czip, czip /)
      qg(17,:) = (/ czip, czip, chaf, czip /)
      qg(18,:) = (/ czip, czip, czip, chaf   /)

      qg(19,:) = (/ cone, -cone, czip, czip /)
      qg(20,:) = (/ cone,  czip,-cone, czip /)
      qg(21,:) = (/ cone,  czip, czip,-cone /)
      qg(22,:) = (/ czip,  cone,-cone, czip /)
      qg(23,:) = (/ czip,  cone, czip,-cone /)
      qg(24,:) = (/ czip,  czip, cone,-cone /)

      qg(25,:) = (/ chaf, cone, czip, czip /)
      qg(26,:) = (/ chaf, czip, cone, czip /)
      qg(27,:) = (/ chaf, czip, czip, cone /)
      qg(28,:) = (/ czip, chaf, cone, czip /)
      qg(29,:) = (/ czip, chaf, czip, cone /)
      qg(30,:) = (/ czip, czip, chaf, cone /)

      qg(31,:) = (/ cone, cone, cone, czip /)
      qg(32,:) = (/ cone, cone, czip, cone /)
      qg(33,:) = (/ cone, czip, cone, cone /)
      qg(34,:) = (/ czip, cone, cone, cone /)

      qg(35,:) = (/ -chaf, czip, czip, czip /)
      qg(36,:) = (/  czip,-chaf, czip, czip /)
      qg(37,:) = (/  czip, czip,-chaf, czip /)
      qg(38,:) = (/  czip, czip, czip,-chaf /)

      qg(39,:) = (/ -chaf, cone, czip, czip /)
      qg(40,:) = (/ -chaf, czip, cone, czip /)
      qg(41,:) = (/ -chaf, czip, czip, cone /)
      qg(42,:) = (/  czip,-chaf, cone, czip /)
      qg(43,:) = (/  czip,-chaf, czip, cone /)
      qg(44,:) = (/  czip, czip,-chaf, cone /)

      qg(45,:) = (/ chaf, -chaf, czip, czip /)
      qg(46,:) = (/ chaf, czip, -chaf, czip /)
      qg(47,:) = (/ chaf, czip, czip, -chaf /)
      qg(48,:) = (/ czip, chaf, -chaf, czip /)
      qg(49,:) = (/ czip, chaf, czip, -chaf /)
      qg(50,:) = (/ czip, czip, chaf, -chaf /)

      qg(51,:) = (/ -cone, -cone, czip, czip /)
      qg(52,:) = (/ -cone,  czip,-cone, czip /)
      qg(53,:) = (/ -cone,  czip, czip,-cone /)
      qg(54,:) = (/  czip, -cone,-cone, czip /)
      qg(55,:) = (/  czip, -cone, czip,-cone /)
      qg(56,:) = (/  czip,  czip,-cone,-cone /)

      qg(57,:) = (/ -cone, cone, cone, czip /)
      qg(58,:) = (/ -cone, cone, czip, cone /)
      qg(59,:) = (/ -cone, czip, cone, cone /)
      qg(60,:) = (/  czip,-cone, cone, cone /)

      qg(61,:) = (/ cone, -cone, cone, czip /)
      qg(62,:) = (/ cone, -cone, czip, cone /) 
      qg(63,:) = (/ cone,  czip,-cone, cone /)
      qg(64,:) = (/ czip,  cone,-cone, cone /)

      qg(65,:) = (/ cone, cone, -cone,  czip /)
      qg(66,:) = (/ cone, cone,  czip, -cone /)
      qg(67,:) = (/ cone, czip,  cone, -cone /)
      qg(68,:) = (/ czip, cone,  cone, -cone /)

      qg(69,:) = (/ cone, cone, cone, cone /)

      qg(70,:) = (/ ctwo, czip, czip, czip /)
      qg(71,:) = (/ czip, ctwo, czip, czip /)
      qg(72,:) = (/ czip, czip, ctwo, czip /)
      qg(73,:) = (/ czip, czip, czip, ctwo  /)

      qg(74,:) = (/ -cone, cone, czip, czip /)
      qg(75,:) = (/ -cone, czip, cone, czip /) 
      qg(76,:) = (/ -cone, czip, czip, cone /) 
      qg(77,:) = (/ czip, -cone, cone, czip /)
      qg(78,:) = (/ czip, -cone, czip, cone /) 
      qg(79,:) = (/ czip, czip, -cone, cone /)

      qg(80,:) = (/ cone, chaf, czip, czip /)
      qg(81,:) = (/ cone, czip, chaf, czip /)
      qg(82,:) = (/ cone, czip, czip, chaf /)
      qg(83,:) = (/ czip, cone, chaf, czip /)
      qg(84,:) = (/ czip, cone, czip, chaf /)
      qg(85,:) = (/ czip, czip, cone, chaf /)

      qg(86,:) = (/ cone, -chaf, czip, czip /)
      qg(87,:) = (/ cone,  czip,-chaf, czip /)
      qg(88,:) = (/ cone,  czip, czip,-chaf /)
      qg(89,:) = (/ czip,  cone,-chaf, czip /)
      qg(90,:) = (/ czip,  cone, czip,-chaf /)
      qg(91,:) = (/ czip,  czip, cone,-chaf /)

      qg(92,:) = (/ -cone, chaf, czip, czip /)
      qg(93,:) = (/ -cone, czip, chaf, czip /)
      qg(94,:) = (/ -cone, czip, czip, chaf /)
      qg(95,:) = (/  czip,-cone, chaf, czip /)
      qg(96,:) = (/  czip,-cone, czip, chaf /)
      qg(97,:) = (/  czip, czip,-cone, chaf /)

      qg(98,:) = (/ -cone, -cone,  cone, czip /)
      qg(99,:) = (/ -cone, -cone,  czip, cone /) 
      qg(100,:) =(/ -cone,  czip, -cone, cone /)
      qg(101,:) =(/  czip, -cone, -cone, cone /)

      qg(102,:) = (/ -cone, cone,-cone,  czip /)
      qg(103,:) = (/ -cone, cone, czip, -cone /)
      qg(104,:) = (/ -cone, czip, cone, -cone /)
      qg(105,:) = (/  czip,-cone, cone, -cone /)

      qg(106,:) = (/ cone, -cone, -cone,  czip /)
      qg(107,:) = (/ cone, -cone,  czip, -cone /)
      qg(108,:) = (/ cone,  czip, -cone, -cone /)
      qg(109,:) = (/ czip,  cone, -cone, -cone /)

      qg(110,:) = (/ chaf, cone, cone, czip /)
      qg(111,:) = (/ chaf, cone, czip, cone /)
      qg(112,:) = (/ chaf, czip, cone, cone /)
      qg(113,:) = (/ czip, chaf, cone, cone /)

      qg(114,:) = (/ cone, chaf, cone, czip /)
      qg(115,:) = (/ cone, chaf, czip, cone /)
      qg(116,:) = (/ cone, czip, chaf, cone /)
      qg(117,:) = (/ czip, cone, chaf, cone /)

      qg(118,:) = (/ cone, cone, chaf, czip /)
      qg(119,:) = (/ cone, cone, czip, chaf /)
      qg(120,:) = (/ cone, czip, cone, chaf /)
      qg(121,:) = (/ czip, cone, cone, chaf /)

      qg(122,:) = (/ -cone, cone, cone, cone /)
      qg(123,:) = (/  cone,-cone, cone, cone /)
      qg(124,:) = (/  cone, cone,-cone, cone /)
      qg(125,:) = (/  cone, cone, cone,-cone /)

      qg(126,:) = (/ -ctwo, czip, czip, czip /)
      qg(127,:) = (/  czip,-ctwo, czip, czip /)
      qg(128,:) = (/  czip, czip,-ctwo, czip /)
      qg(129,:) = (/  czip, czip, czip,-ctwo /)

      qg(130,:) = (/ chaf,-cone, czip,  czip /)
      qg(131,:) = (/ chaf, czip,-cone,  czip  /)
      qg(132,:) = (/ chaf, czip, czip, -cone  /)
      qg(133,:) = (/ czip, chaf,-cone,  czip  /)
      qg(134,:) = (/ czip, chaf, czip, -cone  /)
      qg(135,:) = (/ czip, czip, chaf, -cone  /)

      qg(136,:) = (/ -chaf,-cone, czip,  czip /)
      qg(137,:) = (/ -chaf, czip,-cone,  czip /)
      qg(138,:) = (/ -chaf, czip, czip, -cone /)
      qg(139,:) = (/  czip,-chaf,-cone,  czip /) 
      qg(140,:) = (/  czip,-chaf, czip, -cone /) 
      qg(141,:) = (/  czip, czip,-chaf, -cone /)

      qg(142,:) = (/ -chaf, chaf, czip, czip  /)
      qg(143,:) = (/ -chaf, czip, chaf, czip  /) 
      qg(144,:) = (/ -chaf, czip, czip, chaf  /) 
      qg(145,:) = (/  czip,-chaf, chaf, czip  /)
      qg(146,:) = (/  czip,-chaf, czip, chaf  /)
      qg(147,:) = (/  czip, czip,-chaf, chaf  /)

      qg(148,:) = (/ chaf, ctwo, czip, czip  /)
      qg(149,:) = (/ chaf, czip, ctwo, czip  /)
      qg(150,:) = (/ chaf, czip, czip, ctwo  /)
      qg(151,:) = (/ czip, chaf, ctwo, czip  /)
      qg(152,:) = (/ czip, chaf, czip, ctwo  /)
      qg(153,:) = (/ czip, czip, chaf, ctwo  /)

      qg(154,:) = (/ ctwo, chaf, czip, czip  /)
      qg(155,:) = (/ ctwo, czip, chaf, czip  /)
      qg(156,:) = (/ ctwo, czip, czip, chaf  /)
      qg(157,:) = (/ czip, ctwo, chaf, czip  /)
      qg(158,:) = (/ czip, ctwo, czip, chaf  /)
      qg(159,:) = (/ czip, czip, ctwo, chaf  /)

      qg(160,:) = (/ chaf, -cone, -cone,  czip /)
      qg(161,:) = (/ chaf, -cone,  czip, -cone /)
      qg(162,:) = (/ chaf,  czip, -cone, -cone /)
      qg(163,:) = (/ czip,  chaf, -cone, -cone /)

      qg(164,:) = (/ -chaf, -cone, -cone,  czip /)
      qg(165,:) = (/ -chaf, -cone,  czip, -cone /)
      qg(166,:) = (/ -chaf,  czip, -cone, -cone /)
      qg(167,:) = (/  czip, -chaf, -cone, -cone /)

      qg(168,:) = (/ -chaf, chaf,-cone,  czip /)
      qg(169,:) = (/ -chaf, chaf, czip, -cone /)
      qg(170,:) = (/ -chaf, czip, chaf, -cone /)
      qg(171,:) = (/  czip,-chaf, chaf, -cone /)

      qg(172,:) = (/ chaf, ctwo,-cone,  czip /)
      qg(173,:) = (/ chaf, ctwo, czip, -cone /)
      qg(174,:) = (/ chaf, czip, ctwo, -cone /)
      qg(175,:) = (/ czip, chaf, ctwo, -cone /)

      qg(176,:) = (/ ctwo, chaf,-cone,  czip /)
      qg(177,:) = (/ ctwo, chaf, czip, -cone /)
      qg(178,:) = (/ ctwo, czip, chaf, -cone /)
      qg(179,:) = (/ czip, ctwo, chaf, -cone /)

      qg(180,:) = (/ chaf, -cone, -chaf,  czip /)
      qg(181,:) = (/ chaf, -cone,  czip, -chaf /)
      qg(182,:) = (/ chaf,  czip, -cone, -chaf /)
      qg(183,:) = (/ czip,  chaf, -cone, -chaf /)

      qg(184,:) = (/ -chaf,-cone, -chaf,  czip /)
      qg(185,:) = (/ -chaf,-cone,  czip, -chaf /)
      qg(186,:) = (/ -chaf, czip, -cone, -chaf /)
      qg(187,:) = (/  czip,-chaf, -cone, -chaf /)

      qg(188,:) = (/ -chaf, chaf,-chaf,  czip /)
      qg(189,:) = (/ -chaf, chaf, czip, -chaf /)
      qg(190,:) = (/ -chaf, czip, chaf, -chaf /)
      qg(191,:) = (/  czip,-chaf, chaf, -chaf /)

      qg(192,:) = (/ chaf, ctwo,-chaf,  czip /)
      qg(193,:) = (/ chaf, ctwo, czip, -chaf /)
      qg(194,:) = (/ chaf, czip, ctwo, -chaf /)
      qg(195,:) = (/ czip, chaf, ctwo, -chaf /)

      qg(196,:) = (/ ctwo, chaf,-chaf,  czip /)
      qg(197,:) = (/ ctwo, chaf, czip, -chaf /)
      qg(198,:) = (/ ctwo, czip, chaf, -chaf /)
      qg(199,:) = (/ czip, ctwo, chaf, -chaf /)

      qg(200,:) = (/  chaf, chaf, cone, cone /)
      qg(201,:) = (/ -chaf, cone, cone, cone /)
      qg(202,:) = (/  cone,-chaf, cone, cone /)
      qg(203,:) = (/  cone, cone,-chaf, cone /)
      qg(204,:) = (/  cone, cone, cone,-chaf /)
      qg(205,:) = (/  chaf, cone, ctwo, cone /)
      qg(206,:) = (/  chaf, ctwo, cone, chaf /)
      qg(207,:) = (/  chaf,-ctwo, chaf, cone /)
      qg(208,:) = (/  chaf,-chaf,-ctwo, ctwo /)
      qg(209,:) = (/  ctwo, ctwo, chaf,-chaf /)

   qg_list_is_initialized = .true.
end subroutine init_qg_list

subroutine     init_cg(numeval,rank,xneval,dxneval)
   implicit none

   integer, intent(in) :: rank
   complex(ki), dimension(0:209), intent(out) :: xneval
   complex(ki), dimension(0:9), intent(out) :: dxneval

   integer :: ig, igm

   interface
      function     numeval(ncut, Q, mu2)
        use precision
        implicit none
        integer, intent(in) :: ncut
        complex(ki), dimension(4), intent(in) :: Q
        complex(ki), intent(in) :: mu2
        complex(ki) :: numeval
      end function numeval
   end interface

   select case(rank)
   case(0:1)
      igm = 4
   case(2)
      igm = 14
   case(3)
      igm = 34
   case(4)
      igm = 69
   case(5)
      igm = 125
   case(6)
      igm = 209
   case default
      print*, "In init_qg: rank not yet implemented: rank =", rank
      stop
   end select

   ! Use icut = 9 to avoid confusing numerators which make use of
   ! the icut parameter (like golem-2.0). [TR]
   cg(:)=czip
   cgx(:)=czip
   cg(0) = numeval(9, qg(0,:),czip)
   do ig = 1, igm
      xneval(ig) = numeval(9, qg(ig,:),czip) - cg(0)
!      cg(ig) = czip
   enddo
! initialization of d-xneval   
   dxneval(0)=numeval(9,qg(0,:),cone)- cg(0)
   dxneval(5)=numeval(9,qg(0,:),-cone)- cg(0)
   do ig = 1, 4
      dxneval(ig)=numeval(9,qg(ig,:),cone)- cg(0)
      dxneval(5+ig)=numeval(9,qg(4+ig,:),cone)- cg(0)
   enddo
   
end subroutine init_cg

function numetens(ncut,Q,mu2)
  implicit none
  integer, intent(in) :: ncut
  complex(ki), dimension(4), intent(in) :: Q
  complex(ki), intent(in) :: mu2
  complex(ki) :: numetens
  complex(ki) :: sub1, sub2, sub3, sub4, sub5, sub6, subx
  complex(ki) :: q1,q2,q3,q4,a1,a2,a3
  !integer :: i

  q1 = Q(1)
  q2 = Q(2)
  q3 = Q(3)
  q4 = Q(4)

  a1 = q4*q4
  a2 = a1*a1
  a3 = a1*q4

  sub1=czip
  sub2=czip
  sub3=czip
  sub4=czip
  sub5=czip
  sub6=czip
  subx=czip

 if (myrank.eq.6) goto 1
 if (myrank.eq.5) goto 2
 if (myrank.eq.4) goto 3
 if (myrank.eq.3) goto 4
 if (myrank.eq.2) goto 5
 if (myrank.eq.1) goto 6
 if (myrank.eq.0) then
   numetens=cg(0)
   goto 7
 endif
  
 1 continue

 !rank6
   sub6 = ((((((cg(128)*q3+cg(140)*q4)*q3+cg(152)*a1)*q3+cg(171)*a1*q4)*q3+cg(15&
   &3)*a2)*q3+cg(141)*a2*q4)*q3+(((((cg(138)*q3+cg(164)*q4)*q3+cg(194)*a1)*q3+cg&
   &(195)*a1*q4)*q3+cg(165)*a2)*q3+((((cg(150)*q3+cg(192)*q4)*q3+cg(203)*a1)*q3+&
   &cg(193)*a1*q4)*q3+(((cg(169)*q3+cg(190)*q4)*q3+cg(191)*a1)*q3+((cg(148)*q3+c&
   &g(163)*q4)*q3+(cg(127)*q2+cg(136)*q3+cg(137)*q4)*q2+cg(149)*a1)*q2+cg(170)*a&
   &1*q4)*q2+cg(151)*a2)*q2+cg(139)*a2*q4)*q2+(((((cg(134)*q3+cg(161)*q4)*q3+cg(&
   &188)*a1)*q3+cg(189)*a1*q4)*q3+cg(162)*a2)*q3+((((cg(159)*q3+cg(198)*q4)*q3+c&
   &g(209)*a1)*q3+cg(199)*a1*q4)*q3+(((cg(186)*q3+cg(207)*q4)*q3+cg(208)*a1)*q3+&
   &((cg(184)*q3+cg(197)*q4)*q3+(cg(133)*q2+cg(157)*q3+cg(158)*q4)*q2+cg(185)*a1&
   &)*q2+cg(187)*a1*q4)*q2+cg(160)*a2)*q2+((((cg(146)*q3+cg(182)*q4)*q3+cg(202)*&
   &a1)*q3+cg(183)*a1*q4)*q3+(((cg(180)*q3+cg(205)*q4)*q3+cg(206)*a1)*q3+((cg(20&
   &0)*q3+cg(204)*q4)*q3+(cg(145)*q2+cg(178)*q3+cg(179)*q4)*q2+cg(201)*a1)*q2+cg&
   &(181)*a1*q4)*q2+(((cg(167)*q3+cg(176)*q4)*q3+cg(177)*a1)*q3+((cg(174)*q3+cg(&
   &196)*q4)*q3+(cg(166)*q2+cg(172)*q3+cg(173)*q4)*q2+cg(175)*a1)*q2+((cg(143)*q&
   &3+cg(156)*q4)*q3+(cg(142)*q2+cg(154)*q3+cg(155)*q4)*q2+(cg(126)*q1+cg(130)*q&
   &2+cg(131)*q3+cg(132)*q4)*q1+cg(144)*a1)*q1+cg(168)*a1*q4)*q1+cg(147)*a2)*q1+&
   &cg(135)*a2*q4)*q1+cg(129)*a3*a3)

 2 continue

!rank5
   sub5 = (((((cg(72)*q3+cg(84)*q4)*q3+cg(96)*a1)*q3+cg(97)*a1*q4)*q3+cg(85)*a2)&
   &*q3+((((cg(82)*q3+cg(108)*q4)*q3+cg(121)*a1)*q3+cg(109)*a1*q4)*q3+(((cg(94)*&
   &q3+cg(119)*q4)*q3+cg(120)*a1)*q3+((cg(92)*q3+cg(107)*q4)*q3+(cg(71)*q2+cg(80&
   &)*q3+cg(81)*q4)*q2+cg(93)*a1)*q2+cg(95)*a1*q4)*q2+cg(83)*a2)*q2+((((cg(78)*q&
   &3+cg(105)*q4)*q3+cg(118)*a1)*q3+cg(106)*a1*q4)*q3+(((cg(103)*q3+cg(124)*q4)*&
   &q3+cg(125)*a1)*q3+((cg(116)*q3+cg(123)*q4)*q3+(cg(77)*q2+cg(101)*q3+cg(102)*&
   &q4)*q2+cg(117)*a1)*q2+cg(104)*a1*q4)*q2+(((cg(90)*q3+cg(114)*q4)*q3+cg(115)*&
   &a1)*q3+((cg(112)*q3+cg(122)*q4)*q3+(cg(89)*q2+cg(110)*q3+cg(111)*q4)*q2+cg(1&
   &13)*a1)*q2+((cg(87)*q3+cg(100)*q4)*q3+(cg(86)*q2+cg(98)*q3+cg(99)*q4)*q2+(cg&
   &(70)*q1+cg(74)*q2+cg(75)*q3+cg(76)*q4)*q1+cg(88)*a1)*q1+cg(91)*a1*q4)*q1+cg(&
   &79)*a2)*q1+cg(73)*a2*q4)

 3 continue

!rank4
   sub4 = ((((cg(37)*q3+cg(49)*q4)*q3+cg(56)*a1)*q3+cg(50)*a1*q4)*q3+(((cg(47)*q&
   &3+cg(67)*q4)*q3+cg(68)*a1)*q3+((cg(54)*q3+cg(66)*q4)*q3+(cg(36)*q2+cg(45)*q3&
   &+cg(46)*q4)*q2+cg(55)*a1)*q2+cg(48)*a1*q4)*q2+(((cg(43)*q3+cg(64)*q4)*q3+cg(&
   &65)*a1)*q3+((cg(62)*q3+cg(69)*q4)*q3+(cg(42)*q2+cg(60)*q3+cg(61)*q4)*q2+cg(6&
   &3)*a1)*q2+((cg(52)*q3+cg(59)*q4)*q3+(cg(51)*q2+cg(57)*q3+cg(58)*q4)*q2+(cg(3&
   &5)*q1+cg(39)*q2+cg(40)*q3+cg(41)*q4)*q1+cg(53)*a1)*q1+cg(44)*a1*q4)*q1+cg(38&
   &)*a1*a1)

 4 continue

!rank3
   sub3 = (((cg(17)*q3+cg(29)*q4)*q3+cg(30)*a1)*q3+((cg(27)*q3+cg(34)*q4)*q3+(cg&
   &(16)*q2+cg(25)*q3+cg(26)*q4)*q2+cg(28)*a1)*q2+((cg(23)*q3+cg(33)*q4)*q3+(cg(&
   &22)*q2+cg(31)*q3+cg(32)*q4)*q2+(cg(15)*q1+cg(19)*q2+cg(20)*q3+cg(21)*q4)*q1+&
   &cg(24)*a1)*q1+cg(18)*a1*q4)

 5 continue

!rank2
   sub2 = ((cg(7)*q3+cg(14)*q4)*q3+(cg(6)*q2+cg(12)*q3+cg(13)*q4)*q2+(cg(5)*q1&
        & +cg(9)*q2+cg(10)*q3+cg(11)*q4)*q1+cg(8)*q4*q4)

 6 continue

!rank1
   sub1 = cg(1)*q1+cg(2)*q2+cg(3)*q3+cg(4)*q4

! extra
   subx = cgx(1)*mu2+cgx(2)*mu2*mu2+(cgx(3)*q1+cgx(4)*q2+cgx(5)*q3+cgx(6)*q4)*&
        mu2+(cgx(7)*q1*q1+cgx(8)*q2*q2+cgx(9)*q3*q3+cgx(10)*q4*q4)*mu2
! sum
   numetens=cg(0)+subx+sub1+sub2+sub3+sub4+sub5+sub6

 7 continue

 end function numetens

pure function sub1(is)
   implicit none
   complex(ki) :: sub1
   complex(ki) :: q1,q2,q3,q4
   integer, intent(in):: is

   q1 = qg(is,1)
   q2 = qg(is,2)
   q3 = qg(is,3)
   q4 = qg(is,4)

   sub1 = cg(1)*q1+cg(2)*q2+cg(3)*q3+cg(4)*q4
     
end function sub1

pure function subx(is,mu2)
   implicit none
   complex(ki) :: subx
   complex(ki) :: q1,q2,q3,q4
   complex(ki), intent(in) :: mu2
   integer, intent(in):: is

   q1 = qg(is,1)
   q2 = qg(is,2)
   q3 = qg(is,3)
   q4 = qg(is,4)

   subx = cgx(1)*mu2+cgx(2)*mu2*mu2+(cgx(3)*q1+cgx(4)*q2+cgx(5)*q3+cgx(6)*q4)*mu2+&
   (cgx(7)*q1*q1+cgx(8)*q2*q2+cgx(9)*q3*q3+cgx(10)*q4*q4)*mu2
     
end function subx



pure function sub2(is)
   ! generated: Do, 1 Jul 2010 16:15:15 +0200
   implicit none
   integer, intent(in) :: is
   complex(ki) :: sub2
   complex(ki) :: q1, q2, q3, q4

   q1 = qg(is,1)
   q2 = qg(is,2)
   q3 = qg(is,3)
   q4 = qg(is,4)

   sub2 = ((cg(7)*q3+cg(14)*q4)*q3+(cg(6)*q2+cg(12)*q3+cg(13)*q4)*q2+(cg(5)*q1&
        & +cg(9)*q2+cg(10)*q3+cg(11)*q4)*q1+cg(8)*q4*q4)
end  function sub2

pure function sub3(is)
   ! generated: Do, 1 Jul 2010 16:15:16 +0200
   implicit none
   integer, intent(in) :: is
   complex(ki) :: sub3
   complex(ki) :: q1, q2, q3, q4
   complex(ki) :: a1

   q1 = qg(is,1)
   q2 = qg(is,2)
   q3 = qg(is,3)
   q4 = qg(is,4)

   a1 = q4*q4
   sub3 = (((cg(17)*q3+cg(29)*q4)*q3+cg(30)*a1)*q3+((cg(27)*q3+cg(34)*q4)*q3+(cg&
   &(16)*q2+cg(25)*q3+cg(26)*q4)*q2+cg(28)*a1)*q2+((cg(23)*q3+cg(33)*q4)*q3+(cg(&
   &22)*q2+cg(31)*q3+cg(32)*q4)*q2+(cg(15)*q1+cg(19)*q2+cg(20)*q3+cg(21)*q4)*q1+&
   &cg(24)*a1)*q1+cg(18)*a1*q4)
end function sub3

pure function sub4(is)
   ! generated: Do, 1 Jul 2010 16:15:17 +0200
   implicit none
   integer, intent(in) :: is
   complex(ki) :: sub4
   complex(ki) :: q1, q2, q3, q4
   complex(ki) :: a1

   q1 = qg(is,1)
   q2 = qg(is,2)
   q3 = qg(is,3)
   q4 = qg(is,4)

   a1 = q4*q4
   sub4 = ((((cg(37)*q3+cg(49)*q4)*q3+cg(56)*a1)*q3+cg(50)*a1*q4)*q3+(((cg(47)*q&
   &3+cg(67)*q4)*q3+cg(68)*a1)*q3+((cg(54)*q3+cg(66)*q4)*q3+(cg(36)*q2+cg(45)*q3&
   &+cg(46)*q4)*q2+cg(55)*a1)*q2+cg(48)*a1*q4)*q2+(((cg(43)*q3+cg(64)*q4)*q3+cg(&
   &65)*a1)*q3+((cg(62)*q3+cg(69)*q4)*q3+(cg(42)*q2+cg(60)*q3+cg(61)*q4)*q2+cg(6&
   &3)*a1)*q2+((cg(52)*q3+cg(59)*q4)*q3+(cg(51)*q2+cg(57)*q3+cg(58)*q4)*q2+(cg(3&
   &5)*q1+cg(39)*q2+cg(40)*q3+cg(41)*q4)*q1+cg(53)*a1)*q1+cg(44)*a1*q4)*q1+cg(38&
   &)*a1*a1)
end  function sub4

pure function sub5(is)
   ! generated: Do, 1 Jul 2010 16:15:18 +0200
   implicit none
   integer, intent(in) :: is

   complex(ki) :: sub5
   complex(ki) :: q1, q2, q3, q4


   complex(ki) :: a1
   complex(ki) :: a2

   q1 = qg(is,1)
   q2 = qg(is,2)
   q3 = qg(is,3)
   q4 = qg(is,4)


   a1 = q4*q4
   a2 = a1*a1
   sub5 = (((((cg(72)*q3+cg(84)*q4)*q3+cg(96)*a1)*q3+cg(97)*a1*q4)*q3+cg(85)*a2)&
   &*q3+((((cg(82)*q3+cg(108)*q4)*q3+cg(121)*a1)*q3+cg(109)*a1*q4)*q3+(((cg(94)*&
   &q3+cg(119)*q4)*q3+cg(120)*a1)*q3+((cg(92)*q3+cg(107)*q4)*q3+(cg(71)*q2+cg(80&
   &)*q3+cg(81)*q4)*q2+cg(93)*a1)*q2+cg(95)*a1*q4)*q2+cg(83)*a2)*q2+((((cg(78)*q&
   &3+cg(105)*q4)*q3+cg(118)*a1)*q3+cg(106)*a1*q4)*q3+(((cg(103)*q3+cg(124)*q4)*&
   &q3+cg(125)*a1)*q3+((cg(116)*q3+cg(123)*q4)*q3+(cg(77)*q2+cg(101)*q3+cg(102)*&
   &q4)*q2+cg(117)*a1)*q2+cg(104)*a1*q4)*q2+(((cg(90)*q3+cg(114)*q4)*q3+cg(115)*&
   &a1)*q3+((cg(112)*q3+cg(122)*q4)*q3+(cg(89)*q2+cg(110)*q3+cg(111)*q4)*q2+cg(1&
   &13)*a1)*q2+((cg(87)*q3+cg(100)*q4)*q3+(cg(86)*q2+cg(98)*q3+cg(99)*q4)*q2+(cg&
   &(70)*q1+cg(74)*q2+cg(75)*q3+cg(76)*q4)*q1+cg(88)*a1)*q1+cg(91)*a1*q4)*q1+cg(&
   &79)*a2)*q1+cg(73)*a2*q4)
end  function sub5

pure function sub6(is)
   implicit none
   integer, intent(in) :: is

   complex(ki) :: sub6
   complex(ki) :: q1, q2, q3, q4


   complex(ki) :: a1
   complex(ki) :: a2
   complex(ki) :: a3

   q1 = qg(is,1)
   q2 = qg(is,2)
   q3 = qg(is,3)
   q4 = qg(is,4)


   a1 = q4*q4
   a2 = a1*a1
   a3 = a1*q4
   sub6 = ((((((cg(128)*q3+cg(140)*q4)*q3+cg(152)*a1)*q3+cg(171)*a1*q4)*q3+cg(15&
   &3)*a2)*q3+cg(141)*a2*q4)*q3+(((((cg(138)*q3+cg(164)*q4)*q3+cg(194)*a1)*q3+cg&
   &(195)*a1*q4)*q3+cg(165)*a2)*q3+((((cg(150)*q3+cg(192)*q4)*q3+cg(203)*a1)*q3+&
   &cg(193)*a1*q4)*q3+(((cg(169)*q3+cg(190)*q4)*q3+cg(191)*a1)*q3+((cg(148)*q3+c&
   &g(163)*q4)*q3+(cg(127)*q2+cg(136)*q3+cg(137)*q4)*q2+cg(149)*a1)*q2+cg(170)*a&
   &1*q4)*q2+cg(151)*a2)*q2+cg(139)*a2*q4)*q2+(((((cg(134)*q3+cg(161)*q4)*q3+cg(&
   &188)*a1)*q3+cg(189)*a1*q4)*q3+cg(162)*a2)*q3+((((cg(159)*q3+cg(198)*q4)*q3+c&
   &g(209)*a1)*q3+cg(199)*a1*q4)*q3+(((cg(186)*q3+cg(207)*q4)*q3+cg(208)*a1)*q3+&
   &((cg(184)*q3+cg(197)*q4)*q3+(cg(133)*q2+cg(157)*q3+cg(158)*q4)*q2+cg(185)*a1&
   &)*q2+cg(187)*a1*q4)*q2+cg(160)*a2)*q2+((((cg(146)*q3+cg(182)*q4)*q3+cg(202)*&
   &a1)*q3+cg(183)*a1*q4)*q3+(((cg(180)*q3+cg(205)*q4)*q3+cg(206)*a1)*q3+((cg(20&
   &0)*q3+cg(204)*q4)*q3+(cg(145)*q2+cg(178)*q3+cg(179)*q4)*q2+cg(201)*a1)*q2+cg&
   &(181)*a1*q4)*q2+(((cg(167)*q3+cg(176)*q4)*q3+cg(177)*a1)*q3+((cg(174)*q3+cg(&
   &196)*q4)*q3+(cg(166)*q2+cg(172)*q3+cg(173)*q4)*q2+cg(175)*a1)*q2+((cg(143)*q&
   &3+cg(156)*q4)*q3+(cg(142)*q2+cg(154)*q3+cg(155)*q4)*q2+(cg(126)*q1+cg(130)*q&
   &2+cg(131)*q3+cg(132)*q4)*q1+cg(144)*a1)*q1+cg(168)*a1*q4)*q1+cg(147)*a2)*q1+&
   &cg(135)*a2*q4)*q1+cg(129)*a3*a3)
end  function sub6

end module mtens
