module mgetc1
   use precision, only: ki
   use constants
   use options
   use mfunctions
   use mrestore
   implicit none

   private

   interface getc1
      module procedure getc1_rm
      module procedure getc1_cm
   end interface getc1

   public :: getc1

contains


!CM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine getc1_cm(numeval,nleg,rank,c1,cut1,q1,qt,Vi,msq)
      use mglobal, only: G0c, mu2g, MP12, mu2t, resit, denst, mu2test
      implicit none
      integer, intent(in) :: nleg, rank, cut1
      complex(ki), dimension(0:15), intent(out) :: c1
      complex(ki), dimension(15,4), intent(in) :: q1
      complex(ki), dimension(4), intent(in) :: qt
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      complex(ki), dimension(0:nleg-1), intent(in) :: msq

      integer :: i,m,n,j1,i1,i2,i3,i4,i5
      integer :: dicut5,dicut4,dicut3,dicut2,diff
      complex(ki), dimension(15) :: dens1,dens2,dens3,dens4,dens5,xneval
      complex(ki), dimension(0:14) :: f1
      complex(ki), dimension(15) :: resi5, resi4, resi3, resi2, known
      complex(ki) :: dens2t,dens3t,dens4t,dens5t
      logical evalres
      complex(ki), dimension(15)   :: mu2vec

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

      mu2test(1)=mu2t(1)
      j1=cut1

      resi2(:)=czip
      resi3(:)=czip
      resi4(:)=czip
      resi5(:)=czip
      known(:)=czip
      xneval(:)=czip
      dens1(:)=cone

      !---  for lnntest
      resit(1)=czip
      denst(1)=cone

     !---  for simplified sampling
     diff = nleg-rank

     if (diff.eq.1) then
        c1(1)=czip
        c1(2)=czip
        c1(3)=czip
        c1(4)=czip
         
        if (nleg.eq.5) then
           resi5(1)=res5(1,mu2g(1))
           resit(1)=res5(1,mu2t(1))
           goto 11
        elseif (nleg.eq.4) then
           resi4(1)=Res4(1,q1(1,:),mu2g(1))
           resit(1)=Res4(1,qt,mu2t(1))
           goto 21
        elseif (nleg.eq.3) then
           resi3(1)=Res3(1,q1(1,:),mu2g(1))
           resit(1)=Res3(1,qt,mu2t(1))
           goto 31
        elseif (nleg.eq.2) then
           resi2(1)=Res2(1,q1(1,:),mu2g(1))
           resit(1)=Res2(1,qt,mu2t(1))
           goto 41
         else
            !---#[ Contributo dei pentuple cuts:
            dicut5=1
            do i5=4,nleg-1
               do i4=3,i5-1
                  do i3=2,i4-1
                     do i2=1,i3-1
                        do i1=0,i2-1
                           dens5(1)=cone
                           dens5t=cone
                           evalres=.false.
                           
                           loop_10: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)&
                                      &.and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i).eq.(j1)) then
                                    dens5(1)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_10
                                 else
                                    dens5(1)=dens5(1)*denevalmu2(nleg,i,&
                                            &q1(1,:),Vi,msq,mu2g(1))
                                    dens5t=dens5t*denevalmu2(nleg,i,qt,Vi,msq,&
                                            &mu2t(1))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_10

                           if (evalres) then
                              resi5(1)=resi5(1)+dens5(1)*res5(dicut5,mu2g(1))
                              resit(1)=resit(1)+dens5(1)*res5(dicut5,mu2t(1))
                           endif
         
                           dicut5=dicut5+1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !---#] Contributo dei pentuple cuts:
         endif

 11      continue

         !---#[ Contributo dei quadruple cuts:
         dicut4=1
         do i4=3,nleg-1
            do i3=2,i4-1
               do i2=1,i3-1
                  do i1=0,i2-1
                     dens4(1)=cone
                     dens4t=cone
                     evalres=.false.
      
                     loop_20: do i=0,nleg-1
                        if ((i.ne.i1).and.(i.ne.i2) &
                                &  .and.(i.ne.i3).and.(i.ne.i4)) then
                           if (i.eq.j1) then
                              dens4(1)=czip
                              dens4t=czip
                              evalres=.false.
                              exit loop_20
                           else
                              dens4(1)=dens4(1)*denevalmu2(nleg,i,q1(1,:),&
                                       &Vi,msq,mu2g(1))
                              dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,&
                                       &mu2t(1))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_20

                     if (evalres) then
                        resi4(1)=resi4(1)+dens4(1)*Res4(dicut4,q1(1,:),mu2g(1))
                        resit(1)=resit(1)+dens4t*Res4(dicut4,qt,mu2t(1))
                     endif

                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei quadruple cuts:
 21      continue

         !---#[ Contributo dei Triple cuts:
         dicut3=1
         do i3=2,nleg-1
            do i2=1,i3-1
               do i1=0,i2-1
                  dens3(1)=cone
                  dens3t=cone
                  evalres=.false.
      
                  loop_30: do i=0,nleg-1
                     if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                        if (i.eq.j1) then
                           dens3(1)=czip
                           dens3t=czip
                           evalres=.false.
                           exit loop_30
                        else
                           dens3(1)=dens3(1)*denevalmu2(nleg,i,q1(1,:),&
                              &Vi,msq,mu2g(1))
                           dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                           evalres=.true.
                        endif
                     endif
                  enddo loop_30

                  if (evalres) then
                     resi3(1)=resi3(1)+dens3(1)*Res3(dicut3,q1(1,:),mu2g(1))
                     resit(1)=resit(1)+dens3t*Res3(dicut3,qt,mu2t(1))
                  endif

                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:

 31      continue

         !---#[ Contributo dei Double cuts:
         dicut2=1
         do i2=1,nleg-1
            do i1=0,i2-1
               dens2(1)=cone
               dens2t=cone
               evalres=.false.
      
               loop_40: do i=0,nleg-1
                  if ((i.ne.i1).and.(i.ne.i2)) then
                     if (i.eq.j1) then
                        dens2(1)=czip
                        dens2t=czip
                        evalres=.false.
                        exit loop_40
                     else
                        dens2(1)=dens2(1)*denevalmu2(nleg,i,q1(1,:),&
                             &Vi,msq,mu2g(1))
                        dens2t=dens2t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                        evalres=.true.
                     endif
                  endif
               enddo loop_40

               if (evalres) then
                  resi2(1)=resi2(1)+dens2(1)*Res2(dicut2,q1(1,:),mu2g(1))
                  resit(1)=resit(1)+dens2t*Res2(dicut2,qt,mu2t(1))
               endif

               dicut2=dicut2+1
            enddo
         enddo
         !---#] Contributo dei Double cuts:

 41      continue

         !---
         do i=0,nleg-1
            if (i.ne.j1) then
               dens1(1)=dens1(1)*denevalmu2(nleg,i,q1(1,:),Vi,msq,mu2g(1))
               denst(1)=denst(1)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
            endif
         enddo

         xneval(1)=numeval(cut1,q1(1,:),mu2g(1))

         if     (imeth.eq.'diag') then
            known(1)=(xneval(1)-resi5(1)-resi4(1)-resi3(1)-resi2(1))/dens1(1) 
         elseif (imeth.eq.'tree') then
            known(1)=xneval(1)-(resi5(1)+resi4(1)+resi3(1)+resi2(1))/dens1(1) 
         endif

         c1(0) = known(1)
      else
         !--- Decomposizione standard

         if (nleg.eq.5) then
            resi5(:)=res5(1,mu2g(1))
            resit(1)=res5(1,mu2t(1)) 
            goto 111
         elseif (nleg.eq.4) then
            do n=1,15         
               resi4(n)=Res4(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res4(1,qt,mu2t(1))
            goto 121
         elseif (nleg.eq.3) then
            do n=1,15
               resi3(n)=Res3(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res3(1,qt,mu2t(1))
            goto 131
         elseif (nleg.eq.2) then
            do n=1,15
               resi2(n)=Res2(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res2(1,qt,mu2t(1))
            goto 141
         else
            !---#[ Contributo dei pentuple cuts:
            dicut5=1
            do i5=4,nleg-1
               do i4=3,i5-1
                  do i3=2,i4-1
                     do i2=1,i3-1
                        do i1=0,i2-1
                           dens5(:)=cone
                           dens5t=cone
                           evalres=.false.
      
                           loop_110: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)&
                                     &.and.(i.ne.i4).and.(i.ne.i5)) then
                                 if (i.eq.j1) then
                                    dens5(:)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_110
                                 else
                                    do n=1,15
                                       dens5(n)=dens5(n)*denevalmu2(nleg,i,&
                                           &q1(n,:),Vi,msq,mu2g(1))
                                    enddo
                                    dens5t=dens5t*denevalmu2(nleg,&
                                           &i,qt,Vi,msq,mu2t(1))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_110

                           if (evalres) then
                              resi5(:)=resi5(:)+dens5(:)*res5(dicut5,mu2g(1))
                              resit(1)=resit(1)+dens5t*res5(dicut5,mu2t(1))
                           endif

                           dicut5=dicut5+1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !---#] Contributo dei pentuple cuts:
         endif

 111     continue

         !---#[ Contributo dei quadruple cuts:
         dicut4=1
         do i4=3,nleg-1
            do i3=2,i4-1
               do i2=1,i3-1
                  do i1=0,i2-1
                     dens4(:)=cone
                     dens4t=cone
                     evalres=.false.
      
                     loop_120: do i=0,nleg-1
                        if ((i.ne.i1).and.(i.ne.i2) &
                               &  .and.(i.ne.i3).and.(i.ne.i4)) then
                           if (i.eq.j1) then
                              dens4(:)=czip
                              dens4t=czip
                              evalres=.false.
                              exit loop_120
                           else
                              do n=1,15
                                 dens4(n)=dens4(n)*denevalmu2(nleg,i,&
                                      &q1(n,:),Vi,msq,mu2g(1))
                              enddo
                              dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_120

                     if (evalres) then
                        do n=1,15
                           resi4(n)=resi4(n)+dens4(n)*&
                               &Res4(dicut4,q1(n,:),mu2g(1))
                        enddo
                        resit(1)=resit(1)+dens4t*Res4(dicut4,qt,mu2t(1))
                     endif
                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei quadruple cuts:

 121     continue

         !---#[ Contributo dei Triple cuts:
         dicut3=1
         do i3=2,nleg-1
            do i2=1,i3-1
               do i1=0,i2-1
                  dens3(:)=cone
                  dens3t=cone
                  evalres=.false.
      
                  loop_130: do i=0,nleg-1
                     if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                        if (i.eq.j1) then
                           dens3(:)=czip
                           dens3t=czip
                           evalres=.false.
                           exit loop_130
                        else
                           do n=1,15
                              dens3(n)=dens3(n)*denevalmu2(nleg,i,q1(n,:),&
                                      &Vi,msq,mu2g(1))
                           enddo
                           dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                           evalres=.true.
                        endif
                     endif
                  enddo loop_130

                  if (evalres) then
                     do n=1,15
                        resi3(n)=resi3(n)+dens3(n)*Res3(dicut3,q1(n,:),mu2g(1))
                     enddo
                     resit(1)=resit(1)+dens3t*Res3(dicut3,qt,mu2t(1))
                  endif

                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:

 131     continue

         !---#[ Contributo dei Double cuts:
         dicut2=1
         do i2=1,nleg-1
            do i1=0,i2-1
               dens2(:)=cone
               dens2t=cone
               evalres=.false.
   
               loop_140: do i=0,nleg-1
                  if ((i.ne.i1).and.(i.ne.i2)) then
                     if (i.eq.j1) then
                        dens2(:)=czip
                        dens2t=czip
                        evalres=.false.
                        exit loop_140
                     else
                        do n=1,15
                           dens2(n)=dens2(n)*denevalmu2(nleg,i,q1(n,:),&
                              &Vi,msq,mu2g(1))
                        enddo
                        dens2t=dens2t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                        evalres=.true.
                     endif
                  endif
               enddo loop_140

               if (evalres) then
                  do n=1,15
                     resi2(n)=resi2(n)+dens2(n)*Res2(dicut2,q1(n,:),mu2g(1))
                  enddo
                  resit(1)=resit(1)+dens2t*Res2(dicut2,qt,mu2t(1))
               endif
               dicut2=dicut2+1
            enddo
         enddo

         !---#] Contributo dei Double cuts:
 141     continue

         !---

         do i=0,nleg-1
            if (i.ne.j1) then
               do n=1,15
                  dens1(n)=dens1(n)*denevalmu2(nleg,i,q1(n,:),Vi,msq,mu2vec(n))
               enddo
               denst(1)=denst(1)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
            endif
         enddo

         do n=1,15
            xneval(n)=numeval(cut1,q1(n,:),mu2vec(n))
         enddo

         if     (imeth.eq.'diag') then
            known(:)=(xneval(:)-resi5(:)-resi4(:)-resi3(:)-resi2(:))/dens1(:)
         elseif (imeth.eq.'tree') then
            known(:)=xneval(:)-(resi5(:)+resi4(:)+resi3(:)+resi2(:))/dens1(:)
         endif

         do m=0,4
            f1(m)=known(m+1)
         enddo

         !--- Coefficienti
         c1(0) = (f1(0)+f1(1))/2.0_ki
         c1(1) = (-f1(0)-f1(1)+f1(3)+f1(4))/(2.0_ki*MP12(1))
         c1(2) = -(-2.0_ki*f1(0) + f1(3) + f1(4))/(2.0_ki*G0c*MP12(1))
         c1(3) = (2.0_ki*f1(0)-2.0_ki*f1(2)-f1(3)+f1(4))/(2.0_ki*G0c*MP12(1))
         c1(4) = (f1(0)-f1(2))/MP12(1)
      endif
   end subroutine getc1_cm
!------------------CM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!]






















































subroutine getc1_rm(numeval,nleg,rank,c1,cut1,q1,qt,Vi,msq)
	use mglobal, only: G0, mu2g, MP12, mu2t, resit, denst, mu2test, G0mu
	implicit none
! PARAMETERS===========================================================================================!
	! external
	integer,			    intent(in ) :: nleg, rank, cut1
	complex(ki), dimension(15,4),       intent(in ) :: q1
	complex(ki), dimension(4),          intent(in ) :: qt
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	real(ki),    dimension(0:nleg-1),   intent(in ) :: msq
	complex(ki), dimension(0:15),       intent(out) :: c1 
		!c1(9) does not exist, hence 15 instead of 14
	! internal
	integer 		     :: i, m, n, j1, i1, i2, i3, i4, i5, nold, nsol,nloops
	integer 		     :: dicut5, dicut4, dicut3, dicut2, diff
	complex(ki), dimension(15)   :: dens1, dens2, dens3, dens4, dens5, xneval
	complex(ki), dimension(0:14) :: f1
	complex(ki), dimension(15)   :: resi5, resi4, resi3, resi2, known, mu2vec
	complex(ki) 		     :: dens2t, dens3t, dens4t, dens5t
	logical evalres
! INTERFACE============================================================================================!
	interface
		function     numeval(ncut, Q, mu2)
		use precision
		implicit none
			integer, 		   intent(in) :: ncut
			complex(ki), dimension(4), intent(in) :: Q
			complex(ki), 	     	   intent(in) :: mu2
			complex(ki) 			      :: numeval
		end function numeval
	end interface
! INITIALIZATION=======================================================================================!
		do n=1,14
			mu2vec(n)=czip
		enddo
		mu2vec(15) = mu2g(1)
		mu2test(1) = mu2t(1)
		j1 = cut1
		resi2(:)  = czip
		resi3(:)  = czip
		resi4(:)  = czip
		resi5(:)  = czip
		known(:)  = czip
		xneval(:) = czip
		dens1(:)  = cone
	!---  for lnntest
		resit(1)=czip
		denst(1)=cone
		diff = nleg-rank
		!--
		if(diff .ge. 0) then
		do n=1,5
			mu2vec(n)=mu2g(1)
		enddo
		endif
		!--
		if 	(diff .ge. 1) then ! rank0: 1 coefficients
			nsol =  1
		elseif	(diff .eq. 0) then ! rank1:  5 coefficients
			nsol =  5
		elseif	(diff .eq.-1) then ! rank2: 15 coefficients
			nsol = 15
		else
			print *, 'error in mgetc1.f90'
		endif
! SIMPLIFIED SAMPLING==================================================================================!
! STANDARD DECOMPOSITION===============================================================================!
	if (nleg.eq.5) then
		do n=1,nsol         
			resi5(n)=Res5(1,mu2vec(n))
		enddo
		resit(1)=res5(1,mu2t(1)) 
	elseif (nleg.eq.4) then
		do n=1,nsol         
			resi4(n)=Res4(1,q1(n,:),mu2vec(n))
		enddo
		resit(1)=Res4(1,qt,mu2t(1))
	elseif (nleg.eq.3) then
		do n=1,nsol
			resi3(n)=Res3(1,q1(n,:),mu2vec(n))
		enddo
		resit(1)=Res3(1,qt,mu2t(1))
	elseif (nleg.eq.2) then
		do n=1,nsol
			resi2(n)=Res2(1,q1(n,:),mu2vec(n))
		enddo
		resit(1)=Res2(1,qt,mu2t(1))
	else
! Pentuple cuts contribution---------------------------------------------------------------------------!
		nloops = 5
		call HResL1(resi5,	nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	endif
! Quadrupole cuts contribution-------------------------------------------------------------------------!
	if((nleg .ne. 4) .and. (nleg .ne. 3) .and. (nleg .ne. 2)) then
			nloops = 4
			call HResL1(resi4,	nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	endif
! Triple cuts contribution-----------------------------------------------------------------------------!
	if((nleg .ne. 3) .and. (nleg .ne. 2)) then
			nloops = 3
			call HResL1(resi3,	nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	endif
! Double cuts contribution ----------------------------------------------------------------------------!
	if(nleg .ne. 2) then
			nloops = 2
			call HResL1(resi2,	nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	endif
! Single cut-------------------------------------------------------------------------------------------!
	do i=0,nleg-1
		if (i.ne.j1) then
			do n=1,nsol
			dens1(n)=dens1(n)*denevalmu2(nleg,i,q1(n,:),Vi,msq,mu2vec(n))
			enddo
			denst(1)=denst(1)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
		endif
	enddo
! CALCULATING THE POLYNOMIALS==========================================================================!
	do n=1,nsol
		xneval(n)=numeval(cut1,q1(n,:),mu2vec(n))
	enddo
	if     (imeth.eq.'diag') then
		do n=1,nsol
			known(n)=(xneval(n)-resi5(n)-resi4(n)-resi3(n)-resi2(n))/dens1(n)
		enddo
	elseif (imeth.eq.'tree') then
		do n=1,nsol
			known(n)=xneval(n)-(resi5(n)+resi4(n)+resi3(n)+resi2(n))/dens1(n)
		enddo
	endif
! SELECT BY RANK=======================================================================================!
	select case(diff)
! RANK 0===============================================================================================!
	case(1)
		c1(0) = known(1)
! RANK 1===============================================================================================!
	case(0)
		do m=0,4
			f1(m)=known(m+1)
		enddo
		c1(0) = (f1(0)+f1(1))/2.0_ki
		c1(1) = (-f1(0)-f1(1)+f1(3)+f1(4))/(2.0_ki*MP12(1))
		c1(2) = -(-2.0_ki*f1(0) + f1(3) + f1(4))/(2.0_ki*G0*MP12(1))
		c1(3) = (2.0_ki*f1(0)-2.0_ki*f1(2)-f1(3)+f1(4))/(2.0_ki*G0*MP12(1))
		c1(4) = (f1(0)-f1(2))/MP12(1)
! RANK 2===============================================================================================!
	case(-1)
	nold = 0
	if(abs(G0) .lt. 1.0e-10_ki) then
	! defining fs with G0=0'
		call DFT1(3,	nold,known,f1)
		call DFT1(2,	nold,known,f1)
		call DFT1(2,	nold,known,f1)
		call DFT1(2,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
	!~~~
	elseif(abs(G0) .lt. 1.0e-1_ki) then
		call DFT1(3,	nold,known,f1)
		call DFT1(2,	nold,known,f1)
		call DFT1(3,	nold,known,f1)
		call DFT1(2,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
	!~~~
	else 
	! defining fs with G0 nonzero'
		call DFT1(5,	nold,known,f1)
		call DFT1(5,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
		call DFT1(1,	nold,known,f1)
	endif
! Coefficients~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
if(abs(G0) .lt. 1.0e-10_ki) then 
      c1(0) =  &
     &f1(0) 
      c1(2) =  &
     &f1(1)/MP12(1) 
      c1(6) =  &
     &f1(2)/MP12(1)**2 
      c1(1) =  &
     &f1(4)/MP12(1) 
      c1(5) =  &
     &(-c1(0) + f1(3))/MP12(1)**2 
      c1(4) =  &
     &-(f1(6)/MP12(1)) 
      c1(8) =  &
     &(-c1(0) + f1(5))/MP12(1)**2 
      c1(3) =  &
     &-(f1(8)/MP12(1)) 
      c1(7) =  &
     &(-c1(0) + f1(7))/MP12(1)**2 
      c1(10) =  &
     &        (c1(0) - f1(9) + c1(1)*MP12(1) - c1(3)*MP12(1) +  &
     &    c1(5)*MP12(1)**2 + c1(7)*MP12(1)**2)/MP12(1)**2 
      c1(11) =  &
     &        (c1(0) - f1(10) + c1(1)*MP12(1) - c1(4)*MP12(1) +  &
     &    c1(5)*MP12(1)**2 + c1(8)*MP12(1)**2)/MP12(1)**2 
      c1(12) =  &
     &        (c1(0) - f1(11) + c1(2)*MP12(1) - c1(3)*MP12(1) +  &
     &    c1(6)*MP12(1)**2 + c1(7)*MP12(1)**2)/MP12(1)**2 
      c1(13) =  &
     &        (c1(0) - f1(12) + c1(2)*MP12(1) - c1(4)*MP12(1) +  &
     &    c1(6)*MP12(1)**2 + c1(8)*MP12(1)**2)/MP12(1)**2 
      c1(14) =  &
     &     (-c1(0) + f1(14) - G0mu*c1(1)*MP12(1) - c1(2)*MP12(1) -  &
     &  G0mu**2*c1(5)*MP12(1)**2 - c1(6)*MP12(1)**2)/mu2g(1)
      c1(15) =  &
     &        (-c1(0) + f1(13) - c1(1)*MP12(1) - c1(2)*MP12(1) +  &
     &    c1(3)*MP12(1) + c1(4)*MP12(1) - c1(5)*MP12(1)**2 -  &
     &    c1(6)*MP12(1)**2 - c1(7)*MP12(1)**2 -  &
     &    c1(8)*MP12(1)**2 + c1(10)*MP12(1)**2 +  &
     &    c1(11)*MP12(1)**2 + c1(12)*MP12(1)**2 +  &
     &    c1(13)*MP12(1)**2)/MP12(1)**2 
!~~~~
elseif(abs(G0) .lt. 1.0e-1_ki) then
! print *, 'mgetc1: using the new version with G0 < 0.1'  
      c1(0) =  &
     &f1(0) 
      c1(1) =  &
     &        (G0*f1(1) + G0**5*f1(2) + G0**3*(f1(0) - f1(3)) - f1(4))/ &
     &  ((-one + G0**6)*MP12(1)) 
      c1(2) =  &
     &        -((f1(1) + G0**2*(f1(0) - f1(3) +  &
     &         G0**2*(f1(2) - G0*f1(4))))/((-one + G0**6)*MP12(1))) 
      c1(5) =  &
     &        (f1(0) - f1(3) + G0**2*(G0**2*f1(1) + f1(2) - G0*f1(4)))/ &
     &  ((-one + G0**6)*MP12(1)**2) 
      c1(6) =  &
     &        (-f1(2) + G0*(-(G0*f1(1)) + G0**3*(-f1(0) + f1(3)) +  &
     &       f1(4)))/((-one + G0**6)*MP12(1)**2) 
      c1(3) =  &
     &        (G0*f1(6) + G0**5*f1(7) + G0**3*(f1(5) - f1(8)) + f1(9))/ &
     &  ((-one + G0**6)*MP12(1)) 
      c1(4) =  &
     &        (f1(6) + G0**2*(f1(5) - f1(8) +  &
     &       G0**2*(f1(7) + G0*f1(9))))/((-one + G0**6)*MP12(1)) 
      c1(7) =  &
     &        (f1(5) - f1(8) + G0**2*(f1(7) + G0*(G0*f1(6) + f1(9))))/ &
     &  ((-one + G0**6)*MP12(1)**2) 
      c1(8) =  &
     &        -((f1(7) + G0*(G0*f1(6) + G0**3*(f1(5) - f1(8)) + f1(9)))/ &
     &    ((-one + G0**6)*MP12(1)**2)) 
      c1(15) =  &
     &(c1(0) - f1(5))/(G0*MP12(1)**2) 
      c1(10) =  &
     &        ((-one + G0)*c1(0) - G0*f1(12) + f1(13) +  &
     &    (-one + G0)*((one + G0)*c1(1) - c1(3))*MP12(1) +  &
     &    (-one + G0)*((one + G0 + G0**2)*c1(5) - G0*c1(6) + c1(7))* &
     &     MP12(1)**2)/((-one + G0**2)*MP12(1)**2) 
      c1(11) =  &
     &        ((-one + G0)*c1(0) - G0*f1(10) + f1(11) +  &
     &    (-one + G0)*((one + G0)*c1(1) - c1(4))*MP12(1) +  &
     &    (-one + G0)*((one + G0 + G0**2)*c1(5) - G0*c1(6) + c1(8))* &
     &     MP12(1)**2)/((-one + G0**2)*MP12(1)**2) 
      c1(12) =  &
     &        ((-one + G0)*c1(0) + f1(12) - G0*f1(13) +  &
     &    (-one + G0)*(c1(2) + G0*c1(2) - c1(3))*MP12(1) +  &
     &    (-one + G0)*(c1(6) + G0*(-c1(5) + c1(6) + G0*c1(6)) +  &
     &       c1(7))*MP12(1)**2)/((-one + G0**2)*MP12(1)**2) 
      c1(13) =  &
     &        ((-one + G0)*c1(0) + f1(10) - G0*f1(11) +  &
     &    (-one + G0)*(c1(2) + G0*c1(2) - c1(4))*MP12(1) +  &
     &    (-one + G0)*(c1(6) + G0*(-c1(5) + c1(6) + G0*c1(6)) +  &
     &       c1(8))*MP12(1)**2)/((-one + G0**2)*MP12(1)**2) 
      c1(14) =  &
     &        (-c1(0) + f1(14) - G0mu*c1(1)*MP12(1) - c1(2)*MP12(1) -  &
     &    G0mu**2*c1(5)*MP12(1)**2 - c1(6)*MP12(1)**2)/mu2g(1) 
!~~~~
else
!print *, 'mgetc1.f90, coefficients G0 not 0'  
      c1(0) =  &
     &f1(0) 
      c1(1) =  &
     &f1(4)/(G0*MP12(1)) 
      c1(2) =  &
     &f1(1)/MP12(1) 
      c1(5) =  &
     &f1(3)/(G0**2*MP12(1)**2) 
      c1(6) =  &
     &f1(2)/MP12(1)**2 
      c1(3) =  &
     &f1(9)/(G0*MP12(1)) 
      c1(4) =  &
     &-(f1(6)/MP12(1)) 
      c1(7) =  &
     &f1(8)/(G0**2*MP12(1)**2) 
      c1(8) =  &
     &f1(7)/MP12(1)**2 
      c1(15) =  &
     &(c1(0) - f1(5))/(G0*MP12(1)**2)
      c1(10) = &
     &        (four*c1(0) - eight*f1(12) + four*f1(13) + &
     &    MP12(1)*(six*G0*c1(1) - four*c1(3) + &
     &       (seven*G0**2*c1(5) - eight*c1(6) + four*c1(7))*MP12(1)))/ &
     &  (six*G0*MP12(1)**2)
      c1(11) = &
     &        (four*c1(0) - eight*f1(10) + four*f1(11) + &
     &    MP12(1)*(six*G0*c1(1) - four*c1(4) + &
     &       (seven*G0**2*c1(5) - eight*c1(6) + four*c1(8))*MP12(1)))/&
     &  (six*G0*MP12(1)**2)
      c1(12) = &
     &        (two*c1(0) + two*f1(12) - four*f1(13) + &
     &    MP12(1)*(six*c1(2) - two*c1(3) + &
     &       (-(G0**2*c1(5)) + 14.0_ki*c1(6) + two*c1(7))*MP12(1)))/&
     &  (six*MP12(1)**2)
      c1(13) = &
     &        (two*c1(0) + two*f1(10) - four*f1(11) + &
     &    MP12(1)*(six*c1(2) - two*c1(4) + &
     &       (-(G0**2*c1(5)) + 14.0_ki*c1(6) + two*c1(8))*MP12(1)))/&
     &  (six*MP12(1)**2)
      c1(14) =  &
     &      (-c1(0) + f1(14) - G0mu*c1(1)*MP12(1) - c1(2)*MP12(1) -  &
     &  G0mu**2*c1(5)*MP12(1)**2 - c1(6)*MP12(1)**2)/mu2g(1) 
endif  
! DEFAULT==========================================================================!
	case default
		print *, 'error in mgetc1.f90'
	end select
		c1(9)=czip ! There is no c1(9)
	if (diff.ge.0) then
		do i=5,15
		c1(i)=czip
		enddo
	if (diff.ge.1) then
		do i=1,4
		c1(i)=czip
		enddo
	endif             
	endif
end subroutine getc1_rm

subroutine DFT1(n,nold,known,f1)
! PARAMETERS=======================================================================!
	! external paramters
		complex(ki),intent(in   ) :: known(15)
		integer,    intent(in   ) :: n
		integer,    intent(inout) :: nold
		complex(ki),intent(inout) :: f1(0:14)
	! internal parameters
		complex(ki) :: summation
		real(ki)    :: theta
		integer     :: k,i
! DOING THE DFT N TIMES AND FILLING F2==============================================!
	theta=twopi/n
	!loop over number of number of f1's to fill
	do i=1,n
		summation=czip
		! DFT summation
		do k=0,n-1
			summation = summation + known(k+nold+1)*&
			&(cos(theta*real((i-1)*k,ki))+im*sin(theta*real((i-1)*k,ki)))
		enddo
		f1(i+nold-1) = summation/n
	enddo
	nold=nold+n		
end subroutine DFT1

subroutine HResL1(resin,nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	use mfunctions
	use precision
	use mglobal, only: denst,resit,mu2t
! PARAMETERS============================================================================!
	! external parameters
		complex(ki), dimension(15),	    intent(inout) :: resin


		integer, 			    intent(in   ) :: nsol, nloops, nleg, cut1
		complex(ki), dimension(15,4), 	    intent(in   ) :: q1
		complex(ki), dimension(15),	    intent(in   ) :: mu2vec
		complex(ki), dimension(4), 	    intent(in   ) :: qt
		real(ki),    dimension(0:nleg-1,4), intent(in   ) :: Vi
		real(ki),    dimension(0:nleg-1),   intent(in   ) :: msq
	! internal parameters
		integer			   :: i, j1, n
		integer			   :: i1, i2, i3, i4, i5
		logical			   :: evalres
		complex(ki), dimension(15) :: densn
		complex(ki) 		   :: densnt
		logical 		   :: notcut
		integer,     dimension(5)  :: iarr
		integer 		   :: dicutn
! INITIALIZATION========================================================================!
	! cuts
		j1 = cut1
		dicutn = 1
! LOOP==================================================================================!
	select case(nloops)
	case(5)
	do i5=4,nleg-1
	do i4=3,i5-1
	do i3=2,i4-1
	do i2=1,i3-1
	do i1=0,i2-1
		iarr(5)=i5
		iarr(4)=i4
		iarr(3)=i3
		iarr(2)=i2
		iarr(1)=i1
	! init
		densn(:) =  cone
		densnt	 =  cone
		evalres  = .false.

	l5: do i=0, nleg-1
		notcut   = .true.
		do n=1,nloops
			if(i.eq.iarr(n)) notcut = .false.
		enddo
		if(notcut) then
			if (i.ne.j1) then
				do n=1,nsol
				densn(n) = densn(n)*denevalmu2(nleg,i,q1(n,:),Vi,msq,mu2vec(n))
				enddo
				densnt   = densnt*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
				evalres  = .true.
			else
				densn(:) = czip
				densnt   = czip
				evalres  = .false.
				exit l5
			endif
		endif
	enddo l5
	if (evalres) then
		do n=1,nsol
			resin(n)=resin(n)+densn(n)*res5(dicutn,mu2vec(n))
		enddo
		resit(1)=resit(1)+densnt*Res5(dicutn,mu2t(1))
	endif
	dicutn = dicutn+1
	enddo
	enddo
	enddo
	enddo
	enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	case(4)
	do i4=3,nleg-1
	do i3=2,i4-1
	do i2=1,i3-1
	do i1=0,i2-1
		iarr(4)=i4
		iarr(3)=i3
		iarr(2)=i2
		iarr(1)=i1
	! init
		densn(:) =  cone
		densnt	 =  cone
		evalres  = .false.

	l4: do i=0, nleg-1
		notcut   = .true.
		do n=1,nloops
			if(i.eq.iarr(n)) notcut = .false.
		enddo
		if(notcut) then
			if (i.ne.j1) then
				do n=1,nsol
				densn(n) = densn(n)*denevalmu2(nleg,i,q1(n,:),Vi,msq,mu2vec(n))
				enddo
				densnt   = densnt*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
				evalres  = .true.
			else
				densn(:) = czip
				densnt   = czip
				evalres  = .false.
				exit l4
			endif
		endif
	enddo l4
	if (evalres) then
		do n=1,nsol
			resin(n)=resin(n)+densn(n)*Res4(dicutn,q1(n,:),mu2vec(n))
		enddo
		resit(1)=resit(1)+densnt*Res4(dicutn,qt,mu2t(1))	
	endif
	dicutn = dicutn+1
	enddo
	enddo
	enddo
	enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	case(3)
	do i3=2,nleg-1
	do i2=1,i3-1
	do i1=0,i2-1
		iarr(3)=i3
		iarr(2)=i2
		iarr(1)=i1
	! init
		densn(:) =  cone
		densnt	 =  cone
		evalres  = .false.

	l3: do i=0, nleg-1
		notcut   = .true.
		do n=1,nloops
			if(i.eq.iarr(n)) notcut = .false.
		enddo
		if(notcut) then
			if (i.ne.j1) then
				do n=1,nsol
				densn(n) = densn(n)*denevalmu2(nleg,i,q1(n,:),Vi,msq,mu2vec(n))
				enddo
				densnt   = densnt*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
				evalres  = .true.
			else
				densn(:) = czip
				densnt   = czip
				evalres  = .false.
				exit l3
			endif
		endif
	enddo l3
	if (evalres) then
		do n=1,nsol
			resin(n)=resin(n)+densn(n)*Res3(dicutn,q1(n,:),mu2vec(n))
		enddo
		resit(1)=resit(1)+densnt*Res3(dicutn,qt,mu2t(1))	
	endif
	dicutn = dicutn+1
	enddo
	enddo
	enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	case(2)
	do i2=1,nleg-1
	do i1=0,i2-1
		iarr(2)=i2
		iarr(1)=i1
	! init
		densn(:) =  cone
		densnt	 =  cone
		evalres  = .false.

	l2: do i=0, nleg-1
		notcut   = .true.
		do n=1,nloops
			if(i.eq.iarr(n)) notcut = .false.
		enddo
		if(notcut) then
			if (i.ne.j1) then
				do n=1,nsol
				densn(n) = densn(n)*denevalmu2(nleg,i,q1(n,:),Vi,msq,mu2vec(n))
				enddo
				densnt   = densnt*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
				evalres  = .true.
			else
				densn(:) = czip
				densnt   = czip
				evalres  = .false.
				exit l2
			endif
		endif
	enddo l2
	if (evalres) then
		do n=1,nsol
			resin(n)=resin(n)+densn(n)*Res2(dicutn,q1(n,:),mu2vec(n))
		enddo
		resit(1)=resit(1)+densnt*Res2(dicutn,qt,mu2t(1))	
	endif
	dicutn = dicutn+1
	enddo
	enddo
	end select
!========================================
end subroutine HResL1


end module mgetc1

