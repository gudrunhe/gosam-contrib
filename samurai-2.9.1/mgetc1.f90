module mgetc1
   use precision, only: ki
   use constants
   use options
   use mfunctions
   use mrestore
   implicit none

   public :: getc1

contains

subroutine getc1(numeval,nleg,rank,c1,cut1,q1,qt,Vi,msq)
	use mglobal, only: G0, mu2g, MP12, mu2t, resit, denst, mu2test, G0mu
	implicit none
! PARAMETERS===========================================================================================!
	! external
	integer,			    intent(in ) :: nleg, rank, cut1
	complex(ki), dimension(15,4),       intent(in ) :: q1
	complex(ki), dimension(4),          intent(in ) :: qt
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	complex(ki),    dimension(0:nleg-1),   intent(in ) :: msq
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
			resi4(n)=Res4(diff,1,q1(n,:),mu2vec(n))
		enddo
		resit(1)=Res4(diff,1,qt,mu2t(1))
	elseif (nleg.eq.3) then
		do n=1,nsol
			resi3(n)=Res3(diff,1,q1(n,:),mu2vec(n))
		enddo
		resit(1)=Res3(diff,1,qt,mu2t(1))
	elseif (nleg.eq.2) then
		do n=1,nsol
			resi2(n)=Res2(diff,1,q1(n,:),mu2vec(n))
		enddo
		resit(1)=Res2(diff,1,qt,mu2t(1))
	else
! Pentuple cuts contribution---------------------------------------------------------------------------!
		nloops = 5
		call HResL1(diff,resi5,	nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	endif
! Quadrupole cuts contribution-------------------------------------------------------------------------!
	if((nleg .ne. 4) .and. (nleg .ne. 3) .and. (nleg .ne. 2)) then
			nloops = 4
			call HResL1(diff,resi4,	nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	endif
! Triple cuts contribution-----------------------------------------------------------------------------!
	if((nleg .ne. 3) .and. (nleg .ne. 2)) then
			nloops = 3
			call HResL1(diff,resi3,	nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	endif
! Double cuts contribution ----------------------------------------------------------------------------!
	if(nleg .ne. 2) then
			nloops = 2
			call HResL1(diff,resi2,	nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
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
end subroutine getc1

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

subroutine HResL1(diff,resin,nsol,nloops,nleg,cut1,q1,qt,Vi,msq,mu2vec)
	use mfunctions
	use precision
	use mglobal, only: denst,resit,mu2t
! PARAMETERS============================================================================!
	! external parameters
		integer, intent(in) :: diff
		complex(ki), dimension(15),	    intent(inout) :: resin


		integer, 			    intent(in   ) :: nsol, nloops, nleg, cut1
		complex(ki), dimension(15,4), 	    intent(in   ) :: q1
		complex(ki), dimension(15),	    intent(in   ) :: mu2vec
		complex(ki), dimension(4), 	    intent(in   ) :: qt
		real(ki),    dimension(0:nleg-1,4), intent(in   ) :: Vi
		complex(ki),    dimension(0:nleg-1),   intent(in   ) :: msq
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
			resin(n)=resin(n)+densn(n)*Res4(diff,dicutn,q1(n,:),mu2vec(n))
		enddo
		resit(1)=resit(1)+densnt*Res4(diff,dicutn,qt,mu2t(1))	
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
			resin(n)=resin(n)+densn(n)*Res3(diff,dicutn,q1(n,:),mu2vec(n))
		enddo
		resit(1)=resit(1)+densnt*Res3(diff,dicutn,qt,mu2t(1))	
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
			resin(n)=resin(n)+densn(n)*Res2(diff,dicutn,q1(n,:),mu2vec(n))
		enddo
		resit(1)=resit(1)+densnt*Res2(diff,dicutn,qt,mu2t(1))	
	endif
	dicutn = dicutn+1
	enddo
	enddo
	end select
!========================================
end subroutine HResL1


end module mgetc1

