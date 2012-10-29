module mgetqs
	use precision
	use constants
	use mfunctions
	implicit none

	private

	interface getq5
		module procedure getq5_rm
		module procedure getq5_cm
	end interface getq5

	interface getq4
		module procedure getq4_rm
		module procedure getq4_cm
	end interface getq4

	interface getq3
		module procedure getq3_rm
		module procedure getq3_cm
	end interface getq3

	interface getq2
		module procedure getq2_rm
		module procedure getq2_cm
	end interface getq2

	interface getq1
		module procedure getq1_rm
		module procedure getq1_cm
	end interface getq1

	public :: getq1, getq2, getq3, getq4, getq5

contains

subroutine getq5_rm(nleg,cut5,e1,e2,e3,e4,p0,Vi,msq,r1,r2,q5,mu2)
	implicit none
!PARAMETERS========================================================================!
	! external
	integer, 			    intent(in ) :: nleg, cut5
	real(ki),    dimension(4), 	    intent(in ) :: e1, e2, p0
	complex(ki), dimension(4), 	    intent(in ) :: e3, e4
	real(ki), 			    intent(in ) :: r1, r2
	real(ki),    dimension(0:nleg-1),   intent(in ) :: msq
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	complex(ki), dimension(4),          intent(out) :: q5
	complex(ki), 		  	    intent(out) :: mu2
	! internal
	integer 	       :: j1, j2, j3, j4, j5
	real(ki) 	       :: x1, x2, var2, var3
	real(ki) 	       :: MP12, MP1v2, MP1v3, MP2v2, MP2v3, MPv22, MPv33
	real(ki), dimension(4) :: v2, v3
	complex(ki) 	       :: MP3v2, MP3v3, MP4v2, MP4v3, x3, x4, tmu2, den
! INITIALIZATION===================================================================!
	! cuts
	j5 =  cut5/10000
	j4 = (cut5-j5*10000)/1000
	j3 = (cut5-j5*10000-j4*1000)/100
	j2 = (cut5-j5*10000-j4*1000-j3*100)/10
	j1 =  cut5-j5*10000-j4*1000-j3*100-j2*10
	! internal momenta
	v2(:) = Vi(j3,:) - p0(:)
	v3(:) = Vi(j4,:) - p0(:)
	! dotproducts
	MP12  = sdot(e1,e2)
	MP1v2 = sdot(v2,e1)
	MP1v3 = sdot(v3,e1)
	MP2v2 = sdot(v2,e2)
	MP2v3 = sdot(v3,e2)
	MP3v2 = sdot(v2,e3)
	MP3v3 = sdot(v3,e3)
	MP4v2 = sdot(v2,e4)
	MP4v3 = sdot(v3,e4)
	MPv22 = sdot(v2,v2)
	MPv33 = sdot(v3,v3)
! CALCULATING XS===================================================================!
	den  = -(MP3v3*MP4v2*two) + MP3v2*MP4v3*two
	x1   = (0.5_ki*(  -(MP12*(1.0_ki+r1)*r2*two)-(1.0_ki+r2)*msq(j1) &
		& +r2*msq(j2) + msq(j5)))/(MP12*(-1.0_ki + r1*r2))
	x2   = -((0.5_ki*(-(MP12*r1*(1.0_ki+r2)*two)-(1.0_ki+r1)*msq(j1) &
		& +msq(j2) + r1*msq(j5)))/(MP12*(-1.0_ki + r1*r2)))
	var2 = MPv22 + MP1v2*two*x1 + MP2v2*two*x2+msq(j1)-msq(j3)
	var3 = MPv33 + MP1v3*two*x1 + MP2v3*two*x2+msq(j1)-msq(j4)
	x3   = (-MP4v3*var2+MP4v2*var3)/den
	x4   = ( MP3v3*var2-MP3v2*var3)/den
	! test values
	tmu2 = two*MP12*(x1*x2-x3*x4)-msq(j1)
	mu2  = real(tmu2, ki)
! SETTING OUTPUT===================================================================!
	q5(:) = -p0(:)+x1*e1(:)+x2*e2(:)+x3*e3(:)+x4*e4(:)
end subroutine getq5_rm


subroutine getq4_rm(nleg,cut4,e1,e2,e3,e4,p0,k1,k2,k3,L3,r1,r2,q4,qt,msq)
	use mglobal, only: MP12, mu2g, mu2t, dx
	implicit none
! PARAMETERS=======================================================================!
	! external
	integer, 			  intent(in ) :: nleg, cut4
	real(ki),    dimension(4), 	  intent(in ) :: k1, k2, k3, L3, e1, e2, p0
	complex(ki), dimension(4),        intent(in ) :: e3, e4
	real(ki),    dimension(0:nleg-1), intent(in ) :: msq
	real(ki), 			  intent(in ) :: r1,r2
	complex(ki), dimension(6,4),	  intent(out) :: q4
	complex(ki), dimension(4),	  intent(out) :: qt
	! internal
	complex(ki), dimension(6) :: x3,x4
	complex(ki) 		  :: x3t, x4t
	complex(ki) 		  :: ME33, ME34, ML33, ML34
	complex(ki) 		  :: B0, B1, B2, rtdel, A3, A4
	real(ki)    		  :: beta, x1, x2, A1, A2
	real(ki)    		  :: ME31, ME32, MK11, MK22, MK33
	integer     		  :: j1, j2, j3, j4, j
! INITIALIZATION===================================================================!
	! cuts
	j4 =  cut4/1000
	j3 = (cut4-j4*1000)/100
	j2 = (cut4-j4*1000-j3*100)/10
	j1 =  cut4-j4*1000-j3*100-j2*10
	! dotproducts
	MP12(4) = sdot(e1,e2)
	ML33    = sdot(L3,e3)
	ML34    = sdot(L3,e4)
	ME31    = sdot(k3,e1)
	ME32    = sdot(k3,e2)
	ME33    = sdot(k3,e3)
	ME34    = sdot(k3,e4)
	MK11    = sdot(k1,k1)
	MK22    = sdot(k2,k2)
	MK33    = sdot(k3,k3)
! SELECTING MU2====================================================================!
	mu2g(4) = max(msq(j1),msq(j2),msq(j3),msq(j4),&
		& abs(MK11),abs(MK22),abs(MK33),&
		& abs(two*sdot(k1,k2)),abs(two*sdot(k2,k3)),abs(two*sdot(k1,k3)),&
		& abs(two*MP12(4)))
	if (abs(mu2g(4)).lt.1.0e-10_ki) mu2g(4)=one
! CALCULATING XS===================================================================!
	! all x1 and x2 are the same
	beta    =  one/(one-r1*r2)
	A1      = (msq(j2)-msq(j1)-MK11)/two/MP12(4)
	A2      = (msq(j1)-msq(j4)+MK22)/two/MP12(4)
	x1      =  beta*(A2-r2*A1)
	x2      =  beta*(A1-r1*A2)
	A3      = (msq(j3)-msq(j1)-MK33-two*x1*ME31-two*x2*ME32)/two/ME33
	A4      = -ME34/ME33
	B1      = -two*MP12(4)*A3
	B2      = -two*MP12(4)*A4
	! x(1) and x(2), mu2=0
	B0      = (two*x1*x2*MP12(4)-msq(j1))*cone
	rtdel   =  sqrt(B1**2-four*B0*B2)
	x4(1)   = (-B1+rtdel)/two/B2
	x4(2)   = (-B1-rtdel)/two/B2
	x3(1)   =  x4(1)*A4+A3
	x3(2)   =  x4(2)*A4+A3
	! x(3) and x(4), mu2=mu2 
	B0      = (two*x1*x2*MP12(4)-msq(j1)-mu2g(4))*cone
	rtdel   =  sqrt(B1**2-four*B0*B2)
	x4(3)   = (-B1+rtdel)/two/B2
	x4(4)   = (-B1-rtdel)/two/B2
	x3(3)   =  x4(3)*A4+A3
	x3(4)   =  x4(4)*A4+A3
	! x(5) and x(6), mu2=-mu2
	B0      = (two*x1*x2*MP12(4)-msq(j1)+mu2g(4))*cone
	rtdel   =  sqrt(B1**2-four*B0*B2)
	x4(5)   = (-B1+rtdel)/two/B2
	x4(6)   = (-B1-rtdel)/two/B2
	x3(5)   =  x4(5)*A4+A3
	x3(6)   =  x4(6)*A4+A3
	! xt
	mu2t(4) =  half
	B0      = (two*x1*x2*MP12(4)-msq(j1)-mu2t(4))*cone
	rtdel   =  sqrt(B1**2-four*B0*B2)
	x4t     = (-B1-rtdel)/two/B2
	x3t     =  x4t*A4+A3
! SETTING OUTPUT===================================================================!
	do j=1,6
		q4(j,:) = -p0(:)+x1*e1(:)+x2*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
	enddo
		qt(:)   = -p0(:)+x1*e1(:)+x2*e2(:)+x3t*e3(:)+x4t*e4(:)
	! dx's (mglobal)
		dx(1) = ML33*x3(1)-ML34*x4(1)
		dx(2) = ML33*x3(2)-ML34*x4(2)
		dx(3) = ML33*x3(3)-ML34*x4(3)
		dx(4) = ML33*x3(4)-ML34*x4(4)
		dx(5) = ML33*x3(5)-ML34*x4(5)
end subroutine getq4_rm



subroutine getq3_rm(nleg,irank,cut3,e1,e2,e3,e4,p0,k1,k2,msq,r1,r2,q3,qt)
	use mglobal, only: C0,C1,Cm1,mu2g,MP12,KK,Kmu,mu2t
	use options, only: C0_thrs
	implicit none
! PARAMETERS=======================================================================!
	! external
	integer, 			  intent(in ) :: nleg, cut3, irank
	real(ki),    dimension(4), 	  intent(in ) :: p0, k1, k2, e1, e2
	complex(ki), dimension(4), 	  intent(in ) :: e3, e4
	real(ki),    dimension(0:nleg-1), intent(in ) :: msq
	complex(ki), dimension(15,4), 	  intent(out) :: q3
	complex(ki), dimension(4), 	  intent(out) :: qt
	! internal
	integer 		   :: j1,j2,j3,j,ndiff
	real(ki) 		   :: r1,r2,beta,teta,x1,x2,A1,A2,C1t
	real(ki) 		   :: MK11,MK12,MK22
	complex(ki) 		   :: x3t, x4t
	complex(ki), dimension(15) :: x3, x4
! INITIALIZATION===================================================================!
	ndiff  = nleg-irank
	KK(3)  = one
	Kmu(3) = one
	! cuts
	j3 =  cut3/100
	j2 = (cut3-j3*100)/10
	j1 =  cut3-j3*100-j2*10
	! dotproducts
	MP12(3) = sdot(e1,e2)
	MK11    = sdot(k1,k1)
	MK12    = sdot(k1,k2)
	MK22    = sdot(k2,k2)
	! mu2
	mu2g(3) = max(msq(j1),msq(j2),msq(j3),abs(MK11),abs(MK22),abs(two*MK12))
	if (abs(mu2g(3)).lt.1d-10) mu2g(3)=one

! CALCULATING XS===================================================================!
	! x1 and x2
	beta = 1d0/(1d0-r1*r2)
	A1   = (msq(j2)-msq(j1)-MK11)/two/MP12(3)
	A2   = (msq(j1)-msq(j3)+MK22)/two/MP12(3)
	x1   = beta*(A2-r2*A1)
	x2   = beta*(A1-r1*A2)
	! C0 and C1
	C0   = x1*x2-msq(j1)/two/MP12(3)
	C1   = x1*x2-(msq(j1)+mu2g(3))/two/MP12(3)
! C0=1 (safe around C0=1)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	if (abs(C0-1.0_ki) .lt. C0_thrs) then
	select case(ndiff)
	case(-1)	
		teta = twopi/nine
		do j=1,9
			x3(j)=cos(teta*real(j-1,ki))-im*sin(teta*real(j-1,ki))
			x4(j)=C0*cone/x3(j)
		enddo
		teta = twopi/five
		do j=10,14
			x4(j)=cos(teta*real(j-10,ki))+im*sin(teta*real(j-10,ki))
			x3(j)=C1*cone/x4(j)
		enddo
		Cm1 = x1*x2-(msq(j1)-mu2g(3))/two/MP12(3) !-mu2 instead of +mu2
			x3(15) = cone
			x4(15) = Cm1*cone
		do j=1,15
			q3(j,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
		enddo
	case default
		teta=twopi/seven
		do j=1,7
			x3(j)=cos(teta*real(j-1,ki))-im*sin(teta*real(j-1,ki))
			x4(j)=C0*cone/x3(j)
		enddo
                teta=twopi/three
		do j=8,10
			x4(j)=cos(teta*real(j-8,ki))+im*sin(teta*real(j-8,ki))
			x3(j)=C1*cone/x4(j)
		enddo
		do j=1,10
			q3(j,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
		enddo
	end select
! C0 unequal 1 (safe around C0=0)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	else
	! Branching according to rank
	select case(ndiff)
	case(3)
		! rank1: 1 coefficient, c3(0)--------------------------------------!
		x3(1) = cone
		x4(1) = C0	
		q3(1,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(1)*e3(:)+x4(1)*e4(:)
		x3(2:15) = czip
	case(1) !------------------------------------------------------------------!
		! rank2 c-system: 6  coefficients
		teta = twopi/three
		do j=1,3
			x3(j)=cos(teta*real(j-1,ki))-im*sin(teta*real(j-1,ki))
			x4(j)=C0*cone/x3(j)
		enddo
		teta = twopi/two
		do j=4,5
			x3(j)=cos(teta*real(j-4,ki))+im*sin(teta*real(j-4,ki))
			x4(j)=cone/x3(j)
		enddo 
		teta = twopi
		do j=6,6
			x3(j)=cos(teta*real(j-6,ki))+im*sin(teta*real(j-6,ki))
			x4(j)=cone/x3(j)
		enddo
			x3(7:15) = czip
			x4(7:15) = czip 
		! here are the qs for rank 2
		q3(1,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(1)*e3(:)+x4(1)*e4(:)
		q3(2,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(2)*e3(:)+x4(2)*e4(:)
		q3(3,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(3)*e3(:)+x4(3)*e4(:)
		q3(4,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C0*x3(4)*e3(:)+x4(4)*e4(:)
		q3(5,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C0*x3(5)*e3(:)+x4(5)*e4(:)
		q3(6,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C1*x3(6)*e3(:)+x4(6)*e4(:)
	case(2) !------------------------------------------------------------------!
		! rank1 c-system: 3 coefficients      
			x3(1)    =  cone
			x4(1)    =  C0*cone
			x3(2)    = -cone
			x4(2)    = -C0*cone         
			x3(3)    =  C0*cone
			x4(3)    =  cone
			x3(4:15) =  czip
			x4(4:15) =  czip
		do j=1,3
			q3(j,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
		enddo
	case(0) !------------------------------------------------------------------!
		! 10 coefficients
		teta=twopi/four
		do j=1,4
			x3(j)=cos(teta*real(j-1,ki))-im*sin(teta*real(j-1,ki))
			x4(j)=C0*cone/x3(j)
		enddo
		teta=twopi/three
		do j=5,7
			x4(j)=cos(teta*real(j-5,ki))+im*sin(teta*real(j-5,ki))
			x3(j)=C0*cone/x4(j)
		enddo
		teta=twopi/three
		do j=8,10
			x4(j)=cos(teta*real(j-8,ki))+im*sin(teta*real(j-8,ki))
			x3(j)=C1*cone/x4(j)
		enddo
		do j=1,10
			q3(j,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
		enddo
	case(-1) !-----------------------------------------------------------------!
		! complete c-system: 15 coefficients
		teta=twopi/five
		do j=1,5
			x3(j)=cos(teta*real(j-1,ki))-im*sin(teta*real(j-1,ki))
			x4(j)=C0*cone/x3(j)
		enddo
		teta=twopi/four
		do j=6,9
			x4(j)=cos(teta*real(j-6,ki))+im*sin(teta*real(j-6,ki))
			x3(j)=C0*cone/x4(j)
		enddo        
		teta=twopi/five
		do j=10,14
			x4(j)=cos(teta*real(j-10,ki))+im*sin(teta*real(j-10,ki))
			x3(j)=C1*cone/x4(j)
		enddo
		Cm1=x1*x2-(msq(j1)-mu2g(3))/two/MP12(3) !-mu2 instead of +mu2
			x4(15) = one
			x3(15) = Cm1/x4(15)
		do j=1,15
			q3(j,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
		enddo
	case default !-------------------------------------------------------------!
		print *, 'error in mgetqs.f90, getq3_rm'
	end select !---------------------------------------------------------------!
	end if
	! test values   
	mu2t(3)=half
	C1t=x1*x2-(msq(j1)+mu2t(3))/two/MP12(3)
	x3t=32.4_ki
	x4t=C1t*cone/x3t
	qt(:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3t*e3(:)+x4t*e4(:)
end subroutine getq3_rm



subroutine getq2_rm(nleg,irank,cut2,e1,e2,e3,e4,p0,k1,msq,q2,qt)
	use mglobal, only: mu2g,MP12,mu2t, y1,y2,y3,x1s1,x1s2, Fz,Fp,Fm,F1z,KK,KB,Kmu, FL2
	implicit none
! PARAMETERS===================================================================================!
	!external parameters
	integer, 			  intent( in) :: nleg, cut2, irank
	real(ki),    dimension(4), 	  intent( in) :: p0, k1, e1, e2
	complex(ki), dimension(4), 	  intent( in) :: e3, e4
	real(ki),    dimension(0:nleg-1), intent( in) :: msq
	complex(ki), dimension(20,4), 	  intent(out) :: q2
	complex(ki), dimension(4), 	  intent(out) :: qt
	!internal parameters
	real(ki)		   :: K1SQ,x1t,x2t,Ft
	complex(ki)		   :: x3t,x4t
	real(ki),    dimension(20) :: x1,x2
	complex(ki), dimension(20) :: x3,x4
	integer 		   :: j1,j2,j, ndiff
	real(ki) 		   :: A2,B2,D2,Ar,Br,Dr
	complex(ki)		   :: C2,C2mu,Cr
	integer 		   :: Fstep, nold
	logical 		   :: mu2zero,mu2nonzero,samplex3,samplex4
	real(ki)		   :: X2z,X2p,X2m,teta
! INITIALIZATION===============================================================================!
	j2 = cut2/10
	j1 = cut2-j2*10
	ndiff   = nleg - irank
	K1SQ    = sdot(k1,k1)
	MP12(2) = sdot(e1,e2)
	! Selecting mu2 dynamically
	mu2g(2) = max(msq(j1),msq(j2),K1SQ)
	if (abs(mu2g(2)) .lt. 1.0e-10_ki) mu2g(2)=one
	! Defining some boolians to increase readibility:
		mu2zero    = .false.
		mu2nonzero = .true.
		samplex3   = .true.
		samplex4   = .false.
	! Initialize to zero
		x1s1 = zip
		x1s2 = zip
	! For rank 1 and 2
		KK(2)  = one
		KB     = one
		Kmu(2) = one
		X2z = (msq(j2)-msq(j1)-K1SQ)/two/MP12(2)
		Fz  = -msq(j1)/two/MP12(2)
		X2p =  X2z-one*K1SQ/two/MP12(2)
		Fp  =  one*X2p-msq(j1)/two/MP12(2)
		X2m =  X2z+one*K1SQ/two/MP12(2)
		Fm  = -one*X2m-msq(j1)/two/MP12(2)
		F1z = -(msq(j1)+mu2g(2))/two/MP12(2)
	select case(ndiff)
! RANK 0=======================================================================================!
	case(2)
		x1(1)=zip
		x2(1)=X2z
		! Changed this into x3<->x4 to compare with old version case(2) equivalent
		x4(1)=cone
		x3(1)=Fz
		q2(1,:)=-p0(:)+x1(1)*e1(:)+x2(1)*e2(:)+x3(1)*e3(:)+x4(1)*e4(:)
! RANK 1=======================================================================================!
	case(1)
	teta=twopi/two
	do j=1,2
		x1(j)=zip
		x2(j)=X2z
		x4(j)=cos(teta*real(j-1,ki))+im*sin(teta*real(j-1,ki))
		x3(j)=Fz/x4(j)
	enddo
	do j=3,4
		x1(j)=one
		x2(j)=X2p
		x3(j)=cos(teta*real(j-3,ki))-im*sin(teta*real(j-3,ki))
		x4(j)=Fp/x3(j)
	enddo
	do j=5,10
		x3(j)=czip
		x4(j)=czip
	enddo
	do j=1,4
		q2(j,:)=-p0(:)+x1(j)*e1(:)+x2(j)*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
	enddo
! RANK 2=======================================================================================!
	case(0)
	teta=twopi/three
	do j=1,3
		x1(j)=zip
		x2(j)=X2z
		x4(j)=cos(teta*real(j-1,ki))+im*sin(teta*real(j-1,ki))
		x3(j)=Fz/x4(j)
	enddo
	teta=twopi/two
	do j=4,5
		x1(j)=zip
		x2(j)=X2z
		x3(j)=cos(teta*real(j-4,ki))-im*sin(teta*real(j-4,ki))
		x4(j)=Fz/x3(j)
	enddo
	teta=twopi/two
	do j=6,7
		x1(j)=one
		x2(j)=X2p
		x4(j)=cos(teta*real(j-6,ki))+im*sin(teta*real(j-6,ki))
		x3(j)=Fp/x4(j)
	enddo
		x1(8)=one
		x2(8)=X2p
		x3(8)=cone
		x4(8)=Fp*cone

		x1(9)=-one
		x2(9)=X2m
		x3(9)=Fm*cone
		x4(9)=cone

		x1(10)=zip
		x2(10)=X2z
		x3(10)=F1z*cone
		x4(10)=cone
	do j=1,10
		q2(j,:)=-p0(:)+x1(j)*e1(:)+x2(j)*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
	enddo
! RANK 3=======================================================================================!
	case(-1)
	! Level 2 has in general a second degree polynomial in x1:
	! x3*x4 = F(x1) = A2*x1^2 + B2*x1 + C2 with:
		A2 = -K1SQ/two/MP12(2)
		B2 = (msq(j2)-msq(j1)-K1SQ)/two/MP12(2)
		C2 = -msq(j1)/two/MP12(2) 		!-(msq(j1)+mu2)/(2*MP12) in general
		C2mu = -(msq(j1)+mu2g(2))/two/MP12(2)
	! Round the values close to zero:
		Ar = A2; Br = B2; Cr = C2
		if(abs(A2).lt. 1.0e-10_ki) Ar=zip
		if(abs(B2).lt. 1.0e-10_ki) Br=zip
		if(abs(C2).lt. 1.0e-10_ki) Cr=zip
	! Calculate the (always rounded) discriminant:
		D2 = B2**2 - 4*A2*C2
		Dr = D2
		if(abs(Dr).lt. 1.0e-10_ki) Dr=zip
	! Doing the DFT for level 2, defined as:
	! SampleL2(mu2 nonzero, x1, number of samplings, sampling x3 (vs. x4), standard io)
	! nold is an integer used to count how many times SampleL2 was called before
		nold = 0
	! Three sampling values for x1: y1, y2, y3. 
	! Notice that the following values are forbidden:
	! y1 = y2, y1^2=y2^2, y=0, y^2=0, F(0,y)=0, Same for y2->y3
		y1 = one/sqrt(seven)
		y2 = one/sqrt(five)
		y3 = one/sqrt(three)
	! Sampling is in six steps, labeled by Fstep:
	! x1=0; x1=x1s1 or y1; x1=x1s2 or y2; x1=y3; mu2=mu2,x1=0; mu2=mu2,x1=1
	!---------------------------------------------------------------------------------------!
	Fstep = 1 !First 7 samplings: sampling at x1=0: F(0,0)                                 	!
			call FillFL2(zip,	A2,B2,C2,Fstep)                                	!
		if((Cr==0) .or. (abs(FL2(1)) .lt. 0.1_ki)) then 				!
			! F(0,0) zero or ~0 (x1=0 is a solution)        			!
			call SampleL2(mu2zero,zip, 4,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
			call SampleL2(mu2zero,zip, 3,samplex4,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
		else 	! F(0,0) nonzero (x1=0 is not a solution)                               !
			call SampleL2(mu2zero,zip, 7,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
		endif                                                                           !
	!---------------------------------------------------------------------------------------!
	Fstep = 2 !Next 5 samplings: sampling if possible at the 1st solution x1s1: F(0,x1s1)=0 !
		if (Ar==0 .and. Br==0 .and. Cr==0) then						!
			!print *, 'mgetqs.f90: Ar=Br=Cr=0, x1s1=y1'				!
			x1s1 = y1								!
			call FillFL2(y1,	A2,B2,C2,Fstep)					!
			call SampleL2(mu2zero,y1,3,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
			call SampleL2(mu2zero,y1,2,samplex4,	nold,x1,x2,x3,x4,Fstep,A2,B2)	!
		elseif(	(Ar==0 .and. Br==0).or.(Ar==0 .and. Cr==0) &				!
		&	.or.(Br==0 .and. Cr==0) .or. (Dr .lt. 0)) then				!
			!F(0,x1)=0 has no non-zero solution: Use y1                             !
			call FillFL2(y1,	A2,B2,C2,Fstep)                                 !
			call SampleL2(mu2zero,y1,  5,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
		else 	                                                                        !
			!F(0,x1)=0 has at least one non-zero solution: x1s1                     !
			if(Ar==0) then								!
				x1s1 = -C2/B2				               		!
			elseif(Cr==0) then							!
				x1s1 = -B2/A2							!
			else									!
			!print *, 'Fstep2 else,else (A,C nonzero)'				!
				x1s1 = (-B2-sign(1.0_ki,B2)*sqrt(D2))/two/A2			!
			endif									!
			 !print *, 'mgetqs.f90: using 3+2 x1s1=',x1s1				!
			call FillFL2(x1s1,	A2,B2,C2,Fstep)                                 !
			call SampleL2(mu2zero,x1s1,3,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
			call SampleL2(mu2zero,x1s1,2,samplex4,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
		endif                                                                           !
	!---------------------------------------------------------------------------------------!	
	Fstep = 3 !Next 3 samplings: sampling if possible at the 2nd solution x1s2: F(0,x1s2)=0 !
		if (Ar==0 .and. Br==0 .and. Cr==0) then						!
			x1s2=y2									!
			call FillFL2(x1s2,	A2,B2,C2,Fstep)                                 !
			call SampleL2(mu2zero,x1s2,2,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
			call SampleL2(mu2zero,x1s2,1,samplex4,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
		elseif(Ar==0 .or. Cr==0 .or. (Dr .le. 0) .or. (abs(Cr)/abs(Ar) .lt. 1e-20)) then!
			! F(0,x1)=0 has no second solution: Use y2                              !
			!print *, 'mgetqs.f90: Using y2'					!
			call FillFL2(y2,	A2,B2,C2,Fstep)                                 !
			call SampleL2(mu2zero,y2,  3,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
		else 	                                                                        !
			! F(0,x1)=0 has a second solution: x1s2                                 !
			x1s2 = (-B2+sign(1.0_ki,B2)*sqrt(D2))/two/A2				!
			call FillFL2(x1s2,	A2,B2,C2,Fstep)                                 !
			call SampleL2(mu2zero,x1s2,2,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
			call SampleL2(mu2zero,x1s2,1,samplex4,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
		endif                                                                           !
	!---------------------------------------------------------------------------------------!
	Fstep = 4 !Next sampling: sampling at a value that is not a solution: F(0,y3)=nonzero   !
			call FillFL2(y3,	A2,B2,C2,Fstep)                                 !
			call SampleL2(mu2zero,y3,  1,samplex3,	nold,x1,x2,x3,x4,Fstep,A2,B2)   !
	!---------------------------------------------------------------------------------------!
	Fstep = 5 !Next 3 samplings: sampling at mu2 nonzero, x1=0! 			        !
			call FillFL2(zip,	A2,B2,C2mu,Fstep)			        !
		if(abs(C2mu) .lt. 1.0e-10_ki) then					        !
			!F(mu2,0)=0                         				        !
			call SampleL2(mu2nonzero,zip,2,samplex3, nold,x1,x2,x3,x4,Fstep,A2,B2)  !
			call SampleL2(mu2nonzero,zip,1,samplex4, nold,x1,x2,x3,x4,Fstep,A2,B2)  !
		else 	                                                                        !
			!F(mu2,0)=nonzero				                        !
			call SampleL2(mu2nonzero,zip,3,samplex3, nold,x1,x2,x3,x4,Fstep,A2,B2)  !
		endif                                                                           !
	!---------------------------------------------------------------------------------------!
	Fstep = 6 !Last sampling: sampling at mu2 and x1 nonzero			        !
			call FillFL2(one,	A2,B2,C2mu,Fstep)                               !
			call SampleL2(mu2nonzero,one,1,samplex3, nold,x1,x2,x3,x4,Fstep,A2,B2)  !
	!---------------------------------------------------------------------------------------!
! SETTING OUTPUT================================================================================!
	!Defining q as Sum_i x_i*e_i
		do j=1,20
		q2(j,:) = -p0(:)+x1(j)*e1(:)+x2(j)*e2(:)+x3(j)*e3(:)+x4(j)*e4(:)
		enddo
! DEFAULT CASE==================================================================================!
	case default
		print *, 'error in mgetqs.f90, getq2_rm'
	end select
! TEST VARIABLE=================================================================================!
		mu2t(2) =  half
		x1t     =  2.3_ki
		B2 = (msq(j2)-msq(j1)-K1SQ)/two/MP12(2)
		x2t     =  B2-x1t*K1SQ/two/MP12(2)
		Ft      =  x1t*x2t-(msq(j1)+mu2t(2))/two/MP12(2)
		x3t     = (2.5_ki,0.3_ki)
		x4t     =  Ft/x3t
		qt(:)   = -p0(:)+x1t*e1(:)+x2t*e2(:)+x3t*e3(:)+x4t*e4(:)
end subroutine getq2_rm



subroutine getq1_rm(nleg,irank,cut1,e1,e2,e3,e4,p0,msq,q1,qt)
	use mglobal, only: G0,KK,mu2g,MP12,mu2t,G0mu
	implicit none
! PARAMETERS===================================================================================!
	! external
	integer,		          intent(in ) :: nleg, irank, cut1
	real(ki),    dimension(0:nleg-1), intent(in ) :: msq
	real(ki),    dimension(4),        intent(in ) :: p0, e1, e2
	complex(ki), dimension(4),        intent(in ) :: e3, e4
	complex(ki), dimension(15,4),     intent(out) :: q1
	complex(ki), dimension(4),        intent(out) :: qt
	! internal
	complex(ki), dimension(15) :: x1, x2 ! note: changed to complex to be able to sample
	complex(ki), dimension(15) :: x3, x4
	real(ki) 		  :: varx,x1t,x2t,x3t,x4t,theta
	integer 		  :: j1, i,nold,n, ndiff
! INITIALIZATION===============================================================================!
	j1      = cut1
	mu2g(1) = msq(j1)
	if (abs(mu2g(1)).lt. 1.0e-10_ki) mu2g(1)=one
	MP12(1) = sdot(e1,e2)
	G0      = msq(j1)/two/MP12(1) !(msq(j1)+mu2g(1))/two/MP12(1)
	nold    = 0
	ndiff   = nleg - irank
	select case(ndiff)
! RANK 0=======================================================================================!
	case(1)
		!---
		mu2g(1)=zip
		if (msq(j1).eq.zip) then
		mu2g(1)=one
		endif
		G0 = (msq(j1)+mu2g(1))/two/MP12(1)
		!---

		x1(1)=G0
		x2(1)=one
		x3(1) =  zip
		x4(1) =  zip
		q1(1,:)=-p0(:)+x1(1)*e1(:)+x2(1)*e2(:)+x3(1)*e3(:)+x4(1)*e4(:)
! RANK 1=======================================================================================!
	case(0)
		!---
		mu2g(1)=zip
		if (msq(j1).eq.zip) then
		mu2g(1)=one
		endif
		G0 = (msq(j1)+mu2g(1))/two/MP12(1)
		!---

		x1(1) =  G0
		x2(1) =  one
		x3(1) =  zip
		x4(1) =  zip

		x1(2) = -G0
		x2(2) = -one
		x3(2) =  zip
		x4(2) =  zip

		x1(3) =  G0
		x2(3) =  one
		x3(3) =  one
		x4(3) =  zip

		x1(4) =  zip
		x2(4) =  one
		x3(4) = -one
		x4(4) =  G0

		x1(5) =  zip
		x2(5) =  one
		x3(5) =  one
		x4(5) = -G0
	do i=1,5
		q1(i,:)=-p0(:)+x1(i)*e1(:)+x2(i)*e2(:)+x3(i)*e3(:)+x4(i)*e4(:)
	enddo
! RANK 2=======================================================================================!
	case(-1)
	if(abs(G0) .lt. 1.0e-10_ki) then
		n = 3 ! 
		theta = twopi/n
		do i=1,n
			x1(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x2(i+nold) = czip
			x3(i+nold) = czip
			x4(i+nold) = czip
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
		n = 2 ! 
		theta = twopi/n
		do i=1,n
			x1(i+nold) = czip
			x2(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x3(i+nold) = czip
			x4(i+nold) = czip
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
		n = 2
		theta = twopi/n
		do i=1,n
			x1(i+nold) = czip
			x2(i+nold) = czip
			x3(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x4(i+nold) = czip
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
		n = 2
		theta = twopi/n
		do i=1,n
			x1(i+nold) = czip
			x2(i+nold) = czip
			x3(i+nold) = czip
			x4(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
			x1(10) = czip
			x2(10) = cone
			x3(10) = czip
			x4(10) = cone

			x1(11) = czip
			x2(11) = cone
			x3(11) = cone
			x4(11) = czip

			x1(12) = cone
			x2(12) = czip
			x3(12) = czip
			x4(12) = cone

			x1(13) = cone
			x2(13) = czip
			x3(13) = cone
			x4(13) = czip

			x1(14) = cone
			x2(14) = cone
			x3(14) = cone
			x4(14) = cone
	elseif(abs(G0) .lt. 1.0e-1_ki) then
		n = 3 ! 
		theta = twopi/n
		do i=1,n
			x1(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x2(i+nold) = G0/x1(i+nold)
			x3(i+nold) = czip
			x4(i+nold) = czip
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
		n = 2 ! 
		theta = twopi/n
		do i=1,n
			x2(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x1(i+nold) = G0/x2(i+nold)
			x3(i+nold) = czip
			x4(i+nold) = czip
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
		n = 3
		theta = twopi/n
		do i=1,n
			x1(i+nold) = czip
			x2(i+nold) = czip
			x3(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x4(i+nold) = -G0/x3(i+nold)
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
		n = 2
		theta = twopi/n
		do i=1,n
			x1(i+nold) = czip
			x2(i+nold) = czip
			x4(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x3(i+nold) = -G0/x4(i+nold)
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
			x1(11) = cone
			x2(11) = G0
			x3(11) = cone
			x4(11) = czip

			x1(12) = G0
			x2(12) = cone
			x3(12) = cone
			x4(12) = czip

			x1(13) = cone
			x2(13) = G0
			x3(13) = czip
			x4(13) = cone

			x1(14) = G0
			x2(14) = cone
			x3(14) = czip
			x4(14) = cone
	else !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
		n = 5 ! 
		theta = twopi/n
		do i=1,n
			x1(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x2(i+nold) = G0/x1(i+nold)
			x3(i+nold) = czip
			x4(i+nold) = czip
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
		n = 5
		theta = twopi/n
		do i=1,n
			x1(i+nold) = czip
			x2(i+nold) = czip
			x3(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x4(i+nold) = -G0/x3(i+nold)
		enddo
		nold = nold + n
		!------------------------------------------------------------------------------!
			x1(11) = cone
			x2(11) = G0
			x3(11) = cone
			x4(11) = czip

			x1(12) = two
			x2(12) = G0/two
			x3(12) = cone
			x4(12) = czip

			x1(13) = cone
			x2(13) = G0
			x3(13) = czip
			x4(13) = cone

			x1(14) = two
			x2(14) = G0/two
			x3(14) = czip
			x4(14) = cone
	endif
		!------------------------------------------------------------------------------!
		G0mu = (msq(j1)+mu2g(1))/two/MP12(1)
		x1(15) = cone
		x2(15) = G0mu
		x3(15) = czip
		x4(15) = czip
		nold = nold + n
	do i=1,15
		q1(i,:)=-p0(:)+x1(i)*e1(:)+x2(i)*e2(:)+x3(i)*e3(:)+x4(i)*e4(:)
	enddo
! DEFAULT CASE=================================================================================!
	case default
		print *, 'error in mgetqs.f90, getq1s'
	end select
! TEST VARIABLE================================================================================!
	!Defining the test variable qt
	mu2t(1) =   half
	varx    = -(msq(j1)+mu2t(1))/two/MP12(1)
	x1t     =   3.2_ki
	x2t     =   1.2_ki
	x3t     =   4.2_ki
	x4t     =  (varx+x1t*x2t)/x3t
	qt(:)   =  -p0(:)+x1t*e1(:)+x2t*e2(:)+x3t*e3(:)+x4t*e4(:)
!==============================================================================================!
  end subroutine getq1_rm



























!CM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine getq5_cm(nleg,cut5,e1,e2,e3,e4,p0,Vi,msq,r1,r2,q5,mu2)
      implicit none
      integer, intent(in) :: nleg, cut5
      real(ki), dimension(4), intent(in) :: e1, e2, p0
      complex(ki), dimension(4), intent(in) :: e3, e4
      real(ki), intent(in) :: r1, r2
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      complex(ki), dimension(4), intent(out) :: q5
      complex(ki), intent(out) :: mu2

      integer :: j1,j2,j3,j4,j5
      complex(ki) :: x1, x2,var2, var3
      real(ki) :: MP12, MP1v2, MP1v3, MP2v2, MP2v3, MPv22, MPv33
      real(ki), dimension(4) :: v2, v3
      complex(ki) :: MP3v2, MP3v3, MP4v2, MP4v3, x3, x4, tmu2, den

      j5=cut5/10000
      j4=(cut5-j5*10000)/1000
      j3=(cut5-j5*10000-j4*1000)/100
      j2=(cut5-j5*10000-j4*1000-j3*100)/10
      j1= cut5-j5*10000-j4*1000-j3*100-j2*10

      v2(:)=Vi(j3,:)-p0(:)
      v3(:)=Vi(j4,:)-p0(:)

      MP12 =sdot(e1,e2)
      MP1v2=sdot(v2,e1)
      MP1v3=sdot(v3,e1)
      MP2v2=sdot(v2,e2)
      MP2v3=sdot(v3,e2)
      MP3v2=sdot(v2,e3)
      MP3v3=sdot(v3,e3)
      MP4v2=sdot(v2,e4)
      MP4v3=sdot(v3,e4)
      MPv22=sdot(v2,v2)
      MPv33=sdot(v3,v3)

      den = -(MP3v3*MP4v2*two) + MP3v2*MP4v3*two

      x1 = (0.5_ki*(  -(MP12*(1.0_ki+r1)*r2*two)-(1.0_ki+r2)*msq(j1) &
     & +r2*msq(j2) + msq(j5)))/(MP12*(-1.0_ki + r1*r2))

      x2 = -((0.5_ki*(-(MP12*r1*(1.0_ki+r2)*two)-(1.0_ki+r1)*msq(j1) &
     & +msq(j2) + r1*msq(j5)))/(MP12*(-1.0_ki + r1*r2)))

      var2 = MPv22 + MP1v2*two*x1 + MP2v2*two*x2+msq(j1)-msq(j3)
      var3 = MPv33 + MP1v3*two*x1 + MP2v3*two*x2+msq(j1)-msq(j4)

      x3 = (-MP4v3*var2+MP4v2*var3)/den

      x4 = ( MP3v3*var2-MP3v2*var3)/den

      tmu2=two*MP12*(x1*x2-x3*x4)-msq(j1)
      if (aimag(tmu2)/real(tmu2).lt.1e-10_ki) then
      mu2=real(tmu2,ki)*cone
      else
      mu2=tmu2
      endif
      q5(:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3*e3(:)+x4*e4(:)
   end subroutine getq5_cm

   subroutine getq4_cm(nleg,cut4,e1,e2,e3,e4,p0,k1,k2,k3,L3,r1,r2,q4,qt,msq)
      use mglobal, only: MP12, mu2g, mu2t, dx
      implicit none
      integer, intent(in) :: nleg, cut4
      real(ki), dimension(4) :: k1, k2, k3, L3, e1, e2, p0
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), intent(in) :: r1,r2
      complex(ki), dimension(6,4), intent(out) :: q4
      complex(ki), dimension(4), intent(out) :: qt

      real(ki) :: beta,mu2x
      complex(ki) :: x1,x2,A1,A2
      real(ki) :: ME31, ME32, MK11, MK22, MK33
      complex(ki) :: x31,x32,x33,x34,x35,x36,x3t,x41,x42,x43,x44,x45,x46,x4t
      complex(ki) :: ME33,ME34,ML33,ML34
      complex(ki) :: B0,B1,B2,rtdel,A3,A4
      integer :: j1,j2,j3,j4

      j4=cut4/1000
      j3=(cut4-j4*1000)/100
      j2=(cut4-j4*1000-j3*100)/10
      j1= cut4-j4*1000-j3*100-j2*10

      MP12(4)=sdot(e1,e2)
      ML33=sdot(L3,e3)
      ML34=sdot(L3,e4)
      ME31=sdot(k3,e1)
      ME32=sdot(k3,e2)
      ME33=sdot(k3,e3)
      ME34=sdot(k3,e4)
      MK11=sdot(k1,k1)
      MK22=sdot(k2,k2)
      MK33=sdot(k3,k3)
      beta=one/(one-r1*r2)
      A1=(msq(j2)-msq(j1)-MK11)/two/MP12(4)
      A2=(msq(j1)-msq(j4)+MK22)/two/MP12(4)
      x1=beta*(A2-r2*A1)
      x2=beta*(A1-r1*A2)
      A3=(msq(j3)-msq(j1)-MK33-2d0*x1*ME31-2d0*x2*ME32)/2d0/ME33
      A4=-ME34/ME33
      B1=-two*MP12(4)*A3
      B2=-two*MP12(4)*A4

      mu2g(4)=zip
      B0=(two*x1*x2*MP12(4)-msq(j1)-mu2g(4))*cone
      rtdel=sqrt(B1**2-four*B0*B2)
      x41=(-B1+rtdel)/two/B2
      x42=(-B1-rtdel)/two/B2
      x31=x41*A4+A3
      x32=x42*A4+A3
      q4(1,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x31*e3(:)+x41*e4(:)
      q4(2,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x32*e3(:)+x42*e4(:)

!--- scelta dinamica
      mu2g(4)=max(abs(msq(j1)),abs(msq(j2)),abs(msq(j3)),abs(msq(j4)),&
     & abs(MK11),abs(MK22),abs(MK33),&
     & abs(two*sdot(k1,k2)),abs(two*sdot(k2,k3)),abs(two*sdot(k1,k3)),&
     & abs(two*MP12(4)))

      if (abs(mu2g(4)).lt.1.0e-10_ki) mu2g(4)=one

!--- scelta statica
!      mu2g(4)=one

      B0=(two*x1*x2*MP12(4)-msq(j1)-mu2g(4))*cone
      rtdel=sqrt(B1**2-four*B0*B2)
      x43=(-B1+rtdel)/two/B2
      x44=(-B1-rtdel)/two/B2
      x33=x43*A4+A3
      x34=x44*A4+A3
      q4(3,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x33*e3(:)+x43*e4(:)
      q4(4,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x34*e3(:)+x44*e4(:)

      mu2x=-mu2g(4)
      B0=(two*x1*x2*MP12(4)-msq(j1)-mu2x)*cone
      rtdel=sqrt(B1**2-four*B0*B2)
      x45=(-B1+rtdel)/two/B2
!      x45=(-B1-rtdel)/two/B2
      x35=x45*A4+A3
!      x35=x4t*A4+A3
      q4(5,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x35*e3(:)+x45*e4(:)
!EDITED BY HvD
      x46=(-B1-rtdel)/two/B2
      x36=x46*A4+A3
      q4(6,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x36*e3(:)+x46*e4(:)

      mu2t(4)=half
      B0=(two*x1*x2*MP12(4)-msq(j1)-mu2t(4))*cone
      rtdel=sqrt(B1**2-four*B0*B2)
!      x4t=(-B1+rtdel)/two/B2
      x4t=(-B1-rtdel)/two/B2
!      x3t=x45*A4+A3
      x3t=x4t*A4+A3
      qt(:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3t*e3(:)+x4t*e4(:)

      dx(1)=ML33*x31-ML34*x41
      dx(2)=ML33*x32-ML34*x42
      dx(3)=ML33*x33-ML34*x43
      dx(4)=ML33*x34-ML34*x44
      dx(5)=ML33*x35-ML34*x45
   end subroutine getq4_cm

   subroutine getq3_cm(nleg,irank,cut3,e1,e2,e3,e4,p0,k1,k2,msq,r1,r2,q3,qt)
      use mglobal, only: C0c,C1c,mu2g,MP12,KK,Kmu,mu2t
      use options, only: C0_thrs
      implicit none
      integer, intent(in) :: nleg, cut3, irank
      real(ki), dimension(4), intent(in) :: p0, k1, k2, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(10,4), intent(out) :: q3
      complex(ki), dimension(4), intent(out) :: qt

      integer :: j1,j2,j3,j,ndiff
      real(ki) :: r1,r2,beta,teta
      complex(ki) :: x1,x2,A1,A2,C1t
      real(ki) :: MK11,MK12,MK22
      complex(ki) :: x3t, x4t
      complex(ki), dimension(10) :: x3, x4

      ndiff=nleg-irank
      KK(3)=one
      Kmu(3)=one

      j3=cut3/100
      j2=(cut3-j3*100)/10
      j1=cut3-j3*100-j2*10

      MP12(3)=sdot(e1,e2)
      MK11=sdot(k1,k1)
      MK12=sdot(k1,k2)
      MK22=sdot(k2,k2)
      beta=1d0/(1d0-r1*r2)
      A1=(msq(j2)-msq(j1)-MK11)/two/MP12(3)
      A2=(msq(j1)-msq(j3)+MK22)/two/MP12(3)
      x1=beta*(A2-r2*A1)
      x2=beta*(A1-r1*A2)

      C0c=x1*x2-msq(j1)/two/MP12(3)

      !---#[ scelta dinamica:
      mu2g(3)=max(abs(msq(j1)),abs(msq(j2)),abs(msq(j3)),abs(MK11), &
                & abs(MK22),abs(two*MK12))
      if (abs(mu2g(3)).lt.1d-10) mu2g(3)=one
      !---#] scelta dinamica:
      !---#[ scelta statica:
      ! mu2g(3)=one
      !---#] scelta statica:

      C1c=x1*x2-(msq(j1)+mu2g(3))/two/MP12(3)

      if (abs(C0c-1.0_ki) .lt. C0_thrs) then
         !---#[ New Sampling:
         ! The new sampling is the one that is safe around C0=1
         ! but not around C0=0

         ! complete c-system: 10 coefficients

         teta=twopi/seven
         do j=1,7
            x3(j)=cos(teta*real(j-1,ki))-im*sin(teta*real(j-1,ki))
            x4(j)=C0c/x3(j)
         enddo
         
         teta=twopi/three
         do j=8,10
            x4(j)=cos(teta*real(j-8,ki))+im*sin(teta*real(j-8,ki))
            x3(j)=C1c/x4(j)
         enddo
   
         mu2t(3)=half
         C1t=x1*x2-(msq(j1)+mu2t(3))/two/MP12(3)
         x3t=32.4_ki
         x4t=C1t*cone/x3t

         q3(1,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(1)*e3(:)+x4(1)*e4(:)
         q3(2,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(2)*e3(:)+x4(2)*e4(:)
         q3(3,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(3)*e3(:)+x4(3)*e4(:)
         q3(4,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(4)*e3(:)+x4(4)*e4(:)
         q3(5,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(5)*e3(:)+x4(5)*e4(:)
         q3(6,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(6)*e3(:)+x4(6)*e4(:)
         q3(7,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(7)*e3(:)+x4(7)*e4(:)
         q3(8,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(8)*e3(:)+x4(8)*e4(:)
         q3(9,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(9)*e3(:)+x4(9)*e4(:)
         q3(10,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(10)*e3(:)+x4(10)*e4(:)
         !---#] New Sampling:
      else
         !---#[ Old Sampling:
         ! The old sampling is the one that is safe around C0=0
         ! but not around C0=1

         !!!! HERE WE BRANCH ACCORDING TO THE RANK
         select case(ndiff)
         case(1)
            ! rank2 c-system: 6  coefficients
            teta=twopi/three
            do j=1,3
               x3(j)=cos(teta*real(j-1,ki))-im*sin(teta*real(j-1,ki))
               x4(j)=C0c/x3(j)
            enddo
         
            teta=twopi/two
            do j=4,5
               x3(j)=cos(teta*real(j-4,ki))+im*sin(teta*real(j-4,ki))
               x4(j)=cone/x3(j)
            enddo
         
            teta=twopi
            do j=6,6
               x3(j)=cos(teta*real(j-6,ki))+im*sin(teta*real(j-6,ki))
               x4(j)=cone/x3(j)
            enddo

            x3(7:10)=czip
            x4(7:10)=czip
   
            ! those numbers are for the test!!
            mu2t(3)=half
            C1t=x1*x2-(msq(j1)+mu2t(3))/two/MP12(3)
            x3t=32.4_ki
            x4t=C1t*cone/x3t

            ! here are the qs for system 3.2
            q3(1,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(1)*e3(:)+x4(1)*e4(:)
            q3(2,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(2)*e3(:)+x4(2)*e4(:)
            q3(3,:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3(3)*e3(:)+x4(3)*e4(:)
            q3(4,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C0c*x3(4)*e3(:)+x4(4)*e4(:)
            q3(5,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C0c*x3(5)*e3(:)+x4(5)*e4(:)
            q3(6,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C1c*x3(6)*e3(:)+x4(6)*e4(:)
         case(2)
            ! rank1 c-system: 3 coefficients
       
            x3(1)=cone
            x4(1)=C0c
        
            x3(2)=-cone
            x4(2)=-C0c
         
            x3(3)=C0c
            x4(3)=cone

            x3(4:10)=czip
            x4(4:10)=czip
           
            mu2t(3)=half
            C1t=x1*x2-(msq(j1)+mu2t(3))/two/MP12(3)
            x3t=32.4_ki
            x4t=C1t*cone/x3t

            q3(1,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(1)*e3(:)+x4(1)*e4(:)/KK(3)
            q3(2,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(2)*e3(:)+x4(2)*e4(:)/KK(3)
            q3(3,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(3)*e3(:)+x4(3)*e4(:)/KK(3)
         case default
            ! complete c-system: 10 coefficients

            teta=twopi/four
            do j=1,4
               x3(j)=cos(teta*real(j-1,ki))-im*sin(teta*real(j-1,ki))
               x4(j)=C0c/x3(j)
            enddo

            teta=twopi/three
            do j=5,7
               x4(j)=cos(teta*real(j-5,ki))+im*sin(teta*real(j-5,ki))
               x3(j)=C0c/x4(j)
            enddo
         
            teta=twopi/three
            do j=8,10
               x4(j)=cos(teta*real(j-8,ki))+im*sin(teta*real(j-8,ki))
               x3(j)=C1c/x4(j)
            enddo
   
            mu2t(3)=half
            C1t=x1*x2-(msq(j1)+mu2t(3))/two/MP12(3)
            x3t=32.4_ki
            x4t=C1t*cone/x3t

            q3(1,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(1)*e3(:)+x4(1)*e4(:)/KK(3)
            q3(2,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(2)*e3(:)+x4(2)*e4(:)/KK(3)
            q3(3,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(3)*e3(:)+x4(3)*e4(:)/KK(3)
            q3(4,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(4)*e3(:)+x4(4)*e4(:)/KK(3)
            q3(5,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(5)*e3(:)+x4(5)*e4(:)/KK(3)
            q3(6,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(6)*e3(:)+x4(6)*e4(:)/KK(3)
            q3(7,:)=-p0(:)+x1*e1(:)+x2*e2(:)+KK(3)*x3(7)*e3(:)+x4(7)*e4(:)/KK(3)
            q3(8,:)=-p0(:)+x1*e1(:)+x2*e2(:)+Kmu(3)*x3(8)*e3(:)+x4(8)*e4(:)/Kmu(3)
            q3(9,:)=-p0(:)+x1*e1(:)+x2*e2(:)+Kmu(3)*x3(9)*e3(:)+x4(9)*e4(:)/Kmu(3)
            q3(10,:)=-p0(:)+x1*e1(:)+x2*e2(:)+Kmu(3)*x3(10)*e3(:)+x4(10)*e4(:)/Kmu(3)
         end select
         !---#] Old Sampling:
      end if
      qt(:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3t*e3(:)+x4t*e4(:)
   end subroutine getq3_cm


   subroutine getq2_cm(nleg,irank,cut2,e1,e2,e3,e4,p0,k1,msq,q2,qt)
      use mglobal, only: Fpc,Fzc,Fmc,F1zc,KB,mu2g,KK,Kmu,MP12,mu2t,Fp2c,F11c
      implicit none
      integer, intent(in) :: nleg, cut2, irank
      real(ki), dimension(4), intent(in) :: p0, k1, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(20,4), intent(out) :: q2
      complex(ki), dimension(4), intent(out) :: qt

      complex(ki), dimension(20) :: x1, x2
      real(ki) :: teta,K1SQ
      complex(ki) :: A1,x1t,x2t,Ft
      complex(ki) :: X1p,X2p,X1z,X2z,X1m,X2m, X1p2,X2p2
      complex(ki) :: x3t, x4t
      complex(ki), dimension(20) :: x3, x4
      integer :: j1,j2,j,ndiff
      
      ndiff=nleg-irank
      KK(2)=one
      KB=one
      Kmu(2)=one

      j2=cut2/10
      j1=cut2-j2*10

      K1SQ=sdot(k1,k1)
      MP12(2)=sdot(e1,e2)
      A1=(msq(j2)-msq(j1)-K1SQ)/two/MP12(2)

      X1z=zip
      X2z=A1
      Fzc=-msq(j1)/two/MP12(2)

      X1p=KB
      X2p=A1-KB*K1SQ/two/MP12(2)
      Fpc=KB*X2p-msq(j1)/two/MP12(2)

      X1m=-KB
      X2m=A1+KB*K1SQ/two/MP12(2)
      Fmc=-KB*X2m-msq(j1)/two/MP12(2)

!Added by HvD
      X1p2=two
      X2p2=A1-two*KB*K1SQ/two/MP12(2)
      Fp2c=two*KB*X2p-msq(j1)/two/MP12(2)
      

!--- scelta dinamica

      mu2g(2)=max(abs(msq(j1)),abs(msq(j2)),K1SQ)
      if (abs(mu2g(2)).lt.1.0e-10_ki) mu2g(2)=one
 
!--- scelta statica
!      mu2g(2)=one

      F1zc=-(msq(j1)+mu2g(2))/two/MP12(2)
!Added by HvD
      F11c=KB*X2p-(msq(j1)+mu2g(2))/two/MP12(2)

!!!! HERE WE BRANCH ACCORDING TO THE RANK

      if (ndiff.eq.1) then
! rank 1
          
         teta=twopi/two

         do j=1,2
            x1(j)=X1z
            x2(j)=X2z
            x4(j)=cos(teta*real(j-1,ki))+im*sin(teta*real(j-1,ki))
            x3(j)=Fzc/x4(j)
         enddo

         do j=3,4
            x1(j)=X1p
            x2(j)=X2p
            x3(j)=cos(teta*real(j-3,ki))-im*sin(teta*real(j-3,ki))
            x4(j)=Fpc/x3(j)
         enddo
         
       
      mu2t(2)=half
      x1t=2.3_ki
      x2t=A1-x1t*K1SQ/two/MP12(2)
      Ft=x1t*x2t-(msq(j1)+mu2t(2))/two/MP12(2)
      x3t=(2.5_ki,0.3_ki)
      x4t=Ft/x3t

      do j=5,10
         x3(j)=czip
         x4(j)=czip
      enddo
         
      q2(1,:)=-p0(:)+x1(1)*e1(:)+x2(1)*e2(:)+KK(2)*x3(1)*e3(:)+x4(1)*e4(:)/KK(2)
      q2(2,:)=-p0(:)+x1(2)*e1(:)+x2(2)*e2(:)+KK(2)*x3(2)*e3(:)+x4(2)*e4(:)/KK(2)
      q2(3,:)=-p0(:)+x1(3)*e1(:)+x2(3)*e2(:)+KK(2)*x3(3)*e3(:)+x4(3)*e4(:)/KK(2)
      q2(4,:)=-p0(:)+x1(4)*e1(:)+x2(4)*e2(:)+KK(2)*x3(4)*e3(:)+x4(4)*e4(:)/KK(2)
 !     q2(5,:)=-p0(:)+x1(5)*e1(:)+x2(5)*e2(:)+KK(2)*x3(5)*e3(:)+x4(5)*e4(:)/KK(2)
 !     q2(6,:)=-p0(:)+x1(6)*e1(:)+x2(6)*e2(:)+KK(2)*x3(6)*e3(:)+x4(6)*e4(:)/KK(2)
!      q2(7,:)=-p0(:)+x1(7)*e1(:)+x2(7)*e2(:)+KK(2)*x3(7)*e3(:)+x4(7)*e4(:)/KK(2)
!      q2(8,:)=-p0(:)+x1(8)*e1(:)+x2(8)*e2(:)+KK(2)*x3(8)*e3(:)+x4(8)*e4(:)/KK(2)
!      q2(9,:)=-p0(:)+x1(9)*e1(:)+x2(9)*e2(:)+KK(2)*x3(9)*e3(:)+x4(9)*e4(:)/KK(2)
!      q2(10,:)=-p0(:)+x1(10)*e1(:)+x2(10)*e2(:)+KK(2)*x3(10)*e3(:) &
!     & +x4(10)*e4(:)/KK(2)
      qt(:)=-p0(:)+x1t*e1(:)+x2t*e2(:)+x3t*e3(:)+x4t*e4(:)
   

      else
! standard (rank 2)

      teta=twopi/four
      do j=1,4
      x1(j)=X1z
      x2(j)=X2z
      x4(j)=cos(teta*real(j-1,ki))+im*sin(teta*real(j-1,ki))
      x3(j)=Fzc/x4(j)
      enddo
      teta=twopi/three
      do j=5,7
      x1(j)=X1z
      x2(j)=X2z
      x3(j)=cos(teta*real(j-5,ki))-im*sin(teta*real(j-5,ki))
      x4(j)=Fzc/x3(j)
      enddo

      teta=twopi/five
      do j=8,12
      x1(j)=X1p
      x2(j)=X2p
      x4(j)=cos(teta*real(j-8,ki))+im*sin(teta*real(j-8,ki))
      x3(j)=Fpc/x4(j)
      enddo

      teta=twopi/three
      do j=13,15
      x1(j)=X1m
      x2(j)=X2m
      x4(j)=cos(teta*real(j-13,ki))+im*sin(teta*real(j-13,ki))
      x3(j)=Fmc/x4(j)
      enddo
      
      x1(16)=X1p2
      x2(16)=X2p2
      x3(16)=Fp2c
      x4(16)=cone

      teta=twopi/three
      do j=17,19
      x1(j)=X1z
      x2(j)=X2z
      x4(j)=cos(teta*real(j-17,ki))+im*sin(teta*real(j-17,ki))
      x3(j)=F1zc/x4(j)
      enddo


      x1(20)=X1p
      x2(20)=X2p
      x4(20)=cone
      x3(20)=F11c/x4(20)






      mu2t(2)=half
      x1t=2.3_ki
      x2t=A1-x1t*K1SQ/two/MP12(2)
      Ft=x1t*x2t-(msq(j1)+mu2t(2))/two/MP12(2)
      x3t=(2.5_ki,0.3_ki)
      x4t=Ft/x3t

      q2(1,:)=-p0(:)+x1(1)*e1(:)+x2(1)*e2(:)+KK(2)*x3(1)*e3(:)+x4(1)*e4(:)/KK(2)
      q2(2,:)=-p0(:)+x1(2)*e1(:)+x2(2)*e2(:)+KK(2)*x3(2)*e3(:)+x4(2)*e4(:)/KK(2)
      q2(3,:)=-p0(:)+x1(3)*e1(:)+x2(3)*e2(:)+KK(2)*x3(3)*e3(:)+x4(3)*e4(:)/KK(2)
      q2(4,:)=-p0(:)+x1(4)*e1(:)+x2(4)*e2(:)+KK(2)*x3(4)*e3(:)+x4(4)*e4(:)/KK(2)
      q2(5,:)=-p0(:)+x1(5)*e1(:)+x2(5)*e2(:)+KK(2)*x3(5)*e3(:)+x4(5)*e4(:)/KK(2)
      q2(6,:)=-p0(:)+x1(6)*e1(:)+x2(6)*e2(:)+KK(2)*x3(6)*e3(:)+x4(6)*e4(:)/KK(2)
      q2(7,:)=-p0(:)+x1(7)*e1(:)+x2(7)*e2(:)+KK(2)*x3(7)*e3(:)+x4(7)*e4(:)/KK(2)
      q2(8,:)=-p0(:)+x1(8)*e1(:)+x2(8)*e2(:)+KK(2)*x3(8)*e3(:)+x4(8)*e4(:)/KK(2)
      q2(9,:)=-p0(:)+x1(9)*e1(:)+x2(9)*e2(:)+KK(2)*x3(9)*e3(:)+x4(9)*e4(:)/KK(2)
      q2(10,:)=-p0(:)+x1(10)*e1(:)+x2(10)*e2(:)+KK(2)*x3(10)*e3(:)+x4(10)*e4(:)/KK(2)
      q2(11,:)=-p0(:)+x1(11)*e1(:)+x2(11)*e2(:)+KK(2)*x3(11)*e3(:)+x4(11)*e4(:)/KK(2)
      q2(12,:)=-p0(:)+x1(12)*e1(:)+x2(12)*e2(:)+KK(2)*x3(12)*e3(:)+x4(12)*e4(:)/KK(2)
      q2(13,:)=-p0(:)+x1(13)*e1(:)+x2(13)*e2(:)+KK(2)*x3(13)*e3(:)+x4(13)*e4(:)/KK(2)
      q2(14,:)=-p0(:)+x1(14)*e1(:)+x2(14)*e2(:)+KK(2)*x3(14)*e3(:)+x4(14)*e4(:)/KK(2)
      q2(15,:)=-p0(:)+x1(15)*e1(:)+x2(15)*e2(:)+KK(2)*x3(15)*e3(:)+x4(15)*e4(:)/KK(2)
      q2(16,:)=-p0(:)+x1(16)*e1(:)+x2(16)*e2(:)+KK(2)*x3(16)*e3(:)+x4(16)*e4(:)/KK(2)
      q2(17,:)=-p0(:)+x1(17)*e1(:)+x2(17)*e2(:)+KK(2)*x3(17)*e3(:)+x4(17)*e4(:)/KK(2)
      q2(18,:)=-p0(:)+x1(18)*e1(:)+x2(18)*e2(:)+KK(2)*x3(18)*e3(:)+x4(18)*e4(:)/KK(2)
      q2(19,:)=-p0(:)+x1(19)*e1(:)+x2(19)*e2(:)+KK(2)*x3(19)*e3(:)+x4(19)*e4(:)/KK(2)
      q2(20,:)=-p0(:)+x1(20)*e1(:)+x2(20)*e2(:)+KK(2)*x3(20)*e3(:)+x4(20)*e4(:)/KK(2)
      qt(:)=-p0(:)+x1t*e1(:)+x2t*e2(:)+x3t*e3(:)+x4t*e4(:)
      
   endif

   end subroutine getq2_cm

   subroutine getq1_cm(nleg,rank,cut1,e1,e2,e3,e4,p0,msq,q1,qt)
     use mglobal, only: G0c,KK,mu2g,MP12,mu2t
      implicit none
      integer, intent(in) :: nleg, rank, cut1
      real(ki), dimension(4), intent(in) :: p0, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(5,4), intent(out) :: q1
      complex(ki), dimension(4), intent(out) :: qt
      complex(ki), dimension(0:nleg-1), intent(in) :: msq

      complex(ki), dimension(5) :: x1, x2
      complex(ki), dimension(5) :: x3, x4
      complex(ki) :: varx,x1t,x2t,x3t,x4t
      integer :: j1

      KK(1)=one
      j1=cut1

      mu2g(1)=zip
      if (abs(msq(j1)).eq.zip) then
      mu2g(1)=one
      endif

      MP12(1)=sdot(e1,e2)
      G0c=(msq(j1)+mu2g(1))/two/MP12(1)

      x1(1)=G0c
      x2(1)=cone
      x3(1)=czip
      x4(1)=czip

      x1(2)=-G0c
      x2(2)=-cone
      x3(2)=czip
      x4(2)=czip

      x1(3)=G0c
      x2(3)=cone
      x3(3)=cone
      x4(3)=czip

      x1(4)=czip
      x2(4)=cone
      x3(4)=-cone
      x4(4)=G0c

      x1(5)=czip
      x2(5)=cone
      x3(5)=cone
      x4(5)=-G0c


      mu2t(1)=half
      varx=-(msq(j1)+mu2t(1))/two/MP12(1)
      x1t=3.2_ki
      x2t=1.2_ki
      x3t=4.2_ki
      x4t=(varx+x1t*x2t)/x3t

      q1(1,:)=-p0(:)+x1(1)*e1(:)+x2(1)*e2(:)+x3(1)*e3(:)+x4(1)*e4(:)
      q1(2,:)=-p0(:)+x1(2)*e1(:)+x2(2)*e2(:)+x3(2)*e3(:)+x4(2)*e4(:)
      q1(3,:)=-p0(:)+x1(3)*e1(:)+x2(3)*e2(:)+x3(3)*e3(:)+x4(3)*e4(:)
      q1(4,:)=-p0(:)+x1(4)*e1(:)+x2(4)*e2(:)+x3(4)*e3(:)+x4(4)*e4(:)
      q1(5,:)=-p0(:)+x1(5)*e1(:)+x2(5)*e2(:)+x3(5)*e3(:)+x4(5)*e4(:)
      qt(:)=-p0(:)+x1t  *e1(:)+x2t  *e2(:)+x3t  *e3(:)+x4t  *e4(:)
  end subroutine getq1_cm



































!Added by HvD
subroutine SampleL2(mu2nonzero,x1samplevalue,n,samplex3,nold,x1,x2,x3,x4,Fstep,A2,B2)
!---------------------------------------------------------------!
!Author:	Hans van Deurzen				!
!Date:		25-05-2012					!
!Description:	Sampling for level 2				!
!---------------------------------------------------------------!
	use mglobal, only: FL2,mu2g,mu2vec2
	implicit none
! PARAMETERS=====================================================================================!
	!external parameters
		logical,     intent(in   )   :: mu2nonzero    !mu2=0 or mu2!=0
		logical,     intent(in   )   :: samplex3      !sample x3 or x4
		real(ki),    intent(in   )   :: x1samplevalue !
		integer,     intent(in   )   :: n	      !number of samplings
		integer,     intent(in   )   :: Fstep	      !
		real(ki),    intent(in   )   :: A2,B2	      !
		integer,     intent(inout)   :: nold	      !number of samplings before these
		real(ki),    intent(inout)   :: x1(:), x2(:)  !
		complex(ki), intent(inout)   :: x3(:), x4(:)  !
!		real(ki),    dimension(20), intent(inout) :: x1,x2
!		complex(ki), dimension(20), intent(inout) :: x3,x4
	!internal parameters
		real(ki)	:: Finternal, x1internal, x2internal, theta
		complex(ki) 	:: mu2internal
		integer 	:: i
! INITIALIZATION=================================================================================!
	x1internal = x1samplevalue
	Finternal  = FL2(Fstep)
	!Determine mu2
	if(mu2nonzero) then
		mu2internal=mu2g(2)
	else !mu2==0
		mu2internal=czip		
	endif
! SAMPLE X3 OR X4 AND ASSIGN ALL MU2, X1 AND X2 THE VALUES DETERMINED ABOVE======================!
	theta=twopi/n
	do i=1,n
		x1(nold+i)=x1internal
		x2(nold+i)=A2*x1internal+B2
		if   (samplex3) then
			x3(i+nold) = cos(theta*real(i-1,ki))-im*sin(theta*real(i-1,ki))
			x4(i+nold) = Finternal/x3(i+nold)
		else !sample x4
			x4(i+nold) = cos(theta*real(i-1,ki))+im*sin(theta*real(i-1,ki))
			x3(i+nold) = Finternal/x4(i+nold)
			mu2vec2(i+nold) = mu2internal
		endif
	enddo
	nold = nold + n
end subroutine SampleL2

subroutine FillFL2(x1,A,B,C,Fstep)
! Product x3*x4 F(mu2,x1) in level2
	use mglobal, only: FL2
	implicit none
	!external parameters
	real(ki), intent(in)    :: x1
	real(ki), intent(in)    :: A,B
	complex(ki), intent(in) :: C
	integer,  intent(in)    :: Fstep
	FL2(Fstep) = A*x1**2 + B*x1 + C
end subroutine FillFL2

end module mgetqs

