module mgetc3
   use precision, only: ki
   use constants
   use options
   use mfunctions
   use mrestore
   implicit none

   public :: getc3

contains

subroutine getc3(numeval,nleg,rank,c3,cut3,q3,qt,Vi,msq)
	use mglobal, only: C0,C1,Cm1,mu2g,MP12,KK,mu2t,resit,denst,mu2test
	use options, only: C0_thrs
	implicit none
! PARAMETERS=======================================================================!
	! external
		integer, 			    intent(in ) :: nleg, rank, cut3
		complex(ki), dimension(0:14), 	    intent(out) :: c3
		complex(ki), dimension(15,4), 	    intent(in ) :: q3
		complex(ki), dimension(4), 	    intent(in ) :: qt
		real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
		complex(ki),    dimension(0:nleg-1),   intent(in ) :: msq
		complex(ki), dimension(15) 		        :: mu2vec
	! internal
		integer :: i,m,n,i1,i2,i3,i4,i5,dicut5,dicut4,mx,j1,j2,j3,diff
		integer :: nsol, acc, nloops
		complex(ki), dimension(15)   :: dens3,dens4,dens5,xneval
		complex(ki), dimension(15)   :: resi5,resi4,known
		complex(ki), dimension(0:14) :: f3
		complex(ki) 		     :: dens4t,dens5t
		logical 		     :: evalres
		integer,     dimension(5)    :: iarr
	! haggies !!! TR: I have used haggies to rewrite the systems.
		complex(ki) :: t1
		complex(ki) :: t2
		complex(ki) :: t3
		complex(ki) :: t4
		complex(ki) :: t5
		complex(ki) :: t6
		complex(ki) :: t7
		complex(ki) :: t8
		complex(ki) :: t9
		complex(ki) :: t10
		complex(ki) :: t11
		complex(ki) :: t12
		complex(ki) :: t13
		complex(ki) :: t14
		complex(ki) :: t15
		complex(ki) :: t16
		complex(ki) :: t17
		complex(ki) :: t18
		complex(ki) :: t19
		complex(ki) :: t20
		complex(ki) :: t21
		complex(ki) :: t22
		complex(ki) :: t23
		complex(ki) :: t24
		complex(ki) :: t25
		complex(ki) :: t26
		complex(ki) :: t27
! INTERFACE========================================================================!
	interface
		function     numeval(ncut, Q, mu2)
			use precision
			implicit none
			integer, 		   intent(in) :: ncut
			complex(ki), dimension(4), intent(in) :: Q
			complex(ki), 		   intent(in) :: mu2
			complex(ki) 			      :: numeval
		end function numeval
	end interface
! INITIALIZATION===================================================================!
	! mu2
		mu2test(3) = mu2t(3)
	! cuts
		j3  =  cut3/100
		acc =  j3*100
		j2  = (cut3-acc)/10
		j1  =  cut3-acc-j2*10
	! ---
		resi5(:)  = czip
		resi4(:)  = czip
		known(:)  = czip
		xneval(:) = czip
		dens3(:)  = cone
	! for lnntest
		resit(3)  = czip
		denst(3)  = cone
	! for simplified sampling
		diff      = nleg-rank
! CALCULATING THE MU2 VECTOR=======================================================!
	do n=1,15
		mu2vec(n)=czip
	enddo
	if 	(diff .ge. 3) then ! rank0:  1 coefficients               
		nsol=1
	elseif  (diff .eq. 2) then ! rank1:  3 coefficients               
		nsol=3
	elseif  (diff .eq. 1) then ! rank2:  6 coefficients
		nsol=6
		mu2vec(6)  =  mu2g(3)
	elseif  (diff .eq. 0) then ! rank3: 10 coefficients
		nsol=10
		mu2vec(8)  =  mu2g(3)
		mu2vec(9)  =  mu2g(3)
		mu2vec(10) =  mu2g(3)
	elseif  (diff .eq.-1) then ! rank4: 15 coefficients	
		nsol=15
		do n=10,14
		mu2vec(n)  =  mu2g(3)
		enddo
		mu2vec(15) = -mu2g(3)
	else
		print *, 'error in mgetc3.f90'
	endif
! CALCULATING THE RESIDUES=========================================================!
! DIFF >= 3: only c3(0)
! DIFF < 3=========================================================================!
	if(nleg.eq.5) then
		do n=1,nsol
			resi5(n)=Res5(1,mu2vec(n))
		enddo
		resit(3)=Res5(1,mu2t(3))         
	elseif(nleg.eq.4) then
		do n=1,nsol
			resi4(n)=Res4(diff,1,q3(n,:),mu2vec(n))
		enddo
		resit(3)=Res4(diff,1,qt,mu2t(3))
        elseif(nleg.eq.3) then
		continue
        else
		!-------------------
		nloops = 5
		call HResL3(diff,resi5,	nsol,nloops,nleg,cut3,q3,qt,Vi,msq,mu2vec)
		!------------------
        endif

	if((nleg .ne. 4) .and. (nleg .ne. 3)) then
		!------------------
		nloops=4
		call HResL3(diff,resi4,	nsol,nloops,nleg,cut3,q3,qt,Vi,msq,mu2vec)
		!------------------
	endif
	if(nleg .ne. 3) then
		do i=0,nleg-1
			if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
				do n=1,nsol
				dens3(n)=dens3(n)*denevalmu2(nleg,i,q3(n,:),Vi,msq,mu2vec(n))
				enddo
				denst(3)=denst(3)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
			endif
		enddo
	endif
! CALCULATING THE POLYNOMIALS======================================================!
	do n=1,nsol
		xneval(n)=numeval(cut3,q3(n,:),mu2vec(n))
	enddo
	if     (imeth.eq.'diag') then
		do n=1,nsol
			known(n)=(xneval(n)-resi5(n)-resi4(n))/dens3(n)
		enddo
	elseif (imeth.eq.'tree') then
		do n=1,nsol
			known(n)=xneval(n)-(resi5(n)+resi4(n))/dens3(n)
		enddo
	endif
! C0=1=============================================================================!
! This sampling is the one that is safe around C0=1 , but not around C0=0
	if (abs(C0-1.0_ki) .lt. C0_thrs) then
	select case(diff)
	case(-1)
		do m=0,8
			f3(m)=effe(known,1,9,m)
		enddo
		do m=9,13
			mx=m-9
			f3(m)=effe(known,10,5,mx)
		enddo
		f3(14)=effe(known,15,1,0)
! CACULATING C3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      c3(0) =  &
     &f3(0) 
      c3(1) =  &
     &-(f3(8)/(C0*MP12(3))) 
      c3(2) =  &
     &f3(7)/(C0**2*MP12(3)**2) 
      c3(3) =  &
     &-(f3(6)/(C0**3*MP12(3)**3)) 
      c3(4) =  &
     &-(f3(1)/MP12(3)) 
      c3(5) =  &
     &f3(2)/MP12(3)**2 
      c3(6) =  &
     &-(f3(3)/MP12(3)**3) 
      c3(12) =  &
     &f3(5)/(C0**4*MP12(3)**4) 
      c3(13) =  &
     &f3(4)/MP12(3)**4 
      c3(7) =  &
     &        (Cm1*f3(8))/(C0*mu2g(3)) +  &
     &  (half*(C0**4*C1**2*two*f3(1)*mu2g(3)**4 +  &
     &       C0**4*C1**2*two*f3(2)*mu2g(3)**4 +  &
     &       C0**4*C1**2*f3(3)*mu2g(3)**4 +  &
     &       C0**4*C1**5*Cm1**2*f3(3)*mu2g(3)**4 +  &
     &       C0**4*C1**2*f3(4)*mu2g(3)**4 +  &
     &       C0**4*C1**6*Cm1*f3(4)*mu2g(3)**4 +  &
     &       C1*f3(5)*mu2g(3)**4 +  &
     &       C1**2*Cm1**4*f3(5)*mu2g(3)**4 +  &
     &       C0*f3(6)*mu2g(3)**4 +  &
     &       C0*C1**2*Cm1**3*f3(6)*mu2g(3)**4 +  &
     &       C0**2*C1**2*Cm1**2*two*f3(7)*mu2g(3)**4 +  &
     &       C0**4*C1**2*f3(9)*mu2g(3)**4 -  &
     &       C0**4*C1*f3(10)*mu2g(3)**4 -  &
     &       C0**4*f3(11)*mu2g(3)**4 -  &
     &       C0**4*C1**2*Cm1**2*f3(12)*mu2g(3)**4 -  &
     &       C0**4*C1**2*Cm1*f3(13)*mu2g(3)**4 -  &
     &       C0**4*C1**2*f3(14)*mu2g(3)**4))/ &
     &   (C0**4*C1**2*mu2g(3)**5) 
      c3(8) =  &
     &-((-(C0*C1**4*f3(4)) - f3(8) + C0*f3(13))/(C0*MP12(3)*mu2g(3))) 
      c3(9) =  &
     &        -((-(C0**4*C1*f3(1)) - f3(5) + C0**4*f3(10))/ &
     &    (C0**4*C1*MP12(3)*mu2g(3))) 
      c3(10) =  &
     &        -((C0**2*C1**3*f3(3) + f3(7) - C0**2*f3(12))/ &
     &    (C0**2*MP12(3)**2*mu2g(3))) 
      c3(11) =  &
     &        -((C0**3*C1**2*f3(2) + f3(6) - C0**3*f3(11))/ &
     &    (C0**3*C1**2*MP12(3)**2*mu2g(3))) 
      c3(14) =  &
     &        -(f3(0)/mu2g(3)**2) - (Cm1*f3(8))/(C0*mu2g(3)**2) +  &
     &  (half*(-(C0**4*C1**2*two*f3(1)) - C0**4*C1**2*two*f3(2) -  &
     &       C0**4*C1**2*f3(3) - C0**4*C1**5*Cm1**2*f3(3) -  &
     &       C0**4*C1**2*f3(4) - C0**4*C1**6*Cm1*f3(4) -  &
     &       C1*f3(5) - C1**2*Cm1**4*f3(5) - C0*f3(6) -  &
     &       C0*C1**2*Cm1**3*f3(6) -  &
     &       C0**2*C1**2*Cm1**2*two*f3(7) + C0**4*C1**2*f3(9) +  &
     &       C0**4*C1*f3(10) + C0**4*f3(11) +  &
     &       C0**4*C1**2*Cm1**2*f3(12) + C0**4*C1**2*Cm1*f3(13) +  &
     &       C0**4*C1**2*f3(14)))/(C0**4*C1**2*mu2g(3)**2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! 
	case default
            do m=0,6
               f3(m)=effe(known,1,7,m)
            enddo

            do m=7,9
               mx=m-7
               f3(m)=effe(known,8,3,mx)
            enddo
            c3(0) =   f3(0)
            c3(1) = -(f3(6)/(C0*MP12(3)))
            c3(2) =   f3(5)/(C0**2*MP12(3)**2)
            c3(3) = -(f3(4)/(C0**3*MP12(3)**3))
            c3(4) = -(f3(1)/MP12(3))
            c3(5) =   f3(2)/MP12(3)**2
            c3(6) = -(f3(3)/MP12(3)**3)
            c3(7) = (-f3(0) + MP12(3)**3*(-((C1**3*f3(3))/MP12(3)**3) -&
            &f3(4)/(C0**3*MP12(3)**3)) + f3(7))/mu2g(3)
            c3(8) = (C1**2*f3(2) + f3(6)/C0 - f3(9))/(MP12(3)*mu2g(3))
            c3(9) = (C1*f3(1) + f3(5)/C0**2 - f3(8))/(C1*MP12(3)*mu2g(3))
	end select
! C0 IS NOT 1======================================================================!
! The old sampling is the one that is safe around C0=0
! but not around C0=1
	else
	select case(diff)
! CALCULATING THE POLYNOMIALS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	case(3)
		c3(0)=known(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	case(2)
               f3(0:2)=known(1:3)
               t1 = f3(0)
               t2 = f3(1)
               c3(0) = ((t1+t2)*0.5_ki)
               t3 = KK(3)
               t4 = MP12(3)
               c3(1) = (-1.0_ki)*(t2+((1.0_ki)+C0)*t1-(t2*C0+(2.0_ki)*f3(2)))/((C&
               &0*C0-(1.0_ki))*t4)*t3*0.5_ki
               c3(4) = (-1.0_ki)*(t3*c3(0)+t4*C0*c3(1)-t2*t3)/(t3*t3*t4)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	case(1)
               do m=0,2
                  f3(m)=effe(known,1,3,m)
               enddo
               
               do m=3,4
                  mx=m-3
                  f3(m)=effe(known,4,2,mx)
               enddo
               
               do m=5,5
                  mx=m-5
                  f3(m)=effe(known,6,1,mx)
               enddo
               t1 = f3(0)
               c3(0) = t1
               t2 = C0*C0
               t3 = t2*C0
   
               t4 = f3(2)
               t5 = f3(1)
               t6 = f3(4)
               t7 = MP12(3)
               t8 = t7*t7
               c3(2) = ((t1+(t4+t2*t5-t6*C0)*t2-f3(3))/((t3*t3-one)*t8))
               c3(4) = (t2*t7*c3(2)-t5/t7)
               t1 = t7*C0
               c3(1) = (-1.0_ki)*(t6+t1*c3(4))/t7
               c3(5) = ((t4+t1*c3(1))/t8)
               t1 = t7*C1
               c3(7) = ((f3(5)+t1*c3(4)+t7*c3(1)-(t8*c3(2)+t1*t1*c3(5)+c3(0)))/(m&
               &u2g(3)))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	case(0)
               do m=0,3
                  f3(m)=effe(known,1,4,m)
               enddo
      
               do m=4,6
                  mx=m-4
                  f3(m)=effe(known,5,3,mx)
               enddo
   
               do m=7,9
                  mx=m-7
                  f3(m)=effe(known,8,3,mx)
               enddo
               t1 = f3(0)
               c3(0) = t1
               t2 = C0*C0
               t3 = t2*C0
               t4 = t3*t3
               t5 = f3(6)
               t6 = f3(4)
               t7 = t1-t6
               t8 = t2*t2
               t9 = t8*t8
               t10 = f3(1)
               t11 = t10*C0
               t12 = f3(2)
               t13 = f3(3)
               t14 = t13*C0
               t15 = t8*C0
               t15 = t15*t15
               t16 = f3(5)
               t17 = KK(3)
               t18 = MP12(3)
               t19 = (t11*t8+t12*t2+t14*t15+t7*t9-(t16*t8+t5))*t17
               t3 = (t3-(1.0_ki))*((1.0_ki)+t3)*((1.0_ki)+t4)
               t20 = t18*t3
               c3(1) = (-1.0_ki)*t19/t20
               t21 = t17*t17
               t22 = (t11+t12*t15+t14*t4+t7*t8-(t5*t9+t16))*t21
               t23 = t18*t18
               c3(2) = (t22/(t23*t3))
               t24 = (t7+t11*t9+t12*t4+t14*t2-(t5*t8+t16*t9))*t17*t21
               c3(3) = (-1.0_ki)*t24/(t20*t23)
               t7 = t7*C0
               t20 = t12*C0
               t23 = t16*C0
               t25 = t5*C0
               t26 = t10+t13*t4+t2*t7+t20*t9-(t25*t4+t15*t23)
               t27 = t17*t18*t3
               c3(4) = (t26/t27)
               t5 = (t6-t1)*t4+t15*t5+t16*t2-(t14*t9+t11*t2+t12)
               t6 = t17*t18
               t11 = t6*t6
               c3(5) = (t5/(t11*t3))
               t2 = t13+t10*t4+t2*t20+t7*t9-(t23*t8+t25)
               c3(6) = (t2/(t11*t27))
               t4 = t17*C1
               t7 = mu2g(3)
               c3(7) = (((f3(7)-t1)*t17*t21+t17*t2/t3*t4*t4*C1-t24/t3)/(t17*t21*t&
               &7))
               c3(8) = (-1.0_ki)*(f3(9)*t17-(t17/t3*t5*C1*C1+t19/t3))/(t18*t7)
               c3(9) = ((t22/t3-(t21*t26/t3*C1+f3(8)*t21))/(t21*t6*t7*C1))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	case(-1)
		do m=0,4
			f3(m)=effe(known,1,5,m)
		enddo
		do m=5,8
			mx=m-5
			f3(m)=effe(known,6,4,mx)
		enddo  
		do m=9,13
			mx=m-9
			f3(m)=effe(known,10,5,mx)
		enddo
			f3(14)=effe(known,15,1,0)
! CALCULATING C3-------------------------------------------------------------------!
      c3(0) = &
     &f3(0) 
      c3(1) = &
     &        (C0**11*f3(1) + C0**7*f3(2) + C0**3*f3(3) +  &
     &    C0**19*f3(4) + C0**15*(f3(0) - f3(5)) - C0**10*f3(6) -  &
     &    C0**5*f3(7) - f3(8))/((1 - C0**20)*MP12(3)) 
      c3(2) = &
     &        (-(C0**6*f3(1)) - C0**2*f3(2) - C0**18*f3(3) -  &
     &    C0**14*f3(4) + C0**10*(-f3(0) + f3(5)) + C0**5*f3(6) +  &
     &    f3(7) + C0**15*f3(8))/((1 - C0**20)*MP12(3)**2) 
      c3(3) = &
     &        (-f3(6) + C0*(f3(1) +  &
     &       C0**4*(f3(0) - f3(5) +  &
     &          C0**4*(C0**8*f3(2) + C0**4*f3(3) + f3(4) -  &
     &             C0**6*f3(7) - C0*f3(8)))))/ &
     &  ((1 - C0**20)*MP12(3)**3) 
      c3(4) = &
     &        (-f3(1) + C0**4*(-f3(0) + f3(5) +  &
     &       C0**4*(-(C0**8*f3(2)) - C0**4*f3(3) - f3(4) +  &
     &          C0**11*f3(6) + C0**6*f3(7) + C0*f3(8))))/ &
     &  ((1 - C0**20)*MP12(3))
      c3(5) = &
     &        (f3(2) + C0**3*(-f3(6) +  &
     &       C0*(f3(1) + C0**4* &
     &           (f3(0) - f3(5) +  &
     &             C0**4*(C0**4*f3(3) + f3(4) - C0**6*f3(7) -  &
     &                C0*f3(8))))))/((1 - C0**20)*MP12(3)**2) 
      c3(6) = &
     &        (-(C0**8*f3(1)) - C0**4*f3(2) - f3(3) - C0**16*f3(4) +  &
     &    C0**12*(-f3(0) + f3(5)) + C0**7*f3(6) + C0**2*f3(7) +  &
     &    C0**17*f3(8))/((1 - C0**20)*MP12(3)**3) 
      c3(7) = &
     &        ((C1**6*half*(C0**12*f3(1) + C0**8*f3(2) + C0**4*f3(3) +  &
     &         f3(4) + C0**16*(f3(0) - f3(5)) - C0**11*f3(6) -  &
     &         C0**6*f3(7) - C0*f3(8)))/(1 - C0**20) -  &
     &    (C1**5*half*(-(C0**8*f3(1)) - C0**4*f3(2) - f3(3) -  &
     &         C0**16*f3(4) + C0**12*(-f3(0) + f3(5)) +  &
     &         C0**7*f3(6) + C0**2*f3(7) + C0**17*f3(8)))/ &
     &     (1 - C0**20) + C1*Cm1* &
     &     ((half*(-f3(0) + f3(5) +  &
     &            C0**4*(-(C0**12*f3(1)) - C0**8*f3(2) -  &
     &               C0**4*f3(3) - f3(4) + C0**11*f3(6) +  &
     &               C0**6*f3(7) + C0*f3(8))))/(1 - C0**20) -  &
     &       half*f3(10)) +  &
     &    Cm1**2*(-((half*(-f3(6) +  &
     &              C0*(f3(1) +  &
     &                 C0**4* &
     &                  (f3(0) - f3(5) +  &
     &                    C0**4* &
     &                     (C0**8*f3(2) + C0**4*f3(3) + f3(4) -  &
     &                       C0**6*f3(7) - C0*f3(8))))))/ &
     &          (1 - C0**20)) - half*f3(11)) +  &
     &    C1**2*(half*f3(9) - half*f3(12) - half*f3(13) -  &
     &       half*f3(14) +  &
     &       MP12(3)*(-((C0**11*f3(1) + C0**7*f3(2) +  &
     &               C0**3*f3(3) + C0**19*f3(4) +  &
     &               C0**15*(f3(0) - f3(5)) - C0**10*f3(6) -  &
     &               C0**5*f3(7) - f3(8))/((1 - C0**20)*MP12(3))) &
     &           - (Cm1*(-f3(1) +  &
     &               C0**4* &
     &                (-f3(0) + f3(5) +  &
     &                  C0**4* &
     &                   (-(C0**8*f3(2)) - C0**4*f3(3) - f3(4) +  &
     &                     C0**11*f3(6) + C0**6*f3(7) + C0*f3(8))) &
     &               ))/((1 - C0**20)*MP12(3)) +  &
     &          ((-(C0**6*f3(1)) - C0**2*f3(2) - C0**18*f3(3) -  &
     &                C0**14*f3(4) + C0**10*(-f3(0) + f3(5)) +  &
     &                C0**5*f3(6) + f3(7) + C0**15*f3(8))/ &
     &              ((1 - C0**20)*MP12(3)**2) +  &
     &             (Cm1**2* &
     &                (f3(2) +  &
     &                  C0**3* &
     &                   (-f3(6) +  &
     &                     C0* &
     &                      (f3(1) +  &
     &                       C0**4* &
     &                       (f3(0) - f3(5) +  &
     &                       C0**4* &
     &                       (C0**4*f3(3) + f3(4) - C0**6*f3(7) -  &
     &                       C0*f3(8)))))))/ &
     &              ((1 - C0**20)*MP12(3)**2))*MP12(3) +  &
     &          (-((Cm1**3*half* &
     &                  (-(C0**8*f3(1)) - C0**4*f3(2) - f3(3) -  &
     &                    C0**16*f3(4) +  &
     &                    C0**12*(-f3(0) + f3(5)) + C0**7*f3(6) +  &
     &                    C0**2*f3(7) + C0**17*f3(8)))/ &
     &                ((1 - C0**20)*MP12(3)**3)) -  &
     &             (half*(-f3(6) +  &
     &                  C0* &
     &                   (f3(1) +  &
     &                     C0**4* &
     &                      (f3(0) - f3(5) +  &
     &                       C0**4* &
     &                       (C0**8*f3(2) + C0**4*f3(3) + f3(4) -  &
     &                       C0**6*f3(7) - C0*f3(8))))))/ &
     &              ((1 - C0**20)*MP12(3)**3))*MP12(3)**2 +  &
     &          ((Cm1**4*half* &
     &                (C0**12*f3(1) + C0**8*f3(2) + C0**4*f3(3) +  &
     &                  f3(4) + C0**16*(f3(0) - f3(5)) -  &
     &                  C0**11*f3(6) - C0**6*f3(7) - C0*f3(8)))/ &
     &              ((1 - C0**20)*MP12(3)**4) +  &
     &             (half*(-f3(0) + f3(5) +  &
     &                  C0**4* &
     &                   (-(C0**12*f3(1)) - C0**8*f3(2) -  &
     &                     C0**4*f3(3) - f3(4) + C0**11*f3(6) +  &
     &                     C0**6*f3(7) + C0*f3(8))))/ &
     &              ((1 - C0**20)*MP12(3)**4))*MP12(3)**3)))/ &
     &  (C1**2*mu2g(3)) 
      c3(8) = &
     &        -(((C0**11*f3(1) + C0**7*f3(2) + C0**3*f3(3) +  &
     &         C0**19*f3(4) + C0**15*(f3(0) - f3(5)) -  &
     &         C0**10*f3(6) - C0**5*f3(7) - f3(8))/ &
     &       ((1 - C0**20)*MP12(3)) -  &
     &      (C1**4*(C0**12*f3(1) + C0**8*f3(2) + C0**4*f3(3) +  &
     &           f3(4) + C0**16*(f3(0) - f3(5)) - C0**11*f3(6) -  &
     &           C0**6*f3(7) - C0*f3(8)))/((1 - C0**20)*MP12(3)) &
     &       + f3(13)/MP12(3))/mu2g(3)) 
      c3(9) = &
     &        ((-f3(0) + f3(5) +  &
     &       C0**4*(-(C0**12*f3(1)) - C0**8*f3(2) - C0**4*f3(3) -  &
     &          f3(4) + C0**11*f3(6) + C0**6*f3(7) + C0*f3(8)))/ &
     &     (1 - C0**20) - (C1* &
     &       (-f3(1) + C0**4* &
     &          (-f3(0) + f3(5) +  &
     &            C0**4*(-(C0**8*f3(2)) - C0**4*f3(3) - f3(4) +  &
     &               C0**11*f3(6) + C0**6*f3(7) + C0*f3(8)))))/ &
     &     (1 - C0**20) - f3(10))/(C1*MP12(3)*mu2g(3)) 
      c3(10) = &
     &        (-((-(C0**6*f3(1)) - C0**2*f3(2) - C0**18*f3(3) -  &
     &         C0**14*f3(4) + C0**10*(-f3(0) + f3(5)) +  &
     &         C0**5*f3(6) + f3(7) + C0**15*f3(8))/ &
     &       ((1 - C0**20)*MP12(3)**2)) +  &
     &    (C1**3*(-(C0**8*f3(1)) - C0**4*f3(2) - f3(3) -  &
     &         C0**16*f3(4) + C0**12*(-f3(0) + f3(5)) +  &
     &         C0**7*f3(6) + C0**2*f3(7) + C0**17*f3(8)))/ &
     &     ((1 - C0**20)*MP12(3)**2) + f3(12)/MP12(3)**2)/mu2g(3) 
      c3(11) = &
     &        (-((f3(2) + C0**3*(-f3(6) +  &
     &            C0*(f3(1) +  &
     &               C0**4* &
     &                (f3(0) - f3(5) +  &
     &                  C0**4* &
     &                   (C0**4*f3(3) + f3(4) - C0**6*f3(7) -  &
     &                     C0*f3(8))))))/((1 - C0**20)*MP12(3)**2) &
     &       ) + ((-f3(6) +  &
     &          C0*(f3(1) +  &
     &             C0**4*(f3(0) - f3(5) +  &
     &                C0**4* &
     &                 (C0**8*f3(2) + C0**4*f3(3) + f3(4) -  &
     &                   C0**6*f3(7) - C0*f3(8)))))/(1 - C0**20) &
     &        + f3(11))/(C1**2*MP12(3)**2))/mu2g(3) 
      c3(12) = &
     &        (-f3(0) + f3(5) + C0**4* &
     &     (-(C0**12*f3(1)) - C0**8*f3(2) - C0**4*f3(3) - f3(4) +  &
     &       C0**11*f3(6) + C0**6*f3(7) + C0*f3(8)))/ &
     &  ((1 - C0**20)*MP12(3)**4) 
      c3(13) = &
     &        (C0**12*f3(1) + C0**8*f3(2) + C0**4*f3(3) + f3(4) +  &
     &    C0**16*(f3(0) - f3(5)) - C0**11*f3(6) - C0**6*f3(7) -  &
     &    C0*f3(8))/((1 - C0**20)*MP12(3)**4) 
      c3(14) = &
     &        (-((C1**6*half*(C0**12*f3(1) + C0**8*f3(2) +  &
     &           C0**4*f3(3) + f3(4) + C0**16*(f3(0) - f3(5)) -  &
     &           C0**11*f3(6) - C0**6*f3(7) - C0*f3(8)))/ &
     &       (1 - C0**20)) +  &
     &    (C1**5*half*(-(C0**8*f3(1)) - C0**4*f3(2) - f3(3) -  &
     &         C0**16*f3(4) + C0**12*(-f3(0) + f3(5)) +  &
     &         C0**7*f3(6) + C0**2*f3(7) + C0**17*f3(8)))/ &
     &     (1 - C0**20) + C1*Cm1* &
     &     (-((half*(-f3(0) + f3(5) +  &
     &              C0**4*(-(C0**12*f3(1)) - C0**8*f3(2) -  &
     &                 C0**4*f3(3) - f3(4) + C0**11*f3(6) +  &
     &                 C0**6*f3(7) + C0*f3(8))))/(1 - C0**20)) +  &
     &       half*f3(10)) +  &
     &    Cm1**2*((half*(-f3(6) +  &
     &            C0*(f3(1) +  &
     &               C0**4* &
     &                (f3(0) - f3(5) +  &
     &                  C0**4* &
     &                   (C0**8*f3(2) + C0**4*f3(3) + f3(4) -  &
     &                     C0**6*f3(7) - C0*f3(8))))))/ &
     &        (1 - C0**20) + half*f3(11)) +  &
     &    C1**2*(-f3(0) + half*f3(9) + half*f3(12) +  &
     &       half*f3(13) + half*f3(14) +  &
     &       MP12(3)*((C0**11*f3(1) + C0**7*f3(2) + C0**3*f3(3) +  &
     &             C0**19*f3(4) + C0**15*(f3(0) - f3(5)) -  &
     &             C0**10*f3(6) - C0**5*f3(7) - f3(8))/ &
     &           ((1 - C0**20)*MP12(3)) +  &
     &          (Cm1*(-f3(1) +  &
     &               C0**4* &
     &                (-f3(0) + f3(5) +  &
     &                  C0**4* &
     &                   (-(C0**8*f3(2)) - C0**4*f3(3) - f3(4) +  &
     &                     C0**11*f3(6) + C0**6*f3(7) + C0*f3(8))) &
     &               ))/((1 - C0**20)*MP12(3)) +  &
     &          (-((-(C0**6*f3(1)) - C0**2*f3(2) - C0**18*f3(3) -  &
     &                  C0**14*f3(4) + C0**10*(-f3(0) + f3(5)) +  &
     &                  C0**5*f3(6) + f3(7) + C0**15*f3(8))/ &
     &                ((1 - C0**20)*MP12(3)**2)) -  &
     &             (Cm1**2* &
     &                (f3(2) +  &
     &                  C0**3* &
     &                   (-f3(6) +  &
     &                     C0* &
     &                      (f3(1) +  &
     &                       C0**4* &
     &                       (f3(0) - f3(5) +  &
     &                       C0**4* &
     &                       (C0**4*f3(3) + f3(4) - C0**6*f3(7) -  &
     &                       C0*f3(8)))))))/ &
     &              ((1 - C0**20)*MP12(3)**2))*MP12(3) +  &
     &          ((Cm1**3*half* &
     &                (-(C0**8*f3(1)) - C0**4*f3(2) - f3(3) -  &
     &                  C0**16*f3(4) + C0**12*(-f3(0) + f3(5)) +  &
     &                  C0**7*f3(6) + C0**2*f3(7) + C0**17*f3(8))) &
     &               /((1 - C0**20)*MP12(3)**3) +  &
     &             (half*(-f3(6) +  &
     &                  C0* &
     &                   (f3(1) +  &
     &                     C0**4* &
     &                      (f3(0) - f3(5) +  &
     &                       C0**4* &
     &                       (C0**8*f3(2) + C0**4*f3(3) + f3(4) -  &
     &                       C0**6*f3(7) - C0*f3(8))))))/ &
     &              ((1 - C0**20)*MP12(3)**3))*MP12(3)**2 +  &
     &          (-((Cm1**4*half* &
     &                  (C0**12*f3(1) + C0**8*f3(2) +  &
     &                    C0**4*f3(3) + f3(4) +  &
     &                    C0**16*(f3(0) - f3(5)) - C0**11*f3(6) -  &
     &                    C0**6*f3(7) - C0*f3(8)))/ &
     &                ((1 - C0**20)*MP12(3)**4)) -  &
     &             (half*(-f3(0) + f3(5) +  &
     &                  C0**4* &
     &                   (-(C0**12*f3(1)) - C0**8*f3(2) -  &
     &                     C0**4*f3(3) - f3(4) + C0**11*f3(6) +  &
     &                     C0**6*f3(7) + C0*f3(8))))/ &
     &              ((1 - C0**20)*MP12(3)**4))*MP12(3)**3)))/ &
     &  (C1**2*mu2g(3)**2) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	end select	
	endif
! LOWER RANK CS====================================================================!
	if (diff.ge.0) then
		c3(10)=czip
		c3(11)=czip
		c3(12)=czip
		c3(13)=czip
		c3(14)=czip
	if (diff.ge.1) then
		c3(3)=czip
		c3(6)=czip
		c3(8)=czip
		c3(9)=czip 
	if (diff.ge.2) then
		c3(2)=czip
		c3(5)=czip
		c3(7)=czip
	if (diff.ge.3) then
		c3(1)=czip
		c3(4)=czip
	endif
	endif 
	endif 
	endif
end subroutine getc3

subroutine HResL3(diff,resin,nsol,nloops,nleg,cut3,q3,qt,Vi,msq,mu2vec)
	use mfunctions
	use precision
	use mglobal, only: denst,resit,mu2t
! PARAMETERS============================================================================!
	! external parameters
		integer, intent(in) :: diff
		complex(ki), dimension(15),	    intent(inout) :: resin


		integer, 			    intent(in   ) :: nsol, nloops, nleg, cut3
		complex(ki), dimension(15,4), 	    intent(in   ) :: q3
		complex(ki), dimension(15),	    intent(in   ) :: mu2vec
		complex(ki), dimension(4), 	    intent(in   ) :: qt
		real(ki),    dimension(0:nleg-1,4), intent(in   ) :: Vi
		complex(ki),    dimension(0:nleg-1),   intent(in   ) :: msq
	! internal parameters
		integer			   :: i, j1, j2, j3, acc, n
		integer			   :: i1, i2, i3, i4, i5
		logical			   :: evalres
		complex(ki), dimension(15) :: densn
		complex(ki) 		   :: densnt
		logical 		   :: notcut
		integer,     dimension(5)  :: iarr
		integer 		   :: dicutn
! INITIALIZATION========================================================================!
	! cuts
		j3     =  cut3/100
		acc    =  j3*100
		j2     = (cut3-acc)/10
		j1     =  cut3-acc-j2*10
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
			if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
				do n=1,nsol
				densn(n) = densn(n)*denevalmu2(nleg,i,q3(n,:),Vi,msq,mu2vec(n))
				enddo
				densnt   = densnt*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
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
		resit(3)=resit(3)+densnt*Res5(dicutn,mu2t(3))
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
			if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
				do n=1,nsol
				densn(n) = densn(n)*denevalmu2(nleg,i,q3(n,:),Vi,msq,mu2vec(n))
				enddo
				densnt   = densnt*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
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
			resin(n)=resin(n)+densn(n)*Res4(diff,dicutn,q3(n,:),mu2vec(n))
		enddo
		resit(3)=resit(3)+densnt*Res4(diff,dicutn,qt,mu2t(3))	
	endif
	dicutn = dicutn+1
	enddo
	enddo
	enddo
	enddo
	end select
!========================================
end subroutine HResL3

end module mgetc3

