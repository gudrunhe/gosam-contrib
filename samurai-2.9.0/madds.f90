module madds
   use precision, only: ki
   use constants
   use options
   implicit none

   private

	public :: add4, add3, add2, add1

contains

subroutine add4(c4,MI4,tot4,tot4r)
	implicit none
	
	complex(ki), dimension(0:5),        intent(in ) :: c4
	complex(ki), dimension(-2:0),       intent(in ) :: MI4
	complex(ki), dimension(-2:0),       intent(out) :: tot4
	complex(ki),                        intent(out) :: tot4r
	
	tot4(:)=c4(0)*MI4(:)
	tot4r=-c4(4)/six
	tot4(0)=tot4(0) + tot4r
end subroutine add4



subroutine add3(V,m,c3,MI3,tot3,tot3r)
	implicit none
	
	real(ki),    dimension(3),     intent(in ) :: V
	complex(ki), dimension( 0:2 ), intent(in ) :: m
	complex(ki), dimension( 0:14), intent(in ) :: c3
	complex(ki), dimension(-2:0 ), intent(in ) :: MI3
	complex(ki), dimension(-2:0 ), intent(out) :: tot3
	complex(ki),                   intent(out) :: tot3r
	
	logical     :: highrank
	integer     :: ep

	highrank=(c3(14) /= 0._ki)
	
	tot3(:)=c3(0)*MI3(:)
	tot3r=+c3(7)/two
	
	if (highrank) then
		tot3(0) = tot3(0) + c3(14) * ((V(1)+V(2)+V(3))/four-m(0)-m(1)-m(2))/six
	end if
	tot3(0)=tot3(0) + tot3r

end subroutine add3

subroutine add2(K11,K12,m,c2,MI2a,MI2b,MI2c,MI2d,MI2e,tot2,tot2r,MI2J111)
	implicit none
	
	real(ki),                     intent(in )           :: K11
	real(ki),                     intent(in )           :: K12
	complex(ki), dimension(0:1),  intent(in )           :: m
	complex(ki), dimension(0:19), intent(in )           :: c2
	complex(ki), dimension(-2:0), intent(in )           :: MI2a, MI2b, MI2c, MI2d, MI2e 
	complex(ki), dimension(-2:0), intent(out)           :: tot2
	complex(ki),                  intent(out)           :: tot2r
	complex(ki), dimension(-2:0), intent(in ), optional :: MI2J111 
	
	complex(ki) :: m0, m1, B06
	integer     :: ep
	logical     :: highrank
	complex(ki), dimension(-2:0) :: MI2J0, MI2J1, MI2J00, MI2J01, MI2J11 
! Added by HvD
!------------

	m0=m(0)
	m1=m(1)
	
	highrank = ((c2(13) /= 0._ki) .or. (c2(10) /= 0._ki))
	
	
	B06=-(K11-three*(m0+m1))/six
	
	tot2r = B06*c2(9)
	
	if (isca.eq.1 .or. isca.eq.3 .or. isca.eq.4) then
	MI2J0(:)  = MI2a(:)
	MI2J1(:)  = MI2b(:)
	MI2J00(:) = MI2c(:)
	MI2J01(:) = MI2d(:)
	MI2J11(:) = MI2e(:)
		if (abs(K11).gt.zip1) then
			do ep=-2,0
				tot2(ep)=-(K12*(two*K12*(m0 - m1)*c2(2) +  &
				& K11*(-three*c2(1) + two*K12*c2(2)))*MI2J0(ep))/(six*K11**2) &
				& + ((two*K12**2*(m0 - m1)**2*c2(2) +  &
				& K11*K12*(-three*m0*c2(1) + three*m1*c2(1) +  &
				& two*K12*m0*c2(2) - four*K12*m1*c2(2)) +  &
				& K11**2*(six*c2(0) + K12*(-three*c2(1) + two*K12*c2(2))))* &
				& MI2J01(ep))/(six*K11**2) +  &
				& (K12*(two*K12*(m0 - m1)*c2(2) +  &
				& K11*(-three*c2(1) + four*K12*c2(2)))*MI2J1(ep))/(six*K11**2)
			enddo
			tot2(0)=tot2(0)+(K12**2*c2(2))/18.0_ki &
			& - (K12**2*m0*c2(2))/(six*K11) -  &
			& (K12**2*m1*c2(2))/(six*K11) + B06*c2(9)
		else
			if (m1.eq.m0) then
				do ep=-2,0
					tot2(ep)=(c2(0) + (K12*(-three*c2(1) &
					& + two*K12*c2(2)))/six)*MI2J00(ep)
				enddo
				tot2(0)=tot2(0)+B06*c2(9)
			else
			        do ep=-2,0
					tot2(ep)=(K12*m0**2*(-three*m0*c2(1) + three*m1*c2(1) &
					& + two*K12*m0*c2(2))* &
					& MI2J00(ep))/(six*(m0 - m1)**3) + c2(0)*MI2J01(ep) - &
					& (K12*m1*(m0*m1*(nine*c2(1) - six*K12*c2(2)) - &
					& six*m0**2*(c2(1) - K12*c2(2)) +  &
					& m1**2*(-three*c2(1) + two*K12*c2(2)))*MI2J11(ep))/ &
					& (six*(m0 - m1)**3)
				enddo
				tot2(0)=tot2(0) + (-three*K12*m0*c2(1))/(four*(m0 - m1)) + &
				& (K12*m1*c2(1))/(four*m0 - four*m1) + &
				& (11.0_ki*K12**2*m0**2*c2(2))/(18.0_ki*(m0 - m1)**2) - &
				& (7.0_ki*K12**2*m0*m1*c2(2))/(18.0_ki*(m0 - m1)**2) + &
				& (K12**2*m1**2*c2(2))/(9.0_ki*(m0 - m1)**2) + B06*c2(9)
			 endif
		endif
	elseif (isca.eq.2) then
		tot2(:) = c2(0)*MI2a(:)+c2(1)*K12*MI2b(:)+c2(2)*K12*K12*MI2c(:)
		tot2(0) = tot2(0) + c2(9)*B06
	else
		print *,'error in madds, add2: isca value not allowed'
	end if
	if (highrank) then
		tot2(0) = tot2(0) + c2(13)*(K12**3)*MI2J111(0) + c2(10)*one/12.0_ki* K12 * &
		& (K11 - two*m0-four*m1)
		tot2(-1) = tot2(-1) + c2(13)*(K12**3)*MI2J111(-1)
		tot2(-2) = tot2(-2) + c2(13)*(K12**3)*MI2J111(-2)
	end if
	
 903   format(a5,I2,a1,I2,a4,2(D24.15))

end subroutine add2



subroutine add1(m0,e3e4,c1,MI1,tot1)

	implicit none 
	complex(ki),                  intent(in ) :: m0
	complex(ki),                  intent(in ) :: e3e4
	complex(ki), dimension(0:15), intent(in ) :: c1
	complex(ki), dimension(-2:0), intent(in ) :: MI1
	complex(ki), dimension(-2:0), intent(out) :: tot1

	tot1(:)=c1(0)*MI1(:)

	if ((c1(14) /= 0) .or. (c1(15) /= 0) ) then
		tot1(0) = tot1(0) + (c1(14) + e3e4 * c1(15)/four)*m0*m0/two
	end if

end subroutine add1


end module madds

