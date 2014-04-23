module mmishighrank
use precision
use constants
use options
implicit none

public :: HJ111

contains 

subroutine HJ111(B31, K11,m0,m1,B0p12,B0z11,B0z22)
        implicit none
! PARAMETERS=============================================================================!
        ! external parameters
        real(ki),                       intent(in   ) :: K11
        complex(ki),                    intent(in   ) :: m0,m1
        complex(ki), dimension(-2:0),   intent(in   ) :: B0p12,B0z11,B0z22
        complex(ki), dimension(-2:0),   intent(inout) :: B31
        ! internal parameters
        complex(ki)                     :: A1,A2,A3,A4,B1,B2,B3
        integer :: B31case
! INITIALIZATION=========================================================================!
        ! determine branch
        if(abs(K11) .gt. 1.0e-10) then
                B31case = 0
        elseif(m0 .ne. m1) then
                B31case = 1
        else
                B31case = 2
        endif
! CALCULATE B31(:)=======================================================================!
        select case(B31case)
        case(0)
                A1 = &
                     &        -(K11**3 - K11**2*m0 + K11*m0**2 - m0**3 + three*K11**2*m1 -  &
                     &     four*K11*m0*m1 + three*m0**2*m1 + three*K11*m1**2 - three*m0*m1**2 +  &
                     &     m1**3)/(four*K11**3)
                A2 = &
                     &        -(m0*(K11**2 - K11*m0 + m0**2 + two*K11*m1 - two*m0*m1 +  &
                     &       m1**2))/(four*K11**3)
                A3 = &
                     &        (m1*(three*K11**2 - two*K11*m0 + m0**2 + three*K11*m1 - two*m0*m1 +  &
                     &      m1**2))/(four*K11**3)
                A4 = &
                     &        -(two*K11**3 + 10.0_ki*K11**2*m0 - nine*K11*m0**2 + six*m0**3 -  &
                     &     10.0_ki*K11**2*m1 + 24.0_ki*K11*m0*m1 - 18.0_ki*m0**2*m1 -  &
                     &     15.0_ki*K11*m1**2 + 18.0_ki*m0*m1**2 - six*m1**3)/(24.0_ki*K11**3)
                B31(0) = A1*B0p12(0)+A2*B0z11(0)+A3*B0z22(0)+A4
                B31(-1) = A1*B0p12(-1)+A2*B0z11(-1)+A3*B0z22(-1)
                B31(-2) = A1*B0p12(-2)+A2*B0z11(-2)+A3*B0z22(-2)
        case(1)
                B1 = &
                     &-m0**4/(two*(m0 - m1)**4)
                B2 = &
                     &        (m1*(four*m0**3 - six*m0**2*m1 + four*m0*m1**2 - m1**3))/ &
                     &  (four*(m0 - m1)**4)
                B3 = &
                     &        -(25.0_ki*m0**3 - 23.0_ki*m0**2*m1 + 13.0_ki*m0*m1**2 - three*m1**3)/ &
                     &  (48.0_ki*(m0 - m1)**3)
                B31(0)  = B1*B0z11(0) + B2*B0z22(0) + B3
                B31(-1) = B1*B0z11(-1) + B2*B0z22(-1)
                B31(-2) = B1*B0z11(-2) + B2*B0z22(-2)
        case(2)
                B31(:) = -B0z11(:)/four
        end select
end subroutine HJ111

end module
