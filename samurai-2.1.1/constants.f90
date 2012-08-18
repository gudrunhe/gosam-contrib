module constants
use precision
implicit none

      private :: ki

      real(ki), parameter :: pi = 3.14159265358979323846&
      &2643383279502884197169399375105820974944592307816&
      &4062862089986280348253421170679821480865132823066&
      &47093844609550582231725359408128_ki
      real(ki),parameter :: twopi = 2.0_ki*pi
!-----------------------------------------------------
      real(ki),parameter :: zip=0.0_ki
      real(ki),parameter :: half=0.5_ki
      real(ki),parameter :: one=1.0_ki
      real(ki),parameter :: two=2.0_ki
      real(ki),parameter :: three=3.0_ki
      real(ki),parameter :: four=4.0_ki
      real(ki),parameter :: six=6.0_ki
      real(ki),parameter :: seven=7.0_ki
      real(ki),parameter :: eight=8.0_ki
      real(ki),parameter :: rt2 = 1.414213562373095048801&
     &688724209698078569671875376948073176679737990732478&
     &4621070388503875343276415727350138462309122970249248361_ki
      real(ki),parameter :: rt3 = 1.732050807568877293527&
     &446341505872366942805253810380628055806979451933016&
     &9088000370811461867572485756756261414154067030299699451_ki
!-----------------------------------------------------
      complex(ki),parameter :: chaf=(0.5_ki,0.0_ki)
      complex(ki),parameter :: im=(0.0_ki,1.0_ki)
      complex(ki),parameter :: impi =(0.0_ki,pi)
      complex(ki),parameter :: czip=(0.0_ki,0.0_ki)
      complex(ki),parameter :: cone=(1.0_ki,0.0_ki)
      complex(ki),parameter :: ctwo=(2.0_ki,0.0_ki)
!-----------------------------------------------------
!                 for nleg max =         ! 8! 7! 6! 5! 4!
!                                        !--------------!
      integer,parameter :: max5=  56     !56!21! 6! 1! 0!
      integer,parameter :: max4=  70     !70!35!15! 5! 1!
      integer,parameter :: max3=  56     !56!35!20!10! 4!
      integer,parameter :: max2=  28     !28!21!15!10! 6!
      integer,parameter :: max1=   8     ! 8! 7! 6! 5! 4!
      integer,parameter :: maxleg= 8     ! 8! 7! 6! 5! 4!
!                                        !--------------!              
      integer,parameter :: nq5=1
      integer,parameter :: nq4=5
      integer,parameter :: nq3=10
      integer,parameter :: nq2=10
      integer,parameter :: nq1=5

!-----------------------------------------------------
      real(ki),parameter ::      zip1 = 1.0e-14_ki

end module constants


