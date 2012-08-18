module save
   use precision, only: ki
   use constants, only: max1, max2, max3, max4, max5
   implicit none
   save

   private :: max1, max2, max3, max4, max5, ki


   real(ki) :: sav5p0(max5,4)
   real(ki) :: sav5e1(max5,4)
   real(ki) :: sav5e2(max5,4)
   complex(ki) :: sav5e3(max5,4)
   complex(ki) :: sav5e4(max5,4)
   complex(ki) :: savc5(max5)
   integer :: savcut5(max5)

   real(ki) :: savL3(max4,4)
   real(ki) :: sav4p0(max4,4)
   real(ki) :: sav4e1(max4,4)
   real(ki) :: sav4e2(max4,4)
   complex(ki) :: sav4e3(max4,4)
   complex(ki) :: sav4e4(max4,4)
   complex(ki) :: savc4(max4,0:4)
   integer :: savcut4(max4)

   real(ki) :: sav3p0(max3,4)
   real(ki) :: sav3e1(max3,4)
   real(ki) :: sav3e2(max3,4)
   complex(ki) :: sav3e3(max3,4)
   complex(ki) :: sav3e4(max3,4)
   complex(ki) :: savc3(max3,0:9)
   integer :: savcut3(max3)

   real(ki) :: sav2p0(max2,4)
   real(ki) :: sav2e1(max2,4)
   real(ki) :: sav2e2(max2,4)
   complex(ki) :: sav2e3(max2,4)
   complex(ki) :: sav2e4(max2,4)
   complex(ki) :: savc2(max2,0:9)
   integer :: savcut2(max2)

   real(ki) :: sav1p0(max1,4)
   real(ki) :: sav1e1(max1,4)
   real(ki) :: sav1e2(max1,4)
   complex(ki) :: sav1e3(max1,4)
   complex(ki) :: sav1e4(max1,4)
   complex(ki) :: savc1(max1,0:4)
   integer :: savcut1(max1)

end module save

