module mgetbase
   use precision
   use constants
   use mfunctions, only: sdot
   use kinematic, only: epsi
   implicit none

   private

   public :: getbase

contains

      subroutine getbase(p1,p2,r1,r2,e1,e2,e3,e4)
      implicit none

      real(ki), dimension(4), intent(in) :: p1, p2
      real(ki), dimension(4), intent(out) :: e1, e2

      real(ki) :: gamma, MP12, MP11, MP22, den
      real(ki), intent(out) :: r1, r2
      complex(ki), intent(out) :: e3(4), e4(4)

      MP11=sdot(p1,p1)
      MP12=sdot(p1,p2)
      MP22=sdot(p2,p2)

      gamma=MP12+sign(+1.0_ki,MP12)*sqrt(MP12**2-MP11*MP22)

!--- Warning: if p1 ~ p2 and massless this formula brakes down
      if ( abs(gamma) .lt. zip1 ) then
         write(6,*) 'getbase: small gamma'
         write(6,*) 'gamma=', gamma
         write(6,*) ' MP11=', MP11
         write(6,*) ' MP12=', MP12
         write(6,*) ' MP22=', MP22
      endif

      r1=MP11/gamma
      r2=MP22/gamma
      den=1.0_ki-r1*r2

      e1(:)=(p1(:)-r1*p2(:))/den
      e2(:)=(p2(:)-r2*p1(:))/den
      e3 = sqrt(abs(sdot(e1,e2))) * epsi(+1, e1, e2)
      e4 = conjg(e3)*sign(1.0_ki,e1(4)*e2(4))
  end subroutine getbase
end module mgetbase

