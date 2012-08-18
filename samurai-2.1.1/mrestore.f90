module mrestore
   use precision
   use constants
   use save
   use mfunctions
   implicit none

   private

   public :: store1, store2, store3, store4, store5
   public :: res1, res2, res3, res4, res5

contains
  subroutine store5(icut5,cut5,p0,e1,e2,e3,e4,c5)
      implicit none
      integer, intent(in) :: icut5,cut5
      real(ki), dimension(4), intent(in) :: p0, e1, e2
      complex(ki), intent(in) :: c5
      complex(ki), dimension(4), intent(in) ::  e3, e4
      sav5p0(icut5,:)=p0(:)
      sav5e1(icut5,:)=e1(:)
      sav5e2(icut5,:)=e2(:)
      sav5e3(icut5,:)=e3(:)
      sav5e4(icut5,:)=e4(:)
      savc5(icut5)=c5
      savcut5(icut5)=cut5
  end subroutine store5

  subroutine store4(icut4,cut4,L3,p0,e1,e2,e3,e4,c4)
      implicit none
      real(ki), dimension(4), intent(in) :: L3, p0, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:4), intent(in) :: c4
      integer, intent(in) :: icut4,cut4
      savL3(icut4,:)=L3(:)
      sav4p0(icut4,:)=p0(:)
      sav4e1(icut4,:)=e1(:)
      sav4e2(icut4,:)=e2(:)
      sav4e3(icut4,:)=e3(:)
      sav4e4(icut4,:)=e4(:)
      savc4(icut4,:)=c4(:)
      savcut4(icut4)=cut4
  end subroutine store4

  subroutine store3(icut3,cut3,p0,e1,e2,e3,e4,c3)
      implicit none
      real(ki), dimension(4), intent(in) :: p0, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:9), intent(in) :: c3
      integer, intent(in) :: icut3,cut3
      sav3p0(icut3,:)=p0(:)
      sav3e1(icut3,:)=e1(:)
      sav3e2(icut3,:)=e2(:)
      sav3e3(icut3,:)=e3(:)
      sav3e4(icut3,:)=e4(:)
      savc3(icut3,:)=c3(:)
      savcut3(icut3)=cut3
  end subroutine store3

  subroutine store2(icut2,cut2,p0,e1,e2,e3,e4,c2)
      implicit none
      real(ki), dimension(4), intent(in) :: p0, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:9), intent(in) :: c2
      integer, intent(in) :: icut2,cut2
      sav2p0(icut2,:)=p0(:)
      sav2e1(icut2,:)=e1(:)
      sav2e2(icut2,:)=e2(:)
      sav2e3(icut2,:)=e3(:)
      sav2e4(icut2,:)=e4(:)
      savc2(icut2,:)=c2(:)
      savcut2(icut2)=cut2
  end subroutine store2

  subroutine store1(icut1,cut1,p0,e1,e2,e3,e4,c1)
      implicit none
      real(ki), dimension(4), intent(in) :: p0, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:4), intent(in) :: c1
      integer, intent(in) :: icut1,cut1

      sav1p0(icut1,:)=p0(:)
      sav1e1(icut1,:)=e1(:)
      sav1e2(icut1,:)=e2(:)
      sav1e3(icut1,:)=e3(:)
      sav1e4(icut1,:)=e4(:)
      savc1(icut1,:)=c1(:)
      savcut1(icut1)=cut1
  end subroutine store1

   pure function res5(icut5,mu2)
      implicit none
      integer, intent(in) :: icut5
      complex(ki), intent(in) :: mu2
      complex(ki) :: res5

      res5=savc5(icut5)*mu2
      ! res5=savc5(icut5)*mu2**2
   end  function res5

   pure function res4(icut4,q,mu2)
      implicit none
      integer, intent(in) :: icut4
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(4), intent(in) :: q
      complex(ki) :: res4

      real(ki), dimension(4) :: L3
      complex(ki), dimension(4) :: e3, e4, pm
      complex(ki), dimension(0:4) :: c4

      L3(:)=savL3(icut4,:)
      pm(:)=sav4p0(icut4,:)+q(:)
      e3(:)=sav4e3(icut4,:)
      e4(:)=sav4e4(icut4,:)
      c4(:)=savc4(icut4,:)

      res4=poly4(c4,pm,mu2,L3,e3,e4)
   end  function res4

   pure function res3(icut3,q,mu2)
      implicit none
      integer, intent(in) :: icut3
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(4), intent(in) :: q
      complex(ki) :: res3

      complex(ki), dimension(4) :: e3, e4, pm
      complex(ki) :: c3(0:9)

      pm(:)=sav3p0(icut3,:)+q(:)
      e3(:)=sav3e3(icut3,:)
      e4(:)=sav3e4(icut3,:)
      c3(:)=savc3(icut3,:)

      res3=poly3(c3,pm,mu2,e3,e4)
   end  function res3

   pure function res2(icut2,q,mu2)
      implicit none
      integer, intent(in) :: icut2
      complex(ki), dimension(4), intent(in) :: q
      complex(ki), intent(in) :: mu2
      complex(ki) :: res2

      real(ki), dimension(4) :: e2
      complex(ki), dimension(4) :: e3,e4,pm
      complex(ki), dimension(0:9) :: c2

      pm(:)=sav2p0(icut2,:)+q(:)
      e2(:)=sav2e2(icut2,:)
      e3(:)=sav2e3(icut2,:)
      e4(:)=sav2e4(icut2,:)
      c2(:)=savc2(icut2,:)
      res2=poly2(c2,pm,mu2,e2,e3,e4)
   end  function res2

   pure function res1(icut1,q)
      implicit none
      integer, intent(in) :: icut1
      complex(ki), dimension(4), intent(in) :: q
      complex(ki) :: res1

      real(ki), dimension(4) ::  e1, e2
      complex(ki), dimension(4) :: e3, e4, pm
      complex(ki), dimension(0:4) :: c1

      pm(:)=sav1p0(icut1,:)+q(:)
      e1(:)=sav1e1(icut1,:)
      e2(:)=sav1e2(icut1,:)
      e3(:)=sav1e3(icut1,:)
      e4(:)=sav1e4(icut1,:)
      c1(:)=savc1(icut1,:)

      res1=poly1(c1,pm,e1,e2,e3,e4)
   end function res1

end module mrestore
