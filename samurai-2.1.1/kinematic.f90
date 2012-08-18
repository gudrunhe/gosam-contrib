module     kinematic
   use precision, only: ki
   implicit none
   save

   private

   interface dotproduct
      module procedure dotproduct_rr
      module procedure dotproduct_rc
      module procedure dotproduct_cr
      module procedure dotproduct_cc
   end interface dotproduct


   interface zb
      module procedure zb_rr
      module procedure zb_rc
      module procedure zb_cr
      module procedure zb_cc
   end interface zb

   interface za
      module procedure za_rr
      module procedure za_rc
      module procedure za_cr
      module procedure za_cc
   end interface za

   interface zab
      module procedure zab_rcr
      module procedure zab_ccc
      module procedure zab_ccr
      module procedure zab_rcc
      module procedure zab_rrr
      module procedure zab_crc
      module procedure zab_crr
      module procedure zab_rrc
   end interface zab

   interface zba
      module procedure zba_rcr
      module procedure zba_ccc
      module procedure zba_ccr
      module procedure zba_rcc
      module procedure zba_rrr
      module procedure zba_crc
      module procedure zba_crr
      module procedure zba_rrc
   end interface zba

   interface zbb
      module procedure zbb_rcrr
   end interface zbb

   public :: dotproduct, zb, za, zab, zbb, zba
   public :: inspect_Vi, epsi, epso, gamma_6

contains

   subroutine inspect_Vi(Vi)
      implicit none
      real(ki), dimension(8,4), intent(in) :: Vi
      integer :: k, i

      do k = 1,8
      do i = 1,4
       write(*,'(A3,I1,A1,I1,A4,F24.15)') "Vi(", k, ",", i, ") = ", &
     & Vi(k, i)
      end do
      end do
   end subroutine


! blocco epsi epso

   pure function epsi(pol, vec, aux)
      implicit none
      integer, intent(in) :: pol
      real(ki), dimension(4), intent(in) :: vec, aux
      complex(ki), dimension(4) :: epsi

      real(ki), dimension(4) :: e1, e2, p, k
      real(ki) :: r, s, n
      integer, dimension(1) :: dir

      p(1:3) = aux(1:3) / aux(4)
      p(4) = 1.0_ki

      k(1:3) = vec(1:3) / vec(4)
      k(4) = 1.0_ki

      r = 1.0_ki - dotproduct(k, p)
      s = sign(1.0_ki, real(pol, ki))
              
      e1(4) = 0.0_ki
      if (abs(r + 1.0_ki) > 1.0E+03_ki * epsilon(1.0_ki)) then
              n = 1.0_ki / sqrt(2.0_ki * abs(1.0_ki - r*r))

              e1(1) = k(2) * p(3) - k(3) * p(2)
              e1(2) = k(3) * p(1) - k(1) * p(3)
              e1(3) = k(1) * p(2) - k(2) * p(1)
              e1(1:3) = n * e1(1:3)

              e2(1:3) = k(1:3) + p(1:3)
              e2(4) = 1.0_ki + r
              e2(:) = n * e2(:)
      else
              dir = minloc(abs(k(1:3)))
              if (dir(1) .eq. 1) then
                      n = sqrt(0.5_ki / (k(2)*k(2) + k(3)*k(3)))
                      e1(1) = 0.0_ki
                      e1(2) = - k(3) * n
                      e1(3) =   k(2) * n
              elseif (dir(1) .eq. 2) then
                      n = sqrt(0.5_ki / (k(1)*k(1) + k(3)*k(3)))
                      e1(2) = 0.0_ki
                      e1(3) = - k(1) * n
                      e1(1) =   k(3) * n
              else
                      n = sqrt(0.5_ki / (k(1)*k(1) + k(2)*k(2)))
                      e1(3) = 0.0_ki
                      e1(1) = - k(2) * n
                      e1(2) =   k(1) * n
              endif

              e2(1) = k(2) * e1(3) - k(3) * e1(2)
              e2(2) = k(3) * e1(1) - k(1) * e1(3)
              e2(3) = k(1) * e1(2) - k(2) * e1(1)
              e2(4) = 0.0_ki
      end if
      epsi(:) = e1(:) + cmplx(0.0_ki, s, ki) * e2(:)
   end  function epsi

   pure function epso(pol, vec, aux)
      implicit none
      integer, intent(in) :: pol
      real(ki), dimension(4), intent(in) :: vec, aux
      complex(ki), dimension(4) :: epso

      epso = conjg(epsi(pol, vec, aux))
   end  function epso


! blocco dotproduct

   pure function dotproduct_rr(p, q)
      implicit none
      real(ki), dimension(4), intent(in) :: p, q
      real(ki) :: dotproduct_rr
      dotproduct_rr = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
   end  function dotproduct_rr

   pure function dotproduct_cc(p, q)
      implicit none
      complex(ki), dimension(4), intent(in) :: p, q
      complex(ki) :: dotproduct_cc
      dotproduct_cc = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
   end  function dotproduct_cc

   pure function dotproduct_rc(p, q)
      implicit none
      real(ki), dimension(4), intent(in) :: p
      complex(ki), dimension(4), intent(in) :: q
      complex(ki) :: dotproduct_rc
      dotproduct_rc = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
   end  function dotproduct_rc

   pure function dotproduct_cr(p, q)
      implicit none
      complex(ki), dimension(4), intent(in) :: p
      real(ki), dimension(4), intent(in) :: q
      complex(ki) :: dotproduct_cr
      dotproduct_cr = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
   end  function dotproduct_cr


! blocco zab

   pure function zab_ccc(p1, k, p2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1
      complex(ki), dimension(4), intent(in) :: p2
      complex(ki), dimension(4), intent(in) :: k
      complex(ki) :: zab_ccc

      complex(ki) :: kp, km, im
      complex(ki) :: kr, kl, pr1, pr2, pl1, pl2, rt1, rt2

      im=(0.0_ki,1.0_ki)

      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=+k(3)-im*k(2)
      kl=+k(3)+im*k(2)

      rt1=sqrt((p1(4)+p1(1)))
      rt2=sqrt((p2(4)+p2(1)))
      pr1=p1(3)-im*p1(2)
      pr2=p2(3)-im*p2(2)
      pl1=p1(3)+im*p1(2)
      pl2=p2(3)+im*p2(2)

      zab_ccc=&
     &   (+pr1*pl2*kp/(rt1*rt2)&
     &    -pr1*kl*rt2/rt1&
     &    -rt1/rt2*kr*pl2+rt1*rt2*km)

   end  function zab_ccc


   pure function zab_ccr(p1, k, k2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1
      real(ki), dimension(4), intent(in) :: k2
      complex(ki), dimension(4), intent(in) :: k
      complex(ki) :: zab_ccr

      complex(ki) :: kp, km, im
      complex(ki) :: kr, kl, pr1, pr2, pl1, pl2, rt1, f2
      real(ki) :: flip2, rt2

      im=(0.0_ki,1.0_ki)

      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=+k(3)-im*k(2)
      kl=+k(3)+im*k(2)

      rt1=sqrt((p1(4)+p1(1)))
      pr1=p1(3)-im*p1(2)
      pl1=p1(3)+im*p1(2)

      if (k2(4) .gt. 0.0_ki) then
         flip2=1.0_ki
         f2=(1.0_ki, 0.0_ki)
      else
         flip2=-1.0_ki
         f2=(0.0_ki, 1.0_ki)
      endif
      rt2=sqrt(flip2*(k2(4)+k2(1)))
      pr2=cmplx(flip2*k2(3),-flip2*k2(2), ki)
      pl2=conjg(pr2)

      zab_ccr=f2*&
     &   (+pr1*pl2*kp/(rt1*rt2)&
     &    -pr1*kl*rt2/rt1&
     &    -rt1/rt2*kr*pl2+rt1*rt2*km)

   end  function zab_ccr

   pure function zab_rcc(k1, k, p2)
      implicit none
      real(ki), dimension(4), intent(in) :: k1
      complex(ki), dimension(4), intent(in) :: p2
      complex(ki), dimension(4), intent(in) :: k
      complex(ki) :: zab_rcc

      complex(ki) :: kp, km, im
      complex(ki) :: kr, kl, pr1, pr2, pl1, pl2, rt2, f1
      real(ki) :: flip1, rt1

      im=(0.0_ki,1.0_ki)

      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=+k(3)-im*k(2)
      kl=+k(3)+im*k(2)

      if (k1(4) .gt. 0.0_ki) then
         flip1=1.0_ki
         f1=1.0_ki
      else
         flip1=-1.0_ki
         f1=(0.0_ki, 1.0_ki)
      endif
      rt1=sqrt(flip1*(k1(4)+k1(1)))
      pr1=cmplx(flip1*k1(3),-flip1*k1(2), ki)
      pl1=conjg(pr1)

      rt2=sqrt((p2(4)+p2(1)))
      pr2=p2(3)-im*p2(2)
      pl2=p2(3)+im*p2(2)

      zab_rcc=f1*&
     &   (+pr1*pl2*kp/(rt1*rt2)&
     &    -pr1*kl*rt2/rt1&
     &    -rt1/rt2*kr*pl2+rt1*rt2*km)

   end  function zab_rcc



   pure function zab_crc(p1, k, p2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1
      complex(ki), dimension(4), intent(in) :: p2
      real(ki), dimension(4), intent(in) :: k
      complex(ki) :: zab_crc

      real(ki) :: kp, km
      complex(ki) :: kr, kl, pr1, pr2, pl1, pl2, rt1, rt2, im

      im=(0.0_ki,1.0_ki)

      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=+k(3)-im*k(2)
      kl=+k(3)+im*k(2)

      rt1=sqrt((p1(4)+p1(1)))
      rt2=sqrt((p2(4)+p2(1)))
      pr1=p1(3)-im*p1(2)
      pr2=p2(3)-im*p2(2)
      pl1=p1(3)+im*p1(2)
      pl2=p2(3)+im*p2(2)

      zab_crc=&
     &   (+pr1*pl2*kp/(rt1*rt2)&
     &    -pr1*kl*rt2/rt1&
     &    -rt1/rt2*kr*pl2+rt1*rt2*km)

   end  function zab_crc



   pure   function zab_rrr(p1, Q, p2)
      implicit none
      real(ki), dimension(4), intent(in) :: p1, p2, Q
      real(ki), dimension(4) :: q1
      real(ki) :: r2
      complex(ki) :: zab_rrr

!---- decomposing Q along p2

      r2=dotproduct(Q,Q)/(2.0_ki*dotproduct(Q,p2))

      q1(:) = Q(:) - r2*p2(:)

      zab_rrr = za(p1,q1)*zb(q1,p2)
   end  function zab_rrr


   pure   function zabccc(p1, Q, p2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1, p2, Q
      complex(ki), dimension(4) :: q1
      complex(ki) :: r2
      complex(ki) :: zabccc

!---- decomposing Q along p2

      r2=dotproduct(Q,Q)/(2.0_ki*dotproduct(Q,p2))

      q1(:) = Q(:) - r2*p2(:)

      zabccc = za(p1,q1)*zb(q1,p2)
   end  function zabccc



   pure function PMzab_rcr(p1, Q, p2)
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      real(ki), dimension(4), intent(in) :: p1, p2
      complex(ki), dimension(4) :: q1
      complex(ki) :: r2
      complex(ki) :: PMzab_rcr

!---- decomposing Q along p2

      r2=dotproduct(Q,Q)/(2.0_ki*dotproduct(Q,p2))

      q1(:) = Q(:) - r2*p2(:)

      PMzab_rcr = za(p1,q1)*zb(q1,p2)
   end  function PMzab_rcr



   pure   function zab_rcr(k1, k, k2)
      implicit none
      real(ki), dimension(4), intent(in) :: k1
      real(ki), dimension(4), intent(in) :: k2
      complex(ki), dimension(4), intent(in) :: k
      complex(ki) :: zab_rcr
      complex(ki) :: kp, km
      complex(ki) :: kr, kl
      complex(ki) :: pr1, pr2, pl1, pl2, im, f1, f2
      real(ki) :: flip1, flip2, rt1, rt2

      im=(0.0_ki,1.0_ki)

      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=+k(3)-im*k(2)
      kl=+k(3)+im*k(2)

      if (k1(4) .gt. 0.0_ki) then
         flip1=1.0_ki
         f1=1.0_ki
      else
         flip1=-1.0_ki
         f1=(0.0_ki, 1.0_ki)
      endif
      rt1=sqrt(flip1*(k1(4)+k1(1)))
      pr1=cmplx(flip1*k1(3),-flip1*k1(2), ki)
      pl1=conjg(pr1)

      if (k2(4) .gt. 0.0_ki) then
         flip2=1.0_ki
         f2=(1.0_ki, 0.0_ki)
      else
         flip2=-1.0_ki
         f2=(0.0_ki, 1.0_ki)
      endif
      rt2=sqrt(flip2*(k2(4)+k2(1)))
      pr2=cmplx(flip2*k2(3),-flip2*k2(2), ki)
      pl2=conjg(pr2)

      zab_rcr=f1*f2*(&
     &   (+pr1*pl2*kp/(rt1*rt2)&
     &    -pr1*kl*rt2/rt1&
     &    -rt1/rt2*kr*pl2+rt1*rt2*km))

   end  function zab_rcr


   pure function zab_crr(p1, k, k2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1
      real(ki), dimension(4), intent(in) :: k2
      real(ki), dimension(4), intent(in) :: k
      complex(ki) :: zab_crr

      real(ki) :: kp, km
      complex(ki) :: kr, kl, pr1, pr2, pl1, pl2, rt1, f2, im
      real(ki) :: flip2, rt2

      im=(0.0_ki,1.0_ki)

      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=+k(3)-im*k(2)
      kl=+k(3)+im*k(2)

      rt1=sqrt((p1(4)+p1(1)))
      pr1=p1(3)-im*p1(2)
      pl1=p1(3)+im*p1(2)

      if (k2(4) .gt. 0.0_ki) then
         flip2=1.0_ki
         f2=(1.0_ki, 0.0_ki)
      else
         flip2=-1.0_ki
         f2=(0.0_ki, 1.0_ki)
      endif
      rt2=sqrt(flip2*(k2(4)+k2(1)))
      pr2=cmplx(flip2*k2(3),-flip2*k2(2), ki)
      pl2=conjg(pr2)

      zab_crr=f2*&
     &   (+pr1*pl2*kp/(rt1*rt2)&
     &    -pr1*kl*rt2/rt1&
     &    -rt1/rt2*kr*pl2+rt1*rt2*km)

   end  function zab_crr

   pure function zab_rrc(k1, k, p2)
      implicit none
      real(ki), dimension(4), intent(in) :: k1
      complex(ki), dimension(4), intent(in) :: p2
      real(ki), dimension(4), intent(in) :: k
      complex(ki) :: zab_rrc

      real(ki) :: kp, km
      complex(ki) :: kr, kl, pr1, pr2, pl1, pl2, rt2, im, f1
      real(ki) :: flip1, rt1

      im=(0.0_ki,1.0_ki)

      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=+k(3)-im*k(2)
      kl=+k(3)+im*k(2)

      if (k1(4) .gt. 0.0_ki) then
         flip1=1.0_ki
         f1=1.0_ki
      else
         flip1=-1.0_ki
         f1=(0.0_ki, 1.0_ki)
      endif
      rt1=sqrt(flip1*(k1(4)+k1(1)))
      pr1=cmplx(flip1*k1(3),-flip1*k1(2), ki)
      pl1=conjg(pr1)

      rt2=sqrt((p2(4)+p2(1)))
      pr2=p2(3)-im*p2(2)
      pl2=p2(3)+im*p2(2)

      zab_rrc=f1*&
     &   (+pr1*pl2*kp/(rt1*rt2)&
     &    -pr1*kl*rt2/rt1&
     &    -rt1/rt2*kr*pl2+rt1*rt2*km)

   end  function zab_rrc


! blocco zba

   pure function zba_ccc(p,k,q) result(zba)
      implicit none
      complex(ki), dimension(4), intent(in) :: p,k,q
      complex(ki) :: zba
      zba = zab(q,k,p)
   end  function zba_ccc

   pure function zba_rcc(p,k,q) result(zba)
      implicit none
      real(ki), dimension(4), intent(in) :: p
      complex(ki), dimension(4), intent(in) :: k,q
      complex(ki) :: zba
      zba = zab(q,k,p)
   end  function zba_rcc

   pure function zba_crc(p,k,q) result(zba)
      implicit none
      real(ki), dimension(4), intent(in) :: k
      complex(ki), dimension(4), intent(in) :: p,q
      complex(ki) :: zba
      zba = zab(q,k,p)
   end  function zba_crc

   pure function zba_ccr(p,k,q) result(zba)
      implicit none
      real(ki), dimension(4), intent(in) :: q
      complex(ki), dimension(4), intent(in) :: p,k
      complex(ki) :: zba
      zba = zab(q,k,p)
   end  function zba_ccr

   pure function zba_rrc(p,k,q) result(zba)
      implicit none
      real(ki), dimension(4), intent(in) :: p,k
      complex(ki), dimension(4), intent(in) :: q
      complex(ki) :: zba
      zba = zab(q,k,p)
   end  function zba_rrc

   pure function zba_rcr(p,k,q) result(zba)
      implicit none
      real(ki), dimension(4), intent(in) :: p,q
      complex(ki), dimension(4), intent(in) :: k
      complex(ki) :: zba
      zba = zab(q,k,p)
   end  function zba_rcr

   pure function zba_crr(p,k,q) result(zba)
      implicit none
      real(ki), dimension(4), intent(in) :: k,q
      complex(ki), dimension(4), intent(in) :: p
      complex(ki) :: zba
      zba = zab(q,k,p)
   end  function zba_crr

   pure function zba_rrr(p,k,q) result(zba)
      implicit none
      real(ki), dimension(4), intent(in) :: p,k,q
      complex(ki) :: zba
      zba = zab(q,k,p)
   end  function zba_rrr

! blocco zb

   pure function zb_rr(p, q)
      implicit none
      real(ki), dimension(4), intent(in) :: p, q
      complex(ki) :: zb_rr
      zb_rr = 2.0_ki*dotproduct(p, q)/za(q, p)
   end  function zb_rr


! blocco za

   pure function za_rr(k1, k2)
      implicit none
      real(ki), dimension(4), intent(in) :: k1, k2
      complex(ki) :: za_rr

      real(ki) :: rt1, rt2
      complex(ki) :: c231, c232, f1, f2
!-----positive energy case
         if (k1(4) .gt. 0.0_ki) then
            rt1=sqrt(k1(4)+k1(1))
            c231=cmplx(k1(3),-k1(2), ki)
            f1=1.0_ki
         else
!-----negative energy case
            rt1=sqrt(-k1(4)-k1(1))
            c231=cmplx(-k1(3),k1(2), ki)
            f1=(0.0_ki, 1.0_ki)
         endif
!-----positive energy case
         if (k2(4) .gt. 0.0_ki) then
            rt2=sqrt(k2(4)+k2(1))
            c232=cmplx(k2(3),-k2(2), ki)
            f2=1.0_ki
         else
!-----negative energy case
            rt2=sqrt(-k2(4)-k2(1))
            c232=cmplx(-k2(3),k2(2), ki)
            f2=(0.0_ki, 1.0_ki)
         endif

         za_rr = -f2*f1*(c232*rt1/rt2-c231*rt2/rt1)

   end  function za_rr

   pure function za_cc(p1, p2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1, p2
      complex(ki) :: za_cc

      complex(ki) :: rt1, rt2, im
      complex(ki) :: c231, c232

      im=(0.0_ki,1.0_ki)

      rt1=sqrt((p1(4)+p1(1)))
      c231=p1(3)-im*p1(2)

      rt2=sqrt((p2(4)+p2(1)))
      c232=p2(3)-im*p2(2)

      za_cc = -(c232*rt1/rt2-c231*rt2/rt1)

   end  function za_cc


   pure function zb_cc(p1, p2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1, p2
      complex(ki) :: zb_cc

      complex(ki) :: rt1, rt2, im
      complex(ki) :: c231, c232

      im=(0.0_ki,1.0_ki)

      rt1=sqrt((p1(4)+p1(1)))
      c231=p1(3)+im*p1(2)

      rt2=sqrt((p2(4)+p2(1)))
      c232=p2(3)+im*p2(2)

      zb_cc = c232*rt1/rt2-c231*rt2/rt1

   end  function zb_cc


   pure function za_rc(k1, p2)
      implicit none
      real(ki), dimension(4), intent(in) :: k1
      complex(ki), dimension(4), intent(in) :: p2
      complex(ki) :: za_rc
      real(ki) :: rt1
      complex(ki) :: c231, c232, f1, im, rt2

      im=(0.0_ki,1.0_ki)

!-----positive energy case
         if (k1(4) .gt. 0.0_ki) then
            rt1=sqrt(k1(4)+k1(1))
            c231=cmplx(k1(3),-k1(2), ki)
            f1=1.0_ki
         else
!-----negative energy case
            rt1=sqrt(-k1(4)-k1(1))
            c231=cmplx(-k1(3),k1(2), ki)
            f1=(0.0_ki, 1.0_ki)
         endif

      rt2=sqrt((p2(4)+p2(1)))
      c232=p2(3)-im*p2(2)

      za_rc = -f1*(c232*rt1/rt2-c231*rt2/rt1)

   end  function za_rc


   pure function zb_rc(k1, p2)
      implicit none
      real(ki), dimension(4), intent(in) :: k1
      complex(ki), dimension(4), intent(in) :: p2
      complex(ki) :: zb_rc
      real(ki) :: rt1
      complex(ki) :: c231, c232, f1, im, rt2

      im=(0.0_ki,1.0_ki)

!-----positive energy case
         if (k1(4) .gt. 0.0_ki) then
            rt1=sqrt(k1(4)+k1(1))
            c231=cmplx(k1(3),k1(2), ki)
            f1=1.0_ki
         else
!-----negative energy case
            rt1=sqrt(-k1(4)-k1(1))
            c231=cmplx(-k1(3),-k1(2), ki)
            f1=(0.0_ki, 1.0_ki)
         endif

      rt2=sqrt((p2(4)+p2(1)))
      c232=p2(3)+im*p2(2)

      zb_rc = f1*(c232*rt1/rt2-c231*rt2/rt1)

   end  function zb_rc


   pure function za_cr(p1, k2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1
      real(ki), dimension(4), intent(in) :: k2
      complex(ki) :: za_cr

      real(ki) :: rt2
      complex(ki) :: c231, c232, f2, rt1, im

      im=(0.0_ki,1.0_ki)
      rt1=sqrt((p1(4)+p1(1)))
      c231=p1(3)-im*p1(2)

!-----positive energy case
         if (k2(4) .gt. 0.0_ki) then
            rt2=sqrt(k2(4)+k2(1))
            c232=cmplx(k2(3),-k2(2), ki)
            f2=1.0_ki
         else
!-----negative energy case
            rt2=sqrt(-k2(4)-k2(1))
            c232=cmplx(-k2(3),k2(2), ki)
            f2=(0.0_ki, 1.0_ki)
         endif

         za_cr = -f2*(c232*rt1/rt2-c231*rt2/rt1)

   end  function za_cr

   pure function zb_cr(p1, k2)
      implicit none
      complex(ki), dimension(4), intent(in) :: p1
      real(ki), dimension(4), intent(in) :: k2
      complex(ki) :: zb_cr

      real(ki) :: rt2
      complex(ki) :: c231, c232, f2, rt1, im

      im=(0.0_ki,1.0_ki)
      rt1=sqrt((p1(4)+p1(1)))
      c231=p1(3)+im*p1(2)

!-----positive energy case
         if (k2(4) .gt. 0.0_ki) then
            rt2=sqrt(k2(4)+k2(1))
            c232=cmplx(k2(3),k2(2), ki)
            f2=1.0_ki
         else
!-----negative energy case
            rt2=sqrt(-k2(4)-k2(1))
            c232=cmplx(-k2(3),-k2(2), ki)
            f2=(0.0_ki, 1.0_ki)
         endif

         zb_cr = f2*(c232*rt1/rt2-c231*rt2/rt1)

   end  function zb_cr


   function zbb_rcrr(p1,P, Q,p2)
      implicit none
      complex(ki), dimension(4), intent(in) :: P
      real(ki), dimension(4), intent(in) :: p1, p2, Q
      real(ki), dimension(4) :: q1
      real(ki) :: r2
      complex(ki) :: zbb_rcrr

!---- decomposing Q along p2

      r2=dotproduct(Q,Q)/(2.0_ki*dotproduct(Q,p2))

      q1(:) = Q(:) - r2*p2(:)

      zbb_rcrr = -zb(p2,q1)*zab(q1,P,p1)
   end  function zbb_rcrr

  function gamma_6(a,b,c,d,e,f)
  implicit none
  complex(ki), dimension(4), intent(in) :: a,b,c,d,e,f
  complex(ki) gamma_6
     gamma_6 = 4.0_ki*(&
    &  + dotproduct(a,b)*dotproduct(c,d)*dotproduct(e,f)&
    &  - dotproduct(a,b)*dotproduct(c,e)*dotproduct(d,f)&
    &  + dotproduct(a,b)*dotproduct(c,f)*dotproduct(d,e)&
    &  - dotproduct(a,c)*dotproduct(b,d)*dotproduct(e,f)&
    &  + dotproduct(a,c)*dotproduct(b,e)*dotproduct(d,f)&
    &  - dotproduct(a,c)*dotproduct(b,f)*dotproduct(d,e)&
    &  + dotproduct(a,d)*dotproduct(b,c)*dotproduct(e,f)&
    &  - dotproduct(a,d)*dotproduct(b,e)*dotproduct(c,f)&
    &  + dotproduct(a,d)*dotproduct(b,f)*dotproduct(c,e)&
    &  - dotproduct(a,e)*dotproduct(b,c)*dotproduct(d,f)&
    &  + dotproduct(a,e)*dotproduct(b,d)*dotproduct(c,f)&
    &  - dotproduct(a,e)*dotproduct(b,f)*dotproduct(c,d)&
    &  + dotproduct(a,f)*dotproduct(b,c)*dotproduct(d,e)&
    &  - dotproduct(a,f)*dotproduct(b,d)*dotproduct(c,e)&
    &  + dotproduct(a,f)*dotproduct(b,e)*dotproduct(c,d))
   end function gamma_6


end module kinematic


