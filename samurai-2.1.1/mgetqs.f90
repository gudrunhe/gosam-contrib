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

   subroutine getq5_rm(nleg,cut5,e1,e2,e3,e4,p0,Vi,msq,r1,r2,q5,mu2)
      implicit none
      integer, intent(in) :: nleg, cut5
      real(ki), dimension(4), intent(in) :: e1, e2, p0
      complex(ki), dimension(4), intent(in) :: e3, e4
      real(ki), intent(in) :: r1, r2
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      complex(ki), dimension(4), intent(out) :: q5
      complex(ki), intent(out) :: mu2

      integer :: j1,j2,j3,j4,j5
      real(ki) :: x1, x2,var2, var3
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

      mu2=real(tmu2, ki)
      q5(:)=-p0(:)+x1*e1(:)+x2*e2(:)+x3*e3(:)+x4*e4(:)
   end subroutine getq5_rm

   subroutine getq4_cm(nleg,cut4,e1,e2,e3,e4,p0,k1,k2,k3,L3,r1,r2,q4,qt,msq)
      use mglobal, only: MP12, mu2g, mu2t, dx
      implicit none
      integer, intent(in) :: nleg, cut4
      real(ki), dimension(4) :: k1, k2, k3, L3, e1, e2, p0
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), intent(in) :: r1,r2
      complex(ki), dimension(5,4), intent(out) :: q4
      complex(ki), dimension(4), intent(out) :: qt

      real(ki) :: beta,mu2x
      complex(ki) :: x1,x2,A1,A2
      real(ki) :: ME31, ME32, MK11, MK22, MK33
      complex(ki) :: x31,x32,x33,x34,x35,x3t,x41,x42,x43,x44,x45,x4t
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

   subroutine getq4_rm(nleg,cut4,e1,e2,e3,e4,p0,k1,k2,k3,L3,r1,r2,q4,qt,msq)
      use mglobal, only: MP12, mu2g, mu2t, dx
      implicit none
      integer, intent(in) :: nleg, cut4
      real(ki), dimension(4) :: k1, k2, k3, L3, e1, e2, p0
      complex(ki), dimension(4), intent(in) :: e3, e4
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), intent(in) :: r1,r2
      complex(ki), dimension(5,4), intent(out) :: q4
      complex(ki), dimension(4), intent(out) :: qt

      real(ki) :: beta,mu2x,x1,x2,A1,A2
      real(ki) :: ME31, ME32, MK11, MK22, MK33
      complex(ki) :: x31,x32,x33,x34,x35,x3t,x41,x42,x43,x44,x45,x4t
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
      mu2g(4)=max(msq(j1),msq(j2),msq(j3),msq(j4),&
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
   end subroutine getq4_rm

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

   subroutine getq3_rm(nleg,irank,cut3,e1,e2,e3,e4,p0,k1,k2,msq,r1,r2,q3,qt)
      use mglobal, only: C0,C1,mu2g,MP12,KK,Kmu,mu2t
      use options, only: C0_thrs
      implicit none
      integer, intent(in) :: nleg, cut3, irank
      real(ki), dimension(4), intent(in) :: p0, k1, k2, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(10,4), intent(out) :: q3
      complex(ki), dimension(4), intent(out) :: qt

      integer :: j1,j2,j3,j,ndiff
      real(ki) :: r1,r2,beta,teta,x1,x2,A1,A2,C1t
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

      C0=x1*x2-msq(j1)/two/MP12(3)

      !---#[ scelta dinamica:
      mu2g(3)=max(msq(j1),msq(j2),msq(j3),abs(MK11),abs(MK22),abs(two*MK12))
      if (abs(mu2g(3)).lt.1d-10) mu2g(3)=one
      !---#] scelta dinamica:
      !---#[ scelta statica:
      ! mu2g(3)=one
      !---#] scelta statica:

      C1=x1*x2-(msq(j1)+mu2g(3))/two/MP12(3)

      if (abs(C0-1.0_ki) .lt. C0_thrs) then
         !---#[ New Sampling:
         ! The new sampling is the one that is safe around C0=1
         ! but not around C0=0

         ! complete c-system: 10 coefficients

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
               x4(j)=C0*cone/x3(j)
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
            q3(4,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C0*x3(4)*e3(:)+x4(4)*e4(:)
            q3(5,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C0*x3(5)*e3(:)+x4(5)*e4(:)
            q3(6,:)=-p0(:)+x1*e1(:)+x2*e2(:)+C1*x3(6)*e3(:)+x4(6)*e4(:)
         case(2)
            ! rank1 c-system: 3 coefficients
       
            x3(1)=cone
            x4(1)=C0*cone
        
            x3(2)=-cone
            x4(2)=-C0*cone
         
            x3(3)=C0*cone
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
   end subroutine getq3_rm

   subroutine getq2_cm(nleg,irank,cut2,e1,e2,e3,e4,p0,k1,msq,q2,qt)
      use mglobal, only: Fpc,Fzc,Fmc,F1zc,KB,mu2g,KK,Kmu,MP12,mu2t
      implicit none
      integer, intent(in) :: nleg, cut2, irank
      real(ki), dimension(4), intent(in) :: p0, k1, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(10,4), intent(out) :: q2
      complex(ki), dimension(4), intent(out) :: qt

      complex(ki), dimension(10) :: x1, x2
      real(ki) :: teta,K1SQ
      complex(ki) :: A1,x1t,x2t,Ft
      complex(ki) :: X1p,X2p,X1z,X2z,X1m,X2m
      complex(ki) :: x3t, x4t
      complex(ki), dimension(10) :: x3, x4
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

!--- scelta dinamica

      mu2g(2)=max(abs(msq(j1)),abs(msq(j2)),K1SQ)
      if (abs(mu2g(2)).lt.1.0e-10_ki) mu2g(2)=one
 
!--- scelta statica
!      mu2g(2)=one

      F1zc=-(msq(j1)+mu2g(2))/two/MP12(2)

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

      teta=twopi/three
      do j=1,3
      x1(j)=X1z
      x2(j)=X2z
      x4(j)=cos(teta*real(j-1,ki))+im*sin(teta*real(j-1,ki))
      x3(j)=Fzc/x4(j)
      enddo
      teta=twopi/two
      do j=4,5
      x1(j)=X1z
      x2(j)=X2z
      x3(j)=cos(teta*real(j-4,ki))-im*sin(teta*real(j-4,ki))
      x4(j)=Fzc/x3(j)
      enddo

      teta=twopi/two
      do j=6,7
      x1(j)=X1p
      x2(j)=X2p
      x4(j)=cos(teta*real(j-6,ki))+im*sin(teta*real(j-6,ki))
      x3(j)=Fpc/x4(j)
      enddo
      x1(8)=X1p
      x2(8)=X2p
      x3(8)=cone
      x4(8)=Fpc

      x1(9)=X1m
      x2(9)=X2m
      x3(9)=Fmc
      x4(9)=cone

      x1(10)=X1z
      x2(10)=X2z
      x3(10)=F1zc
      x4(10)=cone

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
      q2(10,:)=-p0(:)+x1(10)*e1(:)+x2(10)*e2(:)+KK(2)*x3(10)*e3(:) &
     & +x4(10)*e4(:)/KK(2)
      qt(:)=-p0(:)+x1t*e1(:)+x2t*e2(:)+x3t*e3(:)+x4t*e4(:)
      
   endif

   end subroutine getq2_cm

   subroutine getq2_rm(nleg,irank,cut2,e1,e2,e3,e4,p0,k1,msq,q2,qt)
      use mglobal, only: Fp,Fz,Fm,F1z,KB,mu2g,KK,Kmu,MP12,mu2t
      implicit none
      integer, intent(in) :: nleg, cut2, irank
      real(ki), dimension(4), intent(in) :: p0, k1, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(10,4), intent(out) :: q2
      complex(ki), dimension(4), intent(out) :: qt

      real(ki), dimension(10) :: x1, x2
      real(ki) :: teta,A1,K1SQ,x1t,x2t,Ft
      real(ki) :: X1p,X2p,X1z,X2z,X1m,X2m
      complex(ki) :: x3t, x4t
      complex(ki), dimension(10) :: x3, x4
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
      Fz=-msq(j1)/two/MP12(2)

      X1p=KB
      X2p=A1-KB*K1SQ/two/MP12(2)
      Fp=KB*X2p-msq(j1)/two/MP12(2)

      X1m=-KB
      X2m=A1+KB*K1SQ/two/MP12(2)
      Fm=-KB*X2m-msq(j1)/two/MP12(2)

!--- scelta dinamica
      mu2g(2)=max(msq(j1),msq(j2),K1SQ)
      if (abs(mu2g(2)).lt.1.0e-10_ki) mu2g(2)=one
 
!--- scelta statica
!      mu2g(2)=one

      F1z=-(msq(j1)+mu2g(2))/two/MP12(2)

!!!! HERE WE BRANCH ACCORDING TO THE RANK

      if (ndiff.eq.1) then
! rank 1
          
         teta=twopi/two

         do j=1,2
            x1(j)=X1z
            x2(j)=X2z
            x4(j)=cos(teta*real(j-1,ki))+im*sin(teta*real(j-1,ki))
            x3(j)=Fz/x4(j)
         enddo

         do j=3,4
            x1(j)=X1p
            x2(j)=X2p
            x3(j)=cos(teta*real(j-3,ki))-im*sin(teta*real(j-3,ki))
            x4(j)=Fp/x3(j)
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

      teta=twopi/three
      do j=1,3
      x1(j)=X1z
      x2(j)=X2z
      x4(j)=cos(teta*real(j-1,ki))+im*sin(teta*real(j-1,ki))
      x3(j)=Fz/x4(j)
      enddo
      teta=twopi/two
      do j=4,5
      x1(j)=X1z
      x2(j)=X2z
      x3(j)=cos(teta*real(j-4,ki))-im*sin(teta*real(j-4,ki))
      x4(j)=Fz/x3(j)
      enddo

      teta=twopi/two
      do j=6,7
      x1(j)=X1p
      x2(j)=X2p
      x4(j)=cos(teta*real(j-6,ki))+im*sin(teta*real(j-6,ki))
      x3(j)=Fp/x4(j)
      enddo
      x1(8)=X1p
      x2(8)=X2p
      x3(8)=cone
      x4(8)=Fp*cone

      x1(9)=X1m
      x2(9)=X2m
      x3(9)=Fm*cone
      x4(9)=cone

      x1(10)=X1z
      x2(10)=X2z
      x3(10)=F1z*cone
      x4(10)=cone

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
      q2(10,:)=-p0(:)+x1(10)*e1(:)+x2(10)*e2(:)+KK(2)*x3(10)*e3(:) &
     & +x4(10)*e4(:)/KK(2)
      qt(:)=-p0(:)+x1t*e1(:)+x2t*e2(:)+x3t*e3(:)+x4t*e4(:)
      
   endif

   end subroutine getq2_rm

   subroutine getq1_cm(nleg,cut1,e1,e2,e3,e4,p0,msq,q1,qt)
     use mglobal, only: G0c,KK,mu2g,MP12,mu2t
      implicit none
      integer, intent(in) :: nleg, cut1
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

   subroutine getq1_rm(nleg,cut1,e1,e2,e3,e4,p0,msq,q1,qt)
     use mglobal, only: G0,KK,mu2g,MP12,mu2t
      implicit none
      integer, intent(in) :: nleg, cut1
      real(ki), dimension(4), intent(in) :: p0, e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4
      complex(ki), dimension(5,4), intent(out) :: q1
      complex(ki), dimension(4), intent(out) :: qt
      real(ki), dimension(0:nleg-1), intent(in) :: msq

      real(ki), dimension(5) :: x1, x2
      complex(ki), dimension(5) :: x3, x4
      real(ki) :: varx,x1t,x2t,x3t,x4t
      integer :: j1

      KK(1)=one
      j1=cut1

      mu2g(1)=zip
      if (msq(j1).eq.zip) then
      mu2g(1)=one
      endif

      MP12(1)=sdot(e1,e2)
      G0=(msq(j1)+mu2g(1))/two/MP12(1)

      x1(1)=G0
      x2(1)=one
      x3(1)=zip
      x4(1)=zip

      x1(2)=-G0
      x2(2)=-one
      x3(2)=zip
      x4(2)=zip

      x1(3)=G0
      x2(3)=one
      x3(3)=one
      x4(3)=zip

      x1(4)=zip
      x2(4)=one
      x3(4)=-one
      x4(4)=G0

      x1(5)=zip
      x2(5)=one
      x3(5)=one
      x4(5)=-G0


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
  end subroutine getq1_rm

end module mgetqs

