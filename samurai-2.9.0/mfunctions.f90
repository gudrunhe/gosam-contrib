module     mfunctions
   use precision
   use constants
   implicit none

   private

   interface sdot
      module procedure sdot_rr
      module procedure sdot_rc
      module procedure sdot_cr
      module procedure sdot_cc
   end interface sdot

   interface denevalmu2
      module procedure denevalmu2_rr
      module procedure denevalmu2_cr
      module procedure denevalmu2_rm
      module procedure denevalmu2_cm
   end interface denevalmu2

   public :: sdot, denevalmu2, effe,effe2, poly1, poly2, poly3, poly4
contains

   pure function sdot_rr(p, q)
      implicit none
      real(ki), dimension(4), intent(in) :: p, q
      real(ki) :: sdot_rr
      sdot_rr = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
   end  function sdot_rr

   pure function sdot_cc(p, q)
      implicit none
      complex(ki), dimension(4), intent(in) :: p, q
      complex(ki) :: sdot_cc
      sdot_cc = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
   end  function sdot_cc

   pure function sdot_rc(p, q)
      implicit none
      real(ki), dimension(4), intent(in) :: p
      complex(ki), dimension(4), intent(in) :: q
      complex(ki) :: sdot_rc
      sdot_rc = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
   end  function sdot_rc

   pure function sdot_cr(p, q)
      implicit none
      complex(ki), dimension(4), intent(in) :: p
      real(ki), dimension(4), intent(in) :: q
      complex(ki) :: sdot_cr
      sdot_cr = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
   end  function sdot_cr

   pure function denevalmu2_cm(nleg,j,q,Vi,msq,mu2)
      implicit none
      integer, intent(in) :: nleg, j
      integer :: i
      complex(ki), intent(in) :: q(4)
      complex(ki) :: L(4), denevalmu2_cm
      complex(ki), intent(in) :: mu2
      complex(ki), intent(in) :: msq(0:nleg-1)
      real(ki), intent(in) :: Vi(0:nleg-1,4)
      do i=1,4
         L(i)=q(i)+Vi(j,i)*cone
      enddo
      denevalmu2_cm=sdot(L,L)-msq(j)-mu2
   end  function denevalmu2_cm

   pure function denevalmu2_rm(nleg,j,q,Vi,msq,mu2)
      implicit none
      integer, intent(in) :: nleg, j
      integer :: i
      complex(ki), intent(in) :: q(4)
      complex(ki) :: L(4), denevalmu2_rm
      complex(ki), intent(in) :: mu2
      real(ki), intent(in) :: msq(0:nleg-1)
      real(ki), intent(in) :: Vi(0:nleg-1,4)
      do i=1,4
         L(i)=q(i)+Vi(j,i)*cone
      enddo
      denevalmu2_rm=sdot(L,L)-msq(j)-mu2*cone
   end  function denevalmu2_rm

   pure function denevalmu2_cr(nleg,j,q,Vi,msq,mu2)
      implicit none
      integer, intent(in) :: nleg, j
      integer :: i
      complex(ki), intent(in) :: q(4)
      complex(ki) :: L(4), denevalmu2_cr
      real(ki), intent(in) :: mu2
      complex(ki), intent(in) :: msq(0:nleg-1)
      real(ki), intent(in) :: Vi(0:nleg-1,4)
      do i=1,4
         L(i)=q(i)+Vi(j,i)*cone
      enddo
      denevalmu2_cr=sdot(L,L)-msq(j)-mu2
   end  function denevalmu2_cr

   pure function denevalmu2_rr(nleg,j,q,Vi,msq,mu2)
      implicit none
      integer, intent(in) :: nleg, j
      integer :: i
      complex(ki), intent(in) :: q(4)
      complex(ki) :: L(4), denevalmu2_rr
      real(ki), intent(in) :: mu2
      real(ki), intent(in) :: msq(0:nleg-1)
      real(ki), intent(in) :: Vi(0:nleg-1,4)
      do i=1,4
         L(i)=q(i)+Vi(j,i)*cone
      enddo
      denevalmu2_rr=sdot(L,L)-msq(j)-mu2*cone
   end  function denevalmu2_rr

   pure function effe(known,nk,ns,m)
      implicit none
      complex(ki),intent(in) :: known(15)
      complex(ki) :: temp, effe
      real(ki) :: teta
      integer, intent(in) :: nk,ns,m
      integer :: k
      teta=twopi/ns
      temp=czip
      do k=0,ns-1
      temp=temp+known(k+nk)*&
     &    (cos(teta*real(m*k, ki))+im*sin(teta*real(m*k, ki)))
      enddo
      effe=temp/ns
   end function effe

   pure function effe2(known,nk,ns,m)
      implicit none
      complex(ki),intent(in) :: known(20)
      complex(ki) :: temp, effe2
      real(ki) :: teta
      integer, intent(in) :: nk,ns,m
      integer :: k
      teta=twopi/ns
      temp=czip
      do k=0,ns-1
      temp=temp+known(k+nk)*&
     &    (cos(teta*real(m*k, ki))+im*sin(teta*real(m*k, ki)))
      enddo
      effe2=temp/ns
   end function effe2


   pure function poly4(diff,c4,pm,mu2,L3,e3,e4)
      implicit none
	integer, intent(in) :: diff
      complex(ki), dimension(0:5), intent(in) :: c4
      complex(ki), dimension(4), intent(in) :: pm, e3, e4
      real(ki), dimension(4), intent(in) :: L3
      complex(ki), intent(in) :: mu2
      complex(ki) :: poly4

	if(diff.ge.4) then
	      poly4=c4(0) 
	elseif(diff.eq.3) then
	      poly4=c4(0) &
	     &     +c4(1) &
	     &      *(+sdot(pm,e3)*sdot(L3,e4) &
	     &        -sdot(pm,e4)*sdot(L3,e3))
	elseif(diff.eq.2) then
	      poly4=c4(0) &
	     &     +mu2*c4(2) &
	     &     +c4(1) &
	     &      *(+sdot(pm,e3)*sdot(L3,e4) &
	     &        -sdot(pm,e4)*sdot(L3,e3))
	elseif(diff.eq.1) then
	      poly4=c4(0) &
	     &     +mu2*c4(2) &
	     &     +(c4(1)+c4(3)*mu2) &
	     &      *(+sdot(pm,e3)*sdot(L3,e4) &
	     &        -sdot(pm,e4)*sdot(L3,e3))
	elseif(diff.eq.0) then
	      poly4=c4(0) &
	     &     +mu2*(c4(2)+c4(4)*mu2) &
	     &     +(c4(1)+c4(3)*mu2) &
	     &      *(+sdot(pm,e3)*sdot(L3,e4) &
	     &        -sdot(pm,e4)*sdot(L3,e3))
	else
	      poly4=c4(0) &
	     &     +mu2*(c4(2)+c4(4)*mu2) &
	     &     +(c4(1)+c4(3)*mu2+c4(5)*mu2**2) &
	     &      *(+sdot(pm,e3)*sdot(L3,e4) &
	     &        -sdot(pm,e4)*sdot(L3,e3))
	endif
   end  function poly4

   pure function poly3(diff,c3,pm,mu2,e3,e4)
      implicit none
	integer, intent(in) :: diff
      complex(ki), dimension(0:14), intent(in) :: c3
      complex(ki), dimension(4), intent(in) :: pm, e3, e4
      complex(ki), intent(in) :: mu2
      complex(ki) :: poly3

      complex(ki) :: pme3, pme4

      pme3=sdot(pm,e3)
      pme4=sdot(pm,e4)

	if(diff.ge.3) then
	      poly3=+c3(0) 
	elseif(diff.eq.2) then
	      poly3=+c3(0) &
	     &      +pme3*c3(1) &
	     &      +pme4*c3(4) 
	elseif(diff.eq.1) then
	      poly3=+c3(0) &
	     &      +pme3*(c3(1)+pme3*c3(2)) &
	     &      +pme4*(c3(4)+pme4*c3(5)) &
	     &      +mu2*c3(7)
	elseif(diff.eq.0) then
	      poly3=+c3(0) &
	     &      +pme3*(c3(1)+pme3*(c3(2)+c3(3)*pme3)) &
	     &      +pme4*(c3(4)+pme4*(c3(5)+c3(6)*pme4)) &
	     &      +mu2*(c3(7)+c3(8)*pme3+c3(9)*pme4)
	else
	      poly3=+c3(0) &
	     &      +pme3*(c3(1)+pme3*(c3(2)+c3(3)*pme3)) &
	     &      +pme4*(c3(4)+pme4*(c3(5)+c3(6)*pme4)) &
	     &      +mu2*( c3(7)+c3(8)*pme3+c3(9)*pme4     &
	     &            +c3(10)*pme3**2+c3(11)*pme4**2)       &
	     &      +c3(12)*pme3**4+c3(13)*pme4**4        &
	     &      +c3(14)*mu2**2
	endif
   end  function poly3

   pure function poly2(diff,c2,pm,mu2,e2,e3,e4)
      implicit none
	integer, intent(in) :: diff
      complex(ki), dimension(0:19), intent(in) :: c2
      real(ki), dimension(4), intent(in) :: e2
      complex(ki), dimension(4), intent(in) :: e3,e4,pm
      complex(ki), intent(in) :: mu2
      complex(ki) :: poly2

      complex(ki) :: pme2, pme3, pme4

      pme2=sdot(e2,pm)
      pme3=sdot(pm,e3)
      pme4=sdot(pm,e4)
	if(diff.ge.2) then
	      poly2=+c2(0) 
	elseif(diff.eq.1) then
	      poly2=+c2(0) &
	     &      +pme2*c2(1) &
	     &      +pme3*c2(3) &
	     &      +pme4*c2(5) 
	elseif(diff.eq.0) then
	      poly2=+c2(0) &
	     &      +pme2*(c2(1)+c2(2)*pme2+c2(7)*pme3+c2(8)*pme4) &
	     &      +pme3*(c2(3)+c2(4)*pme3) &
	     &      +pme4*(c2(5)+c2(6)*pme4) &
	     &      +c2(9)*mu2
	else
	      poly2=+c2(0) &
	     &      +pme2*(c2(1)+c2(2)*pme2+c2(7)*pme3+c2(8)*pme4) &
	     &      +pme3*(c2(3)+c2(4)*pme3) &
	     &      +pme4*(c2(5)+c2(6)*pme4) &
	     &      +c2(9)*mu2 &
	     &      +mu2*(c2(10)*pme2+c2(11)*pme3+c2(12)*pme4) &
	     &      +c2(13)*pme2**3+c2(14)*pme3**3+c2(15)*pme4**3 &
	     &      +pme2**2*(c2(16)*pme3+c2(17)*pme4) &
	     &      +pme2*(c2(18)*pme3**2+c2(19)*pme4**2)
	endif

  end function poly2


   pure function poly1(diff,c1,pm,mu2,e1,e2,e3,e4)
      implicit none
	integer, intent(in) :: diff
      complex(ki), dimension(0:15), intent(in) :: c1
      real(ki), dimension(4), intent(in) ::  e1, e2
      complex(ki), dimension(4), intent(in) :: e3, e4, pm
      complex(ki), intent(in) :: mu2
      complex(ki) :: poly1

      complex(ki) :: pme1,pme2,pme3,pme4

      pme1=sdot(e1,pm)
      pme2=sdot(e2,pm)
      pme3=sdot(pm,e3)
      pme4=sdot(pm,e4)

	if(diff.ge.1) then
	      poly1=+c1(0) 
	elseif(diff.eq.0) then
	      poly1=+c1(0) &
	     &      +c1(1)*pme1 &
	     &      +c1(2)*pme2 &
	     &      +c1(3)*pme3 &
	     &      +c1(4)*pme4 
	else
	      poly1=+c1(0) &
	     &      +c1(1)*pme1 &
	     &      +c1(2)*pme2 &
	     &      +c1(3)*pme3 &
	     &      +c1(4)*pme4 &
		& + c1(5)*pme1**2 &
		& + c1(6)*pme2**2 &
		& + c1(7)*pme3**2 &
		& + c1(8)*pme4**2 &
		& + c1(10)*pme1*pme3 &
		& + c1(11)*pme1*pme4 &
		& + c1(12)*pme2*pme3 &
		& + c1(13)*pme2*pme4 &
		& + c1(14)*mu2 	     &
		& + c1(15)*pme3*pme4
	endif
  end  function poly1

end module mfunctions

