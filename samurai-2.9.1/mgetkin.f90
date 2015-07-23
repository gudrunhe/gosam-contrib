module mgetkin
	use precision
	use mfunctions, only: sdot
implicit none

	complex(ki), dimension(:,:), allocatable, public :: s_mat
	public :: checksmatalloc, getV4,getV4smat,getV3,getV3smat,get2smat

contains

subroutine checksmatalloc(smatallocated)
	implicit none
	logical, intent(out) :: smatallocated
	if (allocated(s_mat)) then
		smatallocated = .true.
	else
		smatallocated = .false.
	endif
end subroutine

subroutine getV4(nleg,cut4,Vi,V)
	implicit none
	integer,                         intent(in ) :: nleg,cut4
	real(ki), dimension(0:nleg-1,4), intent(in ) :: Vi
	real(ki), dimension(6),          intent(out) :: V

	real(ki), dimension(4)  :: Vi1, Vi2, Vi3, Vi21, Vi31, Vi32
	integer 		:: j1,j2,j3,j4

	j4=cut4/1000
	j3=(cut4-j4*1000)/100
	j2=(cut4-j4*1000-j3*100)/10
	j1=cut4-j4*1000-j3*100-j2*10

		Vi1(:)=Vi(j2,:)-Vi(j1,:)
		Vi2(:)=Vi(j3,:)-Vi(j1,:)
		Vi3(:)=Vi(j1,:)-Vi(j4,:)
		Vi21(:)=Vi(j3,:)-Vi(j2,:)
		Vi31(:)=Vi(j4,:)-Vi(j2,:)
		Vi32(:)=Vi(j4,:)-Vi(j3,:)
		
		V(1)=sdot(Vi1,Vi1)
		V(2)=sdot(Vi2,Vi2)
		V(3)=sdot(Vi3,Vi3)
		V(4)=sdot(Vi21,Vi21)
		V(5)=sdot(Vi31,Vi31)
		V(6)=sdot(Vi32,Vi32)
end subroutine

subroutine getV4smat(cut4,m,V)
implicit none
	integer,                  intent(in ) :: cut4
	real(ki), dimension(0:3), intent(in ) :: m
	real(ki), dimension(6),   intent(out) :: V

	real(ki) :: V1, V2, V3, V21, V31, V32
	integer  :: j1,j2,j3,j4

	j4=cut4/1000
	j3=(cut4-j4*1000)/100
	j2=(cut4-j4*1000-j3*100)/10
	j1=cut4-j4*1000-j3*100-j2*10

	V(1) = s_mat(j2+1, j1+1) + m(1) + m(0)
	V(2) = s_mat(j3+1, j1+1) + m(2) + m(0)
	V(3) = s_mat(j4+1, j1+1) + m(3) + m(0)
	V(4) = s_mat(j3+1, j2+1) + m(2) + m(1)
	V(5) = s_mat(j4+1, j2+1) + m(3) + m(1)
	V(6) = s_mat(j4+1, j3+1) + m(3) + m(2)
end subroutine

subroutine getV3(nleg,cut3,Vi,V)
implicit none
	integer,                         intent(in ) :: nleg,cut3
	real(ki), dimension(0:nleg-1,4), intent(in ) :: Vi
	real(ki), dimension(3),          intent(out) :: V

	integer                :: j1,j2,j3
	real(ki), dimension(4) :: Vi1, Vi2, Vi3

	j3=cut3/100
	j2=(cut3-j3*100)/10
	j1=cut3-j3*100-j2*10

	Vi1(:)=Vi(j2,:)-Vi(j1,:)
	Vi2(:)=Vi(j3,:)-Vi(j2,:)
	Vi3(:)=Vi(j1,:)-Vi(j3,:)
	
	V(1)=sdot(Vi1,Vi1)
	V(2)=sdot(Vi2,Vi2)
	V(3)=sdot(Vi3,Vi3)
end subroutine

subroutine getV3smat(cut3,m,V)
implicit none
	integer,                  intent(in ) :: cut3
	real(ki), dimension(0:2), intent(in ) :: m
	real(ki), dimension(3),   intent(out) :: V

	integer                :: j1,j2,j3
	real(ki)               :: V1, V2, V3

	j3=cut3/100
	j2=(cut3-j3*100)/10
	j1=cut3-j3*100-j2*10

	V(1) = s_mat(j2+1, j1+1) + m(1) + m(0)
	V(2) = s_mat(j3+1, j2+1) + m(2) + m(1)
	V(3) = s_mat(j3+1, j1+1) + m(2) + m(0)
end subroutine


subroutine get2smat(cut2,m,K11)
	integer,                  intent(in ) :: cut2
	real(ki), dimension(0:1), intent(in ) :: m
	real(ki),                 intent(out) :: K11

	integer :: j1,j2

	j2=cut2/10
	j1=cut2-j2*10
	
	K11 = s_mat(j2+1, j1+1) + m(0) + m(1)
end subroutine


end module
