module mgetc5
   use precision, only: ki
   use constants
   use options
   use mfunctions
   implicit none

   public :: getc5

contains

subroutine getc5(numeval,nleg,c5,cut5,Vi,msq,q5,mu2)
      implicit none

      interface
         function     numeval(ncut, Q, mu2)
           use precision
           implicit none
           integer, intent(in) :: ncut
           complex(ki), dimension(4), intent(in) :: Q
           complex(ki), intent(in) :: mu2
           complex(ki) :: numeval
         end function numeval
      end interface
! PARAMETERS========================================================!
	! external
	integer, 			     intent(in   ) :: nleg
	complex(ki), 			     intent(  out) :: c5
	integer, 			     intent(in   ) :: cut5
	real(ki),    dimension(0:nleg-1, 4), intent(in   ) :: Vi
	complex(ki),    dimension(0:nleg-1),    intent(in   ) :: msq
	complex(ki), dimension(4), 	     intent(inout) :: q5
	complex(ki), 			     intent(in   ) :: mu2
	! internal
	integer     :: i,j1,j2,j3,j4,j5,acc
	complex(ki) :: denoms
! INITIALIZATION====================================================!
	! cuts
	if(meth_is_diag) then
		j5  =  cut5/10000
		acc =  j5*10000
		j4  = (cut5-acc)/1000
		acc =  acc + j4*1000
		j3  = (cut5-acc)/100
		acc =  acc + j3*100
		j2  = (cut5-acc)/10
		j1  =  cut5-acc-j2*10
! CALCULATING C5====================================================!
	denoms=cone
	do i=0,nleg-1
		if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)&
		& .and.(i.ne.j4).and.(i.ne.j5)) then
		    denoms = denoms*denevalmu2(nleg,i,q5,Vi,msq,mu2)
		endif
	enddo
		c5 = numeval(cut5,q5,mu2)/denoms/mu2
	elseif (meth_is_tree) then
		c5=numeval(cut5,q5,mu2)/mu2
	else
		print*, "imeth must be 'tree' or 'diag'"
		stop
	endif
end subroutine getc5

end module mgetc5




