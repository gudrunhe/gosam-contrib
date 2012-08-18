module options
   use precision, only: ki
   implicit none

   private :: ki

   integer :: isca, verbosity, itest, iresc
   character(len=4) :: imeth
   logical :: meth_is_tree
   logical :: meth_is_diag

   integer :: iout = 20
   integer :: ibad = 30

   logical :: use_maccu = .false.

   ! value of C0 at which we switch between the two different samplings
   real(ki), parameter :: C0_thrs = 1.0E-10_ki

end module options

