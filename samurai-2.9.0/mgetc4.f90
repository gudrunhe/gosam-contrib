module mgetc4
   use precision, only: ki
   use constants
   use options
   use mfunctions
   use mrestore
   implicit none

   public :: getc4

contains

subroutine getc4(numeval,nleg,rank,c4,cut4,q4,qt,p0,Vi,msq)
	use mglobal, only: MP12,mu2g,mu2t,resit,denst,mu2test,dx
	implicit none
! PARAMETERS=======================================================================!
	! external
	integer, 			    intent(in ) :: nleg, rank, cut4
	complex(ki), dimension(0:5), 	    intent(out) :: c4
	complex(ki), dimension(6,4), 	    intent(in ) :: q4
	complex(ki), dimension(4), 	    intent(in ) :: qt
	real(ki),    dimension(4), 	    intent(in ) :: p0
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	complex(ki),    dimension(0:nleg-1),   intent(in ) :: msq
	! internal
	complex(ki) 		  :: dens5t
	integer 		  :: i,m,i1,i2,i3,i4,i5,j1,j2,j3,j4
	integer 		  :: dicut5,diff, acc
	complex(ki), dimension(6) :: dens5, dens4, xneval, resi5, f4
	complex(ki), dimension(6) :: mu2vec
	logical evalres
! INTERFACE========================================================================!
	interface
	    function numeval(ncut, Q, mu2)
	            use precision
	            implicit none
	            integer, 		       intent(in) :: ncut
	            complex(ki), dimension(4), intent(in) :: Q
	            complex(ki), 	       intent(in) :: mu2
	            complex(ki) 			  :: numeval
	    end function numeval
	end interface
! INITIALIZATION===================================================================!
	dens4(:)  = cone
	resi5(:)  = czip
	xneval(:) = czip
	! mu2
		mu2test(4) = mu2t(4)
		mu2vec = (/ czip, czip, mu2g(4), mu2g(4),-mu2g(4),-mu2g(4) /)
	! for lnntest
		resit(4) = czip
		denst(4) = cone
	! for simplified sampling
	diff = nleg-rank
	! cuts
	j4  =  cut4/1000
	acc =  j4*1000
	j3  = (cut4-acc)/100
	acc =  acc + j3*100
	j2  = (cut4-acc)/10
	j1  =  cut4-acc-j2*10
! CALCULATING THE RESIDUES=========================================================!
      	nleg_ne_4: if (nleg.ne.4) then
          nleg_eq_5: if (nleg.eq.5) then
            do m=1,6
               resi5(m)=Res5(1,mu2vec(m))
            enddo
            resit(4)=Res5(1,mu2t(4))
	  else
            dicut5=1
            loop_i5: do i5=4,nleg-1
               loop_i4: do i4=3,i5-1
                  loop_i3: do i3=2,i4-1
                     loop_i2: do i2=1,i3-1
                        loop_i1: do i1=0,i2-1
                           dens5(:)=cone
                           dens5t=cone
                           evalres=.false.
                           loop110: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)&
                                   &  .and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i.eq.j1).or.(i.eq.j2)& 
                                       &  .or.(i.eq.j3).or.(i.eq.j4)) then
                                    dens5(:)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop110
                                 else
                                    do m=1,6
                                       dens5(m)=dens5(m)*denevalmu2(&
                                             &nleg,i,q4(m,:),Vi,msq,mu2vec(m))
                                    enddo
                                    dens5t=dens5t*denevalmu2(&
                                       &nleg,i,qt,Vi,msq,mu2t(4))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop110
                           if (evalres) then
                              do m=1,6
                                 resi5(m)=resi5(m)+dens5(m)*&
                                     &  Res5(dicut5,mu2vec(m))
                              enddo
                              resit=resit+dens5t*Res5(dicut5,mu2t(4))
                           endif
                           dicut5=dicut5+1
                        enddo loop_i1
                     enddo loop_i2
                  enddo loop_i3
               enddo loop_i4
            enddo loop_i5
	  end if nleg_eq_5
          ! 111 continue
          do i=0,nleg-1
            if ((i.ne.j1).and.(i.ne.j2) &
            &  .and.(i.ne.j3).and.(i.ne.j4)) then
               do m=1,6
                  dens4(m)=dens4(m)*denevalmu2(nleg,i,q4(m,:),Vi,msq,mu2vec(m))
               enddo
               denst=denst*denevalmu2(nleg,i,qt,Vi,msq,mu2t(4))
            endif
          enddo
	end if nleg_ne_4
      	! 116 continue
! CALCULATING THE POLYNOMIALS======================================================!
      	do m=1,6
         xneval(m)=numeval(cut4,q4(m,:),mu2vec(m))
      	enddo
	if (imeth.eq.'diag') then
		f4(:) = (xneval(:)-resi5(:))/dens4(:)
	elseif (imeth.eq.'tree') then
		f4(:) = xneval(:)-resi5(:)/dens4(:)
	endif
! CALCULATING THE C4===============================================================!
	select case(diff)
	case(-1)
      c4(0) = & 
     &half*(f4(1) + f4(2))
      c4(1) = & 
     &(half*(f4(1) - one*f4(2)))/(dx(1)*MP12(4))
      c4(2) = & 
     &(one*(f4(3) + f4(4) - one*f4(5) - one*f4(6)))/(four*mu2g(4))
      c4(3) = & 
     &        (one*(dx(5)*(f4(3) - one*f4(4)) + &
     &      dx(3)*(-(one*f4(5)) + f4(6))))/&
     &  (four*dx(3)*dx(5)*MP12(4)*mu2g(4))
      c4(4) = & 
     &        (one*(-(two*f4(1)) - two*f4(2) + f4(3) + f4(4) + f4(5) + &
     &      f4(6)))/(four*mu2g(4)**2)
      c4(5) = & 
     &        (one*(dx(1)*dx(5)*(f4(3) - one*f4(4)) + &
     &      dx(3)*(two*dx(5)*(-(one*f4(1)) + f4(2)) + &
     &         dx(1)*(f4(5) - one*f4(6)))))/&
     &  (four*dx(1)*dx(3)*dx(5)*MP12(4)*mu2g(4)**2)
! RANK 0 ==========================================================================!
	case default
      c4(0)= (f4(1)+f4(2))/two
      c4(1)= (f4(1)-f4(2))/(two*dx(1)*MP12(4))
      c4(3)= -(two*dx(3)*MP12(4)*c4(1)-f4(3)+f4(4))/(two*dx(3)*MP12(4)*mu2g(4))
      c4(2)= ((-dx(3)+dx(5))*MP12(4)*c4(1)-(dx(3)+dx(5))*MP12(4)*mu2g(4)*c4(3) &
        & +f4(3)-f4(5))/(two*mu2g(4))
      c4(4)= -((c4(0)+mu2g(4)*c4(2)+dx(3)*MP12(4)*(c4(1)+mu2g(4)*c4(3)) &
        & -f4(3))/mu2g(4)**2)
	end select
! LOWER RANK CS====================================================================!
	if (diff.ge.0) then
		c4(5)=czip
	if (diff.ge.1) then
		c4(4)=czip
	if (diff.ge.2) then
		c4(3)=czip
	if (diff.ge.3) then
		c4(2)=czip
	endif 
	endif
	endif
	endif
end subroutine getc4

end module mgetc4




