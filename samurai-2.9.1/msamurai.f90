module msamurai
   use precision
   use constants
   use options
   use ncuts
   use notfirst
   use ltest
   use madds
   use mgetqs
   use mgetbase
   use mrestore
!   use mgetc5
!   use mgetc4
!   use mgetc3
!   use mgetc2
!   use mgetc1
   use mtests
   use mtens,  only: tensor_reconstruction, numetens
	use mglobal, only: mu2g
   use mmasters
   use mgetkin
   implicit none

   private

   public :: initsamurai, exitsamurai, InitDenominators, samurai

contains

subroutine initsamurai(ameth, asca, averb, atest, aresc)
	implicit none
	integer,          intent(in)           :: asca, averb, atest
	character(len=4), intent(in)           :: ameth
	integer,          intent(in), optional :: aresc
	
	call banner
	
	
	if (ameth.eq.'diag') then
		imeth = ameth
		meth_is_tree = .false.
		meth_is_diag = .true.
	elseif (ameth.eq.'tree') then
		imeth = ameth
		meth_is_tree = .true.
		meth_is_diag = .false.
	else
		write(6,*) 'incompatible value for imeth'
		write(6,*) 'imeth =','ameth'
		stop
	endif
	
	if ((asca.ge.1).and.(asca.le.4)) then
		isca = asca
	else
		write(6,*) 'incompatible value for isca'
		write(6,*) 'isca =',isca
		stop
	endif
	
	if ((averb.ge.0).or.(averb.le.3)) then
		verbosity = averb
	else
		write(6,*) 'incompatible value for verbosity'
		write(6,*) 'verbosity =',averb
		stop
	endif
	
	if ((atest.ge.0).or.(atest.le.3)) then
		itest = atest
	else
		write(6,*) 'incompatible value for itest'
		write(6,*) 'itest =',atest
		stop
	endif
	
	if (present(aresc)) then
		if ((0.le.aresc).and.(aresc.le.3)) then
			iresc = aresc
		else
			write(6,*) 'incompatible value for iresc'
			write(6,*) 'iresc =',aresc
			stop
		endif
	else
		iresc=0
	endif
	
	if (imeth.eq.'tree') itest=2
	
	call initmasters()
	
	901  format(a40,a4,a10,I1,a15,I1,a11,I1)
	904  format(a68,I1)
	
	if ( (verbosity.gt.0) .and. (notfirstp.eqv.(.false.)) ) then
		open(unit=iout,file='output.dat',status='unknown')
		write(iout,*) '---------------------------------------------------------------------------------'
		write(iout,901) ' SAMURAI called with arguments: imeth = ',imeth, &
		&               ' ; isca = ',isca,' ; verbosity = ',verbosity,' ; itest = ',itest
		write(iout,*) '---------------------------------------------------------------------------------'
		notfirstp = .TRUE.
	endif
	if ( (itest.gt.0) .and. (notfirstd.eqv.(.false.)) .and. (ibad.gt.0) .and. (verbosity.gt.0)  ) then
		open(unit=ibad,file='bad.points',status='unknown')
		write(ibad,*) '--------------------------------------------------------------------'
		write(ibad,904) ' Points that have been discarded by SAMURAI because failing itest = ',itest
		write(ibad,*) '--------------------------------------------------------------------'
		notfirstd = .true.
	endif
	
	call rtlimit

end subroutine initsamurai



subroutine rtlimit
	implicit none
	real    :: tpwlimit,tnnlimit,tnnlimit4,tnnlimit3,tnnlimit2,tnnlimit1
	integer :: ierr
	open(unit=10,file='ltest.dat',status='old',iostat=ierr)
	
	if (ierr.eq.0) then
		read(10,*)
		read(10,*)
		read(10,*)  tpwlimit
		pwlimit   = real(tpwlimit, ki)
		read(10,*)  tnnlimit
		nnlimit   = real(tnnlimit, ki)
		read(10,*)  tnnlimit4
		lnnlimit4 = real(tnnlimit4, ki)
		read(10,*)  tnnlimit3
		lnnlimit3 = real(tnnlimit3, ki)
		read(10,*)  tnnlimit2
		lnnlimit2 = real(tnnlimit2, ki)
		read(10,*)  tnnlimit1
		lnnlimit1 = real(tnnlimit1, ki)
		
		close(10)
	else
		pwlimit   = 1.0E-03_ki
		nnlimit   = 1.0E-03_ki
		lnnlimit4 = 1.0E-02_ki
		lnnlimit3 = 1.0E-02_ki
		lnnlimit2 = 1.0E-02_ki
		lnnlimit1 = 1.0E+01_ki
	end if
end subroutine rtlimit



subroutine exitsamurai()
	implicit none
	
	if (verbosity.gt.0) close(iout)
	if (itest.gt.0 .and. ibad.gt.0 .and. verbosity.gt.0) close(ibad)
	if (isca.eq.4) call exitmasters
end subroutine exitsamurai



subroutine banner
        implicit none
        
        print *
        print *, ' ********************************************************************'
        print *, ' ********************** SAMURAI - version 2.9.1'
        print *, ' ********************************************************************'
        print *, ' *                                                                  *'
        print *, ' *                                                                  *'
        print *, ' * Authors: H. van Deurzen, P. Mastrolia,                           *'
        print *, ' *          E. Mirabella, G. Ossola, and F. Tramontano              *'
        print *, ' *                                                                  *'
        print *, ' * hdeurzen@mpp.mpg.de                                              *'
        print *, ' * pierpaolo.mastrolia@cern.ch                                      *'
        print *, ' * mirabell@mpp.mpg.de                                              *'
        print *, ' * gossola@citytech.cuny.edu                                        *'
        print *, ' * francesco.tramontano@cern.ch                                     *'
        print *, ' *                                                                  *'
        print *, ' *  Based on:                                                       *'
        print *, ' *    - Mastrolia, Ossola, Reiter, Tramontano,                      *'
        print *, ' *      JHEP 1008 (2010) 080, arXiv:1006.0710                       *'
        print *, ' *    - van Deurzen, Acta Phys.Polon. B44 (2013) 11, 2223-2230      *'
        print *, ' *                                                                  *'
        print *, ' *  On the web:  http://cern.ch/samurai                             *'
        print *, ' *                                                                  *'
        print *, ' ********************************************************************'
        print *, ' *                                                                  *'
        print *, ' * output files: <output.log>   [ for verbosity.gt.0 ]              *'
        print *, ' *                                                                  *'
        print *, ' *               <bad.points>   [ for itest.gt.0     ]              *'
        print *, ' *                                                                  *'
        print *, ' ********************************************************************'

end subroutine banner



pure subroutine InitDenominators(nleg,Vi,msq,Q0,m0,Q1,m1,Q2,m2,Q3,m3,Q4,m4,Q5,m5,Q6,m6,Q7,m7)
      implicit none
      integer,                            intent(in ) :: nleg
      real(ki),    dimension(0:nleg-1,4), intent(out) :: Vi
      complex(ki), dimension(0:nleg-1),   intent(out) :: msq

      real(ki), dimension(4), intent(in), optional :: Q0
      complex(ki),            intent(in), optional :: m0
      real(ki), dimension(4), intent(in), optional :: Q1
      complex(ki),            intent(in), optional :: m1
      real(ki), dimension(4), intent(in), optional :: Q2
      complex(ki),            intent(in), optional :: m2
      real(ki), dimension(4), intent(in), optional :: Q3
      complex(ki),            intent(in), optional :: m3
      real(ki), dimension(4), intent(in), optional :: Q4
      complex(ki),            intent(in), optional :: m4
      real(ki), dimension(4), intent(in), optional :: Q5
      complex(ki),            intent(in), optional :: m5
      real(ki), dimension(4), intent(in), optional :: Q6
      complex(ki),            intent(in), optional :: m6
      real(ki), dimension(4), intent(in), optional :: Q7
      complex(ki),            intent(in), optional :: m7

      if(present(Q0) .and. present(m0)) then
         Vi(0,:) = Q0
         msq(0)  = m0*m0
      end if
      if(present(Q1) .and. present(m1)) then
         Vi(1,:) = Q1
         msq(1)  = m1*m1
      end if
      if(present(Q2) .and. present(m2)) then
         Vi(2,:) = Q2
         msq(2)  = m2*m2
      end if
      if(present(Q3) .and. present(m3)) then
         Vi(3,:) = Q3
         msq(3)  = m3*m3
      end if
      if(present(Q4) .and. present(m4)) then
         Vi(4,:) = Q4
         msq(4)  = m4*m4
      end if
      if(present(Q5) .and. present(m5)) then
         Vi(5,:) = Q5
         msq(5)  = m5*m5
      end if
      if(present(Q6) .and. present(m6)) then
         Vi(6,:) = Q6
         msq(6)  = m6*m6
      end if
      if(present(Q7) .and. present(m7)) then
         Vi(7,:) = Q7
         msq(7)  = m7*m7
      end if

end  subroutine InitDenominators















subroutine samurai(numeval,tot,totr,Vi,msq,nleg,rank,istop,scale2,ok,cache_flag, scalar_cache,rank_numeval)
	use options, only: use_maccu
	use maccu
	implicit none

	integer,                            intent(in ) :: nleg, rank, istop
	real(ki),                           intent(in ) :: scale2
	complex(ki),                        intent(out) :: totr
	complex(ki), dimension(0:nleg-1),   intent(in ) :: msq
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	complex(ki), dimension(-2:0),       intent(out) :: tot
	logical,                            intent(out) :: ok
	logical,     intent(inout), optional                    :: cache_flag
	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
	optional :: rank_numeval
	
	type(accumulator_type), dimension(-2:0) :: acc_re, acc_im
	type(accumulator_type) :: accr_re, accr_im
	
	integer :: i,j,k
	integer :: ep
	integer :: diff 
	integer :: cache_offset
	
	complex(ki) :: mu2, mu2test
	
	complex(ki) :: totr2sum,totr3sum,totr4sum
	complex(ki), dimension(4)    :: qt, qtest
	complex(ki), dimension(-2:0) :: tot4sum, tot3sum, tot2sum, tot1sum
	
	
	interface
	   function     numeval(ncut, Q, mu2)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     complex(ki), dimension(4), intent(in) :: Q
	     complex(ki),               intent(in) :: mu2
	     complex(ki) :: numeval
	   end function numeval
	end interface
	
	interface
	   function     rank_numeval(ncut)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     integer :: rank_numeval
	   end function rank_numeval
	end interface

	if (nleg.gt.maxleg) then
		write(6,*) 'reduction called for nleg.gt.maxleg'
		write(6,*) 'maxleg, nleg ',maxleg,nleg
		write(6,*) 'please change the values max1,..,max5,maxleg in constants.f90'
		write(6,*) 'and compile again the library'
		stop
	endif

	if(rank.lt.0 .and. .not. present(rank_numeval)) then
		write (6, *) 'dynamical rank (rank<0) needs rank_numeval parameter'
		stop
	endif
	
	cache_offset = 0
	diff = nleg-rank
	
	if(verbosity.gt.0)then
		write(iout,*)
		write(iout,*)
		write(iout,*) 'Denominators: '
		do k=0,nleg-1
			write(iout,902) ' Pi(',k,') = ', Vi(k,:)
			write(iout,903) 'msq(',k,') = ', msq(k)
			write(iout,*)
		enddo
		write(iout,*)' '
	endif
	
	totr       = czip
	totr4sum   = czip
	totr3sum   = czip
	totr2sum   = czip
	
	tot(:)     = czip
	tot4sum(:) = czip
	tot3sum(:) = czip
	tot2sum(:) = czip
	tot1sum(:) = czip
	
	
	if (use_maccu) then
		do ep=-2,0
			acc_re(ep)%a(:) = 0.0_ki
			acc_im(ep)%a(:) = 0.0_ki
		end do
		accr_re%a(:) = 0.0_ki
		accr_im%a(:) = 0.0_ki
	end if
	
	notfirsti = .false.
	
	if (iresc.eq.1) call tensor_reconstruction(numeval,nleg,rank) 
	
	!=========5ple cut
	if(nleg.ge.5 .and. istop.le.5) then
		if(verbosity.gt.0) write(iout,*) 'Pentagon coefficients: '
		call samurai5plecut(numeval,nleg, rank, Vi, msq, scale2)
	endif
	!=========4ple cut
	if(nleg.ge.4 .and. istop.le.4) then
		if(verbosity.gt.0) write(iout,*) 'Box coefficients: '
		call samurai4plecut(numeval,nleg,rank,Vi,msq,scale2,tot4sum,totr4sum,ok,cache_flag,scalar_cache,cache_offset,rank_numeval)
	endif
	!-========3ple cut
	if(nleg.ge.3 .and. istop.le.3) then
		if(verbosity.gt.0) write(iout,*) 'Triangle coefficients: '
		call samurai3plecut(numeval,nleg,rank,Vi,msq,scale2,tot3sum,totr3sum,ok,cache_flag,scalar_cache,cache_offset,rank_numeval)
	endif
	!=========2ple cut
	if(nleg.ge.2 .and. istop.le.2) then
		if(verbosity.gt.0) write(iout,*) 'Bubble coefficients: '
		call samurai2plecut(numeval,nleg,rank,Vi,msq,scale2,tot2sum,totr2sum,ok,cache_flag,scalar_cache,cache_offset,rank_numeval)
	endif
	!========1ple cut
	if(nleg.ge.1 .and. istop.le.1) then
		if(verbosity.gt.0) write(iout,*) 'Tadpole coefficients: '
		call samurai1plecut(numeval,nleg,rank,Vi,msq,scale2,tot1sum,ok,cache_flag,scalar_cache,cache_offset,rank_numeval)
	endif
	!===============
	
	tot(:) = tot4sum(:) + tot3sum(:) + tot2sum(:) + tot1sum(:)
	totr   = totr4sum   + totr3sum   + totr2sum  
	
	if (use_maccu) then
		tot = reduce_accu(acc_re, acc_im)
		totr = reduce_accu(accr_re, accr_im)
	end if
	
	
	! TEST N=N -------------------------------------------
	if (itest.eq.1) then
		qtest=10.3_ki*cone
		mu2test=13.0_ki
		call nntest(diff,numeval,qtest,mu2test,nleg,Vi,msq,ok)
	endif
	!-----------------------------------------------------         
	
	! POWER test ----------------------------------
	if (itest.eq.3) call pwtest(nleg,rank,ok)
	! ---------------------------------------------   
	
	
	if (ok) then
		if(verbosity.gt.0)then
			write(iout,*)
			write(iout,*)' Result: '
			write(iout,*)' Double   Pole = ', tot(-2)
			write(iout,*)' Single   Pole = ', tot(-1)
			write(iout,*)' Finite   Part = ', tot(0)
			write(iout,*)
			write(iout,*)'[Rational Part = ', totr,']'
			write(iout,*)
			write(iout,*)
		endif
		
	else    
		
		if(ibad.gt.0.and.verbosity.gt.0)then
			write(ibad,*) 'Denominators: '
			do k=0,nleg-1
				write(ibad,902) ' Pi(',k,') = ', Vi(k,:)
				write(ibad,903) 'msq(',k,') = ', msq(k)
				write(ibad,*)
			enddo
			write(ibad,*)'---------------------------------------- '
			write(ibad,*)' '
		endif
	endif
	
	if (present(cache_flag)) cache_flag = .true.
	
	902  format(a4,I1,a4,4(D24.15))
	903  format(a4,I1,a4,1(D24.15))
	
	9005 format(A3,I5,A6,D24.15,A1,D24.15,A3)
	9004 format(A3,I4,A1,I1,A5,D24.15,A1,D24.15,A3)
	9003 format(A3,I3,A1,I2,A6,D24.15,A1,D24.15,A3)
	9002 format(A3,I2,A1,I2,A7,D24.15,A1,D24.15,A3)
	9001 format(A3,I2,A1,I2,A7,D24.15,A1,D24.15,A3)
	
end subroutine samurai
!---#] subroutine samurai:


subroutine samurai5plecut(numeval,nleg, rank, Vi, msq, scale2)
	use mgetc5
	implicit none
	
	integer,                            intent(in) :: nleg, rank
	real(ki),                           intent(in) :: scale2
	complex(ki), dimension(0:nleg-1),   intent(in) :: msq
	real(ki),    dimension(0:nleg-1,4), intent(in) :: Vi
	
	integer     :: j1,j2,j3,j4,j5,icut5
	integer     :: cut5,n1
	integer     :: diff 
	complex(ki) :: mu2, mu2test
	real(ki)    :: r1, r2
	complex(ki) :: c5
	real(ki),    dimension(4) :: k1, k2, k3, e1, e2, p0, L3
	complex(ki), dimension(4) :: q5, e3, e4


	interface
		function     numeval(ncut, Q, mu2)
			use precision
			implicit none
			integer,                   intent(in) :: ncut
			complex(ki), dimension(4), intent(in) :: Q
			complex(ki),               intent(in) :: mu2
			complex(ki) :: numeval
		end function numeval
	end interface

	9005 format(A3,I5,A6,D24.15,A1,D24.15,A3)

	icut5=1
	do j5=4,nleg-1
	do j4=3,j5-1
	do j3=2,j4-1
	do j2=1,j3-1
	do j1=0,j2-1
		cut5=j5*10000+j4*1000+j3*100+j2*10+j1
		do n1=1,4
			k1(n1)=Vi(j2,n1)-Vi(j1,n1)
			k2(n1)=Vi(j1,n1)-Vi(j5,n1)
			p0(n1)=Vi(j1,n1)
		enddo
		
		call getbase(k1,k2,r1,r2,e1,e2,e3,e4)
		call getq5(nleg,cut5,e1,e2,e3,e4,p0,Vi,msq,r1,r2,q5,mu2)
		
		if (iresc.eq.1) then
			call getc5(numetens,nleg,c5,cut5,Vi,msq,q5,mu2)
		else
			call getc5(numeval,nleg,c5,cut5,Vi,msq,q5,mu2)
		endif
		
		call store5(icut5,cut5,p0,e1,e2,e3,e4,c5)
		
		if(verbosity.gt.0)then
			write(iout,9005) 'c5(', cut5,')  = (',real(c5),',',aimag(c5),'  )'
			write(iout,*)
		endif
		
		icut5=icut5+1
	enddo
	enddo
	enddo
	enddo
	enddo
	nc5=icut5-1
end subroutine





subroutine samurai4plecut(numeval,nleg,rnk,Vi,msq,scale2,tot4sum,totr4sum,ok,cache_flag,scalar_cache,cache_offset,rank_numeval)
	use options, only: use_maccu
	use maccu
	use mgetc4
	implicit none
	
	integer,                            intent(in ) :: nleg, rnk
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	complex(ki), dimension(0:nleg-1),   intent(in ) :: msq
	real(ki),                           intent(in ) :: scale2
	complex(ki),                        intent(out) :: totr4sum
	complex(ki), dimension(-2:0),       intent(out) :: tot4sum
	logical,                            intent(out) :: ok
	logical,     intent(inout), optional                    :: cache_flag
	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
	integer,     intent(inout), optional                    :: cache_offset
	optional :: rank_numeval
	
	type(accumulator_type), dimension(-2:0) :: acc_re,  acc_im
	type(accumulator_type)                  :: accr_re, accr_im
	
	integer :: i,j1,j2,j3,j4,icut4
	integer :: cut4,n1,ep
	integer :: diff 
	integer :: rank
	
	complex(ki) :: mu2, mu2test
	real(ki)    :: r1, r2
	real(ki), dimension(4) :: k1, k2, k3, e1, e2, p0, L3
	
	logical     :: smatallocated
	complex(ki) :: totr4tmp
	complex(ki), dimension(6,4)  :: q4
	complex(ki), dimension(4)    :: e3, e4, qt, qtest
	complex(ki), dimension(0:5)  :: c4
	complex(ki), dimension(-2:0) :: tot4tmp
	
	complex(ki), dimension(0:3)  :: m
	real(ki),    dimension(6)    :: V
	complex(ki), dimension(-2:0) :: MI4
	
	
	interface
	   function     numeval(ncut, Q, mu2)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     complex(ki), dimension(4), intent(in) :: Q
	     complex(ki),               intent(in) :: mu2
	     complex(ki) :: numeval
	   end function numeval
	end interface

	interface
	   function     rank_numeval(ncut)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     integer :: rank_numeval
	   end function rank_numeval
	end interface

	9041 format(A3,I4,A1,I1,A5,D24.15,A1,D24.15,A3)
	9042 format(A3,I4,A1,I2,A5,D24.15,A1,D24.15,A3)


	diff=nleg-rank
	call checksmatalloc(smatallocated)
	tot4sum  = czip
	totr4sum = czip
	icut4=1
	do j4=3,nleg-1
	do j3=2,j4-1
	do j2=1,j3-1
	do j1=0,j2-1
		cut4=j4*1000+j3*100+j2*10+j1
		
		m(0)=msq(j1)
		m(1)=msq(j2)
		m(2)=msq(j3)
		m(3)=msq(j4)
		
		do n1=1,4
			k1(n1)=Vi(j2,n1)-Vi(j1,n1)
			k2(n1)=Vi(j1,n1)-Vi(j4,n1)
			k3(n1)=Vi(j3,n1)-Vi(j1,n1)
			p0(n1)=Vi(j1,n1)
			L3(n1)=Vi(j4,n1)-Vi(j3,n1)
		enddo
		
		call getbase(k1,k2,r1,r2,e1,e2,e3,e4)
		call getq4(nleg,cut4,e1,e2,e3,e4,p0,k1,k2,k3,L3,r1,r2,q4,qt,msq)

                if (rnk.lt.0) then
                  rank=rank_numeval(cut4)
                  if (rank.lt.0) cycle
                else
                  rank=rnk
                endif

		if (iresc.eq.1) then
			call getc4(numetens,nleg,rank,c4,cut4,q4,qt,p0,Vi,msq)
		else
			call getc4(numeval,nleg,rank,c4,cut4,q4,qt,p0,Vi,msq)
		endif
		
		if (itest.eq.2) call lnntest4(diff,numeval,cut4,c4,qt,p0,L3,e3,e4,ok)
		call store4(icut4,cut4,L3,p0,e1,e2,e3,e4,c4)
		
		if (smatallocated) then
			call getV4smat(cut4,real(m),V)
		else
			call getV4(nleg,cut4,Vi,V)
		endif
		
		if (present(cache_flag)) then
			call getMI4(V,m,scale2,MI4,cache_flag, cache_offset, scalar_cache)
		else
			call getMI4(V,m,scale2,MI4)
		end if
		
		call add4(c4,MI4,tot4tmp,totr4tmp)
		
		if (verbosity.ge.2) then 
			do ep=-2,0
			write(iout,9042) 'I4(',cut4,',',ep,') = (',real(MI4(ep)),',',aimag(MI4(ep)),'  )'
			enddo
		endif
		
		if(verbosity.gt.0)then
			do i=0,5
			write(iout,9041) 'c4(', cut4,',',i,') = (',real(c4(i)),',',aimag(c4(i)),'  )'
			enddo
		write(iout,*)
		endif
		
		if (use_maccu) then
			call add_accu(accr_re, accr_im, totr4sum)
			call add_accu(acc_re(:), acc_im(:), tot4sum(:))
		else 
			tot4sum(:)=tot4sum(:)+tot4tmp(:)
			totr4sum=totr4sum+totr4tmp
		end if
		
		icut4=icut4+1
	enddo
	enddo
	enddo
	enddo
	nc4=icut4-1
end subroutine samurai4plecut

subroutine samurai3plecut(numeval,nleg,rnk,Vi,msq,scale2,tot3sum,totr3sum,ok,cache_flag,scalar_cache, cache_offset,rank_numeval)
	use options, only: use_maccu
	use maccu
	use mgetc3
	implicit none
	integer,                            intent(in ) :: nleg, rnk
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	complex(ki), dimension(0:nleg-1),   intent(in ) :: msq
	real(ki),                           intent(in ) :: scale2
	complex(ki), dimension(-2:0),       intent(out) :: tot3sum
	complex(ki),                        intent(out) :: totr3sum
	logical,                            intent(out) :: ok
	logical,     intent(inout), optional                    :: cache_flag
	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
	integer,     intent(inout), optional                    :: cache_offset
	optional :: rank_numeval
	
	type(accumulator_type), dimension(-2:0) :: acc_re, acc_im
	type(accumulator_type) :: accr_re, accr_im
	
	integer :: rank
	integer :: i,j,ep,j1,j2,j3,icut3
	integer :: cut3,n1
	integer :: diff 
	
	logical     :: smatallocated
	complex(ki) :: mu2, mu2test
	real(ki)    :: r1, r2
	real(ki), dimension(4):: k1, k2, k3, e1, e2, p0, L3
	
	complex(ki) :: totr3tmp
	complex(ki), dimension(15,4) :: q3
	complex(ki), dimension(4)    :: e3, e4, qt, qtest
	complex(ki), dimension(0:14) :: c3
	complex(ki), dimension(-2:0) :: tot3tmp
	
	complex(ki), dimension(0:2)  :: m
	real(ki),    dimension(3)    :: V
	complex(ki), dimension(-2:0) :: MI3
	
	interface
	   function     numeval(ncut, Q, mu2)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     complex(ki), dimension(4), intent(in) :: Q
	     complex(ki),               intent(in) :: mu2
	     complex(ki) :: numeval
	   end function numeval
	end interface

	interface
	   function     rank_numeval(ncut)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     integer :: rank_numeval
	   end function rank_numeval
	end interface

	9031 format(A3,I3,A1,I2,A5,D24.15,A1,D24.15,A3)
	9032 format(A3,I3,A1,I2,A6,D24.15,A1,D24.15,A3)

	call checksmatalloc(smatallocated)

        if (rnk.lt.0) goto 5

	diff=nleg-rnk
	if (diff.ge.4 .and. itest.ne.1) then
		do j=0,14
		c3(j)=czip
		enddo
                goto 10
        endif

 5      continue


		tot3sum  = czip
		totr3sum = czip
		icut3=1
		do j3=2,nleg-1
		do j2=1,j3-1
		do j1=0,j2-1
		
			cut3=j3*100+j2*10+j1
			m(0)=msq(j1)
			m(1)=msq(j2)
			m(2)=msq(j3)
			
			do n1=1,4
				k1(n1)=Vi(j2,n1)-Vi(j1,n1)
				k2(n1)=Vi(j1,n1)-Vi(j3,n1)
				p0(n1)=Vi(j1,n1)
			enddo

                if (rnk.lt.0) then
                  rank=rank_numeval(cut3)
                  if (rank.lt.0) cycle
                else
                  rank=rnk
                endif

			call getbase(k1,k2,r1,r2,e1,e2,e3,e4)
			call getq3(nleg,rank,cut3,e1,e2,e3,e4,p0,k1,k2,msq,r1,r2,q3,qt)
			
			if (iresc.eq.1) then
				call getc3(numetens,nleg,rank,c3,cut3,q3,qt,Vi,msq)
			else 
				call getc3(numeval,nleg,rank,c3,cut3,q3,qt,Vi,msq)
			endif
			
			if (itest.eq.2) call lnntest3(diff,numeval,cut3,c3,qt,p0,e3,e4,ok)
			call store3(icut3,cut3,p0,e1,e2,e3,e4,c3)
			
			if (smatallocated) then
				call getV3smat(cut3,real(m),V)
			else
				call getV3(nleg,cut3,Vi,V)
			endif
			
			if (present(cache_flag)) then
				call getMI3(V,m,scale2,MI3,cache_flag, cache_offset, scalar_cache)
			else
				call getMI3(V,m,scale2,MI3)
			end if
			
			call add3(V,m,c3,MI3,tot3tmp,totr3tmp)
			
			if (verbosity.ge.2) then
				do ep=-2,0
				write(iout,9031) 'I3(',cut3,',',ep,') = (',real(MI3(ep)),',',aimag(MI3(ep)),'  )'
				enddo
			endif
			if(verbosity.gt.0.and.diff.le.3)then
				do i=0,14
				write(iout,9032) 'c3(',cut3,',',i,')  = (',real(c3(i)),',',aimag(c3(i)),'  )'
				enddo
			write(iout,*)
			endif
			
			if (use_maccu) then
				call add_accu(accr_re, accr_im, totr3sum)
				call add_accu(acc_re(:), acc_im(:), tot3sum(:))
			else
				tot3sum(:) = tot3sum(:) + tot3tmp(:)
				totr3sum   = totr3sum   + totr3tmp
			end if
			
			icut3=icut3+1
		enddo
		enddo
		enddo
		nc3=icut3-1

 10	continue
end subroutine

subroutine samurai2plecut(numeval,nleg,rnk,Vi,msq,scale2,tot2sum,totr2sum,ok,cache_flag,scalar_cache, cache_offset,rank_numeval)
	use options, only: use_maccu
	use maccu
	use mgetc2
	implicit none
	
	integer,                            intent(in ) :: nleg, rnk
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	complex(ki), dimension(0:nleg-1),   intent(in ) :: msq
	real(ki),                           intent(in ) :: scale2
	complex(ki),                        intent(out) :: totr2sum
	complex(ki), dimension(-2:0),       intent(out) :: tot2sum
	logical,                            intent(out) :: ok
	logical,     intent(inout), optional                    :: cache_flag
	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
	integer,     intent(inout), optional                    :: cache_offset
	optional :: rank_numeval
	
	type(accumulator_type), dimension(-2:0) :: acc_re, acc_im
	type(accumulator_type) :: accr_re, accr_im

        integer :: rank	
	integer :: i,j,k,j1,j2,icut2
	integer :: cut2,n1,ep
	integer :: diff 
	
	complex(ki) :: mu2, mu2test
	real(ki)    :: r1, r2, factor
	real(ki), dimension(4) :: k1, k2, k3, e1, e2, p0, L3
	
	logical     :: highrank, smatallocated
	complex(ki) :: totr2tmp
	complex(ki), dimension(20,4) :: q2
	complex(ki), dimension(4)    :: e3, e4, qt, qtest
	complex(ki), dimension(0:19) :: c2
	complex(ki), dimension(-2:0) :: tot2tmp
	
	complex(ki), dimension(0:1)  :: m
	real(ki)                     :: K11, K12
	complex(ki), dimension(-2:0) :: MI2a, MI2b, MI2c, MI2d, MI2e, MI2J111
!	complex(ki), dimension(-2:0) :: B0p12, B0z11, B0z22
	
	interface
	   function     numeval(ncut, Q, mu2)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     complex(ki), dimension(4), intent(in) :: Q
	     complex(ki),               intent(in) :: mu2
	     complex(ki) :: numeval
	   end function numeval
	end interface

	interface
	   function     rank_numeval(ncut)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     integer :: rank_numeval
	   end function rank_numeval
	end interface
	

	1903 format(A4,I2,A1,I2,A3,1(D24.15),1(D24.15))
	9002 format(A3,I2,A1,I2,A7,D24.15,A1,D24.15,A3)
	call checksmatalloc(smatallocated)

	tot2sum=czip
	totr2sum=czip

        if (rnk.lt.0) goto 5


	diff=nleg-rnk
	if (diff.ge.3) then 
		do j=0,19
		c2(j)=czip
		enddo
                goto 10
        endif
	


 5      continue

		icut2=1
		do j2=1,nleg-1
		do j1=0,j2-1
			
			cut2=j2*10+j1
			m(0)=msq(j1)
			m(1)=msq(j2)
			
			do n1=1,4
				k1(n1)=Vi(j2,n1)-Vi(j1,n1)
				p0(n1)=Vi(j1,n1)
			enddo
			
			if (abs(abs(k1(4))-1.0_ki) .lt. 0.1_ki) then
				factor = 0.345_ki
			else
				factor = one
			end if
			
			k2(1)=-sign(factor/rt3,k1(1))
			k2(2)=-sign(factor/rt3,k1(2))
			k2(3)=-sign(factor/rt3,k1(3))
			k2(4)= sign(factor,k1(4))
			
			call getbase(k1,k2,r1,r2,e1,e2,e3,e4)
			
                        if (rnk.lt.0) then
                           rank=rank_numeval(cut2)
                           if (rank.lt.0.or.(nleg-rank.ge.3)) cycle
                        else
                           rank=rnk
                        endif

			call getq2(nleg,rank,cut2,e1,e2,e3,e4,p0,k1,msq,q2,qt)
			
			if (iresc.eq.1) then
				call getc2(numetens,nleg,rank,c2,cut2,q2,qt,Vi,msq)
			else
				call getc2(numeval,nleg,rank,c2,cut2,q2,qt,Vi,msq)
			endif
			
			highrank = ((c2(13) /= 0._ki) .or. (c2(10) /= 0._ki))
			
			if (itest.eq.2) call lnntest2(diff,numeval,cut2,c2,qt,p0,e2,e3,e4,ok)
			call store2(icut2,cut2,p0,e1,e2,e3,e4,c2)
			
			if (smatallocated) then
				call get2smat(cut2,real(m),K11)
			else
				K11 = sdot(k1,k1)
			end if
			K12=sdot(K1,K2)
			
			if (present(cache_flag)) then
				call getMI2(K11,m,scale2,MI2a,MI2b,MI2c,MI2d,MI2e,cache_flag, cache_offset, scalar_cache)
				if (highrank) call getMI2hr(K11,m,scale2,MI2J111)!,cache_flag, cache_offset, scalar_cache)
			else
				call getMI2(K11,m,scale2,MI2a,MI2b,MI2c,MI2d,MI2e)
				if (highrank) call getMI2hr(K11,m,scale2,MI2J111)
			end if
			
			if (highrank) then
				call add2(K11,K12,m,c2,MI2a,MI2b,MI2c,MI2d,MI2e,tot2tmp,totr2tmp,MI2J111)
			else
				call add2(K11,K12,m,c2,MI2a,MI2b,MI2c,MI2d,MI2e,tot2tmp,totr2tmp)
			endif
			
			if (verbosity.ge.2) then
				do ep=0,2
					write(iout,1903) 'B0 (',cut2,',',-ep,') =',real(MI2a(-ep)),aimag(MI2a(-ep))
					write(iout,1903) 'B1 (',cut2,',',-ep,') =',real(MI2b(-ep)),aimag(MI2b(-ep))
					!write(iout,1903) 'B00(',cut2,',',-ep,') =',real(MI2c(-ep)),aimag(MI2c(-ep))
					write(iout,1903) 'B00(',cut2,',',-ep,') =',real(MI2c(-ep)),aimag(MI2c(-ep))
					write(iout,1903) 'B01(',cut2,',',-ep,') =',real(MI2d(-ep)),aimag(MI2d(-ep))
					write(iout,1903) 'B11(',cut2,',',-ep,') =',real(MI2e(-ep)),aimag(MI2e(-ep))
					if (highrank) then
						write(iout,1903) 'B31(',cut2,',',-ep,') =',real(MI2J111(-ep)),aimag(MI2J111(-ep))
					end if
				end do
				if (highrank) then
					write(iout,1903) 'B06(',cut2,',',-1,') =',real(one/24.0_ki* (K11 - two*m(0)-four*m(1))),0.0_ki
				end if
			endif
	
			if(verbosity.gt.0)then
				do i=0,19
				write(iout,9002) 'c2(',cut2,',',i,')   = (',real(c2(i)),',',aimag(c2(i)),'  )'
				enddo
				write(iout,*)
			endif
			
			if (use_maccu) then
				call add_accu(accr_re, accr_im, totr2sum)
				call add_accu(acc_re(:), acc_im(:), tot2sum(:))
			else 
				tot2sum(:)=tot2sum(:)+tot2tmp(:)
				totr2sum=totr2sum+totr2tmp
			end if
			
			icut2=icut2+1
		enddo
		enddo
		nc2=icut2-1
 10 continue
end subroutine

subroutine samurai1plecut(numeval,nleg,rnk,Vi,msq,scale2,tot1sum,ok,cache_flag,scalar_cache, cache_offset,rank_numeval)
	use options, only: use_maccu
	use maccu
	use mgetc1
	implicit none
	
	integer,                            intent(in ) :: nleg, rnk
	real(ki),    dimension(0:nleg-1,4), intent(in ) :: Vi
	complex(ki), dimension(0:nleg-1),   intent(in ) :: msq
	real(ki),                           intent(in ) :: scale2
	complex(ki), dimension(-2:0),       intent(out) :: tot1sum
	logical,                            intent(out) :: ok
	logical,     intent(inout), optional                    :: cache_flag
	complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache
	integer,     intent(inout), optional                    :: cache_offset
	optional :: rank_numeval
	
	type(accumulator_type), dimension(-2:0) :: acc_re, acc_im
	type(accumulator_type) :: accr_re, accr_im

        integer :: rank	
	integer :: i,j,k,ep,j1,icut1
	integer :: cut1,n1
	integer :: diff 
	
	complex(ki) :: mu2, mu2test
	real(ki)    :: r1, r2 
	real(ki), dimension(4):: k1, k2, k3, e1, e2, p0, L3
	
	complex(ki), dimension(15,4) :: q1
	complex(ki), dimension(4)    :: e3, e4, qt, qtest
	complex(ki), dimension(0:15) :: c1
	complex(ki), dimension(-2:0) :: tot1tmp
	
	complex(ki) :: m0
	complex(ki) :: e3e4 
	
	interface
	   function     numeval(ncut, Q, mu2)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     complex(ki), dimension(4), intent(in) :: Q
	     complex(ki),               intent(in) :: mu2
	     complex(ki) :: numeval
	   end function numeval
	end interface
	
	interface
	   function     rank_numeval(ncut)
	     use precision
	     implicit none
	     integer,                   intent(in) :: ncut
	     integer :: rank_numeval
	   end function rank_numeval
	end interface

	complex(ki), dimension(-2:0) :: MI1

	9011 format(A3,I1,A1,I2,A5,D24.15,A1,D24.15,A3)
	9012 format(A3,I2,A1,I2,A7,D24.15,A1,D24.15,A3)

	tot1sum  = czip

        if (rnk.lt.0) goto 5


	diff=nleg-rnk
	if (diff.ge.2) then
		do j=0,15
		c1(j)=czip
		enddo
                goto 10
        endif
	


 5      continue



		icut1=1
		do j1=0,nleg-1
			
			cut1=j1
			m0=msq(j1)
			
			k1(1)=+one/rt3
			k1(2)=-one/rt3
			k1(3)=+one/rt3
			k1(4)=+two
			
			k2(1)=+one/rt2
			k2(2)=-one/rt2
			k2(3)=+one/rt2
			k2(4)=+rt3
			
			e3e4=sdot(e3,e4)
			
			do n1=1,4
				p0(n1)=Vi(j1,n1)
			enddo
			
			call getbase(k1,k2,r1,r2,e1,e2,e3,e4)    
			
                        if (rnk.lt.0) then
                           rank=rank_numeval(cut1)
                           if (rank.lt.0.or.(nleg-rank.ge.2)) cycle
                        else
                           rank=rnk
                        endif

			call getq1(nleg,rank,cut1,e1,e2,e3,e4,p0,msq,q1,qt)
			
			if (iresc.eq.1) then
				call getc1(numetens,nleg,rank,c1,cut1,q1,qt,Vi,msq)
			else
				call getc1(numeval,nleg,rank,c1,cut1,q1,qt,Vi,msq)
			endif
			
			if (itest.eq.2) call lnntest1(diff,numeval,cut1,c1,qt,p0,e1,e2,e3,e4,ok)
			call store1(icut1,cut1,p0,e1,e2,e3,e4,c1)
			
			if (present(cache_flag)) then
				call getMI1(m0,scale2,MI1,cache_flag, cache_offset, scalar_cache)
			else
				call getMI1(m0,scale2,MI1)
			end if
			
			call add1(m0,e3e4,c1,MI1,tot1tmp)
			
			if (verbosity.ge.2) then 
				do ep=-2,0
				write(iout,9011) 'I1(',cut1,',',ep,') = (',real(MI1(ep)),',',aimag(MI1(ep)),'  )'
				enddo
			endif

			if(verbosity.gt.0)then
				do i=0,15
					if(i .eq. 9) cycle
					write(iout,9012) 'c1(',cut1,',',i,')   = (',real(c1(i)),',',&
					& aimag(c1(i)),'  )'
				enddo
				write(iout,*)
			endif
			
			if (use_maccu) then
				call add_accu(acc_re, acc_im, tot1sum)
			else 
				tot1sum(:)=tot1sum(:)+tot1tmp(:)
			end if
			
			icut1=icut1+1
		enddo
		nc1=icut1-1
 10      continue
end subroutine samurai1plecut




end module msamurai

