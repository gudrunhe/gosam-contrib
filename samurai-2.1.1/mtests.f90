module mtests
   use precision, only: ki
   use constants
   use options
   use mfunctions
   use mrestore
   use save
   use ltest
   use mglobal, only: resit, denst
   implicit none

   private
   public :: pwtest, nntest, lnntest1, lnntest2, lnntest3, lnntest4

   real(ki), parameter :: nprec = epsilon(1.0_ki) * 1.0E+03_ki

   interface nntest
      module procedure nntest_rm
      module procedure nntest_cm
   end interface nntest

contains

   subroutine pwtest(nleg,rank,ok)
      implicit none
      integer, intent(in) :: nleg, rank
      logical, intent(out) :: ok

      integer :: i,j
      integer :: diff
      real(ki), dimension(4) :: e1, e2
      real(ki) :: pwval
      complex(ki), dimension(4) :: e3, e4
      complex(ki) :: pwt
      complex(ki), dimension(4) :: pwtv

      diff = nleg-rank
      pwt=czip
 
      if(verbosity.ge.3)then
          write(iout,*) '------------------------------------'
          write(iout,*) 'PoWer test'
      endif

      if (diff.eq.0) then
         !---#[ maximum rank :
         pwtv=(/ czip, czip, czip, czip /)
         
         do j =1,max1
            do i=1,4
               e1(i)=sav1e1(j,i)
               e2(i)=sav1e2(j,i)
               e3(i)=sav1e3(j,i)
               e4(i)=sav1e4(j,i)
            enddo
            do i=1,4
               pwtv(i)=pwtv(i)+savc1(j,1)*e1(i)+savc1(j,2)*e2(i) &
              &               +savc1(j,3)*e3(i)+savc1(j,4)*e4(i)
            enddo
         enddo
 
         pwval =  maxval(abs(pwtv))
         !---#] maximum rank :
      elseif (diff.eq.1) then
         !---#[ nleg-rank = 1:
         pwt=czip
         do i =1,6
            pwt=pwt+savc1(i,0)
         enddo
          
         pwval =  abs(pwt)
         !---#] nleg-rank = 1:
      elseif (diff.eq.2) then
         !---#[ nleg-rank = 2:
         pwt=czip
         do i =1,max2
            pwt=pwt+savc2(i,0)
         enddo
         pwval =  abs(pwt)
         !---#] nleg-rank = 2:
      else
         !---#[ nleg-rank > 2:
         pwt=czip
         do i =1,max3
            pwt=pwt+savc3(i,0)
         enddo
         pwval =  abs(pwt)
         !---#] nleg-rank > 2:
       endif

       if(verbosity.ge.3) write(iout,*) 'pwtest = ', pwval
       if (pwval.ge.pwlimit) then 
          ok=.false.
          if(verbosity.gt.0) then
             write(iout,*) ' POWERTEST FAILED! ' 
             write(iout,*) 'point above discarded'
          endif
       else
          ok=.true.
          if(verbosity.ge.3) write(iout,*) ' PowerTest passed! ' 
       endif
       if(verbosity.ge.3) write(iout,*) '------------------------------------'
   end subroutine

   subroutine nntest_rm(numeval,q1,mu2,nleg,Vi,msq,ok)
      implicit none

      complex(ki), dimension(4), intent(in) :: q1
      complex(ki), intent(in) :: mu2
      integer, intent(in) :: nleg
      real(ki), dimension(0:nleg-1) :: msq
      real(ki), dimension(0:nleg-1,4) :: Vi
      logical, intent(out) :: ok

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

      integer :: i, i1, i2, i3, i4, i5, ncut
      integer :: dicut5, dicut4, dicut3, dicut2, dicut1
      complex(ki) :: dens1,dens2,dens3,dens4,dens5,xneval
      complex(ki) :: resi5, resi4, resi3, resi2 ,resi1, resitot, ztest

      resi1=czip
      resi2=czip
      resi3=czip
      resi4=czip
      resi5=czip

      ztest=czip
      xneval=czip

      if (nleg.gt.5) then
         !---#[ Contributo dei pentuple cuts:
         dicut5=1

         do i5=4,nleg-1
            do i4=3,i5-1
               do i3=2,i4-1
                  do i2=1,i3-1
                     do i1=0,i2-1
                        dens5=cone
                        do i=0,nleg-1 
                           if ((i.ne.i1).and.(i.ne.i2) &
                           & .and.(i.ne.i3).and.(i.ne.i4)&
                           & .and.(i.ne.i5)) then 
                              dens5=dens5*denevalmu2(nleg,i,q1,Vi,msq,mu2)
                           endif
                        enddo
                        resi5=resi5+dens5*res5(dicut5,mu2)
                        dicut5=dicut5+1
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei pentuple cuts:
      elseif (nleg.eq.5) then
         resi5=res5(1,mu2)
      end if

      if (nleg.ge.5) then
         !---#[ Contributo dei quadruple cuts:
         dicut4=1
         do i4=3,nleg-1
            do i3=2,i4-1
               do i2=1,i3-1
                  do i1=0,i2-1
                     dens4=cone
                     do i=0,nleg-1
                        if ((i.ne.i1).and.(i.ne.i2) &
                        &.and.(i.ne.i3).and.(i.ne.i4)) then
                           dens4=dens4*denevalmu2(nleg,i,q1,Vi,msq,mu2)
                        endif
                     enddo
                     
                     resi4=resi4+dens4*Res4(dicut4,q1,mu2)
                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei quadruple cuts:
      elseif (nleg.eq.4) then
         resi4=Res4(1,q1,mu2)
      end if

      if (nleg.ge.4) then
         !---#[ Contributo dei Triple cuts:
         dicut3=1
         do i3=2,nleg-1
            do i2=1,i3-1
               do i1=0,i2-1
                  dens3=cone
                  do i=0,nleg-1
                     if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                        dens3=dens3*denevalmu2(nleg,i,q1,Vi,msq,mu2)
                     endif
                  enddo

                  resi3=resi3+dens3*Res3(dicut3,q1,mu2)
                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:
      elseif (nleg.eq.3) then
         resi3=Res3(1,q1,mu2)
      end if
       
      if (nleg.ge.3) then
         !---#[ Contributo dei Double cuts:
         dicut2=1
         do i2=1,nleg-1
            do i1=0,i2-1
               dens2=cone
               do i=0,nleg-1
                  if ((i.ne.i1).and.(i.ne.i2)) then
                     dens2=dens2*denevalmu2(nleg,i,q1,Vi,msq,mu2)
                  endif
               enddo

               resi2=resi2+dens2*Res2(dicut2,q1,mu2)
               dicut2=dicut2+1
            enddo
         enddo
         !---#] Contributo dei Double cuts:
      elseif (nleg.eq.2) then
         resi2=Res2(1,q1,mu2)
      endif

      !---#[ Contribution of the single cut:
      dicut1=1
      do i1=0,nleg-1
         dens1=cone
         ! ---> I replaced this loop
         !do i=0,nleg-1
         !   if (i.ne.i1) then
         !      dens1=dens1*denevalmu2(nleg,i,q1,Vi,msq,mu2)
         !   endif
         !enddo
         ! ---> because we can do this without an if statement inside:
         do i=0,i1-1
            dens1=dens1*denevalmu2(nleg,i,q1,Vi,msq,mu2)
         enddo
         do i=i1+1,nleg-1
            dens1=dens1*denevalmu2(nleg,i,q1,Vi,msq,mu2)
         enddo

         resi1=resi1+dens1*Res1(dicut1,q1)
         dicut1=dicut1+1
      enddo
      !---#] Contribution of the single cut:

      xneval=numeval(ncut,q1,mu2)
      resitot=resi5+resi4+resi3+resi2+resi1
      ztest=xneval-resitot

      if(verbosity.ge.3)then
         write(iout,*) '--------------------------------------------'
         write(iout,*) 'N=N test for'
         write(iout,*) 'q(0) = ', q1(4)
         write(iout,*) 'q(1) = ', q1(1)
         write(iout,*) 'q(2) = ', q1(2)
         write(iout,*) 'q(3) = ', q1(3)
         write(iout,*) 'and mu2 = ', mu2
         write(iout,*) 'N    calculated', xneval
         write(iout,*) 'N reconstructed', resitot
         write(iout,*) ' '
         write(iout,*) 'difference = ', ztest
         write(iout,*) 'rel.diff.  = ', ztest/xneval
         write(iout,*) '--------------------------------------------'
      endif
      
      if (abs(ztest/xneval).gt.nnlimit) then
         ok=.false.
         if (verbosity.gt.0) write(iout,*) 'N=N test FAILED'
      else
         ok=.true.
      endif
   end subroutine nntest_rm

   subroutine nntest_cm(numeval,q1,mu2,nleg,Vi,msq,ok)
      implicit none

      complex(ki), dimension(4), intent(in) :: q1
      complex(ki), intent(in) :: mu2
      integer, intent(in) :: nleg
      complex(ki), dimension(0:nleg-1) :: msq
      real(ki), dimension(0:nleg-1,4) :: Vi
      logical, intent(out) :: ok

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

      integer :: i, i1, i2, i3, i4, i5, ncut
      integer :: dicut5, dicut4, dicut3, dicut2, dicut1
      complex(ki) :: dens1,dens2,dens3,dens4,dens5,xneval
      complex(ki) :: resi5, resi4, resi3, resi2 ,resi1, resitot, ztest

      resi1=czip
      resi2=czip
      resi3=czip
      resi4=czip
      resi5=czip

      ztest=czip
      xneval=czip

      if (nleg.gt.5) then
         !---#[ Contributo dei pentuple cuts:
         dicut5=1

         do i5=4,nleg-1
            do i4=3,i5-1
               do i3=2,i4-1
                  do i2=1,i3-1
                     do i1=0,i2-1
                        dens5=cone
                        do i=0,nleg-1 
                           if ((i.ne.i1).and.(i.ne.i2) &
                           & .and.(i.ne.i3).and.(i.ne.i4)&
                           & .and.(i.ne.i5)) then 
                              dens5=dens5*denevalmu2(nleg,i,q1,Vi,msq,mu2)
                           endif
                        enddo
                        resi5=resi5+dens5*res5(dicut5,mu2)
                        dicut5=dicut5+1
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei pentuple cuts:
      elseif (nleg.eq.5) then
         resi5=res5(1,mu2)
      end if

      if (nleg.ge.5) then
         !---#[ Contributo dei quadruple cuts:
         dicut4=1
         do i4=3,nleg-1
            do i3=2,i4-1
               do i2=1,i3-1
                  do i1=0,i2-1
                     dens4=cone
                     do i=0,nleg-1
                        if ((i.ne.i1).and.(i.ne.i2) &
                        &.and.(i.ne.i3).and.(i.ne.i4)) then
                           dens4=dens4*denevalmu2(nleg,i,q1,Vi,msq,mu2)
                        endif
                     enddo
                     
                     resi4=resi4+dens4*Res4(dicut4,q1,mu2)
                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei quadruple cuts:
      elseif (nleg.eq.4) then
         resi4=Res4(1,q1,mu2)
      end if

      if (nleg.ge.4) then
         !---#[ Contributo dei Triple cuts:
         dicut3=1
         do i3=2,nleg-1
            do i2=1,i3-1
               do i1=0,i2-1
                  dens3=cone
                  do i=0,nleg-1
                     if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                        dens3=dens3*denevalmu2(nleg,i,q1,Vi,msq,mu2)
                     endif
                  enddo

                  resi3=resi3+dens3*Res3(dicut3,q1,mu2)
                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:
      elseif (nleg.eq.3) then
         resi3=Res3(1,q1,mu2)
      end if
       
      if (nleg.ge.3) then
         !---#[ Contributo dei Double cuts:
         dicut2=1
         do i2=1,nleg-1
            do i1=0,i2-1
               dens2=cone
               do i=0,nleg-1
                  if ((i.ne.i1).and.(i.ne.i2)) then
                     dens2=dens2*denevalmu2(nleg,i,q1,Vi,msq,mu2)
                  endif
               enddo

               resi2=resi2+dens2*Res2(dicut2,q1,mu2)
               dicut2=dicut2+1
            enddo
         enddo
         !---#] Contributo dei Double cuts:
      elseif (nleg.eq.2) then
         resi2=Res2(1,q1,mu2)
      endif

      !---#[ Contribution of the single cut:
      dicut1=1
      do i1=0,nleg-1
         dens1=cone
         ! ---> I replaced this loop
         !do i=0,nleg-1
         !   if (i.ne.i1) then
         !      dens1=dens1*denevalmu2(nleg,i,q1,Vi,msq,mu2)
         !   endif
         !enddo
         ! ---> because we can do this without an if statement inside:
         do i=0,i1-1
            dens1=dens1*denevalmu2(nleg,i,q1,Vi,msq,mu2)
         enddo
         do i=i1+1,nleg-1
            dens1=dens1*denevalmu2(nleg,i,q1,Vi,msq,mu2)
         enddo

         resi1=resi1+dens1*Res1(dicut1,q1)
         dicut1=dicut1+1
      enddo
      !---#] Contribution of the single cut:

      xneval=numeval(ncut,q1,mu2)
      resitot=resi5+resi4+resi3+resi2+resi1
      ztest=xneval-resitot

      if(verbosity.ge.3)then
         write(iout,*) '--------------------------------------------'
         write(iout,*) 'N=N test for'
         write(iout,*) 'q(0) = ', q1(4)
         write(iout,*) 'q(1) = ', q1(1)
         write(iout,*) 'q(2) = ', q1(2)
         write(iout,*) 'q(3) = ', q1(3)
         write(iout,*) 'and mu2 = ', mu2
         write(iout,*) 'N    calculated', xneval
         write(iout,*) 'N reconstructed', resitot
         write(iout,*) ' '
         write(iout,*) 'difference = ', ztest
         write(iout,*) 'rel.diff.  = ', ztest/xneval
         write(iout,*) '--------------------------------------------'
      endif
      
      if (abs(ztest/xneval).gt.nnlimit) then
         ok=.false.
         if (verbosity.gt.0) write(iout,*) 'N=N test FAILED'
      else
         ok=.true.
      endif
   end subroutine nntest_cm

   subroutine lnntest4(numeval,cut4,c4,qt,p0,k3,e3,e4,ok)
      implicit none
      integer, intent(in) :: cut4
      complex(ki), dimension(0:4), intent(in) :: c4
      complex(ki), dimension(4), intent(in) :: qt, e3, e4
      real(ki), dimension(4), intent(in) :: p0, k3
      logical, intent(out) :: ok

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

      complex(ki), dimension(4) :: pm
      complex(ki) :: test4, poli4, reldif

      pm(:)=qt(:)+p0(:)
      if     (imeth.eq.'diag') then
         test4=(numeval(cut4,qt,chaf)-resit(4))/denst(4)
      elseif (imeth.eq.'tree') then
         test4= numeval(cut4,qt,chaf)-resit(4) /denst(4)
      endif
      poli4=poly4(c4,pm,chaf,k3,e3,e4)

      if (abs(poli4).lt.nprec .and. abs(test4).lt.nprec) then
         reldif = nprec
      else
         reldif=(test4-poli4)/poli4
      endif

      if(verbosity.ge.3)then
         write(iout,*) ' LOCAL N=N TEST - BOX'
         write(iout,*) 'cut4    =',cut4
         write(iout,*) 'test4   =',test4
         write(iout,*) 'poli4   =',poli4
         write(iout,*) 'diff4   =',test4-poli4
         write(iout,*) 'rel.dif4=',abs(reldif)
      endif

      if (abs(reldif).gt.lnnlimit4) then
         ok = .false.
         if (verbosity.gt.0) &
            & write(iout,*) 'LOCAL N=N test FAILED for the BOX',cut4
      else
         ok = .true.
      endif
   end subroutine lnntest4

   subroutine lnntest3(numeval,cut3,c3,qt,p0,e3,e4,ok)
      implicit none
      integer, intent(in) :: cut3
      complex(ki), dimension(0:9), intent(in) :: c3
      complex(ki), dimension(4), intent(in) :: qt, e3, e4
      real(ki), dimension(4), intent(in) :: p0 
      logical, intent(out) :: ok  

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

      complex(ki), dimension(4) :: pm
      complex(ki) ::  test3, poli3, reldif

      pm(:)=qt(:)+p0(:)

      if     (imeth.eq.'diag') then
         test3=(numeval(cut3,qt,chaf)-resit(3)) / denst(3)
      elseif (imeth.eq.'tree') then
         test3= numeval(cut3,qt,chaf)-resit(3) / denst(3)
      endif
      poli3=poly3(c3,pm,chaf,e3,e4)

      if (abs(poli3).lt.nprec .and. abs(test3).lt.nprec) then
         reldif = nprec
      else
         reldif=(test3-poli3)/poli3
      endif

      if(verbosity.ge.3)then
         write(iout,*) ' LOCAL N=N TEST - triangle'
         write(iout,*) " cut3   =",cut3
         write(iout,*) "test3   =",test3
         write(iout,*) "poli3   =",poli3
         write(iout,*) "diff3   =",test3-poli3
         write(iout,*) "rel.dif3=",abs(reldif)
      endif
         
      if (abs(reldif).gt.lnnlimit3) then
         ok = .false.
         if (verbosity.gt.0) write(iout,*) &
                & 'LOCAL N=N test FAILED for the 3-cut',cut3
      else
         ok = .true.
      endif
   end subroutine lnntest3

   subroutine lnntest2(numeval,cut2,c2,qt,p0,e2,e3,e4,ok)
      implicit none
      integer, intent(in) :: cut2
      complex(ki), dimension(0:9), intent(in) :: c2
      complex(ki), dimension(4), intent(in) :: qt, e3, e4
      real(ki), dimension(4), intent(in) :: p0, e2
      logical, intent(out) :: ok

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

      complex(ki), dimension(4) :: pm
      complex(ki) :: test2, poli2, reldif

      pm(:)=qt(:)+p0(:)
      if     (imeth.eq.'diag') then
         test2=(numeval(cut2,qt,chaf)-resit(2))/denst(2)
      elseif (imeth.eq.'tree') then
         test2= numeval(cut2,qt,chaf)-resit(2) /denst(2)
      endif
      poli2=poly2(c2,pm,chaf,e2,e3,e4)

      if (abs(poli2).lt.nprec .and. abs(test2).lt.nprec) then
         reldif = nprec
      else
         reldif=(test2-poli2)/poli2
      endif
      
      if(verbosity.ge.3)then
         write(iout,*) ' LOCAL N=N TEST - Bubble'
         write(iout,*) " cut2   =",cut2
         write(iout,*) "test2   =",test2
         write(iout,*) "poli2   =",poli2
         write(iout,*) "diff2   =",test2-poli2
         write(iout,*) "rel.dif2=",abs(reldif)
      endif

      if (abs(reldif).gt.lnnlimit2) then
         ok = .false.
         if (verbosity.gt.0) write(iout,*) &
              &'LOCAL N=N test FAILED for the 2-cut',cut2
      else
         ok = .true.
      endif
   end subroutine lnntest2


   subroutine lnntest1(numeval,cut1,c1,qt,p0,e1,e2,e3,e4,ok)
      implicit none
      integer, intent(in) :: cut1
      complex(ki), dimension(0:4), intent(in) :: c1
      complex(ki), dimension(4), intent(in) :: qt, e3, e4
      real(ki), dimension(4), intent(in) :: p0, e1, e2
      logical, intent(out) :: ok

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

      complex(ki), dimension(4) :: pm
      complex(ki) :: test1, poli1, reldif
     
      pm(:)=qt(:)+p0(:)
      if     (imeth.eq.'diag') then
         test1=(numeval(cut1,qt,chaf)-resit(1)) / denst(1)
      elseif (imeth.eq.'tree') then
         test1= numeval(cut1,qt,chaf)-resit(1) /denst(1)
      endif
      poli1=poly1(c1,pm,e1,e2,e3,e4)

      if (abs(poli1).lt.nprec .and. abs(test1).lt.nprec) then
         reldif = nprec
      else
         reldif=(test1-poli1)/poli1
      endif

      if (verbosity.ge.3) then
         write(iout,*) ' LOCAL N=N TEST - Tadpole'         
         write(iout,*) " cut1   =",cut1
         write(iout,*) "test1   =",test1
         write(iout,*) "poli1   =",poli1
         write(iout,*) "diff1   =",test1-poli1
         write(iout,*) "rel.dif1=",abs(reldif)
      endif

      if (abs(reldif).gt.lnnlimit1) then
         ok = .false.
         if (verbosity.gt.0) write(iout,*) &
              &'LOCAL N=N test FAILED for the 1-cut',cut1
      else
         ok = .true.
      endif
   end subroutine lnntest1

end module mtests
