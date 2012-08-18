module mgetc1
   use precision, only: ki
   use constants
   use options
   use mfunctions
   use mrestore
   implicit none

   private

   interface getc1
      module procedure getc1_rm
      module procedure getc1_cm
   end interface getc1

   public :: getc1

contains

   subroutine getc1_cm(numeval,nleg,rank,c1,cut1,q1,qt,Vi,msq)
      use mglobal, only: G0c, mu2g, MP12, mu2t, resit, denst, mu2test
      implicit none
      integer, intent(in) :: nleg, rank, cut1
      complex(ki), dimension(0:4), intent(out) :: c1
      complex(ki), dimension(5,4), intent(in) :: q1
      complex(ki), dimension(4), intent(in) :: qt
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      complex(ki), dimension(0:nleg-1), intent(in) :: msq

      integer :: i,m,n,j1,i1,i2,i3,i4,i5
      integer :: dicut5,dicut4,dicut3,dicut2,diff
      complex(ki), dimension(5) :: dens1,dens2,dens3,dens4,dens5,xneval
      complex(ki), dimension(0:4) :: f1
      complex(ki), dimension(5) :: resi5, resi4, resi3, resi2, known
      complex(ki) :: dens2t,dens3t,dens4t,dens5t
      logical evalres

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

      mu2test(1)=mu2t(1)

      j1=cut1

      resi2(:)=czip
      resi3(:)=czip
      resi4(:)=czip
      resi5(:)=czip
      known(:)=czip
      xneval(:)=czip
      dens1(:)=cone

      !---  for lnntest
      resit(1)=czip
      denst(1)=cone

     !---  for simplified sampling
     diff = nleg-rank

     if (diff.eq.1) then
        c1(1)=czip
        c1(2)=czip
        c1(3)=czip
        c1(4)=czip
         
        if (nleg.eq.5) then
           resi5(1)=res5(1,mu2g(1))
           resit(1)=res5(1,mu2t(1))
           goto 11
        elseif (nleg.eq.4) then
           resi4(1)=Res4(1,q1(1,:),mu2g(1))
           resit(1)=Res4(1,qt,mu2t(1))
           goto 21
        elseif (nleg.eq.3) then
           resi3(1)=Res3(1,q1(1,:),mu2g(1))
           resit(1)=Res3(1,qt,mu2t(1))
           goto 31
        elseif (nleg.eq.2) then
           resi2(1)=Res2(1,q1(1,:),mu2g(1))
           resit(1)=Res2(1,qt,mu2t(1))
           goto 41
         else
            !---#[ Contributo dei pentuple cuts:
            dicut5=1
            do i5=4,nleg-1
               do i4=3,i5-1
                  do i3=2,i4-1
                     do i2=1,i3-1
                        do i1=0,i2-1
                           dens5(1)=cone
                           dens5t=cone
                           evalres=.false.
                           
                           loop_10: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)&
                                      &.and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i).eq.(j1)) then
                                    dens5(1)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_10
                                 else
                                    dens5(1)=dens5(1)*denevalmu2(nleg,i,&
                                            &q1(1,:),Vi,msq,mu2g(1))
                                    dens5t=dens5t*denevalmu2(nleg,i,qt,Vi,msq,&
                                            &mu2t(1))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_10

                           if (evalres) then
                              resi5(1)=resi5(1)+dens5(1)*res5(dicut5,mu2g(1))
                              resit(1)=resit(1)+dens5(1)*res5(dicut5,mu2t(1))
                           endif
         
                           dicut5=dicut5+1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !---#] Contributo dei pentuple cuts:
         endif

 11      continue

         !---#[ Contributo dei quadruple cuts:
         dicut4=1
         do i4=3,nleg-1
            do i3=2,i4-1
               do i2=1,i3-1
                  do i1=0,i2-1
                     dens4(1)=cone
                     dens4t=cone
                     evalres=.false.
      
                     loop_20: do i=0,nleg-1
                        if ((i.ne.i1).and.(i.ne.i2) &
                                &  .and.(i.ne.i3).and.(i.ne.i4)) then
                           if (i.eq.j1) then
                              dens4(1)=czip
                              dens4t=czip
                              evalres=.false.
                              exit loop_20
                           else
                              dens4(1)=dens4(1)*denevalmu2(nleg,i,q1(1,:),&
                                       &Vi,msq,mu2g(1))
                              dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,&
                                       &mu2t(1))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_20

                     if (evalres) then
                        resi4(1)=resi4(1)+dens4(1)*Res4(dicut4,q1(1,:),mu2g(1))
                        resit(1)=resit(1)+dens4t*Res4(dicut4,qt,mu2t(1))
                     endif

                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei quadruple cuts:
 21      continue

         !---#[ Contributo dei Triple cuts:
         dicut3=1
         do i3=2,nleg-1
            do i2=1,i3-1
               do i1=0,i2-1
                  dens3(1)=cone
                  dens3t=cone
                  evalres=.false.
      
                  loop_30: do i=0,nleg-1
                     if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                        if (i.eq.j1) then
                           dens3(1)=czip
                           dens3t=czip
                           evalres=.false.
                           exit loop_30
                        else
                           dens3(1)=dens3(1)*denevalmu2(nleg,i,q1(1,:),&
                              &Vi,msq,mu2g(1))
                           dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                           evalres=.true.
                        endif
                     endif
                  enddo loop_30

                  if (evalres) then
                     resi3(1)=resi3(1)+dens3(1)*Res3(dicut3,q1(1,:),mu2g(1))
                     resit(1)=resit(1)+dens3t*Res3(dicut3,qt,mu2t(1))
                  endif

                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:

 31      continue

         !---#[ Contributo dei Double cuts:
         dicut2=1
         do i2=1,nleg-1
            do i1=0,i2-1
               dens2(1)=cone
               dens2t=cone
               evalres=.false.
      
               loop_40: do i=0,nleg-1
                  if ((i.ne.i1).and.(i.ne.i2)) then
                     if (i.eq.j1) then
                        dens2(1)=czip
                        dens2t=czip
                        evalres=.false.
                        exit loop_40
                     else
                        dens2(1)=dens2(1)*denevalmu2(nleg,i,q1(1,:),&
                             &Vi,msq,mu2g(1))
                        dens2t=dens2t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                        evalres=.true.
                     endif
                  endif
               enddo loop_40

               if (evalres) then
                  resi2(1)=resi2(1)+dens2(1)*Res2(dicut2,q1(1,:),mu2g(1))
                  resit(1)=resit(1)+dens2t*Res2(dicut2,qt,mu2t(1))
               endif

               dicut2=dicut2+1
            enddo
         enddo
         !---#] Contributo dei Double cuts:

 41      continue

         !---
         do i=0,nleg-1
            if (i.ne.j1) then
               dens1(1)=dens1(1)*denevalmu2(nleg,i,q1(1,:),Vi,msq,mu2g(1))
               denst(1)=denst(1)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
            endif
         enddo

         xneval(1)=numeval(cut1,q1(1,:),mu2g(1))

         if     (imeth.eq.'diag') then
            known(1)=(xneval(1)-resi5(1)-resi4(1)-resi3(1)-resi2(1))/dens1(1) 
         elseif (imeth.eq.'tree') then
            known(1)=xneval(1)-(resi5(1)+resi4(1)+resi3(1)+resi2(1))/dens1(1) 
         endif

         c1(0) = known(1)
      else
         !--- Decomposizione standard

         if (nleg.eq.5) then
            resi5(:)=res5(1,mu2g(1))
            resit(1)=res5(1,mu2t(1)) 
            goto 111
         elseif (nleg.eq.4) then
            do n=1,5         
               resi4(n)=Res4(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res4(1,qt,mu2t(1))
            goto 121
         elseif (nleg.eq.3) then
            do n=1,5
               resi3(n)=Res3(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res3(1,qt,mu2t(1))
            goto 131
         elseif (nleg.eq.2) then
            do n=1,5
               resi2(n)=Res2(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res2(1,qt,mu2t(1))
            goto 141
         else
            !---#[ Contributo dei pentuple cuts:
            dicut5=1
            do i5=4,nleg-1
               do i4=3,i5-1
                  do i3=2,i4-1
                     do i2=1,i3-1
                        do i1=0,i2-1
                           dens5(:)=cone
                           dens5t=cone
                           evalres=.false.
      
                           loop_110: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)&
                                     &.and.(i.ne.i4).and.(i.ne.i5)) then
                                 if (i.eq.j1) then
                                    dens5(:)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_110
                                 else
                                    do n=1,5
                                       dens5(n)=dens5(n)*denevalmu2(nleg,i,&
                                           &q1(n,:),Vi,msq,mu2g(1))
                                    enddo
                                    dens5t=dens5t*denevalmu2(nleg,&
                                           &i,qt,Vi,msq,mu2t(1))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_110

                           if (evalres) then
                              resi5(:)=resi5(:)+dens5(:)*res5(dicut5,mu2g(1))
                              resit(1)=resit(1)+dens5t*res5(dicut5,mu2t(1))
                           endif

                           dicut5=dicut5+1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !---#] Contributo dei pentuple cuts:
         endif

 111     continue

         !---#[ Contributo dei quadruple cuts:
         dicut4=1
         do i4=3,nleg-1
            do i3=2,i4-1
               do i2=1,i3-1
                  do i1=0,i2-1
                     dens4(:)=cone
                     dens4t=cone
                     evalres=.false.
      
                     loop_120: do i=0,nleg-1
                        if ((i.ne.i1).and.(i.ne.i2) &
                               &  .and.(i.ne.i3).and.(i.ne.i4)) then
                           if (i.eq.j1) then
                              dens4(:)=czip
                              dens4t=czip
                              evalres=.false.
                              exit loop_120
                           else
                              do n=1,5
                                 dens4(n)=dens4(n)*denevalmu2(nleg,i,&
                                      &q1(n,:),Vi,msq,mu2g(1))
                              enddo
                              dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_120

                     if (evalres) then
                        do n=1,5
                           resi4(n)=resi4(n)+dens4(n)*&
                               &Res4(dicut4,q1(n,:),mu2g(1))
                        enddo
                        resit(1)=resit(1)+dens4t*Res4(dicut4,qt,mu2t(1))
                     endif
                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei quadruple cuts:

 121     continue

         !---#[ Contributo dei Triple cuts:
         dicut3=1
         do i3=2,nleg-1
            do i2=1,i3-1
               do i1=0,i2-1
                  dens3(:)=cone
                  dens3t=cone
                  evalres=.false.
      
                  loop_130: do i=0,nleg-1
                     if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                        if (i.eq.j1) then
                           dens3(:)=czip
                           dens3t=czip
                           evalres=.false.
                           exit loop_130
                        else
                           do n=1,5
                              dens3(n)=dens3(n)*denevalmu2(nleg,i,q1(n,:),&
                                      &Vi,msq,mu2g(1))
                           enddo
                           dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                           evalres=.true.
                        endif
                     endif
                  enddo loop_130

                  if (evalres) then
                     do n=1,5
                        resi3(n)=resi3(n)+dens3(n)*Res3(dicut3,q1(n,:),mu2g(1))
                     enddo
                     resit(1)=resit(1)+dens3t*Res3(dicut3,qt,mu2t(1))
                  endif

                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:

 131     continue

         !---#[ Contributo dei Double cuts:
         dicut2=1
         do i2=1,nleg-1
            do i1=0,i2-1
               dens2(:)=cone
               dens2t=cone
               evalres=.false.
   
               loop_140: do i=0,nleg-1
                  if ((i.ne.i1).and.(i.ne.i2)) then
                     if (i.eq.j1) then
                        dens2(:)=czip
                        dens2t=czip
                        evalres=.false.
                        exit loop_140
                     else
                        do n=1,5
                           dens2(n)=dens2(n)*denevalmu2(nleg,i,q1(n,:),&
                              &Vi,msq,mu2g(1))
                        enddo
                        dens2t=dens2t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                        evalres=.true.
                     endif
                  endif
               enddo loop_140

               if (evalres) then
                  do n=1,5
                     resi2(n)=resi2(n)+dens2(n)*Res2(dicut2,q1(n,:),mu2g(1))
                  enddo
                  resit(1)=resit(1)+dens2t*Res2(dicut2,qt,mu2t(1))
               endif
               dicut2=dicut2+1
            enddo
         enddo

         !---#] Contributo dei Double cuts:
 141     continue

         !---

         do i=0,nleg-1
            if (i.ne.j1) then
               do n=1,5
                  dens1(n)=dens1(n)*denevalmu2(nleg,i,q1(n,:),Vi,msq,mu2g(1))
               enddo
               denst(1)=denst(1)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
            endif
         enddo

         do n=1,5
            xneval(n)=numeval(cut1,q1(n,:),mu2g(1))
         enddo

         if     (imeth.eq.'diag') then
            known(:)=(xneval(:)-resi5(:)-resi4(:)-resi3(:)-resi2(:))/dens1(:)
         elseif (imeth.eq.'tree') then
            known(:)=xneval(:)-(resi5(:)+resi4(:)+resi3(:)+resi2(:))/dens1(:)
         endif

         do m=0,4
            f1(m)=known(m+1)
         enddo

         !--- Coefficienti
         c1(0) = (f1(0)+f1(1))/2.0_ki
         c1(1) = (-f1(0)-f1(1)+f1(3)+f1(4))/(2.0_ki*MP12(1))
         c1(2) = -(-2.0_ki*f1(0) + f1(3) + f1(4))/(2.0_ki*G0c*MP12(1))
         c1(3) = (2.0_ki*f1(0)-2.0_ki*f1(2)-f1(3)+f1(4))/(2.0_ki*G0c*MP12(1))
         c1(4) = (f1(0)-f1(2))/MP12(1)
      endif
   end subroutine getc1_cm

   subroutine getc1_rm(numeval,nleg,rank,c1,cut1,q1,qt,Vi,msq)
      use mglobal, only: G0, mu2g, MP12, mu2t, resit, denst, mu2test
      implicit none
      integer, intent(in) :: nleg, rank, cut1
      complex(ki), dimension(0:4), intent(out) :: c1
      complex(ki), dimension(5,4), intent(in) :: q1
      complex(ki), dimension(4), intent(in) :: qt
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      real(ki), dimension(0:nleg-1), intent(in) :: msq

      integer :: i,m,n,j1,i1,i2,i3,i4,i5
      integer :: dicut5,dicut4,dicut3,dicut2,diff
      complex(ki), dimension(5) :: dens1,dens2,dens3,dens4,dens5,xneval
      complex(ki), dimension(0:4) :: f1
      complex(ki), dimension(5) :: resi5, resi4, resi3, resi2, known
      complex(ki) :: dens2t,dens3t,dens4t,dens5t
      logical evalres

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

      mu2test(1)=mu2t(1)

      j1=cut1

      resi2(:)=czip
      resi3(:)=czip
      resi4(:)=czip
      resi5(:)=czip
      known(:)=czip
      xneval(:)=czip
      dens1(:)=cone

      !---  for lnntest
      resit(1)=czip
      denst(1)=cone

     !---  for simplified sampling
     diff = nleg-rank

     if (diff.eq.1) then
        c1(1)=czip
        c1(2)=czip
        c1(3)=czip
        c1(4)=czip
         
        if (nleg.eq.5) then
           resi5(1)=res5(1,mu2g(1))
           resit(1)=res5(1,mu2t(1))
           goto 11
        elseif (nleg.eq.4) then
           resi4(1)=Res4(1,q1(1,:),mu2g(1))
           resit(1)=Res4(1,qt,mu2t(1))
           goto 21
        elseif (nleg.eq.3) then
           resi3(1)=Res3(1,q1(1,:),mu2g(1))
           resit(1)=Res3(1,qt,mu2t(1))
           goto 31
        elseif (nleg.eq.2) then
           resi2(1)=Res2(1,q1(1,:),mu2g(1))
           resit(1)=Res2(1,qt,mu2t(1))
           goto 41
         else
            !---#[ Contributo dei pentuple cuts:
            dicut5=1
            do i5=4,nleg-1
               do i4=3,i5-1
                  do i3=2,i4-1
                     do i2=1,i3-1
                        do i1=0,i2-1
                           dens5(1)=cone
                           dens5t=cone
                           evalres=.false.
                           
                           loop_10: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)&
                                      &.and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i).eq.(j1)) then
                                    dens5(1)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_10
                                 else
                                    dens5(1)=dens5(1)*denevalmu2(nleg,i,&
                                            &q1(1,:),Vi,msq,mu2g(1))
                                    dens5t=dens5t*denevalmu2(nleg,i,qt,Vi,msq,&
                                            &mu2t(1))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_10

                           if (evalres) then
                              resi5(1)=resi5(1)+dens5(1)*res5(dicut5,mu2g(1))
                              resit(1)=resit(1)+dens5(1)*res5(dicut5,mu2t(1))
                           endif
         
                           dicut5=dicut5+1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !---#] Contributo dei pentuple cuts:
         endif

 11      continue

         !---#[ Contributo dei quadruple cuts:
         dicut4=1
         do i4=3,nleg-1
            do i3=2,i4-1
               do i2=1,i3-1
                  do i1=0,i2-1
                     dens4(1)=cone
                     dens4t=cone
                     evalres=.false.
      
                     loop_20: do i=0,nleg-1
                        if ((i.ne.i1).and.(i.ne.i2) &
                                &  .and.(i.ne.i3).and.(i.ne.i4)) then
                           if (i.eq.j1) then
                              dens4(1)=czip
                              dens4t=czip
                              evalres=.false.
                              exit loop_20
                           else
                              dens4(1)=dens4(1)*denevalmu2(nleg,i,q1(1,:),&
                                       &Vi,msq,mu2g(1))
                              dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,&
                                       &mu2t(1))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_20

                     if (evalres) then
                        resi4(1)=resi4(1)+dens4(1)*Res4(dicut4,q1(1,:),mu2g(1))
                        resit(1)=resit(1)+dens4t*Res4(dicut4,qt,mu2t(1))
                     endif

                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei quadruple cuts:
 21      continue

         !---#[ Contributo dei Triple cuts:
         dicut3=1
         do i3=2,nleg-1
            do i2=1,i3-1
               do i1=0,i2-1
                  dens3(1)=cone
                  dens3t=cone
                  evalres=.false.
      
                  loop_30: do i=0,nleg-1
                     if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                        if (i.eq.j1) then
                           dens3(1)=czip
                           dens3t=czip
                           evalres=.false.
                           exit loop_30
                        else
                           dens3(1)=dens3(1)*denevalmu2(nleg,i,q1(1,:),&
                              &Vi,msq,mu2g(1))
                           dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                           evalres=.true.
                        endif
                     endif
                  enddo loop_30

                  if (evalres) then
                     resi3(1)=resi3(1)+dens3(1)*Res3(dicut3,q1(1,:),mu2g(1))
                     resit(1)=resit(1)+dens3t*Res3(dicut3,qt,mu2t(1))
                  endif

                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:

 31      continue

         !---#[ Contributo dei Double cuts:
         dicut2=1
         do i2=1,nleg-1
            do i1=0,i2-1
               dens2(1)=cone
               dens2t=cone
               evalres=.false.
      
               loop_40: do i=0,nleg-1
                  if ((i.ne.i1).and.(i.ne.i2)) then
                     if (i.eq.j1) then
                        dens2(1)=czip
                        dens2t=czip
                        evalres=.false.
                        exit loop_40
                     else
                        dens2(1)=dens2(1)*denevalmu2(nleg,i,q1(1,:),&
                             &Vi,msq,mu2g(1))
                        dens2t=dens2t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                        evalres=.true.
                     endif
                  endif
               enddo loop_40

               if (evalres) then
                  resi2(1)=resi2(1)+dens2(1)*Res2(dicut2,q1(1,:),mu2g(1))
                  resit(1)=resit(1)+dens2t*Res2(dicut2,qt,mu2t(1))
               endif

               dicut2=dicut2+1
            enddo
         enddo
         !---#] Contributo dei Double cuts:

 41      continue

         !---
         do i=0,nleg-1
            if (i.ne.j1) then
               dens1(1)=dens1(1)*denevalmu2(nleg,i,q1(1,:),Vi,msq,mu2g(1))
               denst(1)=denst(1)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
            endif
         enddo

         xneval(1)=numeval(cut1,q1(1,:),mu2g(1))

         if     (imeth.eq.'diag') then
            known(1)=(xneval(1)-resi5(1)-resi4(1)-resi3(1)-resi2(1))/dens1(1) 
         elseif (imeth.eq.'tree') then
            known(1)=xneval(1)-(resi5(1)+resi4(1)+resi3(1)+resi2(1))/dens1(1) 
         endif

         c1(0) = known(1)
      else
         !--- Decomposizione standard

         if (nleg.eq.5) then
            resi5(:)=res5(1,mu2g(1))
            resit(1)=res5(1,mu2t(1)) 
            goto 111
         elseif (nleg.eq.4) then
            do n=1,5         
               resi4(n)=Res4(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res4(1,qt,mu2t(1))
            goto 121
         elseif (nleg.eq.3) then
            do n=1,5
               resi3(n)=Res3(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res3(1,qt,mu2t(1))
            goto 131
         elseif (nleg.eq.2) then
            do n=1,5
               resi2(n)=Res2(1,q1(n,:),mu2g(1))
            enddo
            resit(1)=Res2(1,qt,mu2t(1))
            goto 141
         else
            !---#[ Contributo dei pentuple cuts:
            dicut5=1
            do i5=4,nleg-1
               do i4=3,i5-1
                  do i3=2,i4-1
                     do i2=1,i3-1
                        do i1=0,i2-1
                           dens5(:)=cone
                           dens5t=cone
                           evalres=.false.
      
                           loop_110: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)&
                                     &.and.(i.ne.i4).and.(i.ne.i5)) then
                                 if (i.eq.j1) then
                                    dens5(:)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_110
                                 else
                                    do n=1,5
                                       dens5(n)=dens5(n)*denevalmu2(nleg,i,&
                                           &q1(n,:),Vi,msq,mu2g(1))
                                    enddo
                                    dens5t=dens5t*denevalmu2(nleg,&
                                           &i,qt,Vi,msq,mu2t(1))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_110

                           if (evalres) then
                              resi5(:)=resi5(:)+dens5(:)*res5(dicut5,mu2g(1))
                              resit(1)=resit(1)+dens5t*res5(dicut5,mu2t(1))
                           endif

                           dicut5=dicut5+1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !---#] Contributo dei pentuple cuts:
         endif

 111     continue

         !---#[ Contributo dei quadruple cuts:
         dicut4=1
         do i4=3,nleg-1
            do i3=2,i4-1
               do i2=1,i3-1
                  do i1=0,i2-1
                     dens4(:)=cone
                     dens4t=cone
                     evalres=.false.
      
                     loop_120: do i=0,nleg-1
                        if ((i.ne.i1).and.(i.ne.i2) &
                               &  .and.(i.ne.i3).and.(i.ne.i4)) then
                           if (i.eq.j1) then
                              dens4(:)=czip
                              dens4t=czip
                              evalres=.false.
                              exit loop_120
                           else
                              do n=1,5
                                 dens4(n)=dens4(n)*denevalmu2(nleg,i,&
                                      &q1(n,:),Vi,msq,mu2g(1))
                              enddo
                              dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_120

                     if (evalres) then
                        do n=1,5
                           resi4(n)=resi4(n)+dens4(n)*&
                               &Res4(dicut4,q1(n,:),mu2g(1))
                        enddo
                        resit(1)=resit(1)+dens4t*Res4(dicut4,qt,mu2t(1))
                     endif
                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo
         !---#] Contributo dei quadruple cuts:

 121     continue

         !---#[ Contributo dei Triple cuts:
         dicut3=1
         do i3=2,nleg-1
            do i2=1,i3-1
               do i1=0,i2-1
                  dens3(:)=cone
                  dens3t=cone
                  evalres=.false.
      
                  loop_130: do i=0,nleg-1
                     if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                        if (i.eq.j1) then
                           dens3(:)=czip
                           dens3t=czip
                           evalres=.false.
                           exit loop_130
                        else
                           do n=1,5
                              dens3(n)=dens3(n)*denevalmu2(nleg,i,q1(n,:),&
                                      &Vi,msq,mu2g(1))
                           enddo
                           dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                           evalres=.true.
                        endif
                     endif
                  enddo loop_130

                  if (evalres) then
                     do n=1,5
                        resi3(n)=resi3(n)+dens3(n)*Res3(dicut3,q1(n,:),mu2g(1))
                     enddo
                     resit(1)=resit(1)+dens3t*Res3(dicut3,qt,mu2t(1))
                  endif

                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:

 131     continue

         !---#[ Contributo dei Double cuts:
         dicut2=1
         do i2=1,nleg-1
            do i1=0,i2-1
               dens2(:)=cone
               dens2t=cone
               evalres=.false.
   
               loop_140: do i=0,nleg-1
                  if ((i.ne.i1).and.(i.ne.i2)) then
                     if (i.eq.j1) then
                        dens2(:)=czip
                        dens2t=czip
                        evalres=.false.
                        exit loop_140
                     else
                        do n=1,5
                           dens2(n)=dens2(n)*denevalmu2(nleg,i,q1(n,:),&
                              &Vi,msq,mu2g(1))
                        enddo
                        dens2t=dens2t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
                        evalres=.true.
                     endif
                  endif
               enddo loop_140

               if (evalres) then
                  do n=1,5
                     resi2(n)=resi2(n)+dens2(n)*Res2(dicut2,q1(n,:),mu2g(1))
                  enddo
                  resit(1)=resit(1)+dens2t*Res2(dicut2,qt,mu2t(1))
               endif
               dicut2=dicut2+1
            enddo
         enddo

         !---#] Contributo dei Double cuts:
 141     continue

         !---

         do i=0,nleg-1
            if (i.ne.j1) then
               do n=1,5
                  dens1(n)=dens1(n)*denevalmu2(nleg,i,q1(n,:),Vi,msq,mu2g(1))
               enddo
               denst(1)=denst(1)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(1))
            endif
         enddo

         do n=1,5
            xneval(n)=numeval(cut1,q1(n,:),mu2g(1))
         enddo

         if     (imeth.eq.'diag') then
            known(:)=(xneval(:)-resi5(:)-resi4(:)-resi3(:)-resi2(:))/dens1(:)
         elseif (imeth.eq.'tree') then
            known(:)=xneval(:)-(resi5(:)+resi4(:)+resi3(:)+resi2(:))/dens1(:)
         endif

         do m=0,4
            f1(m)=known(m+1)
         enddo

         !--- Coefficienti
         c1(0) = (f1(0)+f1(1))/2.0_ki
         c1(1) = (-f1(0)-f1(1)+f1(3)+f1(4))/(2.0_ki*MP12(1))
         c1(2) = -(-2.0_ki*f1(0) + f1(3) + f1(4))/(2.0_ki*G0*MP12(1))
         c1(3) = (2.0_ki*f1(0)-2.0_ki*f1(2)-f1(3)+f1(4))/(2.0_ki*G0*MP12(1))
         c1(4) = (f1(0)-f1(2))/MP12(1)
      endif
   end subroutine getc1_rm

end module mgetc1

