module mgetc2
   use precision, only: ki
   use constants
   use options
   use mfunctions
   use mrestore
   implicit none

   private

   interface getc2
      module procedure getc2_rm
      module procedure getc2_cm
   end interface getc2

   public :: getc2

contains

   subroutine getc2_rm(numeval,nleg,rank,c2,cut2,q2,qt,Vi,msq)
      use mglobal, only: Fp,Fz,Fm,F1z,KB,mu2g,KK,MP12,mu2t,&
                       & resit,denst,mu2test
      implicit none
      integer, intent(in) :: nleg,rank,cut2
      complex(ki), dimension(0:9), intent(out) :: c2
      complex(ki), dimension(10,4), intent(in) :: q2
      complex(ki), dimension(4), intent(in) :: qt
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi

      integer :: i,m,n,j1,j2,i1,i2,i3,i4,i5
      integer :: dicut5,dicut4,dicut3,mx,diff,nsol
      complex(ki), dimension(10) :: dens2,dens3,dens4,dens5,xneval
      complex(ki), dimension(10) :: resi5,resi4,resi3,known
      complex(ki), dimension(0:9) :: f2
      complex(ki) :: dens3t,dens4t,dens5t
      logical evalres
      !real(ki) :: Fppow2, Fppow3, Fppow4
      !real(ki) :: Fzpow2, Fzpow3, Fzpow4, Fzpow5, Fzpow6, tmp1
      !---#[ HAGGIES:
      complex(ki) :: t1
      complex(ki) :: t2
      complex(ki) :: t3
      complex(ki) :: t4
      complex(ki) :: t5
      complex(ki) :: t6
      complex(ki) :: t7
      complex(ki) :: t8
      complex(ki) :: t9
      complex(ki) :: t10
      complex(ki) :: t11
      complex(ki) :: t12
      complex(ki) :: t13
      complex(ki) :: t14
      complex(ki) :: t15
      complex(ki) :: t16
      complex(ki) :: t17
      complex(ki) :: t18
      complex(ki) :: t19
      complex(ki) :: t20
      complex(ki) :: t21
      complex(ki) :: t22
      complex(ki) :: t23
      complex(ki) :: t24
      complex(ki) :: t25
      complex(ki) :: t26
      complex(ki) :: t27
      complex(ki) :: t28
      complex(ki) :: t29
      complex(ki) :: t30
      complex(ki) :: t31
      complex(ki) :: t32
      complex(ki) :: t33
      complex(ki) :: t34
      complex(ki) :: t35
      complex(ki) :: t36
      complex(ki) :: t37
      complex(ki) :: t38
      complex(ki) :: t39
      complex(ki) :: t40
      complex(ki) :: t41
      complex(ki) :: t42
      complex(ki) :: t43
      complex(ki) :: t44
      complex(ki) :: t45
      !---#] HAGGIES:
      complex(ki), dimension(10) :: mu2vec

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

      mu2test(2) = mu2t(2)
     
      j2=cut2/10
      j1=cut2-j2*10

      resi3(:)=czip
      resi4(:)=czip
      resi5(:)=czip
      known(:)=czip
      xneval(:)=czip
      dens2(:)=cone

!---  for lnntest
      resit(2)=czip
      denst(2)=cone


!---  for simplified sampling
      diff = nleg-rank
      
      if (diff.eq.2) then
         !---#[ simplified sampling -- only c2(0):
         if     (nleg.eq.5) then
            resi5(1)=res5(1,czip)
            resit(2)=res5(1,mu2t(2))
            goto 11
         elseif (nleg.eq.4) then
            resi4(1)=Res4(1,q2(1,:),czip)    
            resit(2)   =res4(1,qt,mu2t(2))
            goto 21
         elseif (nleg.eq.3) then
            resi3(1)=Res3(1,q2(1,:),czip)
            resit(2)=Res3(1,qt,mu2t(2))
            goto 31
         elseif (nleg.eq.2) then
            goto 36
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
                                      &  .and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i.eq.j1).or.(i.eq.j2)) then
                                    dens5(1)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_10
                                 else
                                    dens5(1)=dens5(1)&
                                       &*denevalmu2(nleg,i,q2(1,:),Vi,msq,czip)
                                    dens5t=dens5t&
                                       &*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_10

                           if (evalres) then
                              resi5(1)=resi5(1)+dens5(1)*res5(dicut5,czip)
                              resit(2)=resit(2)+dens5t*res5(dicut5,mu2t(2))
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
                        if ((i.ne.i1).and.(i.ne.i2).and.&
                                 &(i.ne.i3).and.(i.ne.i4)) then
                           if ((i.eq.j1).or.(i.eq.j2)) then
                              dens4(1)=czip
                              dens4t=czip
                              evalres=.false.
                              exit loop_20
                           else
                              dens4(1)=dens4(1)*denevalmu2(nleg,i,q2(1,:),&
                                        &Vi,msq,czip)
                              dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_20
   
   
                     if (evalres) then
                        resi4(1)=resi4(1)+dens4(1)*Res4(dicut4,q2(1,:),czip)    
                        resit(2)=resit(2)+dens4t*Res4(dicut4,qt,mu2t(2))    
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
                        if ((i.eq.j1).or.(i.eq.j2)) then
                           dens3(1)=czip
                           dens3t=czip
                           evalres=.false.
                           exit loop_30
                        else
                           dens3(1)=dens3(1)*denevalmu2(nleg,i,q2(1,:),&
                                  &Vi,msq,czip)
                           dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
                           evalres=.true.
                        endif
                     endif
                  enddo loop_30
   
                  if (evalres) then
                     resi3(1)=resi3(1)+dens3(1)*Res3(dicut3,q2(1,:),czip)
                     resit(2)=resit(2)+dens3t*Res3(dicut3,qt,mu2t(2))
                  endif
   
                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:

 31      continue

         do i=0,nleg-1
            if ((i.ne.j1).and.(i.ne.j2)) then
               dens2(1)=dens2(1)*denevalmu2(nleg,i,q2(1,:),Vi,msq,czip)
               denst(2)=denst(2)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
            endif
         enddo

 36      continue
   
         xneval(1)=numeval(cut2,q2(1,:),czip)
       
         if     (imeth.eq.'diag') then
            known(1)=(xneval(1)-resi5(1)-resi4(1)-resi3(1))/dens2(1)
         elseif (imeth.eq.'tree') then
            known(1)=xneval(1)-(resi5(1)+resi4(1)+resi3(1))/dens2(1)
         endif
        
         c2(0)=known(1)
         do m=1,9
            c2(m)=czip
         enddo
   
         !---#] simplified sampling -- only c2(0):
      else
         !---#[ Decomposizione standard:
         if (diff.eq.1) then
            ! rank1 c-system: 4 coefficients
            mu2vec=(/czip,czip,czip,czip,czip,czip,czip,czip,czip,czip/)
            nsol=4
            
         else
           ! traditional system
            mu2vec=(/czip,czip,czip,czip,czip,czip,czip,czip,czip,mu2g(2)/)
            nsol=10
         endif
       
         if     (nleg.eq.5) then
            do n=1,nsol
               resi5(n)=Res5(1,mu2vec(n))
            enddo
            resit(2)=Res5(1,mu2t(2))
            goto 111
         elseif (nleg.eq.4) then
            do n=1,nsol
               resi4(n)=Res4(1,q2(n,:),mu2vec(n))
            enddo
            resit(2)=Res4(1,qt,mu2t(2))
            goto 121
         elseif (nleg.eq.3) then
            do n=1,nsol
               resi3(n)=Res3(1,q2(n,:),mu2vec(n))
            enddo
            resit(2)=Res3(1,qt,mu2t(2))
            goto 131
         elseif (nleg.eq.2) then
            goto 136
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
                              if ((i.ne.i1).and.(i.ne.i2) &
                                & .and.(i.ne.i3).and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i.eq.j1).or.(i.eq.j2)) then
                                    dens5(:)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_110
                                 else
                                    do n=1,nsol
                                       dens5(n)=dens5(n)*denevalmu2(nleg,i,&
                                              &q2(n,:),Vi,msq,mu2vec(n))
                                    enddo
                                    dens5t=dens5t*denevalmu2(nleg,i,qt,Vi,msq,&
                                             &mu2t(2))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_110
      
                           if (evalres) then
                              do n=1,nsol
                                 resi5(n)=resi5(n)+dens5(n)*res5(dicut5,mu2vec(n))
                              enddo
                              resit(2)=resit(2)+dens5t*res5(dicut5,mu2t(2))
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
                         if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3).and.&
                                  &(i.ne.i4)) then
                            if ((i.eq.j1).or.(i.eq.j2)) then
                               dens4(:)=czip
                               dens4t=czip
                               evalres=.false.
                               exit loop_120
                            else
                               do n=1,nsol
                                  dens4(n)=dens4(n)*denevalmu2(nleg,i,&
                                      &q2(n,:),Vi,msq,mu2vec(n))
                               enddo
                               dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,&
                                      &mu2t(2))
                               evalres=.true.
                            endif
                         endif
                      enddo loop_120
      
                      if (evalres) then
                         do n=1,nsol
                            resi4(n)=resi4(n)+dens4(n)*Res4(dicut4,q2(n,:),&
                                &mu2vec(n))
                         enddo
                         resit(2)=resit(2)+dens4t*Res4(dicut4,qt,mu2t(2))
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
          loop_dicut3: do i3=2,nleg-1
             do i2=1,i3-1
                do i1=0,i2-1
                   dens3(:)=cone
                   dens3t=cone
                   evalres=.false.
                   loop_130: do i=0,nleg-1
                      if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                         if ((i.eq.j1).or.(i.eq.j2)) then
                            dens3(:)=czip
                            dens3t=czip
                            evalres=.false.
                            exit loop_130
                         else
                            do n=1,nsol
                               dens3(n)=dens3(n)&
                               & *denevalmu2(nleg,i,q2(n,:),Vi,msq,mu2vec(n))
                            enddo
                            dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
                            evalres=.true.
                         endif
                      endif
                   enddo loop_130
   
                   if (evalres) then
                      do n=1,nsol
                         resi3(n)=resi3(n)+dens3(n)*Res3(dicut3,q2(n,:),&
                            &mu2vec(n))
                      enddo
                      resit(2)=resit(2)+dens3t*Res3(dicut3,qt,mu2t(2))
                   endif
   
                   dicut3=dicut3+1
                enddo
             enddo
          enddo loop_dicut3
          !---#] Contributo dei Triple cuts:

 131     continue
   
         do i=0,nleg-1
            if ((i.ne.j1).and.(i.ne.j2)) then
               do n=1,nsol
                  dens2(n)=dens2(n)*denevalmu2(nleg,i,q2(n,:),Vi,msq,mu2vec(n))
               enddo
               denst(2)=denst(2)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
            endif
         enddo

 136     continue

         do n=1,nsol
            xneval(n)=numeval(cut2,q2(n,:),mu2vec(n))
         enddo

         if     (imeth.eq.'diag') then
            known(:)=(xneval(:)-resi5(:)-resi4(:)-resi3(:))/dens2(:)
         elseif (imeth.eq.'tree') then
            known(:)=xneval(:)-(resi5(:)+resi4(:)+resi3(:))/dens2(:)
         endif

         if (diff.eq.1) then
            ! rank1 c-system: 4 coefficients
        
            do m=0,1
               f2(m)=effe(known,1,2,m)
            enddo
  
            do m=2,3
               mx=m-2
               f2(m)=effe(known,3,2,mx)
            enddo
  
            !---#[ getc2_S1:
            !----#[ original code:
!           c2(0) = f2(0)
!           c2(1) = (-f2(0) + f2(2))/(KB*MP12(2))
!           c2(3) = (KK(2)*(f2(1) - Fz*f2(3)))/((-one + Fp*Fz)*MP12(2))
!           c2(5) = (-(Fp*MP12(2)*c2(3)) - KK(2)*f2(3))/(KK(2)**2*MP12(2))
            !----#] original code:
            !----#[ HAGGIES:
            t1 = f2(0)
            c2(0) = t1
            c2(1) = ((f2(2)-t1)/(KB*MP12(2)))
            t1 = f2(3)
            c2(3) = ((f2(1)-t1*Fz)/((Fp*Fz-one)*MP12(2))*KK(2))
            c2(5) = ((-(Fp*MP12(2)*c2(3)+t1*KK(2)))/(KK(2)*KK(2)*MP12(2)))
            !----#] HAGGIES:
            !---#] getc2_S1:
         else
            ! traditional system
            do m=0,2
              f2(m)=effe(known,1,3,m)
            enddo
            do m=3,4
               mx=m-3
               f2(m)=effe(known,4,2,mx)
            enddo
            do m=5,6
               mx=m-5
               f2(m)=effe(known,6,2,mx)
            enddo
            f2(7)=known(8)
            f2(8)=known(9)
            f2(9)=known(10)
  
            !---#[ getc2_S2:
            !----#[ original code:
!           !--- coefficienti
!           Fppow2 = Fp*Fp
!           Fppow3 = Fppow2*Fp
!           Fppow4 = Fppow3*Fp
!           Fzpow2 = Fz*Fz
!           Fzpow3 = Fzpow2*Fz
!           Fzpow4 = Fzpow3*Fz
!           Fzpow5 = Fzpow4*Fz
!           Fzpow6 = Fzpow5*Fz
!           tmp1 = one - Fppow2
!  c2(2) = (0.5_ki*(Fppow3*((one + Fzpow4)*f2(0) + f2(1) - f2(3) + Fz*&
!    &(Fz*(f2(1) + (one + Fzpow2)*f2(2) + Fz*(-Fz*f2(3) - f2(4))) - f2(4)))&
!    &+ Fppow4*(-f2(0) - Fzpow2*f2(1) - Fzpow4*f2(2) + f2(3) + Fzpow3*f2(4))&
!    &+ Fm**2*(tmp1*f2(0) + tmp1*Fzpow2*f2(1) + tmp1*Fzpow4*f2(2) + (-tmp1)&
!    &*f2(3) + (-tmp1)*Fzpow3*f2(4)) + two*((-one - Fzpow2 - Fzpow4 + Fzpow6)&
!    &*f2(0) + (-one - Fzpow4)*f2(1) - f2(2) + Fz*(-Fz*f2(2) + Fz*f2(3)&
!    &+ Fzpow3*f2(3) + f2(4) + Fzpow4*f2(4))) + f2(5) + f2(6) + Fm*((one&
!    &+ Fzpow4 + Fzpow3*two + Fppow2*(-one - Fzpow4 - Fzpow3*two))*f2(0)&
!    &+ (one + Fzpow2 + Fzpow5*two + Fppow2*(-one - Fzpow2 - Fzpow5*two))*f2(1)&
!    &+ tmp1*Fzpow2*f2(2) + (-tmp1)*f2(3) + Fzpow4*(tmp1*f2(2) + (-tmp1)*f2(3))&
!    &- two*f2(4) + Fppow2*two*f2(4) - tmp1*Fzpow3*(two*f2(3) + f2(4)) + Fz&
!    &*(tmp1*two*f2(2) - tmp1*f2(4)) - f2(5) - Fp*f2(6) + Fzpow6*(f2(5) + Fp&
!    &*f2(6) - f2(7)) + f2(7)) + Fp*((-one - Fzpow4)*f2(0) - f2(1) + f2(3)&
!    &+ f2(5) - f2(7) + Fz*(f2(4) + Fz*(-f2(1) - (one + Fzpow2)*f2(2) + Fz&
!    &*(Fz*f2(3) + f2(4) + Fzpow3*(-f2(5) + f2(7)))))) + Fzpow6*(-f2(5) - f2(6)&
!    &- f2(8)) + f2(8) + Fppow2*((three + Fzpow2*(one + Fzpow2 - Fzpow4)*two)&
!    &*f2(0) + two*f2(1) + two*f2(2) - f2(3) + Fzpow4*(two*f2(1) + f2(2) - two&
!    &*f2(3)) + Fzpow2*(f2(1) + two*f2(2) - two*f2(3)) - Fzpow3*f2(4) - Fz*two&
!    &*f2(4) - Fzpow5*two*f2(4) - f2(5) - f2(8) + Fzpow6*(f2(5) + f2(8)))))/&
!    &((-tmp1)*(-one + Fzpow6)*MP12(2)**2)
!  c2(7) = (KK(2)*(f2(2) + Fppow2*(-(Fzpow4*f2(1)) - f2(2) + Fzpow2*(-f2(0)&
!    &+ f2(3)) + Fzpow5*f2(4)) + Fppow3*(-((one + Fzpow4)*f2(0))&
!    &- (one + Fzpow2)*(f2(1) + Fzpow2*f2(2)) + (one + Fzpow4)*f2(3) &
!    &+ (Fz + Fzpow3)*f2(4)) - f2(6) + Fzpow2*(f2(0) - f2(3) +&
!    &Fzpow2*(f2(1) + Fz*(-f2(4) + Fz*f2(6)))) +&
!    &Fp*((one + Fzpow4)*f2(0) + f2(1) - f2(3) - f2(5) +&
!    &Fz*(Fz*(f2(1) + f2(2)) + Fzpow3*(f2(2) - f2(3)) -&
!    &f2(4) - Fzpow2*f2(4) + Fzpow5*(f2(5) - f2(7))) +&
!    &f2(7))))/((-tmp1)*(-one + Fzpow6)*MP12(2)**2)
!  c2(8) = ((-tmp1)*(one + Fzpow3 + Fzpow4)*f2(0) +&
!    &(-tmp1)*(one + Fzpow2 + Fzpow5)*f2(1) + (-tmp1)*Fzpow2*f2(2) +&
!    &(-tmp1)*Fzpow4*(f2(2) - f2(3)) + f2(3) + (-tmp1)*Fz*(f2(2) - f2(4))&
!    &+ f2(4) - (-tmp1)*Fzpow3*(f2(3) + f2(4)) + f2(5) +&
!    &Fp*(-(Fp*(f2(3) + f2(4))) + f2(6)) - Fzpow6*(f2(5) + Fp*f2(6) - f2(7))&
!    &- f2(7))/ ((-tmp1)*(-one + Fzpow6)*KK(2)*MP12(2)**2)
!  c2(4) = -((KK(2)**2*(f2(1) + Fzpow2*f2(2) + Fzpow4*(f2(0) - f2(3)) -&
!    &Fz*f2(4)))/((-one + Fzpow6)*MP12(2)**2))
!  c2(5) = (Fzpow3*f2(0) + Fzpow5*f2(1) + Fz*f2(2) - Fzpow3*f2(3) -&
!    &f2(4))/(KK(2)*MP12(2) - Fzpow6*KK(2)*MP12(2))
!  c2(3) = (KK(2)*(f2(2) + Fzpow2*(f2(0) - f2(3) + Fzpow2*(f2(1)&
!  - Fz*f2(4)))))/((-one + Fzpow6)*MP12(2))
!  c2(6) = (f2(0) - f2(3) + Fzpow2*(f2(1) + Fz*(Fz*f2(2) - f2(4))))/&
!    &((-one + Fzpow6)*KK(2)**2*MP12(2)**2)
!  c2(0) = f2(0)
!  c2(9) = (KK(2)*MP12(2)*c2(3) - MP12(2)**2*c2(4)&
!    &+ F1z*KK(2)**3*MP12(2)*c2(5) - F1z**2*KK(2)**4*MP12(2)**2*c2(6) &
!    &+ KK(2)**2*(-c2(0) + f2(9)))/(KK(2)**2*mu2g(2))
!  c2(1) = -(c2(3)/KK(2)) - Fm*KK(2)*c2(5) + (MP12(2)*(c2(4) + KK(2)*(c2(7)&
!    &+ KK(2)*(c2(2) + Fm*KK(2)*(Fm*KK(2)*c2(6) + c2(8))))))/KK(2)**2&
!    &+ (c2(0) - f2(8))/MP12(2)
            !----#] original code:
            !----#[ HAGGIES:
            t1 = (1.0_ki)+Fp
            t2 = Fz*Fz
            t3 = t2*Fz
            t4 = t3-(1.0_ki)
            t5 = (1.0_ki)+t3
            t6 = f2(2)
            t7 = t2*t2
            t8 = -(t7+(1.0_ki))
            t9 = f2(1)
            t10 = t3*t3
            t11 = f2(0)
            t12 = f2(4)
            t13 = t6*Fz
            t14 = f2(3)
            t15 = t14*Fz
            t16 = t15*t2
            t17 = f2(5)
            t18 = f2(6)
            t19 = f2(8)
            t20 = t9-t14
            t21 = (1.0_ki)+t7
            t22 = t11*t21
            t23 = (1.0_ki)+t2
            t24 = t23*t6
            t25 = Fp*Fp
            t26 = t14-t11
            t27 = t12*Fz
            t28 = t2*t27
            t29 = ((1.0_ki)-Fp)*t1
            t30 = t14*t29
            t31 = t29*t6
            t32 = t12*t29
            t33 = f2(7)
            t34 = t6+t9
            t35 = t27*t7
            t36 = t6-t14
            t37 = t17-t33
            t38 = t18*Fp
            t39 = (t37+t38)*t10
            t40 = t7*Fz
            t41 = (2.0_ki)*t40
            t42 = (2.0_ki)*t3
            t43 = t2*t31
            t44 = MP12(2)*MP12(2)
            t45 = t29*t4*t44*t5
            c2(2) = (-1.0_ki)*((2.0_ki)*((t10-(t7+t2+(1.0_ki)))*t11+(t12+t15+t&
            &16+t12*t7-t13)*Fz+t8*t9-t6)+t17+t18+t19+(-(t19+t18+t17))*t10+(t26&
            &+t28-(t6*t7+t2*t9))*t25*t25+(t11*t29+t31*t7+t2*t29*t9-(t2*t32*Fz+&
            &t30))*Fm*Fm+(t14+t17+(t12+((t12+t15+(t33-t17)*t2*Fz)*Fz-(t9+t24))&
            &*Fz)*Fz+t11*t8-(t9+t33))*Fp+((2.0_ki)*(t34-(t35+t27))+((3.0_ki)+(&
            &2.0_ki)*((1.0_ki)+t2-t7)*t2)*t11+(t17+t19)*t10+((2.0_ki)*t20+t6)*&
            &t7+((2.0_ki)*t36+t9)*t2-(t28+t19+t17+t14))*t25+((2.0_ki)*(t12*t25&
            &-t12)+t33+t39+t43+(t31-t30)*t7+((2.0_ki)*t31-t32)*Fz+((1.0_ki)+t2&
            &+t41+(-(t41+t2+(1.0_ki)))*t25)*t9+((1.0_ki)+t42+t7+(-(t7+t42+(1.0&
            &_ki)))*t25)*t11-((t12+(2.0_ki)*t14)*t2*t29*Fz+t38+t30+t17))*Fm+(t&
            &20+t22+((t24+t9+(-(t15+t12))*Fz)*Fz-t12)*Fz)*t25*Fp)/t45*0.5_ki
            t8 = t11-t14
            t10 = t2*t6
            c2(7) = (-1.0_ki)*(t6+(t8+(t9+(t18*Fz-t12)*Fz)*t2)*t2+(t35+t2*t26-&
            &(t7*t9+t6))*t25+(t20+t22+t33+(t34*Fz+t2*t36*Fz+t37*t7*Fz-(t12*t2+&
            &t12))*Fz-t17)*Fp+((t3+Fz)*t12+t14*t21-((t10+t9)*t23+t22))*t25*Fp-&
            &t18)/t45*KK(2)
            t14 = t12+t14
            c2(8) = ((t14+t17+(t18-t14*Fp)*Fp+t14*t2*t29*Fz-(t29*t36*t7+((1.0_&
            &ki)+t3+t7)*t11*t29+((1.0_ki)+t2+t40)*t29*t9+(t6-t12)*t29*Fz+t43+t&
            &39+t33))/((Fp-(1.0_ki))*t1*t4*t44*t5*KK(2)))
            t1 = KK(2)*KK(2)
            t3 = t4*t5
            c2(4) = (-1.0_ki)*(t10+t9+t7*t8-t27)/(t3*t44)*t1
            c2(5) = (-1.0_ki)*(t13+t11*t2*Fz+t7*t9*Fz-(t16+t12))/(t3*KK(2)*MP1&
            &2(2))
            c2(3) = ((t6+(t8+(t9-t27)*t2)*t2)/(t3*MP12(2))*KK(2))
            t4 = KK(2)*MP12(2)
            c2(6) = ((t8+(t9+(t13-t12)*Fz)*t2)/(t3*t4*t4))
            c2(0) = t11
            t2 = t1*F1z*MP12(2)
            c2(9) = (((f2(9)-c2(0))*t1+t4*c2(3)+t1*F1z*KK(2)*MP12(2)*c2(5)-(t4&
            &4*c2(4)+t2*t2*c2(6)))/(t1*mu2g(2)))
            t2 = Fm*KK(2)
            c2(1) = ((c2(0)-t19)/MP12(2)+(c2(4)+(c2(7)+(c2(2)+(c2(8)+t2*c2(6))&
            &*Fm*KK(2))*KK(2))*KK(2))/t1*MP12(2)-(c2(3)/KK(2)+t2*c2(5)))
            !----#] HAGGIES:
            !---#] getc2_S2:
         endif
         !---#] Decomposizione standard:
      endif

      if (diff.ge.1) then
         c2(2)=czip
         c2(4)=czip
         c2(6)=czip
         c2(7)=czip
         c2(8)=czip
         c2(9)=czip
         if (diff.ge.2) then
            c2(1)=czip
            c2(3)=czip
            c2(5)=czip
         endif             
      endif
   end subroutine getc2_rm

   subroutine getc2_cm(numeval,nleg,rank,c2,cut2,q2,qt,Vi,msq)
      use mglobal, only: Fpc,Fzc,Fmc,F1zc,KB,mu2g,KK,MP12,mu2t,&
                       & resit,denst,mu2test
      implicit none
      integer, intent(in) :: nleg,rank,cut2
      complex(ki), dimension(0:9), intent(out) :: c2
      complex(ki), dimension(10,4), intent(in) :: q2
      complex(ki), dimension(4), intent(in) :: qt
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi

      integer :: i,m,n,j1,j2,i1,i2,i3,i4,i5
      integer :: dicut5,dicut4,dicut3,mx,diff,nsol
      complex(ki), dimension(10) :: dens2,dens3,dens4,dens5,xneval
      complex(ki), dimension(10) :: resi5,resi4,resi3,known
      complex(ki), dimension(0:9) :: f2
      complex(ki) :: dens3t,dens4t,dens5t
      logical evalres
      !real(ki) :: Fpcpow2, Fpcpow3, Fpcpow4
      !real(ki) :: Fzcpow2, Fzcpow3, Fzcpow4, Fzcpow5, Fzcpow6, tmp1
      !---#[ HAGGIES:
      complex(ki) :: t1
      complex(ki) :: t2
      complex(ki) :: t3
      complex(ki) :: t4
      complex(ki) :: t5
      complex(ki) :: t6
      complex(ki) :: t7
      complex(ki) :: t8
      complex(ki) :: t9
      complex(ki) :: t10
      complex(ki) :: t11
      complex(ki) :: t12
      complex(ki) :: t13
      complex(ki) :: t14
      complex(ki) :: t15
      complex(ki) :: t16
      complex(ki) :: t17
      complex(ki) :: t18
      complex(ki) :: t19
      complex(ki) :: t20
      complex(ki) :: t21
      complex(ki) :: t22
      complex(ki) :: t23
      complex(ki) :: t24
      complex(ki) :: t25
      complex(ki) :: t26
      complex(ki) :: t27
      complex(ki) :: t28
      complex(ki) :: t29
      complex(ki) :: t30
      complex(ki) :: t31
      complex(ki) :: t32
      complex(ki) :: t33
      complex(ki) :: t34
      complex(ki) :: t35
      complex(ki) :: t36
      complex(ki) :: t37
      complex(ki) :: t38
      complex(ki) :: t39
      complex(ki) :: t40
      complex(ki) :: t41
      complex(ki) :: t42
      complex(ki) :: t43
      complex(ki) :: t44
      complex(ki) :: t45
      !---#] HAGGIES:
      complex(ki), dimension(10) :: mu2vec

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

      mu2test(2) = mu2t(2)
     
      j2=cut2/10
      j1=cut2-j2*10

      resi3(:)=czip
      resi4(:)=czip
      resi5(:)=czip
      known(:)=czip
      xneval(:)=czip
      dens2(:)=cone

!---  for lnntest
      resit(2)=czip
      denst(2)=cone


!---  for simplified sampling
      diff = nleg-rank
      
      if (diff.eq.2) then
         !---#[ simplified sampling -- only c2(0):
         if     (nleg.eq.5) then
            resi5(1)=res5(1,czip)
            resit(2)=res5(1,mu2t(2))
            goto 11
         elseif (nleg.eq.4) then
            resi4(1)=Res4(1,q2(1,:),czip)    
            resit(2)   =res4(1,qt,mu2t(2))
            goto 21
         elseif (nleg.eq.3) then
            resi3(1)=Res3(1,q2(1,:),czip)
            resit(2)=Res3(1,qt,mu2t(2))
            goto 31
         elseif (nleg.eq.2) then
            goto 36
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
                                      &  .and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i.eq.j1).or.(i.eq.j2)) then
                                    dens5(1)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_10
                                 else
                                    dens5(1)=dens5(1)&
                                       &*denevalmu2(nleg,i,q2(1,:),Vi,msq,czip)
                                    dens5t=dens5t&
                                       &*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_10

                           if (evalres) then
                              resi5(1)=resi5(1)+dens5(1)*res5(dicut5,czip)
                              resit(2)=resit(2)+dens5t*res5(dicut5,mu2t(2))
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
                        if ((i.ne.i1).and.(i.ne.i2).and.&
                                 &(i.ne.i3).and.(i.ne.i4)) then
                           if ((i.eq.j1).or.(i.eq.j2)) then
                              dens4(1)=czip
                              dens4t=czip
                              evalres=.false.
                              exit loop_20
                           else
                              dens4(1)=dens4(1)*denevalmu2(nleg,i,q2(1,:),&
                                        &Vi,msq,czip)
                              dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_20
   
   
                     if (evalres) then
                        resi4(1)=resi4(1)+dens4(1)*Res4(dicut4,q2(1,:),czip)    
                        resit(2)=resit(2)+dens4t*Res4(dicut4,qt,mu2t(2))    
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
                        if ((i.eq.j1).or.(i.eq.j2)) then
                           dens3(1)=czip
                           dens3t=czip
                           evalres=.false.
                           exit loop_30
                        else
                           dens3(1)=dens3(1)*denevalmu2(nleg,i,q2(1,:),&
                                  &Vi,msq,czip)
                           dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
                           evalres=.true.
                        endif
                     endif
                  enddo loop_30
   
                  if (evalres) then
                     resi3(1)=resi3(1)+dens3(1)*Res3(dicut3,q2(1,:),czip)
                     resit(2)=resit(2)+dens3t*Res3(dicut3,qt,mu2t(2))
                  endif
   
                  dicut3=dicut3+1
               enddo
            enddo
         enddo
         !---#] Contributo dei Triple cuts:

 31      continue

         do i=0,nleg-1
            if ((i.ne.j1).and.(i.ne.j2)) then
               dens2(1)=dens2(1)*denevalmu2(nleg,i,q2(1,:),Vi,msq,czip)
               denst(2)=denst(2)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
            endif
         enddo

 36      continue
   
         xneval(1)=numeval(cut2,q2(1,:),czip)
       
         if     (imeth.eq.'diag') then
            known(1)=(xneval(1)-resi5(1)-resi4(1)-resi3(1))/dens2(1)
         elseif (imeth.eq.'tree') then
            known(1)=xneval(1)-(resi5(1)+resi4(1)+resi3(1))/dens2(1)
         endif
        
         c2(0)=known(1)
         do m=1,9
            c2(m)=czip
         enddo
   
         !---#] simplified sampling -- only c2(0):
      else
         !---#[ Decomposizione standard:
         if (diff.eq.1) then
            ! rank1 c-system: 4 coefficients
            mu2vec=(/czip,czip,czip,czip,czip,czip,czip,czip,czip,czip/)
            nsol=4
            
         else
           ! traditional system
            mu2vec=(/czip,czip,czip,czip,czip,czip,czip,czip,czip,mu2g(2)/)
            nsol=10
         endif
       
         if     (nleg.eq.5) then
            do n=1,nsol
               resi5(n)=Res5(1,mu2vec(n))
            enddo
            resit(2)=Res5(1,mu2t(2))
            goto 111
         elseif (nleg.eq.4) then
            do n=1,nsol
               resi4(n)=Res4(1,q2(n,:),mu2vec(n))
            enddo
            resit(2)=Res4(1,qt,mu2t(2))
            goto 121
         elseif (nleg.eq.3) then
            do n=1,nsol
               resi3(n)=Res3(1,q2(n,:),mu2vec(n))
            enddo
            resit(2)=Res3(1,qt,mu2t(2))
            goto 131
         elseif (nleg.eq.2) then
            goto 136
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
                              if ((i.ne.i1).and.(i.ne.i2) &
                                & .and.(i.ne.i3).and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i.eq.j1).or.(i.eq.j2)) then
                                    dens5(:)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_110
                                 else
                                    do n=1,nsol
                                       dens5(n)=dens5(n)*denevalmu2(nleg,i,&
                                              &q2(n,:),Vi,msq,mu2vec(n))
                                    enddo
                                    dens5t=dens5t*denevalmu2(nleg,i,qt,Vi,msq,&
                                             &mu2t(2))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_110
      
                           if (evalres) then
                              do n=1,nsol
                                 resi5(n)=resi5(n)+dens5(n)*res5(dicut5,mu2vec(n))
                              enddo
                              resit(2)=resit(2)+dens5t*res5(dicut5,mu2t(2))
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
                         if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3).and.&
                                  &(i.ne.i4)) then
                            if ((i.eq.j1).or.(i.eq.j2)) then
                               dens4(:)=czip
                               dens4t=czip
                               evalres=.false.
                               exit loop_120
                            else
                               do n=1,nsol
                                  dens4(n)=dens4(n)*denevalmu2(nleg,i,&
                                      &q2(n,:),Vi,msq,mu2vec(n))
                               enddo
                               dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,&
                                      &mu2t(2))
                               evalres=.true.
                            endif
                         endif
                      enddo loop_120
      
                      if (evalres) then
                         do n=1,nsol
                            resi4(n)=resi4(n)+dens4(n)*Res4(dicut4,q2(n,:),&
                                &mu2vec(n))
                         enddo
                         resit(2)=resit(2)+dens4t*Res4(dicut4,qt,mu2t(2))
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
          loop_dicut3: do i3=2,nleg-1
             do i2=1,i3-1
                do i1=0,i2-1
                   dens3(:)=cone
                   dens3t=cone
                   evalres=.false.
                   loop_130: do i=0,nleg-1
                      if ((i.ne.i1).and.(i.ne.i2).and.(i.ne.i3)) then
                         if ((i.eq.j1).or.(i.eq.j2)) then
                            dens3(:)=czip
                            dens3t=czip
                            evalres=.false.
                            exit loop_130
                         else
                            do n=1,nsol
                               dens3(n)=dens3(n)&
                               & *denevalmu2(nleg,i,q2(n,:),Vi,msq,mu2vec(n))
                            enddo
                            dens3t=dens3t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
                            evalres=.true.
                         endif
                      endif
                   enddo loop_130
   
                   if (evalres) then
                      do n=1,nsol
                         resi3(n)=resi3(n)+dens3(n)*Res3(dicut3,q2(n,:),&
                            &mu2vec(n))
                      enddo
                      resit(2)=resit(2)+dens3t*Res3(dicut3,qt,mu2t(2))
                   endif
   
                   dicut3=dicut3+1
                enddo
             enddo
          enddo loop_dicut3
          !---#] Contributo dei Triple cuts:

 131     continue
   
         do i=0,nleg-1
            if ((i.ne.j1).and.(i.ne.j2)) then
               do n=1,nsol
                  dens2(n)=dens2(n)*denevalmu2(nleg,i,q2(n,:),Vi,msq,mu2vec(n))
               enddo
               denst(2)=denst(2)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(2))
            endif
         enddo

 136     continue

         do n=1,nsol
            xneval(n)=numeval(cut2,q2(n,:),mu2vec(n))
         enddo
        
         if     (imeth.eq.'diag') then
            known(:)=(xneval(:)-resi5(:)-resi4(:)-resi3(:))/dens2(:)
         elseif (imeth.eq.'tree') then
            known(:)=xneval(:)-(resi5(:)+resi4(:)+resi3(:))/dens2(:)
         endif

         if (diff.eq.1) then
            ! rank1 c-system: 4 coefficients
        
            do m=0,1
               f2(m)=effe(known,1,2,m)
            enddo
  
            do m=2,3
               mx=m-2
               f2(m)=effe(known,3,2,mx)
            enddo
  
            !---#[ getc2_S1:
            !----#[ original code:
!           c2(0) = f2(0)
!           c2(1) = (-f2(0) + f2(2))/(KB*MP12(2))
!           c2(3) = (KK(2)*(f2(1) - Fzc*f2(3)))/((-one + Fpc*Fzc)*MP12(2))
!           c2(5) = (-(Fpc*MP12(2)*c2(3)) - KK(2)*f2(3))/(KK(2)**2*MP12(2))
            !----#] original code:
            !----#[ HAGGIES:
            t1 = f2(0)
            c2(0) = t1
            c2(1) = ((f2(2)-t1)/(KB*MP12(2)))
            t1 = f2(3)
            c2(3) = ((f2(1)-t1*Fzc)/((Fpc*Fzc-one)*MP12(2))*KK(2))
            c2(5) = ((-(Fpc*MP12(2)*c2(3)+t1*KK(2)))/(KK(2)*KK(2)*MP12(2)))
            !----#] HAGGIES:
            !---#] getc2_S1:
         else
            ! traditional system
            do m=0,2
              f2(m)=effe(known,1,3,m)
            enddo
            do m=3,4
               mx=m-3
               f2(m)=effe(known,4,2,mx)
            enddo
            do m=5,6
               mx=m-5
               f2(m)=effe(known,6,2,mx)
            enddo
            f2(7)=known(8)
            f2(8)=known(9)
            f2(9)=known(10)
  
            !---#[ getc2_S2:
            !----#[ original code:
!           !--- coefficienti
!           Fpcpow2 = Fpc*Fpc
!           Fpcpow3 = Fpcpow2*Fpc
!           Fpcpow4 = Fpcpow3*Fpc
!           Fzcpow2 = Fzc*Fzc
!           Fzcpow3 = Fzcpow2*Fzc
!           Fzcpow4 = Fzcpow3*Fzc
!           Fzcpow5 = Fzcpow4*Fzc
!           Fzcpow6 = Fzcpow5*Fzc
!           tmp1 = one - Fpcpow2
!  c2(2) = (0.5_ki*(Fpcpow3*((one + Fzcpow4)*f2(0) + f2(1) - f2(3) + Fzc*&
!    &(Fzc*(f2(1) + (one + Fzcpow2)*f2(2) + Fzc*(-Fzc*f2(3) - f2(4))) - f2(4)))&
!    &+ Fpcpow4*(-f2(0) - Fzcpow2*f2(1) - Fzcpow4*f2(2) + f2(3) + Fzcpow3*f2(4))&
!    &+ Fmc**2*(tmp1*f2(0) + tmp1*Fzcpow2*f2(1) + tmp1*Fzcpow4*f2(2) + (-tmp1)&
!    &*f2(3) + (-tmp1)*Fzcpow3*f2(4)) + two*((-one - Fzcpow2 - Fzcpow4 + Fzcpow6)&
!    &*f2(0) + (-one - Fzcpow4)*f2(1) - f2(2) + Fzc*(-Fzc*f2(2) + Fzc*f2(3)&
!    &+ Fzcpow3*f2(3) + f2(4) + Fzcpow4*f2(4))) + f2(5) + f2(6) + Fmc*((one&
!    &+ Fzcpow4 + Fzcpow3*two + Fpcpow2*(-one - Fzcpow4 - Fzcpow3*two))*f2(0)&
!    &+ (one + Fzcpow2 + Fzcpow5*two + Fpcpow2*(-one - Fzcpow2 - Fzcpow5*two))*f2(1)&
!    &+ tmp1*Fzcpow2*f2(2) + (-tmp1)*f2(3) + Fzcpow4*(tmp1*f2(2) + (-tmp1)*f2(3))&
!    &- two*f2(4) + Fpcpow2*two*f2(4) - tmp1*Fzcpow3*(two*f2(3) + f2(4)) + Fzc&
!    &*(tmp1*two*f2(2) - tmp1*f2(4)) - f2(5) - Fpc*f2(6) + Fzcpow6*(f2(5) + Fpc&
!    &*f2(6) - f2(7)) + f2(7)) + Fpc*((-one - Fzcpow4)*f2(0) - f2(1) + f2(3)&
!    &+ f2(5) - f2(7) + Fzc*(f2(4) + Fzc*(-f2(1) - (one + Fzcpow2)*f2(2) + Fzc&
!    &*(Fzc*f2(3) + f2(4) + Fzcpow3*(-f2(5) + f2(7)))))) + Fzcpow6*(-f2(5) - f2(6)&
!    &- f2(8)) + f2(8) + Fpcpow2*((three + Fzcpow2*(one + Fzcpow2 - Fzcpow4)*two)&
!    &*f2(0) + two*f2(1) + two*f2(2) - f2(3) + Fzcpow4*(two*f2(1) + f2(2) - two&
!    &*f2(3)) + Fzcpow2*(f2(1) + two*f2(2) - two*f2(3)) - Fzcpow3*f2(4) - Fzc*two&
!    &*f2(4) - Fzcpow5*two*f2(4) - f2(5) - f2(8) + Fzcpow6*(f2(5) + f2(8)))))/&
!    &((-tmp1)*(-one + Fzcpow6)*MP12(2)**2)
!  c2(7) = (KK(2)*(f2(2) + Fpcpow2*(-(Fzcpow4*f2(1)) - f2(2) + Fzcpow2*(-f2(0)&
!    &+ f2(3)) + Fzcpow5*f2(4)) + Fpcpow3*(-((one + Fzcpow4)*f2(0))&
!    &- (one + Fzcpow2)*(f2(1) + Fzcpow2*f2(2)) + (one + Fzcpow4)*f2(3) &
!    &+ (Fzc + Fzcpow3)*f2(4)) - f2(6) + Fzcpow2*(f2(0) - f2(3) +&
!    &Fzcpow2*(f2(1) + Fzc*(-f2(4) + Fzc*f2(6)))) +&
!    &Fpc*((one + Fzcpow4)*f2(0) + f2(1) - f2(3) - f2(5) +&
!    &Fzc*(Fzc*(f2(1) + f2(2)) + Fzcpow3*(f2(2) - f2(3)) -&
!    &f2(4) - Fzcpow2*f2(4) + Fzcpow5*(f2(5) - f2(7))) +&
!    &f2(7))))/((-tmp1)*(-one + Fzcpow6)*MP12(2)**2)
!  c2(8) = ((-tmp1)*(one + Fzcpow3 + Fzcpow4)*f2(0) +&
!    &(-tmp1)*(one + Fzcpow2 + Fzcpow5)*f2(1) + (-tmp1)*Fzcpow2*f2(2) +&
!    &(-tmp1)*Fzcpow4*(f2(2) - f2(3)) + f2(3) + (-tmp1)*Fzc*(f2(2) - f2(4))&
!    &+ f2(4) - (-tmp1)*Fzcpow3*(f2(3) + f2(4)) + f2(5) +&
!    &Fpc*(-(Fpc*(f2(3) + f2(4))) + f2(6)) - Fzcpow6*(f2(5) + Fpc*f2(6) - f2(7))&
!    &- f2(7))/ ((-tmp1)*(-one + Fzcpow6)*KK(2)*MP12(2)**2)
!  c2(4) = -((KK(2)**2*(f2(1) + Fzcpow2*f2(2) + Fzcpow4*(f2(0) - f2(3)) -&
!    &Fzc*f2(4)))/((-one + Fzcpow6)*MP12(2)**2))
!  c2(5) = (Fzcpow3*f2(0) + Fzcpow5*f2(1) + Fzc*f2(2) - Fzcpow3*f2(3) -&
!    &f2(4))/(KK(2)*MP12(2) - Fzcpow6*KK(2)*MP12(2))
!  c2(3) = (KK(2)*(f2(2) + Fzcpow2*(f2(0) - f2(3) + Fzcpow2*(f2(1)&
!  - Fzc*f2(4)))))/((-one + Fzcpow6)*MP12(2))
!  c2(6) = (f2(0) - f2(3) + Fzcpow2*(f2(1) + Fzc*(Fzc*f2(2) - f2(4))))/&
!    &((-one + Fzcpow6)*KK(2)**2*MP12(2)**2)
!  c2(0) = f2(0)
!  c2(9) = (KK(2)*MP12(2)*c2(3) - MP12(2)**2*c2(4)&
!    &+ F1zc*KK(2)**3*MP12(2)*c2(5) - F1zc**2*KK(2)**4*MP12(2)**2*c2(6) &
!    &+ KK(2)**2*(-c2(0) + f2(9)))/(KK(2)**2*mu2g(2))
!  c2(1) = -(c2(3)/KK(2)) - Fmc*KK(2)*c2(5) + (MP12(2)*(c2(4) + KK(2)*(c2(7)&
!    &+ KK(2)*(c2(2) + Fmc*KK(2)*(Fmc*KK(2)*c2(6) + c2(8))))))/KK(2)**2&
!    &+ (c2(0) - f2(8))/MP12(2)
            !----#] original code:
            !----#[ HAGGIES:
            t1 = (1.0_ki)+Fpc
            t2 = Fzc*Fzc
            t3 = t2*Fzc
            t4 = t3-(1.0_ki)
            t5 = (1.0_ki)+t3
            t6 = f2(2)
            t7 = t2*t2
            t8 = -(t7+(1.0_ki))
            t9 = f2(1)
            t10 = t3*t3
            t11 = f2(0)
            t12 = f2(4)
            t13 = t6*Fzc
            t14 = f2(3)
            t15 = t14*Fzc
            t16 = t15*t2
            t17 = f2(5)
            t18 = f2(6)
            t19 = f2(8)
            t20 = t9-t14
            t21 = (1.0_ki)+t7
            t22 = t11*t21
            t23 = (1.0_ki)+t2
            t24 = t23*t6
            t25 = Fpc*Fpc
            t26 = t14-t11
            t27 = t12*Fzc
            t28 = t2*t27
            t29 = ((1.0_ki)-Fpc)*t1
            t30 = t14*t29
            t31 = t29*t6
            t32 = t12*t29
            t33 = f2(7)
            t34 = t6+t9
            t35 = t27*t7
            t36 = t6-t14
            t37 = t17-t33
            t38 = t18*Fpc
            t39 = (t37+t38)*t10
            t40 = t7*Fzc
            t41 = (2.0_ki)*t40
            t42 = (2.0_ki)*t3
            t43 = t2*t31
            t44 = MP12(2)*MP12(2)
            t45 = t29*t4*t44*t5
            c2(2) = (-1.0_ki)*((2.0_ki)*((t10-(t7+t2+(1.0_ki)))*t11+(t12+t15+t&
            &16+t12*t7-t13)*Fzc+t8*t9-t6)+t17+t18+t19+(-(t19+t18+t17))*t10+(t26&
            &+t28-(t6*t7+t2*t9))*t25*t25+(t11*t29+t31*t7+t2*t29*t9-(t2*t32*Fzc+&
            &t30))*Fmc*Fmc+(t14+t17+(t12+((t12+t15+(t33-t17)*t2*Fzc)*Fzc-(t9+t24))&
            &*Fzc)*Fzc+t11*t8-(t9+t33))*Fpc+((2.0_ki)*(t34-(t35+t27))+((3.0_ki)+(&
            &2.0_ki)*((1.0_ki)+t2-t7)*t2)*t11+(t17+t19)*t10+((2.0_ki)*t20+t6)*&
            &t7+((2.0_ki)*t36+t9)*t2-(t28+t19+t17+t14))*t25+((2.0_ki)*(t12*t25&
            &-t12)+t33+t39+t43+(t31-t30)*t7+((2.0_ki)*t31-t32)*Fzc+((1.0_ki)+t2&
            &+t41+(-(t41+t2+(1.0_ki)))*t25)*t9+((1.0_ki)+t42+t7+(-(t7+t42+(1.0&
            &_ki)))*t25)*t11-((t12+(2.0_ki)*t14)*t2*t29*Fzc+t38+t30+t17))*Fmc+(t&
            &20+t22+((t24+t9+(-(t15+t12))*Fzc)*Fzc-t12)*Fzc)*t25*Fpc)/t45*0.5_ki
            t8 = t11-t14
            t10 = t2*t6
            c2(7) = (-1.0_ki)*(t6+(t8+(t9+(t18*Fzc-t12)*Fzc)*t2)*t2+(t35+t2*t26-&
            &(t7*t9+t6))*t25+(t20+t22+t33+(t34*Fzc+t2*t36*Fzc+t37*t7*Fzc-(t12*t2+&
            &t12))*Fzc-t17)*Fpc+((t3+Fzc)*t12+t14*t21-((t10+t9)*t23+t22))*t25*Fpc-&
            &t18)/t45*KK(2)
            t14 = t12+t14
            c2(8) = ((t14+t17+(t18-t14*Fpc)*Fpc+t14*t2*t29*Fzc-(t29*t36*t7+((1.0_&
            &ki)+t3+t7)*t11*t29+((1.0_ki)+t2+t40)*t29*t9+(t6-t12)*t29*Fzc+t43+t&
            &39+t33))/((Fpc-(1.0_ki))*t1*t4*t44*t5*KK(2)))
            t1 = KK(2)*KK(2)
            t3 = t4*t5
            c2(4) = (-1.0_ki)*(t10+t9+t7*t8-t27)/(t3*t44)*t1
            c2(5) = (-1.0_ki)*(t13+t11*t2*Fzc+t7*t9*Fzc-(t16+t12))/(t3*KK(2)*MP1&
            &2(2))
            c2(3) = ((t6+(t8+(t9-t27)*t2)*t2)/(t3*MP12(2))*KK(2))
            t4 = KK(2)*MP12(2)
            c2(6) = ((t8+(t9+(t13-t12)*Fzc)*t2)/(t3*t4*t4))
            c2(0) = t11
            t2 = t1*F1zc*MP12(2)
            c2(9) = (((f2(9)-c2(0))*t1+t4*c2(3)+t1*F1zc*KK(2)*MP12(2)*c2(5)-(t4&
            &4*c2(4)+t2*t2*c2(6)))/(t1*mu2g(2)))
            t2 = Fmc*KK(2)
            c2(1) = ((c2(0)-t19)/MP12(2)+(c2(4)+(c2(7)+(c2(2)+(c2(8)+t2*c2(6))&
            &*Fmc*KK(2))*KK(2))*KK(2))/t1*MP12(2)-(c2(3)/KK(2)+t2*c2(5)))
            !----#] HAGGIES:
            !---#] getc2_S2:
         endif
         !---#] Decomposizione standard:
      endif

      if (diff.ge.1) then
         c2(2)=czip
         c2(4)=czip
         c2(6)=czip
         c2(7)=czip
         c2(8)=czip
         c2(9)=czip
         if (diff.ge.2) then
            c2(1)=czip
            c2(3)=czip
            c2(5)=czip
         endif             
      endif
   end subroutine getc2_cm

end module mgetc2

