module mgetc3
   use precision, only: ki
   use constants
   use options
   use mfunctions
   use mrestore
   implicit none

   private

   interface getc3
      module procedure getc3_rm
      module procedure getc3_cm
   end interface getc3


   public :: getc3

contains

   subroutine getc3_cm(numeval,nleg,rank,c3,cut3,q3,qt,Vi,msq)
      use mglobal, only: C0c,C1c,mu2g,MP12,KK,mu2t,resit,denst,mu2test
      implicit none
      integer, intent(in) :: nleg, rank, cut3
      complex(ki), dimension(0:9), intent(out) :: c3
      complex(ki), dimension(10,4), intent(in) :: q3
      complex(ki), dimension(4), intent(in) :: qt
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(10) :: mu2vec

      integer :: i,m,n,i1,i2,i3,i4,i5,dicut5,dicut4,mx,j1,j2,j3,diff
      integer :: nsol, acc
      complex(ki), dimension(10) :: dens3,dens4,dens5,xneval
      complex(ki), dimension(10) :: resi5,resi4,known
      complex(ki), dimension(0:9) :: f3
      complex(ki) :: dens4t,dens5t
      logical evalres 

      !!! TR: I have used haggies to rewrite the systems.
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
      !---#] HAGGIES:

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

      mu2test(3) = mu2t(3)

      j3=cut3/100
      acc = j3*100
      j2=(cut3-acc)/10
      j1=cut3-acc-j2*10

      resi5(:)=czip
      resi4(:)=czip
      known(:)=czip
      xneval(:)=czip
      dens3(:)=cone

      !---  for lnntest
      resit(3)=czip
      denst(3)=cone

      !---  for simplified sampling
      diff = nleg-rank

      if_diff: if (diff.ge.3) then
         ! simplified sampling -- only c3(0)

         select case(nleg)
         case(5)
            resi5(1)=res5(1,czip)
            resit(3)=res5(1,mu2t(3))
            goto 11
         case(4)
            resi4(1)=Res4(1,q3(1,:),czip)
            resit(3)   =res4(1,qt,mu2t(3))
            goto 21
         case(3)
            goto 26
         case default
            dicut5=1
            loop_i5: do i5=4,nleg-1
               loop_i4: do i4=3,i5-1
                  loop_i3: do i3=2,i4-1
                     loop_i2: do i2=1,i3-1
                        loop_i1: do i1=0,i2-1

                           dens5(1)=cone
                           dens5t=cone
                        
                           evalres=.false.
      
                           loop_10: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2) &
                                   & .and.(i.ne.i3)&
                                   & .and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i.eq.j1).or.(i.eq.j2).or.(i.eq.j3)) then
                                    dens5(1)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_10
                                 else
                                    dens5(1)=dens5(1)&
                                      &*denevalmu2(nleg,i,q3(1,:),Vi,msq,czip)
                                    dens5t=dens5t&
                                      &*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_10

                           if (evalres) then
                              resi5(1)=resi5(1)+dens5(1)*res5(dicut5,czip)
                              resit(3)=resit(3)+dens5t*Res5(dicut5,mu2t(3))
                           endif

                         dicut5=dicut5+1
                      enddo loop_i1
                   enddo loop_i2
                enddo loop_i3
             enddo loop_i4
          enddo loop_i5

       end select

 11    continue

       dicut4=1
       loop_21: do i4=3,nleg-1
          loop_21_i3: do i3=2,i4-1
             loop_21_i2: do i2=1,i3-1
                loop_21_i1: do i1=0,i2-1

                   dens4(1)=cone
                   dens4t=cone
                   evalres=.false.
      
                   loop_20: do i=0,nleg-1
                      if ((i.ne.i1).and.(i.ne.i2) &
                         &  .and.(i.ne.i3).and.(i.ne.i4)) then
                         if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
                            dens4(1)=dens4(1)&
                                &*denevalmu2(nleg,i,q3(1,:),Vi,msq,czip)
                            dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
                            evalres=.true.
                         else
                            dens4(1)=czip
                            dens4t=czip
                            evalres=.false.
                            exit loop_20
                         endif
                      endif
                   enddo loop_20

                   if (evalres) then
                      resi4(1)=resi4(1)+dens4(1)*Res4(dicut4,q3(1,:),czip)
                      resit(3)=resit(3)+dens4t*Res4(dicut4,qt,mu2t(3))
                   endif

                   dicut4=dicut4+1
                enddo loop_21_i1
              enddo loop_21_i2
           enddo loop_21_i3
        enddo loop_21

 21     continue

        loop_26: do i=0,nleg-1
           if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
             dens3(1)=dens3(1)*denevalmu2(nleg,i,q3(1,:),Vi,msq,czip)
             denst(3)=denst(3)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
           endif
        enddo loop_26

 26     continue

        xneval(1)=numeval(cut3,q3(1,:),czip)
        if     (imeth.eq.'diag') then
           known(1)=(xneval(1)-resi5(1)-resi4(1))/dens3(1)
        elseif (imeth.eq.'tree') then
           known(1)=xneval(1)-(resi5(1)+resi4(1))/dens3(1)
        endif
 
        c3(0)=known(1)
        do m=1,9
           c3(m)=czip
        enddo

     else
         if (abs(C0c-1.0_ki) .lt. C0_thrs) then
            !---#[ New Sampling:
            ! The new sampling is the one that is safe around C0=1
            ! but not around C0=0

            ! traditional system
            mu2vec=(/czip,czip,czip,czip,czip,czip,czip,mu2g(3),mu2g(3),mu2g(3)/)
            nsol=10
            !---#] New Sampling:
         else
            !---#[ Old Sampling:
            ! The old sampling is the one that is safe around C0=0
            ! but not around C0=1
            if (diff.eq.2) then
               ! rank1 c-system: 3 coefficients
               mu2vec=(/czip,czip,czip,czip,czip,czip,czip,czip,czip,czip/)
               nsol=3
            elseif (diff.eq.1) then
               ! rank1 c-system: 3 coefficients
               mu2vec=(/czip,czip,czip,czip,czip,mu2g(3),czip,czip,czip,czip/)
               nsol=6
            else
               ! traditional system
               mu2vec=(/czip,czip,czip,czip,czip,czip,czip,mu2g(3),mu2g(3),mu2g(3)/)
               nsol=10
            endif
            !---#] Old Sampling:
         endif
      

        if_nleg: if     (nleg.eq.5) then
           do n=1,nsol
              resi5(n)=Res5(1,mu2vec(n))
           enddo
           resit(3)=Res5(1,mu2t(3))         
           goto 111
        elseif (nleg.eq.4) then
           do n=1,nsol
              resi4(n)=Res4(1,q3(n,:),mu2vec(n))
           enddo
           resit(3)=Res4(1,qt,mu2t(3))         
           goto 121
        elseif (nleg.eq.3) then
           goto 126
        else
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
                                  & .and.(i.ne.i3).and.&
                                  &(i.ne.i4).and.(i.ne.i5)) then
                                if ((i.eq.j1).or.(i.eq.j2).or.(i.eq.j3)) then
                                   dens5(:)=czip
                                   dens5t=czip
                                   evalres=.false.
                                   exit loop_110
                                else
                                   do n=1,nsol
                                      dens5(n)=dens5(n)*denevalmu2(nleg,i,&
                                            &q3(n,:),Vi,msq,mu2vec(n))
                                   enddo
                                   dens5t=dens5t&
                                     &*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
                                   evalres=.true.
                                endif
                             endif
                          enddo loop_110

                          if (evalres) then
                             do n=1,nsol
                                resi5(n)=resi5(n)&
                                   &+dens5(n)*res5(dicut5,mu2vec(n))
                             enddo
                             resit(3)=resit(3)+dens5t*res5(dicut5,mu2t(3))
                          endif

                          dicut5=dicut5+1
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        endif if_nleg

 111    continue

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
                           if ((i.eq.j1).or.(i.eq.j2).or.(i.eq.j3)) then
                              dens4(:)=czip
                              dens4t=czip
  
                              evalres=.false.
                              exit loop_120
                           else
                              do n=1,nsol
                                 dens4(n)=dens4(n)&
                                   &*denevalmu2(nleg,i,q3(n,:),Vi,msq,mu2vec(n))
                              enddo
                              dens4t=dens4t&
                                  &*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_120

                     if (evalres) then
                        do n=1,nsol
                           resi4(n)=resi4(n)&
                                & +dens4(n)*Res4(dicut4,q3(n,:),mu2vec(n))
                        enddo
                        resit(3)=resit(3)+dens4t*Res4(dicut4,qt,mu2t(3))
                     endif

                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo

 121     continue

         do i=0,nleg-1
            if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
               do n=1,nsol
                  dens3(n)=dens3(n)*denevalmu2(nleg,i,q3(n,:),Vi,msq,mu2vec(n))
               enddo
               denst(3)=denst(3)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
            endif
         enddo

 126     continue

         do n=1,nsol
            xneval(n)=numeval(cut3,q3(n,:),mu2vec(n))
         enddo

         if     (imeth.eq.'diag') then
            known(:)=(xneval(:)-resi5(:)-resi4(:))/dens3(:)
         elseif (imeth.eq.'tree') then
            known(:)=xneval(:)-(resi5(:)+resi4(:))/dens3(:)
         endif

         if (abs(C0c-1.0_ki) .lt. C0_thrs) then
            !---#[ New Sampling:
            ! The new sampling is the one that is safe around C0=1
            ! but not around C0=0

            ! traditional system
        
            do m=0,6
               f3(m)=effe(known,1,7,m)
            enddo

            do m=7,9
               mx=m-7
               f3(m)=effe(known,8,3,mx)
            enddo

            !---#[ getc3_S3:
            c3(0) =   f3(0)
            c3(1) = -(f3(6)/(C0c*MP12(3)))
            c3(2) =   f3(5)/(C0c**2*MP12(3)**2)
            c3(3) = -(f3(4)/(C0c**3*MP12(3)**3))
            c3(4) = -(f3(1)/MP12(3))
            c3(5) =   f3(2)/MP12(3)**2
            c3(6) = -(f3(3)/MP12(3)**3)
            c3(7) = (-f3(0) + MP12(3)**3*(-((C1c**3*f3(3))/MP12(3)**3) -&
            &f3(4)/(C0c**3*MP12(3)**3)) + f3(7))/mu2g(3)
            c3(8) = (C1c**2*f3(2) + f3(6)/C0c - f3(9))/(MP12(3)*mu2g(3))
            c3(9) = (C1c*f3(1) + f3(5)/C0c**2 - f3(8))/(C1c*MP12(3)*mu2g(3))
            !---#] getc3_S3:  

            !---#] New Sampling:
         else
            !---#[ Old Sampling:
            ! The old sampling is the one that is safe around C0=0
            ! but not around C0=1
            select case(diff)
            case(2)
               ! rank1 c-system: 3 coefficients
               f3(0:2)=known(1:3)
               !---#[ getc3_S1:  
               !----#[ original code:
      !      c3(0) = (f3(0) + f3(1))/two
      !      c3(1) = -(KK(3)*((one + C0c)*f3(0) + f3(1) - C0c*f3(1) - two*f3(2)))/ &
      !        (two*(-one + C0c**2)*MP12(3))
      !      c3(4) = &
      !       -((KK(3)*c3(0) + C0c*MP12(3)*c3(1) - KK(3)*f3(1))/(KK(3)**2*MP12(3)))
               !----#] original code:
               !----#[ HAGGIES:
               t1 = f3(0)
               t2 = f3(1)
               c3(0) = ((t1+t2)*0.5_ki)
               t3 = KK(3)
               t4 = MP12(3)
               c3(1) = (-1.0_ki)*(t2+((1.0_ki)+C0c)*t1-(t2*C0c+(2.0_ki)*f3(2)))/((C&
               &0c*C0c-(1.0_ki))*t4)*t3*0.5_ki
               c3(4) = (-1.0_ki)*(t3*c3(0)+t4*C0c*c3(1)-t2*t3)/(t3*t3*t4)
               !----#] HAGGIES:
               !---#] getc3_S1:  
            case(1)
               ! rank1 c-system: 3 coefficients
               do m=0,2
                  f3(m)=effe(known,1,3,m)
               enddo
               
               do m=3,4
                  mx=m-3
                  f3(m)=effe(known,4,2,mx)
               enddo
               
               do m=5,5
                  mx=m-5
                  f3(m)=effe(known,6,1,mx)
               enddo
               !---#[ getc3_S2:
               !----#[ original code:
      !!$      c3(0) = f3(0)
      !!$      c3(2) = &
      !!$     &        (KK(3)**2*(f3(0) - f3(3) + &
      !!$     &      C0c**2*(C0c**2*f3(1) + f3(2) - C0c*f3(4))))/&
      !!$     &  ((-one + C0c**6)*MP12(3)**2)
      !!$      c3(4) = (C0c**2*MP12(3)*c3(2))/KK(3)**3 - f3(1)/(KK(3)*MP12(3))
      !!$      c3(1) = -((KK(3)*(C0c*KK(3)*MP12(3)*c3(4) + f3(4)))/MP12(3))
      !!$      c3(5) = (C0c*MP12(3)*c3(1) + KK(3)*f3(2))/(KK(3)**3*MP12(3)**2)
      !!$      c3(7) = &
      !!$     &  (MP12(3)*c3(1) - MP12(3)**2*c3(2) + C1c*MP12(3)*c3(4)-& 
      !!$     &    C1c**2*MP12(3)**2*c3(5) + (-c3(0) + f3(5)))/mu2g(3)
      
      !        c3(0) = f3(0)
      !     c3(2) = &
      !    &        (f3(0) - f3(3) + &
      !    &      C0c**2*(C0c**2*f3(1) + f3(2) - C0c*f3(4)))/&
      !    &  ((-one + C0c**6)*MP12(3)**2)
      !     c3(4) = C0c**2*MP12(3)*c3(2) - f3(1)/MP12(3)
      !     c3(1) = -(C0c*MP12(3)*c3(4) + f3(4))/MP12(3)
      !     c3(5) = (C0c*MP12(3)*c3(1) + f3(2))/MP12(3)**2
      !     c3(7) = &
      !    &  (MP12(3)*c3(1) - MP12(3)**2*c3(2) + C1c*MP12(3)*c3(4)-& 
      !    &    C1c**2*MP12(3)**2*c3(5) + (-c3(0) + f3(5)))/mu2g(3)
                  !----#] original code:
                  !----#[ HAGGIES:
               t1 = f3(0)
               c3(0) = t1
               t2 = C0c*C0c
               t3 = t2*C0c
   
               t4 = f3(2)
               t5 = f3(1)
               t6 = f3(4)
               t7 = MP12(3)
               t8 = t7*t7
               c3(2) = ((t1+(t4+t2*t5-t6*C0c)*t2-f3(3))/((t3*t3-one)*t8))
               c3(4) = (t2*t7*c3(2)-t5/t7)
               t1 = t7*C0c
               c3(1) = (-1.0_ki)*(t6+t1*c3(4))/t7
               c3(5) = ((t4+t1*c3(1))/t8)
               t1 = t7*C1c
               c3(7) = ((f3(5)+t1*c3(4)+t7*c3(1)-(t8*c3(2)+t1*t1*c3(5)+c3(0)))/(m&
               &u2g(3)))
               !----#] HAGGIES:
               !---#] getc3_S2:  
            case default 
               ! traditional system
              
               do m=0,3
                  f3(m)=effe(known,1,4,m)
               enddo
      
               do m=4,6
                  mx=m-4
                  f3(m)=effe(known,5,3,mx)
               enddo
   
               do m=7,9
                  mx=m-7
                  f3(m)=effe(known,8,3,mx)
               enddo
      
               !---#[ getc3_S3:
               !----#[ original code:
      !     c3(0) =  &
      !    &f3(0)
      !     c3(1) =  &
      !    &        -((KK(3)*(C0c**5*f3(1) + C0c**2*f3(2) + C0c**11*f3(3) +  &
      !    &        C0c**8*(f3(0) - f3(4)) - C0c**4*f3(5) - f3(6)))/ &
      !    &    ((-one + C0c**12)*MP12(3)))
      !     c3(2) =  &
      !    &        (KK(3)**2*(C0c*f3(1) + C0c**10*f3(2) + C0c**7*f3(3) +  &
      !    &      C0c**4*(f3(0) - f3(4)) - f3(5) - C0c**8*f3(6)))/ &
      !    &  ((-one + C0c**12)*MP12(3)**2)
      !     c3(3) =  &
      !    &        -((KK(3)**3*(f3(0) + C0c**9*f3(1) + C0c**6*f3(2) +  &
      !    &        C0c**3*f3(3) - f3(4) - C0c**8*f3(5) - C0c**4*f3(6)))/ &
      !    &    ((-one + C0c**12)*MP12(3)**3))
      !     c3(4) =  &
      !    &        (f3(1) + C0c**9*f3(2) + C0c**6*f3(3) +  &
      !    &    C0c**3*(f3(0) - f3(4)) - C0c**11*f3(5) - C0c**7*f3(6))/ &
      !    &  ((-one + C0c**12)*KK(3)*MP12(3))
      !     c3(5) =  &
      !    &        (-(C0c**3*f3(1)) - f3(2) - C0c**9*f3(3) +  &
      !    &    C0c**6*(-f3(0) + f3(4)) + C0c**2*f3(5) + C0c**10*f3(6))/ &
      !    &  ((-one + C0c**12)*KK(3)**2*MP12(3)**2)
      !     c3(6) =  &
      !    &        (C0c**6*f3(1) + C0c**3*f3(2) + f3(3) +  &
      !    &    C0c**9*(f3(0) - f3(4)) - C0c**5*f3(5) - C0c*f3(6))/ &
      !    &  ((-one + C0c**12)*KK(3)**3*MP12(3)**3)
      !     c3(7) =  &
      !    &        ((C1c**3*KK(3)**3*(C0c**6*f3(1) + C0c**3*f3(2) + f3(3) +  &
      !    &         C0c**9*(f3(0) - f3(4)) - C0c**5*f3(5) - C0c*f3(6)))/ &
      !    &     (-one + C0c**12) -  &
      !    &    (KK(3)**3*(f3(0) + C0c**9*f3(1) + C0c**6*f3(2) +  &
      !    &         C0c**3*f3(3) - f3(4) - C0c**8*f3(5) - C0c**4*f3(6)))/ &
      !    &     (-one + C0c**12) + KK(3)**3*(-f3(0) + f3(7)))/(KK(3)**3*mu2g(3))
      !     c3(8) =  &
      !    &        -((-((KK(3)*(C0c**5*f3(1) + C0c**2*f3(2) + C0c**11*f3(3) +  &
      !    &             C0c**8*(f3(0) - f3(4)) - C0c**4*f3(5) - f3(6)))/ &
      !    &         (-one + C0c**12)) -  &
      !    &      (C1c**2*KK(3)*(-(C0c**3*f3(1)) - f3(2) - C0c**9*f3(3) +  &
      !    &           C0c**6*(-f3(0) + f3(4)) + C0c**2*f3(5) +  &
      !    &           C0c**10*f3(6)))/(-one + C0c**12) + KK(3)*f3(9))/ &
      !    &    (MP12(3)*mu2g(3)))
      !     c3(9) =  &
      !    &        (-((C1c*KK(3)**2*(f3(1) + C0c**9*f3(2) + C0c**6*f3(3) +  &
      !    &           C0c**3*(f3(0) - f3(4)) - C0c**11*f3(5) -  &
      !    &           C0c**7*f3(6)))/(-one + C0c**12)) +  &
      !    &    (KK(3)**2*(C0c*f3(1) + C0c**10*f3(2) + C0c**7*f3(3) +  &
      !    &         C0c**4*(f3(0) - f3(4)) - f3(5) - C0c**8*f3(6)))/ &
      !    &     (-one + C0c**12) - KK(3)**2*f3(8))/(C1c*KK(3)**3*MP12(3)*mu2g(3))
               !----#] original code:
               !----#[ HAGGIES:
               t1 = f3(0)
               c3(0) = t1
               t2 = C0c*C0c
               t3 = t2*C0c
               t4 = t3*t3
               t5 = f3(6)
               t6 = f3(4)
               t7 = t1-t6
               t8 = t2*t2
               t9 = t8*t8
               t10 = f3(1)
               t11 = t10*C0c
               t12 = f3(2)
               t13 = f3(3)
               t14 = t13*C0c
               t15 = t8*C0c
               t15 = t15*t15
               t16 = f3(5)
               t17 = KK(3)
               t18 = MP12(3)
               t19 = (t11*t8+t12*t2+t14*t15+t7*t9-(t16*t8+t5))*t17
               t3 = (t3-(1.0_ki))*((1.0_ki)+t3)*((1.0_ki)+t4)
               t20 = t18*t3
               c3(1) = (-1.0_ki)*t19/t20
               t21 = t17*t17
               t22 = (t11+t12*t15+t14*t4+t7*t8-(t5*t9+t16))*t21
               t23 = t18*t18
               c3(2) = (t22/(t23*t3))
               t24 = (t7+t11*t9+t12*t4+t14*t2-(t5*t8+t16*t9))*t17*t21
               c3(3) = (-1.0_ki)*t24/(t20*t23)
               t7 = t7*C0c
               t20 = t12*C0c
               t23 = t16*C0c
               t25 = t5*C0c
               t26 = t10+t13*t4+t2*t7+t20*t9-(t25*t4+t15*t23)
               t27 = t17*t18*t3
               c3(4) = (t26/t27)
               t5 = (t6-t1)*t4+t15*t5+t16*t2-(t14*t9+t11*t2+t12)
               t6 = t17*t18
               t11 = t6*t6
               c3(5) = (t5/(t11*t3))
               t2 = t13+t10*t4+t2*t20+t7*t9-(t23*t8+t25)
               c3(6) = (t2/(t11*t27))
               t4 = t17*C1c
               t7 = mu2g(3)
               c3(7) = (((f3(7)-t1)*t17*t21+t17*t2/t3*t4*t4*C1c-t24/t3)/(t17*t21*t&
               &7))
               c3(8) = (-1.0_ki)*(f3(9)*t17-(t17/t3*t5*C1c*C1c+t19/t3))/(t18*t7)
               c3(9) = ((t22/t3-(t21*t26/t3*C1c+f3(8)*t21))/(t21*t6*t7*C1c))
               !----#] HAGGIES:
               !---#] getc3_S3:  
            end select
            !---#] Old Sampling:
         endif
      end if if_diff

      if (diff.ge.1) then
         c3(3)=czip
         c3(6)=czip
         c3(8)=czip
         c3(9)=czip 
         if (diff.ge.2) then
            c3(2)=czip
            c3(5)=czip
            c3(7)=czip
            if (diff.ge.3) then
               c3(1)=czip
               c3(4)=czip
            endif 
         endif 
      endif
   end subroutine getc3_cm

   subroutine getc3_rm(numeval,nleg,rank,c3,cut3,q3,qt,Vi,msq)
      use mglobal, only: C0,C1,mu2g,MP12,KK,mu2t,resit,denst,mu2test
      use options, only: C0_thrs
      implicit none
      integer, intent(in) :: nleg, rank, cut3
      complex(ki), dimension(0:9), intent(out) :: c3
      complex(ki), dimension(10,4), intent(in) :: q3
      complex(ki), dimension(4), intent(in) :: qt
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(10) :: mu2vec

      integer :: i,m,n,i1,i2,i3,i4,i5,dicut5,dicut4,mx,j1,j2,j3,diff
      integer :: nsol, acc
      complex(ki), dimension(10) :: dens3,dens4,dens5,xneval
      complex(ki), dimension(10) :: resi5,resi4,known
      complex(ki), dimension(0:9) :: f3
      complex(ki) :: dens4t,dens5t
      logical evalres 

      !!! TR: I have used haggies to rewrite the systems.
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
      !---#] HAGGIES:

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

      mu2test(3) = mu2t(3)

      j3=cut3/100
      acc = j3*100
      j2=(cut3-acc)/10
      j1=cut3-acc-j2*10

      resi5(:)=czip
      resi4(:)=czip
      known(:)=czip
      xneval(:)=czip
      dens3(:)=cone

      !---  for lnntest
      resit(3)=czip
      denst(3)=cone

      !---  for simplified sampling
      diff = nleg-rank


      if_diff: if (diff.ge.3) then
         ! simplified sampling -- only c3(0)

         select case(nleg)
         case(5)
            resi5(1)=res5(1,czip)
            resit(3)=res5(1,mu2t(3))
            goto 11
         case(4)
            resi4(1)=Res4(1,q3(1,:),czip)
            resit(3)   =res4(1,qt,mu2t(3))
            goto 21
         case(3)
            goto 26
         case default
            dicut5=1
            loop_i5: do i5=4,nleg-1
               loop_i4: do i4=3,i5-1
                  loop_i3: do i3=2,i4-1
                     loop_i2: do i2=1,i3-1
                        loop_i1: do i1=0,i2-1

                           dens5(1)=cone
                           dens5t=cone
                        
                           evalres=.false.
      
                           loop_10: do i=0,nleg-1
                              if ((i.ne.i1).and.(i.ne.i2) &
                                   & .and.(i.ne.i3)&
                                   & .and.(i.ne.i4).and.(i.ne.i5)) then
                                 if ((i.eq.j1).or.(i.eq.j2).or.(i.eq.j3)) then
                                    dens5(1)=czip
                                    dens5t=czip
                                    evalres=.false.
                                    exit loop_10
                                 else
                                    dens5(1)=dens5(1)&
                                      &*denevalmu2(nleg,i,q3(1,:),Vi,msq,czip)
                                    dens5t=dens5t&
                                      &*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
                                    evalres=.true.
                                 endif
                              endif
                           enddo loop_10

                           if (evalres) then
                              resi5(1)=resi5(1)+dens5(1)*res5(dicut5,czip)
                              resit(3)=resit(3)+dens5t*Res5(dicut5,mu2t(3))
                           endif

                         dicut5=dicut5+1
                      enddo loop_i1
                   enddo loop_i2
                enddo loop_i3
             enddo loop_i4
          enddo loop_i5

       end select

 11    continue

       dicut4=1
       loop_21: do i4=3,nleg-1
          loop_21_i3: do i3=2,i4-1
             loop_21_i2: do i2=1,i3-1
                loop_21_i1: do i1=0,i2-1

                   dens4(1)=cone
                   dens4t=cone
                   evalres=.false.
      
                   loop_20: do i=0,nleg-1
                      if ((i.ne.i1).and.(i.ne.i2) &
                         &  .and.(i.ne.i3).and.(i.ne.i4)) then
                         if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
                            dens4(1)=dens4(1)&
                                &*denevalmu2(nleg,i,q3(1,:),Vi,msq,czip)
                            dens4t=dens4t*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
                            evalres=.true.
                         else
                            dens4(1)=czip
                            dens4t=czip
                            evalres=.false.
                            exit loop_20
                         endif
                      endif
                   enddo loop_20

                   if (evalres) then
                      resi4(1)=resi4(1)+dens4(1)*Res4(dicut4,q3(1,:),czip)
                      resit(3)=resit(3)+dens4t*Res4(dicut4,qt,mu2t(3))
                   endif

                   dicut4=dicut4+1
                enddo loop_21_i1
              enddo loop_21_i2
           enddo loop_21_i3
        enddo loop_21

 21     continue

        loop_26: do i=0,nleg-1
           if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
             dens3(1)=dens3(1)*denevalmu2(nleg,i,q3(1,:),Vi,msq,czip)
             denst(3)=denst(3)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
           endif
        enddo loop_26

 26     continue

        xneval(1)=numeval(cut3,q3(1,:),czip)
        if     (imeth.eq.'diag') then
           known(1)=(xneval(1)-resi5(1)-resi4(1))/dens3(1)
        elseif (imeth.eq.'tree') then
           known(1)=xneval(1)-(resi5(1)+resi4(1))/dens3(1)
        endif
 
        c3(0)=known(1)
        do m=1,9
           c3(m)=czip
        enddo

     else
         if (abs(C0-1.0_ki) .lt. C0_thrs) then
            !---#[ New Sampling:
            ! The new sampling is the one that is safe around C0=1
            ! but not around C0=0

            ! traditional system
            mu2vec=(/czip,czip,czip,czip,czip,czip,czip,mu2g(3),mu2g(3),mu2g(3)/)
            nsol=10
            !---#] New Sampling:
         else
            !---#[ Old Sampling:
            ! The old sampling is the one that is safe around C0=0
            ! but not around C0=1
            if (diff.eq.2) then
               ! rank1 c-system: 3 coefficients
               mu2vec=(/czip,czip,czip,czip,czip,czip,czip,czip,czip,czip/)
               nsol=3
            elseif (diff.eq.1) then
               ! rank1 c-system: 3 coefficients
               mu2vec=(/czip,czip,czip,czip,czip,mu2g(3),czip,czip,czip,czip/)
               nsol=6
            else
               ! traditional system
               mu2vec=(/czip,czip,czip,czip,czip,czip,czip,mu2g(3),mu2g(3),mu2g(3)/)
               nsol=10
            endif
            !---#] Old Sampling:
         endif
      

        if_nleg: if     (nleg.eq.5) then
           do n=1,nsol
              resi5(n)=Res5(1,mu2vec(n))
           enddo
           resit(3)=Res5(1,mu2t(3))         
           goto 111
        elseif (nleg.eq.4) then
           do n=1,nsol
              resi4(n)=Res4(1,q3(n,:),mu2vec(n))
           enddo
           resit(3)=Res4(1,qt,mu2t(3))         
           goto 121
        elseif (nleg.eq.3) then
           goto 126
        else
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
                                  & .and.(i.ne.i3).and.&
                                  &(i.ne.i4).and.(i.ne.i5)) then
                                if ((i.eq.j1).or.(i.eq.j2).or.(i.eq.j3)) then
                                   dens5(:)=czip
                                   dens5t=czip
                                   evalres=.false.
                                   exit loop_110
                                else
                                   do n=1,nsol
                                      dens5(n)=dens5(n)*denevalmu2(nleg,i,&
                                            &q3(n,:),Vi,msq,mu2vec(n))
                                   enddo
                                   dens5t=dens5t&
                                     &*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
                                   evalres=.true.
                                endif
                             endif
                          enddo loop_110

                          if (evalres) then
                             do n=1,nsol
                                resi5(n)=resi5(n)&
                                   &+dens5(n)*res5(dicut5,mu2vec(n))
                             enddo
                             resit(3)=resit(3)+dens5t*res5(dicut5,mu2t(3))
                          endif

                          dicut5=dicut5+1
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        endif if_nleg

 111    continue

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
                           if ((i.eq.j1).or.(i.eq.j2).or.(i.eq.j3)) then
                              dens4(:)=czip
                              dens4t=czip
  
                              evalres=.false.
                              exit loop_120
                           else
                              do n=1,nsol
                                 dens4(n)=dens4(n)&
                                   &*denevalmu2(nleg,i,q3(n,:),Vi,msq,mu2vec(n))
                              enddo
                              dens4t=dens4t&
                                  &*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
                              evalres=.true.
                           endif
                        endif
                     enddo loop_120

                     if (evalres) then
                        do n=1,nsol
                           resi4(n)=resi4(n)&
                                & +dens4(n)*Res4(dicut4,q3(n,:),mu2vec(n))
                        enddo
                        resit(3)=resit(3)+dens4t*Res4(dicut4,qt,mu2t(3))
                     endif

                     dicut4=dicut4+1
                  enddo
               enddo
            enddo
         enddo

 121     continue

         do i=0,nleg-1
            if ((i.ne.j1).and.(i.ne.j2).and.(i.ne.j3)) then
               do n=1,nsol
                  dens3(n)=dens3(n)*denevalmu2(nleg,i,q3(n,:),Vi,msq,mu2vec(n))
               enddo
               denst(3)=denst(3)*denevalmu2(nleg,i,qt,Vi,msq,mu2t(3))
            endif
         enddo

 126     continue

         do n=1,nsol
            xneval(n)=numeval(cut3,q3(n,:),mu2vec(n))
         enddo

         if     (imeth.eq.'diag') then
            known(:)=(xneval(:)-resi5(:)-resi4(:))/dens3(:)
         elseif (imeth.eq.'tree') then
            known(:)=xneval(:)-(resi5(:)+resi4(:))/dens3(:)
         endif

         if (abs(C0-1.0_ki) .lt. C0_thrs) then
            !---#[ New Sampling:
            ! The new sampling is the one that is safe around C0=1
            ! but not around C0=0

            ! traditional system
        
            do m=0,6
               f3(m)=effe(known,1,7,m)
            enddo

            do m=7,9
               mx=m-7
               f3(m)=effe(known,8,3,mx)
            enddo

            !---#[ getc3_S3:
            c3(0) =   f3(0)
            c3(1) = -(f3(6)/(C0*MP12(3)))
            c3(2) =   f3(5)/(C0**2*MP12(3)**2)
            c3(3) = -(f3(4)/(C0**3*MP12(3)**3))
            c3(4) = -(f3(1)/MP12(3))
            c3(5) =   f3(2)/MP12(3)**2
            c3(6) = -(f3(3)/MP12(3)**3)
            c3(7) = (-f3(0) + MP12(3)**3*(-((C1**3*f3(3))/MP12(3)**3) -&
            &f3(4)/(C0**3*MP12(3)**3)) + f3(7))/mu2g(3)
            c3(8) = (C1**2*f3(2) + f3(6)/C0 - f3(9))/(MP12(3)*mu2g(3))
            c3(9) = (C1*f3(1) + f3(5)/C0**2 - f3(8))/(C1*MP12(3)*mu2g(3))
            !---#] getc3_S3:  

            !---#] New Sampling:
         else
            !---#[ Old Sampling:
            ! The old sampling is the one that is safe around C0=0
            ! but not around C0=1
            select case(diff)
            case(2)
               ! rank1 c-system: 3 coefficients
               f3(0:2)=known(1:3)
               !---#[ getc3_S1:  
               !----#[ original code:
      !      c3(0) = (f3(0) + f3(1))/two
      !      c3(1) = -(KK(3)*((one + C0)*f3(0) + f3(1) - C0*f3(1) - two*f3(2)))/ &
      !        (two*(-one + C0**2)*MP12(3))
      !      c3(4) = &
      !       -((KK(3)*c3(0) + C0*MP12(3)*c3(1) - KK(3)*f3(1))/(KK(3)**2*MP12(3)))
               !----#] original code:
               !----#[ HAGGIES:
               t1 = f3(0)
               t2 = f3(1)
               c3(0) = ((t1+t2)*0.5_ki)
               t3 = KK(3)
               t4 = MP12(3)
               c3(1) = (-1.0_ki)*(t2+((1.0_ki)+C0)*t1-(t2*C0+(2.0_ki)*f3(2)))/((C&
               &0*C0-(1.0_ki))*t4)*t3*0.5_ki
               c3(4) = (-1.0_ki)*(t3*c3(0)+t4*C0*c3(1)-t2*t3)/(t3*t3*t4)
               !----#] HAGGIES:
               !---#] getc3_S1:  
            case(1)
               ! rank1 c-system: 3 coefficients
               do m=0,2
                  f3(m)=effe(known,1,3,m)
               enddo
               
               do m=3,4
                  mx=m-3
                  f3(m)=effe(known,4,2,mx)
               enddo
               
               do m=5,5
                  mx=m-5
                  f3(m)=effe(known,6,1,mx)
               enddo
               !---#[ getc3_S2:
               !----#[ original code:
      !!$      c3(0) = f3(0)
      !!$      c3(2) = &
      !!$     &        (KK(3)**2*(f3(0) - f3(3) + &
      !!$     &      C0**2*(C0**2*f3(1) + f3(2) - C0*f3(4))))/&
      !!$     &  ((-one + C0**6)*MP12(3)**2)
      !!$      c3(4) = (C0**2*MP12(3)*c3(2))/KK(3)**3 - f3(1)/(KK(3)*MP12(3))
      !!$      c3(1) = -((KK(3)*(C0*KK(3)*MP12(3)*c3(4) + f3(4)))/MP12(3))
      !!$      c3(5) = (C0*MP12(3)*c3(1) + KK(3)*f3(2))/(KK(3)**3*MP12(3)**2)
      !!$      c3(7) = &
      !!$     &  (MP12(3)*c3(1) - MP12(3)**2*c3(2) + C1*MP12(3)*c3(4)-& 
      !!$     &    C1**2*MP12(3)**2*c3(5) + (-c3(0) + f3(5)))/mu2g(3)
      
      !        c3(0) = f3(0)
      !     c3(2) = &
      !    &        (f3(0) - f3(3) + &
      !    &      C0**2*(C0**2*f3(1) + f3(2) - C0*f3(4)))/&
      !    &  ((-one + C0**6)*MP12(3)**2)
      !     c3(4) = C0**2*MP12(3)*c3(2) - f3(1)/MP12(3)
      !     c3(1) = -(C0*MP12(3)*c3(4) + f3(4))/MP12(3)
      !     c3(5) = (C0*MP12(3)*c3(1) + f3(2))/MP12(3)**2
      !     c3(7) = &
      !    &  (MP12(3)*c3(1) - MP12(3)**2*c3(2) + C1*MP12(3)*c3(4)-& 
      !    &    C1**2*MP12(3)**2*c3(5) + (-c3(0) + f3(5)))/mu2g(3)
                  !----#] original code:
                  !----#[ HAGGIES:
               t1 = f3(0)
               c3(0) = t1
               t2 = C0*C0
               t3 = t2*C0
   
               t4 = f3(2)
               t5 = f3(1)
               t6 = f3(4)
               t7 = MP12(3)
               t8 = t7*t7
               c3(2) = ((t1+(t4+t2*t5-t6*C0)*t2-f3(3))/((t3*t3-one)*t8))
               c3(4) = (t2*t7*c3(2)-t5/t7)
               t1 = t7*C0
               c3(1) = (-1.0_ki)*(t6+t1*c3(4))/t7
               c3(5) = ((t4+t1*c3(1))/t8)
               t1 = t7*C1
               c3(7) = ((f3(5)+t1*c3(4)+t7*c3(1)-(t8*c3(2)+t1*t1*c3(5)+c3(0)))/(m&
               &u2g(3)))
               !----#] HAGGIES:
               !---#] getc3_S2:  
            case default 
               ! traditional system
              
               do m=0,3
                  f3(m)=effe(known,1,4,m)
               enddo
      
               do m=4,6
                  mx=m-4
                  f3(m)=effe(known,5,3,mx)
               enddo
   
               do m=7,9
                  mx=m-7
                  f3(m)=effe(known,8,3,mx)
               enddo
      
               !---#[ getc3_S3:
               !----#[ original code:
      !     c3(0) =  &
      !    &f3(0)
      !     c3(1) =  &
      !    &        -((KK(3)*(C0**5*f3(1) + C0**2*f3(2) + C0**11*f3(3) +  &
      !    &        C0**8*(f3(0) - f3(4)) - C0**4*f3(5) - f3(6)))/ &
      !    &    ((-one + C0**12)*MP12(3)))
      !     c3(2) =  &
      !    &        (KK(3)**2*(C0*f3(1) + C0**10*f3(2) + C0**7*f3(3) +  &
      !    &      C0**4*(f3(0) - f3(4)) - f3(5) - C0**8*f3(6)))/ &
      !    &  ((-one + C0**12)*MP12(3)**2)
      !     c3(3) =  &
      !    &        -((KK(3)**3*(f3(0) + C0**9*f3(1) + C0**6*f3(2) +  &
      !    &        C0**3*f3(3) - f3(4) - C0**8*f3(5) - C0**4*f3(6)))/ &
      !    &    ((-one + C0**12)*MP12(3)**3))
      !     c3(4) =  &
      !    &        (f3(1) + C0**9*f3(2) + C0**6*f3(3) +  &
      !    &    C0**3*(f3(0) - f3(4)) - C0**11*f3(5) - C0**7*f3(6))/ &
      !    &  ((-one + C0**12)*KK(3)*MP12(3))
      !     c3(5) =  &
      !    &        (-(C0**3*f3(1)) - f3(2) - C0**9*f3(3) +  &
      !    &    C0**6*(-f3(0) + f3(4)) + C0**2*f3(5) + C0**10*f3(6))/ &
      !    &  ((-one + C0**12)*KK(3)**2*MP12(3)**2)
      !     c3(6) =  &
      !    &        (C0**6*f3(1) + C0**3*f3(2) + f3(3) +  &
      !    &    C0**9*(f3(0) - f3(4)) - C0**5*f3(5) - C0*f3(6))/ &
      !    &  ((-one + C0**12)*KK(3)**3*MP12(3)**3)
      !     c3(7) =  &
      !    &        ((C1**3*KK(3)**3*(C0**6*f3(1) + C0**3*f3(2) + f3(3) +  &
      !    &         C0**9*(f3(0) - f3(4)) - C0**5*f3(5) - C0*f3(6)))/ &
      !    &     (-one + C0**12) -  &
      !    &    (KK(3)**3*(f3(0) + C0**9*f3(1) + C0**6*f3(2) +  &
      !    &         C0**3*f3(3) - f3(4) - C0**8*f3(5) - C0**4*f3(6)))/ &
      !    &     (-one + C0**12) + KK(3)**3*(-f3(0) + f3(7)))/(KK(3)**3*mu2g(3))
      !     c3(8) =  &
      !    &        -((-((KK(3)*(C0**5*f3(1) + C0**2*f3(2) + C0**11*f3(3) +  &
      !    &             C0**8*(f3(0) - f3(4)) - C0**4*f3(5) - f3(6)))/ &
      !    &         (-one + C0**12)) -  &
      !    &      (C1**2*KK(3)*(-(C0**3*f3(1)) - f3(2) - C0**9*f3(3) +  &
      !    &           C0**6*(-f3(0) + f3(4)) + C0**2*f3(5) +  &
      !    &           C0**10*f3(6)))/(-one + C0**12) + KK(3)*f3(9))/ &
      !    &    (MP12(3)*mu2g(3)))
      !     c3(9) =  &
      !    &        (-((C1*KK(3)**2*(f3(1) + C0**9*f3(2) + C0**6*f3(3) +  &
      !    &           C0**3*(f3(0) - f3(4)) - C0**11*f3(5) -  &
      !    &           C0**7*f3(6)))/(-one + C0**12)) +  &
      !    &    (KK(3)**2*(C0*f3(1) + C0**10*f3(2) + C0**7*f3(3) +  &
      !    &         C0**4*(f3(0) - f3(4)) - f3(5) - C0**8*f3(6)))/ &
      !    &     (-one + C0**12) - KK(3)**2*f3(8))/(C1*KK(3)**3*MP12(3)*mu2g(3))
               !----#] original code:
               !----#[ HAGGIES:
               t1 = f3(0)
               c3(0) = t1
               t2 = C0*C0
               t3 = t2*C0
               t4 = t3*t3
               t5 = f3(6)
               t6 = f3(4)
               t7 = t1-t6
               t8 = t2*t2
               t9 = t8*t8
               t10 = f3(1)
               t11 = t10*C0
               t12 = f3(2)
               t13 = f3(3)
               t14 = t13*C0
               t15 = t8*C0
               t15 = t15*t15
               t16 = f3(5)
               t17 = KK(3)
               t18 = MP12(3)
               t19 = (t11*t8+t12*t2+t14*t15+t7*t9-(t16*t8+t5))*t17
               t3 = (t3-(1.0_ki))*((1.0_ki)+t3)*((1.0_ki)+t4)
               t20 = t18*t3
               c3(1) = (-1.0_ki)*t19/t20
               t21 = t17*t17
               t22 = (t11+t12*t15+t14*t4+t7*t8-(t5*t9+t16))*t21
               t23 = t18*t18
               c3(2) = (t22/(t23*t3))
               t24 = (t7+t11*t9+t12*t4+t14*t2-(t5*t8+t16*t9))*t17*t21
               c3(3) = (-1.0_ki)*t24/(t20*t23)
               t7 = t7*C0
               t20 = t12*C0
               t23 = t16*C0
               t25 = t5*C0
               t26 = t10+t13*t4+t2*t7+t20*t9-(t25*t4+t15*t23)
               t27 = t17*t18*t3
               c3(4) = (t26/t27)
               t5 = (t6-t1)*t4+t15*t5+t16*t2-(t14*t9+t11*t2+t12)
               t6 = t17*t18
               t11 = t6*t6
               c3(5) = (t5/(t11*t3))
               t2 = t13+t10*t4+t2*t20+t7*t9-(t23*t8+t25)
               c3(6) = (t2/(t11*t27))
               t4 = t17*C1
               t7 = mu2g(3)
               c3(7) = (((f3(7)-t1)*t17*t21+t17*t2/t3*t4*t4*C1-t24/t3)/(t17*t21*t&
               &7))
               c3(8) = (-1.0_ki)*(f3(9)*t17-(t17/t3*t5*C1*C1+t19/t3))/(t18*t7)
               c3(9) = ((t22/t3-(t21*t26/t3*C1+f3(8)*t21))/(t21*t6*t7*C1))
               !----#] HAGGIES:
               !---#] getc3_S3:  
            end select
            !---#] Old Sampling:
         endif
      end if if_diff

      if (diff.ge.1) then
         c3(3)=czip
         c3(6)=czip
         c3(8)=czip
         c3(9)=czip 
         if (diff.ge.2) then
            c3(2)=czip
            c3(5)=czip
            c3(7)=czip
            if (diff.ge.3) then
               c3(1)=czip
               c3(4)=czip
            endif 
         endif 
      endif
   end subroutine getc3_rm

end module mgetc3

