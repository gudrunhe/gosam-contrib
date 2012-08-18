module mglobal
   use precision, only: ki
   implicit none

   private :: ki


   !mgetc1.f90:      common/FGzip/G0,KK,mu2,MP12,mu2t
   ! G0   --->    G0
   ! KK           KK(1)
   ! mu2          mu2g(1)
   ! MP12         MP12(1)
   ! mu2t         mu2t(1)
   !mgetc1.f90:      common/t1/resit,dens1t,mu2test
   ! resit ---->  resit(1)
   ! dens1t --->  denst(1)
   ! mu2test -->  mu2test(1)

   !mgetc2.f90:      common/Fzip/Fp,Fz,Fm,F1z,mu2,KK,Kmu,KB,MP12,mu2t
   ! Fp,Fz,Fm,F1z --> Fp,Fz,Fm,F1z
   ! KB   --->    KB
   ! mu2  --->    mu2g(2)
   ! KK   --->    KK(2)
   ! Kmu  -->     Kmu(2)
   ! MP12 --->    MP12(2)
   ! mu2t --->    mu2t(2)
   !mgetc2.f90:      common/t2/resit,dens2t,mu2test
   ! resit --> resit(2)
   ! dens2t --> denst(2)
   ! mu2test --> mu2test(2)


   !mgetc3.f90:      common/C0C1/C0,C1,mu2,MP12,KK,Kmu,mu2t
   ! C0,C1 ---> C0,C1
   ! mu2 ----> mu2g(3)
   ! MP12 ---> MP12(3)
   ! KK  ----> KK(3)
   ! Kmu ---> Kmu(3)
   ! mu2t ---> mu2t(3)
   !mgetc3.f90:      common/t3/resit,dens3t,mu2test
   ! resit ---> resit(3)
   ! dens3t --> denst(3)
   ! mu2test  -> mu2test(3)

   !mgetqs.f90:      common/FGzip/G0,KK,mu2,MP12,mu2t
   ! G0   --->    G0
   ! KK           KK(1)
   ! mu2          mu2g(1)
   ! MP12         MP12(1)
   ! mu2t         mu2t(1)
   !mgetqs.f90:      common/Fzip/Fp,Fz,Fm,F1z,mu2,KK,Kmu,KB,MP12,mu2t
   ! Fp,Fz,Fm,F1z --> Fp,Fz,Fm,F1z
   ! KB   --->    KB
   ! mu2  --->    mu2g(2)
   ! KK   --->    KK(2)
   ! Kmu  -->     Kmu(2)
   ! MP12 --->    MP12(2)
   ! mu2t --->    mu2t(2)
   !mgetqs.f90:      common/ds/dx1,dx2,dx3,dx4,dx5
   ! dxi --> dx(i)
   !mgetqs.f90:      common/C0C1/C0,C1,mu2,MP12,KK,Kmu,mu2t
   ! C0,C1 ---> C0,C1
   ! mu2 ----> mu2g(3)
   ! MP12 ---> MP12(3)
   ! KK  ----> KK(3)
   ! Kmu ---> Kmu(3)
   ! mu2t ---> mu2t(3)
   !mgetqs.f90:      common/mp12mu2/MP12,mu2,mu2t
   ! mp12 ---> MP12(4)
   ! mu2 ---> mu2g(4)
   ! mu2t ---> mu2t(4)

   !mgetc4.f90:      common/ds/dx1,dx2,dx3,dx4,dx5
   ! dxi --> dx(i)
   !mgetc4.f90:      common/mp12mu2/mp12,mu2,mu2t
   ! mp12 ---> MP12(4)
   ! mu2 ---> mu2g(4)
   ! mu2t ---> mu2t(4)
   !mgetc4.f90:      common/t4/resi5t,dens4t,mu2test
   ! resi5t ---> resit(4)
   ! dens4t ---> denst(4)
   ! mu2test --> mu2test(4)

   real(ki), dimension(1:4) :: MP12
   real(ki), dimension(1:3) :: KK
   complex(ki), dimension(1:4) :: mu2g, mu2t, mu2test
   complex(ki), dimension(1:4) :: resit, denst
   real(ki), dimension(2:3) :: Kmu
   complex(ki), dimension(1:5) :: dx
   real(ki) :: Fp,Fz,Fm,F1z,KB,C0,C1,G0
   complex(ki) :: Fpc,Fzc,Fmc,F1zc,C0c,C1c,G0c

end module mglobal
