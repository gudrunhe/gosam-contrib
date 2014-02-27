! 
!****h* src/kinematic/inverse_matrice
! NAME
!
!  Module inverse_matrice
!
! USAGE
!
!  use inverse_matrice
!
! DESCRIPTION
!
! This module provides some routines and tools to inverse a n x n matrix. 
!
! OUTPUT
!
!  This module exports two routines:
!  * inverse -- to inverse a nXn matrix
!  * imprime_mat -- to print a nXn matrix
!
! USES
!
!  * precision (src/module/precision_golem.f90)
!  * equal (src/module/equal.f90)
!  * s_matrix_type (src/module/s_matrix_type.f90)
!  * constante (src/module/constante.f90)
!
!*****
module inverse_matrice
  use precision_golem
  use equal
  use s_matrix_type
  use parametre, only : accuracy_par
  use sortie_erreur, only : tab_erreur_par,catch_exception
  use constante, only:czero
  implicit none
  !
  private
  !
  interface imprime_mat
     !
     module procedure imprime_mat_r, imprime_mat_c
     module procedure imprime_mat_p
     !
  end interface
  !
  interface inverse
     !
     module procedure inverse_r, inverse_c
     module procedure inverse_pr, inverse_pc
     !
  end interface
  !
  interface inverse_rescue
     !
     module procedure inverse_rescue_r, inverse_rescue_c
     !
  end interface
  !  
  interface inverse_true
     !
     module procedure inverse_true_r, inverse_true_c
     !
  end interface
  !
  interface lu_decomp
     !
     module procedure lu_decomp_r, lu_decomp_c
     !
  end interface
  !
  interface inverse_triangular
     !
     module procedure inverse_triangular_r, inverse_triangular_c
     !
  end interface
    !
  interface inverse_greville
     !
     module procedure inverse_greville_r, inverse_greville_c
     !
  end interface
  !
  interface compt
     !
     module procedure compt_r, compt_c
     !
  end interface
    !
  interface verif
     !
     module procedure verif_r, verif_c
     !
  end interface
  !
  public :: inverse, imprime_mat
  real(ki) :: glob_eps = 1.e-12_ki        ! valeur en de ca duquelle on
                                          ! passe au cas singulier
  contains
    !
    !****f* src/kinematic/inverse_matrice/inverse
    ! NAME
    !
    !  Subroutine inverse
    !
    ! USAGE
    !
    !  call inverse(mat,inv_mat,error,pinch1,pinch2,pinch3,pinch4,pinch5)
    !
    ! DESCRIPTION
    !
    ! This routine first tries the Gauss method with partial pivoting strategy. 
    ! If the error returned is too large (greater than the global variable accuracy_par),
    ! then it switches to another method : the Greville method.
    ! In the case of the Gauss method, if some reduced matrices need to be inverted, a new matrix is built
    ! by removing the row(s) and column(s) pinch1,pinch2, etc. then the inverse is computed and the result returned
    ! is a nXn matrix where the column(s) and row(s) pinch1, pinch2, etc. are filled by 0, the other elements
    ! are those of the inverse computed. In the Greville method, the reduce matrix which is a nXn matrix
    ! where the column(s) and row(s) pinch1, pinch2, etc. are filled by 0 is directly inverted.
    ! Note that the error is computed in the following way:
    ! first the matrix is rescaled : i. e. divided by the greatest (in absolute value) element
    ! then the inverse is computed and the two matrices abs(1 - A^(-1) A) and abs(1 - A A^(-1)) are computed
    ! the error is the greatest element of these two matrices.
    ! In the case of the Greville method, the Moore_Penrose conditions are also tested
    !
    ! INPUTS
    !
    !  * mat -- a real/complex (type ki) array of rank 2, or an s_matrix_poly type.
    !  * pinch1 -- an integer (optional), specified a pinch
    !  * pinch2 -- an integer (optional), specified a pinch
    !  * pinch3 -- an integer (optional), specified a pinch
    !  * pinch4 -- an integer (optional), specified a pinch
    !  * pinch5 -- an integer (optional), specified a pinch
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * inv_mat -- a real (type ki) array of rank 2, same shape, the inverse 
    !               of the matrix mat
    !  * error -- a real (type ki), the estimation of the error of the numerical inversion
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine inverse_pr(mat_p,inv_mat_r,error,pinch1,pinch2,pinch3,pinch4,pinch5)
      type(s_matrix_poly), intent(in) :: mat_p
      real(ki), intent(out), dimension(size(mat_p%pt_real,1),size(mat_p%pt_real,1)) :: inv_mat_r
      real(ki), intent(out) :: error     
      integer, optional :: pinch1, pinch2,pinch3,pinch4,pinch5
      !
      if (associated(mat_p%pt_real)) then
         call inverse_r(mat_p%pt_real,inv_mat_r,error, pinch1=pinch1, pinch2=pinch2,pinch3=pinch3,pinch4=pinch4,pinch5=pinch5)
      end if
      !
    end subroutine inverse_pr
    !
    subroutine inverse_pc(mat_p,inv_mat_c,error,pinch1,pinch2,pinch3,pinch4,pinch5)
      type(s_matrix_poly), intent(in) :: mat_p
      complex(ki), intent(out), dimension(size(mat_p%pt_cmplx,1),size(mat_p%pt_cmplx,1)) :: inv_mat_c
      real(ki), intent(out) :: error     
      integer, optional :: pinch1, pinch2,pinch3,pinch4,pinch5
      !
      if (associated(mat_p%pt_cmplx)) then
         call inverse_c(mat_p%pt_cmplx,inv_mat_c,error,pinch1=pinch1,pinch2=pinch2,pinch3=pinch3,pinch4=pinch4,pinch5=pinch5)
      end if
      !
    end subroutine inverse_pc
    !
    subroutine inverse_r(mat_r,inv_mat_r,error,pinch1,pinch2,pinch3,pinch4,pinch5)
      !
      real(ki), intent(in), dimension(:,:) :: mat_r
      real(ki), intent(out), dimension(size(mat_r,1),size(mat_r,1)) :: inv_mat_r
      real(ki), intent(out) :: error
      integer, optional :: pinch1,pinch2,pinch3,pinch4,pinch5
      logical :: invertible
      !
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: norm_mat_r, mat_greville
      real(ki) :: plus_grand,g_error,o_error
      integer :: pin1,pin2,pin3,pin4,pin5
      !
      pin1 = -1
      pin2 = -1
      pin3 = -1
      pin4 = -1
      pin5 = -1
      !
      if (present(pinch1)) pin1 = pinch1
      if (present(pinch2)) pin2 = pinch2
      if (present(pinch3)) pin3 = pinch3
      if (present(pinch4)) pin4 = pinch4
      if (present(pinch5)) pin5 = pinch5
      !
      g_error = 1._ki
      o_error = 1._ki
      !
      !
      ! First we rescale the matrix
      !
      plus_grand = maxval(array=abs(mat_r))
      if (equal_real(plus_grand,0._ki)) then
              plus_grand = 1._ki ! to prevent divison by zero
      end if
      norm_mat_r = mat_r/plus_grand
      !
      ! We first use the Gauss method
      !
      call inverse_rescue(norm_mat_r,inv_mat_r,o_error,invertible,pin1,pin2,pin3,pin4,pin5)
      !
      if ((o_error >= accuracy_par) .and. invertible) then
        !
        call inverse_greville(norm_mat_r,mat_greville,g_error)
        !
        if (g_error .lt. o_error) then
           inv_mat_r = mat_greville
        end if
        !
      end if
      !
      error = min(g_error,o_error)
      !
      if ((error >= accuracy_par) .and. invertible) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'the Greville method failed'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the Gauss method failed too'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'the error returned %f0'
        tab_erreur_par(3)%arg_real = error
        call catch_exception(1)
      else if (.not. invertible) then
        tab_erreur_par(1)%chaine = 'The Gauss method failed. Matrix not invertible.'
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'The error returned %f0'
        tab_erreur_par(2)%arg_real = error
        call catch_exception(1)
        !
      end if
      !
      inv_mat_r = inv_mat_r/plus_grand
      !
    end subroutine inverse_r
    !
    subroutine inverse_c(mat_c,inv_mat_c,error,pinch1,pinch2,pinch3,pinch4,pinch5)
      !
      complex(ki), intent(in), dimension(:,:) :: mat_c
      complex(ki), intent(out), dimension(size(mat_c,1),size(mat_c,1)) :: inv_mat_c
      real(ki), intent(out) :: error
      integer, optional :: pinch1,pinch2,pinch3,pinch4,pinch5
      !
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: norm_mat_c, mat_greville_c
      real(ki) :: plus_grand,g_error,o_error
      integer :: pin1,pin2,pin3,pin4,pin5
      logical :: invertible
      !
      pin1 = -1
      pin2 = -1
      pin3 = -1
      pin4 = -1
      pin5 = -1
      !
      if (present(pinch1)) pin1 = pinch1
      if (present(pinch2)) pin2 = pinch2
      if (present(pinch2)) pin3 = pinch3
      if (present(pinch2)) pin4 = pinch4
      if (present(pinch2)) pin5 = pinch5
      !
      g_error = 1._ki
      o_error = 1._ki
      !
      !
      ! First we rescale the matrix
      !
      plus_grand = max(maxval( array=abs( real(mat_c,ki) ) ), maxval( array=abs( aimag(mat_c) ) ) ) 
      if (equal_real(plus_grand,0._ki)) then
              plus_grand = 1._ki ! to prevent divison by zero
      end if
      norm_mat_c = mat_c/cmplx(plus_grand,0._ki,ki)
      !
      ! We first use the Gauss method
      !
      call inverse_rescue(norm_mat_c,inv_mat_c,o_error,invertible,pin1,pin2,pin3,pin4,pin5)
      !
      if ((o_error >= accuracy_par) .and. invertible) then
        !
        call inverse_greville(norm_mat_c,mat_greville_c,g_error)
        !
        if (g_error .lt. o_error) then
           inv_mat_c = mat_greville_c
        end if
      end if
      !
      error = min(g_error,o_error)
      !
      if ((error >= accuracy_par) .and. invertible) then
        !
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(1)%chaine = 'the Greville method failed'
        tab_erreur_par(2)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'the Gauss method failed too'
        tab_erreur_par(3)%a_imprimer = .true.
        tab_erreur_par(3)%chaine = 'the error returned %f0'
        tab_erreur_par(3)%arg_real = error
        call catch_exception(1)
        !
      else if (.not. invertible) then
        tab_erreur_par(1)%chaine = 'The Gauss method failed. Matrix not invertible.'
        tab_erreur_par(1)%a_imprimer = .true.
        tab_erreur_par(2)%chaine = 'The error returned %f0'
        tab_erreur_par(2)%arg_real = error
        call catch_exception(1)
        !
      end if
      !
      inv_mat_c = inv_mat_c/plus_grand
      !
    end subroutine inverse_c
    !****if* src/kinematic/inverse_matrice/inverse_rescue
    ! NAME
    !
    !  Subroutine inverse_rescue
    !
    ! USAGE
    !
    !  call inverse_rescue(mat,inv_mat,error,pinch1,pinch2)
    !
    ! DESCRIPTION
    !
    ! The role of this routine is just to reduce the size the input 
    ! matrix mat if pinch1 and/or pinch2 etc. is/are present
    !
    ! INPUTS
    !
    !  * mat -- a real/complex (type ki) array of rank 2
    !  * pinch1 -- an integer (optional), specified a pinch
    !  * pinch2 -- an integer (optional), specified a pinch
    !  * pinch3 -- an integer (optional), specified a pinch
    !  * pinch4 -- an integer (optional), specified a pinch
    !  * pinch5 -- an integer (optional), specified a pinch
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * inv_mat -- a real/compelx (type ki) array of rank 2, same shape, the inverse 
    !               of the matrix mat
    !  * error -- a real (type ki), the estimation of the error of the numerical inversion
    !  * invertible -- logical, true if the matrix is invertible
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine inverse_rescue_r(mat_r,inv_mat_r,error,invertible,pin1,pin2,pin3,pin4,pin5)
      !
      real(ki), intent(in), dimension(:,:) :: mat_r
      real(ki), intent(out), dimension(size(mat_r,1),size(mat_r,1)) :: inv_mat_r
      real(ki), intent(out) :: error
      logical, intent(out) :: invertible
      integer, intent(in) :: pin1,pin2,pin3,pin4,pin5
      !
      integer :: i,j,n_dim
      logical, dimension(size(mat_r,1),size(mat_r,1)) :: masque
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: f2
      real(ki), dimension(:,:), allocatable :: true_mat_r,true_inv_mat_r
      real(ki), dimension(:), allocatable :: true_vect_r
      integer, dimension(5) :: tab_pinch
      integer :: nb_pinch
      !
      tab_pinch = (/pin1,pin2,pin3,pin4,pin5/)
      nb_pinch = count(mask=tab_pinch /= -1)
      masque = .true.
      !
      n_dim = size(mat_r,1)-nb_pinch
      !
      allocate(true_mat_r(n_dim,n_dim))
      allocate(true_inv_mat_r(n_dim,n_dim))
      allocate(true_vect_r(n_dim*n_dim))
      !
      select case (nb_pinch)
      !
      case(0)
        !
        true_mat_r = mat_r
        !
      case default ! nb_pinch > 0
        !
        do i = 1, size(tab_pinch)
          j=tab_pinch(i)
          if (j /= -1) then
                masque(j,:) = .false.
                masque(:,j) = .false.
          end if
        end do
        true_vect_r = pack(mat_r,mask=masque)
        true_mat_r = reshape(source=true_vect_r,shape=(/n_dim,n_dim/))
      end select
      !
      call inverse_true(true_mat_r,true_inv_mat_r,error,invertible)
      !
      f2(:,:) = 0._ki
      true_vect_r = pack(true_inv_mat_r,.true.)
      inv_mat_r = unpack(true_vect_r,masque,f2)
      !
      deallocate(true_mat_r)
      deallocate(true_inv_mat_r)
      deallocate(true_vect_r)
      !
    end subroutine inverse_rescue_r
    !
   subroutine inverse_rescue_c(mat_c,inv_mat_c,error,invertible,pin1,pin2,pin3,pin4,pin5)
      !
      complex(ki), intent(in), dimension(:,:) :: mat_c
      complex(ki), intent(out), dimension(size(mat_c,1),size(mat_c,1)) :: inv_mat_c
      real(ki), intent(out) :: error
      logical, intent(out) :: invertible
      integer, intent(in) :: pin1,pin2,pin3,pin4,pin5
      !
      integer :: i,j,n_dim
      logical, dimension(size(mat_c,1),size(mat_c,1)) :: masque
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: f2
      complex(ki), dimension(:,:), allocatable :: true_mat_c,true_inv_mat_c
      complex(ki), dimension(:), allocatable :: true_vect_c
      integer, dimension(5) :: tab_pinch
      integer :: nb_pinch
      !
      tab_pinch = (/pin1,pin2,pin3,pin4,pin5/)
      nb_pinch = count(mask=tab_pinch /= -1)
      masque = .true.
      !
      n_dim = size(mat_c,1)-nb_pinch
      !
      allocate(true_mat_c(n_dim,n_dim))
      allocate(true_inv_mat_c(n_dim,n_dim))
      allocate(true_vect_c(n_dim*n_dim))
      !
      select case (nb_pinch)
      !
      case(0)
        !
        true_mat_c = mat_c
        !
      case default ! nb_pinch > 0
        !
        do i = 1, size(tab_pinch)
          j=tab_pinch(i)
          if (j /= -1) then
                masque(j,:) = .false.
                masque(:,j) = .false.
          end if
        end do
        true_vect_c = pack(mat_c,mask=masque)
        true_mat_c = reshape(source=true_vect_c,shape=(/n_dim,n_dim/))
      end select
      !
      call inverse_true(true_mat_c,true_inv_mat_c,error,invertible)
      !
      f2(:,:) = czero
      true_vect_c = pack(true_inv_mat_c,.true.)
      inv_mat_c = unpack(true_vect_c,masque,f2)
      !
      deallocate(true_mat_c)
      deallocate(true_inv_mat_c)
      deallocate(true_vect_c)
      !
    end subroutine inverse_rescue_c
    !
    !****if* src/kinematic/inverse_matrice/inverse_true
    ! NAME
    !
    !  Subroutine inverse_true
    !
    ! USAGE
    !
    !  call inverse_true(mat,inv_mat,error,invertible)
    !
    ! DESCRIPTION
    !
    ! This routine uses a Gauss pivot method with partial pivoting to inverse M a nXn matrix.
    ! It returns an estimation of the error by computing the two matrices abs(1 - M^(-1) M) 
    ! and abs(1 - M M^(-1)). The error is the greatest element of these two matrices.
    !
    ! INPUTS
    !
    !  * mat -- a real/complex (type ki) array of rank 2
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * inv_mat -- a real/complex (type ki) array of rank 2, same shape, the inverse 
    !               of the matrix mat
    !  * error -- a real (type ki), the estimation of the error of the numerical inversion
    !  * invertible -- logical, true if mat is invertible
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine inverse_true_r(mat_r,inv_mat_r,error,invertible)
      !
      real(ki), intent(in), dimension(:,:) :: mat_r
      real(ki), intent(out), dimension(size(mat_r,1),size(mat_r,1)) :: inv_mat_r
      real(ki), intent(out) :: error
      logical, intent(out) :: invertible
      !
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: mat1,mat2,unit_mat
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: p_mat,l_mat,u_mat
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: inv_l_mat,inv_u_mat
      integer :: i,n_dim
      real(ki) :: max1,max2
      !integer :: errorflag
      !
      n_dim = size(mat_r,1)               ! dimension de la matrice
      unit_mat(:,:) = 0._ki
      !
      do i=1,n_dim
        !
        unit_mat(i,i) = 1._ki
        !
      end do
      !
      call lu_decomp(mat_r,p_mat,l_mat,u_mat,invertible)
      if (.not. invertible) then
              inv_mat_r = 0._ki
              error = 1._ki
              return
      end if
      call inverse_triangular(l_mat,'inf',inv_l_mat)
      call inverse_triangular(u_mat,'sup',inv_u_mat)
      inv_mat_r = matmul(inv_u_mat,inv_l_mat)
      inv_mat_r = matmul(inv_mat_r,transpose(p_mat))
      !
      mat1 = matmul(inv_mat_r,mat_r)
      mat2 = matmul(mat_r,inv_mat_r)
      mat1 = abs(mat1-unit_mat)
      mat2 = abs(mat2-unit_mat)
      !
      max1 = maxval(mat1)
      max2 = maxval(mat2)
      !
      error = max(max1,max2)
      !
    end subroutine inverse_true_r
    !
    subroutine inverse_true_c(mat_c,inv_mat_c,error,invertible)
      !
      complex(ki), intent(in), dimension(:,:) :: mat_c
      complex(ki), intent(out), dimension(size(mat_c,1),size(mat_c,1)) :: inv_mat_c
      real(ki), intent(out) :: error
      logical, intent(out) :: invertible
      !
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: mat1c,mat2c,unit_mat
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: p_mat,l_mat,u_mat
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: inv_l_mat,inv_u_mat
      real(ki), dimension(size(mat_c,1),size(mat_c,1)) :: mat1, mat2
      integer :: i,n_dim
      real(ki) :: max1,max2
      !integer :: errorflag
      !
      n_dim = size(mat_c,1)               ! dimension de la matrice
      unit_mat(:,:) = czero
      !
      do i=1,n_dim
        !
        unit_mat(i,i) = cmplx(1._ki,0._ki,ki)
        !
      end do
      !
      call lu_decomp(mat_c,p_mat,l_mat,u_mat,invertible)
      if (.not. invertible) then
              inv_mat_c = 0._ki
              error = 1._ki
              return
      end if
      call inverse_triangular(l_mat,'inf',inv_l_mat)
      call inverse_triangular(u_mat,'sup',inv_u_mat)
      inv_mat_c = matmul(inv_u_mat,inv_l_mat)
      inv_mat_c = matmul(inv_mat_c,transpose(p_mat))
      !
      mat1c = matmul(inv_mat_c,mat_c)
      mat2c = matmul(mat_c,inv_mat_c)
      mat1 = abs(mat1c-unit_mat)
      mat2 = abs(mat2c-unit_mat)
      !
      max1 = maxval(mat1)
      max2 = maxval(mat2)
      !
      error = max(max1,max2)
      !
    end subroutine inverse_true_c
    !
    !****if* src/kinematic/inverse_matrice/inverse_greville
    ! NAME
    !
    !  Subroutine inverse_greville
    !
    ! USAGE
    !
    !  call inverse_greville(mat,inv_mat,error)
    !
    ! DESCRIPTION
    !
    !  This routine inverses a n x n matrix using the Greville method. 
    !  This method enables to invert singular matrix using the Moore-Penrose definition.
    !  One builds an iterative process, a rectengular matrix A is defined with the 
    !  first row of the n x n matrix, the the other rows are added step by step
    !  At each step, the pseudo inverse (Moore-Penrose) A' is built such that
    !  A A' = transpose(A A')
    !  A' A = transpose(A' A)
    !  A' A A' = A'
    !  A A' A = A
    !
    ! INPUTS
    !
    !  * mat -- a real/complex (type ki) array of rank 2
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * inv_mat -- a real/complex (type ki) array of rank 2, same shape, the inverse 
    !               of the matrix mat
    !  * error -- a real (type ki), the estimation of the error of the numerical inversion
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    ! On construit un processus iteratif : on definit une matrice A qui au 
    ! debut est construit en partant de la premiere colonne
    ! puis en ajoutant les colonnes une a une. 
    ! Pour chaque etape, on construit le pseudo inverse 
    ! (au sens de Moore-Penrose) A' tel que :
    ! A A' = transpose(A A')
    ! A' A = transpose(A' A)
    ! A' A A' = A'
    ! A A' A = A
    ! cette routine retourne une estimation de l'erreur
    !
    subroutine inverse_greville_r(mat_r,inv_mat_r,error)
      !
      real(ki), intent(in), dimension(:,:) :: mat_r
      real(ki), intent(out), dimension(size(mat_r,1),size(mat_r,1)) :: inv_mat_r
      real(ki), intent(out) :: error
      !
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: id_mat, mat1, mat2
      real(ki), dimension(size(mat_r,1),1) :: a
      real(ki), dimension(1,size(mat_r,1)) :: at
      real(ki), dimension(:,:), allocatable :: ga,ga_prime
      real(ki), dimension(:,:), allocatable :: gab,gab_prime
      integer :: k,n_dim
      integer :: res
      real(ki) :: denom, max1, max2
      !
      n_dim = size(mat_r,1)               ! dimension de la matrice
      !
      !
      id_mat = 0._ki
      do k=1,n_dim
         id_mat(k,k) = 1._ki
      end do
      !
      ! premiere iteration
      !
      allocate(ga(n_dim,1),stat=res)
      allocate(ga_prime(1,n_dim),stat=res)
      !
      a(:,1) = mat_r(:,1)         ! premiere colonne de mat
      ga = a                         ! A est la 1ere colonne de mat
      at = transpose(a)
      denom = sum(matmul(at,a))         ! on utilise sum pour rendre scalaire
                                  ! matmul, sinon tableau (1,1)
      if (sqrt(denom) >= glob_eps) then
        !
        ga_prime = at/denom             ! cas standard A' pseudo-inverse de A
        !
      else
        !
        ga_prime(:,:) = 0._ki                 ! si mat est singuliere
        !
      end if
      !
      ! autres iterations
      !
      do k=2,n_dim
        !
        allocate(gab(n_dim,k),stat=res)
        allocate(gab_prime(k,n_dim),stat=res)
        !
        a(:,1) = mat_r(:,k)
        call compt(n_dim,k-1,a,ga,ga_prime,gab,gab_prime)
!        call verif(n_dim,k,gab,gab_prime,error)
        !
        deallocate(ga,stat=res)
        deallocate(ga_prime,stat=res)
        !
        allocate(ga(n_dim,k),stat=res)
        allocate(ga_prime(k,n_dim),stat=res)
        !
        ga = gab
        ga_prime = gab_prime
        !
        deallocate(gab,stat=res)
        deallocate(gab_prime,stat=res)
        !
      end do
      !
      inv_mat_r = ga_prime
      !
      mat1 = matmul(inv_mat_r,mat_r)
      mat2 = matmul(mat_r,inv_mat_r)
      mat1 = abs(mat1-id_mat)
      mat2 = abs(mat2-id_mat)
      !
      max1 = maxval(mat1)
      max2 = maxval(mat2)
      !
      error = max(max1,max2)
      !
    end subroutine inverse_greville_r
    !
    subroutine inverse_greville_c(mat_c,inv_mat_c,error)
      !
      complex(ki), intent(in), dimension(:,:) :: mat_c
      complex(ki), intent(out), dimension(size(mat_c,1),size(mat_c,1)) :: inv_mat_c
      real(ki), intent(out) :: error
      !
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: id_mat, mat1, mat2
      complex(ki), dimension(size(mat_c,1),1) :: a
      complex(ki), dimension(1,size(mat_c,1)) :: at
      complex(ki), dimension(:,:), allocatable :: ga,ga_prime
      complex(ki), dimension(:,:), allocatable :: gab,gab_prime
      integer :: k,n_dim
      integer :: res
      real(ki) :: denom, max1, max2
      !
      n_dim = size(mat_c,1)               ! dimension de la matrice
      !
      id_mat = czero
      do k=1,n_dim
         id_mat(k,k) = cmplx(1._ki,0._ki,ki)
      end do
      !
      ! premiere iteration
      !
      allocate(ga(n_dim,1),stat=res)
      allocate(ga_prime(1,n_dim),stat=res)
      !
      a(:,1) = mat_c(:,1)         ! premiere colonne de mat
      ga = a                         ! A est la 1ere colonne de mat
      at = conjg(transpose(a))
      denom = real(sum(matmul(at,a)),ki)       ! on utilise sum pour rendre scalaire
                                  ! matmul, sinon tableau (1,1)
      if (sqrt(denom) >= glob_eps) then
        !
        ga_prime = at/denom             ! cas standard A' pseudo-inverse de A
        !
      else
        !
        ga_prime(:,:) = czero               ! si mat est singuliere
        !
      end if
      !
      ! autres iterations
      !
      do k=2, n_dim
        !
        allocate(gab(n_dim,k),stat=res)
        allocate(gab_prime(k,n_dim),stat=res)
        !
        a(:,1) = mat_c(:,k)
        call compt(n_dim,k-1,a,ga,ga_prime,gab,gab_prime)
!        call verif(n_dim,k,gab,gab_prime,error)
        !
        deallocate(ga,stat=res)
        deallocate(ga_prime,stat=res)
        !
        allocate(ga(n_dim,k),stat=res)
        allocate(ga_prime(k,n_dim),stat=res)
        !
        ga = gab
        ga_prime = gab_prime
        !
        deallocate(gab,stat=res)
        deallocate(gab_prime,stat=res)
        !
      end do
      !
      inv_mat_c = ga_prime
      !
      mat1 = matmul(inv_mat_c,mat_c)
      mat2 = matmul(mat_c,inv_mat_c)
      mat1 = abs(mat1-id_mat)
      mat2 = abs(mat2-id_mat)
      !
      max1 = maxval(real(mat1,ki))
      max2 = maxval(real(mat2,ki))
      !
      error = max(max1,max2)
    end subroutine inverse_greville_c
    !
    !****if* src/kinematic/inverse_matrice/compt
    ! NAME
    !
    !  Subroutine compt
    !
    ! USAGE
    !
    !  call compt(n,k,a,ga,ga_prime,ga_plus,ga_prime_plus)
    !
    ! DESCRIPTION
    !
    !  This routine is used to iterate the inversion processus
    !
    ! INPUTS
    !
    !  * n -- an integer, the shape of the matrix to inverse
    !  * k -- an integer, number of rows
    !  * a -- a real/complex (type ki) array of rank 2,  its shape is (n,1)
    !  * ga -- a real/complex (type ki) array of rank 2, its shapes is (n,k) 
    !  * ga_prime -- a real/complex (type ki) array of rank 2, its shapes is (k,n)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * ga_plus -- a real/complex(type ki) array of rank 2, its shape is (n,k+1)
    !  * ga_prime_plus -- a real/complex(type ki) array of rank 2, its shape is (k+1,n)
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    ! cette routine itere le processus, elle recoit en argument les dimensions 
    ! de A at A' ainsi que A et A' et retourne A+ et A'+, A+ a une colonne de 
    ! plus que A et A'+ une ligne de plus que A'
    !
    subroutine compt_r(n,k,a,ga,ga_prime,ga_plus,ga_prime_plus)
      !
      integer, intent(in) :: n,k
      real(ki), intent(in), dimension(n,1) :: a
      real(ki), intent(in), dimension(n,k) :: ga
      real(ki), intent(in), dimension(k,n) :: ga_prime
      real(ki), intent(out), dimension(n,k+1) :: ga_plus
      real(ki), intent(out), dimension(k+1,n) :: ga_prime_plus
      !
      real(ki), dimension(k,1) :: d
      real(ki), dimension(1,k) :: dt
      real(ki), dimension(n,1) :: c
      real(ki), dimension(1,n) :: ct
      real(ki) :: denom
      !
      !
      d = matmul(ga_prime,a)
      dt = transpose(d)
      c = a - matmul(ga,d)
      ct = transpose(c)                         ! cas standard
      denom = sum(matmul(ct,c))         ! cas standard
      !
      if (sqrt(denom) < glob_eps) then
        !
        ct = matmul(dt,ga_prime)        ! cas singulier
        denom = 1._ki + sum(matmul(dt,d))       ! cas singulier
        !
      end if
      !
      ga_prime_plus(1:k,:) = ga_prime - matmul(d,ct)/denom
      ga_prime_plus(k+1,:) = ct(1,:)/denom
      ga_plus(:,1:k) = ga
      ga_plus(:,k+1) = a(:,1)
      !
    end subroutine compt_r

    subroutine compt_c(n,k,a,ga,ga_prime,ga_plus,ga_prime_plus)
      !
      integer, intent(in) :: n,k
      complex(ki), intent(in), dimension(n,1) :: a
      complex(ki), intent(in), dimension(n,k) :: ga
      complex(ki), intent(in), dimension(k,n) :: ga_prime
      complex(ki), intent(out), dimension(n,k+1) :: ga_plus
      complex(ki), intent(out), dimension(k+1,n) :: ga_prime_plus
      !
      complex(ki), dimension(k,1) :: d
      complex(ki), dimension(1,k) :: dt
      complex(ki), dimension(n,1) :: c
      complex(ki), dimension(1,n) :: ct
      real(ki) :: denom
      !
      d = matmul(ga_prime,a)
      dt = conjg(transpose(d))
      c = a - matmul(ga,d)
      ct = conjg(transpose(c))                        ! cas standard
      denom = real(sum(matmul(ct,c)),ki)         ! cas standard
      !
      if (sqrt(abs(denom)) < glob_eps) then
        !
        ct = matmul(dt,ga_prime)        ! cas singulier
        denom = 1._ki + real(sum(matmul(dt,d)),ki)       ! cas singulier
        !
      end if
      !
      ga_prime_plus(1:k,:) = ga_prime - matmul(d,ct)/denom
      ga_prime_plus(k+1,:) = ct(1,:)/denom
      ga_plus(:,1:k) = ga
      ga_plus(:,k+1) = a(:,1)
      !
    end subroutine compt_c

    !
    !****if* src/kinematic/inverse_matrice/verif
    ! NAME
    !
    !  Subroutine verif
    !
    ! USAGE
    !
    !  call verif(n,k,ga,ga_prime,error)
    !
    ! DESCRIPTION
    !
    !  This routine verifies the Moor-Penrose conditions, it returns the 
    !  maximum error obtained
    !
    ! INPUTS
    !
    !  * n -- an integer, the shape of the matrix to inverse
    !  * k -- an integer, number of rows
    !  * ga -- a real/complex(type ki) array of rank 2, its shapes is (n,k) 
    !  * ga_prime -- a real/complex(type ki) array of rank 2, its shapes is (k,n)
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * error -- a real(type ki), the maximum error obtained to fulfill the 
    !             Moor-Penrose conditions
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    ! cette routine verifie les differentes conditions de Moor-Penrose, 
    ! elle retourne une erreur qui est le maximum des nombres trouves 
    ! au lieu de 0
    !
    subroutine verif_r(n,k,ga,ga_prime,error)
      !
      integer, intent(in) :: n,k
      real(ki), intent(in), dimension(n,k) :: ga
      real(ki), intent(in), dimension(k,n) :: ga_prime
      real(ki), intent (out) :: error
      !
      real(ki), dimension(n,n) :: mat1
      real(ki), dimension(k,k) :: mat2
      real(ki), dimension(n,k) :: mat3
      real(ki), dimension(k,n) :: mat4
      real(ki) :: max1,max2,max3,max4
      !
      mat1 = matmul(ga,ga_prime)
      mat2 = matmul(ga_prime,ga)
      mat3 = abs(matmul(mat1,ga)-ga)
      mat4 = abs(matmul(mat2,ga_prime)-ga_prime)
      !
      mat1 = abs(mat1-transpose(mat1))
      mat2 = abs(mat2-transpose(mat2))
      !
      max1 = maxval(mat1)
      max2 = maxval(mat2)
      max3 = maxval(mat3)
      max4 = maxval(mat4)
      !
      error = max(max1,max2,max3,max4)
      !
    end subroutine verif_r

    subroutine verif_c(n,k,ga,ga_prime,error)
      !
      integer, intent(in) :: n,k
      complex(ki), intent(in), dimension(n,k) :: ga
      complex(ki), intent(in), dimension(k,n) :: ga_prime
      real(ki), intent (out) :: error
      !
      complex(ki), dimension(n,n) :: mat1
      complex(ki), dimension(k,k) :: mat2
      real(ki), dimension(n,n) :: mat1a
      real(ki), dimension(k,k) :: mat2a
      real(ki), dimension(n,k) :: mat3
      real(ki), dimension(k,n) :: mat4
      real(ki) :: max1,max2,max3,max4
      !
      mat1 = matmul(ga,ga_prime)
      mat2 = matmul(ga_prime,ga)
      mat3 = abs(matmul(mat1,ga)-ga)
      mat4 = abs(matmul(mat2,ga_prime)-ga_prime)
      !
      mat1a = abs(mat1-transpose(mat1))
      mat2a = abs(mat2-transpose(mat2))
      !
      max1 = maxval(mat1a)
      max2 = maxval(mat2a)
      max3 = maxval(mat3)
      max4 = maxval(mat4)
      !
      error = max(max1,max2,max3,max4)
      !
    end subroutine verif_c
    !
    !****f* src/kinematic/inverse_matrice/imprime_mat
    ! NAME
    !
    !  Subroutine imprime_mat
    !
    ! USAGE
    !
    !  call imprime_mat(mat)
    !
    ! DESCRIPTION
    !
    !  This routine prints a n x n matrix
    !
    ! INPUTS
    !
    !  * mat -- a real/complex (type ki) array of rank 2
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  No value returned
    !
    ! EXAMPLE
    !
    ! WARNING: swapped lines and columns! mat(line, column)
    !
    !*****
    ! cette routine sert a imprimer une matrice carree

    subroutine imprime_mat_p(mat_p)
      type(s_matrix_poly), intent(inout) :: mat_p
      !
      if (associated(mat_p%pt_cmplx)) then
         call imprime_mat(mat_p%pt_cmplx)
      elseif (associated(mat_p%pt_real)) then
         call imprime_mat(mat_p%pt_real)
      end if
      !
    end subroutine imprime_mat_p
    !
    !
    subroutine imprime_mat_r(mat)
      !
      real(ki), intent(in), dimension(:,:) :: mat
      !
      character(len=11*(size(mat,2)-1)+8+14) :: form 
      integer :: i
      integer, dimension(2) :: dim
      !
      dim = shape(mat)
      form = '(1x,"[",'//repeat('E17.10,TR2,', dim(2)-1)//'E17.10,"]",1x)'

      do i=1,dim(1)
        !
        write (*, fmt=form ) mat(i,:)
        !
      end do
      !
    end subroutine imprime_mat_r
    !
    subroutine imprime_mat_c(matc)
      !
      complex(ki), intent(in), dimension(:,:) :: matc
      !
      character(len=32*(size(matc,2)-1)+30) :: form 
      integer :: i
      integer, dimension(2) :: dim
      !
      dim = shape(matc)
      form = ""
      !
      do i=1,dim(2)-1
        !
        form = trim(form)//'"(",e16.10,1x,"I*",e16.10,")",2x'
        !
      end do
      !
      form = trim(form)//'"(",e16.10,1x,"I*",e16.10,")"'
      !
      do i=1,dim(1)
        !
        write (*,'(1x,"[",'//form//',"]")') matc(i,:)
        !
      end do
      !
    end subroutine imprime_mat_c
    !
    !****if* src/kinematic/inverse_matrice/inverse_triangular
    ! NAME
    !
    !  Subroutine inverse_triangular
    !
    ! USAGE
    !
    !  call inverse_triangular(mat,inf_or_sup,inv_mat)
    !
    ! DESCRIPTION
    !
    !  This routine inverses an upper or lower triangular matrix
    !  The program assumes that the matrix is lower triangular, if it is upper
    !  the program works with the transposed matrix
    !
    ! INPUTS
    !
    !  * mat -- a real/complex (type ki) array of rank 2
    !  * inf_or_sup -- a string of 3 characters : inf or sup
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * inv_mat -- a real/complex (type ki) array of rank 2, the inverse of the triangular matrix 
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine inverse_triangular_r(mat_r,inf_or_sup,inv_mat_r)
      !
      real(ki), intent(in), dimension(:,:) :: mat_r
      character(len=3), intent(in) :: inf_or_sup
      real(ki), intent(out), dimension(size(mat_r,1),size(mat_r,1)) :: inv_mat_r
      !
      logical :: inversible,triangular_sup,triangular_inf
      integer :: n_dim,i,j,k
      real(ki) :: somme
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: mat1
      !
      n_dim = size(mat_r,1) ! dimension de la matrice
      inversible = .true.
      triangular_inf = .false.
      triangular_sup = .false.
      inv_mat_r(:,:) = 0._ki
      !
      if (inf_or_sup == 'inf') triangular_inf = .true.
      if (inf_or_sup == 'sup') triangular_sup = .true.
      !
      do i=1,n_dim
        !
        inversible = inversible .and. .not.(equal_real(mat_r(i,i),0._ki))
        !
      end do
      !
      if ( (inversible) .and. (triangular_inf .or. triangular_sup) ) then
        !
        if (triangular_inf) then
          !
          mat1 = mat_r
          !
        else
          !
          mat1 = transpose(mat_r)
          !
        end if
        !
        do i=1,n_dim
          !
          do j=1,i
            !
            if (j == i) then
              !
              inv_mat_r(i,j) = 1._ki/mat1(i,i)
              !
            else
              !
              somme = 0._ki
              !
              do k=1,i-1
                !
                somme = somme + mat1(i,k)*inv_mat_r(k,j)
                !
              end do
              !
              inv_mat_r(i,j) = -somme/mat1(i,i)
              !
            end if
            !
          end do
          !
        end do
        !
        if (triangular_sup) inv_mat_r = transpose(inv_mat_r)
        !
      else
        !
        !~ if (.not.(inversible)) write(*,*) 'matrice pas inversible'
        if (.not.(inversible)) then
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'alerte, internal error'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'In the LU decomposition'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'One triangular matrix is not invertible'
          call catch_exception(1)
        end if
        !~ if (triangular_inf .and. triangular_sup) write(*,*) 'matrice pas diagonal'
        if (triangular_inf .and. triangular_sup) then
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'alerte, internal error'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'In the LU decomposition'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'One matrix is not diagonal'
          call catch_exception(1)
        end if
        !
      end if
      !
    end subroutine inverse_triangular_r
    !
    subroutine inverse_triangular_c(mat_c,inf_or_sup,inv_mat_c)
      !
      complex(ki), intent(in), dimension(:,:) :: mat_c
      character(len=3), intent(in) :: inf_or_sup
      complex(ki), intent(out), dimension(size(mat_c,1),size(mat_c,1)) :: inv_mat_c
      !
      logical :: inversible,triangular_sup,triangular_inf
      integer :: n_dim,i,j,k
      complex(ki) :: somme
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: mat1
      !
      n_dim = size(mat_c,1) ! dimension de la matrice
      inversible = .true.
      triangular_inf = .false.
      triangular_sup = .false.
      inv_mat_c(:,:) = czero
      !
      if (inf_or_sup == 'inf') triangular_inf = .true.
      if (inf_or_sup == 'sup') triangular_sup = .true.
      !
      do i=1,n_dim
        !
        inversible = inversible .and. .not.(equal_real(abs(mat_c(i,i)),0._ki))
        !
      end do
      !
      if ( (inversible) .and. (triangular_inf .or. triangular_sup) ) then
        !
        if (triangular_inf) then
          !
          mat1 = mat_c
          !
        else
          !
          mat1 = transpose(mat_c)
          !
        end if
        !
        do i=1,n_dim
          !
          do j=1,i
            !
            if (j == i) then
              !
              inv_mat_c(i,j) = cmplx(1._ki,0._ki,ki)/mat1(i,i)
              !
            else
              !
              somme = czero
              !
              do k=1,i-1
                !
                somme = somme + mat1(i,k)*inv_mat_c(k,j)
                !
              end do
              !
              inv_mat_c(i,j) = -somme/mat1(i,i)
              !
            end if
            !
          end do
          !
        end do
        !
        if (triangular_sup) inv_mat_c = transpose(inv_mat_c)
        !
      else
        !
        !~ if (.not.(inversible)) write(*,*) 'matrice pas inversible'
        if (.not.(inversible)) then
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'alerte, internal error'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'In the LU decomposition'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'One triangular matrix is not invertible'
          call catch_exception(1)
        end if
        !~ if (triangular_inf .and. triangular_sup) write(*,*) 'matrice pas diagonal'
        if (triangular_inf .and. triangular_sup) then
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'alerte, internal error'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'In the LU decomposition'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'One matrix is not diagonal'
          call catch_exception(1)
        end if
        !
      end if
      !
    end subroutine inverse_triangular_c
    !
    !****if* src/kinematic/inverse_matrice/inverse_diagonal
    ! NAME
    !
    !  Subroutine inverse_diagonal
    !
    ! USAGE
    !
    !  call inverse_diagonal(mat,inv_mat)
    !
    ! DESCRIPTION
    !
    !  This routine inverses a diagonal matrix
    !
    ! INPUTS
    !
    !  * mat -- a real (type ki) array of rank 2
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * inv_mat -- a real (type ki) array of rank 2, the inverse of the triangular matrix 
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine inverse_diagonal(mat,inv_mat)
      !
      real(ki), intent(in), dimension(:,:) :: mat
      real(ki), intent(out), dimension(size(mat,1),size(mat,1)) :: inv_mat
      !
      logical :: inversible
      integer :: n_dim,i
      !
      n_dim = size(mat,1)  ! dimension de la matrice
      inversible = .true.
      inv_mat(:,:) = 0._ki
      !
      do i=1,n_dim
        !
        inversible = inversible .and. (mat(i,i) /= 0._ki)
        !
      end do
      !
      if (inversible) then
        !
        do i=1,n_dim
          !
          inv_mat(i,i) = 1._ki/mat(i,i)
          !
        end do
        !
      else
        !
        !~ write(*,*) 'matrice pas inversible'
          tab_erreur_par(1)%a_imprimer = .true.
          tab_erreur_par(1)%chaine = 'alerte, internal error'
          tab_erreur_par(2)%a_imprimer = .true.
          tab_erreur_par(2)%chaine = 'In the LU decomposition'
          tab_erreur_par(3)%a_imprimer = .true.
          tab_erreur_par(3)%chaine = 'One diagonal matrix is not invertible'
          call catch_exception(1)
        !
      end if
      !
    end subroutine inverse_diagonal
    !
    !
    !****if* src/kinematic/inverse_matrice/lu_decomp
    ! NAME
    !
    !  Subroutine lu_decomp
    !
    ! USAGE
    !
    !  call lu_decomp(mat,p_mat,l_mat,u_mat)
    !
    ! DESCRIPTION
    !
    !  This routine writes a symmetric matrix as P L U where P is a
    !  permutation matrix, L is a lower triangular matrix and  
    !  U is an upper triangular matrix
    !
    ! INPUTS
    !
    !  * mat -- a real (type ki) array of rank 2
    !
    ! SIDE EFFECTS
    !
    !  No side effect
    !
    ! RETURN VALUE
    !
    !  * p_mat -- a real (type ki) array of rank 2, the permutation matrix 
    !  * l_mat -- a real (type ki) array of rank 2, the lower triangular matrix 
    !  * u_mat -- a real (type ki) array of rank 2, the upper triangular matrix 
    !  * succeed -- logical, true if the decomposition succeeded
    !
    ! EXAMPLE
    !
    !
    !
    !*****
    !
    subroutine lu_decomp_r(mat_r,p_mat,l_mat,u_mat, succeeded)
      !
      real(ki), intent(in), dimension(:,:) :: mat_r
      real(ki), intent(out), dimension(size(mat_r,1),size(mat_r,1)) :: l_mat,u_mat,p_mat
      logical, intent(out) :: succeeded
      !
      integer :: n_dim,i,j,k
      integer, dimension(1) :: loc_m
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: mat_b
      real(ki), dimension(size(mat_r,1),size(mat_r,1)) :: temp_mat,id_mat
      real(ki), dimension(size(mat_r,1)) :: temp_vec
      !
      n_dim = size(mat_r,1)  ! dimension de la matrice
      u_mat(:,:) = 0._ki
      id_mat(:,:) = 0._ki
      succeeded = .true.
      !
      do i=1,n_dim
        !
        id_mat(i,i) = 1._ki
        !
      end do
      !
      mat_b = mat_r
      l_mat(:,:) = 0._ki
      p_mat = id_mat
      !
      do k=1,n_dim-1
        !
        ! plus grand element de la colonne k entre les lignes k et n_dim
        !
        loc_m = maxloc(abs(mat_b(k:n_dim,k))) + k - 1 
        !
        ! si ce plus grand element n'est pas sur la ligne k ou bien on ne 
        ! traite pas la derniere ligne, on permute les lignes.
        ! p_mat garde trace de la permutation
        !
        if ( (loc_m(1) /= k) .and.(k /= n_dim) ) then
          !
          temp_vec = mat_b(k,:)
          mat_b(k,:) = mat_b(loc_m(1),:)
          mat_b(loc_m(1),:) = temp_vec
          !
          temp_vec = l_mat(k,:)
          l_mat(k,:) = l_mat(loc_m(1),:)
          l_mat(loc_m(1),:) = temp_vec
          !
          temp_mat = id_mat
          temp_mat(k,:) = id_mat(loc_m(1),:)
          temp_mat(loc_m(1),:) = id_mat(k,:)
          p_mat = matmul(p_mat,temp_mat)
          !
        end if
        !
        l_mat(k,k) = 1._ki
        if (equal_real(mat_b(k,k),0._ki)) then
                l_mat(k+1:n_dim,k) = 0._ki ! matrix not invertible
                succeeded=.false.
        else
                l_mat(k+1:n_dim,k) = mat_b(k+1:n_dim,k)/mat_b(k,k)
        end if
        !
        do i=k+1,n_dim
          !
          mat_b(i,1:k) = 0._ki
          !
          do j=k+1,n_dim
            !
            mat_b(i,j) = mat_b(i,j) - l_mat(i,k)*mat_b(k,j)
            !
          end do
          !
        end do
        !
      end do
      !
      l_mat(n_dim,n_dim) = 1._ki
      u_mat = mat_b
      !
    end subroutine lu_decomp_r
    !
    !
    subroutine lu_decomp_c(mat_c,p_mat,l_mat,u_mat,succeeded)
      !
      complex(ki), intent(in), dimension(:,:) :: mat_c
      complex(ki), intent(out), dimension(size(mat_c,1),size(mat_c,1)) :: l_mat,u_mat,p_mat
      logical, intent(out) :: succeeded
      !
      integer :: n_dim,i,j,k
      integer, dimension(1) :: loc_m
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: mat_b
      complex(ki), dimension(size(mat_c,1),size(mat_c,1)) :: temp_mat,id_mat
      complex(ki), dimension(size(mat_c,1)) :: temp_vec
      !
      n_dim = size(mat_c,1)  ! dimension de la matrice
      u_mat(:,:) = czero
      id_mat(:,:) = czero
      succeeded=.true.
      !
      do i=1,n_dim
        !
        id_mat(i,i) = cmplx(1._ki,0._ki,ki)
        !
      end do
      !
      mat_b = mat_c
      l_mat(:,:) = czero
      p_mat = id_mat
      !
      do k=1,n_dim-1
        !
        ! plus grand element de la colonne k entre les lignes k et n_dim
        !
        loc_m = maxloc(abs(mat_b(k:n_dim,k))) + k - 1 
        !
        ! si ce plus grand element n'est pas sur la ligne k ou bien on ne 
        ! traite pas la derniere ligne, on permute les lignes.
        ! p_mat garde trace de la permutation
        !
        if ( (loc_m(1) /= k) .and.(k /= n_dim) ) then
          !
          temp_vec = mat_b(k,:)
          mat_b(k,:) = mat_b(loc_m(1),:)
          mat_b(loc_m(1),:) = temp_vec
          !
          temp_vec = l_mat(k,:)
          l_mat(k,:) = l_mat(loc_m(1),:)
          l_mat(loc_m(1),:) = temp_vec
          !
          temp_mat = id_mat
          temp_mat(k,:) = id_mat(loc_m(1),:)
          temp_mat(loc_m(1),:) = id_mat(k,:)
          p_mat = matmul(p_mat,temp_mat)
          !
        end if
        !
        l_mat(k,k) = cmplx(1._ki,0._ki,ki)
        if (equal_real(abs(mat_b(k,k)),0._ki)) then
                l_mat(k+1:n_dim,k) = 0._ki ! matrix not invertible
                succeeded=.false.
        else
                l_mat(k+1:n_dim,k) = mat_b(k+1:n_dim,k)/mat_b(k,k)
        end if
        !
        do i=k+1,n_dim
          !
          mat_b(i,1:k) = czero
          !
          do j=k+1,n_dim
            !
            mat_b(i,j) = mat_b(i,j) - l_mat(i,k)*mat_b(k,j)
            !
          end do
          !
        end do
        !
      end do
      !
      l_mat(n_dim,n_dim) = cmplx(1._ki,0._ki,ki)
      u_mat = mat_b
      !
    end subroutine lu_decomp_c
    !
end module inverse_matrice
