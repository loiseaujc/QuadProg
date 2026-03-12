submodule(quadprog) lstsq_variants
   use quadprog_constants, only: dp
   use stdlib_linalg, only: eye
   use stdlib_linalg_lapack, only: gemm
   implicit none(type, external)
contains

   !----------------------------------------------
   !-----     Non-Negative Least-Squares     -----
   !----------------------------------------------

   module procedure nnls
   ! Problem's dimensions.
   integer :: m, n
   ! QR factorization.
   real(dp), allocatable    :: Q(:, :), R(:, :)
   ! Constraint matrix.
   real(dp), allocatable    :: C(:, :), d(:)
   integer, allocatable     :: icmat(:, :)

   ! Miscellaneous.
   integer :: i, info

   !> Get problem's dimensions.
   m = size(A, 1); n = size(A, 2)

   !> Sanity checks.
   call assert(assertion=m > n, &
               description="m < n : Problem is not strictly convex.")
   call assert(assertion=size(b) == m, &
               description="A and b have inconsistent number of rows.")

   !> Allocate matrices.
   allocate (x(n), source=0.0_dp)

   !> Constraint matrix.
   allocate (C(1, n), source=1.0_dp)
   allocate (d(n), source=0.0_dp)
   allocate (icmat(2, n)); icmat(1, :) = 1; icmat(2, :) = [(i, i=1, n)]

   !> QR decomposition of the data matrix.
   allocate (Q(m, n), source=0.0_dp)
   allocate (R(n, n), source=0.0_dp)
   call qr(A, Q, R)

   !> Compute inv(R).
   call trtri("u", "n", n, R, n, info)

   !> Solve the corresponding Quadratic Program.
   block
      real(dp) :: qvec(n), y(n), obj
      real(dp), allocatable :: work(:)
      integer :: neq, ncons, nact, iter(2), iact(n), lwork
      ! Prepare variables.
      neq = 0; ncons = n; lwork = 2*n + n*(n + 5)/2 + 2*n + 1
      allocate (work(lwork), source=0.0_dp)
      ! Solve QP.
      info = 1; qvec = matmul(transpose(A), b)
      call qpgen1(R, qvec, n, n, x, y, obj, C, icmat, d, 1, ncons, &
                  neq, iact, nact, iter, work, info)
   end block
   end procedure nnls

   !---------------------------------------------------
   !-----     Bounded-Variables Least-Squares     -----
   !---------------------------------------------------

   module procedure bvls
   ! Problem's dimensions.
   integer :: m, n
   ! QR factorization.
   real(dp), allocatable :: Q(:, :), R(:, :)
   ! Constraint matrix.
   real(dp), allocatable :: C(:, :), d(:)
   integer, allocatable  :: icmat(:, :)

   ! Miscellaneous.
   integer :: i, info, ncons
   logical :: is_upper_bounded, is_lower_bounded

   !> Get problem's dimensions.
   m = size(A, 1); n = size(A, 2)

   !----- Sanity checks -----
   call assert(assertion=m > n, &
               description=" m < n : Problem is not strictly convex.")
   call assert(assertion=size(b) == m, &
               description="A and b have inconsistent number of rows.")

   is_upper_bounded = .false.
   if (present(ub)) then
      call assert(assertion=size(ub) == n, &
                  description="The number of upper bounds is incompatible with the number of variables.")
      is_upper_bounded = .true.
   end if

   if (present(lb)) then
      call assert(assertion=size(lb) == n, &
                  description="The number of lower bounds is incompatible with the number of variables.")
      is_lower_bounded = .true.
   end if

   if (is_lower_bounded .and. is_upper_bounded) then
      call assert(assertion=all(lb <= ub), &
                  description="Some lower bounds are larger than upper ones. Problem is infeasible.")
   end if

   !----- Sets up the problem -----
   allocate (x(n), source=0.0_dp)

   !> Constraint matrix.
   if (is_lower_bounded .and. .not. is_upper_bounded) then
      allocate (C(1, n), source=1.0_dp); allocate (d, source=lb)
      allocate (icmat(2, n)); icmat(1, :) = 1; icmat(2, :) = [(i, i=1, n)]
      ncons = n
   else if (is_upper_bounded .and. .not. is_lower_bounded) then
      allocate (C(1, n), source=1.0_dp); allocate (d, source=ub)
      allocate (icmat(2, n)); icmat(1, :) = 1; icmat(2, :) = [(i, i=1, n)]
      ncons = n
   else if (is_upper_bounded .and. is_lower_bounded) then
      allocate (C(1, 2*n), d(2*n), icmat(2, 2*n))
      icmat(1, :) = 1
      do i = 1, n
         icmat(2, i) = i; icmat(2, i + n) = i
         C(1, i) = 1.0_dp; C(1, i + n) = -1.0_dp
         d(i) = lb(i); d(i + n) = ub(i)
      end do
      ncons = 2*n
   else
      allocate (C(1, 1), icmat(1, 1), d(1))
      ncons = 0
   end if

   !> QR decomposition of the data matrix.
   allocate (Q(m, n), source=0.0_dp)
   allocate (R(n, n), source=0.0_dp)
   call qr(A, Q, R)

   !> Compute inv(R).
   call trtri("u", "n", n, R, n, info)

   !> Solve the corresponding Quadratic Program.
   block
      real(dp) :: qvec(n), y(ncons), obj
      real(dp), allocatable :: work(:)
      integer :: neq, nact, iter(2), iact(ncons), lwork, r_
      ! Prepare variables.
      neq = 0; ncons = n; r_ = min(n, ncons); lwork = 2*n + r_*(r_ + 5)/2 + 2*ncons + 1
      allocate (work(lwork), source=0.0_dp)
      ! Solve QP.
      info = 1; qvec = matmul(transpose(A), b)
      call qpgen1(R, qvec, n, n, x, y, obj, C, icmat, d, 1, ncons, &
                  neq, iact, nact, iter, work, info)
   end block

   end procedure bvls

   !------------------------------------------------------------------------
   !-----  Least Absolute Shrinkage and Selection Operator (LASSO)     -----
   !------------------------------------------------------------------------

   module procedure lasso
   !> Optional arguments.
   real(dp) :: tol_, rho_
   integer :: maxiter_
   !> Miscellaneous.
   integer  :: m, n, info
   real(dp), allocatable :: Atb(:), z(:), v(:)
   real(dp), allocatable :: L(:, :)

   !> Problem's dimensions.
   m = size(A, 1); n = size(A, 2)

   !> Optional arguments.
   maxiter_ = optval(maxiter, 100)
   tol_ = optval(tol, sqrt(epsilon(1.0_dp)))
   rho_ = optval(rho, 2.0_dp)

   !> Sanity checks.
   call assert(assertion=m >= n, &
               description="A needs to be a tall matrix (m >= n).")
   call assert(assertion=size(b) == m, &
               description="A and b have inconsistent number of rows.")
   call assert(assertion=lambda > 0, &
               description="lambda needs to be strictly positive.")
   call assert(assertion=rho_ > 0, &
               description="rho needs to be strictly positive.")
   call assert(assertion=maxiter_ > 0, &
               description="Maximum number of iterations needs to be positive.")
   call assert(assertion=tol_ > 0, &
               description="Tolerance needs to be strictly positive.")

   !> Allocate arrays.
   allocate (x(n), z(n), v(n), source=0.0_dp)
   !> Cholesky factorization of (rho*I + A^T @ A).
   L = rho_*eye(n) + matmul(transpose(A), A)
   call cholesky(L)
   end procedure lasso
end submodule lstsq_variants
