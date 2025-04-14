submodule(quadprog) lstsq_variants
   use quadprog_constants, only: dp
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
   integer :: i, j

   !> Get problem's dimensions.
   m = size(A, 1); n = size(A, 2)

   !> Sanity checks.
   if (m < n) error stop "For m < n, the problem is not strictly convex and cannot be solved with QuadProg."
   if (size(b) /= m) error stop "A and b have inconsistent number of rows."
   if (present(lambda)) then
      if (lambda < 0) error stop "Regularization parameter lambda needs to be positive."
      error stop "Regularization not yet supported."
   end if

   !> Allocate matrices.
   Q = A; allocate (R(n, n)); R = 0.0_dp
   allocate (x(n)); x = 0.0_dp

   !> QR decomposition of the data matrix.
   block
      real(dp) :: tau(n)
      real(dp), allocatable :: work(:)
      integer :: lwork, info
      !----- Factorization -----
      !> Workspace query
      lwork = -1; allocate (work(1)); call dgeqrf(m, n, Q, m, tau, work, lwork, info)
      !> Factorization.
      lwork = work(1); deallocate (work); allocate (work(lwork)); work = 0.0_dp
      call dgeqrf(m, n, Q, m, tau, work, lwork, info)

      !> Extract the R matrix
      do concurrent(i=1:n, j=1:n)
         if (j >= i) R(i, j) = Q(i, j)
      end do

      !> Compute inv(R).
      call dtrtri("u", "n", n, R, n, info)

      !----- Extract the Q matrix -----
      !> Workspace query.
      lwork = -1; call dorgqr(m, n, n, Q, m, tau, work, lwork, info)
      !> Compute Q.
      lwork = work(1); deallocate (work); allocate (work(lwork)); work = 0.0_dp
      call dorgqr(m, n, n, Q, m, tau, work, lwork, info)
   end block

   !> Constraint matrix.
   allocate (C(1, n), icmat(2, n), d(n)); C = 1.0_dp; d = 0.0_dp
   icmat(1, :) = 1; icmat(2, :) = [(i, i=1, n)]

   !> Solve the corresponding Quadratic Program.
   block
      real(dp) :: qvec(n), y(n), obj
      real(dp), allocatable :: work(:)
      integer :: neq, ncons, nact, iter(2), iact(n), info, lwork
      ! Prepare variables.
      neq = 0; ncons = n; lwork = 2*n + n*(n + 5)/2 + 2*n + 1
      allocate (work(lwork)); work = 0.0_dp
      ! Solve QP.
      info = 1; qvec = matmul(transpose(A), b)
      call qpgen1(R, qvec, n, n, x, y, obj, C, icmat, d, 1, ncons, neq, iact, nact, iter, work, info)
   end block
   end procedure

end submodule
