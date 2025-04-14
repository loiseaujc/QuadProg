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
   integer :: i, j, info

   !> Get problem's dimensions.
   m = size(A, 1); n = size(A, 2)

   !> Sanity checks.
   if (m < n) error stop "For m < n, the problem is not strictly convex and cannot be solved with QuadProg."
   if (size(b) /= m) error stop "A and b have inconsistent number of rows."

   !> Allocate matrices.
   allocate (x(n)); x = 0.0_dp

   !> Constraint matrix.
   allocate (C(1, n), icmat(2, n), d(n)); C = 1.0_dp; d = 0.0_dp
   icmat(1, :) = 1; icmat(2, :) = [(i, i=1, n)]

   !> QR decomposition of the data matrix.
   call qr(A, Q, R)

   !> Compute inv(R).
   call dtrtri("u", "n", n, R, n, info)

   !> Solve the corresponding Quadratic Program.
   block
      real(dp) :: qvec(n), y(n), obj
      real(dp), allocatable :: work(:)
      integer :: neq, ncons, nact, iter(2), iact(n), lwork
      ! Prepare variables.
      neq = 0; ncons = n; lwork = 2*n + n*(n + 5)/2 + 2*n + 1
      allocate (work(lwork)); work = 0.0_dp
      ! Solve QP.
      info = 1; qvec = matmul(transpose(A), b)
      call qpgen1(R, qvec, n, n, x, y, obj, C, icmat, d, 1, ncons, neq, iact, nact, iter, work, info)
   end block
   end procedure

end submodule
