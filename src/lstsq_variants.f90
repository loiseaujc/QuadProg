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
   integer :: i, j, info, ncons
   logical :: is_upper_bounded, is_lower_bounded

   !> Get problem's dimensions.
   m = size(A, 1); n = size(A, 2)

   !----- Sanity checks -----

   !> Sanity checks.
   if (m < n) error stop "For m < n, the problem is not strictly convex and cannot be solved with QuadProg."
   if (size(b) /= m) error stop "A and b have inconsistent number of rows."

   is_upper_bounded = .false.
   if (present(ub)) then
      if (size(ub) /= n) error stop "The number of upper bounds is incompatible with the number of variables."
      is_upper_bounded = .true.
   end if

   if (present(lb)) then
      if (size(lb) /= n) error stop "The number of lower bounds is incompatible with the number of variables."
      is_lower_bounded = .true.
   end if

   if (is_lower_bounded .and. is_upper_bounded) then
      do i = 1, n
         if (lb(i) > ub(i)) error stop "Some lower bounds are actually larger than the provided upper bounds."
      end do
   end if

   !----- Sets up the problem -----
   allocate (x(n)); x = 0.0_dp

   !> Constraint matrix.
   if (is_lower_bounded .and. .not. is_upper_bounded) then
      allocate (C(1, n)); C = 1.0_dp; allocate (d, source=lb)
      allocate (icmat(2, n)); icmat(1, :) = 1; icmat(2, :) = [(i, i=1, n)]
      ncons = n
   else if (is_upper_bounded .and. .not. is_lower_bounded) then
      allocate (C(1, n)); C = 1.0_dp; allocate (d, source=ub)
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
   call qr(A, Q, R)

   !> Compute inv(R).
   call dtrtri("u", "n", n, R, n, info)

   !> Solve the corresponding Quadratic Program.
   block
      real(dp) :: qvec(n), y(ncons), obj
      real(dp), allocatable :: work(:)
      integer :: neq, nact, iter(2), iact(ncons), lwork, r_
      ! Prepare variables.
      neq = 0; ncons = n; r_ = min(n, ncons); lwork = 2*n + r_*(r_ + 5)/2 + 2*ncons + 1
      allocate (work(lwork)); work = 0.0_dp
      ! Solve QP.
      info = 1; qvec = matmul(transpose(A), b)
      call qpgen1(R, qvec, n, n, x, y, obj, C, icmat, d, 1, ncons, neq, iact, nact, iter, work, info)
   end block

   end procedure
end submodule
