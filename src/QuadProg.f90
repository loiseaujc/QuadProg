module QuadProg
   use quadprog_constants, only: dp
   implicit none
   private

   public :: dp
   public :: solve
   public :: nnls, bvls
   public :: qpgen1, qpgen2

   !---------------------------------
   !-----     DERIVED-TYPES     -----
   !---------------------------------

   type, public :: OptimizeResult
      real(dp), allocatable :: x(:)
      !! Solution of the optimization problem.
      real(dp), allocatable :: y(:)
      !! Lagrange multipliers associated to each constraint.
      real(dp)              :: obj
      !! Objective function evaluated at the minimizer.
      logical               :: success
      !! Whether the problem has been successfully solved or not.
   end type OptimizeResult

   type, public :: qp_problem
      real(dp), allocatable :: P(:, :), q(:)
      !! Symmetric positive definite matrix and associated vector defining the
      !! quadratic form 1/2 x.T @ P @ x - x.T @ q.
      real(dp), allocatable :: A(:, :), b(:)
      !! Matrix and vector defining the lineear equality constraints A @ x = b.
      real(dp), allocatable :: C(:, :), d(:)
      !! Matrix and vector defining the inequality constraints C @ x >= d.
      integer               :: neq, ncons
   end type qp_problem
   interface qp_problem
      module procedure initialize_qp_problem
   end interface

   type, public :: compact_qp_problem
      real(dp), allocatable :: P(:, :), q(:)
      !! Symmetric positive definite matrix and associated vector defining the
      !! quadratic form 1/2 x.T @ P @ x - x.T @ q.
      real(dp), allocatable :: A(:, :), b(:)
      integer, allocatable :: iamat(:, :)
      !! Matrix and vector defining the linear equality constraints A @ x = b using
      !! its compact representation.
      real(dp), allocatable :: C(:, :), d(:)
      integer, allocatable :: icmat(:, :)
      !! Matrix and vector defining the inequality constraints C @ x >= d using
      !! its compact representation.
      integer               :: neq, ncons
   end type compact_qp_problem
   interface compact_qp_problem
      module procedure initialize_compact_qp_problem
   end interface

   !---------------------------------------------------
   !-----     INTERFACE FOR SOLVE(QP_PROBLEM)     -----
   !---------------------------------------------------

   interface solve
      module procedure solve_standard_qp
      module procedure solve_compact_qp
   end interface

   !------------------------------------------------------------------
   !-----     INTERFACES FOR THE QUADPROG MODERNIZED DRIVERS     -----
   !------------------------------------------------------------------

   interface
      module subroutine qpgen1(dmat, dvec, fddmat, n, sol, lagr, crval, amat, iamat, bvec, &
                               fdamat, q, meq, iact, nact, iter, work, ierr)
         implicit none
         integer, intent(in)     :: fddmat, n
         !! Dimensions of the symmetric positive definit matrix Dmat.
         integer, intent(in)     :: fdamat, q
         !! Dimensions of the constraints matrix Amat.
         integer, intent(in)     :: meq
         !! Number of equality constraints.
         integer, intent(out)    :: iact(*), nact
         !! Indices and number of active constraints at the optimum.
         integer, intent(out)    :: iter(*)
         !! Number of iterations.
         integer, intent(inout)  :: ierr
         !! Information flag.
         real(dp), intent(inout) :: dmat(fddmat, *), dvec(*)
         !! Sym. pos. def. matrix and vector defining the quadratic cost.
         real(dp), intent(out)   :: lagr(*), sol(*)
         !! Lagrange multipliers and solution vector.
         integer, intent(in)     :: iamat(fdamat + 1, *)
         real(dp), intent(inout) :: amat(fdamat, *), bvec(*)
         !! Matrix and vector defining the (in-)equality constraints.
         real(dp), intent(inout) :: work(*)
         !! Workspace.
         real(dp), intent(out)   :: crval
         !! Cost function at the optimum.
      end subroutine qpgen1

      module subroutine qpgen2(dmat, dvec, fddmat, n, sol, lagr, crval, amat, bvec, fdamat, q, &
                               meq, iact, nact, iter, work, ierr)
         implicit none
         integer, intent(in)     :: fddmat, n
         !! Dimensions of the symmetric positive definite matrix Dmat.
         integer, intent(in)     :: fdamat, q
         !! Dimensions of the constraint matrix Amat
         integer, intent(in)     :: meq
         !! Number of equality constraints.
         integer, intent(out)    :: iact(*), nact
         !! Indices and number of active constraints at the optimum.
         integer, intent(out)    :: iter(*)
         !! Number of iterations.
         integer, intent(inout)  :: ierr
         !! Information flag.
         real(dp), intent(inout) :: dmat(fddmat, *), dvec(*)
         !! Sym. pos. def. matrix and vector defining the quadratic cost.
         real(dp), intent(out)   :: lagr(*), sol(*)
         !! Lagrange multipliers and solution vector.
         real(dp), intent(inout) :: amat(fdamat, *), bvec(*)
         !! Matrix and vector defining the (in-)equality constraints.
         real(dp), intent(inout) :: work(*)
         !! Workspace.
         real(dp), intent(out)   :: crval
         !! Cost function at the optimum.
      end subroutine qpgen2
   end interface

   !------------------------------------------------------------------------
   !-----     INTERFACES FOR THE VARIANTS OF LEAST-SQUARES SOLVERS     -----
   !------------------------------------------------------------------------

   interface
      module function nnls(A, b) result(x)
         implicit none
         real(dp), intent(in)           :: A(:, :)
         real(dp), intent(in)           :: b(:)
         real(dp), allocatable          :: x(:)
      end function nnls

      module function bvls(A, b, ub, lb) result(x)
         implicit none
         real(dp), intent(in)           :: A(:, :)
         real(dp), intent(in)           :: b(:)
         real(dp), optional, intent(in) :: ub(:)
         real(dp), optional, intent(in) :: lb(:)
         real(dp), allocatable          :: x(:)
      end function bvls
   end interface

   !---------------------------------------------------------------
   !-----     INTERFACES FOR THE LINALG UTILITY FUNCTIONS     -----
   !---------------------------------------------------------------

   interface
      module subroutine qr(A, Q, R)
         implicit none
         real(dp), intent(in) :: A(:, :)
         real(dp), allocatable, intent(out) :: Q(:, :), R(:, :)
      end subroutine qr
   end interface

contains

   !---------------------------------------
   !-----                             -----
   !-----     STANDARD QP PROBLEM     -----
   !-----                             -----
   !---------------------------------------

   module type(qp_problem) function initialize_qp_problem(P, q, A, b, C, d) result(prob)
      implicit none
      real(dp), intent(in)           :: P(:, :), q(:)
      real(dp), optional, intent(in) :: A(:, :), b(:)
      real(dp), optional, intent(in) :: C(:, :), d(:)
      integer :: info, n

      prob%neq = 0; prob%ncons = 0

      !> Sanity checks for the quadratic form.
      if (size(P, 1) /= size(P, 2)) error stop "Matrix P is not square."
      if (size(P, 1) /= size(q)) error stop "Matrix P and vector q have incompatible dimensions."

      !> Quadratic cost.
      prob%P = P; prob%q = q; n = size(P, 1)

      !> Pre-factorize the symmetric positive definite matrix.
      call dpotrf("u", n, prob%P, n, info)
      call dtrtri("u", "n", n, prob%P, n, info)

      !> Sanity checks for the equality constraints.
      if (present(A) .and. .not. present(b)) error stop "Right-hand side vector b for the equality constraints is missing."
      if (.not. present(A) .and. present(b)) error stop "Matrix A for the equality constraints is missing."
      if (present(A) .and. present(b)) then
         if (size(P, 2) /= size(A, 2)) error stop "Matrices P and A have incompatible number of columns."
         if (size(A, 1) /= size(b)) error stop "Matrix A and vector b have incompatible dimensions."
         prob%A = A; prob%b = b; prob%neq = size(b); prob%ncons = size(b)
      end if

      !> Sanity checks for the inequality constraints.
      if (present(C) .and. .not. present(d)) error stop "Right-hand side vector d for the inequality constraints is missing."
      if (.not. present(C) .and. present(d)) error stop "Matrix C for the inequality constraints is missing."
      if (present(C) .and. present(d)) then
         if (size(P, 2) /= size(C, 2)) error stop "Matrices P and C have incompatible number of columns."
         if (size(C, 1) /= size(d)) error stop "Matrix C and vector d have incompatible dimensions."
         prob%C = C; prob%d = d; prob%ncons = prob%neq + size(d)
      end if

      return
   end function initialize_qp_problem

   module subroutine get_constraints_matrix(prob, G, h)
      implicit none
      type(qp_problem), intent(in) :: prob
      !! Quadratic Problem to be solved.
      real(dp), allocatable, intent(out) :: G(:, :)
      !! Constraints matrix expected by qpgen2.
      real(dp), allocatable, intent(out) :: h(:)
      !! Constraints vector expected by qpgen2.

      ! Internal variables.
      logical :: is_constrained
      integer :: i

      associate (n => size(prob%P, 1), neq => prob%neq, ncons => prob%ncons)
         is_constrained = allocated(prob%A) .or. allocated(prob%C)
         if (is_constrained) then
            allocate (G(n, ncons), h(ncons), source=0.0_dp)
            !> Linear equality constraints.
            if (allocated(prob%A)) then
               do i = 1, neq
                  G(:, i) = prob%A(i, :); h(i) = prob%b(i)
               end do
            end if
            !> Linear inequality constraints.
            if (allocated(prob%C)) then
               do i = neq + 1, ncons
                  G(:, i) = prob%C(i - neq, :); h(i) = prob%d(i - neq)
               end do
            end if
         else
            allocate (G(1, 1), h(1), source=0.0_dp)
         end if
      end associate

      return
   end subroutine get_constraints_matrix

   module type(OptimizeResult) function solve_standard_qp(problem, legacy) result(result)
      implicit none
      type(qp_problem), intent(in) :: problem
      logical, optional, intent(in) :: legacy
      logical :: legacy_
      real(dp), allocatable :: P(:, :), q(:)
      real(dp), allocatable :: G(:, :), h(:)
      real(dp), allocatable :: work(:)
      integer               :: n, neq, ncons, r, lwork, nact, iter(2), info
      integer, allocatable  :: iact(:)

      n = size(problem%P, 1); neq = problem%neq; ncons = problem%ncons
      legacy_ = .false.; if (present(legacy)) legacy_ = legacy
      !> Allocate data.
      allocate (iact(ncons))
      allocate (P, source=problem%P); allocate (q, source=problem%q)
      allocate (result%x, mold=q); result%x = 0.0_dp
      allocate (result%y(ncons), source=0.0_dp)
      !> Allocate workspace
      r = min(n, ncons); lwork = 2*n + r*(r + 5)/2 + 2*ncons + 1
      allocate (work(lwork), source=0.0_dp)
      !> Get the constraints matrix and vector.
      call get_constraints_matrix(problem, G, h)
      !> Solve the QP problem.
      info = 1 ! P is already factorized when defining the QP.
      if (legacy_) then
         call legacy_qpgen2(P, q, n, n, result%x, result%y, result%obj, G, h, n, &
                            ncons, neq, iact, nact, iter, work, info)
      else
         call qpgen2(P, q, n, n, result%x, result%y, result%obj, G, h, n, &
                     ncons, neq, iact, nact, iter, work, info)
      end if
      !> Success?
      result%success = (info == 0)
      return
   end function solve_standard_qp

   !--------------------------------------
   !-----                            -----
   !-----     COMPACT QP PROBLEM     -----
   !-----                            -----
   !--------------------------------------

   module type(compact_qp_problem) function initialize_compact_qp_problem(P, q, A, iamat, b, &
                                                                          C, icmat, d) result(prob)
      implicit none
      real(dp), intent(in)           :: P(:, :), q(:)
      real(dp), optional, intent(in) :: A(:, :), b(:)
      integer, optional, intent(in)  :: iamat(:, :)
      real(dp), optional, intent(in) :: C(:, :), d(:)
      integer, optional, intent(in)  :: icmat(:, :)
      integer :: info, n

      prob%neq = 0; prob%ncons = 0

      !> Sanity checks for the quadratic form.
      if (size(P, 1) /= size(P, 2)) error stop "Matrix P is not square."
      if (size(P, 1) /= size(q)) error stop "Matrix P and vector q have incompatible dimensions."

      !> Quadratic cost.
      prob%P = P; prob%q = q; n = size(P, 1)

      !> Pre-factorize the symmetric positive definite matrix.
      call dpotrf("u", n, prob%P, n, info)
      call dtrtri("u", "n", n, prob%P, n, info)

      !> Sanity checks for the equality constraints.
      if (present(A) .and. .not. present(iamat)) error stop "Matrix A is provided but not iamat."
      if (.not. present(A) .and. present(iamat)) error stop "iamat is provided but not matrix A."
      if (present(A) .and. .not. present(b)) error stop "Right-hand side vector b for the equality constraints is missing."
      if (.not. present(A) .and. present(b)) error stop "Matrix A for the equality constraints is missing."

      if (present(A) .and. present(b) .and. present(iamat)) then
         if (size(iamat, 1) /= size(A, 1) + 1) error stop "Matrix A and index iamat have incompatible dimensions."
         if (size(iamat, 2) /= size(A, 2)) error stop "Matrix A and index iamat have incompatible dimensions."
         prob%A = A; prob%iamat = iamat; prob%b = b; prob%neq = size(b); prob%ncons = size(b)
      end if

      !> Sanity checks for the inequality constraints.
      if (present(C) .and. .not. present(icmat)) error stop "Matrix C is provided but not icmat."
      if (.not. present(C) .and. present(icmat)) error stop "icmat is provided but not matrix C."
      if (present(C) .and. .not. present(d)) error stop "Right-hand side vector d for the inequality constraints is missing."
      if (.not. present(C) .and. present(d)) error stop "Matrix C for the inequality constraints is missing."

      if (present(C) .and. present(d) .and. present(icmat)) then
         if (size(icmat, 1) /= size(C, 1) + 1) error stop "Matrix C and index icmat have incompatible dimensions."
         if (size(icmat, 2) /= size(C, 2)) error stop "Matrix C and index icmat have incompatible dimensions."
         prob%C = C; prob%icmat = icmat; prob%d = d; prob%ncons = prob%neq + size(d)
      end if
   end function initialize_compact_qp_problem

   module subroutine get_compact_constraints_matrix(prob, G, igmat, h)
      implicit none
      type(compact_qp_problem), intent(in) :: prob
      !! Quadratic Problem to be solved.
      real(dp), allocatable, intent(out)   :: G(:, :)
      integer, allocatable, intent(out)    :: igmat(:, :)
      !! Constraints matrix expected by qpgen1 (compact representation).
      real(dp), allocatable, intent(out)   :: h(:)
      !! Constraints vector expected by qpgen1.

      ! Internal variables.
      logical :: is_constrained
      integer :: ma, na, mc, nc

      associate (n => size(prob%P, 1), neq => prob%neq, ncons => prob%ncons)
         is_constrained = allocated(prob%A) .or. allocated(prob%C)
         if (is_constrained) then
            !> Linear equality constraints.
            if (allocated(prob%A)) then
               G = prob%A; igmat = prob%iamat; h = prob%b
            end if
            !> Linear inequality constraints.
            if (allocated(prob%C)) then
               if (allocated(prob%A)) then
                  !> Both equality and inequality constraints.
                  ma = size(prob%A, 1); na = size(prob%A, 2)
                  mc = size(prob%C, 1)

                  allocate (G(ma + mc, n))
                  allocate (igmat(ma + mc + 1, n))

                  !> Construct the constraint matrix.
                  G(:na, :) = prob%A
                  G(na + 1:, :) = prob%C

                  igmat(1, :) = prob%iamat(1, :) + prob%icmat(1, :)
                  igmat(2:na + 1, :) = prob%iamat
                  igmat(na + 2:, :) = prob%icmat
               else
                  !> Only inequality constraints.
                  G = prob%C; igmat = prob%icmat; h = prob%d
               end if
            end if
         else
            allocate (G(1, n), h(1), source=0.0_dp)
            allocate (igmat(2, n), source=0)
         end if
      end associate

      return
   end subroutine get_compact_constraints_matrix

   module type(OptimizeResult) function solve_compact_qp(problem) result(result)
      implicit none
      type(compact_qp_problem), intent(in) :: problem
      real(dp), allocatable :: P(:, :), q(:)
      real(dp), allocatable :: G(:, :), h(:)
      integer, allocatable  :: igmat(:, :)
      real(dp), allocatable :: work(:)
      integer               :: n, neq, ncons, r, lwork, nact, iter(2), info
      integer, allocatable  :: iact(:)

      n = size(problem%P, 1); neq = problem%neq; ncons = problem%ncons

      !> Allocate data.
      allocate (iact(ncons))
      allocate (P, source=problem%P); allocate (q, source=problem%q)
      allocate (result%x, mold=q); result%x = 0.0_dp
      allocate (result%y(ncons), source=0.0_dp)

      !> Allocate workspace
      r = min(n, ncons); lwork = 2*n + r*(r + 5)/2 + 2*ncons + 1
      allocate (work(lwork), source=0.0_dp)

      !> Get the constraints matrix and vector.
      call get_compact_constraints_matrix(problem, G, igmat, h)

      !> Solve the QP problem.
      info = 1 ! P is already factorized when defining the QP.
      call qpgen1(P, q, n, n, result%x, result%y, result%obj, G, igmat, h, size(G, 1), &
                  ncons, neq, iact, nact, iter, work, info)
      !> Success?
      result%success = (info == 0)
      return
   end function solve_compact_qp

end module QuadProg
