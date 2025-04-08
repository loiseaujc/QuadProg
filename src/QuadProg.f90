module QuadProg
   use quadprog_constants, only: dp
   implicit none
   private

   public :: dp
   public :: solve
   public :: qpgen2

   type, public :: OptimizeResult
      real(dp), allocatable :: x(:)
      !! Solution of the optimization problem.
      real(dp), allocatable :: y(:)
      !! Lagrange multipliers associated to each constraint.
      real(dp)              :: obj
      !! Objective function evaluated at the minimizer.
      logical               :: success
      !! Whether the problem has been successfully solved or not.
   end type

   type, public :: qp_problem
      real(dp), allocatable :: P(:, :), q(:)
      !! Symmetric positive definite matrix and associated vector defining the
      !! quadratic form 1/2 x.T @ P @ x - x.T @ q.
      real(dp), allocatable :: A(:, :), b(:)
      !! Matrix and vector defining the lineear equality constraints A @ x = b.
      real(dp), allocatable :: C(:, :), d(:)
      !! Matrix and vector defining the inequality constraints C @ x >= d.
      integer               :: neq, ncons
   end type
   interface qp_problem
      module procedure initialize_qp_problem
   end interface

   interface
      pure module subroutine qpgen1(dmat, dvec, fddmat, n, sol, lagr, crval, amat, iamat, bvec, fdamat, q, &
                                    meq, iact, nact, iter, work, ierr)
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
      end subroutine

      pure module subroutine qpgen2(dmat, dvec, fddmat, n, sol, lagr, crval, amat, bvec, fdamat, q, &
                                    meq, iact, nact, iter, work, ierr)
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
      end subroutine
   end interface

   interface
      pure module subroutine dpofa(A, lda, n, info)
         integer, intent(in) :: lda
         real(dp), intent(out) :: A(lda, *)
         integer, intent(in) :: n
         integer, intent(out) :: info
      end subroutine

      pure module subroutine dpori(A, lda, n)
         integer, intent(in) :: lda
         real(dp), intent(out) :: A(lda, *)
         integer, intent(in) :: n
      end subroutine
   end interface

contains

   type(qp_problem) function initialize_qp_problem(P, q, A, b, C, d) result(prob)
      real(dp), intent(in)           :: P(:, :), q(:)
      real(dp), optional, intent(in) :: A(:, :), b(:)
      real(dp), optional, intent(in) :: C(:, :), d(:)
      integer :: info

      prob%neq = 0; prob%ncons = 0

      !> Sanity checks for the quadratic form.
      if (size(P, 1) /= size(P, 2)) error stop "Matrix P is not square."
      if (size(P, 1) /= size(q)) error stop "Matrix P and vector q have incompatible dimensions."

      !> Quadratic cost.
      prob%P = P; prob%q = q

      !> Pre-factorize the symmetric positive definite matrix.
      call dpofa(prob%P, size(P, 1), size(P, 2), info)
      call dpori(prob%P, size(P, 1), size(P, 2))

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
   end function

   subroutine get_constraints_matrix(prob, G, h)
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
            allocate (G(n, ncons), h(ncons)); G = 0.0_dp; h = 0.0_dp
            !> Linear equality constraints.
            if (allocated(prob%A)) then
               do i = 1, neq
                  G(:, i) = prob%A(i, :); h(i) = prob%b(i)
               end do
            end if
            !> Linear inequality constraints.
            if (allocated(prob%C)) then
               do i = neq + 1, ncons
                  G(:, i) = prob%C(i, :); h(i) = prob%d(i)
               end do
            end if
         else
            allocate (G(1, 1), h(1)); G = 0.0_dp; h = 0.0_dp
         end if
      end associate

      return
   end subroutine

   type(OptimizeResult) function solve(problem) result(result)
      type(qp_problem), intent(in) :: problem
      real(dp), allocatable :: P(:, :), q(:)
      real(dp), allocatable :: G(:, :), h(:)
      real(dp), allocatable :: work(:)
      integer               :: n, neq, ncons, r, lwork, nact, iter(2), info
      integer, allocatable  :: iact(:)

      n = size(problem%P, 1); neq = problem%neq; ncons = problem%ncons
      !> Allocate data.
      allocate (iact(ncons))
      allocate (P, source=problem%P); allocate (q, source=problem%q)
      allocate (result%x, source=q); result%x = 0.0_dp
      allocate (result%y(ncons)); result%y = 0.0_dp
      !> Allocate workspace
      r = min(n, ncons); lwork = 2*n + r*(r + 5)/2 + 2*ncons + 1
      allocate (work(lwork)); work = 0.0_dp
      !> Get the constraints matrix and vector.
      call get_constraints_matrix(problem, G, h)
      !> Solve the QP problem.
      info = 1 ! P is already factorized when defining the QP.
      call qpgen2(P, q, n, n, result%x, result%y, result%obj, G, h, n, ncons, neq, iact, nact, iter, work, info)
      !> Success?
      result%success = (info == 0)
      return
   end function

end module QuadProg
