module QuadProg
   use quadprog_constants, only: dp
   implicit none
   private

   public :: solve_qp
   public :: qpgen2

   !----- Interface to the legacy (yet modernized) Fortran code for solving a QP problem -----

   interface
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

   !----- Stdlib-like interface for solving a QP problem -----

   interface
    !!  #### Description
    !!
    !!  ...
    !!
    !!  #### Syntax
    !!
    !!  - To solve an unconstrained quadratic program:
    !!
    !!  ```fortran
    !!  call solve_qp(P, q, x)
    !!  ```
    !!
    !!  #### Arguments
    !!
      module function solve_qp(P, q, Aeq, beq, C, d, y, obj, &
                               factorized, overwrite_p, info) result(x)
         real(dp), intent(inout), contiguous, target  :: P(:, :)
      !! n x n Symmetric positive-definite matrix defining the quadratic form.
         real(dp), intent(in)                         :: q(:)
      !! Vector defining the linear term of the quadratic cost.
         real(dp), optional, intent(in)               :: Aeq(:, :)
      !! Matrix defining the linear equality constraints Aeq @ x = beq.
         real(dp), optional, intent(in)               :: beq(:)
      !! Right-hand side vector of the linear equality constraints Aeq @ x = beq.
         real(dp), optional, intent(in)               :: C(:, :)
      !! Matrix defining the linear inequality constraints C @ x >= d.
         real(dp), optional, intent(in)               :: d(:)
      !! Right-hand side vector of the linear inequality constraints C @ x >= d.
         real(dp), optional, allocatable, target, intent(out) :: y(:)
      !! Vector of Lagrange multipliers at the optimum.
         real(dp), optional, intent(out), target      :: obj
      !! Cost function at the optimum.
         logical, optional, intent(inout)            :: factorized
      !! Whether P has already been factorized using ??? (default .false.)
         logical, optional, intent(in)               :: overwrite_p
      !! Whether P can be overwritten (default .false.)
         integer, intent(out)              :: info
      !! Information flag.

         real(dp), allocatable :: x(:)
      !! Vector corresponding to the minimizer.
      end function
   end interface

contains

   pure elemental function optval(x, default) result(y)
      logical, intent(in), optional :: x
      logical, intent(in)           :: default
      logical                       :: y
      if (present(x)) then
         y = x
      else
         y = default
      end if
   end function

   module procedure solve_qp
   integer                       :: i, n, ncons, neq, r, lwork
   logical                       :: is_factorized, can_overwrite, is_constrained
   real(dp), target              :: obj_
   real(dp), allocatable         :: work(:), q_(:), G(:, :), h(:)
   real(dp), allocatable, target :: P_(:, :)
   real(dp), pointer             :: Pmat(:, :)
   real(dp), allocatable, target :: y_(:)
   real(dp), pointer             :: yp(:)

   !> Sanity checks.
   if (size(P, 1) /= size(P, 2)) error stop "Matrix P is not square."
   if (size(P, 1) /= size(q)) error stop "Matrix P and vector q have incompatible dimensions."

   if (present(Aeq) .and. .not. present(beq)) error stop "Right-hand side beq vector for the equality constraints is missing."
   if (.not. present(Aeq) .and. present(beq)) error stop "Matrix Aeq for the equality constraints is missing."
   if (present(Aeq) .and. present(beq)) then
      if (size(P, 2) /= size(Aeq, 2)) error stop "Matrices P and Aeq have incompatible number of columns."
      if (size(Aeq, 1) /= size(beq)) error stop "Matrix Aeq and vector beq have incompatible dimensions."
   end if

   if (present(C) .and. .not. present(d)) error stop "Right-hand side d vector for the inequality constraints is missing."
   if (.not. present(C) .and. present(d)) error stop "Matrix C for the inequality constraints is missing."
   if (present(C) .and. present(d)) then
      if (size(C, 2) /= size(P, 2)) error stop "Matrices P and C have incompatible number of columns."
      if (size(C, 1) /= size(d)) error stop "Matrix C and vector d have incompatible dimensions."
   end if

   !> Sets up problem's dimensions.
   n = size(P, 1)              ! Dimension of the problem.
   if (present(Aeq)) then      ! Number of linear equality constraints.
      neq = size(Aeq, 1)
   else
      neq = 0
   end if
   if (present(C)) then        ! Total number of constraints.
      ncons = neq + size(C, 1)
   else
      ncons = neq
   end if

   !> Allocate workspace and solution vector.
   allocate (x, mold=q); x = 0.0_dp
   r = min(n, ncons); lwork = 2*n + r*(r + 5)/2 + 2*ncons + 1
   allocate (work(lwork)); work = 0.0_dp

   !> Optional boolean arguments.
   is_factorized = optval(factorized, .false.)
   can_overwrite = optval(overwrite_p, .false.)

   !> Setup the problem's data.
   if (can_overwrite) then
      Pmat(1:n, 1:n) => P
   else
      P_ = P; Pmat(1:n, 1:n) => P_
   end if
   q_ = q

   !> Setup the matrix of constraints.
   is_constrained = present(Aeq) .or. present(C)
   if (is_constrained) then
      allocate (G(n, ncons)); G = 0.0_dp
      allocate (h(ncons)); h = 0.0_dp
      !> Linear equality constraints.
      if (present(Aeq)) then
         do i = 1, neq
            G(:, i) = Aeq(i, :); h(i) = beq(i)
         end do
      end if

      !> Linear inequality constraints.
      if (present(C)) then
         do i = neq + 1, ncons
            G(:, i) = C(i, :); h(i) = d(i)
         end do
      end if
   else
      allocate (G(1, 1)); G = 0.0_dp
      allocate (h(1)); h = 0.0_dp
   end if

   !> Setup the vector of Lagrange multipliers.
   if (present(y)) then
      if (allocated(y) .and. (size(y) < ncons)) then
         error stop "Dimension of vector y is too small."
      else
         allocate (y(max(1, ncons))); y = 0.0_dp
      end if
      yp(1:max(1, ncons)) => y
   else
      allocate (y_(max(1, ncons))); y_ = 0.0_dp
      yp(1:max(1, ncons)) => y_
   end if

   !> Solve the QP problem.
   block
      integer           :: iact(2), nact, iter(2)
      real(dp), pointer :: crval
      if (present(obj)) then
         crval => obj
      else
         crval => obj_
      end if
      info = merge(1, 0, is_factorized)
      call qpgen2(Pmat, q_, n, n, x, yp, crval, &
                  G, h, n, ncons, neq, iact, nact, iter, work, info)
   end block
   return
   end procedure

end module QuadProg
