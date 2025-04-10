module TestProblems
   use iso_fortran_env, only: output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use QuadProg
   implicit none
   private

   real(dp), parameter, private :: atol = 10.0_dp**(-precision(1.0_dp))
   real(dp), parameter, private :: rtol = sqrt(atol)

   public :: collect_test_problems

contains

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   subroutine collect_test_problems(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Unconstrained QP problem", test_problem_1), &
                  new_unittest("(Legacy) Unconstrained QP problem", test_legacy_problem_1), &
                  new_unittest("Inequality Constrained QP problem", test_problem_2), &
                  new_unittest("Equality Constrained QP problem", test_problem_3) &
                  ]
   end subroutine

   !---------------------------------------------
   !-----                                   -----
   !-----     UNCONSTRAINED QR PROBLEMS     -----
   !-----                                   -----
   !---------------------------------------------
   subroutine test_legacy_problem_1(error)
      type(error_type), allocatable, intent(out) :: error
      ! Size of the problem.
      integer, parameter :: n = 3
      ! Quadratic Problem.
      real(dp) :: P(n, n), q(n), obj
      real(dp), allocatable :: x(:), y(:), G(:, :), h(:)
      real(dp), allocatable :: work(:)
      integer               :: i, neq, ncons, r, lwork, nact, iter(2), info
      integer, allocatable  :: iact(:)

      ! Setup the quadratic form.
      P = 0.0_dp; forall (i=1:n) P(i, i) = 1.0_dp; q = 1.0_dp

      ! Setup the constraints.
      neq = 0; ncons = 0
      !> Allocate data.
      allocate (iact(ncons))
      allocate (x, source=q); x = 0.0_dp
      allocate (y(ncons)); y = 0.0_dp
      !> Allocate workspace
      r = min(n, ncons); lwork = 2*n + r*(r + 5)/2 + 2*ncons + 1
      allocate (work(lwork)); work = 0.0_dp
      !> Get the constraints matrix and vector.
      allocate (G(1, 1), h(1)); G = 0.0_dp; h = 0.0_dp
      !> Solve the QP problem.
      info = 1 ! P is already factorized when defining the QP.
      call qpgen2(P, q, n, n, x, y, obj, G, h, n, ncons, neq, iact, nact, iter, work, info)

      !> Check QuadProg info message.
      call check(error, info == 0, &
                 "QP solver did not converge.")

      !> Check solution correctness.
      call check(error, norm2(matmul(P, x) - q) < rtol, &
                 "Unconstrained QP solution is not accurate.")
      if (allocated(error)) return
   end subroutine

   subroutine test_problem_1(error)
      type(error_type), allocatable, intent(out) :: error
      ! Size of the problem.
      integer, parameter :: n = 3
      ! Quadratic Problem.
      real(dp) :: P(n, n), q(n)
      type(qp_problem) :: problem
      type(OptimizeResult) :: result
      integer :: i, info

      ! Setup the quadratic form.
      P = 0.0_dp; forall (i=1:n) P(i, i) = 1.0_dp; q = 1.0_dp
      problem = qp_problem(P, q)

      ! Solve the QP problem.
      result = solve(problem)

      !> Check QuadProg info message.
      call check(error, result%success, &
                 "QP solver did not converge.")

      !> Check solution correctness.
      call check(error, norm2(matmul(P, result%x) - q) < rtol, &
                 "Unconstrained QP solution is not accurate.")
      if (allocated(error)) return
   end subroutine

   !-------------------------------------------
   !-----                                 -----
   !-----     CONSTRAINED QP PROBLEMS     -----
   !-----                                 -----
   !-------------------------------------------

   subroutine test_problem_2(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 3
      ! Size of the problem.
      type(qp_problem) :: problem
      type(OptimizeResult) :: result
      real(dp) :: P(n, n), q(n), C(n, n), d(n)
      integer :: i, info

      ! Setup problem.
      P = 0.0_dp; forall (i=1:n) P(i, i) = 1.0_dp
      q = [0.0_dp, 5.0_dp, 0.0_dp]
      C(1, :) = [-4, 2, 0]
      C(2, :) = [-3, 1, -2]
      C(3, :) = [0, 0, 1]
      C = transpose(C); d = [-8, 2, 0]

      problem = qp_problem(P, q, C=C, d=d)
      result = solve(problem)

      !> Check QuadProg info message.
      call check(error, result%success, &
                 "QP solver did not converge.")
      if (allocated(error)) return
      !> Check solution correctness.
      block
         real(dp) :: xref(n)
         ! Reference solution.
         real(dp) :: yref(n), obj_ref
         ! Reference Lagrange multipliers, reference cost.

         xref = [0.4761905_dp, 1.0476190_dp, 2.0952381_dp]
         yref = [0.0_dp, 0.2380952_dp, 2.0952381_dp]
         obj_ref = -2.380952380952381_dp

         !> Constrained solution.
         call check(error, maxval(abs(xref - result%x)) < 1e-6_dp, &
                    "Constrained solution is not correct.")
         if (allocated(error)) return
         !> Lagrange multipliers.
         call check(error, maxval(abs(yref - result%y)) < 1e-6_dp, &
                    "Lagrange mutilpliers are not correct.")
         if (allocated(error)) return
         !> Objective value.
         call check(error, abs(obj_ref - result%obj) < 1e-6_dp, &
                    "Minimum cost is not correct.")
         if (allocated(error)) return
      end block
   end subroutine

   ! Test case from https://github.com/quadprog/quadprog/issues/2#issue-443570242
   subroutine test_problem_3(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 5, m = 3
      ! Size of the problem.
      real(dp) :: P(n, n), q(n)
      real(dp), allocatable :: A(:, :), b(:)
      type(qp_problem) :: problem
      type(OptimizeResult) :: result
      integer :: i, info
      real(dp) :: xref(n), yref(m), obj_ref
      ! Reference solution, Lagrange multipliers, reference cost.

      ! Setup problem.
      P = 0.0_dp
      do concurrent(i=1:n)
         P(i, i) = 1.0_dp
      end do
      q = [0.73727161_dp, 0.75526241_dp, 0.04741426_dp, -0.11260887_dp, -0.11260887_dp]
      allocate (A(n, m))
      A(1, :) = [3.6_dp, 0.0_dp, -9.72_dp]
      A(2, :) = [-3.4_dp, -1.9_dp, -8.67_dp]
      A(3, :) = [-3.8_dp, -1.7_dp, 0.0_dp]
      A(4, :) = [1.6_dp, -4.0_dp, 0.0_dp]
      A(5, :) = [1.6_dp, -4.0_dp, 0.0_dp]
      A = transpose(A)
      b = [1.02_dp, 0.03_dp, 0.081_dp]

      problem = qp_problem(P, q, A=A, b=b)
      result = solve(problem)

      !> Check QuadProg info message.
      call check(error, result%success, &
                 "QP solver did not converge.")
      if (allocated(error)) return

      !> Check solution correctness.
      xref = [0.07313507_dp, -0.09133482_dp, -0.08677699_dp, 0.03638213_dp, 0.03638213_dp]
      yref = [0.0440876_dp, 0.01961271_dp, 0.08465554_dp]
      obj_ref = 0.0393038880729888_dp

      !> Constrained solution.
      call check(error, norm2(xref - result%x) < 1e-6_dp, "Constrained solution is not correct.")
      if (allocated(error)) return
      !> Lagrange multipliers.
      call check(error, norm2(yref - result%y) < 1e-6_dp, "Lagrange mutilpliers are not correct.")
      if (allocated(error)) return
      !> Objective value.
      call check(error, abs(obj_ref - result%obj) < 1e-6_dp, "Minimum cost is not correct.")
      if (allocated(error)) return
      return
   end subroutine

end module
