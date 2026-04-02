module TestLstsqVariants
   use iso_fortran_env, only: output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_linalg, only: norm
   use QuadProg, only: dp, bvls, lasso, nnls
   implicit none(type, external)
   private

   real(dp), parameter, private :: atol = 10.0_dp**(-precision(1.0_dp))
   real(dp), parameter, private :: rtol = sqrt(atol)

   public :: collect_lstsq_problems

contains

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   subroutine collect_lstsq_problems(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Non-negative Least-Squares", test_nnls), &
                  new_unittest("Bounded-Variables Least-Squares", test_bvls), &
                  new_unittest("LASSO-ADMM", test_lasso) &
                  ]
   end subroutine collect_lstsq_problems

   subroutine test_nnls(error)
      type(error_type), allocatable, intent(out) :: error
      ! Size of the problem.
      integer, parameter :: m = 4, n = 3
      ! Least-Squares problem.
      real(dp) :: A(m, n), b(m)
      real(dp) :: x(n), xref(n)

      ! Setup the least-squares problem.
      A(1, :) = [1.0_dp, 1.0_dp, 2.0_dp]
      A(2, :) = [10.0_dp, 11.0_dp, -9.0_dp]
      A(3, :) = [-1.0_dp, 0.0_dp, 0.0_dp]
      A(4, :) = [-5.0_dp, 6.0_dp, -7.0_dp]

      b = [-1.0_dp, 11.0_dp, 0.0_dp, 1.0_dp]

      ! Solve the nnls problem.
      x = nnls(A, b)

      !> Check solution.
      xref = [0.4610_dp, 0.5611_dp, 0.0_dp]
      call check(error, maxval(x - xref) < 1e-4_dp, &
                 "NNLS did not find the correct solution.")
      if (allocated(error)) return
   end subroutine test_nnls

   subroutine test_bvls(error)
      type(error_type), allocatable, intent(out) :: error
      ! Size of the problem.
      integer, parameter :: m = 4, n = 3
      ! Least-Squares problem.
      real(dp) :: A(m, n), b(m)
      real(dp) :: x(n), xref(n)
      real(dp) :: ub(n), lb(n)

      ! Setup the least-squares problem.
      A(1, :) = [1.0_dp, 1.0_dp, 2.0_dp]
      A(2, :) = [10.0_dp, 11.0_dp, -9.0_dp]
      A(3, :) = [-1.0_dp, 0.0_dp, 0.0_dp]
      A(4, :) = [-5.0_dp, 6.0_dp, -7.0_dp]

      b = [-1.0_dp, 11.0_dp, 0.0_dp, 1.0_dp]

      ! Solve the nnls problem.
      ub = huge(1.0_dp); lb = 0.0_dp
      x = bvls(A, b, ub=ub, lb=lb)
      !> Check solution.
      xref = [0.4610_dp, 0.5611_dp, 0.0_dp]
      call check(error, maxval(x - xref) < 1e-4_dp, &
                 "BVLS did not find the correct solution.")
      if (allocated(error)) return

      ! Solve the nnls problem.
      x = bvls(A, b, lb=lb)
      !> Check solution.
      xref = [0.4610_dp, 0.5611_dp, 0.0_dp]
      call check(error, maxval(x - xref) < 1e-4_dp, &
                 "BVLS did not find the correct solution.")
      if (allocated(error)) return
   end subroutine test_bvls

   subroutine test_lasso(error)
      type(error_type), allocatable, intent(out) :: error
      ! Size of the problem.
      integer, parameter :: m = 4, n = 3
      ! Least-Squares problem.
      real(dp) :: A(m, n), b(m)
      real(dp) :: x(n), xref(n)
      real(dp) :: obj, obj_ref

      ! Setup the least-squares problem.
      A(1, :) = [1.0_dp, 1.0_dp, 2.0_dp]
      A(2, :) = [10.0_dp, 11.0_dp, -9.0_dp]
      A(3, :) = [-1.0_dp, 0.0_dp, 0.0_dp]
      A(4, :) = [-5.0_dp, 6.0_dp, -7.0_dp]

      b = [-1.0_dp, 11.0_dp, 0.0_dp, 1.0_dp]

      ! Solve the lasso problem.
      x = lasso(A, b, lambda=1.0_dp, tol=1e-14_dp)

      ! Check againt solution obtained with cvxpy.
      xref = [0.572242945_dp, 0.0_dp, -0.572172567_dp]
      obj_ref = 1.418256034907453_dp
      obj = 0.5_dp*norm(matmul(A, x) - b, 2)**2 + norm(x, 1)
      call check(error, maxval(x - xref) < 1e-8_dp, &
                 "LASSO did not find the correct solution.")
      if (allocated(error)) return
   end subroutine test_lasso

end module TestLstsqVariants
