module test_lstsq_variants
   use iso_fortran_env, only: output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use QuadProg
   implicit none
   private

   real(dp), parameter, private :: atol = 10.0_dp**(-precision(1.0_dp))
   real(dp), parameter, private :: rtol = sqrt(atol)

   public :: collect_suite

contains

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   subroutine collect_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Non-negative Least-Squares", test_nnls), &
                  new_unittest("Bounded-Variables Least-Squares n°1", test_bvls_1) &
                  ]
   end subroutine

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
   end subroutine

   subroutine test_bvls_1(error)
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
   end subroutine

end module

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_lstsq_variants, only: collect_suite
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("Simple problems", collect_suite) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program
