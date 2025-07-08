module TestBenchmarkProblems
   use iso_fortran_env, only: output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use QuadProg
   use quadprog_benchmark
   implicit none
   private

   real(dp), parameter, private :: atol = 10.0_dp**(-precision(1.0_dp))
   real(dp), parameter, private :: rtol = sqrt(atol)

   public :: collect_test_benchmark_problems

contains

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   subroutine collect_test_benchmark_problems(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Maros-Mezaros Benchmarks", test_benchmark_problems)]
   end subroutine

   !---------------------------------------------
   !-----                                   -----
   !-----     UNCONSTRAINED QR PROBLEMS     -----
   !-----                                   -----
   !---------------------------------------------
   subroutine test_benchmark_problems(error)
      type(error_type), allocatable, intent(out) :: error
      type(qp_problem) :: problem
      type(OptimizeResult) :: solution
      integer :: i

      do i = 1, size(problems)
         problem = get_problem_data(problems(i))

         !> Modern driver.
         solution = solve(problem, legacy=.false.)
         call check(error, solution%success, .true.)
         if (allocated(error)) return

         !> Legacy driver.
         solution = solve(problem, legacy=.true.)
         call check(error, solution%success, .true.)
         if (allocated(error)) return
      end do
   end subroutine

end module
