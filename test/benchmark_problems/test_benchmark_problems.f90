module benchmark_utils
   use stdlib_io_npy, only: load_npy
   use quadprog, only: qp_problem, dp
   implicit none
   private

   character(len=255), dimension(19), parameter, public :: problems = [character(len=255) :: &
                                                                       "data/HS21", &
                                                                       "data/QPTEST", &
                                                                       "data/HS35MOD", &
                                                                       "data/HS35", &
                                                                       "data/HS76", &
                                                                       "data/HS268", &
                                                                       "data/S268", &
                                                                       "data/DUALC5", &
                                                                       "data/DUALC1", &
                                                                       "data/HS118", &
                                                                       "data/KSIP", &
                                                                       "data/DUAL4", &
                                                                       "data/QPCBLEND", &
                                                                       "data/DUAL1", &
                                                                       "data/DUAL2", &
                                                                       "data/DUAL3", &
                                                                       "data/QPCSTAIR", &
                                                                       "data/MOSARQP2", &
                                                                       "data/LASER" &
                                                                       ]

   public :: get_problem_data, problem_description
contains

   integer function problem_description(problem, pname) result(ierr)
      type(qp_problem), intent(in) :: problem
      character(len=*), intent(in) :: pname

      print *
      print *, "--------------------"
      print *, "     * Problem ID :                           ", pname(6:)
      print *, "         - Number of variables              : ", size(problem%P, 1)
      print *, "         - Number of equality constraints   : ", problem%neq
      print *, "         - Number of inequality constraints : ", problem%ncons - problem%neq
      print *
      return
   end function

   type(qp_problem) function get_problem_data(fname) result(problem)
      character(len=*), intent(in) :: fname
      real(dp), allocatable :: P(:, :), q(:)
      real(dp), allocatable :: A(:, :), b(:)
      real(dp), allocatable :: C(:, :), d(:)

      !> Load P and q for the quadratic form.
      call load_npy(trim(fname)//"/P_matrix.npy", P)
      call load_npy(trim(fname)//"/q_vector.npy", q)

      !> Load A and b for the equality constraints.
      call load_npy(trim(fname)//"/A_matrix.npy", A)
      call load_npy(trim(fname)//"/b_vector.npy", b)

      !> Load A and b for the inequality constraints.
      call load_npy(trim(fname)//"/C_matrix.npy", C)
      call load_npy(trim(fname)//"/d_vector.npy", d)

      if (size(A, 2) == 1) A = transpose(A)

      if ((size(A, 2) == 1) .and. (size(C, 2) == 1)) then
         problem = qp_problem(P, -q)
      else if ((size(A, 2) == 1) .and. (size(C, 2) /= 1)) then
         problem = qp_problem(P, -q, C=C, d=d)
      else if ((size(A, 2) /= 1) .and. (size(C, 2) == 1)) then
         problem = qp_problem(P, -q, A=A, b=b)
      else
         problem = qp_problem(P, -q, A=A, b=b, C=C, d=d)
      end if

      return
   end function
end module benchmark_utils

module test_benchmark_problems
   use iso_fortran_env, only: output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use QuadProg
   use benchmark_utils
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
      testsuite = [new_unittest("Maros-Mezaros Benchmarks", test_problems)]
   end subroutine

   !---------------------------------------------
   !-----                                   -----
   !-----     UNCONSTRAINED QR PROBLEMS     -----
   !-----                                   -----
   !---------------------------------------------
   subroutine test_problems(error)
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

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_benchmark_problems, only: collect_suite
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
