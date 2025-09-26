program check
   ! Fortran Standard.
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   ! Unit-test utilities.
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   ! Dense solver for Convex Quadratic Programs.
   use QuadProg
   ! Collection of test problems.
   use TestProblems
   use TestBenchmarkProblems
   use TestLstsqVariants
   implicit none

   ! Unit-test related.
   integer :: status, i
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   ! Collection of test suites.
   status = 0
   testsuites = [ &
                new_testsuite("Simple problems", collect_test_problems), &
                new_testsuite("Benchmark problems", collect_test_benchmark_problems), &
                new_testsuite("Least-squares variants", collect_lstsq_problems) &
                ]

   ! Run all the test suites.
   do i = 1, size(testsuites)
      write (*, *) "-------------------------------"
      write (error_unit, fmt) "Testing :", testsuites(i)%name
      write (*, *) "-------------------------------"
      write (*, *)
      call run_testsuite(testsuites(i)%collect, error_unit, status)
   end do

   if (status > 0) then
      write (error_unit, '(i0, 1x, a)') status, "test(s) failed!"
   else if (status == 0) then
      write (output_unit, *) "All tests succesfully passed!", new_line('a')
   end if
end program check
