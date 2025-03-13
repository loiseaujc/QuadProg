module TestProblems
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_test_problems

contains
    subroutine collect_test_problems(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
    end subroutine

    subroutine test_problem_1(error)
        type(error_type), allocatable, intent(out) :: error
    end subroutine
end module
