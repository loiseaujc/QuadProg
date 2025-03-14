module TestProblems
    use iso_fortran_env, only: output_unit
    use stdlib_math, only: is_close, all_close
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use QuadProg
    implicit none
    private

    integer , parameter, private :: dp = selected_real_kind(15, 307)
    real(dp), parameter, private :: atol = 10.0_dp ** (-precision(1.0_dp))
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
                    new_unittest("Inequality Constrained QP problem", test_problem_2), &
                    new_unittest("Equality Constrained QP problem", test_problem_3) &
                    ]
    end subroutine

    !---------------------------------------------
    !-----                                   -----
    !-----     UNCONSTRAINED QR PROBLEMS     -----
    !-----                                   -----
    !---------------------------------------------

    subroutine test_problem_1(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 3
        ! Size of the problem.
        real(dp) :: P(n, n), q(n), x(n)
        ! Quadratic form f(x) = 0.5 * x.T @ P @ x - x.T @ q
        integer :: i
        
        ! Setup the quadratic form.
        P = 0.0_dp ; forall(i=1:n) P(i, i) = 1.0_dp
        q = 1.0_dp ; x = 0.0_dp

        block
        real(dp) :: A(1, 1), b(1), qvec(n), Pmat(n, n)
        integer :: ncons, meq
        ! Constraints
        real(dp) :: y(1), obj
        ! Lagrange multiplier + Minimum cost.
        integer :: iact(1), nact
        ! Indices of active constraints + number of active constraints.
        integer :: iter(2), ierr
        ! Number of iterations + Error code.
        real(dp), allocatable :: work(:)
        integer :: lwork, r
        ! Workspace.

        !> Allocate workspace.
        ncons = 0 ; meq = 0 ; r = min(n, ncons)
        lwork = 2*n + r*(r+5)/2 + 2*ncons + 1
        allocate(work(lwork)) ; work = 0.0_dp

        ierr = 0 ; qvec = q ; Pmat = P

        !> Solve the unconstrained QP problem.
        call qpgen2(Pmat, qvec, n, n, x, y, obj, A, b, 0, ncons, meq, iact, nact, iter, work, ierr)

        !> Check QuadProg info message.
        call check(error, ierr == 0, &
                   "QP solver did not converge.")
        if (allocated(error)) then
            select case(ierr)
            case (1)
                write(output_unit, *) "Unconstrained QP has no solution."
            case (2)
                write(output_unit, *) "Cholesky decomposition of P failed."
            end select
            return
        endif

        !> Check solution correctness.
        call check(error, norm2(matmul(P, x) - q) < rtol, &
                   "Unconstrained QP solution is not accurate.")
        if (allocated(error)) return
        end block
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
        real(dp) :: P(n, n), q(n), x(n)
        ! Quadratic form f(x) = 0.5 * x.T @ P @ x - x.T @ q
        real(dp) :: C(n, n), b(n)
        ! Constrained C.T @ x >= b
        integer :: i

        ! Setup problem.
        P = 0.0_dp ; forall(i=1:n) P(i, i) = 1.0_dp
        q = [0.0_dp, 5.0_dp, 0.0_dp]
        C(1, :) = [-4, 2, 0] ; C(2, :) = [-3, 1, -2] ; C(3, :) = [0, 0, 1]
        b = [-8, 2, 0]

        block
        integer :: ncons = 3, meq = 0
        ! Constraints.
        real(dp) :: y(n), obj
        ! Lagrange multipliers + Minimum cost.
        integer :: iact(n), nact
        ! Indices of active constraints + number of active constraints.
        integer :: iter(2), ierr
        ! Number of iterations + Error code.
        real(dp), allocatable :: work(:)
        integer :: lwork, r
        ! Workspace.

        !> Allocate workspace.
        ncons = size(C, 2) ; meq = 0 ; r = min(n, ncons)
        lwork = 2*n + r*(r+5)/2 + 2*ncons + 1
        allocate(work(lwork)) ; work = 0.0_dp

        ierr = 0

        !> Solve the inequality constrained QR problem.
        call qpgen2(P, q, n, n, x, y, obj, C, b, n, ncons, meq, iact, nact, iter, work, ierr)

        !> Check QuadProg info message.
        call check(error, ierr == 0, &
                   "QP solver did not converge.")
        if (allocated(error)) then
            select case(ierr)
            case (1)
                write(output_unit, *) "Unconstrained QP has no solution."
            case (2)
                write(output_unit, *) "Cholesky decomposition of P failed."
            end select
            return
        endif

        !> Check solution correctness.
        block
        real(dp) :: xref(n), xunc(n)
        ! Reference solution, Unconstrained solution.
        real(dp) :: yref(n), obj_ref
        ! Reference Lagrange multipliers, reference cost.

        xref = [0.4761905_dp, 1.0476190_dp, 2.0952381_dp]
        xunc = [0.0_dp, 5.0_dp, 0.0_dp]
        yref = [0.0_dp, 0.2380952_dp, 2.0952381_dp]
        obj_ref = -2.380952380952381_dp

        !> Unconstrained solution.
        call check(error, all_close(xunc, q, abs_tol=1e-6_dp), &
                   "Unconstrained solution is not correct.")
        !> Constrained solution.
        call check(error, all_close(xref, x, abs_tol=1e-6_dp), &
                   "Constrained solution is not correct.")
        !> Lagrange multipliers.
        call check(error, all_close(yref, y, abs_tol=1e-6_dp), &
                   "Lagrange mutilpliers are not correct.")
        !> Objective value.
        call check(error, is_close(obj_ref, obj, abs_tol=1e-6_dp), &
                   "Minimum cost is not correct.")
        end block
        end block
    end subroutine

    ! Test case from https://github.com/quadprog/quadprog/issues/2#issue-443570242
    subroutine test_problem_3(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 5, m = 3
        ! Size of the problem.
        real(dp) :: P(n, n), q(n), x(n)
        ! Quadratic form f(x) = 0.5 * x.T @ P @ x - x.T @ q
        real(dp) :: C(n, m), b(m)
        ! Constrained C.T @ x = b
        integer :: i

        ! Setup problem.
        P = 0.0_dp ; forall(i=1:n) P(i, i) = 1.0_dp
        q = [0.73727161_dp, 0.75526241_dp, 0.04741426_dp, -0.11260887_dp, -0.11260887_dp]
        C(1, :) = [3.6_dp, 0.0_dp, -9.72_dp]
        C(2, :) = [-3.4_dp, -1.9_dp, -8.67_dp]
        C(3, :) = [-3.8_dp, -1.7_dp, 0.0_dp]
        C(4, :) = [1.6_dp, -4.0_dp, 0.0_dp]
        C(5, :) = [1.6_dp, -4.0_dp, 0.0_dp]
        b = [1.02_dp, 0.03_dp, 0.081_dp]

        block
        integer :: ncons, meq
        ! Constraints.
        real(dp) :: y(m), obj
        ! Lagrange multipliers + Minimum cost.
        integer :: iact(n), nact
        ! Indices of active constraints + number of active constraints.
        integer :: iter(2), ierr
        ! Number of iterations + Error code.
        real(dp), allocatable :: work(:)
        integer :: lwork, r
        ! Workspace.

        !> Allocate workspace.
        ncons = m ; meq = ncons ; r = min(n, ncons)
        lwork = 2*n + r*(r+5)/2 + 2*ncons + 1
        allocate(work(lwork)) ; work = 0.0_dp

        ierr = 0

        !> Solve the inequality constrained QR problem.
        call qpgen2(P, q, n, n, x, y, obj, C, b, n, ncons, meq, iact, nact, iter, work, ierr)

        !> Check QuadProg info message.
        call check(error, ierr == 0, &
                   "QP solver did not converge.")
        if (allocated(error)) then
            select case(ierr)
            case (1)
                write(output_unit, *) "Unconstrained QP has no solution."
            case (2)
                write(output_unit, *) "Cholesky decomposition of P failed."
            end select
            return
        endif

        !> Check solution correctness.
        block
        real(dp) :: xref(n), xunc(n)
        ! Reference solution, Unconstrained solution.
        real(dp) :: yref(m), obj_ref
        ! Reference Lagrange multipliers, reference cost.

        xref = [0.07313507_dp, -0.09133482_dp, -0.08677699_dp, 0.03638213_dp, 0.03638213_dp]
        yref = [0.0440876_dp, 0.01961271_dp, 0.08465554_dp]
        obj_ref = 0.0393038880729888_dp

        !> Constrained solution.
        call check(error, all_close(xref, x, abs_tol=1e-6_dp), &
                   "Constrained solution is not correct.")
        !> Lagrange multipliers.
        call check(error, all_close(yref, y, abs_tol=1e-6_dp), &
                   "Lagrange mutilpliers are not correct.")
        !> Objective value.
        call check(error, is_close(obj_ref, obj, abs_tol=1e-6_dp), &
                   "Minimum cost is not correct.")
        end block
        end block
    end subroutine
        
end module
