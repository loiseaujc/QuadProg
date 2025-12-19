![](img/quadprog_logo.png)

---
Modern Fortran Edition of the `quadprog` solver

### Status

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![GitHub release](https://img.shields.io/github/release/loiseaujc/quadprog.svg)](https://github.com/loiseaujc/quadprog/releases/latest)
[![Build Status](https://github.com/loiseaujc/quadprog/actions/workflows/ci.yml/badge.svg)](https://github.com/loiseaujc/quadprog/actions)
[![codecov](https://codecov.io/gh/loiseaujc/quadprog/branch/main/graph/badge.svg)](https://codecov.io/gh/loiseaujc/quadprog)
[![last-commit](https://img.shields.io/github/last-commit/loiseaujc/quadprog)](https://github.com/loiseaujc/quadprog/commits/main)

### Description

This is an updated version of the `quadprog` solver initially written by Berwin A. Turlach in Fortran 77. It can be used to solve strictly convex quadratic programs of the form

$$
\begin{aligned}
\mathrm{minimize}   &   \quad   \dfrac12 \mathbf{x}^\top \mathbf{Px} - \mathbf{x}^\top \mathbf{q} \\
\mathrm{subject~to} &   \quad   \mathbf{Ax} = \mathbf{b}    \\
                    &   \quad   \mathbf{Cx} \geq \mathbf{d}.
\end{aligned}
$$

using an [active set method](https://en.wikipedia.org/wiki/Active-set_method). It is most efficient for small to moderate sized QP described using dense matrices.

**Updates to the original code include:**
 - Sources have been translated from FORTRAN 77 fixed-form to Fortran 90 free-form.
 - All obsolescent features (`goto`, etc) have been removed. It is now 100% standard compliant (Fortran 2018).
 - It makes use of derived-type and easy to use interfaces. The `qp_problem` class is used to defined the quadratic program and `solve` to compute its solution.
 - Calls to `blas` functions and subroutines now replace some hand-crafted implementations for improved performances.
 - Calls to `lapack` subroutines now replace the original functionalities provided by [`linpack`](https://www.netlib.org/linpack/).
 - Utility solvers for non-negative least-squares (`nnls`) and bounded-variables least-squares (`bvls`) are provided.


### Building QuadProg

#### Fortran Package Manager

The library can be built with the [Fortran Package Manager](https://github.com/fortran-lang/fpm) `fpm` using the provided `fpm.toml` file like so:

```bash
fpm build --release
```

Only double precision (`real64`) is currently supported.

To use `QuadProg` within your `fpm` project, add the following to your `fpm.toml` file:

```toml
[dependencies]
QuadProg = { git="https://github.com/loiseaujc/QuadProg.git"}
```

### Example

The following program solves the QP

$$
\begin{aligned}
\mathrm{minimize}   &   \quad   \dfrac12 \left( x_1^2 + x_2^2 + x_3^2 \right) - 5 x_2   \\
\mathrm{subject~to} &   \quad   -4x_1 + 2x_2 \geq -8    \\
                    &   \quad   -3x_1 + x_2 - 2x_3 \geq -2  \\
                    &   \quad   x_3 \geq 0.
\end{aligned}
$$

```fortran
program example
    use quadprog
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    ! Size of the problem.
    integer, parameter :: n = 3
    ! Quadratic cost.
    real(dp) :: P(n, n), q(n)
    ! Inequality constraints.
    real(dp) :: C(n, n), d(n)
    ! Convenience types.
    type(qp_problem) :: prob
    type(OptimizeResult) :: solution
    ! Miscellaneous
    integer :: i

    !> Setup the quadratic function..
    P = 0.0_dp ; q = [0.0_dp, 5.0_dp, 0.0_dp]
    do i = 1, n
        P(i, i) = 1.0_dp
    enddo
    
    !> Setup the inequality constraints.
    C(:, 1) = [-4.0_dp, 2.0_dp, 0.0_dp]
    C(:, 2) = [-3.0_dp, 1.0_dp, -2.0_dp]
    C(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
    d = [-8.0_dp, -2.0_dp, 0.0_dp]

    !> Solve the inequality constrained QP.
    prob = qp_problem(P, q, C=C, d=d)
    solution = solve(prob)

    if (solution%success) then
        print *, "x   =", solution%x ! Solution of the QP.
        print *, "y   =", solution%y ! Lagrange multipliers.
        print *, "obj =", solution%obj ! Objective function.
    endif
end program
```

### Licence

- The original source code was released under the [GNU General Public Licence version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html) (GPL-2).

### Development
 
- Development continues on [Github](https://github.com/loiseaujc/QuadProg).

### Similar projects

- [`quadprog`](https://github.com/albertosantini/quadprog), an equivalent project written in JavaScript.
- [`quadprog`](https://docs.rs/quadprog/latest/quadprog/), an equivalent project written in Rust.
- [`quadprog`](https://rdrr.io/cran/quadprog/), some R bindings to the original `quadprog` solver.
- [`GoldfarbIdnaniSolver.jl`](https://github.com/fabienlefloch/GoldfarbIdnaniSolver.jl), a port of the original `quadprog` code in Julia.

### References

- D. Goldfarb and A. Idnani (1983). A numerically stable dual method for solving strictly convex quadratic programs. Mathematical Programming, 27, 1-33.
