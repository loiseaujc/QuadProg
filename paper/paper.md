---
title: '`Modern QuadProg`: a modernized implementation of the classic `quadprog` Fortran solver'
tags:
  - Fortran
  - Quadratic Programs
  - Convex optimization
authors:
  - name: Jean-Christophe Loiseau
    orcid: 0000-0002-7244-8416
    corresponding: true
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Arts et Métiers Institute of Technology
   index: 1
date: 19 December 2025
bibliography: paper.bib
---

# Summary

`Modern QuadProg` is a modernized implementation of the original Fortran 77 `quadprog` solver written by Berwin A. Turlach to solve strictly convex quadratic programs (QP) of the form

$$
\begin{aligned}
  \mathrm{minimize} \quad & \dfrac12 \mathbf{x}^T \mathbf{Px} - \mathbf{x}^T \mathbf{q} \\
  \mathrm{subject~to} \quad & \mathbf{Ax} = \mathbf{b} \\
                             & \mathbf{Cx} \geq \mathbf{d},
\end{aligned}
$$

where $\mathbf{P} \in \mathbb{R}^{n \times n}$ is a symmetric positive-definite matrix.
Based on the method by @goldfarb-idnani, the linear constraints are handled using an *active set method*.
The solver is most efficient for small to moderate sized QP described using dense matrices.
A specialized implementation is also provided when the constraints are described by a set of sparse equations.

# Statement of need

Many problems in science and engineering can be formulated as convex quadratic problems. A non-exhaustive list includes: support vector machines in machine learning, Markowitz portfolio optimization in financial mathematics, or linear model predictive control in system engineering. Additionally, solving convex QPs forms the computational bottleneck of many optimization algorithms, e.g. *Newton's method* or *sequential linear-quadratic programming* for nonlinear optimization problems with twice differentiable objective and constraints.

## A modernized implementation

Among the many algorithms proposed to solve convex QPs, the one by @goldfarb-idnani has proven to be efficient, numerically stable and accurate. A popular implementation is `quadprog` by Berwin Turlach, interfaced with the `R` programming language as early as 1997 [@turlach2007quadprog]. Since then, `quadprog` has been ported to many different languages, including [JavaScript](https://github.com/albertosantini/quadprog), [Rust](https://docs.rs/quadprog/latest/quadprog/), or [Julia](https://github.com/fabienlefloch/GoldfarbIdnaniSolver.jl).
Yet, little effort within the Fortran community has been devoted to modernizing the original source code. This contribution is a step in this direction. It is part of a wider community-driven effort aiming at modernizing the Fortran ecosystem at large.

Written in FORTRAN 77, the original `quadprog` implementation makes use of language features now considered as obsolete. Moreover, `blas` and `lapack` being not as well established back then as they are today, many vector and matrix-vector operations relied on simple implementations, potentially hindering the use of modern CPU instructions or hardware acceleration. In our modernization effort, the most important updates to the original code include:

- Sources have been translated from FORTRAN 77 fixed-form to Fortran 90 free-from.
- All obsolescent features (`goto`, `continue`, etc) have been removed and the code is now fully compliant with the Fortran 2018 standard.
- Calls to appropriate `blas` and `lapack` functions now replace most hand-crafted or `linpack` implementations for improved performances.

We also provide modern object-oriented interfaces (see `qp_problem`) as well as utility functions to solve non-negative least-squares (`nnls`) and bounded-variables least-squares (`bvls`).

## A modern object-oriented interface

A notable introduction in `Modern QuadProg` are object-oriented interfaces. Given the datum $\mathbf{P}$, $\mathbf{q}$, $\mathbf{A}$, $\mathbf{b}$, $\mathbf{C}$ and $\mathbf{d}$ defining the quadratic problem, a `qp_problem` instance can be created as follows

```fortran
problem = qp_problem(P, q, A=A, b=b, C=C, d=d)
```

where `A`, `b`, `C` and `d` are optional arguments. While we do check the dimensions of the different matrices and vectors are consistent, it is left to the user to make sure $\mathbf{P}$ is indeed symmetric positive definite as its factorization relies on `lapack` and makes use only of its upper triangular part. Once instantiated, this problem can be solved with

```fortran
solution = solve(problem)
```

where `solution` is a derived-type storing the solution of the constrained QP, the vector of Lagrange multipliers and the value of the objective function evaluated at the constrained solution.

# Example

The following program illustrates how to use `Modern QuadProg` to solve the constrained QP

$$
\begin{aligned}
  \mathrm{minimize} \quad & \dfrac12 \left( x_1^2 + x_2^2 + x_3^2 \right) - 5 x_2 \\
  \mathrm{subject~to} \quad & -4 x_1 + 2 x_2 \geq -8 \\
                            & -3 x_1 + x_2 - 2x_3 \geq -2 \\
                            & x_3 \geq 0
\end{aligned}
$$

```fortran
program example
    use stdlib_linalg, only: eye
    use quadprog
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    ! Size of the problem.
    integer, parameter :: n = 3
    ! Quadratic cost and inequality constraints.
    real(dp) :: P(n, n), q(n), C(n, n), d(n)
    ! Convenience types.
    type(qp_problem) :: prob
    type(OptimizeResult) :: solution
    ! Miscellaneous
    integer :: i

    !> Setup the quadratic function..
    P = eye(n) ; q = [0.0_dp, 5.0_dp, 0.0_dp]
    
    !> Setup the inequality constraints.
    C(:, 1) = [-4.0_dp, 2.0_dp, 0.0_dp]  ; d(1) = -8.0_dp
    C(:, 2) = [-3.0_dp, 1.0_dp, -2.0_dp] ; d(2) = -2.0_dp
    C(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]   ; d(3) = 0.0_dp

    !> Solve the inequality constrained QP.
    prob = qp_problem(P, q, C=C, d=d)
    solution = solve(prob) ! pass legacy=.true. to use the original F77 driver.

    if (solution%success) then
        print *, "x   =", solution%x ! Solution of the QP.
        print *, "y   =", solution%y ! Lagrange multipliers.
        print *, "obj =", solution%obj ! Objective function.
    endif
end program
```

More examples can be found in the dedicated folder [here](https://github.com/loiseaujc/QuadProg/tree/main/example). These include the construction of a linear MPC controller with bounded actuation, as well as a Markowitz portfolio optimization problem.

# Performance considerations

Beyond the source code translation from Fortran 77 to modern Fortran, computational performances have been improved by making explicit calls to the appropriate `blas` functions wherever appropriate.
Similarly, calls to deprecated `linpack` functions have been replaced by their modern `lapack` equivalent.

| Problem ID  | Number of variables | Number of constraints | Legacy | Modern | Speed-up|
|:-----------:|:-------------------:|:---------------------:|:------:|:-------:|:--------:|
|       HS118 | 15                  | 64                    | 24µs   | 31µs | 0.8x |
|   LASER     | 1002                | 4004                  | 8.71s   | 2.44s | 3.6x |
|   AUG3DCQP  | 3873                | 8746                  | 82.2s    | 36.2s | 2.3x |

The table above reports the computational time needed by the legacy and modernized implementations to solve three representative problems from the Maros-Meszaros test suite [@maros-meszaros].
This test suite contains 138 convex quadratic problems.
Following the methodology in @qpbenchmark, we extracted a subset of 25 of them corresponding to problems having fewer than 4000 optimization variables and 10 000 constraints.
All computations have been run on a computer equipped with an 11th Gen Intel Core i7-11850H @ 2.50 GHz.
Both the legacy and modernized solvers have been compiled with `gfortran 15.2.0` and the same compilation options.
The `blas`/`lapack` backend used is `openblas 0.3.29` installed using `conda`.
Both solvers providing the option to use a pre-factorized matrix $\mathbf{P}$, we restrict the timings to the actual solve only.
The modernized implementation outperforms the legacy one, with an average speed-up of 2 to 3 for problems having roughly 50 optimization variables or more and similar performances for smaller problems (in the tenths of microseconds range).
More details about this benchmark can be found in the [quadprog_benchmark](https://github.com/loiseaujc/quadprog_benchmark) Github repository.

# Limitations and perspectives

**Strict convexity :** `Modern QuadProg` (and its ancestor) is limited to strictly convex QP.
When the problem is not strictly convex, the symmetric positive semi-definite matrix $\mathbf{P}$ can be replaced with $\mathbf{P} + \varepsilon \mathbf{I}$ (with $\varepsilon > 0$) at the expense of solving a slightly perturbed (albeit now strictly convex) problem.
In most applications, this small regularization might hardly change the result of the optimizer while robustifying the solution process.
Another alternative would be to implement the extension of the Goldfarb & Idnani algorithm for non-strictly convex QP by @boland1996dual. 

**Lack of interfaces with other languages :** We do not currently provide bindings to other languages.
Interfacing `Fortran` with `Python` can however be done relatively easily using utilities such as `f2py` [@f2py] or [`f90wrap`](https://github.com/jameskermode/f90wrap) [@f90wrap].
Similar packages likely exist to interface with other languages (e.g. `R` or `Julia`).
Note moreover that the latest `Fortran` standards have introduced many features to facilitate interoperability with the `C` language.

# Acknowledgements

We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-33-CE46-0008-CONMAN grant agreement.

# References
