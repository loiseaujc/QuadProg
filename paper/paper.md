---
title: '[Modern QuadProg]{.sc}: a modern, open-source Fortran implementation of the Goldfarb-Idnani algorithm for convex quadratic programming'
tags:
  - Fortran
  - Quadratic Programming
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

`Modern QuadProg` is a modernized implementation of the original FORTRAN 77 `quadprog` solver written by Berwin A. Turlach to solve strictly convex quadratic programs (QP) of the form

$$
\begin{aligned}
  \mathrm{minimize}   \quad & \dfrac12 x^\top Px - x^\top q \\
  \mathrm{subject~to} \quad & Ax = b \\
                            & Cx \geq d,
\end{aligned}
$$

where $P \in \mathbb{R}^{n \times n}$ is a symmetric positive-definite matrix.
It is based on the *active set method* by @goldfarb-idnani.
The solver is most efficient for small to moderate-sized QP described using dense matrices for which high-accuracy solutions are required.
A specialized implementation is also provided when the constraints are described by a set of sparse equations.

# Statement of need

Many problems in science and engineering can be formulated as convex quadratic problems, e.g. support vector machines in machine learning, Markowitz portfolio optimization in financial mathematics, or linear model predictive control in system engineering.
Solving convex QPs also forms the computational bottleneck of many optimization algorithms, e.g. *Newton's method* or *sequential linear-quadratic programming* for nonlinear optimization problems.
Among the many algorithms proposed to solve convex QPs, the one by @goldfarb-idnani has proven to be efficient, numerically stable and accurate.
A popular implementation is `quadprog` by Berwin Turlach, interfaced with the `R` programming language as early as 1997 [@turlach2007quadprog].
Since then, `quadprog` has been ported to many different languages, including [JavaScript](https://github.com/albertosantini/quadprog), [Rust](https://docs.rs/quadprog/latest/quadprog/), or [Julia](https://github.com/fabienlefloch/GoldfarbIdnaniSolver.jl).
Yet, little effort within the Fortran community has been devoted to modernizing the original source code.
This contribution is part of a wider community-driven effort aiming at modernizing the Fortran ecosystem at large.

# State of the field

Open-source state-of-the-art solvers for generic convex programs often are primal-dual interior point methods based on the standard conic form of said programs [@boyd:book:2004].
These include [`SCS`](https://github.com/cvxgrp/scs) [@ocpb:16; @scs] and [`ECOS`](https://github.com/embotech/ecos) [@bib:Domahidi2013ecos], both written in C, or [`Clarabel`](https://github.com/oxfordcontrol/Clarabel.jl) [@Clarabel_2024] in Julia.
They are accessible via domain specific languages provided by packages such `cvx` [@cvx; @gb08] in MATLAB, `cvxpy` [@diamond2016cvxpy; @agrawal2018rewriting] in Python, or `Convex.jl` [@convexjl] in Julia.
Solvers specialized for convex quadratic programs include [`HiGHS`](https://github.com/ERGO-Code/HiGHS) [@huangfu2018parallelizing] (active set method in C++), [`OSQP`](https://osqp.org/) [@osqp] (Douglas-Ratchford splitting in C/C++), or [`ProxSuite`](https://github.com/Simple-Robotics/proxsuite) [@bambade2022prox] (augmented Lagrangian method in C++) to name just a few.
A thorough comparison of their performances has been performed by @qpbenchmark.
It needs to be emphasized however that, despite the long-standing history of Fortran in scientific computing, none of these solvers are written in it.
This is a gap we wish to (partially) fill with this contribution to facilitate the use of convex programming techniques for Fortran programmers.

# Software design

## A modernized implementation

Written in FORTRAN 77, the original `quadprog` implementation makes use of language features now considered obsolete and calls to the outdated `linpack` library, potentially hindering the use of modern CPU instructions or hardware acceleration.
In our modernization effort, the most important updates to the original code include:

- Sources have been translated from FORTRAN 77 fixed-form to Fortran 90 free-from.
- All obsolescent features (`goto`, `continue`, etc) have been removed and the code is now fully compliant with the Fortran 2018 standard.
- Calls to appropriate `blas` and `lapack` functions now replace hand-crafted or `linpack` implementations for improved performances.
- Build system based on the Fortran Package Manager `fpm`.

We also provide modern object-oriented interfaces (see `qp_problem` and `solve`) as well as utility functions for non-negative least-squares (`nnls`), bounded-variables least-squares (`bvls`) and $\ell_1$ regularized least-squares (`lasso`).

**High-level interface -** Given the datum $P$, $q$, $A$, $b$, $C$ and $d$ defining the quadratic problem, a `qp_problem` instance can be created as follows

```fortran
problem = qp_problem(P, q, A=A, b=b, C=C, d=d)
```

where `A`, `b`, `C` and `d` are optional arguments.
Once instantiated, this problem can be solved with

```fortran
solution = solve(problem)
```

where `solution` is a derived-type storing the solution of the constrained QP, the vector of Lagrange multipliers and the value of the objective function evaluated at the constrained solution.

**Example -** The following program illustrates how to use `Modern QuadProg` to solve the constrained QP

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
    use iso_fortran_env, only: dp => real64
    use stdlib_linalg, only: eye
    use quadprog, only: OptimizeResult, qp_problem, solve
    implicit none
    integer, parameter :: n = 3                 ! Number of variables.
    real(dp) :: P(n, n), q(n), C(n, n), d(n)    ! Quadratic cost and constraints.
    type(qp_problem) :: prob
    type(OptimizeResult) :: solution
    !> Quadratic cost function.
    P = eye(n) ; q = [0.0_dp, 5.0_dp, 0.0_dp]
    !> Linear inequality constraints.
    C(1, :) = [-4.0_dp, 2.0_dp, 0.0_dp]  ; d(1) = -8.0_dp
    C(2, :) = [-3.0_dp, 1.0_dp, -2.0_dp] ; d(2) = -2.0_dp
    C(3, :) = [0.0_dp, 0.0_dp, 1.0_dp]   ; d(3) = 0.0_dp
    !> Solve the inequality constrained QP.
    prob = qp_problem(P, q, C=C, d=d)
    solution = solve(prob)
    if (solution%success) then
        print *, "x   =", solution%x   ! Solution of the QP.
        print *, "y   =", solution%y   ! Lagrange multipliers.
        print *, "obj =", solution%obj ! Objective function.
    endif
end program
```

More examples can be found in the dedicated folder [here](https://github.com/loiseaujc/QuadProg/tree/main/example).
These include a linear MPC controller with bounded actuation, and a Markowitz portfolio optimization problem with long-only positions.

## Performance considerations


| Problem ID  | Number of variables | Number of constraints | Legacy | Modern | Speed-up |
|:-----------:|:-------------------:|:---------------------:|:------:|:-------:|:--------:|
|       HS118 | 15                  | 64                    | 26µs   | 41µs | 0.63x |
|   LASER     | 1002                | 4004                  | 9.27s   | 2.96s | 3.1x |
|   AUG3DCQP  | 3873                | 8746                  | 85.1s    | 31.3s | 2.7x |

The table above reports the wall-clock time needed by the legacy and modernized implementations to solve three representative problems from the Maros-Meszaros test suite [@maros-meszaros].
This test suite contains 138 convex quadratic problems.
Following the methodology in @qpbenchmark, we extracted a subset of 25 of them corresponding to strictly convex problems having fewer than 4000 optimization variables and 10 000 constraints.
All computations were performed on a computer equipped with an 11th Gen Intel Core i7-11850H @ 2.50 GHz.
Both the legacy and modernized solvers have been compiled with `gfortran 14.3.0` and the standard options `-march=native -mtune=native -03`.
The `blas`/`lapack` backend used is Intel's Math Kernel Library `mkl`.
Both solvers providing the option to use a pre-factorized matrix $P$, we restrict the timings to the actual solve only.
The modernized implementation outperforms the legacy one, with an average speed-up of 2 to 3 for problems having roughly 50 optimization variables or more, and has similar performances for smaller problems (in the tenths of microseconds range).
Similar results have been observed using `ifx` or `flang-new` as compilers.
More details about this benchmark can be found in the [quadprog_benchmark](https://github.com/loiseaujc/quadprog_benchmark) Github repository.

## Limitations

**Strict convexity :** `Modern QuadProg` (and its ancestor) is limited to strictly convex QP.
If the matrix $P$ is only symmetric positive semi-definite, it can be replaced with $P + \varepsilon I$ (with $\varepsilon > 0$) at the expense of solving a slightly perturbed (albeit now strictly convex) problem.
In most applications, this small regularization hardly changes the result of the optimizer while robustifying the solution process.
Another alternative would be to implement the extension of the Goldfarb & Idnani algorithm for non-strictly convex QP by @boland1996dual. 

**Lack of bindings with other languages :** We do not currently provide binding to other languages.
Interfacing `Fortran` with `Python` can however be done relatively easily using utilities such as `f2py` [@f2py] or [`f90wrap`](https://github.com/jameskermode/f90wrap) [@f90wrap].
Similar packages exist for other languages (e.g. `R` or `Julia`).
Moreover that the latest `Fortran` standards have introduced many features to facilitate interoperability with the `C` language.

# Research impact statement

# Acknowledgements

No AI tools were used for software creation, documentation, nor paper authoring.
We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-33-CE46-0008-CONMAN grant agreement.

# References
