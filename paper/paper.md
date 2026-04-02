---
title: 'Modern Quadprog: a modern, open-source Fortran implementation of the Goldfarb-Idnani algorithm for convex quadratic programming'
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

`Modern Quadprog` is a modernized implementation of the original FORTRAN 77 `quadprog` solver written by Berwin A. Turlach to solve strictly convex quadratic programs (QP) of the form

$$
\begin{aligned}
  \mathrm{minimize}   \quad & \dfrac12 x^\top Px - x^\top q \\
  \mathrm{subject~to} \quad & Ax = b \\
                            & Cx \geq d,
\end{aligned}
$$

where $P \in \mathbb{R}^{n \times n}$ is a symmetric positive-definite matrix.
It is based on the *active set method* by @goldfarb-idnani, excelling on small-to-moderate dense QPs demanding high-accuracy solutions.
A specialized implementation is also provided when the constraints are described by a set of sparse equations.
This software is released under the GNU GPL-v3 license.

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

**Example -** The following program illustrates how to use `Modern Quadprog` to solve the constrained QP

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

## Limitations

**Strict convexity:** `Modern Quadprog` (and its ancestor) is limited to strictly convex QP.
If $P$ is only positive semi-definite, replacing it with $P + \varepsilon I$ (with $\varepsilon > 0$) yields a slightly perturbed, albeit strictly convex problem with negligible effect on the solution.
Another alternative would be to implement the extension of the Goldfarb & Idnani algorithm for non-strictly convex QP by @boland1996dual. 

**Lack of bindings with other languages:** We do not currently provide binding to other languages.

# Research impact statement

We follow the methodology in @qpbenchmark to evaluate the performances of the modernized implementation.
We restrict ourselves to problems from the Maros-Meszaros test suite [@maros-meszaros] with fewer than 4000 variables and 10 000 constraints.
All computations were performed on a computer equipped with an 11th Gen Intel Core i7-11850H @ 2.50 GHz.
More details about this benchmark can be found in the [quadprog_benchmark](https://github.com/loiseaujc/quadprog_benchmark) Github repository.

## Legacy vs modernized solver


| Problem ID  | Number of variables | Number of constraints | Legacy | Modern | Speed-up |
|:-----------:|:-------------------:|:---------------------:|:------:|:-------:|:--------:|
|   HS118     | 15                  | 64                    | 26µs   | 41µs | 0.63x |
|   LASER     | 1002                | 4004                  | 9.27s   | 2.96s | 3.1x |
|   AUG3DCQP  | 3873                | 8746                  | 85.1s    | 31.3s | 2.7x |

The table above reports the wall-clock time needed by the legacy and modernized implementations to solve three representative problems from the Maros-Meszaros test suite [@maros-meszaros] restricted to strictly convex quadratic programs.
Both the legacy and modernized solvers have been compiled with `gfortran 14.3.0` and the standard options `-march=native -mtune=native -03`.
The `blas`/`lapack` backend used is Intel's Math Kernel Library `mkl`, although it can easily be switched to another backend thanks to the Fortran package manager `fpm`.
Both solvers providing the option to use a pre-factorized matrix $P$, we restrict the timings to the actual solve only.
The modernized implementation outperforms the legacy one, with an average speed-up of 2 to 3 for problems having roughly 50 optimization variables or more, and has similar performances for smaller problems (in the tenths of microseconds range).
Similar results have been observed using `ifx` or `flang-new` as compilers.

## Benchmarking against other QP solvers

We now benchmark the modernized solver against a large collection of QP solvers (see the dedicated [repository](https://github.com/loiseaujc/quadprog_benchmark) for more details).
The settings used for each solver follows the high-accuracy settings given by @qp_benchmark.
For fair comparison across all solvers considered, all non-strictly convex QP in this test suite are regularized by replacing $P$ with $P + \varepsilon I$ where $\varepsilon = 10^{-8}$.

# Conclusion

`Modern Quadprog` delivers a robust, high‑performance Fortran implementation of the Goldfarb–Idnani algorithm, effective for small-to-moderate dense convex quadratic programs.
While the current release supports only strictly convex problems, a simple regularisation of the objective matrix suffices to treat non‑strictly convex cases.
Leveraging the enhanced Fortran-to-C interoperability in recent standards, the next phase of development will add bindings for Python, R, and Julia to broaden adoption and maximize research impact.

# Acknowledgements

No AI tools were used for software creation, documentation, nor paper authoring.
We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-33-CE46-0008-CONMAN grant agreement.
We are grateful to the fortran-lang community for the development of the Fortran standard library `stdlib` [@federico_perini_2026_18346789] and the Fortran Package manager `fpm`, two critical components in this modernization effort.

# References
