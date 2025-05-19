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
date: 13 August 2017
bibliography: paper.bib
---

# Summary

`Modern QuadProg` is a modernized implementation of the original Fortran 77 `quadprog` solver written by Berwin A. Turlach.
It can be used to solve strictly convex quadratic programs (QP) of the form

$$
\begin{aligned}
  \mathrm{minimize} \quad & \dfrac12 \mathbf{x}^T \mathbf{Px} - \mathbf{x}^T \mathbf{q} \\
  \mathrm{subject~to} \quad & \mathbf{Ax} = \mathbf{b} \\
                             & \mathbf{Cx} \geq \mathbf{d},
\end{aligned}
$$

where $\mathbf{P} \in \mathbb{R}^{n \times n}$ is a symmetric positive-definite matrix.
Based on the method originally proposed by @goldfarb-idnani, the linear equality and inequality constraints are handled using an *active set method*.
The solver is most efficient for small to moderate sized QP described using dense matrices.
A specialized implementation is also provided when the constraints are described by a set of sparse equations.

# Statement of need

Many problems in science and engineering can be formulated as convex quadratic problems. A non-exhaustive list includes: support vector machines in classical machine learning, Markowitz portfolio optimization in financial mathematics, or linear model predictive control in system engineering. Likewise, solving convex QPs forms the computational bottleneck of many optimization algorithms, e.g. *Newton's method* in convex optimization or *sequential linear-quadratic programming* for nonlinear optimization problems with twice differentiable objective and constraints.

## A modernized implementation

Among the many algorithms proposed over the years to solve convex QPs, the method by @goldfarb-idnani has proven to be one of the most numerically stable and accurate of all. A popular implementation of this algorithm is `quadprog` by Berwin Turlach, interfaced with the `R` programming language as early as 1997 by Andreas Weingessel [@turlach2007quadprog]. Since then, `quadprog` has been ported to many different languages, including [JavaScript](https://github.com/albertosantini/quadprog), [Rust](https://docs.rs/quadprog/latest/quadprog/), or [Julia](https://github.com/fabienlefloch/GoldfarbIdnaniSolver.jl).
Yet, very little effort within the Fortran community has been devoted to modernizing the Fortran source code itself. This contribution is a step in this direction. It is part of a wider community-driven effort aiming at modernizing the overall Fortran ecosystem.

Written in FORTRAN 77, the original `quadprog` implementation makes use of language features now considered as obsolete. Moreover, `blas` and `lapack` being not as well established back then as they are today, many vector and matrix-vector operations relied on simple implementations, potentially hindering the use of modern CPU instructions or hardware acceleration. In our modernization effort, the most important updates to the original code include:

- Sources have been translated from FORTRAN 77 fixed-form to Fortran 90 free-from.
- All obsolescent features (`goto`, `continue`, etc) have been removed and the code base now is fully compliant with the Fortran 2018 standard.
- Calls to appropriate `blas` functions now replace most hand-crafted implementations for improved performances.
- Calls to appropriate `lapack` now replace the functionalities originally provided by `linpack`.

While we retained the definition of the original interfaces (see `qpgen1` and `qpgen2`), we also provide modern object-oriented interfaces (see next section) as well as utility functions so solve non-negative least-squares (`nnls`) and bounded-variables least-squares (`bvls`).

## A modern object-oriented interface

A notable introduction in `Modern QuadProg` is the definition of modern object-oriented interfaces. Given the datum $\mathbf{P}$, $\mathbf{q}$, $\mathbf{A}$, $\mathbf{b}$, $\mathbf{C}$ and $\mathbf{d}$ defining the quadratic problem, a `qp_problem` instance can be created as follows

```fortran
problem = qp_problem(P, q, A=A, b=b, C=C, d=d)
```

where `A`, `b`, `C` and `d` are optional arguments. It needs to be noted that, while we do check that the dimensions of the different matrices and vectors are consistent, it is left to the user to make sure $\mathbf{P}$ is indeed symmetric as its factorization relies on `lapack` and makes use only of the upper triangular part of $\mathbf{P}$.
Once defined, this problem can be solved with

```fortran
solution = solve(problem)
```

where `solution` is a derived-type with the following attributes

- `solution%x` : Solution of the QP.
- `solution%y` : Corresponding vector of Lagrange multipliers.
- `solution%obj` : Minimum of the objective function evaluated at the constrained solution.
- `solution%success` : Boolean flag determining whether the solver terminated successfully (`solution%success = .true.`) or if the problem is unfeasible (`solution%success = .false.`).

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
    P = 0.0_dp ; q = [0.0_dp, 5.0_dp, 0.0_dp]
    do i = 1, n
        P(i, i) = 1.0_dp
    enddo
    
    !> Setup the inequality constraints.
    C(:, 1) = [-4.0_dp, 2.0_dp, 0.0_dp]
    C(:, 2) = [-3.0_dp, 1.0_dp, -2.0_dp]
    C(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
    d = [-8.0_dp, 2.0_dp, 0.0_dp]

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

Additionally, `Modern QuadProg` exposes the following specialized interfaces:

- `x = nnls(A, b)` : solve a non-negative least-squares problem.
- `x = bvls(A, b, ub, lb)` : solve the bounded least-squares problem.
- `x = lasso(A, b, lambda)`: solve the $\ell_1$-regularized least-squares problem.
- `x = elastic_net(A, b, l1, l2)` : solve the mixed $\ell_1$-$\ell_2$ regularized least-squares problem.

More examples can be found in the dedicated folder [here](https://github.com/loiseaujc/QuadProg/tree/main/example). These include the construction of a linear MPC controller with bounded actuation, an SVM classifier, as well as a Markowitz portfolio optimization problem.

# Performance considerations

## Comparison with the legacy implementation

Beyond the source code translation from Fortran 77 to modern Fortran, computational performances have been improved by making explicit calls to the appropriate `blas` functions wherever appropriate.
Similarly, the calls to the deprecated `linpack` functions have been replaced by their modern `lapack` equivalent.

| Problem n° | Number of variables | Number of constraints | Legacy | Modern QuadProg |
|:----------:|:-------------------:|:---------------------:|:------:|:---------------:|
|       # 1 |
|       # 2 |
|       # 3 |

The table above reports the computational time needed by the legacy and modernized implementations to solve three representative problems for the ??? test-suite.
The platform considered is a something-something computer with something-something CPU.
Both codes have been compiled with `gfortran 14` along with the following options: `-03 -march=native -mtune=native`.
In all cases, the modernized implementation outperforms the legacy one, with speed-up reaching almost ??x on the largest problem considered.

## Comparison with other QP solvers

We make use of the [`qpbenchmark`](https://github.com/qpsolvers/qpbenchmark/tree/main) [@qpbenchmark] utility `python` package to compare the performances of the modernized `QuadProg` implementation against a fairly complete set of alternatives.
Only the subset of dense problems from the [Maros-Mesaros](https://www.cuter.rl.ac.uk/Problems/marmes.html) test suite is being considered.

# Limitations and perspectives

**Strict convexity :** Making explicit usage of the Cholesky decomposition, `Modern QuadProg` (and its legacy ancestor) is limited to strictly convex QP (i.e. problems for which $\mathbf{P}$ is symmetric positive-definite).
Note however that, when the problem is not strictly convex, the symmetric positive semi-definite matrix $\mathbf{P}$ can be replaced with $\mathbf{P} + \varepsilon \mathbf{I}$ at the expense of solving a slightly perturbed (albeit now strictly convex) problem. In most applications, this small regularization might hardly change the result of the optimizer while robustifying the solution process.

**Lack of interfaces with other languages :**

**Build system :**


# Acknowledgements

We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-33-CE46-0008-CONMAN grant agreement.

# References
