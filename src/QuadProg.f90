module QuadProg
  use quadprog_constants
  implicit none
  private

  public :: qpgen2

  interface
!  copyright (c) 1995-2010 berwin a. turlach <berwin.turlach@gmail.com>

!  this program is free software; you can redistribute it and/or modify
!  it under the terms of the gnu general public license as published by
!  the free software foundation; either version 2 of the license, or
!  (at your option) any later version.

!  this program is distributed in the hope that it will be useful,
!  but without any warranty; without even the implied warranty of
!  merchantability or fitness for a particular purpose.  see the
!  gnu general public license for more details.

!  you should have received a copy of the gnu general public license
!  along with this program; if not, write to the free software
!  foundation, inc., 59 temple place - suite 330, boston, ma 02111-1307,
!  usa.

!  this routine uses the goldfarb/idnani algorithm to solve the
!  following minimization problem:

!        minimize  -d^t x + 1/2 *  x^t d x
!        where   a1^t x  = b1
!                a2^t x >= b2

!  the matrix d is assumed to be positive definite.  especially,
!  w.l.o.g. d is assumed to be symmetric.

!  input parameter:
!  dmat   nxn matrix, the matrix d from above (dp)
!         *** will be destroyed on exit ***
!         the user has two possibilities:
!         a) give d (ierr=0), in this case we use routines from linpack
!            to decompose d.
!         b) to get the algorithm started we need r^-1, where d=r^tr.
!            so if it is cheaper to calculate r^-1 in another way (d may
!            be a band matrix) then with the general routine, the user
!            may pass r^{-1}.  indicated by ierr not equal to zero.
!  dvec   nx1 vector, the vector d from above (dp)
!         *** will be destroyed on exit ***
!         contains on exit the solution to the initial, i.e.,
!         unconstrained problem
!  fddmat scalar, the leading dimension of the matrix dmat
!  n      the dimension of dmat and dvec (int)
!  amat   nxq matrix, the matrix a from above (dp) [ a=(a1 a2)^t ]
!         *** entries corresponding to equality constraints may have
!             changed signes on exit ***
!  bvec   qx1 vector, the vector of constants b in the constraints (dp)
!         [ b = (b1^t b2^t)^t ]
!         *** entries corresponding to equality constraints may have
!             changed signes on exit ***
!  fdamat the first dimension of amat as declared in the calling program.
!         fdamat >= n !!
!  q      integer, the number of constraints.
!  meq    integer, the number of equality constraints, 0 <= meq <= q.
!  ierr   integer, code for the status of the matrix d:
!            ierr =  0, we have to decompose d
!            ierr != 0, d is already decomposed into d=r^tr and we were
!                       given r^{-1}.

!  output parameter:
!  sol   nx1 the final solution (x in the notation above)
!  lagr  qx1 the final lagrange multipliers
!  crval scalar, the value of the criterion at the minimum
!  iact  qx1 vector, the constraints which are active in the final
!        fit (int)
!  nact  scalar, the number of constraints active in the final fit (int)
!  iter  2x1 vector, first component gives the number of "main"
!        iterations, the second one says how many constraints were
!        deleted after they became active
!  ierr  integer, error code on exit, if
!           ierr = 0, no problems
!           ierr = 1, the minimization problem has no solution
!           ierr = 2, problems with decomposing d, in this case sol
!                     contains garbage!!

!  working space:
!  work  vector with length at least 2*n+r*(r+5)/2 + 2*q +1
!        where r=min(n,q)
 
    module subroutine qpgen2(dmat, dvec, fddmat, n, sol, lagr, crval, amat, bvec, fdamat, q, meq, iact, nact, iter, work, ierr)
      integer :: fddmat, n, fdamat, q, meq, iact(*), nact, iter(*), ierr
      real(dp) :: dmat(fddmat, *), dvec(*), lagr(*), sol(*), bvec(*), work(*), amat(fdamat, *), crval
    end subroutine
  end interface

end module QuadProg
