!  Copyright (c) 1995-2010 Berwin A. Turlach <berwin.turlach@gmail.com>

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
submodule(quadprog) quadprog_legacy
   use quadprog_constants, only: dp
contains
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

   module procedure qpgen2
   integer  :: i, j, l, l1, info, it1, iwzv, iwrv, iwrm, iwsv, iwuv, nvl, r, iwnbv
   real(dp) :: temp, sum, t1, tt, gc, gs, nu, vsmall, tmpa, tmpb
   logical  :: t1inf, t2min

   r = min(n, q); l = 2*n + (r*(r + 5))/2 + 2*q + 1   ! Workspace size.
   vsmall = epsilon(1.0_dp)                           ! Machine precision.

   !> Store the initial dvec and initialize some arrays.
   do concurrent(i=1:n)
      work(i) = dvec(i)
   end do
   do concurrent(i=n + 1:l)
      work(i) = 0.0_dp
   end do
   do concurrent(i=1:q)
      iact(i) = 0; lagr(i) = 0.0_dp
   end do

   !------------------------------------------------------------
   !-----     SOLVE UNCONSTRAINED QP AS STARTING POINT     -----
   !------------------------------------------------------------

   if (ierr == 0) then
      !> Matrix has not been factorized yet.
      call dpofa(dmat, fddmat, n, info)   ! Cholesky factorization.
      if (info /= 0) then
         ierr = 2; return
      end if
      call dposl(dmat, fddmat, n, dvec)   ! Solve Ax = b with pre-computed Chol. fact.
      call dpori(dmat, fddmat, n)         ! Compute inv(A) from pre-computed Chol. fact.
   else
      !> Matrix has already been pre-factorized by calling dpofa and dpori.
      !> Multiply by inv(R).T
      do j = 1, n
         sol(j) = 0.0_dp
         do i = 1, j
            sol(j) = sol(j) + dmat(i, j)*dvec(i)
         end do
      end do
      !> Multiply by inv(R).
      do j = 1, n
         dvec(j) = 0.0_dp
         do i = j, n
            dvec(j) = dvec(j) + dmat(j, i)*sol(i)
         end do
      end do
   end if

   !> Book-keeping:
   !>    - set low triangular part of dmat to zero,
   !>    - store dvec in sol,
   !>    - calculate value of the criterion at unconstrained minima.
   crval = 0.0_dp
   do j = 1, n
      sol(j) = dvec(j)
      crval = crval + work(j)*sol(j)
      work(j) = 0.0_dp
      do concurrent(i=j + 1:n)
         dmat(i, j) = 0.0_dp
      end do
   end do
   crval = -crval/2.0_dp
   ierr = 0

   !> Calculate some constants, i.e., from which index on the different
   !> quantities are stored in the work matrix
   iwzv = n
   iwrv = iwzv + n
   iwuv = iwrv + r
   iwrm = iwuv + r + 1
   iwsv = iwrm + (r*(r + 1))/2
   iwnbv = iwsv + q

   !> Calculate the norm of each column of the constraint matrix
   do i = 1, q
      sum = 0.0_dp
      do j = 1, n
         sum = sum + amat(j, i)*amat(j, i)
      end do
      work(iwnbv + i) = sqrt(sum)
   end do
   nact = 0
   iter(1) = 0
   iter(2) = 0

   !-------------------------------------------------------------------------
   !-----     ACTIVE SET METHOD FOR SOLVING THE CONSTRAINED PROBLEM     -----
   !-------------------------------------------------------------------------

   loop50: do
      !> Update iteration counter.
      iter(1) = iter(1) + 1

      !> Verify all constraints.
      !>    - Check which are being violated.
      !>    - For equality ones, the normal vector may have to be negated, bvec also.
      l = iwsv
      do i = 1, q
         l = l + 1
         sum = -bvec(i)
         do j = 1, n
            sum = sum + amat(j, i)*sol(j)
         end do
         if (abs(sum) < vsmall) then
            sum = 0.0_dp
         end if
         if (i > meq) then
            work(l) = sum
         else
            work(l) = -abs(sum)
            if (sum > 0.0_dp) then
               do j = 1, n
                  amat(j, i) = -amat(j, i)
               end do
               bvec(i) = -bvec(i)
            end if
         end if
      end do

      !> Active constraints are set explicitly to zero.
      do i = 1, nact
         work(iwsv + iact(i)) = 0.0_dp
      end do

      !> Select which violated constraint to consider.
      !>    - Each violation is weighted by the number of non-zero elements in
      !>      the corresponding row of A.
      !>    - Selected violated constraint is the one with maximal absolute value.
      nvl = 0
      temp = 0.0_dp
      do i = 1, q
         if (work(iwsv + i) < temp*work(iwnbv + i)) then
            nvl = i
            temp = work(iwsv + i)/work(iwnbv + i)
         end if
      end do
      if (nvl == 0) then
         do i = 1, nact
            lagr(iact(i)) = work(iwuv + i)
         end do
         return
      end if

      !> Compute d = J.T @ n+
      !>    -  n+ is the normal vector of the violated constraint.
      !>    -  J is stored in dmat in this implementation.
      !>    -  If a constraint is dropped, we come back here.
      loop55: do
         block700: block

            do i = 1, n
               sum = 0.0_dp
               do j = 1, n
                  sum = sum + dmat(j, i)*amat(j, nvl)
               end do
               work(i) = sum
            end do

            !> Compute z = J_2 @ d_2
            l1 = iwzv
            do i = 1, n
               work(l1 + i) = 0.0_dp
            end do
            do j = nact + 1, n
               do i = 1, n
                  work(l1 + i) = work(l1 + i) + dmat(i, j)*work(j)
               end do
            end do

            !> Compute r = inv(R) @ d_1
            !>    -  Check if r has positive entries (among entries corresponding
            !>       to inequality constraints).
            t1inf = .true.
            do i = nact, 1, -1
               sum = work(i)
               l = iwrm + (i*(i + 3))/2
               l1 = l - i
               do j = i + 1, nact
                  sum = sum - work(l)*work(iwrv + j)
                  l = l + j
               end do
               sum = sum/work(l1)
               work(iwrv + i) = sum
               if (iact(i) <= meq) cycle
               if (sum <= 0.0_dp) cycle
               t1inf = .false.
               it1 = i
            end do

            !> If r has positive elements:
            !>    -  Find partial step length t1 (maximum step in dual space without
            !>       violating dual feasibility).
            !>    -  it1 stores wich component t1, the min of u/r, occurs.
            if (.not. t1inf) then
               t1 = work(iwuv + it1)/work(iwrv + it1)
               do i = 1, nact
                  if (iact(i) <= meq) cycle
                  if (work(iwrv + i) <= 0.0_dp) cycle
                  temp = work(iwuv + i)/work(iwrv + i)
                  if (temp < t1) then
                     t1 = temp
                     it1 = i
                  end if
               end do
            end if

            !> Test if the z vector is equal to zero.
            sum = 0.0_dp
            do i = iwzv + 1, iwzv + n
               sum = sum + work(i)*work(i)
            end do
            if (abs(sum) <= vsmall) then
               !> Not step in primal space such that the new constraint becomes feasible.
               !>    -  Take step in dual space and drop a constraint.
               if (t1inf) then
                  !> No step in dual possible is possible. Problem is infeasible.
                  ierr = 1; return
               else
                  !> Take partial step in dual space.
                  !>    -  Drop the it1-th active constraint.
                  !>    -  Continue at step 2(a) (marked by label 55).
                  do i = 1, nact
                     work(iwuv + i) = work(iwuv + i) - t1*work(iwrv + i)
                  end do
                  work(iwuv + nact + 1) = work(iwuv + nact + 1) + t1
                  exit block700
               end if
            else
               !> Minimum step in primal space such that the constraint becomes feasible.
               !>    -  Full step length t2.
               !>    -  Keep sum = z.T @ n+ to update crval below.
               sum = 0.0_dp
               do i = 1, n
                  sum = sum + work(iwzv + i)*amat(i, nvl)
               end do
               tt = -work(iwsv + nvl)/sum
               t2min = .true.
               if (.not. t1inf) then
                  if (t1 < tt) then
                     tt = t1
                     t2min = .false.
                  end if
               end if

               !> Take step in primal and dual space.
               do i = 1, n
                  sol(i) = sol(i) + tt*work(iwzv + i)
               end do
               crval = crval + tt*sum*(tt/2.0_dp + work(iwuv + nact + 1))
               do i = 1, nact
                  work(iwuv + i) = work(iwuv + i) - tt*work(iwrv + i)
               end do
               work(iwuv + nact + 1) = work(iwuv + nact + 1) + tt

               !> If a full step has been taken:
               !>    -  Check whether further constraints are violated.
               !>    -  If not, drop current constraint and iterate once more.
               if (t2min) then
                  !> Full step was performed.
                  !>    -  Add constraint nvl to the list of active constraints.
                  !>    -  Update J and r.
                  nact = nact + 1
                  iact(nact) = nvl
                  !> To update r:
                  !>    -  Put the first nact-1 components of the d vector into column
                  !>       nact of r.
                  l = iwrm + ((nact - 1)*nact)/2 + 1
                  do i = 1, nact - 1
                     work(l) = work(i)
                     l = l + 1
                  end do

                  !> If nact = n:
                  !>    -  Add the last element to the new row or r.
                  !> Else:
                  !>    -  Use Given rotations to turn the vector d(nact:n) into a multiple
                  !>       of the first unit vector.
                  !>    -  Multiple goes into the last element of the new row of r.
                  !>    -  J is updated accordingly to the Givens rotations.
                  if (nact == n) then
                     work(l) = work(n)
                  else
                     do i = n, nact + 1, -1
                        !> Find the Givens rotation reducing the l1-th element of d to zero.
                        !> If it is already zero, do nothing except for decreasing l1.
                        if (work(i) == 0.0_dp) cycle
                        gc = max(abs(work(i - 1)), abs(work(i)))
                        gs = min(abs(work(i - 1)), abs(work(i)))
                        temp = sign(gc*sqrt(1 + (gs/gc)*(gs/gc)), work(i - 1))
                        gc = work(i - 1)/temp
                        gs = work(i)/temp

                        !> Givens rotation is done with the matrix [gc gs, gs -gc].
                        !>    -  If gc = 1, i-th element of d is zero compared with element
                        !>       (l1-1). Hence, do nothing.
                        !>    -  If gc = 0, switch column i and column i-1 of J. Since we only
                        !>       switch columns in J, need to be careful how d is updated
                        !>       depending on the sign of gs.
                        !>    -  Else, Givens rotation is applied to the columns. The i-1
                        !>       element of d has to be updated to temp.
                        if (gc == 1.0_dp) cycle
                        if (gc == 0.0_dp) then
                           work(i - 1) = gs*temp
                           do j = 1, n
                              temp = dmat(j, i - 1)
                              dmat(j, i - 1) = dmat(j, i)
                              dmat(j, i) = temp
                           end do
                        else
                           work(i - 1) = temp
                           nu = gs/(1.0_dp + gc)
                           do j = 1, n
                              temp = gc*dmat(j, i - 1) + gs*dmat(j, i)
                              dmat(j, i) = nu*(dmat(j, i - 1) + temp) - dmat(j, i)
                              dmat(j, i - 1) = temp
                           end do
                        end if
                     end do

                     !> l is still pointing to element (nact,nact) of the matrix r.
                     !> Store d(nact) in r(nact,nact)
                     work(l) = work(nact)
                  end if
               else
                  !> A partial step in dual space was taken.
                  !>    -  Drop the it1-th active constraint.
                  !>    -  Continue to step 2(a) (marked by label 55).
                  !> Since fit changed, we need to recalculate by "how much" the chosen
                  !> constraint is now violated.
                  sum = -bvec(nvl)
                  do j = 1, n
                     sum = sum + sol(j)*amat(j, nvl)
                  end do
                  if (nvl > meq) then
                     work(iwsv + nvl) = sum
                  else
                     work(iwsv + nvl) = -abs(sum)
                     if (sum > 0.0_dp) then
                        do j = 1, n
                           amat(j, nvl) = -amat(j, nvl)
                        end do
                        bvec(nvl) = -bvec(nvl)
                     end if
                  end if
                  exit block700
               end if
            end if
            cycle loop50
            !> Drop the it1-th constraint.
         end block block700

         !> if it1 = nact
         !>    -  it is only necessary to update the vector u and nact
         if (it1 /= nact) then ! go to 799
            !> After updating one row of r (column of j), come back here.
            loop797: do
               !> Find the Givens rotation reducing the (it1+1, it+1)-th element of R
               !> to zero.
               !> If it is already zero:
               !>    - do nothing except updating u, iact and shifting column (it+1) of R
               !>      to column (it1).
               !> l will point ot element (1,it1+1) of R.
               !> l1 will point to element (it+1, it1+1) of R.
               ! we have to find the givens rotation which will reduce the element
               ! (it1+1,it1+1) of r to zero.
               l = iwrm + (it1*(it1 + 1))/2 + 1
               l1 = l + it1
               if (work(l1) /= 0.0_dp) then ! first go to 798
                  gc = max(abs(work(l1 - 1)), abs(work(l1)))
                  gs = min(abs(work(l1 - 1)), abs(work(l1)))
                  temp = sign(gc*sqrt(1 + (gs/gc)*(gs/gc)), work(l1 - 1))
                  gc = work(l1 - 1)/temp
                  gs = work(l1)/temp

                  !> Givens rotation is done with the matrix [gc gc ; gs -gc].
                  !> If gc = 0:
                  !>    - Switch row (it1) and (it1+1) of R
                  !>    - Switch column (it1= and (it+1) of J.
                  !>    - Since we switch rows in R and columns R, we can ignore sign(gs).
                  !> Else:
                  !>    - Apply Givens rotation to these rows/columns.
                  if (gc /= 1.0_dp) then ! second go to 798
                     if (gc == 0.0_dp) then
                        do i = it1 + 1, nact
                           temp = work(l1 - 1)
                           work(l1 - 1) = work(l1)
                           work(l1) = temp
                           l1 = l1 + i
                        end do
                        do i = 1, n
                           temp = dmat(i, it1)
                           dmat(i, it1) = dmat(i, it1 + 1)
                           dmat(i, it1 + 1) = temp
                        end do
                     else
                        nu = gs/(1.0_dp + gc)
                        do i = it1 + 1, nact
                           temp = gc*work(l1 - 1) + gs*work(l1)
                           work(l1) = nu*(work(l1 - 1) + temp) - work(l1)
                           work(l1 - 1) = temp
                           l1 = l1 + i
                        end do
                        do i = 1, n
                           temp = gc*dmat(i, it1) + gs*dmat(i, it1 + 1)
                           dmat(i, it1 + 1) = nu*(dmat(i, it1) + temp) - dmat(i, it1 + 1)
                           dmat(i, it1) = temp
                        end do
                     end if
                     ! shift column (it1+1) of r to column (it1) (that is, the first it1
                     ! elements). the posit1on of element (1,it1+1) of r was calculated above
                     ! and stored in l.
                  end if
               end if

               l1 = l - it1
               do i = 1, it1
                  work(l1) = work(l)
                  l = l + 1
                  l1 = l1 + 1
               end do

               !> Update vector u and iact as necessary and continue
               !> with updating the matrices J and R.
               work(iwuv + it1) = work(iwuv + it1 + 1)
               iact(it1) = iact(it1 + 1)
               it1 = it1 + 1
               if (it1 >= nact) exit loop797
            end do loop797
         end if ! go to 799

         work(iwuv + nact) = work(iwuv + nact + 1)
         work(iwuv + nact + 1) = 0.0_dp
         iact(nact) = 0
         nact = nact - 1
         iter(2) = iter(2) + 1
      end do loop55
   end do loop50

   return
   end procedure

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
!  amat   lxq matrix (dp)
!         *** entries corresponding to equality constraints may have
!             changed signes on exit ***
!  iamat  (l+1)xq matrix (int)
!         these two matrices store the matrix a in compact form. the format
!         is: [ a=(a1 a2)^t ]
!           iamat(1,i) is the number of non-zero elements in column i of a
!           iamat(k,i) for k>=2, is equal to j if the (k-1)-th non-zero
!                      element in column i of a is a(i,j)
!            amat(k,i) for k>=1, is equal to the k-th non-zero element
!                      in column i of a.

!  bvec   qx1 vector, the vector of constants b in the constraints (dp)
!         [ b = (b1^t b2^t)^t ]
!         *** entries corresponding to equality constraints may have
!             changed signes on exit ***
!  fdamat the first dimension of amat as declared in the calling program.
!         fdamat >= n (and iamat must have fdamat+1 as first dimension)
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

   module procedure qpgen1
   integer :: i, j, l, l1, info, it1, iwzv, iwrv, iwrm, iwsv, iwuv, nvl, r, iwnbv
   real(dp) :: temp, sum, t1, tt, gc, gs, nu, vsmall, tmpa, tmpb
   logical :: t1inf, t2min

   r = min(n, q)
   l = 2*n + (r*(r + 5))/2 + 2*q + 1
   vsmall = epsilon(1.0_dp)

   !> Store the initial dvec and initialize some arrays.
   do i = 1, n
      work(i) = dvec(i)
   end do
   do i = n + 1, l
      work(i) = 0.d0
   end do
   do i = 1, q
      iact(i) = 0
      lagr(i) = 0.d0
   end do

   !------------------------------------------------------------
   !-----     SOLVE UNCONSTRAINED QP AS STARTING POINT     -----
   !------------------------------------------------------------

   if (ierr == 0) then
      !> Matrix has not been factorized yet.
      call dpofa(dmat, fddmat, n, info)   ! Cholesky factorization.
      if (info /= 0) then
         ierr = 2; return
      end if
      call dposl(dmat, fddmat, n, dvec)   ! Solve Ax = b with pre-computed Chol. fact.
      call dpori(dmat, fddmat, n)         ! Compute inv(A) from pre-computed Chol. fact.
   else
      !> Matrix has already been pre-factorized by calling dpofa and dpori.
      !> Multiply by inv(R).T.
      do j = 1, n
         sol(j) = 0.d0
         do i = 1, j
            sol(j) = sol(j) + dmat(i, j)*dvec(i)
         end do
      end do
      !> Multiply by inv(R).
      do j = 1, n
         dvec(j) = 0.d0
         do i = j, n
            dvec(j) = dvec(j) + dmat(j, i)*sol(i)
         end do
      end do
   end if

   !> Book-keeping:
   !>    - set lower triangular part of dmat to zero,
   !>    - store dvec in sol,
   !>    - calculate value of the criterion at unconstrained minima
   crval = 0.d0
   do j = 1, n
      sol(j) = dvec(j)
      crval = crval + work(j)*sol(j)
      work(j) = 0.d0
      do i = j + 1, n
         dmat(i, j) = 0.d0
      end do
   end do
   crval = -crval/2.d0
   ierr = 0

   !> Calculate some constants, i.e., from which index on the different
   !> quantities are stored in the work matrix
   iwzv = n
   iwrv = iwzv + n
   iwuv = iwrv + r
   iwrm = iwuv + r + 1
   iwsv = iwrm + (r*(r + 1))/2
   iwnbv = iwsv + q

   !> Calculate the norm of each column of the constraint matrix
   do i = 1, q
      sum = 0.d0
      do j = 1, iamat(1, i)
         sum = sum + amat(j, i)*amat(j, i)
      end do
      work(iwnbv + i) = sqrt(sum)
   end do
   nact = 0
   iter(1) = 0
   iter(2) = 0

50 continue

! start a new iteration

   iter(1) = iter(1) + 1

! calculate all constraints and check which are still violated
! for the equality constraints we have to check whether the normal
! vector has to be negated (as well as bvec in that case)

   l = iwsv
   do i = 1, q
      l = l + 1
      sum = -bvec(i)
      do j = 1, iamat(1, i)
         sum = sum + amat(j, i)*sol(iamat(j + 1, i))
      end do
      if (abs(sum) < vsmall) then
         sum = 0.0d0
      end if
      if (i > meq) then
         work(l) = sum
      else
         work(l) = -abs(sum)
         if (sum > 0.d0) then
            do j = 1, iamat(1, i)
               amat(j, i) = -amat(j, i)
            end do
            bvec(i) = -bvec(i)
         end if
      end if
   end do

! as safeguard against rounding errors set already active constraints
! explicitly to zero

   do i = 1, nact
      work(iwsv + iact(i)) = 0.d0
   end do

! we weight each violation by the number of non-zero elements in the
! corresponding row of a. then we choose the violated constraint which
! has maximal absolute value, i.e., the minimum.
! by obvious commenting and uncommenting we can choose the strategy to
! take always the first constraint which is violated. ;-)

   nvl = 0
   temp = 0.d0
   do i = 1, q
      if (work(iwsv + i) < temp*work(iwnbv + i)) then
         nvl = i
         temp = work(iwsv + i)/work(iwnbv + i)
      end if
   end do
   if (nvl == 0) then
      do i = 1, nact
         lagr(iact(i)) = work(iwuv + i)
      end do
      return
   end if

! calculate d=j^tn^+ where n^+ is the normal vector of the violated
! constraint. j is stored in dmat in this implementation!!
! if we drop a constraint, we have to jump back here.

55 continue
   do i = 1, n
      sum = 0.d0
      do j = 1, iamat(1, nvl)
         sum = sum + dmat(iamat(j + 1, nvl), i)*amat(j, nvl)
      end do
      work(i) = sum
   end do

! now calculate z = j_2 d_2

   l1 = iwzv
   do i = 1, n
      work(l1 + i) = 0.d0
   end do
   do j = nact + 1, n
      do i = 1, n
         work(l1 + i) = work(l1 + i) + dmat(i, j)*work(j)
      end do
   end do

! and r = r^{-1} d_1, check also if r has positive elements (among the
! entries corresponding to inequalities constraints).

   t1inf = .true.
   do i = nact, 1, -1
      sum = work(i)
      l = iwrm + (i*(i + 3))/2
      l1 = l - i
      do j = i + 1, nact
         sum = sum - work(l)*work(iwrv + j)
         l = l + j
      end do
      sum = sum/work(l1)
      work(iwrv + i) = sum
      if (iact(i) <= meq) cycle
      if (sum <= 0.d0) cycle
      t1inf = .false.
      it1 = i
   end do

! if r has positive elements, find the partial step length t1, which is
! the maximum step in dual space without violating dual feasibility.
! it1  stores in which component t1, the min of u/r, occurs.

   if (.not. t1inf) then
      t1 = work(iwuv + it1)/work(iwrv + it1)
      do i = 1, nact
         if (iact(i) <= meq) cycle
         if (work(iwrv + i) <= 0.d0) cycle
         temp = work(iwuv + i)/work(iwrv + i)
         if (temp < t1) then
            t1 = temp
            it1 = i
         end if
      end do
   end if

! test if the z vector is equal to zero

   sum = 0.d0
   do i = iwzv + 1, iwzv + n
      sum = sum + work(i)*work(i)
   end do
   if (abs(sum) <= vsmall) then

! no step in primal space such that the new constraint becomes
! feasible. take step in dual space and drop a constant.

      if (t1inf) then
         ! no step in dual space possible either, problem is not solvable
         ierr = 1; return
      else

! we take a partial step in dual space and drop constraint it1,
! that is, we drop the it1-th active constraint.
! then we continue at step 2(a) (marked by label 55)

         do i = 1, nact
            work(iwuv + i) = work(iwuv + i) - t1*work(iwrv + i)
         end do
         work(iwuv + nact + 1) = work(iwuv + nact + 1) + t1
         go to 700
      end if
   else

! compute full step length t2, minimum step in primal space such that
! the constraint becomes feasible.
! keep sum (which is z^tn^+) to update crval below!

      sum = 0.d0
      do i = 1, iamat(1, nvl)
         sum = sum + work(iwzv + iamat(i + 1, nvl))*amat(i, nvl)
      end do
      tt = -work(iwsv + nvl)/sum
      t2min = .true.
      if (.not. t1inf) then
         if (t1 < tt) then
            tt = t1
            t2min = .false.
         end if
      end if

! take step in primal and dual space

      do i = 1, n
         sol(i) = sol(i) + tt*work(iwzv + i)
      end do
      crval = crval + tt*sum*(tt/2.d0 + work(iwuv + nact + 1))
      do i = 1, nact
         work(iwuv + i) = work(iwuv + i) - tt*work(iwrv + i)
      end do
      work(iwuv + nact + 1) = work(iwuv + nact + 1) + tt

! if it was a full step, then we check wheter further constraints are
! violated otherwise we can drop the current constraint and iterate once
! more
      if (t2min) then

! we took a full step. thus add constraint nvl to the list of active
! constraints and update j and r

         nact = nact + 1
         iact(nact) = nvl

! to update r we have to put the first nact-1 components of the d vector
! into column (nact) of r

         l = iwrm + ((nact - 1)*nact)/2 + 1
         do i = 1, nact - 1
            work(l) = work(i)
            l = l + 1
         end do

! if now nact=n, then we just have to add the last element to the new
! row of r.
! otherwise we use givens transformations to turn the vector d(nact:n)
! into a multiple of the first unit vector. that multiple goes into the
! last element of the new row of r and j is accordingly updated by the
! givens transformations.

         if (nact == n) then
            work(l) = work(n)
         else
            do i = n, nact + 1, -1

! we have to find the givens rotation which will reduce the element
! (l1) of d to zero.
! if it is already zero we don't have to do anything, except of
! decreasing l1

               if (work(i) == 0.d0) cycle
               gc = max(abs(work(i - 1)), abs(work(i)))
               gs = min(abs(work(i - 1)), abs(work(i)))
               temp = sign(gc*sqrt(1 + (gs/gc)*(gs/gc)), work(i - 1))
               gc = work(i - 1)/temp
               gs = work(i)/temp

! the givens rotation is done with the matrix (gc gs, gs -gc).
! if gc is one, then element (i) of d is zero compared with element
! (l1-1). hence we don't have to do anything.
! if gc is zero, then we just have to switch column (i) and column (i-1)
! of j. since we only switch columns in j, we have to be careful how we
! update d depending on the sign of gs.
! otherwise we have to apply the givens rotation to these columns.
! the i-1 element of d has to be updated to temp.

               if (gc == 1.d0) cycle
               if (gc == 0.d0) then
                  work(i - 1) = gs*temp
                  do j = 1, n
                     temp = dmat(j, i - 1)
                     dmat(j, i - 1) = dmat(j, i)
                     dmat(j, i) = temp
                  end do
               else
                  work(i - 1) = temp
                  nu = gs/(1.d0 + gc)
                  do j = 1, n
                     temp = gc*dmat(j, i - 1) + gs*dmat(j, i)
                     dmat(j, i) = nu*(dmat(j, i - 1) + temp) - dmat(j, i)
                     dmat(j, i - 1) = temp
                  end do
               end if
            end do

! l is still pointing to element (nact,nact) of the matrix r.
! so store d(nact) in r(nact,nact)
            work(l) = work(nact)
         end if
      else

! we took a partial step in dual space. thus drop constraint it1,
! that is, we drop the it1-th active constraint.
! then we continue at step 2(a) (marked by label 55)
! but since the fit changed, we have to recalculate now "how much"
! the fit violates the chosen constraint now.

         sum = -bvec(nvl)
         do j = 1, iamat(1, nvl)
            sum = sum + sol(iamat(j + 1, nvl))*amat(j, nvl)
         end do
         if (nvl > meq) then
            work(iwsv + nvl) = sum
         else
            work(iwsv + nvl) = -abs(sum)
            if (sum > 0.d0) then
               do j = 1, iamat(1, nvl)
                  amat(j, nvl) = -amat(j, nvl)
               end do
               bvec(nvl) = -bvec(nvl)
            end if
         end if
         go to 700
      end if
   end if
   go to 50

! drop constraint it1

700 continue

! if it1 = nact it is only necessary to update the vector u and nact

   if (it1 == nact) go to 799

! after updating one row of r (column of j) we will also come back here

797 continue

! we have to find the givens rotation which will reduce the element
! (it1+1,it1+1) of r to zero.
! if it is already zero we don't have to do anything except of updating
! u, iact, and shifting column (it1+1) of r to column (it1)
! l  will point to element (1,it1+1) of r
! l1 will point to element (it1+1,it1+1) of r

   l = iwrm + (it1*(it1 + 1))/2 + 1
   l1 = l + it1
   if (work(l1) == 0.d0) go to 798
   gc = max(abs(work(l1 - 1)), abs(work(l1)))
   gs = min(abs(work(l1 - 1)), abs(work(l1)))
   temp = sign(gc*sqrt(1 + (gs/gc)*(gs/gc)), work(l1 - 1))
   gc = work(l1 - 1)/temp
   gs = work(l1)/temp

! the givens rotatin is done with the matrix (gc gs, gs -gc).
! if gc is one, then element (it1+1,it1+1) of r is zero compared with
! element (it1,it1+1). hence we don't have to do anything.
! if gc is zero, then we just have to switch row (it1) and row (it1+1)
! of r and column (it1) and column (it1+1) of j. since we swithc rows in
! r and columns in j, we can ignore the sign of gs.
! otherwise we have to apply the givens rotation to these rows/columns.

   if (gc == 1.d0) go to 798
   if (gc == 0.d0) then
      do i = it1 + 1, nact
         temp = work(l1 - 1)
         work(l1 - 1) = work(l1)
         work(l1) = temp
         l1 = l1 + i
      end do
      do i = 1, n
         temp = dmat(i, it1)
         dmat(i, it1) = dmat(i, it1 + 1)
         dmat(i, it1 + 1) = temp
      end do
   else
      nu = gs/(1.d0 + gc)
      do i = it1 + 1, nact
         temp = gc*work(l1 - 1) + gs*work(l1)
         work(l1) = nu*(work(l1 - 1) + temp) - work(l1)
         work(l1 - 1) = temp
         l1 = l1 + i
      end do
      do i = 1, n
         temp = gc*dmat(i, it1) + gs*dmat(i, it1 + 1)
         dmat(i, it1 + 1) = nu*(dmat(i, it1) + temp) - dmat(i, it1 + 1)
         dmat(i, it1) = temp
      end do
   end if

! shift column (it1+1) of r to column (it1) (that is, the first it1
! elements). the posit1on of element (1,it1+1) of r was calculated above
! and stored in l.

798 continue
   l1 = l - it1
   do i = 1, it1
      work(l1) = work(l)
      l = l + 1
      l1 = l1 + 1
   end do

! update vector u and iact as necessary
! continue with updating the matrices j and r

   work(iwuv + it1) = work(iwuv + it1 + 1)
   iact(it1) = iact(it1 + 1)
   it1 = it1 + 1
   if (it1 < nact) go to 797
799 work(iwuv + nact) = work(iwuv + nact + 1)
   work(iwuv + nact + 1) = 0.d0
   iact(nact) = 0
   nact = nact - 1
   iter(2) = iter(2) + 1
   go to 55
999 continue
   return
   end procedure qpgen1

!     dpofa factors a double precision symmetric positive definite
!     matrix.
!     dpofa is usually called by dpoco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
!     on entry
!        a       double precision(lda, n)
!                the symmetric matrix to be factored.  only the
!                diagonal and upper triangle are used.
!        lda     integer
!                the leading dimension of the array  a .
!        n       integer
!                the order of the matrix  a .
!     on return
!        a       an upper triangular matrix  r  so that  a = trans(r)*r
!                where  trans(r)  is the transpose.
!                the strict lower triangle is unaltered.
!                if  info .ne. 0 , the factorization is not complete.
!        info    integer
!                = 0  for normal return.
!                = k  signals an error condition.  the leading minor
!                     of order  k  is not positive definite.
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
   module procedure dpofa
   real(dp) :: t, s
   integer  :: j, k

   do j = 1, n
      info = j; s = 0.0_dp
      do k = 1, j - 1
         t = A(k, j) - dot_product(A(:k - 1, k), A(:k - 1, j))
         t = t/A(k, k); A(k, j) = t; s = s + t*t
      end do
      s = A(j, j) - s
      if (s <= 0.0_dp) return
      A(j, j) = sqrt(s)
   end do
   info = 0
   return
   end procedure

!     dpori computes the inverse of the factor of a
!     double precision symmetric positive definite matrix
!     using the factors computed by dpofa.
!
!     modification of dpodi by bat 05/11/95
!
!     on entry
!
!        a       double precision(lda, n)
!                the output  a  from dpofa
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       if dpofa was used to factor  a  then
!                dpodi produces the upper half of inverse(a) .
!                elements of  a  below the diagonal are unchanged.
!
!     error condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if dpoco or dpofa has set info .eq. 0 .
!
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!     modified by berwin a. turlach 05/11/95
   module procedure dpori
   real(dp) :: t
   integer :: j, k

   !> Compute the inverse.
   do k = 1, n
      A(k, k) = 1.0d0/A(k, k); t = -A(k, k)
      A(:k - 1, k) = t*A(:k - 1, k)
      do j = k + 1, n
         t = A(k, j); A(k, j) = 0.0d0
         A(:k, j) = A(:k, j) + t*A(:k, k)
      end do
   end do
   return
   end procedure

!     dposl solves the double precision symmetric positive definite
!     system a * x = b
!     using the factors computed by dpoco or dpofa.
!     on entry
!        a       double precision(lda, n)
!                the output from dpoco or dpofa.
!        lda     integer
!                the leading dimension of the array  a .
!        n       integer
!                the order of the matrix  a .
!        b       double precision(n)
!                the right hand side vector.
!     on return
!        b       the solution vector  x .
!     error condition
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal.  technically this indicates
!        singularity but it is usually caused by improper subroutine
!        arguments.  it will not occur if the subroutines are called
!        correctly and  info .eq. 0 .
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dpoco(a,lda,n,rcond,z,info)
!           if (rcond is too small .or. info .ne. 0) go to ...
!           do 10 j = 1, p
!              call dposl(a,lda,n,c(1,j))
!        10 continue
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
   pure subroutine dposl(A, lda, n, b)
!      use quadprog_constants, only: dp
      integer, intent(in)  :: lda
  !! Leading dimension of the matrix A (i.e. number of rows).
      integer, intent(in)  :: n
  !! Number of columns of the matrix A.
      real(dp), intent(in) :: A(lda, *)
  !! Matrix A already factorized.
      real(dp), intent(inout) :: b(*)
  !! On entry, right-hand side vector.
  !! On exit, solution vector.

      !> Internal variables.
      real(dp) :: t
      integer  :: k, kb

      ! Solve trans(r)*y = b
      do k = 1, n
         t = dot_product(A(:k - 1, k), b(:k - 1))
         b(k) = (b(k) - t)/A(k, k)
      end do

      ! Solve r*x = y
      do kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/A(k, k)
         t = -b(k)
         b(:k - 1) = t*A(:k - 1, k) + b(:k - 1)
      end do
      return
   end subroutine dposl
end submodule
