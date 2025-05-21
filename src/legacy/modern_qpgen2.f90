subroutine modern_qpgen2(dmat, dvec, fddmat, n, sol, lagr, crval, &
                         amat, bvec, fdamat, q, meq, iact, nact, iter, work, ierr)

   implicit none
   double precision, intent(in out)         :: dmat(fddmat, *)
   double precision, intent(in out)         :: dvec(*)
   integer, intent(in out)                  :: fddmat
   integer, intent(in)                      :: n
   double precision, intent(out)            :: sol(*)
   double precision, intent(out)            :: lagr(*)
   double precision, intent(out)            :: crval
   double precision, intent(in out)         :: amat(fdamat, *)
   double precision, intent(in out)         :: bvec(*)
   integer, intent(in out)                  :: fdamat
   integer, intent(in)                      :: q
   integer, intent(in)                      :: meq
   integer, intent(out)                     :: iact(*)
   integer, intent(out)                     :: nact
   integer, intent(out)                     :: iter(*)
   double precision, intent(out)            :: work(*)
   integer, intent(in out)                  :: ierr
   integer :: i, j, l, l1, info, it1, &
              iwzv, iwrv, iwrm, iwsv, iwuv, nvl, r, iwnbv
   double precision :: temp, sum, t1, tt, gc, gs, nu, vsmall, tmpa, tmpb
   logical :: t1inf, t2min

   r = min(n, q)
   l = 2*n + (r*(r + 5))/2 + 2*q + 1

!     code gleaned from powell's zqpcvx routine to determine a small
!     number  that can be assumed to be an upper bound on the relative
!     precision of the computer arithmetic.

   vsmall = 1.0d-60
1  vsmall = vsmall + vsmall
   tmpa = 1.0d0 + 0.1d0*vsmall
   tmpb = 1.0d0 + 0.2d0*vsmall
   if (tmpa <= 1.0d0) go to 1
   if (tmpb <= 1.0d0) go to 1

! store the initial dvec to calculate below the unconstrained minima of
! the critical value.

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

! get the initial solution

   if (ierr == 0) then
      call dpofa(dmat, fddmat, n, info)
      if (info /= 0) then
         ierr = 2
         go to 999
      end if
      call dposl(dmat, fddmat, n, dvec)
      call dpori(dmat, fddmat, n)
   else

! matrix d is already factorized, so we have to multiply d first with
! r^-t and then with r^-1.  r^-1 is stored in the upper half of the
! array dmat.

      do j = 1, n
         sol(j) = 0.d0
         do i = 1, j
            sol(j) = sol(j) + dmat(i, j)*dvec(i)
         end do
      end do
      do j = 1, n
         dvec(j) = 0.d0
         do i = j, n
            dvec(j) = dvec(j) + dmat(j, i)*sol(i)
         end do
      end do
   end if

! set lower triangular of dmat to zero, store dvec in sol and
! calculate value of the criterion at unconstrained minima

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

! calculate some constants, i.e., from which index on the different
! quantities are stored in the work matrix

   iwzv = n
   iwrv = iwzv + n
   iwuv = iwrv + r
   iwrm = iwuv + r + 1
   iwsv = iwrm + (r*(r + 1))/2
   iwnbv = iwsv + q

! calculate the norm of each column of the a matrix

   do i = 1, q
      sum = 0.d0
      do j = 1, n
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
      do j = 1, n
         sum = sum + amat(j, i)*sol(j)
      end do
      if (abs(sum) < vsmall) then
         sum = 0.0d0
      end if
      if (i > meq) then
         work(l) = sum
      else
         work(l) = -abs(sum)
         if (sum > 0.d0) then
            do j = 1, n
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
!         if (work(iwsv+i) .lt. 0.d0) then
!            nvl = i
!            goto 72
!         endif
   end do
72 if (nvl == 0) then
      do i = 1, nact
         lagr(iact(i)) = work(iwuv + i)
      end do
      go to 999
   end if

! calculate d=j^tn^+ where n^+ is the normal vector of the violated
! constraint. j is stored in dmat in this implementation!!
! if we drop a constraint, we have to jump back here.

55 continue
   do i = 1, n
      sum = 0.d0
      do j = 1, n
         sum = sum + dmat(j, i)*amat(j, nvl)
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
7     t1inf = .false.
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

         ierr = 1
         go to 999
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
         do j = 1, n
            sum = sum + sol(j)*amat(j, nvl)
         end do
         if (nvl > meq) then
            work(iwsv + nvl) = sum
         else
            work(iwsv + nvl) = -abs(sum)
            if (sum > 0.d0) then
               do j = 1, n
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
eND SUBROUTINE modern_qpgen2
