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
subroutine dpori(A, lda, n)
   use quadprog_constants, only: dp
   integer, intent(in) :: lda
  !! Leading dimension of A.
   real(dp), intent(out) :: A(lda, *)
  !! Matrix A already factorized.
   integer, intent(in) :: n
  !! Number of columns of A.

   !> Internal variables.
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
end subroutine dpori

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
subroutine dposl(A, lda, n, b)
   use quadprog_constants, only: dp
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
