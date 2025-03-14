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
subroutine dpofa(a, lda, n, info)
   use quadprog_constants, only: dp
   integer, intent(in) :: lda
  !! Leading dimension of the matrix A (i.e. number of rows).
   integer, intent(in) :: n
  !! Number of columns of the matrix A.
   real(dp), intent(inout) :: A(lda, *)
  !! Matrix to be factorized.
   integer, intent(out) :: info
  !! Information flag.

   ! internal variables.
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
end subroutine dpofa
