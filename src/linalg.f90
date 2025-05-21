submodule(quadprog) linalg
   use quadprog_constants, only: dp
contains

   module procedure qr
   integer :: m, n, i, j
   real(dp), allocatable :: tau(:), work(:)
   integer :: lwork, info

   !> Allocate matrices.
   m = size(A, 1); n = size(A, 2)
   Q = A; allocate (R(n, n)); R = 0.0_dp
   allocate (tau(n)); tau = 0.0_dp

   !> Workspace query.
   allocate (work(1)); lwork = -1
   call dgeqrf(m, n, Q, m, tau, work, lwork, info)

   !> QR factorization.
   lwork = int(work(1)); deallocate (work); allocate (work(lwork))
   call dgeqrf(m, n, Q, m, tau, work, lwork, info)

   !> Extract the R matrix.
   do concurrent(i=1:n, j=1:n)
      if (j >= i) R(i, j) = Q(i, j)
   end do

   !> Workspace query for reconstructing Q.
   lwork = -1; deallocate (work); allocate (work(1))
   call dorgqr(m, n, n, Q, m, tau, work, lwork, info)

   !> Extract Q matrix.
   lwork = work(1); deallocate (work); allocate (work(lwork))
   call dorgqr(m, n, n, Q, m, tau, work, lwork, info)
   end procedure

end submodule
