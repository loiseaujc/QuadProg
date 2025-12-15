submodule(quadprog) linalg
   use quadprog_constants, only: dp
   implicit none
contains

   module procedure cho_solve
   integer :: n, info
   real(dp), pointer :: x(:, :)
   n = size(A, 1)
   x(1:n, 1:1) => b
   call potrs("u", n, 1, A, n, x, n, info)
   end procedure cho_solve

end submodule linalg
