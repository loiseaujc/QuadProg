module markowitz_utilities
   use QuadProg
   implicit none
   private

   public :: dp
   public :: equicorr, efficient_frontier

contains

   !---------------------------------------------
   !-----     Assets' covariance matrix     -----
   !---------------------------------------------

   pure function equicorr(n, r) result(xmat)
      ! return a correlation matrix with equal off-diagonal elements
      integer, intent(in)  :: n ! dimension of correlation matrix
      real(kind=dp), intent(in)  :: r ! off-diagonal correlation
      real(kind=dp), allocatable :: xmat(:, :)
      integer                    :: i
      allocate (xmat(n, n), source=r)
      do i = 1, n
         xmat(i, i) = 1.0_dp
      end do
   end function equicorr

   !------------------------------------------------
   !-----     Markowitz efficient frontier     -----
   !------------------------------------------------

   function efficient_frontier(mu, cov, gamma) result(output)
      real(dp), intent(in) :: mu(:), cov(:, :)
      ! Assets' expected returns and covariance matrix.
      real(dp), intent(in) :: gamma(:)
      ! Set of risk-aversion parameters.
      real(dp), allocatable :: output(:, :)
      ! Expected returns and standard deviation of the efficient portfolios.

      !----- Internal variables -----
      integer :: i, n, npts
      real(dp), allocatable :: A(:, :), b(:), C(:, :), d(:)
      type(qp_problem) :: problem
      type(OptimizeResult) :: solution

      !> Allocate variables.
      npts = size(gamma); allocate (output(npts, 2), source=0.0_dp)

      !> Long-only constraints.
      n = size(mu); allocate (C(n, n), d(n), source=0.0_dp)
      do i = 1, n
         C(i, i) = 1.0_dp
      end do

      !> Full allocation constraint.
      allocate (A(1, n), b(1), source=1.0_dp)

      !> Initial QP.
      problem = qp_problem(cov, mu, A=A, b=b, C=C, d=d)

      !> Compute efficient frontier.
      do i = 1, npts
         ! Update right-hand side vector.
         problem%q = (1.0_dp/gamma(i))*mu
         ! Solve QP.
         solution = solve(problem)
         ! Compute expected return and standard deviation.
         associate (x => solution%x)
            if (solution%success) then
               ! Expected return.
               output(i, 1) = dot_product(mu, x)
               ! Standard deviation.
               output(i, 2) = sqrt(dot_product(x, matmul(cov, x)))
            else
               error stop
            end if
         end associate
      end do

   end function efficient_frontier

end module markowitz_utilities

program markowitz_example
   use markowitz_utilities
   implicit none

   integer, parameter :: n = 3, npts = 101
   integer :: i, io
   real(dp) :: r, cov(n, n), mu(n), portfolio(n)
   real(dp) :: gamma(npts), exponents(npts), output(npts, 2)
   real(dp) :: expected_return, standard_deviation

   !> Expected returns of the different assets.
   mu = 1.0_dp; mu(1) = 2.0_dp ! Stock n°1 has higher expected return.

   !> Range of risk-aversion parameters.
   exponents = [(-4.0_dp + 0.08_dp*i, i=1, npts)]
   gamma = 10.0_dp**exponents

   !---------------------------------------
   !-----     Uncorrelated assets     -----
   !---------------------------------------

   !> Covariance matrix.
   r = 0.0_dp; cov = equicorr(n, r)

   !> Markowitz efficient frontier.
   output = efficient_frontier(mu, cov, gamma)

   !-----------------------------------------
   !-----     0.5-correlated assets     -----
   !-----------------------------------------

   !> Covariance matrix.
   r = 0.5_dp; cov = equicorr(n, r)

   !> Markowitz efficient frontier.
   output = efficient_frontier(mu, cov, gamma)
   !-----------------------------------------
   !-----     0.8-correlated assets     -----
   !-----------------------------------------

   !> Covariance matrix.
   r = 0.8_dp; cov = equicorr(n, r)

   !> Markowitz efficient frontier.
   output = efficient_frontier(mu, cov, gamma)

end program markowitz_example
