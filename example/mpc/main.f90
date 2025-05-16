module mpc_utilities
   use QuadProg
   implicit none
   private

   public :: dp
   public :: lti

   !---------------------------------------------------------------------------
   !-----     Derived-type for representing discrete-time LTI systems     -----
   !---------------------------------------------------------------------------

   type, public :: lti_system
      real(dp), allocatable :: A(:, :)
      real(dp), allocatable :: B(:, :)
      real(dp), allocatable :: C(:, :)
      real(dp), allocatable :: D(:, :)
   end type
   interface lti
      module procedure initialize_system
   end interface

   !----------------------------------------------------------------
   !-----     Derived-type for representing MPC controller     -----
   !----------------------------------------------------------------

   type, public :: mpc_controller
      type(lti_system)      :: plant
      !! Plant to be controlled.
      integer               :: horizon
      !! Horizon of the MPC problem.
      type(qp_problem) :: prob
      real(dp), allocatable :: P(:, :), q(:, :), x0(:)
      !! Condensed quadratic cost function + initial condition
      real(dp), allocatable :: C(:, :), d(:)
      !! Inequality constraints.
   contains
      procedure, pass(self), public :: compute
   end type
   interface mpc_controller
      module procedure initialize_mpc_controller
   end interface

contains

   !-------------------------------------
   !-----     Type constructors     -----
   !-------------------------------------

   pure type(lti_system) function initialize_system(A, B, C, D) result(plant)
      real(dp), intent(in) :: A(:, :)
      real(dp), intent(in) :: B(:, :)
      real(dp), intent(in) :: C(:, :)
      real(dp), optional, intent(in) :: D(:, :)
      integer :: nstates, ninputs, noutputs

      !> Problem's dimensions.
      nstates = size(A, 1); ninputs = size(B, 2); noutputs = size(C, 1)

      !> Sanity checks.
      if (size(A, 1) /= size(A, 2)) error stop "Dynamics matrix A is not square."
      if (size(B, 1) /= size(A, 1)) error stop "Input matrix B and dynamics matrix A have incompatible dimensions."
      if (size(C, 2) /= size(A, 1)) error stop "Measurement matrix C and dynamics matrix A have incompatible dimensions."

      !> Initialize the LTI system.
      plant%A = A; plant%B = B; plant%C = C
      if (present(D)) then
         plant%D = D
      else
         allocate (plant%D(noutputs, ninputs), source=0.0_dp)
      end if
      return
   end function

   type(mpc_controller) function initialize_mpc_controller(plant, horizon, Q, R, u_ub, u_lb) result(controller)
      type(lti_system), intent(in) :: plant
      !! Plant to be controller.
      integer, intent(in)  :: horizon
      !! Horizon for the MPC optimizer.
      real(dp), intent(in) :: Q(:, :)
      !! State cost (assumed constant over time).
      real(dp), intent(in) :: R(:, :)
      !! Input cost (assumed constant over time).
      real(dp), optional, intent(in) :: u_ub(:), u_lb(:)
      !! Upper and lower-bounds for the input (assumed constant over time.)

      integer :: nstates, ninputs ! Number of states and inputs
      integer :: i, j, ncons
      real(dp), allocatable :: blkdiag_Q(:, :), blkdiag_R(:, :)
      real(dp), allocatable :: O(:, :), blktoeplitz_B(:, :)

      !> Plant's dimensions.
      nstates = size(plant%A, 1); ninputs = size(plant%B, 2)
      controller%plant = plant

      !> Sanity checks for the state cost.
      if (size(Q, 1) /= size(Q, 2)) error stop "Matrix Q is not square."
      if (size(Q, 1) /= nstates) error stop "Dimensions of Q are inconsistent with the state's dimension."

      !> Santity checks for the input cost.
      if (size(R, 1) /= size(R, 2)) error stop "Matrix R is not square."
      if (size(R, 1) /= ninputs) error stop "Dimensions of R are inconsistent with the number of inputs."

      !> Sanity checks for the input's box constraints.
      if (present(u_ub)) then
         if (size(u_ub) /= ninputs*horizon) error stop "Dimension of u_ub is inconsistent with the number of inputs."
      end if
      if (present(u_lb)) then
         if (size(u_lb) /= ninputs*horizon) error stop "Dimension of u_lb is inconsistent with the number of inputs."
      end if
      if (present(u_lb) .and. present(u_ub)) then
         do i = 1, size(u_lb)
            if (u_lb(i) > u_ub(i)) error stop "Upper and lower bounds are incompatible."
         end do
      end if

      !> Sanity check for the MPC optimization horizon.
      if (horizon <= 0) error stop "Horizon needs to be strictly positive."

      !------------------------------------------------------
      !-----     Construction of the quadratic cost     -----
      !------------------------------------------------------

      blkdiag_Q = blkdiag(Q, horizon)
      blkdiag_R = blkdiag(R, horizon)
      O = obs(plant%A, horizon)
      blktoeplitz_B = blktoeplitz(plant%A, plant%B, horizon)

      controller%P = matmul(matmul(transpose(blktoeplitz_B), blkdiag_Q), blktoeplitz_B) + blkdiag_R
      controller%q = matmul(matmul(transpose(blktoeplitz_B), blkdiag_Q), O)

      !---------------------------------------------------------------
      !-----     Construction of the input's box constraints     -----
      !---------------------------------------------------------------

      ncons = 0
      if (present(u_ub)) ncons = ncons + ninputs*horizon
      if (present(u_lb)) ncons = ncons + ninputs*horizon

      if (ncons /= 0) then
         allocate (controller%C(ncons, ninputs*horizon), source=0.0_dp)
         allocate (controller%d(ncons), source=0.0_dp)

         block
            integer :: i_start
            if (present(u_ub)) then
               i_start = 1
               do concurrent(i=i_start:ninputs*horizon)
                  controller%C(i, i) = -1.0_dp
                  controller%d(i) = -u_ub(i)
               end do
            end if

            if (present(u_lb)) then
               i_start = 1; if (present(u_ub)) i_start = i_start + ninputs*horizon
               do concurrent(i=i_start:(i_start - 1) + ninputs*horizon)
                  controller%C(i, i - i_start + 1) = 1.0_dp
                  controller%d(i) = u_lb(i - i_start + 1)
               end do
            end if
         end block
      end if

      !> Initialize the resulting QP problem.
      if (allocated(controller%C)) then
         controller%prob = qp_problem(controller%P, controller%q(:, 1), C=controller%C, d=controller%d)
      else
         controller%prob = qp_problem(controller%P, controller%q(:, 1))
      end if
   end function

   !-----------------------------------------
   !-----     Type-bound procedures     -----
   !-----------------------------------------

   function compute(self, x0) result(u)
      class(mpc_controller), intent(inout) :: self
      real(dp), intent(in) :: x0(:)
      real(dp) :: u(size(self%plant%B, 2))
      type(OptimizeResult) :: opt
      self%prob%q = -matmul(self%q, x0)
      opt = solve(self%prob)
      if (.not. opt%success) error stop "Quadratic solver failed to find a solution."
      u = opt%x(1)
      return
   end function

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   pure function matrix_power(A, pow) result(B)
      real(dp), intent(in) :: A(:, :)
      integer, intent(in)  :: pow
      real(dp), allocatable :: B(:, :)
      integer :: i, n

      !> Sanity checks.
      if (size(A, 1) /= size(A, 2)) error stop "Matrix A is not square."
      if (pow < 0) error stop "Power needs to be positive."

      !> Allocate matrix.
      n = size(A, 1); allocate (B, mold=A); B = 0.0_dp
      do concurrent(i=1:n)
         B(i, i) = 1.0_dp
      end do

      !> Treat the special case.
      if (pow == 0) return

      !> General case.
      do i = 1, pow
         B = matmul(A, B)
      end do
   end function

   pure function blkdiag(A, nblk) result(B)
      real(dp), intent(in) :: A(:, :)
      integer, intent(in)  :: nblk
      real(dp), allocatable :: B(:, :)
      integer :: i, n
      n = size(A, 1); allocate (B(n*nblk, n*nblk), source=0.0_dp)
      do concurrent(i=1:nblk)
         B((i - 1)*n + 1:i*n, (i - 1)*n + 1:i*n) = A
      end do
      return
   end function

   pure function obs(A, n) result(O)
      real(dp), intent(in) :: A(:, :)
      integer, intent(in) :: n
      real(dp), allocatable :: O(:, :)
      integer :: i, m
      m = size(A, 1); allocate (O(m*n, m), source=0.0_dp)

      do concurrent(i=1:n)
         O((i - 1)*m + 1:i*m, :) = matrix_power(A, i)
      end do
      return
   end function

   pure function blktoeplitz(A, B, n) result(T)
      real(dp), intent(in) :: A(:, :), B(:, :)
      integer, intent(in) :: n
      real(dp), allocatable :: T(:, :)
      integer :: nstate, ninput, i, j

      nstate = size(A, 1); ninput = size(B, 2)
      allocate (T(nstate*n, ninput*n), source=0.0_dp)

      do i = 1, n
         do j = 1, i
            T((i - 1)*nstate + 1:i*nstate, (j - 1)*ninput + 1:j*ninput) = matmul(matrix_power(A, i - j), B)
         end do
      end do
      return
   end function

end module

program mpc_examples
   use mpc_utilities
   implicit none
   integer :: i, j

   block
      !---------------------------------------------------
      !-----     Example n°1 : Double Integrator     -----
      !---------------------------------------------------

      real(dp) :: A(2, 2), B(2, 1), C(1, 2), D(1, 1)
      !! State-space model.
      type(lti_system) :: plant
      !! Corresponding LTI system.

      real(dp), parameter :: T = 25.0_dp, dt = 0.01_dp
      integer, parameter :: nt = int(T/dt)
      real(dp) :: x(2, 0:nt), y(1, 0:nt), u(1, 0:nt)

      integer, parameter :: horizon = int(5_dp/dt)
      !! Horizon for the MPC controller.

      real(dp) :: Q(2, 2), R(1, 1), cost
      !! State and input cost.
      real(dp), parameter :: u_ub(horizon) = 1_dp, u_lb(horizon) = -1_dp
      type(mpc_controller) :: controller
      !! Controller.
      real(dp) :: start_time, end_time

      ! Dynamics matrix.
      A(1, 1) = 1.0_dp; A(1, 2) = dt
      A(2, 1) = 0.0_dp; A(2, 2) = 1.0_dp

      ! Input-to-state matrix.
      B(1, 1) = 0.0_dp; B(2, 1) = dt

      ! State-to-output matrix.
      C(1, 1) = 1.0_dp; C(1, 2) = 0.0_dp

      ! Feedthrough matrix.
      D(1, 1) = 0.0_dp

      ! LTI system.
      plant = lti(A, B, C, D)

      ! Controller.
      Q = 0.0_dp; Q(1, 1) = 1.0_dp; Q(2, 2) = 1.0_dp; Q = Q*dt
      R = 0.0_dp; R(1, 1) = 1_dp; R = R*dt
      controller = mpc_controller(plant=plant, horizon=horizon, Q=Q, R=R, u_ub=u_ub, u_lb=u_lb)

      ! Initial condition.
      u = 0.0_dp; x(:, 0) = [-5.0_dp, 0.0_dp]; y(:, 0) = matmul(plant%C, x(:, 0))

      do i = 1, nt
         ! Compute control signal.
         u(:, i - 1) = controller%compute(x(:, i - 1))
         ! Integrate system.
         x(:, i) = matmul(plant%A, x(:, i - 1)) + matmul(plant%B, u(:, i - 1))
         ! Take measurement.
         y(:, i) = matmul(plant%C, x(:, i)) + matmul(plant%D, u(:, i))
      end do

      cost = 0.0_dp
      do i = 0, nt - 1
         cost = cost + 0.5_dp*dot_product(x(:, i), matmul(Q, x(:, i)))
         cost = cost + 0.5_dp*dot_product(u(:, i), matmul(R, u(:, i)))
      end do
      cost = cost + dot_product(x(:, nt), matmul(Q, x(:, nt)))

      print *, "Total cost :", cost

      open (unit=1234, file="example/mpc/mpc_response.dat")
      do i = 0, nt
         write (1234, *) i*dt, x(1, i), x(2, i), u(1, i)
      end do
      close (1234)

   end block

end program mpc_examples
