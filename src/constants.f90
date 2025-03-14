module quadprog_constants
   implicit none

   integer, parameter, public :: dp = selected_real_kind(15, 307)
   real(dp), parameter, public :: atol = 10.0_dp**(-precision(1.0_dp))
   real(dp), parameter, public :: rtol = sqrt(atol)
   real(dp), parameter, public :: eps = epsilon(1.0_dp)

end module
