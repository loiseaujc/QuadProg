module quadprog_benchmark
   use stdlib_io_npy, only: load_npy
   use quadprog, only: qp_problem, dp
   implicit none
   private

   character(len=255), dimension(19), parameter, public :: problems = [character(len=255) :: &
                                                                       "test/data/HS21", &
                                                                       "test/data/QPTEST", &
                                                                       "test/data/HS35MOD", &
                                                                       "test/data/HS35", &
                                                                       "test/data/HS76", &
                                                                       "test/data/HS268", &
                                                                       "test/data/S268", &
                                                                       "test/data/DUALC5", &
                                                                       "test/data/DUALC1", &
                                                                       "test/data/HS118", &
                                                                       "test/data/KSIP", &
                                                                       "test/data/DUAL4", &
                                                                       "test/data/QPCBLEND", &
                                                                       "test/data/DUAL1", &
                                                                       "test/data/DUAL2", &
                                                                       "test/data/DUAL3", &
                                                                       "test/data/QPCSTAIR", &
                                                                       "test/data/MOSARQP2", &
                                                                       "test/data/LASER" &
                                                                       ]

   public :: get_problem_data, problem_description
contains

   integer function problem_description(problem, pname) result(ierr)
      type(qp_problem), intent(in) :: problem
      character(len=*), intent(in) :: pname

      print *
      print *, "--------------------"
      print *, "     * Problem ID :                                    ", pname(6:)
      print *, "         - Number of variables              : ", size(problem%P, 1)
      print *, "         - Number of equality constraints   : ", problem%neq
      print *, "         - Number of inequality constraints : ", problem%ncons - problem%neq
      print *
      return
   end function

   type(qp_problem) function get_problem_data(fname) result(problem)
      character(len=*), intent(in) :: fname
      real(dp), allocatable :: P(:, :), q(:)
      real(dp), allocatable :: A(:, :), b(:)
      real(dp), allocatable :: C(:, :), d(:)
      integer :: i

      !> Load P and q for the quadratic form.
      call load_npy(trim(fname)//"/P_matrix.npy", P)
      call load_npy(trim(fname)//"/q_vector.npy", q)

      !> Load A and b for the equality constraints.
      call load_npy(trim(fname)//"/A_matrix.npy", A)
      call load_npy(trim(fname)//"/b_vector.npy", b)

      !> Load A and b for the inequality constraints.
      call load_npy(trim(fname)//"/C_matrix.npy", C)
      call load_npy(trim(fname)//"/d_vector.npy", d)

      if (size(A, 2) == 1) A = transpose(A)

      if ((size(A, 2) == 1) .and. (size(C, 2) == 1)) then
         problem = qp_problem(P, -q)
      else if ((size(A, 2) == 1) .and. (size(C, 2) /= 1)) then
         problem = qp_problem(P, -q, C=C, d=d)
      else if ((size(A, 2) /= 1) .and. (size(C, 2) == 1)) then
         problem = qp_problem(P, -q, A=A, b=b)
      else
         problem = qp_problem(P, -q, A=A, b=b, C=C, d=d)
      end if

      return
   end function
end module quadprog_benchmark
