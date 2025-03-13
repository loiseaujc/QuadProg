module QuadProg
  implicit none
  private

  public :: qpgen2

  interface
    subroutine qpgen2(dmat, dvec, fddmat, n, sol, lagr, crval, amat, bvec, fdamat, q, meq, iact, nact, iter, work, ierr)
      integer :: fddmat, n, fdamat, q, meq, iact(*), nact, iter(*), ierr
      double precision :: dmat(fddmat, *), dvec(*), lagr(*), sol(*), bvec(*), work(*), amat(fdamat, *), crval
    end subroutine
  end interface

end module QuadProg
