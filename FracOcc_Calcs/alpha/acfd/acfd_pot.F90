!
! Module containing potentials and fock matrices
!
module acfd_pot
  use iso_fortran_env, only: dp => real64
  implicit none

  real(dp), pointer :: fock(:), old_fock(:) ! fock matrix and previous fock matrix
  real(dp), pointer :: vx(:) ! exchange potential
  real(dp), pointer :: vsol(:) ! solution vector, i.e, rhs in unprocessed basis for EXX-OEP

  real(dp), pointer :: focka(:), fockb(:), old_focka(:), old_fockb(:)
  real(dp), pointer :: vxa(:), vxb(:), vxb1(:), vxb2(:)
  real(dp), pointer :: vsola(:), vsolb(:)

  real(dp), pointer :: vj(:) ! Coulomb potential
  real(dp), pointer :: h0(:) ! one electron potential

end module acfd_pot
