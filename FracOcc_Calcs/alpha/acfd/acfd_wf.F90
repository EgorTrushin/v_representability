!
! module containing orbitals, eigenvalues and densities
!
module acfd_wf
  use iso_fortran_env, only: dp => real64
  implicit none

  integer :: no, nv ! number of occupied and virtual orbitals
  real(dp), pointer :: orb(:), eig(:), den(:), occ(:) ! orbitals, eigenvalues, denisty, occupations

  integer :: noa, nob, nva, nvb
  real(dp), pointer :: orba(:), orbb(:), orbb_frac(:)
  real(dp), pointer :: eiga(:), eigb(:)
  real(dp), pointer :: dena(:), denb(:), denb_frac(:)
  real(dp), pointer :: occa(:), occb(:), occb1(:), occb2(:)

  real(dp), pointer :: smh(:), sh_a(:), int_all(:, :, :)

end module acfd_wf
