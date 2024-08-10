!
! contains transformation matrices and their dimensions
!
module acfd_auxiliary_basis
  use acfd_derived_datatypes
  use iso_fortran_env, only: dp => real64
  implicit none

  integer :: naux_overlap ! number of basis functions after reduction (applying trans_mat_overlap)
  type(matrix_type) :: trans_mat_overlap ! transformation matrix for throwing out functions with small overlap

  integer :: naux_constraint ! number of basis functions after reduction and constraints
  integer :: naux_charge ! number of basis functions after reduction and charge constraint
  type(matrix_type) :: trans_mat_constraint ! transformation matrix for reduction and constraints
  type(matrix_type) :: trans_mat_charge ! transformation matrix for reduction and charge constraint

  type(matrix_type) :: trans_mat_charge_ri ! transformation matrix for reduction and charge constraint
  integer :: naux_ri ! number of basis functions in RI basis set
  integer :: naux_charge_ri ! number of remaining basis functions after reduction and charge constraint

  integer :: naux_constraint_a, naux_constraint_b ! number of basis functions after reduction and constraints
  type(matrix_type) :: trans_mat_constraint_a ! transformation matrix for reduction and constraints
  type(matrix_type) :: trans_mat_constraint_b ! transformation matrix for reduction and constraints

  integer :: naux_oep ! number of OEP basis functions

end module acfd_auxiliary_basis
