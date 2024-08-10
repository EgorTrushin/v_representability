!
! Derived datatypes for passing not allocated matrices and vectors as arguments.
!
module acfd_derived_datatypes
  use iso_fortran_env, only: dp => real64
  implicit none

  private

  type, public :: matrix_type
    real(dp), allocatable :: mat(:, :)
    integer :: n, m

  contains

    procedure :: init => init_mat
    procedure :: free => free_mat

  end type matrix_type

  type, public :: vector_type
    real(dp), allocatable :: mat(:)
    integer :: n

  contains

    procedure :: init => init_vec
    procedure :: free => free_vec

  end type vector_type

contains

  subroutine init_mat(mat, n, m)
    implicit none
    class(matrix_type) :: mat
    integer :: n, m

    if (allocated(mat%mat)) deallocate (mat%mat)
    allocate (mat%mat(n, m))
    mat%mat = 0d0
    mat%n = n
    mat%m = m
  end subroutine init_mat

  subroutine init_vec(mat, n)
    implicit none
    class(vector_type) :: mat
    integer :: n

    if (allocated(mat%mat)) deallocate (mat%mat)
    allocate (mat%mat(n))
    mat%mat = 0d0
    mat%n = n
  end subroutine init_vec

  subroutine free_vec(mat)
    implicit none
    class(vector_type) :: mat
    if (allocated(mat%mat)) deallocate (mat%mat)
  end subroutine free_vec

  subroutine free_mat(mat)
    implicit none
    class(matrix_type) :: mat
    if (allocated(mat%mat)) deallocate (mat%mat)
  end subroutine free_mat

end module acfd_derived_datatypes
