module acfd_quadrature
  private
  public :: quadrature_driver

contains

  ! Driver for generating quadrature points and weights
  !
  ! \param [out] omega quadrature points
  ! \param [out] weight quadrature weights
  ! \param [in] nquad number of points
  ! \param [in] eig array ?
  ! \param [in] no
  ! \param [in] nv
  subroutine quadrature_driver(omega, weight, nquad)
    use acfd_parameters, only: quad_w0
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: nquad
    real(dp), intent(out) :: omega(nquad), weight(nquad)

    call quad_ratio(omega, weight, nquad, quad_w0)

  end subroutine quadrature_driver

  ! Gauss-Legendre quadrature algorithm using a rational function
  ! for distributing supporting nodes on the interval [a; +inf)
  !
  ! \param [out] omega quadrature points
  ! \param [out] weight quadrature weights
  ! \param [in] nquad number of quadrature points
  ! \param [in] w0 scaling factor for quadrature points
  ! \param [in] a starting point of integration interval
  subroutine quad_rationaltail(omega, weight, nquad, w0, a)
    use memory
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: nquad
    real(dp), intent(in) :: w0, a
    real(dp), intent(out) :: omega(nquad), weight(nquad)

    real(dp), pointer, dimension(:) :: omega_raw, weight_raw
    integer :: i, j

    omega = 0d0
    weight = 0d0

    omega_raw => memory_allocate(nquad)
    omega_raw = 0d0
    weight_raw => memory_allocate(nquad)
    weight_raw = 0d0
    call Gauss_Legendre_points(-1d0, 1d0, omega_raw, weight_raw, nquad)

    j = nquad
    do i = 1, nquad
      omega(i) = w0*(1d0-omega_raw(j)+a)/(1d0+omega_raw(j))
      weight(i) = w0*weight_raw(j)*(a+2d0)/((1d0+omega_raw(j))**2)
      j = j-1
    end do
    call memory_release(omega_raw)
  end subroutine quad_rationaltail

  ! Gauss-Legendre quadrature algorithm using a rational function
  ! for distributing supporting nodes on the interval [a; +inf)
  !
  ! \param [out] omega quadrature points
  ! \param [out] weight quadrature weights
  ! \param [in] nquad number of quadrature points
  ! \param [in] w0 scaling factor for quadrature points
  ! \param [in] a starting point of integration interval
  subroutine quad_rationaltail2(omega, weight, nquad, w0, a)
    use memory
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: nquad
    real(dp), intent(in) :: w0, a
    real(dp), intent(out) :: omega(nquad), weight(nquad)

    real(dp), pointer, dimension(:) :: omega_raw, weight_raw
    integer :: i, j

    omega = 0d0
    weight = 0d0

    omega_raw => memory_allocate(nquad)
    omega_raw = 0d0
    weight_raw => memory_allocate(nquad)
    weight_raw = 0d0
    call Gauss_Legendre_points(-1d0, 1d0, omega_raw, weight_raw, nquad)

    j = nquad
    do i = 1, nquad
      omega(i) = w0*(1d0-omega_raw(j))/(1d0+omega_raw(j))+a
      weight(i) = w0*weight_raw(j)*(2d0)/((1d0+omega_raw(j))**2)
      j = j-1
    end do
    call memory_release(omega_raw)
  end subroutine quad_rationaltail2

  ! Gauss-Legendre quadrature algorithm using a rational function
  ! for distributing supporting nodes on the interval [0; +inf)
  !
  ! \param [out] omega quadrature points
  ! \param [out] weight quadrature weights
  ! \param [in] nquad number of quadrature points
  ! \param [in] w0 scaling factor for quadrature points
  subroutine quad_ratio(omega, weight, nquad, w0)
    use acfd_parameters, only: nquadint
    use iso_fortran_env, only: dp => real64
    implicit none

    !integer,parameter :: n_int =5
    integer :: i

    real(dp) :: a(nquadint)!=0d0

    integer, intent(in) :: nquad
    real(dp), intent(in) :: w0
    real(dp), intent(out) :: omega(nquad), weight(nquad)

    if (nquadint .le. 0) then
      write (*, *) "something nasty is going on"
    end if

    do i = 1, nquadint
      a(i) = 0d0
    end do

    if (nquadint .eq. 1) then
      call quad_rationaltail(omega, weight, nquad, w0, a(nquadint))
    else
      a(2) = 10.0d0**(-(nquadint-1))

      do i = 1, nquadint-1
        call Gauss_Legendre_points(a(i), a(i+1), &
                                   omega(((i-1)*nquad/nquadint+1):(i*nquad/nquadint)), &
                                   weight(((i-1)*nquad/nquadint+1):(i*nquad/nquadint)), &
                                   nquad/nquadint)

        a(i+2) = a(i+1)*10.0d0

      end do

      call quad_rationaltail2(omega((nquadint-1)*nquad/nquadint+1:nquadint*nquad/nquadint), &
                weight((nquadint-1)*nquad/nquadint+1:nquadint*nquad/nquadint), &
                              nquad/nquadint, w0, a(nquadint))

    end if

  end subroutine quad_ratio

end module acfd_quadrature
