module acfd_plot
  use acfd_derived_datatypes, only: vector_type
  use iso_fortran_env, only: dp => real64
  implicit none

  type(vector_type) :: vref_aux_upa, vref_aux_upb, &
                       vref_aux_up, vref_aux_dummy ! reference potential, saved for plotting
  integer, dimension(3), parameter :: ngrid = (/71, 71, 71/) ! number of grid points per dimension
  real(dp), dimension(3), parameter :: dist = (/10d0, 10d0, 10d0/) ! size of cube
  real(dp), dimension(3), parameter :: &
    origin = (/-dist(1)/2d0, -dist(2)/2d0, -dist(3)/2d0/) ! origin of cube (lower edge)

contains

  ! plots exchange and correlation potentials along the x,y and z axis
  ! CAUTON: only works if a hf calculation has been done before, something with the initialization of the grid
  ! output format for 1d plot: x,y,z,r,pot,ref,pot+ref (pot+ref is nonsense for vc plots, as is ref)
  subroutine plot_potential(orb, den, oep_solution, naux, itc, basis, na, &
                            potfil, vref_aux_unproc)
    use common_cbas, only: ntdg, ntqg
    use molpro_options, only: molpro_pwd
    use acfd_derived_datatypes, only: vector_type
    use acfd_parameters, only: l_plot_x_line, l_plot_y_line, l_plot_z_line
    use iso_fortran_env, only: dp => real64
    implicit none
    character(len=*), intent(in) :: itc ! iteration as character
    character(len=*), intent(in) :: basis ! name of auxiliary basis set (usually OEP)
    character(len=*), intent(in) :: potfil ! name of the plot file
    integer, intent(in) :: naux ! dimension of oep_solution
    integer, intent(in) :: na ! number of occupied orbitals
    real(dp), intent(in)  :: orb(ntqg) ! orbitals
    real(dp), intent(in)  :: den(ntdg) ! density
    real(dp), intent(in)  :: oep_solution(naux) ! solution vector of the OEP equation in unprocessed OEP basis
    type(vector_type) :: vref_aux_unproc ! reference potential in the unprocessed auxiliary basis set
    real(dp) :: r1(3), r2(3)
    real(dp), pointer :: dum(:, :)

    call coulomb_basis(basis, 0)

    allocate (dum(naux, 0))

    if (l_plot_x_line) then
      r1(:) = (/-1d0, 0d0, 0d0/)
      r2(:) = (/1d0, 0d0, 0d0/)
      call oep_potential(orb, den, oep_solution, vref_aux_unproc%mat, dum, &
                         r1, r2, na, naux, 0, basis, 'BASIS2', &
                         molpro_pwd//trim(potfil)//'-'//itc//'.x', 'J')
    end if

    if (l_plot_y_line) then
      r1(:) = (/0d0, -1d0, 0d0/)
      r2(:) = (/0d0, 1d0, 0d0/)
      call oep_potential(orb, den, oep_solution, vref_aux_unproc%mat, dum, &
                         r1, r2, na, naux, 0, basis, 'BASIS2', &
                         molpro_pwd//trim(potfil)//'-'//itc//'.y', 'J')
    end if

    if (l_plot_z_line) then
      r1(:) = (/0d0, 0d0, -1d0/)
      r2(:) = (/0d0, 0d0, 1d0/)

      call oep_potential(orb, den, oep_solution, vref_aux_unproc%mat, dum, &
                         r1, r2, na, naux, 0, basis, 'BASIS2', &
                         molpro_pwd//trim(potfil)//'-'//itc//'.z', 'J')
    end if

  end subroutine plot_potential

end module acfd_plot
