module ksinv_plot
  implicit none
  private
  public :: plot_driver
contains
  subroutine plot_driver(naux, iter_tag, type_tag, vxc_sol, vx_sol, ref_den, &
                         den, v_ref)
    use common_cbas, only: ntg, ntdg
    use memory, only: memory_allocate, memory_release
    use ksinv_mod, only: l_plot_z_line, l_plot_y_line, l_plot_x_line, &
                         i_gridsize, r0_plot_range
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: naux
    real(dp), intent(in) :: v_ref(naux), vxc_sol(naux), vx_sol(naux)
    real(dp), intent(in) :: ref_den(ntdg), den(ntdg)
    character(len=*) :: iter_tag, type_tag
    character(5) :: grid_tag
    real(dp)  :: r1(3), r2(3)
    real(dp), pointer  ::  grid(:, :)

    grid => memory_allocate(3, i_gridsize)

    if (l_plot_z_line) then
      ! make grid z
      grid_tag = 'z'
      r1(:) = (/0d0, 0d0, -1d0/)
      r2(:) = (/0d0, 0d0, 1d0/)
      call generate_line_grid(r1, r2, i_gridsize, r0_plot_range, grid)
      call plot_subdriver(naux, ntdg, iter_tag, type_tag, grid_tag, &
                          vxc_sol, vx_sol, ref_den, den, v_ref, &
                          ntg, grid, &
                          i_gridsize)
    end if

    if (l_plot_y_line) then
      ! make grid y
      grid_tag = 'y'
      r1(:) = (/0d0, -1d0, 0d0/)
      r2(:) = (/0d0, 1d0, 0d0/)
      call generate_line_grid(r1, r2, i_gridsize, r0_plot_range, grid)
      call plot_subdriver(naux, ntdg, iter_tag, type_tag, grid_tag, &
                          vxc_sol, vx_sol, ref_den, den, v_ref, &
                          ntg, grid, &
                          i_gridsize)
    end if

    if (l_plot_x_line) then
      ! make grid x
      grid_tag = 'x'
      r1(:) = (/-1d0, 0d0, 0d0/)
      r2(:) = (/1d0, 0d0, 0d0/)
      call generate_line_grid(r1, r2, i_gridsize, r0_plot_range, grid)
      call plot_subdriver(naux, ntdg, iter_tag, type_tag, grid_tag, &
                          vxc_sol, vx_sol, ref_den, den, v_ref, &
                          ntg, grid, &
                          i_gridsize)
    end if

  end subroutine plot_driver
  ! This subroutine is supposed to order all the different cases and allow for
  ! plotting without the clutter of having it show up somwhere in the code
  subroutine plot_subdriver(naux, ntdg, iter_tag, type_tag, grid_tag, vxc_sol, &
                            vx_sol, ref_den, den, v_ref, ntg, grid, i_gridsize)
    use memory, only: memory_allocate, memory_release
    use ksinv_mod, only: l_plot_vx, l_plot_vc, l_plot_vref, l_plot_vxc, &
                         l_plot_rho_ks, l_plot_rho_ref, l_plot_rho_diff
    use iso_fortran_env, only: dp => real64
    implicit none

    ! Input
    integer, intent(in) :: naux, ntdg, ntg, i_gridsize
    real(dp), intent(in) ::  grid(3, i_gridsize)
    character(len=*) :: iter_tag, type_tag
    character(5) :: grid_tag
    real(dp), intent(in) :: v_ref(naux), vxc_sol(naux), vx_sol(naux)
    real(dp), intent(in) :: ref_den(ntdg), den(ntdg)

    ! Temporary
    real(dp) :: P_delta(ntdg)
    real(dp), pointer  ::  vc_sol(:)

    P_delta = ref_den-den
    vc_sol => memory_allocate(naux)
    vc_sol = vxc_sol

    call daxpy_x(naux, -1d0, vx_sol, 1, vc_sol, 1)
    call coulomb_basis('OEP', 0) ! u need to do this, otherwise you provoke a crash

    ! grid is now input, grid needs directionyl nametag !
    if (l_plot_vx) call plot_potentials(vx_sol, naux, i_gridsize, grid, 'vx', &
                                        iter_tag, grid_tag, type_tag)
    if (l_plot_vc) call plot_potentials(vc_sol, naux, i_gridsize, grid, 'vc', &
                                        iter_tag, grid_tag, type_tag)
    if (l_plot_vxc) call plot_potentials(vxc_sol, naux, i_gridsize, grid, &
                                         'vxc', iter_tag, grid_tag, type_tag)
    if (l_plot_vref) call plot_potentials(v_ref, naux, i_gridsize, grid, &
                                          'vref', iter_tag, grid_tag, type_tag)

    if (l_plot_rho_ks) call plot_densities(ntg, ntdg, den, i_gridsize, grid, &
                                           'rho_ks', iter_tag, grid_tag, &
                                           type_tag)
    if (l_plot_rho_ref) call plot_densities(ntg, ntdg, ref_den, i_gridsize, &
                                            grid, 'rho_ref', iter_tag, &
                                            grid_tag, type_tag)
    if (l_plot_rho_diff) call plot_densities(ntg, ntdg, P_delta, i_gridsize, &
                                             grid, 'rho_diff', iter_tag, &
                                             grid_tag, type_tag)

    call memory_release(vc_sol)

  end subroutine plot_subdriver

  ! Simple Rewrite of Plotting routines
  ! Original routines are buggy messy and unreadable, therefor this rewrite

  ! Takes a upper triangular version of a density matrix and writes it into a file onto a real space grid
  ! ao_matrix goes, in gridsized R^3 grid is made, ao_matrix is projected onto grid, prjection is written to file
  subroutine plot_densities(ntg, ntdg, ao_matrix, gridsize, grid, name_, &
                            iter_tag, grid_tag, type_tag)
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: ntg, ntdg, gridsize
    real(dp), intent(in) :: ao_matrix(ntdg)
    character(len=*), intent(in) :: name_, iter_tag, grid_tag, type_tag
    real(dp), intent(in)  :: grid(3, gridsize) ! real space grid for plotting
    real(dp)  :: rs_vector(gridsize)

    call get_rho(grid, gridsize, ao_matrix, ntdg, ntg, rs_vector)
    call print_to_file(grid, gridsize, rs_vector, name_, iter_tag, grid_tag, &
                       type_tag)

  end subroutine plot_densities

  subroutine get_rho(grid, gridsize, ao_matrix, ntdg, ntg, rs_vector)
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: gridsize, ntg, ntdg
    real(dp), intent(in) :: grid(3, gridsize)
    real(dp), intent(in) :: ao_matrix(ntdg)
    real(dp), intent(out) :: rs_vector(gridsize)

    integer :: map_(ntg), nmap, map_inverse(ntg)
    real(dp) :: orbval(gridsize, 1, ntg)
    real(dp)  :: tolorb
    integer :: imap

    integer  :: icentres(ntg)
    real(dp) :: zk(gridsize)

    icentres = 0
    call get_inpf('TOLORB', 'KSOPT', tolorb)
    if (tolorb .lt. 0d0) tolorb = 1d-11
    call grid_orbital_value(gridsize, 0, 1, grid, tolorb, map_, nmap, orbval, &
                            icentres)
    map_inverse(:) = -1
    do imap = 1, nmap
      map_inverse(map_(imap)) = imap
    end do
    call dft_rho_matrix_0(gridsize, orbval, 1, map_inverse, nmap, ao_matrix, &
                          rs_vector, zk)

  end subroutine get_rho

  subroutine plot_potentials(oep_vector, naux, gridsize, grid, name_, &
                             iter_tag, grid_tag, type_tag)
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: naux, gridsize
    real(dp), intent(in) :: oep_vector(naux)
    character(len=*), intent(in) :: name_, iter_tag, grid_tag, type_tag
    real(dp), intent(in)  :: grid(3, gridsize) ! real space grid for plotting
    real(dp)  :: rs_vector(gridsize)

    call auxbasisspace_to_realspace(oep_vector, naux, grid, gridsize, rs_vector)
    call print_to_file(grid, gridsize, rs_vector, name_, iter_tag, grid_tag, &
                       type_tag)

  end subroutine plot_potentials

  ! calculates from the oep vector, the representation of this vector on a real space grid
  ! i.e. v_x(aux-basis) -> v_x(r)
  ! needs grid and coefficient-vector for oep-basis as input
  subroutine auxbasisspace_to_realspace(oep_vector, naux, grid, gridsize, &
                                        rs_vector)
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: naux, gridsize
    real(dp), intent(in) :: oep_vector(naux)
    real(dp), intent(in)  :: grid(3, gridsize)
    real(dp), intent(out)  :: rs_vector(gridsize)
    integer :: ichunk, manifold, rest
    real(dp)  :: auxpot(naux, gridsize)

    ! Functions
    real(dp), external :: ddot_x

    ! this is the one component that is needed auxpot is only allocated, has no values that come in here
    ! call coulomb_auxbasis_pot(npt,pt,auxpot,naux) ! -> issue npt is some block, i dont care for blocks
    rs_vector = 0d0
    ! coulomb_auxbasis_pot has a max accepted gridsize of 1024, so larger grids must be subdivided into blocks
    ! This is why here batchwise auxpot is made
    if (gridsize .gt. 1024) then
      manifold = int(gridsize/1024)
      rest = mod(gridsize, 1024)
      do ichunk = 1, manifold
        call coulomb_auxbasis_pot(1024, &
                                  grid(:, (ichunk-1)*1024+1:ichunk*1024), &
                                  auxpot(:, (ichunk-1)*1024+1:ichunk*1024), &
                                  naux)
      end do
      call coulomb_auxbasis_pot(rest, &
                                grid(:, manifold*1024+1:manifold*1024+rest), &
                                auxpot(:, manifold*1024+1:manifold*1024+rest), &
                                naux)
    else
      call coulomb_auxbasis_pot(gridsize, grid, auxpot, naux)
    end if

    call dgemm_x('T', 'N', gridsize, 1, naux, -1d0, &
                 auxpot, naux, &
                 oep_vector, naux, &
                 0d0, rs_vector, gridsize)

  end subroutine auxbasisspace_to_realspace
  ! creates file in runnig directory and writes values of grid there
  subroutine print_to_file(grid, gridsize, rs_vector, name_, iter_tag, &
                           grid_tag, type_tag)
    use molpro_options, only: molpro_pwd
    use iso_fortran_env, only: dp => real64
    implicit none
    !
    ! Input
    ! name supposed to contain type like vx or vc may contain anything
    ! tag  contains current iteration or "final" if final iteration
    ! rs_vector contains values for "result" axis in grid
    ! grid contains R^3-grid
    ! gridsize number of gridpoints
    !
    character(len=*), intent(in) :: name_, iter_tag, grid_tag, type_tag
    integer, intent(in) :: gridsize
    real(dp), intent(in)  :: rs_vector(gridsize), grid(3, gridsize)

    integer :: iunit, ipoint

    open (newunit=iunit, &
          file=molpro_pwd//trim(name_)//'-'//trim(adjustl(type_tag)) &
          //'-'//trim(iter_tag)//'.'//trim(grid_tag), &
          status='replace', form='formatted')
    do ipoint = 1, gridsize
      write (iunit, '(4(1x,ES21.14))') grid(:, ipoint), rs_vector(ipoint)
    end do
    close (iunit)

  end subroutine print_to_file

  subroutine generate_line_grid(r1, r2, gridsize, dist, grid)
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in)  :: gridsize
    real(dp), intent(inout) :: grid(3, gridsize)
    real(dp), intent(in)  :: r1(3), r2(3)
    real(dp), intent(in)  :: dist

    ! internal
    integer :: ipt, i, npt_half
    real(dp) :: r12(3), d12, r12h(3), d1(3), d2(3), step, d0(3)

    !...generate grid between points r1 and r2
    npt_half = gridsize/2
    r12(:) = r1(:)-r2(:)
    d12 = sqrt(dot_product(r12, r12))
    r12h(:) = r1(:)-r12(:)/2d0
    d1(:) = r1(:)-r12h(:)
    d2(:) = r2(:)-r12h(:)
    d1(:) = d1(:)/sqrt(dot_product(d1, d1))
    d2(:) = d2(:)/sqrt(dot_product(d2, d2))
    step = dist/dble(npt_half)
    ipt = 1
    do i = npt_half, 1, -1
      grid(:, ipt) = r12h(:)+dble(i)*d1(:)*dist/real(npt_half, dp)
      d0(:) = grid(:, ipt)-r12h(:)
      ipt = ipt+1
    end do
    do i = 1, npt_half
      grid(:, ipt) = r12h(:)+dble(i)*d2(:)*dist/real(npt_half, dp)
      d0(:) = grid(:, ipt)-r12h(:)
      ipt = ipt+1
    end do

  end subroutine generate_line_grid

end module ksinv_plot

