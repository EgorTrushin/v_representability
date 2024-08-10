!
! Module with input parameters and its parsing and processing as well as
! with subroutine for initialization in the beginning of calculation
!
module ksinv_mod
  use acfd_derived_datatypes, only: matrix_type
  use iso_fortran_env, only: dp => real64
  implicit none

  ! Group of input parameters, i.e., parameters parsed from input file

  integer :: irec_den, ifil_den
  integer :: irec_orb, ifil_orb
  integer :: irec_orb_out, ifil_orb_out
  real(dp) :: thr_scf_inversion ! convergence threshold for the potential inversion
  integer :: maxit
  integer :: minit
  real(dp) :: mixing_a
  integer :: mix_switch_iter
  real(dp) :: reference_energy ! total energy of calculation which delivered density matrix to be inverted
  integer :: iverbose
  logical :: l_backfilter

  real(dp) :: thr_symmetry ! threshold for symmetrizing the OEP basis
  real(dp) :: thr_overlap_oep  ! threshold for processing auxbas by throwing away small EV in overlap matrix
  real(dp) :: thr_fai_oep  ! threshold for processing auxbas by throwing away small EV in D_{III}-matrix
  logical :: l_vref_fa
  logical :: l_vh_oep

  real(dp) :: thr_oep  ! threshold for eigenvalue difference in lambda calculation for OEP
  real(dp) :: thr_solve ! threshold for diagonalization in eigenproblem solver for OEP equation
  character(32) :: csolvemeth ! eigenproblem solver for OEP equation

  logical :: l_plot_always ! plot in the end (F) or at every iteration (T)

  logical  :: l_plot_vx, l_plot_vc, l_plot_vref, l_plot_vxc
  logical  :: l_plot_rho_ks, l_plot_rho_ref, l_plot_rho_diff
  logical  :: l_plot_x_line, l_plot_y_line, l_plot_z_line
  integer  :: i_gridsize
  real(dp) :: r0_plot_range

  ! Arrays and variables which are not input parameters

  integer :: no, nv ! number of occupied and virtual orbitals
  integer :: no_a, no_b ! number of occupied orbitals in the alpha and beta channel respectively
  real(dp), pointer :: orb(:), eig(:), den(:), occ(:) ! orbitals, eigenvalues, denisty, occupations
  real(dp), pointer :: orb_a(:), eig_a(:), den_a(:), &
                       orb_b(:), eig_b(:), den_b(:) ! orbitals, eigenvalues, denisty, occupations for alpha and beta spin channel
  real(dp), pointer :: smh(:) ! S^{-1/2}
  real(dp), pointer :: smat(:, :), den_sqr(:, :), UT_smat(:) ! S, P

  real(dp), pointer :: vx(:) ! exchange potential
  real(dp), pointer :: vx_a(:), vx_b(:) ! UKSINV: exchange potential of alpha and beta charge density
  real(dp), pointer :: vxc(:) ! exchange-correlation potential
  real(dp), pointer :: vxc_a(:), vxc_b(:) ! UKSINV: exchange-correlation potential of alpha and beta charge density
  real(dp), pointer :: vx_sol(:)
  real(dp), pointer :: vxc_sol(:) ! representation of exchange correlation potential in auxiliary basis
  real(dp), pointer :: vxc_sol_a(:), vxc_sol_b(:) ! representation of exchange correlation potential in auxiliary basis set
  real(dp), pointer :: vx_sol_a(:), vx_sol_b(:)   ! representation of exchange potential in auxiliary basis set
  real(dp), pointer :: vj(:) ! Coulomb potential
  real(dp), pointer :: h0(:) ! one electron potential

  integer :: naux_constraint ! number of basis functions after reduction and constraints
  integer :: naux_charge ! number of basis functions after reduction and charge constraint
  type(matrix_type) :: trans_mat_constraint ! transformation matrix for reduction and constraints
  type(matrix_type) :: trans_mat_constraint_inv
  integer :: naux_oep ! number of OEP basis functions

  real(dp), pointer :: ref_den(:) ! the reference density matrix to which the inversion will converge
  real(dp), pointer :: ref_den_a(:), ref_den_b(:) ! UKSINV: majority (alpha, a) and minority (beta, b) channel reference charge densities
  real(dp), pointer :: ref_dens(:) ! UKSINV: the reference spin-density to which the unrestricted inversion will converge
  real(dp), pointer :: vj_frozen(:) ! made to get density from fci density, then freeze it in the calculation

  real(dp) :: c_vsol_delta ! sets up logistic funtion for mixing, in current implementation, inverse of initial |t|
  real(dp) :: c_vsol_delta_a, c_vsol_delta_b  ! sets up logistic funtion for mixing, in current implementation, inverse of initial |t| of alpha and beta channel
  real(dp) :: inversion_rhs ! contains |t^{\Delta}| and can be converged
  real(dp) :: inversion_rhs_a, inversion_rhs_b ! contains |t^{\Delta}| of alpha and beta channel
  real(dp), pointer :: int_fkl(:, :, :)
  real(dp), pointer :: u_frozen(:) ! the u-vector needed for the generation of the reference potential
  real(dp), pointer :: ref_v_ref_den(:) ! refernce potential made from reference density matrix
  real(dp), pointer :: v_kin(:) ! kinetic energy operator
  logical :: l_stopped = .false.

  logical :: l_swap, l_swap2
  integer :: iswap1, iswap2
  logical :: l_frac_occ
  character(32) :: system

  real(dp), pointer :: orbb_frac(:), denb_frac(:), vxb1(:), vxb2(:), occb1(:), occb2(:)

contains

  ! read input parameters
  subroutine get_input_parameters()
    use common_cinput, only: ncol
    use common_tapes, only: iout
    use acfd_parameters, only: iverbose_ => iverbose, &
                               thr_oep_ => thr_oep, &
                               thr_solve_ => thr_solve, &
                               csolvemeth_ => csolvemeth
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, parameter :: MAXCOL = 20

    character(32) :: str, namx
    integer :: icol
    real(dp) :: ain(MAXCOL+1)

    call set_default()  ! set default values

    ! read input values
    do icol = 2, ncol
      call geta(icol, namx, str, ain, 1)
      if (namx .eq. 'REFDEN') then
        irec_den = int(ain(1))
        ifil_den = int((ain(1)-dble(irec_den))*10.1d0)
      elseif (namx .eq. 'ORB') then
        irec_orb = int(ain(1))
        ifil_orb = int((ain(1)-dble(irec_orb))*10.1d0)
      elseif (namx .eq. 'SAVE') then
        irec_orb_out = int(ain(1))
        ifil_orb_out = int((ain(1)-dble(irec_orb_out))*10.1d0)
      elseif (namx .eq. 'E_REF') then
        reference_energy = ain(1)
      elseif (namx .eq. 'THR_INVERSION') then
        thr_scf_inversion = ain(1)
      elseif (namx .eq. 'MAXIT') then
        maxit = int(ain(1))
      elseif (namx .eq. 'MINIT') then
        minit = int(ain(1))
      elseif (namx .eq. 'MIXING_A') then
        mixing_a = ain(1)
      elseif (namx .eq. 'MIX_SWITCH_ITER') then
        mix_switch_iter = int(ain(1))
      elseif (namx .eq. 'BACKFILTER') then
        l_backfilter = ain(1) .ne. 0
      elseif (namx .eq. 'THR_SYM') then
        thr_symmetry = ain(1)
      elseif (namx .eq. 'THR_FAI_OEP') then
        thr_fai_oep = ain(1)
      elseif (namx .eq. 'THR_OVERLAP_OEP') then
        thr_overlap_oep = ain(1)
      elseif (namx .eq. 'VREF_FA') then
        l_vref_fa = ain(1) .ne. 0
      elseif (namx .eq. 'VHOEP') then
        l_vh_oep = ain(1) .ne. 0
      elseif (namx .eq. 'THR_OEP') then
        thr_oep = ain(1)
      elseif (namx .eq. 'THR_SOLVE') then
        thr_solve = ain(1)
      elseif (namx .eq. 'SOLVE') then
        csolvemeth = str
      elseif (namx .eq. 'PLOT_ALWAYS') then
        l_plot_always = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_VC') then
        l_plot_vc = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_VX') then
        l_plot_vx = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_VXC') then
        l_plot_vxc = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_VREF') then
        l_plot_vref = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_RHO_KS') then
        l_plot_rho_ks = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_RHO_REF') then
        l_plot_rho_ref = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_RHO_DIFF') then
        l_plot_rho_diff = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_X') then
        l_plot_x_line = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_Y') then
        l_plot_y_line = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_Z') then
        l_plot_z_line = ain(1) .ne. 0
      elseif (namx .eq. 'GRIDSIZE') then
        i_gridsize = int(ain(1))
      elseif (namx .eq. 'PLOTRANGE') then
        r0_plot_range = ain(1)
      elseif (namx .eq. 'VERB') then
        iverbose = int(ain(1))
      elseif (namx .eq. 'NOA') then
        no_a = int(ain(1))
      elseif (namx .eq. 'NOB') then
        no_b = int(ain(1))
      elseif (namx .eq. 'SWAP') then
        l_swap = ain(1) .ne. 0
      elseif (namx .eq. 'ISWAP1') then
        iswap1 = int(ain(1))
      elseif (namx .eq. 'ISWAP2') then
        iswap2 = int(ain(1))
      elseif (namx .eq. 'FRACOCC') then
        l_frac_occ = ain(1) .ne. 0
      elseif (namx .eq. 'SYSTEM') then
        system = str
      elseif (namx .eq. 'SWAP2') then
        l_swap2 = ain(1) .ne. 0
      end if
    end do

    ! provide parameters required for EXX OEP to acfd_parameters
    iverbose_ = iverbose
    thr_oep_ = thr_oep
    thr_solve_ = thr_solve
    csolvemeth_ = csolvemeth

    ! print used parameters
    call print_parameters()

  contains

    ! sets default values for parameters
    subroutine set_default()
      implicit none

      irec_den = 0
      ifil_den = 0

      irec_orb = 2100
      ifil_orb = 2
      irec_orb_out = 2101
      ifil_orb_out = 2

      maxit = 30
      minit = 3

      thr_scf_inversion = 1d-10
      mixing_a = 0.1d0
      mix_switch_iter = 3
      l_backfilter = .true.

      thr_overlap_oep = 1d-99
      thr_fai_oep = 1d-14
      thr_symmetry = 0d0
      l_vref_fa = .true.
      l_vh_oep = .false.

      thr_oep = 1d-6
      thr_solve = 1d-99
      csolvemeth = 'GTSVD'

      l_plot_always = .false.
      l_plot_vc = .false.
      l_plot_vx = .false.
      l_plot_vxc = .false.
      l_plot_vref = .false.
      l_plot_rho_ks = .false.
      l_plot_rho_ref = .false.
      l_plot_rho_diff = .false.
      l_plot_x_line = .false.
      l_plot_y_line = .false.
      l_plot_z_line = .false.
      i_gridsize = 2048
      r0_plot_range = 20d0

      iverbose = 0

      no_a = -1
      no_b = -1

      l_swap = .false.
      iswap1 = 3
      iswap2 = 5
      l_frac_occ = .false.
      system = 'O'
      l_swap2 = .false.

    end subroutine set_default

    ! prints parameters
    subroutine print_parameters()
      implicit none

      write (iout, '(1x,a,i4,a,i1)') &
        'Read Reference Density from ', irec_den, '.', ifil_den
      write (iout, '(1x,a,i4,a,i1)') &
        'Read Orbital Records from ', irec_orb, '.', ifil_orb
      write (iout, '(1x,a,i4,a,i1)') &
        'Write Orbital Records to ', irec_orb_out, '.', ifil_orb_out
      write (iout, '(1x,a,1x,f20.12)') 'Reference energy: ', reference_energy
      write (iout, '(1x,a,1x,es12.3)') &
        'Threshold for inversion convergence: ', thr_scf_inversion
      write (iout, '(1x,a,i4)') 'Maximum number of iterations:', maxit
      write (iout, '(1x,a,i4)') 'Minimum number of iterations:', minit
      write (iout, '(1x,a,1x,es12.3)') 'mixing_a:', mixing_a
      write (iout, '(1x,a,1x,i3)') 'mix_switch_iter:', mix_switch_iter
      write (iout, '(1x,a,l2)') 'Use backfiltering:', l_backfilter
      write (iout, '(1x,a)') 'Thresholds for processing of OEP basis set:'
      write (iout, '(1x,a,es12.3)') 'thr_symmetry:', thr_symmetry
      write (iout, '(1x,a,1x,es12.3)') 'thr_overlap_oep: ', thr_overlap_oep
      write (iout, '(1x,a,1x,es12.3)') 'thr_fai_oep: ', thr_fai_oep
      write (iout, '(1x,a,l2)') &
        'Fermi-Amaldi as reference potential:', l_vref_fa
      write (iout, '(1x,a,l2)') &
        'Construct the Hartree potential in the OEP basis:', l_vh_oep
      write (iout, '(1x,a,1x,es12.3)') &
        'Threshold for eigenvalue difference in calculation of Lambda:', thr_oep
      write (iout, '(1x,a,1x,es12.3)') &
        'Threshold for eigenproblem solver of OEP equation:', thr_solve
      write (iout, '(1x,a)') 'Solver for eigenproblems: '//trim(csolvemeth)

      if ((l_plot_x_line .or. l_plot_y_line .or. l_plot_z_line) .and. &
          (l_plot_vc .or. l_plot_vx .or. l_plot_vxc .or. l_plot_vref .or. &
           l_plot_rho_ks .or. l_plot_rho_ref .or. l_plot_rho_diff)) then
        write (iout, '(1x,a)') 'Plotting parameters:'
        write (iout, '(1x,a,l2)') 'l_plot_always:', l_plot_always
        write (iout, '(1x,a,l2)') 'l_plot_vc:', l_plot_vc
        write (iout, '(1x,a,l2)') 'l_plot_vx:', l_plot_vx
        write (iout, '(1x,a,l2)') 'l_plot_vxc:', l_plot_vxc
        write (iout, '(1x,a,l2)') 'l_plot_vref:', l_plot_vref
        write (iout, '(1x,a,l2)') 'l_plot_rho_ks:', l_plot_rho_ks
        write (iout, '(1x,a,l2)') 'l_plot_rho_ref:', l_plot_rho_ref
        write (iout, '(1x,a,l2)') 'l_plot_rho_diff:', l_plot_rho_diff
        write (iout, '(1x,a,l2)') 'l_plot_x_line:', l_plot_x_line
        write (iout, '(1x,a,l2)') 'l_plot_y_line:', l_plot_y_line
        write (iout, '(1x,a,l2)') 'l_plot_z_line:', l_plot_z_line
        write (iout, '(1x,a,i6)') 'i_gridsize:', i_gridsize
        write (iout, '(1x,a,f12.4)') 'r0_plot_range:', r0_plot_range
      end if

      if ((no_a > -1) .and. (no_b > -1)) then
        write (iout, '(1x,a,i4)') 'noa:', no_a
        write (iout, '(1x,a,i4)') 'nob:', no_b
      end if

      write (iout, '(1x,a,l2)') 'l_swap:', l_swap
      write (iout, '(1x,a,l2)') 'l_swap2:', l_swap2
      if ((l_swap) .or. (l_swap2)) write (iout, '(1x,a,1x,i3,i3)') "orbital to swap:", iswap1, iswap2

      write (iout, '(1x,a,l2)') 'l_frac_occ:', l_frac_occ
      if (l_frac_occ) write (iout, '(1x,a,1x,a)') "system:", system

      write (iout, *)

    end subroutine print_parameters
  end subroutine get_input_parameters

  ! allocates variables on stack and gets values needed for calculation from records
  subroutine ksinv_init()
    use common_cbas, only: ntg, ntdg, ntqg
    use common_code, only: ekin, ho1, rs
    use memory, only: memory_allocate
    use iso_fortran_env, only: dp => real64
    implicit none

    ! auxiliary variables
    real(dp) :: r_number_of_electrons
    character(16) :: c_dummy
    integer :: i_dummy1, i_dummy2

    real(dp), external :: dsum

    call basis_number('OEP', naux_oep)
    ! initialize variables for a spin-restricted potential inversion calculation
    orb => memory_allocate(ntqg)
    orb = 0d0
    eig => memory_allocate(ntg)
    eig = 0d0
    den => memory_allocate(ntdg)
    den = 0d0
    ref_den => memory_allocate(ntdg)
    ref_den = 0d0
    occ => memory_allocate(ntg)
    occ = 0d0
    vx => memory_allocate(ntdg)
    vx = 0d0
    vxc => memory_allocate(ntdg)
    vxc = 0d0
    vxc_sol => memory_allocate(naux_oep)
    vxc_sol = 0d0
    vx_sol => memory_allocate(naux_oep)
    vx_sol = 0d0
    int_fkl => memory_allocate(naux_oep, ntg, ntg)
    int_fkl = 0d0
    vj_frozen => memory_allocate(ntdg)
    vj_frozen = 0d0
    u_frozen => memory_allocate(naux_oep)
    u_frozen = 0d0
    ref_v_ref_den => memory_allocate(ntdg)
    ref_v_ref_den = 0d0

    ! get values from the record
    call get_den(ref_den, irec_den, ifil_den, 1, 1)
    call get_occ(occ, irec_orb, ifil_orb, 1, 1)
    call getvar('NELEC', r_number_of_electrons, &
                c_dummy, i_dummy1, i_dummy2, 1, 1)

    v_kin => memory_allocate(ntdg)
    v_kin = 0d0
    call les(v_kin, ntdg, 1, ekin)

    no = int(r_number_of_electrons/2d0)

    nv = ntg-no

    vj => memory_allocate(ntdg)
    vj = 0d0
    h0 => memory_allocate(ntdg)
    h0 = 0d0

    call les(h0, ntdg, 1, ho1)

    smh => memory_allocate(ntqg)
    smh = 0d0

    call les(smh, ntdg, 1, rs)
    call expan(smh, smh, 1, 1, 1)

  end subroutine ksinv_init
  ! allocates variables on stack and gets values needed for calculation from records
  subroutine uksinv_init()
    use common_cbas, only: ntg, ntdg, ntqg
    use common_code, only: ekin, ho1, rs
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    !use ksinv_utils, only: UT_smat
    use iso_fortran_env, only: dp => real64
    implicit none

    ! auxiliary variables
    real(dp) :: r_number_of_electrons
    character(16) :: c_dummy
    integer :: i_dummy1, i_dummy2

    real(dp), external :: dsum, spur

    call basis_number('OEP', naux_oep)
    ! initialize variables for a spin-restricted potential inversion calculation
    orb => memory_allocate(ntqg)
    orb = 0d0
    orb_a => memory_allocate(ntqg)
    orb_a = 0d0
    orb_b => memory_allocate(ntqg)
    orb_b = 0d0
    eig => memory_allocate(ntg)
    eig = 0d0
    eig_a => memory_allocate(ntg)
    eig_a = 0d0
    eig_b => memory_allocate(ntg)
    eig_b = 0d0
    den => memory_allocate(ntdg)
    den = 0d0
    den_a => memory_allocate(ntdg)
    den_a = 0d0
    den_b => memory_allocate(ntdg)
    den_b = 0d0
    ref_den => memory_allocate(ntdg)
    ref_den = 0d0
    ref_dens => memory_allocate(ntdg)
    ref_dens = 0d0
    ref_den_a => memory_allocate(ntdg)
    ref_den_a = 0d0
    ref_den_b => memory_allocate(ntdg)
    ref_den_b = 0d0
    occ => memory_allocate(ntg)
    occ = 0d0
    vx => memory_allocate(ntdg)
    vx = 0d0
    vx_a => memory_allocate(ntdg)
    vx_a = 0d0
    vx_b => memory_allocate(ntdg)
    vx_b = 0d0
    vxc => memory_allocate(ntdg)
    vxc = 0d0
    vxc_a => memory_allocate(ntdg)
    vxc_a = 0d0
    vxc_b => memory_allocate(ntdg)
    vxc_b = 0d0
    vxc_sol => memory_allocate(naux_oep)
    vxc_sol = 0d0
    vxc_sol_a => memory_allocate(naux_oep)
    vxc_sol_a = 0d0
    vxc_sol_b => memory_allocate(naux_oep)
    vxc_sol_b = 0d0
    vx_sol => memory_allocate(naux_oep)
    vx_sol = 0d0
    vx_sol_a => memory_allocate(naux_oep)
    vx_sol_a = 0d0
    vx_sol_b => memory_allocate(naux_oep)
    vx_sol_b = 0d0
    int_fkl => memory_allocate(naux_oep, ntg, ntg)
    int_fkl = 0d0
    vj_frozen => memory_allocate(ntdg)
    vj_frozen = 0d0
    u_frozen => memory_allocate(naux_oep)
    u_frozen = 0d0
    ref_v_ref_den => memory_allocate(ntdg)
    ref_v_ref_den = 0d0

    ! get values from the record
    call get_den(ref_den, irec_den, ifil_den, 1, 1)
    call get_den(ref_dens, irec_den, ifil_den, 2, 1)
    call get_occ(occ, irec_orb, ifil_orb, 1, 1)
    call getvar('NELEC', r_number_of_electrons, &
                c_dummy, i_dummy1, i_dummy2, 1, 1)

    v_kin => memory_allocate(ntdg)
    v_kin = 0d0
    call les(v_kin, ntdg, 1, ekin)

    no = int(r_number_of_electrons/2d0)

    nv = ntg-no

    vj => memory_allocate(ntdg)
    vj = 0d0
    h0 => memory_allocate(ntdg)
    h0 = 0d0

    call les(h0, ntdg, 1, ho1)

    smh => memory_allocate(ntqg)
    smh = 0d0

    ref_den_a = 0.5d0*(ref_den+ref_dens)
    ref_den_b = 0.5d0*(ref_den-ref_dens)

    smat => memory_allocate(ntg, ntg)
    UT_smat => memory_allocate(ntdg)
    call basis_sqr("ORBITAL", 'S', smat, 0d0)
    call reduc(smat, UT_smat, 1, 1, 1d0, 0)
    no_a = NINT(spur(UT_smat, ref_den_a))
    no_b = NINT(spur(UT_smat, ref_den_b))
    write (iout, '(1x,a,i4)') 'noa:', no_a
    write (iout, '(1x,a,i4)') 'nob:', no_b

    call memory_release(smat)

    call les(smh, ntdg, 1, rs)
    call expan(smh, smh, 1, 1, 1)

  end subroutine uksinv_init

  subroutine deallocate_arrays()
    use memory, only: memory_release
    implicit none
    if (allocated(trans_mat_constraint%mat)) &
      deallocate (trans_mat_constraint%mat)
    if (allocated(trans_mat_constraint_inv%mat)) &
      deallocate (trans_mat_constraint_inv%mat)
  end subroutine deallocate_arrays
end module ksinv_mod
