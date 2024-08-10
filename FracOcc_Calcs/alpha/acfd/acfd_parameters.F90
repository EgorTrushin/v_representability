!
! This module contains the parameters for calculations
!
module acfd_parameters
  use iso_fortran_env, only: dp => real64
  implicit none

  logical :: ldfit  ! if .true., density fitting is used
  integer :: iverbose  ! govern level of verbosity in output
  real(dp) :: vc_scaling  ! scaling factor for Coulomb potential

  logical :: l_sigma  ! govern usage of sigma-functional
  character(32) :: sigma_param  ! determine parametrization for sigma-functional
  logical :: l_write_sigma  ! govern writing dat-file with sigma-values

  ! record number of the orbital record used as input
  integer :: irec_orb
  integer :: ifil_orb

  integer :: nquad  ! number of quadrature points per interval for frequency integration
  integer :: nquadint  ! number of logarithmic intervals for frequency integration
  real(dp) :: quad_w0  ! scaling factor for frequency integration

  real(dp) :: thr_overlap_ri  ! threshold for processing auxbas by throwing away small EV in overlap matrix
  real(dp) :: thr_fai_ri  ! threshold for processing auxbas by throwing away small EV in D_{III}-matrix

  real(dp) :: thr_rpa  ! threshold for eigenvalue difference in X0 construction (RPA)
  real(dp) :: thr_oep  ! threshold for eigenvalue difference in lambda calculation for OEP
  real(dp) :: thr_solve ! threshold for diagonalization in eigenproblem solver for OEP equation

  integer :: irec_orb_out ! record number of the orbital record used as output
  integer :: ifil_orb_out ! file number for the orbital record used as output

  integer :: maxit ! maximum number of interations in SCF calculations
  integer :: minit ! minimum amount of iterations in SCF calculations
  integer :: maxdim_diis ! maximum size of DIIS matrix
  logical :: l_fixmix ! mix old fock matrix into new fock matrix at a fixed rate instead of DIIS
  real(dp) :: fixmix_rate ! amount of old fock matrix in new fock matrix for l_fixmix=.true.

  real(dp) :: thr_symmetry ! threshold for symmetrizing the OEP basis
  real(dp) :: thr_overlap_oep  ! threshold for processing auxbas by throwing away small EV in overlap matrix
  real(dp) :: thr_fai_oep  ! threshold for processing auxbas by throwing away small EV in D_{III}-matrix

  logical :: l_vref_fa ! if true, Fermi-Amaldi reference potential is constructed
  logical :: l_homo ! if true, HOMO condition is used

  real(dp) :: thr_scf_energy ! convergence threshold for the energy in SCF calculations
  real(dp) :: thr_scf_density ! convergence threshold for the density in SCF calculations

  logical :: l_spin_sym ! if .true., apply spin-symmetrization
  logical :: l_vref_fa_sameab
  logical :: l_vh_oep ! if .true., Hartree potential is constructed in OEP basis

  logical :: l_test_pot ! if .true., do test for potentials after SCF

  logical :: l_plot_x_line, l_plot_y_line, l_plot_z_line
  logical :: l_plot_always
  character(32) :: csolvemeth ! eigenproblem solver for OEP equation

  logical :: l_swap, l_swap2
  integer :: iswap1, iswap2
  logical :: l_frac_occ
  character(32) :: system

  public :: ldfit, iverbose, vc_scaling, &
            l_sigma, sigma_param, l_write_sigma, &
            irec_orb, ifil_orb, nquad, nquadint, quad_w0, &
            thr_overlap_ri, thr_fai_ri, thr_rpa, thr_oep, thr_solve, &
            get_parameters_rirpa, irec_orb_out, ifil_orb_out, &
            maxit, minit, maxdim_diis, l_fixmix, fixmix_rate, &
            thr_symmetry, thr_overlap_oep, thr_fai_oep, &
            thr_scf_energy, thr_scf_density, l_vref_fa, l_homo, &
            l_spin_sym, l_vref_fa_sameab, l_vh_oep, &
            csolvemeth, l_test_pot, l_plot_always, &
            l_plot_x_line, l_plot_y_line, l_plot_z_line, &
            get_parameters_scexx, &
            l_swap, l_swap2, iswap1, iswap2, l_frac_occ, system

  private

contains

  ! Read parameters for the RIRPA and URIRPA calculations
  subroutine get_parameters_rirpa(method)
    use common_cinput, only: ncol
    use common_tapes, only: iout
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: MAXCOL = 20

    character(*), intent(in) :: method

    character(32) :: str, namx
    integer :: icol
    real(dp) :: ain(MAXCOL+1)

    call set_default(method)  ! set default values

    ! read input values
    do icol = 2, ncol
      call geta(icol, namx, str, ain, 1)
      if (namx .eq. 'ORB') then
        irec_orb = int(ain(1))
        ifil_orb = int((ain(1)-dble(irec_orb))*10.1d0)
      elseif (namx .eq. 'DFIT') then
        ldfit = ain(1) .ne. 0
      elseif (namx .eq. 'SIGMA') then
        l_sigma = ain(1) .ne. 0
      elseif (namx .eq. 'SIGMA_PARAM') then
        sigma_param = str
      elseif (namx .eq. 'WRITE_SIGMA') then
        l_write_sigma = ain(1) .ne. 0
      elseif (namx .eq. 'THR_FAI_RI') then
        thr_fai_ri = ain(1)
      elseif (namx .eq. 'THR_OVERLAP_RI') then
        thr_overlap_ri = ain(1)
      elseif (namx .eq. 'THR_RPA') then
        thr_rpa = ain(1)
      elseif (namx .eq. 'NQUAD') then
        nquad = int(ain(1))
      elseif (namx .eq. 'NQUADINT') then
        nquadint = int(ain(1))
        nquad = nquad*nquadint
      elseif (namx .eq. 'W0') then
        quad_w0 = ain(1)
      elseif (namx .eq. 'VC_SCAL') then
        vc_scaling = ain(1)
      elseif (namx .eq. 'VERB') then
        iverbose = int(ain(1))
      end if
    end do

    ! print used parameters
    call print_parameters(method)

  contains

    ! Sets default values for parameters depending on the method
    subroutine set_default(method)
      implicit none
      character(*), intent(in) :: method

      logical :: l_unrestricted

      l_unrestricted = (Trim(method) == 'URIRPA')

      ! Parameters that are the same in all methods
      iverbose = 0
      vc_scaling = 1d0
      ldfit = .true.

      thr_overlap_ri = 1d-99
      thr_fai_ri = 1d-14

      thr_rpa = 1d-6

      nquad = 50
      quad_w0 = 2.5d0
      nquadint = 1

      l_sigma = .true.
      sigma_param = "PBE_S2"
      l_write_sigma = .false.

      ! Parameters that are not the same in restricted and unrestricted methods
      if (l_unrestricted) then
        irec_orb = 2200
        ifil_orb = 2
      else
        irec_orb = 2100
        ifil_orb = 2
      end if

    end subroutine set_default

    ! Prints parameters depending on the method
    subroutine print_parameters(method)
      implicit none

      character(*), intent(in) :: method
      logical :: l_unrestricted

      l_unrestricted = (Trim(method) == 'URIRPA')

      write (iout, '(1x,a,i4,a,i1)') &
        'Read Orbital Records from ', irec_orb, '.', ifil_orb
      write (iout, '(1x,a,1x,es12.3)') &
        'Threshold for eigenvalue difference in construction of X0: ', thr_rpa
      write (iout, '(1x,a,1x,es12.3)') &
        'Threshold for processing RI basis thr_overlap_ri: ', thr_overlap_ri
      write (iout, '(1x,a,1x,es12.3)') &
        'Threshold for processing RI basis thr_fai_ri: ', thr_fai_ri
      write (iout, '(1x,a,l2)') &
        'Use density fitting for evaluation of reference xc energy:', ldfit
      write (iout, '(1x,a)') &
        'Gauss-Legendre integration for frequency with rational function:'
      write (iout, '(1x,a,i4)') &
        'Number of logarithmic integration intervals: ', nquadint
      write (iout, '(1x,a,i4)') &
        'Number of quadrature points per interval:', nquad
      write (iout, '(1x,a,es10.3)') 'Scaling factor: ', quad_w0
      if (abs(vc_scaling-1d0) > 1d-10) write (iout, '(/,1x,a,es10.3)') &
        'Scaling factor for correlation potential:', vc_scaling
      if (l_sigma) then
        write (iout, '(1x,a)') "Sigma-functional will be used"
        write (iout, '(1x,a,1x,a)') &
          "Sigma-functional parametrization:", sigma_param
        if (l_write_sigma) &
          write (iout, '(1x,a)') 'Sigma-values will be stored into sigma.dat'
      end if

      write (iout, '(a)') ''

    end subroutine print_parameters
  end subroutine get_parameters_rirpa

  ! Read parameters for the SCEXX and USCEXX calculations
  subroutine get_parameters_scexx(method)
    use common_cinput, only: ncol
    use common_tapes, only: iout
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: MAXCOL = 20

    character(*), intent(in) :: method

    character(32) :: str, namx
    integer :: icol
    real(dp) :: ain(MAXCOL+1)

    call set_default(method)  ! set default values

    ! read input values
    do icol = 2, ncol
      call geta(icol, namx, str, ain, 1)
      if (namx .eq. 'ORB') then
        irec_orb = int(ain(1))
        ifil_orb = int((ain(1)-dble(irec_orb))*10.1d0)
        irec_orb_out = irec_orb+1
        ifil_orb_out = ifil_orb
      elseif (namx .eq. 'SAVE') then
        irec_orb_out = int(ain(1))
        ifil_orb_out = int((ain(1)-dble(irec_orb_out))*10.1d0)
      elseif (namx .eq. 'DFIT') then
        ldfit = ain(1) .ne. 0
      elseif (namx .eq. 'MAXIT') then
        maxit = int(ain(1))
      elseif (namx .eq. 'MINIT') then
        minit = int(ain(1))
      elseif (namx .eq. 'MAXDIIS') then
        maxdim_diis = int(ain(1))
        maxdim_diis = max(4, maxdim_diis)
      elseif (namx .eq. 'FIXMIX') then
        l_fixmix = ain(1) .ne. 0
      elseif (namx .eq. 'MIXRATE') then
        fixmix_rate = ain(1)
      elseif (namx .eq. 'ENERGY') then
        thr_scf_energy = ain(1)
      elseif (namx .eq. 'DENSITY') then
        thr_scf_density = ain(1)
      elseif (namx .eq. 'THR_SYM') then
        thr_symmetry = ain(1)
      elseif (namx .eq. 'THR_FAI_OEP') then
        thr_fai_oep = ain(1)
      elseif (namx .eq. 'THR_OVERLAP_OEP') then
        thr_overlap_oep = ain(1)
      elseif (namx .eq. 'THR_OEP') then
        thr_oep = ain(1)
      elseif (namx .eq. 'THR_SOLVE') then
        thr_solve = ain(1)
      elseif (namx .eq. 'SOLVE') then
        csolvemeth = str
      elseif (namx .eq. 'SPIN_SYM') then
        l_spin_sym = ain(1) .ne. 0
      elseif (namx .eq. 'VHOEP') then
        l_vh_oep = ain(1) .ne. 0
      elseif (namx .eq. 'VREF_FA') then
        l_vref_fa = ain(1) .ne. 0
      elseif (namx .eq. 'HOMO') then
        l_homo = ain(1) .ne. 0
      elseif (namx .eq. 'VREF_FA_SAMEAB') then
        l_vref_fa_sameab = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_X') then
        l_plot_x_line = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_Y') then
        l_plot_y_line = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_Z') then
        l_plot_z_line = ain(1) .ne. 0
      elseif (namx .eq. 'PLOT_ALWAYS') then
        l_plot_always = ain(1) .ne. 0
      elseif (namx .eq. 'TEST_POT') then
        l_test_pot = ain(1) .ne. 0
      elseif (namx .eq. 'VERB') then
        iverbose = int(ain(1))
      elseif (namx .eq. 'SWAP') then
        l_swap = ain(1) .ne. 0
      elseif (namx .eq. 'SWAP2') then
        l_swap2 = ain(1) .ne. 0
      elseif (namx .eq. 'ISWAP1') then
        iswap1 = int(ain(1))
      elseif (namx .eq. 'ISWAP2') then
        iswap2 = int(ain(1))
      elseif (namx .eq. 'FRACOCC') then
        l_frac_occ = ain(1) .ne. 0
      elseif (namx .eq. 'SYSTEM') then
          system = str
      end if
    end do

    ! print used parameters
    call print_parameters(method)

  contains

    ! Sets default values for parameters depending on the method
    subroutine set_default(method)
      implicit none
      character(*), intent(in) :: method

      logical :: l_unrestricted

      l_unrestricted = (Trim(method) == 'USCEXX')

      maxit = 30
      minit = 3
      maxdim_diis = 10
      fixmix_rate = 0.95d0
      l_fixmix = .false.

      ! Parameters that are the same in all methods
      iverbose = 0
      ldfit = .true.

      thr_overlap_oep = 1d-99
      thr_fai_oep = 1d-14
      thr_symmetry = 0d0

      thr_scf_energy = 1d-8
      thr_scf_density = 0d0

      thr_oep = 1d-6
      thr_solve = 1d-99
      csolvemeth = 'GTSVD'

      l_vref_fa = .true.
      l_homo = .false.
      l_vref_fa_sameab = .false.
      l_vh_oep = .false.

      l_plot_x_line = .false.
      l_plot_y_line = .false.
      l_plot_z_line = .false.

      l_plot_always = .false.

      l_test_pot = .false.

      ! Parameters that are not the same in restricted and unrestricted methods
      if (l_unrestricted) then
        irec_orb = 2200
        ifil_orb = 2
        irec_orb_out = 2201
        ifil_orb_out = 2
        l_spin_sym = .false.
      else
        irec_orb = 2100
        ifil_orb = 2
        irec_orb_out = 2101
        ifil_orb_out = 2
      end if

      l_swap = .false.
      l_swap2 = .false.
      iswap1 = 3
      iswap2 = 5
      l_frac_occ = .false.
      system = 'O'

    end subroutine set_default

    ! Prints parameters depending on the method
    subroutine print_parameters(method)
      implicit none

      character(*), intent(in) :: method
      logical :: l_unrestricted

      l_unrestricted = (Trim(method) == 'USCEXX')

      write (iout, '(1x,a,i4,a,i1)') &
        'Read Orbital Records from ', irec_orb, '.', ifil_orb
      write (iout, '(1x,a,i4,a,i1)') &
        'Write Orbital Records to ', irec_orb_out, '.', ifil_orb_out
      write (iout, '(1x,a,l2)') &
        'Use density fitting to evaluate the Coulomb and exchange integrals:', &
        ldfit
      write (iout, '(1x,a,i4)') 'Maximum number of iterations:', maxit
      write (iout, '(1x,a,i4)') 'Minimum number of iterations:', minit
      if (l_fixmix) then
        write (iout, '(1x,a)') 'Use a fixed mixing scheme instead of DIIS:'
        write (iout, '(1x,a, es10.3)') 'Ratio of old fock matrix: ', fixmix_rate
      else
        write (iout, '(1x,a)') 'DIIS is used:'
        write (iout, '(1x,a,i4)') 'Maximum size of DIIS matrix:', maxdim_diis
      end if
      write (iout, '(1x,a,es12.3)') &
        'Threshold for energy convergence:', thr_scf_energy
      write (iout, '(1x,a,es12.3)') &
        'Threshold for density convergence:', thr_scf_density
      write (iout, '(1x,a)') 'Thresholds for processing of OEP basis set:'
      write (iout, '(1x,a,es12.3)') 'thr_symmetry:', thr_symmetry
      write (iout, '(1x,a,1x,es12.3)') 'thr_overlap_oep: ', thr_overlap_oep
      write (iout, '(1x,a,1x,es12.3)') 'thr_fai_oep: ', thr_fai_oep
      write (iout, '(1x,a,1x,es12.3)') &
        'Threshold for eigenvalue difference in calculation of Lambda:', thr_oep
      write (iout, '(1x,a,1x,es12.3)') &
        'Threshold for eigenproblem solver of OEP equation:', thr_solve
      write (iout, '(1x,a)') 'Solver for eigenproblems: '//trim(csolvemeth)
      if (l_unrestricted) then
        write (iout, '(1x,a,l2)') &
          'Spin-symmetrization:', l_spin_sym
      end if
      write (iout, '(1x,a,l2)') &
        'Constucting Hartree potential in OEP basis:', l_vh_oep
      write (iout, '(1x,a,l2)') &
        'Fermi-Amaldi as reference potential:', l_vref_fa
      write (iout, '(1x,a,l2)') 'HOMO condition:', l_homo
      if ((l_unrestricted) .and. (l_vref_fa)) &
        write (iout, '(1x,a,l2)') &
        'The same reference potential for alpha and beta:', l_vref_fa_sameab

      if (l_plot_x_line) &
        write (iout, '(1x,a)') 'Plot exchange potential along x-axis'
      if (l_plot_y_line) &
        write (iout, '(1x,a)') 'Plot exchange potential along y-axis'
      if (l_plot_z_line) &
        write (iout, '(1x,a)') 'Plot exchange potential along z-axis'
      if ((l_plot_x_line) .or. (l_plot_y_line) .or. (l_plot_z_line)) &
        write (iout, '(1x,a,1x,l2)') 'Plot every iteration: ', l_plot_always

      write (iout, '(1x,a,l2)') 'l_swap:', l_swap
      write (iout, '(1x,a,l2)') 'l_swap2:', l_swap2
      if (l_swap) write (iout, '(1x,a,i3,i3)') "orbitals to swap:", iswap1, iswap2

      write (iout, '(1x,a,l2)') 'l_frac_occ:', l_frac_occ
      if (l_frac_occ) write (iout, '(1x,a,1x,a)') "system:", system

      write (iout, '(a)') ''

    end subroutine print_parameters
  end subroutine get_parameters_scexx

end module acfd_parameters
