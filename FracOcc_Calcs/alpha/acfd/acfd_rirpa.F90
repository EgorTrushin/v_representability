!
! Spin-restricted post-SCF dRPA (and sigma-functional) calculations
!
! Last refactoring: E. Trushin <e.v.trushin@gmail.com>
!
! The code is formatted with fprettify:
! fprettify -i 2 -l 80 -w 1 -s
!
subroutine acfd_rirpa
  use memory, only: memory_allocate, memory_release, memory_allocate_integer
  use common_cbas, only: ntg, ntdg, ntqg
  use common_cmpp, only: nprocs, mpp_state, global_sum, &
                         initask, nextask_or_serial, endtask
  use common_tapes, only: iout
  use molpro_options, only: molpro_pwd
  use outputResult, only: output_result
  use iso_fortran_env, only: dp => real64
  use acfd_auxiliary_basis, only: naux_ri, naux_charge_ri, &
                                  trans_mat_charge_ri
  use acfd_basis_process, only: basis_process_rirpa
  use acfd_integrals, only: build_fai
  use acfd_parameters, only: ldfit, iverbose, nquad, l_sigma, l_write_sigma, &
                             vc_scaling, irec_orb, ifil_orb, &
                             get_parameters_rirpa
  use acfd_quadrature, only: quadrature_driver
  use acfd_sigma, only: cubic_spline_integr
  use acfd_utils, only: dfit_init, closed_e_ref, get_x0, matrix_spectral, &
                        deallocate_arrays
  use acfd_wf, only: no, nv, orb, den, eig
  implicit none
  real(dp), parameter :: pi = 4d0*atan(1d0)

  integer, pointer :: ibase(:) ! pointer for stack management

  real(dp), pointer :: int_fai_ri(:, :, :) ! aux,virt,occ integrals in RI basis
  real(dp), pointer :: int_fai(:, :, :) ! aux,virt,occ integrals

  integer :: ifreq, iaux, imxs, iunit ! auxiliary variables

  real(dp) :: e_ref ! reference energy, i.e., EXX energy
  real(dp) :: e_corr_rpa, e_corr_sigma ! correlation energy
  real(dp) :: dedw_rpa, dedw_sigma ! auxiliary variable used in correlation energy calculation
  real(dp) :: e_tot_rpa, e_tot_sigma ! total energy
  real(dp) :: r_nelec

  integer :: nov ! number of occupied times number of virtual orbitals

  real(dp), pointer :: omega(:), weight(:) ! nodes and weights for frequency integration
  real(dp), pointer :: x0(:, :), evec(:, :), sigma(:, :) ! response matrix, eigenvectors and eigenvalues

  ! for measurements of timings
  real(dp) :: ts
  real(dp) :: tf_prep
  real(dp) :: ts_3idx, tf_3idx
  real(dp) :: ts_ri_prep, tf_ri_prep
  real(dp) :: ts_3idx_ri, tf_3idx_ri
  real(dp) :: ts_omega_loop, tf_omega_loop
  real(dp), external :: wallcl

  ! parallelization
  real(dp) :: tasks
  integer :: ngrp, minbt, maxbt

  ! dummy/auxiliary
  integer :: i_dummy1, i_dummy2
  character(16) :: c_dummy

  ts = wallcl()

  ! allocate integer as the starting point of memory stack used by present subroutine
  ibase => memory_allocate_integer(1)

  ! print information about the code and authors
  write (iout, *) ''
  write (iout, *) 'PROGRAM * RIRPA  (restricted resolution of identity random-phase approximation)'
  write (iout, *)
  write (iout, *) 'Authors:'
  write (iout, *) 'A. Hesselmann, P. Bleiziffer, D. Schmidtel, J. Erhard, A. Thierbach, E. Trushin'
  write (iout, *)

  ! set and display parameters
  call get_parameters_rirpa('RIRPA')

  ! allocate required arrays
  orb => memory_allocate(ntqg)
  eig => memory_allocate(ntg)
  den => memory_allocate(ntdg)

  ! get values from the record
  call get_orb(orb, irec_orb, ifil_orb, 1, 1)
  call get_eig(eig, irec_orb, ifil_orb, 1, 1)
  call get_den(den, irec_orb, ifil_orb, 1, 1)
  call getvar('NELEC', r_nelec, c_dummy, i_dummy1, i_dummy2, 1, 1)

  no = int(r_nelec/2d0)
  nv = ntg-no

  call basis_number('RI', naux_ri)
  nov = no*nv

  ! initialize density fitting, if required
  if (ldfit) call dfit_init()

  ! get reference_energy
  e_ref = closed_e_ref(den, orb, ldfit)

  tf_prep = wallcl()

  ! construct 3-index integral (basis,virt,occ)
  ts_3idx = wallcl()
  int_fai => memory_allocate(naux_ri, nv, no)
  call build_fai('RI', 'J', orb, no, nv, naux_ri, int_fai)
  tf_3idx = wallcl()

  ! preprocessing of RI basis set
  ts_ri_prep = wallcl()
  call basis_process_rirpa(trans_mat_charge_ri, int_fai)
  tf_ri_prep = wallcl()

  naux_charge_ri = trans_mat_charge_ri%m

  ! if high level of verbosity specified, print information
  if (iverbose > 1) write (iout, *) 'naux_ri =', naux_ri
  if (iverbose > 1) write (iout, *) 'naux_charge_ri =', naux_charge_ri

  ! construct 3-index integral (basis,virt,occ) and transform
  ts_3idx_ri = wallcl()
  int_fai_ri => memory_allocate(naux_charge_ri, nv, no)
  call dgemm_x('T', 'N', naux_charge_ri, nv*no, naux_ri, 1d0, &
               trans_mat_charge_ri%mat, naux_ri, int_fai, naux_ri, &
               0d0, int_fai_ri, naux_charge_ri)
  tf_3idx_ri = wallcl()

  ! allocate required arrays before frequency integration
  x0 => memory_allocate(naux_charge_ri, naux_charge_ri)
  sigma => memory_allocate(nquad, naux_charge_ri)
  sigma = 0d0  ! mandatory for mpi/mpp
  evec => memory_allocate(naux_charge_ri, naux_charge_ri)

  ! allocate arrays for nodes and weights of frequency integration
  omega => memory_allocate(nquad)
  weight => memory_allocate(nquad)

  ! generate grid for frequency integration
  call quadrature_driver(omega, weight, nquad)

  ! print header for additional information which will be printed inside
  ! frequency loop, if high level of verbosity was specified
  if (iverbose > 0) then
    if (l_sigma) then
      write (iout, '(/,4x,a,4x,a,11x,a,5x,a,3x,a,4x,a)') &
        'ifreq', 'omega', 'dE/domega RPA', 'Ecorr RPA', &
        'dE/domega SIGMA', 'Ecorr SIGMA'
    else
      write (iout, '(/,4x,a,4x,a,11x,a,5x,a)') &
        'ifreq', 'omega', 'dE/domega RPA', 'Ecorr RPA'
    end if
  end if

  ! mpi/mpp init
  if (nprocs > 1) then
    tasks = dble(nquad)
    call initask('Omega integration', 1, 0, tasks, ngrp, minbt, maxbt)
    mpp_state = 2
  end if

  ! initialize variables for correlation energy
  e_corr_rpa = 0d0
  e_corr_sigma = 0d0

  ! loop over frequencies
  ts_omega_loop = wallcl()
  omega_loop: do ifreq = 1, nquad
    if (.not. nextask_or_serial()) cycle ! mpi/mpp

    ! calculate response function x0
!    call get_x0(x0, eig, omega(ifreq), no, nv, naux_charge_ri, int_fai_ri, &
!                'RIRPA')

    ! diagonalize x0
    call matrix_spectral(x0, evec, sigma(ifreq, :), naux_charge_ri)

    dedw_rpa = 0d0
    dedw_sigma = 0d0
    do iaux = 1, naux_charge_ri

      dedw_rpa = dedw_rpa+ &
            log(1d0+vc_scaling*sigma(ifreq, iaux))/vc_scaling-sigma(ifreq, iaux)

      ! if sigma-functional calculation is active, calculate spline-correction
      if ((l_sigma) .and. (sigma(ifreq, iaux) > 0d0)) &
        dedw_sigma = dedw_sigma-cubic_spline_integr(sigma(ifreq, iaux))

    end do

    e_corr_rpa = e_corr_rpa+weight(ifreq)*dedw_rpa/(2d0*pi)
    if (l_sigma) e_corr_sigma = e_corr_sigma &
              +weight(ifreq)*dedw_rpa/(2d0*pi)+weight(ifreq)*dedw_sigma/(2d0*pi)

    ! print additional information, if high level of verbosity was specified
    if (iverbose > 0) then
      if (l_sigma) then
        write (iout, '(1x,i5,f14.6,4x,5f16.8)') &
          ifreq, omega(ifreq), dedw_rpa, e_corr_rpa, &
          dedw_rpa+dedw_sigma, e_corr_sigma
      else
        write (iout, '(1x,i5,f14.6,4x,2f16.8)') &
          ifreq, omega(ifreq), dedw_rpa, e_corr_rpa
      end if
    end if

  end do omega_loop
  tf_omega_loop = wallcl()

  ! mpi/mpp term
  if (nprocs > 1) then
    call endtask
    mpp_state = 1
    call global_sum(e_corr_rpa, 1, '+')
    if (l_sigma) call global_sum(e_corr_sigma, 1, '+')
    call global_sum(sigma, nquad*naux_charge_ri, '+')
  end if

  ! write sigma.dat file, if required
  if (l_write_sigma) then
    open (newunit=iunit, file=molpro_pwd//'sigma.dat', status='new')
    write (iunit, '(es30.15)') e_ref
    write (iunit, *) nquad, naux_charge_ri
    do ifreq = 1, nquad
      write (iunit, '(es30.15)') weight(ifreq)
    end do
    do ifreq = 1, nquad
      do iaux = 1, naux_charge_ri
        write (iunit, '(es30.15)') sigma(ifreq, iaux)
      end do
    end do
    close (iunit)
  end if

  ! calculate and print total energy
  ! provide information about calculated energy to main program
  e_tot_rpa = e_ref+e_corr_rpa
  write (iout, *)
  call output_result('Reference', 'Energy', &
                     e_ref, showstate=.false., principal=.false.)
  if (l_sigma) then
    call output_result('RPA', 'Correlation Energy', &
                       e_corr_rpa, showstate=.false., principal=.false.)
    call output_result('RPA', 'Total Energy', &
                       e_tot_rpa, showstate=.false., principal=.false.)
    e_tot_sigma = e_ref+e_corr_sigma
    call output_result('SIGMA', 'Correlation Energy', &
                       e_corr_sigma, showstate=.false., principal=.false.)
    call output_result('SIGMA', 'Total Energy', &
                       e_tot_sigma, showstate=.false., principal=.true.)
    call setvar('ENERGY', e_tot_sigma, 'AU', 1, 1, imxs, 0)
    call output_energy(e_tot_sigma, 'SIGMA')
  else
    call output_result('RPA', 'Correlation Energy', &
                       e_corr_rpa, showstate=.false., principal=.false.)
    call output_result('RPA', 'Total Energy', &
                       e_tot_rpa, showstate=.false., principal=.true.)
    call setvar('ENERGY', e_tot_rpa, 'AU', 1, 1, imxs, 0)
    call output_energy(e_tot_rpa, 'RPA')
  end if

  ! free memory
  call memory_release(ibase)
  call deallocate_arrays()

  if (iverbose > 0) then
    write (iout, '(/,1x,a,f12.2)') &
      "Time for preparations (sec):", tf_prep-ts
    write (iout, '(1x,a,f12.2)') &
      "Time for construction of 3-idx integrals (sec):", tf_3idx-ts_3idx
    write (iout, '(1x,a,f12.2)') &
      "Time for RI basis preprocessing (sec):", tf_ri_prep-ts_ri_prep
    write (iout, '(1x,a,f12.2)') &
     "Time for construction of RI 3-idx integrals (sec):", tf_3idx_ri-ts_3idx_ri
    write (iout, '(1x,a,f12.2)') &
      "Time for frequency loop (sec):", tf_omega_loop-ts_omega_loop
    write (iout, '(1x,a,f12.2)') &
      "Total execution time for RIRPA (sec):", wallcl()-ts
  end if

end subroutine acfd_rirpa
