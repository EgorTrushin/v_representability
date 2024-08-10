!
! Spin-unrestricted post-SCF dRPA (and sigma-functional) calculations
!
! Last refactoring: E. Trushin <e.v.trushin@gmail.com>
!
! The code is formatted with fprettify:
! fprettify -i 2 -l 80 -w 1 -s
!
subroutine acfd_urirpa
  use memory, only: memory_allocate, memory_release, memory_allocate_integer
  use common_cbas, only: ntg, ntdg, ntqg
  use common_tapes, only: iout
  use common_cmpp, only: nprocs, mpp_state, global_sum, &
                         initask, nextask_or_serial, endtask
  use molpro_options, only: molpro_pwd
  use outputResult, only: output_result
  use iso_fortran_env, only: dp => real64
  use acfd_wf, only: noa, nob, nva, nvb, dena, denb, orba, orbb, &
                     eiga, eigb, occa, occb
  use acfd_parameters, only: ldfit, nquad, iverbose, l_sigma, l_write_sigma, &
                             irec_orb, ifil_orb, get_parameters_rirpa
  use acfd_auxiliary_basis, only: naux_ri, naux_charge_ri, trans_mat_charge_ri
  use acfd_basis_process, only: basis_process_rirpa
  use acfd_integrals, only: build_fai
  use acfd_quadrature, only: quadrature_driver
  use acfd_sigma, only: cubic_spline_integr
  use acfd_utils, only: dfit_init, open_e_ref, get_x0, matrix_spectral, &
                        deallocate_arrays
  implicit none
  real(dp), parameter :: pi = 4d0*atan(1d0)
  integer, pointer :: ibase(:) ! pointers for memory stack management
  real(dp), pointer :: int_fai_ri_a(:, :, :), int_fai_ri_b(:, :, :)
  real(dp), pointer :: int_fai_a(:, :, :), int_fai_b(:, :, :)

  integer :: ifreq, iaux, imxs, iunit ! auxiliary variables

  real(dp) :: e_ref ! reference energy, i.e., EXX energy
  real(dp) :: e_corr_rpa, e_corr_sigma ! correlation energy
  real(dp) :: dedw_rpa, dedw_sigma ! auxiliary variable used in correlation energy calculation
  real(dp) :: e_tot_rpa, e_tot_sigma ! total energy
  real(dp), pointer :: omega(:), weight(:) ! nodes and weights for frequency integration
  real(dp), pointer :: x0(:, :), x0_aux(:, :) ! response matrix
  real(dp), pointer :: sigma(:, :) ! eigenvalues of response matrix
  real(dp), pointer :: evec(:, :) ! eigenvectors

  real(dp), pointer :: denc(:), dens(:) ! auxiliary arrays used during reading of densities

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

  real(dp), external :: dsum

  ts = wallcl()

  ! allocate integer as the starting point of memory stack for the subroutine
  ibase => memory_allocate_integer(1)

  ! print information about the code and authors
  write (iout, *) ''
  write (iout, *) 'PROGRAM * URIRPA'
  write (iout, *) '(unrestricted resolution of identity random-phase approximation)'
  write (iout, *)
  write (iout, *) 'Authors:'
  write (iout, *) 'A. Hesselmann, P. Bleiziffer, D. Schmidtel, J. Erhard, A. Thierbach, E. Trushin'
  write (iout, *)

  ! set and display parameters
  call get_parameters_rirpa("URIRPA")

  ! allocate required arrays
  orba => memory_allocate(ntqg)
  orbb => memory_allocate(ntqg)
  eiga => memory_allocate(ntg)
  eigb => memory_allocate(ntg)
  dena => memory_allocate(ntdg)
  denb => memory_allocate(ntdg)
  occa => memory_allocate(ntg)
  occb => memory_allocate(ntg)

  denc => memory_allocate(ntdg)
  dens => memory_allocate(ntdg)

  ! get values from the record
  call get_orb(orba, irec_orb, ifil_orb, 1, 1)
  call get_eig(eiga, irec_orb, ifil_orb, 1, 1)
  call get_den(denc, irec_orb, ifil_orb, 1, 1)
  call get_orb(orbb, irec_orb, ifil_orb, 2, 1)
  call get_eig(eigb, irec_orb, ifil_orb, 2, 1)
  call get_den(dens, irec_orb, ifil_orb, 2, 1)
  call get_occ(occa, irec_orb, ifil_orb, 1, 1)
  call get_occ(occb, irec_orb, ifil_orb, 2, 1)
  noa = int(dsum(ntg, occa, 1))
  nob = int(dsum(ntg, occb, 1))
  nva = ntg-noa
  nvb = ntg-nob
  dena = 0.5d0*(denc+dens)
  denb = 0.5d0*(denc-dens)
  call memory_release(denc)

  call basis_number('RI', naux_ri)

  ! initialize density fitting, if required
  if (ldfit) call dfit_init()

  ! get reference_energy
  e_ref = open_e_ref(dena, denb, orba, orbb, ldfit)

  tf_prep = wallcl()

  ! construct 3-index integral (basis,virt,occ) for alpha spin
  ts_3idx = wallcl()
  int_fai_a => memory_allocate(naux_ri, nva, noa)
  call build_fai('RI', 'J', orba, noa, nva, naux_ri, int_fai_a)

  if (nob /= 0) then

    ! construct 3-index integral (basis,virt,occ) for beta spin
    int_fai_b => memory_allocate(naux_ri, nvb, nob)
    call build_fai('RI', 'J', orbb, nob, nvb, naux_ri, int_fai_b)
    tf_3idx = wallcl()

    ! process basis set for alpha and beta spin simultaniously
    ts_ri_prep = wallcl()
    call basis_process_rirpa(trans_mat_charge_ri, int_fai_a, int_fai_b)
    tf_ri_prep = wallcl()

  else

    tf_3idx = wallcl()
    ! process basis set for alpha spin only
    ts_ri_prep = wallcl()
    call basis_process_rirpa(trans_mat_charge_ri, int_fai_a)
    tf_ri_prep = wallcl()

  end if

  naux_charge_ri = trans_mat_charge_ri%m

  ! if high level of verbosity specified, print information
  if (iverbose > 1) write (iout, *) 'naux_charge_ri:', naux_charge_ri

  ts_3idx_ri = wallcl()
  ! construct 3-index integral (basis,virt,occ) for alpha spin and transform
  int_fai_ri_a => memory_allocate(naux_charge_ri, nva, noa)
  call dgemm_x('T', 'N', naux_charge_ri, nva*noa, naux_ri, 1d0, &
               trans_mat_charge_ri%mat, naux_ri, int_fai_a, naux_ri, &
               0d0, int_fai_ri_a, naux_charge_ri)

  ! construct 3-index integral (basis,virt,occ) for beta spin and transform
  if (nob /= 0) then
    int_fai_ri_b => memory_allocate(naux_charge_ri, nvb, nob)
    call dgemm_x('T', 'N', naux_charge_ri, nob*nvb, naux_ri, 1d0, &
                 trans_mat_charge_ri%mat, naux_ri, int_fai_b, naux_ri, &
                 0d0, int_fai_ri_b, naux_charge_ri)
  end if
  tf_3idx_ri = wallcl()

  ! preparation before frequency loop
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

    ! calculate response function
    ! for RPA MPI parallelization over virtual orbitals is used
    x0 = 0d0
!    call get_x0(x0, eiga, omega(ifreq), noa, nva, naux_charge_ri, &
!                int_fai_ri_a, 'RIRPA')
    if (nob /= 0) then
      x0_aux => memory_allocate(naux_charge_ri, naux_charge_ri)
      x0_aux = 0d0
!      call get_x0(x0_aux, eigb, omega(ifreq), nob, nvb, naux_charge_ri, &
!                  int_fai_ri_b, 'RIRPA')
      x0 = x0+x0_aux
      call memory_release(x0_aux)
    end if

    x0 = x0/2d0

    call matrix_spectral(x0, evec, sigma(ifreq, :), naux_charge_ri)

    dedw_rpa = 0d0
    dedw_sigma = 0d0
    do iaux = 1, naux_charge_ri

      dedw_rpa = dedw_rpa+log(1d0+sigma(ifreq, iaux))-sigma(ifreq, iaux)

      ! if sigma-functional calculation is active, calculate spline-correction
      if ((l_sigma) .and. (sigma(ifreq, iaux) > 0d0)) &
        dedw_sigma = dedw_sigma-cubic_spline_integr(sigma(ifreq, iaux))

    end do

    ! calculate dRPA correlation energy
    e_corr_rpa = e_corr_rpa+weight(ifreq)*dedw_rpa/(2d0*pi)
    if (l_sigma) e_corr_sigma = e_corr_sigma &
                                +weight(ifreq)*dedw_rpa/(2d0*pi) &
                                +weight(ifreq)*dedw_sigma/(2d0*pi)

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
    open (newunit=iunit, file=molpro_pwd//'sigma.dat',status='new')
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

end subroutine acfd_urirpa
