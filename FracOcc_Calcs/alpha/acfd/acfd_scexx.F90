!
! Driver for spin-restricted optimized effective potential
! exact-exchange calculations
!
! Last refactoring: E. Trushin <e.v.trushin@gmail.com>
!
! The code is formatted with fprettify:
! fprettify -i 2 -l 80 -w 1 -s
!
subroutine acfd_scexx
  use common_cbas, only: ntdg, ntqg, ntg
  use common_code, only: ho1, rs
  use common_tapes, only: iout
  use common_molen, only: ekern
  use common_scfopt, only: maxdis
  use memory, only: memory_allocate, memory_release, memory_allocate_integer
  use outputResult, only: output_result
  use acfd_opm, only: scexx_opm
  use acfd_parameters, only: ldfit, maxit, minit, maxdim_diis, irec_orb, &
                             ifil_orb, thr_scf_energy, thr_scf_density, &
                             l_fixmix, fixmix_rate, &
                             irec_orb_out, ifil_orb_out, iverbose, &
                             irec_orb, ifil_orb, l_test_pot, &
                             get_parameters_scexx
  use acfd_plot, only: plot_potential, vref_aux_up
  use acfd_pot, only: vj, vx, fock, old_fock, h0, vsol
  use acfd_pot_test, only: scexx_pot_test
  use acfd_utils, only: dfit_init, print_eig, deallocate_arrays
  use acfd_wf, only: den, orb, eig, smh, no, nv, occ
  use iso_fortran_env, only: dp => real64
  implicit none

  integer, pointer :: ibase(:), ibase2(:) ! Pointers for stack management

  integer :: iter ! loop counter: number of current iteration

  real(dp) :: twoel_energy ! two electron energy
  real(dp) :: h0_energy ! single electron energy
  real(dp) :: exchange_energy
  real(dp) :: coulomb_energy
  real(dp) :: total_energy ! total energy
  real(dp) :: total_energy_old ! total energy of last iteration
  real(dp) :: diff_energy ! energy difference between two iterations
  real(dp) :: diff_density ! density difference between two iterations

  real(dp), pointer :: denc_old(:) ! old charge density for diff_density calculation

  ! variables for timing
  real(dp) :: ts

  ! DIIS parameters
  integer :: ndim, ndel
  real(dp) :: bfak

  ! variables for density fitting
  integer :: icoulfil = 4, icoulrec = 5000

  ! auxiliary variables
  real(dp), pointer :: r1_tmp1(:)
  real(dp), pointer :: fock_tri(:)
  real(dp), pointer :: dipole(:)
  integer :: imxs
  real(dp) :: r_number_of_electrons
  character(16) :: c_dummy
  integer :: i_dummy1, i_dummy2
  integer :: naux_oep

  real(dp), external :: spur ! trace of the matrix product of to triagonal matrices
  real(dp), external :: wallcl ! timing measurment
  real(dp), external :: ddot_x ! scalar product

  logical, external :: exx_eris_exist ! check whether ERIS record exists or not

  ts = wallcl()

  ibase => memory_allocate_integer(1) ! create a pointer to the bottom of the stack used in this driver

  write (iout, *)
  write (iout, *) 'PROGRAM * SCEXX      (restricted self-consistent exact-exchange)'
  write (iout, *)
  write (iout, *) 'Authors:'
  write (iout, *) 'A. Hesselmann, P. Bleiziffer, D. Schmidtel, J. Erhard, A. Thierbach, E. Trushin'
  write (iout, *)

  ! set and display parameters
  call get_parameters_scexx("SCEXX")

  if ((.not. ldfit) .and. (.not. exx_eris_exist(1))) then
    write (iout, *) &
      'Warning! Calculations with ldfit=0 and GDIRECT are not supported'
    write (iout, *) &
      'Do not worry about this warning if you are doing parallel calculation'
    write (iout, *)
  end if

  ! set DIIS paramters
  maxdis = maxdim_diis
  ndim = 1
  ndel = 0
  bfak = 1d0

  ! allocate arrays for calculations
  orb => memory_allocate(ntqg)
  eig => memory_allocate(ntg)
  den => memory_allocate(ntdg)
  occ => memory_allocate(ntg)
  fock => memory_allocate(ntqg)
  old_fock => memory_allocate(ntdg)
  old_fock = 0d0
  vx => memory_allocate(ntdg)
  vj => memory_allocate(ntdg)
  vj = 0d0
  call basis_number('OEP', naux_oep)
  vsol => memory_allocate(naux_oep)
  h0 => memory_allocate(ntdg)
  h0 = 0d0
  call les(h0, ntdg, 1, ho1) ! get h0
  smh => memory_allocate(ntqg)
  call les(smh, ntdg, 1, rs) ! get smh
  call expan(smh, smh, 1, 1, 1) ! expand smh

  ! get values from the record
  call get_orb(orb, irec_orb, ifil_orb, 1, 1)
  call get_eig(eig, irec_orb, ifil_orb, 1, 1)
  call get_den(den, irec_orb, ifil_orb, 1, 1)
  call get_occ(occ, irec_orb, ifil_orb, 1, 1)

  call getvar('NELEC', r_number_of_electrons, c_dummy, i_dummy1, i_dummy2, 1, 1)
  no = int(r_number_of_electrons/2d0)
  nv = ntg-no

  ! get initial charge density for diff_density calculation
  denc_old => memory_allocate(ntdg)
  denc_old = den

  if (ldfit) call dfit_init() ! initialize density fitting if required

  write (iout, '(/,3x,a,8x,a,12x,a,7x,a,7x,a,9x,a)') &
    'ITER', 'ENERGY', 'EDIFF', 'DDIFF', 'DIIS', 'TIME'

  ! start scf cycle
  scf_loop: do iter = 0, maxit

    ibase2 => memory_allocate_integer(1) ! create pointer to start of scf loop

    ! calculate nonlocal exchange and coulomb matrix
    vj = 0d0
    vx = 0d0
    if (.not. ldfit) then
      r1_tmp1 => memory_allocate(1)
      call hfma(den, den, vj, vx, 1, 1, 0, .true., .true., &
                r1_tmp1, .false., .false.)
      call memory_release(r1_tmp1)
    else
      r1_tmp1 => memory_allocate(1)
      call reservem(ntdg, 0, icoulfil, icoulrec, 0, 0, 'JOP')
      call dfexch(den, vx, vj, orb, .true., .false., &
                  r1_tmp1(1), icoulrec, icoulfil)
      call lesw(vj, ntdg, icoulfil, icoulrec, 0)
      call truncatem(icoulfil, icoulrec, 0)
      call daxpy_x(ntdg, -1d0, vj, 1, vx, 1)
      call memory_release(r1_tmp1)
    end if

    ! calculate two-electron energy
    twoel_energy = 0.5d0*(spur(den, vx)+spur(den, vj))
    h0_energy = spur(den, h0)
    total_energy = h0_energy+twoel_energy+ekern
    exchange_energy = 0.5d0*spur(den, vx)
    coulomb_energy = 0.5d0*spur(den, vj)
    if (iverbose > 2) then
      write (iout, '(1x,a,f24.12)') 'E_kin  =', h0_energy
      write (iout, '(1x,a,f24.12)') 'E_x    =', exchange_energy
      write (iout, '(1x,a,f24.12)') 'E_Coul =', coulomb_energy
    end if

    call scexx_opm(iter) ! generate local exchange potential

    ! print output
    if (iter .eq. 0) then
      diff_energy = total_energy
      write (iout, '(1x,i5,2x,f18.10,38x,f10.2)') &
      & iter, total_energy, wallcl()-ts
    else
      diff_energy = total_energy-total_energy_old
      write (iout, '(1x,i5,2x,f18.10,2x,2e12.3,2x,i5,5x,f10.2)') &
            & iter, total_energy, diff_energy, diff_density, ndim-1, wallcl()-ts
    end if

    ! print eigenvalues
    if (iverbose > 0) then
      call print_eig(eig, min(ntg, no+10), 'SCEXX eigenvalues')
      write (iout, '(1x,a,1x,e12.5)') 'LUMO-HOMO:', (eig(no+1)-eig(no))
      write (iout, *)
    end if

    total_energy_old = total_energy

    ! check for convergence
    if (iter .ge. minit) then
      if ((thr_scf_energy .le. 0d0) .or. &
          abs(diff_energy) .lt. thr_scf_energy) then
        if ((thr_scf_density .le. 0d0) .or. &
            abs(diff_density) .lt. thr_scf_density) then
          exit scf_loop
        end if
      end if
    end if

    if (iter .eq. maxit) exit scf_loop

    ! build fock matrix
    fock_tri => memory_allocate(ntdg)
    fock_tri = 0d0
    fock_tri = h0+vj+vx

    if (l_fixmix) then
      ! fixed mixing of fock matrices instead of DIIS
      if (iter > 0) then
        fock_tri = (1d0-fixmix_rate)*fock_tri+old_fock*fixmix_rate
      else
        fock_tri = h0+vj+vx
      end if
      old_fock = fock_tri
    else

      ! DIIS for accelerating convergence
      call diiss(den, fock_tri, ndim, ndel, bfak, 1)
      if (ndim .gt. maxdim_diis) then
        ndim = 1
        ndel = 0
        bfak = 1d0
      end if

    end if

    ! calculate orbitals and eigenvalues, expand fock matrix
    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp1 = 0d0
    fock = 0d0
    call expan(fock_tri, fock, 1, 1, 1)
    call mmult(smh, fock, r1_tmp1, 1, 1, 0, -1, 0, 0)
    call mmult(r1_tmp1, smh, fock, 1, 1, 0, 0, 0, 0)
    call diag(fock, orb, eig, 0)
    call mmult(smh, orb, r1_tmp1, 1, 1, 0, -1, 0, 0)
    orb = r1_tmp1

    ! calculate densities
    r1_tmp1 = 0d0
    call dgemm_x('N', 'T', ntg, ntg, no, 1d0, &
                 orb, ntg, &
                 orb, ntg, &
                 0d0, r1_tmp1, ntg)
    call reduc(r1_tmp1, den, 1, 1, 1d0, 1)

    call memory_release(fock_tri)

    ! calculate density difference for convergence criteria
    r1_tmp1 => memory_allocate(ntdg)
    r1_tmp1 = 0d0
    r1_tmp1 = den-denc_old
    diff_density = sqrt(ddot_x(ntdg, r1_tmp1, 1, r1_tmp1, 1)/dble(ntdg))
    denc_old = den

    call memory_release(ibase2)

  end do scf_loop

  if (iter .lt. maxit) then
    write (iout, '(/,1x,a,i5,a)') 'SCF converged in ', iter, ' iterations'
  else
    write (iout, '(/,1x,a,i5,a)') 'SCF NOT converged in ', iter, ' iterations'
  end if

  ! test potential, if required
  if (l_test_pot) then
    if (ldfit) then
      write (iout, '(/,1x, a)') &
        'WARNING: test for potential is skipped. It works only with ldfit=0.'
    else
      call scexx_pot_test(6, 1d0)
    end if
  end if

  ! plot/write final potential, if required (cases inside)
  call plot_potential(orb, den, vsol, naux_oep, 'final', 'OEP', no, 'vx', &
                      vref_aux_up)

  ! create SCEXX records
  irec_orb = irec_orb_out
  ifil_orb = ifil_orb_out
  call reserve_dump(irec_orb, ifil_orb, 'SCEXX', ntqg*2+ntg*2+ntdg)
  call write_orb(orb, 1, 'CANONICAL')
  call write_den(den, 1, 'CHARGE')
  call write_eig(eig, 1, 'CANONICAL')
  call write_occ(occ, 1, 'OCC')
  call write_fock(fock, 1, 'TOTAL')
  call flush_dump
  write (iout, '(/,1x,a,i8,a,i1)') &
    'SCEXX orbitals written to record', irec_orb, '.', ifil_orb

  ! print energy
  write (iout, *) ''
  call output_result('SCEXX', 'Total energy', &
                     total_energy, showstate=.false., principal=.true.)
  call output_result('SCEXX', 'One-electron energy', &
                     h0_energy, showstate=.false., principal=.false.)
  call output_result('SCEXX', 'Two-electron energy', &
                     twoel_energy, showstate=.false., principal=.false.)
  call output_result('SCEXX', 'Nuclear energy', &
                     ekern, showstate=.false., principal=.false.)
  call output_result('SCEXX', 'Coulomb energy', &
                     coulomb_energy, showstate=.false., principal=.false.)
  call output_result('SCEXX', 'Exchange energy', &
                     exchange_energy, showstate=.false., principal=.false.)

  ! calculate and print dipole moment
  dipole => memory_allocate(3)
  call dipole_moment(dipole, den)
  write (iout, *) ''
  call output_result('SCEXX', 'Dipole moment', dipole(1:3), &
       & showstate=.false., numberformat='3f15.8', debye=.TRUE.)
  call memory_release(dipole)

  ! print eigenvalues and LUMO-HOMO
  call print_eig(eig, min(ntg, no+10), 'SCEXX eigenvalues')
  write (iout, '(1x,a,1x,e12.5)') 'HOMO:', eig(no)
  write (iout, '(1x,a,1x,e12.5)') 'LUMO:', eig(no+1)
  write (iout, '(1x,a,1x,e12.5)') 'LUMO-HOMO:', (eig(no+1)-eig(no))

  ! set energy variables
  call setvar('ENERGY', total_energy, 'AU', 1, 1, imxs, 0)
  call output_energy(total_energy, 'SCEXX')

  ! release memory
  call memory_release(ibase)
  call deallocate_arrays()

end subroutine acfd_scexx
