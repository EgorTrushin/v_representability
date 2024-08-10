!
! Driver for spin-unrestricted optimized effective potential
! exact-exchange calculations
!
! Last refactoring: E. Trushin <e.v.trushin@gmail.com>
!
! The code is formatted with fprettify:
! fprettify -i 2 -l 80 -w 1 -s
!
subroutine acfd_uscexx
  use common_cbas, only: ntdg, ntqg, ntg
  use common_code, only: ho1, rs, sm
  use common_tapes, only: iout
  use common_molen, only: ekern
  use common_scfopt, only: maxdis
  use memory, only: memory_allocate, memory_release, memory_allocate_integer
  use outputResult, only: output_result
  use acfd_opm, only: uscexx_opm, uscexx_opm_spin_symmetrized
  use acfd_parameters, only: ldfit, maxit, minit, maxdim_diis, irec_orb, &
                             ifil_orb, thr_scf_energy, thr_scf_density, &
                             l_fixmix, fixmix_rate, &
                             irec_orb_out, ifil_orb_out, iverbose, &
                             l_spin_sym, l_test_pot, get_parameters_scexx, &
                             l_swap, l_swap2, iswap1, iswap2, l_frac_occ, &
                             system
  use acfd_plot, only: plot_potential, vref_aux_upa, vref_aux_upb, vref_aux_up
  use acfd_pot, only: vj, h0, vxa, vxb, focka, fockb, old_focka, old_fockb, &
                      vsol, vsola, vsolb, vxb1, vxb2
  use acfd_pot_test, only: uscexx_pot_test
  use acfd_utils, only: dfit_init, print_eig, ssquared, deallocate_arrays, uscexx_pots
  use acfd_wf, only: dena, denb, orba, orbb, eiga, eigb, noa, nob, nva, nvb, &
                     occa, occb, smh, sh_a, orbb_frac, denb_frac, int_all, occb1, occb2
  use ksinv_utils, only: den_from_orb
  use iso_fortran_env, only: dp => real64
  implicit none

  integer, pointer :: ibase(:), ibase2(:) ! pointers for stack management

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
  real(dp) ::  ts

  ! DIIS parameters
  integer :: ndim, ndel
  real(dp) :: bfak

  real(dp), external :: spur ! trace of the matrix product of to triagonal matrices
  real(dp), external :: wallcl ! timing measurment
  real(dp), external :: ddot_x ! scalar product
  real(dp), external :: dsum

  ! auxiliaray variables
  real(dp), pointer :: focka_tri(:), fockb_tri(:)
  real(dp), pointer :: r1_tmp1(:)
  real(dp), pointer :: dipole(:)
  integer :: imxs
  integer :: naux_oep
  real(dp), pointer :: denc(:), dens(:)
  logical, external :: exx_eris_exist ! check whether ERIS record exists or not

  character(32) :: system_aux

  ts = wallcl()

  ! create a pointer to the bottom of the stack used in this driver
  ibase => memory_allocate_integer(1)

  write (iout, *)
  write (iout, *) 'PROGRAM * USCEXX      (unrestricted self-consistent exact-exchange)'
  write (iout, *)
  write (iout, *) 'Authors:'
  write (iout, *) 'A. Hesselmann, P. Bleiziffer, D. Schmidtel, J. Erhard, A. Thierbach, E. Trushin'
  write (iout, *)

  ! set and display parameters
  call get_parameters_scexx('USCEXX')

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

  if (l_frac_occ) then
    orbb_frac => memory_allocate(ntqg)
    denb_frac => memory_allocate(ntdg)
    vxb1 => memory_allocate(ntdg)
    vxb2 => memory_allocate(ntdg)
    occb1 => memory_allocate(ntg)
    occb2 => memory_allocate(ntg)
  end if

  call basis_number('OEP', naux_oep)
  int_all => memory_allocate(naux_oep, ntg, ntg)

  ! allocate arrays for calculations
  orba => memory_allocate(ntqg)
  orbb => memory_allocate(ntqg)
  eiga => memory_allocate(ntg)
  eigb => memory_allocate(ntg)
  dena => memory_allocate(ntdg)
  denb => memory_allocate(ntdg)
  occa => memory_allocate(ntg)
  occb => memory_allocate(ntg)
  focka => memory_allocate(ntqg)
  fockb => memory_allocate(ntqg)
  old_focka => memory_allocate(ntdg)
  old_fockb => memory_allocate(ntdg)
  vxa => memory_allocate(ntdg)
  vxb => memory_allocate(ntdg)
  vj => memory_allocate(ntdg)
  call basis_number('OEP', naux_oep)
  if (l_spin_sym) then
    vsol => memory_allocate(naux_oep)
    if (l_test_pot) then
      vsola => memory_allocate(naux_oep)
      vsolb => memory_allocate(naux_oep)
    end if
  else
    vsola => memory_allocate(naux_oep)
    vsolb => memory_allocate(naux_oep)
  end if
  h0 => memory_allocate(ntdg)
  call les(h0, ntdg, 1, ho1) ! get h0
  ! calculate smh and sh_a
  smh => memory_allocate(ntqg)
  sh_a => memory_allocate(ntqg)
  r1_tmp1 => memory_allocate(ntqg)
  call les(smh, ntdg, 1, rs)
  call expan(smh, smh, 1, 1, 1)
  call les(r1_tmp1, ntdg, 1, sm)
  call expan(r1_tmp1, r1_tmp1, 1, 1, 1)
  call mulq(r1_tmp1, smh, sh_a)
  call memory_release(r1_tmp1)

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

  ! get initial charge density for diff_density calculation
  denc_old => memory_allocate(ntdg)
  denc_old = dena+denb

  ! initialize density fitting, if required
  if (ldfit) call dfit_init()

  write (iout, '(/,3x,a,8x,a,12x,a,7x,a,7x,a,9x,a)') &
    'ITER', 'ENERGY', 'EDIFF', 'DDIFF', 'DIIS', 'TIME'

  ! start scf cycle
  scf_loop: do iter = 0, maxit

    ! Create pointer to start of scf loop
    ibase2 => memory_allocate_integer(1)

    ! Swap orbitals and eigenvalues, if required
    if (l_swap) then
      call uscexx_swap(orba, orbb, eiga, eigb, dena, denb, ntqg, ntdg, &
                       ntg, noa, nob, iswap1, iswap2)
    end if

    call uscexx_pots(dena, denb, vj, vxa, vxb)

    ! calculate two-electron energy
    twoel_energy = 0.5d0*(spur(dena, vj)+spur(denb, vj) &
                          +spur(dena, vxa)+spur(denb, vxb))
    h0_energy = spur(dena, h0)+spur(denb, h0)
    total_energy = h0_energy+twoel_energy+ekern
    coulomb_energy = 0.5d0*(spur(dena, vj)+spur(denb, vj))
    exchange_energy = 0.5d0*(spur(dena, vxa)+spur(denb, vxb))

    ! print output
    if (iter .eq. 0) then
      diff_energy = total_energy
      write (iout, '(1x,i5,2x,f18.10,38x,f10.2)') &
        iter, total_energy, wallcl()-ts
    else
      diff_energy = total_energy-total_energy_old
      write (iout, '(1x,i5,2x,f18.10,2x,2e12.3,2x,i5,5x,f10.2)') &
        iter, total_energy, diff_energy, diff_density, &
        ndim-1, wallcl()-ts
    end if

    if (l_swap2) then
      call uscexx_swap(orba, orbb, eiga, eigb, dena, denb, ntqg, ntdg, &
                       ntg, noa, nob, iswap1, iswap2)

      call uscexx_pots(dena, denb, vj, vxa, vxb)

      twoel_energy = 0.5d0*(spur(dena, vj)+spur(denb, vj) &
                   +spur(dena, vxa)+spur(denb, vxb))
      h0_energy = spur(dena, h0)+spur(denb, h0)
      total_energy = h0_energy+twoel_energy+ekern
      coulomb_energy = 0.5d0*(spur(dena, vj)+spur(denb, vj))
      exchange_energy = 0.5d0*(spur(dena, vxa)+spur(denb, vxb))
  
      ! print output
      if (iter .eq. 0) then
        diff_energy = total_energy
        write (iout, '(1x,i5,2x,f18.10,38x,f10.2)') &
          iter, total_energy, wallcl()-ts
      else
        diff_energy = total_energy-total_energy_old
        write (iout, '(1x,i5,2x,f18.10,2x,2e12.3,2x,i5,5x,f10.2)') &
          iter, total_energy, diff_energy, diff_density, &
          ndim-1, wallcl()-ts
      end if

      call uscexx_swap(orba, orbb, eiga, eigb, dena, denb, ntqg, ntdg, &
                       ntg, noa, nob, iswap1, iswap2)

      call uscexx_pots(dena, denb, vj, vxa, vxb)

    end if

    ! Create orbitals and density matrix with fractional
    ! occupation numbers, if required
    if (l_frac_occ) then

!      call get_orbb_frac(orbb, occb, ntg, ntqg, system, orbb_frac)
!      call den_from_orb(orbb_frac, noa, denb_frac)
!      call uscexx_pots(dena, denb_frac, vj_out=vj, vxb_out=vxb)

      system_aux = trim(system)//"1"
      occb1 = occa
      call get_orbb_frac(orba, occb1, ntg, ntqg, system_aux, orbb_frac)
      call den_from_orb(orbb_frac, 9, denb_frac)
      call uscexx_pots(denb_frac, denb, vxa_out=vxb1)

      system_aux = trim(system)//"2" 
      occb2 = occa
      call get_orbb_frac(orba, occb2, ntg, ntqg, system_aux, orbb_frac)
      call den_from_orb(orbb_frac, 9, denb_frac)
      call uscexx_pots(denb_frac, denb, vxa_out=vxb2)

      call get_orbb_frac(orba, occa, ntg, ntqg, system, orbb_frac)
      call den_from_orb(orbb_frac, 9, denb_frac)
      call uscexx_pots(denb_frac, denb, vj_out=vj, vxa_out=vxa)
    end if

    ! generate local exchange potentials
    if (l_spin_sym) then
      call uscexx_opm_spin_symmetrized(iter)
    else
      call uscexx_opm(iter)
    end if

    ! print eigenvalues
    if (iverbose > 0) then
      call print_eig(eiga, min(ntg, max(noa, nob)+10), &
                     'USCEXX eigenvalues (alpha)')
      write (iout, '(1x,a,1x,e12.5)') &
        'LUMO-HOMO (alpha):', (eiga(noa+1)-eiga(noa))
      call print_eig(eigb, min(ntg, max(noa, nob)+10), &
                     'USCEXX eigenvalues (beta)')
      write (iout, '(1x,a,1x,e12.5)') &
        'LUMO-HOMO (beta):', (eigb(nob+1)-eigb(nob))
      write (iout, *) ''
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
    focka_tri => memory_allocate(ntdg)
    focka_tri = 0d0
    fockb_tri => memory_allocate(ntdg)
    fockb_tri = 0d0
    focka_tri = h0+vj+vxa
    fockb_tri = h0+vj+vxb

    if (l_fixmix) then
      ! fixed mixing of fock matrices instead of DIIS
      if (iter > 0) then
        focka_tri = (1d0-fixmix_rate)*focka_tri+old_focka*fixmix_rate
        fockb_tri = (1d0-fixmix_rate)*fockb_tri+old_fockb*fixmix_rate
      else
        focka_tri = h0+vj+vxa
        fockb_tri = h0+vj+vxb
      end if
      old_focka = focka_tri
      old_fockb = fockb_tri
    else

      ! DIIS for accelerating convergence
      call udiis(dena, denb, focka_tri, fockb_tri, smh, sh_a, ndim, ndel, bfak)
      if (ndim .gt. maxdim_diis) then
        ndim = 1
        ndel = 0
        bfak = 1d0
      end if

    end if

    ! calculate orbitals and eigenvalues, expand fock matrix
    r1_tmp1 => memory_allocate(ntqg)
    focka = 0d0
    r1_tmp1 = 0d0
    call expan(focka_tri, focka, 1, 1, 1)
    call mmult(smh, focka, r1_tmp1, 1, 1, 0, -1, 0, 0)
    call mmult(r1_tmp1, smh, focka, 1, 1, 0, 0, 0, 0)
    call diag(focka, orba, eiga, 0)
    call mmult(smh, orba, r1_tmp1, 1, 1, 0, -1, 0, 0)
    orba = r1_tmp1

    fockb = 0d0
    r1_tmp1 = 0d0
    call expan(fockb_tri, fockb, 1, 1, 1)
    call mmult(smh, fockb, r1_tmp1, 1, 1, 0, -1, 0, 0)
    call mmult(r1_tmp1, smh, fockb, 1, 1, 0, 0, 0, 0)
    call diag(fockb, orbb, eigb, 0)
    call mmult(smh, orbb, r1_tmp1, 1, 1, 0, -1, 0, 0)
    orbb = r1_tmp1

    ! calculate densities
    r1_tmp1 = 0d0
    call dgemm_x('N', 'T', ntg, ntg, noa, 1d0, orba, ntg, &
                 orba, ntg, 0d0, r1_tmp1, ntg)
    call reduc(r1_tmp1, dena, 1, 1, 0.5d0, 1)

    r1_tmp1 = 0d0
    call dgemm_x('N', 'T', ntg, ntg, nob, 1d0, orbb, ntg, &
                 orbb, ntg, 0d0, r1_tmp1, ntg)
    call reduc(r1_tmp1, denb, 1, 1, 0.5d0, 1)

    call memory_release(focka_tri)

    ! calculate density difference for convergence criteria
    r1_tmp1 => memory_allocate(ntdg)
    r1_tmp1 = 0d0
    r1_tmp1 = dena+denb-denc_old
    diff_density = sqrt(ddot_x(ntdg, r1_tmp1, 1, r1_tmp1, 1)/dble(ntdg))
    denc_old = dena+denb

    call memory_release(ibase2)

  end do scf_loop

  if (iter .lt. maxit) then
    write (iout, '(/,1x,a,i5,a)') 'SCF converged in ', iter, ' iterations'
  else
    write (iout, '(/,1x,a,i5,a)') 'SCF NOT converged in ', iter, ' iterations'
  end if

  ! test potentials, if required
  if (l_test_pot) then
    if (ldfit) then
      write (iout, '(/,1x, a)') &
        'WARNING: test for potential is skipped. It works only with ldfit=0.'
    else
      call uscexx_pot_test(6, 1d-0)
    end if
  end if

  ! plot/write final potentials, if required (cases inside)
  if (l_spin_sym) then
    call plot_potential(orba, dena, vsol, naux_oep, 'final', 'OEP', noa, &
                        'vxa', vref_aux_up)
    call plot_potential(orbb, denb, vsol, naux_oep, 'final', 'OEP', nob, &
                        'vxb', vref_aux_up)
  else
    call plot_potential(orba, dena, vsola, naux_oep, 'final', 'OEP', noa, &
                        'vxa', vref_aux_upa)
    call plot_potential(orbb, denb, vsolb, naux_oep, 'final', 'OEP', nob, &
                        'vxb', vref_aux_upb)
  end if

  ! create USCEXX records
  irec_orb = irec_orb_out
  ifil_orb = ifil_orb_out

  call reserve_dump(irec_orb, ifil_orb, 'USCEXX', ntqg*4+ntg*4+ntdg*2)
  call write_orb(orba, 1, 'ALPHA')
  call write_orb(orbb, 2, 'BETA')
  call write_eig(eiga, 1, 'ALPHA')
  call write_eig(eigb, 2, 'BETA')
  if (l_frac_occ) then
    call write_den(dena+denb_frac, 1, 'CHARGE')
    call write_den(dena-denb_frac, 2, 'SPIN')
  else
    call write_den(dena+denb, 1, 'CHARGE')
    call write_den(dena-denb, 2, 'SPIN')
  end if
  call write_occ(occa, 1, 'ALPHA')
  call write_occ(occb, 2, 'BETA')
  call write_fock(focka+fockb, 1, 'TOTAL')
  call write_fock(focka-fockb, 2, 'OPEN')
  call flush_dump
  write (iout, '(/,1x,a,i8,a,i1)') &
    'USCEXX orbitals written to record', irec_orb, '.', ifil_orb

  ! print energy
  write (iout, *) ''
  call output_result('USCEXX', 'Total energy', &
                     total_energy, showstate=.false., principal=.true.)
  call output_result('USCEXX', 'One-electron energy', &
                     h0_energy, showstate=.false., principal=.false.)
  call output_result('USCEXX', 'Two-electron energy', &
                     twoel_energy, showstate=.false., principal=.false.)
  call output_result('USCEXX', 'Nuclear energy', &
                     ekern, showstate=.false., principal=.false.)
  call output_result('USCEXX', 'Coulomb energy', coulomb_energy, &
                     showstate=.false., principal=.false.)
  call output_result('USCEXX', 'Exchange energy', exchange_energy, &
                     showstate=.false., principal=.false.)

  ! calculate and print dipole moment
  dipole => memory_allocate(3)
  call dipole_moment(dipole, dena+denb)
  write (iout, *) ''
  call output_result('USCEXX', 'Dipole moment', dipole(1:3), &
       & showstate=.false., numberformat='3f15.8', debye=.TRUE.)
  call memory_release(dipole)

  ! calculate and print S**2
  r1_tmp1 => memory_allocate(ntqg)
  call les(r1_tmp1, ntdg, 1, sm)
  call expan(r1_tmp1, r1_tmp1, 1, 1, 1)
  call ssquared(dena+denb, dena-denb, r1_tmp1)
  call memory_release(r1_tmp1)

  ! print eigenvalues and LUMO-HOMO
  call print_eig(eiga, min(ntg, max(noa, nob)+10), &
                 'USCEXX eigenvalues (alpha)')
  write (iout, '(1x,a,1x,e12.5)') 'HOMO (alpha):', eiga(noa)
  write (iout, '(1x,a,1x,e12.5)') 'LUMO (alpha):', eiga(noa+1)
  write (iout, '(1x,a,1x,e12.5)') &
    'LUMO-HOMO (alpha):', (eiga(noa+1)-eiga(noa))
  call print_eig(eigb, min(ntg, max(noa, nob)+10), &
                 'USCEXX eigenvalues (beta)')
  write (iout, '(1x,a,1x,e12.5)') 'HOMO (beta):', eigb(nob)
  write (iout, '(1x,a,1x,e12.5)') 'LUMO (beta):', eigb(nob+1)
  write (iout, '(1x,a,1x,e12.5)') &
    'LUMO-HOMO (beta):', (eigb(nob+1)-eigb(nob))

  call print_eig(occa, min(ntg, max(noa, nob)+10), &
                 'USCEXX occupation numbers (alpha)')
  call print_eig(occb, min(ntg, max(noa, nob)+10), &
                 'USCEXX occupation numbers (beta)')


  if (l_swap2) then
    call uscexx_swap(orba, orbb, eiga, eigb, dena, denb, ntqg, ntdg, &
                     ntg, noa, nob, iswap1, iswap2)

    call uscexx_pots(dena, denb, vj, vxa, vxb)

    twoel_energy = 0.5d0*(spur(dena, vj)+spur(denb, vj) &
                 +spur(dena, vxa)+spur(denb, vxb))
    h0_energy = spur(dena, h0)+spur(denb, h0)
    total_energy = h0_energy+twoel_energy+ekern
    coulomb_energy = 0.5d0*(spur(dena, vj)+spur(denb, vj))
    exchange_energy = 0.5d0*(spur(dena, vxa)+spur(denb, vxb))

    ! print energy
    write (iout, *) ''
    call output_result('USCEXX', 'Total energy', &
                       total_energy, showstate=.false., principal=.true.)
    call output_result('USCEXX', 'One-electron energy', &
                       h0_energy, showstate=.false., principal=.false.)
    call output_result('USCEXX', 'Two-electron energy', &
                       twoel_energy, showstate=.false., principal=.false.)
    call output_result('USCEXX', 'Nuclear energy', &
                       ekern, showstate=.false., principal=.false.)
    call output_result('USCEXX', 'Coulomb energy', coulomb_energy, &
                       showstate=.false., principal=.false.)
    call output_result('USCEXX', 'Exchange energy', exchange_energy, &
                       showstate=.false., principal=.false.)
  end if

  ! set energy variables
  call setvar('ENERGY', total_energy, 'AU', 1, 1, imxs, 0)
  call output_energy(total_energy, 'USCEXX')

  ! release memory
  call memory_release(ibase)
  call deallocate_arrays()

end subroutine acfd_uscexx

subroutine uscexx_swap(orba, orbb, eiga, eigb, dena, denb, ntqg, ntdg, &
                       ntg, noa, nob, iswap1, iswap2)
  use memory, only: memory_allocate, memory_release
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: ntg, ntdg, ntqg, iswap1, iswap2, noa, nob
  real(dp), intent(inout) :: orba(ntqg), orbb(ntqg), eiga(ntg), eigb(ntg)
  real(dp), intent(out) :: dena(ntdg), denb(ntdg)

  real(dp), pointer :: r1_tmp1(:), r2_tmp1(:, :), r2_tmp2(:, :)

  r2_tmp1 => memory_allocate(ntg, ntg)
  r2_tmp2 => memory_allocate(ntg, ntg)

  r2_tmp1 = reshape(orba, (/ntg, ntg/))
  r2_tmp2 = reshape(orba, (/ntg, ntg/))

  r2_tmp2(:, iswap1) = r2_tmp1(:, iswap2)
  r2_tmp2(:, iswap2) = r2_tmp1(:, iswap1)

  orba = reshape(r2_tmp2, (/ntqg/))

  r2_tmp1 = reshape(orbb, (/ntg, ntg/))
  r2_tmp2 = reshape(orbb, (/ntg, ntg/))

  r2_tmp2(:, iswap1) = r2_tmp1(:, iswap2)
  r2_tmp2(:, iswap2) = r2_tmp1(:, iswap1)

  orbb = reshape(r2_tmp2, (/ntqg/))

  eiga(iswap1) = eiga(iswap2)
  eiga(iswap2) = eigb(iswap1)
  eigb = eiga

  call memory_release(r2_tmp1)

  ! Calculate densities from orbitals
  r1_tmp1 => memory_allocate(ntqg)
  r1_tmp1 = 0d0 
  call dgemm_x('N', 'T', ntg, ntg, noa, 1d0, &
               orba, ntg, &
               orba, ntg, &
               0d0, r1_tmp1, ntg)
  call reduc(r1_tmp1, dena, 1, 1, 1d0, 0)
  r1_tmp1 = 0d0 
  call dgemm_x('N', 'T', ntg, ntg, nob, 1d0, &
               orbb, ntg, &
               orbb, ntg, &
               0d0, r1_tmp1, ntg)
  call reduc(r1_tmp1, denb, 1, 1, 1d0, 0)
  call memory_release(r1_tmp1)

end subroutine uscexx_swap

subroutine get_orbb_frac(orbb, occb, ntg, ntqg, system, orbb_frac)
  use memory, only: memory_allocate, memory_release
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: ntg, ntqg
  real(dp), intent(in) :: orbb(ntqg)
  real(dp), intent(inout) :: occb(ntg)
  character(32), intent(in) :: system
  real(dp), intent(out) :: orbb_frac(ntqg)

  real(dp), pointer :: r2_tmp1(:, :)

  r2_tmp1 => memory_allocate(ntg, ntg)
  r2_tmp1 = reshape(orbb, (/ntg, ntg/))

  if (trim(system) == "B") then
    occb(3) = 0.5d0
    occb(4) = 0.5d0
    occb(5) = 0d0
    r2_tmp1(:, 3) = sqrt(occb(3))*r2_tmp1(:, 3)
    r2_tmp1(:, 4) = sqrt(occb(4))*r2_tmp1(:, 4)
    r2_tmp1(:, 5) = sqrt(occb(5))*r2_tmp1(:, 5)
  end if

  if (trim(system) == "B1") then
    occb(3) = 0d0
    occb(4) = 1d0
    occb(5) = 0d0
    r2_tmp1(:, 3) = sqrt(occb(3))*r2_tmp1(:, 3)
    r2_tmp1(:, 4) = sqrt(occb(4))*r2_tmp1(:, 4)
    r2_tmp1(:, 5) = sqrt(occb(5))*r2_tmp1(:, 5)
  end if

  if (trim(system) == "B2") then
    occb(3) = 1d0
    occb(4) = 0d0
    occb(5) = 0d0
    r2_tmp1(:, 3) = sqrt(occb(3))*r2_tmp1(:, 3)
    r2_tmp1(:, 4) = sqrt(occb(4))*r2_tmp1(:, 4)
    r2_tmp1(:, 5) = sqrt(occb(5))*r2_tmp1(:, 5)
  end if

  if (trim(system) == "C") then
    occb(3) = 1d0
    occb(4) = 0.5d0
    occb(5) = 0.5d0 
    r2_tmp1(:, 3) = sqrt(occb(3))*r2_tmp1(:, 3)
    r2_tmp1(:, 4) = sqrt(occb(4))*r2_tmp1(:, 4)
    r2_tmp1(:, 5) = sqrt(occb(5))*r2_tmp1(:, 5)
  end if

  if (trim(system) == "C1") then
    occb(3) = 1d0 
    occb(4) = 1d0 
    occb(5) = 0d0 
    r2_tmp1(:, 3) = sqrt(occb(3))*r2_tmp1(:, 3)
    r2_tmp1(:, 4) = sqrt(occb(4))*r2_tmp1(:, 4)
    r2_tmp1(:, 5) = sqrt(occb(5))*r2_tmp1(:, 5)
  end if

  if (trim(system) == "C2") then
    occb(3) = 1d0 
    occb(4) = 0d0 
    occb(5) = 1d0 
    r2_tmp1(:, 3) = sqrt(occb(3))*r2_tmp1(:, 3)
    r2_tmp1(:, 4) = sqrt(occb(4))*r2_tmp1(:, 4)
    r2_tmp1(:, 5) = sqrt(occb(5))*r2_tmp1(:, 5)
  end if

  if (trim(system) == "Al") then
    occb(7) = 1d0/2d0
    occb(8) = 1d0/2d0
    occb(9) = 0d0
    r2_tmp1(:, 7) = sqrt(occb(7))*r2_tmp1(:, 7)
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "Al1") then
    occb(7) = 0d0 
    occb(8) = 1d0
    occb(9) = 0d0
    r2_tmp1(:, 7) = sqrt(occb(7))*r2_tmp1(:, 7)
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "Al2") then
    occb(7) = 1d0 
    occb(8) = 0d0
    occb(9) = 0d0
    r2_tmp1(:, 7) = sqrt(occb(7))*r2_tmp1(:, 7)
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "Si") then
    occb(7) = 1.0d0
    occb(8) = 0.5d0
    occb(9) = 0.5d0
    r2_tmp1(:, 7) = sqrt(occb(7))*r2_tmp1(:, 7)
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "Si1") then
    occb(7) = 1d0
    occb(8) = 0d0
    occb(9) = 1d0 
    r2_tmp1(:, 7) = sqrt(occb(7))*r2_tmp1(:, 7)
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "Si2") then
    occb(7) = 1d0
    occb(8) = 1d0
    occb(9) = 0d0 
    r2_tmp1(:, 7) = sqrt(occb(7))*r2_tmp1(:, 7)
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "Cl_eq") then
    occb(7) = 2d0/3d0
    occb(8) = 2d0/3d0
    occb(9) = 2d0/3d0
    r2_tmp1(:, 7) = sqrt(occb(7))*r2_tmp1(:, 7)
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "HS") then
    occb(8) = 0.5d0
    occb(9) = 0.5d0
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "HS1") then
    occb(8) = 1d0
    occb(9) = 0d0
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  if (trim(system) == "HS2") then
    occb(8) = 0d0
    occb(9) = 1d0
    r2_tmp1(:, 8) = sqrt(occb(8))*r2_tmp1(:, 8)
    r2_tmp1(:, 9) = sqrt(occb(9))*r2_tmp1(:, 9)
  end if

  orbb_frac = reshape(r2_tmp1, (/ntqg/))

  call memory_release(r2_tmp1)

end subroutine get_orbb_frac
