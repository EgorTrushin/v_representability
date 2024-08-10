!
! Driver for spin-restricted (possibly open-shell) KS inversion calculations
!
! J. Chem. Phys. 156, 204124 (2022); https://doi.org/10.1063/5.0087356
!
! Last refactoring: E. Trushin <e.v.trushin@gmail.com>
!
! The code is formatted with fprettify:
! fprettify -i 2 -l 80 -w 1 -s
!
subroutine ksinv_driver(cmd)
  use common_cbas, only: ntg, ntdg, ntqg
  use common_molen, only: ekern
  use common_tapes, only: iout
  use memory, only: memory_allocate, memory_release, memory_allocate_integer
  use outputResult, only: output_result
  use acfd_derived_datatypes, only: vector_type
  use acfd_plot, only: vref_aux_up
  use acfd_utils, only: print_eig
  use acfd_wf, only: occa, occb
  use ksinv_basis_process, only: basis_process
  use ksinv_integrals, only: build_fkl
  use ksinv_mod, only: irec_orb, ifil_orb, irec_orb_out, ifil_orb_out, &
                       maxit, minit, iverbose, den, orb, eig, occ, no, smh, &
                       vx, vxc, vxc_sol, vx_sol, vj, h0, naux_oep, &
                       inversion_rhs, c_vsol_Delta, ref_den, vj_frozen, &
                       ref_v_ref_den, int_fkl, v_kin, reference_energy, &
                       thr_scf_inversion, l_plot_always, l_stopped, &
                       get_input_parameters, ksinv_init, deallocate_arrays, &
                       no_a, no_b, den_a, den_b, vx_a, vx_b, &
                       l_swap, l_swap2, iswap1, iswap2, l_frac_occ, system, & 
                       orbb_frac, denb_frac, vxb1, vxb2, occb1, occb2
  use ksinv_plot, only: plot_driver
  use ksinv_utils, only: get_u_frozen, get_orb_and_den, linear_response, &
                         get_energies, den_from_orb, get_potentials
  use iso_fortran_env, only: dp => real64
  implicit none
  ! identifies the algorithm by its internal name (here KSINV)
  character(64), intent(in) :: cmd

  ! auxiliary variables
  real(dp), pointer :: r1_tmp1(:)
  real(dp), pointer :: dipole(:)
  integer :: imxs

  integer, pointer :: ibase(:), ibase2(:) ! pointers for stack management
  integer :: iter ! number of current iteration
  character(5) :: tag ! for differentiation of potentials and densities (iter, final)

  real(dp) :: twoel_energy ! two-electron energy
  real(dp) :: h0_energy ! single electron energy
  real(dp) :: total_energy ! total energy
  real(dp) :: E0 ! total energy without correlation energy
  real(dp) :: total_energy_old ! total energy of previous iteration

  ! other energies
  real(dp) :: E_X_Ref, E_H_Ref, E_KIN_Ref, E_ext_Ref, &
              E_X, E_H, E_KIN, E_ext, E_ee_Ref, &
              T_c, V_c, E_c

  real(dp) :: ts ! variables for timings

  logical :: l_converged

  type(vector_type) :: vref_ao

  real(dp), external :: spur ! trace of the matrix product of two triagonal matrices
  real(dp), external :: ddot_x ! scalar product
  real(dp), external ::  wallcl ! timing function
  real(dp), external :: trace

  ! create a pointer to the bottom of the stack used in this driver
  ibase => memory_allocate_integer(1)

  write (iout, *)
  write (iout, *) 'PROGRAM * KSINV     (Kohn-Sham Inversion)'
  write (iout, *) 'Authors: J. Erhard, E. Trushin'
  write (iout, *)

  ! parse and print input parameters
  call get_input_parameters()

  ! zget orbitals, eigenvalues, density, etc
  call ksinv_init()

  if (l_frac_occ) then
    orbb_frac => memory_allocate(ntqg)
    denb_frac => memory_allocate(ntdg)
    vxb1 => memory_allocate(ntdg)
    vxb2 => memory_allocate(ntdg)
    occb1 => memory_allocate(ntg)
    occb2 => memory_allocate(ntg)
  end if

  occa => memory_allocate(ntg)
  occb => memory_allocate(ntg)
  occa = 0d0
  occb = 0d0
  occa(1:no_a) = 1d0
  occb(1:no_b) = 1d0

  ts = wallcl()

  ! print header
  write (iout, '(/,3x,a,8x,a,13x,a,10x,a)') 'ITER', 'ENERGY', '|t|', 'TIME'

  ! start scf cycle
  scf_loop: do iter = 0, maxit
    if (iter == 0) then

      ! Use the reference densities to get vj and vx^{nl}
      ! Calculate the correspontding energies
      r1_tmp1 => memory_allocate(1)
      call hfma(ref_den, ref_den, vj_frozen, vx, 1, 1, 0, &
                .true., .true., r1_tmp1, .false., .false.)
      call memory_release(r1_tmp1)

      ! calculate energy components
      call get_energies(ref_den, h0, v_kin, vx, vj_frozen, reference_energy, &
                        ekern, h0_energy, E_kin_Ref, E_x_Ref, E_H_Ref, &
                        E_ee_Ref, E_ext_Ref, ntdg)

      ! calculate three-index integrals
      call build_fkl('OEP', 'J', naux_oep, int_fkl)

      ! u is generated, needed for the construction of the reference potential
      call get_u_frozen()

      ! AO-Reference potential from FCI density
      call basis_process(naux_oep, vref_aux_up, vref_ao)
      ref_v_ref_den = vref_ao%mat

      ! diagonalise fock matrix to initial orbitals, eigenvalues and density
      if (no_a > 0) then
        call get_orb_and_den(h0+vj_frozen+ref_v_ref_den, &
                             orb, den, smh, eig, no_a, cmd, no_b)
      else
        call get_orb_and_den(h0+vj_frozen+ref_v_ref_den, &
                             orb, den, smh, eig, no, cmd)
      end if

      ! Sum up all contributionj into total energy of reference-density-matrix
      total_energy = E_X_Ref+E_H_Ref+E_KIN_Ref+E_ext_Ref+ekern

      ! print output
      write (iout, '(1x,i5,36x,f10.2)') iter, wallcl()-ts

      ! print additional info, if required
      if (iverbose > 0) then
        if (no_a > 0) then
          write (iout, '(/,1x,a,1x,f18.10)') 'HOMO (alpha):', eig(no_a)
          write (iout, '(1x,a,1x,f18.10)') 'LUMO (alpha):', eig(no_a+1)
          write (iout, '(/,1x,a,1x,f18.10)') 'HOMO (beta):', eig(no_b)
          write (iout, '(1x,a,1x,f18.10)') 'LUMO (beta):', eig(no_b+1)
        else
          write (iout, '(/,1x,a10,1x,f18.10)') 'HOMO:', eig(no)
          write (iout, '(1x,a10,1x,f18.10)') 'LUMO:', eig(no+1)
        end if
        write (iout, '(1x,a10,1x,f18.10)') 'E_X_Ref:', E_X_Ref
        write (iout, '(1x,a10,1x,f18.10)') 'E_H_Ref:', E_H_Ref
        write (iout, '(1x,a10,1x,f18.10)') 'E_KIN_Ref:', E_KIN_Ref
        write (iout, '(1x,a10,1x,f18.10)') 'E_ext_Ref:', E_ext_Ref
        write (iout, '(1x,a10,1x,f18.10)') 'E_ee_Ref:', E_ee_Ref
      end if

    else

      ! Create pointer to start of scf loop
      ibase2 => memory_allocate_integer(1)

      ! ###### Swap orbitals and eigenvalues, if required
      if (l_swap) then
        call ksinv_swap(orb, den, eig, ntqg, ntdg, ntg, no_a, no_b, &
                        iswap1, iswap2)
      end if

      ! calculate nonlocal exchange and coulomb matrix
      vj = 0d0
      vx = 0d0
      r1_tmp1 => memory_allocate(1)
      call hfma(den, den, vj, vx, 1, 1, 0, .true., .true., &
                r1_tmp1, .false., .false.)
      call memory_release(r1_tmp1)

      ! calculate two-electron energy
      call get_energies(den, h0, v_kin, vx, vj, reference_energy, ekern, &
                        h0_energy, E_kin, E_x, E_H, twoel_energy, E_ext, ntdg)

      twoel_energy = E_x+E_H
      E0 = h0_energy+twoel_energy+ekern

      T_c = E_KIN_Ref-E_KIN
      V_c = E_ee_Ref-E_H-E_X
      E_c = T_c+V_c

      if (l_frac_occ) then
        call get_orbb_frac(orb, occb, ntg, ntqg, system, orbb_frac)
        call den_from_orb(orbb_frac, no_a, denb_frac)
        den_a => memory_allocate(ntdg)
        call den_from_orb(orb, no_a, den_a) ! get alpha density
        den = den_a + denb_frac
        call memory_release(den_a)
!        call uscexx_pots(dena, denb_frac, vj_out=vj, vxb_out=vxb)
      end if

      ! generate local exchange potential, correlation potential
      if (no_a > 0) then
        call linear_response(eig, den, ref_den, orb, vxc, vxc_sol, vx_sol, &
                             c_vsol_Delta, inversion_rhs, no_a, iter, cmd, no_b)
      else
        call linear_response(eig, den, ref_den, orb, vxc, vxc_sol, vx_sol, &
                             c_vsol_Delta, inversion_rhs, no, iter, cmd)
      end if

      ! plot potentials and densities
      write (tag, '(I5.5)') iter
      l_converged = (inversion_rhs .lt. thr_scf_inversion)
      l_stopped = (iter .eq. maxit)
      if (l_converged .or. l_stopped) then
        tag = 'final'
      end if
      if (l_stopped .or. l_converged .or. l_plot_always) then
        call plot_driver(naux_oep, tag, 'total', vxc_sol, vx_sol, &
                         ref_den, den, vref_aux_up%mat)
      end if

      ! calculate total energy
      total_energy = E0+E_c

      ! print output
      write (iout, '(1x,i5,2x,f18.10,2x,e12.3,2x,f10.2)') &
        iter, total_energy, inversion_rhs, wallcl()-ts

      ! print additional info, if required
      if (iverbose > 0) then
        write (iout, '(/,1x,a6,1x,f18.10)') 'E_X:', E_X
        write (iout, '(1x,a6,1x,f18.10)') 'E_H:', E_H
        write (iout, '(1x,a6,1x,f18.10)') 'E_KIN:', E_KIN
        write (iout, '(1x,a6,1x,f18.10)') 'E_ext:', E_ext
        write (iout, '(1x,a6,1x,f18.10)') 'T_c:', T_c
        write (iout, '(1x,a6,1x,f18.10)') 'V_c:', V_c
        write (iout, '(1x,a6,1x,f18.10)') 'E_c:', E_c
      end if

      total_energy_old = total_energy

      ! check for convergence
      if (iter .ge. minit) then
        if ((thr_scf_inversion .le. 0d0) .or. &
            abs(inversion_rhs) .lt. thr_scf_inversion) then
          write (iout, '(/,1x,a,3x,f18.12)') 'Diff E_x:', E_X_Ref-E_X
          write (iout, '(1x,a,3x,f18.12)') 'Diff E_H:', E_H_Ref-E_H
          write (iout, '(1x,a,1x,f18.12)') 'Diff E_kin:', E_KIN_Ref-E_KIN
          write (iout, '(1x,a,1x,f18.12)') 'Diff E_ext:', E_ext_Ref-E_ext
          exit scf_loop
        end if
      end if

      if (iter .eq. maxit) exit scf_loop

      ! Takes potentials and overlap matrix (S^{-1/2}) and generates orbitals, eigenavlues and density
      if (no_a > 0) then
        call get_orb_and_den(h0+vj_frozen+vxc+ref_v_ref_den, orb, &
                             den, smh, eig, no_a, cmd, no_b)
      else
        call get_orb_and_den(h0+vj_frozen+vxc+ref_v_ref_den, orb, &
                             den, smh, eig, no, cmd)
      end if

      ! If required, print accuracy of Hartree energy and external energy in inversion, HOMO and LUMO
      if (iverbose .gt. 0) then
        if (no_a > 0) then
          write (iout, '(/,1x,a,1x,f18.10)') 'HOMO (alpha):', eig(no_a)
          write (iout, '(1x,a,1x,f18.10)') 'LUMO (alpha):', eig(no_a+1)
          write (iout, '(/,1x,a,1x,f18.10)') 'HOMO (beta):', eig(no_b)
          write (iout, '(1x,a,1x,f18.10)') 'LUMO (beta):', eig(no_b+1)
        else
          write (iout, '(/,1x,a10,1x,f18.10)') 'HOMO:', eig(no)
          write (iout, '(1x,a10,1x,f18.10)') 'LUMO:', eig(no+1)
        end if
        write (iout, '(1x,a,1x,es21.14)') 'Diff E_H:', E_H_Ref-E_H
        write (iout, '(1x,a,1x,es21.14)') 'Diff E_ext:', E_ext_Ref-E_ext
      end if

      call memory_release(ibase2)
    end if
  end do scf_loop

  if (iter .lt. maxit) then
    write (iout, '(/,1x,a,i5,a)') 'SCF converged in ', iter, ' iterations'
  else
    write (iout, '(/,1x,a,i5,a)') 'SCF NOT converged in ', iter, ' iterations'
  end if

!  if (l_swap2) then
!    call ksinv_swap(orb, den, eig, ntqg, ntdg, ntg, no_a, no_b, &
!                    iswap1, iswap2)
!  end if

  if (no_a > 0) then
    ! Calculate correct exchange energy and energies dependent on it
    den_a => memory_allocate(ntdg)
    den_b => memory_allocate(ntdg)
    vx_a => memory_allocate(ntdg)
    vx_b => memory_allocate(ntdg)
    call den_from_orb(orb, no_a, den_a) ! get alpha density
    call den_from_orb(orb, no_b, den_b) ! get beta density
    call get_potentials(den_a, den_b, vx_a, vx_b) ! get exchange potentials
!    call get_energies(den, h0, v_kin, vx, vj, reference_energy, ekern, &
!                      h0_energy, E_kin, E_x, E_H, twoel_energy, E_ext, ntdg)
    twoel_energy = 0.5d0*(spur(den_a, vx_a)+spur(den_b, vx_b)+spur(den, vj))
    E0 = h0_energy+twoel_energy+ekern
    E_X = 0.5*(spur(den_a, vx_a)+spur(den_b, vx_b))
    T_c = E_KIN_Ref-E_KIN
    V_c = E_ee_Ref-E_H-E_X
    E_c = T_c+V_c
  end if

  ! create KSINV records
  irec_orb = irec_orb_out
  ifil_orb = ifil_orb_out
  call reserve_dump(irec_orb, ifil_orb, 'KSINV', ntqg*2+ntg*2+ntdg)
  call write_orb(orb, 1, 'CANONICAL')
  call write_den(den, 1, 'CHARGE')
  call write_eig(eig, 1, 'CANONICAL')
  call write_occ(occ, 1, 'OCC')
  call flush_dump
  write (iout, '(/,1x,a,i8,a,i1)') &
    'KSINV orbitals written to record', irec_orb, '.', ifil_orb

  ! print energy
  write (iout, *) ''
  call output_result('KSINV', 'Total energy', &
                     total_energy, showstate=.false., principal=.true.)
  call output_result('KSINV', 'One-electron energy', &
                     h0_energy, showstate=.false., principal=.false.)
  call output_result('KSINV', 'Two-electron energy', &
                     E_X+E_H, showstate=.false., principal=.false.)
  call output_result('KSINV', 'Nuclear energy', &
                     ekern, showstate=.false., principal=.false.)
  call output_result('KSINV', 'Coulomb energy', &
                     E_H, showstate=.false., principal=.false.)
  call output_result('KSINV', 'Exchange energy', &
                     E_X, showstate=.false., principal=.false.)
  call output_result('KSINV', 'Correlation energy', &
                     E_c, showstate=.false., principal=.false.)

  ! calculate and print dipole moment
  dipole => memory_allocate(3)
  call dipole_moment(dipole, den)
  write (iout, *) ''
  call output_result('KSINV', 'Dipole moment', dipole(1:3), &
       & showstate=.false., numberformat='3f15.8', debye=.TRUE.)
  call memory_release(dipole)

  ! print eigenvalues and LUMO-HOMO
  call print_eig(eig, min(ntg, no+10), 'KSINV eigenvalues')

  if (no_a > 0) then
    write (iout, '(1x,a,1x,e12.5)') 'HOMO (alpha):', eig(no_a)
    write (iout, '(1x,a,1x,e12.5)') 'LUMO (alpha):', eig(no_a+1)
    write (iout, '(1x,a,1x,e12.5)') 'LUMO-HOMO (alpha):', &
      (eig(no_a+1)-eig(no_a))
    write (iout, '(1x,a,1x,e12.5)') 'HOMO (beta):', eig(no_b)
    write (iout, '(1x,a,1x,e12.5)') 'LUMO (beta):', eig(no_b+1)
    write (iout, '(1x,a,1x,e12.5)') 'LUMO-HOMO (beta):', (eig(no_b+1)-eig(no_b))
  else
    write (iout, '(1x,a,1x,e12.5)') 'HOMO:', eig(no)
    write (iout, '(1x,a,1x,e12.5)') 'LUMO:', eig(no+1)
    write (iout, '(1x,a,1x,e12.5)') 'LUMO-HOMO:', (eig(no+1)-eig(no))
  end if

  call print_eig(occa, min(ntg, max(no_a, no_b)+10), &
                 'Occupation numbers (alpha)')
  call print_eig(occb, min(ntg, max(no_a, no_b)+10), &
                 'Occupation numbers (beta)')

  ! set energy variables
  call setvar('ENERGY', total_energy, 'AU', 1, 1, imxs, 0)
  call output_energy(total_energy, 'KSINV')

  ! release memory
  call memory_release(ibase)
  call deallocate_arrays()

end subroutine ksinv_driver

subroutine ksinv_swap(orb, den, eig, ntqg, ntdg, ntg, noa, nob, iswap1, iswap2)
  use memory, only: memory_allocate, memory_release
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: ntg, ntdg, ntqg, iswap1, iswap2, noa, nob
  real(dp), intent(inout) :: orb(ntqg), eig(ntg)
  real(dp), intent(out) :: den(ntdg)

  real(dp), pointer :: r1_tmp1(:), r1_tmp2(:), r2_tmp1(:, :), r2_tmp2(:, :)
  real(dp) :: rtmp
  
  r2_tmp1 => memory_allocate(ntg, ntg)
  r2_tmp2 => memory_allocate(ntg, ntg)

  r2_tmp1 = reshape(orb, (/ntg, ntg/))
  r2_tmp2 = reshape(orb, (/ntg, ntg/))

  r2_tmp2(:, iswap1) = r2_tmp1(:, iswap2)
  r2_tmp2(:, iswap2) = r2_tmp1(:, iswap1)

  orb = reshape(r2_tmp2, (/ntqg/))

  rtmp = eig(iswap1)
  eig(iswap1) = eig(iswap2)
  eig(iswap2) = rtmp

  call memory_release(r2_tmp1)

  ! Calculate densities from orbitals
  r1_tmp1 => memory_allocate(ntqg)
  r1_tmp1 = 0d0
  call dgemm_x('N', 'T', ntg, ntg, noa, 1d0, &
               orb, ntg, &
               orb, ntg, &
               0d0, r1_tmp1, ntg)
  r1_tmp2 => memory_allocate(ntqg)
  r1_tmp2 = 0d0
  call dgemm_x('N', 'T', ntg, ntg, nob, 1d0, &
               orb, ntg, &
               orb, ntg, &
               0d0, r1_tmp2, ntg)
  r1_tmp1 = r1_tmp1 + r1_tmp2
  call reduc(r1_tmp1, den, 1, 1, 1d0, 0)
  call memory_release(r1_tmp1)

end subroutine ksinv_swap
