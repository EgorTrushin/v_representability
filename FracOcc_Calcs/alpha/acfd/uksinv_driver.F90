!
! Driver for spin-unrestricted KS inversion calculations
!
! J. Chem. Phys. 156, 204124 (2022); https://doi.org/10.1063/5.0087356
!
! Last refactoring: E. Trushin <e.v.trushin@gmail.com>
!
! The code is formatted with fprettify:
! fprettify -i 2 -l 80 -w 1 -s
!
subroutine uksinv_driver(cmd)
  use common_cbas, only: ntg, ntdg, ntqg
  use common_molen, only: ekern
  use common_tapes, only: iout
  use memory, only: memory_allocate, memory_release, memory_allocate_integer
  use outputResult, only: output_result
  use acfd_derived_datatypes, only: vector_type
  use acfd_plot, only: vref_aux_up
  use acfd_utils, only: print_eig
  use ksinv_basis_process, only: basis_process
  use ksinv_integrals, only: build_fkl
  use ksinv_mod, only: irec_orb, ifil_orb, irec_orb_out, ifil_orb_out, &
                       maxit, minit, iverbose, den, den_a, den_b, orb_a, &
                       orb_b, eig_a, eig_b, no_a, no_b, smh, vx, vx_a, vx_b, &
                       vxc_a, vxc_b, vxc_sol_a, vxc_sol_b, vx_sol_a, vx_sol_b, &
                       vj, h0, naux_oep, inversion_rhs_a, inversion_rhs_b, &
                       c_vsol_Delta_a, c_vsol_Delta_b, ref_den, &
                       ref_den_a, ref_den_b, vj_frozen, ref_v_ref_den, &
                       int_fkl, naux_oep, v_kin, reference_energy, occ, &
                       thr_scf_inversion, l_plot_always, l_stopped, &
                       get_input_parameters, uksinv_init, deallocate_arrays
  use ksinv_plot, only: plot_driver
  use ksinv_utils, only: get_u_frozen, get_orb_and_den, linear_response, &
                         get_energies
  use iso_fortran_env, only: dp => real64
  implicit none
  ! identifies the algorithm by its internal name (here UKSINV)
  character(64), intent(in) :: cmd

  ! auxiliary variables
  real(dp), pointer :: r1_tmp1(:)
  real(dp), pointer :: dipole(:)
  integer :: imxs

  integer, pointer :: ibase(:), ibase2(:) ! pointers for stack management
  integer :: iter ! number of current iteration
  character(5) :: tag ! for differentiation of potentials and densities (iter, final)

  real(dp) :: twoel_energy ! two-electron energy
  real(dp) :: h0_energy, h0_energy_a, h0_energy_b ! single electron energy
  real(dp) :: total_energy ! total energy
  real(dp) :: E0 ! total energy without correlation energy
  real(dp) :: total_energy_old ! total energy of previous iteration

  ! other energies
  real(dp) :: E_X_Ref, E_H_Ref, E_KIN_Ref, E_ext_Ref, &
              E_X, E_H, E_KIN, E_ext, E_ee_Ref, &
              T_c, V_c, E_c
  real(dp) :: E_X_Ref_a, E_H_Ref_a, E_KIN_Ref_a, E_ext_Ref_a
  real(dp) :: E_X_Ref_b, E_H_Ref_b, E_KIN_Ref_b, E_ext_Ref_b
  real(dp) :: E_X_a, E_H_a, E_KIN_a, E_ext_a
  real(dp) :: E_X_b, E_H_b, E_KIN_b, E_ext_b

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
  write (iout, *) 'PROGRAM * UKSINV     (Kohn-Sham Inversion)'
  write (iout, *) 'Authors: J. Erhard, E. Trushin'
  write (iout, *)

  ! parse and print input parameters
  call get_input_parameters()

  ! get orbitals, eigenvalues, density, etc
  call uksinv_init()

  ts = wallcl()

  ! print header
  write (iout, '(/,3x,a,8x,a,12x,a,9x,a,10x,a)') &
    'ITER', 'ENERGY', '|t_a|', '|t_b|', 'TIME'

  ! start scf cycle
  scf_loop: do iter = 0, maxit
    if (iter == 0) then

      ! Use the reference densities to get vj and vx^{nl}
      ! Calculate the correspontding energies
      r1_tmp1 => memory_allocate(1)
      call hfma(ref_den_a, ref_den_a, vj_frozen, vx_a, 1, 1, 0, &
                .true., .true., r1_tmp1, .false., .false.)
      call hfma(ref_den_b, ref_den_b, vj_frozen, vx_b, 1, 1, 0, &
                .true., .true., r1_tmp1, .false., .false.)
      call hfma(ref_den, ref_den, vj_frozen, vx, 1, 1, 0, &
                .true., .true., r1_tmp1, .false., .false.)
      vx_a = 2d0*vx_a
      vx_b = 2d0*vx_b
      call memory_release(r1_tmp1)

      ! calculate energy contributions from reference density matrix
      call get_energies(ref_den_a, h0, v_kin, vx_a, vj_frozen, &
                        reference_energy, ekern, h0_energy_a, E_kin_Ref_a, &
                        E_x_Ref_a, E_H_Ref_a, E_ee_Ref, E_ext_Ref_a, ntdg)
      call get_energies(ref_den_b, h0, v_kin, vx_b, vj_frozen, &
                        reference_energy, ekern, h0_energy_b, E_kin_Ref_b, &
                        E_x_Ref_b, E_H_Ref_b, E_ee_Ref, E_ext_Ref_b, ntdg)
      call get_energies(ref_den, h0, v_kin, vx, vj_frozen, reference_energy, &
                        ekern, h0_energy, E_kin_Ref, E_x_Ref, E_H_Ref, &
                        E_ee_Ref, E_ext_Ref, ntdg)
      E_x_Ref = E_x_Ref_a+E_x_Ref_b

      ! calculate three-index integrals
      call build_fkl('OEP', 'J', naux_oep, int_fkl)

      ! u is generated, needed for the construction of the reference potential
      call get_u_frozen()

      ! AO-Reference potential from FCI density
      call basis_process(naux_oep, vref_aux_up, vref_ao)
      ref_v_ref_den = vref_ao%mat

      ! diagonalise fock matrix to initial orbitals, eigenvalues and density
      call get_orb_and_den(h0+vj_frozen+ref_v_ref_den, orb_a, den_a, smh, &
                           eig_a, no_a, cmd)
      call get_orb_and_den(h0+vj_frozen+ref_v_ref_den, orb_b, den_b, smh, &
                           eig_b, no_b, cmd)
      den = den_a+den_b

      ! Sum up all contributionj into total energy of reference-density-matrix
      total_energy = E_X_Ref+E_H_Ref+E_KIN_Ref+E_ext_Ref+ekern

      ! print output
      write (iout, '(1x,i5,50x,f10.2)') iter, wallcl()-ts

      ! print additional info, if required
      if (iverbose > 0) then
        write (iout, '(/,1x,a10,19x,2(1x,f18.10))') &
          'HOMO:', eig_a(no_a), eig_b(no_b)
        write (iout, '(1x,a10,19x,2(1x,f18.10))') &
          'LUMO:', eig_a(no_a+1), eig_b(no_b+1)
        write (iout, '(1x,a10,3(1x,f18.10))') &
          'E_X_Ref:', E_X_Ref, E_X_Ref_a, E_X_Ref_b
        write (iout, '(1x,a10,3(1x,f18.10))') &
          'E_H_Ref:', E_H_Ref, E_H_Ref_a, E_H_Ref_b
        write (iout, '(1x,a10,3(1x,f18.10))') &
          'E_KIN_Ref:', E_KIN_Ref, E_KIN_Ref_a, E_KIN_Ref_b
        write (iout, '(1x,a10,3(1x,f18.10))') &
          'E_ext_Ref:', E_ext_Ref, E_ext_Ref_a, E_ext_Ref_b
        write (iout, '(1x,a10,1x,f18.10)') 'E_ee_Ref:', E_ee_Ref
      end if

    else

      ! Create pointer to start of scf loop
      ibase2 => memory_allocate_integer(1)

      ! calculate nonlocal exchange and coulomb matrix
      vj = 0d0
      vx = 0d0
      r1_tmp1 => memory_allocate(1)
      call hfma(den_a, den_a, vj, vx_a, -1, 1, 0, .true., .true., r1_tmp1, &
                .false., .false.)
      call hfma(den_a+den_b, den_a-den_b, vj, vx_b, -1, 1, 0, &
                .true., .true., r1_tmp1, .false., .false.)
      vx_b = vx_a-vx_b
      vj = (vj-(vx_a+vx_b))
      vx_a = 2d0*vx_a
      vx_b = 2d0*vx_b
      call memory_release(r1_tmp1)

      ! calculate two-electron energy
      call get_energies(den, h0, v_kin, vx, vj, reference_energy, ekern, &
                        h0_energy, E_kin, E_x, E_H, twoel_energy, E_ext, ntdg)
      call get_energies(den_a, h0, v_kin, vx_a, vj, reference_energy, &
                        ekern, h0_energy_a, E_kin_a, E_x_a, E_H_a, &
                        twoel_energy, E_ext_a, ntdg)
      call get_energies(den_b, h0, v_kin, vx_b, vj, reference_energy, &
                        ekern, h0_energy_b, E_kin_b, E_x_b, E_H_b, &
                        twoel_energy, E_ext_b, ntdg)

      E_x = E_x_a+E_x_b
      twoel_energy = E_x+E_H
      E0 = h0_energy+twoel_energy+ekern
      T_c = E_KIN_Ref-E_KIN
      V_c = E_ee_Ref-E_H-E_X
      E_c = T_c+V_c

      ! print additional info, if required
      if (iverbose > 0) then
        write (iout, '(/,1x,a6,3(1x,f18.10))') 'E_X:', E_X, E_x_a, E_x_b
        write (iout, '(1x,a6,3(1x,f18.10))') 'E_H:', E_H, E_H_a, E_H_b
        write (iout, '(1x,a6,3(1x,f18.10))') 'E_KIN:', E_KIN, E_kin_a, E_kin_b
        write (iout, '(1x,a6,3(1x,f18.10))') 'E_ext:', E_ext, E_ext_a, E_ext_b
        write (iout, '(1x,a6,1x,f18.10)') 'T_c:', T_c
        write (iout, '(1x,a6,1x,f18.10)') 'V_c:', V_c
        write (iout, '(1x,a6,1x,f18.10)') 'E_c:', E_c
      end if

      ! generate local exchange potential, correlation potential and correlation energy
      if (abs(inversion_rhs_a) .gt. thr_scf_inversion .or. iter .lt. 2) then
        vx = vx_a ! hook to make exx_oep working, required for plotting
        call linear_response(eig_a, den_a, ref_den_a, orb_a, vxc_a, vxc_sol_a, &
                             vx_sol_a, c_vsol_Delta_a, inversion_rhs_a, no_a, &
                             iter, cmd)
      end if
      if (abs(inversion_rhs_b) .gt. thr_scf_inversion .or. iter .lt. 2) then
        if (no_b /= 0) then
          vx = vx_b ! hook to make exx_oep working, required for plotting
          call linear_response(eig_b, den_b, ref_den_b, orb_b, vxc_b, &
                               vxc_sol_b, vx_sol_b, c_vsol_Delta_b, &
                               inversion_rhs_b, no_b, iter, cmd)
        end if
      end if

      ! plot potentials and densities
      write (tag, '(I5.5)') iter
      l_converged = (2d0*inversion_rhs_a .lt. thr_scf_inversion &
                     .and. 2d0*inversion_rhs_b .lt. thr_scf_inversion)
      l_stopped = (iter .eq. maxit)
      if (l_converged .or. l_stopped) then
        tag = 'final'
      end if
      if (l_stopped .or. l_converged .or. l_plot_always) then
        call plot_driver(naux_oep, tag, 'alpha', vxc_sol_a, vx_sol_a, &
                         ref_den_a, den_a, vref_aux_up%mat)
        if (no_b /= 0) then
          call plot_driver(naux_oep, tag, 'beta', vxc_sol_b, vx_sol_b, &
                           ref_den_b, den_b, vref_aux_up%mat)
        end if
      end if

      ! calculate total energy
      total_energy = E0+E_c

      ! print output
      write (iout, '(1x,i5,2x,f18.10,2(2x,e12.3),2x,f10.2)') &
        iter, total_energy, inversion_rhs_a, inversion_rhs_b, wallcl()-ts

      total_energy_old = total_energy

      ! check for convergence
      if (iter .ge. minit) then
        if ((thr_scf_inversion .le. 0d0) .or. &
            (abs(inversion_rhs_a) .lt. thr_scf_inversion &
             .and. abs(inversion_rhs_b) .lt. thr_scf_inversion)) then
          write (iout, '(/,1x,a,3x,2(f18.12))') &
            'Diff E_x:', E_X_Ref_a-E_X_a, E_X_Ref_b-E_X_b
          write (iout, '(1x,a,3x,2(f18.12))') &
            'Diff E_H:', E_H_Ref_a-E_H_a, E_H_Ref_b-E_H_b
          write (iout, '(1x,a,1x,2(f18.12))') &
            'Diff E_kin:', E_KIN_Ref_a-E_KIN_a, E_KIN_Ref_b-E_KIN_b
          write (iout, '(1x,a,1x,2(f18.12))') &
            'Diff E_ext:', E_ext_Ref_a-E_ext_a, E_ext_Ref_b-E_ext_b
          exit scf_loop
        end if
      end if

      if (iter .eq. maxit) exit scf_loop

      ! Takes potentials and overlap matrix (S^{-1/2}) and generates orbitals, eigenavlues and density
      if (abs(inversion_rhs_a) .gt. thr_scf_inversion) then
        call get_orb_and_den(h0+vj_frozen+vxc_a+ref_v_ref_den, orb_a, den_a, &
                             smh, eig_a, no_a, cmd)
      end if
      if (abs(inversion_rhs_b) .gt. thr_scf_inversion) then
        if (no_b /= 0) then
          call get_orb_and_den(h0+vj_frozen+vxc_b+ref_v_ref_den, orb_b, den_b, &
                               smh, eig_b, no_b, cmd)
        else
          den_b = 0d0
        end if
      end if

      den = den_a+den_b

      ! If required, print accuracy of Hartree energy and external energy in inversion, HOMO and LUMO
      if (iverbose .gt. 0) then
        write (iout, '(1x,a11,2x,a6,14x,2x,a5,15x)') '', 'Alpha:', 'Beta:'
        write (iout, '(1x,a11,2(1x,es21.14))') &
          'HOMO:', eig_a(no_a), eig_b(no_b)
        write (iout, '(1x,a11,2(1x,es21.14))') &
          'LUMO:', eig_a(no_a+1), eig_b(no_b+1)
        write (iout, '(1x,a11,2(1x,es21.14))') &
          'Diff E_H:', E_H_Ref_b-E_H_b, E_H_Ref_a-E_H_a
        write (iout, '(1x,a11,2(1x,es21.14))') &
          'Diff E_ext:', E_ext_Ref_a-E_ext_a, E_ext_Ref_b-E_ext_b
      end if

      call memory_release(ibase2)
    end if
  end do scf_loop

  if (iter .lt. maxit) then
    write (iout, '(/,1x,a,i5,a)') 'SCF converged in ', iter, ' iterations'
  else
    write (iout, '(/,1x,a,i5,a)') 'SCF NOT converged in ', iter, ' iterations'
  end if

  ! create UKSINV records
  irec_orb = irec_orb_out
  ifil_orb = ifil_orb_out
  call reserve_dump(irec_orb, ifil_orb, 'UKSINV', 2*(ntqg*2+ntg*2+ntdg))
  call write_orb(orb_a, 1, 'ALPHA')
  call write_orb(orb_b, 2, 'BETA')
  call write_eig(eig_a, 1, 'ALPHA')
  call write_eig(eig_b, 2, 'BETA')
  call write_den(den_a+den_b, 1, 'CHARGE')
  call write_den(den_a-den_b, 2, 'SPIN')
  occ = 0d0
  occ(1:no_a) = 1d0
  call write_occ(occ, 1, 'ALPHA')
  occ = 0d0
  occ(1:no_b) = 1d0
  call write_occ(occ, 2, 'BETA')
  call flush_dump
  write (iout, '(/,1x,a,i8,a,i1)') &
    'UKSINV orbitals written to record', irec_orb, '.', ifil_orb

  ! print energy
  write (iout, *) ''
  call output_result('UKSINV', 'Total energy', &
                     total_energy, showstate=.false., principal=.true.)
  call output_result('UKSINV', 'One-electron energy', &
                     h0_energy, showstate=.false., principal=.false.)
  call output_result('UKSINV', 'Two-electron energy', &
                     E_X+E_H, showstate=.false., principal=.false.)
  call output_result('UKSINV', 'Nuclear energy', &
                     ekern, showstate=.false., principal=.false.)
  call output_result('UKSINV', 'Coulomb energy', &
                     E_H, showstate=.false., principal=.false.)
  call output_result('UKSINV', 'Exchange energy', &
                     E_X, showstate=.false., principal=.false.)
  call output_result('UKSINV', 'Correlation energy', &
                     E_c, showstate=.false., principal=.false.)

  ! calculate and print dipole moment
  dipole => memory_allocate(3)
  call dipole_moment(dipole, den)
  write (iout, *) ''
  call output_result('UKSINV', 'Dipole moment', dipole(1:3), &
       & showstate=.false., numberformat='3f15.8', debye=.TRUE.)
  call memory_release(dipole)

  ! print eigenvalues and LUMO-HOMO
  call print_eig(eig_a, min(ntg, no_a+10), 'alpha - UKSINV eigenvalues')
  call print_eig(eig_b, min(ntg, no_b+10), 'beta - UKSINV eigenvalues')

  write (iout, '(1x,a10,2x,a6,5x,2x,a5,4x)') '', 'Alpha:', 'Beta:'
  write (iout, '(1x,a10,2(1x,e12.5))') &
    'HOMO:', eig_a(no_a), eig_b(no_b)
  write (iout, '(1x,a10,2(1x,e12.5))') &
    'LUMO:', eig_a(no_a+1), eig_b(no_b+1)
  write (iout, '(1x,a10,2(1x,e12.5))') &
    'LUMO-HOMO:', (eig_a(no_a+1)-eig_a(no_a)), (eig_b(no_b+1)-eig_b(no_b))

  ! set energy variables
  call setvar('ENERGY', total_energy, 'AU', 1, 1, imxs, 0)
  call output_energy(total_energy, 'UKSINV')

  ! release memory
  call memory_release(ibase)
  call deallocate_arrays()

end subroutine uksinv_driver
