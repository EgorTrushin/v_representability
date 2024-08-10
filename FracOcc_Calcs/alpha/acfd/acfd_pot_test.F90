module acfd_pot_test
  implicit none
  private
  public :: scexx_pot_test, uscexx_pot_test
contains
  subroutine scexx_pot_test(ntest, eps0)
    use common_tapes, only: iout
    use common_cbas, only: ntqg, ntdg, ntg
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_opm, only: scexx_opm
    use acfd_pot, only: fock, vx, vj
    use acfd_wf, only: orb, den, eig
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: ntest
    real(dp), intent(in) :: eps0

    real(dp), pointer :: orb_save(:), den_save(:), eig_save(:), fock_save(:), &
                         vx_save(:), vj_save(:)

    real(dp), pointer :: vxnl0(:), vxl0(:), vj0(:)
    real(dp), pointer :: r1_tmp1(:)
    integer, pointer :: ibase(:)

    integer :: itest
    real(dp) :: eps
    real(dp), pointer :: diff_Ex(:), diff_rhox(:), diff_diffx(:), diff_rhsx(:)
    real(dp), pointer :: vpertx(:)

    ibase => memory_allocate_integer(1)

    ! save properties
    orb_save => memory_allocate(ntqg)
    orb_save = orb
    den_save => memory_allocate(ntdg)
    den_save = den
    eig_save => memory_allocate(ntg)
    eig_save = eig
    fock_save => memory_allocate(ntqg)
    fock_save = fock
    vx_save => memory_allocate(ntdg)
    vx_save = vx
    vj_save => memory_allocate(ntdg)
    vj_save = vj

    ! calculate VxNL and Vj
    vj = 0d0
    vx = 0d0
    vxnl0 => memory_allocate(ntdg)
    vj0 => memory_allocate(ntdg)
    r1_tmp1 => memory_allocate(1)
    call hfma(den, den, vj, vx, 1, 1, 0, .true., .true., &
              r1_tmp1, .false., .false.)
    vxnl0 = vx
    vj0 = vj
    call memory_release(r1_tmp1)

    ! calculate VxL
    vxl0 => memory_allocate(ntdg)
    call scexx_opm(999)
    vxl0 = vx

    ! Start Test
    diff_Ex => memory_allocate(ntest)
    diff_rhox => memory_allocate(ntest)
    diff_diffx => memory_allocate(ntest)
    diff_rhsx => memory_allocate(ntest)
    vpertx => memory_allocate(ntdg)

    eps = eps0
    do itest = 1, ntest
      eps = eps*1d-1
      vx = vxnl0
      call scexx_opm_test(eps, vpertx, diff_rhsx(itest))
      call scexx_pert_exch(den, fock, vxnl0, vxl0, vpertx, &
                           diff_Ex(itest), diff_rhox(itest))
    end do
    diff_diffx = diff_Ex-diff_rhox

    ! load properties
    orb = orb_save
    den = den_save
    eig = eig_save
    fock = fock_save
    vx = vx_save
    vj = vj_save

    write (iout, '(/,1x,a,7x,a,9x,a,7x,a,5x,a,2x,a,3x,a)') &
      'Perturbation', 'd_Ex', 'd_rho_vx', 'd_rhs_x', &
      'd_Ex-d_rho_vx', 'd_Ex-d_rhs_x', 'd_rho_vx-d_rhs_x'
    eps = eps0
    do itest = 1, ntest
      eps = eps*1d-1
      write (iout, '(1x,7(es10.3,5x))') eps, diff_Ex(itest), diff_rhox(itest), &
        diff_rhsx(itest), diff_diffx(itest), &
        diff_Ex(itest)-diff_rhsx(itest), &
        diff_rhox(itest)-diff_rhsx(itest)
    end do

    call memory_release(ibase)
  end subroutine scexx_pot_test

  subroutine scexx_opm_test(eps, vpertx, drhsx)
    use common_cbas, only: ntqg, ntdg
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_auxiliary_basis, only: trans_mat_constraint, naux_oep, &
                                    naux_constraint
    use acfd_basis_process, only: basis_process_scexx
    use acfd_derived_datatypes, only: vector_type
    use acfd_plot, only: vref_aux_dummy
    use acfd_pot, only: vx
    use acfd_integrals, only: build_fai, build_fij
    use acfd_utils, only: vxao_to_vxmo
    use acfd_wf, only: orb, eig, no, nv
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, pointer :: ibase(:)

    ! Test Stuff
    real(dp), intent(in) :: eps
    real(dp), intent(out) :: vpertx(ntdg), drhsx

    real(dp), pointer :: int_fij(:, :, :) ! aux,occ,occ integrals
    real(dp), pointer :: int_fai(:, :, :) ! aux,virt,occ integrals

    type(vector_type) :: vref_ao
    real(dp), pointer :: auxcharg(:), homo(:)
    real(dp), pointer :: vxmo(:), vsol(:)
    real(dp) :: vxl(ntdg)

    ibase => memory_allocate_integer(1)

    vxl = vx

    ! get number of auxiliary basis functions
    call basis_number('OEP', naux_oep)

    ! calculate occ,occ and virt,occ 3-index-integrals
    int_fij => memory_allocate(naux_oep, no, no)
    call build_fij('OEP', 'J', orb, no, nv, naux_oep, int_fij)
    int_fai => memory_allocate(naux_oep, nv, no)
    call build_fai('OEP', 'J', orb, no, nv, naux_oep, int_fai)

    ! transform nonlocal exchange to MO basis
    vxmo => memory_allocate(ntqg)
    call vxao_to_vxmo(vx, orb, ntdg, ntqg, vxmo)

    auxcharg => memory_allocate(naux_oep)
    homo => memory_allocate(naux_oep) ! dummy variable for basis_process_scexx
    vsol => memory_allocate(naux_oep)

    ! process auxiliary basis set
    call basis_process_scexx(trans_mat_constraint, int_fai, auxcharg, &
                             int_fij, vref_ao, vref_aux_dummy, homo, vxmo, eig)

    naux_constraint = trans_mat_constraint%m

    ! calculate vx local
    call exx_oep_pert(orb, eig, vxl, vxmo, no, naux_oep, &
                      naux_constraint, trans_mat_constraint, vref_ao, &
                      vsol, auxcharg, int_fai, eps, vpertx, drhsx)

    call memory_release(int_fij)

    call memory_release(ibase)

  end subroutine scexx_opm_test

  subroutine exx_oep_pert(orb, eig, vx, vxmo, noo, naux, &
                          naux_constraint, trans_mat, vref_ao, vsol, &
                          auxcharg, int_fai, eps, vpertx, drhsx)
    use common_cbas, only: ntdg, ntqg, ntg
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_derived_datatypes, only: matrix_type, vector_type
    use acfd_parameters, only: thr_solve
    use acfd_utils, only: subtract_vref, x0_exx, rhs_exx, oep_solve, &
                          test_orthogonality
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp), intent(in) :: eps
    real(dp), intent(out) :: drhsx, vpertx(ntdg)
    integer, intent(in) :: noo, naux, naux_constraint
    real(dp), intent(in) :: orb(ntqg), eig(ntg)
    real(dp), intent(inout) :: vx(ntdg)
    real(dp), intent(inout) :: vxmo(ntg, ntg)
    type(matrix_type), intent(in) :: trans_mat
    type(vector_type), intent(in) :: vref_ao
    real(dp), intent(in) :: int_fai(naux, ntg-noo, noo)
    real(dp), intent(out) :: vsol(naux)
    real(dp), intent(in) :: auxcharg(naux)

    ! Pointers for stack management
    integer, pointer :: ibase(:)

    ! Temporary variables
    real(dp), pointer :: r1_tmp1(:), r1_tmp2(:)

    integer :: nvo ! number of virtual orbitals
    real(dp), pointer :: rX0(:, :) ! response matrix and 1/(eig(virt)-eig(occ))
    real(dp), pointer :: rhs(:) ! right hand side
    real(dp), external :: ddot_x

    ! create a pointer to the bottom of the stack used in this driver
    ibase => memory_allocate_integer(1)

    nvo = ntg-noo

    ! subtract vref from local exchange
    call subtract_vref(vref_ao, orb, vxmo)

    ! construct and transform to processed basis response matrix
    rX0 => memory_allocate(naux_constraint, naux_constraint)
!    call x0_exx(eig, int_fai, trans_mat, noo, nvo, naux, naux_constraint, &
!                rX0)

    ! construct and transform to processed basis right hand side
    rhs => memory_allocate(naux_constraint)
!    call rhs_exx(eig, vxmo, trans_mat, int_fai, nvo, noo, naux, &
!                 naux_constraint, rhs)

    ! Perturbation Stuff for exchange
    r1_tmp1 => memory_allocate(naux_constraint)
    r1_tmp1 = eps
    r1_tmp2 => memory_allocate(naux_constraint)
    r1_tmp2 = -4d0*rhs
    drhsx = ddot_x(naux_constraint, r1_tmp2, 1, r1_tmp1, 1)
    call memory_release(r1_tmp2)
    r1_tmp2 => memory_allocate(naux)
    call dgemv_x('N', naux, naux_constraint, 1d0, trans_mat%mat, naux, &
                 r1_tmp1, 1, 0d0, r1_tmp2, 1)
    call basis_vec_contract_a('OEP', 'ORBITAL', 'J', r1_tmp2, vpertx, 0d0)
    call memory_release(r1_tmp1)

    ! solve oep equation
    call oep_solve(rX0, rhs, naux_constraint, thr_solve)

    ! backtransformation of solution vector
    call dgemv_x('N', naux, naux_constraint, 1d0, trans_mat%mat, naux, &
                 rhs, 1, 0d0, vsol, 1)

    ! transform local potential to AO basis
    r1_tmp1 => memory_allocate(ntdg)
    call basis_vec_contract_a('OEP', 'ORBITAL', 'J', vsol, r1_tmp1, 0d0)

    ! test orthogonality
    call test_orthogonality(auxcharg, vsol, naux, 'charge orthogonality')

    ! add reference potential
    call daxpy_x(ntdg, 1d0, vref_ao%mat, 1, r1_tmp1, 1)

    ! replace non-local exchange with local exchange
    vx = r1_tmp1

    call memory_release(ibase)

  end subroutine exx_oep_pert

  subroutine scexx_pert_exch(den0, fock0, vxnl0, vxl0, vpertx, dEx, drhox)
    use common_cbas, only: ntdg, ntqg, ntg
    use memory, only: memory_allocate, memory_release
    use acfd_wf, only: smh, no
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp), intent(in) :: den0(ntdg), fock0(ntqg)
    real(dp), intent(in) :: vxnl0(ntdg), vxl0(ntdg), vpertx(ntdg)
    real(dp), intent(out) :: dEx, drhox

    real(dp) :: E_xnl0, E_xnl, E_xl0, E_xl
    real(dp), pointer :: fockP(:)
    real(dp), pointer :: r1_tmp1(:), r1_tmp2(:)
    real(dp), pointer :: orbP(:), eigP(:), vxnlP(:), vjP(:), denP(:)
    real(dp), external :: spur

    fockP => memory_allocate(ntqg)
    fockP = fock0

    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp2 => memory_allocate(ntqg)
    call expan(vpertx, r1_tmp1, 1, 1, 1)
    call mmult(smh, r1_tmp1, r1_tmp2, 1, 1, 0, -1, 0, 0)
    call mmult(r1_tmp2, smh, r1_tmp1, 1, 1, 0, 0, 0, 0)
    fockP = fockP+r1_tmp1
    call memory_release(r1_tmp1)

    E_xnl0 = 0.5d0*spur(den0, vxnl0)
    E_xl0 = spur(den0, vxl0)

    orbP => memory_allocate(ntqg)
    eigP => memory_allocate(ntg)
    denP => memory_allocate(ntdg)
    r1_tmp1 => memory_allocate(ntqg)
    call diag(fockP, r1_tmp1, eigP, 0)
    call mmult(smh, r1_tmp1, orbP, 1, 1, 0, -1, 0, 0)

    call dgemm_x('N', 'T', ntg, ntg, no, 2d0, orbP, ntg, &
                 orbP, ntg, 0d0, r1_tmp1, ntg)
    call reduc(r1_tmp1, denP, 1, 1, 1d0, 0)
    call memory_release(r1_tmp1)

    vxnlP => memory_allocate(ntdg)
    vxnlP = 0d0
    vjP => memory_allocate(ntdg)
    vjP = 0d0
    r1_tmp1 => memory_allocate(1)
    call hfma(denP, denP, vjP, vxnlP, 1, 1, 0, .true., .true., &
              r1_tmp1, .false., .false.)
    call memory_release(r1_tmp1)

    E_xnl = 0.5d0*spur(denP, vxnlP)
    E_xl = spur(denP, vxl0)

    dEx = E_xnl-E_xnl0
    drhox = E_xl-E_xl0

    call memory_release(fockP)

  end subroutine scexx_pert_exch

  subroutine uscexx_pot_test(ntest, eps0)
    use common_cbas, only: ntdg, ntqg, ntg
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_opm, only: uscexx_opm
    use acfd_pot, only: focka, fockb, vxa, vxb, vj
    use acfd_wf, only: orba, dena, eiga, orbb, eigb, denb
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: ntest
    real(dp), intent(in) :: eps0

    real(dp), pointer :: orba_save(:), dena_save(:), eiga_save(:), &
                         focka_save(:), vxa_save(:), vj_save(:)
    real(dp), pointer :: orbb_save(:), denb_save(:), eigb_save(:), &
                         fockb_save(:), vxb_save(:)
    real(dp), pointer :: denc(:), dens(:)

    real(dp), pointer :: vxnla0(:), vxla0(:), vj0(:)
    real(dp), pointer :: vxnlb0(:), vxlb0(:)
    real(dp), pointer :: r1_tmp1(:)
    integer, pointer :: ibase(:)

    integer :: itest
    real(dp) :: eps
    real(dp), pointer :: diff_Exa(:), diff_rhoxa(:), diff_diffxa(:), &
                         diff_rhsxa(:)
    real(dp), pointer :: vec_perta(:), vpertxa(:)

    real(dp), pointer :: diff_Exb(:), diff_rhoxb(:), diff_diffxb(:), &
                         diff_rhsxb(:)
    real(dp), pointer :: vec_pertb(:), vpertxb(:)

    ibase => memory_allocate_integer(1)

    ! save properties
    orba_save => memory_allocate(ntqg)
    orba_save = orba
    dena_save => memory_allocate(ntdg)
    dena_save = dena
    eiga_save => memory_allocate(ntg)
    eiga_save = eiga
    focka_save => memory_allocate(ntqg)
    focka_save = focka
    vxa_save => memory_allocate(ntdg)
    vxa_save = vxa

    orbb_save => memory_allocate(ntqg)
    orbb_save = orbb
    denb_save => memory_allocate(ntdg)
    denb_save = denb
    eigb_save => memory_allocate(ntg)
    eigb_save = eigb
    fockb_save => memory_allocate(ntqg)
    fockb_save = fockb
    vxb_save => memory_allocate(ntdg)
    vxb_save = vxb

    vj_save => memory_allocate(ntdg)
    vj_save = vj

    ! calculate VxNL and Vj
    vj = 0d0
    vxa = 0d0
    vxb = 0d0
    vxnla0 => memory_allocate(ntdg)
    vxnlb0 => memory_allocate(ntdg)
    vj0 => memory_allocate(ntdg)

    denc => memory_allocate(ntdg)
    denc = dena+denb
    dens => memory_allocate(ntdg)
    dens = dena-denb
    r1_tmp1 => memory_allocate(ntqg*2)
    r1_tmp1(1:ntqg) = orba
    call hfma(dena, dena, vj, vxa, -1, 1, 0, .true., .true., &
              r1_tmp1, .false., .false.)
    r1_tmp1(1:ntqg) = orba
    r1_tmp1(1+ntqg:2*ntqg) = orbb
    call hfma(denc, dens, vj, vxb, -1, 1, 0, .true., .true., &
              r1_tmp1, .false., .false.)
    call memory_release(denc)
    vxb = vxa-vxb
    vj = (vj-(vxa+vxb))

    vxa = 2d0*vxa
    vxb = 2d0*vxb

    vxnla0 = vxa
    vxnlb0 = vxb
    vj0 = vj

    ! calculate VxL, Vc and ecorr0
    vxla0 => memory_allocate(ntdg)
    vxlb0 => memory_allocate(ntdg)
    call uscexx_opm(999)
    vxla0 = vxa
    vxlb0 = vxb

    ! Start Test
    diff_Exa => memory_allocate(ntest)
    diff_rhoxa => memory_allocate(ntest)
    diff_diffxa => memory_allocate(ntest)
    diff_rhsxa => memory_allocate(ntest)

    diff_Exb => memory_allocate(ntest)
    diff_rhoxb => memory_allocate(ntest)
    diff_diffxb => memory_allocate(ntest)
    diff_rhsxb => memory_allocate(ntest)

    vec_perta => memory_allocate(ntest)
    vpertxa => memory_allocate(ntdg)

    vec_pertb => memory_allocate(ntest)
    vpertxb => memory_allocate(ntdg)

    eps = eps0
    do itest = 1, ntest
      eps = eps*1d-1
      vxa = vxnla0
      vxb = vxnlb0
      call uscexx_opm_test(eps, vpertxa, vpertxb, diff_rhsxa(itest), &
                           diff_rhsxb(itest), 'both')
      call uscexx_pert_exch(dena, denb, focka, fockb, vxnla0, vxnlb0, &
                            vxla0, vxlb0, vpertxa, vpertxb, &
                            diff_Exa(itest), diff_rhoxa(itest))
    end do
    diff_rhsxa = 0.5d0*(diff_rhsxa+diff_rhsxb)

    diff_diffxa = diff_Exa-diff_rhoxa

    ! load properties
    orba = orba_save
    dena = dena_save
    eiga = eiga_save
    focka = focka_save
    vxa = vxa_save

    orbb = orbb_save
    denb = denb_save
    eigb = eigb_save
    fockb = fockb_save
    vxb = vxb_save

    vj = vj_save

    write (iout, '(/,1x,a,7x,a,9x,a,7x,a,5x,a,2x,a,3x,a)') &
      'Perturbation', 'd_Ex', 'd_rho_vx', 'd_rhs_x', &
      'd_Ex-d_rho_vx', 'd_Ex-d_rhs_x', 'd_rho_vx-d_rhs_x'
    eps = eps0
    do itest = 1, ntest
      eps = eps*1d-1
      write (iout, '(1x,7(es10.3,5x))') eps, diff_Exa(itest), &
        diff_rhoxa(itest), &
        diff_rhsxa(itest), diff_diffxa(itest), &
        diff_Exa(itest)-diff_rhsxa(itest), &
        diff_rhoxa(itest)-diff_rhsxa(itest)
    end do

    call memory_release(ibase)
  end subroutine uscexx_pot_test

  subroutine uscexx_opm_test(eps, vpertxa, vpertxb, drhsxa, drhsxb, spin)
    use common_cbas, only: ntqg, ntdg
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_auxiliary_basis, only: naux_oep, naux_constraint_a, &
                                    naux_constraint_b, trans_mat_constraint_a, &
                                    trans_mat_constraint_b
    use acfd_basis_process, only: basis_process_scexx
    use acfd_derived_datatypes, only: vector_type
    use acfd_integrals, only: build_fai, build_fij
    use acfd_pot, only: vxa, vxb
    use acfd_plot, only: vref_aux_dummy
    use acfd_utils, only: vxao_to_vxmo
    use acfd_wf, only: orba, orbb, eiga, eigb, noa, nob, nva, nvb
    use iso_fortran_env, only: dp => real64
    implicit none

    ! Test Stuff
    real(dp), intent(in) :: eps
    character(*), intent(in) :: spin
    real(dp), intent(out) :: vpertxa(ntdg), drhsxa
    real(dp), intent(out) :: vpertxb(ntdg), drhsxb

    integer, pointer :: ibase(:)

    ! int_fij_a         aux,occ,occ integrals for alpha spin
    ! int_fij_b         aux,occ,occ integrals for beta spin
    ! int_fai_a         aux,virt,occ integrals for alpha spin
    ! int_fai_b         aux,virt,occ integrals for beta spin
    real(dp), pointer :: int_fij_a(:, :, :), int_fij_b(:, :, :)
    real(dp), pointer :: int_fai_a(:, :, :), int_fai_b(:, :, :)

    ! vxmo        non-local exchange potential in MO basis
    ! vref_ao     reference potential in AO basis
    ! vsol_a      solution vector for alpha spin
    ! vsol_b      solution vector for beta spin
    type(vector_type) :: vref_ao
    real(dp), pointer :: auxcharg(:), homo(:)
    real(dp), pointer :: vxmo(:), vsol_a(:), vsol_b(:)
    real(dp) :: vxla(ntdg)
    real(dp) :: vxlb(ntdg)

    ibase => memory_allocate_integer(1)

    vxla = vxa
    vxlb = vxb

    ! get number of auxiliary basis functions
    call basis_number('OEP', naux_oep)

    ! calculate occ,occ and virt,occ 3-index-integrals
    int_fij_a => memory_allocate(naux_oep, noa, noa)
    int_fij_a = 0d0
    call build_fij('OEP', 'J', orba, noa, nva, naux_oep, int_fij_a)
    int_fai_a => memory_allocate(naux_oep, nva, noa)
    int_fai_a = 0d0
    call build_fai('OEP', 'J', orba, noa, nva, naux_oep, int_fai_a)

    !   processing of the auxiliary basis and
    !   calculation of the local exchange potential
    !   for ALPHA spin
    if (noa > 0) then

      ! transform nonlocal exchange to MO basis
      vxmo => memory_allocate(ntqg)
      call vxao_to_vxmo(vxa, orba, ntdg, ntqg, vxmo)

      auxcharg => memory_allocate(naux_oep)
      homo => memory_allocate(naux_oep)
      vsol_a => memory_allocate(naux_oep)

      ! process auxiliary basis set
      call basis_process_scexx(trans_mat_constraint_a, int_fai_a, auxcharg, &
                               int_fij_a, vref_ao, vref_aux_dummy, homo, &
                               vxmo, eiga)

      naux_constraint_a = trans_mat_constraint_a%m

      ! calculate vx local for alpha spin
      if (trim(spin) == 'alpha' .or. trim(spin) == 'both') then
        call exx_oep_pert(orba, eiga, vxla, vxmo, noa, &
                          naux_oep, naux_constraint_a, &
                          trans_mat_constraint_a, vref_ao, vsol_a, &
                          auxcharg, int_fai_a, &
                          eps, vpertxa, drhsxa)
      else
        call exx_oep_pert(orba, eiga, vxla, vxmo, noa, &
                          naux_oep, naux_constraint_a, &
                          trans_mat_constraint_a, vref_ao, vsol_a, &
                          auxcharg, int_fai_a, &
                          0d0, vpertxa, drhsxa)
      end if

      call memory_release(vxmo)

    else
      vxa = 0d0
    end if

    !   processing of the auxiliary basis and
    !   calculation of the local exchange potential
    !   for BETA spin
    if (nob > 0) then

      ! calculate occ,occ and virt,occ 3-index-integrals
      int_fij_b => memory_allocate(naux_oep, nob, nob)
      call build_fij('OEP', 'J', orbb, nob, nvb, naux_oep, int_fij_b)
      int_fai_b => memory_allocate(naux_oep, nvb, nob)
      call build_fai('OEP', 'J', orbb, nob, nvb, naux_oep, int_fai_b)

      ! transform nonlocal exchange to MO basis
      vxmo => memory_allocate(ntqg)
      call vxao_to_vxmo(vxb, orbb, ntdg, ntqg, vxmo)

      auxcharg => memory_allocate(naux_oep)
      homo => memory_allocate(naux_oep)
      vsol_b => memory_allocate(naux_oep)

      ! process auxiliary basis set
      call basis_process_scexx(trans_mat_constraint_b, int_fai_b, &
                               auxcharg, int_fij_b, vref_ao, vref_aux_dummy, &
                               homo, vxmo, eigb)

      naux_constraint_b = trans_mat_constraint_b%m

      ! calculate vx local for beta spin
      if (trim(spin) == 'beta' .or. trim(spin) == 'both') then
        call exx_oep_pert(orbb, eigb, vxlb, vxmo, nob, &
                          naux_oep, naux_constraint_b, &
                          trans_mat_constraint_b, vref_ao, vsol_b, &
                          auxcharg, int_fai_b, &
                          eps, vpertxb, drhsxb)
      else
        call exx_oep_pert(orbb, eigb, vxlb, vxmo, nob, &
                          naux_oep, naux_constraint_b, &
                          trans_mat_constraint_b, vref_ao, vsol_b, &
                          auxcharg, int_fai_b, &
                          0d0, vpertxb, drhsxb)
      end if

      call memory_release(vxmo)

    else
      vxb = 0d0
    end if

    call memory_release(int_fij_a)

    call memory_release(ibase)

  end subroutine uscexx_opm_test

  subroutine uscexx_pert_exch(dena0, denb0, focka0, fockb0, &
                              vxnla0, vxnlb0, vxla0, vxlb0, &
                              vpertxa, vpertxb, dEx, drhox)
    use common_cbas, only: ntdg, ntqg, ntg
    use memory, only: memory_allocate, memory_release
    use acfd_wf, only: smh, noa, nob
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp), intent(in) :: dena0(ntdg), focka0(ntqg)
    real(dp), intent(in) :: vxnla0(ntdg), vxla0(ntdg), vpertxa(ntdg)
    real(dp), intent(in) :: denb0(ntdg), fockb0(ntqg)
    real(dp), intent(in) :: vxnlb0(ntdg), vxlb0(ntdg), vpertxb(ntdg)
    real(dp), intent(out) :: dEx, drhox

    real(dp) :: E_xnl0, E_xnl, E_xl0, E_xl
    real(dp), pointer :: fockaP(:), fockbP(:)
    real(dp), pointer :: r1_tmp1(:), r1_tmp2(:)
    real(dp), pointer :: orbaP(:), eigaP(:), vxnlaP(:), vjP(:), denaP(:)
    real(dp), pointer :: orbbP(:), eigbP(:), vxnlbP(:), denbP(:)
    real(dp), pointer :: dens(:), denc(:)
    real(dp), external :: spur

    fockaP => memory_allocate(ntqg)
    fockaP = focka0
    fockbP => memory_allocate(ntqg)
    fockbP = fockb0

    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp2 => memory_allocate(ntqg)
    call expan(vpertxa, r1_tmp1, 1, 1, 1)
    call mmult(smh, r1_tmp1, r1_tmp2, 1, 1, 0, -1, 0, 0)
    call mmult(r1_tmp2, smh, r1_tmp1, 1, 1, 0, 0, 0, 0)
    fockaP = fockaP+r1_tmp1
    call memory_release(r1_tmp1)

    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp2 => memory_allocate(ntqg)
    call expan(vpertxb, r1_tmp1, 1, 1, 1)
    call mmult(smh, r1_tmp1, r1_tmp2, 1, 1, 0, -1, 0, 0)
    call mmult(r1_tmp2, smh, r1_tmp1, 1, 1, 0, 0, 0, 0)
    fockbP = fockbP+r1_tmp1
    call memory_release(r1_tmp1)

    E_xnl0 = 0.5d0*(spur(dena0, vxnla0)+spur(denb0, vxnlb0))
    E_xl0 = spur(dena0, vxla0)+spur(denb0, vxlb0)

    orbaP => memory_allocate(ntqg)
    eigaP => memory_allocate(ntg)
    denaP => memory_allocate(ntdg)
    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp1 = 0d0
    call diag(fockaP, r1_tmp1, eigaP, 0)
    call mmult(smh, r1_tmp1, orbaP, 1, 1, 0, -1, 0, 0)

    r1_tmp1 = 0d0
    call dgemm_x('N', 'T', ntg, ntg, noa, 2d0, orbaP, ntg, &
                 orbaP, ntg, 0d0, r1_tmp1, ntg)
    call reduc(r1_tmp1, denaP, 1, 1, 1d0, 0)
    call memory_release(r1_tmp1)

    orbbP => memory_allocate(ntqg)
    eigbP => memory_allocate(ntg)
    denbP => memory_allocate(ntdg)
    r1_tmp1 => memory_allocate(ntqg)
    call diag(fockbP, r1_tmp1, eigbP, 0)
    call mmult(smh, r1_tmp1, orbbP, 1, 1, 0, -1, 0, 0)

    r1_tmp1 = 0d0
    call dgemm_x('N', 'T', ntg, ntg, nob, 2d0, orbbP, ntg, &
                 orbbP, ntg, 0d0, r1_tmp1, ntg)
    call reduc(r1_tmp1, denbP, 1, 1, 1d0, 0)
    call memory_release(r1_tmp1)

    vxnlaP => memory_allocate(ntdg)
    vxnlaP = 0d0
    vxnlbP => memory_allocate(ntdg)
    vxnlbP = 0d0
    vjP => memory_allocate(ntdg)
    vjP = 0d0

    denc => memory_allocate(ntdg)
    denc = denaP+denbP
    dens => memory_allocate(ntdg)
    dens = denaP-denbP
    r1_tmp1 => memory_allocate(ntqg*2)
    r1_tmp1 = 0d0
    r1_tmp1(1:ntqg) = orbaP
    call hfma(denaP, denaP, vjP, vxnlaP, -1, 1, 0, .true., .true., &
              r1_tmp1, .false., .false.)
    r1_tmp1(1:ntqg) = orbaP
    r1_tmp1(1+ntqg:2*ntqg) = orbbP
    call hfma(denc, dens, vjP, vxnlbP, -1, 1, 0, .true., .true., &
              r1_tmp1, .false., .false.)
    call memory_release(denc)
    vxnlbP = vxnlaP-vxnlbP
    vjP = (vjP-(vxnlaP+vxnlbP))

    E_xnl = 0.25d0*(spur(denaP, vxnlaP)+spur(denbP, vxnlbP))
    E_xl = 0.5d0*(spur(denaP, vxla0)+spur(denbP, vxlb0))
    dEx = E_xnl-E_xnl0
    drhox = E_xl-E_xl0

    call memory_release(fockaP)

  end subroutine uscexx_pert_exch

end module acfd_pot_test
