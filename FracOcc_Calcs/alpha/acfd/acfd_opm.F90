module acfd_opm

contains

  ! EXX OEP calculation in spin-restricted case
  subroutine scexx_opm(iter)
    use common_cbas, only: ntqg, ntdg
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use acfd_auxiliary_basis, only: naux_oep, naux_constraint, &
                                    trans_mat_constraint
    use acfd_basis_process, only: basis_process_scexx
    use acfd_derived_datatypes, only: vector_type
    use acfd_integrals, only: build_fai, build_fij
    use acfd_plot, only: vref_aux_up
    use acfd_pot, only: vx, vsol
    use acfd_utils, only: vxao_to_vxmo, exx_oep
    use acfd_wf, only: no, nv, orb, eig, den
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: iter ! number of current iteration

    integer, pointer :: ibase(:) ! pointer for stack management

    real(dp), pointer :: int_fij(:, :, :) ! 3-index integrals (aux,occ,occ)
    real(dp), pointer :: int_fai(:, :, :) ! 3-index integrals (aux,virt,occ)
    character(3) :: c_iter ! number of current iteration as character for filenames
    type(vector_type) :: vref_ao ! reference potential in AO basis
    real(dp), pointer :: auxcharg(:), homo(:)
    real(dp), pointer :: vxmo(:) ! non-local exchange potential in MO basis

    ! allocate integer as the starting point of memory stack used by present subroutine
    ibase => memory_allocate_integer(1)

    ! convert number of iteration to character for use in filenames
    write (c_iter, '(I3.3)') iter

    ! determine number of basis functions in OEP basis
    call basis_number('OEP', naux_oep)

    ! calculate occ,occ and virt,occ 3-index-integrals
    int_fij => memory_allocate(naux_oep, no, no)
    call build_fij('OEP', 'J', orb, no, nv, naux_oep, int_fij)
    int_fai => memory_allocate(naux_oep, nv, no)
    call build_fai('OEP', 'J', orb, no, nv, naux_oep, int_fai)

    ! transform nonlocal exchange potential to MO basis
    vxmo => memory_allocate(ntqg)
    call vxao_to_vxmo(vx, orb, ntdg, ntqg, vxmo)

    auxcharg => memory_allocate(naux_oep)
    homo => memory_allocate(naux_oep)

    ! process OEP basis set
    call basis_process_scexx(trans_mat_constraint, int_fai, auxcharg, int_fij, &
                             vref_ao, vref_aux_up, homo, vxmo, eig)
    naux_constraint = trans_mat_constraint%m

    ! calculate vx local
    call exx_oep(vx, vsol, no, naux_oep, naux_constraint, orb, eig, vxmo, &
                 trans_mat_constraint, vref_ao, auxcharg, homo, int_fai, '', &
                 vref_aux_up, den, c_iter)

    call memory_release(int_fij)

    call memory_release(ibase)

  end subroutine scexx_opm

  ! EXX OEP calculation in spin-unrestricted case
  subroutine uscexx_opm(iter)
    use common_cbas, only: ntqg, ntdg
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use acfd_auxiliary_basis, only: naux_constraint_a, naux_constraint_b, &
                                    trans_mat_constraint_a, &
                                    trans_mat_constraint_b, &
                                    naux_oep
    use acfd_derived_datatypes, only: vector_type
    use acfd_basis_process, only: basis_process_scexx
    use acfd_integrals, only: build_fai, build_fij
    use acfd_plot, only: vref_aux_upa, vref_aux_upb
    use acfd_pot, only: vxa, vxb, vsola, vsolb
    use acfd_utils, only: vxao_to_vxmo, exx_oep
    use acfd_wf, only: noa, nob, nva, nvb, orba, orbb, eiga, eigb, dena, denb
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: iter

    integer, pointer :: ibase(:)

    real(dp), pointer :: int_fij_a(:, :, :), int_fij_b(:, :, :) ! aux,occ,occ integrals for alpha and beta spin
    real(dp), pointer :: int_fai_a(:, :, :), int_fai_b(:, :, :) ! aux,virt,occ integrals for alpha and beta spin

    character(3) :: c_iter ! number of current iteration as character for filenames

    type(vector_type) :: vref_ao ! reference potential in AO basis
    real(dp), pointer :: auxcharg(:), homo(:)
    real(dp), pointer :: vxmo(:) ! non-local exchange potential in MO basis

    ibase => memory_allocate_integer(1)

    ! convert number of current iteration to character for use in filenames
    write (c_iter, '(I3.3)') iter

    ! get number of functions in OEP basis set
    call basis_number('OEP', naux_oep)

    ! calculate occ,occ and virt,occ 3-index-integrals
    int_fij_a => memory_allocate(naux_oep, noa, noa)
    call build_fij('OEP', 'J', orba, noa, nva, naux_oep, int_fij_a)
    int_fai_a => memory_allocate(naux_oep, nva, noa)
    call build_fai('OEP', 'J', orba, noa, nva, naux_oep, int_fai_a)

    if (nob > 0) then
      int_fij_b => memory_allocate(naux_oep, nob, nob)
      call build_fij('OEP', 'J', orbb, nob, nvb, naux_oep, int_fij_b)
      int_fai_b => memory_allocate(naux_oep, nvb, nob)
      call build_fai('OEP', 'J', orbb, nob, nvb, naux_oep, int_fai_b)
    end if

    ! processing of OEP basis and calculation of the local exchange potential
    ! for alpha spin
    if (noa > 0) then

      ! transform nonlocal exchange to MO basis
      vxmo => memory_allocate(ntqg)
      call vxao_to_vxmo(vxa, orba, ntdg, ntqg, vxmo)

      auxcharg => memory_allocate(naux_oep)
      homo => memory_allocate(naux_oep)

      ! process auxiliary basis set
      if (nob > 0) then
        call basis_process_scexx(trans_mat_constraint_a, int_fai_a, &
                                 auxcharg, int_fij_a, vref_ao, vref_aux_upa, &
                                 homo, vxmo, eiga, int_fij_b=int_fij_b)
      else
        call basis_process_scexx(trans_mat_constraint_a, int_fai_a, auxcharg, &
                                 int_fij_a, vref_ao, vref_aux_upa, homo, &
                                 vxmo, eiga)
      end if

      naux_constraint_a = trans_mat_constraint_a%m

      ! calculate vx local for alpha spin
      call exx_oep(vxa, vsola, noa, naux_oep, naux_constraint_a, orba, eiga, &
                   vxmo, trans_mat_constraint_a, vref_ao, auxcharg, homo, &
                   int_fai_a, 'a', vref_aux_upa, dena, c_iter)

      call memory_release(vxmo)

    else
      vxa = 0d0
    end if

    ! processing of OEP basis and calculation of the local exchange potential
    ! for beta spin
    if (nob > 0) then

      ! transform nonlocal exchange to MO basis
      vxmo => memory_allocate(ntqg)
      call vxao_to_vxmo(vxb, orbb, ntdg, ntqg, vxmo)

      auxcharg => memory_allocate(naux_oep)
      homo => memory_allocate(naux_oep)

      ! process OEP basis set
      call basis_process_scexx(trans_mat_constraint_b, int_fai_b, auxcharg, &
                               int_fij_b, vref_ao, vref_aux_upb, homo, vxmo, &
                               eigb, int_fij_b=int_fij_a)

      naux_constraint_b = trans_mat_constraint_b%m

      ! calculate vx local for beta spin
      call exx_oep(vxb, vsolb, nob, naux_oep, naux_constraint_b, orbb, eigb, &
                   vxmo, trans_mat_constraint_b, vref_ao, auxcharg, homo, &
                   int_fai_b, 'b', vref_aux_upb, denb, c_iter)

      call memory_release(vxmo)

    else
      vxb = 0d0
    end if

    call memory_release(int_fij_a)

    call memory_release(ibase)

  end subroutine uscexx_opm

  ! EXX OEP calculation in spin-unrestricted case with averaged spin
  subroutine uscexx_opm_spin_symmetrized(iter)
    use common_cbas, only: ntqg, ntdg, ntg
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use acfd_auxiliary_basis, only: naux_constraint_a, trans_mat_constraint_a, &
                                    naux_oep
    use acfd_derived_datatypes, only: vector_type
    use acfd_basis_process, only: basis_process_scexx
    use acfd_integrals, only: build_fai, build_fij, build_all
    use acfd_plot, only: vref_aux_up
    use acfd_pot, only: vxa, vxb, vsol, vxb1, vxb2
    use acfd_utils, only: vxao_to_vxmo, exx_oep_spin_symmetrized
    use acfd_wf, only: noa, nob, nva, nvb, orba, orbb, eiga, eigb, dena, denb, int_all
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: iter
    integer, pointer :: ibase(:)

    real(dp), pointer :: int_fij_a(:, :, :), int_fij_b(:, :, :) ! aux,occ,occ integrals for alpha and beta spin

    character(3) :: c_iter ! number of current iteration as character for filenames

    type(vector_type) :: vref_ao
    real(dp), pointer :: auxcharg(:), homo(:)
    real(dp), pointer :: vxmoa(:), vxmob(:), vxmob1(:), vxmob2(:) ! non-local exchange potential in MO basis

    ibase => memory_allocate_integer(1)

    ! convert number of iteration to character for use in filenames
    write (c_iter, '(I3.3)') iter

    ! get number of auxiliary basis functions
    call basis_number('OEP', naux_oep)

    ! calculate occ,occ and virt,occ 3-index-integrals
    int_fij_a => memory_allocate(naux_oep, noa, noa)
    call build_fij('OEP', 'J', orba, noa, nva, naux_oep, int_fij_a)
    int_fij_b => memory_allocate(naux_oep, nob, nob)
    call build_fij('OEP', 'J', orbb, nob, nvb, naux_oep, int_fij_b)

    call build_all('OEP', 'J', orba, naux_oep, int_all)

    ! transform nonlocal exchange to MO basis
    vxmoa => memory_allocate(ntqg)
    call vxao_to_vxmo(vxa, orba, ntdg, ntqg, vxmoa)
    vxmob => memory_allocate(ntqg)
    vxmob1 => memory_allocate(ntqg)
    vxmob2 => memory_allocate(ntqg)
    call vxao_to_vxmo(vxb, orbb, ntdg, ntqg, vxmob)
    call vxao_to_vxmo(vxb1, orba, ntdg, ntqg, vxmob1)
    call vxao_to_vxmo(vxb2, orba, ntdg, ntqg, vxmob2)
!    vxmob = 0.5d0*vxmob1 + 0.5d0*vxmob2

    auxcharg => memory_allocate(naux_oep)
    homo => memory_allocate(naux_oep)

    ! process auxiliary basis set
    call basis_process_scexx(trans_mat_constraint_a, int_all, auxcharg, &
                             int_fij_a, vref_ao, vref_aux_up, homo, vxmoa, &
                             eiga, int_fai_b=int_all, int_fij_b=int_fij_b)

    naux_constraint_a = trans_mat_constraint_a%m

    ! calculate vx local
    call exx_oep_spin_symmetrized(vxa, vxb, vsol, noa, nob, naux_oep, &
                                  naux_constraint_a, orba, orbb, eiga, eigb, &
                                  vxmob1, vxmob, trans_mat_constraint_a, &
                                  vref_ao, auxcharg, homo, vref_aux_up, dena, denb, c_iter, vxmob2)

    call memory_release(vxmoa)

    call memory_release(int_fij_a)

    call memory_release(ibase)

  end subroutine uscexx_opm_spin_symmetrized

end module acfd_opm
