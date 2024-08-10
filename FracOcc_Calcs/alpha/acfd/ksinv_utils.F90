module ksinv_utils
  implicit none
  private
  public :: get_u_frozen, get_orb_and_den, linear_response, get_energies, &
            den_from_orb, get_potentials

contains

  ! Inversion and EXX OEP (if required), EXX OEP is needed to plot v_c
  subroutine linear_response(eig, den, ref_den, orb, vxc, vxc_sol, vx_sol, &
                             c_vsol_Delta, inversion_rhs, no, iter, cmd, nob)
    use common_cbas, only: ntg, ntqg, ntdg
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use acfd_derived_datatypes, only: vector_type
    use acfd_integrals, only: build_fai, build_all
    use acfd_plot, only: vref_aux_up
    use acfd_utils, only: exx_oep, exx_oep_spin_symmetrized, vxao_to_vxmo
    use acfd_wf, only: int_all
    use ksinv_basis_process, only: basis_process
    use ksinv_mod, only: vx, naux_oep, naux_constraint, &
                         trans_mat_constraint, &
                         trans_mat_constraint_inv, &
                         int_fkl, l_plot_always, l_stopped, &
                         thr_scf_inversion
    use iso_fortran_env, only: dp => real64
    implicit none
    ! identifies the algorithm which called the subroutine (UKSINV or KSINV)
    character(64), intent(in) :: cmd
    integer, intent(in) :: iter ! number of current iteration
    integer, intent(in) :: no ! number of occupied orbitals
    integer, intent(in), optional :: nob ! number of occupied orbitals for second spin-channel
    integer, pointer :: ibase(:) ! pointer for stack management
    real(dp), intent(inout) :: vxc(ntdg) ! exchange correlation potential in AO basis
    real(dp), intent(inout) :: vxc_sol(naux_oep), vx_sol(naux_oep) ! exchange and exchange correlation potential in auxiliary basis
    real(dp), intent(in) :: ref_den(ntdg), den(ntdg) !  RDM1 in AO basis
    real(dp), intent(in) :: eig(ntg), orb(ntdg)
    real(dp), intent(inout) :: inversion_rhs, c_vsol_Delta
    real(dp), pointer :: int_fai(:, :, :), int_fai_b(:, :, :) ! 3-index integrals (aux,virt,occ)
    type(vector_type) :: vref_ao ! reference potential in AO basis
    real(dp), pointer :: auxcharg(:)
    real(dp), pointer :: vxmo(:), vxmoa(:), vxmob(:) ! non-local exchange potential in MO basis
    real(dp), pointer :: homo(:)
    real(dp), pointer :: dena(:), denb(:), vxa(:), vxb(:)

    ibase => memory_allocate_integer(1)

    auxcharg => memory_allocate(naux_oep)
    homo => memory_allocate(naux_oep)

    int_all => memory_allocate(naux_oep, ntg, ntg)
    call build_all('OEP', 'J', orb, naux_oep, int_all)

    ! calculate occ,occ and virt,occ 3-index-integrals
    int_fai => memory_allocate(naux_oep, ntg-no, no)
    call build_fai('OEP', 'J', orb, no, ntg-no, naux_oep, int_fai)

    if (present(nob)) then
      int_fai_b => memory_allocate(naux_oep, ntg-nob, nob)
      call build_fai('OEP', 'J', orb, nob, ntg-nob, naux_oep, int_fai_b)

      ! calculate the basis processing matrices
      call basis_process(naux_oep, vref_aux_up, vref_ao, &
                         auxcharg=auxcharg, &
                         trans_mat_constraint=trans_mat_constraint, &
                         trans_mat_constraint_inv=trans_mat_constraint_inv, &
                         int_fai=int_all, int_fai_b=int_all, eig=eig)
      naux_constraint = trans_mat_constraint%m

      ! solve linear response equation for difference density
      call potential_inversion(vxc, vxc_sol, naux_oep, naux_constraint, eig, &
                               trans_mat_constraint, int_all, int_fkl, &
                               den, ref_den, iter, trans_mat_constraint_inv, &
                               inversion_rhs, c_vsol_Delta, cmd, &
                               int_fai_b=int_all, noa=no, nob=nob)

    else
      ! calculate the basis processing matrices
      call basis_process(naux_oep, vref_aux_up, vref_ao, &
                         auxcharg=auxcharg, &
                         trans_mat_constraint=trans_mat_constraint, &
                         trans_mat_constraint_inv=trans_mat_constraint_inv, &
                         int_fai=int_all, eig=eig)

      naux_constraint = trans_mat_constraint%m

      ! solve linear response equation for difference density
      call potential_inversion(vxc, vxc_sol, naux_oep, naux_constraint, eig, &
                               trans_mat_constraint, int_fai, int_fkl, &
                               den, ref_den, iter, trans_mat_constraint_inv, &
                               inversion_rhs, c_vsol_Delta, cmd, noa=no)
    end if

    ! if needed, calculate vx-local for plotting
    if (l_plot_always .or. (inversion_rhs < thr_scf_inversion) .or. &
        l_stopped) then

      if (present(nob)) then

        dena => memory_allocate(ntdg)
        denb => memory_allocate(ntdg)
        vxa => memory_allocate(ntdg)
        vxb => memory_allocate(ntdg)
        call den_from_orb(orb, no, dena) ! get alpha density
        call den_from_orb(orb, nob, denb) ! get beta density
        call get_potentials(dena, denb, vxa, vxb) ! get exchange potentials

        vxmoa => memory_allocate(ntqg)
        call vxao_to_vxmo(vxa, orb, ntdg, ntqg, vxmoa)
        vxmob => memory_allocate(ntqg)
        call vxao_to_vxmo(vxb, orb, ntdg, ntqg, vxmob)

!        call exx_oep_spin_symmetrized(vxa, vxb, vx_sol, &
!                                      no, nob, naux_oep, naux_constraint, &
!                                      orb, orb, eig, eig, vxmoa, vxmob, &
!                                      trans_mat_constraint, vref_ao, auxcharg, &
!                                      homo, int_fai, int_fai_b, &
!                                      vref_aux_up, dena, denb, '')
        call memory_release(dena)

      else

        vxmo => memory_allocate(ntqg)
        call vxao_to_vxmo(vx, orb, ntdg, ntqg, vxmo)
        call exx_oep(vx, vx_sol, no, naux_oep, naux_constraint, orb, eig, &
                     vxmo, trans_mat_constraint, vref_ao, auxcharg, homo, &
                     int_fai, '', vref_aux_up, den, '')
      end if
    end if

    call memory_release(ibase)

  end subroutine linear_response

  ! Inverts single particle density matrix to calculate v_{xc}
  subroutine potential_inversion(vxc_acc_ao, vxc_acc_oep, naux, &
                                 naux_constraint, eig, trans_mat, int_fai, &
                                 int_fkl, den, ref_den, iter, trans_mat_inv, &
                                 inversion_rhs, c_vsol_Delta, cmd, int_fai_b, noa, nob)
    use common_cbas, only: ntg, ntdg, ntqg
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_derived_datatypes, only: matrix_type
    use acfd_utils, only: oep_solve, x0_exx
    use ksinv_mod, only: thr_solve, iverbose, &
                         l_backfilter, mixing_a, mix_switch_iter
    use iso_fortran_env, only: dp => real64
    implicit none
    ! identifies the algorithm which called the subroutine initially (UKSINV or KSINV)
    character(64), intent(in) :: cmd
    integer, intent(in) :: naux ! number of auxiliary basis functions in unproccessed OEP basis
    integer, intent(in) :: naux_constraint ! number of auxiliary basis functions in processed OEP basis
    integer, intent(in) :: iter ! current iteration
    integer, intent(in), optional :: noa, nob
    real(dp), intent(in) :: eig(ntg) ! eigenvalues
    real(dp), intent(in) :: den(ntdg) ! density
    real(dp), intent(in) :: ref_den(ntdg) ! reference density
    real(dp), intent(out) :: vxc_acc_ao(ntdg) ! accumulated exchange correlation potential in AO basis
    real(dp), intent(inout) :: inversion_rhs, c_vsol_Delta ! abusolute vaulue of to be minimized vector |\Delta \rho|
    real(dp), intent(inout) :: vxc_acc_oep(naux) ! accumulated exchange correlation potential in OEP basis
    type(matrix_type), intent(in) :: trans_mat ! transformation matrix
    type(matrix_type), intent(in) :: trans_mat_inv ! inverse transformation matrix used for backfiltering
    real(dp), intent(in) :: int_fai(:, :, :) ! 3-index integrals (oep-function | orbital(virtual), orbial(occupied)
    real(dp), intent(in), optional :: int_fai_b(:, :, :)
    real(dp), intent(in) :: int_fkl(naux, ntg, ntg) ! 3-index integrals (oep-function | ao-function, ao-function

    integer, pointer :: ibase(:) ! pointer for stack management

    integer :: iaux ! loop index
    real(dp), pointer :: rX0(:, :), rX0_b(:, :) ! response matrix
    real(dp), pointer :: rhs(:) ! right hand side
    real(dp), pointer :: vxc_delta_oep(:) ! incremental change of potential
    real(dp) :: mixing
    real(dp) :: P_delta(ntdg) ! current KS density-matrux minus reference density-matrix, upper triangle
    real(dp) :: P_delta_sqr(ntqg) ! expanded P_delta, square matrix

    ! auxiliary variables
    real(dp), pointer :: r1_tmp1(:)
    real(dp), pointer :: r2_tmp1(:, :)

    real(dp), external :: trace ! function that returns trace of a matrix

    ! create a pointer to the bottom of the stack used in this driver
    ibase => memory_allocate_integer(1)

    ! Generate P_Delta and expand
    P_delta = den-ref_den
    call expan(P_delta, P_Delta_sqr, 1, 1, 0)

    !
    ! Generate right hand side
    !
    ! P^{/Delta} is transferred into the unprocessed OEP basis
    ! P^{/Delta}_p = \sum_s \sum_t (f_p|\chi_{\kappa}\chi_{\lambda}) P^{\Delta}_{\kappa,\lambda}
    !
    rhs => memory_allocate(naux_constraint)
    rhs = 0d0
    r1_tmp1 => memory_allocate(naux)
    r1_tmp1 = 0d0
    r2_tmp1 => memory_allocate(ntg, ntg)
    r2_tmp1 = 0d0
    do iaux = 1, naux
      call dgemm_x('N', 'T', ntg, ntg, ntg, 1d0, &
                   P_Delta_sqr, ntg, &
                   int_fkl(iaux, :, :), ntg, &
                   0d0, r2_tmp1, ntg)
      r1_tmp1(iaux) = trace(r2_tmp1, ntg)
    end do

    ! transform right hand side to processed oep basis
    call dgemv_x('T', naux, naux_constraint, 1d0,&
                &   trans_mat%mat, naux, r1_tmp1, 1, 0d0, rhs, 1)
    inversion_rhs = sqrt(dot_product(rhs, rhs))
    call memory_release(r1_tmp1)

    ! construct and transform to processed basis response matrix
    rX0 => memory_allocate(naux_constraint, naux_constraint)
    call x0_exx(eig, trans_mat, noa, size(int_fai, 2), & 
                naux, naux_constraint, rX0)

    if (present(int_fai_b)) then
      rX0_b => memory_allocate(naux_constraint, naux_constraint)
      call x0_exx(eig, trans_mat, nob, &
                  size(int_fai_b, 2), naux, naux_constraint, rX0_b)
      rX0 = rX0+rX0_b
      call memory_release(rX0_b)
    end if

    if (trim(cmd) .eq. 'UKSINV') rX0 = 2d0*rX0
    if (trim(cmd) .eq. 'KSINV') then
      if (present(int_fai_b)) then
        rX0 = 2d0*rX0
      else
        rX0 = 4d0*rX0
      end if
    end if

    ! solve oep equation in oep processed basis
    call oep_solve(rx0, rhs, naux_constraint, thr_solve)

    ! backtransformation of solution vector to unprocessed oep basis
    vxc_delta_oep => memory_allocate(naux)

    call dgemv_x('N', naux, naux_constraint, 1d0, trans_mat%mat, naux, rhs, 1, &
                 0d0, vxc_delta_oep, 1)

    ! determine t_i
    if (iter == mix_switch_iter) c_vsol_Delta = inversion_rhs

    if (iter < mix_switch_iter) then
      mixing = mixing_a
    else
      mixing = 2d0/(1.0+(2d0/mixing_a-1d0)**(inversion_rhs/c_vsol_Delta))
    end if

    if (iverbose .gt. 0) write (iout, '(1x,a,1x,es20.12,/)') "mixing:", mixing

    vxc_acc_oep = vxc_acc_oep+mixing*vxc_delta_oep

    ! backfiltering: W x W^{-1} v_{xc}, removes remenants of now unused basis functions
    if (l_backfilter) then

      r1_tmp1 => memory_allocate(naux_constraint)

      call dgemm_x('N', 'N', naux_constraint, 1, naux, 1d0, &
                   trans_mat_inv%mat, naux_constraint, &
                   vxc_acc_oep, naux, 0d0, &
                   r1_tmp1, naux_constraint)

      call dgemm_x('N', 'N', naux, 1, naux_constraint, 1d0, &
                   trans_mat%mat, naux, &
                   r1_tmp1, naux_constraint, 0d0, &
                   vxc_acc_oep, naux)

      call memory_release(r1_tmp1)

    end if

    ! transform local potential to AO basis
    call basis_vec_contract_A('OEP', 'ORBITAL', 'J', vxc_acc_oep, vxc_acc_ao, &
                              0d0)

    call memory_release(ibase)

  end subroutine potential_inversion

  ! generate u_frozen corresponding to reference denisty
  subroutine get_u_frozen()
    use common_cbas, only: ntg
    use memory, only: memory_allocate, memory_release
    use ksinv_mod, only: ref_den, int_fkl, u_frozen, naux_oep
    use iso_fortran_env, only: dp => real64
    implicit none
    integer :: iaux
    real(dp), pointer :: r2_tmp1(:, :), r2_tmp2(:, :)
    real(dp), external :: trace

    ! prepare u_frozen
    r2_tmp1 => memory_allocate(ntg, ntg)
    r2_tmp1 = 0d0
    r2_tmp2 => memory_allocate(ntg, ntg)
    r2_tmp2 = 0d0
    call expan(ref_den, r2_tmp1, 1, 1, 0)
    do iaux = 1, naux_oep
      call dgemm_x('N', 'T', ntg, ntg, ntg, 1d0, &
                   r2_tmp1, ntg, &
                   int_fkl(iaux, :, :), ntg, &
                   0d0, r2_tmp2, ntg)
      u_frozen(iaux) = trace(r2_tmp2, ntg)
    end do
    call memory_release(r2_tmp1)

  end subroutine get_u_frozen

  subroutine den_from_orb(orb, no, den)
    use common_cbas, only: ntg, ntdg, ntqg
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: no
    real(dp), intent(in) :: orb(ntdg)
    real(dp), intent(out) :: den(ntdg)
    real(dp), pointer :: r1_tmp1(:)

    den = 0d0
    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp1 = 0d0
    call dgemm_x('N', 'T', ntg, ntg, no, 1d0, &
                 orb, ntg, &
                 orb, ntg, &
                 0d0, r1_tmp1, ntg)

    call reduc(r1_tmp1, den, 1, 1, 1d0, 0)
    call memory_release(r1_tmp1)

  end subroutine den_from_orb

  !
  ! generates density matrix and orbitals and eigenvalues from fock_tri, which is the upper triangular matrix
  ! of the current fock matrix
  !
  ! Steps are:
  !
  ! Expands fock matrix
  ! Translates Generalized Eigenvalue Problem (GEP) into Eigenvalue Problem (EP)
  ! calculates eigenvectors and eigenvalues of EP
  ! Transforms eigenvectors of EP  back to orbitals GEP
  ! Uses Eigenvectors to generate Density-Matrix
  !
  subroutine get_orb_and_den(fock_tri, orb, den, smh, eig, no, cmd, nob)
    use common_cbas, only: ntg, ntdg, ntqg
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    character(64), intent(in) :: cmd
    integer, intent(in)   :: no
    integer, intent(in), optional   :: nob
    real(dp), intent(in)  :: fock_tri(ntdg)
    real(dp), intent(in)  :: smh(ntqg)
    real(dp), intent(out) :: eig(ntg)
    real(dp), intent(out) :: orb(ntqg)
    real(dp), intent(out) :: den(ntdg)
    real(dp), pointer :: r1_tmp1(:)
    real(dp), pointer :: fock(:)

    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp1 = 0d0

    ! Expand tridiagonal matrix into full Fock matrix
    fock => memory_allocate(ntqg)
    fock = 0d0
    call expan(fock_tri, fock, 1, 1, 1)

    ! S^{-1/2} F S^{-1/2}
    call dgemm_x('N', 'N', ntg, ntg, ntg, 1d0, &
                 smh, ntg, &
                 fock, ntg, &
                 0d0, r1_tmp1, ntg)

    call dgemm_x('N', 'N', ntg, ntg, ntg, 1d0, &
                 r1_tmp1, ntg, &
                 smh, ntg, &
                 0d0, fock, ntg)

    ! solve S^{-1/2} F S^{-1/2} S^{1/2} c_i = e_i S^{1/2} c_i
    call diag(fock, orb, eig, 0)

    ! S^{-1/2} S^{1/2} c_i = c_i
    call dgemm_x('N', 'N', ntg, ntg, ntg, 1d0, &
                 smh, ntg, &
                 orb, ntg, &
                 0d0, r1_tmp1, ntg)

    orb = r1_tmp1

    call memory_release(r1_tmp1)

    ! calculate density
    call den_from_orb(orb, no, den)

    if (present(nob)) then
      r1_tmp1 => memory_allocate(ntdg)
      call den_from_orb(orb, nob, r1_tmp1) ! get beta density
      den = den+r1_tmp1
    else if (cmd .eq. 'KSINV') then
      den = 2d0*den
    end if

  end subroutine get_orb_and_den

  subroutine get_energies(den, h0, v_kin, vx, vj, E0, E_K, E_h0, E_kin, E_x, &
                          E_j, E_ee, E_ext, ntdg)
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: ntdg !size of a upper triangular matrix in the AO basis set
    real(dp), intent(in) :: den(ntdg), h0(ntdg), v_kin(ntdg), vx(ntdg), &
                            vj(ntdg), E0, E_K !density matrix and matrix representations of operators
    real(dp), intent(out) :: E_h0, E_kin, E_x, E_j, E_ee, E_ext ! energies
    real(dp), external :: spur ! trace of the matrix product of two triagonal matrices

    E_h0 = spur(den, h0)
    E_KIN = spur(den, v_kin)
    E_X = 0.5d0*spur(den, vx)
    E_j = 0.5d0*spur(den, vj)

    E_ext = E_h0-E_KIN
    E_ee = E0-E_h0-E_K

  end subroutine get_energies

  subroutine get_potentials(dena, denb, vxa, vxb)
    use common_cbas, only: ntqg
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: dena(:), denb(:)
    real(dp), intent(out) :: vxa(:), vxb(:)
    real(dp), pointer :: vj_aux(:), r1_tmp1(:)

    vj_aux => memory_allocate(size(vxa, 1))
    vj_aux = 0d0
    vxa = 0d0
    vxb = 0d0
    r1_tmp1 => memory_allocate(ntqg*2)
    r1_tmp1 = 0d0
    call hfma(dena, dena, vj_aux, vxa, -1, 1, 0, .true., .true., &
              r1_tmp1, .false., .false.)
    call hfma(dena+denb, dena-denb, vj_aux, vxb, -1, 1, 0, .true., .true., &
              r1_tmp1, .false., .false.)
    vxb = vxa-vxb
    vj_aux = (vj_aux-(vxa+vxb))
    vxa = 2d0*vxa
    vxb = 2d0*vxb

    call memory_release(vj_aux)

  end subroutine get_potentials

end module ksinv_utils
