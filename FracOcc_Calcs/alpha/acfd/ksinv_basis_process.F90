module ksinv_basis_process
  implicit none
  private
  public :: basis_process
contains
  !
  ! Processing of auxiliary basis set for KS inversion calculations
  !
  ! Some ideas behind the content can be found in
  ! J. Chem. Phys. 155, 054109 (2021); https://doi.org/10.1063/5.0056431
  ! J. Chem. Phys. 156, 204124 (2022); https://doi.org/10.1063/5.0087356
  !
  ! Parameters:
  ! naux [integer, in]: number of function in unprocessed auxiliary basis set
  ! vref_ao [vector_type, inout]: reference potential in AO basis
  ! vref_aux_up [vector_type, inout]: reference potential
  ! auxcharg [real(dp), out, optional]: charge constraint vector in unprocessed
  !     basis
  ! trans_mat_constraint [matrix_type, out, optional]: transformation matrix
  ! trans_mat_constraint_inv [matrix_type, out, optional]: transformation
  !     used for backfiltering during KS inversion
  ! int_fai [real(dp), in, optional]: 3-index integrals (basis,virt,occ)
  ! int_fai_b [real(dp), in, optional]: 3-index integrals (basis,virt,occ)
  !
  subroutine basis_process(naux, vref_aux_up, vref_ao, auxcharg, &
                           trans_mat_constraint, trans_mat_constraint_inv, &
                           int_fai, int_fai_b, eig)
    use common_cbas, only: ntdg
    use common_tapes, only: iout
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use acfd_basis_process, only: calc_amat, aevec_put_zeros, vref_const, &
                                  reduce_trans_mat_w_aevec, basis_reduc, proj, &
                                  charge_constr_vec, overlap_diag_unit
    use acfd_derived_datatypes, only: matrix_type, vector_type
    use acfd_utils, only: matrix_spectral
    use acfd_wf, only: occa, occb
    use ksinv_mod, only: thr_overlap_oep, thr_fai_oep, thr_symmetry, &
                         l_vref_fa, l_backfilter, iverbose
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: naux
    type(vector_type), intent(inout) :: vref_aux_up
    type(vector_type), intent(inout) :: vref_ao
    real(dp), intent(out), optional :: auxcharg(naux)
    type(matrix_type), intent(out), optional :: trans_mat_constraint
    type(matrix_type), intent(out), optional :: trans_mat_constraint_inv
    real(dp), intent(in), optional :: int_fai(:, :, :), int_fai_b(:, :, :)
    real(dp), intent(in), optional :: eig(:)

    integer, pointer :: ibase(:) ! pointer for stack management

    integer :: i ! counters
    integer :: naux_overlap ! size of auxiliary basis after reduction with respect to thr_overlap (STEP 1)
    integer :: naux_constraint
    integer :: nsing ! ! number of eigenvalues of A-matrix or B^T*B-matrix which are lower than specified threshold thr_fai_oep
    real(dp), pointer :: smat(:, :), smh(:, :), ediag_mat(:, :), ediag(:) ! variables for overlap matrix S
    real(dp), pointer :: amat(:, :), aevec(:, :), aeval(:), amat_b(:, :) ! variables for A-matrix and its construction
    real(dp), pointer :: smh_inv(:, :), ediag_mat_inv(:, :), &
                         mat_unity_test(:, :)

    ! transformation matrix to remove functions with small overlap (STEP 1)
    real(dp), pointer :: trans_mat_overlap(:, :), trans_mat_overlap_inv(:, :), &
                         trans_mat_aux(:, :)
    real(dp), pointer :: auxcharg_overlap(:) ! charge constraint vector
    real(dp), pointer :: proj_constraint(:, :) ! projection matrix for constraints only

    type(matrix_type) :: trans_mat_constraint_ ! auxiliary counterpart for optional trans_mat_constraint
    real(dp) :: auxcharg_(naux) ! auxiliary counterpart for optional auxcharg

    ! if high level of verbosity is chosen, print detailed information
    if (iverbose > 2) then
      write (iout, '(/,1x,a,es14.7)') 'thr_overlap_oep:', thr_overlap_oep
      write (iout, '(1x,a,es14.7)') 'thr_symmetry:', thr_symmetry
      write (iout, '(1x,a,es14.7)') 'thr_fai_oep:', thr_fai_oep
      write (iout, '(1x,a,l2)') 'l_vref_fa:', l_vref_fa
    end if

    ! create a pointer to the bottom of the stack used in this driver
    ibase => memory_allocate_integer(1)

    ! ###### STEP 1
    ! construct the overlap matrix S of the basis functions from basis set
    ! library with respect to the Coulomb norm
    smat => memory_allocate(naux, naux)
    call basis_sqr('OEP', 'J', smat, 0d0)

    ! generate the overlap matrix S^I (will be rewriten to smat) with a diagonal
    ! of ones and store corresponding transformation matrix (in ediag_mat)
    ediag => memory_allocate(naux)
    ediag_mat => memory_allocate(naux, naux)
    call overlap_diag_unit(smat, ediag, ediag_mat)

    ! eliminate eigenvectors with small eigenvalues in overlap matrix S^I
    ! and determine corresponding transformation matrix (in smh = \nu^{-1/2}*P)
    smh => memory_allocate(naux, naux)
    call basis_reduc(smat, smh, naux_overlap, thr_overlap_oep, thr_symmetry)

    ! determine combination of transformation matrices
    trans_mat_overlap => memory_allocate(naux, naux_overlap)
    ! mu^{-1/2} * \nu^{-1/2}*P
    call dgemm_x('N', 'N', naux, naux_overlap, naux, 1d0, &
                 ediag_mat, naux, smh(:, 1:naux_overlap), naux, &
                 0d0, trans_mat_overlap, naux)

    ! calculate transformation matrix for backfiltering at this step, if required
    if (present(trans_mat_constraint_inv) .and. (l_backfilter)) then

      smat = 0d0
      call basis_sqr('OEP', 'J', smat, 0d0)
      call overlap_diag_unit(smat, ediag, ediag_mat)

      ediag_mat_inv => memory_allocate(naux, naux)
      ediag_mat_inv = 0d0
      do i = 1, naux
        ediag_mat_inv(i, i) = 1d0/ediag_mat(i, i)
      end do

      smh_inv => memory_allocate(naux, naux)
      ! get nu^{1/2}*P^T := smh_inv
      call basis_reduc_inv(smat, smh_inv, naux_overlap, &
                           thr_overlap_oep, thr_symmetry)

      trans_mat_overlap_inv => memory_allocate(naux_overlap, naux)

      call dgemm_x('N', 'N', naux_overlap, naux, naux, 1d0, &
                   smh_inv(1:naux_overlap, :), naux_overlap, &
                   ediag_mat_inv, naux, 0d0, &
                   trans_mat_overlap_inv, naux_overlap)

      mat_unity_test => memory_allocate(naux_overlap, naux_overlap)

      call dgemm_x('N', 'N', naux_overlap, naux_overlap, naux, 1d0, &
                   trans_mat_overlap_inv, naux_overlap, &
                   trans_mat_overlap, naux, 0d0, &
                   mat_unity_test, naux_overlap)

      call matrix_unity_test(mat_unity_test, &
                             "mu^{1/2} nu^{1/2}*P^T * mu^{-1/2} nu^{-1/2}*P")

      call memory_release(mat_unity_test)

    end if

    ! ###### STEP 2
    ! calculate charge constraint vector
    auxcharg_overlap => memory_allocate(naux_overlap)
    call charge_constr_vec('OEP', trans_mat_overlap, auxcharg_, &
                           auxcharg_overlap)

    ! determine dimension of naux_constraint
    naux_constraint = naux_overlap-1

    ! determine transformation matrix which remove contributions from a constant
    ! function from the auxiliary basis set
    proj_constraint => memory_allocate(naux_overlap, naux_constraint)
    call proj(auxcharg_overlap, proj_constraint)

    ! calculate combined transformation matrix from steps 1 and 2
    call trans_mat_constraint_%init(naux, naux_constraint)
    call dgemm_x('N', 'N', naux, naux_constraint, naux_overlap, 1d0, &
                 trans_mat_overlap, naux, &
                 proj_constraint, naux_overlap, &
                 0d0, trans_mat_constraint_%mat, naux)

    ! calculate transformation matrix for backfiltering at this step, if required
    if (present(trans_mat_constraint_inv) .and. (l_backfilter)) then

      ! calculate combined transformation matrix from steps 1 and 2
      call trans_mat_constraint_inv%init(naux_constraint, naux)

      call dgemm_x('T', 'N', naux_constraint, naux, naux_overlap, 1d0, &
                   proj_constraint, naux_overlap, &
                   trans_mat_overlap_inv, naux_overlap, 0d0, &
                   trans_mat_constraint_inv%mat, naux_constraint)

      mat_unity_test => memory_allocate(naux_constraint, naux_constraint)

      call dgemm_x('N', 'N', naux_constraint, naux_constraint, naux, 1d0, &
                   trans_mat_constraint_inv%mat, naux_constraint, &
                   trans_mat_constraint_%mat, naux, 0d0, &
                   mat_unity_test, naux_constraint)

      call matrix_unity_test(mat_unity_test, &
                             "(W_1 W_2 W_3)^{-1} W_1 W_2 W_3")

      call memory_release(mat_unity_test)

    end if

    call memory_release(proj_constraint)

    ! ###### STEP 3
    ! calculate reference potentials
    call vref_ao%init(ntdg)
    call vref_aux_up%init(naux)
    if (l_vref_fa) then
      call vref_fa(trans_mat_overlap, auxcharg_, vref_ao%mat, vref_aux_up%mat)
    else
      call vref_const(trans_mat_overlap, auxcharg_overlap, &
                      vref_ao%mat, vref_aux_up%mat)
    end if

    call memory_release(trans_mat_overlap)

    if (present(int_fai)) then
      ! ###### STEP 4
      ! remove linear combinations of basis functions not coupling to products of
      ! occupied and unoccupied orbitals; construct resulting transformation matrix
      amat => memory_allocate(naux_constraint, naux_constraint)
      aevec => memory_allocate(naux_constraint, naux_constraint)
      aeval => memory_allocate(naux_constraint)

      ! construct A^{III}-matrix
      call calc_amat(trans_mat_constraint_%mat, int_fai, amat, eig, occa)
      if (present(int_fai_b)) then
        amat_b => memory_allocate(naux_constraint, naux_constraint)
        call calc_amat(trans_mat_constraint_%mat, int_fai_b, amat_b, eig, occb)
        amat = (amat+amat_b)/2.0
        call memory_release(amat_b)
      end if

      ! get eigenvectors and eigenvalues of A-matrix (or scaled A-matrix)
      call matrix_spectral(amat, aevec, aeval, naux_constraint)

      if (iverbose > 2) then
        write (iout, '(1x,a)') '----------------------------------------------'
        write (iout, '(1x,a)') '        Eigenvalues of A^{III}-matrix'
        write (iout, '(1x,a)') '----------------------------------------------'
      end if

      ! put zeros for eigenvectors of A-matrix corresponding to eigenvalues which
      ! are lower than specified threshold, count the number of such eigenvectors
      call aevec_put_zeros(aevec, aeval, thr_fai_oep, nsing)

      ! if nsing>0, reduce the size of transformation matrix ,
      ! or calculate new transformation matrix of the same size with zeros
      ! depending on specified l_fai_reduc
      if (nsing > 0) then
        if (iverbose > 0) write (iout, '(1x,a,i5)') 'nsing (A-matrix):', nsing
        call reduce_trans_mat_w_aevec(trans_mat_constraint_, aevec, naux, nsing)
      end if

      ! calculate transformation matrix for backfiltering at this step, if required
      if (present(trans_mat_constraint_inv) .and. (l_backfilter)) then
        if (nsing > 0) then

          trans_mat_aux => memory_allocate(naux_constraint, naux)
          trans_mat_aux = trans_mat_constraint_inv%mat

          call trans_mat_constraint_inv%init(naux_constraint-nsing, naux)

          call dgemm_x('N', 'N', naux_constraint-nsing, naux, &
                       naux_constraint, 1d0, &
                       transpose(aevec(:, 1+nsing:naux_constraint)), &
                       naux_constraint-nsing, &
                       trans_mat_aux, naux_constraint, 0d0, &
                       trans_mat_constraint_inv%mat, naux_constraint-nsing)

          call memory_release(trans_mat_aux)

          naux_constraint = naux_constraint-nsing

          mat_unity_test => memory_allocate(naux_constraint, naux_constraint)

          call dgemm_x('N', 'N', naux_constraint, naux_constraint, naux, 1d0, &
                       trans_mat_constraint_inv%mat, naux_constraint, &
                       trans_mat_constraint_%mat, naux, 0d0, &
                       mat_unity_test, naux_constraint)

          call matrix_unity_test(mat_unity_test, &
                                 "(W_1 W_2 W_3 W_4)^{-1} (W_1 W_2 W_3 W_4)")

          call memory_release(mat_unity_test)

        end if
      end if

      call memory_release(amat)

    end if

    if (iverbose .gt. 0) then
      write (iout, '(1x,a,i5)') &
        'Number of nonsingular components:', trans_mat_constraint_%m
    end if

    if (present(trans_mat_constraint)) then
      call trans_mat_constraint%init(trans_mat_constraint_%n, &
                                     trans_mat_constraint_%m)
      trans_mat_constraint%mat = trans_mat_constraint_%mat
    end if

    if (present(auxcharg)) auxcharg = auxcharg_

    if (allocated(trans_mat_constraint_%mat)) &
      deallocate (trans_mat_constraint_%mat)
    call memory_release(ibase)

  end subroutine basis_process

  ! Calculates Fermi-Amaldi potential
  subroutine vref_fa(trans_mat, y, vref_ao, vref_aux_up)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    use ksinv_mod, only: u_frozen, iverbose, vj_frozen, l_vh_oep
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat(:, :)
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: vref_ao(:)
    real(dp), intent(out) :: vref_aux_up(:)
    real(dp), external :: ddot_x
    integer :: n, m
    real(dp), pointer :: vc2(:), vc(:)
    real(dp) :: nelec

    n = size(trans_mat, 1)
    m = size(trans_mat, 2)

    vc2 => memory_allocate(m)
    ! v_C^{II} = W^T \cdot u
    call dgemv_x('T', n, m, 1d0, trans_mat, n, u_frozen, 1, 0d0, vc2, 1)

    vc => memory_allocate(n)
    ! v_C = W \cdot v_c^{II}
    call dgemv_x('N', n, m, 1d0, trans_mat, n, vc2, 1, 0d0, vc, 1)

    nelec = ddot_x(n, y, 1, vc, 1)
    ! v_{ref} = -1/N*v_C
    vref_aux_up = -1d0/nelec*vc
    ! Get reference potential in atomic basis
    call basis_vec_contract_a('OEP', 'ORBITAL', 'J', vref_aux_up, vref_ao, 0d0)

    if (l_vh_oep) then
      call basis_vec_contract_a('OEP', 'ORBITAL', 'J', &
                                real(nint(nelec), dp)/nelec*vc, vj_frozen, 0d0)
    end if

    if (iverbose > 2) then
      write (iout, '(a,f24.12)') 'E_Coul(OEP):', 0.5d0*ddot_x(m, vc2, 1, vc2, 1)
      write (iout, '(a,f24.12)') 'N:', ddot_x(n, y, 1, vc, 1)
    end if

    call memory_release(vc2)

  end subroutine vref_fa

  ! Eliminates eigenvectors with small eigenvalues in overlap matrix (Inversion)
  subroutine basis_reduc_inv(matrix, evec, n_out, thr, thr_sym)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    use acfd_utils, only: matrix_spectral
    use ksinv_mod, only: iverbose
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(out) :: n_out ! size of the new overlap matrix
    real(dp), intent(in) :: matrix(:, :) ! overlap matrix
    real(dp), intent(in) :: thr ! threshold for the removal of the eigenvectors
    real(dp), intent(in) :: thr_sym ! threshold for the removal of the eigenvectors
    real(dp), intent(out) :: evec(:, :) ! remaining eigenvectors after the removal

    integer :: i, cnt, n_in
    real(dp), pointer :: eval(:) ! eigenvalues
    real(dp), pointer :: r2_tmp(:, :)
    real(dp), pointer :: symvec(:), r1_tmp(:)

    real(dp), external :: ddot_x

    n_out = 0
    n_in = size(matrix, 1)

    eval => memory_allocate(n_in)
    r2_tmp => memory_allocate(n_in, n_in)
    r2_tmp = 0d0
    symvec => memory_allocate(n_in)
    symvec = 0d0
    r1_tmp => memory_allocate(n_in)
    r1_tmp = 0d0

    call matrix_spectral(matrix, evec, eval, n_in)

    if (thr_sym .gt. 0d0) then
      call auxbasis_charge(symvec, r1_tmp, n_in, 'OEP')
      where (abs(symvec) .gt. thr_sym) symvec = 1d0
      do i = 1, n_in
        if (abs(ddot_x(n_in, evec(:, i), 1, symvec, 1)) .lt. thr_sym) then
          eval(i) = 0d0
        end if
      end do
    end if

    if (iverbose > 2) then
      write (iout, *) '------------------------------------------------'
      write (iout, *) '           Eigenvalues of S^I-matrix (inv)'
      write (iout, *) '------------------------------------------------'
    end if
    do i = 1, n_in
      if (eval(i) < thr) then
        if (iverbose > 2) write (iout, '(i4,es30.15,3x,a)') i, eval(i), 'F'
      else
        if (iverbose > 2) write (iout, '(i4,es30.15,3x,a)') i, eval(i), 'T'
      end if
    end do

    cnt = 0
    do i = 1, n_in
      if (eval(i) .lt. thr) then
        cnt = cnt+1
        cycle
      end if
      n_out = n_out+1
      r2_tmp(:, n_out) = evec(:, i)
      call dscal_x(n_in, sqrt(eval(i)), r2_tmp(:, n_out), 1)
    end do
    evec = 0d0
    evec(:, 1:n_out) = r2_tmp(:, 1:n_out)
    evec = transpose(evec) ! this transpose is different from original soubroutine

    if ((cnt > 0) .and. (iverbose > 0)) then
      write (iout, '(a,i5)') 'nsing (S-matrix):', cnt
    end if

    call memory_release(eval)

  end subroutine basis_reduc_inv

  ! check whether an input matrix is unity matrix or not
  subroutine matrix_unity_test(inmatrix, designation)
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    ! size of tested matrix, is square, as is unity
    real(dp), intent(in) :: inmatrix(:, :)
    character(*) :: designation
    ! indices of dimensions of matrix
    integer :: i, j
    logical :: second_test_necessary = .false.

    do i = 1, size(inmatrix, 1)
      if (abs(inmatrix(i, i)-1d0) > 1d-13) then
        write (*, '(3a,2i4,a,ES21.14)') "Warning ", designation, &
          " is deviating from unity: Deviation of 1-Element ", &
          i, i, " is ", inmatrix(i, i)-1d0
        second_test_necessary = .True.
      end if
    end do
    if (second_test_necessary) then
      do i = 1, size(inmatrix, 1)
        do j = 1, size(inmatrix, 1)
          if (abs(inmatrix(i, j)) > 1d-13 .and. i .NE. j) then
            write (*, '(3a,2i4,a,es21.14)') "Warning ", designation, &
              " is deviating from unity: 0-Element ", i, j, " is ", &
              inmatrix(i, j)
          end if
        end do
      end do
    end if

  end subroutine matrix_unity_test

end module ksinv_basis_process
