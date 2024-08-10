!
! Module with subroutines for processing of auxiliary basis sets
!
! Some ideas behind the content can be found in
! J. Chem. Phys. 127, 054102 (2007); https://doi.org/10.1063/1.2751159
! J. Chem. Phys. 136, 134102 (2012); https://doi.org/10.1063/1.3697845
! J. Chem. Phys. 155, 054109 (2021); https://doi.org/10.1063/5.0056431
!
! Notations follow J. Chem. Phys. 155, 054109 (2021) where is possible
!
module acfd_basis_process
  implicit none
  private
  public :: basis_process_rirpa, basis_process_scexx, &
            calc_amat, aevec_put_zeros, proj, vref_const, basis_reduc, & ! for ksinv
            reduce_trans_mat_w_aevec, charge_constr_vec, overlap_diag_unit ! for ksinv
contains
  !
  ! Processing of auxiliary basis set for RIRPA and URIRPA calculations
  !
  ! Parameters:
  ! trans_mat_constraint [matrix_type, out]: transformation matrix
  ! int_fai [real(dp), in]: 3-index integrals (basis,virt,occ)
  ! int_fai_b [real(dp), in, optional]: 3-index integrals (basis,virt,occ) for
  !     second spin-channel
  !
  subroutine basis_process_rirpa(trans_mat_constraint, int_fai, int_fai_b)
    use common_tapes, only: iout
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use acfd_parameters, only: thr_overlap_ri, thr_fai_ri, iverbose
    use acfd_derived_datatypes, only: matrix_type
    use acfd_utils, only: matrix_spectral
    use iso_fortran_env, only: dp => real64
    implicit none

    type(matrix_type), intent(out) :: trans_mat_constraint
    real(dp), intent(in) :: int_fai(:, :, :)
    real(dp), intent(in), optional :: int_fai_b(:, :, :)

    integer, pointer :: ibase(:) ! pointer for stack management

    integer :: naux
    integer :: naux_overlap ! size of auxiliary basis after reduction with respect to thr_overlap
    integer :: naux_constraint ! size of auxiliary basis after all reductions

    ! number of eigenvalues of A-matrix or B^T*B-matrix which are lower
    ! than specified threshold
    integer :: nsing

    ! variables for overlap matrix S
    real(dp), pointer :: smat(:, :), smh(:, :), ediag_mat(:, :), ediag(:)

    ! variables for A-matrix, etc
    real(dp), pointer :: amat(:, :), aevec(:, :), aeval(:), &
                         amat_b(:, :)

    ! transformation matrix to remove functions with small overlap (STEP 1)
    real(dp), pointer :: trans_mat_overlap(:, :)

    ! charge constraint vector
    real(dp), pointer :: auxcharg_overlap(:)

    ! projection matrix for constraints only
    real(dp), pointer :: proj_constraint(:, :)

    ! charge constraint vector in unprocessed basis
    real(dp), pointer :: auxcharg_aux(:)

    ! if high level of verbosity is chosen, print detailed information
    if (iverbose > 2) then
      write (iout, '(1x,a,es14.7)') 'thr_overlap_ri:', thr_overlap_ri
      write (iout, '(1x,a,es14.7)') 'thr_fai_ri:', thr_fai_ri
    end if

    ! create a pointer to the bottom of the stack used in this driver
    ibase => memory_allocate_integer(1)

    naux = size(int_fai, 1)

    ! allocate auxcharg_aux
    auxcharg_aux => memory_allocate(naux)

    ! construct the overlap matrix S of the basis functions from basis set
    ! library with respect to the Coulomb norm
    smat => memory_allocate(naux, naux)
    call basis_sqr('RI', 'J', smat, 0d0)

    ! generate the overlap matrix S^I (will be rewriten to smat) with a diagonal
    ! of ones and store corresponding transformation matrix (in ediag_mat)
    ediag => memory_allocate(naux)
    ediag_mat => memory_allocate(naux, naux)
    call overlap_diag_unit(smat, ediag, ediag_mat)

    ! eliminate eigenvectors with small eigenvalues in overlap matrix S^I
    ! and determine corresponding transformation matrix (in smh = \nu^{-1/2}*P)
    smh => memory_allocate(naux, naux)
    call basis_reduc(smat, smh, naux_overlap, thr_overlap_ri, 0d0)

    ! determine combination of transformation matrices
    trans_mat_overlap => memory_allocate(naux, naux_overlap)
    ! mu^{-1/2} * \nu^{-1/2}*P
    call dgemm_x('N', 'N', naux, naux_overlap, naux, 1d0, &
                 ediag_mat, naux, smh(:, 1:naux_overlap), naux, &
                 0d0, trans_mat_overlap, naux)

    ! calculate charge constraint vector, if required
    auxcharg_overlap => memory_allocate(naux_overlap)
    call charge_constr_vec('RI', trans_mat_overlap, auxcharg_aux, &
                           auxcharg_overlap)

    ! determine dimension of naux_constraint
    naux_constraint = naux_overlap-1

    ! determine transformation matrix which remove contributions from a constant
    ! function from the auxiliary basis set
    proj_constraint => memory_allocate(naux_overlap, naux_constraint)
    call proj(auxcharg_overlap, proj_constraint)

    ! calculate combined transformation matrix from steps 1 and 2
    call trans_mat_constraint%init(naux, naux_constraint)
    call dgemm_x('N', 'N', naux, naux_constraint, naux_overlap, 1d0, &
                 trans_mat_overlap, naux, &
                 proj_constraint, naux_overlap, &
                 0d0, trans_mat_constraint%mat, naux)

    call memory_release(trans_mat_overlap)

    ! remove linear combinations of basis functions not coupling to products of
    ! occupied and unoccupied orbitals; construct resulting transformation matrix
    amat => memory_allocate(naux_constraint, naux_constraint)
    aevec => memory_allocate(naux_constraint, naux_constraint)
    aeval => memory_allocate(naux_constraint)

    ! construct A^{III}-matrix or A^{scIII}-matrix
!    call calc_amat(trans_mat_constraint%mat, int_fai, amat)
    if (present(int_fai_b)) then
      amat_b => memory_allocate(naux_constraint, naux_constraint)
!      call calc_amat(trans_mat_constraint%mat, int_fai_b, amat_b)
      amat = (amat+amat_b)/2.0
      call memory_release(amat_b)
    end if

    ! get eigenvectors and eigenvalues of A-matrix (or scaled A-matrix)
    call matrix_spectral(amat, aevec, aeval, naux_constraint)

    if (iverbose > 2) then
      write (iout, *) '------------------------------------------------'
      write (iout, *) '         Eigenvalues of A^{III}-matrix'
      write (iout, *) '------------------------------------------------'
    end if

    ! put zeros for eigenvectors of A-matrix corresponding to eigenvalues which
    ! are lower than specified threshold, count the number of such eigenvectors
    call aevec_put_zeros(aevec, aeval, thr_fai_ri, nsing)

    ! if nsing>0, reduce the size of transformation matrix ,
    if (nsing > 0) then
      if (iverbose > 0) write (iout, '(1x,a,i5)') 'nsing (A-matrix):', nsing
      call reduce_trans_mat_w_aevec(trans_mat_constraint, aevec, naux, nsing)
    end if

    call memory_release(amat)

    call memory_release(ibase)

  end subroutine basis_process_rirpa

  !
  ! Processing of auxiliary basis set for SCEXX and USCEXX calculations
  !
  ! Parameters:
  ! trans_mat_constraint [matrix_type, out]: transformation matrix
  ! int_fai [real(dp), in]: 3-index integrals (basis,virt,occ)
  ! auxcharg [real(dp), out]: charge constraint vector in unprocessed basis
  ! int_fij [real(dp), in]: 3-index integrals (basis,occ,occ)
  ! vref_ao [vector_type, inout]: reference potential in AO basis
  ! vref_aux_up [vector_type, inout]: reference potential
  ! homo [real(dp), out]: homo constrained vector in unprocessed basis
  ! vxmo [real(dp), in]: non-local exchange potential in MO basis
  ! eig [real(dp), in]: orbital eigenvalues
  ! int_fai_b [real(dp), in, optional]: 3-index integrals (basis,virt,occ)
  !     for second spin-channel
  ! int_fij_b [real(dp), in, optional]: 3-index integrals (basis,occ,occ) for
  !     second spin-channel
  !
  subroutine basis_process_scexx(trans_mat_constraint, int_fai, auxcharg, &
                                 int_fij, vref_ao, vref_aux_up, homo, vxmo, &
                                 eig, int_fai_b, int_fij_b)
    use common_tapes, only: iout
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use acfd_derived_datatypes, only: matrix_type, vector_type
    use acfd_parameters, only: thr_overlap_oep, thr_fai_oep, thr_symmetry, &
                               l_vref_fa, l_homo, l_vh_oep, iverbose
    use acfd_utils, only: matrix_spectral
    use acfd_wf, only: occa, occb, int_all
    use iso_fortran_env, only: dp => real64
    implicit none

    type(vector_type), intent(inout) :: vref_ao, vref_aux_up
    real(dp), intent(in) :: int_fai(:, :, :)
    real(dp), intent(in) :: int_fij(:, :, :)
    real(dp), intent(out) :: auxcharg(:)

    real(dp), intent(out) :: homo(:)
    real(dp), intent(in) :: eig(:)
    real(dp), intent(in) :: vxmo(size(eig), size(eig))
    type(matrix_type), intent(out) :: trans_mat_constraint
    real(dp), intent(in), optional :: int_fai_b(:, :, :)
    real(dp), intent(in), optional :: int_fij_b(:, :, :)

    integer, pointer :: ibase(:) ! pointer for stack management

    integer :: naux
    integer :: ntg ! nvo+noo
    integer :: ntdg ! size of triagonal matrix of (ntg,ntg)
    integer :: naux_overlap ! size of auxiliary basis after reduction with respect to thr_overlap (STEP 1)
    integer :: naux_constraint ! size of auxiliary basis after reduction and constraints

    ! number of eigenvalues of A-matrix or B^T*B-matrix which are lower
    ! than specified threshold (thr_fai_oep)
    integer :: nsing

    ! variables for overlap matrix S
    real(dp), pointer :: smat(:, :), smh(:, :), ediag_mat(:, :), ediag(:)

    ! variables for A-matrix, B-matrix, B^T*B-matrix, etc
    real(dp), pointer :: amat(:, :), aevec(:, :), aeval(:), amat_b(:, :)

    ! transformation matrix to remove functions with small overlap (STEP 1)
    real(dp), pointer :: trans_mat_overlap(:, :)

    ! charge constraint vector
    real(dp), pointer :: auxcharg_overlap(:)

    ! HOMO constraint vector
    real(dp), pointer :: homo_overlap(:)

    ! projection matrix for constraints only
    real(dp), pointer :: proj_constraint(:, :)

    ! if high level of verbosity is chosen, print detailed information
    if (iverbose > 2) then
      write (iout, '(1x,a,es14.7)') 'thr_overlap_oep:', thr_overlap_oep
      write (iout, '(1x,a,es14.7)') 'thr_symmetry:', thr_symmetry
      write (iout, '(1x,a,es14.7)') 'thr_fai_oep:', thr_fai_oep
      write (iout, '(1x,a,l2)') 'l_vref_fa:', l_vref_fa
      write (iout, '(1x,a,l2)') 'l_homo:', l_homo
    end if

    ! create a pointer to the bottom of the stack used in this driver
    ibase => memory_allocate_integer(1)

    naux = size(int_fai, 1)

    ntg = size(int_fai, 2)!+size(int_fai, 3) ! determine number of orbital basis functions
    ntdg = ntg*(ntg+1)/2 ! determine size of triagonal matrix of (ntg,ntg)

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

    ! ###### STEP 2
    ! calculate charge constraint vector
    auxcharg_overlap => memory_allocate(naux_overlap)
    call charge_constr_vec('OEP', trans_mat_overlap, auxcharg, auxcharg_overlap)

    ! calculate HOMO constraint vector, if required
    if (l_homo) then
      homo_overlap => memory_allocate(naux_overlap)
      call homo_constr_vec(trans_mat_overlap, homo, vxmo, int_fij, eig, &
                           homo_overlap)
    end if

    ! determine dimension naux_constraint
    if (l_homo) then
      naux_constraint = naux_overlap-2
    else
      naux_constraint = naux_overlap-1
    end if

    ! determine transformation matrix which remove contributions from a constant
    ! function from the auxiliary basis set
    proj_constraint => memory_allocate(naux_overlap, naux_constraint)
    if (l_homo) then
      call proj_charge_and_homo(auxcharg_overlap, homo_overlap, proj_constraint)
    else
      call proj(auxcharg_overlap, proj_constraint)
    end if

    ! calculate combined transformation matrix from steps 1 and 2
    call trans_mat_constraint%init(naux, naux_constraint)
    call dgemm_x('N', 'N', naux, naux_constraint, naux_overlap, 1d0, &
                 trans_mat_overlap, naux, &
                 proj_constraint, naux_overlap, &
                 0d0, trans_mat_constraint%mat, naux)

    ! Constract Hartree potential in OEP basis, if required
    if (l_vh_oep) then
      if (.not. (present(int_fij_b))) then
        call vh_oep(trans_mat_overlap, int_fij, auxcharg)
      else
        call vh_oep(trans_mat_overlap, int_all, auxcharg, int_all)
      end if
    end if

    ! ###### STEP 3
    ! calculate reference potentials
    call vref_ao%init(ntdg)
    call vref_aux_up%init(naux)

    if (l_homo) then
      if (l_vref_fa) then
        if (present(int_fij_b)) then
          ! spin case
          call vref_fa_and_homo_spin(trans_mat_overlap, homo_overlap, &
                                     auxcharg_overlap, int_fij, int_fij_b, &
                                     vxmo(size(int_fai, 3), size(int_fai, 3)), & ! size incorrect now
                                     auxcharg, vref_ao%mat, vref_aux_up%mat)
        else
          call vref_fa_and_homo(trans_mat_overlap, homo_overlap, &
                                auxcharg_overlap, int_fij, &
                                vxmo(size(int_fai, 3), size(int_fai, 3)), & ! size incorrect now
                                auxcharg, vref_ao%mat, vref_aux_up%mat)
        end if
      else
        call vref_const_and_homo(trans_mat_overlap, &
                                 homo_overlap, auxcharg_overlap, &
                                 vxmo(size(int_fai, 3), size(int_fai, 3)), &
                                 auxcharg, vref_ao%mat, vref_aux_up%mat)
      end if
    else
      if (l_vref_fa) then
        if (present(int_fij_b)) then
          ! spin case
          call vref_fa_spin(trans_mat_overlap, int_fij, int_fij_b, &
                            auxcharg, vref_ao%mat, vref_aux_up%mat)
        else
          call vref_fa(trans_mat_overlap, int_fij, auxcharg, &
                       vref_ao%mat, vref_aux_up%mat)
        end if
      else
        call vref_const(trans_mat_overlap, auxcharg_overlap, vref_ao%mat, &
                        vref_aux_up%mat)
      end if
    end if

    call memory_release(trans_mat_overlap)

    ! ###### STEP 4
    ! remove linear combinations of basis functions not coupling to products of
    ! occupied and unoccupied orbitals; construct resulting transformation matrix
    amat => memory_allocate(naux_constraint, naux_constraint)
    aevec => memory_allocate(naux_constraint, naux_constraint)
    aeval => memory_allocate(naux_constraint)

    call calc_amat(trans_mat_constraint%mat, int_fai, amat, eig, occa)
    if (present(int_fai_b)) then
      amat_b => memory_allocate(naux_constraint, naux_constraint)
      call calc_amat(trans_mat_constraint%mat, int_fai_b, amat_b, eig, occb)
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
    if (nsing > 0) then
      if (iverbose > 1) write (iout, '(1x, a,i5)') 'nsing (A-matrix):', nsing
      call reduce_trans_mat_w_aevec(trans_mat_constraint, aevec, naux, nsing)
    end if

    call memory_release(amat)

    if (iverbose .gt. 1) then
      write (iout, '(1x,a,i5)') &
        'Number of nonsingular components:', trans_mat_constraint%m
    end if

    call memory_release(ibase)

  end subroutine basis_process_scexx

  ! Calculates A^{III}-matrix as $A^{III} = {D^{III}}^T \cdot D^{III}$
  ! with D^{III} = D \cdot W^{III}
  subroutine calc_amat(trans_mat, int_fkl, amat, eig, occ)
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat(:, :), int_fkl(:, :, :), eig(:), occ(:)
    real(dp), intent(out) :: amat(:, :)
    integer :: n, m, ntg, i, j
    real(dp), pointer :: dmat(:, :, :), dmat1(:, :, :)

    n = size(trans_mat, 1)
    m = size(trans_mat, 2)
    ntg = size(int_fkl, 2)

    dmat => memory_allocate(m, ntg, ntg)
    dmat1 => memory_allocate(m, ntg, ntg)

    ! D^{III} = D \cdot W^{III}
    call dgemm_x('T', 'N', m, ntg*ntg, n, 1d0, trans_mat, n, int_fkl, n, 0d0, &
                 dmat, m)

    ! scale D^{III}
    do i = 1, ntg
      do j = 1, ntg
        if (i /= j) then
          dmat1(:, i, j) = occ(i)*dmat(:, i, j)*sign(1d0, eig(j)-eig(i))
        else
          dmat1(:, i, j) = 0d0
        end if
      end do
    end do

    ! A^{III} = {D^{III}}^T \cdot D^{III}
    call dgemm_x('N', 'T', m, m, ntg*ntg, 1d0, dmat1, m, dmat, m, 0d0, amat, m)
    call memory_release(dmat)

  end subroutine calc_amat


  ! Determines the number of eigenvectors of A^{III}-matrix which correspond to
  ! eigenvalues lower than specified threshold and puts zeros in corresponding
  ! columns
  subroutine aevec_put_zeros(evec, eval, thr, nsing)
    use common_tapes, only: iout
    use acfd_parameters, only: iverbose
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(out) :: nsing
    real(dp), intent(in) :: thr, eval(:)
    real(dp), intent(inout) :: evec(:, :)
    integer :: i

    nsing = 0
    do i = 1, size(eval)
      if (eval(i) < thr) then
        nsing = nsing+1
        evec(:, i) = 0d0
        if (iverbose > 2) write (iout, '(1x,i4,es30.15,3x,a)') i, eval(i), 'F'
      else
        if (iverbose > 2) write (iout, '(1x,i4,es30.15,3x,a)') i, eval(i), 'T'
      end if
    end do

  end subroutine aevec_put_zeros

  ! Calculates new transformation matrix W as product of previous transformation
  ! matrix W^{III} and eigenvectors of A^{III}-matrix
  subroutine reduce_trans_mat_w_aevec(trans_mat, aevec, n, nsing)
    use memory, only: memory_allocate, memory_release
    use acfd_derived_datatypes, only: matrix_type
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: n, nsing
    real(dp), intent(in) :: aevec(:, :)
    type(matrix_type), intent(inout) :: trans_mat
    real(dp), pointer :: trans_mat_aux(:, :)
    integer :: m

    m = size(aevec, 1)

    trans_mat_aux => memory_allocate(n, m)
    trans_mat_aux = trans_mat%mat

    if (nsing > 0) call trans_mat%init(n, m-nsing)

    ! W = W^{III} \cdot U
    call dgemm_x('N', 'N', n, m-nsing, m, 1d0, trans_mat_aux, n, &
                 aevec(:, 1+nsing:m), m, 0d0, trans_mat%mat, n)

    call memory_release(trans_mat_aux)

  end subroutine reduce_trans_mat_w_aevec

  ! Calculates charge constraint vector in original and transformed basis
  subroutine charge_constr_vec(basis, trans_mat, y, y_transformed)
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    character(*), intent(in) :: basis
    real(dp), intent(in) :: trans_mat(:, :)
    real(dp), intent(out) :: y(:)
    real(dp), intent(out) :: y_transformed(:)
    real(dp), pointer :: vec_aux(:)

    vec_aux => memory_allocate(size(trans_mat, 1))
    ! Get charge constraint vector y in original basis
    call auxbasis_charge(y, vec_aux, size(trans_mat, 1), basis)
    call memory_release(vec_aux)
    ! Transform charge constraint vector into current basis
    ! y_{transformed} = W^T \cdot y
    call dgemv_x('T', size(trans_mat, 1), size(trans_mat, 2), 1d0, trans_mat, &
                 size(trans_mat, 1), y, 1, 0d0, y_transformed, 1)

  end subroutine charge_constr_vec

  ! generate overlap matrix with a diagonal of ones
  subroutine overlap_diag_unit(s, diag, diag_mat)
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(inout) :: s(:, :) ! overlap matrix
    real(dp), intent(out) :: diag(:) ! diagonal entries of diag_mat
    real(dp), intent(out) :: diag_mat(:, :) ! diagonal matrix used to transform the overlap matrix

    integer :: i

    diag_mat = 0d0
    do i = 1, size(diag)
      diag(i) = 1d0/sqrt(s(i, i))
    end do
    do i = 1, size(diag)
      diag_mat(i, i) = 1d0/sqrt(s(i, i))
    end do

    call matrix_times_diag(s, diag, size(diag), 'r')
    call matrix_times_diag(s, diag, size(diag), 'l')

  end subroutine overlap_diag_unit

  ! Eliminates eigenvectors with small eigenvalues in overlap matrix
  subroutine basis_reduc(matrix, evec, n_out, thr, thr_sym)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    use acfd_parameters, only: iverbose
    use acfd_utils, only: matrix_spectral
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
    real(dp), pointer  :: symvec(:), r1_tmp(:)

    real(dp), external :: ddot_x

    n_in = size(matrix, 1)

    eval => memory_allocate(n_in)
    r2_tmp => memory_allocate(n_in, n_in)
    symvec => memory_allocate(n_in)
    r1_tmp => memory_allocate(n_in)

    symvec = 0d0
    r1_tmp = 0d0
    call matrix_spectral(matrix, evec, eval, n_in)

    r2_tmp = 0d0
    n_out = 0

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
      write (iout, *) '           Eigenvalues of S^I-matrix'
      write (iout, *) '------------------------------------------------'
    end if
    do i = 1, n_in
      if (eval(i) < thr) then
        if (iverbose > 2) write (iout, '(1x,i4,es30.15,3x,a)') i, eval(i), 'F'
      else
        if (iverbose > 2) write (iout, '(1x,i4,es30.15,3x,a)') i, eval(i), 'T'
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
      call dscal_x(n_in, 1d0/sqrt(eval(i)), r2_tmp(:, n_out), 1)
    end do
    evec = 0d0
    evec(:, 1:n_out) = r2_tmp(:, 1:n_out)

    if ((cnt > 0) .and. (iverbose > 0)) then
      write (iout, '(1x,a,i5)') 'nsing (S-matrix):', cnt
    end if

    call memory_release(eval)

  end subroutine basis_reduc

  ! Determines projection/transformation matrix which eliminates contributions
  ! from a constant function from the auxiliary basis set for only charge case
  subroutine proj(vector, proj_mat)
    use memory, only: memory_allocate, memory_release
    use acfd_utils, only: matrix_spectral, identity_mat
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp), intent(in) :: vector(:) ! constraint vector
    real(dp), intent(out) :: proj_mat(:, :) ! projection matrix

    real(dp), pointer :: eval(:)
    real(dp), pointer :: evec(:, :)
    real(dp) :: charge_norm
    real(dp) :: min_val
    integer :: min_pos

    real(dp), pointer :: r2_tmp1(:, :), r2_tmp2(:, :)
    integer :: iaux, n
    integer :: naux
    real(dp), external :: ddot_x

    naux = size(vector)

    proj_mat = 0d0
    r2_tmp1 => memory_allocate(naux, naux)
    call identity_mat(r2_tmp1, naux)
    charge_norm = ddot_x(naux, vector, 1, vector, 1)
    do iaux = 1, naux
      r2_tmp1(:, iaux) = r2_tmp1(:, iaux)-vector(iaux)*vector(:)/charge_norm
    end do

    r2_tmp2 => memory_allocate(naux, naux)
    r2_tmp2 = 0d0
    evec => memory_allocate(naux, naux)
    eval => memory_allocate(naux)
    call dgemm_x('T', 'N', naux, naux, naux, 1d0, r2_tmp1, naux, &
                 r2_tmp1, naux, 0d0, r2_tmp2, naux)
    call matrix_spectral(r2_tmp2, evec, eval, naux)

    min_pos = 0
    min_val = 1d60
    do iaux = 1, naux
      if (eval(iaux) .lt. min_val) then
        min_pos = iaux
        min_val = eval(iaux)
      end if
    end do
    if (min_pos == 0) call error('min_pos=0', 'proj')
    n = 0
    do iaux = 1, naux
      if (iaux == min_pos) cycle
      n = n+1
      proj_mat(:, n) = evec(:, iaux)
      call dscal_x(naux, 1d0/sqrt(eval(iaux)), proj_mat(:, n), 1)
    end do

    call memory_release(r2_tmp1)

  end subroutine proj

  ! Calculates reference potential for only charge constraint case when
  ! Fermi-Amaldi option is not active
  subroutine vref_const(trans_mat_overlap, auxcharg_overlap, &
                        vref_ao, vref_aux_up)
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat_overlap(:, :)
    real(dp), intent(inout) :: auxcharg_overlap(:)
    real(dp), intent(out) :: vref_ao(:)
    real(dp), intent(out) :: vref_aux_up(:)
    integer :: naux, naux_overlap
    real(dp), external :: ddot_x
    real(dp) :: r0_tmp1
    real(dp), pointer :: vref_aux(:), r1_tmp1(:)

    naux = size(trans_mat_overlap, 1)
    naux_overlap = size(trans_mat_overlap, 2)

    r1_tmp1 => memory_allocate(naux)
    r1_tmp1 = 0d0
    vref_aux => memory_allocate(naux_overlap)
    vref_aux = 0d0
    r0_tmp1 = ddot_x(naux_overlap, auxcharg_overlap, 1, auxcharg_overlap, 1)
    call dscal_x(naux_overlap, -1d0/r0_tmp1, auxcharg_overlap, 1)
    vref_aux = auxcharg_overlap
    call dgemv_x('N', naux, naux_overlap, 1d0, &
                 trans_mat_overlap, naux, &
                 vref_aux, 1, 0d0, r1_tmp1, 1)
    call basis_vec_contract_a('OEP', 'ORBITAL', 'J', r1_tmp1, vref_ao, 0d0)
    vref_aux_up = r1_tmp1
    call memory_release(r1_tmp1)
  end subroutine vref_const

  ! Calculates Fermi-Amaldi potential
  subroutine vref_fa(trans_mat, int_fij, y, vref_ao, vref_aux_up)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    use acfd_parameters, only: iverbose
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat(:, :)
    real(dp), intent(in) :: int_fij(:, :, :)
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: vref_ao(:)
    real(dp), intent(out) :: vref_aux_up(:)
    real(dp), external :: ddot_x
    integer :: i, j, n, m
    real(dp), pointer :: u(:), vc2(:), vc(:)
    real(dp) :: nelec

    n = size(trans_mat, 1)
    m = size(trans_mat, 2)

    u => memory_allocate(n)
    u = 0d0

    do i = 1, n
      do j = 1, size(int_fij, 2)
        u(i) = u(i)+2d0*int_fij(i, j, j)
      end do
    end do

    vc2 => memory_allocate(m)
    ! v_C^{II} = W^T \cdot u
    call dgemv_x('T', n, m, 1d0, trans_mat, n, u, 1, 0d0, vc2, 1)

    vc => memory_allocate(n)
    ! v_C = W \cdot v_c^{II}
    call dgemv_x('N', n, m, 1d0, trans_mat, n, vc2, 1, 0d0, vc, 1)

    nelec = ddot_x(n, y, 1, vc, 1)
    ! v_{ref} = -1/N*v_C
    vref_aux_up = -1d0/nelec*vc
    ! Get reference potential in atomic basis
    call basis_vec_contract_a('OEP', 'ORBITAL', 'J', vref_aux_up, vref_ao, 0d0)

    if (iverbose > 2) then
      write (iout, '(1x,a,f24.12)') 'E_Coul (OEP):', &
        0.5d0*ddot_x(m, vc2, 1, vc2, 1)
      write (iout, '(1x,a,f24.12)') 'N:', ddot_x(n, y, 1, vc, 1)
    end if

    call memory_release(u)

  end subroutine vref_fa

  ! Calculates Fermi-Amaldi potential for spin case
  subroutine vref_fa_spin(trans_mat, int_fij, int_fij_b, y, vref_ao, &
                          vref_aux_up)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    use acfd_parameters, only: iverbose, l_vref_fa_sameab, l_spin_sym
    use acfd_wf, only: int_all, occa, occb
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat(:, :)
    real(dp), intent(in) :: int_fij(:, :, :), int_fij_b(:, :, :)
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: vref_ao(:)
    real(dp), intent(out) :: vref_aux_up(:)
    real(dp), external :: ddot_x
    integer :: i, j, n, m
    real(dp), pointer :: u(:), vc2(:), vc(:)
    real(dp) :: nelec

    n = size(trans_mat, 1)
    m = size(trans_mat, 2)

    u => memory_allocate(n)
    u = 0d0

    if ((l_vref_fa_sameab) .or. (l_spin_sym)) then
      do i = 1, n
        do j = 1, size(int_all, 2)
          u(i) = u(i)+occa(j)*int_all(i, j, j)
        end do
      end do
      do i = 1, n
        do j = 1, size(int_all, 2)
          u(i) = u(i)+occb(j)*int_all(i, j, j)
        end do
      end do
      u = 0.5d0*u
    else
      do i = 1, n
        do j = 1, size(int_fij, 2)
          u(i) = u(i)+int_fij(i, j, j)
        end do
      end do
    end if

    vc2 => memory_allocate(m)
    ! v_C^{II} = W^T \cdot u
    call dgemv_x('T', n, m, 1d0, trans_mat, n, u, 1, 0d0, vc2, 1)

    vc => memory_allocate(n)
    ! v_C = W \cdot v_c^{II}
    call dgemv_x('N', n, m, 1d0, trans_mat, n, vc2, 1, 0d0, vc, 1)

    nelec = ddot_x(n, y, 1, vc, 1)
    ! v_{ref} = -1/N*v_C
    vref_aux_up = -1d0/nelec*vc
    ! Get reference potential in atomic basis
    call basis_vec_contract_a('OEP', 'ORBITAL', 'J', vref_aux_up, vref_ao, 0d0)

    if (iverbose > 2) then
      write (iout, '(1x,a,f24.12)') 'E_Coul (OEP):', &
        0.5d0*ddot_x(m, vc2, 1, vc2, 1)
      write (iout, '(1x,a,f24.12)') 'N:', ddot_x(n, y, 1, vc, 1)
    end if

    call memory_release(u)

  end subroutine vref_fa_spin

  ! Calculates Hartree potnetial in OEP basis
  subroutine vh_oep(trans_mat, int_fij, y, int_fij_b)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    use acfd_parameters, only: iverbose
    use acfd_pot, only: vj
    use acfd_wf, only: occa, occb
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat(:, :)
    real(dp), intent(in) :: int_fij(:, :, :)
    real(dp), intent(in), optional :: int_fij_b(:, :, :)
    real(dp), intent(in) :: y(:)
    real(dp), external :: ddot_x
    integer :: i, j, noo, n, m, noo_b, ntg
    real(dp), pointer :: u(:), vc2(:), vc(:)
    real(dp) :: nelec

    n = size(trans_mat, 1)
    m = size(trans_mat, 2)
    noo = nint(sum(occa))
    noo_b = nint(sum(occb))
    ntg = size(int_fij, 2)
!    if (present(int_fij_b)) noo_b = size(int_fij_b, 2)

    u => memory_allocate(n)
    u = 0d0

    if (.not. (present(int_fij_b))) then
      do i = 1, n
        do j = 1, ntg
          u(i) = u(i)+2d0*int_fij(i, j, j)
        end do
      end do
    else
      do i = 1, n
        do j = 1, ntg
          u(i) = u(i)+occa(j)*int_fij(i, j, j)
        end do
      end do
      do i = 1, n
        do j = 1, ntg
          u(i) = u(i)+occb(j)*int_fij_b(i, j, j)
        end do
      end do
    end if

    vc2 => memory_allocate(m)
    ! v_C^{II} = W^T \cdot u
    call dgemv_x('T', n, m, 1d0, trans_mat, n, u, 1, 0d0, vc2, 1)

    vc => memory_allocate(n)
    ! v_C = W \cdot v_c^{II}
    call dgemv_x('N', n, m, 1d0, trans_mat, n, vc2, 1, 0d0, vc, 1)

    nelec = ddot_x(n, y, 1, vc, 1)

    if (.not. (present(int_fij_b))) then
      call basis_vec_contract_a('OEP', 'ORBITAL', 'J', &
                                2d0*real(noo, dp)/nelec*vc, vj, 0d0)
    else
      call basis_vec_contract_a('OEP', 'ORBITAL', 'J', &
                                real(noo+noo_b, dp)/nelec*vc, vj, 0d0)
    end if

    if (iverbose > 1) then
      write (iout, '(1x,a,f24.12)') 'E_Coul (OEP):', &
        0.5d0*ddot_x(m, vc2, 1, vc2, 1)
      write (iout, '(1x,a,f24.12)') 'N:', ddot_x(n, y, 1, vc, 1)
    end if

    call memory_release(u)

  end subroutine vh_oep

  ! Calculates HOMO constraint vector
  subroutine homo_constr_vec(trans_mat_overlap, homo, vxmo, int_fij, eig, &
                             homo_overlap)
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat_overlap(:, :)
    real(dp), intent(in) :: int_fij(:, :, :)
    real(dp), intent(in) :: vxmo(:, :)
    real(dp), intent(in) :: eig(:)
    real(dp), intent(inout) :: homo(:)
    real(dp), intent(out) :: homo_overlap(:)
    real(dp) :: r0_tmp1
    real(dp), pointer :: r1_tmp1(:)

    integer :: naux, naux_overlap, ntg, noo

    naux = size(trans_mat_overlap, 1)
    naux_overlap = size(trans_mat_overlap, 2)
    ntg = size(eig)
    noo = size(int_fij, 2)

    homo_overlap = 0d0

    r1_tmp1 => memory_allocate(naux)
    call orbital_constraint(homo, r1_tmp1, r0_tmp1, naux, vxmo, eig, ntg, &
                                 noo, noo, int_fij)
    call memory_release(r1_tmp1)
    call dgemv_x('T', naux, naux_overlap, 1d0, trans_mat_overlap, naux, homo, &
                 1, 0d0, homo_overlap, 1)
  end subroutine homo_constr_vec

  ! Calculates HOMO constraint vector in original (untransformed) basis set
  subroutine orbital_constraint(vector, potential, val, naux, vxnl, eig, &
                                ntg, no, iorb, int_fij)
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: thr_eig = 1d-6
    integer, intent(in) :: naux ! number of auxiliary basis functions
    integer, intent(in) :: no ! number of occupied orbitals
    integer, intent(in) :: ntg ! number of orbitals
    integer, intent(in) :: iorb ! orbital used for constaint
    real(dp), intent(in) :: vxnl(ntg, ntg) ! non-local exchange potential
    real(dp), intent(in) :: eig(ntg) ! orbital eigenvalues
    real(dp), intent(in) :: int_fij(naux, no, no) ! 3-index integrals (auxbasis,occ,occ)
    real(dp), intent(out) :: vector(naux), potential(naux), val

    integer :: iaux, ndgen, i
    real(dp) :: zz, fac, diff
    integer :: idorb(10), j

    integer, external :: idgeneracy_homo
    real(dp), external :: ddot_x

    val = vxnl(iorb, iorb)
    ndgen = idgeneracy_homo(eig, iorb, thr_eig)
    call fzero(vector, naux)

    if (ndgen == 1) then
      do iaux = 1, naux
        vector(iaux) = vector(iaux)+int_fij(iaux, iorb, iorb)
      end do
    else
      call izero(idorb, 10)
      idorb(1) = iorb
      j = 2
      do i = 1, ndgen
        if ((iorb-i) .gt. 0) then
          diff = abs(eig(iorb-i)-eig(iorb))
          if (diff .lt. thr_eig) then
            idorb(j) = iorb-i
            j = j+1
          end if
        end if
        if ((iorb+i) .gt. no) cycle ! no constraints for virtuals
        if ((iorb+i) .lt. ntg) then
          diff = abs(eig(iorb+i)-eig(iorb))
          if (diff .lt. thr_eig) then
            idorb(j) = iorb+i
            j = j+1
          end if
        end if
      end do
      if ((j-1) /= ndgen) then
        call error('degeneracy check', 'orbital_constraint')
      end if
      do iaux = 1, naux
        do i = 1, ndgen
          j = idorb(i)
          vector(iaux) = vector(iaux)+int_fij(iaux, j, j)
        end do
      end do
      call dscal_x(naux, 1d0/dble(ndgen), vector, 1)
    end if

    zz = ddot_x(naux, vector, 1, vector, 1)
    fac = val/zz
    call fmove(vector, potential, naux)
    call dscal_x(naux, fac, potential, 1)

  end subroutine orbital_constraint

  ! Determines projection/transformation matrix which eliminates contributions
  ! from a constant function from the auxiliary basis set for charge and HOMO case
  subroutine proj_charge_and_homo(auxcharg_overlap, homo_overlap, &
                                  proj_constraint)
    use acfd_utils, only: identity_mat
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp), intent(in) :: auxcharg_overlap(:)
    real(dp), intent(in) :: homo_overlap(:)
    real(dp), intent(out) :: proj_constraint(:, :)

    integer :: info, i, iaux, naux_overlap, naux_constraint

    real(dp), pointer :: constraint_vector(:, :)
    real(dp), pointer :: overlap_cv(:, :) ! overlap matrix of constraint_vector
    real(dp), pointer :: evec_overlap_cv(:, :) ! eigenvectors of overlap_cv
    real(dp), pointer :: eval_overlap_cv(:) ! eigenvalues of overlap_cv
    real(dp), pointer :: eval_ocv_sqrtinv(:, :) ! inverse square root of eval_ocv
    real(dp), pointer :: orth_constraint(:, :) ! matrix of orthonormalized constraints
    real(dp), pointer :: proj_matrix(:, :) ! intermediate projection matrix
    real(dp), pointer :: evec_oov(:, :) ! eigenvectors of orth_overlap
    real(dp), pointer :: eval_oov(:) ! eigenvalues of orth_overlap
    real(dp), pointer :: orth_overlap(:, :) ! overlap matrix of auxiliary basis made orthogonal to the constraints

    real(dp), pointer :: r1_tmp1(:), r1_tmp2(:)

    integer, parameter :: ncon = 2
    real(dp), parameter :: zero = 1d-12

    naux_overlap = size(auxcharg_overlap)
    naux_constraint = size(proj_constraint, 2)

    ! collect constraint vectors
    constraint_vector => memory_allocate(naux_overlap, ncon)
    constraint_vector(:, 1) = auxcharg_overlap
    constraint_vector(:, 2) = homo_overlap

    ! calculate overlap of constraint vector
    overlap_cv => memory_allocate(ncon, ncon)
    call dgemm_x('T', 'N', ncon, ncon, naux_overlap, 1d0, &
                 constraint_vector, naux_overlap, &
                 constraint_vector, naux_overlap, &
                 0d0, overlap_cv, ncon)

    ! diagonalize overlap of constraint vector
    evec_overlap_cv => memory_allocate(ncon, ncon)
    eval_overlap_cv => memory_allocate(ncon)
    r1_tmp1 => memory_allocate(6*ncon)
    evec_overlap_cv = overlap_cv
    info = 0
    call dsyev_x('V', 'U', ncon, evec_overlap_cv, ncon, eval_overlap_cv, &
                 r1_tmp1, ncon*6, info)
    if (info /= 0) &
      call error('error in dsyev_x', 'proj_constraint_charge_and_homo')
    call memory_release(r1_tmp1)
    eval_ocv_sqrtinv => memory_allocate(ncon, ncon)
    eval_ocv_sqrtinv = 0d0
    do i = 1, ncon
      eval_ocv_sqrtinv(i, i) = 1d0/sqrt(eval_overlap_cv(i))
    end do

    ! construction of orthonormalized constraints
    orth_constraint => memory_allocate(naux_overlap, ncon)
    r1_tmp1 => memory_allocate(ncon*ncon)
    call dgemm_x('N', 'N', ncon, ncon, ncon, 1d0, evec_overlap_cv, ncon, &
                 eval_ocv_sqrtinv, ncon, 0d0, r1_tmp1, ncon)
    call dgemm_x('N', 'N', naux_overlap, ncon, ncon, 1d0, &
                 constraint_vector, naux_overlap, &
                 r1_tmp1, ncon, 0d0, orth_constraint, naux_overlap)
    call memory_release(r1_tmp1)

    ! intermediate projection matrix
    proj_matrix => memory_allocate(naux_overlap, naux_overlap)
    call identity_mat(proj_matrix, naux_overlap)
    do iaux = 1, naux_overlap
      do i = 1, ncon
        call daxpy_x(naux_overlap, -orth_constraint(iaux, i), &
                     orth_constraint(:, i), 1, proj_matrix(:, iaux), 1)
      end do
    end do

    ! construct and diagonalize overlap matrix for auxiliary functions
    ! made orthogonal to the constraints
    orth_overlap => memory_allocate(naux_overlap, naux_overlap)
    evec_oov => memory_allocate(naux_overlap, naux_overlap)
    eval_oov => memory_allocate(naux_overlap)
    r1_tmp1 => memory_allocate(naux_overlap**2)
    call dgemm_x('T', 'N', naux_overlap, naux_overlap, naux_overlap, 1d0, &
                 proj_matrix, naux_overlap, &
                 proj_matrix, naux_overlap, &
                 0d0, orth_overlap, naux_overlap)
    evec_oov = orth_overlap
    r1_tmp2 => memory_allocate(naux_overlap**2)
    call dgesvd_x('A', 'A', naux_overlap, naux_overlap, orth_overlap, &
                  naux_overlap, eval_oov, evec_oov, naux_overlap, r1_tmp1, &
                  naux_overlap, r1_tmp2, naux_overlap**2, info)
    if (info /= 0) &
      call error('problem in dgesvd', 'proj_constraint_charge_and_homo')
    call memory_release(r1_tmp1)

    ! remove columns of evec_oov with zero eigenvalues
    i = 0
    proj_constraint = 0d0
    do iaux = 1, naux_overlap
      if (abs(eval_oov(iaux)) .gt. zero) then
        i = i+1
        proj_constraint(:, i) = evec_oov(:, iaux)
      end if
    end do
    if (i /= naux_constraint) &
      call error('projection matrix has wrong number of zero eigenvalues', &
                 'proj_constraint_charge_and_homo')

    call memory_release(constraint_vector)

  end subroutine proj_charge_and_homo

  ! Calculates reference potential for the case of charge and homo condition
  ! without Fermi-Amaldi treatment of charge part
  subroutine vref_const_and_homo(trans_mat_overlap, homo_overlap, &
                                 auxcharg_overlap, vxmonoonoo, auxcharg, &
                                 vref_ao, vref_aux_up)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_parameters, only: iverbose
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat_overlap(:, :)
    real(dp), intent(inout) :: auxcharg_overlap(:)
    real(dp), intent(inout) :: homo_overlap(:)
    real(dp), intent(in) :: vxmonoonoo
    real(dp), intent(inout) :: auxcharg(:)
    real(dp), intent(out) :: vref_ao(:)
    real(dp), intent(out) :: vref_aux_up(:)
    real(dp), external :: ddot_x
    integer :: info, naux, naux_overlap
    integer, pointer :: ipiv(:)
    real(dp), pointer :: vref_aux(:), r1_tmp1(:), r1_tmp2(:), r2_tmp1(:, :)

    naux = size(trans_mat_overlap, 1)
    naux_overlap = size(trans_mat_overlap, 2)

    r1_tmp1 => memory_allocate(naux)
    r1_tmp1 = 0d0
    vref_aux => memory_allocate(naux_overlap)
    vref_aux = 0d0
    r2_tmp1 => memory_allocate(2, 2)
    r2_tmp1 = 0d0
    r1_tmp2 => memory_allocate(2)
    r1_tmp2 = 0d0
    r2_tmp1(1, 1) = ddot_x(naux_overlap, auxcharg_overlap, 1, &
                           auxcharg_overlap, 1)
    r2_tmp1(1, 2) = ddot_x(naux_overlap, auxcharg_overlap, 1, homo_overlap, 1)
    r2_tmp1(2, 1) = ddot_x(naux_overlap, homo_overlap, 1, auxcharg_overlap, 1)
    r2_tmp1(2, 2) = ddot_x(naux_overlap, homo_overlap, 1, homo_overlap, 1)
    r1_tmp2(1) = -1d0
    r1_tmp2(2) = vxmonoonoo
    info = 0
    ipiv => memory_allocate_integer(2)
    call dgesv_x(2, 1, r2_tmp1, 2, ipiv, r1_tmp2, 2, info)
    if (info /= 0) call error('Error in dgesv', 'exx_process')
    vref_aux(1:naux_overlap) = r1_tmp2(1)*auxcharg_overlap(:) &
                               +r1_tmp2(2)*homo_overlap(:)
    call memory_release(r2_tmp1)
    ! transform back
    call dgemv_x('N', naux, naux_overlap, 1d0, &
                 trans_mat_overlap, naux, &
                 vref_aux, 1, 0d0, r1_tmp1, 1)
    if (iverbose > 1) then
      call opm_test_charge(r1_tmp1, auxcharg, naux, 'Vref2')
      write (iout, *) ''
    end if
    call basis_vec_contract_A('OEP', 'ORBITAL', 'J', r1_tmp1, vref_ao, 0d0)
    vref_aux_up = r1_tmp1
    call memory_release(r1_tmp1)

  end subroutine vref_const_and_homo

  ! Calculates reference potential for the case of charge and homo condition
  ! with Fermi-Amaldi treatment of charge part
  subroutine vref_fa_and_homo(trans_mat_overlap, homo_overlap, &
                              auxcharg_overlap, int_fij, vxmonoonoo, auxcharg, &
                              vref_ao, vref_aux_up)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_parameters, only: iverbose
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat_overlap(:, :)
    real(dp), intent(in) :: auxcharg_overlap(:)
    real(dp), intent(in) :: homo_overlap(:)
    real(dp), intent(in) :: vxmonoonoo
    real(dp), intent(in) :: auxcharg(:)
    real(dp), intent(out) :: vref_ao(:)
    real(dp), intent(out) :: vref_aux_up(:)
    real(dp), intent(in) :: int_fij(:, :, :)
    real(dp), external :: ddot_x
    integer :: info, i, j, naux, naux_overlap
    integer, pointer :: ipiv(:)
    real(dp), pointer :: vref_aux(:), r1_tmp1(:), r1_tmp2(:), r2_tmp1(:, :)
    real(dp), pointer :: u(:), vc2(:)

    naux = size(trans_mat_overlap, 1)
    naux_overlap = size(trans_mat_overlap, 2)

    u => memory_allocate(naux)
    u = 0d0
    do i = 1, naux
      do j = 1, size(int_fij, 2)
        u(i) = u(i)+2d0*int_fij(i, j, j)
      end do
    end do

    vc2 => memory_allocate(naux_overlap)
    vc2 = 0d0
    call dgemv_x('T', naux, naux_overlap, 1d0, trans_mat_overlap, naux, &
                 u, 1, 0d0, vc2, 1)

    r1_tmp1 => memory_allocate(naux)
    r1_tmp1 = 0d0
    vref_aux => memory_allocate(naux_overlap)
    vref_aux = 0d0
    r2_tmp1 => memory_allocate(2, 2)
    r2_tmp1 = 0d0
    r1_tmp2 => memory_allocate(2)
    r1_tmp2 = 0d0
    r2_tmp1(1, 1) = ddot_x(naux_overlap, auxcharg_overlap, 1, &
                           auxcharg_overlap, 1)
    r2_tmp1(1, 2) = ddot_x(naux_overlap, auxcharg_overlap, 1, vc2, 1)
    r2_tmp1(2, 1) = ddot_x(naux_overlap, homo_overlap, 1, auxcharg_overlap, 1)
    r2_tmp1(2, 2) = ddot_x(naux_overlap, homo_overlap, 1, vc2, 1)
    r1_tmp2(1) = -1d0
    r1_tmp2(2) = vxmonoonoo
    info = 0
    ipiv => memory_allocate_integer(2)
    call dgesv_x(2, 1, r2_tmp1, 2, ipiv, r1_tmp2, 2, info)
    if (info /= 0) call error('Error in dgesv', 'exx_process')
    vref_aux(1:naux_overlap) = r1_tmp2(1)*auxcharg_overlap(:) &
                               +r1_tmp2(2)*vc2(:)
    call memory_release(r2_tmp1)

    ! transform back
    call dgemv_x('N', naux, naux_overlap, 1d0, &
                 trans_mat_overlap, naux, &
                 vref_aux, 1, 0d0, r1_tmp1, 1)
    if (iverbose > 1) then
      call opm_test_charge(r1_tmp1, auxcharg, naux, 'Vref2')
      write (iout, *) ''
    end if
    call basis_vec_contract_A('OEP', 'ORBITAL', 'J', r1_tmp1, vref_ao, 0d0)
    vref_aux_up = r1_tmp1
    call memory_release(r1_tmp1)

  end subroutine vref_fa_and_homo

  ! Calculates reference potential for the case of charge and homo condition
  ! with Fermi-Amaldi treatment of charge part for spin case
  subroutine vref_fa_and_homo_spin(trans_mat_overlap, homo_overlap, &
                                   auxcharg_overlap, int_fij, int_fij_b, &
                                   vxmonoonoo, auxcharg, vref_ao, vref_aux_up)
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_parameters, only: iverbose, l_vref_fa_sameab, l_spin_sym
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: trans_mat_overlap(:, :)
    real(dp), intent(in) :: auxcharg_overlap(:)
    real(dp), intent(in) :: homo_overlap(:)
    real(dp), intent(in) :: vxmonoonoo
    real(dp), intent(in) :: auxcharg(:)
    real(dp), intent(out) :: vref_ao(:)
    real(dp), intent(out) :: vref_aux_up(:)
    real(dp), intent(in) :: int_fij(:, :, :), int_fij_b(:, :, :)
    real(dp), external :: ddot_x
    integer :: info, i, j, naux, naux_overlap
    integer, pointer :: ipiv(:)
    real(dp), pointer :: vref_aux(:), r1_tmp1(:), r1_tmp2(:), r2_tmp1(:, :)
    real(dp), pointer :: u(:), vc2(:)

    naux = size(trans_mat_overlap, 1)
    naux_overlap = size(trans_mat_overlap, 2)

    u => memory_allocate(naux)
    u = 0d0
    if ((l_vref_fa_sameab) .or. (l_spin_sym)) then
      do i = 1, naux
        do j = 1, size(int_fij, 3)
          u(i) = u(i)+int_fij(i, j, j)
        end do
      end do
      do i = 1, naux
        do j = 1, size(int_fij_b, 3)
          u(i) = u(i)+int_fij_b(i, j, j)
        end do
      end do
    else
      do i = 1, naux
        do j = 1, size(int_fij, 3)
          u(i) = u(i)+int_fij(i, j, j)
        end do
      end do
    end if

    vc2 => memory_allocate(naux_overlap)
    vc2 = 0d0
    call dgemv_x('T', naux, naux_overlap, 1d0, trans_mat_overlap, naux, &
                 u, 1, 0d0, vc2, 1)

    r1_tmp1 => memory_allocate(naux)
    r1_tmp1 = 0d0
    vref_aux => memory_allocate(naux_overlap)
    vref_aux = 0d0
    r2_tmp1 => memory_allocate(2, 2)
    r2_tmp1 = 0d0
    r1_tmp2 => memory_allocate(2)
    r1_tmp2 = 0d0
    r2_tmp1(1, 1) = ddot_x(naux_overlap, auxcharg_overlap, 1, &
                           auxcharg_overlap, 1)
    r2_tmp1(1, 2) = ddot_x(naux_overlap, auxcharg_overlap, 1, vc2, 1)
    r2_tmp1(2, 1) = ddot_x(naux_overlap, homo_overlap, 1, auxcharg_overlap, 1)
    r2_tmp1(2, 2) = ddot_x(naux_overlap, homo_overlap, 1, vc2, 1)
    r1_tmp2(1) = -1d0
    r1_tmp2(2) = vxmonoonoo
    info = 0
    ipiv => memory_allocate_integer(2)
    call dgesv_x(2, 1, r2_tmp1, 2, ipiv, r1_tmp2, 2, info)
    if (info /= 0) call error('Error in dgesv', 'exx_process')
    vref_aux(1:naux_overlap) = r1_tmp2(1)*auxcharg_overlap(:) &
                               +r1_tmp2(2)*vc2(:)
    call memory_release(r2_tmp1)

    ! transform back
    call dgemv_x('N', naux, naux_overlap, 1d0, &
                 trans_mat_overlap, naux, &
                 vref_aux, 1, 0d0, r1_tmp1, 1)
    if (iverbose > 1) then
      call opm_test_charge(r1_tmp1, auxcharg, naux, 'Vref2')
      write (iout, *) ''
    end if
    call basis_vec_contract_A('OEP', 'ORBITAL', 'J', r1_tmp1, vref_ao, 0d0)
    vref_aux_up = r1_tmp1
    call memory_release(r1_tmp1)

  end subroutine vref_fa_and_homo_spin

end module acfd_basis_process
