module acfd_utils
  implicit none
  private
  public :: dfit_init, closed_e_ref, open_e_ref, identity_mat, vxao_to_vxmo, &
            matrix_spectral, exx_oep, exx_oep_spin_symmetrized, &
            print_eig, deallocate_arrays, x0_exx, rhs_exx, get_x0, &
            test_orthogonality, oep_solve, ssquared, subtract_vref, uscexx_pots
contains

subroutine uscexx_pots(dena, denb, vj_out, vxa_out, vxb_out)
  use common_cbas, only: ntqg, ntg, ntdg
  use memory, only: memory_allocate, memory_release
  use iso_fortran_env, only: dp => real64
  implicit none
  real(dp), intent(in) :: dena(ntqg), denb(ntqg)
  real(dp), optional, intent(out) :: vj_out(ntdg), vxa_out(ntdg), vxb_out(ntdg)
  real(dp), pointer :: r1_tmp1(:)
  real(dp) :: vj(ntdg), vxa(ntdg), vxb(ntdg)
  
  ! calculate nonlocal exchange and coulomb matrix
  vj = 0d0
  vxa = 0d0
  vxb = 0d0
  r1_tmp1 => memory_allocate(ntqg*2)
  call hfma(dena, dena, vj, vxa, -1, 1, 0, .true., .true., &
            r1_tmp1, .false., .false.)
  call hfma(dena+denb, dena-denb, vj, vxb, -1, 1, 0, .true., .true., &
            r1_tmp1, .false., .false.)
  vxb = vxa-vxb
  vj = (vj-(vxa+vxb))
  vxa = 2d0*vxa
  vxb = 2d0*vxb
  
  call memory_release(r1_tmp1)
  
  if (present(vj_out)) vj_out = vj
  if (present(vxa_out)) vxa_out = vxa
  if (present(vxb_out)) vxb_out = vxb
  
end subroutine uscexx_pots

  ! initializes density fitting
  subroutine dfit_init()
    implicit none
    integer :: icfit, nbas
    logical :: lfitcoul, lfitexch, lcoarse

    icfit = 3
    call set_inpi('CFIT_SCF', 8, icfit, 'CFIT')
    call dfhf_init(icfit, lfitcoul, lfitexch, nbas, 1, lcoarse, 0)

  end subroutine dfit_init

  ! calculates reference energy for spin-restricted calculation
  function closed_e_ref(den, orb, ldfit) result(energy)
    use common_cbas, only: ntdg, ntqg
    use common_code, only: ho1
    use common_molen, only: ekern
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: den(ntdg) ! density
    real(dp), intent(in) :: orb(ntqg) ! orbitals
    real(dp) :: energy
    logical, intent(in)  :: ldfit ! density fitting is enabled/disabled

    real(dp), pointer :: h0(:), vj(:), vx(:) ! potentials
    real(dp), external :: spur ! function

    h0 => memory_allocate(ntdg)
    h0 = 0d0
    vj => memory_allocate(ntdg)
    vj = 0d0
    vx => memory_allocate(ntdg)
    vx = 0d0

    ! hfma with ...1,0,0... provides vj which contains all two-electron part,
    ! vx is zero
    call hfma(den, den, vj, vx, 1, 0, 0, .true., .true., orb, ldfit, ldfit)
    call les(h0, ntdg, 1, ho1)

    energy = ekern+spur(den, h0)+0.5d0*(spur(den, vj)+spur(den, vx))

    call memory_release(h0)

  end function closed_e_ref

  ! calculates reference energy for spin-unrestricted calculations
  function open_e_ref(dena, denb, orba, orbb, ldfit) result(energy)
    use common_cbas, only: ntdg, ntqg
    use common_code, only: ho1
    use common_molen, only: ekern
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: dena(ntdg), denb(ntdg) ! density for alpha/beta spin
    real(dp), intent(in) :: orba(ntqg), orbb(ntqg) ! orbitals for alpha/beta spin
    real(dp) :: energy
    logical, intent(in) :: ldfit ! density fitting is enabled/disabled

    real(dp), pointer :: h0(:), vj(:), vxa(:), vxb(:) ! potentials
    real(dp), pointer :: denc(:), dens(:), orb(:) ! auxiliary variable
    real(dp), external :: spur ! function

    h0 => memory_allocate(ntdg)
    h0 = 0d0
    vj => memory_allocate(ntdg)
    vj = 0d0
    vxa => memory_allocate(ntdg)
    vxa = 0d0
    vxb => memory_allocate(ntdg)
    vxb = 0d0

    denc => memory_allocate(ntdg)
    denc = dena+denb
    dens => memory_allocate(ntdg)
    dens = dena-denb
    orb => memory_allocate(ntqg*2)
    orb = 0d0
    orb(1:ntqg) = orba
    call hfma(dena, dena, vj, vxa, -1, 1, 0, .true., .true., orb, ldfit, ldfit)
    orb(1:ntqg) = orba
    orb(1+ntqg:2*ntqg) = orbb
    call hfma(denc, dens, vj, vxb, -1, 1, 0, .true., .true., orb, ldfit, ldfit)
    vxb = vxa-vxb
    vj = (vj-(vxa+vxb))

    call les(h0, ntdg, 1, ho1)
    vxa = 2d0*vxa
    vxb = 2d0*vxb

    energy = ekern+spur(dena, h0)+spur(denb, h0) &
             +0.5d0*(spur(dena, vj)+spur(denb, vj) &
                     +spur(dena, vxa) &
                     +spur(denb, vxb))

    call memory_release(h0)

  end function open_e_ref

  ! Calculates the response matrix with different prefactors
  !
  ! \param [out] rX0 response matrix
  ! \param [in] eig orbital eigenvalues
  ! \param [in] omega frequency
  ! \param [in] noo number of occupied orbitals
  ! \param [in] nvo number of virtual orbitals
  ! \param [in] naux number of auxiliary basis functions
  ! \param [in] int_fai 3-index integrals (basis,virtual,occupied)
  ! \param [in] cmode determines prefactor and threshold
  subroutine get_x0(rX0, eig, omega, noo, ntg, naux, cmode)
    use memory, only: memory_allocate, memory_release
    use acfd_wf, only: int_all
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: noo, ntg, naux
    real(dp), intent(out) :: rX0(naux, naux)
    real(dp), intent(in) :: eig(ntg), omega
    character(*), intent(in) :: cmode

    integer :: i, a
    real(dp), pointer :: lambda(:, :)
    real(dp), pointer :: r3_tmp(:, :, :)

    lambda => memory_allocate(ntg, ntg)
    call get_lambda(lambda, eig, omega, noo, ntg, cmode)
    r3_tmp => memory_allocate(naux, ntg, ntg)
    r3_tmp = 0d0
    r3_tmp = int_all
    do i = 1, ntg
      do a = 1, ntg
        call dscal_x(naux, lambda(a, i), r3_tmp(:, a, i), 1)
      end do
    end do

    call dgemm_x('N', 'T', naux, naux, ntg*ntg, 1d0, &
                 int_all, naux, r3_tmp, naux, 0d0, rX0, naux)

    call memory_release(lambda)

  end subroutine get_x0

  ! This subroutine creates unit matrix a of given sizen n
  subroutine identity_mat(a, n)
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: n
    real(dp), intent(out):: a(n, n)
    integer :: i

    call fzero(a, n*n)
    do i = 1, n
      a(i, i) = 1d0
    end do

  end subroutine identity_mat

  ! diagonalizes a symmetric matrix
  subroutine matrix_spectral(matrix, evec, eval, n)
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: n ! size of the matrix
    real(dp), intent(in) :: matrix(n, n) ! matrix
    real(dp), intent(out) :: evec(n, n) ! eigenvectors
    real(dp), intent(out) :: eval(n) ! eigenvalues

    integer :: iwork, info
    real(dp) :: work(3*n)

    iwork = 3*n
    info = 0
    evec = matrix
    call dsyev_x('V', 'U', n, evec, n, eval, work, iwork, info)
    if (info .ne. 0) call error('error in dsyev', 'matrix_spectral')

  end subroutine matrix_spectral

  ! calculates lambda with different prefactors
  ! lambda = $\frac{\Delta\epsilon}{\Delta\epsilon^2+omega^2}$
  subroutine get_lambda(lambda, eig, omega, noo, ntg, cmode)
    use acfd_parameters, only: thr_rpa, thr_oep
    use acfd_wf, only: occa, occb
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: noo, ntg ! number of occupied and virtual orbitals
    real(dp), intent(out) :: lambda(ntg, ntg) ! lambda
    real(dp), intent(in) :: eig(ntg) ! orbital eigenvalues
    real(dp), intent(in) ::  omega ! frequency
    character(*), intent(in) :: cmode ! determines prefactor and threshold

    integer :: i, a
    real(dp) :: eia, thr, fac, occ(ntg)

    thr = 1d99 ! to suppress warning: 'thr' may be used uninitialized
    fac = 0d0 ! to suppress warning:  'fac' may be used uninitialized

    ! choose prefactor and threshold according to method
    if (trim(cmode) == 'EXX') then
      fac = 1d0
      thr = thr_oep
    elseif (trim(cmode) == 'RPA') then
      fac = 4d0
      thr = thr_oep
    elseif (trim(cmode) == 'RIRPA') then
      fac = 4d0
      thr = thr_rpa
    else
      call error('unknown cmode', 'get_lambda')
    end if

    occ = 0d0
    if (noo == nint(sum(occa))) then
      occ = occa
    else if (noo == nint(sum(occb))) then
      occ = occb
    else
      write(*,*) "something is wrong"
    end if

    lambda = 0d0
    do i = 1, ntg
      do a = 1, ntg
        eia = eig(a)-eig(i)
        if (abs(eia) .lt. thr) eia = sign(thr, eia)
        if (i /= a) lambda(a, i) = occ(i)*fac*eia/(eia**2+omega**2)
      end do
    end do
  end subroutine get_lambda

  subroutine get_lambda2(lambda, eig, omega, noo, ntg, cmode, occ)
    use acfd_parameters, only: thr_rpa, thr_oep
!    use acfd_wf, only: occa, occb
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: noo, ntg ! number of occupied and virtual orbitals
    real(dp), intent(out) :: lambda(ntg, ntg) ! lambda
    real(dp), intent(in) :: eig(ntg), occ(ntg) ! orbital eigenvalues
    real(dp), intent(in) ::  omega ! frequency
    character(*), intent(in) :: cmode ! determines prefactor and threshold

    integer :: i, a
    real(dp) :: eia, thr, fac!, occ(ntg)

    thr = 1d99 ! to suppress warning: 'thr' may be used uninitialized
    fac = 0d0 ! to suppress warning:  'fac' may be used uninitialized

    ! choose prefactor and threshold according to method
    if (trim(cmode) == 'EXX') then
      fac = 1d0 
      thr = thr_oep
    elseif (trim(cmode) == 'RPA') then
      fac = 4d0 
      thr = thr_oep
    elseif (trim(cmode) == 'RIRPA') then
      fac = 4d0 
      thr = thr_rpa
    else
      call error('unknown cmode', 'get_lambda')
    end if

!    occ = 0d0 
!    if (noo == nint(sum(occa))) then
!      occ = occa
!    else if (noo == nint(sum(occb))) then
!      occ = occb
!    else
!      write(*,*) "something is wrong"
!    end if

    lambda = 0d0 
    do i = 1, ntg 
      do a = 1, ntg 
        eia = eig(a)-eig(i)
        if (abs(eia) .lt. thr) eia = sign(thr, eia)
        if (i /= a) lambda(a, i) = occ(i)*fac*eia/(eia**2+omega**2)
      end do
    end do
  end subroutine get_lambda2

  ! convert local exchange potential for AO to MO representation
  subroutine vxao_to_vxmo(vxao, orb, ntdg, ntqg, vxmo)
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), pointer :: r1_tmp1(:)
    integer, intent(in) :: ntdg, ntqg
    real(dp), intent(in) :: vxao(ntdg), orb(ntqg)
    real(dp), intent(out) :: vxmo(ntqg)

    vxmo = 0d0
    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp1 = 0d0
    call expan(vxao, vxmo, 1, 1, 1)
    call mmult(orb, vxmo, r1_tmp1, 1, 1, 0, -1, 0, 0)
    call mmult(r1_tmp1, orb, vxmo, 1, 1, 0, 0, 0, 0)
    call memory_release(r1_tmp1)
  end subroutine

  ! calculates EXX exchange potential
  subroutine exx_oep(vx, vsol, noo, naux, naux_constraint, orb, eig, vxmo, &
                     trans_mat, vref_ao, auxcharg, homo, int_fai, spin, &
                     vref_aux_up, den, itc)
    use common_cbas, only: ntqg, ntdg, ntg
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_derived_datatypes, only: matrix_type, vector_type
    use acfd_parameters, only: thr_solve, l_plot_always, l_homo, iverbose
    use acfd_plot, only: plot_potential
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: noo ! number of occupied orbitals
    integer, intent(in) :: naux ! number of auxiliary basis functions
    integer, intent(in) :: naux_constraint ! number of auxiliary basis functions of processed basis
    real(dp), intent(in) :: orb(ntqg) ! orbitals
    real(dp), intent(in) :: eig(ntg) ! eigenvalues
    real(dp), intent(in) :: den(ntdg) ! density
    real(dp), intent(inout) :: vx(ntdg) ! non-local exchange potential
    real(dp), intent(inout) :: vxmo(ntg, ntg) ! non-local exchange potential in MO basis
    type(matrix_type), intent(in) :: trans_mat ! transformation matrix for processed auxiliary basis
    type(vector_type), intent(in) :: vref_ao ! reference potential in AO basis
    type(vector_type), intent(in) :: vref_aux_up
    real(dp), intent(in) :: int_fai(naux, ntg-noo, noo) ! integrals (basis,virtual,occupied)
    real(dp), intent(out) :: vsol(naux) ! solution vector
    real(dp), intent(in) :: auxcharg(naux) ! charge constraint vector
    real(dp), intent(in) :: homo(naux) ! HOMO constraint vector
    character(*) :: spin, itc ! spin and number of current iteration

    integer, pointer :: ibase(:) ! pointers for stack management

    integer :: nvo ! number of virtual orbitals
    real(dp), pointer :: rX0(:, :) ! response matrix
    real(dp), pointer :: rhs(:) ! right hand side

    ! auxiliary variables
    real(dp), pointer :: r1_tmp1(:)

    ! create a pointer to the bottom of the stack used in this driver
    ibase => memory_allocate_integer(1)

    nvo = ntg-noo

    ! subtract vref from local exchange
    call subtract_vref(vref_ao, orb, vxmo)

    ! construct and transform to processed basis response matrix
    rX0 => memory_allocate(naux_constraint, naux_constraint)
!    call x0_exx(eig, int_fai, trans_mat, noo, ntg, naux, naux_constraint, &
!                rX0)

    ! construct and transform to processed basis right hand side
    rhs => memory_allocate(naux_constraint)
!    call rhs_exx(eig, vxmo, trans_mat, int_fai, nvo, ntg, naux, &
!                 naux_constraint, rhs)

    ! solve oep equation
    call oep_solve(rX0, rhs, naux_constraint, thr_solve)

    ! backtransformation of solution vector
    call dgemv_x('N', naux, naux_constraint, 1d0, &
                 trans_mat%mat, naux, rhs, 1, 0d0, vsol, 1)

    ! transform local potential to AO basis
    r1_tmp1 => memory_allocate(ntdg)
    r1_tmp1 = 0d0
    call basis_vec_contract_A('OEP', 'ORBITAL', 'J', vsol, r1_tmp1, 0d0)

    ! test orthogonality
    call test_orthogonality(auxcharg, vsol, naux, &
                                 'charge orthogonality of EXX')
    if (l_homo) then
      call test_orthogonality(homo, vsol, naux, &
                                   'HOMO orthogonality of EXX')
    end if

    ! add reference potential
    call daxpy_x(ntdg, 1d0, vref_ao%mat, 1, r1_tmp1, 1)

    if ((l_homo) .and. (iverbose > 1)) then
      write (iout, *) ''
      call test_homo(r1_tmp1, orb, noo, 'Vx(L)')
      call test_homo(vx, orb, noo, 'Vx(NL)')
      call test_homo(vref_ao%mat, orb, noo, 'Vref')
      write (iout, *) ''
    end if

    ! replace non-local exchange with local exchange
    vx = r1_tmp1

    if (l_plot_always) then
      ! plot exchange potential
      call plot_potential(orb, den, vsol, naux, itc, 'OEP', noo, &
                               'vx'//trim(spin), vref_aux_up)
    end if

    call memory_release(ibase)

  end subroutine exx_oep

  ! Calculates EXX exchange potential with one OEP equation for alpha and beta spin combined
  subroutine exx_oep_spin_symmetrized(vxa, vxb, vsol, noa, nob, naux, &
                                      naux_constraint, orba, orbb, eiga, eigb, &
                                      vxmoa, vxmob, trans_mat, vref_ao, &
                                      auxcharg, homo, &
                                      vref_aux_up, dena, denb, itc, vxmob2)
    use common_cbas, only: ntqg, ntdg, ntg
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use acfd_derived_datatypes, only: matrix_type, vector_type
    use acfd_parameters, only: thr_solve, l_plot_always, l_homo, iverbose
    use acfd_plot, only: plot_potential
    use acfd_wf, only: occa, occb1, occb2, occb
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: noa, nob ! number of occupied orbitals
    integer, intent(in) :: naux ! number of auxiliary basis functions
    integer, intent(in) :: naux_constraint ! number of auxiliary basis functions of processed basis
    real(dp), intent(in) :: orba(ntqg), eiga(ntg), dena(ntdg) ! orbitals, eigenvalues, densities
    real(dp), intent(in) :: orbb(ntqg), eigb(ntg), denb(ntdg)
    real(dp), intent(inout) :: vxa(ntdg), vxb(ntdg) ! non-local exchange potential
    real(dp), intent(inout) :: vxmoa(ntg, ntg), vxmob(ntg, ntg), vxmob2(ntg, ntg) ! non-local exchange potential in MO basis
    type(matrix_type), intent(in) :: trans_mat ! transformation matrix for processed auxiliary basis
    type(vector_type), intent(inout) :: vref_ao ! reference potential in AO basis
    type(vector_type), intent(in) :: vref_aux_up
    real(dp), intent(out) :: vsol(naux) ! solution vector
    real(dp), intent(in) :: auxcharg(naux) ! charge constraint vector
    real(dp), intent(in) :: homo(naux) ! HOMO constraint vector
    character(*) :: itc

    integer, pointer :: ibase(:) ! pointers for stack management

    ! auxiliary variables
    real(dp), pointer :: r1_tmp1(:)

    !integer :: nva, nvb ! number of virtual orbitals
    real(dp), pointer :: rX0a(:, :), rX0b(:, :) ! response matrix
    real(dp), pointer :: rhsa(:), rhsb(:), rhsb2(:) ! right hand side

    real(dp) :: aux(ntdg)

    ! create a pointer to the bottom of the stack used in this driver
    ibase => memory_allocate_integer(1)

    !nva = ntg-noa
    !nvb = ntg-nob

    ! subtract vref from local exchange
    call subtract_vref(vref_ao, orba, vxmoa) ! alpha
    call subtract_vref(vref_ao, orba, vxmob2)

!    call subtract_vref(vref_ao, orbb, vxmob) ! beta
    call subtract_vref(vref_ao, orbb, vxmob)

    ! construct and transform to processed basis X0 and rhs
    ! alpha
    rX0a => memory_allocate(naux_constraint, naux_constraint)
    call x0_exx(eiga, trans_mat, noa, ntg, naux, naux_constraint, &
                rX0a)
    rhsa => memory_allocate(naux_constraint)
    call rhs_exx(eiga, vxmoa, trans_mat, ntg, noa, naux, & !??? noa is not important
                 naux_constraint, rhsa, occb1)

    rhsb2 => memory_allocate(naux_constraint)
    call rhs_exx(eiga, vxmob2, trans_mat, ntg, noa, naux, &
                 naux_constraint, rhsb2, occb2) !??? noa is not important

    ! beta
    rX0b => memory_allocate(naux_constraint, naux_constraint)
    call x0_exx(eigb, trans_mat, nob, ntg, naux, naux_constraint, &
                rX0b)
    rhsb => memory_allocate(naux_constraint)
    call rhs_exx(eigb, vxmob, trans_mat, ntg, nob, naux, &
                 naux_constraint, rhsb, occb)

    ! combine alpha and beta rhs and lhs
    rX0a = rX0a+rX0b
    rhsa = 0.5d0*rhsa+rhsb+0.5d0*rhsb2
    call memory_release(rX0b)

    ! solve oep equation
    call oep_solve(rX0a, rhsa, naux_constraint, 2d0*thr_solve)

    ! backtransformation of solution vector
    call dgemv_x('N', naux, naux_constraint, 1d0, trans_mat%mat, naux, &
                 rhsa, 1, 0d0, vsol, 1)

    ! transform local potential to AO basis
    r1_tmp1 => memory_allocate(ntdg)
    call basis_vec_contract_A('OEP', 'ORBITAL', 'J', vsol, r1_tmp1, 0d0)

    ! test orthogonality
    call test_orthogonality(auxcharg, vsol, naux, 'charge orthogonality')
    if (l_homo) then
      call test_orthogonality(homo, vsol, naux, &
                                   'HOMO orthogonality of EXX')
    end if

    ! add reference potential
    call daxpy_x(ntdg, 1d0, vref_ao%mat, 1, r1_tmp1, 1)

    if ((l_homo) .and. (iverbose > 1)) then
      write (iout, *) ''
      call test_homo(r1_tmp1, orba, noa, 'Vx(L)')
      call test_homo(vxa, orba, noa, 'Vx(NL)')
      call test_homo(vref_ao%mat, orba, noa, 'Vref')
      write (iout, *) ''
    end if

    ! replace non-local exchange with local exchange
    vxa = r1_tmp1
    vxb = r1_tmp1

    if (l_plot_always) then
      ! plot exchange potential
      call plot_potential(orba, dena, vsol, naux, itc, 'OEP', noa, 'vxa', &
                          vref_aux_up)
      call plot_potential(orbb, denb, vsol, naux, itc, 'OEP', nob, 'vxb', &
                          vref_aux_up)
    end if

    call memory_release(ibase)

  end subroutine exx_oep_spin_symmetrized

  ! subtract vref from local exchange
  subroutine subtract_vref(vref_ao, orb, vxmo)
    use common_cbas, only: ntqg, ntg
    use memory, only: memory_allocate, memory_release
    use acfd_derived_datatypes, only: vector_type
    use iso_fortran_env, only: dp => real64
    implicit none
    type(vector_type), intent(in) :: vref_ao ! reference potential in AO basis
    real(dp), intent(in) :: orb(ntqg) ! orbitals
    real(dp), intent(inout) :: vxmo(ntg, ntg) ! non-local exchange potential in MO basis

    real(dp), pointer :: r1_tmp1(:), r1_tmp2(:) ! auxiliary variables

    r1_tmp1 => memory_allocate(ntqg)
    r1_tmp2 => memory_allocate(ntqg)
    call expan(vref_ao%mat, r1_tmp1, 1, 1, 1)
    call mmult(orb, r1_tmp1, r1_tmp2, 1, 1, 0, -1, 0, 0)
    call mmult(r1_tmp2, orb, r1_tmp1, 1, 1, 0, 0, 0, 0)
    call daxpy_x(ntqg, -1d0, r1_tmp1, 1, vxmo, 1)
    call memory_release(r1_tmp1)
  end subroutine subtract_vref

  ! construct and transform to processed basis rhs for EXX OEP calculations
  subroutine rhs_exx(eig, vxmo, trans_mat, ntg, noo, naux, &
                     naux_constraint, rhs, occ)
!    use common_cbas, only: ntg
    use memory, only: memory_allocate, memory_release
    use acfd_wf, only: int_all
    use acfd_derived_datatypes, only: matrix_type
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: ntg, noo, naux, naux_constraint
    real(dp), intent(in) :: eig(ntg), occ(ntg) ! eigenvalues
    real(dp), intent(in) :: vxmo(ntg, ntg) ! non-local exchange potential in MO basis
    type(matrix_type), intent(in) :: trans_mat
    real(dp), intent(out) :: rhs(naux_constraint)

    integer :: i, a, iaux
    real(dp), pointer :: lambda(:, :) ! 1/(eig(virt)-eig(occ))
    real(dp), pointer :: r1_tmp(:)

    rhs = 0d0

    ! right hand side in unprocessed auxbasis
    r1_tmp => memory_allocate(naux)
    r1_tmp = 0d0
    lambda => memory_allocate(ntg, ntg)
    lambda = 0d0
    call get_lambda2(lambda, eig, 0d0, noo, ntg, 'EXX', occ)
    do i = 1, ntg
      do a = 1, ntg
        do iaux = 1, naux
          r1_tmp(iaux) = r1_tmp(iaux) &
                         +vxmo(a, i)*lambda(a, i)*int_all(iaux, a, i)
        end do
      end do
    end do
    call memory_release(lambda)

    ! transform right hand side to processed basis
    call dgemv_x('T', naux, naux_constraint, 1d0, &
                 trans_mat%mat, naux, r1_tmp, 1, 0d0, rhs, 1)

    call memory_release(r1_tmp)
  end subroutine rhs_exx

  ! construct and transform to processed basis response function in
  ! EXX OEP calculations
  subroutine x0_exx(eig, trans_mat, noo, ntg, naux, naux_constraint, &
                    rX0)
    !use common_cbas, only: ntg
    use memory, only: memory_allocate, memory_release
    use acfd_derived_datatypes, only: matrix_type
    use iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: ntg, noo, naux, naux_constraint
    real(dp), intent(in) :: eig(ntg) ! eigenvalues
    type(matrix_type), intent(in) :: trans_mat
    real(dp), intent(out) :: rX0(naux_constraint, naux_constraint)

    real(dp), pointer :: r1_tmp(:)
    real(dp), pointer :: r2_tmp(:, :)

    rX0 = 0d0

    ! response matrix in unprocessed auxbasis
    r2_tmp => memory_allocate(naux, naux)
    r2_tmp = 0d0
    call get_x0(r2_tmp, eig, 0d0, noo, ntg, naux, 'EXX')

    ! transform response matrix and right hand side to processed basis
    r1_tmp => memory_allocate(naux*naux_constraint)
    r1_tmp = 0d0
    call dgemm_x('T', 'N', naux_constraint, naux, naux, 1d0, &
                 trans_mat%mat, naux, &
                 r2_tmp, naux, &
                 0d0, r1_tmp, naux_constraint)
    call dgemm_x('N', 'N', naux_constraint, naux_constraint, naux, 1d0, &
                 r1_tmp, naux_constraint, &
                 trans_mat%mat, naux, &
                 0d0, rX0, naux_constraint)
    call memory_release(r2_tmp)
  end subroutine x0_exx

  ! solves OEP equation
  subroutine oep_solve(rX0, rhs, naux, thr)
    use common_tapes, only: iout
    use acfd_parameters, only: csolvemeth, iverbose
    use memory, only: memory_allocate, memory_release, memory_allocate_integer
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: naux ! number of auxiliary basis functions
    real(dp), intent(inout) :: rX0(naux, naux) ! response matrix
    real(dp), intent(inout) :: rhs(naux) ! right hand side (input) and solution vector
    real(dp), intent(in) :: thr

    ! auxiliary variables
    integer :: info, iaux
    integer, pointer :: ibase(:), ipiv(:)
    real(dp), pointer :: work(:), r1_tmp1(:)
    real(dp), pointer :: r2_tmp1(:, :), r2_tmp2(:, :)

    ibase => memory_allocate_integer(1)

    ! use GESV
    if (trim(csolvemeth) .eq. 'GESV') then
      r2_tmp1 => memory_allocate(naux, naux)
      r2_tmp1 = rX0
      ipiv => memory_allocate_integer(naux)
      ipiv = 0
      info = 0
      call dgesv_x(naux, 1, r2_tmp1, naux, ipiv, rhs, naux, info)
      if (info .ne. 0) call error('error in dgesv', 'oep_solve')

      ! use TSVD (previously SVD)
    elseif (trim(csolvemeth) .eq. 'TSVD') then
      r2_tmp1 => memory_allocate(naux, naux)
      r2_tmp1 = rX0
      r1_tmp1 => memory_allocate(naux)
      r1_tmp1 = 0d0
      work => memory_allocate(3*naux)
      work = 0d0
      info = 0
      call dsyev_x('V', 'U', naux, r2_tmp1, naux, r1_tmp1, work, 3*naux, info)
      if (info .ne. 0) call error('error in dsyev', 'oep_solve')
      r2_tmp2 => memory_allocate(naux, naux)
      r2_tmp2 = r2_tmp1
      do iaux = 1, naux
        if (abs(r1_tmp1(iaux)) .gt. thr) then
          call dscal_x(naux, 1d0/r1_tmp1(iaux), r2_tmp2(:, iaux), 1)
        else
          if (iverbose > 0) write (iout, '(1x,a,i5,e16.8)') &
            'SVD filtered out:', iaux, r1_tmp1(iaux)
          r2_tmp2(:, iaux) = 0d0
        end if
      end do
      call dgemm_x('N', 'T', naux, naux, naux, 1d0, r2_tmp2, naux, &
                   r2_tmp1, naux, 0d0, rX0, naux)
      r1_tmp1 = rhs
      call dgemv_x('N', naux, naux, 1d0, rX0, naux, r1_tmp1, 1, 0d0, rhs, 1)

      ! use GTSVD (previously TIKH2)
    elseif (trim(csolvemeth) .eq. 'GTSVD') then
      r2_tmp1 => memory_allocate(naux, naux)
      r2_tmp1 = rX0
      r1_tmp1 => memory_allocate(naux)
      r1_tmp1 = 0d0
      work => memory_allocate(3*naux)
      work = 0d0
      info = 0
      call dsyev_x('V', 'U', naux, r2_tmp1, naux, r1_tmp1, work, 3*naux, info)
      if (info .ne. 0) call error('error in dsyev', 'oep_solve')
      r2_tmp2 => memory_allocate(naux, naux)
      r2_tmp2 = r2_tmp1
      r1_tmp1 = r1_tmp1+thr
      do iaux = 1, naux
        if (abs(r1_tmp1(iaux)) .gt. thr) then
          call dscal_x(naux, 1d0/r1_tmp1(iaux), r2_tmp2(:, iaux), 1)
        else
          if (iverbose > 0) write (iout, '(1x,a,i5,e16.8)') &
            'SVD filtered out:', iaux, r1_tmp1(iaux)
          r2_tmp2(:, iaux) = 0d0
        end if
      end do
      call dgemm_x('N', 'T', naux, naux, naux, 1d0, r2_tmp2, naux, &
                   r2_tmp1, naux, 0d0, rX0, naux)
      r1_tmp1 = rhs
      call dgemv_x('N', naux, naux, 1d0, rX0, naux, r1_tmp1, 1, 0d0, rhs, 1)
    else
      call error('unknown csolvemeth', 'oep_solve')
    end if

    call memory_release(ibase)
  end subroutine oep_solve

  ! checks if two vectors are orthogonal
  subroutine test_orthogonality(vec1, vec2, n, cname)
    use common_tapes, only: iout
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: thr = 1d-8

    character(*), intent(in) :: cname ! string for output
    integer, intent(in) :: n ! length of vectors
    real(dp), intent(in) :: vec1(n), vec2(n) ! vector
    real(dp) :: overlap
    real(dp), external :: ddot_x

    overlap = 0d0
    overlap = ddot_x(n, vec1, 1, vec2, 1)
    if (abs(overlap) .gt. thr) &
      write (iout, '(1x,a,e16.8,2x,a,a,a)') &
      'Warning! Vectors not orthogonal: ', overlap, '(', cname, ')'
  end subroutine test_orthogonality

  ! print array in formatted form and with specified name
  subroutine print_eig(array, asize, arname)
    use common_tapes, only: iout
    use iso_fortran_env, only: dp => real64
    integer :: asize
    real(dp) :: array(asize)
    character(*) :: arname
    integer :: i, i1, i2

10  format(i4, '-', i4, 1x, 10(1x, e12.5))

    write (iout, *)
    write (iout, *) arname
    i2 = 0
    do
      i1 = i2+1
      i2 = i2+10
      if (i2 > asize) i2 = asize
      write (iout, 10) i1, i2, (array(i), i=i1, i2)
      if (i2 == asize) exit
    end do
    write (iout, *)

  end subroutine print_eig

  ! calculates and prints S**2 operator
  subroutine ssquared(dc, ds, s)
    use common_cref, only: ms2, nelec
    use common_cbas, only: ntdg, ntqg
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: dc(*), ds(*), s(*)
    integer :: i
    real(dp) :: ssq
    real(dp), pointer :: da(:), db(:), q3(:)
    real(dp), external :: tracem
    da => memory_allocate(ntqg)
    db => memory_allocate(ntqg)
    q3 => memory_allocate(2*ntqg)
    do i = 1, ntdg
      da(i) = 0.5d0*(dc(i)+ds(i))
      db(i) = 0.5d0*(dc(i)-ds(i))
    end do
    call tranop_mpp(db, s, q3, 1, 1, 1, 0)
    ssq = dble(ms2)*dble(ms2)*0.25d0+dble(nelec)*0.5d0-tracem(da, db, 1, 1)
    write (iout, "(/1x,'EXPECTATION VALUE OF S**2:',t36,f15.8)") ssq
    call memory_release(da)
  end subroutine ssquared

  subroutine test_homo(pot, orb, noo, cname)
    use common_cbas, only: ntqg, ntdg, ntg
    use common_tapes, only: iout
    use memory, only: memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none

    character(*), intent(in) :: cname ! string for output
    integer, intent(in) :: noo ! number of occupied orbitals
    real(dp), intent(in) :: orb(ntqg), pot(ntdg) ! potential and orbitals

    real(dp), pointer :: r2_tmp1(:, :), r2_tmp2(:, :)
    real(dp) :: v_homo

    r2_tmp1 => memory_allocate(ntg, ntg)
    r2_tmp1 = 0d0
    r2_tmp2 => memory_allocate(ntg, ntg)
    r2_tmp2 = 0d0
    call expan(pot, r2_tmp1, 1, 1, 1)
    call mmult(orb, r2_tmp1, r2_tmp2, 1, 1, 0, -1, 0, 0)
    call mmult(r2_tmp2, orb, r2_tmp1, 1, 1, 0, 0, 0, 0)
    v_homo = r2_tmp1(noo, noo)
    write (iout, '(1x,a,f16.8,2x,a,a,a)') 'v(HOMO) = ', v_homo, '(', cname, ')'
    call memory_release(r2_tmp1)
  end subroutine test_homo

  ! deallocates arrays on heap
  subroutine deallocate_arrays
    use acfd_auxiliary_basis

    if (allocated(trans_mat_overlap%mat)) deallocate (trans_mat_overlap%mat)

    ! arrays for RIRPA
    if (allocated(trans_mat_constraint%mat)) &
      deallocate (trans_mat_constraint%mat)
    if (allocated(trans_mat_charge%mat)) deallocate (trans_mat_charge%mat)

    ! array for SCEXX
    if (allocated(trans_mat_constraint_a%mat)) &
      deallocate (trans_mat_constraint_a%mat)
    if (allocated(trans_mat_constraint_b%mat)) &
      deallocate (trans_mat_constraint_b%mat)

  end subroutine deallocate_arrays
end module acfd_utils
