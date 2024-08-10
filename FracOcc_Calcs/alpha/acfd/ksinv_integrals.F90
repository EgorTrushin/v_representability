module ksinv_integrals
  implicit none
  private
  public :: build_fkl
contains

  subroutine ao_3idx( &
    nfun1, ngrp1, infg1, exp1, cgr1, sph1, &
    nfun2, ngrp2, infg2, exp2, cgr2, sph2, &
    nfun3, ngrp3, infg3, exp3, cgr3, sph3, &
    nclass, ints, int_fxx)
    use common_cbascode
    use common_big
    use memory
    use iso_fortran_env, only: dp => real64
    implicit none
    double precision, external :: second_x

    ! basis 1
    integer, intent(in) :: nfun1, ngrp1
    integer, dimension(ig_ndim, *), intent(in) :: infg1
    real(dp), dimension(*), intent(in) :: exp1, cgr1
    logical, intent(in) :: sph1
    ! basis 2
    integer, intent(in) :: nfun2, ngrp2
    integer, dimension(ig_ndim, *), intent(in) :: infg2
    real(dp), dimension(*), intent(in) :: exp2, cgr2
    logical, intent(in) :: sph2
    ! basis 3
    integer, intent(in) :: nfun3, ngrp3
    integer, dimension(ig_ndim, *), intent(in) :: infg3
    real(dp), dimension(*), intent(in) :: exp3, cgr3
    logical, intent(in) :: sph3

    integer, intent(in) :: nclass
    character(len=*), intent(in) :: ints(nclass)

    real(dp), intent(out) :: int_fxx(nfun3, nfun2, nfun1)

    real(dp) :: tim
    integer, pointer, dimension(:) :: ibase
    integer :: ifil_half, irec_half
    real(dp), pointer, dimension(:) :: zmx
    integer, pointer, dimension(:) :: i1lst, i2lst, i3lst, i3lst_out, i3cls_out
    integer :: ic, ib, ioff
    real(dp), pointer, dimension(:) :: abc, ibc
    real(dp) :: offr, thr
    integer :: ipr
    integer :: ngrp_batch_2, ngrp_batch_3
    integer :: nfun_batch_2, nfun_batch_3
    integer :: ngrpmax1, ngrpmax2, ngrpmax3
    integer :: igrp_2, igrp_3
    integer :: igrp_2_first, igrp_3_first
    integer :: nbatch2, nbatch3, nbatch3_out
    integer :: ngrp_actual_2, ngrp_actual_3
    integer :: mem_reserve, memory_available
    real(dp) :: dum
    integer :: idum
    real(dp) :: zcalc, znonz, ztot
    integer :: i, i1, i2, i3, ii3
    integer :: n23, n3read, nwritten
    character(8) :: string

    int_fxx(:, :, :) = 0.d0

    tim = second_x()
    ibase => memory_allocate_integer(1)

    if (nclass .gt. 1) call error('nclass.gt.1', 'ao_3idx')

    ! locations for half transformed integrals
    ifil_half = 5
    irec_half = 4001

    ! temporary threshold
    call get_inpf('THRINT', 'EXPLICIT', thr)

    ! initialize integral batcher
    call basis_3idx_batchs_init2( &
      infg1, sph1, ngrp1, &
      infg2, sph2, ngrp2, &
      infg3, sph3, ngrp3, &
      thr, nclass, ints, mem_reserve)
    mem_reserve = max(2000000, mem_reserve)   !HJW returned mem_reserve is incorrect.

    ! get extra space needed for batcher arguments
    zmx => memory_allocate(nfun3)
    i1lst => memory_allocate_integer(ngrp1)
    i2lst => memory_allocate_integer(ngrp2)
    i3lst => memory_allocate_integer(ngrp3)
    i3lst_out => memory_allocate_integer(nfun3)
    i3cls_out => memory_allocate_integer(nfun3)

    !
    !cccccccccccccccccccccccccc FIRST TRANSFORMATION ccccccccccccccccccccccccccc
    !

    ! some memory management
    memory_available = memory_remaining()-mem_reserve-10000

    ! get maximum group sizes
    ngrpmax1 = 0
    do i = 1, ngrp1
      ngrpmax1 = max(ngrpmax1, infg1(ig_func, i))
      i1lst(i) = i
    end do
    ngrpmax2 = 0
    do i = 1, ngrp2
      ngrpmax2 = max(ngrpmax2, infg2(ig_func, i))
    end do
    ngrpmax3 = 0
    do i = 1, ngrp3
      ngrpmax3 = max(ngrpmax3, infg3(ig_func, i))
    end do

    ipr = 0
    IF (ipr .GT. 1) THEN
      WRITE (6, 1) TRIM(ints(1))
      WRITE (6, 3) 'Total memory left', dble(memory_available)*1D-6
      WRITE (6, 3) 'Min memory for in-core transformation', &
        dble((nfun1+nfun1)*nfun2*nfun3)*1D-6
      WRITE (6, 3) 'Min memory for paged transformation', &
        dble((nfun1+nfun1)*ngrpmax2*ngrpmax3)*1D-6
1     FORMAT(/, 1x, 'Building 3-index integrals of class ', a, /, 50('-'))
2     FORMAT(1x, a, ':', t50, i8)
3     FORMAT(1x, a, ':', t50, f5.2, ' MW')
4     FORMAT(1x, a, ':', t50, f5.2, ' second')
    END IF

    if ((nfun1+nfun1)*nfun2*nfun3 .lt. memory_available) then
      ngrp_batch_2 = ngrp2
      ngrp_batch_3 = ngrp3
      nfun_batch_2 = nfun2
      nfun_batch_3 = nfun3
    else if ((nfun1+nfun1)*nfun2*ngrpmax3 .lt. memory_available) then
      ngrp_batch_2 = ngrp2
      nfun_batch_2 = nfun2
      nfun_batch_3 = memory_available/((nfun1+nfun1)*nfun2)
      ngrp_batch_3 = nfun_batch_3/ngrpmax3
    else if ((nfun1+nfun1)*ngrpmax2*ngrpmax3 .lt. memory_available) then
      nfun_batch_3 = ngrpmax3
      ngrp_batch_3 = 1
      nfun_batch_2 = memory_available/((nfun1+nfun1)*ngrpmax3)
      ngrp_batch_2 = nfun_batch_2/ngrpmax2
    else
      CALL error('Not enough memory for 3idx first tran', 'ao_3idx')
    end if

    if (ipr .gt. 0) then
      WRITE (6, 2) 'Number of batchs in 1st tran', &
        ngrp2*ngrp3/(ngrp_batch_2*ngrp_batch_3)
    end if

    ! build half transformed integrals
    abc => memory_allocate(nfun1*nfun_batch_2*nfun_batch_3)
    call fzero(abc, nfun1*nfun_batch_2*nfun_batch_3)
    ibc => memory_allocate(nfun1*nfun_batch_2*nfun_batch_3)

    zcalc = 0d0
    znonz = 0d0
    ztot = 0d0
    do igrp_3_first = 1, ngrp3, ngrp_batch_3
      ic = 1
      nbatch3 = 0
      ngrp_actual_3 = 0
      do igrp_3 = igrp_3_first, min(igrp_3_first+ngrp_batch_3-1, ngrp3)
        nbatch3 = nbatch3+infg3(ig_func, igrp_3)
        i3lst(ic) = igrp_3
        ic = ic+1
        ngrp_actual_3 = ngrp_actual_3+1
      end do
      do igrp_2_first = 1, ngrp2, ngrp_batch_2
        ib = 1
        nbatch2 = 0
        ngrp_actual_2 = 0
        do igrp_2 = igrp_2_first, min(igrp_2_first+ngrp_batch_2-1, ngrp2)
          nbatch2 = nbatch2+infg2(ig_func, igrp_2)
          i2lst(ib) = igrp_2
          ib = ib+1
          ngrp_actual_2 = ngrp_actual_2+1
        end do

        call basis_3idx_batchs2( &
          infg1, exp1, cgr1, sph1, &
          infg2, exp2, cgr2, sph2, &
          infg3, exp3, cgr3, sph3, &
          nclass, ints, &
          abc, nfun1, nbatch2, zmx, &
          i1lst, ngrp1, &
          i2lst, ngrp_actual_2, &
          i3lst, ngrp_actual_3, &
          i3lst_out, nbatch3_out, i3cls_out, &
          .FALSE., .FALSE., dum, idum, &
          zcalc, znonz, ztot)

        call ibc_from_abc(abc, ibc, nfun2, nfun3)

        ! scatter integrals onto disk
        ioff = 1
        do i1 = 1, nfun1
          do i3 = 1, nbatch3_out
            ii3 = i3lst_out(i3) ! actual index in full basis3
            ! index 2 is fast and should be fast on disk -- we therefore write a block
            ! offset in q(iabc) compressed, but not in ioffr!!

            offr = dble(infg2(ig_coff, igrp_2_first)) &
               & +dble(nfun2)*dble(ii3-1+nfun3*(i1-1))
            call writem_big(ibc(ioff:ioff+nbatch2-1), nbatch2, &
                            ifil_half, irec_half, offr, '3IDXSCR')
            ioff = ioff+nbatch2
          end do
        end do
      end do
    end do

    ! terminate integral batcher
    call basis_3idx_batchs_term2
    ! clean up
    call memory_release(ibase)

    !
    !cccccccccccccccccccccccccc SECOND TRANSFORMATION ccccccccccccccccccccccccccc
    !
    ! start memory accounting again
    ibase => memory_allocate_integer(1)

    ! decide how big the batches can be (again)
    memory_available = memory_remaining()-10000

    ! treat occupied orbitals one at a time, and assume we can read all b at once
    nfun_batch_3 = min(nfun3, memory_available/(nfun2+nfun2))
    if (nfun_batch_3 .lt. 1) &
      call error('Not enough memory for 2nd transfmation', 'ao_3idx')
    ibc => memory_allocate(nfun2*nfun_batch_3)

    if (ipr .gt. 0) &
      WRITE (6, 2) 'Number of batchs in 2nd tran', nfun1*nfun3/nfun_batch_3

    ! second transformation
    nwritten = 0
    do i1 = 1, nfun1
      do i3 = 1, nfun3, nfun_batch_3
        n3read = min(nfun_batch_3, nfun3-i3+1)
        n23 = nfun2*n3read
        offr = dble(nfun2)*dble(i3-1+nfun3*(i1-1))
        call readm_big(ibc, n23, ifil_half, irec_half, offr, string)
        ioff = 1
        do i2 = 1, n3read
          call fmove(ibc(ioff:nfun2+ioff-1), int_fxx(i2, :, i1), nfun2)
          ioff = ioff+nfun2
        end do
        nwritten = nwritten+nfun2*n3read
      end do
    end do

    ! trash half-transformed integrals and clean up
    call deletem(ifil_half, irec_half, 0)
    call memory_release(ibase)

    if (ipr .gt. 0) then
      WRITE (6, 2) 'Number of 3-idx integrals', nwritten
      WRITE (6, 4) 'Time', second_x()-tim
      WRITE (6, *)
    end if

    return

  contains

    subroutine ibc_from_abc(abc, ibc, nfun2, nfun3)
      use iso_fortran_env, only: dp => real64
      implicit none
      integer, intent(in) :: nfun2, nfun3
      real(dp), intent(in) :: abc(nfun2, nfun2, nfun3)
      real(dp), intent(out) :: ibc(nfun2, nfun3, nfun2)

      integer :: i

      ibc = 0d0

      do i = 1, nfun2
        ibc(i, :, :) = transpose(abc(i, :, :))
      end do

    end subroutine ibc_from_abc

  end subroutine ao_3idx

  ! Calculate 3-index integrals (f_p|\chi_{\kappa} \chi_{\lambda})
  ! Of a fitbasis {f} with all possible binary products of the AO basis (\chi)
  !
  ! [in] vbas Name of the Fitbasis, i.e. ri_basis or oep_basis
  ! [in] fit  Nirmalisation of integrals, either 'J' for Coulumb or 'S' for overlap
  ! [in] naux number of auxbas functions
  ! [out] int_fkl 3-index integrals
  subroutine build_fkl(vbas, fit, naux, int_fkl)
    use common_cbas, only: ntg
    use common_cbascode, only: ig_ndim
    use memory, only: memory_allocate_integer, memory_allocate, memory_release
    use iso_fortran_env, only: dp => real64
    implicit none

    logical :: sph_ao, sph_au
    integer, intent(in) :: naux
    real(dp), intent(out) :: int_fkl(naux, ntg, ntg)
    character(len=*) :: vbas
    character(len=*) :: fit
    integer, pointer, dimension(:) :: ibase

    integer, pointer, dimension(:) :: nfg_au
    real(dp), pointer, dimension(:) :: exp_au, cgr_au
    real(dp), pointer, dimension(:, :, :) :: tmp
    !
    ! Related to AO basis set
    !
    ! Variables describing the AO basis Set
    ! ngrp_ao  number of groups
    ! nshl_ao  total number of distinct exponents on dostinct centres
    ! lenc     total number of contraction coefficients
    ! nprim_ao total number of primitives in AO basis
    ! sph_ao logical, Spherical (T) or Cartesian(F)
    !
    integer :: ngrp_ao, nshl_ao, lenc
    integer :: nprim_ao
    integer :: ngrp_au, nshl_au, nfit
    !
    ! Variables Containing the AO basis Set
    ! nfg_ao Group information
    ! exp_ao Exponents
    ! cgr_ao contraction coefficients
    !
    integer, pointer, dimension(:) :: nfg_ao
    real(dp), pointer, dimension(:) :: exp_ao, cgr_ao

    ibase => memory_allocate_integer(1)

    !... load AO basis
    call basis_ngroup('ORBITAL', ngrp_ao)
    call basis_nprim_sh('ORBITAL', nshl_ao)
    call basis_ncc('ORBITAL', lenc)
    call basis_nprim('ORBITAL', nprim_ao)
    nfg_ao => memory_allocate_integer(ig_ndim*ngrp_ao)
    exp_ao => memory_allocate(nshl_ao)
    cgr_ao => memory_allocate(lenc)
    call basis_all_info('ORBITAL', nfg_ao, exp_ao, cgr_ao, sph_ao)

    !... load aux basis
    call basis_ngroup(vbas, ngrp_au)
    call basis_nprim_sh(vbas, nshl_au)
    call basis_number(vbas, nfit)
    call basis_ncc(vbas, lenc)
    nfg_au => memory_allocate_integer(ig_ndim*ngrp_au)
    exp_au => memory_allocate(nshl_au)
    cgr_au => memory_allocate(lenc)
    call basis_all_info(vbas, nfg_au, exp_au, cgr_au, sph_au)

    !... get the ints
    tmp => memory_allocate(naux, ntg, ntg)
    tmp = 0d0
    call ao_3idx( &
      ntg, ngrp_ao, nfg_ao, exp_ao, cgr_ao, sph_ao, &
      ntg, ngrp_ao, nfg_ao, exp_ao, cgr_ao, sph_ao, &
      nfit, ngrp_au, nfg_au, exp_au, cgr_au, sph_au, &
      1, fit, int_fkl)

    call memory_release(ibase)
    return
  end subroutine build_fkl

end module ksinv_integrals
