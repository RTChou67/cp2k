!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program for atomistic simulations                                              !
!   Copyright 2000-2026 CP2K developers group                                                      !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Complete E1 and M1 operators for electronic circular dichroism.
! **************************************************************************************************
MODULE qs_ecd_operator_utils
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_mom_nl,&
                                              build_com_nl_mag
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_desymmetrize, dbcsr_finalize, dbcsr_init_p, dbcsr_p_type, &
        dbcsr_reserve_blocks, dbcsr_set, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_generate_filename,&
                                              cp_logger_type
   USE input_constants,                 ONLY: tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity
   USE kinds,                           ONLY: dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_local_magmom_matrix,&
                                              build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_ang_mom_matrix,&
                                              build_lin_mom_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PRIVATE       :: ncart = 3
   REAL(KIND=dp), PARAMETER, PRIVATE :: velocity_degeneracy_threshold = 1.0E-8_dp

   TYPE, PUBLIC :: ecd_mo_operators_type
      COMPLEX(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE :: e1, m1
   END TYPE ecd_mo_operators_type

   PUBLIC :: ecd_build_e1_m1_mo
   PUBLIC :: ecd_build_length_ao
   PUBLIC :: ecd_open_transition_files
   PUBLIC :: ecd_write_transition_operators

CONTAINS

! **************************************************************************************************
!> \brief Build complete E1 and M1 MO matrices for every spin channel.
!>
!> LENGTH uses direct position and magnetic AO operators.  VELOCITY constructs
!> the complete local plus non-local velocity internally, then uses it for both
!> the off-diagonal E1 reconstruction and M1.  Velocity matrices never leave
!> this routine.
! **************************************************************************************************
   SUBROUTINE ecd_build_e1_m1_mo(qs_env, gauge, pseudopotential_commutator, reference_point, &
                                 operators, reject_degenerate_velocity)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: gauge
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: reference_point
      TYPE(ecd_mo_operators_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: operators
      LOGICAL, INTENT(IN), OPTIONAL                      :: reject_degenerate_velocity

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: velocity_mo
      INTEGER                                            :: column, iatom, idir, inode, ispin, &
                                                            jatom, last_jatom, nao, &
                                                            ncommutator_block, nmo, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: commutator_cols, commutator_rows
      LOGICAL                                            :: include_nonlocal_velocity, &
                                                            reject_degenerate, use_velocity
      REAL(KIND=dp)                                      :: gap
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: commutator_matrix, magnetic_matrix, &
                                                            operator_matrix
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: coeff_struct
      TYPE(cp_fm_type)                                   :: operator_coeff, operator_mo
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: commutator_ao, magnetic_ao, matrix_s, &
                                                            operator_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (cell, coeff_struct, commutator_ao, dft_control, eigenvalues, magnetic_ao, matrix_s, mo_coeff, mos, &
               nl_iterator, operator_ao, particle_set, qs_kind_set, sab_orb, sap_ppnl)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_s=matrix_s, mos=mos, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set, sab_orb=sab_orb, sap_ppnl=sap_ppnl)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(mos))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(sab_orb))
      IF (gauge /= tddfpt_dipole_length .AND. gauge /= tddfpt_dipole_velocity) THEN
         CPABORT("ECD| E1/M1 construction received an unsupported gauge")
      END IF

      use_velocity = gauge == tddfpt_dipole_velocity
      reject_degenerate = .FALSE.
      IF (PRESENT(reject_degenerate_velocity)) reject_degenerate = reject_degenerate_velocity
      IF (use_velocity) THEN
         CALL dbcsr_allocate_matrix_set(operator_ao, ncart)
         DO idir = 1, ncart
            CALL dbcsr_init_p(operator_ao(idir)%matrix)
            IF (ALL(cell%perd == 0)) THEN
               CALL dbcsr_create(operator_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(operator_ao(idir)%matrix, sab_orb)
            ELSE
               CALL dbcsr_create(operator_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_desymmetrize(matrix_s(1)%matrix, operator_ao(idir)%matrix)
            END IF
            CALL dbcsr_set(operator_ao(idir)%matrix, 0.0_dp)
         END DO
         CALL build_lin_mom_matrix(qs_env, operator_ao)

         include_nonlocal_velocity = pseudopotential_commutator .AND. ASSOCIATED(sap_ppnl)
         IF (include_nonlocal_velocity) THEN
            CALL dbcsr_allocate_matrix_set(commutator_ao, ncart)
            DO idir = 1, ncart
               CALL dbcsr_init_p(commutator_ao(idir)%matrix)
               CALL dbcsr_create(commutator_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               IF (ALL(cell%perd == 0)) THEN
                  CALL cp_dbcsr_alloc_block_from_nbl(commutator_ao(idir)%matrix, sab_orb)
                  CALL dbcsr_set(commutator_ao(idir)%matrix, 0.0_dp)
               END IF
            END DO
            IF (ALL(cell%perd == 1)) THEN
               ncommutator_block = 0
               last_jatom = 0
               CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
               DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
                  CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, inode=inode)
                  IF (inode == 1) last_jatom = 0
                  IF (jatom == last_jatom) CYCLE
                  last_jatom = jatom
                  ncommutator_block = ncommutator_block + 1
               END DO
               CALL neighbor_list_iterator_release(nl_iterator)
               ALLOCATE (commutator_cols(MAX(1, ncommutator_block)), &
                         commutator_rows(MAX(1, ncommutator_block)))
               ncommutator_block = 0
               last_jatom = 0
               CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
               DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
                  CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, inode=inode)
                  IF (inode == 1) last_jatom = 0
                  IF (jatom == last_jatom) CYCLE
                  last_jatom = jatom
                  ncommutator_block = ncommutator_block + 1
                  commutator_rows(ncommutator_block) = MIN(iatom, jatom)
                  commutator_cols(ncommutator_block) = MAX(iatom, jatom)
               END DO
               CALL neighbor_list_iterator_release(nl_iterator)
               DO idir = 1, ncart
                  CALL dbcsr_finalize(commutator_ao(idir)%matrix)
                  CALL dbcsr_reserve_blocks(commutator_ao(idir)%matrix, &
                                            commutator_rows(1:ncommutator_block), &
                                            commutator_cols(1:ncommutator_block))
                  CALL dbcsr_finalize(commutator_ao(idir)%matrix)
                  CALL dbcsr_set(commutator_ao(idir)%matrix, 0.0_dp)
               END DO
               DEALLOCATE (commutator_cols, commutator_rows)
            END IF
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                                  particle_set, cell, matrix_rv=commutator_ao)
         END IF
      ELSE
         CALL ecd_build_length_ao(qs_env, pseudopotential_commutator, reference_point, operator_ao, magnetic_ao)
         include_nonlocal_velocity = .FALSE.
      END IF

      ALLOCATE (operators(SIZE(mos)))
      DO ispin = 1, SIZE(mos)
         NULLIFY (coeff_struct, eigenvalues, mo_coeff)
         CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, mo_coeff=mo_coeff, nmo=nmo)
         CPASSERT(ASSOCIATED(eigenvalues))
         CPASSERT(ASSOCIATED(mo_coeff))
         IF (SIZE(eigenvalues) /= nmo) CPABORT("ECD| E1/M1 eigenvalue and MO dimensions are inconsistent")
         CALL cp_fm_get_info(mo_coeff, matrix_struct=coeff_struct, nrow_global=nao, ncol_global=nmo)
         ALLOCATE (operators(ispin)%e1(ncart, nmo, nmo), operators(ispin)%m1(ncart, nmo, nmo), &
                   operator_matrix(nmo, nmo))
         operators(ispin)%e1 = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         operators(ispin)%m1 = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         IF (use_velocity) THEN
            ALLOCATE (velocity_mo(ncart, nmo, nmo))
            IF (include_nonlocal_velocity) ALLOCATE (commutator_matrix(nmo, nmo))
         ELSE
            ALLOCATE (magnetic_matrix(nmo, nmo))
         END IF

         CALL cp_fm_create(operator_coeff, coeff_struct, name="ECD operator C")
         CALL cp_fm_create(operator_mo, coeff_struct, nrow=nmo, ncol=nmo, name="ECD C^T operator C")
         DO idir = 1, ncart
            CALL cp_dbcsr_sm_fm_multiply(operator_ao(idir)%matrix, mo_coeff, operator_coeff, ncol=nmo, &
                                         alpha=1.0_dp, beta=0.0_dp)
            CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                               1.0_dp, mo_coeff, operator_coeff, 0.0_dp, operator_mo)
            CALL cp_fm_get_submatrix(operator_mo, operator_matrix)
            IF (use_velocity) THEN
               IF (include_nonlocal_velocity) THEN
                  CALL cp_dbcsr_sm_fm_multiply(commutator_ao(idir)%matrix, mo_coeff, operator_coeff, ncol=nmo, &
                                               alpha=1.0_dp, beta=0.0_dp)
                  CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                                     1.0_dp, mo_coeff, operator_coeff, 0.0_dp, operator_mo)
                  CALL cp_fm_get_submatrix(operator_mo, commutator_matrix)
                  operator_matrix = operator_matrix + commutator_matrix
               END IF
               velocity_mo(idir, :, :) = CMPLX(0.0_dp, -operator_matrix, KIND=dp)
            ELSE
               operators(ispin)%e1(idir, :, :) = CMPLX(-operator_matrix, 0.0_dp, KIND=dp)
               CALL cp_dbcsr_sm_fm_multiply(magnetic_ao(idir)%matrix, mo_coeff, operator_coeff, ncol=nmo, &
                                            alpha=1.0_dp, beta=0.0_dp)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                                  1.0_dp, mo_coeff, operator_coeff, 0.0_dp, operator_mo)
               CALL cp_fm_get_submatrix(operator_mo, magnetic_matrix)
               operators(ispin)%m1(idir, :, :) = CMPLX(0.0_dp, 0.5_dp*magnetic_matrix, KIND=dp)
            END IF
         END DO
         CALL cp_fm_release(operator_coeff)
         CALL cp_fm_release(operator_mo)

         IF (use_velocity) THEN
            DO column = 1, nmo
               DO row = 1, nmo
                  IF (row == column) CYCLE
                  gap = eigenvalues(row) - eigenvalues(column)
                  IF (ABS(gap) <= velocity_degeneracy_threshold) THEN
                     IF (reject_degenerate) CPABORT("ECD| VELOCITY gauge requires non-degenerate MO energies")
                     CYCLE
                  END IF
                  operators(ispin)%e1(:, row, column) = &
                     CMPLX(0.0_dp, 1.0_dp, KIND=dp)*velocity_mo(:, row, column)/gap
               END DO
            END DO
            operators(ispin)%m1(1, :, :) = 0.5_dp*( &
                                           MATMUL(operators(ispin)%e1(2, :, :), velocity_mo(3, :, :)) - &
                                           MATMUL(operators(ispin)%e1(3, :, :), velocity_mo(2, :, :)))
            operators(ispin)%m1(2, :, :) = 0.5_dp*( &
                                           MATMUL(operators(ispin)%e1(3, :, :), velocity_mo(1, :, :)) - &
                                           MATMUL(operators(ispin)%e1(1, :, :), velocity_mo(3, :, :)))
            operators(ispin)%m1(3, :, :) = 0.5_dp*( &
                                           MATMUL(operators(ispin)%e1(1, :, :), velocity_mo(2, :, :)) - &
                                           MATMUL(operators(ispin)%e1(2, :, :), velocity_mo(1, :, :)))
            DEALLOCATE (velocity_mo)
            IF (ALLOCATED(commutator_matrix)) DEALLOCATE (commutator_matrix)
         ELSE
            DEALLOCATE (magnetic_matrix)
         END IF
         DEALLOCATE (operator_matrix)
      END DO

      CALL dbcsr_deallocate_matrix_set(operator_ao)
      IF (ASSOCIATED(magnetic_ao)) CALL dbcsr_deallocate_matrix_set(magnetic_ao)
      IF (ASSOCIATED(commutator_ao)) CALL dbcsr_deallocate_matrix_set(commutator_ao)

   END SUBROUTINE ecd_build_e1_m1_mo

! **************************************************************************************************
!> \brief Build direct length-gauge E1 and M1 AO operators.
! **************************************************************************************************
   SUBROUTINE ecd_build_length_ao(qs_env, pseudopotential_commutator, reference_point, position_ao, magnetic_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: reference_point
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: position_ao, magnetic_ao

      INTEGER                                            :: idir
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_nl, magnetic_nl_full, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (cell, dft_control, magnetic_nl, magnetic_nl_full, matrix_s, particle_set, qs_kind_set, sab_all, sab_orb, &
               sap_ppnl)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_s=matrix_s, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set, sab_orb=sab_orb, &
                      sab_all=sab_all, sap_ppnl=sap_ppnl)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(sab_all))
      CPASSERT(ASSOCIATED(sab_orb))

      CALL dbcsr_allocate_matrix_set(position_ao, ncart)
      CALL dbcsr_allocate_matrix_set(magnetic_ao, ncart)
      DO idir = 1, ncart
         CALL dbcsr_init_p(position_ao(idir)%matrix)
         CALL dbcsr_create(position_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(position_ao(idir)%matrix, sab_orb)
         CALL dbcsr_set(position_ao(idir)%matrix, 0.0_dp)
         CALL dbcsr_init_p(magnetic_ao(idir)%matrix)
         IF (ALL(cell%perd == 1)) THEN
            CALL dbcsr_create(magnetic_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(magnetic_ao(idir)%matrix, sab_all)
         ELSE
            CALL dbcsr_create(magnetic_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(magnetic_ao(idir)%matrix, sab_orb)
         END IF
         CALL dbcsr_set(magnetic_ao(idir)%matrix, 0.0_dp)
      END DO

      CALL build_local_moment_matrix(qs_env, position_ao, 1, ref_point=reference_point)
      IF (ALL(cell%perd == 1)) THEN
         CALL build_ang_mom_matrix(qs_env, magnetic_ao, reference_point)
      ELSE
         CALL build_local_magmom_matrix(qs_env, magnetic_ao, 1, ref_point=reference_point)
      END IF

      IF (pseudopotential_commutator .AND. ASSOCIATED(sap_ppnl)) THEN
         CALL dbcsr_allocate_matrix_set(magnetic_nl, ncart)
         DO idir = 1, ncart
            CALL dbcsr_init_p(magnetic_nl(idir)%matrix)
            CALL dbcsr_create(magnetic_nl(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(magnetic_nl(idir)%matrix, sab_orb)
            CALL dbcsr_set(magnetic_nl(idir)%matrix, 0.0_dp)
         END DO
         IF (ALL(cell%perd == 1)) THEN
            CALL build_com_nl_mag(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                                  particle_set, magnetic_nl, reference_point, cell)
            CALL dbcsr_allocate_matrix_set(magnetic_nl_full, ncart)
            DO idir = 1, ncart
               CALL dbcsr_init_p(magnetic_nl_full(idir)%matrix)
               CALL dbcsr_create(magnetic_nl_full(idir)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_desymmetrize(magnetic_nl(idir)%matrix, magnetic_nl_full(idir)%matrix)
               CALL dbcsr_add(magnetic_ao(idir)%matrix, magnetic_nl_full(idir)%matrix, 1.0_dp, -1.0_dp)
            END DO
            CALL dbcsr_deallocate_matrix_set(magnetic_nl_full)
         ELSE
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                                  particle_set, cell, matrix_rxrv=magnetic_nl, ref_point=reference_point)
            DO idir = 1, ncart
               CALL dbcsr_add(magnetic_ao(idir)%matrix, magnetic_nl(idir)%matrix, 1.0_dp, 1.0_dp)
            END DO
         END IF
         CALL dbcsr_deallocate_matrix_set(magnetic_nl)
      END IF

   END SUBROUTINE ecd_build_length_ao

! **************************************************************************************************
!> \brief Open the shared transition-level E1 and M1 output files.
! **************************************************************************************************
   SUBROUTINE ecd_open_transition_files(logger, label, write_matrices, e1_unit, m1_unit, &
                                        e1_filename, m1_filename)
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      LOGICAL, INTENT(OUT)                               :: write_matrices
      INTEGER, INTENT(OUT)                               :: e1_unit, m1_unit
      CHARACTER(LEN=*), INTENT(INOUT)                    :: e1_filename, m1_filename

      CALL cp_logger_generate_filename(logger, e1_filename, logger%iter_info%project_name, &
                                       "-ecd-e1-transition-matrix.dat", local=.FALSE.)
      CALL cp_logger_generate_filename(logger, m1_filename, logger%iter_info%project_name, &
                                       "-ecd-m1-transition-matrix.dat", local=.FALSE.)
      write_matrices = logger%para_env%is_source()
      IF (.NOT. write_matrices) RETURN

      CALL open_file(e1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=e1_unit)
      CALL open_file(m1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=m1_unit)
      WRITE (e1_unit, '(A)') "# "//TRIM(label)//" E1 transition matrix"
      WRITE (m1_unit, '(A)') "# "//TRIM(label)//" M1 transition matrix"
      WRITE (e1_unit, '(A)') "# state energy_eV spin occupied virtual component e1_real e1_imag occupation_difference"
      WRITE (m1_unit, '(A)') "# state energy_eV spin occupied virtual component m1_real m1_imag occupation_difference"

   END SUBROUTINE ecd_open_transition_files

! **************************************************************************************************
!> \brief Write all Cartesian E1 and M1 components for one transition.
! **************************************************************************************************
   SUBROUTINE ecd_write_transition_operators(e1_unit, m1_unit, state, energy, ispin, iocc, ivirt, &
                                             occupation_difference, e1, m1)
      INTEGER, INTENT(IN)                                :: e1_unit, m1_unit, state
      REAL(KIND=dp), INTENT(IN)                          :: energy
      INTEGER, INTENT(IN)                                :: ispin, iocc, ivirt
      REAL(KIND=dp), INTENT(IN)                          :: occupation_difference
      COMPLEX(KIND=dp), DIMENSION(ncart), INTENT(IN)     :: e1, m1

      INTEGER                                            :: idir

      DO idir = 1, ncart
         WRITE (e1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
            state, energy*evolt, ispin, iocc, ivirt, idir, REAL(e1(idir), KIND=dp), AIMAG(e1(idir)), &
            occupation_difference
         WRITE (m1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
            state, energy*evolt, ispin, iocc, ivirt, idir, REAL(m1(idir), KIND=dp), AIMAG(m1(idir)), &
            occupation_difference
      END DO

   END SUBROUTINE ecd_write_transition_operators

END MODULE qs_ecd_operator_utils
