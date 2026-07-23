!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Folded implicit-Gamma zero-order operators for periodic ECD.
! **************************************************************************************************
MODULE qs_ecd_gamma_operators
   USE commutator_rpnl,                 ONLY: build_com_mom_nl
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE kinds,                           ONLY: dp
   USE qs_ecd_gamma_data,               ONLY: ecd_gamma_context_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_moments,                      ONLY: build_local_magmom_matrix,&
                                              build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_gamma_operators'
   INTEGER, PARAMETER, PRIVATE          :: ncart = 3

   TYPE, PUBLIC :: ecd_gamma_zero_order_operators_type
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: position => NULL(), local_current => NULL(), &
                                                    static_magnetic_diagnostic => NULL()
   END TYPE ecd_gamma_zero_order_operators_type

   PUBLIC :: ecd_gamma_zero_order_operators_build, &
             ecd_gamma_zero_order_operators_get, &
             ecd_gamma_zero_order_operators_release

CONTAINS

! **************************************************************************************************
!> \brief Build folded local and diagnostic matrices at implicit Gamma.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_zero_order_operators_build(operators, qs_env, context)
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(INOUT) :: operators
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context

      INTEGER                                            :: idir
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gth_current_symmetric, gth_magnetic

      IF (ASSOCIATED(operators%position) .OR. ASSOCIATED(operators%local_current) .OR. &
          ASSOCIATED(operators%static_magnetic_diagnostic)) THEN
         CPABORT("ECD| Implicit-Gamma zero-order operators are already allocated")
      END IF
      IF (.NOT. ASSOCIATED(context%matrix_s) .OR. .NOT. ASSOCIATED(context%sab_orb) .OR. &
          .NOT. ASSOCIATED(context%sab_all)) THEN
         CPABORT("ECD| Implicit-Gamma zero-order operator context is incomplete")
      END IF

      NULLIFY (gth_current_symmetric, gth_magnetic)
      CALL ecd_gamma_ao_vector_create(operators%position, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD implicit Gamma position', dbcsr_type_symmetric)
      CALL ecd_gamma_ao_vector_create(operators%local_current, context%matrix_s(1)%matrix, context%sab_all, &
                                      'ECD implicit Gamma local current')
      CALL ecd_gamma_ao_vector_create(operators%static_magnetic_diagnostic, context%matrix_s(1)%matrix, &
                                      context%sab_orb, 'ECD implicit Gamma static magnetic diagnostic', &
                                      dbcsr_type_antisymmetric)

      CALL build_local_moment_matrix(qs_env, operators%position, 1, ref_point=context%ref_point)
      CALL build_lin_mom_matrix(qs_env, operators%local_current)
      CALL build_local_magmom_matrix(qs_env, operators%static_magnetic_diagnostic, 1, ref_point=context%ref_point)
      IF (ASSOCIATED(context%sap_ppnl)) THEN
         CALL ecd_gamma_ao_vector_create(gth_current_symmetric, context%matrix_s(1)%matrix, context%sab_orb, &
                                         'ECD implicit Gamma GTH current symmetric', dbcsr_type_antisymmetric)
         CALL ecd_gamma_ao_vector_create(gth_magnetic, context%matrix_s(1)%matrix, context%sab_orb, &
                                         'ECD implicit Gamma static magnetic diagnostic', dbcsr_type_antisymmetric)
         CALL build_com_mom_nl(context%qs_kind_set, context%sab_orb, context%sap_ppnl, &
                               context%dft_control%qs_control%eps_ppnl, context%particle_set, context%cell, &
                               matrix_rv=gth_current_symmetric, matrix_rxrv=gth_magnetic)
         CALL ecd_gamma_gth_magnetic_reference_translate(gth_magnetic, gth_current_symmetric, context%ref_point)
         DO idir = 1, ncart
            CALL dbcsr_add(operators%static_magnetic_diagnostic(idir)%matrix, gth_magnetic(idir)%matrix, &
                           1.0_dp, 1.0_dp)
         END DO
         CALL dbcsr_deallocate_matrix_set(gth_current_symmetric)
         CALL dbcsr_deallocate_matrix_set(gth_magnetic)
      END IF

   END SUBROUTINE ecd_gamma_zero_order_operators_build

! **************************************************************************************************
!> \brief Borrow zero-order folded local and diagnostic matrices without transferring ownership.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_zero_order_operators_get(operators, position, local_current, static_magnetic_diagnostic)
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: operators
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER, OPTIONAL :: position, local_current
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER, OPTIONAL :: static_magnetic_diagnostic

      IF (.NOT. ASSOCIATED(operators%position) .OR. .NOT. ASSOCIATED(operators%local_current) .OR. &
          .NOT. ASSOCIATED(operators%static_magnetic_diagnostic)) THEN
         CPABORT("ECD| Implicit-Gamma zero-order operators are unavailable")
      END IF
      IF (PRESENT(position)) position => operators%position
      IF (PRESENT(local_current)) local_current => operators%local_current
      IF (PRESENT(static_magnetic_diagnostic)) static_magnetic_diagnostic => operators%static_magnetic_diagnostic

   END SUBROUTINE ecd_gamma_zero_order_operators_get

! **************************************************************************************************
!> \brief Release ECD-owned folded zero-order operator matrices.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_zero_order_operators_release(operators)
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(INOUT) :: operators

      IF (ASSOCIATED(operators%position)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%position)
         NULLIFY (operators%position)
      END IF
      IF (ASSOCIATED(operators%local_current)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%local_current)
         NULLIFY (operators%local_current)
      END IF
      IF (ASSOCIATED(operators%static_magnetic_diagnostic)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%static_magnetic_diagnostic)
         NULLIFY (operators%static_magnetic_diagnostic)
      END IF
   END SUBROUTINE ecd_gamma_zero_order_operators_release

! **************************************************************************************************
!> \brief Shift an origin-defined GTH static magnetic diagnostic to the explicit reference point.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_gth_magnetic_reference_translate(magnetic, current, ref_point)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: current, magnetic
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: ref_point

      IF (SIZE(magnetic) /= ncart .OR. SIZE(current) /= ncart) THEN
         CPABORT('ECD| GTH static magnetic reference translation received incompatible Cartesian storage')
      END IF
      CALL dbcsr_add(magnetic(1)%matrix, current(3)%matrix, 1.0_dp, -ref_point(2))
      CALL dbcsr_add(magnetic(1)%matrix, current(2)%matrix, 1.0_dp, ref_point(3))
      CALL dbcsr_add(magnetic(2)%matrix, current(1)%matrix, 1.0_dp, -ref_point(3))
      CALL dbcsr_add(magnetic(2)%matrix, current(3)%matrix, 1.0_dp, ref_point(1))
      CALL dbcsr_add(magnetic(3)%matrix, current(2)%matrix, 1.0_dp, -ref_point(1))
      CALL dbcsr_add(magnetic(3)%matrix, current(1)%matrix, 1.0_dp, ref_point(2))

   END SUBROUTINE ecd_gamma_gth_magnetic_reference_translate

! **************************************************************************************************
!> \brief Allocate a Cartesian folded AO vector owned by the implicit-Gamma provider.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_ao_vector_create(matrices, template, sab_nl, name, matrix_type)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices
      TYPE(dbcsr_type), POINTER                          :: template
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_nl
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL             :: matrix_type

      INTEGER                                            :: idir
      CHARACTER(LEN=1)                                   :: my_matrix_type

      my_matrix_type = dbcsr_type_no_symmetry
      IF (PRESENT(matrix_type)) my_matrix_type = matrix_type
      CALL dbcsr_allocate_matrix_set(matrices, ncart)
      DO idir = 1, ncart
         CALL dbcsr_init_p(matrices(idir)%matrix)
         CALL dbcsr_create(matrices(idir)%matrix, template=template, matrix_type=my_matrix_type, name=name)
         CALL cp_dbcsr_alloc_block_from_nbl(matrices(idir)%matrix, sab_nl)
         CALL dbcsr_set(matrices(idir)%matrix, 0.0_dp)
      END DO

   END SUBROUTINE ecd_gamma_ao_vector_create

END MODULE qs_ecd_gamma_operators
