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
                                              dbcsr_desymmetrize,&
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
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_gamma_operators'
   INTEGER, PARAMETER, PRIVATE          :: ncart = 3

   TYPE, PUBLIC :: ecd_gamma_zero_order_operators_type
      LOGICAL                                      :: gth_available = .FALSE.
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: position => NULL(), local_current => NULL(), &
                                                    gth_current => NULL(), complete_current => NULL()
   END TYPE ecd_gamma_zero_order_operators_type

   PUBLIC :: ecd_gamma_zero_order_operators_build, &
             ecd_gamma_zero_order_operators_get, &
             ecd_gamma_zero_order_operators_release

CONTAINS

! **************************************************************************************************
!> \brief Build folded length and complete-current matrices at implicit Gamma.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_zero_order_operators_build(operators, qs_env, context, pseudopotential_commutator)
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(INOUT) :: operators
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator

      INTEGER                                            :: idir
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gth_current_symmetric

      IF (ASSOCIATED(operators%position) .OR. ASSOCIATED(operators%local_current) .OR. &
          ASSOCIATED(operators%gth_current) .OR. ASSOCIATED(operators%complete_current)) THEN
         CPABORT("ECD| Implicit-Gamma zero-order operators are already allocated")
      END IF
      IF (.NOT. ASSOCIATED(context%matrix_s) .OR. .NOT. ASSOCIATED(context%sab_orb) .OR. &
          .NOT. ASSOCIATED(context%sab_all)) THEN
         CPABORT("ECD| Implicit-Gamma zero-order operator context is incomplete")
      END IF

      NULLIFY (gth_current_symmetric)
      CALL ecd_gamma_ao_vector_create(operators%position, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD implicit Gamma position', dbcsr_type_symmetric)
      CALL ecd_gamma_ao_vector_create(operators%local_current, context%matrix_s(1)%matrix, context%sab_all, &
                                      'ECD implicit Gamma local current')
      CALL ecd_gamma_ao_vector_create(operators%gth_current, context%matrix_s(1)%matrix, context%sab_all, &
                                      'ECD implicit Gamma GTH current')
      CALL ecd_gamma_ao_vector_create(operators%complete_current, context%matrix_s(1)%matrix, context%sab_all, &
                                      'ECD implicit Gamma complete current')

      CALL build_local_moment_matrix(qs_env, operators%position, 1, ref_point=context%ref_point)
      CALL build_lin_mom_matrix(qs_env, operators%local_current)
      DO idir = 1, ncart
         CALL dbcsr_add(operators%complete_current(idir)%matrix, operators%local_current(idir)%matrix, &
                        0.0_dp, 1.0_dp)
      END DO

      IF (pseudopotential_commutator .AND. ASSOCIATED(context%sap_ppnl)) THEN
         CALL ecd_gamma_ao_vector_create(gth_current_symmetric, context%matrix_s(1)%matrix, context%sab_orb, &
                                         'ECD implicit Gamma GTH current symmetric', dbcsr_type_antisymmetric)
         CALL build_com_mom_nl(context%qs_kind_set, context%sab_orb, context%sap_ppnl, &
                               context%dft_control%qs_control%eps_ppnl, context%particle_set, context%cell, &
                               matrix_rv=gth_current_symmetric)
         DO idir = 1, ncart
            CALL dbcsr_desymmetrize(gth_current_symmetric(idir)%matrix, operators%gth_current(idir)%matrix)
            CALL dbcsr_add(operators%complete_current(idir)%matrix, operators%gth_current(idir)%matrix, &
                           1.0_dp, 1.0_dp)
         END DO
         CALL dbcsr_deallocate_matrix_set(gth_current_symmetric)
         operators%gth_available = .TRUE.
      END IF

   END SUBROUTINE ecd_gamma_zero_order_operators_build

! **************************************************************************************************
!> \brief Borrow zero-order folded operator matrices without transferring ownership.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_zero_order_operators_get(operators, position, local_current, gth_current, complete_current, &
                                                  gth_available)
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: operators
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: position, local_current, gth_current, complete_current
      LOGICAL, INTENT(OUT), OPTIONAL                     :: gth_available

      IF (.NOT. ASSOCIATED(operators%position) .OR. .NOT. ASSOCIATED(operators%local_current) .OR. &
          .NOT. ASSOCIATED(operators%gth_current) .OR. .NOT. ASSOCIATED(operators%complete_current)) THEN
         CPABORT("ECD| Implicit-Gamma zero-order operators are unavailable")
      END IF
      position => operators%position
      local_current => operators%local_current
      gth_current => operators%gth_current
      complete_current => operators%complete_current
      IF (PRESENT(gth_available)) gth_available = operators%gth_available

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
      IF (ASSOCIATED(operators%gth_current)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%gth_current)
         NULLIFY (operators%gth_current)
      END IF
      IF (ASSOCIATED(operators%complete_current)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%complete_current)
         NULLIFY (operators%complete_current)
      END IF
      operators%gth_available = .FALSE.

   END SUBROUTINE ecd_gamma_zero_order_operators_release

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
