!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Periodic ECD AO-operator ownership and construction helpers.
! **************************************************************************************************
MODULE qs_ecd_operators
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_current_nl_path
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_moments,                      ONLY: build_local_moments_der_matrix_rs_img
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix_rs_img
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_operators'

   TYPE, PUBLIC :: ecd_periodic_local_operators_type
      INTEGER, PRIVATE                                      :: moment_order = -1
      INTEGER, DIMENSION(:, :, :), POINTER, PRIVATE         :: cell_to_index => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER, PRIVATE :: momentum_rs => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER, &
         PRIVATE                                            :: moments_der_rs => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER, PRIVATE :: moments_rs => NULL()
   END TYPE ecd_periodic_local_operators_type

   TYPE, PUBLIC :: ecd_periodic_gth_operators_type
      INTEGER, PRIVATE                                      :: max_order = -1
      INTEGER, DIMENSION(:, :, :), POINTER, PRIVATE         :: cell_to_index => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER, &
         PRIVATE                                            :: current_path_rs => NULL()
   END TYPE ecd_periodic_gth_operators_type

   PUBLIC :: ecd_periodic_local_operators_build, &
             ecd_periodic_local_operators_get, &
             ecd_periodic_local_operators_release, &
             ecd_periodic_gth_operators_build, &
             ecd_periodic_gth_operators_get, &
             ecd_periodic_gth_operators_release, &
             ecd_project_real_operator_to_mo

CONTAINS

! **************************************************************************************************
!> \brief Project one real sparse AO operator into the full real MO basis.
!> \param operator_ao AO-space sparse operator
!> \param mo_coeff Real MO coefficient matrix
!> \param nao Number of AOs
!> \param nmo Number of MOs
!> \param operator_mo Full MO-space operator
! **************************************************************************************************
   SUBROUTINE ecd_project_real_operator_to_mo(operator_ao, mo_coeff, nao, nmo, operator_mo)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: operator_ao
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER, INTENT(IN)                                :: nao, nmo
      TYPE(cp_fm_type), INTENT(INOUT)                    :: operator_mo

      TYPE(cp_fm_struct_type), POINTER                   :: coeff_struct, mo_struct
      TYPE(cp_fm_type)                                   :: operator_coeff

      NULLIFY (coeff_struct, mo_struct)

      CALL cp_fm_get_info(mo_coeff, matrix_struct=coeff_struct)
      CALL cp_fm_create(operator_coeff, coeff_struct, name="ECD operator C", set_zero=.TRUE.)
      CALL cp_fm_struct_create(fmstruct=mo_struct, template_fmstruct=coeff_struct, &
                               nrow_global=nmo, ncol_global=nmo)
      CALL cp_fm_create(operator_mo, mo_struct, name="ECD C^T operator C", set_zero=.TRUE.)
      CALL cp_fm_struct_release(mo_struct)

      CALL cp_dbcsr_sm_fm_multiply(operator_ao%matrix, mo_coeff, operator_coeff, ncol=nmo, &
                                   alpha=1.0_dp, beta=0.0_dp)
      CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mo_coeff, operator_coeff, &
                         0.0_dp, operator_mo)

      CALL cp_fm_release(operator_coeff)

   END SUBROUTINE ecd_project_real_operator_to_mo

! **************************************************************************************************
!> \brief Build image-resolved local-current primitives for periodic ECD.
!> \param operators ECD-owned image-resolved AO operators
!> \param qs_env Quickstep environment with an explicit k-point layout
!> \param ref_point Reference point for every local multipole factor
!> \param moment_order Highest local multipole order required by the Taylor path
! **************************************************************************************************
   SUBROUTINE ecd_periodic_local_operators_build(operators, qs_env, ref_point, moment_order)
      TYPE(ecd_periodic_local_operators_type), &
         INTENT(INOUT)                                   :: operators
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ref_point
      INTEGER, INTENT(IN)                                :: moment_order
      INTEGER, DIMENSION(:, :, :), POINTER               :: moments_cell_to_index

      IF (moment_order < 1) THEN
         CPABORT("ECD| Periodic local-current moments require a positive Taylor order")
      END IF

      IF (ASSOCIATED(operators%cell_to_index) .OR. ASSOCIATED(operators%momentum_rs) .OR. &
          ASSOCIATED(operators%moments_der_rs) .OR. ASSOCIATED(operators%moments_rs)) THEN
         CPABORT("ECD| Periodic local-current operators are already allocated")
      END IF

      NULLIFY (moments_cell_to_index)
      CALL build_lin_mom_matrix_rs_img(qs_env, operators%momentum_rs, operators%cell_to_index)
      CALL build_local_moments_der_matrix_rs_img(qs_env, operators%moments_der_rs, moment_order, ref_point, &
                                                 moments_rs=operators%moments_rs, &
                                                 cell_to_index_out=moments_cell_to_index)
      IF (ANY(LBOUND(operators%cell_to_index) /= LBOUND(moments_cell_to_index)) .OR. &
          ANY(UBOUND(operators%cell_to_index) /= UBOUND(moments_cell_to_index))) THEN
         CPABORT("ECD| Local momentum and multipole image-cell maps are inconsistent")
      END IF
      IF (ANY(operators%cell_to_index /= moments_cell_to_index)) THEN
         CPABORT("ECD| Local momentum and multipole image-cell maps are inconsistent")
      END IF
      DEALLOCATE (moments_cell_to_index)
      operators%moment_order = moment_order

   END SUBROUTINE ecd_periodic_local_operators_build

! **************************************************************************************************
!> \brief Return ECD-owned image-resolved local-current primitives without transferring ownership.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_local_operators_get(operators, cell_to_index, momentum_rs, moments_der_rs, moments_rs)
      TYPE(ecd_periodic_local_operators_type), &
         INTENT(IN), TARGET                              :: operators
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: momentum_rs
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: moments_der_rs
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                          :: moments_rs

      IF (.NOT. ASSOCIATED(operators%cell_to_index) .OR. .NOT. ASSOCIATED(operators%momentum_rs) .OR. &
          .NOT. ASSOCIATED(operators%moments_der_rs)) THEN
         CPABORT("ECD| Periodic local-current operators are unavailable")
      END IF
      cell_to_index => operators%cell_to_index
      momentum_rs => operators%momentum_rs
      moments_der_rs => operators%moments_der_rs
      IF (PRESENT(moments_rs)) THEN
         IF (.NOT. ASSOCIATED(operators%moments_rs)) THEN
            CPABORT("ECD| Periodic local coordinate moments are unavailable")
         END IF
         moments_rs => operators%moments_rs
      END IF

   END SUBROUTINE ecd_periodic_local_operators_get

! **************************************************************************************************
!> \brief Release ECD-owned image-resolved local-current primitives.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_local_operators_release(operators)
      TYPE(ecd_periodic_local_operators_type), &
         INTENT(INOUT)                                   :: operators

      IF (ASSOCIATED(operators%cell_to_index)) THEN
         DEALLOCATE (operators%cell_to_index)
         NULLIFY (operators%cell_to_index)
      END IF
      IF (ASSOCIATED(operators%momentum_rs)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%momentum_rs)
         NULLIFY (operators%momentum_rs)
      END IF
      IF (ASSOCIATED(operators%moments_der_rs)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%moments_der_rs)
         NULLIFY (operators%moments_der_rs)
      END IF
      IF (ASSOCIATED(operators%moments_rs)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%moments_rs)
         NULLIFY (operators%moments_rs)
      END IF
      operators%moment_order = -1

   END SUBROUTINE ecd_periodic_local_operators_release

! **************************************************************************************************
!> \brief Build real image-resolved Pickard--Mauri GTH path kernels.
!> \param operators ECD-owned GTH raw path kernels
!> \param qs_env Quickstep environment with explicit-Gamma image matrices
!> \param direction Unit propagation direction
!> \param ref_point Common position and path reference point
!> \param max_order Highest normalized spatial Taylor order
! **************************************************************************************************
   SUBROUTINE ecd_periodic_gth_operators_build(operators, qs_env, direction, ref_point, max_order)
      TYPE(ecd_periodic_gth_operators_type), &
         INTENT(INOUT)                                   :: operators
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: direction, ref_point
      INTEGER, INTENT(IN)                                :: max_order

      INTEGER                                            :: component, iimage, nimages, order
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints_all
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      IF (max_order < 0) THEN
         CPABORT("ECD| Periodic GTH path moments require a non-negative Taylor order")
      END IF
      IF (ASSOCIATED(operators%cell_to_index) .OR. ASSOCIATED(operators%current_path_rs)) THEN
         CPABORT("ECD| Periodic GTH path operators are already allocated")
      END IF

      NULLIFY (cell, cell_to_index, dft_control, kpoints_all, matrix_ks_kp, para_env, particle_set, &
               qs_kind_set, sab_all, sap_ppnl)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_ks_kp=matrix_ks_kp, &
                      para_env=para_env, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      sab_all=sab_all, sap_ppnl=sap_ppnl)
      IF (.NOT. ASSOCIATED(cell) .OR. .NOT. ASSOCIATED(dft_control) .OR. .NOT. ASSOCIATED(matrix_ks_kp) .OR. &
          .NOT. ASSOCIATED(para_env) .OR. .NOT. ASSOCIATED(particle_set) .OR. .NOT. ASSOCIATED(qs_kind_set) .OR. &
          .NOT. ASSOCIATED(sab_all)) THEN
         CPABORT("ECD| Explicit-Gamma data needed for GTH path moments are unavailable")
      END IF

      CALL kpoint_create(kpoints_all)
      CALL kpoint_init_cell_index(kpoints_all, sab_all, para_env, nimages)
      CALL get_kpoint_info(kpoints_all, cell_to_index=cell_to_index)
      ALLOCATE (operators%cell_to_index(LBOUND(cell_to_index, 1):UBOUND(cell_to_index, 1), &
                                        LBOUND(cell_to_index, 2):UBOUND(cell_to_index, 2), &
                                        LBOUND(cell_to_index, 3):UBOUND(cell_to_index, 3)))
      operators%cell_to_index = cell_to_index
      ALLOCATE (operators%current_path_rs(3, max_order + 1, nimages))
      DO iimage = 1, nimages
         DO order = 0, max_order
            DO component = 1, 3
               NULLIFY (operators%current_path_rs(component, order + 1, iimage)%matrix)
               ALLOCATE (operators%current_path_rs(component, order + 1, iimage)%matrix)
               CALL dbcsr_create(operators%current_path_rs(component, order + 1, iimage)%matrix, &
                                 template=matrix_ks_kp(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry, &
                                 name="ECD GTH PATH CURRENT")
               CALL cp_dbcsr_alloc_block_from_nbl(operators%current_path_rs(component, order + 1, iimage)%matrix, &
                                                  sab_all)
               CALL dbcsr_set(operators%current_path_rs(component, order + 1, iimage)%matrix, 0.0_dp)
            END DO
         END DO
      END DO

      CALL build_com_current_nl_path(qs_kind_set, sab_all, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                                     particle_set, cell, cell_to_index, direction, ref_point, max_order, &
                                     operators%current_path_rs)
      CALL kpoint_release(kpoints_all)
      operators%max_order = max_order

   END SUBROUTINE ecd_periodic_gth_operators_build

! **************************************************************************************************
!> \brief Return ECD-owned GTH path kernels without transferring ownership.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_gth_operators_get(operators, cell_to_index, current_path_rs)
      TYPE(ecd_periodic_gth_operators_type), &
         INTENT(IN), TARGET                              :: operators
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: current_path_rs

      IF (.NOT. ASSOCIATED(operators%cell_to_index) .OR. .NOT. ASSOCIATED(operators%current_path_rs)) THEN
         CPABORT("ECD| Periodic GTH path operators are unavailable")
      END IF
      cell_to_index => operators%cell_to_index
      current_path_rs => operators%current_path_rs

   END SUBROUTINE ecd_periodic_gth_operators_get

! **************************************************************************************************
!> \brief Release ECD-owned GTH path kernels.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_gth_operators_release(operators)
      TYPE(ecd_periodic_gth_operators_type), &
         INTENT(INOUT)                                   :: operators

      IF (ASSOCIATED(operators%cell_to_index)) THEN
         DEALLOCATE (operators%cell_to_index)
         NULLIFY (operators%cell_to_index)
      END IF
      IF (ASSOCIATED(operators%current_path_rs)) THEN
         CALL dbcsr_deallocate_matrix_set(operators%current_path_rs)
         NULLIFY (operators%current_path_rs)
      END IF
      operators%max_order = -1

   END SUBROUTINE ecd_periodic_gth_operators_release

END MODULE qs_ecd_operators
