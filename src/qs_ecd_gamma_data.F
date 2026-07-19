!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Read-only implicit-Gamma data view for periodic ECD.
! **************************************************************************************************
MODULE qs_ecd_gamma_data
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_init_cell_index,&
                                              kpoint_initialize,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_ecd_subspace,                 ONLY: ecd_band_manifold_type,&
                                              ecd_band_manifolds_build,&
                                              ecd_band_manifolds_release
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_core_matrices,                ONLY: core_matrices,&
                                              kinetic_energy_matrix
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_gamma_data'

   TYPE, PUBLIC :: ecd_gamma_spin_type
      INTEGER                                      :: nao = 0, nmo = 0
      REAL(KIND=dp), DIMENSION(:), POINTER        :: eigenvalues => NULL(), eigenvalue_first => NULL(), &
                                                       eigenvalue_reference => NULL(), occupations => NULL()
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: band_rotation, coefficient_first
      TYPE(cp_fm_type), POINTER                    :: mo_coeff => NULL()
      TYPE(ecd_band_manifold_type), DIMENSION(:), &
         POINTER                                   :: manifolds => NULL()
   END TYPE ecd_gamma_spin_type

   TYPE, PUBLIC :: ecd_gamma_translation_layout_type
      INTEGER                                           :: ntranslations = 0
      INTEGER, DIMENSION(:, :, :), POINTER              :: cell_to_index => NULL()
      TYPE(kpoint_type), POINTER                        :: external_kpoints => NULL()
   END TYPE ecd_gamma_translation_layout_type

   TYPE, PUBLIC :: ecd_gamma_context_type
      INTEGER                                                  :: nspin = 0, spatial_order = 0
      REAL(KIND=dp), DIMENSION(3)                            :: direction = 0.0_dp, ref_point = 0.0_dp
      TYPE(cell_type), POINTER                                :: cell => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: matrix_h => NULL(), matrix_ks => NULL(), &
                                                                  matrix_s => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER            :: core_images => NULL(), &
                                                                  local_potential_images => NULL(), overlap_images => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER            :: current_taylor_imag => NULL(), &
                                                                  current_taylor_real => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: ks_first_imag => NULL(), ks_first_real => NULL(), &
                                                                  overlap_first_imag => NULL(), overlap_first_real => NULL()
      TYPE(dft_control_type), POINTER                         :: dft_control => NULL()
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: sab_orb => NULL(), sab_all => NULL(), &
                                                                  sap_ppnl => NULL()
      TYPE(particle_type), DIMENSION(:), POINTER              :: particle_set => NULL()
      TYPE(qs_kind_type), DIMENSION(:), POINTER               :: qs_kind_set => NULL()
      TYPE(ecd_gamma_translation_layout_type)                 :: translation_layout
      TYPE(ecd_gamma_spin_type), ALLOCATABLE, DIMENSION(:)    :: spin
   END TYPE ecd_gamma_context_type

   PUBLIC :: ecd_gamma_context_create_implicit, &
             ecd_gamma_context_core_images_get, &
             ecd_gamma_context_local_potential_images_get, &
             ecd_gamma_context_overlap_images_get, &
             ecd_gamma_context_overlap_images_status, &
             ecd_gamma_context_release

CONTAINS

! **************************************************************************************************
!> \brief Borrow field-free implicit-Gamma data from a converged Quickstep environment.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_create_implicit(context, qs_env, direction, ref_point, spatial_order, &
                                                eps_degeneracy)
      TYPE(ecd_gamma_context_type), INTENT(OUT)         :: context
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: direction, ref_point
      INTEGER, INTENT(IN)                                :: spatial_order
      REAL(KIND=dp), INTENT(IN)                          :: eps_degeneracy

      INTEGER                                            :: ispin, nao, nmo, nrow, ncol
      LOGICAL                                            :: do_kpoints
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupations
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(cp_blacs_env_type), POINTER                    :: blacs_env
      TYPE(mp_para_env_type), POINTER                     :: para_env
      TYPE(qs_ks_env_type), POINTER                       :: ks_env

      NULLIFY (context%core_images, context%current_taylor_imag, context%current_taylor_real, context%ks_first_imag, &
               context%ks_first_real, context%local_potential_images, &
               context%overlap_first_imag, context%overlap_first_real, &
               context%translation_layout%external_kpoints, context%translation_layout%cell_to_index, &
               context%overlap_images)
      context%translation_layout%ntranslations = 0
      NULLIFY (blacs_env, eigenvalues, ks_env, mo_coeff, mos, occupations, para_env)
      CALL get_qs_env(qs_env, cell=context%cell, do_kpoints=do_kpoints, mos=mos, &
                      ks_env=ks_env, matrix_h=context%matrix_h, matrix_ks=context%matrix_ks, &
                      matrix_s=context%matrix_s, &
                      dft_control=context%dft_control, sab_orb=context%sab_orb, sab_all=context%sab_all, &
                      blacs_env=blacs_env, para_env=para_env, sap_ppnl=context%sap_ppnl, &
                      particle_set=context%particle_set, &
                      qs_kind_set=context%qs_kind_set)

      IF (do_kpoints) THEN
         CPABORT("ECD| Implicit-Gamma data adapter requires the field-free Gamma SCF layout")
      END IF
      IF (.NOT. ASSOCIATED(context%cell) .OR. .NOT. ASSOCIATED(ks_env) .OR. .NOT. ASSOCIATED(mos) .OR. &
          .NOT. ASSOCIATED(context%matrix_h) .OR. .NOT. ASSOCIATED(context%matrix_ks) .OR. &
          .NOT. ASSOCIATED(context%matrix_s) .OR. &
          .NOT. ASSOCIATED(context%dft_control) .OR. .NOT. ASSOCIATED(context%sab_orb) .OR. &
          .NOT. ASSOCIATED(context%sab_all) .OR. .NOT. ASSOCIATED(context%particle_set) .OR. &
          .NOT. ASSOCIATED(context%qs_kind_set)) THEN
         CPABORT("ECD| Implicit-Gamma SCF data are incomplete")
      END IF
      IF (SIZE(mos) < 1 .OR. SIZE(context%matrix_ks) < 1 .OR. SIZE(context%matrix_s) < 1) THEN
         CPABORT("ECD| Implicit-Gamma SCF data contain an empty MO or AO matrix set")
      END IF
      IF (eps_degeneracy <= 0.0_dp) THEN
         CPABORT("ECD| EPS_DEGENERACY must be positive for the implicit-Gamma data adapter")
      END IF
      IF (spatial_order < 0) THEN
         CPABORT("ECD| SPATIAL_ORDER must be non-negative for the implicit-Gamma data adapter")
      END IF

      context%nspin = SIZE(mos)
      IF (SIZE(context%matrix_ks) /= context%nspin) THEN
         CPABORT("ECD| Implicit-Gamma KS spin channels are inconsistent with the MO channels")
      END IF
      context%direction = direction
      context%ref_point = ref_point
      context%spatial_order = spatial_order
      IF (spatial_order >= 1) THEN
         IF (.NOT. ASSOCIATED(blacs_env) .OR. .NOT. ASSOCIATED(para_env)) THEN
            CPABORT("ECD| Implicit-Gamma image snapshots require the parallel and BLACS environments")
         END IF
         CALL ecd_gamma_context_overlap_images_build(context, ks_env, mos, para_env, blacs_env)
         CALL ecd_gamma_context_core_images_build(context, qs_env)
      END IF
      ALLOCATE (context%spin(context%nspin))

      DO ispin = 1, context%nspin
         NULLIFY (eigenvalues, mo_coeff, occupations)
         CALL get_mo_set(mos(ispin), nao=nao, nmo=nmo, eigenvalues=eigenvalues, &
                         occupation_numbers=occupations, mo_coeff=mo_coeff)
         IF (nao < 1 .OR. nmo < 1 .OR. .NOT. ASSOCIATED(eigenvalues) .OR. &
             .NOT. ASSOCIATED(occupations) .OR. .NOT. ASSOCIATED(mo_coeff)) THEN
            CPABORT("ECD| Implicit-Gamma SCF did not provide complete molecular-orbital data")
         END IF
         IF (SIZE(eigenvalues) /= nmo .OR. SIZE(occupations) /= nmo) THEN
            CPABORT("ECD| Implicit-Gamma MO eigenvalue or occupation dimensions are inconsistent")
         END IF
         CALL cp_fm_get_info(mo_coeff, nrow_global=nrow, ncol_global=ncol)
         IF (nrow /= nao .OR. ncol /= nmo) THEN
            CPABORT("ECD| Implicit-Gamma MO coefficient dimensions are inconsistent")
         END IF

         context%spin(ispin)%nao = nao
         context%spin(ispin)%nmo = nmo
         context%spin(ispin)%eigenvalues => eigenvalues
         context%spin(ispin)%occupations => occupations
         context%spin(ispin)%mo_coeff => mo_coeff
         CALL ecd_band_manifolds_build(context%spin(ispin)%manifolds, eigenvalues, eps_degeneracy)
      END DO
   END SUBROUTINE ecd_gamma_context_create_implicit

! **************************************************************************************************
!> \brief Borrow the ECD-owned image-resolved core Hamiltonian and its cell map.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_core_images_get(context, core_images, cell_to_index)
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET   :: context
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: core_images
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      IF (.NOT. ASSOCIATED(context%core_images) .OR. .NOT. ASSOCIATED(context%translation_layout%cell_to_index)) THEN
         CPABORT("ECD| Implicit-Gamma core image snapshot is unavailable")
      END IF
      core_images => context%core_images
      cell_to_index => context%translation_layout%cell_to_index

   END SUBROUTINE ecd_gamma_context_core_images_get

! **************************************************************************************************
!> \brief Borrow frozen local-potential translation matrices and their shared cell map.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_local_potential_images_get(context, local_potential_images, cell_to_index)
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET   :: context
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: local_potential_images
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      IF (.NOT. ASSOCIATED(context%local_potential_images) .OR. &
          .NOT. ASSOCIATED(context%translation_layout%cell_to_index)) THEN
         CPABORT("ECD| Implicit-Gamma local-potential translation snapshot is unavailable")
      END IF
      IF (SIZE(context%local_potential_images, 1) /= context%nspin .OR. &
          SIZE(context%local_potential_images, 2) /= context%translation_layout%ntranslations) THEN
         CPABORT("ECD| Implicit-Gamma local-potential translation snapshot dimensions are inconsistent")
      END IF
      local_potential_images => context%local_potential_images
      cell_to_index => context%translation_layout%cell_to_index

   END SUBROUTINE ecd_gamma_context_local_potential_images_get

! **************************************************************************************************
!> \brief Borrow the ECD-owned image-resolved overlap matrices and their cell map.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_overlap_images_get(context, overlap_images, cell_to_index)
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET   :: context
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: overlap_images
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      IF (.NOT. ASSOCIATED(context%overlap_images) .OR. .NOT. ASSOCIATED(context%translation_layout%cell_to_index)) THEN
         CPABORT("ECD| Implicit-Gamma overlap image snapshot is unavailable")
      END IF
      overlap_images => context%overlap_images
      cell_to_index => context%translation_layout%cell_to_index

   END SUBROUTINE ecd_gamma_context_overlap_images_get

! **************************************************************************************************
!> \brief Report whether the implicit-Gamma context owns an image-resolved overlap snapshot.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_overlap_images_status(context, available, nimages)
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      LOGICAL, INTENT(OUT)                               :: available
      INTEGER, INTENT(OUT)                               :: nimages

      available = .FALSE.
      nimages = 0
      IF (.NOT. ASSOCIATED(context%overlap_images) .OR. .NOT. ASSOCIATED(context%translation_layout%cell_to_index)) THEN
         RETURN
      END IF
      IF (SIZE(context%overlap_images, 1) /= 1) THEN
         RETURN
      END IF
      nimages = context%translation_layout%ntranslations
      IF (nimages < 1 .OR. SIZE(context%overlap_images, 2) /= nimages) THEN
         RETURN
      END IF
      available = .TRUE.

   END SUBROUTINE ecd_gamma_context_overlap_images_status

! **************************************************************************************************
!> \brief Build a read-only image-resolved overlap snapshot for the implicit-Gamma ECD expansion.
!>
!> The standard Gamma SCF matrix is folded over lattice images.  The overlap
!> Taylor coefficients need the individual image blocks, which are generated
!> here without changing the matrices used by the SCF calculation.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_overlap_images_build(context, ks_env, mos, para_env, blacs_env)
      TYPE(ecd_gamma_context_type), INTENT(INOUT)        :: context
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      INTEGER                                            :: nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      IF (ASSOCIATED(context%overlap_images) .OR. ASSOCIATED(context%translation_layout%cell_to_index)) THEN
         CPABORT("ECD| Implicit-Gamma overlap image snapshot is already allocated")
      END IF
      IF (.NOT. ASSOCIATED(context%sab_orb)) THEN
         CPABORT("ECD| Implicit-Gamma overlap image snapshot requires the orbital neighbour list")
      END IF

      NULLIFY (cell_to_index)
      CALL kpoint_create(context%translation_layout%external_kpoints)
      context%translation_layout%external_kpoints%kp_scheme = "GAMMA"
      context%translation_layout%external_kpoints%use_real_wfn = .TRUE.
      CALL kpoint_initialize(context%translation_layout%external_kpoints, context%particle_set, context%cell)
      CALL kpoint_env_initialize(context%translation_layout%external_kpoints, para_env, blacs_env)
      CALL kpoint_initialize_mos(context%translation_layout%external_kpoints, mos)
      CALL kpoint_initialize_mo_set(context%translation_layout%external_kpoints)
      CALL kpoint_init_cell_index(context%translation_layout%external_kpoints, context%sab_orb, para_env, nimages)
      CALL get_kpoint_info(context%translation_layout%external_kpoints, cell_to_index=cell_to_index)
      IF (nimages < 1 .OR. .NOT. ASSOCIATED(cell_to_index)) THEN
         CALL kpoint_release(context%translation_layout%external_kpoints)
         CPABORT("ECD| Implicit-Gamma overlap image snapshot has no image cells")
      END IF

      context%translation_layout%ntranslations = nimages
      ALLOCATE (context%translation_layout%cell_to_index(LBOUND(cell_to_index, 1):UBOUND(cell_to_index, 1), &
                                                          LBOUND(cell_to_index, 2):UBOUND(cell_to_index, 2), &
                                                          LBOUND(cell_to_index, 3):UBOUND(cell_to_index, 3)))
      context%translation_layout%cell_to_index = cell_to_index
      CALL build_overlap_matrix(ks_env, matrixkp_s=context%overlap_images, &
                                matrix_name="ECD IMPLICIT GAMMA OVERLAP IMAGES", basis_type_a="ORB", &
                                basis_type_b="ORB", sab_nl=context%sab_orb, &
                                ext_kpoints=context%translation_layout%external_kpoints)

      IF (.NOT. ASSOCIATED(context%overlap_images) .OR. SIZE(context%overlap_images, 1) /= 1 .OR. &
          SIZE(context%overlap_images, 2) /= nimages .OR. context%translation_layout%ntranslations /= nimages) THEN
         CPABORT("ECD| Implicit-Gamma overlap image snapshot dimensions are inconsistent")
      END IF

   END SUBROUTINE ecd_gamma_context_overlap_images_build

! **************************************************************************************************
!> \brief Build the field-free core Hamiltonian in the overlap-image representation.
!>
!> The core matrix has an external image-kpoint construction path.  The density
!> argument is present for shared integral interfaces and remains unused because
!> this snapshot does not request forces.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_core_images_build(context, qs_env)
      TYPE(ecd_gamma_context_type), INTENT(INOUT)        :: context
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: iimage, ispin, nimages, ppl_method
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_dummy

      IF (ASSOCIATED(context%core_images)) THEN
         CPABORT("ECD| Implicit-Gamma core image snapshot is already allocated")
      END IF
      IF (.NOT. ASSOCIATED(context%translation_layout%external_kpoints) .OR. &
          .NOT. ASSOCIATED(context%overlap_images)) THEN
         CPABORT("ECD| Implicit-Gamma core image snapshot requires the overlap image environment")
      END IF
      nimages = context%translation_layout%ntranslations
      IF (nimages < 1 .OR. SIZE(context%overlap_images, 2) /= nimages) THEN
         CPABORT("ECD| Implicit-Gamma core image snapshot has an inconsistent translation layout")
      END IF
      IF (SIZE(context%matrix_h) /= 1) THEN
         CPABORT("ECD| Implicit-Gamma folded core matrix has an unsupported shape")
      END IF

      NULLIFY (density_dummy)
      CALL dbcsr_allocate_matrix_set(density_dummy, context%nspin, nimages)
      DO iimage = 1, nimages
         DO ispin = 1, context%nspin
            CALL dbcsr_init_p(density_dummy(ispin, iimage)%matrix)
            CALL dbcsr_create(density_dummy(ispin, iimage)%matrix, template=context%matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, name="ECD CORE SNAPSHOT DENSITY")
            CALL cp_dbcsr_alloc_block_from_nbl(density_dummy(ispin, iimage)%matrix, context%sab_orb)
            CALL dbcsr_set(density_dummy(ispin, iimage)%matrix, 0.0_dp)
         END DO
      END DO

      CALL kinetic_energy_matrix(qs_env, matrixkp_t=context%core_images, matrix_p=density_dummy, &
                                 ext_kpoints=context%translation_layout%external_kpoints, matrix_name="ECD CORE KINETIC", &
                                 calculate_forces=.FALSE., sab_orb=context%sab_orb)
      ppl_method = context%dft_control%qs_control%do_ppl_method
      CALL core_matrices(qs_env, context%core_images, density_dummy, calculate_forces=.FALSE., nder=0, &
                         ext_kpoints=context%translation_layout%external_kpoints)
      context%dft_control%qs_control%do_ppl_method = ppl_method
      CALL dbcsr_deallocate_matrix_set(density_dummy)

      IF (.NOT. ASSOCIATED(context%core_images) .OR. SIZE(context%core_images, 1) /= context%nspin .OR. &
          SIZE(context%core_images, 2) /= nimages) THEN
         CPABORT("ECD| Implicit-Gamma core image snapshot dimensions are inconsistent")
      END IF

   END SUBROUTINE ecd_gamma_context_core_images_build

! **************************************************************************************************
!> \brief Release ECD-owned metadata while preserving borrowed Quickstep objects.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_release(context)
      TYPE(ecd_gamma_context_type), INTENT(INOUT)       :: context

      INTEGER                                            :: ispin

      IF (ALLOCATED(context%spin)) THEN
         DO ispin = 1, SIZE(context%spin)
            CALL ecd_band_manifolds_release(context%spin(ispin)%manifolds)
            IF (ASSOCIATED(context%spin(ispin)%eigenvalue_first)) DEALLOCATE (context%spin(ispin)%eigenvalue_first)
            IF (ASSOCIATED(context%spin(ispin)%eigenvalue_reference)) THEN
               DEALLOCATE (context%spin(ispin)%eigenvalue_reference)
            END IF
            IF (ALLOCATED(context%spin(ispin)%band_rotation)) DEALLOCATE (context%spin(ispin)%band_rotation)
            IF (ALLOCATED(context%spin(ispin)%coefficient_first)) DEALLOCATE (context%spin(ispin)%coefficient_first)
            NULLIFY (context%spin(ispin)%eigenvalues, context%spin(ispin)%eigenvalue_first, &
                     context%spin(ispin)%eigenvalue_reference, context%spin(ispin)%occupations, &
                     context%spin(ispin)%mo_coeff)
         END DO
         DEALLOCATE (context%spin)
      END IF
      IF (ASSOCIATED(context%overlap_images)) THEN
         CALL dbcsr_deallocate_matrix_set(context%overlap_images)
         NULLIFY (context%overlap_images)
      END IF
      IF (ASSOCIATED(context%core_images)) THEN
         CALL dbcsr_deallocate_matrix_set(context%core_images)
         NULLIFY (context%core_images)
      END IF
      IF (ASSOCIATED(context%local_potential_images)) THEN
         CALL dbcsr_deallocate_matrix_set(context%local_potential_images)
         NULLIFY (context%local_potential_images)
      END IF
      IF (ASSOCIATED(context%current_taylor_imag)) THEN
         CALL dbcsr_deallocate_matrix_set(context%current_taylor_imag)
         NULLIFY (context%current_taylor_imag)
      END IF
      IF (ASSOCIATED(context%current_taylor_real)) THEN
         CALL dbcsr_deallocate_matrix_set(context%current_taylor_real)
         NULLIFY (context%current_taylor_real)
      END IF
      IF (ASSOCIATED(context%ks_first_real)) THEN
         CALL dbcsr_deallocate_matrix_set(context%ks_first_real)
         NULLIFY (context%ks_first_real)
      END IF
      IF (ASSOCIATED(context%ks_first_imag)) THEN
         CALL dbcsr_deallocate_matrix_set(context%ks_first_imag)
         NULLIFY (context%ks_first_imag)
      END IF
      IF (ASSOCIATED(context%overlap_first_real)) THEN
         CALL dbcsr_deallocate_matrix_set(context%overlap_first_real)
         NULLIFY (context%overlap_first_real)
      END IF
      IF (ASSOCIATED(context%overlap_first_imag)) THEN
         CALL dbcsr_deallocate_matrix_set(context%overlap_first_imag)
         NULLIFY (context%overlap_first_imag)
      END IF
      IF (ASSOCIATED(context%translation_layout%external_kpoints)) THEN
         CALL kpoint_release(context%translation_layout%external_kpoints)
      END IF
      IF (ASSOCIATED(context%translation_layout%cell_to_index)) THEN
         DEALLOCATE (context%translation_layout%cell_to_index)
         NULLIFY (context%translation_layout%cell_to_index)
      END IF
      context%translation_layout%ntranslations = 0
      NULLIFY (context%cell, context%matrix_h, context%matrix_ks, context%matrix_s, context%dft_control, &
               context%sab_orb, &
               context%sab_all, context%sap_ppnl, context%particle_set, context%qs_kind_set)
      context%nspin = 0
      context%spatial_order = 0
      context%direction = 0.0_dp
      context%ref_point = 0.0_dp

   END SUBROUTINE ecd_gamma_context_release

END MODULE qs_ecd_gamma_data
