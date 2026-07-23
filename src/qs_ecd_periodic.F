!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Periodic electronic circular dichroism drivers.
! **************************************************************************************************
MODULE qs_ecd_periodic
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_reserve_blocks,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE eigenvalueproblems,               ONLY: singular_values
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: rskp_directional_taylor_transform,&
                                              rskp_transform
   USE mathconstants,                   ONLY: twopi
   USE orbital_pointers,                 ONLY: coset
   USE qs_ecd_subspace,                 ONLY: ecd_band_first_order_response,&
                                              ecd_band_manifold_type,&
                                              ecd_band_manifolds_validate_window,&
                                              ecd_taylor_retarded_response_pair,&
                                              ecd_taylor_retarded_response
   USE qs_ecd_operators,                ONLY: ecd_periodic_gth_operators_build,&
                                              ecd_periodic_gth_operators_get,&
                                              ecd_periodic_gth_operators_release,&
                                              ecd_periodic_gth_operators_type,&
                                              ecd_periodic_local_operators_build,&
                                              ecd_periodic_local_operators_get,&
                                              ecd_periodic_local_operators_release,&
                                              ecd_periodic_local_operators_type,&
                                              ecd_project_real_operator_to_mo
   USE qs_ecd_gamma_data,               ONLY: ecd_gamma_context_create_implicit,&
                                              ecd_gamma_context_core_images_get,&
                                              ecd_gamma_context_local_potential_images_get,&
                                              ecd_gamma_context_overlap_images_get,&
                                              ecd_gamma_context_overlap_images_status,&
                                              ecd_gamma_context_release,&
                                              ecd_gamma_spin_type,&
                                              ecd_gamma_context_type
   USE qs_ecd_gamma_operators,          ONLY: ecd_gamma_zero_order_operators_build,&
                                              ecd_gamma_zero_order_operators_get,&
                                              ecd_gamma_zero_order_operators_release,&
                                              ecd_gamma_zero_order_operators_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE rt_propagation_velocity_gauge,    ONLY: velocity_gauge_nl_term
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list,&
                                              task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_periodic'
   INTEGER, PARAMETER, PRIVATE          :: ecd_broadening_lorentzian = 3, &
                                           ncart = 3

   TYPE, PRIVATE :: ecd_periodic_transition_selection_type
      INTEGER                                      :: candidate_count = 0, requested_count = 0, selected_count = 0
      INTEGER, DIMENSION(2)                         :: band_range = (/ -1, -1 /)
      REAL(KIND=dp)                                :: cutoff_gap = 0.0_dp
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: iocc, ispin, ivirt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: gap, occupation_difference
   END TYPE ecd_periodic_transition_selection_type

   PUBLIC :: ecd_periodic_ipa_response

CONTAINS

! **************************************************************************************************
!> \brief Validate and dispatch the periodic independent-particle ECD response.
!> \param qs_env Quickstep environment
!> \param ecd_section ECD input section
!> \param ipa_section ECD/IPA input section
!> \param ntransitions Number of transition-resolved records requested
!> \param transition_band_range Band range included in the response
!> \param pseudopotential_commutator Select complete-current or local-current primary output
!> \param ref_point Spatial-expansion reference point
!> \param broadening Spectrum broadening selector
!> \param fwhm Spectrum full width at half maximum
! **************************************************************************************************
   SUBROUTINE ecd_periodic_ipa_response(qs_env, ecd_section, ipa_section, ntransitions, &
                                        transition_band_range, pseudopotential_commutator, &
                                        ref_point, broadening, fwhm)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: ecd_section, ipa_section
      INTEGER, INTENT(IN)                                :: ntransitions
      INTEGER, DIMENSION(2), INTENT(IN)                  :: transition_band_range
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: ref_point
      INTEGER, INTENT(IN)                                :: broadening
      REAL(KIND=dp), INTENT(IN)                          :: fwhm

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: velocity_first_direction, velocity_response_zero, &
                                                            velocity_response_zero_direction, velocity_transition_direction
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: velocity_first_cartesian, velocity_transition_cartesian
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: direct_position_transition, static_magnetic_transition_direction
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: velocity_complete_first_direction, &
                                                            velocity_complete_zero, velocity_complete_zero_direction, &
                                                            velocity_cross_first_direction, velocity_cross_zero, &
                                                            velocity_cross_zero_direction, velocity_gth_transition_direction, &
                                                            velocity_local_first_direction, velocity_local_zero, &
                                                            velocity_local_zero_direction, velocity_nonlocal_first_direction, &
                                                            velocity_nonlocal_zero, velocity_nonlocal_zero_direction
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: velocity_complete_first_cartesian, &
                                                            velocity_cross_first_cartesian, velocity_gth_transition_cartesian, &
                                                            velocity_local_first_cartesian, velocity_nonlocal_first_cartesian
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: static_magnetic_transition_cartesian
      INTEGER                                            :: idir, itransition, npoints, spatial_order
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: selected_iocc, selected_ispin, selected_ivirt
      LOGICAL                                            :: do_kpoints, energy_range_explicit, gth_available, &
                                                            gth_available_reference
      REAL(KIND=dp)                                      :: eps_degeneracy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gap_first_cartesian
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy_range
      REAL(KIND=dp), DIMENSION(ncart)                    :: direction
      REAL(KIND=dp), DIMENSION(ncart, ncart)             :: diamagnetic_kernel, diamagnetic_kernel_direction
      TYPE(cell_type), POINTER                           :: cell
      TYPE(ecd_gamma_context_type)                       :: gamma_context
      TYPE(ecd_gamma_zero_order_operators_type)          :: gamma_operators
      TYPE(ecd_periodic_transition_selection_type)       :: transition_selection
      TYPE(section_vals_type), POINTER                   :: response_section

      NULLIFY (cell, energy_range, response_section)

      CALL get_qs_env(qs_env, cell=cell, do_kpoints=do_kpoints)
      CPASSERT(ASSOCIATED(cell))

      IF (.NOT. ALL(cell%perd == 1)) THEN
         CPABORT("ECD| Periodic IPA currently requires CELL%PERIODIC XYZ")
      END IF
      IF (do_kpoints) THEN
         CPABORT("ECD| Periodic IPA production path requires the implicit-Gamma SCF layout")
      END IF

      CALL section_vals_val_get(ecd_section, "EPS_DEGENERACY", r_val=eps_degeneracy)
      CALL section_vals_val_get(ecd_section, "SPATIAL_ORDER", i_val=spatial_order)

      IF (eps_degeneracy <= 0.0_dp) THEN
         CPABORT("ECD| EPS_DEGENERACY must be positive")
      END IF
      IF (spatial_order /= 1) THEN
         CPABORT("ECD| Periodic ECD response currently requires SPATIAL_ORDER 1")
      END IF
      IF (broadening /= ecd_broadening_lorentzian) THEN
         CPABORT("ECD| Periodic response currently requires BROADENING LORENTZIAN")
      END IF
      IF (fwhm <= 0.0_dp) THEN
         CPABORT("ECD| Periodic response requires positive FWHM")
      END IF
      IF (ntransitions == 0 .OR. ntransitions < -1) THEN
         CPABORT("ECD| Periodic IPA requires NTRANSITIONS -1 or a positive value")
      END IF
      IF (transition_band_range(1) == -1 .OR. transition_band_range(2) == -1) THEN
         IF (ANY(transition_band_range /= -1)) THEN
            CPABORT("ECD| TRANSITION_BAND_RANGE uses -1 -1 or two positive bounds")
         END IF
      ELSE IF (transition_band_range(1) < 1 .OR. &
               transition_band_range(1) > transition_band_range(2)) THEN
         CPABORT("ECD| Invalid periodic TRANSITION_BAND_RANGE")
      END IF

      response_section => section_vals_get_subs_vals(ecd_section, "PRINT%RESPONSE_TENSOR")
      CALL section_vals_val_get(response_section, "ENERGY_RANGE", explicit=energy_range_explicit)
      IF (energy_range_explicit) THEN
         CALL section_vals_val_get(response_section, "ENERGY_RANGE", r_vals=energy_range)
      END IF
      CALL section_vals_val_get(response_section, "NPOINTS", i_val=npoints)
      IF (.NOT. energy_range_explicit) THEN
         CPABORT("ECD| Periodic response requires PRINT%RESPONSE_TENSOR%ENERGY_RANGE")
      END IF
      IF (.NOT. ASSOCIATED(energy_range)) THEN
         CPABORT("ECD| ENERGY_RANGE data are unavailable")
      END IF
      IF (SIZE(energy_range) /= 2) THEN
         CPABORT("ECD| ENERGY_RANGE expects lower and upper photon energies")
      END IF
      IF (energy_range(1) < 0.0_dp .OR. energy_range(1) >= energy_range(2)) THEN
         CPABORT("ECD| ENERGY_RANGE requires 0 <= lower < upper")
      END IF
      IF (energy_range(1) <= 0.0_dp) THEN
         CPABORT("ECD| Periodic ECD response requires ENERGY_RANGE lower > 0")
      END IF
      IF (npoints < 2) THEN
         CPABORT("ECD| Periodic response requires NPOINTS of at least two")
      END IF

      DO idir = 1, ncart
         direction = 0.0_dp
         direction(idir) = 1.0_dp
         CALL ecd_gamma_context_create_implicit(gamma_context, qs_env, direction, ref_point, spatial_order, &
                                                eps_degeneracy)
         CALL ecd_gamma_zero_order_operators_build(gamma_operators, qs_env, gamma_context)
         CALL ecd_periodic_implicit_gamma_overlap_g3(gamma_context)
         CALL ecd_periodic_implicit_gamma_current_g4(qs_env, gamma_context, gamma_operators)
         CALL ecd_periodic_implicit_gamma_ks_g6(qs_env, gamma_context)
         CALL ecd_periodic_implicit_gamma_ks_phase_g7(gamma_context)
         CALL ecd_periodic_implicit_gamma_state_g8(gamma_context)
         CALL ecd_periodic_transition_selection_build(transition_selection, gamma_context, ntransitions, &
                                                      transition_band_range, eps_degeneracy)
         IF (idir == 1) THEN
            ALLOCATE (gap_first_cartesian(transition_selection%selected_count, ncart), &
                      selected_iocc(transition_selection%selected_count), &
                      selected_ispin(transition_selection%selected_count), &
                      selected_ivirt(transition_selection%selected_count), &
                      velocity_complete_first_cartesian(ncart, ncart, ncart, npoints), &
                      velocity_cross_first_cartesian(ncart, ncart, ncart, npoints), &
                      velocity_first_cartesian(ncart, ncart, ncart, npoints), &
                      velocity_gth_transition_cartesian(ncart, 0:1, transition_selection%selected_count, ncart), &
                      velocity_local_first_cartesian(ncart, ncart, ncart, npoints), &
                      velocity_nonlocal_first_cartesian(ncart, ncart, ncart, npoints), &
                      velocity_transition_cartesian(ncart, 0:1, transition_selection%selected_count, ncart), &
                      static_magnetic_transition_cartesian(ncart, transition_selection%selected_count, ncart))
            selected_iocc = transition_selection%iocc
            selected_ispin = transition_selection%ispin
            selected_ivirt = transition_selection%ivirt
         ELSE
            IF (transition_selection%selected_count /= SIZE(selected_iocc)) THEN
               CPABORT('ECD| Cartesian Gamma directions selected incompatible occupied-empty transition counts')
            END IF
            IF (ANY(transition_selection%iocc /= selected_iocc) .OR. &
                ANY(transition_selection%ispin /= selected_ispin) .OR. &
                ANY(transition_selection%ivirt /= selected_ivirt)) THEN
               CPABORT('ECD| Cartesian Gamma directions selected inconsistent occupied-empty transitions')
            END IF
         END IF
         DO itransition = 1, transition_selection%selected_count
            gap_first_cartesian(itransition, idir) = 0.5_dp*( &
               gamma_context%spin(transition_selection%ispin(itransition))% &
                  eigenvalue_first(transition_selection%ivirt(itransition)) + &
               gamma_context%spin(transition_selection%ispin(itransition))% &
                  eigenvalue_first(transition_selection%iocc(itransition)))
         END DO
         IF (idir == 1) THEN
            CALL ecd_periodic_implicit_gamma_velocity_response_g10(qs_env, gamma_context, gamma_operators, &
                                                                    transition_selection, pseudopotential_commutator, &
                                                                    energy_range, npoints, fwhm, &
                                                                    velocity_response_zero_direction, velocity_first_direction, &
                                                                    velocity_complete_zero_direction, &
                                                                    velocity_complete_first_direction, &
                                                                    velocity_local_zero_direction, velocity_local_first_direction, &
                                                                    velocity_cross_zero_direction, velocity_cross_first_direction, &
                                                                    velocity_nonlocal_zero_direction, &
                                                                    velocity_nonlocal_first_direction, &
                                                                    diamagnetic_kernel_direction, velocity_transition_direction, &
                                                                    velocity_gth_transition_direction, &
                                                                    static_magnetic_transition_direction, &
                                                                    gth_available)
            diamagnetic_kernel = diamagnetic_kernel_direction
            gth_available_reference = gth_available
         ELSE
            CALL ecd_periodic_implicit_gamma_velocity_response_g10(qs_env, gamma_context, gamma_operators, &
                                                                    transition_selection, pseudopotential_commutator, &
                                                                    energy_range, npoints, fwhm, &
                                                                    velocity_response_zero_direction, velocity_first_direction, &
                                                                    velocity_complete_zero_direction, &
                                                                    velocity_complete_first_direction, &
                                                                    velocity_local_zero_direction, velocity_local_first_direction, &
                                                                    velocity_cross_zero_direction, velocity_cross_first_direction, &
                                                                    velocity_nonlocal_zero_direction, &
                                                                    velocity_nonlocal_first_direction, &
                                                                    diamagnetic_kernel_direction, velocity_transition_direction, &
                                                                    velocity_gth_transition_direction, &
                                                                    static_magnetic_transition_direction, &
                                                                    gth_available, diamagnetic_kernel)
            IF (gth_available .NEQV. gth_available_reference) THEN
               CPABORT('ECD| Cartesian Gamma directions reported inconsistent GTH-current availability')
            END IF
         END IF
         IF (idir == 1) THEN
            CALL MOVE_ALLOC(velocity_response_zero_direction, velocity_response_zero)
         ELSE
            CALL ecd_periodic_response_zero_check(velocity_response_zero, velocity_response_zero_direction, 'velocity')
            DEALLOCATE (velocity_response_zero_direction)
         END IF
         IF (idir == 1) THEN
            CALL MOVE_ALLOC(velocity_complete_zero_direction, velocity_complete_zero)
            CALL MOVE_ALLOC(velocity_local_zero_direction, velocity_local_zero)
            CALL MOVE_ALLOC(velocity_cross_zero_direction, velocity_cross_zero)
            CALL MOVE_ALLOC(velocity_nonlocal_zero_direction, velocity_nonlocal_zero)
         ELSE
            CALL ecd_periodic_response_zero_check(velocity_complete_zero, velocity_complete_zero_direction, &
                                                  'complete velocity')
            CALL ecd_periodic_response_zero_check(velocity_local_zero, velocity_local_zero_direction, &
                                                  'local velocity')
            CALL ecd_periodic_response_zero_check(velocity_cross_zero, velocity_cross_zero_direction, &
                                                  'cross velocity')
            CALL ecd_periodic_response_zero_check(velocity_nonlocal_zero, velocity_nonlocal_zero_direction, &
                                                  'nonlocal velocity')
            DEALLOCATE (velocity_complete_zero_direction, velocity_local_zero_direction, &
                        velocity_cross_zero_direction, velocity_nonlocal_zero_direction)
         END IF
         velocity_first_cartesian(:, :, idir, :) = velocity_first_direction
         velocity_complete_first_cartesian(:, :, idir, :) = velocity_complete_first_direction
         velocity_local_first_cartesian(:, :, idir, :) = velocity_local_first_direction
         velocity_cross_first_cartesian(:, :, idir, :) = velocity_cross_first_direction
         velocity_nonlocal_first_cartesian(:, :, idir, :) = velocity_nonlocal_first_direction
         velocity_transition_cartesian(:, :, :, idir) = velocity_transition_direction
         velocity_gth_transition_cartesian(:, :, :, idir) = velocity_gth_transition_direction
         static_magnetic_transition_cartesian(:, :, idir) = static_magnetic_transition_direction
         DEALLOCATE (velocity_first_direction, &
                     velocity_complete_first_direction, velocity_cross_first_direction, &
                     velocity_gth_transition_direction, velocity_local_first_direction, &
                     velocity_nonlocal_first_direction, velocity_transition_direction, &
                     static_magnetic_transition_direction)
         IF (idir < ncart) THEN
            CALL ecd_periodic_transition_selection_release(transition_selection)
            CALL ecd_gamma_zero_order_operators_release(gamma_operators)
            CALL ecd_gamma_context_release(gamma_context)
         END IF
      END DO
      CALL ecd_periodic_implicit_gamma_direct_position_transitions(gamma_context, gamma_operators, &
                                                                    transition_selection, direct_position_transition)
      CALL ecd_periodic_print_cartesian_response(ecd_section, energy_range, velocity_response_zero, &
                                                 velocity_first_cartesian, gamma_context, gth_available, &
                                                 pseudopotential_commutator, fwhm, &
                                                 transition_selection)
      CALL ecd_periodic_print_cartesian_transition_matrix_elements(ecd_section, gamma_context, transition_selection, &
                                                                    pseudopotential_commutator, &
                                                                    gap_first_cartesian, &
                                                                    velocity_transition_cartesian, &
                                                                    velocity_gth_transition_cartesian, &
                                                                    direct_position_transition, &
                                                                    static_magnetic_transition_cartesian)
      CALL ecd_periodic_print_cartesian_velocity_decomposition(ecd_section, energy_range, velocity_local_zero, &
                                                                velocity_local_first_cartesian, velocity_cross_zero, &
                                                                velocity_cross_first_cartesian, velocity_nonlocal_zero, &
                                                                velocity_nonlocal_first_cartesian, velocity_complete_zero, &
                                                                velocity_complete_first_cartesian, gamma_context, &
                                                                gth_available, fwhm, transition_selection)
      DEALLOCATE (direct_position_transition, gap_first_cartesian, &
                  selected_iocc, selected_ispin, selected_ivirt, velocity_first_cartesian, &
                  velocity_response_zero, velocity_transition_cartesian, velocity_complete_first_cartesian, &
                  velocity_complete_zero, velocity_cross_first_cartesian, velocity_cross_zero, &
                  velocity_gth_transition_cartesian, &
                  velocity_local_first_cartesian, velocity_local_zero, &
                  velocity_nonlocal_first_cartesian, velocity_nonlocal_zero, &
                  static_magnetic_transition_cartesian)
      CALL ecd_periodic_transition_selection_release(transition_selection)
      CALL ecd_gamma_zero_order_operators_release(gamma_operators)
      CALL ecd_gamma_context_release(gamma_context)

   END SUBROUTINE ecd_periodic_ipa_response

! **************************************************************************************************
!> \brief Build the immutable energy-ordered periodic IPA transition selection.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_transition_selection_build(selection, context, nrequested, transition_band_range, &
                                                       eps_degeneracy)
      TYPE(ecd_periodic_transition_selection_type), INTENT(INOUT) :: selection
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      INTEGER, INTENT(IN)                                :: nrequested
      INTEGER, DIMENSION(2), INTENT(IN)                  :: transition_band_range
      REAL(KIND=dp), INTENT(IN)                          :: eps_degeneracy

      INTEGER                                            :: first_band, icandidate, ispin, ivirt, iocc, last_band, &
                                                            nbase, ncandidates, nmo
      REAL(KIND=dp), PARAMETER                           :: occupation_threshold = 1.0E-8_dp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: candidate_iocc, candidate_ispin, candidate_ivirt, &
                                                            candidate_occ_manifold, candidate_virt_manifold, order
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: candidate_selected
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: candidate_gap, candidate_occupation_difference

      CALL ecd_periodic_transition_selection_release(selection)
      IF (nrequested == 0 .OR. nrequested < -1 .OR. eps_degeneracy <= 0.0_dp) THEN
         CPABORT('ECD| Periodic IPA transition selection received invalid controls')
      END IF

      ncandidates = 0
      DO ispin = 1, context%nspin
         nmo = context%spin(ispin)%nmo
         IF (.NOT. ASSOCIATED(context%spin(ispin)%eigenvalue_reference) .OR. &
             .NOT. ASSOCIATED(context%spin(ispin)%occupations) .OR. &
             .NOT. ASSOCIATED(context%spin(ispin)%manifolds)) THEN
            CPABORT('ECD| Periodic IPA transition selection requires retained Gamma state data')
         END IF
         CALL ecd_periodic_transition_window(context%spin(ispin), transition_band_range, first_band, last_band)
         DO iocc = first_band, last_band
            IF (context%spin(ispin)%occupations(iocc) <= occupation_threshold) CYCLE
            DO ivirt = first_band, last_band
               IF (context%spin(ispin)%occupations(ivirt) > occupation_threshold) CYCLE
               IF (context%spin(ispin)%eigenvalue_reference(ivirt) - &
                   context%spin(ispin)%eigenvalue_reference(iocc) <= SQRT(TINY(1.0_dp))) CYCLE
               IF (context%spin(ispin)%occupations(iocc) - &
                   context%spin(ispin)%occupations(ivirt) <= occupation_threshold) CYCLE
               ncandidates = ncandidates + 1
            END DO
         END DO
      END DO
      IF (ncandidates < 1) THEN
         CPABORT('ECD| Periodic IPA transition selection found no occupied-empty candidates')
      END IF

      ALLOCATE (candidate_gap(ncandidates), candidate_iocc(ncandidates), candidate_ispin(ncandidates), &
                candidate_ivirt(ncandidates), candidate_occ_manifold(ncandidates), &
                candidate_virt_manifold(ncandidates), candidate_occupation_difference(ncandidates), order(ncandidates))
      icandidate = 0
      DO ispin = 1, context%nspin
         CALL ecd_periodic_transition_window(context%spin(ispin), transition_band_range, first_band, last_band)
         DO iocc = first_band, last_band
            IF (context%spin(ispin)%occupations(iocc) <= occupation_threshold) CYCLE
            DO ivirt = first_band, last_band
               IF (context%spin(ispin)%occupations(ivirt) > occupation_threshold) CYCLE
               IF (context%spin(ispin)%eigenvalue_reference(ivirt) - &
                   context%spin(ispin)%eigenvalue_reference(iocc) <= SQRT(TINY(1.0_dp))) CYCLE
               IF (context%spin(ispin)%occupations(iocc) - &
                   context%spin(ispin)%occupations(ivirt) <= occupation_threshold) CYCLE
               icandidate = icandidate + 1
               candidate_gap(icandidate) = context%spin(ispin)%eigenvalue_reference(ivirt) - &
                                            context%spin(ispin)%eigenvalue_reference(iocc)
               candidate_iocc(icandidate) = iocc
               candidate_ispin(icandidate) = ispin
               candidate_ivirt(icandidate) = ivirt
               candidate_occ_manifold(icandidate) = ecd_periodic_band_manifold_index( &
                  context%spin(ispin)%manifolds, iocc)
               candidate_virt_manifold(icandidate) = ecd_periodic_band_manifold_index( &
                  context%spin(ispin)%manifolds, ivirt)
               candidate_occupation_difference(icandidate) = context%spin(ispin)%occupations(iocc) - &
                                                               context%spin(ispin)%occupations(ivirt)
            END DO
         END DO
      END DO
      CPASSERT(icandidate == ncandidates)

      CALL ecd_periodic_transition_order_sort(candidate_gap, candidate_ispin, candidate_iocc, candidate_ivirt, order)
      IF (nrequested == -1) THEN
         nbase = ncandidates
      ELSE
         nbase = MIN(nrequested, ncandidates)
      END IF
      selection%cutoff_gap = candidate_gap(order(nbase))
      ALLOCATE (candidate_selected(ncandidates))
      candidate_selected = .FALSE.
      IF (nrequested == -1) THEN
         candidate_selected = .TRUE.
      ELSE
         DO icandidate = 1, nbase
            candidate_selected(order(icandidate)) = .TRUE.
         END DO
         DO icandidate = 1, ncandidates
            IF (.NOT. candidate_selected(icandidate)) CYCLE
            DO nbase = 1, ncandidates
               IF (candidate_ispin(nbase) /= candidate_ispin(icandidate)) CYCLE
               IF (candidate_occ_manifold(nbase) /= candidate_occ_manifold(icandidate)) CYCLE
               IF (candidate_virt_manifold(nbase) /= candidate_virt_manifold(icandidate)) CYCLE
               candidate_selected(nbase) = .TRUE.
            END DO
         END DO
      END IF
      selection%selected_count = COUNT(candidate_selected)
      IF (nrequested > ncandidates) THEN
         CPWARN('ECD| Requested more periodic IPA transitions than available')
      END IF
      selection%candidate_count = ncandidates
      selection%requested_count = nrequested
      selection%band_range = transition_band_range
      ALLOCATE (selection%gap(selection%selected_count), selection%iocc(selection%selected_count), &
                selection%ispin(selection%selected_count), selection%ivirt(selection%selected_count), &
                selection%occupation_difference(selection%selected_count))
      nbase = 0
      DO icandidate = 1, ncandidates
         IF (.NOT. candidate_selected(order(icandidate))) CYCLE
         nbase = nbase + 1
         selection%gap(nbase) = candidate_gap(order(icandidate))
         selection%iocc(nbase) = candidate_iocc(order(icandidate))
         selection%ispin(nbase) = candidate_ispin(order(icandidate))
         selection%ivirt(nbase) = candidate_ivirt(order(icandidate))
         selection%occupation_difference(nbase) = candidate_occupation_difference(order(icandidate))
      END DO
      CPASSERT(nbase == selection%selected_count)
      DEALLOCATE (candidate_gap, candidate_iocc, candidate_ispin, candidate_ivirt, candidate_occ_manifold, &
                  candidate_virt_manifold, candidate_occupation_difference, candidate_selected, order)

   END SUBROUTINE ecd_periodic_transition_selection_build

! **************************************************************************************************
!> \brief Return the complete degenerate-manifold index containing a band.
! **************************************************************************************************
   FUNCTION ecd_periodic_band_manifold_index(manifolds, iband) RESULT(imanifold)
      TYPE(ecd_band_manifold_type), DIMENSION(:), POINTER :: manifolds
      INTEGER, INTENT(IN)                                :: iband
      INTEGER                                            :: imanifold

      IF (.NOT. ASSOCIATED(manifolds)) THEN
         CPABORT('ECD| Periodic transition selection requires degenerate-manifold metadata')
      END IF
      DO imanifold = 1, SIZE(manifolds)
         IF (iband >= manifolds(imanifold)%first_band .AND. iband <= manifolds(imanifold)%last_band) RETURN
      END DO
      CPABORT('ECD| Periodic transition selection encountered a band outside the manifold partition')

   END FUNCTION ecd_periodic_band_manifold_index

! **************************************************************************************************
!> \brief Resolve a complete-manifold periodic IPA transition window for one spin channel.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_transition_window(spin, transition_band_range, first_band, last_band)
      TYPE(ecd_gamma_spin_type), INTENT(IN)              :: spin
      INTEGER, DIMENSION(2), INTENT(IN)                  :: transition_band_range
      INTEGER, INTENT(OUT)                               :: first_band, last_band

      IF (ALL(transition_band_range == -1)) THEN
         first_band = 1
         last_band = spin%nmo
      ELSE
         first_band = transition_band_range(1)
         last_band = transition_band_range(2)
         CALL ecd_band_manifolds_validate_window(spin%manifolds, transition_band_range)
      END IF
      IF (first_band < 1 .OR. last_band > spin%nmo .OR. first_band > last_band) THEN
         CPABORT('ECD| Periodic IPA transition band range is inconsistent with the Gamma MO space')
      END IF

   END SUBROUTINE ecd_periodic_transition_window

! **************************************************************************************************
!> \brief Deterministically order periodic IPA candidates by gap, spin, occupied band, and virtual band.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_transition_order_sort(gap, ispin, iocc, ivirt, order)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: gap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: iocc, ispin, ivirt
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: order

      INTEGER                                            :: icandidate, insert_at, key

      IF (SIZE(gap) /= SIZE(ispin) .OR. SIZE(gap) /= SIZE(iocc) .OR. SIZE(gap) /= SIZE(ivirt) .OR. &
          SIZE(gap) /= SIZE(order)) THEN
         CPABORT('ECD| Periodic IPA transition-order arrays are inconsistent')
      END IF
      DO icandidate = 1, SIZE(order)
         order(icandidate) = icandidate
      END DO
      DO icandidate = 2, SIZE(order)
         key = order(icandidate)
         insert_at = icandidate - 1
         DO WHILE (insert_at >= 1)
            IF (.NOT. ecd_periodic_transition_precedes(key, order(insert_at), gap, ispin, iocc, ivirt)) EXIT
            order(insert_at + 1) = order(insert_at)
            insert_at = insert_at - 1
         END DO
         order(insert_at + 1) = key
      END DO

   END SUBROUTINE ecd_periodic_transition_order_sort

! **************************************************************************************************
!> \brief Compare two periodic IPA candidate indices using deterministic lexical keys.
! **************************************************************************************************
   LOGICAL FUNCTION ecd_periodic_transition_precedes(left, right, gap, ispin, iocc, ivirt)
      INTEGER, INTENT(IN)                                :: left, right
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: gap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: iocc, ispin, ivirt

      ecd_periodic_transition_precedes = .FALSE.
      IF (gap(left) < gap(right)) THEN
         ecd_periodic_transition_precedes = .TRUE.
      ELSE IF (gap(left) == gap(right)) THEN
         IF (ispin(left) < ispin(right)) THEN
            ecd_periodic_transition_precedes = .TRUE.
         ELSE IF (ispin(left) == ispin(right) .AND. iocc(left) < iocc(right)) THEN
            ecd_periodic_transition_precedes = .TRUE.
         ELSE IF (ispin(left) == ispin(right) .AND. iocc(left) == iocc(right) .AND. &
                  ivirt(left) < ivirt(right)) THEN
            ecd_periodic_transition_precedes = .TRUE.
         END IF
      END IF

   END FUNCTION ecd_periodic_transition_precedes

! **************************************************************************************************
!> \brief Release periodic IPA transition-selection storage.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_transition_selection_release(selection)
      TYPE(ecd_periodic_transition_selection_type), INTENT(INOUT) :: selection

      IF (ALLOCATED(selection%gap)) DEALLOCATE (selection%gap)
      IF (ALLOCATED(selection%iocc)) DEALLOCATE (selection%iocc)
      IF (ALLOCATED(selection%ispin)) DEALLOCATE (selection%ispin)
      IF (ALLOCATED(selection%ivirt)) DEALLOCATE (selection%ivirt)
      IF (ALLOCATED(selection%occupation_difference)) DEALLOCATE (selection%occupation_difference)
      selection%candidate_count = 0
      selection%requested_count = 0
      selection%selected_count = 0
      selection%band_range = (/ -1, -1 /)
      selection%cutoff_gap = 0.0_dp

   END SUBROUTINE ecd_periodic_transition_selection_release

! **************************************************************************************************
!> \brief Validate the ECD-owned implicit-Gamma overlap image snapshot through first order.
!>
!> This provider establishes the image-resolved overlap and its directional
!> Taylor coefficient before KS image assembly and state-response construction.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_overlap_g3(context)
      TYPE(ecd_gamma_context_type), INTENT(INOUT), TARGET :: context

      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-5_dp

      INTEGER                                            :: nimages
      LOGICAL                                            :: available
      REAL(KIND=dp)                                      :: derivative_norm_sq, fd_residual_norm_sq, &
                                                            fold_norm_sq, fold_residual_norm_sq
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_minus, xkp_plus, xkp_step, xkp_zero
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: derivative_imag, derivative_real, fd_imag, &
                                                            fold_residual, fold_zero, minus_imag, minus_real, &
                                                            plus_imag, plus_real
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: overlap_images
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      NULLIFY (cell_to_index, derivative_imag, derivative_real, fd_imag, fold_residual, fold_zero, logger, &
               minus_imag, minus_real, overlap_images, plus_imag, plus_real)
      CALL ecd_gamma_context_overlap_images_status(context, available, nimages)
      IF (.NOT. available) THEN
         CPABORT("ECD| Stage G3 requires an image-resolved overlap snapshot")
      END IF
      CALL ecd_gamma_context_overlap_images_get(context, overlap_images, cell_to_index)
      IF (SIZE(overlap_images, 2) /= nimages) THEN
         CPABORT("ECD| Stage G3 overlap image count is inconsistent")
      END IF
      IF (ASSOCIATED(context%overlap_first_real) .OR. ASSOCIATED(context%overlap_first_imag)) THEN
         CPABORT("ECD| Stage G3 overlap Taylor coefficients are already allocated")
      END IF

      xkp_zero = 0.0_dp
      xkp_step = fd_step*MATMUL(TRANSPOSE(context%cell%hmat), context%direction)/twopi
      xkp_plus = xkp_zero + xkp_step
      xkp_minus = xkp_zero - xkp_step

      CALL ecd_periodic_ao_set_create(fold_zero, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap fold', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(fold_residual, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap fold residual', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(context%overlap_first_real, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap derivative real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(context%overlap_first_imag, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap derivative imag', dbcsr_type_antisymmetric)
      derivative_real => context%overlap_first_real
      derivative_imag => context%overlap_first_imag
      CALL ecd_periodic_ao_set_create(plus_real, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap positive real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(plus_imag, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap positive imag', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(minus_real, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap negative real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(minus_imag, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap negative imag', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(fd_imag, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G3 overlap finite difference imag', dbcsr_type_antisymmetric)

      CALL rskp_transform(fold_zero(1)%matrix, rsmat=overlap_images, ispin=1, xkp=xkp_zero, &
                          cell_to_index=cell_to_index, sab_nl=context%sab_orb)
      CALL dbcsr_set(fold_residual(1)%matrix, 0.0_dp)
      CALL dbcsr_add(fold_residual(1)%matrix, fold_zero(1)%matrix, 0.0_dp, 1.0_dp)
      CALL dbcsr_add(fold_residual(1)%matrix, context%matrix_s(1)%matrix, 1.0_dp, -1.0_dp)
      CALL dbcsr_dot(fold_zero(1)%matrix, fold_zero(1)%matrix, fold_norm_sq)
      CALL dbcsr_dot(fold_residual(1)%matrix, fold_residual(1)%matrix, fold_residual_norm_sq)

      CALL rskp_directional_taylor_transform(derivative_real(1)%matrix, derivative_imag(1)%matrix, &
                                             overlap_images, 1, xkp_zero, cell_to_index, context%sab_orb, &
                                             context%cell%hmat, context%direction, 1)
      CALL rskp_transform(plus_real(1)%matrix, plus_imag(1)%matrix, overlap_images, 1, xkp_plus, &
                          cell_to_index, context%sab_orb)
      CALL rskp_transform(minus_real(1)%matrix, minus_imag(1)%matrix, overlap_images, 1, xkp_minus, &
                          cell_to_index, context%sab_orb)
      CALL dbcsr_set(fd_imag(1)%matrix, 0.0_dp)
      CALL dbcsr_add(fd_imag(1)%matrix, plus_imag(1)%matrix, 0.0_dp, 0.5_dp/fd_step)
      CALL dbcsr_add(fd_imag(1)%matrix, minus_imag(1)%matrix, 1.0_dp, -0.5_dp/fd_step)
      CALL dbcsr_add(fd_imag(1)%matrix, derivative_imag(1)%matrix, 1.0_dp, -1.0_dp)
      CALL dbcsr_dot(derivative_imag(1)%matrix, derivative_imag(1)%matrix, derivative_norm_sq)
      CALL dbcsr_dot(fd_imag(1)%matrix, fd_imag(1)%matrix, fd_residual_norm_sq)

      logger => cp_get_default_logger()
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,I8)') &
            'ECD| Stage G3 overlap translation count', nimages
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G3 overlap fold norm', SQRT(ABS(fold_norm_sq))
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G3 overlap fold residual', SQRT(ABS(fold_residual_norm_sq))
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G3 overlap Taylor derivative norm', SQRT(ABS(derivative_norm_sq))
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G3 overlap finite-difference residual', SQRT(ABS(fd_residual_norm_sq))
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G3 overlap finite-difference step [bohr^-1]', fd_step
      END IF

      NULLIFY (derivative_imag, derivative_real)
      CALL dbcsr_deallocate_matrix_set(fd_imag)
      CALL dbcsr_deallocate_matrix_set(fold_residual)
      CALL dbcsr_deallocate_matrix_set(fold_zero)
      CALL dbcsr_deallocate_matrix_set(minus_imag)
      CALL dbcsr_deallocate_matrix_set(minus_real)
      CALL dbcsr_deallocate_matrix_set(plus_imag)
      CALL dbcsr_deallocate_matrix_set(plus_real)

   END SUBROUTINE ecd_periodic_implicit_gamma_overlap_g3

! **************************************************************************************************
!> \brief Reconstruct and validate the frozen local-potential KS image snapshot.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_ks_g6(qs_env, context)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(INOUT), TARGET :: context

      INTEGER                                            :: ispin, nimages
      REAL(KIND=dp)                                      :: exc_dummy, fold_norm_sq, residual_norm_sq
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_zero
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: core_images, local_images
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: core_fold, ks_fold, ks_residual, local_fold, &
                                                            local_image_spin
      TYPE(dft_control_type), POINTER                     :: dft_control
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(pw_env_type), POINTER                          :: pw_env
      TYPE(pw_pool_type), POINTER                         :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), POINTER                       :: v_hartree_rspace
      TYPE(pw_r3d_rs_type)                                :: local_potential
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER         :: vxc_rho, vxc_tau
      TYPE(qs_ks_env_type), POINTER                       :: ks_env
      TYPE(qs_rho_type), POINTER                          :: rho
      TYPE(section_vals_type), POINTER                    :: input, xc_section
      TYPE(task_list_type), POINTER                       :: image_task_list

      NULLIFY (auxbas_pw_pool, cell_to_index, core_fold, core_images, dft_control, image_task_list, input, ks_env, &
               ks_fold, ks_residual, local_fold, local_images, local_image_spin, logger, pw_env, rho, &
               v_hartree_rspace, vxc_rho, vxc_tau, xc_section)
      CALL ecd_gamma_context_core_images_get(context, core_images, cell_to_index)
      IF (.NOT. ASSOCIATED(context%translation_layout%external_kpoints)) THEN
         CPABORT("ECD| Stage G6 requires the implicit-Gamma translation layout")
      END IF
      nimages = SIZE(core_images, 2)
      IF (nimages < 1) THEN
         CPABORT("ECD| Stage G6 requires at least one KS image")
      END IF
      IF (SIZE(core_images, 1) /= context%nspin .OR. SIZE(context%matrix_ks) /= context%nspin) THEN
         CPABORT("ECD| Stage G6 KS image spin structure is inconsistent")
      END IF
      IF (ASSOCIATED(context%local_potential_images)) THEN
         CPABORT("ECD| Stage G6 local-potential translation snapshot is already allocated")
      END IF

      CALL get_qs_env(qs_env, dft_control=dft_control, input=input, ks_env=ks_env, pw_env=pw_env, rho=rho, &
                      v_hartree_rspace=v_hartree_rspace)
      IF (.NOT. ASSOCIATED(dft_control) .OR. .NOT. ASSOCIATED(input) .OR. .NOT. ASSOCIATED(ks_env) .OR. &
          .NOT. ASSOCIATED(pw_env) .OR. .NOT. ASSOCIATED(rho) .OR. .NOT. ASSOCIATED(v_hartree_rspace)) THEN
         CPABORT("ECD| Stage G6 frozen local-potential data are incomplete")
      END IF
      IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc .OR. dft_control%do_admm .OR. &
          dft_control%qs_control%rigpw .OR. dft_control%qs_control%lrigpw) THEN
         CPABORT("ECD| Stage G6 supports the ordinary GPW local-potential SCF path")
      END IF

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      IF (.NOT. ASSOCIATED(auxbas_pw_pool)) THEN
         CPABORT("ECD| Stage G6 cannot obtain the auxiliary plane-wave pool")
      END IF
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      CALL qs_vxc_create(ks_env, rho, xc_section, vxc_rho, vxc_tau, exc_dummy)
      IF (.NOT. ASSOCIATED(vxc_rho) .OR. SIZE(vxc_rho) /= context%nspin) THEN
         CPABORT("ECD| Stage G6 local XC potential has an inconsistent spin structure")
      END IF
      IF (ASSOCIATED(vxc_tau)) THEN
         CPABORT("ECD| Stage G6 meta-GGA kinetic-density potentials are outside the current scope")
      END IF

      CALL allocate_task_list(image_task_list)
      CALL generate_qs_task_list(ks_env, image_task_list, basis_type="ORB", reorder_rs_grid_ranks=.FALSE., &
                                 skip_load_balance_distributed=dft_control%qs_control%skip_load_balance_distributed, &
                                 sab_orb_external=context%sab_orb, &
                                 ext_kpoints=context%translation_layout%external_kpoints)
      CALL ecd_periodic_task_image_matrix_set_create(local_images, context%matrix_ks, context%nspin, nimages, &
                                                      image_task_list, 'ECD G6 LOCAL POTENTIAL IMAGE')

      DO ispin = 1, context%nspin
         CALL auxbas_pw_pool%create_pw(local_potential)
         CALL pw_copy(v_hartree_rspace, local_potential)
         CALL pw_axpy(vxc_rho(ispin), local_potential, alpha=vxc_rho(ispin)%pw_grid%dvol)
         local_image_spin => local_images(ispin, :)
         CALL integrate_v_rspace(v_rspace=local_potential, hmat_kp=local_image_spin, qs_env=qs_env, &
                                 calculate_forces=.FALSE., task_list_external=image_task_list, &
                                 ext_kpoints=context%translation_layout%external_kpoints)
         CALL auxbas_pw_pool%give_back_pw(local_potential)
         NULLIFY (local_image_spin)
      END DO
      DO ispin = 1, context%nspin
         CALL auxbas_pw_pool%give_back_pw(vxc_rho(ispin))
      END DO
      DEALLOCATE (vxc_rho)
      context%local_potential_images => local_images

      xkp_zero = 0.0_dp
      CALL ecd_periodic_ao_set_create(core_fold, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G6 core fold', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(local_fold, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G6 local-potential fold', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(ks_fold, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G6 KS fold', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(ks_residual, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G6 KS fold residual', dbcsr_type_symmetric)

      logger => cp_get_default_logger()
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,I8)') &
            'ECD| Stage G6 KS translation count', nimages
      END IF
      DO ispin = 1, context%nspin
         CALL rskp_transform(core_fold(ispin)%matrix, rsmat=core_images, ispin=ispin, xkp=xkp_zero, &
                             cell_to_index=cell_to_index, sab_nl=context%sab_orb)
         CALL rskp_transform(local_fold(ispin)%matrix, rsmat=local_images, ispin=ispin, xkp=xkp_zero, &
                             cell_to_index=cell_to_index, sab_nl=context%sab_orb)
         CALL dbcsr_set(ks_fold(ispin)%matrix, 0.0_dp)
         CALL dbcsr_add(ks_fold(ispin)%matrix, core_fold(ispin)%matrix, 0.0_dp, 1.0_dp)
         CALL dbcsr_add(ks_fold(ispin)%matrix, local_fold(ispin)%matrix, 1.0_dp, 1.0_dp)
         CALL dbcsr_set(ks_residual(ispin)%matrix, 0.0_dp)
         CALL dbcsr_add(ks_residual(ispin)%matrix, ks_fold(ispin)%matrix, 0.0_dp, 1.0_dp)
         CALL dbcsr_add(ks_residual(ispin)%matrix, context%matrix_ks(ispin)%matrix, 1.0_dp, -1.0_dp)
         CALL dbcsr_dot(ks_fold(ispin)%matrix, ks_fold(ispin)%matrix, fold_norm_sq)
         CALL dbcsr_dot(ks_residual(ispin)%matrix, ks_residual(ispin)%matrix, residual_norm_sq)
         IF (cp_logger_get_default_io_unit(logger) > 0) THEN
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G6 KS fold norm spin', ispin, SQRT(ABS(fold_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G6 KS fold residual spin', ispin, SQRT(ABS(residual_norm_sq))
         END IF
      END DO
      CALL dbcsr_deallocate_matrix_set(core_fold)
      CALL dbcsr_deallocate_matrix_set(ks_fold)
      CALL dbcsr_deallocate_matrix_set(ks_residual)
      CALL dbcsr_deallocate_matrix_set(local_fold)
      NULLIFY (local_images)
      CALL deallocate_task_list(image_task_list)

   END SUBROUTINE ecd_periodic_implicit_gamma_ks_g6

! **************************************************************************************************
!> \brief Build and validate the first directional KS Taylor coefficient from frozen translations.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_ks_phase_g7(context)
      TYPE(ecd_gamma_context_type), INTENT(INOUT), TARGET :: context

      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-5_dp

      INTEGER                                            :: ispin, nimages
      REAL(KIND=dp)                                      :: fd_residual_norm_sq, first_imag_norm_sq, first_real_norm_sq
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_minus, xkp_plus, xkp_step, xkp_zero
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: core_images, local_potential_images
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: core_first_imag, core_first_real, &
                                                            local_first_imag, local_first_real, fd_imag, &
                                                            minus_imag, minus_real, plus_imag, plus_real
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      NULLIFY (cell_to_index, core_first_imag, core_first_real, core_images, fd_imag, local_first_imag, &
               local_first_real, local_potential_images, logger, minus_imag, minus_real, plus_imag, plus_real)
      CALL ecd_gamma_context_core_images_get(context, core_images, cell_to_index)
      CALL ecd_gamma_context_local_potential_images_get(context, local_potential_images, cell_to_index)
      nimages = context%translation_layout%ntranslations
      IF (SIZE(core_images, 1) /= context%nspin .OR. SIZE(local_potential_images, 1) /= context%nspin .OR. &
          SIZE(core_images, 2) /= nimages .OR. SIZE(local_potential_images, 2) /= nimages) THEN
         CPABORT("ECD| Stage G7 KS translation snapshots have inconsistent dimensions")
      END IF
      IF (ASSOCIATED(context%ks_first_real) .OR. ASSOCIATED(context%ks_first_imag)) THEN
         CPABORT("ECD| Stage G7 KS Taylor coefficients are already allocated")
      END IF

      xkp_zero = 0.0_dp
      xkp_step = fd_step*MATMUL(TRANSPOSE(context%cell%hmat), context%direction)/twopi
      xkp_plus = xkp_zero + xkp_step
      xkp_minus = xkp_zero - xkp_step
      CALL ecd_periodic_ao_set_create(core_first_real, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 core first real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(core_first_imag, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 core first imaginary', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(local_first_real, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 local first real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(local_first_imag, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 local first imaginary', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(context%ks_first_real, context%nspin, context%matrix_ks(1)%matrix, &
                                      context%sab_orb, 'ECD G7 KS first real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(context%ks_first_imag, context%nspin, context%matrix_ks(1)%matrix, &
                                      context%sab_orb, 'ECD G7 KS first imaginary', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(plus_real, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 KS positive real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(plus_imag, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 KS positive imaginary', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(minus_real, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 KS negative real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(minus_imag, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 KS negative imaginary', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(fd_imag, context%nspin, context%matrix_ks(1)%matrix, context%sab_orb, &
                                      'ECD G7 KS finite difference imaginary', dbcsr_type_antisymmetric)

      DO ispin = 1, context%nspin
         CALL rskp_directional_taylor_transform(core_first_real(ispin)%matrix, core_first_imag(ispin)%matrix, &
                                                core_images, ispin, xkp_zero, cell_to_index, context%sab_orb, &
                                                context%cell%hmat, context%direction, 1)
         CALL rskp_directional_taylor_transform(local_first_real(ispin)%matrix, local_first_imag(ispin)%matrix, &
                                                local_potential_images, ispin, xkp_zero, cell_to_index, &
                                                context%sab_orb, context%cell%hmat, context%direction, 1)
         CALL dbcsr_set(context%ks_first_real(ispin)%matrix, 0.0_dp)
         CALL dbcsr_add(context%ks_first_real(ispin)%matrix, core_first_real(ispin)%matrix, 0.0_dp, 1.0_dp)
         CALL dbcsr_add(context%ks_first_real(ispin)%matrix, local_first_real(ispin)%matrix, 1.0_dp, 1.0_dp)
         CALL dbcsr_set(context%ks_first_imag(ispin)%matrix, 0.0_dp)
         CALL dbcsr_add(context%ks_first_imag(ispin)%matrix, core_first_imag(ispin)%matrix, 0.0_dp, 1.0_dp)
         CALL dbcsr_add(context%ks_first_imag(ispin)%matrix, local_first_imag(ispin)%matrix, 1.0_dp, 1.0_dp)

         CALL rskp_transform(plus_real(ispin)%matrix, plus_imag(ispin)%matrix, core_images, ispin, xkp_plus, &
                             cell_to_index, context%sab_orb)
         CALL rskp_transform(plus_real(ispin)%matrix, plus_imag(ispin)%matrix, local_potential_images, ispin, &
                             xkp_plus, cell_to_index, context%sab_orb)
         CALL rskp_transform(minus_real(ispin)%matrix, minus_imag(ispin)%matrix, core_images, ispin, xkp_minus, &
                             cell_to_index, context%sab_orb)
         CALL rskp_transform(minus_real(ispin)%matrix, minus_imag(ispin)%matrix, local_potential_images, ispin, &
                             xkp_minus, cell_to_index, context%sab_orb)
         CALL dbcsr_set(fd_imag(ispin)%matrix, 0.0_dp)
         CALL dbcsr_add(fd_imag(ispin)%matrix, plus_imag(ispin)%matrix, 0.0_dp, 0.5_dp/fd_step)
         CALL dbcsr_add(fd_imag(ispin)%matrix, minus_imag(ispin)%matrix, 1.0_dp, -0.5_dp/fd_step)
         CALL dbcsr_add(fd_imag(ispin)%matrix, context%ks_first_imag(ispin)%matrix, 1.0_dp, -1.0_dp)
         CALL dbcsr_dot(context%ks_first_real(ispin)%matrix, context%ks_first_real(ispin)%matrix, &
                        first_real_norm_sq)
         CALL dbcsr_dot(context%ks_first_imag(ispin)%matrix, context%ks_first_imag(ispin)%matrix, &
                        first_imag_norm_sq)
         CALL dbcsr_dot(fd_imag(ispin)%matrix, fd_imag(ispin)%matrix, fd_residual_norm_sq)
         logger => cp_get_default_logger()
         IF (cp_logger_get_default_io_unit(logger) > 0) THEN
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G7 KS first real norm spin', ispin, SQRT(ABS(first_real_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G7 KS first imaginary norm spin', ispin, SQRT(ABS(first_imag_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G7 KS first finite-difference residual spin', ispin, SQRT(ABS(fd_residual_norm_sq))
         END IF
      END DO
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G7 KS finite-difference step [bohr^-1]', fd_step
      END IF
      CALL dbcsr_deallocate_matrix_set(core_first_imag)
      CALL dbcsr_deallocate_matrix_set(core_first_real)
      CALL dbcsr_deallocate_matrix_set(fd_imag)
      CALL dbcsr_deallocate_matrix_set(local_first_imag)
      CALL dbcsr_deallocate_matrix_set(local_first_real)
      CALL dbcsr_deallocate_matrix_set(minus_imag)
      CALL dbcsr_deallocate_matrix_set(minus_real)
      CALL dbcsr_deallocate_matrix_set(plus_imag)
      CALL dbcsr_deallocate_matrix_set(plus_real)

   END SUBROUTINE ecd_periodic_implicit_gamma_ks_phase_g7

! **************************************************************************************************
!> \brief Solve the first Gamma generalized-eigenproblem response from retained phase coefficients.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_state_g8(context)
      TYPE(ecd_gamma_context_type), INTENT(INOUT), TARGET :: context

      INTEGER                                            :: first_band, iband, imanifold, ispin, last_band, nao, nimages, nmo
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-5_dp
      REAL(KIND=dp)                                      :: direct_energy_residual_max, direct_energy_residual_norm, &
                                                            eigenvector_derivative_norm, equation_residual_norm, &
                                                            manifold_residual_norm, normalization_residual_norm, &
                                                            rotation_non_degenerate_residual, rotation_unitarity_residual, &
                                                            state_derivative_residual_max, state_derivative_residual_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalue_derivative, eigenvalue_reference, &
                                                            minus_eigenvalues, plus_eigenvalues
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_plus, xkp_step, xkp_zero
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: rotation_gram
      TYPE(cp_cfm_type)                                  :: hamiltonian_derivative_mo, hamiltonian_minus, hamiltonian_plus, &
                                                            overlap_derivative_mo, overlap_minus, overlap_plus, &
                                                            minus_eigenvectors, minus_work, plus_eigenvectors, plus_work
      TYPE(cp_fm_type)                                   :: hamiltonian_derivative_raw, hamiltonian_minus_imag, &
                                                            hamiltonian_minus_real, hamiltonian_plus_imag, hamiltonian_plus_real, &
                                                            overlap_derivative_raw, overlap_minus_imag, overlap_minus_real, &
                                                            overlap_plus_imag, overlap_plus_real
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: core_images, local_images, overlap_images
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: hamiltonian_minus_imag_ao, hamiltonian_minus_real_ao, &
                                                            hamiltonian_plus_imag_ao, hamiltonian_plus_real_ao, &
                                                            overlap_minus_imag_ao, overlap_minus_real_ao, &
                                                            overlap_plus_imag_ao, overlap_plus_real_ao
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      NULLIFY (cell_to_index, core_images, hamiltonian_minus_imag_ao, hamiltonian_minus_real_ao, &
               hamiltonian_plus_imag_ao, hamiltonian_plus_real_ao, local_images, logger, overlap_images, &
               overlap_minus_imag_ao, overlap_minus_real_ao, overlap_plus_imag_ao, overlap_plus_real_ao)
      IF (.NOT. ASSOCIATED(context%ks_first_imag) .OR. .NOT. ASSOCIATED(context%overlap_first_imag)) THEN
         CPABORT("ECD| Stage G8 requires retained first-order KS and overlap coefficients")
      END IF
      IF (SIZE(context%ks_first_imag) /= context%nspin .OR. SIZE(context%overlap_first_imag) /= 1) THEN
         CPABORT("ECD| Stage G8 first-order coefficient spin dimensions are inconsistent")
      END IF
      CALL ecd_gamma_context_core_images_get(context, core_images, cell_to_index)
      CALL ecd_gamma_context_local_potential_images_get(context, local_images, cell_to_index)
      CALL ecd_gamma_context_overlap_images_get(context, overlap_images, cell_to_index)
      nimages = context%translation_layout%ntranslations
      IF (nimages < 1 .OR. SIZE(core_images, 2) /= nimages .OR. SIZE(local_images, 2) /= nimages .OR. &
          SIZE(overlap_images, 2) /= nimages) THEN
         CPABORT("ECD| Stage G8 translation snapshot dimensions are inconsistent")
      END IF
      IF (SIZE(core_images, 1) /= context%nspin .OR. SIZE(local_images, 1) /= context%nspin .OR. &
          SIZE(overlap_images, 1) /= 1) THEN
         CPABORT("ECD| Stage G8 translation snapshot spin dimensions are inconsistent")
      END IF

      xkp_zero = 0.0_dp
      xkp_step = fd_step*MATMUL(TRANSPOSE(context%cell%hmat), context%direction)/twopi
      xkp_plus = xkp_zero + xkp_step
      CALL ecd_periodic_ao_set_create(hamiltonian_plus_real_ao, context%nspin, context%matrix_ks(1)%matrix, &
                                      context%sab_orb, 'ECD G8 positive-q KS real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(hamiltonian_plus_imag_ao, context%nspin, context%matrix_ks(1)%matrix, &
                                      context%sab_orb, 'ECD G8 positive-q KS imaginary', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(hamiltonian_minus_real_ao, context%nspin, context%matrix_ks(1)%matrix, &
                                      context%sab_orb, 'ECD G8 negative-q KS real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(hamiltonian_minus_imag_ao, context%nspin, context%matrix_ks(1)%matrix, &
                                      context%sab_orb, 'ECD G8 negative-q KS imaginary', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(overlap_plus_real_ao, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G8 positive-q overlap real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(overlap_plus_imag_ao, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G8 positive-q overlap imaginary', dbcsr_type_antisymmetric)
      CALL ecd_periodic_ao_set_create(overlap_minus_real_ao, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G8 negative-q overlap real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(overlap_minus_imag_ao, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G8 negative-q overlap imaginary', dbcsr_type_antisymmetric)
      CALL rskp_transform(overlap_plus_real_ao(1)%matrix, overlap_plus_imag_ao(1)%matrix, overlap_images, 1, &
                          xkp_plus, cell_to_index, context%sab_orb)
      CALL rskp_transform(overlap_minus_real_ao(1)%matrix, overlap_minus_imag_ao(1)%matrix, overlap_images, 1, &
                          -xkp_step, cell_to_index, context%sab_orb)
      logger => cp_get_default_logger()
      DO ispin = 1, context%nspin
         nmo = context%spin(ispin)%nmo
         nao = context%spin(ispin)%nao
         IF (nao < 1 .OR. nmo < 1 .OR. .NOT. ASSOCIATED(context%spin(ispin)%mo_coeff) .OR. &
             .NOT. ASSOCIATED(context%spin(ispin)%eigenvalues) .OR. &
             .NOT. ASSOCIATED(context%spin(ispin)%manifolds)) THEN
            CPABORT("ECD| Stage G8 Gamma MO data are incomplete")
         END IF
         IF (ASSOCIATED(context%spin(ispin)%eigenvalue_first) .OR. &
             ASSOCIATED(context%spin(ispin)%eigenvalue_reference) .OR. &
             ALLOCATED(context%spin(ispin)%band_rotation) .OR. &
             ALLOCATED(context%spin(ispin)%coefficient_first)) THEN
            CPABORT("ECD| Stage G8 state-response data are already allocated")
         END IF

         CALL ecd_project_real_operator_to_mo(context%ks_first_imag(ispin), context%spin(ispin)%mo_coeff, &
                                              context%spin(ispin)%nao, nmo, hamiltonian_derivative_raw)
         CALL ecd_project_real_operator_to_mo(context%overlap_first_imag(1), context%spin(ispin)%mo_coeff, &
                                              context%spin(ispin)%nao, nmo, overlap_derivative_raw)
         CALL cp_cfm_create(hamiltonian_derivative_mo, hamiltonian_derivative_raw%matrix_struct, &
                            name='ECD G8 KS derivative in MO basis', set_zero=.TRUE.)
         CALL cp_cfm_create(overlap_derivative_mo, overlap_derivative_raw%matrix_struct, &
                            name='ECD G8 overlap derivative in MO basis', set_zero=.TRUE.)
         CALL cp_fm_to_cfm(msourcei=hamiltonian_derivative_raw, mtarget=hamiltonian_derivative_mo)
         CALL cp_fm_to_cfm(msourcei=overlap_derivative_raw, mtarget=overlap_derivative_mo)
         CALL cp_fm_release(hamiltonian_derivative_raw)
         CALL cp_fm_release(overlap_derivative_raw)

         CALL ecd_band_first_order_response(eigenvalue_reference, eigenvalue_derivative, &
                                            equation_residual_norm, normalization_residual_norm, &
                                            manifold_residual_norm, eigenvector_derivative_norm, &
                                            hamiltonian_derivative_mo, overlap_derivative_mo, &
                                            context%spin(ispin)%eigenvalues, context%spin(ispin)%manifolds, &
                                            context%spin(ispin)%band_rotation, context%spin(ispin)%coefficient_first)
         ALLOCATE (context%spin(ispin)%eigenvalue_first(nmo))
         context%spin(ispin)%eigenvalue_first = eigenvalue_derivative
         ALLOCATE (context%spin(ispin)%eigenvalue_reference(nmo))
         context%spin(ispin)%eigenvalue_reference = eigenvalue_reference
         ALLOCATE (rotation_gram(nmo, nmo))
         rotation_gram = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                context%spin(ispin)%band_rotation)
         DO iband = 1, nmo
            rotation_gram(iband, iband) = rotation_gram(iband, iband) - CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         END DO
         rotation_unitarity_residual = SQRT(SUM(ABS(rotation_gram)**2))
         rotation_non_degenerate_residual = 0.0_dp
         DO imanifold = 1, SIZE(context%spin(ispin)%manifolds)
            first_band = context%spin(ispin)%manifolds(imanifold)%first_band
            last_band = context%spin(ispin)%manifolds(imanifold)%last_band
            IF (first_band == last_band) THEN
               rotation_non_degenerate_residual = MAX(rotation_non_degenerate_residual, &
                                                      ABS(context%spin(ispin)%band_rotation(first_band, first_band) - &
                                                          CMPLX(1.0_dp, 0.0_dp, KIND=dp)))
            END IF
         END DO
         DEALLOCATE (rotation_gram)

         CALL rskp_transform(hamiltonian_plus_real_ao(ispin)%matrix, hamiltonian_plus_imag_ao(ispin)%matrix, &
                             core_images, ispin, xkp_plus, cell_to_index, context%sab_orb)
         CALL rskp_transform(hamiltonian_plus_real_ao(ispin)%matrix, hamiltonian_plus_imag_ao(ispin)%matrix, &
                             local_images, ispin, xkp_plus, cell_to_index, context%sab_orb)
         CALL rskp_transform(hamiltonian_minus_real_ao(ispin)%matrix, hamiltonian_minus_imag_ao(ispin)%matrix, &
                             core_images, ispin, -xkp_step, cell_to_index, context%sab_orb)
         CALL rskp_transform(hamiltonian_minus_real_ao(ispin)%matrix, hamiltonian_minus_imag_ao(ispin)%matrix, &
                             local_images, ispin, -xkp_step, cell_to_index, context%sab_orb)
         CALL cp_fm_create(hamiltonian_plus_real, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                           set_zero=.TRUE.)
         CALL cp_fm_create(hamiltonian_plus_imag, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                           set_zero=.TRUE.)
         CALL cp_fm_create(overlap_plus_real, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                           set_zero=.TRUE.)
         CALL cp_fm_create(overlap_plus_imag, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                           set_zero=.TRUE.)
         CALL copy_dbcsr_to_fm(hamiltonian_plus_real_ao(ispin)%matrix, hamiltonian_plus_real)
         CALL copy_dbcsr_to_fm(hamiltonian_plus_imag_ao(ispin)%matrix, hamiltonian_plus_imag)
         CALL copy_dbcsr_to_fm(overlap_plus_real_ao(1)%matrix, overlap_plus_real)
         CALL copy_dbcsr_to_fm(overlap_plus_imag_ao(1)%matrix, overlap_plus_imag)
         CALL cp_cfm_create(hamiltonian_plus, hamiltonian_plus_real%matrix_struct, &
                            name='ECD G8 positive-q KS', set_zero=.TRUE.)
         CALL cp_cfm_create(overlap_plus, overlap_plus_real%matrix_struct, &
                            name='ECD G8 positive-q overlap', set_zero=.TRUE.)
         CALL cp_fm_to_cfm(hamiltonian_plus_real, hamiltonian_plus_imag, hamiltonian_plus)
         CALL cp_fm_to_cfm(overlap_plus_real, overlap_plus_imag, overlap_plus)
         CALL cp_fm_release(hamiltonian_plus_imag)
         CALL cp_fm_release(hamiltonian_plus_real)
         CALL cp_fm_release(overlap_plus_imag)
         CALL cp_fm_release(overlap_plus_real)
         CALL cp_fm_create(hamiltonian_minus_real, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                           set_zero=.TRUE.)
         CALL cp_fm_create(hamiltonian_minus_imag, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                           set_zero=.TRUE.)
         CALL cp_fm_create(overlap_minus_real, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                           set_zero=.TRUE.)
         CALL cp_fm_create(overlap_minus_imag, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                           set_zero=.TRUE.)
         CALL copy_dbcsr_to_fm(hamiltonian_minus_real_ao(ispin)%matrix, hamiltonian_minus_real)
         CALL copy_dbcsr_to_fm(hamiltonian_minus_imag_ao(ispin)%matrix, hamiltonian_minus_imag)
         CALL copy_dbcsr_to_fm(overlap_minus_real_ao(1)%matrix, overlap_minus_real)
         CALL copy_dbcsr_to_fm(overlap_minus_imag_ao(1)%matrix, overlap_minus_imag)
         CALL cp_cfm_create(hamiltonian_minus, hamiltonian_minus_real%matrix_struct, &
                            name='ECD G8 negative-q KS', set_zero=.TRUE.)
         CALL cp_cfm_create(overlap_minus, overlap_minus_real%matrix_struct, &
                            name='ECD G8 negative-q overlap', set_zero=.TRUE.)
         CALL cp_fm_to_cfm(hamiltonian_minus_real, hamiltonian_minus_imag, hamiltonian_minus)
         CALL cp_fm_to_cfm(overlap_minus_real, overlap_minus_imag, overlap_minus)
         CALL cp_fm_release(hamiltonian_minus_imag)
         CALL cp_fm_release(hamiltonian_minus_real)
         CALL cp_fm_release(overlap_minus_imag)
         CALL cp_fm_release(overlap_minus_real)
         CALL cp_cfm_create(plus_eigenvectors, context%spin(ispin)%mo_coeff%matrix_struct, &
                            name='ECD G8 positive-q eigenvectors', set_zero=.TRUE.)
         CALL cp_cfm_create(plus_work, hamiltonian_plus%matrix_struct, &
                            name='ECD G8 positive-q diagonalization work', set_zero=.TRUE.)
         ALLOCATE (plus_eigenvalues(nmo))
         CALL cp_cfm_geeig(hamiltonian_plus, overlap_plus, plus_eigenvectors, plus_eigenvalues, plus_work)
         CALL cp_cfm_create(minus_eigenvectors, context%spin(ispin)%mo_coeff%matrix_struct, &
                            name='ECD G8 negative-q eigenvectors', set_zero=.TRUE.)
         CALL cp_cfm_create(minus_work, hamiltonian_minus%matrix_struct, &
                            name='ECD G8 negative-q diagonalization work', set_zero=.TRUE.)
         ALLOCATE (minus_eigenvalues(nmo))
         CALL cp_cfm_geeig(hamiltonian_minus, overlap_minus, minus_eigenvectors, minus_eigenvalues, minus_work)
         direct_energy_residual_norm = SQRT(SUM((plus_eigenvalues - &
                                                 (eigenvalue_reference + fd_step*eigenvalue_derivative))**2))
         direct_energy_residual_max = MAXVAL(ABS(plus_eigenvalues - &
                                                 (eigenvalue_reference + fd_step*eigenvalue_derivative)))
         CALL ecd_periodic_state_derivative_validate(context, ispin, plus_eigenvectors, minus_eigenvectors, fd_step, &
                                                     state_derivative_residual_norm, state_derivative_residual_max)
         IF (cp_logger_get_default_io_unit(logger) > 0) THEN
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 band-energy derivative norm spin', ispin, SQRT(SUM(eigenvalue_derivative**2))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 manifold residual spin', ispin, manifold_residual_norm
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 eigenproblem residual spin', ispin, equation_residual_norm
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 normalization residual spin', ispin, normalization_residual_norm
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 eigenvector derivative norm spin', ispin, eigenvector_derivative_norm
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 rotation unitarity residual spin', ispin, rotation_unitarity_residual
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 non-degenerate rotation residual spin', ispin, rotation_non_degenerate_residual
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 positive-q energy residual spin', ispin, direct_energy_residual_norm
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 positive-q maximum residual spin', ispin, direct_energy_residual_max
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 central state-derivative residual spin', ispin, state_derivative_residual_norm
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 central state-derivative maximum residual spin', ispin, state_derivative_residual_max
         END IF
         DEALLOCATE (eigenvalue_derivative, eigenvalue_reference, minus_eigenvalues, plus_eigenvalues)
         CALL cp_cfm_release(hamiltonian_derivative_mo)
         CALL cp_cfm_release(overlap_derivative_mo)
         CALL cp_cfm_release(hamiltonian_minus)
         CALL cp_cfm_release(hamiltonian_plus)
         CALL cp_cfm_release(overlap_minus)
         CALL cp_cfm_release(overlap_plus)
         CALL cp_cfm_release(minus_eigenvectors)
         CALL cp_cfm_release(minus_work)
         CALL cp_cfm_release(plus_eigenvectors)
         CALL cp_cfm_release(plus_work)
      END DO
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G8 central-q step [bohr^-1]', fd_step
      END IF
      CALL dbcsr_deallocate_matrix_set(hamiltonian_minus_imag_ao)
      CALL dbcsr_deallocate_matrix_set(hamiltonian_minus_real_ao)
      CALL dbcsr_deallocate_matrix_set(hamiltonian_plus_imag_ao)
      CALL dbcsr_deallocate_matrix_set(hamiltonian_plus_real_ao)
      CALL dbcsr_deallocate_matrix_set(overlap_minus_imag_ao)
      CALL dbcsr_deallocate_matrix_set(overlap_minus_real_ao)
      CALL dbcsr_deallocate_matrix_set(overlap_plus_imag_ao)
      CALL dbcsr_deallocate_matrix_set(overlap_plus_real_ao)

   END SUBROUTINE ecd_periodic_implicit_gamma_state_g8

! **************************************************************************************************
!> \brief Validate retained first-order state coefficients by a centered q difference.
!>
!> Each finite-q generalized eigensolver result is aligned to the reference
!> eigenspace by a manifold-local polar rotation before the centered difference
!> is compared with the analytic G8 coefficient derivative.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_state_derivative_validate(context, ispin, plus_eigenvectors, minus_eigenvectors, fd_step, &
                                                     residual_norm, residual_max)
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET  :: context
      INTEGER, INTENT(IN)                               :: ispin
      TYPE(cp_cfm_type), INTENT(IN)                     :: minus_eigenvectors, plus_eigenvectors
      REAL(KIND=dp), INTENT(IN)                         :: fd_step
      REAL(KIND=dp), INTENT(OUT)                        :: residual_norm, residual_max

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)    :: aligned_overlap_minus, aligned_overlap_plus, &
                                                           finite_difference_derivative, reference_coefficients
      INTEGER                                            :: nao, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap_zero, zero_coefficients
      TYPE(cp_fm_type)                                  :: overlap_zero_fm

      IF (ispin < 1 .OR. ispin > context%nspin .OR. fd_step <= 0.0_dp) THEN
         CPABORT('ECD| Stage G8 state-derivative validation received invalid controls')
      END IF
      nao = context%spin(ispin)%nao
      nmo = context%spin(ispin)%nmo
      IF (nao < 1 .OR. nmo < 1 .OR. .NOT. ASSOCIATED(context%spin(ispin)%mo_coeff) .OR. &
          .NOT. ASSOCIATED(context%spin(ispin)%manifolds) .OR. &
          .NOT. ALLOCATED(context%spin(ispin)%band_rotation) .OR. &
          .NOT. ALLOCATED(context%spin(ispin)%coefficient_first)) THEN
         CPABORT('ECD| Stage G8 state-derivative validation requires retained state data')
      END IF
      IF (SIZE(context%spin(ispin)%band_rotation, 1) /= nmo .OR. &
          SIZE(context%spin(ispin)%band_rotation, 2) /= nmo .OR. &
          SIZE(context%spin(ispin)%coefficient_first, 1) /= nmo .OR. &
          SIZE(context%spin(ispin)%coefficient_first, 2) /= nmo) THEN
         CPABORT('ECD| Stage G8 state-derivative validation received inconsistent state dimensions')
      END IF
      IF (plus_eigenvectors%matrix_struct%nrow_global /= nao .OR. &
          plus_eigenvectors%matrix_struct%ncol_global /= nmo .OR. &
          minus_eigenvectors%matrix_struct%nrow_global /= nao .OR. &
          minus_eigenvectors%matrix_struct%ncol_global /= nmo) THEN
         CPABORT('ECD| Stage G8 central-q eigenvector dimensions are inconsistent')
      END IF

      ALLOCATE (aligned_overlap_minus(nmo, nmo), aligned_overlap_plus(nmo, nmo), &
                finite_difference_derivative(nmo, nmo), reference_coefficients(nao, nmo), &
                overlap_zero(nao, nao), zero_coefficients(nao, nmo))
      CALL cp_fm_get_submatrix(context%spin(ispin)%mo_coeff, zero_coefficients)
      reference_coefficients = MATMUL(CMPLX(zero_coefficients, 0.0_dp, KIND=dp), &
                                      context%spin(ispin)%band_rotation)
      CALL cp_fm_create(overlap_zero_fm, context%spin(ispin)%mo_coeff%matrix_struct, nrow=nao, ncol=nao, &
                        set_zero=.TRUE.)
      CALL copy_dbcsr_to_fm(context%matrix_s(1)%matrix, overlap_zero_fm)
      CALL cp_fm_get_submatrix(overlap_zero_fm, overlap_zero)
      CALL cp_fm_release(overlap_zero_fm)

      CALL ecd_periodic_state_eigenspace_align(context, ispin, plus_eigenvectors, reference_coefficients, &
                                               overlap_zero, aligned_overlap_plus)
      CALL ecd_periodic_state_eigenspace_align(context, ispin, minus_eigenvectors, reference_coefficients, &
                                               overlap_zero, aligned_overlap_minus)
      finite_difference_derivative = (aligned_overlap_plus - aligned_overlap_minus)/(2.0_dp*fd_step)
      residual_norm = SQRT(SUM(ABS(finite_difference_derivative - context%spin(ispin)%coefficient_first)**2))
      residual_max = MAXVAL(ABS(finite_difference_derivative - context%spin(ispin)%coefficient_first))
      DEALLOCATE (aligned_overlap_minus, aligned_overlap_plus, finite_difference_derivative, reference_coefficients, &
                  overlap_zero, zero_coefficients)

   END SUBROUTINE ecd_periodic_state_derivative_validate

! **************************************************************************************************
!> \brief Align one finite-q eigenspace to the retained Gamma reference eigenspace.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_state_eigenspace_align(context, ispin, eigenvectors, reference_coefficients, &
                                                  overlap_zero, aligned_overlap)
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET  :: context
      INTEGER, INTENT(IN)                               :: ispin
      TYPE(cp_cfm_type), INTENT(IN)                     :: eigenvectors
      COMPLEX(KIND=dp), INTENT(IN)                      :: reference_coefficients(:, :)
      REAL(KIND=dp), INTENT(IN)                         :: overlap_zero(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                     :: aligned_overlap(:, :)

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)    :: eigenvector_coefficients, left_singular, &
                                                           manifold_overlap, right_singular_adjoint, rotation
      INTEGER                                            :: first_band, imanifold, last_band, nao, nmanifold, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: singular_value

      nao = context%spin(ispin)%nao
      nmo = context%spin(ispin)%nmo
      IF (SIZE(reference_coefficients, 1) /= nao .OR. SIZE(reference_coefficients, 2) /= nmo .OR. &
          SIZE(overlap_zero, 1) /= nao .OR. SIZE(overlap_zero, 2) /= nao .OR. &
          SIZE(aligned_overlap, 1) /= nmo .OR. SIZE(aligned_overlap, 2) /= nmo) THEN
         CPABORT('ECD| Stage G8 eigenspace alignment received inconsistent matrix dimensions')
      END IF
      IF (eigenvectors%matrix_struct%nrow_global /= nao .OR. eigenvectors%matrix_struct%ncol_global /= nmo) THEN
         CPABORT('ECD| Stage G8 eigenspace alignment received inconsistent eigenvector dimensions')
      END IF

      ALLOCATE (eigenvector_coefficients(nao, nmo))
      CALL cp_cfm_get_submatrix(eigenvectors, eigenvector_coefficients)
      aligned_overlap = MATMUL(CONJG(TRANSPOSE(reference_coefficients)), &
                               MATMUL(CMPLX(overlap_zero, 0.0_dp, KIND=dp), eigenvector_coefficients))
      DO imanifold = 1, SIZE(context%spin(ispin)%manifolds)
         first_band = context%spin(ispin)%manifolds(imanifold)%first_band
         last_band = context%spin(ispin)%manifolds(imanifold)%last_band
         nmanifold = last_band - first_band + 1
         IF (first_band < 1 .OR. last_band > nmo .OR. nmanifold < 1) THEN
            CPABORT('ECD| Stage G8 eigenspace alignment encountered invalid manifold bounds')
         END IF
         ALLOCATE (left_singular(nmanifold, nmanifold), manifold_overlap(nmanifold, nmanifold), &
                   right_singular_adjoint(nmanifold, nmanifold), rotation(nmanifold, nmanifold), &
                   singular_value(nmanifold))
         manifold_overlap = aligned_overlap(first_band:last_band, first_band:last_band)
         CALL singular_values(manifold_overlap, singular_value, nmanifold, nmanifold, left_singular, &
                              right_singular_adjoint)
         rotation = MATMUL(CONJG(TRANSPOSE(right_singular_adjoint)), CONJG(TRANSPOSE(left_singular)))
         aligned_overlap(:, first_band:last_band) = MATMUL(aligned_overlap(:, first_band:last_band), rotation)
         DEALLOCATE (left_singular, manifold_overlap, right_singular_adjoint, rotation, singular_value)
      END DO
      DEALLOCATE (eigenvector_coefficients)

   END SUBROUTINE ecd_periodic_state_eigenspace_align

! **************************************************************************************************
!> \brief Allocate external-translation matrices from the exact real-space task-list block contract.
!>
!> The distributed gather path addresses atom_pair_send while the replicated-grid
!> path addresses atom_pair_recv.  Reserving their union gives every target
!> translation matrix the blocks that either path can write.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_task_image_matrix_set_create(matrices, templates, nspin, nimages, task_list, name)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrices
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: templates
      INTEGER, INTENT(IN)                                :: nspin, nimages
      TYPE(task_list_type), POINTER                       :: task_list
      CHARACTER(LEN=*), INTENT(IN)                       :: name

      INTEGER                                            :: iimage, imatrix, ipair, npairs, nreserve
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: reserve_cols, reserve_rows

      IF (.NOT. ASSOCIATED(task_list)) THEN
         CPABORT("ECD| External translation matrix allocation requires a task list")
      END IF
      IF (nspin < 1 .OR. nimages < 1 .OR. SIZE(templates) /= nspin) THEN
         CPABORT("ECD| External translation matrix allocation has inconsistent dimensions")
      END IF
      npairs = 0
      IF (ASSOCIATED(task_list%atom_pair_recv)) npairs = npairs + SIZE(task_list%atom_pair_recv)
      IF (ASSOCIATED(task_list%atom_pair_send)) npairs = npairs + SIZE(task_list%atom_pair_send)
      ! A compact AO support can leave every external translation block empty.
      ! The corresponding local-potential images are then exact zero matrices.
      IF (npairs > 0) ALLOCATE (reserve_rows(npairs), reserve_cols(npairs))
      CALL dbcsr_allocate_matrix_set(matrices, nspin, nimages)
      DO iimage = 1, nimages
         nreserve = 0
         IF (ASSOCIATED(task_list%atom_pair_recv)) THEN
            DO ipair = 1, SIZE(task_list%atom_pair_recv)
               IF (task_list%atom_pair_recv(ipair)%image < 1 .OR. &
                   task_list%atom_pair_recv(ipair)%image > nimages) THEN
                  CPABORT("ECD| External translation receive pair has an invalid translation index")
               END IF
               IF (task_list%atom_pair_recv(ipair)%image == iimage) THEN
                  nreserve = nreserve + 1
                  reserve_rows(nreserve) = task_list%atom_pair_recv(ipair)%row
                  reserve_cols(nreserve) = task_list%atom_pair_recv(ipair)%col
               END IF
            END DO
         END IF
         IF (ASSOCIATED(task_list%atom_pair_send)) THEN
            DO ipair = 1, SIZE(task_list%atom_pair_send)
               IF (task_list%atom_pair_send(ipair)%image < 1 .OR. &
                   task_list%atom_pair_send(ipair)%image > nimages) THEN
                  CPABORT("ECD| External translation send pair has an invalid translation index")
               END IF
               IF (task_list%atom_pair_send(ipair)%image == iimage) THEN
                  nreserve = nreserve + 1
                  reserve_rows(nreserve) = task_list%atom_pair_send(ipair)%row
                  reserve_cols(nreserve) = task_list%atom_pair_send(ipair)%col
               END IF
            END DO
         END IF
         DO imatrix = 1, nspin
            CALL dbcsr_init_p(matrices(imatrix, iimage)%matrix)
            CALL dbcsr_create(matrices(imatrix, iimage)%matrix, template=templates(imatrix)%matrix, &
                              matrix_type=dbcsr_type_symmetric, name=name)
            IF (nreserve > 0) THEN
               CALL dbcsr_reserve_blocks(matrices(imatrix, iimage)%matrix, reserve_rows(1:nreserve), &
                                         reserve_cols(1:nreserve))
            END IF
            CALL dbcsr_set(matrices(imatrix, iimage)%matrix, 0.0_dp)
         END DO
      END DO
      IF (npairs > 0) DEALLOCATE (reserve_cols, reserve_rows)

   END SUBROUTINE ecd_periodic_task_image_matrix_set_create

! **************************************************************************************************
!> \brief Build implicit-Gamma current Taylor coefficients through first spatial order.
!>
!> The local and GTH image kernels are ECD-owned.  Their zeroth coefficient is
!> checked against the folded complete current before the first coefficient is
!> made available to the finite-q transition provider.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_current_g4(qs_env, context, gamma_operators)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(INOUT), TARGET :: context
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: gamma_operators

      INTEGER                                            :: idir
      LOGICAL                                            :: gth_available, gth_enabled
      REAL(KIND=dp)                                      :: complete_norm_sq, first_imag_norm_sq, first_norm_sq, gth_norm_sq, &
                                                            gth_reference_norm_sq, gth_relative_residual, gth_residual_norm_sq, &
                                                            local_norm_sq, local_reference_norm_sq, local_relative_residual, &
                                                            local_residual_norm_sq, residual_norm_sq, zero_norm_sq, &
                                                            zero_relative_residual
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_zero
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: complete_current, current_residual, &
                                                            gth_current, local_current, position
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: current_taylor_imag, current_taylor_real, &
                                                            gth_current_taylor_imag, gth_current_taylor_real, &
                                                            local_moments_rs, local_momentum_rs, overlap_images
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: gth_diamagnetic_taylor_imag, &
                                                            gth_diamagnetic_taylor_real, gth_current_path_rs, &
                                                            local_moments_der_rs
      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: gth_diamagnetic_path_rs
      TYPE(ecd_periodic_gth_operators_type)              :: gth_operators
      TYPE(ecd_periodic_local_operators_type)            :: local_operators
      INTEGER, DIMENSION(:, :, :), POINTER               :: gth_cell_to_index, local_cell_to_index, &
                                                            translation_cell_to_index

      NULLIFY (complete_current, current_residual, current_taylor_imag, current_taylor_real, &
               gth_cell_to_index, gth_current, gth_current_path_rs, gth_current_taylor_imag, &
               gth_current_taylor_real, gth_diamagnetic_path_rs, gth_diamagnetic_taylor_imag, &
               gth_diamagnetic_taylor_real, local_cell_to_index, local_current, local_moments_der_rs, &
               local_moments_rs, local_momentum_rs, logger, translation_cell_to_index, overlap_images, position)
      xkp_zero = 0.0_dp
      CALL ecd_gamma_context_overlap_images_get(context, overlap_images, translation_cell_to_index)
      IF (ASSOCIATED(context%current_taylor_real) .OR. ASSOCIATED(context%current_taylor_imag) .OR. &
          ASSOCIATED(context%gth_current_taylor_real) .OR. ASSOCIATED(context%gth_current_taylor_imag)) THEN
         CPABORT("ECD| Stage G4 current Taylor coefficients are already allocated")
      END IF
      CALL ecd_gamma_zero_order_operators_get(gamma_operators, position, local_current, gth_current, complete_current, &
                                               gth_available)

      CALL ecd_periodic_local_operators_build(local_operators, qs_env, context%ref_point, 1)
      CALL ecd_periodic_local_operators_get(local_operators, local_cell_to_index, local_momentum_rs, &
                                            local_moments_der_rs, local_moments_rs)
      gth_enabled = gth_available
      IF (gth_enabled) THEN
          CALL ecd_periodic_gth_operators_build(gth_operators, qs_env, context%direction, context%ref_point, 1)
         CALL ecd_periodic_gth_operators_get(gth_operators, gth_cell_to_index, gth_current_path_rs, &
                                              gth_diamagnetic_path_rs)
      END IF

      CALL ecd_periodic_current_taylor_build(current_taylor_real, current_taylor_imag, &
                                             gth_current_taylor_real, gth_current_taylor_imag, &
                                             local_momentum_rs, local_moments_der_rs, local_moments_rs, &
                                             gth_current_path_rs, gth_enabled, overlap_images, xkp_zero, &
                                             translation_cell_to_index, local_cell_to_index, gth_cell_to_index, &
                                             context%sab_all, context%cell%hmat, context%direction, 1)
      IF (gth_enabled) THEN
         CALL ecd_periodic_gth_diamagnetic_taylor_build(gth_diamagnetic_taylor_real, &
                                                         gth_diamagnetic_taylor_imag, gth_diamagnetic_path_rs, &
                                                         xkp_zero, gth_cell_to_index, context%sab_all)
      END IF
      CALL ecd_periodic_ao_vector_create(current_residual, context%matrix_s(1)%matrix, context%sab_all, &
                                         'ECD G4 current zeroth residual', dbcsr_type_no_symmetry)

      logger => cp_get_default_logger()
      DO idir = 1, ncart
         ! Verify each production contribution separately so component errors cannot cancel in the total current.
         CALL dbcsr_set(current_residual(idir)%matrix, 0.0_dp)
         CALL dbcsr_add(current_residual(idir)%matrix, current_taylor_real(idir, 1)%matrix, 0.0_dp, 1.0_dp)
         CALL dbcsr_add(current_residual(idir)%matrix, complete_current(idir)%matrix, 1.0_dp, -1.0_dp)
         CALL dbcsr_dot(current_taylor_real(idir, 1)%matrix, current_taylor_real(idir, 1)%matrix, zero_norm_sq)
         CALL dbcsr_dot(complete_current(idir)%matrix, complete_current(idir)%matrix, complete_norm_sq)
         CALL dbcsr_dot(current_taylor_real(idir, 2)%matrix, current_taylor_real(idir, 2)%matrix, first_norm_sq)
         CALL dbcsr_dot(current_taylor_imag(idir, 2)%matrix, current_taylor_imag(idir, 2)%matrix, &
                        first_imag_norm_sq)
         CALL dbcsr_dot(current_residual(idir)%matrix, current_residual(idir)%matrix, residual_norm_sq)
         zero_relative_residual = 2.0_dp*SQRT(ABS(residual_norm_sq))/ &
                                  MAX(SQRT(TINY(1.0_dp)), SQRT(ABS(zero_norm_sq)) + SQRT(ABS(complete_norm_sq)))
         CALL dbcsr_set(current_residual(idir)%matrix, 0.0_dp)
         CALL dbcsr_add(current_residual(idir)%matrix, current_taylor_real(idir, 1)%matrix, 0.0_dp, 1.0_dp)
         CALL dbcsr_add(current_residual(idir)%matrix, gth_current_taylor_real(idir, 1)%matrix, 1.0_dp, -1.0_dp)
         CALL dbcsr_dot(current_residual(idir)%matrix, current_residual(idir)%matrix, local_norm_sq)
         CALL dbcsr_dot(local_current(idir)%matrix, local_current(idir)%matrix, local_reference_norm_sq)
         CALL dbcsr_add(current_residual(idir)%matrix, local_current(idir)%matrix, 1.0_dp, -1.0_dp)
         CALL dbcsr_dot(current_residual(idir)%matrix, current_residual(idir)%matrix, local_residual_norm_sq)
         local_relative_residual = 2.0_dp*SQRT(ABS(local_residual_norm_sq))/ &
                                   MAX(SQRT(TINY(1.0_dp)), SQRT(ABS(local_norm_sq)) + SQRT(ABS(local_reference_norm_sq)))
         gth_norm_sq = 0.0_dp
         gth_reference_norm_sq = 0.0_dp
         gth_residual_norm_sq = 0.0_dp
         gth_relative_residual = 0.0_dp
         IF (gth_enabled) THEN
            CALL dbcsr_set(current_residual(idir)%matrix, 0.0_dp)
            CALL dbcsr_add(current_residual(idir)%matrix, gth_current_taylor_real(idir, 1)%matrix, 0.0_dp, 1.0_dp)
            CALL dbcsr_dot(current_residual(idir)%matrix, current_residual(idir)%matrix, gth_norm_sq)
            CALL dbcsr_dot(gth_current(idir)%matrix, gth_current(idir)%matrix, gth_reference_norm_sq)
            CALL dbcsr_add(current_residual(idir)%matrix, gth_current(idir)%matrix, 1.0_dp, -1.0_dp)
            CALL dbcsr_dot(current_residual(idir)%matrix, current_residual(idir)%matrix, gth_residual_norm_sq)
            gth_relative_residual = 2.0_dp*SQRT(ABS(gth_residual_norm_sq))/ &
                                    MAX(SQRT(TINY(1.0_dp)), SQRT(ABS(gth_norm_sq)) + &
                                        SQRT(ABS(gth_reference_norm_sq)))
         END IF
         IF (cp_logger_get_default_io_unit(logger) > 0) THEN
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current zeroth norm component', idir, SQRT(ABS(zero_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current zeroth residual component', idir, SQRT(ABS(residual_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current zeroth relative residual component', idir, zero_relative_residual
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 local current zeroth residual component', idir, SQRT(ABS(local_residual_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 local current zeroth relative residual component', idir, local_relative_residual
            IF (gth_enabled) THEN
               WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
                  'ECD| Stage G4 GTH current zeroth residual component', idir, SQRT(ABS(gth_residual_norm_sq))
               WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
                  'ECD| Stage G4 GTH current zeroth relative residual component', idir, gth_relative_residual
            END IF
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current first real norm component', idir, SQRT(ABS(first_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current first imaginary norm component', idir, SQRT(ABS(first_imag_norm_sq))
         END IF
      END DO

      CALL dbcsr_deallocate_matrix_set(current_residual)
      context%current_taylor_imag => current_taylor_imag
      context%current_taylor_real => current_taylor_real
      context%gth_current_taylor_imag => gth_current_taylor_imag
      context%gth_current_taylor_real => gth_current_taylor_real
      context%gth_diamagnetic_taylor_imag => gth_diamagnetic_taylor_imag
      context%gth_diamagnetic_taylor_real => gth_diamagnetic_taylor_real
      NULLIFY (current_taylor_imag, current_taylor_real, gth_current_taylor_imag, gth_current_taylor_real, &
               gth_diamagnetic_taylor_imag, gth_diamagnetic_taylor_real)
      CALL ecd_periodic_gth_operators_release(gth_operators)
      CALL ecd_periodic_local_operators_release(local_operators)

   END SUBROUTINE ecd_periodic_implicit_gamma_current_g4

! **************************************************************************************************
!> \brief Assemble the first-order velocity-gauge causal response from retained Gamma providers.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_velocity_response_g10(qs_env, context, operators, transition_selection, &
                                                                 pseudopotential_commutator, energy_range, npoints, fwhm, &
                                                                 response_zero_output, response_first_output, &
                                                                 complete_zero_output, complete_first_output, &
                                                                 local_zero_output, local_first_output, &
                                                                 cross_zero_output, cross_first_output, &
                                                                 nonlocal_zero_output, nonlocal_first_output, &
                                                                 diamagnetic_kernel_output, transition_output, &
                                                                 gth_transition_output, &
                                                                 static_magnetic_transition_output, &
                                                                 gth_available_output, &
                                                                 diamagnetic_kernel_reference)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET   :: context
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: operators
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: transition_selection
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      INTEGER, INTENT(IN)                                :: npoints
      REAL(KIND=dp), INTENT(IN)                          :: fwhm
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: response_first_output, response_zero_output
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: complete_first_output, complete_zero_output, &
                                                                          cross_first_output, cross_zero_output, &
                                                                          local_first_output, local_zero_output, &
                                                                          nonlocal_first_output, nonlocal_zero_output
      REAL(KIND=dp), DIMENSION(ncart, ncart), INTENT(OUT) :: diamagnetic_kernel_output
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: transition_output
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: gth_transition_output
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: static_magnetic_transition_output
      LOGICAL, INTENT(OUT)                               :: gth_available_output
      REAL(KIND=dp), DIMENSION(ncart, ncart), INTENT(IN), OPTIONAL :: diamagnetic_kernel_reference

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: current_first_mo, current_first_rotated, &
                                                            current_zero_mo, current_zero_rotated, &
                                                            gth_first_mo, gth_first_rotated, gth_zero_mo, &
                                                            gth_zero_rotated, state_first, &
                                                            static_magnetic_rotated, transition_first
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: current_first_components, current_zero_components, &
                                                            cross_counter, cross_first_counter, cross_first_resonant, &
                                                            cross_resonant, local_counter, local_first_counter, &
                                                            local_first_resonant, local_resonant, nonlocal_counter, &
                                                            nonlocal_first_counter, nonlocal_first_resonant, &
                                                            nonlocal_resonant, gth_first_components, &
                                                            gth_zero_components, static_magnetic_components
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: gth_transition_taylor, local_transition_taylor, &
                                                            transition_taylor
      COMPLEX(KIND=dp), DIMENSION(ncart, ncart, 0:1)     :: counter_taylor, resonant_taylor
      INTEGER                                            :: iband, idir, ipoint, ispin, itransition, ivirt, nmo, &
                                                            ntransitions, ntransitions_total
      LOGICAL                                            :: gth_available
      REAL(KIND=dp)                                      :: antisymmetric_first_norm, cell_volume, gap, gap_derivative, &
                                                            occupation_difference, photon_energy, prefactor, &
                                                            response_counter_norm, response_first_counter_norm, &
                                                            response_first_norm, response_first_resonant_norm, &
                                                            response_norm, gth_diamagnetic_zero_relative_residual, &
                                                            gth_diamagnetic_zero_residual, response_resonant_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: gap_taylor
      REAL(KIND=dp), DIMENSION(ncart, ncart)             :: gth_diamagnetic_kernel, local_diamagnetic_kernel
      COMPLEX(KIND=dp), DIMENSION(ncart, ncart)          :: gth_diamagnetic_first
      REAL(KIND=dp), DIMENSION(ncart, ncart)             :: gth_diamagnetic_path_zero
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_fm_type)                                   :: static_magnetic_mo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: static_magnetic_buffer
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: complete_current, gth_current, local_current, position, &
                                                            static_magnetic_diagnostic

      NULLIFY (complete_current, gth_current, local_current, logger, position, static_magnetic_diagnostic)
      IF (.NOT. ASSOCIATED(context%current_taylor_real) .OR. .NOT. ASSOCIATED(context%current_taylor_imag) .OR. &
          .NOT. ASSOCIATED(context%gth_current_taylor_real) .OR. &
          .NOT. ASSOCIATED(context%gth_current_taylor_imag)) THEN
         CPABORT("ECD| Stage G10 requires retained current Taylor coefficients")
      END IF
      IF (SIZE(context%current_taylor_real, 1) /= ncart .OR. SIZE(context%current_taylor_imag, 1) /= ncart .OR. &
          SIZE(context%gth_current_taylor_real, 1) /= ncart .OR. &
          SIZE(context%gth_current_taylor_imag, 1) /= ncart .OR. &
          SIZE(context%current_taylor_real, 2) < 2 .OR. SIZE(context%current_taylor_imag, 2) < 2 .OR. &
          SIZE(context%gth_current_taylor_real, 2) < 2 .OR. SIZE(context%gth_current_taylor_imag, 2) < 2) THEN
         CPABORT("ECD| Stage G10 current Taylor coefficient dimensions are inconsistent")
      END IF
      IF (fwhm <= 0.0_dp) THEN
         CPABORT("ECD| Stage G10 requires a positive FWHM")
      END IF
      cell_volume = ABS(context%cell%deth)
      IF (cell_volume <= 0.0_dp) THEN
         CPABORT("ECD| Stage G10 requires a positive cell volume")
      END IF
      CALL ecd_gamma_zero_order_operators_get(operators, position, local_current, gth_current, complete_current, &
                                               gth_available, static_magnetic_diagnostic)
      gth_available_output = gth_available

      ALLOCATE (cross_counter(ncart, ncart, npoints), cross_first_counter(ncart, ncart, npoints), &
                cross_first_resonant(ncart, ncart, npoints), cross_resonant(ncart, ncart, npoints), &
                local_counter(ncart, ncart, npoints), local_first_counter(ncart, ncart, npoints), &
                local_first_resonant(ncart, ncart, npoints), local_resonant(ncart, ncart, npoints), &
                nonlocal_counter(ncart, ncart, npoints), nonlocal_first_counter(ncart, ncart, npoints), &
                nonlocal_first_resonant(ncart, ncart, npoints), nonlocal_resonant(ncart, ncart, npoints))
      cross_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      cross_first_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      cross_first_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      cross_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      local_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      local_first_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      local_first_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      local_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      nonlocal_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      nonlocal_first_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      nonlocal_first_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      nonlocal_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ALLOCATE (transition_output(ncart, 0:1, transition_selection%selected_count))
      transition_output = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ALLOCATE (gth_transition_output(ncart, 0:1, transition_selection%selected_count))
      gth_transition_output = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ALLOCATE (static_magnetic_transition_output(ncart, transition_selection%selected_count))
      static_magnetic_transition_output = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      logger => cp_get_default_logger()
      ntransitions_total = 0
      DO ispin = 1, context%nspin
         nmo = context%spin(ispin)%nmo
         IF (nmo < 1 .OR. .NOT. ASSOCIATED(context%spin(ispin)%eigenvalue_first) .OR. &
             .NOT. ASSOCIATED(context%spin(ispin)%eigenvalue_reference) .OR. &
             .NOT. ASSOCIATED(context%spin(ispin)%mo_coeff) .OR. &
             .NOT. ASSOCIATED(context%spin(ispin)%manifolds) .OR. &
             .NOT. ALLOCATED(context%spin(ispin)%band_rotation) .OR. &
             .NOT. ALLOCATED(context%spin(ispin)%coefficient_first)) THEN
            CPABORT("ECD| Stage G10 state-response data are incomplete")
         END IF
         IF (SIZE(context%spin(ispin)%eigenvalue_first) /= nmo .OR. &
             SIZE(context%spin(ispin)%eigenvalue_reference) /= nmo .OR. &
             SIZE(context%spin(ispin)%band_rotation, 1) /= nmo .OR. &
             SIZE(context%spin(ispin)%band_rotation, 2) /= nmo .OR. &
             SIZE(context%spin(ispin)%coefficient_first, 1) /= nmo .OR. &
             SIZE(context%spin(ispin)%coefficient_first, 2) /= nmo) THEN
            CPABORT("ECD| Stage G10 retained state-response dimensions are inconsistent")
         END IF
         ALLOCATE (current_first_components(nmo, nmo, ncart), current_zero_components(nmo, nmo, ncart), &
                   gth_first_components(nmo, nmo, ncart), gth_zero_components(nmo, nmo, ncart), &
                   static_magnetic_components(nmo, nmo, ncart))
         DO idir = 1, ncart
            CALL ecd_project_complex_operator_to_mo(context%current_taylor_real(idir, 1), &
                                                    context%current_taylor_imag(idir, 1), &
                                                    context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                    current_zero_mo)
            CALL ecd_project_complex_operator_to_mo(context%current_taylor_real(idir, 2), &
                                                    context%current_taylor_imag(idir, 2), &
                                                    context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                    current_first_mo)
            CALL ecd_project_complex_operator_to_mo(context%gth_current_taylor_real(idir, 1), &
                                                    context%gth_current_taylor_imag(idir, 1), &
                                                    context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                    gth_zero_mo)
            CALL ecd_project_complex_operator_to_mo(context%gth_current_taylor_real(idir, 2), &
                                                    context%gth_current_taylor_imag(idir, 2), &
                                                    context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                    gth_first_mo)
            ALLOCATE (current_first_rotated(nmo, nmo), current_zero_rotated(nmo, nmo), &
                      gth_first_rotated(nmo, nmo), gth_zero_rotated(nmo, nmo), &
                      state_first(nmo, nmo), transition_first(nmo, nmo))
            current_zero_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                           MATMUL(current_zero_mo, context%spin(ispin)%band_rotation))
            current_first_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                            MATMUL(current_first_mo, context%spin(ispin)%band_rotation))
            gth_zero_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                      MATMUL(gth_zero_mo, context%spin(ispin)%band_rotation))
            gth_first_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                       MATMUL(gth_first_mo, context%spin(ispin)%band_rotation))
            state_first = 0.5_dp*(MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%coefficient_first)), &
                                          current_zero_rotated) - &
                                    MATMUL(current_zero_rotated, context%spin(ispin)%coefficient_first))
            transition_first = current_first_rotated + state_first
            current_zero_components(:, :, idir) = current_zero_rotated
            current_first_components(:, :, idir) = transition_first
            state_first = 0.5_dp*(MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%coefficient_first)), &
                                         gth_zero_rotated) - &
                                   MATMUL(gth_zero_rotated, context%spin(ispin)%coefficient_first))
            gth_first_components(:, :, idir) = gth_first_rotated + state_first
            gth_zero_components(:, :, idir) = gth_zero_rotated
            CALL ecd_project_real_operator_to_mo(static_magnetic_diagnostic(idir), context%spin(ispin)%mo_coeff, &
                                                 context%spin(ispin)%nao, nmo, static_magnetic_mo)
            ALLOCATE (static_magnetic_buffer(nmo, nmo), static_magnetic_rotated(nmo, nmo))
            CALL cp_fm_get_submatrix(static_magnetic_mo, static_magnetic_buffer)
            static_magnetic_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                               MATMUL(CMPLX(static_magnetic_buffer, 0.0_dp, KIND=dp), &
                                                      context%spin(ispin)%band_rotation))
            static_magnetic_components(:, :, idir) = static_magnetic_rotated
            DEALLOCATE (static_magnetic_buffer, static_magnetic_rotated)
            CALL cp_fm_release(static_magnetic_mo)
            DEALLOCATE (current_first_mo, current_first_rotated, current_zero_mo, current_zero_rotated, &
                        gth_first_mo, gth_first_rotated, gth_zero_mo, gth_zero_rotated, state_first, &
                        transition_first)
         END DO
         ALLOCATE (gap_taylor(0:1), gth_transition_taylor(ncart, 0:1), &
                   local_transition_taylor(ncart, 0:1), transition_taylor(ncart, 0:1))
         ntransitions = 0
         DO itransition = 1, transition_selection%selected_count
            IF (transition_selection%ispin(itransition) /= ispin) CYCLE
            iband = transition_selection%iocc(itransition)
            ivirt = transition_selection%ivirt(itransition)
            gap = transition_selection%gap(itransition)
            occupation_difference = transition_selection%occupation_difference(itransition)
            gap_derivative = 0.5_dp*(context%spin(ispin)%eigenvalue_first(ivirt) + &
                                     context%spin(ispin)%eigenvalue_first(iband))
            gap_taylor(0) = gap
            gap_taylor(1) = gap_derivative
            DO idir = 1, ncart
               transition_taylor(idir, 0) = current_zero_components(ivirt, iband, idir)
               transition_taylor(idir, 1) = current_first_components(ivirt, iband, idir)
               gth_transition_taylor(idir, 0) = gth_zero_components(ivirt, iband, idir)
               gth_transition_taylor(idir, 1) = gth_first_components(ivirt, iband, idir)
            END DO
            local_transition_taylor = transition_taylor - gth_transition_taylor
            transition_output(:, :, itransition) = transition_taylor
            gth_transition_output(:, :, itransition) = gth_transition_taylor
            static_magnetic_transition_output(:, itransition) = &
               0.5_dp*static_magnetic_components(ivirt, iband, :)
            ntransitions = ntransitions + 1
            CALL ecd_periodic_response_add_transition_pair(local_transition_taylor, local_transition_taylor, &
                                                            gap_taylor, occupation_difference, cell_volume, &
                                                            energy_range, fwhm, local_resonant, local_counter, &
                                                            local_first_resonant, local_first_counter, &
                                                            resonant_taylor, counter_taylor)
            CALL ecd_periodic_response_add_transition_pair(local_transition_taylor, gth_transition_taylor, &
                                                            gap_taylor, occupation_difference, cell_volume, &
                                                            energy_range, fwhm, cross_resonant, cross_counter, &
                                                            cross_first_resonant, cross_first_counter, &
                                                            resonant_taylor, counter_taylor)
            CALL ecd_periodic_response_add_transition_pair(gth_transition_taylor, local_transition_taylor, &
                                                            gap_taylor, occupation_difference, cell_volume, &
                                                            energy_range, fwhm, cross_resonant, cross_counter, &
                                                            cross_first_resonant, cross_first_counter, &
                                                            resonant_taylor, counter_taylor)
            CALL ecd_periodic_response_add_transition_pair(gth_transition_taylor, gth_transition_taylor, &
                                                            gap_taylor, occupation_difference, cell_volume, &
                                                            energy_range, fwhm, nonlocal_resonant, nonlocal_counter, &
                                                            nonlocal_first_resonant, nonlocal_first_counter, &
                                                            resonant_taylor, counter_taylor)
         END DO
         ntransitions_total = ntransitions_total + ntransitions
         DEALLOCATE (current_first_components, current_zero_components, gth_first_components, &
                     gth_zero_components, static_magnetic_components, gap_taylor, gth_transition_taylor, &
                     local_transition_taylor, transition_taylor)
      END DO
      IF (ntransitions_total /= transition_selection%selected_count) THEN
         CPABORT('ECD| Stage G10 transition selection is inconsistent with the Gamma spin space')
      END IF

      IF (PRESENT(diamagnetic_kernel_reference)) THEN
         diamagnetic_kernel_output = diamagnetic_kernel_reference
         local_diamagnetic_kernel = 0.0_dp
         DO ispin = 1, context%nspin
            DO idir = 1, ncart
               local_diamagnetic_kernel(idir, idir) = local_diamagnetic_kernel(idir, idir) - &
                                                      SUM(context%spin(ispin)%occupations)
            END DO
         END DO
         gth_diamagnetic_kernel = diamagnetic_kernel_output - local_diamagnetic_kernel
         IF (.NOT. gth_available .AND. MAXVAL(ABS(gth_diamagnetic_kernel)) > 100.0_dp*EPSILON(1.0_dp)) THEN
            CPABORT('ECD| Cached diamagnetic kernel contains a GTH contribution without GTH current data')
         END IF
      ELSE
         CALL ecd_periodic_implicit_gamma_diamagnetic_kernel(qs_env, context, gth_available, &
                                                              local_diamagnetic_kernel, gth_diamagnetic_kernel, &
                                                              diamagnetic_kernel_output)
      END IF
      gth_diamagnetic_first = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      gth_diamagnetic_path_zero = 0.0_dp
      gth_diamagnetic_zero_residual = 0.0_dp
      gth_diamagnetic_zero_relative_residual = 0.0_dp
      IF (gth_available) THEN
         CALL ecd_periodic_implicit_gamma_gth_diamagnetic_taylor(context, gth_diamagnetic_path_zero, &
                                                                   gth_diamagnetic_first)
         gth_diamagnetic_zero_residual = SQRT(SUM((gth_diamagnetic_path_zero - gth_diamagnetic_kernel)**2))
         gth_diamagnetic_zero_relative_residual = 2.0_dp*gth_diamagnetic_zero_residual/ &
                                                   MAX(SQRT(TINY(1.0_dp)), SQRT(SUM(gth_diamagnetic_path_zero**2)) + &
                                                       SQRT(SUM(gth_diamagnetic_kernel**2)))
      END IF
      CALL ecd_periodic_current_fsum_closure_g11(transition_selection, transition_output, &
                                                  diamagnetic_kernel_output)
      ALLOCATE (complete_zero_output(ncart, ncart, npoints), complete_first_output(ncart, ncart, npoints), &
                local_zero_output(ncart, ncart, npoints), local_first_output(ncart, ncart, npoints), &
                cross_zero_output(ncart, ncart, npoints), cross_first_output(ncart, ncart, npoints), &
                nonlocal_zero_output(ncart, ncart, npoints), nonlocal_first_output(ncart, ncart, npoints))
      local_zero_output = local_resonant + local_counter
      local_first_output = local_first_resonant + local_first_counter
      cross_zero_output = cross_resonant + cross_counter
      cross_first_output = cross_first_resonant + cross_first_counter
      nonlocal_zero_output = nonlocal_resonant + nonlocal_counter
      nonlocal_first_output = nonlocal_first_resonant + nonlocal_first_counter
      DO ipoint = 1, npoints
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                         (energy_range(2) - energy_range(1))/REAL(npoints - 1, KIND=dp)
         prefactor = 2.0_dp*twopi/(cell_volume*photon_energy**2)
         local_zero_output(:, :, ipoint) = local_zero_output(:, :, ipoint) + &
                                           CMPLX(prefactor*local_diamagnetic_kernel, 0.0_dp, KIND=dp)
         nonlocal_zero_output(:, :, ipoint) = nonlocal_zero_output(:, :, ipoint) + &
                                              CMPLX(prefactor*gth_diamagnetic_kernel, 0.0_dp, KIND=dp)
         nonlocal_first_output(:, :, ipoint) = nonlocal_first_output(:, :, ipoint) + &
                                               prefactor*gth_diamagnetic_first
      END DO
      complete_zero_output = local_zero_output + cross_zero_output + nonlocal_zero_output
      complete_first_output = local_first_output + cross_first_output + nonlocal_first_output

      ALLOCATE (response_zero_output(ncart, ncart, npoints), response_first_output(ncart, ncart, npoints))
      response_zero_output = complete_zero_output
      response_first_output = complete_first_output
      IF (.NOT. pseudopotential_commutator) THEN
         response_zero_output = local_zero_output
         response_first_output = local_first_output
      END IF

      ipoint = (npoints + 1)/2
      photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                     (energy_range(2) - energy_range(1))/REAL(npoints - 1, KIND=dp)
      response_resonant_norm = SQRT(SUM(ABS(local_resonant(:, :, ipoint) + cross_resonant(:, :, ipoint) + &
                                             nonlocal_resonant(:, :, ipoint))**2))
      response_counter_norm = SQRT(SUM(ABS(local_counter(:, :, ipoint) + cross_counter(:, :, ipoint) + &
                                            nonlocal_counter(:, :, ipoint))**2))
      response_first_resonant_norm = SQRT(SUM(ABS(local_first_resonant(:, :, ipoint) + &
                                                   cross_first_resonant(:, :, ipoint) + &
                                                   nonlocal_first_resonant(:, :, ipoint))**2))
      response_first_counter_norm = SQRT(SUM(ABS(local_first_counter(:, :, ipoint) + &
                                                  cross_first_counter(:, :, ipoint) + &
                                                  nonlocal_first_counter(:, :, ipoint))**2))
      response_norm = SQRT(SUM(ABS(complete_zero_output(:, :, ipoint))**2))
      response_first_norm = SQRT(SUM(ABS(complete_first_output(:, :, ipoint))**2))
      antisymmetric_first_norm = 0.5_dp*SQRT(SUM(ABS(complete_first_output(:, :, ipoint) - &
                                                       TRANSPOSE(complete_first_output(:, :, ipoint)))**2))
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,I10)') &
            'ECD| Stage G10 velocity response grid points', npoints
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 midpoint photon energy [hartree]', photon_energy
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 midpoint resonant zeroth-order norm', response_resonant_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 midpoint counter zeroth-order norm', response_counter_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 midpoint full zeroth-order norm', response_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 midpoint resonant first-order norm', response_first_resonant_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 midpoint counter first-order norm', response_first_counter_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 midpoint full first-order norm', response_first_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 midpoint antisymmetric first-order norm', antisymmetric_first_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 GTH diamagnetic zero-order residual', gth_diamagnetic_zero_residual
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 GTH diamagnetic zero-order relative residual', gth_diamagnetic_zero_relative_residual
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G10 GTH diamagnetic first-order norm', SQRT(SUM(ABS(gth_diamagnetic_first)**2))
      END IF
      DEALLOCATE (cross_counter, cross_first_counter, cross_first_resonant, cross_resonant, local_counter, &
                  local_first_counter, local_first_resonant, local_resonant, nonlocal_counter, &
                  nonlocal_first_counter, nonlocal_first_resonant, nonlocal_resonant)

   END SUBROUTINE ecd_periodic_implicit_gamma_velocity_response_g10

! **************************************************************************************************
!> \brief Diagnose the q=0 complete-current f-sum closure against the diamagnetic kernel.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_current_fsum_closure_g11(transition_selection, transition_output, &
                                                     diamagnetic_kernel)
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: transition_selection
      COMPLEX(KIND=dp), DIMENSION(:, 0:, :), INTENT(IN) :: transition_output
      REAL(KIND=dp), DIMENSION(ncart, ncart), INTENT(IN) :: diamagnetic_kernel

      COMPLEX(KIND=dp), DIMENSION(ncart, ncart) :: closure, fsum, hermiticity
      INTEGER :: ia, ib, itransition
      LOGICAL :: complete_window
      REAL(KIND=dp) :: closure_norm, diamagnetic_norm, fsum_norm, hermiticity_norm, relative_residual
      TYPE(cp_logger_type), POINTER :: logger

      IF (SIZE(transition_output, 1) /= ncart .OR. SIZE(transition_output, 2) < 2 .OR. &
          SIZE(transition_output, 3) /= transition_selection%selected_count) THEN
         CPABORT('ECD| Stage G11 complete-current transition storage is inconsistent')
      END IF
      IF (ANY(transition_selection%gap <= SQRT(TINY(1.0_dp)))) THEN
         CPABORT('ECD| Stage G11 current f-sum received a nonpositive transition gap')
      END IF

      fsum = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO itransition = 1, transition_selection%selected_count
         DO ia = 1, ncart
            DO ib = 1, ncart
               fsum(ia, ib) = fsum(ia, ib) + 2.0_dp*transition_selection%occupation_difference(itransition)* &
                              CONJG(transition_output(ia, 0, itransition))* &
                              transition_output(ib, 0, itransition)/transition_selection%gap(itransition)
            END DO
         END DO
      END DO
      closure = fsum + CMPLX(diamagnetic_kernel, 0.0_dp, KIND=dp)
      hermiticity = fsum - CONJG(TRANSPOSE(fsum))
      fsum_norm = SQRT(SUM(ABS(fsum)**2))
      diamagnetic_norm = SQRT(SUM(diamagnetic_kernel**2))
      closure_norm = SQRT(SUM(ABS(closure)**2))
      hermiticity_norm = SQRT(SUM(ABS(hermiticity)**2))
      relative_residual = 2.0_dp*closure_norm/MAX(SQRT(TINY(1.0_dp)), fsum_norm + diamagnetic_norm)
      complete_window = transition_selection%requested_count == -1 .AND. &
                        ALL(transition_selection%band_range == -1) .AND. &
                        transition_selection%selected_count == transition_selection%candidate_count

      logger => cp_get_default_logger()
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,I10)') &
            'ECD| Stage G11 current f-sum selected transitions', transition_selection%selected_count
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,L1)') &
            'ECD| Stage G11 current f-sum complete window', complete_window
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G11 current f-sum norm', fsum_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G11 negative diamagnetic norm', diamagnetic_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G11 current f-sum closure residual', closure_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G11 current f-sum relative closure residual', relative_residual
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G11 current f-sum Hermiticity residual', hermiticity_norm
      END IF

   END SUBROUTINE ecd_periodic_current_fsum_closure_g11

! **************************************************************************************************
!> \brief Build direct AO position transitions for isolated-molecule-limit diagnostics.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_direct_position_transitions(context, operators, transition_selection, &
                                                                       transition_output)
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET   :: context
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: operators
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: transition_selection
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: transition_output

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: length_components
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: position_rotated
      INTEGER                                            :: iband, idir, ispin, itransition, ivirt, nmo, &
                                                            ntransitions, ntransitions_total
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: position_buffer
      TYPE(cp_fm_type)                                   :: position_mo
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: complete_current, gth_current, local_current, position
      NULLIFY (complete_current, gth_current, local_current, position)
      CALL ecd_gamma_zero_order_operators_get(operators, position, local_current, gth_current, complete_current)
      ALLOCATE (transition_output(ncart, transition_selection%selected_count))
      transition_output = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      ntransitions_total = 0
      DO ispin = 1, context%nspin
         nmo = context%spin(ispin)%nmo
         IF (nmo < 1 .OR. .NOT. ASSOCIATED(context%spin(ispin)%mo_coeff) .OR. &
             .NOT. ALLOCATED(context%spin(ispin)%band_rotation)) THEN
            CPABORT('ECD| Direct-position transition data are incomplete')
         END IF
         IF (SIZE(context%spin(ispin)%band_rotation, 1) /= nmo .OR. &
             SIZE(context%spin(ispin)%band_rotation, 2) /= nmo) THEN
            CPABORT('ECD| Direct-position transition rotation dimensions are inconsistent')
         END IF
         ALLOCATE (length_components(nmo, nmo, ncart))
         DO idir = 1, ncart
            CALL ecd_project_real_operator_to_mo(position(idir), context%spin(ispin)%mo_coeff, &
                                                 context%spin(ispin)%nao, nmo, position_mo)
            ALLOCATE (position_buffer(nmo, nmo), position_rotated(nmo, nmo))
            CALL cp_fm_get_submatrix(position_mo, position_buffer)
            position_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                      MATMUL(CMPLX(position_buffer, 0.0_dp, KIND=dp), &
                                             context%spin(ispin)%band_rotation))
            length_components(:, :, idir) = position_rotated
            DEALLOCATE (position_buffer, position_rotated)
            CALL cp_fm_release(position_mo)
         END DO
         ntransitions = 0
         DO itransition = 1, transition_selection%selected_count
            IF (transition_selection%ispin(itransition) /= ispin) CYCLE
            iband = transition_selection%iocc(itransition)
            ivirt = transition_selection%ivirt(itransition)
            DO idir = 1, ncart
               transition_output(idir, itransition) = length_components(ivirt, iband, idir)
            END DO
            ntransitions = ntransitions + 1
         END DO
         ntransitions_total = ntransitions_total + ntransitions
         DEALLOCATE (length_components)
      END DO
      IF (ntransitions_total /= transition_selection%selected_count) THEN
         CPABORT('ECD| Direct-position transition selection is inconsistent with the Gamma spin space')
      END IF

   END SUBROUTINE ecd_periodic_implicit_gamma_direct_position_transitions

! **************************************************************************************************
!> \brief Add one transition Taylor series to a first-order retarded response grid.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_response_add_transition(transition_taylor, gap_taylor, occupation_difference, &
                                                   cell_volume, energy_range, fwhm, response_resonant, &
                                                   response_counter, response_first_resonant, response_first_counter, &
                                                   resonant_taylor, counter_taylor)
      COMPLEX(KIND=dp), DIMENSION(:, 0:), INTENT(IN)     :: transition_taylor
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT) :: response_counter, response_first_counter, &
                                                             response_first_resonant, response_resonant
      COMPLEX(KIND=dp), DIMENSION(:, :, 0:), INTENT(OUT) :: counter_taylor, resonant_taylor
      INTEGER                                            :: ipoint, ncart_local, npoints
      REAL(KIND=dp), DIMENSION(0:), INTENT(IN)           :: gap_taylor
      REAL(KIND=dp), INTENT(IN)                          :: cell_volume, fwhm, occupation_difference
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      REAL(KIND=dp)                                      :: photon_energy, prefactor

      ncart_local = SIZE(transition_taylor, 1)
      npoints = SIZE(response_resonant, 3)
      IF (ncart_local /= ncart .OR. LBOUND(transition_taylor, 2) /= 0 .OR. &
          UBOUND(transition_taylor, 2) /= 1 .OR. LBOUND(gap_taylor, 1) /= 0 .OR. &
          UBOUND(gap_taylor, 1) /= 1) THEN
         CPABORT("ECD| Retarded response addition requires first-order transition and gap Taylor coefficients")
      END IF
      IF (npoints < 2 .OR. cell_volume <= 0.0_dp .OR. fwhm <= 0.0_dp .OR. &
          energy_range(1) <= 0.0_dp .OR. energy_range(1) >= energy_range(2)) THEN
         CPABORT("ECD| Retarded response addition received an invalid response grid")
      END IF
      IF (SIZE(response_resonant, 1) /= ncart .OR. SIZE(response_resonant, 2) /= ncart .OR. &
          SIZE(response_counter, 1) /= ncart .OR. SIZE(response_counter, 2) /= ncart .OR. &
          SIZE(response_first_resonant, 1) /= ncart .OR. SIZE(response_first_resonant, 2) /= ncart .OR. &
          SIZE(response_first_counter, 1) /= ncart .OR. SIZE(response_first_counter, 2) /= ncart .OR. &
          SIZE(response_counter, 3) /= npoints .OR. SIZE(response_first_resonant, 3) /= npoints .OR. &
          SIZE(response_first_counter, 3) /= npoints) THEN
         CPABORT("ECD| Retarded response addition received incompatible response tensors")
      END IF
      IF (SIZE(resonant_taylor, 1) /= ncart .OR. SIZE(resonant_taylor, 2) /= ncart .OR. &
          SIZE(counter_taylor, 1) /= ncart .OR. SIZE(counter_taylor, 2) /= ncart .OR. &
          LBOUND(resonant_taylor, 3) /= 0 .OR. UBOUND(resonant_taylor, 3) /= 1 .OR. &
          LBOUND(counter_taylor, 3) /= 0 .OR. UBOUND(counter_taylor, 3) /= 1) THEN
         CPABORT("ECD| Retarded response addition received incompatible Taylor work tensors")
      END IF

      DO ipoint = 1, npoints
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                         (energy_range(2) - energy_range(1))/REAL(npoints - 1, KIND=dp)
         prefactor = -2.0_dp*twopi/(cell_volume*photon_energy**2)
         CALL ecd_taylor_retarded_response(transition_taylor, gap_taylor, photon_energy, 0.5_dp*fwhm, &
                                            resonant_taylor, counter_taylor)
         response_resonant(:, :, ipoint) = response_resonant(:, :, ipoint) + &
                                            prefactor*occupation_difference*resonant_taylor(:, :, 0)
         response_counter(:, :, ipoint) = response_counter(:, :, ipoint) + &
                                           prefactor*occupation_difference*counter_taylor(:, :, 0)
         response_first_resonant(:, :, ipoint) = response_first_resonant(:, :, ipoint) + &
                                                  prefactor*occupation_difference*resonant_taylor(:, :, 1)
         response_first_counter(:, :, ipoint) = response_first_counter(:, :, ipoint) + &
                                                 prefactor*occupation_difference*counter_taylor(:, :, 1)
      END DO

   END SUBROUTINE ecd_periodic_response_add_transition

! **************************************************************************************************
!> \brief Add one bilinear current-transition Taylor series to a retarded dielectric grid.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_response_add_transition_pair(left_transition_taylor, right_transition_taylor, &
                                                         gap_taylor, occupation_difference, cell_volume, &
                                                         energy_range, fwhm, response_resonant, response_counter, &
                                                         response_first_resonant, response_first_counter, &
                                                         resonant_taylor, counter_taylor)
      COMPLEX(KIND=dp), DIMENSION(:, 0:), INTENT(IN)     :: left_transition_taylor, right_transition_taylor
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT) :: response_counter, response_first_counter, &
                                                             response_first_resonant, response_resonant
      COMPLEX(KIND=dp), DIMENSION(:, :, 0:), INTENT(OUT) :: counter_taylor, resonant_taylor
      INTEGER                                            :: ipoint, npoints
      REAL(KIND=dp), DIMENSION(0:), INTENT(IN)           :: gap_taylor
      REAL(KIND=dp), INTENT(IN)                          :: cell_volume, fwhm, occupation_difference
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      REAL(KIND=dp)                                      :: photon_energy, prefactor

      npoints = SIZE(response_resonant, 3)
      IF (SIZE(left_transition_taylor, 1) /= ncart .OR. SIZE(right_transition_taylor, 1) /= ncart .OR. &
          LBOUND(left_transition_taylor, 2) /= 0 .OR. LBOUND(right_transition_taylor, 2) /= 0 .OR. &
          UBOUND(left_transition_taylor, 2) /= 1 .OR. UBOUND(right_transition_taylor, 2) /= 1 .OR. &
          LBOUND(gap_taylor, 1) /= 0 .OR. UBOUND(gap_taylor, 1) /= 1) THEN
         CPABORT('ECD| Bilinear retarded response requires first-order transition and gap Taylor coefficients')
      END IF
      IF (npoints < 2 .OR. cell_volume <= 0.0_dp .OR. fwhm <= 0.0_dp .OR. &
          energy_range(1) <= 0.0_dp .OR. energy_range(1) >= energy_range(2)) THEN
         CPABORT('ECD| Bilinear retarded response received an invalid response grid')
      END IF
      IF (SIZE(response_resonant, 1) /= ncart .OR. SIZE(response_resonant, 2) /= ncart .OR. &
          SIZE(response_counter, 1) /= ncart .OR. SIZE(response_counter, 2) /= ncart .OR. &
          SIZE(response_first_resonant, 1) /= ncart .OR. SIZE(response_first_resonant, 2) /= ncart .OR. &
          SIZE(response_first_counter, 1) /= ncart .OR. SIZE(response_first_counter, 2) /= ncart .OR. &
          SIZE(response_counter, 3) /= npoints .OR. SIZE(response_first_resonant, 3) /= npoints .OR. &
          SIZE(response_first_counter, 3) /= npoints) THEN
         CPABORT('ECD| Bilinear retarded response received incompatible response tensors')
      END IF
      IF (SIZE(resonant_taylor, 1) /= ncart .OR. SIZE(resonant_taylor, 2) /= ncart .OR. &
          SIZE(counter_taylor, 1) /= ncart .OR. SIZE(counter_taylor, 2) /= ncart .OR. &
          LBOUND(resonant_taylor, 3) /= 0 .OR. UBOUND(resonant_taylor, 3) /= 1 .OR. &
          LBOUND(counter_taylor, 3) /= 0 .OR. UBOUND(counter_taylor, 3) /= 1) THEN
         CPABORT('ECD| Bilinear retarded response received incompatible Taylor work tensors')
      END IF

      DO ipoint = 1, npoints
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                         (energy_range(2) - energy_range(1))/REAL(npoints - 1, KIND=dp)
         prefactor = -2.0_dp*twopi/(cell_volume*photon_energy**2)
         CALL ecd_taylor_retarded_response_pair(left_transition_taylor, right_transition_taylor, gap_taylor, &
                                                 photon_energy, 0.5_dp*fwhm, resonant_taylor, counter_taylor)
         response_resonant(:, :, ipoint) = response_resonant(:, :, ipoint) + &
                                            prefactor*occupation_difference*resonant_taylor(:, :, 0)
         response_counter(:, :, ipoint) = response_counter(:, :, ipoint) + &
                                           prefactor*occupation_difference*counter_taylor(:, :, 0)
         response_first_resonant(:, :, ipoint) = response_first_resonant(:, :, ipoint) + &
                                                  prefactor*occupation_difference*resonant_taylor(:, :, 1)
         response_first_counter(:, :, ipoint) = response_first_counter(:, :, ipoint) + &
                                                 prefactor*occupation_difference*counter_taylor(:, :, 1)
      END DO

   END SUBROUTINE ecd_periodic_response_add_transition_pair

! **************************************************************************************************
!> \brief Validate that the q=0 response is independent of the Cartesian q-derivative direction.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_response_zero_check(reference_response, candidate_response, representation_label)
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: candidate_response, reference_response
      CHARACTER(LEN=*), INTENT(IN)                       :: representation_label

      REAL(KIND=dp), PARAMETER                           :: relative_tolerance = 1.0E-8_dp

      REAL(KIND=dp)                                      :: candidate_norm, reference_norm, relative_residual, residual_norm

      IF (SIZE(reference_response, 1) /= ncart .OR. SIZE(reference_response, 2) /= ncart .OR. &
          SIZE(candidate_response, 1) /= ncart .OR. SIZE(candidate_response, 2) /= ncart .OR. &
          SIZE(reference_response, 3) < 2 .OR. SIZE(candidate_response, 3) /= SIZE(reference_response, 3)) THEN
         CPABORT('ECD| Cartesian zero-order response check received incompatible tensor dimensions')
      END IF
      reference_norm = SQRT(SUM(ABS(reference_response)**2))
      candidate_norm = SQRT(SUM(ABS(candidate_response)**2))
      residual_norm = SQRT(SUM(ABS(candidate_response - reference_response)**2))
      relative_residual = 2.0_dp*residual_norm/MAX(SQRT(TINY(1.0_dp)), reference_norm + candidate_norm)
      IF (relative_residual > relative_tolerance) THEN
         CPABORT('ECD| Cartesian q directions produced inconsistent '//TRIM(representation_label)//' zero-order responses')
      END IF

   END SUBROUTINE ecd_periodic_response_zero_check

! **************************************************************************************************
!> \brief Assemble the folded local plus GTH diamagnetic kernel at implicit Gamma.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_diamagnetic_kernel(qs_env, context, gth_available, &
                                                              local_kernel_output, gth_kernel_output, &
                                                              diamagnetic_kernel)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      LOGICAL, INTENT(IN)                                :: gth_available
      REAL(KIND=dp), DIMENSION(ncart, ncart), INTENT(OUT) :: local_kernel_output, gth_kernel_output
      REAL(KIND=dp), DIMENSION(ncart, ncart), INTENT(OUT) :: diamagnetic_kernel

      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-4_dp

      INTEGER                                            :: ia, ib, iband, ispin, nmo
      REAL(KIND=dp)                                      :: gth_component
      REAL(KIND=dp), DIMENSION(ncart, ncart)             :: gth_kernel, local_kernel
      REAL(KIND=dp), DIMENSION(ncart)                    :: vector_potential
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: curvature_buffer
      TYPE(cp_fm_type)                                   :: curvature_mo
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: curvature, nl_mm, nl_mp, nl_pm, nl_pp, nl_zero

      NULLIFY (curvature, logger, nl_mm, nl_mp, nl_pm, nl_pp, nl_zero)
      local_kernel = 0.0_dp
      DO ispin = 1, context%nspin
         DO ia = 1, ncart
            local_kernel(ia, ia) = local_kernel(ia, ia) - SUM(context%spin(ispin)%occupations)
         END DO
      END DO
      gth_kernel = 0.0_dp

      IF (gth_available) THEN
         CALL ecd_periodic_ao_vector_create(curvature, context%matrix_s(1)%matrix, context%sab_orb, &
                                             'ECD implicit Gamma diamagnetic GTH curvature', dbcsr_type_symmetric)
         CALL ecd_periodic_gauge_nl_pair_create(nl_zero, context%matrix_s(1)%matrix, context%sab_orb, &
                                                'ECD implicit Gamma diamagnetic gauge zero')
         CALL ecd_periodic_gauge_nl_pair_create(nl_pp, context%matrix_s(1)%matrix, context%sab_orb, &
                                                'ECD implicit Gamma diamagnetic gauge pp')
         CALL ecd_periodic_gauge_nl_pair_create(nl_pm, context%matrix_s(1)%matrix, context%sab_orb, &
                                                'ECD implicit Gamma diamagnetic gauge pm')
         CALL ecd_periodic_gauge_nl_pair_create(nl_mp, context%matrix_s(1)%matrix, context%sab_orb, &
                                                'ECD implicit Gamma diamagnetic gauge mp')
         CALL ecd_periodic_gauge_nl_pair_create(nl_mm, context%matrix_s(1)%matrix, context%sab_orb, &
                                                'ECD implicit Gamma diamagnetic gauge mm')
         CALL ecd_periodic_gauge_nl_pair_reset(nl_zero)
         vector_potential = 0.0_dp
         CALL velocity_gauge_nl_term(qs_env, nl_zero, vector_potential)

         DO ia = 1, ncart
            DO ib = ia, ncart
               CALL ecd_periodic_gauge_nl_curvature_build(qs_env, nl_zero, nl_pp, nl_pm, nl_mp, nl_mm, &
                                                           fd_step, ia, ib, curvature(1))
               gth_component = 0.0_dp
               DO ispin = 1, context%nspin
                  nmo = context%spin(ispin)%nmo
                  CALL ecd_project_real_operator_to_mo(curvature(1), context%spin(ispin)%mo_coeff, &
                                                       context%spin(ispin)%nao, nmo, curvature_mo)
                  ALLOCATE (curvature_buffer(nmo, nmo))
                  CALL cp_fm_get_submatrix(curvature_mo, curvature_buffer)
                  DO iband = 1, nmo
                     gth_component = gth_component - context%spin(ispin)%occupations(iband)* &
                                     curvature_buffer(iband, iband)
                  END DO
                  DEALLOCATE (curvature_buffer)
                  CALL cp_fm_release(curvature_mo)
               END DO
               gth_kernel(ia, ib) = gth_component
               gth_kernel(ib, ia) = gth_component
            END DO
         END DO

         CALL dbcsr_deallocate_matrix_set(curvature)
         CALL dbcsr_deallocate_matrix_set(nl_mm)
         CALL dbcsr_deallocate_matrix_set(nl_mp)
         CALL dbcsr_deallocate_matrix_set(nl_pm)
         CALL dbcsr_deallocate_matrix_set(nl_pp)
         CALL dbcsr_deallocate_matrix_set(nl_zero)
      END IF

      diamagnetic_kernel = local_kernel + gth_kernel
      local_kernel_output = local_kernel
      gth_kernel_output = gth_kernel
      logger => cp_get_default_logger()
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         DO ia = 1, ncart
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G2 local diamagnetic kernel component', ia, local_kernel(ia, ia)
            DO ib = ia, ncart
               WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,A,I1,T55,ES14.6)') &
                  'ECD| Stage G2 GTH diamagnetic kernel component', ia, ' ', ib, gth_kernel(ia, ib)
               WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,A,I1,T55,ES14.6)') &
                  'ECD| Stage G2 total diamagnetic kernel component', ia, ' ', ib, diamagnetic_kernel(ia, ib)
            END DO
         END DO
      END IF

   END SUBROUTINE ecd_periodic_implicit_gamma_diamagnetic_kernel

! **************************************************************************************************
!> \brief Contract folded GTH two-field Taylor matrices with the occupied Gamma density.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_gth_diamagnetic_taylor(context, zero_kernel, first_kernel)
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET   :: context
      REAL(KIND=dp), DIMENSION(ncart, ncart), INTENT(OUT) :: zero_kernel
      COMPLEX(KIND=dp), DIMENSION(ncart, ncart), INTENT(OUT) :: first_kernel

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: operator_mo
      INTEGER                                            :: ia, ib, iband, ispin, nmo

      IF (.NOT. ASSOCIATED(context%gth_diamagnetic_taylor_real) .OR. &
          .NOT. ASSOCIATED(context%gth_diamagnetic_taylor_imag)) THEN
         CPABORT('ECD| GTH diamagnetic Taylor coefficients are unavailable')
      END IF
      IF (SIZE(context%gth_diamagnetic_taylor_real, 1) /= ncart .OR. &
          SIZE(context%gth_diamagnetic_taylor_real, 2) /= ncart .OR. &
          SIZE(context%gth_diamagnetic_taylor_real, 3) /= 2 .OR. &
          SIZE(context%gth_diamagnetic_taylor_imag, 1) /= ncart .OR. &
          SIZE(context%gth_diamagnetic_taylor_imag, 2) /= ncart .OR. &
          SIZE(context%gth_diamagnetic_taylor_imag, 3) /= 2) THEN
         CPABORT('ECD| GTH diamagnetic Taylor coefficient dimensions are inconsistent')
      END IF
      zero_kernel = 0.0_dp
      first_kernel = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO ia = 1, ncart
         DO ib = 1, ncart
            DO ispin = 1, context%nspin
               nmo = context%spin(ispin)%nmo
               CALL ecd_project_complex_operator_to_mo(context%gth_diamagnetic_taylor_real(ia, ib, 1), &
                                                       context%gth_diamagnetic_taylor_imag(ia, ib, 1), &
                                                       context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                       operator_mo)
               DO iband = 1, nmo
                  zero_kernel(ia, ib) = zero_kernel(ia, ib) - context%spin(ispin)%occupations(iband)* &
                                         REAL(operator_mo(iband, iband), KIND=dp)
               END DO
               DEALLOCATE (operator_mo)
               CALL ecd_project_complex_operator_to_mo(context%gth_diamagnetic_taylor_real(ia, ib, 2), &
                                                       context%gth_diamagnetic_taylor_imag(ia, ib, 2), &
                                                       context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                       operator_mo)
               DO iband = 1, nmo
                  first_kernel(ia, ib) = first_kernel(ia, ib) - context%spin(ispin)%occupations(iband)* &
                                          operator_mo(iband, iband)
               END DO
               DEALLOCATE (operator_mo)
            END DO
         END DO
      END DO

   END SUBROUTINE ecd_periodic_implicit_gamma_gth_diamagnetic_taylor

! **************************************************************************************************
!> \brief Write complete Cartesian first-order periodic IPA response tensors.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_print_cartesian_response(ecd_section, energy_range, velocity_zero, velocity_first, &
                                                    context, gth_available, &
                                                    pseudopotential_commutator, fwhm, selection)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: velocity_zero
      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN) :: velocity_first
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      LOGICAL, INTENT(IN)                                :: gth_available
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), INTENT(IN)                          :: fwhm
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: selection

      CHARACTER(LEN=default_path_length)                 :: response_filename
      CHARACTER(LEN=64)                                  :: column_label
      CHARACTER(LEN=1)                                   :: q_label
      COMPLEX(KIND=dp)                                   :: response_value
      INTEGER                                            :: ia, ib, ipoint, iq, response_unit
      REAL(KIND=dp)                                      :: photon_energy
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tensor_print_key

      IF (SIZE(velocity_zero, 1) /= ncart .OR. SIZE(velocity_zero, 2) /= ncart .OR. &
          SIZE(velocity_zero, 3) < 2 .OR. &
          SIZE(velocity_first, 1) /= ncart .OR. SIZE(velocity_first, 2) /= ncart .OR. &
          SIZE(velocity_first, 3) /= ncart .OR. SIZE(velocity_first, 4) /= SIZE(velocity_zero, 3)) THEN
         CPABORT('ECD| Cartesian periodic response output received incompatible tensor dimensions')
      END IF

      NULLIFY (logger, tensor_print_key)
      logger => cp_get_default_logger()
      tensor_print_key => section_vals_get_subs_vals(ecd_section, 'PRINT%RESPONSE_TENSOR')
      response_filename = cp_print_key_generate_filename(logger, tensor_print_key, &
                                                         extension='.ecd_response', my_local=.FALSE.)
      response_unit = cp_print_key_unit_nr(logger, ecd_section, 'PRINT%RESPONSE_TENSOR', &
                                            extension='.ecd_response', file_position='REWIND', &
                                            file_action='WRITE', file_form='FORMATTED')
      IF (response_unit <= 0) RETURN

      WRITE (response_unit, '(A)') '# CP2K periodic ECD Cartesian dielectric response'
      WRITE (response_unit, '(A)') '# Format_version 9'
      WRITE (response_unit, '(A)') '# Gamma_layout IMPLICIT'
      WRITE (response_unit, '(A)') '# Electric_representation COMPLETE_CURRENT_VELOCITY'
      WRITE (response_unit, '(A,1X,L1)') '# GTH_current_available', gth_available
      WRITE (response_unit, '(A,1X,ES18.10E3)') '# FWHM_hartree', fwhm
      WRITE (response_unit, '(A,3(1X,ES18.10E3))') '# Reference_point_bohr', context%ref_point
      WRITE (response_unit, '(A,2(1X,I0))') '# Transition_band_range', selection%band_range
      WRITE (response_unit, '(A,1X,I0)') '# Requested_occupied_empty_pairs', selection%requested_count
      WRITE (response_unit, '(A,1X,I0)') '# Candidate_occupied_empty_pairs', selection%candidate_count
      WRITE (response_unit, '(A,1X,I0)') '# Selected_occupied_empty_pairs', selection%selected_count
      WRITE (response_unit, '(A,1X,ES18.10E3)') '# Transition_cutoff_gap_hartree', selection%cutoff_gap
      WRITE (response_unit, '(A,2(1X,ES18.10E3))') '# Energy_range_hartree', energy_range
      WRITE (response_unit, '(A,1X,I0)') '# NPOINTS', SIZE(velocity_zero, 3)
      IF (pseudopotential_commutator) THEN
         WRITE (response_unit, '(A)') '# Velocity epsilon uses the complete current plus diamagnetic kernel.'
      ELSE
         WRITE (response_unit, '(A)') '# Velocity epsilon uses the local-current ablation plus local diamagnetic kernel.'
      END IF
      WRITE (response_unit, '(A)') '# epsilon1_abq stores the q_x, q_y, q_z Cartesian derivative components.'
      IF (pseudopotential_commutator) THEN
         WRITE (response_unit, '(A)') '# Velocity epsilon1 includes the GTH finite-q two-field diamagnetic coefficient.'
      ELSE
         WRITE (response_unit, '(A)') &
            '# Velocity epsilon1 uses the local-current ablation and excludes the GTH two-field coefficient.'
      END IF
      WRITE (response_unit, '(A)', ADVANCE='NO') '# Columns: Energy_hartree Energy_eV'
      DO ia = 1, ncart
         DO ib = 1, ncart
            WRITE (column_label, '(A,I0,A,I0,A)') 'velocity_epsilon0_', ia, '_', ib, '_re'
            WRITE (response_unit, '(1X,A)', ADVANCE='NO') TRIM(column_label)
            WRITE (column_label, '(A,I0,A,I0,A)') 'velocity_epsilon0_', ia, '_', ib, '_im'
            WRITE (response_unit, '(1X,A)', ADVANCE='NO') TRIM(column_label)
         END DO
      END DO
      DO ia = 1, ncart
         DO ib = 1, ncart
            DO iq = 1, ncart
               q_label = ACHAR(IACHAR('w') + iq)
               WRITE (column_label, '(A,I0,A,I0,A,A,A)') 'velocity_epsilon1_', ia, '_', ib, '_q', q_label, '_re'
               WRITE (response_unit, '(1X,A)', ADVANCE='NO') TRIM(column_label)
               WRITE (column_label, '(A,I0,A,I0,A,A,A)') 'velocity_epsilon1_', ia, '_', ib, '_q', q_label, '_im'
               WRITE (response_unit, '(1X,A)', ADVANCE='NO') TRIM(column_label)
            END DO
         END DO
      END DO
      WRITE (response_unit, '()')

      DO ipoint = 1, SIZE(velocity_zero, 3)
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                          (energy_range(2) - energy_range(1))/REAL(SIZE(velocity_zero, 3) - 1, KIND=dp)
         WRITE (response_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') photon_energy, photon_energy*evolt
         DO ia = 1, ncart
            DO ib = 1, ncart
               response_value = velocity_zero(ia, ib, ipoint)
               IF (ia == ib) response_value = response_value + CMPLX(1.0_dp, 0.0_dp, KIND=dp)
               WRITE (response_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') &
                  REAL(response_value, KIND=dp), AIMAG(response_value)
            END DO
         END DO
         DO ia = 1, ncart
            DO ib = 1, ncart
               DO iq = 1, ncart
                  response_value = velocity_first(ia, ib, iq, ipoint)
                  WRITE (response_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') &
                     REAL(response_value, KIND=dp), AIMAG(response_value)
               END DO
            END DO
         END DO
         WRITE (response_unit, '()')
      END DO
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,1X,A)') &
            'ECD| Periodic Cartesian response tensor file', TRIM(response_filename)
      END IF
      CALL cp_print_key_finished_output(response_unit, logger, ecd_section, 'PRINT%RESPONSE_TENSOR')

   END SUBROUTINE ecd_periodic_print_cartesian_response

! **************************************************************************************************
!> \brief Write additive local, cross, nonlocal, and complete velocity susceptibilities.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_print_cartesian_velocity_decomposition(ecd_section, energy_range, local_zero, local_first, &
                                                                   cross_zero, cross_first, nonlocal_zero, nonlocal_first, &
                                                                   complete_zero, complete_first, context, gth_available, &
                                                                   fwhm, selection)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: complete_zero, cross_zero, local_zero, nonlocal_zero
      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN) :: complete_first, cross_first, local_first, nonlocal_first
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      LOGICAL, INTENT(IN)                                :: gth_available
      REAL(KIND=dp), INTENT(IN)                          :: fwhm
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: selection

      CHARACTER(LEN=default_path_length)                 :: decomposition_filename
      COMPLEX(KIND=dp)                                   :: response_value
      INTEGER                                            :: ia, ib, icomponent, idir, ipoint, decomposition_unit
      REAL(KIND=dp)                                      :: closure_first, closure_zero, photon_energy
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tensor_print_key

      IF (SIZE(complete_zero, 1) /= ncart .OR. SIZE(complete_zero, 2) /= ncart .OR. &
          SIZE(complete_zero, 3) < 2 .OR. ANY(SHAPE(local_zero) /= SHAPE(complete_zero)) .OR. &
          ANY(SHAPE(cross_zero) /= SHAPE(complete_zero)) .OR. ANY(SHAPE(nonlocal_zero) /= SHAPE(complete_zero)) .OR. &
          SIZE(complete_first, 1) /= ncart .OR. SIZE(complete_first, 2) /= ncart .OR. &
          SIZE(complete_first, 3) /= ncart .OR. SIZE(complete_first, 4) /= SIZE(complete_zero, 3) .OR. &
          ANY(SHAPE(local_first) /= SHAPE(complete_first)) .OR. &
          ANY(SHAPE(cross_first) /= SHAPE(complete_first)) .OR. &
          ANY(SHAPE(nonlocal_first) /= SHAPE(complete_first))) THEN
         CPABORT('ECD| Velocity decomposition output received incompatible response tensors')
      END IF

      NULLIFY (logger, tensor_print_key)
      logger => cp_get_default_logger()
      tensor_print_key => section_vals_get_subs_vals(ecd_section, 'PRINT%RESPONSE_TENSOR')
      decomposition_filename = cp_print_key_generate_filename(logger, tensor_print_key, &
                                                               extension='.ecd_velocity_decomposition', my_local=.FALSE.)
      decomposition_unit = cp_print_key_unit_nr(logger, ecd_section, 'PRINT%RESPONSE_TENSOR', &
                                                 extension='.ecd_velocity_decomposition', file_position='REWIND', &
                                                 file_action='WRITE', file_form='FORMATTED')
      IF (decomposition_unit <= 0) RETURN

      WRITE (decomposition_unit, '(A)') '# CP2K periodic ECD velocity susceptibility decomposition'
      WRITE (decomposition_unit, '(A)') '# Format_version 1'
      WRITE (decomposition_unit, '(A)') '# Gamma_layout IMPLICIT'
      WRITE (decomposition_unit, '(A,1X,L1)') '# GTH_current_available', gth_available
      WRITE (decomposition_unit, '(A,1X,ES18.10E3)') '# FWHM_hartree', fwhm
      WRITE (decomposition_unit, '(A,3(1X,ES18.10E3))') '# Reference_point_bohr', context%ref_point
      WRITE (decomposition_unit, '(A,1X,I0)') '# Selected_occupied_empty_pairs', selection%selected_count
      WRITE (decomposition_unit, '(A,2(1X,ES18.10E3))') '# Energy_range_hartree', energy_range
      WRITE (decomposition_unit, '(A)') '# Components LOCAL CROSS NONLOCAL COMPLETE'
      WRITE (decomposition_unit, '(A)') '# LOCAL includes the local diamagnetic kernel.'
      WRITE (decomposition_unit, '(A)') '# NONLOCAL includes the GTH diamagnetic kernel.'
      WRITE (decomposition_unit, '(A)') '# Per component: chi0_ab_re chi0_ab_im for a,b=X,Y,Z, then chi1_abq_re chi1_abq_im.'
      WRITE (decomposition_unit, '(A)') '# Final columns: zero_order_closure_norm first_order_closure_norm.'

      DO ipoint = 1, SIZE(complete_zero, 3)
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                         (energy_range(2) - energy_range(1))/REAL(SIZE(complete_zero, 3) - 1, KIND=dp)
         closure_zero = SQRT(SUM(ABS(complete_zero(:, :, ipoint) - local_zero(:, :, ipoint) - &
                                     cross_zero(:, :, ipoint) - nonlocal_zero(:, :, ipoint))**2))
         closure_first = SQRT(SUM(ABS(complete_first(:, :, :, ipoint) - local_first(:, :, :, ipoint) - &
                                      cross_first(:, :, :, ipoint) - nonlocal_first(:, :, :, ipoint))**2))
         WRITE (decomposition_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') photon_energy, photon_energy*evolt
         DO icomponent = 1, 4
            DO ia = 1, ncart
               DO ib = 1, ncart
                  SELECT CASE (icomponent)
                  CASE (1)
                     response_value = local_zero(ia, ib, ipoint)
                  CASE (2)
                     response_value = cross_zero(ia, ib, ipoint)
                  CASE (3)
                     response_value = nonlocal_zero(ia, ib, ipoint)
                  CASE DEFAULT
                     response_value = complete_zero(ia, ib, ipoint)
                  END SELECT
                  WRITE (decomposition_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') &
                     REAL(response_value, KIND=dp), AIMAG(response_value)
               END DO
            END DO
            DO ia = 1, ncart
               DO ib = 1, ncart
                  DO idir = 1, ncart
                     SELECT CASE (icomponent)
                     CASE (1)
                        response_value = local_first(ia, ib, idir, ipoint)
                     CASE (2)
                        response_value = cross_first(ia, ib, idir, ipoint)
                     CASE (3)
                        response_value = nonlocal_first(ia, ib, idir, ipoint)
                     CASE DEFAULT
                        response_value = complete_first(ia, ib, idir, ipoint)
                     END SELECT
                     WRITE (decomposition_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') &
                        REAL(response_value, KIND=dp), AIMAG(response_value)
                  END DO
               END DO
            END DO
         END DO
         WRITE (decomposition_unit, '(2(1X,ES18.10E3))') closure_zero, closure_first
      END DO
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,1X,A)') &
            'ECD| Periodic velocity decomposition file', TRIM(decomposition_filename)
      END IF
      CALL cp_print_key_finished_output(decomposition_unit, logger, ecd_section, 'PRINT%RESPONSE_TENSOR')

   END SUBROUTINE ecd_periodic_print_cartesian_velocity_decomposition

! **************************************************************************************************
!> \brief Write complete periodic IPA transition coefficients for both electric representations.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_print_cartesian_transition_matrix_elements(ecd_section, context, transition_selection, &
                                                                       gth_requested, gap_first, velocity_transition, &
                                                                       gth_transition, length_transition, &
                                                                       static_magnetic_transition)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: transition_selection
      LOGICAL, INTENT(IN)                                :: gth_requested
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)          :: gap_first
      COMPLEX(KIND=dp), DIMENSION(:, 0:, :, :), INTENT(IN) :: gth_transition, velocity_transition
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: length_transition
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: static_magnetic_transition

      CHARACTER(LEN=default_path_length)                 :: matrix_filename
      INTEGER                                            :: idir, iocc, ispin, iq, itransition, ivirt, matrix_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: matrix_print_key

      IF (SIZE(velocity_transition, 1) /= ncart .OR. SIZE(gth_transition, 1) /= ncart .OR. &
          SIZE(length_transition, 1) /= ncart .OR. SIZE(static_magnetic_transition, 1) /= ncart .OR. &
          LBOUND(velocity_transition, 2) /= 0 .OR. &
          LBOUND(gth_transition, 2) /= 0 .OR. &
          UBOUND(velocity_transition, 2) < 1 .OR. UBOUND(gth_transition, 2) < 1 .OR. &
          SIZE(velocity_transition, 3) /= transition_selection%selected_count .OR. &
          SIZE(gth_transition, 3) /= transition_selection%selected_count .OR. &
          SIZE(length_transition, 2) /= transition_selection%selected_count .OR. &
          SIZE(static_magnetic_transition, 2) /= transition_selection%selected_count .OR. &
          SIZE(velocity_transition, 4) /= ncart .OR. SIZE(gth_transition, 4) /= ncart .OR. &
          SIZE(static_magnetic_transition, 3) /= ncart .OR. &
          SIZE(gap_first, 1) /= transition_selection%selected_count .OR. SIZE(gap_first, 2) /= ncart) THEN
         CPABORT('ECD| Periodic transition-matrix output received inconsistent transition coefficients')
      END IF
      IF (.NOT. ALLOCATED(context%spin)) THEN
         CPABORT('ECD| Periodic transition-matrix output requires retained Gamma spin data')
      END IF

      NULLIFY (logger, matrix_print_key)
      logger => cp_get_default_logger()
      matrix_print_key => section_vals_get_subs_vals(ecd_section, 'PRINT%MATRIX_ELEMENTS')
      matrix_filename = cp_print_key_generate_filename(logger, matrix_print_key, &
                                                       extension='.ecd_matrix', my_local=.FALSE.)
      matrix_unit = cp_print_key_unit_nr(logger, ecd_section, 'PRINT%MATRIX_ELEMENTS', &
                                         extension='.ecd_matrix', file_position='REWIND', &
                                         file_action='WRITE', file_form='FORMATTED')
      IF (matrix_unit <= 0) RETURN

      WRITE (matrix_unit, '(A)') '# CP2K periodic ECD IPA transition matrix elements'
      WRITE (matrix_unit, '(A)') '# Format_version 7'
      WRITE (matrix_unit, '(A)') '# Gamma_layout IMPLICIT'
      WRITE (matrix_unit, '(A,1X,L1)') '# Pseudopotential_commutator_selects_complete', gth_requested
      WRITE (matrix_unit, '(A)') '# GTH current columns are retained for both selector values.'
      WRITE (matrix_unit, '(A,3(1X,ES18.10E3))') '# Reference_point_bohr', context%ref_point
      WRITE (matrix_unit, '(A)') '# Q_derivative_components X Y Z'
      WRITE (matrix_unit, '(A,2(1X,I0))') '# Transition_band_range', transition_selection%band_range
      WRITE (matrix_unit, '(A,1X,I0)') '# Requested_occupied_empty_pairs', transition_selection%requested_count
      WRITE (matrix_unit, '(A,1X,I0)') '# Candidate_occupied_empty_pairs', transition_selection%candidate_count
      WRITE (matrix_unit, '(A,1X,I0)') '# Selected_occupied_empty_pairs', transition_selection%selected_count
      WRITE (matrix_unit, '(A,1X,ES18.10E3)') '# Transition_cutoff_gap_hartree', transition_selection%cutoff_gap
      WRITE (matrix_unit, '(A)') '# Current coefficients use normalized q-Taylor orders 0 and 1.'
      WRITE (matrix_unit, '(A)') '# Per Cartesian component: complete_current0_re complete_current0_im '// &
                                  'complete_current1_re complete_current1_im local_current0_re local_current0_im '// &
                                  'local_current1_re local_current1_im gth_current0_re gth_current0_im '// &
                                  'gth_current1_re gth_current1_im direct_position0_re direct_position0_im'
      WRITE (matrix_unit, '(A)') '# Current coefficients define the periodic finite-q observable and its first derivative.'
      WRITE (matrix_unit, '(A)') '# Direct position0 is retained for isolated-molecule-limit diagnostics only.'
      WRITE (matrix_unit, '(A)') '# Static magnetic diagnostic is not used for finite-q ECD assembly.'
      WRITE (matrix_unit, '(A)') '# It supports molecular-limit rotatory-strength validation only.'
      WRITE (matrix_unit, '(A)') '# Columns: Transition_index Q_derivative_component Spin Occupied_band '// &
                                  'Virtual_band Occupation_difference '// &
                                  'Gap0_hartree Gap1_hartree followed by X, Y, and Z component groups and '// &
                                  'static_magnetic_x_re static_magnetic_x_im static_magnetic_y_re '// &
                                  'static_magnetic_y_im static_magnetic_z_re static_magnetic_z_im.'

      DO itransition = 1, transition_selection%selected_count
         ispin = transition_selection%ispin(itransition)
         iocc = transition_selection%iocc(itransition)
         ivirt = transition_selection%ivirt(itransition)
         IF (ispin < 1 .OR. ispin > context%nspin) THEN
            CPABORT('ECD| Periodic transition-matrix output found an invalid spin index')
         END IF
         IF (iocc < 1 .OR. ivirt < 1 .OR. iocc > context%spin(ispin)%nmo .OR. &
             ivirt > context%spin(ispin)%nmo .OR. &
             .NOT. ASSOCIATED(context%spin(ispin)%eigenvalue_first)) THEN
            CPABORT('ECD| Periodic transition-matrix output found incomplete first-order state data')
         END IF
         DO iq = 1, ncart
            WRITE (matrix_unit, '(5(I0,1X),3(ES18.10E3,1X))', ADVANCE='NO') itransition, iq, ispin, iocc, ivirt, &
               transition_selection%occupation_difference(itransition), transition_selection%gap(itransition), &
               gap_first(itransition, iq)
            DO idir = 1, ncart
               WRITE (matrix_unit, '(14(ES18.10E3,1X))', ADVANCE='NO') &
                  REAL(velocity_transition(idir, 0, itransition, iq), KIND=dp), &
                  AIMAG(velocity_transition(idir, 0, itransition, iq)), &
                  REAL(velocity_transition(idir, 1, itransition, iq), KIND=dp), &
                  AIMAG(velocity_transition(idir, 1, itransition, iq)), &
                  REAL(velocity_transition(idir, 0, itransition, iq) - &
                       gth_transition(idir, 0, itransition, iq), KIND=dp), &
                  AIMAG(velocity_transition(idir, 0, itransition, iq) - &
                        gth_transition(idir, 0, itransition, iq)), &
                  REAL(velocity_transition(idir, 1, itransition, iq) - &
                       gth_transition(idir, 1, itransition, iq), KIND=dp), &
                  AIMAG(velocity_transition(idir, 1, itransition, iq) - &
                        gth_transition(idir, 1, itransition, iq)), &
                  REAL(gth_transition(idir, 0, itransition, iq), KIND=dp), &
                  AIMAG(gth_transition(idir, 0, itransition, iq)), &
                  REAL(gth_transition(idir, 1, itransition, iq), KIND=dp), &
                  AIMAG(gth_transition(idir, 1, itransition, iq)), &
                  REAL(length_transition(idir, itransition), KIND=dp), &
                  AIMAG(length_transition(idir, itransition))
            END DO
            DO idir = 1, ncart
               WRITE (matrix_unit, '(2(ES18.10E3,1X))', ADVANCE='NO') &
                  REAL(static_magnetic_transition(idir, itransition, iq), KIND=dp), &
                  AIMAG(static_magnetic_transition(idir, itransition, iq))
            END DO
            WRITE (matrix_unit, '()')
         END DO
      END DO

      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,1X,A)') &
            'ECD| Periodic transition matrix file', TRIM(matrix_filename)
      END IF
      CALL cp_print_key_finished_output(matrix_unit, logger, ecd_section, 'PRINT%MATRIX_ELEMENTS')

   END SUBROUTINE ecd_periodic_print_cartesian_transition_matrix_elements

! **************************************************************************************************
!> \brief Allocate a Cartesian vector of real AO matrices on an explicit-Gamma neighbour list.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_ao_vector_create(matrices, template, sab_nl, name, matrix_type)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices
      TYPE(dbcsr_type), POINTER                          :: template
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL             :: matrix_type

      INTEGER                                            :: idir
      CHARACTER(LEN=1)                                   :: my_matrix_type

      my_matrix_type = dbcsr_type_no_symmetry
      IF (PRESENT(matrix_type)) my_matrix_type = matrix_type
      CALL dbcsr_allocate_matrix_set(matrices, ncart)
      DO idir = 1, ncart
         CALL dbcsr_init_p(matrices(idir)%matrix)
         CALL dbcsr_create(matrices(idir)%matrix, template=template, &
                           matrix_type=my_matrix_type, name=name)
         CALL cp_dbcsr_alloc_block_from_nbl(matrices(idir)%matrix, sab_nl)
         CALL dbcsr_set(matrices(idir)%matrix, 0.0_dp)
      END DO

   END SUBROUTINE ecd_periodic_ao_vector_create

! **************************************************************************************************
!> \brief Allocate a set of real AO matrices on an explicit-Gamma neighbour list.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_ao_set_create(matrices, nmatrices, template, sab_nl, name, matrix_type)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices
      INTEGER, INTENT(IN)                                :: nmatrices
      TYPE(dbcsr_type), POINTER                          :: template
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      CHARACTER(LEN=1), INTENT(IN)                       :: matrix_type

      INTEGER                                            :: imatrix

      CALL dbcsr_allocate_matrix_set(matrices, nmatrices)
      DO imatrix = 1, nmatrices
         CALL dbcsr_init_p(matrices(imatrix)%matrix)
         CALL dbcsr_create(matrices(imatrix)%matrix, template=template, &
                           matrix_type=matrix_type, name=name)
         CALL cp_dbcsr_alloc_block_from_nbl(matrices(imatrix)%matrix, sab_nl)
         CALL dbcsr_set(matrices(imatrix)%matrix, 0.0_dp)
      END DO

   END SUBROUTINE ecd_periodic_ao_set_create

! **************************************************************************************************
!> \brief Build complete local and GTH current Taylor coefficients in the explicit-Gamma AO basis.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_current_taylor_build(current_real, current_imag, gth_real, gth_imag, &
                                                local_momentum_rs, local_moments_der_rs, local_moments_rs, &
                                                gth_current_path_rs, include_gth_path, matrix_s, xkp, &
                                                field_cell_to_index, local_cell_to_index, gth_cell_to_index, sab_nl, &
                                                hmat, direction, &
                                                max_order)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: current_real, current_imag, gth_real, gth_imag
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: local_momentum_rs, local_moments_rs, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: gth_current_path_rs, local_moments_der_rs
      LOGICAL, INTENT(IN)                                 :: include_gth_path
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)             :: xkp, direction
      INTEGER, DIMENSION(:, :, :), POINTER                :: field_cell_to_index, local_cell_to_index, &
                                                             gth_cell_to_index, source_cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                          :: sab_nl
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)          :: hmat
      INTEGER, INTENT(IN)                                 :: max_order

      COMPLEX(KIND=dp)                                    :: scale
      INTEGER                                             :: idir, ix, iy, iz, lattice_order, &
                                                             moment_index, operator_order, total_order
      REAL(KIND=dp)                                       :: directional_weight
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: source_rs
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: work_imag, work_real

      IF (max_order < 0) THEN
         CPABORT("ECD| Current Taylor construction requires a non-negative order")
      END IF
      IF (SIZE(local_momentum_rs, 1) /= ncart .OR. SIZE(matrix_s, 1) /= 1) THEN
         CPABORT("ECD| Current Taylor construction received incompatible local AO storage")
      END IF
      IF (SIZE(local_moments_der_rs, 2) /= ncart .OR. SIZE(local_moments_rs, 2) /= SIZE(local_momentum_rs, 2)) THEN
         CPABORT("ECD| Current Taylor construction received incompatible local multipole storage")
      END IF
      IF (include_gth_path) THEN
         IF (SIZE(gth_current_path_rs, 1) /= ncart .OR. SIZE(gth_current_path_rs, 2) < max_order + 1) THEN
            CPABORT("ECD| Current Taylor construction received incompatible GTH path storage")
         END IF
      END IF

      CALL dbcsr_allocate_matrix_set(current_real, ncart, max_order + 1)
      CALL dbcsr_allocate_matrix_set(current_imag, ncart, max_order + 1)
      CALL dbcsr_allocate_matrix_set(gth_real, ncart, max_order + 1)
      CALL dbcsr_allocate_matrix_set(gth_imag, ncart, max_order + 1)
      CALL ecd_periodic_ao_set_create(work_real, 1, matrix_s(1, 1)%matrix, sab_nl, &
                                      'ECD CURRENT TAYLOR WORK REAL', dbcsr_type_no_symmetry)
      CALL ecd_periodic_ao_set_create(work_imag, 1, matrix_s(1, 1)%matrix, sab_nl, &
                                      'ECD CURRENT TAYLOR WORK IMAG', dbcsr_type_no_symmetry)
      DO total_order = 0, max_order
         DO idir = 1, ncart
            CALL dbcsr_init_p(current_real(idir, total_order + 1)%matrix)
            CALL dbcsr_create(current_real(idir, total_order + 1)%matrix, template=matrix_s(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, name="ECD CURRENT TAYLOR REAL")
            CALL cp_dbcsr_alloc_block_from_nbl(current_real(idir, total_order + 1)%matrix, sab_nl)
            CALL dbcsr_set(current_real(idir, total_order + 1)%matrix, 0.0_dp)
            CALL dbcsr_init_p(current_imag(idir, total_order + 1)%matrix)
            CALL dbcsr_create(current_imag(idir, total_order + 1)%matrix, template=matrix_s(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, name="ECD CURRENT TAYLOR IMAG")
            CALL cp_dbcsr_alloc_block_from_nbl(current_imag(idir, total_order + 1)%matrix, sab_nl)
            CALL dbcsr_set(current_imag(idir, total_order + 1)%matrix, 0.0_dp)
            CALL dbcsr_init_p(gth_real(idir, total_order + 1)%matrix)
            CALL dbcsr_create(gth_real(idir, total_order + 1)%matrix, template=matrix_s(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, name="ECD GTH CURRENT TAYLOR REAL")
            CALL cp_dbcsr_alloc_block_from_nbl(gth_real(idir, total_order + 1)%matrix, sab_nl)
            CALL dbcsr_set(gth_real(idir, total_order + 1)%matrix, 0.0_dp)
            CALL dbcsr_init_p(gth_imag(idir, total_order + 1)%matrix)
            CALL dbcsr_create(gth_imag(idir, total_order + 1)%matrix, template=matrix_s(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, name="ECD GTH CURRENT TAYLOR IMAG")
            CALL cp_dbcsr_alloc_block_from_nbl(gth_imag(idir, total_order + 1)%matrix, sab_nl)
            CALL dbcsr_set(gth_imag(idir, total_order + 1)%matrix, 0.0_dp)
         END DO
      END DO

      DO total_order = 0, max_order
         DO idir = 1, ncart
            DO operator_order = 0, total_order
               lattice_order = total_order - operator_order
               IF (operator_order == 0) THEN
                  CALL rskp_directional_taylor_transform(work_real(1)%matrix, work_imag(1)%matrix, &
                                                         local_momentum_rs, idir, xkp, local_cell_to_index, sab_nl, &
                                                         hmat, direction, lattice_order, shift_scale=-0.5_dp)
                  CALL ecd_periodic_ao_add_complex_scaled(current_real(idir, total_order + 1), &
                                                          current_imag(idir, total_order + 1), work_real(1), &
                                                          work_imag(1), CMPLX(1.0_dp, 0.0_dp, KIND=dp))
               ELSE
                  DO ix = 0, operator_order
                     DO iy = 0, operator_order - ix
                        iz = operator_order - ix - iy
                        ! moments_der_rs omits the scalar Cartesian moment at coset index one.
                        moment_index = coset(ix, iy, iz) - 1
                        IF (moment_index < 1 .OR. moment_index > SIZE(local_moments_der_rs, 1)) THEN
                           CPABORT("ECD| Current Taylor derivative-moment index is outside local storage")
                        END IF
                        source_rs => local_moments_der_rs(moment_index:moment_index, idir, :)
                        source_cell_to_index => local_cell_to_index
                        CALL rskp_directional_taylor_transform(work_real(1)%matrix, work_imag(1)%matrix, &
                                                               source_rs, 1, xkp, source_cell_to_index, sab_nl, hmat, &
                                                               direction, lattice_order, shift_scale=-0.5_dp)
                        directional_weight = ecd_periodic_factorial(operator_order)/ &
                           (ecd_periodic_factorial(ix)*ecd_periodic_factorial(iy)*ecd_periodic_factorial(iz))
                        IF (ix > 0) directional_weight = directional_weight*direction(1)**ix
                        IF (iy > 0) directional_weight = directional_weight*direction(2)**iy
                        IF (iz > 0) directional_weight = directional_weight*direction(3)**iz
                        scale = ecd_periodic_i_power(operator_order)*directional_weight/ &
                                ecd_periodic_factorial(operator_order)
                        CALL ecd_periodic_ao_add_complex_scaled(current_real(idir, total_order + 1), &
                                                                current_imag(idir, total_order + 1), work_real(1), &
                                                                work_imag(1), scale)
                        NULLIFY (source_rs)
                     END DO
                  END DO
                  DO ix = 0, operator_order - 1
                     DO iy = 0, operator_order - 1 - ix
                        iz = operator_order - 1 - ix - iy
                        directional_weight = ecd_periodic_factorial(operator_order - 1)/ &
                           (ecd_periodic_factorial(ix)*ecd_periodic_factorial(iy)*ecd_periodic_factorial(iz))
                        IF (ix > 0) directional_weight = directional_weight*direction(1)**ix
                        IF (iy > 0) directional_weight = directional_weight*direction(2)**iy
                        IF (iz > 0) directional_weight = directional_weight*direction(3)**iz
                        IF (operator_order == 1) THEN
                           source_rs => matrix_s
                           source_cell_to_index => field_cell_to_index
                        ELSE
                           ! moments_rs uses the same scalar-free Cartesian-moment convention.
                           moment_index = coset(ix, iy, iz) - 1
                           IF (moment_index < 1 .OR. moment_index > SIZE(local_moments_rs, 1)) THEN
                              CPABORT("ECD| Current Taylor coordinate-moment index is outside local storage")
                           END IF
                           source_rs => local_moments_rs(moment_index:moment_index, :)
                           source_cell_to_index => local_cell_to_index
                        END IF
                        CALL rskp_directional_taylor_transform(work_real(1)%matrix, work_imag(1)%matrix, &
                                                               source_rs, 1, xkp, source_cell_to_index, sab_nl, hmat, &
                                                               direction, lattice_order, shift_scale=-0.5_dp)
                        scale = ecd_periodic_i_power(operator_order)*direction(idir)*directional_weight/ &
                                (2.0_dp*ecd_periodic_factorial(operator_order - 1))
                        CALL ecd_periodic_ao_add_complex_scaled(current_real(idir, total_order + 1), &
                                                                current_imag(idir, total_order + 1), work_real(1), &
                                                                work_imag(1), scale)
                        NULLIFY (source_rs)
                        NULLIFY (source_cell_to_index)
                     END DO
                  END DO
               END IF
               IF (include_gth_path) THEN
                  source_rs => gth_current_path_rs(idir:idir, operator_order + 1, :)
                  source_cell_to_index => gth_cell_to_index
                  CALL rskp_directional_taylor_transform(work_real(1)%matrix, work_imag(1)%matrix, source_rs, &
                                                         1, xkp, source_cell_to_index, sab_nl, hmat, direction, &
                                                         lattice_order, shift_scale=-0.5_dp)
                  ! Match the finite-q Fourier convention used by the local current vertex.
                  scale = CONJG(ecd_periodic_i_power(operator_order))
                  CALL ecd_periodic_ao_add_complex_scaled(current_real(idir, total_order + 1), &
                                                          current_imag(idir, total_order + 1), work_real(1), &
                                                          work_imag(1), scale)
                  CALL ecd_periodic_ao_add_complex_scaled(gth_real(idir, total_order + 1), &
                                                          gth_imag(idir, total_order + 1), work_real(1), &
                                                          work_imag(1), scale)
                  NULLIFY (source_rs)
                  NULLIFY (source_cell_to_index)
               END IF
            END DO
         END DO
      END DO

      CALL dbcsr_deallocate_matrix_set(work_imag)
      CALL dbcsr_deallocate_matrix_set(work_real)

   END SUBROUTINE ecd_periodic_current_taylor_build

! **************************************************************************************************
!> \brief Fold Pickard--Mauri two-field path coefficients into the implicit-Gamma AO basis.
!>
!> The raw first path coefficient is real.  Its physical Hamiltonian
!> coefficient carries the explicit imaginary Taylor factor.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_gth_diamagnetic_taylor_build(diamagnetic_real, diamagnetic_imag, path_rs, xkp_zero, &
                                                         cell_to_index, sab_nl)
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: diamagnetic_imag, diamagnetic_real
      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: path_rs
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_zero
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_nl

      INTEGER                                            :: ia, ib, iorder
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: source_rs

      IF (SIZE(path_rs, 1) /= ncart .OR. SIZE(path_rs, 2) /= ncart .OR. SIZE(path_rs, 3) /= 2) THEN
         CPABORT('ECD| GTH diamagnetic Taylor construction received incompatible path storage')
      END IF
      CALL dbcsr_allocate_matrix_set(diamagnetic_real, ncart, ncart, 2)
      CALL dbcsr_allocate_matrix_set(diamagnetic_imag, ncart, ncart, 2)
      DO iorder = 1, 2
         DO ia = 1, ncart
            DO ib = 1, ncart
               CALL dbcsr_init_p(diamagnetic_real(ia, ib, iorder)%matrix)
               CALL dbcsr_create(diamagnetic_real(ia, ib, iorder)%matrix, &
                                 template=path_rs(ia, ib, iorder, 1)%matrix, matrix_type=dbcsr_type_no_symmetry, &
                                 name='ECD GTH DIAMAGNETIC TAYLOR REAL')
               CALL cp_dbcsr_alloc_block_from_nbl(diamagnetic_real(ia, ib, iorder)%matrix, sab_nl)
               CALL dbcsr_set(diamagnetic_real(ia, ib, iorder)%matrix, 0.0_dp)
               CALL dbcsr_init_p(diamagnetic_imag(ia, ib, iorder)%matrix)
               CALL dbcsr_create(diamagnetic_imag(ia, ib, iorder)%matrix, &
                                 template=path_rs(ia, ib, iorder, 1)%matrix, matrix_type=dbcsr_type_no_symmetry, &
                                 name='ECD GTH DIAMAGNETIC TAYLOR IMAG')
               CALL cp_dbcsr_alloc_block_from_nbl(diamagnetic_imag(ia, ib, iorder)%matrix, sab_nl)
               CALL dbcsr_set(diamagnetic_imag(ia, ib, iorder)%matrix, 0.0_dp)
               source_rs => path_rs(ia:ia, ib, iorder, :)
               CALL rskp_transform(diamagnetic_real(ia, ib, iorder)%matrix, &
                                   diamagnetic_imag(ia, ib, iorder)%matrix, source_rs, 1, xkp_zero, &
                                   cell_to_index, sab_nl)
               IF (iorder == 2) THEN
                  CALL dbcsr_add(diamagnetic_imag(ia, ib, iorder)%matrix, &
                                 diamagnetic_real(ia, ib, iorder)%matrix, 0.0_dp, -1.0_dp)
                  CALL dbcsr_set(diamagnetic_real(ia, ib, iorder)%matrix, 0.0_dp)
               END IF
               NULLIFY (source_rs)
            END DO
         END DO
      END DO

   END SUBROUTINE ecd_periodic_gth_diamagnetic_taylor_build

! **************************************************************************************************
!> \brief Add a complex scaled AO matrix pair to a complex AO matrix pair.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_ao_add_complex_scaled(target_real, target_imag, source_real, source_imag, scale)
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: target_real, target_imag
      TYPE(dbcsr_p_type), INTENT(IN)                     :: source_real, source_imag
      COMPLEX(KIND=dp), INTENT(IN)                       :: scale

      CALL dbcsr_add(target_real%matrix, source_real%matrix, 1.0_dp, REAL(scale, KIND=dp))
      CALL dbcsr_add(target_real%matrix, source_imag%matrix, 1.0_dp, -AIMAG(scale))
      CALL dbcsr_add(target_imag%matrix, source_imag%matrix, 1.0_dp, REAL(scale, KIND=dp))
      CALL dbcsr_add(target_imag%matrix, source_real%matrix, 1.0_dp, AIMAG(scale))

   END SUBROUTINE ecd_periodic_ao_add_complex_scaled

! **************************************************************************************************
!> \brief Return i raised to a non-negative integer power.
! **************************************************************************************************
   PURE FUNCTION ecd_periodic_i_power(order) RESULT(value)
      INTEGER, INTENT(IN)                                :: order
      COMPLEX(KIND=dp)                                   :: value

      SELECT CASE (MODULO(order, 4))
      CASE (0)
         value = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      CASE (1)
         value = CMPLX(0.0_dp, 1.0_dp, KIND=dp)
      CASE (2)
         value = CMPLX(-1.0_dp, 0.0_dp, KIND=dp)
      CASE DEFAULT
         value = CMPLX(0.0_dp, -1.0_dp, KIND=dp)
      END SELECT

   END FUNCTION ecd_periodic_i_power

! **************************************************************************************************
!> \brief Return n factorial as a real scalar.
! **************************************************************************************************
   PURE FUNCTION ecd_periodic_factorial(order) RESULT(value)
      INTEGER, INTENT(IN)                                :: order
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: i

      IF (order < 0) THEN
         value = 0.0_dp
         RETURN
      END IF
      value = 1.0_dp
      DO i = 2, order
         value = value*REAL(i, KIND=dp)
      END DO

   END FUNCTION ecd_periodic_factorial

! **************************************************************************************************
!> \brief Allocate the real symmetric and imaginary antisymmetric parts of one GTH gauge term.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_gauge_nl_pair_create(matrices, template, sab_nl, name)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices
      TYPE(dbcsr_type), POINTER                          :: template
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      CHARACTER(LEN=*), INTENT(IN)                       :: name

      CALL dbcsr_allocate_matrix_set(matrices, 2)
      CALL dbcsr_init_p(matrices(1)%matrix)
      CALL dbcsr_create(matrices(1)%matrix, template=template, &
                        matrix_type=dbcsr_type_symmetric, name=name//" real")
      CALL cp_dbcsr_alloc_block_from_nbl(matrices(1)%matrix, sab_nl)
      CALL dbcsr_set(matrices(1)%matrix, 0.0_dp)
      CALL dbcsr_init_p(matrices(2)%matrix)
      CALL dbcsr_create(matrices(2)%matrix, template=template, &
                        matrix_type=dbcsr_type_antisymmetric, name=name//" imaginary")
      CALL cp_dbcsr_alloc_block_from_nbl(matrices(2)%matrix, sab_nl)
      CALL dbcsr_set(matrices(2)%matrix, 0.0_dp)

   END SUBROUTINE ecd_periodic_gauge_nl_pair_create

! **************************************************************************************************
!> \brief Reset one GTH gauge-term pair before a fresh vector-potential evaluation.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_gauge_nl_pair_reset(matrices)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices

      CPASSERT(SIZE(matrices) == 2)
      CALL dbcsr_set(matrices(1)%matrix, 0.0_dp)
      CALL dbcsr_set(matrices(2)%matrix, 0.0_dp)

   END SUBROUTINE ecd_periodic_gauge_nl_pair_reset

! **************************************************************************************************
!> \brief Form one central finite-difference GTH uniform-gauge curvature matrix.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_gauge_nl_curvature_build(qs_env, nl_zero, nl_pp, nl_pm, nl_mp, nl_mm, &
                                                     step, ia, ib, curvature)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: nl_zero, nl_pp, nl_pm, nl_mp, nl_mm
      REAL(KIND=dp), INTENT(IN)                          :: step
      INTEGER, INTENT(IN)                                :: ia, ib
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: curvature

      REAL(KIND=dp), DIMENSION(ncart)                    :: vector_potential

      CPASSERT(SIZE(nl_zero) == 2)
      CPASSERT(SIZE(nl_pp) == 2)
      CPASSERT(SIZE(nl_pm) == 2)
      CPASSERT(SIZE(nl_mp) == 2)
      CPASSERT(SIZE(nl_mm) == 2)
      CPASSERT(ia >= 1 .AND. ia <= ncart)
      CPASSERT(ib >= 1 .AND. ib <= ncart)
      CPASSERT(step > 0.0_dp)

      CALL dbcsr_set(curvature%matrix, 0.0_dp)
      CALL ecd_periodic_gauge_nl_pair_reset(nl_pp)
      CALL ecd_periodic_gauge_nl_pair_reset(nl_pm)
      CALL ecd_periodic_gauge_nl_pair_reset(nl_mp)
      CALL ecd_periodic_gauge_nl_pair_reset(nl_mm)
      IF (ia == ib) THEN
         vector_potential = 0.0_dp
         vector_potential(ia) = step
         CALL velocity_gauge_nl_term(qs_env, nl_pp, vector_potential)
         vector_potential(ia) = -step
         CALL velocity_gauge_nl_term(qs_env, nl_mm, vector_potential)
         CALL dbcsr_add(curvature%matrix, nl_pp(1)%matrix, 0.0_dp, 1.0_dp/step**2)
         CALL dbcsr_add(curvature%matrix, nl_mm(1)%matrix, 1.0_dp, 1.0_dp/step**2)
         CALL dbcsr_add(curvature%matrix, nl_zero(1)%matrix, 1.0_dp, -2.0_dp/step**2)
      ELSE
         vector_potential = 0.0_dp
         vector_potential(ia) = step
         vector_potential(ib) = step
         CALL velocity_gauge_nl_term(qs_env, nl_pp, vector_potential)
         vector_potential(ib) = -step
         CALL velocity_gauge_nl_term(qs_env, nl_pm, vector_potential)
         vector_potential(ia) = -step
         vector_potential(ib) = step
         CALL velocity_gauge_nl_term(qs_env, nl_mp, vector_potential)
         vector_potential(ib) = -step
         CALL velocity_gauge_nl_term(qs_env, nl_mm, vector_potential)
         CALL dbcsr_add(curvature%matrix, nl_pp(1)%matrix, 0.0_dp, 0.25_dp/step**2)
         CALL dbcsr_add(curvature%matrix, nl_pm(1)%matrix, 1.0_dp, -0.25_dp/step**2)
         CALL dbcsr_add(curvature%matrix, nl_mp(1)%matrix, 1.0_dp, -0.25_dp/step**2)
         CALL dbcsr_add(curvature%matrix, nl_mm(1)%matrix, 1.0_dp, 0.25_dp/step**2)
      END IF

   END SUBROUTINE ecd_periodic_gauge_nl_curvature_build

! **************************************************************************************************
!> \brief Project one complex AO operator, represented by real sparse components, to the MO basis.
! **************************************************************************************************
   SUBROUTINE ecd_project_complex_operator_to_mo(operator_real, operator_imag, mo_coeff, nao, nmo, operator_mo)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: operator_real, operator_imag
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER, INTENT(IN)                                :: nao, nmo
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                    :: operator_mo

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: operator_imag_buffer, operator_real_buffer
      TYPE(cp_fm_type)                                   :: operator_imag_mo, operator_real_mo

      CALL ecd_project_real_operator_to_mo(operator_real, mo_coeff, nao, nmo, operator_real_mo)
      CALL ecd_project_real_operator_to_mo(operator_imag, mo_coeff, nao, nmo, operator_imag_mo)
      ALLOCATE (operator_imag_buffer(nmo, nmo), operator_mo(nmo, nmo), operator_real_buffer(nmo, nmo))
      CALL cp_fm_get_submatrix(operator_real_mo, operator_real_buffer)
      CALL cp_fm_get_submatrix(operator_imag_mo, operator_imag_buffer)
      operator_mo = CMPLX(operator_real_buffer, operator_imag_buffer, KIND=dp)
      DEALLOCATE (operator_imag_buffer, operator_real_buffer)
      CALL cp_fm_release(operator_imag_mo)
      CALL cp_fm_release(operator_real_mo)

   END SUBROUTINE ecd_project_complex_operator_to_mo

END MODULE qs_ecd_periodic
