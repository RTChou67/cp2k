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
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_desymmetrize,&
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
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
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
   USE mathlib,                         ONLY: diag_complex
   USE orbital_pointers,                 ONLY: coset, ncoset
   USE qs_ecd_subspace,                 ONLY: ecd_band_first_order_response,&
                                              ecd_band_manifolds_validate_window,&
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
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE physcon,                         ONLY: c_light_au,&
                                              evolt
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
                                           ecd_gauge_length = 1, &
                                           ecd_gauge_velocity = 2, &
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
!> \param pseudopotential_commutator Include the GTH non-local current
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

      CHARACTER(LEN=default_string_length)               :: gauge_label
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: velocity_response_first, velocity_response_zero
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: length_transition, velocity_transition
      INTEGER                                            :: electric_gauge, npoints, spatial_order
      LOGICAL                                            :: direction_explicit, do_kpoints, &
                                                            energy_range_explicit
      REAL(KIND=dp)                                      :: direction_norm, eps_degeneracy
      REAL(KIND=dp), DIMENSION(:), POINTER               :: direction_ptr, energy_range
      REAL(KIND=dp), DIMENSION(ncart)                    :: direction
      REAL(KIND=dp), DIMENSION(ncart, ncart)             :: diamagnetic_kernel
      TYPE(cell_type), POINTER                           :: cell
      TYPE(ecd_gamma_context_type)                       :: gamma_context
      TYPE(ecd_gamma_zero_order_operators_type)          :: gamma_operators
      TYPE(ecd_periodic_transition_selection_type)       :: transition_selection
      TYPE(section_vals_type), POINTER                   :: print_section, spectrum_section

      NULLIFY (cell, direction_ptr, energy_range, print_section, spectrum_section)

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
      CALL section_vals_val_get(ecd_section, "PROPAGATION_DIRECTION", explicit=direction_explicit)
      IF (direction_explicit) THEN
         CALL section_vals_val_get(ecd_section, "PROPAGATION_DIRECTION", r_vals=direction_ptr)
      END IF
      CALL section_vals_val_get(ipa_section, "ELECTRIC_GAUGE", i_val=electric_gauge)

      IF (eps_degeneracy <= 0.0_dp) THEN
         CPABORT("ECD| EPS_DEGENERACY must be positive")
      END IF
      IF (spatial_order /= 1) THEN
         CPABORT("ECD| Periodic ECD response currently requires SPATIAL_ORDER 1")
      END IF
      IF (.NOT. direction_explicit) THEN
         CPABORT("ECD| Periodic ECD requires PROPAGATION_DIRECTION")
      END IF
      IF (.NOT. ASSOCIATED(direction_ptr)) THEN
         CPABORT("ECD| PROPAGATION_DIRECTION data are unavailable")
      END IF
      IF (SIZE(direction_ptr) /= ncart) THEN
         CPABORT("ECD| PROPAGATION_DIRECTION expects three Cartesian components")
      END IF
      direction = direction_ptr
      direction_norm = SQRT(SUM(direction**2))
      IF (direction_norm <= SQRT(TINY(1.0_dp))) THEN
         CPABORT("ECD| PROPAGATION_DIRECTION must have non-zero length")
      END IF
      direction = direction/direction_norm

      SELECT CASE (electric_gauge)
      CASE (ecd_gauge_length)
         gauge_label = "LENGTH"
      CASE (ecd_gauge_velocity)
         gauge_label = "VELOCITY"
      CASE DEFAULT
         CPABORT("ECD| Unknown periodic ELECTRIC_GAUGE")
      END SELECT

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

      print_section => section_vals_get_subs_vals(ecd_section, "PRINT")
      spectrum_section => section_vals_get_subs_vals(print_section, "SPECTRUM")
      CALL section_vals_val_get(spectrum_section, "ENERGY_RANGE", explicit=energy_range_explicit)
      IF (energy_range_explicit) THEN
         CALL section_vals_val_get(spectrum_section, "ENERGY_RANGE", r_vals=energy_range)
      END IF
      CALL section_vals_val_get(spectrum_section, "NPOINTS", i_val=npoints)
      IF (.NOT. energy_range_explicit) THEN
         CPABORT("ECD| Periodic response requires PRINT%SPECTRUM%ENERGY_RANGE")
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

      CALL ecd_gamma_context_create_implicit(gamma_context, qs_env, direction, ref_point, spatial_order, &
                                             eps_degeneracy)
      CALL ecd_gamma_zero_order_operators_build(gamma_operators, qs_env, gamma_context, &
                                                 pseudopotential_commutator)
      IF (spatial_order >= 1) THEN
         CALL ecd_periodic_implicit_gamma_overlap_g3(gamma_context)
         CALL ecd_periodic_implicit_gamma_current_g4(qs_env, gamma_context, gamma_operators, &
                                                      pseudopotential_commutator)
         CALL ecd_periodic_implicit_gamma_length_g11(qs_env, gamma_context, gamma_operators)
         CALL ecd_periodic_implicit_gamma_ks_g6(qs_env, gamma_context)
         CALL ecd_periodic_implicit_gamma_ks_phase_g7(gamma_context)
         CALL ecd_periodic_implicit_gamma_state_g8(gamma_context)
      END IF
      CALL ecd_periodic_transition_selection_build(transition_selection, gamma_context, ntransitions, &
                                                   transition_band_range, eps_degeneracy)
      CALL ecd_periodic_implicit_gamma_velocity_response_g10(qs_env, ecd_section, gamma_context, gamma_operators, &
                                                              transition_selection, energy_range, npoints, fwhm, &
                                                              electric_gauge == ecd_gauge_velocity, velocity_response_zero, &
                                                              velocity_response_first, diamagnetic_kernel, velocity_transition)
      CALL ecd_periodic_implicit_gamma_length_response_g13(ecd_section, gamma_context, gamma_operators, &
                                                            transition_selection, energy_range, npoints, fwhm, &
                                                            gauge_label, electric_gauge == ecd_gauge_length, &
                                                            velocity_response_zero, velocity_response_first, &
                                                            diamagnetic_kernel, length_transition)
      CALL ecd_periodic_print_transition_matrix_elements(ecd_section, gamma_context, transition_selection, &
                                                          gauge_label, pseudopotential_commutator, &
                                                          velocity_transition, length_transition)
      DEALLOCATE (length_transition, velocity_transition)
      DEALLOCATE (velocity_response_first, velocity_response_zero)
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
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: candidate_iocc, candidate_ispin, candidate_ivirt, order
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
                candidate_ivirt(ncandidates), candidate_occupation_difference(ncandidates), order(ncandidates))
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
      selection%selected_count = nbase
      DO WHILE (selection%selected_count < ncandidates)
         icandidate = order(selection%selected_count + 1)
         IF (ABS(candidate_gap(icandidate) - selection%cutoff_gap) > eps_degeneracy) EXIT
         selection%selected_count = selection%selected_count + 1
      END DO
      IF (nrequested > ncandidates) THEN
         CPWARN('ECD| Requested more periodic IPA transitions than available')
      END IF
      selection%candidate_count = ncandidates
      selection%requested_count = nrequested
      selection%band_range = transition_band_range
      ALLOCATE (selection%gap(selection%selected_count), selection%iocc(selection%selected_count), &
                selection%ispin(selection%selected_count), selection%ivirt(selection%selected_count), &
                selection%occupation_difference(selection%selected_count))
      DO icandidate = 1, selection%selected_count
         nbase = order(icandidate)
         selection%gap(icandidate) = candidate_gap(nbase)
         selection%iocc(icandidate) = candidate_iocc(nbase)
         selection%ispin(icandidate) = candidate_ispin(nbase)
         selection%ivirt(icandidate) = candidate_ivirt(nbase)
         selection%occupation_difference(icandidate) = candidate_occupation_difference(nbase)
      END DO
      DEALLOCATE (candidate_gap, candidate_iocc, candidate_ispin, candidate_ivirt, &
                  candidate_occupation_difference, order)

   END SUBROUTINE ecd_periodic_transition_selection_build

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

      INTEGER                                            :: ispin, nao, nimages, nmo
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-5_dp
      REAL(KIND=dp)                                      :: direct_energy_residual_max, direct_energy_residual_norm, &
                                                            eigenvector_derivative_norm, equation_residual_norm, &
                                                            manifold_residual_norm, normalization_residual_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalue_derivative, eigenvalue_reference, plus_eigenvalues
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_plus, xkp_step, xkp_zero
      TYPE(cp_cfm_type)                                  :: hamiltonian_derivative_mo, hamiltonian_plus, &
                                                            overlap_derivative_mo, overlap_plus, plus_eigenvectors, plus_work
      TYPE(cp_fm_type)                                   :: hamiltonian_derivative_raw, hamiltonian_plus_imag, &
                                                            hamiltonian_plus_real, overlap_derivative_raw, overlap_plus_imag, &
                                                            overlap_plus_real
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: core_images, local_images, overlap_images
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: hamiltonian_plus_imag_ao, hamiltonian_plus_real_ao, &
                                                            overlap_plus_imag_ao, overlap_plus_real_ao
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      NULLIFY (cell_to_index, core_images, hamiltonian_plus_imag_ao, hamiltonian_plus_real_ao, local_images, logger, &
               overlap_images, overlap_plus_imag_ao, overlap_plus_real_ao)
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
      CALL ecd_periodic_ao_set_create(overlap_plus_real_ao, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G8 positive-q overlap real', dbcsr_type_symmetric)
      CALL ecd_periodic_ao_set_create(overlap_plus_imag_ao, 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G8 positive-q overlap imaginary', dbcsr_type_antisymmetric)
      CALL rskp_transform(overlap_plus_real_ao(1)%matrix, overlap_plus_imag_ao(1)%matrix, overlap_images, 1, &
                          xkp_plus, cell_to_index, context%sab_orb)
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

         CALL rskp_transform(hamiltonian_plus_real_ao(ispin)%matrix, hamiltonian_plus_imag_ao(ispin)%matrix, &
                             core_images, ispin, xkp_plus, cell_to_index, context%sab_orb)
         CALL rskp_transform(hamiltonian_plus_real_ao(ispin)%matrix, hamiltonian_plus_imag_ao(ispin)%matrix, &
                             local_images, ispin, xkp_plus, cell_to_index, context%sab_orb)
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
         CALL cp_cfm_create(plus_eigenvectors, context%spin(ispin)%mo_coeff%matrix_struct, &
                            name='ECD G8 positive-q eigenvectors', set_zero=.TRUE.)
         CALL cp_cfm_create(plus_work, hamiltonian_plus%matrix_struct, &
                            name='ECD G8 positive-q diagonalization work', set_zero=.TRUE.)
         ALLOCATE (plus_eigenvalues(nmo))
         CALL cp_cfm_geeig(hamiltonian_plus, overlap_plus, plus_eigenvectors, plus_eigenvalues, plus_work)
         direct_energy_residual_norm = SQRT(SUM((plus_eigenvalues - &
                                                 (eigenvalue_reference + fd_step*eigenvalue_derivative))**2))
         direct_energy_residual_max = MAXVAL(ABS(plus_eigenvalues - &
                                                 (eigenvalue_reference + fd_step*eigenvalue_derivative)))
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
               'ECD| Stage G8 positive-q energy residual spin', ispin, direct_energy_residual_norm
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I2,T55,ES14.6)') &
               'ECD| Stage G8 positive-q maximum residual spin', ispin, direct_energy_residual_max
         END IF
         DEALLOCATE (eigenvalue_derivative, eigenvalue_reference, plus_eigenvalues)
         CALL cp_cfm_release(hamiltonian_derivative_mo)
         CALL cp_cfm_release(overlap_derivative_mo)
         CALL cp_cfm_release(hamiltonian_plus)
         CALL cp_cfm_release(overlap_plus)
         CALL cp_cfm_release(plus_eigenvectors)
         CALL cp_cfm_release(plus_work)
      END DO
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G8 positive-q step [bohr^-1]', fd_step
      END IF
      CALL dbcsr_deallocate_matrix_set(hamiltonian_plus_imag_ao)
      CALL dbcsr_deallocate_matrix_set(hamiltonian_plus_real_ao)
      CALL dbcsr_deallocate_matrix_set(overlap_plus_imag_ao)
      CALL dbcsr_deallocate_matrix_set(overlap_plus_real_ao)

   END SUBROUTINE ecd_periodic_implicit_gamma_state_g8

! **************************************************************************************************
!> \brief Project first-order length vertices into the rotated Gamma MO frame.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_gamma_length_components_build(context, ispin, length_zero_components, &
                                                         length_first_operator_components, &
                                                         length_first_state_components)
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET  :: context
      INTEGER, INTENT(IN)                                :: ispin
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: length_first_operator_components, &
                                                                           length_first_state_components, &
                                                                           length_zero_components

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)    :: length_first_mo, length_first_rotated, &
                                                           length_zero_mo, length_zero_rotated, state_first
      INTEGER                                           :: idir, nmo

      IF (ispin < 1 .OR. ispin > context%nspin) THEN
         CPABORT("ECD| Length transition component construction received an invalid spin index")
      END IF
      IF (.NOT. ASSOCIATED(context%length_taylor_real) .OR. .NOT. ASSOCIATED(context%length_taylor_imag)) THEN
         CPABORT("ECD| Length transition component construction requires retained length Taylor coefficients")
      END IF
      IF (SIZE(context%length_taylor_real, 1) /= ncart .OR. SIZE(context%length_taylor_imag, 1) /= ncart .OR. &
          SIZE(context%length_taylor_real, 2) < 2 .OR. SIZE(context%length_taylor_imag, 2) < 2) THEN
         CPABORT("ECD| Length transition component construction received inconsistent Taylor coefficients")
      END IF
      nmo = context%spin(ispin)%nmo
      IF (nmo < 1 .OR. .NOT. ASSOCIATED(context%spin(ispin)%eigenvalue_first) .OR. &
          .NOT. ASSOCIATED(context%spin(ispin)%eigenvalue_reference) .OR. &
          .NOT. ASSOCIATED(context%spin(ispin)%mo_coeff) .OR. &
          .NOT. ALLOCATED(context%spin(ispin)%band_rotation) .OR. &
          .NOT. ALLOCATED(context%spin(ispin)%coefficient_first)) THEN
         CPABORT("ECD| Length transition component construction requires first-order state data")
      END IF
      IF (SIZE(context%spin(ispin)%eigenvalue_first) /= nmo .OR. &
          SIZE(context%spin(ispin)%eigenvalue_reference) /= nmo .OR. &
          SIZE(context%spin(ispin)%band_rotation, 1) /= nmo .OR. &
          SIZE(context%spin(ispin)%band_rotation, 2) /= nmo .OR. &
          SIZE(context%spin(ispin)%coefficient_first, 1) /= nmo .OR. &
          SIZE(context%spin(ispin)%coefficient_first, 2) /= nmo) THEN
         CPABORT("ECD| Length transition component construction received inconsistent state dimensions")
      END IF

      ALLOCATE (length_first_operator_components(nmo, nmo, ncart), &
                length_first_state_components(nmo, nmo, ncart), length_zero_components(nmo, nmo, ncart))
      DO idir = 1, ncart
         CALL ecd_project_complex_operator_to_mo(context%length_taylor_real(idir, 1), &
                                                 context%length_taylor_imag(idir, 1), &
                                                 context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                 length_zero_mo)
         CALL ecd_project_complex_operator_to_mo(context%length_taylor_real(idir, 2), &
                                                 context%length_taylor_imag(idir, 2), &
                                                 context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                 length_first_mo)
         ALLOCATE (length_first_rotated(nmo, nmo), length_zero_rotated(nmo, nmo), state_first(nmo, nmo))
         length_zero_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                       MATMUL(length_zero_mo, context%spin(ispin)%band_rotation))
         length_first_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                        MATMUL(length_first_mo, context%spin(ispin)%band_rotation))
         state_first = 0.5_dp*(MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%coefficient_first)), &
                                       length_zero_rotated) - &
                                 MATMUL(length_zero_rotated, context%spin(ispin)%coefficient_first))
         length_zero_components(:, :, idir) = length_zero_rotated
         length_first_operator_components(:, :, idir) = length_first_rotated
         length_first_state_components(:, :, idir) = state_first
         DEALLOCATE (length_first_mo, length_first_rotated, length_zero_mo, length_zero_rotated, state_first)
      END DO

   END SUBROUTINE ecd_periodic_gamma_length_components_build

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
      IF (npairs < 1) THEN
         CPABORT("ECD| External translation task list contains no atom pairs")
      END IF

      ALLOCATE (reserve_rows(npairs), reserve_cols(npairs))
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
      DEALLOCATE (reserve_cols, reserve_rows)

   END SUBROUTINE ecd_periodic_task_image_matrix_set_create

! **************************************************************************************************
!> \brief Build implicit-Gamma current Taylor coefficients through first spatial order.
!>
!> The local and GTH image kernels are ECD-owned.  Their zeroth coefficient is
!> checked against the folded complete current before the first coefficient is
!> made available to the finite-q transition provider.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_current_g4(qs_env, context, gamma_operators, &
                                                       pseudopotential_commutator)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(INOUT), TARGET :: context
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: gamma_operators
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator

      INTEGER                                            :: idir
      LOGICAL                                            :: gth_available, gth_enabled
      REAL(KIND=dp)                                      :: first_imag_norm_sq, first_norm_sq, residual_norm_sq, zero_norm_sq
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_zero
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: complete_current, current_residual, &
                                                            gth_current, local_current, position
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: current_taylor_imag, current_taylor_real, &
                                                            gth_current_taylor_imag, gth_current_taylor_real, &
                                                            local_moments_rs, local_momentum_rs, overlap_images
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: gth_current_path_rs, local_moments_der_rs
      TYPE(ecd_periodic_gth_operators_type)              :: gth_operators
      TYPE(ecd_periodic_local_operators_type)            :: local_operators
      INTEGER, DIMENSION(:, :, :), POINTER               :: gth_cell_to_index, local_cell_to_index, &
                                                            translation_cell_to_index

      NULLIFY (complete_current, current_residual, current_taylor_imag, current_taylor_real, &
               gth_cell_to_index, gth_current, gth_current_path_rs, gth_current_taylor_imag, &
               gth_current_taylor_real, local_cell_to_index, local_current, local_moments_der_rs, &
               local_moments_rs, local_momentum_rs, logger, translation_cell_to_index, overlap_images, position)
      xkp_zero = 0.0_dp
      CALL ecd_gamma_context_overlap_images_get(context, overlap_images, translation_cell_to_index)
      IF (ASSOCIATED(context%current_taylor_real) .OR. ASSOCIATED(context%current_taylor_imag)) THEN
         CPABORT("ECD| Stage G4 current Taylor coefficients are already allocated")
      END IF
      CALL ecd_gamma_zero_order_operators_get(gamma_operators, position, local_current, gth_current, complete_current, &
                                               gth_available)

      CALL ecd_periodic_local_operators_build(local_operators, qs_env, context%ref_point, 1)
      CALL ecd_periodic_local_operators_get(local_operators, local_cell_to_index, local_momentum_rs, &
                                            local_moments_der_rs, local_moments_rs)
      gth_enabled = pseudopotential_commutator .AND. gth_available
      IF (gth_enabled) THEN
         CALL ecd_periodic_gth_operators_build(gth_operators, qs_env, context%direction, context%ref_point, 1)
         CALL ecd_periodic_gth_operators_get(gth_operators, gth_cell_to_index, gth_current_path_rs)
      END IF

      CALL ecd_periodic_current_taylor_build(current_taylor_real, current_taylor_imag, &
                                             gth_current_taylor_real, gth_current_taylor_imag, &
                                             local_momentum_rs, local_moments_der_rs, local_moments_rs, &
                                             gth_current_path_rs, gth_enabled, overlap_images, xkp_zero, &
                                             translation_cell_to_index, local_cell_to_index, gth_cell_to_index, &
                                             context%sab_all, context%cell%hmat, context%direction, 1)
      CALL ecd_periodic_ao_vector_create(current_residual, context%matrix_s(1)%matrix, context%sab_all, &
                                         'ECD G4 current zeroth residual', dbcsr_type_no_symmetry)

      logger => cp_get_default_logger()
      DO idir = 1, ncart
         CALL dbcsr_set(current_residual(idir)%matrix, 0.0_dp)
         CALL dbcsr_add(current_residual(idir)%matrix, current_taylor_real(idir, 1)%matrix, 0.0_dp, 1.0_dp)
         CALL dbcsr_add(current_residual(idir)%matrix, complete_current(idir)%matrix, 1.0_dp, -1.0_dp)
         CALL dbcsr_dot(current_taylor_real(idir, 1)%matrix, current_taylor_real(idir, 1)%matrix, zero_norm_sq)
         CALL dbcsr_dot(current_taylor_real(idir, 2)%matrix, current_taylor_real(idir, 2)%matrix, first_norm_sq)
         CALL dbcsr_dot(current_taylor_imag(idir, 2)%matrix, current_taylor_imag(idir, 2)%matrix, &
                        first_imag_norm_sq)
         CALL dbcsr_dot(current_residual(idir)%matrix, current_residual(idir)%matrix, residual_norm_sq)
         IF (cp_logger_get_default_io_unit(logger) > 0) THEN
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current zeroth norm component', idir, SQRT(ABS(zero_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current zeroth residual component', idir, SQRT(ABS(residual_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current first real norm component', idir, SQRT(ABS(first_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G4 current first imaginary norm component', idir, SQRT(ABS(first_imag_norm_sq))
         END IF
      END DO

      CALL dbcsr_deallocate_matrix_set(current_residual)
      context%current_taylor_imag => current_taylor_imag
      context%current_taylor_real => current_taylor_real
      NULLIFY (current_taylor_imag, current_taylor_real)
      CALL dbcsr_deallocate_matrix_set(gth_current_taylor_imag)
      CALL dbcsr_deallocate_matrix_set(gth_current_taylor_real)
      CALL ecd_periodic_gth_operators_release(gth_operators)
      CALL ecd_periodic_local_operators_release(local_operators)

   END SUBROUTINE ecd_periodic_implicit_gamma_current_g4

! **************************************************************************************************
!> \brief Build and validate the first-order length-gauge position vertex at implicit Gamma.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_length_g11(qs_env, context, gamma_operators)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(INOUT), TARGET :: context
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: gamma_operators

      INTEGER                                            :: idir
      REAL(KIND=dp)                                      :: first_imag_norm_sq, first_real_norm_sq, &
                                                            residual_norm_sq, zero_norm_sq
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: complete_current, gth_current, local_current, &
                                                            position, position_general, residual, length_moments

      NULLIFY (complete_current, gth_current, local_current, length_moments, logger, position, position_general, &
               residual)
      IF (ASSOCIATED(context%length_taylor_real) .OR. ASSOCIATED(context%length_taylor_imag)) THEN
         CPABORT("ECD| Stage G11 length Taylor coefficients are already allocated")
      END IF
      CALL ecd_gamma_zero_order_operators_get(gamma_operators, position, local_current, gth_current, complete_current)
      CALL ecd_periodic_ao_set_create(length_moments, ncoset(2) - 1, context%matrix_s(1)%matrix, context%sab_orb, &
                                      'ECD G11 length moments', dbcsr_type_symmetric)
      CALL build_local_moment_matrix(qs_env, length_moments, 2, ref_point=context%ref_point)
      CALL ecd_periodic_length_taylor_build(context%length_taylor_real, context%length_taylor_imag, length_moments, &
                                            context%sab_all, context%matrix_s(1)%matrix, context%direction)

      CALL ecd_periodic_ao_vector_create(position_general, context%matrix_s(1)%matrix, context%sab_all, &
                                         'ECD G11 position general', dbcsr_type_no_symmetry)
      CALL ecd_periodic_ao_vector_create(residual, context%matrix_s(1)%matrix, context%sab_all, &
                                         'ECD G11 length zero residual', dbcsr_type_no_symmetry)
      logger => cp_get_default_logger()
      DO idir = 1, ncart
         CALL dbcsr_desymmetrize(position(idir)%matrix, position_general(idir)%matrix)
         CALL dbcsr_set(residual(idir)%matrix, 0.0_dp)
         CALL dbcsr_add(residual(idir)%matrix, context%length_taylor_real(idir, 1)%matrix, 0.0_dp, 1.0_dp)
         CALL dbcsr_add(residual(idir)%matrix, position_general(idir)%matrix, 1.0_dp, -1.0_dp)
         CALL dbcsr_dot(context%length_taylor_real(idir, 1)%matrix, &
                        context%length_taylor_real(idir, 1)%matrix, zero_norm_sq)
         CALL dbcsr_dot(context%length_taylor_real(idir, 2)%matrix, &
                        context%length_taylor_real(idir, 2)%matrix, first_real_norm_sq)
         CALL dbcsr_dot(context%length_taylor_imag(idir, 2)%matrix, &
                        context%length_taylor_imag(idir, 2)%matrix, first_imag_norm_sq)
         CALL dbcsr_dot(residual(idir)%matrix, residual(idir)%matrix, residual_norm_sq)
         IF (cp_logger_get_default_io_unit(logger) > 0) THEN
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G11 length zero norm component', idir, SQRT(ABS(zero_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G11 length zero residual component', idir, SQRT(ABS(residual_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G11 length first real norm component', idir, SQRT(ABS(first_real_norm_sq))
            WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,I1,T55,ES14.6)') &
               'ECD| Stage G11 length first imaginary norm component', idir, SQRT(ABS(first_imag_norm_sq))
         END IF
      END DO
      CALL dbcsr_deallocate_matrix_set(position_general)
      CALL dbcsr_deallocate_matrix_set(residual)
      CALL dbcsr_deallocate_matrix_set(length_moments)

   END SUBROUTINE ecd_periodic_implicit_gamma_length_g11

! **************************************************************************************************
!> \brief Assemble the first-order velocity-gauge causal response from retained Gamma providers.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_velocity_response_g10(qs_env, ecd_section, context, operators, &
                                                                 transition_selection, energy_range, npoints, fwhm, &
                                                                 write_response, response_zero_output, response_first_output, &
                                                                 diamagnetic_kernel_output, transition_output)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: ecd_section
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET   :: context
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: operators
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: transition_selection
      LOGICAL, INTENT(IN)                                :: write_response
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      INTEGER, INTENT(IN)                                :: npoints
      REAL(KIND=dp), INTENT(IN)                          :: fwhm
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: response_first_output, response_zero_output
      REAL(KIND=dp), DIMENSION(ncart, ncart), INTENT(OUT) :: diamagnetic_kernel_output
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: transition_output

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: current_first_mo, current_first_rotated, &
                                                            current_zero_mo, current_zero_rotated, &
                                                            state_first, transition_first
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: current_first_components, current_zero_components, &
                                                            response_counter, &
                                                            response_first, response_first_counter, response_resonant, &
                                                            response_zero
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: transition_taylor
      COMPLEX(KIND=dp), DIMENSION(ncart, ncart, 0:1)     :: counter_taylor, resonant_taylor
      INTEGER                                            :: iband, idir, ipoint, ispin, itransition, ivirt, nmo, &
                                                            ntransitions, ntransitions_total
      LOGICAL                                            :: gth_available
      REAL(KIND=dp)                                      :: antisymmetric_first_norm, cell_volume, gap, gap_derivative, &
                                                            occupation_difference, photon_energy, prefactor, response_counter_norm, &
                                                            response_first_counter_norm, response_first_norm, &
                                                            response_first_resonant_norm, response_norm, &
                                                            response_resonant_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: gap_taylor
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: complete_current, gth_current, local_current, position

      NULLIFY (complete_current, gth_current, local_current, logger, position)
      IF (.NOT. ASSOCIATED(context%current_taylor_real) .OR. .NOT. ASSOCIATED(context%current_taylor_imag)) THEN
         CPABORT("ECD| Stage G10 requires retained current Taylor coefficients")
      END IF
      IF (SIZE(context%current_taylor_real, 1) /= ncart .OR. SIZE(context%current_taylor_imag, 1) /= ncart .OR. &
          SIZE(context%current_taylor_real, 2) < 2 .OR. SIZE(context%current_taylor_imag, 2) < 2) THEN
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
                                               gth_available)

      ALLOCATE (response_counter(ncart, ncart, npoints), &
                response_first(ncart, ncart, npoints), response_first_counter(ncart, ncart, npoints), &
                response_resonant(ncart, ncart, npoints), response_zero(ncart, ncart, npoints))
      response_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      response_first = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      response_first_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      response_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ALLOCATE (transition_output(ncart, 0:1, transition_selection%selected_count))
      transition_output = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

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
         ALLOCATE (current_first_components(nmo, nmo, ncart), current_zero_components(nmo, nmo, ncart))
         DO idir = 1, ncart
            CALL ecd_project_complex_operator_to_mo(context%current_taylor_real(idir, 1), &
                                                    context%current_taylor_imag(idir, 1), &
                                                    context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                    current_zero_mo)
            CALL ecd_project_complex_operator_to_mo(context%current_taylor_real(idir, 2), &
                                                    context%current_taylor_imag(idir, 2), &
                                                    context%spin(ispin)%mo_coeff, context%spin(ispin)%nao, nmo, &
                                                    current_first_mo)
            ALLOCATE (current_first_rotated(nmo, nmo), current_zero_rotated(nmo, nmo), &
                      state_first(nmo, nmo), transition_first(nmo, nmo))
            current_zero_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                           MATMUL(current_zero_mo, context%spin(ispin)%band_rotation))
            current_first_rotated = MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%band_rotation)), &
                                            MATMUL(current_first_mo, context%spin(ispin)%band_rotation))
            state_first = 0.5_dp*(MATMUL(CONJG(TRANSPOSE(context%spin(ispin)%coefficient_first)), &
                                          current_zero_rotated) - &
                                    MATMUL(current_zero_rotated, context%spin(ispin)%coefficient_first))
            transition_first = current_first_rotated + state_first
            current_zero_components(:, :, idir) = current_zero_rotated
            current_first_components(:, :, idir) = transition_first
            DEALLOCATE (current_first_mo, current_first_rotated, current_zero_mo, current_zero_rotated, &
                        state_first, transition_first)
         END DO

         ALLOCATE (gap_taylor(0:1), transition_taylor(ncart, 0:1))
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
            END DO
            transition_output(:, :, itransition) = transition_taylor
            ntransitions = ntransitions + 1
            CALL ecd_periodic_response_add_transition(transition_taylor, gap_taylor, occupation_difference, &
                                                      cell_volume, energy_range, fwhm, response_resonant, &
                                                      response_counter, response_first, response_first_counter, &
                                                      resonant_taylor, counter_taylor)
         END DO
         ntransitions_total = ntransitions_total + ntransitions
         DEALLOCATE (current_first_components, current_zero_components, gap_taylor, transition_taylor)
      END DO
      IF (ntransitions_total /= transition_selection%selected_count) THEN
         CPABORT('ECD| Stage G10 transition selection is inconsistent with the Gamma spin space')
      END IF

      CALL ecd_periodic_implicit_gamma_diamagnetic_kernel(qs_env, context, gth_available, diamagnetic_kernel_output)
      response_zero = response_resonant + response_counter
      DO ipoint = 1, npoints
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                         (energy_range(2) - energy_range(1))/REAL(npoints - 1, KIND=dp)
         prefactor = 2.0_dp*twopi/(cell_volume*photon_energy**2)
         response_zero(:, :, ipoint) = response_zero(:, :, ipoint) + &
                                        CMPLX(prefactor*diamagnetic_kernel_output, 0.0_dp, KIND=dp)
      END DO
      response_first = response_first + response_first_counter

      ipoint = (npoints + 1)/2
      photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                     (energy_range(2) - energy_range(1))/REAL(npoints - 1, KIND=dp)
      response_resonant_norm = SQRT(SUM(ABS(response_resonant(:, :, ipoint))**2))
      response_counter_norm = SQRT(SUM(ABS(response_counter(:, :, ipoint))**2))
      response_first_resonant_norm = SQRT(SUM(ABS(response_first(:, :, ipoint) - &
                                                       response_first_counter(:, :, ipoint))**2))
      response_first_counter_norm = SQRT(SUM(ABS(response_first_counter(:, :, ipoint))**2))
      response_norm = SQRT(SUM(ABS(response_zero(:, :, ipoint))**2))
      response_first_norm = SQRT(SUM(ABS(response_first(:, :, ipoint))**2))
      antisymmetric_first_norm = 0.5_dp*SQRT(SUM(ABS(response_first(:, :, ipoint) - &
                                                       TRANSPOSE(response_first(:, :, ipoint)))**2))
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
      END IF
      ALLOCATE (response_zero_output(ncart, ncart, npoints), response_first_output(ncart, ncart, npoints))
      response_zero_output = response_zero
      response_first_output = response_first
      IF (write_response) THEN
         CALL ecd_periodic_print_response_tensor(ecd_section, energy_range, response_zero, response_first, &
                                                  context, 'VELOCITY', gth_available, fwhm, transition_selection%band_range, &
                                                  ntransitions_total, candidate_count=transition_selection%candidate_count, &
                                                  requested_count=transition_selection%requested_count, &
                                                  cutoff_gap=transition_selection%cutoff_gap)
         CALL ecd_periodic_print_maxwell_spectrum(ecd_section, energy_range, context%direction, response_zero, &
                                                  response_first, context, 'VELOCITY', gth_available, fwhm, &
                                                  transition_selection%band_range, ntransitions_total, &
                                                  candidate_count=transition_selection%candidate_count, &
                                                  requested_count=transition_selection%requested_count, &
                                                  cutoff_gap=transition_selection%cutoff_gap)
      END IF

      DEALLOCATE (response_counter, response_first, response_first_counter, response_resonant, response_zero)

   END SUBROUTINE ecd_periodic_implicit_gamma_velocity_response_g10

! **************************************************************************************************
!> \brief Assemble the first-order length-gauge causal response and compare it with the velocity tensor.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_length_response_g13(ecd_section, context, operators, &
                                                               transition_selection, energy_range, npoints, fwhm, &
                                                               gauge_label, write_response, velocity_response_zero, &
                                                               velocity_response_first, diamagnetic_kernel, transition_output)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      TYPE(ecd_gamma_context_type), INTENT(IN), TARGET   :: context
      TYPE(ecd_gamma_zero_order_operators_type), INTENT(IN), TARGET :: operators
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: transition_selection
      CHARACTER(LEN=*), INTENT(IN)                       :: gauge_label
      LOGICAL, INTENT(IN)                                :: write_response
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      INTEGER, INTENT(IN)                                :: npoints
      REAL(KIND=dp), INTENT(IN)                          :: fwhm
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: velocity_response_first, velocity_response_zero
      REAL(KIND=dp), DIMENSION(ncart, ncart), INTENT(IN) :: diamagnetic_kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT) :: transition_output

      REAL(KIND=dp), PARAMETER                           :: occupation_threshold = 1.0E-8_dp

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: response_counter, response_first, &
                                                            response_first_counter, response_first_resonant, &
                                                            response_resonant, response_zero
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: length_first_operator_components, &
                                                            length_first_state_components, length_zero_components
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: transition_taylor
      COMPLEX(KIND=dp), DIMENSION(ncart, ncart, 0:1)     :: counter_taylor, resonant_taylor
      INTEGER                                            :: iband, idir, ipoint, ispin, itransition, ivirt, nmo, &
                                                            ntransitions, ntransitions_total
      LOGICAL                                            :: gth_available
      REAL(KIND=dp)                                      :: antisymmetric_first_norm, cell_volume, gap, gap_derivative, &
                                                            gauge_first_residual_max, gauge_first_residual_norm, &
                                                            gauge_first_relative_max, gauge_first_relative_norm, &
                                                            gauge_zero_residual_max, gauge_zero_residual_norm, &
                                                            gauge_zero_relative_max, gauge_zero_relative_norm, &
                                                            occupation_difference, photon_energy, prefactor, &
                                                            response_counter_norm, response_first_counter_norm, &
                                                            response_first_norm, response_first_resonant_norm, &
                                                            response_norm, response_resonant_norm
      REAL(KIND=dp)                                      :: velocity_first_norm, velocity_zero_norm
      REAL(KIND=dp), DIMENSION(0:1)                      :: gap_taylor
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: complete_current, gth_current, local_current, position

      NULLIFY (complete_current, gth_current, local_current, logger, position)
      IF (fwhm <= 0.0_dp) THEN
         CPABORT("ECD| Stage G13 requires a positive FWHM")
      END IF
      IF (SIZE(velocity_response_zero, 1) /= ncart .OR. SIZE(velocity_response_zero, 2) /= ncart .OR. &
          SIZE(velocity_response_first, 1) /= ncart .OR. SIZE(velocity_response_first, 2) /= ncart .OR. &
          SIZE(velocity_response_zero, 3) /= npoints .OR. SIZE(velocity_response_first, 3) /= npoints) THEN
         CPABORT("ECD| Stage G13 received incompatible velocity response tensors")
      END IF
      cell_volume = ABS(context%cell%deth)
      IF (cell_volume <= 0.0_dp) THEN
         CPABORT("ECD| Stage G13 requires a positive cell volume")
      END IF
      CALL ecd_gamma_zero_order_operators_get(operators, position, local_current, gth_current, complete_current, &
                                               gth_available)

      ALLOCATE (response_counter(ncart, ncart, npoints), response_first(ncart, ncart, npoints), &
                response_first_counter(ncart, ncart, npoints), response_first_resonant(ncart, ncart, npoints), &
                response_resonant(ncart, ncart, npoints), response_zero(ncart, ncart, npoints))
      response_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      response_first = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      response_first_counter = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      response_first_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      response_resonant = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ALLOCATE (transition_output(ncart, 0:1, transition_selection%selected_count))
      transition_output = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      logger => cp_get_default_logger()
      ntransitions_total = 0
      DO ispin = 1, context%nspin
         nmo = context%spin(ispin)%nmo
         CALL ecd_periodic_gamma_length_components_build(context, ispin, length_zero_components, &
                                                         length_first_operator_components, &
                                                         length_first_state_components)

         ALLOCATE (transition_taylor(ncart, 0:1))
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
               transition_taylor(idir, 0) = -gap*length_zero_components(ivirt, iband, idir)
               transition_taylor(idir, 1) = -gap_derivative*length_zero_components(ivirt, iband, idir) - &
                                             gap*(length_first_operator_components(ivirt, iband, idir) + &
                                                  length_first_state_components(ivirt, iband, idir))
            END DO
            transition_output(:, :, itransition) = transition_taylor
            ntransitions = ntransitions + 1
            CALL ecd_periodic_response_add_transition(transition_taylor, gap_taylor, occupation_difference, &
                                                      cell_volume, energy_range, fwhm, response_resonant, &
                                                      response_counter, response_first_resonant, &
                                                      response_first_counter, resonant_taylor, counter_taylor)
         END DO
         ntransitions_total = ntransitions_total + ntransitions
         DEALLOCATE (length_first_operator_components, length_first_state_components, length_zero_components, &
                     transition_taylor)
      END DO
      IF (ntransitions_total /= transition_selection%selected_count) THEN
         CPABORT('ECD| Stage G13 transition selection is inconsistent with the Gamma spin space')
      END IF

      response_zero = response_resonant + response_counter
      response_first = response_first_resonant + response_first_counter
      DO ipoint = 1, npoints
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                         (energy_range(2) - energy_range(1))/REAL(npoints - 1, KIND=dp)
         prefactor = 2.0_dp*twopi/(cell_volume*photon_energy**2)
         response_zero(:, :, ipoint) = response_zero(:, :, ipoint) + &
                                        CMPLX(prefactor*diamagnetic_kernel, 0.0_dp, KIND=dp)
      END DO

      gauge_zero_residual_max = 0.0_dp
      gauge_first_residual_max = 0.0_dp
      gauge_zero_relative_max = 0.0_dp
      gauge_first_relative_max = 0.0_dp
      DO ipoint = 1, npoints
         gauge_zero_residual_norm = SQRT(SUM(ABS(response_zero(:, :, ipoint) - &
                                                  velocity_response_zero(:, :, ipoint))**2))
         gauge_first_residual_norm = SQRT(SUM(ABS(response_first(:, :, ipoint) - &
                                                   velocity_response_first(:, :, ipoint))**2))
         velocity_zero_norm = SQRT(SUM(ABS(velocity_response_zero(:, :, ipoint))**2))
         velocity_first_norm = SQRT(SUM(ABS(velocity_response_first(:, :, ipoint))**2))
         response_norm = SQRT(SUM(ABS(response_zero(:, :, ipoint))**2))
         response_first_norm = SQRT(SUM(ABS(response_first(:, :, ipoint))**2))
         gauge_zero_relative_norm = 2.0_dp*gauge_zero_residual_norm/MAX(SQRT(TINY(1.0_dp)), &
                                                                          velocity_zero_norm + response_norm)
         gauge_first_relative_norm = 2.0_dp*gauge_first_residual_norm/MAX(SQRT(TINY(1.0_dp)), &
                                                                            velocity_first_norm + response_first_norm)
         gauge_zero_residual_max = MAX(gauge_zero_residual_max, gauge_zero_residual_norm)
         gauge_first_residual_max = MAX(gauge_first_residual_max, gauge_first_residual_norm)
         gauge_zero_relative_max = MAX(gauge_zero_relative_max, gauge_zero_relative_norm)
         gauge_first_relative_max = MAX(gauge_first_relative_max, gauge_first_relative_norm)
      END DO
      ipoint = (npoints + 1)/2
      photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                      (energy_range(2) - energy_range(1))/REAL(npoints - 1, KIND=dp)
      gauge_zero_residual_norm = SQRT(SUM(ABS(response_zero(:, :, ipoint) - &
                                               velocity_response_zero(:, :, ipoint))**2))
      gauge_first_residual_norm = SQRT(SUM(ABS(response_first(:, :, ipoint) - &
                                                velocity_response_first(:, :, ipoint))**2))
      velocity_zero_norm = SQRT(SUM(ABS(velocity_response_zero(:, :, ipoint))**2))
      velocity_first_norm = SQRT(SUM(ABS(velocity_response_first(:, :, ipoint))**2))
      response_resonant_norm = SQRT(SUM(ABS(response_resonant(:, :, ipoint))**2))
      response_counter_norm = SQRT(SUM(ABS(response_counter(:, :, ipoint))**2))
      response_norm = SQRT(SUM(ABS(response_zero(:, :, ipoint))**2))
      response_first_resonant_norm = SQRT(SUM(ABS(response_first_resonant(:, :, ipoint))**2))
      response_first_counter_norm = SQRT(SUM(ABS(response_first_counter(:, :, ipoint))**2))
      response_first_norm = SQRT(SUM(ABS(response_first(:, :, ipoint))**2))
      gauge_zero_relative_norm = 2.0_dp*gauge_zero_residual_norm/MAX(SQRT(TINY(1.0_dp)), &
                                                                       velocity_zero_norm + response_norm)
      gauge_first_relative_norm = 2.0_dp*gauge_first_residual_norm/MAX(SQRT(TINY(1.0_dp)), &
                                                                         velocity_first_norm + response_first_norm)
      antisymmetric_first_norm = 0.5_dp*SQRT(SUM(ABS(response_first(:, :, ipoint) - &
                                                       TRANSPOSE(response_first(:, :, ipoint)))**2))
      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,I10)') &
            'ECD| Stage G13 length response grid points', npoints
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint photon energy [hartree]', photon_energy
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint resonant zeroth-order norm', response_resonant_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint counter zeroth-order norm', response_counter_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint full zeroth-order norm', response_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint resonant first-order norm', response_first_resonant_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint counter first-order norm', response_first_counter_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint full first-order norm', response_first_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint antisymmetric first-order norm', antisymmetric_first_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint zero-order gauge residual norm', gauge_zero_residual_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint first-order gauge residual norm', gauge_first_residual_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 maximum zero-order gauge residual norm', gauge_zero_residual_max
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 maximum first-order gauge residual norm', gauge_first_residual_max
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint zero-order gauge relative residual', gauge_zero_relative_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 midpoint first-order gauge relative residual', gauge_first_relative_norm
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 maximum zero-order gauge relative residual', gauge_zero_relative_max
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
            'ECD| Stage G13 maximum first-order gauge relative residual', gauge_first_relative_max
      END IF
      CALL ecd_periodic_print_gauge_diagnostics(ecd_section, context, gauge_label, gth_available, fwhm, &
                                                transition_selection%band_range, ntransitions_total, &
                                                transition_selection%candidate_count, transition_selection%requested_count, &
                                                transition_selection%cutoff_gap, energy_range, &
                                                velocity_response_zero, velocity_response_first, response_zero, &
                                                response_first)
      IF (write_response) THEN
         CALL ecd_periodic_print_response_tensor(ecd_section, energy_range, response_zero, response_first, &
                                                  context, gauge_label, gth_available, fwhm, transition_selection%band_range, &
                                                  ntransitions_total, candidate_count=transition_selection%candidate_count, &
                                                  requested_count=transition_selection%requested_count, &
                                                  cutoff_gap=transition_selection%cutoff_gap)
         CALL ecd_periodic_print_maxwell_spectrum(ecd_section, energy_range, context%direction, response_zero, &
                                                  response_first, context, gauge_label, gth_available, fwhm, &
                                                  transition_selection%band_range, ntransitions_total, &
                                                  candidate_count=transition_selection%candidate_count, &
                                                  requested_count=transition_selection%requested_count, &
                                                  cutoff_gap=transition_selection%cutoff_gap)
      END IF
      CALL ecd_periodic_print_output_summary(context, gauge_label, gth_available, ntransitions_total, &
                                             transition_selection%candidate_count, transition_selection%requested_count, &
                                             transition_selection%cutoff_gap, &
                                             gauge_zero_relative_norm, gauge_first_relative_norm, &
                                             gauge_zero_relative_max, gauge_first_relative_max)

      DEALLOCATE (response_counter, response_first, response_first_counter, response_first_resonant, &
                  response_resonant, response_zero)

   END SUBROUTINE ecd_periodic_implicit_gamma_length_response_g13

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
!> \brief Assemble the folded local plus GTH diamagnetic kernel at implicit Gamma.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_implicit_gamma_diamagnetic_kernel(qs_env, context, gth_available, &
                                                              diamagnetic_kernel)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      LOGICAL, INTENT(IN)                                :: gth_available
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
!> \brief Write complete periodic IPA transition coefficients for both electric representations.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_print_transition_matrix_elements(ecd_section, context, transition_selection, &
                                                             gauge_label, gth_requested, &
                                                             velocity_transition, length_transition)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      TYPE(ecd_periodic_transition_selection_type), INTENT(IN) :: transition_selection
      CHARACTER(LEN=*), INTENT(IN)                       :: gauge_label
      LOGICAL, INTENT(IN)                                :: gth_requested
      COMPLEX(KIND=dp), DIMENSION(:, 0:, :), INTENT(IN)  :: length_transition, velocity_transition

      CHARACTER(LEN=default_path_length)                 :: matrix_filename
      INTEGER                                            :: idir, iocc, ispin, itransition, ivirt, matrix_unit
      REAL(KIND=dp)                                      :: gap_derivative
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: matrix_print_key

      IF (SIZE(velocity_transition, 1) /= ncart .OR. SIZE(length_transition, 1) /= ncart .OR. &
          LBOUND(velocity_transition, 2) /= 0 .OR. LBOUND(length_transition, 2) /= 0 .OR. &
          UBOUND(velocity_transition, 2) < 1 .OR. UBOUND(length_transition, 2) < 1 .OR. &
          SIZE(velocity_transition, 3) /= transition_selection%selected_count .OR. &
          SIZE(length_transition, 3) /= transition_selection%selected_count) THEN
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
      WRITE (matrix_unit, '(A)') '# Format_version 1'
      WRITE (matrix_unit, '(A)') '# Gamma_layout IMPLICIT'
      WRITE (matrix_unit, '(A,1X,A)') '# Primary_electric_gauge', TRIM(gauge_label)
      WRITE (matrix_unit, '(A,1X,L1)') '# GTH_nonlocal_current_requested', gth_requested
      WRITE (matrix_unit, '(A,3(1X,ES18.10E3))') '# Reference_point_bohr', context%ref_point
      WRITE (matrix_unit, '(A,3(1X,ES18.10E3))') '# Propagation_direction', context%direction
      WRITE (matrix_unit, '(A,2(1X,I0))') '# Transition_band_range', transition_selection%band_range
      WRITE (matrix_unit, '(A,1X,I0)') '# Requested_occupied_empty_pairs', transition_selection%requested_count
      WRITE (matrix_unit, '(A,1X,I0)') '# Candidate_occupied_empty_pairs', transition_selection%candidate_count
      WRITE (matrix_unit, '(A,1X,I0)') '# Selected_occupied_empty_pairs', transition_selection%selected_count
      WRITE (matrix_unit, '(A,1X,ES18.10E3)') '# Transition_cutoff_gap_hartree', transition_selection%cutoff_gap
      WRITE (matrix_unit, '(A)') '# Velocity and length coefficients use normalized q-Taylor orders 0 and 1.'
      WRITE (matrix_unit, '(A)') '# Length coefficients are reconstructed as velocity-equivalent transitions.'
      WRITE (matrix_unit, '(A)') '# Per Cartesian component: velocity0_re velocity0_im velocity1_re velocity1_im '// &
                                  'length0_re length0_im length1_re length1_im'
      WRITE (matrix_unit, '(A)') '# Columns: Transition_index Spin Occupied_band Virtual_band Occupation_difference '// &
                                  'Gap0_hartree Gap1_hartree followed by X, Y, and Z component groups.'

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
         gap_derivative = 0.5_dp*(context%spin(ispin)%eigenvalue_first(ivirt) + &
                                  context%spin(ispin)%eigenvalue_first(iocc))
         WRITE (matrix_unit, '(4(I0,1X),3(ES18.10E3,1X))', ADVANCE='NO') itransition, ispin, iocc, ivirt, &
            transition_selection%occupation_difference(itransition), transition_selection%gap(itransition), &
            gap_derivative
         DO idir = 1, ncart
            WRITE (matrix_unit, '(8(ES18.10E3,1X))', ADVANCE='NO') &
               REAL(velocity_transition(idir, 0, itransition), KIND=dp), &
               AIMAG(velocity_transition(idir, 0, itransition)), &
               REAL(velocity_transition(idir, 1, itransition), KIND=dp), &
               AIMAG(velocity_transition(idir, 1, itransition)), &
               REAL(length_transition(idir, 0, itransition), KIND=dp), &
               AIMAG(length_transition(idir, 0, itransition)), &
               REAL(length_transition(idir, 1, itransition), KIND=dp), &
               AIMAG(length_transition(idir, 1, itransition))
         END DO
         WRITE (matrix_unit, '()')
      END DO

      IF (cp_logger_get_default_io_unit(logger) > 0) THEN
         WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,1X,A)') &
            'ECD| Periodic transition matrix file', TRIM(matrix_filename)
      END IF
      CALL cp_print_key_finished_output(matrix_unit, logger, ecd_section, 'PRINT%MATRIX_ELEMENTS')

   END SUBROUTINE ecd_periodic_print_transition_matrix_elements

! **************************************************************************************************
!> \brief Write common reproducibility metadata for implicit-Gamma periodic ECD files.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_write_output_metadata(output_unit, context, gauge_label, gth_available, fwhm, &
                                                 transition_band_range, ntransitions, candidate_count, requested_count, &
                                                 cutoff_gap, energy_range, npoints, transition_source)
      INTEGER, INTENT(IN)                                :: output_unit
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      CHARACTER(LEN=*), INTENT(IN)                       :: gauge_label
      LOGICAL, INTENT(IN)                                :: gth_available
      INTEGER, DIMENSION(2), INTENT(IN)                  :: transition_band_range
      INTEGER, INTENT(IN)                                :: candidate_count, ntransitions, requested_count
      INTEGER, INTENT(IN)                                :: npoints
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_gap, fwhm
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: transition_source

      CHARACTER(LEN=12)                                  :: source_label

      IF (output_unit <= 0 .OR. .NOT. ASSOCIATED(context%cell)) THEN
         CPABORT('ECD| Periodic output metadata require an open unit and a cell')
      END IF
      source_label = 'IPA'
      IF (PRESENT(transition_source)) source_label = transition_source
      WRITE (output_unit, '(A,1X,A)') '# Gamma_layout', 'IMPLICIT'
      WRITE (output_unit, '(A,1X,A)') '# Electric_gauge', TRIM(gauge_label)
      WRITE (output_unit, '(A,1X,L1)') '# GTH_current_available', gth_available
      WRITE (output_unit, '(A,1X,A)') '# Broadening', 'LORENTZIAN'
      WRITE (output_unit, '(A,1X,ES18.10E3)') '# FWHM_hartree', fwhm
      WRITE (output_unit, '(A,1X,ES18.10E3)') '# FWHM_eV', fwhm*evolt
      WRITE (output_unit, '(A,1X,I0)') '# Spatial_order', context%spatial_order
      WRITE (output_unit, '(A,3(1X,ES18.10E3))') '# Propagation_direction', context%direction
      WRITE (output_unit, '(A,3(1X,ES18.10E3))') '# Reference_point_bohr', context%ref_point
      WRITE (output_unit, '(A,1X,ES18.10E3)') '# Cell_volume_bohr3', ABS(context%cell%deth)
      WRITE (output_unit, '(A,1X,I0)') '# Spin_channels', context%nspin
      WRITE (output_unit, '(A,2(1X,I0))') '# Transition_band_range', transition_band_range
      SELECT CASE (TRIM(source_label))
      CASE ('IPA')
         WRITE (output_unit, '(A,1X,I0)') '# Requested_occupied_empty_pairs', requested_count
         WRITE (output_unit, '(A,1X,I0)') '# Candidate_occupied_empty_pairs', candidate_count
         WRITE (output_unit, '(A,1X,I0)') '# Selected_occupied_empty_pairs', ntransitions
         WRITE (output_unit, '(A,1X,ES18.10E3)') '# Transition_cutoff_gap_hartree', cutoff_gap
      CASE ('TDDFPT_TDA')
         WRITE (output_unit, '(A,1X,I0)') '# Selected_TDA_states', ntransitions
      CASE DEFAULT
         CPABORT('ECD| Unsupported periodic output transition source')
      END SELECT
      WRITE (output_unit, '(A,2(1X,ES18.10E3))') '# Energy_range_hartree', energy_range
      WRITE (output_unit, '(A,1X,I0)') '# NPOINTS', npoints

   END SUBROUTINE ecd_periodic_write_output_metadata

! **************************************************************************************************
!> \brief Write energy-resolved finite-AO length--velocity residual diagnostics.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_print_gauge_diagnostics(ecd_section, context, gauge_label, gth_available, fwhm, &
                                                    transition_band_range, ntransitions, candidate_count, requested_count, &
                                                    cutoff_gap, energy_range, &
                                                    velocity_zero, velocity_first, length_zero, length_first)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      CHARACTER(LEN=*), INTENT(IN)                       :: gauge_label
      LOGICAL, INTENT(IN)                                :: gth_available
      INTEGER, DIMENSION(2), INTENT(IN)                  :: transition_band_range
      INTEGER, INTENT(IN)                                :: candidate_count, ntransitions, requested_count
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_gap, fwhm
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: length_first, length_zero, velocity_first, velocity_zero

      CHARACTER(LEN=default_path_length)                 :: diagnostics_filename
      INTEGER                                            :: diagnostics_unit, ipoint, iounit
      REAL(KIND=dp)                                      :: absolute_first, absolute_zero, energy, length_first_norm, &
                                                            length_zero_norm, relative_first, relative_zero, &
                                                            velocity_first_norm, velocity_zero_norm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tensor_print_key

      IF (SIZE(velocity_zero, 1) /= ncart .OR. SIZE(velocity_zero, 2) /= ncart .OR. &
          SIZE(velocity_first, 1) /= ncart .OR. SIZE(velocity_first, 2) /= ncart .OR. &
          SIZE(length_zero, 1) /= ncart .OR. SIZE(length_zero, 2) /= ncart .OR. &
          SIZE(length_first, 1) /= ncart .OR. SIZE(length_first, 2) /= ncart .OR. &
          SIZE(velocity_zero, 3) /= SIZE(velocity_first, 3) .OR. &
          SIZE(velocity_zero, 3) /= SIZE(length_zero, 3) .OR. &
          SIZE(velocity_zero, 3) /= SIZE(length_first, 3) .OR. SIZE(velocity_zero, 3) < 2) THEN
         CPABORT('ECD| Gauge diagnostics received incompatible response tensors')
      END IF
      NULLIFY (logger, tensor_print_key)
      logger => cp_get_default_logger()
      tensor_print_key => section_vals_get_subs_vals(ecd_section, 'PRINT%RESPONSE_TENSOR')
      diagnostics_filename = cp_print_key_generate_filename(logger, tensor_print_key, &
                                                            extension='.ecd_gauge', my_local=.FALSE.)
      diagnostics_unit = cp_print_key_unit_nr(logger, ecd_section, 'PRINT%RESPONSE_TENSOR', &
                                               extension='.ecd_gauge', file_position='REWIND', &
                                               file_action='WRITE', file_form='FORMATTED')
      IF (diagnostics_unit <= 0) RETURN

      iounit = cp_logger_get_default_io_unit(logger)
      IF (iounit > 0) THEN
         WRITE (iounit, '(T2,A,1X,A)') 'ECD| Periodic gauge diagnostics file', TRIM(diagnostics_filename)
      END IF
      WRITE (diagnostics_unit, '(A)') '# CP2K periodic ECD gauge diagnostics'
      WRITE (diagnostics_unit, '(A)') '# Format_version 2'
      CALL ecd_periodic_write_output_metadata(diagnostics_unit, context, gauge_label, gth_available, fwhm, &
                                              transition_band_range, ntransitions, candidate_count, requested_count, &
                                              cutoff_gap, energy_range, SIZE(velocity_zero, 3))
      WRITE (diagnostics_unit, '(A)') '# Relative residual: 2*norm(velocity-length)/(norm(velocity)+norm(length))'
      WRITE (diagnostics_unit, '(A)') '# Columns: Energy_hartree Energy_eV velocity0_norm length0_norm '// &
                                     'absolute0_residual relative0_residual velocity1_norm length1_norm '// &
                                     'absolute1_residual relative1_residual'
      DO ipoint = 1, SIZE(velocity_zero, 3)
         energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                  (energy_range(2) - energy_range(1))/REAL(SIZE(velocity_zero, 3) - 1, KIND=dp)
         velocity_zero_norm = SQRT(SUM(ABS(velocity_zero(:, :, ipoint))**2))
         length_zero_norm = SQRT(SUM(ABS(length_zero(:, :, ipoint))**2))
         absolute_zero = SQRT(SUM(ABS(velocity_zero(:, :, ipoint) - length_zero(:, :, ipoint))**2))
         relative_zero = 2.0_dp*absolute_zero/MAX(SQRT(TINY(1.0_dp)), velocity_zero_norm + length_zero_norm)
         velocity_first_norm = SQRT(SUM(ABS(velocity_first(:, :, ipoint))**2))
         length_first_norm = SQRT(SUM(ABS(length_first(:, :, ipoint))**2))
         absolute_first = SQRT(SUM(ABS(velocity_first(:, :, ipoint) - length_first(:, :, ipoint))**2))
         relative_first = 2.0_dp*absolute_first/MAX(SQRT(TINY(1.0_dp)), velocity_first_norm + length_first_norm)
         WRITE (diagnostics_unit, '(10(1X,ES18.10E3))') energy, energy*evolt, velocity_zero_norm, length_zero_norm, &
            absolute_zero, relative_zero, velocity_first_norm, length_first_norm, absolute_first, relative_first
      END DO
      CALL cp_print_key_finished_output(diagnostics_unit, logger, ecd_section, 'PRINT%RESPONSE_TENSOR')

   END SUBROUTINE ecd_periodic_print_gauge_diagnostics

! **************************************************************************************************
!> \brief Print the compact final implicit-Gamma ECD validation summary.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_print_output_summary(context, gauge_label, gth_available, ntransitions, candidate_count, &
                                                requested_count, cutoff_gap, &
                                                midpoint_zero_relative, midpoint_first_relative, &
                                                maximum_zero_relative, maximum_first_relative)
      TYPE(ecd_gamma_context_type), INTENT(IN)           :: context
      CHARACTER(LEN=*), INTENT(IN)                       :: gauge_label
      LOGICAL, INTENT(IN)                                :: gth_available
      INTEGER, INTENT(IN)                                :: candidate_count, ntransitions, requested_count
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_gap, maximum_first_relative, maximum_zero_relative, &
                                                            midpoint_first_relative, midpoint_zero_relative

      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (cp_logger_get_default_io_unit(logger) <= 0) RETURN
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A)') 'ECD| Periodic output summary'
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,A)') &
         'ECD| Primary electric gauge', TRIM(gauge_label)
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,L1)') &
         'ECD| GTH current available', gth_available
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,I10)') &
         'ECD| Requested occupied-empty pairs', requested_count
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,I10)') &
         'ECD| Candidate occupied-empty pairs', candidate_count
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,I10)') &
         'ECD| Selected occupied-empty pairs', ntransitions
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
         'ECD| Transition cutoff gap [hartree]', cutoff_gap
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
         'ECD| Midpoint zero-order gauge relative residual', midpoint_zero_relative
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
         'ECD| Midpoint first-order gauge relative residual', midpoint_first_relative
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
         'ECD| Maximum zero-order gauge relative residual', maximum_zero_relative
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,ES14.6)') &
         'ECD| Maximum first-order gauge relative residual', maximum_first_relative
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,3F12.6)') &
         'ECD| Reference point [bohr]', context%ref_point
      WRITE (cp_logger_get_default_io_unit(logger), '(T2,A,T55,3F12.6)') &
         'ECD| Propagation direction', context%direction

   END SUBROUTINE ecd_periodic_print_output_summary

! **************************************************************************************************
!> \brief Write the explicit-Gamma zero- and first-order dielectric response tensors.
!> \param ecd_section ECD input section supplying PRINT%RESPONSE_TENSOR
!> \param energy_range Photon-energy range in hartree
!> \param response_zero Full zero-order response contribution, excluding vacuum identity
!> \param response_first Directional first-order dielectric coefficient
! **************************************************************************************************
   SUBROUTINE ecd_periodic_print_response_tensor(ecd_section, energy_range, response_zero, response_first, &
                                                 context, gauge_label, gth_available, fwhm, transition_band_range, &
                                                 ntransitions, transition_source, candidate_count, requested_count, cutoff_gap)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: response_first, response_zero
      TYPE(ecd_gamma_context_type), INTENT(IN), OPTIONAL :: context
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: gauge_label
      LOGICAL, INTENT(IN), OPTIONAL                      :: gth_available
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: transition_band_range
      INTEGER, INTENT(IN), OPTIONAL                      :: candidate_count, ntransitions, requested_count
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff_gap, fwhm
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: transition_source

      CHARACTER(LEN=default_path_length)                 :: response_filename
      CHARACTER(LEN=12)                                  :: source_label
      COMPLEX(KIND=dp)                                   :: epsilon_zero
      INTEGER                                            :: ia, ib, ipoint, iounit, response_unit
      REAL(KIND=dp)                                      :: photon_energy
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tensor_print_key

      CPASSERT(SIZE(response_zero, 1) == ncart)
      CPASSERT(SIZE(response_zero, 2) == ncart)
      CPASSERT(SIZE(response_first, 1) == ncart)
      CPASSERT(SIZE(response_first, 2) == ncart)
      CPASSERT(SIZE(response_zero, 3) == SIZE(response_first, 3))
      CPASSERT(SIZE(response_zero, 3) >= 2)

      NULLIFY (logger, tensor_print_key)
      source_label = 'IPA'
      IF (PRESENT(transition_source)) source_label = transition_source
      logger => cp_get_default_logger()
      tensor_print_key => section_vals_get_subs_vals(ecd_section, "PRINT%RESPONSE_TENSOR")
      response_filename = cp_print_key_generate_filename(logger, tensor_print_key, &
                                                         extension=".ecd_response", my_local=.FALSE.)
      response_unit = cp_print_key_unit_nr(logger, ecd_section, "PRINT%RESPONSE_TENSOR", &
                                            extension=".ecd_response", file_position="REWIND", &
                                            file_action="WRITE", file_form="FORMATTED")
      IF (response_unit <= 0) RETURN

      iounit = cp_logger_get_default_io_unit(logger)
      IF (iounit > 0) THEN
         WRITE (iounit, '(T2,A,1X,A)') 'ECD| Periodic response tensor file', TRIM(response_filename)
      END IF
      WRITE (response_unit, '(A)') '# CP2K periodic ECD dielectric response'
      IF (PRESENT(context)) THEN
         IF (.NOT. PRESENT(gauge_label) .OR. .NOT. PRESENT(gth_available) .OR. .NOT. PRESENT(fwhm) .OR. &
             .NOT. PRESENT(transition_band_range) .OR. .NOT. PRESENT(ntransitions) .OR. &
             .NOT. PRESENT(candidate_count) .OR. .NOT. PRESENT(requested_count) .OR. .NOT. PRESENT(cutoff_gap)) THEN
            CPABORT('ECD| Periodic response metadata are incomplete')
         END IF
         WRITE (response_unit, '(A)') '# Format_version 3'
         CALL ecd_periodic_write_output_metadata(response_unit, context, gauge_label, gth_available, fwhm, &
                                                 transition_band_range, ntransitions, candidate_count, requested_count, &
                                                 cutoff_gap, energy_range, SIZE(response_zero, 3), &
                                                 source_label)
      END IF
      WRITE (response_unit, '(A)') '# epsilon0 includes the vacuum identity and diamagnetic baseline'
      WRITE (response_unit, '(A)') '# epsilon1 is the directional Taylor coefficient in bohr'
      WRITE (response_unit, '(A)', ADVANCE='NO') '# Columns: Energy_hartree Energy_eV'
      DO ia = 1, ncart
         DO ib = 1, ncart
            WRITE (response_unit, '(1X,A,I1,I1,A)', ADVANCE='NO') 'epsilon0_', ia, ib, '_re'
            WRITE (response_unit, '(1X,A,I1,I1,A)', ADVANCE='NO') 'epsilon0_', ia, ib, '_im'
         END DO
      END DO
      DO ia = 1, ncart
         DO ib = 1, ncart
            WRITE (response_unit, '(1X,A,I1,I1,A)', ADVANCE='NO') 'epsilon1_', ia, ib, '_re'
            WRITE (response_unit, '(1X,A,I1,I1,A)', ADVANCE='NO') 'epsilon1_', ia, ib, '_im'
         END DO
      END DO
      WRITE (response_unit, '()')

      DO ipoint = 1, SIZE(response_zero, 3)
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                          (energy_range(2) - energy_range(1))/REAL(SIZE(response_zero, 3) - 1, KIND=dp)
         WRITE (response_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') photon_energy, photon_energy*evolt
         DO ia = 1, ncart
            DO ib = 1, ncart
               epsilon_zero = response_zero(ia, ib, ipoint)
               IF (ia == ib) epsilon_zero = epsilon_zero + CMPLX(1.0_dp, 0.0_dp, KIND=dp)
               WRITE (response_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') &
                  REAL(epsilon_zero, KIND=dp), AIMAG(epsilon_zero)
            END DO
         END DO
         DO ia = 1, ncart
            DO ib = 1, ncart
               WRITE (response_unit, '(2(1X,ES18.10E3))', ADVANCE='NO') &
                  REAL(response_first(ia, ib, ipoint), KIND=dp), AIMAG(response_first(ia, ib, ipoint))
            END DO
         END DO
         WRITE (response_unit, '()')
      END DO

      CALL cp_print_key_finished_output(response_unit, logger, ecd_section, "PRINT%RESPONSE_TENSOR")

   END SUBROUTINE ecd_periodic_print_response_tensor

! **************************************************************************************************
!> \brief Solve the first-order transverse Maxwell equation and write a directional ECD spectrum.
!> \param ecd_section ECD input section supplying PRINT%SPECTRUM
!> \param energy_range Photon-energy range in hartree
!> \param direction Normalized propagation direction
!> \param response_zero Zero-order response contribution, excluding vacuum identity
!> \param response_first Directional first-order dielectric coefficient
! **************************************************************************************************
   SUBROUTINE ecd_periodic_print_maxwell_spectrum(ecd_section, energy_range, direction, response_zero, &
                                                   response_first, context, gauge_label, gth_available, fwhm, &
                                                   transition_band_range, ntransitions, transition_source, candidate_count, &
                                                   requested_count, cutoff_gap)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: energy_range
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: direction
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: response_first, response_zero
      TYPE(ecd_gamma_context_type), INTENT(IN), OPTIONAL :: context
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: gauge_label
      LOGICAL, INTENT(IN), OPTIONAL                      :: gth_available
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: transition_band_range
      INTEGER, INTENT(IN), OPTIONAL                      :: candidate_count, ntransitions, requested_count
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff_gap, fwhm
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: transition_source

      CHARACTER(LEN=default_path_length)                 :: spectrum_filename
      CHARACTER(LEN=12)                                  :: source_label
      COMPLEX(KIND=dp), DIMENSION(2, 2)                  :: epsilon_first_transverse, &
                                                            epsilon_zero_transverse, maxwell_matrix, &
                                                            polarizations
      COMPLEX(KIND=dp), DIMENSION(4)                     :: roots
      COMPLEX(KIND=dp)                                   :: q_l, q_r
      INTEGER, DIMENSION(2)                              :: forward_roots
      INTEGER                                            :: imode, ipoint, iounit, left_index, right_index, spectrum_unit
      REAL(KIND=dp)                                      :: alpha_l, alpha_r, beta, label_separation, &
                                                            max_residual, min_label_separation, photon_energy, &
                                                            residuals(2)
      REAL(KIND=dp), DIMENSION(ncart, 2)                 :: transverse_basis
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: spectrum_print_key

      CPASSERT(SIZE(response_zero, 1) == ncart)
      CPASSERT(SIZE(response_zero, 2) == ncart)
      CPASSERT(SIZE(response_first, 1) == ncart)
      CPASSERT(SIZE(response_first, 2) == ncart)
      CPASSERT(SIZE(response_zero, 3) == SIZE(response_first, 3))
      CPASSERT(SIZE(response_zero, 3) >= 2)

      CALL ecd_periodic_transverse_basis(direction, transverse_basis)
      NULLIFY (logger, spectrum_print_key)
      source_label = 'IPA'
      IF (PRESENT(transition_source)) source_label = transition_source
      logger => cp_get_default_logger()
      spectrum_print_key => section_vals_get_subs_vals(ecd_section, "PRINT%SPECTRUM")
      spectrum_filename = cp_print_key_generate_filename(logger, spectrum_print_key, &
                                                         extension=".ecd", my_local=.FALSE.)
      spectrum_unit = cp_print_key_unit_nr(logger, ecd_section, "PRINT%SPECTRUM", &
                                            extension=".ecd", file_position="REWIND", &
                                            file_action="WRITE", file_form="FORMATTED")
      IF (spectrum_unit <= 0) RETURN

      iounit = cp_logger_get_default_io_unit(logger)
      WRITE (spectrum_unit, '(A)') '# CP2K periodic ECD Maxwell spectrum'
      IF (PRESENT(context)) THEN
         IF (.NOT. PRESENT(gauge_label) .OR. .NOT. PRESENT(gth_available) .OR. .NOT. PRESENT(fwhm) .OR. &
             .NOT. PRESENT(transition_band_range) .OR. .NOT. PRESENT(ntransitions) .OR. &
             .NOT. PRESENT(candidate_count) .OR. .NOT. PRESENT(requested_count) .OR. .NOT. PRESENT(cutoff_gap)) THEN
            CPABORT('ECD| Maxwell-spectrum metadata are incomplete')
         END IF
         WRITE (spectrum_unit, '(A)') '# Format_version 3'
         CALL ecd_periodic_write_output_metadata(spectrum_unit, context, gauge_label, gth_available, fwhm, &
                                                 transition_band_range, ntransitions, candidate_count, requested_count, &
                                                 cutoff_gap, energy_range, SIZE(response_zero, 3), &
                                                 source_label)
      END IF
      WRITE (spectrum_unit, '(A,3(1X,ES18.10E3))') '# Propagation_direction', direction
      WRITE (spectrum_unit, '(A,3(1X,ES18.10E3))') '# Transverse_e1', transverse_basis(:, 1)
      WRITE (spectrum_unit, '(A,3(1X,ES18.10E3))') '# Transverse_e2', transverse_basis(:, 2)
      WRITE (spectrum_unit, '(A)') '# Circular convention: eL=(e1+i*e2)/sqrt(2), eR=(e1-i*e2)/sqrt(2)'
      WRITE (spectrum_unit, '(A)') '# Columns: Energy_hartree Energy_eV q_L_re q_L_im q_R_re q_R_im '// &
                                   'alpha_L_bohr-1 alpha_R_bohr-1 alpha_mean_bohr-1 Delta_alpha_bohr-1 '// &
                                   'Maxwell_residual_L Maxwell_residual_R L_label_separation'

      max_residual = 0.0_dp
      min_label_separation = HUGE(1.0_dp)
      DO ipoint = 1, SIZE(response_zero, 3)
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                         (energy_range(2) - energy_range(1))/REAL(SIZE(response_zero, 3) - 1, KIND=dp)
         beta = (photon_energy/c_light_au)**2
         CALL ecd_periodic_project_transverse(response_zero(:, :, ipoint), transverse_basis, &
                                              epsilon_zero_transverse)
         epsilon_zero_transverse(1, 1) = epsilon_zero_transverse(1, 1) + CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         epsilon_zero_transverse(2, 2) = epsilon_zero_transverse(2, 2) + CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         CALL ecd_periodic_project_transverse(response_first(:, :, ipoint), transverse_basis, &
                                              epsilon_first_transverse)
         CALL ecd_periodic_maxwell_roots(beta*epsilon_zero_transverse, beta*epsilon_first_transverse, roots)
         CALL ecd_periodic_select_forward_roots(roots, forward_roots)
         DO imode = 1, SIZE(forward_roots)
            maxwell_matrix = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            maxwell_matrix(1, 1) = roots(forward_roots(imode))**2
            maxwell_matrix(2, 2) = roots(forward_roots(imode))**2
            maxwell_matrix = maxwell_matrix - roots(forward_roots(imode))*beta*epsilon_first_transverse - &
                             beta*epsilon_zero_transverse
            CALL ecd_periodic_maxwell_nullvector(maxwell_matrix, polarizations(:, imode), residuals(imode))
         END DO
         CALL ecd_periodic_assign_circular_modes(polarizations, left_index, right_index, label_separation)
         q_l = roots(forward_roots(left_index))
         q_r = roots(forward_roots(right_index))
         alpha_l = 2.0_dp*AIMAG(q_l)
         alpha_r = 2.0_dp*AIMAG(q_r)
         max_residual = MAX(max_residual, MAXVAL(residuals))
         min_label_separation = MIN(min_label_separation, label_separation)
         WRITE (spectrum_unit, '(13(1X,ES18.10E3))') photon_energy, photon_energy*evolt, &
            REAL(q_l, KIND=dp), AIMAG(q_l), REAL(q_r, KIND=dp), AIMAG(q_r), alpha_l, alpha_r, &
            0.5_dp*(alpha_l + alpha_r), alpha_l - alpha_r, residuals(left_index), residuals(right_index), &
            label_separation
      END DO

      IF (iounit > 0) THEN
         ipoint = (SIZE(response_zero, 3) + 1)/2
         photon_energy = energy_range(1) + REAL(ipoint - 1, KIND=dp)* &
                         (energy_range(2) - energy_range(1))/REAL(SIZE(response_zero, 3) - 1, KIND=dp)
         WRITE (iounit, '(T2,A,T55,ES14.6)') 'ECD| Maxwell midpoint photon energy [hartree]', photon_energy
         WRITE (iounit, '(T2,A,T55,ES14.6)') 'ECD| Maxwell maximum residual', max_residual
         WRITE (iounit, '(T2,A,T55,ES14.6)') 'ECD| Maxwell minimum L-label separation', min_label_separation
         WRITE (iounit, '(T2,A,1X,A)') 'ECD| Periodic Maxwell spectrum file', TRIM(spectrum_filename)
      END IF
      CALL cp_print_key_finished_output(spectrum_unit, logger, ecd_section, "PRINT%SPECTRUM")

   END SUBROUTINE ecd_periodic_print_maxwell_spectrum

! **************************************************************************************************
!> \brief Construct a right-handed real transverse basis for a propagation direction.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_transverse_basis(direction, transverse_basis)
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: direction
      REAL(KIND=dp), DIMENSION(ncart, 2), INTENT(OUT)    :: transverse_basis

      INTEGER                                            :: seed_index
      REAL(KIND=dp)                                      :: norm_e1, projection
      REAL(KIND=dp), DIMENSION(ncart)                    :: seed

      seed_index = MINLOC(ABS(direction), DIM=1)
      seed = 0.0_dp
      seed(seed_index) = 1.0_dp
      projection = DOT_PRODUCT(seed, direction)
      transverse_basis(:, 1) = seed - projection*direction
      norm_e1 = SQRT(SUM(transverse_basis(:, 1)**2))
      IF (norm_e1 <= SQRT(TINY(1.0_dp))) THEN
         CPABORT('ECD| Unable to construct a transverse Maxwell basis')
      END IF
      transverse_basis(:, 1) = transverse_basis(:, 1)/norm_e1
      transverse_basis(1, 2) = direction(2)*transverse_basis(3, 1) - direction(3)*transverse_basis(2, 1)
      transverse_basis(2, 2) = direction(3)*transverse_basis(1, 1) - direction(1)*transverse_basis(3, 1)
      transverse_basis(3, 2) = direction(1)*transverse_basis(2, 1) - direction(2)*transverse_basis(1, 1)

   END SUBROUTINE ecd_periodic_transverse_basis

! **************************************************************************************************
!> \brief Project a Cartesian tensor into a real two-dimensional transverse basis.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_project_transverse(cartesian_tensor, transverse_basis, transverse_tensor)
      COMPLEX(KIND=dp), DIMENSION(ncart, ncart), INTENT(IN) :: cartesian_tensor
      REAL(KIND=dp), DIMENSION(ncart, 2), INTENT(IN)     :: transverse_basis
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(OUT)     :: transverse_tensor

      INTEGER                                            :: ia, ib, it, jt

      transverse_tensor = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO it = 1, 2
         DO jt = 1, 2
            DO ia = 1, ncart
               DO ib = 1, ncart
                  transverse_tensor(it, jt) = transverse_tensor(it, jt) + &
                                               transverse_basis(ia, it)*cartesian_tensor(ia, ib)* &
                                               transverse_basis(ib, jt)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE ecd_periodic_project_transverse

! **************************************************************************************************
!> \brief Solve the quartic first-order transverse Maxwell dispersion polynomial.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_maxwell_roots(maxwell_zero, maxwell_first, roots)
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(IN)      :: maxwell_first, maxwell_zero
      COMPLEX(KIND=dp), DIMENSION(4), INTENT(OUT)        :: roots

      COMPLEX(KIND=dp)                                   :: a0, a1, a2, a3
      COMPLEX(KIND=dp), DIMENSION(4, 4)                  :: companion
      COMPLEX(KIND=dp), DIMENSION(1, 1)                  :: left_dummy, right_dummy
      COMPLEX(KIND=dp), DIMENSION(8)                     :: work
      INTEGER                                            :: info
      REAL(KIND=dp), DIMENSION(8)                        :: rwork

      EXTERNAL                                           :: zgeev

      a3 = -(maxwell_first(1, 1) + maxwell_first(2, 2))
      a2 = maxwell_first(1, 1)*maxwell_first(2, 2) - maxwell_first(1, 2)*maxwell_first(2, 1) - &
           maxwell_zero(1, 1) - maxwell_zero(2, 2)
      a1 = maxwell_first(1, 1)*maxwell_zero(2, 2) + maxwell_first(2, 2)*maxwell_zero(1, 1) - &
           maxwell_first(1, 2)*maxwell_zero(2, 1) - maxwell_first(2, 1)*maxwell_zero(1, 2)
      a0 = maxwell_zero(1, 1)*maxwell_zero(2, 2) - maxwell_zero(1, 2)*maxwell_zero(2, 1)

      companion = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      companion(1, 4) = -a0
      companion(2, 1) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      companion(2, 4) = -a1
      companion(3, 2) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      companion(3, 4) = -a2
      companion(4, 3) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      companion(4, 4) = -a3
      CALL zgeev('N', 'N', 4, companion, 4, roots, left_dummy, 1, right_dummy, 1, work, SIZE(work), rwork, info)
      IF (info /= 0) THEN
         CPABORT('ECD| Complex Maxwell root solver failed')
      END IF

   END SUBROUTINE ecd_periodic_maxwell_roots

! **************************************************************************************************
!> \brief Identify the two roots that propagate along the requested positive direction.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_select_forward_roots(roots, forward_roots)
      COMPLEX(KIND=dp), DIMENSION(4), INTENT(IN)         :: roots
      INTEGER, DIMENSION(2), INTENT(OUT)                 :: forward_roots

      INTEGER                                            :: iroot, nforward
      REAL(KIND=dp)                                      :: root_scale, tolerance

      root_scale = MAX(1.0_dp, MAXVAL(ABS(roots)))
      tolerance = 100.0_dp*EPSILON(1.0_dp)*root_scale
      nforward = 0
      DO iroot = 1, SIZE(roots)
         IF (REAL(roots(iroot), KIND=dp) > tolerance .OR. &
             (ABS(REAL(roots(iroot), KIND=dp)) <= tolerance .AND. AIMAG(roots(iroot)) > tolerance)) THEN
            nforward = nforward + 1
            IF (nforward <= SIZE(forward_roots)) forward_roots(nforward) = iroot
         END IF
      END DO
      IF (nforward /= SIZE(forward_roots)) THEN
         CPABORT('ECD| Maxwell dispersion did not yield two forward propagation roots')
      END IF

   END SUBROUTINE ecd_periodic_select_forward_roots

! **************************************************************************************************
!> \brief Construct a normalized null vector for a two-dimensional Maxwell matrix.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_maxwell_nullvector(maxwell_matrix, polarization, residual)
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(IN)      :: maxwell_matrix
      COMPLEX(KIND=dp), DIMENSION(2), INTENT(OUT)        :: polarization
      REAL(KIND=dp), INTENT(OUT)                         :: residual

      COMPLEX(KIND=dp), DIMENSION(2)                     :: candidate
      REAL(KIND=dp)                                      :: first_norm, second_norm

      first_norm = SQRT(SUM(ABS(maxwell_matrix(1, :))**2))
      second_norm = SQRT(SUM(ABS(maxwell_matrix(2, :))**2))
      IF (first_norm >= second_norm) THEN
         candidate = [-maxwell_matrix(1, 2), maxwell_matrix(1, 1)]
      ELSE
         candidate = [-maxwell_matrix(2, 2), maxwell_matrix(2, 1)]
      END IF
      IF (SQRT(SUM(ABS(candidate)**2)) <= SQRT(TINY(1.0_dp))) THEN
         CPABORT('ECD| Maxwell matrix has no stable transverse null vector')
      END IF
      polarization = candidate/SQRT(SUM(ABS(candidate)**2))
      residual = SQRT(SUM(ABS(MATMUL(maxwell_matrix, polarization))**2))

   END SUBROUTINE ecd_periodic_maxwell_nullvector

! **************************************************************************************************
!> \brief Assign the two forward transverse modes by circular-polarization overlap.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_assign_circular_modes(polarizations, left_index, right_index, label_separation)
      COMPLEX(KIND=dp), DIMENSION(2, 2), INTENT(IN)      :: polarizations
      INTEGER, INTENT(OUT)                               :: left_index, right_index
      REAL(KIND=dp), INTENT(OUT)                         :: label_separation

      COMPLEX(KIND=dp), DIMENSION(2)                     :: circular_left
      REAL(KIND=dp), DIMENSION(2)                        :: left_overlap

      circular_left(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp, KIND=dp)
      circular_left(2) = CMPLX(0.0_dp, 1.0_dp/SQRT(2.0_dp), KIND=dp)
      left_overlap(1) = ABS(DOT_PRODUCT(circular_left, polarizations(:, 1)))**2
      left_overlap(2) = ABS(DOT_PRODUCT(circular_left, polarizations(:, 2)))**2
      IF (left_overlap(1) >= left_overlap(2)) THEN
         left_index = 1
         right_index = 2
      ELSE
         left_index = 2
         right_index = 1
      END IF
      label_separation = ABS(left_overlap(1) - left_overlap(2))

   END SUBROUTINE ecd_periodic_assign_circular_modes

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
                                                gth_current_path_rs, pseudopotential_commutator, matrix_s, xkp, &
                                                field_cell_to_index, local_cell_to_index, gth_cell_to_index, sab_nl, &
                                                hmat, direction, &
                                                max_order)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: current_real, current_imag, gth_real, gth_imag
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: local_momentum_rs, local_moments_rs, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: gth_current_path_rs, local_moments_der_rs
      LOGICAL, INTENT(IN)                                 :: pseudopotential_commutator
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
      IF (pseudopotential_commutator) THEN
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
               IF (pseudopotential_commutator) THEN
                  source_rs => gth_current_path_rs(idir:idir, operator_order + 1, :)
                  source_cell_to_index => gth_cell_to_index
                  CALL rskp_directional_taylor_transform(work_real(1)%matrix, work_imag(1)%matrix, source_rs, &
                                                         1, xkp, source_cell_to_index, sab_nl, hmat, direction, &
                                                         lattice_order, shift_scale=-0.5_dp)
                  scale = ecd_periodic_i_power(operator_order)
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
!> \brief Build directional Taylor coefficients of the length-gauge position vertex.
! **************************************************************************************************
   SUBROUTINE ecd_periodic_length_taylor_build(length_real, length_imag, moments, sab_nl, template, direction)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: length_real, length_imag
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                          :: sab_nl
      TYPE(dbcsr_type), POINTER                           :: template
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: direction

      INTEGER                                             :: idir, ix, iy, iz, jdir, moment_index
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: work

      NULLIFY (work)
      IF (SIZE(moments) < ncoset(2) - 1) THEN
         CPABORT("ECD| Length Taylor construction requires coordinate moments through second order")
      END IF
      CALL dbcsr_allocate_matrix_set(length_real, ncart, 2)
      CALL dbcsr_allocate_matrix_set(length_imag, ncart, 2)
      DO idir = 1, ncart
         DO ix = 1, 2
            CALL dbcsr_init_p(length_real(idir, ix)%matrix)
            CALL dbcsr_create(length_real(idir, ix)%matrix, template=template, &
                              matrix_type=dbcsr_type_no_symmetry, name='ECD LENGTH TAYLOR REAL')
            CALL cp_dbcsr_alloc_block_from_nbl(length_real(idir, ix)%matrix, sab_nl)
            CALL dbcsr_set(length_real(idir, ix)%matrix, 0.0_dp)
            CALL dbcsr_init_p(length_imag(idir, ix)%matrix)
            CALL dbcsr_create(length_imag(idir, ix)%matrix, template=template, &
                              matrix_type=dbcsr_type_no_symmetry, name='ECD LENGTH TAYLOR IMAG')
            CALL cp_dbcsr_alloc_block_from_nbl(length_imag(idir, ix)%matrix, sab_nl)
            CALL dbcsr_set(length_imag(idir, ix)%matrix, 0.0_dp)
         END DO
      END DO
      CALL ecd_periodic_ao_set_create(work, 1, template, sab_nl, &
                                      'ECD LENGTH TAYLOR WORK', dbcsr_type_no_symmetry)

      DO idir = 1, ncart
         ix = 0
         iy = 0
         iz = 0
         SELECT CASE (idir)
         CASE (1)
            ix = 1
         CASE (2)
            iy = 1
         CASE (3)
            iz = 1
         END SELECT
         moment_index = coset(ix, iy, iz) - 1
         CALL dbcsr_desymmetrize(moments(moment_index)%matrix, length_real(idir, 1)%matrix)

         DO jdir = 1, ncart
            SELECT CASE (idir)
            CASE (1)
               ix = 1
               iy = 0
               iz = 0
            CASE (2)
               ix = 0
               iy = 1
               iz = 0
            CASE (3)
               ix = 0
               iy = 0
               iz = 1
            END SELECT
            SELECT CASE (jdir)
            CASE (1)
               ix = ix + 1
            CASE (2)
               iy = iy + 1
            CASE (3)
               iz = iz + 1
            END SELECT
            moment_index = coset(ix, iy, iz) - 1
            CALL dbcsr_desymmetrize(moments(moment_index)%matrix, work(1)%matrix)
            CALL dbcsr_add(length_imag(idir, 2)%matrix, work(1)%matrix, 1.0_dp, direction(jdir))
         END DO
      END DO

      CALL dbcsr_deallocate_matrix_set(work)

   END SUBROUTINE ecd_periodic_length_taylor_build

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
