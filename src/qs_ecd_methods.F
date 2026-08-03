!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Electronic circular dichroism post-SCF driver.
! **************************************************************************************************
MODULE qs_ecd_methods
   USE cell_types, ONLY: cell_type, &
                         use_perd_none
   USE commutator_rpnl, ONLY: build_com_mom_nl
   USE cp_files, ONLY: close_file, &
                       open_file
   USE cp_control_types, ONLY: dft_control_type
   USE cp_dbcsr_api, ONLY: dbcsr_add, &
                           dbcsr_create, &
                           dbcsr_init_p, &
                           dbcsr_p_type, &
                           dbcsr_set, &
                           dbcsr_type_antisymmetric, &
                           dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link, ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations, ONLY: cp_dbcsr_sm_fm_multiply, &
                                  dbcsr_allocate_matrix_set, &
                                  dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg, ONLY: cp_fm_trace
   USE cp_fm_struct, ONLY: cp_fm_struct_type
   USE cp_fm_types, ONLY: cp_fm_create, &
                          cp_fm_get_element, &
                          cp_fm_get_info, &
                          cp_fm_get_submatrix, &
                          cp_fm_release, &
                          cp_fm_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_generate_filename, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE input_section_types, ONLY: section_vals_get, &
                                  section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE input_constants, ONLY: no_sf_tddfpt, &
                              tddfpt_dipole_length, &
                              tddfpt_dipole_velocity, &
                              tddfpt_kernel_full, &
                              tddfpt_kernel_none, &
                              use_mom_ref_user
   USE kinds, ONLY: default_path_length, &
                    dp
   USE moments_utils, ONLY: get_reference_point
   USE physcon, ONLY: evolt
   USE qs_ecd_operator_utils, ONLY: ecd_project_operator_to_mo
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_ecd_spectra_xyz, ONLY: ecd_spectra_xyz
   USE qs_kind_types, ONLY: qs_kind_type
   USE qs_mo_types, ONLY: get_mo_set, &
                          mo_set_type
   USE qs_moments, ONLY: build_local_magmom_matrix, &
                         build_local_moment_matrix
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type
   USE qs_operators_ao, ONLY: build_lin_mom_matrix
   USE qs_tddfpt2_types, ONLY: tddfpt_ground_state_mos
   USE particle_types, ONLY: particle_type
   USE scf_control_types, ONLY: scf_control_type
   USE util, ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_methods'
   INTEGER, PARAMETER, PRIVATE          :: ncart = 3
   REAL(KIND=dp), PARAMETER, PRIVATE    :: velocity_degeneracy_threshold = 1.0E-8_dp

   PUBLIC :: ecd_post_scf, ecd_tddfpt

CONTAINS

! **************************************************************************************************
!> \brief Dispatch direct post-SCF ECD according to the cell boundary condition.
! **************************************************************************************************
   SUBROUTINE ecd_post_scf(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_post_scf'

      INTEGER                                            :: handle, poisson_periodic, spinflip, tddfpt_kernel
      LOGICAL                                            :: do_ecd, do_hfx, do_kpoints, do_soc, pseudopotential_commutator, &
                                                            tddfpt_enabled
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: ecd_section, hfx_section, poisson_section, soc_section, &
                                                            tddfpt_section

      CALL timeset(routineN, handle)

      NULLIFY (cell, ecd_section, hfx_section, poisson_section, soc_section, tddfpt_section)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")
      CALL section_vals_val_get(ecd_section, "DO_ECD", l_val=do_ecd)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF
      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CPABORT("ECD| Current ECD implementation does not support HFX")
      END IF

      CALL get_qs_env(qs_env, cell=cell, do_kpoints=do_kpoints)
      CPASSERT(ASSOCIATED(cell))
      IF (ALL(cell%perd == 0)) THEN
         poisson_section => section_vals_get_subs_vals(qs_env%input, "DFT%POISSON")
         CALL section_vals_val_get(poisson_section, "PERIODIC", i_val=poisson_periodic)
         IF (poisson_periodic /= use_perd_none) THEN
            CPABORT("ECD| PERIODIC NONE ECD requires POISSON%PERIODIC NONE")
         END IF
         IF (do_kpoints) THEN
            CPABORT("ECD| PERIODIC NONE ECD requires a Gamma-only SCF layout")
         END IF
      ELSEIF (ALL(cell%perd == 1)) THEN
         CALL section_vals_val_get(ecd_section, "PSEUDOPOTENTIAL_COMMUTATOR", &
                                   l_val=pseudopotential_commutator)
         IF (do_kpoints) THEN
            CPABORT("ECD| PERIODIC XYZ ECD currently requires an implicit Gamma SCF layout")
         END IF
      ELSE
         CPABORT("ECD| Current implementation requires CELL%PERIODIC NONE or XYZ")
      END IF

      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")
      CALL section_vals_get(tddfpt_section, explicit=tddfpt_enabled)
      IF (tddfpt_enabled) THEN
         CALL section_vals_val_get(tddfpt_section, "KERNEL", i_val=tddfpt_kernel)
         CALL section_vals_val_get(tddfpt_section, "SPINFLIP", i_val=spinflip)
         soc_section => section_vals_get_subs_vals(tddfpt_section, "SOC")
         CALL section_vals_get(soc_section, explicit=do_soc)
         IF (do_soc) CPABORT("ECD| TDDFPT ECD does not support TDDFPT SOC")
         IF (spinflip /= no_sf_tddfpt) CPABORT("ECD| TDDFPT ECD does not support spin-flip TDDFPT")
         IF (tddfpt_kernel /= tddfpt_kernel_full .AND. tddfpt_kernel /= tddfpt_kernel_none) THEN
            CPABORT("ECD| TDDFPT ECD supports TDDFPT KERNEL FULL and KERNEL NONE")
         END IF
         IF (ALL(cell%perd == 1) .AND. tddfpt_kernel == tddfpt_kernel_full) THEN
            CPABORT("ECD| PERIODIC XYZ ECD currently supports independent-particle transitions only")
         END IF
         CALL timestop(handle)
         RETURN
      END IF

      IF (ALL(cell%perd == 0)) THEN
         CALL ecd_spectra_none(qs_env)
      ELSE
         CALL ecd_spectra_xyz(qs_env, pseudopotential_commutator)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_post_scf

! **************************************************************************************************
!> \brief Dispatch ECD while TDDFPT response data are available.
! **************************************************************************************************
   SUBROUTINE ecd_tddfpt(qs_env, evects, evals, gs_mos, tddfpt_kernel, mult)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: evects
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), POINTER, INTENT(IN) :: gs_mos
      INTEGER, INTENT(IN)                                :: mult, tddfpt_kernel

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_tddfpt'

      INTEGER                                            :: handle
      LOGICAL                                            :: do_ecd, pseudopotential_commutator
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: ecd_section

      CALL timeset(routineN, handle)

      NULLIFY (cell, ecd_section)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")
      CALL section_vals_val_get(ecd_section, "DO_ECD", l_val=do_ecd)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL get_qs_env(qs_env, cell=cell)
      CPASSERT(ASSOCIATED(cell))
      IF (ALL(cell%perd == 0)) THEN
         IF (tddfpt_kernel == tddfpt_kernel_full) THEN
            CALL ecd_spectra_none(qs_env, evects, evals, gs_mos, tddfpt_kernel, mult)
         ELSE
            CALL ecd_spectra_none(qs_env)
         END IF
      ELSEIF (ALL(cell%perd == 1)) THEN
         IF (tddfpt_kernel == tddfpt_kernel_full) CPABORT("ECD| PERIODIC XYZ TDDFPT ECD is unsupported")
         CALL section_vals_val_get(ecd_section, "PSEUDOPOTENTIAL_COMMUTATOR", &
                                   l_val=pseudopotential_commutator)
         CALL ecd_spectra_xyz(qs_env, pseudopotential_commutator)
      ELSE
         CPABORT("ECD| Current implementation requires CELL%PERIODIC NONE or XYZ")
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_tddfpt

! **************************************************************************************************
!> \brief Assemble PERIODIC NONE ECD sticks from direct or TDDFPT transition data.
!> \param qs_env Quickstep environment
! **************************************************************************************************
   SUBROUTINE ecd_spectra_none(qs_env, evects, evals, gs_mos, tddfpt_kernel, mult)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:, :), OPTIONAL, INTENT(IN) :: evects
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, INTENT(IN)  :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), POINTER, OPTIONAL :: gs_mos
      INTEGER, OPTIONAL, INTENT(IN)                      :: mult, tddfpt_kernel

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_spectra_none'

      INTEGER                                            :: column, e1_unit, ecd_gauge, handle, homo, idir, iocc, iounit, ispin, &
                                                            istate, &
                                                            itrans, ivirt, m1_unit, nao, nactive, nmo, nspins, nstates, &
                                                            ntrans_total, output_index, reference, velocity_e1_unit, &
                                                            velocity_m1_unit, velocity_sticks_unit, row
      CHARACTER(LEN=default_path_length)                 :: e1_filename, m1_filename, velocity_e1_filename, &
                                                            velocity_m1_filename, velocity_sticks_filename
      LOGICAL                                            :: pseudopotential_commutator, use_direct_transitions, use_velocity_gauge
      LOGICAL                                            :: ref_point_explicit, write_transition_matrices
      REAL(KIND=dp), DIMENSION(ncart)                    :: reference_point
      REAL(KIND=dp), DIMENSION(ncart)                    :: transition_magnetic, transition_position
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: transition_occupied, transition_order, transition_spin, &
                                                            transition_virtual
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE     :: magnetic_spin, position_spin
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: transition_magnetic_all, transition_position_all, &
                                                            transition_velocity_all
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: rotatory_strength, rotatory_strength_velocity, &
                                                            transition_energy, transition_weight_all
      COMPLEX(KIND=dp), DIMENSION(:, :), ALLOCATABLE     :: transition_magnetic_velocity_all
      COMPLEX(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE  :: magnetic_velocity_matrix, position_velocity_matrix, &
                                                            velocity_matrix
      COMPLEX(KIND=dp), DIMENSION(ncart)                 :: electric_dipole_velocity
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: velocity_raw_matrix
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: magnetic_active_action, position_active_action
      TYPE(cp_fm_type), DIMENSION(ncart)                 :: magnetic_mo, position_mo, velocity_mo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_ao, position_ao, velocity_ao
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: ecd_section

      CALL timeset(routineN, handle)

      NULLIFY (ecd_section, eigenvalues, fm_struct, logger, magnetic_ao, mo_coeff, mos, occupation_numbers, &
               position_ao, ref_point, scf_control, velocity_ao)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")

      use_direct_transitions = .TRUE.
      IF (PRESENT(evects)) THEN
         IF (.NOT. (PRESENT(evals) .AND. PRESENT(gs_mos) .AND. PRESENT(tddfpt_kernel) .AND. PRESENT(mult))) THEN
            CPABORT("ECD| PERIODIC NONE TDDFPT ECD requires complete response data")
         END IF
         IF (tddfpt_kernel /= tddfpt_kernel_full) THEN
            CPABORT("ECD| PERIODIC NONE response ECD requires TDDFPT KERNEL FULL")
         END IF
         use_direct_transitions = .FALSE.
      END IF

      CALL section_vals_val_get(ecd_section, "GAUGE", i_val=ecd_gauge)
      IF (ecd_gauge /= tddfpt_dipole_length .AND. ecd_gauge /= tddfpt_dipole_velocity) THEN
         CPABORT("ECD| PERIODIC NONE ECD supports GAUGE LENGTH and VELOCITY")
      END IF
      IF (.NOT. use_direct_transitions .AND. ecd_gauge /= tddfpt_dipole_length) THEN
         CPABORT("ECD| TDDFPT ECD currently supports GAUGE LENGTH")
      END IF
      use_velocity_gauge = use_direct_transitions .AND. ecd_gauge == tddfpt_dipole_velocity

      CALL section_vals_val_get(ecd_section, "PSEUDOPOTENTIAL_COMMUTATOR", &
                                l_val=pseudopotential_commutator)

      CALL section_vals_val_get(ecd_section, "REFERENCE", i_val=reference)
      CALL section_vals_val_get(ecd_section, "REF_POINT", explicit=ref_point_explicit)
      IF (ref_point_explicit) THEN
         CALL section_vals_val_get(ecd_section, "REF_POINT", r_vals=ref_point)
      END IF
      IF (reference == use_mom_ref_user) THEN
         IF (.NOT. ref_point_explicit) THEN
            CPABORT("ECD| REFERENCE USER_DEFINED requires REF_POINT")
         END IF
         IF (.NOT. ASSOCIATED(ref_point) .OR. SIZE(ref_point) /= ncart) THEN
            CPABORT("ECD| REF_POINT expects three coordinates")
         END IF
         CALL get_reference_point(reference_point, qs_env=qs_env, reference=reference, ref_point=ref_point)
      ELSE
         IF (ref_point_explicit) THEN
            CPWARN("ECD| REF_POINT is used only with REFERENCE USER_DEFINED")
         END IF
         CALL get_reference_point(reference_point, qs_env=qs_env, reference=reference)
      END IF

      IF (use_direct_transitions) THEN
         ! Direct occupied-virtual transitions are assembled independently.
         CALL get_qs_env(qs_env, mos=mos, scf_control=scf_control)
         CPASSERT(ASSOCIATED(mos))
         CPASSERT(ASSOCIATED(scf_control))
         IF (scf_control%use_ot) THEN
            CPABORT("ECD| IPA currently requires diagonalization SCF with available virtual MOs")
         END IF

         nspins = SIZE(mos)
         ntrans_total = 0
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), homo=homo, nmo=nmo)
            IF (homo > 0 .AND. nmo > homo) ntrans_total = ntrans_total + homo*(nmo - homo)
         END DO
         IF (ntrans_total <= 0) CPABORT("ECD| No independent-particle transitions available")

         ALLOCATE (transition_magnetic_all(ncart, ntrans_total), transition_position_all(ncart, ntrans_total), &
                   transition_occupied(ntrans_total), transition_order(ntrans_total), transition_spin(ntrans_total), &
                   transition_energy(ntrans_total), transition_virtual(ntrans_total), rotatory_strength(ntrans_total), &
                   transition_weight_all(ntrans_total))
         IF (use_velocity_gauge) THEN
            ALLOCATE (transition_magnetic_velocity_all(ncart, ntrans_total), &
                      transition_velocity_all(ncart, ntrans_total), rotatory_strength_velocity(ntrans_total))
         END IF
         logger => cp_get_default_logger()
         iounit = cp_logger_get_default_io_unit(logger)
         CALL cp_logger_generate_filename(logger, e1_filename, logger%iter_info%project_name, &
                                          "-ecd-e1-transition-matrix.dat", local=.FALSE.)
         CALL cp_logger_generate_filename(logger, m1_filename, logger%iter_info%project_name, &
                                          "-ecd-m1-transition-matrix.dat", local=.FALSE.)
         IF (use_velocity_gauge) THEN
            CALL cp_logger_generate_filename(logger, velocity_e1_filename, logger%iter_info%project_name, &
                                             "-ecd-velocity-e1-transition-matrix.dat", local=.FALSE.)
            CALL cp_logger_generate_filename(logger, velocity_m1_filename, logger%iter_info%project_name, &
                                             "-ecd-velocity-m1-transition-matrix.dat", local=.FALSE.)
            CALL cp_logger_generate_filename(logger, velocity_sticks_filename, logger%iter_info%project_name, &
                                             "-ecd-velocity-ipa-rotatory-strengths.dat", local=.FALSE.)
         END IF
         write_transition_matrices = logger%para_env%is_source()
         IF (write_transition_matrices) THEN
            CALL open_file(e1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=e1_unit)
            WRITE (e1_unit, '(A)') "# PERIODIC NONE length-gauge E1 transition matrix"
            WRITE (e1_unit, '(A)') "# branch is zero for direct molecular transitions"
            WRITE (e1_unit, '(A)') "# current = -gap*r; dipole = -r"
            WRITE (e1_unit, '(A)') "# transition energy_eV spin occupied virtual branch component "// &
                                   "current_real current_imag dipole_real dipole_imag occupation_difference"
            CALL open_file(m1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=m1_unit)
            WRITE (m1_unit, '(A)') "# PERIODIC NONE M1 transition angular and magnetic-dipole vertices"
            WRITE (m1_unit, '(A)') "# lambda is the complete local plus GTH angular matrix; magnetic = i lambda/2"
            WRITE (m1_unit, '(A)') "# transition energy_eV spin occupied virtual component "// &
                                   "lambda_real lambda_imag magnetic_real magnetic_imag occupation_difference"
            IF (use_velocity_gauge) THEN
               CALL open_file(velocity_e1_filename, file_status="UNKNOWN", file_action="WRITE", &
                               unit_number=velocity_e1_unit)
               WRITE (velocity_e1_unit, '(A)') "# PERIODIC NONE E1 transition matrix in the velocity representation"
               WRITE (velocity_e1_unit, '(A)') "# velocity = -i*(nabla + [r,V_nl]); velocity dipole = i*velocity/gap"
               WRITE (velocity_e1_unit, '(A)') "# transition energy_eV spin occupied virtual branch component "// &
                                               "velocity_real velocity_imag dipole_real dipole_imag occupation_difference"
               CALL open_file(velocity_m1_filename, file_status="UNKNOWN", file_action="WRITE", &
                              unit_number=velocity_m1_unit)
               WRITE (velocity_m1_unit, '(A)') "# PERIODIC NONE M1 transition matrix in the velocity representation"
               WRITE (velocity_m1_unit, '(A)') "# angular = r_velocity cross velocity; magnetic = -angular/2"
               WRITE (velocity_m1_unit, '(A)') "# transition energy_eV spin occupied virtual component "// &
                                               "angular_real angular_imag magnetic_real magnetic_imag occupation_difference"
               CALL open_file(velocity_sticks_filename, file_status="UNKNOWN", file_action="WRITE", &
                              unit_number=velocity_sticks_unit)
               WRITE (velocity_sticks_unit, '(A)') &
                  "# PERIODIC NONE IPA direct length and complete velocity-gauge rotatory strengths"
               WRITE (velocity_sticks_unit, '(A)') "# state energy_eV spin occupied virtual occupation_difference "// &
                                                   "rotatory_strength_length rotatory_strength_velocity"
            END IF
         END IF
         itrans = 0
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, homo=homo, mo_coeff=mo_coeff, nao=nao, &
                            nmo=nmo, occupation_numbers=occupation_numbers)
            CPASSERT(ASSOCIATED(eigenvalues))
            CPASSERT(ASSOCIATED(mo_coeff))
            CPASSERT(ASSOCIATED(occupation_numbers))
            IF (homo <= 0 .OR. nmo <= homo) CYCLE

            IF (use_velocity_gauge) THEN
               CALL ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, reference_point, &
                                                          position_ao, magnetic_ao, velocity_ao)
            ELSE
               CALL ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, reference_point, &
                                                          position_ao, magnetic_ao)
            END IF
            DO idir = 1, ncart
               CALL ecd_project_operator_to_mo(position_ao(idir), mo_coeff, position_mo(idir))
               CALL ecd_project_operator_to_mo(magnetic_ao(idir), mo_coeff, magnetic_mo(idir))
               IF (use_velocity_gauge) THEN
                  CALL ecd_project_operator_to_mo(velocity_ao(idir), mo_coeff, velocity_mo(idir))
               END IF
            END DO
            IF (use_velocity_gauge) THEN
               ! Convert raw D=nabla+[r,V_nl] to physical v=-iD, then set
               ! r^V_pq = v_pq/[i(e_p-e_q)] in the non-degenerate MO subspace.
               ! The diagonal is fixed to zero; its orbital-gauge contribution drops out of
               ! Im[(mu^V_ai)^* dot m^V_ai] for the real-orbital IPA implementation.
               ALLOCATE (magnetic_velocity_matrix(ncart, nmo, nmo), &
                         position_velocity_matrix(ncart, nmo, nmo), velocity_matrix(ncart, nmo, nmo), &
                         velocity_raw_matrix(nmo, nmo))
               DO idir = 1, ncart
                  CALL cp_fm_get_submatrix(velocity_mo(idir), velocity_raw_matrix)
                  velocity_matrix(idir, :, :) = CMPLX(0.0_dp, -velocity_raw_matrix, KIND=dp)
               END DO
               DEALLOCATE (velocity_raw_matrix)
               position_velocity_matrix = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
               DO column = 1, nmo
                  DO row = 1, nmo
                     IF (row == column) CYCLE
                     IF (ABS(eigenvalues(row) - eigenvalues(column)) <= velocity_degeneracy_threshold) THEN
                        CPABORT("ECD| VELOCITY-gauge M1 requires non-degenerate molecular orbital energies")
                     END IF
                     position_velocity_matrix(:, row, column) = velocity_matrix(:, row, column)/ &
                                                                 CMPLX(0.0_dp, eigenvalues(row) - eigenvalues(column), KIND=dp)
                  END DO
               END DO
               magnetic_velocity_matrix(1, :, :) = -0.5_dp*( &
                  MATMUL(position_velocity_matrix(2, :, :), velocity_matrix(3, :, :)) - &
                  MATMUL(position_velocity_matrix(3, :, :), velocity_matrix(2, :, :)))
               magnetic_velocity_matrix(2, :, :) = -0.5_dp*( &
                  MATMUL(position_velocity_matrix(3, :, :), velocity_matrix(1, :, :)) - &
                  MATMUL(position_velocity_matrix(1, :, :), velocity_matrix(3, :, :)))
               magnetic_velocity_matrix(3, :, :) = -0.5_dp*( &
                  MATMUL(position_velocity_matrix(1, :, :), velocity_matrix(2, :, :)) - &
                  MATMUL(position_velocity_matrix(2, :, :), velocity_matrix(1, :, :)))
            END IF
            DO iocc = 1, homo
               DO ivirt = homo + 1, nmo
                  itrans = itrans + 1
                  transition_energy(itrans) = eigenvalues(ivirt) - eigenvalues(iocc)
                  IF (transition_energy(itrans) <= 0.0_dp) CPABORT("ECD| IPA requires positive occupied-virtual gaps")
                  transition_weight_all(itrans) = occupation_numbers(iocc) - occupation_numbers(ivirt)
                  IF (transition_weight_all(itrans) <= 0.0_dp) THEN
                     CPABORT("ECD| IPA requires a positive occupied-virtual occupation difference")
                  END IF
                  transition_occupied(itrans) = iocc
                  transition_spin(itrans) = ispin
                  transition_virtual(itrans) = ivirt
                  DO idir = 1, ncart
                     CALL cp_fm_get_element(position_mo(idir), ivirt, iocc, &
                                             transition_position_all(idir, itrans))
                     CALL cp_fm_get_element(magnetic_mo(idir), ivirt, iocc, &
                                             transition_magnetic_all(idir, itrans))
                     IF (use_velocity_gauge) THEN
                        CALL cp_fm_get_element(velocity_mo(idir), ivirt, iocc, &
                                                transition_velocity_all(idir, itrans))
                        transition_magnetic_velocity_all(idir, itrans) = &
                           magnetic_velocity_matrix(idir, ivirt, iocc)
                     END IF
                  END DO
                  ! For real Gamma-point orbitals without SOC,
                  ! R_ia = -w_ia*sum_lambda(r_lambda,ai*Lambda_lambda,ai)/2.
                  rotatory_strength(itrans) = -0.5_dp*transition_weight_all(itrans)* &
                                               SUM(transition_position_all(:, itrans)*transition_magnetic_all(:, itrans))
                  IF (use_velocity_gauge) THEN
                     ! mu^V_ci = i*v_ci/E_ci and m^V_ci = -[r^V x v]_ci/2.
                     electric_dipole_velocity = CMPLX(0.0_dp, 1.0_dp, KIND=dp)* &
                                                 velocity_matrix(:, ivirt, iocc)/transition_energy(itrans)
                     rotatory_strength_velocity(itrans) = transition_weight_all(itrans)*AIMAG( &
                        SUM(CONJG(electric_dipole_velocity)*transition_magnetic_velocity_all(:, itrans)))
                  END IF
               END DO
            END DO
            DO idir = 1, ncart
               CALL cp_fm_release(magnetic_mo(idir))
               CALL cp_fm_release(position_mo(idir))
               IF (use_velocity_gauge) CALL cp_fm_release(velocity_mo(idir))
            END DO
            IF (use_velocity_gauge) THEN
               DEALLOCATE (magnetic_velocity_matrix, position_velocity_matrix, velocity_matrix)
            END IF
            CALL dbcsr_deallocate_matrix_set(magnetic_ao)
            CALL dbcsr_deallocate_matrix_set(position_ao)
            IF (use_velocity_gauge) CALL dbcsr_deallocate_matrix_set(velocity_ao)
         END DO
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T2,A)') "ECD| PERIODIC NONE E1 transition matrix written to "//TRIM(e1_filename)
            WRITE (iounit, '(T2,A)') "ECD| PERIODIC NONE M1 transition matrix written to "//TRIM(m1_filename)
            IF (use_velocity_gauge) THEN
               WRITE (iounit, '(T2,A)') "ECD| PERIODIC NONE velocity E1 matrix written to "//TRIM(velocity_e1_filename)
               WRITE (iounit, '(T2,A)') "ECD| PERIODIC NONE velocity M1 matrix written to "//TRIM(velocity_m1_filename)
               WRITE (iounit, '(T2,A)') "ECD| PERIODIC NONE length/velocity rotatory strengths written to "// &
                                         TRIM(velocity_sticks_filename)
            END IF
         END IF
         CALL sort(transition_energy, ntrans_total, transition_order)
         IF (write_transition_matrices) THEN
            WRITE (e1_unit, '(A)') "# records are ordered by increasing excitation energy"
            WRITE (m1_unit, '(A)') "# records are ordered by increasing excitation energy"
            IF (use_velocity_gauge) WRITE (velocity_e1_unit, '(A)') "# records are ordered by increasing excitation energy"
            IF (use_velocity_gauge) WRITE (velocity_m1_unit, '(A)') "# records are ordered by increasing excitation energy"
            DO output_index = 1, ntrans_total
               itrans = transition_order(output_index)
               DO idir = 1, ncart
                  WRITE (e1_unit, '(I0,1X,ES24.15E3,5(1X,I0),5(1X,ES24.15E3))') &
                     output_index, transition_energy(output_index)*evolt, transition_spin(itrans), &
                     transition_occupied(itrans), transition_virtual(itrans), 0, idir, &
                     -transition_energy(output_index)*transition_position_all(idir, itrans), 0.0_dp, &
                     -transition_position_all(idir, itrans), 0.0_dp, &
                     transition_weight_all(itrans)
                  WRITE (m1_unit, '(I0,1X,ES24.15E3,4(1X,I0),5(1X,ES24.15E3))') &
                     output_index, transition_energy(output_index)*evolt, transition_spin(itrans), &
                     transition_occupied(itrans), transition_virtual(itrans), idir, &
                     transition_magnetic_all(idir, itrans), 0.0_dp, 0.0_dp, &
                     0.5_dp*transition_magnetic_all(idir, itrans), &
                     transition_weight_all(itrans)
                  IF (use_velocity_gauge) THEN
                     WRITE (velocity_e1_unit, '(I0,1X,ES24.15E3,5(1X,I0),5(1X,ES24.15E3))') &
                        output_index, transition_energy(output_index)*evolt, transition_spin(itrans), &
                        transition_occupied(itrans), transition_virtual(itrans), 0, idir, &
                        0.0_dp, -transition_velocity_all(idir, itrans), &
                        transition_velocity_all(idir, itrans)/transition_energy(itrans), 0.0_dp, &
                        transition_weight_all(itrans)
                     WRITE (velocity_m1_unit, '(I0,1X,ES24.15E3,4(1X,I0),5(1X,ES24.15E3))') &
                        output_index, transition_energy(output_index)*evolt, transition_spin(itrans), &
                        transition_occupied(itrans), transition_virtual(itrans), idir, &
                        REAL(-2.0_dp*transition_magnetic_velocity_all(idir, itrans), KIND=dp), &
                        AIMAG(-2.0_dp*transition_magnetic_velocity_all(idir, itrans)), &
                        REAL(transition_magnetic_velocity_all(idir, itrans), KIND=dp), &
                        AIMAG(transition_magnetic_velocity_all(idir, itrans)), transition_weight_all(itrans)
                  END IF
               END DO
               IF (use_velocity_gauge) THEN
                  WRITE (velocity_sticks_unit, '(I0,1X,ES24.15E3,3(1X,I0),3(1X,ES24.15E3))') &
                     output_index, transition_energy(output_index)*evolt, transition_spin(itrans), &
                     transition_occupied(itrans), transition_virtual(itrans), transition_weight_all(itrans), &
                     rotatory_strength(itrans), rotatory_strength_velocity(itrans)
               END IF
            END DO
            CALL close_file(e1_unit)
            CALL close_file(m1_unit)
            IF (use_velocity_gauge) THEN
               CALL close_file(velocity_e1_unit)
               CALL close_file(velocity_m1_unit)
               CALL close_file(velocity_sticks_unit)
            END IF
         END IF
         IF (use_velocity_gauge) THEN
            CALL ecd_print_sticks(transition_energy, rotatory_strength_velocity, transition_order)
         ELSE
            CALL ecd_print_sticks(transition_energy, rotatory_strength, transition_order)
         END IF
         DEALLOCATE (transition_magnetic_all, transition_occupied, transition_order, transition_position_all, &
                     transition_spin, transition_energy, transition_virtual, rotatory_strength, transition_weight_all)
         IF (use_velocity_gauge) THEN
            DEALLOCATE (rotatory_strength_velocity, transition_magnetic_velocity_all, transition_velocity_all)
         END IF
      ELSE
         ! TDDFPT KERNEL FULL contracts the coupled response amplitudes with the operator actions.
         nspins = SIZE(gs_mos)
         nstates = SIZE(evals)
         IF (SIZE(evects, 1) /= nspins .OR. SIZE(evects, 2) /= nstates) THEN
            CPABORT("ECD| TDDFPT eigenvector dimensions are inconsistent with TDDFPT response data")
         END IF
         IF (nspins == 1 .AND. mult /= 1) THEN
            CPABORT("ECD| TDDFPT ECD supports restricted singlet TDDFPT states")
         END IF
         IF (nspins > 2) CPABORT("ECD| TDDFPT ECD supports at most two spin channels")

         ALLOCATE (position_active_action(ncart, nspins), magnetic_active_action(ncart, nspins))
         CALL ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, reference_point, &
                                                    position_ao, magnetic_ao)
         DO ispin = 1, nspins
            CPASSERT(ASSOCIATED(gs_mos(ispin)%mos_active))
            NULLIFY (fm_struct)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_active, matrix_struct=fm_struct, ncol_global=nactive)
            DO idir = 1, ncart
               CALL cp_fm_create(position_active_action(idir, ispin), fm_struct)
               CALL cp_fm_create(magnetic_active_action(idir, ispin), fm_struct)
               CALL cp_dbcsr_sm_fm_multiply(position_ao(idir)%matrix, gs_mos(ispin)%mos_active, &
                                             position_active_action(idir, ispin), ncol=nactive, alpha=1.0_dp, beta=0.0_dp)
               CALL cp_dbcsr_sm_fm_multiply(magnetic_ao(idir)%matrix, gs_mos(ispin)%mos_active, &
                                             magnetic_active_action(idir, ispin), ncol=nactive, alpha=1.0_dp, beta=0.0_dp)
            END DO
         END DO
         NULLIFY (fm_struct)
         CALL dbcsr_deallocate_matrix_set(magnetic_ao)
         CALL dbcsr_deallocate_matrix_set(position_ao)
         ALLOCATE (position_spin(nstates, ncart, nspins), magnetic_spin(nstates, ncart, nspins))
         position_spin = 0.0_dp
         magnetic_spin = 0.0_dp
         DO ispin = 1, nspins
            CALL cp_fm_get_info(evects(ispin, 1), ncol_global=nactive)
            IF (nactive /= gs_mos(ispin)%nmo_active) THEN
               CPABORT("ECD| TDDFPT response amplitudes and active occupied orbitals are inconsistent")
            END IF
            DO idir = 1, ncart
               CALL cp_fm_trace(evects(ispin, :), position_active_action(idir, ispin), &
                                position_spin(:, idir, ispin))
               CALL cp_fm_trace(evects(ispin, :), magnetic_active_action(idir, ispin), &
                                magnetic_spin(:, idir, ispin))
            END DO
         END DO

         ALLOCATE (rotatory_strength(nstates))
         DO istate = 1, nstates
            IF (nspins == 1) THEN
               transition_position = SQRT(2.0_dp)*position_spin(istate, :, 1)
               transition_magnetic = SQRT(2.0_dp)*magnetic_spin(istate, :, 1)
            ELSE
               transition_position = position_spin(istate, :, 1) + position_spin(istate, :, 2)
               transition_magnetic = magnetic_spin(istate, :, 1) + magnetic_spin(istate, :, 2)
            END IF
            ! TDDFPT amplitudes already provide the state transition moments, so w=1.
            ! The same real-orbital reduction gives R = -sum_lambda(r_lambda*Lambda_lambda)/2.
            rotatory_strength(istate) = -0.5_dp*SUM(transition_position*transition_magnetic)
         END DO
         CALL ecd_print_sticks(evals, rotatory_strength)

         DO ispin = 1, nspins
            DO idir = 1, ncart
               CALL cp_fm_release(position_active_action(idir, ispin))
               CALL cp_fm_release(magnetic_active_action(idir, ispin))
            END DO
         END DO
         DEALLOCATE (magnetic_active_action, position_active_action, magnetic_spin, position_spin, rotatory_strength)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_spectra_none

! **************************************************************************************************
!> \brief Print molecular ECD excitation energies and rotatory strengths.
! **************************************************************************************************
   SUBROUTINE ecd_print_sticks(excitation_energy, rotatory_strength, order)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: excitation_energy, rotatory_strength
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)        :: order

      INTEGER                                            :: istate, iounit, source_state
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (SIZE(excitation_energy) /= SIZE(rotatory_strength)) THEN
         CPABORT("ECD| Excitation-energy and rotatory-strength dimensions are inconsistent")
      END IF
      IF (PRESENT(order)) THEN
         IF (SIZE(order) /= SIZE(excitation_energy)) THEN
            CPABORT("ECD| ECD stick ordering dimensions are inconsistent")
         END IF
      END IF

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      IF (iounit <= 0) RETURN

      WRITE (iounit, '(/,T8,A,T20,A,T43,A)') "State", "Excitation energy", "Rotatory strength"
      WRITE (iounit, '(T8,A,T20,A,T43,A)') "number", "eV", "au"
      WRITE (iounit, '(T8,62("-"))')

      DO istate = 1, SIZE(excitation_energy)
         source_state = istate
         IF (PRESENT(order)) source_state = order(istate)
         WRITE (iounit, '(1X,A,T11,I7,T20,F11.5,T43,ES12.5E2)') &
            "ECD|", istate, excitation_energy(istate)*evolt, rotatory_strength(source_state)
      END DO
   END SUBROUTINE ecd_print_sticks

! **************************************************************************************************
!> \brief Build direct length, complete velocity, and orbital magnetic AO operators for PERIODIC NONE.
!> \param qs_env Quickstep environment
!> \param pseudopotential_commutator Include the complete GTH non-local velocity and magnetic corrections
!> \param reference_point Common reference point for all molecular operators
!> \param position_ao Three direct AO-space position matrices
!> \param magnetic_ao Three direct AO-space angular-momentum matrices
!> \param velocity_ao Optional three real antisymmetric AO velocity matrices
! **************************************************************************************************
   SUBROUTINE ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, reference_point, &
                                                    position_ao, magnetic_ao, velocity_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: reference_point
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER, &
         INTENT(INOUT)                                   :: magnetic_ao, position_ao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER, &
         OPTIONAL, INTENT(INOUT)                          :: velocity_ao

      INTEGER                                            :: idir
      LOGICAL                                            :: build_velocity
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_nl, matrix_s, velocity_nl
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (cell, dft_control, magnetic_nl, matrix_s, particle_set, qs_kind_set, sab_orb, sap_ppnl, velocity_nl)
      build_velocity = PRESENT(velocity_ao)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_s=matrix_s, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set, sab_orb=sab_orb, &
                      sap_ppnl=sap_ppnl)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(sab_orb))

      CALL dbcsr_allocate_matrix_set(position_ao, ncart)
      CALL dbcsr_allocate_matrix_set(magnetic_ao, ncart)
      IF (build_velocity) CALL dbcsr_allocate_matrix_set(velocity_ao, ncart)
      DO idir = 1, ncart
         CALL dbcsr_init_p(position_ao(idir)%matrix)
         CALL dbcsr_create(position_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(position_ao(idir)%matrix, sab_orb)
         CALL dbcsr_set(position_ao(idir)%matrix, 0.0_dp)

         CALL dbcsr_init_p(magnetic_ao(idir)%matrix)
         CALL dbcsr_create(magnetic_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(magnetic_ao(idir)%matrix, sab_orb)
         CALL dbcsr_set(magnetic_ao(idir)%matrix, 0.0_dp)

         IF (build_velocity) THEN
            CALL dbcsr_init_p(velocity_ao(idir)%matrix)
            CALL dbcsr_create(velocity_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(velocity_ao(idir)%matrix, sab_orb)
            CALL dbcsr_set(velocity_ao(idir)%matrix, 0.0_dp)
         END IF
      END DO

      CALL build_local_moment_matrix(qs_env, position_ao, 1, ref_point=reference_point)
      CALL build_local_magmom_matrix(qs_env, magnetic_ao, 1, ref_point=reference_point)
      IF (build_velocity) CALL build_lin_mom_matrix(qs_env, velocity_ao)

      IF (pseudopotential_commutator .AND. ASSOCIATED(sap_ppnl)) THEN
         CALL dbcsr_allocate_matrix_set(magnetic_nl, ncart)
         IF (build_velocity) CALL dbcsr_allocate_matrix_set(velocity_nl, ncart)
         DO idir = 1, ncart
            CALL dbcsr_init_p(magnetic_nl(idir)%matrix)
            CALL dbcsr_create(magnetic_nl(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(magnetic_nl(idir)%matrix, sab_orb)
            CALL dbcsr_set(magnetic_nl(idir)%matrix, 0.0_dp)

            IF (build_velocity) THEN
               CALL dbcsr_init_p(velocity_nl(idir)%matrix)
               CALL dbcsr_create(velocity_nl(idir)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(velocity_nl(idir)%matrix, sab_orb)
               CALL dbcsr_set(velocity_nl(idir)%matrix, 0.0_dp)
            END IF
         END DO
         ! r_R0 x [r_R0, V_nl] is the non-local contribution required by
         ! the velocity-defined orbital magnetic dipole.
         IF (build_velocity) THEN
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                                  particle_set, cell, matrix_rv=velocity_nl, matrix_rxrv=magnetic_nl, &
                                  ref_point=reference_point)
         ELSE
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                                  particle_set, cell, matrix_rxrv=magnetic_nl, ref_point=reference_point)
         END IF
         DO idir = 1, ncart
            CALL dbcsr_add(magnetic_ao(idir)%matrix, magnetic_nl(idir)%matrix, 1.0_dp, 1.0_dp)
            IF (build_velocity) THEN
               CALL dbcsr_add(velocity_ao(idir)%matrix, velocity_nl(idir)%matrix, 1.0_dp, 1.0_dp)
            END IF
         END DO
         CALL dbcsr_deallocate_matrix_set(magnetic_nl)
         IF (build_velocity) CALL dbcsr_deallocate_matrix_set(velocity_nl)
      END IF

   END SUBROUTINE ecd_build_molecular_operators_ao_none

END MODULE qs_ecd_methods
