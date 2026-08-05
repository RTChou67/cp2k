!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Electronic circular dichroism input validation and spectrum drivers.
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
   USE qs_ecd_spectra_xyz, ONLY: ecd_spectra_xyz_ipa
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

   PUBLIC :: ecd_dispatch, ecd_spectra_none_tddfpt, ecd_validate_input

CONTAINS

! **************************************************************************************************
!> \brief Validate ECD input before any optional TDDFPT calculation.
! **************************************************************************************************
   SUBROUTINE ecd_validate_input(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_validate_input'

      INTEGER                                            :: handle, poisson_periodic, spinflip, tddfpt_kernel
      LOGICAL                                            :: do_ecd, do_hfx, do_kpoints, do_soc, &
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
         IF (tddfpt_kernel == tddfpt_kernel_full) THEN
            IF (ALL(cell%perd == 1)) THEN
               CPABORT("ECD| PERIODIC XYZ ECD currently supports independent-particle transitions only")
            END IF
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_validate_input

! **************************************************************************************************
!> \brief Dispatch independent-particle ECD after optional TDDFPT has completed.
! **************************************************************************************************
   SUBROUTINE ecd_dispatch(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_dispatch'

      INTEGER                                            :: handle, tddfpt_kernel
      LOGICAL                                            :: do_ecd, pseudopotential_commutator, tddfpt_enabled
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: ecd_section, tddfpt_section

      CALL timeset(routineN, handle)

      NULLIFY (cell, ecd_section, tddfpt_section)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")
      CALL section_vals_val_get(ecd_section, "DO_ECD", l_val=do_ecd)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF

      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")
      CALL section_vals_get(tddfpt_section, explicit=tddfpt_enabled)
      IF (tddfpt_enabled) THEN
         CALL section_vals_val_get(tddfpt_section, "KERNEL", i_val=tddfpt_kernel)
         IF (tddfpt_kernel == tddfpt_kernel_full) THEN
            CALL timestop(handle)
            RETURN
         END IF
      END IF

      CALL get_qs_env(qs_env, cell=cell)
      CPASSERT(ASSOCIATED(cell))
      IF (ALL(cell%perd == 0)) THEN
         CALL ecd_spectra_none_ipa(qs_env)
      ELSEIF (ALL(cell%perd == 1)) THEN
         CALL section_vals_val_get(ecd_section, "PSEUDOPOTENTIAL_COMMUTATOR", &
                                   l_val=pseudopotential_commutator)
         CALL ecd_spectra_xyz_ipa(qs_env, pseudopotential_commutator)
      ELSE
         CPABORT("ECD| Current implementation requires CELL%PERIODIC NONE or XYZ")
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_dispatch

! **************************************************************************************************
!> \brief Assemble PERIODIC NONE IPA ECD sticks from occupied--virtual transitions.
!> \param qs_env Quickstep environment
! **************************************************************************************************
   SUBROUTINE ecd_spectra_none_ipa(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_spectra_none_ipa'

      INTEGER                                            :: column, e1_unit, ecd_gauge, handle, homo, idir, iocc, iounit, ispin, &
                                                            itrans, ivirt, m1_unit, nmo, nspins, &
                                                            ntrans_total, reference, row
      CHARACTER(LEN=default_path_length)                 :: e1_filename, m1_filename
      LOGICAL                                            :: pseudopotential_commutator, use_velocity_gauge
      LOGICAL                                            :: ref_point_explicit, write_transition_matrices
      REAL(KIND=dp), DIMENSION(ncart)                    :: reference_point
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: transition_order
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: rotatory_strength, transition_energy
      COMPLEX(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE  :: e1_mo, m1_mo, velocity_mo_matrix
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: operator_mo_matrix
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_fm_type)                                   :: operator_mo
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_ao, position_ao, velocity_ao
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: ecd_section

      CALL timeset(routineN, handle)

      NULLIFY (ecd_section, eigenvalues, logger, magnetic_ao, mo_coeff, mos, occupation_numbers, &
               position_ao, ref_point, scf_control, velocity_ao)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")

      CALL section_vals_val_get(ecd_section, "GAUGE", i_val=ecd_gauge)
      IF (ecd_gauge /= tddfpt_dipole_length .AND. ecd_gauge /= tddfpt_dipole_velocity) THEN
         CPABORT("ECD| PERIODIC NONE ECD supports GAUGE LENGTH and VELOCITY")
      END IF
      use_velocity_gauge = ecd_gauge == tddfpt_dipole_velocity

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

      ALLOCATE (transition_order(ntrans_total), transition_energy(ntrans_total), rotatory_strength(ntrans_total))
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      CALL cp_logger_generate_filename(logger, e1_filename, logger%iter_info%project_name, &
                                       "-ecd-e1-transition-matrix.dat", local=.FALSE.)
      CALL cp_logger_generate_filename(logger, m1_filename, logger%iter_info%project_name, &
                                       "-ecd-m1-transition-matrix.dat", local=.FALSE.)
      write_transition_matrices = logger%para_env%is_source()
      IF (write_transition_matrices) THEN
         CALL open_file(e1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=e1_unit)
         WRITE (e1_unit, '(A)') "# PERIODIC NONE IPA E1 transition matrix"
         IF (use_velocity_gauge) THEN
            WRITE (e1_unit, '(A)') "# mu_ai = i*v_ai/(epsilon_a-epsilon_i)"
         ELSE
            WRITE (e1_unit, '(A)') "# mu_ai = -r_ai"
         END IF
         WRITE (e1_unit, '(A)') "# records follow the occupied-empty MO traversal order"
         WRITE (e1_unit, '(A)') "# state energy_eV spin occupied virtual component e1_real e1_imag occupation_difference"
         CALL open_file(m1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=m1_unit)
         WRITE (m1_unit, '(A)') "# PERIODIC NONE IPA M1 transition matrix"
         IF (use_velocity_gauge) THEN
            WRITE (m1_unit, '(A)') "# m_ai = -[r_velocity cross v]_ai/2"
         ELSE
            WRITE (m1_unit, '(A)') "# m_ai = i*Lambda_ai/2"
         END IF
         WRITE (m1_unit, '(A)') "# records follow the occupied-empty MO traversal order"
         WRITE (m1_unit, '(A)') "# state energy_eV spin occupied virtual component m1_real m1_imag occupation_difference"
      END IF
      IF (use_velocity_gauge) THEN
         CALL ecd_build_molecular_velocity_ao_none(qs_env, pseudopotential_commutator, velocity_ao)
      ELSE
         CALL ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, reference_point, &
                                                    position_ao, magnetic_ao)
      END IF

      itrans = 0
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, homo=homo, mo_coeff=mo_coeff, &
                         nmo=nmo, occupation_numbers=occupation_numbers)
         CPASSERT(ASSOCIATED(eigenvalues))
         CPASSERT(ASSOCIATED(mo_coeff))
         CPASSERT(ASSOCIATED(occupation_numbers))
         IF (homo <= 0 .OR. nmo <= homo) CYCLE

         IF (use_velocity_gauge) THEN
            ! Build the physical E1 and M1 MO matrices before selecting
            ! occupied-empty elements.  v is only an intermediate matrix.
            ALLOCATE (e1_mo(ncart, nmo, nmo), m1_mo(ncart, nmo, nmo), &
                      velocity_mo_matrix(ncart, nmo, nmo), operator_mo_matrix(nmo, nmo))
            DO idir = 1, ncart
               CALL ecd_project_operator_to_mo(velocity_ao(idir), mo_coeff, operator_mo)
               CALL cp_fm_get_submatrix(operator_mo, operator_mo_matrix)
               CALL cp_fm_release(operator_mo)
               velocity_mo_matrix(idir, :, :) = CMPLX(0.0_dp, -operator_mo_matrix, KIND=dp)
            END DO
            DEALLOCATE (operator_mo_matrix)
            e1_mo = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            DO column = 1, nmo
               DO row = 1, nmo
                  IF (row == column) CYCLE
                  IF (ABS(eigenvalues(row) - eigenvalues(column)) <= velocity_degeneracy_threshold) THEN
                     CPABORT("ECD| VELOCITY-gauge M1 requires non-degenerate molecular orbital energies")
                  END IF
                  e1_mo(:, row, column) = CMPLX(0.0_dp, 1.0_dp, KIND=dp)*velocity_mo_matrix(:, row, column)/ &
                                          (eigenvalues(row) - eigenvalues(column))
               END DO
            END DO
            ! m^V = -r^V x v/2 = E1^V x v/2 because E1^V = -r^V.
            m1_mo(1, :, :) = 0.5_dp*(MATMUL(e1_mo(2, :, :), velocity_mo_matrix(3, :, :)) - &
                                     MATMUL(e1_mo(3, :, :), velocity_mo_matrix(2, :, :)))
            m1_mo(2, :, :) = 0.5_dp*(MATMUL(e1_mo(3, :, :), velocity_mo_matrix(1, :, :)) - &
                                     MATMUL(e1_mo(1, :, :), velocity_mo_matrix(3, :, :)))
            m1_mo(3, :, :) = 0.5_dp*(MATMUL(e1_mo(1, :, :), velocity_mo_matrix(2, :, :)) - &
                                     MATMUL(e1_mo(2, :, :), velocity_mo_matrix(1, :, :)))
         ELSE
            ALLOCATE (e1_mo(ncart, nmo, nmo), m1_mo(ncart, nmo, nmo), operator_mo_matrix(nmo, nmo))
            DO idir = 1, ncart
               CALL ecd_project_operator_to_mo(position_ao(idir), mo_coeff, operator_mo)
               CALL cp_fm_get_submatrix(operator_mo, operator_mo_matrix)
               CALL cp_fm_release(operator_mo)
               e1_mo(idir, :, :) = CMPLX(-operator_mo_matrix, 0.0_dp, KIND=dp)
               CALL ecd_project_operator_to_mo(magnetic_ao(idir), mo_coeff, operator_mo)
               CALL cp_fm_get_submatrix(operator_mo, operator_mo_matrix)
               CALL cp_fm_release(operator_mo)
               m1_mo(idir, :, :) = CMPLX(0.0_dp, 0.5_dp*operator_mo_matrix, KIND=dp)
            END DO
            DEALLOCATE (operator_mo_matrix)
         END IF
         DO iocc = 1, homo
            DO ivirt = homo + 1, nmo
               itrans = itrans + 1
               transition_energy(itrans) = eigenvalues(ivirt) - eigenvalues(iocc)
               IF (transition_energy(itrans) <= 0.0_dp) CPABORT("ECD| IPA requires positive occupied-virtual gaps")
               IF (occupation_numbers(iocc) - occupation_numbers(ivirt) <= 0.0_dp) THEN
                  CPABORT("ECD| IPA requires a positive occupied-virtual occupation difference")
               END IF
               rotatory_strength(itrans) = (occupation_numbers(iocc) - occupation_numbers(ivirt))*AIMAG( &
                                           SUM(CONJG(e1_mo(:, ivirt, iocc))*m1_mo(:, ivirt, iocc)))
               IF (write_transition_matrices) THEN
                  DO idir = 1, ncart
                     WRITE (e1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
                        itrans, transition_energy(itrans)*evolt, ispin, iocc, ivirt, idir, &
                        REAL(e1_mo(idir, ivirt, iocc), KIND=dp), AIMAG(e1_mo(idir, ivirt, iocc)), &
                        occupation_numbers(iocc) - occupation_numbers(ivirt)
                     WRITE (m1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
                        itrans, transition_energy(itrans)*evolt, ispin, iocc, ivirt, idir, &
                        REAL(m1_mo(idir, ivirt, iocc), KIND=dp), AIMAG(m1_mo(idir, ivirt, iocc)), &
                        occupation_numbers(iocc) - occupation_numbers(ivirt)
                  END DO
               END IF
            END DO
         END DO
         IF (use_velocity_gauge) THEN
            DEALLOCATE (velocity_mo_matrix)
         END IF
         DEALLOCATE (e1_mo, m1_mo)
      END DO
      IF (use_velocity_gauge) THEN
         CALL dbcsr_deallocate_matrix_set(velocity_ao)
      ELSE
         CALL dbcsr_deallocate_matrix_set(magnetic_ao)
         CALL dbcsr_deallocate_matrix_set(position_ao)
      END IF
      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T2,A)') "ECD| PERIODIC NONE E1 transition matrix written to "//TRIM(e1_filename)
         WRITE (iounit, '(T2,A)') "ECD| PERIODIC NONE M1 transition matrix written to "//TRIM(m1_filename)
      END IF
      CALL sort(transition_energy, ntrans_total, transition_order)
      IF (write_transition_matrices) THEN
         CALL close_file(e1_unit)
         CALL close_file(m1_unit)
      END IF
      CALL ecd_print_sticks(transition_energy, rotatory_strength, transition_order)
      DEALLOCATE (transition_order, transition_energy, rotatory_strength)
      CALL timestop(handle)

   END SUBROUTINE ecd_spectra_none_ipa

! **************************************************************************************************
!> \brief Assemble PERIODIC NONE TDDFPT ECD sticks from response roots.
!> \param qs_env Quickstep environment
!> \param evects TDDFPT response eigenvectors
!> \param evals TDDFPT excitation energies
!> \param gs_mos Active ground-state orbitals used by TDDFPT
!> \param mult TDDFPT multiplicity
! **************************************************************************************************
   SUBROUTINE ecd_spectra_none_tddfpt(qs_env, evects, evals, gs_mos, mult)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: evects
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), POINTER, INTENT(IN) :: gs_mos
      INTEGER, INTENT(IN)                                :: mult

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_spectra_none_tddfpt'

      INTEGER                                            :: ecd_gauge, handle, idir, ispin, istate, nactive, nspins, &
                                                            nstates, reference
      LOGICAL                                            :: do_ecd, pseudopotential_commutator, ref_point_explicit
      REAL(KIND=dp), DIMENSION(ncart)                    :: reference_point, transition_magnetic, transition_position
      COMPLEX(KIND=dp), DIMENSION(ncart)                 :: state_e1, state_m1
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE     :: magnetic_spin, position_spin
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: rotatory_strength
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: magnetic_active_action, position_active_action
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_ao, position_ao
      TYPE(section_vals_type), POINTER                   :: ecd_section

      CALL timeset(routineN, handle)

      NULLIFY (ecd_section, fm_struct, magnetic_ao, position_ao, ref_point)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")
      CALL section_vals_val_get(ecd_section, "DO_ECD", l_val=do_ecd)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL section_vals_val_get(ecd_section, "GAUGE", i_val=ecd_gauge)
      IF (ecd_gauge /= tddfpt_dipole_length .AND. ecd_gauge /= tddfpt_dipole_velocity) THEN
         CPABORT("ECD| PERIODIC NONE ECD supports GAUGE LENGTH and VELOCITY")
      END IF
      IF (ecd_gauge == tddfpt_dipole_velocity) CPABORT("ECD| TDDFPT ECD currently supports GAUGE LENGTH")

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
         IF (ref_point_explicit) CPWARN("ECD| REF_POINT is used only with REFERENCE USER_DEFINED")
         CALL get_reference_point(reference_point, qs_env=qs_env, reference=reference)
      END IF

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
         ! TDDFPT amplitudes already include the state transition weight.
         state_e1 = CMPLX(-transition_position, 0.0_dp, KIND=dp)
         state_m1 = CMPLX(0.0_dp, 0.5_dp*transition_magnetic, KIND=dp)
         rotatory_strength(istate) = AIMAG(SUM(CONJG(state_e1)*state_m1))
      END DO
      CALL ecd_print_sticks(evals, rotatory_strength)

      DO ispin = 1, nspins
         DO idir = 1, ncart
            CALL cp_fm_release(position_active_action(idir, ispin))
            CALL cp_fm_release(magnetic_active_action(idir, ispin))
         END DO
      END DO
      DEALLOCATE (magnetic_active_action, position_active_action, magnetic_spin, position_spin, rotatory_strength)

      CALL timestop(handle)

   END SUBROUTINE ecd_spectra_none_tddfpt

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
!> \brief Build direct length and orbital magnetic AO operators for PERIODIC NONE.
!> \param qs_env Quickstep environment
!> \param pseudopotential_commutator Include the complete GTH non-local magnetic correction
!> \param reference_point Common reference point for all molecular operators
!> \param position_ao Three direct AO-space position matrices
!> \param magnetic_ao Three direct AO-space angular-momentum matrices
! **************************************************************************************************
   SUBROUTINE ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, reference_point, &
                                                    position_ao, magnetic_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: reference_point
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER, &
         INTENT(INOUT)                                   :: magnetic_ao, position_ao

      INTEGER                                            :: idir
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_nl, matrix_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (cell, dft_control, magnetic_nl, matrix_s, particle_set, qs_kind_set, sab_orb, sap_ppnl)
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

      END DO

      CALL build_local_moment_matrix(qs_env, position_ao, 1, ref_point=reference_point)
      CALL build_local_magmom_matrix(qs_env, magnetic_ao, 1, ref_point=reference_point)

      IF (pseudopotential_commutator .AND. ASSOCIATED(sap_ppnl)) THEN
         CALL dbcsr_allocate_matrix_set(magnetic_nl, ncart)
         DO idir = 1, ncart
            CALL dbcsr_init_p(magnetic_nl(idir)%matrix)
            CALL dbcsr_create(magnetic_nl(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(magnetic_nl(idir)%matrix, sab_orb)
            CALL dbcsr_set(magnetic_nl(idir)%matrix, 0.0_dp)

         END DO
         ! r_R0 x [r_R0, V_nl] is the non-local contribution required by
         ! the orbital magnetic dipole.
         CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                               particle_set, cell, matrix_rxrv=magnetic_nl, ref_point=reference_point)
         DO idir = 1, ncart
            CALL dbcsr_add(magnetic_ao(idir)%matrix, magnetic_nl(idir)%matrix, 1.0_dp, 1.0_dp)
         END DO
         CALL dbcsr_deallocate_matrix_set(magnetic_nl)
      END IF

   END SUBROUTINE ecd_build_molecular_operators_ao_none

! **************************************************************************************************
!> \brief Build the raw molecular velocity-generator AO operators for PERIODIC NONE.
!> \param qs_env Quickstep environment
!> \param pseudopotential_commutator Include the GTH non-local velocity commutator
!> \param velocity_ao Three real antisymmetric AO derivative matrices
! **************************************************************************************************
   SUBROUTINE ecd_build_molecular_velocity_ao_none(qs_env, pseudopotential_commutator, velocity_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER, INTENT(INOUT) :: velocity_ao

      INTEGER                                            :: idir
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, velocity_nl
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (cell, dft_control, matrix_s, particle_set, qs_kind_set, sab_orb, sap_ppnl, velocity_nl)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_s=matrix_s, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set, sab_orb=sab_orb, &
                      sap_ppnl=sap_ppnl)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(sab_orb))

      CALL dbcsr_allocate_matrix_set(velocity_ao, ncart)
      DO idir = 1, ncart
         CALL dbcsr_init_p(velocity_ao(idir)%matrix)
         CALL dbcsr_create(velocity_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(velocity_ao(idir)%matrix, sab_orb)
         CALL dbcsr_set(velocity_ao(idir)%matrix, 0.0_dp)
      END DO
      CALL build_lin_mom_matrix(qs_env, velocity_ao)

      IF (pseudopotential_commutator .AND. ASSOCIATED(sap_ppnl)) THEN
         CALL dbcsr_allocate_matrix_set(velocity_nl, ncart)
         DO idir = 1, ncart
            CALL dbcsr_init_p(velocity_nl(idir)%matrix)
            CALL dbcsr_create(velocity_nl(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(velocity_nl(idir)%matrix, sab_orb)
            CALL dbcsr_set(velocity_nl(idir)%matrix, 0.0_dp)
         END DO
         CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                               particle_set, cell, matrix_rv=velocity_nl)
         DO idir = 1, ncart
            CALL dbcsr_add(velocity_ao(idir)%matrix, velocity_nl(idir)%matrix, 1.0_dp, 1.0_dp)
         END DO
         CALL dbcsr_deallocate_matrix_set(velocity_nl)
      END IF

   END SUBROUTINE ecd_build_molecular_velocity_ao_none

END MODULE qs_ecd_methods
