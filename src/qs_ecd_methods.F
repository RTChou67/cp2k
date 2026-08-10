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
   USE cell_types,                      ONLY: cell_type,&
                                              use_perd_none
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_generate_filename,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE external_potential_types,        ONLY: all_potential_type
   USE input_constants,                 ONLY: do_method_gapw,&
                                              do_method_gpw,&
                                              no_sf_tddfpt,&
                                              tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity,&
                                              tddfpt_kernel_full,&
                                              tddfpt_kernel_none,&
                                              use_mom_ref_user
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_ecd_operator_utils,           ONLY: ecd_build_e1_m1_mo,&
                                              ecd_build_length_ao,&
                                              ecd_mo_operators_type
   USE qs_ecd_spectra_xyz,              ONLY: ecd_spectra_xyz_ipa
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_methods'
   INTEGER, PARAMETER, PRIVATE          :: ncart = 3

   PUBLIC :: ecd_spectra, ecd_spectra_none_tddfpt

CONTAINS

! **************************************************************************************************
!> \brief Read and validate the ECD configuration for the active calculation path.
! **************************************************************************************************
   SUBROUTINE ecd_read_input(qs_env, do_ecd, ecd_gauge, ecd_gauge_explicit, &
                             pseudopotential_commutator, reference_point)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: do_ecd
      INTEGER, INTENT(OUT)                               :: ecd_gauge
      LOGICAL, INTENT(OUT)                               :: ecd_gauge_explicit, &
                                                            pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(OUT)       :: reference_point

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_read_input'

      INTEGER                                            :: handle, ikind, poisson_periodic, &
                                                            reference, spinflip, tddfpt_kernel
      LOGICAL                                            :: do_hfx, do_kpoints, do_soc, &
                                                            ref_point_explicit, tddfpt_enabled
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: ecd_section, hfx_section, &
                                                            poisson_section, soc_section, &
                                                            tddfpt_section

      CALL timeset(routineN, handle)

      NULLIFY (all_potential, cell, dft_control, ecd_section, hfx_section, particle_set, qs_kind_set, ref_point, &
               soc_section, tddfpt_section)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")
      CALL section_vals_val_get(ecd_section, "DO_ECD", l_val=do_ecd)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL section_vals_val_get(ecd_section, "GAUGE", i_val=ecd_gauge, explicit=ecd_gauge_explicit)
      IF (ecd_gauge /= tddfpt_dipole_length .AND. ecd_gauge /= tddfpt_dipole_velocity) THEN
         CPABORT("ECD| ECD supports GAUGE LENGTH and VELOCITY")
      END IF
      CALL section_vals_val_get(ecd_section, "PSEUDOPOTENTIAL_COMMUTATOR", &
                                l_val=pseudopotential_commutator)
      CALL section_vals_val_get(ecd_section, "REFERENCE", i_val=reference)
      CALL section_vals_val_get(ecd_section, "REF_POINT", explicit=ref_point_explicit)
      IF (ref_point_explicit) CALL section_vals_val_get(ecd_section, "REF_POINT", r_vals=ref_point)
      IF (reference == use_mom_ref_user) THEN
         IF (.NOT. ref_point_explicit) CPABORT("ECD| REFERENCE USER_DEFINED requires REF_POINT")
         IF (.NOT. ASSOCIATED(ref_point)) CPABORT("ECD| REF_POINT expects three coordinates")
         IF (SIZE(ref_point) /= ncart) CPABORT("ECD| REF_POINT expects three coordinates")
      ELSEIF (ref_point_explicit) THEN
         CPWARN("ECD| REF_POINT is used only with REFERENCE USER_DEFINED")
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
         IF (ecd_gauge_explicit .AND. ecd_gauge == tddfpt_dipole_length) THEN
            NULLIFY (dft_control, qs_kind_set)
            CALL get_qs_env(qs_env, dft_control=dft_control, qs_kind_set=qs_kind_set)
            CPASSERT(ASSOCIATED(dft_control))
            CPASSERT(ASSOCIATED(qs_kind_set))
            SELECT CASE (dft_control%qs_control%method_id)
            CASE (do_method_gapw)
               DO ikind = 1, SIZE(qs_kind_set)
                  NULLIFY (all_potential)
                  CALL get_qs_kind(qs_kind_set(ikind), all_potential=all_potential)
                  IF (.NOT. ASSOCIATED(all_potential)) THEN
                     CPABORT("ECD| PERIODIC XYZ GAPW GAUGE LENGTH requires POTENTIAL ALL for every KIND")
                  END IF
               END DO
            CASE (do_method_gpw)
               CALL section_vals_val_get(ecd_section, "PSEUDOPOTENTIAL_COMMUTATOR", &
                                         l_val=pseudopotential_commutator)
               IF (.NOT. pseudopotential_commutator) THEN
                  CPABORT("ECD| PERIODIC XYZ GPW GAUGE LENGTH requires PSEUDOPOTENTIAL_COMMUTATOR T")
               END IF
            CASE DEFAULT
               CPABORT("ECD| PERIODIC XYZ GAUGE LENGTH supports METHOD GAPW and METHOD GPW")
            END SELECT
         END IF
      ELSE
         CPABORT("ECD| Current implementation requires CELL%PERIODIC NONE or XYZ")
      END IF

      IF (ALL(cell%perd == 1)) THEN
         CALL get_qs_env(qs_env, particle_set=particle_set)
         CPASSERT(ASSOCIATED(particle_set))
         IF (reference == use_mom_ref_user) THEN
            CALL get_reference_point(reference_point, qs_env=qs_env, reference=reference, ref_point=ref_point, &
                                     ifirst=1, ilast=SIZE(particle_set))
         ELSE
            CALL get_reference_point(reference_point, qs_env=qs_env, reference=reference, &
                                     ifirst=1, ilast=SIZE(particle_set))
         END IF
      ELSEIF (reference == use_mom_ref_user) THEN
         CALL get_reference_point(reference_point, qs_env=qs_env, reference=reference, ref_point=ref_point)
      ELSE
         CALL get_reference_point(reference_point, qs_env=qs_env, reference=reference)
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
            IF (ecd_gauge == tddfpt_dipole_velocity) THEN
               CPABORT("ECD| TDDFPT ECD currently supports GAUGE LENGTH")
            END IF
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_read_input

! **************************************************************************************************
!> \brief Run independent-particle ECD after optional TDDFPT has completed.
! **************************************************************************************************
   SUBROUTINE ecd_spectra(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_spectra'

      INTEGER                                            :: ecd_gauge, handle
      LOGICAL                                            :: do_ecd, ecd_gauge_explicit, &
                                                            pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart)                    :: reference_point
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(dft_control))
      CALL ecd_read_input(qs_env, do_ecd, ecd_gauge, ecd_gauge_explicit, &
                          pseudopotential_commutator, reference_point)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF
      IF (dft_control%tddfpt2_control%enabled .AND. &
          dft_control%tddfpt2_control%kernel == tddfpt_kernel_full) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (ALL(cell%perd == 0) .OR. &
          (ecd_gauge_explicit .AND. ecd_gauge == tddfpt_dipole_length)) THEN
         CALL ecd_spectra_molecular_ipa(qs_env, ecd_gauge, pseudopotential_commutator, reference_point)
      ELSEIF (ALL(cell%perd == 1)) THEN
         CALL ecd_spectra_xyz_ipa(qs_env, pseudopotential_commutator, reference_point)
      ELSE
         CPABORT("ECD| Current implementation requires CELL%PERIODIC NONE or XYZ")
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_spectra

! **************************************************************************************************
!> \brief Assemble molecular IPA ECD sticks from open-boundary or periodic-supercell orbitals.
!> \param qs_env Quickstep environment
!> \param ecd_gauge E1/M1 gauge selected by the ECD input
!> \param pseudopotential_commutator Include the non-local pseudopotential contribution
!> \param reference_point Common E1/M1 reference point
! **************************************************************************************************
   SUBROUTINE ecd_spectra_molecular_ipa(qs_env, ecd_gauge, pseudopotential_commutator, &
                                        reference_point)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ecd_gauge
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: reference_point

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ecd_spectra_molecular_ipa'

      CHARACTER(LEN=default_path_length)                 :: e1_filename, m1_filename
      INTEGER                                            :: e1_unit, handle, homo, idir, iocc, &
                                                            iounit, ispin, itrans, ivirt, m1_unit, &
                                                            nmo, nspins, ntrans_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: transition_order
      LOGICAL                                            :: periodic_supercell, use_velocity_gauge, &
                                                            write_transition_matrices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rotatory_strength, transition_energy
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(ecd_mo_operators_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: operators
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (cell, eigenvalues, logger, mos, occupation_numbers, scf_control)
      CALL get_qs_env(qs_env, cell=cell)
      CPASSERT(ASSOCIATED(cell))
      periodic_supercell = ALL(cell%perd == 1)
      use_velocity_gauge = ecd_gauge == tddfpt_dipole_velocity
      IF (periodic_supercell .AND. use_velocity_gauge) THEN
         CPABORT("ECD| PERIODIC XYZ molecular-supercell ECD uses GAUGE LENGTH")
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
         IF (periodic_supercell) THEN
            WRITE (e1_unit, '(A)') "# PERIODIC XYZ molecular-supercell length-gauge IPA E1 transition matrix"
         ELSE
            WRITE (e1_unit, '(A)') "# PERIODIC NONE IPA E1 transition matrix"
         END IF
         IF (use_velocity_gauge) THEN
            WRITE (e1_unit, '(A)') "# mu_ai = i*v_ai/(epsilon_a-epsilon_i)"
         ELSE
            WRITE (e1_unit, '(A)') "# mu_ai = -r_ai"
         END IF
         WRITE (e1_unit, '(A)') "# records follow the occupied-empty MO traversal order"
         WRITE (e1_unit, '(A)') "# state energy_eV spin occupied virtual component e1_real e1_imag occupation_difference"
         CALL open_file(m1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=m1_unit)
         IF (periodic_supercell) THEN
            WRITE (m1_unit, '(A)') "# PERIODIC XYZ molecular-supercell length-gauge IPA M1 transition matrix"
         ELSE
            WRITE (m1_unit, '(A)') "# PERIODIC NONE IPA M1 transition matrix"
         END IF
         IF (use_velocity_gauge) THEN
            WRITE (m1_unit, '(A)') "# m_ai = -[r_velocity cross v]_ai/2"
         ELSE
            WRITE (m1_unit, '(A)') "# m_ai = i*Lambda_ai/2"
         END IF
         WRITE (m1_unit, '(A)') "# records follow the occupied-empty MO traversal order"
         WRITE (m1_unit, '(A)') "# state energy_eV spin occupied virtual component m1_real m1_imag occupation_difference"
      END IF
      CALL ecd_build_e1_m1_mo(qs_env, ecd_gauge, pseudopotential_commutator, reference_point, &
                              operators, reject_degenerate_velocity=.TRUE.)
      IF (SIZE(operators) /= nspins) CPABORT("ECD| E1/M1 spin dimensions are inconsistent")

      itrans = 0
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, homo=homo, nmo=nmo, &
                         occupation_numbers=occupation_numbers)
         CPASSERT(ASSOCIATED(eigenvalues))
         CPASSERT(ASSOCIATED(occupation_numbers))
         IF (homo <= 0 .OR. nmo <= homo) CYCLE
         DO iocc = 1, homo
            DO ivirt = homo + 1, nmo
               itrans = itrans + 1
               transition_energy(itrans) = eigenvalues(ivirt) - eigenvalues(iocc)
               IF (transition_energy(itrans) <= 0.0_dp) CPABORT("ECD| IPA requires positive occupied-virtual gaps")
               IF (occupation_numbers(iocc) - occupation_numbers(ivirt) <= 0.0_dp) THEN
                  CPABORT("ECD| IPA requires a positive occupied-virtual occupation difference")
               END IF
               rotatory_strength(itrans) = (occupation_numbers(iocc) - occupation_numbers(ivirt))*AIMAG( &
                                           SUM(CONJG(operators(ispin)%e1(:, ivirt, iocc))*operators(ispin)%m1(:, ivirt, iocc)))
               IF (write_transition_matrices) THEN
                  DO idir = 1, ncart
                     WRITE (e1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
                        itrans, transition_energy(itrans)*evolt, ispin, iocc, ivirt, idir, &
                        REAL(operators(ispin)%e1(idir, ivirt, iocc), KIND=dp), &
                        AIMAG(operators(ispin)%e1(idir, ivirt, iocc)), &
                        occupation_numbers(iocc) - occupation_numbers(ivirt)
                     WRITE (m1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
                        itrans, transition_energy(itrans)*evolt, ispin, iocc, ivirt, idir, &
                        REAL(operators(ispin)%m1(idir, ivirt, iocc), KIND=dp), &
                        AIMAG(operators(ispin)%m1(idir, ivirt, iocc)), &
                        occupation_numbers(iocc) - occupation_numbers(ivirt)
                  END DO
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE (operators)
      IF (iounit > 0) THEN
         IF (periodic_supercell) THEN
            WRITE (iounit, '(/,T2,A)') "ECD| PERIODIC XYZ length-gauge E1 transition matrix written to "//TRIM(e1_filename)
            WRITE (iounit, '(T2,A)') "ECD| PERIODIC XYZ length-gauge M1 transition matrix written to "//TRIM(m1_filename)
         ELSE
            WRITE (iounit, '(/,T2,A)') "ECD| PERIODIC NONE E1 transition matrix written to "//TRIM(e1_filename)
            WRITE (iounit, '(T2,A)') "ECD| PERIODIC NONE M1 transition matrix written to "//TRIM(m1_filename)
         END IF
      END IF
      CALL sort(transition_energy, ntrans_total, transition_order)
      IF (write_transition_matrices) THEN
         CALL close_file(e1_unit)
         CALL close_file(m1_unit)
      END IF
      CALL ecd_print_sticks(transition_energy, rotatory_strength, transition_order)
      DEALLOCATE (transition_order, transition_energy, rotatory_strength)
      CALL timestop(handle)

   END SUBROUTINE ecd_spectra_molecular_ipa

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
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(IN), POINTER                             :: gs_mos
      INTEGER, INTENT(IN)                                :: mult

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ecd_spectra_none_tddfpt'

      COMPLEX(KIND=dp), DIMENSION(ncart)                 :: state_e1, state_m1
      INTEGER                                            :: ecd_gauge, handle, idir, ispin, istate, &
                                                            nactive, nspins, nstates
      LOGICAL                                            :: do_ecd, ecd_gauge_explicit, &
                                                            pseudopotential_commutator
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rotatory_strength
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: magnetic_spin, position_spin
      REAL(KIND=dp), DIMENSION(ncart)                    :: reference_point, transition_magnetic, &
                                                            transition_position
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: magnetic_active_action, &
                                                            position_active_action
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_ao, position_ao

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct, magnetic_ao, position_ao)
      CALL ecd_read_input(qs_env, do_ecd, ecd_gauge, ecd_gauge_explicit, &
                          pseudopotential_commutator, reference_point)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF
      IF (ecd_gauge == tddfpt_dipole_velocity) CPABORT("ECD| TDDFPT ECD currently supports GAUGE LENGTH")

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
      CALL ecd_build_length_ao(qs_env, pseudopotential_commutator, reference_point, &
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
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: order

      INTEGER                                            :: iounit, istate, source_state
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

END MODULE qs_ecd_methods
