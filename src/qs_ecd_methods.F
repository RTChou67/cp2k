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
   USE cell_types,                      ONLY: cell_type,&
                                              use_perd_none
   USE commutator_rpnl,                 ONLY: build_com_mom_nl
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE input_constants,                 ONLY: no_sf_tddfpt,&
                                              tddfpt_kernel_full,&
                                              tddfpt_kernel_none,&
                                              use_mom_ref_user
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE moments_utils,                   ONLY: get_reference_point
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_local_magmom_matrix,&
                                              build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE particle_types,                  ONLY: particle_type
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_methods'
   INTEGER, PARAMETER, PRIVATE          :: ecd_broadening_none = 1, &
                                           ecd_broadening_gaussian = 2, &
                                           ecd_broadening_lorentzian = 3, &
                                           ecd_spectrum_x_ev = 1, &
                                           ecd_spectrum_x_nm = 2, &
                                           ecd_spectrum_x_cm1 = 3, &
                                           ecd_spectrum_y_au = 1, &
                                           ecd_spectrum_y_cgs_1e40 = 2
   INTEGER, PARAMETER, PRIVATE          :: ecd_spectrum_points = 2001, ncart = 3
   REAL(KIND=dp), PARAMETER, PRIVATE    :: ecd_ev_to_cm1 = 8065.54394_dp, &
                                           ecd_hc_ev_nm = 1239.84198_dp, &
                                           ecd_rot_au_to_cgs_1e40 = 471.4436_dp

   PUBLIC :: ecd_lr_tddfpt_spectra_none, ecd_spectra

CONTAINS

! **************************************************************************************************
!> \brief Parse the ECD input section and dispatch the selected ECD data source.
!> \param qs_env Quickstep environment
!> \param ecd_section ECD input section
! **************************************************************************************************
   SUBROUTINE ecd_spectra(qs_env, ecd_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: ecd_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_spectra'

      CHARACTER(LEN=default_string_length)               :: broadening_label
      INTEGER                                            :: broadening, handle, iounit, ntransitions, reference
      INTEGER, DIMENSION(2)                              :: identity_band_range, transition_band_range
      INTEGER, DIMENSION(:), POINTER                     :: identity_band_range_ptr, &
                                                            transition_band_range_ptr
      LOGICAL                                            :: do_ecd, do_hfx, do_ipa, do_lr_tddft, pseudopotential_commutator
      LOGICAL                                            :: ref_point_explicit
      REAL(KIND=dp)                                      :: fwhm
      REAL(KIND=dp), DIMENSION(3)                        :: rcc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: hfx_section, ipa_section, lr_tddft_section, &
                                                            tddfpt_section

      CALL timeset(routineN, handle)

      NULLIFY (hfx_section, identity_band_range_ptr, logger, ipa_section, lr_tddft_section, ref_point, &
               tddfpt_section, transition_band_range_ptr)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL section_vals_val_get(ecd_section, "_SECTION_PARAMETERS_", l_val=do_ecd)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ipa_section => section_vals_get_subs_vals(ecd_section, "IPA")
      lr_tddft_section => section_vals_get_subs_vals(ecd_section, "LR_TDDFT")

      CALL section_vals_get(ipa_section, explicit=do_ipa)
      CALL section_vals_get(lr_tddft_section, explicit=do_lr_tddft)

      IF (do_ipa .AND. do_lr_tddft) THEN
         CPABORT("ECD| IPA and LR_TDDFT sources are mutually exclusive")
      END IF
      IF (.NOT. do_ipa .AND. .NOT. do_lr_tddft) do_ipa = .TRUE.

      CALL section_vals_val_get(ecd_section, "PSEUDOPOTENTIAL_COMMUTATOR", &
                                l_val=pseudopotential_commutator)
      CALL section_vals_val_get(ecd_section, "IDENTITY_BAND_RANGE", i_vals=identity_band_range_ptr)
      CALL section_vals_val_get(ecd_section, "REFERENCE", i_val=reference)
      CALL section_vals_val_get(ecd_section, "REF_POINT", explicit=ref_point_explicit)
      IF (ref_point_explicit) THEN
         CALL section_vals_val_get(ecd_section, "REF_POINT", r_vals=ref_point)
      END IF
      CALL section_vals_val_get(ecd_section, "BROADENING", i_val=broadening)
      CALL section_vals_val_get(ecd_section, "FWHM", r_val=fwhm)

      CALL section_vals_val_get(ipa_section, "NTRANSITIONS", i_val=ntransitions)
      CALL section_vals_val_get(ipa_section, "TRANSITION_BAND_RANGE", &
                                i_vals=transition_band_range_ptr)

      IF (SIZE(identity_band_range_ptr) /= 2) THEN
         CPABORT("ECD| IDENTITY_BAND_RANGE expects two integers")
      END IF
      IF (SIZE(transition_band_range_ptr) /= 2) THEN
         CPABORT("ECD| TRANSITION_BAND_RANGE expects two integers")
      END IF
      identity_band_range = identity_band_range_ptr(1:2)
      transition_band_range = transition_band_range_ptr(1:2)

      SELECT CASE (broadening)
      CASE (ecd_broadening_none)
         broadening_label = "NONE"
      CASE (ecd_broadening_gaussian)
         broadening_label = "GAUSSIAN"
      CASE (ecd_broadening_lorentzian)
         broadening_label = "LORENTZIAN"
      CASE DEFAULT
         CPABORT("ECD| Unknown BROADENING")
      END SELECT

      IF (reference == use_mom_ref_user) THEN
         IF (.NOT. ref_point_explicit) THEN
            CPABORT("ECD| REFERENCE USER_DEFINED requires REF_POINT")
         END IF
         IF (.NOT. ASSOCIATED(ref_point) .OR. SIZE(ref_point) /= ncart) THEN
            CPABORT("ECD| REF_POINT expects three coordinates")
         END IF
         CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
      ELSE
         IF (ref_point_explicit) THEN
            CPWARN("ECD| REF_POINT is used only with REFERENCE USER_DEFINED")
         END IF
         CALL get_reference_point(rcc, qs_env=qs_env, reference=reference)
      END IF

      IF (iounit > 0 .AND. .NOT. do_lr_tddft) THEN
         WRITE (iounit, '(/,T2,A)') "ECD| Electronic circular dichroism spectra"
         WRITE (iounit, '(T2,A,T55,L1)') "ECD| IPA source", do_ipa
         WRITE (iounit, '(T2,A,T55,L1)') "ECD| LR-TDDFT source", do_lr_tddft
         WRITE (iounit, '(T2,A,T55,L1)') "ECD| GTH non-local magnetic correction", &
            pseudopotential_commutator
         WRITE (iounit, '(T2,A,T55,2I8)') "ECD| Identity band range", identity_band_range
         WRITE (iounit, '(T2,A,T55,3F12.6)') "ECD| Molecular reference point [bohr]", rcc
         WRITE (iounit, '(T2,A,T55,A)') "ECD| Broadening", TRIM(broadening_label)
         WRITE (iounit, '(T2,A,T55,ES14.6)') "ECD| FWHM [a.u.]", fwhm
         WRITE (iounit, '(T2,A,T55,I8)') "ECD| IPA transitions", ntransitions
         WRITE (iounit, '(T2,A,T55,2I8)') "ECD| IPA transition band range", &
            transition_band_range
      END IF

      CALL ecd_require_nonperiodic(qs_env)

      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CPABORT("ECD| Molecular ECD with HFX requires a non-local magnetic correction")
      END IF

      IF (do_lr_tddft) THEN
         tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")
         CALL section_vals_get(tddfpt_section, explicit=do_lr_tddft)
         IF (.NOT. do_lr_tddft) THEN
            CPABORT("ECD| LR_TDDFT source requires the same-level TDDFPT section")
         END IF
         CALL timestop(handle)
         RETURN
      ELSE
         CALL ecd_ipa_spectra_none(qs_env, ntransitions, transition_band_range, &
                                   identity_band_range, pseudopotential_commutator, &
                                   rcc, broadening, fwhm, ecd_section)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_spectra

! **************************************************************************************************
!> \brief Enforce the currently supported non-periodic ECD calculation domain.
!> \param qs_env Quickstep environment
! **************************************************************************************************
   SUBROUTINE ecd_require_nonperiodic(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_require_nonperiodic'

      INTEGER                                            :: handle, poisson_periodic
      LOGICAL                                            :: do_kpoints
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: input, poisson_section

      CALL timeset(routineN, handle)

      NULLIFY (cell, input, poisson_section)
      CALL get_qs_env(qs_env, cell=cell, do_kpoints=do_kpoints, input=input)

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(input))

      IF (ANY(cell%perd /= 0)) THEN
         CPABORT("ECD| Current implementation supports only CELL%PERIODIC NONE")
      END IF

      poisson_section => section_vals_get_subs_vals(input, "DFT%POISSON")
      CALL section_vals_val_get(poisson_section, "PERIODIC", i_val=poisson_periodic)
      IF (poisson_periodic /= use_perd_none) THEN
         CPABORT("ECD| Current implementation supports only POISSON%PERIODIC NONE")
      END IF

      IF (do_kpoints) THEN
         CPABORT("ECD| Current implementation supports only Gamma-only PERIODIC NONE")
      END IF

      CALL timestop(handle)

   END SUBROUTINE ecd_require_nonperiodic

! **************************************************************************************************
!> \brief Independent-particle ECD entry point for non-periodic molecular calculations.
!> \param qs_env Quickstep environment
!> \param ntransitions Number of requested transitions
!> \param transition_band_range Occupied and virtual band range selector
!> \param identity_band_range Periodic-only identity window, ignored for molecular operators
!> \param pseudopotential_commutator Explicit GTH non-local magnetic correction switch
!> \param rcc Molecular operator reference point
! **************************************************************************************************
   SUBROUTINE ecd_ipa_spectra_none(qs_env, ntransitions, transition_band_range, &
                                   identity_band_range, pseudopotential_commutator, &
                                   rcc, broadening, fwhm, ecd_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ntransitions
      INTEGER, DIMENSION(2), INTENT(IN)                  :: transition_band_range, &
                                                            identity_band_range
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      INTEGER, INTENT(IN)                                :: broadening
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: rcc
      REAL(KIND=dp), INTENT(IN)                          :: fwhm
      TYPE(section_vals_type), POINTER                   :: ecd_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_ipa_spectra_none'

      INTEGER                                            :: handle, homo, iocc, &
                                                            ispin, istate, itrans, ivirt, nao, nmo, nspins, &
                                                            nstates, ntrans_total, occ_first, occ_last, &
                                                            stick_unit, virt_first, virt_last
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds, trans_iocc, trans_ispin, trans_ivirt
      LOGICAL                                            :: added_mos_explicit
      REAL(KIND=dp)                                      :: energy_checksum, fsum, oscillator_strength, &
                                                            rotatory_strength, rsum, transition_weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: oscillator_stick, &
                                                            rotatory_stick, spectrum_energy, &
                                                            transition_energy
      REAL(KIND=dp), DIMENSION(ncart)                    :: magnetic, position
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: magnetic_mo, position_mo
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: input, scf_section

      CALL timeset(routineN, handle)

      NULLIFY (eigenvalues, input, logger, mo_coeff, mos, occupation_numbers, scf_control, scf_section)

      IF (ntransitions < 1) THEN
         CPABORT("ECD| IPA requires NTRANSITIONS larger than zero")
      END IF
      IF (transition_band_range(1) > transition_band_range(2)) THEN
         CPABORT("ECD| IPA TRANSITION_BAND_RANGE lower bound exceeds upper bound")
      END IF
      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, input=input, mos=mos, scf_control=scf_control)
      CPASSERT(ASSOCIATED(input))
      CPASSERT(ASSOCIATED(mos))
      CPASSERT(ASSOCIATED(scf_control))

      IF (scf_control%use_ot) THEN
         CPABORT("ECD| IPA currently requires diagonalization SCF with explicit SCF%ADDED_MOS")
      END IF

      scf_section => section_vals_get_subs_vals(input, "DFT%SCF")
      CALL section_vals_val_get(scf_section, "ADDED_MOS", explicit=added_mos_explicit)
      IF (.NOT. added_mos_explicit) THEN
         CPABORT("ECD| IPA requires explicit SCF%ADDED_MOS to provide virtual MOs")
      END IF

      nspins = SIZE(mos)
      ntrans_total = 0
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), homo=homo, nao=nao, nmo=nmo)
         IF (homo <= 0) THEN
            CPABORT("ECD| IPA requires occupied MOs in every spin channel")
         END IF
         IF (nmo <= homo) THEN
            CPABORT("ECD| IPA requires virtual MOs from SCF%ADDED_MOS")
         END IF
         IF (nao <= 0) THEN
            CPABORT("ECD| IPA requires AO-basis MO coefficients")
         END IF
         CALL ecd_decode_transition_range(transition_band_range, nmo, homo, occ_first, &
                                          occ_last, virt_first, virt_last)
         ntrans_total = ntrans_total + (occ_last - occ_first + 1)*(virt_last - virt_first + 1)
      END DO
      IF (ntrans_total <= 0) THEN
         CPABORT("ECD| No independent-particle transitions available")
      END IF

      ALLOCATE (transition_energy(ntrans_total), trans_ispin(ntrans_total), &
                trans_iocc(ntrans_total), trans_ivirt(ntrans_total), inds(ntrans_total))
      transition_energy = 0.0_dp

      itrans = 0
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), homo=homo, nmo=nmo, eigenvalues=eigenvalues)
         CPASSERT(ASSOCIATED(eigenvalues))
         CALL ecd_decode_transition_range(transition_band_range, nmo, homo, occ_first, &
                                          occ_last, virt_first, virt_last)
         DO iocc = occ_first, occ_last
            DO ivirt = virt_first, virt_last
               itrans = itrans + 1
               transition_energy(itrans) = eigenvalues(ivirt) - eigenvalues(iocc)
               IF (transition_energy(itrans) <= 0.0_dp) THEN
                  CPABORT("ECD| IPA requires positive occupied-virtual gaps")
               END IF
               trans_ispin(itrans) = ispin
               trans_iocc(itrans) = iocc
               trans_ivirt(itrans) = ivirt
            END DO
         END DO
      END DO

      CALL sort(transition_energy, ntrans_total, inds)
      nstates = MIN(ntransitions, ntrans_total)
      IF (ntransitions > ntrans_total) THEN
         CPWARN("ECD| Requested more IPA transitions than available")
      END IF
      ALLOCATE (spectrum_energy(nstates), oscillator_stick(nstates), rotatory_stick(nstates))
      spectrum_energy = transition_energy(1:nstates)
      oscillator_stick = 0.0_dp
      rotatory_stick = 0.0_dp

      IF (ANY(identity_band_range /= -1)) THEN
         CPWARN("ECD| IDENTITY_BAND_RANGE is ignored for PERIODIC NONE direct molecular operators")
      END IF

      ALLOCATE (magnetic_mo(ncart, nspins), position_mo(ncart, nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, nao=nao, nmo=nmo)
         CPASSERT(ASSOCIATED(mo_coeff))
         CALL ecd_build_molecular_operators_none(qs_env, mo_coeff, nao, nmo, &
                                                  pseudopotential_commutator, rcc, &
                                                  position_mo(:, ispin), magnetic_mo(:, ispin))
      END DO

      logger => cp_get_default_logger()
      energy_checksum = 0.0_dp
      fsum = 0.0_dp
      rsum = 0.0_dp
      stick_unit = cp_print_key_unit_nr(logger, ecd_section, "PRINT%ROTATORY_STRENGTHS", &
                                        extension=".Log")
      IF (stick_unit > 0) CALL ecd_print_ipa_header(stick_unit, pseudopotential_commutator)

      DO istate = 1, nstates
         itrans = inds(istate)
         ispin = trans_ispin(itrans)
         iocc = trans_iocc(itrans)
         ivirt = trans_ivirt(itrans)
         CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, nmo=nmo, &
                         occupation_numbers=occupation_numbers)
         CPASSERT(ASSOCIATED(eigenvalues))
         CPASSERT(ASSOCIATED(occupation_numbers))
         transition_weight = occupation_numbers(iocc) - occupation_numbers(ivirt)
         IF (transition_weight <= 0.0_dp) THEN
            CPABORT("ECD| IPA requires a positive occupied-virtual occupation difference")
         END IF
         CALL ecd_transition_quantities(position_mo(:, ispin), magnetic_mo(:, ispin), &
                                        eigenvalues, iocc, ivirt, transition_weight, position, magnetic, &
                                        rotatory_strength, oscillator_strength)
         rotatory_stick(istate) = rotatory_strength
         oscillator_stick(istate) = oscillator_strength
         energy_checksum = energy_checksum + transition_energy(istate)**2
         fsum = fsum + oscillator_strength**2
         rsum = rsum + rotatory_strength**2

         IF (stick_unit > 0) THEN
            CALL ecd_print_ipa_transition(stick_unit, istate, &
                                          transition_energy(istate), ispin, iocc, ivirt, &
                                          transition_weight, position, magnetic, &
                                          rotatory_strength, oscillator_strength)
         END IF
      END DO

      IF (stick_unit > 0) THEN
         CALL ecd_print_ipa_footer(stick_unit, SQRT(energy_checksum), &
                                   SQRT(fsum), SQRT(rsum))
         CALL cp_print_key_finished_output(stick_unit, logger, ecd_section, "PRINT%ROTATORY_STRENGTHS")
      END IF

      CALL ecd_print_ipa_spectrum(ecd_section, spectrum_energy, rotatory_stick, oscillator_stick, &
                                  broadening, fwhm)

      DO ispin = 1, nspins
         DO iocc = 1, ncart
            CALL cp_fm_release(magnetic_mo(iocc, ispin))
            CALL cp_fm_release(position_mo(iocc, ispin))
         END DO
      END DO
      DEALLOCATE (inds, oscillator_stick, rotatory_stick, spectrum_energy, trans_iocc, &
                  trans_ispin, trans_ivirt, transition_energy, magnetic_mo, position_mo)

      CALL timestop(handle)

   END SUBROUTINE ecd_ipa_spectra_none

! **************************************************************************************************
!> \brief Assemble PERIODIC NONE ECD sticks from converged TDDFPT response vectors.
!> \param qs_env Quickstep environment
!> \param evects TDDFPT response orbitals [spin, state]
!> \param evals TDDFPT excitation energies
!> \param gs_mos Ground-state occupied and active MO data
!> \param tddfpt_kernel TDDFPT kernel selector
!> \param mult TDDFPT excitation multiplicity
!> \param do_soc Whether the TDDFPT SOC route is active
!> \param spinflip TDDFPT spin-flip selector
! **************************************************************************************************
   SUBROUTINE ecd_lr_tddfpt_spectra_none(qs_env, evects, evals, gs_mos, tddfpt_kernel, mult, do_soc, spinflip)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: evects
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      INTEGER, INTENT(IN)                                :: mult, spinflip, tddfpt_kernel
      LOGICAL, INTENT(IN)                                :: do_soc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_lr_tddfpt_spectra_none'

      INTEGER                                            :: broadening, handle, idir, ispin, istate, iout, &
                                                            nactive, nocc, nspins, nstates, nstates_out, &
                                                            reference, state
      INTEGER, DIMENSION(:), POINTER                     :: state_list
      LOGICAL                                            :: do_ecd, do_hfx, do_ipa, do_lr_tddft, pseudopotential_commutator
      LOGICAL, DIMENSION(:), ALLOCATABLE                 :: selected
      LOGICAL                                            :: ref_point_explicit
      REAL(KIND=dp)                                      :: energy_checksum, fsum, fwhm, rsum
      REAL(KIND=dp), DIMENSION(ncart)                    :: rcc
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE     :: magnetic_spin, position_spin
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: magnetic_state, position_state
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: oscillator_stick, rotatory_stick, spectrum_energy
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: magnetic_active, position_active
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: magnetic_occ, position_occ
      TYPE(section_vals_type), POINTER                   :: ecd_section, hfx_section, ipa_section, &
                                                            lr_tddfpt_section

      CALL timeset(routineN, handle)

      NULLIFY (ecd_section, fm_struct, hfx_section, ipa_section, lr_tddfpt_section, ref_point, state_list)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")
      CALL section_vals_val_get(ecd_section, "_SECTION_PARAMETERS_", l_val=do_ecd)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ipa_section => section_vals_get_subs_vals(ecd_section, "IPA")
      lr_tddfpt_section => section_vals_get_subs_vals(ecd_section, "LR_TDDFT")
      CALL section_vals_get(ipa_section, explicit=do_ipa)
      CALL section_vals_get(lr_tddfpt_section, explicit=do_lr_tddft)
      IF (do_ipa .AND. do_lr_tddft) CPABORT("ECD| IPA and LR_TDDFT sources are mutually exclusive")
      IF (.NOT. do_lr_tddft) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (do_soc) CPABORT("ECD| LR_TDDFT ECD does not support TDDFPT SOC")
      IF (spinflip /= no_sf_tddfpt) CPABORT("ECD| LR_TDDFT ECD does not support spin-flip TDDFPT")
      IF (tddfpt_kernel /= tddfpt_kernel_full .AND. tddfpt_kernel /= tddfpt_kernel_none) THEN
         CPABORT("ECD| LR_TDDFT ECD supports TDDFPT KERNEL FULL and KERNEL NONE")
      END IF

      CALL ecd_require_nonperiodic(qs_env)
      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CPABORT("ECD| Molecular ECD with HFX requires a non-local magnetic correction")
      END IF

      CALL section_vals_val_get(ecd_section, "PSEUDOPOTENTIAL_COMMUTATOR", &
                                l_val=pseudopotential_commutator)
      CALL section_vals_val_get(ecd_section, "BROADENING", i_val=broadening)
      CALL section_vals_val_get(ecd_section, "FWHM", r_val=fwhm)
      CALL section_vals_val_get(ecd_section, "REFERENCE", i_val=reference)
      CALL section_vals_val_get(ecd_section, "REF_POINT", explicit=ref_point_explicit)
      IF (reference == use_mom_ref_user) THEN
         IF (.NOT. ref_point_explicit) CPABORT("ECD| REFERENCE USER_DEFINED requires REF_POINT")
         CALL section_vals_val_get(ecd_section, "REF_POINT", r_vals=ref_point)
         IF (SIZE(ref_point) /= ncart) CPABORT("ECD| REF_POINT expects three coordinates")
         CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
      ELSE
         IF (ref_point_explicit) CPWARN("ECD| REF_POINT is used only with REFERENCE USER_DEFINED")
         CALL get_reference_point(rcc, qs_env=qs_env, reference=reference)
      END IF

      nspins = SIZE(gs_mos)
      nstates = SIZE(evals)
      IF (SIZE(evects, 1) /= nspins .OR. SIZE(evects, 2) /= nstates) THEN
         CPABORT("ECD| TDDFPT eigenvector dimensions are inconsistent with the requested ECD states")
      END IF
      IF (nspins == 1 .AND. mult /= 1) THEN
         CPABORT("ECD| LR_TDDFT ECD supports restricted singlet TDDFPT states")
      END IF
      IF (nspins > 2) CPABORT("ECD| LR_TDDFT ECD supports at most two spin channels")

      CALL section_vals_val_get(lr_tddfpt_section, "STATE_LIST", i_vals=state_list)
      ALLOCATE (selected(nstates))
      selected = .FALSE.
      IF (SIZE(state_list) == 1 .AND. state_list(1) == -1) THEN
         selected = .TRUE.
      ELSE
         IF (ANY(state_list < 1)) THEN
            CPABORT("ECD| LR_TDDFT STATE_LIST uses -1 alone or positive state numbers")
         END IF
         DO iout = 1, SIZE(state_list)
            state = state_list(iout)
            IF (state > nstates) CPABORT("ECD| LR_TDDFT STATE_LIST exceeds TDDFPT NSTATES")
            selected(state) = .TRUE.
         END DO
      END IF
      nstates_out = COUNT(selected)
      IF (nstates_out == 0) CPABORT("ECD| LR_TDDFT STATE_LIST selected no TDDFPT states")

      ALLOCATE (position_occ(ncart, nspins), magnetic_occ(ncart, nspins))
      CALL ecd_build_lr_operator_actions(qs_env, gs_mos, pseudopotential_commutator, rcc, &
                                         position_occ, magnetic_occ)
      ALLOCATE (position_spin(nstates, ncart, nspins), magnetic_spin(nstates, ncart, nspins))
      position_spin = 0.0_dp
      magnetic_spin = 0.0_dp

      DO ispin = 1, nspins
         CALL cp_fm_get_info(position_occ(1, ispin), ncol_global=nocc)
         CALL cp_fm_get_info(evects(ispin, 1), ncol_global=nactive)
         IF (nocc == nactive) THEN
            DO idir = 1, ncart
               CALL cp_fm_trace(evects(ispin, :), position_occ(idir, ispin), &
                                position_spin(:, idir, ispin))
               CALL cp_fm_trace(evects(ispin, :), magnetic_occ(idir, ispin), &
                                magnetic_spin(:, idir, ispin))
            END DO
         ELSE
            CALL cp_fm_get_info(evects(ispin, 1), matrix_struct=fm_struct)
            ALLOCATE (position_active(ncart), magnetic_active(ncart))
            DO idir = 1, ncart
               CALL cp_fm_create(position_active(idir), fm_struct)
               CALL cp_fm_create(magnetic_active(idir), fm_struct)
               DO iout = 1, nactive
                  state = gs_mos(ispin)%index_active(iout)
                  CALL cp_fm_to_fm(position_occ(idir, ispin), position_active(idir), &
                                   ncol=1, source_start=state, target_start=iout)
                  CALL cp_fm_to_fm(magnetic_occ(idir, ispin), magnetic_active(idir), &
                                   ncol=1, source_start=state, target_start=iout)
               END DO
               CALL cp_fm_trace(evects(ispin, :), position_active(idir), &
                                position_spin(:, idir, ispin))
               CALL cp_fm_trace(evects(ispin, :), magnetic_active(idir), &
                                magnetic_spin(:, idir, ispin))
               CALL cp_fm_release(position_active(idir))
               CALL cp_fm_release(magnetic_active(idir))
            END DO
            DEALLOCATE (position_active, magnetic_active)
            NULLIFY (fm_struct)
         END IF
      END DO

      ALLOCATE (position_state(nstates, ncart), magnetic_state(nstates, ncart))
      IF (nspins == 1) THEN
         position_state = SQRT(2.0_dp)*position_spin(:, :, 1)
         magnetic_state = 0.5_dp*SQRT(2.0_dp)*magnetic_spin(:, :, 1)
      ELSE
         position_state = position_spin(:, :, 1) + position_spin(:, :, 2)
         magnetic_state = 0.5_dp*(magnetic_spin(:, :, 1) + magnetic_spin(:, :, 2))
      END IF

      ALLOCATE (spectrum_energy(nstates_out), oscillator_stick(nstates_out), rotatory_stick(nstates_out))
      energy_checksum = 0.0_dp
      fsum = 0.0_dp
      rsum = 0.0_dp
      iout = 0
      DO istate = 1, nstates
         IF (.NOT. selected(istate)) CYCLE
         iout = iout + 1
         spectrum_energy(iout) = evals(istate)
         rotatory_stick(iout) = -SUM(position_state(istate, :)*magnetic_state(istate, :))
         oscillator_stick(iout) = 2.0_dp/3.0_dp*evals(istate)*SUM(position_state(istate, :)**2)
         energy_checksum = energy_checksum + evals(istate)**2
         fsum = fsum + oscillator_stick(iout)**2
         rsum = rsum + rotatory_stick(iout)**2
      END DO

      CALL ecd_print_lr_tddft(ecd_section, pseudopotential_commutator, rcc, evals, selected, &
                              position_state, magnetic_state, rotatory_stick, oscillator_stick, &
                              SQRT(energy_checksum), SQRT(fsum), SQRT(rsum))
      CALL ecd_print_ipa_spectrum(ecd_section, spectrum_energy, rotatory_stick, oscillator_stick, broadening, fwhm)

      DO ispin = 1, nspins
         DO idir = 1, ncart
            CALL cp_fm_release(position_occ(idir, ispin))
            CALL cp_fm_release(magnetic_occ(idir, ispin))
         END DO
      END DO
      DEALLOCATE (magnetic_occ, position_occ, magnetic_spin, position_spin, magnetic_state, position_state, &
                  oscillator_stick, rotatory_stick, selected, spectrum_energy)

      CALL timestop(handle)

   END SUBROUTINE ecd_lr_tddfpt_spectra_none

! **************************************************************************************************
!> \brief Form direct AO operator actions on the occupied TDDFPT ground-state orbitals.
!> \param qs_env Quickstep environment
!> \param gs_mos Ground-state occupied MO data
!> \param pseudopotential_commutator Include GTH non-local magnetic correction
!> \param rcc Common molecular reference point
!> \param position_occ Position operator actions X*C_occ [direction, spin]
!> \param magnetic_occ Angular-momentum operator actions Lambda*C_occ [direction, spin]
! **************************************************************************************************
   SUBROUTINE ecd_build_lr_operator_actions(qs_env, gs_mos, pseudopotential_commutator, rcc, &
                                            position_occ, magnetic_occ)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: rcc
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(INOUT)   :: magnetic_occ, position_occ

      INTEGER                                            :: idir, ispin, nocc
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_ao, position_ao

      NULLIFY (fm_struct, magnetic_ao, position_ao)
      CALL ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, rcc, &
                                                  position_ao, magnetic_ao)
      DO ispin = 1, SIZE(gs_mos)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct, ncol_global=nocc)
         DO idir = 1, ncart
            CALL cp_fm_create(position_occ(idir, ispin), fm_struct)
            CALL cp_fm_create(magnetic_occ(idir, ispin), fm_struct)
            CALL cp_dbcsr_sm_fm_multiply(position_ao(idir)%matrix, gs_mos(ispin)%mos_occ, &
                                         position_occ(idir, ispin), ncol=nocc, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(magnetic_ao(idir)%matrix, gs_mos(ispin)%mos_occ, &
                                         magnetic_occ(idir, ispin), ncol=nocc, alpha=1.0_dp, beta=0.0_dp)
         END DO
         NULLIFY (fm_struct)
      END DO
      CALL dbcsr_deallocate_matrix_set(magnetic_ao)
      CALL dbcsr_deallocate_matrix_set(position_ao)

   END SUBROUTINE ecd_build_lr_operator_actions

! **************************************************************************************************
!> \brief Print PERIODIC NONE LR-TDDFT ECD sticks and checksums.
! **************************************************************************************************
   SUBROUTINE ecd_print_lr_tddft(ecd_section, pseudopotential_commutator, rcc, evals, selected, &
                                  position_state, magnetic_state, rotatory_stick, oscillator_stick, &
                                  energy_checksum, oscillator_checksum, rotatory_checksum)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: rcc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals, oscillator_stick, rotatory_stick
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: selected
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: magnetic_state, position_state
      REAL(KIND=dp), INTENT(IN)                          :: energy_checksum, oscillator_checksum, &
                                                            rotatory_checksum

      INTEGER                                            :: iout, istate, iounit
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      iounit = cp_print_key_unit_nr(logger, ecd_section, "PRINT%ROTATORY_STRENGTHS", extension=".Log")
      IF (iounit <= 0) RETURN

      WRITE (iounit, "(1X,A)") "", &
         "-------------------------------------------------------------------------------", &
         "-                  ECD LR-TDDFT transition states                              -", &
         "-------------------------------------------------------------------------------"
      WRITE (iounit, "(1X,A)") "ECD| PERIODIC NONE LR-TDDFT sticks use direct AO length and magnetic operators."
      WRITE (iounit, "(1X,A,L1)") "ECD| GTH non-local magnetic correction included: ", &
         pseudopotential_commutator
      WRITE (iounit, "(1X,A,3F12.6)") "ECD| Molecular reference point [bohr]", rcc
      WRITE (iounit, '(/,T8,A,T17,A,T31,A,T52,A,T69,A)') &
         "State", "Energy", "Rotatory", "Oscillator", "Source"
      WRITE (iounit, '(T8,A,T17,A,T31,A,T52,A,T69,A)') &
         "number", "eV", "strength", "strength", "TDDFPT"
      WRITE (iounit, '(T8,80("-"))')

      iout = 0
      DO istate = 1, SIZE(evals)
         IF (.NOT. selected(istate)) CYCLE
         iout = iout + 1
         WRITE (iounit, '(1X,A,T11,I7,T20,F11.5,T32,ES12.5E2,T53,ES12.5E2,T69,A)') &
            "ECD_LR|", istate, evals(istate)*evolt, rotatory_stick(iout), oscillator_stick(iout), "TDDFPT"
         WRITE (iounit, '(1X,A,T17,3(1X,ES13.6E2))') "ECD_LR_R|", position_state(istate, :)
         WRITE (iounit, '(1X,A,T17,3(1X,ES13.6E2))') "ECD_LR_M|", magnetic_state(istate, :)
      END DO
      WRITE (iounit, '(/,T2,A,E16.8)') 'ECD LR-TDDFT : CheckSum E = ', energy_checksum
      WRITE (iounit, '(T2,A,E16.8)') 'ECD LR-TDDFT : CheckSum F = ', oscillator_checksum
      WRITE (iounit, '(T2,A,E16.8)') 'ECD LR-TDDFT : CheckSum R = ', rotatory_checksum
      WRITE (iounit, "(1X,A)") "-------------------------------------------------------------------------------"
      CALL cp_print_key_finished_output(iounit, logger, ecd_section, "PRINT%ROTATORY_STRENGTHS")

   END SUBROUTINE ecd_print_lr_tddft

! **************************************************************************************************
!> \brief Decode the occupied and virtual windows used for IPA transition enumeration.
!> \param band_range Input pair
!> \param nmo Number of available MOs
!> \param homo Highest occupied MO index
!> \param occ_first First occupied MO
!> \param occ_last Last occupied MO
!> \param virt_first First virtual MO
!> \param virt_last Last virtual MO
! **************************************************************************************************
   SUBROUTINE ecd_decode_transition_range(band_range, nmo, homo, occ_first, occ_last, &
                                          virt_first, virt_last)
      INTEGER, DIMENSION(2), INTENT(IN)                  :: band_range
      INTEGER, INTENT(IN)                                :: homo, nmo
      INTEGER, INTENT(OUT)                               :: occ_first, occ_last, virt_first, &
                                                            virt_last

      IF (ALL(band_range == -1)) THEN
         occ_first = 1
         occ_last = homo
         virt_first = homo + 1
         virt_last = nmo
      ELSE
         IF (ANY(band_range < 1)) THEN
            CPABORT("ECD| IPA TRANSITION_BAND_RANGE must be -1 -1 or positive bounds")
         END IF
         IF (band_range(2) > nmo) THEN
            CPABORT("ECD| IPA TRANSITION_BAND_RANGE exceeds available MOs")
         END IF
         occ_first = band_range(1)
         occ_last = MIN(homo, band_range(2))
         virt_first = MAX(homo + 1, band_range(1))
         virt_last = band_range(2)
      END IF

      IF (occ_first > occ_last .OR. virt_first > virt_last) THEN
         CPABORT("ECD| IPA TRANSITION_BAND_RANGE excludes occupied-virtual pairs")
      END IF

   END SUBROUTINE ecd_decode_transition_range

! **************************************************************************************************
!> \brief Build direct length-form and orbital magnetic AO operators for PERIODIC NONE.
!> \param qs_env Quickstep environment
!> \param mo_coeff MO coefficient matrix
!> \param nao Number of AOs
!> \param nmo Number of MOs
!> \param pseudopotential_commutator Include the complete GTH non-local magnetic correction
!> \param rcc Common reference point for all molecular operators
!> \param position_mo Three direct MO-space position matrices
!> \param magnetic_mo Three MO-space angular-momentum matrices
! **************************************************************************************************
   SUBROUTINE ecd_build_molecular_operators_none(qs_env, mo_coeff, nao, nmo, &
                                                  pseudopotential_commutator, rcc, &
                                                  position_mo, magnetic_mo)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER, INTENT(IN)                                :: nao, nmo
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: rcc
      TYPE(cp_fm_type), DIMENSION(ncart), INTENT(INOUT)  :: magnetic_mo, position_mo

      INTEGER                                            :: idir
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_ao, position_ao

      NULLIFY (magnetic_ao, position_ao)
      CALL ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, rcc, &
                                                  position_ao, magnetic_ao)

      DO idir = 1, ncart
         CALL ecd_project_operator_to_mo(position_ao(idir), mo_coeff, nao, nmo, position_mo(idir))
         CALL ecd_project_operator_to_mo(magnetic_ao(idir), mo_coeff, nao, nmo, magnetic_mo(idir))
      END DO
      CALL dbcsr_deallocate_matrix_set(magnetic_ao)
      CALL dbcsr_deallocate_matrix_set(position_ao)

   END SUBROUTINE ecd_build_molecular_operators_none

! **************************************************************************************************
!> \brief Build direct length-form and orbital magnetic AO operators for PERIODIC NONE.
!> \param qs_env Quickstep environment
!> \param pseudopotential_commutator Include the complete GTH non-local magnetic correction
!> \param rcc Common reference point for all molecular operators
!> \param position_ao Three direct AO-space position matrices
!> \param magnetic_ao Three direct AO-space angular-momentum matrices
! **************************************************************************************************
   SUBROUTINE ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, rcc, &
                                                     position_ao, magnetic_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: rcc
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

      CALL build_local_moment_matrix(qs_env, position_ao, 1, ref_point=rcc)
      CALL build_local_magmom_matrix(qs_env, magnetic_ao, 1, ref_point=rcc)

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
         ! the velocity-defined orbital magnetic dipole.
         CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                               particle_set, cell, matrix_rxrv=magnetic_nl, ref_point=rcc)
         DO idir = 1, ncart
            CALL dbcsr_add(magnetic_ao(idir)%matrix, magnetic_nl(idir)%matrix, 1.0_dp, 1.0_dp)
         END DO
         CALL dbcsr_deallocate_matrix_set(magnetic_nl)
      END IF

   END SUBROUTINE ecd_build_molecular_operators_ao_none

! **************************************************************************************************
!> \brief Project a sparse AO operator to the full MO space.
!> \param operator_ao AO-space sparse operator
!> \param mo_coeff MO coefficient matrix
!> \param nao Number of AOs
!> \param nmo Number of MOs
!> \param operator_mo Full MO-space operator
! **************************************************************************************************
   SUBROUTINE ecd_project_operator_to_mo(operator_ao, mo_coeff, nao, nmo, operator_mo)
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

   END SUBROUTINE ecd_project_operator_to_mo

! **************************************************************************************************
!> \brief Assemble direct molecular ECD transition quantities.
!> \param position_mo Three direct MO-space position matrices
!> \param magnetic_mo Three direct MO-space angular-momentum matrices
!> \param eigenvalues MO eigenvalues
!> \param iocc Occupied MO index
!> \param ivirt Virtual MO index
!> \param transition_weight Occupation-number difference
!> \param position Electric position transition dipole
!> \param magnetic Real coefficient of the orbital magnetic transition dipole
!> \param rotatory_strength Isotropic IPA rotatory strength
!> \param oscillator_strength Electric-dipole oscillator strength
! **************************************************************************************************
   SUBROUTINE ecd_transition_quantities(position_mo, magnetic_mo, eigenvalues, iocc, ivirt, &
                                        transition_weight, position, magnetic, &
                                        rotatory_strength, oscillator_strength)
      TYPE(cp_fm_type), DIMENSION(ncart), INTENT(IN)     :: position_mo, magnetic_mo
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      INTEGER, INTENT(IN)                                :: iocc, ivirt
      REAL(KIND=dp), INTENT(IN)                          :: transition_weight
      REAL(KIND=dp), DIMENSION(ncart), INTENT(OUT)       :: magnetic, position
      REAL(KIND=dp), INTENT(OUT)                         :: oscillator_strength, &
                                                            rotatory_strength

      INTEGER                                            :: idir
      REAL(KIND=dp)                                      :: gap

      gap = eigenvalues(ivirt) - eigenvalues(iocc)
      position = 0.0_dp
      magnetic = 0.0_dp
      DO idir = 1, ncart
         CALL cp_fm_get_element(position_mo(idir), ivirt, iocc, position(idir))
         CALL cp_fm_get_element(magnetic_mo(idir), ivirt, iocc, magnetic(idir))
         magnetic(idir) = 0.5_dp*magnetic(idir)
      END DO

      ! mu_ic = -position and m_ci = i*magnetic for real molecular orbitals.
      rotatory_strength = -transition_weight*SUM(position*magnetic)
      oscillator_strength = 2.0_dp/3.0_dp*gap*transition_weight*SUM(position**2)

   END SUBROUTINE ecd_transition_quantities

! **************************************************************************************************
!> \brief Print the IPA ECD stick-table header.
!> \param iounit Output unit
!> \param pseudopotential_commutator GTH non-local magnetic correction switch
! **************************************************************************************************
   SUBROUTINE ecd_print_ipa_header(iounit, pseudopotential_commutator)
      INTEGER, INTENT(IN)                                :: iounit
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator

      WRITE (iounit, "(1X,A)") "", &
         "-------------------------------------------------------------------------------", &
         "-                  ECD IPA transition candidates                              -", &
         "-------------------------------------------------------------------------------"
      WRITE (iounit, "(1X,A)") &
         "ECD| PERIODIC NONE IPA sticks use Gamma-only MO transitions."
      WRITE (iounit, "(1X,A,L1)") &
         "ECD| GTH non-local magnetic correction included: ", &
         pseudopotential_commutator
      WRITE (iounit, "(1X,A)") &
         "ECD| ECD_IPA_R is the direct length-form position matrix element."
      WRITE (iounit, "(1X,A)") &
         "ECD| ECD_IPA_M is the real coefficient of the physical orbital magnetic dipole."
      WRITE (iounit, "(1X,A)") &
         "ECD| Magnetic matrix uses direct local angular and GTH non-local AO operators; spin is omitted."
      WRITE (iounit, '(/,T8,A,T17,A,T29,A,T40,A,T50,A,T60,A,T70,A,T84,A,T99,A)') &
         "State", "Energy", "k-point", "Spin", "Occ", "Virt", "Weight", "Rotatory", &
         "Oscillator"
      WRITE (iounit, '(T8,A,T17,A,T29,A,T40,A,T50,A,T60,A,T70,A,T84,A,T99,A)') &
         "number", "eV", "index", "index", "band", "band", "factor", "strength", &
         "strength"
      WRITE (iounit, '(T8,104("-"))')

   END SUBROUTINE ecd_print_ipa_header

! **************************************************************************************************
!> \brief Print one IPA ECD transition and its vector matrix elements.
! **************************************************************************************************
   SUBROUTINE ecd_print_ipa_transition(iounit, istate, energy, ispin, iocc, ivirt, &
                                       transition_weight, position, magnetic, &
                                       rotatory_strength, oscillator_strength)
      INTEGER, INTENT(IN)                                :: iocc, iounit, ispin, istate, ivirt
      REAL(KIND=dp), INTENT(IN)                          :: energy, oscillator_strength, &
                                                            rotatory_strength, transition_weight
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: magnetic, position

      WRITE (iounit, '(1X,A,T11,I7,T20,F11.5,T32,I7,T43,I7,T53,I7,T63,I7,'// &
                      'T73,F8.3,T84,ES12.5E2,T99,ES12.5E2)') &
         "ECD_IPA|", istate, energy*evolt, 1, ispin, iocc, ivirt, transition_weight, &
         rotatory_strength, oscillator_strength
      WRITE (iounit, '(1X,A,T17,3(1X,ES13.6E2))') "ECD_IPA_R|", position
      WRITE (iounit, '(1X,A,T17,3(1X,ES13.6E2))') "ECD_IPA_M|", magnetic

   END SUBROUTINE ecd_print_ipa_transition

! **************************************************************************************************
!> \brief Print broadened IPA ECD and absorption spectra on an energy grid.
!> \param ecd_section ECD input section
!> \param transition_energy Sorted transition energies in hartree
!> \param rotatory_stick Rotatory strengths in atomic units
!> \param oscillator_stick Oscillator strengths
!> \param broadening Broadening kernel selector
!> \param fwhm Full width at half maximum in hartree
! **************************************************************************************************
   SUBROUTINE ecd_print_ipa_spectrum(ecd_section, transition_energy, rotatory_stick, &
                                     oscillator_stick, broadening, fwhm)
      TYPE(section_vals_type), POINTER                   :: ecd_section
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: oscillator_stick, rotatory_stick, &
                                                            transition_energy
      INTEGER, INTENT(IN)                                :: broadening
      REAL(KIND=dp), INTENT(IN)                          :: fwhm

      CHARACTER(LEN=default_path_length)                 :: spectrum_filename
      CHARACTER(LEN=default_string_length)               :: x_column_label, x_unit_label, &
                                                            y_column_label, y_unit_label
      INTEGER                                            :: ipoint, istate, spectrum_unit
      INTEGER                                            :: iounit, x_unit, y_unit
      REAL(KIND=dp)                                      :: absorption_density, delta_e, emax, emin, &
                                                            energy, energy_ev, jacobian, margin, &
                                                            rotatory_density, step, x_value, y_scale
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: spectrum_print_key

      IF (broadening == ecd_broadening_none) RETURN
      IF (fwhm <= 0.0_dp) THEN
         CPABORT("ECD| FWHM must be positive when BROADENING is active")
      END IF
      CPASSERT(SIZE(transition_energy) == SIZE(rotatory_stick))
      CPASSERT(SIZE(transition_energy) == SIZE(oscillator_stick))

      margin = 8.0_dp*fwhm
      emin = MAX(0.0_dp, MINVAL(transition_energy) - margin)
      emax = MAXVAL(transition_energy) + margin
      IF (emax <= emin) emax = emin + fwhm

      logger => cp_get_default_logger()
      spectrum_print_key => section_vals_get_subs_vals(ecd_section, "PRINT%SPECTRUM")
      CALL section_vals_val_get(ecd_section, "PRINT%SPECTRUM%X_UNIT", i_val=x_unit)
      CALL section_vals_val_get(ecd_section, "PRINT%SPECTRUM%Y_UNIT", i_val=y_unit)

      SELECT CASE (x_unit)
      CASE (ecd_spectrum_x_ev)
         x_unit_label = "eV"
         x_column_label = "Energy_eV"
      CASE (ecd_spectrum_x_nm)
         x_unit_label = "nm"
         x_column_label = "Wavelength_nm"
      CASE (ecd_spectrum_x_cm1)
         x_unit_label = "cm-1"
         x_column_label = "Wavenumber_cm-1"
      CASE DEFAULT
         CPABORT("ECD| Unknown SPECTRUM X_UNIT")
      END SELECT

      SELECT CASE (y_unit)
      CASE (ecd_spectrum_y_au)
         y_scale = 1.0_dp
         y_unit_label = "au"
         y_column_label = "ECD_R_au"
      CASE (ecd_spectrum_y_cgs_1e40)
         y_scale = ecd_rot_au_to_cgs_1e40
         y_unit_label = "1e-40_cgs"
         y_column_label = "ECD_R_1e-40_cgs"
      CASE DEFAULT
         CPABORT("ECD| Unknown SPECTRUM Y_UNIT")
      END SELECT

      IF (x_unit == ecd_spectrum_x_nm .AND. emin <= 0.0_dp) THEN
         emin = MAX(0.1_dp*MINVAL(transition_energy), TINY(emin))
      END IF
      step = (emax - emin)/REAL(ecd_spectrum_points - 1, KIND=dp)

      spectrum_filename = cp_print_key_generate_filename(logger, spectrum_print_key, &
                                                         extension=".ecd", my_local=.FALSE.)
      spectrum_unit = cp_print_key_unit_nr(logger, ecd_section, "PRINT%SPECTRUM", &
                                           extension=".ecd", file_position="REWIND", &
                                           file_action="WRITE", file_form="FORMATTED")
      IF (spectrum_unit <= 0) RETURN

      iounit = cp_logger_get_default_io_unit(logger)
      IF (iounit > 0) THEN
         WRITE (iounit, '(T2,A,1X,A)') "ECD| Broadened spectrum file", &
            TRIM(spectrum_filename)
      END IF

      WRITE (spectrum_unit, '(A)') "# CP2K ECD broadened spectrum"
      SELECT CASE (broadening)
      CASE (ecd_broadening_gaussian)
         WRITE (spectrum_unit, '(A,1X,A)') "# Broadening", "GAUSSIAN"
      CASE (ecd_broadening_lorentzian)
         WRITE (spectrum_unit, '(A,1X,A)') "# Broadening", "LORENTZIAN"
      END SELECT
      WRITE (spectrum_unit, '(A,1X,ES16.8E3,1X,A,1X,F12.6)') "# FWHM", fwhm, "hartree", &
         fwhm*evolt
      WRITE (spectrum_unit, '(A,1X,A)') "# X unit", TRIM(x_unit_label)
      WRITE (spectrum_unit, '(A,1X,A)') "# ECD R unit", TRIM(y_unit_label)
      WRITE (spectrum_unit, '(A,1X,ES16.8E3)') "# R au to 1e-40 cgs", &
         ecd_rot_au_to_cgs_1e40
      WRITE (spectrum_unit, '(A,1X,A,1X,A,1X,A)') "# Columns:", TRIM(x_column_label), &
         TRIM(y_column_label)//"_per_"//TRIM(x_unit_label), &
         "Oscillator_per_"//TRIM(x_unit_label)

      DO ipoint = 1, ecd_spectrum_points
         IF (x_unit == ecd_spectrum_x_nm) THEN
            energy = emax - REAL(ipoint - 1, KIND=dp)*step
         ELSE
            energy = emin + REAL(ipoint - 1, KIND=dp)*step
         END IF
         rotatory_density = 0.0_dp
         absorption_density = 0.0_dp
         DO istate = 1, SIZE(transition_energy)
            delta_e = energy - transition_energy(istate)
            rotatory_density = rotatory_density + rotatory_stick(istate)* &
                               ecd_broadening_function(delta_e, broadening, fwhm)
            absorption_density = absorption_density + oscillator_stick(istate)* &
                                 ecd_broadening_function(delta_e, broadening, fwhm)
         END DO
         energy_ev = energy*evolt
         SELECT CASE (x_unit)
         CASE (ecd_spectrum_x_ev)
            x_value = energy_ev
            jacobian = 1.0_dp/evolt
         CASE (ecd_spectrum_x_nm)
            x_value = ecd_hc_ev_nm/energy_ev
            jacobian = ecd_hc_ev_nm/(evolt*x_value**2)
         CASE (ecd_spectrum_x_cm1)
            x_value = energy_ev*ecd_ev_to_cm1
            jacobian = 1.0_dp/(evolt*ecd_ev_to_cm1)
         END SELECT
         WRITE (spectrum_unit, '(3(1X,ES18.10E3))') x_value, &
            y_scale*rotatory_density*jacobian, absorption_density*jacobian
      END DO

      CALL cp_print_key_finished_output(spectrum_unit, logger, ecd_section, "PRINT%SPECTRUM")

   END SUBROUTINE ecd_print_ipa_spectrum

! **************************************************************************************************
!> \brief Normalized ECD broadening kernel. The width is FWHM in hartree.
!> \param delta_e Energy offset in hartree
!> \param broadening Broadening kernel selector
!> \param fwhm Full width at half maximum in hartree
!> \return Kernel value in inverse hartree
! **************************************************************************************************
   PURE FUNCTION ecd_broadening_function(delta_e, broadening, fwhm) RESULT(value)
      REAL(KIND=dp), INTENT(IN)                          :: delta_e, fwhm
      INTEGER, INTENT(IN)                                :: broadening
      REAL(KIND=dp)                                      :: value

      REAL(KIND=dp)                                      :: gamma, sigma

      SELECT CASE (broadening)
      CASE (ecd_broadening_gaussian)
         sigma = fwhm/(2.0_dp*SQRT(2.0_dp*LOG(2.0_dp)))
         value = EXP(-0.5_dp*(delta_e/sigma)**2)/(sigma*SQRT(2.0_dp*pi))
      CASE (ecd_broadening_lorentzian)
         gamma = 0.5_dp*fwhm
         value = gamma/(pi*(delta_e**2 + gamma**2))
      CASE DEFAULT
         value = 0.0_dp
      END SELECT

   END FUNCTION ecd_broadening_function

! **************************************************************************************************
!> \brief Print IPA ECD stick-table checksums.
! **************************************************************************************************
   SUBROUTINE ecd_print_ipa_footer(iounit, energy_checksum, oscillator_checksum, &
                                   rotatory_checksum)
      INTEGER, INTENT(IN)                                :: iounit
      REAL(KIND=dp), INTENT(IN)                          :: energy_checksum, oscillator_checksum, &
                                                            rotatory_checksum

      WRITE (iounit, '(/,T2,A,E16.8)') 'ECD IPA : CheckSum E = ', energy_checksum
      WRITE (iounit, '(T2,A,E16.8)') 'ECD IPA : CheckSum F = ', oscillator_checksum
      WRITE (iounit, '(T2,A,E16.8)') 'ECD IPA : CheckSum R = ', rotatory_checksum
      WRITE (iounit, "(1X,A)") &
         "-------------------------------------------------------------------------------"

   END SUBROUTINE ecd_print_ipa_footer

END MODULE qs_ecd_methods
