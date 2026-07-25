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
   USE cp_fm_struct, ONLY: cp_fm_struct_create, &
                           cp_fm_struct_type
   USE cp_fm_types, ONLY: cp_fm_create, &
                          cp_fm_get_element, &
                          cp_fm_get_info, &
                          cp_fm_release, &
                          cp_fm_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE input_section_types, ONLY: section_vals_get, &
                                  section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE input_constants, ONLY: no_sf_tddfpt, &
                              tddfpt_kernel_full, &
                              tddfpt_kernel_none, &
                              use_mom_ref_user
   USE kinds, ONLY: dp
   USE moments_utils, ONLY: get_reference_point
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE physcon, ONLY: evolt
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_kind_types, ONLY: qs_kind_type
   USE qs_mo_types, ONLY: get_mo_set, &
                          mo_set_type
   USE qs_moments, ONLY: build_local_magmom_matrix, &
                         build_local_moment_matrix
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type
   USE qs_tddfpt2_types, ONLY: tddfpt_ground_state_mos
   USE particle_types, ONLY: particle_type
   USE scf_control_types, ONLY: scf_control_type
   USE util, ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_methods'
   INTEGER, PARAMETER, PRIVATE          :: ncart = 3

   PUBLIC :: ecd_spectra

CONTAINS

! **************************************************************************************************
!> \brief Assemble molecular ECD sticks from direct transitions or TDDFPT response amplitudes.
!> \param qs_env Quickstep environment
! **************************************************************************************************
   SUBROUTINE ecd_spectra(qs_env, evects, evals, gs_mos, tddfpt_kernel, mult, do_soc, spinflip)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:, :), OPTIONAL, INTENT(IN) :: evects
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, INTENT(IN)  :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), POINTER, OPTIONAL :: gs_mos
      INTEGER, OPTIONAL, INTENT(IN)                      :: mult, spinflip, tddfpt_kernel
      LOGICAL, OPTIONAL, INTENT(IN)                      :: do_soc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ecd_spectra'

      INTEGER                                            :: handle, homo, idir, iocc, ispin, istate, itrans, ivirt, &
                                                            nao, nactive, nmo, nspins, nstates, ntrans_total, reference
      LOGICAL                                            :: do_ecd, do_hfx, pseudopotential_commutator, tddfpt_enabled, &
                                                            use_direct_transitions
      LOGICAL                                            :: ref_point_explicit
      REAL(KIND=dp), DIMENSION(ncart)                    :: reference_point
      REAL(KIND=dp), DIMENSION(ncart)                    :: transition_magnetic, transition_position
      REAL(KIND=dp)                                      :: transition_weight
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: transition_order
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE     :: magnetic_spin, position_spin
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: rotatory_strength, transition_energy
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: magnetic_active_action, position_active_action
      TYPE(cp_fm_type), DIMENSION(ncart)                 :: magnetic_mo, position_mo
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_ao, position_ao
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: ecd_section, hfx_section, tddfpt_section

      CALL timeset(routineN, handle)

      NULLIFY (ecd_section, eigenvalues, hfx_section, magnetic_ao, mo_coeff, mos, occupation_numbers, &
               position_ao, ref_point, scf_control, tddfpt_section)
      ecd_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ECD")

      CALL section_vals_val_get(ecd_section, "DO_ECD", l_val=do_ecd)
      IF (.NOT. do_ecd) THEN
         CALL timestop(handle)
         RETURN
      END IF

      use_direct_transitions = .TRUE.
      IF (PRESENT(evects)) THEN
         IF (.NOT. (PRESENT(evals) .AND. PRESENT(gs_mos) .AND. PRESENT(tddfpt_kernel) .AND. &
                    PRESENT(mult) .AND. PRESENT(do_soc) .AND. PRESENT(spinflip))) THEN
            CPABORT("ECD| TDDFPT callback requires complete response data")
         END IF
         IF (do_soc) CPABORT("ECD| TDDFPT ECD does not support TDDFPT SOC")
         IF (spinflip /= no_sf_tddfpt) CPABORT("ECD| TDDFPT ECD does not support spin-flip TDDFPT")
         ! FULL supplies coupled response amplitudes; NONE uses direct occupied-virtual transitions.
         SELECT CASE (tddfpt_kernel)
         CASE (tddfpt_kernel_full)
            use_direct_transitions = .FALSE.
         CASE (tddfpt_kernel_none)
            use_direct_transitions = .TRUE.
         CASE DEFAULT
            CPABORT("ECD| TDDFPT ECD supports TDDFPT KERNEL FULL and KERNEL NONE")
         END SELECT
      END IF

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

      CALL ecd_require_nonperiodic(qs_env)

      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CPABORT("ECD| Molecular ECD currently does not support HFX")
      END IF

      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")
      CALL section_vals_get(tddfpt_section, explicit=tddfpt_enabled)
      ! The TDDFPT callback emits ECD while its response data are available.
      IF (.NOT. PRESENT(evects) .AND. tddfpt_enabled) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (use_direct_transitions) THEN
         ! IPA and TDDFPT KERNEL NONE use individual occupied-virtual transitions.
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
         itrans = 0
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, homo=homo, mo_coeff=mo_coeff, nao=nao, &
                            nmo=nmo, occupation_numbers=occupation_numbers)
            CPASSERT(ASSOCIATED(eigenvalues))
            CPASSERT(ASSOCIATED(mo_coeff))
            CPASSERT(ASSOCIATED(occupation_numbers))
            IF (homo <= 0 .OR. nmo <= homo) CYCLE

            CALL ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, reference_point, &
                                                       position_ao, magnetic_ao)
            DO idir = 1, ncart
               CALL ecd_project_operator_to_mo(position_ao(idir), mo_coeff, position_mo(idir))
               CALL ecd_project_operator_to_mo(magnetic_ao(idir), mo_coeff, magnetic_mo(idir))
            END DO
            DO iocc = 1, homo
               DO ivirt = homo + 1, nmo
                  itrans = itrans + 1
                  transition_energy(itrans) = eigenvalues(ivirt) - eigenvalues(iocc)
                  IF (transition_energy(itrans) <= 0.0_dp) CPABORT("ECD| IPA requires positive occupied-virtual gaps")
                  transition_weight = occupation_numbers(iocc) - occupation_numbers(ivirt)
                  IF (transition_weight <= 0.0_dp) THEN
                     CPABORT("ECD| IPA requires a positive occupied-virtual occupation difference")
                  END IF
                  DO idir = 1, ncart
                     CALL cp_fm_get_element(position_mo(idir), ivirt, iocc, transition_position(idir))
                     CALL cp_fm_get_element(magnetic_mo(idir), ivirt, iocc, transition_magnetic(idir))
                  END DO
                  CALL ecd_transition_quantities(transition_position, transition_magnetic, transition_weight, &
                                                 rotatory_strength(itrans))
               END DO
            END DO
            DO idir = 1, ncart
               CALL cp_fm_release(magnetic_mo(idir))
               CALL cp_fm_release(position_mo(idir))
            END DO
            CALL dbcsr_deallocate_matrix_set(magnetic_ao)
            CALL dbcsr_deallocate_matrix_set(position_ao)
         END DO
         CALL sort(transition_energy, ntrans_total, transition_order)
         CALL ecd_print_sticks(transition_energy, rotatory_strength, transition_order)
         DEALLOCATE (transition_order, rotatory_strength, transition_energy)
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
         CALL ecd_build_tddfpt_operator_actions(qs_env, gs_mos, pseudopotential_commutator, reference_point, &
                                                position_active_action, magnetic_active_action)
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
            CALL ecd_transition_quantities(transition_position, transition_magnetic, 1.0_dp, &
                                           rotatory_strength(istate))
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
!> \brief Form direct AO operator actions on the active occupied TDDFPT ground-state orbitals.
!> \param qs_env Quickstep environment
!> \param gs_mos Ground-state active occupied MO data
!> \param pseudopotential_commutator Include GTH non-local magnetic correction
!> \param reference_point Common molecular reference point
!> \param position_active_action Position operator actions X*C_active [direction, spin]
!> \param magnetic_active_action Angular-momentum operator actions Lambda*C_active [direction, spin]
! **************************************************************************************************
   SUBROUTINE ecd_build_tddfpt_operator_actions(qs_env, gs_mos, pseudopotential_commutator, reference_point, &
                                                position_active_action, magnetic_active_action)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: reference_point
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(INOUT)   :: magnetic_active_action, position_active_action

      INTEGER                                            :: idir, ispin, nactive
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magnetic_ao, position_ao

      NULLIFY (fm_struct, magnetic_ao, position_ao)
      CALL ecd_build_molecular_operators_ao_none(qs_env, pseudopotential_commutator, reference_point, &
                                                 position_ao, magnetic_ao)
      DO ispin = 1, SIZE(gs_mos)
         CPASSERT(ASSOCIATED(gs_mos(ispin)%mos_active))
         CALL cp_fm_get_info(gs_mos(ispin)%mos_active, matrix_struct=fm_struct, ncol_global=nactive)
         DO idir = 1, ncart
            CALL cp_fm_create(position_active_action(idir, ispin), fm_struct)
            CALL cp_fm_create(magnetic_active_action(idir, ispin), fm_struct)
            CALL cp_dbcsr_sm_fm_multiply(position_ao(idir)%matrix, gs_mos(ispin)%mos_active, &
                                         position_active_action(idir, ispin), ncol=nactive, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(magnetic_ao(idir)%matrix, gs_mos(ispin)%mos_active, &
                                         magnetic_active_action(idir, ispin), ncol=nactive, alpha=1.0_dp, beta=0.0_dp)
         END DO
         NULLIFY (fm_struct)
      END DO
      CALL dbcsr_deallocate_matrix_set(magnetic_ao)
      CALL dbcsr_deallocate_matrix_set(position_ao)

   END SUBROUTINE ecd_build_tddfpt_operator_actions

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
!> \brief Build direct length-form and orbital magnetic AO operators for PERIODIC NONE.
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
         ! the velocity-defined orbital magnetic dipole.
         CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                               particle_set, cell, matrix_rxrv=magnetic_nl, ref_point=reference_point)
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
!> \param operator_mo Full MO-space operator
! **************************************************************************************************
   SUBROUTINE ecd_project_operator_to_mo(operator_ao, mo_coeff, operator_mo)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: operator_ao
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_fm_type), INTENT(INOUT)                    :: operator_mo

      INTEGER                                            :: nao, nmo
      TYPE(cp_fm_struct_type), POINTER                   :: coeff_struct
      TYPE(cp_fm_type)                                   :: operator_coeff

      NULLIFY (coeff_struct)

      CALL cp_fm_get_info(mo_coeff, matrix_struct=coeff_struct, nrow_global=nao, ncol_global=nmo)
      CALL cp_fm_create(operator_coeff, coeff_struct, name="ECD operator C")
      CALL cp_fm_create(operator_mo, coeff_struct, nrow=nmo, ncol=nmo, name="ECD C^T operator C")

      CALL cp_dbcsr_sm_fm_multiply(operator_ao%matrix, mo_coeff, operator_coeff, ncol=nmo, &
                                   alpha=1.0_dp, beta=0.0_dp)
      CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mo_coeff, operator_coeff, &
                         0.0_dp, operator_mo)

      CALL cp_fm_release(operator_coeff)

   END SUBROUTINE ecd_project_operator_to_mo

! **************************************************************************************************
!> \brief Assemble an isotropic ECD rotatory strength from Cartesian transition components.
!> \param transition_position Position transition components
!> \param transition_magnetic Angular-momentum transition components
!> \param transition_weight Occupation or state weight
!> \param rotatory_strength Isotropic rotatory strength
! **************************************************************************************************
   SUBROUTINE ecd_transition_quantities(transition_position, transition_magnetic, transition_weight, &
                                        rotatory_strength)
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: transition_position, transition_magnetic
      REAL(KIND=dp), INTENT(IN)                          :: transition_weight
      REAL(KIND=dp), INTENT(OUT)                         :: rotatory_strength

      ! mu_ic = -r_ic and m_ci = i*Lambda_ci/2 for real molecular orbitals.
      rotatory_strength = -0.5_dp*transition_weight*SUM(transition_position*transition_magnetic)

   END SUBROUTINE ecd_transition_quantities

END MODULE qs_ecd_methods
