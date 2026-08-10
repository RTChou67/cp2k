!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program for atomistic simulations                                              !
!   Copyright 2000-2026 CP2K developers group                                                      !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gamma-point PERIODIC XYZ IPA ECD spectra.
! **************************************************************************************************
MODULE qs_ecd_spectra_xyz
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_generate_filename,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_constants,                 ONLY: tddfpt_dipole_velocity
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE physcon,                         ONLY: evolt
   USE qs_ecd_operator_utils,           ONLY: ecd_build_e1_m1_mo,&
                                              ecd_mo_operators_type,&
                                              ecd_open_transition_files,&
                                              ecd_write_transition_operators
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PRIVATE       :: ncart = 3
   REAL(KIND=dp), PARAMETER, PRIVATE :: occupation_threshold = 1.0E-8_dp

   PUBLIC :: ecd_spectra_xyz_ipa

CONTAINS

! **************************************************************************************************
!> \brief Assemble the Gamma-point periodic E1--M1 tensor from complete MO operators.
! **************************************************************************************************
   SUBROUTINE ecd_spectra_xyz_ipa(qs_env, pseudopotential_commutator, reference_point)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: pseudopotential_commutator
      REAL(KIND=dp), DIMENSION(ncart), INTENT(IN)        :: reference_point

      CHARACTER(LEN=default_path_length)                 :: e1_filename, m1_filename, tensor_filename
      COMPLEX(KIND=dp), DIMENSION(ncart)                 :: e1_ai, m1_ai
      INTEGER                                            :: e1_unit, idir, iocc, iounit, ispin, &
                                                            itransition, ivirt, jdir, m1_unit, &
                                                            nmo, nspin, ntransition, output_index, &
                                                            tensor_unit
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: occupied, spin_channel, &
                                                            transition_order, virtual
      LOGICAL                                            :: write_transition_matrices
      REAL(KIND=dp)                                      :: cell_volume, e1_m1_imag, gap, &
                                                            occupation_difference
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: transition_energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: spectral_tensor
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupations
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(ecd_mo_operators_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: operators
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      NULLIFY (cell, eigenvalues, logger, mos, occupations)
      CALL get_qs_env(qs_env, cell=cell, mos=mos)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(mos))
      CALL get_cell(cell, deth=cell_volume)
      cell_volume = ABS(cell_volume)
      IF (cell_volume <= 0.0_dp) CPABORT("ECD| Periodic E1/M1 tensor requires a positive cell volume")

      nspin = SIZE(mos)
      ntransition = 0
      DO ispin = 1, nspin
         NULLIFY (eigenvalues, occupations)
         CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, nmo=nmo, occupation_numbers=occupations)
         IF (nmo < 1 .OR. .NOT. ASSOCIATED(eigenvalues) .OR. .NOT. ASSOCIATED(occupations)) THEN
            CPABORT("ECD| Periodic velocity E1/M1 output received incomplete MO data")
         END IF
         DO iocc = 1, nmo
            IF (occupations(iocc) <= occupation_threshold) CYCLE
            DO ivirt = 1, nmo
               IF (occupations(ivirt) > occupation_threshold) CYCLE
               gap = eigenvalues(ivirt) - eigenvalues(iocc)
               occupation_difference = occupations(iocc) - occupations(ivirt)
               IF (gap > 0.0_dp .AND. occupation_difference > 0.0_dp) ntransition = ntransition + 1
            END DO
         END DO
      END DO
      IF (ntransition < 1) CPABORT("ECD| Periodic velocity E1/M1 output found no occupied-empty transitions")

      ALLOCATE (occupied(ntransition), spin_channel(ntransition), spectral_tensor(ncart, ncart, ntransition), &
                transition_energy(ntransition), transition_order(ntransition), virtual(ntransition))
      spectral_tensor = 0.0_dp

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      CALL cp_logger_generate_filename(logger, tensor_filename, logger%iter_info%project_name, &
                                       "-periodic-ipa-e1m1-tensor.dat", local=.FALSE.)
      CALL ecd_open_transition_files(logger, "PERIODIC XYZ Gamma velocity-gauge IPA", write_transition_matrices, &
                                     e1_unit, m1_unit, e1_filename, m1_filename)

      CALL ecd_build_e1_m1_mo(qs_env, tddfpt_dipole_velocity, pseudopotential_commutator, reference_point, operators)
      IF (SIZE(operators) /= nspin) CPABORT("ECD| Periodic E1/M1 spin dimensions are inconsistent")

      itransition = 0
      DO ispin = 1, nspin
         NULLIFY (eigenvalues, occupations)
         CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, nmo=nmo, occupation_numbers=occupations)
         DO iocc = 1, nmo
            IF (occupations(iocc) <= occupation_threshold) CYCLE
            DO ivirt = 1, nmo
               IF (occupations(ivirt) > occupation_threshold) CYCLE
               gap = eigenvalues(ivirt) - eigenvalues(iocc)
               occupation_difference = occupations(iocc) - occupations(ivirt)
               IF (gap <= 0.0_dp .OR. occupation_difference <= 0.0_dp) CYCLE
               itransition = itransition + 1
               transition_energy(itransition) = gap
               occupied(itransition) = iocc
               spin_channel(itransition) = ispin
               virtual(itransition) = ivirt

               e1_ai = operators(ispin)%e1(:, ivirt, iocc)
               m1_ai = operators(ispin)%m1(:, ivirt, iocc)
               e1_m1_imag = AIMAG(SUM(CONJG(e1_ai)*m1_ai))
               DO idir = 1, ncart
                  DO jdir = 1, ncart
                     spectral_tensor(idir, jdir, itransition) = 4.0_dp*pi*pi*occupation_difference/ &
                                                                (gap*cell_volume)*( &
                                                                2.0_dp*MERGE(e1_m1_imag, 0.0_dp, idir == jdir) - &
                                                                AIMAG(CONJG(e1_ai(idir))*m1_ai(jdir) + &
                                                                      CONJG(e1_ai(jdir))*m1_ai(idir)))
                  END DO
               END DO
               IF (write_transition_matrices) THEN
                  CALL ecd_write_transition_operators( &
                     e1_unit, m1_unit, itransition, gap, ispin, iocc, ivirt, occupation_difference, &
                     operators(ispin)%e1(:, ivirt, iocc), operators(ispin)%m1(:, ivirt, iocc))
               END IF
            END DO
         END DO
      END DO
      IF (itransition /= ntransition) CPABORT("ECD| Periodic velocity E1/M1 transition count is inconsistent")
      DEALLOCATE (operators)
      IF (write_transition_matrices) THEN
         CALL close_file(e1_unit)
         CALL close_file(m1_unit)
      END IF

      CALL sort(transition_energy, ntransition, transition_order)
      IF (write_transition_matrices) THEN
         CALL open_file(tensor_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=tensor_unit)
         WRITE (tensor_unit, '(A)') "# Periodic Gamma IPA CD tensor coefficient"
         WRITE (tensor_unit, '(A)') "# E1_M1_ONLY: D = D_M1 and D_E2 = 0"
         WRITE (tensor_unit, '(A)') "# Delta_alpha(mu,nu;omega) = sum_t coefficient(t,mu,nu)*delta(gap_t-omega)"
         WRITE (tensor_unit, '(A)') "# state energy_eV spin occupied virtual alpha beta eq17_coefficient_au"
         DO output_index = 1, ntransition
            itransition = transition_order(output_index)
            DO idir = 1, ncart
               DO jdir = 1, ncart
                  WRITE (tensor_unit, '(I0,1X,ES24.15E3,5(1X,I0),1X,ES24.15E3)') &
                     output_index, transition_energy(output_index)*evolt, spin_channel(itransition), &
                     occupied(itransition), virtual(itransition), idir, jdir, &
                     spectral_tensor(idir, jdir, itransition)
               END DO
            END DO
         END DO
         CALL close_file(tensor_unit)
      END IF

      DEALLOCATE (occupied, spectral_tensor, spin_channel, transition_energy, transition_order, virtual)
      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T2,A)') "ECD| Periodic Gamma E1 velocity matrix written to "//TRIM(e1_filename)
         WRITE (iounit, '(T2,A)') "ECD| Periodic Gamma M1 velocity matrix written to "//TRIM(m1_filename)
         WRITE (iounit, '(T2,A)') "ECD| Periodic Gamma E1--M1 tensor written to "//TRIM(tensor_filename)
         WRITE (iounit, '(T2,A,I0)') "ECD| Transition count ", ntransition
      END IF

   END SUBROUTINE ecd_spectra_xyz_ipa

END MODULE qs_ecd_spectra_xyz
