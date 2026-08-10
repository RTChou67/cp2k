!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program for atomistic simulations                                              !
!   Copyright 2000-2026 CP2K developers group                                                     !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gamma-point PERIODIC XYZ IPA ECD spectra.
! **************************************************************************************************
MODULE qs_ecd_spectra_xyz
   USE cell_types, ONLY: cell_type, &
                         get_cell
   USE commutator_rpnl, ONLY: build_com_mom_nl
   USE cp_control_types, ONLY: dft_control_type
   USE cp_dbcsr_api, ONLY: dbcsr_create, &
                           dbcsr_desymmetrize, &
                           dbcsr_finalize, &
                           dbcsr_init_p, &
                           dbcsr_p_type, &
                           dbcsr_reserve_blocks, &
                           dbcsr_release_p, &
                           dbcsr_set, &
                           dbcsr_type_antisymmetric, &
                           dbcsr_type_no_symmetry
   USE cp_files, ONLY: close_file, &
                       open_file
   USE cp_fm_types, ONLY: cp_fm_get_submatrix, &
                          cp_fm_release, &
                          cp_fm_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_generate_filename, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE kinds, ONLY: default_path_length, &
                    dp
   USE mathconstants, ONLY: pi
   USE physcon, ONLY: evolt
   USE qs_ecd_operator_utils, ONLY: ecd_project_operator_to_mo
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_kind_types, ONLY: qs_kind_type
   USE qs_mo_types, ONLY: get_mo_set, &
                          mo_set_type
   USE qs_neighbor_list_types, ONLY: get_iterator_info, &
                                     neighbor_list_iterate, &
                                     neighbor_list_iterator_create, &
                                     neighbor_list_iterator_p_type, &
                                     neighbor_list_iterator_release, &
                                     neighbor_list_set_p_type
   USE qs_operators_ao, ONLY: build_lin_mom_matrix
   USE particle_types, ONLY: particle_type
   USE util, ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PRIVATE :: ncart = 3
   REAL(KIND=dp), PARAMETER, PRIVATE :: occupation_threshold = 1.0E-8_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: energy_degeneracy_threshold = 1.0E-8_dp
   PUBLIC :: ecd_spectra_xyz_ipa

CONTAINS

! **************************************************************************************************
!> \brief Build Gamma E1 and M1 from the complete MO velocity matrix and write transitions.
!>
!> E1 and M1 use the complete MO velocity matrix.  The off-diagonal velocity
!> relation first reconstructs r^V, then m^V=-r^V x v/2.  The E1--M1-only
!> coefficient of the Gamma CD tensor is written for every transition.
!> \param pseudopotential_commutator Include the GTH non-local velocity commutator
! **************************************************************************************************
   SUBROUTINE ecd_spectra_xyz_ipa(qs_env, pseudopotential_commutator)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                 :: pseudopotential_commutator

      INTEGER                                            :: e1_unit, idir, iatom, iocc, inode, iounit, ispin, &
                                                            itransition, ivirt, jatom, jdir, last_jatom, m1_unit, &
                                                            ncommutator_block, nmo, nspin, ntransition, output_index, &
                                                            tensor_unit
      CHARACTER(LEN=default_path_length)                 :: e1_filename, m1_filename, tensor_filename
      LOGICAL                                            :: include_nonlocal_velocity, write_transition_matrices
      REAL(KIND=dp)                                      :: cell_volume, gap, occupation_difference
      COMPLEX(KIND=dp), DIMENSION(ncart)                 :: e1_ai, m1_ai
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: e1_mo, m1_mo, velocity_mo_matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: momentum_matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: commutator_matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: spectral_tensor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: transition_energy
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: commutator_cols, commutator_rows, occupied, spin_channel, &
                                                            transition_order, virtual
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupations
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: commutator_ao, linear_momentum_ao, matrix_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb, sap_ppnl
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: nl_iterator
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cp_fm_type)                                   :: commutator_mo, linear_momentum_mo
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (cell, commutator_ao, dft_control, eigenvalues, linear_momentum_ao, logger, matrix_s, mo_coeff, mos, &
               nl_iterator, occupations, particle_set, qs_kind_set, sab_orb, sap_ppnl)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_s=matrix_s, mos=mos, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set, sab_orb=sab_orb, &
                      sap_ppnl=sap_ppnl)
      IF (.NOT. ASSOCIATED(cell) .OR. .NOT. ASSOCIATED(matrix_s) .OR. SIZE(matrix_s) /= 1 .OR. &
          .NOT. ASSOCIATED(mos)) THEN
         CPABORT("ECD| Periodic velocity E1/M1 output received incomplete Gamma SCF data")
      END IF
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

      ALLOCATE (transition_energy(ntransition), &
                occupied(ntransition), spin_channel(ntransition), transition_order(ntransition), virtual(ntransition))

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      CALL cp_logger_generate_filename(logger, e1_filename, logger%iter_info%project_name, &
                                       "-ecd-e1-transition-matrix.dat", local=.FALSE.)
      CALL cp_logger_generate_filename(logger, m1_filename, logger%iter_info%project_name, &
                                       "-ecd-m1-transition-matrix.dat", local=.FALSE.)
      CALL cp_logger_generate_filename(logger, tensor_filename, logger%iter_info%project_name, &
                                       "-periodic-ipa-e1m1-tensor.dat", local=.FALSE.)
      write_transition_matrices = logger%para_env%is_source()
      IF (write_transition_matrices) THEN
         CALL open_file(e1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=e1_unit)
         WRITE (e1_unit, '(A)') "# Periodic Gamma E1 matrix in the velocity representation"
         WRITE (e1_unit, '(A)') "# E1_ai = i*v_ai/(epsilon_a-epsilon_i)"
         WRITE (e1_unit, '(A)') "# records follow the occupied-empty MO traversal order"
         WRITE (e1_unit, '(A)') "# state energy_eV spin occupied virtual component e1_real e1_imag occupation_difference"
         CALL open_file(m1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=m1_unit)
         WRITE (m1_unit, '(A)') "# Periodic Gamma M1 matrix in the velocity representation"
         WRITE (m1_unit, '(A)') "# M1 = E1 cross v/2"
         WRITE (m1_unit, '(A)') "# records follow the occupied-empty MO traversal order"
         WRITE (m1_unit, '(A)') "# state energy_eV spin occupied virtual component m1_real m1_imag occupation_difference"
      END IF

      ! At Gamma, every real-space image carries unit phase.  Build the folded
      ! local momentum p_lambda^loc,AO directly instead of retaining image matrices.
      ALLOCATE (linear_momentum_ao(ncart))
      DO idir = 1, ncart
         CALL dbcsr_init_p(linear_momentum_ao(idir)%matrix)
         CALL dbcsr_create(linear_momentum_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, linear_momentum_ao(idir)%matrix)
         CALL dbcsr_set(linear_momentum_ao(idir)%matrix, 0.0_dp)
      END DO

      CALL build_lin_mom_matrix(qs_env, linear_momentum_ao)

      ! Keep [r_lambda,V_nl] in an antisymmetric DBCSR matrix.  build_com_mom_nl
      ! fills the canonical atom-pair block; the matrix type supplies its
      ! transposed negative counterpart.
      include_nonlocal_velocity = pseudopotential_commutator .AND. ASSOCIATED(sap_ppnl)
      IF (include_nonlocal_velocity) THEN
         CPASSERT(ASSOCIATED(dft_control))
         CPASSERT(ASSOCIATED(particle_set))
         CPASSERT(ASSOCIATED(qs_kind_set))
         CPASSERT(ASSOCIATED(sab_orb))
         ALLOCATE (commutator_ao(ncart))
         DO idir = 1, ncart
            CALL dbcsr_init_p(commutator_ao(idir)%matrix)
            CALL dbcsr_create(commutator_ao(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
         END DO

         ! Reserve the same canonical atom-pair blocks used by the old
         ! antisymmetric commutator path, using only public DBCSR APIs.
         ncommutator_block = 0
         last_jatom = 0
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, inode=inode)
            IF (inode == 1) last_jatom = 0
            IF (jatom == last_jatom) CYCLE
            last_jatom = jatom
            ncommutator_block = ncommutator_block + 1
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         ALLOCATE (commutator_cols(MAX(1, ncommutator_block)), commutator_rows(MAX(1, ncommutator_block)))
         ncommutator_block = 0
         last_jatom = 0
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, inode=inode)
            IF (inode == 1) last_jatom = 0
            IF (jatom == last_jatom) CYCLE
            last_jatom = jatom
            ncommutator_block = ncommutator_block + 1
            IF (iatom <= jatom) THEN
               commutator_rows(ncommutator_block) = iatom
               commutator_cols(ncommutator_block) = jatom
            ELSE
               commutator_rows(ncommutator_block) = jatom
               commutator_cols(ncommutator_block) = iatom
            END IF
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)
         DO idir = 1, ncart
            CALL dbcsr_finalize(commutator_ao(idir)%matrix)
            CALL dbcsr_reserve_blocks(commutator_ao(idir)%matrix, commutator_rows(1:ncommutator_block), &
                                      commutator_cols(1:ncommutator_block))
            CALL dbcsr_finalize(commutator_ao(idir)%matrix)
            CALL dbcsr_set(commutator_ao(idir)%matrix, 0.0_dp)
         END DO
         DEALLOCATE (commutator_cols, commutator_rows)
         CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                               particle_set, cell, matrix_rv=commutator_ao)
      END IF

      ! cp_fm_get_submatrix replicates each distributed MO block on every rank.
      ! Every rank therefore assembles the same complete E1 and M1 MO matrices.
      ALLOCATE (spectral_tensor(ncart, ncart, ntransition))
      spectral_tensor = 0.0_dp
      itransition = 0
      DO ispin = 1, nspin
         NULLIFY (eigenvalues, mo_coeff, occupations)
         CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, mo_coeff=mo_coeff, nmo=nmo, &
                         occupation_numbers=occupations)
         IF (.NOT. ASSOCIATED(eigenvalues) .OR. .NOT. ASSOCIATED(mo_coeff) .OR. .NOT. ASSOCIATED(occupations)) THEN
            CPABORT("ECD| Periodic velocity E1/M1 output received incomplete Gamma MO data")
         END IF
         ALLOCATE (e1_mo(ncart, nmo, nmo), m1_mo(ncart, nmo, nmo), momentum_matrix(nmo, nmo), &
                   velocity_mo_matrix(ncart, nmo, nmo))
         IF (include_nonlocal_velocity) ALLOCATE (commutator_matrix(nmo, nmo))
         DO idir = 1, ncart
            ! Project, collect, and release each Cartesian MO operator before the next one.
            CALL ecd_project_operator_to_mo(linear_momentum_ao(idir), mo_coeff, linear_momentum_mo)
            ! cp_fm_get_submatrix replicates the distributed MO matrix on every rank.
            CALL cp_fm_get_submatrix(linear_momentum_mo, momentum_matrix, 1, 1, nmo, nmo)
            IF (include_nonlocal_velocity) THEN
               CALL ecd_project_operator_to_mo(commutator_ao(idir), mo_coeff, commutator_mo)
               CALL cp_fm_get_submatrix(commutator_mo, commutator_matrix, 1, 1, nmo, nmo)
               momentum_matrix = momentum_matrix + commutator_matrix
               CALL cp_fm_release(commutator_mo)
            END IF
            ! The real antisymmetric AO convention becomes a Hermitian velocity through -i.
            velocity_mo_matrix(idir, :, :) = CMPLX(0.0_dp, -momentum_matrix, KIND=dp)
            CALL cp_fm_release(linear_momentum_mo)
         END DO
         DEALLOCATE (momentum_matrix)
         IF (include_nonlocal_velocity) DEALLOCATE (commutator_matrix)
         ! Construct physical E1 and M1 matrices before selecting occupied-empty elements.
         ! The diagonal E1 element is the omitted Bloch connection and remains zero.
         e1_mo = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO ivirt = 1, nmo
            DO iocc = 1, nmo
               IF (ivirt == iocc) CYCLE
               IF (ABS(eigenvalues(ivirt) - eigenvalues(iocc)) <= energy_degeneracy_threshold) CYCLE
               e1_mo(:, ivirt, iocc) = CMPLX(0.0_dp, 1.0_dp, KIND=dp)* &
                                       velocity_mo_matrix(:, ivirt, iocc)/(eigenvalues(ivirt) - eigenvalues(iocc))
            END DO
         END DO
         ! M1 = -r^V x v/2 = E1 x v/2 because E1 = -r^V.
         m1_mo(1, :, :) = 0.5_dp*(MATMUL(e1_mo(2, :, :), velocity_mo_matrix(3, :, :)) - &
                                  MATMUL(e1_mo(3, :, :), velocity_mo_matrix(2, :, :)))
         m1_mo(2, :, :) = 0.5_dp*(MATMUL(e1_mo(3, :, :), velocity_mo_matrix(1, :, :)) - &
                                  MATMUL(e1_mo(1, :, :), velocity_mo_matrix(3, :, :)))
         m1_mo(3, :, :) = 0.5_dp*(MATMUL(e1_mo(1, :, :), velocity_mo_matrix(2, :, :)) - &
                                  MATMUL(e1_mo(2, :, :), velocity_mo_matrix(1, :, :)))
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

               e1_ai = e1_mo(:, ivirt, iocc)
               m1_ai = m1_mo(:, ivirt, iocc)
               ! Direct form of the Gamma E1--M1 tensor coefficient:
               ! Delta_alpha(mu,nu)_ai = 4*pi^2*w_ia/(gap*Omega) *
               !   [2 delta_mu,nu Im(mu^* dot m) - Im(mu_mu^*m_nu + mu_nu^*m_mu)].
               DO idir = 1, ncart
                  DO jdir = 1, ncart
                     spectral_tensor(idir, jdir, itransition) = 4.0_dp*pi*pi*occupation_difference/ &
                        (gap*cell_volume)*(2.0_dp*MERGE(AIMAG(SUM(CONJG(e1_ai)*m1_ai)), 0.0_dp, idir == jdir) - &
                        AIMAG(CONJG(e1_ai(idir))*m1_ai(jdir) + CONJG(e1_ai(jdir))*m1_ai(idir)))
                  END DO
               END DO
               IF (write_transition_matrices) THEN
                  DO idir = 1, ncart
                     WRITE (e1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
                        itransition, gap*evolt, ispin, iocc, ivirt, idir, REAL(e1_mo(idir, ivirt, iocc), KIND=dp), &
                        AIMAG(e1_mo(idir, ivirt, iocc)), occupation_difference
                     WRITE (m1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
                        itransition, gap*evolt, ispin, iocc, ivirt, idir, REAL(m1_mo(idir, ivirt, iocc), KIND=dp), &
                        AIMAG(m1_mo(idir, ivirt, iocc)), occupation_difference
                  END DO
               END IF
            END DO
         END DO
         DEALLOCATE (e1_mo, m1_mo, velocity_mo_matrix)
      END DO
      IF (itransition /= ntransition) CPABORT("ECD| Periodic velocity E1/M1 transition count is inconsistent")
      DO idir = 1, ncart
         CALL dbcsr_release_p(linear_momentum_ao(idir)%matrix)
         IF (include_nonlocal_velocity) CALL dbcsr_release_p(commutator_ao(idir)%matrix)
      END DO
      DEALLOCATE (linear_momentum_ao)
      IF (include_nonlocal_velocity) DEALLOCATE (commutator_ao)
      IF (write_transition_matrices) THEN
         CALL close_file(e1_unit)
         CALL close_file(m1_unit)
      END IF

      ! All E1, M1, and tensor denominators above use the original transition index.
      ! Sort only after assembly, so transition_order maps the energy-ordered output to it.
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

      DEALLOCATE (occupied, spectral_tensor, spin_channel, &
                  transition_energy, transition_order, virtual)

      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T2,A)') "ECD| Periodic Gamma E1 velocity matrix written to "//TRIM(e1_filename)
         WRITE (iounit, '(T2,A)') "ECD| Periodic Gamma M1 velocity matrix written to "//TRIM(m1_filename)
         WRITE (iounit, '(T2,A)') "ECD| Periodic Gamma E1--M1 tensor written to "//TRIM(tensor_filename)
         WRITE (iounit, '(T2,A,I0)') "ECD| Transition count ", ntransition
      END IF

   END SUBROUTINE ecd_spectra_xyz_ipa

END MODULE qs_ecd_spectra_xyz
