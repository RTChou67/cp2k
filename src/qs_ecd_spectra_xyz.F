!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program for atomistic simulations                                              !
!   Copyright 2000-2026 CP2K developers group                                                     !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gamma-point PERIODIC XYZ ECD spectra.
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
   REAL(KIND=dp), DIMENSION(ncart, ncart), PARAMETER, PRIVATE :: cartesian_identity = &
      RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [ncart, ncart])

   PUBLIC :: ecd_spectra_xyz

CONTAINS

! **************************************************************************************************
!> \brief Build Gamma E1 and M1 from the complete MO velocity matrix and write transitions.
!>
!> E1 uses the hypervirial velocity relation. M1 uses the Eq. (14) intermediate-state
!> sum, including its band-space Hermitian-conjugate term. The E1--M1-only coefficient
!> of the Gamma CD tensor is written for every transition.
!> \param pseudopotential_commutator Include the GTH non-local velocity commutator
! **************************************************************************************************
   SUBROUTINE ecd_spectra_xyz(qs_env, pseudopotential_commutator)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                 :: pseudopotential_commutator

      INTEGER                                            :: e1_unit, idir, iatom, iocc, inode, iounit, ispin, &
                                                            itransition, ivirt, jatom, jdir, last_jatom, m1_unit, &
                                                            ncommutator_block, nmo, nspin, ntransition, output_index, &
                                                            tensor_unit
      CHARACTER(LEN=default_path_length)                 :: e1_filename, m1_filename, tensor_filename
      LOGICAL                                            :: include_nonlocal_velocity
      REAL(KIND=dp)                                      :: cell_volume, gap, occupation_difference, velocity_angular_momentum_dot
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: angular_momentum_transition, electric_dipole_transition
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: velocity_matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: momentum_matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: commutator_matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: spectral_tensor
      COMPLEX(KIND=dp), DIMENSION(ncart, ncart)           :: x_matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: intermediate_weight_i, intermediate_weight_a, &
                                                            occupation_weight, transition_energy
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

      ALLOCATE (occupation_weight(ntransition), transition_energy(ntransition), &
                occupied(ntransition), spin_channel(ntransition), transition_order(ntransition), virtual(ntransition))
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
               occupation_weight(itransition) = occupation_difference
               occupied(itransition) = iocc
               spin_channel(itransition) = ispin
               virtual(itransition) = ivirt
            END DO
         END DO
      END DO
      IF (itransition /= ntransition) CPABORT("ECD| Periodic velocity E1/M1 transition count is inconsistent")

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
      ! Every rank therefore assembles the same complete velocity matrix.
      ALLOCATE (angular_momentum_transition(ncart, ntransition), electric_dipole_transition(ncart, ntransition), &
                spectral_tensor(ncart, ncart, ntransition))
      spectral_tensor = 0.0_dp
      DO ispin = 1, nspin
         NULLIFY (eigenvalues, mo_coeff)
         CALL get_mo_set(mos(ispin), eigenvalues=eigenvalues, mo_coeff=mo_coeff, nmo=nmo)
         IF (.NOT. ASSOCIATED(eigenvalues) .OR. .NOT. ASSOCIATED(mo_coeff)) THEN
            CPABORT("ECD| Periodic velocity E1/M1 output received incomplete Gamma MO data")
         END IF
         ALLOCATE (intermediate_weight_a(nmo), intermediate_weight_i(nmo), momentum_matrix(nmo, nmo), &
                   velocity_matrix(ncart, nmo, nmo))
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
            velocity_matrix(idir, :, :) = CMPLX(0.0_dp, -momentum_matrix, KIND=dp)
            CALL cp_fm_release(linear_momentum_mo)
         END DO
         DEALLOCATE (momentum_matrix)
         IF (include_nonlocal_velocity) DEALLOCATE (commutator_matrix)
         DO itransition = 1, ntransition
            IF (spin_channel(itransition) /= ispin) CYCLE
            ! E1 in the periodic velocity representation:
            ! mu_lambda,ai = i*pi_lambda,ai/(epsilon_a-epsilon_i).
            electric_dipole_transition(:, itransition) = CMPLX(0.0_dp, 1.0_dp, KIND=dp)* &
                                                         velocity_matrix(:, virtual(itransition), occupied(itransition))/ &
                                                         transition_energy(itransition)
            ! Eq. (14) velocity reconstruction of X_lambda,eta,ai.  In atomic units,
            ! A_lambda,eta,ai = (1/(2i)) sum_b v_eta,ab*v_lambda,bi/
            ! (epsilon_b-epsilon_i), and X = A_ai + A_ia^*.  The Cartesian order is
            ! transposed so X_lambda,eta yields r_lambda*p_eta in the local limit.
            intermediate_weight_i = 0.0_dp
            WHERE (ABS(eigenvalues - eigenvalues(occupied(itransition))) > energy_degeneracy_threshold)
               intermediate_weight_i = 1.0_dp/(eigenvalues - eigenvalues(occupied(itransition)))
            END WHERE
            intermediate_weight_a = 0.0_dp
            WHERE (ABS(eigenvalues - eigenvalues(virtual(itransition))) > energy_degeneracy_threshold)
               intermediate_weight_a = 1.0_dp/(eigenvalues - eigenvalues(virtual(itransition)))
            END WHERE
            x_matrix = CMPLX(0.0_dp, -0.5_dp, KIND=dp)* &
                       MATMUL(velocity_matrix(:, :, occupied(itransition))* &
                              SPREAD(intermediate_weight_i, DIM=1, NCOPIES=ncart), &
                              TRANSPOSE(velocity_matrix(:, virtual(itransition), :))) + &
                       CONJG(CMPLX(0.0_dp, -0.5_dp, KIND=dp)* &
                             MATMUL(velocity_matrix(:, :, virtual(itransition))* &
                                    SPREAD(intermediate_weight_a, DIM=1, NCOPIES=ncart), &
                                    TRANSPOSE(velocity_matrix(:, occupied(itransition), :))))
            ! L_ai = (X_yz-X_zy, X_zx-X_xz, X_xy-X_yx)_ai and, without SOC,
            ! b_ai = -L_ai/2.  No periodic position operator or reference point enters.
            angular_momentum_transition(1, itransition) = x_matrix(2, 3) - x_matrix(3, 2)
            angular_momentum_transition(2, itransition) = x_matrix(3, 1) - x_matrix(1, 3)
            angular_momentum_transition(3, itransition) = x_matrix(1, 2) - x_matrix(2, 1)
            ! Assemble the E1--M1 coefficient C^M1_ai of the Gamma CD tensor.
            ! With b=-L/2, off-diagonal elements are Re[pi_mu^* b_nu+pi_nu^* b_mu],
            ! while a diagonal element is -2 Re[sum_{lambda/=mu} pi_lambda^* b_lambda].
            spectral_tensor(:, :, itransition) = -0.5_dp*REAL( &
               MATMUL(RESHAPE(CONJG(velocity_matrix(:, virtual(itransition), occupied(itransition))), [ncart, 1]), &
                      RESHAPE(angular_momentum_transition(:, itransition), [1, ncart])), KIND=dp)
            spectral_tensor(:, :, itransition) = spectral_tensor(:, :, itransition) + &
                                                 TRANSPOSE(spectral_tensor(:, :, itransition))
            velocity_angular_momentum_dot = REAL( &
               SUM(CONJG(velocity_matrix(:, virtual(itransition), occupied(itransition)))* &
                   angular_momentum_transition(:, itransition)), KIND=dp)
            spectral_tensor(:, :, itransition) = spectral_tensor(:, :, itransition) + &
                                                 velocity_angular_momentum_dot*cartesian_identity
            ! Delta_alpha^Gamma_ai = 4*pi^2*w_ia*C^M1_ai /
            ! (Delta_epsilon_ia^2*Omega) delta(Delta_epsilon_ia-omega).
            spectral_tensor(:, :, itransition) = 4.0_dp*pi*pi*occupation_weight(itransition)* &
                                                 spectral_tensor(:, :, itransition)/ &
                                                 (transition_energy(itransition)**2*cell_volume)
         END DO
         DEALLOCATE (intermediate_weight_a, intermediate_weight_i, velocity_matrix)
      END DO
      DO idir = 1, ncart
         CALL dbcsr_release_p(linear_momentum_ao(idir)%matrix)
         IF (include_nonlocal_velocity) CALL dbcsr_release_p(commutator_ao(idir)%matrix)
      END DO
      DEALLOCATE (linear_momentum_ao)
      IF (include_nonlocal_velocity) DEALLOCATE (commutator_ao)

      ! All E1, M1, and tensor denominators above use the original transition index.
      ! Sort only after assembly, so transition_order maps the energy-ordered output to it.
      CALL sort(transition_energy, ntransition, transition_order)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      CALL cp_logger_generate_filename(logger, e1_filename, logger%iter_info%project_name, &
                                       "-ecd-e1-transition-matrix.dat", local=.FALSE.)
      CALL cp_logger_generate_filename(logger, m1_filename, logger%iter_info%project_name, &
                                       "-ecd-m1-transition-matrix.dat", local=.FALSE.)
      CALL cp_logger_generate_filename(logger, tensor_filename, logger%iter_info%project_name, &
                                       "-periodic-ipa-e1m1-tensor.dat", local=.FALSE.)
      IF (logger%para_env%is_source()) THEN
         CALL open_file(e1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=e1_unit)
         WRITE (e1_unit, '(A)') "# Periodic Gamma E1 transition matrix in the velocity representation"
         WRITE (e1_unit, '(A)') "# mu(a,i) = i*v(a,i)/(e(a)-e(i))"
         WRITE (e1_unit, '(A)') "# state energy_eV spin occupied virtual component dipole_real dipole_imag "// &
                                "occupation_difference"

         CALL open_file(m1_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=m1_unit)
         WRITE (m1_unit, '(A)') "# Periodic Gamma M1 matrix from the Eq. (14) velocity sum"
         WRITE (m1_unit, '(A)') "# X(lambda,eta) uses the Cartesian-transposed Eq. (14) product"
         WRITE (m1_unit, '(A)') "# X(lambda,eta; a,i) = A(lambda,eta; a,i) + A(lambda,eta; i,a)^*"
         WRITE (m1_unit, '(A)') "# L = (X(y,z)-X(z,y), X(z,x)-X(x,z), X(x,y)-X(y,x))"
         WRITE (m1_unit, '(A)') "# magnetic_dipole = -L/2; no periodic position reference is used"
         WRITE (m1_unit, '(A)') "# state energy_eV spin occupied virtual component angular_momentum_real "// &
                                "angular_momentum_imag "// &
                                "magnetic_dipole_real magnetic_dipole_imag occupation_difference"

         CALL open_file(tensor_filename, file_status="UNKNOWN", file_action="WRITE", unit_number=tensor_unit)
         WRITE (tensor_unit, '(A)') "# Periodic Gamma IPA CD tensor coefficient"
         WRITE (tensor_unit, '(A)') "# E1_M1_ONLY: D = D_M1 and D_E2 = 0"
         WRITE (tensor_unit, '(A)') "# Delta_alpha(mu,nu;omega) = sum_t coefficient(t,mu,nu)*delta(gap_t-omega)"
         WRITE (tensor_unit, '(A)') "# state energy_eV spin occupied virtual alpha beta eq17_coefficient_au"

         DO output_index = 1, ntransition
            itransition = transition_order(output_index)
            DO idir = 1, ncart
               WRITE (e1_unit, '(I0,1X,ES24.15E3,4(1X,I0),3(1X,ES24.15E3))') &
                  output_index, transition_energy(output_index)*evolt, spin_channel(itransition), &
                  occupied(itransition), virtual(itransition), idir, &
                  REAL(electric_dipole_transition(idir, itransition), KIND=dp), &
                  AIMAG(electric_dipole_transition(idir, itransition)), occupation_weight(itransition)
               WRITE (m1_unit, '(I0,1X,ES24.15E3,4(1X,I0),5(1X,ES24.15E3))') &
                  output_index, transition_energy(output_index)*evolt, spin_channel(itransition), &
                  occupied(itransition), virtual(itransition), idir, &
                  REAL(angular_momentum_transition(idir, itransition), KIND=dp), &
                  AIMAG(angular_momentum_transition(idir, itransition)), &
                  REAL(-0.5_dp*angular_momentum_transition(idir, itransition), KIND=dp), &
                  AIMAG(-0.5_dp*angular_momentum_transition(idir, itransition)), occupation_weight(itransition)
               DO jdir = 1, ncart
                  WRITE (tensor_unit, '(I0,1X,ES24.15E3,5(1X,I0),1X,ES24.15E3)') &
                     output_index, transition_energy(output_index)*evolt, spin_channel(itransition), &
                     occupied(itransition), virtual(itransition), idir, jdir, &
                     spectral_tensor(idir, jdir, itransition)
               END DO
            END DO
         END DO
         CALL close_file(e1_unit)
         CALL close_file(m1_unit)
         CALL close_file(tensor_unit)
      END IF

      DEALLOCATE (angular_momentum_transition, electric_dipole_transition, occupation_weight, &
                  occupied, spectral_tensor, spin_channel, &
                  transition_energy, transition_order, virtual)

      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T2,A)') "ECD| Periodic Gamma E1 velocity matrix written to "//TRIM(e1_filename)
         WRITE (iounit, '(T2,A)') "ECD| Periodic Gamma M1 velocity matrix written to "//TRIM(m1_filename)
         WRITE (iounit, '(T2,A)') "ECD| Periodic Gamma E1--M1 tensor written to "//TRIM(tensor_filename)
         WRITE (iounit, '(T2,A,I0)') "ECD| Transition count ", ntransition
      END IF

   END SUBROUTINE ecd_spectra_xyz

END MODULE qs_ecd_spectra_xyz
