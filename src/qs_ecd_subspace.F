!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Degenerate-band manifold utilities for periodic ECD.
! **************************************************************************************************
MODULE qs_ecd_subspace
   USE cp_cfm_diag,                    ONLY: cp_cfm_heevd
   USE cp_cfm_types,                   ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_subspace'

   TYPE, PUBLIC :: ecd_band_manifold_type
      INTEGER :: first_band = 0
      INTEGER :: last_band = -1
   END TYPE ecd_band_manifold_type

   PUBLIC :: ecd_band_manifolds_build, &
             ecd_band_manifolds_release, &
             ecd_band_manifolds_validate_window, &
             ecd_band_first_order_response, &
             ecd_band_taylor_non_degenerate_response, &
             ecd_taylor_retarded_response

CONTAINS

! **************************************************************************************************
!> \brief Group consecutive eigenvalues into complete EPS_DEGENERACY manifolds.
! **************************************************************************************************
   SUBROUTINE ecd_band_manifolds_build(manifolds, eigenvalues, eps_degeneracy)
      TYPE(ecd_band_manifold_type), DIMENSION(:), &
         POINTER                                         :: manifolds
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      REAL(KIND=dp), INTENT(IN)                          :: eps_degeneracy

      INTEGER                                            :: iband, imanifold, nmanifolds, nmo

      IF (eps_degeneracy <= 0.0_dp) THEN
         CPABORT("ECD| EPS_DEGENERACY must be positive when forming band manifolds")
      END IF
      nmo = SIZE(eigenvalues)
      IF (nmo < 1) THEN
         CPABORT("ECD| Cannot form band manifolds from an empty eigenvalue list")
      END IF

      nmanifolds = 1
      DO iband = 2, nmo
         IF (ABS(eigenvalues(iband) - eigenvalues(iband - 1)) > eps_degeneracy) nmanifolds = nmanifolds + 1
      END DO
      ALLOCATE (manifolds(nmanifolds))

      imanifold = 1
      manifolds(imanifold)%first_band = 1
      DO iband = 2, nmo
         IF (ABS(eigenvalues(iband) - eigenvalues(iband - 1)) > eps_degeneracy) THEN
            manifolds(imanifold)%last_band = iband - 1
            imanifold = imanifold + 1
            manifolds(imanifold)%first_band = iband
         END IF
      END DO
      manifolds(imanifold)%last_band = nmo

   END SUBROUTINE ecd_band_manifolds_build

! **************************************************************************************************
!> \brief Release a manifold list.
! **************************************************************************************************
   SUBROUTINE ecd_band_manifolds_release(manifolds)
      TYPE(ecd_band_manifold_type), DIMENSION(:), &
         POINTER                                         :: manifolds

      IF (ASSOCIATED(manifolds)) DEALLOCATE (manifolds)
      NULLIFY (manifolds)

   END SUBROUTINE ecd_band_manifolds_release

! **************************************************************************************************
!> \brief Require a selected band interval to contain every intersected manifold.
! **************************************************************************************************
   SUBROUTINE ecd_band_manifolds_validate_window(manifolds, band_range)
      TYPE(ecd_band_manifold_type), DIMENSION(:), &
         POINTER                                         :: manifolds
      INTEGER, DIMENSION(2), INTENT(IN)                  :: band_range

      INTEGER                                            :: imanifold

      IF (.NOT. ASSOCIATED(manifolds)) THEN
         CPABORT("ECD| Band-manifold metadata are unavailable")
      END IF
      IF (band_range(1) < 1 .OR. band_range(1) > band_range(2)) THEN
         CPABORT("ECD| Invalid band range for manifold validation")
      END IF
      IF (band_range(2) > manifolds(SIZE(manifolds))%last_band) THEN
         CPABORT("ECD| Transition band range exceeds the available explicit-Gamma MOs")
      END IF
      DO imanifold = 1, SIZE(manifolds)
         IF (manifolds(imanifold)%last_band < band_range(1)) CYCLE
         IF (manifolds(imanifold)%first_band > band_range(2)) CYCLE
         IF (manifolds(imanifold)%first_band < band_range(1) .OR. &
             manifolds(imanifold)%last_band > band_range(2)) THEN
            CPABORT("ECD| A selected transition band range cuts through a degenerate manifold")
         END IF
      END DO

   END SUBROUTINE ecd_band_manifolds_validate_window

! **************************************************************************************************
!> \brief Build the first directional generalized-eigenproblem response in the MO basis.
!> \param eigenvalue_reference Degenerate-manifold reference energies
!> \param eigenvalue_derivative First directional band-energy coefficients
!> \param equation_residual_norm First-order projected eigenproblem residual norm
!> \param normalization_residual_norm First-order MO normalization residual norm
!> \param manifold_residual_norm Projected manifold diagonalization residual norm
!> \param eigenvector_derivative_norm Parallel-transport eigenvector derivative norm
!> \param hamiltonian_derivative First directional Hamiltonian coefficient in the MO basis
!> \param overlap_derivative First directional overlap coefficient in the MO basis
!> \param eigenvalues Zero-order band energies
!> \param manifolds Complete degenerate band manifolds
!> \param band_rotation Optional manifold-resolved Gamma band rotation
!> \param coefficient_derivative_output Optional parallel-transport band response
! **************************************************************************************************
   SUBROUTINE ecd_band_first_order_response(eigenvalue_reference, eigenvalue_derivative, &
                                             equation_residual_norm, normalization_residual_norm, &
                                             manifold_residual_norm, eigenvector_derivative_norm, &
                                             hamiltonian_derivative, overlap_derivative, eigenvalues, manifolds, &
                                             band_rotation, coefficient_derivative_output)
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: eigenvalue_reference, &
                                                               eigenvalue_derivative
      REAL(KIND=dp), INTENT(OUT)                          :: equation_residual_norm, &
                                                            normalization_residual_norm, &
                                                            manifold_residual_norm, &
                                                            eigenvector_derivative_norm
      TYPE(cp_cfm_type), INTENT(IN)                       :: hamiltonian_derivative, overlap_derivative
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)             :: eigenvalues
      TYPE(ecd_band_manifold_type), DIMENSION(:), &
         INTENT(IN)                                       :: manifolds
      COMPLEX(KIND=dp), DIMENSION(:, :), ALLOCATABLE, &
         OPTIONAL, INTENT(OUT)                            :: band_rotation, coefficient_derivative_output

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)      :: coefficient_derivative, equation_residual, &
                                                            hamiltonian_mo, hamiltonian_rotated, &
                                                            manifold_vectors, normalization_residual, &
                                                            overlap_mo, overlap_rotated, rotation, &
                                                            work_block
      INTEGER                                             :: first_band, iband, imanifold, jband, last_band, &
                                                             nband, nmanifold
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: manifold_index
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: manifold_eigenvalues
      TYPE(cp_cfm_type)                                   :: manifold_matrix, manifold_rotation

      nband = SIZE(eigenvalues)
      IF (nband < 1) THEN
         CPABORT("ECD| First-order band response requires at least one band")
      END IF
      IF (.NOT. ASSOCIATED(hamiltonian_derivative%local_data) .OR. &
          .NOT. ASSOCIATED(overlap_derivative%local_data)) THEN
         CPABORT("ECD| First-order band response requires projected H and S coefficients")
      END IF
      IF (hamiltonian_derivative%matrix_struct%nrow_global /= nband .OR. &
          hamiltonian_derivative%matrix_struct%ncol_global /= nband .OR. &
          overlap_derivative%matrix_struct%nrow_global /= nband .OR. &
          overlap_derivative%matrix_struct%ncol_global /= nband) THEN
         CPABORT("ECD| First-order band response received inconsistent MO matrix dimensions")
      END IF
      IF (SIZE(manifolds) < 1) THEN
         CPABORT("ECD| First-order band response requires degenerate-manifold metadata")
      END IF

      ALLOCATE (eigenvalue_reference(nband), eigenvalue_derivative(nband), manifold_index(nband), &
                hamiltonian_mo(nband, nband), hamiltonian_rotated(nband, nband), &
                overlap_mo(nband, nband), overlap_rotated(nband, nband), &
                rotation(nband, nband), coefficient_derivative(nband, nband), &
                equation_residual(nband, nband), normalization_residual(nband, nband))
      CALL cp_cfm_get_submatrix(hamiltonian_derivative, hamiltonian_mo)
      CALL cp_cfm_get_submatrix(overlap_derivative, overlap_mo)

      rotation(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      manifold_index(:) = 0
      eigenvalue_reference(:) = 0.0_dp
      eigenvalue_derivative(:) = 0.0_dp
      DO iband = 1, nband
         rotation(iband, iband) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      END DO

      DO imanifold = 1, SIZE(manifolds)
         first_band = manifolds(imanifold)%first_band
         last_band = manifolds(imanifold)%last_band
         IF (first_band < 1 .OR. last_band < first_band .OR. last_band > nband) THEN
            CPABORT("ECD| First-order band response encountered invalid manifold bounds")
         END IF
         nmanifold = last_band - first_band + 1
         eigenvalue_reference(first_band:last_band) = &
            SUM(eigenvalues(first_band:last_band))/REAL(nmanifold, KIND=dp)
         manifold_index(first_band:last_band) = imanifold

         ALLOCATE (work_block(nmanifold, nmanifold), manifold_vectors(nmanifold, nmanifold), &
                   manifold_eigenvalues(nmanifold))
         work_block(:, :) = hamiltonian_mo(first_band:last_band, first_band:last_band) - &
                              eigenvalue_reference(first_band)* &
                              overlap_mo(first_band:last_band, first_band:last_band)
         CALL cp_cfm_create(manifold_matrix, hamiltonian_derivative%matrix_struct, &
                            nrow=nmanifold, ncol=nmanifold, set_zero=.TRUE.)
         CALL cp_cfm_create(manifold_rotation, hamiltonian_derivative%matrix_struct, &
                            nrow=nmanifold, ncol=nmanifold, set_zero=.TRUE.)
         CALL cp_cfm_set_submatrix(manifold_matrix, work_block)
         CALL cp_cfm_heevd(manifold_matrix, manifold_rotation, manifold_eigenvalues)
         CALL cp_cfm_get_submatrix(manifold_rotation, manifold_vectors)
         rotation(first_band:last_band, first_band:last_band) = manifold_vectors
         eigenvalue_derivative(first_band:last_band) = manifold_eigenvalues
         CALL cp_cfm_release(manifold_matrix)
         CALL cp_cfm_release(manifold_rotation)
         DEALLOCATE (manifold_eigenvalues, manifold_vectors, work_block)
      END DO

      hamiltonian_rotated = MATMUL(CONJG(TRANSPOSE(rotation)), MATMUL(hamiltonian_mo, rotation))
      overlap_rotated = MATMUL(CONJG(TRANSPOSE(rotation)), MATMUL(overlap_mo, rotation))
      coefficient_derivative(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO imanifold = 1, SIZE(manifolds)
         first_band = manifolds(imanifold)%first_band
         last_band = manifolds(imanifold)%last_band
         coefficient_derivative(first_band:last_band, first_band:last_band) = -0.5_dp* &
            overlap_rotated(first_band:last_band, first_band:last_band)
      END DO
      DO iband = 1, nband
         DO jband = 1, nband
            IF (manifold_index(iband) == manifold_index(jband)) CYCLE
            coefficient_derivative(iband, jband) = &
               (eigenvalue_reference(jband)*overlap_rotated(iband, jband) - &
                hamiltonian_rotated(iband, jband))/ &
               (eigenvalue_reference(iband) - eigenvalue_reference(jband))
         END DO
      END DO

      equation_residual(:, :) = hamiltonian_rotated(:, :)
      DO iband = 1, nband
         DO jband = 1, nband
            equation_residual(iband, jband) = equation_residual(iband, jband) + &
               (eigenvalue_reference(iband) - eigenvalue_reference(jband))* &
               coefficient_derivative(iband, jband) - &
               eigenvalue_reference(jband)*overlap_rotated(iband, jband)
         END DO
         equation_residual(iband, iband) = equation_residual(iband, iband) - &
                                             eigenvalue_derivative(iband)
      END DO
      normalization_residual = coefficient_derivative + CONJG(TRANSPOSE(coefficient_derivative)) + &
                               overlap_rotated

      manifold_residual_norm = 0.0_dp
      DO imanifold = 1, SIZE(manifolds)
         first_band = manifolds(imanifold)%first_band
         last_band = manifolds(imanifold)%last_band
         nmanifold = last_band - first_band + 1
         ALLOCATE (work_block(nmanifold, nmanifold))
         work_block = hamiltonian_rotated(first_band:last_band, first_band:last_band) - &
                      eigenvalue_reference(first_band)*overlap_rotated(first_band:last_band, first_band:last_band)
         DO iband = 1, nmanifold
            work_block(iband, iband) = work_block(iband, iband) - &
                                        eigenvalue_derivative(first_band + iband - 1)
         END DO
         manifold_residual_norm = manifold_residual_norm + SUM(ABS(work_block)**2)
         DEALLOCATE (work_block)
      END DO
      manifold_residual_norm = SQRT(manifold_residual_norm)
      equation_residual_norm = SQRT(SUM(ABS(equation_residual)**2))
      normalization_residual_norm = SQRT(SUM(ABS(normalization_residual)**2))
      eigenvector_derivative_norm = SQRT(SUM(ABS(coefficient_derivative)**2))

      IF (PRESENT(band_rotation)) THEN
         ALLOCATE (band_rotation(nband, nband))
         band_rotation = rotation
      END IF
      IF (PRESENT(coefficient_derivative_output)) THEN
         ALLOCATE (coefficient_derivative_output(nband, nband))
         coefficient_derivative_output = coefficient_derivative
      END IF

      DEALLOCATE (coefficient_derivative, equation_residual, hamiltonian_mo, hamiltonian_rotated, &
                  manifold_index, normalization_residual, overlap_mo, overlap_rotated, rotation)

   END SUBROUTINE ecd_band_first_order_response

! **************************************************************************************************
!> \brief Solve normalized Taylor coefficients of a non-degenerate generalized eigenproblem.
!> \param hamiltonian_coeff Normalized MO Hamiltonian coefficients, including order zero
!> \param overlap_coeff Normalized MO overlap coefficients, including order zero
!> \param eigenvalues Zero-order band energies
!> \param manifolds Complete zero-order band-manifold metadata
!> \param energy_coeff Band-energy Taylor coefficients
!> \param coefficient_coeff MO eigenvector Taylor coefficients in the order-zero MO frame
!> \param eigenproblem_residual_norm Coefficient-wise generalized-eigenproblem residual norms
!> \param normalization_residual_norm Coefficient-wise MO-normalization residual norms
! **************************************************************************************************
   SUBROUTINE ecd_band_taylor_non_degenerate_response(hamiltonian_coeff, overlap_coeff, eigenvalues, manifolds, &
                                                       energy_coeff, coefficient_coeff, &
                                                       eigenproblem_residual_norm, normalization_residual_norm)
      COMPLEX(KIND=dp), DIMENSION(:, :, 0:), INTENT(IN)  :: hamiltonian_coeff, overlap_coeff
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      TYPE(ecd_band_manifold_type), DIMENSION(:), &
         INTENT(IN)                                      :: manifolds
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: energy_coeff
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: coefficient_coeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: eigenproblem_residual_norm, normalization_residual_norm

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: eigenproblem_residual, normalization_known, &
                                                            normalization_residual, rhs
      INTEGER                                            :: ia, ib, ic, iorder, imanifold, jorder, &
                                                            nband, max_order
      REAL(KIND=dp)                                      :: energy_denominator

      nband = SIZE(eigenvalues)
      max_order = UBOUND(hamiltonian_coeff, 3)
      IF (nband < 1) THEN
         CPABORT("ECD| Taylor band response requires at least one band")
      END IF
      IF (LBOUND(hamiltonian_coeff, 3) /= 0 .OR. LBOUND(overlap_coeff, 3) /= 0) THEN
         CPABORT("ECD| Taylor band-response coefficients must start at order zero")
      END IF
      IF (SIZE(hamiltonian_coeff, 1) /= nband .OR. SIZE(hamiltonian_coeff, 2) /= nband .OR. &
          SIZE(overlap_coeff, 1) /= nband .OR. SIZE(overlap_coeff, 2) /= nband .OR. &
          UBOUND(overlap_coeff, 3) /= max_order) THEN
         CPABORT("ECD| Taylor band response received inconsistent MO coefficient dimensions")
      END IF
      IF (max_order < 1) THEN
         CPABORT("ECD| Taylor band response requires at least first order")
      END IF
      DO imanifold = 1, SIZE(manifolds)
         IF (manifolds(imanifold)%first_band /= manifolds(imanifold)%last_band) THEN
            CPABORT("ECD| Higher-order response requires a covariant degenerate-manifold solver")
         END IF
      END DO

      ALLOCATE (coefficient_coeff(nband, nband, 0:max_order), energy_coeff(nband, 0:max_order), &
                eigenproblem_residual_norm(0:max_order), normalization_residual_norm(0:max_order), &
                eigenproblem_residual(nband, nband), normalization_known(nband, nband), &
                normalization_residual(nband, nband), rhs(nband, nband))
      coefficient_coeff = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      energy_coeff = 0.0_dp
      eigenproblem_residual_norm = 0.0_dp
      normalization_residual_norm = 0.0_dp
      DO ia = 1, nband
         coefficient_coeff(ia, ia, 0) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         energy_coeff(ia, 0) = eigenvalues(ia)
      END DO

      DO iorder = 1, max_order
         rhs = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO ia = 0, iorder
            DO ib = 0, iorder - ia
               ic = iorder - ia - ib
               IF (ia == 0 .AND. ib == iorder .AND. ic == 0) CYCLE
               IF (ia == 0 .AND. ib == 0 .AND. ic == iorder) CYCLE
               DO jorder = 1, nband
                  rhs(:, jorder) = rhs(:, jorder) + energy_coeff(jorder, ic)* &
                                    MATMUL(overlap_coeff(:, :, ia), coefficient_coeff(:, jorder, ib))
               END DO
            END DO
         END DO
         DO ia = 1, iorder
            rhs = rhs - MATMUL(hamiltonian_coeff(:, :, ia), coefficient_coeff(:, :, iorder - ia))
         END DO
         DO ia = 1, nband
            energy_coeff(ia, iorder) = -REAL(rhs(ia, ia), KIND=dp)
            rhs(ia, ia) = rhs(ia, ia) + CMPLX(energy_coeff(ia, iorder), 0.0_dp, KIND=dp)
         END DO
         DO ib = 1, nband
            DO ia = 1, nband
               IF (ia == ib) CYCLE
               energy_denominator = eigenvalues(ia) - eigenvalues(ib)
               IF (ABS(energy_denominator) <= SQRT(TINY(1.0_dp))) THEN
                  CPABORT("ECD| Higher-order response encountered an unresolved energy degeneracy")
               END IF
               coefficient_coeff(ia, ib, iorder) = rhs(ia, ib)/energy_denominator
            END DO
         END DO
         normalization_known = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO ia = 0, iorder
            DO ib = 0, iorder - ia
               ic = iorder - ia - ib
               IF (ia == iorder .AND. ib == 0 .AND. ic == 0) CYCLE
               IF (ia == 0 .AND. ib == 0 .AND. ic == iorder) CYCLE
               normalization_known = normalization_known + &
                  MATMUL(CONJG(TRANSPOSE(coefficient_coeff(:, :, ia))), &
                         MATMUL(overlap_coeff(:, :, ib), coefficient_coeff(:, :, ic)))
            END DO
         END DO
         DO ia = 1, nband
            coefficient_coeff(ia, ia, iorder) = -0.5_dp*normalization_known(ia, ia)
         END DO

         eigenproblem_residual = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         normalization_residual = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO ia = 0, iorder
            eigenproblem_residual = eigenproblem_residual + &
               MATMUL(hamiltonian_coeff(:, :, ia), coefficient_coeff(:, :, iorder - ia))
            DO ib = 0, iorder - ia
               ic = iorder - ia - ib
               DO jorder = 1, nband
                  eigenproblem_residual(:, jorder) = eigenproblem_residual(:, jorder) - &
                     energy_coeff(jorder, ic)*MATMUL(overlap_coeff(:, :, ia), &
                                                       coefficient_coeff(:, jorder, ib))
               END DO
            END DO
         END DO
         DO ia = 0, iorder
            DO ib = 0, iorder - ia
               ic = iorder - ia - ib
               normalization_residual = normalization_residual + &
                  MATMUL(CONJG(TRANSPOSE(coefficient_coeff(:, :, ia))), &
                         MATMUL(overlap_coeff(:, :, ib), coefficient_coeff(:, :, ic)))
            END DO
         END DO
         eigenproblem_residual_norm(iorder) = SQRT(SUM(ABS(eigenproblem_residual)**2))
         normalization_residual_norm(iorder) = SQRT(SUM(ABS(normalization_residual)**2))
      END DO

      DEALLOCATE (eigenproblem_residual, normalization_known, normalization_residual, rhs)

   END SUBROUTINE ecd_band_taylor_non_degenerate_response

! **************************************************************************************************
!> \brief Assemble arbitrary-order resonant and counter-rotating IPA response coefficients.
!> \param transition_coeff Current-vertex Taylor coefficients indexed by Cartesian component and order
!> \param gap_coeff Positive transition-gap Taylor coefficients
!> \param photon_energy Positive photon energy in hartree
!> \param half_width Lorentzian half width in hartree
!> \param resonant_coeff Coefficients of K(q)/[omega-gap(q)+i eta]
!> \param counter_coeff Coefficients of -K_ba(-q)/[omega+gap(-q)+i eta]
! **************************************************************************************************
   SUBROUTINE ecd_taylor_retarded_response(transition_coeff, gap_coeff, photon_energy, half_width, &
                                            resonant_coeff, counter_coeff)
      COMPLEX(KIND=dp), DIMENSION(:, 0:), INTENT(IN)     :: transition_coeff
      COMPLEX(KIND=dp), DIMENSION(:, :, 0:), INTENT(OUT) :: resonant_coeff, counter_coeff
      REAL(KIND=dp), DIMENSION(0:), INTENT(IN)           :: gap_coeff
      REAL(KIND=dp), INTENT(IN)                          :: photon_energy, half_width

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: counter_inverse, counter_denominator, &
                                                            resonant_inverse, resonant_denominator
      INTEGER                                            :: ia, ib, iorder, jorder, max_order, ncart

      ncart = SIZE(transition_coeff, 1)
      max_order = UBOUND(transition_coeff, 2)
      IF (ncart < 1) THEN
         CPABORT("ECD| Taylor response requires at least one Cartesian component")
      END IF
      IF (LBOUND(transition_coeff, 2) /= 0 .OR. LBOUND(gap_coeff, 1) /= 0) THEN
         CPABORT("ECD| Taylor response coefficients must start at order zero")
      END IF
      IF (UBOUND(gap_coeff, 1) /= max_order) THEN
         CPABORT("ECD| Taylor response received inconsistent transition and gap orders")
      END IF
      IF (SIZE(resonant_coeff, 1) /= ncart .OR. SIZE(resonant_coeff, 2) /= ncart .OR. &
          LBOUND(resonant_coeff, 3) /= 0 .OR. UBOUND(resonant_coeff, 3) /= max_order) THEN
         CPABORT("ECD| Taylor response received incompatible resonant-coefficient storage")
      END IF
      IF (SIZE(counter_coeff, 1) /= ncart .OR. SIZE(counter_coeff, 2) /= ncart .OR. &
          LBOUND(counter_coeff, 3) /= 0 .OR. UBOUND(counter_coeff, 3) /= max_order) THEN
         CPABORT("ECD| Taylor response received incompatible counter-coefficient storage")
      END IF
      IF (photon_energy <= 0.0_dp .OR. half_width <= 0.0_dp) THEN
         CPABORT("ECD| Taylor response requires positive photon energy and half width")
      END IF

      ALLOCATE (kernel(ncart*ncart, 0:max_order), counter_denominator(0:max_order), &
                counter_inverse(0:max_order), resonant_denominator(0:max_order), resonant_inverse(0:max_order))
      kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO iorder = 0, max_order
         DO jorder = 0, iorder
            DO ia = 1, ncart
               DO ib = 1, ncart
                  kernel((ia - 1)*ncart + ib, iorder) = kernel((ia - 1)*ncart + ib, iorder) + &
                     CONJG(transition_coeff(ia, jorder))*transition_coeff(ib, iorder - jorder)
               END DO
            END DO
         END DO
      END DO

      resonant_denominator(0) = CMPLX(photon_energy - gap_coeff(0), half_width, KIND=dp)
      counter_denominator(0) = CMPLX(photon_energy + gap_coeff(0), half_width, KIND=dp)
      DO iorder = 1, max_order
         resonant_denominator(iorder) = CMPLX(-gap_coeff(iorder), 0.0_dp, KIND=dp)
         counter_denominator(iorder) = CMPLX((-1)**iorder*gap_coeff(iorder), 0.0_dp, KIND=dp)
      END DO
      CALL ecd_taylor_scalar_inverse(resonant_denominator, resonant_inverse)
      CALL ecd_taylor_scalar_inverse(counter_denominator, counter_inverse)

      resonant_coeff(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      counter_coeff(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO iorder = 0, max_order
         DO jorder = 0, iorder
            DO ia = 1, ncart
               DO ib = 1, ncart
                  resonant_coeff(ia, ib, iorder) = resonant_coeff(ia, ib, iorder) + &
                     kernel((ia - 1)*ncart + ib, jorder)*resonant_inverse(iorder - jorder)
                  counter_coeff(ia, ib, iorder) = counter_coeff(ia, ib, iorder) - &
                     (-1)**jorder*kernel((ib - 1)*ncart + ia, jorder)*counter_inverse(iorder - jorder)
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (counter_denominator, counter_inverse, kernel, resonant_denominator, resonant_inverse)

   END SUBROUTINE ecd_taylor_retarded_response

! **************************************************************************************************
!> \brief Invert a normalized scalar Taylor series.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_inverse(series, inverse)
      COMPLEX(KIND=dp), DIMENSION(0:), INTENT(IN)        :: series
      COMPLEX(KIND=dp), DIMENSION(0:), INTENT(OUT)       :: inverse

      INTEGER                                            :: iorder, jorder, max_order

      max_order = UBOUND(series, 1)
      IF (UBOUND(inverse, 1) /= max_order) THEN
         CPABORT("ECD| Taylor inverse received inconsistent storage")
      END IF
      IF (ABS(series(0)) <= SQRT(TINY(1.0_dp))) THEN
         CPABORT("ECD| Taylor inverse requires a non-zero zeroth coefficient")
      END IF
      inverse(0) = 1.0_dp/series(0)
      DO iorder = 1, max_order
         inverse(iorder) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO jorder = 1, iorder
            inverse(iorder) = inverse(iorder) - series(jorder)*inverse(iorder - jorder)
         END DO
         inverse(iorder) = inverse(iorder)/series(0)
      END DO

   END SUBROUTINE ecd_taylor_scalar_inverse

END MODULE qs_ecd_subspace
