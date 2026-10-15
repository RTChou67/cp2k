!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Dense screened-interaction algebra for real and complex full matrices.
! **************************************************************************************************
MODULE mbpt_screening
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_add_on_diag, cp_cfm_uplo_to_full
   USE cp_cfm_types, ONLY: cp_cfm_create, cp_cfm_release, cp_cfm_type
   USE cp_fm_basic_linalg, ONLY: cp_fm_add_on_diag, cp_fm_uplo_to_full
   USE cp_fm_types, ONLY: cp_fm_create, cp_fm_release, cp_fm_type
   USE gw_utils_fm, ONLY: fm_invert, cfm_invert, fm_sqrt, cfm_sqrt
   USE kinds, ONLY: dp
   USE mathconstants, ONLY: z_one, z_zero
   USE parallel_gemm_api, ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mbpt_screening'
   PUBLIC :: screened_interaction, prepare_screening_factors

   INTERFACE screened_interaction
      MODULE PROCEDURE screened_interaction_real, screened_interaction_complex
   END INTERFACE screened_interaction

   INTERFACE prepare_screening_factors
      MODULE PROCEDURE prepare_screening_factors_real, prepare_screening_factors_complex
   END INTERFACE prepare_screening_factors

CONTAINS

   #:for suffix, prefix, transpose, one, zero in [('real', 'cp_fm', 'T', '1.0_dp', '0.0_dp'), ('complex', 'cp_cfm', 'C', 'z_one', 'z_zero')]
! **************************************************************************************************
!> \brief Prepare Coulomb factors V = B^H B and F = M^-1 B^H for dielectric screening.
!>        B may be triangular. The input V and metric inverse are retained.
!> \param matrix_V Coulomb matrix in the auxiliary basis
!> \param metric_inverse inverse auxiliary metric
!> \param coulomb_factor output factor B
!> \param metric_coulomb_factor output factor F
!> \param eigenvalue_threshold Coulomb cutoff on the spectral fallback
!> \param unit_nr optional diagnostic output unit
!> \param workspace optional reusable full-matrix workspace
! **************************************************************************************************
      SUBROUTINE prepare_screening_factors_${suffix}$ (matrix_V, metric_inverse, coulomb_factor, &
                                                       metric_coulomb_factor, eigenvalue_threshold, &
                                                       unit_nr, workspace)
         TYPE(${prefix}$_type), INTENT(IN) :: matrix_V, metric_inverse
         TYPE(${prefix}$_type), INTENT(INOUT) :: coulomb_factor, metric_coulomb_factor
         REAL(KIND=dp), INTENT(IN) :: eigenvalue_threshold
         INTEGER, INTENT(IN), OPTIONAL :: unit_nr
         TYPE(${prefix}$_type), INTENT(INOUT), OPTIONAL, TARGET :: workspace

         INTEGER :: n

         n = matrix_V%matrix_struct%nrow_global
         CALL ${prefix[3:]}$_sqrt(matrix_V, coulomb_factor, eigenvalue_threshold, unit_nr, workspace)
         CALL parallel_gemm('N', '${transpose}$', n, n, n, ${one}$, &
                            metric_inverse, coulomb_factor, ${zero}$, metric_coulomb_factor)
      END SUBROUTINE prepare_screening_factors_${suffix}$

! **************************************************************************************************
!> \brief Compute the correlation-only screened interaction from a positive response matrix.
!>        The input is the positive response stored by the imaginary-time paths, i.e. minus
!>        the physical polarizability. For a Coulomb factor satisfying V = B^H B and the
!>        metric-weighted factor F = M^-1 B^H, form epsilon = I + F^H response F and
!>        Wc = B^H (epsilon^-1 - I) B. H denotes T for real matrices and C for complex ones.
!>        The factor B may be triangular; no symmetric-square-root assumption is made.
!>        Cholesky inversion is used first, with a positive-spectrum inverse as fallback.
!> \param response input response; overwritten by Wc when matrix_Wc is absent
!> \param coulomb_factor factor B satisfying V = B^H B
!> \param metric_coulomb_factor F = M^-1 B^H
!> \param eigenvalue_threshold nonnegative eigenvalue cutoff for the inverse fallback
!> \param matrix_Wc optional separate output; retains response when supplied
!> \param unit_nr optional output unit for the inverse fallback diagnostics
! **************************************************************************************************
      SUBROUTINE screened_interaction_${suffix}$ (response, coulomb_factor, metric_coulomb_factor, &
                                                  eigenvalue_threshold, matrix_Wc, unit_nr)
         TYPE(${prefix}$_type), INTENT(INOUT), TARGET :: response
         TYPE(${prefix}$_type), INTENT(IN) :: coulomb_factor, metric_coulomb_factor
         REAL(KIND=dp), INTENT(IN) :: eigenvalue_threshold
         TYPE(${prefix}$_type), INTENT(INOUT), OPTIONAL, TARGET :: matrix_Wc
         INTEGER, INTENT(IN), OPTIONAL :: unit_nr

         CHARACTER(LEN=*), PARAMETER :: routineN = 'screened_interaction_${suffix}$'
         INTEGER :: handle, n
         TYPE(${prefix}$_type) :: work
         TYPE(${prefix}$_type), POINTER :: epsilon

         CALL timeset(routineN, handle)
         n = response%matrix_struct%nrow_global
         CPASSERT(n == response%matrix_struct%ncol_global)
         CPASSERT(eigenvalue_threshold >= 0.0_dp)
         IF (PRESENT(matrix_Wc)) THEN
            epsilon => matrix_Wc
         ELSE
            epsilon => response
         END IF
         CALL ${prefix}$_create(work, response%matrix_struct)

         ! epsilon = I + (M^-1 B^H)^H response (M^-1 B^H).
         CALL parallel_gemm('N', 'N', n, n, n, ${one}$, response, metric_coulomb_factor, ${zero}$, work)
         CALL parallel_gemm('${transpose}$', 'N', n, n, n, ${one}$, &
                            metric_coulomb_factor, work, ${zero}$, epsilon)
         CALL ${prefix}$_add_on_diag(epsilon, ${one}$)
         CALL ${prefix}$_uplo_to_full(epsilon, work)

         ! Retain epsilon for a failed Cholesky decomposition, reusing the same workspace.
         CALL ${prefix[3:]}$_invert(epsilon, eigenvalue_threshold, unit_nr, workspace=work)

         ! Wc = B^H (epsilon^-1 - I) B, including triangular Coulomb factors.
         CALL ${prefix}$_add_on_diag(epsilon, -${one}$)
         CALL parallel_gemm('N', 'N', n, n, n, ${one}$, epsilon, coulomb_factor, ${zero}$, work)
         CALL parallel_gemm('${transpose}$', 'N', n, n, n, ${one}$, coulomb_factor, work, ${zero}$, epsilon)

         CALL ${prefix}$_release(work)
         CALL timestop(handle)

      END SUBROUTINE screened_interaction_${suffix}$
   #:endfor

END MODULE mbpt_screening
