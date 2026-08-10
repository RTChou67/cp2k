!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program for atomistic simulations                                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common AO-to-MO operator projection for electronic circular dichroism.
! **************************************************************************************************
MODULE qs_ecd_operator_utils
   USE cp_dbcsr_api, ONLY: dbcsr_p_type
   USE cp_dbcsr_operations, ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_struct, ONLY: cp_fm_struct_type
   USE cp_fm_types, ONLY: cp_fm_create, &
                          cp_fm_get_info, &
                          cp_fm_release, &
                          cp_fm_type
   USE kinds, ONLY: dp
   USE parallel_gemm_api, ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ecd_project_operator_to_mo

CONTAINS

! **************************************************************************************************
!> \brief Project a real sparse AO operator into the full real Gamma-point MO basis.
!> The two matrix products evaluate O^MO = C^T O^AO C.  For length-gauge E1,
!> O^AO is r_lambda-R0_lambda and mu_lambda=-O^MO_lambda.  The same projection
!> is used for the real AO angular-momentum and velocity-generator matrices.
!> \param operator_ao AO-space sparse operator
!> \param mo_coeff Real MO coefficient matrix
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

END MODULE qs_ecd_operator_utils
