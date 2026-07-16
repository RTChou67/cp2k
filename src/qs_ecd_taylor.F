!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Normalized Taylor coefficient algebra for finite-wave-vector ECD.
! **************************************************************************************************
MODULE qs_ecd_taylor
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_solve
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_taylor'

   TYPE, PUBLIC :: ecd_taylor_scalar_type
      INTEGER                                      :: order = -1
      COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE  :: coeff
   END TYPE ecd_taylor_scalar_type

   TYPE, PUBLIC :: ecd_taylor_cfm_type
      INTEGER                                      :: order = -1
      TYPE(cp_cfm_type), DIMENSION(:), ALLOCATABLE :: coeff
   END TYPE ecd_taylor_cfm_type

   PUBLIC :: ecd_taylor_scalar_add, &
             ecd_taylor_scalar_conjugate, &
             ecd_taylor_scalar_create, &
             ecd_taylor_scalar_exp_iq, &
             ecd_taylor_scalar_inverse, &
             ecd_taylor_scalar_multiply, &
             ecd_taylor_scalar_release, &
             ecd_taylor_scalar_set_constant
   PUBLIC :: ecd_taylor_cfm_add, &
             ecd_taylor_cfm_conjugate, &
             ecd_taylor_cfm_create, &
             ecd_taylor_cfm_multiply, &
             ecd_taylor_cfm_release, &
             ecd_taylor_cfm_set_zero, &
             ecd_taylor_cfm_solve

CONTAINS

! **************************************************************************************************
!> \brief Allocate coefficients with bounds 0:order and initialize them to zero.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_create(series, order)
      TYPE(ecd_taylor_scalar_type), INTENT(OUT)          :: series
      INTEGER, INTENT(IN)                                :: order

      CPASSERT(order >= 0)

      series%order = order
      ALLOCATE (series%coeff(0:order))
      series%coeff = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

   END SUBROUTINE ecd_taylor_scalar_create

! **************************************************************************************************
!> \brief Release coefficient storage.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_release(series)
      TYPE(ecd_taylor_scalar_type), INTENT(INOUT)        :: series

      IF (ALLOCATED(series%coeff)) DEALLOCATE (series%coeff)
      series%order = -1

   END SUBROUTINE ecd_taylor_scalar_release

! **************************************************************************************************
!> \brief Set a q-independent complex scalar.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_set_constant(series, value, order)
      TYPE(ecd_taylor_scalar_type), INTENT(OUT)          :: series
      COMPLEX(KIND=dp), INTENT(IN)                       :: value
      INTEGER, INTENT(IN)                                :: order

      CALL ecd_taylor_scalar_create(series, order)
      series%coeff(0) = value

   END SUBROUTINE ecd_taylor_scalar_set_constant

! **************************************************************************************************
!> \brief Add two normalized Taylor series.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_add(sum_series, left, right)
      TYPE(ecd_taylor_scalar_type), INTENT(OUT)          :: sum_series
      TYPE(ecd_taylor_scalar_type), INTENT(IN)           :: left, right

      CALL ecd_taylor_scalar_assert_compatible(left, right)
      CALL ecd_taylor_scalar_create(sum_series, left%order)
      sum_series%coeff = left%coeff + right%coeff

   END SUBROUTINE ecd_taylor_scalar_add

! **************************************************************************************************
!> \brief Form the Cauchy product of two normalized Taylor series.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_multiply(product, left, right)
      TYPE(ecd_taylor_scalar_type), INTENT(OUT)          :: product
      TYPE(ecd_taylor_scalar_type), INTENT(IN)           :: left, right

      INTEGER                                            :: iorder, jorder

      CALL ecd_taylor_scalar_assert_compatible(left, right)
      CALL ecd_taylor_scalar_create(product, left%order)

      DO iorder = 0, product%order
         DO jorder = 0, iorder
            product%coeff(iorder) = product%coeff(iorder) + &
                                    left%coeff(jorder)*right%coeff(iorder - jorder)
         END DO
      END DO

   END SUBROUTINE ecd_taylor_scalar_multiply

! **************************************************************************************************
!> \brief Conjugate every normalized Taylor coefficient.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_conjugate(conjugated, input)
      TYPE(ecd_taylor_scalar_type), INTENT(OUT)          :: conjugated
      TYPE(ecd_taylor_scalar_type), INTENT(IN)           :: input

      CALL ecd_taylor_scalar_assert_valid(input)
      CALL ecd_taylor_scalar_create(conjugated, input%order)
      conjugated%coeff = CONJG(input%coeff)

   END SUBROUTINE ecd_taylor_scalar_conjugate

! **************************************************************************************************
!> \brief Form the normalized Taylor series for the reciprocal of input.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_inverse(inverse, input)
      TYPE(ecd_taylor_scalar_type), INTENT(OUT)          :: inverse
      TYPE(ecd_taylor_scalar_type), INTENT(IN)           :: input

      INTEGER                                            :: iorder, jorder

      CALL ecd_taylor_scalar_assert_valid(input)
      IF (ABS(input%coeff(0)) <= TINY(1.0_dp)) THEN
         CPABORT("ECD| Cannot invert a Taylor series with a zero constant coefficient")
      END IF

      CALL ecd_taylor_scalar_create(inverse, input%order)
      inverse%coeff(0) = 1.0_dp/input%coeff(0)
      DO iorder = 1, inverse%order
         DO jorder = 1, iorder
            inverse%coeff(iorder) = inverse%coeff(iorder) + &
                                    input%coeff(jorder)*inverse%coeff(iorder - jorder)
         END DO
         inverse%coeff(iorder) = -inverse%coeff(iorder)/input%coeff(0)
      END DO

   END SUBROUTINE ecd_taylor_scalar_inverse

! **************************************************************************************************
!> \brief Form exp(i*q*phase) through the requested normalized Taylor order.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_exp_iq(phase_series, phase, order)
      TYPE(ecd_taylor_scalar_type), INTENT(OUT)          :: phase_series
      REAL(KIND=dp), INTENT(IN)                          :: phase
      INTEGER, INTENT(IN)                                :: order

      COMPLEX(KIND=dp)                                   :: iphase
      INTEGER                                            :: iorder

      CALL ecd_taylor_scalar_create(phase_series, order)
      phase_series%coeff(0) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      iphase = CMPLX(0.0_dp, phase, KIND=dp)
      DO iorder = 1, order
         phase_series%coeff(iorder) = phase_series%coeff(iorder - 1)*iphase/ &
                                      REAL(iorder, KIND=dp)
      END DO

   END SUBROUTINE ecd_taylor_scalar_exp_iq

! **************************************************************************************************
!> \brief Allocate a normalized Taylor series of distributed complex matrices.
!> \param series Matrix coefficient container
!> \param matrix_struct CP2K distributed-matrix structure
!> \param order Highest Taylor order
!> \param nrow Optional global row count
!> \param ncol Optional global column count
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_create(series, matrix_struct, order, nrow, ncol)
      TYPE(ecd_taylor_cfm_type), INTENT(OUT)             :: series
      TYPE(cp_fm_struct_type), INTENT(IN), TARGET        :: matrix_struct
      INTEGER, INTENT(IN)                                :: order
      INTEGER, INTENT(IN), OPTIONAL                      :: nrow, ncol

      INTEGER                                            :: iorder

      CPASSERT(order >= 0)

      series%order = order
      ALLOCATE (series%coeff(0:order))
      DO iorder = 0, order
         IF (PRESENT(nrow)) THEN
            IF (PRESENT(ncol)) THEN
               CALL cp_cfm_create(series%coeff(iorder), matrix_struct, nrow=nrow, ncol=ncol, &
                                  name="ECD Taylor coefficient", set_zero=.TRUE.)
            ELSE
               CALL cp_cfm_create(series%coeff(iorder), matrix_struct, nrow=nrow, &
                                  name="ECD Taylor coefficient", set_zero=.TRUE.)
            END IF
         ELSE IF (PRESENT(ncol)) THEN
            CALL cp_cfm_create(series%coeff(iorder), matrix_struct, ncol=ncol, &
                               name="ECD Taylor coefficient", set_zero=.TRUE.)
         ELSE
            CALL cp_cfm_create(series%coeff(iorder), matrix_struct, &
                               name="ECD Taylor coefficient", set_zero=.TRUE.)
         END IF
      END DO

   END SUBROUTINE ecd_taylor_cfm_create

! **************************************************************************************************
!> \brief Release all distributed complex matrix coefficients.
!> \param series Matrix coefficient container
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_release(series)
      TYPE(ecd_taylor_cfm_type), INTENT(INOUT)           :: series

      INTEGER                                            :: iorder

      IF (ALLOCATED(series%coeff)) THEN
         DO iorder = 0, SIZE(series%coeff) - 1
            IF (ASSOCIATED(series%coeff(iorder)%local_data)) THEN
               CALL cp_cfm_release(series%coeff(iorder))
            END IF
         END DO
         DEALLOCATE (series%coeff)
      END IF
      series%order = -1

   END SUBROUTINE ecd_taylor_cfm_release

! **************************************************************************************************
!> \brief Set every matrix coefficient to zero.
!> \param series Matrix coefficient container
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_set_zero(series)
      TYPE(ecd_taylor_cfm_type), INTENT(INOUT)           :: series

      INTEGER                                            :: iorder

      CALL ecd_taylor_cfm_assert_valid(series)
      DO iorder = 0, series%order
         CALL cp_cfm_set_all(series%coeff(iorder), CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      END DO

   END SUBROUTINE ecd_taylor_cfm_set_zero

! **************************************************************************************************
!> \brief Add two normalized Taylor matrix series.
!> \param sum_series Sum of left and right
!> \param left Left matrix coefficient series
!> \param right Right matrix coefficient series
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_add(sum_series, left, right)
      TYPE(ecd_taylor_cfm_type), INTENT(OUT)             :: sum_series
      TYPE(ecd_taylor_cfm_type), INTENT(IN)              :: left, right

      INTEGER                                            :: iorder

      CALL ecd_taylor_cfm_assert_compatible(left, right)
      CALL ecd_taylor_cfm_create(sum_series, left%coeff(0)%matrix_struct, left%order)
      DO iorder = 0, left%order
         CALL cp_cfm_to_cfm(left%coeff(iorder), sum_series%coeff(iorder))
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, KIND=dp), sum_series%coeff(iorder), &
                                   CMPLX(1.0_dp, 0.0_dp, KIND=dp), right%coeff(iorder))
      END DO

   END SUBROUTINE ecd_taylor_cfm_add

! **************************************************************************************************
!> \brief Conjugate each normalized Taylor matrix coefficient.
!> \param conjugated Conjugated coefficient series
!> \param input Input matrix coefficient series
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_conjugate(conjugated, input)
      TYPE(ecd_taylor_cfm_type), INTENT(OUT)             :: conjugated
      TYPE(ecd_taylor_cfm_type), INTENT(IN)              :: input

      INTEGER                                            :: iorder

      CALL ecd_taylor_cfm_assert_valid(input)
      CALL ecd_taylor_cfm_create(conjugated, input%coeff(0)%matrix_struct, input%order)
      DO iorder = 0, input%order
         conjugated%coeff(iorder)%local_data = CONJG(input%coeff(iorder)%local_data)
      END DO

   END SUBROUTINE ecd_taylor_cfm_conjugate

! **************************************************************************************************
!> \brief Form the Cauchy product of two normalized Taylor matrix series.
!> \param product Product of left and right
!> \param left Left matrix coefficient series
!> \param right Right matrix coefficient series
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_multiply(product, left, right)
      TYPE(ecd_taylor_cfm_type), INTENT(OUT)             :: product
      TYPE(ecd_taylor_cfm_type), INTENT(IN)              :: left, right

      INTEGER                                            :: iorder, jorder, ncol, ninner, nrow, &
                                                            right_nrow

      CALL ecd_taylor_cfm_assert_valid(left)
      CALL ecd_taylor_cfm_assert_valid(right)
      CPASSERT(left%order == right%order)
      CALL cp_cfm_get_info(left%coeff(0), nrow_global=nrow, ncol_global=ninner)
      CALL cp_cfm_get_info(right%coeff(0), nrow_global=right_nrow, ncol_global=ncol)
      CPASSERT(ninner == right_nrow)

      CALL ecd_taylor_cfm_create(product, left%coeff(0)%matrix_struct, left%order, nrow=nrow, ncol=ncol)
      DO iorder = 0, product%order
         DO jorder = 0, iorder
            CALL cp_cfm_gemm("N", "N", nrow, ncol, ninner, CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                             left%coeff(jorder), right%coeff(iorder - jorder), &
                             CMPLX(1.0_dp, 0.0_dp, KIND=dp), product%coeff(iorder))
         END DO
      END DO

   END SUBROUTINE ecd_taylor_cfm_multiply

! **************************************************************************************************
!> \brief Solve a normalized Taylor-series linear system coefficient by coefficient.
!> \param solution Series satisfying coefficient*solution=right_hand_side
!> \param coefficient Square coefficient-matrix series
!> \param right_hand_side Right-hand-side series
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_solve(solution, coefficient, right_hand_side)
      TYPE(ecd_taylor_cfm_type), INTENT(OUT)             :: solution
      TYPE(ecd_taylor_cfm_type), INTENT(IN)              :: coefficient, right_hand_side

      INTEGER                                            :: iorder, jorder, ncol, ninner, nrow, &
                                                            rhs_nrow
      TYPE(cp_cfm_type)                                  :: coefficient_zero

      CALL ecd_taylor_cfm_assert_valid(coefficient)
      CALL ecd_taylor_cfm_assert_valid(right_hand_side)
      CPASSERT(coefficient%order == right_hand_side%order)
      CALL cp_cfm_get_info(coefficient%coeff(0), nrow_global=nrow, ncol_global=ninner)
      CALL cp_cfm_get_info(right_hand_side%coeff(0), nrow_global=rhs_nrow, ncol_global=ncol)
      CPASSERT(nrow == ninner)
      CPASSERT(nrow == rhs_nrow)

      CALL ecd_taylor_cfm_create(solution, coefficient%coeff(0)%matrix_struct, coefficient%order, &
                                 nrow=nrow, ncol=ncol)
      DO iorder = 0, solution%order
         CALL cp_cfm_to_cfm(right_hand_side%coeff(iorder), solution%coeff(iorder))
         DO jorder = 1, iorder
            CALL cp_cfm_gemm("N", "N", nrow, ncol, ninner, CMPLX(-1.0_dp, 0.0_dp, KIND=dp), &
                             coefficient%coeff(jorder), solution%coeff(iorder - jorder), &
                             CMPLX(1.0_dp, 0.0_dp, KIND=dp), solution%coeff(iorder))
         END DO

         CALL cp_cfm_create(coefficient_zero, coefficient%coeff(0)%matrix_struct, &
                            name="ECD Taylor coefficient zero")
         CALL cp_cfm_to_cfm(coefficient%coeff(0), coefficient_zero)
         CALL cp_cfm_solve(coefficient_zero, solution%coeff(iorder))
         CALL cp_cfm_release(coefficient_zero)
      END DO

   END SUBROUTINE ecd_taylor_cfm_solve

! **************************************************************************************************
!> \brief Check allocation and normalized coefficient bounds.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_assert_valid(series)
      TYPE(ecd_taylor_scalar_type), INTENT(IN)           :: series

      CPASSERT(series%order >= 0)
      CPASSERT(ALLOCATED(series%coeff))
      CPASSERT(LBOUND(series%coeff, 1) == 0)
      CPASSERT(UBOUND(series%coeff, 1) == series%order)

   END SUBROUTINE ecd_taylor_scalar_assert_valid

! **************************************************************************************************
!> \brief Check that two series use the same normalized Taylor order.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_scalar_assert_compatible(left, right)
      TYPE(ecd_taylor_scalar_type), INTENT(IN)           :: left, right

      CALL ecd_taylor_scalar_assert_valid(left)
      CALL ecd_taylor_scalar_assert_valid(right)
      CPASSERT(left%order == right%order)

   END SUBROUTINE ecd_taylor_scalar_assert_compatible

! **************************************************************************************************
!> \brief Check allocation, normalized bounds, and matrix ownership.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_assert_valid(series)
      TYPE(ecd_taylor_cfm_type), INTENT(IN)              :: series

      INTEGER                                            :: iorder

      CPASSERT(series%order >= 0)
      CPASSERT(ALLOCATED(series%coeff))
      CPASSERT(LBOUND(series%coeff, 1) == 0)
      CPASSERT(UBOUND(series%coeff, 1) == series%order)
      DO iorder = 0, series%order
         CPASSERT(ASSOCIATED(series%coeff(iorder)%matrix_struct))
         CPASSERT(ASSOCIATED(series%coeff(iorder)%local_data))
      END DO

   END SUBROUTINE ecd_taylor_cfm_assert_valid

! **************************************************************************************************
!> \brief Check equal order and distributed-matrix shape for coefficientwise operations.
! **************************************************************************************************
   SUBROUTINE ecd_taylor_cfm_assert_compatible(left, right)
      TYPE(ecd_taylor_cfm_type), INTENT(IN)              :: left, right

      INTEGER                                            :: iorder, left_ncol, left_nrow, &
                                                            right_ncol, right_nrow

      CALL ecd_taylor_cfm_assert_valid(left)
      CALL ecd_taylor_cfm_assert_valid(right)
      CPASSERT(left%order == right%order)
      DO iorder = 0, left%order
         CALL cp_cfm_get_info(left%coeff(iorder), nrow_global=left_nrow, ncol_global=left_ncol)
         CALL cp_cfm_get_info(right%coeff(iorder), nrow_global=right_nrow, ncol_global=right_ncol)
         CPASSERT(left_nrow == right_nrow)
         CPASSERT(left_ncol == right_ncol)
      END DO

   END SUBROUTINE ecd_taylor_cfm_assert_compatible

END MODULE qs_ecd_taylor
