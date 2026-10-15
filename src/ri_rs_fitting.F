!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared numerical operations for computing the RI-RS matrix Z_lP.
!> \par History
!>      09.2026 created Jan Wilhelm
! **************************************************************************************************
MODULE ri_rs_fitting
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_put_block,&
                                              dbcsr_type
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_solve
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_fitting'
   REAL(KIND=dp), PARAMETER, PRIVATE :: jacobi_floor = 1.0E-16_dp

   PUBLIC :: solve_ri_rs_fit, &
             store_Z_lP_columns

CONTAINS

! **************************************************************************************************
!> \brief Computes the regularized RI-RS fit from an AO collocation matrix and a replicated RHS.
!>
!>          D_ll' = [Σ_μ ϕ_μ(r_l) ϕ_μ(r_l')]²,
!>          S_ll = 1/sqrt(max(D_ll, jacobi_floor)),
!>          (S D S + λ I) Y = S R,   Z = S Y.
!>
!>        The regularization is defined in the Jacobi-scaled space, corresponding to
!>        (D + λ S^-2) Z = R. Both LAPACK and ScaLAPACK use this definition.
!> \param phi_local AO values on the fitting grid, replicated within the solve subgroup
!> \param d_lp replicated RHS R on entry and fitted coefficients Z on return
!> \param tikhonov regularization parameter λ in the Jacobi-scaled space
!> \param para_env_sub optional subgroup communicator for the distributed solve
!> \param blacs_env_sub optional BLACS context, supplied together with para_env_sub
! **************************************************************************************************
   SUBROUTINE solve_ri_rs_fit(phi_local, d_lp, tikhonov, para_env_sub, blacs_env_sub)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_local
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_sub
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_sub

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_ri_rs_fit'

      INTEGER                                            :: handle, handle2, i, info, j, n_ao, &
                                                            n_grid, n_rhs
      REAL(KIND=dp), ALLOCATABLE                         :: gram(:, :), jacobi_scale(:)

      CALL timeset(routineN, handle)

      n_grid = SIZE(phi_local, 1)
      n_ao = SIZE(phi_local, 2)
      n_rhs = SIZE(d_lp, 2)
      CPASSERT(SIZE(d_lp, 1) == n_grid)
      CPASSERT(PRESENT(para_env_sub) .EQV. PRESENT(blacs_env_sub))
      ALLOCATE (jacobi_scale(n_grid))

      IF (PRESENT(para_env_sub)) THEN
         ! Only the row norms are replicated; the distributed solve builds its own Gram tiles.
         !$OMP PARALLEL DO DEFAULT(NONE) &
         !$OMP SHARED(n_grid, phi_local, jacobi_scale) &
         !$OMP PRIVATE(i) SCHEDULE(STATIC)
         DO i = 1, n_grid
            jacobi_scale(i) = DOT_PRODUCT(phi_local(i, :), phi_local(i, :))
         END DO
         !$OMP END PARALLEL DO
      ELSE
         ALLOCATE (gram(n_grid, n_grid))
         CALL timeset(routineN//'_dsyrk', handle2)
         CALL dsyrk('L', 'N', n_grid, n_ao, 1.0_dp, phi_local, n_grid, &
                    0.0_dp, gram, n_grid)
         CALL timestop(handle2)
         DO i = 1, n_grid
            jacobi_scale(i) = gram(i, i)
         END DO
      END IF

      ! Convert the AO row norms squared to S_ll = 1/max(sum_mu phi_lmu**2, sqrt(jacobi_floor)).
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(n_grid, jacobi_scale) &
      !$OMP PRIVATE(i) &
      !$OMP SCHEDULE(STATIC)
      DO i = 1, n_grid
         jacobi_scale(i) = 1.0_dp/MAX(jacobi_scale(i), SQRT(jacobi_floor))
      END DO
      !$OMP END PARALLEL DO

      CALL scale_rows_by_diag(d_lp, jacobi_scale)
      IF (PRESENT(para_env_sub)) THEN
         CALL solve_D_lp_distributed(phi_local, jacobi_scale, d_lp, tikhonov, &
                                     para_env_sub, blacs_env_sub)
      ELSE
         !$OMP PARALLEL DO DEFAULT(NONE) &
         !$OMP SHARED(n_grid, gram, jacobi_scale, tikhonov) &
         !$OMP PRIVATE(j, i) SCHEDULE(DYNAMIC)
         DO j = 1, n_grid
            DO i = j, n_grid
               gram(i, j) = gram(i, j)**2*jacobi_scale(i)*jacobi_scale(j)
            END DO
            gram(j, j) = gram(j, j) + tikhonov
         END DO
         !$OMP END PARALLEL DO
         CALL timeset(routineN//'_dpotrf', handle2)
         CALL dpotrf('L', n_grid, gram, n_grid, info)
         CALL timestop(handle2)
         IF (info /= 0) CPABORT('RI-RS Cholesky factorization failed')
         CALL timeset(routineN//'_dpotrs', handle2)
         CALL dpotrs('L', n_grid, n_rhs, gram, n_grid, d_lp, n_grid, info)
         CALL timestop(handle2)
         IF (info /= 0) CPABORT('RI-RS Cholesky solve failed')
      END IF
      CALL scale_rows_by_diag(d_lp, jacobi_scale)

      CALL timestop(handle)

   END SUBROUTINE solve_ri_rs_fit

! **************************************************************************************************
!> \brief Multiplies every matrix row by the corresponding diagonal entry:
!>        A(l, :) <- d_l A(l, :).
!> \param matrix ...
!> \param diagonal ...
! **************************************************************************************************
   SUBROUTINE scale_rows_by_diag(matrix, diagonal)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: diagonal

      CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_rows_by_diag'

      INTEGER                                            :: handle, i, j

      CALL timeset(routineN, handle)
      CPASSERT(SIZE(diagonal) == SIZE(matrix, 1))

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(matrix, diagonal) &
      !$OMP PRIVATE(j, i) &
      !$OMP SCHEDULE(STATIC)
      DO j = 1, SIZE(matrix, 2)
         DO i = 1, SIZE(matrix, 1)
            matrix(i, j) = matrix(i, j)*diagonal(i)
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE scale_rows_by_diag

! **************************************************************************************************
!> \brief Stores dense Z_lP columns in the distributed block-sparse matrix.
!> \param mat_Z_lP ...
!> \param Z_local ...
!> \param local_grid_idx ...
!> \param n_local_grid ...
!> \param n_loc_ri ...
!> \param atom_P ...
!> \param r_blk_sizes ...
!> \param row_offset ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE store_Z_lP_columns(mat_Z_lP, Z_local, local_grid_idx, n_local_grid, n_loc_ri, &
                                 atom_P, r_blk_sizes, row_offset, eps_filter)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_Z_lP
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Z_local
      INTEGER, DIMENSION(:), INTENT(IN)                  :: local_grid_idx
      INTEGER, INTENT(IN)                                :: n_local_grid, n_loc_ri, atom_P
      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes, row_offset
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'store_Z_lP_columns'

      INTEGER                                            :: current_chunk_size, g_pt, handle, i_blk, &
                                                            loc_ptr, r_end, r_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Z_blk

      CALL timeset(routineN, handle)

      ALLOCATE (Z_blk(MAXVAL(r_blk_sizes), n_loc_ri))
      loc_ptr = 1

      DO i_blk = 1, SIZE(r_blk_sizes)
         r_start = row_offset(i_blk) + 1
         r_end = row_offset(i_blk) + r_blk_sizes(i_blk)
         current_chunk_size = r_blk_sizes(i_blk)
         Z_blk = 0.0_dp

         DO WHILE (loc_ptr <= n_local_grid)
            g_pt = local_grid_idx(loc_ptr)
            IF (g_pt > r_end) EXIT
            Z_blk(g_pt - r_start + 1, 1:n_loc_ri) = Z_local(loc_ptr, 1:n_loc_ri)
            loc_ptr = loc_ptr + 1
         END DO

         IF (MAXVAL(ABS(Z_blk(1:current_chunk_size, 1:n_loc_ri))) > eps_filter) THEN
            CALL dbcsr_put_block(mat_Z_lP, row=i_blk, col=atom_P, &
                                 block=Z_blk(1:current_chunk_size, 1:n_loc_ri))
         END IF
      END DO

      DEALLOCATE (Z_blk)

      CALL timestop(handle)

   END SUBROUTINE store_Z_lP_columns

! **************************************************************************************************
!> \brief Solves the Jacobi-scaled, regularized RI-RS system with ScaLAPACK.
!>
!>        Each rank builds its block-cyclic part of
!>
!>          D'_ll' = d_l [Σ_μ ϕ_μ(r_l) ϕ_μ(r_l')]² d_l' + λ δ_ll'.
!>
!>        The Cholesky solution is gathered into the replicated right-hand side.
!> \param phi_local ...
!> \param d_vec ...
!> \param d_lp ...
!> \param tikhonov ...
!> \param para_env_sub ...
!> \param blacs_env_sub ...
! **************************************************************************************************
   SUBROUTINE solve_D_lp_distributed(phi_local, d_vec, d_lp, tikhonov, para_env_sub, blacs_env_sub)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_local
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: d_vec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_D_lp_distributed'

      INTEGER                                            :: handle, i_loc, ig, info, j_loc, jg, &
                                                            n_ao, n_loc, n_rhs, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_b, fm_struct_D
      TYPE(cp_fm_type)                                   :: fm_b, fm_D

      CALL timeset(routineN, handle)
      n_loc = SIZE(phi_local, 1)
      n_ao = SIZE(phi_local, 2)
      n_rhs = SIZE(d_lp, 2)

      NULLIFY (fm_struct_D, fm_struct_b)
      CALL cp_fm_struct_create(fm_struct_D, para_env=para_env_sub, &
                               context=blacs_env_sub, &
                               nrow_global=n_loc, ncol_global=n_loc)
      CALL cp_fm_struct_create(fm_struct_b, para_env=para_env_sub, &
                               context=blacs_env_sub, &
                               nrow_global=n_loc, ncol_global=n_rhs)
      CALL cp_fm_create(fm_D, fm_struct_D)
      CALL cp_fm_create(fm_b, fm_struct_b)

      CALL cp_fm_get_info(fm_D, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data)

      IF (nrow_local > 0 .AND. ncol_local > 0) THEN
         BLOCK
            INTEGER, PARAMETER :: ntile = 1024
            INTEGER :: ib, ie, jb, je, mb, kb, ti, tj, handle_dgemm
            REAL(KIND=dp), ALLOCATABLE :: gram_t(:, :), phi_cols_t(:, :), phi_rows_t(:, :)
            ALLOCATE (phi_rows_t(ntile, n_ao), phi_cols_t(n_ao, ntile), gram_t(ntile, ntile))
            DO ib = 1, nrow_local, ntile
               ie = MIN(ib + ntile - 1, nrow_local)
               mb = ie - ib + 1
               !$OMP PARALLEL DO DEFAULT(NONE) &
               !$OMP SHARED(mb, n_ao, phi_rows_t, phi_local, row_indices, ib) &
               !$OMP PRIVATE(ti, j_loc) SCHEDULE(STATIC)
               DO j_loc = 1, n_ao
                  DO ti = 1, mb
                     phi_rows_t(ti, j_loc) = phi_local(row_indices(ib + ti - 1), j_loc)
                  END DO
               END DO
               !$OMP END PARALLEL DO
               DO jb = 1, ncol_local, ntile
                  je = MIN(jb + ntile - 1, ncol_local)
                  kb = je - jb + 1
                  !$OMP PARALLEL DO DEFAULT(NONE) &
                  !$OMP SHARED(kb, n_ao, phi_cols_t, phi_local, col_indices, jb) &
                  !$OMP PRIVATE(tj, i_loc) SCHEDULE(STATIC)
                  DO tj = 1, kb
                     DO i_loc = 1, n_ao
                        phi_cols_t(i_loc, tj) = phi_local(col_indices(jb + tj - 1), i_loc)
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
                  CALL timeset(routineN//'_dgemm', handle_dgemm)
                  CALL dgemm('N', 'N', mb, kb, n_ao, &
                             1.0_dp, phi_rows_t, ntile, phi_cols_t, n_ao, &
                             0.0_dp, gram_t, ntile)
                  CALL timestop(handle_dgemm)
                  !$OMP PARALLEL DO DEFAULT(NONE) &
                  !$OMP SHARED(mb, kb, gram_t, d_vec, row_indices, col_indices, ib, jb) &
                  !$OMP SHARED(local_data, tikhonov) &
                  !$OMP PRIVATE(ti, tj, ig, jg) SCHEDULE(STATIC)
                  DO tj = 1, kb
                     jg = col_indices(jb + tj - 1)
                     DO ti = 1, mb
                        ig = row_indices(ib + ti - 1)
                        local_data(ib + ti - 1, jb + tj - 1) = &
                           gram_t(ti, tj)*gram_t(ti, tj)*d_vec(ig)*d_vec(jg)
                        IF (ig == jg) THEN
                           local_data(ib + ti - 1, jb + tj - 1) = &
                              local_data(ib + ti - 1, jb + tj - 1) + tikhonov
                        END IF
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
               END DO
            END DO
            DEALLOCATE (phi_rows_t, phi_cols_t, gram_t)
         END BLOCK
      END IF

      CALL cp_fm_set_submatrix(fm_b, d_lp)

      CALL cp_fm_cholesky_decompose(fm_D, n=n_loc, info_out=info)
      IF (info /= 0) CPABORT("pdpotrf failed in solve_D_lp_distributed")

      CALL cp_fm_cholesky_solve(fm_D, fm_b, n=n_loc, info_out=info)
      IF (info /= 0) CPABORT("pdpotrs failed in solve_D_lp_distributed")

      CALL cp_fm_get_submatrix(fm_b, d_lp)

      CALL cp_fm_release(fm_D)
      CALL cp_fm_release(fm_b)
      CALL cp_fm_struct_release(fm_struct_D)
      CALL cp_fm_struct_release(fm_struct_b)

      CALL timestop(handle)

   END SUBROUTINE solve_D_lp_distributed

END MODULE ri_rs_fitting
