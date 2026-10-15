!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using RI-RS Approximation for molecules
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE gw_ri_rs_large_cell_Gamma
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_binary_read, dbcsr_binary_write, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_info, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output
   USE gw_ri_rs_kernels,                ONLY: compute_Sigma_c,&
                                              compute_Sigma_x,&
                                              get_mat_chi_Gamma_tau
   USE gw_tensor_large_cell_Gamma,      ONLY: compute_QP_energies,&
                                              delete_unnecessary_files,&
                                              get_W_MIC
   USE gw_utils,                        ONLY: de_init_bs_env
   USE gw_utils_compute_integrals,      ONLY: build_3c_integral_block_ctx,&
                                              gw_3c_ctx_create,&
                                              gw_3c_ctx_release,&
                                              gw_3c_ctx_type,&
                                              gw_3c_ws_create,&
                                              gw_3c_ws_release,&
                                              gw_3c_ws_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE ri_rs_fitting,                   ONLY: build_gram_jacobi_blas,&
                                              scale_rows_by_diag,&
                                              store_Z_lP_columns
   USE ri_rs_grid_setup,                ONLY: setup_ri_rs_grid
   USE ri_rs_matrices,                  ONLY: atomic_basis_at_grid_point
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points,&
                                              precompute_ri_rs_radii
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_large_cell_Gamma'

   PUBLIC :: gw_calc_ri_rs_large_cell_Gamma

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism for molecules
!> \param qs_env ...
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************

   SUBROUTINE gw_calc_ri_rs_large_cell_Gamma(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_ri_rs_large_cell_Gamma'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

       !!========================================================================
       !! 0. Precompute AO and RI Radii
       !!    Per-atom cutoff radii from the most diffuse Gaussian primitives in
       !!    the AO and RI auxiliary basis sets. Stored in bs_env%ri_rs%
       !!    radius_ao_per_atom and radius_ri_per_atom, used for sphere-cutoff
       !!    and phi_local screening.
       !!========================================================================
      CALL precompute_ri_rs_radii(bs_env%ri_rs, bs_env%basis_set_AO, bs_env%basis_set_RI, &
                                  bs_env%eps_filter, bs_env%unit_nr)

       !!========================================================================
       !! 1. Grid Generation for RI-RS
       !!    (Modified Lebedev grids from Ivan Duchemin and Xavier Blase)
       !!    Generate flattened 1D array of grid points for RI-RS.
       !!    Equation: r_g(k) = R_A + r_g(A)
       !!========================================================================
      CALL setup_ri_rs_grid(bs_env, bs_env%ri_rs%grid_points)

       !!========================================================================
       !! 2. Atomic Basis Evaluation
       !!    Compute values of spherical atomic basis functions at grid points.
       !!    Expression: ϕ_μl = ϕ_μ(r_l) (mat_phi_mu_l)
       !!========================================================================
      CALL atomic_basis_at_grid_point(bs_env%ri_rs, bs_env%basis_set_AO, &
                                      bs_env%mat_ao_ao%matrix, bs_env%para_env, bs_env%eps_filter, &
                                      bs_env%dbcsr_msg_elem_limit, bs_env%unit_nr, periodic_images=.TRUE.)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, *) 'Done with evaluation of phi'

       !!========================================================================
       !! 3. Compute RI-RS Coefficients (Z_lp)
       !!    Solve the regularized system for each atom P, where the grid domain
       !!    is restricted to r_l within a cutoff distance of atom P:
       !!    a. D_ll' = [ Σ_μ ϕ_μ(r_l) ϕ_μ(r_l') ]^2 (Equation 13)
       !!    b. D_lP  = Σ_{μν} ϕ_μ(r_l) ϕ_ν(r_l) (μν|P) (Equation 15)
       !!    c. Conditioning:
       !!       Dvec_l   = 1 / sqrt(D_ll)  (Diagonal scaling vector)
       !!       D'_ll' = Dvec_l * D_ll' * Dvec_l' + λδ_ll'
       !!       D'_lP  = Dvec_l * D_lP
       !!    d. Solve: Σ_l' D'_ll' * Z'_l'P = D'_lP (Equation 14)
       !!    e. Rescale: Z_lP = Z'_lP * Dvec_l   (Z_lP stored in mat_Z_lP)
       !!========================================================================
      CALL compute_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                        bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)
      bs_env%ri_rs%grid_built = .TRUE.

       !!========================================================================
       !! 4. Compute Independent-Particle Polarizability (χ)
       !!    G^occ_µλ(i|τ|)  = sum_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^vir_µλ(i|τ|)  = sum_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^occ_ll'(i|τ|) = sum_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')
       !!    G^vir_ll'(i|τ|) = sum_µν ϕ_µ(r_l) G^vir_µν ϕ_ν(r_l')
       !!    χ_ll'(iτ)       = G^occ_ll'(i|τ|) * G^vir_ll'(i|τ|)
       !!    χ_PQ(iτ)        = sum_ll' Z_lP χ_ll'(iτ) Z_l'Q
       !!========================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP, cutoff=0.0_dp)

       !!========================================================================
       !! 5. Compute Screened Interaction (W^MIC)
       !!    χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W^MIC_PQ(iτ)
       !!========================================================================
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

       !!========================================================================
       !! 6. Compute Exact Exchange Self-Energy (Σ^x)
       !!    D_µν          = sum_n^occ C_µn C_νn
       !!    D_ll'         = sum_µν ϕ_µ(r_l) D_µν ϕ_ν(r_l')
       !!    V^trunc_ll'   = sum_PQ Z_lP V^trunc_PQ Z_l'Q
       !!    Σ^x_ll'       = D_ll' * V^trunc_ll'
       !!    Σ^x_λσ(k=0)   = -sum_ll' ϕ_λ(r_l) Σ^x_ll' ϕ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma, cutoff=0.0_dp)

       !!========================================================================
       !! 7. Compute Correlation Self-Energy (Σ^c)
       !!    W^MIC_ll'(iτ) =  sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
       !!    Σ^c_ll'(iτ)   = -G^occ_ll'(i|τ|) * W^MIC_ll'(iτ), for τ < 0
       !!    Σ^c_ll'(iτ)   =  G^vir_ll'(i|τ|) * W^MIC_ll'(iτ), for τ > 0
       !!    Σ^c_λσ(iτ)    =  sum_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_c(bs_env, fm_W_time, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time, cutoff=0.0_dp)
      CALL cp_fm_release(fm_W_time)
      CALL delete_unnecessary_files(bs_env)

       !!========================================================================
       !! 8. Compute Quasiparticle Energies
       !!    Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
       !!    ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ) + Σ^x_nn - v^xc_nn
       !!========================================================================
      CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_large_cell_Gamma

! **************************************************************************************************
!> \brief Compute RI-RS Coefficients (Z_lP)
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************

   SUBROUTINE compute_Z_lP(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP)

      ! Arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'compute_Z_lP'

      INTEGER :: atom_P, col_end, col_start, group_handle, handle, handle_dpotrf, handle_dpotrs, &
         i_blk, ikind, info, j, l, loc_idx, max_ao_size, max_loc_ri, n_ao_total, n_grid_total, &
         n_loc_ri, n_local_grid, natom, nkind, npcol_phi, num_grid_chunks, P_loop_atom, ri_atom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_grid_idx, row_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_ri, r_blk_sizes, ri_blk_sizes, &
                                                            row_dist_grid
      REAL(KIND=dp)                                      :: cutoff_ri, cutoff_ri_2, d_sP, dist2_min, &
                                                            r2_threshold, r_c, t1, t2, t3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cutoff_ri_per_atom, cutoff_ri_per_kind, &
                                                            d_vec_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: D_local, d_lp_local, phi_local, &
                                                            sphere_grid
      REAL(KIND=dp), DIMENSION(3)                        :: dist_vec_raw, pos_P
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist_phi, dist_Z
      TYPE(gw_3c_ctx_type)                               :: ctx_3c
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL get_qs_env(qs_env, para_env=para_env, particle_set=particle_set, input=input, cell=cell)
      atomic_kind_set => bs_env%ri_rs%atomic_kind_set

      natom = bs_env%n_atom
      n_ao_total = bs_env%i_ao_end_from_atom(natom)
      n_grid_total = bs_env%ri_rs%n_grid_points
      CPASSERT(SIZE(ri_rs_grid_points, 2) == n_grid_total)

      ! =========================================================================
      ! 1. SETUP DBCSR TOPOLOGY & EXACT OFFSETS
      ! =========================================================================
      CALL dbcsr_get_info(mat_phi_mu_l, row_blk_size=r_blk_sizes, distribution=dist_phi)
      CALL dbcsr_distribution_get(dist_phi, row_dist=row_dist_grid, &
                                  group=group_handle, npcols=npcol_phi)

      num_grid_chunks = SIZE(r_blk_sizes)

      ALLOCATE (row_offset(num_grid_chunks))
      row_offset(1) = 0
      DO i_blk = 2, num_grid_chunks
         row_offset(i_blk) = row_offset(i_blk - 1) + r_blk_sizes(i_blk - 1)
      END DO

      ALLOCATE (ri_blk_sizes(natom), col_dist_ri(natom))
      DO atom_P = 1, natom
         ri_blk_sizes(atom_P) = bs_env%i_RI_end_from_atom(atom_P) - bs_env%i_RI_start_from_atom(atom_P) + 1
         col_dist_ri(atom_P) = MOD(atom_P - 1, npcol_phi)
      END DO

      CALL dbcsr_distribution_new(dist_Z, template=dist_phi, row_dist=row_dist_grid, col_dist=col_dist_ri)

      IF (bs_env%ri_rs%Z_lP_exists) THEN
         CALL dbcsr_binary_read(filepath=TRIM(bs_env%prefix)//"Z_lP.matrix", &
                                distribution=dist_Z, &
                                matrix_new=mat_Z_lP)
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Read Z_lP from file ', ' Execution time', m_walltime() - t1, ' s'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
      ELSE

         CALL dbcsr_create(mat_Z_lP, name="mat_Z_lP", dist=dist_Z, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=r_blk_sizes, col_blk_size=ri_blk_sizes)

         max_ao_size = 0
         DO j = 1, natom
            max_ao_size = MAX(max_ao_size, bs_env%i_ao_end_from_atom(j) - bs_env%i_ao_start_from_atom(j) + 1)
         END DO
         max_loc_ri = MAXVAL(ri_blk_sizes)

         ! Per-atom RI-RS integration sphere:
         !   cutoff_ri(P) = r_c + r_AO(P)
         ! where r_c is the truncated-Coulomb cutoff of the RI metric.
         nkind = SIZE(atomic_kind_set)
         ALLOCATE (cutoff_ri_per_atom(natom))

         r_c = bs_env%ri_metric%cutoff_radius
         DO P_loop_atom = 1, natom
            cutoff_ri_per_atom(P_loop_atom) = &
               r_c + bs_env%ri_rs%radius_ao_per_atom(P_loop_atom)
         END DO

         ALLOCATE (cutoff_ri_per_kind(nkind))
         cutoff_ri_per_kind(:) = 0.0_dp
         IF (bs_env%unit_nr > 0) THEN
            DO P_loop_atom = 1, natom
               ikind = particle_set(P_loop_atom)%atomic_kind%kind_number
               cutoff_ri_per_kind(ikind) = MAX(cutoff_ri_per_kind(ikind), &
                                               cutoff_ri_per_atom(P_loop_atom))
            END DO
            WRITE (bs_env%unit_nr, '(T2,A)') 'Per-kind maximum RI-RS sphere cutoff (Bohr):'
            WRITE (bs_env%unit_nr, '(T4,A4,A14)') 'Kind', 'max cutoff_ri'
            DO ikind = 1, nkind
               WRITE (bs_env%unit_nr, '(T4,A4,F14.4)') &
                  atomic_kind_set(ikind)%element_symbol, &
                  cutoff_ri_per_kind(ikind)
            END DO
            WRITE (bs_env%unit_nr, '(A)') ' '
            DEALLOCATE (cutoff_ri_per_kind)
         END IF

         ! Shared 3c-integral context: hoists libint / t_c_g0 / md_ftable / contracted
         ! sphi tables out of the per-triple call so compute_d_lp threads only allocate
         ! a lightweight per-thread workspace. MPI-collective; must be outside any
         ! OMP region.
         CALL gw_3c_ctx_create(ctx_3c, bs_env, bs_env%ri_metric, &
                               basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO, &
                               basis_i=bs_env%basis_set_RI)

         ! =========================================================================
         ! 2. MPI LOOP OVER ATOMS (Fully independent, no MPI barriers inside)
         ! phi_local for each atom_P's cutoff sphere is built on the fly via
         ! image-summed AO evaluation; no dense replicated phi_global or allreduce.
         ! =========================================================================
         DO atom_P = para_env%mepos + 1, natom, para_env%num_pe

            n_loc_ri = ri_blk_sizes(atom_P)
            pos_P(:) = particle_set(atom_P)%r(:)

            cutoff_ri = cutoff_ri_per_atom(atom_P)
            cutoff_ri_2 = cutoff_ri**2

            ! ---------------------------------------------------------------------
            ! A. Determine Local Grid Domain based on cutoff_ri (PBC distance)
            ! ---------------------------------------------------------------------
            n_local_grid = 0
            DO l = 1, n_grid_total
               dist_vec_raw = pbc(ri_rs_grid_points(1:3, l), pos_P(1:3), cell)
               dist2_min = DOT_PRODUCT(dist_vec_raw, dist_vec_raw)
               IF (dist2_min <= cutoff_ri_2) n_local_grid = n_local_grid + 1
            END DO

            ALLOCATE (local_grid_idx(n_local_grid))

            n_local_grid = 0
            DO l = 1, n_grid_total
               dist_vec_raw = pbc(ri_rs_grid_points(1:3, l), pos_P(1:3), cell)
               dist2_min = DOT_PRODUCT(dist_vec_raw, dist_vec_raw)
               IF (dist2_min <= cutoff_ri_2) THEN
                  n_local_grid = n_local_grid + 1
                  local_grid_idx(n_local_grid) = l
               END IF
            END DO

            ! ---------------------------------------------------------------------
            ! B. Build phi_local on the fly with image-summed AO values.
            !    Only source atoms whose AO basis can reach the cutoff sphere of
            !    atom_P (MIC distance) contribute; the rest are pruned. The
            !    evaluate_ao_basis_on_points sums over (ix, iy, iz) images of
            !    ri_atom internally.
            ! ---------------------------------------------------------------------
            ALLOCATE (sphere_grid(3, n_local_grid))
            DO loc_idx = 1, n_local_grid
               sphere_grid(:, loc_idx) = ri_rs_grid_points(:, local_grid_idx(loc_idx))
            END DO

            ALLOCATE (phi_local(n_local_grid, n_ao_total))
            phi_local = 0.0_dp

            DO ri_atom = 1, natom
               dist_vec_raw = pbc(particle_set(ri_atom)%r(:), pos_P(:), cell)
               d_sP = NORM2(dist_vec_raw)
               IF (d_sP > bs_env%ri_rs%radius_ao_per_atom(ri_atom) + cutoff_ri) CYCLE

               col_start = bs_env%i_ao_start_from_atom(ri_atom)
               col_end = bs_env%i_ao_end_from_atom(ri_atom)
               r2_threshold = bs_env%ri_rs%radius_ao_per_atom(ri_atom)**2

               CALL evaluate_ao_basis_on_points(phi_local(:, col_start:col_end), sphere_grid, &
                                                bs_env%basis_set_AO( &
                                                particle_set(ri_atom)%atomic_kind%kind_number)%gto_basis_set, &
                                                particle_set(ri_atom)%r, cell, &
                                                cutoff_squared=r2_threshold, periodic_images=.TRUE.)
            END DO

            DEALLOCATE (sphere_grid)

            ! ---------------------------------------------------------------------
            ! C. Build Local RHS Matrix (d_lp_local). compute_d_lp does not depend
            !    on D_local or d_vec_local.
            ! ---------------------------------------------------------------------
            ALLOCATE (d_lp_local(n_local_grid, n_loc_ri))
            d_lp_local = 0.0_dp

            t2 = m_walltime()

            CALL compute_d_lp(bs_env, ctx_3c, phi_local, d_lp_local, n_local_grid, &
                              n_loc_ri, atom_P, max_ao_size)

            t3 = m_walltime()

            ! ---------------------------------------------------------------------
            ! D. Build d_vec_local (Jacobi diagonal) + LHS
            ! ---------------------------------------------------------------------
            ALLOCATE (d_vec_local(n_local_grid))

            CALL build_gram_jacobi_blas(phi_local, n_local_grid, n_ao_total, &
                                        bs_env%ri_rs%tikhonov, D_local, d_vec_local)

            ! ---------------------------------------------------------------------
            ! E. Pre-scale d_lp by d_vec
            ! ---------------------------------------------------------------------
            CALL scale_rows_by_diag(d_lp_local, d_vec_local, n_local_grid, n_loc_ri)

            ! ---------------------------------------------------------------------
            ! F. Solve — LAPACK dpotrf/dpotrs
            ! ---------------------------------------------------------------------
            CALL timeset(routineN//"_dpotrf", handle_dpotrf)
            CALL dpotrf('L', n_local_grid, D_local, n_local_grid, info)
            CALL timestop(handle_dpotrf)
            IF (info /= 0) CPABORT("RI-RS Cholesky factorization failed")
            CALL timeset(routineN//"_dpotrs", handle_dpotrs)
            CALL dpotrs('L', n_local_grid, n_loc_ri, D_local, n_local_grid, &
                        d_lp_local, n_local_grid, info)
            CALL timestop(handle_dpotrs)
            IF (info /= 0) CPABORT("RI-RS Cholesky solve failed")
            DEALLOCATE (D_local)

            ! ---------------------------------------------------------------------
            ! G. Post-scale solution by d_vec
            ! ---------------------------------------------------------------------
            CALL scale_rows_by_diag(d_lp_local, d_vec_local, n_local_grid, n_loc_ri)

            ! ---------------------------------------------------------------------
            ! H. Scatter Local Solution Back to Global DBCSR Matrix.
            !    DBCSR routes blocks to their global owner on finalize.
            !    local_grid_idx is ascending (built by the ordered scan above), so
            !    a single walking pointer over chunks works.
            ! ---------------------------------------------------------------------
            CALL store_Z_lP_columns(mat_Z_lP, d_lp_local, local_grid_idx, n_local_grid, &
                                    n_loc_ri, atom_P, r_blk_sizes, row_offset, &
                                    bs_env%eps_filter)

            DEALLOCATE (d_vec_local, d_lp_local)
            DEALLOCATE (local_grid_idx, phi_local)

         END DO

         DEALLOCATE (cutoff_ri_per_atom)
         CALL gw_3c_ctx_release(ctx_3c)

         CALL dbcsr_finalize(mat_Z_lP)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Computed Z_lP ', ' Execution time', m_walltime() - t1, ' s'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF

         logger => cp_get_default_logger()

         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, key), cp_p_file)) THEN
            CALL dbcsr_binary_write(matrix=mat_Z_lP, filepath=TRIM(bs_env%prefix)//"Z_lP.matrix")
         END IF

      END IF

      DEALLOCATE (row_offset, ri_blk_sizes, col_dist_ri)
      CALL dbcsr_distribution_release(dist_Z)

      DEALLOCATE (ri_rs_grid_points)

      CALL timestop(handle)

   END SUBROUTINE compute_Z_lP

! **************************************************************************************************
!> \brief Computes the dense localized RHS d_lp(l,P) = Σ_{μν,R,S} ϕ_μ(r_l)·ϕ_ν(r_l)·(μν|P) for one
!>        RI atom P. OMP-threaded over (atom_j, atom_k) AO-pair blocks: per thread, sweep all
!>        (cell_R, cell_S) periodic images of (atom_j, atom_k) about atom_P at cell (0,0,0); each
!>        3c block is built by build_3c_integral_block_ctx (cached libint / sphi tables in ctx,
!>        kind-radius triangle screen → `screened` short-circuits negligible image triples), and
!>        grid-chunked pair densities are contracted into a private d_lp partial that is reduced
!>        into d_lp at the end of the parallel region.
!> \param bs_env ...
!> \param ctx shared 3c-integral context (gw_3c_ctx_create)
!> \param phi_val ϕ_μ(r_l) on the local-sphere grid (n_grid_total × n_ao)
!> \param d_lp output (n_grid_total × n_loc_ri), zeroed by the caller, accumulated here
!> \param n_grid_total number of local-sphere grid rows
!> \param n_loc_ri number of RI functions of atom_P
!> \param atom_P RI atom (pinned to cell (0,0,0))
!> \param max_ao_size ...
! **************************************************************************************************

   SUBROUTINE compute_d_lp(bs_env, ctx, phi_val, d_lp, n_grid_total, n_loc_ri, atom_P, &
                           max_ao_size)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_val
      INTEGER, INTENT(IN)                                :: n_grid_total, n_loc_ri
      REAL(KIND=dp), INTENT(INOUT)                       :: d_lp(n_grid_total, n_loc_ri)
      INTEGER, INTENT(IN)                                :: atom_P, max_ao_size

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_d_lp'
      INTEGER, PARAMETER                                 :: grid_chunk = 1024

      INTEGER :: atom_j, atom_k, c, handle, handle_dgemm, ix_max, ix_min, ix_R, ix_S, iy_max, &
         iy_min, iy_R, iy_S, iz_max, iz_min, iz_R, iz_S, j, jk_idx, jsize, jstart, k, ksize, &
         kstart, l, l0, natom, ri
      INTEGER, DIMENSION(3)                              :: cell_R_vec, cell_S_vec
      LOGICAL                                            :: any_kept, screened
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_lp_prv, int_2d_prv, rho_chunk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c_prv, int_3c_sum
      TYPE(gw_3c_ws_type)                                :: ws

      CALL timeset(routineN, handle)

      natom = bs_env%n_atom

      IF (bs_env%ri_rs%cell%perd(1) == 1) THEN; ix_min = -1; ix_max = 1; ELSE; ix_min = 0; ix_max = 0
      END IF
      IF (bs_env%ri_rs%cell%perd(2) == 1) THEN; iy_min = -1; iy_max = 1; ELSE; iy_min = 0; iy_max = 0
      END IF
      IF (bs_env%ri_rs%cell%perd(3) == 1) THEN; iz_min = -1; iz_max = 1; ELSE; iz_min = 0; iz_max = 0
      END IF

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(bs_env, ctx, phi_val, d_lp, n_grid_total, n_loc_ri, atom_P, max_ao_size, &
      !$OMP        natom, ix_min, ix_max, iy_min, iy_max, iz_min, iz_max) &
      !$OMP PRIVATE(any_kept, atom_j, atom_k, c, handle_dgemm, j, jk_idx, jsize, jstart, k, &
      !$OMP         ksize, kstart, l, l0, ri, ix_R, iy_R, iz_R, ix_S, iy_S, iz_S, cell_R_vec, &
      !$OMP         cell_S_vec, screened, d_lp_prv, int_2d_prv, rho_chunk, int_3c_prv, int_3c_sum, ws)

      CALL gw_3c_ws_create(ws, ctx)
      ALLOCATE (int_3c_prv(max_ao_size, max_ao_size, n_loc_ri))
      ALLOCATE (int_3c_sum(max_ao_size, max_ao_size, n_loc_ri))
      ALLOCATE (int_2d_prv(max_ao_size*max_ao_size, n_loc_ri))
      ALLOCATE (rho_chunk(grid_chunk, max_ao_size*max_ao_size))
      ALLOCATE (d_lp_prv(n_grid_total, n_loc_ri))
      d_lp_prv(:, :) = 0.0_dp

      ! atom_P pinned at cell (0,0,0); enumerate (atom_j, cell_R) × (atom_k, cell_S). The ctx
      ! integral builder's kind_radius triangle screen sets screened=.TRUE. for the bulk of
      ! image triples (one or both AO atoms beyond the truncated-Coulomb reach of atom_P),
      ! so the 27 × 27 = 729 candidate cells collapse to "adjacent cells" in practice.
      ! The inner atom_k loop carries enough work for DYNAMIC.
      !$OMP DO SCHEDULE(DYNAMIC)
      DO atom_j = 1, natom
         DO atom_k = 1, natom
            jstart = bs_env%i_ao_start_from_atom(atom_j)
            jsize = bs_env%i_ao_end_from_atom(atom_j) - jstart + 1
            kstart = bs_env%i_ao_start_from_atom(atom_k)
            ksize = bs_env%i_ao_end_from_atom(atom_k) - kstart + 1

            int_3c_sum(1:jsize, 1:ksize, 1:n_loc_ri) = 0.0_dp
            any_kept = .FALSE.

            DO ix_R = ix_min, ix_max
            DO iy_R = iy_min, iy_max
            DO iz_R = iz_min, iz_max
               cell_R_vec = [ix_R, iy_R, iz_R]
               DO ix_S = ix_min, ix_max
               DO iy_S = iy_min, iy_max
               DO iz_S = iz_min, iz_max
                  cell_S_vec = [ix_S, iy_S, iz_S]

                  int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri) = 0.0_dp

                  CALL build_3c_integral_block_ctx(int_3c_prv( &
                                                   1:jsize, 1:ksize, 1:n_loc_ri), ctx, ws, &
                                                   atom_j=atom_j, atom_k=atom_k, atom_i=atom_P, &
                                                   cell_j=cell_R_vec, cell_k=cell_S_vec, cell_i=[0, 0, 0], &
                                                   screened=screened)
                  IF (screened) CYCLE

                  any_kept = .TRUE.
                  int_3c_sum(1:jsize, 1:ksize, 1:n_loc_ri) = &
                     int_3c_sum(1:jsize, 1:ksize, 1:n_loc_ri) + &
                     int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri)
               END DO
               END DO
               END DO
            END DO
            END DO
            END DO

            IF (.NOT. any_kept) CYCLE

            ! Flatten 3D B_{μν,P} → 2D B_{(μν),P}
            DO ri = 1, n_loc_ri
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     int_2d_prv(jk_idx, ri) = int_3c_sum(j, k, ri)
                  END DO
               END DO
            END DO

            ! Pair density ρ(l,μν) = ϕ_μ(r_l)ϕ_ν(r_l) in grid chunks, contracted on the fly:
            ! d_{l,P} += ρ(l,μν) B_{(μν),P}  (dgemm runs serially inside the parallel region)
            DO l0 = 1, n_grid_total, grid_chunk
               c = MIN(grid_chunk, n_grid_total - l0 + 1)
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     DO l = 1, c
                        rho_chunk(l, jk_idx) = phi_val(l0 + l - 1, jstart + j - 1)* &
                                               phi_val(l0 + l - 1, kstart + k - 1)
                     END DO
                  END DO
               END DO
               CALL timeset(routineN//"_dgemm", handle_dgemm)
               CALL dgemm("N", "N", c, n_loc_ri, jsize*ksize, &
                          1.0_dp, rho_chunk, grid_chunk, &
                          int_2d_prv, max_ao_size*max_ao_size, &
                          1.0_dp, d_lp_prv(l0, 1), n_grid_total)
               CALL timestop(handle_dgemm)
            END DO
         END DO
      END DO
      !$OMP END DO

      !$OMP CRITICAL (compute_d_lp_reduce)
      d_lp(1:n_grid_total, 1:n_loc_ri) = d_lp(1:n_grid_total, 1:n_loc_ri) + &
                                         d_lp_prv(1:n_grid_total, 1:n_loc_ri)
      !$OMP END CRITICAL (compute_d_lp_reduce)

      DEALLOCATE (int_3c_prv, int_3c_sum, int_2d_prv, rho_chunk, d_lp_prv)
      CALL gw_3c_ws_release(ws)

      !$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE compute_d_lp

END MODULE gw_ri_rs_large_cell_Gamma
