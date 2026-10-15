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

MODULE gw_ri_rs_non_periodic
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_add_on_diag,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_uplo_to_full
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE ri_rs_setup,           ONLY: ensure_ri_rs_representation
   USE gw_ri_rs_kernels,                ONLY: compute_Sigma_c,&
                                              compute_Sigma_x,&
                                              get_mat_chi_Gamma_tau,&
                                              print_matrix_occupation
   USE gw_tensor_large_cell_Gamma,      ONLY: Fourier_transform_w_to_t,&
                                              compute_fm_chi_Gamma_freq,&
                                              create_fm_W_MIC_time,&
                                              delete_unnecessary_files,&
                                              fm_write
   USE gw_utils,                        ONLY: analyt_conti_and_print,&
                                              de_init_bs_env
   USE gw_utils_fm,                     ONLY: fm_contract_ABA,&
                                              fm_invert,&
                                              fm_sqrt
   USE input_constants,                 ONLY: G0W0,&
                                              evGW0,&
                                              rtp_method_bse
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              max_line_length
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_print_mem_per_rank
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: get_all_VBM_CBM_bandgaps
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE ri_rs_matrices,                  ONLY: mask_grid_blocks_near_panel,&
                                              panel_mem_estimate_GB,&
                                              panel_template_elems,&
                                              resolve_grid_panels
   USE time_frequency_grids,            ONLY: negative_time,&
                                              positive_time,&
                                              transform_time_to_frequency
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_non_periodic'

   PUBLIC :: gw_calc_ri_rs_non_periodic

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism for molecules
!> \param qs_env ...
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************
   SUBROUTINE gw_calc_ri_rs_non_periodic(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_ri_rs_non_periodic'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time

      CALL timeset(routineN, handle)

      ! Build the common grid, phi and Z representation with the system's boundary conditions.
      CALL ensure_ri_rs_representation(qs_env, bs_env)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '
      CALL print_matrix_occupation(bs_env%ri_rs%mat_phi_mu_l, 'ϕ_μ(r_l)', bs_env)
      CALL print_ri_rs_memory_estimate(bs_env)

      CALL mp_print_mem_per_rank(bs_env%para_env, bs_env%unit_nr, &
                                 label='Memory per MPI process after computing Z_lP:')

      ! ==============================================================================
      ! 4. Polarizability matrix χ on the imaginary-time grid
      !      G^occ_µλ(i|τ|)  = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^vir_µλ(i|τ|)  = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^occ_ll'(i|τ|) = Σ_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')  (G^vir analogous)
      !      χ_ll'(iτ)       = G^occ_ll'(i|τ|) ∘ G^vir_ll'(i|τ|) (element-wise)
      !      χ_PQ(iτ)        = Σ_ll' Z_lP χ_ll'(iτ) Z_l'Q
      ! ==============================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP, bs_env%ri_rs%cutoff_radius_v_w)

      ! ==============================================================================
      ! 5. Screened Coulomb interaction W (RI basis)
      !      χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W_PQ(iτ)
      ! ==============================================================================
      CALL compute_W(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

      ! ==============================================================================
      ! 6. Exact-exchange self-energy Σ^x
      !      D_µν          = Σ_n^occ C_µn C_νn                (density matrix)
      !      D_ll'         = Σ_µν ϕ_µ(r_l) D_µν ϕ_ν(r_l')
      !      V^tr_ll'      = Σ_PQ Z_lP V^tr_PQ Z_l'Q          (truncated Coulomb)
      !      Σ^x_ll'       = D_ll' ∘ V^tr_ll'
      !      Σ^x_λσ(k=0)   = -Σ_ll' ϕ_λ(r_l) Σ^x_ll' ϕ_σ(r_l')
      ! ==============================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma, bs_env%ri_rs%cutoff_radius_v_w)

      ! ==============================================================================
      ! 7. Correlation self-energy Σ^c and quasiparticle energies, iterated
      !    until eigenvalue self-consistency if the &EVGW0 section is given and
      !    done in a single pass for G0W0.
      !
      !    (a) W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
      !        Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
      !        Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
      !        Σ^c_λσ(iτ)  =  Σ_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
      !    (b) Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
      !        ϵ_n^GW,(i) = ϵ_n^DFT + Σ^c_nn[G^(i-1),W](ϵ_n^GW,(i)) + Σ^x_nn - v^xc_nn
      ! ==============================================================================
      CALL compute_Sigma_c_and_QP_energies(bs_env, fm_W_time, fm_Sigma_x_Gamma)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_non_periodic

! **************************************************************************************************
!> \brief Correlation self-energy and quasiparticle energies, iterated to eigenvalue
!>        self-consistency in G (evGW0).
!>
!>        For G0W0 this runs once with G^(0) built from the DFT eigenvalues. For evGW0 the
!>        Green's function is rebuilt from the quasiparticle energies of all states until the
!>        quasiparticle HOMO, LUMO and HOMO-LUMO gap change by less than EPS_ITER between two
!>        cycles, or MAX_ITER cycles are spent. W, Σ^x and everything computed before this
!>        routine stay frozen: they either do not depend on the eigenvalues at all (Σ^x is G
!>        at τ = 0, a pure density matrix) or are held fixed by construction in GW0. The
!>        imaginary-time grid is fixed as well, since W(iτ) lives on it.
!>
!> \param bs_env ...
!> \param fm_W_time ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_and_QP_energies(bs_env, fm_W_time, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time, fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_c_and_QP_energies'

      INTEGER                                            :: handle, i_iter, n_iter
      LOGICAL                                            :: converged
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2)                        :: e_fermi_dft
      REAL(KIND=dp), DIMENSION(3, 2)                     :: band_prev
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

      converged = .FALSE.

      CALL init_evGW0(bs_env, n_iter, band_prev, eigenval_scf_Gamma_dft, e_fermi_dft)

      ! evGW0 self-consistency loop; for G0W0, loop is terminated after one iteration
      DO i_iter = 1, n_iter

         bs_env%evgw0_i_iter = i_iter

         ! W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         ! Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
         ! Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
         ! Σ^c_λσ(iτ)  =  Σ_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
         CALL compute_Sigma_c(bs_env, fm_W_time, bs_env%ri_rs%mat_phi_mu_l, &
                              bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time, bs_env%ri_rs%cutoff_radius_v_w)

         ! Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
         ! ϵ_n^GW,(i) = ϵ_n^DFT + Σ^c_nn[G^(i-1),W](ϵ_n^GW,(i)) + Σ^x_nn - v^xc_nn
         CALL compute_QP_energies(bs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

         IF (bs_env%gw_flavour == G0W0) EXIT

         CALL print_evGW0_band_edges(bs_env, band_prev, i_iter, n_iter, converged)

         IF (i_iter == n_iter .OR. converged) EXIT

         ! eigenvalues ϵ_n to be updated in G:
         ! G^occ_µλ(i|τ|) = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
         ! G^vir_µλ(i|τ|) = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
         CALL update_eigenvalues_G(bs_env)

      END DO

      CALL cp_fm_release(fm_W_time)

      CALL reset_and_clean_bs_env(bs_env, eigenval_scf_Gamma_dft, e_fermi_dft, fm_Sigma_x_Gamma)

      CALL delete_unnecessary_files(bs_env)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_and_QP_energies

! **************************************************************************************************
!> \brief Sets up the evGW0 eigenvalue self-consistency loop: the cycle count, the eigenvalues
!>        the first Green's function is built from, and the DFT reference that the loop
!>        overwrites. A G0W0 run reduces to a single cycle and needs none of it.
!> \param bs_env ...
!> \param n_iter ...
!> \param band_prev ...
!> \param eigenval_scf_Gamma_dft ...
!> \param e_fermi_dft ...
! **************************************************************************************************
   SUBROUTINE init_evGW0(bs_env, n_iter, band_prev, eigenval_scf_Gamma_dft, e_fermi_dft)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(OUT)                               :: n_iter
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(OUT)        :: band_prev
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: e_fermi_dft

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_evGW0'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      n_iter = 1
      band_prev(:, :) = 0.0_dp
      e_fermi_dft(:) = 0.0_dp

      IF (bs_env%gw_flavour == evGW0) THEN
         n_iter = bs_env%evgw0_iter

         ! eigenvalues currently in G; the first cycle starts is a plain G0W0 step
         bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_scf(:, :, :)
         ! the loop overwrites these; post-GW printing expects the DFT values back
         ALLOCATE (eigenval_scf_Gamma_dft, SOURCE=bs_env%eigenval_scf_Gamma)
         e_fermi_dft(:) = bs_env%e_fermi(:)
      END IF

      CALL timestop(handle)

   END SUBROUTINE init_evGW0

! **************************************************************************************************
!> \brief Feeds the quasiparticle energies of the current evGW0 cycle back into the Green's
!>        function used by the next one, and re-centres the Fermi level between the new
!>        quasiparticle HOMO and LUMO.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE update_eigenvalues_G(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'update_eigenvalues_G'

      INTEGER                                            :: handle, i_mo, ispin

      CALL timeset(routineN, handle)

      ! record this cycle's evGW0 result; it is also what the next Green's function is built from
      bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_GW(:, :, :)

      DO ispin = 1, bs_env%n_spin
         ! Update all physical states; exclude linear-dependency placeholders.
         DO i_mo = 1, bs_env%n_mo_retained
            bs_env%eigenval_scf_Gamma(i_mo, ispin) = bs_env%eigenval_GW(i_mo, 1, ispin)
         END DO
         bs_env%e_fermi(ispin) = &
            0.5_dp*(bs_env%eigenval_GW(bs_env%n_occ(ispin), 1, ispin) + &
                    bs_env%eigenval_GW(bs_env%n_occ(ispin) + 1, 1, ispin))
      END DO

      CALL timestop(handle)

   END SUBROUTINE update_eigenvalues_G

! **************************************************************************************************
!> \brief Restores the DFT reference that the evGW0 loop overwrote and cleanup
!> \param bs_env ...
!> \param eigenval_scf_Gamma_dft ...
!> \param e_fermi_dft ...
!> \param fm_Sigma_x_Gamma ..
! **************************************************************************************************
   SUBROUTINE reset_and_clean_bs_env(bs_env, eigenval_scf_Gamma_dft, e_fermi_dft, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: e_fermi_dft
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reset_and_clean_bs_env'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (bs_env%gw_flavour == evGW0) THEN
         ! Update array with final evGW0 result
         bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_GW(:, :, :)
         ! put the DFT reference back for the post-GW DOS/band-edge printing
         bs_env%eigenval_scf_Gamma(:, :) = eigenval_scf_Gamma_dft(:, :)
         bs_env%e_fermi(:) = e_fermi_dft(:)
         DEALLOCATE (eigenval_scf_Gamma_dft)
         CALL cp_fm_release(fm_Sigma_x_Gamma)
      END IF

      CALL timestop(handle)

   END SUBROUTINE reset_and_clean_bs_env

! **************************************************************************************************
!> \brief Print the quasiparticle HOMO, LUMO and HOMO-LUMO gap of the current evGW0 cycle
!> \param bs_env ...
!> \param band_prev ...
!> \param i_iter ...
!> \param n_iter ...
!> \param converged ...
! **************************************************************************************************
   SUBROUTINE print_evGW0_band_edges(bs_env, band_prev, i_iter, n_iter, converged)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(INOUT)      :: band_prev
      INTEGER, INTENT(IN)                                :: i_iter, n_iter
      LOGICAL, INTENT(OUT)                               :: converged

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_evGW0_band_edges'

      INTEGER                                            :: handle, homo, ispin, u
      REAL(KIND=dp)                                      :: max_delta
      REAL(KIND=dp), DIMENSION(3)                        :: band

      CALL timeset(routineN, handle)

      u = bs_env%unit_nr

      converged = (i_iter > 1)
      max_delta = 0.0_dp

      IF (u > 0) THEN
         WRITE (u, '(A)') ' '
         WRITE (u, '(T2,A)') REPEAT('-', 79)
         WRITE (u, '(T2,A,I4,A,I4)') 'evGW0 cycle', i_iter, ' /', n_iter
         WRITE (u, '(T2,A)') REPEAT('-', 79)
      END IF

      DO ispin = 1, bs_env%n_spin

         homo = bs_env%n_occ(ispin)
         band(1) = bs_env%eigenval_GW(homo, 1, ispin)
         band(2) = bs_env%eigenval_GW(homo + 1, 1, ispin)
         band(3) = band(2) - band(1)

         IF (i_iter > 1) THEN
            max_delta = MAX(max_delta, MAXVAL(ABS(band(:) - band_prev(:, ispin))))
         END IF

         IF (u > 0) THEN
            IF (bs_env%n_spin == 2) WRITE (u, '(T2,A,I0)') 'Spin ', ispin
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 HOMO (eV)', band(1)*evolt
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 LUMO (eV)', band(2)*evolt
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 HOMO-LUMO gap (eV)', band(3)*evolt
         END IF

         band_prev(:, ispin) = band(:)

      END DO

      IF (i_iter > 1) THEN
         converged = (max_delta < bs_env%evgw0_eps_iter)
         IF (u > 0) WRITE (u, '(T2,A,T61,F20.6)') 'Max. change to previous cycle (eV)', &
            max_delta*evolt
      END IF

      IF (u > 0) WRITE (u, '(T2,A)') REPEAT('-', 79)

      IF (converged) THEN
         IF (u > 0) THEN
            WRITE (u, '(A)') ' '
            WRITE (u, '(T2,A,I4,A)') &
               'evGW0 eigenvalue self-consistency reached in', i_iter, ' cycles.'
            WRITE (u, '(A)') ' '
         END IF
      ELSE IF (i_iter == n_iter) THEN
         CALL cp_warn(__LOCATION__, &
                      "The evGW0 eigenvalue self-consistency cycle did not converge "// &
                      "within MAX_ITER cycles. The reported quasiparticle energies are "// &
                      "those of the last cycle.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_evGW0_band_edges

! **************************************************************************************************
!> \brief Estimates and prints per-process memory requirements for the RI-RS GW calculation.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE print_ri_rs_memory_estimate(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_ri_rs_memory_estimate'

      CHARACTER(LEN=max_line_length)                     :: msg
      INTEGER                                            :: handle, ipan, max_n_local_grid, &
                                                            n_grid_total, n_procs, n_RI, pan_rows, &
                                                            pan_width
      INTEGER(KIND=int_8)                                :: nze_tmpl
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: r_blk_sizes
      LOGICAL                                            :: use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used
      REAL(KIND=dp)                                      :: mem_avail_GB, mem_gram_GB, mem_pan_GB, &
                                                            mem_panels_GB, mem_Z_lP_GB, &
                                                            mem_Zlp_peak_GB

      CALL timeset(routineN, handle)

      n_grid_total = bs_env%ri_rs%n_grid_points
      CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, row_blk_size=r_blk_sizes)
      CPASSERT(SUM(r_blk_sizes) == n_grid_total)
      n_RI = bs_env%n_RI
      n_procs = bs_env%para_env%num_pe

      ! Z_lP upper bound: dense n_grid × n_RI, distributed evenly across all ranks.
      ! The actual sparse Z_lP is smaller due to the per-atom locality cutoff.
      mem_Z_lP_GB = REAL(n_grid_total, dp)*REAL(n_RI, dp)*8.0_dp/ &
                    REAL(n_procs, dp)*1.0E-9_dp

      ! Peak panel memory during Σ^c: two G panels (A_occ, A_vir) + one W panel plus the
      ! grid × RI / grid × AO intermediates. With the CUTOFF_RADIUS_RL_W restriction the panel
      ! matrices only allocate the geo-template blocks, so use the same nze-aware model as the
      ! panel planner (panel_mem_estimate_GB); without the cutoff, dense panel_rows × n_grid.
      ! Plus the n_RI × n_RI W_aux matrix. Distributed over n_procs ranks.
      use_cutoff = bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp .AND. &
                   ALLOCATED(bs_env%ri_rs%chunk_centroids)
      CALL resolve_grid_panels(bs_env%ri_rs, bs_env%para_env, bs_env%dbcsr_msg_elem_limit, &
                               pan_first, pan_last, bs_env%ri_rs%cutoff_radius_v_w)
      mem_panels_GB = 0.0_dp
      DO ipan = 1, SIZE(pan_first)
         pan_rows = SUM(r_blk_sizes(pan_first(ipan):pan_last(ipan)))
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(bs_env%ri_rs%chunk_centroids, pan_first(ipan), &
                                             pan_last(ipan), bs_env%ri_rs%cutoff_radius_v_w, &
                                             grid_used)
            pan_width = SUM(r_blk_sizes, MASK=grid_used)
            CALL panel_template_elems(r_blk_sizes, bs_env%ri_rs%chunk_centroids, &
                                      grid_used, pan_first(ipan), pan_last(ipan), &
                                      bs_env%ri_rs%cutoff_radius_v_w, nze_tmpl)
            CALL panel_mem_estimate_GB(nze_tmpl, pan_rows, pan_width, n_grid_total, &
                                       n_RI, bs_env%n_ao, n_procs, mem_pan_GB)
         ELSE
            pan_width = n_grid_total
            mem_pan_GB = (3.0_dp*REAL(pan_rows, dp)*REAL(pan_width, dp) + &
                          2.0_dp*REAL(pan_rows, dp)*REAL(n_RI, dp))* &
                         8.0_dp/REAL(n_procs, dp)*1.0E-9_dp
         END IF
         mem_panels_GB = MAX(mem_panels_GB, mem_pan_GB)
      END DO
      mem_panels_GB = mem_panels_GB + &
                      REAL(n_RI, dp)*REAL(n_RI, dp)*8.0_dp/REAL(n_procs, dp)*1.0E-9_dp

      ! Fitting setup records the same Gram + phi + threaded RHS model used for
      ! solver selection. AUTO-RI records its own common-grid or localized dimensions.
      ! A restart read has no fitting peak and needs no diagnostic geometry scan.
      max_n_local_grid = bs_env%ri_rs%fit_max_grid
      mem_gram_GB = 8.0_dp*REAL(max_n_local_grid, dp)**2*1.0E-9_dp
      mem_Zlp_peak_GB = bs_env%ri_rs%fit_peak_bytes*1.0E-9_dp

      ! Available memory per process = node MemLikelyFree / ranks-per-node, min across ranks
      ! (0 on non-Linux => warnings suppressed below).
      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, mem_avail_GB)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ' '
         WRITE (bs_env%unit_nr, '(T2,A)') 'RI-RS memory estimate per MPI process:'
         WRITE (bs_env%unit_nr, '(T4,A,F37.2,A)') &
            'Available memory per process (system)', mem_avail_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F18.2,A)') &
            'Required for Z_lP (dense upper bound; actual is sparser)', mem_Z_lP_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F25.2,A)') &
            'Required for χ, W, Σ panels (peak per panel step)', mem_panels_GB, ' GB'
         IF (bs_env%ri_rs%Z_lP_exists .AND. .NOT. bs_env%auto_ri%enabled) THEN
            WRITE (bs_env%unit_nr, '(T4,A)') 'Z_lP fitting skipped (matrix read from restart)'
         ELSE
            WRITE (bs_env%unit_nr, '(T4,A,F17.2,A)') &
               'Required for Z_lP fit (Gram+phi+RHS, single-rank bound)', mem_Zlp_peak_GB, ' GB'
            WRITE (bs_env%unit_nr, '(T4,A,T69,I12)') &
               'Largest fitting-grid number of grid points:', max_n_local_grid
            WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
               'Largest dense Gram matrix:', mem_gram_GB, ' GB'
         END IF

      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_Z_lP_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated memory for Z_lP, ", mem_Z_lP_GB, " GB per process, exceeds the "// &
            "available ", mem_avail_GB, " GB. Z_lP (n_grid x n_RI) is distributed across all "// &
            "MPI ranks, so add nodes or use fewer MPI ranks per node."
         CPWARN(TRIM(msg))
      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_panels_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated peak memory of the chi/W/Sigma panels, ", mem_panels_GB, &
            " GB per process, exceeds the available ", mem_avail_GB, " GB. Panel memory "// &
            "scales roughly as 3*panel_size*n_grid/n_procs, so add nodes, use fewer MPI ranks "// &
            "per node."
         CPWARN(TRIM(msg))
      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_Zlp_peak_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated peak memory of the Z_lP solve, ", mem_Zlp_peak_GB, &
            " GB per process, exceeds the available ", mem_avail_GB, &
            " GB. The dense Gram matrix scales as n_local_grid^2; "// &
            "the estimate assumes a single-rank solve. Localized fits may distribute the Gram "// &
            "matrix over a ScaLAPACK subgroup. Use fewer MPI ranks per node if more memory "// &
            "per rank is needed."
         CPWARN(TRIM(msg))
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_ri_rs_memory_estimate

! **************************************************************************************************
!> \brief Computes the screened Coulomb interaction on the imaginary-time grid, entirely in the
!>        RI auxiliary (PQ) basis:
!>          χ_PQ(iω) = Σ_τ w(ω,τ) cos(ωτ) χ_PQ(iτ)                  (cosine transform)
!>          ε(iω)    = Id - V^0.5 M^-1 χ(iω) M^-1 V^0.5             (dielectric function)
!>          W(iω)    = V^0.5 ( ε^-1(iω) - Id ) V^0.5                (correlation part only)
!>          W(iτ)    = Σ_ω w̃(τ,ω) cos(ωτ) W(iω)                    (back transform)
!>          W(iτ)   <- M^-1 W(iτ) M^-1                              (fold in the RI metric)
!>        where V is the bare Coulomb matrix and M the RI metric.
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_time ...
! **************************************************************************************************
   SUBROUTINE compute_W(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W'

      INTEGER                                            :: handle, i_t, j_w
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type)                                   :: fm_M_inv_V_sqrt, fm_V, fm_V_sqrt

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL create_fm_W_MIC_time(bs_env, fm_W_time)

      ! 1. Allocate V and M matrices
      CALL cp_fm_create(fm_V, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_M_inv_V_sqrt, bs_env%fm_RI_RI%matrix_struct)

      ! Compute V and M^-1 * V^0.5
      CALL compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_M_inv_V_sqrt)

      ! 2. Loop over frequencies
      DO j_w = 1, bs_env%num_time_freq_points
         ! Fourier transformation of χ_PQ(iτ) to χ_PQ(iω_j)
         CALL compute_fm_chi_Gamma_freq(bs_env, bs_env%fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)

         ! ε(iω_j) = Id - V^0.5*M^-1*χ(iω_j)*M^-1*V^0.5
         ! W(iω_j) = V^0.5*(ε^-1(iω_j)-Id)*V^0.5
         CALL compute_fm_W_freq(bs_env, bs_env%fm_chi_Gamma_freq, fm_V_sqrt, &
                                fm_M_inv_V_sqrt, bs_env%fm_W_MIC_freq)

         ! Fourier transform from W_PQ^MIC(iω_j) to W_PQ^MIC(iτ)
         CALL Fourier_transform_w_to_t(bs_env, fm_W_time, bs_env%fm_W_MIC_freq, j_w)
      END DO

      ! M^-1(k=0) W^MIC(iτ) M^-1(k=0) -> fm_W_time
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_W_time)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed W(iτ),', ' Execution time', m_walltime() - t1, ' s'
      END IF

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      ! Cleanup
      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_V_sqrt)
      CALL cp_fm_release(fm_M_inv_V_sqrt)

      ! Marek : Fourier transform W^MIC(itau) back to get it at a specific im.frequency point - iomega = 0
      IF (bs_env%rtp_method == rtp_method_bse) THEN
         t1 = m_walltime()
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_W_MIC_freq%matrix_struct)
         ! Set to zero
         CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_zero, 0.0_dp)
         ! Sum over all times
         DO i_t = 1, bs_env%num_time_freq_points
            ! Add the relevant structure with correct weight
            CALL cp_fm_scale_and_add(1.0_dp, bs_env%fm_W_MIC_freq_zero, &
                                     bs_env%time_frequency_grid%time_weights_at_zero_frequency(i_t), fm_W_time(i_t))
         END DO
         ! Done, save to file
         CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
         ! Report calculation
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Computed W(0),', ' Execution time', m_walltime() - t1, ' s'
         END IF
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W

! **************************************************************************************************
!> \brief Computes V, V^0.5, and M^-1 V^0.5 for the RI-RS dielectric function.
!>        The Coulomb matrix V is constructed by the RI-RS k-point path. The inverse metric
!>        M^-1(k=0) is precomputed once in gw_utils and stored in bs_env.
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V Coulomb matrix V(k=0)
!> \param fm_V_sqrt symmetric factor V^0.5
!> \param fm_Minv_Vsqrt product M^-1(k=0) V^0.5
! **************************************************************************************************
   SUBROUTINE compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_Minv_Vsqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_V, fm_V_sqrt, fm_Minv_Vsqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_MinvVsqrt'

      INTEGER                                            :: handle
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (bs_env%auto_ri%enabled) THEN
         ! -------------------------------------------------------------------
         ! 1a. The optimized AB functions span two atoms, so their previously transformed
         !     Coulomb matrix cannot be rebuilt by the atom-local k-point integral routine.
         ! -------------------------------------------------------------------
         CALL cp_fm_to_fm(bs_env%auto_ri%V_pq, fm_V)
      ELSE
         ! -------------------------------------------------------------------
         ! 1b. Build Coulomb Matrix V(k=0) using the kp-routine but only for ikp=1.
         ! -------------------------------------------------------------------
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                         qs_kind_set=qs_kind_set)
         atomic_kind_set => bs_env%ri_rs%atomic_kind_set

         ALLOCATE (mat_V_kp(1:1, 1:2))
         NULLIFY (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
         ALLOCATE (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
         CALL dbcsr_create(mat_V_kp(1, 1)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 1)%matrix)
         CALL dbcsr_set(mat_V_kp(1, 1)%matrix, 0.0_dp)
         CALL dbcsr_create(mat_V_kp(1, 2)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 2)%matrix)
         ! The dummy imaginary part is required only by the k-point routine interface.
         CALL dbcsr_set(mat_V_kp(1, 2)%matrix, 0.0_dp)

         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig
         CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, "RI_AUX", cell, &
                                         particle_set, qs_kind_set, atomic_kind_set, &
                                         bs_env%size_lattice_sum_V, operator_coulomb, 1, 1)
         CALL copy_dbcsr_to_fm(mat_V_kp(1, 1)%matrix, fm_V)

         CALL dbcsr_deallocate_matrix(mat_V_kp(1, 1)%matrix)
         CALL dbcsr_deallocate_matrix(mat_V_kp(1, 2)%matrix)
         DEALLOCATE (mat_V_kp)
      END IF

      ! -----------------------------------------------------------------------
      ! 2. V -> V^0.5.
      ! -----------------------------------------------------------------------
      CALL fm_sqrt(fm_V, fm_V_sqrt, bs_env%eps_eigval_mat_RI, bs_env%unit_nr)

      ! -----------------------------------------------------------------------
      ! 3. M^-1(k=0) V^0.5.
      ! -----------------------------------------------------------------------
      CALL parallel_gemm("N", "T", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI, 1.0_dp, &
                         bs_env%fm_Minv_Gamma, fm_V_sqrt, 0.0_dp, fm_Minv_Vsqrt)

      CALL timestop(handle)

   END SUBROUTINE compute_V_MinvVsqrt

! **************************************************************************************************
!> \brief Computes the screened interaction at one imaginary frequency:
!>          ε(iω_j)  = Id - (M^-1 V^0.5)^T χ(iω_j) (M^-1 V^0.5)
!>          W(iω_j)  = V^0.5^T ( ε^-1(iω_j) - Id ) V^0.5
!>        ε is inverted via Cholesky; if that fails due to conditioning, via
!>        eigendecomposition (cp_fm_power) with eigenvalue filtering.
!> \param bs_env ...
!> \param fm_chi_freq_j ...
!> \param fm_V_sqrt ...
!> \param fm_Minv_Vsqrt ...
!> \param fm_W_freq_j ...
! **************************************************************************************************
   SUBROUTINE compute_fm_W_freq(bs_env, fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt, fm_W_freq_j)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_W_freq_j

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_fm_W_freq'

      INTEGER                                            :: handle, n_RI
      TYPE(cp_fm_type)                                   :: fm_eps_freq_j, fm_work

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      CALL cp_fm_create(fm_eps_freq_j, fm_chi_freq_j%matrix_struct)
      CALL cp_fm_create(fm_work, fm_chi_freq_j%matrix_struct)

      ! -----------------------------------------------------------------------
      ! 1. ε(iω_j) = Id - (M^-1 * V^0.5)^T * χ(iω_j) * (M^-1 * V^0.5)
      ! -----------------------------------------------------------------------
      ! work = χ(iω_j) * (M^-1 * V^0.5)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_chi_freq_j, fm_Minv_Vsqrt, 0.0_dp, fm_work)

      ! eps_work = (M^-1 * V^0.5)^T * work
      CALL parallel_gemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_Minv_Vsqrt, fm_work, 0.0_dp, fm_eps_freq_j)

      ! ε(iω_j) = Id - eps_work  -->  -eps_work + Id
      CALL cp_fm_add_on_diag(fm_eps_freq_j, 1.0_dp)

      ! Force perfect symmetry before Cholesky to avoid info != 0 due to GEMM noise
      CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)

      ! -----------------------------------------------------------------------
      ! 2. W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      ! -----------------------------------------------------------------------

      ! a) Invert ε by Cholesky decomposition or, if that fails, by diagonalization.
      CALL fm_invert(fm_eps_freq_j, bs_env%eps_eigval_mat_RI, bs_env%unit_nr)

      ! b) ε^-1(iω_j) - Id
      CALL cp_fm_add_on_diag(fm_eps_freq_j, -1.0_dp)

      ! c) W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      CALL fm_contract_ABA(fm_V_sqrt, fm_eps_freq_j, fm_W_freq_j)

      ! Cleanup
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_eps_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_freq

! **************************************************************************************************
!> \brief Σ^c_λσ(iτ) -> Σ^c_nn(ϵ) and the quasi-particle levels of the non-periodic RI-RS path,
!>        ϵ_n^GW = ϵ_n^DFT + Σ^c_nn(ϵ_n^GW) + Σ^x_nn - v^xc_nn.
!> \param bs_env ...
!> \param fm_Sigma_x_Gamma ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(bs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: Sigma_c_n_freq
      INTEGER                                            :: handle, ispin, j_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_n, V_xc_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_n_time
      TYPE(cp_fm_type)                                   :: fm_ks, fm_mos, fm_s, fm_work

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_ks, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_s, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_mos, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_work, bs_env%fm_s_Gamma%matrix_struct)

      ALLOCATE (V_xc_n(bs_env%n_ao), Sigma_x_n(bs_env%n_ao))
      ALLOCATE (Sigma_c_n_time(bs_env%n_ao, bs_env%num_time_freq_points, 2))
      ALLOCATE (Sigma_c_n_freq(bs_env%n_ao, bs_env%num_freq_points_fit))

      DO ispin = 1, bs_env%n_spin

         ! 1. Roothaan-Hall H^KS_µν C_νn = S_µν C_νn ϵ_n
         CALL cp_fm_to_fm(bs_env%fm_ks_Gamma(ispin), fm_ks)
         CALL cp_fm_to_fm(bs_env%fm_s_Gamma, fm_s)
         CALL cp_fm_geeig(fm_ks, fm_s, fm_mos, bs_env%eigenval_scf(:, 1, ispin), fm_work)

         ! 2. v^xc_µν -> v^xc_nn  and  Σ^x_µν -> Σ^x_nn
         CALL to_Gamma_and_mo_real(V_xc_n, bs_env%fm_V_xc_Gamma(ispin), fm_mos)
         CALL to_Gamma_and_mo_real(Sigma_x_n, fm_Sigma_x_Gamma(ispin), fm_mos)

         ! 3. Σ^c_µν(+/-i|τ_j|) -> Σ^c_nn(+/-i|τ_j|)
         DO j_t = 1, bs_env%num_time_freq_points
            CALL to_Gamma_and_mo_real(Sigma_c_n_time(:, j_t, positive_time), &
                                      fm_Sigma_c_Gamma_time(j_t, positive_time, ispin), fm_mos)
            CALL to_Gamma_and_mo_real(Sigma_c_n_time(:, j_t, negative_time), &
                                      fm_Sigma_c_Gamma_time(j_t, negative_time, ispin), fm_mos)
         END DO

         ! 4. Σ^c_nn(iτ) -> Σ^c_nn(iω)
         CALL transform_time_to_frequency(bs_env%time_frequency_grid, Sigma_c_n_time, Sigma_c_n_freq)
         Sigma_c_n_freq(1:bs_env%n_occ(ispin), :) = CONJG(Sigma_c_n_freq(1:bs_env%n_occ(ispin), :))

         ! 5. Analytic continuation Σ^c_nn(iω) -> Σ^c_nn(ϵ) and the QP levels
         CALL analyt_conti_and_print(bs_env, Sigma_c_n_freq, Sigma_x_n, V_xc_n, &
                                     bs_env%eigenval_scf(:, 1, ispin), 1, ispin)

      END DO ! ispin

      CALL get_all_VBM_CBM_bandgaps(bs_env)

      IF (bs_env%gw_flavour == G0W0) CALL cp_fm_release(fm_Sigma_x_Gamma)
      CALL cp_fm_release(fm_Sigma_c_Gamma_time)

      CALL cp_fm_release(fm_ks)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_mos)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies

! **************************************************************************************************
!> \brief AO -> MO transform of a Γ-point matrix
!> \param array_n ...
!> \param fm_Gamma ...
!> \param fm_mos ...
! **************************************************************************************************
   SUBROUTINE to_Gamma_and_mo_real(array_n, fm_Gamma, fm_mos)

      REAL(KIND=dp), DIMENSION(:)                        :: array_n
      TYPE(cp_fm_type)                                   :: fm_Gamma, fm_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'to_Gamma_and_mo_real'

      INTEGER                                            :: handle
      TYPE(cp_fm_type)                                   :: fm_mo

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mo, fm_Gamma%matrix_struct)

      ! A_nn' = Σ_μν C_μn A_μν C_νn'
      CALL fm_contract_ABA(fm_mos, fm_Gamma, fm_mo)

      CALL cp_fm_get_diag(fm_mo, array_n)

      CALL cp_fm_release(fm_mo)

      CALL timestop(handle)

   END SUBROUTINE to_Gamma_and_mo_real

END MODULE gw_ri_rs_non_periodic
