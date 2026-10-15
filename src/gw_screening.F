!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Imaginary-time GW screening, minimum-image integration and screened-interaction restart.
! **************************************************************************************************
MODULE gw_screening
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_cfm_diag,                     ONLY: cp_cfm_power
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE gw_utils,                        ONLY: Fourier_transform_w_to_t,&
                                              compute_fm_chi_Gamma_freq,&
                                              create_fm_W_MIC_time,&
                                              fm_read,&
                                              fm_write,&
                                              sum_W_MIC_freq_zero
   USE gw_utils_fm,                     ONLY: cfm_invert,&
                                              fm_contract_ABA
   USE input_constants,                 ONLY: ri_rs_non_periodic,&
                                              rtp_method_bse,&
                                              rtp_method_bse_linearized
   USE kinds,                           ONLY: dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE machine,                         ONLY: m_walltime
   USE mbpt_screening,                  ONLY: prepare_screening_factors,&
                                              screened_interaction
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat_kp
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: MIC_contribution_from_ikp,&
                                              cfm_ikp_from_fm_Gamma
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_screening'
   PUBLIC :: compute_W, get_W_MIC, compute_V_k_by_lattice_sum, compute_MinvVsqrt_Vsqrt, &
             compute_fm_W_MIC_freq_j_from_chi, apply_extrapol_factor

CONTAINS

! **************************************************************************************************
!> \brief Compute the correlation interaction from the Gamma imaginary-time response.
!>        The molecular RI-RS channel uses real factors directly. Other Gamma paths retain
!>        their complex k-point factors, MIC projection and extrapolation. All channels share
!>        frequency integration, the final M^-1 sandwich and static-interaction construction.
!> \param bs_env screening settings, reusable matrices and owned static output
!> \param qs_env electronic-structure environment
!> \param mat_chi_Gamma_tau input response matrices, released after consumption
!> \param fm_W_time allocated output M^-1 Wc(i tau) M^-1; owned by the caller
! **************************************************************************************************
   SUBROUTINE compute_W(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), POINTER                        :: mat_chi_Gamma_tau(:)
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_W_time(:)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W'

      INTEGER                                            :: handle, i_t, ikp, ikp_batch, &
                                                            ikp_in_batch, j_w, n_batches, &
                                                            nkp_batch, nkp_total
      LOGICAL                                            :: keep_static, real_screening, write_static
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_cfm_type)                                  :: cfm_Minv_Vsqrt, cfm_V_sqrt
      TYPE(cp_cfm_type), ALLOCATABLE                     :: cfm_V_kp(:)
      TYPE(cp_fm_type)                                   :: fm_Minv_Vsqrt, fm_V, fm_V_sqrt

      CALL timeset(routineN, handle)
      CALL create_fm_W_MIC_time(bs_env, fm_W_time)

      real_screening = bs_env%gw_implementation == ri_rs_non_periodic
      n_batches = 1
      nkp_batch = 1
      nkp_total = 1
      IF (real_screening) THEN
         CALL cp_fm_create(fm_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
         CALL cp_fm_create(fm_Minv_Vsqrt, bs_env%fm_RI_RI%matrix_struct)
      ELSE
         n_batches = bs_env%num_chi_eps_W_batches
         nkp_batch = bs_env%nkp_chi_eps_W_batch
         nkp_total = bs_env%nkp_chi_eps_W_orig_plus_extra
      END IF

      DO ikp_batch = 1, n_batches
         t1 = m_walltime()
         IF (real_screening) THEN
            CALL cp_fm_create(fm_V, bs_env%fm_RI_RI%matrix_struct)
            CALL compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_Minv_Vsqrt)
            CALL cp_fm_release(fm_V)
         ELSE
            CALL compute_V_k_by_lattice_sum(bs_env, qs_env, cfm_V_kp, ikp_batch)
         END IF

         DO ikp_in_batch = 1, nkp_batch
            ikp = (ikp_batch - 1)*nkp_batch + ikp_in_batch
            IF (ikp > nkp_total) CYCLE
            IF (.NOT. real_screening) THEN
               CALL compute_MinvVsqrt_Vsqrt(bs_env, qs_env, cfm_V_kp, cfm_V_sqrt, cfm_Minv_Vsqrt, ikp)
               CALL bs_env%para_env%sync()
               CALL cp_cfm_release(cfm_V_kp(ikp))
            END IF

            DO j_w = 1, bs_env%num_time_freq_points
               IF (.NOT. real_screening .AND. bs_env%approx_kp_extrapol .AND. j_w > 1 .AND. &
                   ikp > bs_env%nkp_chi_eps_W_orig) CYCLE
               ! Transform the Gamma response once for this channel and frequency.
               CALL compute_fm_chi_Gamma_freq(bs_env, bs_env%fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)
               IF (real_screening) THEN
                  CALL screened_interaction(bs_env%fm_chi_Gamma_freq, fm_V_sqrt, fm_Minv_Vsqrt, &
                                            bs_env%eps_eigval_mat_RI, matrix_Wc=bs_env%fm_W_MIC_freq, &
                                            unit_nr=bs_env%unit_nr)
               ELSE
                  CALL compute_fm_W_MIC_freq_j_from_chi(bs_env, qs_env, bs_env%fm_W_MIC_freq, j_w, ikp, &
                                                        cfm_Minv_Vsqrt, cfm_V_sqrt)
               END IF
               ! Accumulate the real Gamma/MIC contribution into every imaginary-time sample.
               CALL Fourier_transform_w_to_t(bs_env, fm_W_time, bs_env%fm_W_MIC_freq, j_w)
            END DO
         END DO

         IF (.NOT. real_screening) DEALLOCATE (cfm_V_kp)
         IF (bs_env%unit_nr > 0) THEN
            IF (real_screening) THEN
               WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
                  'Computed W(iτ),', ' Execution time', m_walltime() - t1, ' s'
            ELSE
               WRITE (bs_env%unit_nr, '(T2,A,I12,A,I3,A,F10.1,A)') &
                  'Computed W(iτ,k) for k-point batch', ikp_batch, ' /', n_batches, &
                  ', Execution time', m_walltime() - t1, ' s'
            END IF
         END IF
      END DO

      IF (real_screening) THEN
         CALL cp_fm_release(fm_V_sqrt)
         CALL cp_fm_release(fm_Minv_Vsqrt)
      ELSE
         CALL cp_cfm_release(cfm_V_sqrt)
         CALL cp_cfm_release(cfm_Minv_Vsqrt)
         IF (bs_env%approx_kp_extrapol) CALL apply_extrapol_factor(bs_env, fm_W_time)
      END IF

      ! Fold the Gamma metric into the integrated interaction exactly once.
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_W_time)
      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)
      IF (.NOT. real_screening) THEN
         DO i_t = 1, bs_env%num_time_freq_points
            CALL fm_write(fm_W_time(i_t), i_t, bs_env%W_time_name, qs_env)
         END DO
      END IF

      write_static = bs_env%rtp_method == rtp_method_bse .OR. &
                     bs_env%rtp_method == rtp_method_bse_linearized
      keep_static = write_static .OR. bs_env%keep_W_freq_zero
      IF (keep_static) CALL sum_W_MIC_freq_zero(bs_env, fm_W_time)
      IF (write_static) CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '
      CALL timestop(handle)
   END SUBROUTINE compute_W

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE get_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_W_MIC'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (bs_env%all_W_exist) THEN
         CALL read_W_MIC_time(bs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      ELSE
         CALL compute_W(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param cfm_V_kp ...
!> \param ikp_batch ...
! **************************************************************************************************
   SUBROUTINE compute_V_k_by_lattice_sum(bs_env, qs_env, cfm_V_kp, ikp_batch)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_V_kp
      INTEGER                                            :: ikp_batch

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_k_by_lattice_sum'

      INTEGER                                            :: handle, ikp, ikp_end, ikp_start, &
                                                            nkp_chi_eps_W_batch, re_im
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_V_kp_work
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      nkp_chi_eps_W_batch = bs_env%nkp_chi_eps_W_batch

      ikp_start = (ikp_batch - 1)*bs_env%nkp_chi_eps_W_batch + 1
      ikp_end = MIN(ikp_batch*bs_env%nkp_chi_eps_W_batch, bs_env%kpoints_chi_eps_W%nkp)

      NULLIFY (mat_V_kp)
      ALLOCATE (mat_V_kp(ikp_start:ikp_end, 2))

      DO re_im = 1, 2
         DO ikp = ikp_start, ikp_end
            NULLIFY (mat_V_kp(ikp, re_im)%matrix)
            ALLOCATE (mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_create(mat_V_kp(ikp, re_im)%matrix, template=bs_env%mat_RI_RI%matrix)
            CALL dbcsr_reserve_all_blocks(mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_set(mat_V_kp(ikp, re_im)%matrix, 0.0_dp)
         END DO ! ikp
      END DO ! re_im

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      cell=cell, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      IF (ikp_end <= bs_env%nkp_chi_eps_W_orig) THEN

         ! 1. 2c Coulomb integrals for the first "original" k-point grid
         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      ELSE IF (ikp_start > bs_env%nkp_chi_eps_W_orig .AND. &
               ikp_end <= bs_env%nkp_chi_eps_W_orig_plus_extra) THEN

         ! 2. 2c Coulomb integrals for the second "extrapolation" k-point grid
         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_extra

      ELSE

         CPABORT("Error with k-point parallelization.")

      END IF

      CALL build_2c_coulomb_matrix_kp(mat_V_kp, &
                                      bs_env%kpoints_chi_eps_W, &
                                      basis_type="RI_AUX", &
                                      cell=cell, &
                                      particle_set=particle_set, &
                                      qs_kind_set=qs_kind_set, &
                                      atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=bs_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=ikp_start, &
                                      ikp_end=ikp_end)

      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      ALLOCATE (cfm_V_kp(ikp_start:ikp_end))
      CALL cp_fm_create(fm_V_kp_work, bs_env%fm_RI_RI%matrix_struct)
      DO ikp = ikp_start, ikp_end
         CALL cp_cfm_create(cfm_V_kp(ikp), bs_env%fm_RI_RI%matrix_struct)
         CALL cp_fm_set_all(fm_V_kp_work, 0.0_dp)
         CALL copy_dbcsr_to_fm(mat_V_kp(ikp, 1)%matrix, fm_V_kp_work)
         cfm_V_kp(ikp)%local_data = CMPLX(fm_V_kp_work%local_data, 0.0_dp, KIND=dp)

         CALL cp_fm_set_all(fm_V_kp_work, 0.0_dp)
         CALL copy_dbcsr_to_fm(mat_V_kp(ikp, 2)%matrix, fm_V_kp_work)
         cfm_V_kp(ikp)%local_data = CMPLX(REAL(cfm_V_kp(ikp)%local_data, KIND=dp), &
                                          fm_V_kp_work%local_data, KIND=dp)
      END DO
      CALL cp_fm_release(fm_V_kp_work)
      CALL dbcsr_deallocate_matrix_set(mat_V_kp)

      CALL timestop(handle)

   END SUBROUTINE compute_V_k_by_lattice_sum

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param cfm_V_kp ...
!> \param cfm_V_sqrt_ikp ...
!> \param cfm_M_inv_V_sqrt_ikp ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE compute_MinvVsqrt_Vsqrt(bs_env, qs_env, cfm_V_kp, cfm_V_sqrt_ikp, &
                                      cfm_M_inv_V_sqrt_ikp, ikp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_V_kp
      TYPE(cp_cfm_type)                                  :: cfm_V_sqrt_ikp, cfm_M_inv_V_sqrt_ikp
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_MinvVsqrt_Vsqrt'

      INTEGER                                            :: handle, n_RI
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_ikp, cfm_work
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_M_ikp

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      ! get here M(k) and write it to cfm_M_ikp
      CALL RI_2c_integral_mat_kp(qs_env, cfm_M_ikp, cfm_V_kp(ikp)%matrix_struct, n_RI, bs_env%ri_metric, &
                                 bs_env%kpoints_chi_eps_W, regularization_RI=bs_env%regularization_RI, &
                                 ikp_ext=ikp, do_build_cell_index=(ikp == 1))

      IF (ikp == 1) THEN
         CALL cp_cfm_create(cfm_V_sqrt_ikp, cfm_V_kp(ikp)%matrix_struct)
         CALL cp_cfm_create(cfm_M_inv_V_sqrt_ikp, cfm_V_kp(ikp)%matrix_struct)
      END IF
      CALL cp_cfm_create(cfm_M_inv_ikp, cfm_V_kp(ikp)%matrix_struct)

      CALL cp_cfm_to_cfm(cfm_M_ikp(1), cfm_M_inv_ikp)

      CALL cp_cfm_release(cfm_M_ikp(1))
      DEALLOCATE (cfm_M_ikp)

      CALL cp_cfm_create(cfm_work, cfm_V_kp(ikp)%matrix_struct)

      ! M(k) -> M^-1(k), retaining the metric's eigenvalue cutoff on the fallback path.
      CALL cfm_invert(cfm_M_inv_ikp, bs_env%eps_eigval_mat_RI, workspace=cfm_work)

      ! Prepare B^H(k) B(k) = V(k) and F(k) = M^-1(k) B^H(k), retaining all positive Coulomb modes.
      CALL prepare_screening_factors(cfm_V_kp(ikp), cfm_M_inv_ikp, cfm_V_sqrt_ikp, &
                                     cfm_M_inv_V_sqrt_ikp, 0.0_dp, workspace=cfm_work)
      CALL cp_cfm_release(cfm_work)

      CALL cp_cfm_release(cfm_M_inv_ikp)

      CALL timestop(handle)

   END SUBROUTINE compute_MinvVsqrt_Vsqrt

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE read_W_MIC_time(bs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'read_W_MIC_time'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: t1

      CALL timeset(routineN, handle)

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)
      CALL create_fm_W_MIC_time(bs_env, fm_W_MIC_time)

      DO i_t = 1, bs_env%num_time_freq_points

         t1 = m_walltime()

         CALL fm_read(fm_W_MIC_time(i_t), bs_env, bs_env%W_time_name, i_t)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I5,A,I3,A,F10.1,A)') &
               'Read W^MIC(iτ) from file for time point  ', i_t, ' /', bs_env%num_time_freq_points, &
               ', Execution time', m_walltime() - t1, ' s'
         END IF

      END DO

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      ! W^MIC(iω=0) for the post-GW consumers, built here because fm_W_MIC_time does not outlive
      ! the W stage: RT-BSE reads its own restart file, PROPERTIES%BSE sums the W(iτ) just read.
      IF (bs_env%rtp_method == rtp_method_bse .OR. &
          bs_env%rtp_method == rtp_method_bse_linearized) THEN
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_W_MIC_freq%matrix_struct)
         t1 = m_walltime()
         CALL fm_read(bs_env%fm_W_MIC_freq_zero, bs_env, "W_freq_rtp", 0)
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I3,A,I3,A,F10.1,A)') &
               'Read W^MIC(f=0) from file for freq. point  ', 1, ' /', 1, &
               ', Execution time', m_walltime() - t1, ' s'
         END IF
      ELSE IF (bs_env%keep_W_freq_zero) THEN
         ! no RTP file in a LRBSE-only run: W(iτ) was just read from the W_time files
         CALL sum_W_MIC_freq_zero(bs_env, fm_W_MIC_time)
      END IF

      CALL timestop(handle)

   END SUBROUTINE read_W_MIC_time

! **************************************************************************************************
!> \brief W_PQ^MIC(iω_j) from χ_PQ(iω_j,k=0) in bs_env%fm_chi_Gamma_freq
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_MIC_freq_j ...
!> \param j_w ...
!> \param ikp ...
!> \param cfm_M_inv_V_sqrt_ikp ...
!> \param cfm_V_sqrt_ikp ...
! **************************************************************************************************
   SUBROUTINE compute_fm_W_MIC_freq_j_from_chi(bs_env, qs_env, fm_W_MIC_freq_j, j_w, ikp, &
                                               cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      INTEGER                                            :: j_w, ikp
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_W_MIC_freq_j_from_chi'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_chi_eps_W_ikp_freq_j

      CALL timeset(routineN, handle)

      CALL cp_fm_set_all(fm_W_MIC_freq_j, 0.0_dp)

      ! 1. Get χ_PQ(iω_j,k_i) from χ_PQ(iω_j,k=0) using the minimum image convention
      CALL cfm_ikp_from_fm_Gamma(cfm_chi_eps_W_ikp_freq_j, bs_env%fm_chi_Gamma_freq, &
                                 ikp, qs_env, bs_env%kpoints_chi_eps_W, "RI_AUX")

      ! 2. Remove all negative eigenvalues from χ_PQ(iω_j,k_i)
      CALL cp_cfm_power(cfm_chi_eps_W_ikp_freq_j, threshold=0.0_dp, exponent=1.0_dp)

      ! 3. ε(iω_j,k_i) = Id - V^0.5(k_i)*M^-1(k_i)*χ(iω_j,k_i)*M^-1(k_i)*V^0.5(k_i)
      CALL screened_interaction(cfm_chi_eps_W_ikp_freq_j, cfm_V_sqrt_ikp, cfm_M_inv_V_sqrt_ikp, &
                                bs_env%eps_eigval_mat_RI, unit_nr=bs_env%unit_nr)

      ! 4. k-point integration W_PQ(iω_j, k_i) to W_PQ^MIC(iω_j)
      SELECT CASE (bs_env%approx_kp_extrapol)
      CASE (.FALSE.)
         ! default: standard k-point extrapolation
         CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_W_MIC_freq_j, &
                                        cfm_chi_eps_W_ikp_freq_j, ikp, &
                                        bs_env%kpoints_chi_eps_W, "RI_AUX")
      CASE (.TRUE.)
         ! for approximate kpoint extrapolation: get W_PQ^MIC(iω_1) with and without k-point
         ! extrapolation to compute the extrapolation factor f_PQ for every PQ-matrix element,
         ! f_PQ = (W_PQ^MIC(iω_1) with extrapolation) / (W_PQ^MIC(iω_1) without extrapolation)

         ! for ω_1, we compute the k-point extrapolated result using all k-points
         IF (j_w == 1) THEN

            ! k-point extrapolated
            CALL MIC_contribution_from_ikp(bs_env, qs_env, bs_env%fm_W_MIC_freq_1_extra, &
                                           cfm_chi_eps_W_ikp_freq_j, ikp, &
                                           bs_env%kpoints_chi_eps_W, &
                                           "RI_AUX")
            ! non-kpoint extrapolated
            IF (ikp <= bs_env%nkp_chi_eps_W_orig) THEN
               CALL MIC_contribution_from_ikp(bs_env, qs_env, bs_env%fm_W_MIC_freq_1_no_extra, &
                                              cfm_chi_eps_W_ikp_freq_j, ikp, &
                                              bs_env%kpoints_chi_eps_W, &
                                              "RI_AUX", wkp_ext=bs_env%wkp_orig)
            END IF

         END IF

         ! for all ω_j, we need to compute W^MIC without k-point extrpolation
         IF (ikp <= bs_env%nkp_chi_eps_W_orig) THEN
            CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_W_MIC_freq_j, &
                                           cfm_chi_eps_W_ikp_freq_j, &
                                           ikp, bs_env%kpoints_chi_eps_W, "RI_AUX", &
                                           wkp_ext=bs_env%wkp_orig)
         END IF
      END SELECT

      CALL cp_cfm_release(cfm_chi_eps_W_ikp_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_MIC_freq_j_from_chi

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE apply_extrapol_factor(bs_env, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_extrapol_factor'
      REAL(KIND=dp), PARAMETER                           :: eps_W_no_extra_1 = 1.0E-13_dp

      INTEGER                                            :: handle, i, i_t, j, ncol_local, nrow_local
      REAL(KIND=dp)                                      :: extrapol_factor, W_extra_1, W_no_extra_1

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_W_MIC_time(1), nrow_local=nrow_local, ncol_local=ncol_local)

      DO i_t = 1, SIZE(fm_W_MIC_time)
         DO j = 1, ncol_local
            DO i = 1, nrow_local

               W_extra_1 = bs_env%fm_W_MIC_freq_1_extra%local_data(i, j)
               W_no_extra_1 = bs_env%fm_W_MIC_freq_1_no_extra%local_data(i, j)

               IF (ABS(W_no_extra_1) > eps_W_no_extra_1) THEN
                  extrapol_factor = ABS(W_extra_1/W_no_extra_1)
               ELSE
                  extrapol_factor = 1.0_dp
               END IF

               ! reset extrapolation factor if it is very large
               IF (extrapol_factor > 10.0_dp) extrapol_factor = 1.0_dp

               fm_W_MIC_time(i_t)%local_data(i, j) = fm_W_MIC_time(i_t)%local_data(i, j) &
                                                     *extrapol_factor
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE apply_extrapol_factor

! **************************************************************************************************
!> \brief Computes V, V^0.5, and M^-1 V^0.5 for molecular dielectric screening.
!>        The Coulomb matrix V uses the auxiliary basis or AUTO_RI data. The inverse metric
!>        M^-1(k=0) is precomputed once in gw_utils and stored in bs_env.
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V Coulomb matrix V(k=0)
!> \param fm_V_sqrt Coulomb factor satisfying V = V_sqrt^T V_sqrt
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
                         qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

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

      ! Prepare B^T B = V and F = M^-1 B^T with the molecular Coulomb cutoff.
      CALL prepare_screening_factors(fm_V, bs_env%fm_Minv_Gamma, fm_V_sqrt, fm_Minv_Vsqrt, &
                                     bs_env%eps_eigval_mat_RI, unit_nr=bs_env%unit_nr)

      CALL timestop(handle)

   END SUBROUTINE compute_V_MinvVsqrt

END MODULE gw_screening
