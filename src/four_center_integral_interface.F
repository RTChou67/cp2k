!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Four-center integral evaluation with an encapsulated Libint backend.
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      11.2019 Fixed potential_id initial values (A. Bussy)
!> \author Manuel Guidon
! **************************************************************************************************
MODULE four_center_integral_interface

   USE cell_types,                      ONLY: cell_type,&
                                              pbc,&
                                              real_to_scaled
   USE four_center_integral_derivatives,ONLY: add_hermite_derivatives,&
                                              add_ssss_derivatives,&
                                              eri4_hermite_workspace_type,&
                                              expand_ssss_virial,&
                                              extract_hermite_derivatives,&
                                              prepare_hermite_workspace
   USE four_center_integral_preparation,ONLY: min_image_gap_sq,&
                                              prepare_eri4_pair,&
                                              prepare_eri4_primitive,&
                                              prepare_eri4_product
   USE four_center_integral_types,      ONLY: eri4_cell_type,&
                                              eri4_image_group_type,&
                                              eri4_pgf_image_type,&
                                              eri4_pgf_pair_type
   USE hfx_contraction_methods,         ONLY: contract
   USE input_constants,                 ONLY: do_potential_id,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_wrapper,                  ONLY: &
        cp_libint_cleanup_eri, cp_libint_cleanup_eri1, cp_libint_get_derivs, cp_libint_get_eris, &
        cp_libint_get_image_derivs, cp_libint_get_image_eris, cp_libint_image_workspace_type, &
        cp_libint_init_eri, cp_libint_init_eri1, cp_libint_set_contrdepth, &
        cp_libint_set_params_eri, cp_libint_set_params_eri_deriv, cp_libint_t, &
        eri_primitive_data_type, get_eri_permutation, libint_can_contract_images, reorder_eri_block
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: four_center_integral_init, four_center_integral_cleanup, &
             prepare_four_center_batch, prepare_four_center_derivative_batch, evaluate_eri, &
             evaluate_deriv_eri, &
             evaluate_eri_screen

!> \brief Context owning the energy and first-derivative backend workspaces.
!> \note Each thread owns its context. Pass it directly without intrinsic copying.
   TYPE, PUBLIC :: four_center_integral_context_type
      PRIVATE
      TYPE(cp_libint_t) :: eri, deriv
      TYPE(cp_libint_image_workspace_type), ALLOCATABLE :: image_workspace
      TYPE(eri_primitive_data_type), ALLOCATABLE :: products(:)
      TYPE(eri4_image_group_type), ALLOCATABLE :: image_groups(:)
      TYPE(eri4_hermite_workspace_type), ALLOCATABLE :: hermite
      REAL(dp), ALLOCATABLE :: derivative_blocks(:, :)
      LOGICAL :: derivatives_prepared = .FALSE., ssss_prepared = .FALSE.
      ! A, C, and D gradients followed by the lower triangle of the Cartesian strain derivative.
      REAL(dp) :: ssss_derivatives(15) = 0.0_dp
      INTEGER :: nproducts = 0, nimage_groups = 0
   END TYPE four_center_integral_context_type

CONTAINS

! **************************************************************************************************
!> \brief Select image-resolved quartets and prepare their backend cache before angular traversal.
!> \param context owning per-thread geometry and backend workspace
!> \param list1 bra primitive pair with image topology and screening bounds
!> \param list2 ket primitive pair with image topology and screening bounds
!> \param nproducts number of retained image quartets
!> \param log10_pmax logarithmic density bound
!> \param log10_eps_schwarz logarithmic integral screening threshold
!> \param neighbor_cells periodic image translations
!> \param cell simulation cell
!> \param potential_parameter physical radial operator
!> \param m_max maximum radial order including any derivative increment
!> \param do_periodic whether to enumerate periodic translations
!> \param group_images whether the caller can use Libint-contracted groups
!> \param force_image_groups whether the caller needs the groups independently of Libint support
! **************************************************************************************************
   SUBROUTINE prepare_four_center_batch(context, list1, list2, nproducts, &
                                        log10_pmax, log10_eps_schwarz, neighbor_cells, &
                                        cell, potential_parameter, m_max, do_periodic, group_images, force_image_groups)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      TYPE(eri4_pgf_pair_type), INTENT(IN)               :: list1, list2
      INTEGER, INTENT(OUT)                               :: nproducts
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz
      TYPE(eri4_cell_type), DIMENSION(:), POINTER        :: neighbor_cells
      TYPE(cell_type), POINTER                           :: cell
      TYPE(coulomb_operator_type)                        :: potential_parameter
      INTEGER, INTENT(IN)                                :: m_max
      LOGICAL, INTENT(IN)                                :: do_periodic, group_images
      LOGICAL, INTENT(IN), OPTIONAL                      :: force_image_groups

      REAL(dp), PARAMETER                                :: r_cutoff_eps = 1.0E-12_dp

      INTEGER                                            :: first, i, j, k, k_max, nimages1, &
                                                            nimages2, tmp_i4
      LOGICAL                                            :: collect_groups, has_r_cutoff
      REAL(dp) :: C11(3), EtaInv, min_img_delta, P(3), pgf_max_1, pgf_max_2, PQ(3), pq2, Q(3), R1, &
         R2, r_cutoff_sq, r_cutoff_sq_guard, rc(3), rcut_eff, rd(3), rpq2, S1234, S1234a, S1234b, &
         shift(3), temp(3), temp_CC(3), temp_DD(3), tmp_D(3), Zeta_C, Zeta_D
      TYPE(eri4_image_group_type), ALLOCATABLE           :: tmp_groups(:)
      TYPE(eri4_pgf_image_type)                          :: shifted_image
      TYPE(eri_primitive_data_type)                      :: product
      TYPE(eri_primitive_data_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: tmp_product_list

      nimages1 = list1%nimages
      nimages2 = list2%nimages
      nproducts = 0
      context%nproducts = 0
      context%nimage_groups = 0
      context%derivatives_prepared = .FALSE.
      IF (nimages1 == 0 .OR. nimages2 == 0) RETURN
      collect_groups = libint_can_contract_images .AND. group_images
      IF (PRESENT(force_image_groups)) collect_groups = collect_groups .OR. force_image_groups
      CALL prepare_eri4_product(list1, list2, m_max, product)
      EtaInv = list2%ZetaInv
      Zeta_C = list2%zeta
      Zeta_D = list2%zetb
      temp_CC = 0.0_dp
      temp_DD = 0.0_dp
      r_cutoff_sq = 0.0_dp
      r_cutoff_sq_guard = 0.0_dp
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_truncated, do_potential_short, do_potential_mix_cl_trunc)
         has_r_cutoff = .TRUE.
         rcut_eff = potential_parameter%cutoff_radius
      CASE DEFAULT
         has_r_cutoff = .FALSE.
         rcut_eff = 0.0_dp
      END SELECT
      DO i = 1, nimages1
         P = list1%image_list(i)%P
         R1 = list1%image_list(i)%R
         S1234a = list1%image_list(i)%S1234
         pgf_max_1 = list1%image_list(i)%pgf_max
         DO j = 1, nimages2
            pgf_max_2 = list2%image_list(j)%pgf_max
            IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz) CYCLE
            Q = list2%image_list(j)%P
            R2 = list2%image_list(j)%R
            S1234b = list2%image_list(j)%S1234
            rc = list2%image_list(j)%ra
            rd = list2%image_list(j)%rb

            IF (has_r_cutoff) THEN
               r_cutoff_sq = (R1 + R2 + rcut_eff)**2
               r_cutoff_sq_guard = r_cutoff_sq*(1.0_dp + r_cutoff_eps)
            END IF
            k_max = SIZE(neighbor_cells)
            IF (do_periodic) THEN
               temp = P - Q
               PQ = pbc(temp, cell)
               IF (has_r_cutoff .AND. cell%orthorhombic .AND. ALL(cell%perd == 1) .AND. k_max >= 1) THEN
                  IF (ALL(neighbor_cells(1)%cell_r == 0.0_dp)) THEN
                     pq2 = PQ(1)**2 + PQ(2)**2 + PQ(3)**2
                     IF (pq2 > r_cutoff_sq_guard) CYCLE
                     min_img_delta = min_image_gap_sq(PQ, cell%hmat)
                     IF (min_img_delta >= 0.0_dp .AND. pq2 + min_img_delta > r_cutoff_sq_guard) k_max = 1
                  END IF
               END IF
               shift = -PQ + temp
               temp_CC = rc + shift
               temp_DD = rd + shift
            END IF

            S1234 = EXP(S1234a + S1234b)
            first = nproducts + 1

            DO k = 1, k_max
               IF (do_periodic) THEN
                  IF (has_r_cutoff .AND. k_max > 1) THEN
                     IF ((PQ(1) - neighbor_cells(k)%cell_r(1))**2 + &
                         (PQ(2) - neighbor_cells(k)%cell_r(2))**2 + &
                         (PQ(3) - neighbor_cells(k)%cell_r(3))**2 > r_cutoff_sq_guard) CYCLE
                  END IF
                  C11 = temp_CC + neighbor_cells(k)%cell_r(:)
                  tmp_D = temp_DD + neighbor_cells(k)%cell_r(:)
               ELSE
                  C11 = rc
                  tmp_D = rd
               END IF
               Q = (Zeta_C*C11 + Zeta_D*tmp_D)*EtaInv
               rpq2 = (P(1) - Q(1))**2 + (P(2) - Q(2))**2 + (P(3) - Q(3))**2
               IF (has_r_cutoff) THEN
                  IF (rpq2 > r_cutoff_sq) CYCLE
               END IF
               nproducts = nproducts + 1

               ! Retain the enlarged buffer in this thread's owning context.
               IF (nproducts > SIZE(context%products)) THEN
                  tmp_i4 = (3*nproducts + 1)/2
                  ALLOCATE (tmp_product_list(tmp_i4))
                  tmp_product_list(1:SIZE(context%products)) = context%products
                  CALL MOVE_ALLOC(tmp_product_list, context%products)
               END IF

               shifted_image%ra = C11
               shifted_image%rb = tmp_D
               shifted_image%P = Q
               CALL prepare_eri4_primitive(product, list1%image_list(i), shifted_image, S1234, rpq2, potential_parameter, &
                                           context%products(nproducts))
            END DO

            IF (collect_groups .AND. do_periodic .AND. nproducts > first) THEN
               IF (.NOT. ALLOCATED(context%image_groups)) ALLOCATE (context%image_groups(16))
               IF (context%nimage_groups == SIZE(context%image_groups)) THEN
                  ALLOCATE (tmp_groups((3*context%nimage_groups + 2)/2))
                  tmp_groups(1:context%nimage_groups) = context%image_groups
                  CALL MOVE_ALLOC(tmp_groups, context%image_groups)
               END IF
               context%nimage_groups = context%nimage_groups + 1
               ASSOCIATE (group => context%image_groups(context%nimage_groups))
                  group%first = first
                  group%last = nproducts
                  group%displacements(:, 1) = list1%image_list(i)%ra - list1%image_list(i)%rb
                  group%displacements(:, 2) = rc - rd
               END ASSOCIATE
            END IF
         END DO
      END DO

      context%nproducts = nproducts
   END SUBROUTINE prepare_four_center_batch

! **************************************************************************************************
!> \brief Prepare a primitive quartet and its Cartesian derivatives when image sums can be reused.
!> \param context owning per-thread geometry and derivative workspace
!> \param list1 bra primitive pair
!> \param list2 ket primitive pair
!> \param nproducts number of retained image quartets
!> \param log10_pmax logarithmic density bound
!> \param log10_eps_schwarz logarithmic integral screening threshold
!> \param neighbor_cells periodic image translations
!> \param cell simulation cell
!> \param potential_parameter physical radial operator
!> \param lower minimum angular momentum at each center
!> \param upper maximum angular momentum at each center
!> \param do_periodic whether to enumerate periodic translations
!> \param use_virial whether to compute the virial tensor
! **************************************************************************************************
   SUBROUTINE prepare_four_center_derivative_batch(context, list1, list2, nproducts, &
                                                   log10_pmax, log10_eps_schwarz, neighbor_cells, &
                                                   cell, potential_parameter, lower, upper, do_periodic, use_virial)
      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      TYPE(eri4_pgf_pair_type), INTENT(IN)               :: list1, list2
      INTEGER, INTENT(OUT)                               :: nproducts
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz
      TYPE(eri4_cell_type), POINTER                      :: neighbor_cells(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      INTEGER, INTENT(IN)                                :: lower(4), upper(4)
      LOGICAL, INTENT(IN)                                :: do_periodic, use_virial

      INTEGER, PARAMETER                                 :: min_images_per_group = 4

      INTEGER                                            :: group, i, j, ncart, ncomponents, ngroups
      LOGICAL                                            :: hermite, radial_moments, ssss
      REAL(dp)                                           :: abcd(3, 2)

      context%ssss_prepared = .FALSE.
      radial_moments = potential_parameter%potential_type /= do_potential_id
      ssss = radial_moments .AND. ALL(upper == 0)
      hermite = radial_moments .AND. do_periodic .AND. ALL(upper <= 2) .AND. .NOT. ssss
      CALL prepare_four_center_batch(context, list1, list2, nproducts, log10_pmax, log10_eps_schwarz, neighbor_cells, &
                                     cell, potential_parameter, SUM(upper) + 1, do_periodic, &
                                     group_images=.NOT. use_virial, force_image_groups=ssss .OR. hermite)
      IF (nproducts == 0) RETURN

      IF (hermite) THEN
         ngroups = nproducts
         DO group = 1, context%nimage_groups
            ngroups = ngroups - (context%image_groups(group)%last - context%image_groups(group)%first)
         END DO
         hermite = nproducts/ngroups >= min_images_per_group
      END IF
      IF (.NOT. (ssss .OR. hermite)) RETURN

      IF (ssss) THEN
         context%ssss_derivatives = 0.0_dp
      ELSE
         ncart = PRODUCT(ncoset(upper) - ncoset(lower - 1))
         ncomponents = MERGE(21, 12, use_virial)
         IF (ALLOCATED(context%derivative_blocks)) THEN
            IF (SIZE(context%derivative_blocks, 1) < ncart .OR. SIZE(context%derivative_blocks, 2) < ncomponents) THEN
               DEALLOCATE (context%derivative_blocks)
            END IF
         END IF
         IF (.NOT. ALLOCATED(context%derivative_blocks)) ALLOCATE (context%derivative_blocks(ncart, ncomponents))
         context%derivative_blocks(:ncart, :ncomponents) = 0.0_dp
         IF (.NOT. ALLOCATED(context%hermite)) ALLOCATE (context%hermite)
         CALL prepare_hermite_workspace(context%hermite, lower, upper, use_virial)
      END IF

      i = 1
      group = 1
      DO WHILE (i <= nproducts)
         j = i
         IF (group <= context%nimage_groups) THEN
            IF (i == context%image_groups(group)%first) THEN
               j = context%image_groups(group)%last
               abcd = context%image_groups(group)%displacements
               group = group + 1
            END IF
         END IF
         IF (j == i) THEN
            abcd(:, 1) = context%products(i)%centers(:, 1) - context%products(i)%centers(:, 2)
            abcd(:, 2) = context%products(i)%centers(:, 3) - context%products(i)%centers(:, 4)
         END IF
         IF (ssss) THEN
            CALL add_ssss_derivatives(context%products(i:j), abcd, use_virial, context%ssss_derivatives)
         ELSE
            CALL add_hermite_derivatives(context%hermite, context%products(i:j), abcd, cell, context%derivative_blocks, ncart)
         END IF
         i = j + 1
      END DO
      context%ssss_prepared = ssss
      context%derivatives_prepared = .TRUE.
   END SUBROUTINE prepare_four_center_derivative_batch

! **************************************************************************************************
!> \brief Initialize the four-center energy and first-derivative workspaces.
!> \param context owning context; global backend initialization is provided by init_cp2k
!> \param max_am maximum angular momentum on any center
! **************************************************************************************************
   SUBROUTINE four_center_integral_init(context, max_am)
      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      INTEGER, INTENT(IN)                                :: max_am

      CALL cp_libint_init_eri(context%eri, max_am)
      CALL cp_libint_init_eri1(context%deriv, max_am)
      CALL cp_libint_set_contrdepth(context%eri, 1)
      CALL cp_libint_set_contrdepth(context%deriv, 1)
      ALLOCATE (context%products(128))
      IF (libint_can_contract_images) ALLOCATE (context%image_workspace, context%image_groups(16))
   END SUBROUTINE four_center_integral_init

! **************************************************************************************************
!> \brief Release the workspaces of an owning four-center context.
!> \param context original initialized owner
! **************************************************************************************************
   SUBROUTINE four_center_integral_cleanup(context)
      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context

      IF (ALLOCATED(context%image_workspace)) DEALLOCATE (context%image_workspace)
      IF (ALLOCATED(context%image_groups)) DEALLOCATE (context%image_groups)
      IF (ALLOCATED(context%hermite)) DEALLOCATE (context%hermite)
      IF (ALLOCATED(context%derivative_blocks)) DEALLOCATE (context%derivative_blocks)
      context%derivatives_prepared = .FALSE.
      context%ssss_prepared = .FALSE.
      CALL cp_libint_cleanup_eri(context%eri)
      CALL cp_libint_cleanup_eri1(context%deriv)
      DEALLOCATE (context%products)
      context%nproducts = 0
      context%nimage_groups = 0
   END SUBROUTINE four_center_integral_cleanup

! **************************************************************************************************
!> \brief Evaluates derivatives of electron repulsion integrals for a primitive quartet
!> \param context four-center integral context
!> \param n_a ...
!> \param n_b ...
!> \param n_c ...
!> \param n_d ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param ncod ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
!> \param primitives ...
!> \param max_contraction ...
!> \param tmp_max_all ...
!> \param eps_schwarz ...
!> \param neris ...
!> \param s_offset_a ...
!> \param s_offset_b ...
!> \param s_offset_c ...
!> \param s_offset_d ...
!> \param nl_a ...
!> \param nl_b ...
!> \param nl_c ...
!> \param nl_d ...
!> \param nsoa ...
!> \param nsob ...
!> \param nsoc ...
!> \param nsod ...
!> \param angular_a ...
!> \param angular_b ...
!> \param angular_c ...
!> \param angular_d ...
!> \param weight_a radial weights of the selected primitive on A
!> \param weight_b radial weights of the selected primitive on B
!> \param weight_c radial weights of the selected primitive on C
!> \param weight_d radial weights of the selected primitive on D
!> \param work ...
!> \param work2 ...
!> \param work_forces ...
!> \param buffer1 ...
!> \param buffer2 ...
!> \param use_virial ...
!> \param work_virial ...
!> \param work2_virial ...
!> \param primitives_virial ...
!> \param cell ...
!> \param tmp_max_all_virial ...
!> \par History
!>      03.2007 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_deriv_eri(context, &
                                 n_a, n_b, n_c, n_d, &
                                 ncoa, ncob, ncoc, ncod, &
                                 nsgfa, nsgfb, nsgfc, nsgfd, &
                                 primitives, max_contraction, tmp_max_all, &
                                 eps_schwarz, neris, &
                                 s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                                 nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod, &
                                 angular_a, angular_b, angular_c, angular_d, &
                                 weight_a, weight_b, weight_c, weight_d, &
                                 work, work2, work_forces, buffer1, buffer2, &
                                 use_virial, work_virial, work2_virial, &
                                 primitives_virial, cell, tmp_max_all_virial)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d, ncoa, ncob, ncoc, &
                                                            ncod, nsgfa, nsgfb, nsgfc, nsgfd
      REAL(dp), &
         DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd, 12)       :: primitives
      REAL(dp)                                           :: max_contraction, tmp_max_all, eps_schwarz
      INTEGER(int_8)                                     :: neris
      INTEGER                                            :: s_offset_a, s_offset_b, s_offset_c, &
                                                            s_offset_d, nl_a, nl_b, nl_c, nl_d, &
                                                            nsoa, nsob, nsoc, nsod
      REAL(dp), DIMENSION(ncoa, nsoa), INTENT(IN)        :: angular_a
      REAL(dp), DIMENSION(ncob, nsob), INTENT(IN)        :: angular_b
      REAL(dp), DIMENSION(ncoc, nsoc), INTENT(IN)        :: angular_c
      REAL(dp), DIMENSION(ncod, nsod), INTENT(IN)        :: angular_d
      REAL(dp), INTENT(IN)                               :: weight_a(nl_a), weight_b(nl_b), &
                                                            weight_c(nl_c), weight_d(nl_d)
      REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod, 12)       :: work, work2, work_forces
      REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod)           :: buffer1, buffer2
      LOGICAL, INTENT(IN)                                :: use_virial
      REAL(dp) :: work_virial(ncoa*ncob*ncoc*ncod, 3, 3), work2_virial(ncoa*ncob*ncoc*ncod, 3, 3)
      REAL(dp), &
         DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd, 3, 3)     :: primitives_virial
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp)                                           :: tmp_max_all_virial

      INTEGER                                            :: angular_momenta(4), group_index, i, j, &
                                                            k, m, mysize, next_group, &
                                                            permutation(4)
      LOGICAL                                            :: compact_ssss
      REAL(dp)                                           :: cell_derivative(3, 3), cell_scale(3), &
                                                            force_max, origin(3), relative(3, 3), &
                                                            strain(6), stress_max

      mysize = ncoa*ncob*ncoc*ncod
      compact_ssss = context%derivatives_prepared .AND. context%ssss_prepared

      IF (context%derivatives_prepared) THEN
         IF (compact_ssss) THEN
            work2(1, 1:3) = context%ssss_derivatives(1:3)
            work2(1, 7:12) = context%ssss_derivatives(4:9)
            work2(1, 4:6) = -(work2(1, 1:3) + work2(1, 7:9) + work2(1, 10:12))
            force_max = MAXVAL(ABS(work2(1, :)))
            IF (use_virial) THEN
               strain = context%ssss_derivatives(10:15)
               work2_virial(1, :, 1) = strain(1:3)
               work2_virial(1, 2:3, 2) = strain(4:5)
               work2_virial(1, 3, 3) = strain(6)
               ! Screen the full cell derivative and contract only the strain lower triangle.
               IF (cell%orthorhombic) THEN
                  cell_scale = ABS([cell%h_inv(1, 1), cell%h_inv(2, 2), cell%h_inv(3, 3)])
                  stress_max = MAXVAL(ABS(strain)*[cell_scale(1), MAX(cell_scale(1), cell_scale(2)), &
                                                   MAX(cell_scale(1), cell_scale(3)), cell_scale(2), &
                                                   MAX(cell_scale(2), cell_scale(3)), cell_scale(3)])
               ELSE
                  cell_derivative(:, 1) = strain(1:3)
                  cell_derivative(2:3, 2) = strain(4:5)
                  cell_derivative(3, 3) = strain(6)
                  CALL expand_ssss_virial(1, cell_derivative, cell)
                  stress_max = MAXVAL(ABS(cell_derivative))
               END IF
            END IF
         ELSE IF (n_a + n_b + n_c + n_d == 0) THEN
            work2(1, :) = context%derivative_blocks(1, :12)
            force_max = MAXVAL(ABS(work2))
            IF (use_virial) THEN
               work2_virial(1, :, :) = RESHAPE(context%derivative_blocks(1, 13:21), [3, 3])
               stress_max = MAXVAL(ABS(work2_virial))
            END IF
         ELSE
            CALL extract_hermite_derivatives(context%hermite, [n_a, n_b, n_c, n_d], context%derivative_blocks, &
                                             use_virial, work2, work2_virial, force_max, stress_max)
         END IF
         tmp_max_all = MAX(tmp_max_all, max_contraction*force_max)
         IF (use_virial) tmp_max_all_virial = MAX(tmp_max_all_virial, max_contraction*stress_max)
      ELSE
         work = 0.0_dp

         IF (use_virial) THEN
            work_virial = 0.0_dp
         END IF

         CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
         angular_momenta = [n_a, n_b, n_c, n_d]
         angular_momenta = angular_momenta(permutation)

         i = 1
         group_index = 1
         next_group = 0
         IF (libint_can_contract_images .AND. .NOT. use_virial .AND. context%nimage_groups > 0) THEN
            next_group = context%image_groups(1)%first
         END IF
         DO WHILE (i <= context%nproducts)
            IF (i == next_group) THEN
               ASSOCIATE (group => context%image_groups(group_index))
                  CALL cp_libint_get_image_derivs(context%deriv, context%image_workspace, angular_momenta, &
                                                  context%products(i:group%last), &
                                                  permutation, group%displacements, work_forces)
                  i = group%last
               END ASSOCIATE
               group_index = group_index + 1
               next_group = 0
               IF (group_index <= context%nimage_groups) next_group = context%image_groups(group_index)%first
            ELSE
               CALL cp_libint_set_params_eri_deriv(context%deriv, context%products(i), &
                                                   permutation, SUM(angular_momenta) + 1)
               CALL cp_libint_get_derivs(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                         angular_momenta(1), context%deriv, work_forces)
            END IF
            DO k = 1, 12
               DO j = 1, mysize
                  work(j, k) = work(j, k) + work_forces(j, k)
               END DO
            END DO
            IF (use_virial) THEN
               ! The four center derivatives sum to zero. Use the canonical second center as the origin.
               ASSOCIATE (centers => context%products(i)%centers)
                  origin = centers(:, permutation(2))
                  CALL real_to_scaled(relative(:, 1), centers(:, permutation(1)) - origin, cell)
                  CALL real_to_scaled(relative(:, 2), centers(:, permutation(3)) - origin, cell)
                  CALL real_to_scaled(relative(:, 3), centers(:, permutation(4)) - origin, cell)
               END ASSOCIATE
               DO m = 1, 3
                  DO k = 1, 3
                     DO j = 1, mysize
                        work_virial(j, k, m) = work_virial(j, k, m) + work_forces(j, k)*relative(m, 1) + &
                                               work_forces(j, k + 6)*relative(m, 2) + work_forces(j, k + 9)*relative(m, 3)
                     END DO
                  END DO
               END DO
            END IF
            i = i + 1
         END DO

         tmp_max_all = MAX(tmp_max_all, max_contraction*MAXVAL(ABS(work)))
         IF (use_virial) THEN
            tmp_max_all_virial = MAX(tmp_max_all_virial, max_contraction*MAXVAL(ABS(work_virial)))
         END IF

      END IF
      neris = neris + INT(context%nproducts, int_8)*12*mysize

      IF (tmp_max_all >= eps_schwarz) THEN
         IF (.NOT. context%derivatives_prepared) THEN
            CALL reorder_eri_block(work, work2, [ncoa, ncob, ncoc, ncod], permutation, 12, 1)
         END IF
         DO i = 1, 12
            ! forces4 recovers center B after contraction for compact ssss and mixed shell sets.
            IF ((mysize > 1 .OR. context%derivatives_prepared) .AND. i >= 4 .AND. i <= 6) CYCLE
            CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                          n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2(1, i), &
                          angular_a, angular_b, angular_c, angular_d, &
                          weight_a, weight_b, weight_c, weight_d, nsgfa, nsgfb, nsgfc, &
                          primitives(s_offset_a + 1, s_offset_b + 1, s_offset_c + 1, s_offset_d + 1, i), &
                          buffer1, buffer2)
         END DO
      END IF

      IF (use_virial .AND. tmp_max_all_virial >= eps_schwarz) THEN
         IF (.NOT. context%derivatives_prepared) THEN
            CALL reorder_eri_block(work_virial, work2_virial, [ncoa, ncob, ncoc, ncod], permutation, 3, 3)
         END IF
         DO i = 1, 3
            DO m = 1, MERGE(i, 3, compact_ssss)
               CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                             n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2_virial(1, i, m), &
                             angular_a, angular_b, angular_c, angular_d, &
                             weight_a, weight_b, weight_c, weight_d, nsgfa, nsgfb, nsgfc, &
                             primitives_virial(s_offset_a + 1, s_offset_b + 1, s_offset_c + 1, s_offset_d + 1, i, m), &
                             buffer1, buffer2)
            END DO
         END DO
      END IF

   END SUBROUTINE evaluate_deriv_eri

! **************************************************************************************************
!> \brief Evaluates max-abs values of  electron repulsion integrals for a primitive quartet
!> \param context four-center integral context
!> \param A ...
!> \param B ...
!> \param C ...
!> \param D ...
!> \param Zeta_A ...
!> \param Zeta_B ...
!> \param Zeta_C ...
!> \param Zeta_D ...
!> \param n_a ...
!> \param n_b ...
!> \param n_c ...
!> \param n_d ...
!> \param max_val ...
!> \param potential_parameter ...
!> \param R1 ...
!> \param R2 ...
!> \par History
!>      03.2007 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_eri_screen(context, A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
                                  n_a, n_b, n_c, n_d, &
                                  max_val, potential_parameter, R1, R2)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      REAL(dp), INTENT(IN)                               :: A(3), B(3), C(3), D(3), Zeta_A, Zeta_B, &
                                                            Zeta_C, Zeta_D
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d
      REAL(dp), INTENT(INOUT)                            :: max_val
      TYPE(coulomb_operator_type)                        :: potential_parameter
      REAL(dp)                                           :: R1, R2

      INTEGER                                            :: angular_momenta(4), i, m_max, mysize, &
                                                            permutation(4)
      REAL(dp)                                           :: centers(3, 4), exponents(4), PQ(3), PQ2, &
                                                            S1234
      REAL(dp), POINTER                                  :: p_work(:)
      TYPE(eri4_pgf_image_type)                          :: image_a, image_b
      TYPE(eri4_pgf_pair_type)                           :: pair_a, pair_b
      TYPE(eri_primitive_data_type)                      :: product

      m_max = n_a + n_b + n_c + n_d
      mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)

      angular_momenta = [n_a, n_b, n_c, n_d]
      centers(:, 1) = A
      centers(:, 2) = B
      centers(:, 3) = C
      centers(:, 4) = D
      exponents = [Zeta_A, Zeta_B, Zeta_C, Zeta_D]

      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = angular_momenta(permutation)
      centers = centers(:, permutation)
      exponents = exponents(permutation)

      PQ = centers(:, 1) - centers(:, 2)
      CALL prepare_eri4_pair(centers(:, 1), centers(:, 2), exponents(1), exponents(2), &
                             DOT_PRODUCT(PQ, PQ), pair_a, image_a)
      PQ = centers(:, 3) - centers(:, 4)
      CALL prepare_eri4_pair(centers(:, 3), centers(:, 4), exponents(3), exponents(4), &
                             DOT_PRODUCT(PQ, PQ), pair_b, image_b)
      PQ = image_a%P - image_b%P
      PQ2 = DOT_PRODUCT(PQ, PQ)
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_truncated, do_potential_short, do_potential_mix_cl_trunc)
         IF (PQ2 > (R1 + R2 + potential_parameter%cutoff_radius)**2) RETURN
      END SELECT

      CALL prepare_eri4_product(pair_a, pair_b, m_max, product)
      S1234 = EXP(image_a%S1234 + image_b%S1234)
      CALL prepare_eri4_primitive(product, image_a, image_b, S1234, PQ2, potential_parameter, context%products(1))
      IF (m_max == 0) THEN
         max_val = ABS(context%products(1)%Fm(1))
      ELSE
         CALL cp_libint_set_params_eri(context%eri, context%products(1), [1, 2, 3, 4], m_max)
         CALL cp_libint_get_eris(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                 angular_momenta(1), context%eri, p_work)
         DO i = 1, mysize
            max_val = MAX(max_val, ABS(p_work(i)))
         END DO
      END IF

   END SUBROUTINE evaluate_eri_screen

! **************************************************************************************************
!> \brief Evaluate electron repulsion integrals for a primitive quartet
!> \param context four-center integral context
!> \param n_a ...
!> \param n_b ...
!> \param n_c ...
!> \param n_d ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param ncod ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
!> \param primitives ...
!> \param max_contraction ...
!> \param tmp_max ...
!> \param eps_schwarz ...
!> \param neris ...
!> \param s_offset_a ...
!> \param s_offset_b ...
!> \param s_offset_c ...
!> \param s_offset_d ...
!> \param nl_a ...
!> \param nl_b ...
!> \param nl_c ...
!> \param nl_d ...
!> \param nsoa ...
!> \param nsob ...
!> \param nsoc ...
!> \param nsod ...
!> \param angular_a ...
!> \param angular_b ...
!> \param angular_c ...
!> \param angular_d ...
!> \param weight_a radial weights of the selected primitive on A
!> \param weight_b radial weights of the selected primitive on B
!> \param weight_c radial weights of the selected primitive on C
!> \param weight_d radial weights of the selected primitive on D
!> \param work2 ...
!> \param buffer1 ...
!> \param buffer2 ...
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_eri(context, &
                           n_a, n_b, n_c, n_d, &
                           ncoa, ncob, ncoc, ncod, &
                           nsgfa, nsgfb, nsgfc, nsgfd, &
                           primitives, max_contraction, tmp_max, &
                           eps_schwarz, neris, &
                           s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                           nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod, &
                           angular_a, angular_b, angular_c, angular_d, &
                           weight_a, weight_b, weight_c, weight_d, work2, buffer1, buffer2)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d, ncoa, ncob, ncoc, &
                                                            ncod, nsgfa, nsgfb, nsgfc, nsgfd
      REAL(dp), DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd)    :: primitives
      REAL(dp)                                           :: max_contraction, tmp_max, eps_schwarz
      INTEGER(int_8)                                     :: neris
      INTEGER                                            :: s_offset_a, s_offset_b, s_offset_c, &
                                                            s_offset_d, nl_a, nl_b, nl_c, nl_d, &
                                                            nsoa, nsob, nsoc, nsod
      REAL(dp), DIMENSION(ncoa, nsoa), INTENT(IN)        :: angular_a
      REAL(dp), DIMENSION(ncob, nsob), INTENT(IN)        :: angular_b
      REAL(dp), DIMENSION(ncoc, nsoc), INTENT(IN)        :: angular_c
      REAL(dp), DIMENSION(ncod, nsod), INTENT(IN)        :: angular_d
      REAL(dp), INTENT(IN)                               :: weight_a(nl_a), weight_b(nl_b), &
                                                            weight_c(nl_c), weight_d(nl_d)
      REAL(dp)                                           :: work2(ncoa, ncob, ncoc, ncod)
      REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod), TARGET   :: buffer1
      REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod)           :: buffer2

      INTEGER                                            :: angular_momenta(4), group_index, i, &
                                                            m_max, mysize, next_group, &
                                                            permutation(4)
      REAL(dp), POINTER                                  :: p_work(:)

      IF (context%nproducts == 0) RETURN
      mysize = ncoa*ncob*ncoc*ncod
      m_max = n_a + n_b + n_c + n_d
      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = [n_a, n_b, n_c, n_d]
      angular_momenta = angular_momenta(permutation)
      IF (m_max /= 0) THEN
         IF (context%nproducts > 1) buffer1 = 0.0_dp
         i = 1
         group_index = 1
         next_group = 0
         IF (context%nimage_groups > 0) next_group = context%image_groups(1)%first
         DO WHILE (i <= context%nproducts)
            IF (i == next_group) THEN
               ASSOCIATE (group => context%image_groups(group_index))
                  CALL cp_libint_get_image_eris(context%eri, context%image_workspace, angular_momenta, &
                                                context%products(i:group%last), &
                                                permutation, group%displacements, p_work)
                  i = group%last
               END ASSOCIATE
               group_index = group_index + 1
               next_group = 0
               IF (group_index <= context%nimage_groups) next_group = context%image_groups(group_index)%first
            ELSE
               CALL cp_libint_set_params_eri(context%eri, context%products(i), permutation, m_max)
               CALL cp_libint_get_eris(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                       angular_momenta(1), context%eri, p_work)
            END IF
            IF (context%nproducts > 1) buffer1 = buffer1 + p_work(1:mysize)
            i = i + 1
         END DO
         IF (context%nproducts > 1) p_work => buffer1
         tmp_max = MAX(tmp_max, MAXVAL(ABS(p_work)))
      ELSE
         work2 = 0.0_dp
         DO i = 1, context%nproducts
            work2(1, 1, 1, 1) = work2(1, 1, 1, 1) + context%products(i)%Fm(1)
         END DO
         tmp_max = MAX(tmp_max, ABS(work2(1, 1, 1, 1)))
      END IF
      neris = neris + INT(context%nproducts, int_8)*mysize
      tmp_max = tmp_max*max_contraction
      IF (tmp_max < eps_schwarz) RETURN
      IF (m_max /= 0) CALL reorder_eri_block(p_work, work2, [ncoa, ncob, ncoc, ncod], permutation, 1, 1)

      CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                    n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2, &
                    angular_a, angular_b, angular_c, angular_d, &
                    weight_a, weight_b, weight_c, weight_d, nsgfa, nsgfb, nsgfc, &
                    primitives(s_offset_a + 1, s_offset_b + 1, s_offset_c + 1, s_offset_d + 1), &
                    buffer1, buffer2)

   END SUBROUTINE evaluate_eri

END MODULE four_center_integral_interface
