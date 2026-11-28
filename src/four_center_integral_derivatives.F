!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Derivatives from image sums of prepared radial fundamentals.
! **************************************************************************************************
MODULE four_center_integral_derivatives
   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE four_center_integral_types,      ONLY: eri4_contraction_type,&
                                              eri4_density_type,&
                                              eri4_image_group_type
   USE kinds,                           ONLY: dp
   USE libint_wrapper,                  ONLY: eri_primitive_data_type
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset,&
                                              nso
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: prepare_hermite_derivatives, extract_hermite_derivatives, expand_symmetric_virial, &
             prepare_hermite_density, select_hermite_density, ssss_density_weight

   ! Polynomial degrees for a first derivative of (dd|dd) and its coordinate moment.
   INTEGER, PARAMETER :: max_order = 9, max_monomials = 286, max_hermites = 220, max_pair_hermites = 56
   INTEGER, PARAMETER :: max_cartesian = 10**4, max_pair_cartesian = 10**2
   INTEGER, PARAMETER :: max_radial_moments = (max_order + 1)*max_monomials
   INTEGER, PARAMETER :: max_polynomial_terms = 2*max_monomials*(max_order/2 + 2)

!> \brief Reusable polynomial, basis transform and gather data for one angular range.
   TYPE :: eri4_hermite_workspace_type
      PRIVATE
      INTEGER :: lower(4) = -1, upper(4) = -1, order = 0, nherm = 0, nshift = 0, nradial = 0
      LOGICAL :: virial = .FALSE.
      INTEGER :: powers(3, max_monomials) = 0
      INTEGER :: parent(max_monomials) = 0, axis(max_monomials) = 0
      INTEGER :: monomial_start(0:max_order) = 0, radial_offset(0:max_order + 1) = 0
      INTEGER :: term_offset(max_hermites + max_monomials + 1) = 0
      INTEGER :: term_index(max_polynomial_terms) = 0
      REAL(dp) :: term_factor(max_polynomial_terms) = 0.0_dp
      INTEGER :: lower_index(3, max_hermites) = 0, upper_index(3, max_hermites) = 0
      INTEGER :: pair_index(max_pair_hermites, max_pair_hermites) = 0
      INTEGER :: moment_count(max_hermites) = 0
      REAL(dp) :: ket_sign(max_pair_hermites) = 0.0_dp
      INTEGER :: block_start(0:80) = 0, cartesian_index(max_cartesian) = 0
      REAL(dp) :: bra(max_pair_hermites, max_pair_cartesian, 0:3) = 0.0_dp
      REAL(dp) :: ket(max_pair_hermites, max_pair_cartesian, 0:6) = 0.0_dp
      REAL(dp) :: auxiliary(max_pair_hermites, max_pair_hermites) = 0.0_dp
      REAL(dp) :: intermediate(max_pair_hermites*max_pair_cartesian) = 0.0_dp
      REAL(dp) :: force(max_cartesian, 12) = 0.0_dp
   END TYPE eri4_hermite_workspace_type

!> \brief Derivative blocks for a prepared primitive quartet.
   TYPE, PUBLIC :: eri4_hermite_batch_type
      PRIVATE
      TYPE(eri4_hermite_workspace_type), ALLOCATABLE :: workspace
      REAL(dp), ALLOCATABLE :: cartesian(:, :)
   END TYPE eri4_hermite_batch_type

!> \brief Density for a shell quartet, transformed once and reused across primitive quartets.
   TYPE :: eri4_hermite_density_block
      REAL(dp), ALLOCATABLE :: values(:, :, :), norms(:, :, :)
   END TYPE eri4_hermite_density_block

   TYPE :: eri4_density_basis_type
      REAL(dp), ALLOCATABLE :: values(:, :)
      REAL(dp), CONTIGUOUS, POINTER :: radial(:, :) => NULL()
   END TYPE eri4_density_basis_type

   TYPE, PUBLIC :: eri4_hermite_density_type
      PRIVATE
      TYPE(eri4_hermite_density_block) :: pair(4)
      INTEGER :: nterms = 0, lhs(6) = 0, rhs(6) = 0, offset(4) = 0
      REAL(dp) :: factor = 0.0_dp, norm = 0.0_dp, remaining = 0.0_dp
      REAL(dp), ALLOCATABLE :: weight(:), ssss(:, :, :, :)
      LOGICAL :: weights_ready = .FALSE.
   END TYPE eri4_hermite_density_type

CONTAINS

! **************************************************************************************************
!> \brief Accumulate the retained displacement groups with the Hermite or all-s kernel.
!> \param batch reusable derivative buffers
!> \param records prepared quartet images
!> \param groups ranges and displacements of multi-image groups
!> \param ngroups number of active entries in groups
!> \param lower minimum angular momentum per center
!> \param upper maximum angular momentum per center
!> \param symmetric_strain use the all-s specialization and accumulate the lower strain triangle
!> \param cell simulation cell
!> \param use_virial whether to include cell derivatives
!> \param density primitive exchange density and screening budget
!> \param virial density-weighted cell derivative accumulated by Hermite
!> \param all_s independent all-s force and strain derivatives
! **************************************************************************************************
   SUBROUTINE prepare_hermite_derivatives(batch, records, groups, ngroups, lower, upper, symmetric_strain, &
                                          cell, use_virial, density, virial, all_s)
      TYPE(eri4_hermite_batch_type), ALLOCATABLE, &
         INTENT(INOUT)                                   :: batch
      TYPE(eri_primitive_data_type), INTENT(IN)          :: records(:)
      TYPE(eri4_image_group_type), ALLOCATABLE, &
         INTENT(IN)                                      :: groups(:)
      INTEGER, INTENT(IN)                                :: ngroups, lower(4), upper(4)
      LOGICAL, INTENT(IN)                                :: symmetric_strain
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(eri4_hermite_density_type), INTENT(INOUT)     :: density
      REAL(dp), INTENT(INOUT)                            :: virial(3, 3)
      REAL(dp), INTENT(OUT)                              :: all_s(15)

      INTEGER                                            :: group, i, j, ncart
      REAL(dp)                                           :: abcd(3, 2)

      IF (symmetric_strain) THEN
         all_s = 0.0_dp
      ELSE
         IF (.NOT. ALLOCATED(batch)) ALLOCATE (batch)
         ncart = PRODUCT(ncoset(upper) - ncoset(lower - 1))
         IF (ALLOCATED(batch%cartesian)) THEN
            IF (SIZE(batch%cartesian, 1) < ncart) DEALLOCATE (batch%cartesian)
         END IF
         IF (.NOT. ALLOCATED(batch%cartesian)) ALLOCATE (batch%cartesian(ncart, 12))
         batch%cartesian(:ncart, :) = 0.0_dp
         IF (.NOT. ALLOCATED(batch%workspace)) ALLOCATE (batch%workspace)
         CALL prepare_hermite_workspace(batch%workspace, lower, upper, use_virial)
      END IF

      i = 1
      group = 1
      DO WHILE (i <= SIZE(records))
         j = i
         IF (group <= ngroups) THEN
            IF (i == groups(group)%first) THEN
               j = groups(group)%last
               abcd = groups(group)%displacements
               group = group + 1
            END IF
         END IF
         IF (j == i) THEN
            abcd(:, 1) = records(i)%centers(:, 1) - records(i)%centers(:, 2)
            abcd(:, 2) = records(i)%centers(:, 3) - records(i)%centers(:, 4)
         END IF
         IF (symmetric_strain) THEN
            CALL add_ssss_derivatives(records(i:j), abcd, use_virial, all_s)
         ELSE
            CALL add_hermite_derivatives(batch%workspace, records(i:j), abcd, cell, batch%cartesian, &
                                         ncart, density, virial)
         END IF
         i = j + 1
      END DO
   END SUBROUTINE prepare_hermite_derivatives

! **************************************************************************************************
!> \brief Accumulate forces on A, C, and D and the lower triangle of the all-s strain derivative.
!> \param records prepared quartets with common exponents and internal pair displacements
!> \param abcd fixed internal bra and ket displacements
!> \param virial whether to accumulate strain components
!> \param derivatives Cartesian force derivatives for A, C, D followed by the strain lower triangle
! **************************************************************************************************
   SUBROUTINE add_ssss_derivatives(records, abcd, virial, derivatives)
      TYPE(eri_primitive_data_type), INTENT(IN)          :: records(:)
      REAL(dp), INTENT(IN)                               :: abcd(3, 2)
      LOGICAL, INTENT(IN)                                :: virial
      REAL(dp), INTENT(INOUT)                            :: derivatives(15)

      INTEGER                                            :: component, i, j
      REAL(dp)                                           :: ab(3), cd(3), d1(3), f0, f1(3), frac(4), &
                                                            moment(3, 3), rel(3), rho2, x(3)

      CPASSERT(records(1)%prepared_max_order >= 1)
      x = records(1)%P - records(1)%Q
      f0 = records(1)%Fm(1)
      f1 = x*records(1)%Fm(2)
      IF (virial) THEN
         ! Only the product-center separation varies within a displacement group.
         DO j = 1, 3
            moment(j:3, j) = f1(j:3)*x(j)
         END DO
         DO i = 2, SIZE(records)
            ASSOCIATE (p => records(i))
               x = p%P - p%Q
               d1 = x*p%Fm(2)
               f0 = f0 + p%Fm(1)
               f1 = f1 + d1
               DO j = 1, 3
                  moment(j:3, j) = moment(j:3, j) + d1(j:3)*x(j)
               END DO
            END ASSOCIATE
         END DO
      ELSE
         DO i = 2, SIZE(records)
            ASSOCIATE (p => records(i))
               d1 = (p%P - p%Q)*p%Fm(2)
               f0 = f0 + p%Fm(1)
               f1 = f1 + d1
            END ASSOCIATE
         END DO
      END IF

      ASSOCIATE (p => records(1))
         frac = p%exponents*[p%ZetaInv, p%ZetaInv, p%EtaInv, p%EtaInv]
         ab = (-2.0_dp*p%exponents(1)*frac(2))*abcd(:, 1)
         cd = (-2.0_dp*p%exponents(3)*frac(4))*abcd(:, 2)
         rho2 = -2.0_dp*p%Rho
         rel = rho2*f1
         derivatives(1:3) = derivatives(1:3) + (ab*f0 + frac(1)*rel)
         derivatives(4:6) = derivatives(4:6) + (cd*f0 - frac(3)*rel)
         derivatives(7:9) = derivatives(7:9) + (-cd*f0 - frac(4)*rel)

         IF (virial) THEN
            component = 9
            DO j = 1, 3
               DO i = j, 3
                  component = component + 1
                  derivatives(component) = derivatives(component) + &
                                           (f0*(ab(i)*abcd(j, 1) + cd(i)*abcd(j, 2)) + rho2*moment(i, j))
               END DO
            END DO
         END IF
      END ASSOCIATE
   END SUBROUTINE add_ssss_derivatives

! **************************************************************************************************
!> \brief Expand all contracted strain tensors and transform them to cell derivatives.
!> \param n number of contracted AO quartets
!> \param virial lower triangle on input, full cell derivatives on output
!> \param cell simulation cell
! **************************************************************************************************
   SUBROUTINE expand_symmetric_virial(n, virial, cell)
      INTEGER, INTENT(IN)                                :: n
      REAL(dp), INTENT(INOUT)                            :: virial(n, 3, 3)
      TYPE(cell_type), INTENT(IN)                        :: cell

      INTEGER                                            :: i, j, k
      REAL(dp)                                           :: row(3)

      IF (cell%orthorhombic) THEN
         DO j = 1, 3
            DO i = j + 1, 3
               virial(:, j, i) = virial(:, i, j)*cell%h_inv(i, i)
            END DO
            virial(:, j:3, j) = virial(:, j:3, j)*cell%h_inv(j, j)
         END DO
      ELSE
         virial(:, 1, 2) = virial(:, 2, 1)
         virial(:, 1, 3) = virial(:, 3, 1)
         virial(:, 2, 3) = virial(:, 3, 2)
         DO j = 1, 3
!$OMP SIMD PRIVATE(row)
            DO i = 1, n
               row = virial(i, j, :)
               DO k = 1, 3
                  virial(i, j, k) = SUM(row*cell%h_inv(k, :))
               END DO
            END DO
         END DO
      END IF
   END SUBROUTINE expand_symmetric_virial

! **************************************************************************************************
!> \brief Initialize angular and polynomial lookup tables for the current shell ranges.
!> \param work context-owned workspace reused between primitive quartets
!> \param lower minimum angular momentum per center
!> \param upper maximum angular momentum per center
!> \param virial whether to include the shifted radial moments
! **************************************************************************************************
   SUBROUTINE prepare_hermite_workspace(work, lower, upper, virial)
      TYPE(eri4_hermite_workspace_type), INTENT(INOUT)   :: work
      INTEGER, INTENT(IN)                                :: lower(4), upper(4)
      LOGICAL, INTENT(IN)                                :: virial

      INTEGER :: a, aa, axis, b, bb, c, cc, d, dc, degree, dims(4), dst, h, i, j, k, l, &
         lookup(0:max_order + 1, 0:max_order + 1, 0:max_order + 1), m, max_degree, min_degree, nh, &
         nslot, nt, nterms, start(4), stride(4), t(3), term_slot(max_polynomial_terms), u(3), x, &
         y, z
      REAL(dp) :: coefficient(0:max_order + 1, 0:(max_order + 1)/2)

      IF (ALL(lower == work%lower) .AND. ALL(upper == work%upper)) THEN
         IF (virial .EQV. work%virial) RETURN
      END IF
      work%lower = lower
      work%upper = upper
      work%virial = virial
      work%order = SUM(upper) + 1
      ! Degrees of the ordinary radial derivatives and the shifted auxiliaries.
      l = work%order
      work%nherm = (l + 1)*(l + 2)*(l + 3)/6
      nh = l + MERGE(1, 0, virial)
      work%nshift = (nh + 1)*(nh + 2)*(nh + 3)/6

      i = 0
      DO l = 0, nh
         DO x = l, 0, -1
            DO y = l - x, 0, -1
               z = l - x - y
               i = i + 1
               work%powers(:, i) = [x, y, z]
               lookup(x, y, z) = i
               IF (i > 1) THEN
                  t = [x, y, z]
                  axis = FINDLOC(t > 0, .TRUE., DIM=1)
                  t(axis) = t(axis) - 1
                  work%axis(i) = axis
                  work%parent(i) = lookup(t(1), t(2), t(3))
               END IF
            END DO
         END DO
         coefficient(l, 0) = 1.0_dp
         DO j = 1, l/2
            coefficient(l, j) = coefficient(l, j - 1)*REAL((l - 2*j + 2)*(l - 2*j + 1), dp)/REAL(2*j, dp)
         END DO
      END DO

      ! Each Cartesian derivative of a radial seed is a sum of monomials times radial orders.
      ! Store only the distinct (order, monomial) pairs accumulated inside the image loop.
      nterms = 0
      DO k = 0, MERGE(1, 0, virial)
         nt = MERGE(work%nshift, work%nherm, k == 1)
         DO i = k + 1, nt
            h = i + k*work%nherm
            work%term_offset(h) = nterms + 1
            t = work%powers(:, i)
            degree = SUM(t)
            DO z = 0, t(3)/2
               DO y = 0, t(2)/2
                  DO x = 0, t(1)/2
                     u = t - 2*[x, y, z]
                     j = lookup(u(1), u(2), u(3))
                     m = degree - x - y - z - k
                     nterms = nterms + 1
                     term_slot(nterms) = j + m*max_monomials
                     work%term_factor(nterms) = coefficient(t(1), x)*coefficient(t(2), y)*coefficient(t(3), z)
                  END DO
               END DO
            END DO
            work%term_offset(h + 1) = nterms + 1
         END DO
      END DO
      nslot = 0
      DO m = 0, work%order
         ! At radial order m, trace terms give all degrees from MAX(0, 2*m - order) to m.
         ! The shifted virial auxiliary also needs degree m + 1. Degrees are stored consecutively.
         min_degree = MAX(0, 2*m - work%order)
         max_degree = m + MERGE(1, 0, virial)
         work%monomial_start(m) = 1 + min_degree*(min_degree + 1)*(min_degree + 2)/6
         work%radial_offset(m) = nslot + 1
         nslot = nslot + (max_degree + 1)*(max_degree + 2)*(max_degree + 3)/6 - work%monomial_start(m) + 1
      END DO
      work%radial_offset(work%order + 1) = nslot + 1
      work%nradial = nslot
      DO i = 1, nterms
         j = term_slot(i) - 1
         m = j/max_monomials
         work%term_index(i) = work%radial_offset(m) + MOD(j, max_monomials) + 1 - work%monomial_start(m)
      END DO

      IF (virial) THEN
         DO j = 1, work%nherm
            DO axis = 1, 3
               t = work%powers(:, j)
               t(axis) = t(axis) + 1
               work%upper_index(axis, j) = lookup(t(1), t(2), t(3))
               t(axis) = MAX(0, t(axis) - 2)
               work%lower_index(axis, j) = lookup(t(1), t(2), t(3))
            END DO
         END DO
      END IF

      a = SUM(upper(1:2)) + 1
      b = SUM(upper(3:4)) + 1
      IF (virial) work%moment_count(:work%nherm) = 0
      DO j = 1, (b + 1)*(b + 2)*(b + 3)/6
         work%ket_sign(j) = REAL((-1)**SUM(work%powers(:, j)), dp)
         DO i = 1, (a + 1)*(a + 2)*(a + 3)/6
            t = work%powers(:, i) + work%powers(:, j)
            IF (SUM(t) > work%order) CYCLE
            work%pair_index(i, j) = lookup(t(1), t(2), t(3))
            IF (virial .AND. i <= a*(a + 1)*(a + 2)/6) THEN
               dst = work%pair_index(i, j)
               work%moment_count(dst) = work%moment_count(dst) + 1
            END IF
         END DO
      END DO

      dims = ncoset(upper) - ncoset(lower - 1)
      stride = [1, dims(1), dims(1)*dims(2), dims(1)*dims(2)*dims(3)]
      dst = 1
      DO d = lower(4), upper(4)
         DO c = lower(3), upper(3)
            DO b = lower(2), upper(2)
               DO a = lower(1), upper(1)
                  start = ncoset([a, b, c, d] - 1) - ncoset(lower - 1)
                  work%block_start(a + 3*b + 9*c + 27*d) = dst
                  DO dc = 0, nco(d) - 1
                     DO cc = 0, nco(c) - 1
                        DO bb = 0, nco(b) - 1
                           DO aa = 1, nco(a)
                              work%cartesian_index(dst) = SUM((start + [aa, bb, cc, dc])*stride)
                              dst = dst + 1
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE prepare_hermite_workspace

! **************************************************************************************************
!> \brief Construct Hermite expansions of the pair and its Gaussian-center derivatives.
!> \param powers Hermite multi-indices
!> \param lower minimum angular momentum at each center
!> \param upper maximum angular momentum at each center
!> \param exponents the pair's two Gaussian exponents
!> \param displacement displacement from the second center to the first
!> \param expansion pair expansion matrices with the undifferentiated matrix first
! **************************************************************************************************
   SUBROUTINE pair_expansions(powers, lower, upper, exponents, displacement, expansion)
      INTEGER, INTENT(IN)                                :: powers(3, max_monomials), lower(2), &
                                                            upper(2)
      REAL(dp), INTENT(IN)                               :: exponents(2), displacement(3)
      REAL(dp), CONTIGUOUS, INTENT(OUT)                  :: expansion(:, :, 0:)

      INTEGER                                            :: a, b, h, i, j, jcart, kcart, last(2), n, &
                                                            pair, start(2), t, top(2), x, y
      REAL(dp)                                           :: aa, bb, coeff(-1:7, -1:3, -1:3, 3), &
                                                            factors(3), prefac, product, u

      top = upper
      start = ncoset(lower - 1) + 1
      last = ncoset(top)
      aa = exponents(1)
      bb = exponents(2)
      u = 1.0_dp/(aa + bb)
      coeff = 0.0_dp
      ! The Gaussian attenuation is already included in the radial seeds.
      coeff(0, 0, 0, :) = 1.0_dp
      DO j = 1, 3
         DO b = 0, top(2) + 1
            DO a = 0, top(1) + 1
               IF (a > 0) THEN
                  prefac = -bb*u*displacement(j)
                  DO t = 0, a + b
                     coeff(t, a, b, j) = 0.5_dp*u*coeff(t - 1, a - 1, b, j) + &
                                         prefac*coeff(t, a - 1, b, j) + (t + 1)*coeff(t + 1, a - 1, b, j)
                  END DO
               ELSE IF (b > 0) THEN
                  prefac = aa*u*displacement(j)
                  DO t = 0, b
                     coeff(t, a, b, j) = 0.5_dp*u*coeff(t - 1, a, b - 1, j) + &
                                         prefac*coeff(t, a, b - 1, j) + (t + 1)*coeff(t + 1, a, b - 1, j)
                  END DO
               END IF
            END DO
         END DO
      END DO

      n = SUM(top)
      h = (n + 2)*(n + 3)*(n + 4)/6
      pair = 0
      DO kcart = start(2), last(2)
         DO jcart = start(1), last(1)
            pair = pair + 1
            DO i = 1, h
               DO j = 1, 3
                  factors(j) = coeff(powers(j, i), indco(j, jcart), indco(j, kcart), j)
               END DO
               expansion(i, pair, 0) = factors(1)*factors(2)*factors(3)
               DO j = 1, 3
                  t = powers(j, i)
                  x = indco(j, jcart)
                  y = indco(j, kcart)
                  product = factors(MOD(j, 3) + 1)*factors(MOD(j + 1, 3) + 1)
                  expansion(i, pair, j) = product*(2.0_dp*aa*coeff(t, x + 1, y, j) - x*coeff(t, x - 1, y, j))
                  IF (UBOUND(expansion, 3) == 6) THEN
                     expansion(i, pair, j + 3) = product*(2.0_dp*bb*coeff(t, x, y + 1, j) - y*coeff(t, x, y - 1, j))
                  END IF
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE pair_expansions

! **************************************************************************************************
!> \brief Accumulate all Cartesian angular blocks of one displacement group.
!> \param work prepared reusable workspace
!> \param records prepared quartets with common exponents and pair displacements
!> \param abcd fixed pair displacements
!> \param cell cell for fractional-coordinate virial weighting
!> \param derivatives Cartesian center derivatives
!> \param nc active rows in the reusable derivatives buffer
!> \param density primitive exchange density and screening budget
!> \param virial accumulated density-weighted cell derivative
! **************************************************************************************************
   SUBROUTINE add_hermite_derivatives(work, records, abcd, cell, derivatives, nc, density, virial)
      TYPE(eri4_hermite_workspace_type), INTENT(INOUT)   :: work
      TYPE(eri_primitive_data_type), INTENT(IN)          :: records(:)
      REAL(dp), INTENT(IN)                               :: abcd(3, 2)
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), CONTIGUOUS, INTENT(INOUT)                :: derivatives(:, :)
      INTEGER, INTENT(IN)                                :: nc
      TYPE(eri4_hermite_density_type), INTENT(INOUT)     :: density
      REAL(dp), INTENT(INOUT)                            :: virial(3, 3)

      INTEGER                                            :: a, b, dha, dhb, dims(4), first, ha, hb, &
                                                            i, j, k, last, m, na, nb, offset, p, t
      REAL(dp) :: bound, bra_norm, delta(3), f, f_a, f_d, f_norm(3, 2), ket_norm, moment_norm(3), &
         moments(3, max_hermites), mono(max_monomials), projected(max_hermites), &
         radial(max_radial_moments), ref_norm(3, 3), reference(3), relative(3, 2), &
         rho_power(0:max_order), shifted(max_monomials), zero(max_hermites)

      radial(:work%nradial) = 0.0_dp
      rho_power(0) = 1.0_dp
      DO m = 1, work%order
         rho_power(m) = (-2.0_dp*records(1)%Rho)*rho_power(m - 1)
      END DO
      DO i = 1, SIZE(records)
         delta = records(i)%P - records(i)%Q
         mono(1) = 1.0_dp
         DO j = 2, work%nshift
            mono(j) = delta(work%axis(j))*mono(work%parent(j))
         END DO
         DO m = 0, work%order
            f = rho_power(m)*records(i)%Fm(m + 1)
            first = work%radial_offset(m)
            last = work%radial_offset(m + 1) - first - 1
            offset = work%monomial_start(m)
            DO j = 0, last
               radial(first + j) = radial(first + j) + f*mono(offset + j)
            END DO
         END DO
      END DO

      DO j = 1, work%nherm
         zero(j) = 0.0_dp
         DO t = work%term_offset(j), work%term_offset(j + 1) - 1
            zero(j) = zero(j) + work%term_factor(t)*radial(work%term_index(t))
         END DO
      END DO
      IF (work%virial) THEN
         DO j = 2, work%nshift
            p = j + work%nherm
            shifted(j) = 0.0_dp
            DO t = work%term_offset(p), work%term_offset(p + 1) - 1
               shifted(j) = shifted(j) + work%term_factor(t)*radial(work%term_index(t))
            END DO
         END DO
         ! Weight the ket translation by C-B for every image, including the reference image.
         reference = (records(1)%P - records(1)%Q) + (records(1)%centers(:, 3) - records(1)%centers(:, 2))
         DO j = 1, work%nherm
            ! Multiplying a Hermite derivative by P-Q raises its auxiliary order
            ! and subtracts the derivative of the coordinate factor.
            delta = reference*zero(j) + &
                    work%powers(:, j)*zero(work%lower_index(:, j)) - shifted(work%upper_index(:, j))
            CALL real_to_scaled(moments(:, j), delta, cell)
         END DO
         CALL real_to_scaled(relative(:, 1), abcd(:, 1), cell)
         CALL real_to_scaled(relative(:, 2), -abcd(:, 2), cell)
      END IF

      CALL pair_expansions(work%powers, work%lower(1:2), work%upper(1:2), records(1)%exponents(1:2), abcd(:, 1), work%bra)
      CALL pair_expansions(work%powers, work%lower(3:4), work%upper(3:4), records(1)%exponents(3:4), abcd(:, 2), work%ket)
      a = SUM(work%upper(1:2))
      b = SUM(work%upper(3:4))
      ha = (a + 1)*(a + 2)*(a + 3)/6
      hb = (b + 1)*(b + 2)*(b + 3)/6
      dha = (a + 2)*(a + 3)*(a + 4)/6
      dhb = (b + 2)*(b + 3)*(b + 4)/6
      dims = ncoset(work%upper) - ncoset(work%lower - 1)
      na = dims(1)*dims(2)
      nb = dims(3)*dims(4)

      DO j = 1, hb
         DO i = 1, dha
            work%auxiliary(i, j) = work%ket_sign(j)*zero(work%pair_index(i, j))
         END DO
      END DO
      CALL DGEMM('N', 'N', dha, nb, hb, 1.0_dp, work%auxiliary, max_pair_hermites, work%ket(1, 1, 0), &
                 max_pair_hermites, 0.0_dp, work%intermediate, dha)
      DO k = 1, 3
         CALL DGEMM('T', 'N', na, nb, dha, 1.0_dp, work%bra(1, 1, k), max_pair_hermites, work%intermediate, &
                    dha, 0.0_dp, work%force(1, k), na)
      END DO

      DO j = 1, dhb
         DO i = 1, ha
            work%auxiliary(i, j) = work%ket_sign(j)*zero(work%pair_index(i, j))
         END DO
      END DO
      CALL DGEMM('T', 'N', na, dhb, ha, 1.0_dp, work%bra(1, 1, 0), max_pair_hermites, work%auxiliary, &
                 max_pair_hermites, 0.0_dp, work%intermediate, na)
      DO k = 1, 6
         CALL DGEMM('N', 'N', na, nb, dhb, 1.0_dp, work%intermediate, na, work%ket(1, 1, k), &
                    max_pair_hermites, 0.0_dp, work%force(1, 6 + k), na)
      END DO
      DO k = 1, 3
         work%force(:nc, k + 3) = -(work%force(:nc, k) + work%force(:nc, k + 6) + work%force(:nc, k + 9))
      END DO
      derivatives(:nc, :12) = derivatives(:nc, :12) + work%force(:nc, :)

      IF (.NOT. work%virial) RETURN
      IF (density%norm == 0.0_dp) RETURN

      IF (density%remaining > 0.0_dp) THEN
         DO k = 1, 3
            f_norm(k, 1) = NORM2(work%force(:nc, k))
            f_norm(k, 2) = NORM2(work%force(:nc, k + 9))
         END DO
         ref_norm = MATMUL(f_norm, TRANSPOSE(ABS(relative)))
         moment_norm = 0.0_dp
         DO j = 1, work%nherm
            moment_norm = moment_norm + work%moment_count(j)*moments(:, j)**2
         END DO
         moment_norm = SQRT(moment_norm)
         bra_norm = NORM2(work%bra(:ha, :na, 0))
         ! A translation shifts the ket expansion into distinct rows and preserves its norm.
         ket_norm = NORM2(work%ket(:hb, :nb, 0))
         DO m = 1, 3
            ref_norm(:, m) = ref_norm(:, m) + (bra_norm*ket_norm)*moment_norm(m)
         END DO
         bound = MAXVAL(ref_norm)
         ! Include roundoff in the componentwise bound before charging it to the budget.
         bound = bound*density%norm*(1.0_dp + 1.0E-12_dp)
         IF (bound <= density%remaining) THEN
            density%remaining = density%remaining - bound
            RETURN
         END IF
      END IF
      IF (.NOT. density%weights_ready) CALL form_hermite_weights(density, dims)

      DO k = 1, 3
         f_a = DOT_PRODUCT(density%weight(:nc), work%force(:nc, k))
         f_d = DOT_PRODUCT(density%weight(:nc), work%force(:nc, k + 9))
         virial(k, :) = virial(k, :) + relative(:, 1)*f_a + relative(:, 2)*f_d
      END DO

      ! Both ket centers translate together, so the ket derivative is an index shift.
      CALL DGEMM('N', 'N', ha, nb, na, 1.0_dp, work%bra(1, 1, 0), max_pair_hermites, &
                 density%weight, na, 0.0_dp, work%intermediate, ha)
      CALL DGEMM('N', 'T', ha, hb, nb, 1.0_dp, work%intermediate, ha, work%ket(1, 1, 0), &
                 max_pair_hermites, 0.0_dp, work%auxiliary, max_pair_hermites)

      ! All entries with the same total Hermite index use the same translated moment.
      ! Raising the ket index flips its parity.
      p = work%order*(work%order + 1)*(work%order + 2)/6
      projected(:p) = 0.0_dp
      DO j = 1, hb
         DO i = 1, ha
            t = work%pair_index(i, j)
            projected(t) = projected(t) - work%ket_sign(j)*work%auxiliary(i, j)
         END DO
      END DO
      DO k = 1, 3
         DO j = 1, p
            virial(k, :) = virial(k, :) + projected(j)*moments(:, work%upper_index(k, j))
         END DO
      END DO
   END SUBROUTINE add_hermite_derivatives

! **************************************************************************************************
!> \brief Extract one angular block and its screening maxima from the group sums.
!> \param batch the accumulated derivative blocks
!> \param l angular momenta in original center order
!> \param forces force block for the angular quartet
!> \param force_max force maximum for Cartesian screening
! **************************************************************************************************
   SUBROUTINE extract_hermite_derivatives(batch, l, forces, force_max)
      TYPE(eri4_hermite_batch_type), INTENT(IN)          :: batch
      INTEGER, INTENT(IN)                                :: l(4)
      REAL(dp), INTENT(OUT)                              :: forces(:, :), force_max

      INTEGER                                            :: i, k, ncart, start
      REAL(dp)                                           :: value

      IF (ALL(l == 0)) THEN
         forces(1, :) = batch%cartesian(1, :12)
         force_max = MAXVAL(ABS(forces))
         RETURN
      END IF

      ncart = SIZE(forces, 1)
      start = batch%workspace%block_start(l(1) + 3*l(2) + 9*l(3) + 27*l(4)) - 1
      force_max = 0.0_dp
      DO k = 1, 12
!$OMP SIMD REDUCTION(MAX:force_max)
         DO i = 1, ncart
            value = batch%cartesian(batch%workspace%cartesian_index(start + i), k)
            forces(i, k) = value
            force_max = MAX(force_max, ABS(value))
         END DO
      END DO
   END SUBROUTINE extract_hermite_derivatives

! **************************************************************************************************
!> \brief Backtransform each AO density block across the full shell-set angular range.
!> \param density transformed densities and exchange product terms
!> \param ao AO density in BD, BC, AD and AC order
!> \param factor exchange coefficient including quartet symmetry
!> \param lower first angular momentum of each shell set
!> \param upper last angular momentum of each shell set
!> \param np number of primitives per center
!> \param ns number of spherical functions per center
!> \param ca contraction data for A
!> \param cb contraction data for B
!> \param cc contraction data for C
!> \param cd contraction data for D
! **************************************************************************************************
   SUBROUTINE prepare_hermite_density(density, ao, factor, lower, upper, np, ns, ca, cb, cc, cd)
      TYPE(eri4_hermite_density_type), INTENT(OUT)       :: density
      TYPE(eri4_density_type), INTENT(IN)                :: ao
      REAL(dp), INTENT(IN)                               :: factor
      INTEGER, INTENT(IN)                                :: lower(4), upper(4), np(4), ns(4)
      TYPE(eri4_contraction_type), INTENT(IN), TARGET    :: ca(0:), cb(0:), cc(0:), cd(0:)

      INTEGER, PARAMETER                                 :: left(4) = [2, 2, 1, 1], &
                                                            right(4) = [4, 3, 4, 3]

      INTEGER                                            :: a, b, block, ip, jp, nc(4), ncol, ndens, &
                                                            nr, spin, t
      LOGICAL                                            :: all_s
      TYPE(eri4_density_basis_type)                      :: basis(4)

      all_s = ALL(upper == 0)
      nc = ncoset(upper) - ncoset(lower - 1)
      ndens = SIZE(ao%blocks, 3)
      IF (all_s) THEN
         ! The angular transform of a normalized s orbital is the identity.
         basis(1)%radial => ca(0)%weights
         basis(2)%radial => cb(0)%weights
         basis(3)%radial => cc(0)%weights
         basis(4)%radial => cd(0)%weights
         ALLOCATE (density%ssss(MAXVAL(np), MAXVAL(np), 4, ndens))
      ELSE
         ALLOCATE (basis(1)%values(nc(1)*np(1), ns(1)), basis(2)%values(nc(2)*np(2), ns(2)), &
                   basis(3)%values(nc(3)*np(3), ns(3)), basis(4)%values(nc(4)*np(4), ns(4)))
         CALL pack_hermite_basis(lower(1), upper(1), np(1), ca, basis(1)%values)
         CALL pack_hermite_basis(lower(2), upper(2), np(2), cb, basis(2)%values)
         CALL pack_hermite_basis(lower(3), upper(3), np(3), cc, basis(3)%values)
         CALL pack_hermite_basis(lower(4), upper(4), np(4), cd, basis(4)%values)
      END IF

      density%factor = factor
      DO spin = 1, ao%nspins
         t = density%nterms
         IF (.NOT. ao%response_only .OR. ndens == ao%nspins) THEN
            t = t + 1
            density%lhs(t) = spin
            density%rhs(t) = spin
         END IF
         IF (ndens > ao%nspins) THEN
            density%lhs(t + 1:t + 2) = [spin, spin + ao%nspins]
            density%rhs(t + 1:t + 2) = [spin + ao%nspins, spin]
            t = t + 2
         END IF
         density%nterms = t
      END DO

      DO block = 1, 4
         a = left(block)
         b = right(block)
         IF (all_s) THEN
            DO t = 1, ndens
               CALL backtransform_ssss_density(ns(a), ns(b), np(a), np(b), basis(a)%radial, basis(b)%radial, &
                                               ao%blocks(:, block, t), density%ssss(:, :, block, t))
            END DO
            CYCLE
         END IF
         nr = nc(a)*np(a)
         ncol = nc(b)*np(b)
         ASSOCIATE (pair => density%pair(block))
            ALLOCATE (pair%values(nr, ncol, ndens), pair%norms(np(a), np(b), ndens))
            DO t = 1, ndens
               pair%values(:, :, t) = MATMUL(MATMUL(basis(a)%values, &
                                                    RESHAPE(ao%blocks(:ns(a)*ns(b), block, t), [ns(a), ns(b)])), &
                                             TRANSPOSE(basis(b)%values))
               DO jp = 1, np(b)
                  DO ip = 1, np(a)
                     pair%norms(ip, jp, t) = NORM2(pair%values((ip - 1)*nc(a) + 1:ip*nc(a), &
                                                               (jp - 1)*nc(b) + 1:jp*nc(b), t))
                  END DO
               END DO
            END DO
         END ASSOCIATE
      END DO
      IF (.NOT. all_s) ALLOCATE (density%weight(PRODUCT(nc)))
   END SUBROUTINE prepare_hermite_density

! **************************************************************************************************
!> \brief Pack the separable radial and angular coefficients in primitive Cartesian by AO order.
!> \param lower minimum angular momentum
!> \param upper maximum angular momentum
!> \param np number of primitives
!> \param contraction angular and radial coefficients
!> \param basis packed transform
! **************************************************************************************************
   SUBROUTINE pack_hermite_basis(lower, upper, np, contraction, basis)
      INTEGER, INTENT(IN)                                :: lower, upper, np
      TYPE(eri4_contraction_type), INTENT(IN)            :: contraction(0:)
      REAL(dp), INTENT(OUT)                              :: basis(:, :)

      INTEGER                                            :: cfirst, l, nc, ncart, ns, p, q, sfirst

      ncart = ncoset(upper) - ncoset(lower - 1)
      basis = 0.0_dp
      sfirst = 1
      DO l = lower, upper
         nc = nco(l)
         ns = nso(l)
         ASSOCIATE (shell => contraction(l))
            DO q = 1, SIZE(shell%weights, 1)
               DO p = 1, np
                  cfirst = (p - 1)*ncart + ncoset(l - 1) - ncoset(lower - 1) + 1
                  basis(cfirst:cfirst + nc - 1, sfirst:sfirst + ns - 1) = shell%weights(q, p)*shell%angular
               END DO
               sfirst = sfirst + ns
            END DO
         END ASSOCIATE
      END DO
   END SUBROUTINE pack_hermite_basis

! **************************************************************************************************
!> \brief Select primitive density subblocks and bound the exchange weight norm.
!> \param density prepared shell-set density
!> \param primitive indices of the four primitives
!> \param ncart Cartesian size of each primitive
!> \param budget maximum omitted cell derivative from this primitive quartet
! **************************************************************************************************
   SUBROUTINE select_hermite_density(density, primitive, ncart, budget)
      TYPE(eri4_hermite_density_type), INTENT(INOUT)     :: density
      INTEGER, INTENT(IN)                                :: primitive(4), ncart(4)
      REAL(dp), INTENT(IN)                               :: budget

      INTEGER                                            :: a, b, c, d, r, s, t

      a = primitive(1)
      b = primitive(2)
      c = primitive(3)
      d = primitive(4)
      density%offset = (primitive - 1)*ncart
      density%remaining = budget
      density%weights_ready = .FALSE.
      density%norm = 0.0_dp
      DO t = 1, density%nterms
         r = density%lhs(t)
         s = density%rhs(t)
         density%norm = density%norm + density%pair(4)%norms(a, c, r)*density%pair(1)%norms(b, d, s) + &
                        density%pair(3)%norms(a, d, r)*density%pair(2)%norms(b, c, s)
      END DO
      density%norm = ABS(density%factor)*density%norm
   END SUBROUTINE select_hermite_density

! **************************************************************************************************
!> \brief Build the exchange weights only after a displacement group passes screening.
!> \param density selected primitive density
!> \param dims Cartesian block sizes
! **************************************************************************************************
   SUBROUTINE form_hermite_weights(density, dims)
      TYPE(eri4_hermite_density_type), INTENT(INOUT)     :: density
      INTEGER, INTENT(IN)                                :: dims(4)

      INTEGER                                            :: a, b, c, d, first(4), i, last(4), r, s, t
      REAL(dp)                                           :: weight

      first = density%offset + 1
      last = density%offset + dims
      i = 0
      ASSOCIATE (bd => density%pair(1)%values, bc => density%pair(2)%values, &
                 ad => density%pair(3)%values, ac => density%pair(4)%values)
         DO d = first(4), last(4)
            DO c = first(3), last(3)
               DO b = first(2), last(2)
                  DO a = first(1), last(1)
                     weight = 0.0_dp
                     DO t = 1, density%nterms
                        r = density%lhs(t)
                        s = density%rhs(t)
                        weight = weight + ac(a, c, r)*bd(b, d, s) + ad(a, d, r)*bc(b, c, s)
                     END DO
                     i = i + 1
                     density%weight(i) = density%factor*weight
                  END DO
               END DO
            END DO
         END DO
      END ASSOCIATE
      density%weights_ready = .TRUE.
   END SUBROUTINE form_hermite_weights

! **************************************************************************************************
!> \brief Backtransform a two-center s-shell density block using radial coefficients.
!> \param na left AO count
!> \param nb right AO count
!> \param np left primitive count
!> \param nq right primitive count
!> \param wa left radial coefficients
!> \param wb right radial coefficients
!> \param ao AO density block
!> \param primitive resulting primitive density block
! **************************************************************************************************
   SUBROUTINE backtransform_ssss_density(na, nb, np, nq, wa, wb, ao, primitive)
      INTEGER, INTENT(IN)                                :: na, nb, np, nq
      REAL(dp), INTENT(IN)                               :: wa(na, np), wb(nb, nq), ao(na, nb)
      REAL(dp), INTENT(OUT)                              :: primitive(:, :)

      INTEGER                                            :: a, b, i, j
      REAL(dp)                                           :: intermediate(np, nb)

      DO b = 1, nb
         DO i = 1, np
            intermediate(i, b) = DOT_PRODUCT(wa(:, i), ao(:, b))
         END DO
      END DO
      DO j = 1, nq
         primitive(:np, j) = intermediate(:, 1)*wb(1, j)
         DO a = 2, nb
            primitive(:np, j) = primitive(:np, j) + intermediate(:, a)*wb(a, j)
         END DO
      END DO
   END SUBROUTINE backtransform_ssss_density

! **************************************************************************************************
!> \brief Exchange weight of a single all-s primitive quartet.
!> \param density backtransformed shell-set density
!> \param a primitive index at center A
!> \param b primitive index at center B
!> \param c primitive index at center C
!> \param d primitive index at center D
!> \return exchange weight including the exchange coefficient
! **************************************************************************************************
   PURE FUNCTION ssss_density_weight(density, a, b, c, d) RESULT(weight)
      TYPE(eri4_hermite_density_type), INTENT(IN)        :: density
      INTEGER, INTENT(IN)                                :: a, b, c, d
      REAL(dp)                                           :: weight

      INTEGER                                            :: r, s, t

      weight = 0.0_dp
      ASSOCIATE (p => density%ssss)
         DO t = 1, density%nterms
            r = density%lhs(t)
            s = density%rhs(t)
            weight = weight + p(a, c, 4, r)*p(b, d, 1, s) + p(a, d, 3, r)*p(b, c, 2, s)
         END DO
      END ASSOCIATE
      weight = weight*density%factor
   END FUNCTION ssss_density_weight

END MODULE four_center_integral_derivatives
