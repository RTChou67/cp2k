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
   USE kinds,                           ONLY: dp
   USE libint_wrapper,                  ONLY: eri_primitive_data_type
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: add_ssss_derivatives, add_hermite_derivatives, prepare_hermite_workspace, extract_hermite_derivatives

   ! Polynomial degrees for a first derivative of (dd|dd) and its coordinate moment.
   INTEGER, PARAMETER :: max_order = 9, max_monomials = 286, max_hermites = 220, max_pair_hermites = 56
   INTEGER, PARAMETER :: max_cartesian = 10**4, max_pair_cartesian = 10**2
   INTEGER, PARAMETER :: max_radial_moments = (max_order + 1)*max_monomials
   INTEGER, PARAMETER :: max_polynomial_terms = 2*max_monomials*(max_order/2 + 2)

!> \brief Reusable polynomial, basis transform and gather data for one angular range.
   TYPE, PUBLIC :: eri4_hermite_workspace_type
      PRIVATE
      INTEGER :: lower(4) = -1, upper(4) = -1, order = 0, nherm = 0, nshift = 0, nradial = 0
      LOGICAL :: virial = .FALSE.
      INTEGER :: powers(3, max_monomials) = 0
      INTEGER :: parent(max_monomials) = 0, axis(max_monomials) = 0
      INTEGER :: monomial(max_radial_moments) = 0, radial_offset(0:max_order + 1) = 0
      INTEGER :: term_offset(max_hermites + max_monomials + 1) = 0
      INTEGER :: term_index(max_polynomial_terms) = 0
      REAL(dp) :: term_factor(max_polynomial_terms) = 0.0_dp
      INTEGER :: lower_index(3, max_hermites) = 0, upper_index(3, max_hermites) = 0
      INTEGER :: pair_index(max_pair_hermites, max_pair_hermites) = 0
      REAL(dp) :: ket_sign(max_pair_hermites) = 0.0_dp
      INTEGER :: block_start(0:80) = 0, cartesian_index(max_cartesian) = 0
      REAL(dp) :: bra(max_pair_hermites, max_pair_cartesian, 0:3) = 0.0_dp
      REAL(dp) :: ket(max_pair_hermites, max_pair_cartesian, 0:6) = 0.0_dp
      REAL(dp) :: auxiliary(max_pair_hermites, max_pair_hermites) = 0.0_dp
      REAL(dp) :: intermediate(max_pair_hermites*max_pair_cartesian) = 0.0_dp
      REAL(dp) :: force(max_cartesian, 12) = 0.0_dp
   END TYPE eri4_hermite_workspace_type

CONTAINS

! **************************************************************************************************
!> \brief Sum one displacement group of all-s derivatives in Cartesian center order.
!> \param records prepared quartets with common exponents and internal pair displacements
!> \param abcd fixed internal bra and ket displacements
!> \param cell cell for fractional-coordinate virial weighting
!> \param virial whether to accumulate virial components
!> \param derivatives twelve Cartesian derivatives followed by nine virial components
! **************************************************************************************************
   SUBROUTINE add_ssss_derivatives(records, abcd, cell, virial, derivatives)
      TYPE(eri_primitive_data_type), INTENT(IN)          :: records(:)
      REAL(dp), INTENT(IN)                               :: abcd(3, 2)
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: virial
      REAL(dp), CONTIGUOUS, INTENT(INOUT)                :: derivatives(:, :)

      INTEGER                                            :: first, first_order, i, j, k
      REAL(dp)                                           :: ab(3), cd(3), d1(3), dg(3, 4), f0, &
                                                            f1(3), frac(4), moment(3, 3), rel(3), &
                                                            x(3)

      first_order = MIN(1, records(1)%prepared_max_order) + 1
      CPASSERT(first_order == 2)
      x = records(1)%P - records(1)%Q
      f0 = records(1)%Fm(1)
      f1 = x*records(1)%Fm(first_order)
      IF (virial) THEN
         ! Only the product-center separation varies within a displacement group.
         DO j = 1, 3
            moment(j:3, j) = f1(j:3)*x(j)
         END DO
         DO i = 2, SIZE(records)
            ASSOCIATE (p => records(i))
               x = p%P - p%Q
               d1 = x*p%Fm(first_order)
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
               d1 = (p%P - p%Q)*p%Fm(first_order)
               f0 = f0 + p%Fm(1)
               f1 = f1 + d1
            END ASSOCIATE
         END DO
      END IF

      ASSOCIATE (p => records(1))
         ! Differentiate the pair attenuation and the product-center separation in the ssss seed.
         frac = p%exponents*[p%ZetaInv, p%ZetaInv, p%EtaInv, p%EtaInv]
         ab = (-2.0_dp*p%exponents(1)*frac(2))*abcd(:, 1)
         cd = (-2.0_dp*p%exponents(3)*frac(4))*abcd(:, 2)
         rel = (-2.0_dp*p%Rho)*f1
         dg(:, 1) = ab*f0 + frac(1)*rel
         dg(:, 2) = -ab*f0 + frac(2)*rel
         dg(:, 3) = cd*f0 - frac(3)*rel
         dg(:, 4) = -cd*f0 - frac(4)*rel
         DO k = 1, 4
            first = 3*k - 2
            derivatives(1, first:first + 2) = derivatives(1, first:first + 2) + dg(:, k)
         END DO

         IF (virial) THEN
            moment(1, 2) = moment(2, 1)
            moment(1, 3) = moment(3, 1)
            moment(2, 3) = moment(3, 2)
            DO j = 1, 3
               moment(:, j) = f0*(ab*abcd(j, 1) + cd*abcd(j, 2)) - 2.0_dp*p%Rho*moment(:, j)
            END DO
            IF (cell%orthorhombic) THEN
               DO j = 1, 3
                  moment(:, j) = moment(:, j)*cell%h_inv(j, j)
               END DO
            ELSE
               moment = MATMUL(moment, TRANSPOSE(cell%h_inv))
            END IF
            DO j = 1, 3
               first = 13 + 3*(j - 1)
               derivatives(1, first:first + 2) = derivatives(1, first:first + 2) + moment(:, j)
            END DO
         END IF
      END ASSOCIATE
   END SUBROUTINE add_ssss_derivatives

! **************************************************************************************************
!> \brief Initialize angular and polynomial lookup tables for the current shell ranges.
!> \param work context-owned workspace reused between primitive quartets
!> \param lower,upper minimum and maximum angular momentum per center
!> \param virial whether to include the shifted radial moments
! **************************************************************************************************
   SUBROUTINE prepare_hermite_workspace(work, lower, upper, virial)
      TYPE(eri4_hermite_workspace_type), INTENT(INOUT)   :: work
      INTEGER, INTENT(IN)                                :: lower(4), upper(4)
      LOGICAL, INTENT(IN)                                :: virial

      INTEGER :: a, aa, axis, b, bb, c, cc, d, dc, degree, dims(4), dst, h, i, j, k, l, &
         lookup(0:max_order + 1, 0:max_order + 1, 0:max_order + 1), m, nh, nslot, nt, nterms, &
         start(4), stride(4), t(3), term_slot(max_polynomial_terms), u(3), &
         usage(max_monomials, 0:max_order), x, y, z
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
      usage = 0
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
                     usage(j, m) = 1
                  END DO
               END DO
            END DO
            work%term_offset(h + 1) = nterms + 1
         END DO
      END DO
      nslot = 0
      DO m = 0, work%order
         work%radial_offset(m) = nslot + 1
         DO j = 1, work%nshift
            IF (usage(j, m) == 0) CYCLE
            nslot = nslot + 1
            usage(j, m) = nslot
            work%monomial(nslot) = j
         END DO
      END DO
      work%radial_offset(work%order + 1) = nslot + 1
      work%nradial = nslot
      DO i = 1, nterms
         j = term_slot(i) - 1
         work%term_index(i) = usage(MOD(j, max_monomials) + 1, j/max_monomials)
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
      DO j = 1, (b + 1)*(b + 2)*(b + 3)/6
         work%ket_sign(j) = REAL((-1)**SUM(work%powers(:, j)), dp)
         DO i = 1, (a + 1)*(a + 2)*(a + 3)/6
            t = work%powers(:, i) + work%powers(:, j)
            IF (SUM(t) > work%order) CYCLE
            work%pair_index(i, j) = lookup(t(1), t(2), t(3))
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
!> \param lower,upper angular range for the two centers
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
!> \param derivatives twelve Cartesian derivatives followed by nine virial components
!> \param nc active rows in the reusable derivatives buffer
! **************************************************************************************************
   SUBROUTINE add_hermite_derivatives(work, records, abcd, cell, derivatives, nc)
      TYPE(eri4_hermite_workspace_type), INTENT(INOUT)   :: work
      TYPE(eri_primitive_data_type), INTENT(IN)          :: records(:)
      REAL(dp), INTENT(IN)                               :: abcd(3, 2)
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), CONTIGUOUS, INTENT(INOUT)                :: derivatives(:, :)
      INTEGER, INTENT(IN)                                :: nc

      INTEGER                                            :: a, b, dha, dhb, dims(4), first, ha, hb, &
                                                            i, j, k, last, m, na, nb, p, plane, t
      REAL(dp) :: delta(3), f, moments(3, max_hermites), mono(max_monomials), &
         radial(max_radial_moments), relative(3, 3), rho_power(0:max_order), &
         shifted(max_monomials), zero(max_hermites)

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
            DO j = work%radial_offset(m), work%radial_offset(m + 1) - 1
               radial(j) = radial(j) + f*mono(work%monomial(j))
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
         DO j = 1, work%nherm
            ! R = delta(first) - delta(image) for translations relative to the first retained image.
            ! The product of delta and a derivative is an order-raised shifted auxiliary
            ! minus the derivative of the coordinate factor. No additional radial seed is needed.
            delta = (records(1)%P - records(1)%Q)*zero(j) + &
                    work%powers(:, j)*zero(work%lower_index(:, j)) - shifted(work%upper_index(:, j))
            CALL real_to_scaled(moments(:, j), delta, cell)
         END DO
         CALL real_to_scaled(relative(:, 1), abcd(:, 1), cell)
         CALL real_to_scaled(relative(:, 2), -abcd(:, 2), cell)
         CALL real_to_scaled(relative(:, 3), records(1)%centers(:, 3) - records(1)%centers(:, 2), cell)
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

      DO plane = 0, MERGE(3, 0, work%virial)
         DO j = 1, dhb
            IF (plane == 0) THEN
               DO i = 1, ha
                  work%auxiliary(i, j) = work%ket_sign(j)*zero(work%pair_index(i, j))
               END DO
            ELSE
               DO i = 1, ha
                  work%auxiliary(i, j) = work%ket_sign(j)*moments(plane, work%pair_index(i, j))
               END DO
            END IF
         END DO
         CALL DGEMM('T', 'N', na, dhb, ha, 1.0_dp, work%bra(1, 1, 0), max_pair_hermites, work%auxiliary, &
                    max_pair_hermites, 0.0_dp, work%intermediate, na)
         IF (plane == 0) THEN
            DO k = 1, 6
               CALL DGEMM('N', 'N', na, nb, dhb, 1.0_dp, work%intermediate, na, work%ket(1, 1, k), &
                          max_pair_hermites, 0.0_dp, work%force(1, 6 + k), na)
            END DO
            DO k = 1, 3
               work%force(:nc, k + 3) = -(work%force(:nc, k) + work%force(:nc, k + 6) + work%force(:nc, k + 9))
            END DO
            derivatives(:nc, :12) = derivatives(:nc, :12) + work%force(:nc, :)
            IF (work%virial) THEN
               ! With B as the origin, the reference-cell virial uses A-B, D-C and C-B.
               DO m = 1, 3
                  DO k = 1, 3
                     first = 12 + 3*(m - 1) + k
                     derivatives(:nc, first) = derivatives(:nc, first) + relative(m, 1)*work%force(:nc, k) + &
                                               relative(m, 2)*work%force(:nc, k + 9) + &
                                               relative(m, 3)*(work%force(:nc, k + 6) + work%force(:nc, k + 9))
                  END DO
               END DO
               ! The ket pair translates together. Contract the C+D derivative against each image moment.
               work%ket(:dhb, :nb, 1:3) = work%ket(:dhb, :nb, 1:3) + work%ket(:dhb, :nb, 4:6)
            END IF
         ELSE
            DO k = 1, 3
               last = 12 + 3*(plane - 1) + k
               CALL DGEMM('N', 'N', na, nb, dhb, 1.0_dp, work%intermediate, na, work%ket(1, 1, k), &
                          max_pair_hermites, 1.0_dp, derivatives(:nc, last), na)
            END DO
         END IF
      END DO
   END SUBROUTINE add_hermite_derivatives

! **************************************************************************************************
!> \brief Extract one angular block and its screening maxima from the group sums.
!> \param work angular lookup tables
!> \param l angular momenta in original center order
!> \param derivatives twelve Cartesian derivatives followed by nine virial components
!> \param virial whether virial components are present
!> \param forces,stress destination blocks for the angular quartet
!> \param force_max,stress_max maxima for Cartesian screening
! **************************************************************************************************
   SUBROUTINE extract_hermite_derivatives(work, l, derivatives, virial, forces, stress, force_max, stress_max)
      TYPE(eri4_hermite_workspace_type), INTENT(IN)      :: work
      INTEGER, INTENT(IN)                                :: l(4)
      REAL(dp), INTENT(IN)                               :: derivatives(:, :)
      LOGICAL, INTENT(IN)                                :: virial
      REAL(dp), INTENT(OUT)                              :: forces(:, :), stress(:, :, :), &
                                                            force_max, stress_max

      INTEGER                                            :: i, j, k, ncart, start
      REAL(dp)                                           :: value

      ncart = SIZE(forces, 1)
      start = work%block_start(l(1) + 3*l(2) + 9*l(3) + 27*l(4)) - 1
      force_max = 0.0_dp
      stress_max = 0.0_dp
      DO k = 1, 12
!$OMP SIMD REDUCTION(MAX:force_max)
         DO i = 1, ncart
            value = derivatives(work%cartesian_index(start + i), k)
            forces(i, k) = value
            force_max = MAX(force_max, ABS(value))
         END DO
      END DO
      IF (virial) THEN
         DO j = 1, 3
            DO k = 1, 3
!$OMP SIMD REDUCTION(MAX:stress_max)
               DO i = 1, ncart
                  value = derivatives(work%cartesian_index(start + i), 12 + 3*(j - 1) + k)
                  stress(i, k, j) = value
                  stress_max = MAX(stress_max, ABS(value))
               END DO
            END DO
         END DO
      END IF
   END SUBROUTINE extract_hermite_derivatives

END MODULE four_center_integral_derivatives
