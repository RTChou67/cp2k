!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Read-only implicit-Gamma data view for periodic ECD.
! **************************************************************************************************
MODULE qs_ecd_gamma_data
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_ecd_subspace,                 ONLY: ecd_band_manifold_type,&
                                              ecd_band_manifolds_build,&
                                              ecd_band_manifolds_release
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ecd_gamma_data'

   TYPE, PUBLIC :: ecd_gamma_spin_type
      INTEGER                                      :: nao = 0, nmo = 0
      REAL(KIND=dp), DIMENSION(:), POINTER        :: eigenvalues => NULL(), occupations => NULL()
      TYPE(cp_fm_type), POINTER                    :: mo_coeff => NULL()
      TYPE(ecd_band_manifold_type), DIMENSION(:), &
         POINTER                                   :: manifolds => NULL()
   END TYPE ecd_gamma_spin_type

   TYPE, PUBLIC :: ecd_gamma_context_type
      INTEGER                                                  :: nspin = 0, spatial_order = 0
      REAL(KIND=dp), DIMENSION(3)                            :: direction = 0.0_dp, ref_point = 0.0_dp
      TYPE(cell_type), POINTER                                :: cell => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: matrix_ks => NULL(), matrix_s => NULL()
      TYPE(dft_control_type), POINTER                         :: dft_control => NULL()
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: sab_orb => NULL(), sab_all => NULL(), &
                                                                  sap_ppnl => NULL()
      TYPE(particle_type), DIMENSION(:), POINTER              :: particle_set => NULL()
      TYPE(qs_kind_type), DIMENSION(:), POINTER               :: qs_kind_set => NULL()
      TYPE(ecd_gamma_spin_type), ALLOCATABLE, DIMENSION(:)    :: spin
   END TYPE ecd_gamma_context_type

   PUBLIC :: ecd_gamma_context_create_implicit, &
             ecd_gamma_context_release

CONTAINS

! **************************************************************************************************
!> \brief Borrow field-free implicit-Gamma data from a converged Quickstep environment.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_create_implicit(context, qs_env, direction, ref_point, spatial_order, &
                                                eps_degeneracy)
      TYPE(ecd_gamma_context_type), INTENT(OUT)         :: context
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: direction, ref_point
      INTEGER, INTENT(IN)                                :: spatial_order
      REAL(KIND=dp), INTENT(IN)                          :: eps_degeneracy

      INTEGER                                            :: ispin, nao, nmo, nrow, ncol
      LOGICAL                                            :: do_kpoints
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupations
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      NULLIFY (eigenvalues, mo_coeff, mos, occupations)
      CALL get_qs_env(qs_env, cell=context%cell, do_kpoints=do_kpoints, mos=mos, &
                      matrix_ks=context%matrix_ks, matrix_s=context%matrix_s, &
                      dft_control=context%dft_control, sab_orb=context%sab_orb, sab_all=context%sab_all, &
                      sap_ppnl=context%sap_ppnl, particle_set=context%particle_set, &
                      qs_kind_set=context%qs_kind_set)

      IF (do_kpoints) THEN
         CPABORT("ECD| Implicit-Gamma data adapter requires the field-free Gamma SCF layout")
      END IF
      IF (.NOT. ASSOCIATED(context%cell) .OR. .NOT. ASSOCIATED(mos) .OR. &
          .NOT. ASSOCIATED(context%matrix_ks) .OR. .NOT. ASSOCIATED(context%matrix_s) .OR. &
          .NOT. ASSOCIATED(context%dft_control) .OR. .NOT. ASSOCIATED(context%sab_orb) .OR. &
          .NOT. ASSOCIATED(context%sab_all) .OR. .NOT. ASSOCIATED(context%particle_set) .OR. &
          .NOT. ASSOCIATED(context%qs_kind_set)) THEN
         CPABORT("ECD| Implicit-Gamma SCF data are incomplete")
      END IF
      IF (SIZE(mos) < 1 .OR. SIZE(context%matrix_ks) < 1 .OR. SIZE(context%matrix_s) < 1) THEN
         CPABORT("ECD| Implicit-Gamma SCF data contain an empty MO or AO matrix set")
      END IF
      IF (eps_degeneracy <= 0.0_dp) THEN
         CPABORT("ECD| EPS_DEGENERACY must be positive for the implicit-Gamma data adapter")
      END IF
      IF (spatial_order < 0) THEN
         CPABORT("ECD| SPATIAL_ORDER must be non-negative for the implicit-Gamma data adapter")
      END IF

      context%nspin = SIZE(mos)
      IF (SIZE(context%matrix_ks) /= context%nspin) THEN
         CPABORT("ECD| Implicit-Gamma KS spin channels are inconsistent with the MO channels")
      END IF
      context%direction = direction
      context%ref_point = ref_point
      context%spatial_order = spatial_order
      ALLOCATE (context%spin(context%nspin))

      DO ispin = 1, context%nspin
         NULLIFY (eigenvalues, mo_coeff, occupations)
         CALL get_mo_set(mos(ispin), nao=nao, nmo=nmo, eigenvalues=eigenvalues, &
                         occupation_numbers=occupations, mo_coeff=mo_coeff)
         IF (nao < 1 .OR. nmo < 1 .OR. .NOT. ASSOCIATED(eigenvalues) .OR. &
             .NOT. ASSOCIATED(occupations) .OR. .NOT. ASSOCIATED(mo_coeff)) THEN
            CPABORT("ECD| Implicit-Gamma SCF did not provide complete molecular-orbital data")
         END IF
         IF (SIZE(eigenvalues) /= nmo .OR. SIZE(occupations) /= nmo) THEN
            CPABORT("ECD| Implicit-Gamma MO eigenvalue or occupation dimensions are inconsistent")
         END IF
         CALL cp_fm_get_info(mo_coeff, nrow_global=nrow, ncol_global=ncol)
         IF (nrow /= nao .OR. ncol /= nmo) THEN
            CPABORT("ECD| Implicit-Gamma MO coefficient dimensions are inconsistent")
         END IF

         context%spin(ispin)%nao = nao
         context%spin(ispin)%nmo = nmo
         context%spin(ispin)%eigenvalues => eigenvalues
         context%spin(ispin)%occupations => occupations
         context%spin(ispin)%mo_coeff => mo_coeff
         CALL ecd_band_manifolds_build(context%spin(ispin)%manifolds, eigenvalues, eps_degeneracy)
      END DO

   END SUBROUTINE ecd_gamma_context_create_implicit

! **************************************************************************************************
!> \brief Release ECD-owned metadata while preserving borrowed Quickstep objects.
! **************************************************************************************************
   SUBROUTINE ecd_gamma_context_release(context)
      TYPE(ecd_gamma_context_type), INTENT(INOUT)       :: context

      INTEGER                                            :: ispin

      IF (ALLOCATED(context%spin)) THEN
         DO ispin = 1, SIZE(context%spin)
            CALL ecd_band_manifolds_release(context%spin(ispin)%manifolds)
            NULLIFY (context%spin(ispin)%eigenvalues, context%spin(ispin)%occupations, &
                     context%spin(ispin)%mo_coeff)
         END DO
         DEALLOCATE (context%spin)
      END IF
      NULLIFY (context%cell, context%matrix_ks, context%matrix_s, context%dft_control, context%sab_orb, &
               context%sab_all, context%sap_ppnl, context%particle_set, context%qs_kind_set)
      context%nspin = 0
      context%spatial_order = 0
      context%direction = 0.0_dp
      context%ref_point = 0.0_dp

   END SUBROUTINE ecd_gamma_context_release

END MODULE qs_ecd_gamma_data
