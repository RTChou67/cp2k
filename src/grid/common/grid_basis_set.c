/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "grid_basis_set.h"
#include "grid_common.h"

/*******************************************************************************
 * \brief Allocates a basis set which can be passed to grid_create_task_list.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_create_basis_set(const int nset, const int nsgf, const int maxco,
                           const int maxpgf, const int lmin[nset],
                           const int lmax[nset], const int npgf[nset],
                           const int nsgf_set[nset], const int first_sgf[nset],
                           const double sphi[nsgf][maxco],
                           const double zet[nset][maxpgf], const int ncgf,
                           const int ncart[nset], const int *cart,
                           const double *angular, const double *radial,
                           grid_basis_set **basis_set_out) {

  grid_basis_set *basis_set = calloc(1, sizeof(grid_basis_set));
  assert(basis_set != NULL);

  basis_set->nset = nset;
  basis_set->nsgf = nsgf;
  basis_set->maxco = maxco;
  basis_set->maxpgf = maxpgf;

  size_t size = nset * sizeof(int);
  basis_set->lmin = malloc(size);
  assert(basis_set->lmin != NULL || size == 0);
  basis_set->lmax = malloc(size);
  assert(basis_set->lmax != NULL || size == 0);
  basis_set->npgf = malloc(size);
  assert(basis_set->npgf != NULL || size == 0);
  basis_set->nsgf_set = malloc(size);
  assert(basis_set->nsgf_set != NULL || size == 0);
  basis_set->first_sgf = malloc(size);
  assert(basis_set->first_sgf != NULL || size == 0);
  if (size != 0) {
    memcpy(basis_set->lmin, lmin, size);
    memcpy(basis_set->lmax, lmax, size);
    memcpy(basis_set->npgf, npgf, size);
    memcpy(basis_set->nsgf_set, nsgf_set, size);
    memcpy(basis_set->first_sgf, first_sgf, size);
  }

  size = nsgf * maxco * sizeof(double);
  basis_set->sphi = malloc(size);
  assert(basis_set->sphi != NULL || size == 0);
  if (size != 0) {
    memcpy(basis_set->sphi, sphi, size);
  }

  size = nset * maxpgf * sizeof(double);
  basis_set->zet = malloc(size);
  assert(basis_set->zet != NULL || size == 0);
  if (size != 0) {
    memcpy(basis_set->zet, zet, size);
  }

  // Validate the explicit factors before exposing them to any backend.
  // Generic and replay bases retain the original sphi representation.
  if (nset > 0 && nsgf > 0 && maxpgf > 0 && ncgf > 0 && ncart != NULL &&
      cart != NULL && angular != NULL && radial != NULL) {
    int *first_cart = malloc((size_t)nset * sizeof(int));
    assert(first_cart != NULL);
    bool exact = true;
    int offset = 0;
    for (int iset = 0; iset < nset && exact; iset++) {
      first_cart[iset] = offset;
      if (ncart[iset] < 0 || ncart[iset] > ncgf - offset) {
        exact = false;
        break;
      }
      offset += ncart[iset];
    }
    exact = exact && offset == ncgf;

    for (int iset = 0; iset < nset && exact; iset++) {
      const int nco = ncoset(lmax[iset]);
      const int sgf = first_sgf[iset] - 1;
      if (npgf[iset] < 0 || npgf[iset] > maxpgf ||
          (npgf[iset] > 0 && nco > maxco / npgf[iset]) || sgf < 0 ||
          sgf > nsgf || nsgf_set[iset] < 0 || nsgf_set[iset] > nsgf - sgf) {
        exact = false;
        break;
      }
      const int first = first_cart[iset];
      const int last = first + ncart[iset];
      for (int icart = first; icart < last; icart++) {
        if (cart[icart] < 0 || cart[icart] >= nco) {
          exact = false;
          break;
        }
      }
      for (int ipgf = 0; ipgf < npgf[iset] && exact; ipgf++) {
        for (int isgf = sgf; isgf < sgf + nsgf_set[iset] && exact; isgf++) {
          double reconstructed[nco];
          memset(reconstructed, 0, sizeof(reconstructed));
          for (int icart = first; icart < last; icart++) {
            reconstructed[cart[icart]] +=
                radial[(size_t)icart * maxpgf + ipgf] *
                angular[(size_t)isgf * ncgf + icart];
          }
          for (int ico = 0; ico < nco; ico++) {
            if (!isfinite(reconstructed[ico]) ||
                reconstructed[ico] != sphi[isgf][ipgf * nco + ico]) {
              exact = false;
              break;
            }
          }
        }
      }
    }

    if (exact) {
      basis_set->ncgf = ncgf;
      basis_set->first_cart = first_cart;
      size = (size_t)nset * sizeof(int);
      basis_set->ncart = malloc(size);
      assert(basis_set->ncart != NULL);
      memcpy(basis_set->ncart, ncart, size);
      size = (size_t)ncgf * sizeof(int);
      basis_set->cart = malloc(size);
      assert(basis_set->cart != NULL);
      memcpy(basis_set->cart, cart, size);
      size = (size_t)nsgf * ncgf * sizeof(double);
      basis_set->angular = malloc(size);
      assert(basis_set->angular != NULL);
      memcpy(basis_set->angular, angular, size);
      size = (size_t)ncgf * maxpgf * sizeof(double);
      basis_set->radial = malloc(size);
      assert(basis_set->radial != NULL);
      memcpy(basis_set->radial, radial, size);
    } else {
      free(first_cart);
    }
  }

  *basis_set_out = basis_set;
}

/*******************************************************************************
 * \brief Deallocates given basis set.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_basis_set(grid_basis_set *basis_set) {
  free(basis_set->lmin);
  free(basis_set->lmax);
  free(basis_set->npgf);
  free(basis_set->nsgf_set);
  free(basis_set->first_sgf);
  free(basis_set->sphi);
  free(basis_set->zet);
  free(basis_set->ncart);
  free(basis_set->first_cart);
  free(basis_set->cart);
  free(basis_set->angular);
  free(basis_set->radial);
  free(basis_set);
}

// EOF
