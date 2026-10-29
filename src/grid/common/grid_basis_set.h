/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_BASIS_SET_H
#define GRID_BASIS_SET_H

/*******************************************************************************
 * \brief Internal representation of a basis set.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int nset;
  int nsgf;
  int maxco;
  int maxpgf;
  int *lmin;
  int *lmax;
  int *npgf;
  int *nsgf_set;
  int *first_sgf;
  double *sphi;
  double *zet;
  // Exact angular/radial factors; ncart == NULL selects the sphi fallback.
  int ncgf;
  int *ncart;      // [nset], contracted Cartesian counts per set
  int *first_cart; // [nset], zero-based offsets into cart
  int *cart;       // [ncgf], zero-based primitive Cartesian coset indices
  double *angular; // row-major [nsgf][ncgf]
  double *radial;  // row-major [ncgf][maxpgf], includes Cartesian normalization
} grid_basis_set;

#ifndef __cplusplus

/*******************************************************************************
 * \brief Allocates a basis set which can be passed to grid_create_task_list.
 *
 * \param nset            Number of sets this basis is composed of.
 * \param nsgf            Size of contracted spherical basis, ie. the block size
 * \param maxco           Maximum number of Cartesian functions across all sets.
 * \param maxpgf          Maximum number of primitive Gaussians across all sets.
 *
 *      The following params are given for each set:
 *
 * \param lmin            Lowest angular momentum.
 * \param lmax            Highest angular momentum.
 * \param npgf            Number of primitive Gaussians, ie. exponents.
 * \param nsgf_set        Number of spherical basis functions
 * \param first_sgf       Index of first spherical basis function (one based).
 * \param sphi            Transformation matrix for (de-)contracting the basis.
 * \param zet             Exponents of primitive Gaussians.
 *
 *      Optional exact angular/radial factors, copied only when they reconstruct
 *      every used sphi entry exactly. Pass zero and NULL pointers to use sphi.
 *
 * \param ncgf            Total number of contracted Cartesian functions.
 * \param ncart           Contracted Cartesian counts per set, or NULL.
 * \param cart            Primitive coset index per Cartesian function,
 * zero-based.
 * \param angular         Row-major angular transform [nsgf][ncgf].
 * \param radial          Row-major radial coefficients [ncgf][maxpgf].
 *
 * \param basis_set       Handle to the created basis set.
 *
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
                           grid_basis_set **basis_set);

/*******************************************************************************
 * \brief Deallocates given basis set.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_basis_set(grid_basis_set *basis_set);

#endif

#endif

// EOF
