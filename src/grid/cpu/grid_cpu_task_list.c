/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/grid_common.h"
#include "grid_cpu_collocate.h"
#include "grid_cpu_integrate.h"
#include "grid_cpu_task_list.h"
#include "grid_cpu_task_list_internal.h"

/*******************************************************************************
 * \brief Comperator passed to qsort to compare two tasks.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_tasks(const void *a, const void *b) {
  const grid_cpu_task *task_a = a, *task_b = b;
  if (task_a->level != task_b->level) {
    return task_a->level - task_b->level;
  } else if (task_a->block_num != task_b->block_num) {
    return task_a->block_num - task_b->block_num;
  } else if (task_a->iset != task_b->iset) {
    return task_a->iset - task_b->iset;
  } else if (task_a->jset != task_b->jset) {
    return task_a->jset - task_b->jset;
  } else if (task_a->ipgf != task_b->ipgf) {
    return task_a->ipgf - task_b->ipgf;
  } else {
    return task_a->jpgf - task_b->jpgf;
  }
}
/*******************************************************************************
 * \brief Allocates a task list for the cpu backend.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_cpu_create_task_list(
    const bool orthorhombic, const int ntasks, const int nlevels,
    const int natoms, const int nkinds, const int nblocks,
    const int block_offsets[nblocks], const double atom_positions[natoms][3],
    const int atom_kinds[natoms], const grid_basis_set *basis_sets[nkinds],
    const int level_list[ntasks], const int iatom_list[ntasks],
    const int jatom_list[ntasks], const int iset_list[ntasks],
    const int jset_list[ntasks], const int ipgf_list[ntasks],
    const int jpgf_list[ntasks], const int border_mask_list[ntasks],
    const int block_num_list[ntasks], const double radius_list[ntasks],
    const double rab_list[ntasks][3], const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], grid_cpu_task_list **task_list_out) {

  if (*task_list_out != NULL) {
    // This is actually an opportunity to reuse some buffers.
    grid_cpu_free_task_list(*task_list_out);
  }

  grid_cpu_task_list_internal *task_list =
      malloc(sizeof(grid_cpu_task_list_internal));
  assert(task_list != NULL);

  task_list->orthorhombic = orthorhombic;
  task_list->ntasks = ntasks;
  task_list->nlevels = nlevels;
  task_list->natoms = natoms;
  task_list->nkinds = nkinds;
  task_list->nblocks = nblocks;

  size_t size = nblocks * sizeof(int);
  task_list->block_offsets = malloc(size);
  assert(task_list->block_offsets != NULL || size == 0);
  if (size != 0) {
    memcpy(task_list->block_offsets, block_offsets, size);
  }

  size = 3 * natoms * sizeof(double);
  task_list->atom_positions = malloc(size);
  assert(task_list->atom_positions != NULL || size == 0);
  if (size != 0) {
    memcpy(task_list->atom_positions, atom_positions, size);
  }

  size = natoms * sizeof(int);
  task_list->atom_kinds = malloc(size);
  assert(task_list->atom_kinds != NULL || size == 0);
  if (size != 0) {
    memcpy(task_list->atom_kinds, atom_kinds, size);
  }

  size = nkinds * sizeof(grid_basis_set *);
  task_list->basis_sets = malloc(size);
  assert(task_list->basis_sets != NULL || size == 0);
  if (size != 0) {
    memcpy(task_list->basis_sets, basis_sets, size);
  }

  size = ntasks * sizeof(grid_cpu_task);
  task_list->tasks = malloc(size);
  assert(task_list->tasks != NULL || size == 0);
#pragma omp parallel for schedule(static) if (ntasks > GRID_OMP_MIN_ITERATIONS)
  for (int i = 0; i < ntasks; i++) {
    task_list->tasks[i].level = level_list[i];
    task_list->tasks[i].iatom = iatom_list[i];
    task_list->tasks[i].jatom = jatom_list[i];
    task_list->tasks[i].iset = iset_list[i];
    task_list->tasks[i].jset = jset_list[i];
    task_list->tasks[i].ipgf = ipgf_list[i];
    task_list->tasks[i].jpgf = jpgf_list[i];
    task_list->tasks[i].border_mask = border_mask_list[i];
    task_list->tasks[i].block_num = block_num_list[i];
    task_list->tasks[i].radius = radius_list[i];
    task_list->tasks[i].rab[0] = rab_list[i][0];
    task_list->tasks[i].rab[1] = rab_list[i][1];
    task_list->tasks[i].rab[2] = rab_list[i][2];
  }

  // Store grid layouts.
  size = nlevels * sizeof(grid_cpu_layout);
  task_list->layouts = malloc(size);
  assert(task_list->layouts != NULL || size == 0);
  for (int level = 0; level < nlevels; level++) {
    for (int i = 0; i < 3; i++) {
      task_list->layouts[level].npts_global[i] = npts_global[level][i];
      task_list->layouts[level].npts_local[i] = npts_local[level][i];
      task_list->layouts[level].shift_local[i] = shift_local[level][i];
      task_list->layouts[level].border_width[i] = border_width[level][i];
      for (int j = 0; j < 3; j++) {
        task_list->layouts[level].dh[i][j] = dh[level][i][j];
        task_list->layouts[level].dh_inv[i][j] = dh_inv[level][i][j];
      }
    }
  }

  // Group equal primitive pairs as well as level, block_num, iset, and jset.
  qsort(task_list->tasks, ntasks, sizeof(grid_cpu_task), &compare_tasks);

  // Index set-pair groups by grid level and AO block.
  size = nlevels * nblocks * sizeof(int);
  task_list->first_level_block_group = malloc(size);
  assert(task_list->first_level_block_group != NULL || size == 0);
  task_list->last_level_block_group = malloc(size);
  assert(task_list->last_level_block_group != NULL || size == 0);
  const int nlevel_blocks = nlevels * nblocks;
#pragma omp parallel for schedule(static) if (nlevel_blocks >                  \
                                                  GRID_OMP_MIN_ITERATIONS)
  for (int i = 0; i < nlevel_blocks; i++) {
    task_list->first_level_block_group[i] = 0;
    task_list->last_level_block_group[i] = -1; // last < first means no tasks
  }
  int ngroups = 0;
  for (int i = 0; i < ntasks; i++) {
    const grid_cpu_task *t = &task_list->tasks[i];
    if (i == 0 || t->level != t[-1].level || t->block_num != t[-1].block_num ||
        t->iset != t[-1].iset || t->jset != t[-1].jset)
      ngroups++;
  }
  task_list->groups = malloc(ngroups * sizeof(grid_cpu_task_group));
  assert(task_list->groups != NULL || ngroups == 0);
  int first = 0;
  for (int igroup = 0; igroup < ngroups; igroup++) {
    grid_cpu_task *t = &task_list->tasks[first];
    const grid_basis_set *a = basis_sets[atom_kinds[t->iatom - 1] - 1];
    const grid_basis_set *b = basis_sets[atom_kinds[t->jatom - 1] - 1];
    const int iset = t->iset - 1, jset = t->jset - 1;
    const int ca = ncoset(a->lmax[iset]), cb = ncoset(b->lmax[jset]);
    const int na = a->nsgf_set[iset], nb = b->nsgf_set[jset];
    int last = first, npairs = 0, nleft = 0;
    while (last < ntasks) {
      grid_cpu_task *next = &task_list->tasks[last];
      if (next->level != t->level || next->block_num != t->block_num ||
          next->iset != t->iset || next->jset != t->jset)
        break;
      if (last == first || next->ipgf != next[-1].ipgf)
        nleft++;
      if (last == first || next->ipgf != next[-1].ipgf ||
          next->jpgf != next[-1].jpgf)
        npairs++;
      next->pair = npairs - 1;
      last++;
    }
    grid_cpu_task_group *group = &task_list->groups[igroup];
    *group = (grid_cpu_task_group){.first_task = first,
                                   .last_task = last - 1,
                                   .mode = GRID_CONTRACT_DENSE};
    // Compare the two dense GEMMs with shared half-transforms over active
    // primitive pairs. The factored path moves angular work to the smaller
    // contracted Cartesian space and applies sparse radial maps afterwards.
    const int ma = a->npgf[iset] * ca, mb = b->npgf[jset] * cb;
    const double dense_cost = 2.0 * nb * ma * (na + mb);
    double cost = dense_cost;
    const double selective = 2.0 * nb * ca * (na * nleft + cb * npairs);
    // Keep the two dense GEMMs for marginal savings: small GEMM dispatches and
    // sparse-loop overhead need a substantial arithmetic reduction to pay off.
    if (2.0 * selective < dense_cost) {
      group->mode = GRID_CONTRACT_SELECTIVE;
      cost = selective;
    }
    if (a->ncart != NULL && b->ncart != NULL && a->ncart[iset] > 0 &&
        b->ncart[jset] > 0) {
      const double xa = a->ncart[iset], xb = b->ncart[jset];
      const double factored =
          2.0 * nb * xa * (na + xb) + 2.0 * xb * (xa * nleft + ca * npairs);
      if (2.0 * factored < dense_cost && factored < cost)
        group->mode = GRID_CONTRACT_FACTORED;
    }
    group->ncoeff =
        group->mode == GRID_CONTRACT_DENSE ? ma * mb : npairs * ca * cb;
    const int slot = (t->level - 1) * nblocks + t->block_num - 1;
    if (task_list->last_level_block_group[slot] < 0)
      task_list->first_level_block_group[slot] = igroup;
    task_list->last_level_block_group[slot] = igroup;
    first = last;
  }

  // Find largest Cartesian subblock size.
  task_list->maxco = 0;
  for (int i = 0; i < nkinds; i++) {
    task_list->maxco = imax(task_list->maxco, task_list->basis_sets[i]->maxco);
  }

  // Initialize thread-local storage.
  size = omp_get_max_threads() * sizeof(double *);
  task_list->threadlocals = malloc(size);
  assert(task_list->threadlocals != NULL);
  memset(task_list->threadlocals, 0, size);
  size = omp_get_max_threads() * sizeof(size_t);
  task_list->threadlocal_sizes = malloc(size);
  assert(task_list->threadlocal_sizes != NULL);
  memset(task_list->threadlocal_sizes, 0, size);

  *task_list_out = task_list;
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 * \author Ole Schuett
 ******************************************************************************/
void grid_cpu_free_task_list(grid_cpu_task_list *ptr) {
  if (ptr == NULL)
    return;

  grid_cpu_task_list_internal *task_list = (grid_cpu_task_list_internal *)ptr;

  free(task_list->block_offsets);
  free(task_list->atom_positions);
  free(task_list->atom_kinds);
  free(task_list->basis_sets);
  free(task_list->tasks);
  free(task_list->groups);
  free(task_list->layouts);
  free(task_list->first_level_block_group);
  free(task_list->last_level_block_group);
  for (int i = 0; i < omp_get_max_threads(); i++) {
    if (task_list->threadlocals[i] != NULL) {
      free(task_list->threadlocals[i]);
    }
  }
  free(task_list->threadlocals);
  free(task_list->threadlocal_sizes);
  free(task_list);
}

/*******************************************************************************
 * \brief Prototype for BLAS dgemm.
 * \author Ole Schuett
 ******************************************************************************/
void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
            const int *k, const double *alpha, const double *a, const int *lda,
            const double *b, const int *ldb, const double *beta, double *c,
            const int *ldc);

/*******************************************************************************
 * \brief Convenient wrapper to hide Fortran nature of dgemm_, swapping a and b.
 * \author Ole Schuett
 ******************************************************************************/
static void dgemm(const char transa, const char transb, const int m,
                  const int n, const int k, const double alpha, const double *a,
                  const int lda, const double *b, const int ldb,
                  const double beta, double *c, const int ldc) {
  dgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c,
         &ldc);
}

/*******************************************************************************
 * \brief Apply a two-sided basis transform, or its adjoint, to one AO subblock.
 * \author Ole Schuett
 ******************************************************************************/
static void transform_matrix(const bool integrate, const bool transpose,
                             const int na, const int nb, const int ma,
                             const int mb, const double *sa, const int lda,
                             const double *sb, const int ldb, const int ldblock,
                             const double *source, double *target) {
  double work[nb * ma];
  if (!integrate) {
    dgemm(transpose ? 'N' : 'T', 'N', nb, ma, na, 1.0, source, ldblock, sa, lda,
          0.0, work, ma);
    dgemm('T', 'N', mb, ma, nb, 1.0, sb, ldb, work, ma, 0.0, target, ma);
  } else {
    dgemm('N', 'N', nb, ma, mb, 1.0, sb, ldb, source, ma, 0.0, work, ma);
    if (transpose)
      dgemm('N', 'T', nb, na, ma, 1.0, work, ma, sa, lda, 1.0, target, ldblock);
    else
      dgemm('N', 'T', na, nb, ma, 1.0, sa, lda, work, ma, 1.0, target, ldblock);
  }
}

/*******************************************************************************
 * \brief Accumulate a sparse radial transform along one Cartesian dimension.
 *        Each contracted Cartesian function maps to one primitive coset index.
 ******************************************************************************/
static inline void transform_radial(const bool integrate, const bool right,
                                    const int nvec, const int ncart,
                                    const int nco, const int *cart,
                                    const double *radial, const int stride,
                                    const double *source, double *target) {
  for (int c = 0; c < ncart; c++) {
    const double coefficient = radial[c * stride];
    for (int v = 0; v < nvec; v++) {
      const int contracted = right ? v * ncart + c : c * nvec + v;
      const int primitive = right ? v * nco + cart[c] : cart[c] * nvec + v;
      if (integrate)
        target[contracted] += coefficient * source[primitive];
      else
        target[primitive] += coefficient * source[contracted];
    }
  }
}

/*******************************************************************************
 * \brief Transform the unique primitive pairs of one set pair and grid level.
 *        Forward and adjoint share the same basis factors and pair selection.
 ******************************************************************************/
static void transform_group(const bool integrate,
                            const grid_cpu_task_list_internal *list,
                            const grid_cpu_task_group *group,
                            const double *source, double *target) {
  const grid_cpu_task *first = &list->tasks[group->first_task];
  const grid_basis_set *a =
      list->basis_sets[list->atom_kinds[first->iatom - 1] - 1];
  const grid_basis_set *b =
      list->basis_sets[list->atom_kinds[first->jatom - 1] - 1];
  const int iset = first->iset - 1, jset = first->jset - 1;
  const int ca = ncoset(a->lmax[iset]), cb = ncoset(b->lmax[jset]);
  const int na = a->nsgf_set[iset], nb = b->nsgf_set[jset];
  const int ia = a->first_sgf[iset] - 1, ib = b->first_sgf[jset] - 1;
  const bool transpose = first->iatom <= first->jatom;
  const int ld = transpose ? a->nsgf : b->nsgf;
  const int offset = transpose ? ib * ld + ia : ia * ld + ib;
  source += integrate ? 0 : offset;
  target += integrate ? offset : 0;

  if (group->mode == GRID_CONTRACT_DENSE) {
    transform_matrix(integrate, transpose, na, nb, a->npgf[iset] * ca,
                     b->npgf[jset] * cb, &a->sphi[ia * a->maxco], a->maxco,
                     &b->sphi[ib * b->maxco], b->maxco, ld, source, target);
    return;
  }

  const bool factored = group->mode == GRID_CONTRACT_FACTORED;
  const int xa = factored ? a->ncart[iset] : 1;
  const int xb = factored ? b->ncart[jset] : nb;
  const int fa = factored ? a->first_cart[iset] : 0;
  const int fb = factored ? b->first_cart[jset] : 0;
  double q[xa * xb], work[xb * ca];
  if (factored) {
    if (integrate)
      memset(q, 0, sizeof(q));
    else
      transform_matrix(false, transpose, na, nb, xa, xb,
                       &a->angular[ia * a->ncgf + fa], a->ncgf,
                       &b->angular[ib * b->ncgf + fb], b->ncgf, ld, source, q);
  }

  for (int i = group->first_task; i <= group->last_task;) {
    const int p = list->tasks[i].ipgf - 1;
    int end = i + 1;
    while (end <= group->last_task && list->tasks[end].ipgf - 1 == p)
      end++;
    if (integrate || factored)
      memset(work, 0, sizeof(work));
    const double *sa = &a->sphi[ia * a->maxco + p * ca];
    if (!integrate) {
      if (factored)
        transform_radial(false, true, xb, xa, ca, &a->cart[fa],
                         &a->radial[fa * a->maxpgf + p], a->maxpgf, q, work);
      else
        dgemm(transpose ? 'N' : 'T', 'N', nb, ca, na, 1.0, source, ld, sa,
              a->maxco, 0.0, work, ca);
    }
    for (int j = i; j < end;) {
      const grid_cpu_task *task = &list->tasks[j];
      const int pgb = task->jpgf - 1;
      const int pair_offset = task->pair * ca * cb;
      const double *input = integrate ? source + pair_offset : work;
      double *output = integrate ? work : target + pair_offset;
      if (factored) {
        if (!integrate)
          memset(output, 0, ca * cb * sizeof(double));
        transform_radial(integrate, false, ca, xb, cb, &b->cart[fb],
                         &b->radial[fb * b->maxpgf + pgb], b->maxpgf, input,
                         output);
      } else if (integrate) {
        dgemm('N', 'N', nb, ca, cb, 1.0, &b->sphi[ib * b->maxco + pgb * cb],
              b->maxco, input, ca, 1.0, output, ca);
      } else {
        dgemm('T', 'N', cb, ca, nb, 1.0, &b->sphi[ib * b->maxco + pgb * cb],
              b->maxco, input, ca, 0.0, output, ca);
      }
      do {
        j++;
      } while (j < end && list->tasks[j].pair == task->pair);
    }
    if (integrate) {
      if (factored)
        transform_radial(true, true, xb, xa, ca, &a->cart[fa],
                         &a->radial[fa * a->maxpgf + p], a->maxpgf, work, q);
      else if (transpose)
        dgemm('N', 'T', nb, na, ca, 1.0, work, ca, sa, a->maxco, 1.0, target,
              ld);
      else
        dgemm('N', 'T', na, nb, ca, 1.0, sa, a->maxco, work, ca, 1.0, target,
              ld);
    }
    i = end;
  }
  if (factored && integrate)
    transform_matrix(true, transpose, na, nb, xa, xb,
                     &a->angular[ia * a->ncgf + fa], a->ncgf,
                     &b->angular[ib * b->ncgf + fb], b->ncgf, ld, q, target);
}

/*******************************************************************************
 * \brief Collocate a range of tasks which are destined for the same grid level.
 * \author Ole Schuett
 ******************************************************************************/
static void collocate_one_grid_level(
    const grid_cpu_task_list *ptr, const int *first_block_group,
    const int *last_block_group, const enum grid_func func,
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const double *pab_blocks, offload_buffer *grid) {

  if (ptr == NULL)
    return;

  grid_cpu_task_list_internal *task_list = (grid_cpu_task_list_internal *)ptr;

// Using default(shared) because with GCC 9 the behavior around const changed:
// https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared)
  {
    const int ithread = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();

    // Matrix pab is re-used across tasks.
    double pab[imax(task_list->maxco * task_list->maxco, 1)];

    const int npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];
    const size_t grid_size = npts_local_total * sizeof(double);
    // Multiple threads need private grids; a single thread uses the output.
    if (nthreads > 1 && task_list->threadlocal_sizes[ithread] < grid_size) {
      if (task_list->threadlocals[ithread] != NULL) {
        free(task_list->threadlocals[ithread]);
      }
      task_list->threadlocals[ithread] = malloc(grid_size);
      assert(task_list->threadlocals[ithread] != NULL);
      task_list->threadlocal_sizes[ithread] = grid_size;
    }
    double *const my_grid =
        nthreads == 1 ? grid->host_buffer : task_list->threadlocals[ithread];
    memset(my_grid, 0, grid_size);

    // Each thread owns complete AO blocks and reuses coefficients within a
    // group.
    const int chunk_size = imax(1, task_list->nblocks / (nthreads * 50));
#pragma omp for schedule(dynamic, chunk_size)
    for (int block_num = 0; block_num < task_list->nblocks; block_num++) {
      const int first_group = first_block_group[block_num];
      const int last_group = last_block_group[block_num];

      for (int igroup = first_group; igroup <= last_group; igroup++) {
        const grid_cpu_task_group *group = &task_list->groups[igroup];
        const grid_cpu_task *first = &task_list->tasks[group->first_task];
        const int block_offset = task_list->block_offsets[first->block_num - 1];
        transform_group(false, task_list, group, &pab_blocks[block_offset],
                        pab);
        for (int itask = group->first_task; itask <= group->last_task;
             itask++) {
          // Define some convenient aliases.
          const grid_cpu_task *task = &task_list->tasks[itask];
          const int iatom = task->iatom - 1;
          const int jatom = task->jatom - 1;
          const int iset = task->iset - 1;
          const int jset = task->jset - 1;
          const int ipgf = task->ipgf - 1;
          const int jpgf = task->jpgf - 1;
          const int ikind = task_list->atom_kinds[iatom] - 1;
          const int jkind = task_list->atom_kinds[jatom] - 1;
          const grid_basis_set *ibasis = task_list->basis_sets[ikind];
          const grid_basis_set *jbasis = task_list->basis_sets[jkind];
          const double zeta = ibasis->zet[iset * ibasis->maxpgf + ipgf];
          const double zetb = jbasis->zet[jset * jbasis->maxpgf + jpgf];
          const int ncoseta = ncoset(ibasis->lmax[iset]);
          const int ncosetb = ncoset(jbasis->lmax[jset]);
          const bool dense = group->mode == GRID_CONTRACT_DENSE;
          const int ncoa = (dense ? ibasis->npgf[iset] : 1) * ncoseta;
          const int ncob = (dense ? jbasis->npgf[jset] : 1) * ncosetb;
          const double *pair = pab + (dense ? 0 : task->pair * ncoa * ncob);

          grid_cpu_collocate_pgf_product(
              /*orthorhombic=*/task_list->orthorhombic,
              /*border_mask=*/task->border_mask,
              /*func=*/func,
              /*la_max=*/ibasis->lmax[iset],
              /*la_min=*/ibasis->lmin[iset],
              /*lb_max=*/jbasis->lmax[jset],
              /*lb_min=*/jbasis->lmin[jset],
              /*zeta=*/zeta,
              /*zetb=*/zetb,
              /*rscale=*/(iatom == jatom) ? 1 : 2,
              /*dh=*/dh,
              /*dh_inv=*/dh_inv,
              /*ra=*/&task_list->atom_positions[3 * iatom],
              /*rab=*/task->rab,
              /*npts_global=*/npts_global,
              /*npts_local=*/npts_local,
              /*shift_local=*/shift_local,
              /*border_width=*/border_width,
              /*radius=*/task->radius,
              /*o1=*/dense ? ipgf * ncoseta : 0,
              /*o2=*/dense ? jpgf * ncosetb : 0,
              /*n1=*/ncoa,
              /*n2=*/ncob,
              /*pab=*/(const double (*)[ncoa])pair,
              /*grid=*/my_grid);

        } // end of task loop
      } // end of group loop
    } // end of block loop

// While there should be an implicit barrier at the end of the block loop, this
// explicit barrier eliminates occasional seg faults with icc compiled binaries.
#pragma omp barrier

    // Merge thread-local grids via an efficient tree reduction.
    const int nreduction_cycles = ceil(log(nthreads) / log(2)); // tree depth
    for (int icycle = 1; icycle <= nreduction_cycles; icycle++) {
      // Threads are divided into groups, whose size doubles with each cycle.
      // After a cycle the reduced data is stored at first thread of each group.
      const int group_size = 1 << icycle; // 2**icycle
      const int igroup = ithread / group_size;
      const int dest_thread = igroup * group_size;
      const int src_thread = dest_thread + group_size / 2;
      // The last group might actually be a bit smaller.
      const int actual_group_size = imin(group_size, nthreads - dest_thread);
      // Parallelize summation by dividing grid points across group members.
      const int rank = modulo(ithread, group_size); // position within the group
      const int64_t lb = ((int64_t)npts_local_total * rank) / actual_group_size;
      const int64_t ub =
          ((int64_t)npts_local_total * (rank + 1)) / actual_group_size;
      if (src_thread < nthreads) {
        GRID_PRAGMA_SIMD_LOOP
        for (int i = (int)lb; i < (int)ub; i++) {
          task_list->threadlocals[dest_thread][i] +=
              task_list->threadlocals[src_thread][i];
        }
      }
#pragma omp barrier
    }

    if (nthreads > 1) {
      // Copy final result from first thread into shared grid.
      const int64_t lb = ((int64_t)npts_local_total * ithread) / nthreads;
      const int64_t ub = ((int64_t)npts_local_total * (ithread + 1)) / nthreads;
      GRID_PRAGMA_SIMD_LOOP
      for (int i = (int)lb; i < (int)ub; i++) {
        grid->host_buffer[i] = task_list->threadlocals[0][i];
      }
    }

  } // end of omp parallel region
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_cpu_collocate_task_list(const grid_cpu_task_list *ptr,
                                  const enum grid_func func, const int nlevels,
                                  const offload_buffer *pab_blocks,
                                  offload_buffer *grids[nlevels]) {
  if (ptr == NULL)
    return;

  grid_cpu_task_list_internal *task_list = (grid_cpu_task_list_internal *)ptr;
  assert(task_list->nlevels == nlevels);

  for (int level = 0; level < task_list->nlevels; level++) {
    const int idx = level * task_list->nblocks;
    const int *first_block_group = &task_list->first_level_block_group[idx];
    const int *last_block_group = &task_list->last_level_block_group[idx];
    const grid_cpu_layout *layout = &task_list->layouts[level];
    collocate_one_grid_level(task_list, first_block_group, last_block_group,
                             func, layout->npts_global, layout->npts_local,
                             layout->shift_local, layout->border_width,
                             layout->dh, layout->dh_inv,
                             pab_blocks->host_buffer, grids[level]);
  }
}

/*******************************************************************************
 * \brief Integrate a range of tasks that belong to the same grid level.
 * \author Ole Schuett
 ******************************************************************************/
static void integrate_one_grid_level(
    const grid_cpu_task_list *ptr, const int *first_block_group,
    const int *last_block_group, const bool compute_tau, const int natoms,
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const offload_buffer *pab_blocks, const offload_buffer *grid,
    offload_buffer *hab_blocks, double forces[natoms][3], double virial[3][3]) {

  if (ptr == NULL)
    return;

  grid_cpu_task_list_internal *task_list = (grid_cpu_task_list_internal *)ptr;

// Using default(shared) because with GCC 9 the behavior around const changed:
// https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared)
  {
    // Matrix pab and hab are re-used across tasks.
    double pab[imax(task_list->maxco * task_list->maxco, 1)];
    double hab[imax(task_list->maxco * task_list->maxco, 1)];

    // Parallelize over blocks to avoid concurred access to hab_blocks.
    const int nthreads = omp_get_num_threads();
    const int chunk_size = imax(1, task_list->nblocks / (nthreads * 50));
#pragma omp for schedule(dynamic, chunk_size)
    for (int block_num = 0; block_num < task_list->nblocks; block_num++) {
      const int first_group = first_block_group[block_num];
      const int last_group = last_block_group[block_num];
      if (first_group > last_group)
        continue;

      // Accumulate forces per block as it corresponds to a pair of atoms.
      const grid_cpu_task *first =
          &task_list->tasks[task_list->groups[first_group].first_task];
      const int iatom = first->iatom - 1;
      const int jatom = first->jatom - 1;
      const int block_offset = task_list->block_offsets[block_num];
      const bool pab_required = (forces != NULL || virial != NULL);
      double my_forces[2][3] = {0};
      double my_virials[2][3][3] = {0};

      for (int igroup = first_group; igroup <= last_group; igroup++) {
        const grid_cpu_task_group *group = &task_list->groups[igroup];
        if (pab_required)
          transform_group(false, task_list, group,
                          &pab_blocks->host_buffer[block_offset], pab);
        memset(hab, 0, group->ncoeff * sizeof(double));
        for (int itask = group->first_task; itask <= group->last_task;
             itask++) {
          // Define some convenient aliases.
          const grid_cpu_task *task = &task_list->tasks[itask];
          assert(task->block_num - 1 == block_num);
          assert(task->iatom - 1 == iatom && task->jatom - 1 == jatom);
          const int ikind = task_list->atom_kinds[iatom] - 1;
          const int jkind = task_list->atom_kinds[jatom] - 1;
          grid_basis_set *ibasis = task_list->basis_sets[ikind];
          grid_basis_set *jbasis = task_list->basis_sets[jkind];
          const int iset = task->iset - 1;
          const int jset = task->jset - 1;
          const int ipgf = task->ipgf - 1;
          const int jpgf = task->jpgf - 1;
          const double zeta = ibasis->zet[iset * ibasis->maxpgf + ipgf];
          const double zetb = jbasis->zet[jset * jbasis->maxpgf + jpgf];
          const int ncoseta = ncoset(ibasis->lmax[iset]);
          const int ncosetb = ncoset(jbasis->lmax[jset]);
          const bool dense = group->mode == GRID_CONTRACT_DENSE;
          const int ncoa = (dense ? ibasis->npgf[iset] : 1) * ncoseta;
          const int ncob = (dense ? jbasis->npgf[jset] : 1) * ncosetb;
          const int pair_offset = dense ? 0 : task->pair * ncoa * ncob;

          grid_cpu_integrate_pgf_product(
              /*orthorhombic=*/task_list->orthorhombic,
              /*compute_tau=*/compute_tau,
              /*border_mask=*/task->border_mask,
              /*la_max=*/ibasis->lmax[iset],
              /*la_min=*/ibasis->lmin[iset],
              /*lb_max=*/jbasis->lmax[jset],
              /*lb_min=*/jbasis->lmin[jset],
              /*zeta=*/zeta,
              /*zetb=*/zetb,
              /*dh=*/dh,
              /*dh_inv=*/dh_inv,
              /*ra=*/&task_list->atom_positions[3 * iatom],
              /*rab=*/task->rab,
              /*npts_global=*/npts_global,
              /*npts_local=*/npts_local,
              /*shift_local=*/shift_local,
              /*border_width=*/border_width,
              /*radius=*/task->radius,
              /*o1=*/dense ? ipgf * ncoseta : 0,
              /*o2=*/dense ? jpgf * ncosetb : 0,
              /*n1=*/ncoa,
              /*n2=*/ncob,
              /*grid=*/grid->host_buffer,
              /*hab=*/(double (*)[ncoa])(hab + pair_offset),
              /*pab=*/
              (pab_required) ? (const double (*)[ncoa])(pab + pair_offset)
                             : NULL,
              /*forces=*/(forces != NULL) ? my_forces : NULL,
              /*virials=*/(virial != NULL) ? my_virials : NULL,
              /*hdab=*/NULL,
              /*hadb=*/NULL,
              /*a_hdab=*/NULL);

        } // end of task loop
        transform_group(true, task_list, group, hab,
                        &hab_blocks->host_buffer[block_offset]);
      } // end of group loop

      // Merge thread-local forces and virial into shared ones.
      // It does not seem worth the trouble to accumulate them thread-locally.
      const double scalef = (iatom == jatom) ? 1.0 : 2.0;
      if (forces != NULL) {
#pragma omp critical(forces)
        for (int i = 0; i < 3; i++) {
          forces[iatom][i] += scalef * my_forces[0][i];
          forces[jatom][i] += scalef * my_forces[1][i];
        }
      }
      if (virial != NULL) {
#pragma omp critical(virial)
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            virial[i][j] += scalef * my_virials[0][i][j];
            virial[i][j] += scalef * my_virials[1][i][j];
          }
        }
      }

    } // end of block loop

  } // end of omp parallel region
}

/*******************************************************************************
 * \brief Integrate all tasks of in given list from given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_cpu_integrate_task_list(
    const grid_cpu_task_list *ptr, const bool compute_tau, const int natoms,
    const int nlevels, const offload_buffer *pab_blocks,
    const offload_buffer *grids[nlevels], offload_buffer *hab_blocks,
    double forces[natoms][3], double virial[3][3]) {

  if (ptr == NULL)
    return;

  grid_cpu_task_list_internal *task_list = (grid_cpu_task_list_internal *)ptr;

  assert(task_list->nlevels == nlevels);
  assert(task_list->natoms == natoms);
  assert(hab_blocks != NULL);

  // Zero result arrays.
  if (hab_blocks->size != 0) {
    memset(hab_blocks->host_buffer, 0, hab_blocks->size);
  }
  if (forces != NULL) {
    memset(forces, 0, natoms * 3 * sizeof(double));
  }
  if (virial != NULL) {
    memset(virial, 0, 9 * sizeof(double));
  }

  for (int level = 0; level < task_list->nlevels; level++) {
    const int idx = level * task_list->nblocks;
    const int *first_block_group = &task_list->first_level_block_group[idx];
    const int *last_block_group = &task_list->last_level_block_group[idx];
    const grid_cpu_layout *layout = &task_list->layouts[level];
    integrate_one_grid_level(
        task_list, first_block_group, last_block_group, compute_tau, natoms,
        layout->npts_global, layout->npts_local, layout->shift_local,
        layout->border_width, layout->dh, layout->dh_inv, pab_blocks,
        grids[level], hab_blocks, forces, virial);
  }
}

// EOF
