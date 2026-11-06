/*
  MAT and VEC: the solver's dense matrices and vectors.

  These are the parts of the Meschach Library 1.2b (David E. Stewart and
  Zbigniew Leyk) that the solver used, kept with Meschach's names, field
  names, signatures and semantics so that the code written against Meschach
  is unchanged (stage 2 of docs/meschach_removal.md; the library itself was
  vendored in src/MESCHACH/ until stage 3). mat.C says, per function, which
  Meschach routine it reproduces and what that routine does in each case.

  - MAT: m x n, stored row by row in one block, base; me[i] = base + i*n.
    Code may rely on that layout (permcols in linalg.C walks me[0] across
    all rows). max_m, max_n and max_size describe the allocation, as in
    Meschach.
  - VEC: dim entries in ve; max_dim is the allocation. dim may be lowered by
    hand (qpbnd.C does gsort->dim--).
  - m_get/v_get return zeroed storage; m_resize keeps the top-left block and
    zeroes the rest; v_resize keeps the first entries and zeroes new ones;
    m_copy/v_copy resize the target only if it is NULL or too small.
  - Where Meschach called error() (negative size, NULL argument, out of
    memory), these print the function's name and abort.

  The row permutation Meschach's PERM and px_rows did is done with an int
  array and permrows (linalg.C).

  Derived from the Meschach Library, whose copyright notice follows. As it
  requires, this is a notice of the last modification: rewritten for
  Nathan Brixius's QAP solver (github.com/natebrix/qap), keeping only the
  routines above, on 6 November 2026.

**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
**
** This Meschach Library is provided "as is" without any express
** or implied warranty of any kind with respect to this software.
** In particular the authors shall not be liable for any direct,
** indirect, special, incidental or consequential damages arising
** in any way from use of the software.
**
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

#ifndef MAT_H
#define MAT_H

// The C headers that Meschach's matrix.h included and that the solver's files
// have relied on since (printf, sqrt, HUGE_VAL, memcpy, ...). Keeping them
// here keeps every call resolving to the same function as before.
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* vector: dim entries in ve[0..dim-1], storage for max_dim */
typedef struct {
  unsigned int dim, max_dim;
  double *ve;
} VEC;

/* matrix: m x n, me[i] = base + i*n; storage for max_size entries */
typedef struct {
  unsigned int m, n;
  unsigned int max_m, max_n, max_size;
  double **me, *base;
} MAT;

#define VNULL ((VEC *)NULL)
#define MNULL ((MAT *)NULL)

/* m x n matrix, all zeros */
MAT *m_get(int m, int n);
/* A resized to new_m x new_n: top-left block kept, new entries zero;
   A == NULL: m_get(new_m, new_n) */
MAT *m_resize(MAT *A, int new_m, int new_n);
/* frees A; returns -1 (and does nothing) if A is NULL */
int m_free(MAT *A);
/* A <- 0 */
MAT *m_zero(MAT *A);
/* copies rows i0.. and columns j0.. of in into out; out is resized to in's
   size first if it is NULL or smaller than in in either dimension */
MAT *_m_copy(MAT *in, MAT *out, unsigned int i0, unsigned int j0);

/* vector of dimension dim, all zeros */
VEC *v_get(int dim);
/* x resized to new_dim: first entries kept, new entries zero;
   x == NULL: v_get(new_dim) */
VEC *v_resize(VEC *x, int new_dim);
/* frees x; returns -1 (and does nothing) if x is NULL */
int v_free(VEC *x);
/* copies in[i0..] into out[i0..]; out is resized if NULL or shorter */
VEC *_v_copy(VEC *in, VEC *out, unsigned int i0);

/* sum over i0 <= i < min(a->dim, b->dim) of a[i]*b[i], summed in order */
double _in_prod(VEC *a, VEC *b, unsigned int i0);

#define m_copy(in, out) _m_copy(in, out, 0, 0)
#define v_copy(in, out) _v_copy(in, out, 0)
#define in_prod(a, b)   _in_prod(a, b, 0)

/* free and set the pointer to NULL */
#define M_FREE(mat) (m_free(mat), (mat) = (MAT *)NULL)
#define V_FREE(vec) (v_free(vec), (vec) = (VEC *)NULL)

#endif
