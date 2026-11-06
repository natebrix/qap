/*
  MAT and VEC storage (see mat.h): Meschach 1.2b's allocation, resizing and
  copying routines, with the same results in every case, including the
  storage layout and which entries end up zero. Each function names the
  Meschach routine and file (src/MESCHACH/ before stage 3 of
  docs/meschach_removal.md) it reproduces. Meschach's memory accounting
  (mem_info, off unless a program turned it on; the solver never did) and its
  segmented-memory variant are left out.

  Derived from the Meschach Library, whose copyright notice follows. As it
  requires, this is a notice of the last modification: rewritten for
  Nathan Brixius's QAP solver (github.com/natebrix/qap), keeping only the
  routines in mat.h, on 6 November 2026.

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

#include "mat.h"

namespace {

// Meschach's error(): it printed the error and exited. These cases are bugs
// in the caller (or out of memory), so stop with the function's name.
[[noreturn]] void mat_fail(const char *fn, const char *what)
{
  fprintf(stderr, "\n%s: %s\n", fn, what);
  abort();
}

// calloc that fails loudly; n == 0 is allowed (Meschach got a non-NULL
// pointer from glibc's calloc(0) there)
template <class T> T *zalloc(size_t n, const char *fn)
{
  T *p = (T *)calloc(n, sizeof(T));
  if (p == NULL && n > 0) mat_fail(fn, "out of memory");
  return p;
}

// Meschach's RENEW: realloc, or calloc if p is NULL (new space is not zeroed
// by realloc; the callers zero what becomes visible). n > 0 at every call.
template <class T> T *renew(T *p, size_t n, const char *fn)
{
  T *q = p ? (T *)realloc(p, n * sizeof(T)) : (T *)calloc(n, sizeof(T));
  if (q == NULL) mat_fail(fn, "out of memory");
  return q;
}

// __zero__ (machine.c, with CHAR0ISDBL0): memset
inline void zero(double *dp, int len) { memset(dp, '\0', len * sizeof(double)); }

}  // namespace

/* m_get (memory.c) -- m x n matrix, zeroed, rows contiguous */
MAT *m_get(int m, int n)
{
  if (m < 0 || n < 0) mat_fail("m_get", "negative size");

  MAT *matrix = zalloc<MAT>(1, "m_get");
  matrix->m = m;      matrix->n = matrix->max_n = n;
  matrix->max_m = m;  matrix->max_size = m * n;
  matrix->base = zalloc<double>((size_t)m * n, "m_get");
  matrix->me = zalloc<double *>(m, "m_get");
  for (int i = 0; i < m; i++)
    matrix->me[i] = &(matrix->base[i * n]);
  return matrix;
}

/* m_free (memory.c) */
int m_free(MAT *mat)
{
  if (mat == MNULL || (int)(mat->m) < 0 || (int)(mat->n) < 0)
    return -1;  /* don't trust it */
  if (mat->base != NULL) free(mat->base);
  if (mat->me != NULL) free(mat->me);
  free(mat);
  return 0;
}

/* m_resize (memory.c) -- A becomes new_m x new_n. The top-left
   min(m,new_m) x min(n,new_n) block keeps its values, every other entry is
   zero, and the rows stay contiguous (me[i] = base + i*new_n): the data is
   moved within base when n changes. Storage only grows. Same size: A is
   returned untouched. */
MAT *m_resize(MAT *A, int new_m, int new_n)
{
  int i;
  int new_max_m, new_max_n, new_size, old_m, old_n;

  if (new_m < 0 || new_n < 0) mat_fail("m_resize", "negative size");

  if (!A)
    return m_get(new_m, new_n);

  /* nothing was changed */
  if (new_m == (int)A->m && new_n == (int)A->n)
    return A;

  old_m = A->m;  old_n = A->n;
  if (new_m > (int)A->max_m)
    A->me = renew(A->me, new_m, "m_resize");
  new_max_m = new_m >= (int)A->max_m ? new_m : (int)A->max_m;
  new_max_n = new_n >= (int)A->max_n ? new_n : (int)A->max_n;

  new_size = new_max_m * new_max_n;
  if (new_size > (int)A->max_size) {
    A->base = renew(A->base, new_size, "m_resize");
    A->max_size = new_size;
  }

  /* now set up A->me[i] */
  for (i = 0; i < new_m; i++)
    A->me[i] = &(A->base[i * new_n]);

  /* now shift data in matrix */
  const int keep_m = old_m <= new_m ? old_m : new_m;
  if (old_n > new_n) {
    for (i = 1; i < keep_m; i++)
      memmove(&(A->base[i * new_n]), &(A->base[i * old_n]), sizeof(double) * new_n);
  } else if (old_n < new_n) {
    for (i = keep_m - 1; i > 0; i--) {
      /* copy & then zero extra space */
      memmove(&(A->base[i * new_n]), &(A->base[i * old_n]), sizeof(double) * old_n);
      zero(&(A->base[i * new_n + old_n]), new_n - old_n);
    }
    /* the rest of row 0. Meschach also did this for new_m == 0, where row 0
       is not part of the matrix and may lie outside the storage (when A
       never had a row); skipping it then changes no entry of any matrix */
    if (new_m > 0)
      zero(&(A->base[old_n]), new_n - old_n);
    A->max_n = new_n;
  }
  /* zero out the new rows.. */
  for (i = old_m; i < new_m; i++)
    zero(&(A->base[i * new_n]), new_n);

  A->max_m = new_max_m;
  A->max_n = new_max_n;
  A->max_size = A->max_m * A->max_n;
  A->m = new_m;  A->n = new_n;

  return A;
}

/* m_zero (init.c) */
MAT *m_zero(MAT *A)
{
  if (A == MNULL) mat_fail("m_zero", "NULL matrix");
  int A_m = A->m, A_n = A->n;
  double **A_me = A->me;
  for (int i = 0; i < A_m; i++)
    zero(A_me[i], A_n);
  return A;
}

/* _m_copy (copy.c) -- out keeps its size if it is at least in's in both
   dimensions (only the top-left block is written); otherwise it is resized
   to in's size first. in == out: nothing to do. */
MAT *_m_copy(MAT *in, MAT *out, unsigned int i0, unsigned int j0)
{
  if (in == MNULL) mat_fail("_m_copy", "NULL matrix");
  if (in == out)
    return out;
  if (out == MNULL || out->m < in->m || out->n < in->n)
    out = m_resize(out, in->m, in->n);

  for (unsigned int i = i0; i < in->m; i++)
    memmove(&(out->me[i][j0]), &(in->me[i][j0]), (in->n - j0) * sizeof(double));

  return out;
}

/* v_get (memory.c) -- zeroed vector */
VEC *v_get(int size)
{
  if (size < 0) mat_fail("v_get", "negative size");

  VEC *vector = zalloc<VEC>(1, "v_get");
  vector->dim = vector->max_dim = size;
  vector->ve = zalloc<double>(size, "v_get");
  return vector;
}

/* v_free (memory.c) */
int v_free(VEC *vec)
{
  if (vec == VNULL || (int)(vec->dim) < 0)
    return -1;  /* don't trust it */
  if (vec->ve != NULL) free(vec->ve);
  free(vec);
  return 0;
}

/* v_resize (memory.c) -- x becomes new_dim long: entries below
   min(dim, new_dim) kept, new ones zero. Storage only grows. Same dim: x is
   returned untouched. A vector with max_dim == 0 (v_get(0)) is not resized:
   a new vector is returned and x is left as it was, as in Meschach (which
   took such a vector to be a view made by sub_vec). */
VEC *v_resize(VEC *x, int new_dim)
{
  if (new_dim < 0) mat_fail("v_resize", "negative size");

  if (!x)
    return v_get(new_dim);

  /* nothing is changed */
  if (new_dim == (int)x->dim)
    return x;

  if (x->max_dim == 0)
    return v_get(new_dim);

  if (new_dim > (int)x->max_dim) {
    x->ve = renew(x->ve, new_dim, "v_resize");
    x->max_dim = new_dim;
  }

  if (new_dim > (int)x->dim)
    zero(&(x->ve[x->dim]), new_dim - x->dim);
  x->dim = new_dim;

  return x;
}

/* _v_copy (copy.c) -- out keeps its size if it is at least as long as in */
VEC *_v_copy(VEC *in, VEC *out, unsigned int i0)
{
  if (in == VNULL) mat_fail("_v_copy", "NULL vector");
  if (in == out)
    return out;
  if (out == VNULL || out->dim < in->dim)
    out = v_resize(out, in->dim);

  memmove(&(out->ve[i0]), &(in->ve[i0]), (in->dim - i0) * sizeof(double));

  return out;
}

/* _in_prod (vecop.c) and __ip__ (machine.c, without VUNROLL): one running
   sum in index order */
double _in_prod(VEC *a, VEC *b, unsigned int i0)
{
  if (a == VNULL || b == VNULL) mat_fail("_in_prod", "NULL vector");
  unsigned int limit = a->dim <= b->dim ? a->dim : b->dim;
  if (i0 > limit) mat_fail("_in_prod", "index out of bounds");

  const double *dp1 = &(a->ve[i0]), *dp2 = &(b->ve[i0]);
  int len = (int)(limit - i0);
  double sum = 0.0;
  for (int i = 0; i < len; i++)
    sum += dp1[i] * dp2[i];
  return sum;
}
