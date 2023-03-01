/*
============================================================
  Linear-algebra-type stuff goes here.

  Basically, lower-level functions that operate on Meschach
  data structures like MAT and VEC go here.  Exception:
  the singular value decomposition routine and dependencies
  are in svd.C.

============================================================
*/

#ifndef LINALG_H
#define LINALG_H

#include "meschach.h"
#undef catch

// trace a*b'
inline double mm_trace(MAT *a,MAT *b)
{
  register unsigned int i,j;
  register double sum=0.0;
  for(i=0;i<a->m;i++) {
    for(j=0;j<a->n;j++)
      sum += a->me[i][j]*b->me[i][j];
  }
  return sum;
}

/* trace of X'*permmat(p) */
inline double mm_trace_alt(MAT *X,int *p)
{
  double sum=0.0;
  for(int i=0;i<X->n;i++)
    sum += X->me[i][*p++];
  //    sum += X->me[i][p[i]];

  return sum;
}

void nice_m_output(MAT *A);
void my_nice_m_output(MAT *X,char *title,int m,int n);

void my_svd(MAT *A, VEC *d);

// permute rows
void permrows(MAT *A,int *p,MAT *OUT);
void permrows0(MAT *A,int *p,MAT *OUT);

// permute cols
void permcols(MAT *A,int *p,MAT *OUT);

// OUT = A*B;
MAT *my_m_mlt(MAT *A,MAT *B,MAT *OUT);
// OUT = A'*B;
MAT *my_mtrm_mlt(MAT *A,MAT *B,MAT *OUT);
// OUT = A*diag(b)*A';
MAT *mvmt_mlt(MAT *A, VEC *b, MAT *R);
// OUT = A'*diag(b)*A;
MAT *mtvm_mlt(MAT *A, VEC *b, MAT *R);

// for a matrix V whose columns are the nullspace of e^T:
// compute V'*A*V, V*A*V'
void vavt(MAT *A, MAT *OUT);
void vtav(MAT *A, MAT *OUT);

void warmstart(MAT *X, MAT *X1, int di, int dj);
MAT *projectionMatrix(int n);


#endif
