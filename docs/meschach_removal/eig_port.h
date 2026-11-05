/*
  sym_eig_port: Meschach's symmetric eigendecomposition as the solver uses
  it (my_svd in latest/linalg.C: Hfactor, makeHQ, trieig, m_transp), ported
  to plain arrays with the same floating-point operations in the same order,
  so that the results are bit for bit Meschach's when compiled without FMA
  contraction (-ffp-contract=off, as latest/makefile does).

  Differences from Meschach that do not change any value:
  - the Householder vectors are read out of H once, not once per column of Q
    (makeHQ calls get_col n*(n-1) times);
  - Q is kept transposed, so trieig's column rotations are row rotations on
    contiguous memory and the final m_transp is gone;
  - no static work vectors (Meschach's Hfactor/makeHQ/hhtrcols keep some),
    no allocation per call after the first call at a given size.

  A: n row pointers to a symmetric matrix; on return row k holds the k-th
  eigenvector. d: the eigenvalues, in the (unsorted) order trieig leaves them.
*/
#ifndef EIG_PORT_H
#define EIG_PORT_H
void sym_eig_port(int n, double **A, double *d);
#endif
