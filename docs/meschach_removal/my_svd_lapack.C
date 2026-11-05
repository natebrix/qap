/*
  Prototype: my_svd (latest/linalg.C) on LAPACK instead of Meschach's
  Hfactor/makeHQ/trieig. Appended to a scratch copy of linalg.C by
  build_lapack.sh, which renames the original my_svd to my_svd_meschach.

  MES_EIG selects the routine: meschach (default), dsyev (QR, the same
  algorithm family as Meschach's), dsyevd (divide and conquer) or dsyevr
  (MRRR). Same contract as my_svd: on return d holds the eigenvalues and
  the rows of A the eigenvectors. Eigenvalues come back ascending; every
  caller sorts them itself.
*/
#include <vector>
#include <cstring>
#include <cstdlib>

extern "C" {
void scipy_dsyev_(const char *, const char *, const int *, double *, const int *, double *,
                  double *, const int *, int *, size_t, size_t);
void scipy_dsyevd_(const char *, const char *, const int *, double *, const int *, double *,
                   double *, const int *, int *, const int *, int *, size_t, size_t);
void scipy_dsyevr_(const char *, const char *, const char *, const int *, double *, const int *,
                   const double *, const double *, const int *, const int *, const double *,
                   int *, double *, double *, const int *, int *, double *, const int *, int *,
                   const int *, int *, size_t, size_t, size_t);
void scipy_openblas_set_num_threads(int);
}

void my_svd_meschach(MAT *A, VEC *d);

namespace {
int eig_method()
{
  static int m = -1;
  if (m < 0) {
    const char *e = getenv("MES_EIG");
    m = (!e || !strcmp(e, "meschach")) ? 0 : !strcmp(e, "dsyev") ? 1 : !strcmp(e, "dsyevd") ? 2 : 3;
    if (m) scipy_openblas_set_num_threads(1);
  }
  return m;
}
std::vector<double> a_, z_, w_, work_;
std::vector<int> iwork_, isuppz_;
}

void my_svd(MAT *A, VEC *d)
{
  const int meth = eig_method();
  if (meth == 0) { my_svd_meschach(A, d); return; }
  const int n = A->m;
  d = v_resize(d, n);
  if (n == 1) { d->ve[0] = A->me[0][0]; A->me[0][0] = 1.0; return; }
  if ((int)a_.size() < n * n) { a_.resize(n * n); z_.resize(n * n); w_.resize(n); isuppz_.resize(2 * n); }
  // symmetric, so row major A is its own column-major transpose; the
  // eigenvectors come back in the columns of the column-major result,
  // i.e. the rows of the row-major copy
  for (int i = 0; i < n; i++) memcpy(&a_[i * n], A->me[i], n * sizeof(double));
  int N = n, info = 0, lw = -1, liw = -1, iq = 0;
  double q = 0;
  const double *vecs = a_.data();
  if (meth == 1) {
    scipy_dsyev_("V", "L", &N, a_.data(), &N, w_.data(), &q, &lw, &info, 1, 1);
    lw = (int)q; if ((int)work_.size() < lw) work_.resize(lw);
    scipy_dsyev_("V", "L", &N, a_.data(), &N, w_.data(), work_.data(), &lw, &info, 1, 1);
  } else if (meth == 2) {
    scipy_dsyevd_("V", "L", &N, a_.data(), &N, w_.data(), &q, &lw, &iq, &liw, &info, 1, 1);
    lw = (int)q; liw = iq;
    if ((int)work_.size() < lw) work_.resize(lw);
    if ((int)iwork_.size() < liw) iwork_.resize(liw);
    scipy_dsyevd_("V", "L", &N, a_.data(), &N, w_.data(), work_.data(), &lw, iwork_.data(), &liw, &info, 1, 1);
  } else {
    double vl = 0, vu = 0, abstol = 0; int il = 1, iu = n, m = 0;
    scipy_dsyevr_("V", "A", "L", &N, a_.data(), &N, &vl, &vu, &il, &iu, &abstol, &m, w_.data(), z_.data(),
                  &N, isuppz_.data(), &q, &lw, &iq, &liw, &info, 1, 1, 1);
    lw = (int)q; liw = iq;
    if ((int)work_.size() < lw) work_.resize(lw);
    if ((int)iwork_.size() < liw) iwork_.resize(liw);
    scipy_dsyevr_("V", "A", "L", &N, a_.data(), &N, &vl, &vu, &il, &iu, &abstol, &m, w_.data(), z_.data(),
                  &N, isuppz_.data(), work_.data(), &lw, iwork_.data(), &liw, &info, 1, 1, 1);
    vecs = z_.data();
  }
  if (info != 0) { fprintf(stderr, "my_svd: LAPACK info %d (n=%d)\n", info, n); exit(2); }
  for (int k = 0; k < n; k++) {
    d->ve[k] = w_[k];
    memcpy(A->me[k], vecs + k * n, n * sizeof(double));
  }
}
