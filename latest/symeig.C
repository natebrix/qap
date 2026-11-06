// See symeig.h. Each block names the Meschach routine it replaces
// (src/MESCHACH/hessen.c, hsehldr.c, givens.c, symmeig.c, machine.c).
#include "symeig.h"
#include <cmath>
#include <cfloat>
#include <vector>

namespace {
// work space, reused across calls (thread_local: no state shared between threads)
thread_local std::vector<double> H_, Qt_, col_, w_, hh_, beta_, diag_, b_;

// __ip__: sequential inner product
inline double ip(const double *x, const double *y, int len)
{
  double sum = 0.0;
  for (int i = 0; i < len; i++) sum += x[i] * y[i];
  return sum;
}
// __mltadd__: x += s*y
inline void mltadd(double *x, const double *y, double s, int len)
{
  for (int i = 0; i < len; i++) x[i] += s * y[i];
}
inline void givens_p(double x, double y, double *c, double *s)
{
  double norm = sqrt(x * x + y * y);
  if (norm == 0.0) { *c = 1.0; *s = 0.0; }
  else { *c = x / norm; *s = y / norm; }
}
// rot_cols(Q,i,k,c,s,Q) on Q stored transposed: rows i and k of Qt
inline void rot_rows_t(double *Qt, int n, int i, int k, double c, double s)
{
  double *qi = Qt + (size_t)i * n, *qk = Qt + (size_t)k * n;
  for (int j = 0; j < n; j++) {
    double temp = c * qi[j] + s * qk[j];
    qk[j] = -s * qi[j] + c * qk[j];
    qi[j] = temp;
  }
}
inline int sgn_p(double x) { return x >= 0 ? 1 : -1; }
}  // namespace

void sym_eig(int n, double **A, double *d)
{
  if (n == 1) { d[0] = A[0][0]; A[0][0] = 1.0; return; }
  if ((int)col_.size() < n) {
    H_.resize((size_t)n * n); Qt_.resize((size_t)n * n);
    col_.resize(n); w_.resize(n); beta_.resize(n); diag_.resize(n); b_.resize(n);
    hh_.resize((size_t)n * n);
  }
  double *H = H_.data(), *Qt = Qt_.data(), *t = col_.data(), *w = w_.data();
  double *beta = beta_.data(), *diag = diag_.data(), *b = b_.data();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) H[i * n + j] = A[i][j];

  // ---- Hfactor(H, diag, beta) ----
  const int limit = n - 1;
  for (int k = 0; k < limit; k++) {
    for (int i = 0; i < n; i++) t[i] = H[i * n + k];          // get_col
    // hhvec(t, k+1, &beta[k], t, &H[k+1][k])
    const int i0 = k + 1;
    double norm = sqrt(ip(t + i0, t + i0, n - i0));
    if (norm <= 0.0) {
      beta[k] = 0.0;
    } else {
      beta[k] = 1.0 / (norm * (norm + fabs(t[i0])));
      double newval = t[i0] > 0.0 ? -norm : norm;
      H[i0 * n + k] = newval;
      t[i0] -= newval;
    }
    diag[k] = t[k + 1];
    const double bk = beta[k];
    const int j0 = k + 1;
    if (bk != 0.0) {
      // hhtrcols(H, k+1, k+1, t, beta[k])
      for (int j = 0; j < n; j++) w[j] = 0.0;
      for (int i = i0; i < n; i++)
        if (t[i] != 0.0) mltadd(w + j0, H + i * n + j0, t[i], n - j0);
      for (int i = i0; i < n; i++)
        if (t[i] != 0.0) mltadd(H + i * n + j0, w + j0, -bk * t[i], n - j0);
      // hhtrrows(H, 0, k+1, t, beta[k])
      for (int i = 0; i < n; i++) {
        double scale = bk * ip(H + i * n + j0, t + j0, n - j0);
        if (scale == 0.0) continue;
        mltadd(H + i * n + j0, t + j0, -scale, n - j0);
      }
    }
  }

  // ---- makeHQ(H, diag, beta, Q), Q kept as Qt ----
  // Householder vector j (entries j+1..n-1): diag[j], then H[j+2..][j]
  double *hh = hh_.data();
  for (int j = 0; j < limit; j++) {
    double *v = hh + (size_t)j * n;
    for (int r = 0; r < n; r++) v[r] = H[r * n + j];
    v[j + 1] = diag[j];
  }
  for (int i = 0; i < n; i++) {
    double *q = Qt + (size_t)i * n;                            // tmp1 = e_i
    for (int j = 0; j < n; j++) q[j] = 0.0;
    q[i] = 1.0;
    for (int j = limit - 1; j >= 0; j--) {                     // hhtrvec(hh_j, beta[j], j+1, q, q)
      const double *v = hh + (size_t)j * n;
      const int r0 = j + 1;
      double scale = beta[j] * ip(v + r0, q + r0, n - r0);
      mltadd(q + r0, v + r0, -scale, n - r0);
    }
  }

  // ---- my_svd: tridiagonal from H ----
  for (int i = 0; i < n - 1; i++) { d[i] = H[i * n + i]; b[i] = H[i * n + i + 1]; }
  d[n - 1] = H[(n - 1) * n + (n - 1)];

  // ---- trieig(d, b, Q) ----
  const double SQRT2 = 1.4142135623730949;
  int i_min = 0;
  while (i_min < n) {
    int i_max = n - 1;
    for (int i = i_min; i < n - 1; i++)
      if (b[i] == 0.0) { i_max = i; break; }
    if (i_max <= i_min) { i_min = i_max + 1; continue; }
    bool split = false;
    while (!split) {
      double c, s, c2, s2, cs, z, bk, ak1, bk1, ak2, bk2;
      double dd = (d[i_max - 1] - d[i_max]) / 2;
      double b_sqr = b[i_max - 1] * b[i_max - 1];
      double mu = d[i_max] - b_sqr / (dd + sgn_p(dd) * sqrt(dd * dd + b_sqr));
      givens_p(d[i_min] - mu, b[i_min], &c, &s);
      s = -s;
      if (fabs(c) < SQRT2) { c2 = c * c; s2 = 1 - c2; }
      else { s2 = s * s; c2 = 1 - s2; }
      cs = c * s;
      ak1 = c2 * d[i_min] + s2 * d[i_min + 1] - 2 * cs * b[i_min];
      bk1 = cs * (d[i_min] - d[i_min + 1]) + (c2 - s2) * b[i_min];
      ak2 = s2 * d[i_min] + c2 * d[i_min + 1] + 2 * cs * b[i_min];
      bk2 = (i_min < i_max - 1) ? c * b[i_min + 1] : 0.0;
      z = (i_min < i_max - 1) ? -s * b[i_min + 1] : 0.0;
      d[i_min] = ak1;
      d[i_min + 1] = ak2;
      b[i_min] = bk1;
      if (i_min < i_max - 1) b[i_min + 1] = bk2;
      rot_rows_t(Qt, n, i_min, i_min + 1, c, -s);
      for (int i = i_min + 1; i < i_max; i++) {
        givens_p(b[i - 1], z, &c, &s);
        s = -s;
        if (fabs(c) < SQRT2) { c2 = c * c; s2 = 1 - c2; }
        else { s2 = s * s; c2 = 1 - s2; }
        cs = c * s;
        bk = c * b[i - 1] - s * z;
        ak1 = c2 * d[i] + s2 * d[i + 1] - 2 * cs * b[i];
        bk1 = cs * (d[i] - d[i + 1]) + (c2 - s2) * b[i];
        ak2 = s2 * d[i] + c2 * d[i + 1] + 2 * cs * b[i];
        bk2 = (i + 1 < i_max) ? c * b[i + 1] : 0.0;
        z = (i + 1 < i_max) ? -s * b[i + 1] : 0.0;
        d[i] = ak1; d[i + 1] = ak2;
        b[i] = bk1;
        if (i < i_max - 1) b[i + 1] = bk2;
        if (i > i_min) b[i - 1] = bk;
        rot_rows_t(Qt, n, i, i + 1, c, -s);
      }
      for (int i = i_min; i < i_max; i++)
        if (fabs(b[i]) < DBL_EPSILON * (fabs(d[i]) + fabs(d[i + 1]))) { b[i] = 0.0; split = true; }
    }
  }
  // m_transp(Q, A): A = Qt
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) A[i][j] = Qt[(size_t)i * n + j];
}
