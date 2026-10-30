/*
============================================================
 dnn.C: certified DNN lower bound for a QAP subproblem (see dnn.h).

 Port of experiments/dnn_bound/nug20/dnn.py. The math is the same; what
 differs is how it is computed:

 1. Products with V. dnn.py forms V explicitly and multiplies dense
    (m^2+1)-square matrices by it (O(m^6) per product). Here the basis
    Vb of e-perp (m x (m-1)) is the first m-1 columns of the Householder
    reflector H with H e_{m-1} = e/sqrt(m), so V'MV and VRV' are a
    reflection along each of the four tensor modes plus a 2x2 rotation:
    O(m^4) per product. Any orthonormal basis gives the same R-step and
    the same bound, because both are invariant under an orthogonal
    change of basis.

 2. Certificate. For a symmetric multiplier Z, every Y = V R V' feasible
    for the relaxation has trace(R) = trace(Y) = m+1 (the arrow
    constraint follows from the gangster zeros and range(V)), so

      DNN value >= min_{Y in P} <L + Z, Y> - (m+1) lambda_max(V'ZV)

    for any polyhedral P containing the feasible Y. dnn.py takes P = box
    {Y_00 = 1, gangster zeros, 0 <= Y <= 1} and Z projected so that
    V'ZV <= 0, at the cost of an extra eigendecomposition per check.
    Here (default) Z' = Z + beta (Y - V R V') from the R-step is used:
    V'Z'V = beta * (W - P_psd(W)) is negative semidefinite by
    construction, with no extra eigendecomposition. lambda_max(V'Z'V)
    is bounded rigorously by a Cholesky factorization of delta I - V'Z'V
    (Rump's criterion), and the bound pays (m+1) * delta.

    The polyhedral part (DNN_POLY_GL) uses a smaller P that still
    contains every feasible Y: each column of Y lies in range(V), so
    for column (ij) the m x m block Y_(kl),(ij) has all row and column
    sums equal to x_ij = Y_0,(ij), and x is doubly stochastic. The
    minimum over that set is a Gilmore-Lawler computation on M = L + Z:

      M_00 + LAP_m(c),  c_ij = 2 M_0,(ij) + M_(ij),(ij)
                               + LAP_{m-1}(M_(kl),(ij), k != i, l != j)

    It is never below the box value, has the same limit at the dual
    optimum, and in the nug20 tests reaches a given value 1.5-6x
    sooner. Every LAP value used is the objective of a dual-feasible
    point, so an inexact inner solve (DNN_LAP_REDUCE) is still valid.

 3. R-step eigensolver. dnn.py takes a full double eigendecomposition
    of W = V'(Y + Z/beta)V every iteration; it is most of the cost, and
    W has hundreds of eigenvalues of each sign, so a partial spectrum
    does not help. The default (DNN_EIG_SUBSPACE) tracks the negative
    eigenspace of W between iterations by a Chebyshev-filtered subspace
    iteration in single precision (BLAS-3), with a full ssyevd every 20
    iterations. R is formed as W minus the computed negative part, so
    W - R stays negative semidefinite exactly and the certificate in 2
    is unaffected by the eigensolver's accuracy; only the ADMM step is
    inexact, and the tests show the same iteration counts.

 4. Step size. dnn.py's penalty is beta = m/3 (for the objective
    scaled to max |L| = 1). The default here is 0.3 m/3, which
    certifies 2-3x sooner on the nug20/nug30 nodes tested. Prototype
    mode (dnntest -P) restores beta, the eigensolver and the box
    certificate of dnn.py and reproduces its certified values.

 5. Floating point. The certified value subtracts a safety margin
    built from standard error bounds with generous constants: Rump's
    Cholesky criterion, the error of the Householder transforms, and
    the summation error of the polyhedral part (accumulated in long
    double). This is not interval arithmetic, but the margins are
    orders of magnitude above the errors seen (DNNResult::margin
    reports the largest one used; ~1e-7 at m = 28 for nug30).

 Index conventions (as in dnn.py): the lifted index of x_ij (facility i
 at location j) is 1 + i*m + j, index 0 is the constant 1. Matrices of
 the lifted size n1 = m^2+1 are stored dense and row-major.
============================================================
*/

#include "dnn.h"
#include "assign.h"   // LowerSum
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <ctime>
#include <vector>
#include <algorithm>

/*
 LAPACK / BLAS. DNN_BLAS_PREFIX is prepended to the Fortran names
 (scipy-openblas exports scipy_dsyevd_ etc.); DNN_OPENBLAS enables the
 OpenBLAS thread control. Set both from latest/makefile. Character
 arguments are followed by their hidden Fortran lengths.
*/
#ifndef DNN_BLAS_PREFIX
#define DNN_BLAS_PREFIX
#endif
#define DNN_CAT2(a, b) a##b
#define DNN_CAT(a, b) DNN_CAT2(a, b)
#define F77(name) DNN_CAT(DNN_BLAS_PREFIX, name##_)

typedef int lapack_int;

extern "C" {
void F77(dsyevd)(const char *jobz, const char *uplo, const lapack_int *n, double *a,
                 const lapack_int *lda, double *w, double *work, const lapack_int *lwork,
                 lapack_int *iwork, const lapack_int *liwork, lapack_int *info,
                 size_t, size_t);
void F77(dsyevr)(const char *jobz, const char *range, const char *uplo, const lapack_int *n,
                 double *a, const lapack_int *lda, const double *vl, const double *vu,
                 const lapack_int *il, const lapack_int *iu, const double *abstol,
                 lapack_int *m, double *w, double *z, const lapack_int *ldz,
                 lapack_int *isuppz, double *work, const lapack_int *lwork,
                 lapack_int *iwork, const lapack_int *liwork, lapack_int *info,
                 size_t, size_t, size_t);
void F77(ssyevd)(const char *jobz, const char *uplo, const lapack_int *n, float *a,
                 const lapack_int *lda, float *w, float *work, const lapack_int *lwork,
                 lapack_int *iwork, const lapack_int *liwork, lapack_int *info,
                 size_t, size_t);
void F77(dsyrk)(const char *uplo, const char *trans, const lapack_int *n, const lapack_int *k,
                const double *alpha, const double *a, const lapack_int *lda,
                const double *beta, double *c, const lapack_int *ldc, size_t, size_t);
void F77(dpotrf)(const char *uplo, const lapack_int *n, double *a, const lapack_int *lda,
                 lapack_int *info, size_t);
void F77(spotrf)(const char *uplo, const lapack_int *n, float *a, const lapack_int *lda,
                 lapack_int *info, size_t);
void F77(sgemm)(const char *ta, const char *tb, const lapack_int *m, const lapack_int *n,
                const lapack_int *k, const float *alpha, const float *a, const lapack_int *lda,
                const float *b, const lapack_int *ldb, const float *beta, float *c,
                const lapack_int *ldc, size_t, size_t);
void F77(ssyrk)(const char *uplo, const char *trans, const lapack_int *n, const lapack_int *k,
                const float *alpha, const float *a, const lapack_int *lda,
                const float *beta, float *c, const lapack_int *ldc, size_t, size_t);
void F77(strsm)(const char *side, const char *uplo, const char *transa, const char *diag,
                const lapack_int *m, const lapack_int *n, const float *alpha, const float *a,
                const lapack_int *lda, float *b, const lapack_int *ldb,
                size_t, size_t, size_t, size_t);
#ifdef DNN_OPENBLAS
void DNN_CAT(DNN_BLAS_PREFIX, openblas_set_num_threads)(int);
#endif
}

#define DNN_STR2(x) #x
#define DNN_STR(x) DNN_STR2(x)

const char *dnn_lapack_name()
{
#ifdef DNN_OPENBLAS
  return "OpenBLAS (prefix '" DNN_STR(DNN_BLAS_PREFIX) "')";
#elif defined(DNN_ACCELERATE)
  return "Accelerate";
#else
  return "system LAPACK/BLAS";
#endif
}

bool dnn_set_blas_threads(int n)
{
#ifdef DNN_OPENBLAS
  DNN_CAT(DNN_BLAS_PREFIX, openblas_set_num_threads)(n);
  return true;
#else
  (void)n;
  return false;
#endif
}

const char *dnn_status_name(int s)
{
  switch (s) {
  case DNN_FATHOMED: return "fathomed";
  case DNN_CONVERGED: return "converged";
  case DNN_STALLED: return "stalled";
  case DNN_HOPELESS: return "hopeless";
  case DNN_ITER_LIMIT: return "iter_limit";
  case DNN_TIME_LIMIT: return "time_limit";
  default: return "error";
  }
}

DNNParams::DNNParams()
{
  max_iter = 4000;
  time_limit = 0;
  check_every = 10;
  tol = 1e-5;
  gamma = 1.618;
  beta = 0.3;
  poly = DNN_POLY_GL;
  inner_lap = DNN_LAP_AUTO;
  auto_gap = 0;
  prototype_cert = 0;
  eig_mode = DNN_EIG_SUBSPACE;
  refresh = 20;
  cheb_degree = 2;
  stall_window = 0;
  stall_lb = 0.05;
  stall_gap = 0.5;
  hopeless_after = 200;
  L_add = NULL;
  elim = NULL;
  L_add_abs = 0;
  child_bounds = 0;
  verbose = 0;
}

namespace {

const double U_ROUND = DBL_EPSILON / 2;   // unit roundoff 2^-53

double thread_cpu()
{
  timespec ts;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/*
 Dot products with four partial sums: without -ffast-math the compiler
 keeps a single sum in order, which serializes on the add latency. The
 order is fixed, so results stay reproducible.
*/
inline double dot4(const double *x, const double *y, size_t n)
{
  double a0 = 0, a1 = 0, a2 = 0, a3 = 0;
  size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    a0 += x[i] * y[i];
    a1 += x[i + 1] * y[i + 1];
    a2 += x[i + 2] * y[i + 2];
    a3 += x[i + 3] * y[i + 3];
  }
  for (; i < n; i++) a0 += x[i] * y[i];
  return (a0 + a1) + (a2 + a3);
}

/*
 Hungarian method (shortest augmenting paths, O(n^3)) on a dense n x n
 cost matrix. Returns a lower bound on the assignment value that is the
 objective of a dual-feasible point: the row duals of the method, with
 the column duals recomputed as v_j = min_i (c_ij - u_i), so the bound
 does not depend on the method having converged exactly.
*/
struct Hungarian {
  int nmax;
  std::vector<double> u, v, minv;
  std::vector<int> p, way;
  std::vector<char> used;

  void init(int n)
  {
    nmax = n;
    u.assign(n + 1, 0); v.assign(n + 1, 0); minv.assign(n + 1, 0);
    p.assign(n + 1, 0); way.assign(n + 1, 0); used.assign(n + 1, 0);
  }

  // c is row-major with leading dimension ldc
  long double solve(int n, const double *c, int ldc)
  {
    for (int j = 0; j <= n; j++) { u[j] = 0; v[j] = 0; p[j] = 0; way[j] = 0; }
    for (int i = 1; i <= n; i++) {
      p[0] = i;
      int j0 = 0;
      for (int j = 0; j <= n; j++) { minv[j] = HUGE_VAL; used[j] = 0; }
      do {
        used[j0] = 1;
        int i0 = p[j0], j1 = 0;
        double delta = HUGE_VAL;
        const double *row = c + (size_t)(i0 - 1) * ldc;
        for (int j = 1; j <= n; j++) {
          if (used[j]) continue;
          double cur = row[j - 1] - u[i0] - v[j];
          if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
          if (minv[j] < delta) { delta = minv[j]; j1 = j; }
        }
        if (j1 == 0) break;   // only with NaN costs; the bound below stays valid
        for (int j = 0; j <= n; j++) {
          if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
          else minv[j] -= delta;
        }
        j0 = j1;
      } while (p[j0] != 0);
      do { int j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
    }
    // dual-feasible repair and value
    long double lb = 0;
    for (int i = 1; i <= n; i++) lb += u[i];
    for (int j = 0; j < n; j++) {
      double mn = HUGE_VAL;
      for (int i = 0; i < n; i++) mn = std::min(mn, c[(size_t)i * ldc + j] - u[i + 1]);
      lb += mn;
    }
    return lb;
  }
};

// Row then column reduction: the cheapest dual-feasible LAP bound.
long double lap_reduce(int n, const double *c, int ldc, double *colmin)
{
  long double lb = 0;
  for (int j = 0; j < n; j++) colmin[j] = HUGE_VAL;
  for (int i = 0; i < n; i++) {
    const double *row = c + (size_t)i * ldc;
    double mn = HUGE_VAL;
    for (int j = 0; j < n; j++) mn = std::min(mn, row[j]);
    lb += mn;
    for (int j = 0; j < n; j++) colmin[j] = std::min(colmin[j], row[j] - mn);
  }
  for (int j = 0; j < n; j++) lb += colmin[j];
  return lb;
}

/*
 Work: everything one dnn_bound call needs. No state outlives the call.
*/
struct Work {
  const unsigned char *elim = NULL;   // eliminated assignments (DNNParams::elim)
  std::vector<int> elim_idx;          // their lifted indices 1 + i m + j
  std::vector<double> red;   // certify: reduced costs of the outer LAP (child bounds)
  bool want_red = false;
  int m, N, n1, K, ll;
  double tau;                 // H = I - tau u u'
  std::vector<double> hu;     // Householder vector u (length m)
  std::vector<int> kept;      // K kept transformed indices
  std::vector<int> kept_of;   // lifted index -> position in kept, or -1
  std::vector<double> L;      // unscaled symmetric lifted objective
  std::vector<double> Y, Z, VRV, Zc, T;
  std::vector<double> W, W2, F, R, tmp;
  std::vector<double> evals, ework;
  std::vector<float> Wf, evalsf, eworkf;   // single-precision eigensolver
  // subspace tracking of the negative eigenspace of W (DNN_EIG_SUBSPACE)
  std::vector<float> Xs, Ys, Ts, Hs;
  std::vector<double> Wprev;
  int kb, since_refresh;
  bool have_basis;
  double lam_hi, theta_top;
  std::vector<lapack_int> eiwork;
  std::vector<double> sub, cmat, colmin;
  Hungarian hung;
  double prof[8];             // CPU seconds per phase, printed with verbose >= 2
  double delta_rel;           // last relative Cholesky shift that worked (certify)

  void alloc(int m_)
  {
    m = m_; N = m * m; n1 = N + 1; K = (m - 1) * (m - 1) + 1;
    ll = 1 + (m - 1) * m + (m - 1);
    size_t nn = (size_t)n1 * n1, kk = (size_t)K * K;
    L.assign(nn, 0); Y.assign(nn, 0); Z.assign(nn, 0); VRV.assign(nn, 0);
    Zc.assign(nn, 0); T.assign(nn, 0);
    W.assign(kk, 0); W2.assign(kk, 0); F.assign(kk, 0); R.assign(kk, 0);
    tmp.assign(3 * (size_t)n1 + m + 8, 0);
    evals.assign(K, 0);
    // dsyevd workspace for JOBZ='V'
    ework.assign(1 + 6 * (size_t)K + 2 * kk, 0);
    eiwork.assign(3 + 5 * (size_t)K, 0);
    sub.assign((size_t)m * m, 0); cmat.assign((size_t)m * m, 0); colmin.assign(m, 0);
    hung.init(m);
    for (int t = 0; t < 8; t++) prof[t] = 0;
    delta_rel = 1e-13;

    // Householder vector mapping e/sqrt(m) to e_{m-1}: u = e/sqrt(m) - e_{m-1}
    hu.assign(m, 1.0 / sqrt((double)m));
    hu[m - 1] -= 1.0;
    double uu = 0;
    for (int i = 0; i < m; i++) uu += hu[i] * hu[i];
    tau = 2.0 / uu;

    kept.clear();
    kept.push_back(0);
    for (int a = 0; a < m - 1; a++)
      for (int b = 0; b < m - 1; b++) kept.push_back(1 + a * m + b);
    kept_of.assign(n1, -1);
    for (int t = 0; t < K; t++) kept_of[kept[t]] = t;
  }

  void alloc_single()
  {
    size_t kk = (size_t)K * K;
    Wf.assign(kk, 0); evalsf.assign(K, 0); eworkf.assign(1 + 6 * (size_t)K + 2 * kk, 0);
  }

  void alloc_subspace()
  {
    size_t kk = (size_t)K * K;
    alloc_single();
    Xs.assign(kk, 0); Ys.assign(kk, 0); Ts.assign(kk, 0); Hs.assign(kk, 0);
    Wprev.assign(kk, 0);
    kb = since_refresh = 0;
    have_basis = false;
    lam_hi = theta_top = 0;
  }
};

/*
 Transforms. G = Rot * diag(1, H kron H), where Rot rotates the pair of
 lifted coordinates (0, ll) by 45 degrees (ll = index of (m-1, m-1)).
 Then V = G' restricted to the kept columns, so V'MV = (G M G')[kept,
 kept] and V R V' = G' pad(R) G. G is orthogonal and both factors are
 symmetric involutions, so G^-1 = G' = diag(1, H kron H) * Rot.
*/

// X <- H X H for an m x m row-major block
inline void hh_block(Work &w, double *X)
{
  const int m = w.m;
  const double *u = w.hu.data(), tau = w.tau;
  double *s = w.tmp.data();
  for (int k = 0; k < m; k++) {
    double *row = X + k * m, d = tau * dot4(row, u, m);
    for (int l = 0; l < m; l++) row[l] -= d * u[l];
  }
  for (int l = 0; l < m; l++) s[l] = 0;
  for (int k = 0; k < m; k++) {
    const double *row = X + k * m;
    for (int l = 0; l < m; l++) s[l] += u[k] * row[l];
  }
  for (int k = 0; k < m; k++) {
    double *row = X + k * m, f = tau * u[k];
    for (int l = 0; l < m; l++) row[l] -= f * s[l];
  }
}

// rotation of the pair (0, ll) of a row: x <- x Rot
inline void col_rot(double *row, int ll)
{
  const double h = sqrt(0.5);
  double a = row[0], b = row[ll];
  row[0] = h * (a + b);
  row[ll] = h * (a - b);
}

// rows 0 and ll: M <- Rot M
void rows_rot(Work &w, double *M)
{
  const int n1 = w.n1;
  const double h = sqrt(0.5);
  double *r0 = M, *rl = M + (size_t)w.ll * n1;
  for (int c = 0; c < n1; c++) {
    double a = r0[c], b = rl[c];
    r0[c] = h * (a + b);
    rl[c] = h * (a - b);
  }
}

// M <- D M (D = diag(1, H kron H) on the row index). The rows 1 + i*m + j
// form an m x m grid of row vectors; H is applied along j and along i. Done
// in chunks of columns, so that a chunk of all N rows stays in cache and
// the matrix is read and written once.
void rows_D(Work &w, double *M)
{
  const int m = w.m, n1 = w.n1;
  const double *u = w.hu.data(), tau = w.tau;
  double *acc = w.tmp.data();
  const int CB = 48;
  for (int c0 = 0; c0 < n1; c0 += CB) {
    const int nc = std::min(CB, n1 - c0);
    for (int pass = 0; pass < 2; pass++)
      for (int g = 0; g < m; g++) {
        // pass 0: along j for fixed i = g; pass 1: along i for fixed j = g
        for (int c = 0; c < nc; c++) acc[c] = 0;
        for (int t = 0; t < m; t++) {
          int r = pass == 0 ? 1 + g * m + t : 1 + t * m + g;
          const double *row = M + (size_t)r * n1 + c0;
          double ut = u[t];
          for (int c = 0; c < nc; c++) acc[c] += ut * row[c];
        }
        for (int t = 0; t < m; t++) {
          int r = pass == 0 ? 1 + g * m + t : 1 + t * m + g;
          double *row = M + (size_t)r * n1 + c0, f = tau * u[t];
          for (int c = 0; c < nc; c++) row[c] -= f * acc[c];
        }
      }
  }
}

// symmetrize an n x n row-major matrix in place, (X + X')/2, in tiles
void symmetrize(double *X, int n)
{
  const int bs = 32;
  for (int ib = 0; ib < n; ib += bs)
    for (int jb = ib; jb < n; jb += bs)
      for (int i = ib; i < std::min(ib + bs, n); i++)
        for (int j = std::max(jb, i + 1); j < std::min(jb + bs, n); j++) {
          double a = 0.5 * (X[(size_t)i * n + j] + X[(size_t)j * n + i]);
          X[(size_t)i * n + j] = X[(size_t)j * n + i] = a;
        }
}

/*
 W = V' M V (K x K, symmetric) with M = A + coef * B (B may be NULL).
 Works in w.T (A may be w.T itself). Column side row by row (the row is in
 cache), then the row side in one chunked pass.
*/
void forward(Work &w, const double *A, double *Wout, const double *B = NULL,
             double coef = 0)
{
  const int n1 = w.n1, K = w.K, ll = w.ll;
  double *T = w.T.data();
  for (int r = 0; r < n1; r++) {
    double *row = T + (size_t)r * n1;
    const double *a = A + (size_t)r * n1;
    if (B) {
      const double *b = B + (size_t)r * n1;
      for (int c = 0; c < n1; c++) row[c] = a[c] + coef * b[c];
    } else if (A != T) {
      memcpy(row, a, n1 * sizeof(double));
    }
    hh_block(w, row + 1);
    col_rot(row, ll);
  }
  rows_D(w, T);
  rows_rot(w, T);
  for (int s = 0; s < K; s++) {
    const double *row = T + (size_t)w.kept[s] * n1;
    for (int t = 0; t < K; t++) Wout[(size_t)s * K + t] = row[w.kept[t]];
  }
  symmetrize(Wout, K);
}

// Mout = V R V' (n1 x n1, symmetric) for a full symmetric K x K R:
// pad R, then Rot and D on both sides (rows outside kept stay zero under
// the column-side operations)
void expand(Work &w, const double *R, double *Mout)
{
  const int K = w.K, n1 = w.n1, ll = w.ll;
  for (int r = 0; r < n1; r++) {
    double *row = Mout + (size_t)r * n1;
    memset(row, 0, n1 * sizeof(double));
    int s = w.kept_of[r];
    if (s < 0) continue;
    const double *rs = R + (size_t)s * K;
    for (int t = 0; t < K; t++) row[w.kept[t]] = rs[t];
    col_rot(row, ll);
    hh_block(w, row + 1);
  }
  rows_rot(w, Mout);
  rows_D(w, Mout);
  symmetrize(Mout, n1);
}

/*
 R = P_psd(W) for symmetric K x K W (row-major = column-major). Uses
 the eigendecomposition (dsyevd) and sums the smaller of the positive
 and negative parts with dsyrk. With neg_form, R = W - (negative part)
 always, so that W - R is negative semidefinite by construction (a sum
 of -|w| v v'), however inaccurate the eigenvectors: the certificate
 then passes its Cholesky check with a tiny shift. Writes R (full) into
 Rout != Win. Returns false on LAPACK failure. npos receives the
 positive count.
*/
bool proj_psd(Work &w, const double *Win, double *Rout, int *npos, bool single = false,
              bool neg_form = false)
{
  const lapack_int K = w.K;
  const size_t kk = (size_t)K * K;
  double *A = w.W2.data();
  lapack_int liwork = (lapack_int)w.eiwork.size(), info = 0;
  if (!single) {
    memcpy(A, Win, kk * sizeof(double));
    lapack_int lwork = (lapack_int)w.ework.size();
    F77(dsyevd)("V", "L", &K, A, &K, w.evals.data(), w.ework.data(), &lwork,
                w.eiwork.data(), &liwork, &info, 1, 1);
  } else {
    // about twice as fast; eigenvectors good to ~1e-7, enough for the ADMM
    // step (the certificate checks its own semidefiniteness in double)
    float *Af = w.Wf.data();
    for (size_t t = 0; t < kk; t++) Af[t] = (float)Win[t];
    lapack_int lwork = (lapack_int)w.eworkf.size();
    F77(ssyevd)("V", "L", &K, Af, &K, w.evalsf.data(), w.eworkf.data(), &lwork,
                w.eiwork.data(), &liwork, &info, 1, 1);
    for (int i = 0; i < K; i++) w.evals[i] = w.evalsf[i];
    for (size_t t = 0; t < kk; t++) A[t] = Af[t];
  }
  if (info != 0) return false;
  const double *ev = w.evals.data();
  int np = 0;
  for (int i = 0; i < K; i++) if (ev[i] > 0) np++;
  *npos = np;
  const int nneg = K - np;   // eigenvalues ascending: negatives first
  double *F = w.F.data();
  const double one = 1, zero = 0;
  if (np <= nneg && !neg_form) {
    // R = sum_{w>0} w v v'
    for (int c = 0; c < np; c++) {
      int i = nneg + c;
      double s = sqrt(ev[i]);
      const double *v = A + (size_t)i * K;
      double *f = F + (size_t)c * K;
      for (int r = 0; r < K; r++) f[r] = s * v[r];
    }
    lapack_int k = np;
    if (k > 0) F77(dsyrk)("L", "N", &K, &k, &one, F, &K, &zero, Rout, &K, 1, 1);
    else memset(Rout, 0, (size_t)K * K * sizeof(double));
  } else {
    // R = W - sum_{w<=0} w v v' = W + sum |w| v v'
    for (int c = 0; c < nneg; c++) {
      double s = sqrt(-ev[c]);
      const double *v = A + (size_t)c * K;
      double *f = F + (size_t)c * K;
      for (int r = 0; r < K; r++) f[r] = s * v[r];
    }
    memcpy(Rout, Win, (size_t)K * K * sizeof(double));
    lapack_int k = nneg;
    if (k > 0) F77(dsyrk)("L", "N", &K, &k, &one, F, &K, &one, Rout, &K, 1, 1);
  }
  // dsyrk filled the lower triangle (column-major) = upper triangle (row-major)
  // in our row-major view: copy it to the other half
  for (int c = 0; c < K; c++)
    for (int r = c + 1; r < K; r++) Rout[(size_t)r * K + c] = Rout[(size_t)c * K + r];
  return true;
}

/*
 R = P_psd(W) approximately, by tracking the negative eigenspace of W from
 one ADMM iteration to the next (DNN_EIG_SUBSPACE). W changes little per
 iteration, so a block of kb >= #negative eigenvalues + guard vectors from
 the previous iteration is refined by a Chebyshev filter of degree `degree`
 that damps [theta_top, lam_hi] (theta_top = the block's largest Ritz
 value, lam_hi an upper bound on the spectrum kept by Weyl's inequality)
 and a Rayleigh-Ritz step, all in single precision with BLAS-3. Every
 `refresh` iterations, or when the block looks too small or the
 orthonormalization fails, a full ssyevd restarts it.

 R = W - sum_{theta < 0} theta x x' (the negative form): W - R is negative
 semidefinite by construction, so certificates are unaffected by the
 approximation; only the ADMM step is inexact.
*/
bool proj_psd_subspace(Work &w, const double *Win, double *Rout, int *npos, int refresh,
                       int degree, int guard)
{
  const lapack_int K = w.K;
  const size_t kk = (size_t)K * K;
  float *Wf = w.Wf.data();
  const float one = 1, zero = 0;

  // Weyl: lambda_max(W) <= lambda_max(Wprev) + ||W - Wprev||_F
  if (w.have_basis) {
    double *d = w.F.data();   // scratch: F is rebuilt below
    for (size_t t = 0; t < kk; t++) d[t] = Win[t] - w.Wprev[t];
    w.lam_hi += sqrt(dot4(d, d, kk)) * (1 + 1e-6);
  }
  memcpy(w.Wprev.data(), Win, kk * sizeof(double));
  for (size_t t = 0; t < kk; t++) Wf[t] = (float)Win[t];

  bool full = !w.have_basis || w.since_refresh >= refresh || w.theta_top >= w.lam_hi;
  int nneg = 0;
  lapack_int kb = w.kb;
  float *X = w.Xs.data();
  std::vector<float> &theta = w.evalsf;
  while (true) {
    if (full) {
      lapack_int lwork = (lapack_int)w.eworkf.size(), liwork = (lapack_int)w.eiwork.size(),
                 info = 0;
      memcpy(X, Wf, kk * sizeof(float));
      F77(ssyevd)("V", "L", &K, X, &K, theta.data(), w.eworkf.data(), &lwork,
                  w.eiwork.data(), &liwork, &info, 1, 1);
      if (info != 0) return false;
      nneg = 0;
      while (nneg < K && theta[nneg] < 0) nneg++;
      kb = std::min<lapack_int>(K, nneg + guard);
      w.lam_hi = theta[K - 1] + 1e-6 * (fabs(theta[K - 1]) + fabs(theta[0]));
      w.since_refresh = 0;
      w.have_basis = true;
      break;
    }
    // Chebyshev filter: Y_{j+1} = 2 (W - c) Y_j / h - Y_{j-1}, Y_0 = X
    float *Y = w.Ys.data(), *T = w.Ts.data(), *Yp = X;
    const float c = (float)(0.5 * (w.lam_hi + w.theta_top)),
                h = (float)(0.5 * (w.lam_hi - w.theta_top));
    const size_t nb = (size_t)K * kb;
    F77(sgemm)("N", "N", &K, &kb, &K, &one, Wf, &K, Yp, &K, &zero, T, &K, 1, 1);
    for (size_t t = 0; t < nb; t++) Y[t] = (T[t] - c * Yp[t]) / h;
    for (int d = 2; d <= degree; d++) {
      F77(sgemm)("N", "N", &K, &kb, &K, &one, Wf, &K, Y, &K, &zero, T, &K, 1, 1);
      for (size_t t = 0; t < nb; t++) Yp[t] = 2 * (T[t] - c * Y[t]) / h - Yp[t];
      std::swap(Y, Yp);
    }
    // Y now holds the filtered block (in Ys or Xs); orthonormalize by CholQR2
    float *G = w.Hs.data();
    bool ok = true;
    for (int pass = 0; pass < 2 && ok; pass++) {
      F77(ssyrk)("U", "T", &kb, &K, &one, Y, &K, &zero, G, &kb, 1, 1);
      lapack_int info = 0;
      F77(spotrf)("U", &kb, G, &kb, &info, 1);
      if (info != 0) { ok = false; break; }
      F77(strsm)("R", "U", "N", "N", &K, &kb, &one, G, &kb, Y, &K, 1, 1, 1, 1);
    }
    if (!ok) { full = true; continue; }
    // Rayleigh-Ritz: H = Y' W Y, X = Y Q
    F77(sgemm)("N", "N", &K, &kb, &K, &one, Wf, &K, Y, &K, &zero, T, &K, 1, 1);
    F77(sgemm)("T", "N", &kb, &kb, &K, &one, Y, &K, T, &K, &zero, G, &kb, 1, 1);
    for (int i = 0; i < kb; i++)
      for (int j = i + 1; j < kb; j++)
        G[(size_t)i * kb + j] = G[(size_t)j * kb + i] =
            0.5f * (G[(size_t)i * kb + j] + G[(size_t)j * kb + i]);
    {
      lapack_int lwork = (lapack_int)w.eworkf.size(), liwork = (lapack_int)w.eiwork.size(),
                 info = 0;
      F77(ssyevd)("V", "L", &kb, G, &kb, theta.data(), w.eworkf.data(), &lwork,
                  w.eiwork.data(), &liwork, &info, 1, 1);
      if (info != 0) { full = true; continue; }
    }
    float *Xn = (Y == X) ? T : X;
    F77(sgemm)("N", "N", &K, &kb, &kb, &one, Y, &K, G, &kb, &zero, Xn, &K, 1, 1);
    if (Xn != X) memcpy(X, Xn, nb * sizeof(float));
    nneg = 0;
    while (nneg < kb && theta[nneg] < 0) nneg++;
    if (nneg > kb - (guard + 1) / 2) { full = true; continue; }   // block too small
    w.since_refresh++;
    kb = std::min<lapack_int>(kb, nneg + guard);   // keep the lowest Ritz vectors
    break;
  }
  w.kb = kb;
  w.theta_top = theta[kb - 1];
  *npos = K - nneg;

  // R = W - sum_{theta<0} theta x x'
  double *F = w.F.data();
  for (int i = 0; i < nneg; i++) {
    double sq = sqrt(-(double)theta[i]);
    const float *x = X + (size_t)i * K;
    double *f = F + (size_t)i * K;
    for (int r = 0; r < K; r++) f[r] = sq * x[r];
  }
  memcpy(Rout, Win, kk * sizeof(double));
  const double done = 1;
  lapack_int k = nneg;
  if (k > 0) F77(dsyrk)("L", "N", &K, &k, &done, F, &K, &done, Rout, &K, 1, 1);
  for (int c = 0; c < K; c++)
    for (int r = c + 1; r < K; r++) Rout[(size_t)r * K + c] = Rout[(size_t)c * K + r];
  return true;
}

/*
 Certificate. Given a symmetric multiplier Zm (scaled units) with
 V'Zm V (approximately) negative semidefinite, return a certified lower
 bound on the subproblem (unscaled, shift excluded), or -HUGE_VAL when
 the semidefiniteness check fails. Uses w.T, w.W, w.W2.
*/
struct CertInfo {
  double margin;             // total safety margin subtracted
  double sdp_penalty;        // its semidefinite part, (m+1) * scale * lambda bound
  bool chol_ok;
  double poly_reduce, poly;  // diagnostics (unscaled)
};

double certify(Work &w, const double *Zm, double scale, double sumabsL, int poly,
               int inner_lap, double exact_above, double exact_below, CertInfo &ci)
{
  const int m = w.m, n1 = w.n1;
  const lapack_int K = w.K;
  ci.chol_ok = false;
  ci.margin = ci.sdp_penalty = 0;
  ci.poly_reduce = ci.poly = -HUGE_VAL;

  // lambda_max(V'ZmV) <= lam, by Cholesky of delta I - S (Rump 2006: if
  // floating-point Cholesky of A - cI succeeds, A is positive definite)
  double *S = w.W.data();
  forward(w, Zm, S);
  // (1 + 1e-8): rounding of the sums of squares, with room for any m that fits in memory
  double znorm = sqrt(dot4(Zm, Zm, (size_t)n1 * n1)) * (1 + 1e-8);
  double snorm = sqrt(dot4(S, S, (size_t)K * K)) * (1 + 1e-8);
  // try delta = rel * ||S||_F, starting a decade below the last value that
  // worked and growing tenfold, up to 1e-3 (single-precision eigenvectors
  // typically need ~1e-7, double ~1e-13)
  double lam = HUGE_VAL;
  // (an integer exponent: repeated *= 10 can overshoot 1e-3 and skip the last try)
  int e = (int)floor(log10(std::max(1e-14, w.delta_rel / 10)) + 0.5);
  double *A = w.W2.data();
  for (e = std::max(e, -14); e <= -3 && lam == HUGE_VAL; e++) {
    const double rel = pow(10.0, e);
    double delta = rel * snorm + 1e-14 * znorm + 1e-300;
    long double tr = 0;
    double maxdiag = 0;
    for (size_t i = 0; i < (size_t)K * K; i++) A[i] = -S[i];
    for (int i = 0; i < K; i++) {
      A[(size_t)i * K + i] += delta;
      tr += A[(size_t)i * K + i];
      maxdiag = std::max(maxdiag, fabs(A[(size_t)i * K + i]));
    }
    lapack_int info = 0;
    F77(dpotrf)("L", &K, A, &K, &info, 1);
    if (info != 0) continue;
    double g = (K + 1) * U_ROUND / (1 - (K + 1) * U_ROUND);
    // c >= g/(1-g) * trace(A + cI): solve for c, generously rounded up
    double trA = (double)tr * (1 + 1e-12) + 1e-300;
    // 1.25: slack for LAPACK/OpenBLAS scaling by a reciprocal instead of dividing
    double c = (g / (1 - g)) * trA / (1 - K * g / (1 - g)) * 1.25
               + 4.0 * K * (2.0 * (K + 1) + maxdiag) * DBL_TRUE_MIN;
    lam = delta + c + U_ROUND * maxdiag;
    w.delta_rel = rel;
  }
  if (lam == HUGE_VAL) return -HUGE_VAL;
  ci.chol_ok = true;
  // error of the computed transform against an exact orthonormal V
  lam += 64.0 * m * U_ROUND * znorm;

  // polyhedral part on M = L + scale * Zm (unscaled)
  const double *L = w.L.data();
  long double absM = 0;
  long double m00 = L[0] + scale * Zm[0];
  long double val = 0;
  if (poly == DNN_POLY_BOX) {
    for (int r = 0; r < n1; r++) {
      const double *lr = L + (size_t)r * n1, *zr = Zm + (size_t)r * n1;
      int i = (r - 1) / m, j = (r - 1) % m;
      for (int c = 0; c < n1; c++) {
        if (r == 0 && c == 0) continue;
        if (r > 0 && c > 0) {
          int k = (c - 1) / m, l = (c - 1) % m;
          if ((k == i) != (l == j)) continue;   // gangster
        }
        // eliminated assignments: Y is zero there in every solution
        if ((r > 0 && w.elim && w.elim[r - 1]) || (c > 0 && w.elim && w.elim[c - 1])) continue;
        double x = lr[c] + scale * zr[c];
        absM += fabs(x);
        if (x < 0) val += x;
      }
    }
    val += m00;
    ci.poly_reduce = ci.poly = (double)val;
  } else {
    // Gilmore-Lawler on M: inner LAPs per row (ij), outer LAP on c
    double *sub = w.sub.data(), *cm = w.cmat.data();
    bool exact = inner_lap == DNN_LAP_EXACT;
    // Eliminated assignments get a penalty P in the inner and outer LAPs. Any
    // P is valid: permutations of the node never use them, so changing their
    // costs changes no feasible objective. P = 4 (m+1) max|M| keeps them out of
    // the LAP solutions whenever a feasible alternative exists, and is counted
    // in absM like any other entry, so the margin stays rigorous.
    double pen = 0;
    if (w.elim) {
      double mx = 0;
      for (size_t t = 0; t < (size_t)n1 * n1; t++) mx = std::max(mx, fabs(L[t] + scale * Zm[t]));
      pen = 4.0 * (m + 1) * mx;
    }
    for (int pass = 0; pass < 2; pass++) {
      absM = 0;
      for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) {
          int r = 1 + i * m + j;
          const double *lr = L + (size_t)r * n1, *zr = Zm + (size_t)r * n1;
          // sub = M_(ij),(kl) for k != i, l != j, as (m-1) x (m-1)
          int kk = 0;
          for (int k = 0; k < m; k++) {
            if (k == i) continue;
            int cb = 1 + k * m, ss = 0;
            for (int l = 0; l < m; l++) {
              if (l == j) continue;
              double x = lr[cb + l] + scale * zr[cb + l];
              if (w.elim && w.elim[k * m + l]) x += pen;
              sub[kk * (m - 1) + ss++] = x;
              absM += fabs(x);
            }
            kk++;
          }
          long double inner = 0;
          if (m > 1)
            inner = exact ? w.hung.solve(m - 1, sub, m - 1)
                          : lap_reduce(m - 1, sub, m - 1, w.colmin.data());
          double d0 = lr[0] + scale * zr[0], dd = lr[r] + scale * zr[r];
          absM += 2 * fabs(d0) + fabs(dd) + fabsl(inner);
          cm[i * m + j] = (double)(2.0L * d0 + dd + inner);
          if (w.elim && w.elim[i * m + j]) {
            cm[i * m + j] += pen;
            absM += pen;
          }
        }
      val = m00 + w.hung.solve(m, cm, m);
      if (pass == 0 && !exact) ci.poly_reduce = (double)val;
      if (exact || inner_lap != DNN_LAP_AUTO || (double)val <= exact_above ||
          (double)val > exact_below)
        break;
      exact = true;   // AUTO: close to the cutoff, redo with exact inner LAPs
    }
    ci.poly = (double)val;
    for (int t = 0; t < m * m; t++) absM += fabs(cm[t]);
    if (w.want_red) {
      // reduced costs of the outer LAP with the dual-repaired potentials that
      // gave val: v_j = min_i (cm_ij - u_i), r_ij = cm_ij - u_i - v_j >= 0, so
      // that every permutation with sigma(i) = j has polyhedral value at least
      // val + r_ij. Less a rounding allowance for the two subtractions.
      const double *u = w.hung.u.data() + 1;
      w.red.resize((size_t)m * m);
      for (int j = 0; j < m; j++) {
        double vj = HUGE_VAL;
        for (int i = 0; i < m; i++) vj = std::min(vj, cm[i * m + j] - u[i]);
        for (int i = 0; i < m; i++) {
          double c = cm[i * m + j];
          w.red[i * m + j] = ((c - u[i]) - vj)
                             - 4 * U_ROUND * (fabs(c) + fabs(u[i]) + fabs(vj));
        }
      }
    }
  }
  absM += fabsl(m00);
  // summation and rounding errors of the polyhedral part, objective rounding,
  // and the semidefinite term
  ci.sdp_penalty = (m + 1) * scale * lam;
  // (m + 2) u sum|L|: rounding of L itself and of cancellations in a row
  double margin = 8.0 * (m + 2) * U_ROUND * (double)absM + (m + 2) * U_ROUND * sumabsL
                  + ci.sdp_penalty;
  ci.margin = margin;
  // a NaN or infinity anywhere (OpenBLAS's Cholesky accepts NaN) certifies nothing
  if (!std::isfinite(lam) || !std::isfinite((double)val) || !std::isfinite(margin))
    return -HUGE_VAL;
  return (double)val - margin;
}

// whether some permutation avoids every eliminated assignment (elim is m x m,
// 1 = forbidden): augmenting paths (Kuhn), O(m^3)
static bool augment(const unsigned char *elim, int m, int i, std::vector<int> &match,
                    std::vector<char> &seen)
{
  for (int j = 0; j < m; j++) {
    if (elim[i * m + j] || seen[j]) continue;
    seen[j] = 1;
    if (match[j] < 0 || augment(elim, m, match[j], match, seen)) { match[j] = i; return true; }
  }
  return false;
}

bool has_perfect_matching(const unsigned char *elim, int m)
{
  std::vector<int> match(m, -1);
  for (int i = 0; i < m; i++) {
    std::vector<char> seen(m, 0);
    if (!augment(elim, m, i, match, seen)) return false;
  }
  return true;
}

// largest |x[t]|
double maxabs_of(const double *x, size_t n)
{
  double mx = 0;
  for (size_t t = 0; t < n; t++) mx = std::max(mx, fabs(x[t]));
  return mx;
}

// optimum for m <= 1 (the relaxation needs m >= 2), plus shift: the value
// rounded to nearest, and a lower bound on it (LowerSum: a product cancelling
// against c or the shift can lose more than an ulp of the result)
double trivial_value(const MAT *A, const MAT *B, const MAT *C, int m, double shift,
                     double &lower)
{
  LowerSum s;
  s.add(shift);
  if (m == 1) {
    s.add_product(A->me[0][0], B->me[0][0]);
    s.add(C->me[0][0]);
  }
  lower = s.lower();
  return (double)s.sum;
}

} // namespace

/*
 Warm start of a child from its parent's multiplier Z: the child fixes the
 parent's facility i to location j (parent-local indices). Lifted indices
 of the child map to the parent's with row i and column j removed; the
 child's index 0 (the constant 1) stands for the parent's 0 and (ij)
 together, since x_ij = 1 in the child, so Z is mapped like the objective
 (<Z, Y> is the same for every permutation of the child). Only Z is
 carried over: starting Y from the parent's (conditioned on x_ij = 1)
 made the child slower in tests, starting Z alone 2-3x faster.
*/
bool dnn_warm_child(const DNNWarm &par, int i, int j, DNNWarm &child)
{
  const int m = par.m, mc = m - 1;
  if (m < 3 || i < 0 || i >= m || j < 0 || j >= m ||
      par.Z.size() != (size_t)(m * m + 1) * (m * m + 1))
    return false;
  const int n1 = m * m + 1, nc = mc * mc + 1;
  std::vector<int> map(nc);   // child lifted index -> parent lifted index
  map[0] = 0;
  for (int k = 0, kc = 0; k < m; k++) {
    if (k == i) continue;
    for (int l = 0, lc = 0; l < m; l++) {
      if (l == j) continue;
      map[1 + kc * mc + lc] = 1 + k * m + l;
      lc++;
    }
    kc++;
  }
  const int r = 1 + i * m + j;
  const double *Z = par.Z.data();
  child.m = mc;
  child.scale = par.scale;
  child.Z.assign((size_t)nc * nc, 0);
  for (int a = 1; a < nc; a++) {
    int pa = map[a];
    for (int b = 1; b < nc; b++) child.Z[(size_t)a * nc + b] = Z[(size_t)pa * n1 + map[b]];
    child.Z[a] = child.Z[(size_t)a * nc] = Z[pa] + Z[(size_t)r * n1 + pa];
  }
  child.Z[0] = Z[0] + 2 * Z[r] + Z[(size_t)r * n1 + r];
  return true;
}

DNNResult dnn_bound(const MAT *A, const MAT *B, const MAT *C, double shift,
                    double cutoff, const DNNParams &par, const DNNWarm *start,
                    DNNWarm *save)
{
  DNNResult res;
  res.lb = -HUGE_VAL; res.primal = HUGE_VAL; res.status = DNN_ERROR;
  res.iterations = 0; res.cpu = 0; res.certificates = 0; res.cert_failures = 0;
  res.refinements = 0;
  res.margin = 0;
  double t0 = thread_cpu();

  if (!A || !B || !C || A->m != A->n || B->m != A->m || B->n != A->m ||
      C->m != A->m || C->n != A->m)
    return res;
  const int m = A->m;
  if (par.elim && !has_perfect_matching(par.elim, m)) {
    // no permutation avoids the eliminated assignments: no feasible completion
    res.lb = res.primal = HUGE_VAL;
    res.status = DNN_FATHOMED;
    res.cpu = thread_cpu() - t0;
    return res;
  }
  if (m <= 1) {
    res.primal = trivial_value(A, B, C, m, shift, res.lb);
    res.status = res.lb > cutoff ? DNN_FATHOMED : DNN_CONVERGED;
    res.cpu = thread_cpu() - t0;
    return res;
  }

  Work w;
  w.alloc(m);
  w.elim = par.elim;
  if (par.elim)
    for (int a = 0; a < m * m; a++) if (par.elim[a]) w.elim_idx.push_back(1 + a);
  if (par.eig_mode == DNN_EIG_SINGLE) w.alloc_single();
  if (par.eig_mode == DNN_EIG_SUBSPACE) w.alloc_subspace();
  const int n1 = w.n1;

  // lifted objective, symmetrized: L_(ij),(kl) = (a_ik b_jl + a_ki b_lj)/2,
  // L_0,(ij) = c_ij / 2. The margin's "rounding of L" term needs the sum of
  // the terms' magnitudes, not of |L|: when neither A nor B is symmetric the
  // two products can be large and cancel, and their rounding errors do not.
  double *L = w.L.data();
  double maxabs = 0;
  long double sumabsL = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) {
      double *row = L + (size_t)(1 + i * m + j) * n1;
      long double rowabs = 0;
      for (int k = 0; k < m; k++)
        for (int l = 0; l < m; l++) {
          double p1 = A->me[i][k] * B->me[j][l], p2 = A->me[k][i] * B->me[l][j];
          row[1 + k * m + l] = 0.5 * (p1 + p2);
          rowabs += 0.5 * (fabs(p1) + fabs(p2));
        }
      row[0] = L[1 + i * m + j] = 0.5 * C->me[i][j];
      sumabsL += rowabs + fabs(C->me[i][j]);
    }
  if (par.L_add) {
    for (size_t t = 0; t < (size_t)n1 * n1; t++) L[t] += par.L_add[t];
    sumabsL += par.L_add_abs + (long double)maxabs_of(par.L_add, (size_t)n1 * n1);
  }
  for (size_t t = 0; t < (size_t)n1 * n1; t++)
    maxabs = std::max(maxabs, fabs(L[t]));
  // as dnn.py: scale the objective to max |L| = 1, beta = m/3
  const double scale = std::max(1.0, maxabs), iscale = 1.0 / scale;
  const double beta = (par.beta > 0 ? par.beta : 1.0) * m / 3.0, ibeta = 1.0 / beta;
  const double gamma = par.gamma;
  const int check_every = std::max(1, par.check_every);

  double *Y = w.Y.data(), *Z = w.Z.data(), *VRV = w.VRV.data(), *Zc = w.Zc.data();
  Y[0] = 1.0;
  if (start && start->m == m && start->Z.size() == (size_t)n1 * n1 && start->scale > 0) {
    // multipliers are kept in unscaled units across nodes
    const double f = start->scale / scale;
    for (size_t t = 0; t < (size_t)n1 * n1; t++) Z[t] = start->Z[t] * f;
  }
  double best = -HUGE_VAL;
  std::vector<double> red_cur, red_best;   // child bounds (par.child_bounds)
  double lbs_best = -HUGE_VAL;
  std::vector<double> lbhist;   // certified lb at each check (for the stall test)
  std::vector<double> besthist; // running best at each check (for the hopeless test)
  int status = -1, npos = 0;
  double obj = HUGE_VAL;
  const double gap = par.auto_gap > 0 ? par.auto_gap : 1e-3 * std::max(1.0, fabs(cutoff));
  // no finite cutoff: never close to it, so never the exact pass (was NaN)
  const double exact_above = std::isfinite(cutoff) ? cutoff - gap - shift : HUGE_VAL;

  for (int it = 1; it <= par.max_iter; it++) {
    // R-step: W = V'(Y + Z/beta)V, R = P_psd(W), VRV = V R V'
    double tp = thread_cpu(), tq;
    forward(w, Y, w.W.data(), Z, ibeta);
    tq = thread_cpu(); w.prof[0] += tq - tp; tp = tq;
    bool check = it % check_every == 0 || it == par.max_iter;
    bool rstep_cert = check && !par.prototype_cert;
    bool okp = par.eig_mode == DNN_EIG_SUBSPACE
                   ? proj_psd_subspace(w, w.W.data(), w.R.data(), &npos, par.refresh,
                                       par.cheb_degree, std::max(8, w.K / 20))
                   : proj_psd(w, w.W.data(), w.R.data(), &npos,
                              par.eig_mode == DNN_EIG_SINGLE, rstep_cert);
    if (!okp) {
      status = DNN_ERROR;
      break;
    }
    tq = thread_cpu(); w.prof[1] += tq - tp; tp = tq;
    expand(w, w.R.data(), VRV);
    tq = thread_cpu(); w.prof[2] += tq - tp; tp = tq;

    // Y-step, Z-step (and the certificate multiplier Zc = Z + beta (Y - VRV)),
    // a row at a time: clip, then the fixed entries, then the updates
    long double pres2 = 0, dy2 = 0;
    double *yn = w.tmp.data();
    for (int r = 0; r < n1; r++) {
      const double *lr = L + (size_t)r * n1, *vr = VRV + (size_t)r * n1;
      double *yr = Y + (size_t)r * n1, *zr = Z + (size_t)r * n1, *cr = Zc + (size_t)r * n1;
      for (int c = 0; c < n1; c++) {
        double y = vr[c] - (lr[c] * iscale + zr[c]) * ibeta;
        yn[c] = std::min(1.0, std::max(0.0, y));
      }
      if (r == 0) {
        yn[0] = 1;
      } else {
        // gangster: same facility i or same location j, except (ij) itself
        int i = (r - 1) / m, j = (r - 1) % m;
        for (int l = 0; l < m; l++) if (l != j) yn[1 + i * m + l] = 0;
        for (int k = 0; k < m; k++) if (k != i) yn[1 + k * m + j] = 0;
      }
      if (w.elim) {
        // eliminated assignments: their lifted rows and columns are zero
        if (r > 0 && w.elim[r - 1]) {
          for (int c = 0; c < n1; c++) yn[c] = 0;
        } else {
          for (int a : w.elim_idx) yn[a] = 0;
        }
      }
      if (rstep_cert)
        for (int c = 0; c < n1; c++) cr[c] = zr[c] + beta * (yr[c] - vr[c]);
      double *rs = w.tmp.data() + n1, *ds = rs + n1;
      for (int c = 0; c < n1; c++) {
        double y = yn[c], res_ = y - vr[c];
        zr[c] += gamma * beta * res_;
        rs[c] = res_;
        ds[c] = y - yr[c];
        yr[c] = y;
      }
      pres2 += dot4(rs, rs, n1);
      dy2 += dot4(ds, ds, n1);
    }
    res.iterations = it;
    tq = thread_cpu(); w.prof[3] += tq - tp; tp = tq;
    if (!check) continue;

    // certified bound
    CertInfo ci;
    double lbs;
    if (par.prototype_cert) {
      // dnn.py: Zt = Z - V P_psd(V'ZV) V' (extra eigendecomposition)
      forward(w, Z, w.W.data());
      int np2;
      if (!proj_psd(w, w.W.data(), w.R.data(), &np2)) { status = DNN_ERROR; break; }
      expand(w, w.R.data(), Zc);
      for (size_t t = 0; t < (size_t)n1 * n1; t++) Zc[t] = Z[t] - Zc[t];
    }
    const double exact_below = cutoff - shift + 1e-6 * (1 + fabs(cutoff));
    w.want_red = par.child_bounds && par.poly != DNN_POLY_BOX;
    lbs = certify(w, Zc, scale, (double)sumabsL, par.poly, par.inner_lap, exact_above,
                  exact_below, ci);
    if (w.want_red) red_cur = w.red;
    res.certificates++;
    if (!ci.chol_ok) res.cert_failures++;
    if (!par.prototype_cert && par.eig_mode != DNN_EIG_DOUBLE &&
        (!ci.chol_ok || (lbs + shift <= cutoff && lbs + shift + ci.sdp_penalty > cutoff))) {
      // single-precision R-steps leave V'Z'V slightly indefinite; when that
      // costs the fathom, project it in double: Z' <- Z' - V P_psd(V'Z'V) V'
      forward(w, Zc, w.W.data());
      int np2;
      if (!proj_psd(w, w.W.data(), w.R.data(), &np2)) { status = DNN_ERROR; break; }
      expand(w, w.R.data(), VRV);   // VRV is free until the next R-step
      for (size_t t = 0; t < (size_t)n1 * n1; t++) Zc[t] -= VRV[t];
      CertInfo c2;
      double l2 = certify(w, Zc, scale, (double)sumabsL, par.poly, par.inner_lap,
                          exact_above, exact_below, c2);
      res.refinements++;
      if (l2 > lbs) { lbs = l2; ci = c2; if (w.want_red) red_cur = w.red; }
    }
    if (lbs > -HUGE_VAL) {
      double lb = nextafter(lbs + shift, -HUGE_VAL);   // covers the addition's rounding
      if (lb > best) {
        best = lb;
        if (w.want_red) { red_best = red_cur; lbs_best = lbs; }
      }
      res.margin = std::max(res.margin, ci.margin);
    }
    tq = thread_cpu(); w.prof[4] += tq - tp; tp = tq;
    obj = dot4(L, Y, (size_t)n1 * n1) + shift;
    double pres = sqrt((double)pres2), dy = sqrt((double)dy2);
    if (par.verbose)
      printf("it %5d lb %.4f (cert %.4f, margin %.2e%s) primal %.4f pres %.2e dY %.2e npos %d cpu %.2f\n",
             it, best, lbs > -HUGE_VAL ? lbs + shift : -HUGE_VAL, ci.margin,
             ci.chol_ok ? "" : ", chol FAILED", obj, pres, dy, npos, thread_cpu() - t0);
    lbhist.push_back(lbs > -HUGE_VAL ? lbs + shift : -HUGE_VAL);

    if (best > cutoff) { status = DNN_FATHOMED; break; }
    if (pres < par.tol && dy < par.tol) { status = DNN_CONVERGED; break; }
    if (par.time_limit > 0 && thread_cpu() - t0 > par.time_limit) { status = DNN_TIME_LIMIT; break; }
    besthist.push_back(best);
    if (par.hopeless_after > 0 && it >= par.hopeless_after && best > -HUGE_VAL) {
      // extrapolate the running best certified value: with d1, d2 its gains
      // over the last two windows of ~50 iterations and r = d1/d2, a
      // geometric tail would add d1 r/(1-r); allow 5 times that
      int w = std::max(1, (50 + check_every - 1) / check_every);
      int h = (int)besthist.size();
      if (h > 2 * w && besthist[h - 1 - 2 * w] > -HUGE_VAL) {
        double d1 = besthist[h - 1] - besthist[h - 1 - w];
        double d2 = besthist[h - 1 - w] - besthist[h - 1 - 2 * w];
        double r = d2 > 0 ? std::min(d1 / d2, 0.95) : (d1 > 0 ? 0.95 : 0.0);
        if (best + 5 * d1 * r / (1 - r) < cutoff) {
          status = DNN_HOPELESS;
          break;
        }
      }
    }
    if (par.stall_window > 0 && (int)lbhist.size() > par.stall_window) {
      // as dnn.py: the last stall_window certified values (not the running
      // best) span less than stall_lb, and the last is near the primal value
      double mx = -HUGE_VAL, mn = HUGE_VAL;
      for (size_t t = lbhist.size() - par.stall_window; t < lbhist.size(); t++) {
        mx = std::max(mx, lbhist[t]);
        mn = std::min(mn, lbhist[t]);
      }
      if (mx - mn < par.stall_lb && fabs(obj - lbhist.back()) < par.stall_gap) {
        status = DNN_STALLED;
        break;
      }
    }
  }
  if (status < 0) status = DNN_ITER_LIMIT;
  if (par.verbose >= 2)
    printf("profile (s): forward %.3f  eig %.3f  expand %.3f  Y/Z-step %.3f  certify %.3f\n",
           w.prof[0], w.prof[1], w.prof[2], w.prof[3], w.prof[4]);
  res.status = status;
  res.lb = best;
  res.primal = obj;
  if (par.child_bounds && status != DNN_FATHOMED && !red_best.empty()) {
    res.child_lb.resize((size_t)m * m);
    for (size_t t = 0; t < (size_t)m * m; t++) {
      LowerSum c;   // lbs + r_ij + shift, rounded down
      c.add(lbs_best);
      c.add(red_best[t]);
      c.add(shift);
      res.child_lb[t] = c.lower();
    }
  }
  if (save) {
    save->m = m;
    save->scale = scale;
    save->Z.assign(Z, Z + (size_t)n1 * n1);
    save->Y.assign(Y, Y + (size_t)n1 * n1);
  }
  res.cpu = thread_cpu() - t0;
  return res;
}
