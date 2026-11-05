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

#include "amx_gemm.h"   // fused split-bf16 W*Y on Intel AMX (DNN_AMX=1)
#include "dnn.h"
#include "assign.h"   // LowerSum
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <ctime>
#include <vector>
#include <algorithm>
#include <string>
#include <atomic>

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
void F77(sgemv)(const char *trans, const lapack_int *m, const lapack_int *n, const float *alpha,
                const float *a, const lapack_int *lda, const float *x, const lapack_int *incx,
                const float *beta, float *y, const lapack_int *incy, size_t);
void F77(ssyrk)(const char *uplo, const char *trans, const lapack_int *n, const lapack_int *k,
                const float *alpha, const float *a, const lapack_int *lda,
                const float *beta, float *c, const lapack_int *ldc, size_t, size_t);
void F77(strtri)(const char *uplo, const char *diag, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *info, int, int);
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
  gamma1 = 0;
  beta = 0.3;
  poly = DNN_POLY_GL;
  inner_lap = DNN_LAP_AUTO;
  auto_gap = 0;
  prototype_cert = 0;
  eig_mode = DNN_EIG_SUBSPACE;
  refresh = 20;
  cheb_degree = 2;
  tight_lam = 0;
  halpern = 0;
  anderson = 0;
  momentum = 0;
  momentum_cap = 1.0;
  stall_window = 0;
  stall_lb = 0.05;
  stall_gap = 0.5;
  hopeless_after = 200;
  hopeless_gap_ratio = 0;
  hopeless_zone = 0;
  hopeless_delay = 100;
  L_add = NULL;
  elim = NULL;
  L_add_abs = 0;
  warm_primal = 0;
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
 Type-II Anderson acceleration of the fixed-point iteration u -> T(u) on the
 ADMM state u = (Y, Z/beta), upper triangles packed (DNNParams::anderson).
 With g = T(u), f = g - u and the last `mem` differences dF, dG of f and g,
   gamma = argmin |f - dF gamma|  (ridge-regularized normal equations),
   u_next = g - dG gamma.
 Safeguard: when an extrapolated point's residual |f| grew past 2x the last
 one, the memory is cleared and the plain step g is taken.
*/
struct Anderson {
  int mem = 0, n1 = 0, cnt = 0, head = 0;
  size_t P = 0, N = 0;
  std::vector<double> u, fprev, gprev, dF, dG, gram, cur;
  double rprev = -1;
  bool have_prev = false, last_extrap = false;
  int steps = 0, resets = 0;
  void init(int mem_, int n1_) {
    mem = mem_; n1 = n1_;
    P = (size_t)n1 * (n1 + 1) / 2; N = 2 * P;
    u.assign(N, 0); fprev.assign(N, 0); gprev.assign(N, 0); cur.assign(N, 0);
    dF.assign((size_t)mem * N, 0); dG.assign((size_t)mem * N, 0);
    gram.assign((size_t)mem * mem, 0);
  }
  void pack(const double *Y, const double *Z, double ibeta, double *out) const {
    size_t t = 0;
    for (int r = 0; r < n1; r++) {
      const double *yr = Y + (size_t)r * n1;
      for (int c = r; c < n1; c++) out[t++] = yr[c];
    }
    for (int r = 0; r < n1; r++) {
      const double *zr = Z + (size_t)r * n1;
      for (int c = r; c < n1; c++) out[t++] = zr[c] * ibeta;
    }
  }
  void unpack(const double *in, double beta, double *Y, double *Z) const {
    size_t t = 0;
    for (int r = 0; r < n1; r++)
      for (int c = r; c < n1; c++, t++)
        Y[(size_t)r * n1 + c] = Y[(size_t)c * n1 + r] = in[t];
    for (int r = 0; r < n1; r++)
      for (int c = r; c < n1; c++, t++)
        Z[(size_t)r * n1 + c] = Z[(size_t)c * n1 + r] = in[t] * beta;
  }
  // before the iteration: remember its input
  void begin(const double *Y, const double *Z, double ibeta) { pack(Y, Z, ibeta, u.data()); }
  // after the iteration (Y, Z) = T(u): replace them by the accelerated point
  void step(double *Y, double *Z, double beta) {
    double *g = cur.data();
    pack(Y, Z, 1 / beta, g);
    // f = g - u, kept in u
    double *f = u.data();
    for (size_t t = 0; t < N; t++) f[t] = g[t] - f[t];
    const double r = sqrt(dot4(f, f, N));
    if (last_extrap && rprev > 0 && r > 2 * rprev) {
      cnt = 0; head = 0; resets++;
      have_prev = false;
    }
    if (have_prev) {
      double *df = dF.data() + (size_t)head * N, *dg = dG.data() + (size_t)head * N;
      for (size_t t = 0; t < N; t++) { df[t] = f[t] - fprev[t]; dg[t] = g[t] - gprev[t]; }
      cnt = std::min(cnt + 1, mem);
      for (int k = 0; k < cnt; k++) {
        double v = dot4(df, dF.data() + (size_t)k * N, N);
        gram[(size_t)head * mem + k] = gram[(size_t)k * mem + head] = v;
      }
      head = (head + 1) % mem;
    }
    memcpy(fprev.data(), f, N * sizeof(double));
    memcpy(gprev.data(), g, N * sizeof(double));
    have_prev = true;
    rprev = r;
    last_extrap = false;
    if (cnt == 0) return;   // (Y, Z) = g already
    // solve (dF'dF + reg I) gamma = dF'f by Gaussian elimination, pivoting
    const int k = cnt;
    std::vector<double> M((size_t)k * (k + 1));
    double tr = 0;
    for (int a = 0; a < k; a++) tr += gram[(size_t)a * mem + a];
    const double reg = 1e-10 * tr + 1e-300;
    for (int a = 0; a < k; a++) {
      for (int b = 0; b < k; b++) M[(size_t)a * (k + 1) + b] = gram[(size_t)a * mem + b];
      M[(size_t)a * (k + 1) + a] += reg;
      M[(size_t)a * (k + 1) + k] = dot4(dF.data() + (size_t)a * N, f, N);
    }
    for (int p = 0; p < k; p++) {
      int piv = p;
      for (int a = p + 1; a < k; a++)
        if (fabs(M[(size_t)a * (k + 1) + p]) > fabs(M[(size_t)piv * (k + 1) + p])) piv = a;
      if (piv != p)
        for (int b = 0; b <= k; b++)
          std::swap(M[(size_t)p * (k + 1) + b], M[(size_t)piv * (k + 1) + b]);
      const double d = M[(size_t)p * (k + 1) + p];
      if (!(fabs(d) > 0)) return;
      for (int a = p + 1; a < k; a++) {
        const double l = M[(size_t)a * (k + 1) + p] / d;
        for (int b = p; b <= k; b++) M[(size_t)a * (k + 1) + b] -= l * M[(size_t)p * (k + 1) + b];
      }
    }
    std::vector<double> gam(k);
    for (int p = k - 1; p >= 0; p--) {
      double s = M[(size_t)p * (k + 1) + k];
      for (int b = p + 1; b < k; b++) s -= M[(size_t)p * (k + 1) + b] * gam[b];
      gam[p] = s / M[(size_t)p * (k + 1) + p];
      if (!std::isfinite(gam[p])) return;
    }
    for (int a = 0; a < k; a++) {
      const double *dg = dG.data() + (size_t)a * N;
      const double ga = gam[a];
      for (size_t t = 0; t < N; t++) g[t] -= ga * dg[t];
    }
    unpack(g, beta, Y, Z);
    last_extrap = true;
    steps++;
  }
};

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
  std::vector<float> vtop, vtmp;   // top eigenvector estimate (tight_lam)
  std::vector<double> Wprev;
  int kb, since_refresh;
  int kb_ritz = 0;            // leading block columns with Ritz values in evalsf (DNN_CHEB_COL)
  std::vector<int> coldeg;    // per-column filter degrees (DNN_CHEB_COL)
  bool have_basis;
  double lam_hi, theta_top;
  std::vector<lapack_int> eiwork;
  std::vector<double> sub, cmat, colmin;
  Hungarian hung;
  double prof[8];             // CPU seconds per phase, printed with verbose >= 2
  double sp[14];              // subspace-step split, ms, printed with verbose >= 2 (DNN_AMX etc.)
  int n_refresh, sum_kb, sum_nneg, n_sub;   // subspace eigensolver statistics
  int ref_sched, ref_weyl, ref_small, ref_fail;   // why a full refresh ran
  int n_grow;                 // block extended with fresh vectors instead of a refresh
  unsigned grow_seed;
  double delta_rel;           // last relative Cholesky shift that worked (certify)
  // Chebyshev filter and Rayleigh-Ritz products W Y on Intel AMX (amx_gemm.h):
  // opt-in with DNN_AMX=1 when the CPU has AMX-BF16; W is split into bf16 hi
  // and lo parts once per iteration, Y per product, three tile products
  // (about 16 significant bits per operand). Falls back to sgemm.
  amx_ctx amx = {};
  int use_amx = 0;
  bool amx_packed = false;
  // transform fusion (DNN_FUSE=1): the Y/Z step leaves the next forward's row phase
  // in T (t_rows_ready) and forward writes Wf alongside W (wf_ready)
  bool t_rows_ready = false, wf_ready = false;
  ~Work() { if (use_amx) amx_ctx_free(&amx); }

  void alloc(int m_)
  {
    m = m_; N = m * m; n1 = N + 1; K = (m - 1) * (m - 1) + 1;
    ll = 1 + (m - 1) * m + (m - 1);
    size_t nn = (size_t)n1 * n1, kk = (size_t)K * K;
    L.assign(nn, 0); Y.assign(nn, 0); Z.assign(nn, 0); VRV.assign(nn, 0);
    Zc.assign(nn, 0); T.assign(nn, 0);
    W.assign(kk, 0); W2.assign(kk, 0); F.assign(kk, 0); R.assign(kk, 0);
    tmp.assign(4 * (size_t)n1 + m + 8, 0);
    evals.assign(K, 0);
    // dsyevd workspace for JOBZ='V'
    ework.assign(1 + 6 * (size_t)K + 2 * kk, 0);
    eiwork.assign(3 + 5 * (size_t)K, 0);
    sub.assign((size_t)m * m, 0); cmat.assign((size_t)m * m, 0); colmin.assign(m, 0);
    hung.init(m);
    for (int t = 0; t < 8; t++) prof[t] = 0;
    for (int t = 0; t < 14; t++) sp[t] = 0;
    n_refresh = sum_kb = sum_nneg = n_sub = 0;
    ref_sched = ref_weyl = ref_small = ref_fail = 0;
    n_grow = 0; grow_seed = 12345u;
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
    vtop.assign(K, 0); vtmp.assign(K, 0);
    Wprev.assign(kk, 0);
    kb = since_refresh = 0;
    have_basis = false;
    lam_hi = theta_top = 0;
    if (const char *e = getenv("DNN_AMX"))
      if (atoi(e) > 0 && !use_amx && amx_available() && amx_ctx_init(&amx, K, K) == 0) use_amx = atoi(e);   // 1: W*Y; 2: also CholQR and Rayleigh-Ritz
    amx_packed = false;
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
// rows_done: the row phase (A + coef B, Householder block, rotation) is already in w.T
// (left there by the Y/Z step under DNN_FUSE=1). Wf_out: also write W in single
// precision, so the eigensolver's conversion pass is saved. DNN_NOSYMW=1 skips the
// final symmetrisation of W (W is symmetric to rounding by construction; the
// eigensolver reads one triangle or multiplies by the full matrix, and R is
// symmetrised again in expand).
void forward(Work &w, const double *A, double *Wout, const double *B = NULL,
             double coef = 0, bool rows_done = false, float *Wf_out = NULL)
{
  static const int nosym = getenv("DNN_NOSYMW") ? atoi(getenv("DNN_NOSYMW")) : 0;
  const int n1 = w.n1, K = w.K, ll = w.ll;
  double *T = w.T.data();
  if (!rows_done)
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
  w.t_rows_ready = false;   // T is consumed (or overwritten) by this call
  rows_D(w, T);
  rows_rot(w, T);
  if (Wf_out) {
    for (int s = 0; s < K; s++) {
      const double *row = T + (size_t)w.kept[s] * n1;
      double *wo = Wout + (size_t)s * K;
      float *wf = Wf_out + (size_t)s * K;
      for (int t = 0; t < K; t++) { const double v = row[w.kept[t]]; wo[t] = v; wf[t] = (float)v; }
    }
    w.wf_ready = !(!nosym);   // with the symmetrisation below Wf would be stale: only valid when it is skipped
  } else {
    for (int s = 0; s < K; s++) {
      const double *row = T + (size_t)w.kept[s] * n1;
      for (int t = 0; t < K; t++) Wout[(size_t)s * K + t] = row[w.kept[t]];
    }
  }
  if (!nosym) symmetrize(Wout, K);
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
// DNN_DUMP_RR=prefix (diagnostic, off by default; reads nothing back): at the
// subspace steps listed in DNN_DUMP_RR_AT (comma-separated step counts), write
// the step's inputs to prefix_<step>.bin: K, kb, degree, guard (int32), lam_hi,
// theta_top (double), theta[kb], the block X (K x kb) and W (K x K), fp32.
// Used by experiments/wil50_feasibility/speed/cs2/rr_bench.C.
static void dump_subspace_step(const char *prefix, int step, int K, int kb, int degree,
                               int guard, double lam_hi, double theta_top,
                               const float *theta, const float *X, const float *Wf)
{
  static const char *at = getenv("DNN_DUMP_RR_AT");
  if (!at) return;
  bool hit = false;
  for (const char *p = at; *p;) {
    if (atoi(p) == step) hit = true;
    while (*p && *p != ',') p++;
    if (*p) p++;
  }
  if (!hit) return;
  char fn[1024];
  snprintf(fn, sizeof fn, "%s_%d.bin", prefix, step);
  FILE *f = fopen(fn, "wb");
  if (!f) return;
  const int hdr[4] = {K, kb, degree, guard};
  const double lh[2] = {lam_hi, theta_top};
  fwrite(hdr, sizeof hdr, 1, f);
  fwrite(lh, sizeof lh, 1, f);
  fwrite(theta, sizeof(float), kb, f);
  fwrite(X, sizeof(float), (size_t)K * kb, f);
  fwrite(Wf, sizeof(float), (size_t)K * K, f);
  fclose(f);
}

bool proj_psd_subspace(Work &w, const double *Win, double *Rout, int *npos, int refresh,
                       int degree, int guard, int tight_lam)
{
  const lapack_int K = w.K;
  const size_t kk = (size_t)K * K;
  float *Wf = w.Wf.data();
  const float one = 1, zero = 0;

  // Weyl: lambda_max(W) <= lambda_max(Wprev) + ||W - Wprev||_F. With tight_lam
  // the estimate below replaces it, so the difference and the Wprev copy (two
  // passes over W, about 16 ms at K = 2117) are skipped (node_speed.md trim a).
  double tw0 = thread_cpu();
  if (w.have_basis && !(tight_lam > 0)) {
    double *d = w.F.data();   // scratch: F is rebuilt below
    for (size_t t = 0; t < kk; t++) d[t] = Win[t] - w.Wprev[t];
    w.lam_hi += sqrt(dot4(d, d, kk)) * (1 + 1e-6);
  }
  if (!(tight_lam > 0)) memcpy(w.Wprev.data(), Win, kk * sizeof(double));
  w.sp[0] += (thread_cpu() - tw0) * 1e3; tw0 = thread_cpu();
  if (!w.wf_ready) for (size_t t = 0; t < kk; t++) Wf[t] = (float)Win[t];
  w.wf_ready = false;
  w.sp[1] += (thread_cpu() - tw0) * 1e3;
  w.amx_packed = false;   // W changed: the AMX pack is rebuilt at the first product
  if (tight_lam > 0 && w.have_basis) {
    // lam_hi from shifted power steps on the top eigenvector of the last full
    // refresh, carried along. The shift by |most negative Ritz value| makes
    // lambda_max + s the dominant eigenvalue of W + s I (late in a run
    // |lambda_min| exceeds lambda_max). The Rayleigh quotient is a lower bound
    // on lambda_max, so it is padded; an underestimate only weakens the filter
    // for this iteration (the Ritz values above it are not damped), and the
    // theta_top >= lam_hi test below still forces a refresh.
    float *v = w.vtop.data(), *t = w.vtmp.data();
    const lapack_int inc = 1;
    const float sh = (float)(1.05 * std::max(0.0, -(double)w.evalsf[0]));
    double rho = 0;
    for (int it = 0; it < tight_lam; it++) {
      F77(sgemv)("N", &K, &K, &one, Wf, &K, v, &inc, &zero, t, &inc, 1);
      double vt = 0, tt = 0;
      for (int i = 0; i < K; i++) {
        t[i] += sh * v[i];
        vt += (double)v[i] * t[i]; tt += (double)t[i] * t[i];
      }
      rho = vt - sh;   // Rayleigh quotient of W (v is a unit vector)
      const double nrm = sqrt(tt) + 1e-30;
      for (int i = 0; i < K; i++) v[i] = (float)(t[i] / nrm);
    }
    const double est = rho * 1.1 + 1e-3 * fabs(rho) + 1e-6;
    w.lam_hi = std::max(est, w.theta_top + 1e-3 * (fabs(w.theta_top) + 1e-6));
  }

  bool full = !w.have_basis || w.since_refresh >= refresh || w.theta_top >= w.lam_hi;
  if (full) { if (!w.have_basis || w.since_refresh >= refresh) w.ref_sched++; else w.ref_weyl++; }
  int nneg = 0;
  lapack_int kb = w.kb;
  float *X = w.Xs.data();
  std::vector<float> &theta = w.evalsf;
  // T = W Yin (K x kb): the five products of a filter step, on AMX when enabled
  auto WY = [&](const float *Yin, float *Tout) {
    double tg0 = thread_cpu();
    if (w.use_amx) {
      if (!w.amx_packed) { amx_pack_W(&w.amx, Wf); w.amx_packed = true; w.sp[2] += (thread_cpu() - tg0) * 1e3; tg0 = thread_cpu(); }
      amx_gemm(&w.amx, Yin, (int)kb, Tout);
    } else
      F77(sgemm)("N", "N", &K, &kb, &K, &one, Wf, &K, Yin, &K, &zero, Tout, &K, 1, 1);
    w.sp[3] += (thread_cpu() - tg0) * 1e3;
  };
  // the same on the n trailing columns Yin, Tout (DNN_CHEB_COL)
  auto WYn = [&](const float *Yin, float *Tout, lapack_int n) {
    double tg0 = thread_cpu();
    if (w.use_amx) {
      if (!w.amx_packed) { amx_pack_W(&w.amx, Wf); w.amx_packed = true; w.sp[2] += (thread_cpu() - tg0) * 1e3; tg0 = thread_cpu(); }
      amx_gemm(&w.amx, Yin, (int)n, Tout);
    } else
      F77(sgemm)("N", "N", &K, &n, &K, &one, Wf, &K, Yin, &K, &zero, Tout, &K, 1, 1);
    w.sp[3] += (thread_cpu() - tg0) * 1e3;
  };
  while (true) {
    if (full) {
      lapack_int lwork = (lapack_int)w.eworkf.size(), liwork = (lapack_int)w.eiwork.size(),
                 info = 0;
      double tf0 = thread_cpu();
      memcpy(X, Wf, kk * sizeof(float));
      F77(ssyevd)("V", "L", &K, X, &K, theta.data(), w.eworkf.data(), &lwork,
                  w.eiwork.data(), &liwork, &info, 1, 1);
      w.prof[5] += thread_cpu() - tf0; w.n_refresh++;
      if (info != 0) return false;
      nneg = 0;
      while (nneg < K && theta[nneg] < 0) nneg++;
      kb = std::min<lapack_int>(K, nneg + guard);
      w.lam_hi = theta[K - 1] + 1e-6 * (fabs(theta[K - 1]) + fabs(theta[0]));
      if (tight_lam > 0) memcpy(w.vtop.data(), X + (size_t)(K - 1) * K, K * sizeof(float));
      w.since_refresh = 0;
      w.have_basis = true;
      w.kb_ritz = kb;
      break;
    }
    static const char *dump_rr = getenv("DNN_DUMP_RR");
    if (dump_rr) dump_subspace_step(dump_rr, w.n_sub, K, kb, degree, guard, w.lam_hi, w.theta_top,
                                    theta.data(), X, Wf);
    // Chebyshev filter: Y_{j+1} = 2 (W - c) Y_j / h - Y_{j-1}, Y_0 = X
    double ts0 = thread_cpu();
    float *Y = w.Ys.data(), *T = w.Ts.data(), *Yp = X;
    const float c = (float)(0.5 * (w.lam_hi + w.theta_top)),
                h = (float)(0.5 * (w.lam_hi - w.theta_top));
    const size_t nb = (size_t)K * kb;
    double tc0;
    // DNN_CHEB_COL=1 (opt-in, experiments/wil50_feasibility/speed/cs2): per-column degree
    // (projection.md idea 3). Column i of X (ascending Ritz order, value theta_i from the
    // last Rayleigh-Ritz) stops at d_i = min{d : |C_d(t_i)| >= target}, t_i = (theta_i - c)/h;
    // target = |C_degree(t_e)| at the largest negative Ritz value (the memo's rule) or
    // DNN_CHEB_COL_RHO; guard, positive and appended columns keep `degree`. d_i is made
    // non-decreasing, so step d runs on a trailing range. Columns are scaled to unit
    // norm before CholQR. Only the ADMM step changes; the certificate never sees it.
    static const int cheb_col = getenv("DNN_CHEB_COL") ? atoi(getenv("DNN_CHEB_COL")) : 0;
    static const double col_rho = getenv("DNN_CHEB_COL_RHO") ? atof(getenv("DNN_CHEB_COL_RHO")) : 0;
    if (cheb_col > 0 && degree >= 2) {
      std::vector<int> &dg = w.coldeg;
      dg.assign(kb, degree);
      int np = 0;
      while (np < w.kb_ritz && np < kb && theta[np] < 0) np++;
      if (np > 0) {
        auto cabs = [](int dd, double t) { return cosh(dd * acosh(std::max(1.0, fabs(t)))); };
        const double target = col_rho > 0 ? col_rho : cabs(degree, (theta[np - 1] - c) / h);
        for (int i = 0; i < np; i++) {
          const double ti = (theta[i] - c) / h;
          int di = 1;
          while (di < degree && cabs(di, ti) < target) di++;
          dg[i] = di;
        }
        for (int i = 1; i < kb; i++) dg[i] = std::max(dg[i], dg[i - 1]);
      }
      std::vector<lapack_int> jd(degree + 2, kb);   // jd[d] = first column with degree >= d
      for (int d = 1; d <= degree; d++) {
        lapack_int q = 0;
        while (q < kb && dg[q] < d) q++;
        jd[d] = q;
      }
      float *A = Yp, *B = Y;   // A = Y_{d-1}, B = Y_d on the active range
      WY(A, T);
      tc0 = thread_cpu();
      for (size_t t = 0; t < nb; t++) B[t] = (T[t] - c * A[t]) / h;
      w.sp[4] += (thread_cpu() - tc0) * 1e3;
      for (int d = 2; d <= degree; d++) {
        const lapack_int j0 = jd[d], n = kb - j0;
        const size_t o = (size_t)j0 * K, no = (size_t)n * K;
        WYn(B + o, T + o, n);
        tc0 = thread_cpu();
        for (size_t t = 0; t < no; t++) A[o + t] = 2 * (T[o + t] - c * B[o + t]) / h - A[o + t];
        w.sp[4] += (thread_cpu() - tc0) * 1e3;
        std::swap(A, B);
      }
      // columns that stopped at degree d were last written to the buffer B held after
      // step d; those of the other parity than `degree` move into the final B
      tc0 = thread_cpu();
      for (int d = 1; d < degree; d++)
        if ((degree - d) % 2 == 1 && jd[d + 1] > jd[d])
          memcpy(B + (size_t)jd[d] * K, A + (size_t)jd[d] * K, (size_t)(jd[d + 1] - jd[d]) * K * sizeof(float));
      for (lapack_int i = 0; i < kb; i++) {
        float *y = B + (size_t)i * K;
        double ss = 0;
        for (int t = 0; t < K; t++) ss += (double)y[t] * y[t];
        const float f = (float)(1 / sqrt(ss + 1e-300));
        for (int t = 0; t < K; t++) y[t] *= f;
      }
      w.sp[4] += (thread_cpu() - tc0) * 1e3;
      Y = B; Yp = A;
    } else {   // the uniform filter, as before
    WY(Yp, T);
    tc0 = thread_cpu();
    for (size_t t = 0; t < nb; t++) Y[t] = (T[t] - c * Yp[t]) / h;
    w.sp[4] += (thread_cpu() - tc0) * 1e3;
    for (int d = 2; d <= degree; d++) {
      WY(Y, T);
      tc0 = thread_cpu();
      for (size_t t = 0; t < nb; t++) Yp[t] = 2 * (T[t] - c * Y[t]) / h - Yp[t];
      w.sp[4] += (thread_cpu() - tc0) * 1e3;
      std::swap(Y, Yp);
    }
    }
    // Y now holds the filtered block (in Ys or Xs); orthonormalize by CholQR2
    // (DNN_CHOLQR_PASSES: 1 saves a pass, node_speed.md trim d; default 2)
    static const int cholqr_passes = getenv("DNN_CHOLQR_PASSES") ? std::max(1, atoi(getenv("DNN_CHOLQR_PASSES"))) : 2;
    float *G = w.Hs.data();
    bool ok = true;
    tc0 = thread_cpu();
    for (int pass = 0; pass < cholqr_passes && ok; pass++) {
      lapack_int info = 0;
      if (w.use_amx >= 2) {
        // G = Y'Y on AMX (symmetric, so row-major is column-major), Cholesky, then
        // Y <- Y R^-1 as an explicit inverse times Y on AMX (strtri is 1-3 ms)
        amx_gemm_tn(&w.amx, Y, Y, (int)kb, (int)kb, G);
        F77(spotrf)("U", &kb, G, &kb, &info, 1);
        if (info != 0) { ok = false; break; }
        F77(strtri)("U", "N", &kb, G, &kb, &info, 1, 1);
        if (info != 0) { ok = false; break; }
        for (int cc = 0; cc < kb; cc++)   // the strictly lower triangle (column-major) of the inverse is zero
          for (int r = cc + 1; r < kb; r++) G[(size_t)cc * kb + r] = 0;
        amx_gemm_ng(&w.amx, Y, (int)kb, G, (int)kb, (int)kb, Y);
      } else {
        F77(ssyrk)("U", "T", &kb, &K, &one, Y, &K, &zero, G, &kb, 1, 1);
        F77(spotrf)("U", &kb, G, &kb, &info, 1);
        if (info != 0) { ok = false; break; }
        F77(strsm)("R", "U", "N", "N", &K, &kb, &one, G, &kb, Y, &K, 1, 1, 1, 1);
      }
    }
    w.sp[5] += (thread_cpu() - tc0) * 1e3;
    if (!ok) { w.ref_fail++; full = true; continue; }
    // Rayleigh-Ritz: H = Y' W Y, X = Y Q
    WY(Y, T);
    tc0 = thread_cpu();
    if (w.use_amx >= 2) amx_gemm_tn(&w.amx, NULL, T, (int)kb, (int)kb, G);   // Y's tiles are the last product's
    else F77(sgemm)("T", "N", &kb, &kb, &K, &one, Y, &K, T, &K, &zero, G, &kb, 1, 1);
    for (int i = 0; i < kb; i++)
      for (int j = i + 1; j < kb; j++)
        G[(size_t)i * kb + j] = G[(size_t)j * kb + i] =
            0.5f * (G[(size_t)i * kb + j] + G[(size_t)j * kb + i]);
    w.sp[6] += (thread_cpu() - tc0) * 1e3; tc0 = thread_cpu();
    {
      lapack_int lwork = (lapack_int)w.eworkf.size(), liwork = (lapack_int)w.eiwork.size(),
                 info = 0;
      F77(ssyevd)("V", "L", &kb, G, &kb, theta.data(), w.eworkf.data(), &lwork,
                  w.eiwork.data(), &liwork, &info, 1, 1);
      if (info != 0) { w.ref_fail++; full = true; continue; }
    }
    w.sp[7] += (thread_cpu() - tc0) * 1e3; tc0 = thread_cpu();
    float *Xn = (Y == X) ? T : X;
    if (w.use_amx >= 2) amx_gemm_ng(&w.amx, Y, (int)kb, G, (int)kb, (int)kb, Xn);
    else F77(sgemm)("N", "N", &K, &kb, &kb, &one, Y, &K, G, &kb, &zero, Xn, &K, 1, 1);
    if (Xn != X) memcpy(X, Xn, nb * sizeof(float));
    w.sp[8] += (thread_cpu() - tc0) * 1e3;
    nneg = 0;
    while (nneg < kb && theta[nneg] < 0) nneg++;
    w.prof[6] += thread_cpu() - ts0;
    w.since_refresh++; w.n_sub++; w.sum_kb += kb; w.sum_nneg += nneg;
    w.theta_top = theta[kb - 1];
    if (nneg > kb - (guard + 1) / 2) {
      // block too small: the negative count grew into the guard. Rather than
      // a full refresh, append `guard` fresh pseudo-random vectors (a
      // deterministic LCG, so runs stay reproducible); the next filter and
      // CholQR2 orthonormalize them against the block and the next
      // Rayleigh-Ritz gives them Ritz values. A full ssyevd runs only when
      // the block cannot grow.
      lapack_int add = std::min<lapack_int>(guard, K - kb);
      if (add <= 0) { w.ref_small++; full = true; continue; }
      float *xn = X + (size_t)kb * K;
      for (size_t t = 0; t < (size_t)add * K; t++) {
        w.grow_seed = w.grow_seed * 1664525u + 1013904223u;
        xn[t] = (float)((w.grow_seed >> 8) & 0xffffu) / 65536.0f - 0.5f;
      }
      w.kb_ritz = kb;
      kb += add;
      w.n_grow++;
    } else {
      kb = std::min<lapack_int>(kb, nneg + guard);   // keep the lowest Ritz vectors
      w.kb_ritz = kb;
    }
    break;
  }
  w.kb = kb;
  if (full) w.theta_top = theta[kb - 1];
  *npos = K - nneg;

  // R = W - sum_{theta<0} theta x x'. DNN_R32=1 (node_speed.md trim b): the
  // rank-nneg term in fp32 (ssyrk on the fp32 eigenvectors, which is their
  // precision anyway), added to W in fp64; the default keeps dsyrk.
  static const int r32 = getenv("DNN_R32") ? atoi(getenv("DNN_R32")) : 0;
  double tr0 = thread_cpu();
  lapack_int k = nneg;
  if (r32 && k > 0) {
    float *Ff = w.Hs.data(), *Rf = w.Ts.data();   // free after the filter step (K x K each)
    for (int i = 0; i < nneg; i++) {
      const float sq = (float)sqrt(-(double)theta[i]);
      const float *x = X + (size_t)i * K;
      float *f = Ff + (size_t)i * K;
      for (int r = 0; r < K; r++) f[r] = sq * x[r];
    }
    w.sp[9] += (thread_cpu() - tr0) * 1e3; tr0 = thread_cpu();
    F77(ssyrk)("L", "N", &K, &k, &one, Ff, &K, &zero, Rf, &K, 1, 1);
    w.sp[10] += (thread_cpu() - tr0) * 1e3; tr0 = thread_cpu();
    // lower triangle (column-major) of Rf holds F F'; R = W + it, both triangles
    for (int c = 0; c < K; c++) {
      const float *rc = Rf + (size_t)c * K;
      Rout[(size_t)c * K + c] = Win[(size_t)c * K + c] + rc[c];
      for (int r = c + 1; r < K; r++) {
        const double v = Win[(size_t)c * K + r] + rc[r];
        Rout[(size_t)c * K + r] = v; Rout[(size_t)r * K + c] = v;
      }
    }
    w.sp[11] += (thread_cpu() - tr0) * 1e3;
    w.prof[7] += thread_cpu() - tr0;
    return true;
  }
  double *F = w.F.data();
  for (int i = 0; i < nneg; i++) {
    double sq = sqrt(-(double)theta[i]);
    const float *x = X + (size_t)i * K;
    double *f = F + (size_t)i * K;
    for (int r = 0; r < K; r++) f[r] = sq * x[r];
  }
  memcpy(Rout, Win, kk * sizeof(double));
  w.sp[9] += (thread_cpu() - tr0) * 1e3; tr0 = thread_cpu();
  const double done = 1;
  if (k > 0) F77(dsyrk)("L", "N", &K, &k, &done, F, &K, &done, Rout, &K, 1, 1);
  w.sp[10] += (thread_cpu() - tr0) * 1e3; tr0 = thread_cpu();
  for (int c = 0; c < K; c++)
    for (int r = c + 1; r < K; r++) Rout[(size_t)r * K + c] = Rout[(size_t)c * K + r];
  w.sp[11] += (thread_cpu() - tr0) * 1e3;
  w.prof[7] += thread_cpu() - tr0;
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
  double lam = HUGE_VAL;     // the rigorous lambda_max bound used (scaled units; DNN_CERTX)
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
  ci.lam = lam;

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
bool dnn_warm_child(const DNNWarm &par, int i, int j, DNNWarm &child, int ymode)
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
  child.Y.clear();
  if (ymode > 0) {
    // the parent's final primal iterate, for a full-state start (ymode 1:
    // restricted; ymode 2: conditioned on x_ij = 1). Both keep the child's
    // gangster zeros (they are the parent's) and the box; the child's Y-step
    // reimposes both at its first iteration anyway.
    if (par.Y.size() != (size_t)n1 * n1) return false;
    const double *Y = par.Y.data();
    auto clip = [](double v) { return v < 0 ? 0.0 : (v > 1 ? 1.0 : v); };
    child.Y.assign((size_t)nc * nc, 0);
    child.Y[0] = 1;
    if (ymode == 1) {
      for (int a = 1; a < nc; a++) {
        const int pa = map[a];
        child.Y[a] = child.Y[(size_t)a * nc] = clip(Y[pa]);
        for (int b = 1; b < nc; b++) child.Y[(size_t)a * nc + b] = clip(Y[(size_t)pa * n1 + map[b]]);
      }
    } else {
      const double xr = Y[(size_t)r * n1 + r];
      std::vector<double> q(nc, 0);
      for (int a = 1; a < nc; a++) {
        const int pa = map[a];
        q[a] = xr >= 0.05 ? clip(Y[(size_t)r * n1 + pa] / xr) : clip(Y[pa]);
        child.Y[a] = child.Y[(size_t)a * nc] = child.Y[(size_t)a * nc + a] = q[a];
      }
      for (int a = 1; a < nc; a++) {
        const double da = Y[(size_t)map[a] * n1 + map[a]];
        for (int b = 1; b < nc; b++) {
          if (b == a) continue;
          const double db = Y[(size_t)map[b] * n1 + map[b]];
          const double den = da * db;
          child.Y[(size_t)a * nc + b] = den > 1e-12
              ? clip(Y[(size_t)map[a] * n1 + map[b]] * sqrt(q[a] * q[b] / den)) : 0.0;
        }
      }
    }
  }
  return true;
}

/*
 Experimental certificate variants (experiments/wil50_feasibility/roadmap/
 certificate). Opt-in with the environment variable DNN_CERTX=<bitmask>;
 when it is unset none of this runs and no number changes. Every variant is
 a lower bound of the same kind as certify()'s, from the same iterate, with
 the same rounding allowances, and the iterates are never changed: the
 variants are printed ("certx" lines) and the first check at which each
 crosses the cutoff is reported, so that they can be compared at equal
 iterations. Bits:
   1  exact inner LAPs on the standard multiplier Zc
   2  Hahn-Grant skew polishing of M = L + scale Zc (DNN_CERTX_PASSES passes)
   4  line search on Zc + t (Z - Zc), Z the multiplier after the Z-step
      (t on a grid; lambda_max bounded by Lanczos + Rump's Cholesky test)
   8  the average of the last DNN_CERTX_AVG certificate multipliers Zc
  16  dump M = L + scale Zc at the last evaluated check (DNN_CERTX_DUMP=prefix)
  32  the same polishing with row-column reductions for the inner blocks
      (DNN_CERTX_RPASSES passes, O(m^4) each) and one exact pass at the end
  64  warm polishing: the skew matrix K accumulated by earlier checks of this
      call is added to M (GL_exact(M + K) is a certificate for any skew K),
      then DNN_CERTX_WPASSES (2) exact passes from there update K
 128  the post-step multiplier Z with its face part projected in double,
      Z - V P_psd(V'ZV) V' (one dsyevd), with the exact GL
 Evaluated at every DNN_CERTX_EVERY-th check (5), and at every check once the
 running best is within DNN_CERTX_WINDOW (80) of the cutoff, and at the last
 iteration. DNN_CERTX_CALLS=<bitmask> restricts it to those dnn_bound calls of
 the process (bit k = call k; default all). DNN_CERTX_USE=1 lets the best
 variant replace the standard value (changes the run; off by default), and
 the polish then stops as soon as it crosses the cutoff. DNN_CERTX_QUIET=1
 prints nothing (for search runs). A search configuration that uses the
 polish only near the cutoff: DNN_CERTX=32 DNN_CERTX_USE=1 DNN_CERTX_QUIET=1
 DNN_CERTX_EVERY=1000000 DNN_CERTX_WINDOW=<W> [DNN_CERTX_RPASSES=<k>].
 Modes 2, 32 and 34 alone reuse the standard certificate's lambda bound
 (no second Cholesky), as a production version would.
*/
struct CertXEnv {
  int mode = 0, passes = 8, rpasses = 12, wpasses = 2, navg = 4, every = 5, use = 0, quiet = 0;
  long calls = -1;
  double window = 80;
  std::string dump;
  CertXEnv()
  {
    if (const char *e = getenv("DNN_CERTX")) mode = atoi(e);
    if (const char *e = getenv("DNN_CERTX_PASSES")) passes = std::max(1, atoi(e));
    if (const char *e = getenv("DNN_CERTX_RPASSES")) rpasses = std::max(1, atoi(e));
    if (const char *e = getenv("DNN_CERTX_WPASSES")) wpasses = std::max(1, atoi(e));
    if (const char *e = getenv("DNN_CERTX_AVG")) navg = std::max(2, std::min(16, atoi(e)));
    if (const char *e = getenv("DNN_CERTX_EVERY")) every = std::max(1, atoi(e));
    if (const char *e = getenv("DNN_CERTX_WINDOW")) window = atof(e);
    if (const char *e = getenv("DNN_CERTX_CALLS")) calls = atol(e);
    if (const char *e = getenv("DNN_CERTX_USE")) use = atoi(e);
    if (const char *e = getenv("DNN_CERTX_DUMP")) dump = e;
    if (const char *e = getenv("DNN_CERTX_QUIET")) quiet = atoi(e);   // 1: no certx output (search runs)
  }
};
static const CertXEnv &certx_env() { static CertXEnv e; return e; }

// Gilmore-Lawler value of a full lifted matrix M (unscaled, n1 x n1, rows
// used as in certify), and the magnitude sum for certify's margin formula.
static long double gl_full(Work &w, const double *M, bool exact, long double &absM)
{
  const int m = w.m, n1 = w.n1;
  double *sub = w.sub.data(), *cm = w.cmat.data();
  absM = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) {
      const int r = 1 + i * m + j;
      const double *mr = M + (size_t)r * n1;
      int kk = 0;
      for (int k = 0; k < m; k++) {
        if (k == i) continue;
        int cb = 1 + k * m, ss = 0;
        for (int l = 0; l < m; l++) {
          if (l == j) continue;
          const double x = mr[cb + l];
          sub[kk * (m - 1) + ss++] = x;
          absM += fabs(x);
        }
        kk++;
      }
      long double inner = exact ? w.hung.solve(m - 1, sub, m - 1)
                                : lap_reduce(m - 1, sub, m - 1, w.colmin.data());
      const double d0 = mr[0], d0t = M[r], dd = mr[r];   // M_(ij),0, M_0,(ij), diagonal
      absM += fabs(d0) + fabs(d0t) + fabs(dd) + fabsl(inner);
      cm[i * m + j] = (double)((long double)d0 + d0t + dd + inner);
    }
  long double val = M[0] + w.hung.solve(m, cm, m);
  for (int t = 0; t < m * m; t++) absM += fabs(cm[t]);
  absM += fabs(M[0]);
  return val;
}

/*
 Hahn-Grant / Adams-Johnson dual ascent on M (no eliminations). R starts as
 M. Each pass: exact LAP on every block row (ij) of R (k != i, l != j),
 leaving the reduced block (>= 0); outer LAP on c_ij = R_0,ij + R_ij,0 +
 R_ij,ij + LAP_ij; the pass adds R_00 + LAP(c) to the bound B and leaves
 R_00 = linear = diagonal = 0; the reduced outer costs are spread evenly over
 their block row (a lifted permutation with x_ij = 1 uses exactly m-1 entries
 of it); then each non-gangster pair R_ab, R_ba is replaced by its average
 (a skew-symmetric change, which no symmetric Y sees). For every lifted
 permutation y, <M, yy'> = B + <R, yy'> >= B in exact arithmetic of the
 computed numbers; err bounds the rounding of the passes. vals[p] = B after
 pass p+1 (the first equals the exact GL value), errs[p] the allowance.
*/
static void hg_polish(Work &w, const double *M, int passes, std::vector<double> &R,
                      std::vector<long double> &vals, std::vector<double> &errs,
                      long double *absM1 = NULL, int exact_from = 0, double *Kacc = NULL,
                      double stop_val = HUGE_VAL)
{
  // passes p >= exact_from solve the inner LAPs exactly; earlier ones use the
  // row-then-column reduction (also dual feasible, O(m^2) a block).
  // Kacc (optional, n1 x n1): accumulates the skew-symmetric part of the
  // transfers. With exact LAPs, GL is invariant to the other changes (the
  // reductions and the spreading are face terms of the form GL's own LAPs
  // absorb), so GL_exact(M + Kacc) >= the value reached here; a later check
  // can start from M' + Kacc (bit 64 of DNN_CERTX).
  const int m = w.m, n1 = w.n1, m1 = m - 1;
  R.assign(M, M + (size_t)n1 * n1);
  vals.clear(); errs.clear();
  long double absM = 0;   // certify's magnitude sum for pass 1 (= the exact GL of M)
  std::vector<double> c((size_t)m * m), u(m1), v(m1), U(m), Vv(m);
  double *sub = w.sub.data();
  long double B = 0;
  double err = 0;
  for (int p = 0; p < passes; p++) {
    double maxR = 0, maxd = 0, maxRn = 0, maxc = 0, maxUV = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        const int r = 1 + i * m + j;
        double *rr = R.data() + (size_t)r * n1;
        int kk = 0;
        for (int k = 0; k < m; k++) {
          if (k == i) continue;
          int cb = 1 + k * m, ss = 0;
          for (int l = 0; l < m; l++) {
            if (l == j) continue;
            sub[kk * m1 + ss++] = rr[cb + l];
            maxR = std::max(maxR, fabs(rr[cb + l]));
            if (p == 0) absM += fabs(rr[cb + l]);
          }
          kk++;
        }
        const bool ex = p >= exact_from;
        long double lv = 0;
        if (ex) {
          lv = w.hung.solve(m1, sub, m1);   // = sum u + sum v (repaired)
          for (int a = 0; a < m1; a++) u[a] = w.hung.u[a + 1];
        } else {
          // row then column reduction: u = row minima, v below; lv = sum u + sum v
          for (int a = 0; a < m1; a++) {
            double mn = HUGE_VAL;
            for (int b = 0; b < m1; b++) mn = std::min(mn, sub[a * m1 + b]);
            u[a] = mn;
          }
        }
        for (int a = 0; a < m1; a++) maxd = std::max(maxd, fabs(u[a]));
        for (int b = 0; b < m1; b++) {
          double mn = HUGE_VAL;
          for (int a = 0; a < m1; a++) mn = std::min(mn, sub[a * m1 + b] - u[a]);
          v[b] = mn; maxd = std::max(maxd, fabs(mn));
        }
        if (!ex) {
          for (int a = 0; a < m1; a++) lv += u[a];
          for (int b = 0; b < m1; b++) lv += v[b];
        }
        if (p == 0) absM += fabs(rr[0]) + fabs(R[r]) + fabs(rr[r]) + fabsl(lv);
        kk = 0;
        for (int k = 0; k < m; k++) {
          if (k == i) continue;
          int cb = 1 + k * m, ss = 0;
          for (int l = 0; l < m; l++) {
            if (l == j) continue;
            rr[cb + l] = (sub[kk * m1 + ss] - u[kk]) - v[ss];   // >= 0
            ss++;
          }
          kk++;
        }
        c[i * m + j] = (double)((long double)rr[0] + R[r] + rr[r] + lv);
        maxc = std::max(maxc, fabs(c[i * m + j]));
      }
    long double ov = w.hung.solve(m, c.data(), m);
    if (p == 0) {
      for (int t = 0; t < m * m; t++) absM += fabs(c[t]);
      absM += fabs(R[0]);
      if (absM1) *absM1 = absM;
    }
    for (int a = 0; a < m; a++) { U[a] = w.hung.u[a + 1]; maxUV = std::max(maxUV, fabs(U[a])); }
    for (int b = 0; b < m; b++) {
      double mn = HUGE_VAL;
      for (int a = 0; a < m; a++) mn = std::min(mn, c[a * m + b] - U[a]);
      Vv[b] = mn; maxUV = std::max(maxUV, fabs(mn));
    }
    B += (long double)R[0] + ov;
    R[0] = 0;
    for (int r = 1; r < n1; r++) { R[r] = 0; R[(size_t)r * n1] = 0; R[(size_t)r * n1 + r] = 0; }
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        const double cr = (c[i * m + j] - U[i]) - Vv[j];   // >= 0
        const double s = cr / m1;
        if (s == 0) continue;
        double *rr = R.data() + (size_t)(1 + i * m + j) * n1;
        for (int k = 0; k < m; k++) {
          if (k == i) continue;
          for (int l = 0; l < m; l++) if (l != j) rr[1 + k * m + l] += s;
        }
      }
    // average each non-gangster pair (the gangster entries are never read)
    for (int a = 1; a < n1; a++) {
      const int i = (a - 1) / m, j = (a - 1) % m;
      double *ra = R.data() + (size_t)a * n1;
      for (int b = a + 1; b < n1; b++) {
        const int k = (b - 1) / m, l = (b - 1) % m;
        if (k == i || l == j) continue;
        double &x = ra[b], &y = R[(size_t)b * n1 + a];
        const double s = 0.5 * (x + y);
        if (Kacc) {
          // the skew part of this pass: K_ab += (R_ba - R_ab)/2, K_ba = -K_ab exactly
          double &ka = Kacc[(size_t)a * n1 + b];
          ka += 0.5 * (y - x);
          Kacc[(size_t)b * n1 + a] = -ka;
        }
        x = y = s;
        maxRn = std::max(maxRn, fabs(s));
      }
    }
    // rounding of this pass for one lifted permutation: m(m-1) block entries
    // (reduction 2 roundings, spread 1, average 2, each relative to operands
    // bounded by maxR, maxd, maxRn) and m linear terms (c, outer reduction),
    // doubled for safety
    err += 2.0 * U_ROUND * ((double)m * m1 * (2 * maxR + 4 * maxd + 6 * maxRn) +
                            (double)m * (8 * maxc + 4 * maxUV) + 4 * fabsl(B));
    vals.push_back(B);
    errs.push_back(err);
    // early exit once the value net of this pass's margins clears stop_val
    // (the caller's cutoff less the margins that do not depend on the passes)
    if ((double)B - err - 8.0 * (m + 2) * U_ROUND * (double)absM > stop_val) break;
  }
}

// lambda_max(S) <= delta + c by Rump's Cholesky test on delta I - S (as in
// certify); HUGE_VAL when the factorization fails. Clobbers w.W2.
static double rump_bound(Work &w, const double *S, double delta)
{
  const lapack_int K = w.K;
  double *A = w.W2.data();
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
  if (info != 0) return HUGE_VAL;
  double g = (K + 1) * U_ROUND / (1 - (K + 1) * U_ROUND);
  double trA = (double)tr * (1 + 1e-12) + 1e-300;
  double c = (g / (1 - g)) * trA / (1 - K * g / (1 - g)) * 1.25
             + 4.0 * K * (2.0 * (K + 1) + maxdiag) * DBL_TRUE_MIN;
  return delta + c + U_ROUND * maxdiag;
}

// Largest eigenvalue estimate of symmetric S (K x K) by Lanczos with full
// reorthogonalisation (steps), and bisection on the tridiagonal matrix.
static double lanczos_top(int K, const double *S, int steps)
{
  steps = std::min(steps, K);
  std::vector<double> Q((size_t)(steps + 1) * K), al(steps), be(steps + 1, 0), wv(K);
  unsigned s = 987654321u;
  double nrm = 0;
  for (int i = 0; i < K; i++) { s = s * 1664525u + 1013904223u; Q[i] = (double)(s >> 8) / (1 << 24) - 0.5; nrm += Q[i] * Q[i]; }
  nrm = sqrt(nrm);
  for (int i = 0; i < K; i++) Q[i] /= nrm;
  int nst = 0;
  for (int k = 0; k < steps; k++) {
    const double *q = Q.data() + (size_t)k * K;
    for (int i = 0; i < K; i++) wv[i] = dot4(S + (size_t)i * K, q, K);
    al[k] = dot4(wv.data(), q, K);
    for (int r = 0; r < 2; r++)   // full reorthogonalisation, twice
      for (int t = 0; t <= k; t++) {
        const double *qt = Q.data() + (size_t)t * K;
        double h = dot4(wv.data(), qt, K);
        for (int i = 0; i < K; i++) wv[i] -= h * qt[i];
      }
    double b = sqrt(dot4(wv.data(), wv.data(), K));
    nst = k + 1;
    if (b < 1e-14 * (fabs(al[k]) + 1e-300)) break;
    be[k + 1] = b;
    double *qn = Q.data() + (size_t)(k + 1) * K;
    for (int i = 0; i < K; i++) qn[i] = wv[i] / b;
  }
  // largest eigenvalue of the nst x nst tridiagonal (al, be[1..]) by bisection
  double lo = HUGE_VAL, hi = -HUGE_VAL;
  for (int k = 0; k < nst; k++) {
    double r = (k > 0 ? fabs(be[k]) : 0) + (k + 1 < nst ? fabs(be[k + 1]) : 0);
    lo = std::min(lo, al[k] - r); hi = std::max(hi, al[k] + r);
  }
  for (int it = 0; it < 200 && hi - lo > 1e-15 * (fabs(lo) + fabs(hi)); it++) {
    double x = 0.5 * (lo + hi), d = 1;
    int cnt = 0;   // eigenvalues above x
    for (int k = 0; k < nst; k++) {
      d = al[k] - x - (k > 0 ? be[k] * be[k] / d : 0);
      if (d == 0) d = 1e-300;
      if (d > 0) cnt++;
    }
    if (cnt > 0) lo = x; else hi = x;
  }
  return hi;
}

// a rigorous lambda_max bound for V'ZmV (scaled units): Lanczos estimate, then
// Rump's test at growing shifts above it; HUGE_VAL if every try fails
static double cert_lambda(Work &w, const double *Zm)
{
  const int m = w.m;
  const int K = w.K;
  double *S = w.W.data();
  forward(w, Zm, S);
  const double znorm = sqrt(dot4(Zm, Zm, (size_t)w.n1 * w.n1)) * (1 + 1e-8);
  const double snorm = sqrt(dot4(S, S, (size_t)K * K)) * (1 + 1e-8);
  const double th = lanczos_top(K, S, 40);
  double gap = std::max(1e-3 * fabs(th), 1e-9 * snorm) + 1e-14 * znorm + 1e-300;
  for (int tries = 0; tries < 8; tries++, gap *= 8) {
    double lam = rump_bound(w, S, std::max(th, 0.0) + gap);
    if (lam < HUGE_VAL) return lam + 64.0 * m * U_ROUND * znorm;
  }
  return HUGE_VAL;
}

static std::atomic<int> dnn_call_counter(0);   // DNN_CERTX_CALLS: index of this dnn_bound call in the process

DNNResult dnn_bound(const MAT *A, const MAT *B, const MAT *C, double shift,
                    double cutoff, const DNNParams &par, const DNNWarm *start,
                    DNNWarm *save)
{
  DNNResult res;
  const int call_idx = dnn_call_counter++;
  res.lb = -HUGE_VAL; res.primal = HUGE_VAL; res.status = DNN_ERROR;
  res.iterations = 0; res.cpu = 0; res.certificates = 0; res.cert_failures = 0;
  res.refinements = 0;
  res.halpern_restarts = 0;
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
  const double gamma = par.gamma, gamma1 = par.gamma1;
  const int check_every = std::max(1, par.check_every);

  double *Y = w.Y.data(), *Z = w.Z.data(), *VRV = w.VRV.data(), *Zc = w.Zc.data();
  Y[0] = 1.0;
  if (par.warm_primal && start) {
    // Explicit same-node initialization; never infer it from dimension alone.
    const size_t count = (size_t)n1 * n1;
    if (start->m != m || start->Y.size() != count || start->Z.size() != count ||
        !std::isfinite(start->scale) || start->scale <= 0) {
      res.status = DNN_ERROR; res.cpu = thread_cpu() - t0; return res;
    }
    for (size_t t = 0; t < count; t++)
      if (!std::isfinite(start->Y[t]) || !std::isfinite(start->Z[t])) {
        res.status = DNN_ERROR; res.cpu = thread_cpu() - t0; return res;
      }
    std::copy(start->Y.begin(), start->Y.end(), Y);
  }
  if (start && start->m == m && start->Z.size() == (size_t)n1 * n1 && start->scale > 0) {
    // multipliers are kept in unscaled units across nodes
    const double f = start->scale / scale;
    for (size_t t = 0; t < (size_t)n1 * n1; t++) {
      Z[t] = start->Z[t] * f;
      if (par.warm_primal && !std::isfinite(Z[t])) {
        res.status = DNN_ERROR; res.cpu = thread_cpu() - t0; return res;
      }
    }
  }
  double best = -HUGE_VAL;
  std::vector<double> red_cur, red_best;   // child bounds (par.child_bounds)
  double lbs_best = -HUGE_VAL;
  std::vector<double> lbhist;   // certified lb at each check (for the stall test)
  std::vector<double> besthist; // running best at each check (for the hopeless test)
  int confirm_at = -1;          // par.hopeless_zone: a deferred give-up is confirmed from here
  int status = -1, npos = 0;
  double obj = HUGE_VAL;
  const double gap = par.auto_gap > 0 ? par.auto_gap : 1e-3 * std::max(1.0, fabs(cutoff));
  // no finite cutoff: never close to it, so never the exact pass (was NaN)
  const double exact_above = std::isfinite(cutoff) ? cutoff - gap - shift : HUGE_VAL;

  // Halpern anchoring (par.halpern)
  std::vector<double> Y0, Z0;
  int h_k = 0;
  double h_r0 = -1;
  Anderson aa;
  if (par.anderson > 0) aa.init(par.anderson, n1);
  if (par.halpern > 0) {
    Y0.assign(Y, Y + (size_t)n1 * n1);
    Z0.assign(Z, Z + (size_t)n1 * n1);
  }
  // momentum (par.momentum): Y, Z hold the point the next R-step consumes
  // (extrapolated), Yt, Zt the last true iterates
  bool m_on = par.momentum > 0;
  std::vector<double> Yt, Zt;
  double m_a = 1, m_c = -1, m_cap = par.momentum_cap > 0 ? par.momentum_cap : 1.0;
  int m_drops = 0;
  if (m_on) {
    Yt.assign(Y, Y + (size_t)n1 * n1);
    Zt.assign(Z, Z + (size_t)n1 * n1);
  }

  // DNN_CERTX (experimental certificate variants, see above): off unless set
  const CertXEnv &cxe = certx_env();
  const bool cx_on = cxe.mode != 0 && !par.elim && par.poly == DNN_POLY_GL && !par.prototype_cert &&
                     (cxe.calls < 0 || (call_idx < 62 && ((cxe.calls >> call_idx) & 1)));
  std::vector<std::vector<double>> cx_ring;   // the last navg certificate multipliers (bit 8)
  int cx_ring_n = 0, cx_ring_pos = 0;
  std::vector<double> cx_M, cx_M2, cx_M3, cx_Zt, cx_R, cx_K;
  // std, exact, hg, ls, avg, hgr, hgw, proj: first check above the cutoff, best value
  int cx_first[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
  double cx_best[8] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  double cx_cpu = 0;
  if (cx_on && !cxe.quiet) printf("certx call %d m %d mode %d passes %d rpasses %d avg %d every %d window %g use %d\n", call_idx,
                    m, cxe.mode, cxe.passes, cxe.rpasses, cxe.navg, cxe.every, cxe.window, cxe.use);

  // DNN_FUSE=1: the Y/Z step writes the next forward's row phase into T while the row
  // is in cache, and forward writes Wf with W (node_speed.md trim f); not with the
  // accelerations, which rewrite Y and Z after the step
  static const int fuse_env = getenv("DNN_FUSE") ? atoi(getenv("DNN_FUSE")) : 0;
  const bool fuse_ok = fuse_env > 0 && par.anderson == 0 && par.halpern == 0 && !m_on;
  // subspace block guard = K / DNN_GUARD_DIV (20; node_speed.md trim c tries 40)
  static const int guard_div = getenv("DNN_GUARD_DIV") ? std::max(1, atoi(getenv("DNN_GUARD_DIV"))) : 20;
  for (int it = 1; it <= par.max_iter; it++) {
    // R-step: W = V'(Y + Z/beta)V, R = P_psd(W), VRV = V R V'
    double tp = thread_cpu(), tq;
    if (par.anderson > 0) aa.begin(Y, Z, ibeta);
    forward(w, Y, w.W.data(), Z, ibeta, w.t_rows_ready,
            par.eig_mode == DNN_EIG_SUBSPACE ? w.Wf.data() : NULL);   // Wf with W (used when DNN_NOSYMW=1)
    tq = thread_cpu(); w.prof[0] += tq - tp; tp = tq;
    bool check = it % check_every == 0 || it == par.max_iter;
    bool rstep_cert = check && !par.prototype_cert;
    bool okp = par.eig_mode == DNN_EIG_SUBSPACE
                   ? proj_psd_subspace(w, w.W.data(), w.R.data(), &npos, par.refresh,
                                       par.cheb_degree, std::max(8, w.K / guard_div), par.tight_lam)
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
    long double pres2 = 0, dy2 = 0, dz2 = 0;
    double *yn = w.tmp.data();
    for (int r = 0; r < n1; r++) {
      const double *lr = L + (size_t)r * n1, *vr = VRV + (size_t)r * n1;
      double *yr = Y + (size_t)r * n1, *zr = Z + (size_t)r * n1, *cr = Zc + (size_t)r * n1;
      if (rstep_cert)
        for (int c = 0; c < n1; c++) cr[c] = zr[c] + beta * (yr[c] - vr[c]);
      if (gamma1 > 0)
        for (int c = 0; c < n1; c++) zr[c] += gamma1 * beta * (yr[c] - vr[c]);
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
      if (fuse_ok) {
        // the next R-step's forward, row phase: (Y + Z/beta) row, Householder block, rotation
        double *trow = w.T.data() + (size_t)r * n1;
        for (int c = 0; c < n1; c++) trow[c] = yr[c] + ibeta * zr[c];
        hh_block(w, trow + 1);   // scratch w.tmp[0..m): yn is finished with
        col_rot(trow, w.ll);
      }
      if (par.halpern > 0 || m_on) {
        // the Z change of this iteration, in units of Y: gamma1 (Y_old - VRV) +
        // gamma (Y_new - VRV)
        double *dzs = ds + n1;
        for (int c = 0; c < n1; c++)
          dzs[c] = gamma * rs[c] + (gamma1 > 0 ? gamma1 * (rs[c] - ds[c]) : 0);
        dz2 += dot4(dzs, dzs, n1);
      }
    }
    if (fuse_ok) w.t_rows_ready = true;
    if (par.anderson > 0) {
      aa.step(Y, Z, beta);   // timed with the Y/Z-step
      res.halpern_restarts = aa.resets;
    } else if (par.halpern > 0) {
      const double r = sqrt((double)(dy2 + dz2));
      if (h_r0 < 0) h_r0 = r;
      if (h_k >= par.halpern || r <= 0.2 * h_r0) {
        // restart: the anchor is the current point
        memcpy(Y0.data(), Y, (size_t)n1 * n1 * sizeof(double));
        memcpy(Z0.data(), Z, (size_t)n1 * n1 * sizeof(double));
        h_k = 0; h_r0 = -1; res.halpern_restarts++;
      } else {
        const double lam = 1.0 / (h_k + 2);
        for (size_t t = 0; t < (size_t)n1 * n1; t++) {
          Y[t] = lam * Y0[t] + (1 - lam) * Y[t];
          Z[t] = lam * Z0[t] + (1 - lam) * Z[t];
        }
        h_k++;
      }
    } else if (m_on) {
      // fast ADMM: the combined residual decides between a momentum step
      // and a restart at the current (true) point
      const double c = (double)(dy2 + dz2);
      if (m_c < 0 || c < par.momentum * m_c) {
        const double an = 0.5 * (1 + sqrt(1 + 4 * m_a * m_a));
        const double om = std::min((m_a - 1) / an, m_cap);
        m_a = an; m_c = c;
        for (size_t t = 0; t < (size_t)n1 * n1; t++) {
          const double yt = Y[t], zt = Z[t];
          Y[t] = yt + om * (yt - Yt[t]); Yt[t] = yt;
          Z[t] = zt + om * (zt - Zt[t]); Zt[t] = zt;
        }
      } else {
        m_a = 1; m_c = c; res.halpern_restarts++;
        memcpy(Yt.data(), Y, (size_t)n1 * n1 * sizeof(double));
        memcpy(Zt.data(), Z, (size_t)n1 * n1 * sizeof(double));
      }
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
    if (cx_on) {
      // experimental certificate variants at this iterate (DNN_CERTX; see above)
      const double tx0 = thread_cpu();
      const size_t nn = (size_t)n1 * n1;
      if (cxe.mode & 8) {
        if ((int)cx_ring.size() < cxe.navg) cx_ring.emplace_back(nn);
        std::copy(Zc, Zc + nn, cx_ring[cx_ring_pos].begin());
        cx_ring_pos = (cx_ring_pos + 1) % cxe.navg;
        cx_ring_n = std::min(cx_ring_n + 1, cxe.navg);
      }
      const double run_best = std::max(best, lbs > -HUGE_VAL ? lbs + shift : -HUGE_VAL);
      const bool due = res.certificates % cxe.every == 0 || it == par.max_iter ||
                       (std::isfinite(cutoff) && run_best > -HUGE_VAL && cutoff - run_best <= cxe.window);
      if (due) {
        const double dr_save = w.delta_rel;
        double u_[8] = {lbs, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};   // unshifted
        double tm[7] = {0, 0, 0, 0, 0, 0, 0}, t_best = 0;
        std::vector<long double> hv, hv2, hv3;
        std::vector<double> he, he2, he3;
        double mg_hg = 0, mg_hgr = 0, mg_hgw = 0;
        CertInfo cx;
        double tq = thread_cpu();
        // the polishings alone (modes 2, 32, 64 and their sums): reuse the standard
        // certificate's lambda bound (the same Zc), so they cost only their LAP passes
        const bool lean = (cxe.mode & ~(2 | 32 | 64)) == 0;
        if (lean) cx = ci;
        else u_[1] = certify(w, Zc, scale, (double)sumabsL, par.poly, DNN_LAP_EXACT, exact_above, exact_below, cx);
        tm[0] = thread_cpu() - tq;
        if ((cxe.mode & 64) && cx.chol_ok) {
          tq = thread_cpu();
          if (cx_K.empty()) cx_K.assign(nn, 0.0);
          cx_M3.resize(nn);
          double mx3 = 0;
          for (size_t t = 0; t < nn; t++) {
            cx_M3[t] = (L[t] + scale * Zc[t]) + cx_K[t];
            mx3 = std::max(mx3, fabs(cx_M3[t]));
          }
          long double absW = 0;
          static const int wtest = getenv("DNN_CERTX_WTEST") ? atoi(getenv("DNN_CERTX_WTEST")) : 0;
          if (wtest) {
            // check of the invariance claim at one check: K from a cold polish of M,
            // then GL_exact(M + K) against the polish's own value
            std::fill(cx_K.begin(), cx_K.end(), 0.0);
            for (size_t t = 0; t < nn; t++) cx_M3[t] = L[t] + scale * Zc[t];
            std::vector<long double> hvt; std::vector<double> het;
            hg_polish(w, cx_M3.data(), cxe.wpasses, cx_R, hvt, het, NULL, 0, cx_K.data());
            for (size_t t = 0; t < nn; t++) cx_M3[t] += cx_K[t];
            long double ab;
            const long double g = gl_full(w, cx_M3.data(), true, ab);
            printf("certx-wtest it %d polish %.6f GL_exact(M+K) %.6f diff %.3e\n", it,
                   (double)hvt.back() + shift, (double)g + shift, (double)(g - hvt.back()));
            std::fill(cx_K.begin(), cx_K.end(), 0.0);
            for (size_t t = 0; t < nn; t++) cx_M3[t] = L[t] + scale * Zc[t];
          }
          hg_polish(w, cx_M3.data(), cxe.wpasses, cx_R, hv3, he3, &absW, 0, cx_K.data());
          // certify's margin with pass 1's magnitudes (M + K formed), the second
          // rounding of forming M + K for the (m+1)^2 entries a permutation reads
          mg_hgw = 8.0 * (m + 2) * U_ROUND * (double)absW + (m + 2) * U_ROUND * (double)sumabsL +
                   cx.sdp_penalty + 2.0 * U_ROUND * (m + 1) * (m + 1) * mx3;
          for (size_t p = 0; p < hv3.size(); p++) {
            const double v = (double)hv3[p] - he3[p] - mg_hgw;
            if (std::isfinite(v)) u_[6] = std::max(u_[6], v);
          }
          tm[5] = thread_cpu() - tq;
        }
        if (cx.chol_ok && (cxe.mode & (2 | 16 | 32))) {
          cx_M.resize(nn);
          for (size_t t = 0; t < nn; t++) cx_M[t] = L[t] + scale * Zc[t];
        }
        if ((cxe.mode & 2) && cx.chol_ok) {
          tq = thread_cpu();
          long double absM1 = 0;
          // with DNN_CERTX_USE, stop as soon as the cutoff is crossed
          const double stop_val = (cxe.use && std::isfinite(cutoff))
            ? cutoff - shift + (m + 2) * U_ROUND * (double)sumabsL + cx.sdp_penalty + 1e-9 * (1 + fabs(cutoff))
            : HUGE_VAL;
          hg_polish(w, cx_M.data(), cxe.passes, cx_R, hv, he, &absM1, 0, NULL, stop_val);
          // pass 1 is the exact GL value of the same M: certify's margin formula
          // (forming M, the GL sums, L) with pass 1's magnitudes, the semidefinite
          // term, plus the passes' own rounding (he)
          mg_hg = 8.0 * (m + 2) * U_ROUND * (double)absM1 + (m + 2) * U_ROUND * (double)sumabsL +
                  cx.sdp_penalty;
          for (size_t p = 0; p < hv.size(); p++) {
            const double v = (double)hv[p] - he[p] - mg_hg;
            if (std::isfinite(v)) u_[2] = std::max(u_[2], v);
          }
          tm[1] = thread_cpu() - tq;
        }
        if ((cxe.mode & 32) && cx.chol_ok) {
          tq = thread_cpu();
          long double absR = 0;
          const double stop_val = (cxe.use && std::isfinite(cutoff))
            ? cutoff - shift + (m + 2) * U_ROUND * (double)sumabsL + cx.sdp_penalty + 1e-9 * (1 + fabs(cutoff))
            : HUGE_VAL;
          hg_polish(w, cx_M.data(), cxe.rpasses + 1, cx_R, hv2, he2, &absR, cxe.rpasses, NULL, stop_val);
          mg_hgr = 8.0 * (m + 2) * U_ROUND * (double)absR + (m + 2) * U_ROUND * (double)sumabsL +
                   cx.sdp_penalty;
          for (size_t p = 0; p < hv2.size(); p++) {
            const double v = (double)hv2[p] - he2[p] - mg_hgr;
            if (std::isfinite(v)) u_[5] = std::max(u_[5], v);
          }
          tm[4] = thread_cpu() - tq;
        }
        if (cxe.mode & 4) {
          tq = thread_cpu();
          static const double grid[] = {0.01, 0.03, 0.1, 0.3, 1.0};
          cx_Zt.resize(nn); cx_M2.resize(nn);
          const double *Zn = Z;   // the multiplier after this iteration's Z-step
          for (double tt : grid) {
            for (size_t q = 0; q < nn; q++) cx_Zt[q] = Zc[q] + tt * (Zn[q] - Zc[q]);
            const double lam = cert_lambda(w, cx_Zt.data());
            if (lam == HUGE_VAL) continue;
            for (size_t q = 0; q < nn; q++) cx_M2[q] = L[q] + scale * cx_Zt[q];
            long double absM;
            const long double val = gl_full(w, cx_M2.data(), true, absM);
            const double v = (double)val - (8.0 * (m + 2) * U_ROUND * (double)absM +
                                            (m + 2) * U_ROUND * (double)sumabsL + (m + 1) * scale * lam);
            if (par.verbose >= 2)
              printf("certx-ls it %d t %.3f value %.4f poly %.4f sdp %.4f\n", it, tt, v + shift,
                     (double)val + shift, (m + 1) * scale * lam);
            if (std::isfinite(v) && v > u_[3]) { u_[3] = v; t_best = tt; }
          }
          tm[2] = thread_cpu() - tq;
        }
        if ((cxe.mode & 8) && cx_ring_n >= 2) {
          tq = thread_cpu();
          cx_Zt.resize(nn);
          std::fill(cx_Zt.begin(), cx_Zt.end(), 0.0);
          for (int r = 0; r < cx_ring_n; r++) {
            const double *zr = cx_ring[r].data();
            for (size_t q = 0; q < nn; q++) cx_Zt[q] += zr[q];
          }
          const double f = 1.0 / cx_ring_n;
          for (size_t q = 0; q < nn; q++) cx_Zt[q] *= f;
          CertInfo ca;
          u_[4] = certify(w, cx_Zt.data(), scale, (double)sumabsL, par.poly, DNN_LAP_EXACT, exact_above,
                          exact_below, ca);
          tm[3] = thread_cpu() - tq;
        }
        if (cxe.mode & 128) {
          // the post-step multiplier Z (box side exact) with its face part projected
          // onto the negative semidefinite cone in double: Zp = Z - V P_psd(V'ZV) V'
          // (the prototype's certificate, here with the exact GL); w.R and VRV are
          // free until the next R-step, as in the refinement above
          tq = thread_cpu();
          cx_Zt.resize(nn);
          forward(w, Z, w.W.data());
          int np2;
          if (proj_psd(w, w.W.data(), w.R.data(), &np2)) {
            expand(w, w.R.data(), VRV);
            for (size_t q = 0; q < nn; q++) cx_Zt[q] = Z[q] - VRV[q];
            CertInfo cp2;
            u_[7] = certify(w, cx_Zt.data(), scale, (double)sumabsL, par.poly, DNN_LAP_EXACT, exact_above,
                            exact_below, cp2);
          }
          tm[6] = thread_cpu() - tq;
        }
        if ((cxe.mode & 16) && !cxe.dump.empty() && cx.chol_ok) {
          // M = L + scale Zc and what is needed to turn a polyhedral value into a bound
          std::string pm = cxe.dump + ".M.bin", pt = cxe.dump + ".meta";
          if (FILE *f = fopen(pm.c_str(), "wb")) { fwrite(cx_M.data(), sizeof(double), nn, f); fclose(f); }
          if (FILE *f = fopen(pt.c_str(), "w")) {
            fprintf(f, "m %d\nit %d\nscale %.17g\nshift %.17g\ncutoff %.17g\nmargin %.17g\nsdp_penalty %.17g\n"
                    "exact %.17g\nstd %.17g\nsumabsL %.17g\n", m, it, scale, shift, cutoff, cx.margin,
                    cx.sdp_penalty, u_[1], lbs, (double)sumabsL);
            fclose(f);
          }
        }
        w.delta_rel = dr_save;
        if (!cxe.quiet) {
        const double primal = dot4(L, Y, nn) + shift;
        printf("certx it %5d std %.4f exact %.4f hg", it, u_[0] + shift, u_[1] + shift);
        for (size_t p = 0; p < hv.size(); p++)
          if (p == 0 || p == 1 || p == 3 || p + 1 == hv.size())
            printf(" %zu:%.4f", p + 1, (double)hv[p] - he[p] - mg_hg + shift);
        printf(" ls %.4f (t %.2f) avg%d %.4f primal %.4f cutoff %.4f cpu %.2f/%.2f/%.2f/%.2f",
               u_[3] + shift, t_best, cx_ring_n, u_[4] + shift, primal, cutoff, tm[0], tm[1], tm[2], tm[3]);
        if (!hv2.empty()) {
          // reduce passes 1, 2, 4, the last, and the final exact pass
          printf(" hgr");
          for (size_t p = 0; p < hv2.size(); p++)
            if (p == 0 || p == 1 || p == 3 || p + 2 == hv2.size() || p + 1 == hv2.size())
              printf(" %zu:%.4f", p + 1, (double)hv2[p] - he2[p] - mg_hgr + shift);
          printf(" cpu %.2f", tm[4]);
        }
        if (!hv3.empty()) {
          printf(" hgw");
          for (size_t p = 0; p < hv3.size(); p++)
            printf(" %zu:%.4f", p + 1, (double)hv3[p] - he3[p] - mg_hgw + shift);
          printf(" cpuw %.2f", tm[5]);
        }
        if (cxe.mode & 128) printf(" proj %.4f cpup %.2f", u_[7] + shift, tm[6]);
        printf("\n");
        }
        double vmax = -HUGE_VAL;
        for (int v = 0; v < 8; v++) {
          if (u_[v] > -HUGE_VAL) {
            cx_best[v] = std::max(cx_best[v], u_[v] + shift);
            vmax = std::max(vmax, u_[v]);
            if (cx_first[v] < 0 && u_[v] + shift > cutoff) cx_first[v] = it;
          }
        }
        if (cxe.use && vmax > lbs) lbs = vmax;   // DNN_CERTX_USE: the best valid value counts
      }
      cx_cpu += thread_cpu() - tx0;
    }
    if (lbs > -HUGE_VAL) {
      double lb = nextafter(lbs + shift, -HUGE_VAL);   // covers the addition's rounding
      if (m_on && best > -HUGE_VAL &&
          (lb < best - 1e-6 * (1 + fabs(best)) ||
           (std::isfinite(cutoff) && cutoff - std::max(best, lb) < 1e-3 * std::max(1.0, fabs(cutoff))))) {
        // the certificate fell below the running best (the momentum overshot):
        // back to the last true iterate, cap halved; or near the cutoff: plain
        // iterations from the last true iterate for the rest of the node
        memcpy(Y, Yt.data(), (size_t)n1 * n1 * sizeof(double));
        memcpy(Z, Zt.data(), (size_t)n1 * n1 * sizeof(double));
        m_a = 1; m_c = -1;
        if (lb < best - 1e-6 * (1 + fabs(best))) { m_cap *= 0.5; m_drops++; }
        else m_on = false;
      }
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
    if (par.hopeless_after > 0 && it >= par.hopeless_after && best > -HUGE_VAL &&
        it >= confirm_at) {
      // extrapolate the running best certified value: with d1, d2 its gains
      // over the last two windows of ~50 iterations and r = d1/d2, a
      // geometric tail would add d1 r/(1-r); allow 5 times that
      int w = std::max(1, (50 + check_every - 1) / check_every);
      int h = (int)besthist.size();
      bool fire = false;
      if (par.hopeless_gap_ratio > 0) {
        // gap-ratio rule: the gap to the cutoff over the last w checks
        if (h > w && besthist[h - 1 - w] > -HUGE_VAL) {
          double g1 = cutoff - besthist[h - 1], g0 = cutoff - besthist[h - 1 - w];
          if (g0 > 0 && g1 / g0 > par.hopeless_gap_ratio) fire = true;
        }
      } else if (h > 2 * w && besthist[h - 1 - 2 * w] > -HUGE_VAL) {
        double d1 = besthist[h - 1] - besthist[h - 1 - w];
        double d2 = besthist[h - 1 - w] - besthist[h - 1 - 2 * w];
        double r = d2 > 0 ? std::min(d1 / d2, 0.95) : (d1 > 0 ? 0.95 : 0.0);
        if (best + 5 * d1 * r / (1 - r) < cutoff) fire = true;
      }
      if (fire) {
        // near the cutoff a give-up must be confirmed: run hopeless_delay more
        // iterations and give up only if the rule fires again after them
        // (once per node; the rule is not consulted during the deferral)
        if (par.hopeless_zone > 0 && confirm_at < 0 &&
            cutoff - best <= par.hopeless_zone * std::max(1.0, fabs(cutoff))) {
          confirm_at = it + par.hopeless_delay;
          if (par.verbose) printf("it %5d give-up deferred to %d: %.4f within %.2g of %.4f\n",
                                  it, confirm_at, best, par.hopeless_zone, cutoff);
        } else {
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
  if (cx_on && !cxe.quiet) {
    printf("certx-summary call %d status %s iters %d cutoff %.4f first std %d exact %d hg %d ls %d avg %d hgr %d"
           " hgw %d proj %d best std %.4f exact %.4f hg %.4f ls %.4f avg %.4f hgr %.4f hgw %.4f proj %.4f"
           " cpu_certx %.2f\n",
           call_idx, dnn_status_name(status), res.iterations, cutoff, cx_first[0], cx_first[1], cx_first[2],
           cx_first[3], cx_first[4], cx_first[5], cx_first[6], cx_first[7], cx_best[0], cx_best[1], cx_best[2],
           cx_best[3], cx_best[4], cx_best[5], cx_best[6], cx_best[7], cx_cpu);
    fflush(stdout);
  }
  if (par.verbose >= 2)
    printf("profile (s): forward %.3f  eig %.3f  expand %.3f  Y/Z-step %.3f  certify %.3f\n",
           w.prof[0], w.prof[1], w.prof[2], w.prof[3], w.prof[4]);
  if (par.verbose >= 2 && par.eig_mode == DNN_EIG_SUBSPACE)
    printf("eig split (s): full ssyevd %.3f (%d refreshes, K=%d)  filter+RR %.3f (%d steps, mean kb %.0f, mean nneg %.0f)  form R %.3f\n",
           w.prof[5], w.n_refresh, w.K, w.prof[6], w.n_sub,
           w.n_sub ? (double)w.sum_kb / w.n_sub : 0.0, w.n_sub ? (double)w.sum_nneg / w.n_sub : 0.0,
           w.prof[7]);
  if (par.verbose >= 2 && par.eig_mode == DNN_EIG_SUBSPACE && w.n_sub > 0) {
    const double ns = w.n_sub, nr = std::max(1, w.n_sub + w.n_refresh);
    printf("step split (ms per subspace step): weyl %.1f  convert %.1f  pack %.1f  W*Y %.1f  cheb %.1f  cholqr %.1f  rr-gemm %.1f  rr-eig %.1f  x-upd %.1f\n",
           w.sp[0] / nr, w.sp[1] / nr, w.sp[2] / ns, w.sp[3] / ns, w.sp[4] / ns, w.sp[5] / ns, w.sp[6] / ns, w.sp[7] / ns, w.sp[8] / ns);
    printf("form R (ms per iteration): F %.1f  syrk %.1f  symmetrize %.1f\n", w.sp[9] / nr, w.sp[10] / nr, w.sp[11] / nr);
  }
  if (par.verbose >= 2 && par.eig_mode == DNN_EIG_SUBSPACE)
    printf("refresh reasons: scheduled %d  weyl bound %d  block too small %d  factorization failed %d  (block grown %d times)\n",
           w.ref_sched, w.ref_weyl, w.ref_small, w.ref_fail, w.n_grow);
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
