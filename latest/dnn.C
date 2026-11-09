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

 Provenance of the certificate (noted Nov 8 2026). The form of the bound,
    a dual multiplier Z with V'ZV <= 0 and the minimum of <L + Z, Y> over a
    polyhedral set containing the feasible Y (the box), is Oliveira,
    Wolkowicz and Xu, "ADMM for the SDP relaxation of the QAP", Math. Prog.
    Comp. 2018, reused by Graham, Hu, Im, Li and Wolkowicz (INFORMS J.
    Comput. 2022). Bounding lambda_max by a Cholesky test is Rump (2006), in
    the tradition of Jansson, Chaykin and Keil's rigorous SDP bounds (SIAM J.
    Numer. Anal. 2007). The Gilmore-Lawler bound (1962) as a dual-feasible
    point of the level-1 RLT relaxation is Adams and Johnson (1994), and
    Hahn and Grant's dual ascent. What we believe is new here: the
    polyhedral set above (the RLT-1 column-block structure) evaluated by a
    Gilmore-Lawler computation on the residual M = L + Z, with every LAP
    value dual-feasible so inexact inner solves stay valid; the multiplier
    Z' = Z + beta (Y - VRV') that is negative semidefinite on the face by
    construction, with no second eigendecomposition; and the explicit
    floating-point margins in place of interval arithmetic. A literature
    check on Nov 8 (three searches) found no earlier Gilmore-Lawler
    evaluation inside an SDP dual certificate; it was not an exhaustive
    review. The certificate study (experiments/wil50_feasibility/roadmap/
    certificate) shows the remaining freedom at a fixed iterate is a skew
    transfer; a Hahn-Grant polish of it bought about 1.04x in iterations to
    fathom (the DNN_CERTX experiment, removed Nov 9; the report has the numbers).

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
void F77(ssyevd)(const char *jobz, const char *uplo, const lapack_int *n, float *a,
                 const lapack_int *lda, float *w, float *work, const lapack_int *lwork,
                 lapack_int *iwork, const lapack_int *liwork, lapack_int *info,
                 size_t, size_t);
// the pieces of ssyevd, for a refresh that back-transforms only the columns the
// tracker reads (the partial refresh)
void F77(ssytrd)(const char *uplo, const lapack_int *n, float *a, const lapack_int *lda,
                 float *d, float *e, float *tau, float *work, const lapack_int *lwork,
                 lapack_int *info, size_t);
void F77(sstedc)(const char *compz, const lapack_int *n, float *d, float *e, float *z,
                 const lapack_int *ldz, float *work, const lapack_int *lwork, lapack_int *iwork,
                 const lapack_int *liwork, lapack_int *info, size_t);
void F77(sormtr)(const char *side, const char *uplo, const char *trans, const lapack_int *m,
                 const lapack_int *n, const float *a, const lapack_int *lda, const float *tau,
                 float *c, const lapack_int *ldc, float *work, const lapack_int *lwork,
                 lapack_int *info, size_t, size_t, size_t);
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
  cut_round = 0;
  warm_primal = 0;
  child_bounds = 0;
  ycuts = NULL;
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
  std::vector<float> trd_e, trd_tau;   // ssytrd's off-diagonal and reflector scalars (the partial refresh)
  std::vector<double> Wprev;
  int kb, since_refresh;
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
  // CP-ADMM (DNNParams::ycuts): the projected Y, the certificate's objective and
  // cut term; allocated only then
  std::vector<double> Ycp, Lcp, Lcp_add;
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
    trd_e.assign(K, 0); trd_tau.assign(K, 0);
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
// precision, so the eigensolver's conversion pass is saved. W is not symmetrised
// (it is symmetric to rounding by construction; the eigensolver reads one
// triangle or multiplies by the full matrix; unconditional since Nov 9).
// force_sym: symmetrise anyway (the certificate's transform, whose multiplier is
// symmetric only to rounding).
void forward(Work &w, const double *A, double *Wout, const double *B = NULL,
             double coef = 0, bool rows_done = false, float *Wf_out = NULL, bool force_sym = false)
{
  const bool nosym = !force_sym;   // symmetrise only for the certificate (force_sym); the
                                   // iterate is symmetric to rounding (speed/nov6_core, Nov 8)
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
  if (Wf_out && nosym) {   // with the symmetrisation below Wf would be stale, so it is not written then
    for (int s = 0; s < K; s++) {
      const double *row = T + (size_t)w.kept[s] * n1;
      double *wo = Wout + (size_t)s * K;
      float *wf = Wf_out + (size_t)s * K;
      for (int t = 0; t < K; t++) { const double v = row[w.kept[t]]; wo[t] = v; wf[t] = (float)v; }
    }
    w.wf_ready = true;
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
  // No symmetrisation pass (a read and write of the n1 x n1 matrix; measured in
  // experiments/wil50_feasibility/speed/nov6_core, unconditional since Nov 9).
  // The transform of a symmetric R is symmetric to rounding, so Y, Z and the
  // certificate multiplier are symmetric to rounding only; the certificate then reads both arrow
  // entries of M (cert_poly) and symmetrises its own transform (cert_lam), which
  // keeps it a bound for sym(Zc), and <Zc, Y> = <sym(Zc), Y> for every symmetric Y.
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
                       int degree, int guard, int tight_lam, bool cert_step = true)
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
  while (true) {
    if (full) {
      lapack_int lwork = (lapack_int)w.eworkf.size(), liwork = (lapack_int)w.eiwork.size(),
                 info = 0;
      double tf0 = thread_cpu();
      memcpy(X, Wf, kk * sizeof(float));
      // The refresh (experiments/wil50_feasibility/speed/nov6_core, Nov 8):
      // ssyevd's three steps called one by one (ssytrd, sstedc, sormtr: the same
      // eigenvalues and vectors bit for bit), with the back-transformation applied
      // to the nneg + guard lowest eigenvectors and the top one only, which is all
      // the tracker reads; the other K - kb columns of ssyevd's 2 K^3 flops are
      // never formed (about a third of the refresh at K = 2117).
      {
        float *Zt = w.Ts.data(), *e = w.trd_e.data(), *tau = w.trd_tau.data();   // Ts is free until form R
        F77(ssytrd)("L", &K, X, &K, theta.data(), e, tau, w.eworkf.data(), &lwork, &info, 1);
        if (info == 0) F77(sstedc)("I", &K, theta.data(), e, Zt, &K, w.eworkf.data(), &lwork,
                                   w.eiwork.data(), &liwork, &info, 1);
        if (info == 0) {
          nneg = 0;
          while (nneg < K && theta[nneg] < 0) nneg++;
          kb = std::min<lapack_int>(K, nneg + guard);
          F77(sormtr)("L", "L", "N", &K, &kb, X, &K, tau, Zt, &K, w.eworkf.data(), &lwork, &info, 1, 1, 1);
          if (info == 0 && tight_lam > 0 && kb < K) {
            const lapack_int one = 1;
            F77(sormtr)("L", "L", "N", &K, &one, X, &K, tau, Zt + (size_t)(K - 1) * K, &K,
                        w.eworkf.data(), &lwork, &info, 1, 1, 1);
          }
          if (info == 0) {
            if (tight_lam > 0) memcpy(w.vtop.data(), Zt + (size_t)(K - 1) * K, K * sizeof(float));
            memcpy(X, Zt, (size_t)kb * K * sizeof(float));
          }
        }
      }
      w.prof[5] += thread_cpu() - tf0; w.n_refresh++;
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
    double ts0 = thread_cpu();
    float *Y = w.Ys.data(), *T = w.Ts.data(), *Yp = X;
    const float c = (float)(0.5 * (w.lam_hi + w.theta_top)),
                h = (float)(0.5 * (w.lam_hi - w.theta_top));
    const size_t nb = (size_t)K * kb;
    WY(Yp, T);
    double tc0 = thread_cpu();
    for (size_t t = 0; t < nb; t++) Y[t] = (T[t] - c * Yp[t]) / h;
    w.sp[4] += (thread_cpu() - tc0) * 1e3;
    for (int d = 2; d <= degree; d++) {
      WY(Y, T);
      tc0 = thread_cpu();
      for (size_t t = 0; t < nb; t++) Yp[t] = 2 * (T[t] - c * Y[t]) / h - Yp[t];
      w.sp[4] += (thread_cpu() - tc0) * 1e3;
      std::swap(Y, Yp);
    }
    // Y now holds the filtered block (in Ys or Xs); orthonormalize by one CholQR
    // pass (a second pass was measured not worth it: node_speed.md trim d, Nov 8)
    float *G = w.Hs.data();
    bool ok = true;
    tc0 = thread_cpu();
    do {
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
    } while (0);
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
      kb += add;
      w.n_grow++;
    } else {
      kb = std::min<lapack_int>(kb, nneg + guard);   // keep the lowest Ritz vectors
    }
    break;
  }
  w.kb = kb;
  if (full) w.theta_top = theta[kb - 1];
  *npos = K - nneg;

  // R = W - sum_{theta<0} theta x x'. The rank-nneg term is formed in fp32
  // (ssyrk on the fp32 eigenvectors, which is their precision anyway) and added
  // to W in fp64, except on certificate iterations, which keep dsyrk
  // (node_speed.md trim b, speed/r32_noncert; unconditional since Nov 9).
  // DNN_R32=1 forces fp32 on certificate iterations too. The certificate multiplier
  // Zc = Z + beta (Y - VRV') reads its own iteration's R and no earlier one
  // (V'ZcV = beta (W - R)), so only certificate iterations need W - R =
  // -fl64(FF') negative semidefinite; elsewhere R is just the ADMM step.
  static const int r32 = getenv("DNN_R32") ? atoi(getenv("DNN_R32")) : 0;
  double tr0 = thread_cpu();
  lapack_int k = nneg;
  if ((r32 || !cert_step) && k > 0) {   // fp32 except on certificate iterations (speed/r32_noncert)
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
  double lam = HUGE_VAL;     // the rigorous lambda_max bound used (scaled units)
  double margin_poly = 0;    // the margin's rounding terms (polyhedral part and L), without sdp_penalty
};

/*
 The certificate's two parts, so that the polyhedral value can be evaluated
 alone (DNN_CORE_LAZYCERT: it is an upper bound on the certified value, since
 the semidefinite penalty is never negative). certify() below calls both in
 the order it always did, so its result is unchanged bit for bit.
*/
// lambda_max(V'ZmV) <= lam, by Cholesky of delta I - S (Rump 2006: if
// floating-point Cholesky of A - cI succeeds, A is positive definite).
// Sets ci.lam, ci.sdp_penalty and ci.chol_ok; false when no shift up to 1e-3
// relative passes. Uses w.T, w.W, w.W2.
static bool cert_lam(Work &w, const double *Zm, double scale, CertInfo &ci)
{
  const int m = w.m, n1 = w.n1;
  const lapack_int K = w.K;
  ci.chol_ok = false;
  ci.sdp_penalty = 0;
  double *S = w.W.data();
  forward(w, Zm, S, NULL, 0, false, NULL, true);
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
  if (lam == HUGE_VAL) return false;
  ci.chol_ok = true;
  // error of the computed transform against an exact orthonormal V
  lam += 64.0 * m * U_ROUND * znorm;
  ci.lam = lam;
  ci.sdp_penalty = (m + 1) * scale * lam;
  return true;
}

// The polyhedral part on M = L + scale * Zm (unscaled): the minimum over P
// (box or Gilmore-Lawler) and, in ci.margin_poly, the rounding terms of the
// margin that do not depend on the semidefinite part. Returns the value val
// (not yet reduced by any margin); -HUGE_VAL only on a non-finite value.
static double cert_poly(Work &w, const double *Zm, double scale, double sumabsL, int poly,
                        int inner_lap, double exact_above, double exact_below, CertInfo &ci)
{
  const int m = w.m, n1 = w.n1;
  ci.poly_reduce = ci.poly = -HUGE_VAL;
  ci.margin_poly = 0;
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
          // the two arrow entries M_(ij),0 and M_0,(ij) separately: equal (and the
          // sums below exactly those of 2 M_(ij),0) when M is symmetric, and a valid
          // row-by-row evaluation of <M, yy'> when it is symmetric to rounding only
          // (no symmetrisation pass after the transforms)
          const double d0 = lr[0] + scale * zr[0], d0t = L[r] + scale * Zm[r], dd = lr[r] + scale * zr[r];
          absM += fabs(d0) + fabs(d0t) + fabs(dd) + fabsl(inner);
          cm[i * m + j] = (double)((long double)d0 + d0t + dd + inner);
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
  // summation and rounding errors of the polyhedral part and objective rounding;
  // (m + 2) u sum|L|: rounding of L itself and of cancellations in a row
  ci.margin_poly = 8.0 * (m + 2) * U_ROUND * (double)absM + (m + 2) * U_ROUND * sumabsL;
  if (!std::isfinite((double)val) || !std::isfinite(ci.margin_poly)) return -HUGE_VAL;
  return (double)val;
}

double certify(Work &w, const double *Zm, double scale, double sumabsL, int poly,
               int inner_lap, double exact_above, double exact_below, CertInfo &ci)
{
  ci.chol_ok = false;
  ci.margin = ci.sdp_penalty = 0;
  ci.poly_reduce = ci.poly = -HUGE_VAL;
  if (!cert_lam(w, Zm, scale, ci)) return -HUGE_VAL;
  const double val = cert_poly(w, Zm, scale, sumabsL, poly, inner_lap, exact_above, exact_below, ci);
  // the semidefinite term is added last, as before the split: the same sum, bit for bit
  double margin = ci.margin_poly + ci.sdp_penalty;
  ci.margin = margin;
  // a NaN or infinity anywhere (OpenBLAS's Cholesky accepts NaN) certifies nothing
  if (!std::isfinite(ci.lam) || !std::isfinite(val) || !std::isfinite(margin))
    return -HUGE_VAL;
  return val - margin;
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
 (<Z, Y> is the same for every permutation of the child). Z is always
 carried over; ymode (QAP_DNN_WARM_Y) says whether Y is too: 0 no (the
 child's Y starts from scratch), 1 the parent's Y restricted to the child's
 indices, 2 (the default) the parent's Y conditioned on x_ij = 1. Early
 tests favoured Z alone; the conditioned Y won once the rest of the recipe
 was in place (experiments/wil50_feasibility/warm_y).
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
 DNN_BETA_SCHED (experiments/wil50_feasibility/speed/beta_sched): opt-in
 schedules of the ADMM penalty. Unset, beta is par.beta * m/3 for the whole
 call and nothing below runs. beta changes only at the top of an iteration
 whose R-step runs a scheduled full eigendecomposition (subspace mode; every
 iteration in the exact modes), and never after iteration `freeze`: the
 tracker restarts from a full ssyevd of the jumped W, and the give-up rule
 sees a stationary iteration. Z is the unscaled multiplier (Y-step
 Pi_P(VRV' - (L/scale + Z)/beta), R-step P_psd(V'(Y + Z/beta)V), Z-step
 Z += gamma beta (Y - VRV')), so nothing is rescaled when beta changes: the
 fixed points do not depend on beta, a child's mapped Z is beta-free, and the
 certificate multiplier Z' = Z + beta (Y - VRV') satisfies V'Z'V =
 beta (W - R) <= 0 for the beta of its own iteration (W, the steps and Z' all
 use it). Comma-separated, kind first:
   two,f=F,t=T[,f2=F2][,force=1]  beta = F * base until the first scheduled
       refresh at or after iteration T, then F2 * base (F2 = 1); with
       force=1 the change is at iteration T itself, which then refreshes
   rb,c=C[,every=40][,max=4][,lo=0.25][,hi=4]  residual balancing: at each
       check S += log(rhat/shat) - log C, rhat = |Y - VRV'| / max(|Y|, |VRV'|),
       shat = beta |dY| / max(|L/scale|, |Z|) (Frobenius norms of the lifted
       matrices); at a scheduled refresh at least `every` iterations after
       the last change, with fewer than `max` changes made,
       f = clip(exp(S / 2n), 0.5, 2), and beta *= f if f is outside
       [1/1.5, 1.5] (beta kept in [lo, hi] * base), S = n = 0
   freeze=N   no change after iteration N (200)
   nodes=all|cold|warm  the calls it applies to (warm: start given; all)
   w=FW       warm-started calls instead run at FW * base throughout (no
              change), cold ones as the kind says; `w=FW` alone leaves cold
              calls at base. Measured: two,f=1.6667,t=60,w=0.7 (cold
              0.5 -> 0.3 m/3, warm 0.21 m/3; experiments/wil50_feasibility/
              speed/beta_sched); on 24 warm children w=0.25 needs 40% fewer
              iterations than 0.7 (speed/nov6_admm)
   wc=FC      warm-started rounds of a cut phase (par.cut_round, set by
              dnncut.C's rounds) run at FC * base instead of w's factor, so
              that children and cut rounds can differ (unset: w for both). A
              pooled child carries a cut term too (par.L_add) but is a child:
              it takes w (until Nov 6 it took wc, which made the laptop's
              wlowc arm identical to fast; pipeline_ab/README.md)
 Not with the accelerations (Anderson keeps Z/beta; Halpern, momentum): off
 there. Read at every call (tests set and unset it within a process).
*/
struct BetaSched {
  int kind = 0;          // 0 off, 1 two-phase, 2 residual balancing, 3 fixed (set per call: w)
  double f = 1, f2 = 1, c = 1, lo = 0.25, hi = 4, w = 0, wc = 0;   // wc: w for warm calls with a cut term (L_add)
  int t = 60, force = 0, every = 40, maxchg = 4, freeze = 200, nodes = 0;   // nodes: 0 all, 1 cold, 2 warm
  BetaSched()
  {
    const char *e = getenv("DNN_BETA_SCHED");
    if (!e) e = "w=0.36";      // the record recipe's warm-child beta, default since Nov 8 2026; DNN_BETA_SCHED="" turns it off
    if (!*e) return;
    std::string s = e;
    for (size_t pos = 0; pos <= s.size();) {
      size_t end = s.find(',', pos);
      if (end == std::string::npos) end = s.size();
      const std::string tok = s.substr(pos, end - pos);
      pos = end + 1;
      if (tok == "two") kind = 1;
      else if (tok == "rb") kind = 2;
      else if (tok.size() > 2 && tok.find('=') != std::string::npos) {
        const std::string k = tok.substr(0, tok.find('=')), v = tok.substr(tok.find('=') + 1);
        if (k == "f") f = atof(v.c_str());
        else if (k == "f2") f2 = atof(v.c_str());
        else if (k == "t") t = atoi(v.c_str());
        else if (k == "force") force = atoi(v.c_str());
        else if (k == "c") c = atof(v.c_str());
        else if (k == "every") every = atoi(v.c_str());
        else if (k == "max") maxchg = atoi(v.c_str());
        else if (k == "lo") lo = atof(v.c_str());
        else if (k == "hi") hi = atof(v.c_str());
        else if (k == "freeze") freeze = atoi(v.c_str());
        else if (k == "nodes") nodes = v == "cold" ? 1 : v == "warm" ? 2 : 0;
        else if (k == "w") w = atof(v.c_str());
        else if (k == "wc") wc = atof(v.c_str());
      }
    }
    if (!(f > 0) || !(f2 > 0) || !(c > 0) || !(lo > 0) || !(hi >= lo) || !(w >= 0) || !(wc >= 0)) kind = 0, w = 0, wc = 0;
  }
};

DNNResult dnn_bound(const MAT *A, const MAT *B, const MAT *C, double shift,
                    double cutoff, const DNNParams &par, const DNNWarm *start,
                    DNNWarm *save)
{
  DNNResult res;
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
  // beta changes only under DNN_BETA_SCHED (bs below)
  const double beta0 = (par.beta > 0 ? par.beta : 1.0) * m / 3.0;
  double beta = beta0, ibeta = 1.0 / beta;
  const double gamma_par = par.gamma, gamma1_par = par.gamma1;
  const int check_every = std::max(1, par.check_every);
  // DNN_ADMM_CROSS=n (off unless set; experiments/wil50_feasibility/speed/nov6_admm):
  // between two regular certificates, when the running best is below the cutoff
  // and rising, predict the iteration at which it crosses from the last two
  // certified values (their gain per iteration, decayed by the ratio of the last
  // two gains) and certify there when that is before the next regular check; at
  // most n extra certificates per interval. An extra check tests the fathom only:
  // the give-up and stall rules see the regular checks alone, so with the switch
  // a node can only fathom earlier, never decide anything else differently.
  const int cross_max = getenv("DNN_ADMM_CROSS") ? atoi(getenv("DNN_ADMM_CROSS")) : 2;   // 2 since Nov 8 2026 (was 0)
  // DNN_ADMM_CROSS_FIRST=k (with DNN_ADMM_CROSS, warm-started calls only): one more
  // extra certificate at iteration k < check_every, so that the first regular check
  // already has two points and can predict a crossing in its own first interval
  // (in the production search a third of the fathoms come at the second check)
  const int cross_first = cross_max > 0 && getenv("DNN_ADMM_CROSS_FIRST") ? atoi(getenv("DNN_ADMM_CROSS_FIRST")) : 0;
  int cross_at = -1, cross_n = 0, cross_extra = 0;
  std::vector<std::pair<int, double>> cr_hist;   // (iteration, running best) at every check

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

  // DNN_FUSE=1: the Y/Z step writes the next forward's row phase into T while the row
  // is in cache, and forward writes Wf with W (node_speed.md trim f); not with the
  // accelerations, which rewrite Y and Z after the step
  static const int fuse_env = getenv("DNN_FUSE") ? atoi(getenv("DNN_FUSE")) : 0;
  const bool fuse_ok = fuse_env > 0 && par.anderson == 0 && par.halpern == 0 && !m_on && !par.ycuts;
  const int guard_div = 20;   // subspace block guard = K / 20 (node_speed.md trim c: 40 was no better)
  // DNN_BETA_SCHED (see BetaSched): off unless set
  BetaSched bs;
  const bool warm_call = start && start->Z.size() == (size_t)n1 * n1;
  const bool accel = par.anderson > 0 || par.halpern > 0 || par.momentum > 0;
  // CP-ADMM (par.ycuts): cuts projected in the Y-step, certified with their
  // multipliers at every check (see DNNYCuts)
  DNNYCuts *const yc = par.ycuts;
  if (yc && (accel || par.prototype_cert || par.child_bounds)) {
    res.status = DNN_ERROR; res.cpu = thread_cpu() - t0; return res;
  }
  double cp_cpu = 0;   // CPU seconds in the hook's projection
  if (yc) { w.Ycp.assign((size_t)n1 * n1, 0.0); w.Lcp.assign((size_t)n1 * n1, 0.0); w.Lcp_add.assign((size_t)n1 * n1, 0.0); }
  {
    if (accel || (bs.nodes == 1 && warm_call) || (bs.nodes == 2 && !warm_call))
      bs.kind = 0, bs.w = 0;
    if (warm_call && bs.w > 0) bs.kind = 3;   // w: a fixed beta on warm-started calls
    // wc: the factor for the warm rounds of a cut phase (par.cut_round, dnncut.C),
    // so that children and cut rounds can run at different betas
    // (experiments/wil50_feasibility/speed/nov6_admm); unset: w for both. Not
    // par.L_add: a pooled child carries a cut term as well, and is a child
    if (warm_call && par.cut_round && bs.wc > 0) { bs.kind = 3; bs.w = bs.wc; }
  }
  // DNN_ADMM_WARM_BETA=F1:T:F2[:force] (off unless set): a warm-started call runs at
  // F1 x base until the first scheduled refresh at or after iteration T, then at
  // F2 x base (force=1: at T itself, with a full refresh there), whatever
  // DNN_BETA_SCHED says for warm calls; cold calls are untouched. Same mechanism as
  // BetaSched kind 1 (a change only where the R-step runs a full eigendecomposition;
  // Z is the unscaled multiplier, so nothing is rescaled). Not on the rounds of a
  // cut phase (par.cut_round), which keep DNN_BETA_SCHED's w or wc: the schedule is
  // for a node's own solve, low early so that warm children fathom fast, higher
  // after T so that a miss hands its cut phase a state like today's (Nov 6).
  if (warm_call && !accel && !par.cut_round) {
    const char *e = getenv("DNN_ADMM_WARM_BETA");
    double f1 = 0, f2 = 0; int t = 0, force = 0;
    if (e && *e && sscanf(e, "%lf:%d:%lf:%d", &f1, &t, &f2, &force) >= 3 && f1 > 0 && f2 > 0 && t > 0) {
      bs.kind = 1; bs.f = f1; bs.t = t; bs.f2 = f2; bs.force = force; bs.w = 0;
      bs.freeze = std::max(bs.freeze, t + 2 * std::max(1, par.refresh));
    }
  }
  // DNN_ADMM_WARM_DAMP=D:K (off unless set): on a warm-started call both dual steps
  // are scaled by D for the first K iterations (gamma1, gamma -> D gamma1, D gamma):
  // the parent's multiplier is near the child's optimum, and the start's primal
  // residual would otherwise move it far in the first iterations. Certificates
  // are unaffected (the R-step multiplier does not depend on the dual steps).
  double wd_fac = 1; int wd_k = 0;
  if (warm_call && !accel) {
    const char *e = getenv("DNN_ADMM_WARM_DAMP");
    double d = 0; int k = 0;
    if (e && *e && sscanf(e, "%lf:%d", &d, &k) == 2 && d > 0 && k > 0) { wd_fac = d; wd_k = k; }
  }
  if (par.verbose && (cross_max > 0 || wd_k > 0))
    printf("admm cross %d warm_damp %g:%d\n", cross_max, wd_fac, wd_k);
  int bs_phase = 0, bs_last = 0, bs_changes = 0, bs_n = 0;
  double bs_S = 0, bs_Lnorm = 0;
  // DNN_CORE_LAZYCERT=1 (experiments/wil50_feasibility/speed/nov6_core): Rump's test
  // only at checks that can decide something (see the check below). Off by default;
  // not with the prototype certificate, the certificate variants, the stall rule or
  // the child bounds, which read every check's certified value. Under CP-ADMM
  // (par.ycuts) only with DNN_CP_LAZYCERT=1 as well: the skipped checks keep their
  // multipliers in the hook (DNNYCuts::snapshot), and the exit certifies them with
  // their own cut term (roadmap/admm/cpadmm, "In the search")
  static const int lazy_env = getenv("DNN_CORE_LAZYCERT") ? atoi(getenv("DNN_CORE_LAZYCERT")) : 1;   // on since Nov 8 2026
  static const int cp_lazy_env = getenv("DNN_CP_LAZYCERT") ? atoi(getenv("DNN_CP_LAZYCERT")) : 0;
  const bool lazy_ok = lazy_env > 0 && !par.prototype_cert && par.stall_window == 0 &&
                       !par.child_bounds && (!par.ycuts || cp_lazy_env > 0);
  std::vector<double> lz_Z;                 // the stored multiplier with the best polyhedral estimate
  bool lz_have = false, lz_exit = false;    // stored; the run stopped on a check that skipped the test
  double lz_pen = 0;                        // the semidefinite penalty of the last full certificate
  double lz_est = -HUGE_VAL;                // this check's estimate (unshifted)
  double lz_best_proxy = -HUGE_VAL, lz_best_est = -HUGE_VAL;   // running maxima (shifted)
  int lz_skipped = 0;
  if (bs.kind == 1) { beta = bs.f * beta0; ibeta = 1.0 / beta; bs_phase = 1; }
  if (bs.kind == 3) { beta = bs.w * beta0; ibeta = 1.0 / beta; }
  if (bs.kind) bs_Lnorm = sqrt(dot4(L, L, (size_t)n1 * n1)) * iscale;
  if (bs.kind && par.verbose)
    printf("bsched kind %d beta %.6g (base %.6g) f %g t %d f2 %g force %d c %g every %d max %d freeze %d nodes %d w %g\n",
           bs.kind, beta, beta0, bs.f, bs.t, bs.f2, bs.force, bs.c, bs.every, bs.maxchg, bs.freeze,
           bs.nodes, bs.w);
  for (int it = 1; it <= par.max_iter; it++) {
    if (bs.kind && it > 1 && it <= bs.freeze) {
      // a change of beta, only where the R-step will run a scheduled full
      // eigendecomposition (force: make this iteration one)
      const bool sub = par.eig_mode == DNN_EIG_SUBSPACE;
      bool at_refresh = !sub || !w.have_basis || w.since_refresh >= par.refresh;
      double nb = beta;
      if (bs.kind == 1 && bs_phase == 1 && it >= bs.t && (at_refresh || bs.force)) {
        nb = bs.f2 * beta0; bs_phase = 2;
        if (sub && !at_refresh) { w.since_refresh = par.refresh; at_refresh = true; }
      } else if (bs.kind == 2 && at_refresh && bs_n > 0 && it - bs_last >= bs.every &&
                 bs_changes < bs.maxchg) {
        const double fr = std::min(2.0, std::max(0.5, exp(bs_S / (2.0 * bs_n))));
        if (fr > 1.5 || fr < 1 / 1.5) {
          nb = std::min(bs.hi * beta0, std::max(bs.lo * beta0, fr * beta));
          bs_S = 0; bs_n = 0; bs_last = it; bs_changes++;
        }
      }
      if (nb != beta) {
        if (par.verbose) printf("bsched it %d beta %.6g -> %.6g (factor %.4g of base)\n", it, beta, nb, nb / beta0);
        beta = nb; ibeta = 1.0 / beta;
        w.t_rows_ready = false;   // DNN_FUSE's rows were formed with the old 1/beta
      }
    }
    // R-step: W = V'(Y + Z/beta)V, R = P_psd(W), VRV = V R V'
    double tp = thread_cpu(), tq;
    if (par.anderson > 0) aa.begin(Y, Z, ibeta);
    forward(w, Y, w.W.data(), Z, ibeta, w.t_rows_ready,
            par.eig_mode == DNN_EIG_SUBSPACE ? w.Wf.data() : NULL);   // Wf with W (forward writes it unless force_sym)
    tq = thread_cpu(); w.prof[0] += tq - tp; tp = tq;
    const double gamma1 = it <= wd_k ? gamma1_par * wd_fac : gamma1_par;
    const double gamma = it <= wd_k ? gamma_par * wd_fac : gamma_par;
    const bool regular = it % check_every == 0 || it == par.max_iter;
    // cross_at: DNN_ADMM_CROSS's extra certificate; cross_first: its early one on warm calls
    const bool first_extra = warm_call && cross_first > 0 && cross_first < check_every && it == cross_first;
    bool check = regular || it == cross_at || first_extra;
    if (it == cross_at || first_extra) { cross_at = -1; cross_extra++; }
    bool rstep_cert = check && !par.prototype_cert;
    bool okp = par.eig_mode == DNN_EIG_SUBSPACE
                   ? proj_psd_subspace(w, w.W.data(), w.R.data(), &npos, par.refresh,
                                       par.cheb_degree, std::max(8, w.K / guard_div), par.tight_lam,
                                       rstep_cert)
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
    if (yc) {
      // CP-ADMM: Zc and the first dual step as below, X = VRV - (L/scale + Z)/beta
      // (in T, free without DNN_FUSE), its projection onto P and the cuts by the
      // hook, then the Z-step
      double *X = w.T.data(), *Yn = w.Ycp.data();
      for (int r = 0; r < n1; r++) {
        const double *lr = L + (size_t)r * n1, *vr = VRV + (size_t)r * n1;
        double *yr = Y + (size_t)r * n1, *zr = Z + (size_t)r * n1, *cr = Zc + (size_t)r * n1;
        double *xr = X + (size_t)r * n1;
        if (rstep_cert)
          for (int c = 0; c < n1; c++) cr[c] = zr[c] + beta * (yr[c] - vr[c]);
        if (gamma1 > 0)
          for (int c = 0; c < n1; c++) zr[c] += gamma1 * beta * (yr[c] - vr[c]);
        for (int c = 0; c < n1; c++) xr[c] = vr[c] - (lr[c] * iscale + zr[c]) * ibeta;
      }
      const double tc = thread_cpu();
      yc->project(X, Yn, beta, w.elim, it);
      cp_cpu += thread_cpu() - tc;
      double *rs = w.tmp.data() + n1, *ds = rs + n1;
      for (int r = 0; r < n1; r++) {
        const double *vr = VRV + (size_t)r * n1, *ynr = Yn + (size_t)r * n1;
        double *yr = Y + (size_t)r * n1, *zr = Z + (size_t)r * n1;
        for (int c = 0; c < n1; c++) {
          double y = ynr[c], res_ = y - vr[c];
          zr[c] += gamma * beta * res_;
          rs[c] = res_;
          ds[c] = y - yr[c];
          yr[c] = y;
        }
        pres2 += dot4(rs, rs, n1);
        dy2 += dot4(ds, ds, n1);
      }
    } else   // the plain Y/Z step
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
    if (bs.kind && check) {
      // residual balance at this check (VRV is intact until the certificate's refinement)
      const size_t nn = (size_t)n1 * n1;
      const double ny = sqrt(dot4(Y, Y, nn)), nv = sqrt(dot4(VRV, VRV, nn)), nz = sqrt(dot4(Z, Z, nn));
      const double pres = sqrt((double)pres2), dy = sqrt((double)dy2);
      const double rhat = pres / std::max(1e-300, std::max(ny, nv));
      const double shat = beta * dy / std::max(1e-300, std::max(bs_Lnorm > 0 ? bs_Lnorm : 1.0, nz));
      if (bs.kind == 2 && rhat > 0 && shat > 0 && it <= bs.freeze) {
        bs_S += log(rhat / shat) - log(bs.c);
        bs_n++;
      }
      if (par.verbose)
        printf("bsched it %d beta %.6g pres/(beta dY) %.4g rhat %.4g shat %.4g rhat/shat %.4g |Y| %.4g |VRV| %.4g |Z| %.4g\n",
               it, beta, dy > 0 ? pres / (beta * dy) : 0.0, rhat, shat, shat > 0 ? rhat / shat : 0.0, ny, nv, nz);
    }
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
    // CP-ADMM: this check's cut term G'mu goes into the objective the certificate
    // reads (w.L, swapped back below) and h'mu into the shift; else shift_c == shift
    double cp_hsum = 0, sumabsL_c = (double)sumabsL;
    bool cp_swapped = false;
    if (yc) {
      const size_t nn = (size_t)n1 * n1;
      std::fill(w.Lcp_add.begin(), w.Lcp_add.end(), 0.0);
      double abs_add = 0;
      if (yc->cert_term(scale, beta, w.Lcp_add.data(), abs_add, cp_hsum)) {
        for (size_t t = 0; t < nn; t++) w.Lcp[t] = L[t] + w.Lcp_add[t];
        sumabsL_c += abs_add + maxabs_of(w.Lcp_add.data(), nn);
        std::swap(w.L, w.Lcp);   // L (the pointer) still reads the plain objective
        cp_swapped = true;
      } else cp_hsum = 0;
    }
    const double shift_c = yc ? shift - cp_hsum : shift;
    const double exact_above_c = !yc ? exact_above : (std::isfinite(cutoff) ? cutoff - gap - shift_c : HUGE_VAL);
    if (par.prototype_cert) {
      // dnn.py: Zt = Z - V P_psd(V'ZV) V' (extra eigendecomposition)
      forward(w, Z, w.W.data());
      int np2;
      if (!proj_psd(w, w.W.data(), w.R.data(), &np2)) { status = DNN_ERROR; break; }
      expand(w, w.R.data(), Zc);
      for (size_t t = 0; t < (size_t)n1 * n1; t++) Zc[t] = Z[t] - Zc[t];
    }
    const double exact_below = cutoff - shift_c + 1e-6 * (1 + fabs(cutoff));
    w.want_red = par.child_bounds && par.poly != DNN_POLY_BOX;
    bool lazy_this = false;   // DNN_CORE_LAZYCERT: this check skipped Rump's test
    if (lazy_ok) {
      // The polyhedral value first: it is an upper bound on the certified value
      // (the semidefinite penalty is never negative), so when even it cannot clear
      // the cutoff the node cannot fathom at this check and Rump's test (forward +
      // dpotrf, about 130 of the certificate's 143 ms at K = 2117) is skipped. The
      // test runs at the first check, whenever the polyhedral value could fathom,
      // at the last iteration, and (below) at the exit on the check that stopped the
      // run and on the stored multiplier with the best polyhedral value, so the
      // fathom decisions and the returned lb are those of certifying every check.
      // Under CP-ADMM the values are those of the cut-modified objective (w.L swapped
      // above): shift_c, sumabsL_c and exact_above_c, which equal shift, sumabsL and
      // exact_above without it.
      ci.chol_ok = false; ci.margin = ci.sdp_penalty = 0;
      const double val = cert_poly(w, Zc, scale, sumabsL_c, par.poly, par.inner_lap,
                                   exact_above_c, exact_below, ci);
      const double proxy = val > -HUGE_VAL ? val - ci.margin_poly : -HUGE_VAL;
      const bool could_fathom = proxy > -HUGE_VAL && nextafter(proxy + shift_c, -HUGE_VAL) > cutoff;
      if (res.certificates == 0 || it == par.max_iter || could_fathom || proxy == -HUGE_VAL) {
        // decisive: the full certificate, as without the switch
        if (!cert_lam(w, Zc, scale, ci)) lbs = -HUGE_VAL;
        else {
          ci.margin = ci.margin_poly + ci.sdp_penalty;
          lbs = (std::isfinite(ci.lam) && std::isfinite(val) && std::isfinite(ci.margin)) ? val - ci.margin : -HUGE_VAL;
        }
        lz_pen = ci.chol_ok ? ci.sdp_penalty : lz_pen;
      } else {
        lazy_this = true;
        lz_skipped++;
        lbs = -HUGE_VAL;   // no certified value at this check
        // the estimate for the stopping rules, and the stored multiplier with the best
        // polyhedral value so far (verified at the exit if it can raise the result)
        lz_est = proxy - lz_pen - (cp_swapped ? cp_hsum : 0.0);   // of the plain objective (unshifted)
        if (yc) yc->snapshot(0, scale, beta);   // this check's multipliers, for an exit on it
        if (lz_est + shift > lz_best_proxy) {
          lz_best_proxy = lz_est + shift;
          if (lz_Z.empty()) lz_Z.resize((size_t)n1 * n1);
          memcpy(lz_Z.data(), Zc, (size_t)n1 * n1 * sizeof(double));
          if (yc) yc->snapshot(1, scale, beta);
          lz_have = true;
        }
      }
    } else {
      lbs = certify(w, Zc, scale, sumabsL_c, par.poly, par.inner_lap, exact_above_c,
                    exact_below, ci);
    }
    if (w.want_red) red_cur = w.red;
    if (!lazy_this) res.certificates++;
    if (!lazy_this && !ci.chol_ok) res.cert_failures++;
    if (!lazy_this && !par.prototype_cert && par.eig_mode != DNN_EIG_DOUBLE &&
        (!ci.chol_ok || (lbs + shift_c <= cutoff && lbs + shift_c + ci.sdp_penalty > cutoff))) {
      // single-precision R-steps leave V'Z'V slightly indefinite; when that
      // costs the fathom, project it in double: Z' <- Z' - V P_psd(V'Z'V) V'
      forward(w, Zc, w.W.data());
      int np2;
      if (!proj_psd(w, w.W.data(), w.R.data(), &np2)) { status = DNN_ERROR; break; }
      expand(w, w.R.data(), VRV);   // VRV is free until the next R-step
      for (size_t t = 0; t < (size_t)n1 * n1; t++) Zc[t] -= VRV[t];
      CertInfo c2;
      double l2 = certify(w, Zc, scale, sumabsL_c, par.poly, par.inner_lap,
                          exact_above_c, exact_below, c2);
      res.refinements++;
      if (l2 > lbs) { lbs = l2; ci = c2; if (w.want_red) red_cur = w.red; }
    }
    if (cp_swapped) std::swap(w.L, w.Lcp);
    // CP-ADMM: the bound of the cut-modified objective less h'mu (rounded up by
    // the hook), rounded down; lb = lbs + shift below rounds down once more
    if (cp_swapped && lbs > -HUGE_VAL) lbs = nextafter(lbs - cp_hsum, -HUGE_VAL);
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
    // the running best the stopping rules see: the certified best, or under
    // DNN_CORE_LAZYCERT the polyhedral estimate of the checks that skipped Rump's test
    double best_rule = best;
    if (lazy_ok) {
      lz_best_est = std::max(lz_best_est, lazy_this ? (lz_est > -HUGE_VAL ? nextafter(lz_est + shift, -HUGE_VAL) : -HUGE_VAL) : best);
      best_rule = lz_best_est;
    }
    tq = thread_cpu(); w.prof[4] += tq - tp; tp = tq;
    obj = dot4(L, Y, (size_t)n1 * n1) + shift;
    double pres = sqrt((double)pres2), dy = sqrt((double)dy2);
    if (par.verbose)
      printf("it %5d lb %.4f (cert %.4f, margin %.2e%s) primal %.4f pres %.2e dY %.2e npos %d cpu %.2f%s\n",
             it, best_rule, lazy_this ? lz_est + shift : (lbs > -HUGE_VAL ? lbs + shift : -HUGE_VAL), ci.margin,
             lazy_this ? ", lazy" : (ci.chol_ok ? "" : ", chol FAILED"), obj, pres, dy, npos, thread_cpu() - t0,
             regular ? "" : " extra");
    if (par.verbose && yc) {
      const std::string t = yc->trace();
      const double cv = lbs > -HUGE_VAL ? nextafter(lbs + shift, -HUGE_VAL) : -HUGE_VAL;
      printf("cp it %5d cert %.17g hsum %.6g cpcpu %.2f %s\n", it, cv, cp_hsum, cp_cpu, t.c_str());
    }
    // DNN_ADMM_CROSS: the crossing prediction from the (iteration, running best)
    // history of the checks; an extra check at a predicted crossing that did not
    // fathom predicts again from its own value (at most cross_max per interval)
    // It reads best_rule, the running best the stopping rules see: equal to best
    // unless DNN_CORE_LAZYCERT, whose skipped checks leave best unchanged.
    auto cross_predict = [&]() {
      const double bb = best_rule;
      if (cross_max <= 0 || !std::isfinite(cutoff) || bb <= -HUGE_VAL) return;
      cr_hist.push_back(std::make_pair(it, bb));
      if (regular) cross_n = 0;
      if (bb > cutoff || cross_n >= cross_max) return;
      const size_t n = cr_hist.size();
      if (n < 2) return;
      const int i1 = cr_hist[n - 2].first;
      const double b1 = cr_hist[n - 2].second;
      if (!(b1 > -HUGE_VAL) || it <= i1) return;
      double g = (bb - b1) / (it - i1);   // gain per iteration over the last interval
      if (!(g > 0)) return;
      // decayed as the gains decay: by the ratio of the last two gains (clamped to
      // [0.3, 1]) or, with two points only, by 0.5 (on wil50 warm children the gain
      // of the second interval is about 0.36 of the first's; a rate taken at face
      // value predicts the crossing 20-30 iterations early, diag_k1.py)
      double decay = 0.5;
      if (n >= 3 && cr_hist[n - 3].second > -HUGE_VAL && i1 > cr_hist[n - 3].first) {
        const double g0 = (b1 - cr_hist[n - 3].second) / (i1 - cr_hist[n - 3].first);
        if (g0 > 0) decay = std::min(1.0, std::max(0.3, g / g0));
      }
      g *= decay;
      const int next_regular = std::min(par.max_iter, (it / check_every + 1) * check_every);
      const double need = (cutoff - bb) / g;
      if (!(need < check_every)) return;
      const int p = it + std::max(1, (int)ceil(need));
      if (p < next_regular) { cross_at = p; cross_n++; }
    };
    if (!regular) {
      // an extra certificate: the fathom test only, then the next prediction
      if (best > cutoff) { status = DNN_FATHOMED; break; }
      cross_predict();
      continue;
    }
    lbhist.push_back(lazy_this ? lz_est + shift : (lbs > -HUGE_VAL ? lbs + shift : -HUGE_VAL));

    if (best > cutoff) { status = DNN_FATHOMED; break; }
    if (yc && yc->stop(it, best_rule)) { status = DNN_STALLED; break; }   // CP-ADMM: the caller's rule
    if (pres < par.tol && dy < par.tol) { status = DNN_CONVERGED; lz_exit = lazy_this; break; }
    if (par.time_limit > 0 && thread_cpu() - t0 > par.time_limit) { status = DNN_TIME_LIMIT; lz_exit = lazy_this; break; }
    besthist.push_back(best_rule);
    if (par.hopeless_after > 0 && it >= par.hopeless_after && best_rule > -HUGE_VAL &&
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
        if (best_rule + 5 * d1 * r / (1 - r) < cutoff) fire = true;   // best_rule == best unless DNN_CORE_LAZYCERT
      }
      if (fire) {
        // near the cutoff a give-up must be confirmed: run hopeless_delay more
        // iterations and give up only if the rule fires again after them
        // (once per node; the rule is not consulted during the deferral)
        if (par.hopeless_zone > 0 && confirm_at < 0 &&
            cutoff - best_rule <= par.hopeless_zone * std::max(1.0, fabs(cutoff))) {
          confirm_at = it + par.hopeless_delay;
          if (par.verbose) printf("it %5d give-up deferred to %d: %.4f within %.2g of %.4f\n",
                                  it, confirm_at, best_rule, par.hopeless_zone, cutoff);
        } else {
          status = DNN_HOPELESS;
          lz_exit = lazy_this;
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
    cross_predict();
  }
  if (status < 0) status = DNN_ITER_LIMIT;
  if (lazy_ok && status != DNN_FATHOMED && status != DNN_ERROR && (lz_exit || lz_have)) {
    // DNN_CORE_LAZYCERT: the checks that skipped Rump's test left no certified value.
    // Certify the multiplier of the check that stopped the run (still in Zc) and the
    // stored one with the best polyhedral estimate, if either can raise the result.
    for (int pass = 0; pass < 2; pass++) {
      const double *Zx = pass == 0 ? Zc : lz_Z.data();
      if (pass == 0 && !lz_exit) continue;
      if (pass == 1 && (!lz_have || lz_best_proxy <= best)) continue;
      // CP-ADMM: the cut term of that check's multipliers (DNNYCuts::snapshot), as at
      // a regular check; without it the plain objective's values
      double sh = shift, sab = (double)sumabsL, ea = exact_above, hs = 0;
      bool sw = false;
      if (yc) {
        const size_t nn = (size_t)n1 * n1;
        std::fill(w.Lcp_add.begin(), w.Lcp_add.end(), 0.0);
        double abs_add = 0;
        const int got = yc->cert_term_at(pass, w.Lcp_add.data(), abs_add, hs);
        if (got < 0) continue;   // not available: that check stays uncertified
        if (got > 0) {
          for (size_t t = 0; t < nn; t++) w.Lcp[t] = L[t] + w.Lcp_add[t];
          sab += abs_add + maxabs_of(w.Lcp_add.data(), nn);
          std::swap(w.L, w.Lcp);
          sw = true;
          sh = shift - hs;
          ea = std::isfinite(cutoff) ? cutoff - gap - sh : HUGE_VAL;
        } else hs = 0;
      }
      const double exact_below = cutoff - sh + 1e-6 * (1 + fabs(cutoff));
      CertInfo ci;
      w.want_red = par.child_bounds && par.poly != DNN_POLY_BOX;
      double lbs = certify(w, Zx, scale, sab, par.poly, par.inner_lap, ea, exact_below, ci);
      if (sw) std::swap(w.L, w.Lcp);
      if (sw && lbs > -HUGE_VAL) lbs = nextafter(lbs - hs, -HUGE_VAL);
      res.certificates++;
      if (!ci.chol_ok) res.cert_failures++;
      if (lbs > -HUGE_VAL) {
        const double lb = nextafter(lbs + shift, -HUGE_VAL);
        if (lb > best) { best = lb; if (w.want_red) { red_best = w.red; lbs_best = lbs; } }
        res.margin = std::max(res.margin, ci.margin);
      }
      if (par.verbose)
        printf("lazycert exit %s: cert %.4f (margin %.2e%s), best %.4f, %d checks skipped Rump's test\n",
               pass == 0 ? "last check" : "stored best", lbs > -HUGE_VAL ? lbs + shift : -HUGE_VAL,
               ci.margin, ci.chol_ok ? "" : ", chol FAILED", best, lz_skipped);
    }
  }
  if (par.verbose && cross_max > 0)
    printf("admm cross: %d extra certificates\n", cross_extra);
  if (par.verbose >= 2)
    printf("profile (s): forward %.3f  eig %.3f  expand %.3f  Y/Z-step %.3f  certify %.3f\n",
           w.prof[0], w.prof[1], w.prof[2], w.prof[3], w.prof[4]);
  if (par.verbose >= 2 && yc) printf("profile (s): cp projection %.3f (inside Y/Z-step)\n", cp_cpu);
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
