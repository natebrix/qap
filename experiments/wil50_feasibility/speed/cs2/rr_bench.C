// rr_bench.C: microbenchmark of the subspace step of latest/dnn.C
// (proj_psd_subspace) on steps dumped by dnntest (DNN_DUMP_RR=prefix
// DNN_DUMP_RR_AT=s1,s2,...; format in dnn.C, dump_subspace_step): the
// dumped block X (the previous step's Ritz vectors), its Ritz values, W and
// the filter interval are exactly what the step saw.
//
// Rayleigh-Ritz variants, on the same orthonormal block Y (filter of the
// dumped degree, one CholQR pass), repeated -r times, median thread CPU:
//   today: G = Y'T, symmetrize, ssyevd(kb), X = Y G, F = X(:,neg) |theta|^1/2 (fp64)
//   CS2:   G = Y'T, symmetrize, ssyevr(kb, theta > 0 only), S = U+ diag(theta+) U+' - G,
//          spstrf(S), F = Y(:,piv) L (strmm + sgemm, trapezoidal), F to fp64
// Both give R = W + F F' (the dsyrk that follows is common; its rank is nneg
// or r). Accuracy of each against the exact R = P_psd(W) (dsyevd of W in
// fp64), with -x: ||R - R_exact||_F / ||P_-(W)||_F. Also timed: ssytrd alone,
// full ssyevr, ssyevx and ssytrd + sstedc + sormtr for the positive pairs,
// and the "shift" variant (Cholesky of theta_top I - G, no eigensolve).
//
// Per-column filter degree (projection.md idea 3), with -x: the same error
// for the degree rules below, and the filter's cost in product columns.
//
// Build: build_bench.sh. Run pinned: taskset -c 3 ./rr_bench [-r reps] [-x] dump.bin ...
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

#ifndef BLAS_PREFIX
#define BLAS_PREFIX scipy_
#endif
#define CAT2(a, b) a##b
#define CAT(a, b) CAT2(a, b)
#define F77(name) CAT(BLAS_PREFIX, name##_)
typedef int li;
extern "C" {
void F77(ssyevd)(const char *, const char *, const li *, float *, const li *, float *, float *,
                 const li *, li *, const li *, li *, size_t, size_t);
void F77(dsyevd)(const char *, const char *, const li *, double *, const li *, double *, double *,
                 const li *, li *, const li *, li *, size_t, size_t);
void F77(ssyevr)(const char *, const char *, const char *, const li *, float *, const li *,
                 const float *, const float *, const li *, const li *, const float *, li *,
                 float *, float *, const li *, li *, float *, const li *, li *, const li *, li *,
                 size_t, size_t, size_t);
void F77(sgemm)(const char *, const char *, const li *, const li *, const li *, const float *,
                const float *, const li *, const float *, const li *, const float *, float *,
                const li *, size_t, size_t);
void F77(ssytrd)(const char *, const li *, float *, const li *, float *, float *, float *, float *,
                 const li *, li *, size_t);
void F77(sstedc)(const char *, const li *, float *, float *, float *, const li *, float *, const li *,
                 li *, const li *, li *, size_t);
void F77(sormtr)(const char *, const char *, const char *, const li *, const li *, const float *,
                 const li *, const float *, float *, const li *, float *, const li *, li *, size_t,
                 size_t, size_t);
void F77(ssyevx)(const char *, const char *, const char *, const li *, float *, const li *,
                 const float *, const float *, const li *, const li *, const float *, li *, float *,
                 float *, const li *, float *, const li *, li *, li *, li *, size_t, size_t, size_t);
void F77(dgemm)(const char *, const char *, const li *, const li *, const li *, const double *,
                const double *, const li *, const double *, const li *, const double *, double *,
                const li *, size_t, size_t);
void F77(ssyrk)(const char *, const char *, const li *, const li *, const float *, const float *,
                const li *, const float *, float *, const li *, size_t, size_t);
void F77(spotrf)(const char *, const li *, float *, const li *, li *, size_t);
void F77(spstrf)(const char *, const li *, float *, const li *, li *, li *, const float *, float *,
                 li *, size_t);
void F77(strsm)(const char *, const char *, const char *, const char *, const li *, const li *,
                const float *, const float *, const li *, float *, const li *, size_t, size_t,
                size_t, size_t);
void F77(strmm)(const char *, const char *, const char *, const char *, const li *, const li *,
                const float *, const float *, const li *, float *, const li *, size_t, size_t,
                size_t, size_t);
}

static double cpu()
{
  timespec ts;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

struct Dump {
  int K, kb, degree, guard;
  double lam_hi, theta_top;
  std::vector<float> theta, X, W;
};

static bool load(const char *fn, Dump &d)
{
  FILE *f = fopen(fn, "rb");
  if (!f) return false;
  int hdr[4];
  double lh[2];
  bool ok = fread(hdr, sizeof hdr, 1, f) == 1 && fread(lh, sizeof lh, 1, f) == 1;
  if (ok) {
    d.K = hdr[0]; d.kb = hdr[1]; d.degree = hdr[2]; d.guard = hdr[3];
    d.lam_hi = lh[0]; d.theta_top = lh[1];
    d.theta.resize(d.kb); d.X.resize((size_t)d.K * d.kb); d.W.resize((size_t)d.K * d.K);
    ok = fread(d.theta.data(), 4, d.kb, f) == (size_t)d.kb &&
         fread(d.X.data(), 4, d.X.size(), f) == d.X.size() &&
         fread(d.W.data(), 4, d.W.size(), f) == d.W.size();
  }
  fclose(f);
  return ok;
}

static double median(std::vector<double> v)
{
  std::sort(v.begin(), v.end());
  return v.empty() ? 0 : v[v.size() / 2];
}

static double cheb_abs(int d, double t)   // |C_d(t)| for |t| >= 1, else 1
{
  return cosh(d * acosh(std::max(1.0, fabs(t))));
}

static const float one = 1, zero = 0;

// Chebyshev filter with per-column degrees deg[i] (non-decreasing in i, each
// in 1..D): step d runs on the trailing columns with deg >= d. X is the
// input block (overwritten as scratch), the result lands in Yout. Returns
// the number of product columns.
static long filter(const float *W, li K, li kb, const std::vector<int> &deg, float c, float h,
                   float *X, float *Yout, float *T)
{
  const int D = deg.empty() ? 1 : *std::max_element(deg.begin(), deg.end());
  std::vector<li> j(D + 2, kb);   // j[d] = first column with deg >= d
  for (int d = 1; d <= D; d++) {
    li q = 0;
    while (q < kb && deg[q] < d) q++;
    j[d] = q;
  }
  float *A = X, *B = Yout;   // A = Y_{d-1}, B = Y_d
  long cols = kb;
  F77(sgemm)("N", "N", &K, &kb, &K, &one, W, &K, A, &K, &zero, T, &K, 1, 1);
  for (size_t t = 0; t < (size_t)K * kb; t++) B[t] = (T[t] - c * A[t]) / h;
  for (int d = 2; d <= D; d++) {
    const li j0 = j[d], n = kb - j0;
    if (n <= 0) break;
    const size_t o = (size_t)j0 * K, nb = (size_t)n * K;
    F77(sgemm)("N", "N", &K, &n, &K, &one, W, &K, B + o, &K, &zero, T + o, &K, 1, 1);
    for (size_t t = 0; t < nb; t++) A[o + t] = 2 * (T[o + t] - c * B[o + t]) / h - A[o + t];
    std::swap(A, B);
    cols += n;
  }
  // columns finished at degree d sit in B if d has D's parity, else in A
  if (B != Yout) {   // D even: the final columns are in X; move everything to Yout
    for (int d = 1; d <= D; d++) {
      const li j0 = j[d], j1 = (d < D) ? j[d + 1] : kb;
      if (j1 <= j0) continue;
      const float *src = ((D - d) % 2 == 0) ? B : A;
      if (src != Yout) memcpy(Yout + (size_t)j0 * K, src + (size_t)j0 * K, (size_t)(j1 - j0) * K * 4);
    }
  } else {
    for (int d = 1; d < D; d++) {
      const li j0 = j[d], j1 = j[d + 1];
      if (j1 <= j0 || (D - d) % 2 == 0) continue;
      memcpy(Yout + (size_t)j0 * K, A + (size_t)j0 * K, (size_t)(j1 - j0) * K * 4);
    }
  }
  return cols;
}

static bool cholqr(float *Y, li K, li kb, bool normalize)
{
  if (normalize)
    for (int i = 0; i < kb; i++) {
      float *y = Y + (size_t)i * K;
      double s = 0;
      for (int t = 0; t < K; t++) s += (double)y[t] * y[t];
      const float f = (float)(1 / sqrt(s + 1e-300));
      for (int t = 0; t < K; t++) y[t] *= f;
    }
  std::vector<float> G((size_t)kb * kb);
  li info = 0;
  F77(ssyrk)("U", "T", &kb, &K, &one, Y, &K, &zero, G.data(), &kb, 1, 1);
  F77(spotrf)("U", &kb, G.data(), &kb, &info, 1);
  if (info) return false;
  F77(strsm)("R", "U", "N", "N", &K, &kb, &one, G.data(), &kb, Y, &K, 1, 1, 1, 1);
  return true;
}

struct Exact {
  int nneg = 0;
  std::vector<double> Vs;   // K x nneg: v_i |lambda_i|^1/2
  double n2 = 0;            // ||P_-(W)||_F^2
};

// ||F F' - N||_F / ||N||_F, N = -P_-(W) = Vs Vs'
static double err_exact(const Exact &e, const double *F, li K, li r)
{
  if (r <= 0) return 1;
  std::vector<double> FtF((size_t)r * r), FtV((size_t)r * e.nneg);
  const double done = 1, dzero = 0;
  li nn = e.nneg;
  F77(dgemm)("T", "N", &r, &r, &K, &done, F, &K, F, &K, &dzero, FtF.data(), &r, 1, 1);
  F77(dgemm)("T", "N", &r, &nn, &K, &done, F, &K, e.Vs.data(), &K, &dzero, FtV.data(), &r, 1, 1);
  double a = 0, b = 0;
  for (double x : FtF) a += x * x;
  for (double x : FtV) b += x * x;
  return sqrt(std::max(0.0, a + e.n2 - 2 * b) / e.n2);
}

struct RRWork {
  std::vector<float> G, H2, S, Zp, theta, wp, Xn, YP, work_d, work_r, work_p;
  std::vector<li> iwork_d, iwork_r, isuppz, piv;
  li lwork_d, liwork_d, lwork_r = -1, liwork_r = -1;
  void init(li K, li kb)
  {
    size_t k2 = (size_t)kb * kb, nb = (size_t)K * kb;
    G.resize(k2); H2.resize(k2); S.resize(k2); Zp.resize(k2); theta.resize(kb); wp.resize(kb);
    Xn.resize(nb); YP.resize(nb); isuppz.resize(2 * (size_t)kb); piv.resize(kb);
    lwork_d = 1 + 6 * kb + 2 * kb * kb; liwork_d = 3 + 5 * kb;
    work_d.resize(lwork_d); iwork_d.resize(liwork_d);
    float wq; li iq, m, info;
    const float vl = 0, vu = FLT_MAX, abstol = 0;
    const li il = 1, iu = kb;
    lwork_r = liwork_r = -1;
    F77(ssyevr)("V", "V", "L", &kb, H2.data(), &kb, &vl, &vu, &il, &iu, &abstol, &m, wp.data(),
                Zp.data(), &kb, isuppz.data(), &wq, &lwork_r, &iq, &liwork_r, &info, 1, 1, 1);
    lwork_r = (li)wq + 1; liwork_r = iq + 1;
    work_r.resize(lwork_r); iwork_r.resize(liwork_r); work_p.resize(2 * (size_t)kb);
  }
};

static void form_G(const float *Y, const float *T, li K, li kb, float *G)
{
  F77(sgemm)("T", "N", &kb, &kb, &K, &one, Y, &K, T, &K, &zero, G, &kb, 1, 1);
  for (int i = 0; i < kb; i++)
    for (int j = i + 1; j < kb; j++)
      G[(size_t)i * kb + j] = G[(size_t)j * kb + i] = 0.5f * (G[(size_t)i * kb + j] + G[(size_t)j * kb + i]);
}

// today's Rayleigh-Ritz from G (destroyed): F (fp64, K x nneg); times in t[3]
static int rr_today(RRWork &r, const float *Y, li K, li kb, double *F, double *t)
{
  li info;
  double s1 = cpu();
  F77(ssyevd)("V", "L", &kb, r.G.data(), &kb, r.theta.data(), r.work_d.data(), &r.lwork_d,
              r.iwork_d.data(), &r.liwork_d, &info, 1, 1);
  double s2 = cpu();
  F77(sgemm)("N", "N", &K, &kb, &kb, &one, Y, &K, r.G.data(), &kb, &zero, r.Xn.data(), &K, 1, 1);
  double s3 = cpu();
  int nneg = 0;
  while (nneg < kb && r.theta[nneg] < 0) nneg++;
  for (int i = 0; i < nneg; i++) {
    const double sq = sqrt(-(double)r.theta[i]);
    const float *x = r.Xn.data() + (size_t)i * K;
    double *f = F + (size_t)i * K;
    for (int q = 0; q < K; q++) f[q] = sq * x[q];
  }
  double s4 = cpu();
  t[0] = s2 - s1; t[1] = s3 - s2; t[2] = s4 - s3;
  return nneg;
}

// F = Y(:,piv) L from the pivoted Cholesky of S (in r.S, lower); returns the rank
static int factor_F(RRWork &r, const float *Y, li K, li kb, double *F, double *t)
{
  double s0 = cpu();
  const float tol = -1;   // LAPACK default: kb * eps * max diag
  li rk = 0, info = 0;
  F77(spstrf)("L", &kb, r.S.data(), &kb, r.piv.data(), &rk, &tol, r.work_p.data(), &info, 1);
  double s1 = cpu();
  for (int i = 0; i < kb; i++)
    memcpy(r.YP.data() + (size_t)i * K, Y + (size_t)(r.piv[i] - 1) * K, (size_t)K * 4);
  if (rk > 0) {
    F77(strmm)("R", "L", "N", "N", &K, &rk, &one, r.S.data(), &kb, r.YP.data(), &K, 1, 1, 1, 1);
    li rest = kb - rk;
    if (rest > 0)
      F77(sgemm)("N", "N", &K, &rk, &rest, &one, r.YP.data() + (size_t)rk * K, &K, r.S.data() + rk, &kb, &one,
                 r.YP.data(), &K, 1, 1);
  }
  double s2 = cpu();
  for (size_t q = 0; q < (size_t)K * rk; q++) F[q] = r.YP[q];
  double s3 = cpu();
  t[0] = s1 - s0; t[1] = s2 - s1; t[2] = s3 - s2;
  return rk;
}

// CS2 from G (kept): positive Ritz pairs only; times in t[5]; *m = positive count
static int rr_cs2(RRWork &r, const float *Y, li K, li kb, double *F, double *t, int *mpos)
{
  double s1 = cpu();
  memcpy(r.H2.data(), r.G.data(), (size_t)kb * kb * 4);
  li m = 0, info;
  const float vl = 0, vu = FLT_MAX, abstol = 0;
  const li il = 1, iu = kb;
  F77(ssyevr)("V", "V", "L", &kb, r.H2.data(), &kb, &vl, &vu, &il, &iu, &abstol, &m, r.wp.data(),
              r.Zp.data(), &kb, r.isuppz.data(), r.work_r.data(), &r.lwork_r, r.iwork_r.data(),
              &r.liwork_r, &info, 1, 1, 1);
  double s2 = cpu();
  for (int j = 0; j < m; j++) {
    const float sq = sqrtf(r.wp[j]);
    float *z = r.Zp.data() + (size_t)j * kb;
    for (int i = 0; i < kb; i++) z[i] *= sq;
  }
  if (m > 0) F77(ssyrk)("L", "N", &kb, &m, &one, r.Zp.data(), &kb, &zero, r.S.data(), &kb, 1, 1);
  else memset(r.S.data(), 0, (size_t)kb * kb * 4);
  for (int j = 0; j < kb; j++)
    for (int i = j; i < kb; i++) r.S[(size_t)j * kb + i] -= r.G[(size_t)j * kb + i];
  double s3 = cpu();
  int rk = factor_F(r, Y, K, kb, F, t + 2);
  t[0] = s2 - s1; t[1] = s3 - s2;
  *mpos = m;
  return rk;
}

// no eigensolve at all: S = sigma I - G (PSD when sigma >= lambda_max(G)), pivoted
// Cholesky; R then carries a PSD bias Y (sigma I - P_+(G)) Y' of norm <= sigma
static int rr_shift(RRWork &r, const float *Y, li K, li kb, double *F, double *t, float sigma)
{
  double s1 = cpu();
  for (int j = 0; j < kb; j++)
    for (int i = j; i < kb; i++) r.S[(size_t)j * kb + i] = -r.G[(size_t)j * kb + i] + (i == j ? sigma : 0.0f);
  double s2 = cpu();
  int rk = factor_F(r, Y, K, kb, F, t + 1);
  t[0] = s2 - s1;
  return rk;
}

int main(int argc, char **argv)
{
  int reps = 5;
  bool exact = false;
  int a = 1;
  while (a < argc && argv[a][0] == '-') {
    if (!strcmp(argv[a], "-r") && a + 1 < argc) { reps = atoi(argv[a + 1]); a += 2; }
    else if (!strcmp(argv[a], "-x")) { exact = true; a++; }
    else { fprintf(stderr, "usage: rr_bench [-r reps] [-x] dump.bin ...\n"); return 1; }
  }
  for (; a < argc; a++) {
    Dump d;
    if (!load(argv[a], d)) { printf("%s: cannot read\n", argv[a]); continue; }
    const li K = d.K, kb = d.kb;
    const size_t nb = (size_t)K * kb;
    const int D = d.degree;
    const float c = (float)(0.5 * (d.lam_hi + d.theta_top)), h = (float)(0.5 * (d.lam_hi - d.theta_top));
    const float *W = d.W.data();
    int nprev = 0;
    while (nprev < kb && d.theta[nprev] < 0) nprev++;
    printf("%s: K %d kb %d degree %d guard %d lam_hi %.4g theta_top %.4g theta_min %.4g (prev nneg %d)\n",
           argv[a], K, kb, D, d.guard, d.lam_hi, d.theta_top, (double)d.theta[0], nprev);

    Exact ex;
    if (exact) {
      std::vector<double> Wd(d.W.begin(), d.W.end()), lam(K);
      li lw = 1 + 6 * K + 2 * K * K, liw = 3 + 5 * K, info;
      std::vector<double> wk(lw);
      std::vector<li> iwk(liw);
      F77(dsyevd)("V", "L", &K, Wd.data(), &K, lam.data(), wk.data(), &lw, iwk.data(), &liw, &info, 1, 1);
      while (ex.nneg < K && lam[ex.nneg] < 0) ex.nneg++;
      ex.Vs.resize((size_t)K * ex.nneg);
      for (int i = 0; i < ex.nneg; i++) {
        const double s = sqrt(-lam[i]);
        for (int q = 0; q < K; q++) ex.Vs[(size_t)i * K + q] = s * Wd[(size_t)i * K + q];
        ex.n2 += lam[i] * lam[i];
      }
      int nsmall = 0;
      for (int i = 0; i < K; i++) if (fabs(lam[i]) < 1e-3) nsmall++;
      printf("  exact: nneg %d of K, lambda_min %.4g, lambda_max %.4g, %d eigenvalues within 1e-3 of 0\n",
             ex.nneg, lam[0], lam[K - 1], nsmall);
    }

    RRWork rw;
    rw.init(K, kb);
    std::vector<float> Xc(nb), Y(nb), T(nb);
    std::vector<double> FA(nb), FB(nb), FC(nb);

    // ---- the step's Y (uniform degree, as dnn.C) ----
    std::vector<int> degU(kb, D);
    memcpy(Xc.data(), d.X.data(), nb * 4);
    double t0 = cpu();
    filter(W, K, kb, degU, c, h, Xc.data(), Y.data(), T.data());
    double t_filter = cpu() - t0;
    t0 = cpu();
    if (!cholqr(Y.data(), K, kb, false)) { printf("  CholQR failed\n"); continue; }
    double t_cholqr = cpu() - t0;
    t0 = cpu();
    F77(sgemm)("N", "N", &K, &kb, &K, &one, W, &K, Y.data(), &K, &zero, T.data(), &K, 1, 1);
    double t_rrprod = cpu() - t0;

    std::vector<double> tA[3], tB[5], tC[4], tG, tR, tT, tX, tD;
    int mX = 0;
    int nnegA = 0, rankB = 0, mB = 0, rankC = 0;
    const float sigma = (float)std::max(0.0, d.theta_top);
    for (int rep = 0; rep < reps; rep++) {
      double tt[5];
      double g0 = cpu();
      form_G(Y.data(), T.data(), K, kb, rw.G.data());
      tG.push_back(cpu() - g0);
      // CS2 and shift first (they keep G), today last (destroys G)
      rankB = rr_cs2(rw, Y.data(), K, kb, FB.data(), tt, &mB);
      for (int q = 0; q < 5; q++) tB[q].push_back(tt[q]);
      rankC = rr_shift(rw, Y.data(), K, kb, FC.data(), tt, sigma);
      for (int q = 0; q < 4; q++) tC[q].push_back(tt[q]);
      {   // full ssyevr (MRRR, all pairs) as a drop-in for today's ssyevd
        memcpy(rw.H2.data(), rw.G.data(), (size_t)kb * kb * 4);
        li m = 0, info;
        const float vl = 0, vu = 0, abstol = 0;
        const li il = 1, iu = kb;
        double r0 = cpu();
        F77(ssyevr)("V", "A", "L", &kb, rw.H2.data(), &kb, &vl, &vu, &il, &iu, &abstol, &m, rw.wp.data(),
                    rw.Zp.data(), &kb, rw.isuppz.data(), rw.work_r.data(), &rw.lwork_r, rw.iwork_r.data(),
                    &rw.liwork_r, &info, 1, 1, 1);
        tR.push_back(cpu() - r0);
      }
      {   // the tridiagonal reduction alone: the floor of any ssyev* on G
        memcpy(rw.H2.data(), rw.G.data(), (size_t)kb * kb * 4);
        std::vector<float> dd(kb), ee(kb), tau(kb), wk(128 * (size_t)kb);
        li lw = 128 * kb, info;
        double r0 = cpu();
        F77(ssytrd)("L", &kb, rw.H2.data(), &kb, dd.data(), ee.data(), tau.data(), wk.data(), &lw, &info, 1);
        tT.push_back(cpu() - r0);
      }
      {   // other partial eigensolvers for the positive pairs: ssyevx (bisection + inverse
          // iteration), and ssytrd + sstedc (all tridiagonal vectors) + sormtr on the m positive only
        memcpy(rw.H2.data(), rw.G.data(), (size_t)kb * kb * 4);
        li m = 0, info;
        const float vl = 0, vu = FLT_MAX, abstol = 0;
        const li il = 1, iu = kb, lw = 128 * kb;
        std::vector<float> wk(lw);
        std::vector<li> iwk(5 * (size_t)kb), ifail(kb);
        double r0 = cpu();
        F77(ssyevx)("V", "V", "L", &kb, rw.H2.data(), &kb, &vl, &vu, &il, &iu, &abstol, &m, rw.wp.data(),
                    rw.Zp.data(), &kb, wk.data(), &lw, iwk.data(), ifail.data(), &info, 1, 1, 1);
        tX.push_back(cpu() - r0);
        mX = m;
        memcpy(rw.H2.data(), rw.G.data(), (size_t)kb * kb * 4);
        std::vector<float> dd(kb), ee(kb), tau(kb), Zt((size_t)kb * kb);
        r0 = cpu();
        F77(ssytrd)("L", &kb, rw.H2.data(), &kb, dd.data(), ee.data(), tau.data(), wk.data(), &lw, &info, 1);
        li lwc = 1 + 4 * kb + kb * kb, liwc = 3 + 5 * kb;
        std::vector<float> wc(lwc);
        std::vector<li> iwc(liwc);
        F77(sstedc)("I", &kb, dd.data(), ee.data(), Zt.data(), &kb, wc.data(), &lwc, iwc.data(), &liwc, &info, 1);
        li first = 0;
        while (first < kb && dd[first] <= 0) first++;
        li mm = kb - first;
        if (mm > 0)
          F77(sormtr)("L", "L", "N", &kb, &mm, rw.H2.data(), &kb, tau.data(), Zt.data() + (size_t)first * kb, &kb,
                      wk.data(), &lw, &info, 1, 1, 1);
        tD.push_back(cpu() - r0);
      }
      nnegA = rr_today(rw, Y.data(), K, kb, FA.data(), tt);
      for (int q = 0; q < 3; q++) tA[q].push_back(tt[q]);
    }
    auto ms = [](const std::vector<double> &v) { return median(v) * 1e3; };
    double A_tot = ms(tA[0]) + ms(tA[1]) + ms(tA[2]);
    double B_tot = 0, C_tot = 0;
    for (int q = 0; q < 5; q++) B_tot += ms(tB[q]);
    for (int q = 0; q < 4; q++) C_tot += ms(tC[q]);
    printf("  step (ms): filter %.1f  cholqr %.1f  rr product %.1f  G=Y'T+sym %.1f\n", t_filter * 1e3,
           t_cholqr * 1e3, t_rrprod * 1e3, ms(tG));
    printf("  today (ms): ssyevd %.1f  X=YG %.1f  F %.1f  = %.1f   (rank nneg %d)\n", ms(tA[0]), ms(tA[1]),
           ms(tA[2]), A_tot, nnegA);
    printf("  CS2   (ms): ssyevr+ %.1f (m %d)  S %.1f  spstrf %.1f  F=Y(:,p)L %.1f  F64 %.1f  = %.1f   (rank %d)"
           "  saved %.1f\n", ms(tB[0]), mB, ms(tB[1]), ms(tB[2]), ms(tB[3]), ms(tB[4]), B_tot, rankB, A_tot - B_tot);
    printf("  shift (ms): S %.1f  spstrf %.1f  F=Y(:,p)L %.1f  F64 %.1f  = %.1f   (rank %d, sigma %.3g)  saved %.1f"
           " before the rank-%d dsyrk's extra columns\n", ms(tC[0]), ms(tC[1]), ms(tC[2]), ms(tC[3]), C_tot, rankC,
           sigma, A_tot - C_tot, rankC);
    printf("  full ssyevr (all pairs) %.1f ms against ssyevd %.1f; ssytrd alone %.1f\n", ms(tR), ms(tA[0]), ms(tT));
    printf("  positive pairs only: ssyevr %.1f  ssyevx %.1f (m %d)  ssytrd+sstedc+sormtr(m) %.1f\n", ms(tB[0]), ms(tX), mX, ms(tD));
    printf("  SUMMARY kb %d nneg %d m %d today %.1f cs2 %.1f saved %.1f\n", (int)kb, nnegA, mB, A_tot, B_tot, A_tot - B_tot);
    if (exact)
      printf("  ||R - R_exact||_F / ||P_-(W)||_F: today %.3e  CS2 %.3e  shift %.3e\n",
             err_exact(ex, FA.data(), K, nnegA), err_exact(ex, FB.data(), K, rankB),
             err_exact(ex, FC.data(), K, rankC));

    // ---- per-column degree rules ----
    // t_i from the dumped (previous) Ritz values; guard, positive and appended columns keep D
    std::vector<double> tcol(kb);
    for (int i = 0; i < kb; i++) tcol[i] = (d.theta[i] - c) / h;
    const double te = nprev > 0 ? tcol[nprev - 1] : -1;
    struct Rule { std::string name; std::vector<int> deg; };
    std::vector<Rule> rules;
    rules.push_back({"uniform " + std::to_string(D), degU});
    rules.push_back({"uniform " + std::to_string(D - 1), std::vector<int>(kb, std::max(1, D - 1))});
    rules.push_back({"uniform " + std::to_string(D + 1), std::vector<int>(kb, D + 1)});
    {
      // memo: d_i = min{d : |C_d(t_i)| >= |C_D(t_e)|}
      const double target = cheb_abs(D, te);
      std::vector<int> g(kb, D);
      for (int i = 0; i < nprev; i++) {
        int di = 1;
        while (di < D && cheb_abs(di, tcol[i]) < target) di++;
        g[i] = di;
      }
      for (int i = 1; i < kb; i++) g[i] = std::max(g[i], g[i - 1]);
      rules.push_back({"memo (edge damping)", g});
    }
    for (int Dt : {D, D + 1})
      for (double rho : {1.5, 2.0, 3.0, 5.0, 10.0}) {
        // amplification floor: d_i = min{d <= Dt : |C_d(t_i)| >= rho}
        std::vector<int> g(kb, Dt);
        for (int i = 0; i < nprev; i++) {
          int di = 1;
          while (di < Dt && cheb_abs(di, tcol[i]) < rho) di++;
          g[i] = di;
        }
        for (int i = 1; i < kb; i++) g[i] = std::max(g[i], g[i - 1]);
        char nm[64];
        snprintf(nm, sizeof nm, "floor rho=%g top %d", rho, Dt);
        rules.push_back({nm, g});
      }
    for (const Rule &ru : rules) {
      memcpy(Xc.data(), d.X.data(), nb * 4);
      std::vector<float> Y2(nb);
      double f0 = cpu();
      long cols = filter(W, K, kb, ru.deg, c, h, Xc.data(), Y2.data(), T.data());
      double tf = cpu() - f0;
      bool ok = cholqr(Y2.data(), K, kb, ru.name.compare(0, 7, "uniform") != 0);
      std::vector<int> hist(D + 2, 0);
      for (int x : ru.deg) hist[std::min(x, D + 1)]++;
      printf("  %-24s cols %.2f kb (filter %.1f ms)  d1 %d d2 %d d3 %d d4 %d", ru.name.c_str(), (double)cols / kb,
             tf * 1e3, hist[1], hist[2], hist[3], D + 1 >= 4 ? hist[4] : 0);
      if (!ok) { printf("  CholQR failed\n"); continue; }
      if (exact) {
        F77(sgemm)("N", "N", &K, &kb, &K, &one, W, &K, Y2.data(), &K, &zero, T.data(), &K, 1, 1);
        form_G(Y2.data(), T.data(), K, kb, rw.G.data());
        double tt[3];
        int nn = rr_today(rw, Y2.data(), K, kb, FA.data(), tt);
        printf("  err %.3e", err_exact(ex, FA.data(), K, nn));
      }
      printf("\n");
    }
    fflush(stdout);
  }
  return 0;
}
