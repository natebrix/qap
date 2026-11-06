/*
  eigbench: Meschach's symmetric eigendecomposition (the solver's my_svd:
  Hfactor + makeHQ + trieig + m_transp) against LAPACK dsyev / dsyevd /
  dsyevr from the scipy-openblas32 wheel, on matrices dumped from a real
  solver run (build_timer.sh with MES_DUMP=<file>).

    eigbench DUMPFILE [REPS]

  Prints, per matrix size n: count, microseconds per call for each method
  (CLOCK_THREAD_CPUTIME_ID, REPS passes over all matrices of that size), and
  the largest eigenvalue difference from Meschach's relative to max|lambda|,
  and the largest eigen-residual ||A v - lambda v|| / ||A||_max.

  Build (from the root of a checkout of 7b49d09, the last commit before the
  eigen port, which still has Meschach, removed from the branch on Nov 6; after make):
    g++ -O3 -march=native -ffp-contract=off -Isrc/MESCHACH \
      docs/meschach_removal/eigbench.C -o eigbench src/MESCHACH/meschach.a \
      -L$OB -lscipy_openblas -Wl,-rpath,$OB -lm
  where OB=$(python3 -c "import scipy_openblas32 as o; print(o.get_lib_dir())")
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <map>
extern "C" {
#include "matrix.h"
#include "matrix2.h"
}
#undef catch
#undef min
#undef max

#define F77(x) scipy_##x##_
typedef int lapack_int;
extern "C" {
void F77(dsyev)(const char *, const char *, const lapack_int *, double *, const lapack_int *,
                double *, double *, const lapack_int *, lapack_int *, size_t, size_t);
void F77(dsyevd)(const char *, const char *, const lapack_int *, double *, const lapack_int *,
                 double *, double *, const lapack_int *, lapack_int *, const lapack_int *,
                 lapack_int *, size_t, size_t);
void F77(dsyevr)(const char *, const char *, const char *, const lapack_int *, double *,
                 const lapack_int *, const double *, const double *, const lapack_int *,
                 const lapack_int *, const double *, lapack_int *, double *, double *,
                 const lapack_int *, lapack_int *, double *, const lapack_int *, lapack_int *,
                 const lapack_int *, lapack_int *, size_t, size_t, size_t);
void scipy_openblas_set_num_threads(int);
}

static double cpu()
{
  timespec t; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

// the solver's my_svd (latest/linalg.C), on a MAT; eigenvalues in d, eigenvectors as rows of A
static void mes_eig(MAT *A, VEC *d)
{
  const unsigned int n = A->m;
  MAT *H = m_copy(A, MNULL);
  MAT *Q = m_get(n, n);
  VEC *diag = v_get(n), *beta = v_get(n), *b = v_get(n - 1);
  Hfactor(H, diag, beta);
  makeHQ(H, diag, beta, Q);
  for (unsigned int i = 0; i < n - 1; i++) { d->ve[i] = H->me[i][i]; b->ve[i] = H->me[i][i + 1]; }
  d->ve[n - 1] = H->me[n - 1][n - 1];
  trieig(d, b, Q);
  m_transp(Q, A);
  V_FREE(diag); V_FREE(beta); V_FREE(b); M_FREE(H); M_FREE(Q);
}

struct Ws {  // LAPACK workspace, sized once per n (as a solver port would keep it)
  int n; std::vector<double> a, z, w, work; std::vector<lapack_int> iwork, isuppz;
  lapack_int lwork_ev, lwork_evd, liwork_evd, lwork_evr, liwork_evr;
  explicit Ws(int n_) : n(n_), a(n_ * n_), z(n_ * n_), w(n_), isuppz(2 * n_) {
    lapack_int N = n, info, lw = -1, liw = -1, iq;
    double q;
    F77(dsyev)("V", "L", &N, a.data(), &N, w.data(), &q, &lw, &info, 1, 1); lwork_ev = (lapack_int)q;
    F77(dsyevd)("V", "L", &N, a.data(), &N, w.data(), &q, &lw, &iq, &liw, &info, 1, 1);
    lwork_evd = (lapack_int)q; liwork_evd = iq;
    double vl = 0, vu = 0, abstol = 0; lapack_int il = 1, iu = n, m;
    F77(dsyevr)("V", "A", "L", &N, a.data(), &N, &vl, &vu, &il, &iu, &abstol, &m, w.data(), z.data(),
                &N, isuppz.data(), &q, &lw, &iq, &liw, &info, 1, 1, 1);
    lwork_evr = (lapack_int)q; liwork_evr = iq;
    work.resize(std::max(lwork_ev, std::max(lwork_evd, lwork_evr)));
    iwork.resize(std::max(liwork_evd, liwork_evr));
  }
};

enum { M_MES, M_EV, M_EVD, M_EVR, M_N };
static const char *mname[M_N] = {"meschach", "dsyev", "dsyevd", "dsyevr"};

// run method on matrix src (row major n x n); eigenvalues ascending into lam,
// eigenvectors as rows into vec (vec[k*n + i] = component i of vector k)
static void run(int meth, int n, const double *src, Ws &ws, MAT *Am, VEC *dm,
                std::vector<double> &lam, std::vector<double> &vec, bool collect)
{
  lapack_int N = n, info;
  if (meth == M_MES) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Am->me[i][j] = src[i * n + j];
    mes_eig(Am, dm);
    if (collect) {
      std::vector<int> ord(n); for (int i = 0; i < n; i++) ord[i] = i;
      std::sort(ord.begin(), ord.end(), [&](int x, int y) { return dm->ve[x] < dm->ve[y]; });
      for (int k = 0; k < n; k++) { lam[k] = dm->ve[ord[k]]; for (int i = 0; i < n; i++) vec[k * n + i] = Am->me[ord[k]][i]; }
    }
    return;
  }
  // symmetric: row major == column major
  std::copy(src, src + n * n, ws.a.begin());
  if (meth == M_EV) {
    F77(dsyev)("V", "L", &N, ws.a.data(), &N, ws.w.data(), ws.work.data(), &ws.lwork_ev, &info, 1, 1);
    if (collect) for (int k = 0; k < n; k++) { lam[k] = ws.w[k]; for (int i = 0; i < n; i++) vec[k * n + i] = ws.a[k * n + i]; }
  } else if (meth == M_EVD) {
    F77(dsyevd)("V", "L", &N, ws.a.data(), &N, ws.w.data(), ws.work.data(), &ws.lwork_evd,
                ws.iwork.data(), &ws.liwork_evd, &info, 1, 1);
    if (collect) for (int k = 0; k < n; k++) { lam[k] = ws.w[k]; for (int i = 0; i < n; i++) vec[k * n + i] = ws.a[k * n + i]; }
  } else {
    double vl = 0, vu = 0, abstol = 0; lapack_int il = 1, iu = n, m;
    F77(dsyevr)("V", "A", "L", &N, ws.a.data(), &N, &vl, &vu, &il, &iu, &abstol, &m, ws.w.data(),
                ws.z.data(), &N, ws.isuppz.data(), ws.work.data(), &ws.lwork_evr, ws.iwork.data(),
                &ws.liwork_evr, &info, 1, 1, 1);
    if (collect) for (int k = 0; k < n; k++) { lam[k] = ws.w[k]; for (int i = 0; i < n; i++) vec[k * n + i] = ws.z[k * n + i]; }
  }
  if (info != 0) { fprintf(stderr, "%s info %d\n", mname[meth], info); exit(1); }
}

int main(int argc, char **argv)
{
  if (argc < 2) { fprintf(stderr, "usage: eigbench DUMPFILE [REPS]\n"); return 1; }
  int reps = argc > 2 ? atoi(argv[2]) : 5;
  scipy_openblas_set_num_threads(1);
  FILE *f = fopen(argv[1], "rb");
  if (!f) { perror(argv[1]); return 1; }
  std::map<int, std::vector<std::vector<double>>> mats;
  int n;
  while (fread(&n, sizeof n, 1, f) == 1) {
    std::vector<double> a(n * n);
    if (fread(a.data(), sizeof(double), n * n, f) != (size_t)(n * n)) break;
    mats[n].push_back(a);
  }
  fclose(f);
  printf("%4s %6s | %9s %9s %9s %9s | %10s %10s %10s | %10s %10s\n", "n", "count", "mes us", "dsyev",
         "dsyevd", "dsyevr", "dlam ev", "dlam evd", "dlam evr", "res mes", "res evr");
  double tot_t[M_N] = {0}; long tot_c = 0;
  for (auto &kv : mats) {
    n = kv.first;
    auto &list = kv.second;
    if (n < 2) continue;
    Ws ws(n);
    MAT *Am = m_get(n, n); VEC *dm = v_get(n);
    std::vector<double> lam(n), vec(n * n), lam0(n), vec0(n * n);
    double t[M_N];
    for (int meth = 0; meth < M_N; meth++) {
      double t0 = cpu();
      for (int r = 0; r < reps; r++)
        for (auto &a : list) run(meth, n, a.data(), ws, Am, dm, lam, vec, false);
      t[meth] = (cpu() - t0) / (reps * (double)list.size()) * 1e6;
      tot_t[meth] += t[meth] * list.size();
    }
    tot_c += list.size();
    double dl[M_N] = {0}, res_mes = 0, res_evr = 0;
    for (auto &a : list) {
      double amax = 0; for (double x : a) amax = std::max(amax, fabs(x));
      run(M_MES, n, a.data(), ws, Am, dm, lam0, vec0, true);
      double lmax = 0; for (double x : lam0) lmax = std::max(lmax, fabs(x));
      auto resid = [&](const std::vector<double> &L, const std::vector<double> &V) {
        double r = 0;
        for (int k = 0; k < n; k++)
          for (int i = 0; i < n; i++) {
            double s = -L[k] * V[k * n + i];
            for (int j = 0; j < n; j++) s += a[i * n + j] * V[k * n + j];
            r = std::max(r, fabs(s) / (amax > 0 ? amax : 1));
          }
        return r;
      };
      res_mes = std::max(res_mes, resid(lam0, vec0));
      for (int meth = 1; meth < M_N; meth++) {
        run(meth, n, a.data(), ws, Am, dm, lam, vec, true);
        for (int k = 0; k < n; k++) dl[meth] = std::max(dl[meth], fabs(lam[k] - lam0[k]) / (lmax > 0 ? lmax : 1));
        if (meth == M_EVR) res_evr = std::max(res_evr, resid(lam, vec));
      }
    }
    printf("%4d %6zu | %9.2f %9.2f %9.2f %9.2f | %10.2e %10.2e %10.2e | %10.2e %10.2e\n", n, list.size(),
           t[M_MES], t[M_EV], t[M_EVD], t[M_EVR], dl[M_EV], dl[M_EVD], dl[M_EVR], res_mes, res_evr);
    M_FREE(Am); V_FREE(dm);
  }
  printf("weighted by the dump's size mix (%ld matrices): mean us/call  meschach %.2f  dsyev %.2f  dsyevd %.2f  dsyevr %.2f\n",
         tot_c, tot_t[M_MES] / tot_c, tot_t[M_EV] / tot_c, tot_t[M_EVD] / tot_c, tot_t[M_EVR] / tot_c);
  return 0;
}
