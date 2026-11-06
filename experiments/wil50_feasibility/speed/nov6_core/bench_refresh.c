// The subspace tracker's full refresh: ssyevd(K) (all K eigenvectors) against the
// same tridiagonalisation and divide-and-conquer with the back-transformation
// (sormtr) applied to the lowest `nwant` columns and the top one only, which is
// all the tracker reads (nneg + guard negative-side vectors and the top eigenvector
// for tight lambda). Build: see build_bench.sh. Usage: bench_refresh [K] [nwant] [reps]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define F(n) scipy_##n##_
extern void F(ssyevd)(const char*, const char*, const int*, float*, const int*, float*, float*, const int*, int*, const int*, int*, int, int);
extern void F(ssytrd)(const char*, const int*, float*, const int*, float*, float*, float*, float*, const int*, int*, int);
extern void F(sstedc)(const char*, const int*, float*, float*, float*, const int*, float*, const int*, int*, const int*, int*, int);
extern void F(sormtr)(const char*, const char*, const char*, const int*, const int*, const float*, const int*, const float*, float*, const int*, float*, const int*, int*, int, int, int);
extern void F(sstevd)(const char*, const int*, float*, float*, float*, const int*, float*, const int*, int*, const int*, int*, int);
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(int argc, char **argv)
{
  int K = argc > 1 ? atoi(argv[1]) : 2117, nwant = argc > 2 ? atoi(argv[2]) : 700, reps = argc > 3 ? atoi(argv[3]) : 3;
  size_t kk = (size_t)K * K;
  float *W = malloc(kk * 4), *A = malloc(kk * 4), *Z = malloc(kk * 4), *ev = malloc(K * 4), *ev2 = malloc(K * 4);
  float *d = malloc(K * 4), *e = malloc(K * 4), *tau = malloc(K * 4);
  srand(7);
  // an indefinite matrix with a spectrum like the solver's W: about a third negative
  for (int i = 0; i < K; i++) for (int j = 0; j <= i; j++) { float v = ((float)rand() / RAND_MAX - 0.5f) / sqrtf((float)K); if (i == j) v += (i < K / 3 ? -1.0f : 1.0f) * ((float)rand() / RAND_MAX); W[i + (size_t)j * K] = W[j + (size_t)i * K] = v; }
  int lwork = 1 + 6 * K + 2 * K * K, liwork = 3 + 5 * K, info = 0;
  float *work = malloc((size_t)lwork * 4); int *iwork = malloc((size_t)liwork * 4);
  double t_evd = 1e9, t_trd = 1e9, t_edc = 1e9, t_mtr_all = 1e9, t_mtr = 1e9;
  for (int r = 0; r < reps; r++) {
    memcpy(A, W, kk * 4);
    double t0 = now();
    F(ssyevd)("V", "L", &K, A, &K, ev, work, &lwork, iwork, &liwork, &info, 1, 1);
    double t = now() - t0; if (t < t_evd) t_evd = t;
    if (info) { printf("ssyevd info %d\n", info); return 1; }
  }
  for (int r = 0; r < reps; r++) {
    memcpy(A, W, kk * 4);
    double t0 = now();
    F(ssytrd)("L", &K, A, &K, d, e, tau, work, &lwork, &info, 1);
    double t1 = now();
    // eigenvectors of the tridiagonal matrix into Z (K x K), divide and conquer
    F(sstedc)("I", &K, d, e, Z, &K, work, &lwork, iwork, &liwork, &info, 1);
    double t2 = now();
    // back-transform all K columns (what ssyevd does)
    memcpy(ev2, d, K * 4);
    float *Zall = malloc(kk * 4); memcpy(Zall, Z, kk * 4);
    F(sormtr)("L", "L", "N", &K, &K, A, &K, tau, Zall, &K, work, &lwork, &info, 1, 1, 1);
    double t3 = now();
    // back-transform the lowest nwant columns and the top column only
    int nw = nwant, one = 1;
    F(sormtr)("L", "L", "N", &K, &nw, A, &K, tau, Z, &K, work, &lwork, &info, 1, 1, 1);
    F(sormtr)("L", "L", "N", &K, &one, A, &K, tau, Z + (size_t)(K - 1) * K, &K, work, &lwork, &info, 1, 1, 1);
    double t4 = now();
    if (t1 - t0 < t_trd) t_trd = t1 - t0;
    if (t2 - t1 < t_edc) t_edc = t2 - t1;
    if (t3 - t2 < t_mtr_all) t_mtr_all = t3 - t2;
    if (t4 - t3 < t_mtr) t_mtr = t4 - t3;
    // agreement with ssyevd on the wanted columns (sign-invariant)
    if (r == reps - 1) {
      double maxd = 0, maxv = 0;
      for (int i = 0; i < K; i++) maxv = fmax(maxv, fabs(ev[i] - ev2[i]));
      for (int j = 0; j < nwant; j++) {
        double s = 0; for (int i = 0; i < K; i++) s += (double)A[i + (size_t)j * K] * 0; // placeholder (A holds the reflectors now)
        double dp = 0; for (int i = 0; i < K; i++) dp += (double)Zall[i + (size_t)j * K] * Z[i + (size_t)j * K];
        double sg = dp >= 0 ? 1 : -1, dd = 0;
        for (int i = 0; i < K; i++) dd = fmax(dd, fabs(Zall[i + (size_t)j * K] - sg * Z[i + (size_t)j * K]));
        maxd = fmax(maxd, dd); (void)s;
      }
      printf("eigenvalues: max |ssyevd - ssytrd/sstedc| = %.2e; wanted vectors: max |all-column sormtr - partial sormtr| = %.2e\n", maxv, maxd);
    }
    free(Zall);
  }
  printf("K %d nwant %d (min of %d reps)\n", K, nwant, reps);
  printf("ssyevd                 %7.1f ms\n", t_evd * 1e3);
  printf("ssytrd                 %7.1f ms\n", t_trd * 1e3);
  printf("sstedc (T eigvecs)     %7.1f ms\n", t_edc * 1e3);
  printf("sormtr, all K columns  %7.1f ms   -> ssytrd+sstedc+sormtr(all) = %.1f ms\n", t_mtr_all * 1e3, (t_trd + t_edc + t_mtr_all) * 1e3);
  printf("sormtr, %d + 1 columns %6.1f ms   -> ssytrd+sstedc+sormtr(part) = %.1f ms (%.0f%% of ssyevd)\n", nwant, t_mtr * 1e3, (t_trd + t_edc + t_mtr) * 1e3, 100 * (t_trd + t_edc + t_mtr) / t_evd);
  return 0;
}
