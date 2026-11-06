/*
  eigport_test: sym_eig_port (eig_port.C) against Meschach's path (the
  solver's my_svd), bit for bit and in time, on dumped matrices.

    eigport_test DUMPFILE [REPS]

  Build from the root of a checkout of 7b49d09 (it needs Meschach, removed
  from the branch on Nov 6; run make there first):
    g++ -O3 -march=native -ffp-contract=off -Isrc/MESCHACH \
      docs/meschach_removal/eigport_test.C docs/meschach_removal/eig_port.C \
      -o eigport_test src/MESCHACH/meschach.a -lm
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include <map>
extern "C" {
#include "matrix.h"
#include "matrix2.h"
}
#undef catch
#include "eig_port.h"

static double cpu()
{
  timespec t; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

static void mes_eig(MAT *A, VEC *d)  // latest/linalg.C my_svd
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

int main(int argc, char **argv)
{
  if (argc < 2) { fprintf(stderr, "usage: eigport_test DUMPFILE [REPS]\n"); return 1; }
  int reps = argc > 2 ? atoi(argv[2]) : 20;
  FILE *f = fopen(argv[1], "rb");
  if (!f) { perror(argv[1]); return 1; }
  std::map<int, std::vector<std::vector<double>>> mats;
  int n;
  while (fread(&n, sizeof n, 1, f) == 1) {
    std::vector<double> a((size_t)n * n);
    if (fread(a.data(), sizeof(double), a.size(), f) != a.size()) break;
    mats[n].push_back(a);
  }
  fclose(f);
  long total = 0, identical = 0;
  double tm_all = 0, tp_all = 0;
  printf("%4s %6s %10s %10s %8s %10s\n", "n", "count", "mes us", "port us", "speedup", "identical");
  for (auto &kv : mats) {
    n = kv.first;
    MAT *Am = m_get(n, n); VEC *dm = v_get(n);
    std::vector<double> P((size_t)n * n), dp(n);
    std::vector<double *> rows(n);
    for (int i = 0; i < n; i++) rows[i] = &P[(size_t)i * n];
    long ident = 0;
    for (auto &a : kv.second) {
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Am->me[i][j] = P[i * n + j] = a[i * n + j];
      mes_eig(Am, dm);
      sym_eig_port(n, rows.data(), dp.data());
      bool same = memcmp(dm->ve, dp.data(), n * sizeof(double)) == 0;
      for (int i = 0; i < n && same; i++) same = memcmp(Am->me[i], rows[i], n * sizeof(double)) == 0;
      ident += same;
    }
    double t0 = cpu();
    for (int r = 0; r < reps; r++)
      for (auto &a : kv.second) {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Am->me[i][j] = a[i * n + j];
        mes_eig(Am, dm);
      }
    double tm = cpu() - t0;
    t0 = cpu();
    for (int r = 0; r < reps; r++)
      for (auto &a : kv.second) {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) P[i * n + j] = a[i * n + j];
        sym_eig_port(n, rows.data(), dp.data());
      }
    double tp = cpu() - t0;
    size_t c = kv.second.size();
    printf("%4d %6zu %10.2f %10.2f %8.2f %6ld/%zu\n", n, c, tm / (reps * c) * 1e6, tp / (reps * c) * 1e6,
           tm / tp, ident, c);
    total += c; identical += ident; tm_all += tm; tp_all += tp;
    M_FREE(Am); V_FREE(dm);
  }
  printf("all: %ld/%ld bit-identical (eigenvalues and eigenvectors); time port/meschach %.3f\n",
         identical, total, tp_all / tm_all);
  return identical == total ? 0 : 1;
}
