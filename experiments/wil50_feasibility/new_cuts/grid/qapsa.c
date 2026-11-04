/* qapsa: simulated annealing (pair swaps, Connolly-style schedule) for a QAP
   with a linear term:  min sum_{i,j} F[i][j] D[p_i][p_j] + sum_i C[p_i][i]
   (p_i = location of facility i).  Input file (text): m, then D (m x m),
   F (m x m), C (m x m: C[a][i] = cost of facility i at location a).
   Usage: qapsa file runs sweeps seed   -> prints "BEST cost" and the permutation.
   Grid-agent tool (experiments/wil50_feasibility/new_cuts/grid); heuristic
   upper bounds only, never used in a certified value. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int m;
static double *D, *F, *C;
#define d(a,b) D[(a)*m+(b)]
#define f(i,j) F[(i)*m+(j)]
#define c(a,i) C[(a)*m+(i)]

static double cost(const int *p) {
  double s = 0;
  for (int i = 0; i < m; i++) { s += c(p[i], i); for (int j = 0; j < m; j++) s += f(i,j) * d(p[i], p[j]); }
  return s;
}
/* change of cost when facilities r and s swap locations (F, D symmetric) */
static double delta(const int *p, int r, int s) {
  double g = c(p[s], r) + c(p[r], s) - c(p[r], r) - c(p[s], s);
  int pr = p[r], ps = p[s];
  for (int k = 0; k < m; k++) {
    if (k == r || k == s) continue;
    int pk = p[k];
    g += 2.0 * (f(r,k) - f(s,k)) * (d(ps, pk) - d(pr, pk));
  }
  return g;
}
static unsigned long long rs = 88172645463325252ULL;
static double urand(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (rs >> 11) * (1.0 / 9007199254740992.0); }

int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: qapsa file runs sweeps seed\n"); return 1; }
  FILE *fp = fopen(argv[1], "r"); if (!fp) return 1;
  if (fscanf(fp, "%d", &m) != 1) return 1;
  D = malloc(sizeof(double) * m * m); F = malloc(sizeof(double) * m * m); C = malloc(sizeof(double) * m * m);
  for (int k = 0; k < m * m; k++) if (fscanf(fp, "%lf", &D[k]) != 1) return 1;
  for (int k = 0; k < m * m; k++) if (fscanf(fp, "%lf", &F[k]) != 1) return 1;
  for (int k = 0; k < m * m; k++) if (fscanf(fp, "%lf", &C[k]) != 1) return 1;
  fclose(fp);
  int runs = atoi(argv[2]); long sweeps = atol(argv[3]); rs ^= (unsigned long long)atoll(argv[4]) * 2654435761ULL;
  int *p = malloc(sizeof(int) * m), *best = malloc(sizeof(int) * m);
  double bestc = 1e300;
  for (int run = 0; run < runs; run++) {
    for (int i = 0; i < m; i++) p[i] = i;
    for (int i = m - 1; i > 0; i--) { int j = (int)(urand() * (i + 1)); int t = p[i]; p[i] = p[j]; p[j] = t; }
    double cur = cost(p), runbest = cur;
    /* temperature range from a sample of deltas */
    double dmin = 1e300, dmax = 0;
    for (int k = 0; k < 2000; k++) { int r = (int)(urand() * m), s = (int)(urand() * m); if (r == s) continue;
      double g = fabs(delta(p, r, s)); if (g > 0 && g < dmin) dmin = g; if (g > dmax) dmax = g; }
    double T0 = dmin + (dmax - dmin) / 10.0, Tf = dmin;
    long total = sweeps * (long)m * (m - 1) / 2;
    double beta = (T0 - Tf) / (total * T0 * Tf), T = T0;
    for (long it = 0; it < total; it++) {
      int r = (int)(urand() * m), s = (int)(urand() * (m - 1)); if (s >= r) s++;
      double g = delta(p, r, s);
      if (g < 0 || urand() < exp(-g / T)) {
        int t = p[r]; p[r] = p[s]; p[s] = t; cur += g;
        if (cur < runbest - 1e-9) runbest = cur;
        if (cur < bestc - 1e-9) { bestc = cur; memcpy(best, p, sizeof(int) * m); }
      }
      T = T / (1.0 + beta * T);
    }
    /* steepest-descent polish of the run's end point */
    for (int imp = 1; imp;) { imp = 0;
      for (int r = 0; r < m; r++) for (int s = r + 1; s < m; s++) { double g = delta(p, r, s);
        if (g < -1e-9) { int t = p[r]; p[r] = p[s]; p[s] = t; cur += g; imp = 1; } } }
    if (cur < bestc - 1e-9) { bestc = cur; memcpy(best, p, sizeof(int) * m); }
  }
  double chk = cost(best);
  printf("BEST %.1f (check %.1f)\nPERM", bestc, chk);
  for (int i = 0; i < m; i++) printf(" %d", best[i]);
  printf("\n");
  return 0;
}
