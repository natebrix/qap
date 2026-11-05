// Hand check of a beta change under DNN_BETA_SCHED (latest/dnn.C): Z is the
// unscaled multiplier, so a change of beta must leave Z (and Y) as they are.
// On random QAPs of size 4..8 (brute-forced optimum), PRSM:
//   scheduled: DNN_BETA_SCHED=two,f=2,t=21 (beta = 2 x base for iterations
//              1-20, base from 21), 45 iterations;
//   split:     20 iterations at beta = 2 x base, saved, then a same-node
//              restart (warm_primal: Y and Z as saved) at base for 25;
//   scaled:    the same split, but Z multiplied by beta_new/beta_old at the
//              restart, as a scaled dual U = Z/beta would require (negative
//              control: it must differ, or the test cannot see a wrong rescale).
// Exact eigensolvers (dsyevd, ssyevd): scheduled and split must agree to
// 1e-12 in Y and Z and in the best certified bound. Subspace tracking with
// force=1 (refresh at the change): reported. Every bound must be <= the optimum.
// Build from latest/ (after make dnntest):
//   g++ -O3 -march=native -ffp-contract=off -DBNB_STANDALONE -DHAVE_DNN -I. -I../src/MESCHACH \
//     ../experiments/wil50_feasibility/speed/beta_sched/beta_change_check.C qpbnd.o qpbnd2.o lap.o \
//     dlap.o qap.o assign.o link1.o solution.o branch.o bnb.o util.o linalg.o fw.o estimate.o svd.o \
//     sa_qap.o stats.o console.o lgb_train.o state.o dnn.o dnncut.o amx_gemm.o -o beta_change_check \
//     -lm ../src/MESCHACH/meschach.a <the DNN_LAPACK flags of latest/makefile>
#include "dnn.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <vector>

struct Matrix {
  MAT a{}; std::vector<double> values; std::vector<double *> rows;
  explicit Matrix(int n) : values(n * n), rows(n) {
    a.m = a.n = n; for (int i = 0; i < n; i++) rows[i] = values.data() + i * n; a.me = rows.data();
  }
};
static double dist(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size()) return HUGE_VAL;
  double d = 0; for (size_t i = 0; i < a.size(); i++) d = std::max(d, fabs(a[i] - b[i])); return d;
}
static double optimum(const MAT *A, const MAT *B, const MAT *C, double shift) {
  int n = A->m; std::vector<int> p(n); std::iota(p.begin(), p.end(), 0); double best = HUGE_VAL;
  do { double v = shift; for (int i = 0; i < n; i++) { v += C->me[i][p[i]]; for (int j = 0; j < n; j++) v += A->me[i][j] * B->me[p[i]][p[j]]; } best = std::min(best, v); }
  while (std::next_permutation(p.begin(), p.end()));
  return best;
}
int main() {
  dnn_set_blas_threads(1);
  std::mt19937 rng(20261105);
  int cases = 0, unsafe = 0, exact_bad = 0, neg_same = 0;
  double max_exact = 0, max_lb = 0, max_sub = 0, min_neg = HUGE_VAL;
  for (int n = 4; n <= 8; n++) for (int seed = 0; seed < 3; seed++) {
    Matrix A(n), B(n), C(n); double shift = seed - 3;
    for (int i = 0; i < n; i++) for (int j = i; j < n; j++) {
      A.a.me[i][j] = A.a.me[j][i] = (int)(rng() % 16) - 5;
      B.a.me[i][j] = B.a.me[j][i] = (int)(rng() % 13) - 4;
    }
    for (auto &v : C.values) v = (int)(rng() % 11) - 5;
    const double opt = optimum(&A.a, &B.a, &C.a, shift);
    for (int mode = 0; mode <= 2; mode++) {
      DNNParams p; p.eig_mode = mode; p.gamma = 1; p.gamma1 = 1; p.cheb_degree = 3; p.tight_lam = 2;
      p.check_every = 5; p.tol = 0; p.stall_window = 0; p.hopeless_after = 0; p.inner_lap = DNN_LAP_EXACT;
      p.refresh = 40;
      // scheduled
      setenv("DNN_BETA_SCHED", mode == 2 ? "two,f=2,t=21,force=1" : "two,f=2,t=21", 1);
      DNNWarm ws; p.max_iter = 45;
      DNNResult rs = dnn_bound(&A.a, &B.a, &C.a, shift, HUGE_VAL, p, NULL, &ws);
      unsetenv("DNN_BETA_SCHED");
      // split
      DNNWarm w1, w2, w3; DNNParams q = p; q.beta = 2 * p.beta; q.max_iter = 20;
      DNNResult r1 = dnn_bound(&A.a, &B.a, &C.a, shift, HUGE_VAL, q, NULL, &w1);
      q.beta = p.beta; q.max_iter = 25; q.warm_primal = 1;
      DNNResult r2 = dnn_bound(&A.a, &B.a, &C.a, shift, HUGE_VAL, q, &w1, &w2);
      // scaled-dual convention (wrong here)
      DNNWarm w1s = w1; for (double &v : w1s.Z) v *= 0.5;
      DNNResult r3 = dnn_bound(&A.a, &B.a, &C.a, shift, HUGE_VAL, q, &w1s, &w3);
      for (const DNNResult *r : {&rs, &r1, &r2, &r3})
        if (r->status == DNN_ERROR || !(r->lb <= opt + 1e-7)) unsafe++;
      const double d = std::max(dist(ws.Y, w2.Y), dist(ws.Z, w2.Z));
      const double dl = fabs(rs.lb - std::max(r1.lb, r2.lb));
      const double dn = std::max(dist(ws.Y, w3.Y), dist(ws.Z, w3.Z));
      if (mode < 2) {
        max_exact = std::max(max_exact, d); max_lb = std::max(max_lb, dl);
        if (d > 1e-12 || dl > 1e-7) exact_bad++;
      } else max_sub = std::max(max_sub, d);
      min_neg = std::min(min_neg, dn); if (dn <= 1e-6) neg_same++;
      printf("n %d seed %d eig %d opt %.4f lb sched %.6f split %.6f/%.6f scaled %.6f |sched-split| %.3g |sched-scaled| %.3g\n",
             n, seed, mode, opt, rs.lb, r1.lb, r2.lb, r3.lb, d, dn);
      cases++;
    }
  }
  printf("%s cases=%d unsafe=%d exact_mismatch=%d max_exact_state_diff=%.3g max_exact_lb_diff=%.3g "
         "max_subspace_state_diff=%.3g negative_control_min_diff=%.3g negative_control_insensitive=%d\n",
         unsafe == 0 && exact_bad == 0 && neg_same == 0 ? "PASS" : "FAIL", cases, unsafe, exact_bad,
         max_exact, max_lb, max_sub, min_neg, neg_same);
  return unsafe == 0 && exact_bad == 0 && neg_same == 0 ? 0 : 1;
}
