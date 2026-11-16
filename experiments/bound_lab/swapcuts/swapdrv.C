/*
 swapdrv: exchange-optimality cuts on top of the certified DNN bound (README.md here).
 NOT valid for every permutation: valid for every permutation that no pairwise exchange
 improves, hence for every global optimum, hence usable in a search that proves "no
 permutation costs less than the incumbent".

   swapdrv [-x i:a,...] [-u ub | -c cutoff] [-I first iters] [-er rounds] [-ei iters]
           [-ek cuts per round] [-et theta] [-edelta D] [-eq quiet iters] [-fp]
           [-grid R:C -side S -k K -cpsep k:K -cpiters N]  [-saveL file] [-check N]
           [-joint -r rounds -i iters]   (with -grid: dnn_cut_rounds on the regional triangles,
                                          one exchange Polyak step between rounds; -ek 0: triangles alone)
           [-T threads] [-v] file.dat

 Cuts. Cost C(p) = sum A[k][l] B[p(k)][p(l)]. For a facility i fixed at a, a free facility
 j and a free location b, with x_jb = 1 the gain of exchanging i and j is
   delta = c0 + sum_{free k != j, free c != b} g_k(c) x_kc,
   g_k(c) = (A[i][k]-A[j][k]) (B[b][c]-B[a][c]) + (A[k][i]-A[k][j]) (B[c][b]-B[c][a]),
   c0 = (A[i][j]-A[j][i]) (B[b][a]-B[a][b]) + (A[i][i]-A[j][j]) (B[b][b]-B[a][a])
        + sum over the other fixed k of g_k(p(k)),
 so  c0 Y[jb,jb] + sum g_k(c) Y[jb,kc] >= 0   (type P, one per fixed i, free j, free b).
 -rot adds the 3-cycles: fixed i (at a) to the location a2 of a second fixed facility i2, i2 to
 b, the free j (at b) to a: again linear in row jb (type R, d (d-1) m^2 cuts).
 -fp adds, for two fixed facilities i (at a), i2 (at b): l(x) = c0 + sum g_k(c) x_kc >= 0 on
 row 0 (type F0) and l(x) x_ld >= 0 on every row ld (type F1: diagonal c0 + g_l(d)).

 Phase E: Polyak rounds as dnn_cut_rounds runs them (solve warm with the term
 -sum mu_t G_t in the objective, separate the -ek most violated on the last primal iterate,
 step every pooled multiplier toward min(cutoff, best + delta), halve the step after a
 round that lowered the bound). Then optionally a quiet round at the best round's
 multipliers (-eq), and the regional-triangle CP phase of dnncuts (-grid ... -cpsep ...)
 with the exchange term frozen as DNNCutParams::L_ext (needs swapcuts.patch: the CP form
 ignores L_ext in the stock library).
*/
#include "dnncut.h"
#include "qap.h"
#include "assign.h"
#include "util.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <algorithm>

static bool parse_fixes(const std::string &fixes, QAPAssignment &node, int n)
{
  if (fixes == "-") return true;
  for (size_t pos = 0; pos < fixes.size();) {
    size_t end = fixes.find(',', pos);
    if (end == std::string::npos) end = fixes.size();
    int i, j;
    if (sscanf(fixes.substr(pos, end - pos).c_str(), "%d:%d", &i, &j) != 2 ||
        i < 0 || i >= n || j < 0 || j >= n || node.p[i] != NONE || node.pinv[j] != NONE)
      return false;
    node.fix(i, j);
    pos = end + 1;
  }
  return true;
}

struct XCut {
  int type;      // 0: P (fixed i, free j at b), 1: F0 (two fixed, row 0), 2: F1 (two fixed, row ld),
                 // 3: R (rotation: fixed i -> the location of fixed i2, i2 -> b, free j at b -> a)
  int ns, S[3], from[3], to[3];   // the move: facility S[s] goes from from[s] to to[s] (original indices)
  void move2() { ns = 2; S[0] = i; from[0] = a; to[0] = b; S[1] = j; from[1] = b; to[1] = a; }
  int i, a;      // original: the fixed facility and its location
  int j, b;      // original: the other facility and the location it takes (P: free; F: fixed)
  int row;       // lifted row (0 for F0)
  int exk, exc;  // local facility / location excluded from the sum (-1: none)
  double diag;   // coefficient of Y[row,row] (F0: the constant)
  double norm;   // 2-norm of the coefficients
  double mu;
  bool pooled;
};

struct Ctx {
  int n, m, n1;
  std::vector<double> A, B;            // original n x n (as the solver holds them after symmetrize_QAP)
  std::vector<int> p;                  // original facility -> location or -1
  std::vector<int> ff, fl;             // free facilities / locations, original indices in local order
  std::vector<int> fixf;               // fixed facilities
  double a_(int i, int k) const { return A[(size_t)i * n + k]; }
  double b_(int i, int k) const { return B[(size_t)i * n + k]; }
  // the move's gain is  cpair + sum over the facilities k outside the move of g(k, p(k))
  double g(const XCut &t, int k, int c) const {   // original k, c
    double v = 0;
    for (int s = 0; s < t.ns; s++)
      v += a_(t.S[s], k) * (b_(t.to[s], c) - b_(t.from[s], c)) + a_(k, t.S[s]) * (b_(c, t.to[s]) - b_(c, t.from[s]));
    return v;
  }
  double cpair(const XCut &t) const {
    double v = 0;
    for (int s = 0; s < t.ns; s++) for (int u = 0; u < t.ns; u++)
      v += a_(t.S[s], t.S[u]) * (b_(t.to[s], t.to[u]) - b_(t.from[s], t.from[u]));
    return v;
  }
  // the coefficients over local (k, c), row-major m x m (zero at the excluded row and column)
  void coef(const XCut &t, std::vector<double> &G) const {
    G.assign((size_t)m * m, 0.0);
    for (int s = 0; s < t.ns; s++)
      for (int k = 0; k < m; k++) {
        if (k == t.exk) continue;
        const double as = a_(t.S[s], ff[k]), sa = a_(ff[k], t.S[s]);
        if (as == 0 && sa == 0) continue;
        double *row = G.data() + (size_t)k * m;
        for (int c = 0; c < m; c++) {
          if (c == t.exc) continue;
          row[c] += as * (b_(t.to[s], fl[c]) - b_(t.from[s], fl[c])) + sa * (b_(fl[c], t.to[s]) - b_(fl[c], t.from[s]));
        }
      }
  }
  double cost(const std::vector<int> &q) const {
    double s = 0;
    for (int k = 0; k < n; k++) for (int l = 0; l < n; l++) s += a_(k, l) * b_(q[k], q[l]);
    return s;
  }
};

// lhs of the cut at Y (>= 0 when satisfied)
static double lhs(const Ctx &cx, const XCut &t, const std::vector<double> &Y,
                  std::vector<double> &u1, std::vector<double> &v1, std::vector<double> &u2, std::vector<double> &v2)
{
  const int m = cx.m, n1 = cx.n1;
  cx.coef(t, u1);
  const double *row = Y.data() + (size_t)t.row * n1;
  double s = t.diag * (t.row == 0 ? 1.0 : row[t.row]);
  const double *G = u1.data(), *r = row + 1;
  for (int q = 0; q < m * m; q++) s += G[q] * r[q];
  return s;
}

// min-cost assignment (Hungarian, O(n^3)); cost row-major n x n; returns the minimum
static double lap_min(const std::vector<double> &cost, int n)
{
  const double INF = 1e300;
  std::vector<double> u(n + 1, 0), v(n + 1, 0), minv(n + 1);
  std::vector<int> p(n + 1, 0), way(n + 1, 0);
  std::vector<char> used(n + 1);
  for (int i = 1; i <= n; i++) {
    p[0] = i; int j0 = 0;
    std::fill(minv.begin(), minv.end(), INF); std::fill(used.begin(), used.end(), 0);
    do {
      used[j0] = 1;
      int i0 = p[j0], j1 = 0; double delta = INF;
      for (int j = 1; j <= n; j++) if (!used[j]) {
        double cur = cost[(size_t)(i0 - 1) * n + (j - 1)] - u[i0] - v[j];
        if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
        if (minv[j] < delta) { delta = minv[j]; j1 = j; }
      }
      for (int j = 0; j <= n; j++) { if (used[j]) { u[p[j]] += delta; v[j] -= delta; } else minv[j] -= delta; }
      j0 = j1;
    } while (p[j0] != 0);
    do { int j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
  }
  double s = 0;
  for (int j = 1; j <= n; j++) s += cost[(size_t)(p[j] - 1) * n + (j - 1)];
  return s;
}

// the exchange cuts' pool, multipliers and objective term; one Polyak step per update
struct XState {
  const Ctx *cx;
  std::vector<XCut> *cuts;
  std::vector<double> Ladd, best_mu, sg;
  std::vector<double> u1, v1, u2, v2;
  double abs_add = 0, hsum = 0, theta_cur = 1, delta = 50, prev_lb = -HUGE_VAL, my_best = -HUGE_VAL;
  int ek = 0, verbose = 0, topk = 12, calls = 0;
  double ecap = 2.0, alpha_first = 0, emu = 0.5;
  int assemble() {
    const int m = cx->m, n1 = cx->n1;
    std::fill(Ladd.begin(), Ladd.end(), 0.0);
    abs_add = 0; hsum = 0;
    int active = 0;
    for (auto &t : *cuts) {
      if (!(t.mu > 0)) continue;
      active++;
      cx->coef(t, u1);
      if (t.row == 0) hsum += t.mu * t.diag; else Ladd[(size_t)t.row * n1 + t.row] -= t.mu * t.diag;
      abs_add += t.mu * fabs(t.diag);
      for (int k = 0; k < m; k++) {
        for (int c = 0; c < m; c++) {
          const double gv = u1[(size_t)k * m + c];
          if (gv == 0) continue;
          const int col = 1 + k * m + c;
          Ladd[(size_t)t.row * n1 + col] -= 0.5 * t.mu * gv;
          Ladd[(size_t)col * n1 + t.row] -= 0.5 * t.mu * gv;
          abs_add += t.mu * fabs(gv);
        }
      }
    }
    return active;
  }
  void restore_best() {
    if (best_mu.size() == cuts->size()) for (size_t t = 0; t < cuts->size(); t++) (*cuts)[t].mu = best_mu[t];
    assemble();
  }
  // lb_now: the certified bound of the solve that produced Y, run with the current multipliers
  void update(const std::vector<double> &Y, double lb_now, double best_lb, double cutoff, int r, bool step) {
    const int n1 = cx->n1;
    std::vector<XCut> &C = *cuts;
    if (lb_now > my_best) { my_best = lb_now; best_mu.resize(C.size()); for (size_t t = 0; t < C.size(); t++) best_mu[t] = C[t].mu; }
    sg.resize(C.size());
    std::vector<std::pair<double, int>> viol;
    std::map<int, std::pair<int, double>> byfix;   // fixed facility -> (count, sum of -lhs)
    size_t nv = 0, nv1 = 0, nrot = 0; double mass = 0, worst_cond = 0, worst_norm = 0, mrot = 0;
    for (size_t t = 0; t < C.size(); t++) {
      const XCut &c = C[t];
      const double v = lhs(*cx, c, Y, u1, v1, u2, v2);
      sg[t] = -v;
      if (v < -1e-6) {
        const double yd = c.row == 0 ? 1.0 : Y[(size_t)c.row * n1 + c.row];
        nv++; mass += -v;
        if (c.norm > 0) worst_norm = std::max(worst_norm, -v / c.norm);
        if (c.type == 3) { nrot++; mrot += -v; }
        if (yd > 1e-4) { worst_cond = std::min(worst_cond, v / yd); if (v / yd < -1.0) nv1++; }
        auto &e = byfix[c.type == 0 ? c.i : -1 - c.type];
        e.first++; e.second += -v;
        viol.push_back({c.norm > 0 ? -v / c.norm : 0, (int)t});
      }
    }
    std::sort(viol.begin(), viol.end(), [](const std::pair<double, int> &x, const std::pair<double, int> &y) { return x.first > y.first; });
    printf("  exchange separation %d: %zu of %zu cuts violated (%zu with conditional gain below -1), sum of violations %.2f (cost units x probability), "
           "worst conditional gain %.1f, largest normalised violation %.5f\n", r, nv, C.size(), nv1, mass, worst_cond, worst_norm);
    if (nrot) printf("  of which rotations: %zu violated, sum %.2f\n", nrot, mrot);
    if (verbose || calls == 0) {
      printf("  by fixed facility:");
      for (auto &e : byfix) {
        if (e.first >= 0) printf("  %d@%d: %d (%.2f)", e.first, cx->p[e.first], e.second.first, e.second.second);
        else printf("  fixed-pair type %d: %d (%.2f)", -1 - e.first, e.second.first, e.second.second);
      }
      printf("\n");
      for (int s = 0; s < std::min<int>(topk, (int)viol.size()); s++) {
        const XCut &c = C[viol[s].second];
        const double yd = c.row == 0 ? 1.0 : Y[(size_t)c.row * n1 + c.row];
        printf("    top %2d  type %d  fixed %d@%d  with %d@%d  row %d  y %.5f  lhs %.4f  conditional gain %.2f  normalised %.5f\n", s, c.type,
               c.i, c.a, c.j, c.b, c.row, yd, -sg[viol[s].second], yd > 0 ? -sg[viol[s].second] / yd : 0.0, viol[s].first);
      }
    }
    calls++;
    if (!step || ek == 0) return;   // -ek 0: no exchange cut (the reference arms)
    int added = 0;
    for (auto &e : viol) { if (added >= ek) break; if (!C[e.second].pooled) { C[e.second].pooled = true; added++; } }
    if (lb_now < prev_lb - 1e-9) theta_cur *= 0.5;
    prev_lb = lb_now;
    double s2 = 0; size_t pool = 0;
    for (size_t t = 0; t < C.size(); t++)
      if (C[t].pooled) { pool++; if (C[t].mu > 0 || sg[t] > 0) s2 += sg[t] * sg[t]; }
    if (s2 <= 0) { printf("  no violated exchange cut in the pool\n"); return; }
    const double target = std::min(cutoff, best_lb + delta);
    double alpha = theta_cur * std::max(target - lb_now, 1e-6) / s2;
    // the step cap (-ecap F, default 2): at most F times the first step. When the pooled violations
    // get small the Polyak step explodes (d5_d: 5.8e-2 after 6e-3, the bound fell 3,600)
    if (alpha_first <= 0) alpha_first = alpha;
    else if (ecap > 0 && alpha > ecap * alpha_first) alpha = ecap * alpha_first;
    // the multiplier cap (-emu M, default 0.5): no multiplier moves by more than M in one step. A
    // multiplier is dimensionless (1 = the row's cost replaced by the exchanged cost); where the
    // violations are small (wil50 node 13: 15 units in all) the Polyak step toward best + 50 was
    // 4.3 per unit of violation and the bound fell 331 (logs/s3_wil50_n13_jTX_overshoot.log)
    double maxsg = 0;
    for (size_t t = 0; t < C.size(); t++) if (C[t].pooled) maxsg = std::max(maxsg, sg[t]);
    if (emu > 0 && alpha * maxsg > emu) alpha = emu / maxsg;
    for (size_t t = 0; t < C.size(); t++)
      if (C[t].pooled) C[t].mu = std::max(0.0, C[t].mu + alpha * sg[t]);
    const int active = assemble();
    double mumax = 0;
    for (size_t t = 0; t < C.size(); t++) mumax = std::max(mumax, C[t].mu);
    printf("  exchange: added %d, pool %zu, active %d, step %.3e (theta %.3f, target %.1f), largest multiplier %.3f\n", added, pool, active, alpha, theta_cur, target, mumax);
  }
};
static void xhook(DNNExtHook *h, const std::vector<double> &Y, double lb_now, double best_lb, double cutoff, int round)
{
  XState *xs = (XState *)h->ctx;
  if (round < 0) xs->restore_best();
  else xs->update(Y, lb_now, best_lb, cutoff, round, true);
  h->abs_ = xs->abs_add; h->cnst = xs->hsum;
}

int main(int argc, char **argv)
{
  std::string fixes = "-", saveL;
  double ub = HUGE_VAL, cutoff = HUGE_VAL;
  bool have_cutoff = false, use_fp = false;
  int first_iters = 300, erounds = 0, eiters = 200, ek = 500, equiet = 0, threads = 1, verbose = 0, ncheck = 50;
  int degree = 4, topk = 12, jrounds = 8, jiters = 150;
  bool joint = false, use_elim = false, use_rot = false;
  double theta = 1.0, edelta = 50.0, ecap = 2.0, emu = 0.5;
  DNNCutParams cp;
  cp.kcuts = 1000;
  DNNCPParams cpp;
  int cp_iters = 0;
  const char *file = NULL;
  for (int a = 1; a < argc; a++) {
    std::string o = argv[a];
    bool more = a + 1 < argc;
    if (o == "-x" && more) fixes = argv[++a];
    else if (o == "-u" && more) ub = atof(argv[++a]);
    else if (o == "-c" && more) { cutoff = atof(argv[++a]); have_cutoff = true; }
    else if (o == "-I" && more) first_iters = atoi(argv[++a]);
    else if (o == "-er" && more) erounds = atoi(argv[++a]);
    else if (o == "-ei" && more) eiters = atoi(argv[++a]);
    else if (o == "-ek" && more) ek = atoi(argv[++a]);
    else if (o == "-et" && more) theta = atof(argv[++a]);
    else if (o == "-edelta" && more) edelta = atof(argv[++a]);
    else if (o == "-eq" && more) equiet = atoi(argv[++a]);
    else if (o == "-emu" && more) emu = atof(argv[++a]);     // largest multiplier change per step (0.5; 0: none)
    else if (o == "-ecap" && more) ecap = atof(argv[++a]);   // exchange step <= F x the first step (2; 0: no cap)
    else if (o == "-fp") use_fp = true;
    else if (o == "-rot") use_rot = true;                  // 3-cycles through two fixed facilities and one free
    else if (o == "-elim") use_elim = true;                // forbid x_jb when every completion with j at b gains by an exchange
    else if (o == "-delta" && more) cp.delta = atof(argv[++a]);   // the triangles' Polyak target (joint mode)
    else if (o == "-nocap") cp.cap_at_cutoff = 0;
    else if (o == "-joint") joint = true;                 // the library's triangle rounds with the exchange hook
    else if (o == "-r" && more) jrounds = atoi(argv[++a]);
    else if (o == "-i" && more) jiters = atoi(argv[++a]);
    else if (o == "-top" && more) topk = atoi(argv[++a]);
    else if (o == "-check" && more) ncheck = atoi(argv[++a]);
    else if (o == "-saveL" && more) saveL = argv[++a];
    else if (o == "-grid" && more) {
      if (sscanf(argv[++a], "%d:%d", &cp.grid_rows, &cp.grid_cols) != 2) { fprintf(stderr, "bad -grid\n"); return 1; }
    }
    else if (o == "-prefix") cp.prefix = true;
    else if (o == "-side" && more) cp.side = atoi(argv[++a]);
    else if (o == "-k" && more) cp.kcuts = atoi(argv[++a]);
    else if (o == "-cpiters" && more) cp_iters = atoi(argv[++a]);
    else if (o == "-cpsep" && more) {
      if (sscanf(argv[++a], "%d:%d:%d:%d", &cpp.sep_every, &cpp.sep_k, &cpp.pool_cap, &cpp.drop_after) < 2) { fprintf(stderr, "bad -cpsep\n"); return 1; }
    }
    else if (o == "-d" && more) degree = atoi(argv[++a]);
    else if (o == "-T" && more) threads = atoi(argv[++a]);
    else if (o == "-v") verbose++;
    else if (o[0] == '-') { fprintf(stderr, "bad option %s\n", o.c_str()); return 1; }
    else file = argv[a];
  }
  if (!file) { fprintf(stderr, "usage: see the comment at the top of swapdrv.C\n"); return 1; }
  setvbuf(stdout, NULL, _IOLBF, 0);
  dnn_set_blas_threads(threads);
  QAP *q = readQAPLIB((char *)file);
  if (!q) return 1;
  set_fathom_tolerance(q);
  symmetrize_QAP(q);
  if (!have_cutoff) cutoff = ub < HUGE_VAL ? ub - DELTA : HUGE_VAL;
  const int n = q->A->m;
  QAPAssignment node(n);
  if (!parse_fixes(fixes, node, n)) { fprintf(stderr, "bad fixes\n"); return 1; }
  QAP qr;
  reduceQAP(q, &qr, node);
  const int m = n - node.nfix, n1 = m * m + 1;
  const size_t nn = (size_t)n1 * n1;

  Ctx cx;
  cx.n = n; cx.m = m; cx.n1 = n1;
  cx.A.resize((size_t)n * n); cx.B.resize((size_t)n * n);
  for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) { cx.A[(size_t)i * n + k] = q->A->me[i][k]; cx.B[(size_t)i * n + k] = q->B->me[i][k]; }
  cx.p.assign(n, -1);
  std::vector<int> locl(n, -1), facl(n, -1);
  for (int i = 0; i < n; i++) {
    if (node.p[i] != NONE) { cx.p[i] = node.p[i]; cx.fixf.push_back(i); }
    else { facl[i] = (int)cx.ff.size(); cx.ff.push_back(i); }
    if (node.pinv[i] == NONE) { locl[i] = (int)cx.fl.size(); cx.fl.push_back(i); }
  }
  const int d = (int)cx.fixf.size();
  printf("instance %s n=%d node %s m=%d d=%d cutoff %.4f  exchange cuts: %d pairwise%s\n", file, n, fixes.c_str(), m, d, cutoff,
         d * m * m, use_fp ? " + fixed-pair rows" : "");
  if (use_rot) printf("rotations: %d cuts\n", d * (d - 1) * m * m);

  // ---- the candidate cuts
  std::vector<XCut> cuts;
  for (int i : cx.fixf)
    for (int jl = 0; jl < m; jl++)
      for (int bl = 0; bl < m; bl++) {
        XCut t; t.type = 0; t.i = i; t.a = cx.p[i]; t.j = cx.ff[jl]; t.b = cx.fl[bl];
        t.row = 1 + jl * m + bl; t.exk = jl; t.exc = bl; t.mu = 0; t.pooled = false;
        t.move2();
        double c0 = cx.cpair(t);
        for (int k : cx.fixf) if (k != i) c0 += cx.g(t, k, cx.p[k]);
        t.diag = c0;
        cuts.push_back(t);
      }
  if (use_fp)
    for (int s = 0; s < d; s++)
      for (int s2 = s + 1; s2 < d; s2++) {
        XCut t; t.i = cx.fixf[s]; t.a = cx.p[t.i]; t.j = cx.fixf[s2]; t.b = cx.p[t.j]; t.mu = 0; t.pooled = false;
        t.move2();
        double c0 = cx.cpair(t);
        for (int k : cx.fixf) if (k != t.i && k != t.j) c0 += cx.g(t, k, cx.p[k]);
        t.type = 1; t.row = 0; t.exk = t.exc = -1; t.diag = c0;
        cuts.push_back(t);
        for (int l = 0; l < m; l++)
          for (int dl = 0; dl < m; dl++) {
            XCut r = t; r.type = 2; r.row = 1 + l * m + dl; r.exk = l; r.exc = dl;
            r.diag = c0 + cx.g(t, cx.ff[l], cx.fl[dl]);
            cuts.push_back(r);
          }
      }
  if (use_rot)
    // rotations: fixed i (at a) -> a2, fixed i2 (at a2) -> b, free j (at b) -> a; both directions through the ordered pairs
    for (int i : cx.fixf)
      for (int i2 : cx.fixf) {
        if (i2 == i) continue;
        for (int jl = 0; jl < m; jl++)
          for (int bl = 0; bl < m; bl++) {
            XCut t; t.type = 3; t.i = i; t.a = cx.p[i]; t.j = cx.ff[jl]; t.b = cx.fl[bl];
            t.row = 1 + jl * m + bl; t.exk = jl; t.exc = bl; t.mu = 0; t.pooled = false;
            t.ns = 3;
            t.S[0] = i; t.from[0] = t.a; t.to[0] = cx.p[i2];
            t.S[1] = i2; t.from[1] = cx.p[i2]; t.to[1] = t.b;
            t.S[2] = t.j; t.from[2] = t.b; t.to[2] = t.a;
            double c0 = cx.cpair(t);
            for (int k : cx.fixf) if (k != i && k != i2) c0 += cx.g(t, k, cx.p[k]);
            t.diag = c0;
            cuts.push_back(t);
          }
      }
  std::vector<double> u1, v1, u2, v2;
  for (auto &t : cuts) {
    cx.coef(t, u1);
    double s = t.diag * t.diag;
    for (double gv : u1) s += gv * gv;
    t.norm = sqrt(s);
  }
  // ---- guard: every cut's linear form against brute-force exchange gains on random completions
  {
    std::mt19937 rng(4242);
    double worst = 0; long checked = 0;
    std::vector<double> Yp(nn);
    for (int trial = 0; trial < ncheck; trial++) {
      std::vector<int> perm(m);
      for (int k = 0; k < m; k++) perm[k] = k;
      std::shuffle(perm.begin(), perm.end(), rng);
      std::vector<int> full(cx.p);
      for (int k = 0; k < m; k++) full[cx.ff[k]] = cx.fl[perm[k]];
      const double c = cx.cost(full);
      std::vector<double> y(n1, 0.0); y[0] = 1;
      for (int k = 0; k < m; k++) y[1 + k * m + perm[k]] = 1;
      for (int r = 0; r < n1; r++) for (int s = 0; s < n1; s++) Yp[(size_t)r * n1 + s] = y[r] * y[s];
      std::map<std::vector<int>, double> gain;
      for (auto &t : cuts) {
        const double v = lhs(cx, t, Yp, u1, v1, u2, v2);
        double want = 0;
        if (t.row == 0 || y[t.row] == 1) {
          std::vector<int> key;
          for (int s = 0; s < t.ns; s++) { key.push_back(t.S[s]); key.push_back(t.to[s]); }
          auto it = gain.find(key);
          if (it == gain.end()) {
            std::vector<int> sw(full);
            for (int s = 0; s < t.ns; s++) sw[t.S[s]] = t.to[s];
            it = gain.insert({key, cx.cost(sw) - c}).first;
          }
          want = it->second;
        }
        worst = std::max(worst, fabs(v - want)); checked++;
      }
    }
    printf("guard: %ld cut evaluations on %d random completions against brute-force exchange gains: max |err| %.3g\n", checked, ncheck, worst);
    if (worst > 1e-6 * (1 + fabs(q->shift) + 1e6)) { fprintf(stderr, "cut coefficients are wrong\n"); return 1; }
  }

  DNNParams base;
  base.eig_mode = DNN_EIG_SUBSPACE;
  base.cheb_degree = degree;
  base.gamma = 1.0; base.gamma1 = 1.0;   // as dnncuts (PRSM); DNNParams' own defaults are 1.618, 0
  DNNWarm warm, save;
  XState xs;
  xs.cx = &cx; xs.cuts = &cuts; xs.Ladd.assign(nn, 0.0); xs.ek = ek; xs.theta_cur = theta; xs.delta = edelta;
  xs.verbose = verbose; xs.topk = topk; xs.ecap = ecap; xs.emu = emu;
  std::vector<double> &Ladd = xs.Ladd;
  int total_iters = 0;
  double total_cpu = 0;
  std::vector<int> free_orig, other_orig, fixed_at(n, -1);
  for (int a = 0; a < n; a++) {
    if ((cp.side == 0 ? node.p[a] : node.pinv[a]) == NONE) free_orig.push_back(a);
    if ((cp.side == 0 ? node.pinv[a] : node.p[a]) == NONE) other_orig.push_back(a);
    else fixed_at[a] = cp.side == 0 ? node.pinv[a] : node.p[a];
  }
  if (joint) {
    // ---- the library's Polyak rounds on the regional triangles, the exchange term stepped
    // between rounds through the hook (-ek 0: the triangles alone, the reference)
    DNNExtHook hook;
    hook.ctx = &xs; hook.update = xhook; hook.L = Ladd.data(); hook.abs_ = 0; hook.cnst = 0;
    cp.first_iters = first_iters; cp.rounds = jrounds; cp.iters = jiters; cp.freeze_iters = equiet;
    cp.verbose = 1;
    if (ek > 0) cp.ext_hook = &hook;
    DNNCutResult jr = dnn_cut_rounds(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, cp, &warm);
    printf("FINAL best_lb %.6f first_lb %.6f gain %.6f  rounds %d  iterations %d (after the first %d)  status %s\n", jr.best_lb, jr.first_lb,
           jr.best_lb - jr.first_lb, jr.rounds, jr.iterations, first_iters, dnn_status_name(jr.status));
    delete q;
    return 0;
  }
  // ---- eliminations: x_jb = 0 in every exchange-optimal completion when the largest gain
  // c0 + max over assignments of the rest of <g, x> of exchanging some fixed i with j at b is negative
  std::vector<unsigned char> elim((size_t)m * m, 0);
  int nelim = 0;
  if (use_elim && m >= 2) {
    std::vector<double> cst((size_t)(m - 1) * (m - 1));
    double closest = -HUGE_VAL;
    for (auto &t : cuts) {
      if (t.type != 0 && t.type != 3) continue;
      const int jl = t.exk, bl = t.exc;
      if (elim[(size_t)jl * m + bl]) continue;
      cx.coef(t, u1);
      int r = 0;
      for (int k = 0; k < m; k++) {
        if (k == jl) continue;
        int c2 = 0;
        for (int c = 0; c < m; c++) { if (c == bl) continue; cst[(size_t)r * (m - 1) + c2] = -u1[(size_t)k * m + c]; c2++; }
        r++;
      }
      const double mx = t.diag - lap_min(cst, m - 1);   // the largest exchange gain over the completions
      if (mx < -1e-6) { elim[(size_t)jl * m + bl] = 1; nelim++; }
      else closest = std::max(closest, -mx);
    }
    printf("eliminations: %d of %d free assignments can never be exchange-optimal against a fixed facility (LAP test)\n", nelim, m * m);
  }
  int active = 0;
  auto solve = [&](int iters, bool warm_start) {
    active = xs.assemble();
    DNNParams par = base;
    par.max_iter = iters;
    par.hopeless_after = 0; par.hopeless_gap_ratio = 0;
    par.L_add = active ? Ladd.data() : NULL;
    par.L_add_abs = xs.abs_add;
    par.cut_round = 1;
    if (nelim > 0 && warm_start) par.elim = elim.data();
    DNNResult res = dnn_bound(qr.A, qr.B, qr.C, qr.shift - xs.hsum, cutoff, par, warm_start ? &warm : NULL, &save);
    warm = save;
    total_iters += res.iterations; total_cpu += res.cpu;
    return res;
  };

  // ---- round 0: the plain solve
  DNNResult res = solve(first_iters, false);
  const double first_lb = res.lb;
  double best_lb = res.lb;
  DNNWarm best_save = warm;
  {
    // the primal value recomputed from Y (index convention check)
    const std::vector<double> &Y = warm.Y;
    double pv = qr.shift;
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
      const double *row = Y.data() + (size_t)(1 + i * m + j) * n1;
      pv += qr.C->me[i][j] * Y[1 + i * m + j];
      for (int k = 0; k < m; k++) for (int l = 0; l < m; l++) pv += qr.A->me[i][k] * qr.B->me[j][l] * row[1 + k * m + l];
    }
    printf("round  0  lb %.4f  primal %.4f (recomputed from Y %.4f)  status %s  iters %d  cpu %.1f\n", res.lb, res.primal, pv,
           dnn_status_name(res.status), res.iterations, res.cpu);
  }
  if (nelim > 0) {
    double mass = 0;
    for (int k = 0; k < m * m; k++) if (elim[k]) mass += warm.Y[(size_t)(1 + k) * n1 + 1 + k];
    printf("eliminations: probability mass %.4f of %d on the eliminated assignments at the plain iterate\n", mass, m);
  }
  bool fathomed = res.status == DNN_FATHOMED;
  for (int r = 0; r <= erounds && !fathomed; r++) {
    if (r > 0) {
      res = solve(eiters, true);
      if (res.lb > best_lb) { best_lb = res.lb; best_save = warm; }
      printf("round %2d  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f  active %d\n", r, res.lb,
             res.lb - first_lb, best_lb - first_lb, res.primal, dnn_status_name(res.status), res.iterations, res.cpu, active);
      if (res.status == DNN_FATHOMED) { fathomed = true; break; }
    }
    // separation and statistics on the last primal iterate, then the Polyak step
    xs.update(warm.Y, res.lb, best_lb, cutoff, r, r < erounds);
  }
  // ---- back to the best round's multipliers
  if (erounds > 0 && !fathomed) {
    xs.restore_best();
    warm = best_save;
    if (equiet > 0) {
      res = solve(equiet, true);
      best_lb = std::max(best_lb, res.lb);
      printf("quiet round  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f  active %d\n", res.lb,
             res.lb - first_lb, best_lb - first_lb, res.primal, dnn_status_name(res.status), res.iterations, res.cpu, active);
      fathomed = res.status == DNN_FATHOMED;
    }
  }
  printf("EXCHANGE best_lb %.6f first_lb %.6f gain %.6f  iterations %d  cpu %.1f%s\n", best_lb, first_lb, best_lb - first_lb,
         total_iters, total_cpu, fathomed ? "  FATHOMED" : "");
  active = xs.assemble();
  const double abs_add = xs.abs_add, hsum = xs.hsum;
  if (!saveL.empty()) {
    FILE *f = fopen(saveL.c_str(), "wb");
    if (f) { fwrite(Ladd.data(), sizeof(double), nn, f); fclose(f); printf("wrote the exchange term (%d active cuts, const %.6f) to %s\n", active, hsum, saveL.c_str()); }
  }
  // ---- the regional-triangle CP phase of dnncuts, the exchange term frozen as L_ext
  if (!fathomed && cpp.sep_every > 0 && cp_iters > 0 && (cp.grid_rows > 0 || cp.prefix)) {
    if (active) { cp.L_ext = Ladd.data(); cp.L_ext_abs = abs_add; cp.ext_const = hsum; }
    cp.verbose = 1;
    DNNCutParams c0 = cp;
    c0.rounds = 0; c0.first_iters = 0; c0.freeze_iters = 0; c0.freeze_after = 0;
    c0.pool_in = NULL; c0.pool_out = NULL;
    DNNPool p0, pin;
    c0.pool_all_round = 0; c0.pool_all_out = &p0; c0.other_orig = &other_orig;
    dnn_cut_rounds(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, c0, &warm, best_lb, res.primal);
    pin = p0;
    printf("cp pool: the separation at the start, %zu regions, %zu cuts; in-solve separation every %d iterations, %d cuts\n",
           pin.regions.size(), pin.cuts.size(), cpp.sep_every, cpp.sep_k);
    DNNParams p = base;
    p.max_iter = cp_iters;
    p.hopeless_after = 0; p.hopeless_gap_ratio = 0;
    p.warm_primal = 0;
    p.cut_round = 1;
    p.verbose = 1;
    DNNCutParams c1 = cp;
    c1.pool_in = &pin; c1.other_orig = &other_orig; c1.fixed_at = &fixed_at; c1.pool_out = NULL; c1.pool_all_out = NULL;
    DNNWarm sv;
    int kept = 0;
    DNNResult r = dnn_cp_bound(qr.A, qr.B, qr.C, qr.shift, cutoff, p, free_orig, n, c1, cpp, &warm, &sv, &kept);
    total_iters += r.iterations; total_cpu += r.cpu;
    printf("cp solve  lb %.10f (gain %.4f over the plain first solve)  primal %.4f  status %s  iters %d  cpu %.1f  kept %d\n",
           r.lb, r.lb - first_lb, r.primal, dnn_status_name(r.status), r.iterations, r.cpu, kept);
    best_lb = std::max(best_lb, r.lb);
    if (sv.Y.size() == nn) { printf("exchange cuts at the CP solve's final iterate:\n"); xs.verbose = 1; xs.update(sv.Y, r.lb, best_lb, cutoff, 0, false); }
  }
  printf("FINAL best_lb %.6f first_lb %.6f gain %.6f  iterations %d  cpu %.1f\n", best_lb, first_lb, best_lb - first_lb, total_iters, total_cpu);
  delete q;
  return 0;
}
