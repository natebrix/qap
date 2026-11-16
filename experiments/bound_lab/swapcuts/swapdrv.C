/*
 swapdrv: exchange-optimality cuts on top of the certified DNN bound (README.md here).
 NOT valid for every permutation: valid for every permutation that no pairwise exchange
 improves, hence for every global optimum, hence usable in a search that proves "no
 permutation costs less than the incumbent".

   swapdrv [-x i:a,...] [-u ub | -c cutoff] [-I first iters] [-er rounds] [-ei iters]
           [-ek cuts per round] [-et theta] [-edelta D] [-eq quiet iters] [-fp]
           [-grid R:C -side S -k K -cpsep k:K -cpiters N]  [-saveL file] [-check N]
           [-T threads] [-v] file.dat

 Cuts. Cost C(p) = sum A[k][l] B[p(k)][p(l)]. For a facility i fixed at a, a free facility
 j and a free location b, with x_jb = 1 the gain of exchanging i and j is
   delta = c0 + sum_{free k != j, free c != b} g_k(c) x_kc,
   g_k(c) = (A[i][k]-A[j][k]) (B[b][c]-B[a][c]) + (A[k][i]-A[k][j]) (B[c][b]-B[c][a]),
   c0 = (A[i][j]-A[j][i]) (B[b][a]-B[a][b]) + (A[i][i]-A[j][j]) (B[b][b]-B[a][a])
        + sum over the other fixed k of g_k(p(k)),
 so  c0 Y[jb,jb] + sum g_k(c) Y[jb,kc] >= 0   (type P, one per fixed i, free j, free b).
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
  int type;      // 0: P (fixed i, free j at b), 1: F0 (two fixed, row 0), 2: F1 (two fixed, row ld)
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
  double g(const XCut &t, int k, int c) const {   // original k, c
    return (a_(t.i, k) - a_(t.j, k)) * (b_(t.b, c) - b_(t.a, c)) + (a_(k, t.i) - a_(k, t.j)) * (b_(c, t.b) - b_(c, t.a));
  }
  double cpair(const XCut &t) const {
    return (a_(t.i, t.j) - a_(t.j, t.i)) * (b_(t.b, t.a) - b_(t.a, t.b)) + (a_(t.i, t.i) - a_(t.j, t.j)) * (b_(t.b, t.b) - b_(t.a, t.a));
  }
  // the four factor vectors over local indices (zero at the excluded ones)
  void factors(const XCut &t, std::vector<double> &u1, std::vector<double> &v1, std::vector<double> &u2, std::vector<double> &v2) const {
    u1.assign(m, 0); v1.assign(m, 0); u2.assign(m, 0); v2.assign(m, 0);
    for (int k = 0; k < m; k++) {
      if (k == t.exk) continue;
      u1[k] = a_(t.i, ff[k]) - a_(t.j, ff[k]); u2[k] = a_(ff[k], t.i) - a_(ff[k], t.j);
    }
    for (int c = 0; c < m; c++) {
      if (c == t.exc) continue;
      v1[c] = b_(t.b, fl[c]) - b_(t.a, fl[c]); v2[c] = b_(fl[c], t.b) - b_(fl[c], t.a);
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
  cx.factors(t, u1, v1, u2, v2);
  const double *row = Y.data() + (size_t)t.row * n1;
  double s = t.diag * (t.row == 0 ? 1.0 : row[t.row]);
  for (int k = 0; k < m; k++) {
    if (u1[k] == 0 && u2[k] == 0) continue;
    const double *r = row + 1 + (size_t)k * m;
    double s1 = 0, s2 = 0;
    for (int c = 0; c < m; c++) { s1 += v1[c] * r[c]; s2 += v2[c] * r[c]; }
    s += u1[k] * s1 + u2[k] * s2;
  }
  return s;
}

int main(int argc, char **argv)
{
  std::string fixes = "-", saveL;
  double ub = HUGE_VAL, cutoff = HUGE_VAL;
  bool have_cutoff = false, use_fp = false;
  int first_iters = 300, erounds = 0, eiters = 200, ek = 500, equiet = 0, threads = 1, verbose = 0, ncheck = 50;
  int degree = 4, topk = 12;
  double theta = 1.0, edelta = 50.0;
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
    else if (o == "-fp") use_fp = true;
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

  // ---- the candidate cuts
  std::vector<XCut> cuts;
  for (int i : cx.fixf)
    for (int jl = 0; jl < m; jl++)
      for (int bl = 0; bl < m; bl++) {
        XCut t; t.type = 0; t.i = i; t.a = cx.p[i]; t.j = cx.ff[jl]; t.b = cx.fl[bl];
        t.row = 1 + jl * m + bl; t.exk = jl; t.exc = bl; t.mu = 0; t.pooled = false;
        double c0 = cx.cpair(t);
        for (int k : cx.fixf) if (k != i) c0 += cx.g(t, k, cx.p[k]);
        t.diag = c0;
        cuts.push_back(t);
      }
  if (use_fp)
    for (int s = 0; s < d; s++)
      for (int s2 = s + 1; s2 < d; s2++) {
        XCut t; t.i = cx.fixf[s]; t.a = cx.p[t.i]; t.j = cx.fixf[s2]; t.b = cx.p[t.j]; t.mu = 0; t.pooled = false;
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
  std::vector<double> u1, v1, u2, v2;
  for (auto &t : cuts) {
    cx.factors(t, u1, v1, u2, v2);
    double s = t.diag * t.diag;
    for (int k = 0; k < m; k++) for (int c = 0; c < m; c++) { double gv = u1[k] * v1[c] + u2[k] * v2[c]; s += gv * gv; }
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
      std::map<long, double> gain;
      for (auto &t : cuts) {
        const double v = lhs(cx, t, Yp, u1, v1, u2, v2);
        double want = 0;
        if (t.row == 0 || y[t.row] == 1) {
          const long key = (long)t.i * n + t.j;
          auto it = gain.find(key);
          if (it == gain.end()) {
            std::vector<int> sw(full); std::swap(sw[t.i], sw[t.j]);
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
  DNNWarm warm, save;
  std::vector<double> Ladd(nn, 0.0);
  double abs_add = 0, hsum = 0;
  int total_iters = 0;
  double total_cpu = 0;
  auto assemble_x = [&]() {
    std::fill(Ladd.begin(), Ladd.end(), 0.0);
    abs_add = 0; hsum = 0;
    int active = 0;
    for (auto &t : cuts) {
      if (!(t.mu > 0)) continue;
      active++;
      cx.factors(t, u1, v1, u2, v2);
      if (t.row == 0) hsum += t.mu * t.diag; else Ladd[(size_t)t.row * n1 + t.row] -= t.mu * t.diag;
      abs_add += t.mu * fabs(t.diag);
      for (int k = 0; k < m; k++) {
        if (u1[k] == 0 && u2[k] == 0) continue;
        for (int c = 0; c < m; c++) {
          const double gv = u1[k] * v1[c] + u2[k] * v2[c];
          if (gv == 0) continue;
          const int col = 1 + k * m + c;
          Ladd[(size_t)t.row * n1 + col] -= 0.5 * t.mu * gv;
          Ladd[(size_t)col * n1 + t.row] -= 0.5 * t.mu * gv;
          abs_add += t.mu * fabs(gv);
        }
      }
    }
    return active;
  };
  auto solve = [&](int iters, bool warm_start, int &active) {
    active = assemble_x();
    DNNParams par = base;
    par.max_iter = iters;
    par.hopeless_after = 0; par.hopeless_gap_ratio = 0;
    par.L_add = active ? Ladd.data() : NULL;
    par.L_add_abs = abs_add;
    par.cut_round = warm_start ? 1 : 0;
    DNNResult res = dnn_bound(qr.A, qr.B, qr.C, qr.shift - hsum, cutoff, par, warm_start ? &warm : NULL, &save);
    warm = save;
    total_iters += res.iterations; total_cpu += res.cpu;
    return res;
  };

  // ---- round 0: the plain solve
  int active = 0;
  DNNResult res = solve(first_iters, false, active);
  const double first_lb = res.lb;
  double best_lb = res.lb, prev_lb = -HUGE_VAL, theta_cur = theta;
  std::vector<double> best_mu(cuts.size(), 0.0);
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
  bool fathomed = res.status == DNN_FATHOMED;
  for (int r = 0; r <= erounds && !fathomed; r++) {
    if (r > 0) {
      res = solve(eiters, true, active);
      if (res.lb > best_lb) { best_lb = res.lb; for (size_t t = 0; t < cuts.size(); t++) best_mu[t] = cuts[t].mu; best_save = warm; }
      printf("round %2d  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f  active %d\n", r, res.lb,
             res.lb - first_lb, best_lb - first_lb, res.primal, dnn_status_name(res.status), res.iterations, res.cpu, active);
      if (res.status == DNN_FATHOMED) { fathomed = true; break; }
    }
    // ---- separation and statistics on the last primal iterate
    const std::vector<double> &Y = warm.Y;
    std::vector<double> sg(cuts.size());
    std::vector<std::pair<double, int>> viol;
    std::map<int, std::pair<int, double>> byfix;   // fixed facility -> (count, sum of -lhs)
    size_t nv = 0, nv1 = 0; double mass = 0, worst_cond = 0, worst_norm = 0;
    for (size_t t = 0; t < cuts.size(); t++) {
      const XCut &c = cuts[t];
      const double v = lhs(cx, c, Y, u1, v1, u2, v2);
      sg[t] = -v;
      if (v < -1e-6) {
        const double yd = c.row == 0 ? 1.0 : Y[(size_t)c.row * n1 + c.row];
        nv++; mass += -v;
        if (c.norm > 0) worst_norm = std::max(worst_norm, -v / c.norm);
        if (yd > 1e-4) { worst_cond = std::min(worst_cond, v / yd); if (v / yd < -1.0) nv1++; }
        auto &e = byfix[c.type == 0 ? c.i : -1 - c.type];
        e.first++; e.second += -v;
        viol.push_back({c.norm > 0 ? -v / c.norm : 0, (int)t});
      }
    }
    std::sort(viol.begin(), viol.end(), [](const std::pair<double, int> &x, const std::pair<double, int> &y) { return x.first > y.first; });
    printf("  separation %d: %zu of %zu cuts violated (%zu with conditional gain below -1), sum of violations %.2f (cost units x probability), "
           "worst conditional gain %.1f, largest normalised violation %.5f\n", r, nv, cuts.size(), nv1, mass, worst_cond, worst_norm);
    if (verbose || r == 0) {
      printf("  by fixed facility:");
      for (auto &e : byfix) {
        if (e.first >= 0) printf("  %d@%d: %d (%.2f)", e.first, cx.p[e.first], e.second.first, e.second.second);
        else printf("  fixed-pair type %d: %d (%.2f)", -1 - e.first, e.second.first, e.second.second);
      }
      printf("\n");
      for (int s = 0; s < std::min<int>(topk, (int)viol.size()); s++) {
        const XCut &c = cuts[viol[s].second];
        const double yd = c.row == 0 ? 1.0 : Y[(size_t)c.row * n1 + c.row];
        printf("    top %2d  type %d  fixed %d@%d  with %d@%d  row %d  y %.5f  lhs %.4f  conditional gain %.2f  normalised %.5f\n", s, c.type,
               c.i, c.a, c.j, c.b, c.row, yd, -sg[viol[s].second], yd > 0 ? -sg[viol[s].second] / yd : 0.0, viol[s].first);
      }
    }
    if (r == erounds) break;
    int added = 0;
    for (auto &e : viol) { if (added >= ek) break; if (!cuts[e.second].pooled) { cuts[e.second].pooled = true; added++; } }
    // ---- Polyak step on the pool
    const double lb_now = res.lb;
    if (r > 0 && lb_now < prev_lb - 1e-9) theta_cur *= 0.5;
    prev_lb = lb_now;
    double s2 = 0; size_t pool = 0;
    for (size_t t = 0; t < cuts.size(); t++)
      if (cuts[t].pooled) { pool++; if (cuts[t].mu > 0 || sg[t] > 0) s2 += sg[t] * sg[t]; }
    if (s2 <= 0) { printf("  no violated cut in the pool: the exchange phase ends\n"); break; }
    const double target = std::min(cutoff, best_lb + edelta);
    const double alpha = theta_cur * std::max(target - lb_now, 1e-6) / s2;
    for (size_t t = 0; t < cuts.size(); t++)
      if (cuts[t].pooled) cuts[t].mu = std::max(0.0, cuts[t].mu + alpha * sg[t]);
    printf("  added %d, pool %zu, step %.3e (theta %.3f, target %.1f)\n", added, pool, alpha, theta_cur, target);
  }
  // ---- back to the best round's multipliers
  if (erounds > 0 && !fathomed) {
    for (size_t t = 0; t < cuts.size(); t++) cuts[t].mu = best_mu[t];
    warm = best_save;
    if (equiet > 0) {
      res = solve(equiet, true, active);
      best_lb = std::max(best_lb, res.lb);
      printf("quiet round  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f  active %d\n", res.lb,
             res.lb - first_lb, best_lb - first_lb, res.primal, dnn_status_name(res.status), res.iterations, res.cpu, active);
      fathomed = res.status == DNN_FATHOMED;
    }
  }
  printf("EXCHANGE best_lb %.6f first_lb %.6f gain %.6f  iterations %d  cpu %.1f%s\n", best_lb, first_lb, best_lb - first_lb,
         total_iters, total_cpu, fathomed ? "  FATHOMED" : "");
  active = assemble_x();
  if (!saveL.empty()) {
    FILE *f = fopen(saveL.c_str(), "wb");
    if (f) { fwrite(Ladd.data(), sizeof(double), nn, f); fclose(f); printf("wrote the exchange term (%d active cuts, const %.6f) to %s\n", active, hsum, saveL.c_str()); }
  }
  // ---- the regional-triangle CP phase of dnncuts, the exchange term frozen as L_ext
  if (!fathomed && cpp.sep_every > 0 && cp_iters > 0 && (cp.grid_rows > 0 || cp.prefix)) {
    std::vector<int> free_orig, other_orig, fixed_at(n, -1);
    for (int a = 0; a < n; a++) {
      if ((cp.side == 0 ? node.p[a] : node.pinv[a]) == NONE) free_orig.push_back(a);
      if ((cp.side == 0 ? node.pinv[a] : node.p[a]) == NONE) other_orig.push_back(a);
      else fixed_at[a] = cp.side == 0 ? node.pinv[a] : node.p[a];
    }
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
  }
  printf("FINAL best_lb %.6f first_lb %.6f gain %.6f  iterations %d  cpu %.1f\n", best_lb, first_lb, best_lb - first_lb, total_iters, total_cpu);
  delete q;
  return 0;
}
