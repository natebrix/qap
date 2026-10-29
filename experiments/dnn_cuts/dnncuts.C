/*
 dnncuts: the certified DNN bound strengthened by triangle inequalities, in
 Lagrangian form (a prototype; see README.md).

   make -C latest dnncuts
   experiments/dnn_cuts/dnncuts [-x i:j,...] [-u ub | -c cutoff] [-r rounds]
                                [-I first-round iters] [-i iters] [-k cuts] [-t theta] [-T threads] file.dat

 For lifted variables a, b, c (facility-location pairs with distinct
 facilities and distinct locations), every lifted permutation Y = yy' satisfies
   (T1)  Y_ab + Y_ac - Y_bc - Y_aa <= 0      (apex a)
   (T2)  Y_aa + Y_bb + Y_cc - Y_ab - Y_ac - Y_bc <= 1
 but the DNN relaxation does not imply them. With multipliers mu >= 0,
   <L, Y> >= <L + sum mu_t G_t, Y> - sum mu_t h_t
 for every permutation (each g_t(Y) - h_t <= 0), so the certified DNN bound of
 the modified objective, with the shift lowered by sum mu_t h_t, is a valid
 bound for the node: every value printed is certified. Each round solves the
 DNN relaxation (warm-started from the previous multiplier Z), separates the
 most violated triangles in its last primal iterate Y, and updates mu by a
 Polyak subgradient step toward the cutoff:
   mu_t <- max(0, mu_t + theta (cutoff - lb) / |s|^2 * s_t),  s_t = g_t(Y) - h_t.
*/

#include "dnn.h"
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
#include <algorithm>

struct Cut { int type, a, b, c; double mu; };

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

// g_t(Y) - h_t
static double violation(const Cut &t, const std::vector<double> &Y, int n1)
{
  auto y = [&](int r, int c) { return Y[(size_t)r * n1 + c]; };
  if (t.type == 1)
    return y(t.a, t.b) + y(t.a, t.c) - y(t.b, t.c) - y(t.a, t.a);
  return y(t.a, t.a) + y(t.b, t.b) + y(t.c, t.c) - y(t.a, t.b) - y(t.a, t.c) - y(t.b, t.c) - 1;
}

// add mu * G_t to L (symmetric: an off-diagonal coefficient c on Y_ab is c/2 on
// both (a,b) and (b,a)); returns sum |entries added|
static double add_cut(const Cut &t, std::vector<double> &L, int n1, double mu)
{
  auto off = [&](int r, int c, double v) {
    L[(size_t)r * n1 + c] += v / 2; L[(size_t)c * n1 + r] += v / 2; };
  if (t.type == 1) {
    off(t.a, t.b, mu); off(t.a, t.c, mu); off(t.b, t.c, -mu);
    L[(size_t)t.a * n1 + t.a] -= mu;
  } else {
    L[(size_t)t.a * n1 + t.a] += mu; L[(size_t)t.b * n1 + t.b] += mu;
    L[(size_t)t.c * n1 + t.c] += mu;
    off(t.a, t.b, -mu); off(t.a, t.c, -mu); off(t.b, t.c, -mu);
  }
  return 6 * mu;
}

int main(int argc, char **argv)
{
  std::string fixes = "-";
  double ub = HUGE_VAL, cutoff = HUGE_VAL, theta = 1.0;
  bool have_cutoff = false;
  int rounds = 15, iters = 400, first_iters = 0, kcuts = 2000, threads = 1;
  const char *file = NULL;
  for (int a = 1; a < argc; a++) {
    std::string o = argv[a];
    bool more = a + 1 < argc;
    if (o == "-x" && more) fixes = argv[++a];
    else if (o == "-u" && more) ub = atof(argv[++a]);
    else if (o == "-c" && more) { cutoff = atof(argv[++a]); have_cutoff = true; }
    else if (o == "-r" && more) rounds = atoi(argv[++a]);
    else if (o == "-i" && more) iters = atoi(argv[++a]);
    else if (o == "-I" && more) first_iters = atoi(argv[++a]);
    else if (o == "-k" && more) kcuts = atoi(argv[++a]);
    else if (o == "-t" && more) theta = atof(argv[++a]);
    else if (o == "-T" && more) threads = atoi(argv[++a]);
    else if (o[0] == '-') { fprintf(stderr, "bad option %s\n", o.c_str()); return 1; }
    else file = argv[a];
  }
  if (!file) { fprintf(stderr, "usage: see the comment at the top of dnncuts.C\n"); return 1; }
  dnn_set_blas_threads(threads);
  QAP *q = readQAPLIB((char *)file);
  if (!q) return 1;
  set_fathom_tolerance(q);
  symmetrize_QAP(q);
  if (!have_cutoff) cutoff = ub < HUGE_VAL ? ub - DELTA : HUGE_VAL;
  if (!(cutoff < HUGE_VAL)) { fprintf(stderr, "need -u or -c\n"); return 1; }
  const int n = q->A->m;
  QAPAssignment node(n);
  if (!parse_fixes(fixes, node, n)) { fprintf(stderr, "bad fixes\n"); return 1; }
  QAP qr;
  reduceQAP(q, &qr, node);
  const int m = n - node.nfix, n1 = m * m + 1;
  printf("instance %s n=%d node %s m=%d cutoff %.4f\n", file, n, fixes.c_str(), m, cutoff);

  std::vector<Cut> cuts;
  std::map<long long, int> index;
  auto key = [&](int type, int a, int b, int c) {
    return (((long long)type * n1 + a) * n1 + b) * (long long)n1 + c; };
  std::vector<double> Ladd((size_t)n1 * n1, 0.0);
  DNNWarm warm, save;
  double first_lb = -HUGE_VAL, best_lb = -HUGE_VAL;

  for (int r = 0; r <= rounds; r++) {
    // the modified objective for the current multipliers
    std::fill(Ladd.begin(), Ladd.end(), 0.0);
    double abs_add = 0, hsum = 0;
    int active = 0;
    for (const Cut &t : cuts)
      if (t.mu > 0) {
        abs_add += add_cut(t, Ladd, n1, t.mu);
        if (t.type == 2) hsum += t.mu;
        active++;
      }
    DNNParams par;
    par.max_iter = (r == 0 && first_iters > 0) ? first_iters : iters;
    par.hopeless_after = 0;
    par.eig_mode = DNN_EIG_SINGLE;
    par.L_add = active ? Ladd.data() : NULL;
    par.L_add_abs = abs_add;
    DNNResult res = dnn_bound(qr.A, qr.B, qr.C, qr.shift - hsum, cutoff, par,
                              r > 0 ? &warm : NULL, &save);
    warm = save;
    if (r == 0) first_lb = res.lb;
    best_lb = std::max(best_lb, res.lb);
    // separate on the last primal iterate
    const std::vector<double> &Y = save.Y;
    std::vector<std::pair<double, Cut>> found;
    std::vector<int> sup;
    for (int a = 1; a < n1; a++) if (Y[(size_t)a * n1 + a] > 1e-4) sup.push_back(a);
    auto compatible = [&](int a, int b) {
      int ia = (a - 1) / m, ja = (a - 1) % m, ib = (b - 1) / m, jb = (b - 1) % m;
      return ia != ib && ja != jb; };
    for (size_t x = 0; x < sup.size(); x++)
      for (size_t yv = x + 1; yv < sup.size(); yv++) {
        int a = sup[x], b = sup[yv];
        if (!compatible(a, b)) continue;
        for (size_t z = yv + 1; z < sup.size(); z++) {
          int c = sup[z];
          if (!compatible(a, c) || !compatible(b, c)) continue;
          Cut cand[4] = {{1, a, b, c, 0}, {1, b, a, c, 0}, {1, c, a, b, 0}, {2, a, b, c, 0}};
          for (const Cut &t : cand) {
            double v = violation(t, Y, n1);
            if (v > 1e-4) found.push_back({v, t});
          }
          if (found.size() > 4000000) {   // keep memory bounded: the most violated
            size_t keep = std::max<size_t>(100000, 8 * (size_t)kcuts);
            std::nth_element(found.begin(), found.begin() + keep, found.end(),
                             [](const std::pair<double, Cut> &u, const std::pair<double, Cut> &w) {
                               return u.first > w.first; });
            found.resize(keep);
          }
        }
      }
    std::sort(found.begin(), found.end(),
              [](const std::pair<double, Cut> &u, const std::pair<double, Cut> &w) {
                return u.first > w.first; });
    int added = 0;
    for (auto &f : found) {
      if (added >= kcuts) break;
      long long k = key(f.second.type, f.second.a, f.second.b, f.second.c);
      if (index.count(k)) continue;
      index[k] = (int)cuts.size();
      cuts.push_back(f.second);
      added++;
    }
    double maxv = found.empty() ? 0 : found[0].first;
    printf("round %2d  lb %.4f (gain %.4f, best %.4f)  status %s  iters %d  cpu %.1f  active %d  "
           "violated %zu (max %.4f)  added %d\n", r, res.lb, res.lb - first_lb, best_lb - first_lb,
           dnn_status_name(res.status), res.iterations, res.cpu, active, found.size(), maxv,
           added);
    fflush(stdout);
    if (res.status == DNN_FATHOMED || r == rounds) break;
    // Polyak subgradient step on all cuts toward the cutoff
    double s2 = 0;
    std::vector<double> s(cuts.size());
    for (size_t t = 0; t < cuts.size(); t++) {
      s[t] = violation(cuts[t], Y, n1);
      if (cuts[t].mu > 0 || s[t] > 0) s2 += s[t] * s[t];
    }
    if (s2 <= 0) break;
    double alpha = theta * std::max(cutoff - res.lb, 1e-6) / s2;
    for (size_t t = 0; t < cuts.size(); t++)
      cuts[t].mu = std::max(0.0, cuts[t].mu + alpha * s[t]);
  }
  delete q;
  return 0;
}
