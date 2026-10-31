/*
 dnncuts: the certified DNN bound strengthened by triangle inequalities, in
 Lagrangian form (a prototype; see README.md and
 ../wil50_feasibility/stronger_bounds/README.md).

   make -C latest dnncuts
   experiments/dnn_cuts/dnncuts [-x i:j,...] [-u ub | -c cutoff] [-r rounds]
       [-I first-round iters] [-i iters] [-k cuts] [-t theta] [-delta D]
       [-mu0 M] [-grid R:C | -prefix] [-d degree] [-e eigmode] [-T threads] [-v]
       file.dat

 Two cut families, both valid for every permutation and folded into the
 objective with multipliers mu >= 0:

 Atomic (the default; experiments/dnn_cuts/README.md): for lifted variables
 a, b, c (location-facility pairs with distinct locations and distinct
 facilities), the Boolean-quadric triangle inequalities on Y = yy'
   (T1)  Y_ab + Y_ac - Y_bc - Y_aa <= 0      (apex a)
   (T2)  Y_aa + Y_bb + Y_cc - Y_ab - Y_ac - Y_bc <= 1

 Regional (-grid or -prefix; Astra's memo in ../wil50_feasibility/
 stronger_bounds/): for a set S of locations and facilities i, j, k, the
 membership variables z_i = sum_{a in S} x_{ai} are 0/1 with
   p_i = sum_{a in S} Y_(ai),(ai),   W_ij = sum_{a,b in S} Y_(ai),(bj)
 their first and second moments, and the same four triangle inequalities
 hold for (p, W):
   W_ij + W_ik - W_jk - p_i <= 0  (and the apexes j, k)
   p_i + p_j + p_k - W_ij - W_ik - W_jk <= 1.
 These are not implied by the DNN relaxation with the atomic triangles (an
 exact witness is in the memo). -grid R:C takes the regions of a rows x
 columns grid of locations, the row prefixes {a : a / C < t} and column
 prefixes {a : a mod C < t}, which decompose the Manhattan distance exactly
 (wil50: 5:10); the regions are in original location indices and are mapped
 to the free locations of the node. -prefix takes the prefixes of the free
 locations {0..t-1}, valid for any instance (for validation).

 With multipliers mu >= 0,
   <L, Y> >= <L + sum mu_t G_t, Y> - sum mu_t h_t
 for every permutation (each g_t(Y) - h_t <= 0), so the certified DNN bound of
 the modified objective, with the shift lowered by sum mu_t h_t, is a valid
 bound for the node: every value printed is certified. Each round solves the
 DNN relaxation (warm-started from the previous multiplier Z), separates the
 most violated cuts in its last primal iterate Y, and updates mu by a
 Polyak subgradient step toward target = min(cutoff, best + delta):
   mu_t <- max(0, mu_t + theta (target - lb) / |s|^2 * s_t),  s_t = g_t(Y) - h_t.
 New cuts enter with multiplier mu0 (default 0: the step sets them).
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
#include <set>
#include <algorithm>

// s < 0: atomic, a, b, c lifted indices; s >= 0: region s, a, b, c facilities
// (local indices), type 0/1/2 = apex a/b/c, 3 = the sum inequality
struct Cut { int type, a, b, c, s; double mu; };

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

// ---- atomic cuts

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

// ---- regional cuts

struct Region { std::string name; std::vector<int> loc; };   // local location indices

// p_i = sum_{a in S} Y_(ai),(ai),  W_ij = sum_{a,b in S} Y_(ai),(bj)
static void aggregate(const Region &S, const std::vector<double> &Y, int m, int n1,
                      std::vector<double> &p, std::vector<double> &W)
{
  p.assign(m, 0.0); W.assign((size_t)m * m, 0.0);
  for (int a : S.loc)
    for (int i = 0; i < m; i++) {
      const size_t r = 1 + (size_t)a * m + i;
      p[i] += Y[r * n1 + r];
      double *w = &W[(size_t)i * m];
      for (int b : S.loc) {
        const double *yr = &Y[r * n1 + 1 + (size_t)b * m];
        for (int j = 0; j < m; j++) w[j] += yr[j];
      }
    }
}

static double rviolation(const Cut &t, const std::vector<double> &p,
                         const std::vector<double> &W, int m)
{
  auto w = [&](int i, int j) { return W[(size_t)i * m + j]; };
  const int i = t.a, j = t.b, k = t.c;
  switch (t.type) {
    case 0: return w(i, j) + w(i, k) - w(j, k) - p[i];
    case 1: return w(i, j) + w(j, k) - w(i, k) - p[j];
    case 2: return w(i, k) + w(j, k) - w(i, j) - p[k];
    default: return p[i] + p[j] + p[k] - w(i, j) - w(i, k) - w(j, k) - 1;
  }
}

static double add_rcut(const Cut &t, const Region &S, std::vector<double> &L, int m, int n1,
                       double mu)
{
  double abs_ = 0;
  auto pair = [&](int i, int j, double v) {
    for (int a : S.loc)
      for (int b : S.loc) {
        if (a == b) continue;   // gangster zero
        const size_t r = 1 + (size_t)a * m + i, c = 1 + (size_t)b * m + j;
        L[r * n1 + c] += v / 2; L[c * n1 + r] += v / 2; abs_ += fabs(v);
      }
  };
  auto lin = [&](int i, double v) {
    for (int a : S.loc) {
      const size_t r = 1 + (size_t)a * m + i;
      L[r * n1 + r] += v; abs_ += fabs(v);
    }
  };
  const int i = t.a, j = t.b, k = t.c;
  if (t.type < 3) {
    const int apex = t.type == 0 ? i : (t.type == 1 ? j : k);
    const int o1 = t.type == 0 ? j : i, o2 = t.type == 2 ? j : k;
    pair(apex, o1, mu); pair(apex, o2, mu); pair(o1, o2, -mu); lin(apex, -mu);
  } else {
    lin(i, mu); lin(j, mu); lin(k, mu);
    pair(i, j, -mu); pair(i, k, -mu); pair(j, k, -mu);
  }
  return abs_;
}

// the regions of a rows x cols grid in original location indices, mapped to the
// node's free locations (local indices), without empty, full and duplicate or
// complementary sets
static std::vector<Region> grid_regions(int rows, int cols, const QAPAssignment &node, int n,
                                        int m)
{
  std::vector<int> local(n, -1);
  for (int a = 0, ii = 0; a < n; a++) if (node.p[a] == NONE) local[a] = ii++;
  std::vector<std::pair<std::string, std::vector<char>>> cand;
  for (int t = 1; t < rows; t++) {
    std::vector<char> in(n, 0);
    for (int a = 0; a < n; a++) in[a] = a / cols < t;
    cand.push_back({"row<" + std::to_string(t), in});
  }
  for (int t = 1; t < cols; t++) {
    std::vector<char> in(n, 0);
    for (int a = 0; a < n; a++) in[a] = a % cols < t;
    cand.push_back({"col<" + std::to_string(t), in});
  }
  std::vector<Region> out;
  std::set<std::vector<int>> seen;
  for (auto &c : cand) {
    std::vector<int> loc, comp;
    for (int a = 0; a < n; a++)
      if (local[a] >= 0) (c.second[a] ? loc : comp).push_back(local[a]);
    if (loc.empty() || comp.empty() || seen.count(loc) || seen.count(comp)) continue;
    seen.insert(loc);
    out.push_back({c.first, loc});
  }
  (void)m;
  return out;
}

static std::vector<Region> prefix_regions(int m)
{
  std::vector<Region> out;
  for (int t = 1; t < m; t++) {
    Region r; r.name = "loc<" + std::to_string(t);
    for (int a = 0; a < t; a++) r.loc.push_back(a);
    out.push_back(r);
  }
  return out;
}

int main(int argc, char **argv)
{
  std::string fixes = "-";
  double ub = HUGE_VAL, cutoff = HUGE_VAL, theta = 1.0, delta = 50.0, mu0 = 0.0;
  bool have_cutoff = false, prefix = false;
  int rounds = 15, iters = 400, first_iters = 0, kcuts = 2000, threads = 1, grid_r = 0,
      grid_c = 0, degree = 4, eig_mode = DNN_EIG_SUBSPACE, verbose = 0;
  double gamma = 1.0, gamma1 = 1.0;
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
    else if (o == "-delta" && more) delta = atof(argv[++a]);
    else if (o == "-mu0" && more) mu0 = atof(argv[++a]);
    else if (o == "-grid" && more) {
      if (sscanf(argv[++a], "%d:%d", &grid_r, &grid_c) != 2 || grid_r < 1 || grid_c < 1) {
        fprintf(stderr, "bad -grid R:C\n"); return 1; }
    }
    else if (o == "-prefix") prefix = true;
    else if (o == "-d" && more) degree = atoi(argv[++a]);
    else if (o == "-e" && more) eig_mode = atoi(argv[++a]);
    else if (o == "-g" && more) gamma = atof(argv[++a]);
    else if (o == "-G" && more) gamma1 = atof(argv[++a]);
    else if (o == "-T" && more) threads = atoi(argv[++a]);
    else if (o == "-v") verbose++;
    else if (o[0] == '-') { fprintf(stderr, "bad option %s\n", o.c_str()); return 1; }
    else file = argv[a];
  }
  if (!file) { fprintf(stderr, "usage: see the comment at the top of dnncuts.C\n"); return 1; }
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
  if (grid_r > 0 && grid_r * grid_c != n) { fprintf(stderr, "-grid R:C needs R*C = n\n"); return 1; }
  QAP qr;
  reduceQAP(q, &qr, node);
  const int m = n - node.nfix, n1 = m * m + 1;
  const bool regional = grid_r > 0 || prefix;
  std::vector<Region> regions;
  if (grid_r > 0) regions = grid_regions(grid_r, grid_c, node, n, m);
  else if (prefix) regions = prefix_regions(m);
  printf("instance %s n=%d node %s m=%d cutoff %.4f  cuts %s", file, n, fixes.c_str(), m, cutoff,
         regional ? (grid_r > 0 ? "regional (grid)" : "regional (prefixes)") : "atomic");
  if (regional) printf(", %zu regions", regions.size());
  printf("\n");
  if (verbose && regional)
    for (auto &R : regions) printf("  region %-8s capacity %zu\n", R.name.c_str(), R.loc.size());

  std::vector<Cut> cuts;
  std::map<std::vector<int>, int> index;
  std::vector<double> Ladd((size_t)n1 * n1, 0.0);
  std::vector<double> p, W;                       // regional aggregates of the last Y
  std::vector<std::vector<double>> P, WW;         // per region
  DNNWarm warm, save;
  double first_lb = -HUGE_VAL, best_lb = -HUGE_VAL, prev_lb = -HUGE_VAL, theta_cur = theta;

  for (int r = 0; r <= rounds; r++) {
    // the modified objective for the current multipliers
    std::fill(Ladd.begin(), Ladd.end(), 0.0);
    double abs_add = 0, hsum = 0;
    int active = 0;
    for (const Cut &t : cuts)
      if (t.mu > 0) {
        abs_add += t.s < 0 ? add_cut(t, Ladd, n1, t.mu) : add_rcut(t, regions[t.s], Ladd, m, n1, t.mu);
        if (t.type == (t.s < 0 ? 2 : 3)) hsum += t.mu;
        active++;
      }
    DNNParams par;
    par.max_iter = (r == 0 && first_iters > 0) ? first_iters : iters;
    par.hopeless_after = 0;
    par.eig_mode = eig_mode;
    par.cheb_degree = degree;
    par.gamma = gamma; par.gamma1 = gamma1;
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
    size_t nviol = 0;
    auto bound_memory = [&]() {
      if (found.size() > 4000000) {   // keep memory bounded: the most violated
        size_t keep = std::max<size_t>(100000, 8 * (size_t)kcuts);
        std::nth_element(found.begin(), found.begin() + keep, found.end(),
                         [](const std::pair<double, Cut> &u, const std::pair<double, Cut> &w) {
                           return u.first > w.first; });
        found.resize(keep);
      }
    };
    if (regional) {
      P.assign(regions.size(), {}); WW.assign(regions.size(), {});
      for (size_t s = 0; s < regions.size(); s++) {
        aggregate(regions[s], Y, m, n1, P[s], WW[s]);
        double rmax = -HUGE_VAL; size_t rn = 0;
        for (int i = 0; i < m; i++)
          for (int j = i + 1; j < m; j++)
            for (int k = j + 1; k < m; k++)
              for (int type = 0; type < 4; type++) {
                Cut t = {type, i, j, k, (int)s, 0.0};
                double v = rviolation(t, P[s], WW[s], m);
                rmax = std::max(rmax, v);
                if (v > 1e-4) { found.push_back({v, t}); rn++; }
              }
        nviol += rn;
        if (verbose > 1)
          printf("  region %-8s violated %zu  max %.4f\n", regions[s].name.c_str(), rn, rmax);
        bound_memory();
      }
    } else {
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
            Cut cand[4] = {{1, a, b, c, -1, 0}, {1, b, a, c, -1, 0}, {1, c, a, b, -1, 0},
                           {2, a, b, c, -1, 0}};
            for (const Cut &t : cand) {
              double v = violation(t, Y, n1);
              if (v > 1e-4) { found.push_back({v, t}); nviol++; }
            }
            bound_memory();
          }
        }
    }
    std::sort(found.begin(), found.end(),
              [](const std::pair<double, Cut> &u, const std::pair<double, Cut> &w) {
                return u.first > w.first; });
    int added = 0;
    for (auto &f : found) {
      if (added >= kcuts) break;
      std::vector<int> key = {f.second.type, f.second.a, f.second.b, f.second.c, f.second.s};
      if (index.count(key)) continue;
      index[key] = (int)cuts.size();
      f.second.mu = mu0;
      cuts.push_back(f.second);
      added++;
    }
    double maxv = found.empty() ? 0 : found[0].first;
    printf("round %2d  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f"
           "  active %d  violated %zu (max %.4f)  added %d\n", r, res.lb, res.lb - first_lb,
           best_lb - first_lb, res.primal, dnn_status_name(res.status), res.iterations, res.cpu,
           active, nviol, maxv, added);
    if (res.status == DNN_FATHOMED || r == rounds) break;
    // Polyak subgradient step on all cuts toward the target; the step factor
    // halves whenever a round lowers the bound (an overshoot)
    if (theta <= 0) continue;
    if (r > 0 && res.lb < prev_lb - 1e-9) theta_cur *= 0.5;
    prev_lb = res.lb;
    double s2 = 0;
    std::vector<double> s(cuts.size());
    for (size_t t = 0; t < cuts.size(); t++) {
      s[t] = cuts[t].s < 0 ? violation(cuts[t], Y, n1)
                           : rviolation(cuts[t], P[cuts[t].s], WW[cuts[t].s], m);
      if (cuts[t].mu > 0 || s[t] > 0) s2 += s[t] * s[t];
    }
    if (s2 <= 0) break;
    const double target = std::min(cutoff, best_lb + delta);
    double alpha = theta_cur * std::max(target - res.lb, 1e-6) / s2;
    for (size_t t = 0; t < cuts.size(); t++)
      cuts[t].mu = std::max(0.0, cuts[t].mu + alpha * s[t]);
  }
  printf("FINAL best_lb %.6f first_lb %.6f gain %.6f\n", best_lb, first_lb, best_lb - first_lb);
  delete q;
  return 0;
}
