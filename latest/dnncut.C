// dnncut.C: triangle cuts on top of the DNN bound (see dnncut.h); the code of
// experiments/dnn_cuts/dnncuts.C made a library, with regions on either side.
#include "dnncut.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <map>
#include <set>
#include <algorithm>

DNNCutParams::DNNCutParams()
{
  rounds = 12; iters = 300; first_iters = 0; kcuts = 1000;
  delta = 50; theta = 1.0; mu0 = 0;
  grid_rows = grid_cols = 0; side = 0; prefix = false; atomic = false; verbose = 0;
  stall_rounds = 0; stall_factor = 3.0;
  cap_at_cutoff = 1; freeze_after = 0; freeze_iters = 0; freeze_if_primal_above = 0;
  freeze_best = 1;
}

namespace {

// s < 0: atomic, a, b, c lifted indices; s >= 0: region s, a, b, c free
// elements of the other side (local indices), type 0/1/2 = apex a/b/c, 3 =
// the sum inequality
struct Cut { int type, a, b, c, s; double mu; };
struct Region { std::string name; std::vector<int> loc; };   // local indices on `side`

// lifted index of the variable with local index s on the region side and t on
// the other side
struct Lift {
  int m, side;
  size_t operator()(int s, int t) const { return side == 0 ? 1 + (size_t)s * m + t : 1 + (size_t)t * m + s; }
};

// ---- atomic cuts

double violation(const Cut &t, const std::vector<double> &Y, int n1)
{
  auto y = [&](int r, int c) { return Y[(size_t)r * n1 + c]; };
  if (t.type == 1)
    return y(t.a, t.b) + y(t.a, t.c) - y(t.b, t.c) - y(t.a, t.a);
  return y(t.a, t.a) + y(t.b, t.b) + y(t.c, t.c) - y(t.a, t.b) - y(t.a, t.c) - y(t.b, t.c) - 1;
}

double add_cut(const Cut &t, std::vector<double> &L, int n1, double mu)
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

// p_i = sum_{a in S} Y_(ai),(ai),  W_ij = sum_{a,b in S} Y_(ai),(bj)
void aggregate(const Region &S, const std::vector<double> &Y, const Lift &lift, int n1,
               std::vector<double> &p, std::vector<double> &W)
{
  const int m = lift.m;
  p.assign(m, 0.0); W.assign((size_t)m * m, 0.0);
  for (int a : S.loc)
    for (int i = 0; i < m; i++) {
      const size_t r = lift(a, i);
      p[i] += Y[r * n1 + r];
      double *w = &W[(size_t)i * m];
      for (int b : S.loc)
        for (int j = 0; j < m; j++) w[j] += Y[r * n1 + lift(b, j)];
    }
}

double rviolation(const Cut &t, const std::vector<double> &p, const std::vector<double> &W, int m)
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

double add_rcut(const Cut &t, const Region &S, std::vector<double> &L, const Lift &lift, int n1,
                double mu)
{
  double abs_ = 0;
  auto pair = [&](int i, int j, double v) {
    for (int a : S.loc)
      for (int b : S.loc) {
        if (a == b) continue;   // gangster zero
        const size_t r = lift(a, i), c = lift(b, j);
        L[r * n1 + c] += v / 2; L[c * n1 + r] += v / 2; abs_ += fabs(v);
      }
  };
  auto lin = [&](int i, double v) {
    for (int a : S.loc) {
      const size_t r = lift(a, i);
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

// the regions of a rows x cols grid in original indices, mapped to the free
// elements (local indices), without empty, full, duplicate or complementary sets
std::vector<Region> grid_regions(int rows, int cols, const std::vector<int> &free_orig, int n)
{
  std::vector<int> local(n, -1);
  for (size_t t = 0; t < free_orig.size(); t++) local[free_orig[t]] = (int)t;
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
  return out;
}

std::vector<Region> prefix_regions(int m)
{
  std::vector<Region> out;
  for (int t = 1; t < m; t++) {
    Region r; r.name = "free<" + std::to_string(t);
    for (int a = 0; a < t; a++) r.loc.push_back(a);
    out.push_back(r);
  }
  return out;
}

}  // namespace

DNNCutResult dnn_cut_rounds(MAT *A, MAT *B, MAT *C, double shift, double cutoff,
                            const DNNParams &base, const std::vector<int> &free_orig,
                            int n_orig, const DNNCutParams &cp, DNNWarm *warm, double first_lb)
{
  const int m = A->m, n1 = m * m + 1;
  const Lift lift = {m, cp.side};
  std::vector<Region> regions;
  if (!cp.atomic) {
    if (cp.grid_rows > 0) regions = grid_regions(cp.grid_rows, cp.grid_cols, free_orig, n_orig);
    else if (cp.prefix) regions = prefix_regions(m);
  }
  const bool regional = !cp.atomic;
  if (cp.verbose >= 2)
    for (auto &R : regions) printf("  region %-8s capacity %zu\n", R.name.c_str(), R.loc.size());

  DNNCutResult out;
  out.first_lb = first_lb; out.best_lb = first_lb; out.status = -1;
  out.rounds = 0; out.iterations = 0; out.cpu = 0; out.active = 0; out.stalled = 0;
  out.nregions = (int)regions.size();
  std::vector<double> best_hist;   // best after each round (index = round)
  if (regional && regions.empty()) return out;

  std::vector<Cut> cuts;
  std::map<std::vector<int>, int> index;
  std::vector<double> Ladd((size_t)n1 * n1, 0.0);
  std::vector<std::vector<double>> P, WW;         // regional aggregates of the last Y
  DNNWarm save;
  double prev_lb = -HUGE_VAL, theta_cur = cp.theta;

  double last_primal = HUGE_VAL;
  std::vector<double> best_mu;   // the multipliers of the best round's solve (freeze_best)
  DNNWarm best_save;
  double best_mu_lb = -HUGE_VAL;
  // one solve of the modified problem with the current multipliers
  auto solve = [&](int r, int iters) {
    std::fill(Ladd.begin(), Ladd.end(), 0.0);
    double abs_add = 0, hsum = 0;
    int active = 0;
    for (const Cut &t : cuts)
      if (t.mu > 0) {
        abs_add += t.s < 0 ? add_cut(t, Ladd, n1, t.mu) : add_rcut(t, regions[t.s], Ladd, lift, n1, t.mu);
        if (t.type == (t.s < 0 ? 2 : 3)) hsum += t.mu;
        active++;
      }
    out.active = active;
    DNNParams par = base;
    par.max_iter = iters;
    par.hopeless_after = 0;
    par.hopeless_gap_ratio = 0;
    par.L_add = active ? Ladd.data() : NULL;
    par.L_add_abs = abs_add;
    DNNResult res = dnn_bound(A, B, C, shift - hsum, cutoff, par, (r > 0 && warm) ? warm : NULL, &save);
    if (warm) *warm = save; else warm = &save;
    out.cpu += res.cpu;
    if (r == 0) { out.first_lb = first_lb = res.lb; }
    else { out.rounds++; out.iterations += res.iterations; }
    out.best_lb = std::max(out.best_lb, res.lb);
    if (cp.freeze_best && res.lb > best_mu_lb) {
      best_mu_lb = res.lb;
      best_mu.resize(cuts.size());
      for (size_t t = 0; t < cuts.size(); t++) best_mu[t] = cuts[t].mu;
      best_save = save;
    }
    out.status = res.status;
    last_primal = res.primal;
    return res;
  };

  // round 0 is the uncut solve: run here (first_iters > 0) or the caller's
  // warm start, whose Y is separated; rounds 1.. solve with the cuts
  const int last_round = cp.freeze_after > 0 ? std::min(cp.rounds, cp.freeze_after) : cp.rounds;
  for (int r = 0; r <= last_round; r++) {
    DNNResult res;
    res.lb = -HUGE_VAL; res.primal = HUGE_VAL; res.status = -1; res.iterations = 0; res.cpu = 0;
    const bool solved = r > 0 || cp.first_iters > 0;
    if (solved) res = solve(r, r == 0 ? cp.first_iters : cp.iters);
    // separate on the last primal iterate
    const std::vector<double> &Y = warm->Y;
    if (Y.size() != (size_t)n1 * n1) break;
    std::vector<std::pair<double, Cut>> found;
    size_t nviol = 0;
    auto bound_memory = [&]() {
      if (found.size() > 4000000) {
        size_t keep = std::max<size_t>(100000, 8 * (size_t)cp.kcuts);
        std::nth_element(found.begin(), found.begin() + keep, found.end(),
                         [](const std::pair<double, Cut> &u, const std::pair<double, Cut> &w) {
                           return u.first > w.first; });
        found.resize(keep);
      }
    };
    if (regional) {
      P.assign(regions.size(), {}); WW.assign(regions.size(), {});
      for (size_t s = 0; s < regions.size(); s++) {
        aggregate(regions[s], Y, lift, n1, P[s], WW[s]);
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
        if (cp.verbose > 1)
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
      if (added >= cp.kcuts) break;
      std::vector<int> key = {f.second.type, f.second.a, f.second.b, f.second.c, f.second.s};
      if (index.count(key)) continue;
      index[key] = (int)cuts.size();
      f.second.mu = cp.mu0;
      cuts.push_back(f.second);
      added++;
    }
    const double maxv = found.empty() ? 0 : found[0].first;
    if (cp.verbose && solved)
      printf("round %2d  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f"
             "  active %d  violated %zu (max %.4f)  added %d\n", r, res.lb, res.lb - first_lb,
             out.best_lb - first_lb, res.primal, dnn_status_name(res.status), res.iterations,
             res.cpu, out.active, nviol, maxv, added);
    else if (cp.verbose)
      printf("round  0  lb %.4f (the node's bound)  violated %zu (max %.4f)  added %d\n",
             first_lb, nviol, maxv, added);
    if (solved && res.status == DNN_FATHOMED) break;
    if (r == last_round) break;
    if (cuts.empty()) break;
    best_hist.push_back(out.best_lb);
    if (cp.stall_rounds > 0 && r >= cp.stall_rounds && std::isfinite(cutoff)) {
      const double gain = out.best_lb - best_hist[r - cp.stall_rounds];
      if (cutoff - out.best_lb > cp.stall_factor * std::max(gain, 0.0)) {
        out.stalled = 1;
        if (cp.verbose) printf("stalled: gap %.4f, gain of the last %d rounds %.4f\n",
                               cutoff - out.best_lb, cp.stall_rounds, gain);
        break;
      }
    }
    // Polyak subgradient step on all cuts toward the target; the step factor
    // halves whenever a round lowers the bound (an overshoot)
    if (cp.theta <= 0) continue;
    const double lb_now = solved ? res.lb : first_lb;
    if (solved && r > 0 && lb_now < prev_lb - 1e-9) theta_cur *= 0.5;
    prev_lb = lb_now;
    double s2 = 0;
    std::vector<double> sg(cuts.size());
    for (size_t t = 0; t < cuts.size(); t++) {
      sg[t] = cuts[t].s < 0 ? violation(cuts[t], Y, n1)
                            : rviolation(cuts[t], P[cuts[t].s], WW[cuts[t].s], m);
      if (cuts[t].mu > 0 || sg[t] > 0) s2 += sg[t] * sg[t];
    }
    if (s2 <= 0) break;
    const double target = cp.cap_at_cutoff ? std::min(cutoff, out.best_lb + cp.delta)
                                           : out.best_lb + cp.delta;
    const double alpha = theta_cur * std::max(target - lb_now, 1e-6) / s2;
    for (size_t t = 0; t < cuts.size(); t++)
      cuts[t].mu = std::max(0.0, cuts[t].mu + alpha * sg[t]);
  }
  // the quiet round: the multipliers frozen, no new cuts, so that the
  // certificate catches up with the Lagrangian value of the final multipliers
  if (cp.freeze_iters > 0 && out.best_lb <= cutoff && !cuts.empty() &&
      (!cp.freeze_if_primal_above || last_primal > cutoff)) {
    if (cp.freeze_best && !best_mu.empty() && best_mu_lb > -HUGE_VAL) {
      // back to the best round's multipliers and iterate (cuts added later keep mu 0)
      for (size_t t = 0; t < cuts.size(); t++) cuts[t].mu = t < best_mu.size() ? best_mu[t] : 0.0;
      *warm = best_save;
    }
    DNNResult res = solve(1, cp.freeze_iters);
    if (cp.verbose)
      printf("quiet round  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f  active %d\n",
             res.lb, res.lb - first_lb, out.best_lb - first_lb, res.primal,
             dnn_status_name(res.status), res.iterations, res.cpu, out.active);
  }
  if (out.best_lb > cutoff) out.status = DNN_FATHOMED;
  return out;
}
