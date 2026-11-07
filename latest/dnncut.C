// dnncut.C: triangle cuts on top of the DNN bound (see dnncut.h); the code of
// experiments/dnn_cuts/dnncuts.C made a library, with regions on either side.
#include "dnncut.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <array>
#include <cfloat>
#include <ctime>
#include <random>

DNNCutParams::DNNCutParams()
{
  warm_primal = 0;
  rounds = 12; iters = 300; first_iters = 0; kcuts = 1000; kcuts5 = 0; add_until = 0;
  kcuts7 = 0; sub7 = 1000;
  kbox = 0; box_qmax = 24; box_starts = 40; box_bmax = 7;
  delta = 50; theta = 1.0; damp1 = 1.0; damp_rest = 1.0; mu0 = 0;
  grid_rows = grid_cols = 0; side = 0; prefix = false; atomic = false; verbose = 0;
  stall_rounds = 0; stall_factor = 3.0;
  pgate_rounds = 0; pgate_margin = 20;
  cap_at_cutoff = 1; freeze_after = 0; freeze_iters = 0; freeze_if_primal_above = 0;
  freeze_best = 1;
  L_ext = NULL; L_ext_abs = 0; ext_const = 0;
  pool_in = NULL; pool_out = NULL; other_orig = NULL; fixed_at = NULL; pool_check = 0;
  pool_all_round = -1; pool_all_out = NULL;
  cp_every = 0; cp_k = 0; cp_iters = 0; cp_sweeps = 1; cp_passes = 1;
}

namespace {

// s < 0: atomic, a, b, c lifted indices; s >= 0: region s, a, b, c free
// elements of the other side (local indices), type 0/1/2 = apex a/b/c, 3 =
// the sum inequality; type 4 = the 5-point count inequality on a < b < c <
// d < e with flip set `mask` (bit x = element x of a..e flipped) and alpha
// type 6 = a local-7 inequality: a is the index into the pool's L7Cut vector;
// type 8 = a box inequality: a is the index into the pool's BoxCut vector
struct Cut { int type, a, b, c, s; double mu; int d, e, mask, alpha; };

// ---- box (count) inequalities (type 8), see dnncut.h:
//   b sum_{i in Q} p_i - sum_{i<j in Q} W_ij - b (b + 1) / 2 <= 0
struct BoxCut { std::vector<int> Q; int b; };
static std::vector<BoxCut> *g_box = NULL;   // the current call's pool of box cuts
double box_value(const BoxCut &c, const std::vector<double> &p, const std::vector<double> &W, int m)
{
  double v = -0.5 * c.b * (c.b + 1);
  for (size_t x = 0; x < c.Q.size(); x++) {
    v += c.b * p[c.Q[x]];
    const double *w = &W[(size_t)c.Q[x] * m];
    for (size_t y = x + 1; y < c.Q.size(); y++) v -= w[c.Q[y]];
  }
  return v;
}
// best-improvement add/drop local search from the set `in` for one b: the
// gain of adding i is b p_i - sum_{k in Q} W_ik, of dropping i its negative
// (W_ii excluded); sizes kept in [3, qmax]
void box_local(const std::vector<double> &p, const std::vector<double> &W, int m, int b, int qmax,
               std::vector<char> &in)
{
  std::vector<double> Wm(m, 0.0);
  int sz = 0;
  for (int k = 0; k < m; k++)
    if (in[k]) {
      sz++;
      const double *w = &W[(size_t)k * m];
      for (int i = 0; i < m; i++) if (i != k) Wm[i] += w[i];
    }
  for (int moves = 0; moves < 20 * m; moves++) {
    double best = 1e-10; int bi = -1;
    for (int i = 0; i < m; i++) {
      const double g = b * p[i] - Wm[i];
      if (in[i]) { if (sz > 3 && -g > best) { best = -g; bi = i; } }
      else if (sz < qmax && g > best) { best = g; bi = i; }
    }
    if (bi < 0) return;
    const double sgn = in[bi] ? -1.0 : 1.0;
    in[bi] = !in[bi]; sz += in[bi] ? 1 : -1;
    const double *w = &W[(size_t)bi * m];
    for (int i = 0; i < m; i++) if (i != bi) Wm[i] += sgn * w[i];
  }
}

// ---- local-7 cuts (type 6), see dnncut.h
struct L7Cut { int f[7]; double ca[7], cb[21], beta; };
static std::vector<L7Cut> *g_l7 = NULL;   // the current call's pool of local-7 cuts
static int P7[21][2];
static void init_p7() { int q = 0; for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++) { P7[q][0] = i; P7[q][1] = j; q++; } }
// the admissible configurations of a region of ksz locations among m, as 28-vectors
struct L7Space { std::vector<std::array<double, 28>> V; };
L7Space l7_space(int ksz, int m)
{
  L7Space sp;
  for (int c = 0; c < 128; c++) {
    int s = __builtin_popcount(c);
    if (s > ksz || 7 - s > m - ksz) continue;
    std::array<double, 28> v;
    for (int i = 0; i < 7; i++) v[i] = (c >> i) & 1;
    for (int q = 0; q < 21; q++) v[7 + q] = v[P7[q][0]] * v[P7[q][1]];
    sp.V.push_back(v);
  }
  return sp;
}
// separate y = (p_T, W_TT) from the hull of sp.V: Frank-Wolfe projection, the
// direction y - x scaled to max |coefficient| 1, beta exact over the vertices
bool l7_separate(const L7Space &sp, const double *y, L7Cut &out, double &viol)
{
  const size_t nv = sp.V.size();
  if (nv == 0) return false;
  std::array<double, 28> x;
  double best = HUGE_VAL; size_t c0 = 0;
  for (size_t c = 0; c < nv; c++) {
    double d2 = 0; for (int k = 0; k < 28; k++) d2 += (sp.V[c][k] - y[k]) * (sp.V[c][k] - y[k]);
    if (d2 < best) { best = d2; c0 = c; }
  }
  // away-step Frank-Wolfe on the simplex weights (linear convergence on a
  // polytope): x = sum lambda_c V_c, exact line search, stop at a tiny gap
  std::vector<double> lam(nv, 0.0);
  lam[c0] = 1.0; x = sp.V[c0];
  for (int it = 0; it < 5000; it++) {
    double gmin = HUGE_VAL, gmax = -HUGE_VAL; size_t cs = 0, ca = 0;
    for (size_t c = 0; c < nv; c++) {
      double g = 0; for (int k = 0; k < 28; k++) g += (x[k] - y[k]) * sp.V[c][k];
      if (g < gmin) { gmin = g; cs = c; }
      if (lam[c] > 0 && g > gmax) { gmax = g; ca = c; }
    }
    double gx = 0; for (int k = 0; k < 28; k++) gx += (x[k] - y[k]) * x[k];
    const double fw_gap = gx - gmin, away_gap = gmax - gx;
    if (fw_gap < 1e-10) break;
    std::array<double, 28> d; double gmax_step;
    bool away = away_gap > fw_gap;
    if (!away) { for (int k = 0; k < 28; k++) d[k] = sp.V[cs][k] - x[k]; gmax_step = 1.0; }
    else { for (int k = 0; k < 28; k++) d[k] = x[k] - sp.V[ca][k]; gmax_step = lam[ca] / (1.0 - lam[ca] + 1e-300); }
    double gd = 0, dd = 0;
    for (int k = 0; k < 28; k++) { gd += (x[k] - y[k]) * d[k]; dd += d[k] * d[k]; }
    if (dd <= 0) break;
    const double step = std::min(gmax_step, std::max(0.0, -gd / dd));
    if (step <= 0) break;
    if (!away) { for (size_t c = 0; c < nv; c++) lam[c] *= (1 - step); lam[cs] += step; }
    else { for (size_t c = 0; c < nv; c++) lam[c] *= (1 + step); lam[ca] -= step; if (lam[ca] < 1e-15) lam[ca] = 0; }
    for (int k = 0; k < 28; k++) x[k] += step * d[k];
  }
  double dmax = 0;
  for (int k = 0; k < 28; k++) dmax = std::max(dmax, fabs(y[k] - x[k]));
  if (dmax < 1e-7) return false;
  for (int k = 0; k < 7; k++) out.ca[k] = (y[k] - x[k]) / dmax;
  for (int q = 0; q < 21; q++) out.cb[q] = (y[7 + q] - x[7 + q]) / dmax;
  double beta = -HUGE_VAL, ay = 0;
  for (size_t c = 0; c < nv; c++) {
    double v = 0;
    for (int k = 0; k < 7; k++) v += out.ca[k] * sp.V[c][k];
    for (int q = 0; q < 21; q++) v += out.cb[q] * sp.V[c][7 + q];
    beta = std::max(beta, v);
  }
  for (int k = 0; k < 7; k++) ay += out.ca[k] * y[k];
  for (int q = 0; q < 21; q++) ay += out.cb[q] * y[7 + q];
  out.beta = beta; viol = ay - beta;
  return viol > 1e-6;
}
double rviolation7(const Cut &t, const std::vector<double> &p, const std::vector<double> &W, int m)
{
  const L7Cut &c = (*g_l7)[t.a];
  double v = -c.beta;
  for (int k = 0; k < 7; k++) v += c.ca[k] * p[c.f[k]];
  for (int q = 0; q < 21; q++) v += c.cb[q] * W[(size_t)c.f[P7[q][0]] * m + c.f[P7[q][1]]];
  return v;
}
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

// ---- 5-point count inequalities (type 4). For T = {a..e}, the flip set F
// and alpha, with z'_i = 1 - z_i on F and q' = sum_T z'_i, the inequality
// (q' - alpha)(q' - alpha - 1) >= 0 reads, in the moments p_i = E z_i and
// W_ij = E z_i z_j,
//   sum_i cp_i p_i + sum_{i<j} cw_ij W_ij + c0 <= 0,
//   cp_i = sg_i (alpha - |F \ {i}|),  cw_ij = -sg_i sg_j,
//   c0 = alpha |F| - C(|F|, 2) - alpha (alpha + 1) / 2,   sg_i = -1 on F, +1 off it
// (expand p'_i = 1 - p_i on F, W'_ij = p_j - W_ij for i in F only, 1 - p_i -
// p_j + W_ij for both). |T| = 3, alpha = 1 gives the four triangles. Pairs
// are indexed (0,1) (0,2) (0,3) (0,4) (1,2) (1,3) (1,4) (2,3) (2,4) (3,4).
struct Five { double cp[5], cw[10], c0; };
void five_coef(int mask, int alpha, Five &f)
{
  int nf = 0;
  for (int i = 0; i < 5; i++) nf += (mask >> i) & 1;
  double sg[5];
  for (int i = 0; i < 5; i++) {
    const int fi = (mask >> i) & 1;
    sg[i] = fi ? -1.0 : 1.0;
    f.cp[i] = sg[i] * (alpha - (nf - fi));
  }
  int t = 0;
  for (int i = 0; i < 5; i++)
    for (int j = i + 1; j < 5; j++, t++) f.cw[t] = -sg[i] * sg[j];
  f.c0 = alpha * nf - nf * (nf - 1) / 2 - alpha * (alpha + 1) / 2;
}

// the 48 distinct patterns: alpha = 1 with every flip set, alpha = 2 with the
// flip sets not containing the first element (its complement with alpha 2 is
// the same inequality; alpha = 3 is the complement of alpha = 1)
struct FivePatterns {
  std::vector<std::pair<int, int>> pat;   // (mask, alpha)
  std::vector<Five> f;
  FivePatterns() {
    for (int alpha = 1; alpha <= 2; alpha++)
      for (int mask = 0; mask < 32; mask++) {
        if (alpha == 2 && (mask & 1)) continue;
        Five x; five_coef(mask, alpha, x);
        pat.push_back({mask, alpha}); f.push_back(x);
      }
  }
};

double rviolation5(const Cut &t, const std::vector<double> &p, const std::vector<double> &W, int m)
{
  Five f; five_coef(t.mask, t.alpha, f);
  const int T[5] = {t.a, t.b, t.c, t.d, t.e};
  double v = f.c0;
  int k = 0;
  for (int i = 0; i < 5; i++) {
    v += f.cp[i] * p[T[i]];
    for (int j = i + 1; j < 5; j++, k++) v += f.cw[k] * W[(size_t)T[i] * m + T[j]];
  }
  return v;
}

// the constant of a cut lin(Y) <= rhs, moved into the shift as mu * rhs
double cut_rhs(const Cut &t)
{
  if (t.s < 0) return t.type == 2 ? 1.0 : 0.0;
  if (t.type == 3) return 1.0;
  if (t.type == 4) { Five f; five_coef(t.mask, t.alpha, f); return -f.c0; }
  if (t.type == 6) return (*g_l7)[t.a].beta;
  if (t.type == 8) { const BoxCut &c = (*g_box)[t.a]; return 0.5 * c.b * (c.b + 1); }
  return 0.0;
}

double rviolation(const Cut &t, const std::vector<double> &p, const std::vector<double> &W, int m)
{
  if (t.type == 4) return rviolation5(t, p, W, m);
  if (t.type == 6) return rviolation7(t, p, W, m);
  if (t.type == 8) return box_value((*g_box)[t.a], p, W, m);
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
  } else if (t.type == 3) {
    lin(i, mu); lin(j, mu); lin(k, mu);
    pair(i, j, -mu); pair(i, k, -mu); pair(j, k, -mu);
  } else if (t.type == 6) {
    const L7Cut &c = (*g_l7)[t.a];
    for (int x = 0; x < 7; x++) if (c.ca[x] != 0) lin(c.f[x], mu * c.ca[x]);
    for (int q = 0; q < 21; q++) if (c.cb[q] != 0) pair(c.f[P7[q][0]], c.f[P7[q][1]], mu * c.cb[q]);
  } else {
    Five f; five_coef(t.mask, t.alpha, f);
    const int T[5] = {t.a, t.b, t.c, t.d, t.e};
    int q = 0;
    for (int x = 0; x < 5; x++) {
      if (f.cp[x] != 0) lin(T[x], mu * f.cp[x]);
      for (int y = x + 1; y < 5; y++, q++) pair(T[x], T[y], mu * f.cw[q]);
    }
  }
  return abs_;
}

// the summed box terms of one region: coefficient pc_i on p_i and wc_ij
// (i < j, upper triangle of an m x m array) on W_ij, written as add_rcut does
double add_region_terms(const Region &S, const std::vector<double> &pc, const std::vector<double> &wc,
                        std::vector<double> &L, const Lift &lift, int n1)
{
  const int m = lift.m;
  double abs_ = 0;
  for (int i = 0; i < m; i++) {
    if (pc[i] != 0)
      for (int a : S.loc) { const size_t r = lift(a, i); L[r * n1 + r] += pc[i]; abs_ += fabs(pc[i]); }
    for (int j = i + 1; j < m; j++) {
      const double v = wc[(size_t)i * m + j];
      if (v == 0) continue;
      for (int a : S.loc)
        for (int b : S.loc) {
          if (a == b) continue;   // gangster zero
          const size_t r = lift(a, i), c = lift(b, j);
          L[r * n1 + c] += v / 2; L[c * n1 + r] += v / 2; abs_ += fabs(v);
        }
    }
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

// the objective term of the cuts with mu > 0 added into L (the box terms
// summed per region first), its sum of |entries| into abs_add, the constant
// sum mu h into hsum; returns the number of active cuts
int assemble(const std::vector<Cut> &cuts, const std::vector<Region> &regions,
             const std::vector<BoxCut> &boxes, const Lift &lift, int n1, std::vector<double> &L,
             double &abs_add, double &hsum)
{
  const int m = lift.m;
  int active = 0;
  std::vector<std::vector<double>> bpc, bwc;   // per region, the summed box terms
  for (const Cut &t : cuts)
    if (t.mu > 0) {
      if (t.type == 8) {
        if (bpc.empty()) { bpc.assign(regions.size(), {}); bwc.assign(regions.size(), {}); }
        if (bpc[t.s].empty()) { bpc[t.s].assign(m, 0.0); bwc[t.s].assign((size_t)m * m, 0.0); }
        const BoxCut &c = boxes[t.a];
        for (size_t x = 0; x < c.Q.size(); x++) {
          bpc[t.s][c.Q[x]] += t.mu * c.b;
          for (size_t y = x + 1; y < c.Q.size(); y++) {
            const int i = std::min(c.Q[x], c.Q[y]), j = std::max(c.Q[x], c.Q[y]);
            bwc[t.s][(size_t)i * m + j] -= t.mu;
          }
        }
      } else
        abs_add += t.s < 0 ? add_cut(t, L, n1, t.mu) : add_rcut(t, regions[t.s], L, lift, n1, t.mu);
      hsum += t.mu * cut_rhs(t);
      active++;
    }
  for (size_t s = 0; s < bpc.size(); s++)
    if (!bpc[s].empty()) abs_add += add_region_terms(regions[s], bpc[s], bwc[s], L, lift, n1);
  return active;
}

// a parent's pool remapped to this node (dnncut.h, DNNPool): regions to the
// node's free region-side elements, matched to `regions` or appended (empty or
// full ones dropped with their cuts); a triangle with a fixed member dropped; a
// box's fixed members leave Q, those fixed inside the region lowering b by one,
// kept while |Q| >= 3 and 1 <= b <= |Q| - 1. Appends to cuts (and boxes) and
// returns the number kept.
int remap_pool(const DNNPool &in, const std::vector<int> &free_orig, int n_orig, const DNNCutParams &cp,
               std::vector<Region> &regions, std::vector<Cut> &cuts, std::vector<BoxCut> &boxes,
               int &dropped, bool keep_zero = false)
{
  const int m = (int)free_orig.size();
  std::vector<int> local(n_orig, -1), olocal(n_orig, -1);
  for (int t = 0; t < m; t++) local[free_orig[t]] = t;
  for (size_t t = 0; t < cp.other_orig->size(); t++) olocal[(*cp.other_orig)[t]] = (int)t;
  std::vector<int> rmap(in.regions.size(), -1);
  std::vector<std::vector<char>> inreg(in.regions.size());
  for (size_t s = 0; s < in.regions.size(); s++) {
    std::vector<int> loc;
    inreg[s].assign(n_orig, 0);
    for (int x : in.regions[s])
      if (x >= 0 && x < n_orig) { inreg[s][x] = 1; if (local[x] >= 0) loc.push_back(local[x]); }
    std::sort(loc.begin(), loc.end());
    if (loc.empty() || (int)loc.size() == m) continue;
    size_t k = 0;
    while (k < regions.size() && regions[k].loc != loc) k++;
    if (k == regions.size()) regions.push_back({"pool" + std::to_string(s), loc});
    rmap[s] = (int)k;
  }
  int kept = 0;
  dropped = 0;
  for (const DNNPoolCut &c : in.cuts) {
    if (c.region < 0 || c.region >= (int)rmap.size() || rmap[c.region] < 0 || !(keep_zero ? c.mu >= 0 : c.mu > 0)) { dropped++; continue; }
    bool ok = true;
    for (int x : c.mem) ok = ok && x >= 0 && x < n_orig;
    if (!ok) { dropped++; continue; }
    if (c.type >= 0 && c.type <= 3 && c.mem.size() == 3) {
      const int a = olocal[c.mem[0]], b = olocal[c.mem[1]], d = olocal[c.mem[2]];
      if (a < 0 || b < 0 || d < 0 || !(a < b && b < d)) { dropped++; continue; }
      Cut t = {c.type, a, b, d, rmap[c.region], c.mu, 0, 0, 0, 0};
      cuts.push_back(t); kept++;
    } else if (c.type == 8) {
      BoxCut bc; bc.b = c.b;
      for (int x : c.mem) {
        if (olocal[x] >= 0) bc.Q.push_back(olocal[x]);
        else { const int at = (*cp.fixed_at)[x]; if (at >= 0 && at < n_orig && inreg[c.region][at]) bc.b--; }
      }
      std::sort(bc.Q.begin(), bc.Q.end());
      if (bc.Q.size() < 3 || bc.b < 1 || bc.b > (int)bc.Q.size() - 1) { dropped++; continue; }
      boxes.push_back(bc);
      Cut t = {8, (int)boxes.size() - 1, 0, 0, rmap[c.region], c.mu, 0, 0, 0, 0};
      cuts.push_back(t); kept++;
    } else dropped++;
  }
  return kept;
}

// the pool of these cuts at multipliers mu (or the cuts' own when NULL), in
// original indices; only regional triangles and boxes with mu > 0
void export_pool(const std::vector<Cut> &cuts, const std::vector<double> *mu, const std::vector<Region> &regions,
                 const std::vector<BoxCut> &boxes, const std::vector<int> &free_orig,
                 const std::vector<int> &other_orig, int side, DNNPool &out, bool keep_zero = false)
{
  out = DNNPool();
  out.side = side;
  std::vector<int> rid(regions.size(), -1);
  for (size_t t = 0; t < cuts.size(); t++) {
    const Cut &c = cuts[t];
    const double u = mu ? (t < mu->size() ? (*mu)[t] : 0.0) : c.mu;
    if (!(keep_zero ? u >= 0 : u > 0) || c.s < 0 || !((c.type >= 0 && c.type <= 3) || c.type == 8)) continue;
    if (rid[c.s] < 0) {
      rid[c.s] = (int)out.regions.size();
      std::vector<int> r;
      for (int a : regions[c.s].loc) r.push_back(free_orig[a]);
      std::sort(r.begin(), r.end());
      out.regions.push_back(r);
    }
    DNNPoolCut pc;
    pc.type = c.type; pc.region = rid[c.s]; pc.mu = u; pc.b = 0;
    if (c.type == 8) {
      const BoxCut &bc = boxes[c.a];
      pc.b = bc.b;
      for (int x : bc.Q) pc.mem.push_back(other_orig[x]);
    } else
      pc.mem = {other_orig[c.a], other_orig[c.b], other_orig[c.c]};
    out.cuts.push_back(pc);
  }
}

// the largest <L, v v'> - hsum over K random permutations of the node
// (v = (1, x), lifted index 1 + (first index) m + (second index)): every cut is
// valid, so a positive value means the term was assembled wrong
double check_term(const std::vector<double> &L, double hsum, int m, int K, unsigned seed)
{
  const int n1 = m * m + 1;
  std::mt19937 rng(seed);
  std::vector<int> pi(m), idx(m + 1);
  for (int i = 0; i < m; i++) pi[i] = i;
  double worst = -HUGE_VAL;
  for (int k = 0; k < K; k++) {
    std::shuffle(pi.begin(), pi.end(), rng);
    idx[0] = 0;
    for (int a = 0; a < m; a++) idx[a + 1] = 1 + a * m + pi[a];
    double v = 0;
    for (int r : idx) for (int c : idx) v += L[(size_t)r * n1 + c];
    worst = std::max(worst, v - hsum);
  }
  return worst;
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

// ---- CP-ADMM (DNNYCuts in dnn.h; experiments/wil50_feasibility/roadmap/admm/cpadmm)
//
// A regional cut is lin(Y) <= h with lin(Y) = sum_i cp_i p_i + sum_{i<j} cw_ij W_ij
// in the moments of its region S (aggregate(); W symmetrised): lin(Y) = <G, Y>
// for the matrix G that add_rcut writes at mu = 1, ||G||^2 = |S| sum cp^2 +
// |S|(|S|-1)/2 sum cw^2. The Y-step projects X onto P ∩ {G_t Y <= h_t} by
// Dykstra's algorithm with the sets P (the box and fixed entries, as dnn.C's
// Y-step) and the half-spaces, in cyclic order. A half-space's correction is
// nu_t G_t with nu_t >= 0, and P comes first in a sweep, so with
// y + p_P + sum_t nu_t G_t = X the P-step is y = P(X - G'nu) whatever p_P was:
// the warm start is nu alone. The half-space steps (Hildreth's
// nu_t <- max(0, nu_t + (<G_t, y> - h_t) / ||G_t||^2), y -= d G_t) run on the
// regions' moments, kept exact under y -= d G_t: region S' gains -d cp_i k on
// p_i and -d cw_ij k(k-1)/2 on W_ij, k = |S ∩ S'|. The dense passes group the
// (a, b) blocks of Y by the set of regions holding both a and b (sigma), so
// each costs one pass over Y. Either side: on side 0 the region side is the
// first index (row 1 + a m + i, the (a, b) block contiguous), on side 1 the
// second (row 1 + i m + a, the (a, b) block the stride-m columns 1 + j m + b);
// only the two dense passes (pstep, dense_term) depend on it. Side 1 visits the
// same (a, b, i, j) in the same order of a and b, so the projector of side 1 on
// a transposed X is side 0's to the bit (cp_selftest (f)).
struct CPTerm { int i, j; double c; };          // j < 0: c p_i; else i < j: c W_ij
struct CPCut { int s, t0, t1, src; double h, g2, absrow; };   // terms [t0, t1); src: index in the pool

// the terms and right-hand side of a regional cut (the sign convention of
// add_rcut and cut_rhs), pairs normalised to i < j and merged
void cp_cut_terms(const Cut &t, std::vector<CPTerm> &out, double &h)
{
  std::map<std::pair<int, int>, double> acc;
  auto lin = [&](int i, double v) { acc[{i, -1}] += v; };
  auto pair = [&](int i, int j, double v) { acc[{std::min(i, j), std::max(i, j)}] += v; };
  const int i = t.a, j = t.b, k = t.c;
  if (t.type < 3) {
    const int apex = t.type == 0 ? i : (t.type == 1 ? j : k);
    const int o1 = t.type == 0 ? j : i, o2 = t.type == 2 ? j : k;
    pair(apex, o1, 1); pair(apex, o2, 1); pair(o1, o2, -1); lin(apex, -1);
  } else if (t.type == 3) {
    lin(i, 1); lin(j, 1); lin(k, 1); pair(i, j, -1); pair(i, k, -1); pair(j, k, -1);
  } else if (t.type == 6) {
    const L7Cut &c = (*g_l7)[t.a];
    for (int x = 0; x < 7; x++) if (c.ca[x] != 0) lin(c.f[x], c.ca[x]);
    for (int q = 0; q < 21; q++) if (c.cb[q] != 0) pair(c.f[P7[q][0]], c.f[P7[q][1]], c.cb[q]);
  } else if (t.type == 8) {
    const BoxCut &c = (*g_box)[t.a];
    for (size_t x = 0; x < c.Q.size(); x++) {
      lin(c.Q[x], c.b);
      for (size_t y = x + 1; y < c.Q.size(); y++) pair(c.Q[x], c.Q[y], -1);
    }
  } else {
    Five f; five_coef(t.mask, t.alpha, f);
    const int T[5] = {t.a, t.b, t.c, t.d, t.e};
    int q = 0;
    for (int x = 0; x < 5; x++) {
      if (f.cp[x] != 0) lin(T[x], f.cp[x]);
      for (int y = x + 1; y < 5; y++, q++) pair(T[x], T[y], f.cw[q]);
    }
  }
  h = cut_rhs(t);
  out.clear();
  for (auto &e : acc) if (e.second != 0) out.push_back({e.first.first, e.first.second, e.second});
}

class CutProjector : public DNNYCuts {
public:
  int m = 0, n1 = 0, R = 0;
  int side = 0;                    // the region side: 0 the first lifted index, 1 the second (Lift)
  int sweeps = 1, passes = 1;      // Dykstra sweeps per Y-step, half-space passes per sweep
  double omega = 1.0;              // over-relaxation of the half-space steps (1: Hildreth/Dykstra)
  std::vector<double> init_mu;     // optional start multipliers (unscaled), applied at the first Y-step
  double init_scale = 0;           // as nu = mu / (scale beta)
  // separation inside the solve (DNNCPParams::sep_every): every sep_every Y-steps, up to
  // sep_k new violated regional triangles on the Y-step's output, ranked as in
  // dnn_cut_rounds, appended with nu = 0; pool_cap bounds the pool (zero-multiplier cuts
  // dropped first, the longest at zero first), drop_after drops cuts whose nu has been 0
  // for that many Y-steps at a separation
  int sep_every = 0, sep_k = 0, pool_cap = 0, drop_after = 0;
  std::vector<Cut> pool;           // the cuts (CPCut::src indexes it)
  std::set<std::vector<int>> keys; // the pool's triangles, as dnn_cut_rounds keys them
  std::vector<int> last_pos;       // per projector cut: the last separation at which nu > 0 (or its arrival)
  int n_sep = 0, last_found = 0, last_added = 0, n_added = 0, n_dropped = 0;
  double last_maxv = 0, t_sep = 0;
  // the stall rule per check (dnn_cut_cp): after stall_window iterations, stop when
  // cutoff - best > stall_factor * (best(t) - best(t - stall_window)), best >= first_lb
  int stall_window = 0;
  double stall_factor = 3.0, stall_cutoff = HUGE_VAL, stall_first = -HUGE_VAL;
  bool stalled = false;
  std::vector<std::pair<int, double>> best_hist;
  int diag_every = 0, diag_sweeps = 30;   // every diag_every Y-steps: distance to a longer run
  std::vector<std::vector<int>> rloc;
  std::vector<CPCut> cuts;
  std::vector<CPTerm> terms;
  std::vector<double> nu;          // Dykstra corrections: nu_t G_t, nu_t >= 0
  std::vector<int> sig2;           // m x m: the sigma of (a, b), -1 when no region holds both
  std::vector<std::vector<int>> sig_reg;     // regions of each sigma
  std::vector<std::vector<std::pair<int, int>>> ovl;   // per region: (region, overlap k >= 1)
  std::vector<double> Pag, Wag;    // moments of the current y, per region (W: i < j used)
  std::vector<double> CPs, CWs;    // sum_t nu_t c_t per region (the P-step's G'nu)
  std::vector<double> Csig, Dsig, Gsig, Gdsig;   // per sigma: 1/2 sum CW, sum CP; block sums of y
  double beta_last = 0;
  // diagnostics of the last Y-step (trace)
  double last_viol = 0, last_dstep = 0, last_rnorm = 0;
  int last_nviol = 0, last_nact = 0, ysteps = 0;
  double diag_dy = -1, diag_dnu = -1, diag_viol = -1, diag_ynorm = -1;
  int diag_it = -1;
  // CPU seconds: P-steps, half-space passes, the diagnostic runs (inside the hook's total)
  double t_p = 0, t_h = 0, t_diag = 0;
  int n_p = 0, n_h = 0;
  static double tcpu() { timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

  static std::vector<int> cut_key(const Cut &c) { return {c.type, c.a, c.b, c.c, c.s, c.d, c.e, c.mask, c.alpha}; }

  // the cut pool[ci] into the projector (terms, ||G||^2, magnitudes); false when it
  // has no terms or is not regional
  bool add_cut(int ci)
  {
    const Cut &c = pool[ci];
    if (c.s < 0 || c.s >= R) return false;   // atomic cuts are not supported
    std::vector<CPTerm> tt;
    double h;
    cp_cut_terms(c, tt, h);
    if (tt.empty()) return false;
    CPCut k;
    k.s = c.s; k.h = h; k.t0 = (int)terms.size(); k.src = ci;
    const double sz = (double)rloc[c.s].size();
    double g2 = 0, ab = 0;
    for (auto &x : tt) {
      terms.push_back(x);
      if (x.j < 0) { g2 += sz * x.c * x.c; ab += fabs(x.c) * sz; }
      else { g2 += 0.5 * sz * (sz - 1) * x.c * x.c; ab += fabs(x.c) * sz * (sz - 1); }
    }
    k.t1 = (int)terms.size(); k.g2 = g2; k.absrow = ab;
    if (!(g2 > 0)) { terms.resize(k.t0); return false; }
    cuts.push_back(k);
    nu.push_back(0.0);
    last_pos.push_back(0);
    if (c.type <= 3) keys.insert(cut_key(c));
    return true;
  }

  CutProjector(int m_, const std::vector<Region> &regions, const std::vector<Cut> &src_pool, int side_ = 0)
  {
    m = m_; n1 = m * m + 1; R = (int)regions.size(); side = side_ == 1 ? 1 : 0;
    for (auto &r : regions) rloc.push_back(r.loc);
    pool = src_pool;
    for (size_t ci = 0; ci < pool.size(); ci++) add_cut((int)ci);
    // region membership signatures, the sigma of every (a, b)
    std::vector<std::vector<char>> in(R, std::vector<char>(m, 0));
    for (int s = 0; s < R; s++) for (int a : rloc[s]) in[s][a] = 1;
    std::map<std::vector<int>, int> sid;
    sig2.assign((size_t)m * m, -1);
    for (int a = 0; a < m; a++)
      for (int b = 0; b < m; b++) {
        std::vector<int> key;
        for (int s = 0; s < R; s++) if (in[s][a] && in[s][b]) key.push_back(s);
        if (key.empty()) continue;
        auto it = sid.find(key);
        if (it == sid.end()) { it = sid.insert({key, (int)sig_reg.size()}).first; sig_reg.push_back(key); }
        sig2[(size_t)a * m + b] = it->second;
      }
    ovl.assign(R, {});
    for (int s = 0; s < R; s++)
      for (int s2 = 0; s2 < R; s2++) {
        int k = 0;
        for (int a = 0; a < m; a++) k += in[s][a] && in[s2][a];
        if (k > 0) ovl[s].push_back({s2, k});
      }
    const size_t mm = (size_t)m * m, ns = sig_reg.size();
    Pag.assign((size_t)R * m, 0); Wag.assign((size_t)R * mm, 0);
    CPs.assign((size_t)R * m, 0); CWs.assign((size_t)R * mm, 0);
    Csig.assign(ns * mm, 0); Dsig.assign(ns * m, 0); Gsig.assign(ns * mm, 0); Gdsig.assign(ns * m, 0);
  }

  size_t nsig() const { return sig_reg.size(); }

  // CPs, CWs from nu, exactly (the incremental sums drift by rounding)
  void rebuild_sums()
  {
    std::fill(CPs.begin(), CPs.end(), 0.0); std::fill(CWs.begin(), CWs.end(), 0.0);
    const size_t mm = (size_t)m * m;
    for (size_t t = 0; t < cuts.size(); t++) {
      if (!(nu[t] > 0)) continue;
      const CPCut &c = cuts[t];
      for (int q = c.t0; q < c.t1; q++) {
        const CPTerm &x = terms[q];
        if (x.j < 0) CPs[(size_t)c.s * m + x.i] += nu[t] * x.c;
        else CWs[c.s * mm + (size_t)x.i * m + x.j] += nu[t] * x.c;
      }
    }
  }

  // the P-step: Y = P(X - G'nu) (P as dnn.C's Y-step: clip to [0, 1], Y_00 = 1, the
  // gangster and eliminated zeros), and the regions' moments of Y
  void pstep(const double *X, double *Y, const unsigned char *elim)
  {
    const size_t mm = (size_t)m * m, ns = nsig();
    std::fill(Csig.begin(), Csig.end(), 0.0); std::fill(Dsig.begin(), Dsig.end(), 0.0);
    for (size_t g = 0; g < ns; g++) {
      double *cs = &Csig[g * mm], *ds = &Dsig[g * m];
      for (int s : sig_reg[g]) {
        const double *cw = &CWs[s * mm], *cpv = &CPs[(size_t)s * m];
        for (int i = 0; i < m; i++) {
          ds[i] += cpv[i];
          for (int j = i + 1; j < m; j++) cs[(size_t)i * m + j] += 0.5 * cw[(size_t)i * m + j];
        }
      }
      for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) cs[(size_t)j * m + i] = cs[(size_t)i * m + j];
    }
    std::fill(Gsig.begin(), Gsig.end(), 0.0); std::fill(Gdsig.begin(), Gdsig.end(), 0.0);
    auto clip = [](double v) { return std::min(1.0, std::max(0.0, v)); };
    Y[0] = 1.0;
    for (int c = 1; c < n1; c++) Y[c] = clip(X[c]);
    if (elim) for (int c = 1; c < n1; c++) if (elim[c - 1]) Y[c] = 0;
    if (side == 1) { pstep_rows1(X, Y, elim); moments_from_groups(); return; }
    for (int a = 0; a < m; a++)
      for (int i = 0; i < m; i++) {
        const size_t r = 1 + (size_t)a * m + i;
        const double *xr = X + r * n1;
        double *yr = Y + r * n1;
        const bool erow = elim && elim[r - 1];
        for (int b = 0; b < m; b++) {
          const int g = sig2[(size_t)a * m + b];
          const double *xb = xr + 1 + (size_t)b * m;
          double *yb = yr + 1 + (size_t)b * m;
          if (b == a) {
            for (int j = 0; j < m; j++) yb[j] = 0;   // gangster: same first index
            yb[i] = clip(g >= 0 ? xb[i] - Dsig[(size_t)g * m + i] : xb[i]);
          } else {
            if (g >= 0) { const double *cs = &Csig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) yb[j] = clip(xb[j] - cs[j]); }
            else for (int j = 0; j < m; j++) yb[j] = clip(xb[j]);
            yb[i] = 0;   // gangster: same second index
          }
        }
        yr[0] = clip(xr[0]);
        if (erow) for (int c = 0; c < n1; c++) yr[c] = 0;
        else if (elim) for (int c = 1; c < n1; c++) if (elim[c - 1]) yr[c] = 0;
        for (int b = 0; b < m; b++) {
          const int g = sig2[(size_t)a * m + b];
          if (g < 0) continue;
          const double *yb = yr + 1 + (size_t)b * m;
          if (b == a) Gdsig[(size_t)g * m + i] += yb[i];
          else { double *gs = &Gsig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) gs[j] += yb[j]; }
        }
      }
    moments_from_groups();
  }

  // the rows of the P-step on side 1: row r = 1 + i m + a (i the other side's element,
  // a the region side's), column 1 + j m + b; the side-0 block (a, b) is the columns of
  // second index b, at stride m. Gangster zeros: same first index (j = i, b != a) and
  // same second index (b = a, j != i). Csig, Dsig and the block sums Gsig, Gdsig are
  // side 0's, and each (i, j, sigma) receives its terms in the same order of (a, b).
  // Each row goes through two row buffers (rowx, rowy), so that memory is read and
  // written in sequence and the stride-m passes stay in L1, and the next row is
  // prefetched meanwhile: 10.3 ms per P-step against 16 ms in place and side 0's
  // 8.8 ms, m = 47, caches flushed (README.md, "Side 1")
  std::vector<double> rowx, rowy;
  void pstep_rows1(const double *X, double *Y, const unsigned char *elim)
  {
    const size_t mm = (size_t)m * m, sm = (size_t)m;
    auto clip = [](double v) { return std::min(1.0, std::max(0.0, v)); };
    rowx.resize(n1); rowy.resize(n1);
    const int per = (n1 + 8 * m - 1) / (8 * m);   // cache lines of the next row to prefetch per b
    for (int i = 0; i < m; i++)
      for (int a = 0; a < m; a++) {
        const size_t r = 1 + (size_t)i * m + a;
        const double *xr = X + r * n1;
        double *yr = Y + r * n1;
        const double *xn = r + 1 < (size_t)n1 ? xr + n1 : xr;
        double *yn = r + 1 < (size_t)n1 ? yr + n1 : yr;
        const bool erow = elim && elim[r - 1];
        std::copy(xr, xr + n1, rowx.begin());
        const double *xx = rowx.data() + 1;
        double *yy = rowy.data() + 1;
        for (int b = 0; b < m; b++) {
          for (int k = 0; k < per; k++) {
            const size_t o = ((size_t)b * per + k) * 8;
            if (o < (size_t)n1) { __builtin_prefetch(xn + o, 0, 3); __builtin_prefetch(yn + o, 1, 3); }
          }
          const int g = sig2[(size_t)a * m + b];
          const double *xc = xx + b;   // column (j, b) at xc[j m]
          double *yc = yy + b;
          if (b == a) {
            for (int j = 0; j < m; j++) yc[j * sm] = 0;   // gangster: same second index
            yc[i * sm] = clip(g >= 0 ? xc[i * sm] - Dsig[(size_t)g * m + i] : xc[i * sm]);
          } else {
            if (g >= 0) { const double *cs = &Csig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) yc[j * sm] = clip(xc[j * sm] - cs[j]); }
            else for (int j = 0; j < m; j++) yc[j * sm] = clip(xc[j * sm]);
            yc[i * sm] = 0;   // gangster: same first index
          }
        }
        rowy[0] = clip(rowx[0]);
        if (erow) std::fill(rowy.begin(), rowy.end(), 0.0);
        else if (elim) for (int c = 1; c < n1; c++) if (elim[c - 1]) rowy[c] = 0;
        for (int b = 0; b < m; b++) {
          const int g = sig2[(size_t)a * m + b];
          if (g < 0) continue;
          const double *yc = yy + b;
          if (b == a) Gdsig[(size_t)g * m + i] += yc[i * sm];
          else { double *gs = &Gsig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) gs[j] += yc[j * sm]; }
        }
        std::copy(rowy.begin(), rowy.end(), yr);
      }
  }

  void moments_from_groups()
  {
    const size_t mm = (size_t)m * m, ns = nsig();
    std::fill(Pag.begin(), Pag.end(), 0.0); std::fill(Wag.begin(), Wag.end(), 0.0);
    for (size_t g = 0; g < ns; g++) {
      const double *gs = &Gsig[g * mm], *gd = &Gdsig[g * m];
      for (int s : sig_reg[g]) {
        double *p = &Pag[(size_t)s * m], *w = &Wag[s * mm];
        for (int i = 0; i < m; i++) {
          p[i] += gd[i];
          for (int j = i + 1; j < m; j++) w[(size_t)i * m + j] += 0.5 * (gs[(size_t)i * m + j] + gs[(size_t)j * m + i]);
        }
      }
    }
  }

  // the cut's value <G_t, y> - h_t on the current moments
  double value(const CPCut &c) const
  {
    const double *p = &Pag[(size_t)c.s * m], *w = &Wag[(size_t)c.s * m * m];
    double v = -c.h;
    for (int q = c.t0; q < c.t1; q++) {
      const CPTerm &x = terms[q];
      v += x.c * (x.j < 0 ? p[x.i] : w[(size_t)x.i * m + x.j]);
    }
    return v;
  }

  // one cyclic pass of half-space steps on the moments; returns max |d| ||G||
  double hpass()
  {
    const size_t mm = (size_t)m * m;
    double dmax = 0;
    for (size_t t = 0; t < cuts.size(); t++) {
      const CPCut &c = cuts[t];
      const double v = value(c);
      const double nn = std::max(0.0, nu[t] + omega * v / c.g2), d = nn - nu[t];
      if (d == 0) continue;
      nu[t] = nn;
      dmax = std::max(dmax, fabs(d) * sqrt(c.g2));
      for (const auto &o : ovl[c.s]) {
        const double fk = o.second, fkk = 0.5 * o.second * (o.second - 1);
        double *p = &Pag[(size_t)o.first * m], *w = &Wag[o.first * mm];
        for (int q = c.t0; q < c.t1; q++) {
          const CPTerm &x = terms[q];
          if (x.j < 0) p[x.i] -= d * x.c * fk;
          else if (fkk > 0) w[(size_t)x.i * m + x.j] -= d * x.c * fkk;
        }
      }
      for (int q = c.t0; q < c.t1; q++) {
        const CPTerm &x = terms[q];
        if (x.j < 0) CPs[(size_t)c.s * m + x.i] += d * x.c;
        else CWs[c.s * mm + (size_t)x.i * m + x.j] += d * x.c;
      }
    }
    return dmax;
  }

  // the output's cut residuals: max violation / ||G|| and the count above 1e-4
  void residuals()
  {
    last_viol = 0; last_nviol = 0; last_nact = 0; last_rnorm = 0;
    for (size_t t = 0; t < cuts.size(); t++) {
      const double v = value(cuts[t]);
      if (v > 1e-4) last_nviol++;
      if (v > 0) { last_viol = std::max(last_viol, v / sqrt(cuts[t].g2)); last_rnorm += v * v / cuts[t].g2; }
      if (nu[t] > 0) last_nact++;
    }
    last_rnorm = sqrt(last_rnorm);
  }

  void run(const double *X, double *Y, const unsigned char *elim, int nsweeps)
  {
    double t0 = tcpu(), t1;
    rebuild_sums();
    pstep(X, Y, elim);
    t1 = tcpu(); t_p += t1 - t0; t0 = t1; n_p++;
    for (int k = 0; k < nsweeps; k++) {
      double d = 0;
      for (int p = 0; p < passes; p++) d = hpass();
      t1 = tcpu(); t_h += t1 - t0; t0 = t1; n_h += passes;
      last_dstep = d;
      pstep(X, Y, elim);
      t1 = tcpu(); t_p += t1 - t0; t0 = t1; n_p++;
    }
  }

  void project(const double *X, double *Y, double beta, const unsigned char *elim, int it) override
  {
    if (!init_mu.empty() && init_scale > 0) {
      for (size_t t = 0; t < nu.size() && t < init_mu.size(); t++) nu[t] = std::max(0.0, init_mu[t] / (init_scale * beta));
      init_mu.clear();
      beta_last = beta;
    }
    // the multipliers mu = scale beta nu stay put when beta changes
    if (beta_last > 0 && beta != beta_last) for (double &v : nu) v *= beta_last / beta;
    beta_last = beta;
    ysteps++;
    const bool diag = diag_every > 0 && it % diag_every == 0;
    const bool sep = sep_every > 0 && sep_k > 0 && it % sep_every == 0;
    std::vector<double> nu0, Psep, Wsep;
    if (diag) nu0 = nu;
    run(X, Y, elim, sweeps);
    residuals();
    if (sep) { Psep = Pag; Wsep = Wag; }   // the output Y's moments (the diagnostic overwrites them)
    if (diag) {
      const double td = tcpu(), tp0 = t_p, th0 = t_h; const int np0 = n_p, nh0 = n_h;
      // the distance of this Y-step's output to a longer run of the same projection
      // from the same start (nu0): ||Y - Y_long||_F, ||G'(nu - nu_long)|| proxy
      std::vector<double> nu_out = nu, Yl((size_t)n1 * n1);
      const double dstep_out = last_dstep;
      nu = nu0;
      run(X, Yl.data(), elim, diag_sweeps);
      last_dstep = dstep_out;
      double dy = 0, yn = 0;
      for (size_t t = 0; t < Yl.size(); t++) { dy += (Y[t] - Yl[t]) * (Y[t] - Yl[t]); yn += Yl[t] * Yl[t]; }
      double dn = 0;
      for (size_t t = 0; t < nu.size(); t++) dn += (nu_out[t] - nu[t]) * (nu_out[t] - nu[t]) * cuts[t].g2;
      double vmax = 0;
      for (size_t t = 0; t < cuts.size(); t++) vmax = std::max(vmax, value(cuts[t]) / sqrt(cuts[t].g2));
      diag_dy = sqrt(dy); diag_ynorm = sqrt(yn); diag_dnu = sqrt(dn); diag_viol = vmax; diag_it = it;
      nu = nu_out;
      rebuild_sums();   // the state of the regular run; Pag/Wag are refreshed by the next pstep
      t_p = tp0; t_h = th0; n_p = np0; n_h = nh0;   // the regular runs' counters exclude the diagnostic
      t_diag += tcpu() - td;
    }
    if (sep) separate(it, Psep, Wsep);
  }

  // keep the projector cuts with keep[t]; the others leave the pool and its keys
  void compact(const std::vector<char> &keep)
  {
    std::vector<CPCut> nc; std::vector<CPTerm> nt; std::vector<double> nn; std::vector<Cut> np; std::vector<int> nl;
    for (size_t t = 0; t < cuts.size(); t++) {
      const Cut &c = pool[cuts[t].src];
      if (!keep[t]) { if (c.type <= 3) keys.erase(cut_key(c)); n_dropped++; continue; }
      CPCut k = cuts[t];
      const int t0 = (int)nt.size();
      for (int q = k.t0; q < k.t1; q++) nt.push_back(terms[q]);
      k.t0 = t0; k.t1 = (int)nt.size(); k.src = (int)np.size();
      np.push_back(c); nc.push_back(k); nn.push_back(nu[t]); nl.push_back(last_pos[t]);
    }
    cuts.swap(nc); terms.swap(nt); nu.swap(nn); pool.swap(np); last_pos.swap(nl);
    rebuild_sums();
    gen++;   // the snapshots' cut indices are stale
  }

  // separation on the moments (Pm, Wm) of the Y-step's output at iteration it: the
  // regional triangles of dnn_cut_rounds (types 0-3 on i < j < k of every region,
  // violation above 1e-4, ranked by violation), the sep_k most violated ones not in the
  // pool appended with nu = 0
  void separate(int it, const std::vector<double> &Pm, const std::vector<double> &Wm)
  {
    const double t0 = tcpu();
    const size_t mm = (size_t)m * m;
    for (size_t t = 0; t < cuts.size(); t++) if (nu[t] > 0) last_pos[t] = it;
    std::vector<std::pair<double, Cut>> found;
    std::vector<double> p(m), W(mm);
    double vmax = 0;
    for (int s = 0; s < R; s++) {
      std::copy(Pm.begin() + (size_t)s * m, Pm.begin() + (size_t)(s + 1) * m, p.begin());
      std::copy(Wm.begin() + s * mm, Wm.begin() + (s + 1) * mm, W.begin());
      for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++)
          for (int k = j + 1; k < m; k++)
            for (int type = 0; type < 4; type++) {
              Cut t = {type, i, j, k, s, 0.0, 0, 0, 0, 0};
              const double v = rviolation(t, p, W, m);
              vmax = std::max(vmax, v);
              if (v > 1e-4) found.push_back({v, t});
            }
    }
    std::sort(found.begin(), found.end(), [](const std::pair<double, Cut> &u, const std::pair<double, Cut> &w) {
      return u.first > w.first; });
    // room: cuts whose nu has been 0 for drop_after Y-steps leave; under the cap the
    // zero-multiplier cuts leave too, the longest at zero first, as far as needed
    std::vector<char> keep(cuts.size(), 1);
    bool any = false;
    if (drop_after > 0)
      for (size_t t = 0; t < cuts.size(); t++)
        if (!(nu[t] > 0) && it - last_pos[t] >= drop_after) { keep[t] = 0; any = true; }
    int kept = 0;
    for (char k : keep) kept += k;
    if (pool_cap > 0 && kept + sep_k > pool_cap) {
      std::vector<std::pair<int, int>> zero;   // (last_pos, t) of the zero-multiplier cuts still kept
      for (size_t t = 0; t < cuts.size(); t++) if (keep[t] && !(nu[t] > 0)) zero.push_back({last_pos[t], (int)t});
      std::sort(zero.begin(), zero.end());
      for (size_t q = 0; q < zero.size() && kept + sep_k > pool_cap; q++) { keep[zero[q].second] = 0; kept--; any = true; }
    }
    if (any) compact(keep);
    const int room = pool_cap > 0 ? std::max(0, pool_cap - (int)cuts.size()) : sep_k;
    int added = 0;
    for (auto &f : found) {
      if (added >= std::min(sep_k, room)) break;
      if (keys.count(cut_key(f.second))) continue;
      pool.push_back(f.second);
      if (add_cut((int)pool.size() - 1)) { last_pos.back() = it; added++; }
      else pool.pop_back();
    }
    n_sep++; n_added += added; last_added = added; last_found = (int)found.size(); last_maxv = vmax;
    t_sep += tcpu() - t0;
  }

  // the dense cut term at mu = f nu (long double sums, one rounding per entry)
  // into L (added), its magnitudes and h'mu rounded up
  bool dense_term(long double f, const std::vector<double> &nuv, double *L, double &abs_add, double &hsum,
                  bool signed_ = false) const
  {
    const size_t mm = (size_t)m * m;
    std::vector<long double> cp((size_t)R * m, 0.0L), cw((size_t)R * mm, 0.0L);
    long double hs = 0, hab = 0, ab = 0;
    bool any = false;
    for (size_t t = 0; t < cuts.size(); t++) {
      if (signed_ ? nuv[t] == 0 : !(nuv[t] > 0)) continue;   // signed_: a test of the update formulas
      any = true;
      const CPCut &c = cuts[t];
      const long double mu = f * (long double)nuv[t];
      for (int q = c.t0; q < c.t1; q++) {
        const CPTerm &x = terms[q];
        if (x.j < 0) cp[(size_t)c.s * m + x.i] += mu * x.c;
        else cw[c.s * mm + (size_t)x.i * m + x.j] += mu * x.c;
      }
      hs += mu * c.h; hab += fabsl(mu * c.h); ab += fabsl(mu) * c.absrow;
    }
    abs_add = (double)ab * (1 + 1e-12);
    hsum = (double)hs;
    hsum = nextafter(hsum + 4.0 * (DBL_EPSILON / 2) * (double)hab + DBL_TRUE_MIN, HUGE_VAL);
    if (!any) { hsum = 0; abs_add = 0; return false; }
    // per sigma: 1/2 sum of the regions' W coefficients (symmetric), sum of the p ones
    const size_t ns = nsig();
    std::vector<long double> cs(mm), ds(m);
    std::vector<double> csd(mm), dsd(m);
    std::vector<std::vector<double>> C(ns), D(ns);
    for (size_t g = 0; g < ns; g++) {
      std::fill(cs.begin(), cs.end(), 0.0L); std::fill(ds.begin(), ds.end(), 0.0L);
      for (int s : sig_reg[g])
        for (int i = 0; i < m; i++) {
          ds[i] += cp[(size_t)s * m + i];
          for (int j = i + 1; j < m; j++) cs[(size_t)i * m + j] += 0.5L * cw[s * mm + (size_t)i * m + j];
        }
      C[g].assign(mm, 0.0); D[g].assign(m, 0.0);
      for (int i = 0; i < m; i++) {
        D[g][i] = (double)ds[i];
        for (int j = i + 1; j < m; j++) C[g][(size_t)i * m + j] = C[g][(size_t)j * m + i] = (double)cs[(size_t)i * m + j];
      }
    }
    if (side == 1) {
      // row 1 + i m + a, column 1 + j m + b (pstep_rows1): the same entries, one addition each
      const size_t sm = (size_t)m;
      for (int i = 0; i < m; i++)
        for (int a = 0; a < m; a++) {
          double *lr = L + (1 + (size_t)i * m + a) * n1;
          for (int b = 0; b < m; b++) {
            const int g = sig2[(size_t)a * m + b];
            if (g < 0) continue;
            double *lc = lr + 1 + b;   // column (j, b) at lc[j m]
            if (b == a) lc[i * sm] += D[g][i];
            else { const double *cg = &C[g][(size_t)i * m]; for (int j = 0; j < m; j++) if (j != i) lc[j * sm] += cg[j]; }
          }
        }
      return true;
    }
    for (int a = 0; a < m; a++)
      for (int i = 0; i < m; i++) {
        double *lr = L + (1 + (size_t)a * m + i) * n1;
        for (int b = 0; b < m; b++) {
          const int g = sig2[(size_t)a * m + b];
          if (g < 0) continue;
          double *lb = lr + 1 + (size_t)b * m;
          if (b == a) lb[i] += D[g][i];
          else { const double *cg = &C[g][(size_t)i * m]; for (int j = 0; j < m; j++) if (j != i) lb[j] += cg[j]; }
        }
      }
    return true;
  }

  bool cert_term(double scale, double beta, double *Ladd, double &abs_add, double &hsum) override
  {
    return dense_term((long double)scale * (long double)beta, nu, Ladd, abs_add, hsum);
  }

  // DNN_CP_LAZYCERT: a check's multipliers, for its certificate at the exit (the cuts
  // separated since enter at zero; a compaction since makes the slot unavailable)
  std::vector<double> snap_nu[2];
  long double snap_f[2] = {0, 0};
  int snap_gen[2] = {-1, -1}, gen = 0;
  void snapshot(int slot, double scale, double beta) override
  {
    if (slot < 0 || slot > 1) return;
    snap_nu[slot] = nu; snap_f[slot] = (long double)scale * (long double)beta; snap_gen[slot] = gen;
  }
  int cert_term_at(int slot, double *Ladd, double &abs_add, double &hsum) override
  {
    if (slot < 0 || slot > 1 || snap_gen[slot] != gen) return -1;
    std::vector<double> v = snap_nu[slot];
    v.resize(cuts.size(), 0.0);
    return dense_term(snap_f[slot], v, Ladd, abs_add, hsum) ? 1 : 0;
  }

  bool stop(int it, double best) override
  {
    if (stall_window <= 0 || !std::isfinite(stall_cutoff)) return false;
    const double b = std::max(best, stall_first);
    best_hist.push_back({it, b});
    if (it < stall_window) return false;
    double b0 = stall_first;   // the running best at it - W (first_lb before the first check)
    for (const auto &h : best_hist) { if (h.first > it - stall_window) break; b0 = h.second; }
    if (stall_cutoff - b > stall_factor * std::max(0.0, b - b0)) { stalled = true; return true; }
    return false;
  }

  std::string trace() const override
  {
    char buf[512];
    double sn = 0; for (double v : nu) sn += v;
    snprintf(buf, sizeof buf, "cuts %zu active %d sum_nu %.6g viol %d (max %.3e/||G||, rnorm %.3e) dstep %.3e",
             cuts.size(), last_nact, sn, last_nviol, last_viol, last_rnorm, last_dstep);
    std::string s = buf;
    snprintf(buf, sizeof buf, " | ms pstep %.2f (%d) hpass %.2f (%d) diag %.1f s",
             n_p ? 1e3 * t_p / n_p : 0.0, n_p, n_h ? 1e3 * t_h / n_h : 0.0, n_h, t_diag);
    s += buf;
    if (sep_every > 0) {
      snprintf(buf, sizeof buf, " | pool %zu sep %d (last: violated %d max %.4f added %d; total added %d dropped %d) sep %.1f s",
               cuts.size(), n_sep, last_found, last_maxv, last_added, n_added, n_dropped, t_sep);
      s += buf;
    }
    if (diag_it > 0) {
      snprintf(buf, sizeof buf, "  diag@%d (%d sweeps): |Y-Ylong| %.3e (|Ylong| %.3e) |G'(nu-nulong)| %.3e long maxviol %.3e",
               diag_it, diag_sweeps, diag_dy, diag_ynorm, diag_dnu, diag_viol);
      s += buf;
    }
    return s;
  }
};

}  // namespace

std::string dnn_pool_encode(const DNNPool &p)
{
  if (p.cuts.empty()) return "-";
  std::string s = "P1/" + std::to_string(p.side) + "/" + std::to_string(p.regions.size()) + "/";
  char buf[64];
  for (const auto &r : p.regions) {
    for (size_t k = 0; k < r.size(); k++) { if (k) s += ','; s += std::to_string(r[k]); }
    s += ';';
  }
  s += "/" + std::to_string(p.cuts.size()) + "/";
  for (const DNNPoolCut &c : p.cuts) {
    snprintf(buf, sizeof buf, "%.17g", c.mu);
    s += std::to_string(c.type) + "," + std::to_string(c.region) + "," + std::to_string(c.b) + "," + buf;
    for (int x : c.mem) s += "," + std::to_string(x);
    s += ';';
  }
  return s;
}

bool dnn_pool_decode(const std::string &s, DNNPool &p)
{
  p = DNNPool();
  if (s.empty() || s == "-") return true;
  auto split = [](const std::string &t, char d) {
    std::vector<std::string> v; size_t a = 0;
    while (a <= t.size()) { size_t b = t.find(d, a); if (b == std::string::npos) b = t.size(); v.push_back(t.substr(a, b - a)); a = b + 1; }
    return v; };
  const std::vector<std::string> f = split(s, '/');
  // P1 / side / regions / r0;r1;... / cuts / c0;c1;...
  if (f.size() != 6 || f[0] != "P1") return false;
  char *e;
  p.side = (int)strtol(f[1].c_str(), &e, 10);
  const long nr = strtol(f[2].c_str(), &e, 10), nc = strtol(f[4].c_str(), &e, 10);
  if ((p.side != 0 && p.side != 1) || nr < 0 || nc < 0) return false;
  std::vector<std::string> rs = split(f[3], ';'), cs = split(f[5], ';');
  if ((long)rs.size() != nr + 1 || (long)cs.size() != nc + 1) return false;
  for (long k = 0; k < nr; k++) {
    std::vector<int> r;
    for (const std::string &x : split(rs[k], ',')) { if (x.empty()) return false; r.push_back(atoi(x.c_str())); }
    p.regions.push_back(r);
  }
  for (long k = 0; k < nc; k++) {
    std::vector<std::string> t = split(cs[k], ',');
    if (t.size() < 4) return false;
    DNNPoolCut c;
    c.type = atoi(t[0].c_str()); c.region = atoi(t[1].c_str()); c.b = atoi(t[2].c_str());
    c.mu = strtod(t[3].c_str(), &e);
    if (!std::isfinite(c.mu) || c.mu < 0 || c.region < 0 || c.region >= nr) return false;
    for (size_t x = 4; x < t.size(); x++) c.mem.push_back(atoi(t[x].c_str()));
    p.cuts.push_back(c);
  }
  return true;
}

DNNCutResult dnn_cut_rounds(MAT *A, MAT *B, MAT *C, double shift, double cutoff,
                            const DNNParams &base, const std::vector<int> &free_orig,
                            int n_orig, const DNNCutParams &cp, DNNWarm *warm, double first_lb,
                            double first_primal)
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
  out.rounds = 0; out.iterations = 0; out.cpu = 0; out.active = 0; out.stalled = 0; out.gated = 0; out.cp_form = 0;
  out.nregions = (int)regions.size();
  std::vector<double> best_hist;   // best after each round (index = round)
  if (regional && regions.empty()) return out;

  std::vector<Cut> cuts;
  std::map<std::vector<int>, int> index;
  std::vector<L7Cut> l7;                // the pool's local-7 cuts (Cut::a indexes it)
  init_p7(); g_l7 = &l7;
  std::vector<BoxCut> boxes;            // the pool's box cuts (Cut::a indexes it)
  g_box = &boxes;
  // an inherited pool: its cuts at their multipliers, before any new cut
  const bool inherited = regional && cp.pool_in && !cp.pool_in->cuts.empty() && cp.other_orig && cp.fixed_at;
  if (inherited) {
    int dropped = 0;
    const int kept = remap_pool(*cp.pool_in, free_orig, n_orig, cp, regions, cuts, boxes, dropped);
    for (size_t k = 0; k < cuts.size(); k++) {
      const Cut &c = cuts[k];
      std::vector<int> key = {c.type, c.a, c.b, c.c, c.s, c.d, c.e, c.mask, c.alpha};
      if (c.type == 8) { key = {8, c.s, boxes[c.a].b}; key.insert(key.end(), boxes[c.a].Q.begin(), boxes[c.a].Q.end()); }
      index[key] = (int)k;
    }
    out.nregions = (int)regions.size();
    if (cp.verbose) printf("inherited pool: kept %d dropped %d, %d regions\n", kept, dropped, out.nregions);
  }
  std::vector<L7Space> l7space;         // per region, the admissible configurations
  if (cp.kcuts7 > 0 && m >= 7)
    for (auto &R : regions) l7space.push_back(l7_space((int)R.loc.size(), m));
  std::vector<double> Ladd((size_t)n1 * n1, 0.0);
  std::vector<std::vector<double>> P, WW;         // regional aggregates of the last Y
  DNNWarm save;
  double prev_lb = -HUGE_VAL, theta_cur = cp.theta;
  bool first_step = true;   // for cp.damp1

  double last_primal = HUGE_VAL, best_primal = -HUGE_VAL;   // best_primal: over the cut rounds
  std::vector<double> best_mu;   // the multipliers of the best round's solve (freeze_best)
  DNNWarm best_save;
  double best_mu_lb = -HUGE_VAL, best_mu_primal = HUGE_VAL;
  if (inherited && cp.first_iters <= 0 && warm && first_lb > -HUGE_VAL) {
    // the caller's solve with the inherited pool is round 0 and the best so far
    best_mu_lb = first_lb; best_mu_primal = first_primal;
    best_mu.resize(cuts.size());
    for (size_t t = 0; t < cuts.size(); t++) best_mu[t] = cuts[t].mu;
    best_save = *warm;
  }
  bool pool_checked = false;
  // one solve of the modified problem with the current multipliers
  auto solve = [&](int r, int iters) {
    if (cp.L_ext) std::copy(cp.L_ext, cp.L_ext + (size_t)n1 * n1, Ladd.begin());
    else std::fill(Ladd.begin(), Ladd.end(), 0.0);
    double abs_add = cp.L_ext_abs, hsum = cp.ext_const;
    const int active = assemble(cuts, regions, boxes, lift, n1, Ladd, abs_add, hsum);
    if (cp.pool_check > 0 && !cp.L_ext && !pool_checked) {
      const double w = check_term(Ladd, hsum - cp.ext_const, m, cp.pool_check, 12345u + r);
      if (w > 1e-7 * (1 + fabs(hsum))) {
        fprintf(stderr, "dnn_cut_rounds: the cut term is positive (%g) at a permutation\n", w); abort(); }
      pool_checked = true;
    }
    out.active = active;
    DNNParams par = base;
    par.warm_primal = cp.warm_primal && r > 0;
    par.max_iter = iters;
    par.hopeless_after = 0;
    par.hopeless_gap_ratio = 0;
    par.L_add = (active || cp.L_ext) ? Ladd.data() : NULL;
    par.L_add_abs = abs_add;
    par.cut_round = 1;   // DNN_BETA_SCHED's wc (dnn.C); a pooled child (dnn_pool_bound) leaves it 0
    DNNResult res = dnn_bound(A, B, C, shift - hsum, cutoff, par, (r > 0 && warm) ? warm : NULL, &save);
    if (warm) *warm = save; else warm = &save;
    out.cpu += res.cpu;
    if (r == 0) { out.first_lb = first_lb = res.lb; }
    else { out.rounds++; out.iterations += res.iterations; }
    out.best_lb = std::max(out.best_lb, res.lb);
    if (cp.freeze_best && res.lb > best_mu_lb) {
      best_mu_lb = res.lb;
      best_mu_primal = res.primal;
      best_mu.resize(cuts.size());
      for (size_t t = 0; t < cuts.size(); t++) best_mu[t] = cuts[t].mu;
      best_save = save;
    }
    out.status = res.status;
    last_primal = res.primal;
    if (r > 0 && std::isfinite(res.primal)) best_primal = std::max(best_primal, res.primal);
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
    std::vector<std::pair<double, Cut>> found, found5, found7, found8;
    std::vector<L7Cut> cand7;           // the local-7 candidates of this round (Cut::a indexes it)
    std::vector<BoxCut> cand8;          // the box candidates of this round (Cut::a indexes it)
    size_t nviol = 0, nviol5 = 0, nviol7 = 0, nviol8 = 0;
    auto by_violation = [](const std::pair<double, Cut> &u, const std::pair<double, Cut> &w) {
      return u.first > w.first; };
    auto trim = [&](std::vector<std::pair<double, Cut>> &v, int k) {
      if (v.size() > 4000000) {
        size_t keep = std::max<size_t>(100000, 8 * (size_t)k);
        std::nth_element(v.begin(), v.begin() + keep, v.end(), by_violation);
        v.resize(keep);
      }
    };
    auto bound_memory = [&]() { trim(found, cp.kcuts); };
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
        // the 5-point count inequalities: every 5-subset of the free elements
        // against the 48 patterns; the part over the first four elements is
        // shared by the innermost loop
        double r5 = -HUGE_VAL; size_t n5 = 0;
        if (cp.kcuts5 > 0 && m >= 5) {
          static const FivePatterns FP;
          const size_t np = FP.f.size();
          const std::vector<double> &p = P[s], &W = WW[s];
          auto w = [&](int i, int j) { return W[(size_t)i * m + j]; };
          std::vector<double> base(np), c4(np), c04(np), c14(np), c24(np), c34(np);
          for (size_t q = 0; q < np; q++) {
            c4[q] = FP.f[q].cp[4]; c04[q] = FP.f[q].cw[3]; c14[q] = FP.f[q].cw[6];
            c24[q] = FP.f[q].cw[8]; c34[q] = FP.f[q].cw[9];
          }
          int T[5];
          for (T[0] = 0; T[0] < m; T[0]++)
            for (T[1] = T[0] + 1; T[1] < m; T[1]++)
              for (T[2] = T[1] + 1; T[2] < m; T[2]++)
                for (T[3] = T[2] + 1; T[3] < m; T[3]++) {
                  for (size_t q = 0; q < np; q++) {
                    const Five &f = FP.f[q];
                    base[q] = f.c0 + f.cp[0] * p[T[0]] + f.cp[1] * p[T[1]] + f.cp[2] * p[T[2]] + f.cp[3] * p[T[3]]
                            + f.cw[0] * w(T[0], T[1]) + f.cw[1] * w(T[0], T[2]) + f.cw[2] * w(T[0], T[3])
                            + f.cw[4] * w(T[1], T[2]) + f.cw[5] * w(T[1], T[3]) + f.cw[7] * w(T[2], T[3]);
                  }
                  for (T[4] = T[3] + 1; T[4] < m; T[4]++) {
                    const double pe = p[T[4]], w04 = w(T[0], T[4]), w14 = w(T[1], T[4]),
                                 w24 = w(T[2], T[4]), w34 = w(T[3], T[4]);
                    for (size_t q = 0; q < np; q++) {
                      const double v = base[q] + c4[q] * pe + c04[q] * w04 + c14[q] * w14 + c24[q] * w24 + c34[q] * w34;
                      if (v > r5) r5 = v;
                      if (v > 1e-3) {
                        Cut t = {4, T[0], T[1], T[2], (int)s, 0.0, T[3], T[4], FP.pat[q].first, FP.pat[q].second};
                        found5.push_back({v, t}); n5++;
                      }
                    }
                  }
                  trim(found5, cp.kcuts5);
                }
          nviol5 += n5;
        }
        // the local-7 inequalities: cp.sub7 random 7-subsets of the free
        // elements plus a third as many from the 15 most fractional ones
        double r7 = -HUGE_VAL; size_t n7 = 0;
        if (cp.kcuts7 > 0 && m >= 7) {
          const std::vector<double> &p = P[s], &W = WW[s];
          std::mt19937 rng(1000003u * (unsigned)(r + 1) + (unsigned)s);
          std::vector<int> frac(m);
          for (int i = 0; i < m; i++) frac[i] = i;
          std::sort(frac.begin(), frac.end(), [&](int i, int j) { return p[i] * (1 - p[i]) > p[j] * (1 - p[j]); });
          const int nf = std::min(m, 15);
          const int nsub = cp.sub7 + cp.sub7 / 3;
          std::set<std::vector<int>> seen;
          for (int t = 0; t < nsub; t++) {
            std::vector<int> T;
            if (t < cp.sub7) { std::vector<int> all(m); for (int i = 0; i < m; i++) all[i] = i;
              std::shuffle(all.begin(), all.end(), rng); T.assign(all.begin(), all.begin() + 7); }
            else { std::vector<int> top(frac.begin(), frac.begin() + nf);
              std::shuffle(top.begin(), top.end(), rng); T.assign(top.begin(), top.begin() + 7); }
            std::sort(T.begin(), T.end());
            if (!seen.insert(T).second) continue;
            double y[28];
            for (int k = 0; k < 7; k++) y[k] = p[T[k]];
            for (int q = 0; q < 21; q++) y[7 + q] = W[(size_t)T[P7[q][0]] * m + T[P7[q][1]]];
            L7Cut c; double v;
            if (!l7_separate(l7space[s], y, c, v)) continue;
            for (int k = 0; k < 7; k++) c.f[k] = T[k];
            r7 = std::max(r7, v);
            if (v > 1e-3) {
              cand7.push_back(c);
              Cut ct = {6, (int)cand7.size() - 1, 0, 0, (int)s, 0.0};
              found7.push_back({v, ct}); n7++;
            }
          }
          nviol7 += n7;
        }
        // the box inequalities: local search over Q for b = 1..box_bmax from
        // the b + 2 largest p_i and from random weighted starts
        double r8 = -HUGE_VAL; size_t n8 = 0;
        if (cp.kbox > 0 && m >= 4) {
          const std::vector<double> &p = P[s], &W = WW[s];
          std::mt19937 rng(7919u * (unsigned)(r + 1) + 104729u * (unsigned)s + 17u);
          std::uniform_real_distribution<double> U(0.0, 1.0);
          const int qmax = std::min(cp.box_qmax, m - 1);
          std::set<std::vector<int>> seen;
          std::vector<int> order(m);
          for (int i = 0; i < m; i++) order[i] = i;
          std::sort(order.begin(), order.end(), [&](int i, int j) { return p[i] > p[j]; });
          for (int b = 1; b <= cp.box_bmax; b++) {
            if (b + 2 > qmax) break;
            for (int st = 0; st < cp.box_starts; st++) {
              std::vector<char> in(m, 0);
              if (st == 0) for (int i = 0; i < b + 2; i++) in[order[i]] = 1;
              else {
                const int hi = std::min(qmax, 3 * b + 6), k0 = b + 1 + (int)(U(rng) * (hi - b));
                // weighted sampling without replacement, weights p(1-p) + 1e-3
                std::vector<std::pair<double, int>> key(m);
                for (int i = 0; i < m; i++)
                  key[i] = {log(std::max(U(rng), 1e-300)) / (p[i] * (1 - p[i]) + 1e-3), i};
                std::partial_sort(key.begin(), key.begin() + std::min(k0, m), key.end(),
                                  [](const std::pair<double, int> &u, const std::pair<double, int> &w) { return u.first > w.first; });
                for (int i = 0; i < std::min(k0, m); i++) in[key[i].second] = 1;
              }
              box_local(p, W, m, b, qmax, in);
              BoxCut c; c.b = b;
              for (int i = 0; i < m; i++) if (in[i]) c.Q.push_back(i);
              std::vector<int> key = c.Q; key.push_back(-b);
              if (!seen.insert(key).second) continue;
              const double v = box_value(c, p, W, m);
              r8 = std::max(r8, v);
              if (v > 1e-3) {
                cand8.push_back(c);
                Cut ct = {8, (int)cand8.size() - 1, 0, 0, (int)s, 0.0};
                found8.push_back({v, ct}); n8++;
              }
            }
          }
          nviol8 += n8;
        }
        if (cp.verbose > 1) {
          printf("  region %-8s violated %zu  max %.4f", regions[s].name.c_str(), rn, rmax);
          if (cp.kcuts5 > 0) printf("  five %zu  max %.4f", n5, r5);
          if (cp.kcuts7 > 0) printf("  seven %zu  max %.4f", n7, r7);
          if (cp.kbox > 0) printf("  box %zu  max %.4f", n8, r8);
          printf("\n");
        }
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
    std::sort(found.begin(), found.end(), by_violation);
    std::sort(found5.begin(), found5.end(), by_violation);
    auto add_from = [&](std::vector<std::pair<double, Cut>> &v, int k) {
      int added = 0;
      for (auto &f : v) {
        if (added >= k) break;
        const Cut &c = f.second;
        std::vector<int> key = {c.type, c.a, c.b, c.c, c.s, c.d, c.e, c.mask, c.alpha};
        if (index.count(key)) continue;
        index[key] = (int)cuts.size();
        f.second.mu = cp.mu0;
        cuts.push_back(f.second);
        added++;
      }
      return added;
    };
    const bool adding = cp.add_until <= 0 || r <= cp.add_until;
    const int added = add_from(found, adding ? cp.kcuts : 0), added5 = add_from(found5, adding ? cp.kcuts5 : 0);
    // local-7 cuts: one per (region, 7-subset) in the pool; the candidate's data moves into l7
    std::sort(found7.begin(), found7.end(), by_violation);
    int added7 = 0;
    for (auto &f : found7) {
      if (added7 >= (adding ? cp.kcuts7 : 0)) break;
      const L7Cut &c = cand7[f.second.a];
      std::vector<int> key = {6, f.second.s, c.f[0], c.f[1], c.f[2], c.f[3], c.f[4], c.f[5], c.f[6]};
      if (index.count(key)) continue;
      index[key] = (int)cuts.size();
      l7.push_back(c);
      Cut ct = f.second; ct.a = (int)l7.size() - 1; ct.mu = cp.mu0;
      cuts.push_back(ct);
      added7++;
    }
    // box cuts: one per (region, b, Q) in the pool; the candidate's data moves into boxes
    std::sort(found8.begin(), found8.end(), by_violation);
    int added8 = 0;
    for (auto &f : found8) {
      if (added8 >= (adding ? cp.kbox : 0)) break;
      const BoxCut &c = cand8[f.second.a];
      std::vector<int> key = {8, f.second.s, c.b};
      key.insert(key.end(), c.Q.begin(), c.Q.end());
      if (index.count(key)) continue;
      index[key] = (int)cuts.size();
      boxes.push_back(c);
      Cut ct = f.second; ct.a = (int)boxes.size() - 1; ct.mu = cp.mu0;
      cuts.push_back(ct);
      added8++;
    }
    // the whole pool after this round's additions, zero multipliers included (CP-ADMM study)
    if (cp.pool_all_out && regional && cp.other_orig && r == cp.pool_all_round)
      export_pool(cuts, NULL, regions, boxes, free_orig, *cp.other_orig, cp.side, *cp.pool_all_out, true);
    const double maxv8 = found8.empty() ? 0 : found8[0].first;
    const double maxv = found.empty() ? 0 : found[0].first;
    const double maxv5 = found5.empty() ? 0 : found5[0].first;
    const double maxv7 = found7.empty() ? 0 : found7[0].first;
    if (cp.verbose && solved)
      printf("round %2d  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f"
             "  active %d  violated %zu (max %.4f)  added %d", r, res.lb, res.lb - first_lb,
             out.best_lb - first_lb, res.primal, dnn_status_name(res.status), res.iterations,
             res.cpu, out.active, nviol, maxv, added);
    else if (cp.verbose)
      printf("round  0  lb %.4f (the node's bound)  violated %zu (max %.4f)  added %d",
             first_lb, nviol, maxv, added);
    if (cp.verbose) {
      if (cp.kcuts5 > 0) printf("  five %zu (max %.4f) added5 %d", nviol5, maxv5, added5);
      if (cp.kcuts7 > 0) printf("  seven %zu (max %.4f) added7 %d", nviol7, maxv7, added7);
      if (cp.kbox > 0) printf("  box %zu (max %.4f) addedB %d", nviol8, maxv8, added8);
      printf("\n");
    }
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
    if (cp.pgate_rounds > 0 && r >= cp.pgate_rounds && std::isfinite(cutoff) &&
        best_primal < cutoff - cp.pgate_margin) {
      out.gated = 1;
      if (cp.verbose) printf("primal gate: best primal %.4f is %.4f below the cutoff after %d rounds\n",
                             best_primal, cutoff - best_primal, r);
      break;
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
    // the first step of an inherited phase moves multipliers that are already
    // near-optimal for this node: optionally damped (cp.damp1, later steps cp.damp_rest)
    const double theta_step = !inherited ? theta_cur : theta_cur * (first_step ? cp.damp1 : cp.damp_rest);
    first_step = false;
    const double alpha = theta_step * std::max(target - lb_now, 1e-6) / s2;
    for (size_t t = 0; t < cuts.size(); t++)
      cuts[t].mu = std::max(0.0, cuts[t].mu + alpha * sg[t]);
  }
  // the quiet round: the multipliers frozen, no new cuts, so that the
  // certificate catches up with the Lagrangian value of the final multipliers
  // an inherited phase tests the primal of the round the quiet round starts from
  const double quiet_primal = (inherited && cp.freeze_best && !best_mu.empty()) ? best_mu_primal : last_primal;
  bool quiet_ran = false;
  double quiet_lb = -HUGE_VAL;
  if (cp.freeze_iters > 0 && out.best_lb <= cutoff && !cuts.empty() && !out.gated &&
      (!cp.freeze_if_primal_above || quiet_primal > cutoff)) {
    quiet_ran = true;
    if (cp.freeze_best && !best_mu.empty() && best_mu_lb > -HUGE_VAL) {
      // back to the best round's multipliers and iterate (cuts added later keep mu 0)
      for (size_t t = 0; t < cuts.size(); t++) cuts[t].mu = t < best_mu.size() ? best_mu[t] : 0.0;
      *warm = best_save;
    }
    DNNResult res = solve(1, cp.freeze_iters);
    quiet_lb = res.lb;
    if (cp.verbose)
      printf("quiet round  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f  active %d\n",
             res.lb, res.lb - first_lb, out.best_lb - first_lb, res.primal,
             dnn_status_name(res.status), res.iterations, res.cpu, out.active);
  }
  if (out.best_lb > cutoff) out.status = DNN_FATHOMED;
  if (cp.pool_out && cp.other_orig && regional && out.status != DNN_FATHOMED) {
    // the children's pool: the quiet round's multipliers, else the best round's
    // (and its state), else the last solve's; a phase that did not inherit
    // hands its pool on only when that solve beat the uncut bound (else zero
    // multipliers are the better start and the children get none)
    *cp.pool_out = DNNPool();
    if (quiet_ran) {
      if (inherited || quiet_lb > first_lb)
        export_pool(cuts, NULL, regions, boxes, free_orig, *cp.other_orig, cp.side, *cp.pool_out);
    } else if (cp.freeze_best && !best_mu.empty() && best_mu_lb > -HUGE_VAL) {
      if (inherited || best_mu_lb > first_lb) {
        export_pool(cuts, &best_mu, regions, boxes, free_orig, *cp.other_orig, cp.side, *cp.pool_out);
        if (warm) *warm = best_save;
      }
    } else if (inherited)
      export_pool(cuts, NULL, regions, boxes, free_orig, *cp.other_orig, cp.side, *cp.pool_out);
  }
  return out;
}

DNNResult dnn_pool_bound(MAT *A, MAT *B, MAT *C, double shift, double cutoff,
                         const DNNParams &par, const std::vector<int> &free_orig, int n_orig,
                         const DNNCutParams &cp, const DNNWarm *start, DNNWarm *save,
                         int *kept, int *dropped)
{
  const int m = A->m, n1 = m * m + 1;
  const Lift lift = {m, cp.side};
  std::vector<Region> regions;
  std::vector<Cut> cuts;
  std::vector<BoxCut> boxes;
  g_box = &boxes;
  int nd = 0, nk = 0;
  if (cp.pool_in && cp.other_orig && cp.fixed_at)
    nk = remap_pool(*cp.pool_in, free_orig, n_orig, cp, regions, cuts, boxes, nd);
  if (kept) *kept = nk;
  if (dropped) *dropped = nd;
  std::vector<double> L;
  DNNParams p = par;
  double abs_add = 0, hsum = 0;
  if (nk > 0) {
    L.assign((size_t)n1 * n1, 0.0);
    assemble(cuts, regions, boxes, lift, n1, L, abs_add, hsum);
    if (cp.pool_check > 0) {
      const double w = check_term(L, hsum, m, cp.pool_check, 54321u);
      if (w > 1e-7 * (1 + fabs(hsum))) {
        fprintf(stderr, "dnn_pool_bound: the pool term is positive (%g) at a permutation\n", w); abort(); }
    }
    p.L_add = L.data();
    p.L_add_abs = abs_add;
  }
  DNNResult r = dnn_bound(A, B, C, shift - hsum, cutoff, p, start, save);
  if (cp.pool_out) export_pool(cuts, NULL, regions, boxes, free_orig, *cp.other_orig, cp.side, *cp.pool_out);
  return r;
}

namespace {
// The CP-ADMM self-test (DNNCPParams::selftest): the cut extraction against the
// code the cut rounds use. (a) the dense term of the projector at multipliers mu
// against assemble() with the same mu: L_add, the magnitude sum and h'mu; (b) for
// a sample of cuts, the explicit sparse row of the matrix add_rcut writes (mu = 1):
// ||G||^2 and <G, Y> against the projector's formulas, on Y from a P-step; (c) the
// P-step's moments against aggregate(); (d) the moments kept by a half-space pass
// against aggregate() on y_P - G'(nu_new - nu_old); (e) the projection's KKT
// conditions after up to 3,000 sweeps on a sample of the cuts, by independent code;
// (f) on side 1, the projector against side 0's on the transposed X, to the bit.
// All of it uses the Lift of the projector's side. Prints one line per check;
// returns false on a mismatch beyond rounding.
bool cp_selftest(CutProjector &pj, std::vector<Cut> cuts, const std::vector<Region> &regions,
                 const std::vector<BoxCut> &boxes, const Lift &lift, int n1)
{
  const int m = lift.m;
  const size_t nn = (size_t)n1 * n1;
  bool ok = true;
  const double tp_keep = pj.t_p, th_keep = pj.t_h;   // (e) and (f) run sweeps: the trace's timings exclude them
  const int np_keep = pj.n_p, nh_keep = pj.n_h;
  std::mt19937 rng(20261107u);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  // (a) random multipliers on about half the cuts
  std::vector<double> mu(pj.cuts.size(), 0.0);
  for (auto &v : mu) if (U(rng) < 0.5) v = U(rng) * 10;
  for (auto &c : cuts) c.mu = 0;
  for (size_t t = 0; t < pj.cuts.size(); t++) cuts[pj.cuts[t].src].mu = mu[t];
  std::vector<double> La(nn, 0.0), Lb(nn, 0.0);
  double absa = 0, hsa = 0, absb = 0, hsb = 0;
  pj.dense_term(1.0L, mu, La.data(), absa, hsa);
  assemble(cuts, regions, boxes, lift, n1, Lb, absb, hsb);
  double dmax = 0, lmax = 0; size_t nza = 0, nzb = 0;
  for (size_t t = 0; t < nn; t++) { dmax = std::max(dmax, fabs(La[t] - Lb[t])); lmax = std::max(lmax, fabs(Lb[t])); nza += La[t] != 0; nzb += Lb[t] != 0; }
  const bool oka = dmax <= 1e-12 * std::max(1.0, lmax) && fabs(absa - absb) <= 1e-9 * absb && fabs(hsa - hsb) <= 1e-9 * (1 + hsb);
  printf("cp-selftest (a) L_add vs assemble: max |diff| %.3e (max |entry| %.3e, nonzeros %zu vs %zu), abs %.10g vs %.10g, h'mu %.10g vs %.10g  %s\n",
         dmax, lmax, nza, nzb, absa, absb, hsa, hsb, oka ? "ok" : "MISMATCH");
  ok = ok && oka;
  // a P-step on a random symmetric X with entries in [-0.2, 1.2] and the multipliers above
  std::vector<double> X(nn), Y(nn);
  for (int r = 0; r < n1; r++) for (int c = r; c < n1; c++) X[(size_t)r * n1 + c] = X[(size_t)c * n1 + r] = -0.2 + 1.4 * U(rng);
  std::vector<double> nu_keep = pj.nu;
  for (size_t t = 0; t < mu.size(); t++) pj.nu[t] = mu[t] * 1e-3;
  pj.rebuild_sums();
  pj.pstep(X.data(), Y.data(), NULL);
  // (c) moments against aggregate()
  double dmom = 0, mmax = 0;
  std::vector<double> p, W;
  for (size_t s = 0; s < regions.size(); s++) {
    aggregate(regions[s], Y, lift, n1, p, W);
    for (int i = 0; i < m; i++) {
      dmom = std::max(dmom, fabs(p[i] - pj.Pag[s * m + i])); mmax = std::max(mmax, fabs(p[i]));
      for (int j = i + 1; j < m; j++) {
        const double w = W[(size_t)i * m + j];
        dmom = std::max(dmom, fabs(w - pj.Wag[s * m * m + (size_t)i * m + j])); mmax = std::max(mmax, fabs(w));
      }
    }
  }
  const bool okc = dmom <= 1e-10 * std::max(1.0, mmax);
  printf("cp-selftest (c) P-step moments vs aggregate(): max |diff| %.3e (max %.3e)  %s\n", dmom, mmax, okc ? "ok" : "MISMATCH");
  ok = ok && okc;
  // (b) explicit sparse rows for a sample of cuts
  std::vector<double> G(nn, 0.0);
  const int nsample = std::min<int>((int)pj.cuts.size(), 60);
  double dg2 = 0, dval = 0, dviol = 0;
  size_t nzmin = (size_t)-1, nzmax = 0;
  for (int q = 0; q < nsample; q++) {
    const size_t t = (size_t)(U(rng) * pj.cuts.size()) % pj.cuts.size();
    const CPCut &k = pj.cuts[t];
    const Cut &c = cuts[k.src];
    std::fill(G.begin(), G.end(), 0.0);
    if (c.type == 8) {
      std::vector<Cut> one(1, c); one[0].mu = 1;
      double a_ = 0, h_ = 0;
      assemble(one, regions, boxes, lift, n1, G, a_, h_);
    } else add_rcut(c, regions[c.s], G, lift, n1, 1.0);
    // the sparse row: (index, value) over the nonzeros of G
    std::vector<std::pair<size_t, double>> row;
    for (size_t e = 0; e < nn; e++) if (G[e] != 0) row.push_back({e, G[e]});
    double g2 = 0, gy = 0;
    for (auto &e : row) { g2 += e.second * e.second; gy += e.second * Y[e.first]; }
    nzmin = std::min(nzmin, row.size()); nzmax = std::max(nzmax, row.size());
    dg2 = std::max(dg2, fabs(g2 - k.g2) / k.g2);
    dval = std::max(dval, fabs(gy - k.h - pj.value(k)) / std::max(1.0, fabs(gy)));
    // and the separation's own violation on the moments
    aggregate(regions[c.s], Y, lift, n1, p, W);
    dviol = std::max(dviol, fabs(rviolation(c, p, W, m) - pj.value(k)) / std::max(1.0, fabs(gy)));
  }
  const bool okb = dg2 <= 1e-12 && dval <= 1e-10 && dviol <= 1e-10;
  printf("cp-selftest (b) %d sparse rows (%zu-%zu nonzeros): ||G||^2 rel diff %.3e, <G,Y> - h vs moments %.3e, vs rviolation %.3e  %s\n",
         nsample, nzmin, nzmax, dg2, dval, dviol, okb ? "ok" : "MISMATCH");
  ok = ok && okb;
  // (d) one half-space pass: the kept moments against aggregate() on y_P - G'(dnu)
  std::vector<double> nu0 = pj.nu;
  pj.hpass();
  std::vector<double> dnu(pj.nu.size());
  for (size_t t = 0; t < dnu.size(); t++) dnu[t] = pj.nu[t] - nu0[t];
  std::vector<double> T(nn, 0.0);
  double a_ = 0, h_ = 0;
  pj.dense_term(1.0L, dnu, T.data(), a_, h_, true);
  for (size_t e = 0; e < nn; e++) Y[e] -= T[e];
  double dh = 0; mmax = 0;
  for (size_t s = 0; s < regions.size(); s++) {
    aggregate(regions[s], Y, lift, n1, p, W);
    for (int i = 0; i < m; i++) {
      dh = std::max(dh, fabs(p[i] - pj.Pag[s * m + i])); mmax = std::max(mmax, fabs(p[i]));
      for (int j = i + 1; j < m; j++) {
        const double w = W[(size_t)i * m + j];
        dh = std::max(dh, fabs(w - pj.Wag[s * m * m + (size_t)i * m + j])); mmax = std::max(mmax, fabs(w));
      }
    }
  }
  size_t moved = 0; for (double d : dnu) moved += d != 0;
  const bool okd = dh <= 1e-9 * std::max(1.0, mmax);
  printf("cp-selftest (d) a half-space pass (%zu multipliers moved): kept moments vs aggregate(y_P - G'dnu) max |diff| %.3e (max %.3e)  %s\n",
         moved, dh, mmax, okd ? "ok" : "MISMATCH");
  ok = ok && okd;
  // (e) the projection's optimality conditions, by code independent of the projector's
  // grouping, on a projector of its own with a random sample of at most 500 of the
  // pool's cuts (the same code; Dykstra from a random X with the whole pool of 13,000
  // is at 2.4e-3 after 1,000 sweeps), from nu = 0 on X, after 1, 10, 100, 1,000 and, if
  // needed, 3,000 sweeps: Y in P (Y_00 = 1, the box, the gangster zeros from the lifted
  // indices decoded), Y = P(X - G'nu) with G'nu from assemble() at mu = nu, and the cut
  // values from aggregate() and rviolation(): feasibility max (<G,Y> - h)+ / ||G|| and
  // complementarity max over nu > 0 of |<G,Y> - h| / ||G||. Y in P and Y = P(X - G'nu)
  // hold at every sweep (the sweep ends with a P-step); with nu >= 0 the two residuals
  // at 0 are the KKT conditions of Y = argmin ||Y - X|| over P and the half-spaces
  {
    std::vector<size_t> pick(pj.cuts.size());
    for (size_t t = 0; t < pick.size(); t++) pick[t] = t;
    if (pick.size() > 500) { std::shuffle(pick.begin(), pick.end(), rng); pick.resize(500); std::sort(pick.begin(), pick.end()); }
    std::vector<Cut> sub;
    for (size_t t : pick) { sub.push_back(cuts[pj.cuts[t].src]); sub.back().mu = 0; }
    CutProjector pe(m, regions, sub, pj.side);
    pe.passes = pj.passes;
    std::vector<double> Ye(nn), Lg(nn);
    std::vector<std::vector<double>> pr(regions.size()), Wr(regions.size());
    const int stages[5] = {1, 10, 100, 1000, 3000};
    int done = 0;
    double feas = 0, comp = 0, dP = 0, dfix = 0;
    size_t nact = 0;
    for (int st = 0; st < 5; st++) {
      if (done >= 1000 && feas <= 1e-6 && comp <= 1e-6) break;
      pe.run(X.data(), Ye.data(), NULL, stages[st] - done);
      done = stages[st];
      auto gangster = [&](int r, int c) {
        if (r == 0 || c == 0) return false;
        const int f1 = (r - 1) / m, s1 = (r - 1) % m, f2 = (c - 1) / m, s2 = (c - 1) % m;
        return (f1 == f2) != (s1 == s2);
      };
      // Y in P
      dP = fabs(Ye[0] - 1.0);
      for (int r = 0; r < n1; r++)
        for (int c = 0; c < n1; c++) {
          const double v = Ye[(size_t)r * n1 + c];
          if (gangster(r, c)) dP = std::max(dP, fabs(v));
          else if (r + c > 0) dP = std::max(dP, std::max(0.0, std::max(-v, v - 1.0)));
        }
      // Y = P(X - G'nu), G'nu by assemble()
      std::vector<Cut> cn = sub;
      for (size_t t = 0; t < pe.cuts.size(); t++) cn[pe.cuts[t].src].mu = pe.nu[t];
      std::fill(Lg.begin(), Lg.end(), 0.0);
      double a_ = 0, h_ = 0;
      assemble(cn, regions, boxes, lift, n1, Lg, a_, h_);
      dfix = 0;
      for (int r = 0; r < n1; r++)
        for (int c = 0; c < n1; c++) {
          if (r + c == 0) continue;
          const size_t e = (size_t)r * n1 + c;
          const double ref = gangster(r, c) ? 0.0 : std::min(1.0, std::max(0.0, X[e] - Lg[e]));
          dfix = std::max(dfix, fabs(Ye[e] - ref));
        }
      // feasibility and complementarity on aggregate()'s moments
      for (size_t s = 0; s < regions.size(); s++) aggregate(regions[s], Ye, lift, n1, pr[s], Wr[s]);
      feas = 0; comp = 0; nact = 0;
      for (size_t t = 0; t < pe.cuts.size(); t++) {
        const Cut &c = sub[pe.cuts[t].src];
        const double v = rviolation(c, pr[c.s], Wr[c.s], m) / sqrt(pe.cuts[t].g2);
        feas = std::max(feas, v);
        if (pe.nu[t] > 0) { comp = std::max(comp, fabs(v)); nact++; }
      }
      printf("cp-selftest (e) %4d sweeps: Y in P %.1e, Y - P(X - G'nu) %.1e, max violation/||G|| %.3e, complementarity %.3e (%zu active of %zu)\n",
             done, dP, dfix, feas, comp, nact, pe.cuts.size());
    }
    const bool oke = dP == 0 && dfix <= 1e-12 && feas <= 1e-6 && comp <= 1e-6;
    printf("cp-selftest (e) the projection's KKT conditions after %d sweeps (%zu of the %zu cuts)  %s\n", done,
           pe.cuts.size(), pj.cuts.size(), oke ? "ok" : "MISMATCH");
    ok = ok && oke;
  }
  // (f) side 1: the projector against side 0's on the transposed problem (lifted index
  // 1 + i m + a <-> 1 + a m + i): the P-step's Y, the moments, a half-space pass, three
  // sweeps and the dense term, all equal to the bit
  if (pj.side == 1) {
    CutProjector p0(m, regions, cuts, 0);
    std::vector<size_t> perm(n1);
    perm[0] = 0;
    for (int i = 0; i < m; i++) for (int a = 0; a < m; a++) perm[1 + (size_t)i * m + a] = 1 + (size_t)a * m + i;
    std::vector<double> Xt(nn), Y1(nn), Y0(nn);
    for (int r = 0; r < n1; r++) for (int c = 0; c < n1; c++) Xt[perm[r] * n1 + perm[c]] = X[(size_t)r * n1 + c];
    auto ydiff = [&]() { double d = 0; for (int r = 0; r < n1; r++) for (int c = 0; c < n1; c++)
                           d = std::max(d, fabs(Y1[(size_t)r * n1 + c] - Y0[perm[r] * n1 + perm[c]])); return d; };
    auto mdiff = [&]() { double d = 0; for (size_t q = 0; q < pj.Pag.size(); q++) d = std::max(d, fabs(pj.Pag[q] - p0.Pag[q]));
                           for (size_t q = 0; q < pj.Wag.size(); q++) d = std::max(d, fabs(pj.Wag[q] - p0.Wag[q])); return d; };
    auto ndiff = [&]() { double d = 0; for (size_t q = 0; q < pj.nu.size(); q++) d = std::max(d, fabs(pj.nu[q] - p0.nu[q])); return d; };
    bool same = p0.cuts.size() == pj.cuts.size() && p0.nsig() == pj.nsig();
    double dy = -1, dm = -1, dn = -1, dy3 = -1, dn3 = -1, dl = -1;
    if (same) {
      for (size_t t = 0; t < mu.size(); t++) pj.nu[t] = p0.nu[t] = mu[t] * 1e-3;
      pj.rebuild_sums(); p0.rebuild_sums();
      pj.pstep(X.data(), Y1.data(), NULL); p0.pstep(Xt.data(), Y0.data(), NULL);
      dy = ydiff(); dm = mdiff();
      pj.hpass(); p0.hpass();
      dn = ndiff();
      pj.run(X.data(), Y1.data(), NULL, 3); p0.run(Xt.data(), Y0.data(), NULL, 3);
      dy3 = ydiff(); dn3 = ndiff();
      std::vector<double> L1(nn, 0.0), L0(nn, 0.0);
      double a1 = 0, h1 = 0, a0 = 0, h0 = 0;
      pj.dense_term(1.0L, mu, L1.data(), a1, h1); p0.dense_term(1.0L, mu, L0.data(), a0, h0);
      dl = fabs(a1 - a0) + fabs(h1 - h0);
      for (int r = 0; r < n1; r++) for (int c = 0; c < n1; c++) dl = std::max(dl, fabs(L1[(size_t)r * n1 + c] - L0[perm[r] * n1 + perm[c]]));
    }
    const bool okf = same && dy == 0 && dm == 0 && dn == 0 && dy3 == 0 && dn3 == 0 && dl == 0;
    printf("cp-selftest (f) side 1 vs side 0 on the transposed X (%zu cuts, %zu classes): P-step Y %.1e, moments %.1e, "
           "half-space pass nu %.1e, 3 sweeps Y %.1e nu %.1e, dense term %.1e  %s\n",
           pj.cuts.size(), pj.nsig(), dy, dm, dn, dy3, dn3, dl, okf ? "ok" : "MISMATCH");
    ok = ok && okf;
  }
  pj.nu = nu_keep;
  pj.rebuild_sums();
  pj.t_p = tp_keep; pj.t_h = th_keep; pj.n_p = np_keep; pj.n_h = nh_keep;
  return ok;
}
}  // namespace

DNNResult dnn_cp_bound(MAT *A, MAT *B, MAT *C, double shift, double cutoff, const DNNParams &par,
                       const std::vector<int> &free_orig, int n_orig, const DNNCutParams &cp,
                       DNNCPParams &cpp, const DNNWarm *start, DNNWarm *save, int *kept, int *dropped)
{
  DNNResult bad;
  bad.lb = -HUGE_VAL; bad.primal = HUGE_VAL; bad.status = DNN_ERROR; bad.iterations = 0; bad.cpu = 0;
  bad.certificates = 0; bad.cert_failures = 0; bad.refinements = 0; bad.halpern_restarts = 0; bad.margin = 0;
  const int m = A->m, n1 = m * m + 1;
  if (cp.atomic || (cp.side != 0 && cp.side != 1)) { fprintf(stderr, "dnn_cp_bound: regional cuts only (side 0 or 1)\n"); return bad; }
  const Lift lift = {m, cp.side};
  std::vector<Region> regions;
  if (cp.grid_rows > 0) regions = grid_regions(cp.grid_rows, cp.grid_cols, free_orig, n_orig);
  else if (cp.prefix) regions = prefix_regions(m);
  std::vector<Cut> cuts;
  std::vector<BoxCut> boxes;
  std::vector<L7Cut> l7;
  init_p7(); g_box = &boxes; g_l7 = &l7;
  int nd = 0, nk = 0;
  if (cp.pool_in && cp.other_orig && cp.fixed_at)
    nk = remap_pool(*cp.pool_in, free_orig, n_orig, cp, regions, cuts, boxes, nd, !cpp.child);
  if (kept) *kept = nk;
  if (dropped) *dropped = nd;
  if (cpp.child && nk == 0) {   // a pooled child whose pool did not survive: the plain solve
    if (cp.pool_out) *cp.pool_out = DNNPool();
    return dnn_bound(A, B, C, shift, cutoff, par, start, save);
  }
  CutProjector pj(m, regions, cuts, cp.side);
  pj.sweeps = std::max(1, cpp.sweeps); pj.passes = std::max(1, cpp.passes);
  pj.diag_every = cpp.diag_every; pj.diag_sweeps = std::max(1, cpp.diag_sweeps);
  if (cpp.omega > 0 && cpp.omega < 2) pj.omega = cpp.omega;
  pj.sep_every = std::max(0, cpp.sep_every); pj.sep_k = std::max(0, cpp.sep_k);
  pj.pool_cap = std::max(0, cpp.pool_cap); pj.drop_after = std::max(0, cpp.drop_after);
  if (cpp.nu_from_pool) {
    // start from the pool's multipliers: nu = mu / (scale beta), scale = max(1, max |L|)
    // of the lifted objective as dnn_bound forms it (no L_add here)
    double mx = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        mx = std::max(mx, fabs(0.5 * C->me[i][j]));
        for (int k = 0; k < m; k++)
          for (int l = 0; l < m; l++)
            mx = std::max(mx, fabs(0.5 * (A->me[i][k] * B->me[j][l] + A->me[k][i] * B->me[l][j])));
      }
    pj.init_scale = std::max(1.0, mx);
    pj.init_mu.assign(pj.cuts.size(), 0.0);
    for (size_t t = 0; t < pj.cuts.size(); t++) pj.init_mu[t] = cuts[pj.cuts[t].src].mu;
  }
  if (cp.verbose) {
    size_t nt = 0; for (auto &c : pj.cuts) nt += c.t1 - c.t0;
    printf("cp pool: %d cuts kept (%d dropped), %zu in the projector, %zu terms, %d regions, %zu region-pair classes; sweeps %d passes %d omega %g\n",
           nk, nd, pj.cuts.size(), nt, (int)regions.size(), pj.nsig(), pj.sweeps, pj.passes, pj.omega);
  }
  if (cpp.selftest && !cp_selftest(pj, cuts, regions, boxes, lift, n1)) {
    fprintf(stderr, "dnn_cp_bound: self-test failed\n");
    return bad;
  }
  DNNParams p = par;
  p.ycuts = &pj;
  p.L_add = NULL; p.L_add_abs = 0;
  DNNWarm own;
  DNNWarm *sv = save ? save : &own;
  DNNResult r = dnn_bound(A, B, C, shift, cutoff, p, start, sv);
  if (cpp.child && r.status == DNN_ERROR && r.iterations == 0)   // a setting the CP solve refuses
    return dnn_pool_bound(A, B, C, shift, cutoff, par, free_orig, n_orig, cp, start, save, kept, dropped);
  if (cpp.child && cp.pool_out) {
    // the children's pool: the cuts with nu > 0 at the final mu_D = scale beta nu
    *cp.pool_out = DNNPool();
    if (r.status != DNN_FATHOMED && sv->scale > 0 && cp.other_orig) {
      const long double f = (long double)sv->scale * (long double)pj.beta_last;
      std::vector<Cut> ex;
      for (size_t t = 0; t < pj.cuts.size(); t++) {
        Cut c = pj.pool[pj.cuts[t].src];
        c.mu = (double)(f * (long double)pj.nu[t]);
        if (c.mu > 0) ex.push_back(c);
      }
      if (cp.pool_check > 0 && !ex.empty()) {
        std::vector<double> L((size_t)n1 * n1, 0.0);
        double abs_add = 0, hsum = 0;
        assemble(ex, regions, boxes, lift, n1, L, abs_add, hsum);
        const double w = check_term(L, hsum, m, cp.pool_check, 13579u);
        if (w > 1e-7 * (1 + fabs(hsum))) {
          fprintf(stderr, "dnn_cp_bound: the pool term is positive (%g) at a permutation\n", w); abort(); }
      }
      if (!ex.empty()) export_pool(ex, NULL, regions, boxes, free_orig, *cp.other_orig, cp.side, *cp.pool_out);
    }
  }
  cpp.pool_final = (int)pj.cuts.size(); cpp.sep_added = pj.n_added; cpp.sep_dropped = pj.n_dropped;
  if (cp.verbose && pj.sep_every > 0)
    printf("cp separation: %d separations, %d cuts added, %d dropped, pool %zu at the end, %.1f s\n", pj.n_sep, pj.n_added,
           pj.n_dropped, pj.cuts.size(), pj.t_sep);
  if (cpp.selftest && pj.sep_every > 0) {
    // the extraction self-test again, on the pool the separations grew
    std::vector<double> nu_end = pj.nu;
    if (!cp_selftest(pj, pj.pool, regions, boxes, lift, n1)) fprintf(stderr, "dnn_cp_bound: self-test after the separations failed\n");
    pj.nu = nu_end; pj.rebuild_sums();
  }
  cpp.quiet_lb = -HUGE_VAL;
  if (cpp.quiet_iters > 0 && r.status != DNN_FATHOMED && r.status != DNN_ERROR && sv->scale > 0) {
    // the quiet round: mu_D = scale beta nu frozen as L_add, the cut rounds' path
    std::vector<double> La((size_t)n1 * n1, 0.0);
    double abs_add = 0, hsum = 0;
    const bool any = pj.dense_term((long double)sv->scale * (long double)pj.beta_last, pj.nu, La.data(), abs_add, hsum);
    DNNParams q = par;
    q.ycuts = NULL;
    q.max_iter = cpp.quiet_iters;
    q.L_add = any ? La.data() : NULL; q.L_add_abs = any ? abs_add : 0;
    if (!any) hsum = 0;
    DNNWarm st = *sv;
    DNNResult rq = dnn_bound(A, B, C, shift - hsum, cutoff, q, &st, sv);
    cpp.quiet_lb = rq.lb;
    if (cp.verbose)
      printf("cp quiet round  lb %.10f  primal %.4f  status %s  iters %d  cpu %.1f  h'mu %.6g\n", rq.lb, rq.primal,
             dnn_status_name(rq.status), rq.iterations, rq.cpu, hsum);
  }
  return r;   // the CP solve's own result; the quiet round's bound is in cpp.quiet_lb
}

// max(1, max |L|) of the lifted objective as dnn_bound forms it, without L_add:
// the scale of a CP solve, for nu = mu / (scale beta)
static double cp_lifted_scale(const MAT *A, const MAT *B, const MAT *C, int m)
{
  double mx = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) {
      mx = std::max(mx, fabs(0.5 * C->me[i][j]));
      for (int k = 0; k < m; k++)
        for (int l = 0; l < m; l++)
          mx = std::max(mx, fabs(0.5 * (A->me[i][k] * B->me[j][l] + A->me[k][i] * B->me[l][j])));
    }
  return std::max(1.0, mx);
}

DNNCutResult dnn_cut_cp(MAT *A, MAT *B, MAT *C, double shift, double cutoff,
                        const DNNParams &base, const std::vector<int> &free_orig,
                        int n_orig, const DNNCutParams &cp, DNNWarm *warm,
                        double first_lb, double first_primal)
{
  const int m = A->m, n1 = m * m + 1;
  const size_t nn = (size_t)n1 * n1, mm = (size_t)m * m;
  // what the CP form does not cover runs as rounds
  if (cp.atomic || (cp.side != 0 && cp.side != 1) || !warm || warm->m != m || warm->Y.size() != nn || warm->Z.size() != nn ||
      !(warm->scale > 0))
    return dnn_cut_rounds(A, B, C, shift, cutoff, base, free_orig, n_orig, cp, warm, first_lb, first_primal);
  const double t0 = CutProjector::tcpu();
  const Lift lift = {m, cp.side};
  std::vector<Region> regions;
  if (cp.grid_rows > 0) regions = grid_regions(cp.grid_rows, cp.grid_cols, free_orig, n_orig);
  else if (cp.prefix) regions = prefix_regions(m);
  DNNCutResult out;
  out.first_lb = first_lb; out.best_lb = first_lb; out.status = -1;
  out.rounds = 0; out.iterations = 0; out.cpu = 0; out.active = 0; out.stalled = 0; out.gated = 0; out.cp_form = 0;
  out.nregions = (int)regions.size();
  out.cp_form = 1;
  if (regions.empty()) return out;
  std::vector<Cut> cuts;
  std::vector<BoxCut> boxes;
  std::vector<L7Cut> l7;
  init_p7(); g_box = &boxes; g_l7 = &l7;
  // the inherited pool at its multipliers, then round 0's separation on the first solve's Y
  const bool inherited = cp.pool_in && !cp.pool_in->cuts.empty() && cp.other_orig && cp.fixed_at;
  if (inherited) {
    int dropped = 0;
    remap_pool(*cp.pool_in, free_orig, n_orig, cp, regions, cuts, boxes, dropped);
  }
  out.nregions = (int)regions.size();
  CutProjector pj(m, regions, cuts, cp.side);
  pj.sweeps = std::max(1, cp.cp_sweeps); pj.passes = std::max(1, cp.cp_passes);
  {
    const int R = (int)regions.size();
    std::vector<double> Pm((size_t)R * m, 0.0), Wm((size_t)R * mm, 0.0), p, W;
    for (int s = 0; s < R; s++) {
      aggregate(regions[s], warm->Y, lift, n1, p, W);
      std::copy(p.begin(), p.end(), Pm.begin() + (size_t)s * m);
      std::copy(W.begin(), W.end(), Wm.begin() + s * mm);
    }
    pj.sep_k = std::max(0, cp.kcuts);
    if (pj.sep_k > 0) pj.separate(0, Pm, Wm);
  }
  if (pj.cuts.empty()) {
    // no cut at all: the rounds stop here too (no solve, status -1, no pool)
    if (cp.pool_out && cp.other_orig) *cp.pool_out = DNNPool();
    out.cpu = CutProjector::tcpu() - t0;
    return out;
  }
  if (inherited) {
    // nu from the inherited multipliers (the separated cuts start at zero)
    pj.init_scale = cp_lifted_scale(A, B, C, m);
    pj.init_mu.assign(pj.cuts.size(), 0.0);
    for (size_t t = 0; t < pj.cuts.size(); t++) pj.init_mu[t] = std::max(0.0, pj.pool[pj.cuts[t].src].mu);
  }
  pj.sep_every = std::max(0, cp.cp_every); pj.sep_k = std::max(0, cp.cp_k);
  pj.stall_window = cp.stall_rounds > 0 ? cp.stall_rounds * cp.iters : 0;
  pj.stall_factor = cp.stall_factor; pj.stall_cutoff = cutoff; pj.stall_first = first_lb;
  DNNParams par = base;
  par.ycuts = &pj;
  par.L_add = NULL; par.L_add_abs = 0;
  par.max_iter = cp.cp_iters > 0 ? cp.cp_iters : std::max(1, cp.rounds * cp.iters + std::max(0, cp.freeze_iters));
  par.hopeless_after = 0; par.hopeless_gap_ratio = 0;   // as the cut rounds
  par.warm_primal = cp.warm_primal;
  par.cut_round = 1;   // DNN_BETA_SCHED's wc, as the rounds
  DNNWarm save;
  DNNResult res = dnn_bound(A, B, C, shift, cutoff, par, warm, &save);
  if (res.status == DNN_ERROR && res.iterations == 0)   // a setting the CP solve refuses (the accelerations)
    return dnn_cut_rounds(A, B, C, shift, cutoff, base, free_orig, n_orig, cp, warm, first_lb, first_primal);
  out.iterations = res.iterations;
  out.rounds = pj.n_sep;
  out.best_lb = std::max(first_lb, res.lb);
  out.stalled = pj.stalled ? 1 : 0;
  out.status = out.best_lb > cutoff ? DNN_FATHOMED : res.status;
  for (double v : pj.nu) out.active += v > 0;
  if (cp.pool_out && cp.other_orig && out.status != DNN_FATHOMED) {
    // the children's pool: the cuts with nu > 0 at the final mu_D = scale beta nu, handed
    // on with the final state; a phase that did not inherit hands on a pool only when
    // its solve beat the uncut bound (as the rounds)
    *cp.pool_out = DNNPool();
    if ((inherited || res.lb > first_lb) && save.scale > 0) {
      const long double f = (long double)save.scale * (long double)pj.beta_last;
      std::vector<Cut> ex;
      for (size_t t = 0; t < pj.cuts.size(); t++) {
        Cut c = pj.pool[pj.cuts[t].src];
        c.mu = (double)(f * (long double)pj.nu[t]);
        if (c.mu > 0) ex.push_back(c);
      }
      if (cp.pool_check > 0 && !ex.empty()) {
        std::vector<double> L(nn, 0.0);
        double abs_add = 0, hsum = 0;
        assemble(ex, regions, boxes, lift, n1, L, abs_add, hsum);
        const double w = check_term(L, hsum, m, cp.pool_check, 24680u);
        if (w > 1e-7 * (1 + fabs(hsum))) {
          fprintf(stderr, "dnn_cut_cp: the pool term is positive (%g) at a permutation\n", w); abort(); }
      }
      export_pool(ex, NULL, regions, boxes, free_orig, *cp.other_orig, cp.side, *cp.pool_out);
    }
  }
  *warm = save;
  out.cpu = CutProjector::tcpu() - t0;
  if (cp.verbose)
    printf("cp phase  lb %.4f (first %.4f)  status %s  iters %d  separations %d  pool %zu  active %d  stalled %d  cpu %.1f\n",
           out.best_lb, first_lb, dnn_status_name(out.status), out.iterations, out.rounds, pj.cuts.size(), out.active,
           out.stalled, out.cpu);
  return out;
}
