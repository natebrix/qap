// exch.C: stage 2, the exact ceiling of the exchange-optimality cuts on small instances.
//   g++ -O2 -o exch exch.C
//   ./exch file.dat [seed=1] [lines=3] [maxdepth=4] [capfrac=0.05]
// Finds a global optimum p* by enumeration (cheap bound), then for depth d = 1..maxdepth and
// `lines` random facility orders: the node prefix = the first d-1 facilities of the order at
// their p* locations (kind G, the good path), or the same with the last prefix facility moved
// to the wrong location of smallest node optimum (kind W, the near-miss sibling); the line =
// the next facility at every free location. For each child:
//   v_node = min cost over its completions,
//   v_ex1  = min over completions that no exchange of a FIXED facility with a FREE one improves,
//   v_ex2  = the same plus exchanges of two fixed facilities,
//   v_2opt = min over completions that no exchange at all improves (not linear in Y; reference).
// Values above cap = v_node + capfrac * opt are reported as "inf" (none below the cap).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <random>
typedef long long ll;
const int N = 16;
int n; ll A[N][N], B[N][N];
int p[N]; bool usedl[N];
int order_[N]; int nfix;
ll lin[N + 1][N][N];
ll best0, best1, best2, best3, thr_cap; int mode; // mode 0: plain minimum; 1: the exchange minima
int bestp[N];
long long leaves;

static ll cost(const int *q) {
  ll s = 0;
  for (int k = 0; k < n; k++) for (int l = 0; l < n; l++) s += A[k][l] * B[q[k]][q[l]];
  return s;
}
// C(q with i, j exchanged) - C(q)
static ll delta(const int *q, int i, int j) {
  int a = q[i], b = q[j];
  ll s = (A[i][j] - A[j][i]) * (B[b][a] - B[a][b]) + (A[i][i] - A[j][j]) * (B[b][b] - B[a][a]);
  for (int k = 0; k < n; k++) {
    if (k == i || k == j) continue;
    int c = q[k];
    s += (A[i][k] - A[j][k]) * (B[b][c] - B[a][c]) + (A[k][i] - A[k][j]) * (B[c][b] - B[c][a]);
  }
  return s;
}
static void leaf(ll c) {
  leaves++;
  if (mode == 0) { if (c < best0) { best0 = c; memcpy(bestp, p, sizeof(int) * n); } return; }
  if (c < best0) best0 = c;
  // fixed facilities are order_[0..nfix-1]
  if (c >= best3) return;
  bool ok1 = true;
  for (int s = 0; s < nfix && ok1; s++)
    for (int t = nfix; t < n; t++)
      if (delta(p, order_[s], order_[t]) < 0) { ok1 = false; break; }
  if (!ok1) return;
  if (c < best1) best1 = c;
  bool ok2 = true;
  for (int s = 0; s < nfix && ok2; s++)
    for (int t = s + 1; t < nfix; t++)
      if (delta(p, order_[s], order_[t]) < 0) { ok2 = false; break; }
  if (!ok2) return;
  if (c < best2) best2 = c;
  for (int s = nfix; s < n; s++)
    for (int t = s + 1; t < n; t++)
      if (delta(p, order_[s], order_[t]) < 0) return;
  best3 = c;
}
static void dfs(int t, ll partial) {
  if (t == n) { leaf(partial); return; }
  int f = order_[t];
  // candidate locations by increasing linear cost
  int cand[N], nc = 0;
  for (int c = 0; c < n; c++) if (!usedl[c]) cand[nc++] = c;
  std::sort(cand, cand + nc, [&](int x, int y) { return lin[t][f][x] < lin[t][f][y]; });
  for (int ci = 0; ci < nc; ci++) {
    int c = cand[ci];
    ll nc_ = partial + lin[t][f][c];
    ll thr = mode == 0 ? best0 : std::min(thr_cap, best3);
    if (nc_ >= thr) break;   // sorted: the rest are no better (data nonnegative)
    usedl[c] = true; p[f] = c;
    ll bound = nc_;
    for (int u = t + 1; u < n && bound < thr; u++) {
      int k = order_[u];
      ll mn = -1;
      for (int c2 = 0; c2 < n; c2++) {
        if (usedl[c2]) continue;
        ll v = lin[t][k][c2] + A[k][f] * B[c2][c] + A[f][k] * B[c][c2];
        lin[t + 1][k][c2] = v;
        if (mn < 0 || v < mn) mn = v;
      }
      bound += mn;
    }
    if (bound < thr) dfs(t + 1, nc_);
    usedl[c] = false;
  }
}
// solve the node with fixed pairs (fac[s] at loc[s]), free facilities in the order given
static void solve(const std::vector<int> &fac, const std::vector<int> &loc, const std::vector<int> &freeorder, int m_, ll cap) {
  mode = m_; thr_cap = cap;
  memset(usedl, 0, sizeof(usedl));
  nfix = fac.size();
  for (int s = 0; s < nfix; s++) { order_[s] = fac[s]; p[fac[s]] = loc[s]; usedl[loc[s]] = true; }
  for (int s = 0; s < (int)freeorder.size(); s++) order_[nfix + s] = freeorder[s];
  ll partial = 0;
  for (int s = 0; s < nfix; s++) for (int t = 0; t < nfix; t++) partial += A[fac[s]][fac[t]] * B[loc[s]][loc[t]];
  for (int u = nfix; u < n; u++) {
    int k = order_[u];
    for (int c = 0; c < n; c++) {
      ll v = A[k][k] * B[c][c];
      for (int s = 0; s < nfix; s++) v += A[k][fac[s]] * B[c][loc[s]] + A[fac[s]][k] * B[loc[s]][c];
      lin[nfix][k][c] = v;
    }
  }
  leaves = 0;
  dfs(nfix, partial);
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r");
  if (!f || fscanf(f, "%d", &n) != 1 || n > N) { fprintf(stderr, "bad file\n"); return 1; }
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (fscanf(f, "%lld", &A[i][j]) != 1) return 1;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (fscanf(f, "%lld", &B[i][j]) != 1) return 1;
  fclose(f);
  int seed = argc > 2 ? atoi(argv[2]) : 1, lines = argc > 3 ? atoi(argv[3]) : 3, maxd = argc > 4 ? atoi(argv[4]) : 4;
  double capfrac = argc > 5 ? atof(argv[5]) : 0.05;
  std::mt19937 rng(seed);
  // facility order for the search: decreasing total flow
  std::vector<int> byflow(n);
  for (int i = 0; i < n; i++) byflow[i] = i;
  std::vector<ll> tot(n, 0);
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) tot[i] += A[i][j] + A[j][i];
  std::sort(byflow.begin(), byflow.end(), [&](int x, int y) { return tot[x] > tot[y]; });
  // incumbent by 2-opt restarts, then the exact root solve
  int q[N]; ll ub = -1; int ubp[N];
  for (int r = 0; r < 200; r++) {
    for (int i = 0; i < n; i++) q[i] = i;
    std::shuffle(q, q + n, rng);
    bool imp = true;
    while (imp) { imp = false; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (delta(q, i, j) < 0) { std::swap(q[i], q[j]); imp = true; } }
    ll c = cost(q);
    if (ub < 0 || c < ub) { ub = c; memcpy(ubp, q, sizeof(int) * n); }
  }
  best0 = ub + 1; memcpy(bestp, ubp, sizeof(int) * n);
  solve({}, {}, byflow, 0, 0);
  ll opt = best0; int popt[N]; memcpy(popt, bestp, sizeof(int) * n);
  printf("# %s n=%d optimum %lld (root leaves %lld) p* =", argv[1], n, opt, leaves);
  for (int i = 0; i < n; i++) printf(" %d", popt[i]);
  printf("\n# cap = v_node + %.3f * opt; kind G = good-path prefix, W = near-miss sibling prefix\n", capfrac);
  printf("# kind depth line prefix fac loc v_node v_ex1 v_ex2 v_2opt onpath\n");
  // exchange-optimality of p* itself (a global optimum must be)
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (delta(popt, i, j) < 0) { printf("BUG: p* not exchange-optimal\n"); return 1; }
  ll capadd = (ll)(capfrac * opt);
  for (int ln = 0; ln < lines; ln++) {
    std::vector<int> ord(n);
    for (int i = 0; i < n; i++) ord[i] = i;
    std::shuffle(ord.begin(), ord.end(), rng);
    for (int kind = 0; kind < 2; kind++) {
      for (int d = (kind ? 2 : 1); d <= maxd; d++) {
        std::vector<int> fac(ord.begin(), ord.begin() + d - 1), loc;
        for (int s = 0; s < d - 1; s++) loc.push_back(popt[fac[s]]);
        if (kind == 1) {
          // move the last prefix facility to the wrong location with the smallest node optimum
          int fl = fac[d - 2]; ll bv = -1; int bl = -1;
          std::vector<bool> used(n, false);
          for (int s = 0; s < d - 2; s++) used[loc[s]] = true;
          for (int c = 0; c < n; c++) {
            if (used[c] || c == popt[fl]) continue;
            std::vector<int> l2(loc); l2[d - 2] = c;
            std::vector<int> fo; for (int x : byflow) if (std::find(fac.begin(), fac.end(), x) == fac.end()) fo.push_back(x);
            best0 = (ll)4e18; solve(fac, l2, fo, 0, 0);
            if (bv < 0 || best0 < bv) { bv = best0; bl = c; }
          }
          loc[d - 2] = bl;
        }
        int j = ord[d - 1];
        std::vector<bool> used(n, false);
        for (int s = 0; s < d - 1; s++) used[loc[s]] = true;
        for (int b = 0; b < n; b++) {
          if (used[b]) continue;
          std::vector<int> f2(fac), l2(loc); f2.push_back(j); l2.push_back(b);
          std::vector<int> fo; for (int x : byflow) if (std::find(f2.begin(), f2.end(), x) == f2.end()) fo.push_back(x);
          best0 = (ll)4e18; solve(f2, l2, fo, 0, 0);
          ll vnode = best0;
          ll cap = vnode + capadd + 1;
          best0 = best1 = best2 = best3 = (ll)4e18;
          solve(f2, l2, fo, 1, cap);
          printf("%c %d %d ", kind ? 'W' : 'G', d, ln);
          for (int s = 0; s < d - 1; s++) printf("%s%d:%d", s ? "," : "", fac[s], loc[s]);
          if (d == 1) printf("-");
          printf(" %d %d %lld ", j, b, vnode);
          ll v[3] = {best1, best2, best3};
          for (int t = 0; t < 3; t++) { if (v[t] > (ll)3e18) printf("inf "); else printf("%lld ", v[t]); }
          bool onp = kind == 0 && b == popt[j];
          printf("%d\n", onp ? 1 : 0);
          fflush(stdout);
        }
      }
    }
  }
  return 0;
}
