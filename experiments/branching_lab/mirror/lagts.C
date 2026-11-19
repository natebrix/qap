// lagts: robust tabu search for min over the completions of a node of
//   f(p) = cost(p) + sum_j mu_j * Delta_j(p),   Delta_j(p) = cost(p) - cost(p with rows i_j, k_j exchanged)
// cost(p) = sum_{u,v} A[u][v] B[p(u)][p(v)], A the file's first matrix; a fix r:c means p(r) = c.
// Any completion's f is >= min_p f(p) >= any relaxation's bound on f, so the best f found is a CEILING
// on what a Lagrangian bound with these multipliers can certify at this node.
// usage: lagts FILE.dat FIXES ("r:c,r:c,...") PAIRS ("i:k:mu,i:k:mu,..." or "-") ITERS RESTARTS SEED [MODE 1=SA (default), 0=tabu] [T0 factor] [T1/T0]
// prints: best f, its cost, each Delta_j, and the permutation (0-based); and every restart's local best.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>
typedef double D;
int n;
std::vector<D> A, B;
#define a(u, v) A[(size_t)(u) * n + (v)]
#define b(u, v) B[(size_t)(u) * n + (v)]
D cost(const std::vector<int> &p) {
  D s = 0;
  for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) s += a(u, v) * b(p[u], p[v]);
  return s;
}
int main(int argc, char **argv) {
  if (argc < 7) { fprintf(stderr, "usage\n"); return 1; }
  FILE *f = fopen(argv[1], "r");
  if (fscanf(f, "%d", &n) != 1) return 1;
  A.resize(n * n); B.resize(n * n);
  for (auto &x : A) { long v; if (fscanf(f, "%ld", &v) != 1) return 1; x = v; }
  for (auto &x : B) { long v; if (fscanf(f, "%ld", &v) != 1) return 1; x = v; }
  fclose(f);
  std::vector<int> fixc(n, -1);
  { std::string s = argv[2]; if (s != "-") { char *t = strtok(&s[0], ","); while (t) { int r, c; sscanf(t, "%d:%d", &r, &c); fixc[r] = c; t = strtok(NULL, ","); } } }
  struct Pair { int i, k; D mu; };
  std::vector<Pair> pr;
  { std::string s = argv[3]; if (s != "-") { char *t = strtok(&s[0], ","); while (t) { Pair q; sscanf(t, "%d:%d:%lf", &q.i, &q.k, &q.mu); pr.push_back(q); t = strtok(NULL, ","); } } }
  long iters = atol(argv[4]); int restarts = atoi(argv[5]); unsigned seed = atoi(argv[6]);
  int mode = argc > 7 ? atoi(argv[7]) : 1; D t0f = argc > 8 ? atof(argv[8]) : 0.3, t1f = argc > 9 ? atof(argv[9]) : 1e-4;
  std::vector<int> fr, fc; std::vector<bool> used(n, false);
  for (int r = 0; r < n; r++) if (fixc[r] < 0) fr.push_back(r); else used[fixc[r]] = true;
  for (int c = 0; c < n; c++) if (!used[c]) fc.push_back(c);
  const int m = fr.size();
  // linear term L[l][c] over free rows l (index into fr) and columns c (original), and constant K
  std::vector<D> L((size_t)m * n, 0.0); D K = 0;
  for (int u = 0; u < n; u++) if (fixc[u] >= 0) for (int v = 0; v < n; v++) if (fixc[v] >= 0) K += a(u, v) * b(fixc[u], fixc[v]);
  for (int x = 0; x < m; x++) { int l = fr[x];
    for (int c : fc) { D s = 0;
      for (int u = 0; u < n; u++) if (fixc[u] >= 0) s += a(l, u) * b(c, fixc[u]) + a(u, l) * b(fixc[u], c);
      s += a(l, l) * b(c, c);  // the diagonal term of a free row
      L[(size_t)x * n + c] = s; } }
  for (auto &q : pr) {
    int i = q.i, k = q.k, aa = fixc[i], bb = fixc[k];
    if (aa < 0 || bb < 0) { fprintf(stderr, "pair %d:%d not both fixed\n", i, k); return 1; }
    D d0 = (a(i, i) - a(k, k)) * (b(aa, aa) - b(bb, bb)) + (a(i, k) - a(k, i)) * (b(aa, bb) - b(bb, aa));
    for (int l = 0; l < n; l++) if (fixc[l] >= 0 && l != i && l != k) { int c = fixc[l];
      d0 += (a(i, l) - a(k, l)) * (b(aa, c) - b(bb, c)) + (a(l, i) - a(l, k)) * (b(c, aa) - b(c, bb)); }
    K += q.mu * d0;
    for (int x = 0; x < m; x++) { int l = fr[x];
      for (int c : fc) L[(size_t)x * n + c] += q.mu * ((a(i, l) - a(k, l)) * (b(aa, c) - b(bb, c)) + (a(l, i) - a(l, k)) * (b(c, aa) - b(c, bb))); }
  }
  // the reduced problem over free rows: Q(x) = sum_{x,y free} A[fr x][fr y] B[c_x][c_y] (x != y) + L + K
  std::vector<D> Af((size_t)m * m);
  for (int x = 0; x < m; x++) for (int y = 0; y < m; y++) Af[(size_t)x * m + y] = (x == y) ? 0 : a(fr[x], fr[y]);
  auto full = [&](const std::vector<int> &q) {  // q[x] = column of free row x
    D s = K; for (int x = 0; x < m; x++) { s += L[(size_t)x * n + q[x]]; for (int y = 0; y < m; y++) s += Af[(size_t)x * m + y] * b(q[x], q[y]); } return s; };
  auto swapd = [&](const std::vector<int> &q, int r, int s) {
    int pr_ = q[r], ps = q[s]; D d = L[(size_t)r * n + ps] + L[(size_t)s * n + pr_] - L[(size_t)r * n + pr_] - L[(size_t)s * n + ps];
    d += Af[(size_t)r * m + s] * (b(ps, pr_) - b(pr_, ps)) + Af[(size_t)s * m + r] * (b(pr_, ps) - b(ps, pr_));
    for (int k2 = 0; k2 < m; k2++) if (k2 != r && k2 != s) { int pk = q[k2];
      d += (Af[(size_t)k2 * m + r] - Af[(size_t)k2 * m + s]) * (b(pk, ps) - b(pk, pr_)) + (Af[(size_t)r * m + k2] - Af[(size_t)s * m + k2]) * (b(ps, pk) - b(pr_, pk)); }
    return d; };
  std::mt19937 rng(seed);
  std::vector<int> best; D bestf = 1e300;
  std::vector<D> dl((size_t)m * m);
  std::vector<long> tabu((size_t)m * n);
  for (int rs = 0; rs < restarts; rs++) {
    std::vector<int> q(fc.begin(), fc.end()); std::shuffle(q.begin(), q.end(), rng);
    D fq = full(q); std::vector<int> lb = q; D lbf = fq;
    if (mode == 1) {  // simulated annealing, geometric cooling, then the best point found
      D sum = 0; for (int t = 0; t < 2000; t++) { int r = rng() % m, s2 = rng() % m; if (r == s2) continue; sum += std::abs(swapd(q, r, s2)); }
      D T0 = t0f * sum / 2000, T1 = T0 * t1f, T = T0, g = std::pow(T1 / T0, 1.0 / iters);
      std::uniform_real_distribution<D> U(0, 1);
      for (long it = 0; it < iters; it++, T *= g) {
        int r = rng() % m, s2 = rng() % m; if (r == s2) continue;
        D d = swapd(q, r, s2);
        if (d <= 0 || U(rng) < std::exp(-d / T)) { std::swap(q[r], q[s2]); fq += d; if (fq < lbf - 1e-9) { lbf = fq; lb = q; } }
      }
      q = lb; fq = lbf;
      for (bool imp = true; imp;) { imp = false;  // descent to a 2-exchange local optimum
        for (int r = 0; r < m; r++) for (int s2 = r + 1; s2 < m; s2++) { D d = swapd(q, r, s2); if (d < -1e-9) { std::swap(q[r], q[s2]); fq += d; imp = true; } } }
      lb = q; lbf = fq;
    } else {
    std::fill(tabu.begin(), tabu.end(), -1);
    int tmin = std::max(1, (int)(0.9 * m)), tmax = std::max(2, (int)(1.1 * m));
    for (int r = 0; r < m; r++) for (int s = r + 1; s < m; s++) dl[(size_t)r * m + s] = swapd(q, r, s);
    for (long it = 0; it < iters; it++) {
      int br = -1, bs = -1; D bd = 1e300; int cnt = 0;
      for (int r = 0; r < m; r++) for (int s = r + 1; s < m; s++) {
        D d = dl[(size_t)r * m + s];
        bool tb = tabu[(size_t)r * n + q[s]] >= it && tabu[(size_t)s * n + q[r]] >= it;
        bool asp = fq + d < lbf - 1e-9;
        if (tb && !asp) continue;
        if (d < bd - 1e-9) { bd = d; br = r; bs = s; cnt = 1; }
        else if (d < bd + 1e-9 && (rng() % (++cnt)) == 0) { br = r; bs = s; }
      }
      if (br < 0) continue;
      tabu[(size_t)br * n + q[br]] = it + tmin + rng() % (tmax - tmin + 1);
      tabu[(size_t)bs * n + q[bs]] = it + tmin + rng() % (tmax - tmin + 1);
      std::swap(q[br], q[bs]); fq += bd;
      if (fq < lbf - 1e-9) { lbf = fq; lb = q; }
      for (int r = 0; r < m; r++) for (int s = r + 1; s < m; s++) dl[(size_t)r * m + s] = swapd(q, r, s);
    }
    }
    D drift = fq - full(q);
    lbf = full(lb);
    // report this restart's best as a full permutation, with its cost and Deltas
    std::vector<int> p(n); for (int r = 0; r < n; r++) p[r] = fixc[r];
    for (int x = 0; x < m; x++) p[fr[x]] = lb[x];
    D c0 = cost(p); printf("sol f %.1f cost %.1f", lbf, c0);
    D chk = c0;
    for (auto &qq : pr) { std::vector<int> p2 = p; std::swap(p2[qq.i], p2[qq.k]); D dd = c0 - cost(p2); chk += qq.mu * dd; printf(" delta %.1f", dd); }
    printf(" chk %.1f drift %.1f p", chk - lbf, drift);
    for (int r = 0; r < n; r++) printf("%c%d", r ? ',' : ' ', p[r]);
    printf("\n");
    if (lbf < bestf) { bestf = lbf; best = p; }
  }
  printf("best f %.1f\n", bestf);
  return 0;
}
