// brute-force optimum (and number of optimal permutations) of a QAPLIB file, n <= 11
// cost(p) = sum_uv A[u][v] B[p(u)][p(v)] (+ C[u][p(u)] if present), as bnbdrv's QAPobjective
#include <cstdio>
#include <vector>
#include <algorithm>
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); int n; if (fscanf(f, "%d", &n) != 1) return 1;
  std::vector<long long> A(n*n), B(n*n), C(n*n, 0);
  for (auto &x : A) { double v; fscanf(f, "%lf", &v); x = (long long)v; }
  for (auto &x : B) { double v; fscanf(f, "%lf", &v); x = (long long)v; }
  { double v; int k = 0; while (k < n*n && fscanf(f, "%lf", &v) == 1) C[k++] = (long long)v; }
  std::vector<int> p(n); for (int i = 0; i < n; i++) p[i] = i;
  long long best = -1, cnt = 0;
  do {
    long long c = 0;
    for (int u = 0; u < n; u++) { c += C[u*n+p[u]]; for (int v = 0; v < n; v++) c += A[u*n+v]*B[p[u]*n+p[v]]; }
    if (best < 0 || c < best) { best = c; cnt = 1; } else if (c == best) cnt++;
  } while (std::next_permutation(p.begin(), p.end()));
  printf("%lld %lld\n", best, cnt);
}
