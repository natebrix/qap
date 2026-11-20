// t0.C: the C landscape signal (latest/lsla.C of the experiment build) at one node, standalone.
//   t0 DAT SLN FIXES SLACK KICKS[:KICKS3] KICKS2 XHI SEED [LINES]
// KICKS3 > 0: a random-start search per cell (KICKS3 kicks) and the harvest (lsla.h).
// FIXES r:v,r:v,... (p[r] = v, 0-based); LINES e.g. c54,c47 (default: none). Prints the node's own excess,
// each free line's predicted survivors (ex <= SLACK, after base kicks and after refinement), every child of
// LINES with ex/SLACK, and the CPU. Build: g++ -O3 -march=native -std=gnu++17 -I$BUILD/latest t0.C $BUILD/latest/lsla.C
#include "lsla.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

int main(int argc, char **argv) {
  if (argc < 9) { fprintf(stderr, "usage: t0 DAT SLN FIXES SLACK KICKS KICKS2 XHI SEED [LINES]\n"); return 2; }
  FILE *f = fopen(argv[1], "r");
  LSInst I;
  if (fscanf(f, "%d", &I.n) != 1) return 1;
  const int n = I.n;
  I.a.resize(n * n); I.b.resize(n * n);
  for (auto &x : I.a) if (fscanf(f, "%lf", &x) != 1) return 1;
  for (auto &x : I.b) if (fscanf(f, "%lf", &x) != 1) return 1;
  fclose(f);
  f = fopen(argv[2], "r");
  int nn; double obj;
  if (fscanf(f, "%d %lf", &nn, &obj) != 2 || nn != n) return 1;
  std::vector<int> p(n);
  for (auto &x : p) { if (fscanf(f, "%d", &x) != 1) return 1; x--; }
  fclose(f);
  std::vector<int> fixed(n, -1);
  for (char *t = strtok(argv[3], ","); t; t = strtok(NULL, ",")) { int r, v; sscanf(t, "%d:%d", &r, &v); fixed[r] = v; }
  const double s = atof(argv[4]);
  LSParams P; P.kicks = atoi(argv[5]); if (const char *k3 = strchr(argv[5], ':')) P.kicks3 = atoi(k3 + 1); P.kicks2 = atoi(argv[6]); P.xhi = atof(argv[7]);
  const unsigned seed = atoi(argv[8]);
  LSTable T;
  ls_node(I, p.data(), fixed.data(), s, P, false, T, seed);
  printf("n %d m %d z %.0f (sln %.0f) ex_start %.0f ex_node %.0f s %.0f kicks %d kicks2 %d xhi %.2f seed %u\n", n, T.m, T.z, obj,
	 T.ex_start, T.ex_node, s, P.kicks, P.kicks2, P.xhi, seed);
  printf("cpu gate %.2f table %.2f refine %.2f (%ld cells refined) moves %ld\n", T.cpu_gate, T.cpu_table, T.cpu_refine,
	 T.refined, T.moves);
  if (P.kicks3 > 0) printf("kicks3 %d: random searches %.2f s (in table), %ld cells lowered by the harvest\n", P.kicks3, T.cpu_random, T.harvested);
  // every line: predicted survivors after base and after refinement
  struct L { std::string name; int cb, c; double off; };
  std::vector<L> Ls;
  LSTable B = T; B.ex = T.ex_base;
  for (int side = 0; side < 2; side++)
    for (int k = 0; k < T.m; k++) {
      int cb, c; double ob, o;
      ls_line_score(B, side == 0, k, s, 0.01 * s, cb, ob);
      ls_line_score(T, side == 0, k, s, 0.01 * s, c, o);
      Ls.push_back({(side == 0 ? "r" : "c") + std::to_string(side == 0 ? T.frow[k] : T.fcol[k]), cb, c, o});
    }
  std::sort(Ls.begin(), Ls.end(), [](const L &x, const L &y) { return x.c != y.c ? x.c < y.c : x.off < y.off; });
  printf("lines (name base_cnt cnt off), by cnt then off:");
  for (auto &l : Ls) printf(" %s:%d:%d:%.2f", l.name.c_str(), l.cb, l.c, l.off);
  printf("\n");
  if (argc > 9)
    for (char *t = strtok(argv[9], ","); t; t = strtok(NULL, ",")) {
      const bool row = t[0] == 'r';
      const int idx = atoi(t + 1);
      const std::vector<int> &v = row ? T.frow : T.fcol;
      const int k = std::find(v.begin(), v.end(), idx) - v.begin();
      std::vector<std::pair<double, std::string>> ch;
      for (int u = 0; u < T.m; u++) {
	const int a = row ? k : u, b = row ? u : k;
	ch.push_back({T.ex[a * T.m + b], std::to_string(T.frow[a]) + ":" + std::to_string(T.fcol[b]) + " base_x " +
				       std::to_string(T.ex_base[a * T.m + b] / s)});
      }
      std::sort(ch.begin(), ch.end());
      printf("line %s:", t);
      for (int u = 0; u < std::min(12, (int)ch.size()); u++) printf("  %s x %.2f", ch[u].second.c_str(), ch[u].first / s);
      printf("\n");
    }
  return 0;
}
