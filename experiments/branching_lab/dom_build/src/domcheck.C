// domcheck: the solver's own pair exchange test (latest/dominance.C) on every
// node of a node trace, as bnbdrv sets it up (readQAPLIB, set_fathom_tolerance,
// symmetrize_QAP). Prints per node: the fixes, 4 x the best pair's min Delta
// (an exact integer), the pair, and whether the node is dominated.
//   domcheck instance.dat trace[.gz via zcat] > out
#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>
#include "qap.h"
#include "assign.h"
#include "dominance.h"
int main(int argc, char **argv)
{
  QAP *q = readQAPLIB(argv[1]);
  set_fathom_tolerance(q);
  symmetrize_QAP(q);
  PairDominance pd;
  std::string why;
  if (!pd.init(q, why)) { std::cerr << "unavailable: " << why << std::endl; return 2; }
  const int n = q->A->m;
  char line[1 << 16];
  FILE *in = argc > 2 ? fopen(argv[2], "r") : stdin;
  while (fgets(line, sizeof line, in)) {
    if (strncmp(line, "node ", 5)) continue;
    char fixes[1 << 15];
    int depth;
    if (sscanf(line, "node %d %s", &depth, fixes) != 2) continue;
    QAPAssignment a(n);
    if (strcmp(fixes, "-")) {
      char *s = fixes;
      while (*s) {
        int r, c, k;
        if (sscanf(s, "%d:%d%n", &r, &c, &k) != 2) break;
        a.fix(r, c);
        s += k;
        if (*s == ',') s++;
      }
    }
    int pi, pk; long long b4;
    bool dom = pd.dominated(a, &pi, &pk, &b4, true);
    if (pi < 0) printf("%s - - 0\n", fixes);
    else printf("%s %lld %d,%d %d\n", fixes, b4, pi, pk, dom ? 1 : 0);
  }
  return 0;
}
