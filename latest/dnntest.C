/*
============================================================
 dnntest: run the DNN bound (dnn.h) on a QAPLIB instance or on one
 node of it, and print the certified bound.

   ./dnntest [options] file.dat

   -x i:j,i:j,...  fix facility i to location j (0-based) before bounding
   -K i:j          after the node, bound its child fixing i to j, cold and
                   warm-started from the node's ADMM state
   -N file         bound many nodes: one per line, "<i:j,... or -> [cutoff]"
   -u X            incumbent; the cutoff is X - (fathoming tolerance), the
                   tolerance being the solver's (1.95 for even objectives)
   -c X            cutoff given directly (overrides -u)
   -i N            iteration limit (4000)
   -t X            CPU time limit in seconds (none)
   -k N            certify every N iterations (10)
   -b              box polyhedral part (as dnn.py) instead of Gilmore-Lawler
   -l N            inner LAPs of the GL part: 0 reduce, 1 exact, 2 auto (default: exact
                   only within 1e-3 |cutoff| below the cutoff)
   -P              prototype mode: dnn.py's certificate (extra eigendecomposition,
                   box part), beta factor 1, check every 50, stall window 6, no
                   hopeless test (each overridable after -P)
   -s N            stall window in checks (0 = off)
   -H N            give up (status hopeless) after N iterations when the certified
                   bound, extrapolated from its recent gains, stays below the
                   cutoff (default 200, 0 = off)
   -R rho          with -H: give up instead when the gap to the cutoff shrank by
                   less than the factor rho over the last 50 iterations
   -E i:j,...      eliminated assignments (x_ij = 0; free at the node): the bound is
                   for the node with them forbidden (DNNParams::elim)
   -X              print the primal assignment marginals x_ij (XBAR lines)
   -C 1|2          child bounds from the best certificate (DNNResult::child_lb):
                   print how many children it fathoms; 2: also every bound
   -e N            R-step eigensolver: 0 dsyevd (as dnn.py), 1 ssyevd, 2 subspace
                   tracking (default) with a full ssyevd every -r N iterations
                   (20) and a Chebyshev filter of degree -d N (2)
   -A N            Halpern anchoring with restarts at most N iterations apart
                   (DNNParams::halpern; 0 = off)
   -W N            Anderson acceleration with memory N (DNNParams::anderson; 0 = off)
   -L N            subspace tracking: the filter's upper end from N power steps
                   instead of the Weyl bound (DNNParams::tight_lam; 0 = Weyl)
   -f              same as -e 1
   -B X            ADMM penalty factor (beta = X * m/3; default 0.3, dnn.py 1)
   -g X            dual step length gamma (1.618)
   -T N            BLAS threads (1)
   -v              print the certified bound at every check

 The last line is "RESULT key=value ...", for scripts.
============================================================
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include "qap.h"
#include "assign.h"
#include "util.h"
#include "dnn.h"

static void usage()
{
  printf("usage: dnntest [-x i:j,... | -N nodefile] [-u ub | -c cutoff] [-i iters] [-t secs]\n"
         "               [-k every] [-b] [-l 0|1|2] [-P] [-f] [-s window] [-H iters] [-T threads]\n"
         "               [-R rho] [-C 1|2] [-E i:j,...] [-X] [-v] file.dat\n");
  exit(1);
}

// Parse "i:j,i:j,..." ("-" or "" for none) into node. False on bad input.
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

// Bound one node; returns false on bad input or a DNN error. With child
// ("i:j", original indices), the node's final ADMM state then bounds that
// child twice, from scratch and warm-started (results 1 and 2).
static bool print_xbar = false;   // -X: print the primal assignment marginals

static bool run_node(QAP *q, const std::string &fixes, double cutoff, const DNNParams &par0,
                     int index, const std::string &child = "",
                     const std::string &elims = "-")
{
  const int n = q->A->m;
  QAPAssignment node(n);
  if (!parse_fixes(fixes, node, n)) {
    printf("bad assignment list '%s'\n", fixes.c_str());
    return false;
  }
  QAP qr;
  reduceQAP(q, &qr, node);
  // eliminated assignments (original indices), as an m x m mask in local ones
  DNNParams par = par0;
  std::vector<unsigned char> mask;
  if (elims != "-" && !elims.empty()) {
    const int m = n - node.nfix;
    std::vector<int> lr(n, -1), lc(n, -1);
    for (int t = 0, k = 0; t < n; t++) if (node.p[t] == NONE) lr[t] = k++;
    for (int t = 0, k = 0; t < n; t++) if (node.pinv[t] == NONE) lc[t] = k++;
    mask.assign((size_t)m * m, 0);
    for (size_t pos = 0; pos < elims.size();) {
      size_t end = elims.find(',', pos);
      if (end == std::string::npos) end = elims.size();
      int i, j;
      if (sscanf(elims.substr(pos, end - pos).c_str(), "%d:%d", &i, &j) != 2 ||
          i < 0 || i >= n || j < 0 || j >= n || lr[i] < 0 || lc[j] < 0) {
        printf("bad elimination list '%s' (pairs must be free at the node)\n", elims.c_str());
        return false;
      }
      mask[(size_t)lr[i] * m + lc[j]] = 1;
      pos = end + 1;
    }
    par.elim = mask.data();
  }
  if (par.verbose)
    printf("node %d fixed=%d m=%d shift=%.6f cutoff=%.6f\n", index, node.nfix,
           n - node.nfix, qr.shift, cutoff);

  DNNWarm warm;
  DNNResult r = dnn_bound(qr.A, qr.B, qr.C, qr.shift, cutoff, par, NULL,
                          (child.empty() && !print_xbar) ? NULL : &warm);

  printf("RESULT node=%d m=%d lb=%.17g primal=%.17g status=%s iters=%d cpu=%.3f"
         " ms_per_iter=%.2f certs=%d cert_fail=%d refine=%d margin=%.3g fathomed=%d"
         " restarts=%d\n",
         index, n - node.nfix, r.lb, r.primal, dnn_status_name(r.status), r.iterations,
         r.cpu, r.iterations > 0 ? 1e3 * r.cpu / r.iterations : 0.0, r.certificates,
         r.cert_failures, r.refinements, r.margin, r.status == DNN_FATHOMED ? 1 : 0,
         r.halpern_restarts);
  if (!r.child_lb.empty()) {
    // child bounds (-C): how many children the node's certificate fathoms;
    // with -C 2 every bound, as "CHILD i j lb" in original indices
    const int m = n - node.nfix;
    std::vector<int> fr, fc;
    for (int t = 0; t < n; t++) if (node.p[t] == NONE) fr.push_back(t);
    for (int t = 0; t < n; t++) if (node.pinv[t] == NONE) fc.push_back(t);
    int over = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        if (r.child_lb[i * m + j] > cutoff) over++;
        if (par.child_bounds >= 2)
          printf("CHILD %d %d %.17g\n", fr[i], fc[j], r.child_lb[i * m + j]);
      }
    printf("CHILDREN node=%d fathomed=%d of=%d\n", index, over, m * m);
  }
  if (print_xbar && !warm.Y.empty()) {
    // -X: x_ij = Y_0,(ij) of the last primal iterate (a relaxation solution,
    // not a bound), as "XBAR i j value" in original indices
    const int m = n - node.nfix, n1 = m * m + 1;
    std::vector<int> fr, fc;
    for (int t = 0; t < n; t++) if (node.p[t] == NONE) fr.push_back(t);
    for (int t = 0; t < n; t++) if (node.pinv[t] == NONE) fc.push_back(t);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++)
        printf("XBAR %d %d %.6f\n", fr[i], fc[j], warm.Y[1 + i * m + j]);
    (void)n1;
  }
  fflush(stdout);
  if (child.empty() || r.status == DNN_ERROR) return r.status != DNN_ERROR;

  int ci, cj;
  if (sscanf(child.c_str(), "%d:%d", &ci, &cj) != 2 || ci < 0 || ci >= n || cj < 0 ||
      cj >= n || node.p[ci] != NONE || node.pinv[cj] != NONE) {
    printf("bad child '%s'\n", child.c_str());
    return false;
  }
  int il = 0, jl = 0;   // local indices in the node's reduced problem
  for (int t = 0; t < ci; t++) if (node.p[t] == NONE) il++;
  for (int t = 0; t < cj; t++) if (node.pinv[t] == NONE) jl++;
  DNNWarm cw;
  if (!dnn_warm_child(warm, il, jl, cw)) { printf("no warm start\n"); return false; }
  node.fix(ci, cj);
  QAP qc;
  reduceQAP(q, &qc, node);
  for (int pass = 1; pass <= 2; pass++) {
    DNNResult rc = dnn_bound(qc.A, qc.B, qc.C, qc.shift, cutoff, par, pass == 2 ? &cw : NULL);
    printf("RESULT node=%d m=%d lb=%.17g primal=%.17g status=%s iters=%d cpu=%.3f"
           " warm=%d\n", pass, n - node.nfix, rc.lb, rc.primal, dnn_status_name(rc.status),
           rc.iterations, rc.cpu, pass == 2);
  }
  return true;
}

int main(int argc, char **argv)
{
  DNNParams par;
  std::string fixes = "-", child, elims = "-";
  const char *nodefile = NULL;
  double ub = HUGE_VAL, cutoff = HUGE_VAL;
  bool have_cutoff = false;
  int threads = 1;
  char *file = NULL;

  for (int a = 1; a < argc; a++) {
    std::string o = argv[a];
    bool more = a + 1 < argc;
    if (o == "-x" && more) fixes = argv[++a];
    else if (o == "-K" && more) child = argv[++a];
    else if (o == "-N" && more) nodefile = argv[++a];
    else if (o == "-u" && more) ub = atof(argv[++a]);
    else if (o == "-c" && more) { cutoff = atof(argv[++a]); have_cutoff = true; }
    else if (o == "-i" && more) par.max_iter = atoi(argv[++a]);
    else if (o == "-t" && more) par.time_limit = atof(argv[++a]);
    else if (o == "-k" && more) par.check_every = atoi(argv[++a]);
    else if (o == "-b") par.poly = DNN_POLY_BOX;
    else if (o == "-l" && more) par.inner_lap = atoi(argv[++a]);
    else if (o == "-P") {
      par.prototype_cert = 1; par.poly = DNN_POLY_BOX;
      par.check_every = 50; par.stall_window = 6; par.beta = 1; par.hopeless_after = 0;
      par.eig_mode = DNN_EIG_DOUBLE;
    }
    else if (o == "-s" && more) par.stall_window = atoi(argv[++a]);
    else if (o == "-H" && more) par.hopeless_after = atoi(argv[++a]);
    else if (o == "-R" && more) par.hopeless_gap_ratio = atof(argv[++a]);
    else if (o == "-T" && more) threads = atoi(argv[++a]);
    else if (o == "-f") par.eig_mode = DNN_EIG_SINGLE;
    else if (o == "-e" && more) par.eig_mode = atoi(argv[++a]);
    else if (o == "-r" && more) par.refresh = atoi(argv[++a]);
    else if (o == "-d" && more) par.cheb_degree = atoi(argv[++a]);
    else if (o == "-L" && more) par.tight_lam = atoi(argv[++a]);
    else if (o == "-A" && more) par.halpern = atoi(argv[++a]);
    else if (o == "-W" && more) par.anderson = atoi(argv[++a]);
    else if (o == "-B" && more) par.beta = atof(argv[++a]);
    else if (o == "-g" && more) par.gamma = atof(argv[++a]);
    else if (o == "-G" && more) par.gamma1 = atof(argv[++a]);
    else if (o == "-v") par.verbose++;
    else if (o == "-C" && more) par.child_bounds = atoi(argv[++a]);
    else if (o == "-E" && more) elims = argv[++a];
    else if (o == "-X") print_xbar = true;
    else if (o[0] == '-') usage();
    else file = argv[a];
  }
  if (!file) usage();
  dnn_set_blas_threads(threads);

  QAP *q = readQAPLIB(file);
  if (!q) return 1;
  set_fathom_tolerance(q);   // sets DELTA, as the solver does
  symmetrize_QAP(q);
  if (!have_cutoff) cutoff = ub < HUGE_VAL ? ub - DELTA : HUGE_VAL;
  printf("instance %s n=%d lapack=%s\n", file, (int)q->A->m, dnn_lapack_name());

  bool ok = true;
  if (!nodefile) {
    ok = run_node(q, fixes, cutoff, par, 0, child, elims);
  } else {
    // one node per line: "<fixes or -> [cutoff [eliminations]]"; a missing cutoff means the
    // command-line one; '#' starts a comment
    FILE *f = fopen(nodefile, "r");
    if (!f) { printf("cannot open %s\n", nodefile); return 1; }
    char line[65536];
    int index = 0;
    while (fgets(line, sizeof line, f)) {
      char fx[65536];
      double c = cutoff;
      if (line[0] == '#' || sscanf(line, "%s", fx) != 1) continue;
      char cs[64];
      static char el[65536];
      if (sscanf(line, "%*s %63s", cs) == 1) c = atof(cs);
      // optional third field: eliminated assignments "i:j,..." or "-"
      std::string e = sscanf(line, "%*s %*s %65535s", el) == 1 ? el : "-";
      ok = run_node(q, fx, c, par, index++, "", e) && ok;
    }
    fclose(f);
  }
  delete q;
  return ok ? 0 : 2;
}
