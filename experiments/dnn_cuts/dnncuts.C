/*
 dnncuts: the certified DNN bound strengthened by triangle inequalities, in
 Lagrangian form (see README.md, ../wil50_feasibility/stronger_bounds/README.md
 and latest/dnncut.h, which holds the math and the code; this is the driver).

   make -C latest dnncuts
   experiments/dnn_cuts/dnncuts [-x i:j,...] [-u ub | -c cutoff] [-r rounds]
       [-I first-round iters] [-i iters] [-k cuts] [-t theta] [-delta D]
       [-mu0 M] [-grid R:C [-side 0|1] | -prefix | -atomic] [-five K] [-d degree]
       [-Ladd file[:const]] [-warmY] [-dumpY prefix] [-addto R] [-seven K[:S]]
       [-e eigmode] [-g gamma] [-G gamma1] [-T threads] [-v] file.dat

 -grid R:C takes the row and column prefixes of an R x C grid over the
 original indices of one side (default side 0, the first matrix: wil50's
 locations; -side 1 the second matrix: nug30's 5 x 6 grid), mapped to the
 node's free elements. -prefix takes the prefixes of the free elements
 (validation on any instance); -atomic the single-variable triangles of
 README.md; -five K adds K of the regional 5-point count inequalities per
 round on top of the triangles (dnncut.h). -Ladd file[:const] adds an external
 (m^2+1)^2 row-major double matrix (e.g. stronger_bounds/triples/make_ladd.py)
 to every round's objective, const subtracted from the shift; -dumpY prefix
 writes the final primal iterate Y to prefix.Y.bin (raw doubles, row-major).
 Each round: solve (warm-started), separate the k most violated
 cuts in the last primal iterate, Polyak step toward min(cutoff, best + delta);
 stop at the first certified bound above the cutoff. Every value printed is
 certified.
*/

#include "dnncut.h"
#include "qap.h"
#include "assign.h"
#include "util.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

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

int main(int argc, char **argv)
{
  std::string fixes = "-";
  double ub = HUGE_VAL, cutoff = HUGE_VAL;
  bool have_cutoff = false;
  DNNCutParams cp;
  cp.rounds = 15; cp.iters = 400; cp.kcuts = 2000;
  int threads = 1, degree = 4, eig_mode = DNN_EIG_SUBSPACE, verbose = 0;
  double gamma = 1.0, gamma1 = 1.0;
  const char *file = NULL;
  std::string ladd_spec, dumpy;
  for (int a = 1; a < argc; a++) {
    std::string o = argv[a];
    bool more = a + 1 < argc;
    if (o == "-x" && more) fixes = argv[++a];
    else if (o == "-u" && more) ub = atof(argv[++a]);
    else if (o == "-c" && more) { cutoff = atof(argv[++a]); have_cutoff = true; }
    else if (o == "-r" && more) cp.rounds = atoi(argv[++a]);
    else if (o == "-i" && more) cp.iters = atoi(argv[++a]);
    else if (o == "-I" && more) cp.first_iters = atoi(argv[++a]);
    else if (o == "-k" && more) cp.kcuts = atoi(argv[++a]);
    else if (o == "-five" && more) cp.kcuts5 = atoi(argv[++a]);               // 5-point count inequalities per round
    else if (o == "-addto" && more) cp.add_until = atoi(argv[++a]);           // add cuts only in rounds <= R
    else if (o == "-seven" && more) {                                          // -seven K[:S]: K local-7 cuts a round, S random subsets per region
      const char *z = argv[++a], *c = strchr(z, ':');
      cp.kcuts7 = atoi(z); if (c) cp.sub7 = atoi(c + 1);
    }
    else if (o == "-t" && more) cp.theta = atof(argv[++a]);
    else if (o == "-delta" && more) cp.delta = atof(argv[++a]);
    else if (o == "-mu0" && more) cp.mu0 = atof(argv[++a]);
    else if (o == "-grid" && more) {
      if (sscanf(argv[++a], "%d:%d", &cp.grid_rows, &cp.grid_cols) != 2 || cp.grid_rows < 1 || cp.grid_cols < 1) {
        fprintf(stderr, "bad -grid R:C\n"); return 1; }
    }
    else if (o == "-side" && more) cp.side = atoi(argv[++a]);
    else if (o == "-stall" && more) cp.stall_rounds = atoi(argv[++a]);   // stop rule (dnncut.h), 0 = off
    else if (o == "-nocap") cp.cap_at_cutoff = 0;                          // Polyak target = best + delta
    else if (o == "-freeze" && more) {                                     // -freeze R:N: quiet round of N iterations after round R (0 = the last)
      const char *z = argv[++a], *c = strchr(z, ':');
      cp.freeze_after = atoi(z);
      cp.freeze_iters = c ? atoi(c + 1) : 600;
    }
    else if (o == "-freeze-last") cp.freeze_best = 0;                     // quiet round from the last round's multipliers, not the best round's
    else if (o == "-warmY") cp.warm_primal = 1;
    else if (o == "-Ladd" && more) ladd_spec = argv[++a];
    else if (o == "-dumpY" && more) dumpy = argv[++a];
    else if (o == "-prefix") cp.prefix = true;
    else if (o == "-atomic") cp.atomic = true;
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
  if (cp.grid_rows == 0 && !cp.prefix) cp.atomic = true;   // as before: no regions means atomic
  if (cp.first_iters <= 0) cp.first_iters = cp.iters;       // the driver always runs the uncut round
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
  if (cp.grid_rows > 0 && cp.grid_rows * cp.grid_cols != n) { fprintf(stderr, "-grid R:C needs R*C = n\n"); return 1; }
  QAP qr;
  reduceQAP(q, &qr, node);
  const int m = n - node.nfix;
  // the free elements of the region side, original indices in local order
  std::vector<int> free_orig;
  for (int a = 0; a < n; a++)
    if ((cp.side == 0 ? node.p[a] : node.pinv[a]) == NONE) free_orig.push_back(a);
  std::vector<double> Lext;
  if (!ladd_spec.empty()) {
    size_t c = ladd_spec.rfind(':');
    std::string path = ladd_spec;
    if (c != std::string::npos && c > 0 && ladd_spec.find('/', c) == std::string::npos) {
      cp.ext_const = atof(ladd_spec.substr(c + 1).c_str()); path = ladd_spec.substr(0, c); }
    const size_t n1 = (size_t)m * m + 1;
    Lext.resize(n1 * n1);
    FILE *f = fopen(path.c_str(), "rb");
    if (!f || fread(Lext.data(), sizeof(double), n1 * n1, f) != n1 * n1) {
      fprintf(stderr, "cannot read %zu doubles from %s\n", n1 * n1, path.c_str()); return 1; }
    fclose(f);
    double abs_ = 0; size_t nz = 0;
    for (double v : Lext) { abs_ += fabs(v); nz += v != 0; }
    cp.L_ext = Lext.data(); cp.L_ext_abs = abs_;
    printf("external term %s: %zu nonzeros, |L_ext| = %.6f, const %.6f\n", path.c_str(), nz, abs_, cp.ext_const);
  }
  printf("instance %s n=%d node %s m=%d cutoff %.4f  cuts %s", file, n, fixes.c_str(), m, cutoff,
         cp.atomic ? "atomic" : (cp.grid_rows > 0 ? "regional (grid)" : "regional (prefixes)"));
  if (cp.grid_rows > 0) printf(" side %d", cp.side);
  if (cp.warm_primal) printf(" warmY 1");
  DNNParams base;
  base.eig_mode = eig_mode;
  base.cheb_degree = degree;
  base.gamma = gamma; base.gamma1 = gamma1;
  cp.verbose = std::max(1, verbose);
  DNNWarm warm;
  if (!cp.atomic) {
    // count the regions for the banner (the library builds them again)
    DNNCutParams probe = cp; probe.rounds = -1; probe.first_iters = 0; probe.verbose = 0;
    DNNCutResult pr = dnn_cut_rounds(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, probe, &warm);
    printf(", %d regions", pr.nregions);
  }
  printf("\n");
  DNNCutResult res = dnn_cut_rounds(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, cp, &warm);
  printf("FINAL best_lb %.6f first_lb %.6f gain %.6f\n", res.best_lb, res.first_lb,
         res.best_lb - res.first_lb);
  if (!dumpy.empty() && !warm.Y.empty()) {
    std::string path = dumpy + ".Y.bin";
    FILE *f = fopen(path.c_str(), "wb");
    if (f) { fwrite(warm.Y.data(), sizeof(double), warm.Y.size(), f); fclose(f);
      printf("dumped Y (%zu doubles, m=%d) to %s\n", warm.Y.size(), warm.m, path.c_str()); }
  }
  delete q;
  return 0;
}
