/*
 dnncuts: the certified DNN bound strengthened by triangle inequalities, in
 Lagrangian form (see README.md, ../wil50_feasibility/stronger_bounds/README.md
 and latest/dnncut.h, which holds the math and the code; this is the driver).

   make -C latest dnncuts
   experiments/dnn_cuts/dnncuts [-x i:j,...] [-u ub | -c cutoff] [-r rounds]
       [-I first-round iters] [-i iters] [-k cuts] [-t theta] [-delta D]
       [-mu0 M] [-grid R:C [-side 0|1] | -prefix | -atomic] [-five K] [-d degree]
       [-Ladd file[:const]] [-warmY] [-dumpY prefix] [-addto R] [-seven K[:S]]
       [-stall R] [-pgate R:M] [-box K[:QMAX[:STARTS]]] [-loadpool F] [-savepool F]
       [-poolcheck K] [-e eigmode] [-g gamma] [-G gamma1] [-T threads] [-v] file.dat
       [-savepoolall R:F] [-cpadmm F [-cpiters N] [-cpsweeps S] [-cppasses K]
       [-cpdiag E[:S]] [-cpcheck]]

 -grid R:C takes the row and column prefixes of an R x C grid over the
 original indices of one side (default side 0, the first matrix: wil50's
 locations; -side 1 the second matrix: nug30's 5 x 6 grid), mapped to the
 node's free elements. -prefix takes the prefixes of the free elements
 (validation on any instance); -atomic the single-variable triangles of
 README.md; -five K adds K of the regional 5-point count inequalities per
 round on top of the triangles (dnncut.h); -box K adds K box (count)
 inequalities per round (largest |Q| QMAX, default 24; STARTS local-search
 starts per region and b, default 40). The search's per-node settings
 QAP_DNN_CHEB, QAP_DNN_TIGHT_LAM and QAP_DNN_REFRESH are read from the
 environment as bnb does (QAP_DNN_CHEB overrides -d). -loadpool F: a
 parent's cut pool (dnn_pool_encode's token, or Astra's
 QAP_REGIONAL_TRIANGLE_POOL v1 file, experiments/pool_inheritance), remapped
 to the node: the first solve (-I iterations) uses it frozen, as the search's
 QAP_DNN_POOL does, then -r rounds continue from it (-r 0: that solve only).
 -savepool F writes the pool the node would hand to its children. -Ladd file[:const] adds an external
 (m^2+1)^2 row-major double matrix (e.g. stronger_bounds/triples/make_ladd.py)
 to every round's objective, const subtracted from the shift; -dumpY prefix
 writes the final primal iterate Y to prefix.Y.bin (raw doubles, row-major).
 Each round: solve (warm-started), separate the k most violated
 cuts in the last primal iterate, Polyak step toward min(cutoff, best + delta);
 stop at the first certified bound above the cutoff. Every value printed is
 certified.

 CP-ADMM (experiments/wil50_feasibility/roadmap/admm/cpadmm): -savepoolall R:F
 writes the whole pool after round R's new cuts (zero multipliers included) to
 F. -cpadmm F runs the uncut round (-I iterations, as the rounds' round 0: the
 give-up point) and then one solve of -cpiters N iterations (2000) from its
 state with pool F frozen inside the Y-step (dnn_cp_bound, latest/dnncut.h):
 -cpsweeps S Dykstra sweeps per iteration (1), -cppasses K half-space passes per
 sweep (1), -cpdiag E[:S] every E iterations the distance to an S-sweep
 projection (30), -cpcheck the cut extraction self-test first. The per-check
 trace (-vc is implied) is certified at every check.
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
#include <sstream>
#include <cctype>

// a pool file: dnn_pool_encode's token, or Astra's v1 file
// (QAP_REGIONAL_TRIANGLE_POOL 1: n side, the assignment, region-side and
// other-side original index lists, provenance, regions, cuts "type s a b c mu")
static bool read_pool_file(const std::string &path, DNNPool &p)
{
  FILE *f = fopen(path.c_str(), "r");
  if (!f) return false;
  std::string text;
  char buf[1 << 16];
  size_t k;
  while ((k = fread(buf, 1, sizeof buf, f)) > 0) text.append(buf, k);
  fclose(f);
  if (text.compare(0, 26, "QAP_REGIONAL_TRIANGLE_POOL") != 0) {
    while (!text.empty() && isspace((unsigned char)text.back())) text.pop_back();
    return dnn_pool_decode(text, p);
  }
  std::istringstream in(text);
  std::string magic; int version, n, side;
  if (!(in >> magic >> version >> n >> side) || version != 1) return false;
  auto list = [&](std::vector<int> &v) { int c; if (!(in >> c) || c < 0) return false; v.resize(c);
    for (int &x : v) if (!(in >> x)) return false; return true; };
  std::vector<int> asg, orig, other;
  if (!list(asg) || !list(orig) || !list(other)) return false;
  int round, quiet; double value;
  if (!(in >> round >> quiet >> value)) return false;
  int nr; if (!(in >> nr) || nr < 0) return false;
  p = DNNPool(); p.side = side;
  for (int s = 0; s < nr; s++) { std::vector<int> r; if (!list(r)) return false; p.regions.push_back(r); }
  int nc; if (!(in >> nc) || nc < 0) return false;
  for (int c = 0; c < nc; c++) {
    DNNPoolCut pc; int a, b, d;
    if (!(in >> pc.type >> pc.region >> a >> b >> d >> pc.mu)) return false;
    pc.b = 0; pc.mem = {a, b, d};
    p.cuts.push_back(pc);
  }
  std::string end;
  return (bool)(in >> end) && end == "END";
}

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
  int hopeless_after = -1;   // -H N: the give-up of the -loadpool first solve (default: DNNParams')
  int check_trace = 0;       // -vc: per-check trace of every solve (DNNParams::verbose)
  double gamma = 1.0, gamma1 = 1.0;
  const char *file = NULL;
  std::string ladd_spec, dumpy, loadpool, savepool;
  std::string savepoolall, cpadmm;   // CP-ADMM study
  int savepoolall_round = -1, cp_iters = 2000, cp_plain = 0;
  DNNCPParams cpp;
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
    else if (o == "-damp1" && more) {   // F[:G]: first and later steps of an inherited phase
      const char *e = argv[++a];
      cp.damp1 = atof(e);
      if (const char *c = strchr(e, ':')) cp.damp_rest = atof(c + 1);
    }
    else if (o == "-mu0" && more) cp.mu0 = atof(argv[++a]);
    else if (o == "-grid" && more) {
      if (sscanf(argv[++a], "%d:%d", &cp.grid_rows, &cp.grid_cols) != 2 || cp.grid_rows < 1 || cp.grid_cols < 1) {
        fprintf(stderr, "bad -grid R:C\n"); return 1; }
    }
    else if (o == "-side" && more) cp.side = atoi(argv[++a]);
    else if (o == "-stall" && more) cp.stall_rounds = atoi(argv[++a]);   // stop rule (dnncut.h), 0 = off
    else if (o == "-pgate" && more) {                                      // -pgate R:M: primal gate after round R, margin M (dnncut.h)
      const char *z = argv[++a], *c = strchr(z, ':');
      cp.pgate_rounds = atoi(z); if (c) cp.pgate_margin = atof(c + 1);
    }
    else if (o == "-nocap") cp.cap_at_cutoff = 0;                          // Polyak target = best + delta
    else if (o == "-freeze" && more) {                                     // -freeze R:N: quiet round of N iterations after round R (0 = the last)
      const char *z = argv[++a], *c = strchr(z, ':');
      cp.freeze_after = atoi(z);
      cp.freeze_iters = c ? atoi(c + 1) : 600;
    }
    else if (o == "-freeze-last") cp.freeze_best = 0;                     // quiet round from the last round's multipliers, not the best round's
    else if (o == "-warmY") cp.warm_primal = 1;
    else if (o == "-box" && more) {                                        // -box K[:QMAX[:STARTS]]: K box (count) inequalities a round (dnncut.h)
      const char *z = argv[++a], *c = strchr(z, ':');
      cp.kbox = atoi(z);
      if (c) { cp.box_qmax = atoi(c + 1); if (const char *c2 = strchr(c + 1, ':')) cp.box_starts = atoi(c2 + 1); }
    }
    else if (o == "-Ladd" && more) ladd_spec = argv[++a];
    else if (o == "-loadpool" && more) loadpool = argv[++a];              // a parent's cut pool (dnncut.h DNNPool; also Astra's v1 file)
    else if (o == "-savepool" && more) savepool = argv[++a];              // write the pool handed to the children
    else if (o == "-poolcheck" && more) cp.pool_check = atoi(argv[++a]);
    else if (o == "-savepoolall" && more) {                                // R:F the whole pool after round R
      const char *z = argv[++a], *c = strchr(z, ':');
      if (!c) { fprintf(stderr, "bad -savepoolall R:F\n"); return 1; }
      savepoolall_round = atoi(z); savepoolall = c + 1;
    }
    else if (o == "-cpadmm" && more) cpadmm = argv[++a];
    else if (o == "-cpiters" && more) cp_iters = atoi(argv[++a]);
    else if (o == "-cpsweeps" && more) cpp.sweeps = atoi(argv[++a]);
    else if (o == "-cppasses" && more) cpp.passes = atoi(argv[++a]);
    else if (o == "-cpdiag" && more) {
      const char *z = argv[++a], *c = strchr(z, ':');
      cpp.diag_every = atoi(z); if (c) cpp.diag_sweeps = atoi(c + 1);
    }
    else if (o == "-cpcheck") cpp.selftest = 1;
    else if (o == "-cpplain") cp_plain = 1;   // the same solve from the give-up point without the cuts (no hook)
    else if (o == "-dumpY" && more) dumpy = argv[++a];
    else if (o == "-prefix") cp.prefix = true;
    else if (o == "-atomic") cp.atomic = true;
    else if (o == "-d" && more) degree = atoi(argv[++a]);
    else if (o == "-e" && more) eig_mode = atoi(argv[++a]);
    else if (o == "-g" && more) gamma = atof(argv[++a]);
    else if (o == "-G" && more) gamma1 = atof(argv[++a]);
    else if (o == "-T" && more) threads = atoi(argv[++a]);
    else if (o == "-H" && more) hopeless_after = atoi(argv[++a]);    // give-up of the -loadpool solve (the search: 300)
    else if (o == "-vc") check_trace = 1;                             // per-check trace of every solve
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
  if (cp.damp1 != 1.0 || cp.damp_rest != 1.0) printf(" damp1 %g:%g", cp.damp1, cp.damp_rest);
  if (cp.kbox > 0) printf(" box %d qmax %d starts %d", cp.kbox, cp.box_qmax, cp.box_starts);
  DNNParams base;
  base.eig_mode = eig_mode;
  base.cheb_degree = degree;
  base.gamma = gamma; base.gamma1 = gamma1;
  if (const char *e = getenv("QAP_DNN_CHEB")) base.cheb_degree = atoi(e);
  if (const char *e = getenv("QAP_DNN_TIGHT_LAM")) base.tight_lam = atoi(e);
  if (const char *e = getenv("QAP_DNN_REFRESH")) base.refresh = atoi(e);
  if (hopeless_after >= 0) base.hopeless_after = hopeless_after;
  base.verbose = check_trace;   // -vc: the certified bound at every check of every solve
  cp.verbose = std::max(1, verbose);
  DNNWarm warm;
  if (!cp.atomic) {
    // count the regions for the banner (the library builds them again)
    DNNCutParams probe = cp; probe.rounds = -1; probe.first_iters = 0; probe.verbose = 0;
    DNNCutResult pr = dnn_cut_rounds(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, probe, &warm);
    printf(", %d regions", pr.nregions);
  }
  printf("\n");
  // the other side's free elements and where its fixed ones sit (pools)
  std::vector<int> other_orig, fixed_at(n, -1);
  for (int a = 0; a < n; a++) {
    if ((cp.side == 0 ? node.pinv[a] : node.p[a]) == NONE) other_orig.push_back(a);
    else fixed_at[a] = cp.side == 0 ? node.pinv[a] : node.p[a];
  }
  DNNPool pin, pout, pall;
  if (!savepool.empty()) { cp.pool_out = &pout; cp.other_orig = &other_orig; }
  if (!savepoolall.empty()) { cp.pool_all_round = savepoolall_round; cp.pool_all_out = &pall; cp.other_orig = &other_orig; }
  if (!cpadmm.empty() || cp_plain) {
    // CP-ADMM: the give-up point (the uncut round 0 of the rounds), then one solve
    // with the frozen pool inside the Y-step (-cpplain: the same solve, plain Y-step)
    if (!cp_plain && (!read_pool_file(cpadmm, pin) || pin.side != cp.side)) { fprintf(stderr, "cannot read pool %s (or wrong side)\n", cpadmm.c_str()); return 1; }
    size_t pos = 0; for (const auto &c : pin.cuts) pos += c.mu > 0;
    printf("cp pool %s: %zu regions, %zu cuts (%zu with mu > 0)\n", cpadmm.c_str(), pin.regions.size(), pin.cuts.size(), pos);
    DNNCutParams c0 = cp;
    c0.rounds = 0; c0.freeze_iters = 0; c0.freeze_after = 0;
    c0.pool_in = NULL; c0.pool_out = NULL; c0.pool_all_out = NULL;
    DNNCutResult r0 = dnn_cut_rounds(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, c0, &warm);
    DNNParams p = base;
    p.max_iter = cp_iters;
    p.hopeless_after = 0; p.hopeless_gap_ratio = 0;   // as the cut rounds
    p.warm_primal = 0;
    p.cut_round = 1;
    p.verbose = std::max(1, check_trace);
    if (verbose >= 2) p.verbose = 2;
    DNNCutParams c1 = cp;
    c1.pool_in = &pin; c1.other_orig = &other_orig; c1.fixed_at = &fixed_at; c1.pool_out = NULL; c1.pool_all_out = NULL;
    DNNWarm save;
    int kept = 0;
    DNNResult r = cp_plain ? dnn_bound(qr.A, qr.B, qr.C, qr.shift, cutoff, p, &warm, &save)
                           : dnn_cp_bound(qr.A, qr.B, qr.C, qr.shift, cutoff, p, free_orig, n, c1, cpp, &warm, &save, &kept);
    printf("cp solve  lb %.10f (gain %.4f)  primal %.4f  status %s  iters %d  cpu %.1f  certs %d failures %d  kept %d\n",
           r.lb, r.lb - r0.first_lb, r.primal, dnn_status_name(r.status), r.iterations, r.cpu, r.certificates,
           r.cert_failures, kept);
    printf("FINAL best_lb %.6f first_lb %.6f gain %.6f\n", std::max(r.lb, r0.first_lb), r0.first_lb,
           std::max(r.lb, r0.first_lb) - r0.first_lb);
    delete q;
    return 0;
  }
  if (!loadpool.empty()) {
    if (!read_pool_file(loadpool, pin) || pin.side != cp.side) { fprintf(stderr, "cannot read pool %s (or wrong side)\n", loadpool.c_str()); return 1; }
    size_t pos = 0; for (const auto &c : pin.cuts) pos += c.mu > 0;
    printf("pool %s: %zu regions, %zu cuts (%zu with mu > 0)\n", loadpool.c_str(), pin.regions.size(), pin.cuts.size(), pos);
    cp.pool_in = &pin; cp.other_orig = &other_orig; cp.fixed_at = &fixed_at;
    DNNParams p0 = base;
    p0.max_iter = cp.first_iters;
    int kept = 0, dropped = 0;
    DNNCutParams pc = cp;
    if (cp.rounds > 0) pc.pool_out = NULL;
    DNNResult r0 = dnn_pool_bound(qr.A, qr.B, qr.C, qr.shift, cutoff, p0, free_orig, n, pc, NULL, &warm, &kept, &dropped);
    printf("pool solve  lb %.4f  primal %.4f  status %s  iters %d  cpu %.1f  kept %d dropped %d\n", r0.lb, r0.primal,
           dnn_status_name(r0.status), r0.iterations, r0.cpu, kept, dropped);
    if (cp.rounds <= 0 || r0.status == DNN_FATHOMED) {
      printf("FINAL best_lb %.6f first_lb %.6f gain 0.000000\n", r0.lb, r0.lb);
      if (!savepool.empty()) { FILE *f = fopen(savepool.c_str(), "w"); if (f) { fprintf(f, "%s\n", dnn_pool_encode(pout).c_str()); fclose(f); } }
      delete q;
      return 0;
    }
    cp.first_iters = 0;
    DNNCutResult res = dnn_cut_rounds(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, cp, &warm, r0.lb, r0.primal);
    printf("FINAL best_lb %.6f first_lb %.6f gain %.6f\n", std::max(res.best_lb, r0.lb), r0.lb, std::max(res.best_lb, r0.lb) - r0.lb);
    if (!savepool.empty()) { FILE *f = fopen(savepool.c_str(), "w"); if (f) { fprintf(f, "%s\n", dnn_pool_encode(pout).c_str()); fclose(f); } }
    delete q;
    return 0;
  }
  DNNCutResult res = dnn_cut_rounds(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, cp, &warm);
  if (!savepool.empty()) { FILE *f = fopen(savepool.c_str(), "w"); if (f) { fprintf(f, "%s\n", dnn_pool_encode(pout).c_str()); fclose(f); } }
  if (!savepoolall.empty()) {
    FILE *f = fopen(savepoolall.c_str(), "w");
    if (f) { fprintf(f, "%s\n", dnn_pool_encode(pall).c_str()); fclose(f); }
    printf("whole pool after round %d: %zu cuts to %s\n", savepoolall_round, pall.cuts.size(), savepoolall.c_str());
  }
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
