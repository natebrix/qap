/*
============================================================
 dnntest: run the DNN bound (dnn.h) on a QAPLIB instance or on one
 node of it, and print the certified bound.

   ./dnntest [options] file.dat

   -x i:j,i:j,...  fix facility i to location j (0-based) before bounding
   -K i:j          after the node, bound its child fixing i to j, cold and
                   warm-started from the node's ADMM state
   -Q list         with -K, the child's starts, comma-separated (default c,z):
                   c cold, z the parent's Z, a Z and the parent's Y restricted,
                   b Z and the parent's Y conditioned on the fixed assignment
                   (-K also takes several children, i:j,i:j,... or i:j;i:j;..., each started
                   from the same parent state)
   -J N            with -K: the parent runs N iterations with no give-up (or
                   until it fathoms); the children keep -H and -i (default:
                   the parent runs as the children)
   -Wr 16|32       with -K: round the parent's saved Y and Z to bf16 (16) or fp32
                   (32) before the children start from them
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
   -Z z[:d]        a give-up within fraction z of the cutoff is deferred d (100)
                   iterations and must then be confirmed (DNNParams::hopeless_zone)
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
   -M eta[:cap]    fast ADMM: Nesterov momentum on (Y, Z), restarted when the combined
                   residual grows by more than eta (1.05), capped at cap (1.0);
                   DNNParams::momentum, momentum_cap (0 = off)
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
#include <cstdint>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include "qap.h"
#include "assign.h"
#include "util.h"
#include "dnn.h"
#include "settings.h"
#include <iostream>

static void usage()
{
  printf("usage: dnntest [-x i:j,... | -N nodefile] [-u ub | -c cutoff] [-i iters] [-t secs]\n"
         "               [-p params.par: its [settings] dnn.* keys, under any flag given]\n"
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

static std::string ladd_spec, dumpy;   // -Ladd file[:const], -dumpY prefix
static std::string child_passes = "c,z";   // -Q: the child's starts
static int parent_iters = 0;   // -J: the parent's own iteration budget, no give-up (0: as the children)
static int round_warm = 0;     // -Wr 16|32: round the parent's saved Y and Z to bf16 or fp32
                               // before the children start from them (a shipped warm start)

// round to bf16 (8 significant bits, round to nearest even) or to fp32
static double round_to(double x, int bits)
{
  float f = (float)x;
  if (bits == 32 || !std::isfinite(f)) return f;
  uint32_t u;
  memcpy(&u, &f, 4);
  u += 0x7fff + ((u >> 16) & 1);
  u &= 0xffff0000u;
  memcpy(&f, &u, 4);
  return f;
}

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
  // an external objective term (externally separated cuts at fixed multipliers)
  std::vector<double> Lext;
  double ext_const = 0;
  if (!ladd_spec.empty()) {
    const int m = n - node.nfix;
    const size_t n1 = (size_t)m * m + 1;
    size_t c = ladd_spec.rfind(':');
    std::string path = ladd_spec;
    if (c != std::string::npos && c > 0 && ladd_spec.find('/', c) == std::string::npos) {
      ext_const = atof(ladd_spec.substr(c + 1).c_str()); path = ladd_spec.substr(0, c); }
    Lext.resize(n1 * n1);
    FILE *f = fopen(path.c_str(), "rb");
    if (!f || fread(Lext.data(), sizeof(double), n1 * n1, f) != n1 * n1) {
      printf("cannot read %zu doubles from %s\n", n1 * n1, path.c_str()); return false; }
    fclose(f);
    double abs_ = 0; size_t nz = 0;
    for (double v : Lext) { abs_ += fabs(v); nz += v != 0; }
    par.L_add = Lext.data(); par.L_add_abs = abs_;
    printf("external term %s: %zu nonzeros, |L_add| = %.6f, const %.6f\n", path.c_str(), nz, abs_, ext_const);
  }

  DNNWarm warm;
  // -J N (with -K): the parent runs exactly N iterations with no give-up (or
  // until it fathoms); the children keep -H / -i (roadmap/convergence2)
  DNNParams ppar = par;
  if (parent_iters > 0 && !child.empty()) { ppar.max_iter = parent_iters; ppar.hopeless_after = 0; }
  DNNResult r = dnn_bound(qr.A, qr.B, qr.C, qr.shift - ext_const, cutoff, ppar, NULL,
                          (child.empty() && !print_xbar && dumpy.empty()) ? NULL : &warm);
  if (!dumpy.empty() && !warm.Y.empty()) {
    std::string path = dumpy + ".Y.bin";
    FILE *f = fopen(path.c_str(), "wb");
    if (f) { fwrite(warm.Y.data(), sizeof(double), warm.Y.size(), f); fclose(f);
      printf("dumped Y (%zu doubles, m=%d) to %s\n", warm.Y.size(), warm.m, path.c_str()); }
  }

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

  if (round_warm == 16 || round_warm == 32) {
    for (double &v : warm.Y) v = round_to(v, round_warm);
    for (double &v : warm.Z) v = round_to(v, round_warm);
    printf("WARM-ROUNDED to %s\n", round_warm == 16 ? "bf16" : "fp32");
  }
  // -K i:j,i:j,... (or i:j;i:j;...): several children of the same parent
  // state, in turn
  int pass = 0;
  const QAPAssignment parent_node = node;
  for (size_t cpos = 0; cpos < child.size();) {
  size_t cend = child.find_first_of(";,", cpos);
  if (cend == std::string::npos) cend = child.size();
  const std::string one = child.substr(cpos, cend - cpos);
  cpos = cend + 1;
  node = parent_node;
  int ci, cj;
  if (sscanf(one.c_str(), "%d:%d", &ci, &cj) != 2 || ci < 0 || ci >= n || cj < 0 ||
      cj >= n || node.p[ci] != NONE || node.pinv[cj] != NONE) {
    printf("bad child '%s'\n", one.c_str());
    return false;
  }
  printf("CHILD-SPEC %s\n", one.c_str());
  int il = 0, jl = 0;   // local indices in the node's reduced problem
  for (int t = 0; t < ci; t++) if (node.p[t] == NONE) il++;
  for (int t = 0; t < cj; t++) if (node.pinv[t] == NONE) jl++;
  node.fix(ci, cj);
  QAP qc;
  reduceQAP(q, &qc, node);
  printf("CHILDNODE %s\n", one.c_str());
  for (size_t t = 0; t < child_passes.size(); t++) {
    const char mode = child_passes[t];
    if (mode == ',') continue;
    const int ymode = mode == 'a' ? 1 : mode == 'b' ? 2 : 0;
    DNNWarm cw;
    if (mode != 'c' && !dnn_warm_child(warm, il, jl, cw, ymode)) { printf("no warm start (%c)\n", mode); return false; }
    DNNParams pc = par;
    pc.warm_primal = ymode > 0;
    pass++;
    printf("CHILD-PASS %d start=%c\n", pass, mode);
    fflush(stdout);
    DNNResult rc = dnn_bound(qc.A, qc.B, qc.C, qc.shift, cutoff, pc, mode == 'c' ? NULL : &cw);
    printf("RESULT node=%d m=%d lb=%.17g primal=%.17g status=%s iters=%d cpu=%.3f"
           " warm=%c\n", pass, n - node.nfix, rc.lb, rc.primal, dnn_status_name(rc.status),
           rc.iterations, rc.cpu, mode);
    fflush(stdout);
  }
  }
  return true;
}

// dnntest's DNNParams defaults differ from the solver's; the registry's
// entries say so (settings_override_default), so the dump and a -p file read
// against them. A flag that has a key is registered as the flag's value, so
// the file does not override it; the file's values then go into par.
static std::string num_text(double v)
{
  char b[64];   // the shortest text that reads back as v
  for (int prec = 6; prec <= 17; prec++) {
    snprintf(b, sizeof b, "%.*g", prec, v);
    if (strtod(b, NULL) == v) break;
  }
  return b;
}

static void par_to_settings(const DNNParams &par, void (*set)(const char *, const std::string &))
{
  set("dnn.giveup", num_text(par.hopeless_after));
  set("dnn.gap_ratio", num_text(par.hopeless_gap_ratio));
  set("dnn.momentum", num_text(par.momentum));
  set("dnn.gu_zone", num_text(par.hopeless_zone));
  set("dnn.refresh", num_text(par.refresh));
  set("dnn.cheb", num_text(par.cheb_degree));
  set("dnn.tight_lam", num_text(par.tight_lam));
  set("dnn.beta", num_text(par.beta));
  set("dnn.gamma", num_text(par.gamma));
  set("dnn.gamma1", num_text(par.gamma1));
}

static void set_default(const char *k, const std::string &v) { settings_override_default(k, v); }
static void set_flag(const char *k, const std::string &v)
{
  if (setting(k).text != v) settings_set(k, v, SRC_FLAG);   // only what a flag changed
}

static void settings_to_par(DNNParams &par)
{
  if (setting("dnn.giveup").src == SRC_PAR) par.hopeless_after = (int)setting_i("dnn.giveup");
  if (setting("dnn.gap_ratio").src == SRC_PAR) par.hopeless_gap_ratio = setting_d("dnn.gap_ratio");
  if (setting("dnn.momentum").src == SRC_PAR) par.momentum = setting_d("dnn.momentum");
  if (setting("dnn.gu_zone").src == SRC_PAR) par.hopeless_zone = setting_d("dnn.gu_zone");
  if (setting("dnn.refresh").src == SRC_PAR) par.refresh = (int)setting_i("dnn.refresh");
  if (setting("dnn.cheb").src == SRC_PAR) par.cheb_degree = (int)setting_i("dnn.cheb");
  if (setting("dnn.tight_lam").src == SRC_PAR) par.tight_lam = (int)setting_i("dnn.tight_lam");
  if (setting("dnn.beta").src == SRC_PAR) par.beta = setting_d("dnn.beta");
  if (setting("dnn.gamma").src == SRC_PAR) par.gamma = setting_d("dnn.gamma");
  if (setting("dnn.gamma1").src == SRC_PAR) par.gamma1 = setting_d("dnn.gamma1");
}

int main(int argc, char **argv)
{
  DNNParams par;
  par_to_settings(par, set_default);   // the registry's defaults are dnntest's
  const char *par_file = NULL;
  bool print_only = false;   // --print-settings: the resolved settings, then exit
  std::string fixes = "-", child, elims = "-";
  const char *nodefile = NULL;
  double ub = HUGE_VAL, cutoff = HUGE_VAL;
  bool have_cutoff = false;
  int threads = 1;
  char *file = NULL;

  for (int a = 1; a < argc; a++) {
    std::string o = argv[a];
    bool more = a + 1 < argc;
    if (o == "--print-settings") { print_only = true; continue; }
    if (o == "-p" && more) { par_file = argv[++a]; continue; }
    if (o == "-x" && more) fixes = argv[++a];
    else if (o == "-K" && more) child = argv[++a];
    else if (o == "-Q" && more) child_passes = argv[++a];
    else if (o == "-J" && more) parent_iters = atoi(argv[++a]);
    else if (o == "-Wr" && more) round_warm = atoi(argv[++a]);
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
    else if (o == "-M" && more) {   // -M eta[:cap]
      const char *z = argv[++a], *c = strchr(z, ':');
      par.momentum = atof(z);
      if (c) par.momentum_cap = atof(c + 1);
    }
    else if (o == "-Z" && more) {   // -Z zone[:delay]: confirm a give-up near the cutoff
      const char *z = argv[++a], *c = strchr(z, ':');
      par.hopeless_zone = atof(z);
      if (c) par.hopeless_delay = atoi(c + 1);
    }
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
    else if (o == "-Ladd" && more) ladd_spec = argv[++a];   // external objective term file[:const]
    else if (o == "-dumpY" && more) dumpy = argv[++a];      // write the final iterate to prefix.Y.bin
    else if (o[0] == '-') usage();
    else file = argv[a];
  }
  par_to_settings(par, set_flag);   // the flags, so the file does not override them
  if (par_file) {
    settings_read_par(par_file);
    settings_to_par(par);
  }
  if (print_only) { settings_dump(std::cout); return 0; }
  if (!file) usage();
  settings_dump(std::cout);   // the settings this run resolved (settings.h)
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
