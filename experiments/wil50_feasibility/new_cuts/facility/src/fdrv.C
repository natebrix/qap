/*
 fdrv: driver for the facility-side experiments (new_cuts/facility). Links the
 shared solver objects of latest/ read-only and a renamed copy of dnncut.C
 (dnncut_f.C, made by make_dnncut_f.py) that can export the triangle phase's
 Lagrangian term and ADMM state, so a phase can be split across processes
 (each under the container's restart interval) and new families can be tested
 on top of the frozen triangle pool. Every printed lb is certified: the
 objective is L + sum of valid-inequality terms with nonnegative multipliers,
 the shift lowered by their right-hand sides, and the certificate of dnn.C.

 Modes (common: -x fixes -u ub | -c cutoff, -d degree, -L tight, -r refresh, file.dat):

   -phase P   the triangle protocol of dnncuts (-grid R:C -side S -k K -I first
              -i iters -rounds R -delta D -nocap -t theta -freeze N -warmY),
              stopping right before the quiet round; writes P.Y0.bin (the
              round-0 iterate), P.tri.bin (the pool's term at the quiet round's
              multipliers), P.sY.bin, P.sZ.bin (the state the quiet round
              starts from) and P.meta (m scale hsum abs active best first).
   -cont P    loads P.*; runs -iters N iterations of dnn_bound from that state
              with L + P.tri [+ -extra F:const ...] (-notri: without P.tri;
              -wy: also start from the saved Y); -dumpY Q writes Q.Y.bin,
              -save Q writes the final state as Q.sY/sZ/meta (same tri term).
   -pool P    as -cont, plus external cuts (-cuts F.cb, sparse format below)
              with their own multipliers: -rounds R of -i I iterations, Polyak
              toward best + delta, the step halved after a drop; an optional
              separator command (-sep "cmd {Y} {C}") is run after every round
              on the dumped iterate and its new cuts are added; then a quiet
              round of -freeze N iterations from the best round's multipliers.

 Sparse cut file (.cb): int32 ncuts; per cut int64 key, double h, int32 nnz,
 nnz x (int32 r, int32 c, double g): the inequality sum g Y_rc <= h on the
 lifted matrix (r, c lifted indices, each unordered pair once, Y symmetric).
*/
#include "dnncut.h"
#include "fdrv.h"
#include "qap.h"
#include "assign.h"
#include "util.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

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

static bool write_vec(const std::string &path, const std::vector<double> &v)
{
  FILE *f = fopen(path.c_str(), "wb");
  if (!f) return false;
  bool ok = fwrite(v.data(), sizeof(double), v.size(), f) == v.size();
  fclose(f);
  return ok;
}
static bool read_vec(const std::string &path, std::vector<double> &v, size_t count)
{
  v.resize(count);
  FILE *f = fopen(path.c_str(), "rb");
  if (!f) return false;
  bool ok = fread(v.data(), sizeof(double), count, f) == count;
  fclose(f);
  return ok;
}

struct Meta { int m = 0; double scale = 0, hsum = 0, abs_ = 0, best = 0, first = 0; int active = 0; };
static bool write_meta(const std::string &p, const Meta &mt)
{
  FILE *f = fopen((p + ".meta").c_str(), "w");
  if (!f) return false;
  fprintf(f, "%d %.17g %.17g %.17g %d %.17g %.17g\n", mt.m, mt.scale, mt.hsum, mt.abs_, mt.active, mt.best, mt.first);
  fclose(f);
  return true;
}
static bool read_meta(const std::string &p, Meta &mt)
{
  FILE *f = fopen((p + ".meta").c_str(), "r");
  if (!f) return false;
  bool ok = fscanf(f, "%d %lf %lf %lf %d %lf %lf", &mt.m, &mt.scale, &mt.hsum, &mt.abs_, &mt.active, &mt.best, &mt.first) == 7;
  fclose(f);
  return ok;
}

struct SCut { int64_t key; double h; std::vector<int32_t> r, c; std::vector<double> g; double mu; };
static bool read_cuts(const std::string &path, std::vector<SCut> &out, int n1)
{
  FILE *f = fopen(path.c_str(), "rb");
  if (!f) return false;
  int32_t nc;
  if (fread(&nc, 4, 1, f) != 1) { fclose(f); return false; }
  for (int k = 0; k < nc; k++) {
    SCut s; int32_t nnz;
    if (fread(&s.key, 8, 1, f) != 1 || fread(&s.h, 8, 1, f) != 1 || fread(&nnz, 4, 1, f) != 1) { fclose(f); return false; }
    s.r.resize(nnz); s.c.resize(nnz); s.g.resize(nnz);
    for (int t = 0; t < nnz; t++) {
      if (fread(&s.r[t], 4, 1, f) != 1 || fread(&s.c[t], 4, 1, f) != 1 || fread(&s.g[t], 8, 1, f) != 1) { fclose(f); return false; }
      if (s.r[t] < 0 || s.r[t] >= n1 || s.c[t] < 0 || s.c[t] >= n1) { fclose(f); return false; }
    }
    s.mu = 0;
    out.push_back(s);
  }
  fclose(f);
  return true;
}
static double cut_value(const SCut &s, const std::vector<double> &Y, int n1)
{
  long double v = 0;
  for (size_t t = 0; t < s.g.size(); t++) v += s.g[t] * Y[(size_t)s.r[t] * n1 + s.c[t]];
  return (double)v - s.h;
}
static double add_scut(const SCut &s, std::vector<double> &L, int n1, double mu)
{
  double a = 0;
  for (size_t t = 0; t < s.g.size(); t++) {
    const size_t r = s.r[t], c = s.c[t];
    const double v = mu * s.g[t];
    if (r == c) L[r * n1 + r] += v;
    else { L[r * n1 + c] += v / 2; L[c * n1 + r] += v / 2; }
    a += fabs(v);
  }
  return a;
}

int main(int argc, char **argv)
{
  std::string fixes = "-", phase, cont, pool, dumpy, savep, sepcmd, regfile;
  std::vector<std::string> extras, cutfiles;
  double ub = HUGE_VAL, cutoff = HUGE_VAL;
  bool have_cutoff = false, notri = false, wy = false;
  int childb = 0;
  int iters_cont = 1500, threads = 1, degree = 4, tight = 0, refresh = 20, eig_mode = DNN_EIG_SUBSPACE;
  double gamma = 1.0, gamma1 = 1.0;
  DNNCutParams cp;
  cp.rounds = 12; cp.iters = 300; cp.kcuts = 1000; cp.first_iters = 400; cp.delta = 50;
  int pool_rounds = 8, pool_iters = 300, pool_freeze = 1500;
  double pool_delta = 50, pool_theta = 1.0, pool_mu0 = 0;
  const char *file = NULL;
  for (int a = 1; a < argc; a++) {
    std::string o = argv[a];
    bool more = a + 1 < argc;
    if (o == "-x" && more) fixes = argv[++a];
    else if (o == "-u" && more) ub = atof(argv[++a]);
    else if (o == "-c" && more) { cutoff = atof(argv[++a]); have_cutoff = true; }
    else if (o == "-d" && more) degree = atoi(argv[++a]);
    else if (o == "-L" && more) tight = atoi(argv[++a]);
    else if (o == "-r" && more) refresh = atoi(argv[++a]);
    else if (o == "-e" && more) eig_mode = atoi(argv[++a]);
    else if (o == "-g" && more) gamma = atof(argv[++a]);
    else if (o == "-G" && more) gamma1 = atof(argv[++a]);
    else if (o == "-T" && more) threads = atoi(argv[++a]);
    else if (o == "-phase" && more) phase = argv[++a];
    else if (o == "-cont" && more) cont = argv[++a];
    else if (o == "-pool" && more) pool = argv[++a];
    else if (o == "-grid" && more) sscanf(argv[++a], "%d:%d", &cp.grid_rows, &cp.grid_cols);
    else if (o == "-side" && more) cp.side = atoi(argv[++a]);
    else if (o == "-k" && more) cp.kcuts = atoi(argv[++a]);
    else if (o == "-I" && more) cp.first_iters = atoi(argv[++a]);
    else if (o == "-i" && more) { cp.iters = atoi(argv[++a]); pool_iters = cp.iters; }
    else if (o == "-rounds" && more) { cp.rounds = atoi(argv[++a]); pool_rounds = cp.rounds; }
    else if (o == "-delta" && more) { cp.delta = atof(argv[++a]); pool_delta = cp.delta; }
    else if (o == "-t" && more) { cp.theta = atof(argv[++a]); pool_theta = cp.theta; }
    else if (o == "-mu0" && more) pool_mu0 = atof(argv[++a]);
    else if (o == "-nocap") cp.cap_at_cutoff = 0;
    else if (o == "-freeze" && more) { cp.freeze_iters = atoi(argv[++a]); pool_freeze = cp.freeze_iters; }
    else if (o == "-warmY") cp.warm_primal = 1;
    else if (o == "-iters" && more) iters_cont = atoi(argv[++a]);
    else if (o == "-extra" && more) extras.push_back(argv[++a]);
    else if (o == "-cuts" && more) cutfiles.push_back(argv[++a]);
    else if (o == "-sep" && more) sepcmd = argv[++a];
    else if (o == "-notri") notri = true;
    else if (o == "-regions" && more) regfile = argv[++a];
    else if (o == "-nogrid") g_fdrv_regions.nogrid = 1;
    else if (o == "-C") childb = 1;
    else if (o == "-wy") wy = true;
    else if (o == "-dumpY" && more) dumpy = argv[++a];
    else if (o == "-save" && more) savep = argv[++a];
    else if (o[0] == '-') { fprintf(stderr, "bad option %s\n", o.c_str()); return 1; }
    else file = argv[a];
  }
  if (!file || (phase.empty() + cont.empty() + pool.empty()) != 2) {
    fprintf(stderr, "usage: see the comment at the top of fdrv.C\n"); return 1; }
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
  QAP qr;
  reduceQAP(q, &qr, node);
  const int m = n - node.nfix;
  const size_t n1 = (size_t)m * m + 1, N = n1 * n1;
  DNNParams base;
  base.eig_mode = eig_mode; base.cheb_degree = degree; base.tight_lam = tight; base.refresh = refresh;
  base.gamma = gamma; base.gamma1 = gamma1;
  printf("instance %s n=%d node %s m=%d cutoff %.4f\n", file, n, fixes.c_str(), m, cutoff);

  for (int a = 0; a < n; a++) {
    if (node.p[a] == NONE) g_fdrv_regions.free_orig[0].push_back(a);
    if (node.pinv[a] == NONE) g_fdrv_regions.free_orig[1].push_back(a);
  }
  if (!regfile.empty()) {
    // lines: side name i1 i2 ... (original indices on that side)
    FILE *f = fopen(regfile.c_str(), "r");
    if (!f) { fprintf(stderr, "cannot read %s\n", regfile.c_str()); return 1; }
    char line[65536];
    while (fgets(line, sizeof line, f)) {
      if (line[0] == '#' || line[0] == '\n') continue;
      FdrvRegion e; char name[256]; int off = 0;
      if (sscanf(line, "%d %255s %n", &e.side, name, &off) < 2) continue;
      e.name = name;
      for (char *c = line + off;;) { char *end; long v = strtol(c, &end, 10); if (end == c) break; e.orig.push_back((int)v); c = end; }
      g_fdrv_regions.extra.push_back(e);
    }
    fclose(f);
    printf("extra regions from %s: %zu%s\n", regfile.c_str(), g_fdrv_regions.extra.size(), g_fdrv_regions.nogrid ? " (grid regions off)" : "");
  }
  if (!phase.empty()) {
    std::vector<int> free_orig;
    for (int a = 0; a < n; a++)
      if ((cp.side == 0 ? node.p[a] : node.pinv[a]) == NONE) free_orig.push_back(a);
    cp.verbose = getenv("FDRV_VERBOSE") ? atoi(getenv("FDRV_VERBOSE")) : 1;
    if (cp.freeze_iters <= 0) cp.freeze_iters = 1500;   // the hook sits at the quiet round
    CutExport ex; ex.stop_before_quiet = 1;
    DNNWarm warm;
    DNNCutResult res = dnn_cut_rounds_f(qr.A, qr.B, qr.C, qr.shift, cutoff, base, free_orig, n, cp, &warm, -HUGE_VAL, &ex);
    printf("PHASE best_lb %.6f first_lb %.6f gain %.6f rounds %d cpu %.1f fathomed %d\n", res.best_lb, res.first_lb,
           res.best_lb - res.first_lb, res.rounds, res.cpu, res.status == DNN_FATHOMED);
    bool ok = write_vec(phase + ".Y0.bin", ex.Y0);
    if (ex.have) {
      Meta mt; mt.m = m; mt.scale = ex.state.scale; mt.hsum = ex.hsum_tri; mt.abs_ = ex.abs_tri;
      mt.active = ex.active; mt.best = ex.best_mu_lb; mt.first = res.first_lb;
      ok = ok && write_vec(phase + ".tri.bin", ex.Ltri) && write_vec(phase + ".sY.bin", ex.state.Y) &&
           write_vec(phase + ".sZ.bin", ex.state.Z) && write_meta(phase, mt);
      printf("saved %s: tri term %d active cuts, hsum %.6f, |L| %.3f, state scale %.6f, best round lb %.6f\n",
             phase.c_str(), ex.active, ex.hsum_tri, ex.abs_tri, ex.state.scale, ex.best_mu_lb);
    }
    printf("FINAL phase %s ok %d\n", phase.c_str(), (int)ok);
    delete q;
    return 0;
  }

  // -cont / -pool: load the saved phase
  const std::string P = cont.empty() ? pool : cont;
  Meta mt;
  if (!read_meta(P, mt) || mt.m != m) { fprintf(stderr, "cannot read %s.meta for m=%d\n", P.c_str(), m); return 1; }
  DNNWarm state; state.m = m; state.scale = mt.scale;
  std::vector<double> Lbase(N, 0.0);
  double hconst = 0, absL = 0;
  if (!read_vec(P + ".sY.bin", state.Y, N) || !read_vec(P + ".sZ.bin", state.Z, N)) {
    fprintf(stderr, "cannot read the state of %s\n", P.c_str()); return 1; }
  if (!notri) {
    if (!read_vec(P + ".tri.bin", Lbase, N)) { fprintf(stderr, "cannot read %s.tri.bin\n", P.c_str()); return 1; }
    hconst += mt.hsum; absL += mt.abs_;
  }
  for (const std::string &spec : extras) {
    size_t c = spec.rfind(':');
    std::string path = spec; double k = 0;
    if (c != std::string::npos && c > 0 && spec.find('/', c) == std::string::npos) { k = atof(spec.substr(c + 1).c_str()); path = spec.substr(0, c); }
    std::vector<double> E;
    if (!read_vec(path, E, N)) { fprintf(stderr, "cannot read %s\n", path.c_str()); return 1; }
    double a_ = 0; size_t nz = 0;
    for (size_t t = 0; t < N; t++) { Lbase[t] += E[t]; a_ += fabs(E[t]); nz += E[t] != 0; }
    absL += a_; hconst += k;
    printf("extra term %s: %zu nonzeros, |E| %.6f, const %.6f\n", path.c_str(), nz, a_, k);
  }
  printf("base term: tri %s (hsum %.6f, %d cuts, phase best %.4f, first %.4f), total const %.6f\n",
         notri ? "off" : "on", mt.hsum, mt.active, mt.best, mt.first, hconst);
  DNNParams par = base;
  par.hopeless_after = 0; par.hopeless_gap_ratio = 0;
  par.warm_primal = wy ? 1 : 0;

  if (!cont.empty()) {
    par.max_iter = iters_cont;
    par.L_add = Lbase.data(); par.L_add_abs = absL;
    par.child_bounds = childb;
    DNNWarm save;
    DNNResult r = dnn_bound(qr.A, qr.B, qr.C, qr.shift - hconst, cutoff, par, &state, &save);
    printf("RESULT lb=%.6f primal=%.6f status=%s iters=%d cpu=%.1f gain_over_phase_first=%.4f\n", r.lb, r.primal,
           dnn_status_name(r.status), r.iterations, r.cpu, r.lb - mt.first);
    if (childb && !r.child_lb.empty()) {
      // child bounds from the certificate (DNNResult::child_lb): lb + reduced cost
      std::vector<double> ex;
      int fath = 0;
      for (double c : r.child_lb) if (std::isfinite(c)) { ex.push_back(c - r.lb); fath += c > cutoff; }
      std::sort(ex.begin(), ex.end());
      auto q = [&](double f) { return ex.empty() ? 0.0 : ex[(size_t)(f * (ex.size() - 1))]; };
      printf("CHILD %zu children: fathomed by the certificate %d; excess over lb: median %.3f, 90%% %.3f, 99%% %.3f, max %.3f; children within 1 of lb %d\n",
             ex.size(), fath, q(0.5), q(0.9), q(0.99), ex.empty() ? 0.0 : ex.back(),
             (int)std::count_if(ex.begin(), ex.end(), [](double v) { return v < 1.0; }));
      if (!dumpy.empty()) write_vec(dumpy + ".child.bin", r.child_lb);
    }
    if (!dumpy.empty()) write_vec(dumpy + ".Y.bin", save.Y);
    if (!savep.empty()) {
      Meta m2 = mt; m2.scale = save.scale;
      write_vec(savep + ".sY.bin", save.Y); write_vec(savep + ".sZ.bin", save.Z); write_meta(savep, m2);
      if (!notri) write_vec(savep + ".tri.bin", std::vector<double>(Lbase));
    }
    printf("FINAL lb %.6f\n", r.lb);
    delete q;
    return 0;
  }

  // -pool: Polyak on external cuts on top of the base term
  std::vector<SCut> cuts;
  std::set<int64_t> keys;
  auto load = [&](const std::string &f) {
    std::vector<SCut> nc;
    if (!read_cuts(f, nc, (int)n1)) { fprintf(stderr, "cannot read cuts %s\n", f.c_str()); exit(1); }
    int added = 0;
    for (auto &s : nc) if (keys.insert(s.key).second) { s.mu = pool_mu0; cuts.push_back(s); added++; }
    return added;
  };
  for (auto &f : cutfiles) printf("cuts %s: %d loaded\n", f.c_str(), load(f));
  std::vector<double> L(N);
  double best = -HUGE_VAL, prev = -HUGE_VAL, theta = pool_theta, first = -HUGE_VAL;
  std::vector<double> best_mu;
  DNNWarm best_state = state, cur = state;
  auto solve = [&](int iters, const std::vector<double> &mu, bool first_round) {
    std::copy(Lbase.begin(), Lbase.end(), L.begin());
    double a_ = absL, h = hconst;
    int act = 0;
    for (size_t t = 0; t < cuts.size(); t++)
      if (mu[t] > 0) { a_ += add_scut(cuts[t], L, (int)n1, mu[t]); h += mu[t] * cuts[t].h; act++; }
    DNNParams p2 = par;
    p2.max_iter = iters; p2.L_add = L.data(); p2.L_add_abs = a_;
    p2.warm_primal = first_round ? par.warm_primal : 1;
    DNNWarm save;
    DNNResult r = dnn_bound(qr.A, qr.B, qr.C, qr.shift - h, cutoff, p2, &cur, &save);
    cur = save;
    return std::make_pair(r, act);
  };
  for (int rd = 0; rd <= pool_rounds; rd++) {
    std::vector<double> mu(cuts.size());
    for (size_t t = 0; t < cuts.size(); t++) mu[t] = cuts[t].mu;
    auto rr = solve(rd == 0 ? iters_cont : pool_iters, mu, rd == 0);
    DNNResult &r = rr.first;
    if (rd == 0) first = r.lb;
    if (r.lb > best) { best = r.lb; best_mu = mu; best_state = cur; }
    // separation on the new iterate
    int added = 0;
    if (!sepcmd.empty() && rd < pool_rounds) {
      std::string yp = pool + ".poolY.bin", cpth = pool + ".poolsep.cb", cmd = sepcmd;
      write_vec(yp, cur.Y);
      for (size_t at; (at = cmd.find("{Y}")) != std::string::npos;) cmd.replace(at, 3, yp);
      for (size_t at; (at = cmd.find("{C}")) != std::string::npos;) cmd.replace(at, 3, cpth);
      remove(cpth.c_str());
      if (system(cmd.c_str()) == 0) added = load(cpth);
    }
    double s2 = 0, vmax = 0; int nviol = 0;
    std::vector<double> sg(cuts.size());
    for (size_t t = 0; t < cuts.size(); t++) {
      sg[t] = cut_value(cuts[t], cur.Y, (int)n1);
      if (sg[t] > 1e-6) { nviol++; vmax = std::max(vmax, sg[t]); }
      if (cuts[t].mu > 0 || sg[t] > 0) s2 += sg[t] * sg[t];
    }
    printf("pool round %2d  lb %.4f (gain %.4f, best %.4f)  primal %.4f  status %s  iters %d  cpu %.1f  active %d  cuts %zu  violated %d (max %.4f)  added %d\n",
           rd, r.lb, r.lb - first, best - first, r.primal, dnn_status_name(r.status), r.iterations, r.cpu, rr.second,
           cuts.size(), nviol, vmax, added);
    if (r.status == DNN_FATHOMED || rd == pool_rounds || s2 <= 0) break;
    if (rd > 0 && r.lb < prev - 1e-9) theta *= 0.5;
    prev = r.lb;
    const double target = best + pool_delta;
    const double alpha = theta * std::max(target - r.lb, 1e-6) / s2;
    for (size_t t = 0; t < cuts.size(); t++) cuts[t].mu = std::max(0.0, cuts[t].mu + alpha * sg[t]);
  }
  double final_lb = best;
  if (pool_freeze > 0 && best <= cutoff) {
    best_mu.resize(cuts.size(), 0.0);
    cur = best_state;
    auto rr = solve(pool_freeze, best_mu, false);
    printf("pool quiet  lb %.4f (gain %.4f)  primal %.4f  iters %d  cpu %.1f  active %d\n", rr.first.lb,
           rr.first.lb - first, rr.first.primal, rr.first.iterations, rr.first.cpu, rr.second);
    final_lb = std::max(final_lb, rr.first.lb);
    if (!dumpy.empty()) write_vec(dumpy + ".Y.bin", cur.Y);
  }
  printf("FINAL pool best_lb %.6f first %.6f gain %.6f\n", final_lb, first, final_lb - first);
  delete q;
  return 0;
}
