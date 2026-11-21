#include "estimate.h"
#include "util.h"

#include <cmath>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <string>
#include <ctime>
#include <chrono>
#include "settings.h"
#include "state.h"
#include <cstdio>
#include <unistd.h>
#ifdef HAVE_DNN
#include "dnn.h"
#include "dnncut.h"
#endif

// One line per printed dive: depth, node estimate, time estimate.
// QAP_EST_HIST overrides the file name (default kept for old scripts).
double est_resume_root_bound = NAN;   // set by the driver for -e -r
static const char *dive_hist_file() { return setting_s("QAP_EST_HIST"); }

void standard_errs(double n_mean,double n_sumsq,
		   double t_mean,double t_sumsq,
		   double trials);


void print_dive(std::ostream &outs,
		QAPStatistics *top_stats,
		QAPStatistics *path_stats,
		double time_scale)
{
  // terminal depth, node estimate, and time estimate for each dive

  int depth=0;
  // 64-bit, as the counters are (tai50b dives reach 1e10 nodes)
  QAPStatistics::ctr_type node;
  double time;

  node = top_stats->node + path_stats->node;
  time = top_stats->total_time + path_stats->total_time;

  // calculate depth of dive
  while(path_stats->n_level[depth++]==0); 
  while(path_stats->n_level[depth++]>0);
  depth -= 2;

  //  path_stats->print();

  outs << depth << " ";
  outs << node << " ";
  outs << time_scale*time;

  outs << std::endl;

}

/*
 Survivor dives (QAP_EST_SURV=b; docs/estimates_plan.md, experiments/estimator_replay):

 Every node of the start stack (a state file's open nodes with -e -r F) is a stratum. Each is
 bounded once and so are its children, exactly: the nodes that do not branch are counted as they
 are. A dive from a stratum takes one of its children that branched (a survivor), with probability
 p = (1 - mix) w / sum w + mix / count, w = exp(b x), x = the survivor's shortfall (incumbent - its
 bound) over the median shortfall of the branched nodes this process has bounded at that depth
 (mix = est.xbar_mix; b = 0: uniform). It processes the survivor again to get its children, bounds
 all of them, counts the ones that do not branch with weight d = 1 / (p1 p2 ...), and goes on into
 one that branches, until none does. Each dive is an unbiased estimate of the stratum's subtree
 below its exact part (Purdom's partial backtracking; the dives never end on a child the search
 fathoms in one bound, as Knuth's do in a DNN tree: 93-99% of a depth's nodes are fathomed).

 Dives go to the strata by a pilot (QAP_EST_MIN_DIVES rounds whatever the budget, default 1, and a
 second while the budget lasts), then by s_h / sqrt(c_h) / n_h (Neyman at fixed cost, s_h floored
 at the stratum's mean; equal allocation for Knuth's dives), until QAP_EST_BUDGET seconds of wall
 time have passed since the strata's setup ended (checked between dives: a dive in progress
 finishes, so no expensive dive is cut) or, without a budget, the
 estimate.par trials. QAP_EST_DIVES=file: a JSON line per stratum (its exact part) and per dive
 (its estimate, by depth, its deepest node and its own CPU cost); coord's `qapctl estimate` adds
 the strata up. Node and CPU figures are the search's: a node's CPU is what bnb() spent on it.
*/
namespace {
struct Rec {                       // one node as bounded: its own statistics, how it ended
  QAPAssignment a;
  double cpu = 0.0;
  double bound = NAN;
  int nodes = 0;                   // 1, or 0 when it was discarded on its inherited bound
  bool branched = false;
  // the search's own counts at its depth (QAPStatistics f_level, l_level, e_level): fathomed (0 or 1),
  // locations its branching considered and eliminated, for the RESULTS table's %fathom and %c.elim
  double fathom = 0.0, locs = 0.0, elim = 0.0;
};

struct Est {                       // an estimate by depth: nodes and CPU; fathoms, locations, eliminated
  std::map<int, double> n, c, f, l, x;
  double nodes = 0.0, cpu = 0.0;
  void add(int d, double nn, double cc) { n[d] += nn; c[d] += cc; nodes += nn; cpu += cc; }
  void add(const Rec &r, double w) {   // a bounded node, weight w (1/p)
    const int d = r.a.depth;
    add(d, w * r.nodes, w * r.cpu);
    f[d] += w * r.fathom; l[d] += w * r.locs; x[d] += w * r.elim;
  }
};

static std::string fixes_of(const QAPAssignment &a)
{
  std::string f;
  char b[32];
  for (int k = 0; k < a.nfix; k++) {
    snprintf(b, sizeof b, "%s%d:%d", k ? "," : "", (int)a.history[k], (int)a.p[a.history[k]]);
    f += b;
  }
  return f;
}

// bound one node with the search's own code: bnb() on it alone with max_depth = its depth, so its
// children come back unprocessed (set aside, not counted). Their bound is the node's, as bnb leaves it.
static Node *bound_one(QAP *q, const QAPAssignment &a, QAPSolution *s, BNBParameters *params, int nlevel,
                       double root_bound, Rec &r)
{
  Node *one = new Node;
  one->data = a;
  one->link = NULL;
  QAPStatistics st(s->n, nlevel);
  st.root_bound = root_bound;
  st.stack_size = 1;
  bnb(q, one, s, &st, 0, 0, params, nlevel, a.depth, false, -1);
  r.a = a;
  r.cpu = st.total_time;
  r.nodes = (int)st.node;
  r.branched = one != NULL;
  const int dd = std::min((int)a.depth, s->n);
  r.fathom = (double)st.f_level[dd];
  r.locs = (double)st.l_level[dd];
  r.elim = (double)st.e_level[dd];
  r.bound = one ? one->data.bound : NAN;
  return one;
}

static void free_list(Node *&l) { while (l) list_head_remove(l); }

struct Stratum {
  std::string key;
  int depth = 0;
  double pi = 1.0;                 // its probability of being sampled (QAP_EST_SAMPLE): weight 1/pi
  std::vector<Rec> surv;           // its children that branched (survivor dives)
  std::vector<QAPAssignment> kids; // its children, unbounded (Knuth dives), and their probabilities
  std::vector<double> kp;
  Est exact;
  double setup_cpu = 0.0;
  // its children's warm start (the stratum's own final state, shared by them) and where it is kept: in
  // memory ('m'), spilled to a file and read back around each dive ('s'), or dropped: cold ('c')
  WarmStart *ws = NULL;
  std::string spill;
  char warm = 'm';
  int dives = 0;
  double s1 = 0.0, s2 = 0.0, cost = 0.0;   // CPU estimate sum, sum of squares; cost per dive (EMA)
};

static std::map<int, std::vector<double>> short_seen;   // shortfalls of branched nodes, by depth

static double median_short(int d, const std::vector<Rec> &sib, double inc)
{
  std::vector<double> v;
  auto it = short_seen.find(d);
  if (it != short_seen.end() && it->second.size() >= 3) v = it->second;
  else for (const Rec &r : sib) v.push_back(std::max(inc - r.bound, 0.0));
  if (v.empty()) return 1.0;
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return std::max(v[v.size() / 2], 1e-9);
}

static size_t draw_survivor(const std::vector<Rec> &surv, double b, double mix, double inc, double &prob)
{
  const int d = surv[0].a.depth;
  const double med = median_short(d, surv, inc);
  std::vector<double> w(surv.size());
  double sw = 0.0;
  for (size_t k = 0; k < surv.size(); k++) {
    const double x = std::max(inc - surv[k].bound, 0.0) / med;
    w[k] = b == 0.0 ? 1.0 : exp(std::min(b * x, 50.0));
    sw += w[k];
  }
  std::vector<double> pr(surv.size());
  for (size_t k = 0; k < surv.size(); k++)
    pr[k] = sw > 0.0 ? (1.0 - mix) * w[k] / sw + mix / surv.size() : 1.0 / surv.size();
  double u = drand(), acc = 0.0;
  size_t k = 0;
  for (; k + 1 < surv.size(); k++) { acc += pr[k]; if (u <= acc) break; }
  prob = pr[k];
  return k;
}

// Knuth's step: probabilities of a node's children from their own QPB bounds, g^q (g the relative
// gap against the root's QPB bound, as QAP_EST_CHILD_QPB), (1 - mix) w / sum + mix / count
static std::vector<double> qpb_probs(QAP *q, const std::vector<QAPAssignment> &kids, BNBParameters *params,
                                     int nlevel, double root_bound, double inc, double qexp, double mix)
{
  std::vector<double> w(kids.size(), 0.0);
  double sw = 0.0;
  for (size_t k = 0; k < kids.size(); k++) {
    int level;
    double rel_gap;
    BNBParameters *pp = pick_param(kids[k], root_bound, inc, params, level, rel_gap, nlevel);
    const double bq = node_qpb_bound(q, kids[k], pp);
    const double g = bq < inc ? (inc - bq) / std::max(inc - root_bound, 1e-9) : 0.0;
    w[k] = g > 0.0 ? pow(g, qexp) : 0.0;
    sw += w[k];
  }
  for (double &x : w) x = sw > 0.0 ? (1.0 - mix) * x / sw + mix / kids.size() : 1.0 / kids.size();
  return w;
}

// QAP_EST_CHILD_DNN=CAP (experiments/estimates_nov19, "Predictor-weighted dives"): Knuth's dives with
// the thesis's design (ch. 4.3: children drawn by g^q of a cheap bound computed for the purpose) and a
// predictor fit for a DNN tree, where QPB's gaps barely differ between children (estimator_replay: QPB's
// shortfall adds nothing to a subtree's size). Each child is solved briefly, CAP ADMM iterations warm
// from its parent's state as its own first solve would start (the parent's Z and conditioned Y, under
// the parent's cut pool where it hands one on: dnn_pool_bound), no give-up, no cut phase; its shortfall
// (cutoff - lb, 0 when it fathoms) is g's numerator. In experiments/line_signals (tho30's three depth-5
// nodes, every child of 8 lines) 30 warm iterations correlate 0.92-0.98 with the cold 300-iteration bound
// and 0.80-0.90 with the child's subtree, at about 1 CPU s a child at m = 28. NAN: no parent DNN state
// (a QPB-only parent, or one past QAP_EST_WARM_KEEP), and the caller falls back to QPB's weights.
static double child_probe(QAP *q, const QAPAssignment &kid, double inc, int cap)
{
#ifdef HAVE_DNN
  if (kid.bound >= inc - DELTA) return 0.0;
  if (!kid.warm || !kid.warm->dnn) return NAN;
  const int n = kid.n, m = n - kid.nfix;
  const DNNWarm &pw = *(const DNNWarm *)kid.warm->dnn;
  int di, dj;
  if (m < 2 || pw.m != m + 1 || !warm_child_offsets(kid, di, dj)) return NAN;
  QAP qc(1);
  reduceQAP(q, &qc, kid);
  DNNParams dp;
  dp.gamma1 = setting_d("dnn.gamma1");
  dp.gamma = setting_d("dnn.gamma");
  dp.cheb_degree = (int)setting_i("dnn.cheb");
  dp.tight_lam = (int)setting_i("dnn.tight_lam");
  dp.momentum = setting_d("dnn.momentum");
  dp.beta = setting_d("dnn.beta");
  dp.refresh = (int)setting_i("dnn.refresh");
  dp.hopeless_after = 0;
  dp.hopeless_after_warm = 0;
  dp.max_iter = cap;
  static const int warm_y = (int)setting_i("dnn.warm_y");
  DNNWarm start;
  const DNNWarm *sp = NULL;
  if (warm_y > 0 && dnn_warm_child(pw, di, dj, start, warm_y) && !start.Y.empty()) {
    sp = &start;
    dp.warm_primal = 1;
  } else if (dnn_warm_child(pw, di, dj, start)) {
    sp = &start;
  }
  const double cutoff = inc - DELTA;
  static const int pool_on = (int)setting_i("pool.on");
  static const int cut_side = (int)setting_i("cut.side");
  DNNPool pool_in, pool_mid;
  DNNResult r;
  if (pool_on && !kid.warm->dnn_pool.empty() && dnn_pool_decode(kid.warm->dnn_pool, pool_in) &&
      pool_in.side == cut_side && !pool_in.cuts.empty()) {
    std::vector<int> free_orig, other_orig, fixed_at(n, -1);
    for (int a = 0; a < n; a++) {
      if ((cut_side == 0 ? kid.p[a] : kid.pinv[a]) == NONE) free_orig.push_back(a);
      if ((cut_side == 0 ? kid.pinv[a] : kid.p[a]) == NONE) other_orig.push_back(a);
      else fixed_at[a] = cut_side == 0 ? kid.pinv[a] : kid.p[a];
    }
    DNNCutParams pc;
    pc.side = cut_side; pc.pool_in = &pool_in; pc.pool_out = &pool_mid;
    pc.other_orig = &other_orig; pc.fixed_at = &fixed_at; pc.pool_check = 0;
    int kept = 0, dropped = 0;
    r = dnn_pool_bound(qc.A, qc.B, qc.C, qc.shift, cutoff, dp, free_orig, n, pc, sp, NULL, &kept, &dropped);
  } else
    r = dnn_bound(qc.A, qc.B, qc.C, qc.shift, cutoff, dp, sp, NULL);
  return r.status == DNN_FATHOMED ? 0.0 : std::max(0.0, cutoff - r.lb);
#else
  (void)q; (void)kid; (void)inc; (void)cap;
  return NAN;
#endif
}

static long child_probes = 0, child_probe_fallbacks = 0;

// a node's children's probabilities for a Knuth step: by their short DNN solves (QAP_EST_CHILD_DNN),
// g^q of the shortfall as qpb_probs takes QPB's, else (or with no parent state) by their QPB bounds
static std::vector<double> child_probs(QAP *q, const std::vector<QAPAssignment> &kids, BNBParameters *params,
                                       int nlevel, double root_bound, double inc, double qexp, double mix)
{
  static const int cap = (int)setting_i("QAP_EST_CHILD_DNN");
  if (cap > 0 && !kids.empty()) {
    std::vector<double> w(kids.size(), 0.0);
    bool ok = true;
    for (size_t k = 0; k < kids.size() && ok; k++) {
      const double sh = child_probe(q, kids[k], inc, cap);
      if (std::isnan(sh)) { ok = false; break; }
      const double g = sh / std::max(inc - root_bound, 1e-9);
      w[k] = g > 0.0 ? pow(g, qexp) : 0.0;
    }
    if (ok) {
      child_probes += kids.size();
      double sw = 0.0;
      for (double x : w) sw += x;
      for (double &x : w) x = sw > 0.0 ? (1.0 - mix) * x / sw + mix / kids.size() : 1.0 / kids.size();
      return w;
    }
    child_probe_fallbacks++;
  }
  return qpb_probs(q, kids, params, nlevel, root_bound, inc, qexp, mix);
}

static size_t draw_index(const std::vector<double> &p)
{
  double u = drand(), acc = 0.0;
  size_t k = 0;
  for (; k + 1 < p.size(); k++) { acc += p[k]; if (u <= acc) break; }
  return k;
}

static void json_by_depth(std::ostream &o, const char *name, const std::map<int, double> &m)
{
  o << "\"" << name << "\":{";
  bool first = true;
  for (auto &kv : m) {
    if (kv.second == 0.0) continue;
    o << (first ? "" : ",") << "\"" << kv.first << "\":" << kv.second;
    first = false;
  }
  o << "}";
}
}  // namespace

// the spilled strata's states (Stratum::spill): files beside QAP_EST_DIVES (else in TMPDIR), removed at
// the end; QAP_EST_SPILL_MB in all at most
static std::vector<std::string> spill_files;
static double spill_bytes = 0.0;
static long spill_reads = 0;

static bool spill_save(WarmStart *w, Stratum &h, size_t index)
{
  static const double cap = setting_d("QAP_EST_SPILL_MB") * 1e6;
  const std::string tok = dnn_state_encode(w, 32);
  if (tok == "-" || spill_bytes + tok.size() > cap) return false;
  std::string dir = setting_s("QAP_EST_DIVES");
  const size_t slash = dir.rfind('/');
  dir = slash == std::string::npos ? std::string(".") : dir.substr(0, slash);
  if (!*setting_s("QAP_EST_DIVES")) dir = getenv("TMPDIR") && *getenv("TMPDIR") ? getenv("TMPDIR") : "/tmp";
  const std::string path = dir + "/qap-est-spill-" + std::to_string((long)getpid()) + "-" + std::to_string(index) + ".dnn";
  FILE *f = fopen(path.c_str(), "wb");
  if (!f) return false;
  const bool ok = fwrite(tok.data(), 1, tok.size(), f) == tok.size();
  if (fclose(f) != 0 || !ok) { unlink(path.c_str()); return false; }
  spill_files.push_back(path);
  spill_bytes += tok.size();
  h.spill = path;
  w->dnn_free(w->dnn);
  w->dnn = NULL;
  return true;
}

static void spill_load(Stratum &h)
{
  if (h.spill.empty() || !h.ws || h.ws->dnn) return;
  FILE *f = fopen(h.spill.c_str(), "rb");
  if (!f) return;
  std::string tok;
  char buf[1 << 16];
  size_t k;
  while ((k = fread(buf, 1, sizeof buf, f)) > 0) tok.append(buf, k);
  fclose(f);
  void *d = dnn_state_decode(tok);
  if (d) { h.ws->dnn = d; h.ws->dnn_free = dnn_state_free; spill_reads++; }
}

static void spill_drop(Stratum &h)
{
  if (h.spill.empty() || !h.ws || !h.ws->dnn) return;
  h.ws->dnn_free(h.ws->dnn);
  h.ws->dnn = NULL;
}

static QAPStatistics *estimate_surv(QAP *q, Node *&stack, QAPSolution *s, BNBParameters *params, int nlevel,
                                    EstimateParameters *ep)
{
  // QAP_EST_SURV=b: survivor dives, weights exp(b x); =kQ: Knuth's dives inside each stratum (one
  // child a level, g^Q of its own QPB bound), for bushy trees where bounding every child costs the tree
  const char *sv = setting_s("QAP_EST_SURV");
  const bool knuth = sv[0] == 'k';
  const double b = atof(knuth ? sv + 1 : sv);
  const double mix = setting_d("est.xbar_mix");
  const double sample = setting_d("QAP_EST_SAMPLE");
  const bool sample_qpb = std::string(setting_s("QAP_EST_SAMPLE_W")) == "qpb";
  const double budget = setting_d("QAP_EST_BUDGET");
  const int min_dives = setting_i("QAP_EST_MIN_DIVES");
  // QAP_EST_SAMPLE_DEPTHS=D:F,...: a sampling rate by depth (experiments/tho40_est_diag: the open nodes at
  // the depth that holds most of the tree, all of them; the deep ones, which hold little, a few percent)
  std::map<int, double> depth_rate;
  for (std::string spec = setting_s("QAP_EST_SAMPLE_DEPTHS"); !spec.empty();) {
    const size_t comma = spec.find(',');
    const std::string item = spec.substr(0, comma);
    const size_t colon = item.find(':');
    if (colon == std::string::npos) { std::cerr << "QAP_EST_SAMPLE_DEPTHS: " << item << " is not D:F" << std::endl; exit(2); }
    depth_rate[atoi(item.substr(0, colon).c_str())] = std::min(1.0, std::max(0.0, atof(item.substr(colon + 1).c_str())));
    spec = comma == std::string::npos ? "" : spec.substr(comma + 1);
  }
  const char *dives_file = setting_s("QAP_EST_DIVES");
  if (const char *e = setting_s("QAP_EST_SEED"); *e) srand((unsigned)strtoul(e, NULL, 10));
  // the relative gaps that pick each node's strategy are against the root's QPB bound: a state
  // file's (-e -r), else computed here as a -F run does (driver.C)
  const double root_bound = !std::isnan(est_resume_root_bound) ? est_resume_root_bound : root_qpb_bound(q, params);
  const double inc = s->getObjective();
  std::ofstream out;
  if (*dives_file) out.open(dives_file, std::ofstream::app);
  out.precision(10);
  std::cout << (knuth ? "Knuth dives in strata: q = " : "survivor dives: b = ") << b
            << (knuth && setting_i("QAP_EST_CHILD_DNN") > 0
                ? ", children weighted by a " + std::to_string(setting_i("QAP_EST_CHILD_DNN")) + "-iteration warm DNN solve"
                : std::string()) << ", mix = " << mix
            << ", budget " << budget << " s, sample " << sample << (sample_qpb ? " (by QPB)" : "") << ", root bound "
            << root_bound << ", incumbent " << inc << std::endl;

  // QAP_EST_SAMPLE=F: each start node is a stratum with probability pi = min(1, F n w / sum w)
  // (w = 1, or g^2 of its own QPB bound with QAP_EST_SAMPLE_W=qpb), the rest are not bounded;
  // a stratum's estimate is weighted 1/pi (Horvitz-Thompson: unbiased). F = 1 keeps every node.
  std::vector<QAPAssignment> starts;
  while (stack != NULL) { starts.push_back(stack->data); list_head_remove(stack); }
  std::vector<double> pis(starts.size(), 1.0);
  if (sample < 1.0 && !starts.empty()) {
    std::vector<double> w(starts.size(), 1.0);
    std::vector<bool> capped(starts.size(), false);
    if (sample_qpb) {
      // a node its own QPB bound fathoms is kept (pi = 1: it costs one QPB bound); the rest weighted
      // g^2, with est.xbar_mix's uniform share mixed in so that no pi is tiny (a near-zero pi is
      // unbiased in expectation and low in practice: nug15 read 0.91 without it, Nov 18)
      double sw = 0.0;
      int nfree = 0;
      for (size_t k = 0; k < starts.size(); k++) {
        int level;
        double rel_gap;
        BNBParameters *pp = pick_param(starts[k], root_bound, inc, params, level, rel_gap, nlevel);
        const double bq = node_qpb_bound(q, starts[k], pp);
        if (bq > inc - DELTA) { capped[k] = true; w[k] = 0.0; continue; }
        const double g = (inc - bq) / std::max(inc - root_bound, 1e-9);
        w[k] = g * g;
        sw += w[k];
        nfree++;
      }
      for (size_t k = 0; k < starts.size(); k++)
        if (!capped[k]) w[k] = sw > 0.0 ? (1.0 - mix) * w[k] / sw + mix / nfree : 1.0;
    }
    int nqpb = 0;
    for (size_t k = 0; k < starts.size(); k++) nqpb += capped[k];
    const double m = std::max(1.0, sample * (starts.size() - nqpb)) + nqpb;
    for (int pass = 0; pass < 10; pass++) {     // inclusion probabilities, the capped ones' excess spread
      double sw = 0.0, mm = m;
      for (size_t k = 0; k < starts.size(); k++) { if (capped[k]) mm -= 1.0; else sw += w[k]; }
      bool more = false;
      for (size_t k = 0; k < starts.size(); k++) {
        if (capped[k]) { pis[k] = 1.0; continue; }
        pis[k] = sw > 0.0 && mm > 0.0 ? mm * w[k] / sw : 0.0;
        if (pis[k] >= 1.0) { capped[k] = true; more = true; }
      }
      if (!more) break;
    }
  }
  for (size_t k = 0; k < starts.size(); k++) {      // a depth's own rate wins over the share
    auto it = depth_rate.find(starts[k].depth);
    if (it != depth_rate.end()) pis[k] = it->second;
  }
  int kept = 0;
  for (size_t k = 0; k < starts.size(); k++) {
    if (pis[k] < 1.0 && !(drand() < pis[k])) { pis[k] = 0.0; continue; }
    kept++;
  }
  for (size_t k = starts.size(); k-- > 0;)
    if (pis[k] > 0.0) list_head_insert(stack, starts[k]);
  if (out.is_open())
    out << "{\"type\":\"sample\",\"open\":" << starts.size() << ",\"kept\":" << kept << ",\"fraction\":" << sample
        << ",\"weights\":\"" << (sample_qpb ? "qpb" : "unif") << "\",\"depths\":\"" << setting_s("QAP_EST_SAMPLE_DEPTHS")
        << "\"}" << std::endl;
  std::vector<double> kept_pi;
  for (size_t k = 0; k < starts.size(); k++) if (pis[k] > 0.0) kept_pi.push_back(pis[k]);
  size_t next_pi = 0;

  // the strata: each start node bounded, and its children
  std::vector<Stratum> strata;
  Est exact_all;
  double spent = 0.0;
  while (stack != NULL) {
    QAPAssignment a = stack->data;
    list_head_remove(stack);
    Stratum h;
    h.key = fixes_of(a);
    h.depth = a.depth;
    h.pi = next_pi < kept_pi.size() ? kept_pi[next_pi++] : 1.0;
    Rec r;
    Node *kids = bound_one(q, a, s, params, nlevel, root_bound, r);
    h.exact.add(r, 1.0);
    h.setup_cpu += r.cpu;
    if (knuth && kids && r.branched) {
      const double c0 = clock();
      for (Node *k = kids; k != NULL; k = k->link) h.kids.push_back(k->data);
      h.kp = child_probs(q, h.kids, params, nlevel, root_bound, inc, b, mix);
      h.setup_cpu += (clock() - c0) / CLOCKS_PER_SEC;
    } else if (kids && r.branched) {
      short_seen[a.depth].push_back(std::max(inc - r.bound, 0.0));
      for (Node *k = kids; k != NULL; k = k->link) {
        Rec c;
        Node *gk = bound_one(q, k->data, s, params, nlevel, root_bound, c);
        h.setup_cpu += c.cpu;
        if (c.branched) {
          short_seen[c.a.depth].push_back(std::max(inc - c.bound, 0.0));
          h.surv.push_back(c);     // counted by the dives (weight 1/p), not here
        } else
          h.exact.add(c, 1.0);
        free_list(gk);
      }
    }
    free_list(kids);
    // A survivor, or a Knuth child, holds a reference to its parent's (the stratum's) warm start, for
    // the dives: the parent's DNN state, its cut pool and X. Past QAP_EST_WARM_KEEP strata the DNN state
    // (tens of MB at n = 60) goes to a file, fp32 as a state file ships it (ship.bits 32), and is read
    // back around each dive from the stratum (spill_load). Until Nov 19 it was dropped, with the pool:
    // the first bound of a Knuth dive from such a stratum ran cold, without the inherited cuts, and the
    // estimate counted it as the tree's: in tho40-est2-261119, 44% of the strata, whose first children
    // cost 2.1-4.6x the CPU of the others' and branched more (experiments/estimates_nov19). Past
    // QAP_EST_SPILL_MB of files in all the old way stands, and the stratum's line says so ("warm":"c").
    static const long warm_keep = setting_i("QAP_EST_WARM_KEEP");
    if ((long)strata.size() >= warm_keep) {
      WarmStart *w = !h.kids.empty() ? h.kids[0].warm : (!h.surv.empty() ? h.surv[0].a.warm : NULL);
      bool shared = w != NULL;
      for (Rec &c : h.surv) shared = shared && c.a.warm == w;
      for (QAPAssignment &c : h.kids) shared = shared && c.warm == w;
      if (shared && w->dnn && spill_save(w, h, strata.size())) {
        h.ws = w;
        h.warm = 's';
      } else if (w) {
        for (Rec &c : h.surv) c.a.set_warm(NULL);
        for (QAPAssignment &c : h.kids) c.set_warm(NULL);
        h.warm = 'c';
      }
    }
    spent += h.setup_cpu;
    if (out.is_open()) {
      out << "{\"type\":\"stratum\",\"key\":\"" << h.key << "\",\"depth\":" << h.depth << ",\"exact_nodes\":"
          << h.exact.nodes << ",\"exact_cpu\":" << h.exact.cpu << ",";
      json_by_depth(out, "exact_by_depth", h.exact.n);
      out << ",";
      json_by_depth(out, "exact_cpu_by_depth", h.exact.c);
      out << ",";
      json_by_depth(out, "exact_fathom_by_depth", h.exact.f);
      out << ",";
      json_by_depth(out, "exact_locs_by_depth", h.exact.l);
      out << ",";
      json_by_depth(out, "exact_elim_by_depth", h.exact.x);
      out << ",\"survivors\":" << (knuth ? h.kids.size() : h.surv.size()) << ",\"setup_cpu\":" << h.setup_cpu
          << ",\"pi\":" << h.pi << ",\"mode\":\"" << (knuth ? "knuth" : "surv") << "\",\"warm\":\"" << h.warm
          << "\"}" << std::endl;
    }
    strata.push_back(h);
  }

  // one Knuth dive from stratum h: one child a level, by its own QPB bound
  auto knuth_dive = [&](Stratum &h, Est &e, int &deepest, double &cost) {
    deepest = h.depth + 1;
    cost = 0.0;
    size_t k = draw_index(h.kp);
    QAPAssignment cur = h.kids[k];
    double d = 1.0 / h.kp[k];
    while (true) {
      Rec r;
      Node *kids = bound_one(q, cur, s, params, nlevel, root_bound, r);
      cost += r.cpu;
      e.add(r, d);
      deepest = std::max(deepest, (int)cur.depth);
      if (!kids || !r.branched) { free_list(kids); break; }
      std::vector<QAPAssignment> ks;
      for (Node *t = kids; t != NULL; t = t->link) ks.push_back(t->data);
      free_list(kids);
      const double c0 = clock();
      std::vector<double> kp = child_probs(q, ks, params, nlevel, root_bound, inc, b, mix);
      cost += (clock() - c0) / CLOCKS_PER_SEC;
      k = draw_index(kp);
      cur = ks[k];
      d /= kp[k];
    }
  };

  // one dive from stratum h: its estimate below the stratum's exact part
  auto one_dive = [&](Stratum &h, Est &e, int &deepest, double &cost) {
    if (knuth) { knuth_dive(h, e, deepest, cost); return; }
    deepest = h.depth + 1;
    cost = 0.0;
    double pk = 1.0;
    size_t k = draw_survivor(h.surv, b, mix, inc, pk);
    Rec cur = h.surv[k];
    double d = 1.0 / pk;
    while (true) {
      e.add(cur, d);
      deepest = std::max(deepest, cur.a.depth);
      Rec again;
      Node *kids = bound_one(q, cur.a, s, params, nlevel, root_bound, again);   // its children, again
      cost += again.cpu;
      std::vector<Rec> surv;
      for (Node *t = kids; t != NULL; t = t->link) {
        Rec c;
        Node *gk = bound_one(q, t->data, s, params, nlevel, root_bound, c);
        cost += c.cpu;
        deepest = std::max(deepest, c.a.depth);
        if (c.branched) {
          short_seen[c.a.depth].push_back(std::max(inc - c.bound, 0.0));
          surv.push_back(c);
        } else
          e.add(c, d);
        free_list(gk);
      }
      free_list(kids);
      if (surv.empty()) break;
      k = draw_survivor(surv, b, mix, inc, pk);
      cur = surv[k];
      d /= pk;
    }
  };

  std::vector<size_t> live;
  for (size_t i = 0; i < strata.size(); i++) if (!strata[i].surv.empty() || !strata[i].kids.empty()) live.push_back(i);
  int total = 0;
  // the budget is the dives': it starts here, after the strata's setup (bounding the sampled start nodes
  // and their children took 68 of tho40-est's 95 estimator CPU hours, and left 828 of its 1,108 dived
  // strata one dive each: experiments/tho40_est_diag)
  const auto t_dives = std::chrono::steady_clock::now();   // not time(): a budget under a second counts
  auto over = [&]() {
    if (budget > 0) return std::chrono::duration<double>(std::chrono::steady_clock::now() - t_dives).count() >= budget;
    return total >= ep->trials;
  };
  auto run = [&](size_t i) {
    Stratum &h = strata[i];
    Est e;
    int deepest;
    double cost;
    spill_load(h);                    // a spilled stratum's state, for this dive only
    one_dive(h, e, deepest, cost);
    spill_drop(h);
    h.dives++;
    h.s1 += e.cpu;
    h.s2 += e.cpu * e.cpu;
    h.cost = h.dives == 1 ? cost : 0.7 * h.cost + 0.3 * cost;
    spent += cost;
    total++;
    if (out.is_open()) {
      out << "{\"type\":\"dive\",\"key\":\"" << h.key << "\",\"nodes\":" << e.nodes << ",\"cpu\":" << e.cpu << ",";
      json_by_depth(out, "by_depth", e.n);
      out << ",";
      json_by_depth(out, "cpu_by_depth", e.c);
      out << ",";
      json_by_depth(out, "fathom_by_depth", e.f);
      out << ",";
      json_by_depth(out, "locs_by_depth", e.l);
      out << ",";
      json_by_depth(out, "elim_by_depth", e.x);
      out << ",\"deepest\":" << deepest << ",\"cost\":" << cost << "}" << std::endl;
    }
  };
  // the pilot: QAP_EST_MIN_DIVES rounds whatever the budget (default 1), and a second while it lasts
  for (int pilot = 0; pilot < std::max(2, min_dives); pilot++)
    for (size_t i : live) {
      if (pilot >= min_dives && over()) break;
      run(i);
    }
  // QAP_EST_ALLOC=equal: the dives go round the strata in turn. The default, Neyman's allocation from
  // the dives so far, puts few dives where the first ones were small: with Knuth's heavy tails that
  // reads low (nug15: -6% with Neyman, Nov 18), so a fixed allocation is the unbiased choice there
  // (auto: equal for Knuth's dives, Neyman for survivor dives, whose estimates measured unbiased under it)
  static const std::string alloc = setting_s("QAP_EST_ALLOC");
  const bool equal_alloc = alloc == "equal" || (alloc == "auto" && knuth);
  size_t rr = 0;
  while (!live.empty() && !over() && equal_alloc) run(live[rr++ % live.size()]);
  while (!live.empty() && !over()) {
    double best = -1.0;
    size_t pick = live[0];
    for (size_t i : live) {
      const Stratum &h = strata[i];
      const double mu = h.s1 / h.dives;
      double sd = sqrt(std::max(h.s2 / h.dives - mu * mu, 0.0));
      sd = std::max(sd, mu);
      const double w = sd / h.pi / sqrt(std::max(h.cost, 1.0)) / h.dives;
      if (w > best) { best = w; pick = i; }
    }
    run(pick);
  }

  // the estimate of the start nodes' subtrees: each stratum's exact part plus its dives' mean
  QAPStatistics *res = new QAPStatistics(s->n, nlevel);
  double en = 0.0, ec = 0.0;
  for (const Stratum &h : strata) {
    en += h.exact.nodes / h.pi;
    ec += h.exact.cpu / h.pi;
    if (h.dives > 0) ec += h.s1 / h.dives / h.pi;
  }
  std::cout << "estimator: " << kept << " of " << starts.size() << " start nodes sampled, "
            << strata.size() << " strata (" << live.size() << " with children to dive), " << total
            << " dives, " << spent << " CPU s spent; estimate below the start nodes: exact part " << en
            << " nodes, CPU " << ec << " s (exact part and the dives' mean)" << std::endl;
  if (setting_i("QAP_EST_CHILD_DNN") > 0)
    std::cout << "estimator: " << child_probes << " children weighted by their short DNN solves, "
              << child_probe_fallbacks << " Knuth steps by QPB (no parent DNN state)" << std::endl;
  {
    int ns = 0, nc = 0;
    for (const Stratum &h : strata) { ns += h.warm == 's'; nc += h.warm == 'c'; }
    if (ns || nc)
      std::cout << "estimator: past QAP_EST_WARM_KEEP, " << ns << " strata's states spilled (" << spill_bytes / 1e6
                << " MB, read back " << spill_reads << " times), " << nc << " strata cold (QAP_EST_SPILL_MB)" << std::endl;
  }
  for (const std::string &f : spill_files) unlink(f.c_str());
  spill_files.clear();
  res->total_time = ec;
  return res;
}

QAPStatistics *estimate(QAP *q,Node *&stack,QAPSolution *s,
			BNBParameters *params,int nlevel,
			EstimateParameters *ep)
{
  // QAP_EST_SURV=b: survivor dives (estimate_surv below) in place of Knuth's
  if (*setting_s("QAP_EST_SURV")) return estimate_surv(q, stack, s, params, nlevel, ep);

  std::ofstream dive_hist;

  Node *p, *path;
  QAPStatistics *path_stats, *top_stats, *total_stats;
  QAPStatistics *estimated_stats;
  int too_deep = 0;
  double sum_of_bounds, sum_so_far, importance, p_score;
  double val, root_bound, opt, size;
  double n_predicted=0.0, t_predicted=0.0; // node, time sum squared
  double n_sumsq=0.0, t_sumsq=0.0; // node, time sum squared

  top_stats  = new QAPStatistics(s->n,nlevel);
  path_stats = new QAPStatistics(s->n,nlevel);
  total_stats = new QAPStatistics(s->n,nlevel);
  estimated_stats = NULL;

  opt = s->getObjective();

  // QAP_EST_SEED: seed rand() so that parallel estimator processes take
  // different dives; unset keeps the old fixed sequence.
  if (const char *e = setting_s("QAP_EST_SEED"); *e) {
    srand((unsigned)strtoul(e, NULL, 10));
    std::cout << "Estimator seed = " << e << std::endl;
  }
  std::cout << "Dive history file = " << dive_hist_file() << std::endl;

  params->score_exponent = ep->score_exponent;
  std::cout << "Scoring exponent = " << params->score_exponent << std::endl;

  // step 1: run branch and bound code on the given problem,
  // returning all nodes of level start_dive_level.
  // With -e -r F the nodes step 1 bounds (F's nodes shallower than start_dive_level) take their
  // strategy (pick_param) and the DNN gap filter from the root bound, which bnb sets only at a
  // depth-0 node: F's must be in place before the call, not after it (review finding 2, Nov 16;
  // before, those nodes were bounded with root bound 0, so every gap looked small and they got the
  // deepest level's strategy: nug15 from a depth-0 state, strategy 2 at all six depth-1 nodes where
  // the run uses 0 and 1, weaker bounds and other branching lines; tests/estimate_test.py checks it).
  if (!std::isnan(est_resume_root_bound)) top_stats->root_bound = est_resume_root_bound;

  bnb(q,stack,s,top_stats,0,0,params,nlevel,ep->start_dive_level-1,false,-1);

  root_bound = top_stats->root_bound;
  // tai50b_size: dives resumed from a state file (-e -r F): the top of the tree was
  // bounded by the run that wrote F, so the root bound is F's, not this run's (unset)
  if (!std::isnan(est_resume_root_bound)) {
    root_bound = est_resume_root_bound;
    top_stats->root_bound = root_bound;
    std::cout << "estimator resumed from a state file: root bound " << root_bound << std::endl;
  }

  top_stats->print();

  // yuck.
  int sdl = ep->start_dive_level;
  top_stats->node = top_stats->node - top_stats->n_level[sdl+1];
  top_stats->n_level[sdl+1] = 0;
  top_stats->f_level[sdl+1] = 0;
  top_stats->t_level[sdl+1] = 0;
  top_stats->l_level[sdl+1] = 0;
  top_stats->e_level[sdl+1] = 0;
  for(int k=0;k<top_stats->m;k++)
    top_stats->n_strat[sdl+1][k] = 0;

  // step 2: each node's score, and their sum. QAP_EST_CHILD_QPB=1 (as for a
  // dive's children, bnb.C): by each node's own QPB bound, (1-mix) g^q / sum +
  // mix / count with mix = est.xbar_mix, the nodes themselves unchanged;
  // outside the top's and the dives' clocks. Otherwise by the predicted bound.
  std::vector<double> first;
  for(p=stack;p!=NULL;p=p->link)
    first.push_back(score(opt,root_bound,p->data.predicted_bound,
			  params->score_exponent));
  if (setting_i("QAP_EST_CHILD_QPB") && first.size() > 1) {
    const double mix = setting_d("est.xbar_mix");
    int level, k = 0;
    double rel_gap, wsum = 0.0;
    for (p = stack; p != NULL; p = p->link, k++) {
      BNBParameters *pp = pick_param(p->data, root_bound, opt, params, level, rel_gap, nlevel);
      const double b = node_qpb_bound(q, p->data, pp);
      first[k] = b < opt ? score(opt, root_bound, b, params->score_exponent) : 0.0;
      if (!(first[k] > 0.0)) first[k] = 0.0;
      wsum += first[k];
    }
    for (double &w : first)
      w = wsum > 0.0 ? (1.0 - mix) * w / wsum + mix / first.size() : 1.0 / first.size();
    std::cout << "QAP_EST_CHILD_QPB: the first draw's " << first.size() << " nodes scored by their own QPB bound" << std::endl;
  }
  sum_of_bounds = 0.0;
  for (double w : first) sum_of_bounds += w;
  //  sum_of_bounds += 1.0;

  std::cout << "sum = " << sum_of_bounds << std::endl;
  std::cout << "stack = " << top_stats->stack_size << std::endl;

  size = double(top_stats->stack_size);
  if(top_stats->stack_size == 0) {
    total_stats->merge_stats(top_stats);
    delete path_stats;
    delete top_stats;
    return total_stats;
  } 

  for(int i=1;i<=ep->trials;i++) {
    path_stats->clear();

    // step 3: get a random number between 0 and sum_of_bounds.
    //         select a node for estimation based on this value.

    val = sum_of_bounds*drand();
    
    sum_so_far = 0.0; p = stack;
    size_t pk = 0;
    do {
      sum_so_far += first[pk];
      if(sum_so_far >= val) break;
      p = p->link; pk++;
    } while(p->link!=NULL);

    // p now points to a node we want to estimate.

    path = new Node;
    path->data = p->data;
    path->link = NULL;

    p_score = first[pk];
    importance  = p_score / sum_of_bounds;
    
    // step 4: random walk down path of tree.
    path_stats->root_bound = root_bound;
    params->est_factor = 1.0/importance;
    bnb(q,path,s,path_stats,0,0,params,nlevel,-1,true,-1);

    if(dive_depth(path_stats) >= ep->ignore_depth) {
      //      cout << "too deep: " << i << std::endl;
      //      path_stats->print();
      too_deep++;
      i--; continue;
    }
    
    n_predicted = top_stats->node + path_stats->node;
    t_predicted = ep->time_scale*(top_stats->total_time + 
				  path_stats->total_time);
    n_sumsq += n_predicted*n_predicted;
    t_sumsq += t_predicted*t_predicted;

    total_stats->merge_stats(path_stats);
    total_stats->merge_stats(top_stats);

    if(i % ep->print_every == 0) {

      dive_hist.open(dive_hist_file(),std::ofstream::app);
      print_dive(dive_hist,top_stats,path_stats,ep->time_scale);
      dive_hist.close();
    }
    if((i % ep->print_everything_every == 0)&&(i<ep->trials)) {
      print_time();
      printf("[ %5d dives  %14.0f nodes  %12.2f seconds %d ignored]\n",
		i,double(total_stats->node)/i,
		ep->time_scale*total_stats->total_time/i,too_deep);

      estimated_stats = new QAPStatistics(*total_stats);
      estimated_stats->scale(1.0/i);
      estimated_stats->scaleTimes(ep->time_scale);
      estimated_stats->print();
      delete estimated_stats;
    }
  }

  total_stats->scale(1.0/ep->trials);
  total_stats->scaleTimes(ep->time_scale);
  total_stats->print();

  double t_mean,n_mean;

  t_mean = total_stats->total_time;
  n_mean = total_stats->node;
  standard_errs(n_mean,n_sumsq,t_mean,t_sumsq,ep->trials);

  std::cout << too_deep << " dives with depth >= " << ep->ignore_depth 
       << " were thrown out." << std::endl;

  // free up memory...
  delete path_stats;
  delete top_stats;
  // estimated_stats already freed.

  return total_stats;
}

void standard_errs(double n_mean,double n_sumsq,
		   double t_mean,double t_sumsq,
		   double trials)
{
  double t_var, n_var;
  double t_stddev, n_stddev;

  std::cout << "tmean = " << t_mean <<  "\t nmean = " << n_mean << std::endl;
  std::cout << "tssq  = " << t_sumsq << "\t nssq  = " << n_sumsq << std::endl;

  t_var  = t_sumsq/(trials) - t_mean*t_mean;
  n_var  = n_sumsq/(trials) - n_mean*n_mean;

  t_stddev = sqrt(t_var);
  n_stddev = sqrt(n_var);

  std::cout << "tvar  = " << t_var << "\t nvar  = " << n_var << std::endl;
  std::cout << "time stddev = " << t_stddev << ", err = " << t_stddev/sqrt(trials) << std::endl;
  std::cout << "node stddev = " << n_stddev << ", err = " << n_stddev/sqrt(trials) << std::endl;
}

void readEstimateParameters(const char *filename,EstimateParameters *&ep)
{
  if(!ep) ep = new EstimateParameters();

  std::ifstream f;

  f.open(filename);

  if(!f)
    {
      std::cout << "sorry chief, couldn't find [" << filename << "].\n";
      f.close();
      return;
    }

  if(f)
    f >> ep->start_dive_level;
  if(f)
    f >> ep->trials;
  if(f)
    f >> ep->print_every;
  if(f)
    f >> ep->print_everything_every;
  if(f)
    f >> ep->score_exponent;
  if(f)
    f >> ep->ignore_depth;

  if(!f)
    {
      std::cout << "warning: couldn't read estimator parameters, using defaults" 
	   << std::endl;
    }
  f.close();

  printf("start dives at %d, %d trials, print every %d/%d, score exp = %f\n",
	    ep->start_dive_level,ep->trials,ep->print_every,
	    ep->print_everything_every,ep->score_exponent);
  printf("throw out dives with depth >= %d\n",ep->ignore_depth);
}

/*
  return depth of deepest node in s.
 */
int dive_depth(QAPStatistics *s)
{
  int depth = 0;

  // find first level with positive number of nodes.
  while((depth <= s->n)&&(s->n_level[depth] == 0)) depth++;

  if(depth > s->n) return depth;

  // loop until we find a level with no nodes.
  while((depth <= s->n)&&(s->n_level[depth] != 0)) depth++;
  
  if(depth > s->n) return depth;

  // depth-1 is the last level reached.
  return depth - 1;
}



