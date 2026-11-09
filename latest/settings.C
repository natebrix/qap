// The settings registry: see settings.h.
#include "settings.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cerrno>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

extern char **environ;

static const double INF = 1e300;

// The table. Defaults are the solver's behaviour with nothing set (the nug30
// record recipe of Nov 8 2026 for the DNN bound). A policy entry has a key
// and the environment name it had before Phase 2 (an error naming the key
// now); experiment, run and diagnostic entries have an environment name.
#define POLICY(key, old, type, def, lo, hi, desc) {key, NULL, old, type, def, lo, hi, SET_POLICY, desc}
#define ENVVAR(env, type, def, lo, hi, cls, desc) {NULL, env, NULL, type, def, lo, hi, cls, desc}
static Setting TABLE[] = {
  // ---- the DNN bound (bnb.C -x, dnn.C)
  POLICY("dnn.depths", NULL, 's', "", 0, 0, "-x dmin:dmax: the DNN bound's depth range (unset: no DNN bound)"),
  POLICY("dnn.maxgap", NULL, 'd', "1", 0, INF, "-x third field: relative gap (inc - QPB)/(inc - root QPB) at most this"),
  POLICY("dnn.giveup", NULL, 'i', "0", 0, 100000000, "-x fourth field: give up on a hopeless node after this many iterations (0: never)"),
  POLICY("dnn.warm", NULL, 'i', "0", 0, 1, "-x fifth field: 1: a child's ADMM starts from its parent's final state"),
  POLICY("dnn.ext", "QAP_DNN_EXT", 'i', "20", 0, 1000, "levels below the last DNN depth at which a child of a DNN miss keeps the bound"),
  POLICY("dnn.gamma1", "QAP_DNN_GAMMA1", 'd', "1.0", 0, 2, "Peaceman-Rachford first dual step (0: none)"),
  POLICY("dnn.gamma", "QAP_DNN_GAMMA", 'd', "1.0", 0, 2, "ADMM dual step length"),
  POLICY("dnn.cheb", "QAP_DNN_CHEB", 'i', "3", 1, 12, "Chebyshev filter degree of the subspace eigensolver"),
  POLICY("dnn.tight_lam", "QAP_DNN_TIGHT_LAM", 'i', "2", 0, 10, "the filter's upper end lam_hi from this many Lanczos steps (0: off)"),
  POLICY("dnn.gap_ratio", "QAP_DNN_GAP_RATIO", 'd', "0", 0, INF, "gap-ratio give-up rule in place of the extrapolation (0: off)"),
  POLICY("dnn.gu_zone", "QAP_DNN_GU_ZONE", 'd', "0.003", 0, 1,
         "a give-up within this fraction of the cutoff is confirmed 100 iterations later (0: off)"),
  POLICY("dnn.momentum", "QAP_DNN_MOMENTUM", 'd', "0", 0, 10, "fast ADMM with Goux's safeguards, restarted when the residual grows by this factor (1.05; 0: off)"),
  POLICY("dnn.beta", "QAP_DNN_BETA", 'd', "0.21", 1e-6, 100, "ADMM penalty, as a multiple of dnn.py's m/3"),
  POLICY("dnn.refresh", "QAP_DNN_REFRESH", 'i', "120", 1, 1000000, "full eigendecomposition every this many iterations"),
  POLICY("dnn.warm_y", "QAP_DNN_WARM_Y", 'i', "2", 0, 2,
         "a child starts from its parent's Y too: 0 Z only, 1 restricted, 2 conditioned on its assignment"),
  POLICY("dnn.beta_sched", "DNN_BETA_SCHED", 's', "w=0.36", 0, 0,
         "ADMM penalty schedule, e.g. two,f=1.6667,t=60,w=0.25,wc=0.7 (\"\": off)"),
  POLICY("dnn.cross", "DNN_ADMM_CROSS", 'i', "2", 0, 100, "extra certificates at the predicted cutoff crossing (0: off)"),
  POLICY("dnn.lazycert", "DNN_CORE_LAZYCERT", 'i', "1", 0, 1, "1: Rump's test only at checks that can decide something"),
  POLICY("dnn.amx", "DNN_AMX", 'i', "0", 0, 2, "AMX bf16 kernels: 1 W*Y, 2 also CholQR and Rayleigh-Ritz (no effect without AMX)"),
  POLICY("ship.bits", "QAP_DNN_SHIP", 'i', "32", 0, 64, "the parent's DNN state in state files: 0 none, 16 bf16, 32 fp32, 64 fp64"),
  // ---- the cut phase and the pool (bnb.C, dnncut.C)
  POLICY("cut.zone", "QAP_DNN_CUT_ZONE", 'd', "0", 0, 1,
         "run the cut phase when the uncut bound ends within this fraction of |incumbent| of the cutoff (0: off)"),
  POLICY("cut.maxdepth", "QAP_DNN_CUT_MAXDEPTH", 'i', "1073741824", 0, 1073741824, "the cut phase only at nodes of at most this depth"),
  POLICY("cut.side", "QAP_DNN_CUT_SIDE", 'i', "0", 0, 1, "regions over the first (0) or the second (1) matrix"),
  POLICY("cut.grid", "QAP_DNN_CUT_GRID", 's', "", 0, 0, "R:C grid whose prefixes are the regions (unset: prefixes of the free elements)"),
  POLICY("cut.rounds", "QAP_DNN_CUT_ROUNDS", 'i', "12", 1, 1000, "cut rounds per phase"),
  POLICY("cut.iters", "QAP_DNN_CUT_ITERS", 'i', "300", 1, 1000000, "ADMM iterations per cut round"),
  POLICY("cut.k", "QAP_DNN_CUT_K", 'i', "1000", 0, 1000000, "new triangle cuts per round"),
  POLICY("cut.k5", "QAP_DNN_CUT_K5", 'i', "0", 0, 1000000, "5-point count inequalities per round (0: off; rounds form only)"),
  POLICY("cut.k7", "QAP_DNN_CUT_K7", 'i', "0", 0, 1000000, "local-7 inequalities per round (0: off; rounds form only)"),
  POLICY("cut.sub7", "QAP_DNN_CUT_SUB7", 'i', "1000", 1, 10000000, "random 7-subsets sampled per region and round (rounds form only)"),
  POLICY("cut.kbox", "QAP_DNN_CUT_KBOX", 's', "", 0, 0, "box inequalities per round, K[:QMAX] (unset: none; rounds form only)"),
  POLICY("cut.delta", "QAP_DNN_CUT_DELTA", 'd', "50", 0, INF, "Polyak target: best + delta (rounds form only)"),
  POLICY("cut.stall", "QAP_DNN_CUT_STALL", 'i', "3", 0, 1000,
         "after this many rounds, stop when the gap exceeds 3x the gain of the last that many (0: off)"),
  POLICY("cut.pgate", "QAP_DNN_CUT_PGATE", 's', "", 0, 0, "primal gate R:M: stop after round R when the best primal is more than M below the cutoff (rounds form only)"),
  POLICY("cut.nocap", "QAP_DNN_CUT_NOCAP", 'i', "1", 0, 1, "1: Polyak target uncapped at the cutoff; 0: capped (rounds form only)"),
  POLICY("cut.freeze", "QAP_DNN_CUT_FREEZE", 'i', "600", 0, 1000000,
         "iterations of the quiet round with frozen multipliers when a phase ends unfathomed with its primal above the cutoff (0: off)"),
  POLICY("cut.warmy", "QAP_DNN_CUT_WARMY", 'i', "0", 0, 1, "1: cut rounds restart from the previous round's Y as well as its Z"),
  POLICY("cut.damp1", "QAP_DNN_CUT_DAMP1", 's', "", 0, 0, "F[:G]: factor on the first (and later) Polyak steps of a phase that starts from a pool (rounds form only)"),
  POLICY("cut.equip", "QAP_DNN_CUT_EQUIP", 's', "", 0, 0, "G:R: a phase more than G below the cutoff runs at most R rounds, to equip the children"),
  POLICY("cut.cp", "QAP_DNN_CUT_CP", 's', "", 0, 0, "k:K: the phase as one CP-ADMM solve separating K triangles every k iterations"),
  POLICY("cut.cp_iters", "QAP_DNN_CUT_CP_ITERS", 'i', "0", 0, 100000000, "CP-ADMM iteration limit (0: rounds x iters + freeze)"),
  POLICY("cut.cp_sweeps", "QAP_DNN_CUT_CP_SWEEPS", 's', "1:1", 0, 0, "S[:P]: Dykstra sweeps and passes per CP-ADMM iteration"),
  POLICY("cut.cp_lazycert", "DNN_CP_LAZYCERT", 'i', "0", 0, 1, "1: the lazy certificate in CP-ADMM solves too"),
  POLICY("pool.on", "QAP_DNN_POOL", 'i', "0", 0, 1, "1: a node hands its cut pool to its children"),
  POLICY("pool.cp", "QAP_DNN_POOL_CP", 'i', "0", 0, 1, "1: a pooled child's first solve in CP form"),
  // ---- branching (branch.C), Frank-Wolfe (fw.C), the estimator (estimate.C)
  POLICY("branch.weighted", NULL, 'i', "1", 0, 1, "rules 1, 3 and 4: weight the U sums by the row's gap fraction"),
  POLICY("branch.weight", NULL, 'd', "0.5", 0, 1, "the weight of that fraction"),
  POLICY("branch.xbar_lines", NULL, 'i', "1", 1, 64, "rule 9: candidate rows and columns per side from the DNN marginals"),
  POLICY("fw.patience", NULL, 'i', "0", 0, 100000000, "-w: Frank-Wolfe patience at lookahead/strong levels (0: off)"),
  POLICY("est.xbar", "QAP_EST_XBAR", 'd', "0", 0, INF, "estimator: child weights x_ij^a from the DNN marginals (0: QPB scores)"),
  POLICY("est.xbar_mix", "QAP_EST_XBAR_MIX", 'd', "0.1", 0, 1, "estimator: uniform share mixed into the x_ij^a weights"),
  POLICY("par.hash", NULL, 's', "", 0, 0, "FNV-1a of the params file's bytes"),
  // ---- experiments in flight (environment)
  ENVVAR("DNN_R32", 'i', "0", 0, 1, SET_EXPERIMENT, "1: R's rank term in fp32 on certificate iterations too"),
  ENVVAR("DNN_ADMM_CROSS_FIRST", 'i', "0", 0, 1000, SET_EXPERIMENT, "one more extra certificate at this iteration of a warm call (0: off)"),
  ENVVAR("DNN_FUSE", 'i', "0", 0, 1, SET_EXPERIMENT, "1: fuse the Y/Z step with the next forward's row phase"),
  ENVVAR("DNN_ADMM_WARM_BETA", 's', "", 0, 0, SET_EXPERIMENT, "F1:T:F2[:force]: a warm call's own penalty schedule"),
  ENVVAR("DNN_ADMM_WARM_DAMP", 's', "", 0, 0, SET_EXPERIMENT, "D:K: both dual steps scaled by D for the first K iterations of a warm call"),
  // ---- this run (environment, -p)
  ENVVAR("QAP_EST_SEED", 's', "", 0, 0, SET_RUN, "estimator: seed for rand() (unset: the fixed sequence)"),
  ENVVAR("QAP_DNN_WARM_DIR", 's', "", 0, 0, SET_RUN, "parent-state store: write each parent's DNN state to this directory and reference it from state files (unset: inline)"),
  ENVVAR("QAP_ALLOW_UNKNOWN", 'i', "0", 0, 1, SET_RUN, "1: an unknown or pre-Phase-2 QAP_DNN_/QAP_EST_/QAP_XBAR_/DNN_ variable is a warning, not an error"),
  {"par.file", NULL, NULL, 's', "", 0, 0, SET_RUN, "-p: the params file"},
  // ---- diagnostics (environment)
  ENVVAR("QAP_DNN_TRACE", 's', "", 0, 0, SET_DIAGNOSTIC, "append a line per DNN node to this file"),
  ENVVAR("QAP_DNN_POOL_CHECK", 'i', "0", 0, 1000000, SET_DIAGNOSTIC, "check each assembled pool term on this many random permutations"),
  ENVVAR("QAP_EST_HIST", 's', "dive_history_nug25_2.txt", 0, 0, SET_DIAGNOSTIC, "estimator: the per-dive history file"),
  ENVVAR("QAP_XBAR_DEBUG", 'i', "0", 0, 1, SET_DIAGNOSTIC, "1: print the rule 8/9 line choice per node to stderr"),
};
#undef POLICY
#undef ENVVAR
static const int NTABLE = (int)(sizeof TABLE / sizeof TABLE[0]);

static const char *PREFIXES[] = {"QAP_DNN_", "QAP_EST_", "QAP_XBAR_", "DNN_", NULL};
static const char *CLASS_NAMES[] = {"policy", "experiment", "run", "diagnostic"};
static const char *SOURCE_NAMES[] = {"default", "env", "flag", "par"};

// Names that used to be settings and were retired outright (Nov 9 2026,
// docs/settings_plan.md step 1: the winning value is the only behaviour, or
// the experiment was removed). Errors like any unknown name, with a message.
struct Retired { const char *name; const char *what; };
static const Retired RETIRED[] = {
  {"DNN_CHOLQR_PASSES", "one CholQR pass is the only behaviour"},
  {"DNN_R32_NONCERT", "R's rank term in fp32 off certificate iterations is the only behaviour (DNN_R32 remains)"},
  {"DNN_CORE_NOSYMV", "no symmetrisation pass is the only behaviour"},
  {"DNN_NOSYMW", "no symmetrisation of W is the only behaviour"},
  {"DNN_CORE_REFRESH_PARTIAL", "the partial refresh is the only behaviour"},
  {"DNN_GUARD_DIV", "the subspace guard is K/20"},
  {"DNN_CERTX", "the certificate experiment was removed (its report: experiments/wil50_feasibility/roadmap/certificate)"},
};

static Setting *find(const char *name)
{
  for (int k = 0; k < NTABLE; k++)
    if ((TABLE[k].key && strcmp(TABLE[k].key, name) == 0) || (TABLE[k].env && strcmp(TABLE[k].env, name) == 0))
      return &TABLE[k];
  return NULL;
}

static Setting *find_old(const std::string &name)
{
  for (int k = 0; k < NTABLE; k++)
    if (TABLE[k].old_env && name == TABLE[k].old_env) return &TABLE[k];
  return NULL;
}

static const char *retired_note(const std::string &name)
{
  for (size_t k = 0; k < sizeof RETIRED / sizeof RETIRED[0]; k++)
    if (name == RETIRED[k].name || (std::string(RETIRED[k].name) == "DNN_CERTX" && name.compare(0, 9, "DNN_CERTX") == 0))
      return RETIRED[k].what;
  return NULL;
}

static void die(const std::string &msg)
{
  std::cerr << "settings: " << msg << std::endl;
  exit(2);
}

static void internal(const char *name)
{
  std::cerr << "settings: internal error: no entry named " << name << std::endl;
  abort();
}

// Parse `value` into `s` or exit 2.
static void assign(Setting &s, const std::string &value, SettingSource src)
{
  if (s.type == 'i') {
    errno = 0;
    char *end = NULL;
    const long v = strtol(value.c_str(), &end, 10);
    if (value.empty() || *end || errno) die(std::string(s.name()) + "=" + value + " is not an integer");
    if (v < s.lo || v > s.hi) {
      std::ostringstream o;
      o << s.name() << "=" << value << " is outside " << (long)s.lo << ".." << (long)s.hi;
      die(o.str());
    }
    s.i = v; s.d = (double)v;
  } else if (s.type == 'd') {
    errno = 0;
    char *end = NULL;
    const double v = strtod(value.c_str(), &end);
    if (value.empty() || *end || errno || !std::isfinite(v)) die(std::string(s.name()) + "=" + value + " is not a number");
    if (v < s.lo || v > s.hi) {
      std::ostringstream o;
      o << s.name() << "=" << value << " is outside " << s.lo << ".." << s.hi;
      die(o.str());
    }
    s.d = v; s.i = (long)v;
  } else {
    if (value.find_first_of(" \t\n\r") != std::string::npos)
      die(std::string(s.name()) + " must not contain whitespace");
    s.i = 0; s.d = 0;
  }
  s.text = value;
  s.src = src;
}

static bool loaded = false;

void settings_load()
{
  if (loaded) return;
  loaded = true;
  for (int k = 0; k < NTABLE; k++) assign(TABLE[k], TABLE[k].def, SRC_DEFAULT);
  const char *allow = getenv("QAP_ALLOW_UNKNOWN");
  const bool lenient = allow && atoi(allow) > 0;
  std::vector<std::string> unknown;
  for (char **e = environ; e && *e; e++) {
    const char *eq = strchr(*e, '=');
    if (!eq) continue;
    const std::string name(*e, eq - *e);
    bool ours = false;
    for (const char **p = PREFIXES; *p; p++)
      if (name.compare(0, strlen(*p), *p) == 0) ours = true;
    if (!ours && name != "QAP_ALLOW_UNKNOWN") continue;
    Setting *s = find(name.c_str());
    if (s && s->env) { assign(*s, eq + 1, SRC_ENV); continue; }
    if (Setting *o = find_old(name)) {
      unknown.push_back(name + " (now the params-file key " + o->key + ": write `" + o->key + " = " + (eq + 1) +
                        "` in a [settings] section; bench/envpar.py converts a recipe)");
      continue;
    }
    const char *r = retired_note(name);
    unknown.push_back(r ? name + " (retired: " + r + "; unset it)" : name);
  }
  if (!unknown.empty()) {
    std::string list;
    for (size_t k = 0; k < unknown.size(); k++) list += (k ? "; " : "") + unknown[k];
    if (lenient)
      std::cerr << "settings: warning: variable(s) ignored (QAP_ALLOW_UNKNOWN=1): " << list << std::endl;
    else
      die("variable(s) not accepted: " + list + " (bnbdrv --print-settings lists the settings;"
          " QAP_ALLOW_UNKNOWN=1 to run anyway)");
  }
}

static std::string trim(const std::string &s)
{
  const size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

void settings_read_par(const char *path)
{
  settings_load();
  std::ifstream f(path);
  if (!f) return;   // qap_read_params reports a missing file
  std::string line;
  bool in = false;
  int lineno = 0;
  std::vector<std::string> bad;
  while (std::getline(f, line)) {
    lineno++;
    const std::string t = trim(line);
    if (!in) {
      if (t == "[settings]") in = true;
      else if (t.size() > 1 && t[0] == '[' && t.compare(0, 6, "[block") == 0)
        die(std::string(path) + ":" + std::to_string(lineno) + ": per-block settings ([block k]) are not supported yet");
      continue;
    }
    if (t.empty() || t[0] == '#') continue;
    if (t[0] == '[') {
      if (t.compare(0, 6, "[block") == 0)
        die(std::string(path) + ":" + std::to_string(lineno) + ": per-block settings ([block k]) are not supported yet");
      die(std::string(path) + ":" + std::to_string(lineno) + ": unknown section " + t);
    }
    const size_t eq = t.find('=');
    if (eq == std::string::npos) die(std::string(path) + ":" + std::to_string(lineno) + ": expected `key = value`, got `" + t + "`");
    std::string key = trim(t.substr(0, eq)), value = trim(t.substr(eq + 1));
    const size_t hash = value.find(" #");   // a trailing comment
    if (hash != std::string::npos) value = trim(value.substr(0, hash));
    Setting *s = find(key.c_str());
    if (!s || !s->key) { bad.push_back(key); continue; }
    if (s->src == SRC_FLAG) continue;   // the flag wins
    assign(*s, value, SRC_PAR);
  }
  if (!bad.empty()) {
    std::string list;
    for (size_t k = 0; k < bad.size(); k++) list += (k ? ", " : "") + bad[k];
    die(std::string(path) + ": unknown [settings] key(s) " + list + " (bnbdrv --print-settings lists the keys)");
  }
}

const Setting &setting(const char *name)
{
  settings_load();
  Setting *s = find(name);
  if (!s) internal(name);
  return *s;
}

long setting_i(const char *name) { return setting(name).i; }
double setting_d(const char *name) { return setting(name).d; }
const char *setting_s(const char *name) { return setting(name).text.c_str(); }
bool setting_given(const char *name) { return setting(name).src != SRC_DEFAULT; }

void settings_set(const char *name, const std::string &value, SettingSource src)
{
  settings_load();
  Setting *s = find(name);
  if (!s) internal(name);
  if (src == SRC_PAR && s->src == SRC_FLAG) return;
  assign(*s, value, src);
}

void settings_override_default(const char *name, const std::string &value)
{
  settings_load();
  Setting *s = find(name);
  if (!s) internal(name);
  if (s->src != SRC_DEFAULT) die(std::string(name) + ": cannot change the default of a setting that is set");
  assign(*s, value, SRC_DEFAULT);
}

static unsigned long long fnv1a(const std::string &bytes, unsigned long long h = 1469598103934665603ULL)
{
  for (size_t k = 0; k < bytes.size(); k++) {
    h ^= (unsigned char)bytes[k];
    h *= 1099511628211ULL;
  }
  return h;
}

static std::string hex64(unsigned long long h)
{
  char b[20];
  snprintf(b, sizeof b, "%016llx", h);
  return b;
}

static bool recorded(const Setting &s) { return s.cls == SET_POLICY || s.cls == SET_EXPERIMENT; }

static std::string hash_text()
{
  settings_load();
  unsigned long long h = 1469598103934665603ULL;
  for (int k = 0; k < NTABLE; k++)
    if (recorded(TABLE[k])) h = fnv1a(std::string(TABLE[k].name()) + "=" + TABLE[k].text + "\n", h);
  return hex64(h);
}

void settings_dump(std::ostream &os)
{
  settings_load();
  os << "settings " << NTABLE << " hash " << hash_text() << '\n';
  for (int k = 0; k < NTABLE; k++) {
    const Setting &s = TABLE[k];
    os << "setting " << s.name() << ' ' << (s.text.empty() ? "\"\"" : s.text) << ' '
       << SOURCE_NAMES[s.src] << ' ' << CLASS_NAMES[s.cls] << '\n';
  }
  os.flush();
}

// The two listings read the table only (not the environment), so they work
// under any environment, an old recipe's included (bench/envpar.py runs them).
void settings_table(std::ostream &os)
{
  os << "| setting | where | type | default | range | class | what |\n|---|---|---|---|---|---|---|\n";
  for (int k = 0; k < NTABLE; k++) {
    const Setting &s = TABLE[k];
    std::string where = s.key ? (strcmp(s.key, "par.file") == 0 ? "`-p`" : strcmp(s.key, "par.hash") == 0 ? "the `-p` file" : "par key")
                              : "environment";
    if (s.key && strncmp(s.key, "dnn.", 4) == 0 && (strcmp(s.key, "dnn.depths") == 0 || strcmp(s.key, "dnn.maxgap") == 0 ||
                                                      strcmp(s.key, "dnn.giveup") == 0 || strcmp(s.key, "dnn.warm") == 0))
      where = "par key or `-x`";
    if (s.key && strcmp(s.key, "fw.patience") == 0) where = "par key or `-w`";
    std::string range;
    if (s.type == 'i') { std::ostringstream o; o << (long)s.lo << ".." << (s.hi >= 1e9 ? std::string("inf") : std::to_string((long)s.hi)); range = o.str(); }
    else if (s.type == 'd') { std::ostringstream o; o << s.lo << ".."; if (s.hi >= 1e299) o << "inf"; else o << s.hi; range = o.str(); }
    else range = "text";
    os << "| `" << s.name() << "` | " << where << " | " << (s.type == 'i' ? "integer" : s.type == 'd' ? "number" : "text")
       << " | " << (s.def[0] ? std::string("`") + s.def + "`" : "unset") << " | " << range << " | " << CLASS_NAMES[s.cls]
       << " | " << s.desc << " |\n";
  }
  os.flush();
}

void settings_env_map(std::ostream &os)
{
  for (int k = 0; k < NTABLE; k++)
    if (TABLE[k].old_env) os << TABLE[k].old_env << ' ' << TABLE[k].key << '\n';
  os.flush();
}

std::string settings_record()
{
  settings_load();
  std::string r = hash_text();
  for (int k = 0; k < NTABLE; k++)
    if (recorded(TABLE[k]) && TABLE[k].src != SRC_DEFAULT)
      r += std::string(" ") + TABLE[k].name() + "=" + TABLE[k].text;
  return r;
}

bool settings_check_record(const std::string &record, std::string &why)
{
  const std::string mine = settings_record();
  const std::string h_file = record.substr(0, record.find(' '));
  const std::string h_mine = mine.substr(0, mine.find(' '));
  if (h_file == h_mine) return true;
  why = "the state file was written with settings [" + record + "], this run has [" + mine +
        "] (the hash covers every policy and experiment setting; the list shows those set away from"
        " their defaults, so a changed default shows as a hash change alone)";
  return false;
}

std::string settings_hash_file(const char *path)
{
  FILE *f = fopen(path, "rb");
  if (!f) return "";
  std::string bytes;
  char buf[65536];
  size_t n;
  while ((n = fread(buf, 1, sizeof buf, f)) > 0) bytes.append(buf, n);
  fclose(f);
  return hex64(fnv1a(bytes));
}
