// The settings registry: see settings.h.
#include "settings.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cerrno>
#include <climits>
#include <iostream>
#include <sstream>
#include <vector>

extern char **environ;

static const double INF = 1e300;

// The table. Defaults are the solver's behaviour with the variable unset
// (the nug30 record recipe of Nov 8 2026 for the DNN bound). Names that
// are not environment variables are the policy flags (-x, -w) and the
// params file, which driver.C fills in with settings_set.
static Setting TABLE[] = {
  // ---- the search (bnb.C): the DNN bound
  {"QAP_DNN_EXT", 'i', "20", 0, 1000, SET_POLICY,
   "levels below the last DNN depth at which a child of a DNN miss keeps the bound"},
  {"QAP_DNN_GAMMA1", 'd', "1.0", 0, 2, SET_POLICY, "Peaceman-Rachford first dual step (0: none)"},
  {"QAP_DNN_GAMMA", 'd', "1.0", 0, 2, SET_POLICY, "ADMM dual step length"},
  {"QAP_DNN_CHEB", 'i', "3", 1, 12, SET_POLICY, "Chebyshev filter degree of the subspace eigensolver"},
  {"QAP_DNN_TIGHT_LAM", 'i', "2", 0, 10, SET_POLICY, "the filter's upper end lam_hi from this many Lanczos steps (0: off)"},
  {"QAP_DNN_GAP_RATIO", 'd', "0", 0, INF, SET_POLICY, "gap-ratio give-up rule in place of the extrapolation (0: off)"},
  {"QAP_DNN_GU_ZONE", 'd', "0.003", 0, 1, SET_POLICY,
   "a give-up within this fraction of the cutoff is confirmed 100 iterations later (0: off)"},
  {"QAP_DNN_MOMENTUM", 'd', "0", 0, 1, SET_POLICY, "fast ADMM momentum with Goux's safeguards (0: off)"},
  {"QAP_DNN_BETA", 'd', "0.21", 1e-6, 100, SET_POLICY, "ADMM penalty, as a multiple of dnn.py's m/3"},
  {"QAP_DNN_REFRESH", 'i', "120", 1, 1000000, SET_POLICY, "full eigendecomposition every this many iterations"},
  {"QAP_DNN_WARM_Y", 'i', "2", 0, 2, SET_POLICY,
   "a child starts from its parent's Y too: 0 Z only, 1 restricted, 2 conditioned on its assignment"},
  {"QAP_DNN_SHIP", 'i', "32", 0, 64, SET_POLICY,
   "the parent's DNN state in state files: 0 none, 16 bf16, 32 fp32, 64 fp64"},
  // ---- the cut phase and the pool (bnb.C)
  {"QAP_DNN_CUT_ZONE", 'd', "0", 0, 1, SET_POLICY,
   "run the cut phase when the uncut bound ends within this fraction of |incumbent| of the cutoff (0: off)"},
  {"QAP_DNN_CUT_MAXDEPTH", 'i', "1073741824", 0, 1073741824, SET_POLICY, "the cut phase only at nodes of at most this depth"},
  {"QAP_DNN_CUT_SIDE", 'i', "0", 0, 1, SET_POLICY, "regions over the first (0) or the second (1) matrix"},
  {"QAP_DNN_CUT_GRID", 's', "", 0, 0, SET_POLICY, "R:C grid whose prefixes are the regions (unset: prefixes of the free elements)"},
  {"QAP_DNN_CUT_ROUNDS", 'i', "12", 1, 1000, SET_POLICY, "cut rounds per phase"},
  {"QAP_DNN_CUT_ITERS", 'i', "300", 1, 1000000, SET_POLICY, "ADMM iterations per cut round"},
  {"QAP_DNN_CUT_K", 'i', "1000", 0, 1000000, SET_POLICY, "new triangle cuts per round"},
  {"QAP_DNN_CUT_K5", 'i', "0", 0, 1000000, SET_POLICY, "5-point count inequalities per round (0: off)"},
  {"QAP_DNN_CUT_K7", 'i', "0", 0, 1000000, SET_POLICY, "local-7 inequalities per round (0: off)"},
  {"QAP_DNN_CUT_SUB7", 'i', "1000", 1, 10000000, SET_POLICY, "random 7-subsets sampled per region and round"},
  {"QAP_DNN_CUT_KBOX", 's', "", 0, 0, SET_POLICY, "box inequalities per round, K[:QMAX] (unset: none)"},
  {"QAP_DNN_CUT_DELTA", 'd', "50", 0, INF, SET_POLICY, "Polyak target: best + delta"},
  {"QAP_DNN_CUT_STALL", 'i', "3", 0, 1000, SET_POLICY,
   "after this many rounds, stop when the gap exceeds 3x the gain of the last that many (0: off)"},
  {"QAP_DNN_CUT_PGATE", 's', "", 0, 0, SET_POLICY, "primal gate R:M: stop after round R when the best primal is more than M below the cutoff"},
  {"QAP_DNN_CUT_NOCAP", 'i', "1", 0, 1, SET_POLICY, "1: Polyak target uncapped at the cutoff; 0: capped"},
  {"QAP_DNN_CUT_FREEZE", 'i', "600", 0, 1000000, SET_POLICY,
   "iterations of the quiet round with frozen multipliers when a phase ends unfathomed with its primal above the cutoff (0: off)"},
  {"QAP_DNN_CUT_WARMY", 'i', "0", 0, 1, SET_POLICY, "1: cut rounds restart from the previous round's Y as well as its Z"},
  {"QAP_DNN_CUT_DAMP1", 's', "", 0, 0, SET_POLICY, "F[:G]: factor on the first (and later) Polyak steps of a phase that starts from a pool"},
  {"QAP_DNN_CUT_EQUIP", 's', "", 0, 0, SET_POLICY, "G:R: a phase more than G below the cutoff runs at most R rounds, to equip the children"},
  {"QAP_DNN_CUT_CP", 's', "", 0, 0, SET_POLICY, "k:K: the phase as one CP-ADMM solve separating K triangles every k iterations"},
  {"QAP_DNN_CUT_CP_ITERS", 'i', "0", 0, 100000000, SET_POLICY, "CP-ADMM iteration limit (0: rounds x iters + freeze)"},
  {"QAP_DNN_CUT_CP_SWEEPS", 's', "1:1", 0, 0, SET_POLICY, "S[:P]: Dykstra sweeps and passes per CP-ADMM iteration"},
  {"QAP_DNN_POOL", 'i', "0", 0, 1, SET_POLICY, "1: a node hands its cut pool to its children"},
  {"QAP_DNN_POOL_CP", 'i', "0", 0, 1, SET_POLICY, "1: a pooled child's first solve in CP form"},
  {"QAP_DNN_POOL_CHECK", 'i', "0", 0, 1000000, SET_DIAGNOSTIC, "check each assembled pool term on this many random permutations"},
  {"QAP_DNN_TRACE", 's', "", 0, 0, SET_DIAGNOSTIC, "append a line per DNN node to this file"},
  // ---- dnn.C
  {"DNN_AMX", 'i', "0", 0, 2, SET_POLICY, "AMX bf16 kernels: 1 W*Y, 2 also CholQR and Rayleigh-Ritz (AMX hosts only)"},
  {"DNN_R32", 'i', "0", 0, 1, SET_EXPERIMENT, "1: R's rank term in fp32 on certificate iterations too"},
  {"DNN_BETA_SCHED", 's', "w=0.36", 0, 0, SET_POLICY,
   "ADMM penalty schedule, e.g. two,f=1.6667,t=60,w=0.25,wc=0.7 (\"\": off)"},
  {"DNN_ADMM_CROSS", 'i', "2", 0, 100, SET_POLICY, "extra certificates at the predicted cutoff crossing (0: off)"},
  {"DNN_ADMM_CROSS_FIRST", 'i', "0", 0, 1000, SET_EXPERIMENT, "one more extra certificate at this iteration of a warm call (0: off)"},
  {"DNN_FUSE", 'i', "0", 0, 1, SET_EXPERIMENT, "1: fuse the Y/Z step with the next forward's row phase"},
  {"DNN_ADMM_WARM_BETA", 's', "", 0, 0, SET_EXPERIMENT, "F1:T:F2[:force]: a warm call's own penalty schedule"},
  {"DNN_ADMM_WARM_DAMP", 's', "", 0, 0, SET_EXPERIMENT, "D:K: both dual steps scaled by D for the first K iterations of a warm call"},
  {"DNN_CORE_LAZYCERT", 'i', "1", 0, 1, SET_POLICY, "1: Rump's test only at checks that can decide something"},
  {"DNN_CP_LAZYCERT", 'i', "0", 0, 1, SET_POLICY, "1: the lazy certificate in CP-ADMM solves too"},
  // ---- the estimator (estimate.C) and branching diagnostics (branch.C)
  {"QAP_EST_XBAR", 'd', "0", 0, INF, SET_POLICY, "estimator: child weights x_ij^a from the DNN marginals (0: QPB scores)"},
  {"QAP_EST_XBAR_MIX", 'd', "0.1", 0, 1, SET_POLICY, "estimator: uniform share mixed into the x_ij^a weights"},
  {"QAP_EST_SEED", 's', "", 0, 0, SET_RUN, "estimator: seed for rand() (unset: the fixed sequence)"},
  {"QAP_EST_HIST", 's', "dive_history_nug25_2.txt", 0, 0, SET_DIAGNOSTIC, "estimator: the per-dive history file"},
  {"QAP_XBAR_DEBUG", 'i', "0", 0, 1, SET_DIAGNOSTIC, "1: print the rule 8/9 line choice per node to stderr"},
  {"QAP_ALLOW_UNKNOWN", 'i', "0", 0, 1, SET_RUN, "1: an unknown QAP_DNN_/QAP_EST_/QAP_XBAR_/DNN_ variable is a warning, not an error"},
  // ---- the policy flags (driver.C) and the params file
  {"dnn.depths", 's', "", 0, 0, SET_POLICY, "-x dmin:dmax: the DNN bound's depth range (unset: no DNN bound)"},
  {"dnn.maxgap", 'd', "1", 0, INF, SET_POLICY, "-x third field: relative gap at most this for the DNN bound"},
  {"dnn.giveup", 'i', "0", 0, 100000000, SET_POLICY, "-x fourth field: give up on a hopeless node after this many iterations (0: never)"},
  {"dnn.warm", 'i', "0", 0, 1, SET_POLICY, "-x fifth field: 1: a child's ADMM starts from its parent's final state"},
  {"fw.patience", 'i', "0", 0, 100000000, SET_POLICY, "-w: Frank-Wolfe patience at lookahead/strong levels (0: off)"},
  {"par.file", 's', "", 0, 0, SET_RUN, "-p: the params file"},
  {"par.hash", 's', "", 0, 0, SET_POLICY, "FNV-1a of the params file's bytes"},
};
static const int NTABLE = (int)(sizeof TABLE / sizeof TABLE[0]);

static const char *PREFIXES[] = {"QAP_DNN_", "QAP_EST_", "QAP_XBAR_", "DNN_", NULL};
static const char *CLASS_NAMES[] = {"policy", "experiment", "run", "diagnostic"};
static const char *SOURCE_NAMES[] = {"default", "env", "flag", "par"};

static Setting *find(const char *name)
{
  for (int k = 0; k < NTABLE; k++)
    if (strcmp(TABLE[k].name, name) == 0) return &TABLE[k];
  return NULL;
}

static void die(const std::string &msg)
{
  std::cerr << "settings: " << msg << std::endl;
  exit(2);
}

// Parse `value` into `s` or exit 2.
static void assign(Setting &s, const std::string &value, SettingSource src)
{
  if (s.type == 'i') {
    errno = 0;
    char *end = NULL;
    const long v = strtol(value.c_str(), &end, 10);
    if (value.empty() || *end || errno) die(std::string(s.name) + "=" + value + " is not an integer");
    if (v < s.lo || v > s.hi) {
      std::ostringstream o;
      o << s.name << "=" << value << " is outside " << (long)s.lo << ".." << (long)s.hi;
      die(o.str());
    }
    s.i = v; s.d = (double)v;
  } else if (s.type == 'd') {
    errno = 0;
    char *end = NULL;
    const double v = strtod(value.c_str(), &end);
    if (value.empty() || *end || errno || !std::isfinite(v)) die(std::string(s.name) + "=" + value + " is not a number");
    if (v < s.lo || v > s.hi) {
      std::ostringstream o;
      o << s.name << "=" << value << " is outside " << s.lo << ".." << s.hi;
      die(o.str());
    }
    s.d = v; s.i = (long)v;
  } else {
    if (value.find_first_of(" \t\n\r") != std::string::npos)
      die(std::string(s.name) + " must not contain whitespace");
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
    if (!s) { unknown.push_back(name); continue; }
    assign(*s, eq + 1, SRC_ENV);
  }
  if (!unknown.empty()) {
    std::string list;
    for (size_t k = 0; k < unknown.size(); k++) list += (k ? ", " : "") + unknown[k];
    if (lenient)
      std::cerr << "settings: warning: unknown variable(s) ignored (QAP_ALLOW_UNKNOWN=1): " << list << std::endl;
    else
      die("unknown variable(s) " + list + ": not a setting of this solver (bnbdrv --print-settings lists them;"
          " QAP_ALLOW_UNKNOWN=1 to run anyway)");
  }
}

const Setting &setting(const char *name)
{
  settings_load();
  Setting *s = find(name);
  if (!s) {
    std::cerr << "settings: internal error: no entry named " << name << std::endl;
    abort();
  }
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
  if (!s) {
    std::cerr << "settings: internal error: no entry named " << name << std::endl;
    abort();
  }
  assign(*s, value, src);
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
    if (recorded(TABLE[k])) h = fnv1a(std::string(TABLE[k].name) + "=" + TABLE[k].text + "\n", h);
  return hex64(h);
}

void settings_dump(std::ostream &os)
{
  settings_load();
  os << "settings " << NTABLE << " hash " << hash_text() << '\n';
  for (int k = 0; k < NTABLE; k++) {
    const Setting &s = TABLE[k];
    os << "setting " << s.name << ' ' << (s.text.empty() ? "\"\"" : s.text) << ' '
       << SOURCE_NAMES[s.src] << ' ' << CLASS_NAMES[s.cls] << '\n';
  }
  os.flush();
}

std::string settings_record()
{
  settings_load();
  std::string r = hash_text();
  for (int k = 0; k < NTABLE; k++)
    if (recorded(TABLE[k]) && TABLE[k].src != SRC_DEFAULT)
      r += std::string(" ") + TABLE[k].name + "=" + TABLE[k].text;
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
