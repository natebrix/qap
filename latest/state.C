#include "state.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

static const char *STATE_MAGIC = "QAP-STATE";
static const int STATE_VERSION = 9;   // 1: no warm starts, 2: no DNN counts, 3: no dnn_ext on the warm starts, 4: no cut pools, 5: no DNN state, 6: no per-level DNN counts, 7: no settings record, 8: no @path DNN references (all still read)
bool state_allow_settings_change = false;   // driver.C --allow-settings-change: a settings mismatch on -r is a warning

// version 7: by level, the DNN bound, its cut phase and the nodes discarded
// when popped (stats.h), in this order after dnn_time
static const char *const DNN_LEVEL_CTRS[] = {"pop_f_level", "dnn_n_level", "dnn_f_level", "dnn_g_level", "dnn_i_level"};
static const char *const CUT_LEVEL_CTRS[] = {"cut_n_level", "cut_f_level"};

#ifdef HAVE_DNN
#include "dnn.h"
#include <cstdint>
#include <unistd.h>
#include <cstdlib>
#include "settings.h"

/*
 Version 6: a warm start may carry the parent's final DNN state (DNNWarm:
 Z and Y, symmetric, (m^2+1)^2 each), so that its children start their ADMM
 from it when they are bounded in another process (a resumed run, a
 coordinator unit; roadmap/convergence2 F2). Written only with
 QAP_DNN_SHIP=16 (bf16), 32 (fp32) or 64 (the exact doubles, for tests:
 a resumed run is then the same run); otherwise, and when the parent has
 none, the token is "-". Token:
   dnn<bits>:<m>:<scale>:<hasY>:<base64 of Z then Y: the upper triangles at 16 and 32 bits, the whole matrices at 64>
 It only moves where a child's ADMM starts: every bound stays certified.
*/
static int dnn_ship_bits()
{
  static const int b = (int)setting_i("ship.bits");
  if (b != 0 && b != 16 && b != 32 && b != 64) {
    std::cerr << "settings: QAP_DNN_SHIP=" << b << " is not 0, 16, 32 or 64" << std::endl;
    exit(2);
  }
  return b;
}

static void free_dnn_state(void *p) { delete (DNNWarm *)p; }

static const char B64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static std::string b64_encode(const std::vector<unsigned char> &in)
{
  std::string out;
  out.reserve((in.size() + 2) / 3 * 4);
  size_t i = 0;
  for (; i + 2 < in.size(); i += 3) {
    const uint32_t v = (uint32_t)in[i] << 16 | (uint32_t)in[i + 1] << 8 | in[i + 2];
    out += B64[v >> 18]; out += B64[(v >> 12) & 63]; out += B64[(v >> 6) & 63]; out += B64[v & 63];
  }
  if (i < in.size()) {
    uint32_t v = (uint32_t)in[i] << 16 | (i + 1 < in.size() ? (uint32_t)in[i + 1] << 8 : 0);
    out += B64[v >> 18]; out += B64[(v >> 12) & 63];
    out += i + 1 < in.size() ? B64[(v >> 6) & 63] : '=';
    out += '=';
  }
  return out;
}

static bool b64_decode(const char *s, size_t len, std::vector<unsigned char> &out)
{
  static int8_t dec[256];
  static bool init = false;
  if (!init) {
    for (int k = 0; k < 256; k++) dec[k] = -1;
    for (int k = 0; k < 64; k++) dec[(unsigned char)B64[k]] = (int8_t)k;
    init = true;
  }
  if (len % 4) return false;
  out.clear();
  out.reserve(len / 4 * 3);
  for (size_t i = 0; i < len; i += 4) {
    int c[4];
    for (int k = 0; k < 4; k++) c[k] = s[i + k] == '=' ? 0 : dec[(unsigned char)s[i + k]];
    if (c[0] < 0 || c[1] < 0 || c[2] < 0 || c[3] < 0) return false;
    const uint32_t v = (uint32_t)c[0] << 18 | (uint32_t)c[1] << 12 | (uint32_t)c[2] << 6 | (uint32_t)c[3];
    out.push_back((unsigned char)(v >> 16));
    if (s[i + 2] != '=') out.push_back((unsigned char)(v >> 8));
    if (s[i + 3] != '=') out.push_back((unsigned char)v);
  }
  return true;
}

static void put_value(std::vector<unsigned char> &b, double x, int bits)
{
  unsigned char tmp[8];
  if (bits == 64) {
    memcpy(tmp, &x, 8);
    b.insert(b.end(), tmp, tmp + 8);
    return;
  }
  const float f = (float)x;
  uint32_t u;
  memcpy(&u, &f, 4);
  if (bits == 32) {
    memcpy(tmp, &u, 4);
    b.insert(b.end(), tmp, tmp + 4);
    return;
  }
  const uint16_t h = (uint16_t)((u + 0x7FFFu + ((u >> 16) & 1u)) >> 16);   // bf16, round to nearest even
  memcpy(tmp, &h, 2);
  b.insert(b.end(), tmp, tmp + 2);
}

static double get_value(const unsigned char *p, int bits)
{
  if (bits == 64) { double x; memcpy(&x, p, 8); return x; }
  uint32_t u;
  if (bits == 32) memcpy(&u, p, 4);
  else { uint16_t h; memcpy(&h, p, 2); u = (uint32_t)h << 16; }
  float f;
  memcpy(&f, &u, 4);
  return f;
}

static std::string dnn_token(const WarmStart *w)
{
  const int bits = dnn_ship_bits();
  if (!bits || !w->dnn) return "-";
  const DNNWarm &d = *(const DNNWarm *)w->dnn;
  const size_t n1 = (size_t)d.m * d.m + 1;
  if (d.m < 1 || d.Z.size() != n1 * n1 || !(d.scale > 0)) return "-";
  const bool hasY = d.Y.size() == n1 * n1;
  // 64 bits ships the whole matrix: since DNN_CORE_NOSYMV (the default since
  // Nov 8 2026) Z and Y are symmetric to rounding only, and the child's warm
  // start reads both triangles, so an upper triangle mirrored on read is not
  // the parent's state (the ship64 test caught it). 16 and 32 bits are lossy
  // anyway and keep the triangle. The reader tells the layouts apart by length.
  const bool full = bits == 64;
  std::vector<unsigned char> b;
  b.reserve((full ? n1 * n1 : n1 * (n1 + 1) / 2) * (hasY ? 2 : 1) * (bits / 8));
  for (int pass = 0; pass < (hasY ? 2 : 1); pass++) {
    const std::vector<double> &M = pass == 0 ? d.Z : d.Y;
    for (size_t i = 0; i < n1; i++)
      for (size_t j = full ? 0 : i; j < n1; j++) put_value(b, M[i * n1 + j], bits);
  }
  char head[96];
  snprintf(head, sizeof head, "dnn%d:%d:%.17g:%d:", bits, d.m, d.scale, hasY ? 1 : 0);
  return head + b64_encode(b);
}

// ---- the parent-state store (QAP_DNN_WARM_DIR; docs/parent_state_store.md)
// With the variable set, a parent's DNN token goes to <dir>/<fnv64>.dnn once
// and the warm line carries "@<path>" (relative to the state file's directory
// when the store is under it, absolute otherwise); a reader takes the token
// from that file. Two writers of the same parent write the same bytes to the
// same name (a temporary file renamed into place), so the second is a no-op;
// a hash collision only costs a wrong warm start's worth of iterations, never
// a wrong bound (every bound is certified from wherever the ADMM starts).

static std::string dir_of(const std::string &path)
{
  size_t k = path.find_last_of('/');
  return k == std::string::npos ? std::string(".") : (k == 0 ? std::string("/") : path.substr(0, k));
}

static std::string absolute(const std::string &path)
{
  if (!path.empty() && path[0] == '/') return path;
  char cwd[4096];
  if (!getcwd(cwd, sizeof cwd)) return path;
  return std::string(cwd) + "/" + path;
}

static std::string fnv64_hex(const std::string &s)
{
  uint64_t h = 1469598103934665603ULL;
  for (size_t i = 0; i < s.size(); i++) { h ^= (unsigned char)s[i]; h *= 1099511628211ULL; }
  char b[24];
  snprintf(b, sizeof b, "%016llx", (unsigned long long)h);
  return b;
}

// the token field for the warm line of a state file at `state_path`: a
// reference into the store, or the inline token when there is no store
static std::string dnn_field(const WarmStart *w, const std::string &state_path)
{
  const char *dir = setting_s("QAP_DNN_WARM_DIR");
  if (!*dir || !w->dnn) { w->dnn_ref.clear(); return dnn_token(w); }
  if (w->dnn_ref.empty()) {
    const std::string tok = dnn_token(w);
    if (tok == "-") return tok;
    const std::string store = absolute(dir);
    const std::string path = store + "/" + fnv64_hex(tok) + ".dnn";
    std::ifstream probe(path.c_str());
    if (!probe) {
      char tmp[4600];
      snprintf(tmp, sizeof tmp, "%s.%d.tmp", path.c_str(), (int)getpid());
      std::ofstream o(tmp);
      o << tok << '\n';
      o.close();
      if (!o || rename(tmp, path.c_str()) != 0) {
        std::cerr << "write_state: cannot write " << path << " (QAP_DNN_WARM_DIR); the token goes inline" << std::endl;
        remove(tmp);
        return tok;
      }
    }
    w->dnn_ref = path;
  }
  // relative to the state file's directory when the store is under it
  const std::string sdir = absolute(dir_of(state_path));
  const std::string base = sdir == "/" ? sdir : sdir + "/";
  if (w->dnn_ref.compare(0, base.size(), base) == 0) return "@" + w->dnn_ref.substr(base.size());
  return "@" + w->dnn_ref;
}

// the token behind a warm line's field, reading the store when it is a
// reference; sets `ref` to the file's absolute path then, "" otherwise
static bool dnn_field_token(const std::string &field, const std::string &state_path, std::string &tok, std::string &ref)
{
  ref.clear();
  if (field.empty() || field[0] != '@') { tok = field; return true; }
  std::string path = field.substr(1);
  if (path.empty()) return false;
  if (path[0] != '/') path = absolute(dir_of(state_path)) + "/" + path;
  std::ifstream f(path.c_str());
  if (!f) { std::cerr << "read_state: cannot open the DNN state file " << path << std::endl; return false; }
  std::getline(f, tok);
  if (tok.empty()) return false;
  ref = path;
  return true;
}

static DNNWarm *dnn_from_token(const std::string &t)
{
  int bits, m, hasY, used = 0;
  double scale;
  if (sscanf(t.c_str(), "dnn%d:%d:%lg:%d:%n", &bits, &m, &scale, &hasY, &used) != 4 || used <= 0 ||
      (bits != 16 && bits != 32 && bits != 64) || m < 1 || m > 1000 || !(scale > 0))
    return NULL;
  std::vector<unsigned char> b;
  if (!b64_decode(t.c_str() + used, t.size() - used, b)) return NULL;
  const size_t n1 = (size_t)m * m + 1, tri = n1 * (n1 + 1) / 2, bytes = (size_t)bits / 8;
  const size_t nm = hasY ? 2 : 1;
  const bool full = b.size() == n1 * n1 * nm * bytes;   // the whole matrix (64 bits) or the upper triangle
  if (!full && b.size() != tri * nm * bytes) return NULL;
  DNNWarm *d = new DNNWarm;
  d->m = m;
  d->scale = scale;
  const unsigned char *p = b.data();
  for (int pass = 0; pass < (int)nm; pass++) {
    std::vector<double> &M = pass == 0 ? d->Z : d->Y;
    M.assign(n1 * n1, 0.0);
    if (full) {
      for (size_t t = 0; t < n1 * n1; t++, p += bytes) M[t] = get_value(p, bits);
    } else {
      for (size_t i = 0; i < n1; i++)
        for (size_t j = i; j < n1; j++, p += bytes) M[i * n1 + j] = M[j * n1 + i] = get_value(p, bits);
    }
  }
  return d;
}
#endif

template <class T>
static void write_array(std::ostream &f, const char *tag, const T *a, int len)
{
  f << tag;
  for (int i = 0; i < len; i++) f << ' ' << a[i];
  f << '\n';
}

template <class T>
static bool read_array(std::istream &f, const char *tag, T *a, int len)
{
  std::string t;
  if (!(f >> t) || t != tag) return false;
  for (int i = 0; i < len; i++)
    if (!(f >> a[i])) return false;
  return true;
}

bool write_state(const char *file, Node *stack, QAPSolution *s,
		 QAPStatistics *stats)
{
  std::string tmp = std::string(file) + ".tmp";
  std::ofstream f(tmp.c_str());
  if (!f) {
    std::cerr << "write_state: cannot write [" << tmp << "]" << std::endl;
    return false;
  }
  f << std::setprecision(17);

  const int n = s->n, L = stats->n + 1, m = stats->m;
  f << STATE_MAGIC << ' ' << STATE_VERSION << '\n';
  f << "n " << n << " levels " << stats->n << " strategies " << m << '\n';
  f << "incumbent " << s->obj << ' ' << s->ifs_found;
  for (int i = 0; i < n; i++) f << ' ' << s->perm[i];
  f << '\n';
  f << "root_bound " << stats->root_bound << '\n';
  f << "counters " << stats->node << ' ' << stats->fathom << ' '
    << stats->enumerate << ' ' << stats->fw_iter << '\n';
  f << "times " << stats->total_time << ' ' << stats->bound_time << '\n';
  f << "settings " << settings_record() << '\n';   // version 8: the hash of the policy settings and those set away from their defaults
  write_array(f, "n_level", stats->n_level, L);
  write_array(f, "f_level", stats->f_level, L);
  write_array(f, "l_level", stats->l_level, L);
  write_array(f, "e_level", stats->e_level, L);
  write_array(f, "t_level", stats->t_level, L);
  f << "n_strat";
  for (int i = 0; i < L; i++)
    for (int j = 0; j < m; j++) f << ' ' << stats->n_strat[i][j];
  f << '\n';
  write_array(f, "t_strat", stats->t_strat, m);
#ifdef RECORD_STDDEV
  write_array(f, "rg_sum", stats->rg_sum, L);
  write_array(f, "rg_sumsq", stats->rg_sumsq, L);
#endif
  f << "dnn_counts " << stats->dnn_nodes << ' ' << stats->dnn_fathoms << '\n';
  f << "dnn_time " << stats->dnn_time << '\n';
  {
    QAPStatistics::ctr_type *dc[] = {stats->pop_f_level, stats->dnn_n_level, stats->dnn_f_level,
                                     stats->dnn_g_level, stats->dnn_i_level};
    for (int k = 0; k < 5; k++) write_array(f, DNN_LEVEL_CTRS[k], dc[k], L);
    write_array(f, "dnn_t_level", stats->dnn_t_level, L);
    write_array(f, CUT_LEVEL_CTRS[0], stats->cut_n_level, L);
    write_array(f, CUT_LEVEL_CTRS[1], stats->cut_f_level, L);
    write_array(f, "cut_t_level", stats->cut_t_level, L);
  }

  // the parents' X the open nodes start from, each written once
  std::map<const WarmStart *, int> warm_id;
  std::vector<const WarmStart *> warms;
  for (Node *p = stack; p != NULL; p = p->link) {
    const WarmStart *w = p->data.warm;
    if (w && warm_id.find(w) == warm_id.end()) {
      warm_id[w] = (int)warms.size();
      warms.push_back(w);
    }
  }
  f << "warm_starts " << warms.size() << '\n';
  for (size_t k = 0; k < warms.size(); k++) {
    const WarmStart *w = warms[k];
    f << "warm " << w->X->m;
    for (int i = 0; i < n; i++) f << ' ' << w->p[i];
    for (unsigned int i = 0; i < w->X->m; i++)
      for (unsigned int j = 0; j < w->X->n; j++) f << ' ' << w->X->me[i][j];
    f << ' ' << w->dnn_ext;   // version 4: QAP_DNN_EXT levels left below a last-level miss
    f << ' ' << (w->dnn_pool.empty() ? std::string("-") : w->dnn_pool);   // version 5: the cut pool
#ifdef HAVE_DNN
    f << ' ' << dnn_field(w, file);   // version 6: the parent's DNN state (ship.bits); version 9: or @path into the store
#else
    f << " -";
#endif
    f << '\n';
  }

  f << "nodes " << list_length(stack) << '\n';
  for (Node *p = stack; p != NULL; p = p->link) {
    const QAPAssignment &a = p->data;
    f << "node " << a.depth << ' ' << a.bound << ' ' << a.predicted_bound
      << ' ' << a.nfix;
    for (int i = 0; i < n; i++) f << ' ' << a.p[i];
    for (int i = 0; i < n; i++) f << ' ' << a.history[i];
    f << ' ' << (a.warm ? warm_id[a.warm] : -1);
    f << '\n';
  }
  f << "end\n";
  f.close();
  if (!f) {
    std::cerr << "write_state: error writing [" << tmp << "]" << std::endl;
    return false;
  }
  if (std::rename(tmp.c_str(), file) != 0) {
    std::cerr << "write_state: cannot rename to [" << file << "]" << std::endl;
    return false;
  }
  return true;
}

static bool fail(const char *file, const char *what)
{
  std::cerr << "read_state: [" << file << "] " << what << std::endl;
  return false;
}

bool read_state(const char *file, int n, Node *&stack, QAPSolution *s,
		QAPStatistics *stats)
{
  std::ifstream f(file);
  if (!f) return fail(file, "cannot open");

  std::string tag;
  int version, fn, fL, fm;
  if (!(f >> tag >> version) || tag != STATE_MAGIC || version < 1 ||
      version > STATE_VERSION)
    return fail(file, "not a version 1-8 state file");
  std::string t1, t2, t3;
  if (!(f >> t1 >> fn >> t2 >> fL >> t3 >> fm) || t1 != "n")
    return fail(file, "bad header");
  if (fn != n) return fail(file, "was written for a different instance size");
  if (fL != stats->n || fm != stats->m)
    return fail(file, "was written with a different number of levels or strategies");

  const int L = stats->n + 1, m = stats->m;
  if (!(f >> tag >> s->obj >> s->ifs_found) || tag != "incumbent")
    return fail(file, "bad incumbent");
  for (int i = 0; i < n; i++)
    if (!(f >> s->perm[i])) return fail(file, "bad incumbent");
  if (!(f >> tag >> stats->root_bound) || tag != "root_bound")
    return fail(file, "bad root bound");
  if (!(f >> tag >> stats->node >> stats->fathom >> stats->enumerate
	>> stats->fw_iter) || tag != "counters")
    return fail(file, "bad counters");
  if (!(f >> tag >> stats->total_time >> stats->bound_time) || tag != "times")
    return fail(file, "bad times");
  if (version >= 8) {
    std::string rec;
    if (!(f >> tag) || tag != "settings" || !std::getline(f, rec)) return fail(file, "bad settings record");
    rec.erase(0, rec.find_first_not_of(" \t"));
    std::string why;
    if (!settings_check_record(rec, why)) {
      if (!state_allow_settings_change) return fail(file, (why + "; --allow-settings-change to resume anyway").c_str());
      std::cerr << "read_state: [" << file << "] warning: " << why << " (allowed)" << std::endl;
    }
  }
  if (!read_array(f, "n_level", stats->n_level, L) ||
      !read_array(f, "f_level", stats->f_level, L) ||
      !read_array(f, "l_level", stats->l_level, L) ||
      !read_array(f, "e_level", stats->e_level, L) ||
      !read_array(f, "t_level", stats->t_level, L))
    return fail(file, "bad per-level statistics");
  if (!(f >> tag) || tag != "n_strat") return fail(file, "bad n_strat");
  for (int i = 0; i < L; i++)
    for (int j = 0; j < m; j++)
      if (!(f >> stats->n_strat[i][j])) return fail(file, "bad n_strat");
  if (!read_array(f, "t_strat", stats->t_strat, m))
    return fail(file, "bad t_strat");
#ifdef RECORD_STDDEV
  if (!read_array(f, "rg_sum", stats->rg_sum, L) ||
      !read_array(f, "rg_sumsq", stats->rg_sumsq, L))
    return fail(file, "bad gap statistics");
#endif
  if (version >= 3 &&
      (!(f >> tag >> stats->dnn_nodes >> stats->dnn_fathoms) || tag != "dnn_counts" ||
       !(f >> tag >> stats->dnn_time) || tag != "dnn_time"))
    return fail(file, "bad DNN counts");
  {
    QAPStatistics::ctr_type *dc[] = {stats->pop_f_level, stats->dnn_n_level, stats->dnn_f_level,
                                     stats->dnn_g_level, stats->dnn_i_level};
    if (version >= 7) {
      for (int k = 0; k < 5; k++)
	if (!read_array(f, DNN_LEVEL_CTRS[k], dc[k], L)) return fail(file, "bad per-level DNN counts");
      if (!read_array(f, "dnn_t_level", stats->dnn_t_level, L) ||
	  !read_array(f, CUT_LEVEL_CTRS[0], stats->cut_n_level, L) ||
	  !read_array(f, CUT_LEVEL_CTRS[1], stats->cut_f_level, L) ||
	  !read_array(f, "cut_t_level", stats->cut_t_level, L))
	return fail(file, "bad per-level DNN counts");
      // the cut phase's totals were not in earlier versions: rebuilt here
      stats->dnn_cut_nodes = stats->dnn_cut_fathoms = 0;
      stats->dnn_cut_time = 0.0;
      for (int i = 0; i < L; i++) {
	stats->dnn_cut_nodes += stats->cut_n_level[i];
	stats->dnn_cut_fathoms += stats->cut_f_level[i];
	stats->dnn_cut_time += stats->cut_t_level[i];
      }
    } else {            // earlier versions: no causes, the arrays start at 0
      for (int k = 0; k < 5; k++)
	for (int i = 0; i < L; i++) dc[k][i] = 0;
      for (int i = 0; i < L; i++) {
	stats->dnn_t_level[i] = stats->cut_t_level[i] = 0.0;
	stats->cut_n_level[i] = stats->cut_f_level[i] = 0;
      }
    }
  }

  std::vector<WarmStart *> warms;
  if (version >= 2) {
    int nw;
    if (!(f >> tag >> nw) || tag != "warm_starts" || nw < 0)
      return fail(file, "bad warm start count");
    for (int k = 0; k < nw; k++) {
      int wm;
      if (!(f >> tag >> wm) || tag != "warm" || wm < 1 || wm > n)
	return fail(file, "bad warm start");
      Index *wp = new Index[n];
      MAT *X = m_get(wm, wm);
      bool ok = true;
      for (int i = 0; ok && i < n; i++) ok = (bool)(f >> wp[i]);
      for (int i = 0; ok && i < wm; i++)
	for (int j = 0; ok && j < wm; j++) ok = (bool)(f >> X->me[i][j]);
      int ext = 0;
      std::string pool;
      if (ok && version >= 4) ok = (bool)(f >> ext);
      if (ok && version >= 5) ok = (bool)(f >> pool);
      std::string dnn = "-";
      if (ok && version >= 6) ok = (bool)(f >> dnn);
#ifdef HAVE_DNN
      DNNWarm *dw = NULL;
      std::string ref;
      if (ok && dnn != "-") {
	std::string tok;
	ok = dnn_field_token(dnn, file, tok, ref) && (dw = dnn_from_token(tok)) != NULL;
      }
#endif
      if (ok) {
	warms.push_back(new_warm_start(wp, n, X));
	warms.back()->dnn_ext = ext;
	if (pool != "-") warms.back()->dnn_pool = pool;
#ifdef HAVE_DNN
	if (dw) { warms.back()->dnn = dw; warms.back()->dnn_free = free_dnn_state; warms.back()->dnn_ref = ref; }
#endif
      }
      delete [] wp;
      M_FREE(X);
      if (!ok) {
	for (size_t i = 0; i < warms.size(); i++) warm_start_release(warms[i]);
	return fail(file, "bad warm start");
      }
    }
  }

  int count;
  if (!(f >> tag >> count) || tag != "nodes" || count < 0)
    count = -1;

  // read in file order, then insert from the back so the stack order is kept
  const char *err = (count < 0) ? "bad node count" : NULL;
  QAPAssignment *nodes = new QAPAssignment[count > 0 ? count : 1];
  for (int k = 0; !err && k < count; k++) {
    QAPAssignment a(n);
    if (!(f >> tag >> a.depth >> a.bound >> a.predicted_bound >> a.nfix) ||
	tag != "node") {
      err = "bad node";
      break;
    }
    for (int j = 0; j < n; j++) a.pinv[j] = NONE;
    int fixed = 0;
    for (int i = 0; !err && i < n; i++) {
      if (!(f >> a.p[i]) || a.p[i] < NONE || a.p[i] >= n)
	err = "bad node assignment";
      else if (a.p[i] != NONE) { a.pinv[a.p[i]] = i; fixed++; }
    }
    for (int i = 0; !err && i < n; i++)
      if (!(f >> a.history[i])) err = "bad node history";
    if (!err && fixed != a.nfix) err = "node nfix does not match";
    if (!err && version >= 2) {
      int w;
      if (!(f >> w) || w < -1 || w >= (int)warms.size())
	err = "bad node warm start";
      else if (w >= 0)
	a.set_warm(warms[w]);
    }
    if (!err) nodes[k] = a;
  }
  if (!err && (!(f >> tag) || tag != "end")) err = "missing end";

  if (!err) {
    for (int k = count - 1; k >= 0; k--)
      list_head_insert(stack, nodes[k]);
    stats->stack_size = count;
  }
  delete [] nodes;
  // the nodes hold their own references now
  for (size_t k = 0; k < warms.size(); k++) warm_start_release(warms[k]);
  return err ? fail(file, err) : true;
}
