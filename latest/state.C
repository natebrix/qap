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
static const int STATE_VERSION = 2;   // 1: no warm starts (still read)

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
    return fail(file, "not a version 1 or 2 state file");
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
      if (ok) warms.push_back(new_warm_start(wp, n, X));
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
