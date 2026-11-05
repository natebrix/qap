/*
  Timing wrappers for the Meschach functions the solver calls.

  Linked into a scratch copy of bnbdrv with -Wl,--wrap=<symbol> for each
  symbol below (see build_timer.sh). Each wrapper reads the TSC around the
  real call; only the outermost Meschach call is timed (a Meschach routine
  that calls another wrapped one is not counted twice). my_svd (latest/
  linalg.C, the symmetric eigendecomposition built on Hfactor/makeHQ/trieig)
  is timed separately via mes_svd_begin/mes_svd_end, with a histogram by n.

  With MES_DUMP=<file>, every MES_DUMP_EVERY-th (default 500) my_svd input
  matrix is appended to <file> as: int32 n, then n*n doubles (row major).
  At exit a table is printed to stderr: calls, seconds, ns/call per symbol,
  and the my_svd histogram. Seconds come from TSC ticks calibrated against
  CLOCK_MONOTONIC over the whole run.
*/
#include <x86intrin.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include "meschach.h"
#undef catch

namespace {
struct Slot { const char *name; unsigned long long calls, ticks; };
enum { S_M_GET, S_M_FREE, S_M_RESIZE, S_M_COPY, S_M_ZERO, S_M_TRANSP,
       S_V_GET, S_V_FREE, S_V_RESIZE, S_V_COPY, S_IN_PROD,
       S_PX_GET, S_PX_FREE, S_PX_ROWS, S_HFACTOR, S_MAKEHQ, S_TRIEIG, S_N };
Slot slots[S_N] = {
  {"m_get"}, {"m_free"}, {"m_resize"}, {"_m_copy"}, {"m_zero"}, {"m_transp"},
  {"v_get"}, {"v_free"}, {"v_resize"}, {"_v_copy"}, {"_in_prod"},
  {"px_get"}, {"px_free"}, {"px_rows"}, {"Hfactor"}, {"makeHQ"}, {"trieig"}};
int depth = 0;
const int NMAX = 128;
unsigned long long svd_calls[NMAX], svd_ticks[NMAX];
unsigned long long svd_t0;
double mono0; unsigned long long tsc0;
FILE *dumpf = nullptr; long dump_every = 500; unsigned long long svd_seen = 0;

double mono() { timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

struct Init {
  Init() {
    mono0 = mono(); tsc0 = __rdtsc();
    if (const char *f = getenv("MES_DUMP")) dumpf = fopen(f, "wb");
    if (const char *e = getenv("MES_DUMP_EVERY")) dump_every = atol(e);
  }
  ~Init() {
    double secs = mono() - mono0; unsigned long long ticks = __rdtsc() - tsc0;
    double sec_per_tick = secs / (double)ticks;
    unsigned long long tot = 0;
    fprintf(stderr, "\n==== Meschach timing (TSC %.3f GHz, run %.2f s wall) ====\n", 1e-9 / sec_per_tick, secs);
    fprintf(stderr, "%-10s %14s %12s %10s\n", "symbol", "calls", "seconds", "ns/call");
    for (auto &s : slots) {
      if (!s.calls) continue;
      tot += s.ticks;
      fprintf(stderr, "%-10s %14llu %12.3f %10.1f\n", s.name, s.calls, s.ticks * sec_per_tick,
              s.ticks * sec_per_tick * 1e9 / s.calls);
    }
    fprintf(stderr, "%-10s %14s %12.3f\n", "TOTAL", "", tot * sec_per_tick);
    unsigned long long sc = 0, st = 0;
    for (int n = 0; n < NMAX; n++) { sc += svd_calls[n]; st += svd_ticks[n]; }
    fprintf(stderr, "my_svd (eigendecomposition incl. its Meschach calls and copies): %llu calls, %.3f s\n",
            sc, st * sec_per_tick);
    fprintf(stderr, "%4s %12s %10s %10s\n", "n", "calls", "seconds", "us/call");
    for (int n = 0; n < NMAX; n++)
      if (svd_calls[n])
        fprintf(stderr, "%4d %12llu %10.3f %10.2f\n", n, svd_calls[n], svd_ticks[n] * sec_per_tick,
                svd_ticks[n] * sec_per_tick * 1e6 / svd_calls[n]);
    if (dumpf) fclose(dumpf);
  }
} init_;

struct Timer {
  int s; unsigned long long t0;
  explicit Timer(int s_) : s(s_) { if (depth++ == 0) t0 = __rdtsc(); }
  ~Timer() { if (--depth == 0) { slots[s].ticks += __rdtsc() - t0; slots[s].calls++; } }
};
}  // namespace

extern "C" {
MAT *__real_m_get(int, int);
int __real_m_free(MAT *);
MAT *__real_m_resize(MAT *, int, int);
MAT *__real__m_copy(MAT *, MAT *, u_int, u_int);
MAT *__real_m_zero(MAT *);
MAT *__real_m_transp(MAT *, MAT *);
VEC *__real_v_get(int);
int __real_v_free(VEC *);
VEC *__real_v_resize(VEC *, int);
VEC *__real__v_copy(VEC *, VEC *, u_int);
double __real__in_prod(VEC *, VEC *, u_int);
PERM *__real_px_get(int);
int __real_px_free(PERM *);
MAT *__real_px_rows(PERM *, MAT *, MAT *);
MAT *__real_Hfactor(MAT *, VEC *, VEC *);
MAT *__real_makeHQ(MAT *, VEC *, VEC *, MAT *);
VEC *__real_trieig(VEC *, VEC *, MAT *);

MAT *__wrap_m_get(int m, int n) { Timer t(S_M_GET); return __real_m_get(m, n); }
int __wrap_m_free(MAT *a) { Timer t(S_M_FREE); return __real_m_free(a); }
MAT *__wrap_m_resize(MAT *a, int m, int n) { Timer t(S_M_RESIZE); return __real_m_resize(a, m, n); }
MAT *__wrap__m_copy(MAT *a, MAT *b, u_int i, u_int j) { Timer t(S_M_COPY); return __real__m_copy(a, b, i, j); }
MAT *__wrap_m_zero(MAT *a) { Timer t(S_M_ZERO); return __real_m_zero(a); }
MAT *__wrap_m_transp(MAT *a, MAT *b) { Timer t(S_M_TRANSP); return __real_m_transp(a, b); }
VEC *__wrap_v_get(int n) { Timer t(S_V_GET); return __real_v_get(n); }
int __wrap_v_free(VEC *v) { Timer t(S_V_FREE); return __real_v_free(v); }
VEC *__wrap_v_resize(VEC *v, int n) { Timer t(S_V_RESIZE); return __real_v_resize(v, n); }
VEC *__wrap__v_copy(VEC *a, VEC *b, u_int i) { Timer t(S_V_COPY); return __real__v_copy(a, b, i); }
double __wrap__in_prod(VEC *a, VEC *b, u_int i) { Timer t(S_IN_PROD); return __real__in_prod(a, b, i); }
PERM *__wrap_px_get(int n) { Timer t(S_PX_GET); return __real_px_get(n); }
int __wrap_px_free(PERM *p) { Timer t(S_PX_FREE); return __real_px_free(p); }
MAT *__wrap_px_rows(PERM *p, MAT *a, MAT *o) { Timer t(S_PX_ROWS); return __real_px_rows(p, a, o); }
MAT *__wrap_Hfactor(MAT *a, VEC *d, VEC *b) { Timer t(S_HFACTOR); return __real_Hfactor(a, d, b); }
MAT *__wrap_makeHQ(MAT *h, VEC *d, VEC *b, MAT *q) { Timer t(S_MAKEHQ); return __real_makeHQ(h, d, b, q); }
VEC *__wrap_trieig(VEC *a, VEC *b, MAT *q) { Timer t(S_TRIEIG); return __real_trieig(a, b, q); }
}

// called from my_svd (linalg.C in the scratch copy)
void mes_svd_begin(MAT *A)
{
  if (dumpf && (svd_seen++ % dump_every) == 0) {
    int n = A->m;
    fwrite(&n, sizeof n, 1, dumpf);
    for (int i = 0; i < n; i++) fwrite(A->me[i], sizeof(double), n, dumpf);
  }
  svd_t0 = __rdtsc();
}
void mes_svd_end(int n)
{
  unsigned long long dt = __rdtsc() - svd_t0;
  if (n < NMAX) { svd_calls[n]++; svd_ticks[n] += dt; }
}
