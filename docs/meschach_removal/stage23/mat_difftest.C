// Differential test of latest/mat.C against Meschach (stage 2 of
// docs/meschach_removal.md). Runs the same random sequence of m_get,
// m_resize, m_copy/_m_copy, m_zero, M_FREE, v_get, v_resize, v_copy, in_prod,
// V_FREE and hand-lowered dims on pools of matrices and vectors, and prints the
// whole visible state (sizes, max_* fields, row layout, every entry in hex)
// after each operation. Build it twice and diff the outputs:
//
//   M=<a checkout with src/MESCHACH built, e.g. afe99db>
//   g++ -O2 -DUSE_MESCHACH -I$M/src/MESCHACH mat_difftest.C $M/src/MESCHACH/meschach.a -lm -o t_mes
//   g++ -O2 -I latest mat_difftest.C latest/mat.C -o t_mat
//   ./t_mes 1 200000 > a.txt; ./t_mat 1 200000 > b.txt; cmp a.txt b.txt
//
// Entries beyond a matrix's or vector's visible size are not compared: after
// a realloc they are uninitialized in both.
//
// One case is skipped in both builds: resizing a matrix that never had a row
// (max_m == 0) to 0 rows and more columns. Meschach's m_resize then zeroes
// "the rest of row 0" outside its storage (valgrind: invalid write after a
// block of size 0) and corrupts the heap; mat.C skips that write.
#ifdef USE_MESCHACH
extern "C" {
#include "matrix.h"
}
#else
#include "mat.h"
#endif
#include <cstdio>
#include <cstdlib>

static unsigned long long rng_state;
static unsigned rnd(unsigned k)  // 0..k-1
{
  rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
  return (unsigned)((rng_state >> 33) % k);
}
static double rval() { return ((int)rnd(2001) - 1000) / 7.0; }

const int NM = 6, NV = 6, MAXD = 7;
static MAT *M[NM];
static VEC *V[NV];

// would Meschach's m_resize(A, new_m, new_n) write outside A's storage?
static bool meschach_overrun(MAT *A, int new_m, int new_n)
{
  return A && A->max_m == 0 && new_m == 0 && new_n > (int)A->n;
}

static void dump_m(int k)
{
  MAT *A = M[k];
  if (!A) { printf("M%d null\n", k); return; }
  int layout = 1;
  for (unsigned i = 0; i < A->m; i++)
    if (A->me[i] != A->base + i * A->n) layout = 0;
  printf("M%d %ux%u max %u %u %u layout %d:", k, A->m, A->n, A->max_m, A->max_n, A->max_size, layout);
  for (unsigned i = 0; i < A->m; i++)
    for (unsigned j = 0; j < A->n; j++) printf(" %a", A->me[i][j]);
  printf("\n");
}
static void dump_v(int k)
{
  VEC *x = V[k];
  if (!x) { printf("V%d null\n", k); return; }
  printf("V%d %u max %u:", k, x->dim, x->max_dim);
  for (unsigned i = 0; i < x->dim; i++) printf(" %a", x->ve[i]);
  printf("\n");
}

int main(int argc, char **argv)
{
  rng_state = argc > 1 ? strtoull(argv[1], 0, 10) : 1;
  long steps = argc > 2 ? atol(argv[2]) : 100000;
  for (long s = 0; s < steps; s++) {
    int op = rnd(15), i = rnd(NM), j = rnd(NM), r = rnd(MAXD), c = rnd(MAXD);
    printf("%ld op %d i %d j %d r %d c %d\n", s, op, i, j, r, c);
    switch (op) {
      case 0: if (M[i]) M_FREE(M[i]); M[i] = m_get(r, c); dump_m(i); break;
      case 1: case 2:
        if (meschach_overrun(M[i], r, c)) { printf("skip\n"); break; }
        M[i] = m_resize(M[i], r, c); dump_m(i); break;
      case 3:
        if (M[i]) for (unsigned a = 0; a < M[i]->m; a++)
          for (unsigned b = 0; b < M[i]->n; b++) M[i]->me[a][b] = rval();
        dump_m(i); break;
      case 4:
        if (M[i] && M[i] != M[j] && meschach_overrun(M[j], M[i]->m, M[i]->n)) { printf("skip\n"); break; }
        if (M[i]) { M[j] = m_copy(M[i], M[j]); dump_m(j); }
        break;
      case 5:
        if (M[i] && M[i] != M[j] && meschach_overrun(M[j], M[i]->m, M[i]->n)) { printf("skip\n"); break; }
        if (M[i]) {
          unsigned i0 = rnd(M[i]->m + 1), j0 = rnd(M[i]->n + 1);
          printf("i0 %u j0 %u\n", i0, j0);
          M[j] = _m_copy(M[i], M[j], i0, j0); dump_m(j);
        }
        break;
      case 6: if (M[i]) { m_zero(M[i]); dump_m(i); } break;
      case 7: if (rnd(3) == 0) { M_FREE(M[i]); dump_m(i); } break;
      // vectors (i, j taken mod NV)
      case 8: i %= NV; if (V[i]) V_FREE(V[i]); V[i] = v_get(r); dump_v(i); break;
      case 9: case 10: {
        i %= NV;
        VEC *old = V[i], *x = v_resize(V[i], r);
        // a v_get(0) vector is not resized: a new one comes back
        printf("same %d\n", x == old);
        if (old && x != old) v_free(old);
        V[i] = x; dump_v(i); break;
      }
      case 11:
        i %= NV;
        if (V[i]) for (unsigned a = 0; a < V[i]->dim; a++) V[i]->ve[a] = rval();
        dump_v(i); break;
      case 12:
        i %= NV; j %= NV;
        if (V[i]) {
          if (rnd(2)) V[j] = v_copy(V[i], V[j]);
          else { unsigned i0 = rnd(V[i]->dim + 1); printf("i0 %u\n", i0); V[j] = _v_copy(V[i], V[j], i0); }
          dump_v(j);
        }
        break;
      case 13:
        i %= NV; j %= NV;
        if (V[i] && V[j]) {
          unsigned lim = V[i]->dim < V[j]->dim ? V[i]->dim : V[j]->dim;
          unsigned i0 = rnd(2) ? 0 : rnd(lim + 1);
          printf("in_prod %u %a\n", i0, _in_prod(V[i], V[j], i0));
        }
        break;
      case 14:
        i %= NV;
        if (V[i] && V[i]->dim > 0 && rnd(2)) V[i]->dim--;  // as qpbnd.C does
        else if (rnd(3) == 0) V_FREE(V[i]);
        dump_v(i); break;
    }
  }
  for (int k = 0; k < NM; k++) M_FREE(M[k]);
  for (int k = 0; k < NV; k++) V_FREE(V[k]);
  printf("done\n");
  return 0;
}
