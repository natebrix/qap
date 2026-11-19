/* ls.c: 2-exchange descent for the QAP, cost = sum_{i,k} a[i][k] b[p[i]][p[k]] (general matrices).
   Delta evaluation as in Taillard's robust taboo search: the full delta table once (O(n^3)), then
   O(1) per pair not touching the exchanged rows after each move, O(n) for those that do: O(n^2) a move.
   Best-improvement descent over the free rows (free[i] != 0). Build: gcc -O3 -march=native -shared -fPIC.
   Exported:
     double qap_cost(int n, const double *a, const double *b, const int *p)
     double descent(int n, a, b, int *p, const int *fr, long *moves)  -> final cost; p updated in place
     double multistart(int n, a, b, int *best, const int *fr, int nstart, unsigned seed, long *moves)
       -> descent from best (the given start), then nstart random completions of the free rows
          (each a random permutation of the free rows' values), the best kept in best. */
#include <stdlib.h>
#include <string.h>
#define A(i, j) a[(i) * n + (j)]
#define B(i, j) b[(i) * n + (j)]

double qap_cost(int n, const double *a, const double *b, const int *p) {
    double c = 0;
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) c += A(i, k) * B(p[i], p[k]);
    return c;
}

static double delta_full(int n, const double *a, const double *b, const int *p, int i, int j) {
    int pi = p[i], pj = p[j];
    double d = (A(i, i) - A(j, j)) * (B(pj, pj) - B(pi, pi)) + (A(i, j) - A(j, i)) * (B(pj, pi) - B(pi, pj));
    for (int k = 0; k < n; k++) {
        if (k == i || k == j) continue;
        int pk = p[k];
        d += (A(k, i) - A(k, j)) * (B(pk, pj) - B(pk, pi)) + (A(i, k) - A(j, k)) * (B(pj, pk) - B(pi, pk));
    }
    return d;
}

/* after rows r, s were exchanged (p already updated): the delta of pair (i, j), disjoint from {r, s} */
static double delta_part(int n, const double *a, const double *b, const int *p, const double *D, int i, int j, int r,
                         int s) {
    int pi = p[i], pj = p[j], pr = p[r], ps = p[s];
    return D[i * n + j] + (A(r, i) - A(r, j) + A(s, j) - A(s, i)) * (B(ps, pi) - B(ps, pj) + B(pr, pj) - B(pr, pi)) +
           (A(i, r) - A(j, r) + A(j, s) - A(i, s)) * (B(pi, ps) - B(pj, ps) + B(pj, pr) - B(pi, pr));
}

double descent(int n, const double *a, const double *b, int *p, const int *fr, long *moves) {
    int f[512], m = 0;
    for (int i = 0; i < n; i++)
        if (fr[i]) f[m++] = i;
    double c = qap_cost(n, a, b, p);
    if (m < 2) return c;
    double *D = malloc(sizeof(double) * n * n);
    for (int x = 0; x < m; x++)
        for (int y = x + 1; y < m; y++) D[f[x] * n + f[y]] = delta_full(n, a, b, p, f[x], f[y]);
    for (;;) {
        double best = -0.5;
        int bi = -1, bj = -1;
        for (int x = 0; x < m; x++)
            for (int y = x + 1; y < m; y++) {
                double d = D[f[x] * n + f[y]];
                if (d < best) best = d, bi = f[x], bj = f[y];
            }
        if (bi < 0) break;
        int t = p[bi]; p[bi] = p[bj]; p[bj] = t;
        c += best;
        if (moves) (*moves)++;
        for (int x = 0; x < m; x++)
            for (int y = x + 1; y < m; y++) {
                int i = f[x], j = f[y];
                if (i == bi || i == bj || j == bi || j == bj) D[i * n + j] = delta_full(n, a, b, p, i, j);
                else D[i * n + j] = delta_part(n, a, b, p, D, i, j, bi, bj);
            }
    }
    free(D);
    return c;
}

static unsigned rng_state;
static unsigned rnd(void) { rng_state = rng_state * 1664525u + 1013904223u; return rng_state >> 8; }

double multistart(int n, const double *a, const double *b, int *best, const int *fr, int nstart, unsigned seed,
                  long *moves) {
    double cb = descent(n, a, b, best, fr, moves);
    int *p = malloc(sizeof(int) * n), f[512], m = 0;
    for (int i = 0; i < n; i++)
        if (fr[i]) f[m++] = i;
    rng_state = seed * 2654435761u + 12345u;
    for (int s = 0; s < nstart; s++) {
        memcpy(p, best, sizeof(int) * n);
        for (int x = m - 1; x > 0; x--) {
            int y = rnd() % (x + 1), t = p[f[x]];
            p[f[x]] = p[f[y]]; p[f[y]] = t;
        }
        double c = descent(n, a, b, p, fr, moves);
        if (c < cb) cb = c, memcpy(best, p, sizeof(int) * n);
    }
    free(p);
    return cb;
}
