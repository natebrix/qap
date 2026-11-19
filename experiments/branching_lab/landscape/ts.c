/* ts.c: a robust taboo search (after Taillard 1991) over the free rows, cost = sum a[i][k] b[p[i]][p[k]].
   double rots(int n, a, b, int *p, const int *fr, int iters, unsigned seed) -> best cost; p <- best permutation.
   Delta table as in ls.c (O(n^2) a move). Tabu: facility i may not return to location l before tabu[i][l];
   tenure drawn in [0.9 m, 1.1 m] every 2 m iterations; aspiration when a move beats the best so far, or when
   both locations have not been used by those facilities for 5 m^2 iterations. Separate from ls.c so a rebuild
   never touches the ls.so a running process has mapped. */
#include <stdlib.h>
#include <string.h>
#define A(i, j) a[(i) * n + (j)]
#define B(i, j) b[(i) * n + (j)]

static double cost(int n, const double *a, const double *b, const int *p) {
    double c = 0;
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) c += A(i, k) * B(p[i], p[k]);
    return c;
}
static double dfull(int n, const double *a, const double *b, const int *p, int i, int j) {
    int pi = p[i], pj = p[j];
    double d = (A(i, i) - A(j, j)) * (B(pj, pj) - B(pi, pi)) + (A(i, j) - A(j, i)) * (B(pj, pi) - B(pi, pj));
    for (int k = 0; k < n; k++) {
        if (k == i || k == j) continue;
        int pk = p[k];
        d += (A(k, i) - A(k, j)) * (B(pk, pj) - B(pk, pi)) + (A(i, k) - A(j, k)) * (B(pj, pk) - B(pi, pk));
    }
    return d;
}
static double dpart(int n, const double *a, const double *b, const int *p, double d, int i, int j, int r, int s) {
    int pi = p[i], pj = p[j], pr = p[r], ps = p[s];
    return d + (A(r, i) - A(r, j) + A(s, j) - A(s, i)) * (B(ps, pi) - B(ps, pj) + B(pr, pj) - B(pr, pi)) +
           (A(i, r) - A(j, r) + A(j, s) - A(i, s)) * (B(pi, ps) - B(pj, ps) + B(pj, pr) - B(pi, pr));
}
static unsigned st;
static unsigned rnd(void) { st = st * 1664525u + 1013904223u; return st >> 8; }

double rots(int n, const double *a, const double *b, int *p, const int *fr, int iters, unsigned seed) {
    int f[512], m = 0;
    for (int i = 0; i < n; i++)
        if (fr[i]) f[m++] = i;
    double c = cost(n, a, b, p);
    if (m < 2) return c;
    st = seed * 2654435761u + 7u;
    double *D = malloc(sizeof(double) * n * n);
    long *tabu = calloc((size_t)n * n, sizeof(long));
    int *bp = malloc(sizeof(int) * n);
    memcpy(bp, p, sizeof(int) * n);
    double bc = c;
    for (int x = 0; x < m; x++)
        for (int y = x + 1; y < m; y++) D[f[x] * n + f[y]] = dfull(n, a, b, p, f[x], f[y]);
    int lo = (int)(0.9 * m), hi = (int)(1.1 * m) + 1, ten = lo + rnd() % (hi - lo + 1);
    long asp = 5L * m * m;
    for (long it = 1; it <= iters; it++) {
        if (it % (2 * m) == 0) ten = lo + rnd() % (hi - lo + 1);
        double bd = 1e300; int bi = -1, bj = -1, forced = 0;
        for (int x = 0; x < m; x++)
            for (int y = x + 1; y < m; y++) {
                int i = f[x], j = f[y];
                double d = D[i * n + j];
                int istabu = tabu[i * n + p[j]] >= it && tabu[j * n + p[i]] >= it;
                int aspired = c + d < bc - 0.5 ||
                              (tabu[i * n + p[j]] < it - asp && tabu[j * n + p[i]] < it - asp);
                if (aspired && !forced) { forced = 1; bd = d; bi = i; bj = j; continue; }
                if ((aspired && forced && d < bd) || (!forced && !istabu && d < bd)) bd = d, bi = i, bj = j;
            }
        if (bi < 0) continue;
        int t = p[bi]; p[bi] = p[bj]; p[bj] = t;
        c += bd;
        tabu[bi * n + p[bj]] = it + ten;      /* bi may not go back to its old location (now p[bj]) */
        tabu[bj * n + p[bi]] = it + ten;
        if (c < bc - 0.5) bc = c, memcpy(bp, p, sizeof(int) * n);
        for (int x = 0; x < m; x++)
            for (int y = x + 1; y < m; y++) {
                int i = f[x], j = f[y];
                if (i == bi || i == bj || j == bi || j == bj) D[i * n + j] = dfull(n, a, b, p, i, j);
                else D[i * n + j] = dpart(n, a, b, p, D[i * n + j], i, j, bi, bj);
            }
    }
    memcpy(p, bp, sizeof(int) * n);
    free(D); free(tabu); free(bp);
    return bc;
}
