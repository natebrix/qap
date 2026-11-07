// P-step and dense-term timing, side 0 against side 1 (README.md, "Side 1"), on node
// 13's region layout (wil50's 5 x 10 grid, m = 47), 1,000 random triangles, the caches
// flushed before every call. Includes dnncut.C for the projector (anonymous namespace).
// Arms: 0 side 0; 1 side 1 as built (pstep_rows1: row buffers and the next row
// prefetched); 2-6 the side-1 loops tried before it, all with the same Y and moments:
// var 1 row buffers without the prefetch (12.4 ms), var 2 each row's m x m block
// transposed to [b][j] and side 0's contiguous loop (15.3 ms), var 3 the buffers and
// the prefetch (as built), var 4 the buffers with the moment gather fused into the
// strided pass (19.8 ms), var 5 both. The first build's in-place loop (16.2 ms) is
// pstep_rows1 at commit 02655b9.
//   bash build_pstep_bench.sh && ./pstep_bench [reps]
#include "dnncut.C"


// side-1 row loops, variants: 1 = B (stage the row contiguously, compute strided in L1),
// 2 = T (transpose the row's m x m block to [b][j], side 0's contiguous inner loop)
static void rows1_variant(CutProjector &pj, const double *X, double *Y, const unsigned char *elim, int var)
{
  const int m = pj.m, n1 = pj.n1;
  const size_t mm = (size_t)m * m, sm = (size_t)m;
  auto clip = [](double v) { return std::min(1.0, std::max(0.0, v)); };
  static std::vector<double> bx, by;
  bx.resize(n1); by.resize(n1);
  for (int i = 0; i < m; i++)
    for (int a = 0; a < m; a++) {
      const size_t r = 1 + (size_t)i * m + a;
      const double *xr = X + r * n1;
      double *yr = Y + r * n1;
      if (var >= 3) {
        const bool pf = var == 3 || var == 5, fuse = var >= 4;
        const double *xn = (r + 1 < (size_t)n1) ? X + (r + 1) * n1 : xr;
        double *yn = (r + 1 < (size_t)n1) ? Y + (r + 1) * n1 : yr;
        std::copy(xr, xr + n1, bx.begin());
        const double *xx = bx.data() + 1; double *yy = by.data() + 1;
        const int per = (n1 + 8 * m - 1) / (8 * m);   // cache lines of the next row per b
        for (int b = 0; b < m; b++) {
          if (pf) for (int k = 0; k < per; k++) { const size_t o = ((size_t)b * per + k) * 8; if (o < (size_t)n1) { __builtin_prefetch(xn + o, 0, 3); __builtin_prefetch(yn + o, 1, 3); } }
          const int g = pj.sig2[(size_t)a * m + b];
          const double *xc = xx + b; double *yc = yy + b;
          if (b == a) {
            for (int j = 0; j < m; j++) yc[j * sm] = 0;
            yc[i * sm] = clip(g >= 0 ? xc[i * sm] - pj.Dsig[(size_t)g * m + i] : xc[i * sm]);
            if (fuse && g >= 0) pj.Gdsig[(size_t)g * m + i] += yc[i * sm];
          } else {
            if (g >= 0) {
              const double *cs = &pj.Csig[g * mm + (size_t)i * m];
              if (fuse) { double *gs = &pj.Gsig[g * mm + (size_t)i * m];
                for (int j = 0; j < m; j++) { const double v = j == i ? 0.0 : clip(xc[j * sm] - cs[j]); yc[j * sm] = v; gs[j] += v; } }
              else for (int j = 0; j < m; j++) yc[j * sm] = clip(xc[j * sm] - cs[j]);
            }
            else for (int j = 0; j < m; j++) yc[j * sm] = clip(xc[j * sm]);
            yc[i * sm] = 0;
          }
        }
        by[0] = clip(bx[0]);
        if (!fuse)
          for (int b = 0; b < m; b++) {
            const int g = pj.sig2[(size_t)a * m + b];
            if (g < 0) continue;
            const double *yc = yy + b;
            if (b == a) pj.Gdsig[(size_t)g * m + i] += yc[i * sm];
            else { double *gs = &pj.Gsig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) gs[j] += yc[j * sm]; }
          }
        std::copy(by.begin(), by.end(), yr);
      } else if (var == 1) {
        std::copy(xr, xr + n1, bx.begin());
        const double *xx = bx.data() + 1; double *yy = by.data() + 1;
        for (int b = 0; b < m; b++) {
          const int g = pj.sig2[(size_t)a * m + b];
          const double *xc = xx + b; double *yc = yy + b;
          if (b == a) {
            for (int j = 0; j < m; j++) yc[j * sm] = 0;
            yc[i * sm] = clip(g >= 0 ? xc[i * sm] - pj.Dsig[(size_t)g * m + i] : xc[i * sm]);
          } else {
            if (g >= 0) { const double *cs = &pj.Csig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) yc[j * sm] = clip(xc[j * sm] - cs[j]); }
            else for (int j = 0; j < m; j++) yc[j * sm] = clip(xc[j * sm]);
            yc[i * sm] = 0;
          }
        }
        by[0] = clip(bx[0]);
        for (int b = 0; b < m; b++) {
          const int g = pj.sig2[(size_t)a * m + b];
          if (g < 0) continue;
          const double *yc = yy + b;
          if (b == a) pj.Gdsig[(size_t)g * m + i] += yc[i * sm];
          else { double *gs = &pj.Gsig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) gs[j] += yc[j * sm]; }
        }
        std::copy(by.begin(), by.end(), yr);
      } else {
        // bx[b m + j] = xr[1 + j m + b]
        for (int j = 0; j < m; j++) { const double *xj = xr + 1 + j * sm; for (int b = 0; b < m; b++) bx[b * sm + j] = xj[b]; }
        for (int b = 0; b < m; b++) {
          const int g = pj.sig2[(size_t)a * m + b];
          const double *xb = bx.data() + b * sm; double *yb = by.data() + b * sm;
          if (b == a) {
            for (int j = 0; j < m; j++) yb[j] = 0;
            yb[i] = clip(g >= 0 ? xb[i] - pj.Dsig[(size_t)g * m + i] : xb[i]);
          } else {
            if (g >= 0) { const double *cs = &pj.Csig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) yb[j] = clip(xb[j] - cs[j]); }
            else for (int j = 0; j < m; j++) yb[j] = clip(xb[j]);
            yb[i] = 0;
          }
        }
        for (int b = 0; b < m; b++) {
          const int g = pj.sig2[(size_t)a * m + b];
          if (g < 0) continue;
          const double *yb = by.data() + b * sm;
          if (b == a) pj.Gdsig[(size_t)g * m + i] += yb[i];
          else { double *gs = &pj.Gsig[g * mm + (size_t)i * m]; for (int j = 0; j < m; j++) gs[j] += yb[j]; }
        }
        yr[0] = clip(xr[0]);
        for (int j = 0; j < m; j++) { double *yj = yr + 1 + j * sm; for (int b = 0; b < m; b++) yj[b] = by[b * sm + j]; }
      }
    }
}
// pstep with a side-1 variant (no elim): pstep's preamble, then the rows
static void pstep_var(CutProjector &pj, const double *X, double *Y, int var)
{
  const int m = pj.m, n1 = pj.n1;
  const size_t mm = (size_t)m * m, ns = pj.nsig();
  std::fill(pj.Csig.begin(), pj.Csig.end(), 0.0); std::fill(pj.Dsig.begin(), pj.Dsig.end(), 0.0);
  for (size_t g = 0; g < ns; g++) {
    double *cs = &pj.Csig[g * mm], *ds = &pj.Dsig[g * m];
    for (int s : pj.sig_reg[g]) {
      const double *cw = &pj.CWs[s * mm], *cpv = &pj.CPs[(size_t)s * m];
      for (int i = 0; i < m; i++) {
        ds[i] += cpv[i];
        for (int j = i + 1; j < m; j++) cs[(size_t)i * m + j] += 0.5 * cw[(size_t)i * m + j];
      }
    }
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) cs[(size_t)j * m + i] = cs[(size_t)i * m + j];
  }
  std::fill(pj.Gsig.begin(), pj.Gsig.end(), 0.0); std::fill(pj.Gdsig.begin(), pj.Gdsig.end(), 0.0);
  auto clip = [](double v) { return std::min(1.0, std::max(0.0, v)); };
  Y[0] = 1.0;
  for (int c = 1; c < n1; c++) Y[c] = clip(X[c]);
  rows1_variant(pj, X, Y, NULL, var);
  pj.moments_from_groups();
}

int main(int argc, char **argv)
{
  const int n = 50, reps = argc > 1 ? atoi(argv[1]) : 20;
  // node 13 on wil50: facilities 21, 8, 28 fixed to locations 28, 41, 19; side 0 regions
  // over the free locations? the region side is the first index (facilities on wil50.dat)
  std::vector<int> free_orig;
  for (int a = 0; a < n; a++) if (a != 21 && a != 8 && a != 28) free_orig.push_back(a);
  const int m = (int)free_orig.size(), n1 = m * m + 1;
  std::vector<Region> regions = grid_regions(5, 10, free_orig, n);
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::vector<Cut> pool;
  for (int q = 0; q < 1000; q++) {
    int i = rng() % m, j = rng() % m, k = rng() % m;
    if (i == j || j == k || i == k) continue;
    int v[3] = {i, j, k}; std::sort(v, v + 3);
    Cut t = {(int)(rng() % 4), v[0], v[1], v[2], (int)(rng() % regions.size()), 0.0, 0, 0, 0, 0};
    pool.push_back(t);
  }
  const size_t nn = (size_t)n1 * n1;
  std::vector<double> X(nn), Y(nn);
  for (int r = 0; r < n1; r++) for (int c = r; c < n1; c++) X[(size_t)r * n1 + c] = X[(size_t)c * n1 + r] = (-0.2 + 1.4 * U(rng)) / m;
  std::vector<double> Yref;
  for (int arm = 0; arm < 7; arm++) {
    const int side = arm == 0 ? 0 : 1, var = arm <= 1 ? 0 : arm - 1;   // var 1 B, 2 T, 3 B+prefetch, 4 B+fused, 5 both
    CutProjector pj(m, regions, pool, side);
    std::mt19937 r2(11);
    for (size_t t = 0; t < pj.nu.size(); t++) pj.nu[t] = U(r2) < 0.5 ? 1e-3 * U(r2) : 0;
    pj.rebuild_sums();
    pj.pstep(X.data(), Y.data(), NULL);
    double best = 1e9, tot = 0;
    static std::vector<double> junk(64u << 20, 1.0);   // 512 MB: flush the caches before each call
    for (int r = 0; r < reps; r++) {
      double js = 0; for (size_t q = 0; q < junk.size(); q += 8) { junk[q] += 1e-9; js += junk[q]; }
      if (js == 0) printf("x");
      const double t0 = CutProjector::tcpu();
      if (var == 0) pj.pstep(X.data(), Y.data(), NULL); else pstep_var(pj, X.data(), Y.data(), var);
      const double dt = CutProjector::tcpu() - t0;
      best = std::min(best, dt); tot += dt;
    }
    double cs = 0; for (size_t t = 0; t < nn; t += 977) cs += Y[t];
    double ps = 0; for (double v : pj.Pag) ps += v;
    if (var == 0) {   // the dense term at the checks, on this side
      std::vector<double> L(nn);
      double bd = 1e9;
      for (int r = 0; r < 3; r++) {
        std::fill(L.begin(), L.end(), 0.0);
        double js = 0; for (size_t q = 0; q < junk.size(); q += 8) { junk[q] += 1e-9; js += junk[q]; }
        if (js == 0) printf("x");
        double ab = 0, hs = 0;
        const double t0 = CutProjector::tcpu();
        pj.dense_term(1.0L, pj.nu, L.data(), ab, hs);
        bd = std::min(bd, CutProjector::tcpu() - t0);
      }
      printf("  side %d dense_term best %.2f ms\n", side, 1e3 * bd);
    }
    if (arm == 1) { Yref = Y; }
    double dmax = 0; if (arm >= 2) for (size_t t = 0; t < nn; t++) dmax = std::max(dmax, fabs(Y[t] - Yref[t]));
    printf("arm %d (var %d, |Y - side-1 ref| %.1e) side %d: m %d regions %zu classes %zu: pstep best %.2f ms mean %.2f ms (checksum %.12g %.12g)\n", arm, var, dmax, side, m,
           regions.size(), pj.nsig(), 1e3 * best, 1e3 * tot / reps, cs, ps);
  }
  return 0;
}
