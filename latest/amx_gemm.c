// Fused split-bf16 GEMM on Intel AMX for the solver's filter product T = W Y,
// W symmetric K x K (fp32, column-major), Y K x kb (fp32, column-major), T fp32.
// W = Wh + Wl and Y = Yh + Yl in bf16 (round to nearest even); T ~ Wh Yh + Wh Yl + Wl Yh
// with fp32 tile accumulation: about 16 significant bits per operand, against 8 for
// one-pass bf16 and 24 for fp32. Computed as the row-major product C^T = Y^T W (W is
// symmetric): A tiles are read straight from the bf16 split of Y (16 columns of Y x
// 32 k, stride Kp), B tiles from W packed once per W into the VNNI pair layout, and
// C tiles land directly in T's column-major storage (edge tiles through a scratch
// block). Blocking: an i-slab of W (Kp x IB, hi and lo) stays in L2 across the j
// blocks; inside, a k-block of the two j-tiles' Y data (KB x 32 x 4 bytes) stays in
// L1 while the i-tiles of the slab sweep over it, with C partial sums in a small
// local block, so the L2 stream is the W tiles alone (2 KB per 6 tile products).
#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/mman.h>
#include "amx_gemm.h"
#include <time.h>
#ifndef AMX_IB
#define AMX_IB 128   // i-slab width: W tiles Kp x IB x 4 bytes in L2
#endif
#ifndef AMX_KB
#define AMX_KB 8     // k-block: KB x 32 k per C partial sum
#endif
#if defined(__x86_64__) && defined(__linux__) && defined(__AMX_BF16__) && defined(__AVX512BF16__)
#include <immintrin.h>
double amx_split_ms = 0, amx_kernel_ms = 0;   // of the last amx_gemm call
static double now_ms(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec*1e3+1e-6*t.tv_nsec; }

struct __attribute__((packed, aligned(64))) tilecfg { uint8_t palette; uint8_t start_row; uint8_t res[14]; uint16_t colsb[16]; uint8_t rows[16]; };

int amx_available(void)
{
  static int ok = -1;
  if (ok < 0) {
    unsigned a, b, c, d;
    __asm__ volatile("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) : "a"(7), "c"(0));
    int tile = (d >> 24) & 1, bf = (d >> 22) & 1, vbf = (c >> 5) & 1 /* AVX512_BF16 is leaf 7 subleaf 1 eax bit 5 */;
    unsigned a1, b1, c1, d1;
    __asm__ volatile("cpuid" : "=a"(a1), "=b"(b1), "=c"(c1), "=d"(d1) : "a"(7), "c"(1));
    vbf = (a1 >> 5) & 1;
    ok = (tile && bf && vbf && syscall(SYS_arch_prctl, 0x1023 /*ARCH_REQ_XCOMP_PERM*/, 18 /*XTILEDATA*/) == 0) ? 1 : 0;
  }
  return ok;
}

static void *al64(size_t n)
{
  void *p = NULL;
  size_t a = n >= (2u << 20) ? (2u << 20) : 64;          // 2 MB alignment for the large buffers
  if (posix_memalign(&p, a, n ? n : 64)) return NULL;
#ifdef MADV_HUGEPAGE
  if (a > 64) madvise(p, n, MADV_HUGEPAGE);
#endif
  memset(p, 0, n); return p;
}

int amx_ctx_init(amx_ctx *c, int K, int kbmax)
{
  memset(c, 0, sizeof *c);
  c->K = K; c->Kp = (K + 31) / 32 * 32; c->kbp = (kbmax + 31) / 32 * 32;
  size_t nw = (size_t)c->Kp * c->Kp, ny = (size_t)c->Kp * c->kbp;
  c->Wh = al64(nw * 2); c->Wl = al64(nw * 2); c->Yh = al64(ny * 2); c->Yl = al64(ny * 2);
  c->Ah = al64(ny * 2); c->Al = al64(ny * 2); c->Bh = al64(ny * 2); c->Bl = al64(ny * 2);
  c->Gh = al64((size_t)c->kbp * c->kbp * 2); c->Gl = al64((size_t)c->kbp * c->kbp * 2);
  size_t scr = (size_t)(AMX_IB / 16) * 512 * 4;             // edge-tile scratch for one i-slab
  if ((size_t)4 * c->Kp * 2 > scr) scr = (size_t)4 * c->Kp * 2;  // and four bf16 rows for the pack
  c->Cp = al64(scr);
  return (c->Wh && c->Wl && c->Yh && c->Yl && c->Ah && c->Al && c->Bh && c->Bl && c->Gh && c->Gl && c->Cp) ? 0 : -1;
}
void amx_ctx_free(amx_ctx *c)
{
  free(c->Wh); free(c->Wl); free(c->Yh); free(c->Yl); free(c->Ah); free(c->Al); free(c->Bh); free(c->Bl); free(c->Gh); free(c->Gl); free(c->Cp);
  memset(c, 0, sizeof *c);
}

// 16 floats -> bf16 hi (round to nearest even) and lo = bf16(x - hi)
static inline void split16(const float *x, uint16_t *h, uint16_t *l)
{
  __m512 v = _mm512_loadu_ps(x);
  __m256bh hb = _mm512_cvtneps_pbh(v);
  __m256i hi = (__m256i)hb;
  __m512 hf = _mm512_castsi512_ps(_mm512_slli_epi32(_mm512_cvtepu16_epi32(hi), 16));
  __m256bh lb = _mm512_cvtneps_pbh(_mm512_sub_ps(v, hf));
  _mm256_storeu_si256((__m256i *)h, hi);
  _mm256_storeu_si256((__m256i *)l, (__m256i)lb);
}
static inline void split_n(const float *x, uint16_t *h, uint16_t *l, int n)
{
  int k = 0;
  for (; k + 16 <= n; k += 16) split16(x + k, h + k, l + k);
  if (k < n) {   // tail: scalar
    for (; k < n; k++) {
      uint32_t u; float f = x[k]; memcpy(&u, &f, 4); u += 0x7FFF + ((u >> 16) & 1); uint16_t hi = (uint16_t)(u >> 16);
      uint32_t uh = (uint32_t)hi << 16; float fh; memcpy(&fh, &uh, 4);
      float r = f - fh; memcpy(&u, &r, 4); u += 0x7FFF + ((u >> 16) & 1);
      h[k] = hi; l[k] = (uint16_t)(u >> 16);
    }
  }
}

// W (symmetric, column-major, ld K) -> VNNI tiles of the bf16 hi and lo parts.
// Tile (k32, i16) at ((i16 * Kp/32) + k32) * 512, so an i-tile's k sweep reads a contiguous
// stream (1.2x over the (k32, i16) order); row r, pair n: [2n+p] = W[k0+2r+p][i0+n].
// Row k of W is contiguous (W[k][i] = W[i][k] = W + k*K): split it with the vector
// routine, then interleave the pairs.
void amx_pack_W(amx_ctx *c, const float *W)
{
  const int K = c->K, Kp = c->Kp, nI = Kp / 16;
  uint16_t *rh0 = (uint16_t *)c->Cp, *rl0 = rh0 + Kp, *rh1 = rl0 + Kp, *rl1 = rh1 + Kp;   // 4 rows of Kp bf16 fit the scratch (8 KB each at Kp 2144)
  for (int k2 = 0; k2 < Kp / 2; k2++) {
    const int k0 = 2 * k2, k1 = k0 + 1;
    if (k0 < K) { split_n(W + (size_t)k0 * K, rh0, rl0, K); } else { memset(rh0, 0, K * 2); memset(rl0, 0, K * 2); }
    if (k1 < K) { split_n(W + (size_t)k1 * K, rh1, rl1, K); } else { memset(rh1, 0, K * 2); memset(rl1, 0, K * 2); }
    for (int i = K; i < Kp; i++) { rh0[i] = rl0[i] = rh1[i] = rl1[i] = 0; }
    const int k32 = k0 / 32, r = (k0 % 32) / 2;
    for (int i16 = 0; i16 < nI; i16++) {
      uint16_t *th = c->Wh + ((size_t)i16 * (Kp / 32) + k32) * 512 + r * 32, *tl = c->Wl + ((size_t)i16 * (Kp / 32) + k32) * 512 + r * 32;
      // interleave 16 pairs: (rh0[i], rh1[i]) for i in the i16 block
      __m256i a = _mm256_loadu_si256((const __m256i *)(rh0 + i16 * 16)), b = _mm256_loadu_si256((const __m256i *)(rh1 + i16 * 16));
      __m512i w = _mm512_inserti64x4(_mm512_castsi256_si512(a), b, 1);
      // permute 16-bit lanes: out[2n] = a[n], out[2n+1] = b[n]
      static const uint16_t idx[32] = {0,16,1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31};
      __m512i perm = _mm512_loadu_si512(idx);
      _mm512_storeu_si512(th, _mm512_permutexvar_epi16(perm, w));
      a = _mm256_loadu_si256((const __m256i *)(rl0 + i16 * 16)); b = _mm256_loadu_si256((const __m256i *)(rl1 + i16 * 16));
      w = _mm512_inserti64x4(_mm512_castsi256_si512(a), b, 1);
      _mm512_storeu_si512(tl, _mm512_permutexvar_epi16(perm, w));
    }
  }
}

// Columns of X (n rows, ncols columns, column-major, ld ldx) as A tiles: tile (c16, k32) at
// (c16 * nK32 + k32) * 512 holds 16 columns x 32 rows in bf16 (hi in Th, lo in Tl), rows
// n..np-1 and columns ncols..(ncols rounded to 32)-1 zero.
static void split_cols(const float *X, int ldx, int n, int ncols, uint16_t *Th, uint16_t *Tl)
{
  const int np = (n + 31) / 32 * 32, nK32 = np / 32, ncp = (ncols + 31) / 32 * 32;
  for (int j = 0; j < ncp; j++) {
    const int j16 = j / 16, r = j % 16;
    uint16_t *th = Th + (size_t)j16 * nK32 * 512 + r * 32, *tl = Tl + (size_t)j16 * nK32 * 512 + r * 32;
    if (j >= ncols) { for (int k32 = 0; k32 < nK32; k32++) { memset(th + (size_t)k32 * 512, 0, 64); memset(tl + (size_t)k32 * 512, 0, 64); } continue; }
    const float *x = X + (size_t)j * ldx;
    int k32 = 0;
    for (; (k32 + 1) * 32 <= n; k32++) {
      split16(x + k32 * 32, th + (size_t)k32 * 512, tl + (size_t)k32 * 512);
      split16(x + k32 * 32 + 16, th + (size_t)k32 * 512 + 16, tl + (size_t)k32 * 512 + 16);
    }
    if (k32 < nK32) {
      uint16_t hb[32], lb[32]; memset(hb, 0, 64); memset(lb, 0, 64);
      split_n(x + k32 * 32, hb, lb, n - k32 * 32);
      memcpy(th + (size_t)k32 * 512, hb, 64); memcpy(tl + (size_t)k32 * 512, lb, 64);
    }
  }
}

// Y (n rows x ncols, column-major, ld ldy) as VNNI B tiles with the inner dimension along
// the ROWS: tile (k32, c16) at (c16 * nK32 + k32) * 512, row r, pair n: [2n+p] = Y[k0+2r+p][c0+n]
// (gathered: the pair is adjacent in memory, the 16 columns are ld apart).
static void pack_vnni_rows(const float *Y, int ldy, int n, int ncols, uint16_t *Bh, uint16_t *Bl)
{
  const int np = (n + 31) / 32 * 32, nK32 = np / 32, ncp = (ncols + 31) / 32 * 32, nC16 = ncp / 16;
  static const uint16_t idx[32] = {0,16,1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31};
  const __m512i perm = _mm512_loadu_si512(idx);
  __m512i cidx = _mm512_mullo_epi32(_mm512_set_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0), _mm512_set1_epi32(ldy));
  for (int c16 = 0; c16 < nC16; c16++)
    for (int k32 = 0; k32 < nK32; k32++) {
      uint16_t *th = Bh + ((size_t)c16 * nK32 + k32) * 512, *tl = Bl + ((size_t)c16 * nK32 + k32) * 512;
      const int ncol = c16 * 16 + 16 <= ncols ? 16 : (ncols - c16 * 16 > 0 ? ncols - c16 * 16 : 0);
      const __mmask16 cm = (__mmask16)((1u << ncol) - 1);
      for (int r = 0; r < 16; r++) {
        const int k0 = k32 * 32 + 2 * r;
        __m512 a = _mm512_setzero_ps(), b = _mm512_setzero_ps();
        if (k0 < n) a = _mm512_mask_i32gather_ps(a, cm, cidx, Y + (size_t)c16 * 16 * ldy + k0, 4);
        if (k0 + 1 < n) b = _mm512_mask_i32gather_ps(b, cm, cidx, Y + (size_t)c16 * 16 * ldy + k0 + 1, 4);
        __m256i ah = (__m256i)_mm512_cvtneps_pbh(a), bh = (__m256i)_mm512_cvtneps_pbh(b);
        __m512 af = _mm512_castsi512_ps(_mm512_slli_epi32(_mm512_cvtepu16_epi32(ah), 16));
        __m512 bf = _mm512_castsi512_ps(_mm512_slli_epi32(_mm512_cvtepu16_epi32(bh), 16));
        __m256i al = (__m256i)_mm512_cvtneps_pbh(_mm512_sub_ps(a, af)), bl = (__m256i)_mm512_cvtneps_pbh(_mm512_sub_ps(b, bf));
        _mm512_storeu_si512(th + r * 32, _mm512_permutexvar_epi16(perm, _mm512_inserti64x4(_mm512_castsi256_si512(ah), bh, 1)));
        _mm512_storeu_si512(tl + r * 32, _mm512_permutexvar_epi16(perm, _mm512_inserti64x4(_mm512_castsi256_si512(al), bl, 1)));
      }
    }
}

// Y (n rows x ncols, column-major, ld ldy) as VNNI B tiles with the inner dimension along
// the COLUMNS (B = Y^T): tile (j32, k16) at (k16 * nJ32 + j32) * 512, row r, pair n:
// [2n+p] = Y[k0+n][j0+2r+p] (two columns, 16 contiguous rows each, interleaved).
static void pack_vnni_cols(const float *Y, int ldy, int n, int ncols, uint16_t *Bh, uint16_t *Bl)
{
  const int np = (n + 31) / 32 * 32, nK16 = np / 16, ncp = (ncols + 31) / 32 * 32, nJ32 = ncp / 32;
  static const uint16_t idx[32] = {0,16,1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31};
  const __m512i perm = _mm512_loadu_si512(idx);
  for (int k16 = 0; k16 < nK16; k16++) {
    const int nrow = k16 * 16 + 16 <= n ? 16 : (n - k16 * 16 > 0 ? n - k16 * 16 : 0);
    const __mmask16 rm = (__mmask16)((1u << nrow) - 1);
    for (int j32 = 0; j32 < nJ32; j32++) {
      uint16_t *th = Bh + ((size_t)k16 * nJ32 + j32) * 512, *tl = Bl + ((size_t)k16 * nJ32 + j32) * 512;
      for (int r = 0; r < 16; r++) {
        const int j0 = j32 * 32 + 2 * r;
        __m512 a = _mm512_setzero_ps(), b = _mm512_setzero_ps();
        if (j0 < ncols) a = _mm512_maskz_loadu_ps(rm, Y + (size_t)j0 * ldy + k16 * 16);
        if (j0 + 1 < ncols) b = _mm512_maskz_loadu_ps(rm, Y + (size_t)(j0 + 1) * ldy + k16 * 16);
        __m256i ah = (__m256i)_mm512_cvtneps_pbh(a), bh = (__m256i)_mm512_cvtneps_pbh(b);
        __m512 af = _mm512_castsi512_ps(_mm512_slli_epi32(_mm512_cvtepu16_epi32(ah), 16));
        __m512 bf = _mm512_castsi512_ps(_mm512_slli_epi32(_mm512_cvtepu16_epi32(bh), 16));
        __m256i al = (__m256i)_mm512_cvtneps_pbh(_mm512_sub_ps(a, af)), bl = (__m256i)_mm512_cvtneps_pbh(_mm512_sub_ps(b, bf));
        _mm512_storeu_si512(th + r * 32, _mm512_permutexvar_epi16(perm, _mm512_inserti64x4(_mm512_castsi256_si512(ah), bh, 1)));
        _mm512_storeu_si512(tl + r * 32, _mm512_permutexvar_epi16(perm, _mm512_inserti64x4(_mm512_castsi256_si512(al), bl, 1)));
      }
    }
  }
}

// The fused three-product tile loop: C (rows = A tiles ia, columns = B tiles jb) over nK32
// k-steps; A tiles (ia, k32) at (ia * nK32 + k32) * 512 in Ah/Al, B tiles (jb, k32) at
// (jb * nK32 + k32) * 512 in Bh/Bl. Stores C tile (ia, jb) through store(ia, jb, tile).
#define AMX_TILECFG() do { struct tilecfg cfg; memset(&cfg, 0, sizeof cfg); cfg.palette = 1; \
  for (int t = 0; t < 8; t++) { cfg.colsb[t] = 64; cfg.rows[t] = 16; } _tile_loadconfig(&cfg); } while (0)

// M (ka x kb, row-major, ld kb) = X^T Y for X (n x ka) and Y (n x kb), both column-major
// with leading dimension n. X == NULL: X's tiles are the ones of the last amx_gemm's Y.
void amx_gemm_tn(amx_ctx *c, const float *X, const float *Y, int ka, int kb, float *M)
{
  const int n = c->K, np = c->Kp, nK32 = np / 32;
  const int kap = (ka + 31) / 32 * 32, kbp = (kb + 31) / 32 * 32;
  const uint16_t *Ah = c->Yh, *Al = c->Yl;
  if (X) { split_cols(X, n, n, ka, c->Ah, c->Al); Ah = c->Ah; Al = c->Al; }
  pack_vnni_rows(Y, n, n, kb, c->Bh, c->Bl);
  AMX_TILECFG();
  float *cs = c->Cp;
  for (int ia = 0; ia < kap / 16; ia++) {
    const uint16_t *ah = Ah + (size_t)ia * nK32 * 512, *al = Al + (size_t)ia * nK32 * 512;
    for (int jb = 0; jb < kbp / 16; jb += 2) {
      const uint16_t *bh1 = c->Bh + (size_t)jb * nK32 * 512, *bl1 = c->Bl + (size_t)jb * nK32 * 512;
      const uint16_t *bh2 = bh1 + (size_t)nK32 * 512, *bl2 = bl1 + (size_t)nK32 * 512;
      _tile_zero(0); _tile_zero(1);
      for (int k32 = 0; k32 < nK32; k32++) {
        _tile_loadd(6, ah + (size_t)k32 * 512, 64);
        _tile_loadd(7, al + (size_t)k32 * 512, 64);
        _tile_loadd(2, bh1 + (size_t)k32 * 512, 64);
        _tile_loadd(3, bl1 + (size_t)k32 * 512, 64);
        _tile_dpbf16ps(0, 6, 2);
        _tile_dpbf16ps(0, 6, 3);
        _tile_dpbf16ps(0, 7, 2);
        _tile_loadd(4, bh2 + (size_t)k32 * 512, 64);
        _tile_loadd(5, bl2 + (size_t)k32 * 512, 64);
        _tile_dpbf16ps(1, 6, 4);
        _tile_dpbf16ps(1, 6, 5);
        _tile_dpbf16ps(1, 7, 4);
      }
      // C tile rows = 16 a-columns (ia), columns = 16 b-columns (jb): M[ia*16 + r][jb*16 + n]
      const int ra = ia * 16, cb1 = jb * 16, cb2 = cb1 + 16;
      const int in1 = ra + 16 <= ka && cb1 + 16 <= kb, in2 = ra + 16 <= ka && cb2 + 16 <= kb;
      if (in1) _tile_stored(0, M + (size_t)ra * kb + cb1, (long)kb * 4); else _tile_stored(0, cs, 64);
      if (in2) _tile_stored(1, M + (size_t)ra * kb + cb2, (long)kb * 4); else _tile_stored(1, cs + 256, 64);
      const int nr = ra + 16 <= ka ? 16 : ka - ra;
      if (!in1) for (int r = 0; r < nr; r++) for (int q = 0; q < 16 && cb1 + q < kb; q++) M[(size_t)(ra + r) * kb + cb1 + q] = cs[r * 16 + q];
      if (!in2) for (int r = 0; r < nr; r++) for (int q = 0; q < 16 && cb2 + q < kb; q++) M[(size_t)(ra + r) * kb + cb2 + q] = cs[256 + r * 16 + q];
    }
  }
  _tile_release();
}

// Xn (n x kout, column-major, ld n) = Y G for Y (n x kin, column-major, ld n) and G
// (kin x kout, column-major, ld ldg). Computed as Xn^T = G^T Y^T: A tiles from the
// columns of G, B tiles from Y^T (pairs along Y's columns). Xn may alias Y.
void amx_gemm_ng(amx_ctx *c, const float *Y, int kin, const float *G, int ldg, int kout, float *Xn)
{
  const int n = c->K, np = c->Kp, nK16 = np / 16;
  const int kinp = (kin + 31) / 32 * 32, nJ32 = kinp / 32, koutp = (kout + 31) / 32 * 32;
  split_cols(G, ldg, kin, kout, c->Gh, c->Gl);          // tile (o16, j32): 16 output columns x 32 inner j
  pack_vnni_cols(Y, n, n, kin, c->Bh, c->Bl);           // tile (j32, k16) at (k16 * nJ32 + j32)
  AMX_TILECFG();
  float *cs = c->Cp;
  for (int o16 = 0; o16 < koutp / 16; o16++) {
    const uint16_t *gh = c->Gh + (size_t)o16 * nJ32 * 512, *gl = c->Gl + (size_t)o16 * nJ32 * 512;
    for (int k16 = 0; k16 < nK16; k16 += 2) {
      const uint16_t *bh1 = c->Bh + (size_t)k16 * nJ32 * 512, *bl1 = c->Bl + (size_t)k16 * nJ32 * 512;
      const uint16_t *bh2 = bh1 + (size_t)nJ32 * 512, *bl2 = bl1 + (size_t)nJ32 * 512;
      const int two = k16 + 1 < nK16;
      _tile_zero(0); _tile_zero(1);
      for (int j32 = 0; j32 < nJ32; j32++) {
        _tile_loadd(6, gh + (size_t)j32 * 512, 64);
        _tile_loadd(7, gl + (size_t)j32 * 512, 64);
        _tile_loadd(2, bh1 + (size_t)j32 * 512, 64);
        _tile_loadd(3, bl1 + (size_t)j32 * 512, 64);
        _tile_dpbf16ps(0, 6, 2);
        _tile_dpbf16ps(0, 6, 3);
        _tile_dpbf16ps(0, 7, 2);
        if (two) {
          _tile_loadd(4, bh2 + (size_t)j32 * 512, 64);
          _tile_loadd(5, bl2 + (size_t)j32 * 512, 64);
          _tile_dpbf16ps(1, 6, 4);
          _tile_dpbf16ps(1, 6, 5);
          _tile_dpbf16ps(1, 7, 4);
        }
      }
      // C tile rows = 16 output columns (o16), columns = 16 rows k of Xn: Xn[k16*16 + q][o16*16 + r]
      const int co = o16 * 16, k1 = k16 * 16, k2 = k1 + 16;
      const int in1 = co + 16 <= kout && k1 + 16 <= n, in2 = two && co + 16 <= kout && k2 + 16 <= n;
      if (in1) _tile_stored(0, Xn + (size_t)co * n + k1, (long)n * 4); else _tile_stored(0, cs, 64);
      if (two) { if (in2) _tile_stored(1, Xn + (size_t)co * n + k2, (long)n * 4); else _tile_stored(1, cs + 256, 64); }
      if (!in1) for (int r = 0; r < 16 && co + r < kout; r++) for (int q = 0; q < 16 && k1 + q < n; q++) Xn[(size_t)(co + r) * n + k1 + q] = cs[r * 16 + q];
      if (two && !in2) for (int r = 0; r < 16 && co + r < kout; r++) for (int q = 0; q < 16 && k2 + q < n; q++) Xn[(size_t)(co + r) * n + k2 + q] = cs[256 + r * 16 + q];
    }
  }
  _tile_release();
}

// T (K x kb, column-major, ld K) = W Y for Y (K x kb, column-major, ld K)
void amx_gemm(amx_ctx *c, const float *Y, int kb, float *T)
{
  const int K = c->K, Kp = c->Kp, nI = Kp / 16, nK32 = Kp / 32;
  const int kbp = (kb + 31) / 32 * 32;
  // split Y straight into contiguous A tiles: tile (j16, k32) at (j16 * nK32 + k32) * 512,
  // row r = column j16*16 + r of Y, 32 bf16 for k = k32*32 .. +31 (zero beyond K and kb)
  double t0 = now_ms();
  for (int j = 0; j < kbp; j++) {
    const int j16 = j / 16, r = j % 16;
    uint16_t *th = c->Yh + (size_t)j16 * nK32 * 512 + r * 32, *tl = c->Yl + (size_t)j16 * nK32 * 512 + r * 32;
    if (j >= kb) { for (int k32 = 0; k32 < nK32; k32++) { memset(th + (size_t)k32 * 512, 0, 64); memset(tl + (size_t)k32 * 512, 0, 64); } continue; }
    const float *y = Y + (size_t)j * K;
    int k32 = 0;
    for (; (k32 + 1) * 32 <= K; k32++) {
      split16(y + k32 * 32, th + (size_t)k32 * 512, tl + (size_t)k32 * 512);
      split16(y + k32 * 32 + 16, th + (size_t)k32 * 512 + 16, tl + (size_t)k32 * 512 + 16);
    }
    if (k32 < nK32) {   // the last, partial k tile
      uint16_t hb[32], lb[32]; memset(hb, 0, 64); memset(lb, 0, 64);
      split_n(y + k32 * 32, hb, lb, K - k32 * 32);
      memcpy(th + (size_t)k32 * 512, hb, 64); memcpy(tl + (size_t)k32 * 512, lb, 64);
    }
  }
  amx_split_ms = now_ms() - t0; t0 = now_ms();
  struct tilecfg cfg; memset(&cfg, 0, sizeof cfg); cfg.palette = 1;
  for (int t = 0; t < 8; t++) { cfg.colsb[t] = 64; cfg.rows[t] = 16; }
  _tile_loadconfig(&cfg);
  const long sT = (long)K * 4;
  float *cs = c->Cp;   // scratch C block: 32 j x IB i, column-major with ld IB... we use per-tile 16x16 blocks
  const int IB = AMX_IB;          // i-slab: W tiles Kp x IB x 2 x 2 bytes (1.1 MB at IB 128) in L2
  const int KB = AMX_KB;          // k-block of KB x 32: the two j-tiles' Y data KB x 32 x 4 x 32 bytes (32 KB at KB 8) in L1
  for (int i0 = 0; i0 < Kp; i0 += IB) {
    const int i1 = i0 + IB < Kp ? i0 + IB : Kp;
    for (int j0 = 0; j0 < kbp; j0 += 32) {
      const uint16_t *yh1 = c->Yh + (size_t)(j0 / 16) * nK32 * 512, *yl1 = c->Yl + (size_t)(j0 / 16) * nK32 * 512;
      const uint16_t *yh2 = yh1 + (size_t)nK32 * 512, *yl2 = yl1 + (size_t)nK32 * 512;
      // C block for (j0..j0+31) x (i0..i1): tiles live in T when fully inside, else in scratch
      // partial sums over k-blocks: first k-block zeroes, later ones load the tile back
      for (int kb0 = 0; kb0 < nK32; kb0 += KB) {
        const int kb1 = kb0 + KB < nK32 ? kb0 + KB : nK32;
        for (int i = i0; i < i1; i += 16) {
          const int i16 = i / 16;
          const int inside1 = (i + 16 <= K) && (j0 + 16 <= kb), inside2 = (i + 16 <= K) && (j0 + 32 <= kb);
          float *c0 = inside1 ? T + (size_t)j0 * K + i : cs + (size_t)((i - i0) / 16) * 512;
          float *c1 = inside2 ? T + (size_t)(j0 + 16) * K + i : cs + (size_t)((i - i0) / 16) * 512 + 256;
          const long s0 = inside1 ? sT : 64, s1 = inside2 ? sT : 64;
          if (kb0 == 0) { _tile_zero(0); _tile_zero(1); }
          else { _tile_loadd(0, c0, s0); _tile_loadd(1, c1, s1); }
          for (int k32 = kb0; k32 < kb1; k32++) {
            const size_t toff = ((size_t)i16 * nK32 + k32) * 512;
            _tile_stream_loadd(6, c->Wh + toff, 64);   // W streams: the non-temporal hint keeps Y's k-block in L1
            _tile_stream_loadd(7, c->Wl + toff, 64);
            _tile_loadd(2, yh1 + (size_t)k32 * 512, 64);
            _tile_loadd(3, yl1 + (size_t)k32 * 512, 64);
            _tile_dpbf16ps(0, 2, 6);
            _tile_dpbf16ps(0, 2, 7);
            _tile_dpbf16ps(0, 3, 6);
            _tile_loadd(4, yh2 + (size_t)k32 * 512, 64);
            _tile_loadd(5, yl2 + (size_t)k32 * 512, 64);
            _tile_dpbf16ps(1, 4, 6);
            _tile_dpbf16ps(1, 4, 7);
            _tile_dpbf16ps(1, 5, 6);
          }
          _tile_stored(0, c0, s0);
          _tile_stored(1, c1, s1);
        }
      }
      // copy the edge tiles out of the scratch
      for (int i = i0; i < i1 && i < K; i += 16) {   // padded row tiles (i >= K) are discarded
        const int inside1 = (i + 16 <= K) && (j0 + 16 <= kb), inside2 = (i + 16 <= K) && (j0 + 32 <= kb);
        const float *b0 = cs + (size_t)((i - i0) / 16) * 512, *b1 = b0 + 256;
        const int ni = i + 16 <= K ? 16 : K - i;
        if (!inside1) for (int r = 0; r < 16 && j0 + r < kb; r++) memcpy(T + (size_t)(j0 + r) * K + i, b0 + r * 16, ni * 4);
        if (!inside2) for (int r = 0; r < 16 && j0 + 16 + r < kb; r++) memcpy(T + (size_t)(j0 + 16 + r) * K + i, b1 + r * 16, ni * 4);
      }
    }
  }
  _tile_release();
  amx_kernel_ms = now_ms() - t0;
}

#else   /* no AMX-BF16 compile support on this target: stubs, the solver falls back to sgemm */
double amx_split_ms = 0, amx_kernel_ms = 0;
int amx_available(void) { return 0; }
int amx_ctx_init(amx_ctx *c, int K, int kbmax) { (void)K; (void)kbmax; memset(c, 0, sizeof *c); return -1; }
void amx_ctx_free(amx_ctx *c) { memset(c, 0, sizeof *c); }
void amx_pack_W(amx_ctx *c, const float *W) { (void)c; (void)W; }
void amx_gemm(amx_ctx *c, const float *Y, int kb, float *T) { (void)c; (void)Y; (void)kb; (void)T; }
void amx_gemm_tn(amx_ctx *c, const float *X, const float *Y, int ka, int kb, float *M) { (void)c; (void)X; (void)Y; (void)ka; (void)kb; (void)M; }
void amx_gemm_ng(amx_ctx *c, const float *Y, int kin, const float *G, int ldg, int kout, float *Xn) { (void)c; (void)Y; (void)kin; (void)G; (void)ldg; (void)kout; (void)Xn; }
#endif
