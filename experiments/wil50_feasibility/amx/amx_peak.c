// Raw AMX bf16 throughput on this core: tdpbf16ps with operands held in tile registers
// (pure TMUL rate), then with tile loads from an L1-resident buffer each step.
#define _GNU_SOURCE
#include <immintrin.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/syscall.h>
#define ARCH_REQ_XCOMP_PERM 0x1023
#define XFEATURE_XTILEDATA 18
struct __attribute__((packed, aligned(64))) tilecfg { uint8_t palette; uint8_t start_row; uint8_t res[14]; uint16_t colsb[16]; uint8_t rows[16]; };
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(){
  if (syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA)) { perror("arch_prctl"); return 1; }
  struct tilecfg c; memset(&c,0,sizeof c); c.palette=1;
  for (int i=0;i<8;i++){ c.colsb[i]=64; c.rows[i]=16; }
  _tile_loadconfig(&c);
  static uint16_t A[16*32] __attribute__((aligned(64))), B[16*32] __attribute__((aligned(64)));
  static float C[16*16] __attribute__((aligned(64)));
  for (int i=0;i<16*32;i++){ A[i]=0x3f80; B[i]=0x3f80; }   // 1.0 in bf16
  _tile_loadd(4, A, 64); _tile_loadd(5, A, 64); _tile_loadd(6, B, 64); _tile_loadd(7, B, 64);
  _tile_zero(0); _tile_zero(1); _tile_zero(2); _tile_zero(3);
  long n = 200000;
  double t0=now();
  for (long i=0;i<n;i++){ _tile_dpbf16ps(0,4,6); _tile_dpbf16ps(1,4,7); _tile_dpbf16ps(2,5,6); _tile_dpbf16ps(3,5,7); }
  double dt=now()-t0; double flops = (double)n*4*16*16*32*2;
  _tile_stored(0, C, 64);
  printf("TMUL only (operands in tiles): %.0f GF/s bf16  (C[0]=%.0f)\n", flops/dt/1e9, C[0]);
  // with loads: each step loads 2 A tiles and 2 B tiles from L1 (same addresses), 4 dot products
  t0=now();
  for (long i=0;i<n;i++){ _tile_loadd(4, A, 64); _tile_loadd(5, A, 64); _tile_loadd(6, B, 64); _tile_loadd(7, B, 64);
    _tile_dpbf16ps(0,4,6); _tile_dpbf16ps(1,4,7); _tile_dpbf16ps(2,5,6); _tile_dpbf16ps(3,5,7); }
  dt=now()-t0; printf("TMUL + 4 tile loads from L1 per 4 products: %.0f GF/s bf16\n", flops/dt/1e9);
  // the fused split pattern: per step load Wh, Wl, Yh, Yl (4 loads), 3 products into one C, x2 column blocks
  t0=now();
  for (long i=0;i<n;i++){ _tile_loadd(4, A, 64); _tile_loadd(5, A, 64); _tile_loadd(6, B, 64); _tile_loadd(7, B, 64);
    _tile_dpbf16ps(0,4,6); _tile_dpbf16ps(0,4,7); _tile_dpbf16ps(0,5,6);
    _tile_loadd(6, B+16, 64); _tile_loadd(7, B+16, 64);
    _tile_dpbf16ps(1,4,6); _tile_dpbf16ps(1,4,7); _tile_dpbf16ps(1,5,6); }
  dt=now()-t0; printf("fused split pattern (6 loads, 6 products): %.0f GF/s bf16 = %.0f GF/s fp32-equivalent\n", (double)n*6*16*16*32*2/dt/1e9, (double)n*2*16*16*32*2/dt/1e9);
  _tile_release();
  return 0; }
