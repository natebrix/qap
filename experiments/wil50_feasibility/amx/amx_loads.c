// Cost of tile loads by source: the fused pattern (6 loads, 6 products) with the operands
// (a) in L1, contiguous; (b) in L1, A tiles strided 4288 B; (c) streaming from a 1 MB buffer
// (L2); (d) streaming from a 32 MB buffer (L3/DRAM).
#define _GNU_SOURCE
#include <immintrin.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/syscall.h>
struct __attribute__((packed, aligned(64))) tilecfg { uint8_t palette; uint8_t start_row; uint8_t res[14]; uint16_t colsb[16]; uint8_t rows[16]; };
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static void run(const char *name, const uint16_t *A, long sa, size_t amask, const uint16_t *B, size_t bmask, long n){
  double t0=now();
  for (long i=0;i<n;i++){
    size_t ao=(i*2048)&amask, bo=(i*2048)&bmask;
    _tile_loadd(6, B+bo, 64); _tile_loadd(7, B+bo+512, 64);
    _tile_loadd(2, A+ao, sa); _tile_loadd(3, A+ao+32, sa);
    _tile_dpbf16ps(0,2,6); _tile_dpbf16ps(0,2,7); _tile_dpbf16ps(0,3,6);
    _tile_loadd(4, A+ao+64, sa); _tile_loadd(5, A+ao+96, sa);
    _tile_dpbf16ps(1,4,6); _tile_dpbf16ps(1,4,7); _tile_dpbf16ps(1,5,6);
  }
  double dt=now()-t0; printf("%-44s %5.0f GF/s bf16 (%4.0f fp32-equivalent)\n", name, (double)n*6*16*16*32*2/dt/1e9, (double)n*2*16*16*32*2/dt/1e9);
}
int main(){
  syscall(SYS_arch_prctl, 0x1023, 18);
  struct tilecfg c; memset(&c,0,sizeof c); c.palette=1; for (int i=0;i<8;i++){ c.colsb[i]=64; c.rows[i]=16; }
  _tile_loadconfig(&c);
  size_t big = 32u<<20;
  uint16_t *A = aligned_alloc(4096, big), *B = aligned_alloc(4096, big);
  memset(A, 0x3f, big); memset(B, 0x3f, big);
  _tile_zero(0); _tile_zero(1);
  long n = 100000;
  run("L1, contiguous tiles (stride 64)", A, 64, 0, B, 0, n);
  run("L1, A tiles strided 4288 B (16 pages)", A, 4288, 0, B, 0, n);
  run("L2 stream, 1 MB buffers, contiguous", A, 64, (1u<<20)/2-1 & ~1023u, B, (1u<<20)/2-1 & ~1023u, n);
  run("L2 stream, A strided 4288 B", A, 4288, (1u<<20)/2-1 & ~1023u, B, (1u<<20)/2-1 & ~1023u, n);
  run("L3/DRAM stream, 32 MB buffers, contiguous", A, 64, big/2-1 & ~1023u, B, big/2-1 & ~1023u, n);
  run("B (W) from 32 MB, A (Y) from L1", A, 64, 0, B, big/2-1 & ~1023u, n);
  run("B (W) from 1 MB (L2), A (Y) from L1", A, 64, 0, B, (1u<<20)/2-1 & ~1023u, n);
  _tile_release(); return 0; }
