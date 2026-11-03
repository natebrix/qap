#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../../../latest/amx_gemm.h"
extern void scipy_sgemm_(const char*, const char*, const int*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(int argc, char **argv){
  int K = argc>1? atoi(argv[1]) : 2117, kb = argc>2? atoi(argv[2]) : 600, reps = 5;
  if (!amx_available()) { printf("no AMX\n"); return 1; }
  size_t nw=(size_t)K*K, ny=(size_t)K*kb;
  float *W=malloc(nw*4), *Y=malloc(ny*4), *T=malloc(ny*4), *T32=malloc(ny*4); double *Td=calloc(ny,8);
  srand(1);
  for (int i=0;i<K;i++) for (int j=0;j<=i;j++){ float v=(float)rand()/RAND_MAX-0.5f; W[i+(size_t)j*K]=v; W[j+(size_t)i*K]=v; }
  for(size_t i=0;i<ny;i++) Y[i]=(float)rand()/RAND_MAX-0.5f;
  for(int j=0;j<kb;j++) for(int k=0;k<K;k++){ double y=Y[k+(size_t)j*K]; const float *w=W+(size_t)k*K; double *t=Td+(size_t)j*K; for(int i=0;i<K;i++) t[i]+=(double)w[i]*y; }
  double tmax=0; for(size_t i=0;i<ny;i++) tmax=fmax(tmax,fabs(Td[i]));
  amx_ctx c; if (amx_ctx_init(&c, K, kb)) return 1;
  double t0=now(); amx_pack_W(&c, W); double tp=now()-t0;
  amx_gemm(&c, Y, kb, T);
  t0=now(); for(int r=0;r<reps;r++) amx_gemm(&c, Y, kb, T); double dt=(now()-t0)/reps;
  double e=0; for(size_t i=0;i<ny;i++) e=fmax(e,fabs(T[i]-Td[i]));
  printf("AMX fused split: pack W %.1f ms (once per iteration); gemm %.1f ms (Y split %.1f, kernel %.1f) = %.0f GF/s-equivalent; max err rel %.1e\n", tp*1e3, dt*1e3, amx_split_ms, amx_kernel_ms, 2.0*K*K*kb/dt/1e9, e/tmax);
  float one=1, zero=0; scipy_sgemm_("N","N",&K,&kb,&K,&one,W,&K,Y,&K,&zero,T32,&K);
  t0=now(); for(int r=0;r<reps;r++) scipy_sgemm_("N","N",&K,&kb,&K,&one,W,&K,Y,&K,&zero,T32,&K); double d32=(now()-t0)/reps;
  e=0; for(size_t i=0;i<ny;i++) e=fmax(e,fabs(T32[i]-Td[i]));
  printf("OpenBLAS sgemm fp32:            %.1f ms = %.0f GF/s; max err rel %.1e;  speedup %.2fx\n", d32*1e3, 2.0*K*K*kb/d32/1e9, e/tmax, d32/dt);
  amx_ctx_free(&c); return 0; }
