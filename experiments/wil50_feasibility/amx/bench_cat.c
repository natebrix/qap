// Split bf16 as ONE GEMM with a tripled inner dimension: T = [Wh | Wh | Wl] (K x 3K) * [Yh; Yl; Yh] (3K x kb).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <dlfcn.h>
typedef void (*bgemm_t)(int,int,int,int,int,int,float,const uint16_t*,int,const uint16_t*,int,float,float*,int);
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static inline uint16_t bf16(float f){ uint32_t u; memcpy(&u,&f,4); u += 0x7FFF + ((u>>16)&1); return (uint16_t)(u>>16); }
static inline float bf16f(uint16_t h){ uint32_t u=(uint32_t)h<<16; float f; memcpy(&f,&u,4); return f; }
int main(int argc, char **argv){
  int K = argc>1? atoi(argv[1]) : 2117, kb = argc>2? atoi(argv[2]) : 600, reps = 5;
  void *h = dlopen("/usr/local/lib/libmkl_rt.so.3", RTLD_NOW); bgemm_t bgemm = (bgemm_t)dlsym(h,"cblas_gemm_bf16bf16f32");
  size_t nw=(size_t)K*K, ny=(size_t)K*kb;
  float *W=malloc(nw*4), *Y=malloc(ny*4), *T=malloc(ny*4); double *Wd=malloc(nw*8), *Yd=malloc(ny*8), *Td=calloc(ny,8);
  srand(1); for(size_t i=0;i<nw;i++){ W[i]=(float)rand()/RAND_MAX-0.5f; Wd[i]=W[i]; }
  for(size_t i=0;i<ny;i++){ Y[i]=(float)rand()/RAND_MAX-0.5f; Yd[i]=Y[i]; }
  for(int j=0;j<kb;j++) for(int k=0;k<K;k++){ double y=Yd[k+(size_t)j*K]; const double *w=Wd+(size_t)k*K; double *t=Td+(size_t)j*K; for(int i=0;i<K;i++) t[i]+=w[i]*y; }
  double tmax=0; for(size_t i=0;i<ny;i++) tmax=fmax(tmax,fabs(Td[i]));
  // A' column-major K x 3K: columns 0..K-1 = Wh, K..2K-1 = Wh, 2K..3K-1 = Wl (each K x K block contiguous since lda = K)
  uint16_t *A=malloc(nw*3*2), *B=malloc(ny*3*2);
  double t0=now();
  for(size_t i=0;i<nw;i++){ uint16_t hi=bf16(W[i]); A[i]=hi; A[nw+i]=hi; A[2*nw+i]=bf16(W[i]-bf16f(hi)); }
  double tsw=now()-t0;
  // B' column-major 3K x kb: column j = [Yh(:,j); Yl(:,j); Yh(:,j)]
  int K3=3*K;
  t0=now();
  for(int j=0;j<kb;j++){ const float *y=Y+(size_t)j*K; uint16_t *b=B+(size_t)j*K3;
    for(int i=0;i<K;i++){ uint16_t hi=bf16(y[i]); b[i]=hi; b[K+i]=bf16(y[i]-bf16f(hi)); b[2*K+i]=hi; } }
  double tsy=now()-t0;
  bgemm(102,111,111,K,kb,K3,1,A,K,B,K3,0,T,K);
  t0=now(); for(int r=0;r<reps;r++) bgemm(102,111,111,K,kb,K3,1,A,K,B,K3,0,T,K);
  double dt=(now()-t0)/reps;
  double e=0; for(size_t i=0;i<ny;i++) e=fmax(e,fabs(T[i]-Td[i]));
  printf("split W (x3 layout) %.1f ms, split Y (x3 layout) %.1f ms\n", tsw*1e3, tsy*1e3);
  printf("one GEMM K'=3K:  %.1f ms  (%.0f GF/s bf16 actual, %.0f GF/s-equivalent fp32)  max err rel %.1e\n", dt*1e3, 2.0*K*K3*kb/dt/1e9, 2.0*K*K*kb/dt/1e9, e/tmax);
  printf("per filter step incl. Y split: %.1f ms\n", (dt+tsy)*1e3);
  return 0; }
