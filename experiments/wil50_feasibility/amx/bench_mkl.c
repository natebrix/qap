// The filter GEMM T = W Y (W: K x K fp32 symmetric-like, Y: K x kb) through oneMKL,
// loaded at run time: fp32 sgemm, one-pass bf16 (gemm_bf16bf16f32), and split bf16
// (W = Wh + Wl, Y = Yh + Yl in bf16; T ~ Wh Yh + Wh Yl + Wl Yh, fp32 accumulate).
// Errors against an fp64 reference. One thread.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <dlfcn.h>
typedef void (*sgemm_t)(int,int,int,int,int,int,float,const float*,int,const float*,int,float,float*,int);
typedef void (*bgemm_t)(int,int,int,int,int,int,float,const uint16_t*,int,const uint16_t*,int,float,float*,int);
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static inline uint16_t bf16(float f){ uint32_t u; memcpy(&u,&f,4); u += 0x7FFF + ((u>>16)&1); return (uint16_t)(u>>16); }
static inline float bf16f(uint16_t h){ uint32_t u=(uint32_t)h<<16; float f; memcpy(&f,&u,4); return f; }
static void split(const float *a, size_t n, uint16_t *hi, uint16_t *lo){
  for (size_t i=0;i<n;i++){ hi[i]=bf16(a[i]); lo[i]=bf16(a[i]-bf16f(hi[i])); } }
int main(int argc, char **argv){
  int K = argc>1? atoi(argv[1]) : 2117, kb = argc>2? atoi(argv[2]) : 600, reps = 5;
  void *h = dlopen("/usr/local/lib/libmkl_rt.so.3", RTLD_NOW); if(!h){ printf("dlopen: %s\n", dlerror()); return 1; }
  sgemm_t sgemm = (sgemm_t)dlsym(h,"cblas_sgemm"); bgemm_t bgemm = (bgemm_t)dlsym(h,"cblas_gemm_bf16bf16f32");
  if(!sgemm||!bgemm){ printf("dlsym failed\n"); return 1; }
  const int ColMajor=102, NoTrans=111;
  size_t nw=(size_t)K*K, ny=(size_t)K*kb;
  float *W=malloc(nw*4), *Y=malloc(ny*4), *T=malloc(ny*4); double *Wd=malloc(nw*8), *Yd=malloc(ny*8), *Td=calloc(ny,8);
  srand(1); for(size_t i=0;i<nw;i++){ W[i]=(float)rand()/RAND_MAX-0.5f; Wd[i]=W[i]; }
  for(size_t i=0;i<ny;i++){ Y[i]=(float)rand()/RAND_MAX-0.5f; Yd[i]=Y[i]; }
  // fp64 reference (plain loops, column-major: T[i + j*K] = sum_k W[i + k*K] Y[k + j*K])
  for(int j=0;j<kb;j++) for(int k=0;k<K;k++){ double y=Yd[k+(size_t)j*K]; const double *w=Wd+(size_t)k*K; double *t=Td+(size_t)j*K; for(int i=0;i<K;i++) t[i]+=w[i]*y; }
  double tmax=0; for(size_t i=0;i<ny;i++) tmax=fmax(tmax,fabs(Td[i]));
  #define ERR(name) do{ double e=0; for(size_t i=0;i<ny;i++) e=fmax(e,fabs(T[i]-Td[i])); printf("  %-28s max err %.2e (relative to max|T| %.2e: %.1e)\n", name, e, tmax, e/tmax);}while(0)
  // fp32
  sgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,W,K,Y,K,0,T,K);
  double t0=now(); for(int r=0;r<reps;r++) sgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,W,K,Y,K,0,T,K);
  double dt=(now()-t0)/reps; printf("MKL sgemm fp32:           %7.1f ms  %6.0f GF/s\n", dt*1e3, 2.0*K*K*kb/dt/1e9); ERR("fp32");
  // bf16 one pass
  uint16_t *Wh=malloc(nw*2),*Wl=malloc(nw*2),*Yh=malloc(ny*2),*Yl=malloc(ny*2);
  t0=now(); split(W,nw,Wh,Wl); double tsw=now()-t0; t0=now(); split(Y,ny,Yh,Yl); double tsy=now()-t0;
  printf("split W %.1f ms, split Y %.1f ms\n", tsw*1e3, tsy*1e3);
  bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wh,K,Yh,K,0,T,K);
  t0=now(); for(int r=0;r<reps;r++) bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wh,K,Yh,K,0,T,K);
  dt=(now()-t0)/reps; printf("MKL bf16 one pass:        %7.1f ms  %6.0f GF/s\n", dt*1e3, 2.0*K*K*kb/dt/1e9); ERR("bf16 x1");
  // split bf16, three passes (Y split each time, W split once per iteration in the solver)
  t0=now(); for(int r=0;r<reps;r++){ split(Y,ny,Yh,Yl);
    bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wh,K,Yh,K,0,T,K);
    bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wh,K,Yl,K,1,T,K);
    bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wl,K,Yh,K,1,T,K); }
  dt=(now()-t0)/reps; printf("MKL split bf16 x3 (+Y split): %4.1f ms  %6.0f GF/s-equivalent\n", dt*1e3, 2.0*K*K*kb/dt/1e9); ERR("split bf16 x3");
  // two passes (drop Wl Yh): W at 8 bits, Y at 16
  t0=now(); for(int r=0;r<reps;r++){ split(Y,ny,Yh,Yl);
    bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wh,K,Yh,K,0,T,K);
    bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wh,K,Yl,K,1,T,K); }
  dt=(now()-t0)/reps; printf("MKL split bf16 x2 (W hi only): %3.1f ms  %6.0f GF/s-equivalent\n", dt*1e3, 2.0*K*K*kb/dt/1e9); ERR("split bf16 x2");
  return 0; }
