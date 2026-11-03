// Packed split bf16: W's hi and lo parts packed once (A operand), then per GEMM only Y is converted.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <dlfcn.h>
typedef size_t (*getsz_t)(int,int,int,int);
typedef void (*pack_t)(int,int,int,int,int,int,const uint16_t*,int,uint16_t*);
typedef void (*compute_t)(int,int,int,int,int,int,float,const uint16_t*,int,const uint16_t*,int,float,float*,int);
typedef void (*bgemm_t)(int,int,int,int,int,int,float,const uint16_t*,int,const uint16_t*,int,float,float*,int);
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static inline uint16_t bf16(float f){ uint32_t u; memcpy(&u,&f,4); u += 0x7FFF + ((u>>16)&1); return (uint16_t)(u>>16); }
static inline float bf16f(uint16_t h){ uint32_t u=(uint32_t)h<<16; float f; memcpy(&f,&u,4); return f; }
int main(int argc, char **argv){
  int K = argc>1? atoi(argv[1]) : 2117, kb = argc>2? atoi(argv[2]) : 600, reps = 5;
  void *h = dlopen("/usr/local/lib/libmkl_rt.so.3", RTLD_NOW);
  getsz_t getsz=(getsz_t)dlsym(h,"cblas_gemm_bf16bf16f32_pack_get_size"); pack_t pack=(pack_t)dlsym(h,"cblas_gemm_bf16bf16f32_pack");
  compute_t compute=(compute_t)dlsym(h,"cblas_gemm_bf16bf16f32_compute"); bgemm_t bgemm=(bgemm_t)dlsym(h,"cblas_gemm_bf16bf16f32");
  if(!getsz||!pack||!compute){ printf("dlsym failed\n"); return 1; }
  const int ColMajor=102, NoTrans=111, Packed=151, AMatrix=161;
  size_t nw=(size_t)K*K, ny=(size_t)K*kb;
  float *W=malloc(nw*4), *Y=malloc(ny*4), *T=malloc(ny*4), *T0=malloc(ny*4);
  srand(1); for(size_t i=0;i<nw;i++) W[i]=(float)rand()/RAND_MAX-0.5f;
  for(size_t i=0;i<ny;i++) Y[i]=(float)rand()/RAND_MAX-0.5f;
  uint16_t *Wh=malloc(nw*2),*Wl=malloc(nw*2),*Yh=malloc(ny*2),*Yl=malloc(ny*2);
  for(size_t i=0;i<nw;i++){ Wh[i]=bf16(W[i]); Wl[i]=bf16(W[i]-bf16f(Wh[i])); }
  for(size_t i=0;i<ny;i++){ Yh[i]=bf16(Y[i]); Yl[i]=bf16(Y[i]-bf16f(Yh[i])); }
  // reference: the unpacked 3-pass result
  bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wh,K,Yh,K,0,T0,K); bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wh,K,Yl,K,1,T0,K); bgemm(ColMajor,NoTrans,NoTrans,K,kb,K,1,Wl,K,Yh,K,1,T0,K);
  size_t psz = getsz(AMatrix,K,kb,K); printf("packed A size %.1f MB (bf16 K x K = %.1f MB)\n", psz/1e6, nw*2/1e6);
  uint16_t *Ph=malloc(psz), *Pl=malloc(psz);
  double t0=now(); pack(ColMajor,AMatrix,NoTrans,K,kb,K,Wh,K,Ph); pack(ColMajor,AMatrix,NoTrans,K,kb,K,Wl,K,Pl); double tp=now()-t0;
  printf("pack Wh+Wl once: %.1f ms\n", tp*1e3);
  compute(ColMajor,Packed,NoTrans,K,kb,K,1,Ph,K,Yh,K,0,T,K);
  t0=now(); for(int r=0;r<reps;r++) compute(ColMajor,Packed,NoTrans,K,kb,K,1,Ph,K,Yh,K,0,T,K);
  double d1=(now()-t0)/reps; printf("one pass, packed A:  %.1f ms  %.0f GF/s\n", d1*1e3, 2.0*K*K*kb/d1/1e9);
  t0=now(); for(int r=0;r<reps;r++){ for(size_t i=0;i<ny;i++){ Yh[i]=bf16(Y[i]); Yl[i]=bf16(Y[i]-bf16f(Yh[i])); }
    compute(ColMajor,Packed,NoTrans,K,kb,K,1,Ph,K,Yh,K,0,T,K);
    compute(ColMajor,Packed,NoTrans,K,kb,K,1,Ph,K,Yl,K,1,T,K);
    compute(ColMajor,Packed,NoTrans,K,kb,K,1,Pl,K,Yh,K,1,T,K); }
  double d3=(now()-t0)/reps; double e=0; for(size_t i=0;i<ny;i++) e=fmax(e,fabs(T[i]-T0[i]));
  printf("three passes, packed A, + Y split: %.1f ms  %.0f GF/s-equivalent  (diff vs unpacked %.1e)\n", d3*1e3, 2.0*K*K*kb/d3/1e9, e);
  return 0; }
