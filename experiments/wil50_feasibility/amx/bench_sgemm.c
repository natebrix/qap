// Baseline: the solver's filter GEMM, W (K x K) times Y (K x kb), fp32, one thread,
// through the same scipy_openblas32 the solver links (symbols prefixed scipy_).
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
extern void scipy_sgemm_(const char*, const char*, const int*, const int*, const int*, const float*,
                         const float*, const int*, const float*, const int*, const float*, float*, const int*);
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(int argc, char **argv){
  int K = argc>1? atoi(argv[1]) : 2117, kb = argc>2? atoi(argv[2]) : 600, reps = 5;
  float *W = malloc((size_t)K*K*4), *Y = malloc((size_t)K*kb*4), *T = malloc((size_t)K*kb*4);
  srand(1); for (size_t i=0;i<(size_t)K*K;i++) W[i] = (float)rand()/RAND_MAX-0.5f;
  for (size_t i=0;i<(size_t)K*kb;i++) Y[i] = (float)rand()/RAND_MAX-0.5f;
  float one=1, zero=0;
  scipy_sgemm_("N","N",&K,&kb,&K,&one,W,&K,Y,&K,&zero,T,&K);
  double t0=now(); for(int r=0;r<reps;r++) scipy_sgemm_("N","N",&K,&kb,&K,&one,W,&K,Y,&K,&zero,T,&K);
  double dt=(now()-t0)/reps; printf("OpenBLAS sgemm %dx%dx%d: %.1f ms, %.1f GF/s\n", K,K,kb, dt*1e3, 2.0*K*K*kb/dt/1e9);
  return 0; }
