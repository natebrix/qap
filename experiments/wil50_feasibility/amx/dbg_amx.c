#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../../../latest/amx_gemm.h"
static inline uint16_t bf16r(float f){ uint32_t u; memcpy(&u,&f,4); u += 0x7FFF + ((u>>16)&1); return (uint16_t)(u>>16); }
static inline float bf16f(uint16_t h){ uint32_t u=(uint32_t)h<<16; float f; memcpy(&f,&u,4); return f; }
int main(int argc, char **argv){
  int K = argc>1? atoi(argv[1]) : 96, kb = argc>2? atoi(argv[2]) : 40;
  if (!amx_available()) { printf("no AMX\n"); return 1; }
  size_t nw=(size_t)K*K, ny=(size_t)K*kb;
  float *W=malloc(nw*4), *Y=malloc(ny*4), *T=calloc(ny,4); double *Td=calloc(ny,8);
  srand(1);
  for (int i=0;i<K;i++) for (int j=0;j<=i;j++){ float v=(float)rand()/RAND_MAX-0.5f; W[i+(size_t)j*K]=v; W[j+(size_t)i*K]=v; }
  for(size_t i=0;i<ny;i++) Y[i]=(float)rand()/RAND_MAX-0.5f;
  for(int j=0;j<kb;j++) for(int k=0;k<K;k++){ double y=Y[k+(size_t)j*K]; const float *w=W+(size_t)k*K; double *t=Td+(size_t)j*K; for(int i=0;i<K;i++) t[i]+=(double)w[i]*y; }
  amx_ctx c; amx_ctx_init(&c, K, kb);
  amx_pack_W(&c, W);
  // check the pack against the definition
  int bad=0; const int Kp=c.Kp, nI=Kp/16;
  for (int k32=0;k32<Kp/32;k32++) for (int i16=0;i16<nI;i16++) for (int r=0;r<16;r++) for (int p=0;p<2;p++) for (int n=0;n<16;n++){
    int k=k32*32+2*r+p, i=i16*16+n; uint16_t h=c.Wh[((size_t)k32*nI+i16)*512+r*32+2*n+p], l=c.Wl[((size_t)k32*nI+i16)*512+r*32+2*n+p];
    uint16_t eh = (k<K&&i<K)? bf16r(W[k+(size_t)i*K]) : 0; uint16_t el = (k<K&&i<K)? bf16r(W[k+(size_t)i*K]-bf16f(eh)) : 0;
    if (h!=eh||l!=el){ if(bad<5) printf("pack mismatch k=%d i=%d: got %04x/%04x want %04x/%04x\n",k,i,h,l,eh,el); bad++; } }
  printf("pack mismatches: %d\n", bad);
  amx_gemm(&c, Y, kb, T);
  // check the Y split
  bad=0; for (int j=0;j<kb;j++) for (int k=0;k<K;k++){ uint16_t h=c.Yh[k+(size_t)j*Kp], l=c.Yl[k+(size_t)j*Kp]; uint16_t eh=bf16r(Y[k+(size_t)j*K]), el=bf16r(Y[k+(size_t)j*K]-bf16f(eh)); if(h!=eh||l!=el){ if(bad<5) printf("Y split mismatch k=%d j=%d: %04x/%04x want %04x/%04x\n",k,j,h,l,eh,el); bad++; } }
  printf("Y split mismatches: %d\n", bad);
  double tmax=0; for(size_t i=0;i<ny;i++) tmax=fmax(tmax,fabs(Td[i]));
  printf("per-tile max relative error (rows = i tiles, cols = j tiles):\n");
  for (int i=0;i<K;i+=16){ for (int j=0;j<kb;j+=16){ double e=0; for(int a=i;a<i+16&&a<K;a++) for(int b=j;b<j+16&&b<kb;b++) e=fmax(e,fabs(T[a+(size_t)b*K]-Td[a+(size_t)b*K])); printf(" %8.1e", e/tmax);} printf("\n"); }
  return 0; }
