#ifndef AMX_GEMM_H
#define AMX_GEMM_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { int K, Kp, kbp; uint16_t *Wh, *Wl, *Yh, *Yl; float *Cp; } amx_ctx;
int amx_available(void);                         // AMX-BF16 present and permitted for this process
int amx_ctx_init(amx_ctx *c, int K, int kbmax);  // buffers for W (K x K) and Y (K x <= kbmax)
void amx_ctx_free(amx_ctx *c);
void amx_pack_W(amx_ctx *c, const float *W);     // once per W (symmetric, column-major fp32)
void amx_gemm(amx_ctx *c, const float *Y, int kb, float *T);   // T = W Y, column-major fp32
extern double amx_split_ms, amx_kernel_ms;        // timing of the last amx_gemm call
#ifdef __cplusplus
}
#endif
#endif
