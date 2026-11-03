# The filter products on Intel AMX (Nov 3, cloud)

Item 1 of `../roadmap/node_speed.md`: the Chebyshev filter and Rayleigh-Ritz
products `T = W Y` (W the 2117 x 2117 reduced-face matrix of the iterate, Y a
block of 400-600 vectors) are five fp32 GEMMs per iteration, about 33 GF at
about 140 GF/s on this Xeon: 120-320 ms of a 410-530 ms iteration. The Xeon
has AMX-BF16 (`amx_bf16 amx_tile` in `/proc/cpuinfo`, kernel 6.18) that
OpenBLAS cannot use because it only runs in bf16.

## Microbenchmarks (`bench_*.c`, `amx_*.c`; one thread, 2117 x 2117 x 600)

| path | ms | rate | error vs fp64 (relative to max entry) |
|---|---|---|---|
| OpenBLAS sgemm fp32 (`bench_sgemm.c`, the solver's library) | 36.8-40.9 | 132-146 GF/s | 5e-7 |
| oneMKL sgemm fp32 (`bench_mkl.c`, pip `mkl`, AVX-512 only) | 40-43 | 126-134 GF/s | 5e-7 |
| oneMKL bf16 one pass, AMX (`cblas_gemm_bf16bf16f32`) | 7.4-7.8 | 685-730 GF/s | **2e-3** (the level that drifted 7 units in the memo) |
| oneMKL bf16 three passes, W = Wh + Wl, Y = Yh + Yl, drop Wl Yl | 23.6 | 228 GF/s-equivalent | 4e-6 |
| the same with W packed once (`bench_pack.c`) or as one GEMM with tripled inner dimension (`bench_cat.c`) | 23.6-23.7 | 227 | 4e-6 |
| raw AMX, operands in tile registers (`amx_peak.c`) | | **2,100-2,280 GF/s bf16** | |
| raw AMX, the fused split pattern (6 tile loads, 6 products, L1) | | 1,840-1,870 GF/s = 613-622 GF/s-equivalent | |
| our fused kernel (`latest/amx_gemm.c`, `test_amx_gemm.c`) | 17.4-19.8 | 272-310 GF/s-equivalent | 3.4e-6 |

MKL's bf16 kernel sits at a third of this core's AMX peak whatever the call
pattern, so three passes through it gain only 1.6x over OpenBLAS on the GEMM.
A fused kernel (the three products share tile loads and one accumulator)
runs at 2.0-2.3x; its pattern can run at the full peak when both operands
are in L1 or L2 (`amx_loads.c`: W from a 1 MB buffer, Y from L1, 2,100
GF/s in a good state), so the remaining gap is cache behaviour inside the
kernel, open for later tuning. Measurements on this VM swing 2-3x between
identical runs (`amx_peak.c` 1,089 against 2,128 GF/s back to back): AMX
frequency licensing or a neighbour; pin and repeat before trusting a number.

Two-pass splits (one operand at 8 bits) land at 1.5e-3, the one-pass level,
so three passes are the minimum; AMX-FP16 (one pass at 11 bits) is Granite
Rapids only, not on this chip.

## The kernel (`latest/amx_gemm.c`)

W is symmetric, so `T = W Y` is computed as the row-major product `C^T = Y^T W`:
the A tiles are the bf16 split of Y written straight into contiguous 1 KB
tiles, the B tiles are W packed once per iteration into the VNNI pair layout
(hi and lo parts, 18 MB at K = 2117, 5-7 ms), and the C tiles are stored
directly into T's column-major storage (edge tiles through a scratch block).
Tiles: two accumulators (two j-tiles of one i-tile), Yh/Yl for both j-tiles,
Wh/Wl for the i-tile; per 32-k step six tile loads and six `tdpbf16ps`:
`C += Yh Wh + Yh Wl + Yl Wh`. Blocking: an i-slab of W (Kp x 128, 1.1 MB)
in L2 across the j blocks, k-blocks of 8 so the two j-tiles' Y data (32 KB)
sit in L1 while the slab's i-tiles sweep over it, C partial sums reloaded
per k-block. Correctness: `dbg_amx.c` checks the pack and the split element
by element and prints the error per 16 x 16 tile; the first version failed
on padded row tiles (a negative-length copy) and AddressSanitizer found it.

## In the solver

`latest/dnn.C`: `DNN_AMX=1` (environment; off by default) routes the five
products of `proj_psd_subspace` through `amx_gemm` when `amx_available()`;
W is packed at the first product after each new W. Everything else (CholQR2,
the small Rayleigh-Ritz products, ssyevd, forming R in fp64, the certificate)
is unchanged, so the certified bound's validity does not depend on the
kernel: only the ADMM step is inexact, as with the fp32 filter already.
Build: `amx_gemm.o` is plain C compiled with the AMX flags on x86-64 Linux
(`latest/makefile`); other architectures get stubs and the sgemm path.

| check | sgemm | AMX | note |
|---|---|---|---|
| nug12 root, 60 iterations, certified lb | 563.27895 | 563.27899 | 0.54 / 0.70 ms per iteration (tiny K, overhead) |
| nug30 root, 100 iterations, certified lb | 5857.9953 | 5857.9948 | 52.97 / 43.12 ms per iteration: 1.23x |
| wil50 node 10, 300 iterations | pending | pending | |
| wil50 depth-4 node 3 (fathoms uncut at 510 iterations) | pending | pending | |
