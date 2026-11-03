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
| wil50 node 10 (m = 47), 300 iterations, certified lb | 48,582.776 | 48,582.745 | **521 / 443 ms per iteration: 1.18x**; filter+RR 300 -> 217 ms per step at kb 605 (`logs/val_n10_*`) |
| wil50 depth-4 node 3 (m = 46), fathoms uncut | fathomed at iteration 510 | fathomed at iteration 510 | 48,814.559 / 48,814.563; **478 / 397 ms per iteration: 1.20x** (`logs/val_d4n3_*`) |
| `tests/dnn_test.py`, `make continuationcheck` with `DNN_AMX=1` | | pass | |

The four wil50 runs shared the four cores (one process per core), so the
per-iteration times are the loaded figures the search sees, not the idle
409 ms of the single-process profile.

## Where the iteration stands with AMX, and what is left

Node 10 with AMX, 300 iterations: eigensolver 95 s of 133 (71%): filter+RR
62 s (217 ms per step: the five products about 100 ms including the 6 ms
pack, the rest CholQR2, the two small Rayleigh-Ritz GEMMs, the 605 x 605
ssyevd and the element-wise Chebyshev updates), full ssyevd refresh 11 s
(37 ms per iteration), forming R 16 s (53 ms: dsyrk in fp64 plus the
symmetrisation copy); forward 12 s, expand 9 s, Y/Z-step 9 s, certify 4 s.
The AMX item alone is therefore about 1.2x, not the memo's 1.4-1.7x; the
rest of that range needs (a) the kernel at its pattern's measured ceiling
(17 -> about 10 ms per product, another 35 ms per iteration), (b) the
small BLAS-3 of CholQR2 and Rayleigh-Ritz on AMX too (about 25 ms), and (c)
the CPU trims of the memo (R in fp32, the Weyl copy, the refresh cadence,
about 60-80 ms). All three together are worth about 1.45-1.5x on the
iteration over today's sgemm path.

## Trims on top of the kernel (Nov 3, afternoon)

Step-level timers (`dnntest -v -v`, `step split` and `form R` lines) priced
every piece of a subspace step; the trims below were then measured on node 10
(300 iterations) and on the depth-4 node 3 (fathoms uncut; degree 4 fathoms
it at iteration 510), four processes on the four cores (`logs/trim_*`), and
the survivors on the accel study's 32-node harness
(`../solvers/accel/run_all32.sh`, `summarize.py d4 d3L2 amxA amxB`).

| trim | switch | node 10, ms per iteration | depth-4 node 3 | verdict |
|---|---|---|---|---|
| AMX for the five W*Y products | `DNN_AMX=1` | 521 -> 443 | fathoms at 510, as fp32 | in |
| degree 3 with tight lambda (the accel study's recommendation, never adopted by the search; with it the Weyl pass is skipped) | `-d 3 -L 2` (search: `QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2`) | 443 -> 409 | 520 | in |
| one-pass CholQR | `DNN_CHOLQR_PASSES=1` | 409 -> 390 (CholQR 44 -> 22 ms) | 520, no factorisation failure, same 4 block grows | in |
| block guard K/40 instead of K/20 | `DNN_GUARD_DIV=40` | 395 -> 371 | 520, 11 grows instead of 4, bound at 300 lower by 0.3 | out (3% for 0.8 of bound on the 32 nodes' misses) |
| R's rank term in fp32 | `DNN_R32=1` | 409 -> 395 | 520 | out: the certificate's rigorous margin grows from 5e-6 to 0.08 per certificate (W - R is no longer exactly negative semidefinite), a systematic bound loss for 3% |
| the small products of CholQR and Rayleigh-Ritz on AMX | `DNN_AMX=2` | 378 -> 387 | 520 | out: at inner dimension 550 the tile loop is store-bound and the Gram matrix's VNNI pack is a gather; CholQR 22 -> 31 ms, the others break even; code kept, default off |
| full eigendecomposition every 40 iterations (with tight lambda; the accel study found 120 too coarse at m = 47) | `-r 40` (search: `QAP_DNN_REFRESH=40`) | | | in: 32-node run below, -8% CPU, two nodes fathom earlier, none later |
| the transforms: the Y/Z step writes the next forward's row phase (row fusion) | `DNN_FUSE=1` | 375 -> 382: forward -9 ms, Y/Z step +11 (a sixth row stream in an already memory-bound loop); bound bit-identical | 520 | out |
| skip W's symmetrisation (symmetric to rounding by construction; the eigensolver reads one triangle or multiplies by the full matrix) and write Wf inside forward's gather | `DNN_NOSYMW=1` | 375 -> 356: forward -6, convert -4 | 520 | in (bundle E below); bound 48,582.59 against 48,582.61 |

**The 32-node harness** (`runs/amxA_*`, `runs/amxB_*`; d4 = degree 4 fp32, the
search's setting; d3L2 = degree 3 tight lambda fp32):

| configuration | fathomed | ms per iteration | median iterations | total iterations | CPU, fathoms | CPU, misses | misses' lb at 200 vs d4 |
|---|---|---|---|---|---|---|---|
| d4 | 18 | 562 | 220 | 4,240 | 2,181 s | 2,889 s | 0 |
| d3L2 | 18 | 480 | 220 | 4,250 | 1,889 s | 2,330 s | +0.81 |
| **amxA: AMX, degree 3 tight, one-pass CholQR** | **18** | **354** | **220** | **4,250** | **1,394 s (1.56x)** | **1,720 s (1.68x)** | **+0.88** |
| amxB: amxA with guard K/40 | 18 | 341 | 220 | 4,270 | 1,352 s | 1,643 s | +0.11 (min -0.37) |
| **amxD: amxA with the refresh every 40** | **18** | **330** | **220** | **4,220** | **1,277 s (1.71x)** | **1,653 s (1.75x)** | **+0.03 (min -0.41)** |

Bundle A fathoms the same 18 nodes as degree 4 at the same iteration on 15 of
them and 10 iterations off on 3 (one check interval), with the misses' bound
ahead. Bundle D adds the refresh every 40: the same 18 fathoms, two of them
10-20 iterations earlier and none later, the misses' bound at 200 back at
the degree-4 level, 8% less CPU again. **D is the production setting from
here**: `DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2
QAP_DNN_REFRESH=40` in the search (the first two are read by `dnn.C`, the
rest by `bnb.C`): **1.71x less CPU on the fathoms and 1.75x on the misses
than the search's current setting**, with the same decisions.

Where the 354 ms go (node 10, bundle A, loaded): the five products 81 + pack 7,
CholQR 22, Rayleigh-Ritz 13.5 + eigensolve 27.5 + update 13, convert and
Chebyshev updates 8 (eigensolver step about 172); the full ssyevd refresh 35;
forming R 55 (fp64 dsyrk 35, symmetrise 13, F 7); the memory-bound transforms
and the Y/Z step 98; the certificate 14. The remaining levers, in order: the
transforms (fuse the Y/Z step with the next forward's row pass, skip the
redundant symmetrisation of W, fuse gather, convert and pack: about -25 ms),
the refresh cadence (-17 if `-r 40` holds on the 32 nodes), the kernel's last
third (the pattern's measured ceiling is 2 TF against about 1 TF achieved),
and R in fp32 if a 0.08 margin is acceptable (-14).
