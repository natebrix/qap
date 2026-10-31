# Per-iteration cost of the DNN ADMM/PRSM at m = 47: where it goes and what would cut it

Research note, Oct 31 2026. Measurements on this cloud machine (Intel Xeon 2.1 GHz,
Sapphire Rapids core: AVX-512, AMX-bf16; OpenBLAS 0.3.34 via scipy_openblas32,
SkylakeX kernels, one thread), wil50 depth-3 node `34:40,28:22,14:29`, m = 47,
K = 2117, n1 = 2210, `dnntest -T 1 -d 4 -g 1.0 -G 1.0 -H 0 -i 200 -u 48816 -v -v`.
Logs are in the session scratchpad; the experimental tracker changes used below are
in `research_iteration_tracker.patch` (against `latest/dnn.C`, env-switched, not
meant to be applied as is).

## 1. Where the 530 ms go (degree 4, refresh 20, 200 iterations)

| phase | ms/iter | share | arithmetic / traffic | note |
|---|---|---|---|---|
| filter + Rayleigh-Ritz | 299 | 56% | 5 sgemm K x K x kb (kb ~ 564): 5 x 5.1 GF = 25 GF; CholQR2 ~4 GF; ssyevd(kb) + Y G ~3 GF | one sgemm K x K x 564 measured at 49 ms = 103 GF/s, ~75% of the core's fp32 peak (square sgemm reaches 146 GF/s). Compute-bound. |
| form R (dsyrk rank nneg ~460, double) | 53 | 10% | 2.1 GF fp64 + K^2 memcpy + triangle copy | ~40 GF/s fp64 |
| ssyevd refresh (K = 2117, single) | 37 | 7% | 734 ms per refresh, every 20 | |
| Weyl difference, Wprev copy, float conversion | 24 | 5% | 3 passes over K^2 doubles | memory-bound |
| forward V'(Y+Z/beta)V | 40 | 8% | ~50 MFLOP, ~250 MB traffic (T, rows_D, gather, symmetrize) | 6 GB/s: memory/latency-bound, not flop-bound |
| expand V R V' | 32 | 6% | same | |
| Y/Z-step | 30 | 6% | 5 arrays of n1^2 doubles, ~200 MB | memory-bound |
| certificate (every 10) | 14 | 3% | 143 ms each: forward 40 + dpotrf K (3.2 GF) ~90 + GL LAPs ~10 | |

Total arithmetic is only ~45 GFLOP per iteration, 34 of it fp32 GEMM. The O(m^4)
Householder transforms are already negligible in flops; they and the Y/Z-step are
pure memory streams (~130 ms, 25%). Apple M2 Pro's 3x comes from its AMX units
running sgemm at ~1.4 TF/s per cluster, which is exactly the compute-bound part.

## 2. Eigen-step findings (measured)

**The Weyl bound on lambda_max is 8x too loose and it costs iterations.** The
Chebyshev filter damps [theta_top, lam_hi], with lam_hi grown by ||W_k - W_{k-1}||_F
each iteration. At the first refresh lam_hi = 9.9 against a true lambda_max of 1.24;
at the second, 12.4 against 1.46. A filter interval 8x too wide weakens the
polynomial badly (T_4 amplification of an eigenvalue at -0.01 is 1.10 instead of 1.9,
at -0.18 it is 1.7 instead of 10), so early in a run, when W moves fast, the tracker
finds a fraction of the negative eigenvalues (22 of 536 at iteration 20). The ADMM
self-corrects, but late. Against an exact run (`-f`, full ssyevd each iteration):

| variant | lb at it 50 | it 100 | it 150 | it 200 | pres at 200 | ms/iter |
|---|---|---|---|---|---|---|
| exact `-f` | 37119.3 | 47171.3 | 48162.7 | 48287.6 | 3.7e-2 | 987 |
| default: degree 4, Weyl lam_hi | 36752.8 | 46969.2 | 48122.7 | 48279.4 | 5.7e-2 | 530 |
| degree 3, Weyl | 36601.1 | 46731.6 | 48009.1 | 48265.4 | 9.4e-2 | 474 |
| degree 2, tight lam_hi | 36747.1 | 46965.3 | 48072.4 | 48266.9 | 7.3e-2 | 450 |
| degree 3, tight lam_hi | 36925.9 | 47092.1 | 48134.3 | 48280.6 | 4.7e-2 | 520 |
| degree 4, tight lam_hi | 37017.3 | 47134.7 | 48151.0 | 48284.2 | 4.1e-2 | 558 |
| degree 3, tight, filter operands rounded to 10-bit mantissa (TF32/fp16 emulation) | 36925.7 | 47092.7 | 48133.8 | 48280.5 | 4.7e-2 | 522 |
| degree 4, Weyl, filter operands rounded to bf16 (7-bit mantissa) | 36753.2 | 46970.1 | 48119.3 | 48272.3 | 6.2e-2 | 539 |

"Tight lam_hi" = 3 power steps per iteration on the persistent top eigenvector of
W + sI (s = 1.05 |smallest Ritz value|, needed because |lambda_min| > lambda_max after
~60 iterations; without the shift the estimate collapses), lam_hi = 1.1 x estimate;
cost 3 sgemv ~ 6 ms. Conclusions:

- Degree 3 with a tight bound beats the current degree 4 at every checkpoint for 2%
  less time; degree 4 tight is within 3.4 of the exact trajectory at iteration 200
  (default: 8.2) and within 37 at iteration 100 (default: 202). The ms/iter does not
  drop by the 1/5 a dropped sgemm suggests because the block finds more negative
  eigenvalues (kb 598 vs 564). The gain is in iterations, not ms: at the it-100 level
  the exact run is ~10 iterations ahead of the default; for nodes that fathom at
  150-300 iterations that is a 5-10% saving, free. Needs the 800/2500-iteration
  protocol of experiments/eigensolver/README.md to confirm the floor is unchanged.
- TF32/fp16 mantissa is safe for the filter GEMMs (indistinguishable from fp32);
  bf16 drifts (7 below at 200 iterations, pres 6.2e-2 vs 5.7e-2), the same kind of
  floor as the degree-2 filter. This decides the hardware options below: tensor-core
  TF32/fp16 yes, Intel AMX (bf16 only) doubtful.
- The certificate is untouched by any of this: R = W - (computed negative part) keeps
  W - R negative semidefinite by construction; the margins printed (4.5e-6) did not
  change in any run.

## 3. Ideas, with the arithmetic at m = 47 (530 ms baseline)

1. **Tight lambda_max + degree 3** (above). Saving 2% per iteration, 5-10% of the
   iterations early. Effort: 40 lines (done in the patch). Confidence high on the
   measured 200 iterations, medium on the converged floor.
2. **Eigen-step trims.** Form R with ssyrk instead of dsyrk (the single-precision
   error ~1e-6 ||F||^2 is absorbed by the certificate's Cholesky shift, which already
   runs at 1e-7..1e-13 relative): -28 ms. Fuse the Weyl difference, Wprev copy and
   float conversion into forward's gather: -20 ms. Guard K/40 instead of K/20 (kb 564
   -> 512): -25 ms. One CholQR pass with a shifted Cholesky instead of two: -15 ms.
   Refresh every 40 once lam_hi is tight: -18 ms. Together ~105 ms, 1.25x. Effort
   1-2 days. Confidence high. (A sixth sgemm cannot be saved: the Rayleigh-Ritz needs
   W times the orthonormalized block, d+1 products are the floor for degree d.)
3. **Memory-lean transforms and Y/Z-step** (130 ms at 6 GB/s): a fused single pass
   over the n1^2 arrays and float storage for VRV and T halves the traffic: -50 ms,
   and it is what makes 4-8 workers per socket scale (NOTES section 9 measured a
   1.4x slowdown at 8 workers, which is this bandwidth). Effort 2-3 days, medium.
   Items 2+3 together: 530 -> ~375 ms, 1.4x. That is the CPU ceiling without new
   hardware; NOTES section 6's "a few percent each" was right.
4. **Eigenvalue-free projections: no.** Composite polynomial filters (Kang, Han,
   Groudiev, Yang 2025, arXiv 2507.09165; Polar Express, Amsel et al. 2025) need
   22-31 K x K x K GEMMs = 420-590 GF per projection against our 34 GF warm-started
   subspace step, for ~5e-4 relative accuracy, which section 2 shows is already a
   visible floor. Newton-Schulz is worse (slow start at |lambda|/||W|| ~ 1e-3). They
   only pay when a full-K eigendecomposition is the alternative, i.e. for the
   refresh, which is 7% here. They could keep the certificate valid if the negative
   part is formed as -QQ' (PSD by construction), at the price of one more K^3 syrk.
   Randomized/Lanczos partial projections: neither side is small (nneg 460-720,
   npos 1480-1650 over a run; at the optimum nneg = rank of V'Z*V ~ 600 and
   rank R* ~ 1500 by complementarity), and the smaller side is already the tracked one.
5. **Structure: nothing left to exploit.** The facility/location transposition is
   not a symmetry (A != B), the gangster zeros are 4% of Y and vanish under V, and
   R has rank ~1500 (70% of K), so no low-rank form of Y or Z exists and no fixed
   basis makes W nearly block-diagonal. The Kronecker form of V does make forward
   and expand expressible as GEMMs ((n1 m) x m by m x m, 2 m^5 ~ 0.5 GF each), which
   matters only for a GPU port.
6. **Intel AMX bf16** (this machine has it; the linked OpenBLAS has no sbgemm, oneDNN
   and MKL do): 1024 flop/cycle/core ~ 2 TF/s peak, ~16x AVX-512 fp32. Filter 196 ms
   -> ~25 ms, 530 -> ~360 ms (1.5x) if bf16 were acceptable, and section 2 says it
   drifts. Granite Rapids adds AMX-fp16, which would be safe. Effort small (convert
   operands, call gemm_bf16bf16f32). Confidence: speed high, accuracy low.
7. **GPU, many nodes per device: the only 5-20x per iteration.** Everything maps to
   cuBLAS/cuSOLVER plus two elementwise kernels. Per node-iteration on an A100:
   filter 4 x 5.1 GF in TF32 (156 TF peak, say 60 TF achieved) ~0.4 ms; RR sgemm,
   CholQR, small syevd fp32 ~1 ms; form R dsyrk 2.1 GF fp64 (9.7 TF) ~0.3 ms;
   transforms + Y/Z as GEMM/elementwise, ~1.5 GB traffic at 1.5 TB/s ~1 ms;
   certificate dpotrf 3.2 GF fp64 ~0.5 ms every 10 (GL LAPs stay on the host,
   ~10 ms per certificate, overlapped). The refresh is the risk: cusolver syevd at
   K = 2117 runs at <3% of peak (Zhang et al. 2024 report this for A100/H100), so
   100-200 ms each, 5-10 ms/iter amortized, more than the rest; do it asynchronously
   on the host or replace it by two extra filter steps once lam_hi is tight.
   Estimate 4-8 ms per node-iteration against 530 ms: 65-130x per device, provided
   8-16 nodes run concurrently (a 2117 x 564 GEMM alone under-fills 108 SMs; state is
   ~470 MB per node, so a 40 GB card holds ~80). Per dollar: A100 at $0.6-1.4/hr
   (neo-cloud spot/on-demand, 2026 price surveys; hyperscalers ~$3.4) equals 65-130
   Xeon cores at $0.03-0.05 per core-hour = $2-6/hr: 2-8x per dollar, H100 1.5-2x
   more. Effort 3-6 weeks including validation (bitwise reproducibility across
   CPU/GPU is lost; the certificate stays rigorous because the Cholesky check is a
   posteriori and can run on the host in double). Confidence medium.
8. **Apple silicon per dollar.** A 12-core M4 Pro Mac mini (~$1.5-2k) at ~130-180
   ms/iter per worker gives ~65 iter/s for ~$0.1/hr amortized over two years, about
   10x a cloud Xeon core per dollar, but it is not a 100-core cluster and the AMX
   unit is shared per cluster (the 1.4x slowdown at 8 workers).

## 4. Ranking and first experiment

1. GPU port with 8-16 nodes per device (item 7): 65-130x per device, 2-8x per
   dollar; the only route to 5-20x per iteration. TF32 for the filter is cleared by
   section 2.
2. Tight lambda_max with degree 3 (item 1): free, exact trajectory, 2% per iteration
   and 5-10% of iterations; also removes the refresh-frequency dependence (item 2's
   refresh-40).
3. CPU trims, items 2+3: 1.4x, a few days, certain.

**First experiment (cheap, one afternoon):** port the patch's `DNN_TIGHT_LAM`
power-iteration bound into `proj_psd_subspace` properly (shift by the smallest Ritz
value, 1.1x safety, reset from the top eigenvector at each refresh), then run the
README's node-2 protocol (`18:43,11:29,46:44`, `-H 0 -i 2500`) for degree 3 tight
and degree 4 tight against the recorded degree-4 Weyl (202.3 short at 2500) and
exact values, plus `-r 40`. If degree 3 tight matches 202.3 at 2500 and is ahead at
800, make it the default and re-time node 11 (fathom at 290 iterations) for the
iteration saving. The same run with `DNN_ROUND_BITS=13` costs nothing extra and
confirms the TF32 clearance at convergence before any GPU work starts.

Sources consulted (search snippets only; arXiv not fetched directly): Kang et al.,
arXiv 2507.09165; Amsel et al., arXiv 2505.16932; Zhang et al., arXiv 2410.02170
(GPU EVD efficiency); Groudiev et al., cuADMM, arXiv 2406.05846; Intel Architecture
Day 2021 AMX figures; 2026 GPU price surveys (jarvislabs, intuitionlabs, spheron).
