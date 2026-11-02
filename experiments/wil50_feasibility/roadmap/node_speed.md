# wil50 in a week: per-node speed, cut-phase cost, parallel efficiency, hardware

Research memo, Nov 2 2026 (one of three; the others cover bound strength and
tree size). Nothing was run for it. **M** = measured (file given), **E** =
estimated (arithmetic on M, or a price or spec from memory to be checked).

## 0. Two corrections to the inputs

- **Per iteration at m = 46-47 the 2.1 GHz Xeon takes 530-560 ms, not 150.**
  530 ms alone at degree 4 (M, `solvers/research_iteration.md` §1), 562 ms four
  at a time, 473-480 ms at degree 3 or degree 3 with tight λ (M,
  `solvers/accel/README.md`), 520-580 ms at m = 46 (M,
  `../sko49_feasibility/NOTES.md` §2). The M2 Pro takes 129-132 ms alone (M,
  `NOTES.md` §9). So at m = 47 the M2 Pro is **about 4x** a Xeon core alone and
  about 3x with 8 workers (1.4x contention, M). The 2.4x figure is m = 28
  (30 against 73 ms, M, `RUNNING.md`).
- The per-node figures hold anyway, because nodes are short. The run_ab3 base arm
  averaged **63 Xeon-s per DNN node** (139,054 s over 2,193 nodes, M,
  `stronger_bounds/search/wil50n10ext_base.log`). On all32, fathoms averaged
  105 s and give-ups 95 s at H = 200 or 175 s at H = 300 (M,
  `tuning/README.md`).

## 1. Where a dominant-class CPU-second goes

The run_ab3 cut arm (M, `wil50n10ext_cut.log`) spent 76,755 s in all. Of that,
**41,796 s (54%) went to plain DNN bounds over 365 nodes and 34,959 s (46%) to
39 cut phases, about 900 s each.** Of the DNN part, I estimate about 7,000 s
for the 39 give-ups and about 35,000 s for 326 fathoms (E). About 320 of the
365 nodes are children of the 8 phases that failed (E, from the tree shape:
1 + 46 + 8 x 45).

Inside a fathoming phase (M, `stronger_bounds/laptop/README.md` §2):

- 400 uncut iterations;
- 3 to 12 rounds of 300 iterations;
- a quiet round of 960-1,570 iterations, which is 25-50% of the phase.

Fewer rounds is cheaper as long as the quiet round starts from the best round's
multipliers. On d1, 3 rounds cost 689 s and 12 rounds cost 1,215 s (M).

Inside an iteration (M, `research_iteration.md` §1):

| part of the iteration | share | note |
|---|---|---|
| Chebyshev filter and Rayleigh-Ritz | 56% | sgemm at 103 GF/s, about 75% of AVX-512 fp32 peak, compute-bound |
| forming R (dsyrk) | 10% | |
| ssyevd refresh | 7% | |
| memory-bound transforms and Y/Z-step | 25% | |
| certificate | 3% | |

## 2. Ranked levers

| # | idea | factor on the wil50 total | effort (days) | risk | evidence | cheapest first test, and the result that kills it |
|---|---|---|---|---|---|---|
| 1 | **Run the filter and Rayleigh-Ritz GEMMs on Intel AMX.** On Sapphire Rapids, use split-bf16 (two or three bf16 passes, fp32 accumulate) through oneDNN or oneMKL. On Granite Rapids, use one-pass AMX-FP16. | 1.4-1.7x per node on AMX x86 (E) | 2-4 | medium | A 10-bit mantissa in the filter operands is indistinguishable from fp32 to iteration 200; one-pass bf16 drifts 7 units (M, research_iteration §2). This container has `amx_bf16` (M, lscpu). OpenBLAS has no bf16 GEMM. AMX peaks near 2 TF/s per core against the measured 103 GF/s (E). | One hour: time a 2117x2117x600 product with OpenBLAS sgemm and with 2-pass and 3-pass split-bf16. Then run node 2 to 2,500 iterations, then all32. **Kill:** effective throughput below 250 GF/s; node 2's converged 48,610.1 moves by more than 1; or any all32 fathom/miss flips. |
| 2 | **CPU trims.** (a) With tight λ, the Weyl difference and the `Wprev` copy at the top of `proj_psd_subspace` are computed and thrown away: about 16 ms, free. (b) ssyrk for R: -28 ms. (c) guard K/40: -25 ms. (d) One-pass CholQR: -15 ms. (e) Refresh every 40-60 under tight λ: -18 to -25 ms. This is untested: `-r 120` was measured only with the Weyl bound. (f) Fused, memory-lean transforms: -50 ms. | 1.3-1.4x per iteration (E). Together with #1: 530 → 220-280 ms, **1.8-2.2x per node** (E). | 3-5 | low | Phase breakdown (M, research_iteration §1, §3) | `-v -v` profile of node 2 before and after, then all32. **Kill:** less than 15% per node. |
| 3 | **GPU port, 8-16 nodes per device.** TF32 filter. Certificate in fp64 on the device, or on host cores. | Per device: 60-160 Xeon cores (A100), 100-240 (H100), 40-80 (RTX 4090) (E). **Per dollar: 1-5x x86 Spot at neo-cloud prices; at best even at hyperscaler on-demand prices.** | 15-30, plus validation | med-high | The GEMM count and traffic per iteration (M, research_iteration §1). GPU timings are E. | A one-day CuPy or PyTorch prototype of one eigen-step: K = 2117, block 600, 8-16 streams, on a rented A100, H100 or 4090 (about $10). **Kill:** above 10 ms per node-iteration in aggregate on an A100 (below 50 core-equivalents, where per-dollar parity is lost). |
| 4 | **Cut-phase trims.** Fewer rounds or 200-iteration rounds, with the quiet round from the best multipliers; stop the quiet round when the certificate's gap over 100 iterations shrinks by less than about 0.8; start the rounds at the gap-ratio point instead of iteration 400. | Phases 1.2-1.5x faster, so **1.08-1.2x total** (E) | 2-3, plus 2 of runs | medium | 3 rounds against 12 on d1: 689 s against 1,215 s (M). 150-iteration rounds plus a 600-iteration quiet round fathomed 4 of 4 near misses (M, stronger_bounds step 6). | Laptop protocol on d1, d4, d6, d11, 6, 9, 15, 14 with rounds of 200 or 300, at most 3-6 of them, and the quiet-round abort. **Kill:** any fathom is lost, or the median CPU to fathom falls by less than 15%. |
| 5 | **Children of a failed phase start from its cut pool and best multipliers.** A cut whose three facilities stay free remains valid in the child. The child gets a new L_add, and the certificate is unchanged. | 1.1-1.25x on the dominant class (E). These children are about 40-45% of run_ab3's CPU. | 3-5 | medium | Children today keep the *uncut* multiplier (`bnb.C`). | `-F` search under node 16 (its uncapped phase failed, M): compare the iterations the children need to fathom. **Kill:** less than a 30% cut. |
| 6 | **Parallel machinery.** No global checkpoint barrier; per-unit leases (plan.md); keep the DNN warm start inside units. | Moves cores doing real work from 70-85% to 93-97% (E): 1.1-1.3x against today's runner. The brief's 26-83 days already assume 100%. | 1 (barrier) to 25-40 (coordinator) | low | Busy cores (CPU/wall) were 3.90 of 4 in the base arm and 3.50 of 4 in the cut arm, which had 2-hour checkpoints (M, run_ab3 logs). 10-minute checkpoints cost about 25% (M, stronger_bounds step 6). | Engineering, no kill test. Phase 5 of plan.md accepts above 90% for 24 h. |
| 7 | **Fleet choice by plan.md's Phase-0 benchmark** | ±1.5x in dollars (E) | 3-4 | low | §5 below | Below |
| 8 | Low: H = 200 (-26% on all32, -8% on nug30, M; give-ups are about 9% of cut-arm CPU, so about 1.04x); sibling warm starts; adaptive β (β 0.21 was neutral on wil50, M); certificate cadence (3%); **batching nodes per core (no: the GEMMs are already at 75% of peak and nodes share no data)** | ≤1.1x each | - | - | - | - |

## 3. The top three

**Item 1 + 2: AMX plus trims (1.8-2.2x per node on Sapphire or Granite Rapids
cores, about a week).** The M2 Pro's 4x lead at m = 47 is its AMX sgemm. The
Xeon running these studies has AMX too, and leaves it idle because OpenBLAS
cannot use bf16. Precision is already cleared to iteration 200: 10-bit
operands track fp32 exactly, and only one-pass bf16 (7 bits) drifts. Two-pass
split-bf16 carries about 16 bits, so it should be safe, but this is unmeasured
at convergence; the node-2 run settles it. Neither the filter's precision nor
how R is formed can make a bound invalid: R = W − FF′, and the certificate
recomputes V′Z′V in double from Z′ and checks it with Rump's Cholesky
(`dnn.C` header §2). Lower precision can only slow the iteration, or widen the
Cholesky shift slightly when R is formed in fp32. The double re-projection
path already handles that edge. Rule: any precision for the filter and the
Rayleigh-Ritz, fp32 for R, and fp64 for the Y/Z state, the certificate's
transform, dpotrf and the Gilmore-Lawler sums. Once the GEMMs are fast, the
memory-bound 25% becomes about half the iteration. Trim (f) is what keeps 32-64
processes per socket from sharing a bandwidth ceiling (measured 6 GB/s per
process in the transforms, so about 200-380 GB/s per host, E).

**Item 3: the GPU port, which is a wall-clock lever more than a cost lever.**
The A100 numbers below are estimates (E).

- **GEMMs:** about 28 GF of fp32 GEMM per node-iteration (degree 3 plus
  Rayleigh-Ritz plus CholQR2) takes 0.3-0.6 ms in TF32.
- **Transforms and Y/Z-step:** about 0.4 GB of HBM traffic each, a few tenths
  of a ms.
- **The risk is the dense eigendecompositions:** the 600x600 Rayleigh-Ritz
  ssyevd every iteration, and the 2,117 refresh every 20 iterations. The
  refresh takes 734 ms on a Xeon core (M); cuSOLVER runs at a few percent of
  peak at that size (E). Plan to put both on the host's 12-24 vCPUs per GPU,
  or to replace the refresh with extra filter steps once λ is tight.
- **Memory:** about 0.57 GB per node at m = 47 (100 m⁴ bytes plus the warm
  start, `dnn.h`). An 80 GB card holds about 130 nodes, a 24 GB card about
  40. Saturating a device takes 8-16.
- **Speed:** 3-8 ms per node-iteration on an A100 and 2-5 ms on an H100, so a
  device does 125-500 node-iterations per second against 2 for a Xeon core.

**So is a 1,000-core week better bought as 100 GPUs for a day?**

| | compute | cost |
|---|---|---|
| 1,000 Xeon cores for a week | 168,000 core-hours | $3.7k-7.6k on x86 Spot, $1.7k-2.5k on Graviton Spot (at parity) |
| 100 A100s for a day | 144,000-384,000 core-equivalent-hours | $1.9k-9.6k |

The money is about the same; the wall time is 7x shorter. That is "yes"
whenever the week's wall time is the constraint, provided the port (3-6 weeks)
passes the kill test. Consumer RTX 4090s at $0.3-0.6 an hour (E) are the
cheapest per core-equivalent, with an fp32-only device path and the fp64
certificate on the host (143 ms every 10 iterations, M, so 14 ms of host time
per node-iteration).

**Items 4 + 5: fewer iterations per decision in the cut search (1.2-1.5x on
the dominant class).** About 46% of the dominant class's CPU is cut phases,
and a fathoming phase spends 25-50% of it in the quiet round. The measured
shortcuts all point the same way: 3 rounds beat 12 when the quiet round starts
from the best multipliers, and short rounds work once a quiet round follows.
The 8 failed phases and their roughly 320 children (about 40-45% of the arm's
CPU) are where item 5 applies. Today those children start from the parent's *uncut*
multiplier, 10-60 short. Started from the parent's cut multipliers, which sit about
one level higher, they should fathom in far fewer iterations. Every proposed stop
rule is a give-up, so each must be scored on the fathom/miss set, because a
false failure costs 45 children.

## 4. Parallel efficiency at 1,000 cores

`bench/parallel.py` loses time in three places:

1. **Checkpoints:** a barrier makes 1,000 workers wait for the slowest node,
   up to a 30-minute phase. At a 2-hour cadence that costs about 15-25% (E).
2. **Steals:** they wait for the victim's node to end. This matters only at
   the end.
3. **Lost warm starts:** stolen nodes start cold, because the 39-78 MB DNN
   warm start is not in state files. That costs 13-16% more DNN time (M,
   plan.md §2.1).

plan.md's coordinator removes all three: leases instead of a barrier, units of
30 minutes to 4 hours, and a pool kept above 2x the fleet. The tail is the
depth-3-to-7 critical path, about 3-5 hours or 2-3% of a week. Spot
interruptions lose about 230 s per slot each, which is negligible (E).

**Real work: 93-97% of the cores with the coordinator, 70-85% with today's
runner (E).** The open risk is memory bandwidth with 32-64 processes per host;
plan.md's Phase 0 measures it.

## 5. Hardware

| hardware | per-core speed on this solver | price and verdict |
|---|---|---|
| **Sapphire / Granite Rapids** (c7i, m7i; c8i, GCP C4; E) | AMX on every core (Apple's AMX is shared per cluster); with item 1 about 2x | Fastest week per core |
| **Graviton4** | About 90 GF/s fp32 peak per core, close to the measured Xeon sgemm (E). No safe matrix engine: one-pass BF16 MMLA drifts. | Lowest $ per core-hour; best per dollar without AMX |
| **AMD Zen 5** | Full-width AVX-512 (E) | Worth a Phase-0 row |
| **M-series Macs** | 4x a Xeon core per iteration (M) | About $0.06 per core-equivalent-hour on AWS mac2 hosts, worse than x86 Spot (E). Parity with 1,000 cores takes about 40 owned Macs. Not the fleet. |

## 6. Dollars and days

Basis: 2.2e9-7e9 reference-core-seconds today with the cut phase (the brief's
figure). Spot prices are from plan.md §10.1; GPU prices are E.

| scenario | 1,000 cores | x86 Spot | Graviton Spot (parity) | 100 GPUs |
|---|---|---|---|---|
| today, `parallel.py` (75-85% busy) | 30-108 d | - | - | - |
| today, coordinator (95%) | **27-85 d** | **$14k-92k** | $6k-31k | A100: 1.8-15 d, $3k-54k (neo-cloud) or $13k-144k (hyperscaler); H100: 1.2-9 d, $4k-65k |
| items 1, 2, 4, 5 and 6 on AMX x86 (2.1-3.3x, the product of the item ranges) | **8-40 d** | $4k-43k | items 2, 4 and 5 only (1.5-2.1x): $3k-20k | - |
| GPU plus items 4 and 5 (1.2-1.5x) | - | - | - | H100: **0.8-7.6 d**; A100: 1.2-12.6 d |

The target of 6e8 core-seconds is 7.3 days on 1,000 cores at 95% busy. After
my CPU items, the bound-strength and tree angles still have to supply 1.1-5.5x.
With 100 H100-class GPUs instead, today's tree fits in a week at its low and
middle estimates, and needs only up to about 1.3x from the other angles at its
high end.

## 7. Interactions with the other two angles

- **Bound strength via cuts moves CPU into long phases.** Phases are 46% of the
  CPU today and would be more if cuts run on more nodes. So items 1-3 apply in
  full, and items 4-5 grow in value. A new cut family that is more expensive
  per round makes the GPU's GEMM advantage bigger; separation stays on the
  host.
- **A stronger bound at depths 3-4 shifts work to larger m**, roughly
  (m−1)⁶ per iteration, so a depth-3 node costs about 1.3x a depth-5 one. It
  also makes the granularity coarser: the depth-3 frontier (2.9e4 nodes) still
  feeds 1,000 workers, and the tail stays at hours.
- **Every iteration-saving stop rule is a tree risk.** A false give-up or a
  false phase failure costs 45-47 children, so items 4, 5 and 8 must be scored
  on all32 plus a depth-4 sample. The gap-ratio rule at 0.65 already failed on
  nug30 (M).
- **Rounding differences across AMX, GPU and MKL** can flip decisions within
  about 1e-6 of the cutoff. That changes node counts, not validity; plan.md's
  single-family policy A applies.
- **The GPU port must return the primal marginals x̄** (row 0 of Y) for branching
  rules 8 and 9, the tree angle's free lever.
