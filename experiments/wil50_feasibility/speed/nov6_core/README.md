# The core calculation and the dual solution: a cheaper DNN iteration, a certificate sooner (Nov 6, cloud)

One of three Nov 6 studies before the 1,000-core run; this one covers the cost of
a DNN iteration and the certified bound (the two siblings cover ADMM convergence
and nonlinear cuts). **[m]** = measured here (log named), **[e]** = estimate.
Host: the 4-vCPU 2.1 GHz Sapphire Rapids VM with AMX, everything on CPU 3
(`taskset -c 3`, one BLAS thread), the two sibling agents on CPUs 1-2 throughout
(load 1.5-3), so every millisecond below is a loaded figure, and paired
(interleaved) runs decide every comparison. Production configuration = the
`fast` arm of `../../pipeline_ab/run_pair.sh` as far as `dnntest` expresses it:
degree 3, tight lambda, PRSM (1, 1), refresh 40, one-pass CholQR, AMX, R's rank
term in fp32 off certificate iterations, the two-phase beta.

## Verdict

- **Three opt-in switches, together 1.11x per iteration on the 32-node sample
  with identical decisions on every node [m]; every default bit unchanged**
  (`make test`, `make dnncheck`, `make shipcheck`, `make continuationcheck`,
  `bench.py check` exact; the old and new binaries print identical traces
  with the switches unset, `logs/check_identity.txt`):
  - `DNN_CORE_LAZYCERT=1`: Rump's Cholesky test (forward + dpotrf, about 130
    of the certificate's 141 ms) only at checks that can decide something.
    **The same fathom decisions and the same returned lb by construction**;
    **-12 ms of 304 per iteration on all32 [m]** (the profile's estimate),
    the same 18 fathoms at the same iterations, the same give-ups, the
    misses' bounds to the printed digit. Every CPU family.
  - `DNN_CORE_NOSYMV=1` (with the existing `DNN_NOSYMW=1`): no symmetrisation
    pass after `expand`; the certificate reads both arrow entries of M and
    symmetrises its own transform, so it stays a bound for sym(Zc). The
    iterate is symmetric to rounding only, so the trajectory changes at the
    rounding level (not bit-identical): **-14 ms per iteration [m]**, the
    same decisions on all 32 nodes, the misses' bounds within 0.06.
  - `DNN_CORE_REFRESH_PARTIAL=1`: the full refresh as ssytrd + sstedc + a
    back-transformation of the nneg + guard + 1 columns the tracker reads
    instead of all K. **Bit-identical eigenvectors** (`logs/bench_refresh.txt`:
    0.00e+00 difference; traces identical, `logs/check_identity.txt`), the
    refresh at 66-69% of ssyevd [m]: **-4 ms per iteration on all32 [m]**
    (-7 to -9 estimated; more on short warm solves, which refresh 3 times in
    100 iterations). Under the 10 ms line, kept because it is free (nothing
    can change but the clock).
  - On the production search (nodes mostly warm 100-iteration solves plus
    cut phases, both certified every 10 iterations and refreshed every 40)
    the same arithmetic gives **about 1.10-1.12x in DNN CPU [e]**; measured
    only cold, on all32. `QAP_DNN_REFRESH=80` adds 2.6% more on all32 at the
    price of trajectory changes (section 3) and needs the search A/B.
- **AMX kernel knobs** (`DNN_CORE_AMX_IB/KB/NT/ORDER`, bit-identical T for
  any setting): no setting beats today's blocking beyond this host's noise;
  a 256-wide i-slab is 1.5-2x slower (the W slab leaves L2). The kernel's
  "last third" stays open and needs a quiet machine (section 4).
- **Closed here** (section 5): GEMM-only projections (sign/polar/composite
  polynomials, 2025-2026 versions included), the a-priori certificate without
  Cholesky, AMX for the form-R syrk, CholQR by explicit inverse, mixed-
  precision Jacobi for the Rayleigh-Ritz eigensolve, shipping the parent's
  eigenbasis, on-the-fly L, fp32 transforms (estimated, deferred), structure
  (Kronecker, grid symmetry, facial basis), batching siblings, a tighter
  lambda in the certificate (lambda_max(V'ZcV) is exactly 0 on R's range).

## 1. Profile of a production iteration (node 10 at m = 47 and friends)

`profile.sh`, `dnntest -v -v`, 300 iterations each; `logs/prof_*.log`,
table by `profile_table.py` (`logs/profile_table.txt`). ms per iteration, loaded:

| piece | cold, sgemm (no AMX) | cold, AMX | cold, production | warm child (Z + conditioned Y), production | depth-4 fathom (all32 node 19, m = 46) |
|---|---|---|---|---|---|
| W*Y products (+ pack) | 155.3 | 74.8 | 76.4 | 83.8 | 52.4 |
| CholQR (one pass) | 20.2 | 20.6 | 20.7 | 23.9 | 11.6 |
| Rayleigh-Ritz: Y'T, ssyevd(kb), X = YG | 48.9 | 49.9 | 49.1 | 57.7 | 25.6 |
| convert W to fp32 + Chebyshev updates | 6.2 | 6.3 | 6.3 | 6.8 | 4.9 |
| full ssyevd refresh | 19.7 | 20.6 | 20.2 | 19.2 | 19.1 |
| form R (F, syrk, add-and-mirror) | 50.3 | 51.2 | 33.0 | 35.6 | 24.1 |
| forward V'(Y + Z/beta)V | 39.5 | 39.5 | 39.8 | 40.7 | 36.0 |
| expand V R V' | 28.8 | 28.7 | 27.2 | 27.1 | 23.9 |
| Y/Z-step | 23.2 | 23.3 | 24.6 | 23.6 | 22.4 |
| certificate (every 10) | 13.5 | 13.7 | 14.1 | 13.1 | 12.2 |
| **total** | **416.3** | **332.1** | **318.7** | **335.1** | **237.2** |
| iterations, K, mean kb / nneg | 300, 2117, 610 / 510 | 300, 2117, 610 / 510 | 300, 2117, 613 / 511 | 300, 2117, 692 / 586 | 100, 2026, 434 / 349 |

Shares of the production iteration: the W*Y products 24%, the rest of the
projection (CholQR, Rayleigh-Ritz, refresh, convert) 30%, form R 10%, the
transforms and the Y/Z-step 29%, the certificate 4.4%. The warm child (the
search's common case) tracks a larger block (kb 692 against 613: the
conditioned Y starts near the final rank structure), so its projection is 10%
dearer per iteration while it needs fewer iterations. The depth-4 fathom runs
3 refreshes in 100 iterations (8% of its CPU) and 10 certificates (5%).

Achieved rates (production cold run; peaks at the nominal 2.1 GHz, which the
core exceeds under turbo, so the fp32 figure reads above 100%):

- one W*Y product on AMX: 18.4 ms for 5.5 GF fp32-equivalent = 299 GF/s-equivalent = 896 GF/s of bf16 tile products, **42% of the tile peak** (2,150 GF/s at 2.1 GHz; `../../amx/amx_peak.c` measured 2,100-2,280). With OpenBLAS sgemm: 39.9 ms = 138 GF/s, the fp32 peak.
- CholQR (ssyrk + spotrf + strsm, 1.6 GF): 21.3 ms = 75 GF/s (strsm is the slow half).
- Rayleigh-Ritz Y'T and X = YG (1.6 GF each): 13 ms each = 121 GF/s. ssyevd(kb = 613) (about 0.77 GF): 24.7 ms = 31 GF/s (half Level-2).
- full ssyevd(K = 2117) (about 41 GF): 758 ms = 54 GF/s; 8 refreshes in 300 iterations.
- form R's rank term: fp32 ssyrk (2.3 GF) 20.4 ms = 112 GF/s; the fp64 dsyrk of the default 36.0 ms = 64 GF/s, the fp64 peak.
- the memory streams: forward about 357 MB in 39.8 ms (9.0 GB/s), expand about 195 MB in 27.2 ms (7.2 GB/s), the Y/Z-step about 238 MB in 24.6 ms (9.7 GB/s). One core of this VM streams about 9-10 GB/s under the siblings' load.
- the certificate: 141 ms each (30 per 300 iterations): forward 40 + dpotrf(K) (3.2 GF fp64, about 85) + the Gilmore-Lawler LAPs and norms (about 15).

So the compute-bound pieces run at or near their units' peaks except the AMX
kernel (42%) and the two small LAPACK eigensolves (Level-2 bound); the
memory-bound third runs at the VM's streaming rate. What is left to take is
(a) arithmetic that need not be done (the certificate's Cholesky at checks
that cannot fathom, the refresh's unread eigenvectors, the symmetrisation
passes), (b) the kernel's other 58%, and (c) memory traffic.

## 2. Candidates, ranked (arithmetic at m = 47, production iteration 319 ms)

| # | idea | ms per iteration | CPU families | status |
|---|---|---|---|---|
| 1 | **Lazy certificate** (`DNN_CORE_LAZYCERT=1`): the polyhedral value first; Rump's test only when it could fathom, at the first and last check, and at the exit on the stopping check and the stored best | -12 [e]: 130 ms x (27 of 30 checks skipped on a 300-iteration miss; 8 of 10 on a 100-iteration warm fathom) / 10. **all32: -12 ms, -3.9% CPU, identical decisions [m]** | all | built, measured, passes |
| 2 | **NOSYMV + NOSYMW**: no symmetrisation of VRV (78 MB of traffic) and of W (72 MB), Wf written in forward's gather | -10 and -10 [e] (the amx study measured NOSYMW at -6 -4); **all32: -14 ms, -4.6% CPU on top of 1, identical decisions, misses' lb within 0.06 [m]** | all; more under load (bandwidth) | built, measured, passes |
| 3 | **Partial refresh** (`DNN_CORE_REFRESH_PARTIAL=1`): sormtr on kb + 1 columns instead of K | -7 [e per refresh: 758 -> about 520 ms [m]; 8 refreshes per 300 iterations]; -9 on a 100-iteration warm solve (3 refreshes); **all32: -4 ms, -1.6% CPU on top of 1 + 2 [m]** | all | built; bit-identical, so only the clock can change |
| 4 | AMX kernel: the other 58% of the tile peak | up to -40 if the pattern's measured ceiling (2.1 TF with W in L2, Y in L1, `../../amx/amx_loads.c`) were reached; knobs measured here: none beyond noise | Sapphire/Granite Rapids only | open; needs a quiet machine |
| 5 | fp32 storage for VRV and T (the transforms' scratch): about -115 MB of the 790 MB per iteration | -14 [e]; the iterate changes at 1e-7 relative, like R32 | all; most under load | not built: touches forward, expand and the Y/Z-step, which this study was to leave to the siblings' merge |
| 6 | Y/Z-step on the upper triangle with tiled transposed writes (-78 MB) | -8 [e] | all | killed by the 10 ms line |
| 7 | form R's syrk on AMX (split bf16; the ADMM step only) | -8 [e]: 2.3 GF at 300 GF/s-equivalent + packs against 20.4 ms | AMX only | killed by the 10 ms line |
| 8 | CholQR with the explicit inverse (strtri + sgemm instead of strsm) | -5 [e] | all | killed by the 10 ms line |
| 9 | on-the-fly L in the Y/Z-step (bit-identical; L is 39 MB read per iteration; not in cut phases, whose L_add is dense) | -4 [e] | all | killed by the 10 ms line |
| 10 | Rayleigh-Ritz ssyevd(kb) by a mixed-precision or warm-started Jacobi sweep (H is nearly diagonal after a warm start; arXiv 2509.00495, 2607.23187) | -5 to -10 [e]: one cyclic sweep is 2 kb^3 = 0.46 GF; two are needed near the zero cluster | all | killed: at best the 10 ms line, a day of work, and the cluster at zero (`../cs2/README.md`) slows Jacobi too |
| 11 | refresh cadence 80 instead of 40 under tight lambda (untested since `../../solvers/research_iteration.md`; a flag, no code) | -10 cold, -8 on a 100-iteration warm solve [e] if no iteration is lost; **all32: -11 ms, -2.6% CPU on top of the bundle, same 18 fathoms, one fathom a check earlier, two give-ups 30-70 iterations later, misses' lb at 200 down to -0.90 [m]** | all | passes the all32 gate narrowly; the search A/B decides |

Ideas 1-3 together: about -29 ms (-9%) [e] at m = 47; **measured -30 ms (-9.9%
per iteration, -10.0% CPU) on all32** (section 3).

**Kill rules applied.** The study's line was 10 ms per iteration unless the
change is free (bit-identical). Killed by it: ideas 6-10 (the triangle
Y/Z-step, the form-R syrk on AMX, CholQR by explicit inverse, on-the-fly L, a
Jacobi Rayleigh-Ritz), each 4-8 ms [e] and each a trajectory or code-path
change. Kept under the line: idea 3, bit-identical. Killed by measurement:
every AMX kernel knob (section 4: within the host's noise, or slower). Deferred
rather than killed: idea 5 (fp32 VRV and T, -14 ms [e]), because it rewrites
the transforms and the Y/Z-step, which this study left alone so that its
switches merge with the siblings' work. Decision gates for the switches:
the same fathom set, the same iteration on every fathom, the misses' lb at 100
and 200 within 1 of base's (all three switches: identical, or within 0.06),
`make dnncheck` and brute force with the switch exported (section 6 and
`logs/gates/`).

## 3. Measurements on single nodes and on all32

### Identity and single-node checks (`check_identity.sh`, `logs/check_identity.txt`, `logs/ident/`)

Old binary (ab937ca) against the new one, production flags, on CPU 3 [m]:

| node | old vs new, switches unset | `REFRESH_PARTIAL` vs old | `LAZYCERT`: lb, outcome, Rump tests | `LAZYCERT + NOSYMV + NOSYMW` |
|---|---|---|---|---|
| nug20 root, 200 iterations | identical trace and RESULT | identical | same lb to 17 digits (2504.7209480473393), 2 tests instead of 20 | lb 2504.72128 (+1.6e-4: another trajectory) |
| wil50 all32 node 19 (m = 46, fathoms) | identical | identical | fathomed at 130 as before, same lb (48830.419745245599), 2 tests instead of 13; 32.0 s against 32.8 | fathomed at 130, lb 48830.434, 30.1 s |
| wil50 node 10 (m = 47, give-up at 300) | identical | identical | same lb (48582.462771451275), hopeless at 300, 2 tests instead of 30; 91.1 s against 97.2 | hopeless at 300, lb 48582.424 (-0.04), 87.4 s |

The lazy certificate's returned lb equals the full certificate's to the last
digit on both wil50 nodes: the stopping check (the fathom at 130, the give-up
at 300) is decisive and certified in full, and the stored best never beat it.
Single-run CPU figures carry this host's +-5% noise; the all32 run decides.

**What the gate caught.** The first build passed brute force but failed
`make dnncheck` part C: two nug20 nodes the default fathoms at 330 gave up at
200 under the switch (`logs/gates/dnncheck_lazycert_v1_failed.log`). The give-up rule's
extrapolation still started from `best`, which under the switch is the last
*certified* value (the first check's), not the running estimate; the fix is one
identifier (`best_rule`, equal to `best` without the switch), after which the
node fathoms at 330 with the same lb to 17 digits and 2 Rump tests instead of
33 (`logs/gates/dnncheck_lazycert_v2.log`). The all32 arms with the switch were
restarted from the fixed binary (the `base` arm's one finished node was kept:
the default path is bit-identical).

### all32 A/B (`run_all32.sh`, `summarize.py`; `runs/`, `runs/summary.txt`)

The 32 wil50 nodes of `../../tuning/giveup/all32.txt` (16 at m = 47, 16 at
m = 46), cold, production flags (`-d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 2500
-u 48816`, `DNN_AMX=1 DNN_CHOLQR_PASSES=1 DNN_R32_NONCERT=1
DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7`), one process at a time on CPU 3, the
five configurations alternating node by node (base, lazy, lazysym, bundle,
bundle_r80: see `run_all32.sh`), a private copy of the binary, the siblings on
CPUs 1-2 throughout (load about 2). ms per iteration is CPU over iterations
summed over the nodes; the fathom set, the iterations to fathom and the
misses' lb at 100 and 200 are the decision checks.

| configuration | fathomed | ms per iteration | iterations to fathom (18 nodes) | CPU, all 32 | CPU, fathoms | CPU, misses | misses' lb at 100 / 200 against base (median [min, max]) |
|---|---|---|---|---|---|---|---|
| base (production) | 18 | 304 | 3,690, median 190 | 2,449 s | 1,059 s | 1,390 s | - |
| lazy: `DNN_CORE_LAZYCERT=1` | 18, the same nodes | **292 (-3.9%)** | 3,690: every node at the same iteration | **2,354 s (-3.9%)** | 1,006 s (-5.0%) | 1,349 s (-2.9%) | +0.00 [+0.00, +0.00] / +0.00 [+0.00, +0.00] |
| lazysym: + `DNN_CORE_NOSYMV=1 DNN_NOSYMW=1` | 18, the same | **278 (-8.6%)** | 3,690: every node at the same iteration | **2,240 s (-8.5%)** | 957 s (-9.6%) | 1,283 s (-7.7%) | -0.00 [-0.06, +0.00] / -0.00 [-0.04, +0.06] |
| bundle: + `DNN_CORE_REFRESH_PARTIAL=1` | 18, the same | **274 (-9.9%)** | 3,690: every node at the same iteration | **2,205 s (-10.0%)** | 943 s (-11.0%) | 1,262 s (-9.2%) | as lazysym (bit-identical trajectory) |
| bundle_r80: bundle with `-r 80` | 18, the same | 263 (-13.5%) | 3,680: node 31 at 210 instead of 220, the other 17 the same | 2,147 s (-12.3%) | 906 s (-14.4%) | 1,242 s (-10.6%) | +1.16 [-4.98, +2.20] / -0.16 [-0.90, +0.27]; two give-ups later (node 9: 400 instead of 330, node 16: 400 instead of 370) |

- **Lazy certificate: -12 ms per iteration [m], the profile's estimate to the
  millisecond; decisions identical on all 32 nodes** (fathom set, every
  fathom's iteration, every give-up's iteration, the misses' lb at 100 and
  200 to the printed digit). Of 141 ms per certificate it saves about 95:
  forward and dpotrf, less the 39 MB copy of the multiplier at checks whose
  estimate makes a new maximum (nearly every check while the bound rises).
- **NOSYMV + NOSYMW: -14 ms per iteration [m]** (estimate -20), decisions
  identical; the misses' bounds move by at most 0.06 (the rounding-level
  trajectory change). The fathoms gain most (-9.6% against -7.7%): their
  iterations are shorter and the transforms a larger share.
- **Partial refresh: -4 ms per iteration [m]** (estimate -7; fewer refreshes
  per iteration on the fathoms than the 8 per 300 of a miss), bit-identical.
- **The three together: 1.11x per iteration and 1.11x in CPU on all32 with
  the same decisions on every node.** Every CPU family: none of the three
  touches the AMX path (the gain is a larger fraction without AMX, whose
  iteration is 100 ms longer for the same certificate and transforms).
- **Refresh 80 under tight lambda: -11 ms more [m] but not free**: node 31
  fathoms a check earlier, two misses give up 30-70 iterations later (their
  bound at 200 is within 0.3 of base's, so the extrapolation rule sees a
  flatter tail), the misses' lb at 200 is lower by up to 0.90 (the kill line
  is 1). It passes the all32 gate narrowly; its effect on the search's give-up
  pattern needs the node-10 subtree A/B before it is adopted. It is a flag
  (`QAP_DNN_REFRESH=80`), not code.

## 4. The AMX kernel: knobs measured, what is left

`scan_amx.sh` (`logs/scan_amx.txt`) and `pair_amx.sh` (`logs/pair_amx.txt`)
time one W*Y product at K = 2117, kb = 610 with `../../amx/test_amx_gemm.c`
against the kernel's new runtime knobs (`latest/amx_gemm.c`; every setting
gives the same T bit for bit, since each accumulator tile sees the same k
order):

| knob | setting | ms per product, min / median of 3 (base interleaved: 13.8 / 20.8 over 33 runs) |
|---|---|---|
| `DNN_CORE_AMX_NT=0` | W tiles with plain loads instead of the non-temporal hint | 15.8 / 21.1 |
| `DNN_CORE_AMX_ORDER=1` | the two accumulators' six products interleaved | 12.9 / 14.8 (paired run below) |
| `DNN_CORE_AMX_IB=64` | i-slab 64 (W slab 0.55 MB in L2) | 16.9 / 20.1 |
| `DNN_CORE_AMX_IB=256` | i-slab 256 (2.2 MB: larger than L2) | 29.2 / 29.7 |
| `DNN_CORE_AMX_KB=4`, `16`, `67` | k-block 4, 16, or the whole K (C tiles stored once) | 21.6 / 24.7, 23.2 / 29.9, 21.1 / 22.8 |

The run-to-run swing of a product on this VM is 14-30 ms (the amx study saw
2-3x on `amx_peak.c`), so only the slab size has a clear effect: the W slab
must stay inside L2. The paired base/ORDER=1 run (`pair_amx.sh`, 15
alternations at load 3.0): base min 13.9 / median 17.8 ms, ORDER=1 min 13.1 /
median 17.7, paired difference median -0.6 ms, the variant faster in 9 of 15
[m]: **within noise** (at most -3%, 2 ms of a 76 ms piece). The knobs stay in
the kernel for a quieter machine; the defaults are unchanged.

Why the kernel sits at 42% of the tile peak, and what the literature says:
the fused pattern (6 tile loads, 6 `tdpbf16ps`, two accumulators) runs at
1.84-1.87 TF with both operands in L1 and 2.1 TF with W from a 1 MB L2 buffer
(`../../amx/amx_loads.c`), so the loss is in the real access pattern: C tiles
stored and reloaded at stride K per k-block (about 90 MB of strided traffic
per product), the Y k-block and the W tiles competing in a 48 KB L1, and TLB
reach over the 18 MB packed W. The Sept 2026 AMX emulation paper (arXiv
2609.04663: three bf16 components, six products, Xeon Platinum 8462Y+,
oneMKL 2026.1) reaches 1.14-2.56x over oneMKL SGEMM depending on size (1.39x
at n = 2048) and names the eight-tile register file (operand panels and fp32
accumulators at once) as the limit of the microkernel: the published state
of the art is at or below our fused three-product kernel (2.0-2.3x over
sgemm in the microbenchmark, 2.17x in the solver's product), so no known
kernel pattern gives the other 58% for free. A quiet Sapphire Rapids core
(the laptop cannot; an AWS c7i can) is where the kernel's cache behaviour
can be tuned: the test is `pair_amx.sh` with `perf stat` on L1/L2 misses.

## 5. Closed here, with the arithmetic

- **GEMM-only PSD projection (matrix sign, polar, composite polynomials;
  Newton-Schulz, Polar Express, Zolotarev/QDWH).** Checked against the 2025-
  2026 literature: Polar Express (arXiv 2505.16932v5, 2025), the Muon
  Newton-Schulz line (2609.39595, 2610.26742 "Polar Switch", 2610.34158,
  2026), the factorisation-free composite-polynomial PSD projection (Kang et
  al., arXiv 2507.09165v2, 2025). All need 20-30 K x K x K products for a
  relative accuracy of 1e-3 to 1e-4 on a spectrum with hundreds of eigenvalues
  within 1e-3 of zero (`../cs2/README.md`): 20-30 x 19 GF = 380-570 GF per
  projection against our warm-started step's 22 GF (four products of 5.5 GF).
  On AMX at 300 GF/s-equivalent that is 1.3-1.9 s per iteration against 76
  ms. They pay only where a full eigendecomposition is the alternative, which
  is the refresh (6% here), and even there one refresh at 520-760 ms is
  cheaper than 1.3 s. QDWH-eig is fp64 QR/Cholesky-based (a few 2K x K
  factorisations: about 50 GF) and no cheaper. **Dead, as in
  `../../solvers/research_iteration.md` item 4, now against the 2026
  versions.**
- **Mixed precision that keeps the certificate exact.** fp32 everywhere the
  certificate does not read is already the design (fp32 eigenvectors, bf16
  tile products, fp32 R off certificate iterations). What remains in fp64 is
  read by the certificate (Zc, the transforms of a certificate iteration,
  dpotrf) or is the ADMM state (Y, Z, VRV, T). fp32 VRV/T is idea 5 above
  (-14 ms [e], not built here). **An a-priori Cholesky-free lambda bound**
  (lambda_max(V'ZcV) <= the rounding of R = W + fl64(FF') and of the
  transforms, about 1e-5 units) would remove the certificate's dpotrf at
  every check, but it certifies the code path rather than the matrix: a
  later change to how R is formed (e.g. `DNN_R32=1`) silently invalidates it,
  where Rump's test fails loudly. Rejected for robustness; the lazy
  certificate keeps Rump's test on every decision.
- **A tighter lambda in the certificate.** lambda_max(V'ZcV) = beta
  lambda_max(W - R) = 0 exactly: W - R is the negative part of W, whose
  complement is its null space. No negative lambda exists to gain
  (m+1) scale beta |theta| from, and Z -> Z + sI is an exactly flat direction
  of the Gilmore-Lawler certificate (`../rprsm/README.md`). With
  `../../roadmap/certificate/REPORT.md` (polish 1.04x, child bounds flat,
  free multipliers neutral), the dual side offers nothing but cost.
- **Structure.** The Kronecker form of L (L_(ij),(kl) = (a_ik b_jl + a_ki
  b_lj)/2) makes L recomputable on the fly (idea 9, 4 ms, killed); it does
  not shrink any product (W and R are dense, rank R* about 1,500). wil50's
  grid symmetry acts on about 1% of depth-3 nodes with a group of order 2
  (`../../roadmap/admm/projection.md`). The facial basis V is already the
  cheapest (Householder, O(m^4) per transform at 9 GB/s, memory-bound).
- **Reusing work across siblings.** Siblings of one parent share K but not
  W; the one shareable object is the parent's eigenbasis, which would save
  the child's first refresh (758 ms of a 33 s warm solve, 2.3%, 7.6 ms per
  iteration): under the line, and it needs the basis in the state files
  (5 MB fp32 per parent). Batched GEMMs across nodes need one process per
  several nodes, which the search's one-node-per-core design does not have,
  and the GEMMs are compute-bound already.
- **Memory traffic under load.** The transforms and the Y/Z-step (790 MB per
  iteration, 29%) run at 7-10 GB/s per core with three cores busy; at 32-64
  processes per socket this is the piece that shares a ceiling
  (`../../roadmap/node_speed.md` trim f). NOSYMV + NOSYMW remove 150 MB
  (idea 2); fp32 VRV/T another 115 MB (idea 5); the triangle Y/Z-step 78 MB
  (idea 6). Phase 0's loaded/alone ratio measures the ceiling; the switches
  here reduce the exposure.
- **Less frequent certificates** (every 20 instead of 10): the lazy
  certificate makes the Cholesky part conditional instead, keeping the
  10-iteration resolution of the fathom and give-up decisions; a coarser
  cadence would save the remaining 11 ms of GL per check at the price of 5
  iterations of resolution on every fathom (about 2.5% of a 200-iteration
  node), a net loss.

## 6. Validity of the switches

- **Lazy certificate.** The certified value at a check is val - margin_poly
  - sdp_penalty with sdp_penalty = (m+1) scale lambda > 0 (lambda = delta +
  c + u maxdiag, all positive). So val - margin_poly (the polyhedral value
  alone, `cert_poly`) is an upper bound on the certified value, and when
  nextafter(val - margin_poly + shift) <= cutoff no certificate at this check
  can fathom. At such checks Rump's test (`cert_lam`) is skipped; the stopping
  rules read the estimate val - margin_poly - (the last full certificate's
  penalty), which differs from the certified value by the penalty's change
  between checks (about 1e-9 units: `logs/ident/w10_*.log`), and `best`, the
  value the node returns, is updated only by full certificates: the first
  check, every check that could fathom (today's path, including the double
  re-projection when the penalty alone decides), the last iteration, and at
  the exit the check that stopped the run plus the stored multiplier with the
  highest estimate if it can raise `best` (a 39 MB copy when a lazy check's
  estimate makes a new maximum). So a fathom is decided exactly as before
  (the full certificate runs whenever one is possible), and the returned lb
  is the maximum over all checks of certified values as before, up to the
  argmax over estimates against certified values (1e-9). The give-up rule's
  inputs differ by 1e-9, so a decision can move only within 1e-9 of its
  threshold. `certify()` itself is the two parts called in the old order and
  sums the margin in the old association: bit-identical with the switch off.
- **Partial refresh.** ssyevd is ssytrd, sstedc('I') and sormtr; calling them
  with the same workspace gives the same eigenvalues and the same chosen
  columns (`bench_refresh.c`: differences 0.00e+00; solver traces identical).
  Only columns the tracker never reads are left unformed.
- **NOSYMV.** Proposition 2.2 of `docs/dnn_bound/dnn_bound.tex` holds for
  sym(Zc): <Zc, Y> = <sym(Zc), Y> for every symmetric Y. The Gilmore-Lawler
  value is evaluated row by row of the full M (both arrow entries,
  `cert_poly`), which is exactly <M, yy'> = <sym(M), yy'> for every lifted
  permutation and so a valid bound for asymmetric M; the Cholesky test runs
  on the symmetrised transform (`forward(..., force_sym)`), which is
  V'sym(Zc)V to rounding, covered by the existing transform allowance. With
  a symmetric Zc the new arrow terms equal the old ones bit for bit (d0 ==
  d0t, and x + x == 2x exactly), so the default is unchanged.

### Gates on the committed binary (`gates.sh`, `logs/gates/summary.txt`, 06:52-07:22 UTC, CPU 3)

| gate | result |
|---|---|
| `make test` (brute force, resume, parallel, coord, regression against `bench/reference/`) | PASS; the regression check exact (7 instances, tolerance 0.0%) |
| `make dnncheck` (default), `make continuationcheck`, `make shipcheck` | PASS, PASS, PASS |
| `make dnncheck` with `DNN_CORE_LAZYCERT=1`; with `DNN_CORE_REFRESH_PARTIAL=1`; with `DNN_CORE_NOSYMV=1 DNN_NOSYMW=1`; with all three | PASS x 4 (part A: max(lb - optimum) <= 0 over 2,592 node bounds; part C: the study's fathom set) |
| `tests/brute_force.py --solver-args "-x 0:20" --jobs 1` with each of the four settings exported | PASS x 4 (4,536 of 4,536 each) |
| old binary (ab937ca) against new, switches unset: nug20 root, wil50 nodes 19 and 10 | identical traces and RESULT (`logs/check_identity.txt`) |

## 7. References (checked through arXiv's API from this box, Nov 6 2026; abstracts and, for the AMX paper, the HTML text)

- arXiv 2609.04663 (Sept 2026): BF16 component-product emulation of FP32 and FP64 GEMM on Intel AMX. Three bf16 components, six products, Xeon Platinum 8462Y+, oneAPI 2026.1: 1.14-2.56x over oneMKL SGEMM by size (1.39x at n = 2048, 1.14x at 4096); the eight-tile register file named as the microkernel's limit. Our fused three-product kernel: 2.0-2.3x over OpenBLAS sgemm at K = 2117 (`../../amx/README.md`), 2.17x in the solver (section 1).
- arXiv 2505.16932v5 (2025): The Polar Express, optimal polynomial matrix-sign iterations (Muon); arXiv 2610.26742 (Oct 2026): Polar Switch; arXiv 2609.39595, 2610.34158 (2026): convergence of practical Muon with finite Newton-Schulz steps. All GEMM-only sign iterations of 5-8 degree-5 steps at 1e-2 to 1e-3 relative accuracy, more near a spectral cluster at zero; costed in section 5.
- arXiv 2507.09165v2 (Kang, Han, Groudiev, Yang, 2025): factorization-free PSD projection by composite polynomial filtering, evaluated by matrix products only (GPU ADMM); 22-31 K x K x K products at our accuracy (`../../solvers/research_iteration.md` item 4).
- arXiv 2603.10514 (Mar 2026): the condition number of Chebyshev-filtered vectors in ChASE (the one-pass CholQR question; our one pass has not failed on wil50); arXiv 1805.10121 (ChASE), 1404.4161, 1305.5120: Chebyshev subspace iteration over sequences of eigenproblems, the tracker's family.
- arXiv 2509.00495 (2025), 2211.03339 (2022), 2607.23187 (2026): mixed-precision and gap-aware Jacobi eigensolvers (idea 10, killed by cost).
- arXiv 2407.15049v3, 2403.09133v3 (LoRADS, cuLoRADS, 2024-2025): low-rank ADMM splitting for SDP on GPUs; not applicable at rank R* about 1,500 of 2,117 (`../../roadmap/admm/projection.md` idea 9).
- Rump, Verification of positive definiteness, BIT 46 (2006): the Cholesky test the certificate uses, kept on every decision by the lazy certificate.
- In-repo: `../../amx/README.md`, `../../roadmap/admm/projection.md`, `../../solvers/research_iteration.md`, `../../roadmap/node_speed.md`, `../cs2/README.md`, `../r32_noncert/README.md`, `../rprsm/README.md`, `../../roadmap/certificate/REPORT.md`, `docs/dnn_bound/dnn_bound.tex` section 2.

## 8. Files

- `profile.sh`, `profile_table.py`: the profile (section 1); `logs/prof_*.log`, `logs/profile_table.txt`.
- `scan_amx.sh`, `pair_amx.sh`, `build_bench.sh`: the kernel knobs (section 4); `logs/scan_amx.txt`, `logs/pair_amx.txt`.
- `bench_refresh.c`: the refresh microbenchmark; `logs/bench_refresh.txt`.
- `check_identity.sh`: old against new binary, switches off and on; `logs/check_identity.txt`, `logs/ident/`.
- `run_all32.sh`, `summarize.py`: the all32 A/B (section 3); `runs/`.
- `gates.sh`: the repository's checks with each switch; `gates2.sh`: the two dnnchecks re-run after the give-up-rule fix, then the all32 start; `logs/gates/` (`summary.txt`).
- Code: `latest/dnn.C` (`cert_lam`/`cert_poly` split of `certify`, `DNN_CORE_LAZYCERT`, `DNN_CORE_REFRESH_PARTIAL`, `DNN_CORE_NOSYMV`), `latest/amx_gemm.c` (`DNN_CORE_AMX_*`).
