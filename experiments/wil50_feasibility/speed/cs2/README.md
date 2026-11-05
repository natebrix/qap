# Complementary-subspace Rayleigh-Ritz and a per-column filter degree: both fail their gates

Item 5 of `../../roadmap/admm/README.md` (ideas 2 and 3 of
`../../roadmap/admm/projection.md`), Nov 5 2026. **[m]** = measured here,
**[e]** = estimate. Host: 4-core 2.1 GHz Xeon (KVM), everything on core 3
(`taskset -c 3`), one BLAS thread, scipy-openblas32, the sgemm path (no
`DNN_AMX`); the other three cores ran other jobs throughout, so these are
loaded figures. Production-like flags (bundle D of `../../amx/README.md`
without AMX): `-d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i 2500 -u 48816`,
`DNN_CHOLQR_PASSES=1`.

**Verdict: both dead.**

- **CS2 fails the memo's own microbenchmark gate**: it saves **7.0 ms per
  subspace step** [m] on a random sample of 10 steps of a wil50 m = 46 node
  (range 2.5 to 11.9), against the 10 ms kill line and the memo's 15-30 ms
  estimate. That is about 2% of the 356 ms iteration [m]. On nug30 (m = 28)
  it saves 0.6 ms of 27 [m].
- **The per-column degree costs iterations.** With the memo's rule it saves
  11 ms per iteration [m], but node 19 fathoms 10 iterations later (+4% CPU)
  and node 25 30 later (+16% CPU). Its bound trails the baseline all the
  way: by about 475 at iteration 50, 155-171 at 100, and 82-97 at 130-140
  [m].
  - The accurate rule (ρ = 3, below) keeps the single-step accuracy but
    removes only 2% of the filter's products, about 2 ms [m: column counts;
    e: ms].
  - ρ = 1.5 removes 5-6%, about 6 ms [e], for up to +10% error [m].

So the all32 A/B and the nug30 search check were not run: the
microbenchmark gate is the stop. Neither change can be made to pay at this
spectrum. The reason is the same for both, below.

## What the spectrum looks like (dumped steps)

`dnntest` with `DNN_DUMP_RR=prefix DNN_DUMP_RR_AT=s1,s2,...` (a diagnostic, off by default,
in `dnn_cheb_col.patch`) writes a subspace step's inputs: W (fp32), the
block X (the previous step's Ritz vectors), their Ritz values and the filter
interval. `rr_bench.C` replays the step from them. Dumps: wil50 node 22
(`26:27,8:38,34:28,37:35`, steps 60/120/180/240/300), node 29
(`20:12,40:28,17:25,43:29`, 10 steps drawn at random from 1-330) and nug30
`4:10,6:8` (8 random steps). Late in a wil50 m = 46 run [m]:

- The block has kb ≈ 640 columns with about 540 negative Ritz values, from
  θ_min ≈ −3.5 up to just below 0. The kb − nneg ≈ 100 positive Ritz values
  all lie in (0, θ_top] with **θ_top = 3e-4 to 1.4e-3**. W itself has
  **130-520 eigenvalues within 1e-3 of zero**.
- The filter interval is [θ_top, λ_hi] with λ_hi ≈ 4.5-5, so c ≈ h ≈ 2.2-2.5.
  Then t_i = (θ_i − c)/h lies in [−1.5, −1] for all but 15-37 of the ~540
  negative columns.

## 1. CS2: the partial eigensolve is not cheaper than ssyevd

Implemented in `rr_bench.C` (`rr_cs2`, `factor_F`) exactly as
`projection.md` §2 specifies. After G = Y'WY: `ssyevr` for θ > 0 only, then
S = U₊Θ₊U₊' − G, then `spstrf` (LAPACK default tolerance), then
F = Y(:,piv)·L as `strmm` + `sgemm` on the trapezoid, then F to fp64. The
dsyrk that follows is common to both paths, at rank r instead of nneg. Today's
path is G, then `ssyevd(kb)`, then X = YG, then the fp64 scaling.

Node 29, the 10 random steps, median of 5 repetitions, ms [m]:

| step | kb | nneg | m (θ > 0) | today: ssyevd + X=YG + F | CS2: ssyevr + S + spstrf + F=YPL + F64 | saved |
|---|---|---|---|---|---|---|
| 39 | 373 | 324 | 49 | 6.8 + 3.5 + 0.7 = 11.0 | 5.0 + 0.1 + 0.5 + 2.2 + 0.7 = 8.5 | 2.5 |
| 45 | 532 | 437 | 95 | 15.7 + 8.1 + 1.0 = 24.8 | 15.0 + 0.4 + 1.2 + 4.7 + 1.0 = 22.2 | 2.6 |
| 148 | 641 | 539 | 102 | 23.2 + 11.3 + 1.2 = 35.7 | 18.3 + 0.5 + 1.7 + 6.6 + 1.2 = 28.4 | 7.3 |
| 177 | 646 | 544 | 102 | 34.5 | 27.6 | 6.9 |
| 178 | 645 | 545 | 100 | 34.1 | 27.4 | 6.7 |
| 202 | 663 | 579 | 84 | 38.4 | 28.6 | 9.8 |
| 262 | 638 | 539 | 99 | 35.0 | 27.7 | 7.4 |
| 281 | 651 | 594 | 57 | 35.8 | 23.9 | 11.9 |
| 307 | 635 | 537 | 98 | 33.8 | 27.3 | 6.5 |
| 311 | 635 | 536 | 100 | 34.6 | 26.5 | 8.2 |
| **mean** | | | | | | **7.0** |

On node 22 the five fixed steps saved 5.9-7.0 / 1.0-1.4 / 7.1-8.3 /
−6.9 to −4.7 / 7.5-10.1 ms (60 / 120 / 180 / 240 / 300, ranges over three
benchmark runs). Steps 120 and 240 follow a refresh, with kb ≈ 790 and
m = 190-224. Repeated runs of the same dump vary by about ±1.5 ms. On nug30
(kb 91-276, m 27-37) the saving is −0.1 to 1.0 ms per step, against 26.9 ms
per iteration [m].

**Accuracy.** ‖R − R_exact‖_F / ‖P₋(W)‖_F was measured against an fp64
`dsyevd` of the dumped W. CS2 equals today's path to 3 digits on every dump,
for example 1.401e-3 vs 1.401e-3 and 2.547e-3 vs 2.547e-3. The largest gap
was 4.415e-4 vs 4.364e-4 at step 281, where spstrf's default tolerance
dropped 63 negligible directions. So CS2 is a valid drop-in, just not a
faster one.

**Why the memo's estimate was off.** At kb ≈ 640 [m]:

- `ssytrd` alone costs 7.5-9 ms.
- Every way of getting the ~100 positive eigenvectors of the reduced matrix
  costs 17-20 ms, against ssyevd's 21-25 ms for all of them:
  - `ssyevr` (MRRR): 17-19 ms.
  - `ssyevx` (bisection and inverse iteration): 18-20 ms.
  - `ssytrd`, then `sstedc` for all tridiagonal vectors, then `sormtr` on the
    positive ones only: 17-19 ms.
- MRRR for all pairs takes 36-48 ms, twice ssyevd.

The positive pairs sit in the dense cluster at zero, and resolving
eigenvectors inside a cluster is what costs. The saving left is about 5 ms in
the eigensolve and 5 ms from the trapezoidal F instead of X = YG, minus
2.3 ms for S and spstrf.

The same cluster rules out **CS2's second stage** (an inner Chebyshev filter
on H for its top pairs, instead of ssytrd). Separating (0, 3e-4] from the
negative near-zero eigenvalues, against a spread of 4, needs a degree in the
hundreds [e].

The mean of 7.0 ms is also an upper bound for the solver. CS2 filters Y
without the Ritz rotation, so the block cannot be trimmed to nneg + guard on
CS2 steps. It would carry up to 50 extra columns between full Rayleigh-Ritz
steps, at about 0.2 ms per column in the filter [e].

**Shift variant, no eigensolve at all.** S = θ_top·I − G, then Cholesky.
This is not CS2: it adds a PSD bias of norm θ_top to R. It saves 16-28 ms
per step [m], but its R error is 32-49% above today's [m]. That is worse than
uniform degree 2. Dead.

## 2. Per-column degree: the strongly negative columns are too few

The memo's rule (projection.md §3) is d_i = min{d : |C_d(t_i)| ≥ |C_3(t_e)|},
with t_e at the largest negative Ritz value. Because θ_top ≈ 0, t_e ≈ −1.0001
and |C_3(t_e)| ≈ 1. The rule therefore gives degree 1 to almost every negative
column. It "passes" the memo's histogram gate (70-84% of columns drop a
degree; products 1.35-2.42 kb instead of 3 kb), but the drop costs accuracy.

`rr_bench -x`, single-step R error relative to uniform degree 3, on the
wil50 dumps [m]:

| rule | filter products | R error vs degree 3 |
|---|---|---|
| uniform 4 | 4.00 kb | −4 to −15% |
| uniform 2 (known to cost iterations, `projection.md` table row 5) | 2.00 kb | +4 to +19% |
| memo rule | 1.35-2.42 kb | **+7 to +36%** |
| amplification floor ρ = 1.5 (d_i = min{d : \|C_d(t_i)\| ≥ 1.5}) | 2.71-2.86 kb | +0.7 to +10% |
| floor ρ = 3 | 2.91-2.95 kb | +0.0 to +0.7% |
| floor ρ = 1.5 up to degree 4 (node 22 dumps) | 3.35-3.69 kb | −1 to −4% |

The memo rule is worse than uniform degree 2, and every rule trades
accuracy for products at the rate of a uniform degree change or worse.
Only 15-37 negative Ritz values have |t| ≥ 1.5. The bulk sits in
t ∈ [−1.5, −1], where C_1, C_2 and C_3 differ the most in relative terms.
On nug30 the memo rule is +42 to +59% on the six late steps (+12 to +15% on
the two early ones), against +15 to +25% for degree 2.

**In the solver** (`DNN_CHEB_COL=1`, opt-in, in `dnn_cheb_col.patch`), two m = 46
nodes against the baseline, same binary [m]:

| node | fathom it | CPU s | ms/it | W*Y ms/step | lb − baseline at 50 / 100 / 130-140 |
|---|---|---|---|---|---|
| 19 (`26:27,8:38,34:28,37:21`) base | 130 | 41.9 | 322.4 | 110.7 | |
| 19 memo rule | 140 | 43.6 | 311.6 | 94.8 | −475 / −171 / −97 (130) |
| 25 (`20:12,40:28,17:25,43:18`) base | 150 | 48.6 | 324.0 | 116.6 | |
| 25 memo rule | 180 | 56.3 | 312.8 | 92.8 | −473 / −155 / −82 (140) |

So the rule saves 11 ms per iteration and loses on both nodes: by the kill
criterion, it fathoms later and needs more iterations. Error-neutral rules
(ρ ≥ 3) save about 2% of the products, about 2 ms [e].

## Switches (all off by default)

These are not in the branch's `latest/dnn.C`. Both ideas are dead, so the
solver change is kept as `dnn_cheb_col.patch`, the diff against a955e95
(`git apply` it to a955e95's `latest/dnn.C` to rebuild the dump or the
per-column filter).

| switch | what |
|---|---|
| `DNN_DUMP_RR=prefix`, `DNN_DUMP_RR_AT=s1,s2,...` | write the subspace step's inputs at those step counts to prefix_<step>.bin (diagnostic; reads nothing back) |
| `DNN_CHEB_COL=1` | per-column filter degree, the memo's rule; columns scaled to unit norm before CholQR |
| `DNN_CHEB_COL_RHO=x` | with `DNN_CHEB_COL=1`: the amplification-floor rule, target x instead of \|C_D(t_e)\| |

CS2 lives only in `rr_bench.C`. It was not put in the solver because the gate
failed.

## Validity and the default path

- **Default path unchanged.** `identity_check.sh` runs node 21
  (`26:27,8:38,34:28,37:2`, 150 iterations, production-like flags) with the
  binary built from `origin/cognac` a955e95 before any change and with the
  new one (switches off). All 16 lines (15 trace lines and RESULT), without
  the CPU fields, are byte-identical (lb = 48777.317240424942) [m]
  (`logs/ident_*.trace`). `make dnncheck` with default settings passes
  (`logs/dnncheck_default.log`). Node 29, run with the dump on, fathoms at
  340, as bundle D did on the AMX host.
- **The certificate never sees the filter.** R = W + FF' with F from the
  computed vectors in both new paths (CS2 in the bench, the per-column filter
  in the solver), so W − R is negative semidefinite by construction.
  - `make dnncheck` (tests/dnn_test.py) with `DNN_CHEB_COL=1`
    (`logs/dnncheck_chebcol.log`; dnn_test uses degree 2, so d_i is 1 or 2):
    validity holds. Part A gives
    max(lb − optimum) = 0 over 2,592 node bounds, and part B agrees with the
    prototype to 8.6e-8. Part C fails on strength only: nug20 `[1:6,19:8]`,
    which the default fathoms at iteration 330, gives up at 400. Its other
    four fathoms all come 10-20 iterations later than the default's. This is
    the same verdict as on wil50.
  - `tests/brute_force.py --solver-args "-x 0:20" --jobs 1` with
    `DNN_CHEB_COL=1` and the production DNN settings (`QAP_DNN_CHEB=3
    QAP_DNN_TIGHT_LAM=2 QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0
    QAP_DNN_REFRESH=40 DNN_CHOLQR_PASSES=1`) on core 3: **4536/4536 passed**
    (`logs/bruteforce_chebcol.log`).

## Files

- `rr_bench.C`, `build_bench.sh`: the microbenchmark. Usage:
  `taskset -c 3 ./rr_bench [-r reps] [-x] dump.bin ...`. With `-x` it also
  runs the exact comparison and the per-column rules.
- `run_nodes.sh`: sequential A/B on all32 nodes, pinned. It is
  `../../solvers/accel/run_all32.sh` with `-P 1` and `taskset -c 3`.
- `identity_check.sh`: the default-path identity check.
- `dnn_cheb_col.patch`: the solver change (the dump and `DNN_CHEB_COL`),
  against a955e95.
- `logs/`:
  - `bench_*.txt`: the benchmark outputs.
  - `*_dump.log`: the dumping runs, with step splits.
  - `ident_*`: the identity check.
  - `dnncheck_*.log`, `bruteforce_chebcol.log`: the validity checks.
- `runs/`: node 19 and node 25, baseline and memo rule.
