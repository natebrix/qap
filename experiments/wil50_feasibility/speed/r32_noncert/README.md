# R's rank term in fp32 except on certificate iterations (Nov 5, cloud)

Item 1 of `../../roadmap/admm/projection.md`, row 3 of
`../../roadmap/admm/README.md`. Switch: `DNN_R32_NONCERT=1` (default off).

**Verdict: passes; adopt it as an opt-in production setting.** On the 32 wil50
nodes it costs **4.2% less per iteration** (376.2 -> 360.5 ms, measured). The
fathom set is the same: 18 nodes, 17 at the same iteration and one a check
earlier. The certificate margin is the same as the default's on all 891
certificate iterations (median 4.0e-6, max 5.1e-6). On nug30 depth-2 nodes
(m = 28) the gain is smaller: **about 2.5%** (1.7% and 3.3% in two passes,
2.2% from the form-R timers; measured), with the same outcomes and
iterations. None of the three kill criteria fires. The nug30 figure is close
to the 2% line, so on nug30 this is a small gain, not a lever.

## The idea

The subspace R-step forms R = W + FF', where F = X·diag(sqrt(-theta)) holds
the K x nneg scaled negative Ritz vectors (fp32). It uses an fp64 `dsyrk`
and a symmetrising copy. `DNN_R32=1` (the trims table in `../../amx/README.md`)
did the rank term in fp32 on every iteration: 409 -> 395 ms per iteration on
the AMX host, with the same fathoms. It was dropped because the
certificate's rigorous margin grew from 5e-6 to 0.08. With W - R =
-fl32(FF'), W - R is no longer exactly negative semidefinite, so the
Cholesky test needs a larger shift.

The certificate of iteration k uses Zc = Z + beta (Y - VRV'), formed with
iteration k's R, so V'ZcV = beta (W - R). This is "The multiplier" in
`docs/dnn_bound/dnn_bound.tex` Section 2.2. The proposition there holds for
any Z, so the R of earlier iterations only shapes the trajectory. Only the R
of a certificate iteration has to be W + fl64(FF'). Every other R is just the
ADMM step, which is already inexact at the fp32 level because the
eigenvectors are fp32. Validity never depends on R in any case: the
Cholesky test certifies whatever Zc it is given. Only the size of the
margin does.

**Change** (`latest/dnn.C`, `proj_psd_subspace`, commit 89a04c5):

- A new `cert_step` argument, passed `rstep_cert` from the main loop. It is
  true on certifying iterations: every `check_every` = 10th and the last.
  With the prototype certificate (`prototype_cert`) it is never true, since
  that certificate re-projects V'ZV in double and does not read R.
- `DNN_R32_NONCERT=1` takes the existing `DNN_R32=1` path when `!cert_step`:
  ssyrk on F in fp32, added to W in fp64, with both triangles written in one
  pass.
- Certificate iterations keep today's fp64 `dsyrk`.

The default path is untouched: with the switch unset, the new condition
reduces to the old `if (r32 && k > 0)`.

## Method

- **all32.** `run_all32.sh` runs the 32 wil50 nodes of
  `../../tuning/giveup/all32.txt`: 16 depth-3 nodes at m = 47 (lines 1-16)
  and 16 depth-4 nodes at m = 46 (lines 17-32). It is adapted from
  `../../solvers/accel/run_all32.sh`.
  - Flags: degree 3, `-L 2`, PRSM `-g 1.0 -G 1.0`, `-r 40`, `-H 300`,
    `-i 2500 -u 48816`, and `DNN_CHOLQR_PASSES=1`. This is bundle D of
    `../../amx/README.md` without AMX: this host is a 4-core 2.1 GHz Xeon
    used without AMX.
  - One process at a time on core 0: `taskset -c 0`,
    `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, `-T 1`. Three other
    experiments ran on cores 1-3 throughout.
  - **base** (default) and **nc** (`DNN_R32_NONCERT=1`) alternate node by
    node, so drifts in the host's load hit both alike.
  - Tables: `summarize.py` (`summary.txt`). Logs: `runs/`.
- **Timing.** ms per iteration is thread CPU over iterations, summed over
  nodes. The other experiments' memory traffic moves the W·Y products by up
  to ±4% between consecutive runs. Per-node ratios are therefore noisy (one
  node, 23, shows 0.2%). The `-v -v` form-R timers (F, syrk, symmetrise)
  isolate the step that changed. A one-minute compile at nice 19 shared
  core 0 at the start of `base_6` (ratio 0.977, not an outlier).
- **nug30** (m = 28). `run_nug30.sh` runs eight depth-2 nodes, cold: rows 3,
  25, 75, 100, 125, 0 and 175 of `../../../dnn_bound/predictor/nug30_depth2.csv`
  and the eigensolver study's `4:10,6:8`.
  - Flags: the record setting's DNN flags (PRSM 1/1, degree 3, `-L 2`,
    beta 0.21, refresh 120, one-pass CholQR), give-up after 200, ub 6124.
  - Two passes: base then nc in the listed order, then (`PASS=2`) nc then
    base in reverse order (`nug30/*2_*.log`).
- Extra runs through `run_all32.sh`: `CONFIGS="r32:-d_3_-L_2_-r_40:DNN_R32=1"
  NODES=19` gives the `DNN_R32=1` reference, and `CONFIGS="base2:-d_3_-L_2_-r_40:
  ncb:-d_3_-L_2_-r_40:DNN_R32_NONCERT=1" NODES="19 23 2"` gives the tiled
  variant below (built from the dropped commit).
  - Tables: `summarize_nug30.py` (`summary_nug30.txt`).

## Results

### all32 (wil50, m = 46-47; measured)

| | base | nc (`DNN_R32_NONCERT=1`) |
|---|---|---|
| fathomed | 18 | 18, the same nodes |
| iterations to fathom (18 nodes) | 4,230, median 220 | 4,220, median 220: 17 nodes the same; node 22 at 390 instead of 400 |
| **ms per iteration (all 32)** | **376.2** | **360.5 (-4.2%)** |
| m = 47 nodes / m = 46 nodes | 403.0 / 337.7 | 385.9 / 323.8 (-4.2% / -4.1%) |
| per-node ratio nc/base | | median 0.957, range 0.929-0.998 |
| CPU, fathoms / misses / all | 1,468 / 1,888 / 3,356 s | 1,404 / 1,808 / 3,212 s (-4.3%) |
| form R per iteration: F + copy, syrk, symmetrise | 6.6 + 29.6 + 12.0 = 48.2 ms | 1.2 + 17.3 + 14.3 = 32.8 ms (-15.4) |
| misses' lb at 200, nc - base (14 nodes) | | median +0.005, range -0.053 to +0.081 |
| misses' lb at 300, nc - base (14 nodes) | | median -0.011, range -0.056 to +0.069 |
| certificate margin, all certificate iterations | median 4.03e-6, max 5.10e-6 | identical on all 891 paired certificate iterations |
| certificate failures, refinements | 0, 0 | 0, 0 |

- The whole-iteration saving (15.7 ms) matches the form-R timers (15.4 ms).
  The rest of the iteration is unchanged on average: W·Y is 150.6 against
  150.7 ms per step at m = 47 and 120.4 against 120.7 at m = 46.
- The fp32 iterations save about 17 ms of form R:
  - ssyrk 15-16 ms against dsyrk 29-33.
  - No 36 MB copy of W.
  - The fused add-and-mirror pass is about 2.5 ms slower than the fp64
    mirror copy.
- Node 22 fathoms one check earlier because base's bound at iteration 390
  was 48,814.046, 0.004 below the cutoff 48,814.05, and nc's was 48,814.059.
  By iteration 200 the two trajectories' bounds differ by up to ±0.08 (the
  misses), so this is a coin flip at the cutoff, not a systematic gain.
  The fathom set does not change, and no fathom comes later.
- For reference, `DNN_R32=1` (fp32 on every iteration) on node 19
  (`runs/r32_19.log`) fathoms at the same iteration 130 as base. Its margin
  grows from 3e-4 at iteration 10 to 0.077 at 130, against 4.4e-6 for both
  base and nc. This reproduces the amx study's 0.08. The switch tested here
  removes that cost completely.

### nug30 (m = 28; measured)

| node | outcome, iterations (base = nc) | lb nc - base | form R base / nc (ms per iteration) |
|---|---|---|---|
| 0:0,23:1 | give-up, 280 | +0.0014 | 2.8 / 2.1 |
| 0:0,23:26 | fathom, 100 | -0.0004 | 2.5 / 2.1 |
| 0:0,23:4 | fathom, 120 | -0.0005 | 2.6 / 2.2 |
| 0:12,23:1 | give-up, 230 | -0.0004 | 2.7 / 2.1 |
| 0:2,23:18 | fathom, 180 | +0.0069 | 2.7 / 2.2 |
| 0:6,23:14 | fathom, 140 | +0.0008 | 2.7 / 2.0 |
| 0:7,23:10 | fathom, 210 | -0.0003 | 2.6 / 2.1 |
| 4:10,6:8 | fathom, 330 | -0.0017 | 2.7 / 2.3 |

- ms per iteration over the 1,590 iterations:
  - pass 1: 23.06 -> 22.67 (-1.7%)
  - pass 2 (reverse order, nc first): 23.03 -> 22.26 (-3.3%)
- The form-R timers give 2.65 -> 2.15 ms (-0.5 ms, -2.2%). At K = 730 and
  nneg about 190, the rank term is 11% of the iteration, against 13% on
  wil50.
- Every outcome and iteration count is the same, and the bounds agree to
  0.007.
- Margins are identical: median 1.64e-7, max 1.80e-7.

### Validity (on the final code)

- `make dnncheck`: OK with the switch off and with `DNN_R32_NONCERT=1`
  exported (`checks/dnncheck_*.log`). Part A gives max(lb - optimum) =
  -1.2e-13 on the m >= 2 nodes.
- `python3 tests/brute_force.py --solver-args "-x 0:20" --jobs 1` with
  `DNN_R32_NONCERT=1` exported: 4536 of 4536 passed
  (`checks/brute_force_r32nc.log`).
- **The default path is unchanged.** With the switch off, wil50 node 19 (130
  iterations, 4 full eigendecompositions, 6 block growths) was run with
  the binary from before the change and with the new one, and so was nug30
  `4:10,6:8` (330 iterations). The `-v` logs are identical apart from CPU
  times: same iterations, bounds, margins and RESULT (`checks/ident_*.log`).
  The 3-node rerun below reproduced base's traces exactly.

### Tried and not kept: a tiled add-and-mirror pass

A 64 x 64 tiled version of the fp32 path's add-and-mirror loop was
measured on nodes 2, 19 and 23 (`runs/base2_*`, `runs/ncb_*`,
`summary_tiled.txt`). It gives the same R bit for bit: the traces are
identical to nc's. But its loop timer (mean over all iterations) read
16.0-17.3 ms against 12.1-16.4 ms for the column loop on the same nodes,
and the whole-iteration ratio (0.955 against base2) was inside the noise.
The column loop stays, and the commit was dropped.

## Kill criteria (task and memo)

| criterion | result |
|---|---|
| margin back at about 5e-6 on certificate iterations | **yes**: identical to base on every certificate iteration (max 5.1e-6 wil50, 1.8e-7 nug30) |
| any fathom changes on all32 | **no**: the same 18 nodes; one fathom a check earlier (above) |
| gain under 2% per iteration | **no**: 4.2% on wil50; about 2.5% on nug30 (1.7-3.3%) |
| memo: bound at 300 moves by more than 0.1 | **no**: within ±0.07 |

## Using it

Set `DNN_R32_NONCERT=1` in the environment of `dnntest` or `bnbdrv -x` (it
is read by `dnn.C`, like `DNN_CHOLQR_PASSES`). It stacks with the other
subspace settings. It should also stack with AMX, but that was not tested
here. On the AMX host the memo expects -12 to -13 ms of 330-409 (about 4%) [estimate]. The
fp32 path's add-and-mirror pass (about 14.5 ms at K = 2117 here) is now the
largest piece of form R after the ssyrk. A faster transpose-write there is
the next 2-3 ms [estimate].
