# Two-phase ADMM penalty: beta high early on cold nodes, the tuned value from the first refresh after iteration 60

Item 2 of `../../roadmap/admm/splitting.md` (row 4 of `../../roadmap/admm/README.md`):
a penalty beta that changes only at full eigendecomposition refreshes and is
frozen after iteration 200. The memo's first test is its "10-line" two-phase
check: beta = 0.5 m/3 until the first refresh after iteration 60, then the
tuned 0.3 m/3, with 0.2 then 0.3 as the control. Both that schedule and the
memo's residual-balancing controller are built, behind one opt-in switch.
Figures marked [m] are measured on this host (4-core 2.1 GHz Xeon without AMX,
one process at a time on core 1, one BLAS thread); [e] are estimates.

## Verdict

**The two-phase schedule passes the kill test on cold nodes by a wide margin,
and must not be applied to warm-started children.**

- all32 (cold, production flags), `DNN_BETA_SCHED=two,f=1.6667,t=60` against
  the default [m]: the same 18 fathoms (none lost, none gained); **every one
  of the 18 fathoms 20-50 iterations earlier** (median 30; 4,230 -> 3,690
  iterations, -13%; CPU on the fathoms 1,492 -> 1,238 s, -17%); on the 14
  misses **lb@100 +615** (median; range +237 to +2,575) and **lb@200 +45**
  (range +31 to +83). The kill criterion (lb@200 up by under 20 and fathoms
  under 10 iterations earlier, or a lost fathom) is missed on every count.
  Total CPU over the 32 nodes 3,387 -> 2,960 s (-13%).
- nug30 depth-2 nodes, the record's setting (base 0.21, `-r 120`; the switch
  is then at iteration 122, and 0.35 m/3 before it) [m]: on 28 depth-2 nodes
  the same 22 fathoms, **each 20-50 iterations earlier** (median 35; 3,190 ->
  2,430 iterations, -24%; 66 -> 50 s), the 6 misses +3 to +12 at iteration
  200; on 8 more nodes (with `4:10,6:8`) 1,160 -> 900 iterations on the 6
  fathoms. The schedule does not need retuning for m = 28 and beta 0.21.
- Warm-started children (`-Q b`, the search's `QAP_DNN_WARM_Y=2`) [m]: the
  same schedule **hurts**: 8 children of two wil50 parents (all32 nodes 9 and
  6) need 620 iterations instead of 450 (+38%; slower on every child). They
  fathom in 30-130 iterations, almost entirely inside the high-beta phase.
  The opposite helps them: beta 0.21 m/3 on warm calls gives 380 (-16%;
  faster or equal on every child).
- So the setting to use is **`DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7`**:
  cold calls 0.5 -> 0.3 m/3, warm-started calls 0.21 m/3 throughout. On the
  node-9 family (parent cold, four warm children) it takes 630 -> 540
  iterations and 236 -> 203 s (-14%), on the node-6 family 620 -> 540 and
  225 -> 201 s (-11%) [m]. `nodes=cold` (warm calls at the default) is the
  conservative alternative: 630 -> 570 iterations on the node-9 family, the
  children unchanged.
- In the search most DNN nodes are warm-started (the nug30 record ships the
  parent's state to every child), so the cold-node gain (13-24% of the
  iterations of a cold node) reaches only the first DNN level and nodes whose
  parent left no state [e]. The warm factor (-16% on 8 children of two
  parents) is the part that would matter for the whole tree, and it is the
  least measured: it needs a search A/B before any default changes.
- The residual-balancing controller is built (`rb`) but not worth running on
  wil50: the normalized ratio it balances does not separate the transient
  from the tail on nug20 (median rhat/shat 21 at iterations 10-40 and 20 at
  100-200 at the default beta) nor on kra30a (29 and 24) [m], and on the wil50
  children it rises after the switch rather than falling. The memo's
  robustness test fails: on the kra30a root no setting moves the bound off
  about 17,365 at 600 iterations (QPB is about 68,600); `rb,c=20` never
  changes beta there, and the two-phase schedule adds only +0.6. The
  controller is closed; the two-phase schedule is the result.
- Fixed beta 0.5 m/3 (no switch back) is worse than the default on wil50:
  fathoms 40 iterations later on the 8-node subset, lb@200 on the misses only
  +16 [m]. The gain is the schedule, not a larger beta.

## The switch

`DNN_BETA_SCHED` (environment, read by `dnn_bound` at every call, so it acts
in `dnntest` and in `bnbdrv -x` alike). Unset or empty: nothing changes.

| value | what it does |
|---|---|
| `two,f=F,t=T` | beta = F x base until the first **scheduled** full refresh at or after iteration T, then base (base = `par.beta` m/3: 0.3 m/3, or `QAP_DNN_BETA`) |
| `...,f2=F2` | the second phase at F2 x base instead of base |
| `...,force=1` | change at iteration T itself and make that iteration's R-step a full refresh (one extra ssyevd) |
| `...,w=FW` | warm-started calls (a `DNNWarm` start with Z) run at FW x base throughout instead; `w=FW` alone leaves cold calls at base |
| `rb,c=C[,every=40][,max=4][,lo=0.25][,hi=4]` | the memo's residual balancing (below) |
| `...,freeze=N` | no change after iteration N (200) |
| `...,nodes=cold\|warm` | only calls without / with a warm start (default: all) |

Measured settings: **`two,f=1.6667,t=60`** (0.5 m/3 until the first
scheduled refresh at or after iteration 60, then 0.3; with `-r 40` refreshes
are at 1, 42, 83, 124, ..., so the switch is at 83; with nug30's `-r 120` it is
at 122), and **`two,f=1.6667,t=60,w=0.7`** (the same on cold calls, 0.21 m/3
on warm ones). f and w are relative to the base, so on nug30 (base 0.21)
f = 1.6667 gives 0.35 m/3 and w = 0.7 gives 0.147 (untested there).

`rb`: at every check (10 iterations), S += log(rhat/shat) - log C with
rhat = |Y - VRV'|/max(|Y|, |VRV'|) and shat = beta |dY|/max(|L/scale|, |Z|)
(Frobenius norms of the lifted matrices, no V projection: three extra passes
over n1^2 doubles per check). At a scheduled refresh at least `every`
iterations after the last change, with fewer than `max` changes made,
f = clip(exp(S/2n), 0.5, 2); if f is outside [1/1.5, 1.5], beta *= f (kept in
[lo, hi] x base) and S = n = 0. With any schedule on, `-v` prints a `bsched`
line at every check (beta, pres/(beta dY), rhat, shat, the norms) and one at
every change; `two,f=1,t=100000` is a diagnostic-only setting (beta never
changes; the trajectory is identical to the default, checked on nug20).

Implementation (`latest/dnn.C`, struct `BetaSched` and the top of the main
loop): beta and 1/beta become variables; a change happens before the
iteration's forward transform, only when `w.since_refresh >= refresh` (the
subspace tracker is about to run its scheduled ssyevd; in the exact
eigensolver modes every iteration qualifies), so W, the R-step, the
certificate multiplier, the Y-step and both Z-steps of that iteration all use
the new beta. `DNN_FUSE`'s precomputed rows (formed with the old 1/beta) are
invalidated. The schedule is off with Anderson, Halpern or momentum (Anderson
stores Z/beta). Cost: none per iteration; three norms per check.

### Rescaling: none is needed, checked by hand

Z in `dnn.C` is the unscaled multiplier: the Y-step is
Pi_P(VRV' - (L/scale + Z)/beta), the R-step P_psd(V'(Y + Z/beta)V), the Z-steps
Z += gamma beta (Y - VRV'). The KKT points do not depend on beta, so a change of
beta keeps Y and Z as they are (a scaled dual U = Z/beta would need
U *= beta_old/beta_new). A child's warm start maps Z entrywise
(`dnn_warm_child`) and rescales it only by the objective scale, never by beta,
so it is valid for whatever beta the child uses (the `w` key relies on this).
The certificate (Section 2 of `docs/dnn_bound/dnn_bound.tex`) is valid for any
symmetric multiplier; its tightness needs V'Z'V = beta (W - R) <= 0, which
holds because W and Z' use the same beta. A stale W would show as a failed
Cholesky test or a large margin at the first check after a change: on the
wil50 runs every check passed Cholesky (cert_fail = 0, no refinements) and
the margins stay at 3-5e-6, as in the default runs.

`beta_change_check.C` (output `beta_change_check.out`) checks this on 15
random QAPs of size 4-8 with brute-forced optima, PRSM, in all three
eigensolver modes [m]:

- **scheduled**: `two,f=2,t=21` (beta 0.6 m/3 for iterations 1-20, then 0.3), 45 iterations;
- **split**: 20 iterations at 0.6 m/3, saved, then a same-node restart
  (`warm_primal`, Y and Z as saved) at 0.3 m/3 for 25;
- **scaled** (negative control): the same split with Z multiplied by 0.3/0.6
  at the restart, as a scaled dual would require.

With the exact eigensolvers (dsyevd, ssyevd) scheduled and split agree
**exactly** (max difference in Y and Z: 0; in the best certified bound: 0) in
all 30 cases; with subspace tracking and `force=1` they agree exactly in 10 of
15 and within 0.004 in the others (the tracker's block-growth seed differs).
The scaled control differs by 0.16-0.48 in Y or Z in every case (its bound is
lower in 39 of 45), so the test would see a wrong rescale. Every bound is at
or below the optimum. PASS.

### The default path is unchanged

- nug20 root, 400 iterations, production flags: the binary before the change
  and after it (both commits) with the switch unset print identical lines at
  all 40 checks (lb, certificate, margin, primal, pres, dY, npos) and the same
  RESULT lb to 17 digits (`nug20/oldbin.log`, `nug20/base.log`,
  `nug20/v2base.log`) [m]. The second commit leaves cold calls of a schedule
  unchanged too (`nug20/v2two05w.log` = `nug20/two05_60.log`).
- wil50 node 17 (all32), 210 iterations to fathom: identical at all 21 checks
  and in the RESULT line (`runs/base_17.log` from the old binary,
  `runs/basenew_17.log` from the new) [m].
- The all32 baseline reproduces the fathom set (18) and the iteration counts
  (4,230 in total, median 220) of `../../solvers/accel/runs/amxE_*` exactly.

### Validity (`validate.sh`, `validate_warm.sh`, logs in `validate/`) [m]

- `make dnncheck`, switch unset: OK (before and after the `w` key).
- `tests/dnn_test.py --parts ACDE` with the switch on, for `two,f=1.6667,t=60`,
  `two,f=1.6667,t=5,force=1`, `two,f=0.6667,t=5,force=1` and `rb,c=8,every=10`
  (the last three change beta on every node of more than 5 or 10
  iterations): OK in all four (max(lb - optimum) = 0 over 2,592 brute-forced
  node bounds; part C fathoms exactly the study's nug20 nodes, 0-20
  iterations sooner). Part B compares `dnntest -P` with the fixed-beta
  prototype `dnn.py`, so it does not apply.
- `tests/brute_force.py --solver-args "-x 0:20" --jobs 1` with the switch on,
  the same four settings: 4536/4536 passed each.
- `tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1` with
  `QAP_DNN_WARM_Y=2` (children of DNN misses start warm, as in the nug30
  record), switch unset and with `two,f=1.6667,t=60,w=0.7`,
  `two,f=1.6667,t=5,force=1,w=0.7` and `w=1.6667`: 4536/4536 passed each.

## Results

### nug20 root (`nug20_curves.txt`) [m]

Cold, production flags, no cutoff. Certified lb:

| it | default 0.3 | 0.5 -> 0.3 at 83 | 0.2 -> 0.3 at 83 (control) | fixed 0.5 | fixed 0.6 |
|---|---|---|---|---|---|
| 10 | 1421.4 | 1599.2 | 1284.0 | 1599.2 | 1723.0 |
| 50 | 2417.9 | 2470.3 | 2245.1 | 2470.3 | 2475.8 |
| 90 | 2490.0 | 2498.5 | 2463.1 | 2497.9 | 2498.4 |
| 150 | 2502.87 | 2505.01 | 2499.86 | 2504.68 | 2504.80 |
| 200 | 2504.72 | 2505.58 | 2503.96 | 2505.68 | 2505.78 |
| 300 | 2505.77 | 2505.98 | 2505.59 | 2506.07 | 2506.12 |
| 400 | 2506.08 | 2506.15 | 2506.02 | 2506.22 | 2506.23 |

The predicted direction holds: higher beta early, lower later; the control is
behind everywhere. The schedule reaches 2505.0 at iteration 150 against about
220 for the default. Unlike wil50, a fixed higher beta is also good on the
nug20 root (beta 0.3 was tuned with gamma = 1.618, before PRSM), and the
switch back costs it 0.07 at iteration 400; on wil50 the fixed 0.5 loses (below).

### all32 (`run_all32.sh`, `summarize.py`, `all32_summary.txt`, logs in `runs/`) [m]

32 wil50 nodes of `../../tuning/giveup/all32.txt`, cold, `-d 3 -L 2 -g 1.0 -G 1.0
-r 40 -H 300 -i 2500 -u 48816`, `DNN_CHOLQR_PASSES=1`, one at a time on core 1.

| config | fathomed | iterations to fathom: median, total | fathoms earlier by >= 10 / later | lb@100 - default, misses: median [min, max] | lb@200 - default, misses | ms/iteration | CPU: fathoms, misses, total |
|---|---|---|---|---|---|---|---|
| default | 18 | 220, 4,230 | | | | 380 | 1,492, 1,895, 3,387 s |
| `two,f=1.6667,t=60` | 18 (same set) | 190, 3,690 | 18 / 0 (each by 20-50, median 30) | +615 [+237, +2,575] | **+45 [+31, +83]** | 365 | 1,238, 1,722, 2,960 s |

Per iteration the cost is not higher: on the nine misses that both run to the
300-iteration give-up, 404 ms (default) against 390 ms (the overall 380 -> 365 is mostly
the mix of nodes). On the 300-iteration misses the bound at the give-up is 14
to 21 higher. Five misses (6, 9, 14, 15, 16) give up 40-100 iterations
earlier under the schedule, at a bound within -8.5 to +2.3 of the default's
later give-up (at equal iterations, 300-370, the schedule is 9-16 ahead): the give-up
rule extrapolates from recent gains, which are smaller once the transient is
over sooner. That saves miss CPU (1,895 -> 1,722 s) but hands the children a
parent state from an earlier iteration; the warm-children runs below include
that effect (with `nodes=cold` the children's iterations are unchanged).

Subset control (`subset_fix05_summary.txt`; nodes 1, 2, 9, 14 missed, 11, 17,
22, 31 fathomed): **fixed 0.5 m/3** fathoms all four 20-140 iterations
**later** than the default (median +40; 1,140 -> 1,380 iterations) and is
only +16 [-8, +27] at iteration 200 on the misses, against -35 and +47 for the
two-phase schedule. The memo's 0.2 -> 0.3 control was run on nug20 only
(behind at every checkpoint, above).

### nug30 (`run_nug30.sh`, `summarize_n.py`, logs in `nug30/`) [m]

Cold, the nug30 record's DNN setting: `-d 3 -L 2 -g 1.0 -G 1.0 -B 0.21 -r 120
-H 300 -i 1000 -u 6124`, `DNN_CHOLQR_PASSES=1`. m = 28, 22 ms per iteration.
`two05` is `two,f=1.6667,t=60` (0.35 m/3, switch at the refresh at 122);
`fix05` 0.35 m/3 throughout; `two05f` the switch forced at iteration 60.

| node set | config | fathomed | iterations on the common fathoms | fathoms earlier / later by >= 10 | misses: lb@100, lb@200, final lb - default | CPU, all nodes |
|---|---|---|---|---|---|---|
| 28 depth-2 nodes (`nug30_nodes_big.txt`, `nug30_big_summary.txt`) | default | 22 | 3,190 | | | 107.0 s |
| | `two05` | 22 (same) | **2,430 (-24%)** | 22 / 0 (by 20-50, median 35) | +18 to +80, **+3 to +12**, +0.6 to +4.3 | 88.6 s |
| | `fix05` | 22 (same) | 2,600 | 21 / 1 (+60) | +18 to +80, -2 to +7, -3.6 to +1.8 | 95.3 s |
| 8 nodes (`nug30_nodes.txt`, `nug30_small_summary.txt`) | default | 6 | 1,160 | | | 40.4 s |
| | `two05` | 6 (same) | 900 | 6 / 0 (by 30-80) | lb@200 +3.5, +3.8 | 34.0 s |
| | `two05f` | 6 (same) | 930 | 6 / 0 | lb@200 +4.7, +5.2 | 35.1 s |

As on wil50, the fixed higher beta gives most of the early gain but loses at
iteration 200 on the misses and one fathom comes 60 iterations later; the
switch back keeps both. Forcing the switch at 60 instead of waiting for the
refresh at 122 is about the same.

### Warm-started children (`run_warm.sh`, logs in `warm/`) [m]

Parent all32 node 9 (`26:27,8:38,34:28`, a miss) runs as in all32; then
children 37:14, 37:37, 37:21, 37:26 (all32 nodes 17-20) start from its final
state with `-Q b` (Z and the conditioned Y, as `QAP_DNN_WARM_Y=2` does in the
search).

| config | parent: iterations, final lb | children: iterations | parent + children: iterations, CPU |
|---|---|---|---|
| default | 430, 48750.8 | 60, 80, 30, 30 = 200 | 630, 236 s |
| `nodes=warm` (children 0.5 -> 0.3, parent default) | 430, 48750.8 | 90, 100, 40, 60 = 290 | 720, 267 s |
| all calls (`two,f=1.6667,t=60`) | 370, 48748.2 | 90, 110, 40, 60 = 300 | 670, 246 s |
| `nodes=cold` (parent 0.5 -> 0.3, children default) | 370, 48748.2 | 60, 80, 30, 30 = 200 | 570, 213 s |
| `two,f=0.7,t=100000,nodes=warm` (children 0.21, parent default) | 430, 48750.8 | 50, 60, 30, 30 = 170 | 600, 227 s |
| **`two,f=1.6667,t=60,w=0.7`** (parent 0.5 -> 0.3, children 0.21) | 370, 48748.2 | 50, 60, 30, 30 = 170 | **540, 203 s** |

Parent all32 node 6 (`20:12,40:28,17:25`, a miss), children 43:18, 43:7,
43:22, 43:19 (all32 nodes 25-28), `SUFFIX=_p6`:

| config | parent: iterations, final lb | children: iterations | parent + children: iterations, CPU |
|---|---|---|---|
| default | 370, 48726.4 | 30, 50, 130, 40 = 250 | 620, 225 s |
| `nodes=warm` (children 0.5 -> 0.3) | 370, 48726.4 | 40, 70, 160, 60 = 330 | 700, 250 s |
| `two,f=0.7,t=100000,nodes=warm` (children 0.21) | 370, 48726.4 | 30, 40, 100, 40 = 210 | 580, 211 s |
| **`two,f=1.6667,t=60,w=0.7`** | 330, 48727.4 | 30, 40, 100, 40 = 210 | **540, 201 s** |

Both families agree: the cold schedule slows every warm child, beta 0.21 on
warm children speeds up or ties every one, and the parent's earlier give-up
under the schedule does not cost its children anything (`nodes=cold` and
`w=0.7` leave the children's counts as with the default parent). This fits
the earlier tuning: 0.21 lost on cold wil50 nodes (`../../solvers/accel`,
-79 at iteration 200), while 0.21 with `-r 120` took 20% off the
warm-started nug30 search (`../../tuning/goux`, `RUNNING.md`).

### kra30a root (robustness) [m]

The memo keeps the controller only if it lifts the kra30a root bound above
QPB (about 68,600; the default stops near 17,368, `dnn_bound/cpp_report.md`).
Root, cold, production flags, no cutoff, 600 iterations (`kra30a/`,
`kra30a_curves.txt`):

| it | default | `two,f=1.6667,t=60` | fixed 0.5 | `rb,c=20` |
|---|---|---|---|---|
| 10 | 7,896 | 8,762 | 8,762 | 7,896 |
| 100 | 17,020 | 17,239 | 17,223 | 17,020 |
| 200 | 17,309 | 17,338 | 17,343 | 17,309 |
| 600 | 17,364.6 | 17,365.3 | 17,366.4 | 17,364.6 |

All settings converge to the same 17,365: kra30a's weak bound is not a beta
scaling problem within a factor 1.7, and the controller sees a ratio near
its target (rhat/shat 24-29 against c = 20, nug20's calibration) and never
acts. Closed.

## Files

- `run_all32.sh`, `summarize.py`, `all32_summary.txt`, `subset_fix05_summary.txt`, `runs/`: all32.
- `nug20_curves.sh` (also runs other instances: `INST`, `ITERS`), `curves.py`, `nug20_curves.txt`, `nug20/`, `kra30a/`.
- `run_nug30.sh`, `nug30_nodes.txt`, `nug30_nodes_big.txt`, `nug30/`; `run_warm.sh`, `warm/`; `summarize_n.py`.
- `ratios.py`: the residual ratios of the `bsched` lines.
- `beta_change_check.C`, `beta_change_check.out`: the hand check of a beta change.
- `validate.sh`, `validate_warm.sh`, `validate/`.

The scripts take the binary from `BIN` (default `latest/dnntest`). The
default runs (`runs/base_*`, `nug20/oldbin.log`) used the build of a955e95,
before the switch; the schedule runs a frozen copy of the build of the first
commit of the switch, and the runs with the `w` key (`warm/two05w07*.log`,
`warm/*_p6.log`, `nug20/v2*.log`) that of the second; `runs/basenew_17.log`
and `nug20/base.log` are the first build with the switch unset.

## Open

- A search A/B (wil50 subtree or nug30) of `two,f=1.6667,t=60,w=0.7` against
  the default: the warm factor is measured on 8 children of 2 parents only, and on nug30 the
  base is already 0.21 (w = 0.7 would give 0.147, untested).
- The switch point: t = 60 (switch at 83) was the memo's choice; t = 40 and
  t = 120 were tried on nug20 only (t = 120 slightly better there).
- The give-up rule fires earlier under the schedule; with `QAP_DNN_GU_ZONE` the
  near-cutoff deferral still applies, untested here.
