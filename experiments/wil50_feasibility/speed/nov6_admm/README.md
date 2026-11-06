# ADMM/PRSM convergence on warm-started solves: what is left to take (Nov 6, 2026)

Question: beyond the Nov 3 memos (`../../roadmap/admm/`) and the Nov 5 speed
work (`../beta_sched/`, `../rprsm/`, `../r32_noncert/`, `../cs2/`), what in the
research literature, or in our own logs, would cut the iterations a node's DNN
solve needs to fathom, above all for the warm-started solves that make up
nearly every node of the production search? arXiv and optimization-online
were reachable from this box, so the papers were read, not their snippets.
**[m]** = measured here (log named) or in a named report; **[e]** = estimate.
Host: 4-core 2.1 GHz Xeon with AMX, one process at a time on core 1, one BLAS
thread, two other agents' jobs on the other cores; iteration counts are
deterministic, milliseconds carry that noise. Production flags throughout
(bundle D: degree 3, `-L 2`, PRSM (1, 1), `-r 40`, `-H 300`, one CholQR pass,
`DNN_AMX=1`) and the `fast` arm's beta setting
`DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7` unless a row says otherwise.

## Verdict

1. **The literature has no new mechanism for a warm-started, early-stopped
   PRSM.** Everything published since the Nov 3 memos on splitting methods
   (the Halpern Peaceman-Rachford line, Nov 2026 reviews included; ergodic
   PRSM; optimal-rate PRS for smooth strongly convex terms; spectral
   penalties) is either closed here by a measurement or needs structure our
   problem does not have (section 2). What it does offer is a frame for the
   crossing certificate below: in branch and bound a node should stop at the first iterate
   whose dual bound crosses the incumbent (Chen, Ning, Goulart 2023), and our
   certificate is exactly such a bound but is evaluated only every 10
   iterations.
2. **Where the production search spends its iterations [m]:** in the `fast`
   arm of `../../pipeline_ab/` (node-10 subtree, every node warm) 61% of all
   DNN iterations are first solves that fathom, with a median of **30**
   iterations (178 of 226 fathom by iteration 50); 10% are first-solve misses
   (6 nodes, 5 of them 18-51 short, inside the cut zone) and 29% cut phases.
   Fathoms are found at a check and so run **3.7-4.1 iterations past the
   point where the certificate crossed the cutoff** (interpolated; 97
   fathoms in four log sets, section 3). That is 7.5% of the fathoming
   iterations and about 4.6% of all iterations of the arm.
3. **The largest finding: the warm-call beta (`w` of `DNN_BETA_SCHED`), the
   knob the Nov 5 gain rests on, is far from its optimum.** The production
   value 0.7 was the only one ever tried; swept over 0.15-0.85 on the 24 warm
   children [m], the response is monotone all the way down: the same 15
   fathoms in 1,550 (0.85), 1,360 (0.7), 1,080 (0.5), 880 (0.35), **780
   (0.25, -43%)** and 700 (0.15, -49%) iterations, every fathom earlier at
   every step down to 0.25, and every one of the 9 misses higher at
   iterations 50, 100 and 200 (at 0.25: median +46, +28, +19). Between 0.25
   and 0.15 the fathoms still gain 10% while two misses lose 1-2 units at
   200, so **0.25 is the value to carry (0.15 the second point)**. On the
   `fast` arm's mix that is about **1.35x in DNN iterations [e]** before any
   effect of the misses' higher bounds on the phases and children. The `w`
   key also governs cut rounds (`fast` ran them at 0.7); the new `wc=` key
   keeps them there while the children run lower, and the arms `wlow`,
   `wlowc` and `best` in `../../pipeline_ab/run_pair.sh` are the laptop test
   (section 7).
4. **Built and measured: `DNN_ADMM_CROSS=2`**, a certificate at the predicted
   crossing iteration between two regular checks. On the 24 warm children
   [m]: the same 15 fathoms in **1,297 iterations against 1,360 (-4.6%, 4.2
   iterations per fathom)**, 11 of 15 earlier by 3-9 iterations, none later,
   for 15 extra certificates (one per fathom; none on the 9 misses, whose
   traces are identical); the fathoms' CPU 456 -> 433 s (**1.052x**). On four
   cold fathoming nodes: 1,130 -> 1,117 iterations, CPU equal within noise.
   Scaled to the `fast` arm's mix: about **1.05x in iterations and 1.03x in
   CPU** [e]. Nothing else can change by construction. The early check at
   iteration 5 for warm calls (`DNN_ADMM_CROSS_FIRST=5`, built so that the
   first regular check has two points) saves the same iterations (1,300)
   for 40 extra certificates instead of 15 and halves the CPU gain (1.031x):
   not worth it.
5. **Two warm-transient probes:** damping both dual steps for a warm child's
   first 10 iterations (`DNN_ADMM_WARM_DAMP=0.5:10`) is **killed** [m]: 2 of
   8 fast children later by 10, none earlier (section 5.3). The two-phase
   warm beta (`DNN_ADMM_WARM_BETA`) was built and not run: with the response
   to beta monotone over 0.25-1.0 it has nothing to add until the optimum is
   bracketed. Both stay as switches, off.
6. **Closed by arithmetic, no run needed:** an earlier give-up for warm misses
   (10% of the arm's iterations, and 5 of the 6 misses are inside the cut
   zone, where a give-up is deferred and the phase needs the converged state;
   section 4).

## 1. The production node mix (the `fast` arm, cloud box, Nov 5) [m]

`tracemix.py` on `../../pipeline_ab/runs_cloud/fast/trace.txt` (232 nodes,
18,820 DNN iterations, all 226 fathoming first solves warm):

| class | nodes | iterations | share |
|---|---|---|---|
| first solves that fathom | 226 | 11,470 | 61% |
| first solves that miss (give-up at 300-330) | 6 (5 warm) | 1,920 | 10% |
| cut phases | 6 | 5,430 | 29% |

Fathom iteration of the 226 first solves: median 30, mean 50.8, quartiles 20
and 50; by check: 10: 7 nodes, **20: 75, 30: 62, 40: 24**, 50: 10, 60: 15,
70: 5, 80: 2, above 100: 23 nodes (4,870 iterations). The `full` arms (cloud
and laptop) have the same shape (median 30, mean 57-58). The six misses end
18, 28, 31, 31, 51 and 214 short.

So the production node is a warm child that fathoms in 20-40 iterations, and
whatever is saved per such solve has to be a fixed number of iterations, not
a rate: the solve is over before any rate matters.

## 2. Literature, current to Nov 2026, against our code and logs

| paper | what it shows | here |
|---|---|---|
| Chen, Sun, Yuan, Zhang, Zhao, *Recent developments in HPR methods for large-scale LP* (arXiv 2611.01872, Nov 1 2026) and *... convex composite QP* (2611.01887); HPR-LP (2408.12179, v2 2025); *PRSM converges ergodically* (2501.07807, 2025) | Halpern Peaceman-Rachford with semi-proximal terms: O(1/k) on the KKT residual; the practical gain comes from **restarts** (sufficient decay, necessary decay and long inactivity of a KKT-residual merit function, checked every 150 iterations) with the **penalty updated only at restarts**; the ergodic PR sequence converges and, restarted, beats pointwise and ergodic DR on LP benchmarks; 2.4-5.7x over PDLP on GPUs at 1e-8 | Closed with the same evidence as Nov 3: our Halpern anchoring with restarts on the fixed-point residual was 1.9x slower to fathom (`../../solvers/accel/`), because PRSM contracts geometrically on fathoms (gap ratio 0.34-0.55 per 50 iterations) and an anchored O(1/k) scheme can only slow a geometric transient; the ergodic (averaged) multiplier as the certificate was below the pointwise one at the fathom check (`../../roadmap/certificate/REPORT.md` #4); the penalty-at-restart rule is the residual-balancing controller whose normalized ratio does not separate transient from tail on wil50 children (`../beta_sched/README.md`). Nothing to re-run. |
| Briceño-Arias, Roldán, *Optimal leveraging of smoothness and strong convexity for PRS* (2601.12190, Jan 2026) | an optimal linear rate for PRS by moving quadratic terms between the two functions | needs two smooth strongly convex terms; ours are an indicator of the PSD cone and an indicator of a polytope with a linear objective. No handle. |
| Chen, Ning, Goulart, *A unified early termination technique for primal-dual algorithms in MICP* (2303.09337) | in branch and bound, stop a node's primal-dual solve at the first iterate whose dual bound exceeds the incumbent; an O(n) check per iteration; "a moderate reduction in total iterations" | **This is candidate 1.** Our certificate is a valid bound at any iterate (`docs/dnn_bound/dnn_bound.tex` Prop. 2.2) but costs 0.4-1.1 iterations (125 ms far from the cutoff, `../r32_noncert/runs/base2_19.log`: 1.62 s over 13; plus the exact LAP pass within 48.8 of the cutoff, about 0.25 s, `../../roadmap/certificate/REPORT.md`), so it runs every 10 iterations and the fathom lands 0-9 iterations late. Their every-iteration check is too dear here; predicting the crossing from the trajectory and certifying once there is the version that pays. |
| Sambharya, Hall, Amos, Stellato, *Learning to warm-start fixed-point optimization algorithms* (2309.07835, JMLR 2024) | for a linearly convergent operator the iterations to a tolerance grow with the log of the warm start's distance to the fixed point | the frame for "what to carry": a better start saves a fixed number of iterations, about log(2)/log(1/rho) per halving of the distance (30-60 iterations at our per-iteration rho of 0.979-0.988). That is what the parent-Y start measured (`../../warm_y/`: -24% on 300-500-iteration fathoms, 10-20 iterations on 60-90-iteration ones). Learned starts need thousands of solved parents per node class; not for this run. |
| Xu, Figueiredo, Goldstein (AADMM, 1605.07246); Xu (2503.06185, 2025) | spectral (Barzilai-Borwein) penalty selection | closed by the Nov 3 chain: "only if residual balancing survives" (`../../roadmap/admm/splitting.md` #5), and residual balancing was killed on wil50 (`../beta_sched/`). |
| Sinjorgo, Sotirov, Vera (2506.08648, 2025) | ADMM for Lasserre-level SDP bounds, warm-started from the level-1 (theta) solution | a warm start across relaxations; ours already carries the parent's Z and Y and the inherited pool's multipliers. Nothing new. |
| de Meijer, Siebenhofer, Sotirov, Wiegele (2410.04997, v4 Nov 2025) | a PRSM for a DNN relaxation with Chvatal-Gomory cuts (quadratic minimum spanning tree) | a bounding paper, same splitting family, no B&B warm start, no iteration counts for our class. |
| QAP-DNN since the Nov 3 memo | 2403.11210 (rank-one DNN reformulations, penalties), 2401.09439 (tai256c as a BQOP), nothing else 2024-2026 on arXiv or optimization-online for PRSM, DNN or QAP (listings `dl/s6, s11, s13, s24, oo1-4` in the session) | nothing to take |
| Qu, Sun, Xu, *Progressive bound strengthening via DNN cutting planes* (2510.02948, Oct 2025) | cuts for nonconvex QP from DNN relaxations | the cuts agent's topic |

Checked and found unchanged in `latest/dnn.C`: the PRSM loop certifies every
`check_every` = 10 iterations (`bnb.C` never changes it); the certificate uses
the R-step multiplier of its own iteration only; the give-up rules read the
running best at the checks (`besthist`, 50-iteration windows); gamma1 = gamma
= 1 is the largest dual step of the He-Ma-Yuan domain (Nov 3); the child warm
start maps Z entrywise with the fixed row folded into row 0 and conditions Y
(`dnn_warm_child`, ymode 2).

## 3. The iterations past the crossing [m]

`crossing.py` (here) reads a `-v` log, takes each fathomed solve's certified
value at the fathom check F and at F - 10, interpolates the iteration at
which it crossed the cutoff, and reports F minus that:

| log set | fathoms | mean past the crossing | median | mean F |
|---|---|---|---|---|
| warm children, `../beta_sched/warm/*.log` (all configurations) | 40 | **3.72** | 3.20 | 58 |
| warm children, `../../warm_y/runs/child_*.log` (three starts) | 21 | 3.89 | 3.79 | 215 |
| cold all32, `../beta_sched/runs/base_*.log` | 18 | 3.84 | 2.90 | 235 |
| cold all32, `../beta_sched/runs/two05_*.log` | 18 | 4.06 | 4.02 | 205 |

The same 3.7-4.1 on every class (a uniform offset in a 10-iteration interval
would give 4.5; the certificate gains decelerate, so the crossing sits late in
the interval slightly more often). On the `fast` arm's 226 first-solve
fathoms that is about 860 iterations, 7.5% of their 11,470 and 4.6% of the
arm's 18,820 [e from m].

## 4. Candidates, the arithmetic, and what kills each

| # | idea | class it acts on | expected | what kills it | status |
|---|---|---|---|---|---|
| 1 | **`DNN_ADMM_CROSS=n`:** at every check with the running best below the cutoff and rising, predict the crossing iteration from the last two certified values (the gain per iteration, decayed by the ratio of the last two gains) and certify there if it falls before the next regular check; at most n extra certificates per interval; an extra check tests the fathom only (give-up and stall rules see the regular checks alone, so nothing else can change) | every fathom (61% of the arm's iterations) | 3.8 iterations per fathom less the extra certificates, about 0.7 iteration each; **1.04-1.05x in iterations, about 1.03x in CPU on the arm [e]** | a fathom later than without it (impossible by construction), misses changed (impossible: the iterates are untouched), or extra certificates costing more than the iterations saved | **built; measured below** |
| 2 | the warm-call beta: `w` = 0.15 to 0.85 against 0.7 | every warm solve | the Nov 5 factor (1.14x in the search) came from the only value tried; +/-5% [e] before the sweep | w = 0.7 best or within noise on both sides | **measured below: monotone, 0.25 gives -43% iterations (section 5.1)** |
| 3 | **`DNN_ADMM_WARM_DAMP=D:K`:** both dual steps scaled by D for the first K iterations of a warm call. Mechanism: the warm child's certificate at iteration 10 is about 200 below the parent's final value (`../beta_sched/warm/base.log`: 48,551 against 48,751) because the start's primal residual (0.106 at iteration 10 against the parent's 0.025) drives the dual steps away from the parent's near-optimal multiplier; the cut phase's first-step overshoot had the same shape and its damping (`QAP_DNN_CUT_DAMP1=0.3`) turned -102 into +16 | warm solves | 0-20% on 30-60-iteration children [e], speculative | children fathom later or lb@50 lower | **built; killed (section 5.3)** |
| 4 | **`DNN_ADMM_WARM_BETA=F1:T:F2[:force]`:** a two-phase beta for warm calls (F1 x base until the first refresh at or after T, then F2 x base), the other way of damping the transient (beta scales the dual step too) | warm solves | as 3 | as 3 | built, not run: superseded by the sweep (5.1) |
| 5 | an earlier give-up for warm misses (`../../tail_probe/README.md` left it open: "25% of the arm in 300-iteration first solves") | warm misses | **nothing**: in `fast` the misses are 10% of the iterations, and 5 of 6 end 18-51 short, inside the cut zone (`QAP_DNN_GU_ZONE=0.003` defers a give-up there and the phase that follows starts from the first solve's state); the sixth (214 short) is one node | - | **closed by the trace** |
| 6 | HPR restarts, ergodic averaging, Halpern/momentum/Anderson | - | - | - | closed Nov 3, confirmed by the Nov 2026 papers (section 2) |
| 7 | spectral/BB penalty, residual balancing at restarts | - | - | - | closed (`../beta_sched/`) |
| 8 | inexact early R-steps (error schedules) | - | - | - | closed by data (`../../roadmap/admm/projection.md` idea 5) |

## 5. Measurements

All on the warm set (`run_warmset.sh`): six wil50 depth-3 parents, each a
clear miss given up at 300 iterations as in the search, and 24 of their
depth-4 children started from the parent's final Z and conditioned Y
(`dnntest -K ... -Q b`, the search's `QAP_DNN_WARM_Y=2`). Two families are
`../beta_sched/run_warm.sh`'s (all32 nodes 9 and 6; their 8 children all
fathom, in 30-100 iterations) and four are `../../warm_y/run.sh`'s (parents
10, 13, 2, 7 of `../../stronger_bounds/depth4/nodes.txt`; 7 of their 16
children fathom in 40-300 iterations, 9 miss). Production flags with AMX
(above); the parents run the production cold schedule
(`two,f=1.6667,t=60`) in every row, so only the children differ. Tables:
`summarize_warm.py prod cross crossf w085 w05 w035 w025 w015 damp --table`.
Iterations are deterministic; CPU seconds carry the shared box's noise
(about +/-3% between runs of the same configuration).

### 5.1 The warm-call beta [m]

| `w` (beta of warm calls, x 0.3 m/3) | fathoms (of 24) | iterations to fathom | earlier / later by >= 10 | misses (9): lb - prod at it 50 / 100 / 200, median [range] | children's CPU |
|---|---|---|---|---|---|
| 0.85 | 15 (same) | 1,550 (+14%) | 0 / 9 | -16 / -7.1 / -5.3 [all 9 lower] | 1,586 s |
| **0.7 (production, `fast`)** | 15 | **1,360** | - | - | 1,486 s |
| 0.5 | 15 (same) | 1,080 (-21%) | 12 / 0 | +22 / +9.9 / +8.5 [all 9 higher] | 1,339 s |
| 0.35 | 15 (same) | 880 (-35%) | 15 / 0 | +36 / +19 / +14.5 [all 9 higher] | 1,236 s |
| 0.25 | 15 (same) | **780 (-43%)** | **15 / 0** (median -20) | **+46 / +28 / +19** [+8.5, +25.5 at 200; all 9 higher] | 1,176 s |
| 0.15 | 15 (same) | **700 (-49%)** | 15 / 0 (median -20) | +53 / +38 / +22 [+6.0, +32; all 9 higher, but two misses 1-2 units below their 0.25 value at 200] | 1,199 s |

- The response is monotone from 1.0 down to 0.25 (1.0 -> 0.7 was -16% on 8 of
  these children in `../beta_sched/`), on every child: no child is later at
  a smaller beta, and the slow fathoms gain most (node 10's child 240 -> 180
  -> 120, node 13's 220 -> 180 -> 120, node 7's 300 -> 220 -> 150). The
  misses gain at every checkpoint, so within 300 iterations the tail does not
  pay for the smaller beta; one miss (node 2's child 8:21) gives up at 330
  instead of 740 iterations, 24 units higher at iteration 200.
- Why: a warm child starts from a near-optimal multiplier Z and a primal Y
  that is inconsistent with the child's new gangster zeros (primal residual
  0.1 at iteration 10 against the parent's 0.025, `prodk1_f9.log`). beta
  scales the dual steps (gamma beta (Y - VRV')): at 0.7 the first 10
  iterations move Z so far that the certificate drops 830 units below the
  parent's final value at iteration 1 and needs 20 iterations to recover
  (section 5.4); a smaller beta keeps Z near the parent's and lets Y
  re-settle with larger primal steps. Damping only the dual steps (section
  5.3) does not do it: the larger primal step matters as well.
- Cold calls are untouched (the cold schedule is the Nov 5 one). The `w`
  key also governs the rounds of a cut phase, which pass a warm start with
  Z; `fast` ran them at 0.7. The warm set has no phases, so the new key
  `wc=` (this commit) lets a search keep the rounds at the measured 0.7 while
  the children run lower: arms `wlow` (0.25 for both) and `wlowc` (children
  0.25, rounds 0.7) in `../../pipeline_ab/run_pair.sh`.
- Scaled to the `fast` arm's mix (section 1: 61% of the iterations in warm
  first-solve fathoms, 10% in warm misses whose 300 iterations do not
  change, 29% in phases, unchanged under `wc=0.7`): about **1.35x in DNN
  iterations on the node-10 subtree [e]**, before any effect of the misses'
  higher bounds on the phases and children. This is the one lever here that
  reaches the roadmap's 1.1x bar by itself, and the laptop A/B is its test.

### 5.2 The crossing-predicted certificate, `DNN_ADMM_CROSS=2` [m]

| configuration | fathoms | iterations to fathom | extra certificates | fathoms' CPU | misses |
|---|---|---|---|---|---|
| `prod` (w=0.7) | 15 | 1,360 | 0 | 456.1 s | - |
| `cross` (+ `DNN_ADMM_CROSS=2`) | 15 (same) | **1,297 (-4.6%, 4.2 per fathom)**; 11 of 15 earlier by 3-9, none later | 15 (one per fathom; 0 on the misses) | **433.4 s (1.052x)** | identical traces |
| `crossf` (+ `DNN_ADMM_CROSS_FIRST=5`) | 15 (same) | 1,300 | 40 | 442.4 s (1.031x) | identical traces |
| cold all32 nodes 11, 17, 22, 31 (`run_cold.sh`), base -> cross | 4 -> 4 | 1,130 -> 1,117 (nodes 11 and 31 by 6-7, 17 and 22 unchanged) | 4 | 356 -> 360 s (noise) | lb at 100 and 200 identical |

- The predictor behaves as designed: it fires once per fathom on average,
  never on a miss (the predicted need is always >= 10 iterations there),
  and the iterate is never touched, so misses and the give-up decisions are
  identical to the iteration. A missed prediction costs one certificate
  (0.13-0.38 s, 0.4-1.1 iterations); a hit saves 3-9 iterations.
- The first version predicted from the raw rate of the last interval and
  was 20-30 iterations early at the first regular check (`diag_k1.py` on
  `prodk1_f9.log`: child 1 predicted 28, true 50; child 2 29 against 57),
  because the certificate's gains fall to about 0.36 of the previous
  interval's; the shipped version decays the rate by the ratio of the last
  two gains, or by 0.5 when only two points exist.
- `DNN_ADMM_CROSS_FIRST=5` (an early check on warm calls so that the first
  regular check has two points) saves no more iterations and spends 25 more
  certificates: not worth it. Kept as a switch, off.
- On the `fast` arm's mix: about **1.05x in iterations, 1.03x in CPU [e]**.
  It stacks with 5.1 (independent mechanisms; with w=0.25 the fathoms are
  shorter, so the 4.2 iterations per fathom are a larger share).

### 5.3 Dual-step damping, `DNN_ADMM_WARM_DAMP=0.5:10` [m]: killed

Families f9 and f6 (8 children that fathom in 30-100 iterations): 190
against 170 and 220 against 210 iterations (2 children later by 10, none
earlier). The remaining families were not run. Smaller dual steps alone slow
the warm child; the gain of 5.1 needs the smaller beta's larger primal step
too. `DNN_ADMM_WARM_BETA` (a two-phase beta for warm calls) was built for the
same transient and not run: with the response to beta monotone over
0.25-1.0, a two-phase schedule has nothing to add until the optimum is
bracketed (the 0.15 row).

### 5.4 The per-iteration certificate of a warm child (`prodk1_f9.log`, `-k 1`) [m]

Parent all32 node 9 (final certificate 48,749.1 after 372 iterations, w=0.7
children):

| child | certificate at iteration 1 / 5 / 10 / 20, relative to the parent's final | crosses the cutoff at |
|---|---|---|
| 37:14 | -829 / -326 / -153 / -28 | 50 |
| 37:37 | -837 / -297 / -131 / -25 | 57 |
| 37:21 | -826 / -279 / -91 / +66 | 20 |
| 37:26 | -801 / -276 / -85 / +42 | 24 |

The first 10 iterations of a warm child are a transient: the certificate
(the R-step multiplier Z + beta (Y - VRV')) is 800 units below the parent's
value at iteration 1 because the residual term is large, and Z itself is
moved by it. The near-cutoff gains are 1-3 units per iteration with runs of
0-gain iterations after a refresh (iterations 42-46), which is why the
crossing predictor needs the decay.

### 5.5 Validity and the default path [m]

`validate.sh` (logs in `validate/`), on the final binary, every switch unset:
`make test` passes (brute force, resume, parallel, coordinator; the
regression check matches `bench/reference/` on 7 instances at tolerance
0.0%); `make dnncheck` OK; `make shipcheck` exit 0; `make continuationcheck`
PASS (72 identity cases, 648 safe certificates, max state error 0);
`python3 bench/bench.py check` exact. With the switches on:
`tests/dnn_test.py --parts ACDE` OK under `DNN_ADMM_CROSS=2` and under
`DNN_ADMM_WARM_DAMP=0.5:10 DNN_ADMM_WARM_BETA=0.5:10:0.7:1` (part A:
max(lb - optimum) = 0 over 2,592 node bounds each); `tests/brute_force.py
--solver-args "-x 0:20:1:200:1" --jobs 1` with `QAP_DNN_WARM_Y=2` (warm
children) 4,536 of 4,536 under the three `DNN_ADMM_*` switches together and
4,536 of 4,536 under `DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7`. The
default path is unchanged: `identity_check.sh` against the binary built
before the change gives identical traces (all fields but the CPU columns) on
the nug20 root (400 iterations), wil50 all32 node 19 (fathomed at 130) and
a warm child through BetaSched (`ident/`, three OK; run again on the final
binary, three OK).

Why validity cannot be affected: the crossing switch only adds evaluations
of the same certificate (Rump's test on the R-step multiplier of its own
iteration) and never touches an iterate, a history or a give-up decision;
the beta keys change only where the R-step runs a full eigendecomposition,
and Z is the unscaled multiplier (`../beta_sched/beta_change_check.C`); the
dual-step damping changes iterates only.

## 6. Switches (all off by default; `latest/dnn.C`, read at every `dnn_bound` call)

| switch | what |
|---|---|
| `DNN_ADMM_CROSS=n` | between two regular certificates, certify at the predicted cutoff crossing (at most n extra certificates per 10-iteration interval); `-v` marks the extra checks' lines with `extra` and prints their count at the end. Measured: n = 2 |
| `DNN_ADMM_CROSS_FIRST=k` | with it, on warm calls one more certificate at iteration k (measured: not worth it) |
| `DNN_BETA_SCHED=...,wc=F` | the warm-call factor for calls that carry a cut term (the rounds of a cut phase), so children (`w`) and cut rounds (`wc`) can differ; unset: `w` for both |
| `DNN_ADMM_WARM_DAMP=D:K` | warm-started calls: gamma1 and gamma scaled by D for the first K iterations |
| `DNN_ADMM_WARM_BETA=F1:T:F2[:force]` | warm-started calls: beta F1 x base until the first scheduled refresh at or after iteration T (force=1: at T itself, with a full refresh there), then F2 x base; overrides `DNN_BETA_SCHED`'s `w` for warm calls, cold calls untouched |

The certificate is untouched by all three: it is the R-step multiplier of its
own iteration, checked by Rump's test, whatever the dual steps or beta did
before (`docs/dnn_bound/dnn_bound.tex`, Section 2). A beta change happens
only where the R-step runs a full eigendecomposition (BetaSched's rule), and
Z is the unscaled multiplier, so nothing is rescaled
(`../beta_sched/beta_change_check.C`).

Default path unchanged [m]: `identity_check.sh` runs the binary built before the
change and the one after it, every switch unset, on the nug20 root (400
iterations, no cutoff), wil50 all32 node 19 (fathomed at 130) and a warm child
(all32 node 6's child 43:18 from its parent's Z and conditioned Y, the warm
path through BetaSched): the traces are identical in every field but the CPU
columns (`ident/`: nug20_root OK, wil50_n19 OK, wil50_warm OK).

## 7. Next test (laptop, `pipeline_ab/run_pair.sh`)

1. `run_pair.sh fast wlowc`: the warm-call beta 0.25 for children, 0.7 for
   cut rounds, against the production candidate, node-10 subtree, same
   machine and hour. Kill line 1.15x in DNN iterations on the matched nodes
   (the warm set says about 1.35x; the laptop's cut phases and deeper
   children are what the warm set does not cover). Then `wlow wlowc` (does
   the low beta help or hurt the cut rounds?) and `PAB_WARM_W=0.15
   run_pair.sh wlowc` as the second point.
2. `run_pair.sh wlowc best`: the crossing certificate on top. Expected about
   1.05x in iterations and 1.03x in CPU; kill line 1.02x in CPU on the
   matched nodes (it cannot lose iterations).
3. The tail class (`../../tail_probe/run_matched.sh` with the same
   environment) if 1 passes: the levers there act at level 5, where the nodes
   are exactly this warm class.

Here (cloud box), if the laptop is busy: the 47 children of node 10 through
`dnntest -K` (`run_warmset.sh` with `FAMILIES=p10` and the full child list)
give the depth-4 answer in about 40 minutes per configuration.

## 8. Files

- `run_warmset.sh`, `summarize_warm.py`, `warm/`: the warm set (six parents,
  24 children) and its tables.
- `run_cold.sh`, `summarize_cold.py`, `cold/`: four cold fathoming nodes of
  all32.
- `identity_check.sh`, `ident/`; `validate.sh`, `validate/`: the default-path
  check and the validity suite (section 5.5).
- `crossing.py`, `tracemix.py`, `diag_k1.py`: the log and trace analyses of
  sections 1, 3 and 5.4.
- `../../pipeline_ab/run_pair.sh`: the new arms `cross`, `wlow`, `wlowc`,
  `best` for the laptop (section 7).
- `latest/dnn.C`: `DNN_ADMM_CROSS`, `DNN_ADMM_CROSS_FIRST`,
  `DNN_ADMM_WARM_DAMP`, `DNN_ADMM_WARM_BETA` and the `wc` key of
  `DNN_BETA_SCHED`, all off by default.
