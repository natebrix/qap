# bound_lab / admm: is the PRSM at its floor on tai50b-type instances? (Nov 15 2026)

Agent "admm". One solver process at a time on the shared 4-core Xeon, so iterations are the
metric; CPU figures are indicative. [m] = measured (log named), [e] = estimated.
Baseline everywhere = tai50b's production keys (`experiments/estimates/tai50b/tai50b.par`:
PRSM (1, 1), `dnn.beta = 0.21`, `dnn.cheb = 3`, `dnn.tight_lam = 2`, `dnn.refresh = 120`,
`dnn.beta_sched = w=0.36`), as `base.par`, cold from zero unless a row says warm.

## Answer

- **Not a robust 1.5x. A lever of 1.2-1.6x in iterations exists on tai-b nodes that are far
  from their plateau, at the same cost per iteration: a diagonal congruence scaling of the
  lifted variable (`dnn.dscale = 0.5`, patch `admm.patch`).** On tai50b: 1.20-1.26x cold on
  depth-3 node 7, **1.58x on its warm depth-4 child (fathom at 100 iterations instead of
  158)** [m]. On the solved siblings: 1.4-1.5x on the tai25b root, 1.2-1.3x on the tai30b
  root and a depth-3 node cold, 1.4-1.5x on a warm depth-3 child [m]. Where the relaxation
  is nearly tight (deep on tai30b's optimal path) it loses a little, and a tai30b subtree
  took 36% more iterations with it; scaled for the first 100 iterations of each call and
  plain after (`dnn.dscale_until = 100`) the same subtree took 21% fewer (section 4). It is
  not the 5-10x that would move tai50b's estimate (floor 1e11 CPU-s) into reach, and its
  effect on tai50b's tree is not measured.
- **The convergence tail on tai-b is sublinear** (gap to the plateau x iteration is constant
  from iteration 400 to 2,000 on the tai25b root), unlike wil50's geometric tail. It is a
  property of the splitting on this data, not of the eigen-tracker, the certificate or beta:
  exact eigendecompositions give the same curve, and nothing in the (gamma1, gamma),
  momentum, Anderson or beta families beats PRSM (1, 1) by more than the scaling does.
- **beta**: cold, the best factor is instance-dependent (0.05 on the tai20b root, whose
  relaxation is tight; 0.2-0.5 on tai25b; 0.5 on the tai30b root, 1.6-2.0x over 0.21 there;
  flat 0.21-0.5 on a tai30b depth-3 node). Warm children: `w` is flat from 0.36 to 2.4 after
  iteration 100 and 0.36 is best before it. No change proposed for warm calls; for a cold
  tai-b root `dnn.beta = 0.35-0.5`.
- **tho40's regime (B)**: the 36 warm children of a give-up parent first cross the cutoff at
  iterations 5-16 (414 in total) and production fathoms them in 537, waiting for the checks
  at 10 and 20. `dnn.cross_first = 5` (not in tho40.par; default 0): 460. With a low beta for
  the first 30 iterations of a warm call as well (`DNN_ADMM_WARM_BETA=0.1:30:0.36:1`):
  **377, 1.42x fewer iterations, 1.2x less CPU on the children**, and the give-up node's
  bound at 200 iterations is unchanged (-14). About 12% of tho40's CPU [e]. Section 6.

## 1. Diagnosis (tai25b root, m = 25, K = 577; optimum 344,355,646)

`logs/tai25b_root_b0.21.log` to 1,500, `logs/tai25b_root_b0.21_long.log` to 12,000 [m].
Plateau taken as 338.00M (337.981M at 12,000, still gaining 7k per 1,000).

| it | certified lb | gap to plateau | gap x it | primal objective | pres | dY | npos |
|---|---|---|---|---|---|---|---|
| 100 | 309,577,523 | 28.42M | 2.8e9 | 264.7M | 8.4e-2 | 4.5e-2 | 410 |
| 200 | 327,182,549 | 10.82M | 2.2e9 | 304.9M | 5.1e-2 | 1.8e-2 | 400 |
| 400 | 333,594,704 | 4.41M | 1.8e9 | 325.6M | 3.4e-2 | 6.5e-3 | 401 |
| 600 | 335,181,559 | 2.82M | 1.7e9 | 329.8M | 2.9e-2 | 3.6e-3 | 404 |
| 1,000 | 336,350,312 | 1.65M | 1.6e9 | 333.1M | 2.2e-2 | 2.2e-3 | 405 |
| 1,500 | 336,914,272 | 1.09M | 1.6e9 | 334.6M | 1.9e-2 | 1.4e-3 | 408 |
| 2,000 | 337,221,567 | 0.78M | 1.6e9 | 335.4M | 1.7e-2 | 1.0e-3 | 410 |
| 4,000 | 337,734,145 | 0.27M | 1.1e9 | 337.0M | 1.3e-2 | 7.6e-4 | 411 |
| 8,000 | 337,937,968 | 0.06M | 5.0e8 | 337.6M | 1.0e-2 | 3.6e-4 | 417 |
| 12,000 | 337,980,955 | 0.02M | 2.3e8 | 337.9M | 8.7e-3 | 4.8e-4 | 418 |

- The gap falls like 1/k from 400 to 2,000 iterations (wil50: a ratio of 0.34-0.55 per 50
  iterations). At the tree's give-up (600) the bound is 2.8M below its own plateau, 30% of
  the 9.2M it is short of the cutoff: on tai-b the convergence speed is part of the bound.
- The certificate does not lag the iterate: the primal objective is further from the plateau
  than the certified value at every row.
- What it is not (all at 1,500 iterations on the same root, certified lb, [m]
  `logs/tai25b_root_*.log`):

| variant | lb at 400 | at 600 | at 1,000 | at 1,500 |
|---|---|---|---|---|
| baseline, PRSM (1, 1), tracker | 333,594,704 | 335,181,559 | 336,350,312 | 336,914,137 |
| exact ssyevd every iteration (`-e 1`) | 333,654,399 | 335,211,034 | 336,364,122 | 336,923,573 |
| exact dsyevd (`-e 0`) | 333,654,404 | 335,210,885 | 336,364,262 | 336,923,576 |
| filter degree 6 (`-d 6`) | 333,647,303 | 335,204,049 | 336,361,628 | 336,921,739 |
| (gamma1, gamma) = (0.9, 0.9) | 333,072,311 | 334,847,589 | 336,169,952 | 336,789,833 |
| (0.5, 1) | 331,795,493 | 334,150,829 | 335,787,704 | 336,538,610 |
| ADMM (0, 1.618) | 331,634,792 | 333,930,152 | 335,585,599 | 336,377,407 |
| ADMM (0, 1) | 327,284,095 | 331,814,599 | 334,577,209 | 335,791,203 |
| ADMM (0, 1) + momentum 1.05 | 332,860,174 | 334,743,544 | 336,093,457 | 336,724,210 |
| PRSM + momentum 1.05 | diverged (pres 3.6e20) | | | |
| PRSM + Anderson 5 | 333,078,790 | 334,854,790 | 336,064,994 | 336,576,682 |
| ADMM (0, 1) + Anderson 10, exact eig | 332,113,617 | 334,872,940 | 336,119,456 | 336,435,195 |

  The eigen-tracker costs nothing in bound here (it only inflates the printed pres: 3.4e-3
  with exact eigendecompositions against 1.9e-2); every damped or accelerated variant ends
  on the same slow curve, below PRSM (1, 1). So candidates (4) and (5) of the brief are
  closed for tai-b too: the tail is sublinear, but momentum and Anderson do not shorten it.
- Certificate-side extrapolation, also closed: with two stored certificate multipliers of one
  run, Z(t) = Z_k2 + t (Z_k2 - Z_k1), certified after projecting out the positive part of
  V'ZV: every t > 0 is far below t = 0 (k1 = 600, k2 = 800: 335.97M at t = 0, 333.81M at
  0.25, 319.7M at 1) [m] `plogs/tai25b_root_extrap.log`. The multipliers do not crawl along
  a line.
- At nodes (a linear term C from the fixed assignments) dY stays at 0.1-0.25 for hundreds of
  iterations (tai50b node 7: 0.17 at iteration 600). 95-99% of dY^2 is in the arrow (row and
  column 0), an undamped Peaceman-Rachford mode: (0.9, 0.9) removes it (dY 0.094 to 0.0046 on
  the tai30b node below) and the bound is 0.3M lower with it, and scaling the homogenizing
  coordinate (`dnn.dscale0` = 2, 4, 8, 0.5, 0.25) is worse at every value [m]
  `logs/tai30b_opt3_ds0_d*.log`, `logs/tai30b_opt3_g*.log`. The mode is harmless; the
  printed dY is not a convergence measure at nodes. (This is `splitting.md`'s idea 4; its
  zero-code test would have said "go" here, and the measurement says no.)

## 2. beta, cold (certified lb; optimum in the heading)

[m] `logs/<inst>_root_b<F>.log`, `dnn.beta = F`, everything else baseline.

| root | F | at 100 | at 200 | at 400 | at 600 | at 1,000 | at 1,500 |
|---|---|---|---|---|---|---|---|
| tai20b (122,455,319) | 0.02 | 58.9M | 78.1M | 98.4M | 110.2M | 121.92M | 122,455,313 |
| | 0.05 | 84.1M | 104.6M | 121.52M | 122.440M | 122.4550M | fathomed at 1,430 |
| | 0.1 | 103.7M | 120.8M | 122.410M | 122.439M | 122.4526M | 122,455,053 |
| | 0.21 | 119.7M | 121.99M | 122.365M | 122.422M | 122.4387M | 122,450,199 |
| | 0.5 | 120.6M | 121.89M | 122.143M | 122.233M | 122.4079M | 122,427,234 |
| | 1.0 | 119.3M | 121.54M | 121.872M | 122.067M | 122.2204M | 122,360,582 |
| tai25b (344,355,646) | 0.05 | 228.7M | 278.7M | 314.7M | 325.2M | 332.21M | 334,735,520 |
| | 0.21 | 309.6M | 327.2M | 333.59M | 335.18M | 336.35M | 336,914,137 |
| | 0.5 | 321.0M | 331.1M | 334.40M | 335.41M | 336.33M | 336,851,520 |
| | 1.0 | 316.3M | 330.5M | 333.28M | 334.70M | 335.69M | 336,295,289 |
| | 2.0 | 283.7M | 321.3M | 331.59M | 332.85M | 334.48M | 335,340,357 |
| tai30b (637,117,113) | 0.05 | 275.2M | 376.0M | 471.0M | 515.4M | 558.5M | 582,006,681 |
| | 0.21 | 470.6M | 542.3M | 583.37M | 592.17M | 595.94M | 597,471,180 |
| | 0.5 | 542.8M | 582.3M | 593.25M | 596.08M | 597.57M | 598,223,990 |

- tai20b's root relaxation is tight (the bound reaches the optimum less 6), and there a
  small beta has a slow start and then converges fast. tai25b and tai30b are not tight, and
  there a small beta is simply slower. One beta does not fit the family.
- tai30b root, 0.5 against 0.21: reaches the 0.21 run's value at iteration k in k/2.0
  (k = 200), k/1.7 (600), k/1.6 (1,500). tai25b: no gain from 0.5. tai30b depth-3 node
  `0:3,1:7,2:10` (on the optimal path): 0.21, 0.35, 0.5, 1.0 give 633.74M, 633.87M,
  633.81M, 633.30M at 600 [m] `logs/tai30b_opt3_b*.log`: flat.
- A residual-balancing controller (candidate 3) was not built: the printed ratio
  pres/(beta dY) is dominated by the tracker's error (root) and by the arrow mode (nodes),
  so it does not point at the better beta, and at nodes beta is flat anyway.
- tai50b node 7 cold: the Nov 8 feasibility settings (beta 0.3, degree 4, refresh 20:
  `../../tai50b_feasibility/depth3/node7.log`) are ahead of the production keys by 4.2M at
  400 iterations and 1.7M at 600 (451.18M and 454.17M against 446.95M and 452.47M) [m].
  Not split into beta and tracker here; it concerns cold solves only (the root).

## 3. The lever: diagonal congruence scaling (`dnn.dscale`)

`splitting.md` closed equilibration because R -> VRV' is an isometry and a non-scalar
penalty breaks the eigen-solution of the R-step. A **congruence** does not break it:

- Yt = D Y D with D = diag(d0, d1 (x) d2) (one weight per facility, one per location).
  Then Yt = Vt Rt Vt' with Vt = orth(D V) and Rt = T R T' psd, so the R-step is still a
  plain PSD projection, in the basis Vt. In the original variables the penalty becomes
  entrywise, beta d_a^2 d_b^2: large where the objective's entries are large.
- The fast transform survives: range(D1 Vb) is the complement of D1^-1 e, so the
  Householder reflector H becomes one reflector per tensor mode, and the 45-degree rotation
  of the pair (0, ll) becomes one more reflection supported on 2m entries. Same O(m^4)
  products; 364-373 ms per iteration against 372-397 at m = 47 [m].
- The scaled problem: min <D^-1 L D^-1, Yt>, 0 <= Yt_ab <= d_a d_b, Yt_00 = d0^2, gangster
  zeros. Weights: d1_i = (r_i / geometric mean)^p with r_i the 2-norm of row i of
  (|A| + |A'|)/2, d2 from B likewise, d0 = 1. p = 0.5 is one pass of symmetric equilibration
  of L. On tai25b d1 lies in [0.93, 1.22] and d2 in [0.58, 1.40]: the flows carry it.
- **The certificate is untouched.** The multiplier is converted to original units as it is
  formed (Zc = D Zc_t D; V'ZcV = T'(Vt' Zc_t Vt)T <= 0 by congruence) and certified by the
  existing code with the plain transform. Start and save states are in original units, so
  warm starts, state shipping and the x-bar branching read what they read today.
- Prototype first (`proto.py`, dense V, exact eigendecompositions, Gilmore-Lawler
  certificate): its unscaled run reproduces dnntest's iterates to every printed digit, and
  the C++ patch reproduces the prototype's scaled iterates to every printed digit
  (`logs/tai25b_root_ds0.5_e0.log` against `plogs/tai25b_root_p0.5_F0.21.log`). With
  `dnn.dscale = 0` the patched binary prints the baseline's lines unchanged.

### 3.1 Results: iterations for the scaled run to reach the baseline's value at k

"ratio" = k / (iteration at which `dnn.dscale = 0.5` reaches the baseline's certified value
at k, interpolated between checks); "delta" = scaled minus baseline at k. Same node, same
start, same beta. `reach.py base.log scaled.log`.

| node | k | baseline at k | ratio | delta at k |
|---|---|---|---|---|
| tai25b root, cold (`tai25b_root_b0.21`, `_ds0.5`) | 200 | 327.18M | 1.52x | +4.58M |
| | 400 | 333.59M | 1.49x | +1.39M |
| | 600 | 335.18M | 1.40x | +0.80M |
| | 1,000 | 336.35M | 1.37x | +0.50M |
| | 1,500 | 336.91M | 1.42x | +0.34M |
| tai30b root, cold (`tai30b_root_b0.21`, `_ds0.5`) | 200 | 542.27M | 1.22x | +18.0M |
| | 600 | 592.17M | 1.26x | +2.28M |
| | 1,500 | 597.47M | 1.33x | +0.61M |
| tai30b depth 3 `0:3,1:7,2:10`, cold (`tai30b_opt3_base`, `_ds0.5_d1`) | 200 | 627.58M | 1.26x | +2.55M |
| | 400 | 632.43M | 1.32x | +0.93M |
| | 600 | 633.74M | 1.29x | +0.58M |
| the same node, **warm** child of `0:3,1:7` (600-iteration parent; `tai30b_P1_w0.36`, `_w0.36_ds0.5`) | 100 | 632.49M | 1.50x | +0.69M |
| | 200 | 633.57M | 1.46x | +0.56M |
| | 400 | 634.47M | 1.47x | +0.40M |
| | 600 | 634.93M | 1.42x | +0.34M |
| **tai50b depth 3 node 7** `0:44,49:49,25:31`, cold (`tai50b_n7_base`, `_ds0.5`) | 200 | 418.68M | 1.20x | +10.6M |
| | 400 | 446.95M | 1.22x | +3.51M |
| | 600 | 452.47M | 1.26x | +1.65M |

Fathoming nodes (iterations to fathom, baseline -> scaled) [m]:

| node | baseline | `dnn.dscale = 0.5` | ratio |
|---|---|---|---|
| **tai50b** depth 4 `0:44,49:49,25:31,1:23`, warm (`tai50b_n7kid_w0.36*.log`) | 158 (54.9 CPU-s) | 100 (34.5 CPU-s) | **1.58x** |
| tai30b depth 3 `0:3,1:7,2:14`, warm | 94 | 80 | 1.18x |
| tai30b depth 3 `0:3,1:7,5:6`, warm | 27 | 22 | 1.23x |
| tai30b depth 3, three children of `4:9,11:20`, warm | 10, 10, 20 | 10, 10, 20 | 1.0x (first or second check) |

- Validity: the tai30b node `0:3,1:7,2:10` is on the optimal path (optimum 637,117,113); its
  best certified value over all runs is 635,288,753, and in the subtree runs of section 4 no
  on-path node is certified above the optimum. No certificate failed in any run
  (`cert_fail=0`; the diverged momentum run aside).
- p and beta with the scaling (tai30b root at 1,500: p = 0.25, 0.5, 0.75, 1.0 give 597.90M,
  598.08M, 598.07M, 597.72M; the prototype on the tai25b root at 800: 336.39M, 336.58M,
  336.60M, 336.34M against 335.97M unscaled). **p = 0.5 to 0.75 is a flat optimum.** Row
  sums or row maxima in place of row 2-norms: no better (`plogs/`).
- The scaling and beta are not additive on the tai30b root: scaled at beta 0.35 against
  unscaled at 0.5 is 1.05-1.09x. On the depth-3 node, against the best unscaled beta (0.35),
  the scaling is still 1.11-1.28x. On warm children it is 1.4-1.5x at every w tried
  (`tai30b_P1_w1_ds0.5`: 1.5-1.8x against the w = 0.36 baseline).
- From a warm state the same factor holds without the transient: continuing the unscaled
  600-iteration state of the tai25b root for 600 more, the scaled iteration reaches the
  unscaled continuation's value at +600 after +390 (prototype, `plogs/tai25b_root_cont600.log`).
- The tail stays sublinear with the scaling; it is a constant factor, not a change of rate.

### 3.2 Warm-call beta with and without the scaling (tai30b, child `2:10` of `0:3,1:7`)

| config | at 10 | at 50 | at 100 | at 200 | at 400 | at 600 |
|---|---|---|---|---|---|---|
| w = 0.36 (production) | 615.5M | 630.50M | 632.49M | 633.57M | 634.47M | 634.93M |
| w = 1 | 591.3M | 630.24M | 632.78M | 633.69M | 634.57M | 634.98M |
| w = 2.4 | 543.8M | 628.53M | 632.57M | 633.63M | 634.38M | 634.87M |
| w = 0.36, dscale 0.5 | 617.8M | 631.91M | 633.18M | 634.13M | 634.88M | 635.27M |
| w = 1, dscale 0.5 | 597.1M | 631.66M | 633.42M | 634.26M | 634.93M | 635.29M |
| w = 1.67, dscale 0.5 | 577.3M | 630.71M | 633.28M | 634.17M | 634.89M | 635.23M |

w = 0.36 stays: best in the first 50 iterations (where the fast fathoms are), level after.

## 4. In the tree

Two subtrees with bnbdrv from the patched worktree, tai50b's settings (`bb_*.par`: tai50b.par
with `dnn.depths` to the bottom, give-up 600, warm starts), the optimum as the incumbent,
per-node traces (`logs/bb_*_F3_*.dnn.txt`). Both are depth-3 subtrees of the optimal path,
the only kind small enough to finish here; `bbsum.py`.

tai30b, `-F 0:3,1:7,2:10` [m]:

| arm | DNN nodes | fathomed | not fathomed | iterations on fathomed nodes | on the others | total | 
|---|---|---|---|---|---|---|
| baseline | 233 | 218 | 15 | 7,079 | 10,520 | 17,599 |
| `dnn.dscale = 0.5` | 267 | 245 | 22 | 8,376 | 15,490 | 23,866 (+36%) |
| `dnn.dscale = 0.5`, `dnn.dscale_until = 200` | 236 | 220 | 16 | 7,298 | 11,500 | 18,798 (+7%) |
| `dnn.dscale = 0.5`, `dnn.dscale_until = 100` | 195 | 184 | 11 | 5,835 | 8,130 | **13,965 (-21%)** |
| `dnn.dscale = 0.25` | 244 | 228 | 16 | 8,322 | 11,320 | 19,642 (+12%) |

The on-path node at each depth, certified value at its give-up minus the baseline's:

| depth | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|
| baseline's shortfall | 3.38M | 957k | 162k | 111k | 82k | 21k | 50k | 5.6k | 17.7k | 13.8k |
| dscale 0.5 | +580k | +185k | -14.5k | -2.0k | -7.1k | -5.2k | -1.1k | -8.9k | -15.7k | -23.0k |
| dscale 0.5 until 200 | +633k | +199k | +8.7k | +4.3k | +0.3k | +1.7k | +19.9k | -35.1k | -6.3k | -7.7k |
| dscale 0.5 until 100 | +556k | +121k | +8.3k | +5.1k | +1.4k | -0.5k | +16.1k | -6.5k | +7.8k | +11.4k |
| dscale 0.25 | +381k | +130k | +0.3k | -1.7k | -5.4k | -7.4k | -18.9k | -15.3k | +2.1k | +2.4k |

- **The scaling helps where the node is far from its plateau (depths 3-4 here, shortfall
  0.15-0.5% of the cutoff) and loses a little where the relaxation is nearly tight (depth 5
  and below, shortfall under 0.03%)**: there the scaled run is 1k-23k lower at the give-up,
  and four nodes the baseline fathomed after 250-330 iterations (or 28) gave up, so the
  spine branched into 34 more nodes. This is the tai20b-root regime of section 2 (tight
  relaxation, where the tail wants a different penalty).
- `dnn.dscale_until = T` (scaled for the first T iterations of a call, plain after; the
  state converts exactly) keeps the gain at depths 3-4 and removes most of the loss below.
  T = 200: 236 nodes against 233, 7% more iterations. **T = 100: 195 nodes, 21% fewer
  iterations**, the spine ending four levels sooner. One subtree, and a spine's length hangs
  on single near-tight nodes, so read it as "no longer a loss", not as 1.26x. T = 100 is the
  arm to carry.
- The fathoming children here take 10-30 iterations (51 at the first check, 93 at the
  second, of 214): the scaling does not change them (6,225 against 6,173 iterations on the
  214 nodes both arms fathom), as on tho40 (section 6).
- tai25b, `-F 0:3,1:14,2:9`: too small to say anything (23 nodes; the subtree's root runs
  3,110 iterations within 5,429 of the cutoff in the baseline and gives up at 1,580 with the
  scaling, 127k lower, and the on-path child then takes one more level): 4,028 iterations
  against 3,655 [m].
- **What this does not measure is tai50b's own regime**: its CPU is at levels 4-5 in nodes
  0.5-3% short that run to the 600-iteration give-up or fathom after 100-600 iterations, the
  regime of the tai30b depth-3/4 rows and of the two tai50b measurements of section 3.1.
  The subtrees above are spines of near-tight nodes and quick fathoms. A tai50b estimate with
  `dnn.dscale = 0.5`, `dnn.dscale_until = 100` or `200` (the dives of `experiments/estimates/tai50b`)
  is the test that would settle it; from the per-node numbers I expect 1.2-1.5x in CPU at
  levels 4-5 [e], which does not change "out of reach".

## 5. Keys per instance

| instance | change | basis |
|---|---|---|
| tai50b, and Taillard type b generally | `dnn.dscale = 0.5`, `dnn.dscale_until = 100` (needs `admm.patch`) | sections 3 and 4; confirm on tai50b dives first |
| tho40 | `dnn.cross_first = 5`; `DNN_ADMM_WARM_BETA=0.1:30:0.36:1` (environment) | section 6; one parent's 36 children, confirm on the tuning set |
| tai-b roots (cold solves) | `dnn.beta = 0.35` to `0.5` | tai30b root only; node-level flat |
| wil50, sko42, sko49, tho40, nug | `dnn.dscale = 0` (default) | tho40 measured: no effect (section 6). The others not measured; their row norms are nearly uniform, so D is close to the identity |
| kra30a | not measured | a candidate for the same scaling (it "does not converge at the default beta") |

The scaling is off with the cut hook (CP-ADMM), the accelerations and the prototype
certificate; the patch does not implement those paths.

## 6. tho40's regime (B)

Checked first: `experiments/estimates/tho40/tho40.par` has no `dnn.cross_first`; the default
is 0 (off); `dnn.cross = 2`, checks every 10 iterations, `w = 0.36`.

The node: tho40t `26:10,8:0,39:32,38:7` (depth 4, m = 36; in tho40ab-32 it gave up 406 short
and its 36 children all fathomed). Reproduced with dnntest: the parent cold for 600 iterations
(240,132.5, 382 short, as in COMMON.md), then the 36 children warm from that state (Z and
conditioned Y, `-Q b`), production keys (`tho40_base.par`). `q2.sh`, `q7.sh` to `q9.sh`,
`tho.py`, `logs/tho40_d4_*.log` [m]. The children's CPU is on the shared machine.

| arm | iterations to fathom, 36 children | histogram (iteration: children) | extra certificates | children's CPU s |
|---|---|---|---|---|
| a certificate at every iteration (`-k 1`): the first crossing | 414 | 5-9: 14, 10: 7, 11-16: 12, 25-33: 3 | | 49.5 |
| production | **537** | 10: 21, 20: 12, 26: 2, 35: 1 | 4 | 52.3 |
| `dnn.cross_first = 5` | 460 (-14%) | 5: 1, 10: 20, 11-16: 12, 26-35: 3 | 51 | 45.3 |
| `dnn.cross_first = 3` | 472 | 10: 21, 11-19: 12, 29-33: 3 | 57 | 48.1 |
| w = 0.15, `-k 1` | 332 | 5-9: 24, 10-13: 9, 18-25: 3 | | 45.5 |
| w = 0.15, `dnn.cross_first = 5` | 387 (-28%) | 5: 5, 10: 20, 11-15: 8, 18-25: 3 | 47 | 44.3 |
| **0.1 x base for 30 iterations, then 0.36; `dnn.cross_first = 5`** | **377 (-30%)** | 5: 5, 10: 25, 11-17: 4, 22-24: 2 | 42 | 43.0 |

First crossing (`-k 1`) against the warm-call beta factor w: 0.5: 467; 0.36: 414; 0.25: 373;
0.15: 332; 0.1: 306; 0.07: 298.

- **Where the iterations go**: a warm child's first certificate is 5,000-14,000 below the
  cutoff (the parent was 382 short), closes a fifth to a third of that per iteration at
  first and 300-700 units per iteration near the cutoff; it crosses at iteration 5-16. Production then waits for the next multiple of 10: 21 children fathom at
  10 (they crossed at 5-10) and the 12 that crossed at 11-16 wait until 20, because the
  crossing prediction needs two certified values and the first check has one. 537
  iterations against 414 crossed: 1.30x is the ceiling of check timing alone.
- `dnn.cross_first = 5` gives the first check its second point: the 12 children of the
  second check fathom at 11-16. -14% iterations, 47 more certificates (most of them the
  cheap polyhedral part under `dnn.lazycert`), children's CPU 45.3 s against 52.3.
  `cross_first = 3` is no better. wil50's adoption carries over.
- **A lower beta on warm calls shortens the recovery**: the crossing moves from 414 to 306
  iterations at w = 0.1. Kept for the whole solve it costs the nodes that do not fathom:
  the same parent run as a warm child of `26:10,8:0,39:32` (`logs/tho40_d3kid_*.log`) is at
  240,068 after 200 iterations at w = 0.36 and at 239,941 at w = 0.15 (127 lower, 28% more
  shortfall). A schedule fixes that: `DNN_ADMM_WARM_BETA=0.1:30:0.36:1` (0.1 x base for the
  first 30 iterations of a warm call, then 0.36 with a refresh; an environment setting of
  the experiment class, already in `dnn.C`) leaves the give-up node at 240,054 (-14) and,
  with `cross_first = 5`, fathoms the 36 children in **377 iterations against 537 (1.42x)**;
  children's CPU 43.0 s against 52.3 (1.2x: a child has a fixed cost of 0.15-0.25 s, the
  first eigendecomposition and certificate, and about 85 ms per iteration).
- Still on the table after that: 25 of the 36 now fathom at the check at 10, and at w = 0.1
  20 children cross at 6-9 (`-k 1`). A certificate every 2-3 iterations between 5 and 20 on
  warm calls (the check interval is `DNNParams::check_every`, not a key) would take up to
  70 more iterations off (377 towards 306) for about 100 cheap certificates [e].
- Size for the estimate: in tho40ab-32 the fathomed nodes are 70% of the CPU (54,952 of
  78,804 s, `infra/phase1/runs/tho40ab-32/nodes.txt.gz`); 1.2x on them is about 12% of the
  tree's CPU [e], 1.42x in iterations would be 21% if the fixed cost per node came down.
  One parent's children, so the size needs the tuning set (`tho40b_analysis`) before a launch.
- `dnn.dscale = 0.5` on the same node: parent 240,145.6 against 240,132.6 at 600, children
  420 against 414 crossing iterations, 538 against 537 under the production checks. No
  effect on tho40, as expected of a grid instance; no harm either.

## 7. What the patch is and is not

- `admm.patch` (against `latest/dnn.C`, `latest/settings.C` at HEAD 4d3cec57): three keys,
  `dnn.dscale` (exponent p, default 0 = off), `dnn.dscale_until` (iterations of a call that
  run scaled, 0 = all) and `dnn.dscale0` (d0, default 1; measured worse at every other
  value, kept for the record), `Work::set_scaling`, the two-reflector
  transform, the scaled Y/Z step, and a verbose-only `split` line (the arrow's share of
  pres^2 and dY^2 every 50 iterations).
- Not done: `tests/dnn_test.py` and the regression suite were not run on the patched
  binary; the cut phase, child bounds under scaling and `DNN_FUSE` are untested; the
  settings goldens would need regenerating (two new keys).
- The error-bound constants of the certificate are unchanged because the certificate runs
  on the plain transform. The one new step in a certified value's path is the entrywise
  product Zc = (D Zc_t D), whose rounding is a relative 2 ulp per entry of a multiplier that
  the certificate then treats as given (any symmetric Z gives a valid bound), so no margin
  term is needed for it.

## Files

- `base.par`, `par_*.par`: dnntest settings (production keys; variants by name).
- `proto.py`: Python prototype of the PRSM with the scaling; `cont.py` same-node
  continuation; `extrap.py` the certificate-side extrapolation.
- `curve.py` (certified value at given iterations), `reach.py` (iterations to reach the
  baseline's value), `kids.py` (dnntest -K logs).
- `q1.sh` to `q9.sh`, `bb*.sh`: the run queues (`bb_q3.sh`'s whole-tree tai25b run was
  stopped unfinished: `logs/bb_tai25b_full_base_killed.log`); `tho.py`, `bbsum.py` their
  summaries; `tho40_*.par`, `bb_*.par` their settings. `logs/` dnntest and bnbdrv logs, `plogs/`
  prototype logs. `logs/tai50b_n3_*.log` are placeholders (those two runs were dropped).
- `admm.patch`: `git diff` of the worktree `/home/user/wt-admm`.
