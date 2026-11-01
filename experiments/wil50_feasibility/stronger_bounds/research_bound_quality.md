# Certified bound quality at wil50 depths 3-5: what to try next, in what order

Research memo, Nov 1 2026. Question: beyond the regional Boolean-triangle
cuts (README.md steps 2-5), which bound-strengthening ideas could add
significant *certified* strength at wil50's depth-3 to depth-5 nodes, and
in what order should they be tried? Bound quality is the lever: +130 at the
median depth-3 clear miss (0.27% of the cutoff 48,814.05) moves fathoming
up one level and divides the tree by about 47 (`../NOTES.md` section 5);
the cuts give +95 at the median, three quarters of a level. Nothing was
run for this memo (four cores are on an A/B); every number below is read
off the existing logs, the code, or the literature, and estimates are
labelled as such.

## 0. Summary

1. The existing logs say the regional-triangle family is **not exhausted and
   not at its dual optimum**: on the nine clear misses the certified bound
   was still rising +4 to +27 per round (median +12) at rounds 7-8; two to
   three rounds per node were lost to Polyak overshoots (-100 to -390); and
   the primal value of the Lagrangian problem sits 5-28 above the certified
   value at 700-iteration rounds, **25-90 at the 300-iteration rounds used at
   depth 4 and in the search**, 87-108 at 150-iteration rounds. The depth-4
   "gains of +64 to +105" of README step 4 are measured against the give-up
   point; against the children's *converged* uncut value they are +6 to +70
   (table in section 1.2). So before any new inequality, the first thing to
   establish is the Lagrangian ceiling of the present family: the protocol
   in section 4 costs one core for an afternoon and decides between "dual
   method" and "new inequalities".
2. Ranked candidates (section 2): (1) a better dual method and inner
   convergence for the present cuts, estimated +40 to +90 over today's +95
   at the median clear miss (a full level if the upper half holds), low
   effort, low risk; (2) 5-point count inequalities (pentagonal /
   Jünger-Kaibel box inequalities on the regional counts), estimated +10 to
   +40, 2-3 days; (3) mixed-region triangles, rectangles and single
   rows/columns, +5 to +20, 1-2 days; (4) cuts inside the certificate's
   polyhedron, +5 to +30, mostly a measurement tool; (5) projected
   triple-consistency cuts, +10 to +50 but uncertain and expensive to
   separate; (6) the regional second-order moment block, +30 to +80 but
   3-6 weeks and a new certificate; it implies the triangles and *not* the
   5-point count inequalities (Laurent 2003 / Grigoriev 2001), so (2) is not
   subsumed by it.
3. Do not bother (section 5): facility-side regions (wil50's flow matrix is
   random 0-9, no decomposition identity), elimination from the certificate
   (measured: 0 of 2,209 on five near misses, `../tuning/README.md` section
   2), atomic triangles, RLT2/RLT3, cheaper relaxations, more symmetry, the
   cutoff half-space on its own, Lagrangian-DNN/Newton bracketing, and
   hypermetrics of support 7 or more before the support-5 ones are measured.

## 1. What the existing logs say (new readings)

### 1.1 The rounds had not plateaued, and the step rule wastes rounds

Gain per round (certified, units of the objective) in the step-2 runs
(`misses/polyak_*.log`, 700 iterations per round, 1,000 new cuts per round,
Polyak step toward best + 100):

| node | short (converged, uncut) | gains per round 1..8 | best | still rising at the end? |
|---|---|---|---|---|
| 1 | 405 | +29 +34 -103 +84 +43 +23 +12 -124 | +121 | +12/round before the round-8 overshoot |
| 2 | 202 | +9 +42 -45 +50 +24 +2 +7 -126 | +89 | +7/round |
| 3 | 241 | +22 +21 -337 +201 +97 +41 +20 +4 | +69 | +20, +4 |
| 4 | 228 | +3 +51 +6 +28 -20 +30 -9 +14 | +103 | +14 |
| 5 | 316 | +9 +43 -33 +44 +24 +4 +5 -24 | +96 | +5 |
| 7 | 203 | +53 +25 -249 +149 +73 +36 +13 -157 | +101 | +13 |
| 8 | 216 | +27 +21 -389 +213 +112 +58 +27 +8 | +77 | +27, +8 |
| 10 | 158 | +23 +25 -238 +135 +76 +39 +15 -50 | +75 | +15 |
| 13 | 182 | +24 +30 -118 +89 +40 +23 +7 -103 | +95 | +7 |
| 14 | 84 | +26 +19 +11 +10 +8 +2 +3 +3 | +82 | +3 (1.75 short) |

Two things follow. The eight-round value is not the family's value: a
geometric tail on the last gains (ratio 0.7-0.8, i.e. 3-5x the last gain)
puts the plateau of the *current* method +35 to +60 higher (estimate). And
the step rule is the crudest possible: the Polyak step assumes the dual
function grows linearly along the subgradient, but the measured gain per
unit of step is 0.3-0.5 (a 50-100 target produces +20 to +50), and the
overshoot rounds cost 2-3 rounds each on nodes 3, 7, 8, 10, 13.

A second artifact, visible on every fathoming node: the target is
`min(cutoff, best + delta)`, so once the bound is within 0.01 of the cutoff
the step `alpha = theta (target - lb) / |g|^2` vanishes. Node 9 sat at
48,814.0x for rounds 3-7 (gains +1, 0, 0, 0, 0), node 16 for rounds 3-7,
depth-4 child 9 (`depth4/cuts_9.log`) for rounds 4-9, each time with the
primal 20-40 above the cutoff and the fathom decided by rounding. That is
4-6 wasted rounds (1,600-2,400 CPU-s) per fathom, and it understates the
family's gain on exactly the near-miss class the search cares about.

### 1.2 The certified value lags the Lagrangian problem's primal, badly at 300 iterations per round

For each run, the primal `<L_mu, Y> + shift_mu` and the certified value at
the best round, and (depth 4) the children's converged uncut value from
`depth4/uncut_*.log`:

| protocol | primal minus certified at the best round |
|---|---|
| uncut, 1,500-2,500 iterations (converged) | -2 to +6 |
| step 2, 700-iteration rounds (depth 3) | 4, 11, 28, 6, 17, 6, 7, 17, 4, 17, 16, 19, 6, 8 (median 11) |
| step 3/4, 300-iteration rounds (depth 3 near misses, depth 4) | 25-90 (median 54) |
| step 3, 150-iteration rounds | 87-108 |

Depth-4 children of clear misses (`depth4/`, 400 uncut iterations then
300-iteration rounds), converged uncut value against the cut phase:

| child | uncut converged | short | give-up (400 it.) | best cut lb | gain vs give-up | **gain vs converged** | primal at best round | primal minus certified |
|---|---|---|---|---|---|---|---|---|
| 1 | 48,751.9 | 62 | 48,696.9 | 48,767.5 | +71 | **+16** | 48,848.5 | 81 |
| 2 | 48,705.3 | 109 | 48,656.9 | 48,732.6 | +76 | **+27** | 48,803.6 | 71 |
| 4 | 48,748.4 | 66 | 48,690.7 | 48,757.3 | +67 | **+9** | 48,845.5 | 88 |
| 6 | 48,715.4 | 99 | 48,679.9 | 48,785.2 | +105 | +70 | 48,839.7 | 55 |
| 7 | 48,707.3 | 107 | 48,676.5 | 48,776.1 | +100 | +69 | 48,830.3 | 54 |
| 11 | 48,747.9 | 66 | 48,690.0 | 48,753.9 | +64 | **+6** | 48,843.3 | 90 |
| 14 | 48,712.6 | 102 | 48,676.6 | 48,753.2 | +77 | +41 | 48,806.9 | 54 |
| 9 (fathomed) | 48,812.6 | 1.5 | 48,782.7 | 48,814.1 | +31 | +1.5 | 48,854.0 | 40 |
| 10 (fathomed) | 48,792.3 | 22 | 48,755.4 | 48,814.1 | +59 | +22 | 48,869.2 | 55 |

The primal is biased upward right after a Polyak step (the step is sized
so that the new Lagrangian evaluated at the old iterate reaches the
target), so the primal is not an estimate of the Lagrangian value L(mu);
but the uncut cold-start pattern (primal 15-20 below the converged value at
400 iterations, certificate 50 below) and the 700-versus-300 comparison
say the certificate is the slow side. Where L(mu) sits in the 55-90 band is
unknown and is the single most valuable number to measure (section 4,
step A): if it is near the primal, children 1, 2, 4, 11 (62-109 short)
already fathom with the present family and the search's cut phase is
losing a level to inner convergence, not to the inequalities. The
150-iteration failure of step 3 (primal 73-87 *above* the cutoff, certificate
20-30 below) is the same effect in its extreme form, and suggests a free
fix: a quiet final round (no new cuts, multipliers frozen, 300-600
iterations) before giving a near miss up.

### 1.3 Where the violations are

Per-region violation counts on node 1's converged uncut iterate (Astra's
`astra/node1.threshold.json`; capacity = free locations in the region):

| region | row<1 | row<2 | row<3 | row<4 | col<1 | col<2 | col<3 | col<4 | col<5 | col<6 | col<7 | col<8 | col<9 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| capacity (of 47) | 10 | 19 | 28 | 37 | 5 | 10 | 15 | 20 | 23 | 28 | 33 | 38 | 42 |
| violated | 138 | 775 | 817 | 123 | 66 | 241 | 614 | 2,395 | 3,716 | 1,576 | 480 | 200 | 65 |
| max violation | .047 | .068 | .078 | .035 | .034 | .055 | .066 | .089 | .091 | .075 | .068 | .055 | .042 |

The five regions with capacity 19-28 (near m/2) carry 83% of the violations
and the largest ones. This is the bisection picture: the objective is
exactly `sum_S cut_F(z^S)` over the 13 thresholds (the memo's identity), a
near-balanced cut carries the most weight, and the DNN's regional marginal
is loosest there. Astra's larger families on the same iterate: rectangles
(770 regions) 21,705 violated, mixed-region triangles (611 observables,
151M inequalities) 389,765 violated, both with the *same* maximum 0.0911,
and the top mixed cuts are near-duplicates of the top same-region ones
(col<4/col<5 variants of the same facility triple at 0.0910 against 0.0911).

The pool after eight rounds: 5,000-7,000 active cuts, 4,400-16,500 still
violated at the final iterate with maximum 0.035-0.075. Neither the pool
(1,000 of about 10,000 violated per round, out of 843,180 candidates at
m = 47) nor the multipliers are saturated; the trajectories above say the
multipliers and the inner solve bind first.

### 1.4 The facility side

`qaplib/wil50.dat`: the first matrix is the 5 x 10 Manhattan distance
(verified entrywise), the second is a symmetric flow matrix with entries
uniform on 0-9 (2,198 of 2,450 off-diagonal entries nonzero, row sums
178-267, no block or grid structure). The regional cuts work on the
location side because the objective is linear in the regional moments
(p^S, W^S) of the 13 thresholds; there is no such identity on the facility
side, and the nug30 `-side 1` success is because nug30's *second* matrix is
the grid. Facility-side regions are in the do-not-bother list.

## 2. Candidates, ranked

Units: certified gain at the median depth-3 clear miss (today +95 after
eight 700-iteration rounds, converged uncut value 185-405 short, median
240; one level = +130). "Cost" is relative to the search's cut phase (12
rounds x 300 iterations, about 2,000 CPU-s at m = 46-47 with four
concurrent workers). All gains are estimates unless marked measured.

| rank | candidate | expected certified gain, depth-3 median | cost per node vs the cut phase | effort | risk | basis |
|---|---|---|---|---|---|---|
| 1 | **Dual method and inner convergence for the present family**: step rule without the cutoff cap and with theta 0.4; a quiet final round; a proximal bundle or level method on the rounds; then multipliers updated inside the ADMM | +40 to +90 (to +135-185); at depth 4 +30 to +60 from inner convergence alone | 1.0x for the step fixes and bundle (same rounds, no overshoots); 1.05-1.1x per iteration for in-ADMM multipliers | 1 day (step fixes, freeze round), 3-5 days (bundle), 1-2 weeks (in-ADMM) | low-medium: the ceiling may be lower than the tail extrapolation | section 1.1-1.2: unplateaued gains, overshoots, 25-90 primal-certificate gap |
| 2 | **5-point count inequalities** on the regional counts (pentagonal / Jünger-Kaibel box inequalities), on the 5 mid regions first | +10 to +40 | +5-10 s separation per round; cut cost 1.7x a triangle | 2-3 days | medium: value only shows once the triangle dual is near its ceiling | section 3.2; max-cut solvers add them after triangles; not implied by M2 |
| 3 | **More regions / mixed regions**: mixed-region triangles, row/column intervals (rectangles), single rows and columns | +5 to +20 | +5-10 s separation per round; same cut cost | 1-2 days | low effort, low value: same max violation, near-duplicate top cuts | Astra's separation counts (section 1.3) |
| 4 | **Cuts in the certificate's polyhedron**: LP over P_GL and the pool for the final Z (once per node), or a subgradient ascent on mu with GL evaluations every k-th certificate | +5 to +30, largely overlapping with 1 (it recovers the inner lag on the polyhedral side) | one LP of ~4.9M variables per node (minutes) or +20-50% per certificate | 3-5 days (expose the GL evaluator; dump L+Z and the pool) | medium | section 3.3 |
| 5 | **Projected triple-consistency cuts** (Astra section 8, RLT2-like) | +10 to +50, uncertain | 1.5-7.7 s LP per triple: 20-50 triples per round = 0.5-2x a round | 3-5 days with HiGHS offline (scipy 1.18 is installed) | medium-high | 10 of 10 separated triples violated at 0.25-0.28 (Astra), never tested in a bound |
| 6 | **Regional M2 block** (second-order moments of z^S, 1,129 x 1,129 at m = 47) on the mid regions | +30 to +80, uncertain | +0.15x per region per iteration (0.6x for 4 regions, 2.0x for 13) plus consensus links | 3-6 weeks, new certificate (trace-aware) | high | implies the triangles, not the 5-point counts; max-cut level-2 experience |
| 7 | **Cutoff-dependent aggregate cuts** (Kim-Kojima successive cutting inequalities on regional counts) | +0 to +20 | one DNN solve per candidate (100-800 s) | 1 week | high | the single-variable case (elimination) measured 0 of 2,209 |
| 8 | Facility-side regions | +0 to +10 | - | 0 (code exists) | - | no identity (section 1.4) |
| - | Elimination from the certificate | 0 (measured) | - | - | - | `../tuning/README.md` section 2 |

A full level needs +130 at the median: candidate 1's upper half alone, or
candidate 1's lower half plus candidate 2. Candidates 5 and 6 are the ones
to reach for if the section-4 protocol shows the linear families capped
below +130.

### 2.1 Why candidate 1 is first

Three separate losses are visible in the logs and each has a cheap remedy:

- **Step size.** Observed gain per unit of Polyak step 0.3-0.5, so the
  linear model over-predicts by 2-3x and the step overshoots. Setting
  `theta = 0.4` (with the halving rule kept), or a level method (target =
  best + delta, but the step projected onto the level set of a cutting-plane
  model of the dual), removes the -100 to -390 rounds. At a fixed eight
  rounds this is +20 to +40 on the nodes that overshot (3, 7, 8, 10, 13),
  estimated from their recovery rounds.
- **Cap at the cutoff.** `target = best + delta` regardless of the cutoff
  (or `cutoff + 5`): the fathoming nodes stop crawling at 48,814.0x. Saves
  4-6 rounds per fathom and gives the true family gain on near misses.
- **Inner lag.** A quiet final round (theta 0, no new cuts, 300-600
  iterations) at the end of a cut phase, or before giving up: the
  certificate catches up with the Lagrangian value of the final multipliers.
  On the step-3 150-iteration runs (primal 48,887-48,901, certificate
  48,784-48,813) this plausibly fathoms all four near misses in the budget
  that failed; at depth 4 it is the +30 to +60 of the table.
- **Dual optimum of the pool.** A proximal bundle method (Fischer, Gruber,
  Rendl, Sotirov 2006 is the precedent for exactly this: dualized cutting
  planes on an SDP, each evaluation one SDP solve) reaches the dual optimum
  of the pool in 10-30 evaluations, monotonically; the current method does
  not reach it in 8-12. The in-ADMM variant (the cuts as a block with their
  own multipliers updated every iteration, as MADAM does for hypermetric
  cuts: slack s = h - GY, a cached sparse factorization of GG' + I, Hrga and
  Povh 2021) reaches the exact value of DNN ∩ pool at ADMM cost; our
  Y = VRV' splitting is not MADAM's, so the cut block would go into the
  polyhedral step as an inexact projection (a few Dykstra/projected-gradient
  passes between the box and the affine cut set), and the certificate uses
  the block's multipliers as mu exactly as now.

Pool size is a free parameter: `-k 5000` (all violated cuts enter at
mu = 0) costs nothing in the ADMM (L_add is dense anyway; 5,000 regional
cuts are 10-25M flops to assemble) and spreads the Polyak step over more
cuts, which should overshoot less; worth one run.

### 2.2 Why candidate 2 is second

The regional triangle type-3 inequality `p_i + p_j + p_k - W_ij - W_ik -
W_jk <= 1` is the case |T| = 3, r = 1 of the integer-count inequality
`(q - r)(q - r - 1) >= 0`, q = sum_{i in T} z^S_i (the memo's "box"
inequalities; Jünger and Kaibel 2001 derive their QAP box inequalities from
the same hypermetric origin and prove many of them facet-defining for the
QAP polytope). In cut-polytope language (s = 1 - 2z, Q_ij = 1 - 2p_i - 2p_j
+ 4W_ij) these are the hypermetric inequalities b'Qb >= 1: |T| = 3 are the
triangles, |T| = 4 adds nothing (CUT_4 = MET_4), |T| = 5 with r = 2 is the
pure pentagonal inequality, and the signed patterns (three +1, two -1) are
its switchings. So "second moments of the counts" and "pentagonal facets"
are the same family, and it is the next facet class of the per-region cut
polytope after the triangles. Two facts fix its place:

- It is **not implied by the regional M2 block**: Laurent (MOR 28(4), 2003)
  shows ⌊n/2⌋ Lasserre rounds do not give CUT_n for odd n, and the witness
  is exactly the function (q - ⌊n/2⌋)(q - ⌊n/2⌋ - 1) (Grigoriev 2001; restated
  in Kurpisz, Leppänen, Mastrolilli, arXiv 1605.03019: no sum of squares of
  degree-⌊n/2⌋ polynomials certifies it). For n = 5 that is (q - 2)(q - 3)
  >= 0 against degree-2 squares, i.e. against M2. So candidate 2 is
  complementary to candidate 6, not a subset of it.
- In the max-cut solvers that have both, pentagonal and heptagonal cuts
  are added only after the triangle violations fall below a threshold
  (BiqBin: r_tri < 0.2, then 300 pentagonal and 200 heptagonal per
  separation; MADAM the same with thresholds 0.2 and 0.4), because until
  then the triangles dominate. Our triangle violations are still 0.05-0.07
  at round 8, so the value of candidate 2 is conditional on candidate 1
  bringing the triangle dual near its ceiling; the free diagnostic is the
  maximum pentagonal violation at a converged iterate (section 3.2).

Separation is exhaustive at this size: C(47,5) = 1.53M facility 5-subsets
per region, 11 sign patterns, 13 regions: 2.2e8 evaluations of a 5 x 5
quadratic form in (p, W), a few seconds in C per round; on the 5 mid regions
2 seconds. Each cut has 10 pair terms and 5 linear terms, 1.7x a triangle in
`L_add`.

Heptagonal (|T| = 7) separation is C(47,7) = 63M subsets x 36 patterns x 13
regions: heuristic only (BiqBin solves it as a small QAP by simulated
annealing). Defer until the support-5 family is measured.

### 2.3 Candidate 3

Mixed-region triangles (p^S_i, p^T_j, p^U_k, W^{ST}_ij, ...) are valid for
any three 0/1 variables and are a one-struct change (`Cut` gets a region
per vertex; `aggregate` already produces every (S, T) cross moment as
`sum_{a in S, b in T} Y_(ai),(bj)`). The separation space is 4 x C(611,3) =
151M at m = 47 (Astra's count), 5-10 s per round in C. Rectangles (all
row and column intervals: 10 + 45 = 55 regions instead of 13) and single
rows/columns (15 regions, capacity 5-10) are a change to `grid_regions`
only. Expected +5 to +20 in total: the maximum violation of all three
families equals the threshold family's, the top cuts are near-duplicates,
and the small regions carry few violations (section 1.3). Cheap enough to
run once on the 14 misses; not where a level comes from.

## 3. The specific questions

### 3.1 More of the same family (question 1)

Covered by candidates 2 and 3 and section 1.3-1.4. The Manhattan reasoning
cuts both ways: the thresholds are the regions that the objective is linear
in, so cuts on them act on the objective directly, which is why they beat
the atomic triangles (+0.2 on nug30); any other region acts only through
the DNN's coupling, so rectangles, single rows/columns, and facility-side
sets should be expected to behave more like the atomic family. Cardinality
indexed inequalities (q <= k_S and its products with z_i) are implied by
Y >= 0 and the lifted row sums (`sum_j W^S_ij = k_S p^S_i` follows from the
range(V) column structure), and the regional first-moment matrix
[1 p'; p W] is a compression of Y and already PSD, as the memo notes.

### 3.2 A better dual method (question 2)

Section 2.1 has the methods. On "pool or multipliers": at round 8 the pool
holds 5,000-7,000 of 843,180 candidates with 4,400-16,500 still violated at
0.035-0.075, so the pool is not exhausted; but the certified bound was
rising +12 per round at the median, 2-3 rounds per node were overshoots,
and the certificate trails the Lagrangian primal by 5-28 even at 700
iterations per round. The multipliers and the inner solve bind first, the
pool second. The protocol of section 4 measures the three pieces
separately on one node.

The first experiment on the 14 misses (ids = line numbers of
`../tuning/giveup/all32.txt`), scorable against the step-2 table:

1. Zero code: the step-2 protocol with `-k 5000` on the nine clear misses,
   and `-t 0.4` on nodes 3, 7, 8, 10, 13 (the overshooters). Score: best
   certified value after 8 rounds against +95 median; count of rounds with
   a drop.
2. Small change to `latest/dnncut.C` (about 15 lines): `DNNCutParams`
   gains `cap_at_cutoff` (default 1; 0 means `target = best_lb + delta`),
   and `freeze_after`/`freeze_iters` (after that round: `theta_cur = 0`,
   `kcuts = 0`, one round of `freeze_iters` iterations). `dnncuts` gets
   `-nocap` and `-freeze R:N`. Run the four near misses (6, 9, 15, 16) and
   node 14 from the give-up point (step-3 protocol, 300-iteration rounds)
   with `-nocap`: score the fathom round and CPU against step 3's 5-10
   rounds; node 14 with `-nocap -delta 50`: does it fathom (it ended 1.75
   short)?
3. The quiet round: step 3's 150-iteration protocol plus `-freeze 8:600`
   on 6, 9, 15, 16: fathoms in the budget that failed?
4. Then the bundle: a proximal bundle on the rounds (each round one
   evaluation of the dual function at the current mu, subgradient = the
   violations of the final iterate; about 200 lines), compared with the
   Polyak run at equal rounds on all 14.

### 3.3 Cuts in the certificate's polyhedron (question 3)

The certificate is `min_{Y in P} <L + Z, Y> - (m+1) lambda_max(V'ZV)` with
P the Gilmore-Lawler structure (per column block: nonnegative, gangster
zeros, row and column sums equal to x_ij; x doubly stochastic). Any
polyhedron containing every permutation lift is valid, so `P ∩ {G Y <= h}`
with the regional cuts is valid too, and by LP duality

    min_{Y in P, GY <= h} <L + Z, Y> = max_{mu >= 0} [ min_{Y in P} <L + Z + G'mu, Y> - h'mu ],

i.e. putting the cuts into P equals optimizing the cut multipliers *for the
polyhedral part* at the current Z. For the ADMM's Z (tuned to L_mu) the
value is sandwiched: L(mu) <= c(Z, mu') <= L(mu'). So this never certifies
less than the present certificate, and it recovers whatever the inner lag
loses on the polyhedral side without touching Z. What it costs:

- Exact LP: variables = the compatible entries of Y (about m^4 = 4.9M at
  m = 47), 2m^3 = 208k block-sum equalities, 2m doubly-stochastic rows, and
  the pool (5-10k cuts with 1,200-5,300 nonzeros each, 10-50M nonzeros).
  Minutes with a barrier solver, so not every 10 iterations (an ADMM
  iteration is 0.18 s alone at m = 47); affordable once per node at the
  end of a cut phase as a last-chance certificate on near misses.
- Subgradient ascent on mu with GL evaluations: the polyhedral part of a
  certificate is cheap (the 56 ms per certificate in the m = 47 profiles,
  `../../eigensolver/node2_*.log`, is dominated by the 2,117 x 2,117
  Cholesky; the GL part is of order 10-20 ms with reduce-mode LAPs), and the
  lambda_max term is fixed for fixed Z, so 20-50 ascent steps per certificate
  cost 0.3-1 s, every fifth certificate. The subgradient `G Y* - h` comes
  from the GL minimizer (a vertex: the outer assignment and one inner
  assignment per column), which is cheap to aggregate per region. The
  resulting mu' can also seed the next round's multipliers, which turns
  this into the "dual ascent on the cuts with the DNN block fixed" of the
  question, updated every 50 iterations instead of every 300-700.

Expected value: +5 to +30 at the median (estimate: bounded above by the
inner-lag plus multiplier loss for a fixed Z, which the protocol measures
in step D). Effort 3-5 days: expose `certify`'s polyhedral part on an
external objective, dump `L + Z_c`, the pool and mu at the end of a run.

### 3.4 Hypermetric and other facets on the regional variables (question 4)

The regional (p^S, W^S) are the first and second moments of a 0/1 vector
z^S of length m with `sum z^S = k_S`, so every valid inequality of the
Boolean quadric polytope QP_m (equivalently of CUT_{m+1} by the covariance
map) applies, and the cardinality constraint adds the k-cardinality-cut
(graph bisection) facets. In order of what is known to add strength on
correlation-matrix relaxations:

- Triangles (support 3): the workhorse; in hand.
- Pentagonal / 5-clique (support 5) and their switchings: candidate 2.
  Heptagonal (support 7): defer.
- Odd-cycle inequalities on a sparse support graph: for dense weights (our
  F is dense) the cycle inequalities beyond triangles are the hypermetrics
  above; nothing separate.
- Clique inequalities with larger |T| and r = 2, 3 (|T| = 6, 7): the
  Jünger-Kaibel 4-box family; heuristic separation; after support 5.
- Bisection-specific facets: Sotirov (Annals of OR 2018) shows a large
  subset of the facet-defining inequalities of the bisection polytope are
  redundant for the strongest SDP relaxation (which the DNN's regional
  marginal is a compression of); nothing to add there beyond the
  hypermetrics.

The memo's exact witness (K_{3,3}, DNN 9 against optimum 10) is closed by a
single type-3 triangle; it says nothing about support 5. The right witness
for candidate 2 is the Grigoriev function on five facilities: a regional
marginal with p = 1/2 and W such that (q - 2)(q - 3) has negative
expectation while every triangle holds; the same construction as the memo's
(Q with -1/2 within groups) on five variables gives one.

### 3.5 Lifting beyond DNN on the regional variables only (question 5)

Sizes at m = 47 (depth 3), against the 2,117 x 2,117 block we pay for:

| object | size | eigendecomposition cost, relative | notes |
|---|---|---|---|
| present facially reduced DNN block | 2,117 | 1 | (m-1)^2 + 1 |
| one regional M2 block, 1 + m + C(m,2) | 1,129 (1,081 after the cardinality facial reduction) | 0.15 (0.13) | moments of z^S to degree 4; trace 1 + k_S + C(k_S, 2) known |
| four mid regions | 4 x 1,129 | 0.6 | col<4, col<5, col<6, row<3 at this node |
| all 13 | 13 x 1,129 | 2.0 | the memo's 1.95 at the root |
| Lovász-Schrijver N+ on the 13 membership blocks | nothing new | 0 | the conditional rows `sum_{a in S} Y_(ai),.` are already in Y |
| "sum of region membership" variables with their own PSD block | [1 p'; p W] per region | 0 | a compression of Y, already PSD (memo section 4) |
| local marginal polytope on 7 regional variables (128 weights) | 128 columns per block | negligible | Sherali-Adams level 7 on a chosen 7-set; linear, goes through L_add or the Y-step |

Per iteration the M2 route is +0.6x to +2.0x on the eigendecompositions,
plus the consensus links (u_i = p_i, u_ij = W_ij: 1,129 linear links per
region, an averaging step in the Y-projection), plus a two-block splitting
that keeps convergence (the memo's warning about Gauss-Seidel on three
blocks stands). What it buys beyond the linear families: every degree-4
consistency of the per-region bisection (all triangles at once, plus the
4-point local polytope with a PSD condition), which on max-cut at n = 50
closes most of what triangles leave; what it does not buy: the 5-point
count inequalities (section 2.2) and any cross-region coupling. Estimate
+30 to +80 at the median, after candidates 1-2; 3-6 weeks including a
trace-aware certificate (the memo section 7: subtract eps_b tau_b per
block). Not first.

The 7-variable local marginal polytope blocks (the memo's 128-weight
option) are the cheap way to test "is higher-order consistency on a few
chosen variables worth anything": pick the 7 facilities most involved in
the violated cuts of a mid region, add the 128 configuration weights as an
LP block linked to (p, W) in the Y-step. If the gain on the 14 misses is
under +10, the M2 route is unlikely to pay either.

### 3.6 Literature 2021-2026 (question 6)

Nothing published computes a relaxation stronger than the DNN on QAP at
n >= 40; the recent work is on solving the same relaxation faster, and on
cutting planes for DNN relaxations of other binary or nonconvex QPs:

- Kim and Kojima, *Generating cutting inequalities successively for
  quadratic optimization problems in binary variables*, arXiv 2107.08665
  (2021). A candidate integer inequality g'x >= alpha + 1 is proved valid
  for the *optimal* set by solving the conic relaxation of the subproblem
  with g'x <= alpha added: if its bound exceeds the incumbent, the cut is
  valid. On 60 BiqMac QUBOs (n = 100-250) ten rounds reached the optimum
  on 70%. For us: candidate 7; each proof costs a DNN solve, and the
  single-variable case (child bounds) measured zero eliminations at depth
  3 (`../tuning/README.md` section 2, reduced costs flat at convergence).
  A cheap version through the certificate's GL structure (the outer
  assignment with one knapsack side constraint, solvable by bisection on a
  scalar multiplier) is the only form worth a look later.
- Hrga and Povh, *MADAM*, COAP 80 (2021), arXiv 2010.07839: ADMM for the
  max-cut SDP with triangle, pentagonal and heptagonal inequalities as a
  block (slack, cached sparse Cholesky of BB' + I, projection onto the
  orthant), a two-block convergence proof via the orthogonality of diag
  and B, and a safe bound by shifting the dual slack by lambda_min. Add
  pentagonals when r_tri < 0.2, heptagonals when r_pent < 0.4, +200 per
  non-pruned node. Gusmeroli et al., *BiqBin*, ACM TOMS (2022), arXiv
  2009.06240: the same cuts with a bundle method (3 to 15 bundle iterations
  per separation), separation of pentagonals/heptagonals as a small QAP
  by simulated annealing, 10n triangles + 300 pentagonal + 200 heptagonal
  per separation. These are the precedents for candidates 1 (in-ADMM
  block, bundle) and 2.
- Qu, Sun, Xu, *Progressive bound strengthening via doubly nonnegative
  cutting planes for nonconvex QPs*, arXiv 2510.02948 (Oct 2025): DNN
  bounds plus valid cuts from a linear SDP at KKT points of the nonconvex
  QP; continuous QPs, cuts of no fixed form; not applicable to a
  permutation polytope as it stands.
- Hou, Tang, Toh, *RiNNAL+*, arXiv 2507.13776 (2025, MPC 2026) and
  *RNNAL*, arXiv 2502.13849 (2025): low-rank Riemannian ALM for DNN / SDP-RLT
  relaxations of mixed-binary QPs, with the proof that the DNN relaxation of
  dimension n + l + 1 equals SDP-RLT of dimension n + 1; QAP (nug21-25) as a
  test case for the constraint-relaxation strategy. Solver speed on the
  same relaxation; no cuts.
- Campos, Misener, Parpas, *Partial Lasserre relaxation for sparse
  max-cut*, Optim. Eng. 24 (2023): level-2 moment blocks on selected
  cliques only; the precedent for candidate 6's "mid regions only".
- Fischer, Gruber, Rendl, Sotirov, Math. Prog. 105 (2006): bundle method
  with dynamically added, dualized inequalities on SDP relaxations of
  max-cut and equipartition; the precedent for candidate 1's bundle.
- Jünger and Kaibel, Math. Prog. 91 (2001): box inequalities for the QAP
  polytope, from hypermetric inequalities on the count of assignments in a
  row-set x column-set box; candidate 2 in QAP-polytope language.
- Laurent, MOR 28(4) 871-883 (2003) and MOR 28(3) 470-496 (2003);
  Grigoriev, *Complexity of positivstellensatz proofs for the knapsack*
  (2001); Kurpisz, Leppänen, Mastrolilli, arXiv 1605.03019 (2016): the
  degree-4 SOS does not certify the 5-point count inequality (section 2.2).
- Sotirov, *Graph bisection revisited*, Annals of OR 265 (2018): the
  order-n SDP equivalent to the strongest bisection relaxation; many
  bisection facets redundant for it.
- Adams, Guignard, Hahn, Hightower (2007) RLT2; Hahn, Roth, Saltzman,
  Guignard (2013) RLT3 at 10 TB for n = 50: as in `../solvers/research_relaxations.md`.

(optimization-online.org is not reachable from this session; the Laurent
and Fischer et al. statements are from the published versions as cited in
the arXiv papers above and from memory of the originals, and should be
checked against the journal text before being quoted.)

## 4. Protocol: the Lagrangian ceiling of the present family on one node

Purpose: split the gap between today's +95 and the family's value into
three measured pieces (inner lag, multiplier suboptimality, pool
incompleteness), so the next investment is chosen on numbers. One core;
nodes 14 (84 short at convergence, 1.75 short after step 2) and 2 (202
short, overshoot-prone) at depth 3, and depth-4 child 1 of `depth4/nodes.txt`
(62 short at convergence, cut phase +16 over converged with primal 81 above
the certificate). Code change: the `freeze_after`/`freeze_iters` and
`cap_at_cutoff` options of section 3.2 (about 15 lines in `dnncut.{h,C}`,
`-freeze R:N -nocap` in `dnncuts`).

A. **Inner lag, L(mu_8).** Repeat the step-2 run (`-I 1500 -i 700 -k 1000
   -delta 100 -grid 5:10`) with `-freeze 8:2000 -nocap`: after round 8 the
   multipliers freeze, no cuts are added, 2,000 iterations run. Accept
   when the final primal and certified value agree within 3 (as the uncut
   runs do). Delta_inner = L(mu_8) - best_8. Prediction: +5 to +25 at depth
   3 (node 2's gap is 11 at the best round), +20 to +60 at depth 4 with the
   300-iteration protocol (`-I 400 -i 300 -delta 50 -freeze 12:2000`).
   Cost: 12,000 iterations, about 40 min per node alone at 0.18-0.2 s per
   iteration (2 h with four concurrent).
B. **Multiplier loss, L*(pool_8).** From the same pool (a dump/reload of
   the pool and mu, or simply the same run continued with `-k 0`): 20
   rounds of 300 iterations with `theta 0.3, delta 30, -nocap` (no
   overshoot expected at that step), then a freeze round of 2,000.
   Delta_mult = L*(pool_8) - L(mu_8). A proximal bundle in place of the
   conservative Polyak step is the better tool once written; the Polyak
   version is available today.
C. **Pool loss, L*(family).** Alternate: re-separate at the converged
   iterate (`-k 3000`), 10 conservative rounds, freeze; repeat until the
   maximum violation at the converged iterate is below 0.01 or a cycle
   gains under 2. Delta_pool = L*(family) - L*(pool_8). Cost: 30-60 rounds,
   2-4 h per node alone.
D. (optional) **Certificate ceiling for the final Z.** Dump `L + Z_c`, the
   pool and mu at the end of C; solve the LP of section 3.3 with HiGHS
   (scipy 1.18) or run the GL subgradient ascent. Delta_cert = LP value -
   certified value: how much a certificate-side fix alone would buy.

Decision rule (the median of the two depth-3 nodes):

- L*(family) - best_8 >= +40 (ceiling >= +135 on node 2): the family is a
  full level and the dual method is the lever; implement candidate 1 fully
  (bundle, then in-ADMM multipliers), candidate 2 after.
- Delta_inner >= 60% of the total: fix the inner convergence and the
  certificate following first (quiet rounds, candidate 4's ascent); it is
  free quality in the search's 300-iteration regime.
- L*(family) - best_8 < +20: the family is capped near +100-115; go to
  candidate 2, then 5 and 6, and run this protocol again with the pentagonal
  family in the pool.

Scoring the follow-ups on the 14 misses: certified best bound per node
against the step-2 table (median +95, fathoms on 6, 9, 15, 16, node 14 at
1.75 short), under the same iteration budget, every value certified.

## 5. Do not bother (and why)

- **Facility-side regions on wil50.** The flow matrix is random 0-9 with no
  structure (section 1.4); the regional identity that makes the location
  side work has no counterpart. nug30 `-side 1` worked because its second
  matrix is the grid. Keep the option for sko/tho/nug (grid on one side).
- **Elimination from the certificate (reduced-cost fixing).** Measured on
  the five depth-3 near misses: child bounds equal lb for 90% of the
  children, largest excess 2.8-6.5 against 43-84 needed, 0 of 2,209
  eliminated (`../tuning/README.md` section 2). The cuts raise lb by 60-120
  but do not change the flatness of the outer assignment problem, so the
  cut-strengthened certificate will not eliminate either; and the SCIT
  generalization (candidate 7) costs a DNN solve per candidate.
- **Atomic triangles, RLT2 as an LP, RLT3.** +0.2 on nug30; weaker than
  DNN in every computed value at n >= 25 and nobody has solved the LP at
  n = 30; 10 TB at n = 50 (`../solvers/research_relaxations.md`).
- **Cheaper relaxations and the pure SDP on the Manhattan matrix.** All
  2-10% weaker; alpha J - E[X'DX] >= 0 is implied by the lifted PSD block
  (memo section 4).
- **The objective cutoff as a constraint on its own.** min <L, Y> subject
  to <L, Y> <= cutoff has the same minimum unless infeasible; it only helps
  as a side condition when proving other cuts (candidate 7).
- **More symmetry, region branching, eigensolver speedups, Anderson /
  Halpern / momentum.** Not bound quality (the last three matter for the
  inner lag only through iteration counts, which the quiet-round fix
  addresses at fixed budget).
- **Lagrangian-DNN / Newton bracketing.** Same relaxation, same value.
- **Heptagonal and larger hypermetrics, 4-box inequalities.** Heuristic
  separation, and the max-cut solvers add them only after the support-5
  family is exhausted; measure support 5 first.
- **Mixed-region M2 or a full second-order assignment moment matrix.**
  3M x 3M before reductions (memo section 6); the regional M2 is the only
  affordable lift and it is candidate 6.

## 6. Appendix: the numbers used

Per-run primal and certified values (from `misses/`, `cost/`, `depth4/`):

| run | round 0 lb / primal | best round: lb / primal | gap |
|---|---|---|---|
| misses/polyak_1 | 48,409.3 / 48,406.9 | r7 48,530.2 / 48,534.9 | 4.6 |
| misses/polyak_2 | 48,611.8 / 48,600.9 | r7 48,700.5 / 48,711.3 | 10.8 |
| misses/polyak_3 | 48,573.3 / 48,573.4 | r8 48,642.1 / 48,669.6 | 27.5 |
| misses/polyak_8 | 48,598.1 / 48,597.3 | r8 48,674.8 / 48,691.6 | 16.9 |
| misses/polyak_14 | 48,729.9 / 48,728.0 | r8 48,812.3 / 48,831.0 | 18.7 |
| cost/g400_9 (300-it rounds) | 48,746.1 / 48,723.8 | r5 48,814.1 / 48,839.2 | 25.1 |
| cost/g400_16 | 48,725.0 / 48,714.3 | r9 48,814.0 / 48,857.7 | 43.6 |
| cost/g400s_6 (150-it rounds) | 48,731.9 / 48,718.6 | r10 48,795.5 / 48,893.2 | 97.8 |
| cost/g400s_16 | 48,725.0 / 48,714.3 | r20 48,783.6 / 48,891.1 | 107.5 |
| depth4/cuts_1 | 48,696.9 / 48,731.1 | r12 48,767.5 / 48,848.5 | 81.1 |
| depth4/cuts_11 | 48,690.0 / 48,727.0 | r9 48,753.9 / 48,843.3 | 89.5 |
| depth4/uncut_1 (1,500 it.) | - | 48,751.9 / 48,758.1 | 6.2 |

Per-iteration cost at m = 47 on the 2.1 GHz Xeon: 0.17-0.18 s alone
(`../../eigensolver/node2_*.log`, 2,500 iterations in 425-443 s, of which
eig 282 s, certify 14 s = 56 ms per certificate), 0.52-0.6 s with four
concurrent workers (`depth4/uncut_*.log`). A 300-iteration round is 160-180
s concurrent; the search's cut phase (12 rounds) about 2,000 s.

Cut family sizes at m = 47: same-region triangles 13 x 4 x C(47,3) =
843,180; rectangles (Astra, 770 regions) about 50M; mixed-region triangles
151M; pentagonal (support 5, 11 patterns) 13 x 11 x C(47,5) = 219M;
heptagonal 13 x 36 x C(47,7) = 2.9e10.
