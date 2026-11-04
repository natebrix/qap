# The facility side and disjunctive strengthening for wil50: a negative result

Agent 3 (facility side), Nov 4, 05:26-09:50 UTC, core 3. Saved by the
coordinator from the agent's hand-back message (the agent could not write
report files). `NOTES.md` in this directory is the agent's running log.

## Summary

Nothing from the facility side, and no disjunctive or lift-and-project
construction, adds certified bound on top of the converged regional-triangle
phase at the nodes tested.

- **The flow matrix is random.** wil50's flow matrix is, statistically, a
  uniform random 0-9 matrix: no clusters, no spectral structure.
- **The facility side is already in the families tried.** The objective is the
  sum, over the 13 grid thresholds, of the flow-weighted cut of the
  facilities' threshold memberships, so every regional triangle is already a
  bisection inequality of the flow graph.
- **The DNN cheats only on the 13 threshold memberships.** On dumped iterates
  of node 13 (give-up point, uncut converged, after the triangle phase) and of
  child d2 (after its triangle phase), no facility-side count, region, parity,
  dispersion or star inequality is violated. Facility-region triangles, with
  the region optimised per location triple, reach only 0.006 uncut and 0.0001
  after the triangles.
- **The one violated facility-side family gives nothing.** The 3-local
  consistency of the heaviest flow triples (projected RLT2) is violated by up
  to 0.114 after the triangles; in a Polyak pool on top of the frozen triangle
  term the best multiplier is 0 and the gain is 0. The location-side analogue
  (compact L-shaped triples, violated by up to 0.104) gives at most +1.5,
  within the drift from extra iterations.
- **Every disjunction is capped by its worst term**, the child of the cell the
  DNN already favours most. At node 10 (166 short), to convergence: any
  disjunction on its best lines at most +42 to +59; one Sherali-Adams / N+
  level on the bottleneck variable +19 to +24; a count split +2.4 to +2.9.
  Each costs DNN solves per term, so branching dominates all of them.
- **The certificate's child bounds stay flat with cuts**: node 13 after the
  phase 0 of 2,209 children fathomed (largest excess 43.5 against 68.4
  needed); d2 (29.8 short) 14 of 2,116.
- **Recommendation: close both the facility side and disjunctive cuts.** The
  DNN's remaining inconsistency is in the objective's own variables: the
  threshold memberships and how they couple across thresholds.

**Conventions** (checked by having `scripts/ynode.py` reproduce the logged
primal value from a dumped Y to within 1e-6): "location" is the first QAP
index, the 5 x 10 grid, cell r*10+c, with A the Manhattan distance; "facility"
is the second index, the flows B; a fix a:b puts facility b on location a; the
lifted index of x_(a,j) is 1 + a*m + j. The repo sometimes calls the first
index "facility"; this report never does.

## (a) What the facility side and disjunctions offer, and why

### A.1 The facility side is already inside every regional family

For a permutation, cost = sum_{j,l} B_jl d(loc j, loc l). The Manhattan
distance is d = sum_S [j and l on different sides of S] over the 13 row and
column thresholds S. Hence objective = shift + linear + sum_S 2 cut_B(z^S) and
Delta_jl := E d(j,l) = sum_S (p_j^S + p_l^S - 2 W_jl^S), where z^S_j =
[facility j lies in S]. Each threshold is a cardinality-constrained bisection
of the flow graph B; the regional triangles, the 5-point family and the
local-7 family are all bisection-polytope inequalities of the flow graph,
"location side" only in the choice of the 13 partitions.

Facility regions can also be made to speak the objective's language: take a
facility set T and w^T_a = [location a holds a facility of T]; for c large
enough B + c(J - I) lies in the interior of the cut cone, so B + c(J - I) =
sum_T w_T delta_T with w_T >= 0, and objective + constant = sum_T w_T
cut_A(w^T), linear in each region's moments. A generic B needs up to C(47,2)
= 1,081 such regions with real weights (Caratheodory) against the location
side's 13 of weight 1; and as A.3 shows, none of these regions has a violated
triangle.

### A.2 wil50's flow matrix has no structure (`logs/flow_analysis.log`)

| statistic | B | random symmetric 0-9 matrix |
|---|---|---|
| top eigenvalues | 223.8, 34.4, 31.3 | 227.4, 31.7, 29.7 |
| top centred eigenvalue (semicircle edge 41.3) | 39.3 | 36.9 |
| triangles with all flows >= 9 / >= 8 | 32 / 168 | 29 / 165 expected |
| mean flow in the densest 4/6/10/15-sets | 9.00 / 8.13 / 7.00 / 6.30 | 8.83 / 8.00 / 6.93 / 6.16 |

The mean row correlation is 0.000. The only structure is first order: flow
degrees (row sums) run from 168 to 267, and the DNN already uses it (at node
13 the correlation between a facility's flow degree and the expected
eccentricity of its location is -0.90; light facilities are more decided,
-0.48). Facility 45 is the lightest (row sum 178) and is the xbar line at all
three lines-study parents (`tuning/README.md` calls it "location 45"; "col"
in the lines code is the second index); the DNN sends it to the corner, x_(49,45)
= 0.47 at node 10. Star inequalities price exactly this first-order
structure, and the DNN implies them (A.4).

### A.3 Where the DNN cheats: dumped iterates of node 13 (21:28,8:41,28:19, m = 47)

Three iterates: the give-up point (400 iterations, `runs/n13/tri.Y0.bin`; not
PSD, smallest eigenvalue of Y -0.0105, so large-support counts there reflect
infeasibility), uncut converged (1,500 iterations, the grid agent's dump), and
after the triangle phase (8 Polyak rounds plus a 1,500-iteration quiet round,
lb 48,745.63, `runs/n13/quiet.Y.bin`). Violations in triangle units, positive
= violated; logs `logs/*_n13_{giveup400,uncut1500,tri}.log`.

| family (each valid for every permutation) | give-up | uncut | after triangles |
|---|---|---|---|
| threshold triangles: max violation / count | 0.03-0.10 / 129-9,251 | 0.03-0.08 / 41-2,530 | 0.015-0.054 / 45-438 |
| signed boxes on thresholds, support <= 15 | infeasibility | 0.03-0.06 | 0.02-0.07 |
| facility-region triangles, T optimised per location triple (QUBO) | 0.008 | 0.006 | 0.0001 |
| facility-region triangles, random T of size 10-23 | < 0 | < 0 | < 0 |
| boxes c = &#124;pi(Q) cap P&#124; with P and Q both free | - | violated only when P is a threshold | none (cancellation artefacts only) |
| parity: min over 16,215 triples of Delta_jk + Delta_jl + Delta_kl (valid >= 4) | 4.81 | 4.79 | 4.87 |
| checkerboard-region triangles | none | none | none |
| dispersion, k = 4..7 (sum of Delta minus the grid minimum) | 2.8-12.3 | 2.7-12.0 | 2.9-11.0 |
| triple distance hull (LP over the 520 achievable triples) | metric facets only | metric facets only (0.49) | metric facets only (0.14) |
| stars (hub j, weights B_j.): slack over the conditional Gilmore-Lawler bound | - | 110-264 | - |
| LMP3_fac (L1), heaviest / random facility triples | - | 0.003-0.027 / <= 0.002 | median 0.016, max 0.114 over 60 / <= 0.013 |
| LMP3_loc (L1), compact L-shaped / random location triples | - | - | median 0.044, max 0.104 over 60 / <= 0.006 |

The metric facets violated in the distance hull are sums of the apex
triangles over the 13 thresholds on the same triple, already in the pool's
cone. Child d2 after its own triangle phase (48,784.25, 29.8 short;
`logs/*_d2_tri.log`) gives the same picture: facility-region triangles at most
0.003, no checkerboard violations, parity minimum 5.02, dispersion slack
3.3-11.4, only metric facets in the hull (0.08), threshold triangles still
violated by 0.010-0.040.

Against a feasible solution: a local search found 49,036 for node 13
(`logs/gapdecomp_n13_uncut1500.log`). The gap of 407.6 is a small difference
of large sums over pairs (+6,424 and -5,807); heavy pairs (flows 8-9) sit at a
mean expected distance of 3.46 in the DNN against 3.62 in the solution, light
pairs at 6.16 against 6.01, and each facility's share of the gap is
uncorrelated with its degree (0.00). The cheating is the diffuse bisection
smoothing that the threshold cuts already target; no facility group carries
it.

### A.4 Facility-subset inequalities

- **Stars are implied.** For a hub facility j and weights w >= 0: for x_aj >
  0, Y_(aj),(.)/x_aj is doubly stochastic on (locations other than a) x
  (facilities other than j), so sum_l w_l Delta_jl >= sum_a x_aj LAP_a(w) >=
  min_a LAP_a(w). With w = B_j. the measured slack is 110-264 per hub.
- **Parity and dispersion are never tight.** The grid is bipartite, so three
  pairwise distances have an even sum >= 4, which follows from the
  checkerboard region's sum triangle with pairwise exactness. Five points need
  a sum >= 16 where the triple inequalities force only 14, so k = 5 is not
  implied, but it is slack by 4.3-4.7. The densest flow cluster {5, 9, 20,
  27} (all six flows 9) has expected pairwise distance sum 13.9-14.0 against a
  grid minimum of 8, with 2.26 expected adjacent pairs of a possible 4: heavy
  groups are not placed too compactly.
- **Flow-weighted sub-QAP cuts on heavy cliques** are slack by 140-1,700
  (measured by the grid agent on the same iterate).
- **Junger-Kaibel boxes and QSAP clique/cut inequalities** ((c - b)(c - b - 1)
  >= 0 for the integer count c = |pi(Q) cap P|): with P a threshold and |Q| =
  3 or 5 they are the regional triangles and the 5-point family; larger signed
  Q on the thresholds are violated at the same normalised level as the
  triangles (the literature agent's box gate gives the same 0.07); every
  search that lets P vary freely ends at a threshold; facility regions (P = 3
  locations, Q = T) are never violated.
- **LMP3_fac** (a facility triple's three pair blocks consistent; the
  projected RLT2 family of step 8) is the only facility-side family the DNN
  violates measurably, concentrated on the heaviest triples (about 10x random
  ones) and peaking inside {5, 9, 20, 27}. Test 4 in (c) shows the DNN repairs
  it at no cost: pair mass moves between location pairs with the same
  distance and the same threshold memberships, directions the objective does
  not see (also why the distance hull is consistent apart from the metric
  facets).

### A.5 Disjunctions and lift-and-project

1. **The DNN already is one full level of lift-and-project**: the
   Lovasz-Schrijver N+ relaxation of the assignment polytope (Burer and
   Vandenbussche computed QAP bounds with N+; de Klerk and Sotirov note the
   coincidence). The next level needs m^2 PSD blocks of child size; only
   selected variables or lines are affordable.
2. **The lift-and-project closure of a line equals its worst child.** For a
   line l with children K_c = K cap {x_c = 1}, P_l = conv(union_c K_c) and min
   over P_l of <L, Y> = min_c V(c). Any coarser disjunction on the line ("j in
   region R or not") has a larger hull and a value at most min_c V(c). N+ on a
   single variable x_c is the hull of K_c and K cap {x_c = 0}, value
   min(V(c), V(not c)).
3. **The closure cuts off the parent's optimum Y\*, but its value is set by
   the favoured child.** If Y\* were in P_l, V(parent) >= sum_c x_c V(c). On
   the lines study's 30 lines (3 depth-3 clear misses, every child at 300 cold
   iterations, `logs/lines_xweighted.log`), sum_c x_c V(c) exceeds the
   parent's converged bound on 26 of 30 lines, by up to +124, so the closure
   does cut Y\* off; but its value is the minimum child, which is the line's
   largest-x cell on 22 of 30 lines and the second largest on 7.
4. **Node 10, to convergence with identical settings** (`dnntest -d 3 -L 2 -r
   40 -i 1500`, `runs/n10c/`):

   | node | certified | primal |
   |---|---|---|
   | parent | 48,648.55 (166 short) | 48,651.14 |
   | child 49:45 (bottleneck cell, x = 0.47; lowest child of both best lines at 300 iterations) | 48,690.99 | 48,710.09 |
   | anti-child (x_(49,45) = 0) | 48,672.21 | 48,670.03 |

   One N+ level on x_(49,45) is worth +19 to +24 for two DNN solves; any
   disjunction on row 49 or column 45 at most +42 to +59. Branching on row 49
   fathoms 33 of its 47 children even at 300 cold iterations (median child
   +159).
5. **Certificate child bounds** (reduced costs of the Gilmore-Lawler outer
   assignment, `dnntest -C`):

   | node | children fathomed | largest excess over lb | children within 1 of lb | other |
   |---|---|---|---|---|
   | node 10, uncut | 0 of 2,209 | 7.6 | 2,195 | Spearman 0.21 against the children's real gains (median +117) |
   | node 13, cut certificate | 0 of 2,209 | 43.5 (68.4 needed) | 2,043 | median excess 0.054, 99th percentile 19.4 |
   | d2, cut certificate | 14 of 2,116 | 59.3 (29.8 needed) | 1,921 | mostly on the heaviest facilities' lines, at peripheral locations |

   d2's 14 are a free saving of at most 7 of the 46 children on any one line,
   and they appear only because the node ends close to the cutoff.
6. **Count splits are cheap for the DNN to satisfy.** A count split (c <= b OR
   c >= b + 1, c a signed count of facilities in a threshold region) is the one
   disjunction not dominated by a single line; its degree-2 relaxation is the
   box inequality. It beats 47-way line branching only if both sides gain more
   than 130/log2(47) = 23.4. Two splits on node 13 from the converged triangle
   state (Y and Z warm, 400 iterations, control 48,746.78):

   | split | E c, Var c | c <= b | c >= b+1 |
   |---|---|---|---|
   | most violated box (col<2, 12 signed facilities, b = 8) | 8.58, 0.10 | +4.7 | +2.4 |
   | mid-region box (col<4, 14 facilities, b = 7) | 7.46, 0.18 | +2.9 | +2.6 |

   Larger multipliers only lowered the side bounds.

## (b) Candidate families, ranked (on top of the converged triangle phase, depths 3-5)

[m] measured on these nodes, [e] estimated.

| # | family | gain | separation cost | validity | effort | kill test | status |
|---|---|---|---|---|---|---|---|
| 1 | LMP3_loc on compact L-shaped location triples | [m] +0.3 to +0.5 over the control (pool best +1.5 over round 0) | LP over 97k placements, ~6 s per triple | exact rhs by enumeration; permutation check | 1 day (C++ separator) | 60-cut Polyak pool on the frozen triangle term; kill if < +10 | killed |
| 2 | LMP3_fac on the heaviest facility triples (projected RLT2) | [m] 0 (best multiplier 0; with the gentle step the multipliers return to ~0) | same | same | same | same | killed |
| 3 | signed boxes / clique inequalities of support 6-15 on the thresholds | [e] 0 to +10 (same normalised violation as the triangles; support 5 was [m] a wash) | local search, ms per region | integrality of a count | 1 day | joint pool on n13 and d2; kill at < +5 | not run |
| 4 | single-variable N+ (cell disjunction) | [m] +19 to +24 at node 10 | 2 DNN solves | hull of child and anti-child | none | needs > +40 | dead: it is a binary branch |
| 5 | one-line lift-and-project closure and coarser line disjunctions | [m] at most the minimum child: +42 to +59 at node 10 | ~47 child-size PSD blocks | convex hull of the children | 1-2 weeks | minimum child must reach +130 | dead: dominated by branching |
| 6 | count-split disjunctions | [m] +2.4 / +2.6 per split | 2 DNN solves | weak duality of one side constraint | done | both sides >= +23.4 | dead |
| 7 | certificate child bounds | [m] 0 fathomed at n10 (uncut) and n13 (cut); 14 of 2,116 at d2 | free | the certificate | none | - | flat |
| 8 | facility-region triangles and boxes (any T, incl. cut-decomposition regions of B + c(J - I)) | [m] 0 (max 0.006 uncut, 0.0001 after triangles) | QUBO per location triple | 0/1 region vector | - | a violation at a converged iterate | dead |
| 9 | parity, dispersion (k <= 7), flow-weighted sub-QAP, triple distance hull, stars | [m] 0: never violated | trivial | A.4 | - | - | dead |

No family reaches the +40 bar, and none shows a credible path to +130.

## (c) Runs on wil50 nodes

Cutoff 48,814.05; core 3; one solver at a time; `OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1`, `DNN_AMX=0` (this host has no
AMX; see the coordinator's note in `../README.md`); solver flags `-d 3 -L 2
-r 40`; 0.55-0.6 s per iteration at m = 47.

Driver `src/fdrv` (built by `src/build.sh`; links `latest/*.o` read-only plus
`src/dnncut_f.C`, a renamed copy of `latest/dnncut.C` with an export hook
generated by `src/make_dnncut_f.py`; cut algebra and certificate unchanged):
`-phase` runs the triangle protocol and saves the pool's Lagrangian term and
the ADMM state, `-cont` continues from it, `-pool` runs Polyak on an external
list of cuts on top of the frozen triangle term. A plumbing check
(`runs/plumb`) reproduces dnncuts' quiet round to every printed digit. All
values reported are certified; exact command lines are in `jobs.txt`
(executed by `runq.sh`).

1. **Node 13 triangle phase.** `fdrv -x 21:28,8:41,28:19 -u 48816 -d 3 -L 2
   -r 40 -grid 5:10 -k 1000 -I 400 -i 300 -rounds 8 -delta 50 -nocap -freeze
   1500 -phase tri wil50.dat` (`runs/n13/phase.log`): round 0 48,597.17, best
   round 7 48,700.75 (laptop 48,699.80). `fdrv ... -cont tri -iters 1500 -C
   -dumpY quiet -save quiet` (`runs/n13/quiet.log`): 48,745.63, 68.4 short
   (12-round laptop reference 48,755.42); certificate child bounds 0 of 2,209.
2. **Diagnostics** on the three node-13 iterates (`logs/*_n13_*.log`):
   `scripts/post.sh` (ydiag, box_sep, facreg_search, sep_dist (hull3),
   disp_check), `box2_search.py`, `star_check.py`, `lmp3.py [--loc]`,
   `gapdecomp.py`, `node_ls.py`.
3. **Node 10 cell disjunction** (`runs/n10c/{parent,child4945,elim4945}.log`):
   `dnntest -x 39:27,30:48,25:49[,49:45] [-E 49:45] -u 48816 -T 1 -d 3 -L 2 -r
   40 -i 1500 -H 0 [-C 2 -X]`; the lines-study reanalysis in
   `logs/lines_xweighted.log` and `logs/child_n10_parent.log`.
4. **LMP3 pools on node 13**, from the saved quiet state with Y and Z warm;
   cuts from `scripts/lmp3_sep.py` and `lmp3_sep_loc.py` on `quiet.Y.bin`, 60
   each, exact right-hand sides, all passing the permutation check.
   - Facility triples, standard step (`fdrv ... -pool quiet -wy -cuts
     lmp3_heavy60.cb -iters 200 -rounds 4 -i 300 -delta 5 -t 0.5 -freeze 800`,
     `runs/n13/lmp3pool.log`): round 0 48,746.66; the first step
     over-satisfies all 60 cuts and drops the bound by 38; rounds 2-4 recover
     to -22.6, -16.4, -12.9 with no cut violated; best multiplier 0; quiet
     48,747.82, the same as the plain continuation.
   - Facility triples, gentle step (delta 2, theta 0.25,
     `runs/n13/lmp3pool_gentle.log`): every cut over-satisfied after the first
     step; by round 4 only 11 multipliers active and all 60 cuts violated
     again; best 48,748.09, quiet 48,747.96: the plain continuation again.
   - Location triples, gentle step (`runs/n13/lmp3locpool.log`): rounds 1-2
     +1.29, +1.54, then two rounds overshoot; best 48,748.21, quiet 48,748.09:
     +0.3 to +0.5 over the plain continuation at equal iterations.
5. **Count splits on node 13** (`scripts/split_jobs.sh`, `pick_split.py`,
   `count_term.py`; `fdrv -cont quiet -wy -iters 400 -extra term:const`):
   results in A.5 item 6. A first control that restarted Y from scratch (Z
   only) certified 48,713.37, 32 below its starting state (`runs/n13/zonly/`);
   all split runs were rerun with Y warm.
6. **Child d2** (39:27,30:48,25:49,26:5), same phase and quiet commands
   (`runs/d2/phase.log`, `runs/d2/quiet.log`): round 0 48,657.14; best round
   5 48,726.31; round 6 overshot by 264; quiet from round 5 48,784.25 (laptop
   12-round 48,789.87); certificate child bounds 14 of 2,116.

## (d) Recommendation

1. **Close the facility side as a source of new inequality families**
   (facility clusters, facility regions including the cut-decomposition
   regions of B + c(J - I), boxes over facility sets, QSAP clique/cut
   inequalities, stars, parity, dispersion, flow-weighted sub-QAP cuts):
   wil50's flows are random, they already enter every regional family as
   bisections of the flow graph, and at converged iterates nothing on the
   facility side is violated except triple consistency, which is repaired for
   free.
2. **Do not build lift-and-project or disjunctive cuts for this relaxation.**
   Each costs DNN solves per term and is capped by its worst term, the child
   of the most-favoured cell (+42 to +59 for node 10's best lines, +19 to +24
   for one N+ level, +2.4 to +2.9 for count splits); branching on the same
   line pays for the same solves and discards the 33 of 47 children that
   fathom; the DNN is already the first N+ level; and the certificate gives
   almost no child bounds.
3. **Retire the triple-consistency family (step 8) on the heaviest flow
   triples too.** Heavy triples raise the violations tenfold but the
   Lagrangian gain stays at 0; the location analogue gives +0.3 to +0.5.
4. **Where a level could still come from:** moments that couple the 13
   threshold memberships of each facility across thresholds (the row-column
   coupling the grid agent measured), not the facility side or disjunctions.
   Until such a family exists, the measured levers stay where the roadmap puts
   them: inheriting the cut pool, warm starts, branching.
5. **By-products for the search.** Restarting a node from Z alone lost 32
   certified units at 400 iterations, while Y and Z gained +1.2: carry Y
   wherever a converged state is resumed (checkpoints, state files,
   children). The bottleneck of a line is its largest-x cell, and the lightest
   facility (45) gives the best line at all three lines-study parents; rules 8
   and 9 already use this.

## References

arxiv.org and kups.ub.uni-koeln.de did not resolve from this session, so no
full texts were read; [s] = from search-result abstracts, [m] = from memory.

- Lovasz, Schrijver (1991). Cones of matrices and set-functions and 0-1 optimization. SIAM J. Optim. 1:166-190 [m].
- Sherali, Adams (1990). A hierarchy of relaxations between the continuous and convex hull representations for zero-one programming problems. SIAM J. Discrete Math. 3:411-430 [m].
- Balas, Ceria, Cornuejols (1993). A lift-and-project cutting plane algorithm for mixed 0-1 programs. Math. Prog. 58:295-324 [m].
- Burer, Vandenbussche (2006). Solving lift-and-project relaxations of binary integer programs. SIAM J. Optim. 16(3):726-750 [s].
- Povh, Rendl (2009). Copositive and semidefinite relaxations of the quadratic assignment problem. Discrete Optim. 6(3):231-241 [s].
- de Klerk, Sotirov (2010). Exploiting group symmetry in semidefinite programming relaxations of the QAP. Math. Prog. 122:225-246 [s].
- Junger, Kaibel (2001). Box-inequalities for quadratic assignment polytopes. Math. Prog. 91:175-197; The QAP-polytope and the star-transformation. Discrete Appl. Math. 111:283-306 [s].
- Saito, Fujie, Matsui, Matuura (2009). A study of the quadratic semi-assignment polytope. Discrete Optim. 6:37-50 [s].
- Padberg (1989). The Boolean quadric polytope. Math. Prog. 45:139-172 [m].
- Deza, Laurent (1997). Geometry of Cuts and Metrics [m].
- Saxena, Bonami, Lee (IPCO 2008; Math. Prog. 124, 2010). Disjunctive cuts for nonconvex MIQCPs [s].
- Lodi, Tanneau, Vielma (2023). Disjunctive cuts in mixed-integer conic optimization. Math. Prog. 199:671-719 [s].
- Gaar, Rendl (2020). A computational study of exact subgraph based SDP bounds for Max-Cut, stable set and coloring. Math. Prog. 183:283-308 [s].
- Mittelmann, Peng (2010). Estimating bounds for QAPs associated with Hamming and Manhattan distance matrices based on SDP. SIAM J. Optim. 20:3408-3426 [s].
