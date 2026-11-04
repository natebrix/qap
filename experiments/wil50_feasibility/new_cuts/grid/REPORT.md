# The location side: what wil50's 5 x 10 grid offers beyond the regional triangles

Agent 2 (grid), Nov 4, 05:26-10:30 UTC, core 2. Saved by the coordinator from
the agent's hand-back message (the agent could not write report files); the
headline split bounds and the validity check were re-read from the logs
(`runs/*_tri5*.log`, `runs/validate_objpart*.log`). `NOTES.md` is the
agent's working log; scripts, instances and logs are in this directory and
`runs/`, and paths below are relative to it.

Tags: [m] measured on wil50 nodes (log named), [e] estimate. Test nodes (spec
`i:j` = cell : facility): node 13 = `21:28,8:41,28:19` (depth 3, m = 47), node
10 = `39:27,30:48,25:49` (depth 3), d2 = `39:27,30:48,25:49,26:5` (depth 4, m
= 46). Cutoff 48,814.05.

## Summary

1. **No new inequality family from the grid's geometry is worth carrying.**
   Tested: parity (checkerboard) triangles; clique, star and flow-weighted
   dispersion (packing); sub-QAP cuts on facility cliques; triangles on
   rectangles, quadrant checkerboards, interior rows and columns; a free
   search over all cell regions; l1, hypermetric and negative-type conditions
   on the expected distance matrix. Each is either not violated at all at the
   uncut and triangle iterates, or dominated by the prefix triangles [m]. The
   DNN does not over-cluster facilities anywhere.
2. **The row/column decomposition is useless as a bound, and the triangle
   phase leaves nothing in it.** The decoupled bound OPT_row + OPT_col is
   below the uncut DNN at node 13 (at most 48,584 against 48,627) [m]. After
   the triangle phase every one of the 13 threshold terms is at or above its
   own optimum, in the objective's direction, on nodes 13 and 10, and so are
   the row and column problems as wholes [m]. This explains why every
   per-region family failed (5-point, local-7, and on the same argument the
   per-region M2 block): they are inequalities of hulls the iterate already
   respects where the objective looks.
3. **The remaining gap is the two-dimensional coupling of the middle row and
   column thresholds, and it is global.** Grouping cells by the five middle
   thresholds (row<2, row<3, col<4, col<5, col<6) gives 12 "super-cells";
   after the triangles the iterate's value of those five terms is below that
   12-super-cell problem's optimum by 117 at node 13 and 143 at node 10 [m].
   Any single threshold, or any row x column pair, is within 0-25 of its own
   optimum. The deficit grows with the number of thresholds coupled (213 and
   260 for nine), and no few-facility family sees it.
4. **The coupling can be harvested, but only by solving a coarse QAP
   exactly.** The cut <L_K, Y> >= OPT_K(node) is linear and valid, with K the
   five middle thresholds and OPT_K the 12-super-cell problem's optimum at the
   node. At multiplier 0.5 it adds **+95 to the uncut bound at node 13, and
   +38.3 / +35.6 / +32.1 on top of the triangles on nodes 13 / 10 / d2**
   (median +35.6) [m], **conditional on the coarse optimum**, which here comes
   from simulated annealing. Proving it is the problem: the coarse problem's
   own bound is 23,950 plain and 23,989 with triangles against an optimum of
   24,250, a 1.0-1.2% gap, larger than the node's own; one level of coarse
   branching raises the certified coarse bound by only about 40 at its worst
   child [m]; breaking even needs about 5 coarse levels.

**Recommendation: stop looking for location-side inequalities for wil50.**
The prefix triangles have used up the grid's structure. What remains is a
coupling that only branching resolves; an exact coarse subproblem is
branching in disguise. No family here meets the bar (+40 on top of converged
triangles at no more than about 2x the phase's CPU), and there is no credible
path to +130 from the location side.

## 1. The structure, and what the triangles already take

- **Threshold identity.** d(a,b) = |r(a)-r(b)| + |c(a)-c(b)| = sum over 13
  thresholds t (row<1..4, col<1..9) of |1[a in S_t] - 1[b in S_t]|, checked
  entry by entry (`gridlib.py`). The objective is sum_t cut_F(z^t): 13
  capacitated bisections that one bijection between facilities and cells must
  realise.
- **Three coupling layers on top of the bisections:** nesting within the rows
  (4 thresholds), nesting within the columns (9), and the row x column
  coupling (each row meets each column in exactly one cell).

Measured split at node 13 [m] (optima are simulated-annealing values from
`bin/qapsa`, upper bounds on the true optima; 2-5 seeds agree):

| quantity | value |
|---|---|
| sum of the 13 bisection optima (fully decoupled) | 48,048 |
| row-only + column-only optima (nesting kept) | 15,902 + 32,682 = 48,584 |
| best known node solution | 49,022 (row 16,248 + column 32,774) |
| uncut DNN, 1,500 it. (`runs/n13_uncut.log`) | 48,627.3 (row part 16,000, column part 32,632) |
| triangle phase, 5 rounds + 600 quiet (`runs/n13_tri5.log`) | 48,713.9 (row 16,027, column 32,688) |

The row x column coupling is worth at least 49,022 - 48,584 = 438 at node 13;
the uncut DNN captures at most 43 of it. At the root, SA gives row-only +
column-only <= 48,184 against a DNN root of 48,135 and an optimum of 48,816.

## 2. Local families: measured and closed

Measured on node 13's uncut and triangle-phase iterates and on the
triangle-phase iterates of node 10 and d2, in the (p, W) units of the
production logs (reference: prefix triangles 0.073-0.080 in the middle
regions uncut, 0.03-0.06 after the triangles; `analyze_Y.py`).

- **Parity.** Three facilities cannot be pairwise at odd distance on a
  bipartite grid (the regional triangle on the checkerboard). Checkerboard,
  row-parity and 2x2-block-checkerboard regions: zero violated triples at all
  four iterates; column parity at most 2 triples at 0.02.
- **Clique dispersion.** sum over pairs in T of d >= Dmin(|T|), Dmin exact
  from row and column margins plus a max-flow check (`dispersion.py --check`
  matches enumeration; k = 2..10: 1, 4, 8, 16, 25, 38, 54, 72, 96, the
  "n-town" numbers). The most clustered k-sets sit above Dmin by 0.8-1.2 (k =
  3) up to 13-20 (k = 8-9): never violated. Star versions are implied by the
  DNN outright (the column block of Y is doubly stochastic given facility i at
  cell a: the Gilmore-Lawler argument the certificate already uses).
- **Flow-weighted dispersion / sub-QAP cuts** (sum over T of F_ij d_ij >=
  Q(T), Q the sub-QAP optimum by SA, optimistic): slack by 140-1,700 at all 50
  greedy high-flow cliques of 4-12 facilities (`runs/n13_uncut.subqap.log`).
- **Other regions.** A 1-flip search over all cell subsets for the most
  violated triangle cannot improve the prefix regions (col<5 stays at 0.0801;
  col<6 climbs to 0.079 by turning into col<5); from the parity regions it
  climbs to 0.02-0.06 (`runs/n13_uncut.regionsearch.log`). Corner rectangles
  and quadrant checkerboards (row<t xor col<s): max 0.067 uncut, 0.045 after
  the triangles, sums of violations 10-30x smaller (`runs/n13_*.xor.log`).
  Interior single rows and columns: at most 0.027. This repeats, for regions,
  Astra's node-1 finding for 770 rectangles and the mixed-region triangles.
- **l1 / hypermetric structure of E[d].** E[d] is of negative type (max
  eigenvalue of P E P is 0 to rounding); the pentagonal hypermetric is not
  violated (-0.52 / -0.66 over 20,000 random 5-sets); the triangle inequality
  on E[d] is violated by up to 0.96 uncut and 0.39 after the triangles, but it
  is the sum of 13 regional triangles with equal multipliers, which the pool
  holds with free multipliers (`runs/n13_*.metric.log`).

## 3. The coupled problems: where the gap is

For a threshold set K, compare the iterate's K-part (those cut terms plus the
node's linear terms) with the optimum of the K-part alone, which couples the
thresholds through the joint cell capacities (K = {row<3, col<5} is a
4-quadrant partition; mid5 is 12 super-cells of 8,2,2,8 / 4,1,1,4 / 8,2,2,8
cells). Negative = the iterate lies outside the coupled hull in the
objective's direction (`coupled_split.py`, `root_vs_node.py`,
`threshold_split.py`).

| K | node 13 uncut | node 13 after tri5 | node 10 after tri5 |
|---|---|---|---|
| single thresholds | -17..+118 | **+0.3..+113** | **+2.9..+125** |
| row problem / column problem | +98 / -50 | +125 / **+6** | +137 / +77 |
| mid row x column pairs | -36..-62 | -16..-25 | -14..+16 |
| mid5 (12 super-cells) | -214 | **-117** | **-143** |
| 7 middle thresholds | -282 | -181 | -190 |
| 9 thresholds | -303 | -213 | -260 |
| mid5 against the root-level optimum | -76 | +21 | -25 |

d2's triangle iterate is inflated (its <L,Y> is 123 above its certified
value: the quiet round's multipliers over-satisfy the cuts); its readings
agree in sign (9 thresholds: -141) but are not used for magnitudes. The
triangles push each threshold, and with it the row and column problems, into
its own hull, leaving nothing for a family acting on one region or one
direction. The deficit is not local in thresholds (pairs 16-25, five 117-143,
nine 213-260) nor in facilities (section 2), and it carries most of the
remaining gap: at node 13 the best known solution's mid5 part is 24,294
against the iterate's 24,133, its other eight thresholds 24,728 against
24,582.

## 4. Objective splitting: harvesting the coupling, and its price

- **The cut.** Every permutation satisfies <L_K, Y> >= OPT_K(node). With
  multiplier alpha the Lagrangian term is -alpha L_K, i.e. the node of the
  instance A' = sD - s alpha D_K (`make_inst.py`, s = 2), and the bound is
  DNN(A')/s + alpha OPT_K. It runs with no solver change, triangles included
  (`dnncuts -grid 5:10` on `inst/wil50_mid5_a1of2.dat`).
- **What it can gain.** The minimum is superadditive, so the gain is at most
  alpha (OPT_K - relaxation(L_K)): only the coarse problem's own relaxation
  gap can be harvested.
- **Validity.** `validate_objpart.py` against enumeration on 20 small grid
  instances (2x3, 2x4, 3x3; random flows, random K, 0-2 fixes; alpha
  0.25/0.5/0.75), 8 with a DNN gap: **80 certified values, 0 above the node
  optimum** (`runs/validate_objpart*.log`). With K = all thresholds the split
  bound tends to the optimum as alpha goes to 1.

Results at alpha = 0.5 [m]; OPT_K by SA: 24,250 for node 13 (2 x 30 runs
agree), 24,230 for node 10 and d2 (5 seeds).

| node, protocol | split bound | same protocol without the cut | gain | break-even OPT_K |
|---|---|---|---|---|
| 13, uncut 1,500 it. (`runs/n13_mid5a1of2.log`) | 36,597.41 + 12,125 = 48,722.4 | 48,627.3 | **+95.1** | |
| 13, tri5 (`runs/n13_tri5_mid5a1of2.log`) | 36,627.14 + 12,125 = 48,752.1 | 48,713.9 | **+38.3** | 24,173 |
| 10, tri5 (`runs/n10_tri5_mid5a1of2.log`) | 36,632.26 + 12,115 = 48,747.3 | 48,711.6 | **+35.6** | 24,159 |
| d2, tri5 (`runs/d2_tri5_mid5a1of2.log`) | 36,666.73 + 12,115 = 48,781.7 | 48,749.6 | **+32.1** | 24,166 |
| 13, tri5, root-level coarse optimum 24,112 | 48,683.1 | 48,713.9 | -30.7 | |
| 13, tri5, best depth-2 ancestor's coarse optimum 24,198 | 48,726.1 | 48,713.9 | +12.3 | |

At the 400-iteration give-up point the cut adds +71 to +105; the triangles
gain less on top of it (+51 against +117 at node 13), so the two overlap.

**The price is the certificate.** Node 13's coarse problem
(`inst/wil50_mid5only.dat`) has a plain converged bound of 23,949.9
(`runs/coarseconv_*.log`) and 23,988.9 after the 5-round triangle phase
(`runs/n13_tri5_mid5only.log`), against an optimum of 24,250. One coarse
level (facility 40 into each of the 12 super-cells; `coarse_bb.py`,
`runs/n13_coarse_level1_f40.log`): children 23,968-24,044 at 400 iterations,
worst child converged 23,989.8, **+39.9 per level**; the break-even (24,173)
is 223 above the plain coarse bound, about 5 coarse levels. Key sharing and
node cost help but not enough: coarse branching is 12-way (cells of a
super-cell are interchangeable for the K-part), a coarse problem is shared by
about 72 depth-3 nodes per key [e], and a symmetry-reduced coarse DNN (de
Klerk-Sotirov) would have a ~507-dimensional main block instead of 2,117,
about 70x cheaper per iteration [e]. Net: 10^4-10^5 coarse nodes per key,
about 5-500 DNN-equivalents per depth-3 node, for +32 to +38, and a node whose
coarse key equals its parent's gains nothing new (d2's coarse optimum equals
node 10's). Killed by the pre-declared criterion (fewer than 60 per coarse
level).

## 5. Ranked candidates

| # | family / bound | gain at depths 3-5 | separation / cost | validity | effort | cheapest test -> kill |
|---|---|---|---|---|---|---|
| 1 | objective splitting with the node's own coarse (mid5) optimum | +95 uncut; **+32 / +36 / +38 on top of tri5** [m, conditional on OPT_K] | an exact coarse solve per coarse key: 12-way B&B, 1.0-1.2% gap, ~5 levels | proven given OPT_K; enumeration check 80/0 | 2-4 weeks (symmetry-reduced coarse DNN, coarse B&B, key cache) | one coarse level at node 13. **Killed: +40 < 60 per level** [m] |
| 2 | the same with 7-9 thresholds | first-order slope 181-260 on top of tri [m]; +50-70 [e] | coarse problem nearly as fine as the node (18-30 super-cells) | as 1 | harder than 1 | not worth running: certifying it is about solving the node |
| 3 | inherited ancestor (depth 1-2) coarse optimum | -31 to +12 on top of tri [m, SA proxies] | one coarse solve per depth-2 key | as 1 | as 1 | killed by the numbers |
| 4 | row/column decomposition (exact 1-D problems, Lagrangian along L_row, L_col) | 0 [m] | exact 5- and 10-part partitions | proven | - | closed |
| 5 | parity triangles / odd cycles | 0 [m]: zero violated at 4 iterates | trivial | regional triangle | 1 h | closed |
| 6 | clique / star / flow-weighted dispersion, sub-QAP cuts | 0 [m]: slack 0.8-1,700 | exact Dmin, small QAPs | proven; star DNN-implied | - | closed |
| 7 | free-region, rectangle, quadrant-XOR, interior row/column triangles | at most the prefix triangles [m] | cheap | regional triangle | 1 day | closed (the bound memo's item 5, +5-15 expected, looks like 0-5) |
| 8 | hypermetric / negative type on E[d] | 0 [m] | eigenvalue / 5-sets | l1 metric | - | closed |
| - | per-region M2 block (bound memo item 3), as a corollary | [e] about 0 | - | - | - | first-order evidence only, consistent with the local-7 and 5-point results |

## 6. What the grid still offers

The structure helps the search rather than the bound; the gap after the
triangles is in how facilities spread jointly over the middle rows and
columns. Two uses [e, untested]: a two-resolution search (a coarse 12-way
tree alongside the fine one, each fine node getting the split cut from the
coarse node with the same coarse fixings; item 1 as a search design, though
the coarse relaxation being weaker than the fine one, 1.0-1.2% against 0.6%,
argues against it); and branching guidance (the facilities whose coarse
marginals are most spread over the middle super-cells carry the deficit; at
node 13 these are 40, 37, 38, 44 and 23; a rule preferring them could be
replayed against rule 9 offline). Bound work belongs to the families already
found: pool inheritance (Astra's +54 median), the dual method and certificate
tuning, and per-node speed.

## 7. Commands and logs

From the repository root with `G=experiments/wil50_feasibility/new_cuts/grid`,
`export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=1 DNN_CHOLQR_PASSES=1`,
everything under `taskset -c 2`. This host has no AMX (0.31-0.62 s per
iteration at m = 46-47), so the triangle reference is the 5-round "tri5"
protocol, not the laptop's 12 rounds: tri5 48,714 / 48,712 / 48,750 on nodes
13 / 10 / d2 against the laptop's 48,755 / 48,744 / 48,790.

- Uncut iterate: `latest/dnntest -x 21:28,8:41,28:19 -u 48816 -T 1 -d 3 -L 2 -r 40 -i 1500 -H 0 -v -dumpY $G/runs/n13_uncut qaplib/wil50.dat` (48,627.29).
- tri5: `experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 5 -delta 50 -nocap -freeze 0:600 -x SPEC -dumpY $G/runs/ID_tri5 qaplib/wil50.dat` (`runs/{n13,n10,d2}_tri5.log`).
- Splitting (after `python3 $G/make_inst.py`): uncut at alpha 0.5, `latest/dnntest -x SPEC -c 1e12 -T 1 -d 3 -L 2 -r 40 -i 1500 -H 0 -v $G/inst/wil50_mid5_a1of2.dat`; on top of the triangles, the tri5 command on `$G/inst/wil50_mid5_a1of2.dat` with `-c 1e12 -delta 100`; the coarse bound, the tri5 command on `$G/inst/wil50_mid5only.dat` with `-c 1e12 -delta 25` (`runs/coarse_conv.sh`); one coarse level, `python3 $G/coarse_bb.py 21:28,8:41,28:19 40 --iters 400`.
- Analyses (seconds each; arguments SPEC Y.bin): `analyze_Y.py`, `threshold_split.py`, `coupled_split.py --pairs`, `root_vs_node.py`, `region_search.py`, `xor_regions.py`, `metric_check.py`, `subqap_viol.py`; also `decomp_ub.py SPEC` and `check_coupled.py SPEC K RUNS SWEEPS`.
- Validity: `python3 $G/dispersion.py --check`; `python3 $G/validate_objpart.py -N 12 --seed 3` and `-N 8 --seed 11 --gap-only`.

## References

- H. Mittelmann, J. Peng, "Estimating bounds for quadratic assignment problems associated with Hamming and Manhattan distance matrices based on semidefinite programming", SIAM J. Optim. 20 (2010) 3408-3426. alpha E - XBX' PSD for small SDPs up to n = 200; their wil50 values (46,526 / 46,467, as quoted in Astra's memo) are far below the DNN, and the constraint is implied by the full moment PSD.
- E. D. Demaine, S. P. Fekete, G. Rote, N. Schweer, D. Schymura, M. Zelke, "Integer point sets minimizing average pairwise L1 distance: what is the optimal shape of a town?", CCCG 2009 (Comput. Geom. 2011). Source for the Dmin values.
- E. de Klerk, R. Sotirov, "Exploiting group symmetry in semidefinite programming relaxations of the quadratic assignment problem", Math. Prog. 122 (2010).
- E. de Klerk, D. Pasechnik, R. Sotirov, C. Dobre, "On semidefinite programming relaxations of maximum k-section", Math. Prog. B 136 (2012) 253-278.
- A. Caprara, A. N. Letchford, J.-J. Salazar-Gonzalez, "Decorous lower bounds for minimum linear arrangement", INFORMS J. Comput. 23(1) (2011) 26-40.
- A. R. S. Amaral, A. N. Letchford, "A polyhedral approach to the single row facility layout problem", Math. Prog. 141 (2013) 453-477.
- M. F. Anjos, A. Fischer, P. Hungerlander, "Improved exact approaches for row layout problems with departments of equal length", EJOR 270(2) (2018) 514-529.
- M. Junger, V. Kaibel, "Box-inequalities for quadratic assignment polytopes", Math. Prog. 91 (2001).
