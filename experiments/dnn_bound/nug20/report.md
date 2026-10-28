# Would a DNN bound fathom the depth-2/3 nodes of nug20 that QPB cannot?

Short answer: **yes, almost all of them (96/98 at depth 2, 200/200 sampled at depth 3), with certified
bounds - but it costs 8.3 CPU s per depth-2 node and 3.8 s per depth-3 node, while the entire QPB run of
nug20 is 160-200 s. At nug20 scale the strong bound loses by ~5x at depth 2 and ~30x at depth 3.** The
bound quality is not the problem; the price per node is.

All code, state files and per-node results are in
`/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn/`
(`dnn.py` bound, `qaputil.py` parsing/reduction, `batch.py` node runs, `verify_exact.py` exact checks,
`results_depth{1,2,3}.jsonl`, `converged_results_*.jsonl`, `verify.log`, `*.log`). Nothing under
/home/user/qap was modified. Every Python run used one core (OMP/OPENBLAS/MKL_NUM_THREADS=1).

## 1. What was computed

**Bound.** The DNN relaxation of Zhao-Karisch-Rendl-Wolkowicz in the facially reduced form of
Oliveira-Wolkowicz-Xu (OWX 2018): lifted matrix Y of size m^2+1 on [1; vec(X)], Y = V R V' with V
orthonormal ((m-1)^2+1 columns, spanning the row/column-sum-feasible subspace), R psd, gangster
constraints Y_J = 0, Y_00 = 1, and 0 <= Y <= 1. Solved by the OWX ADMM (R-step: eigenprojection of
V'(Y+Z/beta)V; Y-step: clip to the polyhedral set; Z-step with gamma = 1.618; beta = m/3 after scaling
L by max|L|).

**Certified lower bound.** The ADMM primal value is *not* used as a bound. From the current multiplier
Z the Lagrangian dual value is evaluated: Zt = Z - V P_psd(V'ZV) V' (so V' Zt V is negative
semidefinite), then lb = (L+Zt)_00 + sum over non-gangster entries of min(0, (L+Zt)_ij). For any such Zt
this is a valid lower bound on the DNN value, hence on the subproblem optimum (permutation matrices are
feasible for the relaxation; checked numerically for random permutations incl. the linear term). It is
valid at every iteration, so the ADMM is stopped as soon as lb > 2568.05 ("fathomed"), or at 3000
iterations / 150 s. The subproblem of a node fixing facilities f to p(f) is
tr(A_FF Y B_LL Y') + <C', Y> + const with C'[i][l] = 2 sum_f A[i][f] B[l][p(f)], verified against the
full objective on random completions for every node processed.

**Nodes.** From the solver (`bnbdrv -p PARAM/default.par -u 2570 -d D -o state`): with -d D the
solver bounds nodes through depth D and writes their surviving children, whose `bound` field is the
parent's effective bound max(QPB at the node, inherited bound, predicted bound). Grouping children by
parent (history[d] = facility branched at depth d) gives every bounded node's QPB. Counts reported by
the solver with incumbent 2570:

| depth | nodes opened | QPB fathomed | all children eliminated by reduced costs | QPB could not fathom (branched) |
|---|---|---|---|---|
| 1 | 6 | 0 | 0 | 6 |
| 2 | 98 | 1 | 1 | 96 |
| 3 | 1474 | 62 | 2 | 1410 |

DNN was run on all 98 depth-2 nodes (m = 18), on a uniform random sample (seed 0) of 200 of the 1410
branched depth-3 nodes (m = 17), and on the 6 depth-1 nodes (m = 19).

## 2. Results

### Root bounds (sanity)

| instance | optimum | QPB root (solver) | DNN certified lb | DNN primal | ADMM iters / CPU s | literature DNN (from memory, unverified) |
|---|---|---|---|---|---|---|
| nug12 | 578 | 478.47 | 567.99 | 567.99 | 900 / 2.0 | 568 |
| nug20 | 2570 | 2243.94 | 2506.33 | 2506.35 | 1600 / 40 | ~2506 |

### Depth 2 of nug20 (all 98 opened nodes, m = 18)

* DNN certified bound > 2568.05: **96 of 98** (98.0%). Of the 96 nodes QPB could not fathom, DNN
  fathoms 94 (97.9%).
* The two misses: node {1->6, 19->12}: DNN value ~2561.6 (certified lb 2560.59, primal 2561.63);
  its exact subproblem optimum is **2570** - it contains the optimal solution, so no bound can fathom
  it. Node {1->6, 19->13}: DNN ~2567.1 (lb 2566.73, primal 2567.14), exact optimum 2588: a genuine
  DNN miss by about 1 unit.
* QPB at these nodes (percentiles 0/10/25/50/75/90/100): 2346.5 2367.9 2377.0 2402.6 2437.3 2474.5 2563.5
* Certified DNN lb at early stop: 2560.6 2568.3 2569.0 2571.3 2579.1 2590.1 2646.5
* lb - QPB: 79.7 115.6 139.4 169.8 193.1 200.9 217.1 (median +170; QPB's median gap of 167 to the
  optimum is closed entirely).
* 2570 - lb: -76.5 -20.1 -9.1 -1.3 1.0 1.7 9.4 (negative = bound above the optimum; these are
  early-stopped values, the true DNN values are much higher: five nodes run to convergence gave DNN
  values 2572.5, 2572.2, ~2640, ~2732, ~2592 where the early-stopped lbs were 2568.05, 2568.10,
  2579.6, 2588.7, 2570.4).
* CPU seconds per node (single core, Python/numpy, 25 ms per ADMM iteration at lifted size 325 with a
  290-column V): percentiles 2.2 3.1 4.4 6.7 9.2 15.0 43.2, **mean 8.3 s, total 818 s** for 98 nodes.
  Iterations to certify: median 500 (150 to 3000); 12 nodes needed > 1000.
* Accuracy needed: because most DNN values are 10-150 above 2570, a crude dual point suffices; the
  nodes whose DNN value is within ~4 of 2570 (about 10% of them) needed 1300-3000 iterations, and the
  smallest certified margin over 2568.05 was 0.001 (node {1->7,19->11}; converged DNN value 2572.5,
  exact optimum 2594).

### Depth 3 of nug20 (200 random branched nodes of 1410, m = 17)

* DNN certified bound > 2568.05: **200 of 200** (95% Clopper-Pearson lower confidence bound of
  98.2% for the fathom rate over the 1410).
* QPB percentiles: 2397.2 2441.3 2463.7 2486.9 2513.5 2546.0 2566.7
* lb (early stop): 2568.1 2569.3 2571.4 2576.5 2583.6 2595.1 2614.2
* lb - QPB: 14.7 45.9 71.4 92.3 113.7 131.9 170.9
* 2570 - lb: -44.2 -25.1 -13.6 -6.5 -1.4 0.7 1.9 (converged DNN values of five nodes: 2580.6,
  2635.6, 2660.3, 2653.5, 2718.4)
* CPU s per node: 1.6 2.3 2.8 3.4 4.4 5.4 13.8, **mean 3.8 s, total 769 s**; iterations median 300
  (150 to 1250), ~20 ms per iteration (lifted size 290, 257-column V). Smallest certified margin 0.013.

### Depth 1 (6 nodes, m = 19), for context

DNN fathoms 1 of 6 ({1->0}: lb 2568.39). The others converge to DNN values 2558.6, 2556.5, 2566.0,
2540.3, 2543.9 (QPB 2302-2355), i.e. the strong bound only "switches on" once two facilities are
fixed. 38-80 CPU s each.

### Validity checks (all passed)

* nug12 subproblems with the linear term folded in: brute force at m = 6, 7, 8, 8, 8 - DNN lb equals the
  optimum (702, 648, 664, 666, 654). m = 10 reduced instance written as a 3-matrix .dat and solved by
  bnbdrv: optimum 614, DNN lb 613.95, QPB root 554.65.
* nug20 nodes solved exactly by bnbdrv on the reduced 3-matrix instance (C' and const folded in):

| node | m | DNN certified lb | exact subproblem optimum | QPB |
|---|---|---|---|---|
| d2 {1->6,19->12} | 18 | 2560.59 | 2570 | 2346.5 |
| d2 {1->6,19->13} | 18 | 2566.73 | 2588 | 2355.8 |
| d2 {1->7,19->11} | 18 | 2568.05 | 2594 | 2353.9 |
| d2 {1->6,19->8} | 18 | 2568.10 | 2594 | 2357.2 |
| d2 {0->19,1->0} | 18 | 2579.56 | 2646 | 2407.5 |
| d3 {1->1,12->19,16->15} | 17 | 2568.06 | 2584 | 2397.2 |
| d3 {1->6,15->4,19->10} | 17 | 2568.07 | 2644 | 2436.1 |
| d3 {1->6,5->16,19->1} | 17 | 2568.24 | 2638 | 2476.8 |
| d3 {1->0,5->19,11->4} | 17 | 2568.35 | 2716 | 2543.2 |
| d3 {1->5,2->19,8->4} | 17 | 2569.65 | 2602 | 2430.2 |

(QPB-only solves of these m = 18 subtrees took 73k-238k nodes, 27-61 s each; m = 17: 24k-74k nodes, 8-16 s.)

## 3. Verdict: would a strong bound at depths <= 3 pay?

Bound strength is not the issue - the DNN bound fathoms ~98% of the depth-2 and essentially all of
the depth-3 nodes that QPB leaves open. The issue is cost per node versus the cost of the QPB subtree
it would remove:

* **nug20** (QPB serial run: ~600k nodes, 160-200 s). The average QPB subtree below a depth-2 node is
  ~180 s / 98 = 1.8 s; below a depth-3 node ~0.12 s. DNN costs 8.3 s and 3.8 s respectively (with the
  cheapest possible stopping rule: stop the moment the certified bound clears the incumbent). DNN at
  depth 2 alone: 818 s + ~60 s of QPB for the two surviving subtrees (102k and 73k nodes) ~ 900 s
  vs ~180 s: **about 5x slower**. DNN at depth 3: 1474 x 3.8 = 5600 s, 30x slower. Break-even at
  depth 2 would need DNN at ~1.5 s per node, i.e. ~5x faster than measured; one 290x290 eigh is ~10 ms
  and the median node needs 500 of them, so even a tight C implementation (maybe 1.5-2x faster than
  numpy, since LAPACK dsyevd + dgemm dominate) does not get there.
* **nug24** (~27M nodes, 5000-6000 CPU s). Per-node DNN cost grows like (m^2)^3: m = 22 lifted size 485
  with a 442-column V, ~90 ms per iteration, 300-600 iterations -> ~30-55 s per depth-2 node (more if
  nug24's depth-2 nodes are harder to certify). A depth-2 QPB subtree costs roughly 5500 s / (number of
  depth-2 nodes, probably 100-300) ~ 20-50 s. So nug24 is near break-even at depth 2 *if* the fathom
  rate carries over, which is unknown (at nug20 the rate went from 1/6 at depth 1 to 96/98 at
  depth 2 - the transition is sharp and may sit deeper for a larger instance). Depth 3 loses again
  (subtree ~2 s vs DNN ~15-30 s).
* **nug30** is where the arithmetic would flip (depth-2 subtrees cost hours-days of QPB; DNN at m = 28,
  lifted size 785, is ~0.4 s per iteration, 200-800 s per node), but whether the DNN bound still
  fathoms depth-2/3 nodes at m = 27-28 is untested here; the depth-1 nug20 result (1 of 6 at m = 19)
  suggests the relaxation needs the subproblem to be small relative to n before it clears the gap.
  That would be a 1-2 hour follow-up experiment (20 nug30 depth-2 nodes at ~5 min each).

So: at nug20/nug24 scale a DNN bound at depths <= 3 does not pay; the direction is only worth
pursuing for nug30-class instances, and only after a nug30 depth-2 fathom-rate test. Note also that
the early-stop fathom test needs the incumbent (2570 here, as the solver had with -u 2570).

## 4. Not done / uncertainties

* RLT1 (Hahn-Grant / Adams-Johnson dual ascent) was not implemented; DNN turned out practical.
* Depth 3 is a 200-node uniform sample of the 1410 branched nodes, not all of them. Depth-3 nodes QPB
  fathomed (62) or child-eliminated (2) were excluded from the sample on purpose.
* Bound validity rests on (i) the Lagrangian-dual derivation above, (ii) numerical checks that
  permutations are feasible for the relaxation, (iii) 16 exact comparisons (lb <= optimum, all passed).
  Floating-point: eigenvalue projection errors are ~1e-9 on these scales; certified margins below
  ~1e-6 (none here; the smallest was 0.001) should not be trusted. Early-stopped lbs understate the
  DNN value; the "obj" (ADMM primal) values in the jsonl files are approximate and not bounds.
* Literature DNN values (nug12 ~568, nug20 ~2506) are from memory and unverified; they agree with the
  computed values but were not looked up.
* The per-node QPB values are the solver's effective bound max(QPB, inherited, predicted); where QPB
  fathomed via child elimination the node's own QPB is not in the state file (2 depth-2 nodes show
  `qpb: null`).
* Timing is Python 3.13 + numpy 2.5 (scipy-openblas, 1 thread); nug24/nug30 costs are extrapolations.
