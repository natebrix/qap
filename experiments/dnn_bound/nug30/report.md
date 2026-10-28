# nug30: would a certified DNN bound at depth 1/2/3 pay off against QPB branch-and-bound?

Short answer: **yes, by one to two orders of magnitude.** At nug30 the certified DNN bound fathoms
**34 of 36 random depth-2 nodes (94%)** and **19 of 20 depth-3 nodes** at a measured cost of **~130 CPU s
per fathomed depth-2 node** (Python/numpy, one core), while the solver's own tree-size estimator puts
the QPB subtree below an *average* depth-2 node at **~22,000-46,000 CPU s** (6-13 CPU-hours) and the
whole QPB-only nug30 run at **30-120 CPU-days** on this machine. The DNN misses are exactly the lowest-QPB
depth-2 nodes (all nodes with QPB below ~5548, 40 of the 225, 18%), and those are also the nodes
that carry most of the QPB work (~70%), so a depth-2 pass alone removes only ~30% of the work; a second
DNN pass at depth 3 under the 40 survivors fathoms ~98% of their 1096 children. The proposed hybrid
(DNN at depth 2, DNN at depth 3 under survivors, QPB or one more DNN level below) is estimated at
**2.5-5 CPU-days with the present Python code, ~1-2 CPU-days with a C/LAPACK implementation**,
against **30-120 CPU-days QPB-only**.

Everything is in `/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn30/`
(`batch30.py` node runs, `analysis.py`/`analysis.log` the cost model, `summarize.py`/`summary.log`
percentiles, `itercost.py`/`itercost.log` per-iteration profile, `results_depth*.jsonl` per-node records
incl. certified-lb histories, `est/` solver estimator logs, `n30_*.state` frontiers, `*.log` solver runs,
`root30.py`/`root30.log` the root DNN run). The bound code is the previous `dnn/dnn.py` + `dnn/qaputil.py`,
unchanged. Nothing under /home/user/qap was modified; every job ran alone on one core
(OMP/OPENBLAS/MKL_NUM_THREADS=1), sequentially, with per-node JSONL output and skip-on-rerun.
Labels: **[measured]** = computed here; **[estimated]** = solver's importance-sampled Knuth estimator
(`bnbdrv -e`) or an extrapolation from measured numbers, basis stated; **[recalled]** = from memory, unverified.

## 1. Setup

* Instance nug30, optimum 6124; objective values are even, so a node is fathomed when the certified
  lower bound exceeds 6122.05 (the solver's own tolerance 1.95, as printed).
* Root sanity [measured]: QPB root bound 5332.33 (`bnbdrv -p PARAM/default.par -d 0`, 0.16 s; rule2 gives
  the same 5332.33). DNN root (m = 30, lifted 901, V 842 columns): section 2.5.
* Frontiers [measured], `bnbdrv -p PARAM/rule2.par -u 6124 -d D -o state` (rule 2 = cheap branching rule,
  150 FW iterations; the solver uses the nug30 symmetry file, so the root has 9 children, not 30):

| depth | nodes opened | QPB fathomed | children eliminated | QPB-open (branched) | own-QPB percentiles 0/10/50/90/100 | pass time |
|---|---|---|---|---|---|---|
| 1 | 9 | 0 | 0 | 9 | 5406 / - / 5478 / - / 5643 | 0.04 s |
| 2 | 225 | 0 | 0 | 225 | 5474 / 5528 / 5621 / 5773 / 5979 | 0.77 s |
| 3 | 6120 | 19 (0.3%) | 0.5% | 6101 | 5534 / 5640 / 5752 / 5922 / 6122 | 19.8 s |

  (depth-3 nodes have 164,007 depth-4 children.) With the real top-of-tree rule (`default.par`, rule 4
  with 30-row lookahead) the depth-2 frontier is also 225 nodes (branching on facility 0 then 23/6
  instead of 6 then 21/4), median QPB 5645, max 5810, and the pass takes 89 s; the rule2 frontier was
  used for the DNN statistics, as instructed. Node QPBs are recovered from the children's inherited
  `bound` field exactly as in the nug20 study.
* DNN bound: facially-reduced DNN relaxation (Zhao-Karisch-Rendl-Wolkowicz / Oliveira-Wolkowicz-Xu)
  solved by the OWX ADMM; the certified bound is the Lagrangian dual value of the current multiplier
  (projection of V'ZV onto the negative semidefinite cone plus the clipped linear part), valid at every
  iteration; ADMM stops as soon as the certified bound exceeds 6122.05, or at convergence
  (primal residual and step < 1e-5), or at 4000 iterations / 900 s. The reduction of a node to an
  m x m subproblem with linear term and constant was verified on random completions for every node.

## 2. DNN results

### 2.1 Depth 2 (m = 28, lifted 785, V 730 columns) [measured]

* **Uniform random sample of 36 of the 225 QPB-open nodes (seed 0): 34 fathomed (94.4%; 95% Clopper-Pearson
  interval 81-99%).** Misses: {4->29, 6->7} QPB 5514.8, DNN value 6095.9 (converged, certified lb 6095.4);
  {4->7, 6->14} QPB 5542.8, DNN value 6109.1 (lb 6108.6).
* **The 8 lowest-QPB nodes of the frontier (QPB 5473.9-5504.2), sampled deliberately: 0 of 8 fathomed**,
  converged DNN values 6034.9-6085.2 (certified lb within 0.5 of the value), 330-500 s each.
* Over all 44 distinct nodes the outcome is a clean function of QPB: every node with QPB <= 5542.8 misses,
  every node with QPB >= 5554.1 is fathomed. Threshold T ~ 5548. **40 of the 225 depth-2 nodes
  (17.8%) have QPB < T, so the predicted depth-2 fathom rate over the whole frontier is 185/225 = 82%**
  (lower than the random sample's 94% because the random sample under-represented the tail).
* DNN value minus QPB: 480-580 at depth 2 (mean 569 at the missed nodes); the DNN closes almost all of
  the QPB gap (QPB is 320-650 below 6124) but not the last 40-90 at the lowest-QPB nodes.
* Cost (one core, 0.135 s per ADMM iteration): fathomed nodes **mean 128 s, median 101 s, range 67-375 s
  (500-2750 iterations)**; missed nodes run to convergence 307-496 s (2300-3450 iterations), mean 382 s.
  Nodes whose DNN value is within ~10 of 6124 need 2000-2750 iterations to certify (the certified bound
  lags the primal objective: the primal stabilises within ~300 iterations, the dual climbs ~5-10 per
  50 iterations near the end). Percentiles over the random 36: cpu 67/75/89/104/144/293/375 s
  (0/10/25/50/75/90/100), iterations 500/575/650/750/1038/2125/2750.
* Full DNN pass over the 225-node frontier [estimated from the above]: 185 x 128 s + 40 x 382 s
  = 39,000 s = **10.8 CPU-h** (Python).

### 2.2 Depth 3 (m = 27, lifted 730, V 677 columns) [measured]

| sample | n | fathomed | QPB range | mean cpu | iterations |
|---|---|---|---|---|---|
| uniform random of the 6101 QPB-open | 6 | 6 | 5619-5981 | 70 s | 350-1050 |
| children of the 2 random-sample misses | 8 | 8 | 5614-5779 | 88 s | 550-1150 |
| children of the 8 lowest-QPB misses | 6 | 5 | 5559-5753 | 215 s | 700-3950 |

* 19 of 20 fathomed. The one miss, {4->5, 6->13, 21->20}, QPB 5558.7, DNN value 6110.0 (converged).
  The lowest fathomed depth-3 QPB is 5574.8 (certified 6122.29 after 1850 iterations, 218 s), so the
  depth-3 threshold is ~5560-5575, again DNN ~ QPB + 500-550.
* The 40 predicted depth-2 survivors have 1096 QPB-open depth-3 children (QPB percentiles
  5534/5612/5649/5691/5845); only **25 of them have QPB < 5567**, i.e. a depth-3 pass under the survivors
  should fathom ~97-98% of them. DNN cost per depth-3 node under survivors is higher than the random
  mean (lower QPB -> smaller margin): ~120 s average over the three samples, 160 s for fathomed
  children of the lowest nodes.
* Depth-3 pass under the 40 survivors [estimated]: 1096 x ~120 s = **~37 CPU-h** (Python). A depth-3
  pass over the entire 6101-node frontier would be 6101 x ~80 s = 136 CPU-h = 5.7 CPU-days, i.e. worse
  than depth-2-then-depth-3 but still far cheaper than QPB-only.

### 2.3 Depth 1 (m = 29, lifted 842, V 785 columns) [measured]

One node ({6->7}, QPB 5459.3): DNN value 6052.2, converged after 2850 iterations, 518 s (0.172 s per
iteration); not fathomed. Depth-1 QPBs are 5406-5643, so with DNN ~ QPB + 560-590 only the single node
with QPB 5643 ({6->0}) could possibly be fathomed at depth 1. Depth 1 is not the place for the bound,
as at nug20 (1 of 6 there).

### 2.4 Validity

Same machinery as the nug20 study (where 16 exact comparisons passed); the reduction to each
subproblem was re-verified numerically on every node. Certified margins over 6122.05 here were at
least 0.08 (node {4->10, 6->8}: lb 6122.13), far above floating-point noise (~1e-9 relative).
Early-stopped certified values understate the true DNN value, often by 30-230.

### 2.5 DNN root value of nug30

[measured] Certified DNN root lower bound **5949.31** (primal value 5949.42) after 2700 ADMM iterations,
604 CPU s at 0.22 s per iteration (stopped by the 600 s cap; the certified bound was still rising by
~0.1 per 500 iterations, so the DNN root value is ~5949.4-5949.5). The gap to 6124 is 174.6 (2.85%),
versus 791.7 for the QPB root. The 2000 solve's root QPB gap was essentially the same 5332 -> 6124.
Literature DNN root for nug30 [recalled, unverified]: roughly 5950 (ZKRW-type SDP/DNN bounds for nug30
are quoted around 5940-5950 in the DNN/ADMM papers); I did not look it up.

## 3. The alternative: QPB subtree cost at nug30

The solver has its own importance-sampled Knuth estimator (`bnbdrv -e`, the one used for nug30 in 2000;
dives are weighted by relative-gap^1.5, exact through depth 1, parameters in `est/estimate.par`). It
was run from a scratch working directory (SYM symlinked) so nothing in the repo was touched. All with
`default.par` and incumbent 6124, i.e. the production configuration.

* **Whole nug30, QPB-only [estimated]:** run A, 400 dives (3.4 min): 2.38e10 nodes (+-1.1e10), 1.04e7 CPU s
  = 120 CPU-days (+-68, one standard error). Run B, 1500 dives (13 min): 1.39e10 nodes (+-0.35e10),
  4.90e6 CPU s = 57 CPU-days (+-13). The two agree within their errors; Knuth estimates of heavy-tailed
  trees are usually *under*-estimates when few dives reach the deep part. Time spent at depth <= 2 is
  ~90 s, at depths 3-4 ~6,000-9,000 s; essentially all the work is at depths 6-11 (`est/full_nug30*.log`).
* **Subtrees of specific depth-2 nodes [estimated]**, estimator run on the reduced m = 28 instance with
  the incumbent 6124 - const (400 dives each, 1.5-2.5 min):

| depth-2 node | QPB | DNN outcome | est. subtree nodes | est. subtree CPU s |
|---|---|---|---|---|
| {4->29, 6->7} | 5514.8 | miss (6096) | 1.24e8 (+-19%) | 56,300 (+-19%) = 15.6 h |
| {4->7, 6->14} | 5542.8 | miss (6109) | 5.7e7 (+-21%) | 27,700 (+-17%) = 7.7 h |
| {4->15, 6->14} | 5594.8 | fathomed | 1.2e7 (+-13%) | 8,000 (+-7%) = 2.2 h |
| {6->6, 21->2} | 5650.7 | fathomed | 7.7e5 (+-17%) | 1,070 = 0.3 h |
| {6->6, 21->24} | 5784.7 | fathomed | 9.8e3 (+-7%) | 74 |

  log10(subtree CPU s) is linear in QPB (slope -0.0109 per unit, residuals < 0.25 dex). Summing the fit
  over the 225 frontier QPBs gives a third whole-tree estimate of 2.9e6 s = **33 CPU-days**, and says
  that **the 40 lowest-QPB nodes (the DNN survivors) hold ~70% of the QPB work** and the 185 nodes the
  DNN fathoms hold ~30% (8.7e5 s = 241 CPU-h). The lowest node (QPB 5474) is predicted at ~1.4e5 s = 40 h.
* Average QPB subtree below a depth-2 node: 4.9e6-1.04e7 s / 225 = **22,000-46,000 s (6-13 CPU-h)**
  [estimated]; below a depth-3 node ~800-1,700 s on average.
* Historical figures [recalled, unverified]: the 2000 solve (Anstreicher, Brixius, Goux, Linderoth)
  is usually quoted at about 11.9 billion nodes and about 11,900 CPU-days on year-2000 hardware
  (the task statement says ~7 billion; I cannot resolve the difference from memory). The node estimates
  here (1.4-2.4e10) are of the same order even though the branching rules differ; per-node speed on this
  machine with this code is ~0.3-0.75 ms at depths >= 8 (from the estimator's per-level times), ~60 ms at
  depths 3-4 (rule 4 lookahead) and ~0.4 s at depth 2.

## 4. Break-even at nug30

Per depth-2 node: DNN costs **128 s** to fathom, and the subtree it removes is worth **thousands of
seconds** for any node above the threshold (8,000 s at QPB 5595, 1,070 s at QPB 5651, 74 s at 5785:
the DNN is dearer than the QPB subtree only above QPB ~5790, the top ~10% of the frontier). Averaged
over the frontier the DNN pass costs 10.8 CPU-h and removes ~240 CPU-h (per-node fit) to
~0.3 x (57-120 CPU-days) (dive estimates) of QPB work. **Break-even is not in question at nug30; the
depth-2 DNN pass pays for itself 20-100 times over.** The real issue is that the survivors are the
expensive nodes, so a depth-2 pass alone cuts the total by only ~30%; the bound has to be applied again
at depth 3 (and possibly 4) under the survivors.

Why the picture flips between nug20 and nug30: DNN cost per node scales like (m^2)^3 (nug20 depth 2:
8 s; nug30 depth 2: 130 s, 16x), while the QPB subtree below a depth-2 node went from ~2 s to
~20,000-50,000 s (10,000x).

## 5. Proposal: DNN-at-depth-2/3 + QPB hybrid for nug30

Work estimate, single core, incumbent 6124 known (all [estimated] from the measured per-node costs and
the subtree estimates above):

| stage | nodes | cost (Python ADMM) | cost (C/LAPACK, sec. 6) |
|---|---|---|---|
| QPB frontier to depth 2 (default rule) | 235 | 90 s | 90 s |
| DNN at depth 2 on all 225 nodes | 225 | 10.8 CPU-h | ~2.5 CPU-h |
| QPB children of the ~40 survivors | 1096 | ~0.5 CPU-h (QPB) | same |
| DNN at depth 3 on those 1096 children | 1096 | ~37 CPU-h | ~8 CPU-h |
| below the ~25 depth-3 survivors: either QPB subtrees | 25 | ~70 CPU-h (if ~10,000 s each, from the depth-2 fit one level down; uncertain by 3x) | same |
| ... or DNN at depth 4 (m = 26, ~650 children, ~70 s each) + QPB below the few depth-4 survivors | 650 | ~13 CPU-h + QPB remainder (small) | ~3 CPU-h + remainder |
| **total** | | **~60-120 CPU-h = 2.5-5 CPU-days** | **~15-30 CPU-h = 0.6-1.3 CPU-days** |
| QPB-only (solver estimator, two runs + per-node fit) | | **33-120 CPU-days** | |

So roughly **10-40x less CPU** than QPB-only with the Python bound, 30-100x with a C implementation,
and the DNN stages are embarrassingly parallel node-by-node (as the QPB tree was in 2000). With a
C implementation the whole hybrid is a day or two on one modern core, or an hour on a workstation.

Recommended depth: **2 first, then 3 under survivors** (not 1: depth-1 DNN values sit 70-120 below the
optimum; not 3 everywhere: 5.7 CPU-days of DNN for the same effect). Practical details: (i) the node
QPB predicts the DNN outcome almost perfectly (threshold ~5548 at depth 2, ~5565 at depth 3), so nodes far
below threshold can be skipped (saving the 380-500 s convergence runs) and nodes far above can use a
lower iteration cap; (ii) the stopping rule must be the certified dual bound, not the primal objective;
(iii) the hybrid needs the incumbent 6124 (an early-stop test), as the QPB solver does.

## 6. Cost of an efficient implementation

Per-iteration profile at m = 28 [measured, `itercost.log`, numpy + OpenBLAS LAPACK, one thread]:
0.135 s = eigh of the 730x730 matrix 0.066 s (49%) + the two products V'(.)V 0.026 s + V(.)V' 0.024 s
+ U diag(w) U' 0.012 s (46% together) + clipping/Z-update/norms 0.007 s. The certified bound costs one
more eigh+products (0.13 s) every 50 iterations (2% overhead). m = 29: 0.172 s; m = 27: 0.113 s.
Iterations needed [measured]: fathomed depth-2 nodes median 750 (500-2750), depth-3 median 600-800;
misses to convergence 2300-3950.

What C/LAPACK buys [estimated]:
* The dense eigh is already LAPACK dsyevd; a plain C port gains little there. Two structural savings:
  (a) V = [[1/sqrt2, 0], [1/(m sqrt2), Vb kron Vb]] with Vb (m x (m-1)): V'MV and VRV' are four mode
  products on an m x m x m x m tensor, 4 m^5 ~ 70 MFlop instead of two 785x785x730 GEMMs (0.9 GFlop):
  the 62 ms of products become ~5-10 ms. (b) Only the positive eigenpairs of W are needed, and their
  number drops to the (small) rank of the DNN solution as the ADMM converges: dsyevr/dsyevx with a
  half-line range still pays the 2/3 n^3 tridiagonalisation (~0.26 GFlop, 10-15 ms) but saves the
  back-transformation; a warm-started block Lanczos/LOBPCG (k ~ 20-50 vectors, O(k n^2) per iteration)
  can bring the R-step to ~5-10 ms once the active rank is small.
* Realistic per-iteration cost at m = 28: **25-40 ms (3-5x faster than numpy)**, i.e. 20-30 s per
  fathomed depth-2 node, 60-140 s per miss run to convergence, 15-25 s per depth-3 node. Numbers in the
  table above use 4x.
* Iteration count is the other lever: the primal objective settles within ~300 iterations while the
  certified dual bound needs 600-2750; a better dual recovery (or a few rounds of dual-ascent polish on
  the current Z) could cut the iteration count by 2-3x. Not tried here.
* The DNN pass is independent per node and memory-light (five 785x785 dense matrices, ~25 MB), so it
  parallelises trivially.

## 7. What is measured vs estimated vs recalled, and loose ends

* Measured: all DNN bounds, fathom outcomes and CPU times (36 random + 8 lowest depth-2 nodes; 20 depth-3
  nodes; 1 depth-1 node; root 5949.3), the frontiers and their QPBs, the per-iteration profile, the solver
  frontier pass times.
* Estimated: every QPB subtree/whole-tree cost (solver's own importance-sampled estimator, standard
  errors quoted; such estimates can be low by a factor of 2-3 when deep dives are rare), the 82%
  whole-frontier fathom rate (threshold model from 44 nodes, applied to the 225 frontier QPBs), the
  depth-3 survivor count (25), the C/LAPACK speed-ups, and all hybrid totals.
* Recalled, unverified: literature DNN root for nug30 (~5950), the 2000 solve's node count and CPU-days.
* The depth-2 frontier used is the rule2 one; the real search (rule 4 at the top) has a different but
  equally sized frontier (225 nodes) whose QPBs are higher on average (median 5645 vs 5621, max 5810), so
  the fathom rate there should be at least as good. Depth-3 and depth-2 cost estimates ignore the
  QPB cost of re-bounding the children (small: 1096 x ~0.07 s).
* Runs that hit the 900 s / 4000-iteration cap: none; the longest node was 518 s (depth 1).
* The solver's `-e` estimator crashes (one-byte overflow in `set_sym_file`, driver.C:278) when the
  instance path has no directory component; absolute paths avoid it. Not fixed (repo untouched).
* Timings: Python 3 + numpy (OpenBLAS, 1 thread) on this container; other benchmarks were running on the
  other cores throughout, so timings may be ~10-20% pessimistic.
