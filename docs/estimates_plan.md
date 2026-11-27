# Settings and estimates before the next large run (Nov 18 2026)

Nathan, Nov 18, after tai60b: tho40 is parked. First (1) learn from the recent runs and adjust the
settings: the laptop's sko42 zone test, rule 2 in a DNN solver, the use of strong branching. Then use
AWS for estimates: (2a) test the settings hypotheses, (2b) estimates we can trust, with the thesis's
breadth-first top levels, dives biased deep, and the fleet's parallelism, through the coordinator if
need be.

The evidence behind this page:
- `experiments/tai60b_postmortem/`: where tai60b's tree went and why its estimate was 11 times low;
- `experiments/line_signals/`: which signals pick a branching line;
- `experiments/rule9_bottom/`: rule 2 at the bottom of the tai-b trees (its A/B, `tai50b-r9`, was
  launched Nov 18 about 14:40 UTC);
- the laptop's sko42 runs at `cut.zone` 0.010 and 0.025 (`infra/laptop/runs/sko42*/`).

## 1. What the runs say about the settings

**The cost is in a few branchings near the top** (tai60b). Its tree is two chains, the best known
solution's path and its mirror. Two chain branchings decide two thirds of the CPU:
- depth 3, rule 9, line r6: the 14 children that branch off it hold 44% of the CPU;
- depth 7, the first depth of rule 2, c54: the 18 children that branch off it hold 35% of the nodes
  and 24% of the CPU.

16 chain nodes at depths 1-8 each have 49,000 to 175,000 nodes below them.

### a. `cut.zone` (sko42, the laptop)

A node's cut phase leaves a pool its children inherit, even when it does not fathom the node:
- children of nodes at gap 100-120 fathom on their first solve 99.5% of the time;
- just outside the zone (gap 158-170), 60%.

At 8.6 h, zone 0.025 against 0.010 had depth-3 first-solve fathoms at 85% against 72%, and fewer
open nodes at depths 2-3. That supports H1; H2 (less CPU in all) is not yet decided.

Two checks to add, both from data we have:
- **tai60b's off-chain children at depths 3 and 7.** Were they inside `tai60b.par`'s zone? Did their
  children inherit pools? (The collected trace has `cut` and `pool_out` per node.)
- **sko42's run 2 against H3.** Watch CPU per bounded node and whether its open nodes drain.

The zone also changes what the estimate must sample: more phases near the top, more fathoms below.

### b. Rule 2

Rule 2 is a QPB rule, and QPB sits 13-29% below the cutoff at every DNN node of these trees. So:
- **Its count is always 0.** It counts children that bound + U_ij fathoms, and none qualify.
- **What actually decides is its tie-break, the largest U column sum.** That ranks lines no better
  than chance (`line_signals`, section 1: Spearman +0.06 against the true subtrees).
- **It branches on the same location at almost every node of a depth:**
  - tai60b: location 54 at all 97 depth-7 branchings;
  - tai50b: location 20 at depth 7 and location 33 at depth 11.
- **On tai60b, the depth-7 chain node is where the bottom block starts, and it leaves more branching
  children than any other chain depth** (9 per chain).

This is suggestive, not proof, since depth is mixed up with the rule. There is also no DNN form of
rule 2: the certificate's child bounds are flat. The replacements, cheapest first:
1. **Rule 9 at every depth.** `tai50b-r9` decides it for tai50b. Rule 9 costs about nothing extra.
2. **A DNN lookahead at the heavy nodes** (c, below).
3. **A cheap DNN signal for the bottom.** None found yet: the marginal score ranks lines at
   Spearman -0.07.

### c. Strong branching

Today's strong branching is rule 4 and rule 9's arbiter:
- QPB bounds the children of one line (rule 4), or of two lines (rule 9);
- the shortlist comes from the DNN marginals.

`line_signals` measured the alternative: solve each candidate line's children with a short warm
DNN solve under the node's cut pool, and branch on the line whose children fall least short.
- It ranks lines by their true subtrees **once the solves fathom some children**: 30 iterations
  at m <= 33, about 150 at m = 39 (tai40b's root).
- With too few iterations it picks wrong lines: tai40b grew 13%.
- Its cost is lines x children x iterations: about 18 CPU-hours per node at n = 60 for 8 lines at
  150 iterations.
- That is far too much for a whole tree, but a fraction of a percent at tai60b's 16 chain nodes
  and its depth-3 and depth-7 branchings.

So the opportunity is **strong branching where the subtree is large**. Two pieces are needed:
- **A test for "large"** that the node can compute. Candidates:
  - depth;
  - the node's shortfall after its phase, against its depth's typical value (the bound lab's
    survive rule: right on 270 of 286 off-path nodes);
  - the shortfall of its parent's other children.
- **A lookahead that fits n:**
  - an iteration cap that grows with n, or a per-child stop (iterate until the child fathoms or
    its bound stalls);
  - 4-6 lines a side;
  - the chosen line's children started from their lookahead solves.

The one untested step is whether a longer cap picks the right line at n = 60. The test for that is
`line_signals`, Verdict 4: tai60b's 43:58, four lines, 150 iterations, about 2.5 h on the cloud box.

## 2. Estimates on AWS

### What the estimator does today (`latest/estimate.C`, thesis ch. 4)

- **Knuth's dives with importance sampling.** A dive takes child i with probability g_i^q / sum g^q,
  where g is the child's relative gap. q is `score_exponent`, 1.5 in `estimate.par`.
- **The child's gap comes from QPB.** By default it is the parent's prediction (bound + U); with
  `QAP_EST_CHILD_QPB=1`, the child's own QPB bound. `est.xbar` switches to weights from the DNN
  marginals.
- **The thesis's breadth-first top.** The top `start_dive_level` levels are run whole, and each dive
  starts at a level node drawn by the same g^q rule. This is the nug30 design: NBFS = 3, q = 1.5,
  10,000 dives.
- **`-e -r STATE` starts the dives from a saved state.** One `QAP_EST_SEED` per process makes the
  parallel processes take different dives.

Two limits matter in a DNN solver:
- **The g^q bias has lost its information.** The QPB gaps of a node's children barely differ, so the
  weights are nearly uniform. Meanwhile 93-99% of the nodes at a depth are fathomed (tai60b's table
  by depth), most of them on their first DNN solve, so a dive almost always ends at the next level. That is exactly what happened in tai60b's dives (none
  went more than 2 levels below its start).
- **Picking the frontier node by g^q mixes strata.** Every frontier node has its own subtree, and
  the variance between frontier nodes counts against every dive.

### 2b. The design to try

1. **The breadth-first top, exact** (the thesis's NBFS; `qapctl init --frontier-depth k` already does
   it):
   - run every node through depth k and keep every node at depth k+1 with its parent's state, as a
     solve's units start;
   - k = 3 on tai60b is 588 bounded nodes and about 18 CPU-hours;
   - the root alone is 1.5 h on one core.
2. **Stratify instead of sampling the frontier.**
   - Every depth-(k+1) node is a stratum. It is bounded once; a fathomed one is exact.
   - Every other one gets its own dives.
   - The whole-tree estimate is the top plus the sum of the strata means. The dives go where
     the variance is (Neyman allocation, from a first round of a few dives each).
3. **Bias toward survivors, exactly.** At each node of a dive, bound all its children, as the search
   would, and continue into a surviving child. Two choices for which one:
   - by g^q over the survivors' own DNN bounds, which are now known;
   - uniformly.

   The estimate stays unbiased (weight = survivors / p). A dive no longer ends on a child that the
   search fathoms in one solve: it ends where the tree does. This is Purdom's partial backtracking
   (1978), and the bound lab's "near/far" strata did a hand version of it. It costs one level's
   children per step instead of one child, mostly cheap first-solve fathoms.
4. **Share work between the dives of a stratum.** Memoize the bounded nodes within a process, so a
   node two dives pass through is bounded once. This is still unbiased, and saves the repeated top
   of the dives.
5. **Check against the truth.**
   - Exact subtrees for a random handful of strata (the solve's own units, under a node budget) on
     each instance.
   - Every estimate reports its dives' depth profile, and how often the deepest sampled level
     reached the strata's deepest.

**Before any of this runs on AWS, choose the design offline on trees we know.** A replay program
samples dives from a finished run's node trace instead of solving nodes. The trace has every node's
children, bounds (QPB and DNN), outcome and CPU, so the replay can apply each design to the real
tree:
- the 2000 design;
- stratified;
- stratified with survivors;
- each at k = 2, 3, 4;
- each at q = 0, 1.5, 3.

For each one it measures the error and the CPU per dive against the known totals. The trees:
- tai60b-261117 (331k nodes);
- tai50b-261116 (8.2k);
- the laptop's tai35b and tai40b trees (`experiments/line_signals/runs/`);
- nug30's EC2 traces (`infra/phase1/runs/nug30-c7a/`, `nug30-trace-a/`).

The replay answers how many dives, and which design, give ±20% on a tree like tai60b's, for no AWS
cost. It can only follow the recorded tree, so it says nothing about settings: that is 2a.

### 2a. Settings hypotheses

A whole-tree estimate's error swamps a 10-20% settings effect, unless both arms are estimated from
the same frontier with many dives. Two ways to compare settings:
- **Paired estimates.** The same depth-(k+1) strata under settings A and B, the dives of each
  stratum under both. The difference is estimated stratum by stratum.
- **Paired exact subtrees, the decisive form.** A stratified random sample of strata, solved to the
  end under A and under B, as in sko42's matched-node A/B and `giveup_ab`. The solve units do this
  today, given two params files.

The hypotheses, with the evidence for each:
- **H-rule2: settled Nov 18.** `tai50b-r9`: -19% nodes, -12% CPU; rule 9 now replaces rule 2 in every recipe's
  bottom block (`experiments/rule9_bottom/README.md`). Was: rule 9 (or the DNN lookahead) in place of rule 2 below depth 6 (1b; `tai50b-r9`
  first).
- **H-zone:** a wider `cut.zone` (1a; sko42 on the laptop).
- **H-strong:** a DNN lookahead at heavy nodes with a cap that grows with n (1c; the tai60b
  43:58 test first).

### 2biii. The coordinator: estimate runs

**Built, Nov 18** (step 2), as the replay revised it:
- the estimate is the start of a solve;
- survivor dives in `latest/estimate.C` (`QAP_EST_SURV`): unbiased on nug15, within 2% for every
  seed and weight;
- dives runs and `qapctl dives`, `estimate`, `estimate-flow` and `worker --follow`
  (`coord/dives.py`, `coord/README.md`);
- `QAP_ESTIMATE` in `node.sh` and start.sh's prompt (`infra/phase1/README.md`);
- `ESTIMATE=` for `infra/laptop/solve.sh`;
- `tests/dives_test.py`: nug17's estimate within 15%, and the paused solve resumed to the
  serial totals.

The design below was written before the replay.

A run kind `estimate` beside `solve`. The design:
- **init:** the frontier pass as today (`--frontier-depth k`). Its nodes become strata, not units.
- **Units:** (stratum, seed, dives, budget T). The worker runs `bnbdrv -e -r stratum.state` with the
  seed and a wall budget, and writes one JSON line per dive: depth reached, and nodes, CPU and
  iterations by depth. Today's dive history has only the totals.
- **Allocation and stop:**
  - a first round of a few dives per stratum, then Neyman allocation from the observed variances;
  - stop at a wall-clock budget, or when the CI reaches a target;
  - `qapctl stop` ends it early.
- **Collect:** `qapctl export` gives dives.jsonl. `experiments/estimates/estimate.py` (new) combines
  the strata: the estimate, a bootstrap CI over dives within strata, the depth profile, and the
  status page's numbers.
- **Anchors:** some strata as ordinary solve units with a node budget, for item 5 above.

Most of the parts exist (the frontier pass, units, budgets, the store, export, Spot). New code:
- the unit kind;
- the allocation loop;
- the per-dive JSON lines;
- survivor dives and the memo in `estimate.C` (`bnb()`'s dive path).

Tests: the moto backend and a local run, as for the solve.

## Order of work, and what each step decides

1. **Replay study** (offline, the cloud box, no AWS). It picks the design and the number of dives.
2. **Estimator changes** (survivor dives, memo, per-dive JSON) and the coordinator's `estimate`
   kind, tested locally and with moto.
3. **A calibration run on AWS against a known tree.** tai50b: its truth is 8,216 nodes and 133
   core-hours, and a 1-2 h estimate on 64-128 slots costs a fraction of that. Or tai60b with a
   2-hour budget, which is the hard case.
4. **The settings A/Bs (2a)**, as paired strata, on the instance the next run is for.
5. **The estimates of record** under the settings chosen in 4 (Nathan, Nov 18): tho40, sko49, wil50
   and tai60b. tai60b is the reference: its truth is known (331,149 nodes, 3,597 core-hours), and the
   settings are expected to improve on it substantially.

Step 1's first results (`experiments/estimator_replay/`, Nov 18):
- today's design gives 5-9% of tai60b's CPU at any number of dives;
- survivor dives with adaptive strata are within 20% three times in four, at 30% of the tree's CPU.
An estimate run should therefore be the start of a solve, breadth-first for a fixed wall-clock budget,
with survivor dives on its open nodes.

Running beside this: the `tai50b-r9` result (rule 2, tonight), the tai60b 43:58 lookahead test (the
cloud box, 2.5 h), and sko42's zone runs on the laptop.
