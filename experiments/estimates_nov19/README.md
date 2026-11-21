# The Nov 19 estimates: tho40-est2 and tai50b-est (worker groups, readouts)

The first two estimates on worker groups (`infra/phase1/README.md`, Starting a run), launched from specs
(`infra/phase1/specs/tho40-est.spec`, `tai50b-est.spec`) and read from their readouts
(`git show origin/phase1-readouts:runs/RUN/readout.md`; `coord/readout.py`). Both ran without a person
in the loop from launch to the end: coordinators on demand, workers in Spot groups that the
coordinators started after the root and stopped at the end. The records with the dives come with
`infra/phase1/collect.sh RUN` (`runs/RUN/dives/`).

## tho40-est2-261119: tho40c, with the first estimate's diagnosis applied

Settings (`experiments/tho40_est_diag`): Knuth dives `k1.5`, 900 s of dives a unit counted after setup,
`QAP_EST_MIN_DIVES=5`, sampling by depth `2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03`; 1 hour of solving on
135 slots (a `c7a.2xlarge` coordinator and two Spot workers).

| | first (tho40-est, Nov 18) | tho40-est2 |
|---|---|---|
| whole tree, CPU | 4,505 h (3,502-5,553) | **5,728 h (5,182-6,255)** |
| c7a CPU-days | 188; 190 (140-250) with the deep tail and the bound lab's comb arm | **239 (216-261)** |
| nodes | 1.87M | 2.07M (1.75-2.45M) |
| jackknife standard error, CPU | 16% | 7% |
| dived strata with one dive | 828 of 1,108 | 0 of 2,575 |
| estimator CPU: setup, dives | 68.3 h, 26.3 h | 156.4 h, 216.8 h |
| the 10 largest dives' share of the rest's CPU | 37% | 13% |
| readout warnings | HEAVY TAIL, FEW DIVES, SETUP-BOUND | none |

Consistent with the first estimate and with the bound lab's comb arm (215 days, se 82), and three
times tighter. Its cost: 117 CPU h of solving and 373 of estimator, 9% of the tree. The solve itself:
about a day on four `c7a.16xlarge`.

**Read high by about a fifth: 196 c7a CPU-days, not 239** (below, "Knuth dives from cold strata"). 44%
of its strata had the first step of their dives bounded cold, and the estimate counted that CPU as the
tree's. With those strata imputed from the warm ones at their depth, the rest is 4,585 CPU h, not 5,610
(nodes 2.06M, as before); the cold strata branched a little more on their own (warm) bound, so the truth
is likely a little above 196. Its 95% band (216-261) did not cover this: the bootstrap and the
jackknife resample the dives, not the way they were made.

## tai50b-est-261119: a known tree, and a poor estimate that said so

tai50b under its current recipe (`recipes/tai50b.par`, rule 9 at the bottom), the tree of `tai50b-r9`
(6,629 nodes, 419,520 CPU s = 116.5 h). Survivor dives (`b = 0.4`), 600 s a unit, 1 hour of solving.

| | estimate | truth (tai50b-r9) |
|---|---|---|
| tree | 4,076 nodes, 77.2 CPU h | 6,629 nodes, 116.5 CPU h |
| ratio | 0.61 of the nodes, 0.66 of the CPU | |
| rest beyond the solve | 2,959 nodes, 44.4 CPU h, jackknife se 13.7 h | 5,512 nodes, 83.7 CPU h |

Its readout said not to trust it: HEAVY TAIL (the largest dive 42% of the rest's nodes), FEW DIVES (7
of 7 dived strata one dive each), SETUP-BOUND, WIDE (se 31%); its 5-95% band collapsed to a point.

- **The shallow part is exact.** The solve's part equals the truth at depths 0-3 (1, 50, 490, 576), and
  setup, which bounds every sampled open node and its children, gives depth 4: 282 (truth 282), depth 5:
  783 (782), depth 6: 663 (630).
- **The deep chain is one sample path.** The tree hangs a long, thin chain off the 17 depth-4 nodes that
  branch, to depth 20 and more (150-185 nodes a level at 13-18). A survivor dive in tai50b bounds every
  child at each level, about 46 nodes of 88 s, an hour a level: a 600 s budget gives each stratum one
  dive, and the 7 dives reached depth 11 at most. Depths 7-11 read 344, 293, 203, 198, 193 against the
  truth's 484, 430, 378, 328, 200; nothing below 11.
- **Part-way it read 0.28.** The last units to merge were the 7 branched depth-4 nodes (282 = 275 merged
  + 7 pending; 17 = 10 + 7): the units that run longest hold the large subtrees. Hence SLOW PENDING
  (`coord/dives.py`, 7042ae7a).
- **CPU per node: mostly the runs, not the estimator.** Where the counts agree, the estimator's CPU per
  node is under the solve's: depth 4, 212 s against 264 (0.80); depth 5, 53 s against 88 (0.60); depth 6,
  38 against 68 (0.60). Node by node (the collected dives against tai50b-r9's `nodes.txt.gz`, matched by
  their fixes): all 282 depth-4 nodes are in both, and the estimator bounded them in 16.6 CPU h against
  r9's 20.7 (0.80; per node a median of 0.96: 139 nodes at 0.9-1.0, 67 at 0.6-0.8, the phase-heavy ones,
  the 7 that branched among them, near 0.7). The search's own part of this run against r9's, the same
  code on the same nodes: 0.87 (0.95 at depth 1, 0.80 at 2, 0.91 at 3). Both ran tai50b's recipe (r9's
  `tai50b_r9.par` differs from `tai50b.par` in one comment) on `c7a.16xlarge`s at 64 slots, the solver
  changed only in `estimate.C`, `settings.C` and the params reader's comment handling between them. So
  the same nodes cost 13-20% less CPU in one run than in the other, and the estimator's own share of the
  gap is at most a few percent. Two runs of one tree differ by this much: an A/B of settings must run its
  arms side by side, on the same instances at the same time, and compare them paired.

**What the next estimate of a tree like this should change:** Knuth dives (one child a level, about 88 s
a level, so many dives in the budget) with `QAP_EST_MIN_DIVES=5` and a longer dive budget; or no
estimate at all: tai50b's whole tree is 117 CPU h, and this estimate cost 70 (solve 32.8, estimator
37.5). For a small tree the solve is the estimate.

## The decision point (Nov 19)

The change agreed after tai50b-est: an estimate is worth its CPU only while it is a small share of the
tree it estimates, so the estimate-flow decides, once, whether to go on diving or to solve
(`coord/dives.py`, `decision()`; `coord/README.md`). A dives run marks a pilot: a random tenth of its
units (10 at least, half of a smaller run's), the first of its shuffle, so they are sent first. When
every pilot unit has merged, the estimate-flow compares the estimator CPU the whole dives run will take
(the pilot's mean per unit, times the units) with the rest of the tree as the pilot alone estimates it
(its open nodes standing for every other unit's, depth by depth). Above 15% (a spec's `ESTIMATE_DECIDE`)
the dives stop and the solve goes on to the end on the same workers; the readout says SOLVING, then
COMPLETE, with the decision's line. Otherwise the dives go on, and the estimate ends with its cost.

**Why a pilot, not the first units to merge.** A first version decided on the merged units, once a tenth
had merged. Replayed on the collected runs (`decide_sim.py`), it would have stopped tho40-est2 after 55
units: the first to merge are the small ones (a unit runs long when its open nodes branch), and they
read a rest of 42 CPU h against the 5,610 the whole dives run found. A pilot fixed in advance and waited
for whole is a random sample of units:

| | the whole dives run | a random pilot of a tenth (400 draws) | decides to solve |
|---|---|---|---|
| tai50b-est (259 units) | 37.5 estimator CPU h, 84% of the rest it found | share 74-116% (10-90%); its CPU 3.1 h (median) | in every draw |
| tho40-est2 (548 units) | 373 estimator CPU h, 6.7% | share 4.8-8.4%; its CPU 37 h | in none |

So tai50b-est would have gone back to solving after about 3 estimator CPU hours, and tho40-est2 would
have run as it did. `tests/dives_test.py --decide` checks the switch end to end on nug17 (the solve
finishes with the serial run's nodes on the worker that came back from the dives run).

## Knuth dives from cold strata (tho40-est2)

The estimator bounds every sampled open node (a stratum) and keeps it as its children's warm start for
the dives. Past `QAP_EST_WARM_KEEP` (8) strata a process dropped that warm start, for memory (a DNN state
is tens of MB at n = 60): for survivor dives that costs only the estimator's CPU (the estimate takes each
survivor's first bound from setup, made warm), but a Knuth dive's first step bounds a child for the
first time, so that bound ran cold, without the parent's DNN state and without its cut pool, and the
estimate counted it as the tree's. tho40-est2's units held a median of 10 strata, up to 35: 2,795 of its
6,307 strata (44%) were past the 8th of their unit. Their dives' first children, against the warm
strata's at the same depth:

| stratum depth | first child's CPU, warm strata | cold strata | ratio | the first child branched: warm, cold |
|---|---|---|---|---|
| 2 | 79 s (174 dives) | 164 s (216) | 2.07 | 40%, 42% |
| 3 | 19 s (7,975) | 66 s (7,301) | 3.45 | 8.2%, 12.8% |
| 4 | 8.3 s (8,197) | 38 s (149) | 4.55 | 1.3%, 4.0% |

The strata themselves are alike: their own bound, made warm in both groups, took 100 s against 106 at
depth 3 and branched 45.5% against 50.6%. Imputing the cold strata from the warm ones at their depth:

| stratum depth | rest as estimated: nodes, CPU h | the cold strata from the warm ones |
|---|---|---|
| 2 (85 strata, 36 warm) | 317,542; 881 | 496,334; 1,071 |
| 3 (4,737, 2,052 warm) | 1,594,071; 4,358 | 1,423,376; 3,152 |
| 4 (1,377, 1,316 warm) | 136,134; 359 | 137,405; 349 |
| all | 2,053,682; 5,611 | 2,063,050; 4,585 (0.82) |

With the solve's 117 CPU h: **196 c7a CPU-days, not 239**, in line with the first estimate (188; the same
flaw, on 16% of its 3,479 strata, so a smaller one; and rule 2 at the bottom then) and the bound lab's
comb arm (215, se 82). Depth 2's
correction rests on 36 strata and goes the other way; the depth-3 strata carry it.

**The fix** (`latest/estimate.C`, `spill_save`): past `QAP_EST_WARM_KEEP` a stratum's DNN state goes to a
file beside the dives' output, fp32 as a state file ships it, and is read back around each dive from
that stratum; its pool and X stay in memory. Past `QAP_EST_SPILL_MB` (1,000 per process; tho40's
largest unit needs about 350) a stratum falls back to the old way, and every stratum's line says which
it had (`"warm"`: `m`, `s` or `c`). On nug24 under the nug30 recipe, from a state saved at 12 CPU s, one
seed three ways (`tests/estimate_test.py`, spill): every state in memory, every one spilled, and the old
drop give the same 175.7 nodes for the rest, and a first child's CPU of 0.16, 0.15 and 0.43 s; the rest's
CPU 46.6, 45.6 and 81.4 s (the truth about 52).

## Predictor-weighted dives (`QAP_EST_CHILD_DNN`)

The second change agreed after tai50b-est: Knuth's dives as the thesis made them (ch. 4.3: one child a
level, drawn by g^q of a cheap bound computed for the purpose) with a predictor fit for a DNN tree. With
`QAP_EST_CHILD_DNN=CAP`, each child a dive may take is solved for CAP ADMM iterations warm from its
parent's state, as its own first solve would start (the parent's Z and conditioned Y, under its pool),
with no give-up and no cut phase; its shortfall (cutoff - lb, 0 if it fathoms) gives g, and the weights
are g^q with the run's uniform share mixed in, so the estimate stays unbiased. Its CPU is the
estimator's. A parent with no DNN state (a QPB-only depth) falls back to QPB's weights.

What it costs: CAP iterations for every child at every level of a dive, against one full bound of the
child taken. A short warm solve carries a fixed part (the allocation and the first full
eigendecomposition), so it is cheap against a full bound where those are long (the top levels, nodes
with a cut phase) and no cheaper where a warm child fathoms in a few seconds: tho40's depth-5 and 6
nodes take 10 and 3.5 s. It is meant for the chain trees (tai), where Knuth's dives by QPB do not
converge and survivor dives bound every child.

First runs: nug24 under the nug30 recipe (249 nodes, 65 CPU s), one stratum (the root), 60 s of dives,
two seeds an arm: by QPB 97 and 146 nodes (25-28 dives), by a 30-iteration probe 195 and 205 (9-11
dives), both noisy at this budget; the tree is four levels deep, too shallow to show what the probe is
for.

**tai35b: the probe finds the branching children one level down and no further.** tai35b's tree
(`experiments/line_signals/runs/tai35b_base*`: 1,017 nodes on two runs, 6,931 and 7,390 CPU s, to
depth 28) under its recipe (`experiments/line_signals/tai/tai35b.par`, a cut phase at every
unfathomed node), Knuth's dives `k1.5` from the root, an hour of dives after the root's bound (611 CPU
s), two seeds an arm, on this 4-core box:

| | dives | nodes / truth | CPU / truth | the dive went on past depth 1 | past depth 2 |
|---|---|---|---|---|---|
| by QPB | 116, 131 | 0.16, 0.15 (se 0.03) | 0.25, 0.23 | 28 of 247 (11%) | 0 of 28 |
| by a 30-iteration probe | 22, 28 | 0.16, 0.14 (se 0.08, 0.07) | 0.43, 0.28 | 20 of 50 (40%) | 1 of 20 |

The probe cost 58 CPU s on the root's 35 children (1.7 s a child at m = 34, with the pool; a depth-1
node's full bound averages 25 s). It draws the 3 depth-1 nodes that branch four times as often as QPB
does, but at depth 2 it picks the one child in 34 that branches no better than chance, and both arms
see 2 of the tree's 382 nodes at depths 3-5 or none, after an estimator spend of about 60% of the
whole tree's CPU. The likely reason is the recipe's cut phase: the probe measures the uncut shortfall,
and here the cut phase decides which children fathom. So `QAP_EST_CHILD_DNN` stays in, off by
default, as a negative result on a chain tree with cuts; the decision point above is what keeps such
an estimate from costing more than the solve.

## Tuning on paired subtrees (Nov 19)

Two estimates of a tree, one per setting, cannot separate a 10% effect: each is known to about 7%
(tho40-est2's jackknife). Solving the same open nodes to the end under each setting, side by side, and
comparing node by node takes the between-node spread out of the comparison. Built: a run from a sample
of a paused run's open nodes (`coord/sample.py`; `qapctl init --sample-from`; a spec's `FROM=` and
`SAMPLE=`, with `ARMS` and `ARMS_ENV`), and `paired.py`, which reads the arms' samples and node traces.

**How to draw the sample.** tho40's subtrees are uneven: of tho40-est2's warm strata at depth 3, a
median subtree of 0.03 CPU h, a mean of 0.40, the largest 1% holding 35% of the CPU (depth 4: 0.002,
0.05, 43%). A uniform sample of 200 depth-3 open nodes estimates their total to 29%. For a ratio
between arms, if an arm's effect varies by a share sigma from node to node, the ratio's relative
standard error is about sigma times the total's, and the cost of a node is its size, so the variance
at a fixed cost is least with pi ~ sqrt(size). With est2's warm strata as the population
(`sample_design.py`; the cost is the expected CPU an arm):

| depth | draw | n | CPU an arm | the total's rel. se | the ratio's rel. se, sigma 20% / 40% |
|---|---|---|---|---|---|
| 3 | uniform | 200 | 80 h | 28.6% | 5.7% / 11.4% |
| 3 | pi ~ sqrt(size) | 50 | 113 h | 18.8% | 3.8% / 7.5% |
| 3 | pi ~ sqrt(size) | 100 | 225 h | 11.5% | 2.3% / 4.6% |
| 3 | pi ~ size | 50 | 297 h | 11.3% | 2.3% / 4.5% |
| 4 | uniform | 200 | 10 h | 42.9% | 8.6% / 17.2% |
| 4 | pi ~ sqrt(size) | 100 | 35 h | 7.7% | 1.5% / 3.1% |

In proportion to the size itself the sample takes every large subtree and costs more for the same
precision (at depth 3, 297 CPU h for what sqrt gives at 225). The plan in `tho40-tune.spec`: 50 at
depth 3 and 100 at depth 4 by sqrt, about 150 CPU h an arm, the ratio to about +-4-8% at depth 3 and
+-2-3% at depth 4. The open nodes at depths 3 and 4 hold about three quarters of tho40's rest; depth
2 (85 open nodes, a mean subtree of 12.6 CPU h) is left out.

**The arms** (`infra/phase1/recipes/`, each tho40c.par with one change): `base`, tho40c;
`xbar4`, `branch.xbar_lines = 4` for 2; `cross5`, `dnn.cross_first = 5` with
`DNN_ADMM_WARM_BETA=0.1:30:0.36:1` in its environment (bound_lab/admm, section 6: 537 to 377
iterations on one parent's 36 children; the variable is an experiment setting no params file holds,
and making it one would change every run's settings hash, so paused runs could not resume); `z03`,
`cut.zone = 0.03` for 0.06. One worker each, side by side. Read them with
`python3 experiments/estimates_nov19/paired.py infra/phase1/runs/tho40-tune-261119-{base,xbar4,cross5,z03}`
after `collect.sh` on each.

**The re-dive** (`tho40-est3.spec`): tho40-est2's paused frontier dived again with the estimator fixed
(the spill), est2's settings otherwise, so that its difference from est2 is the fix: a direct check of
the 196 c7a CPU-days the correction gives.


**The arms' result (Nov 20; `paired.py` on the four records).** 137 sampled nodes (48 at depth 3, 89 at
depth 4), each arm to the end, weighted by the source's open nodes each stands for:

| arm | CPU, against base (5-95%) | nodes |
|---|---|---|
| base (tho40c) | 5,808 CPU h | 2,576,459 |
| xbar4 | 1.04x (0.98-1.09) | 0.96x |
| cross5 | **1.61x** (1.57-1.67) | 1.01x |
| z03 | 1.05x (1.03-1.08) | 1.00x: the same trees |

- **z03 builds base's trees** node for node: below depth 2 every node its first solve does not fathom
  already runs a phase at 0.03 (as on sko49 at 0.018 and sko42 at 0.025). Its 1.05x is the run-to-run
  noise between two arms here.
- **xbar4 is neutral** (4% fewer nodes, CPU within the noise). **cross5 costs 1.6x** for the same nodes:
  bound_lab's 537 to 377 iterations on one parent's 36 children does not hold over a tree. tho40c stands.
- **The solved sample against est3**, open nodes by depth (base's weighted totals, bootstrap 5-95%, against
  est3's strata at the same depth): depth 3 5,534 CPU h (3,882-7,308) against 4,481; depth 4 275 (205-354)
  against 803. At depth 3 the dives' estimate is inside the band (nodes 2.45 against 2.58 million), at depth
  4 it is 2.9x the solved sample's, outside it. With est3's depth 2 (1,426 CPU h, 85 open nodes, dives only)
  and depths 5-6 (16), the rest of tho40 is about 7,250 CPU h from the solved sample against est3's 6,727.

## The estimator against a known tree: sko42 (Nov 21)

`sko42-est-261121b` (`infra/phase1/specs/sko42-est.spec`): sko42 under `sko42.par` estimated as tho40 was (a 0.3 h
solve, Knuth dives k1.5, 900 s a unit, 5 dives at least per branched open node, sampling by depth 2:1, 3:0.6,
4:0.2, deeper 0.03, the decision point off), against `sko42-261120b`, the same tree solved: **324 CPU h (5-95%:
273-404) and 101,840 nodes estimated, 539.3 CPU h and 128,431 nodes solved**. CPU 0.60, outside the band; nodes
0.79. The estimator cost 70 CPU h.

`est_vs_solve.py`: each of the 2,246 sampled open nodes' estimate against its true subtree in the solve
(the trees are the same node for node: the same params and solver), weighted by 1/pi:

| open nodes at depth | strata (dived) | estimate CPU h | true CPU h | ratio | estimate nodes | true nodes | ratio |
|---|---|---|---|---|---|---|---|
| 2 | 169 (118) | 106.4 | 252.6 | 0.42 | 32,667 | 59,982 | 0.54 |
| 3 | 1,297 (241) | 143.4 | 227.3 | 0.63 | 49,553 | 58,570 | 0.85 |
| 4 | 750 (22) | 23.4 | 29.4 | 0.80 | 9,703 | 9,590 | 1.01 |
| 5-6 | 26 (1) | 4.3 | 4.9 | 0.88 | 2,104 | 2,100 | 1.00 |
| all | | 277.5 | 514.1 | 0.54 | 94,027 | 130,242 | 0.72 |

- **The sample of open nodes is not the problem.** Their true subtrees, weighted, come to 514 CPU h against the
  solve's 493 after the estimate's 46.7 done (+4%).
- **A node bounded by the estimator costs 0.79 of the same node in the solve** (their own bound, 28.7 against
  36.2 CPU h; 0.85 for DNN-fathomed, 0.79 cut-phase-fathomed, 0.78 branched), although bound_one is the solve's
  own bnb() and the twelve depth-1 nodes, solved in both, cost 2.05 against 2.09 CPU h. Whether the work or its
  timing differs, the DNN iterations would say (the two runs' dnn.txt.gz, in their exports). This is the bias
  seen on tai50b (0.80) and on tho40's depth-3 open nodes (est3 0.81 of the solved sample).
- **The dives miss the large subtrees.** Where the node counts are right (depths 4-6) the CPU is that 0.8-0.9;
  from depth-3 open nodes the dives find 0.85 of the nodes, from depth-2 open nodes 0.54, and no dive went below
  depth 6 although the tree reaches 9 and more. With a few dives an open node, Knuth's estimate of a large
  subtree is usually low and now and then very high: the mean is right, the typical draw is not.
- **What estimates well: a solved sample.** The open nodes sized by the dives, then a sample of them drawn in
  proportion to the square root of that size and solved to the end (`coord/sample.py`, as the tuning arms did):
  its weighted total is unbiased by construction. On tho40 the tuning arms' base gave 5,808 CPU h for the
  depth-3 and depth-4 open nodes (bootstrap 5-95% at depth 3: 3,882-7,308), where est3's dives gave 5,284.

**The solved sample, simulated on the known tree (`sample_sim.py`, Nov 21).** The strata of sko42-est-261121b's
dives, drawn as `coord/sample.py` draws (pi ~ the dive estimate^0.5, mix 0.2), each sampled node priced at its
true subtree CPU in sko42-261120b, 2,000 draws a design, against the true rest of 493 CPU h:

| sample at depths 2 / 3 / 4 | total, mean | 5-95% | the sample's own CPU |
|---|---|---|---|
| 20 / 50 / 50 | 513 | 356-697 | 92 h |
| 30 / 60 / 60 | 513 | 380-653 | 131 h |
| 40 / 80 / 80 | 511 | 413-626 | 173 h |
| 60 / 100 / 100 | 513 | 436-601 | 233 h |

It centres on the strata's own weighted truth (514; the 4% over 493 is which open nodes the dives run sampled,
the first stage), and is within +-20% for about 170 CPU h, where the dives alone gave 0.56 of the truth. A
sampled run on sko42 would add one draw of this; it was not run.

**The cost per node: what it is not** (Nov 21, on this container, two sko42 subtrees of 116 and 40 nodes cold
from `-F`, side by side): not the node and DNN traces (1,143 against 1,134 CPU s; 705 against 706), not the
solve units' checkpoints and warm-state shipping (`-o -k 60`, `QAP_DNN_WARM_DIR`: 1,112 against 1,124; 675
against 675), and not the solve's cutoff (it found 15,812 so late that every stratum was bounded at 15,812.05 in
both). What is left is the work itself (the estimator's warm start for a stratum) or the machines under load
(64 processes a worker). `QAP_EST_DNN_TRACE=1` keeps the DNN trace in a dives run's units, so that the
iterations of the same nodes can be compared: `infra/phase1/specs/sko42-est2.spec`.
