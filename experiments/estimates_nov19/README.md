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
- **CPU per node, an open question.** Where the counts agree, the estimator's CPU per node is under the
  solve's: depth 4, 212 s against 264 (0.80); depth 5, 53 s against 88 (0.60); depth 6, 38 against 68
  (0.60). The solve's own part of this run against tai50b-r9's: 0.95 at depth 1, 0.80 at 2, 0.91 at 3.
  A node bounded in setup or a dive warm-starts from its parent's state in memory (64-bit), a node of
  the solve often from a shipped one (`ship.bits`), after its unit was split; or the runs differ by
  more than this. Settled node by node with the dives (collect) and tai50b-r9's node trace (its export's
  `nodes.txt.gz`): the same 282 depth-4 nodes in both, by their fixes.

**What the next estimate of a tree like this should change:** Knuth dives (one child a level, about 88 s
a level, so many dives in the budget) with `QAP_EST_MIN_DIVES=5` and a longer dive budget; or no
estimate at all: tai50b's whole tree is 117 CPU h, and this estimate cost 70 (solve 32.8, estimator
37.5). For a small tree the solve is the estimate.
