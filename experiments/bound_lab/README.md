# Bound lab, Nov 15-16 2026: new cuts and ADMM changes for the open instances

Question (Nathan, Nov 15): are there new cuts that would significantly improve the bounds, or
significant improvements to the ADMM, that would materially change the estimates for the open
QAPLIB instances? Then: pull on the most promising threads, looking for one big swing on tho40 and
two on sko49 and wil50.

This page is the index and the summary. Each directory has its own README with the runs behind the
numbers. **[m]** measured (in the directory named), **[e]** estimated. Work was done by agents in
this session and in cloud sessions, each on its own branch `claude/cool-babbage-ilnogj-*`; their
reports are merged here. Status as of Nov 16, 09:30 UTC: two cloud sessions still running (tai50b
level by level, `tai50b_bfs` on `claude/cool-babbage-ilnogj-tai50b`; tai60b dives).

## Summary

| instance | before | now | what did it | how sure |
|---|---|---|---|---|
| **tai50b** | at least 1.1e6 core-hours, no cut phase | **about 100 core-hours** (60-250) [e] | cluster-region cuts | 275 dives; tai35b and tai40b solved end to end with the same code; levels 4-5 thinly sampled |
| **tho40** | 1,308 c7a CPU-days (tho40-b as running) | **226 CPU-days** with the cut phase at zone 0.06 [m, dives]; **about 110-215** with hand-over and the phase penalty [e] | the grid cut phase used widely, then a rule for when to stop a phase | 39-51 dives per arm over the whole tree; the last factor of two is paired on 34 dives, noisy in total |
| **nug30** | 3,030 CPU-s, 3,688 nodes (record recipe, cloud box) | **1,667 CPU-s, 1,302 nodes** (1.82x) [m] | long phases at depth 1, short below; phase penalty; early entry | whole search, repeated, path checks pass |
| **tai60b** | no estimate | root 605,315,696 certified (published 592,829,061; best known 608,215,054); tree **500-2,500 core-hours** (400 to more than 10,000) [e] | cluster-region cuts | a model from samples on the path to depth 18 and off it at depths 1-5; no dives yet |
| **tai35b, tai40b** | solved by others in 2.6 million and 278,000 nodes | **solved in 922 and 733 nodes**, about 1.4 CPU-hours each [m] | cluster-region cuts | two region families each |
| **wil50** | 7,500-18,600 core-hours | about 1.15-1.2x less [e]; 1.09x and 1.35x measured on two subtrees | early entry, hand-over rule refitted | two subtrees |
| **sko49** | 37,000-41,000 core-hours | about 1.2-1.35x less [m, two subtrees] | the same, plus `wc=1.4` | two subtrees |

- **Two swings, on tho40 and on the tai-b family.** No second swing was found for wil50 and sko49:
  the cuts are closed there (the Nov 4 study stands), and what transfers from tho40 is a policy worth
  1.1-1.35x.
- New best lower bounds at the root, all certified: tho40 230,452.99 (published 226,522); tai50b
  453,552,676 (431,245,612); tai60b 605,315,696 (592,829,061); tai80b 798,336,758 (786,900,529)
  after 1,000 phase iterations and still rising.

## What was found

### 1. Cluster-region cuts for the tai-b instances (`regions/`, `tai50b_size/`, `tai_solves/`, `tai60b/`)

The regional triangle inequalities of `latest/dnncut.h` are valid for any regions. Taillard's type-b
instances have locations in a few tight, well-separated clusters; with regions that separate the
clusters (prefixes along the principal MDS axis, unions of the top clusters, between-cluster
half-planes; `regions/gen_regions.py`) the cut phase closes most of the root gap: 96% on tai25b, 97%
on tai30b, 71-81% on tai35b to tai60b. Solver change: `cut.regions = FILE` (`regions/regions.patch`).

- **The tree becomes a path with a fringe.** Off-path children fathom in their first solve; what
  survives is the best known solution's path and nodes holding a solution nearly as good (swaps
  inside one cluster). tai35b and tai40b solve in 922 and 733 nodes from an incumbent one above the
  optimum, each by two region families, with no bound on the optimum's path above the optimum
  (`tai_solves/`).
- **tai50b**: 275 Knuth dives and a model of the open nodes put it at 5,000-13,000 nodes and about
  100 core-hours on the cloud box, about 1.4x that on c7a (`tai50b_size/`). An independent sample
  below one branching sequence (`regions/`, "Follow-up") read 2,000-7,000 core-hours if the fringe
  compounds for three levels; the dives say two. A level-by-level run is measuring it.
- **tai60b** (`tai60b/`): the root closes 81% of the published gap; the on-path bound is 0.45% short
  at depth 1 and 0.01% at depth 18. Off the path the fringe of same-cluster swaps compounds for at
  least two levels. A model of the surviving prefixes gives 500-2,500 core-hours; the slack at depths
  17-30 decides between 400 and 10,000. tai80b: the root only (8 s per iteration).
- **A fault to know about** (`tai_solves/`, section 5): with `wc=2.88` a child's first solve on its
  parent's pool sometimes fails outright (a certificate of -3e11, or none). The bound stays valid but
  the subtree blows up; `tai_solves/pool_retry.patch` redoes the solve without the pool. Not yet in
  the integration build.

### 2. tho40: the cut phase as a policy (`tho40_search/`, `tho40_nodes/`, `cpfast/`, `tho40_dives/`)

The grid cut phase (`cut.grid = 5:8`, side 1) was known; tho40-b runs without it. Three things make
it pay:

1. **A wide zone** (`cut.zone = 0.06`): a phase at every give-up within 14,400 of the cutoff. By
   Knuth dives over the whole tree: 226 c7a CPU-days (se 39) and 1.1 million nodes, against 1,308 and
   34 million [m, `tho40_dives/`].
2. **A hand-over rule** (`cut.handover`, `tho40_search/stop/`): stop a phase when finishing it costs
   more than handing its pool to the children; with early entry (`dnn.giveup_warm = 100`,
   `dnn.gu_zone = 0`).
3. **The phase's own ADMM penalty** (`dnn.beta_sched = w=0.36,wc=2.88`, `cpfast/`): the phase had
   been running at the warm children's penalty; far phases get 2.4-4.6x shorter.

With 2 and 3, the dives take 0.49 of zone 0.06's CPU on the same paths (34 pairs), which gives about
110 CPU-days; the combined arm's own total is 215 (se 82). On the A/B subtree `26:10,8:0` the four
depth-3 classes measured 7.8x under production (`tho40_search/`, "With the faster phase").

Recommended block, on `experiments/estimates/tho40/tho40.par` with `dnn.gu_zone = 0` and
`dnn.beta_sched = w=0.36,wc=2.88`:

```
dnn.giveup_warm = 100
pool.on = 1
cut.grid = 5:8
cut.side = 1
cut.cp = 100:1000
cut.warmy = 1
cut.rounds = 12
cut.iters = 200
cut.freeze = 1000
cut.zone = 0.06
cut.handover = 20:0:0.22:1.5:100:3:1.25
cut.cp_wc_gap = 300
```

The first ten lines need no patch and are the 5.8x. The last two need the integration build.

### 3. nug30: the same three levers on a whole search (`nug30_phase/`)

Zones alone lose on nug30 (a give-up is cheap to branch on there). What wins is the pool: nine long
phases at depth 1 whose cuts are handed down, phases below stopped at 100 iterations unless nearly
done (`cut.giveup = 100:5,300:2`, `cut.giveup_depth = 2`), the phase penalty (`wc=1.44`) and early
entry. 1,667 CPU-s against 3,030; on the laptop's 3.0-minute record that would be about 1 min 40 s
[e, not run there].

### 4. wil50 and sko49 (`wil50_policy/`)

The wide zone is already there. Early entry and a hand-over rule that is less patient than tho40's
(`cut.handover = 32:0:25` on wil50, `27:0:5` with `wc=1.4` on sko49) give 1.1-1.35x. `wc` hurts
wil50.

### 5. ADMM (`admm/`, `cpfast/`)

- The phase penalty above is the one large ADMM finding.
- Diagonal congruence scaling (`dnn.dscale = 0.5`, `dnn.dscale_until = 100`): 1.2-1.6x on tai-b,
  nothing on the grid instances.
- tho40 children: `dnn.cross_first = 5` with a lower warm penalty at the start takes 36 children
  from 537 to 377 iterations (about 12% of the tree's CPU [e]); not in the combined arm.

### 6. Tested and closed

- **Exchange-optimality cuts** (`swapcuts/`): inequalities valid only for optimal permutations
  (no pairwise exchange improves). Valid, and the exact ceiling is high, but the DNN realises a tenth
  of it: +1.4 to +2.7 on wil50, about +135 on tho40. By-product: wil50's best known solution has a
  zero-gain exchange, so there is a second permutation of value 48,816.
- **5-point, local-7 and box families on tho40** after the triangle phase: flat (`tho40_nodes/`).
- **Rule 2 at the bottom of tho40's tree** (`rule2_ab/`): the same trees as rule 9 throughout,
  0.5-1.8% more iterations.
- **Literature 2021-2026** (`literature/`): nothing published solves tho40, sko49, wil50 or tai50b;
  the regional-triangle idea was not found in print; no relaxation above the DNN is computed at
  n >= 30.

## The solver changes (`integration/`)

One tree with all the patches: commit **b86352e** on `claude/cool-babbage-ilnogj-integration`
(`integration/all.patch` against `cognac` e1596e22). Nine new keys, all default-off: `cut.handover`,
`cut.cp_wc_gap`, `cut.cp_wc_deficit`, `cut.cp_extrap`, `cut.regions`, `cut.giveup_depth`,
`dnn.dscale`, `dnn.dscale0`, `dnn.dscale_until`; plus a fix to `bnbdrv -e -r`.

- Builds clean; `make test`, `nodetrace`, `shipcheck`, `dnncheck` pass; `continuationcheck` fails
  identically before and after the merge.
- With no new key set, node and DNN traces are identical line for line on six runs.
- **The settings hash changes** (71 to 80 entries): a state file written before is refused without
  `--allow-settings-change`. A run in flight has to change binaries on all workers at once.
- 18 review findings, none a bug. Before a long run: `cut.regions` records the path, not the file's
  contents; `cut.handover` silently replaces `cut.stall`; two keys accept bad values silently.
- Not in it: `tai_solves/pool_retry.patch`, and the combination of `cut.giveup_depth` with
  `cut.handover`, which nobody has run.

## What is not measured

- The tho40 combined arm's absolute total (uncertain by about 2x either way); the patches in a run
  through the coordinator; a switch of tho40-b's settings mid-run.
- tai50b below depth 5 other than on the path; the failed pooled solves' cause.
- tai60b by dives; tai80b's tree at all.
- wil50 and sko49 whole-tree effect of the policy (two subtrees each).
- nug30 on the laptop.

## Housekeeping

Early commits on this branch carried large state files (five `*.state` of 44 MB, seven `*.warm`);
they were untracked in later commits but remain in the branch's history. Squash on merge.
