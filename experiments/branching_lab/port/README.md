# DNN strong branching in cognac: the port

Nov 22 2026. [m] = measured here, [e] = estimate or inference.

**The rule of `../root` is in a commit on top of cognac 21dd774a (`0001-dnn-strong-branching.patch`), as `sb.*`
settings in the params file, off by default. With it off the tree is node for node what it was; with it on, in-line
through the coordinator, nug30 closes in 2,980 nodes against production's 4,285 to 4,365 (0.69) and tho30 in 9,918
against 14,815 (0.67), both inside the ranges written beforehand, and its cost is in the run's own statistics.**

## What went in, and what did not

| in | the setting |
|---|---|
| the root row from every row's children | `sb.depths` from 0, `sb.all = 0` (needs `dnn.depths` from 0) |
| the shortlist at depth 1 (a row's real children only) | `sb.depths`, `sb.k = 4`, `sb.iters = 150` |
| depth 2 only where the node is far from fathoming | `sb.gap = 2:G`, the relative gap of the node's own DNN bound |
| the score | `sb.power = 0:8,1:1`: sum of shortfall^8 at the root, the sum of shortfalls below |

Not in: the tail rule, the line-score variants (`QAP_LS*`), the marginals dump, forcing a line
(`QAP_FORCE_LINE`), the exponential score. They stay in the experiment build of `../root`.

Terms: "strong branching" is bounding one level down and choosing the line from those bounds (the thesis's
term). `../root` and `../lsla` call the same thing "the lookahead" and "root scoring".

## The settings (RUNNING.md, "DNN strong branching (sb.*)")

```
dnn.depths = 0:28        # the root needs its own DNN state
sb.depths = 0:2          # where it may choose the line
sb.all = 0               # through this depth every free line is a candidate (rows only under symmetry)
sb.k = 4                 # below: this many rows, and columns where no symmetry is left, by marginal score
sb.iters = 150           # warm ADMM iterations per child
sb.power = 0:8,1:1       # from depth D on, the exponent P of the score sum of shortfall^P
sb.gap = 2:0.0582        # from depth D on, only where (inc - the node's DNN bound) / (inc - root QPB) >= G
```

- **The gap is the params blocks' quantity** ("use this rule when gap is worse than"), with the node's DNN bound in
  place of QPB's and tested after the node is bounded. nug30's 0.0582 is the 45.6 of `../root` over the 783.7
  between the incumbent and the root's QPB bound; tho30's 0.084 the 2,184 over 26,016.
- **The score needs no calibration** [m]: with the 8th power the root's first row is the best one run at nug30
  (r23), production's at tho30 (r5), and at tho40 r2 for every power from 2 to 12 (`../root/README.md`).

## The cost is the run's own

A child solve is not a node, so node counts do not show it. The RESULTS head and a table do, kept in state files
(version 10) and summed by `bench/parallel.py` and the coordinator:

```
 strong branch: 17 nodes, 3226 child solves, 212333 iterations, 6352.28 s

Strong branching by level:
level    nodes	  solves	      iters	    time
 0           1	     270	      40063	 1441.38
 1           9	    1500	     109995	 3358.32
 2           7	    1456	      62275	 1552.58
```

## nug30 (`run.sh nug30 sb nug30_sb.par`; 4 slots, the nug30 recipe with `dnn.depths` from 0, incumbent 6126)

| | nodes | CPU s: the tree | strong branching | total |
|---|---|---|---|---|
| production (`../root/runs/nug30_prod`, `_prod_tr`) | 4,285 and 4,365 | 3,314 and 3,421 | | 3,314 and 3,421 |
| the experiment build: r23 forced, the sum at depth 1 and gated depth 2 (`../root/runs/nug30_r23_d2_gate`) | 3,020 | 2,281 | 5,001 + 1,726 for the root's rows | 9,008 |
| the port, 8th power at every depth (`runs/nug30_sb_p8`) | 3,374 | 2,651 | 6,466 | 9,117 |
| **the port, its defaults (`runs/nug30_sb`)** | **2,980** | **2,265** | **6,352** | **8,617** |

- **The root chooses r23 in-line** [m]: 270 child solves, 1,430 to 1,441 CPU s.
- **The first run missed the prediction, and that fixed the default** [m]: with the 8th power at depths 1 and 2 the
  three largest depth-1 nodes branched on r4 where the sum had chosen c24, and two of those subtrees were 15 and
  18% larger (870 and 924 trace lines against 674 and 674). Below the root many children fathom or nearly and
  their number counts; the 8th power sees only the worst. So the exponent is by depth, and the default is what
  was measured in `../root`: the 8th power at the root, the sum below.
- **With the default the port reproduces the experiment build** [m]: the same line at eight of the nine depth-1
  nodes (c5 for c24 at 23:8, 58 trace lines either way), 3,010 trace lines against 3,033. Predicted 2,900 to 3,250
  (`predictions.txt`).
- **On a tree this small the rule costs more than it saves** [m]: 8,617 CPU s in all against 3,400. The tree itself
  is a third smaller. `../budget` is the case for the open grids: the same firings are about 1% of a run there.

## tho30 (`run.sh tho30 sb tho30_sb.par`)

4 slots, the nug30 recipe with `dnn.depths` from 0, incumbent 149938.

| | nodes | CPU s: the tree | strong branching | total |
|---|---|---|---|---|
| production (`../grid/runs/tho30_prod`) | 14,815 | 9,375 | | 9,375 |
| the experiment build: the sum at depth 1 and gated depth 2 (`../root/runs/tho30_d2_gate`) | 9,713 | 7,041 | 7,772 | 14,813 |
| **the port, its defaults (`runs/tho30_sb`)** | **9,918** | **7,870** | **8,876** | **16,746** |

```
 strong branch: 17 nodes, 2736 child solves, 279927 iterations, 8876.30 s
 0           1	     300	      44901	 1664.98
 1          10	    1460	     158904	 5160.15
 2           6	     976	      76122	 2051.18
```

- **The root chooses r5 in-line** [m], production's row: 300 child solves, 1,665 CPU s. Predicted: r5, 9,200 to
  10,300 nodes.
- **The same lines as the experiment build at nine of the ten depth-1 nodes** [m] (c29 for c9 at 5:0, 811 trace
  lines against 808), and at each of the six depth-2 nodes that fired.
- **Six depth-2 nodes fired, not nine** [m]: the gate is on the node's own DNN bound, and with the root's DNN state
  the depth-1 and depth-2 solves start warm and end a little higher. Three nodes that fired in the experiment build
  sit at 0.0803 to 0.0809 here, under the 0.084. That is most of the 2% between the two trees (10,018 trace lines
  against 9,778) [e]: the gate is a threshold read off one run, and a node near it can fall either side.
- **Here too the rule costs more than it saves** [m]: the tree's CPU is 0.84 of production's, the total 1.79.
  n = 30 is the wrong size to pay for 2,736 child solves; see `../budget`.

## Off is off

- nug20 and nug22 with the nug30 recipe, one process, node traces of cognac 21dd774a and of the port: identical but
  for the CPU column (261 and 202 lines) [m].
- cognac's suites on the port: brute force, resume, parallel, coordinator, settings table, regression check
  (7 instances, tolerance 0), node trace, ship: pass. The estimator test passes or fails here by its CPU check
  alone (the estimated CPU within 15% of the run's): 1.05 to 1.21 over six runs on the port, 1.14 to 1.17 over
  three on cognac 21dd774a (two of the three fail), each while an acceptance run had the machine. A timing check
  at its threshold on this box, not the port's; its node estimate is the same in every run (1.031) [m].
- `make sb` (`tests/sb_test.py`, new): the statistics against the trace, the progress file, `bench/parallel.py`'s
  totals against serial, the gate, refused settings.

## Review

An independent read of the diff (a subagent, the diff and no conclusions of mine) found one defect, fixed: a child
solve that ends without a certificate (-inf, or tai35b's -3e11 under a pool) made the score NaN or flattened every
other line's. A child's bound is now never taken below QPB's for it. Also fixed from it: child solves are kept in
memory for the process (the estimator's dives pass the same top nodes many times), the progress file's key carries
the settings hash and the root bound, long trace lines. Left as it is and documented: under `pool.cp = 1` a pooled
child is bounded in the pool form, not the CP form of its own first solve.

## The review on the PR (Nov 22, `0002-review-fixes.patch`, cd8d5b4e on the PR's branch)

Four comments, all taken:

- **The settings hash** covered the six new keys at their defaults, so every paused run would have needed
  `--allow-settings-change`. Now the `sb.*` keys are hashed only when `sb.depths` is set. The five golden hashes
  are cognac's again; cognac and the port resume each other's nug20 checkpoint (217 nodes either way) [m].
- **The state version** was 10 for every file. Now 10 only once strong branching chose a line, else 9, the file
  cognac writes [m: `tests/sb_test.py`].
- **A failed child solve** scored at QPB's bound. Now a child's bound is never below the node's own DNN bound, and
  with `dnn.retry = 1` a failed child is solved again without the pool, then cold. The trace counts the retried
  and the failed. Under the first test (1% below the inherited bound) tai15b's root at 40 iterations had 194 of
  222 children fail, 16 improve cold, the same tree [m]; the test changed in the second round, below.
- **The marginals** for the shortlist were written into `dnn_xbar`, which `est.xbar` reads. Now their own vector.

A second round (`0003-child-failure-test.patch`, 6d1dfb14): the reviewer proposed that a child solve fail only
below QPB's bound for the child, since a capped solve that is merely short of converged ends above it. On
tai15b's root 52 to 196 of 222 children fail the 1% test at 150 to 1 iterations and none ends below QPB's bound
[m], so that is the test now, and the floor under a child's bound includes `node.bound`. With it no test case has
a failed child: the cold retry was seen running only under the 1% test, the retry without the pool never.

**Merged into cognac Nov 22 (7c468e25, PR 6), the three commits.**

With strong branching on, nug20 and nug22 give the same node traces and strong-branching traces before and after
[m]. The nug30 and tho30 runs above are on the first commit and were not repeated: in them no child's bound was
below its node's, so neither the floor nor a retry would have fired.

## What to know before merging

1. **A run without `sb.depths` keeps its settings hash and its state version** (after the review's fixes). The
   five `settings.golden` of `experiments/estimates` change in their listing only.
2. **A state file is version 10 once strong branching ran.** The website's parser reads version 10 from
   brixius-web 39764ee on (branch `claude/cool-babbage-ilnogj`): deploy it before the first run that uses strong
   branching publishes its totals, or the live pages refuse the file.
3. **A strong-branching node interrupted by a lost instance starts again.** `QAP_SB_PROGRESS=file` avoids that
   where the file outlives the process; the coordinator does not manage it (the root's pass is in the
   coordinator's root cache as before). At tho40 a depth-1 node is about 2.6 core-hours.
4. **The first nodes are one process each.** The root's rows, then each depth-1 node, run on one core while the
   other slots wait for work: at tho40 about 4 core-hours for the root, then 12 depth-1 nodes side by side.
5. **The estimator.** Dives pass through `bnb()`, so a dive that starts above `sb.depths`' last depth does the
   same strong branching (once per process and node). An estimate whose strata start below pays nothing.

## To apply

```
git checkout cognac && git am experiments/branching_lab/port/0001-dnn-strong-branching.patch
make && make test
```

The patch is against 21dd774a. Nothing here was pushed to cognac. The same commit, rebased onto cognac 6f614567
(three commits of run records later, no file in common), is branch `claude/dnn-strong-branching`:
https://github.com/natebrix/qap/pull/6.

## Files

`0001-dnn-strong-branching.patch` (the commit), `0002-review-fixes.patch` and `0003-child-failure-test.patch`
(on top of it; cognac has all three since 7c468e25); `run.sh`, `nug30_sb.par`, `tho30_sb.par`, `nug30_sb_p8.par` (the
acceptance runs); `predictions.txt` (written before each run); `runs/` (results, the node trace, `sb.txt`: a line
per strong-branching node, `sb_progress.txt`: its child solves).
