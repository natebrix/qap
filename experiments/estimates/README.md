# Tree estimates for tho40, wil50, sko42, sko49 and tai50b (Nov 10)

One place that says, for each instance, which settings are current and how to
get an estimate of the whole run under them. Before this, only wil50 had a
settings file of its own (`../wil50_feasibility/fullrun/fullrun.par`); the
other four estimates were made under settings that have since changed, and
sko42 had none.

## Current settings

Each instance has a directory here with `INSTANCE.par` (the rule blocks of
`latest/PARAM/nug_deep_xbar9_all.par` plus a `[settings]` section with every
solver setting), `settings.sh` (instance, symmetry file, incumbent; source it,
then `fullrun_check`) and `settings.golden` (the resolved settings from
HEAD's `bnbdrv`, checked before every run). wil50 keeps the full run's files.

| instance | settings | what they are | not measured on this instance |
|---|---|---|---|
| wil50 | `../wil50_feasibility/fullrun/` | the full run: CP cut phase, `cut.zone = 0.0075`, rule 9, DNN at depths 1-4 then 20 levels below a miss, give-up 300; `wil50t.dat` with its symmetry | |
| tho40 | `tho40/` | the Nov 8 nug30 record recipe (rule 9, DNN at depths 1-7, give-up 200, no cut phase) with `branch.xbar_lines = 2`, the tho40 A/B's setting (`../xbar_lines_sweep/tho40_ab/`); `tho40t.dat` with its symmetry | 2 lines was measured on three depth-3 subtrees only |
| sko42 | `sko42/` | the nug30 record recipe as it stands (rule 9, 1 line); `sko42t.dat` with its symmetry; incumbent 15,812, the known optimum | everything: the one run (Nov 8) used rule 8 |
| sko49 | `sko49/` | wil50's full run with the grid 7 x 7; `sko49t.dat` (new in `qaplib/`, the matrices of `sko49.dat` swapped) with the 7 x 7 grid's symmetry (`latest/SYM/sko49t.sym`, 10 depth-1 orbits), so `cut.side = 1` | the CP phase and the 0.0075 zone (wil50's Nov 5 levers were 2.2x and 4.2x on two sko49 subtrees, `../sko49_feasibility/vs_wil50.md`) |
| tai50b | `tai50b/` | the Nov 8 tuned estimate's values as keys (DNN at depths 1-9, give-up 600, the nug30 recipe otherwise) with rule 9, which its notes name; no symmetry, no cut phase (no grid); the dive sampler `est.xbar = 1`, `est.xbar_mix = 0.3` | rule 9: the Nov 8 runs used rule 8 (their log header), not rule 9 as the notes say |

A setting changes only in a commit that means to change it, with the golden
regenerated: `(. experiments/estimates/INSTANCE/settings.sh; fullrun_dump latest/bnbdrv > "$FULLRUN_GOLDEN")`.

## The runs (laptop, from the repository root)

`git pull origin cognac && make solver` first. Every script checks its
settings against the golden, uses all performance cores, keeps the Mac awake,
resumes when run again, and drains on ^C. Each prints the `git add` line for
its results at the end. Times are for the M2 Pro with 6 performance cores and
are estimates [e].

### wil50: three more measured subtrees (`wil50/`)

```
experiments/estimates/wil50/run_anchors.sh          # a220, a251, a325 in turn, then the refit
experiments/estimates/wil50/run_anchors.sh a325     # or one of them
```

- **Why.** The record's estimate reads each frontier node's subtree size off
  three measured subtrees by shortfall. The laptop's subtree of `13:1,21:0,45:2`
  (Nov 10) sits 4.7x above that curve. Fitting all four measured subtrees
  gives 8.4e5 DNN nodes (4.4e5-1.3e6) and, with their scatter (a factor 2.3
  either way), a mean of 1.2e6 (6.4e5-1.8e6), against the record's 4.9e5
  (`wil50/refit.py` with no arguments). More measured subtrees pin that curve.
- **The nodes**, from stage B's sample (wil50t indices, converged shortfall):
  `a220` 21:22,31:23,40:15 (220.2), `a251` 8:9,21:11,45:40 (251.1),
  `a325` 6:49,21:0,45:1 (324.6). The 210-300 class is 30% of the record's
  total and the tail (300 and beyond) 44%.
- **How.** `subtree.sh` runs one node's whole subtree through the coordinator
  with `fullrun.par` itself (not `run_pair.sh`'s arm settings), cold from
  `-F`, like every other measured subtree.
- **Time.** By the fit, `a220` and `a251` are about 20 and 40 minutes; `a325`
  is 2-3 h by the fit and about 9 h by the per-level model (at the laptop's
  12.7 s per DNN node of the Nov 10 subtree).
- **Result.** `wil50/subtrees/NAME/summary.txt` and `wil50/refit.txt`: the
  estimate with every measured subtree. The cold start probably overstates a
  frontier node's subtree (a frontier node starts warm, with its parent's
  pool), as for the record's anchors.

### tho40 and sko42: traced search, depth-2 sweep, fit (`trace_estimate.sh`)

```
experiments/estimates/trace_estimate.sh tho40 5
experiments/estimates/trace_estimate.sh sko42 6
```

- **The method** is tho40's of Nov 8 (`../tho40_estimate/README.md`):
  1. the production search under `parallel.py`, traced, for HOURS of wall;
  2. the depth-2 sweep: the root's children from a depth-0 pass, then each
     depth-1 node with its children (`bnbdrv -F node -d 2`);
  3. `estimate.py` fits log10(subtree size) to the gap per depth on the
     finished subtrees and applies it to the open nodes. Three calibrations:
     the instance's own finished subtrees and the solved tho30 and nug30
     trees.
- **Time.** tho40 about 5.5 h. sko42 about 6.5 h: its nodes cost 3.4 s on
  the laptop against tho40's 1.4.
- **Result.** `INSTANCE/trace/estimate.txt`, with the traces next to it (the
  checkpoint and its `warm/` store are ignored by git). To search longer
  after it has stopped: delete `INSTANCE/trace/search.done` and run again
  with more hours; the search resumes from its checkpoint, and the sweep and
  the fit are redone if their outputs are deleted too.
- **Checked here.** Fed the Nov 8 traces, step 3 gives exactly the record's
  three figures (6.1e5, 1.7e6 and 2.9e6 nodes). Steps 1 and 2 ran on tho40
  on the cloud box (`EST_OUT` and `EST_SWEEP_LIMIT` keep a smoke test out of
  the repo):
  - a two-minute search stopped at the wall limit with its checkpoint;
  - the depth-0 pass listed the 12 depth-1 nodes, and the root line is still
    facility 26 at 2 lines;
  - the first depth-1 node's sweep bounded 33 of its 39 children before the
    smoke test's own 25-minute limit stopped it (one process on the cloud
    box, which shared its cores with two other smoke tests). The Nov 8 sweep
    took 17 minutes of wall for all 12 depth-1 nodes on the laptop, 8 at a
    time.
- **sko42 doubles as a check.** Its optimum is known and Fujii et al. solved
  it with 6.0e6 nodes (rule M). If the estimate puts it in reach, solving it
  tests the estimate and the solver.

### sko49: the wil50 sampling plan (`sko49/run_stages.sh`)

```
experiments/estimates/sko49/run_stages.sh a          # stage A, about 3 h
experiments/estimates/sko49/run_stages.sh b          # stage B, about 3 h
experiments/estimates/sko49/run_stages.sh anchors    # 3 complete subtrees, hours each, then the estimate
```

- **Why this method.** sko49's estimate of Nov 6 (6.4e8 laptop CPU s) scaled
  a random sample of 12 depth-3 nodes, all of them clear misses. wil50's
  random-line sample said the same kind of thing (85-91% of the cost in the
  tail), and the production frontier then showed 74% of depth-3 nodes
  fathoming in their first solve and a tail of 1%.
  - **Stage A** bounds every node at depths 0-2 under the production search
    and saves the depth-3 frontier.
  - **Stage B** takes a seeded sample of 400 frontier nodes to their give-up
    point and converges the 80 farthest.
  - **The anchors** are three complete subtrees at the 50th, 80th and 95th
    percentile of the converged misses, the size curve for the refit.
- **Shared code changed.** `../wil50_feasibility/sampling/common.py` reads
  `SAMPLING_SETTINGS` and `SAMPLING_PARAMS` (unset: wil50's full run, as
  before). It also writes one-node state files from HEAD's version-9 state
  files (the "state version 9, expected 6-8" failure of Nov 9), and pins CPUs
  only where `taskset` exists, so the plan runs on macOS.
- **Checked here.** Stage A's root and one depth-1 task ran under sko49's
  settings: depth-1 node 41:0 is 450.8 short after 300 iterations (132 s on
  the cloud box) and has 27 children.
- **Disk.** Stage A keeps the DNN state of the parents of stage B's
  candidates: about 7 GB on wil50, 10-15 GB here.

### tai50b: Knuth dives (`tai50b/run_dives.sh`)

```
experiments/estimates/tai50b/run_dives.sh 8          # one estimator per core, 8 h
experiments/estimates/tai50b/run_dives.sh 8 6 201    # a second night, new seeds; the pool takes both
```

- **The method** is the Nov 8 tuned run's (`../tai50b_feasibility/NOTES.md`
  section 5): the estimator's dives, sampled by the DNN marginals, from level
  1, one seed per process, pooled by `../tai50b_feasibility/estimate/pool.py`
  into `tai50b/dives/pool.txt`.
- **What changes:**
  - the settings (rule 9, HEAD's QPB bound, which matters because QPB
    bounds every level below 9);
  - the number of dives. That run had 10, one of them 77% of the total.
    Six processes for 8 h should give on the order of 100.
- **An exact top of the tree** (the record's "dives from an exact depth-2
  frontier"): `start_dive_level` (the first number in `estimate.par`) at 2
  bounds the root's 50 children before the first dive, about 5 h in one
  process. At 3 it bounds the 2,450 depth-2 nodes, about 10 CPU-days in one
  process. Neither is set here: more dives buy more than an exact level 1.
- **Checked here.** Four minutes of one estimator on the cloud box. It reads
  the settings (rule 9, depths 1-9, give-up 600, the sampler keys) and starts
  its dives over the root's 50 children. No dive finished in that time; they
  took 10-25 minutes on Nov 8. HEAD's root QPB bound is 2.8858e8 against
  Nov 8's 2.8944e8: the Nov 8 change to QPB's parametric update moves the
  bound either way by instance.
- **Cloud box windows** (seeds 301 and up, four estimators for 1.9 h each;
  CPU-s are the cloud box's, with all four cores busy). The pool takes the
  laptop's seeds too: their nodes pool as they are, their CPU-s are in the
  laptop's units.

  | window | dives | nodes | CPU-s | largest dive | without it |
  |---|---|---|---|---|---|
  | 301-304 | 21 | 7.6e9 ± 7.2e9 | 3.1e10 ± 2.7e10 (360 days on 1,000 cores) | 86% of the CPU (seed 303 dive 3, to level 10, 1.5e11 nodes) | 4.1e8 nodes, 4.5e9 CPU-s (52 days) |
  | 305-308 | 19 | 8.2e9 ± 7.8e9 | 1.0e11 ± 9.9e10 (1,200 days) | 96% (seed 308 dive 3, to level 7, 1.5e11 nodes, 1.9e12 CPU-s) | |
  | 309-312 | 24 | 2.6e9 ± 2.3e9 | 3.0e11 ± 2.9e11 (3,450 days) | 99% (seed 312 dive 5, to level 9, 5.6e10 nodes, 7.1e12 CPU-s) | |
  | 313-316 | 7 | 3.1e7 | 3.8e8 | | stopped after 33 minutes by a container restart; the finished dives are kept |
  | 317-320 | 23 | 1.9e10 ± 1.8e10 | 1.3e11 ± 1.2e11 (1,520 days) | 90% (seed 317 dive 4, to level 9, 4.0e11 nodes, 2.7e12 CPU-s) | |
  | **pooled** | **94** | **8.7e9 ± 4.8e9** | **1.4e11 ± 8.2e10 (1,575 days)** | 312/5 55%, 317/4 21%, 308/3 15%, 303/3 4% | without the top three: 2.3e9 nodes, 1.3e10 ± 6.6e9 CPU-s (144 days) |

  The median dive is 4.5e6 nodes and 1.4e8 CPU-s. Level 10, which QPB
  bounds alone, costs 1.8 CPU-s a node; level 9, under the DNN, 108-125.
  Every full window has had one dive 3-4 orders of magnitude above the
  median (5.6e11, 1.9e12, 7.1e12, 2.7e12 CPU-s); seed 312's spends 97% of
  its CPU on 5.5e10 level-9 nodes, all fathomed by the DNN.

  **The dives cannot settle the total.** The tail of the per-dive CPU
  estimates has a Hill index of 0.35-0.42 (the largest 5 to 30 of the 94).
  Below 1, the sample mean is carried by the largest dive and keeps rising
  as dives are added, until the sample reaches the part of the tree where
  its finite size cuts the tail off; nothing in 94 dives shows where that
  is. The mean of such a sample is, more often than not, below the true
  mean. So what the dives establish is a floor: **under these settings
  tai50b's tree is at least about 1e11 CPU-s on this Xeon, 3-4 years on
  1,000 cores, and probably more.** More windows of the same dives sharpen
  that floor and the tail index; they do not converge.

  This also revises Nov 8. The tuned settings' 4.1e9 CPU-s against the base
  settings' 1.8e11 (45x) came from 10 dives with no giant among them; with
  94 dives the tuned settings give 1.4e11, the same order as the base.
  Whether the tuning helps is not settled by these dives either.

  The estimator samples each child with probability (1 - 0.3) x_ij / sum +
  0.3 / kids (`latest/bnb.C`, the dive branch) and divides the dive's
  weight by it, so the estimate stays unbiased; a dive that takes a few
  low-marginal children gets up to 3.3 kids times the weight per level.
  The log does not print the per-level probabilities, so whether the giant
  dives are such paths (the marginals misleading the sampler) or high-
  marginal paths into a large subtree is not known from these files.

## Suggested order

wil50's anchors first: they move the instance we are going to run. Then
tai50b overnight, tho40 and sko42 as one evening each, and sko49's stages A
and B as one night.

## Cold against warm: a220 as production would bound it (cloud box, Nov 10)

Every measured subtree behind wil50's estimate started cold from `-F`: no
parent DNN state and no inherited pool. Under `cut.zone = 0.0075`, the
parents of three of the four nodes measured so far (a220, a251, and the
laptop's 13:1,21:0,45:2) are 299-361 short and now run a cut phase, so in the
production search those nodes start warm and pooled.

a220 run as production would, `wil50/warm/` (`run_warm_a220.sh`):
1. its ancestors bounded by stage A's own tasks under `fullrun.par`
   (`stage_a.py --only r,r.12,r.12.22`);
2. the parent 21:22,31:23: first solve 299.2 short, then its cut phase to
   190.6 short (800 iterations), 28 children, a220 among them;
3. a220's one-node state from the parent's state, then its whole subtree with
   `bench/parallel.py --resume` on 4 cores. `--allow-settings-change` is
   needed because stage A records its copy of `fullrun.par`: the settings
   are identical value for value, and only the params file's hash differs.

| a220 | first solve | DNN nodes | iterations | DNN CPU |
|---|---|---|---|---|
| warm, pooled (this run) | 136.4 short, phase to 110.4 | 232 | 18,029 | 5,871 s (cloud box) |
| warm, no pool (stage B, zone 0.005) | 226.8 short | | | |
| cold `-F` (laptop, `wil50/subtrees/a220/`) | | | | |

The cold row comes from the laptop's `run_anchors.sh`.
