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
| sko42 | `sko42/` | the nug30 record recipe as it stands (rule 9, 1 line); `sko42t.dat` with its symmetry; incumbent 15,812, the known optimum. Measured worse than wil50's settings with the 6 x 7 grid, `sko42/tune/B.par` (below) | the recipe: the one run (Nov 8) used rule 8 |
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
  pool), as for the record's anchors. For a220 it does not: see *Cold
  against warm* below.
- **Results (laptop, Nov 10-11).** All three ran to the end with
  `fullrun.par` on the 6 performance cores:

  | node | shortfall | DNN nodes | DNN iterations | DNN CPU | per node | wall |
  |---|---|---|---|---|---|---|
  | a220 | 220.2 | 232 | 18,852 | 3,071 s | 13.2 s | 15 min |
  | a251 | 251.1 | 3,935 | 387,735 | 62,288 s | 15.8 s | 3.6 h |
  | a325 | 324.6 | 11,923 | 1,066,452 | 159,914 s | 13.4 s | about 8 h running (17.7 h end to end) |

  a251 took 3.6 h against the predicted 40 minutes. a325 ran about 8 h
  (6.0 h until the evening drain, then 2.1 h until the disk filled at 22:41,
  `No space left on device` in its `worker.log`), against 2-3 h by the fit
  and 9 h by the per-level model. The resume replayed the
  finished units, and the run went on to the end.

  The refit (`wil50/refit.txt`) uses seven measured subtrees, two of
  them ([e]) estimates themselves: log(size) = 1.742 + 0.0205 s, a factor
  2.9 of scatter either way. In DNN nodes for the whole of wil50, at 22.5
  CPU-s per node:

  | size model | DNN nodes (95%) | core-hours |
  |---|---|---|
  | the record's anchors | 4.9e5 (2.8e5-7.5e5) | 3,070 |
  | the fit | 1.2e6 (6.1e5-1.9e6) | 7,500 |
  | the fit with its scatter (mean) | 2.1e6 (1.1e6-3.3e6) | 13,200 |
  | per level, 47^(s/130) | 3.0e6 (1.1e6-5.5e6) | 18,600 |

  The four-subtree fit gave 8.4e5, with a mean of 1.2e6. With the
  anchors, wil50 is 2.4-4.3 times the record's estimate: about 1.2-2.1
  million DNN nodes and 7,500-13,200 core-hours.

- **Through AWS's queue and store.** `subtree.sh` takes `EST_BACKEND`, a
  coordinator backend other than its local directory: S3 for the state
  files, two SQS queues, two DynamoDB tables (`coord/README.md`). The solver
  still runs on the machine that starts the script. Re-running a220 this way
  checks the AWS backend on real AWS against a known tree (232 DNN nodes);
  until now it had only run against `moto`. It costs a few cents in requests.
  Checked here: a nug20 subtree (`-F 0:5`, `fullrun.par`, 38 DNN nodes) on
  the local backend and on the AWS backend against a local moto server gave
  the same trace line for line, apart from process IDs and CPU times, and a
  rerun of the finished AWS run took the resume path and fetched the same
  result.

  ```
  python3 -c 'import boto3' || python3 -m pip install --user boto3
  ACCT=$(aws sts get-caller-identity --query Account --output text)
  export EST_BACKEND="aws:qap-coord-$ACCT-eu-north-1?region=eu-north-1"
  python3 -m coord.qapctl aws-setup --backend "$EST_BACKEND"     # once: the bucket, queues and tables
  EST_RUNS=experiments/estimates/wil50/subtrees experiments/estimates/subtree.sh \
    experiments/wil50_feasibility/fullrun/settings.sh a220aws 21:22,31:23,40:15
  ```

  **On real AWS (laptop, Nov 11, `wil50/subtrees/a220aws/`): the same
  tree.** The same 232 nodes, each with the same status, iterations,
  cut phase and fathoming as the local run. Four nodes' first bounds differ
  by at most 0.15, all of them fathomed either way: the coordinator split the
  work into 47 units against 51, and where a unit boundary falls depends on
  timing. The solver binary is the same (no change in `latest/` between the
  two commits). The run took 1,013 s of wall time against 926 s, with the
  slots busy 51.5% of the time against 65.3%: on this 15-minute run, the
  queue and store round trips cost about a tenth of the wall time. No errors
  in either log.

  Removing it all: `aws s3 rb s3://qap-coord-$ACCT-eu-north-1 --force`;
  `aws sqs delete-queue --queue-url $(aws sqs get-queue-url --queue-name Q --query QueueUrl --output text)`
  for `Q` in `qap-units-large` and `qap-units-small`; `aws dynamodb
  delete-table --table-name T` for `qap-runs` and `qap-units`.

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
  | 321-324 | 21 | 3.2e9 ± 3.0e9 | 4.2e10 ± 3.2e10 (490 days) | 74% (seed 321 dive 5, to level 10, 6.2e10 nodes, 6.5e11 CPU-s) | 2.5e8 nodes, 1.2e10 CPU-s (134 days) |
  | **pooled** | **115** | **7.7e9 ± 4.0e9** | **1.2e11 ± 6.8e10 (1,375 days)** | 312/5 52%, 317/4 20%, 308/3 14%, 321/5 5%, 303/3 4% | without the top three: 2.5e9 nodes, 1.8e10 ± 8.0e9 CPU-s (209 days) |

  The median dive is 4.7e6 nodes and 1.5e8 CPU-s. Level 10, which QPB
  bounds alone, costs 1.8-7.6 CPU-s a node; level 9, under the DNN, 108-125.
  Every full window has had one dive 3-4 orders of magnitude above the
  median (5.6e11, 1.9e12, 7.1e12, 2.7e12, 6.5e11 CPU-s); seed 312's spends
  97% of its CPU on 5.5e10 level-9 nodes, all fathomed by the DNN.

  **The dives cannot settle the total.** The tail of the per-dive CPU
  estimates has a Hill index of 0.37-0.45 (the largest 5 to 30 of the 115).
  Below 1, the sample mean is carried by the largest dive and keeps rising
  as dives are added, until the sample reaches the part of the tree where
  its finite size cuts the tail off; nothing in 115 dives shows where that
  is. The mean of such a sample is, more often than not, below the true
  mean. So what the dives establish is a floor: **under these settings
  tai50b's tree is at least about 1e11 CPU-s on this Xeon, 3-4 years on
  1,000 cores, and probably more.** More windows of the same dives sharpen
  that floor and the tail index; they do not converge.

  This also revises Nov 8. The tuned settings' 4.1e9 CPU-s against the base
  settings' 1.8e11 (45x) came from 10 dives with no giant among them; with
  115 dives the tuned settings give 1.2e11, the same order as the base.
  Whether the tuning helps is not settled by these dives either.

  The estimator samples each child with probability (1 - 0.3) x_ij / sum +
  0.3 / kids (`latest/bnb.C`, the dive branch) and divides the dive's
  weight by it, so the estimate stays unbiased; a dive that takes a few
  low-marginal children gets up to 3.3 kids times the weight per level.
  The log does not print the per-level probabilities, so whether the giant
  dives are such paths (the marginals misleading the sampler) or high-
  marginal paths into a large subtree is not known from these files.

  **Stopped, Nov 11: tai50b is out of reach under the current method.**
  Six windows (115 dives, about 10 h of four cores) put a floor of about 1e11
  CPU-s on the tree, and the tail says the true total is probably higher;
  more dives of the same kind cannot lower that. Where the pooled CPU goes:

  | level | ≤5 | 6 | 7 | 8 | 9 | ≥10 |
  |---|---|---|---|---|---|---|
  | share of CPU | 0.7% | 4.6% | 13.3% | 11.4% | 64.6% | 5.4% |

  **Taking the DNN to every depth would not help.** Levels 1-9 are
  already under the DNN and hold 95% of the CPU; QPB bounds only levels 10
  and below, which hold 5.4%. The number of level-10 nodes is set by what
  the DNN fathoms at level 9, so a DNN below 9 can save at most that 5.4%,
  and less in practice, since a DNN node at level 9 costs 108-125 CPU-s
  against QPB's 1.8-7.6. The largest dives end at level 9 or 10 (312/5 and
  317/4 at 9, 321/5 and 303/3 at 10, 308/3 at 7). What remains open is bound
  strength at levels 6-9, the levels that hold the cost: DNN bounds run
  closer to convergence than a 600-iteration give-up, at a higher cost per
  node, which Fujii et al.'s 5.5e7 nodes suggests can shrink the tree. These
  dives do not test that, so it remains untested, not ruled out.

## Suggested order

wil50's anchors first: they move the instance we are going to run. Then
tai50b overnight, tho40 and sko42 as one evening each, and sko49's stages A
and B as one night.

## sko42: settings and a preliminary estimate (cloud box, Nov 11-12)

`sko42/tune/README.md` has the details. Three settings were tried on sko42:
- **A**, `sko42.par` as it stands: the nug30 recipe, no cut phase.
- **B**, wil50's full run with `cut.grid = 6:7`.
- **C**, B with `cut.zone = 0.010`.

**The comparison.** At depth 1, B's first solve is 87-175 tighter than A's.
On five matched complete subtrees, B takes 1.15-2.8x less CPU than A. On the
farthest of them, A was stopped with its tree still open, already at 1.9x
B's CPU. C beat B by 24% and 36% on the two nodes tried so far.

**The estimate under B.** The wil50 sampling plan:
- stage A: depths 0-2 bounded, 10,726 frontier nodes;
- stage B: 250 sampled frontier nodes, 49% fathomed in their first solve, the
  rest at most 189 short;
- ten complete subtrees, read by the sampled nodes' give-up shortfalls.

That gives about 1e5 DNN nodes (8.5e4-1.3e5) and 2.2-2.4e6 CPU s:
600-670 core-hours on the cloud box, about 300-335 on the laptop, some 2-2.5
days on its 6 performance cores. Fujii et al. needed 6.0e6 nodes. The
`trace_estimate.sh sko42` run above uses `sko42/settings.sh`, that is A.

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
| cold `-F` (laptop, `wil50/subtrees/a220/`) | 254.1 short, phase to 115.6 | 232 | 18,852 | 3,071 s (laptop) |

**Starting warm and pooled does not change a220's tree.** Both runs give
the same 232 DNN nodes by depth (1, 47, 184), with 227 fathomed and the same
7 cut phases. The warm start nearly halves the root's first shortfall (136
against 254). But its cut phase ends about where the cold one does (110
against 116 short), in 600 iterations against 980. Over the whole subtree,
the warm run takes 4% fewer iterations (18,029 against 18,852): 2% fewer in
first solves and 7% fewer in cut phases. CPU does not compare: the cloud box
spends about twice the laptop's time per iteration (0.33 s against 0.16 s).
So, at least for this node, the cold measurements behind the refit do not
overstate the subtree. a220 is one small subtree, though: whether the same
holds for a deep one such as a325's 11,923 nodes is not measured.
