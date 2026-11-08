# wil50 sampling plan: stages A and B (Nov 7, 2026)

Step 6 of `../roadmap/README.md` section 4 and the plan of
`../roadmap/tree_size.md` section 3. The full run's cost estimate spans 6x,
and 85-91% of it comes from a few far-from-cutoff depth-3 nodes, so the real
depth-3 frontier that the production branching rule builds at depths 0-2 is
the main unknown. Stage A builds that frontier; stage B samples it.

**[m]** = measured here, with its file; **[e]** = estimate.

## Settings

| part | setting | why |
|---|---|---|
| environment | `../fullrun/settings.sh`, sourced by the launchers (`common.fullrun_env`), with `fullrun_check` on the binary; the caller's own `QAP_*`/`DNN_*` variables are dropped | the full run's settings, decided Nov 7 (the `cpcutt` arm) |
| instance | `qaplib/wil50t.dat`, `-s latest/SYM/wil50t.sym`, `-u 48816` (cutoff 48,814.05) | `../sym/README.md` |
| parameter file | `latest/PARAM/nug_deep_xbar9_all.par` | the file of every measured production arm (`../pipeline_ab/run_pair.sh`): rule 9 through depth 6, rule 2 below |
| depth spec | `-x 1:4:1:300:1` | **D1 = 1**: the DNN bound from depth 1 (`tree_size.md` section 4), so that rule 9 has marginals at depths 1-2; the root gets QPB only (its line is a row restricted to the 15 J1 cells, which the marginals cannot change). **D2 = 4**: with `QAP_DNN_EXT=20` every child of a miss keeps the DNN bound, so D2 >= 3 is all that matters (section 4); 4 makes depths 3 and below run exactly the absolute depth spec of the measured node-10 arms (`-x 0:4:1:300:1` from a depth-3 node, `../pipeline_ab/`). **G = 1** (every gap), **H = 300** (the production give-up), **W = 1** (warm start from the parent) as in every production arm |

## Why one bnbdrv task per node, not the coordinator

`QAP_DNN_SHIP=32` puts a branched node's final DNN state into the state file
its children travel in. At depth 2 (m = 48) that is the upper triangles of Z
and Y, (48^2+1)^2 entries each at fp32, 28.4 MB of base64 per parent, and
about 85 MB of doubles in memory while its children are open [m: token size
on nug20, scaled]. A `bnbdrv -d 2` run, or a coordinator unit holding many
depth-2 nodes, keeps every finished depth-2 node's state until it writes the
frontier: 635 x 85 MB = 54 GB of memory and 635 x 28 MB = 18 GB of state
files, against this box's 15 GB of memory and 17 GB of free disk. The
coordinator also re-queues returned units' depth-3 nodes, so each would be
copied with its parent's state several times.

So `stage_a.py` runs each node as its own task, the way a fleet run hands a
node to another process: a state file holding the node alone with its
parent's warm start (QPB point, cut pool, DNN state), resumed with `-r` and
bounded with `-d <its depth>`, so its children stay open in the task's
`out.state`. One process per slot, pinned with `taskset`, one thread each.
The root runs first, then the 15 depth-1 nodes, then the depth-2 nodes as
their parents finish, round-robin over parents so that the first tasks are a
spread sample for the cost projection. Each depth-2 task's children are
written twice: `out.slim.state` without the DNN state (the frontier), and
`out.full.state` with it only when one of the children is a stage-B
candidate (below).

Resuming: a task with `done.json` is never rerun, so the stage is resumable
at task granularity; SIGTERM drains (no new task, running ones finish), a
second signal kills the running solvers (their tasks rerun). The work redone
at a stop is at most the tasks in flight. Sessions may overlap: a task is
claimed atomically before it runs (`tasks/<id>/claim`, stale once its process
is gone), and a slot waits until no other `bnbdrv` is pinned to its CPU, so
a new session can start while an old one drains its long cut phases without
idling the freed CPUs and without exceeding one solver per CPU. The last
depth runs in a seeded random order, and the projection estimates the
remaining tasks from the random-order ones only. Under the parents branched
on one line, the same child index is the same facility at the same location,
so a round-robin by index (the first 57 tasks of the wil50 run) is a
correlated sample.

The stage-B sample is fixed before stage A runs: the 400 frontier nodes
with the smallest `unit_hash(20261107, node)` (SHA-256 of the seed and the
node's assignments), a uniform random sample of the frontier whatever its
size. Stage A keeps a parent's DNN state when one of its children hashes
below tau = 0.01454 (430 expected below it of 29,580), so stage B can start
the sampled nodes from their parents' shipped state without 18 GB of
frontier. A sampled node whose parent was not kept (if the frontier is
smaller than about 27,500) has its parent rebuilt by rerunning its stage-A
task from the depth-1 parent (those 15 states are kept); stage B checks the
rebuilt children against stage A's record.

## Running

```
cd experiments/wil50_feasibility/sampling
python3 validate.py --instance ../../../qaplib/nug20.dat --sym ../../../latest/SYM/nug20.sym \
    --ub 2570 --grid 4:5 --depth 2 --ship 32 --tag nug20_ship32
setsid nohup python3 stage_a.py > runs/stageA.out 2>&1 < /dev/null &     # resumable; kill -TERM drains
python3 summarize_a.py                                                     # results/stageA_summary.md
python3 stage_b.py --limit 4 --check-F 1 --converge 2 --calibrate 49:23,15:37,5:6   # the smoke test
python3 stage_b.py                                                         # 400 nodes, the 80 farthest converged
```

## Validation (cloud box, Nov 7, binary built from 4d8693d's `latest/`)

`validate.py` runs `stage_a.py` (4 slots) and one serial
`bnbdrv -d 2 -o serial.state` with the same private binary and environment
(the settings plus the instance's grid, `QAP_DNN_CUT_GRID`, and the shipping
precision under test), then compares the frontiers node by node, the per-level
statistics (`n_level`, `f_level`, `l_level`, `e_level`, DNN counts) and every
DNN node's trace line. Reports in `results/validation/`.

| instance, depth 2 | shipping | frontier serial / parallel | same node set | per-level statistics | DNN nodes, same iterations / status | largest first-lb difference |
|---|---|---|---|---|---|---|
| nug20 (nug20.sym, -u 2570) | 64 (exact) | 18 / 18 | yes | equal | 36, 36 / 36 | 0.0017 (4 warm nodes) |
| nug20, symmetrisation switches off | 64 | 18 / 18 | yes | equal | 36, 36 / 36 | **0** |
| nug20 | 32 (production) | 18 / 18 | yes | equal | 36, 36 / 36 | 0.15 |
| nug30 (nug30.sym, -u 6124) | 32 (production) | 552 / 552 | yes | equal | 234, 228 / 234 | 0.29 |

[m: `results/validation/*.json`]

- **The frontier is the same node set as the serial run's** on both
  instances, with the production shipping (nug30: 552 nodes under 21
  branched depth-2 parents, 204 of 225 depth-2 nodes fathomed, 39 CP cut
  phases, every depth-2 node warm-started from its parent's shipped state).
- **Bit-identical with exact shipping once the iterates are symmetric.** With
  `QAP_DNN_SHIP=64` the first bounds of 4 warm-started nodes still differed by
  up to 0.0017. Rerun with `DNN_NOSYMW=0 DNN_CORE_NOSYMV=0`, everything is
  equal to the last printed digit, so the difference is the two production
  switches that skip symmetrisation passes. They leave the in-memory Z and Y
  slightly asymmetric, and a state file carries only the upper triangles.
  With fp32 shipping, 6 of 234 nug30 nodes take a different number of
  iterations, none a different outcome. A depth-2 node in the full run
  crosses a process boundary the same way whenever its unit returns, so this
  is the production behaviour, not an artefact of the launcher.
- **Restart against -F (3 frontier nodes each on nug20 and nug30, one on
  wil50t below):** the stage-B state file with no warm start gives exactly
  the first bound and iteration count of the same node started with `-F`
  (nug20: 2,569.7582 / 2,568.8553 / 2,569.1469; nug30: 6,125.6529 /
  6,127.2605 / 6,127.7758). With the parent's warm start, as stage B runs
  it, the same nodes fathom in 10-20 iterations instead of 30-66.
- **Stage B's paths** on nug30's frontier (`results/validation/nug30_stageB_test.md`):
  the give-up step, the converge step (4,000 iterations, 1.1-1.2 units above the
  300-iteration value), the -F check and the calibration path; and the rebuild
  of a parent whose state was not kept: its out.state equals stage A's in
  every line but the CPU times, the DNN state included, and the sampled
  child's give-up result is the same.

## Stage A on wil50t (cloud box, Nov 7 18:39 to Nov 8 12:57 UTC)

Binary: `latest/` at 4d8693d (the private copy's commit.txt says a667727, which
adds only this directory; `latest/` unmodified, no source newer than the
binary). Every number is measured [m: `results/stageA_summary.md`, built by
`summarize_a.py` from the 651 tasks' `done.json`], except where marked [e].

| depth | nodes bounded | fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB = all | cut phases |
|---|---|---|---|---|---|
| 0 | 1 (QPB only) | 0 | - | 1 | 0 |
| 1 | 15 | 0 | 4,500 + 0 | 1,762 + 0 + 13 = 1,775 | 0 |
| 2 | 635 | **331** (57 in the first solve, 274 in the cut phase) | 187,936 + 455,344 | 69,816 + 183,881 + 139 = 253,837 | 463 |
| all | 651 | 331 | 192,436 + 455,344 = 647,780 | **255,612** | 463 |

- **Frontier: 14,192 depth-3 nodes** (284 depth-2 parents x 48 + 20 rows-only
  parents x 28), against the 29,580 depth-3 orbits that `../roadmap/tree_size.md`
  assumes (0.48x). No child was eliminated in branching (`e_level` = 0); the
  difference is the 331 depth-2 fathoms.
- **The depth-2 nodes are near the cutoff, and the cut zone catches most of
  them.** After the first solve (300 iterations, warm from the depth-1 parent's
  shipped state), 463 of 635 were within the 244-unit zone and ran a CP cut
  phase: 274 fathomed (first-solve shortfall 4-171, median 84), 174 stalled
  (134-244 short at the start, 6.5-145 after, median 69), 15 hit the
  3,400-iteration limit (0.1-7.5 short after). Phases are 72% of the stage's CPU.
- **Shortfall histogram of the depth-2 nodes** (after the phase where one ran):
  fathomed 331; 0-50: 75; 50-100: 72; 100-150: 42; 150-244: 0 (a phase either
  fathoms or ends below 150); 244-300: 54; 300-350: 32; 350-400: 15; 400-450: 13;
  450-500: 1.
- **The frontier's inherited gap** (the parent's best bound; an upper bound on
  each node's own shortfall before it is bounded): 0-50: 3,500; 50-100: 3,336;
  100-150: 1,976; 244-300: 2,492; 300-350: 1,496; 350-400: 720; 400-450: 624;
  450-500: 48. So 8,812 frontier nodes (62%) are children of stalled phases and
  start with the parent's cut pool; 5,380 (38%) have parents 244-458 short.
- **Branching lines.** Root: row 21 (facility 21 in the 15 J1 cells; the root has
  QPB only, `-x 1:...`). Depth 1: 8 on a column (location 0, a corner, at
  `21:1`-`21:4`, `21:10`, `21:12`-`21:14`), 7 on a row (facility 45 at `21:0`,
  `21:11`, `21:20`, `21:21`; facility 31 at `21:22`-`21:24`). Depth 2: 183 on a
  column, 121 on a row (per-parent table in the summary; `frontier_stats.tsv`).
  Which rule chose them: see "Root lines" below.
- **Symmetry**: "read symmetry info from .../latest/SYM/wil50t.sym" in all 651
  solver logs. 15 depth-1 nodes; 635 depth-2 nodes before bounding (10 x 49 + 5 x
  29, the five row-2 cells branching rows-only), as `../sym/README.md` counts.
- **Root and depth-1 bounds against earlier values**: the root's QPB bound is
  46,100.957 (`../sym/logs/d1_sym.txt`: 46,101), and its 15 children's
  predicted bounds give a mean relative gap of 0.886607 (sd 0.025193), equal to
  that log's to the 6 printed digits, so the root branched on the same line as
  the Nov 3 run. The depth-1 DNN bounds (333-570 short at 300 iterations) are new;
  the laptop's 1,500-iteration cold bounds of the same 15 children (`2d26817`,
  `tuning/lines/root/children/f21_*`) are 55-62 higher, child for child.
- **Cost**: 255,612 CPU s against `tree_size.md`'s estimate of 1.6e5 (651 x 250
  s): 1.6x, from the phases. 18.3 h of wall on 4 slots; CPU/wall 0.990 per task.
- **Segments and redone work**: five sessions (34 tasks for the cost gate; then
  the switch to the random order; a 3-slot session beside a draining 25-minute
  cut phase; the claim files; tau raised to 0.04). About one slot-hour idled in
  drains before the claim files existed; one task (`r.5.17`) ran twice at once
  for about a minute when the random-order switch raced a drain, and both copies
  were killed and it was rerun. Nothing else was redone.
- **Kept for stage B**: 251 parents' full states (7.2 GB under `runs/stageA/`);
  the 400th smallest hash is 0.0287, so 9 parents of the 400 sampled nodes
  (finished before tau rose to 0.04) are rebuilt by stage B. Committed:
  `results/frontier.state.gz` (the frontier as a state file without DNN states,
  11 MB), `results/frontier.tsv` (every node: parent, inherited and predicted
  bound, hash), the trace lines, the launcher log, config and environment.

**What it changes in the estimate [e].** The tree model multiplies 29,580
depth-3 orbits by a per-orbit cost (1.0e4-6e4 DNN nodes under the Oct 29
settings, 85-91% of it from the far tail). The production frontier is 14,192
nodes. With the per-node mix unchanged that would divide the total by 2.08. But
the 331 fathomed depth-2 parents removed the depth-3 nodes nearest the cutoff,
which were the cheap ones, so the factor on the total is smaller than 2.08. What
it is depends on how many of the 14,192 are tail nodes (300 or more short at
convergence). The 2,888 frontier nodes whose parents were 300-458 short are the
candidates. Stage B's class fractions decide it. Stage A's own cost, 2.6e5 CPU s,
is still under 0.05% of a 6e8-CPU-s week.

## Root lines: forcing the root's facility (`--root-facility`, Nov 8)

The laptop's root study (`2d26817`, `tuning/lines/README.md`, "Root summary")
scored six root facilities by the lines cost of their 15 orbit children (cold,
1,500 iterations): facility 45 best, 21 at 1.2x, 25 (rule 8's pick with root
marginals) at 5.3x. This run's root is facility 21.

**Which rule picked 21.** The root has QPB only (`-x 1:...`), and in
`nug_deep_xbar9_all.par` the root's block is rule 9 (`BRANCH_XBAR_LOOKAHEAD`).
At the root `determine_J3` restricts the columns to J1 (15 cells), so
`can_branch_on_cols` is false and `select()` calls `select_row_lookahead`: the
QPB row lookahead of rule 4 over the 30 best rows by U, scored by each child's
QPB bound. That function never reads the DNN marginals (`dnn_xbar` appears only
in rule 8's selector and in `select_ind_lookahead`'s two-line mode). So the
root's line is a QPB lookahead whether or not the DNN runs at the root:
**D1 = 0 would give the same facility 21**, since branching uses the QPB bound
and U either way. The same holds at every rows-only node: the five depth-1
nodes in row 2 (`21:20`-`21:24`) and the 45 rows-only depth-2 nodes branch by
QPB row lookahead. Only the 10 other depth-1 nodes and the 590 other depth-2
nodes use rule 9's two DNN-marginal lines.

**How the line is forced.** The root task runs as before. Its 15 children (the
rule's facility in the 15 J1 cells) are rewritten to facility F in the same
cells, and the rule's version is kept as `out.rule.state`. Every other field
stays: depth 1, the inherited root bound, the root's warm start, the state's
root bound. The only change is the predicted bound (the root's bound plus U at
the rule's cell, which the state does not hold for F), which becomes the
inherited bound. It enters only `pick_param`'s relative gap, which stays above
0.2, so the first parameter block applies as it does for every depth-1 node.
It also enters a `max()` with the node's own bounds, which the DNN bound
exceeds. A child finds its warm start by comparing its `p` with the parent's
(`warm_child_offsets`), so F at c starts from the root's QPB point cut at
(F, c), like any root child.

Why not `bnbdrv -F F:c`: `stats->root_bound` is set only when a depth-0 node is
bounded, so under `-F` it stays 0 (the node-10 `-F` logs print "root bound: 0").
`pick_param` then computes the relative gap against 48,816 instead of against
48,816 - 46,101, about 0.05 instead of 0.89 at depth 1. That selects the third
parameter block (FW 100/50/25, 5 best rows) instead of the first (150/100/25,
30 best rows), so the QPB bounds and the lookahead differ. At depth 3 the two
coincide (the third block) for nodes within 326 of the cutoff; nodes farther
out get the second block in the full run and the third under `-F`. The tail
probes' `-F` roots (node 1 at 455 short, node 5 at 364) are such nodes. Their
DNN bounds are unaffected; their QPB parameters were not the full run's.

**Symmetry below a forced node.** `determine_J3` looks only at the node's depth
and its assigned locations: depth 0 means J1, and all assigned locations in row
2 means rows only (J3 = rows 0-2). A rewritten depth-1 node (depth 1, one
assignment) is therefore treated exactly like a rule's root child: 29 children
in the row-2 cells, 49 elsewhere, 635 in all. An `-F` node has depth equal to
its number of fixes (`QAPAssignment::fix` increments it), so it is never
treated as a root restricted to J1.

**Validation** (`validate_force.py`, `results/validation/force21.txt`): a run
forcing facility 21, restricted to the root, `r.0` [21:0], `r.4` [21:4], `r.14`
[21:24] (row 2: 29 children) and `r.14.0`, with a byte-identical binary,
reproduces stage A's tasks in every field compared: first lb 48,248.852,
48,411.3768, 48,480.9476 and 48,829.9474, the primal values, iterations,
statuses, the lines (row 45, col 0, row 31), the children's assignments and
their inherited and predicted bounds, and the per-level statistics.

**The samples** (`run_samples.sh`, detached; facility 45, then 25): the root,
the 15 depth-1 nodes of the forced line, and a seeded uniform sample of 200 of
the 635 depth-2 nodes (`--sample-last 200`). No depth-2 task starts until every
depth-1 task is done; then the 200 task ids with the smallest seeded order hash
are bounded, the same ids in every run. The same 200 ids in stage A's full run
form facility 21's matched sample. No parent state is kept (`--tau 0`). A run
drains (exit 4) if its projection exceeds 9 h of wall, and the chain then
stops. `compare_lines.py` writes the comparison (`results/root_lines.md`).

## Stage B (`stage_b.py`): built and smoke-tested, not run

- **The sample.** The 400 frontier nodes with the smallest
  `unit_hash(20261107, node)`, from stage A's records (`sample.tsv`). The seed
  was fixed before stage A ran; the hash makes the sample uniform over the
  frontier and independent of the order in which stage A produced it.
- **The give-up step**, one task per node, the way the full run treats a
  depth-3 node's first solve: a state file holding the node alone with its
  parent's warm start as stage A wrote it (the parent's QPB point, its cut
  pool when its phase ran, and its fp32 DNN state from `out.full.state`, or
  from a rebuilt parent), resumed with `bnbdrv -r node.state -d 3 -o out.state`
  under the settings and `-x 1:4:1:300:1`. So the node is warm-started from
  its parent's shipped state (a pooled node's first solve is on the inherited
  pool), and stops at the production give-up rule (300 iterations when
  hopeless) or fathoms. The cut phase is switched off at depth 3
  (`QAP_DNN_CUT_MAXDEPTH=2`, which gates nothing else): the strata are defined
  on the first-solve shortfall, and the phase (up to 3,400 iterations) would
  multiply the step's cost; `--phase` runs it as production would.
- **The converge step**: the 80 nodes farthest from the cutoff at the give-up
  point, from the same state file with no give-up (`-x 1:4:1:0:1`) and no cut
  phase (`QAP_DNN_CUT_ZONE=0`): the DNN runs to its own limit, 4,000 iterations,
  or its convergence test. On nug30 the 4,000-iteration bound is 1.1-1.2 units
  above the warm 300-iteration one. On wil50's clear misses of
  `../tuning/giveup/all32.txt`, the 1,500-iteration uncut values
  (`../stronger_bounds/README.md` step 2) are 50-74 above the cold
  300-iteration ones (`../tuning/giveup/cloud/H300_*.log`, Oct 31 settings).
- **The report** (`summary.md`): every node's first bound, shortfall, status,
  iterations, CPU, warm and pool flags and children; the converged bounds; the
  class fractions in the plan's strata (fathomed in the first solve, under 150
  short, 150-210, 210-300, 300 and more), each with a 1-sigma Wilson interval,
  classed by the converged shortfall where there is one and by the give-up
  shortfall (an upper bound on it) otherwise; the give-up shortfall quantiles.
- **Checks built in**: `--check-F K` reruns the first K nodes with no warm start
  and with `-F` (equal first bound required); `--calibrate SPEC` runs a `-F`
  node through both steps cold. Node 5 of `../tuning/giveup/all32.txt`
  (`49:23,15:37,5:6` on wil50t) has two recorded values: 48,450.4827 at 300
  iterations cold under these settings (363.6 short;
  `../tail_probe/runs_node5_cpt/cpcutt/L3_0.trace`, the same solver sources),
  and 48,498.4 after 1,500 uncut iterations (315.7 short; Oct 31 settings).
  The "329 short" quoted for it matches no log; `../roadmap/tree_size.md`
  section 3 calls it "the extrapolated 329", so it is model.py's
  extrapolation, not a measurement.

(TO BE COMPLETED: the smoke test and the projection.)
