# Phase 2 gate, cloud box: wil50 depth-2 replays, stage-A binary against HEAD (Nov 9)

`docs/settings_plan.md`, Phase 2, "Gate, cloud box (after stage B)". The 8 depth-2 tasks of stage A
were replayed on the machine stage A ran on, from their kept depth-1 parent states. The depth-1
task r.1 and the root task `r` were replayed too. Records: `phase2_gate_wil50/<side>/<task>/`
(trace.txt, solver.log, done.json) and `phase2_gate_wil50/compare.json`.

## Verdict

**The settings are identical.** On all 8 depth-2 tasks and r.1, HEAD under `fullrun.par`
reproduces stage A's record on every DNN result:
- **Trace lines:** identical field for field, pid and CPU aside. That covers the first bound,
  primal, status, iterations, warm, pool, the cut phase (best bound, rounds, iterations, outcome),
  cutoff, fathomed and pool_out.
- **Children:** the same count, assignments and inherited bounds. The inherited bound is the
  final DNN bound in full double precision (e.g. 48677.80197002722).
- **Branching line:** the same on every branching task.
- **Shipped state:** the cut pool and the 28.3 MB shipped DNN state carried to the children are
  byte-identical.

**Not identical: the QPB values, and the cause is 18ebd58, not the settings.**
- **What differs:**
  - the node's Frank-Wolfe point in the children's warm start (all 2,304 entries, |dX| up to 0.07);
  - the children's predicted bounds (up to 33 apart);
  - at the root, the QPB bound itself: 46100.957 becomes 46131.508.
- **The cause** is 18ebd58, the QPB h-term fix. It postdates a667727, and the plan has it under
  "Out of scope" to be gated on its own.
- **Proof:**
  - 18ebd58, run with stage A's environment and flags, gives HEAD's output byte for byte (traces,
    logs, every node line, every warm line) on `r`, r.1 and r.6.41.
  - Its parent 347ff84 gives stage A's output byte for byte on `r` and r.1.
  - 347ff84 to 18ebd58 changes `latest/qpbnd.C` alone.

Against a binary that carries 18ebd58, the children match too. That was shown byte for byte with
18ebd58 itself on r, r.1 and r.6.41, and the laptop's fe14ac0 replays showed it at depth 1. Against
a667727, the gate's wording ("the same children") fails, and only by the QPB values. All six tasks
with children (the 4 branching depth-2 tasks, r.1 and the root) show the same pattern. Three of
them were attributed by experiment. Before the
full run, the 18ebd58 decision is the one open item in "Re-gate the full run":
- **What a fresh run would change:** under HEAD, a fresh run from the root starts with depth-1
  nodes inheriting 46131.5, not 46101. Every QPB point and predicted bound below them can move.
- **The line choice can move too:** rule 9's lookahead scores lines by the children's QPB bounds.
- **What these replays don't test:** both binaries got the same a667727-made parent state, so they
  don't show whether the new upstream QPB points move the DNN path.
- **What they do show:** the search uses a node's predicted bound only as
  `max(inherited, predicted)`, in `pick_param` (bnb.C:240) and in the node's bound (bnb.C:610).
  - **Depth 2 and 3:** every predicted bound is 1,300 to 2,000 below the inherited one, so its
    change decides nothing there.
  - **Depth 1 (the root's children):** the predicted bound is above the inherited one and moves by
    at most 11.2. These nodes stay in the first parameter block.
  - **What can matter:** the QPB point X warm-starts the child's QPB.
- **The root bound moves the parameter blocks.** This was computed, not run: `pick_param`'s
  rule applied to stage A's recorded bounds.
  - The root bound is the denominator of every node's relative gap, `(opt - bnd)/(opt - root)`.
    So under 18ebd58 the 0.2 threshold moves from bound 48272.99 to 48279.10, and the 0.12
    threshold from 48490.19 to 48493.86.
  - In stage A's tree, all 49 depth-2 children of r.1 (bound 48273.088, gap 0.19996 against
    0.2022) would move from block 1 (FW `150/100/12`, best 3) to block 0 (`150/100/25`, best 30).
    Stage A's logs show block 1 for all 49.
  - 124 depth-3 frontier nodes (the children of r.10.20, r.5.28 and r.7.44, bounds 48492.4 to
    48493.9) would move from block 2 to block 1.
  - The replays don't exercise this, because both sides read the a667727 root bound from the
    input state's header.
- **The 4 branching depth-2 tasks:** the line was unchanged, and so was every DNN result.

The old binary reproduces its own record on this machine. That held for `r`, r.3.36 and r.6.41,
the DNN state included.

## Tasks

Each pick is the cheapest of its class by stage A's CPU. For the stalled and iteration-limit
classes, a row-line case at the same cost replaced the second-cheapest (r.2.26 for r.3.33, r.13.4
for r.0.15). The 8 tasks have 7 depth-1 parents: column-0 parents r.2, r.3, r.4, r.7 and r.9, the
row-45 parent r.6, and the row-31 parent r.13.

| class (stage A count) | tasks |
|---|---|
| fathomed in the first solve (57) | r.9.39, r.4.39 |
| fathomed in the CP cut phase (274) | r.3.36, r.7.35 |
| CP phase stalled (174) | r.6.41, r.2.26 |
| CP phase hit the 3,400-iteration limit (15) | r.3.15, r.13.4 |

**There is no task r.20.** The root has 15 children, r.0 to r.14, and the laptop's r.10 is cell
21:20. In its place I ran r.1 (21:1, a column-0 line, 49 children), because the laptop's three
replays were all row lines. I also ran the root `r` (QPB only, 1 s).

## Commands

- **New binary:** HEAD c378cad, built in this worktree with
  `make solver CF="-O3 -march=native -ffp-contract=off -include cstring"`. This is the makefile's
  `CF` plus a forced `<cstring>`, because plain `make solver` fails with g++ 13.3 (see
  "Surprises"). sha256 `50e3df6d...`, linked against the same scipy_openblas32 as stage A's
  binary. `fullrun_check` passes on it.
- **Old binary:** stage A's `bin/bnbdrv` (a667727, sha256 `e25b0286...`), run in place.
- **Attribution binaries:** `git archive 347ff84 latest` and `git archive 18ebd58 latest` in the
  scratchpad, each built with its own makefile (`make bnb`).
- **Input states:** each task's one-node input was built as stage A built it, with
  `common.write_one(stageA/tasks/<parent>/out.state, <index>, in.state)`, warm start kept. All 14
  inputs are byte-identical to a667727's `write_one`. HEAD reads the version-6 parent states.
- **New side:** the sampling code path, as `stage_a.py` runs a task. The environment is
  `fullrun_env` (the thread variables and `QAP_DNN_POOL_CHECK=20`). The params file is
  `params.par = params_with(fullrun.par, {})`, which holds the same keys, comments dropped. Each run
  goes through `run_bnbdrv` with `taskset -c <cpu>` and `QAP_DNN_TRACE`:
  `bnbdrv -p params.par -s latest/SYM/wil50t.sym -u 48816.0 -d <depth> -o out.state -r in.state qaplib/wil50t.dat`
- **Old, 347ff84 and 18ebd58 sides:** stage A's own record. The environment is `config.json`'s
  `settings_env` (31 variables, neither `QAP_DNN_BETA` nor `QAP_DNN_GU_ZONE`), and the flags are
  stage A's:
  `bnbdrv -p latest/PARAM/nug_deep_xbar9_all.par -x 1:4:1:300:1 -s ... -u 48816.0 -d <depth> -o out.state -r in.state wil50t.dat`
  The root runs without `-r`, at `-d 0`.
- **Exact commands:** each `solver.log`'s first line has the exact command.
- **Comparison:**
  - each replay's trace line, with pid and cpu removed, against stage A's `tasks/<id>/trace.txt`;
  - the log's node, iteration and DNN counters;
  - every node line, and every warm line (QPB point, `dnn_ext`, cut pool, DNN state) of the
    output, against stage A's `out.full.state` (all four branching tasks were retained) or the
    root's or r.1's `out.state`.

## Old vs new

The stage A rows are the task records. The HEAD rows are the replays. The first bound and the cut
phase's best bound are the trace's `%.4f`. The children's bounds come from the state files at full
precision.

| task | class | side | first lb | status, iters | cut phase: outcome, best, iters, rounds | fathomed | children, line | children's inherited bound | children's predicted bounds | CPU s |
|---|---|---|---|---|---|---|---|---|---|---|
| `r.9.39` `21:14,40:0` | first solve | stage A | 48816.2166 | fathomed 23 | - | 1 | 0, - | - | - | 10.2 |
| | | HEAD | 48816.2166 | fathomed 23 | - | 1 | 0, - | - | - | 9.9 |
| `r.4.39` `21:4,40:0` | first solve | stage A | 48814.3576 | fathomed 28 | - | 1 | 0, - | - | - | 11.0 |
| | | HEAD | 48814.3576 | fathomed 28 | - | 1 | 0, - | - | - | 10.6 |
| `r.3.36` `21:3,37:0` | CP fathomed | stage A | 48793.2235 | hopeless 300 | fathomed 48814.3827, 78 it, 1 rd | 1 | 0, - | - | - | 150.0 |
| | | HEAD | 48793.2235 | hopeless 300 | fathomed 48814.3827, 78 it, 1 rd | 1 | 0, - | - | - | 150.4 |
| `r.7.35` `21:12,36:0` | CP fathomed | stage A | 48797.5467 | hopeless 310 | fathomed 48814.5727, 90 it, 1 rd | 1 | 0, - | - | - | 153.3 |
| | | HEAD | 48797.5467 | hopeless 310 | fathomed 48814.5727, 90 it, 1 rd | 1 | 0, - | - | - | 150.2 |
| `r.6.41` `21:11,45:42` | stalled | stage A | 48580.3874 | hopeless 310 | stalled 48677.8020, 860 it, 9 rd | 0 | 48, col 9 | 48677.80197002722 | 46717.28 - 47115.65 | 418.6 |
| | | HEAD | 48580.3874 | hopeless 310 | stalled 48677.8020, 860 it, 9 rd | 0 | 48, col 9 | 48677.80197002722 | 46713.98 - 47100.70 | 429.0 |
| `r.2.26` `21:2,27:0` | stalled | stage A | 48582.2542 | hopeless 300 | stalled 48678.2352, 860 it, 9 rd | 0 | 48, row 45 | 48678.23521765621 | 46738.52 - 47193.46 | 432.4 |
| | | HEAD | 48582.2542 | hopeless 300 | stalled 48678.2352, 860 it, 9 rd | 0 | 48, row 45 | 48678.23521765621 | 46743.33 - 47186.20 | 408.7 |
| `r.3.15` `15:0,21:3` | iter limit | stage A | 48682.9782 | hopeless 300 | iter_limit 48810.7809, 3400 it, 35 rd | 0 | 48, col 9 | 48810.78086665441 | 46856.51 - 47165.75 | 1355.0 |
| | | HEAD | 48682.9782 | hopeless 300 | iter_limit 48810.7809, 3400 it, 35 rd | 0 | 48, col 9 | 48810.78086665441 | 46854.15 - 47167.76 | 1374.2 |
| `r.13.4` `21:23,31:4` | iter limit | stage A | 48653.0738 | hopeless 310 | iter_limit 48806.5861, 3400 it, 35 rd | 0 | 48, row 45 | 48806.5861089415 | 46942.10 - 47186.46 | 1444.0 |
| | | HEAD | 48653.0738 | hopeless 310 | iter_limit 48806.5861, 3400 it, 35 rd | 0 | 48, row 45 | 48806.5861089415 | 46930.61 - 47184.08 | 1410.7 |
| `r.1` `21:1` | depth 1 | stage A | 48273.0884 | hopeless 300 | - | 0 | 49, col 0 | 48273.08840905269 | 46464.04 - 46952.95 | 117.8 |
| | | HEAD | 48273.0884 | hopeless 300 | - | 0 | 49, col 0 | 48273.08840905269 | 46463.84 - 46957.36 | 106.4 |
| `r` `root` | root, QPB | stage A | - | QPB only | - | 0 | 15, row 21 | 46100.957091528995 | 46313.64 - 46551.18 | 1.1 |
| | | HEAD | - | QPB only | - | 0 | 15, row 21 | 46131.50786271595 | 46316.82 - 46552.00 | 1.0 |

Per child, the largest |new - old| predicted bound is:

| task | r.6.41 | r.2.26 | r.3.15 | r.13.4 | r.1 | r |
|---|---|---|---|---|---|---|
| max \|new - old\| | 20.98 | 27.57 | 33.14 | 26.69 | 21.07 | 11.20 |

In each warm start, every entry of the QPB point X differs (2,304 entries, max |dX| 0.029, 0.063,
0.070 and 0.037 for the four depth-2 tasks). The `dnn_ext`, the cut pool (73-85 KB) and the DNN
state (28,348,443 characters) are identical. At depths 2 and 3, every predicted bound is 1,316 to
1,964 below the inherited bound.

| attribution run | task | against | result |
|---|---|---|---|
| a667727 (old) | r, r.3.36, r.6.41 | stage A's record | identical (traces, counters, node and warm lines) |
| 347ff84 (18ebd58's parent) | r, r.1 | stage A's record | identical |
| 18ebd58 | r, r.1, r.6.41 | HEAD's replay | identical |
| 18ebd58 | r, r.1, r.6.41 | stage A's record | the same QPB differences as HEAD |

## Surprises

- **HEAD does not build on this box.** g++ 13.3 (Ubuntu 24.04) stops at
  `util.C:133:13: error: 'strcmp' was not declared in this scope` (`util.C` doesn't include
  `<cstring>`; the laptop's clang gets it transitively). The line came with 0ac0209. Here it was
  built with `-include cstring` and the source unchanged. The full run's Linux machines need the
  include.
- **`par.hash` is the FNV-1a of the params file's bytes**, and `params_with`/`envpar.write_par`
  drops `fullrun.par`'s comment block.
  - So the sampling launchers' `params.par` records `par.hash da0f5610dd2b176e` and settings
    hash `d39d7b7e50ca1315`. The golden has `3f2f813daeeb4a6b` and `ebd9848cb3f7e0d2`. Every value
    is identical (the dumps differ in those two lines only).
  - Consequence: a version-8 state written under one file does not resume under the other.
    Resuming `new/r.9.39`'s output state with `-p fullrun.par` exits 1 with
    `read_state: [...] the state file was written with settings [d39d7b7e50ca1315 ... par.hash=da0f5610dd2b176e], this run has [ebd9848cb3f7e0d2 ... par.hash=3f2f813daeeb4a6b] ...; --allow-settings-change to resume anyway`.
  - With `-p params.par`, the same state resumes.
  - The version-6 stage A and B states have no record and are unaffected. Before the full run,
    decide whether its launcher and the sampling launchers share one file byte for byte.
- **Old and new CPU agree within 10%** (r.1: 106 s against 118 s; r.2.26: 409 s against 432 s).
  The iteration counts are identical, so this is machine noise.

## CPU

The replays took 5,289 CPU s in 18 runs:

| side | CPU s | runs |
|---|---|---|
| HEAD | 4,051 | 10 |
| old binary | 585 | 3 |
| 347ff84 | 112 | 2 |
| 18ebd58 | 541 | 3 |

The three builds took about 10 CPU minutes more. The runs took 28 minutes of wall time on 4
pinned CPUs. The output states (about 30 MB per branching task) were deleted after the
comparison.
