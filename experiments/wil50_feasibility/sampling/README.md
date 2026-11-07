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
at a stop is at most the tasks in flight.

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

(TO BE COMPLETED: stage A results, stage B smoke test and projection.)
