# Pipeline A/B: the Nov 4 candidates on the node-10 subtree

The convergence and certificate studies (`../roadmap/convergence2/REPORT.md`,
`../roadmap/certificate/REPORT.md`) left four candidates. Each is tested
end to end on the node-10 subtree through the coordinator, against the
production configuration on the same machine at the same time.

| arm | what it adds | judged against | kill line (DNN iterations) |
|---|---|---|---|
| `base` | production: rule 9, bundle D, parent-Y warm start, pool inheritance (`../coord_search/run_coord_search.sh`) | - | - |
| `zone` | `QAP_DNN_CUT_ZONE=0.005` (from 0.003): misses up to 244 short run a cut phase, node 10 (232 short) among them, so their children inherit a pool (the "equip phase") | base | 1.15x |
| `warmy` | zone + Y kept between cut rounds (`QAP_DNN_CUT_WARMY=1`), 12 rounds of 200 iterations, quiet round at most 1,000 (Astra's Y12s), first Polyak step damped (0.3) for phases that start from an inherited pool | zone | 1.05x, and the tree no larger |
| `certx` | the certificate polish within 10 units of the cutoff (`DNN_CERTX=32 DNN_CERTX_USE=1 ...`) | base | 1.03x |

Iterations are the deciding measure: they do not depend on the machine.
Running two arms at once on one machine also makes their CPU seconds
comparable. The smoke runs on nug20 (`runs_smoke/`, not committed) checked
the plumbing only.

## Running it (laptop)

```
git pull origin cognac
make solver                      # latest/bnbdrv; on macOS with Accelerate
experiments/wil50_feasibility/pipeline_ab/run_pair.sh base zone
```

- Each arm gets half the performance cores (`SLOTS=N` overrides).
- The script keeps the Mac awake (`caffeinate`) while it runs. Keep the lid
  open, or the power and an external display connected.
- `kill PID` (the script prints it) or ^C stops it cleanly. The workers
  drain, and the same command resumes the arms where they stopped. A
  finished arm is skipped.
- At the end it writes `summary.txt` and prints the `git add/commit/push`
  line for the results: `runs/ARM/` trace, results, logs and hosts. The
  stores, work directories and binaries stay out of git.
- The second pair, `run_pair.sh warmy certx`, needs the first-step damping,
  which comes in a later commit (`git pull` and `make solver` first; the
  script refuses `warmy` on a binary without it).

`summarize.py runs/` can be run at any time; unfinished arms are listed as
such.

## Results

### base vs zone (laptop, Nov 4; `runs/base`, `runs/zone`, `summary.txt`)

Both arms at f6d448e, Apple M2 Pro, 3 slots each, run at the same time
(16:07-17:21 PST). **zone passes**: 1.635x fewer DNN iterations against a
kill line of 1.15x.

| arm | nodes / fathomed | cut phases (fathomed) | DNN iterations (first + phase) | DNN CPU | wall (3 slots) |
|---|---|---|---|---|---|
| base | 454 / 444 | 15 (6) | 83,640 (59,600 + 24,040) | 11,318 s | 4,431 s |
| zone | 366 / 358 | 8 (0) | 51,150 (40,650 + 10,500) | 7,985 s | 3,192 s |

- Iterations 1.635x, DNN CPU 1.417x, wall 1.388x. The CPU ratio is below
  the iteration ratio because the pooled iterations are dearer per step.
- The root's cut phase (zone only: 48,582 -> 48,685 in 1,200 iterations,
  still 129 short) fathoms nothing itself; the gain is in the children.
  Every zone node is a base node (0 nodes only in zone, 88 only in base).
  On the 366 shared nodes the iterations fall 74,390 -> 51,150 (1.454x):
  first solves are shorter from the inherited pool, and only 8 nodes need
  a phase of their own against 15.
- Two shared nodes change outcome, both to the good: a depth-5 node base
  branched after a 2,400-iteration phase (best 48,807.8) and zone fathomed
  in the first solve (48,814.1 at 1,570 iterations); a depth-7 node
  likewise (base phase best 48,813.4; zone first solve 48,814.8 at 190
  iterations). Those two account for the 88 extra base nodes.
- Against the cloud's Nov 3 runs of the same subtree (rule-9 arm 322 nodes,
  54,621 CPU-s on 4 Xeon workers), base here is 454 nodes at 11,456 CPU-s:
  the tree differs because the production configuration now carries the
  parent-Y warm start and pool inheritance, and the laptop is about 2.4x
  the Xeon per core without the AMX kernel.

Next: `run_pair.sh warmy certx` (warmy is judged against zone, kill line
1.05x and the tree no larger).
