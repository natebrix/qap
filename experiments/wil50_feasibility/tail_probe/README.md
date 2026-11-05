# Tail probe: do the Nov 4-5 levers hold on the tail class? (Nov 5)

**Why.** Every Nov 4-5 lever (`../pipeline_ab/`: a wider cut zone, warm-Y cut
rounds with a damped first step, the shipped parent state; together 4.3x
fewer DNN iterations) was measured on the node-10 subtree. Node 10 is in the
233-244-short class of depth-3 nodes. `../roadmap/tree_size.md` puts 85-91%
of the wil50 total in the 329-405-short class, whose subtrees are estimated
at 7e4-4e5 DNN nodes each. That is far too many to finish on a laptop.

**What.** Node 1 of `../tuning/lines/parents.txt` (`34:40,28:22,14:29`, 455
short at 300 iterations, 405 converged) is probed layer by layer
(`tail_probe.py`):

| level | what is bounded |
|---|---|
| 3 | the node itself |
| 4 | all of its children |
| 5 | a seeded sample of 30 grandchildren, stratified over their parents |
| 6 | a sample of 20 of their children |
| below | estimated as open nodes |

Each level's node count is estimated from the level above (mean children
per bounded node), and its iterations and CPU from the bounded sample. Every
task is one node in its own state file and solver process, as in a fleet
run.

| arm | settings | what it shows |
|---|---|---|
| `base` | the Nov 4 production setting | every node starts cold (no parent state across tasks) |
| `baseship` | base + the shipped parent state | every node warm, as in a serial run of base |
| `full` | base + cut zone 0.005 + warm-Y rounds (damped) + shipped state | the production candidate |
| `fullwide` | full with cut zone 0.010 | node 1 (455 short) runs a cut phase too, so its children inherit a pool |

Read the results like this:
- **full against baseship** isolates the cut-phase levers on the tail class.
- **baseship against base** bounds what warm starts are worth there.
- **fullwide against full** asks whether equipping the tail node pays.

The shipping arms ship in bf16 (`QAP_DNN_SHIP=16`), so that each task file
carries about 13 MB of parent state instead of 26. On 8 children of node 10,
bf16, fp32 and exact gave the same iterations (`../ship/`).

## Running it (laptop)

```
git pull origin cognac && make solver
experiments/wil50_feasibility/tail_probe/run_probe.sh base full
experiments/wil50_feasibility/tail_probe/run_probe.sh baseship fullwide
```

- **Slots:** each arm gets half the performance cores (`SLOTS=N` overrides).
- **Interruptions:** `kill` with the PID it prints, or ^C, stops it. Run
  the same command again: finished tasks are kept, and the run resumes.
- **Time:** an estimate of 20-60 minutes per pair. The 47 level-4 children
  dominate, more in the arms whose zone gives them cut phases.
- **Disk:** about 1 GB per shipping arm under `runs/`. The state files and
  logs are ignored by git; the traces and summaries are committed.
- **Output:** the end prints `compare.txt` (estimated iterations per level
  and in total, the open count below the probe, and the ratios against base
  and baseship) and the `git add/commit/push` line.

The probe was smoke-tested on nug20 from a depth-1 node: the layers, the
sampling, resuming, the warm counts (0 of 19 in base, 19 of 19 with shipping)
and the comparison. Its cut phases end in `status error` there because the
test passed nug20's grid for the wrong matrix. That is a smoke-test setting;
no wil50 trace has such a phase (0 of 87).

## Results

(pending)
