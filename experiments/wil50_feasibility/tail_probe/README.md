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

(laptop, Nov 5, 08:35-11:20 PST; binary at 8632a46, 3 slots per arm, `base
full` then `baseship fullwide`; `runs/<arm>/summary.txt`, `runs/compare.txt`,
one `.trace` per bounded node)

### Level 4: all 47 children of node 1, the same nodes in every arm

| arm | fathomed | cut phases | DNN iterations | CPU | iterations per child (first + phase) |
|---|---|---|---|---|---|
| base | 16 / 47 | 17 | 33,780 | 6,600 s | 719 (340 + 379) |
| baseship | 17 / 47 | 25 | 47,590 | 9,197 s | 1,013 (307 + 706) |
| full | 18 / 47 | 41 | 53,900 | 10,515 s | 1,147 (307 + 840) |
| fullwide | 18 / 47 | 38 | 49,790 | 9,420 s | 1,059 (341 + 718) |

This is the exact part of the probe. On the tail class the children of a
455-short node are 150-250 short themselves, so the levers that paid on
node 10 behave differently here:

- **Warm starts cost at level 4.** baseship's children start warm and
  their first solves land higher (48,680 against 48,641 on the flipped
  nodes), which puts 25 of them inside the 0.003 zone instead of 17; each
  phase is 2,000-3,900 iterations and fathoms almost nothing (1 of 25).
  Warm first solves are the same 300 iterations (the 300 are a floor, not
  a convergence test here). Net +41% iterations for one more fathom.
- **The cut levers cost more still.** full runs 41 phases (every child
  within 244 short), 18 fathom; fullwide's root phase (800 iterations,
  48,582 -> about 48,650) gives the children a pool and they run 38
  phases at 1,059 each. Phase iterations dominate level 4 in every arm
  but base: 70-75% of the level's cost.
- **Four marginal children flip** between arms, none consistently: the
  phases end at 48,742-48,814 on them, and which arm's phase crosses the
  cutoff is not predictable from the setting (fullwide fathoms 13:45 and
  5:45 that full does not, full fathoms 43:45 and 6:45 that fullwide does
  not). These are the 49:45-style marginal nodes of `../pipeline_ab/`, now
  at level 4.

### Level 5: 30 sampled grandchildren per arm, mostly different nodes

The sample is stratified over the level-4 nodes that branched, which
differ by arm, so only 5 of the 30 are shared by all four arms. Per-node
means across each arm's own sample:

| arm | fathomed | iterations per node (first + phase) | CPU per node | phases | pooled |
|---|---|---|---|---|---|
| base | 25 / 30 | 537 (267 + 270) | 94.9 s | 6 | 1 |
| baseship | 29 / 30 | 295 (105 + 190) | 52.4 s | 5 | 8 |
| full | 29 / 30 | 186 (88 + 99) | 23.7 s | 4 | 23 |
| fullwide | 29 / 30 | 166 (129 + 37) | 30.4 s | 3 | 30 |

On the 5 shared nodes: base 4 of 5 fathomed in 5,660 iterations, baseship
4 of 5 in 5,040, full 5 of 5 in 2,660, fullwide 5 of 5 in 1,790. The one
flip (`1:45 21:21`) is the warm-Y lever doing exactly what it did on node
10: the same warm first solve (48,728.2 in baseship and full), then the
cold-Y phase reaches 48,773 and the warm-Y damped phase 48,814.1.

Level 5 is where the levers pay on the tail class: warm starts about 1.8x
per node (the first solve 267 -> 105 iterations), the cut levers another
1.6x on top (phases 190 -> 99, and 23 of 30 nodes arrive with a pool),
with the fathom rate 25 -> 29 of 30. The per-node means are on different
node sets, so the ratios are indicative, not matched.

### Level 6 and the totals: not comparable across arms

Each arm's level-6 sample is the children of its own branched level-5
nodes. base branched 5 of 30; the other three arms branched 1 of 30 each,
and not the same one: baseship's is `1:45 21:21` (children `49:x`), full's
and fullwide's is `0:45 10:21` (children `40:x`), a harder node whose
children run 300 iterations to a hopeless 48,705-48,776 before their
phases. So baseship's 20 of 20 at 38 iterations and full's 17 of 20 at 570
are different parents, not the pool hurting at depth 6, and the compare.py
totals (baseship 4.17x over base, full 0.37x of baseship, the 13,207 open
at level 7) are driven by this: they should not be quoted. The probe needs
a matched design below level 4 (sample the level-5 and level-6 nodes once,
from base's branching, and bound the same nodes in every arm) before its
totals mean anything.

### Reading

1. On the tail class the production candidate's gains are at level 5 and
   presumably below, where nodes arrive 50-100 short with a pool and a
   warm start and fathom in 100-200 iterations. They are not at level 4,
   where the children are still 150-250 short: every lever adds cut
   phases there that cost 2,000-4,000 iterations each and fathom one or
   two more of 47. Level 4 is 2-4% of the subtree's estimated cost in any
   arm, so this does not change the verdict, but the zone should not be
   widened further for the tail class (fullwide against full: same
   fathoms, root phase +800 iterations, level-4 phases cheaper by 11%,
   level 5 by 11%; a wash within the sampling noise).
2. The node-10 ratios do not transfer as a single factor. Level 4 is
   0.63x (full against base), level 5 about 2.9x on unmatched samples,
   and the subtree total depends on the level-6 and deeper nodes that
   this probe did not measure comparably.
3. Next: the matched version of the probe (same sampled nodes in every
   arm), and a larger level-5 sample (the 30 of 1,334-1,426 gives one
   flip of 5 shared nodes to argue from).

### How big is node 1's subtree? (Nov 5, a second reading)

The arms cannot be compared below level 4, but each arm's own estimate is
still an estimate of the tail class's size. base, the Nov 4 production
setting with every node started cold, puts node 1's subtree at about
**12,000 DNN nodes and 2.2e6 DNN iterations** (1 + 47 + 1,426 + 10,695),
337,530 laptop CPU s, or about 6-9e5 Xeon s (the Xeon is 1.8x the laptop
on nug30's record run and about 2.7x per wil50 iteration without AMX) [e].

- **The roadmap's tail assumption.** `../roadmap/tree_size.md` built 7e4-4e5
  DNN nodes per 329-405-short node with the Oct 29 settings and anchors
  (node 10 = 2,193 nodes). The roadmap then applies node-10 factors to the
  whole tree; under the Nov 4 setting node 10's subtree is 322-454 nodes
  (`../pipeline_ab/`), 4.8-6.8x fewer, so the roadmap implicitly puts the
  tail at about 1e4-8e4 nodes today.
- **The probe's point estimate sits at the low end of that range.** The tail
  class shrank under the Nov 4 levers at least as much as node 10 did, not
  less.
- **Its uncertainty is wide.** Level 6 rests on the children of the 5 of 30
  level-5 nodes that branched (95% upper bound 35%), and level 7 on 0 of 20
  level-6 nodes branching (95% upper bound 17%). At those bounds the subtree
  would be about 1e5-2e5 nodes, 8-16x the point estimate.

So the first run moves nothing in the roadmap by itself; the matched run
below, with three times the level-5 sample and every level sampled until
nothing is open, measures it.

## The matched probe (`matched.py`, `run_matched.sh`)

**Design.** Levels 3 and 4 are as before, and are reused: `matched.py` runs
in the same `runs/ARM/` directories, so the first run's `L3_0` and `L4_*`
tasks are skipped. Below level 4, of each branched sampled node's children
it bounds the k = 3 with the smallest hash of their assignment, level after
level, until nothing is open or level 9. The arms branch the same level-4
parents on the same line in nearly every case, so a hash pick bounds the
same nodes in every arm where the parent branched. A first-run level-5 or
level-6 task that is the very node picked is reused (hard link).

**Estimates** are Horvitz-Thompson (a sampled node stands for its parent's
weight times the parent's open children over the number picked), with an
ultimate-cluster standard error stratified by level-4 parent. The
comparison prints each arm's totals with their SE and, for every pair of
arms and level, the nodes bounded in both: fathoms and iterations side by
side.

**Running it (laptop).**

```
git pull origin cognac && make solver
experiments/wil50_feasibility/tail_probe/run_matched.sh base full
experiments/wil50_feasibility/tail_probe/run_matched.sh baseship fullwide   (optional)
```

- **Time:** about 90 level-5 nodes per arm, then about 3 per branched node
  below. base's level 5 is about 9,000 CPU s on 3 slots (about 50 minutes),
  full's about a quarter of that [e]. The first pair should take about an
  hour, the second (optional: fullwide was a wash against full) less.
- **Disk:** a parent's state file (the shipped DNN state, about 13 MB per
  file) is cut down to its open count once its picked children are done.
  The first run's level-5 and level-6 state files are no longer needed
  except where reused; `rm runs/*/L5_* runs/*/L6_*` after the first pair
  frees them (the traces are in git).
- **Output:** `runs/ARM/matched.txt` (per level and total),
  `runs/matched_compare.txt` (all arms, matched nodes), and the git line.

Smoke-tested on nug20 from a depth-1 node at a loose cutoff (4 levels, 38
matched nodes at level 3 and 16 at level 4 in both arms), including resuming,
the reuse of an earlier task, and the run script.
