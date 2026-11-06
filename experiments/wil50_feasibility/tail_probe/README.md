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

**Superseded by the matched results below:** with three times the level-5
sample and level 6 branching 2 of 24 (not 0 of 20), base's estimate is
27,000 nodes, not 12,000, and sits at the center of the roadmap's implied
range, not its low end.

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

### Matched results (laptop, Nov 5, 12:06-14:00 PST; binary at 8a9033c + the fp32-R commit, default off)

`base full` then `baseship fullwide`, 3 slots per arm, levels 3-4 reused
from the morning run; `runs/<arm>/matched.txt`, `runs/matched_compare.txt`,
one `M<level>_<hash>.trace` per bounded node.

Horvitz-Thompson estimates of node 1's subtree:

| arm | est. DNN nodes | est. iterations (± 1 SE) | est. CPU s (± 1 SE) | per level: nodes, fathomed/bounded |
|---|---|---|---|---|
| base | 27,234 | 12,572,760 (± 9,292,023) | 1,094,282 (± 742,452) | L4: 47, 16/47; L5: 1,426, 85/93; L6: 5,520, 22/24; L7: 20,240, 6/6 |
| baseship | 27,878 | 9,009,763 (± 7,306,724) | 883,997 (± 665,636) | L4: 47, 17/47; L5: 1,380, 81/90; L6: 6,210, 25/27; L7: 20,240, 6/6 |
| full | 14,262 | 2,552,000 (± 1,778,097) | 353,678 (± 238,944) | L4: 47, 18/47; L5: 1,334, 83/87; L6: 2,760, 11/12; L7: 10,120, 3/3 |
| fullwide | 13,572 | 2,320,683 (± 1,631,517) | 360,723 (± 251,233) | L4: 47, 18/47; L5: 1,334, 84/87; L6: 2,070, 8/9; L7: 10,120, 3/3 |

Matched nodes (the same assignment bounded in both arms), fathomed and
DNN iterations, first arm / second arm:

| arms | L4 (47 shared) | L5 (shared) | L6 (shared) | L7 (shared) |
|---|---|---|---|---|
| base / baseship | 16 / 17, 33,780 / 47,590 (0.71x) | 90: 82 / 81, 42,810 / 35,970 (1.19x) | 24: 22 / 22, 5,580 / 3,780 (1.48x) | 6: 6 / 6, 3,140 / 2,230 (1.41x) |
| base / full | 16 / 18, 33,780 / 53,900 (0.63x) | 87: 79 / 83, 41,360 / 15,450 (**2.68x**) | 12: 10 / 11, 3,340 / 3,230 (1.03x) | 3: 3 / 3, 2,750 / 450 (6.1x) |
| baseship / full | 17 / 18, 47,590 / 53,900 (0.88x) | 87: 78 / 83, 35,780 / 15,450 (**2.32x**) | 12: 10 / 11, 3,070 / 3,230 (0.95x) | 3: 3 / 3, 2,130 / 450 (4.7x) |
| full / fullwide | 18 / 18, 53,900 / 49,790 (1.08x) | 81: 77 / 78, 15,220 / 15,710 (0.97x) | 9: 8 / 8, 3,080 / 2,790 (1.10x) | 3: 3 / 3, 450 / 410 (1.10x) |

Reading, in the README's three questions:

- **full against baseship (the cut-phase levers on the tail class): 2.3x
  per node at level 5 on 87 matched nodes, 83 fathomed against 78.**
  Level 4 costs 12% more (41 phases against 25, one more fathom). Level 6
  is even on 12 nodes and level 7 4.7x on 3; too few to weigh. End to end
  3.5x in estimated iterations and 2.5x in CPU, driven by level 5 and by
  the level-5 branching rate (4 of 87 against 9 of 90), which halves the
  estimated nodes below (14,262 against 27,878).
- **baseship against base (warm starts): 1.2x at level 5, 1.4-1.5x at
  levels 6-7, 0.71x at level 4** (the warm first solves put 25 children in
  the zone instead of 17). The fathom counts are the same at every level.
  End to end 1.4x in iterations, 1.24x in CPU: on the tail class the warm
  start alone is worth less than on node 10 (1.66x), because the children
  are far from the cutoff and the warm start saves iterations only where
  the node fathoms in its first solve.
- **fullwide against full (equipping the tail node): a wash.** The root
  phase costs 800 iterations, level 4 is 8% cheaper, level 5 3% dearer,
  the same fathoms everywhere (18, 78 against 77, 8, 3). The zone stays at
  0.005.

Against base, the production candidate is 4.9x in estimated iterations
and 3.1x in CPU on this tail-class subtree, with the level-5 matched
figure (2.68x on 87 nodes, 83 against 79 fathomed) the solid part. That
is in line with node 10 (4.26x on the whole subtree), so the Nov 4-5 levers
hold on the tail class.

Size of the tail class: base's estimate is 27,000 DNN nodes (SE 74% in
iterations). Level 7 is 20,240 estimated nodes from 6 bounded, none
branching, and carries 84% of base's estimated iterations at 523 per node
(one of the 6 ran a cut phase); level 6 is 5,500 from 24 bounded, 2
branching. The point estimate is 2.3x this morning's 12,000 because the
larger level-5 sample found 8 of 93 branching against 5 of 30, and the
level-6 branching is now 2 of 24 against 0 of 20. Under full the subtree
is about 14,000 nodes and 2.5e6 iterations, 3.5e5 laptop CPU s. The
uncertainty is still the level-6 and level-7 branching rates: with 2 of
24 and 0 of 6, the 95% upper bounds (about 25% and 40%) would put the
subtree several times higher. A deeper sample under full (k = 3 is too
few at levels 6-7) is the next step if the tail class's size matters for
the roadmap; the lever ratios do not need it.

### Matched results, node 5: the 329-short node (laptop, Nov 5, 21:02-21:48 PST; binary at 151aded)

The other member of the tail class is node 5 of `../tuning/giveup/all32.txt`
(`23:49,37:15,6:5`: 329 short at 300 iterations, 315.7 converged, against
node 1's 455 and 404.8). The same matched probe, base against full, 3 slots
each, in `runs_node5/` (`PROBE_RUNS=.../runs_node5 PROBE_ARGS="--node
23:49,37:15,6:5" run_matched.sh base full`). Both arms finished in 46
minutes of wall, against about two hours each for node 1.

| arm | est. DNN nodes | est. iterations (± 1 SE) | est. CPU s (± 1 SE) | per level: nodes, fathomed/bounded |
|---|---|---|---|---|
| base | 1,290 | 381,930 (± 159,454) | 43,955 (± 15,861) | L4: 47, 35/47; L5: 552, 35/36; L6: 690, 3/3 |
| full | 1,290 | 161,253 (± 62,995) | 24,204 (± 9,150) | L4: 47, 35/47; L5: 552, 35/36; L6: 690, 3/3 |

Matched nodes (the same assignment bounded in both arms):

| arms | L4 (47 shared) | L5 (33 shared) | L6 (3 shared) |
|---|---|---|---|
| base / full | 35 / 35, 26,510 / 25,100 (1.06x) | 32 / 32, 16,060 / 3,980 (**4.04x**) | 3 / 3, 440 / 170 (2.59x) |

- **The levers hold here too, more cleanly than on node 1:** 2.37x in
  estimated iterations and 1.82x in CPU on the subtree, 4.04x per node at
  level 5 on 33 matched nodes (every one fathomed in both arms). Level 4 is
  a wash (1.06x) instead of node 1's 0.63x penalty: the 12 branched
  children are the only ones whose phases run to the end in both arms, and
  the 35 fathoms are cheap in both. The identical trees (35 of 47 fathomed
  at level 4, 35 of 36 at level 5, 3 of 3 at level 6, the same nodes) mean
  the whole difference is per-node cost, as on node 10 and node 1.
- **Node 1 does not stand for the class.** Under the same Nov 4 setting
  node 5's subtree is about 1,290 DNN nodes and 44,000 laptop CPU s against
  node 1's 27,000 and 1.09e6: 21x fewer nodes, 25x less CPU (under full:
  11x and 15x). Level 4 tells the story: 35 of 47 children fathom here,
  16 of 47 under node 1. The two converged shortfalls differ by 89 units,
  0.69 of the 130-unit level of `../NOTES.md`, and 47^0.69 is 14x, so the
  gap between the two subtrees is about what the per-level model predicts.
  The tree model's "329-405-short class" (`../roadmap/tree_size.md`, 2 of
  16 sampled depth-3 nodes carrying 85-91% of the total) is therefore not
  one class: its weight has to be spread over the converged shortfall, and
  with node 1 as the representative it is overstated, by up to 2x on these
  two nodes. Which end of the class the other clear misses sit at (nodes
  2-5 and 7-8 converge 202-316 short) is what the sampling plan measures;
  the probe cannot say with two nodes.
- **For the roadmap's center (section 3):** the tail-class CPU per node it
  assumes is anchored on node 1, the worse of the two members by 25x. The
  21-30-day center is conservative on this evidence, not optimistic. The
  lever ratios it uses (full against the Nov 4 setting) are unchanged:
  3.1x in CPU on node 1, 1.8x here, 3.5x on node 10, with the level-5 per-node
  factor 2.7x, 4.0x and 4.3x.

## Equip phases: shorter cut phases where they cannot fathom (`QAP_DNN_CUT_EQUIP`, Nov 5 afternoon)

**Hypothesis.** In the matched probe's full arm, phases that started 100 or
more short consumed 42% of node 1's iterations (38,350 of 90,310) and
fathomed 2 of 30; no phase in 101 traces has gained more than 143. Such a
phase only equips the children with a pool, so capping it at a few Polyak
rounds should save most of its cost at no loss. Switch
`QAP_DNN_CUT_EQUIP=G:R` (`latest/bnb.C`, default off): a phase whose uncut
value is more than G below the cutoff runs at most R rounds. Arms
(`tail_probe.py`): `equip3` = full + `120:3`, `equip150` = full + `150:3`,
`equip150r2` = full + `150:2`. Levels 3-4 exact, level 5 the 87 hash-picked
nodes of full's parents, so the children of capped and uncapped parents can
be compared node for node (`runs/<arm>/`, `runs/matched_compare.txt`).

**Result: the hypothesis fails. The extra rounds are paid back by the
children.**

| arm | L4 fathomed / 47 | L4 iterations | L5 matched (87): fathomed, iterations | L5 children of capped parents: iterations vs full | est. subtree iterations |
|---|---|---|---|---|---|
| full | 18 | 53,900 | 83, 15,450 | — | 2,552,000 |
| equip3 (120:3) | 16 | 35,550 | 83, 15,920 | 84 nodes: 15,370 -> 15,840 (+3%) | 2,834,337 |
| equip150 (150:3) | 18 | 44,900 | 83, 16,670 | 63 nodes: 14,640 -> 15,860 (**+8%**) | 2,662,907 |
| equip150r2 (150:2) | 18 | 41,900 | 84, 19,010 | 63 nodes: 14,640 -> 18,200 (**+24%**) | 2,471,153 |

- The children of uncapped parents are identical to full's to the
  iteration (24 of 24 nodes in both 150 arms, 3 of 3 in equip3), so the
  runs are deterministic and the differences on the capped parents'
  children are the effect, not noise.
- **Level 4 saves what was predicted**: 9,000 iterations (17%) at 150:3,
  12,000 (22%) at 150:2, 18,350 (34%) at 120:3. At 150 no fathom is lost.
  At 120 two are: children 43:45 and 6:45 started 127 and 128 short and
  full's phases closed them in 7 and 8 rounds (48,814.1 and 48,814.2); the
  fathoming band reaches about 130-145 with the warm-Y rounds.
- **Level 5 pays it back.** A child of a 3-round phase costs 8% more than
  a child of the full phase, a child of a 2-round phase 24% more. The
  first solves are unchanged (7,170 -> 6,990 and 7,240 on the 63): the
  weaker pool shows in the children's own phases (per node 86 -> 102 and
  126), which is where a child that does not fathom in the first solve
  spends its iterations. With 1,334 level-5 nodes against 47 at level 4,
  +8% on 70% of level 5 is about +13,000 estimated iterations against
  the 9,000 saved: a wash at 3 rounds and a clear loss at 2. The subtree
  estimates agree in direction (equip150 +4%, SE 70%, so not significant
  on their own; the matched node figures are the evidence).
- The one case that would pay is a parent whose children all fathom in
  their first solves regardless of pool strength; the probe cannot tell
  those parents apart in advance, and the children of the two near-miss
  parents (43:45, 6:45) show the opposite shape: trivial children (48
  iterations each) under a parent that the full phase would have fathomed.

**Reading.** The phases at 100-250 short are not wasted; they are the
cheapest place to buy the pool that the 28x more numerous children spend.
The switch stays in the code, default off, for a later test on a class
where level 5 is cheap. The two remaining items from the Nov 5 list stand:
an earlier give-up for warm nodes that are still far short at 150
iterations (25% of the arm in 300-iteration first solves), and pruning
inactive pool cuts (pooled iterations about 15% dearer).

## Node 5 under the Nov 6 warm betas: `fast`, `wlowc`, `wtwo` (laptop, Nov 6, 12:10-12:58 PST)

Run 3 of the cloud's Nov 6 list (`../pipeline_ab/README.md`, "Next (Nov 6
evening)"): the matched probe on node 5 (`23:49,37:15,6:5`, 316 short) with
three arms at 2 slots each, binary 4f4198c, `runs_node5_w/`. `fast` is the
Nov 5 production candidate plus fp32 R and the 0.7 warm beta; `wlowc` runs
the children at 0.25 and cut rounds at 0.7; `wtwo` runs a child at 0.25
until the first refresh at or after iteration 60, then 0.7
(`DNN_ADMM_WARM_BETA=0.25:60:0.7`; the arm was added to `tail_probe.py` for
this run).

| arm | est. DNN nodes | est. iterations (± 1 SE) | est. CPU s (± 1 SE) | per level: fathomed/bounded |
|---|---|---|---|---|
| fast | 1,290 | 122,513 (± 50,146) | 15,219 (± 4,719) | L4 35/47; L5 35/36; L6 3/3 |
| wlowc | 1,290 | 91,027 (± 36,555) | 14,042 (± 4,989) | L4 35/47; L5 35/36; L6 3/3 |
| wtwo | 1,290 | 95,517 (± 37,028) | 14,389 (± 4,848) | L4 35/47; L5 35/36; L6 3/3 |

Matched nodes (iterations, first arm / second):

| arms | L4 (47 shared) | L5 (36 shared) | L6 (3 shared) |
|---|---|---|---|
| fast / wlowc | 35 / 35, 23,620 / 21,880 (1.08x) | 35 / 35, 4,330 / 3,290 (**1.32x**) | 3 / 3, 140 / 80 (1.75x) |
| fast / wtwo | 35 / 35, 23,620 / 21,770 (1.08x) | 35 / 35, 4,330 / 3,590 (1.21x) | 3 / 3, 140 / 80 (1.75x) |
| wlowc / wtwo | 21,880 / 21,770 (1.01x) | 3,290 / 3,590 (0.92x) | 80 / 80 (1.00x) |

- **No marginal fathom moves on this tail node.** The three trees are
  identical level for level (35 of 47, 35 of 36, 3 of 3), so the
  node-10 phase that went either way (`../pipeline_ab/README.md`, Nov 6
  runs) has no counterpart here, on 47 + 36 + 3 bounded nodes. One node
  of 86 on node 10, none of 86 here: the lost fathom is noise on the
  cutoff, not a pattern at the deep levels.
- **The warm beta pays at every level and most where the children are
  many**: `wlowc` 1.32x on the 36 matched level-5 nodes (first solves 82 to
  53 iterations a node), 1.75x at level 6 (47 to 27), 1.08x at level 4;
  1.35x on the subtree's iterations. `wtwo` sits between: 1.21x at level 5
  (61 a node, its switch to 0.7 at the refresh after iteration 60 lands in
  the middle of these solves), the same as `wlowc` at levels 4 and 6.
- **In CPU the subtree gains only 1.08x (`wlowc`) and 1.06x (`wtwo`)**: a
  53-iteration solve carries the same start-up, shipping and certificate
  cost as an 82-iteration one, so per iteration these children are dearer
  (level 5: 123 ms under `fast`, 151 under `wlowc`); the fixed cost per
  node is now a third of a tail child's time and the next thing to trim if
  the children get any faster.
- **For the production choice:** `wlowc` is the better arm on the tail
  (1.32x against 1.21x at level 5, no fathom lost), `wtwo` on node 10
  (where it kept the one marginal phase). Neither effect is large in CPU;
  `wlowc` as the candidate stands, and `wtwo` is the fallback if the
  sampling plan shows phases flipping against `wlowc` at any rate.
