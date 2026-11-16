# tho40: the grid cut phase in the search (bound_lab agent `tho40_search`, Nov 15 2026)

Question (`../COMMON.md`): tho40 has never run a cut phase. On complete subtrees, end to end, what is
the grid cut phase worth, and which `cut.zone`? **[m]** measured here (run named), **[e]** estimated.

## Findings

1. **The phase works in the search, and far beyond a near miss.** In the cut arms it fathomed nodes
   from a few units up to 3,193 units short (1.3% of the cutoff), and after the few misses the pooled
   children fathom almost at once. The tree all but ends two to three levels higher: 8x to 290x fewer
   nodes on every subtree run.
2. **The cost does not fall with the nodes: it moves into the phases.** Total DNN iterations (first
   solves + phases) under `cut.zone = 0.03` against no phase:
   - cold-rooted matched subtrees [m]: gate 1.6x fewer, t1 4.9x, t2 2.7x, t3 1.7x, t4 4.4x (t4's
     base from production [e]);
   - warm-rooted depth-3 subtrees against production [m, production's iterations e]: `38:32` (287
     nodes) 2.1x fewer iterations and 1.5x less CPU; `18:32` (1,490 nodes, the class that holds half
     of the A/B subtree's CPU) **2.9x fewer iterations, 2.7x less CPU**, with 74 nodes.
3. **Best zone: 0.03**, a phase at every give-up from depth 3 down. ZONE_SUMMARY
4. **Extrapolated [e]: about 2.4x less CPU on the A/B subtree (1.9-2.9x), and if the whole tree
   responds alike, about 550 c7a CPU-days (450-700) instead of 1,300.** One depth-2 subtree of 468,
   with one node standing for 54% of it: a first sizing, not an estimate to launch on.
5. **Surprises.**
   - Phases run long. Two thirds of a cut arm's iterations are phase iterations; a phase from 2,400
     short at depth 3 takes 3,341 iterations to fathom; t3's cold root ran to the 3,400 limit and
     ended 2.9 short; one pooled phase in `w18` spent 2,610 iterations to end 146 short, after which
     its 36 children fathomed in 21 iterations each.
   - First solves linger before the phase: near misses take 320-620 iterations, not 200, because
     `dnn.gu_zone` defers their give-up; the phase then needs 30-140.
   - On t2-t4 a cold root's first solve is 1,000-1,600 looser than the warm one, which moved their roots out
     of the 0.01 zone here; the warm-rooted runs are the ones to read for a run.
   - This box bounds a matched node in 0.70 of production's CPU.
6. **Sanity: PASS.** No arm found a solution under 240,516 and every subtree was fully fathomed; the
   three nodes checked on the path to the best known solution (depths 7, 9, 10, the last with an
   inherited pool) all branch with bounds at most 240,516.0000.

## How it was run

- One `latest/bnbdrv` process per subtree and arm (`run.sh`; two at a time, `queue.sh` over `jobs.txt`),
  `qaplib/tho40t.dat`, `latest/SYM/tho40t.sym`, `-u 240516`, `QAP_NODE_TRACE` and `QAP_DNN_TRACE` on.
  Each `runs/NAME_ARM/` has `cmd.txt` (the command), `log.txt`, `nodes.txt`, `dnn.txt`.
- **Arms** (`par/`): `base` = `experiments/estimates/tho40/tho40.par`, unchanged. `z003`, `z010`, `z030`,
  `zall` = base plus the block at the end of this file with `cut.zone` = 0.003, 0.01, 0.03, 0.06.
- **What the zone does** (`latest/bnb.C`, line 884). After a node's first solve ends unfathomed, the phase
  runs when `cutoff - first-solve bound <= cut.zone x |incumbent|`: 722, 2,405, 7,215 and 14,431 units
  short for the four zones.
  - The first solve does not depend on the zone. It gives up from iteration `dnn.giveup` = 200 when its
    extrapolated bound cannot reach the cutoff; `dnn.gu_zone` = 0.003 makes a give-up within 722 of the
    cutoff wait 100 more iterations, so near misses reach the phase after 320 to 400 iterations, not 200.
  - The phase (CP form) runs at most `cut.rounds x cut.iters + cut.freeze` = 3,400 iterations. It stops
    when it fathoms, or when after 600 iterations (`cut.stall` 3 x `cut.iters`) the gap is more than
    3 x the gain of the last 600 ("stalled").
  - An unfathomed phase hands its pool to the children (`pool.on`); a pooled child's first solve carries
    it (give-up rule unchanged), and its own phase, if it runs, starts from it.
  - `cut.damp1` does nothing in the CP form (the solver says so on stderr); it is in the block only
    because the brief's and sko42's have it.
- **Cold roots: `gate`, `t1`-`t4`.** Every arm starts the subtree's root cold from `-F`, so the arms are
  matched node for node. But a cold root stops at 200 iterations far looser than the same node in a
  run, which starts from its parent (t2: 3,193 short cold, 1,931 in production; t3: 3,462 and 2,427;
  t4: 4,274 and 2,675). So on t2-t4 the root is outside `z010` here, although in a run t2's would be
  inside it and t3's 22 units outside: read `z010` on t2-t4 as "no phase at the root, phases on its children".
- **Warm roots: `w18`, `w38`.** The A/B root's state file (`runs/wroot_ARM/root.state`, `-d 2 -o`), cut
  to one depth-3 child (`cutstate.py`, as `path_check/dive.py --warm`), resumed with `-r`: the node
  starts from its parent's shipped state as in a run. The first solves reproduce the production
  trace's bounds to the digit (A/B root 231,879.8596; `18:32` 4,474.4 short; `38:32` 2,426.7 short).
- **Reference: the production run** of the subtree `-F 26:10,8:0` (`infra/phase1/runs/tho40ab-32/`,
  22,650 nodes, 78,801 c7a CPU s, the base settings, no cut phase). The base arm here gave its node
  counts exactly on t1, t2, t3 (106, 211, 287) and the same outcome at every node both visited
  (`calib.py`). This box's CPU per matched node is 0.70 of production's, and the base arm takes 25.6
  iterations per node (t1-t3, cold roots aside); production's iterations, which its trace does not
  have, are estimated as nodes x 25.6 [e].
- **The box was shared** (load 5-7 on 4 cores; my two processes got about half a core each).
  **Iterations and nodes are the metrics.** CPU s are the solver's own CPU time and indicative. A phase
  iteration costs about 1.25x a first-solve iteration (0.126 against 0.103 s at m = 36, t4 z030).

## Cold-rooted subtrees, arm by arm [m]

Production's (warm-root) figures for the same nodes: gate 36 nodes below the root, 407 short; t1
`26:10,8:0,33:32,7:8` 106 nodes, 796 short; t2 `26:10,8:0,2:32,37:39` 211, 1,931; t3 `26:10,8:0,38:32`
(depth 3) 287, 2,427; t4 `26:10,8:0,35:32,18:39` 688, 2,675. `python3 table.py`:

@@ table.py

- Every phase in these runs but two fathomed its node; the two are the cold roots of t3 and t4 under
  `z030` (next section).
- **t4 base was not rerun** (about 1,400 CPU s, 45 minutes of the shared box): its ratios use
  production's 688 nodes and 688 x 25.6 iterations [e]. gate, t1 and `z030` arms where the root fathoms
  are the same under any zone that covers the root, so they were run once.

## Warm-rooted depth-3 subtrees under `z030` against production [m, production iterations e]

`python3 wtable.py`:

@@ wtable.py

`w18`, child by child (`python3 wchildren.py w18_z030 w18_zall`; children that fathom at once in
production and here are left out; `zall` is the next section's arm):

@@ wchildren.py w18_z030 w18_zall

## The phases [m]

`python3 phases.py` (every phase of every cut arm, the path check aside), by the shortfall of the solve
before the phase:

```
@@ phases.py
```

- **Up to about 1,000 short an unpooled phase fathoms in about 0.3 iterations per unit of shortfall**
  (44 to 280 iterations), against 500 to 2,300 base iterations under the same node
  (`python3 payoff.py`: 4 to 11 base iterations saved per phase iteration).
- **Farther out it still fathoms, but the price per unit rises**: 0.4-0.6 at 1,140-1,540 short
  (depths 4-5), 0.9 at 1,908 (depth 5: 1,717 iterations), 0.65 at 3,193 (t2's cold root, depth 4:
  2,070 iterations), 1.4 at 2,427 on a depth-3 node (`w38`: 3,341 iterations, 59 before the 3,400
  limit).
- **A phase that starts from an inherited pool is slower per unit**: the pool has already taken the
  easy part. In `w18`: 214 short 97 iterations, 382: 227, 493: 387, 517: 500, 567: 630, 648: 1,290,
  849: 1,900; and 1,037 short: stalled after 2,610 iterations, 146 short, after which its 36 children
  fathomed in 21 iterations each.
- **Misses at the roots.** t4's cold root (4,274 short): stalled after 1,810 iterations, 832 short.
  `w18`'s warm root (4,474): stalled after 920, 2,421 short, still gaining 0.8 an iteration. t3's cold root (3,462,
  depth 3): ran to the 3,400 limit and ended **2.9 short**; its 37 pooled children then all fathomed
  in 19 iterations each. That one is the run where `z030` lost to `z010` (4,307 against 3,401
  iterations): a limit 100 iterations longer, or stopping at 1,500 and leaving the rest to the pooled
  children, would both have been cheaper.
- **After a miss the pool carries.** Children of the three missed roots: t3 37 of 37 and t4 34 of 36
  fathom in their pooled first solve, the other 2 in short phases; in `w18` the pooled first solve
  takes about 1,370 off a child's production shortfall (1,891 to 517, 2,214 to 849, 1,006 and 953 to
  fathomed).

## The zone

- **`z003` (722 short) leaves most of it on the table.** ZONE003_TEXT
- **`z030` against `z010`** on the cold-rooted subtrees: t2 2,270 against 3,629 iterations, t4 3,966
  against 8,240, t3 4,307 against 3,401 (the 2.9-short miss above; warm, `w38` fathoms at the root in
  3,541). In t4 the eight phases `z010` runs on unpooled depth-5 children cost 4,042 iterations; under
  `z030` the root's one stalled phase (1,810) lets 34 of the 36 children fathom at once.
- **Where a phase at depth 3 does not pay.** A depth-3 node 1,800-2,400 short whose production subtree
  is 38-146 nodes (1,000-3,700 iterations) would get a phase of some 2,500-3,400 iterations at 1.25x
  the cost per iteration: break-even or a loss. `w38` (287 nodes, 7,300 iterations) wins 2x. Such nodes
  are 7% of the A/B subtree's CPU, so this trims the total by a few percent; it argues for a phase
  limit that depends on what is left, not for a narrower zone, which would also drop the phases that
  equip the children of the far nodes.
- **`zall` (a phase at the depth-2 root too).** ZALL_TEXT

**Best zone measured: 0.03** (every give-up from depth 3 down in this subtree; the largest depth-3
shortfall is 5,039).

## What it means for the tree [e]

**The A/B subtree** (22,650 nodes, 78,801 CPU s), by the shortfall class of its 38 depth-3 nodes
(production CPU shares from the trace, `prodstats.py`):

| depth-3 class (short) | nodes of that class | share of the subtree's CPU | CPU factor under `z030` | from |
|---|---|---|---|---|
| fathomed at once, or <= 1,500 | 12 | 1.2% | about 1.5x | [e] gate, t1 (deeper nodes of that shortfall) |
| 1,500-2,500 | 11 | 7.4% | 1.0-1.6x, take 1.2x | [m] `w38`, the largest of them, 1.6x; the small ones break even (above) |
| 2,500-3,500 | 5 | 12.3% | about 2x | [e] between `w38` and `w18`; not run (`w30` was queued) |
| 3,500-5,000 | 9 | 53.7% | 2.7x | [m] `w18` (4,474 short, 1,490 nodes): 2.70x CPU, 2.93x iterations |
| above 5,000 | 1 (`2:32`, 6,281 nodes) | 25.3% | as `w18` | [e] not run; its depth-4 children are t2/t4-like (2.7x and 4.4x in iterations, cold) |

Weighted by CPU share (1 / sum of share / factor): **about 2.4x less CPU on the A/B subtree under `z030`**,
1.9x if every class is at the low end of what the runs allow (1.0, 1.0, 1.5, 2.2, 2.0) and 2.9x at the
high end (2, 1.6, 2.5, 3.2, 3.5). In iterations about 2.6x; in nodes 20x or more (`w18`: 1,490 to 74).
Measured in it: two warm depth-3 subtrees, of which **one node, `w18`, stands for 54% of the CPU** and,
by assumption, for the 25% above it. The four cold-rooted subtrees (1.6x to 4.9x in iterations, depth-4
roots mostly) say the deeper nodes do at least as well.

`python3 replay.py` is a cruder model on the same trace: every branched node within the zone gets a
phase that fathoms from at most 3,000 short at the iteration costs measured at depths 4-5 (low for
depth 3), and after a miss the children either all fathom (`all`) or get no help (`none`). It gives
`z003` x1.35 less CPU with 2.5x fewer nodes (no pools arise, so this one is firm), `z010` x1.9, and
`z030` x1.9 (`none`) to x5.7 (`all`). The class table's 2.4x sits inside that bracket: the pool's
carry is partial, and depth-3 phases cost more than the model's.

**The whole tree** (3.4e7 nodes, about 1,300 c7a CPU-days by Knuth's estimator, itself "likely low"):
- **If it responds as the A/B subtree: about 550 c7a CPU-days (450-700), from 1,300**, and of the order
  of 1e6 to 2e6 nodes instead of 3.4e7. That is an extrapolation from one depth-2 subtree of 468.
- **Why it may differ.** The A/B subtree (22,650 nodes) is smaller than the average depth-2 subtree
  (73,000), and its root is 8,634 short where the path's depth-2 node is 9,075-10,052. In a larger
  subtree more depth-3 nodes are beyond 5,000, the class not measured here; the phase then does its
  work one level down, on smaller matrices. tho40-b's merged tree at 25 h had the same shape by depth
  as the A/B subtree, shifted about half a level deeper (`experiments/tho40b_analysis/README.md`), and
  its better-bound curve was within 5 points of the A/B subtree's.
- **What the factor is made of.** The tree all but disappears below depth 4-5 (in `w18` 65 of 74 nodes
  fathom in a pooled first solve of 20-140 iterations), and the cost moves into the phases: 66% of
  `w18`'s iterations are phase iterations, at 1.25x the cost each. Bounding fewer, harder nodes for
  longer is also a different load for the coordinator: units of minutes, not seconds.
- **Not in it:** any tuning of the phase for tho40 (next section); the incumbent is the best known
  value, as in the 1,300.

## Sanity

- **Outcomes.** Every run ended with no `New objective value` line and its subtree fully fathomed
  (`summarize.py`, column `newobj` = 0; each `log.txt` ends with the RESULTS block). The base arm and
  production fathom the same subtrees node by node without a better solution, so every node a cut arm
  fathomed has an optimum above the cutoff: no cut-arm fathom contradicts the plain search.
- **The path to the best known solution** (`python3 pathverify.py`; cold path of
  `experiments/tho40b_analysis/path_check/cold.tsv`, incumbent 240,518 as there, arm `z030`). A node on
  the path has optimum 240,516, so no certified value on it may pass 240,516 and it must branch:

```
@@ pathverify.py
```

  The first solves are `cold.tsv`'s to the digit. At depth 7 the node is 3,351 short and the phase stalls
  1,312 below the optimum (t2's root, 3,193 short and with no optimal solution under it, fathomed). At
  depth 9 the phase climbs to 240,494.30 in 3,400 iterations and stops at the limit; its path child,
  with the inherited pool, reaches 240,516.0000 and branches; the 30 siblings fathom. Not covered:
  depths 1-6 and 8, the warm path, pooled nodes above depth 10 (path8, path10 and the children of
  path7 were dropped for time). `../tho40_nodes/` checks the CP solve itself on the path and by
  brute force on small instances.
- `cmd.txt` of `t3_base` and `t3_z010` lack the exit line: `run.sh` was edited while they ran (the
  logs are complete; noted in the files). `run.sh` is now replaced only by `mv`.

## Recommended settings for tho40

`experiments/estimates/tho40/tho40.par` plus (this is `par/z030.par`, less `cut.damp1`, which the CP
form ignores):

```
pool.on = 1
cut.grid = 5:8
cut.side = 1
cut.cp = 100:1000
cut.warmy = 1
cut.rounds = 12
cut.iters = 200
cut.freeze = 1000
cut.zone = 0.03
```

`tho40t.dat`'s second matrix is the 5 x 8 grid (`cut.side = 1`), as the gate test used. Not changed
here and worth an A/B before a long run (next section): `dnn.gu_zone`, the phase's limit and stop rule,
`cut.zone = 0.06`.

## What I would do next

1. **The phase's stop rule, on `w18` warm (13,035 iterations to beat) and t3/t4.** The evidence above
   says a phase should stop and equip its children once it will clearly not fathom soon: t4's root
   stalled at 1,810 and the 36 pooled children cost 1,956 in all, where the phases that insisted cost
   1,900-3,400 each. Candidates: `cut.giveup = T:G`, a shorter `cut.cp_iters`, `cut.stall = 2`. The DNN
   traces written after 23:00 (`w*`, `path7`, `t1_z003`, `t4_z030g`) carry each phase's running best
   every 10 iterations (`QAP_DNN_TRACE_HIST`), so stop rules can be replayed offline first.
2. **`dnn.gu_zone = 0` with the phase** (`par/z030g.par`): GZ_TEXT
3. **The unmeasured classes before any whole-tree number:** `30:32` warm (`w30`, queued, not run),
   the `2:32` node (25% of the A/B subtree) and one depth-2 subtree with a far root; then Knuth dives
   (`bnbdrv -e`) under the cut arm, which is how the 1,300 was made.
4. **sko42's `dnn.*` keys with the phase** (`par/z030c.par`: beta 0.3, the two-phase schedule, refresh 40,
   give-up 300): prepared, not run.
5. **The pool through state files works** (the `w*` runs resume a parent's state with `-r`; `zall`'s
   children received the root's pool that way), so a coordinator run needs nothing new; its units
   become few and long.

## Files

| | |
|---|---|
| `par/*.par` | the arms (`z030c`, `z030g`: prepared, not run) |
| `run.sh`, `queue.sh`, `jobs.txt`, `queue.log`, `path_jobs.txt` | the runner, the job list in the order run, its log |
| `runs/NAME_ARM/` | `cmd.txt`, `log.txt`, `nodes.txt` (node trace), `dnn.txt` (DNN trace; with the phase's history from 23:00 on), `err.txt` |
| `summarize.py`, `table.py`, `wtable.py`, `wchildren.py`, `phases.py`, `payoff.py` | the tables above |
| `calib.py`, `pick.py`, `prodstats.py`, `replay.py` | the production trace: matched nodes, subtree choice, CPU by class, the replay model |
| `cutstate.py`, `pathverify.py`, `mkreadme.py`, `README.tmpl.md` | the warm-root state cut, the path verdict, this file's assembly |
