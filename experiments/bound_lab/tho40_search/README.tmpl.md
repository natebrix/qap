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
3. **Best zone: 0.03 or wider**, a phase at every give-up from depth 3 down. `cut.zone = 0.003` gets 1.4-1.9x
   where 0.03 gets 2.7-4.9x (t1, t2); 0.01 is behind 0.03 on t2 and t4 and ahead on t3. A zone that
   also covers depths 1-2 (0.06) ran the same 74 nodes on `18:32` in 11% fewer iterations: the first
   thing to A/B, and what I would launch with.
4. **Extrapolated [e]: about 2.4x less CPU on the A/B subtree (1.9-2.9x), and if the whole tree
   responds alike, about 550 c7a CPU-days (450-700) instead of 1,300.** One depth-2 subtree of 468,
   with one node standing for 54% of it: a first sizing, not an estimate to launch on.
5. **Surprises.**
   - Phases run long. Half to nine tenths of a `z030` arm's iterations are phase iterations; a phase from 2,400
     short at depth 3 takes 3,341 iterations to fathom; t3's cold root ran to the 3,400 limit and
     ended 2.9 short; one pooled phase in `w18` spent 2,610 iterations to end 146 short, after which
     its 36 children fathomed in 21 iterations each.
   - First solves linger before the phase: near misses take 320-620 iterations, not 200, mostly
     because `dnn.gu_zone` defers their give-up; the phase then needs 30-140.
   - On t2-t4 a cold root's first solve is 1,000-1,600 looser than the warm one, which moved their
     roots out of the 0.01 zone here; the warm-rooted runs are the ones to read for a run.
   - This box bounds a matched node in 0.70 of production's CPU.
6. **Sanity: PASS.** No arm found a solution under 240,516 and every subtree was fully fathomed; the
   three nodes checked on the path to the best known solution (depths 7, 9, 10, the last with an
   inherited pool) all branch with bounds at most 240,516.0000.

7. **Follow-up, Nov 16 ("The stop rule" below): another 2.1-2.3x against `zall`, 2.4-2.7x against `z030`.** With a hand-over stop rule
   (a new key, patch in `stop/`), `dnn.gu_zone = 0` and `dnn.giveup_warm = 100` under `cut.zone = 0.06`,
   `18:32` takes 5,523 iterations (13,035 above) and `38:32` 1,312 (3,541); the extrapolation becomes
   about 5.8x on the A/B subtree (4.0-6.8x) and **about 225 c7a CPU-days (190-330)** for the tree if
   it responds alike. Items 3 and 4 above are the first study's numbers, kept as written.

8. **Third study, Nov 16 ("With the faster phase" below): with cpfast's `wc=2.88` as well.** `18:32`
   4,322 iterations (38 nodes), `38:32` 410, `30:32` 1,687, and `2:32`, the last unmeasured class,
   20,992 against production's 160,800 [e]. About 7.8x on the A/B subtree (6.2-9.4x) and **about 165
   c7a CPU-days (140-210)** for the tree if it responds alike. One block, one patch
   (`stop/combined.patch`).

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
- **Warm roots: `w18`, `w38`.** The A/B root's state file (`runs/wroot_ARM/root.state`, `-d 2 -o`, cold from `-F` as production's root was), cut
  to one depth-3 child (`cutstate.py`, as `path_check/dive.py --warm`), resumed with `-r`: the node
  starts from its parent's shipped state as in a run. The first solves reproduce the production
  trace's bounds to the digit (A/B root 231,879.8596; `18:32` 4,474.4 short; `38:32` 2,426.7 short).
- **Reference: the production run** of the subtree `-F 26:10,8:0` (`infra/phase1/runs/tho40ab-32/`,
  22,650 nodes, 78,801 c7a CPU s, the base settings, no cut phase). The base arm here gave its node
  counts exactly on t1, t2, t3 (106, 211, 287) and the same outcome at every node both visited
  (`calib.py`). This box's CPU per matched node is 0.70 of production's, and the base arm takes 25.6
  iterations per node (t1-t3, cold roots aside); production's iterations, which its trace does not
  have, are estimated as nodes x 25.6 [e] (counting each node at the base arm's mean for its depth and
  outcome instead gives the same within 8%; t3's measured 7,281 is 287 x 25.4).
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
  `w18`'s warm root (4,474): stalled after 920, 2,421 short, still gaining 0.8 an iteration. t3's cold
  root (3,462, depth 3): ran to the 3,400 limit and ended **2.9 short**; its 37 pooled children then all fathomed
  in 19 iterations each. That one is the run where `z030` lost to `z010` (4,307 against 3,401
  iterations): a limit 100 iterations longer, or stopping at 1,500 and leaving the rest to the pooled
  children, would both have been cheaper.
- **After a miss the pool carries.** Children of the three missed roots: t3 37 of 37 and t4 34 of 36
  fathom in their pooled first solve, the other 2 in short phases; in `w18` the pooled first solve
  takes about 1,370 off a child's production shortfall (1,891 to 517, 2,214 to 849, 1,006 and 953 to
  fathomed).

## The zone

- **`z003` (722 short) leaves most of it on the table.** t1: 1,396 iterations (37 nodes) against 540 under
  `z010` or `z030` and 2,626 base. t2: 4,331 (72 nodes) against 3,629 (`z010`), 2,270 (`z030`) and
  6,187 base. All its phases fathom, in 73 to 140 iterations, but nodes 800-3,000 short still branch.
  On the production trace `replay.py` puts it at x1.35 less CPU for the A/B subtree [e].
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
- **`zall` (a phase at the depth-2 root too).** `zall` differs from `z030` only above depth 3 here.
  The A/B root (8,634 short) ran a phase that stalled after 760 iterations 5,100 short, and its pool
  went to the 38 depth-3 children. `18:32` then started 3,163 short instead of 4,474, its own phase
  stalled after 650 iterations instead of 920 at about the same place (2,298 against 2,421), and the
  whole subtree took the same 74 nodes in **11,565 iterations against 13,035** (`python3 matched.py
  w18_z030 w18_zall`: 0.90 on the 37 children, every hard child a little cheaper). The root's 760
  iterations are 20 per child. One subtree, and a 3.0x CPU factor against production where `z030` has
  2.7x: enough to prefer it, not enough to move the estimate.

**Best zone: 0.03 or wider.** 0.03 is every give-up from depth 3 down in this subtree (the largest
depth-3 shortfall is 5,039) and is what six subtrees ran; 0.06 adds the phases at depths 1-2 (a few
hundred nodes in the whole tree, at most 3,400 iterations each) and was 11% cheaper on the one subtree
that can show the difference.

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
by assumption, for the 25% above it. The five cold-rooted subtrees (1.6x to 4.9x in iterations, depth-4
roots but one) say the deeper nodes do at least as well.

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
  work one level down, on smaller matrices. tho40-b's merged tree at 25 h has its CPU
  deeper than the A/B subtree's (depth 6 and below: 59% against 40%), while its better-bound curve
  was within 5 points of the A/B subtree's (`experiments/tho40b_analysis/README.md`).
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
  depths 1-6 and 8, the warm path, any pooled node but the depth-10 one (path8, path10 and the
  children of path7 were dropped for time). `../tho40_nodes/` checks the CP solve itself on the path and by
  brute force on small instances.
- `cmd.txt` of `t3_base` and `t3_z010` lack the exit line: `run.sh` was edited while they ran (the
  logs are complete; noted in the files). `run.sh` is now replaced only by `mv`.

## Recommended settings for tho40

(The first study's. "Recommended settings, restated" in the stop-rule section adds three keys.)

`experiments/estimates/tho40/tho40.par` plus (this is `par/zall.par`, less `cut.damp1`, which the CP
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
cut.zone = 0.06
```

- `cut.zone = 0.06` is `par/zall.par`: a phase at every give-up within 14,431 units, which is depth 2
  and below and the depth-1 node on the path (14,269 short). `cut.zone = 0.03` (`par/z030.par`) is the
  conservative choice: identical from depth 3 down wherever no pool comes from above, and the setting behind every factor quoted here
  but `w18 zall`'s.
- `tho40t.dat`'s second matrix is the 5 x 8 grid (`cut.side = 1`), as the gate test used.
- Left as they are in `tho40.par`, each worth an A/B before a long run (next section): `dnn.gu_zone`
  (0 looks 7% better), the phase's limit and stop rule.

## What I would do next

(Written after the first study. Items 1 and 2 were then done: "The stop rule" below.)

1. **The phase's stop rule, on `w18` warm (13,035 iterations to beat) and t3/t4.** The evidence above
   says a phase should stop and equip its children once it will clearly not fathom soon: t4's root
   stalled at 1,810 and the 36 pooled children cost 1,956 in all, where the phases that insisted cost
   1,900-3,400 each. Candidates: `cut.giveup = T:G`, a shorter `cut.cp_iters`, `cut.stall = 2`. The DNN
   traces of the runs started after 22:48 (`t4_z010`, `w*`, `path7`, `t1_z003`, `t2_z003`, `t4_z030g`) carry each phase's running best
   every 10 iterations (`QAP_DNN_TRACE_HIST`), so stop rules can be replayed offline first.
2. **`dnn.gu_zone = 0` with the phase** (`par/z030g.par`): run once [m], on t4: the same 37 nodes and
   the same phases, 3,697 iterations against `z030`'s 3,966 (7% fewer; the two pooled near misses hand
   over after 210-260 iterations instead of 350-400). In `w18` seven children waited 300-350 iterations
   before their phase. Worth adopting after one warm-rooted A/B; a lower `dnn.giveup` for pooled
   children is the same idea taken further.
3. **The unmeasured classes before any whole-tree number:** `30:32` warm (`w30`, queued, not run),
   the `2:32` node (25% of the A/B subtree) and one depth-2 subtree with a far root; then Knuth dives
   (`bnbdrv -e`) under the cut arm, which is how the 1,300 was made.
4. **sko42's `dnn.*` keys with the phase** (`par/z030c.par`: beta 0.3, the two-phase schedule, refresh 40,
   give-up 300): prepared, not run.
5. **The pool through state files works** (the `w*` runs resume a parent's state with `-r`; `zall`'s
   children received the root's pool that way), so a coordinator run should need nothing new (none
   was tried here); its units become few and long.

## The stop rule (follow-up, Nov 16)

Question: when should a phase stop and hand over to its pooled children? The first study's stop rule
is the stall rule: after 600 iterations, stop when the gap exceeds 3 x the gain of the last 600. Its
runs showed it is too patient near the cutoff (2,610 iterations to end 146 short; 3,341 to fathom from
2,427) and, as it turned out, too impatient far from it. Result [m]: with a
hand-over rule (a new key, patch below), `dnn.gu_zone = 0` and `dnn.giveup_warm = 100`, all under
`cut.zone = 0.06`, the warm `18:32` subtree takes **5,523 iterations against 11,565 (`zall`) and
13,035 (`z030`)**, and `38:32` **1,312 against 3,012 and 3,541**: 2.1x and 2.3x fewer than `zall`,
2.4x and 2.7x fewer than `z030`, the same trees (74 nodes and 1 node). `30:32`, not run before, takes 2,005 iterations and 38 nodes where
production took 571 nodes.

### Offline: candidate rules replayed on the recorded phases [m histories, e hand-over cost]

`stop/hist.py` reads the phase histories (`QAP_DNN_TRACE_HIST`, the running best every 10 iterations):
37 phases recorded under the stall rule, the path check aside. `stop/replay_rules.py` replays a rule on
each as the solver's `stop()` would see it. A phase a rule stops before its recorded fathom is a lost
fathom, and costs what the node's children then cost.

**The hand-over cost** is measured where a phase missed and its children ran complete, and modelled as
`C(g) = m x (20 + 2 g + 3 g^2)` iterations, m the node's children, g the gap left in 0.1% of the cutoff
(240.5 units): 20 iterations per child when the parent ends at the cutoff (19.1, 20.5, 20.6 measured),
more as the gap grows and some children need a phase of their own.

| missed phase | gap left | children measured | model |
|---|---|---|---|
| t3 `z030` root (first study, no history) | 2.9 | 707 | 741 |
| `w18` child `2:7`, `z030` / `zall` | 146 / 169 | 741 / 739 | 803 / 824 |
| t4 `z030g` root | 832 | 1,687 | 2,262 |
| `w18` root, `zall` / `z030` | 2,298 / 2,421 | 10,715 / 11,915 | 11,582 / 12,730 |

Three nodes; the quadratic term rests on `w18`'s root and t4's. Cost of a rule = phase iterations + C at
its stops (`python3 stop/replay_rules.py`, in `stop/replay_rules.out`; "near" = phases that start within
2,000 of the cutoff, "far" = the 5 that start beyond; x = recorded cost / the rule's):

| rule | near (32 phases): phase it + hand-over, lost fathoms | far (5): the same | all |
|---|---|---|---|
| recorded: `cut.stall 3`, `cut.iters 200` | 21,830 + 1,628 = 23,458, 0 | 7,481 + 80,207 | x1.00 |
| window 400 (`cut.stall 2`) | 21,027 (x1.12), 0 | x0.84 | x0.88 |
| window 300 (`cut.iters 100`) | 20,656 (x1.14), 0 | x0.73 | x0.79 |
| window 200 (`cut.iters 100`, `cut.stall 2`) | 19,726 (x1.19), 2 | x0.62 | x0.69 |
| limit `cut.cp_iters` 1,500 / 1,000 / 600 | x1.03, 4 / x1.09, 6 / x1.10, 8 | x1.01 / x1.00 / x0.89 | x1.01 / x1.02 / x0.93 |
| `cut.giveup 600:200` | 20,253 (x1.16), 5 | x0.89 | x0.94 |
| `cut.giveup 400:400,800:150` | x1.09, 4 | x0.76 | x0.81 |
| rate rule: stop when gap / rate > 20 m (window 200) | 19,790 (x1.19), 1 | x0.64 | x0.71 |
| the same with 40 m | 21,027 (x1.12), 0 | x0.80 | x0.85 |
| **hand-over rule: stop when 1.5 x gap / rate > C(gap)** (window 200) | 13,818 + 5,466 = 19,284 (x1.22), 3 | x1.01, and 3 recorded stalls it would let run on | x1.05 |

- **Near the cutoff every impatient rule saves about the same, 12-22% of the phase iterations**: the
  late phases are few (in each `w18` run two phases of eight carry the saving), and a lost fathom costs
  20-25 iterations a child.
- **Every rule the solver has loses it back on the far phases.** A shorter window, a limit, a give-up
  pair or the plain rate rule all stop a far phase earlier, and by C its children then cost far more
  than the phase saved: by the model's slope a unit of gap closed at 2,300 short is worth about 9 child
  iterations, and `w18`'s root was closing it at 1.2 phase iterations a unit when it stalled.
- **The replay cannot value letting a phase run on** (its history ends where the stall rule stopped it),
  so "x1.01 far" for the hand-over rule is a floor. That part had to be run.
- Sensitivity of the hand-over rule (same table, `replay_rules.out`): the factor 1.5 to 2, the window
  100 to 300, N0 20 to 30 and the quadratic coefficient 1.5 to 6 all give x1.17 to x1.23 near (factor 1:
  x1.13).

### The key the solver lacks: `cut.handover`

No existing key is impatient near the cutoff and patient far from it: the stall rule compares the gap
with the gain of a fixed window (factor 3, not a setting), `cut.giveup` and `cut.cp_iters` only stop.
`stop/handover.patch` (37 lines over `latest/settings.C`, `dnncut.h`, `dnncut.C`, `bnb.C`; built in my
worktree `/home/user/wt-tho40search`, the main checkout untouched) adds

```
cut.handover = N0:A:B[:PHI[:W]]      (run: 20:2:3, PHI 1.5, W 200)
```

In `CutProjector::stop()`, in place of the stall rule: from iteration W, with `rate` the gain of the
running best over the last W iterations per iteration and g the gap in 0.1% of |cutoff|, the phase
stops when `PHI x gap / rate > m x (N0 + A g + B g^2)`: when finishing at its present rate (PHI for
its decay) would cost more iterations than the children would if it handed over now. Unset, the solver
is the main checkout's: the patched binary reproduced the unpatched root phase's history to the digit
over its 690 recorded iterations. The certified values do not depend on where a phase stops.

**Entering the phase early** (the per-node lab's finding, `../tho40_nodes/README.md`) has a key
already: `dnn.giveup_warm = N` lets a warm-started node give up from iteration N instead of
`dnn.giveup` = 200 (the give-up rule needs 110 iterations of history, so 100 means 110), and
`dnn.gu_zone = 0` removes the 100-iteration confirmation near the cutoff. Both are
in the A/B: on `38:32` and on the 8 give-up nodes of `18:32` the first solve ends at 110 iterations
instead of 200-220 and the phase then needs 0 to 25 iterations more, about 90 saved a node, as the
per-node lab found.

### In the search: warm-rooted subtrees [m]

Arms (all `zall`, `cut.zone = 0.06`): `zall` the first study's; `zallhg` = `zall` + `cut.handover = 20:2:3`
+ `dnn.gu_zone = 0`; `zallhge` = `zallhg` + `dnn.giveup_warm = 100`. One `bnbdrv` at a time
(`stop/run2.sh`, `stop/seq.sh` over `stop/jobs.txt`, log `stop/seq.log`), the box at load 6-8.
Each arm's subtree is resumed from that arm's A/B root (`wroot_ARM`; `zallhg` and `zallhge` from
`wroot_zallh`, whose solve their extra keys do not touch). `python3 stop/ab.py 2`:

@@ stop/ab.py 2

`18:32` child by child (`python3 wchildren.py w18_zall w18_zallhg w18_zallhge`):

@@ wchildren.py w18_zall w18_zallhg w18_zallhge

- **The same trees in about half the iterations**: `18:32` 74 nodes in every arm, `38:32` fathomed at its
  root in every arm. Every run ended with its subtree fathomed and no solution under 240,516.
- **What the rule did.**
  - *Far, at the depth-2 root:* its phase ran 3,030 iterations instead of 760 and ended 4,050 short
    instead of 5,100 (the rule stopped it when the rate had fallen to 0.17 an iteration). Those 2,270
    extra iterations are paid once and every depth-3 child starts closer: `18:32` 2,468 short
    instead of 3,163, `38:32` 984 instead of 1,380, `30:32` 1,807 where production had 3,492.
  - *Near:* `18:32`'s child `2:7` handed over after 490-500 iterations, 261 short, where the stall
    rule ran 2,200 to end 169 short; its 36 children then fathomed in 21 iterations each (767; the
    model said 925). `30:32`'s root handed over after 550 iterations, 1,032 short; 36 of its 37
    children fathomed in their first solve and one in a 210-iteration phase.
  - *No fathom was lost to it that the record shows:* the phases it did not stop fathomed in 74 to
    1,202 iterations.
- **Where `18:32`'s gain comes from** (`runs/w18r_zallhg`: the rule and `dnn.gu_zone = 0` below the
  *first study's* root, whose phase had stalled at 760): 11,565 (`zall`) to **7,813** with the rule
  and `gu_zone` alone (1.48x; here the rule let `18:32`'s own phase run 1,320 iterations instead of
  650, to 2,011 short instead of 2,298, and stopped `2:7` at 470), to 6,231 once the depth-2 root is
  under the rule too (1.25x more), to 5,523 with `dnn.giveup_warm = 100` (1.13x more).
- **The two `dnn.*` keys** show in the first-solve iterations of `18:32`: 4,429 (`zall`), 3,500 with
  `dnn.gu_zone = 0` (give-ups at 200-220 instead of 300-350; part of this is also the better pool),
  2,717 with `dnn.giveup_warm = 100` as well, while the phases stay at 2,731 and 2,806.
- **The hand-over model, checked on the new runs.** It overstates what children cost under the new
  settings: `18:32`'s root left 1,888 and its children cost 4,793 (model 8,160); `30:32`'s left 1,032
  and they cost 1,345 (model 3,100). The children's own phases are cheaper now. Refit, the quadratic
  coefficient would be 0.4 to 1.5, not 3, which makes far phases hand over sooner: the first thing to
  A/B next. The constants run are the first study's; nothing was refit to these runs.
- **Not run:** a second rule (the best of the solver's own keys, window 400, is x1.12 near and loses far
  in the replay); the rule without the two `dnn.*` keys; `zall` alone on `30:32`.

### Sanity under the final arm

Every `zallh*` run ended with its subtree fathomed and no `New objective value`. The stop rule changes
only where a phase stops, never a certified value; the path check was rerun all the same with the
patched solver and `par/zallhge.par` (`python3 pathverify.py zallhge`, incumbent 240,518):

```
@@ pathverify.py zallhge
```

All three branch, none above 240,516, and the first solves at depths 7 and 9 are `cold.tsv`'s. The
depth-10 path child, entered after 110 iterations with the inherited pool, approaches its optimum from
below (0.014 under it at 1,400 iterations, 240,516.0000 to four places when the rule stops it at 1,580).
On a path node the gap to the cutoff shrinks about as fast as the rate, so the rule lets the phase run:
one long phase per path node, as under the stall rule (3,400 and 600 iterations there).

### Recommended settings, restated

`experiments/estimates/tho40/tho40.par` with `dnn.gu_zone = 0` in place of its 0.003, plus
(this is `par/zallhge.par`, less `cut.damp1`):

```
pool.on = 1
cut.grid = 5:8
cut.side = 1
cut.cp = 100:1000
cut.warmy = 1
cut.rounds = 12
cut.iters = 200
cut.freeze = 1000
cut.zone = 0.06
cut.handover = 20:2:3
dnn.giveup_warm = 100
```

`cut.handover` needs `stop/handover.patch` in the solver (and its line in RUNNING.md's table and the
settings goldens, which the patch does not touch). Without the patch the block is the first study's
with the two `dnn.*` changes, which the main checkout already reads: not run as such here; by the
first-solve iterations above, worth 10-15% on `zall`.

### The whole tree, restated [e]

`python3 stop/tree.py`: the first study's class table with the factors measured under `zallhge`, and the
A/B root's own longer phase added:

@@ stop/tree.py

- **About 225 c7a CPU-days (190-330) if the whole tree responds as the A/B subtree**, where the first
  study had 550 (450-700) and the base settings 1,300. In nodes, of the order of 1e6.
- **What is measured:** three warm depth-3 subtrees in the classes that hold 73% of the A/B
  subtree's CPU, each 4x to 7x cheaper than production in CPU and 5.6x to 7.3x in iterations.
- **What is not:** the `2:32` node (25%), taken as `18:32` with half of it as the low end; the other
  35 depth-3 nodes; any depth-2 node but this one, and the depth-2 phase is now where much of the gain
  starts; depth 1. The range is the class factors' spread only. The earlier caveats stand: one
  depth-2 subtree of 468, smaller than average, CPU from a shared box (iterations weighted 1.25 for
  phases give 6.1x on `18:32` where CPU gives 6.4x).
- **Before a launch:** the same arm on a second depth-2 subtree with a far root, warm from depth 1,
  and Knuth dives under it.

## With the faster phase (third study, Nov 16): hand-over policy + `wc=2.88`

`../cpfast/README.md` found that the phase was running at the warm children's penalty and that
`dnn.beta_sched = w=0.36,wc=2.88` makes a far phase 2.4 to 4.6 times shorter. This section puts that
together with the second study's arm. Result [m], warm subtrees, one process
at a time: **`18:32` 4,322 iterations and 38 nodes (5,523 and 74 under the second study's arm, 13,035
under `z030`), `38:32` 410 (1,312; 3,541), `30:32` 1,687 (2,005)**, and `2:32`, never run before,
20,992 iterations and 288 nodes where production took 6,281 nodes. Sweeps and extrapolation on top
add nothing in the subtrees. The path check passes and costs
no more than before. Extrapolated, with all four classes now measured: about 7.8x less CPU than
production on the A/B subtree and **about 165 c7a CPU-days (140-210) for the tree if it responds
alike**.

### The solver: `stop/combined.patch`

My worktree `/home/user/wt-tho40search` with `../cpfast/cpfast.patch` and `stop/handover.patch`
applied (one context conflict, in `DNNCutParams`' defaults; `git diff` of the result is
`stop/combined.patch`, against 2f2f4f1b: 239 lines added over 7 files). Every new key defaults to today's behaviour.
`cut.handover` grew two optional fields here:

```
cut.handover = N0:A:B[:PHI[:W[:P[:KAPPA]]]]     run: 20:0:0.22:1.5:100:3:1.25
```

The children's cost is `C = m x (N0 + A g + B g^P)` (P = 2 before). From iteration max(W, 200), with
`rate` the gain over the last W iterations per iteration, the phase stops when
(i) `PHI x gap / rate > C` (finishing costs more than handing over) and, if KAPPA > 0,
(ii) `rate x dC/dgap < KAPPA` (one more phase iteration no longer saves KAPPA child iterations).
`20:2:3` still means what it meant in the second study.

### The refit [m hand-overs, e model]

`python3 stop/replay3.py` (`stop/replay3.out`): every hand-over whose children ran complete, mine and
cpfast's two `w18` runs at `wc=2.88`, in iterations per child, against the second study's model
(`quad3` = 20 + 2 g + 3 g^2), a quadratic with half the coefficient, and a cubic:

@@ stop/replay3.py --table

- **The cost of a hand-over is cubic in the gap, not quadratic**: 20 a child up to 340 left, 36 at
  1,030, 130-170 at 1,900-2,000. `20 + 0.22 g^3` is within 16% of every hand-over run under the second
  study's settings (`zallhg`, `zallhge`, `w18r`); `quad3` is 1.5 to 2.3 times too high between 800 and
  2,000 (the runs it was fitted to had 2,600-iteration child phases). The cubic is low by 30-45% on
  t4's and cpfast's rows, whose give-up children still spent 300-350 iterations in the first solve,
  and the rows under the first study's settings beyond 2,200 are above every model.
- **The `zallhw*` rows were run after the refit** and are a check on it: 33 a child against the
  cubic's 27 at 757 left, 97 against 72 at 1,486. Low by a fifth to a third: under `wc` the pool a
  stopped phase leaves is less converged. Not refit again.
- **With an honest C, condition (i) alone is too impatient far from the cutoff**: it compares finishing
  with handing over, and a far phase cannot finish. What a far phase buys is a better pool, so the
  stop there is the marginal one, (ii). The quadratic's overstatement had been standing in for it.
- **A phase at a new penalty starts flat** (cpfast: 30-50 iterations, 100-150 from a pool), so the rate
  is taken over the last 100 iterations and never before iteration 200: with W = 200 the replay stops
  cpfast's pooled 479-short phase at iteration 200, 550 iterations before it fathoms.
- **Replayed** on the 13 phases of cpfast's `wc` runs, the refit rule stops none of the 11 that fathom
  and neither root early. On the 66 phases recorded at `w` it would stop 13, six of them recorded
  fathoms that took 1,290 to 3,341 iterations: at `w` those were the slow tails.

### In the search: the combined arms [m]

Arms, on top of `zallhge` (the second study's): `zallhw` = `dnn.beta_sched = w=0.36,wc=2.88` and
`cut.handover = 20:0:0.22:1.5:100:3:1.25`; `zallhwg` = `zallhw` + `cut.cp_wc_gap = 300`; `zallhws` =
`zallhwg` + `cut.cp_sweeps = 3:1` + `cut.cp_extrap = 100:1:30:2`. One `bnbdrv` at a time, the box at load
2-5. Each arm's A/B root was run first (`wroot_zallhw`: 200 + 2,550 iterations, 3,510 short left,
where `zallh`'s 3,030 iterations left 4,050); `zallhwg` resumes `zallhw`'s root, whose phase starts
far beyond the gap. `python3 stop/ab.py`:

@@ stop/ab.py

For reference, cpfast's `z030` + `wc=2.88` runs of the same nodes: `18:32` 38 nodes, 5,891 iterations
(5,605 with the gap rule); `38:32` 1,001.

- **Together: `18:32` 4,322 iterations and 38 nodes, `38:32` 410, `30:32` 1,687.** Against `z030`'s 13,035
  and 3,541 that is 3.0x and 8.6x; against the second study's arm 1.28x, 3.2x and 1.19x; against
  cpfast's `z030` + `wc` runs 1.36x (1.30x with their gap rule) and 2.4x.
- **The two results overlap rather than multiply** on `18:32`: each alone took 13,035 to about 5,500
  to 5,900 (2.2x to 2.4x), both together to 4,322 (3.0x). Both act on the same long phases: the
  hand-over rule had already stopped the ones `wc` shortens most. Where one far phase is the whole
  subtree they do multiply: `38:32`'s phase is 300 iterations (3,341 under `z030`, 1,202 under
  `zallhge`, 801 in cpfast's run from the unpooled start).
- **What `wc` changed under the hand-over policy.**
  - The depth-2 root's phase stops at 2,550 iterations 3,510 short (3,030 and 4,050 at `w`), so the
    depth-3 nodes start closer again: `18:32` 2,225 short (2,551), `30:32` 1,662 (1,807).
  - `18:32`'s own phase is handed over at 620 iterations as before but 1,486 short instead of
    1,888; 31 of its 37 children fathom at once (30) and none branches (one did): no depth 5.
  - Its children's phases: the longer ones are shorter (552 short: 300 iterations, where 685 short took
    750; 340: 210, where 540 took 405), the nearest are not (340-430 short: 158-210 iterations at
    `wc`, where 396-471 short took 93-150 at `w`). cpfast's flat start: a pooled phase at the new
    penalty does not move for 100-200 iterations.
  - `30:32`'s root is handed over after 370 iterations, 757 short, and all 37 children fathom in their
    first solve: one phase in the subtree.
- **`cut.cp_wc_gap = 300` changes nothing on these three**: no phase in them starts within 300 of the
  cutoff. `30:32` was run under both arms and is the same run to the iteration; the other two were
  not rerun. It matters on the path (below).
- **A flaw found and guarded on `2:32`.** In the first run of `2:32` a pooled depth-5 phase 353 short,
  at `wc`, had not moved at all by iteration 200, and the rule handed it over on a rate of zero.
  `stop/combined.patch` as delivered waits until the running best has been moving for W iterations
  (and hands over a phase that has not moved by 400). The three subtrees above are unaffected (every
  stop in them came more than W after the first movement: checked on their histories). With the
  guard that phase runs on and fathoms after 1,370 iterations, a geometric tail (gap x 0.7 per 100
  iterations): dearer than the hand-over would have been (about 200 + 35 x 20). A wider
  `cut.cp_wc_gap` (500) would have kept it at `w`. Not tuned further.
- **`2:32`, the class above 5,000 short, measured for the first time** (25% of the A/B subtree's CPU;
  6,281 nodes and 19,941 c7a CPU s in production): **288 nodes, 20,992 iterations, 1,867 CPU s** under
  the combined arm with the guard (`zallhwf`): 22x fewer nodes, 7.7x fewer iterations [production's e],
  7.5x less CPU. Its root starts 2,643 short (5,039 in production) and is handed over at 690
  iterations, 1,936 short; 5 of its 37 children branch: four handed over 400 to 1,060 short after
  370 to 490 iterations, one after 1,610 iterations 15 short (the geometric tail again: the rule
  does not catch it). The tree ends at depth 6. All 5 phases that start beyond 1,000 from a pool hand
  over; of the 16 that start within 1,000, 14 fathom.
- **Three sweeps and the extrapolation (`zallhws`) change nothing in the subtrees**: `18:32` 4,274
  iterations (4,322), `38:32` 400 (410), `30:32` 1,692 (1,687), at 7% more per phase iteration. They
  shorten the one long far phase, the depth-2 root's: 1,990 iterations to the same 3,500 short
  (2,550), 299 CPU s against 366. Under the hand-over policy the phases inside a subtree are too
  short for it to matter. Not in the block; `2:32` was not run under it.

### Where the iterations are now: first solves [m]

`python3 stop/firstsplit.py ...`:

@@ stop/firstsplit.py w18_z030 cpfast/w18_wc2.88 cpfast/w18_wc2.88_gap300 w18_zallhge w18_zallhw w38_zallhw w30_zallhw w2_zallhwf

- **The early entry already removes most of it.** In cpfast's `wc` run of `18:32` the first solves are
  56% of the iterations (3,323 of 5,891): 1,960 at the 6 nodes that give up, 327 each, and 1,363 at the
  32 that fathom. Under `zallhw` a give-up costs **110** (`dnn.giveup_warm = 100`, `dnn.gu_zone = 0`):
  770 at 7 nodes, 18% of the subtree.
- **What is left** in the first solves:
  - the 110 itself. It is the give-up rule's floor (two windows of 50 iterations of history, then
    the first check); the per-node lab started its phases after 50 to 100. Going below needs a
    change in `latest/dnn.C`'s rule (a shorter window for warm nodes), not a key;
  - the first solves of the children that fathom there: 31 nodes, 1,201 iterations in `18:32` (39
    each, 28% of the subtree), and 37 nodes, 1,207 in `30:32` (72% of that subtree, which has one
    phase). No entry policy touches these; a better pool does (fewer iterations to fathom), and so
    would a cheaper iteration.

### The path under the combined arm

`python3 pathverify.py zallhwg` (the patched solver, `par/zallhwg.par`, incumbent 240,518):

```
@@ pathverify.py zallhwg
```

Nothing above 240,516; all three branch; the first solves at depths 7 and 9 are `cold.tsv`'s. **The cost
on the path did not rise**, because the rule hands over where the higher penalty is slow:

| run (iterations: first solves + phases, CPU s) | `z030` | `zallhge` | `zallhwg` |
|---|---|---|---|
| depth 7, the node alone | 200 + 940, 76 | 200 + 880, 71 | 200 + 1,180, 82 |
| depth 9 and its 31 children | 2,208 + 4,000, 347 | 1,286 + 2,230, 156 | 1,283 + 1,810, 121 |

- Depth 7 (3,351 short): `wc` is slower there, as cpfast found: 239,156 after 1,180 phase iterations
  where `w` had 239,180 after 880. 300 iterations more on this node.
- Depth 9 (1,125 short): the phase is stopped after 290 iterations at 239,862 (240,211 after 650 under
  `zallhge`), so its pool is weaker and the path child starts 194 short instead of 47. That child is
  within `cut.cp_wc_gap`, keeps `w`, and reaches 240,515.9998 in 1,520 iterations (1,580 before).
- Not checked: depths 8 and 11-12, where cpfast measured `wc` 2 to 4 times slower on single nodes; a
  path child that starts beyond the 300 gap would run its long phase at `wc`. The gap rule is in the
  block for that reason, not for the subtrees, where it changed nothing.

### Recommended block (replaces the two above)

`experiments/estimates/tho40/tho40.par` with two of its lines changed,

```
dnn.gu_zone = 0                          (was 0.003)
dnn.beta_sched = w=0.36,wc=2.88          (was w=0.36)
```

and these added (this is `par/zallhwf.par`, less `cut.damp1`):

```
dnn.giveup_warm = 100
pool.on = 1
cut.grid = 5:8
cut.side = 1
cut.cp = 100:1000
cut.warmy = 1
cut.rounds = 12
cut.iters = 200
cut.freeze = 1000
cut.zone = 0.06
cut.handover = 20:0:0.22:1.5:100:3:1.25
cut.cp_wc_gap = 300
```

- The last two keys need `stop/combined.patch` (`git apply` on 2f2f4f1b, then `make -C latest`; it also
  carries cpfast's `cut.cp_wc_deficit` and `cut.cp_extrap`, unused here, and their `dnncuts` options).
  A merge also needs the keys' lines in RUNNING.md's table and the settings goldens.
- **Without any patch**: drop the last two lines. That is the second study's block with `wc=2.88`
  under the stall rule, which nobody has run in this form: cpfast ran `z030` + `wc` (5,891 on `18:32`),
  I ran `zall` + rule + early entry without `wc` (5,523).
- Not in the block: `cut.cp_sweeps = 3:1` and `cut.cp_extrap` (1.28x on the depth-2 root's phase, nothing
  on the three subtrees: above), `cut.cp_wc_deficit` (cpfast's rule for the path: the hand-over rule did that job on the
  path nodes checked).

### The whole tree, restated again [e]

`python3 stop/tree.py zallhwf,zallhwg,zallhw wroot_zallhw`: the class table once more, each class from
its subtree under the combined arm, and now with `2:32` measured:

@@ stop/tree.py zallhwf,zallhwg,zallhw wroot_zallhw

- **About 165 c7a CPU-days (140-210) if the whole tree responds as the A/B subtree**: 7.8x (6.2-9.4x)
  on 1,300. The second study had 225 (190-330), the first 550 (450-700).
- **What changed in the estimate's footing:** every class but the 1.2% of near nodes is now measured
  on a warm subtree under the arm itself, four subtrees holding 99% of the A/B subtree's CPU between
  their classes (36% of it in the four nodes themselves). The range is +-20% on each class factor
  for the nodes of the class not run.
- **What it still rests on:** one depth-2 subtree of 468, smaller than average; CPU from a shared box
  against production's c7a seconds x 0.70 (in iterations with phases weighted 1.25, `18:32` is 7.8x
  where its CPU says 8.8x, so the CPU factors may flatter by a tenth); the depth-2 phase, now 2,550
  iterations a node, measured on this one depth-2 node; nothing at depth 1; the incumbent given.
- **Before a launch**, as before: the same arm on a second depth-2 subtree with a far root, warm from
  depth 1, and Knuth dives under it.

## Files

| | |
|---|---|
| `par/*.par` | the arms (`z030g`: `z030` with `dnn.gu_zone = 0`, run on t4; `z030c`: prepared, not run) |
| `run.sh`, `queue.sh`, `jobs.txt`, `queue.log`, `path_jobs.txt` | the runner, the job list as it ended (reordered as results came in; four jobs not run), the log with the order run |
| `runs/NAME_ARM/` | `cmd.txt`, `log.txt`, `nodes.txt` (node trace), `dnn.txt` (DNN trace; with each phase's history in the runs started after 22:48), `err.txt`. The warm roots' state files (44 MB each) were deleted: the `wroot` jobs rewrite them in 26 and 147 CPU s |
| `summarize.py`, `table.py`, `wtable.py`, `wchildren.py`, `matched.py`, `phases.py`, `payoff.py` | the tables above |
| `calib.py`, `pick.py`, `prodstats.py`, `replay.py` | the production trace: matched nodes, subtree choice, CPU by class, the replay model |
| `cutstate.py`, `pathverify.py`, `mkreadme.py`, `README.tmpl.md` | the warm-root state cut, the path verdict, this file's assembly |
| `stop/hist.py`, `stop/replay_rules.py`, `stop/replay_rules.out` | the stop-rule study: the phase histories, the offline replay and its output |
| `stop/handover.patch` | the `cut.handover` key, as built in the worktree `/home/user/wt-tho40search` |
| `stop/run2.sh`, `stop/seq.sh`, `stop/jobs.txt`, `stop/seq.log` | its runner (patched solver for the arms named `*h*`), the one-at-a-time queue, its log |
| `stop/ab.py`, `stop/tree.py` | the A/B table and the restated extrapolation |
| `par/zallh.par`, `par/zallhg.par`, `par/zallhge.par` | `zall` + `cut.handover = 20:2:3`; + `dnn.gu_zone = 0`; + `dnn.giveup_warm = 100` |
| `stop/combined.patch` | third study: cpfast's patch and the hand-over key (with its exponent, marginal condition and flat-start guard), the diff of the worktree against 2f2f4f1b |
| `stop/replay3.py`, `stop/replay3.out`, `stop/firstsplit.py` | third study: the refit and replay on `w` and `wc` histories; where the iterations are |
| `par/zallhw.par`, `zallhwg.par`, `zallhwf.par`, `zallhws.par` | `zallhge` + `wc=2.88` + the refit rule; + `cut.cp_wc_gap = 300`; the same settings run with the guarded solver; + sweeps and extrapolation |
| `runs/.w2_zallhwg_stopped/`, `runs/.old_wroot_zall/` | the traces of a run stopped part-way (before the guard) and of a root later rerun |
