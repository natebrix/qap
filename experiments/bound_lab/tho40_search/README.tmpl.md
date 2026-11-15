# tho40: the grid cut phase in the search (bound_lab agent `tho40_search`, Nov 15 2026)

Question (`../COMMON.md`): tho40 has never run a cut phase. On complete subtrees, end to end, what is
the grid cut phase worth, and which `cut.zone`? **[m]** measured here (run named), **[e]** estimated.

## Findings

SUMMARY_GOES_HERE

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
  t4: 4,274 and 2,675). So on t2-t4 the root is outside `z010` here although in a run it would be
  inside (t2) or near it: read `z010` on t2-t4 as "no phase at the root, phases on its children".
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
- **Farther out it still fathoms, but the price per unit rises**: 0.4-0.6 at 1,250-1,550 short
  (depth 5), 0.65 at 3,193 (t2's cold root, depth 4: 2,070 iterations), 1.4 at 2,427 on a depth-3 node
  (`w38`: 3,341 iterations, 59 short of the 3,400 limit).
- **A phase that starts from an inherited pool is slower per unit**: the pool has already taken the
  easy part. In `w18`: 214 short 97 iterations, 382: 227, 517: 500, 567: 630, 648: 1,290, 849: 1,900.
- **Misses.** t4's cold root (4,274 short): stalled after 1,810 iterations, 832 short. `w18`'s warm
  root (4,474): stalled after 920, 2,421 short, still gaining 0.8 a iteration. t3's cold root (3,462,
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
(production CPU shares, `prodstats.py` and the class table in `replay.py`'s neighbour `pick.py`):

| depth-3 class (short) | nodes of that class | share of the subtree's CPU | CPU factor under `z030` | from |
|---|---|---|---|---|
| fathomed at once, or <= 1,500 | 12 | 1.2% | about 1.5x | [e] gate, t1 (deeper nodes of that shortfall) |
| 1,500-2,500 | 11 | 7.4% | 1.0-1.6x, take 1.2x | [m] `w38`, the largest of them, 1.6x; the small ones break even (above) |
| 2,500-3,500 | 5 | 12.3% | about 2x | [e] between `w38` and `w18`; not run (`w30` was queued) |
| 3,500-5,000 | 9 | 53.7% | W18CPU | [m] `w18` (4,474 short, 1,490 nodes) |
| above 5,000 | 1 (`2:32`, 6,281 nodes) | 25.3% | as `w18` | [e] not run; its depth-4 children are t2/t4-like (2.7x and 4.4x in iterations, cold) |

CLASS_TEXT

`python3 replay.py` is a cruder cross-check on the same trace (a phase fathoms from at most SF short at
measured iteration costs; after a miss the children either all fathom, `all`, or get no help, `none`):
`z003` x1.35 CPU, `z010` x1.9, `z030` x1.6 (none) to x5 (all). The class table sits inside that bracket
because the pool's carry is partial.

**The whole tree** (3.4e7 nodes, about 1,300 c7a CPU-days by Knuth's estimator, itself "likely low"):
WHOLE_TEXT

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

  PATH_TEXT
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

NEXT_TEXT

## Files

| | |
|---|---|
| `par/*.par` | the arms (`z030c`, `z030g`: prepared, not run) |
| `run.sh`, `queue.sh`, `jobs.txt`, `queue.log`, `path_jobs.txt` | the runner, the job list in the order run, its log |
| `runs/NAME_ARM/` | `cmd.txt`, `log.txt`, `nodes.txt` (node trace), `dnn.txt` (DNN trace; with the phase's history from 23:00 on), `err.txt` |
| `summarize.py`, `table.py`, `wtable.py`, `wchildren.py`, `phases.py`, `payoff.py` | the tables above |
| `calib.py`, `pick.py`, `prodstats.py`, `replay.py` | the production trace: matched nodes, subtree choice, CPU by class, the replay model |
| `cutstate.py`, `pathverify.py`, `mkreadme.py`, `README.tmpl.md` | the warm-root state cut, the path verdict, this file's assembly |
