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

| subtree | arm | nodes | branched | fathomed: first solve / phase | phases (fathomed) | iterations: first solves + phases = total | of which the cold root | CPU s | iterations vs base | nodes vs base |
|---|---|---|---|---|---|---|---|---|---|---|
| gate | base | 37 | 1 | 36 / 0 | 0 (0) | 1,127 + 0 = 1,127 | 590 | 106 |  |  |
| gate | z003 | 1 | 0 | 0 / 1 | 1 (1) | 590 + 130 = 720 | 720 | 68 | 1.57x | 37.0x |
| t1 | base | 106 | 3 | 103 / 0 | 0 (0) | 2,626 + 0 = 2,626 | 320 | 221 |  |  |
| t1 | z003 | 37 | 1 | 35 / 1 | 1 (1) | 1,266 + 130 = 1,396 | 320 | 145 | 1.88x | 2.9x |
| t1 | z010 | 1 | 0 | 0 / 1 | 1 (1) | 320 + 220 = 540 | 540 | 53 | 4.86x | 106.0x |
| t2 | base | 211 | 6 | 205 / 0 | 0 (0) | 6,187 + 0 = 6,187 | 200 | 510 |  |  |
| t2 | z010 | 37 | 1 | 32 / 4 | 4 (4) | 2,812 + 817 = 3,629 | 200 | 323 | 1.70x | 5.7x |
| t2 | z030 | 1 | 0 | 0 / 1 | 1 (1) | 200 + 2,070 = 2,270 | 2,270 | 262 | 2.73x | 211.0x |
| t3 | base | 287 | 8 | 279 / 0 | 0 (0) | 7,281 + 0 = 7,281 | 200 | 734 |  |  |
| t3 | z010 | 38 | 1 | 33 / 4 | 4 (4) | 2,249 + 1,152 = 3,401 | 200 | 393 | 2.14x | 7.6x |
| t3 | z030 | 38 | 1 | 37 / 0 | 1 (0) | 907 + 3,400 = 4,307 | 3,600 | 583 | 1.69x | 7.6x |
| t4 | z010 | 37 | 1 | 28 / 8 | 8 (8) | 4,198 + 4,042 = 8,240 | 200 | 765 | 2.14x [e] | 18.6x |
| t4 | z030 | 37 | 1 | 34 / 2 | 3 (2) | 1,836 + 2,130 = 3,966 | 2,010 | 446 | 4.44x [e] | 18.6x |

- Every phase in these runs but two fathomed its node; the two are the cold roots of t3 and t4 under
  `z030` (next section).
- **t4 base was not rerun** (about 1,400 CPU s, 45 minutes of the shared box): its ratios use
  production's 688 nodes and 688 x 25.6 iterations [e]. gate, t1 and `z030` arms where the root fathoms
  are the same under any zone that covers the root, so they were run once.

## Warm-rooted depth-3 subtrees under `z030` against production [m, production iterations e]

`python3 wtable.py`:

| node (warm root) | production: shortfall, nodes, CPU s (x 0.70 here) | arm | nodes | fathomed: first solve / phase | phases (fathomed) | iterations: first + phases = total | CPU s here | nodes vs production | iterations vs production [e] | CPU vs production x 0.70 |
|---|---|---|---|---|---|---|---|---|---|---|
| `26:10,8:0,18:32` | 4474, 1,490, 5,551 (3,886) | z030 | 74 | 65 / 7 | 9 (7) | 4,474 + 8,561 = 13,035 | 1439 | 20x | 2.93x | 2.70x |
| `26:10,8:0,18:32` | 4474, 1,490, 5,551 (3,886) | zall | 10 (stopped part-way) | 7 / 2 | 3 (2) | 1,098 + 2,780 = 3,878 | 469 | - | - | - |
| `26:10,8:0,38:32` | 2427, 287, 1,000 (700) | z030 | 1 | 0 / 1 | 1 (1) | 200 + 3,341 = 3,541 | 455 | 287x | 2.07x | 1.54x |

`w18`, child by child (`python3 wchildren.py w18_z030 w18_zall`; children that fathom at once in
production and here are left out; `zall` is the next section's arm):

| node | production: shortfall, subtree nodes, iterations [e] | w18_z030 | w18_zall |
|---|---|---|---|
| 26:10,8:0,18:32 (root) | 4474, 1490, - | B 200 it (4474 short) + 920 it stalled 2421 short | B 200 it (3163 short) + 650 it stalled 2298 short |
| 39:7 | 1006, 37, 947 | D 69 it | D 61 it |
| 37:7 | 1891, 107, 2,739 | C 340 it (517 short) + 500 it fathomed | C 350 it (437 short) + 380 it fathomed |
| 36:7 | 268, 37, 947 | D 50 it | D 45 it |
| 35:7 | 2214, 177, 4,531 | C 200 it (849 short) + 1900 it fathomed | C 200 it (795 short) + 1750 it fathomed |
| 34:7 | 953, 37, 947 | D 77 it | D 70 it |
| 33:7 | 1203, 37, 947 | D 95 it | D 92 it |
| 30:7 | 2047, 142, 3,635 | C 340 it (648 short) + 1290 it fathomed | - |
| 28:7 | 556, 37, 947 | D 51 it | - |
| 25:7 | 737, 37, 947 | D 64 it | - |
| 24:7 | 1628, 142, 3,635 | C 350 it (382 short) + 227 it fathomed | - |
| 17:7 | 1644, 72, 1,843 | C 340 it (214 short) + 97 it fathomed | - |
| 16:7 | 1821, 107, 2,739 | C 340 it (567 short) + 630 it fathomed | - |
| 11:7 | 500, 37, 947 | D 52 it | - |
| 10:7 | 444, 37, 947 | D 46 it | - |
| 6:7 | 2001, 72, 1,843 | C 300 it (493 short) + 387 it fathomed | - |
| 4:7 | 187, 37, 947 | D 50 it | - |
| 3:7 | 1258, 72, 1,843 | D 140 it | - |
| 2:7 | 2451, 246, 6,298 | B 200 it (1037 short) + 2610 it stalled 146 short | - |
| all nodes traced so far | 1490 nodes, 38,144 it [e] | 74 nodes, 13,035 it | 10 nodes, 3,878 it |

## The phases [m]

`python3 phases.py` (every phase of every cut arm, the path check aside), by the shortfall of the solve
before the phase:

```
short     0.7 d5 pool -> FATH after   44 it (6 s), left     0.0; first solve 400 it  t4_z030
short    29.3 d5      -> FATH after   30 it (3 s), left     0.0; first solve 620 it  t4_z010
short    88.0 d4      -> FATH after   44 it (6 s), left     0.0; first solve 360 it  t3_z010
short   159.2 d5 pool -> FATH after  276 it (28 s), left     0.0; first solve 350 it  t4_z030
short   182.9 d5      -> FATH after   65 it (7 s), left     0.0; first solve 370 it  t4_z010
short   185.4 d4      -> FATH after   59 it (7 s), left     0.0; first solve 350 it  t3_z010
short   213.7 d4 pool -> FATH after   97 it (11 s), left     0.0; first solve 340 it  w18_z030
short   231.5 d5      -> FATH after   73 it (7 s), left     0.0; first solve 350 it  t2_z010
short   236.1 d5      -> FATH after   78 it (7 s), left     0.0; first solve 340 it  t2_z010
short   302.0 d5      -> FATH after  140 it (14 s), left     0.0; first solve 350 it  t4_z010
short   382.5 d4 pool -> FATH after  227 it (26 s), left     0.0; first solve 350 it  w18_z030
short   384.6 d4      -> FATH after  130 it (14 s), left     0.0; first solve 590 it  gate_z003 ROOT(cold)
short   410.7 d5      -> FATH after  140 it (13 s), left     0.0; first solve 340 it  t2_z010
short   437.3 d4 pool -> FATH after  380 it (47 s), left     0.0; first solve 350 it  w18_zall
short   450.4 d5      -> FATH after  130 it (14 s), left     0.0; first solve 300 it  t1_z003
short   493.1 d4 pool -> FATH after  387 it (44 s), left     0.0; first solve 300 it  w18_z030
short   517.3 d4 pool -> FATH after  500 it (57 s), left     0.0; first solve 340 it  w18_z030
short   566.6 d4 pool -> FATH after  630 it (71 s), left     0.0; first solve 340 it  w18_z030
short   577.4 d5      -> FATH after  202 it (20 s), left     0.0; first solve 320 it  t4_z010
short   648.2 d4 pool -> FATH after 1290 it (141 s), left     0.0; first solve 340 it  w18_z030
short   795.5 d4 pool -> FATH after 1750 it (209 s), left     0.0; first solve 200 it  w18_zall
short   822.8 d4      -> FATH after  276 it (33 s), left     0.0; first solve 200 it  t3_z010
short   849.4 d4 pool -> FATH after 1900 it (203 s), left     0.0; first solve 200 it  w18_z030
short   957.3 d4      -> FATH after  220 it (24 s), left     0.0; first solve 320 it  t1_z010 ROOT(cold)
short  1036.8 d4 pool -> miss after 2610 it (299 s), left   146.0; first solve 200 it  w18_z030
short  1138.3 d5      -> FATH after  415 it (41 s), left     0.0; first solve 200 it  t4_z010
short  1246.7 d5      -> FATH after  526 it (51 s), left     0.0; first solve 200 it  t2_z010
short  1284.0 d5      -> FATH after  561 it (56 s), left     0.0; first solve 200 it  t4_z010
short  1442.7 d4      -> FATH after  773 it (92 s), left     0.0; first solve 200 it  t3_z010
short  1538.8 d5      -> FATH after  912 it (88 s), left     0.0; first solve 200 it  t4_z010
short  1907.9 d5      -> FATH after 1717 it (163 s), left     0.0; first solve 200 it  t4_z010
short  2426.7 d3      -> FATH after 3341 it (431 s), left     0.0; first solve 200 it  w38_z030 ROOT(warm)
short  3162.8 d3 pool -> miss after  650 it (90 s), left  2298.2; first solve 200 it  w18_zall ROOT(warm)
short  3192.6 d4      -> FATH after 2070 it (243 s), left     0.0; first solve 200 it  t2_z030 ROOT(cold)
short  3462.2 d3      -> miss after 3400 it (482 s), left     2.9; first solve 200 it  t3_z030 ROOT(cold)
short  4274.3 d4      -> miss after 1810 it (228 s), left   832.1; first solve 200 it  t4_z030 ROOT(cold)
short  4474.4 d3      -> miss after  920 it (120 s), left  2420.8; first solve 200 it  w18_z030 ROOT(warm)
short  8634.2 d2      -> miss after  760 it (122 s), left  5100.2; first solve 200 it  wroot_zall ROOT(warm)

shortfall bin: phases, fathomed, median it (fathomed), median it (missed), median gain of misses
    0-250      9   9     65   nan      nan
  250-500      7   7    140   nan      nan
  500-750      4   4    630   nan      nan
  750-1000     4   4   1750   nan      nan
 1000-1500     5   4    561  2610    890.8
 1500-2000     2   2   1717   nan      nan
 2000-2500     1   1   3341   nan      nan
 2500-3500     3   1   2070  3400   3459.3
 3500-5000     2   0    nan  1810   3442.2
 5000-99999    1   0    nan   760   3534.0
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

`python3 replay.py` is a cruder cross-check on the same trace (a phase fathoms from at most SF short at
measured iteration costs; after a miss the children either all fathom, `all`, or get no help, `none`):
`z003` x1.35 CPU, `z010` x1.9, `z030` x1.6 (none) to x5 (all). The class table sits inside that bracket
because the pool's carry is partial.

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
depth  7 path node  B first solve 237165.4112 (200 it), phase best 239204.0841 after 940 it (stalled)  -> ok
   (the path node alone: -d 7)
depth  9 path node  B first solve 239391.3413 (570 it), phase best 240494.2955 after 3400 it (iter_limit)  -> ok
depth 10 path child B first solve 240515.9974 (800 it, pooled 995 cuts), phase best 240516.0000 after 600 it (stalled)  -> ok
   siblings of the path child: 30 bounded, 30 fathomed; run complete
PASS
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
