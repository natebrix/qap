# Stage A: depths 0-2 of wil50t.dat (IN PROGRESS; 67 tasks)

- binary: commit `c4f4dc4760f6`, latest/ modified: False, sources newer than the binary: none (built 2026-11-07 18:05:31)
- `bnbdrv -p nug_deep_xbar9_all.par -x 1:4:1:300:1 -s latest/SYM/wil50t.sym -u 48816 -d <depth> -r <one node> -o ...`; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.0
- symmetry file loaded in 67 of 67 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 1.1 | 0 | 0 | 15 | 0 |
| 1 | 15 | 0 | 15 | 0 | 4,500 + 0 = 4,500 | 1,711 + 0 + 12 = 1,723 | 114.9 | 0 | 0 | 635 | 0 |
| 2 | 51 | 37 | 51 | 37 | 16,965 + 30,356 = 47,321 | 6,086 + 11,907 + 8 = 18,001 | 353.0 | 51 | 40 | 652 | 0 |
| all | 67 | 37 | 66 | 37 | 21,465 + 30,356 = 51,821 | 7,798 + 11,907 + 21 = 19,726 | | | 40 | | |

Task wall time summed over slots: 19,928 s; CPU/wall 0.990.

Cutoff 48814.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 15 | 300 / 300 / 300 | 343.7 / 456.8 / 574.4 | 114.5 |
| 2 | fathomed | 7 | 91 / 402 / 1330 | -0.4 / -0.1 / -0.0 | 150.8 |
| 2 | hopeless | 44 | 300 / 310 / 370 | 22.4 / 95.4 / 402.4 | 109.9 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 37 | 5 | 3 | 2 | 0 | 0 | 3 | 0 | 0 | 1 | 0 | 0 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 30 | 22.4 / 66.0 / 128.5 | -1.6 / -0.2 / -0.0 | 310 | 6,303 |
| 2 | stalled | 10 | 140.4 / 164.0 / 228.8 | 7.8 / 47.0 / 126.2 | 1380 | 5,604 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 45 (15)
- depth 1: 15 nodes, 15 branched: 3 on a col, 12 on a row
  `45:0` col 49 (49); `45:1` row 21 (49); `45:2` row 21 (49); `45:3` row 21 (49); `45:4` row 21 (49); `45:10` row 21 (49); `45:11` row 21 (49); `45:12` row 21 (49); `45:13` col 0 (49); `45:14` col 0 (49); `45:20` row 21 (29); `45:21` row 21 (29); `45:22` row 21 (29); `45:23` row 31 (29); `45:24` row 31 (29)
- depth 2: 51 nodes, 14 branched: 12 on a col, 2 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `45:0` | 574.4 | col 49 | 49 | 1 | 0 | row 21 x1 | 48 |
| `45:1` | 534.0 | row 21 | 49 | 3 | 1 | col 9 x1, col 49 x1 | 96 |
| `45:2` | 494.1 | row 21 | 49 | 7 | 5 | col 49 x2 | 96 |
| `45:3` | 456.8 | row 21 | 49 | 3 | 2 | col 0 x1 | 48 |
| `45:4` | 432.7 | row 21 | 49 | 6 | 6 |  | 0 |
| `45:10` | 560.3 | row 21 | 49 | 3 | 3 |  | 0 |
| `45:11` | 505.6 | row 21 | 49 | 4 | 3 | col 49 x1 | 48 |
| `45:12` | 449.4 | row 21 | 49 | 9 | 7 | col 49 x2 | 96 |
| `45:13` | 406.0 | col 0 | 49 | 0 | 0 |  | 0 |
| `45:14` | 379.6 | col 0 | 49 | 4 | 4 |  | 0 |
| `45:20` | 549.6 | row 21 | 29 | 4 | 2 | col 0 x1, col 9 x1 | 96 |
| `45:21` | 491.9 | row 21 | 29 | 1 | 1 |  | 0 |
| `45:22` | 427.7 | row 21 | 29 | 3 | 2 | col 9 x1 | 48 |
| `45:23` | 376.2 | row 31 | 29 | 2 | 0 | row 25 x1, col 0 x1 | 76 |
| `45:24` | 343.7 | row 31 | 29 | 1 | 1 |  | 0 |

Children per depth-2 node: 0: 37 nodes, 28: 1 nodes, 48: 13 nodes
**Depth-2 sample: 51 of 635 nodes, 652 children; expected frontier at depth 3: 8,118 nodes (1 SE 1,789)** [e: 635 x the sample's mean children]
Root line forced: facility 45 in the 15 J1 cells (the rule's line was row 21; out.rule.state)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 220 | 144 | 96 | 0 | 0 | 144 | 0 | 0 | 48 | 0 | 0 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 7.8, 10% 21.9, median 85.0, 90% 282.9, max 402.4. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1558.9, median 1888.0, max 2156.4.
Parents kept with their DNN state for stage B: 0 of 51 (tau history [[0.0, 0]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.6152580063397829; parents of the 400 that stage B must rebuild: 14

Sanity checks:
- the root: 15 children (expected 15 orbits)
- depth 2 before bounding: 635 candidate children of the 15 depth-1 nodes (635 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
