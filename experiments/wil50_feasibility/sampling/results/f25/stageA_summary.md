# Stage A: depths 0-2 of wil50t.dat (finished; 216 tasks)

- binary: commit `8dd01bc3bdd2`, latest/ modified: False, sources newer than the binary: none (built 2026-11-07 18:05:31)
- `bnbdrv -p nug_deep_xbar9_all.par -x 1:4:1:300:1 -s latest/SYM/wil50t.sym -u 48816 -d <depth> -r <one node> -o ...`; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.0
- symmetry file loaded in 216 of 216 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 1.1 | 0 | 0 | 15 | 0 |
| 1 | 15 | 0 | 15 | 0 | 4,500 + 0 = 4,500 | 1,623 + 0 + 12 = 1,634 | 108.9 | 0 | 0 | 635 | 0 |
| 2 | 200 | 37 | 200 | 37 | 61,160 + 144,171 = 205,331 | 20,994 + 54,145 + 61 = 75,201 | 376.0 | 200 | 106 | 7,624 | 0 |
| all | 216 | 37 | 215 | 37 | 65,660 + 144,171 = 209,831 | 22,617 + 54,145 + 74 = 76,836 | | | 106 | | |

Task wall time summed over slots: 77,662 s; CPU/wall 0.989.

Cutoff 48814.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 15 | 300 / 300 / 300 | 504.4 / 549.5 / 609.4 | 108.8 |
| 2 | hopeless | 200 | 300 / 310 / 330 | 65.0 / 237.4 / 432.5 | 104.8 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 37 | 25 | 27 | 17 | 0 | 0 | 50 | 32 | 8 | 4 | 0 | 0 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 37 | 65.0 / 122.3 / 156.0 | -0.3 / -0.1 / -0.0 | 730 | 14,568 |
| 2 | iter_limit | 6 | 153.1 / 160.0 / 182.4 | 0.7 / 4.2 / 7.4 | 3400 | 7,649 |
| 2 | stalled | 63 | 156.7 / 204.0 / 242.0 | 8.0 / 79.4 / 135.3 | 1150 | 31,928 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 25 (15)
- depth 1: 15 nodes, 15 branched: 15 on a row
  `25:0` row 45 (49); `25:1` row 45 (49); `25:2` row 45 (49); `25:3` row 45 (49); `25:4` row 45 (49); `25:10` row 45 (49); `25:11` row 45 (49); `25:12` row 45 (49); `25:13` row 45 (49); `25:14` row 45 (49); `25:20` row 21 (29); `25:21` row 21 (29); `25:22` row 21 (29); `25:23` row 21 (29); `25:24` row 21 (29)
- depth 2: 200 nodes, 163 branched: 79 on a col, 84 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `25:0` | 609.4 | row 45 | 49 | 15 | 0 | col 9 x8, col 40 x4, col 49 x2, col 10 x1 | 720 |
| `25:1` | 583.3 | row 45 | 49 | 17 | 0 | col 9 x6, row 21 x5, col 0 x3, col 40 x3 | 816 |
| `25:2` | 553.9 | row 45 | 49 | 20 | 4 | row 21 x9, col 9 x3, col 49 x2, col 0 x1 ... | 768 |
| `25:3` | 546.3 | row 45 | 49 | 13 | 1 | row 21 x8, col 40 x2, col 0 x1, col 49 x1 | 576 |
| `25:4` | 549.5 | row 45 | 49 | 16 | 5 | row 21 x10, col 9 x1 | 528 |
| `25:10` | 602.8 | row 45 | 49 | 15 | 0 | row 21 x6, col 0 x5, col 9 x3, col 40 x1 | 720 |
| `25:11` | 564.2 | row 45 | 49 | 17 | 2 | row 21 x6, col 9 x4, col 0 x2, col 49 x2 ... | 720 |
| `25:12` | 522.7 | row 45 | 49 | 18 | 4 | row 21 x6, col 9 x3, col 0 x2, col 49 x2 ... | 672 |
| `25:13` | 512.0 | row 45 | 49 | 10 | 3 | col 9 x4, row 21 x2, col 40 x1 | 336 |
| `25:14` | 514.8 | row 45 | 49 | 17 | 7 | row 21 x4, col 0 x2, col 49 x2, col 9 x2 | 480 |
| `25:20` | 603.3 | row 21 | 29 | 10 | 0 | row 45 x7, col 0 x2, row 31 x1 | 400 |
| `25:21` | 557.8 | row 21 | 29 | 6 | 2 | row 45 x3, col 0 x1 | 152 |
| `25:22` | 515.7 | row 21 | 29 | 9 | 4 | row 45 x5 | 240 |
| `25:23` | 504.4 | row 21 | 29 | 9 | 3 | row 45 x6 | 248 |
| `25:24` | 506.0 | row 21 | 29 | 8 | 2 | row 45 x6 | 248 |

Children per depth-2 node: 0: 37 nodes, 28: 10 nodes, 48: 153 nodes
**Depth-2 sample: 200 of 635 nodes, 7,624 children; expected frontier at depth 3: 24,206 nodes (1 SE 694)** [e: 635 x the sample's mean children]
Root line forced: facility 25 in the 15 J1 cells (the rule's line was row 21; out.rule.state)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1,140 | 1,296 | 796 | 0 | 0 | 2,320 | 1,516 | 364 | 192 | 0 | 0 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 0.7, 10% 34.2, median 254.7, 90% 342.3, max 432.5. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1586.4, median 1893.5, max 2192.0.
Parents kept with their DNN state for stage B: 0 of 200 (tau history [[0.0, 0]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.05364202802178595; parents of the 400 that stage B must rebuild: 151

Sanity checks:
- the root: 15 children (expected 15 orbits)
- depth 2 before bounding: 635 candidate children of the 15 depth-1 nodes (635 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
