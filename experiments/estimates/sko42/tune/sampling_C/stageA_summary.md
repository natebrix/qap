# Stage A: depths 0-2 of sko42t.dat (finished; 451 tasks)

- binary: commit `dc5e5d761d07`, latest/ modified: False, sources newer than the binary: none (built 2026-11-11 18:22:16)
- `bnbdrv -p params.par -s latest/SYM/sko42t.sym -u 15812 -d <depth> -r <one node> -o ...`; the params file is `experiments/estimates/sko42/tune/C.par` plus the --env settings; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.01454
- symmetry file loaded in 451 of 451 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 0.5 | 0 | 0 | 12 | 0 |
| 1 | 12 | 0 | 12 | 0 | 3,600 + 0 = 3,600 | 558 + 0 + 4 = 562 | 46.9 | 0 | 0 | 438 | 0 |
| 2 | 438 | 165 | 438 | 165 | 127,256 + 232,156 = 359,412 | 19,828 + 38,495 + 60 = 58,383 | 133.3 | 438 | 218 | 10,704 | 0 |
| all | 451 | 165 | 450 | 165 | 130,856 + 232,156 = 363,012 | 20,386 + 38,495 + 65 = 58,946 | | | 218 | | |

Task wall time summed over slots: 59,559 s; CPU/wall 0.990.

Cutoff 15810.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 12 | 300 / 300 / 300 | 231.4 / 349.5 / 393.0 | 46.5 |
| 2 | fathomed | 41 | 20 / 73 / 1391 | -10.8 / -0.2 / -0.0 | 11.9 |
| 2 | hopeless | 397 | 300 / 300 / 330 | 8.7 / 148.8 / 292.3 | 47.0 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 165 | 70 | 24 | 0 | 82 | 68 | 29 | 0 | 0 | 0 | 0 | 0 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 124 | 8.7 / 65.8 / 121.7 | -1.1 / -0.1 / -0.0 | 348 | 13,257 |
| 2 | iter_limit | 10 | 87.9 / 115.6 / 147.2 | 0.6 / 2.5 / 4.3 | 3400 | 5,667 |
| 2 | stalled | 84 | 94.2 / 135.9 / 158.0 | 10.2 / 39.5 / 76.9 | 1325 | 19,571 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 36 (12)
- depth 1: 12 nodes, 12 branched: 8 on a col, 4 on a row
  `36:0` row 22 (41); `36:1` col 0 (41); `36:2` col 0 (41); `36:3` row 22 (23); `36:7` col 0 (41); `36:8` col 0 (41); `36:9` col 0 (41); `36:10` row 22 (23); `36:14` col 0 (41); `36:15` col 0 (41); `36:16` col 0 (41); `36:17` row 22 (23)
- depth 2: 438 nodes, 273 branched: 180 on a col, 93 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `36:0` | 231.4 | row 22 | 41 | 41 | 36 | col 41 x3, col 35 x1, col 6 x1 | 200 |
| `36:1` | 307.5 | col 0 | 41 | 41 | 23 | row 22 x6, col 7 x3, col 6 x3, row 12 x2 ... | 720 |
| `36:2` | 334.5 | col 0 | 41 | 41 | 19 | col 6 x9, row 22 x6, col 41 x3, row 7 x1 ... | 880 |
| `36:3` | 341.7 | row 22 | 23 | 23 | 5 | col 6 x7, col 35 x3, col 41 x3, col 0 x3 ... | 684 |
| `36:7` | 303.2 | col 0 | 41 | 41 | 23 | col 35 x6, col 1 x3, row 22 x3, col 41 x2 ... | 720 |
| `36:8` | 357.3 | col 0 | 41 | 41 | 15 | row 22 x14, col 41 x6, row 12 x2, row 5 x1 ... | 1,040 |
| `36:9` | 379.9 | col 0 | 41 | 41 | 9 | col 41 x16, row 22 x11, col 6 x2, row 5 x1 ... | 1,280 |
| `36:10` | 385.6 | row 22 | 23 | 23 | 0 | col 6 x8, col 35 x5, col 41 x5, row 24 x4 ... | 830 |
| `36:14` | 326.1 | col 0 | 41 | 41 | 19 | col 35 x10, row 22 x7, col 41 x2, row 5 x1 ... | 880 |
| `36:15` | 375.0 | col 0 | 41 | 41 | 10 | col 35 x14, col 41 x13, row 22 x4 | 1,240 |
| `36:16` | 391.0 | col 0 | 41 | 41 | 6 | col 41 x21, row 22 x8, col 35 x6 | 1,400 |
| `36:17` | 393.0 | row 22 | 23 | 23 | 0 | col 35 x6, col 41 x5, col 6 x5, row 24 x4 ... | 830 |

Children per depth-2 node: 0: 165 nodes, 22: 12 nodes, 40: 261 nodes
**Frontier at depth 3: 10,704 nodes** (Burnside orbit count 29,580 for wil50 at depth 3, roadmap/tree_size.md section 1, before any fathoming)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2,728 | 960 | 0 | 3,208 | 2,648 | 1,160 | 0 | 0 | 0 | 0 | 0 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 0.6, 10% 25.4, median 180.9, 90% 246.5, max 292.3. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1002.1, median 1279.2, max 1503.0.
Parents kept with their DNN state for stage B: 125 of 438 (tau history [[0.01454, 0]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.036803292551284505; parents of the 400 that stage B must rebuild: 88

Sanity checks:
- the root: 12 children (expected 15 orbits)
- depth 2 before bounding: 438 candidate children of the 12 depth-1 nodes (438 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
