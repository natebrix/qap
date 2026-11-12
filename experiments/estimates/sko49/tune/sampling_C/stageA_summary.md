# Stage A: depths 0-2 of sko49t.dat (finished; 344 tasks)

- binary: commit `7c950e5f6461`, latest/ modified: False, sources newer than the binary: ['dnn.C', 'dnn.h', 'bnb.C'] (built 2026-11-11 18:22:16)
- `bnbdrv -p params.par -s latest/SYM/sko49t.sym -u 23386 -d <depth> -r <one node> -o ...`; the params file is `experiments/estimates/sko49/tune/C.par` plus the --env settings; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.025
- symmetry file loaded in 344 of 344 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 0.7 | 0 | 0 | 10 | 0 |
| 1 | 10 | 0 | 10 | 0 | 3,000 + 0 = 3,000 | 1,136 + 0 + 13 = 1,150 | 115.0 | 0 | 0 | 333 | 0 |
| 2 | 333 | 20 | 333 | 20 | 100,680 + 68,998 = 169,678 | 38,541 + 28,248 + 157 = 66,946 | 201.0 | 333 | 54 | 13,682 | 0 |
| all | 344 | 20 | 343 | 20 | 103,680 + 68,998 = 172,678 | 39,677 + 28,248 + 171 = 68,096 | | | 54 | | |

Task wall time summed over slots: 68,815 s; CPU/wall 0.990.

Cutoff 23384.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 10 | 300 / 300 / 300 | 450.8 / 608.3 / 642.2 | 114.0 |
| 2 | fathomed | 1 | 260 / 260 / 260 | -0.0 / -0.0 / -0.0 | 98.7 |
| 2 | hopeless | 332 | 300 / 300 / 780 | 2.2 / 352.5 / 511.5 | 114.9 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 14 | 15 | 6 | 0 | 5 | 52 | 54 | 65 | 51 | 46 | 5 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 19 | 2.2 / 97.7 / 155.9 | -0.6 / -0.1 / -0.0 | 479 | 6,672 |
| 2 | iter_limit | 2 | 158.4 / 161.6 / 164.9 | 0.3 / 3.6 / 7.0 | 3400 | 2,757 |
| 2 | stalled | 33 | 162.5 / 205.2 / 233.5 | 11.0 / 82.3 / 128.1 | 1190 | 18,819 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 41 (10)
- depth 1: 10 nodes, 10 branched: 3 on a col, 7 on a row
  `41:0` row 30 (27); `41:1` col 0 (48); `41:2` col 0 (48); `41:3` row 40 (27); `41:8` row 30 (27); `41:9` col 0 (48); `41:10` row 40 (27); `41:16` row 40 (27); `41:17` row 40 (27); `41:24` row 40 (27)
- depth 2: 333 nodes, 313 branched: 188 on a col, 125 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `41:0` | 450.8 | row 30 | 27 | 27 | 8 | row 46 x14, row 9 x2, row 40 x2, col 42 x1 | 788 |
| `41:1` | 533.5 | col 0 | 48 | 48 | 7 | col 7 x25, row 46 x8, col 48 x5, row 18 x1 ... | 1,927 |
| `41:2` | 562.6 | col 0 | 48 | 48 | 3 | col 1 x16, col 6 x9, row 46 x7, col 7 x6 ... | 2,115 |
| `41:3` | 566.8 | row 40 | 27 | 27 | 1 | col 0 x14, row 35 x5, col 42 x5, row 46 x1 ... | 1,096 |
| `41:8` | 596.0 | row 30 | 27 | 27 | 1 | row 46 x15, col 0 x5, row 9 x3, col 48 x1 ... | 1,117 |
| `41:9` | 620.7 | col 0 | 48 | 48 | 0 | col 48 x17, row 46 x9, col 6 x6, col 1 x5 ... | 2,256 |
| `41:10` | 623.7 | row 40 | 27 | 27 | 0 | col 0 x13, row 35 x6, col 42 x5, row 9 x2 ... | 1,143 |
| `41:16` | 640.5 | row 40 | 27 | 27 | 0 | col 6 x8, row 35 x6, row 9 x5, col 0 x4 ... | 1,143 |
| `41:17` | 642.2 | row 40 | 27 | 27 | 0 | col 0 x9, col 42 x7, row 35 x6, row 9 x3 ... | 1,143 |
| `41:24` | 641.0 | row 40 | 27 | 27 | 0 | row 35 x14, col 6 x6, col 48 x3, col 0 x3 ... | 954 |

Children per depth-2 node: 0: 20 nodes, 26: 49 nodes, 47: 264 nodes
**Frontier at depth 3: 13,682 nodes** (Burnside orbit count 29,580 for wil50 at depth 3, roadmap/tree_size.md section 1, before any fathoming)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 637 | 684 | 261 | 0 | 235 | 2,339 | 2,454 | 2,908 | 2,166 | 1,805 | 193 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 0.3, 10% 103.5, median 358.1, 90% 464.3, max 511.5. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1152.3, median 1674.2, max 2027.6.
Parents kept with their DNN state for stage B: 211 of 333 (tau history [[0.025, 0]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.029094128750720396; parents of the 400 that stage B must rebuild: 16

Sanity checks:
- the root: 10 children (expected 15 orbits)
- depth 2 before bounding: 333 candidate children of the 10 depth-1 nodes (333 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
