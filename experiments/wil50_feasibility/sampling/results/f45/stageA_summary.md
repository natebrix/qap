# Stage A: depths 0-2 of wil50t.dat (finished; 216 tasks)

- binary: commit `c4f4dc4760f6`, latest/ modified: False, sources newer than the binary: none (built 2026-11-07 18:05:31)
- `bnbdrv -p nug_deep_xbar9_all.par -x 1:4:1:300:1 -s latest/SYM/wil50t.sym -u 48816 -d <depth> -r <one node> -o ...`; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.0
- symmetry file loaded in 216 of 216 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 1.1 | 0 | 0 | 15 | 0 |
| 1 | 15 | 0 | 15 | 0 | 4,500 + 0 = 4,500 | 1,711 + 0 + 12 = 1,723 | 114.9 | 0 | 0 | 635 | 0 |
| 2 | 200 | 113 | 200 | 113 | 63,576 + 135,615 = 199,191 | 22,767 + 52,920 + 38 = 75,725 | 378.6 | 200 | 137 | 4,036 | 0 |
| all | 216 | 113 | 215 | 113 | 68,076 + 135,615 = 203,691 | 24,479 + 52,920 + 51 = 77,450 | | | 137 | | |

Task wall time summed over slots: 78,183 s; CPU/wall 0.991.

Cutoff 48814.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 15 | 300 / 300 / 300 | 343.7 / 456.8 / 574.4 | 114.5 |
| 2 | fathomed | 30 | 20 / 255 / 1500 | -4.1 / -0.1 / -0.0 | 97.2 |
| 2 | hopeless | 170 | 300 / 310 / 620 | 9.5 / 138.2 / 449.5 | 108.1 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 113 | 27 | 19 | 8 | 0 | 0 | 17 | 9 | 3 | 4 | 0 | 0 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 83 | 9.5 / 79.0 / 139.2 | -1.6 / -0.1 / -0.0 | 362 | 19,401 |
| 2 | iter_limit | 4 | 135.1 / 144.9 / 176.0 | 2.7 / 3.0 / 3.3 | 3400 | 5,443 |
| 2 | stalled | 50 | 127.8 / 177.0 / 233.9 | 6.3 / 53.2 / 136.7 | 1295 | 28,076 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 45 (15)
- depth 1: 15 nodes, 15 branched: 3 on a col, 12 on a row
  `45:0` col 49 (49); `45:1` row 21 (49); `45:2` row 21 (49); `45:3` row 21 (49); `45:4` row 21 (49); `45:10` row 21 (49); `45:11` row 21 (49); `45:12` row 21 (49); `45:13` col 0 (49); `45:14` col 0 (49); `45:20` row 21 (29); `45:21` row 21 (29); `45:22` row 21 (29); `45:23` row 31 (29); `45:24` row 31 (29)
- depth 2: 200 nodes, 87 branched: 68 on a col, 19 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `45:0` | 574.4 | col 49 | 49 | 15 | 1 | row 21 x9, col 40 x4, col 48 x1 | 672 |
| `45:1` | 534.0 | row 21 | 49 | 17 | 8 | col 49 x6, col 0 x2, col 9 x1 | 432 |
| `45:2` | 494.1 | row 21 | 49 | 20 | 13 | col 49 x6, col 1 x1 | 336 |
| `45:3` | 456.8 | row 21 | 49 | 13 | 7 | col 49 x4, col 0 x2 | 288 |
| `45:4` | 432.7 | row 21 | 49 | 16 | 9 | col 40 x2, col 49 x2, col 9 x2, col 0 x1 | 336 |
| `45:10` | 560.3 | row 21 | 49 | 15 | 10 | col 0 x4, col 49 x1 | 240 |
| `45:11` | 505.6 | row 21 | 49 | 17 | 10 | col 49 x5, col 9 x2 | 336 |
| `45:12` | 449.4 | row 21 | 49 | 18 | 11 | col 49 x7 | 336 |
| `45:13` | 406.0 | col 0 | 49 | 10 | 6 | col 49 x3, row 21 x1 | 192 |
| `45:14` | 379.6 | col 0 | 49 | 17 | 14 | col 49 x2, col 40 x1 | 144 |
| `45:20` | 549.6 | row 21 | 29 | 10 | 7 | col 0 x1, col 9 x1, row 31 x1 | 124 |
| `45:21` | 491.9 | row 21 | 29 | 6 | 4 | col 9 x1, row 40 x1 | 76 |
| `45:22` | 427.7 | row 21 | 29 | 9 | 6 | col 9 x2, col 49 x1 | 144 |
| `45:23` | 376.2 | row 31 | 29 | 9 | 2 | col 0 x2, row 21 x2, row 25 x1, col 9 x1 ... | 296 |
| `45:24` | 343.7 | row 31 | 29 | 8 | 5 | row 40 x2, row 25 x1 | 84 |

Children per depth-2 node: 0: 113 nodes, 28: 7 nodes, 48: 80 nodes
**Depth-2 sample: 200 of 635 nodes, 4,036 children; expected frontier at depth 3: 12,814 nodes (1 SE 866)** [e: 635 x the sample's mean children]
Root line forced: facility 45 in the 15 J1 cells (the rule's line was row 21; out.rule.state)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1,216 | 892 | 384 | 0 | 0 | 796 | 412 | 144 | 192 | 0 | 0 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 2.7, 10% 15.9, median 90.8, 90% 347.9, max 449.5. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1456.6, median 1898.8, max 2265.6.
Parents kept with their DNN state for stage B: 0 of 200 (tau history [[0.0, 0]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.09637573071611424; parents of the 400 that stage B must rebuild: 86

Sanity checks:
- the root: 15 children (expected 15 orbits)
- depth 2 before bounding: 635 candidate children of the 15 depth-1 nodes (635 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
