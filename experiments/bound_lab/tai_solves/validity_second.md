- **A second route for each** (`runs/bfs35h`, `runs/bfs40c`; `logs/solve_bfs35h.txt`, `logs/solve_bfs40c.txt`),
  the other region family from its own short root (the 1,500-phase-iteration roots of section 2), so the root's
  pool, every inherited pool and every phase differ:

  | search | regions | root | nodes | DNN iterations | CPU-s | optimum found, no open node | path: nodes, violations, largest bound |
  |---|---|---|---|---|---|---|---|
  | tai35b main | `comb` (42) | 280,789,962 | 922 | 125,232 | 5,275 | yes | 32, 0, 283,315,444.9994 |
  | tai35b second | `combh` (74) | 279,907,632 | 890 | 112,396 | 4,882 | yes | 32, 0, 283,315,444.9999 |
  | tai40b main | `combh` (85) | 636,110,765 | 733 | 75,582 | 4,862 | yes | 37, 0, 637,250,947.9995 |
  | tai40b second | `comb` (48) | 633,078,560 | 726 | 69,745 | 4,551 | yes | 37, 0, 637,250,947.9998 |

  The second routes are no bigger than the main ones although their roots are 0.9M and 3.0M lower: the
  6,000-iteration roots did not pay for themselves on trees this small (tai40b: 1,227 against 359 CPU-s at the root).
