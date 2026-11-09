# The expected-subtree model (depth1_lines_score.py)

E(s) for a depth-2 node of first-solve shortfall s (no cut phase): 0 when the first solve fathoms, else q x 48 x m, q = the fraction of stage A's depth-2 nodes in the bin (by first-solve shortfall) the production cut phase did not fathom, m = the mean anchor-model size (tree_estimate.size_anchor) of stage B's sampled depth-3 nodes whose parents are in the bin. 95% bootstrap intervals over stage A's 635 and stage B's 400.

| first-solve shortfall | stage A nodes | their mean shortfall | q (not fathomed after the phase) | stage B nodes | m (mean child size) | largest child size | E at the mean shortfall (95%) |
|---|---|---|---|---|---|---|---|
| 0-100.00 | 176 | 60.9 | 0.000 | 0 | - | 0.0 | 0 |
| 100-150.00 | 121 | 126.6 | 0.248 | 51 | 1.00 | 1.0 | 11.9 (8.29-15.7) |
| 150-200.00 | 96 | 174.5 | 0.927 | 122 | 1.00 | 1.5 | 44.7 (41.9-47) |
| 200-244.08 | 70 | 220.6 | 1.000 | 81 | 1.28 | 10.3 | 61.2 (50.3-75.7) |
| 244-300.00 | 54 | 270.5 | 1.000 | 72 | 16.97 | 256.0 | 814 (430-1.31e+03) |
| 300-350.00 | 32 | 327.0 | 1.000 | 34 | 94.41 | 723.6 | 4.53e+03 (2.31e+03-7.15e+03) |
| 350-400.00 | 15 | 370.6 | 1.000 | 26 | 113.88 | 1523.0 | 5.47e+03 (1.52e+03-1.21e+04) |
| 400- | 14 | 419.7 | 1.000 | 14 | 439.74 | 2209.7 | 2.11e+04 (5.59e+03-3.91e+04) |

Stage A depth-2 nodes fathomed in the first solve: 57 (E = 0).
