# Experiment: lookahead on a restricted set of candidate lines

At depths 3-4, default.par runs lookahead on the 5 best rows and 5 best
columns by U sum. Each of those lines has m children, and every child gets a
QP bound. `restricted_lookahead.patch` (against 2c1e399, which has the rank
fix) lets the candidate set be chosen differently. All switches are
environment variables and apply from depth QAP_LA_MINFIX on (default 3):

- `QAP_LA_ECC=1`: rank the free columns (locations) by eccentricity, the
  largest B entry to another free location, instead of by U sum.
- `QAP_LA_KC=k`: how many candidate columns to bound (default use_best).
- `QAP_LA_KR=k`: how many candidate rows to bound, still ranked by U sum
  (default use_best).

## Large instances (ub = optimum, default.par at 2c1e399; 3 solvers at a time on 4 cores)

Times are in seconds, with the ratio to the current default in brackets.

| instance | current: 5 rows + 5 cols by U | 2 + 2 by U | 3 eccentric cols only | 3 eccentric cols + 2 rows by U |
|---|---|---|---|---|
| nug20 | 529,427 nodes / 157.8 s | 566,216 / 122.2 (0.77) | 517,179 / 110.8 (0.70) | 525,610 / 127.0 (0.80) |
| nug21 | 856,864 / 202.0 | 947,003 / 175.6 (0.87) | 888,420 / 160.6 (0.79) | 899,878 / 176.4 (0.87) |
| tai17a | 474,858 / 181.6 | 534,431 / 126.2 (0.70) | 629,479 / 137.4 (0.76) | 530,976 / 142.7 (0.79) |
| chr15a | 5,023,059 / 305-307 | 7,836,758 / 368.3 (1.21) | 6,340,761 / 324.3 (1.06) | 6,887,159 / 345.7 (1.13) |
| esc16a | 10,140,245 / 642.4 | 10,226,125 / 756.1 (1.18) | stopped after 2,268 s (> 3.5) | 7,307,144 / 446.1 (0.69) |

The geometric mean of the time ratios is 0.92 for 2 + 2 by U and 0.84 for
3 eccentric + 2 rows. Columns only is a disaster on esc16a, as rule 6 was:
its hypercube distances make eccentricity meaningless, and there it is the
rows (facilities) that matter.

On small instances (nug16a, nug17, had16, chr12a, tai15a, rou15, scr15),
eccentric versus U ranking of 5 columns made no difference. The saving
there comes from bounding fewer lines.
