## Setup

Binary: `de4e0697b8cb` (commit fe31484, `latest/` unmodified: True); diagnostic binary `5182e9fbecaa` (the same sources plus the QAP_XBAR_DEBUG=2 print). Root: QPB bound 46131.508, line row 21, 15 children; CPU 1.0 s.

| node | DNN bound (HEAD) | shortfall | iterations | CPU s | rule 9's line | lookahead score row / col | rule 9's line's rank by rule 8's score | diagnostic build: trace, children, warm start |
|---|---|---|---|---|---|---|---|---|
| r.1 `21:1` | 48273.0884 | 541.0 | 300 | 114.8 | col 0 | col 0: 114,558,810.7 / 114,585,133.4 | 2 | equal, equal, equal |
| r.5 `21:10` | 48243.9954 | 570.1 | 300 | 113.4 | col 0 | col 0: 114,599,870.2 / 114,614,395.8 | 2 | equal, equal, equal |

Rule 8's score (sum of squared DNN marginals) of the ten best lines, rows and columns ranked together (the candidates in bold):

- r.1: **row 45** 0.2019, **col 0** 0.1930, **col 49** 0.1520, **row 25** 0.1517, col 9 0.1512, col 40 0.1337, row 26 0.1215, row 0 0.1185, row 13 0.0995, row 19 0.0914
- r.5: **row 45** 0.2436, **col 0** 0.2312, **col 40** 0.1673, **col 49** 0.1471, row 0 0.1470, col 9 0.1336, row 26 0.1317, row 25 0.1316, col 20 0.1263, row 13 0.1108

Every child run (490) resumed with `--allow-settings-change`; its warning lists cut.maxdepth and par.hash as the only differences from the node's settings record.

## The rewrite check: rule 9's children rewritten against straight from the solver

| node | children | identical first lb, iterations, status | largest abs. lb difference | same depth-2 line | predicted bound: production (min-max) / rewrite |
|---|---|---|---|---|---|
| r.1 | 49 | 49 | 0.0000 | 49 of 49 | 46464.3-46959.7 / 48273.0884 (inherited) |
| r.5 | 49 | 49 | 0.0000 | 49 of 49 | 46476.6-46959.7 / 48243.9954 (inherited) |

## HEAD against stage A (a667727) on rule 9's children (production predicted bounds)

| node | children matched | identical first lb | median / largest abs. difference | HEAD higher | same status | stage A fathomed in the first solve / HEAD |
|---|---|---|---|---|---|---|
| r.1 | 49 | 49 | 0.000 / 0.000 | 0 | 49 | 2 / 2 |
| r.5 | 49 | 49 | 0.000 / 0.000 | 0 | 49 | 2 / 2 |

## Per node and line (49 depth-2 children each, first DNN solve, no cut phase)

Predicted bound of every child: the inherited bound (the node's DNN bound), for every line, rule 9's rewritten included; "prod" is rule 9's children straight from the solver (production's predicted bounds), the rewrite check. Ratios are rule 9's score over the line's: above 1, the line is better. E: the expected-subtree score (model.md) with the step model (bins) and the log-linear one (between the bins' mean shortfalls), 95% bootstrap over stages A and B.

| node | line | rule 8 score (rank) | children | fathomed | in the zone, not fathomed | beyond 244 | fathomed after a phase, modelled (95%) | shortfall min / median / max | lines cost 48^(gap/130) | rule 9 / line | top-3 share | E bins (95%) | rule 9 / line (95%) | E log-linear | rule 9 / line (95%) | CPU s (mean iterations) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r.1 | col 0 (rule 9) | 0.1930 (2) | 49 | 2 | 28 | 19 | 16.0 (15.5-16.6) | 10.1 / 203.0 / 430.9 | 1.10e+06 | 1.00 | 0.68 | 1.21e+05 (5.72e+04-1.93e+05) | 1 | 9.55e+04 | 1 | 5,212 (293) |
| r.1 | row 45 | 0.2019 (1) | 49 | 0 | 36 | 13 | 20.9 (20.2-21.7) | 31.5 / 164.1 / 407.0 | 5.54e+05 | 1.98 | 0.67 | 6.51e+04 (3.97e+04-9.82e+04) | 1.85 (1.08-2.62) | 5.84e+04 | 1.64 (1.13-2.22) | 5,627 (308) |
| r.1 | col 49 | 0.1520 (3) | 49 | 1 | 25 | 23 | 11.7 (11.3-12.1) | 19.2 / 243.3 / 420.0 | 1.21e+06 | 0.91 | 0.58 | 1.37e+05 (8.20e+04-2.02e+05) | 0.88 (0.61-1.08) | 1.29e+05 | 0.74 (0.58-0.82) | 5,441 (296) |
| r.1 | row 25 | 0.1517 (4) | 49 | 0 | 0 | 49 | 0.0 (0.0-0.0) | 258.0 / 306.6 / 398.6 | 9.07e+05 | 1.21 | 0.32 | 1.51e+05 (9.11e+04-2.23e+05) | 0.80 (0.37-1.50) | 1.42e+05 | 0.67 (0.38-1.03) | 5,407 (300) |
| r.1 | col 0 (rule 9, prod) | 0.1930 (2) | 49 | 2 | 28 | 19 | 16.0 (15.5-16.6) | 10.1 / 203.0 / 430.9 | 1.10e+06 | 1.00 | 0.68 | 1.21e+05 (5.72e+04-1.93e+05) | 1.00 (1.00-1.00) | 9.55e+04 | 1.00 (1.00-1.00) | 5,231 (293) |
| r.5 | col 0 (rule 9) | 0.2312 (2) | 49 | 2 | 25 | 22 | 11.5 (11.1-12.0) | 42.7 / 230.6 / 457.8 | 2.35e+06 | 1.00 | 0.71 | 1.40e+05 (7.48e+04-2.15e+05) | 1 | 1.40e+05 | 1 | 5,185 (298) |
| r.5 | row 45 | 0.2436 (1) | 49 | 0 | 34 | 15 | 20.2 (19.3-21.1) | 52.7 / 166.1 / 437.1 | 1.29e+06 | 1.83 | 0.67 | 1.19e+05 (5.60e+04-1.93e+05) | 1.17 (1.08-1.40) | 1.05e+05 | 1.33 (1.24-1.54) | 5,520 (307) |
| r.5 | col 40 | 0.1673 (3) | 49 | 0 | 20 | 29 | 5.8 (5.4-6.2) | 36.8 / 265.7 / 454.0 | 2.03e+06 | 1.16 | 0.65 | 1.64e+05 (9.57e+04-2.44e+05) | 0.85 (0.75-0.92) | 1.54e+05 | 0.91 (0.80-0.99) | 5,260 (301) |
| r.5 | col 49 | 0.1471 (4) | 49 | 1 | 19 | 29 | 8.1 (7.6-8.5) | 46.1 / 271.2 / 449.2 | 2.83e+06 | 0.83 | 0.56 | 2.12e+05 (1.13e+05-3.28e+05) | 0.66 (0.61-0.70) | 1.94e+05 | 0.72 (0.63-0.87) | 5,443 (300) |
| r.5 | col 0 (rule 9, prod) | 0.2312 (2) | 49 | 2 | 25 | 22 | 11.5 (11.1-12.0) | 42.7 / 230.6 / 457.8 | 2.35e+06 | 1.00 | 0.71 | 1.40e+05 (7.48e+04-2.15e+05) | 1.00 (1.00-1.00) | 1.40e+05 | 1.00 (1.00-1.00) | 5,203 (298) |

Shortfall histogram (first solve; stage A's bins; the zone ends at 244.08):

| node | line | fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r.1 | col 0 (rule 9) | 2 | 2 | 7 | 6 | 7 | 6 | 9 | 5 | 1 | 4 | 0 | 0 | 0 |
| r.1 | row 45 | 0 | 2 | 13 | 7 | 9 | 5 | 4 | 4 | 4 | 1 | 0 | 0 | 0 |
| r.1 | col 49 | 1 | 3 | 5 | 3 | 6 | 8 | 6 | 9 | 5 | 3 | 0 | 0 | 0 |
| r.1 | row 25 | 0 | 0 | 0 | 0 | 0 | 0 | 21 | 21 | 7 | 0 | 0 | 0 | 0 |
| r.1 | col 0 (rule 9, prod) | 2 | 2 | 7 | 6 | 7 | 6 | 9 | 5 | 1 | 4 | 0 | 0 | 0 |
| r.5 | col 0 (rule 9) | 2 | 2 | 4 | 4 | 7 | 8 | 8 | 7 | 3 | 3 | 1 | 0 | 0 |
| r.5 | row 45 | 0 | 0 | 12 | 10 | 9 | 3 | 5 | 3 | 3 | 4 | 0 | 0 | 0 |
| r.5 | col 40 | 0 | 1 | 2 | 3 | 7 | 7 | 10 | 12 | 3 | 3 | 1 | 0 | 0 |
| r.5 | col 49 | 1 | 1 | 2 | 5 | 4 | 7 | 7 | 9 | 7 | 6 | 0 | 0 | 0 |
| r.5 | col 0 (rule 9, prod) | 2 | 2 | 4 | 4 | 7 | 8 | 8 | 7 | 3 | 3 | 1 | 0 | 0 |

CPU: 490 child tasks, 53,529 CPU s (first solves 53,377 s), wall per task 111 s on average; prep (root and 4 depth-1 runs) 455 CPU s.
