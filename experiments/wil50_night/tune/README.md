# wil50's tuning round (Nov 25): both candidates hold, 0.80 of the CPU together

`wil50-tune-261125-{both,noxl2,nogap,neither}` (cognac `infra/phase1/specs/wil50-tune.spec`, records in
`infra/phase1/runs/`): four arms, each solving the same 133 depth-4 open nodes of `wil50-est-261124` to the end
(150 asked; drawn in proportion to the square root of the dives' estimates), one Spot worker an arm
(r7a.16xlarge; `nogap` on c7a.16xlarge, the same EPYC 9R14), no interruption. `paired.txt` is
`experiments/estimates_nov19/paired.py` on the four records, `both` the baseline.

| arm | recipe | CPU h, weighted | against `both` (5-95%) | nodes against `both` | the run itself |
|---|---|---|---|---|---|
| `both` | `wil50_d2r10_g150xl2.par` (the estimate's) | 1,791.7 | 1 | 1 | 21,270 nodes, 172.5 CPU h |
| `noxl2` | without `branch.xbar_lines = 2` | 1,920.5 | 1.07x (1.02-1.14) | 1.07x | 23,267 nodes, 188.0 CPU h |
| `nogap` | without `cut.cp_wc_gap = 150`, `cut.cp_extrap = 100:1:30:2` | 2,078.0 | 1.16x (1.13-1.18) | 1.12x | 23,876 nodes, 202.2 CPU h |
| `neither` | `wil50_d2r10.par` | 2,248.5 | 1.25x (1.20-1.32) | 1.22x | 26,409 nodes, 220.5 CPU h |

- **The two effects multiply.** `branch.xbar_lines = 2` is 0.93 of the CPU with the cut phase's keys and 0.92
  without; the keys are 0.86 with it and 0.85 without; together 0.80 (nodes 0.82).
- **Against the screening** (five subtrees): the shortlist read 0.80 there and is 0.93 here; the cut phase's
  keys read 0.93 there and are 0.86 here, and here they also cut nodes (0.89), which the screening's two
  subtrees did not show.
- **Where the gain is:** in the large subtrees. The 20 sampled nodes with over 200 nodes under `both` hold
  1,418 of the 1,792 weighted CPU h; on them `noxl2` is 1.07x, `nogap` 1.18x, `neither` 1.27x. Node by node the
  arms cross often (`noxl2` is cheaper than `both` at 48 nodes and dearer at 33): the result is the weighted sum's.
- **The dives against the solved sample, depth 4:** the same 133 nodes weigh 1,223 CPU h by their dives'
  estimates and 1,792 solved: 1.47x (bootstrap 5-95%: 0.91 to 2.45), so the dives read about 0.68 of the
  truth at depth 4 on wil50 (sko42: 0.80 at depth 4). One node the dives put at 0.6 CPU h took 20.2; one put
  at 71.7 took 27.1. `wil50-est-261124`'s 10,145 CPU h scaled by that ratio is about 15,000 CPU h, 10,000 to
  25,000 by its band; at depths 2 and 3, where 86% of the estimated rest is, the ratio is not measured.
