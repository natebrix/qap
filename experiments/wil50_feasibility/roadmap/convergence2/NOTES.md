# convergence2: working notes (Nov 4, started 20:30 UTC)

Running notes; REPORT.md is the deliverable.

## Facts checked in code (cognac 99097cc)

- Search DNN call (`latest/bnb.C` ~663-900): `hopeless_after` = 300 from `-x 0:4:1:300:1`
  with the extrapolation rule (factor 5), `max_iter` 4000 (default), checks every 10.
- Cut phase in the search: `DNNCutParams::warm_primal` = 0 (default; bnb.C never sets it),
  so **every Polyak round and the quiet round start from Y = e0e0'** with the previous Z.
  `dnncuts -warmY` exists (same_node_continuation), never wired into bnb.C.
- Polyak step (`dnncut.C` ~990): alpha = theta (best + delta - lb_round) / |g|^2 with
  lb_round the round's *certified* bound (lags the Lagrangian value), g = violations of
  the round's *last* primal iterate.
- Rounds run their full 300 iterations (`hopeless_after = 0` inside rounds); only a
  fathom ends a round early.
- Pooled children: `dnn_pool_bound` with dp.warm_primal = 1 (WARM_Y=2), pool frozen.

## Trace breakdown, node-10 subtree through the coordinator (`coord_search/trace.txt`)

`trace_breakdown.py coord_search/trace.txt`: 322 nodes, first solves 12,016 s, cut phases
6,448 s (35%).

| class | n | mean it | CPU s |
|---|---|---|---|
| d4 cold fathom | 38 | 162 | 1,918 |
| d4 cold hopeless | 9 | 533 (310-990) | 1,603 |
| d5 cold fathom | 22 | 339 (100-1510) | 2,187 |
| d5 cold hopeless | 3 | 1,007 | 930 |
| d5 warm fathom (no pool) | 19 | 49 | 297 |
| d5 pooled warm fathom | 138 | 68 | 2,986 |
| d6 pooled cold fathom | 21 | 174 | 885 |
| d6 pooled warm fathom | 69 | 48 | 899 |
| cut phases | 13 | 600-2,560 it | 6,448 |

Cold = no parent state (unit boundaries; DNN state does not travel in state files).
93 of 322 nodes cold, about 7,500 s (40%) of DNN CPU.

## 21:10 findings from existing logs

- **Cut-phase lag is the Y reset** (`round_lag.py`, `logs/round_lag_budget_study.txt`):
  median primal - certificate at the end of a 300-iteration round is 44-67 with Y reset
  (budget study arm B = today's search) and 4-10 with Y restored (Y12); the quiet round
  then adds 37-60 (B) vs 5.5-10 (Y12).
- **Cold siblings** (`sibling_warm_cold.py`, `logs/sibling_warm_cold.txt`): same parent,
  same pool status, warm vs cold (parent state lost at a unit/task boundary):
  fathoms 49 vs 339 it (d5, no pool), 51 vs 174 (d6 pooled), 70 vs 195 and 70 vs 268
  (d5 pooled, parallel.py run); misses give up at 320-330 warm vs 872-1,007 cold.
- 16 depth-4 children cold (stronger_bounds/depth4/uncut) vs warm Z+Y (warm_y):
  7 shared fathoms 2,440 vs 1,360 iterations (1.8x) (`cold_vs_warm.py`).
