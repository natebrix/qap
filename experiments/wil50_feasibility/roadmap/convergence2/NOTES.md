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
- Search phases (`phase_split.py`, `logs/phase_split_coord.txt`): 19,100 phase iterations,
  4,740 (25%) in quiet rounds; near-near misses (3.5-7.6 short after 990-1,190 uncut
  iterations) spend 870-1,370 phase iterations (3+ full rounds) to gain 8-9.

## 22:45 test results (runs/, summarize_ep.py)

- Measurement gate passed 22:12; tests on cores 0-1 from 22:12.
- E0: node 10's phase reproduces the search's round 0 exactly (48,582.4671 at 300);
  +63/+90/+102 in rounds 1-3, overshoot -18 in round 4, stall; 1,200 phase iterations,
  446 CPU-s; pool (best round) 2,392 cuts kept by the children.
- E children (cold, pool frozen, give-up 300): the 4 phase fathoms of the search
  (48, 38, 28, 8 on location 45) fathom on the first solve at 1,100/530/350/450 it
  (search 3,030/2,600/1,860/2,020 incl. phase); phase misses 29, 9 end 44/42 short
  (search 130/106 short, 79/64 after their phase); 49:45 (out of zone in the search,
  its subtree 137 nodes 8,002 CPU-s = 43%) ends 82.6 short (search 185.7); easy
  fathoms 130->120, 200->150.
- E2 (49:45 with the pool, then the search's phase from it, Y reset as today): round 1
  drops 48,731.5 -> 48,615.6 (-116), rounds 2-3 recover to 48,706.8, stall after 3 rounds
  with no gain over the pooled solve (stopped in the quiet round, which the search would
  not run: primal 48,767.8 < cutoff). E2w = the same with -warmY (pending).
- E2w (the same with -warmY): round 1 still drops (-88, to 48,643.4), rounds 2-3 recover to
  48,746.5 (+15 over the pooled solve), stall: the drop is the first Polyak step on
  inherited near-optimal multipliers (target best + 50), not only the Y reset.
- S (parent 13 at 300, Y and Z rounded to bf16): its 4 children take 340/300/310/90
  iterations, exactly as unrounded (warm_y), bounds at 50/100/200 within 0.8.
- C n6: 780 (B 1,260, Y12 860, Y12s 780); C n9: 700 (B 990, Y12 810, Y12s 720).
- C misses (23:39): n13 final 48,741.39 (B 48,754.62: -13.2), n10 48,720.74 (B 48,743.88:
  -23.1) after 2,100 iterations; overshoots -302 (n13 round 7), -88/-51 (n10 rounds 4, 8).
  **Y12q fails the bound gate on misses; 100-iteration rounds are killed.** Y12s (200, quiet
  1,000) passed on both (+5.2, +4.5) in the budget study: carry Y12s.
- The report text is in the session's final message (the harness refused a REPORT.md
  file from this subagent); `logs/C_summary.txt`, `logs/E_P_S_summary.txt` hold the tables.
- C n14 (Y12q: warm Y, 12 x 100, quiet 500): fathoms at 1,680 total iterations
  (B 4,750, Y12 3,640, Y12s 2,110); round 0 48,692.0959 vs the study's 48,692.0828
  (solver changed since; not bit-identical).
- P (parents 13, 7 at 1,000 vs 300 iterations): fathoming children 950 -> 830 it
  (-13%), misses unchanged (give-up 300-310, +4..+12 at 200). Break-even against the
  parent's extra 700. Closed.
