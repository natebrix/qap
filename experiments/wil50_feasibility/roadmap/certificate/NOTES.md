# Certificate study: running notes (Nov 4, from 20:45 UTC)

Crash-recovery log; the findings are in REPORT.md.

## 20:45-21:10 reading, code, smoke tests (box busy until ~22:00)

- Worktree branch was at master 16a72e8; reset to cognac 99097cc (clean).
- Certificate anatomy (latest/dnn.C certify): lb = s + M00 + GL(M) - (m+1) scale
  lambda - margin, M = L + scale Zc, Zc = Z + beta (Y - VRV') (PRSM half-step
  multiplier; V'ZcV = beta * negative part of W exactly). GL inner LAPs: row-col
  reduction, exact Hungarian only within 1e-3 |cutoff| (48.8 units) of the cutoff.
- Margins in the wil50 logs: 3.6e-6 to 4.95e-6 units, refine=0, cert_fail=0
  (amx/logs/val_d4n3_amx1.log, warm_y/runs/child_9.log). The PSD term and the
  rounding margins have nothing to give.
- bnb.C never sets DNNParams::child_bounds: child_lb is unused (confirmed).
- warm_y/runs/child_9.log (uncut, warm child, pass z): primal ~48,813 from
  iteration 150 on, certificate 48,729.7 (150), 48,780.1 (300), 48,804.3 (600),
  48,810.4 (900). The lag is the dual converging slowly, not the primal.
- Traces: every cut-phase node of pool_search / coord_search had its uncut
  primal BELOW the cutoff at give-up (-1.9 to -112): no certificate change on
  the uncut solve could have spared any of the 13 phases.
- coord_search CPU (18,496 s): fathoming first solves ~9,200 s (cold d4/d5
  4,100; pooled 4,770), first-solve misses ~2,700, cut phases ~7,000 (from
  pool_search's split).
- Pooled children first-solve iterations (pool_search): 6 at 10, 47 at 20, 40
  at 30, ..., tail to 500. A free child bound can save at most the first check.

## Instrumentation (latest/dnn.C, DNN_CERTX, opt-in, default path unchanged)

Variants at the same iterate: 1 exact GL, 2 Hahn-Grant skew polishing, 4 line
search Zc + t (Z - Zc), 8 average of last 4 Zc, 16 dump M. dnntest -K takes a
';' list. Smoke: nug12 to convergence all variants agree with the DNN value
(567.990); instrumented and plain runs give bit-identical RESULT lines;
brute-force n = 8 (scripts/validity_bruteforce.py): no variant above the
enumerated optimum on 12 (sub)problems. child_gl.py's exact parent GL equals
dnn.C's to 4 decimals (nug12 dump).

nug20 node 0:1, 100 iterations: HG +4 to +40 over std, exact +0.6 to +7,
line search never above t = 0 (concave, decreasing), average below std while
the certificate rises.

## 21:00-21:20 more variants, validity

- Bit 32 (HGR): HG with row-column reductions in the inner blocks (O(m^4) a
  pass) and one exact pass at the end. nug20 node: 12+1 passes reach the 8-exact
  -pass value (2,524.91 vs 2,524.82 at it 100) at about half the time at m = 19.
- Bit 64 (warm skew K carried across checks): WORSE than plain exact GL
  (2,496.9 vs 2,521.2 at it 100, nug20): a K tuned to one check's M is stale at
  the next. The invariance claim itself holds: at the same check
  GL_exact(M + K) >= polish value (+0.17 to +7.1, DNN_CERTX_WTEST). Dead.
- validity_qaplib.py (all variants, every check, 1,500 its): nug12, had12,
  chr12a, tai12a, tai12b: no value above the optimum. tai12b is a sharp case:
  HG lifts the converged certificate 39,464,306 -> 39,464,921 (optimum
  39,464,925). logs/valid_*.log.

## 21:20-21:30 one more variant, queue ready

- Bit 128 (proj): the post-step Z (box-dual exact) with its face part projected
  in double (one dsyevd), exact GL. nug20 node: +17 over exact GL at it 40, +6
  at 60, +0.12 at 100 (HG keeps +3.6 there).
- Early exit in the polish (USE mode only) once the cutoff is crossed.
- validity_qaplib.py with all variants (239): PASS (proj and warm included).
- queue.sh 2: E1 cold children 47 (with ls), 37, 27, 6, 46 (modes 175/171,
  window 50, every 5th check), then E1w warm children 7, 47.
- queue.sh 3: E4 node 10 dump + child bounds on column 45; E2 [9:45] phase
  (mode 187, dump, child bounds on row 19); E2 [48:45] phase (mode 171).

## 22:12-22:30 first wil50 results (measurement ended 22:12; cores 2, 3)

E1, cold children of node 10, production flags, same iterates (logs/e1_*.log):
| child | std fathom | HG8 | HGR13 | gain at last check before std fathom (HG / HGR) |
|---|---|---|---|---|
| 47:45 (d4, m 46) | 450 | 420 | 430 | +3.19 / +3.03 (lag 6.7) |
| 49:45,0:12 (d5, m 45) | 650 | 610 | 620 | +2.09 / +2.03 (lag 13.4) |
| 7:45 (d4) | 310 | 280 | 280 | +9.28 / +8.68 (lag 9.6) |
| 37:45 (d4) | 290 | 260 | 270 | +9.41 / +8.38 (lag 9.6) |
Line search: below exact from it ~100 on (47:45: -0.7 at the end). Average:
-1 to -6. Projected post-step multiplier: below std from it ~150 (-2 to -9 at
the end). Exact GL = std in the window (AUTO already exact there).
CPU per evaluation (m 45-46): exact GL 0.25 s, HGR 0.5 s, HG8 1.4 s, proj 1.5
s, line search 2.4 s; 10 iterations cost about 3.8 s.
E2 [48:45] round ends (lag primal - std 59-69): HG +6.3 to +9.8 (10-14% of
the lag); round 0 at 470 its: +4.6 (39% of an 11.9 lag). Stalled after round 5
(best +36.8), quiet round running.

## 22:30-22:45 E2 [48:45] done, E3

- E2 [48:45] (logs/e2_n48_45_m171.log): rounds 1-5 end 44-135 short with
  lags 59-69; HG +6.3 to +9.8 (10-14%). Stalled after round 5; quiet round from
  round 5's multipliers fathoms at 1,060 (std); HG/HGR above the cutoff at the
  1,000 check (evaluated every 100 its): about 60 of 3,030 phase iterations.
  In the quiet round the gain shrinks 31.5 (it 100) -> 2.4 (700) -> 1.2 (1,000),
  9-17% of the lag throughout. Search: the same node's phase 2,560 its.
- E3 A/B (DNN_CERTX=32, USE=1, HGR used within 15 of the cutoff), plain vs
  polished, same core:
  7:45: 310 its 96.3 s -> 280 its 87.4 s (2 evaluations, 1.35 s)
  47:45: 450 its 145.9 s -> 430 its 151.0 s (5.9 s of polishing; ms/iter
  324 vs 351: run-to-run drift of 8% with the other agent on cores 0-1)
  Iterations are deterministic, CPU is not: count iterations, price
  evaluations at 0.5-0.6 s.
- run_e4.sh lacked +x: E4 skipped by queue 3; rerun on core 2 after E3.

## 22:45-23:05 E3 depth-5, E4, E5, E2b, E1w

- E3 49:45,0:12: 650 its 188.9 s -> 620 its 179.1 s (7.3 s of polishing).
- E4 node 10 at its 300-it give-up (229 short): exact per-child GL on column
  45: median +45.8, max +68.1 over the parent, 0 of 47 cross (160-210 short);
  -C 1: 0 of 2,209. logs/e4_n10_join_col45.log.
- E5 37:45, more passes: at it 280 HG 2/4/32 passes +4.3/+7.3/+12.5, HGR 48+1
  +10.4; both fathom at 260 (std 290; HG8 260, HGR13 270); HG32 0.12 short at
  250. Cost 4.3 s (HG32), 2.1 s (HGR49) per evaluation.
- E2b [9:45] (logs/e2_n9_45_m187.log): round 2 best 48,750.28 = the search's
  63.77 short exactly; stalled after round 4; quiet round (the search ran
  none: primal below the cutoff) converges to 48,793.9 (20.1 short), primal
  48,796.3. Polish in the quiet round +32 (it 100), +3.2 (500), +1.0-1.4
  (1,000-1,500); average +2.4 at 1,500 (oscillation at convergence).
  Q3 on its row 19 (logs/e2_n9_45_join_row19_*.log): exact per-child GL median
  +0.58 over the parent (flat), max +57.6; 2 of 46 cross (19:40, 19:31), which
  the search fathomed in 10 and 30 iterations: 13-15 s of 1,043-1,585 s.
- E1w warm children (node 10's Z and conditioned Y, -Q b): 7:45 std 120 ->
  HG/HGR 100; 47:45 std 210 -> 170. Gains near the end +4.5 to +13 (30-60% of
  the lag).

## 22:55-23:08 E1b/E1c: 13 more fathoming nodes (logs/e1b_table.md)

Every node crosses earlier with the polish (first check above the cutoff,
std -> HG8 / HGR13): d4 cold 20:45 130->120/120, 45:45 150->140/140, 16:45
160->150/150, 26:45 160->140/140, 36:45 170->150/160, 6:45 180->160/160,
46:45 200->180/180, 27:45 240->220/220, 17:45 260->240/240; d5 cold under
49:45: 0:22 380->350/350, 0:17 470->440/450, 0:18 510->480/490, 0:11
530->500/490. With E1: 17 cold nodes, 5,240 std iterations, HGR13 saves 350
(6.7%), HG8 400 (7.6%); 2 warm children 330 -> 270 (18%).
E6 started 23:08 (sub-subtree of [9:45], plain vs polish W=10, cores 2/3).

## 23:08-23:35 E6: search A/B on the [9:45] sub-subtree (logs/e6_*, e6_summary.txt)

bnbdrv serial, run_pool_search.sh's configuration, -F 39:27,30:48,25:49,9:45.
Same tree in both arms (47 nodes, 46 fathomed, FW 3,492). Total 1,490.7 s
(plain) -> 1,418.7 s (polish, W = 10): 1.051x. Root: identical first solve
(430 its) and phase (1,200 its, same pool size; the polish arm reports best
48,758.2 instead of 48,750.3 because round ends are evaluation points, but the
Polyak steps use best + delta - lb_now and the exported pool is the same).
46 pooled children: 3,400 -> 3,020 iterations (-11.2%), 943.0 -> 868.3 s
(1.086x); 32 children unchanged (fast fathoms), 14 save 10-60 iterations; the
8 children of 100+ iterations 2,060 -> 1,740 its, 571.6 -> 503.1 s.

Effect estimate on the node-10 subtree (coord_search split, 18,465 DNN s):
pooled 4,770 s at 1.086x (-378 s), non-pooled warm d5 297 s at ~15% (-45),
cold fathoms 4,105 s at ~5% net (-205), fathomed phases ~4,000 s at ~2%
(-80): about -700 s, 3.8%, 1.04x [e from measured per-class rates].
