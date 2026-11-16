# cpfast: what limits the tail of the CP cut phase on tho40, and a faster phase (bound_lab agent `cpfast`, Nov 16 2026)

Question (brief; `../COMMON.md`): the regional-triangle phase (CP-ADMM, `latest/dnncut.h`) adds 2,000 to
3,500 units to tho40's certified bound but takes 1,000 to 3,400 iterations to do it. What limits the
tail, and can a given bound be reached in half the iterations at the same cost per iteration?

**[m]** measured here (the log is named), **[e]** estimated. Iterations are the metric: the box was
shared (load 5 to 9 on 4 cores). Settings everywhere: `experiments/estimates/tho40/tho40.par` plus the
cut block of `../tho40_search/README.md` (`cut.cp = 100:1000`, `cut.warmy = 1`), cutoff 240,514.05.

## Answer

1. **What limits the tail: the ADMM penalty of the phase.** tho40.par's `dnn.beta_sched = w=0.36` is
   the warm children's factor (tuned so that they fathom in 20 to 40 iterations). With no `wc` key the
   cut phase runs at the same factor, beta = 0.36 x 0.21 x m/3. At that penalty the cut multipliers,
   which move by penalty x violation / ||G||^2 per iteration, creep: their mass h'mu is 4,473 at phase
   iteration 300, 9,058 at 1,000 and 11,021 at 2,000, while the primal objective of the iterate stays
   2,000 to 3,400 units below the certified bound for the whole phase [m `logs/A_base.log`]. At 8 times
   the penalty h'mu is 11,109 at iteration 300 [m `logs/A_wc288.log`]. The pool is not the limit (it
   stops growing by iteration 200 to 300), nor is the accuracy of the Dykstra projection (worth 1.15x),
   nor the eigen-tracker (exact eigendecompositions: same curve).
2. **`dnn.beta_sched = w=0.36,wc=2.88` (an existing key; phase penalty x 8) makes the phase 2.4 to 4.6
   times shorter, at the same cost per iteration, on far-miss nodes off the path at depths 3 to 5** [m]:
   - `26:10,8:0,38:32` (depth 3, 2,427 short) in the search, warm: fathoms after **801 phase
     iterations instead of 3,341** (4.2x; 119 CPU s for the node instead of 455);
   - `26:10,8:0,18:32` (depth 3, 4,474 short, cannot fathom): today's end value (2,421 short after 920
     iterations) after 300 iterations; 1,605 short when the stall rule stops it at 1,010;
   - a depth-4 node 1,812 short: past the cutoff at 390 instead of about 940 [e, 27.7 short at 900];
     a depth-5 node 2,527 short: today's 800-iteration value at 270 (3.0x).
3. **End to end** [m]: the warm subtree `26:10,8:0,18:32` (1,490 nodes in production; 74 nodes and
   13,035 DNN iterations, 1,437 CPU s under tho40_search's `z030`): **38 nodes, @@R3_ITS@@ iterations,
   @@R3_CPU@@ CPU s** with the settings below (2.2x fewer with `wc=2.88` alone). The tree ends one level
   higher: every depth-4 child fathoms.
4. **Two places where the higher penalty must not go** (both measured, both handled by the patch):
   - the change of penalty costs the first 30 to 150 iterations of a phase (the running best stays at
     the first solve's value), which is most of a near miss's phase;
   - where the first solve's primal objective is not below its bound (pooled children, and nodes on
     the path to the best known solution at depths 8 to 10) a higher penalty gains nothing or is
     slower (path node, depth 8: 239,664 at 500 iterations under `wc=2.88` against 240,054 today).
   The first solve's **primal deficit** (certified bound less the primal objective of its last
   iterate) separates the cases: 0.35 to 1.44% of the cutoff at unpooled give-ups off the path, between
   -0.22% and +0.11% at pooled children and on the path. Two keys in `cpfast.patch` use it:
   `cut.cp_wc_deficit` (take `wc` only above this deficit) and `cut.cp_wc_gap` (keep `w` within this
   of the cutoff).
5. **Factor, against the bar**: on the far-miss phases that hold most of a subtree's iterations today,
   3x to 4.6x in iterations to a given bound at equal cost per iteration (bar: 2x). On the one subtree
   run, 2.2x to @@R3_FAC@@x in total DNN iterations and 2.3x in CPU. Where the first solve has no
   primal deficit (the second regime) the tail is as slow as before: there three sweeps or a guarded
   extrapolation give 1.1 to 1.15x, i.e. closed.

**Settings to run** (tho40_search's block, plus):

```
dnn.beta_sched = w=0.36,wc=2.88
cut.cp_wc_deficit = 0.002      # cpfast.patch; fraction of the cutoff (481 units)
cut.cp_wc_gap = 300            # cpfast.patch
```

Without the patch, `wc=2.88` alone already gives the 2.2x on the subtree (run 1 below), at the price
of slower phases around the path to the optimum, which that subtree does not contain.

## How it was measured

- **Driver**: the worktree's `experiments/dnn_cuts/dnncuts -phase` (added by the patch): a cold first
  solve of 300 iterations with tho40.par's keys (`-p`), then the search's own phase, `dnn_cut_cp`, as
  `latest/bnb.C` calls it (`cut.warmy`, round 0's 1,000 triangles, 1,000 more every 100 iterations),
  no stall rule, cutoff at 241,500 so that nothing stops at a fathom; every check (10 iterations) is
  certified. `scripts/run.sh NAME SPEC options` writes `logs/NAME.log` (first line: the command); the
  first solve's state is cached in `state/` (not kept: 30 MB each). The driver's baseline reproduces
  the search: node A passes the cutoff at phase iteration 3,340 (`logs/A_base.log`), the search's phase
  fathoms at 3,341 (`../tho40_search/runs/w38_z030`), and the two curves agree within 60 units
  from iteration 50 and within 25 from 200 (2,222 / 1,982 / 1,624 / 1,373 / 1,037 / 584 short at 50 /
  100 / 200 / 300 / 500 / 1,000 in the search; 2,182 / 1,926 / 1,601 / 1,363 / 1,034 / 584 in the
  driver). What the driver does not reproduce is the first solve: cold at 0.21 m/3 here, warm at
  `w` in the search, so the start of a phase at `wc` (section 2) and the deficit (section 3) are
  read from search runs.
- **Search** (`scripts/search.sh`, as `../tho40_search/run.sh` with the worktree's `bnbdrv`): the
  depth-2 root's state (`search/wroot_*`: first solve 231,879.8596, as production), cut to one depth-3
  child, resumed with `-r`. `scripts/hist.py` and `scripts/early.py` read the DNN traces.
- **Small instances** (`scripts/small_sweep.sh`, `logs/small_sweep.txt`): the same phase on the roots
  of nug12 (3x4 grid), nug15 (3x5), nug20 (4x5) and one depth-3 node of nug20, seconds each.
- Test nodes: A = `26:10,8:0,38:32` (depth 3, m = 37; 2,427 short in the search, 2,593 after the
  driver's cold first solve), B = `26:10,8:0,18:32` (depth 3; 4,474 short), C = `26:10,8:0,18:32,24:7`
  (depth 4, m = 36; 1,812 short), D5 = `26:10,8:0,2:32,6:39,33:24` (depth 5, m = 35; 2,527 short; the
  "plateau" node of `../tho40_nodes/`), the root, and four nodes on the path to the best known
  solution (depths 8, 10, 12; `../tho40_nodes/scripts/path.py`).

## 1. Diagnosis

### The penalty (`logs/A_*.log`, node A, driver) [m]

`python3 scripts/levels.py logs/A_base.log logs/A_wc144.log ...`: phase iterations to reach a share of
the baseline's 3,400-iteration gain (2,598.5), and the cutoff; `wc` is the phase's penalty in units
of 0.21 x m/3 (today 0.36, `wc` unset).

| run | its | end bound | 25% (238571) | 50% (239221) | 75% (239870) | 90% (240260) | 100% (240520) | cutoff (240514) | s/it |
|---|---|---|---|---|---|---|---|---|---|
| today (`wc` unset) | 3400 | 240519.8 | 100 | 340 | 910 | 1830 | 3400 | 3340 | 0.120 |
| `wc=1.44` | 1000 | 240503.4 | 80 (1.2x) | 170 (2.0x) | 350 (2.6x) | 620 (3.0x) | > 1000 | > 1000 | 0.121 |
| `wc=2.88` | 1000 | 240606.4 | 100 (1.0x) | 170 (2.0x) | 330 (2.8x) | 560 (3.3x) | 840 (4.0x) | 830 (4.0x) | 0.121 |
| `wc=5.76` | 1000 | 240501.1 | 150 (0.7x) | 220 (1.5x) | 390 (2.3x) | 650 (2.8x) | > 1000 | > 1000 | 0.119 |
| `wc=2.88`, `cut.cp_sweeps = 3:1` | 1000 | 240652.0 | 90 (1.1x) | 150 (2.3x) | 280 (3.2x) | 460 (4.0x) | 740 (4.6x) | 720 (4.6x) | 0.130 |

- `rhat/shat` (primal over dual residual, the trace's `bsched` lines) is 18.7 / 16.9 / 15.4 at
  iterations 100 / 300 / 1,000 today, 3.9 at `wc=1.44`, 2.0 at `wc=2.88`, 1.0 at `wc=5.76`: it falls
  as 1 / penalty and the best penalty has it near 2.
- The primal objective of the iterate less the certified value at the same check: today -3,371 at
  iteration 10, -3,130 at 100, -2,365 at 500, -1,966 at 800; `wc=1.44`: -2,818 ... -801; **`wc=2.88`:
  -1,905, -713, -402, -392**; `wc=5.76`: -18, +891, +167, +130. The best penalty keeps the primal
  objective a few hundred units below the bound; too high a penalty puts it above.
- The multipliers' mass h'mu: today 1,711 / 4,473 / 9,058 / 11,021 / 11,352 at iterations 100 / 300 /
  1,000 / 2,000 / 3,400; `wc=2.88`: 6,996 / 11,109 / 9,908 at 100 / 300 / 1,000.
- Same picture on small instances (`logs/small_sweep.txt`): with `wc=2.88` the nug12 root reaches
  today's 1,600-iteration value at about 280 iterations, nug15 at about 450, nug20 at about 285, a
  depth-3 node of nug20 at about 580. The best factor there is 1.44 to 5.76; 11.5 is worse.
- Why it was never set: `w` was swept down for the children (wil50, `speed/nov6_admm`), `wc=0.7` was
  added to keep wil50's phases at 0.7 x 0.3 x m/3 = 0.21 m/3, and the nug30 recipe tho40.par copies
  has no `wc`: tho40's phase runs at 0.076 m/3, a third of wil50's and sko42's.

### What it is not

| candidate | test | result |
|---|---|---|
| the pool's growth | pool size along `logs/A_base.log` | 2,000 cuts at iteration 100, 2,271 at 200, 2,312 at 300, 2,528 at 1,000, 2,799 at 3,400: the separations after 200 find 900 to 1,400 violated triangles (max 0.02 to 0.05) that are almost all already in the pool. At the root (`../tho40_nodes/logs/p2cp_root.log`): 2,822 at 300, 3,568 at 5,490. Separating every 50 iterations instead of 100 at `wc=2.88` (`logs/A_wc288_sep50.log`): +45 / +30 / +10 / +4 units at iterations 100 / 200 / 300 / 500: nothing |
| the projection (one Dykstra sweep) | `cut.cp_sweeps = 3` at `wc=2.88` (`logs/A_wc288_s3.log`) | 75% at 280 instead of 330, cutoff at 720 instead of 830 (1.15x), at 0.130 s instead of 0.121 s per iteration: 1.07x in CPU. Exact projection (20 sweeps) on nug12: the same as 3 |
| the half-space step size | over-relaxation 1.9 (`logs/A_om19.log`); a step scaled to the free entries of Y, on which the cut can act (37 to 40% of ||G||^2; `logs/A_count.log`, `logs/A_wc288_free.log`) | 1.9 diverges after 200 iterations (certified values stay valid, the bound stops at 238,862); the free-entry step is 140 to 190 ahead until iteration 200, then unstable (240,083 at 1,000 against 240,606); stalls on nug15 and nug20 too. Both removed from `cpfast.patch` |
| the eigen-tracker | a full `ssyevd` at every iteration (`-e 1`) at `wc=2.88`, 200 iterations (`logs/A_wc288_exact.log`) | 238,682 / 239,119 / 239,433 at iterations 100 / 150 / 200 against 238,647 / 239,068 / 239,405 with the tracker (refresh 120, degree 3): 28 to 51 units, i.e. 5 to 8 iterations, at twice the CPU per iteration (0.23 s against 0.12). Not the tail |
| the multipliers as such | the oracle of the brief: the phase started from the first solve's Y and Z with the final pool and the final multipliers of the 1,000-iteration `wc=2.88` run, whose bound is 240,606 (`logs/A_orc_mu.log`) | no jump: 218,640 at iteration 10, 238,839 at 100, 239,589 at 200, against 238,647 and 239,405 without the oracle. The right cut multipliers are worth 190 units; the multiplier of the semidefinite constraint has to follow them, which takes the same hundreds of iterations. (Pooled children fathom in 20 iterations because they inherit the parent's Z as well.) So the slow variable is the pair, and moving the cut multipliers alone (extrapolation mode 1) loses on every instance tried |
| extrapolating the state (Z and the cut multipliers) along its trajectory every 100 iterations, guarded | `cut.cp_extrap = 100:1:30:2` | small instances: 1.1 to 1.25x at the tuned penalty (`logs/small_sweep.txt`); the multipliers alone (mode 1): always worse; a factor of 2: no better. Node A at `wc=2.88` (`logs/A_wc288_ex.log`): 7 jumps, none reverted, 92 to 123 units ahead from iteration 200; 90% at 470 instead of 560 and the cutoff at 700 instead of 830 (**1.19x, at the same cost per iteration**). Path node, depth 8, at `w`: 1.1x (section 3). With three sweeps as well (`logs/A_wc288_s3_ex.log`): 90% at 420 and the cutoff at 640 against 560 and 830 for `wc=2.88` alone: **1.3x together** at 0.125 s per iteration against 0.121, and 5.2x against today |

## 2. The change of penalty at the start of the phase

In the search the first solve of a child runs at `w`; a phase at `wc` changes the penalty by wc / w at
its first iteration. Shortfall of the running best of node A's phase, in the search, warm
(`python3 scripts/early.py ...`; `search/w38*`, 400 to 500 phase iterations each) [m]:

| phase penalty | 20 | 40 | 50 | 60 | 80 | 100 | 150 | 200 | 300 | 400 | 500 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| today (`w`) | 2,399 | 2,278 | 2,222 | 2,169 | 2,070 | 1,982 | 1,786 | 1,624 | 1,373 | 1,186 | 1,037 |
| `wc=1.44` | 2,427 | 2,335 | 2,211 | 2,100 | 1,905 | 1,737 | 1,379 | 1,111 | 748 | 525 | |
| `wc=2.0` | 2,427 | 2,427 | 2,298 | 2,174 | 1,945 | 1,757 | 1,377 | 1,077 | 695 | 474 | |
| `wc=2.88` | 2,427 | 2,427 | 2,427 | 2,336 | 2,077 | 1,868 | 1,451 | 1,106 | 709 | 473 | 316 |
| `w` for 100 iterations, then 2.88 | 2,399 | 2,278 | 2,222 | 2,169 | 2,070 | 1,982 | 1,982 | 1,612 | 968 | 640 | 425 |
| `w` for 100, then 2.0 | 2,399 | 2,278 | 2,222 | 2,169 | 2,070 | 1,982 | 1,906 | 1,491 | 939 | 630 | 434 |
| 1.44 for 100, then 2.88 | 2,427 | 2,335 | 2,211 | 2,100 | 1,905 | 1,737 | 1,655 | 1,318 | 820 | 553 | |
| 0.72 / 1.44 / 2.88 from iterations 1 / 21 / 41 | 2,427 | 2,422 | 2,422 | 2,422 | 2,316 | 2,081 | 1,583 | 1,225 | 773 | 516 | |

- A phase at `wc` makes no progress for its first 30 (1.44), 40 (2.0) or 50 (2.88) iterations and is
  ahead of today's from iteration 60 to 100, by which time today's phase has gained 250 to 450 units.
  A node less than about 300 short is therefore better off at `w`: that is `cut.cp_wc_gap`.
- From an inherited pool the cost is higher: in run 1 below, pooled children 25, 133, 137 and 307
  short took 150, 183, 194 and 281 phase iterations at `wc=2.88`; at `w` (runs 2 and 3) 49 to 60, 92,
  142 and 275.
- Stepping the penalty up later or in stages only moves the cost (each step stalls the bound for 20
  to 50 iterations and the staged arms end 70 to 150 iterations behind `wc=2.0`). Not kept in the
  patch. Starting the phase from the complementary pair (psd part, certificate multiplier) of the
  first solve's state, or damping the dual steps for the first 10 to 30 iterations, changes nothing
  (nug20, emulated: the dip is the same to 0.1): the dip is the certificate at the new penalty, not a
  transient of the start.

## 3. Where the higher penalty does not help: the primal deficit

Nodes on the path to the best known solution, driver, certified bound (`logs/V_warm_d08*.log`, depth 8,
m = 32, first solve 238,668.7) [m]:

| phase penalty | 50 | 100 | 200 | 300 | 500 | 750 |
|---|---|---|---|---|---|---|
| `wc=0.18` | 239,109.5 | 239,460.0 | 239,799.5 | 239,956.2 | 240,110.0 | 240,207.5 |
| today (0.36) | 239,000.2 | 239,285.1 | 239,638.1 | 239,835.5 | 240,054.2 | 240,196.1 |
| `wc=0.72` | 238,944.1 | 239,168.0 | 239,490.2 | 239,695.4 | 239,945.0 | 240,110.5 |
| 1.95 (chosen by the residual ratio, `wca=3`) | 238,894.6 | 239,087.9 | 239,376.6 | 239,558.5 | 239,758.7 | 239,904.1 |
| `wc=2.88`, 3 sweeps | 238,744.5 | 239,032.9 | 239,331.3 | 239,488.0 | 239,664.0 | 239,796.5 |
| today, 3 sweeps | 239,031.6 | 239,330.8 | 239,698.8 | 239,900.3 | 240,114.2 | 240,249.9 |
| today, `cut.cp_extrap = 100:1:30:2` (7 jumps, none reverted) | 239,000.2 | 239,285.1 | 239,661.2 | 239,868.0 | 240,092.5 | 240,238.6 |

- Here the order is reversed: the lower the penalty the better, up to iteration 500, and the two
  lowest meet near 750. The depth-10 path node says the same (`logs/V_warm_d10_base.log`,
  `logs/V_warm_d10.log`): 240,346.7 / 240,459.5 / 240,507.4 at 300 / 500 / 1,000 today, 240,037.4 /
  240,126.9 / 240,268.6 under `wc=2.88` with three sweeps: today's 230-iteration value at 1,000. The residual ratio does not tell this node from node A (16 to 17 at `w` on both),
  so a rule on that ratio (`wca`, tried, removed) picks 1.95 here and loses.
- What does tell them apart is the primal objective against the bound. On this node the iterate's
  primal objective is 54 to 600 below the certified value at `w` and 500 to 1,000 above it at
  `wc=2.88`; on node A it is 3,371 below at `w`.
- **The first solve's deficit** (its certified bound less its last primal objective), over every
  phase of the search traces (`../tho40_search/runs/*/dnn.txt`, `search/*/dnn.txt`) [m]:

| class | phases | deficit, units | % of the cutoff |
|---|---|---|---|
| unpooled give-ups off the path, 1,100 to 4,500 short (depths 3 to 5) | 7 | 1,850 to 3,453 | 0.77 to 1.44 |
| unpooled near misses, 29 to 580 short (depth 5) | 7 | 834 to 1,556 | 0.35 to 0.65 |
| pooled children, 25 to 1,037 short (depth 4) | 18 | -274 to +264 | -0.11 to +0.11 |
| path nodes, depths 7, 9, 10 | 3 | -536 to -72 | -0.22 to -0.03 |

  The higher penalty wins in the first class (nodes A, B, C, D5 and the root are of it), costs its
  30 to 50 start-up iterations in the second, and gains nothing in the last two: pooled phases 307 and
  479 short took 281 and 750 iterations at `wc=2.88` (run 1), and at `w` 275 (run 3) and, from 520
  short, 880 at `wc=2.0` (run 2); today's pooled phases of that size took 227 to 630.
  `cut.cp_wc_deficit = 0.002` (481 units) sits between the classes.
- In that second regime the tail is not the penalty. On the depth-8 path node three sweeps are 45
  to 65 units ahead at every mark (1.15x in iterations, 1.07x in CPU) and the guarded extrapolation
  23 to 42 (1.1x): under the brief's 1.15x, closed. A penalty below `w` is 55 to 175 ahead until
  iteration 500 (1.2 to 1.4x early, nothing at 750): one node, not pursued.

## 4. The test nodes [m]

`python3 scripts/levels.py BASE.log VARIANT.log` (driver; shares of the baseline's own gain):

**C** (depth 4): reference gain 1783.9 over the first solve 238702.4; cutoff 240514.05

| run | its | end bound | 25% (239148) | 50% (239594) | 75% (240040) | 90% (240308) | 100% (240486) | cutoff (240514) | s/it |
|---|---|---|---|---|---|---|---|---|---|
| today | 900 | 240486.3 | 80 | 210 | 440 | 670 | 900 | > 900 | 0.108 |
| `wc=2.88` | 700 | 240893.6 | 100 (0.8x) | 160 (1.3x) | 230 (1.9x) | 310 (2.2x) | 380 (2.4x) | 390 | 0.104 |

**D5** (depth 5): reference gain 1583.4 over the first solve 237987.5; cutoff 240514.05

| run | its | end bound | 25% (238383) | 50% (238779) | 75% (239175) | 90% (239413) | 100% (239571) | cutoff (240514) | s/it |
|---|---|---|---|---|---|---|---|---|---|
| today | 800 | 239570.9 | 80 | 210 | 430 | 630 | 800 | > 800 | 0.092 |
| `wc=2.88` | 800 | 240298.9 | 90 (0.9x) | 140 (1.5x) | 180 (2.4x) | 230 (2.7x) | 270 (3.0x) | > 800 | 0.086 |
| 2.25, chosen by the residual ratio (`wca=3`, removed) | 800 | 240277.4 | 80 (1.0x) | 130 (1.6x) | 190 (2.3x) | 240 (2.6x) | 280 (2.9x) | > 800 | 0.092 |

- **B** (`26:10,8:0,18:32`, in the search, `../tho40_search/runs/w18_z030` against `search/w18_wc2.88`;
  shortfall of the phase's running best): today 3,954 / 3,554 / 3,268 / 2,880 / 2,522 at 100 / 200 /
  300 / 500 / 800, stalled at 920 with 2,421; `wc=2.88`: 3,714 / 2,880 / 2,412 / 1,991 / 1,702, stalled
  at 1,010 with 1,605. Today's 500-iteration value at 200 (2.5x), its end value at 300 (3.1x). The
  node does not fathom (the phase has added 2,870 when the stall rule stops it, 126 of them in the
  last 250 iterations), and its children start from a pool worth 800 units more.
- **The root** (`logs/R_base.log`, `logs/R_wc288.log`; 300 cold iterations, then the phase): first solve 224,074.5; today's settings 225,255 / 226,067 / 226,631 / 227,319 / 227,826 at 100 / 200 / 300 / 500 / 750 and 227,904 at 800; `wc=2.88`: 226,295 / 227,643 / 228,470 / 229,251 / 229,702 and **229,763 at 800**. Today's 750-iteration value at about 215 (3.5x). `../tho40_nodes/logs/p2cp_root.log` (dnncuts' defaults, a penalty 1.4 times today's) reads 229,669 at 3,000 and 229,885 at 4,000: the 800-iteration value here is its 3,400-iteration value (4x). @@ROOTLONG@@
- **Cost per iteration**: unchanged by the penalty (0.119 to 0.121 s at m = 37 in the A runs, load
  permitting); three sweeps add 7%.

## 5. One subtree end to end [m]

The warm depth-3 subtree `26:10,8:0,18:32` (4,474 short; 1,490 nodes and about 38,100 DNN iterations
[e] in production; `../tho40_search`'s `w18_z030` is the reference). `python3 scripts/hist.py
search/NAME/dnn.txt`:

| run | settings on top of `z030` | DNN nodes | first solves + phases = iterations | CPU s | phases (fathomed) | iterations vs `z030` |
|---|---|---|---|---|---|---|
| `../tho40_search/runs/w18_z030` | (reference) | 74 | 4,474 + 8,561 = 13,035 | 1,437 | 9 (7) | |
| `search/w18_wc2.88` (run 1) | `wc=2.88` | 38 | 3,323 + 2,568 = 5,891 | 625 | 6 (5) | 2.21x |
| `search/w18_wc2_gap400` (run 2) | `wc=2.0`, `cut.cp_wc_gap = 400` | 38 | 3,381 + 2,496 = 5,877 | 622 | 6 (5) | 2.22x |
| `search/w18_rule` (run 3) | `wc=2.88`, `cut.cp_wc_gap = 300`, `cut.cp_wc_deficit = 0.004` | 38 | 3,323 + 3,099 = 6,422 | 698 | 6 (5) | 2.03x |
| `search/w18_rule2` (run 4) | run 3's with `cut.cp_wc_deficit = 0.002`, `cut.cp_wc_test = 50:0.004` | @@R4_ROW@@ |

- The root's phase ends 800 units higher (1,605 short against 2,421), so the 37 children start from
  a better pool: 32 fathom in their pooled first solve (29 before) and the 5 that need a phase are 25
  to 479 short (8 phases from 214 to 1,037 short before, three of them 1,290 to 2,610 iterations
  long, one of which gave up and branched). No node reaches depth 5 (36 did).
- Run 3 against run 1: the pooled near misses (25, 133, 137, 307 short), left at `w` by the rules,
  fathom in 49, 82, 110 and 275 phase iterations instead of 150, 183, 194 and 281 (392 saved). But the
  pooled phase 479 short, which the deficit rule also left at `w` (its first solve's deficit is -180),
  took 1,573 iterations instead of 750 (823 lost): a pooled phase that long does gain from the
  higher penalty, and the first solve's deficit does not see it (a pooled first solve has converged
  on the pool it inherited; the phase then adds cuts whose multipliers creep). Run 4 adds a test on
  the phase's own iterate for that case (`cut.cp_wc_test`). @@R4_TEXT@@
- Against production's 38,100 iterations [e]: 6.5x for run 1. CPU: 2.3x against `z030`, measured hours apart on
  a shared box (0.106 s per iteration in run 1, 0.110 in the reference).
- One subtree, the class that holds 54% of the A/B subtree's CPU in tho40_search's table. Not run: a
  subtree whose depth-3 root fathoms (`38:32`: one node, 1,001 iterations against 3,541, 3.5x, is
  the probe above), the classes above 5,000 short, any Knuth dive.

## 6. Validity [m]

- The certificate is untouched: every value is the existing certificate at that check's Z and cut
  multipliers; the penalty, the two rules and the extrapolation only move the iteration.
- **Path to the best known solution** (`python3 scripts/validity.py logs/V_*.log`; nodes of
  `../tho40_nodes/` part 1; optimum 240,516 if the best known value is optimal; every check's
  certified value, cutoff at 241,500 so that nothing stops):

@@VALID@@

- Every search run above ended with its subtree fathomed and no `New objective value` line, as the
  reference did (the reference and production fathom the same subtree without a better solution).
- Small instances: no certified value above the optimum in `logs/small_sweep.txt` (nug12 578: at
  most 574.35; nug15 1,150: 1,144.64; nug20 2,570: 2,527.68), the unstable arms included.
- Not rerun: `tests/dnn_test.py`, the brute-force validator of `../tho40_nodes` (the CP path's
  certificate code is not changed by the patch; with every new key at its default the lean build
  reproduces the unpatched phase to the last digit on nug12, 571.194178 at 300 iterations).

## 7. The patch and the keys

`cpfast.patch` (against 8dd9afe1; `git apply`, then `make -C latest && make -C latest dnntest dnncuts`).
Every key defaults to today's behaviour.

| key | default | what |
|---|---|---|
| `dnn.beta_sched = ...,wc=F` | unset | existing. The phase's penalty factor. 2.88 for tho40 (1.44 to 2.88 all give at least 3x at 90% on node A) |
| `cut.cp_wc_deficit = D` | 0 (off) | new. The phase takes `wc` only when the first solve's bound less its primal objective is at least D x the cutoff; otherwise `w`. 0.002 |
| `cut.cp_wc_gap = G` | 0 (off) | new. A phase that starts within G of the cutoff keeps `w`. 300 |
| `cut.cp_sweeps = 3:1` | 1:1 | existing. 1.15x fewer iterations at 7% more per iteration in both regimes: optional, not in the subtree runs |
| `cut.cp_extrap = E:GAMMA[:GUARD[:MODE]]` | unset | new, experimental. Every E iterations Z and the cut multipliers (mode 2; 3: Y too; 1: the multipliers alone) move GAMMA times their change over the last E/2; if the certified value GUARD iterations later is below the value before the jump, the state is restored and GAMMA halved. 1.1x: not recommended as it stands |
| `dnncuts -p PAR -phase [-savewarm F | -loadwarm F] [-firstbeta B] [-loadpool F [-poolscale S]] [-cpextrap ..] [-cpwcgap G] [-cpwcdef D]` | | the driver used here: tho40.par's `dnn.*` keys, then `dnn_cut_cp` as the search calls it |

`cpfast_experiments.patch` is the worktree before the clean-up: it also has the options that were
measured and dropped (`cut.cp_free`, `cut.cp_omega`, `dnn.beta_sched`'s `wcs=` steps and `wca=`), and
is what the logs that use them were made with (`-cpfree`, `-cpomega` through `-phase`, `par/r*.par`,
`par/pr_*.par`, `par/wca*.par`).

## 8. What I would try next

1. **A/B the three keys on the unmeasured classes** before any estimate moves: `2:32` (25% of the A/B
   subtree), one depth-2 subtree with a far root, and the path (`../tho40_search`'s `path7`, `path9`
   jobs) with the rule on, to see that the rule leaves the path's phases as they are and what the
   phases of the path's off-path siblings gain.
2. **The penalty from the deficit, continuously.** The data here are two clusters; a factor like
   2^(deficit / 0.4%) clamped to [0.25, 8] would also lower the penalty where the primal is ahead (the
   depth-8 path node gains 1.2 to 1.4x early from `wc=0.18`), and would let the first solve itself use
   it: a far miss's first solve at `w` is where the deficit comes from.
3. **The first solve of a node that will get a phase.** With the phase 3 to 4 times shorter the 200 to
   370 iterations before it are now 55% of run 1's iterations. tho40_search's `dnn.gu_zone = 0` and an
   earlier hand-over belong here.
4. **The second regime's tail** (pooled children, path nodes): not the penalty, not the projection.
   The stop-rule work (stop and equip the children) is the lever there.
5. The stall rule now stops B at 1,010 iterations with 126 gained in the last 250: with the faster
   phase a shorter window (`cut.stall`, `cut.iters`) is worth a replay on the new traces.

## Files

| | |
|---|---|
| `cpfast.patch`, `cpfast_experiments.patch` | section 7 |
| `scripts/run.sh`, `scripts/queue.sh`, `scripts/jobs.txt` | one driver run; the sequential queue and what was left in it; `logs/queue.log` has the order and times |
| `scripts/search.sh` | one bnbdrv run (warm subtree or a `-d 3` probe) |
| `scripts/curve.py`, `scripts/levels.py`, `scripts/validity.py` | driver logs: bound by iteration, iterations to levels, the path check |
| `scripts/hist.py`, `scripts/early.py` | DNN traces of the search: phases and their running best |
| `scripts/small.sh`, `scripts/small_sweep.sh`, `logs/small_sweep.txt` | the small instances |
| `par/` | `wcF.par`: tho40.par with `wc=F`; `z030_*.par`: tho40_search's `z030` plus the keys of section 5; `pr_*.par`: the 400 to 500 iteration probes of section 2 |
| `logs/A_*`, `B`: see section 4, `C_*`, `D5_*`, `R_*`, `V_*` | driver runs by node |
| `search/` | `wroot_*` (the depth-2 root), `w38_wc2.88`, `w38p_*` (probes), `w18_*` (the subtree runs): `cmd.txt`, `log.txt`, `nodes.txt`, `dnn.txt` |
