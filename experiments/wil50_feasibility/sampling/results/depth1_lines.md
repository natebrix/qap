# Depth-1 lines: rule 9's line against three others at the two weakest depth-1 nodes (Nov 9)

**The question.** Stage A's tree (root line row 21) has two weak depth-1 nodes: r.1 = `21:1` (DNN bound 48,273.09, 541 short) and r.5 = `21:10` (48,244.00, 570 short). Rule 9 branched both on col 0, and stage A fathomed 13 and 11 of their 49 children. Would another line give clearly better depth-2 children? This is anecdotal evidence from two nodes, with each line's 49 children bounded once. It is not an estimate of the full tree.

**Verdict: no clear winner at either node.** The rule set in advance: a line wins only if it beats rule 9's by at least 1.5x on the subtree-type scores and does not lose on the fathomed count.

- **r.1 (`21:1`): close, but no.** Row 45 beats rule 9's col 0 on every subtree-type score:
  - lines cost 1.98x;
  - expected subtree 1.85x (95%: 1.08-2.62) with the step model and 1.64x (1.13-2.22) with the log-linear model;
  - 13 children beyond the 244 zone against 19.

  It fails the fathomed condition: no child fathoms in its first solve, against 2 for col 0. The lower ends of the 1.5x intervals are 1.08 and 1.13.
- **r.5 (`21:10`): no.** Row 45 is 1.83x better on the lines cost, but only 1.17x (1.08-1.40) and 1.33x (1.24-1.54) on the expected subtree. Its fathomed count is 0 against 2.
- **The other alternatives are worse than rule 9's line or no better.** These are rule 8's lines ranked 3 and 4: col 49 and row 25 at r.1, col 40 and col 49 at r.5.

Row 45 is the best alternative at both nodes. It is rule 8's top line, and it is also the row that rule 9's own lookahead considered and rejected. That lookahead scored the two candidates almost equally: col 0 led by 0.023% at r.1 and 0.013% at r.5. Measured by the DNN first solves, the choice between them is not close (lines cost about 2x).

Row 45's children start nearer the cutoff: median 164 short against 203 at r.1, 166 against 231 at r.5. But none of them ends at or below 0 in the first solve. The fathomed-count condition rests on 2 children per node at -3.5/-0.7 and -0.4/-0.2. Stage A's phase outcomes model 21 against 16 fathoms after the cut phase (r.1) and 20 against 11.5 (r.5); see the table and caveats. The optional phase pass would settle this, but it did not fit the budget (below).

## Method

1. **Binary.** HEAD `fe31484`, built with `make solver` in this worktree. The binary's sha256 is `de4e0697b8cb...`, and `latest/` was unmodified when it was built. The run uses a private copy, with the settings of `../fullrun/settings.sh` and `../fullrun/fullrun.par`, the same as stage A: instance `wil50t.dat`, `wil50t.sym`, `-u 48816` (cutoff 48,814.05), one thread, `taskset`.

   `fullrun_check` fails at HEAD. `settings.golden` lacks `QAP_DNN_WARM_DIR`, a run setting that 1d40615 added and that sits outside the settings hash. The dump otherwise matches, hash `ebd9848cb3f7e0d2` included. `depth1_lines.fullrun_env_head` accepts exactly that difference and records the diff in `depth1_lines/config.json`.
2. **The root under HEAD.** `bnbdrv -p fullrun.par -d 0`: QPB bound 46,131.508 (a667727: 46,100.957), line row 21, 15 children. Child 1 is `21:1` and child 5 is `21:10`, numbered as in stage A.
3. **r.1 and r.5 from HEAD's root state**, as stage A runs a task: `write_one(root/out.state, K)`, then `bnbdrv -r in.state -d 1 -o out.state`. The node gets a cold DNN with 300 iterations. Its children carry its warm start: the root's QPB point and the node's own fp32 DNN state.
   - The DNN bounds are 48,273.0884 and 48,243.9954, equal to a667727's to the printed digit (trace lines below).
   - **Rule 9 still picks col 0 at both.**
4. **The diagnostic build** applies `depth1_lines/xbar_debug2.patch`, which is not committed to `latest/`. Under `QAP_XBAR_DEBUG=2` it prints every free row's and column's rule-8 score at a rule-9 node, summed as `dnn_xbar_lines` sums it. Both nodes ran on both binaries:
   - the traces are equal except pid and CPU;
   - the child sets are equal;
   - the warm-start line, DNN state included, is byte-equal (sha256).
5. **Candidates**: rule 9's line, plus the next 3 lines by rule 8's score, with rows and columns ranked together and rule 9's own line left out. At both nodes rule 9's col 0 ranks 2nd, behind row 45.
   - r.1: col 0, row 45 (0.2019), col 49 (0.1520), row 25 (0.1517). Col 9 (0.1512) ranks 5th and was not run.
   - r.5: col 0, row 45 (0.2436), col 40 (0.1673), col 49 (0.1471). Row 0 (0.1470) ranks 5th.
6. **Children by one rewrite for every line** (`depth1_lines.rewrite_line`, generalising `stage_a.force_root_line`). This includes rule 9's own line, so the comparison is paired.
   - A row-f line has f at each of the 49 free locations; a col-c line has each of the 49 free facilities at c.
   - Every other field of the solver's child is kept: depth 2, the inherited bound (the node's DNN bound), and the node's warm start, which `warm_child_offsets` cuts at the child's own assignment.
   - **Predicted bound = the inherited bound for every line**: 48,273.0884 at r.1 and 48,243.9954 at r.5.
   - "prod" is rule 9's children straight from the solver. These carry production's predicted bounds: QPB + lapU + U_ij, 46,464-46,960. They are the rewrite check.
7. **First DNN solve only.** Each child runs `bnbdrv --allow-settings-change -p child.par -r in.state -d 2 -o out.state`, where `child.par` is `fullrun.par` plus `cut.maxdepth = 1`.
   - The solve starts warm from the node's fp32 DNN state and uses the production give-up rule: stop at 300 iterations when hopeless.
   - There is no cut phase.
   - A miss branches by rule 9 at depth 2 (its line is recorded), but no child of it is bounded.
   - Why the flag: HEAD's states carry a settings record (QAP-STATE 9), and `cut.maxdepth` differs from it. Every child's log warning lists `cut.maxdepth` and `par.hash` as the only differences (checked by the score script).
   - Each child is one task pinned to one of CPUs 0-3, and all 490 tasks run in a seeded random order. The big states are deleted after each task.
8. **Scores per line** (`depth1_lines_score.py`). "Gap" means first-solve shortfall.
   - **fathomed**: the first-solve fathoms.
   - **shortfall histogram**: stage A's bins.
   - **beyond 244**: children past the zone `cut.zone` x 48,816 = 244.08, which the phase cannot rescue.
   - **lines cost**: sum of 48^(max(gap,0)/130). It is exact over the 49 children, so it has no interval. "Top-3 share" is how much of it the 3 worst children carry.
   - **expected subtree E**: the sum over children of E(s), the expected total anchor-model size of a depth-2 node's depth-3 children given its first-solve shortfall s. It is fitted from stages A and B (`depth1_lines/model.md`):
     - E(s) = 0 when the first solve fathoms; otherwise E(s) = q(s) x 48 x m(s).
     - q(s) is the fraction of stage A's 635 depth-2 nodes in s's bin, binned by first-solve shortfall, that the production cut phase did not fathom.
     - m(s) is the mean `tree_estimate.size_anchor` of stage B's 400 sampled depth-3 nodes whose parents' first-solve shortfall is in the bin. It uses the converged shortfall where stage B converged the node and the give-up shortfall otherwise; a fathomed node has size 1.
     - Two variants: step bins, and log-linear in s between the bins' mean shortfalls.
     - The 95% bootstrap resamples stage A's 635 and stage B's 400 (4,000 times). The line's own 49 children are the whole line, not a sample.
   - **fathomed after a phase, modelled**: the first-solve fathoms, plus for each child in the zone the fraction of stage A's nodes in its bin that the phase fathomed. This is a model, not a measurement.

**Is the expected-subtree score sound?** It is sound enough to rank lines, with these limits:
- It assumes a depth-2 node's subtree depends on the node only through its first-solve shortfall.
- It inherits the anchor model's error.
- Its last bin is open, starting at 400 short. Stage B has 14 children of parents 400-458 short, with a mean of 440, and none beyond 458. Every child 400 or more short counts about 21,000, so lines with several such children are scored coarsely.
- Stages A and B ran on a667727. HEAD reproduces a667727's depth-2 first solves bit for bit on rule 9's 98 children (below), so the mapping from s to outcome carries over.
- Below 100 short every stage A node fathomed after its phase, so E = 0 there. In 100-150, q = 0.248.

## Commands

```
make solver                                   # HEAD fe31484 -> latest/bnbdrv (sha256 de4e0697...)
git apply experiments/wil50_feasibility/sampling/results/depth1_lines/xbar_debug2.patch && make solver
cp latest/bnbdrv <diag>/bnbdrv && git checkout latest/branch.C latest/settings.C && make solver
cd experiments/wil50_feasibility/sampling
python3 depth1_lines.py prep --binary ../../../latest/bnbdrv --diag-binary <diag>/bnbdrv   # root, r.1, r.5 x 2 binaries
python3 depth1_lines.py kids --max-wall-h 5 [--stop-after-min 7.3]   # 490 child tasks; resumable
python3 depth1_lines_score.py                 # results/depth1_lines/{tables,model}.md, summary.json, traces/, logs/
```

The background-job limit of the tool environment stopped the driver twice, after 30 minutes and after 2 hours. Each time its 4 running solvers died (exit -15, in `depth1_lines/launcher.log`); those 8 tasks were rerun from scratch. The last 171 tasks ran as foreground slices of 7.3 minutes (`--stop-after-min`).

## Results

Setup, from `depth1_lines/tables.md`:

| node | DNN bound, HEAD (= a667727) | shortfall | CPU s | rule 9's line | rule 9's lookahead score, row 45 / col 0 | rank of col 0 by rule 8's score |
|---|---|---|---|---|---|---|
| r.1 `21:1` | 48,273.0884 | 541.0 | 114.8 | col 0 | 114,558,810.7 / 114,585,133.4 | 2 (row 45 first) |
| r.5 `21:10` | 48,243.9954 | 570.1 | 113.4 | col 0 | 114,599,870.2 / 114,614,395.8 | 2 (row 45 first) |

Per node and line: 49 depth-2 children each, first DNN solve, no cut phase. Every line's children have predicted bound = inherited. Ratios are rule 9's score divided by the line's, so above 1 means the line is better. The full table, with the log-linear E values and per-line CPU, is in `depth1_lines/tables.md`.

| node | line (rule 8 rank) | fathomed | in the zone, not fathomed | beyond 244 | fathomed after a phase, modelled (95%) | shortfall min / median / max | lines cost | rule 9 / line | top-3 share | E, step (95%) | rule 9 / line (95%) | rule 9 / line, log-linear E (95%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r.1 | **col 0, rule 9** (2) | **2** | 28 | 19 | 16.0 (15.5-16.6) | 10.1 / 203.0 / 430.9 | 1.10e6 | 1 | 0.68 | 1.21e5 (5.7e4-1.9e5) | 1 | 1 |
| r.1 | row 45 (1) | 0 | 36 | **13** | **20.9** (20.2-21.7) | 31.5 / **164.1** / 407.0 | 5.54e5 | **1.98** | 0.67 | 6.51e4 (4.0e4-9.8e4) | **1.85** (1.08-2.62) | **1.64** (1.13-2.22) |
| r.1 | col 49 (3) | 1 | 25 | 23 | 11.7 (11.3-12.1) | 19.2 / 243.3 / 420.0 | 1.21e6 | 0.91 | 0.58 | 1.37e5 (8.2e4-2.0e5) | 0.88 (0.61-1.08) | 0.74 (0.58-0.82) |
| r.1 | row 25 (4) | 0 | 0 | 49 | 0.0 | 258.0 / 306.6 / 398.6 | 9.07e5 | 1.21 | 0.32 | 1.51e5 (9.1e4-2.2e5) | 0.80 (0.37-1.50) | 0.67 (0.38-1.03) |
| r.1 | col 0, prod (check) | 2 | 28 | 19 | 16.0 | 10.1 / 203.0 / 430.9 | 1.10e6 | 1.00 | 0.68 | 1.21e5 | 1.00 | 1.00 |
| r.5 | **col 0, rule 9** (2) | **2** | 25 | 22 | 11.5 (11.1-12.0) | 42.7 / 230.6 / 457.8 | 2.35e6 | 1 | 0.71 | 1.40e5 (7.5e4-2.2e5) | 1 | 1 |
| r.5 | row 45 (1) | 0 | 34 | **15** | **20.2** (19.3-21.1) | 52.7 / **166.1** / 437.1 | 1.29e6 | **1.83** | 0.67 | 1.19e5 (5.6e4-1.9e5) | 1.17 (1.08-1.40) | 1.33 (1.24-1.54) |
| r.5 | col 40 (3) | 0 | 20 | 29 | 5.8 (5.4-6.2) | 36.8 / 265.7 / 454.0 | 2.03e6 | 1.16 | 0.65 | 1.64e5 (9.6e4-2.4e5) | 0.85 (0.75-0.92) | 0.91 (0.80-0.99) |
| r.5 | col 49 (4) | 1 | 19 | 29 | 8.1 (7.6-8.5) | 46.1 / 271.2 / 449.2 | 2.83e6 | 0.83 | 0.56 | 2.12e5 (1.1e5-3.3e5) | 0.66 (0.61-0.70) | 0.72 (0.63-0.87) |
| r.5 | col 0, prod (check) | 2 | 25 | 22 | 11.5 | 42.7 / 230.6 / 457.8 | 2.35e6 | 1.00 | 0.71 | 1.40e5 | 1.00 | 1.00 |

Shortfall histogram (first solve; stage A's bins; the zone ends at 244.08):

| node | line | fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | >= 500 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r.1 | col 0 (rule 9) | 2 | 2 | 7 | 6 | 7 | 6 | 9 | 5 | 1 | 4 | 0 | 0 |
| r.1 | row 45 | 0 | 2 | 13 | 7 | 9 | 5 | 4 | 4 | 4 | 1 | 0 | 0 |
| r.1 | col 49 | 1 | 3 | 5 | 3 | 6 | 8 | 6 | 9 | 5 | 3 | 0 | 0 |
| r.1 | row 25 | 0 | 0 | 0 | 0 | 0 | 0 | 21 | 21 | 7 | 0 | 0 | 0 |
| r.5 | col 0 (rule 9) | 2 | 2 | 4 | 4 | 7 | 8 | 8 | 7 | 3 | 3 | 1 | 0 |
| r.5 | row 45 | 0 | 0 | 12 | 10 | 9 | 3 | 5 | 3 | 3 | 4 | 0 | 0 |
| r.5 | col 40 | 0 | 1 | 2 | 3 | 7 | 7 | 10 | 12 | 3 | 3 | 1 | 0 |
| r.5 | col 49 | 1 | 1 | 2 | 5 | 4 | 7 | 7 | 9 | 7 | 6 | 0 | 0 |

The expected-subtree model, from `depth1_lines/model.md`:

| first-solve shortfall | stage A nodes | q, not fathomed after the phase | stage B children | m, mean anchor size | E at the bin's mean shortfall (95%) |
|---|---|---|---|---|---|
| 0-100 | 176 | 0.000 | 0 | - | 0 |
| 100-150 | 121 | 0.248 | 51 | 1.00 | 11.9 (8.3-15.7) |
| 150-200 | 96 | 0.927 | 122 | 1.00 | 44.7 (41.9-47.0) |
| 200-244.08 | 70 | 1.000 | 81 | 1.28 | 61.2 (50.3-75.7) |
| 244-300 | 54 | 1.000 | 72 | 16.97 | 814 (430-1,310) |
| 300-350 | 32 | 1.000 | 34 | 94.4 | 4,530 (2,310-7,150) |
| 350-400 | 15 | 1.000 | 26 | 113.9 | 5,470 (1,520-12,100) |
| 400- | 14 | 1.000 | 14 | 439.7 | 21,100 (5,590-39,100) |

**Check of the modelled fathom count against measured outcomes.** Stage A ran production's cut phase on exactly rule 9's children, and its first solves are identical to HEAD's. It fathomed 13 of r.1's children and 11 of r.5's. The model gives 16.0 and 11.5.

## The rewrite check

Rule 9's 49 children, rewritten (predicted bound = inherited) against straight from the solver (production predicted bounds 46,464-46,960), child for child:

| node | first lb, primal, iterations, status identical | largest abs. lb difference | same depth-2 line (rule 9 at depth 2) |
|---|---|---|---|
| r.1 | 49 of 49 | 0.0000 | 49 of 49 |
| r.5 | 49 of 49 | 0.0000 | 49 of 49 |

They are identical, not merely close. The predicted bound enters `pick_param` only as `max(node.bound, node.predicted_bound)`. Production's predicted bounds are about 1,300 below the inherited DNN bound, so both versions select the same parameter block. They also enter the same `max()` after QPB, and the DNN solve does not read them.

## HEAD against stage A on the same children

HEAD's first solves of rule 9's 98 children equal stage A's (a667727) in first lb, primal, iterations and status: 98 of 98, to the trace's 4 decimals. The depth-1 traces are equal too:

```
HEAD:    pid 2633 depth 1 fixes 21:1 | first lb 48273.0884 primal 47946.1493 status hopeless iters 300 cpu 114.54 warm 0 ...
a667727: pid 25284 depth 1 fixes 21:1 | first lb 48273.0884 primal 47946.1493 status hopeless iters 300 cpu 117.56 warm 0 ...
```

The root's QPB bound changed by 30.6, but the depth-1 DNN is cold and the depth-2 DNN is warm from the parent's DNN state, so the change does not reach them.

It does reach the depth-2 branching. On the 41 children that stage A bounded without a phase (so the marginals are the same), HEAD's rule 9 picks a different depth-2 line at 4: r.1 `4:0`, `19:0`, `47:0` (stage A row 45, HEAD col 40 / col 9 / col 9) and r.5 `5:0` (col 40 against row 45). Two things could cause this:
- At r.1 the children's relative gap is (48,816 - 48,273.09) / (48,816 - root). It is 0.2022 under HEAD's root bound and 0.19996 under a667727's, so they straddle the first block's 0.20 threshold. The blocks differ in the lookahead's FW iterations, 25 against 12.
- The QPB code itself changed (18ebd58).

Which of the two it is was not separated. This does not affect the first-solve measurements here. It does mean HEAD's depth-3 frontier would not be exactly stage A's.

## CPU and wall

- **Child tasks:** 490 tasks, 53,529 CPU s (first solves 53,377 s). That is about 5,200-5,600 CPU s per line of 49 children, and 109 s per task on average. Row 45's children take a few more iterations on average: 307-308, against 293-298 for col 0.
- **Prep:** the root plus 4 depth-1 runs, 455 CPU s.
- **Wall:** 4.06 h on 4 CPUs (17:54-21:58), about 8% of it idle in drains and restarts. This is over the 3 h budget for the first pass but under the 5 h cut-off, so no line was dropped.
- **The optional second pass** (the cut phase on the in-zone children of rule 9's line and row 45) was not run, because the first pass did not fit the budget. At stage A's 290 s per phase plus the 110 s first solve, it would take about 400 CPU s per child for 64 children (r.1) or 123 (both nodes). That is 1.8 h or 3.4 h of wall on 4 CPUs. It would replace the modelled fathom counts with measured ones, and it is the obvious next step for row 45 at r.1.

## Caveats

- **Two nodes, one run per child.** These are the two weakest depth-1 nodes, picked for being weak. The lines are compared on the same node, so the comparison is paired at the node, but not child for child: the lines have different children.
- **The first solve only.** The children's fate in the zone (0-244 short) is not measured. It comes from stage A's phase outcomes by bin, and the "fathomed after a phase, modelled" column is a model. On rule 9's line it overpredicts r.1's measured 13 by 3.
- **The subtree-type scores are concentrated.** The 3 worst children carry 56-71% of each line's lines cost (32% for row 25, whose children are all 258-399 short). Rule 9's col 0 at r.1 has 4 children 400-431 short, against row 45's 1 (407). Most of the 2x is those few children.
- **The expected-subtree score is coarse at the far end.** Every child 400 or more short counts about 21,000, with a 95% interval of 5,600-39,000. Its 1.5x margins are not robust within the model's own bootstrap: the lower ends are 1.08-1.24 for row 45.
- **No QPB elimination for the alternative lines.** Production drops a child at branching when its QPB-predicted bound reaches the incumbent. At depth 2 those bounds are about 1,850 below the cutoff, so it would drop none, as for rule 9's line (`e_level` = 0 in stage A).
- **Rule 8's ranking beyond its top two lines is weak here.** Rank 3 and 4 are no better than col 0 at either node. Row 25 at r.1 (0.1517, level with col 49's 0.1520) puts all 49 children 258-399 short, so none is in the zone.

## Surprises

1. **Rule 9's choice is a near tie in its own score, but not in the DNN bounds.** The diagnostic build's logs (`depth1_lines/depth1/r.*_bin_diag_solver.log`):
   ```
   r.1: xbar9: depth 1 col 0 (row 114558810.7 col 114585133.4)
   r.5: xbar9: depth 1 col 0 (row 114599870.2 col 114614395.8)
   ```
   Col 0 wins by 0.023% and 0.013% of the lookahead score, which is dominated by 48 x the children's QPB bounds. Row 45, which rule 8's score ranks first, has the lower lines cost by 1.8-2.0x and about 32% fewer children beyond the zone (13 against 19, 15 against 22).
2. **HEAD and a667727 give identical DNN bounds at depths 1 and 2 on these nodes**, despite the 30.6-unit change in the root's QPB bound. See the trace lines above.
3. **`fullrun_check` fails at HEAD** because of the stale golden:
   ```
   < settings 65 hash ebd9848cb3f7e0d2
   > settings 64 hash ebd9848cb3f7e0d2
   < setting QAP_DNN_WARM_DIR "" default run
   ```
4. **`common.write_one` rejects HEAD's state files.** It accepts versions 6-8, and HEAD writes `QAP-STATE 9`. `stage_a.py` and `stage_b.py` would fail on HEAD's binary as they stand. `depth1_lines.py` reads version 9 itself and leaves `common.py` unchanged. In addition, a child resumed with different settings needs `--allow-settings-change`.

## Files

- `../depth1_lines.py`: the driver (prep, kids).
- `../depth1_lines_score.py`: the scores.
- `depth1_lines/`:
  - `tables.md` and `model.md`;
  - `summary.json` (every child: spec, line, shortfall, status, iterations, CPU, depth-2 line);
  - `traces/<node>_<line>.txt` (the DNN trace lines) and `logs/<node>_<line>.log.gz` (the solver logs);
  - `depth1/` (the root's log and the depth-1 runs on both binaries);
  - `prep.json`, `config.json`, `child.par`, `launcher.log`;
  - `xbar_debug2.patch` (the diagnostic print).
