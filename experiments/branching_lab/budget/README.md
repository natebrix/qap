# What the lookahead costs at the top of each tree, as a share of the run

Nov 21 2026. No solver runs: `budget.py` over figures already in the repository. [m] = read from a run's or a
stage's files, [e] = an estimate or derived here.

**On the three open grids the lookahead at every depth-1 node costs 0.1 to 1.2% of the run; at every depth-2
node 4 to 21%; at depth 3 more than the run. So: always at depth 1, at depth 2 only on the nodes with the
largest shortfall, never below.**

## The cost of one firing

`QAP_LA=4:150` at a node with m free facilities: 4 rows and 4 columns by marginal score, their children bounded
with 150 warm ADMM iterations each (no certificate): 8m - 16 solves at dimension m - 1. At m = 49 that is
56,400 child iterations: the work of about 190 first solves at give-up 300. The table charges every child all
150 iterations and both sides at every node, so it is an upper bound: at tho30's depth 1 the children averaged
110 to 137 iterations, and nodes with symmetry left score rows only [m `../grid/runs/tho30_la150_d1/la.log`].
The seconds per iteration are each estimate's own, on its own machine, at the children's depth.

## The table (`python3 budget.py`; `budget.py 2` for 2 rows and 2 columns, half the cost)

| instance | firing depth | branching nodes | core-hours a firing | core-hours, this depth | cumulative | share of the run, cumulative |
|---|---|---|---|---|---|---|
| tho40 (settings c) | 1 | 12 | 2.6 | 31 | 31 | 0.5-0.6% |
| | 2 | 403 | 2.1 | 842 | 873 | 14-17% |
| | 3 | about 6,500 | 1.8 | 11,500 | 12,300 | about 2 runs |
| | | | | | | of 5,182-6,255 c7a core-hours (estimate of Nov 19) |
| wil50 | 1 | 15 | 5.8 | 87 | 87 | 0.7-1.2% |
| | 2 | 304 | 4.9 | 1,490 | 1,580 | 12-21% |
| | 3 | about 3,650 | 4.3 | 15,700 | 17,300 | 1.3 to 2.3 runs |
| | | | | | | of 7,500-13,200 cloud Xeon core-hours (fit of Nov 11) |
| sko49 | 1 | 10 | 5.9 | 59 | 59 | 0.1-0.2% |
| | 2 | 313 | 5.2 | 1,630 | 1,690 | 4.1-4.6% |
| | 3 | about 10,800 | 4.5 | 48,800 | 50,500 | 1.2 to 1.4 runs |
| | | | | | | of 37,000-41,000 cloud Xeon core-hours (estimate Cw of Nov 13; 95%: 17,000-65,000) |

The two solved tai-b trees as back-tests (the root has marginals there, so depth 0 fires):

| instance | firing depth | branching nodes | core-hours a firing | cumulative | share of the run, cumulative |
|---|---|---|---|---|---|
| tai50b (`tai50b-r9`, 116 c7a core-hours) | 0 | 1 | 11 | 11 | 10% |
| | 1 | 10 | 10 | 110 | 94% |
| tai60b (`tai60b-261117`, 3,597 m7a core-hours) | 0 | 1 | 31 | 31 | 0.8% |
| | 1 | 5 | 27 | 167 | 4.7% |
| | 2 | 4 | 25 | 266 | 7.4% |
| | 3 | 12 | 22 | 530 | 15% |

## What it says

- **Depth 1 is always affordable on the open grids**: 12, 15 and 10 firings, 31, 87 and 59 core-hours, at most
  1.2% of the run. It pays if it saves more than that; on tho30 the same ten firings saved 25% of the nodes
  [m `../grid`], one grid at n = 30.
- **Depth 2 wholesale is not "underuse"**: 4% (sko49) to 21% (wil50) of the run for 300 to 400 firings.
- **But the work is concentrated, so depth 2 can be selective by shortfall** [e, the estimates' own figures]:
  tho40's two largest depth-1 nodes hold 41 to 53% of the estimated CPU; at wil50 the 1% of depth-3 nodes more
  than 300 short hold 44%; at sko49 the 9 of 250 sampled depth-3 nodes beyond 330 short hold 51 to 68%. The
  tenth of the depth-2 nodes with the largest shortfall is 30 to 40 firings: 0.4 to 2% of the run.
- **Depth 3 and below: never wholesale.** One level of firings costs more than the run.
- **The per-node form of the rule**: fire where the firing costs at most a fraction f of the node's predicted
  subtree, the prediction being the estimate's size curve in the shortfall (wil50, depth 3: ln(DNN nodes) =
  1.742 + 0.0205 s, 22.5 s a node [m `../../estimates/wil50/refit.txt`]). At f = 5% a wil50 depth-3 node would
  have to be 380 short (none sampled is beyond 343); at f = 25%, tho30's gain, 301 short: the 1% tail. So
  the rule underuses by construction, and depth and gap enter through one number.
- **tai50b: the tree is too small.** The root firing alone is a tenth of the 116 core-hour tree, and it chose
  rule 9's line (`../lsla/README.md`). The same firing would be noise on a run a hundred times larger.
- **tai60b: the top is thin.** Through depth 3 there are 22 branching nodes in a tree of 331,149; firing at
  all of them is 15% of the run (7.5% with 2 rows and 2 columns), at the root and depth 1 4.7%.

## A reasonable rule, biased to underuse

1. Depth 0 and 1: always (the root where it has marginals).
2. Depth 2: only where the node's shortfall is in the top tenth of its depth (the threshold read off the
   estimate's stage A or sweep, before the run), or by the per-node form with f about 5%.
3. Below: rule 9.
4. A cap on the lookahead's total spend (2% of the estimate) as the backstop.

## Not covered

- **The grids' roots have no marginals** (QPB only at depth 0; the root's line is a QPB row lookahead under
  symmetry [m stage A summaries, `tho40_dives/logs/root/nodes.txt`]). A root lookahead there would be cheap,
  since a row has only 10 to 15 children after symmetry (4 rows: about 1 core-hour [e]), but it needs a DNN
  bound at depth 0 or a shortlist from another score. Not implemented.
- The lookahead is in the experiment build only (`../lsla/lsla.patch`), not in cognac.
- No whole tree at n >= 40 with the lookahead. The gain side of the ledger is tho30 (25%), tai35b (23% fewer
  nodes, more CPU) and one level at tai60b's chain nodes.

## Sources and how each figure was got

- tho40: 12 depth-1 nodes, 403 of 404 depth-2 nodes branching [m `../../tho40b_analysis/knuth_tho40/s*/est.out`,
  settings b; under c 404 nodes]; depth 3: 13,459 nodes [e] x 2,292 of 4,737 sampled open nodes branching
  [cognac `infra/phase1/runs/tho40-est2-261119/`], which over-counts (finished nodes were not sampled). Seconds
  per iteration: 0.161 and 0.135 at depths 2 and 3 on the Xeon box, x 1.29 for c7a
  [`../../bound_lab/tho40_dives/results/analyze.txt`, README]. The run: 5,728 CPU h (5,182-6,255) [e, the
  readout of tho40-est2-261119].
- wil50: 15 and 304 (of 635) [m `../../wil50_feasibility/sampling/results/stageA_summary.md`, cut.zone 0.005;
  fewer branch under 0.0075]; depth 3: 14,192 x 103 of 400 sampled, no cut phase, so an upper bound. 0.37 and
  0.32 s an iteration at depths 2 and 3 [m, the same and stage B]. The run: `../../estimates/wil50/refit.txt`.
- sko49: 10 and 313 (of 333) [m `../../estimates/sko49/tune/sampling_C/stageA_summary.md`, arm C]; depth 3:
  13,682 x 198 of 250. 0.383, 0.347 and 0.307 s an iteration at depths 2, 3 and 4. The run:
  `../../estimates/README.md`, Cw.
- tai50b, tai60b: nodes and fathomed shares per depth, totals [m cognac `infra/phase1/runs/tai50b-r9/`,
  `tai60b-261117/results.txt`]. Seconds per iteration [e]: the root's time over its iteration caps (4,297 s /
  6,100 = 0.70; 5,344 s / 3,100 = 1.72), scaled to the children's dimension by the power 5.2 measured between
  our lookahead logs (0.012, 0.029 and 0.49 s at dimensions 24, 28 and 49 on the cloud box).
- The per-instance figures were collected by four subagents and the ones the table rests on spot-checked here
  against the files named (the two stage A summaries, the tho40 readout, the two results files).
- Mixed machines: each row is in its own estimate's unit, so shares are comparable and core-hours across
  instances are not.
