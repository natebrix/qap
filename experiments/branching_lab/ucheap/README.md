# ucheap: is there a cheap branching rule with the size score? (branching_lab)

Nov 22 2026, Nathan's question after the port of DNN strong branching: the root's row is now chosen from the
children's DNN bounds with a size score; can a transform of QPB's dual matrix U do the same for free, perhaps
with rule 2's fathom count first, and should rule 9's two lines be opened up? [m] = measured here.

**U itself: no. The children's own QPB bounds (rule 4's lookahead, which rule 9 already runs on its two lines):
yes, if they are scored with the size score and not summed, and over every line and not the two by marginal
score. Run as a rule (Nov 23, `trees/`): nug30 closes in 3,420 nodes and 2,392 CPU s against production's 4,285
to 4,365 and 3,314 to 3,421; tho30 in 8,581 and 6,350 against 14,815 and 9,375. On three depth-3 subtrees of
tho40 it gains nothing (1,862 nodes against 1,818; one subtree 38% larger). On eight subtrees of sko42 rule 10 takes
0.61 of rule 9's nodes and 0.78 of its CPU, all of it from the depth-2 nodes' lines; on tai35b neither rule
gains (Nov 23, below).**

## The trees (`trees/run.sh PROBLEM exp_all`; cognac 1d7e7403 + `trees/lasize.patch`, 4 slots, the nug30 recipe)

`branch.la_exp = 32`, `branch.la_all = 1`: wherever the recipe's lookahead rule chooses (rule 9; rule 4 at the
root, which has no DNN bound), every free line is a candidate and the line with the least sum over its children
of exp(32 (s - s_max) / G) is taken, s the child's QPB shortfall, G the node's own QPB gap.

| | nug30 nodes | CPU s, all of it | tho30 nodes | CPU s, all of it |
|---|---|---|---|---|
| production (rule 9, its two lines) | 4,285 and 4,365 | 3,314 and 3,421 | 14,815 | 9,375 |
| DNN strong branching at depths 0 to 2 (`../port`) | 2,980 | 8,617 | 9,918 | 16,746 |
| **this rule** | **3,420** | **2,392** | **8,581** | **6,350** |
| against production | 0.79 | 0.71 | 0.58 | 0.68 |

- **The roots** [m]: r23 at nug30 and r5 at tho30, the DNN score's rows, for under a second each.
- **The rule's cost is small** [m]: everything but the DNN bound (QPB, the lookahead over every line, branching)
  is 49 CPU s of nug30's 2,392 and 113 of tho30's 6,350.
- **Predictions** (`trees/predictions.txt`, committed before the runs): nug30 3,100 to 3,700, met. tho30 11,000
  to 13,500, missed on the good side; and at 5:1 it took c20, not c9 as predicted from the truth node's dump
  (another recipe and incumbent there).
- **With the keys off** the build's nug20 trace is cognac's, node for node [m].
- **On nug30 the score does the work, not the open candidates** [m `trees/runs/nug30_exp`, `nug30_all`]:

  | nug30 | root row | nodes | CPU s |
  |---|---|---|---|
  | production: rule 4's score, rule 9's two lines | r0 | 4,285 and 4,365 | 3,314 and 3,421 |
  | every line, rule 4's score (`all.par`) | r0 | 5,639 | 3,543 |
  | the exponential score, rule 9's two lines (`exp.par`) | r23 | **3,261** | **2,250** |
  | the exponential score, every line (`exp_all.par`) | r23 | 3,420 | 2,392 |

  Opening the lines under the old score makes nug30 30% worse (as 8 lines a side did on Nov 9); under the
  exponential score it costs 5%. At the root both `exp` arms see every row (no marginals there: rule 4's
  candidates, 30 in the recipe).
- **On tho30 it is the other way: the open candidates are needed** [m `trees/runs/tho30_exp`]: the exponential
  score on rule 9's two lines gives 14,303 nodes and 9,443 CPU s, production's tree (14,815); on every line
  8,581 and 6,350. Rule 9's two lines by marginal score do not hold the lines the score wants there (r19, r28,
  r2 at depth 1 against c9, c20 with every line). So every line is the setting that works on both: 5% dearer
  than the best at nug30, 40% better at tho30.
- Not tried: any instance outside these two, and tho40, where its root row (r26) is not the DNN score's (r2).

## tho40: three depth-3 subtrees (`tho40sub/`, Nov 23): no gain

Three depth-3 nodes under 26:10,8:0 from the traced run tho40ab-32, each solved to the end single-threaded from
the node, incumbent 240516; tho40's recipe (`experiments/estimates/tho40/tho40.par`: rule 9, 2 lines a side)
against the same with `branch.la_exp = 32`, `branch.la_all = 1`.

| subtree | production: depth-3 line, branching children at depth 4, nodes, CPU s | the score, every line (`exp_all.par`) | the score, rule 9's two lines a side (`exp.par`) |
|---|---|---|---|
| 3:32 | c7, 9, 781, 1,647 | c39, 9, 607, 1,514 | c39, 677, 1,506 |
| 30:32 | c39, 8, 571, 1,302 | r5, 14, 787, 1,768 | c39, 571, 1,285 |
| 36:32 | c7, 8, 466, 998 | r15, 10, 468, 1,102 | c39, **359**, 762 |
| all | 1,818 nodes, 3,946 s | 1,862 nodes (1.02), 4,384 s (1.11) | **1,607 nodes (0.88), 3,553 s (0.90)** |

- **On rule 9's own lines the score gains 12% and loses nowhere** [m]: the control, run after the first two arms.
  So at tho40's depth 3, as at nug30, the lines to choose among are the DNN's marginal lines and the score to
  choose by is the exponential; at tho30 the marginal shortlist (one line a side there) holds nothing and every
  line is needed. Both are in cognac: rule 10 (rule 4's candidates) and rule 11 (rule 9's), same score.

- **The n = 30 result does not carry to these subtrees** [m]: one 22% smaller, one 38% larger, one the same. The
  prediction written before the runs (0.6 to 0.9 of production, no subtree more than 10% larger) failed.
- **Where it loses, it loses at the node's own choice** [m]: at 30:32 its line r5 leaves 14 branching children at
  depth 4 where production's c39 leaves 8 (504 nodes at depth 5 against 288; the same 245 at depth 6).
- These subtrees are three levels deep: nearly every child is fathomed by its DNN bound, and what counts is how
  many are not. The marginal lines of rule 9 come from the DNN solution; the QPB shortfalls do not see which
  children the DNN bound fathoms. nug30 said the same more mildly (rule 9's lines 5% better than every line).
- Production's sizes are the traced run's exactly (781, 571, 466).

## Nathan's conjecture (Nov 23): DNN strong branching near the root, this rule farther down

`trees/sb1_exp.par`, `sb1_exp_all.par`: `sb.depths = 1:1` (the root's row is the same under both rules at nug30
and tho30, so the root is left to the cheap rule), the exponential score everywhere else.

| | nug30 nodes | CPU s, all of it | tho30 nodes | CPU s, all of it |
|---|---|---|---|---|
| production | 4,285 and 4,365 | 3,314 and 3,421 | 14,815 | 9,375 |
| DNN strong branching at depths 0 to 2, rule 9 below (`../port`) | 2,980 | 8,617 | 9,918 | 16,746 |
| the cheap rule everywhere (nug30: rule 9's lines; tho30: every line) | 3,261 | 2,250 | 8,581 | 6,350 |
| **DNN strong branching at depth 1, the cheap rule below** | **3,140** | 4,794 (2,782 strong branching) | **7,757** | 9,826 (4,169 strong branching) |

- **Both halves hold on both instances** [m]. Near the root: DNN strong branching at depth 1 in place of the
  cheap rule, the cheap rule below in both: 3,261 to 3,140 (-4%) and 8,581 to 7,757 (-10%). Farther down: the
  cheap rule in place of rule 9 below the DNN's depth-1 lines: nug30 3,241 (`../root/runs/nug30_r23_d1_sum`) to
  3,140 (-3%); tho30, node by node under the same depth-1 line, smaller under all nine nodes where the line is
  the same (2,451 to 2,033, 1,017 to 684, 577 to 349 trace lines, ...).
- **At depth 1 the DNN's line is not always the better one** [m, tho30, the cheap rule below both]: better under
  5:11, 5:12 and 5:13 (2,033 against 2,736, 869 against 1,170, 349 against 560), worse under 5:10 (684 against
  339 under the cheap rule's r2).
- **At n = 30 the depth-1 child solves cost more than they save**: the totals above. The case is the open
  instances, where depth 1 is about 1% of a run (`../budget`).
- Predictions (`trees/predictions.txt`): nug30 2,950 to 3,250, met; tho30 8,000 to 9,300 with no clear gain
  expected, wrong: 7,757.

## sko42 and tai35b (Nov 23): where rules 10 and 11 would go next

For Nathan's retries of sko42 and tai50b. Runs on the cloud container, cognac 626230c8 (rules 10 and 11).
Predictions, written before each batch, and the results against them: `sko42sub/predictions.txt`,
`tai35b/predictions.txt`.

**sko42** (`sko42sub/`): eight subtrees of the solved run sko42-261120b, each solved to the end from its node
(single-threaded, a cold first solve, incumbent 15814) under `sko42.par` (rule 9, one line a side) and the same
with rule 11 (one and two lines a side) or rule 10 on every line. Nodes / CPU s.

| subtree (its depth) | rule 9 | rule 11 | rule 11, two lines | rule 10, every line |
|---|---|---|---|---|
| 36:1,7:0,5:7 (3) | 191 / 1,600 | 191 / 1,619 | 191 / 1,595 | 191 / 1,526 |
| 36:16,5:0,22:35 (3) | 338 / 2,524 | 338 / 2,359 | 264 / 2,140 | 192 / 2,064 |
| 36:3,22:0,24:41 (3) | 226 / 1,637 | 226 / 1,627 | 227 / 1,599 | 265 / 1,861 |
| 36:16,23:0 (2) | 197 / 2,064 | | 197 / 1,793 | 41 / 1,782 |
| 36:8,35:0 (2) | 347 / 2,723 | | 273 / 2,696 | 158 / 1,475 |
| 36:9,12:0 (2) | 310 / 2,471 | | 80 / 1,100 | 80 / 1,077 |
| 36:14,7:0 (2, large) | 575 / 5,014 | | | 578 / 5,754 |
| 36:15,27:0 (2, large) | 691 / 5,134 | | | 235 / 2,477 |
| the first six | 1,609 / 13,020 | | 1,232 / 10,923 (0.77 / 0.84) | 927 / 9,785 (0.58 / 0.75) |
| all eight | 2,875 / 23,168 | | | **1,740 / 18,016 (0.61 / 0.78)** |

- **The gain is the depth-2 node's line** [m]. At four of the five depth-2 nodes rule 10 takes another line than
  rule 9 (r1, c35, c6, c35 for r22, c7, c41, c41) and the subtree is 41 nodes for 197, 158 for 347, 80 for 310,
  235 for 691. At the fifth both take c35 and the subtree is the same size, with 15% more CPU under rule 10.
  The line decides how many of the node's 40 children branch at all.
- **From depth 3 down rule 10 is level with rule 9, not ahead** [m]: CPU 0.95, 0.82 and 1.14 on the three
  depth-3 starts and 1.15 below the shared depth-2 line of 36:14,7:0, where its lines leave 37 children to the
  cut phase for rule 9's 23. Rule 11 with two lines a side was never behind there (1.00, 0.85, 0.98).
- **Fewer nodes is not as much less CPU**: a child fathomed by its first DNN solve costs about 4 s here, one that
  needs its cut phase about 50. At 36:16,23:0 rule 10's one branching leaves 27 of 40 children to the cut phase.
- **What the sample is**: five of the run's 269 branching depth-2 nodes and three depth-3 nodes, subtrees of
  191 to 691 nodes; the 27 largest depth-2 subtrees (over 1,150 nodes) hold 61% of the run and none is here.
  Each arm starts cold at its node, so rule 9's sizes are not the traced run's (575 for 923, 691 for 959).
  Rule 10 at depth 1 (12 nodes, 14 to 103 CPU h each in the run) is not measured.
- Suggested files, `suggest/`: `sko42_r10.par` (rule 9 at the root and depth 1 as in the run, so the same 438
  depth-2 nodes; rule 10 on every line below) and `sko42_r10_sb.par` (the same with DNN strong branching at the
  root and depth 1, as `tho40c_sb.par`; unmeasured on sko42).

**tai35b whole** (`tai35b/`, the stand-in for tai50b; line_signals' recipe, 4 slots): no gain.

| | trace nodes | CPU s |
|---|---|---|
| rule 9 (`../lsla/t2/tai35b_prod`) | 1,017 | 6,935 |
| rule 11 | 1,079 | 6,951 |
| rule 10, every line | 1,123 | 9,211 |
| rule 9 with DNN strong branching at the root and depth 1 | 975 | 12,519 |

- QPB is 30% under the cutoff on a tai b instance, so the children's QPB shortfalls rank nothing: rule 10 takes
  c5 at the root for rule 9's r19 and twice as many depth-1 nodes branch. Rule 11 keeps rule 9's lines at the
  top (one line a side leaves the score two candidates).
- DNN strong branching: c12 at the root, a tie with r19 in the third digit (1.789 against 1.819); four firings,
  5,000 CPU s, for 4% of the nodes.
- tai50b's recipe in cognac (`infra/phase1/recipes/tai50b.par`: rule 9 at every depth, 6,629 nodes and 116.5
  core-hours in tai50b-r9) stays the best known.

## What was measured

| what the score reads | where | result |
|---|---|---|
| U (free): sum of (s/s_max)^P over a line's children, s = cutoff - bound - U_ij | 6 line-truth nodes, 8 lines each | Spearman with the solved subtrees +0.14 to +0.18 for P = 1 to 8; its first line 1.9x to 2.3x the best [m `out.md`] |
| the same | 5 grid roots, every row | Spearman with the DNN root score +0.01 to +0.37; first row r6 at nug30 (4,854 nodes; the best tree is 3,623) [m `out.md`] |
| the children's QPB bounds, summed | 6 truth nodes | best line first at 2 of 6, 1.38x [m `truthla_out.md`] |
| the children's QPB bounds, exp or a high power | 6 truth nodes | **best line first at 5 of 6, 1.22x (worst 2.29x)**; rule 9 as it is: 1 of 6, 1.41x; the DNN children summed: 6 of 6 |
| the same, every row | nug30's root | **r23**, the DNN score's first row and the best tree run; Spearman with the six trees +0.94 [m `rootla_out.md`] |
| the same | tho30's root | r5, the DNN score's first row (production's) |
| the same | tho40, sko49, wil50 roots | r26 (DNN rank 16 of 40; the DNN's r2 is its 2nd), r41 (DNN rank 7 of 49), r21 (DNN rank 2 of 50): no trees to say which is right |

- **Why U fails** [m]: every child is 10 to 18% under the cutoff by QPB, so U's shortfalls span 0.80 to 1 of their
  largest and no power changes their order (P = 1 to 8 pick the same row everywhere). U follows the plain sum of
  the DNN children's shortfalls (+0.37 to +0.57 by row), the score that picks the worst roots, and says nothing
  about a row's worst child (+0.01 to +0.19), which is where the tree is. Shifting the shortfalls (to the
  smallest, or onto the node's DNN bound) does not help.
- **Rule 2's count is 0** at every node here (no child's QPB bound + U reaches the cutoff; none of the 841 or
  1,369 lookahead bounds does either), so "fathoms first" is its tie-break at these depths. The exponential
  score is the smooth form of the count: a fathomed child adds 0, one far from the cutoff adds 1.
- **A child's own QPB bound ranks its siblings' subtrees almost as the DNN bound does** [m `siblings.py` on five
  nug30 trees and two tho30 trees]: Spearman within a sibling group +0.91 against +0.94 at the root, +0.69
  against +0.73 at depth 1 (nug30); it names the largest subtree less often (12 of 36 groups against 24 at
  depth 1).
- **The exponential needs no offset** [m]: over 45 depth-1 children of nug30, DNN shortfall = -430 + 0.82 x QPB
  shortfall (Pearson +0.86), and ln(subtree) has slope 0.042 on the QPB shortfall against 0.048 on the DNN's.
  A sum of exp(b s) orders lines the same whatever constant is added to s; a power does not.
- **Cost** [m]: every row's children at a root by rule 4's lookahead: 1 to 3 s (270 to 750 QPB solves), against
  0.5 to 15 core-hours for the DNN's. All 58 lines of a tho30 depth-1 node: 841 solves.
- **The lines must be opened** [m]: at the truth nodes the score's first line over all free lines is 6th to 63rd
  by marginal score, outside rule 9's shortlist. At tho30's three nodes it is col9, the study's best or tied
  best; at tho40's it is a row the study did not solve (row15, row27).
- Earlier, with rule 9's own arbiter (`../../xbar_lines_sweep`, Nov 9): 8 lines a side took tho30 from 14,256
  to 11,205 nodes and nug30 from 3,716 to 4,602.

## What it does not show

- Two trees have been run under the rule, both at n = 30. The truth nodes are six (tho30 depth 1, tho40 depth
  3), their lines eight chosen by marginal score; the roots with trees are two.
- At tho40's root it disagrees with the DNN score (r26 against r2), and the run now going uses r2.
- tai*b: tai35b whole, above: no gain from rule 10, 11 or DNN strong branching. On Nov 18 every cheap signal
  ranked tai60b's four lines backwards (QPB is 29% under the cutoff there).

## Files

`ucheap.py` (U; `out.md`), `rootla.py` (the children's QPB bounds at the roots; `rootla_out.md`), `truthla.py`
(at the truth nodes; `truthla_out.md`), `siblings.py` (within sibling groups, from whole-tree traces).
`dumps/`: `*_root.err` U at the roots, `*_rootla.err` and `tho30_*`, `tho40_*` with `ls-LAcell row col bound`,
the QPB bound rule 4's lookahead gives each child. The print is three additions to the experiment build
(`../root`'s, in `latest/branch.C`): the line after each `U->me[di][dj] = std::max(...)` in
`select_row_lookahead` and `select_ind_lookahead`, and a rows-only call of `select_row_lookahead` in the
`QAP_LINE_SCORES` block when symmetry is left (`lacell.txt`).
