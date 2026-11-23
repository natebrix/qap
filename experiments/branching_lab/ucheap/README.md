# ucheap: is there a cheap branching rule with the size score? (branching_lab)

Nov 22 2026, Nathan's question after the port of DNN strong branching: the root's row is now chosen from the
children's DNN bounds with a size score; can a transform of QPB's dual matrix U do the same for free, perhaps
with rule 2's fathom count first, and should rule 9's two lines be opened up? [m] = measured here.

**U itself: no. The children's own QPB bounds (rule 4's lookahead, which rule 9 already runs on its two lines):
yes, if they are scored with the size score and not summed, and over every line and not the two by marginal
score. Run as a rule (Nov 23, `trees/`): nug30 closes in 3,420 nodes and 2,392 CPU s against production's 4,285
to 4,365 and 3,314 to 3,421; tho30 in 8,581 and 6,350 against 14,815 and 9,375.**

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
- tai*b: not looked at here. On Nov 18 every cheap signal ranked tai60b's four lines backwards (QPB is 29%
  under the cutoff there).

## Files

`ucheap.py` (U; `out.md`), `rootla.py` (the children's QPB bounds at the roots; `rootla_out.md`), `truthla.py`
(at the truth nodes; `truthla_out.md`), `siblings.py` (within sibling groups, from whole-tree traces).
`dumps/`: `*_root.err` U at the roots, `*_rootla.err` and `tho30_*`, `tho40_*` with `ls-LAcell row col bound`,
the QPB bound rule 4's lookahead gives each child. The print is three additions to the experiment build
(`../root`'s, in `latest/branch.C`): the line after each `U->me[di][dj] = std::max(...)` in
`select_row_lookahead` and `select_ind_lookahead`, and a rows-only call of `select_row_lookahead` in the
`QAP_LINE_SCORES` block when symmetry is left (`lacell.txt`).
