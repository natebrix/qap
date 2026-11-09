# Does a line's marginal score predict its subtree? (tho30, Nov 9 2026)

Question: at a node, rule 8 picks the line (row or column) with the largest
sum of squared DNN marginals, and rule 9 looks ahead between the top row and
the top column by their children's QPB bounds. Does either ranking agree with
the truth, the size of the subtree each line would actually produce?

## Setup

Instance tho30, record recipe `../cutd2_cp.par` (cut phase through depth 2,
upper bound 149936). Three depth-1 nodes with the lowest bounds and no
leftover symmetry (facility 5 at locations 2, 1, 3; the middle-row locations
10 to 14 keep a reflection and were skipped).

Candidate lines at each node: rule 9's line, the top line of the other side,
and the next lines by marginal score, 8 in all (`QAP_XBAR_DEBUG=2` prints the
ranking; `cutd2_cp_r9.par` is the recipe with rule 9 for its choice). Every
line gets the same two treatments on each of its 29 children:

- (a) score: the child's first DNN solve from a cold start (`dnntest`, give-up
  300 iterations, no cut phase). The line's score is the sum over its children
  of cutoff minus lb where positive, i.e. how far the children that do not
  fathom are from fathoming.
- (b) truth: the child's subtree solved to completion with `-F` (the node's fix
  plus the child's), nodes and CPU seconds summed over the children plus 1.

Driver `line_truth.py`; per-child data in `results/<node>.json`, the tables in
`results/summary.md`. 24 lines, 696 children, 696 cold solves and 696 subtrees,
about 2.5 h wall on 8 jobs.

## Result

Ranks are 1 = best. "score" is the marginal-score rank, what rule 8 orders by.

| node | rule 9 line | rank by (b) nodes | rule 8 line (top score) | rank by (b) nodes | best line | best / rule 9, nodes | best / rule 9, CPU |
|---|---|---|---|---|---|---|---|
| 5:2 | col0 | 5 of 8 | row28 | 4 of 8 | col9 | 0.68x | 0.90x |
| 5:1 | col0 | 5 of 8 | col0 | 5 of 8 | col9 | 0.81x | 0.91x |
| 5:3 | col29 | 1 of 8 | row28 | 5 of 8 | col29 | 1.00x | 1.00x |

Spearman rank correlation with the truth (b):

| node | (a) shortfall vs nodes | (a) vs CPU | marginal score vs nodes | score vs CPU |
|---|---|---|---|---|
| 5:2 | 0.71 | 0.88 | -0.21 | -0.29 |
| 5:1 | 0.71 | 0.88 | -0.55 | -0.14 |
| 5:3 | 0.93 | 0.95 | -0.22 | -0.26 |

Three findings.

1. **The marginal score does not predict the subtree.** Its correlation with
   the truth is negative at all three nodes. row25 ranks 2nd or 3rd by score at
   every node and is the worst or second-worst line by truth (3,091 / 4,789 /
   1,812 nodes against 864 / 1,850 / 385 for the best). row28, the top score at
   two nodes, is mid-table. Within these 8 candidates, ordering by score is no
   better than random.
2. **The cold first-solve shortfall (a) does predict it.** Spearman 0.71 to
   0.93 against nodes, 0.88 to 0.95 against CPU. The line whose children are
   closest to fathoming at their first solve is the best or second-best line at
   every node, and (a)'s worst line is the truth's worst or second-worst. (a)
   costs 29 cold DNN solves per line, about 3 s each, so about 90 s per line
   at depth 1; the subtrees it ranks cost 225 to 1,421 CPU s.
3. **Columns beat rows here, and rule 9's lookahead gets that right.**
   The five column lines tested are the five best or near-best by truth at
   every node; the row lines are all worse. Rule 9 chose a column at every node
   (the QPB children bounds favour it), landing 1st, 5th and 5th. Rule 8 would
   have chosen row28, row28 and col0: 4th, 5th and 5th. Where rule 9 was 5th,
   the best column (col9 or col29) was 0.68x to 0.81x its nodes and 0.90x its
   CPU. The shortfall (a) is far from perfect: at 5:3 it ranks rule 9's col29
   3rd when it is actually 1st by nodes, with col9 tied at 385 nodes.

So the headline: ranking by (a) matches the ranking by (b) well (0.7 to 0.9),
ranking by the marginal score does not. Rule 9's line is 1st at one node and
5th at two; the best line beats it by 19 to 32 % in nodes and about 10 % in
CPU where it is not already best. The ceiling available from a perfect
depth-1 choice is modest in CPU (10 %) but the marginal score alone would be a
poor way to reach it; the children's first-solve bounds carry the signal.

## Caveats

- Three nodes, 8 lines each, one instance; the margins between neighbouring
  lines (e.g. 864 vs 894 nodes) are within what a different upper bound or
  cut-phase depth would move.
- (b) was run cold per child (`-F` from the root), so each child pays a cold
  root solve and the node's own DNN state is not warm-started; absolute node
  counts are therefore a little above what the full tree would see, but the
  treatment is identical across lines so the comparison is paired.
- (a) ignores the fathomed children's slack (negative shortfall clipped to 0);
  the number of children that fathom cold ("fathom cold" column in
  `results/summary.md`) tells the same story.

## Reproduce

```
cd experiments/tho30_solve/line_truth
export VECLIB_MAXIMUM_THREADS=1
python3 line_truth.py --node 5:2 --node 5:1 --node 5:3 --lines 8 --jobs 8
python3 line_truth.py --summary-only
```

Needs `latest/bnbdrv` and `latest/dnntest` built from a commit with the
`QAP_XBAR_DEBUG=2` ranking print (this one or later).
