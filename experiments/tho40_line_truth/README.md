# tho40: does a line's marginal score predict its subtree? (Nov 9 2026)

The tho30 study (`experiments/tho30_solve/line_truth/`) repeated on tho40,
with one change forced by the instance: tho40's depth-1 subtrees run to
hundreds of thousands of nodes, so the truth measurement sits at depth 3,
where the Nov 8 traced run (`experiments/tho40_estimate/traces/tho40.trace`)
finished subtrees of a few hundred nodes.

## Setup

Instance `qaplib/tho40t.dat` (the transposed copy the symmetry file applies
to), upper bound 240,516, the Nov 8 run's own recipe
(`infra/phase1/recipes/nug30.par`: the record settings with rule 9, DNN at
depths 1 to 7), so rule 9's choice and the subtree truth come from the same
par. Three depth-3 nodes from the trace with finished, unfathomed subtrees of
moderate size and a root fix outside the grid's middle row (no leftover
symmetry):

| node | first-solve gap | subtree in the Nov 8 run |
|---|---|---|
| 26:3, 31:7, 37:0 | 0.94% | 182 nodes |
| 25:32, 26:8, 39:0 | 1.03% | 182 |
| 24:7, 26:1, 39:0 | 1.18% | 359 |

Eight candidate lines per node (rule 9's, the top line of the other side,
then by marginal score), 37 children per line; every child gets (a) a cold
first DNN solve with give-up 300 (its shortfall, cutoff minus lb where
positive) and (b) its subtree to completion with `-F` (the node's three fixes
plus the child's), nodes and CPU summed plus 1. 24 lines, 888 children, about
2 h wall on 8 jobs. Same driver, with its constants made options:

```
cd experiments/tho30_solve/line_truth
python3 line_truth.py --dat ../../../qaplib/tho40t.dat --par ../../../infra/phase1/recipes/nug30.par \
  --ub 240516 --n 40 --rank-x "" --node 26:3,31:7,37:0 --node 25:32,26:8,39:0 --node 24:7,26:1,39:0 \
  --lines 8 --jobs 8 --out ../../tho40_line_truth/results
```

## Result

Full tables in `results/summary.md`. Ranks are 1 = best.

| node | rule 9 line | rank by (b) nodes | top marginal score (rule 8's line) | its rank | best line | best / rule 9, nodes | best / rule 9, CPU |
|---|---|---|---|---|---|---|---|
| 26:3,31:7,37:0 | col39 | 3 of 8 | row2 | 1 | row2 | 0.60x | 0.97x |
| 25:32,26:8,39:0 | col7 | 4 of 8 | row2 | 2 | row35 | 0.60x | 1.08x |
| 24:7,26:1,39:0 | col39 | 2 of 8 | row35 | 3 | col32 | 0.70x | 0.86x |

Spearman rank correlation with the truth (b):

| node | (a) shortfall vs nodes | (a) vs CPU | marginal score vs nodes | score vs CPU |
|---|---|---|---|---|
| 26:3,31:7,37:0 | 0.92 | 0.90 | 0.60 | 0.62 |
| 25:32,26:8,39:0 | 0.77 | 0.07 | -0.23 | -0.19 |
| 24:7,26:1,39:0 | 0.76 | 0.76 | 0.55 | 0.55 |

1. **The cold first-solve shortfall (a) predicts the subtree here too:**
   0.76 to 0.92 against nodes at all three nodes, as on tho30 (0.71 to 0.93).
   (a)'s top line is the truth's best at every node.
2. **The marginal score is a weaker signal than (a) but not useless on
   tho40:** 0.55 to 0.60 at two nodes, -0.23 at the third, against -0.2 to
   -0.55 on tho30. The top-scored line is 1st, 2nd and 3rd by truth. The
   difference from tho30 is plausibly that at depth 3 with gaps near 1% most
   children fathom cold on every line (29 to 35 of 37), so the lines differ
   less.
3. **Rule 9 lands 2nd to 4th of 8; the best line has 0.60x to 0.70x its
   nodes.** Rule 9's lookahead between the top row and the top column picked
   the column every time; the best line was a row at two nodes and another
   column at the third. The ceiling from a perfect choice at these nodes is
   30 to 40% of the nodes.
4. **CPU at depth 3 is not a fair measure in this design.** Each child is a
   separate `-F` run that pays its own first solve (about 15 to 25 CPU s at
   m = 37 when it does not fathom quickly), so the 37 children cost 400 to
   600 CPU s per line before any subtree is searched, and the subtrees here
   are 110 to 1,175 nodes. That fixed cost compresses the CPU ratios (0.86x to
   1.08x) and gives the second node's 0.07 CPU correlation. In a real run a
   child's first solve is warm-started from its parent's state and costs a
   fraction of that; the node counts carry the comparison.

## Against tho30

| | tho30, depth 1 | tho40, depth 3 |
|---|---|---|
| (a) vs (b) nodes | 0.71, 0.71, 0.93 | 0.92, 0.77, 0.76 |
| marginal score vs (b) nodes | -0.21, -0.55, -0.22 | 0.60, -0.23, 0.55 |
| rule 9's rank of 8 | 5, 5, 1 | 3, 4, 2 |
| best / rule 9, nodes | 0.68, 0.81, 1.00 | 0.60, 0.60, 0.70 |

The picture is consistent: the children's first-solve bounds are the
predictive signal, the marginal score alone is not reliable, and rule 9's
two-line lookahead leaves 20 to 40% of the nodes on the table at the nodes
tested. What this does not yet say is whether a lookahead over more lines
pays at scale: on tho40 the top of the tree is where the nodes are (the Nov 8
estimate is 0.6 to 3 million), a depth-1 or depth-2 node has 37 to 39
children at 20 to 25 s each cold, and only the few hundred nodes at depths 1
to 3 would use it. A cheap version would score only the children of the top
3 or 4 lines by marginal score (not 8), warm-started from the node's own
state, and be switched off below depth 3 or so.

## Caveats

- Three nodes, one depth, one instance; neighbouring lines are often tied
  (110 vs 110, 182 vs 182 nodes).
- Rule 9's own choice is among the candidates by construction; the best line
  of the 8 is a lower bound on the best of all 74 free lines.
- Subtree sizes under `-F` differ slightly from the Nov 8 run's (182 vs 182,
  182 vs 182, 359 vs 359 for rule 9's line here: identical at all three, as
  it happens).
