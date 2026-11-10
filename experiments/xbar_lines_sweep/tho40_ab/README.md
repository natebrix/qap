# tho40 A/B: `branch.xbar_lines` 1 vs 2 on three depth-3 subtrees (Nov 9 2026)

The cheap test the sweep's reading asked for: which of the two regimes
(tho30's, where a wider lookahead shortlist helps; nug30's, where it hurts)
tho40 is in. Each of the three depth-3 nodes of the tho40 line-truth study
(`experiments/tho40_line_truth/`) is solved to completion as a single
`bnbdrv -F` run under the Nov 8 recipe (`infra/phase1/recipes/nug30.par`,
upper bound 240,516, `qaplib/tho40t.dat` with its symmetry file), once with
the default 1 line per side and once with 2 (`logs/x2.par` is the recipe
plus `branch.xbar_lines = 2`). Six single-threaded runs at once on the
laptop, `logs/ab.sh`, 6.7 min wall in all.

| node | lines per side | nodes | fathomed | FW iters | CPU s |
|---|---|---|---|---|---|
| 26:3,31:7,37:0 | 1 | 182 | 177 | 17,126 | 236 |
| 26:3,31:7,37:0 | 2 | 110 | 107 | 14,278 | 197 |
| 25:32,26:8,39:0 | 1 | 182 | 177 | 17,126 | 267 |
| 25:32,26:8,39:0 | 2 | 182 | 177 | 24,878 | 271 |
| 24:7,26:1,39:0 | 1 | 359 | 349 | 34,701 | 381 |
| 24:7,26:1,39:0 | 2 | 289 | 281 | 40,628 | 340 |
| total | 1 | 723 | | | 884 |
| total | 2 | 581 | | | 808 |

The 1-line runs reproduce the Nov 8 run's subtree sizes exactly (182, 182,
359). The first two nodes give identical node and FW-iteration counts at
1 line: their subtrees appear to be images of each other under the
instance's symmetry (the same root bound, the same tree), so they are one
data point there, not two; at 2 lines they diverge because the wider
shortlist ranks the candidate lines differently at the two images.

## Reading

- **tho40 is in tho30's regime, at least at depth 3.** 2 lines per side
  gives 20% fewer nodes (723 to 581) and 9% less CPU (884 to 808 s) over the
  three subtrees; two of the three subtrees shrink (0.60x and 0.80x), the
  third is unchanged. No subtree got worse, which is the nug30 failure mode.
- The per-node CPU gain is smaller than the node gain because the second
  line's QPB lookahead costs FW iterations at m = 36 and the subtrees are
  short; in the full tree, where the depth-1 and depth-2 choices decide
  hundreds of thousands of nodes, the node saving is what matters.
- Three subtrees of a few hundred nodes under one incumbent is a small
  sample, and the 0.60x case is one node. It is enough to say 2 lines is
  the setting to carry into any long tho40 run, not enough to size the
  gain; the tho30 sweep (11,912 at 2 vs 14,256 at 1, -16%) is the better
  guide to the order of magnitude.
- Not tried here: 4 and 8 lines (tho30 gained another 6% at 8). The same
  six-run script with `branch.xbar_lines = 8` would say whether tho40
  follows tho30 there too.
