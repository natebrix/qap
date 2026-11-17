# tai60b: which line at the heaviest node? (Nov 17 2026, running: partial results)

Question (Nathan, after tai60b's estimate came in at a third of the tree): would stronger branching
pay? The evidence says where: line choice matters most where the subtrees are largest, at the top
(wil50, `../wil50_feasibility/tuning/lines/`: the best of the top 4 lines 5.3x better than the worst
at the root, 2.4x at depth 1; about 5% per decision at depths 3-5). And the children's first-solve
shortfall ranks lines by their true subtrees where the marginal score does not
(`../tho30_solve/line_truth/`: Spearman 0.71-0.93 against nodes; the score's is negative).

## The node

tai60b-261117's live tree at 16:03 UTC Nov 17: 99% of the 42,305 merged nodes are under the depth-2
node on the optimum's path, 44:45 then 43:58 (19,691 merged and 7,369 still to bound), and its
near-mirror 44:58 then 43:45. Production branched 43:58 on facility 1 (row 1). The bound lab's estimate
(`experiments/bound_lab/tai60b` on `claude/cool-babbage-ilnogj`) had branched the same node on the
column through 45:26 (location 26) and measured its subtree at 14,270 nodes: the estimate sampled a
different tree from depth 3 down. So: how different are this node's lines, by the measure that
predicts subtrees?

## Design

`run.sh` with the experiment build (`force_line.patch` on cognac 8df9573, in a worktree; production's
`latest/` does not have it): `QAP_FORCE_LINE=D:rI|D:cJ` branches every node of depth D on facility I
or location J; `QAP_LINE_RANK=D` prints a depth's free lines ranked by marginal score. Without them
the build is production's, node for node (tai25b to depth 2: the same iterations, nodes and fathoms).

1. `run.sh parent`: 44:45 cold from `-F`, production settings (`infra/phase1/recipes/tai60b.par`,
   incumbent 608,215,055), branched on facility 43 as production did; its children set aside.
2. `run.sh arm NAME LINE`: 43:58 warm from 44:45's state, as in production, branched on LINE; every
   child bounded warm from it with the production settings (first solve, then the phase); the
   grandchildren of those that do not fathom set aside (`-d 3`). Arms: rule 9's own choice, the
   estimate's location 26, production's facility 1 if rule 9 does not pick it, and the next line by
   marginal score. 4 processes on this box's 4 cores, about 2.3 GB each.
3. Score each line by its children: how many fathom, the shortfall (cutoff - lb) of those that do
   not, and the CPU. Calibrate against production's subtrees under facility 1 once tai60b ends.

## Results (partial, 20:03 UTC Nov 17)

`python3 experiments/tai60b_lines/score.py` (from the repository root) prints the table and the r1 check.
Each arm bounds its children in descending order of the free index (location 59 down for a row,
facility 59 down for a column), so at 20:03 the arms had covered different children:

| arm | children bounded | first solve fathomed (D) | cut phase fathomed (C) | branched | branched children's shortfall |
|---|---|---|---|---|---|
| production r1 (all 58) | 58 | 49 | 3 | 6 | path child 1:33 758k; 510k, 509k, 383k, 373k, 319k |
| r1 | 24 | 17 | 5 | 2 | 772k, 632k (path child 1:33 not reached) |
| r45 | 27 | 23 | 1 | 3 | 719k, 665k, 175k (path child 45:26 not reached) |
| r8 | 18 | 14 | 1 | 3 | 999k, 756k, 726k (path child 8:29 not reached) |
| c26 | 17 | 12 | 0 | 5 | path child 45:26 1,029k; 810k, 623k, 611k, 457k |

A child that branches costs about 1,240 CPU s here (bounded, then set aside); one the first solve
fathoms, about 45 s.

**The arms are weaker than production, by the same amount for every line.** Their 43:58 is 1.76M
short against production's 1.11M: their 44:45 came cold from `-F` (4.1M short against production's
1.89M, warm from the root). Yet r1's first 24 children branch or fathom exactly as production's did
(24 of 24). The two that branch, 1:59 and 1:39, are 0.31M and 0.26M further from the cutoff than
production's. Two that production's first solve fathomed (1:57, 1:49, a few thousand over the cutoff)
need the cut phase here. The phase costs 2-3x production's CPU. Every arm starts from one parent
state, so the comparison between lines should hold, but absolute subtree sizes will not transfer.

**So far c26 is the worst.** That is the bound lab estimate's line, and rule 9's choice from this weaker state. It has five
branching children in its first 17, facilities 47 to 41 (except 43 and 44, fixed) one after another,
and its path child is 1.03M short. Production's path child under r1 was 0.76M short (from a parent
0.65M stronger). Too early for a ranking: each arm has 30-40 children to go. Projected finishes: r45
about 21:45 UTC, r1 22:20, r8 23:20, c26 00:35 Nov 18.
