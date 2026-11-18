# tai60b: which line at the heaviest node? (Nov 17-18 2026, running: partial results)

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
   marginal score. 4 processes on this box's 4 cores, about 2.3 GB each to start, plus about 0.4 GB
   for each child that branches: its DNN state stays in memory for the grandchildren set aside.
3. Score each line by its children: how many fathom, the shortfall (cutoff - lb) of those that do
   not, and the CPU. Calibrate against production's subtrees under facility 1 once tai60b ends.

## Results (partial, 02:20 UTC Nov 18: r1, r45 and r8 done, c26 running)

`python3 experiments/tai60b_lines/score.py` (from the repository root) prints the table and the r1 check.
Each arm bounds its children in descending order of the free index (location 59 down for a row,
facility 59 down for a column), so the arms have covered different children. Marginal score at 43:58
(`QAP_LINE_RANK`, this state): c26 0.145, r45 0.142, r8 0.139, r1 0.100 (13th); c26 is also rule 9's choice.

| arm | children bounded (of 58) | first solve fathomed (D) | cut phase fathomed (C) | branched | branched children's shortfall |
|---|---|---|---|---|---|
| production r1 | 58 | 49 | 3 | 6 | path child 1:33 758k; 510k, 509k, 383k, 373k, 319k |
| r1 (done 00:54) | 58 | 46 | 6 | 6 | path child 1:33 1,000k; 772k, 710k, 693k, 632k, 530k |
| r45 (done 02:09) | 58 | 42 | 6 | 10 | path child 45:26 1,029k; 719k, 665k, 662k, 651k, 329k, 222k, 188k, 175k, 11k |
| r8 (done 01:55) | 58 | 44 | 2 | 12 | path child 8:29 1,133k; 999k, 982k, 853k, 757k, 756k, 726k, 697k, 691k, 490k, 416k, 111k |
| c26 | 25 | 12 | 0 | 13 | path child 45:26 1,029k; 879k, 810k, 809k, 727k, 703k, 661k, 623k, 611k, 536k, 479k, 457k, 418k |

A child that branches costs about 1,100 CPU s here (bounded, then set aside); one the first solve
fathoms, about 45; one the cut phase fathoms, 300 to 1,000. In production the path child's subtree
is most of the line's (1:33: 16,406 nodes; the other five 170 to 2,223).

**The arms are weaker than production, by the same amount for every line.** Their 43:58 is 1.76M
short against production's 1.11M: their 44:45 came cold from `-F` (4.1M short against production's
1.89M, warm from the root). Yet r1's children branch or fathom exactly as production's did, all 58,
and the same six branch, each 0.16M to 0.31M further from the cutoff than in production. Five
that production's first solve fathomed (a few thousand over the cutoff) need the cut phase here, at 2-3x
production's CPU. Every arm starts from one parent state, so the comparison between lines should hold,
but absolute subtree sizes will not transfer.

**r1, production's line, leads; the marginal score has the order backwards.** r1 finished at 00:54 UTC
with all 58 children as production's (`58 of 58 children branch or fathom as in production`): 6 branching
children, the fewest of the four lines: r8 finished at 01:55 with 12, r45 at 02:09 with 10, and c26 has 13
after 25.
Its path child (1.00M short) is also the least short. Against r8 and c26 that is dominance: each of their
six worst children is further from the cutoff than r1's child of the same rank, and they have six more
besides. r45 is closer:
its path child is 45:26, the same node as c26's (both lines contain it), 1.03M short against r1's 1.00M.
Its other branching children are mostly less short than r1's, but there are nine of them to r1's five.
c26, the bound lab estimate's line and rule 9's choice from this state, branches on 13 of its first 25
children; rule 9 chose between c26 and r45 (rule 8's two lines here) and took the worse. The four lines' order by children so far, r1, r45, r8, c26, is the reverse of
the marginal score's, as in `../tho30_solve/line_truth/`.

**Rough subtree sizes: r1 smallest, r45 about 1.2x, c26 at least 1.5x, r8 about 3.4x.** Production's
subtrees under r1 give a scale: its five branching children off the path, 319k to 510k short, had 170 to
2,223 nodes (log size on shortfall: an e-fold every 127k, from 5 noisy points), and the path child 1:33
16,406, five times the rest. Weighting each line's branching children by exp(shortfall/127k), relative to
r1's, and its path child the same way: r1 19,600 nodes (production's own), r45 22,900, r8 67,100, c26 30,000
from its first 25 children. The arms' children are 0.16M-0.31M shorter of the cutoff than production's, but
a common shift cancels in these ratios. The path child carries most of it: r45's 45:26 is 29k shorter of the
cutoff than r1's 1:33, a factor 1.26 on this slope, which outweighs r45's lighter other children (0.75x).
That rests on extending the slope to a different node on the optimum's path; r1 against r45 is within this
model's error, r8 and c26 are not.

**Restarts reproduce the arms exactly.** Container restarts at about 20:35 and 20:58 UTC Nov 17 killed
the first attempt (`runs/*/*.attempt1.txt`, started 17:37 and 18:03) and a second one before its first
checkpoint; the arms restarted from scratch at 21:01 with a checkpoint every 300 CPU s. Every child bounded by both attempts (r1 24, c26 18, r45 30, r8 21) has the same
outcome and the same bound values to the last digit printed, at 8% less CPU the second time.

**c26 ran out of memory and resumed trimmed.** The OOM killer stopped c26 at 00:42 UTC (`exit 137`),
with 4 arms at 14.9 GB in all: each branched child's DNN state, about 0.4 GB, stays in memory for its
grandchildren. `run.sh trim c26` kept the 37 children still to bound and the one DNN state they warm-start from
(the checkpoint went from 571 MB to 61 MB; the whole one is `runs/c26/setaside.1.state`, uncommitted), and
it resumed at 00:50 from its 00:24 checkpoint (`resumed: 37 open nodes`). The children are bounded as
before: `-d 3` never bounds the grandchildren. A child bounded after a resume warm-starts from the
checkpoint's copy of 43:58's state, rounded to fp32 (`ship.bits = 32`) as production's warm starts from the
store are. `supervise.sh` keeps the arms going: each minute it trims and resumes an arm that stopped part
way, when the memory leaves room for it.

**A third container restart at about 01:32 UTC Nov 18** killed c26, r45 and r8 within 5 minutes of
their last checkpoints (01:27-01:31). `supervise.sh` trimmed each checkpoint (c26 35 children kept,
r45 10, r8 9, each with the one DNN state they use) and resumed them at 01:32-01:33: 5.4 GB in all,
against 12.6 GB before.

A fourth container restart at about 03:05 stopped c26 5 minutes after its 03:00 checkpoint;
`supervise.sh` trimmed it (30 children kept) and resumed it at 03:06. At that point c26 had 16 branching
children in its first 28.

Projected finish of c26 from its CPU per child so far (658 s; 17 branching children in its first 29): about 08:45 UTC Nov 18.
