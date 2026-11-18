# tai60b: which line at the heaviest node? (Nov 17-18 2026)

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

## Results (all four arms done, 04:36 UTC Nov 18)

`python3 experiments/tai60b_lines/score.py` (from the repository root) prints the table and the r1 check.
Each arm bounds its children in descending order of the free index (location 59 down for a row,
facility 59 down for a column). Marginal score at 43:58
(`QAP_LINE_RANK`, this state): c26 0.145, r45 0.142, r8 0.139, r1 0.100 (13th); c26 is also rule 9's choice.

| arm | first solve fathomed (D) | cut phase fathomed (C) | branched | branched children's shortfall | children's CPU s | rough subtree |
|---|---|---|---|---|---|---|
| production r1 | 49 | 3 | 6 | path child 1:33 758k; 510k, 509k, 383k, 373k, 319k | 10,246 | 19,638 (merged by 16:03 Nov 17) |
| r1 (done 00:54) | 46 | 6 | 6 | path child 1:33 1,000k; 772k, 710k, 693k, 632k, 530k | 12,449 | 19,600 (1.0x) |
| r45 (done 02:09) | 42 | 6 | 10 | path child 45:26 1,029k; 719k, 665k, 662k, 651k, 329k, 222k, 188k, 175k, 11k | 16,615 | 22,900 (1.2x) |
| c26 (done 04:36) | 38 | 0 | 20 | path child 45:26 1,029k; 923k, 879k, 832k, 810k, 809k, 756k, 727k, 703k, 661k, 657k, 656k, 650k, 623k, 623k, 611k, 536k, 479k, 457k, 418k | 23,400 | 38,800 (2.0x) |
| r8 (done 01:55) | 44 | 2 | 12 | path child 8:29 1,133k; 999k, 982k, 853k, 757k, 756k, 726k, 697k, 691k, 490k, 416k, 111k | 15,969 | 67,100 (3.4x) |

Each line has 58 children. The rough subtree is the model below, from production's subtrees under r1 as
merged by 16:03 UTC Nov 17 (7,369 nodes under 43:58 were still to bound then: the sizes are partial,
the ratios between lines are what the model gives).

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
children, the fewest of the four lines: r45 has 10, r8 12 and c26 20.
Its path child (1.00M short) is also the least short. Against r8 and c26 that is dominance: each of their
six worst children is further from the cutoff than r1's child of the same rank, and they have 6 and 14
more besides. r45 is closer:
its path child is 45:26, the same node as c26's (both lines contain it), 1.03M short against r1's 1.00M.
Its other branching children are mostly less short than r1's, but there are nine of them to r1's five.
c26, the bound lab estimate's line and rule 9's choice from this state, branches on 20 of its 58
children, more than any other line; rule 9 chose between c26 and r45 (rule 8's two lines here) and took
the worse. The four lines' order by branching children, r1, r45, r8, c26, is the reverse of the marginal
score's, as in `../tho30_solve/line_truth/`.

**Rough subtree sizes: r1 smallest, r45 about 1.2x, c26 about 2.0x, r8 about 3.4x.** Production's
subtrees under r1 give a scale: its five branching children off the path, 319k to 510k short, had 170 to
2,223 nodes (log size on shortfall: an e-fold every 127k, from 5 noisy points), and the path child 1:33
16,406, five times the rest. Weighting each line's branching children by exp(shortfall/127k), relative to
r1's, and its path child the same way: r1 19,600 nodes (production's own), r45 22,900, c26 38,800,
r8 67,100. c26 and r45 share their path child, 45:26 (20,500 on this model); c26's other 19 branching
children add 18,300 to r45's 2,400. r8 has fewer branching children than c26 but is larger: its path child
is 133k shorter of the cutoff than r1's (46,500) and its other children are further from it. The arms' children are 0.16M-0.31M shorter of the cutoff than production's, but
a common shift cancels in these ratios. The path child carries most of it: r45's 45:26 is 29k shorter of the
cutoff than r1's 1:33, a factor 1.26 on this slope, which outweighs r45's lighter other children (0.75x).
That rests on extending the slope to a different node on the optimum's path; r1 against r45 is within this
model's error, c26 and r8 are not.

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

c26 finished at 04:36 UTC, four hours ahead of a projection made from its first 29 children: its
branching children are consecutive (every free facility from 47 down to 26, 20 in a row), and the last 26
fathomed in the first solve at about 48 CPU s each.

## Verdict

- **Production's line, r1, is the best of the four** by every measure here: the fewest branching
  children (6), the path child closest to the cutoff, the least CPU, the smallest rough subtree.
- **Rule 9's pick depends on the parent's state, and the wrong pick costs about 2x.** Production's rule 9
  chose r1 from production's state (43:58 1.11M short). From the arms' weaker state (1.76M short, the same
  node and settings) rule 9 compared c26 and r45, rule 8's two lines, by their children's QPB, and chose
  c26: 20 branching children to r1's 6, about 2x r1's subtree; r45 would have cost about 1.2x.
- **The marginal score ranks these lines in reverse.** By score c26, r45, r8, then r1 13th; by branching
  children r1, r45, r8, c26; by rough size r1, r45, c26, r8. The children's first-solve shortfall, which
  ranks lines by their subtrees in `../tho30_solve/line_truth/`, picks r1 here as well.
- **Choosing by the children's bounds is affordable at a node this heavy.** Bounding all 58 children of
  one line cost 12,000-23,000 CPU s here; the four lines, 68,000. Production's subtree under 43:58 had used
  1.67M CPU s by 16:03 Nov 17 with 7,369 nodes still to bound, so the four-way choice is at most 4% of it,
  against a 2x loss from the wrong line. Most nodes are far lighter, so this pays only where a subtree is
  known to be large (an estimate's heavy nodes, the top of the tree).

Limits: one node; the sizes rest on a slope fit to five subtrees; the arms' parent is weaker than
production's (the order between lines should hold, the absolute sizes will not). Step 3 of the design,
calibrating against production's subtrees once tai60b ends, refits the slope on complete subtrees; it
covers r1 only, since production branched no other line here.
