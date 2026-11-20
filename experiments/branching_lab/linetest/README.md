# The depth-7 line test at tai60b's chain node 10:39: c47 against c54 (Nov 19 2026)

branching_lab agent "linetest". This tests the landscape signal's prediction
(`../landscape/README.md`, "Rule and the cheapest test"). At tai60b's depth-7 chain node
`44:45,43:58,1:33,6:24,7:19,0:32,10:39`, production branched on column 54 (rule 2), and that line
leaves 9 surviving children. Column 47 would leave 2.

## Answer

**c47 left fewer branching children than c54, 5 against 9, but not the 2 predicted. The kill rule is
met: c47 has 5 or more.** These counts are measured [m] from one parent state whose depth-7 node
reaches production's bound strength (406k short in production, 416k here; the "strong" arms `sc47`
and `sc54`). From the cold state with production's settings throughout (937k short; arms `c47`,
`c54`) the counts are 7 against 11. There c47 also has 8 children fathomed only by the cut phase,
against c54's 0. Step 4 (bounding the survivors' subtrees) was not run, as the kill rule says.

Two findings change how to read the kill rule:

- **The prediction of 2 failed in the local search, not in the survive rule's threshold.** The
  landscape signal used 30 ILS kicks per child. c47's children are hard to complete: with 3,000
  kicks, 29:47 drops from x = 1.74 to 0.69 (x = a child's excess over the node's slack) and 45:47
  from 1.39 to 0.98. In one seed 31:47 drops to 1.99; the best of three seeds puts it at 1.07, just
  over the threshold [m `effort.py`, `effort_all.py`]. With 3,000 kicks the rule predicts c47's
  survivors at 4 of 5, with no false positive, and c54's at 9 of 9 with none. c54's children are
  easy: their excess does not move between 30 and 3,000 kicks. So the signal ranks c47 above c54
  correctly, at an honest 4-5 against 9.
- **c47's survivors look lighter, though that is not measured.** On a log-linear fit of production's
  own subtrees under c54's children against their shortfalls (an e-fold every 88k of shortfall,
  correlation 0.87, 9 points), c47's four off-path survivors would hold about 8.7k nodes. c54's
  eight would hold 55k (production measured 60k) [e `subtrees.py`]. c47's path child 47:47 is
  330k short against c54's path child 14:54 at 383k, which held 54k nodes in production. That fit
  is a regression from one line's children applied to another line's, so it is a hypothesis for
  step 4, not a result.

So the depth-7 saving is not "2 survivors instead of 10" as predicted, and by the kill rule set in
advance this test stops here. The measured gap, 5 survivors against 9 that are also less short of
the cutoff, is still large enough that step 4 is the cheapest next test if Nathan wants it. That is
for him to decide; the kill rule says stop.

## Down the tree under c47

**Answer (Nov 20, 09:22 UTC; stopped at the 10-hour wall-time limit, not by the kill rule).** Depths 9 to 17 under
c47 are bounded. At every one of them c47's branching count stayed far below half of production's under c54: 7 to 28
against 56 to 199. Its largest share is 15% (depth 13, 28 against 185). Summed over depths 9-17: **5,017 nodes and 111
branching under c47, against production's 50,064 nodes and 1,169 branching, a tenth** [m]. CPU: 120,114 s on this box
against production's 1,674,126 s on the cloud boxes. Scaled by this box's measured speed (0.53x production's CPU per
node at depth 8), c47's levels cost about a seventh of production's [e]. The subtree did not close. Depth 18 has 469
open nodes (production bounded 5,203 there); the limit stopped its first 56 (50 D, 4 C, 2 B), and the rest are still
set aside. Do not read past depth 17. Production's branching rose again at depths 19-22 (146, 335, 504, 249), where its
mirror-like structure sits, and nothing here says whether c47's tree does the same. **The path to the best known
solution is intact** (`pathcheck.py`, `out/pathcheck.txt`): every node on it, depths 7-17, branched, with bounds from
607,798,661 to 608,129,286, all below the cutoff 608,215,054 (the incumbent was best known + 1). No memory limit was
near: under 5 GB of 15 with 4 processes, each holding one parent's state. Three container restarts (about 00:00, 03:25
and 04:10 UTC) cost only the chunks in progress, which were rerun from their parent's state. They count against the
10 hours. To go on: `python3 down.py 18` resumes from depth 18's set-aside states, which are in this container only and not
committed (state files). Of depth 18's 44 chunks, 5 had finished; a rerun keeps them.

How it ran: level by level under the strong c47 arm's survivors, with production's settings (`tai60b.par`) and production's own line
rule throughout (no forced lines): depth D is bounded from the states depth D-1 set aside, and depth D+1 is set aside
(`down.py`; a parent's children in pieces of at most 13, each with a copy of the parent's state, 4 processes (depth 10: one chunk per parent); `runs/down/dD/nodes/*.txt` the traces). Kill rules, fixed by the coordinator
before the run: stop if c47's branching count at a depth reaches half of production's there; stop at 10 hours of wall
time, or if memory will not hold a level. Production's rows: the tai60b-261117 trace under `...,10:39`
(`out/production_by_depth.txt`, checked against the coordinator's figures). Every row is measured [m]; CPU is this box's
against production's cloud boxes (this box ran the depth-8 children in 0.53x production's CPU).

<!-- down-table -->
| depth | nodes bounded | D / C / B (other) | CPU s [m, this box] | branching shortfall: sum / largest | production under c54: nodes / B / CPU s (cloud boxes) |
|---|---|---|---|---|---|
| 9 (from above) | 260 | 247 / 4 / 9 | 9,002 | 1,152k / 254k | 520 / 56 / 80,387 |
| 10 | 459 | 444 / 4 / 11 | 11,981 | 1,018k / 248k | 2,856 / 103 / 191,477 |
| 11 | 550 | 525 / 12 / 13 | 14,725 | 1,179k / 186k | 5,150 / 118 / 255,515 |
| 12 | 637 | 612 / 12 / 13 | 12,979 | 1,505k / 328k | 5,782 / 199 / 289,072 |
| 13 | 624 | 570 / 26 / 28 | 22,046 | 2,979k / 234k | 9,552 / 185 / 319,299 |
| 14 | 1,316 | 1,268 / 37 / 11 | 18,650 | 893k / 271k | 8,695 / 194 / 229,430 |
| 15 | 506 | 460 / 37 / 8 (1) | 14,057 | 373k / 112k | 8,924 / 93 / 148,859 |
| 16 | 359 | 340 / 11 / 7 (1) | 7,472 | 676k / 318k | 4,185 / 100 / 82,632 |
| 17 | 306 | 283 / 11 / 11 (1) | 9,202 | 1,123k / 223k | 4,400 / 121 / 77,455 |
<!-- /down-table -->

### One level deeper (depth 9; the coordinator's first follow-up, Nov 19 22:25-23:10 UTC)

**Under c47, depth 9 has 9 branching nodes against production's 56 under c54, so c47 is materially
better one level down. The kill rule, fixed in advance at 28 (half of 56), is not met, and per that
rule I stop here.** [m `d9score.py`, `out/d9score.txt`] The 260 depth-9 nodes under the strong c47
arm's five survivors were bounded from `runs/sc47/out.state`. Production's settings applied
(`tai60b.par`), with production's own line rule at depth 8 (no forced line), and depth 10 set aside
(`-d 9`; `run.sh deep`, `splitstate.py`, 4 processes, `runs/d9_*`). Production's figures are from
the tai60b-261117 trace.

| depth-8 node | line (production's rule) | depth 9: D | C | B | branching children's shortfalls | CPU s |
|---|---|---|---|---|---|---|
| 47:47 (path) | c10 | 48 | 1 | 3 | 254k, 238k, 179k | 2,692 |
| 46:47 | c53 | 48 | 1 | 3 | 214k, 175k, 18k | 2,502 |
| 29:47 | r9 | 50 | 0 | 2 | 54k, 14k | 1,626 |
| 45:47 | c59 | 51 | 0 | 1 | 6k | 1,093 |
| 31:47 | c59 | 50 | 2 | 0 | | 1,088 |
| **c47 total** | | 247 | 4 | **9** | max 254k, sum 1,152k | 9,002 |
| **production c54 total** (10 depth-8 nodes, all branched on c2) | | 442 | 22 | **56** | max 411k, sum 11,036k; 4-7 per node but 21:54 (0) | 80,387 (cloud boxes) |

The 56 have shortfalls of 8k to 411k, and every production depth-8 node except 21:54 has 4 to 7
branching children. c47's 9 have 6k to 254k. Its path child 47:47 has 3, against production's
14:54 with 7 at 117k-411k. CPU is not like for like: this box bounded the strong c54 arm's children
in 0.53x production's CPU (6,460 s against 12,153 s). Scaled by that, c47's depth 9 costs about a
fifth of production's [e]. The comparison rests on one assumption: that the strong state stands in
for production's below depth 8, as it does at depth 8 (c54's children within about 10k of
production's shortfall). The c54 arm was not taken to depth 9 here, so production's trace is the
reference. Counting branched nodes at depths 8 and 9 together: c47 5 + 9, c54 10 + 56 in production.

What that does not show: whether the gap holds deeper. The path child 47:47's subtree is the rest of
the chain, and its three depth-9 survivors at 179k-254k are lighter than 14:54's seven. The deeper
levels are Nathan's decision.

## Child by child, the strong state (10:39 416,393 short; production 405,921)

x = excess / production's slack (405,921). "x 30" is the landscape signal's search (`predict.py`,
p0 and mirror starts, 30 kicks); "x 3000" is `effort_all.py` (p0, mirror and ILS starts, 3,000
kicks, seed 0). The rule predicts that a child branches iff x <= 1. The outcome is the child's in
the strong arm [m `runs/sc47`, `runs/sc54`]: B branched, C fathomed by the cut phase, D fathomed by
the first solve. Shortfall = cutoff - the final bound (negative means fathomed, by that margin).
All other children are D with x 30 > 2.2.

**c47** (prediction: 2 survivors; realised 5 B, 3 C, 45 D)

| child | x 30 | x 3000 | outcome | shortfall |
|---|---|---|---|---|
| 47:47 (path) | 0.00 | 0.00 | B | 329,620 |
| 46:47 | 0.22 | 0.21 | B | 338,205 |
| 45:47 | 1.39 | 0.98 | B | 60,016 |
| 29:47 | 1.79 | 0.69 | B | 157,889 |
| 31:47 | 2.44 | 1.99 (1.07 best of 3 seeds) | B | 92,495 |
| 28:47 | 2.23 | 1.91 | C | -1,646 |
| 26:47 | 2.69 | 2.17 | C | -1,546 |
| 40:47 | 2.96 | 1.79 | C | -4,798 |

**c54** (prediction: 9 survivors; realised 9 B, 1 C, 43 D; production had 10 B, 43 D)

| child | x 30 | x 3000 | outcome | shortfall | production: outcome, shortfall |
|---|---|---|---|---|---|
| 14:54 (path) | 0.00 | 0.00 | B | 382,757 | B 382,605 |
| 18:54 | 0.16 | 0.16 | B | 351,017 | B 357,622 |
| 24:54 | 0.26 | 0.25 | B | 345,675 | B 350,561 |
| 19:54 | 0.41 | 0.38 | B | 373,809 | B 380,655 |
| 17:54 | 0.41 | 0.40 | B | 228,396 | B 225,571 |
| 23:54 | 0.60 | 0.58 | B | 446,081 | B 445,224 |
| 20:54 | 0.66 | 0.64 | B | 161,891 | B 157,271 |
| 16:54 | 0.66 | 0.64 | B | 177,877 | B 174,501 |
| 25:54 | 0.79 | 0.76 | B | 120,670 | B 128,052 |
| 21:54 | 1.13 | 1.11 | C | -345 | B 4,938 |
| 15:54 | 1.48 | 1.42 | D | | D |

The strong c54 arm is production's branching at 10:39 to within about 10k of shortfall per child,
so its state is a fair stand-in for production's. The one difference is 21:54, the soft-edge child
at x 1.1: production branched it with 4,938 to spare, while here the cut phase fathomed it with 345.

## Per arm [m]

| arm | 10:39 shortfall | D | C | B | first-solve iterations | cut iterations | children CPU s | 10:39 CPU s |
|---|---|---|---|---|---|---|---|---|
| production c54 (tai60b-261117) | 405,921 | 43 | 0 | 10 | | | 12,153 (cloud boxes) | 1,165 |
| sc54 (strong) | 416,393 | 43 | 1 | 9 | 2,705 | 9,910 | 6,460 | 4,855 |
| sc47 (strong) | 416,393 | 45 | 3 | 5 | 3,922 | 6,942 | 5,560 | 4,863 |
| c54 (cold, production settings) | 937,130 | 42 | 0 | 11 | 2,398 | 11,000 | 7,066 | 636 |
| c47 (cold, production settings) | 937,130 | 38 | 8 | 7 | 3,178 | 9,898 | 6,883 | 636 |

Each branched child here costs about 550-630 CPU s to bound: a 110-iteration first solve, then a
1,000-iteration cut phase that stops at production's give-up (`cut.giveup = 1000:0`). A C child
costs 140-510 s and a D child 10-50 s. Per ADMM iteration that is about 0.5 s at m = 52 [m]. RSS was
1.2-2.5 GB per process, growing about 0.15 GB per branched child set aside.
`out/score_*.tsv` has every child. The cold arms' tables, and the cold c54 arm against production,
are in `out/score.txt`. Cold c54 branched all 10 of production's children plus 15:54 (x 1.48).

## The cold parent against production

| node | production (warm down the chain): first DNN / cut bound / short | here: first DNN / cut bound / short | CPU s, production / here |
|---|---|---|---|
| depth 6 `...,0:32` (-F, cold) | 607,675,088 / 607,757,721 / 457,333 | 602,040,926 / 606,555,745 / 1,659,309 | 1,253 / 828 |
| depth 7 `10:39` (warm from the cold parent, production settings) | 607,692,496 / 607,809,133 / 405,921 | 606,819,754 / 607,277,924 / 937,130 | 1,165 / 636 |
| depth 7 `10:39`, cut phase to 8,000 iterations (`strong.par`) | | 606,819,754 / 607,798,661 / 416,393 | 4,855 |

Cold from `-F`, the depth-6 node has no warm DNN state and no inherited cut pool. Its first solve
gives up 5.6M below production's warm one. Its 1,000-iteration cut phase recovers to 1.66M short,
3.6x production's shortfall. 10:39 warm from it is 2.3x production's shortfall. In production, each
node down the chain inherited its parent's pool and state, seven 1,000-iteration phases deep. A
single 8,000-iteration phase at 10:39 recovers almost all of that: 416k against 406k.

Under the cold state the rule's own prediction moves with the slack. At 937k, c47's predicted count
(x 30) is 5-6 and c54's 11-12 (`out_predict.txt` recomputed at the arm's slack). The cold arms
realised 7 and 11. A test at a weaker state than production's asks a different question, which is
why the strong pair was added. The strong pair departs from production's settings at exactly one
node, 10:39's own cut phase (`strong.par`: `cut.giveup_depth = 8`, `cut.cp_iters = 8000`). Its
children are bounded with production's settings.

The depth-7 run (`run.sh kids`: every child of the cold depth-6 node, for the comparison with
production's depth-7 outcomes) and the r31 arms (cold `r31`, strong `sr31`) were started late. See
"Not measured".

## Not measured

- **Step 4, the survivors' subtrees**, by the kill rule. The size comparison above is a regression
  [e].
- **r31**: stopped part-way for the depth-9 run (cold `runs/r31`, strong `runs/sr31`); partial traces
  are committed. With 3,000 kicks the signal now puts r31 at 6 survivors (x <= 1) plus two at 1.02
  and 1.10 [m `out/effort_all.tsv`], worse than c47. The landscape README's 1 to 5 came from 30 or
  100 kicks.
- **Every depth-7 child's outcome under the cold parent** (`runs/kids`): stopped part-way for the depth-9 run: 31 of 54 bounded, 30 D and 10:29 B (production: 10:29 B too, among its 3 B and 51 D).
  Production's depth-7 children are 3 B and 51 D.
- **Bit-reproducibility of the arms**: not rechecked here. cognac's tai60b_lines found restarts
  reproduce bounds to the last digit printed.
- The build check (tai25b to depth 2, the experiment build against production's, node for node
  identical apart from CPU, `validate/`) used tai60b's settings without the cut phase, since the
  region file is for n = 60. The patch only changes line choice.

## Files

`run.sh` (parent, kids, arm, strong, resume, trim), `strong.par`, `cutstate.py`, `trimstate.py`
(from cognac), `supervise.sh`, `predict.py` (the landscape search at 10:39, `out/predict.tsv`),
`score.py` (`out/score*.{tsv,txt}`), `effort.py` / `effort_all.py` (`out/effort*.tsv`),
`subtrees.py` (`out/subtrees.txt`), `runs/<arm>/` (cmd, log, err, the node trace `nodes.txt`).
State files and DNN traces are not committed (`.gitignore`). Note: the landscape README's
"`-d 7 -o parent.state`" for the parent is off by one. `bnbdrv -d D` sets aside nodes deeper than D,
so the parent ran with `-d 6`, leaving its depth-7 children open.
