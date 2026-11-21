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

**Correction (Nov 20, 13:05 UTC): below depth 7 these runs used rule 9, not production's rule 2.** tai60b-261117 ran
`tai60b.par` from commit 316361fc (its `recipe.sh`), whose bottom block (depths 7 and down) is rule 2. Commit 4e1e1b6d
(Nov 18, "Rule 9 for rule 2 in every DNN recipe's bottom block") changed it to rule 9, and every run here used the
file at 0e6c2ecf. The trace's rule field is the parameter level (3 = the bottom block), not the rule number, and I did
not check the rule against the run. *Not affected:* the depth-7 arms (c47 5 against c54 9). Their depth-7 line was
forced and their depth-8 children were only bounded, and a child's outcome does not depend on the line it later
branches on (strong c54 matches production child by child). *Affected:* the depth-8 lines under c47 and every level
below. "Down the tree" compares **c47 with rule 9 below it** against **production's c54 with rule 2 below it**.
`experiments/rule9_bottom` (cognac) measured rule 9 at -19% nodes on tai50b, so part of the 10x at depths 9-18 may be
the rule, not the depth-7 line. Separating the two needs c54 + rule 9 (or c47 + rule 2) bounded below depth 7, which
was not done. In the depth-17 line test the "taken" lines are rule 9's. Rule 2's own lines were bounded as a third
arm (below).

**Answer (Nov 20, 11:21 UTC): the kill rule fired at depth 18. Under c47, 31 nodes branched there against
production's 49, which is over half, so by the rule fixed in advance c47 is no longer clearly ahead at depth 18. The
run stopped there; the subtree is not closed.** [m] Through depths 9-17, c47 branched 7 to 28 times per depth
against production's 56 to 199 (at most 15%). Depth 18 is where production's own tree was thinnest: 49 branching
in 5,203 nodes, just before its 19-22 bulge. c47's 31 come mostly from three depth-17 parents: the two ending in 8:49
and 8:57, which branched at 190k and 223k short of the cutoff, give 12 and 6, and the path node 22:31 gives 7.
Their lines are c16 (17 of the 31) and r9 (6).

Totals for depths 9-18 [m]: **c47 5,486 nodes, 142 branching, 143,075 CPU s on this box; production under c54
55,267 nodes, 1,218 branching, 1,713,605 CPU s on the cloud boxes.** That is 10.1x fewer nodes. Scaled to
production's boxes (CPU / 0.53 = 270,000 s), c47 used 16% of production's CPU over the same depths [e]. Production's
whole subtree under 10:39 was about 114,000 nodes and 770 core-hours. Its depths 19 and below (58,700 nodes) are not
compared, because c47's are not bounded.

**Path check** (`pathcheck.py`, `out/pathcheck.txt`): the best known solution's path stays in the tree through
depth 18. Each of its 12 nodes, depths 7-18, branched, with bounds from 607,798,661 to 608,139,870, all below the
cutoff 608,215,054. The incumbent was best known + 1, so every fathom in the bounded part proves no solution better than the best known there.
The "other" outcomes at depths 15, 16 and 17 are one QPB fathom each (`Q`): the QP bound alone was above the cutoff
(608,292,252, 608,513,649 and 608,252,789).

**To resume** (for example, to see whether c47 has its own 19-22 bulge): `python3 down.py 19`, from depth 18's 44
set-aside states in `runs/down/d18/out/` (503 MB, in this container only, not committed). That means overriding
the kill rule, which is Nathan's call. The run ran from 23:22 to 09:22 under the coordinator's 10-hour limit, then from
09:25 under the 15:30 deadline. Three container restarts (about 00:00, 03:25 and 04:10 UTC) cost only the chunks in
progress, which were rerun.

### The 2x2 at depth 9: the line and the rule, separated (Nov 20, 13:10-19:15 UTC; the coordinator's request)

**The two effects are separate and about multiplicative. At depth 9, c47 instead of c54 cuts the branching children to a
third under either rule, and rule 9 instead of rule 2 halves them under either line.** [m] Rows: the depth-7 line
(from the strong state, 10:39 416k short). Columns: the rule that chose the depth-8 lines. Each cell is depth 9's
nodes, branching (B), shortfall sum and CPU:

| | rule 2 below (production's) | rule 9 below (0e6c2ecf's par) |
|---|---|---|
| **c54** | 520 nodes, **56 B**, 80,387 s (production, cloud boxes); depth-8 lines c2 at all 10 | 468 nodes, **27 B**, 4,583k, 24,609 s; lines r3 (6), c53, c35, c59 |
| **c47** | 260 nodes, **19 B**, 2,469k, 21,055 s; lines c54 at all 5 | 260 nodes, **9 B**, 1,152k, 9,002 s; lines c10, c53, c59, c59, r9 |

The line effect is 19/56 = 0.34 under rule 2 and 9/27 = 0.33 under rule 9. The rule effect is 27/56 = 0.48 under c54
and 9/19 = 0.47 under c47. So the 6x between production (56) and the c47 descent (9) at depth 9 is about 3x from the
depth-7 line and 2x from the rule. Per depth-8 node, c47 + rule 2: 47:47 8 B, 46:47 9, 29:47 2, 31:47 0, 45:47 0
(rule 9: 3, 3, 2, 0, 1). Rule 2's c54 costs most under the two survivors farthest from the cutoff.

*Check of c54 + rule 2 against production:* two of c54's depth-8 parents, rebranched from this run's state under
`rule2.par`, pick production's line (c2) and leave production's counts exactly: 25:54 5 B (production 5) and 16:54
4 B (production 4). So production's 56 is a like-for-like cell here, on the evidence of 2 of its 9 parents.

How: c54 + rule 9 bounds depth 9 from the strong c54 arm's set-aside states (`lvl.py c54r9`; the strong arm's depth-8
lines were chosen by rule 9). c47 + rule 2 and the c54 check need the depth-8 nodes branched again, so 10:39 was
rebounded with its depth-8 children set aside (`run.sh strong7 c47|c54`: 10:39's bound 607,798,660.8327, identical
to the strong arms'). Then each depth-8 node was bounded under `rule2.par` with depth 9 bounded and depth 10 set aside
(`d9arm.sh`, `runs/d9x/`). The depth-8 nodes' own bounds differ from the strong arms' by up to a few hundred
(47:47 607,885,595 against 607,885,434). Their warm start came from a checkpoint (fp32, `ship.bits = 32`); in the arm
it was the parent's in-memory state. Depth 10 and below were not compared.

### Rule 2 against rule 9 over whole levels: every branched node at depths 16 and 17 (Nov 20, 20:40-23:05 UTC; the coordinator's request)

**These are the whole levels, not a selection: all 11 branched depth-17 nodes and all 7 branched depth-16 nodes of the
c47 descent.** [m] At depth 17, rule 2's lines leave 11 branching children against rule 9's 31 (shortfall sums 606k
against 5,144k; CPU 12,855 s against 22,961). That total rests on two nodes. At 8:49 and 8:57, rule 9's c30 leaves 18
and rule 2's r2 none. **On the 8 depth-17 nodes that were not chosen before, rule 2 leaves 8 against rule 9's 6:
slightly worse, not better.** So the earlier 3 against 25 was mostly the selection, though not regression to the mean:
rule 9 does badly at a few nodes, and those happen to be where it was looked at. At depth 16 the two are close, rule 2
9 against rule 9 11; on the path, rule 2 3 against 2. Rule 2's count is 0 at all 18 nodes (15 measured here, the
3 depth-17 parents before), so rule 2 chooses by its U-sum tie-break throughout. At the depth-16 node 4:23, rule 2's line
is rule 9's (r8). Its arm reproduces the descent's children exactly (1 branching, 86k), a reproducibility check for free.
Rule 9's numbers are the descent's (`runs/down`). Rule 2's lines were bounded from the parents' states, rebuilt by five
chains (`all11.sh`, `runs/rebuild/R1`-`R5`: each rebuilt node has the run's bound and line exactly), with `rule2.par`
and the `latest_dbg` build (`armD.sh`, `runs/scan/s16_*`, `s17_*`). Score: `all11_score.py`, `out/all11.txt`.

depth 17: 11 branched nodes in the descent (the whole level)

| node | path | rule 9: line, B, shortfall sum, CPU s | rule 2: line, B, shortfall sum, CPU s |
|---|---|---|---|
| 46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:23,8:29 | no | r5: **0**, 0k, 309 | r2: **1**, 1k, 356 |
| 46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:29,8:49 | no | c30: **12**, 3,350k, 4,881 | r2: **0**, 0k, 117 |
| 46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:29,8:57 | no | c30: **6**, 962k, 3,772 | r2: **0**, 0k, 153 |
| 46:47,36:53,32:40,3:35,11:59,31:10,47:38,51:30,50:6,9:55 | no | r13: **0**, 0k, 733 | r8: **0**, 0k, 803 |
| 47:47,31:10,46:38,36:53,27:40,22:31,11:59,3:35,51:30,48:6 | no | c11: **1**, 30k, 1,515 | c54: **2**, 57k, 1,950 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:36,13:59,12:0,9:29 | no | c6: **2**, 265k, 1,550 | c23: **0**, 0k, 709 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:36,13:59,12:0,9:55 | no | r4: **1**, 12k, 677 | r8: **1**, 39k, 573 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,48:6,22:31 | yes | c26: **7**, 391k, 5,914 | c54: **3**, 231k, 2,749 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,48:6,25:31 | no | c43: **1**, 56k, 2,152 | c54: **2**, 88k, 2,908 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,50:6,9:0 | no | c4: **0**, 0k, 712 | r8: **1**, 121k, 1,380 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,50:6,9:55 | no | r13: **1**, 79k, 747 | r8: **1**, 69k, 1,158 |
| **all** (11 nodes with both) | | rule 9: **31**, 5,144k, 22,961 | rule 2: **11**, 606k, 12,855 |
| **path** (1 nodes with both) | | rule 9: **7**, 391k, 5,914 | rule 2: **3**, 231k, 2,749 |
| **off** (10 nodes with both) | | rule 9: **24**, 4,753k, 17,048 | rule 2: **8**, 375k, 10,106 |

depth 16: 7 branched nodes in the descent (the whole level)

| node | path | rule 9: line, B, shortfall sum, CPU s | rule 2: line, B, shortfall sum, CPU s |
|---|---|---|---|
| 46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:23 | no | r8: **1**, 86k, 593 | r8: **1**, 86k, 542 |
| 46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:29 | no | r8: **2**, 414k, 728 | c23: **0**, 0k, 172 |
| 46:47,36:53,32:40,3:35,11:59,31:10,47:38,51:30,50:6 | no | r9: **1**, 24k, 1,276 | c54: **1**, 16k, 1,395 |
| 47:47,31:10,46:38,36:53,27:40,22:31,11:59,3:35,51:30 | no | c6: **1**, 38k, 1,667 | c54: **1**, 54k, 1,326 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:36,13:59,12:0 | no | r9: **2**, 215k, 941 | r8: **1**, 88k, 588 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,48:6 | yes | c31: **2**, 224k, 2,182 | c54: **3**, 259k, 2,536 |
| 47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,50:6 | no | r9: **2**, 122k, 1,816 | c54: **2**, 176k, 2,115 |
| **all** (7 nodes with both) | | rule 9: **11**, 1,123k, 9,202 | rule 2: **9**, 678k, 8,675 |
| **path** (1 nodes with both) | | rule 9: **2**, 224k, 2,182 | rule 2: **3**, 259k, 2,536 |
| **off** (6 nodes with both) | | rule 9: **9**, 898k, 7,020 | rule 2: **6**, 420k, 6,139 |


### Where does rule 2 overtake rule 9? A scan at depths 11, 13, 15 (Nov 20, 18:15-20:50 UTC; the coordinator's request)

**It does not, at these nodes. Over the nine nodes rule 2's lines left 24 branching children and rule 9's left 21.
Node by node rule 2 is worse at 2, better at 1 (d15_a: 1 against 2) and tied at 6. It costs more: 18,704 CPU s
against 11,641, most of it in cut phases. Rule 2's count, the children QPB + U put at the cutoff, is 0 at every scanned
node, and also at the three depth-17 parents where rule 2 won (8:49, 8:57, 22:31; `count17.sh`).** [m] So rule 2 never
chooses by its count in this subtree: its line is the tie-break's, the largest sum of U over the line. Whether rule 2
wins where its count is not zero cannot be tested here. Between depth 8 (rule 9 better: 9 against 19 under c47,
27 against 56 under c54) and depth 17 (rule 2 better: 3 against 25), the scanned depths show neither ahead.

*Selection bias:* at each depth the three nodes are those where rule 9's line did worst in the c47 descent (the most
branching children, ties by shortfall sum; `out/scan_select.txt`). That favours rule 2, and rule 2 still did not
win. The descent's depths 11-15 had few branching nodes (13, 28 and 8), so "worst" means 1 to 3 branching children.
Rule 2 chose c54 at 6 of the 9 nodes and r8 at 2 (its first choice in the U-sum order). Rule 9's lines are the
descent's (`runs/down`). Rule 2's lines were bounded from the same parent states, rebuilt where needed
(`rebuild.py`: every rebuilt node has the run's bound and line exactly), with `rule2.par` and the `latest_dbg` build
(`maxf_debug.patch`: a print of rule 2's counts under `QAP_MAXF_DEBUG=1`, the search unchanged; `armD.sh`,
`scan.sh`, `scan_score.py`, `out/scan.txt`). QPB below cutoff = (cutoff - the node's QPB bound) / cutoff.

| node (depth) | path | QPB below cutoff | rule 9: line, D / C / B, shortfall sum, CPU s | rule 2: line, D / C / B, shortfall sum, CPU s | rule 2 count (chosen side / other side) |
|---|---|---|---|---|---|
| d11_a: 46:47,36:53,32:40,3:35 (11) | no | 0.171 | r11: 46 / 0 / **3**, 500k, 1,488 | c54: 40 / 4 / **5**, 896k, 3,770 | 0 / 0 |
| d11_b: 46:47,36:53,34:40,3:35 (11) | no | 0.171 | r11: 46 / 0 / **3**, 428k, 1,540 | c54: 40 / 6 / **3**, 326k, 3,480 | 0 / 0 |
| d11_c: 47:47,31:10,46:38,36:53 (11) | yes | 0.170 | c40: 46 / 0 / **3**, 398k, 1,742 | c54: 45 / 1 / **3**, 278k, 2,213 | 0 / 0 |
| d13_a: 47:47,31:10,46:38,36:53,32:40,3:35 (13) | yes | 0.160 | r11: 44 / 0 / **3**, 409k, 1,253 | c54: 39 / 4 / **4**, 614k, 2,662 | 0 / 0 |
| d13_b: 47:47,31:10,46:38,36:53,34:40,3:35 (13) | no | 0.160 | r11: 44 / 0 / **3**, 265k, 1,240 | c54: 39 / 5 / **3**, 199k, 2,287 | 0 / 0 |
| d13_c: 46:47,36:53,32:40,3:35,11:36,13:59 (13) | no | 0.148 | r12: 43 / 3 / **1**, 90k, 979 | r8: 45 / 1 / **1**, 116k, 871 | 0 / 0 |
| d15_a: 46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55 (15) | no | 0.133 | r4: 41 / 1 / **2** (1 Q), 346k, 861 | r8: 40 / 3 / **1** (1 Q), 55k, 830 | 0 / 0 |
| d15_b: 47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30 (15) | yes | 0.141 | c6: 42 / 1 / **2**, 201k, 1,659 | c54: 41 / 1 / **3**, 232k, 1,560 | 0 / 0 |
| d15_c: 47:47,31:10,46:38,36:53,27:40,22:31,11:59,3:35 (15) | no | 0.146 | r51: 44 / 0 / **1**, 51k, 879 | c54: 42 / 2 / **1**, 55k, 1,031 | 0 / 0 |

### The depth-17 line test (Nov 20, 11:25-12:58 UTC; the coordinator's follow-up to the depth-18 stop)

**By the decision rule fixed in advance (the signal helps if its three lines leave at most 12 branching children,
against the taken lines' 25), the signal helps: its lines left 8 (5, 1 and 2) [m]. That is a weaker result than it
looks.** At the two parents off the chain (8:49 and 8:57), the signal predicted 0 survivors for every line. Their own
best completions are 3.4 to 3.65 times their slack above the best known, the regime where the landscape README found
the threshold fails. Its "top line" r51 won only on the tie-break: off 0.002 against 0.002 for the taken c30, which
was the signal's third line. So at those two parents the signal did not tell r51 from c30, and r51's 5 and 1 against
c30's 12 and 6 is one draw, not a prediction borne out. At the path node 22:31 the signal made a real call: r13 with 1
predicted survivor against the taken c26's 2. r13 left 2, c26 left 7. The parents are the down run's, rebuilt
(`rebuild.py`: their ancestors bounded one at a time from the depth-9 states). Every rebuilt node, depths 10 to 17,
has the run's first-solve bound, cut bound and line to the last digit printed (`runs/rebuild/*.out`). The arms are
`arm17.sh`: the parent forced to the signal's line (`QAP_FORCE_LINE=17:LINE`), production's settings, its 43
children bounded, depth 19 set aside. Signal: `signal_node.py`, every free line at the parent's own slack, 30 and 300
kicks (`out/signal_*`). Scores: `score17.py`, `out/score17.txt`.

| parent (slack) | line | signal: predicted survivors, 300 kicks | D / C / B | branching shortfall: sum / largest | CPU s |
|---|---|---|---|---|---|
| 8:49 (190,223) | taken c30 (rule 9) | 0 (rank 3 of 86) | 17 / 14 / **12** | 3,350k / 1,073k | 4,881 |
| | signal r51 | 0 (rank 1, on the tie-break) | 33 / 5 / **5** | 543k / 179k | 2,048 |
| 8:57 (223,279) | taken c30 (rule 9) | 0 (rank 3) | 24 / 13 / **6** | 962k / 276k | 3,772 |
| | signal r51 | 0 (rank 1, on the tie-break) | 35 / 7 / **1** | 706k / 706k | 1,731 |
| 22:31, path (85,768) | taken c26 (rule 9) | 2 (rank 33) | 33 / 3 / **7** | 391k / 157k | 5,914 |
| | signal r13 | 1 (rank 1) | 41 / 0 / **2** | 281k / 158k | 1,651 |
| **total** | taken (rule 9) | 2 | **25** | 4,703k | 14,567 |
| | signal | 1 | **8** | 1,530k | 5,430 |
| 8:49 | rule 2: r2 | 0 (rank 5) | 29 / 0 / **0** (14 Q) | 0 | 117 |
| 8:57 | rule 2: r2 | 0 (rank 5) | 37 / 0 / **0** (6 Q) | 0 | 153 |
| 22:31 | rule 2: c54 | 2 (rank 41) | 36 / 4 / **3** | 231k / 117k | 2,749 |
| **total** | rule 2 | 2 | **3** | 231k | 3,019 |

**The control (the coordinator's, Nov 20 13:00): rule 2, production's actual rule at these depths, beats both.** Rule
9's line is the taken line at all three parents, so it was not rerun. Rule 2's lines (`rule2.par`: the run's par file
at 316361fc; `arm17.sh` with LINE `-`) left 3 branching children (0, 0, 3), against the signal's 8 and rule 9's 25, in
3,019 CPU s against 5,430 and 14,567 [m]. At the two off-chain parents rule 2's r2 fathoms everything in the first
solve or by QPB alone (Q). The signal's own excesses show why: r2's children are at x 3.4-3.8 for two of them, then
74 to 155 times the slack. The signal's count is 0 for every line there and its tie-break is near 0 for all, so it
does not see this; the size of the far children does. So at depth 17 **rule 9's line is the bad one, more than the
signal's line is the good one**: rule 2 > signal > rule 9 at both off-chain parents. At the path node the three are
rule 2's c54 3, the signal's r13 2 and rule 9's c26 7 branching children, with shortfall sums 231k, 281k and 391k.
This is three parents, one level.

The signal's top three lines, 300 kicks: at 8:49 and at 8:57, r51, c38, c30 (0 predicted each; every line has 0); at
22:31, r13, c0, c36 (1 each). Only the top line was bounded at each parent. The second and third were not, and since
the taken line was never the signal's top line, no second line was needed under the coordinator's rule. One branching
child of 8:57's r51 is far from the cutoff (51:9, 706k short), so 1 child there is not necessarily a light subtree.
Below depth 18 nothing is measured.

*Earlier status, as of Nov 20 09:22 UTC (the 10-hour limit, depths 9-17; superseded above):* Depths 9 to 17 under
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

How it ran: level by level under the strong c47 arm's survivors, with `tai60b.par` as it stands at 0e6c2ecf, whose bottom block is rule 9
(production ran rule 2 there: see "Correction" above), choosing the line throughout (no forced lines): depth D is bounded from the states depth D-1 set aside, and depth D+1 is set aside
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
| 18 | 469 | 391 / 47 / 31 | 22,961 | 5,144k / 1,073k | 5,203 / 49 / 39,479 |
<!-- /down-table -->

### One level deeper (depth 9; the coordinator's first follow-up, Nov 19 22:25-23:10 UTC)

**Under c47, depth 9 has 9 branching nodes against production's 56 under c54, so c47 is materially
better one level down. The kill rule, fixed in advance at 28 (half of 56), is not met, and per that
rule I stop here.** [m `d9score.py`, `out/d9score.txt`] The 260 depth-9 nodes under the strong c47
arm's five survivors were bounded from `runs/sc47/out.state`. Production's settings applied
(`tai60b.par` at 0e6c2ecf), with its own line rule at depth 8, rule 9 (no forced line; production ran rule 2 there, see "Correction"), and depth 10 set aside
(`-d 9`; `run.sh deep`, `splitstate.py`, 4 processes, `runs/d9_*`). Production's figures are from
the tai60b-261117 trace.

| depth-8 node | line (rule 9) | depth 9: D | C | B | branching children's shortfalls | CPU s |
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

## Rule 9's failures

See [`failures.md`](failures.md) (Nov 20-21): every failure in this subtree is a column whose children are one
tight facility cluster, chosen for a single large entry it shares with a row; the good line was never one of rule 9's two
candidates; the relaxation's mass does not name the surviving cluster but its non-zero support does; "support" fails its
decision rule.

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
