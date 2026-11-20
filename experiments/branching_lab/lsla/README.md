# lsla: the landscape signal shortlists lines, the DNN lookahead chooses (branching_lab)

Nathan's idea (Nov 19): at nodes where much work is at stake, use the landscape signal (local search on
every child) to shortlist lines, and the cheap, sloppy DNN lookahead to choose among them.

## Answer (in progress; updated as arms finish)

**Not yet known.** What is in hand [m]:

- The signal in C costs **35 CPU s per firing at n = 60** with 30 kicks per child, and 355 s with 3,000
  more kicks on the cells near the threshold. The DNN lookahead it would feed costs hours per firing at that
  size, so the signal is cheap next to it.
- **At tai60b's depth-7 chain node the signal ranks the lines differently depending on how hard it
  searches**, and at either effort its top lines are lines nobody has bounded (c10, c59, c35, c38, r22 after
  refinement; r31, c10, c11, c47 at 30 kicks). A K = 2 shortlist is therefore unverified lines plus rule 9's.
- **As a shortlist on tai35b it agrees with the DNN lookahead**: at 4 of the 5 nodes where (c) fired and
  the gate passes, the signal's first line is the lookahead's choice, for 6 CPU s against 450 s. On tai40b
  the slack below the root is too small for it to tell lines apart (T3).
- **The signal alone (mode a) on tai35b: 983 nodes against production's 1,017 and (c)'s 780**, 6,352 CPU s
  gross (production 6,935), 396 s of it the signal. It matches (c) at the root and depth 1 and loses (c)'s
  depth-2 gains, because the refined signal prefers c5 to c23 (T2).
- tai35b's production arm reproduces line_signals' tree exactly, so (c)'s numbers are reused (below).
- T1 (the lookahead on the shortlist at the depth-7 node) and the other tai35b arms are running.

## The rule (`lsla.patch`, against the experiment build: cognac 0e6c2ecf + `experiment_build.patch`)

`QAP_LS=MODE:K:DMIN:DMAX`; unset, nothing changes. At a node of depth DMIN..DMAX whose DNN bound did not
fathom it:

1. **Gate.** s = cutoff - bound. The node's own best completion from the incumbent permutation (a
   2-exchange descent and 300 kicks) has excess ex_node over the incumbent; the rule goes on only if
   ex_node <= s (the node could hold a completion the bound cannot cut). Firings and passes are logged.
2. **Signal.** For every free cell (r, v), ex(r, v) = the best completion with r:v forced, from the node's
   best completion and from the incumbent adjusted to the fixes, each with one exchange to force r:v, a
   descent, then 30 kicks (6 random exchanges + descent, kept if better) on the better of the two. Cells with
   1 < ex/s <= 2.5 get KICKS2 more kicks (`QAP_LS_SEARCH=KICKS:KICKS2:XHI:GATEKICKS`, default 30:0:2.5:300).
   Taillard's delta table: O(n^3) per node, O(m^2) per exchange, so a child costs O(m^2) to set up.
3. **Shortlist.** Lines by predicted survivors (children with ex <= s), ties by the chain model's
   off-chain mass; the top K, plus rule 9's row and column, plus `QAP_LS_ADD` lines, de-duplicated.
4. **Choice.** Mode a: branch on the top line, no lookahead. Mode b: the existing `QAP_LA` lookahead (warm,
   pooled DNN solves of CAP iterations per child, `QAP_LS_CAP=30[,150]`) over the shortlist; ties keep the
   LS order. Mode c is plain `QAP_LA`.

One log line per firing (`QAP_LS_LOG`); `QAP_LS_DUMP` writes every cell's ex.

**Checks [m].** `QAP_LS` unset: tai25b's trace to depth 2 equals linetest's experiment build in every field
but CPU, 98 nodes (`validate/nodes.unset.txt`). `QAP_LA=2:20:0:1` (I changed its code): the same nodes and
the same lookahead log as a plain experiment-build binary (`validate/nodes.la_*.txt`). Modes a and b run
on tai25b; the signal's permutation costs match the `.sln`.

## T0: the C signal at tai60b's depth-7 chain node (`t0.C`, `t0out/`)

Node `44:45,43:58,1:33,6:24,7:19,0:32,10:39`, slack 405,921, n = 60, m = 53.

| search | c54 | c47 | r31 | CPU s | output |
|---|---|---|---|---|---|
| 30 kicks | 9 | 2 | 2 | 34.6 [m] | `t0out/k30_s0.txt` |
| + 1,000 kicks on cells with 1 < x <= 2.5 | 9 | 3 | - | 108 [m] | `t0out/k30_r1000_x2.5_s0.txt` |
| + 3,000 kicks on those cells (284 cells) | 9 | 4 | 6 | 355 [m] | `t0out/k30_r3000_x2.5_s0.txt` |
| truth (linetest, bounded) | 9 | 5 | - | hours | linetest |

- At 30 kicks the x values match the Python signal to two decimals (c54: 0.16 0.26 0.39 0.41 0.61 0.66
  0.66 0.79, then 1.13 1.48; c47: 45:47 1.39, 29:47 1.79, 31:47 2.44). The C is 2x faster than Python (69 s).
- With 3,000 more kicks c47's children reach 29:47 0.71, 45:47 0.74, 31:47 1.07: linetest's best of three.
  c54 stays right at any effort; c47 is under-predicted (4 of a true 5) even at 3,000.
- **Ranking is unstable with effort.** Top lines by predicted count: 30 kicks r31, c10, c11, c47 (2 each);
  after refinement c10, c59, c35, c38, r22 (3 each), then c47 (4). Neither top list has been bounded.

## T1: the lookahead on the shortlist at that node (`run.sh t1`, `runs/t1_cap*/`)

Running (two processes, CAP 30 and CAP 150, each from linetest's strong parent state at depth 6, which
reproduces linetest's cut bound 606,555,744.87 [m]). Shortlist: K = 2 by the signal, rule 9's two lines,
c47 and c54 (truth c47 5, c54 9).

## T2: whole trees (`t2.sh`, `t2queue.sh`, `summ.py` -> `out/t2_*.md`)

tai35b, coordinator local mode, 2 slots, incumbent the optimum + 1, line_signals' `tai35b.par`, windows
depth 0-2 as (c). **Kill rule, fixed before the runs:** (b) has no legs on tai35b unless its gross CPU is
more than 7% (the noise) below production's, or it reaches <= 780 nodes (c's) with lookahead CPU under
1,402 s (half of c's 2,803 s); and it must close at the optimum with no open node.

| arm | nodes | CPU s | rule CPU s | CPU s net of the rule | source |
|---|---|---|---|---|---|
| production | 1,017 | 6,935 | - | 6,935 | [m] `t2/tai35b_prod` |
| (c) `QAP_LA=4:30:0:2` | 780 | 8,900 | 2,803 (6 nodes) | 6,097 | [m] line_signals `runs/tai35b_la4p` |
| (a) top LS line (refined signal) | 983 | 6,352 | 396 (5 firings, 1 gate fail) | 5,956 | [m] `t2/tai35b_a` |
| (b) K = 2, CAP 30 | queued | | | | |
| (b) K = 4, CAP 30 | queued | | | | |
| (b) K = 2, CAP 30, signal without refinement (`SEARCH=30:0:2.5:300`) | queued | | | | |

Production here has the same 1,017 nodes as line_signals' two production runs, with the same (depth,
fixes, bound, outcome, line) at every node [m: sorted trace comparison]; its CPU, 6,935 s, is inside their
6,931-7,390. So line_signals' (c) arm is reused rather than rerun. (The coordinator's own count, 987, is not
the trace's; line_signals saw the same.)

## T3: the signal's cost, and its ranking against the lookahead's (`t3rank.py`, `t3out/`, `out/t3_rank.md`)

**Cost per firing [m]** (t0, one core, the gate's 300 kicks plus the table):

| node | n, m | 30 kicks | + 3,000 kicks near the threshold | for comparison |
|---|---|---|---|---|
| tai35b root | 35, 35 | 6.4 s | +115 s (236 cells) | the root's DNN node 590 s; (c)'s lookahead 468 s |
| tai40b root | 40, 40 | 10.7 s | +91 s (144 cells) | the root's DNN node 1,277 s |
| tai60b depth 7 | 60, 53 | 34.6 s | +320 s (284 cells) | bounding the depth-7 node 10:39 itself under `strong.par` about 4,900 s; a depth-9 node 21-52 s (linetest `d9score.txt`: 52 children per depth-8 node in 1,088-2,692 s) |

**Ranking at the 12 nodes where (c) fired in line_signals [m `out/t3_rank.md`].** The signal (30 kicks)
was run at each node, with that arm's slack, and the 8 lines the lookahead judged put in the signal's order:

- **tai35b, 5 nodes where the gate passes: the signal's first line is the lookahead's choice at 4**, and the
  root's order is the lookahead's order (c26, r19, r18; the DNN's 23M, 43M, 46M short). At the fifth
  (depth 2, 18:12 19:26) it puts c5 (rule 9's line there in production) first and the lookahead's c23
  second, at 4 survivors against 4. About 6 CPU s per node against the lookahead's 430-470 s.
- **tai40b: the signal has nothing to say below the root.** The slack is small (9k-239k below the root, against
  1.0-2.5M on tai35b): the gate fails at 3 of 6 nodes, and where it passes every line has 1 predicted
  survivor (the optimum's own child) and no off-chain mass, so the order is a tie. At the root it ranks
  production's r23 first and (c)'s c27 7th of 8; (c)'s tree was the worse (829 nodes against
  production's 731), so there it sided with the better rule.
- The incumbent permutation the signal starts from is the optimum (`QAP_LS_PERM` = the `.sln`), as in every
  arm here. On a node off the optimum's chain it has to find the good completions itself.

## linetest's depth 9 (fetched from `claude/cool-babbage-ilnogj`, linetest `out/d9score.txt`)

Under c47, depth 9 has 9 branching nodes (260 bounded, 9,002 CPU s on this box) against production's 56
under c54 (520 bounded). c47, which the refined signal (4 survivors) and the strong bound (5) both prefer
to c54 (9), stays the better line a level down. That supports the signal's ranking of those two lines; it
says nothing about the lines the signal ranks above c47 (c10, c59, c35, c38, r22), which no one has bounded.

**Nodes by depth [m]** (nodes / branched; from the node traces):

| arm | 1 | 2 | 3 | 4 | 5 | depth 0-2 lines |
|---|---|---|---|---|---|---|
| production | 35/3 | 102/3 | 99/4 | 128/5 | 155/2 | r19; r18 x3; c5 x3 (rule 9) |
| (a) | 35/2 | 68/3 | 99/4 | 128/5 | 155/2 | c26; r19, r18; c5, c5, and rule 9 where the gate failed |
| (c) | 35/2 | 68/3 | 99/2 | 64/2 | 62/2 | c26; r19, r18; c23, c23, c22 |

Below depth 5 the three trees have the same counts. **Mode (a) gets (c)'s root and depth 1 (34 nodes fewer
than production at depth 2), and none of (c)'s gain below**: at depth 2 it picks c5, rule 9's line, and from
depth 3 its tree is production's. (c)'s 203 further nodes saved come from its depth-2 lines c23, c23, c22.
The refinement is what loses them: at 18:26 19:12 the 30-kick signal ranks c23 first (T3), the refined one
(3,000 more kicks near the threshold) puts c5 first and c23 outside its top 12. At 18:8 19:26 the gate
fails (ex_node 1.78M against s 165k), so no signal fires where (c) picked c22.

Mode (a)'s CPU is 8% below production's (6,352 against 6,935; production's second run in line_signals was
7,390), at the edge of the 7% noise, from one run; its rule cost 396 s, nearly all refinement.

**Prediction for (b), written before its runs finished:** under the refined signal its K = 2 and K = 4
shortlists at 18:26 19:12 lack c23 (top 4: c5, c33, c6, r23), and it does not fire at 18:8 19:26, so (b)
cannot reach (c)'s depth-3-5 savings at two of the three depth-2 nodes. The extra arm without refinement
tests whether the 30-kick shortlist lets the lookahead find c23.

## Not tested yet

- The break-even at tai60b's depth-3 and depth-7 chain nodes (T3's second half).
- tai40b's arms.
- A tai50b partial tree.
