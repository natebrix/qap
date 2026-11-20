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
| (a) top LS line | running | | | | |
| (b) K = 2, CAP 30 | queued | | | | |
| (b) K = 4, CAP 30 | queued | | | | |

Production here has the same 1,017 nodes as line_signals' two production runs, with the same (depth,
fixes, bound, outcome, line) at every node [m: sorted trace comparison]; its CPU, 6,935 s, is inside their
6,931-7,390. So line_signals' (c) arm is reused rather than rerun. (The coordinator's own count, 987, is not
the trace's; line_signals saw the same.)

## Not tested yet

- T3 (cost per firing at n = 35, 40, 60 and the break-even at tai60b's depth-3 and depth-7 chain nodes).
- tai40b's arms.
- A tai50b partial tree.
