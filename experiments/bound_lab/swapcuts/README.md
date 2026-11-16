# Exchange-optimality cuts (bound_lab agent `swapcuts`, Nov 15-16 2026)

See `../COMMON.md`. Marks: **[m]** measured (log named), **[e]** estimated. The box was shared
(load 6-7 on 4 cores): iterations are the metric, CPU seconds are indicative. Nothing under
`latest/` in the main checkout was touched; nothing was committed.

## Answer

@@ANSWER@@

## The idea

A search run with `-u BESTKNOWN` proves that no permutation costs less than the incumbent.
If one did, a global optimum p* would, and no pairwise exchange improves p*. So every node's
relaxation may be restricted to completions that no exchange improves: optimality-based cuts,
which cut into the convex hull of the node's permutations (no inequality valid for all
permutations can). Symmetry pruning stays valid (an automorphism maps optima to optima). With
facility i fixed at a, j free and b a free location, and C(p) = sum A[k][l] B[p(k)][p(l)]:

    x_jb = 1  =>  delta = c0 + sum_{free k != j, free c != b} g_k(c) x_kc >= 0
    g_k(c) = (A[i][k]-A[j][k]) (B[b][c]-B[a][c]) + (A[k][i]-A[k][j]) (B[c][b]-B[c][a])
    c0 = (A[i][j]-A[j][i]) (B[b][a]-B[a][b]) + (A[i][i]-A[j][j]) (B[b][b]-B[a][a])
         + sum over the other fixed k of g_k(p(k))

    c0 Y[jb,jb] + sum g_k(c) Y[jb,kc] >= 0          (d m^2 cuts at depth d, m free)

In the relaxation: given j at b, the conditional expectation of the gain of exchanging j with
the fixed i is nonnegative. Further members, all linear in Y and all in `swapdrv.C`: two fixed
facilities exchanged (`-fp`: on row 0 and, multiplied by x_ld, on every row); 3-cycles through
two fixed facilities and one free (`-rot`, d (d-1) m^2 cuts).

## Files

| | |
|---|---|
| `swaplib.py`, `stage1.py` | the derivation; stage 1 (`logs/stage1.out`) |
| `exch.C`, `stage2.sh`, `stage2_summary.py` | stage 2, exact enumeration (`logs/stage2_*.out`, `logs/stage2_summary.txt`) |
| `lpceil.py`, `lpceil_summary.py` | stage 2b, the cuts' ceiling as linear inequalities (`logs/lpceil.out`, `logs/lpceil_summary.txt`) |
| `swapdrv.C`, `build.sh` | the driver for stage 3 (its header comment lists the options); built against the worktree `/home/user/wt-swap` |
| `swapcuts.patch` | the worktree's diff to `latest/dnncut.C`, `latest/dnncut.h` (below) |
| `validate_small.py`, `validate_path_small.py`, `check_bks.py` | validity (`logs/validate_small.out`, `logs/validate_path_small.out`, `logs/check_bks.out`) |
| `runseq.sh`, `runseq2.sh`, `after*.sh`, `jobs_*.txt` | the job runners and lists as they were used (one process at a time); each log's first line is its command |
| `s3_summary.py` | the joint-rounds tables (`logs/s3_summary.txt`) |
| `logs/s3_*.log` | stage 3 runs |

`swapdrv` = the driver as of the first builds, `swapdrv2` = `swapdrv.C` as it is now (adds
`-rot`, `-elim`, the step caps `-ecap`, `-emu`). Logs whose first line says `./swapdrv2` were
made with the caps' defaults in force at the time; names ending in `m` have both caps, `c` the
first only, none: uncapped.

**What changed in C++** (`swapcuts.patch`, 87 lines, worktree only): (1) `dnn_cp_bound` keeps
`DNNCutParams::L_ext` in the CP solve's objective (the stock CP form drops it); (2) a hook
`DNNExtHook` in `dnn_cut_rounds`: the caller rewrites an external objective term between
rounds (after the separation, before the Polyak step; once more before the quiet round). The
cuts themselves live in `swapdrv.C`: Lagrangian form, objective term -sum mu_t G_t with
mu >= 0, so every value printed is a certificate of `dnn_bound` and is a valid bound on the
node's exchange-optimal completions for any multipliers.

## Stage 1: derivation and validity. PASS [m]

`logs/stage1.out` (`stage1.py`): the formula against brute-force cost differences on 600
random cases with asymmetric matrices, negative entries and nonzero diagonals: max error 0.
On 8 instances with n = 6 to 9 (random asymmetric, a 2 x 4 grid, a tai-b-like clustered one,
nug07, nug08, esc08a with 17,280 optima), every global optimum (at most 6 per instance)
against every cut of 1,499 nodes that contain it: 79,426 pairwise cuts and 8,892
fixed-pair cuts, smallest left-hand side 0 (never negative), and each left-hand side equals
the brute-force exchange gain. The driver repeats the check at every start (its `guard:`
line: every cut's linear form against brute-force gains on 50 random completions; 8.2 million
evaluations with rotations at tho40 depth 12, max error 0).

The best known solutions are exchange-optimal, as they must be if optimal
(`logs/check_bks.out`): smallest exchange gain +118 for tho40, +4 sko42, +6 sko49, +1,534
tai50b, +8 nug30, and **0 for wil50**: one exchange of its best known solution costs the same,
so wil50 has a second permutation of value 48,816 (not an automorphic image by this test alone;
not looked into further).

## Stage 2: the exact ceiling. PASS (the kill rule, under 10% at depths 3-5, is far off) [m]

`exch.C` / `stage2.sh` / `stage2_summary.py`, `logs/stage2_*.out`, table in
`logs/stage2_summary.txt`. Seven n = 12 instances (nug12, had12, scr12, chr12a, tai12a,
tai12b, rou12), four random facility orders each; nodes = all children of a line under a
good-path prefix (G) or under the near-miss sibling of the path (W: the last prefix facility
at the wrong location of smallest node optimum). v_node = the node's optimum; v_ex1 = its
optimum over completions that no exchange of a fixed with a free facility improves (the
pairwise cuts' logical content); "none" = no such completion within 5% of the optimum above
v_node. Off-path nodes only:

| prefix | depth | nodes | v_ex1 > v_node | of which none within 5% | with fixed-fixed exchanges too |
|---|---|---|---|---|---|
| G | 1 | 321 | 35% | 7% | 35% |
| G | 2 | 282 | 62% | 20% | 64% |
| G | 3 | 252 | 71% | 29% | 74% |
| G | 4 | 224 | 88% | 60% | 90% |
| G | 5 | 196 | 89% | 71% | 91% |
| W | 2 | 308 | 57% | 20% | 58% |
| W | 3 | 280 | 75% | 37% | 78% |
| W | 4 | 252 | 83% | 58% | 87% |
| W | 5 | 224 | 90% | 77% | 94% |

Where finite, the rise v_ex1 - v_node has a median of 1 to 4% of the optimum. Exchanges among
free facilities add nothing to "affected" (the node's best completion is already optimal under
them). nug14 was stopped (the enumerator's bound is too weak for 13!).

### Stage 2b: the ceiling of the cuts as linear inequalities on Y [m]

A relaxation sees only the linear inequalities, and a mixture of completions none of which
is exchange-optimal can satisfy them. That ceiling is an LP over a node's completions
(`lpceil.py`: min sum w_p C(p) subject to every cut, over the exact lifted hull): no bound
from these cuts in any relaxation can pass it. `logs/lpceil_summary.txt` (the footer inside
`logs/lpceil.out` miscounts infeasible LPs and is superseded), off-path nodes of tai10a,
rou10, scr10, had10, tai10b (depths 1-2) and of six n = 12 instances (depths 3-5):

| n | depth | m | nodes | logical (v_ex1 > v_node) | linear (v_lp > v_node) | LP infeasible | median rise, % of v_node | share of the logical rise kept |
|---|---|---|---|---|---|---|---|---|
| 10 | 1 | 9 | 90 | 32% | 32% | 1% | 0.61 | 0.55 |
| 10 | 2 | 8 | 80 | 52% | 52% | 4% | 1.51 | 0.69 |
| 12 | 3 | 9 | 35 | 71% | 71% | 9% | 1.13 | 0.61 |
| 12 | 4 | 8 | 70 | 83% | 83% | 29% | 1.80 | 0.60 |
| 12 | 5 | 7 | 68 | 88% | 88% | 56% | 1.22 | 0.70 |

So the linearisation loses little: every node the logical condition affects is affected, a
third to a half of the nodes at depths 4-5 are excluded outright, and where the rise is finite
the LP keeps 55 to 70% of it. Only 0.004 to 2% of a node's completions are exchange-optimal.
Whatever the DNN relaxation fails to deliver below is the relaxation's looseness, not the cuts'.

## Stage 3: in the DNN relaxation

@@STAGE3@@
