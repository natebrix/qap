# Exchange-optimality cuts (bound_lab agent `swapcuts`, Nov 15-16 2026)

STATUS: in progress. Stages 1 and 2 done; stage 3 (the gate) running. See `../COMMON.md`.
Marks: **[m]** measured (log named), **[e]** estimated. The box was shared (load 6-7 on 4
cores): iterations are the metric.

## The idea

A search run with `-u BESTKNOWN` proves that no permutation costs less than the incumbent.
If one did, a global optimum p* would, and no pairwise exchange improves p*. So every node's
relaxation may be restricted to completions that no exchange improves: optimality-based cuts,
which cut into the convex hull of the node's permutations (no inequality valid for all
permutations can). With facility i fixed at a, j free and b a free location, and
C(p) = sum A[k][l] B[p(k)][p(l)]:

    x_jb = 1  =>  delta = c0 + sum_{free k != j, free c != b} g_k(c) x_kc >= 0
    g_k(c) = (A[i][k]-A[j][k]) (B[b][c]-B[a][c]) + (A[k][i]-A[k][j]) (B[c][b]-B[c][a])
    c0 = (A[i][j]-A[j][i]) (B[b][a]-B[a][b]) + (A[i][i]-A[j][j]) (B[b][b]-B[a][a])
         + sum over the other fixed k of g_k(p(k))

    c0 Y[jb,jb] + sum g_k(c) Y[jb,kc] >= 0          (d m^2 cuts at depth d, m free)

i.e. in the relaxation: the conditional expectation, given j at b, of the gain of exchanging
j with the fixed i is nonnegative. `swaplib.py` has the derivation; `swapdrv.C` the C++.

## Stage 1: derivation and validity. PASS [m]

`logs/stage1.out` (`stage1.py`): the formula against brute-force cost differences on 600
random cases with asymmetric matrices, negative entries and nonzero diagonals: max error 0.
On 8 instances with n = 6 to 9 (random asymmetric, a 2 x 4 grid, a tai-b-like clustered one,
nug07, nug08, esc08a with 17,280 optima), every global optimum (at most 6 per instance)
against every cut of 1,499 nodes that contain it: 79,426 pairwise cuts and 8,892
fixed-pair cuts, smallest left-hand side 0 (never negative), and each left-hand side equals
the brute-force exchange gain.

## Stage 2: the exact ceiling. PASS (the kill rule, under 10% at depths 3-5, is far off) [m]

`exch.C` / `stage2.sh` / `stage2_summary.py`, `logs/stage2_*.out`, table in
`logs/stage2_summary.txt`. Seven n = 12 instances (nug12, had12, scr12, chr12a, tai12a,
tai12b, rou12), four random facility orders each; nodes = all children of a line under a
good-path prefix (G) or under the near-miss sibling of the path (W: the last prefix facility
at the wrong location of smallest node optimum). v_node = the node's optimum; v_ex1 = its
optimum over completions that no exchange of a fixed with a free facility improves (the
pairwise cuts' hull); "none" = no such completion within 5% of the optimum above v_node.
Off-path nodes only:

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

Where finite, the rise v_ex1 - v_node has a median of 1 to 4% of the optimum (per instance and
depth), several times a node's relaxation gap at these depths. So at the level of the integer
hull the cuts are very strong: most wrong nodes have no good exchange-optimal completion at
all. Exchanges among free facilities add nothing to "affected" (the node's best completion is
already optimal under them). nug14 was stopped (my enumerator's bound is too weak for 13!).

