# Exchange-optimality cuts (bound_lab agent `swapcuts`, Nov 15-16 2026)

See `../COMMON.md`. Marks: **[m]** measured (log named), **[e]** estimated. The box was shared
(load 6-7 on 4 cores): iterations are the metric, CPU seconds are indicative. Nothing under
`latest/` in the main checkout was touched; nothing was committed.

## Answer

**The cuts are valid and they do cut into the hull, but in the DNN relaxation they are worth
too little: closed as a material lever.** wil50 fails the gate outright; tho40 sits just under it.

| | measured [m] | the gate |
|---|---|---|
| wil50 node 13, on top of the triangle rounds, equal iterations | **+1.4** (scatter 5) | +10 |
| wil50 node 10, the same (3 rounds) | **+2.7** on the best round (-32, -12 in the two rounds after the triangles' own overshoot) | +10 |
| tho40 far-miss nodes d4_d, d5_d, on top of the triangle rounds, equal iterations | **+70 to +290 a round, median +135**; d4_d fathoms at 1,350 iterations where the triangles alone end 196 short | +150 |
| tho40 gate node, exchange cuts alone against the plain solve | +150 (the triangle phase: +1,327 and a fathom) | |

1. **Valid** [m]. The formula matches brute force exactly (asymmetric matrices, diagonals);
   3,248 certified values on n = 12 nodes never pass the enumerated optimum over
   exchange-optimal completions (or the optimum, on path nodes), and no value on the tho40 path
   to the best known solution passes 240,516 (stages 1, 3e). The best known solutions of tho40,
   wil50, sko42, sko49, tai50b and nug30 are exchange-optimal; wil50's has one exchange of gain
   0, i.e. a second permutation of value 48,816.
2. **At the level of the integer hull the idea is very strong** [m]: on n = 12 instances 71 to
   90% of off-path nodes at depths 3-5 lose their best completion and 29 to 77% have no
   exchange-optimal completion within 5% of the optimum (stage 2). The linear inequalities on Y
   keep most of that over the exact lifted hull: the same nodes are affected, 9 to 56% of them
   are excluded outright, and 55 to 70% of the finite rise survives (stage 2b).
3. **The DNN relaxation realises about a tenth of the violation it shows** [m]. 2 to 4% of the
   d m^2 cuts are violated at converged plain iterates, concentrated on one or two fixed facilities
   (88 and 98% of the mass on facility 8 at corner 0 at two of the three tho40 nodes). The sum
   of violations is 810 to 1,413 cost units at tho40 depths 4-5 and the certified gain +110 to
   +217; it is 15 units at wil50 node 13 (gain +1.4) and 68 at node 10. wil50's dense, random
   flows leave almost nothing for an exchange to gain in expectation; tho40's sparse,
   heavy-tailed flows leave something, an order of magnitude less than the regional triangles
   give there (+1,300 to +1,900).
4. **What it would be worth** [e]. wil50 (and, unmeasured, sko42/sko49, whose flows are of the
   same kind): nothing. tho40: on top of a triangle phase, about one 150-iteration round in the
   climb and two or more at the plateau, at far-miss give-up nodes only; nothing at the near
   misses the phase fathoms in about 100 iterations. Of the order of 10% of the cut-phase
   iterations at far-miss nodes, less on the tree: not material. The gain is a level shift,
   so it needs multipliers that cost no iterations of their own (stepped inside the phase, as
   in the joint rounds here, or inherited).
5. **Not done**: stage 4 (children; the gate did not pass), sko49, the cuts as half-spaces
   inside the CP form's Y-step, a full-protocol (5,500-iteration) wil50 pair, the depth-4 node
   d4_e, the gain of the rotation and fixed-pair families on the large nodes (built; validated on
   small nodes and on one tho40 path node; their gain was not measured). tai50b: violation statistics only (3d').

**What I would do next**, in order: (a) nothing more on wil50/sko with this family; (b) if the
triangle phase goes into tho40's settings (`../tho40_nodes/`), look first at the side
observation of 3c (the rounds form with a target of best + 400 is ahead of the CP form at
equal iterations on the two far-miss nodes measured), which is worth more than these cuts;
(c) only after that, the exchange cuts as an add-on there: row half-spaces in the CP
projector plus inheritance (a child's cut is the parent's with the newly fixed facility's
term moved into c0), expected +135 a node at far-miss nodes; (d) the logical form at the
bottom of the tree, where the LAP elimination test starts to fire (16 of 64 assignments at
nug12 depth 4, 27 of 36 at depth 6) and nodes are cheap anyway.

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
| `runseq.sh`, `runseq2.sh`, `after*.sh`, `jobs_*.txt` | the job runners and lists as they were used (one process at a time; queues were stopped and re-ordered several times, so a list may name jobs that never ran: a log exists only for those that did); each log's first line is its command |
| `s3_summary.py` | the joint-rounds tables (`logs/s3_summary.txt`) |
| `logs/s3_*.log` | stage 3 runs |
| `answer_text.md`, `stage3_text.md`, `mkreadme.py` | the pieces this README was assembled from |

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

Driver `swapdrv.C` on dnncuts' ADMM defaults (PRSM, degree 4, `-T 1`), cold first solve of 300
iterations (wil50, tai50b: 400). All bounds below are certified values; "plain" and "T" arms
run in the same driver, same rounds, so the comparisons are at equal iterations.

### 3a. Violations at converged plain iterates [m]

| node | d | m | cuts | violated | sum of violations (cost x probability) | share on one fixed facility | log |
|---|---|---|---|---|---|---|---|
| tho40 gate `26:10,8:0,39:32,38:7` | 4 | 36 | 5,184 | 192 (3.7%) | 1,413 | 52% on 38@7, 30% on 8@0 | `logs/s3_tho40_d4_ex.log` |
| tho40 d4_d `26:10,8:0,39:32,30:7` | 4 | 36 | 5,184 | 181 (3.5%) | 1,346 | 88% on 8@0 | `logs/s3_tho40_d4d_cp1400.log` |
| the same, after 800 / 1,400 CP iterations of the triangle phase | | | | 314 / 299 | 1,505 / 1,433 | 85% on 8@0 | `logs/s3_tho40_d4d_cp.log`, `..._cp1400.log` |
| tho40 d5_d `26:10,8:0,2:32,6:39,33:24` | 5 | 35 | 6,125 | 130 (2.1%) | 810 | 98% on 8@0 | `logs/s3_tho40_d5d_jTX.log` |
| tho40 path, depth 12 (converged, 2,500 its) | 12 | 28 | 9,408 | 333 | 36 | - | `logs/s3_tho40_path_d12_conv_ex.log` |
| wil50 node 13 `21:28,8:41,28:19` | 3 | 47 | 6,627 | 144 (2.2%) | **15** | 60% on 28@19 | `logs/s3_wil50_n13_jTXm.log` |
| wil50 node 10 `39:27,30:48,25:49` | 3 | 47 | 6,627 | 355 (5.4%) | 68 | 79% on 30@48 | `logs/s3_wil50_n10_stats.log` |
| tai50b `7:40,40:39,42:43` (400 its: far from converged, 21.6 million short) | 3 | 47 | 6,627 | 216 (3.3%) | 756,256 | 89% on 42@43 | `logs/s3_tai50b_d3_stats.log` |

A violated cut reads "given j at b (probability y), exchanging j with the fixed i gains on
average". The largest at tho40 d4_d: fixed 8@0 with 2@8, y = 0.17, conditional gain -1,367;
with 35@8, y = 0.09, -1,235; with 24@8, y = 0.14, -598. In the subtree `26:10,8:0` nearly all
of it sits on facility 8 at corner location 0: the relaxation says 8 would rather trade places
with whatever lands next to it. The triangle phase does not remove these violations (they
grow slightly), so the two families are complementary in principle. On wil50 the violations
are a hundred times smaller: 15 units in all against a gap of 219.

The elimination test (`-elim`: x_jb = 0 when even the best assignment of the rest leaves the
exchange gain negative, one LAP per cut) fires on 0 of 1,296 free assignments at the two depth-4
tho40 nodes and on 1 of 1,156 at depth 6: nothing at these depths (it fires on 16 of 64 at
nug12 depth 4).

### 3b. Exchange cuts alone against the plain solve, tho40 gate node [m]

`logs/s3_tho40_d4_plain.log`, `logs/s3_tho40_d4_ex.log` (target best + 50),
`logs/s3_tho40_gate_ex400.log` (best + 400): 300 iterations, then 8 warm rounds of 150.

| total iterations | 300 | 450 | 600 | 750 | 900 | 1,050 | 1,350 | 1,500 |
|---|---|---|---|---|---|---|---|---|
| plain | 239,727.7 | 240,072.7 | 240,076.9 | 240,080.7 | 240,083.5 | 240,084.4 | 240,085.8 | 240,085.8 |
| exchange cuts (best + 50) | | 240,114.9 | 240,155.5 | 240,186.0 | 240,214.5 | 240,225.0 | 240,235.4 | 240,233.9 |
| exchange cuts (best + 400) | | 240,109.7 | 240,127.7 | 240,209.8 | 239,309.2 | 239,850.5 | 240,226.9 | 240,221.5 |

**+150 over the converged plain value** (+140 with the larger target, which overshoots once).
About 10% of the violation mass at the plain iterate. The triangle CP phase gains +1,327 on
the same node and fathoms it 120 iterations after the first solve (`../tho40_nodes/`).

### 3c. With the triangles, tho40 far-miss nodes [m]

*Joint Polyak rounds* (`-joint`: the library's `dnn_cut_rounds` on the regional triangles,
1,000 new triangles a round, target best + 400; arm TX steps the exchange multipliers between
rounds through the hook, at most 0.5 per multiplier per step). Table by `s3_summary.py`
(`logs/s3_summary.txt`); certified bound after each round, shortfall to 240,514.05 in brackets.

| total its | d4_d: T | d4_d: TX | TX - T | d5_d: T | d5_d: TX | TX - T |
|---|---|---|---|---|---|---|
| 300 | 238,237.8 | 238,237.8 | | 237,994.8 | 237,994.8 | |
| 450 | 238,759.8 | 238,902.3 | +142 | 238,551.4 | 238,631.6 | +80 |
| 600 | 239,134.6 | 239,331.7 | +197 | 238,891.2 | 239,040.1 | +149 |
| 750 | 239,511.7 | 239,636.5 | +125 | 239,251.8 | 239,380.9 | +129 |
| 900 | 239,874.9 | 240,024.1 | +149 | 239,600.0 | 239,694.0 | +94 |
| 1,050 | 240,163.9 | 240,309.2 | +145 | 239,904.6 | 240,039.3 | +135 |
| 1,200 | 240,251.0 | 240,380.2 | +129 | 240,052.6 | 240,123.2 | +71 |
| 1,350 | 240,318.1 (196 short) | **240,535.1 fathomed** | +217 | 239,914.4 | 240,203.0 | +289 |
| 1,500 | 240,288.3 | - | | 240,094.2 (420 short) | 240,204.2 (310 short) | +110 |

Logs `logs/s3_tho40_d4d_jTd400.log`, `..._d4d_jTXm.log`, `..._d5d_jT.log`, `..._d5d_jTXm.log`.
**+70 to +290 at equal iterations, median +135**: the exchange cuts keep, on top of the
triangles, the +150 they give alone. At d4_d that turns a miss into a fathom (the triangles
alone end 196 short after 1,500 iterations; the CP form alone is 214 short after 1,700,
`logs/s3_tho40_d4d_cp1400.log`). The same pair without the step cap: d4_d fathoms at 1,350 too
(`..._d4d_jTXd400.log`, 240,522.7); d5_d leads by +105, +93, +194 and then a Polyak step ten
times the others drops the bound 3,600 (`..._d5d_jTX.log`): the multipliers need the cap.

*Exchange rounds first, then the CP form with the exchange term frozen*
(`logs/s3_tho40_d4d_xcp.log` against `logs/s3_tho40_d4d_cp1400.log`), d4_d:

| CP iteration | 100 | 200 | 400 | 600 | 800 | 950 | 1,400 |
|---|---|---|---|---|---|---|---|
| CP alone (300 + CP) | 238,831.0 | 239,172.5 | 239,583.7 | 239,839.1 | 240,008.6 | 240,105.2 | 240,299.8 |
| 3 exchange rounds, then CP (750 + CP) | 239,019.4 | 239,350.1 | 239,753.5 | 239,996.8 | 240,164.3 | 240,256.8 | |
| difference at equal CP iterations | +188 | +178 | +170 | +158 | +156 | +152 | |

At equal CP iterations +150 to +190; at equal total iterations (1,700) -43, because the 450
iterations spent finding the multipliers are worth about +195 to the CP phase alone at that
point. The gain is a level shift of about 150, not a faster climb: it pays only if the
multipliers come without iterations of their own (inside the phase, as in the joint rounds,
or inherited from the parent).

Side observation, not this agent's topic [m]: at these two far-miss nodes the rounds form with
`-delta 400` is ahead of the CP form at equal iterations (d4_d at 1,350 total: 240,318 against
about 240,160; d5_d at 1,200: 240,053 against 239,780 at 1,100 in `../tho40_nodes/`).

### 3d. wil50 [m]

Node 13 (`21:28,8:41,28:19`, `wil50.dat`, grid 5:10 on side 0, `-u 48816`), a shortened form
of the closed studies' protocol (400 iterations, 5 rounds of 200 at best + 50 uncapped, quiet
round of 300; 1,700 iterations per arm), `logs/s3_wil50_n13_jT.log`, `..._jTXm.log`:

| total its | 400 | 600 | 800 | 1,000 | 1,200 | 1,400 | quiet (+300) |
|---|---|---|---|---|---|---|---|
| T (triangles) | 48,596.8 | 48,656.0 | 48,681.5 | 48,697.6 | 48,709.2 | 48,716.3 | 48,729.6 |
| TX (+ exchange cuts) | 48,596.8 | 48,656.8 | 48,685.1 | 48,692.5 | 48,712.4 | 48,714.0 | 48,731.0 |
| TX - T | | +0.9 | +3.5 | -5.1 | +3.1 | -2.3 | **+1.4** |

**+1.4 at equal iterations, inside the round-to-round scatter of 5**, against a bar of +10.
With the uncapped Polyak step the bound falls 331 in the first round
(`logs/s3_wil50_n13_jTX_overshoot.log`: 15 units of violation cannot carry a step toward
best + 50). Node 10 (`39:27,30:48,25:49`; 3 rounds, no quiet round, `logs/s3_wil50_n10_jT.log`,
`..._jTXm.log`): 48,667.7 against 48,670.4 after round 1 (**+2.7**); in round 2 the library's
own triangle step overshoots in both arms (48,501.2 / 48,469.4) and round 3 has not recovered
(48,606.9 / 48,595.4), so the best certified values are round 1's. Node 10 shows 68 units of
violation, 4.5 times node 13's; at the tenth that tho40 and node 13 realise that is about +7
[e], still under the bar.

### 3d'. tai50b: the violations only [m]

`logs/s3_tai50b_d3_stats.log` (node `7:40,40:39,42:43`, 400 cold iterations, bound 437.22
million against 458.82): 216 cuts violated, 756,256 units in all, which is 3.5% of the node's
gap (tho40 d4_d: 59% of its gap; wil50 node 13: 7%, node 10: 33%). At the tenth realised
elsewhere that is of the order of +75,000, 0.35% of the gap [e]: nothing to build on. Not
measured further (a controlled pair needs about 20 minutes at m = 47; an exchange-rounds run
was cut short for the wil50 node-10 pair, `logs/s3_tai50b_d3_ex_cut_short.log`).

### 3e. Validity [m]

- Brute force, pairwise cuts and fixed-pair cuts: 84 nodes of seven n = 12 instances (depths
  2-4), 1,680 certified values against the enumerated optimum over exchange-optimal completions:
  max(certified - v_ex) = 0.0000 (`logs/validate_small.out`). 5 to 10 of every 12 nodes end
  **above the node's own optimum** v_node, which no inequality valid for all permutations can
  do; the relaxation captures a median of 21% (nug12) to 89% (scr12) of the exact rise.
- All families (pairwise, fixed-pair, rotations) on 112 on-path nodes of the same instances,
  two Polyak targets: 1,568 certified values, none above the optimum
  (`logs/validate_path_small.out`).
- tho40, the path to the best known solution (`../tho40_nodes/scripts/path.py`, cutoff 241,500
  so nothing stops early): exchange cuts alone at depths 4 and 8 end at 235,069 and 238,793
  (`logs/s3_tho40_path_d4_ex.log`, `..._d8_ex.log`); joint with the triangles at depth 12:
  240,255.8 (`..._path_d12_jTX.log`); from a converged plain solve at depth 12 (240,515.40 at
  2,500 iterations) six rounds with all families at a target of best + 2 stay at 240,437 to
  240,457, and the pairwise cuts at best + 400 fall to 128,446 and climb back to 202,793
  (`..._path_d12_conv_exrot_d2.log`, `..._path_d12_conv_ex.log`). **No certified value above
  240,516 in any path run.**
