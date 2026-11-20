# The exchange cut on the mirror nodes (branching_lab agent `mirror`, Nov 19 2026)

**Answer.** The DNN relaxation realises the exchange cut as fully as anything can: on tai35b's twin it
delivers 98 to 100% of the rise that the integer problem itself allows, and the twin (160 of tai35b's 1,017 nodes,
22% of its CPU) is fathomed at depth 2 by its first pooled solve in 60 iterations [m]. On **tai60b's mirror `44:58,43:45` it
cannot work, and no relaxation could**. That node holds a completion 139,349 above the cutoff that is a 2-exchange
local optimum with Delta_(44,43) = -14,070 [m]. Mixing it with the optimum's mirror (+2,786, Delta +2,786) caps
min over the node's convex hull of cost with Delta <= 0, and so any Lagrangian or convex-relaxation bound with
this cut, at **cutoff + 22,060** [m]. Adding every fixed-pair cut does not lift the cap: it stays at +22k to +40k at
every depth of the mirror chain down to 30 [m]. Since that completion is 2-exchange optimal, no per-exchange
cut of any kind (the `swapcuts` family included) can lift it past about +25k [m]. The node is 1,087,650 short, and the
measured law (bound with the cut = integer Lagrangian - the node's gap at mu = 0) gives at most +19k. That is 1.8% of the
shortfall, against the kill rule's 25%. Dropping the node would have removed 175,252 nodes and 1,999 CPU h, 56% of
tai60b's run [m, trace]. The cut does not reach it.
tai50b's 3-cluster twins (depth 3) behave like tai35b: their ceilings are 1.3M to 14.5M above the cutoff.
On a cold solve of one of them, mu = 8 raised the bound by 933k, 82% of the integer rise [m]. With production's
warm state that would fathom it (753k short) [e]. The five such nodes held 1,055 nodes and 29 CPU h, 22% of tai50b's run [m, trace].
In a whole tai35b search the cut removes the twin's subtree, 21% of the run's CPU [m]. Applied blindly, at every node
whose warm start shows a positive Delta, it does more harm than good. The search grew 54% and was stopped, because
where the integer Lagrangian falls with mu the bound collapses. The path to the optimum survived, as it must.

So the test that decides a node is cheap and does not need the DNN. Take the twin's cheapest completion that
violates no cut, mix it with the near-optimal completion the cut removes, and compare the ceiling with the node's
gap at mu = 0. On tai-b, whether the cut can pay is a property of the instance's near-optimal solutions, not of the
relaxation.

Marks: **[m]** measured (script and output named), **[e]** estimated. Everything is in this directory. The solver
hook is in `mirror.patch`, applied to a worktree of cognac 0e6c2ecf (`/home/user/wt-mirror`). Nothing in `latest/`
or in `/home/user/wt-cognac` was touched, and nothing was committed.

## 1. The cut and the hook

For a node with rows i and k fixed (i at a, k at b), let p' be p with the columns of i and k exchanged. Then
Delta(p) = cost(p) - cost(p') = D0 + sum over free rows l of D_l(p(l)), with the dominance README's formula. Any optimal p has
Delta <= 0. For mu >= 0, if a lower bound on min over the node of cost + mu Delta is >= the cutoff, then every
completion below the cutoff has Delta > 0, so the exchange improves it and it is not optimal. The node holds no
optimal solution below the cutoff and can be dropped, under the same invariant as exchange dominance.

**Hook** (`mirror.patch`: `latest/mirror_hook.h`, 13 lines in `bnb.C`, experiment build only):
- `QAP_MIRROR="i,k,mu[,a,b];..."`: right after `reduceQAP`, mu D_l(c) is added to the reduced linear term and
  mu D0 to the shift, both rounded down with `LowerSum`. Every bound (QPB, the first DNN solve, the cut phase with
  its pool) then bounds cost + mu Delta. With a,b given, it applies only at nodes where p(i) = a and p(k) = b.
- `QAP_MIRROR_AUTO="F:CAP"`: at a DNN node whose warm start carries Y (the parent's final Y conditioned on the
  new fix), each fixed pair's Delta is evaluated at the warm start's marginals. The pair with the largest positive
  value gets mu = F × gap / Delta(x_warm), with gap = cutoff - the node's inherited bound, rounded down to 1/64. It is
  skipped if mu > CAP. The term goes in after QPB and before the first solve.
- stderr `mirror:` lines give Delta at the first solve's and the phase's final marginals. `QAP_MIRROR_CHECK=1`
  compares the reduced objective with cost + mu Delta computed from the full problem at every node the hook
  touches: on every completion when m <= 8, otherwise on 500 random ones (200 for the automatic rule).

Build: `make bnb` in `/home/user/wt-mirror/latest` finds the scipy-openblas32 BLAS. The binaries are
`bnbdrv.orig` (unpatched), `bnbdrv.hook` (`QAP_MIRROR`), `bnbdrv.hook2` (+ `QAP_MIRROR_AUTO`) and `bnbdrv.hook3`
(+ the CAP skip). `run.sh NAME BIN ARGS` runs one process with node and DNN traces into `runs/NAME/`.

## 2. Validation [m]

| check | what | result |
|---|---|---|
| hook = cost + mu Delta | `mirror-check` at every node the hook touched in 14 small searches (`validate.sh`: tai10b, nug08, tai12b; 5-6 pairs each, mu 0.25 to 8; incumbent optimum + 1 and about 3% above) | 1,234 nodes, all completions (m <= 8): max difference 1.4e-9, double rounding of values ~4e7; the D0 + sum D_l formula equals Delta exactly |
| certified <= true minimum | every certified value (QPB, first DNN solve, phase) at the 1,210 nodes with m <= 8, against the brute-force min over completions of cost + sum mu Delta (`validate.py` → `validate.out`) | **0 above** |
| the claim | min over completions of cost + sum mu Delta <= min cost over the completions with every Delta <= 0, at the same nodes | 0 failures |
| end of search | all 14 searches end at the optimum (brute force / QAPSOLN) | 14 of 14 |
| mu = 0 | the twin `19:12,18:26` bounded from the same state by `bnbdrv.orig` and by the hook with `QAP_MIRROR=19,18,0,12,26` | node trace, DNN trace (but CPU) and the open-node lines of the output state identical |
| full size | tai35b's twin: certified values for mu = 16, 32, 64, 128 against the best cost + mu Delta of the completions `lagts` finds (simulated annealing; `lagts_twin_mu*.out`, `ceiling_tai35b.out`) | never above: 284.20M vs 285.29M, 285.88M vs 286.96M, 288.66M vs 289.68M, 290.06M vs 294.81M |
| the path | the whole tai35b search with the automatic rule on at every node it could touch (section 5) | all 32 path nodes branched, none fathomed; the optimum found from optimum + 1 (the run was stopped before it closed) |
| the automatic rule | `mirror-check (auto)` on every node it touched (twin run; the 836 nodes of section 5's search) | max difference 0; 6e-8 |

The build reproduces production: root 280,789,962.1796 and node `19:12` 282,252,542.0889, as in the coordinator's
trace (`line_signals/runs/tai35b_base`). The twin differs by +723 (282,378,810 against 282,378,087), since the parent's
state shipped through a different store and build.

## 3. tai35b's twin at depth 2: the bound against mu [m]

The twin `19:12,18:26` swaps the first two fixes of the optimum's path (`19:26,18:12`; rows 18 and 19 are tai35b's
two-location cluster). Here the twin is bounded as the search would bound it: from `19:12`'s state and pool (root →
`19:12` → twin, all with `tai35b.par`, the production tai settings with the bound lab's region file). The cutoff is
283,315,445.05. At mu = 0 the shortfall after the phase is 936,635, and a quarter of it is 234,159.
Scale (`scale.py` → `scale.out`): |D_l(c)| <= 176,930 against the lifted objective's max |entry| of 1.3e7, so up to
mu ~ 70 the term does not change the ADMM's scaling. The optimum's mirror lies in the twin at +118,579 with
Delta = +118,579. A rough scale is then shortfall / Delta at good completions = 8. Over all completions Delta lies in
[-66,094, +343,348] (a rank-one LAP).

`sweep.sh` → `sweep_table.py` → `sweep_tai35b.out`. Values are of cost + mu Delta, and each phase runs at most 2,400 iterations:

| mu | first DNN (its) | after the phase (its) | gain over mu=0 | gain / shortfall | Delta(x) first | Delta(x) phase end | fate |
|---|---|---|---|---|---|---|---|
| 0 | 282,120,575 (110) | 282,378,810 (2,400) | 0 | 0 | +126,606 | +124,406 | branched |
| 0.0625 | 282,127,896 (110) | 282,385,969 (2,400) | +7,159 | 0.008 | +126,510 | +124,335 | branched |
| 0.25 | 282,151,522 (110) | 282,408,836 (2,400) | +30,026 | 0.032 | +126,080 | +124,167 | branched |
| 1 | 282,242,176 (110) | 282,499,060 (2,400) | +120,250 | 0.128 | +124,509 | +123,442 | branched |
| 2 | 282,358,099 (110) | 282,617,143 (2,400) | +238,332 | 0.254 | +122,514 | +122,578 | branched |
| 4 | 282,713,442 (210) | 282,863,277 (2,400) | +484,467 | 0.517 | +117,800 | +120,858 | branched |
| **8** | 283,152,061 (210) | **283,315,477 (2,167)** | +936,667 | 1.000 | +110,056 | +116,187 | **fathomed in the phase** |
| **16** | **283,322,434 (60)** | - | | | +98,625 | | **fathomed by the first solve** |
| 32 | 283,366,798 (40) | - | | | +73,056 | | fathomed by the first solve |
| 64 | 283,426,100 (35) | - | | | +41,709 | | fathomed by the first solve |

At mu = 0 the relaxation's Delta(x) is +126,606, close to the mirror's +118,579, so its solution is
mirror-like. Delta(x) hardly moves as mu grows, and the bound rises by mu × Delta(x). The
automatic rule (F = 2) picks mu = 17.2 from the warm start's Delta(x) = 123,541 and fathoms the twin with its first
solve in 3.7 s, against 217 s and a branching without it (`runs/tw_auto2`).

**Where it peaks** (`sweep_hi.sh` → `sweep_hi_tai35b.out`). The state's incumbent was edited to 290,000,000 so
that nothing stops at the real cutoff, and every solve ran the same 110 + 2,400 iterations. The concave curve is still
rising at mu = 64 (+5.34M above the cutoff). At mu = 128 the phase stopped at the raised cutoff after 160
iterations, so the peak lies beyond 64. Against H(mu), the best value of cost + mu Delta over the completions that
`lagts` finds (an upper bound on the true minimum, `ceiling_tai35b.out`; `gaplaw.py` → `gaplaw_tai35b.out`):

| mu | bound - cutoff | H(mu) - cutoff | H - bound (>= the relaxation's gap) | bound's rise / H's rise |
|---|---|---|---|---|
| 0 | -936,635 | +118,579 | 1,055,214 | |
| 8 | -11,463 | +1,051,427 | 1,062,890 | 0.99 |
| 16 | +886,065 | +1,978,339 | 1,092,274 | 0.98 |
| 32 | +2,562,868 | +3,646,681 | 1,083,813 | 0.99 |
| 64 | +5,344,631 | +6,369,465 | 1,024,834 | 1.00 |

**The law**: the relaxation's gap on cost + mu Delta stays at its mu = 0 value, about 1.05M, within 4% over
mu = 0 to 64. So the bound with the cut is H(mu) minus that gap. The DNN realises the cut completely; the question is
only how high the integer problem's own Lagrangian, max_mu min_p cost + mu Delta, stands above the cutoff. On the
twin it is 16.3M above (at mu ~ 190: the cheapest Delta <= 0 completion found costs 19.5M above the cutoff).

## 4. The ceiling, and tai60b and tai50b

**Ceiling.** Any convex relaxation that contains the node's permutations (the DNN does), with the cut Delta <= 0
added or as min cost + mu Delta, certifies at most min over conv(node) ∩ {Delta <= 0} of cost. By LP duality that is
max_mu min_p cost + mu Delta, which is at most max_mu min over a set S of completions. `ceiling.py` builds S from
the optimum's mirror and `lagts` runs over a mu ladder. `ceiling_multi.py` does the same for all fixed pairs at once,
by column generation: it solves the LP over the multipliers, runs `lagts` at the LP's multipliers and repeats until it
finds nothing below the LP value. Each value printed is a valid upper bound on what the cuts can certify.
Fathoming needs ceiling - cutoff > the relaxation's gap. Under the law above, the best bound is the ceiling minus the
mu = 0 gap.

| node | its shortfall (production trace) | ceiling - cutoff, one pair | all fixed pairs | predicted best bound - cutoff | log |
|---|---|---|---|---|---|
| tai35b twin `19:12,18:26` (d2) | 936,635 (this build) | +16.3M | (one pair) | +15.2M: fathoms (it does) | `ceiling_tai35b.out` |
| **tai60b mirror `44:58,43:45` (d2)** | **1,087,650** | **+22,060** (mu 6.9) | (one pair) | **-1.07M** | `ceiling_tai60b.out` |
| tai60b mirror chain, d3 / d5 / d8 / d11 / d15 / d20 | 758,810 / 545,752 / 390,011 / 303,339 / 244,966 / 199,188 | +22k to +25k | +25,357 / +22,084 / +24,017 / +22,084 / +24,017 / +22,084 | far below | `ceiling_tai60b_chain.out`, `ceiling_multi_tai60b_chain.out`, `ceiling_multi_tai60b_d11.out` |
| tai60b mirror chain, d25 / d30 (994 / 736 nodes left) | 36,531 / 30,900 | | +39,998 / +39,998 | below | `ceiling_multi_tai60b_chain.out` |
| tai50b `31:0,33:32,32:7` (d3, pair 31,32) | 753,372 | +1.37M (mu 8.6) | +1.34M | +0.47M: fathoms [e] | `ceiling_tai50b_31-32.out`, `ceiling_multi_tai50b.out` |
| tai50b `31:7,33:0,32:32`, `31:32,33:7,32:0`, `31:0,33:7,32:32`, `31:32,33:0,32:7` (d3) | 729,763 / 652,488 / 579,997 / 578,289 | | +14.0M / +8.1M / +14.5M / +6.5M | fathom [e] | `ceiling_multi_tai50b.out` |

**Why tai60b's mirror resists** (`obstacle.py` → `obstacle_tai60b_d2.out`, `obstacle_tai60b_d11.out`;
`mirror_exchanges.out`):
- The optimum's mirror (+2,786 above the cutoff) has exactly one improving exchange, (43,44) by 2,786.
- The same node holds completions at **+139,349** (Delta_(44,43) = -14,070) and +217,640 (-218) that are
  **2-exchange local optima**: no exchange of any two rows improves them.
- A mixture of 83% mirror and 17% of the first satisfies every per-exchange cut, whether fixed-fixed, fixed-free
  (the `swapcuts` family, x_jb × gain >= 0 in Y) or free-free. Each of these says one exchange's gain is >= 0. The
  local optimum satisfies all of them, the mirror all but (43,44), and on (43,44) the mixture has Delta = 0.
  Its cost is about cutoff + 25k.
- So no per-exchange cuts, linear in x or in the lifted Y, can lift any relaxation of the mirror chain more than
  about 25k above the cutoff. The DNN's gap there is 1.09M at depth 2, and still 199k at depth 20.

In tai35b the twin has no such completion: the cheapest Delta <= 0 completion found is 19.5M above the cutoff.
tai60b's (43,44) is a near tie. The two facilities are almost interchangeable at good solutions (the mirror costs
2,786 more), and good solutions exist on both sides of the tie.

**The DNN solve on tai60b's mirror, once** (cold from `-F 44:58,43:45`, `tai60b_first.par`: production settings,
first solve only; `runs/t60_first_mu0`, `runs/t60_first_mu7`):

| | mu = 0 | mu = 7 (the ceiling's maximiser) |
|---|---|---|
| first DNN solve | 584,195,461 (470 its, hopeless, 571 s) | 584,680,580 (490 its, hopeless, 599 s) |
| Delta(x) | +6,182 | -38,959 |

The cold solve gives up 24M short (production's warm pooled solve: 606.65M, then 607.13M after the phase). Its
Delta(x) of +6.2k is mirror-like (+2,786), as on tai35b, and at mu = 7 the relaxation moves straight to
Delta < 0. Its +485k is the difference between two unconverged iterates 24M from the cutoff, taken at different
iteration counts. It is not a converged gain: at convergence the bound at mu = 7 is at most H(7) = cutoff + 22k. The
tai35b rule (F = 2 from the relaxation's Delta) would give mu = 2 × 1.09M / 6.2k ≈ 350 at production's gap, where
H is more than 178M below the cutoff (H(256) = cutoff - 178.8M). **Kill rule: the best possible gain, about 19k, is under a quarter of the 1.09M
shortfall. Stopped.** If the node had been dropped at depth 2, 175,252 nodes (52.9%) and 1,999 CPU h (55.6%) of
tai60b's tree would have gone (`subtrees.py` → `subtrees_tai60b.out`).

**tai50b, one node, cold** (`-F 31:0,33:32,32:7`, `tai50b_cold.par`: production settings with the first solve's
give-up at 1,000 iterations; `runs/t50_mu0`, `runs/t50_mu8`):

| | mu = 0 | mu = 8 |
|---|---|---|
| first solve (1,000 its) | 434,396,980 | 434,351,786 |
| after the phase (2,400 its) | 456,211,602 (2.61M short) | 457,144,339 (1.68M short) |
| Delta(x) at the phase's end | +201,740 | -72,187 |
| gain over mu = 0 | | **+932,737**, 82% of H(8) - H(0) = 1,141,161 |

The law holds to 82% on a second instance. Production's warm node is 753,372 short after its phase, and the cold gain
of 933k would close it [e]. The five 3-cluster twins at depth 3, all with ceilings far above, held 1,055 of
tai50b's 8,225 nodes and 29.3 of 132.7 CPU h (22%) (`subtrees_tai50b.out`) [m, trace]. Their branched depth-2 parents
`31:0,33:32` and `31:7,33:0` stay.

## 5. A whole tai35b search with the automatic rule [m]

`search35.sh s35_auto` (`bnbdrv.hook3`, `QAP_MIRROR_AUTO=2:64`) ran the whole search from the root's state
(`runs/root`, unpatched build) with the cut at every DNN node where some fixed pair has a positive Delta at the warm
start. `search_table.py` → `search_s35_auto.out` compares it with the production trace (the same settings via the
coordinator, `line_signals/runs/tai35b_base`):

| subtree | production: nodes / CPU s | with the rule: nodes / CPU s |
|---|---|---|
| `19:12` (the twin's parent) | 194 / 1,766 | **35 / 329** |
| `19:8` | 35 / 295 | 35 / 299 |
| `19:26` (the optimum's) | 755 / 4,095 | 1,476 / 9,360, **stopped with 104 nodes open** |
| all, with the root | 1,017 / 6,931 | 1,571 / 10,699 (stopped) |

- **The twin.** The rule picked mu = 17.2 there, and the twin died at depth 2 with its first solve. Its 159 descendants
  and 1,437 CPU-s went: **21% of the run's CPU, 16% of its nodes**. With the cut confined to the twin
  (`QAP_MIRROR=19,18,16,12,26`), nothing else in the search changes: the hook leaves every other node alone, and
  mu = 0 is trace-identical. So this saving is what the targeted cut gives [m for the subtree, e for the whole run].
- **Everywhere else, the blanket rule hurts.** It touched 836 nodes. Of those that also exist in production, 51 fathomed in
  both and 1 (the twin) fathomed only with the rule. Nine production fathoms became branchings (6) or phase fathoms (3).
  Under the optimum's path some bounds collapsed: `19:26,18:12,26:5,23:6,31:34` fell to 225,754,021 and was branched.
  A positive Delta at the warm start does not mean the integer Lagrangian rises with mu (section 4: on tai50b's node
  H falls steeply past mu ~ 8.6). The bound changes also changed the branching lines below depth 4. The search was
  stopped at 10,699 CPU-s, 54% above production's whole run, with 104 nodes open. A rule needs the ceiling test.
- **Path check: passed.** All 32 nodes of the optimum's path were bounded with the rule active and branched, none
  fathomed. The leaf was enumerated and the search found 283,315,445, the optimum, from an incumbent one above.
  (Not "no open node", since the run was stopped.) The `mirror-check (auto)` lines at the 836 touched nodes
  (200 random completions each) show a largest difference of 6e-8.

## What was not tested

- tai60b's mirror under production conditions (warm from `44:58`'s state and pool). The root, `44:58` and the node
  cost about 2.2 CPU h per arm here. The ceiling makes the outcome certain to within the law's accuracy: the bound
  cannot pass cutoff + 22k at any mu, so fathoming needs the relaxation's gap to fall from 1.09M to under 22k.
- tai50b under production conditions, and a tai50b search with the rule. The 82% is a cold solve's.
- Other exchange families as a lift for the tai60b chain: 3-cycles, or k-exchanges with k >= 3. The local optima
  that cap the ceiling are only 2-exchange optimal. A 3-exchange might improve them, and then their mixture with the
  mirror would no longer satisfy all cuts. Not measured; `ceiling_multi.py` would take a 3-cycle family with small changes.
- The rule's mu: F = 2 suits tai35b. On tai50b's node the ceiling falls steeply past mu ~ 8.6 (completions with
  Delta < 0 at +1.6M), so a rule must stay below the ceiling's peak. A cheap SA ceiling at the node gives both the
  go/no-go and the mu.

## Cheapest next step

A production rule needs neither the DNN nor a sweep to decide where to apply the cut. When a node's first solve misses,
run `ceiling.py`'s computation for its best pair: its Delta at the warm marginals, then a short SA on cost + mu Delta,
which takes seconds in C at n = 50. Apply mu (about the ceiling's argmax) only when the ceiling minus the node's
mu = 0 gap is above the cutoff. On tai50b, test that rule on the five 3-cluster twins from production's warm state
(the root state is the expensive part, 4,286 CPU-s). Expect about 22% of the run's CPU [e]. On tai60b it would never
fire on the mirror chain.

## Files

| | |
|---|---|
| `mirror.patch` | the hook (worktree `/home/user/wt-mirror`, cognac 0e6c2ecf) |
| `run.sh`, `tai35b.par`, `tai60b.par`, `tai60b_first.par`, `tai50b_cold.par`, `small.par` | one solver process; the settings (production recipes with absolute region paths; `_first`: no phase; `_cold`: first-solve give-up at 1,000) |
| `validate.sh`, `validate.py` → `validate.out` | section 2 |
| `runs/root`, `runs/d1` | the root and `19:12`, bounded with the unpatched build; `runs/twin.state` = `splitstate.py runs/d1/d1.state runs/twin.state -pick 18:26` (bound lab's `tai_solves/splitstate.py`) |
| `sweep.sh`, `sweep_hi.sh`, `sweep_table.py` → `sweep_tai35b.out`, `sweep_hi_tai35b.out`; `gaplaw.py` → `gaplaw_tai35b.out`; `runs/tw_*`, `runs/hi_*`, `runs/tw_auto2` | section 3 |
| `lagts.C` (build: `g++ -O3 -march=native -o lagts lagts.C`) | simulated annealing (and a tabu mode) on cost + sum mu Delta at a node; prints every restart's best with its cost and Deltas |
| `ceiling.py`, `ceiling_multi.py`, `obstacle.py`, `chain.py`, `subtrees.py`, `scale.py`, `basics.py` and their `.out` | section 4 |
| `runs/t60_*`, `runs/t50_*` | the tai60b and tai50b solves |
| `search35.sh`, `search_table.py` → `search_s35_auto.out`, `runs/s35_auto` | section 5 |

State files and parent-state stores (88 MB) are not in this directory: they are in the session scratchpad
(`/tmp/claude-0/-home-user/92b66c47-7a9f-5ca5-8a65-57e14ac5b075/scratchpad/mirror_states/`, may be cleared).
`runs/root` and `runs/d1` rebuild them: `run.sh root orig -p tai35b.par -u 283315446 -d 0 -o runs/root/root.state ...`,
then `-r` on the `19:12` node with `-d 1`. The commands are in each run's `cmd.txt`.
