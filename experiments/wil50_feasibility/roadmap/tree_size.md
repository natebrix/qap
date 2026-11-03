# wil50 roadmap: tree size (Nov 3 2026)

One of three memos (bound strength, per-node speed, tree size) on proving wil50
(best known 48,816, cutoff 48,814.05) optimal in about a week on 1,000 Xeon-class
cores, i.e. about 6e8 CPU-s. This one covers everything between the bound and the
node cost: symmetry, branching, where the DNN bound and the cut phase run, the
cutoff, and the tree model. Nothing was run except small Python arithmetic.
**[m: file]** = measured, with its source; **[e]** = estimate (the small models are
specified inline).

## Bottom line

1. **The symmetry file is not a bonus; the estimates already assume it.** NOTES.md
   section 4 counts 3e4 depth-3 nodes = 50·49·48/4, but every wil50 search so far
   ran without one ("couldn't find [./SYM/wil50.sym]" [m: stronger_bounds/search/wil50n10ext_*.log]).
   Without the file the run is 3.98x the estimate. Building it is half a day.
2. **The tree estimate is low, and the tail sets it.** NOTES.md's own tail assumption
   (one depth-3 node in nine or ten reaching depth 6, about 1e5 DNN nodes each) gives
   2e8 to 3.7e8 DNN nodes, not "nearer 1e8". The 2 of 16 sampled depth-3 nodes that are
   329 to 405 short account for 85 to 91% of every estimate I can build. My range is
   **3e8 to 1.8e9 DNN nodes** [e], which is 1.0e10 to 6e10 CPU-s with the cut phase,
   **120 to 700 days on 1,000 cores** [e], against the brief's 26 to 83.
3. **The largest tree lever is the branching line at the DNN depths.** "130 units per
   level" is what a random line gains. On the DNN-marginal line, 31 of 47 children of
   a 241-short parent fathom outright, against 2 to 9 of 47 on random lines
   [m: tuning/lines/wil50_lines_m2pro.log]. Today the dominant class branches by rule
   6, 7 or 2, never by 8 or 9.
4. The cut-phase knobs (zone, stall, depths) are within about 10% of their optimum.
   The one rule that must not change is running the phase at every depth from 4 down.
5. Before 1,000 cores are committed, a sampling plan costing 1e7 to 2e7 CPU-s (2 to 3%
   of the week) can pin the total to within 1.5x. The sampled subtrees are pieces of the
   proof, so none of that work is wasted.

## Ranked table

| # | idea | factor on the wil50 total | effort | risk | evidence | cheapest test; result that kills it |
|---|---|---|---|---|---|---|
| 1 | Symmetry file for the transposed instance | 3.98x against today's runs [m: Burnside, below]; 1x against the estimates, which assume it | 0.5 day | low | 15 / 635 / 29,580 orbits at depths 1-3 against 50 / 2,450 / 117,600; `latest/SYM/nug30.sym` has the same structure; sko49 took an hour | `dnntest` on 3 depth-3 nodes and their transposes, plus a bnbdrv run to depth 2. Killed if the certified lb differs by more than 1e-6, or the counts are not 15 and 635 |
| 2 | Rule 9 (or 8) at every DNN depth, with the DNN from depth 1 so that depths 1-2 branch on marginals too | 2-5x [e] | 1 day + A/B | medium: rule 8 cost +37% on nug28, rule 9 +8% [m: tuning/README.md] | xbar line best at 3 of 3 wil50 parents; eccentric 2-4.7x the best line, random 4-6x (lines cost metric) [m: lines log] | A/B on the node-10 subtree (cut arm, 6 h on 4 cores), then node 3 (242 short). Killed if rule 9 costs ≥ 0.9x the current rules' CPU on both |
| 3 | Sampling plan (section 3) | none; narrows a 15x range to 1.5x | 2 days to set up, 1e7-2e7 CPU-s | low | - | - |
| 4 | Defer the cut phase behind an uncut give-up confirmation (`QAP_DNN_GU_ZONE=0.003`, longer deferral), so slow uncut fathoms do not pay for a phase | 1.05-1.15x [e] | 0 (env var) | low | 3 of 7 uncut depth-4 fathoms needed 510-630 iterations; the phase fathomed them in 1,600 iterations, 868-918 s [m: stronger_bounds/depth4]; run_ab3 ran without the flag | Replay on the 16 depth-4 children with `-Z`. Killed if no phase-fathomed node would have fathomed uncut |
| 5 | Primal gate in the phase (stop after round 3 if the best primal is more than 15 below the cutoff), zone 0.0035 | 1.03-1.06x [e: model] | 0.5 day | low | Round 3 primal: 4 successes within 1.4 of the cutoff, 3 failures 34-85 below [m: laptop/dual/*.log] | Replay the rule on every logged phase. Killed by one false stop |
| 6 | Children of a failed phase start from its best multipliers and cut pool | 1.1-1.3x [e] | 2-3 days | medium | 316 of the 365 nodes in run_ab3's cut arm are children of failed phases [m/e, section 4] | One failed node's 45 children, inherited vs not |
| 7 | Longest-first child order, asynchronous checkpoints | 1.0x CPU; 1.1-1.3x wall [e] | 1 day | low | A parallel.py checkpoint idles every worker until the slowest node ends [m: stronger_bounds/README.md] | - |
| 8 | Short-DNN strong branching at depth ≤ 3 | ≤ 1.2x over #2 [e] | 2 days | medium | xbar already finds the best line | - |
| - | Dead: region branching (a block node is as strong as its weakest cell [m: region_branching]); certificate eliminations (0 of 2,209 [m: tuning/README.md section 2]); flow automorphisms (none, below); tighter incumbent; proving 48,800 first | 1x | | | | |

## 1. Symmetry (Q1)

The first matrix of `qaplib/wil50.dat` is exactly the 5x10 Manhattan distance
[m: checked entry by entry]. Its group is {id, r→4−r, c→9−c, 180°}. The row
reflection fixes row 2 (10 cells); nothing else fixes a cell. Burnside gives
(P(50,d) + P(10,d))/4 orbits, against the plain counts:

| depth | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| plain | 50 | 2,450 | 117,600 | 5.53e6 | 2.54e8 |
| orbits | 15 | 635 | 29,580 | 1.38e6 | 6.36e7 |
| ratio | 3.33 | 3.86 | 3.98 | 4.00 | 4.00 |

The solver's J1/J2/J3 mechanism (`latest/branch.C`, `determine_J3`) achieves exactly
these counts with one J2/J3 pair, as `nug30.sym` does for 5x6. The file is
J1 = rows 0-2 × columns 0-4 (15 cells), J2 = row 2 (cells 21-30) and J3 = rows 0-2
(cells 1-30). At depth 2 that gives 10·49 + 5·29 = 635 nodes, and at depth 3
590·48 + 45·28 = 29,580. Like sko49, the file applies to the transposed instance:
the solver's symmetry acts on locations (the second matrix), and wil50's grid is the
first. The transpose has three consequences. `-F` specs swap i:j to j:i; the cut
phase needs `QAP_DNN_CUT_SIDE=1`; and at nodes whose assigned cells all lie in row 2
the branching is restricted to rows, which `select_ind_dnn_xbar(rows_only)` already
handles (45 depth-2 and 360 depth-3 nodes). Swapping A and B permutes the DNN variables,
so the bound and the warm starts should agree to rounding; the test in the table
checks this. The regional cuts are valid at every node, so symmetry does not touch
them.

The flow matrix has no automorphisms. All 50 rows have distinct value multisets, and
colour refinement gives 50 classes [m: Python on wil50.dat]. On the grid the same
refinement gives 15 classes, matching the 15 orbits. The order-4 group is all there
is. Objective parity: values are even (already used); random permutations hit every
residue 0, 2, 4, 6 mod 8 [m: 300 samples], so there is nothing more to gain.

## 2. Branching at the DNN depths (Q3)

`pick_param` (bnb.C) chooses the parameter block from the inherited bound, and that
is the parent's DNN value. In the `-F` runs a child of node 10 sits at a relative gap
of about 0.11 (rule 6, eccentric), and children of depth-4 misses sit at 0.03-0.08
(rule 6 or 2) [e: my reading of the code]. A full run lands in rules 7, 6 and 2.
Rule 2 uses QPB's U, which fathoms nothing at these gaps, so it picks an essentially
arbitrary line. On wil50's clear misses, cost = Σ 47^(gap/130) [m: lines log]:

| parent (short) | best line | lookahead | eccentric | random mean | children fathomed: best vs random |
|---|---|---|---|---|---|
| 1 (405) | 8.5e4 (xbar) | 1.58x | 4.68x | 5.9x | 0 vs 0 (16 vs 0-1 within 130) |
| 3 (241) | 5,143 (xbar) | 1.11x | 1.96x | 3.7x | 31 vs 2-9 of 47 |
| 10 (158) | 486 (xbar = lookahead) | 1.00x | 2.07x | 4.2x | 33-36 vs 4-9 of 47 |

The tree model's "130 per level" is a random-line figure (NOTES.md section 3 drew a
random facility). In the exponential model, subtree ∝ 47^(s/130). "Children 10%
closer to the cutoff" is then worth 47^(0.1·s/130): 1.6x for a 158-short node, 2.0x
at 242, and **2.7-3.3x for the 330-405 tail** [e], which carries the total. The best
line beats lookahead by 1.0-1.6x at depth 3. Below depth 3 the current rules are
2-6x worse than the best line, and that factor applies again at each level of a tail
subtree. Hence the 2-5x for item 2, discounted for nug28's warning. Fujii et al.'s
three rules spread their wil50 estimate 1.9x (2.05e7 to 3.93e7) [m: sko49 NOTES
section 7]: same order. Child order does not change the total, since the cutoff is
fixed and warm starts come from the parent, not siblings. It only matters for load
balance: hardest children first. A short-DNN strong branching (47 children × 4 lines
× about 30 s, about 6e3 s a node) is affordable at depth ≤ 3. It can only beat rule 9
where xbar's two lines are both wrong, which none of the three parents showed.

## 3. The tree model and the sampling plan (Q5)

**NOTES vs model.py: the model has no bug; NOTES's tail arithmetic is short.** A
depth-3 node fathoming at depth 6 costs 47 + 47·46 + 47·46·45 = 99,499 DNN nodes.
One such node in nine clear misses is about 6x the median class, not 2-4x. That gives
2.2e8, and one in ten depth-3 nodes gives 3e8. model.py rounds two nodes (405 and the
extrapolated 329, 2.53 levels) to depth 6 and gets 3.9e8 [m: re-run]. A step function
at s/g ≈ 2.5 swings the total 2x on one sample, so neither number is robust.

**Measured anchors in the search** (run_ab3, Xeon 2.1 GHz, `nug_deep.par`, H = 300,
EXT = 20): near misses 48 nodes [m: laptop/search]. Node 10 (158 short): 2,193 nodes,
63 CPU-s each. Its node and phase counts reconstruct, with lines.py's 33/47 depth-4
fathoms on the root's lookahead line (the logs print no depths), as about 14 depth-4,
33 depth-5 and 1 depth-6 misses (1 + 47 + 14·46 + 33·45 ≈ 2,193) [e]. Building tail
subtrees from the lines data with these anchors (children within 130 at 300
iterations cost about 153 nodes, as node 10's depth-4 misses did; children 130-260
short cost about 2,193) gives per depth-3 node: 5e3-2e4 nodes for the 233-244 class
and 7e4-4e5 for the 329-405 class [e]. That makes **1.0e4-6e4 DNN nodes per
depth-3 orbit, 3e8-1.8e9 in all**, 85-91% of it from the two tail nodes. An i.i.d.
gain model fitted to node 10 (G ~ U[72,200]) gives 1.6e10. It overcounts, because it
ignores the line effect above, but it shows how far the tail can swing.

**The 16-node sample is the random-line tree.** The production tree's depth-3
frontier comes from rule 4 or 9 lines at depths 0-2 and can be much better (see the
fathomed counts in section 2). That is the main upside risk to my range. The plan
therefore samples the frontier the solver actually builds:

| stage | what | CPU-s [e] | output |
|---|---|---|---|
| A | Run depths 0-2 with the production settings (sym file, rule 9, DNN from depth 1); save the frontier | 651 orbits × 250 = 1.6e5 | the real depth-3 frontier |
| B | 400 random frontier nodes to the give-up point; converge the 80 farthest | 6.6e4 + 1e5 | class fractions; tail fraction ±1.6 points |
| C | End to end, production settings: 3 near, 6 at 150-210, 6 at 210-300 short | 2e6-6e6 | T(s) for the bulk |
| D | 10 tail (≥ 300) nodes: 30 Knuth dives each (bound every child on the path, follow one random miss; about 1.2e4 s per dive), plus 3 run end to end at cutoffs c−130 and c−65 (at c only if affordable; the "cutoff ladder" of section 5) | 5e6-1e7 | the tail mean within about ±35% (1σ) |

That is 1e7-2e7 CPU-s in all: 3-6 hours on 1,000 cores. On the 4-core box, stages A
and B take about a day and stage C one to three weeks. The decision rule: commit only if
(frontier fraction × T) summed over the strata, times 63 s / 1.81, comes under the
budget with every stratum's 1σ interval inside 1.5x.

## 4. Where the DNN and the cut phase run (Q2)

- **D1/D2.** Start the DNN at depth 1. The 651 top orbits cost about 1.6e5 CPU-s
  (< 0.1%) and give rule 9 its marginals at depths 1-2. On nug30, DNN from depth 1
  cut depth-3 misses from 59 to 45 [m: tuning/README.md refine]. With
  `QAP_DNN_EXT=20`, D2 does not matter as long as it is ≥ 3: every child of a miss
  gets the DNN bound. Keep EXT, and keep today's state.C fix.
- **The phase must run at every depth ≥ 4, not only at depth 4.** In the run_ab3
  reconstruction, 14 phases ran at depth 4 (10 fathomed) and 25 at depth 5 (21
  fathomed): 1 + 47 + 4·46 + 3·44 = 364 nodes against 365 [m/e]. With depth-4
  phases only, the 25 depth-5 misses would branch: about 1,360 nodes and about 1.2x
  instead of 1.81x [e]. The economics do not change with depth (a phase at about
  900 s [m: (76,755 − 41,796)/39] against 45 children at about 61 s), so depth 5 and
  below pay as well as depth 4.
- **Zone and stall.** I fitted a stochastic model (gain U[72,200], give-up lag 55, a
  phase that fathoms with probability Φ((125 − u)/20) at give-up shortfall u, success
  1,100 s, failure 1,300 s) to both arms of run_ab3: model 2,103 nodes / 133k s and
  363 / 73k s, against 2,193 / 139k and 365 / 77k measured. Its total varies by only
  2% across zones 0.0025-0.0035, and is 10% worse at 0.002 [e]. Round-1 gains do not
  predict failure: failures gain +31 to +40 and successes +28 to +34 [m: laptop/dual
  logs]. The round-3 primal does predict it (item 5). With that gate the best zone
  moves to 0.0035, worth 4% [e]. Keep the stall rule at 3 rounds.
- **The cut-arm per-node DNN cost was 114 s against the base arm's 63 s** [m]: the
  children of failed phases are the hard ones (an estimated 14% of them miss). That is
  why item 6 is worth testing.

## 5. Upper bound and cutoff (Q4)

48,816 is the long-standing best known. A better incumbent search buys nothing if 48,816 is
optimal, and if it is not, the search finds the better solution anyway. Proving "no
solution below 48,800" is only 16 units weaker, a tree 47^(16/130) ≈ 1.6x smaller
[e], and tells us nothing. The useful version is the **cutoff ladder**. The tree at
cutoff c − X is exactly the "+X bound" row of the model. So running sampled subtrees
(or the whole tree) at c − 130 (about 30-47x cheaper [e]), then c − 65, measures the
growth per unit of cutoff directly, and a full-tree pass at c − 130 is a rehearsal at
2-3% of the cost. It would also certify a lower bound of about 48,686.

## 6. The 1,000-core run

| configuration | DNN nodes [e] | CPU-s [e] | 1,000 cores [e] |
|---|---|---|---|
| Today as configured (no sym file; rules 4/7/6/2; phase zone 0.003; EXT) | 1.2e9-7e9 | 4e10-2.5e11 | 1.3-8 years |
| Today + the symmetry file | 3e8-1.8e9 | 1.0e10-6e10 (63 s/node [m], /1.81 [m]) | 120-700 days, center about 290 |
| + items 2, 4, 5 (2.2-6x) | 5e7-8e8 | 1.7e9-2.7e10 | 20-330 days, center about 80 |

On top of the tree items, a week needs about 12x more at the center (3-50x across the
range) from bound strength and speed together. On a tail-dominated tree each +65
units at depths 3-6 is worth about 47^0.5 ≈ 7x.

## 7. Interactions with the other two angles

- **Bound strength.** The tail decides the total, so gains at depths 4-7, inside the
  tail subtrees, count as much as gains at depth 3. The optimal cut zone tracks the
  phase's gain at the give-up point (zone ≈ gain − 1σ). Every +50 of phase gain
  widens it by about 0.001.
- **Speed.** Cheaper DNN iterations scale the phase and the children alike, so the
  zone and depth choices hold. Rule 9's QPB lookahead is negligible next to a DNN
  node. Changes to the give-up rule interact with item 4 (slow fathoms should finish
  uncut before any phase).
- **Parallel.** About 29,580 frontier nodes is ample parallelism. The tail subtrees
  (1e5-4e5 nodes each) need work stealing and longest-first order, and the
  synchronous checkpoint has to go at 1,000 cores.
