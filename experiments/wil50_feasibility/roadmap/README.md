# wil50 in a week on 1,000 cores: the roadmap (Nov 3, 2026)

Three research memos written on Nov 2-3 from the measurements of Oct 28 to
Nov 2: `bound_strength.md`, `node_speed.md`, `tree_size.md`. This page puts
their factors together. Every number is an estimate unless it names a log.

## 1. The baseline is larger than we had been saying

The brief gave 2.2e9-7e9 core-seconds after the cut phase (26-83 days on
1,000 cores). `tree_size.md` re-derives the tree from the one end-to-end
subtree (node 10: 2,193 DNN nodes at 63 Xeon-s each, `stronger_bounds/
search/wil50n10ext_base.log`) and the `tuning/lines` child data, and finds
two things the earlier estimate missed:

- **The symmetry file was never used.** NOTES.md section 4 counts the
  depth-3 nodes modulo the 5x10 grid's order-4 symmetry (29,580 orbits), but
  every wil50 run printed "couldn't find [./SYM/wil50.sym]". As configured
  today the search is 3.98x the estimate. The file is half a day's work on
  the transposed instance (the solver applies symmetry to the second matrix;
  the sko49 session did the same for the 7x7 grid) and is a prerequisite,
  not a bonus.
- **The depth-6 tail sets the total.** Two of the 16 sampled depth-3 nodes,
  the ones 329-405 short, carry 85-91% of the work; the estimate's step
  function at the level boundaries swings the total 2x on one sample. With
  the symmetry file: 3e8-1.8e9 DNN nodes, 1.0e10-6e10 CPU-s with the cut
  phase, **120-700 days on 1,000 Xeon cores, center about 290**.

So a week on 1,000 cores (6.05e8 core-seconds) needs about **41x** at the
center and **17-100x** across the range, after the symmetry file and the cut
phase (corrected Nov 3: this line said 12x and 3-45x, which do not follow
from the range above).

## 2. What the three angles offer (factors on the total)

| angle | item | factor | effort | first test |
|---|---|---|---|---|
| tree | symmetry file for wil50 (transposed instance) | 3.98x against today's runs, already in the estimate | **done** (`sym/README.md`: `latest/SYM/wil50t.sym`, `qaplib/wil50t.dat`) | passed: 15 and 635 children at depths 1 and 2 (651 nodes against 2,501 plain); node 10's bound agrees to 0.02 after a cut round with `-side 1` |
| tree | branch on the DNN marginal line (rule 9) at every DNN depth, DNN from depth 1 | **measured on the node-10 subtree: 1.41x less CPU, 1.13x fewer nodes** (rule 8 gives the same tree at 1.23x; `stronger_bounds/README.md` step 12); the 2-5x estimate was too high because rule 2 already finds decent lines | done | the node-3 subtree would confirm; DNN from depth 1 untested |
| bound | children inherit the parent's cut pool and multipliers (and a phase at clear misses to equip their children) | 1.5-4x on the dominant class (bound memo); 1.1-1.3x (tree memo). **Gate passed Nov 3** (Astra, `research/pool-inheritance`): one frozen solve on the inherited pool gains a median +54 over the child's own uncut value, d1 and s5 fathom from it; continuing with Polyak rounds adds only +10 to +21 before the quiet round, at 3.3x the iterations, so the frozen solve is the policy. **Measured in the search Nov 4: 1.37x end to end on top of the stacked arm** (`pool_search/README.md`): the same 322-node tree in 23,905 CPU-s against 32,740; all 228 children that inherited a pool fathomed on their first solve (101 iterations on average), cut phases fell from about 14,400 to 7,000 CPU-s | done (`QAP_DNN_POOL=1`, the pool travels in state files) | production: `QAP_DNN_POOL=1` |
| bound | children start from the parent's final Y as well as its Z (`warm_y/README.md`) | **measured Nov 3 on 16 depth-4 children: 1.17x per child** (1.29x on fathoms, 24% fewer iterations; +67 at iteration 100 on misses; no outcome changes) | done; in the search as `QAP_DNN_WARM_Y=2` | in the stacked node-10 A/B (`run_ab6.sh`) |
| speed | AMX bf16 for the filter and Rayleigh-Ritz products, plus CPU trims | **measured Nov 3 on the 32-node sample: 1.71x on fathoms, 1.75x on misses, same 18 fathoms, same iterations** (bundle D: AMX kernel for the five W*Y products, degree 3 with tight lambda, one-pass CholQR, refresh every 40; `amx/README.md`). Closed as not worth it: guard K/40, fp32 R (certificate margin), the small products on AMX, the row fusion; neutral: the symmetrisation skip. Left: the kernel's last third (ceiling 2 TF against 1 TF reached) | done | production: `DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40` |
| speed | cut-phase trims: shorter inner solves with Y/Z restoration (Astra's Y12s: `-warmY`, 12 rounds x 200 iterations, quiet cap 1,000) | **measured on the phase: 1.56x** (64% of baseline CPU, bound and fathom gates kept on all six nodes; `same_node_continuation/budget_study/RESULTS.md`); 1.17-1.20x on the total at the 46% phase share | **not pursued** (Nathan, Nov 3: too small); option exists | none |
| speed | coordinator instead of the checkpoint barrier and node-end steals | 1.15-1.35x of the fleet [e]; **built Nov 4** (`coord/`, plan.md phase 1, local and AWS backends, exact totals under worker kills, coordinator restarts and drains); on the 4-slot node-10 subtree the same tree with **1.09x fewer DNN iterations** (warm starts kept: 23 units, no barrier) and the same 81% busy (the subtree's tail) (`coord_search/README.md`) | phase 2 of plan.md (Terraform, Spot handler) to run it on a fleet | the fleet factor needs plan.md phase 2 (100 Spot cores) |
| tree | confirm give-ups uncut before the phase (QAP_DNN_GU_ZONE was off in run_ab3), primal gate in the phase, zone 0.0035 | **replayed, about 1.1x; closed** (`trims/README.md`): the deferral has nothing to save (the extrapolation rule already lets slow uncut fathoms finish), the primal gate 4:20 is in the code (default off, no false stop on 88 logged phases) and worth about 6%, the zone a modelled 4%; the search A/B was abandoned as not worth four hours | done | switch `QAP_DNN_CUT_PGATE=4:20` on in the production configuration; no further test |
| bound | skew-symmetric certificate polishing at non-converged iterates | 1.05-1.15x on the phase | 2-3 days | offline ascent on a dumped residual; killed under 20% of the gap recovered |
| bound | regional M2 block as a Lagrangian decomposition (Astra's 10b) | +20 units (10-15% chance of +65); **downgraded Nov 4 to about 0** (`new_cuts/SYNTHESIS.md`: after the triangle phase every single threshold is already saturated, so per-region blocks have nothing left) | 2-3 weeks | not worth running |
| bound | new cut families, three research agents (literature, grid, facility; `new_cuts/SYNTHESIS.md`) | **none meets +40; no path to +130 through linear inequalities on Y.** The residual gap is the global coupling of the middle thresholds (117-143 at depth 3), which only branching resolves. The one side-bet, large-Q box inequalities in the pool (+10 to +25 estimated), **failed its gate Nov 4: -32.3 on node 10, -22.6 on node 13** (`new_cuts/box_pool/`): no Lagrangian gain inside the Polyak rounds, a wider certificate lag | done | closed |
| speed | GPU port, 8-16 nodes per device | wall clock: 100 H100s = about 7x the 1,000 cores | 3-6 weeks | one-day CuPy prototype; dropped above 10 ms per node-iteration on an A100 |

Closed this week, do not revisit without new evidence: 5-point counts, triple
consistency, the dual method (bundle), the local-7 family as a pool cut type,
Anderson/Halpern, momentum, atomic triangles, eliminations from the
certificate. Another full level from new inequalities at depths 3-5 is rated
under 10% within a month (`bound_strength.md`). Closed Nov 4 by the three-agent study (`new_cuts/`): every facility-side family, the location-side geometric families (parity, dispersion, free regions, rectangles, l1/hypermetric), the row/column decomposition, objective splitting with a coarse optimum, the three local RLT2 projections, RLT2/RLT3/Lasserre-2, disjunctive and lift-and-project cuts, and the per-region moment block; and, after its own test, the large-Q box inequalities (`new_cuts/box_pool/`).

## 3. Multiplying out (revised Nov 3, evening)

The Nov 2 centers were rule-9 branching 3x, inheritance 2x, per-node speed
2x, phase trims and fleet efficiency 1.3x: about 15x. What has been measured
since (`[m]`) and what is still an estimate (`[e]`):

| lever | factor | status |
|---|---|---|
| rule 9 branching | 1.41x [m] | end to end on the node-10 subtree (`stronger_bounds/README.md` step 12) |
| per-node speed, bundle D | 1.71x [m] | per node on the 32-node sample, same decisions (`amx/README.md`) |
| parent-Y child warm start | 1.17x [m] | per child on 16 depth-4 children (`warm_y/README.md`); 25% on the laptop's full nug30 search with one-pass CholQR and two give-up settings |
| **all three in the search** | **2.34x [m]** | **end to end on the node-10 subtree** (`stronger_bounds/README.md` step 13): the same tree as rule 9 at 32,740 CPU-s against 76,755 for the baseline configuration; the per-node levers give 1.67x in the search against 2.0x measured separately (lost warm starts at checkpoints and restarts, and cut phases not measured per node) |
| pool inheritance in the search | 1.37x [m] | end to end on the node-10 subtree on top of the stacked arm (`pool_search/README.md`): same tree, 23,905 against 32,740 CPU-s |
| **all four in the search** | **3.21x [m]** | **against the baseline configuration's 76,755 CPU-s on the same subtree (5.8x against the plain search's 139,149)** |
| coordinator (fleet busy 70-85% -> 93-97%) | 1.15-1.35x [e] of wall clock; 1.09x [m] in iterations on the 4-slot subtree | built (`coord/`); the fleet factor is measured at plan.md phase 2 |
| shorter cut phase (Y12s), phase trims | 1.1-1.2x each [m on the phase] | not pursued |

Against the baseline of section 1 (120-700 days on 1,000 Xeon cores, center
290; a week needs about 41x at the center):

- **Measured end to end (3.21x, Nov 4, with inheritance):** center about
  **90 days**, range about 40-220 (it was 125 days at 2.34x before
  inheritance was measured).
- **Plus the coordinator at 1.25x (about 4.0x):** center about **72 days**,
  range about 30-175. A week is reached only at the low end of the tree
  range. (The Nov 3 projection assumed inheritance at 2x; it measured 1.37x.)
- **100 H100-class GPUs in place of the cores** (about 7x the wall-clock
  throughput at similar dollars, `node_speed.md`): center about 13 days on
  the measured levers alone, about 10 days with the coordinator.

**A week on 1,000 Xeon-class cores is not reachable at the center of the tree
range with the levers we have.** It comes into reach if the tree is at the
small end of its range (the 16 sampled depth-3 nodes branch on random lines,
the production search on rule-9 lines that fathom far more children), or on
GPUs. Caveat on the arithmetic: the 3.21x is measured end to end in the
search, but on one subtree (node 10); how much of the full tree's work sits
in the children of cut-phase misses, where inheritance acts, is for the
sampling plan to show.

The two measurements that decide between these outcomes are cheap:

1. **The sampling plan** (`tree_size.md` section 6): build the production
   frontier at depths 0-2 (1.6e5 CPU-s), take 400 frontier nodes to the
   give-up point, run stratified subtrees end to end, dive into the ten
   worst tail nodes, and run a cutoff ladder. About 1e7-2e7 CPU-s, 3-6 hours
   on 1,000 cores, 2-3% of the week, and every subtree finished is a piece of
   the proof. It pins the total within 1.5x instead of 15x.
2. **The subtree A/Bs** on the cloud box: rule 9 is done (1.41x); the
   stacked arm is done (2.34x end to end); inheritance in the search is
   done (1.37x on top, 3.21x in all, Nov 4).

## 4. Order of work

1. ~~Symmetry file for wil50~~ done Nov 3 (`sym/README.md`): every wil50 run
   from here uses `qaplib/wil50t.dat` with `QAP_DNN_CUT_SIDE=1` and `-F` specs
   as j:i; the node-10 subtree A/Bs are unaffected (no symmetry below node 10).
2. ~~Inheritance test and the rule-9 subtree A/B~~ done Nov 3: rule 9 1.41x;
   frozen-pool inheritance passed its gate; the stacked node-10 A/B gives
   2.34x end to end (Nov 4, step 13); inheritance in the search adds 1.37x
   on top (3.21x, Nov 4, `pool_search/`).
3. ~~Confirm give-ups uncut, the primal gate, zone 0.0035~~ closed Nov 3 by
   replay (`trims/README.md`, about 1.1x; the gate ships default off). The
   Y12s phase budget was not pursued.
4. ~~AMX and the CPU trims~~ done Nov 3 (1.71x per node); ~~the coordinator
   from docs/distributed/plan.md~~ built Nov 4 (phase 1, `coord/`; on the
   node-10 subtree the same tree with 1.09x fewer iterations,
   `coord_search/`). Next for it: plan.md phases 2-3 (Terraform, the Spot
   handler, 100 cores with induced interruptions) once an AWS account is
   connected.
5. The sampling plan on the first 1,000 cores, before the full run.
6. In parallel, gated: the M2 projection test and the GPU prototype.
